#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include "bwtaln.h"
#include "kvec.h"
#include "bntseq.h"
#include "utils.h"
#include "stdaln.h"

typedef struct {
	int n;
	bwtint_t *a;
} poslist_t;

typedef struct {
	double avg, std;
	bwtint_t low, high;
} isize_info_t;

#include "khash.h"
KHASH_MAP_INIT_INT64(64, poslist_t)

#include "ksort.h"
KSORT_INIT_GENERIC(uint64_t)

typedef struct {
	kvec_t(uint64_t) arr;
	kvec_t(uint64_t) pos[2];
	kvec_t(bwt_aln1_t) aln[2];
} pe_data_t;

#define MIN_HASH_WIDTH 1000

static int g_log_n[256];
static kh_64_t *g_hash;

void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);
void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq, bntseq_t *ntbns);
int bwa_approx_mapQ(const bwa_seq_t *p, int mm);
void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2);
bntseq_t *bwa_open_nt(const char *prefix);
void bwa_print_sam_SQ(const bntseq_t *bns);

pe_opt_t *bwa_init_pe_opt()
{
	pe_opt_t *po;
	po = (pe_opt_t*)calloc(1, sizeof(pe_opt_t));
	po->max_isize = 500;
	po->max_occ = 100000;
	po->type = BWA_PET_STD;
	po->is_sw = 1;
	return po;
}

static inline uint64_t hash_64(uint64_t key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}

// for normal distribution, this is about 3std
#define OUTLIER_BOUND 2.0

static int infer_isize(int n_seqs, bwa_seq_t *seqs[2], isize_info_t *ii)
{
	uint64_t x, *isizes;
	int n, i, tot, p25, p75, p50;

	ii->avg = ii->std = -1.0;
	ii->low = ii->high = 0;
	isizes = (uint64_t*)calloc(n_seqs, 8);
	for (i = 0, tot = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		p[0] = seqs[0] + i; p[1] = seqs[1] + i;
		if (p[0]->mapQ >= 20 && p[1]->mapQ >= 20)
			isizes[tot++] = (p[0]->pos < p[1]->pos)? p[1]->pos + p[1]->len - p[0]->pos : p[0]->pos + p[0]->len - p[1]->pos;
	}
	if (tot < 20) {
		fprintf(stderr, "[infer_isize] fail to infer insert size: too few good pairs\n");
		free(isizes);
		return -1;
	}
	ks_introsort(uint64_t, tot, isizes);
	p25 = isizes[(int)(tot*0.25 + 0.5)];
	p50 = isizes[(int)(tot*0.50 + 0.5)];
	p75 = isizes[(int)(tot*0.75 + 0.5)];
	ii->low  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
	ii->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
	for (i = 0, x = n = 0; i < tot; ++i)
		if (isizes[i] >= ii->low && isizes[i] <= ii->high)
			++n, x += isizes[i];
	ii->avg = (double)x / n;
	for (i = 0; i < tot; ++i)
		if (isizes[i] >= ii->low && isizes[i] <= ii->high)
			ii->std += (isizes[i] - ii->avg) * (isizes[i] - ii->avg);
	ii->std = sqrt(ii->std / n);
	free(isizes);
	if (ii->low <= 0 || isnan(ii->std)) {
		ii->low = ii->high = 0; ii->avg = ii->std = -1.0;
		fprintf(stderr, "[infer_isize] fail to infer insert size: weird pairing\n");
		return -1;
	}
	fprintf(stderr, "[infer_isize] (25, 50, 75) percentile: (%d, %d, %d)\n", p25, p50, p75);
	fprintf(stderr, "[infer_isize] low and high boundaries: %d and %d\n", ii->low, ii->high);
	fprintf(stderr, "[infer_isize] inferred external isize from %d pairs: %.3lf +/- %.3lf\n", n, ii->avg, ii->std);
	return 0;
}

static int pairing(bwa_seq_t *p[2], pe_data_t *d, const pe_opt_t *opt, int s_mm, const isize_info_t *ii)
{
	int i, j, o_n, subo_n, cnt_chg = 0;
	uint64_t last_pos[2][2], o_pos[2], subo_score, o_score;

	// here v>=u. When ii is set, we check insert size with ii; otherwise with opt->max_isize
#define __pairing_aux(u,v) do {											\
		bwtint_t l = ((v)>>32) + p[(v)&1]->len - ((u)>>32);				\
		if ((u) != (uint64_t)-1 && (v)>>32 > (u)>>32					\
			&& ((ii->high && l >= ii->low && l <= ii->high) || (ii->high == 0 && l <= opt->max_isize))) \
		{																\
			uint64_t s = d->aln[(v)&1].a[(uint32_t)(v)>>1].score + d->aln[(u)&1].a[(uint32_t)(u)>>1].score; \
			s *= 10;													\
			if (ii->high) s += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * fabs(l - ii->avg) / ii->std)) + .499); \
			s = s<<32 | (uint32_t)hash_64((u)>>32<<32 | (v)>>32);		\
			if (s>>32 == o_score>>32) ++o_n;							\
			else if (s>>32 < o_score>>32) { subo_n += o_n; o_n = 1; }	\
			else ++subo_n;												\
			if (s < o_score) subo_score = o_score, o_score = s, o_pos[(u)&1] = (u), o_pos[(v)&1] = (v); \
			else if (s < subo_score) subo_score = s;					\
		}																\
	} while (0)

#define __pairing_aux2(q, w) do {										\
		const bwt_aln1_t *r = d->aln[(w)&1].a + ((uint32_t)(w)>>1);		\
		(q)->extra_flag |= SAM_FPP;										\
		if ((q)->pos != (w)>>32 || (q)->strand != r->a) {				\
			(q)->n_mm = r->n_mm; (q)->n_gapo = r->n_gapo; (q)->n_gape = r->n_gape; (q)->strand = r->a; \
			(q)->score = r->score; (q)->mapQ = mapQ_p;					\
			(q)->pos = (w)>>32;											\
			if ((q)->mapQ > 0) ++cnt_chg;								\
		}																\
	} while (0)

	o_score = subo_score = (uint64_t)-1;
	o_n = subo_n = 0;
	ks_introsort(uint64_t, d->arr.n, d->arr.a);
	for (j = 0; j < 2; ++j) last_pos[j][0] = last_pos[j][1] = (uint64_t)-1;
	if (opt->type == BWA_PET_STD) {
		for (i = 0; i < d->arr.n; ++i) {
			uint64_t x = d->arr.a[i];
			int strand = d->aln[x&1].a[(uint32_t)x>>1].a;
			if (strand == 1) { // reverse strand, then check
				int y = 1 - (x&1);
				__pairing_aux(last_pos[y][1], x);
				__pairing_aux(last_pos[y][0], x);
			} else { // forward strand, then push
				last_pos[x&1][0] = last_pos[x&1][1];
				last_pos[x&1][1] = x;
			}
		}
	} else if (opt->type == BWA_PET_SOLID) {
		for (i = 0; i < d->arr.n; ++i) {
			uint64_t x = d->arr.a[i];
			int strand = d->aln[x&1].a[(uint32_t)x>>1].a;
			if ((strand^x)&1) { // push
				int y = 1 - (x&1);
				__pairing_aux(last_pos[y][1], x);
				__pairing_aux(last_pos[y][0], x);
			} else { // check
				last_pos[x&1][0] = last_pos[x&1][1];
				last_pos[x&1][1] = x;
			}
		}
	} else {
		fprintf(stderr, "[paring] not implemented yet!\n");
		exit(1);
	}
	// set pairing
	//fprintf(stderr, "[%d, %d, %d, %d]\n", d->arr.n, (int)(o_score>>32), (int)(subo_score>>32), o_n);
	if (o_score != (uint64_t)-1) {
		int mapQ_p = 0; // this is the maximum mapping quality when one end is moved
		//fprintf(stderr, "%d, %d\n", o_n, subo_n);
		if (o_n == 1) {
			if (subo_score == (uint64_t)-1) mapQ_p = 29; // no sub-optimal pair
			else if ((subo_score>>32) - (o_score>>32) > s_mm * 10) mapQ_p = 23; // poor sub-optimal pair
			else {
				int n = subo_n > 255? 255 : subo_n;
				mapQ_p = ((subo_score>>32) - (o_score>>32)) / 2 - g_log_n[n];
				if (mapQ_p < 0) mapQ_p = 0;
			}
		}
		if (p[0]->pos == o_pos[0]>>32 && p[1]->pos == o_pos[1]>>32) { // both ends not moved
			if (p[0]->mapQ > 0 && p[1]->mapQ > 0) {
				int mapQ = p[0]->mapQ + p[1]->mapQ;
				if (mapQ > 60) mapQ = 60;
				p[0]->mapQ = p[1]->mapQ = mapQ;
			} else {
				if (p[0]->mapQ == 0) p[0]->mapQ = (mapQ_p + 7 < p[1]->mapQ)? mapQ_p + 7 : p[1]->mapQ;
				if (p[1]->mapQ == 0) p[1]->mapQ = (mapQ_p + 7 < p[0]->mapQ)? mapQ_p + 7 : p[0]->mapQ;
			}
		} else if (p[0]->pos == o_pos[0]>>32) { // [1] moved
			p[1]->seQ = 0; p[1]->mapQ = p[0]->mapQ;
			if (p[1]->mapQ > mapQ_p) p[1]->mapQ = mapQ_p;
		} else if (p[1]->pos == o_pos[1]>>32) { // [0] moved
			p[1]->seQ = 0; p[0]->mapQ = p[1]->mapQ;
			if (p[0]->mapQ > mapQ_p) p[0]->mapQ = mapQ_p;
		} else { // both ends moved
			p[0]->seQ = p[1]->seQ = 0;
			mapQ_p -= 20;
			if (mapQ_p < 0) mapQ_p = 0;
			p[0]->mapQ = p[1]->mapQ = mapQ_p;
		}
		__pairing_aux2(p[0], o_pos[0]);
		__pairing_aux2(p[1], o_pos[1]);
	}
	return cnt_chg;
}

typedef struct {
	kvec_t(bwt_aln1_t) aln;
} aln_buf_t;

int bwa_cal_pac_pos_pe(const char *prefix, bwt_t *const _bwt[2], int n_seqs, bwa_seq_t *seqs[2], FILE *fp_sa[2], isize_info_t *ii,
					   const pe_opt_t *opt, const gap_opt_t *gopt, const isize_info_t *last_ii)
{
	int i, j, cnt_chg = 0;
	char str[1024];
	bwt_t *bwt[2];
	pe_data_t *d;
	aln_buf_t *buf[2];

	d = (pe_data_t*)calloc(1, sizeof(pe_data_t));
	buf[0] = (aln_buf_t*)calloc(n_seqs, sizeof(aln_buf_t));
	buf[1] = (aln_buf_t*)calloc(n_seqs, sizeof(aln_buf_t));

	if (_bwt[0] == 0) { // load forward SA
		strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt[0]);
		strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt[1]);
	} else bwt[0] = _bwt[0], bwt[1] = _bwt[1];

	// SE
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {
			int n_aln;
			p[j] = seqs[j] + i;
			p[j]->extra_flag |= SAM_FPD | (j == 0? SAM_FR1 : SAM_FR2);
			fread(&n_aln, 4, 1, fp_sa[j]);
			if (n_aln > kv_max(d->aln[j]))
				kv_resize(bwt_aln1_t, d->aln[j], n_aln);
			d->aln[j].n = n_aln;
			fread(d->aln[j].a, sizeof(bwt_aln1_t), n_aln, fp_sa[j]);
			kv_copy(bwt_aln1_t, buf[j][i].aln, d->aln[j]); // backup d->aln[j]
			// generate SE alignment and mapping quality
			bwa_aln2seq(n_aln, d->aln[j].a, p[j]);
			if (p[j]->type == BWA_TYPE_UNIQUE || p[j]->type == BWA_TYPE_REPEAT) {
				int max_diff = gopt->fnr > 0.0? bwa_cal_maxdiff(p[j]->len, BWA_AVG_ERR, gopt->fnr) : gopt->max_diff;
				p[j]->pos = p[j]->strand? bwt_sa(bwt[0], p[j]->sa)
					: bwt[1]->seq_len - (bwt_sa(bwt[1], p[j]->sa) + p[j]->len);
				p[j]->seQ = p[j]->mapQ = bwa_approx_mapQ(p[j], max_diff);
			}
		}
	}

	// infer isize
	infer_isize(n_seqs, seqs, ii);
	if (ii->avg < 0.0 && last_ii->avg > 0.0) *ii = *last_ii;

	// PE
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {
			p[j] = seqs[j] + i;
			kv_copy(bwt_aln1_t, d->aln[j], buf[j][i].aln);
		}
		if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
			&& (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT))
		{ // only when both ends mapped
			uint64_t x;
			int j, k, n_occ[2];
			for (j = 0; j < 2; ++j) {
				n_occ[j] = 0;
				for (k = 0; k < d->aln[j].n; ++k)
					n_occ[j] += d->aln[j].a[k].l - d->aln[j].a[k].k + 1;
			}
			if (n_occ[0] > opt->max_occ || n_occ[1] > opt->max_occ) continue;
			d->arr.n = 0;
			for (j = 0; j < 2; ++j) {
				for (k = 0; k < d->aln[j].n; ++k) {
					bwt_aln1_t *r = d->aln[j].a + k;
					bwtint_t l;
					if (r->l - r->k + 1 >= MIN_HASH_WIDTH) { // then check hash table
						uint64_t key = (uint64_t)r->k<<32 | r->l;
						int ret;
						khint_t iter = kh_put(64, g_hash, key, &ret);
						if (ret) { // not in the hash table; ret must equal 1 as we never remove elements
							poslist_t *z = &kh_val(g_hash, iter);
							z->n = r->l - r->k + 1;
							z->a = (bwtint_t*)malloc(sizeof(bwtint_t) * z->n);
							for (l = r->k; l <= r->l; ++l)
								z->a[l - r->k] = r->a? bwt_sa(bwt[0], l) : bwt[1]->seq_len - (bwt_sa(bwt[1], l) + p[j]->len);
						}
						for (l = 0; l < kh_val(g_hash, iter).n; ++l) {
							x = kh_val(g_hash, iter).a[l];
							x = x<<32 | k<<1 | j;
							kv_push(uint64_t, d->arr, x);
						}
					} else { // then calculate on the fly
						for (l = r->k; l <= r->l; ++l) {
							x = r->a? bwt_sa(bwt[0], l) : bwt[1]->seq_len - (bwt_sa(bwt[1], l) + p[j]->len);
							x = x<<32 | k<<1 | j;
							kv_push(uint64_t, d->arr, x);
						}
					}
				}
			}
			cnt_chg += pairing(p, d, opt, gopt->s_mm, ii);
		}
	}

	// free
	for (i = 0; i < n_seqs; ++i) {
		kv_destroy(buf[0][i].aln);
		kv_destroy(buf[1][i].aln);
	}
	free(buf[0]); free(buf[1]);
	if (_bwt[0] == 0) {
		bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	}
	kv_destroy(d->arr);
	kv_destroy(d->pos[0]); kv_destroy(d->pos[1]);
	kv_destroy(d->aln[0]); kv_destroy(d->aln[1]);
	free(d);
	return cnt_chg;
}

#define SW_MIN_MATCH_LEN 17

// cnt = n_mm<<16 | n_gapo<<8 | n_gape
uint16_t *bwa_sw_core(bwtint_t l_pac, const ubyte_t *pacseq, int len, const ubyte_t *seq, int64_t *beg, int reglen,
					  int *n_cigar, uint32_t *_cnt)
{
	uint16_t *cigar = 0;
	ubyte_t *ref_seq;
	bwtint_t k, x, y, l;
	int path_len, ret;
	AlnParam ap = aln_param_bwa;
	path_t *path, *p;

	// check whether there are too many N's
	if (reglen < SW_MIN_MATCH_LEN || (int64_t)l_pac - *beg < len) return 0;
	for (k = 0, x = 0; k < len; ++k)
		if (seq[k] >= 4) ++x;
	if ((float)x/len >= 0.25 || len - x < SW_MIN_MATCH_LEN) return 0;

	// get reference subsequence
	ref_seq = (ubyte_t*)calloc(reglen, 1);
	for (k = *beg, l = 0; l < reglen && k < l_pac; ++k)
		ref_seq[l++] = pacseq[k>>2] >> ((~k&3)<<1) & 3;
	path = (path_t*)calloc(l+len, sizeof(path_t));

	// do alignment
	ret = aln_local_core(ref_seq, l, (ubyte_t*)seq, len, &ap, path, &path_len, 1, 0);
	if (ret < 0) {
		free(path); free(cigar); *n_cigar = 0;
		return 0;
	}
	cigar = aln_path2cigar(path, path_len, n_cigar);

	// check whether the alignment is good enough
	for (k = 0, x = y = 0; k < *n_cigar; ++k) {
		uint16_t c = cigar[k];
		if (c>>14 == FROM_M) x += c&0x3fff, y += c&0x3fff;
		else if (c>>14 == FROM_D) x += c&0x3fff;
		else y += c&0x3fff;
	}
	if (x < SW_MIN_MATCH_LEN || y < SW_MIN_MATCH_LEN) { // not good enough
		free(path); free(cigar);
		*n_cigar = 0;
		return 0;
	}

	{ // update cigar and coordinate;
		int start, end;
		p = path + path_len - 1;
		*beg += (p->i? p->i : 1) - 1;
		start = (p->j? p->j : 1) - 1;
		end = path->j;
		cigar = (uint16_t*)realloc(cigar, 2 * (*n_cigar + 2));
		if (start) {
			memmove(cigar + 1, cigar, 2 * (*n_cigar));
			cigar[0] = 3<<14 | start;
			++(*n_cigar);
		}
		if (end < len) {
			cigar[*n_cigar] = 3<<14 | (len - end);
			++(*n_cigar);
		}
	}

	{ // set *cnt
		int n_mm, n_gapo, n_gape;
		n_mm = n_gapo = n_gape = 0;
		p = path + path_len - 1;
		x = p->i? p->i - 1 : 0; y = p->j? p->j - 1 : 0;
		for (k = 0; k < *n_cigar; ++k) {
			uint16_t c = cigar[k];
			if (c>>14 == FROM_M) {
				for (l = 0; l < (c&0x3fff); ++l)
					if (ref_seq[x+l] < 4 && seq[y+l] < 4 && ref_seq[x+l] != seq[y+l]) ++n_mm;
				x += c&0x3fff, y += c&0x3fff;
			} else if (c>>14 == FROM_D) {
				x += c&0x3fff, ++n_gapo, n_gape += (c&0x3fff) - 1;
			} else if (c>>14 == FROM_I) {
				y += c&0x3fff, ++n_gapo, n_gape += (c&0x3fff) - 1;
			}
		}
		*_cnt = (uint32_t)n_mm<<16 | n_gapo<<8 | n_gape;
	}
	
	free(ref_seq); free(path);
	return cigar;
}

ubyte_t *bwa_paired_sw(const bntseq_t *bns, const ubyte_t *_pacseq, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii)
{
	ubyte_t *pacseq;
	int i;
	uint64_t x, n;

	// load reference sequence
	if (_pacseq == 0) {
		pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
		rewind(bns->fp_pac);
		fread(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);
	} else pacseq = (ubyte_t*)_pacseq;
	if (!popt->is_sw || ii->avg < 0.0) return pacseq;

	// perform mate alignment
	for (i = 0, x = n = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		int is_first = 1;
		p[0] = seqs[0] + i; p[1] = seqs[1] + i;
		if ((p[0]->mapQ >= 20 && p[1]->type == BWA_TYPE_NO_MATCH) || (p[1]->mapQ >= 20 && p[0]->type == BWA_TYPE_NO_MATCH)) {
			++n;
			if (p[0]->type == BWA_TYPE_NO_MATCH) {
				p[0] = seqs[1] + i; p[1] = seqs[0] + i; // swap s.t p[0] is the mapped read
				is_first = 0;
			}
			if (popt->type == BWA_PET_STD || popt->type == BWA_PET_SOLID) {
				int64_t beg, end; // this is the start and end of the region
				ubyte_t *seq;
				uint32_t cnt;

#define __set_rght_coor(_a, _b) do {									\
					(_a) = p[0]->pos + ii->avg - 3 * ii->std - p[1]->len * 1.5; \
					(_b) = (_a) + 6 * ii->std + 2 * p[1]->len;			\
					if ((_a) < p[0]->pos + p[0]->len) (_a) = p[0]->pos + p[0]->len; \
					if ((_b) > bns->l_pac) (_b) = bns->l_pac;			\
				} while (0)

#define __set_left_coor(_a, _b) do {									\
					(_a) = p[0]->pos + p[0]->len - ii->avg - 3 * ii->std - p[1]->len * 0.5; \
					(_b) = (_a) + 6 * ii->std + 2 * p[1]->len;			\
					if ((_a) < 0) (_a) = 0;								\
					if ((_b) > p[0]->pos) (_b) = p[0]->pos;				\
				} while (0)

				if (popt->type == BWA_PET_STD) {
					if (p[0]->strand == 0) { // the mate is on the reverse strand and has larger coordinate
						__set_rght_coor(beg, end);
						seq = p[1]->rseq;
					} else { // the mate is on forward stand and has smaller coordinate
						__set_left_coor(beg, end);
						seq = p[1]->seq;
						seq_reverse(p[1]->len, seq, 0); // because ->seq is reversed
					}
				} else { // popt->type == BWA_PET_SOLID
					if (p[0]->strand == 0) {
						if (!is_first) __set_left_coor(beg, end);
						else __set_rght_coor(beg, end);
						seq = p[1]->rseq;
						seq_reverse(p[1]->len, seq, 0); // because ->seq is reversed
					} else {
						if (!is_first) __set_rght_coor(beg, end);
						else __set_left_coor(beg, end);
						seq = p[1]->seq;
					}
				}

				p[1]->cigar = bwa_sw_core(bns->l_pac, pacseq, p[1]->len, seq, &beg, end - beg, &p[1]->n_cigar, &cnt);
				if (p[1]->cigar) { // the SW alignment is good enough
					++x;
					p[1]->type = BWA_TYPE_MATESW;
					p[1]->pos = beg;
					p[1]->mapQ = p[0]->mapQ;
					p[1]->seQ = p[0]->seQ;
					p[1]->strand = (popt->type == BWA_PET_STD)? 1 - p[0]->strand : p[0]->strand;
					p[1]->n_mm = cnt>>16; p[1]->n_gapo = cnt>>8&0xff; p[1]->n_gape = cnt&0xff;
					p[1]->extra_flag |= SAM_FPP;
					p[0]->extra_flag |= SAM_FPP;
				}
				if (popt->type == BWA_PET_STD) {
					if (p[0]->strand) seq_reverse(p[1]->len, seq, 0); // reverse it back
				} else {
					if (p[0]->strand == 0) seq_reverse(p[1]->len, seq, 0);
				}
			} else {
				fprintf(stderr, "[bwa_paired_sw] not implemented!\n");
				exit(1);
			}
		}
	}
	fprintf(stderr, "[bwa_paired_sw] %lld reads aligned out of %lld candidates.\n",
			(long long)x, (long long)n);

	return pacseq;
}

void bwa_sai2sam_pe_core(const char *prefix, char *const fn_sa[2], char *const fn_fa[2], pe_opt_t *popt)
{
	int i, j, n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs[2];
	bwa_seqio_t *ks[2];
	clock_t t;
	bntseq_t *bns, *ntbns = 0;
	FILE *fp_sa[2];
	gap_opt_t opt;
	khint_t iter;
	isize_info_t last_ii; // this is for the last batch of reads
	char str[1024];
	bwt_t *bwt[2];
	uint8_t *pac;

	// initialization
	pac = 0; bwt[0] = bwt[1] = 0;
	for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);
	bns = bns_restore(prefix);
	srand48(bns->seed);
	for (i = 0; i < 2; ++i) {
		ks[i] = bwa_seq_open(fn_fa[i]);
		fp_sa[i] = xopen(fn_sa[i], "r");
	}
	g_hash = kh_init(64);
	last_ii.avg = -1.0;

	fread(&opt, sizeof(gap_opt_t), 1, fp_sa[0]);
	fread(&opt, sizeof(gap_opt_t), 1, fp_sa[1]);
	if (!(opt.mode & BWA_MODE_COMPREAD)) {
		popt->type = BWA_PET_SOLID;
		ntbns = bwa_open_nt(prefix);
	} else { // for Illumina alignment only
		if (popt->is_preload) {
			strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
			strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt[0]);
			strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
			strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt[1]);
			pac = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
			rewind(bns->fp_pac);
			fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);
		}
	}

	// core loop
	bwa_print_sam_SQ(bns);
	while ((seqs[0] = bwa_read_seq(ks[0], 0x40000, &n_seqs, opt.mode & BWA_MODE_COMPREAD, opt.trim_qual)) != 0) {
		int cnt_chg;
		isize_info_t ii;
		ubyte_t *pacseq;

		seqs[1] = bwa_read_seq(ks[1], 0x40000, &n_seqs, opt.mode & BWA_MODE_COMPREAD, opt.trim_qual);
		tot_seqs += n_seqs;
		t = clock();

		fprintf(stderr, "[bwa_sai2sam_pe_core] convert to sequence coordinate... \n");
		cnt_chg = bwa_cal_pac_pos_pe(prefix, bwt, n_seqs, seqs, fp_sa, &ii, popt, &opt, &last_ii);
		fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
		fprintf(stderr, "[bwa_sai2sam_pe_core] change of coordinates in %d alignments.\n", cnt_chg);

		fprintf(stderr, "[bwa_sai2sam_pe_core] align unmapped mate...\n");
		pacseq = bwa_paired_sw(bns, pac, n_seqs, seqs, popt, &ii);
		fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_sai2sam_pe_core] refine gapped alignments... ");
		for (j = 0; j < 2; ++j)
			bwa_refine_gapped(bns, n_seqs, seqs[j], pacseq, ntbns);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
		if (pac == 0) free(pacseq);

		fprintf(stderr, "[bwa_sai2sam_pe_core] print alignments... ");
		for (i = 0; i < n_seqs; ++i) {
			bwa_print_sam1(bns, seqs[0] + i, seqs[1] + i, opt.mode, opt.max_top2);
			bwa_print_sam1(bns, seqs[1] + i, seqs[0] + i, opt.mode, opt.max_top2);
		}
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		for (j = 0; j < 2; ++j)
			bwa_free_read_seq(n_seqs, seqs[j]);
		fprintf(stderr, "[bwa_sai2sam_pe_core] %d sequences have been processed.\n", tot_seqs);
		last_ii = ii;
	}

	// destroy
	bns_destroy(bns);
	if (ntbns) bns_destroy(ntbns);
	for (i = 0; i < 2; ++i) {
		bwa_seq_close(ks[i]);
		fclose(fp_sa[i]);
	}
	for (iter = kh_begin(g_hash); iter != kh_end(g_hash); ++iter)
		if (kh_exist(g_hash, iter)) free(kh_val(g_hash, iter).a);
	kh_destroy(64, g_hash);
	if (pac) {
		free(pac); bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	}
}

int bwa_sai2sam_pe(int argc, char *argv[])
{
	int c;
	pe_opt_t *popt;
	popt = bwa_init_pe_opt();
	while ((c = getopt(argc, argv, "a:o:sP")) >= 0) {
		switch (c) {
		case 'a': popt->max_isize = atoi(optarg); break;
		case 'o': popt->max_occ = atoi(optarg); break;
		case 's': popt->is_sw = 0; break;
		case 'P': popt->is_preload = 1; break;
		default: return 1;
		}
	}

	if (optind + 5 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa sampe [options] <prefix> <in1.sai> <in2.sai> <in1.fq> <in2.fq>\n\n");
		fprintf(stderr, "Options: -a INT   maximum insert size [%d]\n", popt->max_isize);
		fprintf(stderr, "         -o INT   maximum occurrences for one end [%d]\n", popt->max_occ);
		fprintf(stderr, "         -P       preload index into memory (for base-space reads only)\n");
		fprintf(stderr, "         -s       disable Smith-Waterman for the unmapped mate\n\n");
		fprintf(stderr, "Notes: 1. For SOLiD read, <in1.fq> corresponds R3 reads and <in2.fq> to F3.\n");
		fprintf(stderr, "       2. For reads shorter than 30bp, applying a smaller -o is recommended to\n");
		fprintf(stderr, "          to get a sensible speed at the cost of pairing accuracy.\n");
		fprintf(stderr, "\n");
		return 1;
	}
	bwa_sai2sam_pe_core(argv[optind], argv + optind + 1, argv + optind+3, popt);
	free(popt);
	return 0;
}
