/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "bntseq.h"
#include "main.h"
#include "utils.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char nst_nt4_table[256] = {
/*0*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
/*16*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
/*32*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
/*48*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
/*64*/	4, 0 /*'A'*/, 4, 1 /*'C'*/,  4, 4, 4, 2 /*'G'*/,  4, 4, 4, 4,  4, 4, 4, 4, 
/*80*/	4, 4, 4, 4,  3 /*'T'*/, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
/*96*/	4, 0 /*'a'*/, 4, 1/*'c'*/,  4, 4, 4, 2/*'g'*/,  4, 4, 4, 4,  4, 4, 4, 4, 
/*112*/	4, 4, 4, 4,  3/*'t'*/, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
/*128*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
/*144*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
/*160*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
/*176*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
/*192*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
/*208*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
/*224*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
/*240*/	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 /*255*/
};

void bns_dump(const bntseq_t *bns, const char *prefix)
{
	char str[1024];
	FILE *fp;
	int i;
	{ // dump .ann
		strcpy(str, prefix); strcat(str, ".ann");
		fp = xopen(str, "w");
		fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->seed); //pac length, # sequences, seed for random generator
		for (i = 0; i != bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			fprintf(fp, "%d %s", p->gi, p->name); //p->gi?  I don't think this is ever set to anything besides 0.  so I have no idea what this is.
			if (p->anno[0]) fprintf(fp, " %s\n", p->anno); //if there's a comment print it
			else fprintf(fp, "\n"); //else just print a newline
			fprintf(fp, "%lld %d %d\n", (long long)p->offset, p->len, p->n_ambs); //now print the offset and the lenghth of each string as well as the n_ambs. . . n_ambs = missing codons
		}
		fclose(fp);
	}
	{ // dump .amb
		strcpy(str, prefix); strcat(str, ".amb");
		fp = xopen(str, "w");
		fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->n_holes); //length of entire pac, num sequences, number of holes (number of places where ACTG = missing)
		for (i = 0; i != bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			fprintf(fp, "%lld %d %c\n", (long long)p->offset, p->len, p->amb); //the offset, the length, and the char that was repeatedly wrong (should be a '-')
		}
		fclose(fp);
	}
}

bntseq_t *bns_restore(const char *prefix)
{
	char str[1024];
	FILE *fp;
	bntseq_t *bns;
	long long xx;
	int i;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		strcpy(str, prefix); strcat(str, ".ann");
		fp = xopen(str, "r");
		fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			fscanf(fp, "%u%s", &p->gi, str);
			p->name = strdup(str);
			// read fasta comments 
			while ((c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			*q = 0;
			if (q - str > 1) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			p->offset = xx;
		}
		fclose(fp);
	}
	{ // read .amb
		int64_t l_pac;
		int32_t n_seqs;
		strcpy(str, prefix); strcat(str, ".amb");
		fp = xopen(str, "r");
		fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		l_pac = xx;
		xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t));
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			p->offset = xx;
			p->amb = str[0];
		}
		fclose(fp);
	}
	{ // open .pac
		strcpy(str, prefix); strcat(str, ".pac");
		bns->fp_pac = xopen(str, "rb");
	}
	return bns;
}

void bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else {
		int i;
		if (bns->fp_pac) fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i) {
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}

// this converts fasta to pac file.  prefix is filename.nt, fp_fa is fasta source
void bns_fasta2bntseq(gzFile fp_fa, const char *prefix)
{
	kseq_t *seq;
	//kseq_t is
	// 4 kstring_t's (2 size_ts and a string (char*))
	// 1 int 
	// 1 *kstreamt (char*; int begin,end,is_eof; gzFile f)
	char name[1024]; //why not do strlen like everywhere else?
	bntseq_t *bns;
	bntamb1_t *q;
	int l_buf;
	unsigned char buf[0x10000]; //2^16 //bit buffer
	int32_t m_seqs, m_holes, l, i;
	FILE *fp;

	// initialization
	seq = kseq_init(fp_fa);
	//creates a kseq_t, puts fp_fa in seq->f->f, sets seq->f->buf to be a 4096-long char*
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	bns->seed = 11; // fixed seed for random generator
	srand48(bns->seed);
	m_seqs = m_holes = 8;
	bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
	bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
	q = bns->ambs;
	l_buf = 0;
	strcpy(name, prefix); strcat(name, ".pac");
	fp = xopen(name, "wb");
	memset(buf, 0, 0x10000); //zeroes out buf.  could have just used calloc?
	// read sequences
	while ((l = kseq_read(seq)) >= 0) {
		bntann1_t *p;
		int lasts;
		if (bns->n_seqs == m_seqs) { //if we're out of space	
			m_seqs <<= 1; //double m_seqs
			bns->anns = (bntann1_t*)realloc(bns->anns, m_seqs * sizeof(bntann1_t)); //double size
		}

		p = bns->anns + bns->n_seqs; //cycling through each of the anns... what's an ann?  dependant on bns's sequences which start unpopulated.  
		p->name = strdup((char*)seq->name.s);
		p->anno = seq->comment.s? strdup((char*)seq->comment.s) : strdup("(null)"); //if we don't have a comment copy in a filler string
		p->gi = 0; p->len = l;  //l = length of sequence seq
		p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len; //first time around, bns->n_seqs == 0
		p->n_ambs = 0;
		for (i = 0, lasts = 0; i < l; ++i) { //i = index of sequence
			int c = nst_nt4_table[(int)seq->seq.s[i]]; //translate char to 0-3
			if (c >= 4) { // N  c >= 4 means it's not A,C,T or G
				if (lasts == seq->seq.s[i]) { // contiguous N
					++q->len;
				} else { //singular hole(hole?)
					if (bns->n_holes == m_holes) {
						m_holes <<= 1; //double memory
						bns->ambs = (bntamb1_t*)realloc(bns->ambs, m_holes * sizeof(bntamb1_t));
					}
					q = bns->ambs + bns->n_holes;
					q->len = 1; //when will len be greater than 1 for amb that can only be 1 char??  don't totally understand how this works
					q->offset = p->offset + i;
					q->amb = seq->seq.s[i];
					++p->n_ambs;
					++bns->n_holes;
				}
			} //N is a hole methinks.  
			lasts = seq->seq.s[i];
			{ // fill buffer
				if (c >= 4) c = lrand48()&0x3; //if it's not ACTG, pick a random one
				if (l_buf == 0x40000) { //buf is only 10000 chars long, but each char contains 4 codons, 
					fwrite(buf, 1, 0x10000, fp);  //write buf
					memset(buf, 0, 0x10000); //reset buf. cool!
					l_buf = 0; //lbuf is the length of buf.
				}
				//loads each 2-bit codon into the buffer from left to right.
				//ACTG =
				//0132 =
				//0b00011110
				//byte = 0
				//len = 1
				//load A into byte by setting byte |= a << 6 (bits 76 = 00)
				//load C into byte by setting byte |= c << 4 (bits 54 = 01)
				//load T into byte by setting byte |= T << 2 (bits 32 = 11)
				//load G into byte by setting byte |= G << 0 (bits 10 = 10)
				
				buf[l_buf>>2] |= c << ( (3 - (l_buf&3)) << 1 ); 
				
				++l_buf;
			}
		}
		++bns->n_seqs;
		bns->l_pac += seq->seq.l; //l_pac = total length of all the sequences summed
	}
	xassert(bns->l_pac, "zero length sequence."); //if bns->l_pac == 0 quit with message "zero..."
	{ // finalize .pac file
		ubyte_t ct;
		fwrite(buf, 1, (l_buf>>2) + ((l_buf&3) == 0? 0 : 1), fp); //write out the rest of the buffer, and 1 more if buf doesn't cleanly contain multiple of 4 code-ons
		// the following codes make the pac file size always (l_pac/4+1+1)
		if (bns->l_pac % 4 == 0) {
			ct = 0;
			fwrite(&ct, 1, 1, fp);
		}
		ct = bns->l_pac % 4;
		fwrite(&ct, 1, 1, fp);
		// close .pac file
		fclose(fp);
	}
	bns_dump(bns, prefix);
	bns_destroy(bns);
	kseq_destroy(seq);
}

int bwa_fa2pac(int argc, char *argv[])
{
	gzFile fp;
	if (argc < 2) {
		fprintf(stderr, "Usage: bwa fa2pac <in.fasta> [<out.prefix>]\n");
		return 1;
	}
	fp = xzopen(argv[1], "r");
	bns_fasta2bntseq(fp, (argc < 3)? argv[1] : argv[2]);
	gzclose(fp);
	return 0;
}

int bns_coor_pac2real(const bntseq_t *bns, int64_t pac_coor, int len, int32_t *real_seq)
{
	int left, mid, right, nn;
	if (pac_coor >= bns->l_pac)
		err_fatal("bns_coor_pac2real", "bug! Coordinate is longer than sequence (%lld>=%lld).", pac_coor, bns->l_pac);
	// binary search for the sequence ID. Note that this is a bit different from the following one...
	left = 0; mid = 0; right = bns->n_seqs;
	while (left < right) {
		mid = (left + right) >> 1;
		if (pac_coor >= bns->anns[mid].offset) {
			if (mid == bns->n_seqs - 1) break;
			if (pac_coor < bns->anns[mid+1].offset) break;
			left = mid + 1;
		} else right = mid;
	}
	*real_seq = mid;
	// binary search for holes
	left = 0; right = bns->n_holes; nn = 0;
	while (left < right) {
		int64_t mid = (left + right) >> 1;
		if (pac_coor >= bns->ambs[mid].offset + bns->ambs[mid].len) left = mid + 1;
		else if (pac_coor + len <= bns->ambs[mid].offset) right = mid;
		else { // overlap
			if (pac_coor >= bns->ambs[mid].offset) {
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pac_coor + len?
					bns->ambs[mid].offset + bns->ambs[mid].len - pac_coor : len;
			} else {
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pac_coor + len?
					bns->ambs[mid].len : len - (bns->ambs[mid].offset - pac_coor);
			}
			break;
		}
	}
	return nn;
}
