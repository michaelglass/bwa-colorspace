#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "bwt.h"
#include "utils.h"

void bwt_dump_bwt(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp); //write the primary index into the first int
	fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp); //write the cum # of As, ACs, ACGs, and ACGTs into the next four ints
	fwrite(bwt->bwt, sizeof(bwtint_t), bwt->bwt_size, fp); //dump the bwt-permuted sequence
	fclose(fp);
}

void bwt_dump_sa(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	fwrite(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	fwrite(&bwt->seq_len, sizeof(bwtint_t), 1, fp);
	fwrite(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
	fclose(fp);
}

void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = xopen(fn, "rb");
	fread(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	fread(skipped, sizeof(bwtint_t), 4, fp); // skip
	fread(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	fread(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	bwt->sa[0] = -1;

	fread(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
	fclose(fp);
}

bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = xopen(fn, "rb");
	fseek(fp, 0, SEEK_END); //this is to get file length, seeking to the beginning
	//not completely sure why this does this, 
	// I'm going ot assume the bwt is padded with 5 bwt_ints and the rest of the data structures are 4 bits each.
	bwt->bwt_size = (ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4); //and then allocating space for the bwt given that size
	fseek(fp, 0, SEEK_SET); //go back to the beginning
	fread(&bwt->primary, sizeof(bwtint_t), 1, fp); //read one bwtint into primary
	fread(bwt->L2+1, sizeof(bwtint_t), 4, fp); //put 4 bwtints into l2[1-5]
	fread(bwt->bwt, 4, bwt->bwt_size, fp);// read the bwt into bwt->bwt
	bwt->seq_len = bwt->L2[4];	//set the seq_len . . . what's l2[0]?
	fclose(fp);
	bwt_gen_cnt_table(bwt);

	return bwt;
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}
