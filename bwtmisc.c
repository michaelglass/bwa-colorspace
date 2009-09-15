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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bntseq.h"
#include "utils.h"
#include "main.h"
#include "bwt.h"

#ifdef _DIVBWT
#include "divsufsort.h"
#endif

int is_bwt(ubyte_t *T, int n);

//returns the sequence length from the packed file.
//sequence is packed with a blank byte if length is divisible by 4
//last byte / blank byte is followed by a length of that byte.  if divisible by 4, the length is 0.
int64_t bwa_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	ubyte_t c;
	fp = xopen(fn_pac, "rb");
	fseek(fp, -1, SEEK_END); //go to second to last char
	pac_len = ftell(fp); //get curr position
	fread(&c, 1, 1, fp); //read 1 byte off the end
	fclose(fp);
	return (pac_len - 1) * 4 + (int)c; //return the length which is the length - 1 * 4 (4 bits per byte) + the last byte's "extra length"
}

//creates bwt for given packed sequence
//default use_is is true
bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is)
{
	bwt_t *bwt;
	ubyte_t *buf, *buf2;
	int i, pac_size;
	FILE *fp;

	// initialization
	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	bwt->seq_len = bwa_seq_len(fn_pac);
	bwt->bwt_size = (bwt->seq_len + 15) >> 4; //divide by 16 because we're using 4-byte slots not 1-byte slots
	//should be able to store 16 nucl per int instead of 4 per byte
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	//pac_size = enough space for the entire sequence, packed.
	// = len >> 2 + (len divisible by 4?) 0 : 1 extra
	pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
	buf2 = (ubyte_t*)calloc(pac_size, 1);
	fread(buf2, 1, pac_size, fp); //read the pac file into buf2
	fclose(fp);
	memset(bwt->L2, 0, 5 * 4); //clear memory in L2
	buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
	for (i = 0; i < bwt->seq_len; ++i) { //unpack buf into buf.  count all the chars occurence (ie, #C)
		buf[i] = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
		++bwt->L2[1+buf[i]]; //count the # of each nucleotide / color space.
	}
	for (i = 2; i <= 4; ++i) bwt->L2[i] += bwt->L2[i-1];
	free(buf2); //get rid of packed sequence

	// Burrows-Wheeler Transform
	if (use_is) {
		//only worrying about this, dont worry about libdivufsort
		bwt->primary = is_bwt(buf, bwt->seq_len); //index of terminator
	} else {
#ifdef _DIVBWT
		bwt->primary = divbwt(buf, buf, 0, bwt->seq_len);
#else
		err_fatal_simple("libdivsufsort is not compiled in.");
#endif
	}
	bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, 4); //allocate bwt_size ints ~= (seq_len / 16) bytes
	for (i = 0; i < bwt->seq_len; ++i)
		bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1); //pack bwt with bwt-transformed buffer
	free(buf);
	return bwt;
	//return the bwt
}

int bwa_pac2bwt(int argc, char *argv[])
{
	bwt_t *bwt;
	int c, use_is = 1;
	while ((c = getopt(argc, argv, "d")) >= 0) {
		switch (c) {
		case 'd': use_is = 0; break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa pac2bwt [-d] <in.pac> <out.bwt>\n");
		return 1;
	}
	bwt = bwt_pac2bwt(argv[optind], use_is);
	bwt_dump_bwt(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

/**
 returns the kth nucleotide from non-updated bwt_t *b
 @param bwt_t *b the source bwt.  cannot be updated
 @param int k the nucleotide-index of the nuceotide (as oppose to int array index)
 @return the nucleotide value of the kth nucleotide (0-3 representing a,c,t,g)
*/
#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

/**
 updates the bwt with char counts interspersed in the bwt's content
 The char counts are added every OCC_INTERVAL sequences, or OCC_INTERVAL/16 bytes.  OCC_INTERVAL defaults to 80 and must be divisible by 16 (to fit cleanly in bytespace).
 the char counts are in the form of four ints in the order A's count,C's count,G's,T's.

 @param bwt_t *bwt the source bwt to update
*/
void bwt_bwtupdate_core(bwt_t *bwt)
{
	bwtint_t i, k, c[4], n_occ;
	uint32_t *buf;

	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1; //default interval = 128
	bwt->bwt_size += n_occ * 4; // the new size
	buf = (uint32_t*)calloc(bwt->bwt_size, 4); // will be the new bwt
	c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = 0; i < bwt->seq_len; ++i) {
		if (i % OCC_INTERVAL == 0) { //starting at 0, every 80 chars, add the char counts (up to that point) to 4 ints in the bwt.
			memcpy(buf + k, c, sizeof(bwtint_t) * 4);
			k += 4;
		}
		if (i % 16 == 0) buf[k++] = bwt->bwt[i/16]; //then every 16, add the next 16 nucleotides from bwt to the buffer.
		++c[bwt_B00(bwt, i)];  //increments the C (not cumulative)
	}
	// then add the final char counts
	memcpy(buf + k, c, sizeof(bwtint_t) * 4);
	xassert(k + 4 == bwt->bwt_size, "inconsistent bwt_size");
	// update bwt
	free(bwt->bwt); bwt->bwt = buf;
}

int bwa_bwtupdate(int argc, char *argv[])
{
	bwt_t *bwt;
	if (argc < 2) {
		fprintf(stderr, "Usage: bwa bwtupdate <the.bwt>\n");
		return 1;
	}
	bwt = bwt_restore_bwt(argv[1]);
	bwt_bwtupdate_core(bwt);
	bwt_dump_bwt(argv[1], bwt);
	bwt_destroy(bwt);
	return 0;
}

void bwa_pac_rev_core(const char *fn, const char *fn_rev)
{
	int64_t seq_len, i;
	bwtint_t pac_len, j;
	ubyte_t *bufin, *bufout, ct;
	FILE *fp;
	seq_len = bwa_seq_len(fn);
	pac_len = (seq_len >> 2) + 1;
	bufin = (ubyte_t*)calloc(pac_len, 1);
	bufout = (ubyte_t*)calloc(pac_len, 1);
	fp = xopen(fn, "rb");
	fread(bufin, 1, pac_len, fp);
	fclose(fp);
	for (i = seq_len - 1, j = 0; i >= 0; --i) {
		int c = bufin[i>>2] >> ((~i&3)<<1) & 3; 
		//i = seq length.  buffer of 8-bit ints is actually
		//4 times smaller, so need to pull 4 chars from each byte.
		//c is last char, starting from the least significant 2 bits moving left 2 bits each time
		bwtint_t j = seq_len - 1 - i; //0 -> seq_len - 1
		bufout[j>>2] |= c << ((~j&3)<<1); //and then throw it in the first two bits, second, etc...
	}
	free(bufin);
	fp = xopen(fn_rev, "wb");
	fwrite(bufout, 1, pac_len, fp);
	ct = seq_len % 4; //pad with extra length
	fwrite(&ct, 1, 1, fp);
	fclose(fp);
	free(bufout);
}

int bwa_pac_rev(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: bwa pac_rev <in.pac> <out.pac>\n");
		return 1;
	}
	bwa_pac_rev_core(argv[1], argv[2]);
	return 0;
}

const int nst_color_space_table[] = { 4, 0, 0, 1, 0, 2, 3, 4, 0, 3, 2, 4, 1, 4, 4, 4}; //zeroes are for impossible values

/* this function is not memory efficient, but this will make life easier
   Ideally we should also change .amb files as one 'N' in the nucleotide
   sequence leads to two ambiguous colors. I may do this later... */
uint8_t *bwa_pac2cspac_core(const bntseq_t *bns)
{
	uint8_t *pac, *cspac; //sometimes it's a char*, ubyte_t*.  Here it's a uint8_t.  Arbitrary but equivalent
	bwtint_t i;
	int c1, c2;
	pac = (uint8_t*)calloc(bns->l_pac/4 + 1, 1);
	cspac = (uint8_t*)calloc(bns->l_pac/4 + 1, 1);
	fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);
	rewind(bns->fp_pac);
	c1 = pac[0]>>6; cspac[0] = c1<<6; //initial colorspace = header nucleotide
	for (i = 1; i < bns->l_pac; ++i) {
		c2 = pac[i>>2] >> (~i&3)*2 & 3; //pull nucleotide out of .pac.  each byte in pac has 4 except the last which may have 1-4 nucleotides.
		cspac[i>>2] |= nst_color_space_table[(1<<c1)|(1<<c2)] << (~i&3)*2; //set the cspac according to this one, last one.  fit 4 in each byte.
		/*
		a    c    g    t
		0    1    2    3		
	a 0 00-1 01-3 02-5 03-9 
	c 1      11-2 12-6 13-10
	g 2           22-4 23-12
	t 3                33-8
		*/
		c1 = c2;
	}
	free(pac);
	return cspac;
}

//convert .pac to colorspace .pac
int bwa_pac2cspac(int argc, char *argv[])
{
	bntseq_t *bns;
	uint8_t *cspac, ct;
	char *str;
	FILE *fp;

	if (argc < 3) {
		fprintf(stderr, "Usage: bwa pac2cspac <in.nt.prefix> <out.cs.prefix>\n");
		return 1;
	}
	bns = bns_restore(argv[1]);
	cspac = bwa_pac2cspac_core(bns);
	bns_dump(bns, argv[2]);
	// now write cspac
	str = (char*)calloc(strlen(argv[2]) + 5, 1);
	strcat(strcpy(str, argv[2]), ".pac");
	fp = xopen(str, "wb");
	fwrite(cspac, 1, bns->l_pac/4 + 1, fp);
	ct = bns->l_pac % 4; //cspac's n-1 char is already taken care of.  just need the last one.
	fwrite(&ct, 1, 1, fp);	
	fclose(fp);
	bns_destroy(bns);
	free(cspac);
	return 0;
}

int bwa_bwt2sa(int argc, char *argv[])
{
	bwt_t *bwt;
	int c, sa_intv = 32;
	while ((c = getopt(argc, argv, "i:")) >= 0) {
		switch (c) {
		case 'i': sa_intv = atoi(optarg); break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa bwt2sa [-i %d] <in.bwt> <out.sa>\n", sa_intv);
		return 1;
	}
	bwt = bwt_restore_bwt(argv[optind]);
	bwt_cal_sa(bwt, sa_intv);
	bwt_dump_sa(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}
