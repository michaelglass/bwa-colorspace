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

#ifndef BWT_BNTSEQ_H
#define BWT_BNTSEQ_H

#include <stdint.h>
#include <zlib.h>

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef uint8_t ubyte_t;
#endif

typedef struct {
	int64_t offset; //offset from beginning of pac methinks (slash beginning of buffer in multiseq BWS)
	int32_t len; //len of this seq
	int32_t n_ambs; //# missed codons (ie, something that's not ACTG)
	uint32_t gi; //???? initialized to 0 and never changed
	char *name, *anno; //name = name of sequence (ie, everything before first whitespace in initial comment)
	//anno is seq's comment (everything after that whitespace before a newline.)
} bntann1_t;

typedef struct {
	int64_t offset; //position of hole relative to beginning
	int32_t len; //length of contiguous hole containing the same wrong char 
	char amb; //the repeatedly wrong char
} bntamb1_t; 

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;

extern unsigned char nst_nt4_table[256];

#ifdef __cplusplus
extern "C" {
#endif

	void bns_dump(const bntseq_t *bns, const char *prefix);
	bntseq_t *bns_restore(const char *prefix);
	void bns_destroy(bntseq_t *bns);
	void bns_fasta2bntseq(gzFile fp_fa, const char *prefix);
	int bns_coor_pac2real(const bntseq_t *bns, int64_t pac_coor, int len, int32_t *real_seq);

#ifdef __cplusplus
}
#endif

#endif
