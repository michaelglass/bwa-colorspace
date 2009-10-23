/* The MIT License

   Copyright (c) 2008, by Heng Li <lh3@sanger.ac.uk>

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

#ifndef AC_KSEQ_H
#define AC_KSEQ_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

//contains 
#define __KS_TYPE(type_t)						\
	typedef struct __kstream_t {				\
		char *buf;								\
		int begin, end, is_eof;					\
		type_t f;								\
	} kstream_t;

#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)
//__bufsize is 4096, defined at the end of this file in KSEQ_INIT.  type_t is 
#define __KS_BASIC(type_t, __bufsize)								\
	static inline kstream_t *ks_init(type_t f)						\
	{																\
		kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));	\
		ks->f = f;													\
		ks->buf = (char*)malloc(__bufsize);							\
		return ks;													\
	}																\
	static inline void ks_destroy(kstream_t *ks)					\
	{																\
		if (ks) {													\
			free(ks->buf);											\
			free(ks);												\
		}															\
	}

//ks_getc
// if it's eof, return -1
// if begin (where we're in the buffer) >= end (the end of the buffer)
// read into the buffer again, set begin = 0.  end = how much we've
// read into the buffer
// else, read one char from the buffer.
#define __KS_GETC(__read, __bufsize)						\
	static inline int ks_getc(kstream_t *ks)				\
	{														\
		if (ks->is_eof && ks->begin >= ks->end) return -1;	\
		if (ks->begin >= ks->end) {							\
			ks->begin = 0;									\
			ks->end = __read(ks->f, ks->buf, __bufsize);	\
			if (ks->end < __bufsize) ks->is_eof = 1;		\
			if (ks->end == 0) return -1;					\
		}													\
		return (int)ks->buf[ks->begin++];					\
	}

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m; //l = length.  m = size (available size)
	char *s; //actual string
} kstring_t;
#endif

//kroundup32 rounds x upto nearest 2^n.  000011 => 000100, 001001 => 010000 etc.  (5 => 8.  46 => 64.  515 => 1024.)
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif
//ks_getuntil.  delim = 0 means delim is whitespace).
#define __KS_GETUNTIL(__read, __bufsize)								\
/* read ks until delimieter identified.  delimieter of 0 means look for whitespace.   
copies all chars upto delimieter into str, returns delimiter into dret if dret passed in 
returns the length of the str (ie, chars until delimeter)*/ \
	static int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret) \
	{																	\
		if (dret) *dret = 0;											\
		str->l = 0;														\
		if (ks->begin >= ks->end && ks->is_eof) return -1;	/*end of file*/			\
		for (;;) {														\
			int i;														\
			if (ks->begin >= ks->end) {									\
				if (!ks->is_eof) {										\
					ks->begin = 0;										\
					ks->end = __read(ks->f, ks->buf, __bufsize);		\
					if (ks->end < __bufsize) ks->is_eof = 1;			\
					if (ks->end == 0) break;							\
				} else break;											\
			}															\
			if (delimiter) {											\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (ks->buf[i] == delimiter) break;					\
			} else {													\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (isspace(ks->buf[i])) break;						\
			}															\
	/*allocating space for str->s.  not sure why it needs this much mem*/	\
			if (str->m - str->l < i - ks->begin + 1) {					\
				str->m = str->l + (i - ks->begin) + 1;					\
				kroundup32(str->m);										\
				str->s = (char*)realloc(str->s, str->m);				\
			}															\
/* all of this just to copy that stringinto the end of the name var.	
   memcpy is copying the buffer (from ks_begin to i ) 
   into the end of str->s (str-> l grows as we run this iteratively) */ \
			memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin); \
			str->l = str->l + (i - ks->begin);							\
			ks->begin = i + 1;											\
/* if i < ks->end, the delimeter was found within the buffer */         \
			if (i < ks->end) {											\
				if (dret) *dret = ks->buf[i];							\
				break;													\
			}															\
		}																\
		str->s[str->l] = '\0'; /*sets the terminator char */	 		\
		return str->l;													\
	}

#define KSTREAM_INIT(type_t, __read, __bufsize) \
	__KS_TYPE(type_t)							\
	__KS_BASIC(type_t, __bufsize)				\
	__KS_GETC(__read, __bufsize)				\
	__KS_GETUNTIL(__read, __bufsize)

//inline is a GCC trick to sub in the function wherever it's called.  makes the function effectively a macro.
//type_t is gzFile
#define __KSEQ_BASIC(type_t)											\
	static inline kseq_t *kseq_init(type_t fd)							\
	{																	\
		kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));					\
		s->f = ks_init(fd);												\
		return s;														\
	}																	\
	static inline void kseq_rewind(kseq_t *ks)							\
	{																	\
		ks->last_char = 0;												\
		ks->f->is_eof = ks->f->begin = ks->f->end = 0;					\
	}																	\
	static inline void kseq_destroy(kseq_t *ks)							\
	{																	\
		if (!ks) return;												\
		free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s); \
		ks_destroy(ks->f);												\
		free(ks);														\
	}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
#define __KSEQ_READ									\
/*  kseq_read reads a single sequence from a seq which contains a stream/file that may contain multiple sequences
	seq->last_char is set to 
 */\
	static int kseq_read(kseq_t *seq)									\
	{																	\
		int c;															\
		kstream_t *ks = seq->f;											\
		if (seq->last_char == 0) { /* then jump to the next header line */ \
			while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@');	\
			if (c == -1) return -1; /* end of file */					\
			seq->last_char = c;	/* set last_char to first char in seq */ \
		} /* the first header char has been read */						\
		seq->comment.l = seq->seq.l = seq->qual.l = 0;					\
		if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1; /* don't know when this would return something LESS than 0 */			\
		if (c != '\n') ks_getuntil(ks, '\n', &seq->comment, 0); /* if the delimeter we got wasn't a  newline, explicitly look for newlines and plug the text into seq's comment */			\
		while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') { \
			if (isgraph(c)) { /* printable non-space character */		\
				if (seq->seq.l + 1 >= seq->seq.m) { /* double the memory if we are out of space. */ \
					seq->seq.m = seq->seq.l + 2;						\
					kroundup32(seq->seq.m); /* here's where it actually doubles seq->m (nearest 2^n-s).  kroundup is a macro hence the lack of references etc.*/ \
					seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); \
				}														\
				seq->seq.s[seq->seq.l++] = (char)c;	/* append char, increment length*/ \
			}															\
		}																\
		if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */	\
		seq->seq.s[seq->seq.l] = 0;	/* null terminated string */		\
		if (c != '+') return seq->seq.l; /* FASTA if + not found, we're done.  I think + is for .fq */					\
		if (seq->qual.m < seq->seq.m) {	/* allocate enough memory */	\
			seq->qual.m = seq->seq.m;									\
			seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m);		\
		}																\
		while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */ \
		if (c == -1) return -2; /* we should not stop here */			\
		while ((c = ks_getc(ks)) != -1 && seq->qual.l < seq->seq.l)		\
			if (c >= 33 && c <= 127) seq->qual.s[seq->qual.l++] = (unsigned char)c;	\
		seq->qual.s[seq->qual.l] = 0; /* null terminated string */		\
		seq->last_char = 0;	/* we have not come to the next header line */ \
		if (seq->seq.l != seq->qual.l) return -2; /* qual string is shorter than seq string */ \
		return seq->seq.l; /* no more header lines */												\
	}

#define __KSEQ_TYPE(type_t)						\
	typedef struct {							\
		kstring_t name, comment, seq, qual;		\
		int last_char;							\
		kstream_t *f;							\
	} kseq_t;
//(type_t = gzfile, __read = gzread)
#define KSEQ_INIT(type_t, __read)				\
	KSTREAM_INIT(type_t, __read, 4096)			\
	__KSEQ_TYPE(type_t)							\
	__KSEQ_BASIC(type_t)						\
	__KSEQ_READ

#endif
