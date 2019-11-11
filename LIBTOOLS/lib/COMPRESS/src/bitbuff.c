#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef VPP
# include <sys/types.h>
typedef __uint64_t WORD;
#else
#ifdef SX5
typedef unsigned long uint64_t;
#else
# include <inttypes.h>
#endif
typedef uint64_t WORD;
#endif

#define WORDSIZE 64

#ifdef NO_UNDERSCORE
# define SET_FILLIDX set_fillidx
# define GET_FILLIDX get_fillidx
# define FILL_BBUFF  fill_bbuff
# define SET_EXTRACTIDX set_extractidx
# define GET_EXTRACTIDX get_extractidx
# define EXTRACT_BBUFF  extract_bbuff
#else
# define SET_FILLIDX set_fillidx_
# define GET_FILLIDX get_fillidx_
# define FILL_BBUFF  fill_bbuff_
# define SET_EXTRACTIDX set_extractidx_
# define GET_EXTRACTIDX get_extractidx_
# define EXTRACT_BBUFF  extract_bbuff_
#endif

int outidx = 0;
int outbrem = WORDSIZE ;

int inidx = 0;
int inbrem = WORDSIZE;

void SET_FILLIDX(unsigned *idx, unsigned *bitoffset){
    inidx  = *idx;
    inidx += (*bitoffset/WORDSIZE);
    inbrem = WORDSIZE - (*bitoffset%WORDSIZE);
}

void GET_FILLIDX(unsigned *idx, unsigned *bitoffset){
    *idx    = inidx;
    *bitoffset = WORDSIZE - inbrem;
}

void FILL_BBUFF(WORD *out, int *n, unsigned *val){
    /* inidx = index of the current buffer elt to fill */
    /* inbrem = number of bits remaining on buffer elt out[idx] */

    /* fill buffer out with n low bits of val */

    if (inbrem >= *n){
	inbrem = inbrem - *n;
	/* turn to 0 the n bits of out */
	out[inidx] &= ~(~(~(WORD)0 << *n) << inbrem);
	/* now set the n bits of out to val */
	out[inidx] |= (*val & ~(~(WORD)0 << *n)) << inbrem;
	return;
    } else {
	int nex = *n - inbrem; /* number of bits that will be filled later */
	if (inbrem != 0){
	    /* turn to 0 the inbrem lower bits of out */
	    out[inidx] &= (~(WORD)0 << inbrem) ;
	    /* now set the inbrem lower bits of out with val */
	    out[inidx] |= ((*val >> nex) & ~(~(WORD)0 << inbrem));
	}
	inidx++;
	inbrem = WORDSIZE;
	FILL_BBUFF(out, &nex, val);
    }

}

void SET_EXTRACTIDX(unsigned *idx, unsigned *bitoffset) {
    outidx = *idx;
    outidx += (*bitoffset/WORDSIZE);
    outbrem = WORDSIZE-(*bitoffset%WORDSIZE);
}

void GET_EXTRACTIDX(unsigned *idx, unsigned *bitoffset){
    *idx = outidx;
    *bitoffset = WORDSIZE - outbrem;
}

    
void extract_bbuff_rec(WORD *buff, int *n, unsigned *val) {
    
    if (outbrem >= *n){
	outbrem = outbrem - *n;
	*val = (*val << *n) | (unsigned)((buff[outidx]>>outbrem) & ~(~(WORD)0 << *n));
	return;
    } else {
	int nex = *n - outbrem;
	if (outbrem != 0){
	    *val = (*val << outbrem)| (unsigned)(buff[outidx] & ~(~(WORD)0 << outbrem));

	}
	outidx++;
	outbrem=WORDSIZE;
	extract_bbuff_rec(buff,&nex,val);
    }
}

void EXTRACT_BBUFF(WORD *buff, int *n, unsigned *val) {
    
    unsigned tmpval;

    tmpval=0;
    extract_bbuff_rec(buff,n,&tmpval);
    *val = tmpval;
}

