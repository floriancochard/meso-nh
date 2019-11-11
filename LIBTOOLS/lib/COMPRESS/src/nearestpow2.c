#include <stdio.h>
#include "ieee754.h"
#include <math.h>

#ifdef NO_UNDERSCORE
# define NEAREST_POW2 nearest_pow2
# define MINBITS_IN_WORD minbits_in_word
# define FMINBITS_IN_WORD fminbits_in_word
#else
# define NEAREST_POW2 nearest_pow2_
# define MINBITS_IN_WORD minbits_in_word_
# define FMINBITS_IN_WORD fminbits_in_word_
#endif

void NEAREST_POW2(union ieee754_double *xval, unsigned int *pow)
{

  if (xval->d != 0.0)
    *pow = xval->ieee.exponent - IEEE754_DOUBLE_BIAS;
  else {
      printf("Warning : NEAREST_POW2 ne traite que des reels > 0.0\n");
      *pow = 0;
  }

}

void MINBITS_IN_WORD(int *nval, unsigned int *nbit)
{
  union ieee754_double xval;
  int ival = *nval;

  /* ne fonctionne qu'avec des entiers non signés */
  if (ival-- < 0){
    printf("Warning : MINBITS_IN_WORD ne traite que des entiers POSITIFS.\n");
    *nbit = -1;
    return;
  } else
    if (ival > 0){
      xval.d = (double)ival;
      NEAREST_POW2(&xval,nbit);
      (*nbit)++;
    } else 
      *nbit = 0 ;
    
}

int FMINBITS_IN_WORD(int *nval)
{
  union ieee754_double xval;
  int ival = *nval;
  unsigned int nbit;

  /* ne fonctionne qu'avec des entiers non signés */
  if (ival < 0){
    printf("Warning : MINBITS_IN_WORD ne traite que des entiers POSITIFS.\n");
    return -1;
  } else {
    if (ival > 0){
      xval.d = (double)ival;
      NEAREST_POW2(&xval,&nbit);
      nbit++;
    } else 
      nbit = 0 ;
    return nbit;
  }
}

/* int main(){ */

/*   double x; */
/*   int i,nbit; */
/*   int exp2; */

/*   printf("Reel : "); */
/*   scanf("%lf",&x); */
  
/*   nearest_pow2_((union ieee754_double*)&x,&exp2); */

/*   printf("2**%d = %lf est la puissance de 2 la plus proche et inferieure à %lf\n", */
/* 	 exp2,pow(2.,exp2),x); */
/*   printf("%lf <= %lf <= %lf\n",pow(2.,(double)exp2),x,pow(2.,(double)exp2+1.)); */
  
/*   printf("Entier positif : "); */
/*   scanf("%d",&i); */
/*   minbits_in_word_(&i,&nbit); */
/*   printf("%d valeurs : %d bits (2**%d = %d).\n",i,nbit,nbit,(1<<nbit)); */
/* } */
