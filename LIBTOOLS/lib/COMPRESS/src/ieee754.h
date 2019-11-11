#undef __BYTE_ORDER

#ifdef BIG_endian
# define __BYTE_ORDER 1234
#endif
#ifdef LITTLE_endian
# define __BYTE_ORDER 4321
#endif
#if !(defined(__BYTE_ORDER))
 #error "ieee754.h : you MUST specify \
-DBIG_endian or -DLITTLE_endian \
in CPPFLAGS of your Makefile."
/* Compiler must throw us out at this point! */
#endif

#define __BIG_ENDIAN    1234
#define __LITTLE_ENDIAN 4321

union ieee754_double
  {
    double d;

    /* This is the IEEE 754 double-precision format.  */
    struct
      {
#if     __BYTE_ORDER == __BIG_ENDIAN
        unsigned int negative:1;
        unsigned int exponent:11;
        /* Together these comprise the mantissa.  */
        unsigned int mantissa0:20;
        unsigned int mantissa1:32;
#endif                          /* Big endian.  */
#if     __BYTE_ORDER == __LITTLE_ENDIAN
        /* Together these comprise the mantissa.  */
        unsigned int mantissa1:32;
        unsigned int mantissa0:20;
        unsigned int exponent:11;
        unsigned int negative:1;
#endif                          /* Little endian.  */
      } ieee;

    /* This format makes it easier to see if a NaN is a signalling NaN.  */
    struct
      {
#if     __BYTE_ORDER == __BIG_ENDIAN
        unsigned int negative:1;
        unsigned int exponent:11;
        unsigned int quiet_nan:1;
        /* Together these comprise the mantissa.  */
        unsigned int mantissa0:19;
        unsigned int mantissa1:32;
#else
        /* Together these comprise the mantissa.  */
        unsigned int mantissa1:32;
        unsigned int mantissa0:19;
        unsigned int quiet_nan:1;
        unsigned int exponent:11;
        unsigned int negative:1;
#endif
      } ieee_nan;
  };

#define IEEE754_DOUBLE_BIAS     0x3ff /* Added to exponent.  */
