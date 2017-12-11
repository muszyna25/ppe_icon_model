#ifndef HAVE_FAST_MATH_LIB
#define EXPON(x,v) v = exp(x)
#define INV_EXPON(x,v) v = exp(-x)
#else
#if (defined HAVE_ACML)
#define EXPON(x,v) call vrda_exp(kproma, x, v)
#define INV_EXPON(x,v) call vrda_exp(kproma, -x, v)
#elif (defined HAVE_MASS)
#define EXPON(x,v) call vrda_exp(kproma, x, v)
#define INV_EXPON(x,v) call vexp(v, -x, kproma)
#elif (defined HAVE_MKL)
#define EXPON(x,v) call vrda_exp(kproma, x, v)
#define INV_EXPON(x,v) call vrda_exp(kproma, -x, v)
#endif
#endif

