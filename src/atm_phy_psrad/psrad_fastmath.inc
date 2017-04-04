#ifndef HAVE_FAST_MATH_LIB
#define EXPON(x,v) v = exp(x)
#define INV_EXPON(x,v) v = exp(-x)
#else
#if (defined HAVE_ACML)
#define EXPON(x,v) call vrda_exp(n, x, v)
#define INV_EXPON(x,v) call vrda_exp(n, -x, v)
#elif (defined HAVE_MASS)
#define EXPON(x,v) call vrda_exp(n, x, v)
#define INV_EXPON(x,v) call vexp(v, -x, n)
#elif (defined HAVE_MKL)
#define EXPON(x,v) call vrda_exp(n, x, v)
#define INV_EXPON(x,v) call vdexp(n, -x, v)
#endif
#endif

