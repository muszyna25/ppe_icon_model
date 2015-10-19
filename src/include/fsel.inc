#ifndef FSEL_INC_INCLUDED
#define FSEL_INC_INCLUDED

#if !(defined __xlC__ && defined _ARCH_PWR6)
#define SWDIV_NOCHK(a,b) ((a)/(b))
#define FSEL(a,b,c) MERGE(b,c,(a).GE.0._wp)
#endif

#endif
