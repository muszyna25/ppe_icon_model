! For runs that check result consistency we fix the different
! contractions that the Intel compiler performs on some loops (at
! least in version 16.0) for the vectorized part and the
! non-vectorized parts
#ifdef IFORT_CONSISTENCY_ENFORCE
#define PREVENT_INCONSISTENT_IFORT_FMA DIR$ NOFMA
#else
#define PREVENT_INCONSISTENT_IFORT_FMA
#endif
