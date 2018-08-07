!
! Cray ftn compilers 8.4 and 8.6 are known to misidentify argument INTENT
! of pointer components, i.e. will disallow changes to an array pointed to
! by a pointer component of a TYPE
! therefore this case needs to be handled specially
!
#if defined _CRAYFTN
#define PTR_INTENT(arg) INTENT(inout)
#else
#define PTR_INTENT(arg) INTENT(arg)
#endif
