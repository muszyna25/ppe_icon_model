#if defined( _OPENACC) && defined(ACC_LOCAL_FLAG)

#define OPENACC_MODE
#define ACC_PREFIX !$ACC

#else 

#undef  OPENACC_MODE
#define ACC_PREFIX !NOACC

#endif

#undef ACC_LOCAL_FLAG
