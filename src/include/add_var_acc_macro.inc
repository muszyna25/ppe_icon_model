#ifdef _OPENACC
    use openacc
    #define __acc_attach(ptr) CALL acc_attach(ptr)
#else
    #define __acc_attach(ptr)
#endif
