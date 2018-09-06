#ifdef __COMPAD__
    REAL(WP)                                   :: goal, dgoal
    !!  ADJOINT BUFFERS  (for ADJ) AND  INPUT COPIES (for ADJ_CP) 
    REAL(WP), ALLOCATABLE, DIMENSION(:,:        ) :: ADJBUF_H,   ADJINP_H
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:      ) :: ADJBUF_VN,  ADJINP_VN
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:,:    ) :: ADJBUF_TR,  ADJINP_TR   
    REAL(WP), ALLOCATABLE, DIMENSION(:,:        ) :: ADJBUF_H0,  ADJINP_H0
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:      ) :: ADJBUF_VN0, ADJINP_VN0
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:,:    ) :: ADJBUF_TR0, ADJINP_TR0
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:      ) :: ADJBUF_GNM1,ADJINP_GNM1
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:      ) :: ADJBUF_DW,  ADJINP_DW
    !!  SHAPES FOR INPUT ALLOCATION
    INTEGER, DIMENSION(2)                      :: szAB_H
    INTEGER, DIMENSION(3)                      :: szAB_VN, szAB_TR, &
         &                                        szAB_GNM1, szAB_DW, szAB_AV
    CHARACTER(len=10)                          :: istr
#endif /* __COMPAD__ */

#if defined(__COMPAD_CPTEST__) || defined(__COMPAD_ADJ_IO_TEST__)
    INTEGER                                    :: nsteps
    CHARACTER(10)                              :: istr1, istr2
#endif

#ifdef __COMPAD_CPTEST__
    INTEGER, DIMENSION(:), ALLOCATABLE         :: iperm
    INTEGER, DIMENSION(3)                      :: iswap
    CHARACTER(100)                             :: tsname
    INTEGER, DIMENSION(:), ALLOCATABLE         :: tscnt
    REAL(wp) :: original_dtime  
#endif /* __COMPAD_CPTEST__ */
