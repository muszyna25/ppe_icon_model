! **********************************************************************
MODULE messy_main_channel_bi
! **********************************************************************

  ! -------------------------------------------------------------------
  !  CHANNEL INTERFACE SHARED BETWEEN DIFFERENT BASEMODELS
  ! -------------------------------------------------------------------

#if defined(ECHAM5)
  USE mo_grib,                  ONLY: L_BM_ORIG_OUTPUT
#endif

#ifndef MESSYTIMER 
  USE messy_main_bmluse_bi,     ONLY: &
#else
  USE messy_main_timer_event,   ONLY: &
#endif
                                      time_event, io_time_event  &
                                    , TIME_INC_MONTHS, TRIG_EXACT

  ! MESSy
  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_channel,       ONLY: NMAXCHANNELS, STRLEN_CHANNEL, NCHANNEL &
                                    , modstr, IOMODE_OUT, IOMODE_RST         &
                                    , AMODE_WRITE, AMODE_READ, AMODE_INIT    &
                                    , REPR_UNDEF, DIMID_UNDEF
  USE messy_main_channel_attributes, ONLY: t_attribute_list
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_tools,         ONLY: PTR_4D_ARRAY, PTR_2D_ARRAY

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  PUBLIC :: DIMID_UNDEF, REPR_UNDEF, IOMODE_OUT, IOMODE_RST 

  ! DIMENSION IDs
  ! - TIME
  INTEGER, SAVE, PUBLIC :: DIMID_TIME = DIMID_UNDEF
  ! - GRIDPOINT
  INTEGER, SAVE, PUBLIC :: DIMID_LEV  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_ILEV = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LON  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LAT  = DIMID_UNDEF
  ! - SPECTRAL
  INTEGER, SAVE, PUBLIC :: DIMID_COMPLEX  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_NSP      = DIMID_UNDEF
  ! - FOURIER (ANTI-SYMMETRIC, SYMMETRIC)
  INTEGER, SAVE, PUBLIC :: DIMID_NHGL     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_NMP1     = DIMID_UNDEF
  ! - SPECIAL
  INTEGER, SAVE, PUBLIC :: DIMID_BELOWSF  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_2LEV     = DIMID_UNDEF
  ! - OPTIONAL
  INTEGER, SAVE, PUBLIC :: DIMID_NCELL_ATTILA = DIMID_UNDEF
  !
  ! um_ak_20080425+ COSMO specific
  INTEGER, SAVE, PUBLIC :: DIMID_SRLON  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SRLAT  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_BNDS   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_DAV    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_PRESS  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_ALTIT  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_FIXHT  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_H2m    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_H10m   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_HTOA   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_WBT13  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SOIL   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SOIL1  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SOIL2  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_TLV    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_TKE    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_CANOPY = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_CANOPY1= DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SAT8   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_SAT32  = DIMID_UNDEF
  ! um_ak_20080425-
#ifdef BLANK
  INTEGER, SAVE, PUBLIC :: DIMID_AL     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_BND= REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: ARRAY        = REPR_UNDEF
#endif
  ! mz_ab_20100227+
  ! CMAT specific
  INTEGER, SAVE, PUBLIC :: DIMID_cmat_ht   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_cmat_lat  = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_cmat_lon  = DIMID_UNDEF    
  ! ECHAM+CMAT vertical levels
  INTEGER, SAVE, PUBLIC :: DIMID_cmat_ec_ht    = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_cmat_hc_rates = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_cmat_p_rates  = DIMID_UNDEF
  ! mz_ab_20100227-
  !
  ! op_bk_20130902+
#ifdef __ICON__
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_LEV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_2LEV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_3LEV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_LEVP1_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_NCELLS_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_NVERTS_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_NEDGES_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_NTRACER_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_NTLV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_UBCP_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_UBCP1_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: DIMID_EXTRA2D_DOM(:)
  INTEGER, SAVE, PUBLIC              :: DIMID_ONE
  INTEGER, SAVE, PUBLIC              :: DIMID_TWO
  INTEGER, SAVE, PUBLIC              :: DIMID_FIVE
#endif
  ! op_bk_20130902-
  !
  ! REPRESENTATION IDs
  ! - GRIDPOINT
  INTEGER, SAVE, PUBLIC :: GP_3D_MID        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_INT        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_HORIZONTAL = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_1LEV       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_1D_LEV        = REPR_UNDEF
  ! - SPECTRAL
  INTEGER, SAVE, PUBLIC :: SP_3D_MID        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: SP_3D_INT        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: SP_2D_HORIZONTAL = REPR_UNDEF
  ! - FOURIER
  INTEGER, SAVE, PUBLIC :: FAS_MID          = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: FAS_INT          = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: FAS_MID_ZM       = REPR_UNDEF
  ! - SPECIAL
  INTEGER, SAVE, PUBLIC :: SCALAR           = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_BELOWSF    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_2LEV       = REPR_UNDEF
  ! um_ak_20080425+ COSMO specific
  INTEGER, SAVE, PUBLIC :: GP_4D_MID_TLV    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_INT_TLV    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_SOIL1_TLV  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_SOIL2_TLV  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_MID_BND    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_4D_INT_BND    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_SLON   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_SLAT   = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_P      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MID_Z      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SOIL1      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SOIL2      = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_HORIZ_T    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_HORIZ_BND  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_HORIZ_DAV  = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_HTOTLON    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_TOT        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_LAT_BND    = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_1D_ILEV       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_1D_LAT        = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SAT8       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_SAT32      = REPR_UNDEF
  ! um_ak_20080425-
  ! - OPTIONAL
  INTEGER, SAVE, PUBLIC :: LG_ATTILA        = REPR_UNDEF
  ! op_sb_20100726+
  INTEGER, SAVE, PUBLIC :: LG_FMM_POS       = REPR_UNDEF
  ! op_sb_20100726-
  ! op_bk_20130902+
#ifdef __ICON__
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_HORIZONTAL_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_HOR_2_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_HOR_5_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_HOR_TRACER_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_HOR_UBCP_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_HOR_UBCP1_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_HOR_EXTRA_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_LEV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_2LEV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_3LEV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_LEV_TRACER_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_LEV_TLV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_LEVP1_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_LEVP1_TRACER_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_CELL_LEVP1_TLV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_VERT_HORIZONTAL_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_VERT_LEV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_VERT_LEVP1_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_EDGE_HORIZONTAL_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_EDGE_LEV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_EDGE_2LEV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_EDGE_3LEV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_EDGE_LEV_TLV_DOM(:)
  INTEGER, SAVE, PUBLIC, ALLOCATABLE :: UNSTRUCTURED_EDGE_LEVP1_DOM(:)
#endif
  ! op_bk_20130902-
  ! mz_ap_20070913+
  ! - MPIOM
  INTEGER, SAVE, PUBLIC :: GP_3D_MPIOM           = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_3D_MPIOM_INT       = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: GP_2D_MPIOM           = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LEV_MPIOM       = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LEV_MPIOM_INT   = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_DEPTH_MPIOM     = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_DEPTH_MPIOM_INT = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LON_MPIOM       = DIMID_UNDEF
  INTEGER, SAVE, PUBLIC :: DIMID_LAT_MPIOM       = DIMID_UNDEF
  ! mz_ap_20070913-

  ! mz_ab_20100118+ CMAT2 specific
  INTEGER, SAVE, PUBLIC :: REPRID_cmat3D         = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPRID_cmat3D_BND     = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPRID_ec3D           = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPRID_cmat2D         = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPRID_cmat2D_BND     = REPR_UNDEF
  INTEGER, SAVE, PUBLIC :: REPRID_cmat4D_hcrates = REPR_UNDEF  
  INTEGER, SAVE, PUBLIC :: REPRID_cmat4D_prates  = REPR_UNDEF
  ! mz_ab_20100118-

  ! DECOMPOSITION TYPES
  INTEGER, PARAMETER, PUBLIC :: DC_BC = 0  ! BROADCAST (IO-PE ONLY)
  INTEGER, PARAMETER, PUBLIC :: DC_GP = 1  ! GRIDPOINT
  INTEGER, PARAMETER, PUBLIC :: DC_SP = 2  ! SPECTRAL
  INTEGER, PARAMETER, PUBLIC :: DC_SA = 3  ! FOURIER (SYM. AND ANTISYM.)
  INTEGER, PARAMETER, PUBLIC :: DC_IX = 4  ! INDEX
  ! mz_ap_20070913+
  INTEGER, PARAMETER, PUBLIC :: DC_GP_MPIOM = 5  ! GRIDPOINT MPIOM  
  ! mz_ap_20070913-
  ! um_ak_20090120+
  INTEGER, PARAMETER, PUBLIC :: DC_I2C    = 6  ! GRIDPOINT INT2COSMO
  INTEGER, PARAMETER, PUBLIC :: DC_I2C_IN = 7  ! GRIDPOINT INT2COSMO INPUT
  ! um_ak_20090120-
  INTEGER, PARAMETER, PUBLIC :: DC_CG = 8      ! CMAT GRIDPOINT ! mz_ab_20100210

  ! GP-DECOMPOSITION INFORMATION
#if defined(ECHAM5)
  INTEGER, PARAMETER, PUBLIC             :: gp_nseg = 2
#endif
#ifdef COSMO
  INTEGER, PARAMETER, PUBLIC             :: gp_nseg = 1
  INTEGER                                :: notl        ! number of time levels
#endif
  INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_start => NULL()
  INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_cnt   => NULL()
  INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_meml  => NULL()
  INTEGER, DIMENSION(:,:), POINTER, PUBLIC :: gp_memu  => NULL()

  ! OUTPUT TIME
  TYPE(time_event), DIMENSION(:), ALLOCATABLE :: OUTPUT_EVENT
  LOGICAL,          DIMENSION(:), ALLOCATABLE, PUBLIC :: LOUTPUT_NOW
  LOGICAL, SAVE                               :: LFORCE_NEW_OUTPUT = .FALSE.
  ! SPECIAL INDICES FOR ECHAM5
  ! NOTE: CHANNEL OUTPUT OF tdiag, tdiag_gp, nudg, nudg_gp
  !       WILL BE SYNCHRONIZED TO THE SHORTEST OUTPUT INTERVAL (IF PRESENT !)
  INTEGER, SAVE :: js_tdiag     = 0
  INTEGER, SAVE :: js_tdiag_gp  = 0
  INTEGER, SAVE :: js_nudg      = 0
  INTEGER, SAVE :: js_nudg_gp   = 0

  ! mz_pj_20080905+
  ! TRIGGER NEW FILE
  TYPE(time_event), DIMENSION(:), ALLOCATABLE :: TNF_EVENT
  LOGICAL,          DIMENSION(:), ALLOCATABLE, PUBLIC :: LTNF_NOW
  ! mz_pj_20080905-

  ! SPECIAL INDICES FOR COSMO
  ! NOTE: CHANNEL OUTPUT OF COSMOm, COSMOp and COSMOz
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC :: js_COSMOm
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC :: js_COSMOz
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC :: js_COSMOp
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC :: js_COSMOs
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC :: js_COSMOc
  
  ! SPECIAL RESTART ATTRIBUTES
  TYPE(t_attribute_list), POINTER, SAVE :: restart_att => NULL()

#if defined(ECHAM5) || defined(MBM_CMAT)
  ! <EXP_NAME> + _YYYYMMDD_HHMM_  = 15 + 15
  INTEGER, PARAMETER :: BASENAME_LEN = 30
  ! YYYYMMDD hhmm
  INTEGER, PARAMETER :: DATE_TIME_STR_LEN = 15
#endif
#if (defined COSMO) || defined(BLANK)
  ! <EXP_NAME> + _YYYYMMDD_HHMMSS_  = 15 + 17
  INTEGER, PARAMETER :: BASENAME_LEN = 32
  ! YYYYMMDD hhmmss
  INTEGER, PARAMETER :: DATE_TIME_STR_LEN = 17
#endif
#if defined(__ICON__)
  ! <EXP_NAME> + _YYYYMMDD_HHMMSS_  = 15 + 17
  INTEGER, PARAMETER :: BASENAME_LEN = 32
  ! YYYYMMDD hhmm
  INTEGER, PARAMETER :: DATE_TIME_STR_LEN = 15
#endif

  ! ====================================================================
  ! FOR CPL-NAMELIST 
  ! 
  TYPE t_channel_timer
     CHARACTER(LEN=STRLEN_CHANNEL) :: cname    = ''       ! CHANNEL NAME
     TYPE(io_time_event)          :: io_event = &
          io_time_event(1, TIME_INC_MONTHS,TRIG_EXACT,0) ! DEFAULT
  END TYPE t_channel_timer
  !
  TYPE(t_channel_timer),                          SAVE :: TIMER_DEFAULT
  TYPE(t_channel_timer), DIMENSION(NMAXCHANNELS), SAVE :: TIMER_CHANNEL
  !
  ! mz_pj_20080905+
  TYPE(t_channel_timer),                          SAVE :: TIMER_TNF_DEFAULT
  TYPE(t_channel_timer), DIMENSION(NMAXCHANNELS), SAVE :: TIMER_TNF_CHANNEL
  ! mz_pj_20080905-
  ! ====================================================================

#if defined(ECHAM5)
  ! ====================================================================
  ! SPECIAL HANDLING FOR ECHAM5-STREAM-ELEMENTS WITH laccu = .TRUE.
  TYPE(PTR_4D_ARRAY), DIMENSION(:), POINTER, SAVE :: ptr_accu
  INTEGER,                                   SAVE :: naccu = 0
  ! ====================================================================
#endif

#ifdef COSMO
  ! ====================================================================
  ! SPECIAL HANDLING FOR ACCUMULATED COSMO ELEMENTS 
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE :: ptr_accu
  INTEGER,                                   SAVE :: naccu = 0
  ! ====================================================================
  !
  ! POINTER FIELDS FOR COSMO OUTPUT VARIABLES per CHANNEL
  TYPE cosmo_out_fields
     CHARACTER(LEN=9)                                :: label
     TYPE(PTR_4D_ARRAY), DIMENSION(:), POINTER       :: vars
  END TYPE cosmo_out_fields
  !
  TYPE cosmo_output_list
     TYPE(cosmo_out_fields)           :: this
     TYPE(cosmo_output_list), POINTER :: next => NULL()
  END TYPE cosmo_output_list
  !
  PUBLIC :: cosmo_output_list,  cosmo_out_fields
  !
  TYPE(cosmo_output_list), POINTER, SAVE, PUBLIC :: COSMOOUT => NULL()
  !
  ! FORCE OUTPUT CALCULATIONS in COSMO
  LOGICAL, PUBLIC :: L_FORCE_calcout = .TRUE.
#endif
#if defined(COSMO) || defined(BLANK) || defined(MBM_CMAT)
  ! SWITCH OFF BASEMODEL OUTPUT AND RESTARTS
  LOGICAL, PUBLIC :: L_BM_ORIG_OUTPUT ! um_ak_20100126
#endif
  ! op_bk_20130820+
#if defined(__ICON__)
  ! SWITCH OFF BASEMODEL OUTPUT AND RESTARTS
  LOGICAL, PUBLIC :: L_BM_ORIG_OUTPUT ! um_ak_20100126
#endif
  ! op_bk_20130820-

  PUBLIC   :: main_channel_initialize
  !                         |-> main_channel_read_nml_ctrl (CORE !)
  !                         |-> main_channel_initialize_gatts
  !                         |-> main_channel_initialize_dims
  !                         |-> main_channel_initialize_reprs
  !                         |-> main_channel_initialize_parallel_io
  PUBLIC   :: main_channel_init_memory
  !                         |-> associate_streams_to_channels
  !                         |-> set_COSMO_ORI_attributes
  PUBLIC   :: main_channel_init_coupling
  !                         |-> messy_channel_cosmo_output
  !                                    |-> make_cosmo_channel
  !                         |-> fixate_channels (CORE !)
  !                         |-> main_channel_init_timer
  !                                    |-> main_channel_read_nml_cpl
  !                         |-> write_attribute
  !                         |-> write_dimension
  !                         |-> write_representation
  !                         |-> write_channel
  PUBLIC   :: main_channel_global_start
  !                         |-> set_channel_output (COSMOc)
! op_pj_20110616+
!!$  PUBLIC   :: main_channel_global_end
!!$  !                         |-> main_channel_update_timer
  PUBLIC   :: main_channel_write_output
  !                         |-> main_channel_update_timer
! op_pj_20110616-
  PUBLIC   :: main_channel_free_memory
  !                         |-> channel_finish_io
  !                         |-> clean_representations (CORE !)
  !                         |-> clean_dimensions (CORE !)
  !                         |-> clean_channels (CORE !)
  !                         |-> write_channel
  !
  PUBLIC   :: messy_channel_write_output
  !                         !-> reset_accu_stream_elements(1)
  !                             (reset accu elements to instantaneous values)
  !                         !-> reset_accu_cosmo_elements(1)
  !                             (reset accu elements to instantaneous values)
  !                         |-> update_channels(1) (CORE !)
  !                             ( ACCUMULATE 2ndary DATA )
  !                         |-> update_channels(2) (CORE !)
  !                             ( PREPARE FOR OUTPUT )
  !                         !-> channel_finish_io
  !                         !-> reset_accu_stream_elements(2)
  !                             (set accu elements to zero)
  !                         !-> reset_accu_cosmo_elements(2)
  !                             (set accu elements to zero)
  !                         |-> update_channels(3) (CORE !)
  !                             ( RESET AFTER OUTPUT )
  !                         !-> initialize_restart_attributes
  PUBLIC   :: messy_channel_init_restart
  !                         !-> initialize_parallel_io
  !                         !-> initialize_restart_attributes
  PUBLIC   :: messy_channel_read_restart
  !                         !-> initialize_restart_attributes
  !
  PUBLIC   :: channel_halt
  !
#ifdef COSMO
  PUBLIC   :: messy_COSMO_create_channel
  !
  PUBLIC   :: messy_COSMO_dealloc_meteofields
#endif
  ! SHARED BETWEEN DIFFERENT BASEMODELS
  !
  !PRIVATE :: initialize_restart_attributes  ! INIT RESTART ATTRIBUTES
  !PRIVATE :: main_channel_init_timer        ! INITIALIZE OUTPUT TIMERS
  !PRIVATE :: main_channel_update_timer
  !PRIVATE :: main_channel_read_nml_cpl
  !PRIVATE :: bi_decompose
  !PRIVATE :: bi_vector
  !
  ! BASEMODEL SPECIFIC: SEE messy_main_channel_xx.inc
  !PRIVATE :: main_channel_initialize_gatts  ! SET GLOBAL ATTRIBUTES
  !PRIVATE :: main_channel_initialize_dims   ! DEFINE DIMENSIONS
  !PRIVATE :: main_channel_initialize_reprs  ! DEFINE REPRESENTATIONS
  !
  ! ... ECHAM5
  !PRIVATE :: associate_streams_to_channels
  !PRIVATE :: reset_accu_stream_elements
  !PRIVATE :: timer_sync
  !
  ! ... COSMO
  !PRIVATE :: set_cosmo_channel_attributes
  !PRIVATE :: messy_channel_cosmo_output
  !PRIVATE :: reset_accu_cosmo_elements
  !PRIVATE :: set_COSMO_ORI_attributes
  !PRIVATE :: set_cosmo_attributes

CONTAINS

#if defined(ECHAM5)
  !INCLUDE 'messy_main_channel_e5.inc'
#include "messy_main_channel_e5.inc"
#endif

#ifdef COSMO
  !INCLUDE 'messy_main_channel_c4.inc'
#include "messy_main_channel_c4.inc"
#endif
#ifdef BLANK
#include "messy_main_channel_blank.inc"
#endif
#ifdef MBM_CMAT
#include "messy_main_channel_cm.inc"
#endif
  ! op_bk_20130821+
#ifdef __ICON__
#include "messy_main_channel_icon.inc"
#endif
  ! op_bk_20130821-

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES (MAIN ENTRY POINTS)
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize

    ! BM/MESSy
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_mpi_bi,   ONLY: p_parallel_io, p_io, p_bcast &
                                 , p_pe, p_all_comm
    ! MESSy
    USE messy_main_channel,  ONLY: main_channel_read_nml_ctrl    &
                                 , NMAXCHANNELS, NMAXOBJECTS     &
                                 , NMAXADDCHANNEL, NMAXADDREF    &
                                 , ADD_CHANNEL, ADD_REF, OUT_DEFAULT &
                                 , OUT_PREC & ! mz_pj_20080118
                                 , OUT_CHANNEL, OUT_OBJECT, EXP_NAME &
                                 , L_FLUSH_IOBUFFER &
                                 , I_VERBOSE_LEVEL ! op_pj_20110803
    USE messy_main_tools,    ONLY: find_next_free_unit

#ifdef MPI
    USE messy_main_channel_io, ONLY: initialize_parallel_io
! op_pj_20091124+
#ifdef PNETCDF
    USE messy_main_channel_pnetcdf, ONLY: ch_pnetcdf_read_nml_ctrl &
                                        , NMAX_MPI_IO_HINTS        &
                                        , MPI_IO_HINT
#endif
! op_pj_20091124-
#endif
#ifndef MESSYTIMER
    USE messy_main_bmluse_bi,    ONLY: p_bcast_event
#else
    USE messy_main_timer_bi,     ONLY: p_bcast_event
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize'
    INTEGER :: status
    INTEGER :: iou
    INTEGER :: i

    CALL start_message_bi(modstr,'INITIALIZE CHANNELS',substr)

    ! READ CTRL-NAMELIST FOR FIXATION ( ... AT END OF INIT_COUPLING)
    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_channel_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(EXP_NAME, p_io)
    CALL p_bcast(L_FLUSH_IOBUFFER, p_io)
    CALL p_bcast(I_VERBOSE_LEVEL, p_io)  ! op_pj_20110803
    DO i=1, NMAXADDCHANNEL
       CALL p_bcast(ADD_CHANNEL(i)%cname, p_io)
    END DO
    !
    DO i=1, NMAXADDREF
       CALL p_bcast(ADD_REF(i)%cname1, p_io)
       CALL p_bcast(ADD_REF(i)%oname1, p_io)
       CALL p_bcast(ADD_REF(i)%cname2, p_io)
       CALL p_bcast(ADD_REF(i)%oname2, p_io)
    END DO
    !
    ! op_bk_20130905+
    ! Workaround for gfortran < 4.8.1, wrong handling of namelist input
    ! with more than one "sub-struct"
#ifdef __GFORTRAN__
    CALL p_bcast(OUT_DEFAULT%cname, p_io)
    CALL p_bcast(OUT_DEFAULT%ftype(:), p_io)
    CALL p_bcast(OUT_DEFAULT%ntpf, p_io)
    CALL p_bcast(OUT_DEFAULT%lrestart, p_io)
    CALL p_bcast(OUT_DEFAULT%lignore, p_io)
    CALL p_bcast(OUT_DEFAULT%lout(:), p_io)
#else
    CALL p_bcast(OUT_DEFAULT%cname, p_io)
    CALL p_bcast(OUT_DEFAULT%cio%ftype(:), p_io)
    CALL p_bcast(OUT_DEFAULT%cio%ntpf, p_io)
    CALL p_bcast(OUT_DEFAULT%oio%lrestart, p_io)
    CALL p_bcast(OUT_DEFAULT%oio%lignore, p_io)
    CALL p_bcast(OUT_DEFAULT%oio%lout(:), p_io)
#endif
    ! op_bk_20130905-
    !
    CALL p_bcast(OUT_PREC(:), p_io) ! mz_pj_20080118
    !
    ! op_bk_20130905+
    ! Workaround for gfortran < 4.8.1, wrong handling of namelist input
    ! with more than one "sub-struct"
#ifdef __GFORTRAN__
    DO i=1, NMAXCHANNELS
       CALL p_bcast(OUT_CHANNEL(i)%cname, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%ftype(:), p_io)
       CALL p_bcast(OUT_CHANNEL(i)%ntpf, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%lrestart, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%lignore, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%lout(:), p_io)
    END DO
#else
    DO i=1, NMAXCHANNELS
       CALL p_bcast(OUT_CHANNEL(i)%cname, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%cio%ftype(:), p_io)
       CALL p_bcast(OUT_CHANNEL(i)%cio%ntpf, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%oio%lrestart, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%oio%lignore, p_io)
       CALL p_bcast(OUT_CHANNEL(i)%oio%lout(:), p_io)
    END DO
#endif
    ! op_bk_20130905-
    ! 
    DO i=1, NMAXOBJECTS
       CALL p_bcast(OUT_OBJECT(i)%cname, p_io)
       CALL p_bcast(OUT_OBJECT(i)%oname, p_io)
       CALL p_bcast(OUT_OBJECT(i)%io%lrestart, p_io)
       CALL p_bcast(OUT_OBJECT(i)%io%lignore, p_io)
       CALL p_bcast(OUT_OBJECT(i)%io%lout(:), p_io)
       CALL p_bcast(OUT_OBJECT(i)%io%range(:), p_io)
    END DO

    ! INTIALIZE GLOBAL ATTRIBUTES
    CALL main_channel_initialize_gatts

    ! INTIALIZE DIMENSIONS
    CALL main_channel_initialize_dims

    ! INTIALIZE REPRESENTATIONS
    CALL main_channel_initialize_reprs

    CALL end_message_bi(modstr,'INITIALIZE CHANNELS',substr)

#ifdef MPI
    CALL start_message_bi(modstr,'INITIALIZE PARALLEL I/O',substr)
! op_pj_20091124+
#ifdef PNETCDF
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL ch_pnetcdf_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    DO i=1, NMAX_MPI_IO_HINTS
       CALL p_bcast(MPI_IO_HINT(i)%hint, p_io)
       CALL p_bcast(MPI_IO_HINT(i)%value, p_io)
    END DO
! op_pj_20091124-
#endif
    CALL initialize_parallel_io(status, p_pe, p_io, p_all_comm)
    CALL channel_halt(substr, status)
    CALL end_message_bi(modstr,'INITIALIZE PARALLEL I/O',substr)
#endif

    ! um_ak_20081104+
    ! moved here from init_coupling; otherwise L_BM_ORIG_OUTPUT is not defined
    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL main_channel_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(L_BM_ORIG_OUTPUT, p_io)
    !
    ! um_ak_20081126+
    ! p_bcast_event not possible in initialize,
    ! as the TIMER manager is not yet initialized in case
    ! of RESTART. In this case the MANAGER is initialized
    ! during read_restart
    ! Therefore this part of code was moved back to init_coupling
    ! (main_channel_init_timer)
!!$    CALL p_bcast(TIMER_DEFAULT%cname, p_io)
!!$    CALL p_bcast_event(TIMER_DEFAULT%io_event, p_io)
!!$    DO i=1, NMAXCHANNELS
!!$       CALL p_bcast(TIMER_CHANNEL(i)%cname, p_io)
!!$       CALL p_bcast_event(TIMER_CHANNEL(i)%io_event, p_io)
!!$    END DO
!!$
!!$    ! mz_pj_20080905+
!!$    CALL p_bcast(TIMER_TNF_DEFAULT%cname, p_io)
!!$    CALL p_bcast_event(TIMER_TNF_DEFAULT%io_event, p_io)
!!$    DO i=1, NMAXCHANNELS
!!$       CALL p_bcast(TIMER_TNF_CHANNEL(i)%cname, p_io)
!!$       CALL p_bcast_event(TIMER_TNF_CHANNEL(i)%io_event, p_io)
!!$    END DO    
    ! um_ak_20081126-
    ! mz_pj_20080905-
    ! um_ak_20081104-

  END SUBROUTINE main_channel_initialize
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_init_memory

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_init_memory'

    CALL start_message_bi(modstr,'INITIALIZE CHANNEL MEMORY',substr)

#if defined(ECHAM5)
    ! ASSOCIATE ECHAM5 STREAMS TO MESSy CHANNELS
    CALL associate_streams_to_channels
#endif

#ifdef COSMO
    CALL set_COSMO_ORI_attributes
#endif

    ! op_bk_20130903+
#ifdef __ICON__
    CALL associate_var_lists_to_channels
#endif
    ! op_bk_20130903-

    CALL end_message_bi(modstr,'INITIALIZE CHANNEL MEMORY',substr)

  END SUBROUTINE main_channel_init_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_init_coupling

    ! BM/MESSy
    USE messy_main_mpi_bi,   ONLY: p_parallel_io, p_bcast
#ifdef COSMO
    USE messy_main_data_bi,  ONLY: ngribout
#endif
    ! MESSy
    USE messy_main_tools,    ONLY: int2str
    ! CHANNEL
    USE messy_main_channel,  ONLY: fixate_channels, write_channel &
                                 , write_attribute, set_channel_output
    USE messy_main_channel_dimensions,  ONLY: write_dimension
    USE messy_main_channel_repr,        ONLY: write_representation

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_init_coupling'
    INTEGER                     :: status

#ifdef COSMO
    CALL start_message_bi(modstr,'COSMO OUTPUT TO CHANNELS',substr)
    ! ASSOCIATE COSMO OUTPUT WITH CHANNELS
    ! MUST BE CALLED BEFORE FIXATE CHANNELS
    CALL messy_channel_cosmo_output
    CALL end_message_bi(modstr,'COSMO OUTPUT TO CHANNELS',substr)
#endif

    CALL start_message_bi(modstr,'FIXATE CHANNELS',substr)

    CALL fixate_channels(status)
    CALL channel_halt(substr, status)

    ! INITIALIZE CHANNEL OUTPUT TIMERS (via CPL-NAMELIST)
    CALL main_channel_init_timer

    IF (p_parallel_io) THEN
       CALL write_attribute(status)
       CALL channel_halt(substr, status)
    END IF

    IF (p_parallel_io) THEN
       CALL write_dimension(status)
       CALL channel_halt(substr, status)
    END IF

    IF (p_parallel_io) THEN
       CALL write_representation(status)
       CALL channel_halt(substr, status)
    END IF

    IF (p_parallel_io) THEN
       CALL write_channel(status)
       CALL channel_halt(substr, status)
    END IF

    CALL end_message_bi(modstr,'FIXATE CHANNELS',substr)

  END SUBROUTINE main_channel_init_coupling
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_global_start
    
#ifdef COSMO
    USE messy_main_data_bi,  ONLY: ngribout
    ! MESSy
    USE messy_main_tools,    ONLY: int2str
    USE messy_main_timer,    ONLY: lfirst_cycle
    USE messy_main_channel,  ONLY: set_channel_output
    IMPLICIT NONE
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='main_channel_global_start'
    CHARACTER(LEN=3) :: str
    INTEGER          :: iCc, status

    IF (lfirst_cycle) THEN
       ! force output of channels with constants at the beginning
       DO iCc = 1, ngribout
          IF (js_COSMOc(iCc) > 0) THEN
             CALL int2str(str,iCc)
             CALL set_channel_output(status,'COSMOc'//str, .TRUE.)
             CALL channel_halt(substr,status)
          ENDIF
       ENDDO
    ENDIF
#endif

  END SUBROUTINE main_channel_global_start
  ! -------------------------------------------------------------------

  ! op_pj_20110616+
!!$  ! -------------------------------------------------------------------
!!$  SUBROUTINE main_channel_global_end
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! LOCAL
!!$    !CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_global_end'
!!$    
!!$     CALL main_channel_update_timer   
!!$
!!$  END SUBROUTINE main_channel_global_end
!!$  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_write_output

    IMPLICIT NONE

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_write_output'
    
     CALL main_channel_update_timer   

  END SUBROUTINE main_channel_write_output
  ! -------------------------------------------------------------------  
  ! op_pj_20110616-

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_free_memory

    ! BM/MESSy
    USE messy_main_mpi_bi,   ONLY: p_parallel_io

    ! MESSy
    USE messy_main_channel_io,          ONLY: channel_finish_io
    USE messy_main_channel_dimensions,  ONLY: clean_dimensions
    USE messy_main_channel_repr,        ONLY: clean_representations
    USE messy_main_channel,             ONLY: clean_channels, write_channel

    IMPLICIT NONE

    INTRINSIC :: ALLOCATED

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_free_memory'
    INTEGER :: status

    CALL start_message_bi(modstr,'FINISH CHANNELS',substr)

    ! ------------------------------------------
    ! CLOSE ALL OUTPUT FILES
    ! ------------------------------------------
    CALL channel_finish_io(status, p_parallel_io, IOMODE_OUT, .TRUE.)
    CALL channel_halt(substr, status)

    CALL clean_representations(status)
    CALL channel_halt(substr, status)

    CALL clean_dimensions(status)
    CALL channel_halt(substr, status)

    CALL clean_channels(status)
    CALL channel_halt(substr, status)

    IF (p_parallel_io) THEN
       CALL write_channel(status)
       CALL channel_halt(substr, status)
    END IF

#ifdef COSMO
    ! CLEAN COSMO output channel indices
    IF (ALLOCATED(js_COSMOm)) DEALLOCATE(js_COSMOm)
    IF (ALLOCATED(js_COSMOp)) DEALLOCATE(js_COSMOp)
    IF (ALLOCATED(js_COSMOz)) DEALLOCATE(js_COSMOz)
    IF (ALLOCATED(js_COSMOs)) DEALLOCATE(js_COSMOs)
    IF (ALLOCATED(js_COSMOc)) DEALLOCATE(js_COSMOc)
#endif

    CALL end_message_bi(modstr,'FINISH CHANNELS',substr)

  END SUBROUTINE main_channel_free_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PUBLIC, BUT NOT VIA ENTRY POINT IN messy_main_control_xx
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE messy_channel_write_output(IOMODE)

    ! called in stepon.f90 in case of ECHAM5

    ! BM/MESSy
    USE messy_main_mpi_bi,   ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
#ifndef MESSYTIMER
    USE messy_main_data_bi,    ONLY: &
#else
    USE messy_main_timer,      ONLY: &
#endif
                                     YEAR, MONTH, DAY, HOUR, MINUTE, SECOND &
                                   , YEAR_NEXT, MONTH_NEXT, DAY_NEXT     &
                                   , HOUR_NEXT, MINUTE_NEXT, SECOND_NEXT &
                                   , YEAR_START, MONTH_START, DAY_START     &
                                   , HOUR_START, MINUTE_START, SECOND_START &
                                   , delta_time, current_time_step &
                                   , time_step_len
    USE messy_main_data_bi,    ONLY: l2tls

    ! MESSy
    USE messy_main_channel_io, ONLY: channel_init_io           &
                                   , channel_write_header      &
                                   , channel_write_time        &
                                   , channel_write_data        &
                                   , channel_finish_io
    USE messy_main_channel,    ONLY: EXP_NAME, update_channels, I_VERBOSE_LEVEL
    USE messy_main_channel_dimensions, ONLY: update_dimension_variable
    USE messy_main_timer,      ONLY: time_span_d  ! mz_pj_20090519

    IMPLICIT NONE

    INTRINSIC :: ABS, ASSOCIATED, LEN, LEN_TRIM, NULL, REAL

    ! I/O
    INTEGER, INTENT(IN)         :: IOMODE  ! OUTPUT MODE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'messy_channel_write_output'
    INTEGER                     :: status
    CHARACTER(LEN=BASENAME_LEN) :: fnamebase = ''
    INTEGER                     :: i
    REAL(DP)                    :: yyyymmdd
    REAL(DP)                    :: now
    LOGICAL                     :: lexit
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: ptr  => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: gptr => NULL()
    INTEGER                               :: reprid
    INTEGER, SAVE    :: nrstcount     = 0
    CHARACTER(LEN=4) :: nrstcount_str = ''
    LOGICAL          :: lp
!    INTEGER          :: dts ! mz_pj_20080327 ! mz_pj_20090519
    INTEGER          :: YEAR_DATE, MONTH_DATE, DAY_DATE     ! um_ak_20091110
    INTEGER          :: HOUR_DATE, MINUTE_DATE, SECOND_DATE ! um_ak_20091110

    IF (I_VERBOSE_LEVEL >= 2) &
         CALL start_message_bi(modstr,'WRITE OUTPUT',substr)

    SELECT CASE(IOMODE)
       !
    CASE(IOMODE_OUT)   ! ### ---------------- OUTPUT -------------------- ###
       !
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! SPECIAL HANDLING FOR ECHAM5-STREAM-ELEMENTS WITH laccu = .TRUE.
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(ECHAM5)
       CALL reset_accu_stream_elements(1)
#endif
#ifdef COSMO
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! SPECIAL HANDLING FOR ACCUMULATED COSMO FIELDS 
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL reset_accu_cosmo_elements(1)
#endif
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       CALL update_channels(status, 1, time_step_len) ! ACCUMULATE 2ndary
       CALL channel_halt(substr, status)
       !
       CALL update_channels(status, 2, time_step_len) ! PREPARE FOR OUTPUT
       CALL channel_halt(substr, status)
       !
       ! UPDATE FILENAME BASE
       fnamebase = EXP_NAME
       DO i = LEN_TRIM(EXP_NAME)+1, LEN(EXP_NAME)
          WRITE(fnamebase(i:i),'(a1)') '_' 
       END DO

       ! um_ak_20091110+
!!$#ifdef MESSYTIMER         ! op_pj_20100212
!!$       IF (l2tls) THEN       ! op_pj_20100212
!!$          YEAR_DATE   = YEAR
!!$          MONTH_DATE  = MONTH
!!$          DAY_DATE    = DAY
!!$          HOUR_DATE   = HOUR
!!$          MINUTE_DATE = MINUTE
!!$          SECOND_DATE = SECOND
!!$#else                     ! op_pj_20100212
!!$       ELSE                  ! op_pj_20100212
          YEAR_DATE   = YEAR_NEXT
          MONTH_DATE  = MONTH_NEXT
          DAY_DATE    = DAY_NEXT
          HOUR_DATE   = HOUR_NEXT
          MINUTE_DATE = MINUTE_NEXT
          SECOND_DATE = SECOND_NEXT
!!$       END IF               ! op_pj_20100212
!!$#endif                   ! op_pj_20100212
       ! um_ak_20091110-
#if defined(ECHAM5) || defined(MBM_CMAT)
       WRITE(fnamebase(LEN(EXP_NAME)+1:),'(i4,i2.2,i2.2,a1,i2.2,i2.2,a1)') &
! um_ak_20091117+
!!$            YEAR_NEXT, MONTH_NEXT, DAY_NEXT, '_', HOUR_NEXT, MINUTE_NEXT, '_'
            YEAR_DATE, MONTH_DATE, DAY_DATE, '_', HOUR_DATE, MINUTE_DATE, '_'
! um_ak_20091117-
#endif
#if defined(COSMO) || defined(BLANK) || defined (__ICON__)
       WRITE(*,*) "MESSY: main_channel_write_output"
       WRITE(*  &
            , '(i4,i2.2,i2.2,a1,i2.2,i2.2,i2.2,a1)')&
            YEAR_DATE, MONTH_DATE,  DAY_DATE, '_' &
            , HOUR_DATE, MINUTE_DATE, SECOND_DATE,'_'
       WRITE(fnamebase(LEN(EXP_NAME)+1:)  &
            , '(i4.4,i2.2,i2.2,a1,i2.2,i2.2,i2.2,a1)')&
! um_ak_20091117+
!!$              YEAR_NEXT, MONTH_NEXT,  DAY_NEXT, '_' &
!!$            , HOUR_NEXT, MINUTE_NEXT, SECOND_NEXT,'_'
            YEAR_DATE, MONTH_DATE,  DAY_DATE, '_' &
            , HOUR_DATE, MINUTE_DATE, SECOND_DATE,'_'
! um_ak_20091117-
#endif
      !
       ! UPDATE TIME STEP - DATA
       ! - TIME
! mz_pj_20090519+
       CALL time_span_d(now   &
            , YEAR_START, MONTH_START, DAY_START &
            , HOUR_START, MINUTE_START, SECOND_START  &
! um_ak_20091110 , YEAR_NEXT, MONTH_NEXT, DAY_NEXT         &
! um_ak_20091110 , HOUR_NEXT, MINUTE_NEXT, SECOND_NEXT)
            , YEAR_DATE, MONTH_DATE, DAY_DATE         &
            , HOUR_DATE, MINUTE_DATE, SECOND_DATE)
! mz_pj_20090519-

       CALL update_dimension_variable(status, 'time', 'time', (/ now /))
       CALL channel_halt(substr, status)
       ! - YYYYMMDD
       yyyymmdd = REAL(ABS(YEAR)*10000 + MONTH*100 +DAY, DP)        &
            + REAL((HOUR*3600 + MINUTE*60 + SECOND), DP)/86400.0_DP
       IF (YEAR<0) yyyymmdd = -yyyymmdd
       CALL update_dimension_variable(status, 'time', 'YYYYMMDD' &
            , (/ yyyymmdd /))
       CALL channel_halt(substr, status)
       ! - DT
       CALL update_dimension_variable(status, 'time', 'dt', (/ delta_time /))
       CALL channel_halt(substr, status)
       ! - TIME STEP
       CALL update_dimension_variable(status, 'time', 'nstep' &
            , (/ REAL(current_time_step, DP) /))
       CALL channel_halt(substr, status)
       !
    CASE(IOMODE_RST)   ! ### ---------------- RESTART -------------------- ###
       !
       ! FORCE NEW OUTPUT FILE IN NEXT STEP
       LFORCE_NEW_OUTPUT = .TRUE.
       ! UPDATE FILENAME BASE
       nrstcount = nrstcount + 1
       write(nrstcount_str,'(i4.4)') nrstcount
       ! UPDATE FILENAME BASE
       fnamebase = 'restart_'//nrstcount_str//'_'
       !
       ! SET RESTART ATTRIBUTES
       CALL initialize_restart_attributes(AMODE_WRITE)
       !
    END SELECT

    ! PREPARE OUTPUT / RESTART FILE
    ! NEW FILE, SAVE I/O-UNITS, FILE-IDs etc.
    CALL channel_init_io(status, p_parallel_io, IOMODE, fnamebase, AMODE_WRITE)
    CALL channel_halt(substr, status)
    !
    ! WRITE HEADER
    SELECT CASE(IOMODE)
    CASE(IOMODE_OUT)
       WRITE(*,*) "MESSY: Write Header of NC Files"
       CALL channel_write_header(status, p_parallel_io, IOMODE, DIMID_TIME)
       WRITE(*,'(A37, I4)') "MESSY: Write Header of NC Files DONE", status
    CASE(IOMODE_RST)
       CALL channel_write_header(status, p_parallel_io, IOMODE, DIMID_TIME &
            , restart_att)
    END SELECT
    CALL channel_halt(substr, status)
    !
    ! WRITE TIME STEP DATA TO OUTPUT FILE
    WRITE(*,*) "MESSY: Write TIME of NC Files"
    CALL channel_write_time(status, p_parallel_io, IOMODE, DIMID_TIME)
    WRITE(*,'(A37, I4)') "MESSY: Write TIME of NC Files DONE", status
    CALL channel_halt(substr, status)

    ! WRITE DATA
    DO
       WRITE(*,*) "MESSY: Write DATA"
       ! NEXT POINTER
       CALL channel_write_data(status, lp, p_parallel_io &
            , IOMODE, lexit, ptr, reprid)
       WRITE(*,'(A37, I4)') "MESSY: Write DATA DONE", status
       CALL channel_halt(substr, status)
       ! NOTE: In case of non-parallel I/O, lexit is only set correctly
       !       on the I/O PE ...
       IF (.NOT. lp) CALL p_bcast(lexit, p_io)
       IF (lexit) EXIT
       IF (lp) THEN ! parallel I/O
          CALL bi_vector(status, -1, reprid, gptr, ptr)
          IF (status /= 0) &
               CALL error_bi('bi_vector reported an error', substr)
       ELSE
          CALL bi_decompose(status, -1, reprid, gptr, ptr)
          IF (status /= 0) &
               CALL error_bi('bi_decompose reported an error',substr)
       END IF
       ! OUTPUT
       ! NOTE: gptr only associated on I/O-PE
       CALL channel_write_data(status, lp, p_parallel_io &
            , IOMODE, lexit, gptr, reprid)
       CALL channel_halt(substr, status)
    END DO

    ! FLUSH (OUTPUT) / CLOSE (RESTART) BUFFER
    CALL channel_finish_io(status, p_parallel_io, IOMODE, (IOMODE==IOMODE_RST))
    CALL channel_halt(substr, status)

    ! CLEAN MEMORY
    IF (ASSOCIATED(gptr)) THEN 
       DEALLOCATE(gptr)
       NULLIFY(gptr)
    END IF

    IF (IOMODE == IOMODE_OUT) THEN
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! SPECIAL HANDLING FOR ECHAM5-STREAM-ELEMENTS WITH laccu = .TRUE.
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(ECHAM5)
       CALL reset_accu_stream_elements(2)
#endif
#ifdef COSMO
       CALL reset_accu_cosmo_elements(2)
#endif
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       CALL update_channels(status, 3, time_step_len) ! RESET AFTER OUTPUT
       CALL channel_halt(substr, status)
       !
    END IF

    IF (I_VERBOSE_LEVEL >= 2) &
         CALL end_message_bi(modstr,'WRITE OUTPUT',substr)

  END SUBROUTINE messy_channel_write_output
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE messy_channel_init_restart

    ! BM/MESSy
    USE messy_main_mpi_bi,  ONLY: p_io, p_pe, p_all_comm
#ifndef MESSYTIMER
    USE messy_main_data_bi, ONLY: lresume
#else
    USE messy_main_timer,   ONLY: lresume
#endif
    ! MESSy
    USE messy_main_channel_io, ONLY: initialize_parallel_io

    ! called in initialize.f90

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'messy_channel_init_restart'
    INTEGER :: status

    IF (.NOT. lresume) RETURN

    CALL initialize_parallel_io(status, p_pe, p_io, p_all_comm)
    CALL channel_halt(substr, status)

    CALL initialize_restart_attributes(AMODE_INIT)
    
  END SUBROUTINE messy_channel_init_restart
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
! um_ak_20110512+
!!$  SUBROUTINE messy_channel_read_restart
  SUBROUTINE messy_channel_read_restart(chname)
! um_ak_20110512-

    ! called in iorestart.f90 in case of ECHAM5

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi

    USE messy_main_channel_io, ONLY: channel_init_io   &
                                   , channel_read_data &
                                   , channel_finish_io
    USE messy_main_channel,    ONLY: t_channel_list, t_channel, GCHANNELLIST

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL

    !um_ak_20110512+
    ! I/O    
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: chname 
    !um_ak_20110512-

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'messy_channel_read_restart'
    INTEGER                      :: status
    CHARACTER(LEN=BASENAME_LEN)  :: fnamebase = ''
    TYPE(t_channel_list),  POINTER        :: ls
    TYPE(t_channel),       POINTER        :: channel
    LOGICAL                               :: lexit
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: ptr  => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: gptr => NULL()
    INTEGER                               :: reprid
    LOGICAL                               :: lok
    LOGICAL                               :: lp

    CALL start_message_bi(modstr,'READ RESTART',substr)

    !um_ak_20110512+
    IF (PRESENT(chname)) THEN
       IF (p_parallel_io) WRITE(*,*) 'chname: ', TRIM(chname)
    ENDIF
    !um_ak_20110512-
    ! UPDATE FILENAME BASE
    fnamebase = 'restart_'

    ! INITIALIZE RESTART ATTRIBUTES
    CALL initialize_restart_attributes(AMODE_READ)

    ! OPEN RESTART FILES AND CHECK HEADER INFORMATION
    ! OPEN FILE FOR READ
    !um_ak_20110512+
!!$    CALL channel_init_io(status, p_parallel_io &
!!$         , IOMODE_RST, fnamebase, AMODE_READ, restart_att)
    CALL channel_init_io(status, p_parallel_io &
         , IOMODE_RST, fnamebase, AMODE_READ, restart_att, chname)
    !um_ak_20110512-
    CALL channel_halt(substr, status)

    ! BROADCAST %tslo IN ALL CHANNELS
    ! this should only be required for for non-parallel I/O, however ...
    ls => GCHANNELLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(ls)) EXIT
       channel => ls%this
       ! ------------------------------------------------
       CALL p_bcast(channel%int%tslo, p_io)
       ! ------------------------------------------------
       ls => ls%next
    END DO channel_loop

    ! READ DATA
    DO
       ! NEXT POINTER
       !um_ak_20110512+
!!$       CALL channel_read_data(status, p_parallel_io &
!!$            , IOMODE_RST, lexit, gptr, reprid, lp)
       CALL channel_read_data(status, p_parallel_io &
            , IOMODE_RST, lexit, gptr, reprid, lp, chname)
       !um_ak_20110512-
       CALL channel_halt(substr, status)
       ! NOTE: In case of non-parallel I/O, lexit is only set correctly
       !       on the I/O PE ...
       IF (.NOT. lp) CALL p_bcast(lexit, p_io)
       IF (lexit) EXIT
       !
       IF (lp) THEN ! parallel I/O
          IF (ASSOCIATED(gptr)) THEN
             CALL bi_vector(status, 1, reprid, gptr, ptr)
             IF (status /= 0) CALL error_bi( &
                  'bi_vector reported an error', substr)
          END IF
       ELSE
          ! NOTE: gptr only associated on I/O-PE
          IF (p_parallel_io) lok = ASSOCIATED(gptr)
          CALL p_bcast(lok, p_io)
          !
          IF (lok) THEN      
             CALL bi_decompose(status, 1, reprid, gptr, ptr)
             IF (status /= 0) CALL error_bi( &
                  'bi_decompose reported an error',substr)
             ! NOTE: ptr now associated on all PEs
             !       or =>NULL if .NOT. lok
          END IF
       END IF
       !
       ! DISTRIBUTE
       CALL channel_read_data(status, p_parallel_io, &
            IOMODE_RST, lexit, ptr, reprid, lp)
       CALL channel_halt(substr, status)
       !
       ! RESET MEMORY
       IF (ASSOCIATED(ptr)) THEN 
          DEALLOCATE(ptr)
          NULLIFY(ptr)
       END IF
       IF (ASSOCIATED(gptr)) THEN 
          DEALLOCATE(gptr)
          NULLIFY(gptr)
       END IF
    END DO

    ! CLOSE ALL RESTART FILES
    !um_ak_20110512+
!!$    CALL channel_finish_io(status, p_parallel_io, IOMODE_RST, .TRUE.)
    CALL channel_finish_io(status, p_parallel_io, IOMODE_RST, .TRUE., chname)
    !um_ak_20110512-
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr,'READ RESTART',substr)

  END SUBROUTINE messy_channel_read_restart
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PUBLIC HELPER ROUTINES
  ! -------------------------------------------------------------------
  SUBROUTINE channel_halt(substr, status)

    ! BM/MESSy
    USE messy_main_blather_bi,     ONLY: error_bi, info_bi
    ! MESSy
    USE messy_main_constants_mem,  ONLY: STRLEN_VLONG
    USE messy_main_channel_error,  ONLY: channel_error_str

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: substr
    INTEGER,          INTENT(IN)  :: status
    ! LOCAL
    CHARACTER(LEN=STRLEN_VLONG)   :: errstr

    IF (status == 0) RETURN

    errstr = channel_error_str(status)

    CALL error_bi(errstr,substr)

  END SUBROUTINE channel_halt
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE initialize_restart_attributes(AMODE)

#if defined(ECHAM5)
    USE messy_main_bmluse_bi,  ONLY: ec_manager_init                 &
                                   , TC_convert, TC_set, time_intern
#ifndef MESSYTIMER
    ! ECHAM5/MESSy
    USE messy_main_bmluse_bi,  ONLY: next_date, get_date_components      &
                                   , resume_date, start_date             &
                                   , time_days, get_time_step, init_step 
    USE messy_main_data_bi,    ONLY: delta_time
#else
    USE messy_main_bmluse_bi,  ONLY: E5_resume_date => resume_date &
                                   , E5_start_date  => start_date  &
                                   , lresume
    USE messy_main_timer,      ONLY: timer_get_lresume
#endif
#endif

    ! BM/MESSY
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
#ifdef COSMO
    USE messy_main_data_bi,    ONLY: ivctype, nnow, nnew, nold, dt &
                                   , vcflat, p0sl, nfltvc, irefatm &
                                   , hnextrad
#endif

    ! MESSy
#ifdef MESSYTIMER
    USE messy_main_timer_bi,   ONLY: get_time_step => timer_get_time_step &
                                   , messy_timer_init_manager             &
#ifdef COSMO
                                   , messy_timer_COSMO_reinit_time &
#endif
                                   , timer_message                        
    USE messy_main_timer,      ONLY: time_days, INIT_STEP                  &
                                   , resume_date, start_date, next_date    &
                                   , timer_get_date,timer_set_date         &
                                   , delta_time                            &
                                   , YEAR_START,MONTH_START,DAY_START      &
                                   , HOUR_START, MINUTE_START, SECOND_START
#endif    
    USE messy_main_channel_io, ONLY: channel_init_restart
    USE messy_main_channel,    ONLY: new_attribute, get_attribute &
                                   , AF_RST_CMP, AF_RST_INP, AF_RST_NONE
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

    IMPLICIT NONE
    INTRINSIC :: INT, MOD, SIGN, TRIM

    ! I/O
    INTEGER, INTENT(IN) :: AMODE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='initialize_restart_attributes'
    CHARACTER(LEN=DATE_TIME_STR_LEN) :: str_date_time = ''  
    !
    INTEGER           :: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
    INTEGER           :: iflag
    INTEGER           :: iflag_dt
    INTEGER           :: yyyymmdd, hhmmss
#if defined(ECHAM5)
    TYPE(time_intern) :: io_date
#endif
    INTEGER, SAVE     :: nstep
    INTEGER           :: timestep   
    INTEGER           :: status
#ifndef MESSYTIMER
    TYPE(time_days)   :: zdate
#endif
    LOGICAL           :: lp = .FALSE. ! um_ak_20100127
    INTEGER           :: kyr, kmo, kdy, khr, kmn, kse, iymd
    CHARACTER(LEN=STRLEN_MEDIUM) :: zdstr = ''

#if defined(ECHAM5)
#ifdef MESSYTIMER
    ! set lresume of original ECHAM5 timer
    CALL timer_get_lresume(lresume)
#endif
#endif
    SELECT CASE(AMODE)
       CASE(AMODE_WRITE)
          iflag    = AF_RST_NONE
          iflag_dt = AF_RST_NONE
#ifndef MESSYTIMER
          zdate    = next_date
#endif
          zdstr    = 'next'
          nstep    = get_time_step()
       CASE(AMODE_READ)
          iflag    = AF_RST_CMP
          iflag_dt = AF_RST_NONE
#ifndef MESSYTIMER
          zdate    = resume_date
#endif
          zdstr    = 'resume'
       CASE(AMODE_INIT)
          iflag    = AF_RST_INP
          iflag_dt = AF_RST_INP
#ifndef MESSYTIMER
          !um_ak_20090625+
          !use start-date here as restart date ist not yet defined (for TIMER)
          ! an will be overwritten anyway
          zdate    = start_date 
          !zdate    = resume_date
          !um_ak_20090625-
#endif
          zdstr    = 'resume'
          nstep    = INIT_STEP
    END SELECT

#ifndef MESSYTIMER
    ! - START DATE AND TIME
    CALL get_date_components(start_date     &
         ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)
    WRITE(str_date_time,'(i4,i2.2,i2.2,a1,3(i2.2))') &
         YEAR, MONTH, DAY,' ',HOUR, MINUTE, SECOND
#else
    ! NOTE: The time information in *_START for the MESSY-TIMER must be
    !       set here, therfore the local variables (YEAR, MONTH, ...)
    !       cannot be used.
    WRITE(str_date_time,'(i4,i2.2,i2.2,a1,3(i2.2))') &
          YEAR_START, MONTH_START,  DAY_START,' '    &
          ,HOUR_START, MINUTE_START, SECOND_START
#endif
    !
    CALL new_attribute(status, restart_att &
         , 'start_date_time', c=TRIM(str_date_time) &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)

    ! - RESTART DATE AND TIME
#ifndef MESSYTIMER
    CALL get_date_components(zdate &
         ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)
#else
    CALL timer_get_date(status, TRIM(zdstr)  &
         ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)
    CALL timer_message(status, substr)
#endif
    WRITE(str_date_time,'(i4,i2.2,i2.2,a1,3(i2.2))') &
         YEAR, MONTH, DAY,' ',HOUR, MINUTE, SECOND
    !
    CALL new_attribute(status, restart_att &
         , 'restart_date_time', c=TRIM(str_date_time) &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)

    ! - CURRENT TIME STEP
    CALL new_attribute(status, restart_att &
         , 'nstep', i=nstep                &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)

    ! - TIME STEP LENGTH
    CALL new_attribute(status, restart_att &
         , 'timestep', i=INT(delta_time)   &
         , loverwrite=.TRUE., iflag = iflag_dt)
    CALL channel_halt(substr, status)
    
#ifdef COSMO
    ! - counter for radiation calculation
    CALL new_attribute(status, restart_att &
         , 'hnextrad', r=hnextrad   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    ! - COORDINATE TYPE
    CALL new_attribute(status, restart_att &
         , 'irefatm', i=irefatm   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    ! - COORDINATE TYPE
    CALL new_attribute(status, restart_att &
         , 'ivctype', i=ivctype   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    ! - current index of time_levels
    CALL new_attribute(status, restart_att &
         , 'nnow', i=nnow   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, restart_att &
         , 'nnew', i=nnew   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, restart_att &
         , 'nold', i=nold   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, restart_att &
         , 'nfltvc', i=nfltvc   &
         , loverwrite=.TRUE., iflag = iflag)
    CALL channel_halt(substr, status)
#endif   
 
    ! CONTINUE ONLY DURING INITIALIZATION AFTER RESTART
    IF (AMODE /= AMODE_INIT) RETURN

#if defined(ECHAM5)
    CALL channel_init_restart(status, lp, p_parallel_io &
         , 'restart_g1a', restart_att)
    CALL channel_halt(substr, status)
#endif

#ifdef COSMO
    CALL channel_init_restart(status, lp, p_parallel_io &
         , 'restart_COSMO_ORI', restart_att)
    CALL channel_halt(substr, status)
#endif

#ifdef BLANK
    CALL channel_init_restart(status, lp, .TRUE. &
         , 'restart_BLANK', restart_att)
    CALL channel_halt(substr, status)
#endif

#ifdef MBM_CMAT
    CALL channel_init_restart(status, lp, p_parallel_io &
         , 'restart_cmat', restart_att)
    CALL channel_halt(substr, status)
#endif

    ! READ DATA FROM ATTRIBUTES:
    ! - START DATE AND TIME
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'start_date_time' &
            , c=str_date_time)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(str_date_time, p_io)
    READ(str_date_time,*) yyyymmdd, hhmmss
#if defined(ECHAM5)
#ifndef MESSYTIMER
    CALL TC_set (yyyymmdd, hhmmss, io_date)
    CALL TC_convert(io_date, start_date) 
#else
    CALL TC_set (yyyymmdd, hhmmss, io_date)
    CALL TC_convert(io_date, E5_start_date) 
#endif
#endif
#ifdef MESSYTIMER
    iymd = SIGN(1,INT(yyyymmdd))*INT(yyyymmdd)
    kyr = INT(iymd/10000)
    kmo = MOD(iymd,10000)/100
    kdy = MOD(iymd,100)
    khr = hhmmss/10000
    kmn = MOD(hhmmss,10000)/100
    kse = MOD(hhmmss,100)

!!$#ifdef COSMO! um_ak_20100630+
!!$    IF (.NOT. lforcedtime) THEN ! um_ak_20100630
!!$       IF (kyr /= YEAR_START .OR. kmo /= MONTH_START .OR. kdy /= DAY_START .OR. &
!!$            khr /= HOUR_START .OR. kmn /= MINUTE_START .OR. kse /= SECOND_START &
!!$         ) THEN
!!$          IF (p_parallel_io) write (*,*) 'COSMO START DATE (', YEAR_START &
!!$               ,MONTH_START, DAY_START,' ',HOUR_START,MINUTE_START,SECOND_START &
!!$               , ') DOES NOT MATCH RESTART START_DATE (',str_date_time,')' 
!!$          CALL channel_halt(substr,3210)
!!$      ENDIF
!!$   ENDIF
!!$#endif! um_ak_20100630-

    CALL timer_set_date(status,'start',kyr, kmo, kdy, khr, kmn, kse)
#endif

    ! - RESTART DATE
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'restart_date_time' &
            , c=str_date_time)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(str_date_time, p_io)
    READ(str_date_time,*) yyyymmdd, hhmmss
#if defined(ECHAM5)
#ifndef MESSYTIMER
    CALL TC_set (yyyymmdd, hhmmss, io_date)
    CALL TC_convert(io_date, resume_date)
#else
    CALL TC_set (yyyymmdd, hhmmss, io_date)
    CALL TC_convert(io_date, E5_resume_date)
#endif
#endif
#ifdef MESSYTIMER
    iymd = SIGN(1,INT(yyyymmdd))*INT(yyyymmdd)
    kyr = INT(iymd/10000)
    kmo = MOD(iymd,10000)/100
    kdy = MOD(iymd,100)
    khr = hhmmss/10000
    kmn = MOD(hhmmss,10000)/100
    kse = MOD(hhmmss,100)
!!$#ifdef COSMO
!!$    ! um_ak_20100630
!!$    IF (kyr /= YEAR .OR. kmo /= MONTH .OR. kdy /= DAY .OR. &
!!$         khr /= HOUR .OR. kmn /= MINUTE .OR. kse /= SECOND &
!!$         ) THEN
!!$       IF (p_parallel_io) write (*,*) 'COSMO START DATE (', YEAR &
!!$            ,MONTH, DAY,' ',HOUR,MINUTE,SECOND &
!!$            , ') DOES NOT MATCH RESTART START_DATE (',str_date_time,')' 
!!$       CALL channel_halt(substr,3210)
!!$    ENDIF
!!$    ! um_ak_20100630-
!!$#endif
    CALL timer_set_date(status,'resume',kyr, kmo, kdy, khr, kmn, kse)
#endif

    ! - CURRENT TIME STEP
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nstep' &
            , i=nstep)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(nstep, p_io)

    ! - TIME STEP LENGTH
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'timestep' &
            , i=timestep)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(timestep, p_io)   

#ifdef COSMO
    IF (timestep /= INT(dt)) THEN
       IF (p_parallel_io) THEN
          write (*,*) 'COSMO TIME STEP (',dt &
               ,') DOES NOT MATCH RESTART TIME STEP (',timestep,')'
          write (*,*) 'CHANGE OF TIMESTEP NOT POSSIBLE IN COSMO'
          CALL channel_halt(substr,3210)
       ENDIF
    ENDIF
    ! radiation trigger
    IF (.NOT. lp) CALL p_bcast(hnextrad, p_io)   
    ! - coordinate type
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'hnextrad' &
            , r=hnextrad)
       CALL channel_halt(substr, status)
    END IF
    !
    ! - coordinate type
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'ivctype' &
            , i=ivctype)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(ivctype, p_io)   
    ! - coordinate type
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'irefatm' &
            , i=irefatm)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(irefatm, p_io)   
    ! - coordinate type 
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nfltvc' &
            , i=nfltvc)
       CALL channel_halt(substr, status)
    END IF
    !
    IF (.NOT. lp) CALL p_bcast(nfltvc, p_io)   

    ! the indices for the timelevels nnow, nnew, nold are rotating
    ! thus the information about the values of these three varaibles
    ! must be saved in the restart files in check in "init_restart"

    ! GET INDEX OF TIME LEVEL NNOW
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nnow' &
            , i=nnow)
       CALL channel_halt(substr, status)
    END IF
    IF (.NOT. lp) CALL p_bcast(nnow, p_io)   

    ! GET INDEX OF TIME LEVEL NNEW
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nnew' &
            , i=nnew)
       CALL channel_halt(substr, status)
    END IF
    IF (.NOT. lp) CALL p_bcast(nnew, p_io)   

    ! GET INDEX OF TIME LEVEL NOLD
    IF (p_parallel_io .OR. lp) THEN
       CALL get_attribute(status, restart_att, 'nold' &
            , i=nold)
       CALL channel_halt(substr, status)
    END IF
    IF (.NOT. lp) CALL p_bcast(nold, p_io)   
#endif

#if defined(ECHAM5)
#ifdef MESSYTIMER
    CALL timer_sync
#endif

    ! (RE-)INITIALIZE ECHAM5 TIME MANAGER
    CALL ec_manager_init(timestep, nstep)
#endif
#ifdef COSMO
    CALL messy_timer_COSMO_reinit_time
#endif

#ifdef MESSYTIMER
    CALL messy_timer_init_manager(timestep, nstep)
#endif

  END SUBROUTINE initialize_restart_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE bi_decompose(status, flag, reprid, gptr, ptr)
    
#if defined(ECHAM5)
    ! ECHAM5
    USE messy_main_mpi_bi,       ONLY: gather_sa,  gather_sp   &  ! ECHAM5
                                     , scatter_sa, scatter_sp  &  ! ECHAM5
                                     , dcg                        ! ECHAM5
    USE messy_main_transform_bi, ONLY: gather_glix, scatter_glix
    ! op_pj_20110205+
    USE messy_main_mpi_bi,       ONLY: gather_field
    ! op_pj_20110205-
#endif
! mz_ab_20100308+
#if defined(ECHAM5) || defined(MBM_CMAT)
    USE messy_main_transform_bi, ONLY: gather_cg, scatter_cg 
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_decomp_bi,    ONLY: dc
#endif
#if defined(ECHAM5) || (defined COSMO) || defined(BLANK)
! mz_ab_20100308-
    ! BM/MESSy
    USE messy_main_mpi_bi,       ONLY: gather_gp, scatter_gp   &  ! E5, C4
                                     , p_parallel_io, p_io, p_bcast
! mz_ab_20100308+
#endif
! mz_ab_20100308-
#if defined(I2CINC)
    USE messy_main_mpi_bi,       ONLY: switch_par_utilities
    USE messy_main_blather_bi,   ONLY: info_bi
#endif
    ! op_bk_20130820+
#if defined(__ICON__)
    USE mo_kind,                 ONLY: wp, sp
    USE mo_mpi,                  ONLY: my_process_is_mpi_workroot, my_process_is_mpi_seq
    ! op_bk_20131210+
!     USE mo_gather_scatter,       ONLY: gather_cells, scatter_cells
    ! op_bk_20131210-
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE mo_model_domain,         ONLY: p_patch, p_phys_patch, t_phys_patch
    USE mo_grid_config,          ONLY: n_phys_dom
    USE mo_parallel_config,      ONLY: nproma
    USE mo_communication,        ONLY: exchange_data, t_comm_gather_pattern, idx_no, blk_no
    USE mo_name_list_output_types,    ONLY:  l_output_phys_patch, t_reorder_info
    USE mo_name_list_output_init,     ONLY:  patch_info
#endif
    ! op_bk_20130820-
    ! MESSy
    USE messy_main_channel_repr, ONLY: t_representation, get_representation

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    INTEGER,                      INTENT(OUT) :: status
    INTEGER,                      INTENT(IN)  :: flag
    INTEGER,                      INTENT(IN)  :: reprid
    REAL(DP), DIMENSION(:,:,:,:), POINTER     :: gptr
    REAL(DP), DIMENSION(:,:,:,:), POINTER     :: ptr

    ! LOCAL
    CHARACTER(LEN=*),       PARAMETER :: substr = 'bi_decompose'
    TYPE(t_representation), POINTER   :: repr
    INTEGER :: stat
    INTEGER :: npoints, nblks, n_own, nlev, n_points, nlevs
    INTEGER :: i, ib, il, jk, jp, n, i_dom
    LOGICAL, ALLOCATABLE :: phys_owner_mask(:) ! owner mask for physical patch
    INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:)
    INTEGER, ALLOCATABLE :: own_idx(:), own_blk(:), own_dst_idx(:), own_dst_blk(:)
    ! op_bk_20131212+
!     REAL(wp), ALLOCATABLE :: r_tmp(:,:,:)
    REAL(wp), ALLOCATABLE :: r_tmp(:,:)
    ! op_bk_20131212-
    REAL(wp), ALLOCATABLE :: r_ptr(:,:,:)
    TYPE(t_phys_patch), POINTER   :: p_phys
    TYPE(t_comm_gather_pattern), POINTER :: p_pat
    TYPE(t_reorder_info),  POINTER :: p_ri

    LOGICAL :: lflipranks

    ! INIT
    CALL get_representation(status, reprid, repr)
    CALL channel_halt(substr, status)

    SELECT CASE(flag) 
    CASE(-1) 
       ! ################ RE-COMPOSE ##################################
       !
       ! INIT
       IF (ASSOCIATED(gptr)) THEN
          DEALLOCATE(gptr)
          NULLIFY(gptr)
       END IF
       !
       IF (p_parallel_io) THEN
          WRITE(*,'(A14, 4I6, A3)') "MESSY: gptr(/ ", repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4), " /)"
          ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          IF (stat /= 0) THEN
             status = 1000 ! memory allocation failed
             CALL channel_halt(substr, status)
          END IF
       END IF
       !
       SELECT CASE(repr%dctype)
       CASE(DC_BC)
          ! NOTE: MUST BE SYNCHRONIZED ON ALL PEs
          IF (p_parallel_io) &
               gptr(:,:,:,:) = ptr(:,:,:,:)
#ifdef COSMO
       CASE(DC_GP)
!qqq      ! NOTE: The following allocation / deallocation on non-IO PEs
          !       is a workaround to bypass Lahey/Fujitsu runtime checks,
          !       which obviously do not allow to pass a nullified pointer
          !       as parameter to a subroutine.
          IF (.NOT.p_parallel_io) THEN 
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          ! um_ak_20110603  lg32 added to enable GP_3D_1LEV  
          CALL gather_gp(gptr, ptr, lg32=(reprid == GP_3D_1LEV)) 
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
#ifdef I2CINC
       CASE(DC_I2C)
          CALL switch_par_utilities (1)
!qqq      ! NOTE: The following allocation / deallocation on non-IO PEs
          !       is a workaround to bypass Lahey/Fujitsu runtime checks,
          !       which obviously do not allow to pass a nullified pointer
          !       as parameter to a subroutine.
          IF (.NOT.p_parallel_io) THEN 
             ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
                  , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          CALL gather_gp(gptr, ptr)
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
          CALL switch_par_utilities (2)
       CASE(DC_I2C_IN)
          IF (ASSOCIATED(gptr)) DEALLOCATE(gptr) 
          NULLIFY(gptr)
          CALL info_bi('OUTPUT FOR I2C INPUT FIELDS NOT IMPLEMENTED' &
               , substr)
#endif
#endif
#if defined(ECHAM5)
       CASE(DC_GP)
! op_pj_20110205+
!!$          CALL gather_gp(gptr, ptr, dcg)
          CALL gather_field(gptr, repr%gdimlen, ptr)
! op_pj_20110205-
       CASE(DC_SP)
          CALL gather_sp(gptr, ptr, dcg)
       CASE(DC_SA)
          CALL gather_sa(gptr, ptr, dcg)
       CASE(DC_IX)
          CALL gather_glix(gptr, ptr, 1)
#endif
! mz_ab_20100210+
#if defined(ECHAM5) || defined(MBM_CMAT)
       CASE(DC_CG)
          CALL gather_cg(gptr, ptr, dc%cgg)          
#endif
! mz_ab_20100210-
          ! op_bk_20130906+
#ifdef  __ICON__
       CASE(DC_GP)
          WRITE(*,*) "MESSY: bi_decompose ICON DC_GP"
          i_dom = repr%patch_id
          WRITE(*,'(A19,I2)') "MESSY: my patch id", i_dom
          IF (i_dom > 0) THEN
             WRITE(*,'(A27,I2)') "MESSY: number phys domains", n_phys_dom
             DO jp=1, n_phys_dom
                IF (p_phys_patch(jp)%logical_id == i_dom) THEN
                   p_phys => p_phys_patch(jp)
                END IF
             END DO
             WRITE(*,*) "MESSY: my repr name " //TRIM(repr%name)
             IF (repr%dim(1)%ptr%id == DIMID_NCELLS_DOM(i_dom)) THEN
                WRITE(*,*) "MESSY: bi_decompose ICON CELLS"
                p_ri => patch_info(i_dom)%cells
                p_pat => patch_info(i_dom)%p_pat_c
!!$                IF(l_output_phys_patch) THEN
!!$                   p_pat => p_phys_patch(i_dom)%comm_pat_gather_c
!!$!                   n_points = p_phys_patch(i_dom)%n_patch_cells
!!$                ELSE
!!$                   p_pat => p_patch(i_dom)%comm_pat_gather_c
!!$!                   n_points = p_patch(i_dom)%n_patch_cells_g
!!$                ENDIF
             ELSE IF (repr%dim(1)%ptr%id == DIMID_NEDGES_DOM(i_dom)) THEN
                WRITE(*,*) "MESSY: bi_decompose ICON EDGES"
                p_ri => patch_info(i_dom)%edges
                p_pat => patch_info(i_dom)%p_pat_e
!!$                IF(l_output_phys_patch) THEN
!!$                   p_pat => p_phys_patch(i_dom)%comm_pat_gather_e
!!$!                   n_points = p_phys_patch(i_dom)%n_patch_edges
!!$                ELSE
!!$                   p_pat => p_patch(i_dom)%comm_pat_gather_e
!!$!                   n_points = p_patch(i_dom)%n_patch_edges_g
!!$                ENDIF
             ELSE IF (repr%dim(1)%ptr%id == DIMID_NVERTS_DOM(i_dom)) THEN
                WRITE(*,*) "MESSY: bi_decompose ICON VERTS"
                p_ri => patch_info(i_dom)%verts
                p_pat => patch_info(i_dom)%p_pat_v
!!$                IF(l_output_phys_patch) THEN
!!$                   p_pat => p_phys_patch(i_dom)%comm_pat_gather_v
!!$!                   n_points = p_phys_patch(i_dom)%n_patch_verts
!!$                ELSE
!!$                   p_pat => p_patch(i_dom)%comm_pat_gather_v
!!$!                   n_points = p_patch(i_dom)%n_patch_verts_g
!!$                ENDIF
             ELSE
                status = 1
                RETURN
!                CALL finish(routine,'unknown grid type')
             END IF

             ! op_bk_20131120+
             lflipranks = .FALSE.
             ! op_bk_20131120-

             n_points = p_ri%n_glb
             nblks = (n_points-1)/nproma + 1
             IF (repr%rank == 2) THEN
                nlevs = 1
             ELSE
                ! op_bk_20131129+
                IF (reprid == UNSTRUCTURED_CELL_2LEV_DOM(i_dom) .OR.          &
                   & reprid == UNSTRUCTURED_EDGE_2LEV_DOM(i_dom)) THEN
                   nlevs = 2
                ELSE IF (reprid == UNSTRUCTURED_CELL_3LEV_DOM(i_dom) .OR.          &
                   & reprid == UNSTRUCTURED_EDGE_3LEV_DOM(i_dom)) THEN
                   nlevs = 3
                ELSE
                   ! op_bk_20131120+
                   IF (ABS(p_patch(i_dom)%nlev - repr%ldimlen(2)) <= 1) THEN
                      nlevs = repr%ldimlen(2)
!                      nlevs = p_patch(i_dom)%nlev
                   ELSE
                      nlevs = repr%ldimlen(3)
                      lflipranks = .TRUE.
                   END IF
                   ! op_bk_20131120-
                END IF
                ! op_bk_20131129-
             END IF

             ! op_bk_20131212+
!              IF(my_process_is_mpi_workroot()) THEN
!                 ALLOCATE(r_tmp(nproma,nlevs,nblks))
!              ELSE
!                 ! Dimensions 1 and 2 of r_tmp must always be nproma and nlevs,
!                 ! otherwise exchange_data doesn't work!
!                 ALLOCATE(r_tmp(nproma,nlevs,1))
!              ENDIF
             ! op_bk_20131212-
             ! op_bk_20131118+
             IF (nlevs == 1) THEN
             ! op_bk_20131118-
                ALLOCATE(r_ptr(repr%ldimlen(1),1,repr%ldimlen(2)))
                WRITE(*,'(a26,i4,a5,i4,a3)') "MESSY: Allocate r_ptr: (/ ", repr%ldimlen(1), ", 1, ", repr%ldimlen(2), " /)"
                WRITE(*,'(a20,l2)') "MESSY: lflipranks: ", lflipranks
                r_ptr(:,1,:) = ptr(:,:,1,1)
             ELSE
                ! op_bk_20131120+
                IF (.NOT. lflipranks) THEN
                   ALLOCATE(r_ptr(repr%ldimlen(1),repr%ldimlen(2),repr%ldimlen(3)))
                   WRITE(*,'(a26,i4,a2,i4,a2,i4,a3)') "MESSY: Allocate r_ptr: (/ ", repr%ldimlen(1), ", ", repr%ldimlen(2) &
                        & , ", ", repr%ldimlen(3)," /)"
                   WRITE(*,'(a20,l2)') "MESSY: lflipranks: ", lflipranks
                   r_ptr = ptr(:,:,:,1)
                ELSE
                   ALLOCATE(r_ptr(repr%ldimlen(1),repr%ldimlen(3),repr%ldimlen(2)))
                   WRITE(*,'(a26,i4,a2,i4,a2,i4,a3)') "MESSY: Allocate r_ptr: (/ ", repr%ldimlen(1), ", ", repr%ldimlen(3) &
                        & , ", ", repr%ldimlen(2)," /)"
                   DO jk = 1, nlevs
                      r_ptr(:,jk,:) = ptr(:,:,jk,1)
                   END DO
                END IF
                ! op_bk_20131120-
             ENDIF
             ! op_bk_20131212+
              ALLOCATE(r_tmp(MERGE(n_points, 0, &
                   &                      my_process_is_mpi_workroot()), nlevs))
             ! r_tmp(:,:,:) = 0._wp
             r_tmp(:,:) = 0._wp
             ! op_bk_20131212-
             ! Gather data on root
             IF (p_parallel_io) THEN
                gptr(:,:,:,:) = 0._wp
             END IF
             ! op_bk_20130909+
             ! op_bk_20131212+
!              IF(my_process_is_mpi_seq()) THEN
!                 DO jk = 1, nlevs
!                    DO i = 1, p_ri%n_own
!                       r_tmp(p_ri%own_dst_idx(i),jk,p_ri%own_dst_blk(i)) = r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i))
!                    ENDDO
!                 ENDDO
!              ELSE
             ! op_bk_20131212-
                ! op_bk_20130909-
                CALL exchange_data(r_ptr(:,:,:), r_tmp(:,:), p_pat)
!                CALL exchange_data(p_pat, RECV=r_tmp, SEND=r_ptr)
                ! op_bk_20130909+
                ! op_bk_20131212+
!              ENDIF
                ! op_bk_20131212-
                ! op_bk_20130909-
                ! op_bk_20130909+
             IF(my_process_is_mpi_workroot()) THEN
                ! op_bk_20130909-
!                IF (p_parallel_io) THEN
                ! op_bk_20131212+
!                 DO jk = 1, nlevs
!                    !                r_out_sp(:,jk) = REAL(RESHAPE(r_tmp(:,jk,:), (/ n_points /)), dp)
!                    gptr(:,jk,1,1) = REAL(RESHAPE(r_tmp(:,jk,:), (/ n_points /)), dp)
!                 ENDDO
                gptr(:,:,1,1) = r_tmp(:,:)
                ! op_bk_20131212-
             END IF
          ELSE
             status = 1
          END IF


          ! op_bk_20130909+
! !!!!!!!!!!!!!!!
!           i_dom = repr%patch_id
!           IF (i_dom > 0) THEN
!              DO jp=1, n_phys_dom
!                 IF (p_phys_patch(jp)%logical_id == i_dom) THEN
!                    p_phys => p_phys_patch(i)
!                 END IF
!              END DO
!              p_pat => p_patch(i_dom)%comm_pat_gather_c
!              npoints = p_phys%n_patch_cells
!              nblks = (npoints-1)/nproma + 1
!              nlev = p_patch(i_dom)%nlev

!              ALLOCATE(phys_owner_mask(npoints))

!              DO i = 1, npoints
!                 il = idx_no(i)
!                 ib = blk_no(i)
!                 phys_owner_mask(i) = p_patch(i_dom)%cells%owner_mask(il,ib)
! !                IF(l_output_phys_patch) &
! !                     phys_owner_mask(i) = phys_owner_mask(i) .AND. (phys_id(il,ib) == phys_patch_id)
!              END DO
!              n_own = COUNT(phys_owner_mask(:))

!              ALLOCATE(own_idx(n_own))
!              ALLOCATE(own_blk(n_own))
!              ALLOCATE(glbidx_own(n_own)) ! Global index of my own points

!              n = 0
!              DO i = 1, npoints
!                 IF(phys_owner_mask(i)) THEN
!                    n = n+1
!                    own_idx(n) = idx_no(i)
!                    own_blk(n) = blk_no(i)
!                    glbidx_own(n) = p_patch(i_dom)%cells%glb_index(i)
!                 ENDIF
!              ENDDO

!              ! set trivial destination indices:
!              IF(my_process_is_mpi_seq()) THEN
!                 ALLOCATE(own_dst_idx(n_own), &
!                      &   own_dst_blk(n_own))
!                 DO i=1,n_own
!                    own_dst_idx(i) = idx_no(i)
!                    own_dst_blk(i) = blk_no(i)
!                 END DO ! i
!              END IF

!              IF (my_process_is_mpi_workroot()) THEN
!                 ALLOCATE(r_tmp(nproma, nlev, nblks))
!              ELSE
!                 ALLOCATE(r_tmp(nproma, nlev, 1))
!              END IF
!              r_ptr => ptr(:,:,:,1)
!              r_tmp(:,:,:) = 0._wp
!              gptr(:,:,:,:) = 0._wp
!              ! Gather data on root
!              IF(my_process_is_mpi_seq()) THEN
!                 DO jk = 1, nlev
!                    DO i = 1, n_own
!                       !                gptr(own_dst_idx(i),jk,own_dst_blk(i),1) = ptr(own_idx(i),jk,own_blk(i),1)
!                       r_tmp(own_dst_idx(i),jk,own_dst_blk(i)) = r_ptr(own_idx(i),jk,own_blk(i))
!                    ENDDO
!                 ENDDO
!              ELSE
!                 CALL exchange_data(p_pat, RECV=r_tmp, SEND=r_ptr)
! !                CALL exchange_data(p_pat, RECV=gptr, SEND=ptr)
!              ENDIF
!           ELSE
!              status = 1
!           END IF
! ! De-blocking
!           IF(my_process_is_mpi_workroot()) THEN
!              DO jk = 1, nlev
          ! op_bk_20130909-
!                 gptr(:,jk,1,1) = REAL(RESHAPE(r_tmp(:,jk,:), (/ npoints /)), dp)
!              ENDDO
!           ENDIF
          ! op_bk_20130906-
#endif
       CASE DEFAULT
          status = 1
       END SELECT
       !
       ! ##############################################################
    CASE(1)
       ! ################ DE-COMPOSE ##################################
       !
       ! INIT
       IF (ASSOCIATED(ptr)) THEN
          DEALLOCATE(ptr)
          NULLIFY(ptr)
       END IF
       !
       ! mz_ab_20100530+
       ! subtract bounds
       !!$ALLOCATE( ptr( repr%ldimlen(1), repr%ldimlen(2) &
       !!$     , repr%ldimlen(3), repr%ldimlen(4) ), STAT=stat )
       ALLOCATE( ptr( repr%ldimlen(1)-2*repr%bounds%nbounds(1) &
                    , repr%ldimlen(2)-2*repr%bounds%nbounds(2) &
                    , repr%ldimlen(3)-2*repr%bounds%nbounds(3) & 
                    , repr%ldimlen(4)-2*repr%bounds%nbounds(4) ), STAT=stat )
       ! mz_ab_20100530-
       IF (stat /= 0) THEN
          status = 1000 ! memory allocation failed
          CALL channel_halt(substr, status)
       END IF
       !
       SELECT CASE(repr%dctype)
       CASE(DC_BC)
          ! BROADCAST FROM IO-PE TO ALL OTHERS
          IF (p_parallel_io) ptr(:,:,:,:) = gptr(:,:,:,:)
          CALL p_bcast(ptr, p_io)
#ifdef COSMO
       CASE(DC_GP)
!qqq      ! NOTE: The following allocation / deallocation on non-IO PEs
          !       is a workaround to bypass Lahey/Fujitsu runtime checks,
          !       which obviously do not allow to pass a nullified pointer
          !       as parameter to a subroutine.
          IF (.NOT.p_parallel_io) THEN 
              ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          ! um_ak_20110603  lg32 added to enable GP_3D_1LEV  
          CALL scatter_gp(gptr, ptr, lg32=(reprid == GP_3D_1LEV)) 
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
#ifdef I2CINC
       CASE (DC_I2C)
          CALL switch_par_utilities (1)
!qqq      ! NOTE: The following allocation / deallocation on non-IO PEs
          !       is a workaround to bypass Lahey/Fujitsu runtime checks,
          !       which obviously do not allow to pass a nullified pointer
          !       as parameter to a subroutine.
          IF (.NOT.p_parallel_io) THEN 
              ALLOCATE( gptr( repr%gdimlen(1), repr%gdimlen(2) &
               , repr%gdimlen(3), repr%gdimlen(4) ), STAT=stat )
          ENDIF
          CALL scatter_gp(gptr, ptr)
          IF (.NOT.p_parallel_io) THEN
            DEALLOCATE(gptr)
            NULLIFY(gptr)
          ENDIF
          CALL switch_par_utilities (2)
          
       CASE (DC_I2C_IN)   
          
          CALL info_bi('SORRY, NO INPUT OF INT2COSMO INPUT FIELDS POSSIBLE' &
               , substr)
          IF (ASSOCIATED(gptr)) DEALLOCATE(gptr)
          NULLIFY(gptr)
#endif
#endif
#if defined(ECHAM5)
       CASE(DC_GP)
          CALL scatter_gp(gptr, ptr, dcg)
       CASE(DC_SP)
          CALL scatter_sp(gptr, ptr, dcg)
       CASE(DC_SA)
          CALL scatter_sa(gptr, ptr, dcg)
       CASE(DC_IX)
          CALL scatter_glix(gptr, ptr, 1)
#endif
! mz_ab_20100210+
#if defined(ECHAM5) || defined(MBM_CMAT)
       CASE(DC_CG)
          CALL scatter_cg(gptr, ptr, dc%cgg)          
#endif
! mz_ab_20100210-
       CASE DEFAULT
          status = 1
       END SELECT
       !
       ! ##############################################################
    END SELECT

    status = 0

  END SUBROUTINE bi_decompose
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE bi_vector(status, flag, reprid, ptr, vptr)

    ! BM/MESSy
    USE messy_main_mpi_bi,       ONLY: dcl, reorder
    ! MESSy
    USE messy_main_channel_repr, ONLY: t_representation, get_representation

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE

    ! I/O
    INTEGER,                      INTENT(OUT) :: status
    INTEGER,                      INTENT(IN)  :: flag
    INTEGER,                      INTENT(IN)  :: reprid
    REAL(DP), DIMENSION(:,:,:,:), POINTER     :: ptr
    REAL(DP), DIMENSION(:,:,:,:), POINTER     :: vptr

    ! LOCAL
    CHARACTER(LEN=*),       PARAMETER :: substr = 'bi_vector'
    TYPE(t_representation), POINTER   :: repr
    INTEGER :: stat
    INTEGER :: i1, i2, n1, n2

    ! INIT
    CALL get_representation(status, reprid, repr)
    CALL channel_halt(substr, status)

    SELECT CASE(flag) 
    CASE(-1) 
       ! ################ DE-VECTOR ##################################
       !
       ! INIT
       IF (ASSOCIATED(ptr)) THEN
          DEALLOCATE(ptr)
          NULLIFY(ptr)
       END IF
       !
       IF (.NOT. repr%pdecomp%lpdecomp) THEN
          ! REPRESENTATION DECOMPOSITION TABLE DOES NOT EXIST
          status = 2025
          CALL channel_halt(substr, status)
       END IF
       !
       ALLOCATE( ptr( repr%pdecomp%shape_mem(1), repr%pdecomp%shape_mem(2) &
            , repr%pdecomp%shape_mem(3), repr%pdecomp%shape_mem(4) ) &
            , STAT=stat )
       IF (stat /= 0) THEN
          status = 1000 ! memory allocation failed
          CALL channel_halt(substr, status)
       END IF
       !
       SELECT CASE(repr%dctype)
       CASE(DC_BC)
!!$          ptr(:,:,:,:) = vptr(:,:,:,:)
          ptr = vptr
       CASE(DC_GP)
          IF (dcl%lreg) THEN
!!$             ptr(:,:,:,:) = vptr(:,:,:,:)
             ptr = vptr
          ELSE
             SELECT CASE (repr%rank)
             CASE(1)
                ! UNKNOWN REPRESENTATION DECOMPOSITION TYPE / RANK
                status = 2022
                CALL channel_halt(substr, status)
             CASE(2)
                IF (repr%id == GP_3D_1LEV) THEN
                   CALL reorder(ptr(:,1,:,1), vptr(:,1,:,1))
                ELSE
                   CALL reorder(ptr(:,:,1,1), vptr(:,:,1,1))
                END IF
             CASE(3)
                n1 = SIZE(vptr, 2)
                DO i1=1, n1
                   CALL reorder(ptr(:,i1,:,1), vptr(:,i1,:,1))
                END DO
             CASE(4)
                n1 = SIZE(vptr, 2)
                n2 = SIZE(vptr, 3)
                DO i1=1, n1
                   DO i2=1, n2
                      CALL reorder(ptr(:,i1,i2,:), vptr(:,i1,i2,:))
                   END DO
                END DO
             END SELECT
          END IF
       CASE(DC_SP)
!!$          ptr(:,:,:,:) = vptr(:,:,:,:)
          ptr = vptr
       CASE(DC_SA)
!!$          ptr(:,:,:,:) = vptr(:,:,:,:)
          ptr = vptr
       CASE(DC_IX)
!!$          ptr(:,:,:,:) = vptr(:,:,:,:)
          ptr = vptr
! mz_ab_20100504+
       CASE(DC_CG)
!!$          ptr(:,:,:,:) = vptr(:,:,:,:)
          ptr = vptr
! mz_ab_20100504-
       CASE DEFAULT
          status = 1
       END SELECT
       !
       ! ##############################################################
    CASE(1)
       ! ################ VECTOR ######################################
       !
       ! INIT
       IF (ASSOCIATED(vptr)) THEN
          DEALLOCATE(vptr)
          NULLIFY(vptr)
       END IF
       !
       ALLOCATE( vptr( repr%ldimlen(1), repr%ldimlen(2) &
            , repr%ldimlen(3), repr%ldimlen(4) ), STAT=stat )
       IF (stat /= 0) THEN
          status = 1000 ! memory allocation failed
          CALL channel_halt(substr, status)
       END IF
       !
       SELECT CASE(repr%dctype)
       CASE(DC_BC)
!!$          vptr(:,:,:,:) = ptr(:,:,:,:)
          vptr = ptr
       CASE(DC_GP)
          IF (dcl%lreg) THEN
!!$             vptr(:,:,:,:) = ptr(:,:,:,:)
             vptr = ptr
          ELSE
             SELECT CASE (repr%rank)
             CASE(1)
                ! UNKNOWN REPRESENTATION DECOMPOSITION TYPE / RANK
                status = 2022
                CALL channel_halt(substr, status)
                ! NOT EXISTENT
             CASE(2)
                IF (repr%id == GP_3D_1LEV) THEN
                   CALL reorder(vptr(:,1,:,1), ptr(:,1,:,1))
                ELSE
                   CALL reorder(vptr(:,:,1,1), ptr(:,:,1,1))
                END IF
             CASE(3)
                n1 = SIZE(ptr, 2)
                DO i1=1, n1
                   CALL reorder(vptr(:,i1,:,1), ptr(:,i1,:,1))
                END DO
             CASE(4)
                n1 = SIZE(ptr, 2)
                n2 = SIZE(ptr, 3)
                DO i1=1, n1
                   DO i2=1, n2
                      CALL reorder(vptr(:,i1,i2,:), ptr(:,i1,i2,:))
                   END DO
                END DO
             END SELECT
          END IF
       CASE(DC_SP)
!!$          vptr(:,:,:,:) = ptr(:,:,:,:)
          vptr = ptr
       CASE(DC_SA)
!!$          vptr(:,:,:,:) = ptr(:,:,:,:)
          vptr = ptr
       CASE(DC_IX)
!!$          vptr(:,:,:,:) = ptr(:,:,:,:)
          vptr = ptr
! mz_ab_20100504+
       CASE(DC_CG)
!!$          vptr(:,:,:,:) = ptr(:,:,:,:)
          vptr = ptr
! mz_ab_20100504-
       CASE DEFAULT
          status = 1
       END SELECT
       !
       ! ##############################################################
    END SELECT

    status = 0

  END SUBROUTINE bi_vector
  ! -------------------------------------------------------------------


  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_init_timer

    ! BM/MESSy
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,   ONLY: error_bi
#ifdef COSMO
    USE messy_main_data_bi,      ONLY: ngribout
#endif
    USE messy_main_data_bi,      ONLY: l2tls ! um_ak_20091117
#ifndef MESSYTIMER
    USE messy_main_bmluse_bi,    ONLY: echam_ev_init, p_bcast_event
#else
    USE messy_main_timer_bi,     ONLY: timer_event_init, p_bcast_event
#endif
    ! MESSy
    USE messy_main_channel,      ONLY: get_channel_name
    USE messy_main_tools,        ONLY: match_wild, find_next_free_unit

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_init_timer'
    INTEGER :: status
    INTEGER :: iou 
    INTEGER :: i, js
    CHARACTER(LEN=STRLEN_CHANNEL) :: cname
    LOGICAL :: lexplicit
    CHARACTER(LEN=8) :: evaldate     ! um_ak_20091117
#ifdef COSMO
    INTEGER :: iCm, iCp, iCz, iCs, iCc
    
    ! dimension indices for COSMO output channels
    ALLOCATE(js_COSMOm(ngribout)); js_COSMOm(:) = -1 
    ALLOCATE(js_COSMOp(ngribout)); js_COSMOp(:) = -1 
    ALLOCATE(js_COSMOz(ngribout)); js_COSMOz(:) = -1 
    ALLOCATE(js_COSMOs(ngribout)); js_COSMOs(:) = -1 
    ALLOCATE(js_COSMOc(ngribout)); js_COSMOc(:) = -1 
    iCm = 0
    iCz = 0
    iCp = 0
    iCs = 0
    iCc = 0
#endif

       ! op_bk_20130904+
       WRITE(*, '(A24)') "MESSY: Pos. 1 reached."
       ! op_bk_20130904-
!!$    ! INITIALIZE CPL
!!$    IF (p_parallel_io) THEN
!!$       iou = find_next_free_unit(100,200)
!!$       CALL main_channel_read_nml_cpl(status, iou)
!!$       IF (status /= 0) CALL error_bi(' ',substr)
!!$    END IF
!!$    ! BROADCAST RESULTS
!!$    CALL p_bcast(L_BM_ORIG_OUTPUT, p_io)
!!$    !
    ! um_ak_20081125+
    ! moved back here because of TIMER MANAGER initialization
    ! in case of RESTART
    CALL p_bcast(TIMER_DEFAULT%cname, p_io)
    CALL p_bcast_event(TIMER_DEFAULT%io_event, p_io)
    DO i=1, NMAXCHANNELS
       CALL p_bcast(TIMER_CHANNEL(i)%cname, p_io)
       CALL p_bcast_event(TIMER_CHANNEL(i)%io_event, p_io)
    END DO

    ! mz_pj_20080905+
    CALL p_bcast(TIMER_TNF_DEFAULT%cname, p_io)
    CALL p_bcast_event(TIMER_TNF_DEFAULT%io_event, p_io)
    DO i=1, NMAXCHANNELS
       CALL p_bcast(TIMER_TNF_CHANNEL(i)%cname, p_io)
       CALL p_bcast_event(TIMER_TNF_CHANNEL(i)%io_event, p_io)
    END DO  
    ! mz_pj_20080905-
    ! um_ak_20081125-

    ! SPACE FOR ACTIVE OUTPUT EVENTS (ONE PER CHANNEL)
    ALLOCATE(OUTPUT_EVENT(NCHANNEL))
    ALLOCATE(LOUTPUT_NOW(NCHANNEL))
    LOUTPUT_NOW(:) = .FALSE.

    ! mz_pj_20080905+
    ALLOCATE(TNF_EVENT(NCHANNEL))
    ALLOCATE(LTNF_NOW(NCHANNEL))
    LTNF_NOW(:) = .FALSE.
    ! mz_pj_20080905-
    
    ! PRESET NON-EXPLICIT TO DEFAULT
    channel_loop: DO js=1, NCHANNEL
       ! op_bk_20130904+
       WRITE(*, '(A24, I4)') "MESSY: Pos. 2 reached:", js
       ! op_bk_20130904-

       CALL get_channel_name(status, js, cname)
       CALL channel_halt(substr, status)

       ! SET EXPLICIT OR DEFAULT
       lexplicit = .FALSE.
       event_loop: DO i=1, NMAXCHANNELS
          IF (TRIM(ADJUSTL(TIMER_CHANNEL(i)%cname)) == '') CYCLE
          IF (match_wild(TRIM(ADJUSTL(TIMER_CHANNEL(i)%cname)),&
               TRIM(cname))) THEN
             lexplicit  = .TRUE.
             EXIT
          END IF
       END DO event_loop

       ! op_bk_20130904+
       WRITE(*, '(A24, L2, I4)') "MESSY: Pos. 3 reached:", lexplicit, js
       ! op_bk_20130904-
       ! INITIALIZE EVENT
#ifndef MESSYTIMER
#if defined(ECHAM5)
       IF (lexplicit) THEN
          CALL echam_ev_init(OUTPUT_EVENT(js), TIMER_CHANNEL(i)%io_event, &
               TRIM(cname), 'next')
       ELSE
          CALL echam_ev_init(OUTPUT_EVENT(js), TIMER_DEFAULT%io_event,   &
               TRIM(cname), 'next')
       END IF
#endif
#else
       ! Note: to trigger the output to exactly the time choosen in the namelist
       !       the evaluation date depends on the time integration scheme
       ! um_ak_20091106+
!!$       IF (l2tls) THEN
!!$          evaldate =  'present'
!!$       ELSE
          evaldate =  'next'
!!$       ENDIF
       ! um_ak_20091106-
       IF (lexplicit) THEN
          CALL timer_event_init(OUTPUT_EVENT(js), TIMER_CHANNEL(i)%io_event, &
               TRIM(cname), TRIM(evaldate))  ! um_ak_20091106 next replaced
       ELSE
          CALL timer_event_init(OUTPUT_EVENT(js), TIMER_DEFAULT%io_event,   &
               TRIM(cname),  TRIM(evaldate)) ! um_ak_20091106 next replaced
       END IF
#endif

#if defined(ECHAM5)
       ! SPECIAL FOR ECHAM5
       ! NOTE: CHANNEL OUTPUT OF tdiag, tdiag_gp, nudg, nudg_gp
       !       WILL BE SYNCHRONIZED TO THE SHORTEST OUTPUT INTERVAL
       !       (IF PRESENT !)
       SELECT CASE(TRIM(cname))
       CASE('tdiag')
          js_tdiag = js
       CASE('tdiag_gp')
          js_tdiag_gp = js
       CASE('nudg')
          js_nudg = js
       CASE('nudg_gp')
          js_nudg_gp = js
       CASE DEFAULT
          ! DO NOTHING
       END SELECT
#endif

#ifdef COSMO
       ! SPECIAL FOR COSMO
       ! calculation of OUTPUT fields in COSMO will be triggered,
       ! if ANY of the output channels is .TRUE.
       IF (cname(1:5) == 'COSMO') THEN
          SELECT CASE(cname(6:6))
          CASE('m')
             iCm =iCm + 1
             js_COSMOm(iCm) = js
          CASE('p')
             iCp = iCp + 1 
             js_COSMOp(iCp) = js
          CASE('z')
             iCz = iCz +1
             js_COSMOz(iCz) = js
          CASE('s')
             iCs = iCs + 1
             js_COSMOs(iCs) = js
          CASE('c')
             iCc = iCc + 1
             js_COSMOs(iCc) = js
          CASE DEFAULT
          ! DO NOTHING
          END SELECT
       ENDIF
#endif

       ! op_bk_20130904+
       WRITE(*, '(A24, L2, I4)') "MESSY: Pos. 4 reached:", lexplicit, js
       ! op_bk_20130904-

! mz_pj_20080905+
       ! SET EXPLICIT OR DEFAULT
       lexplicit = .FALSE.
       event_loop_tnf: DO i=1, NMAXCHANNELS
          IF (TRIM(ADJUSTL(TIMER_TNF_CHANNEL(i)%cname)) == '') CYCLE
          IF (match_wild(TRIM(ADJUSTL(TIMER_TNF_CHANNEL(i)%cname)),&
               TRIM(cname))) THEN
             lexplicit  = .TRUE.
             EXIT
          END IF
       END DO event_loop_tnf

       ! op_bk_20130904+
       WRITE(*, '(A24, L2, I4)') "MESSY: Pos. 5 reached:", lexplicit, js
       ! op_bk_20130904-
       ! INITIALIZE EVENT
#ifndef MESSYTIMER
#if defined(ECHAM5)
       IF (lexplicit) THEN
          CALL echam_ev_init(TNF_EVENT(js), TIMER_TNF_CHANNEL(i)%io_event, &
               TRIM(cname), 'next')
       ELSE
          CALL echam_ev_init(TNF_EVENT(js), TIMER_TNF_DEFAULT%io_event,   &
               TRIM(cname), 'next')
       END IF
#endif
#else
       ! op_bk_20130904+
       WRITE(*, '(A24, L2, I4)') "MESSY: Pos. 6 reached:", lexplicit, js
       ! op_bk_20130904-

       IF (lexplicit) THEN
          CALL timer_event_init(TNF_EVENT(js), TIMER_TNF_CHANNEL(i)%io_event, &
               TRIM(cname), 'next')
       ELSE
          CALL timer_event_init(TNF_EVENT(js), TIMER_TNF_DEFAULT%io_event,   &
               TRIM(cname), 'next')
       END IF
#endif
! mz_pj_20080905-

    END DO channel_loop

  END SUBROUTINE main_channel_init_timer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_update_timer

#if defined(ECHAM5)
    ! ECHAM5
    USE mo_diag_tendency,       ONLY: dio_index
    USE mo_nudging_buffer,      ONLY: nio_index
    USE messy_main_bmluse_bi,   ONLY: l_putdata
#ifndef MESSYTIMER
    ! ECHAM5/MESSy
    USE messy_main_bmluse_bi,   ONLY: event_state, next_date, current_date &
                                    , time_days ! um_ak_20091106
#endif
#endif

#ifdef COSMO
    ! COSMO/MESSY
    USE messy_main_data_bi,     ONLY: ngribout
#endif

#ifdef MESSYTIMER
    USE messy_main_timer_bi,    ONLY: event_state
    ! MESSy
    USE messy_main_timer,       ONLY: next_date, current_date, time_days ! um_ak_20091106
#endif
    USE messy_main_data_bi,     ONLY: l2tls ! um_ak_20091106
    USE messy_main_channel,     ONLY: trigger_channel_output &
                                    , set_channel_output
    USE messy_main_tools,       ONLY: int2str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_update_timer'
    INTEGER :: status
    INTEGER :: js
    LOGICAL :: lsync
    CHARACTER(LEN=3) :: str
    INTEGER          :: iCc
    TYPE(time_days)  :: date ! um_ak_20091106

    ! um_ak_20091106+
    ! Note: to trigger the output to exactly the time choosen in the namelist
    !       the evaluation date depends on the time integration scheme
!!$    IF (l2tls) THEN
!!$       date = current_date
!!$    ELSE
       date = next_date
!!$    ENDIF
    ! um_ak_20091106-
    DO js = 1, NCHANNEL
!       LOUTPUT_NOW(js) = event_state(OUTPUT_EVENT(js), next_date)
       LOUTPUT_NOW(js) = event_state(OUTPUT_EVENT(js), date) ! um_ak_20091106
    END DO

    ! op_bk_20130904+
!    WRITE(*,'(A22, 3I3)') "MESSY: channel output", date%day, date%second, NCHANNEL
!    WRITE(*,'(A3, L2, A3)') "(/", LOUTPUT_NOW(:), "/)"
!    LOUTPUT_NOW(:) = .TRUE.
    ! op_bk_20130904-

    ! mz_pj_20080905+
    DO js = 1, NCHANNEL
!       LTNF_NOW(js) = event_state(TNF_EVENT(js), next_date)
       LTNF_NOW(js) = event_state(TNF_EVENT(js), date) ! um_ak_20091106
    END DO    
    ! mz_pj_20080905-

#ifdef COSMO
    ! SPECIAL FOR COSMO
    ! set L_FORCE_calcout to TRUE, if any COSMO output channel is true
    L_FORCE_calcout = .FALSE.
    DO js =1, ngribout
       IF (js_COSMOm(js) > 0 ) &
            L_FORCE_calcout = L_FORCE_calcout .OR. LOUTPUT_NOW(js_COSMOm(js))
       IF (js_COSMOp(js) > 0 ) &
            L_FORCE_calcout = L_FORCE_calcout .OR. LOUTPUT_NOW(js_COSMOp(js))
       IF (js_COSMOz(js) > 0 ) &
            L_FORCE_calcout = L_FORCE_calcout .OR. LOUTPUT_NOW(js_COSMOz(js))
       IF (js_COSMOs(js) > 0 ) &
            L_FORCE_calcout = L_FORCE_calcout .OR. LOUTPUT_NOW(js_COSMOs(js))
    ENDDO
#endif

#if defined(ECHAM5)
    ! SPECIAL FOR ECHAM5
    ! NOTE: CHANNEL OUTPUT OF tdiag, tdiag_gp, nudg, nudg_gp
    !       WILL BE SYNCHRONIZED TO THE SHORTEST OUTPUT INTERVAL (IF PRESENT !)
    lsync = .FALSE.
    IF (js_tdiag > 0)    lsync = lsync .OR. LOUTPUT_NOW(js_tdiag)
    IF (js_tdiag_gp > 0) lsync = lsync .OR. LOUTPUT_NOW(js_tdiag_gp)
    IF (js_nudg > 0)     lsync = lsync .OR. LOUTPUT_NOW(js_nudg)
    IF (js_nudg_gp > 0)  lsync = lsync .OR. LOUTPUT_NOW(js_nudg_gp)
    !
    IF (js_tdiag + js_tdiag_gp > 0) l_putdata(dio_index) = lsync
    IF (js_nudg  + js_nudg_gp  > 0) l_putdata(nio_index) = lsync
    !
    IF (js_tdiag > 0)    LOUTPUT_NOW(js_tdiag)    = lsync
    IF (js_tdiag_gp > 0) LOUTPUT_NOW(js_tdiag_gp) = lsync
    IF (js_nudg > 0)     LOUTPUT_NOW(js_nudg)     = lsync
    IF (js_nudg_gp > 0)  LOUTPUT_NOW(js_nudg_gp)  = lsync
#endif

! mz_pj_20080905+
!!$    CALL trigger_channel_output(status, LOUTPUT_NOW, LFORCE_NEW_OUTPUT)
    CALL trigger_channel_output(status, LOUTPUT_NOW, LTNF_NOW &
         , LFORCE_NEW_OUTPUT)
! mz_pj_20080905-
    CALL channel_halt(substr, status)
    LFORCE_NEW_OUTPUT = .FALSE.

#ifdef COSMO
    ! suppress output of channels with constants at the beginning
    DO iCc = 1, ngribout
       IF (js_COSMOc(iCc) > 0) THEN
          CALL int2str(str,iCc)
          CALL set_channel_output(status,'COSMOc'//str, .FALSE.)
          CALL channel_halt(substr,status)
       ENDIF
    ENDDO
#endif
  END SUBROUTINE main_channel_update_timer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_read_nml_cpl(status, iou)

    ! MODULE ROUTINE (SMIL)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2004

    USE messy_main_tools,   ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_channel, ONLY: modstr

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CPL/ L_BM_ORIG_OUTPUT, TIMER_DEFAULT, TIMER_CHANNEL &
         , TIMER_TNF_DEFAULT, TIMER_TNF_CHANNEL ! mz_pj_20080905

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_read_nml_cpl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! DO NOT ENABLE ADDITIONAL ECHAM5 STANDARD STREAM OUTPUT PER DEFAULT
    ! (only needed for GRIB-template generation)
    L_BM_ORIG_OUTPUT = .FALSE.

    ! -> OTHER DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_channel_read_nml_cpl
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_main_channel_bi
! **********************************************************************
