! **********************************************************************+
MODULE messy_main_tracer_mem_bi
! **********************************************************************+

  ! MESSy
  USE messy_main_tracer, ONLY: dp, t_trinfo, t_trinfo_tp, t_trinfo_list &
       , modstr, ON, OFF, I_FORCE_COL, I_ADVECT, I_VDIFF, I_INTEGRATE, I_MIX &
       , I_HDIFF, I_RELAX ! um_ak_20081030 for COSMO

  ! SPECIAL FOR LAGRANGIAN TRACERS / CHANNEL OBJECTS
  USE messy_main_tools,  ONLY: PTR_4D_ARRAY  ! um_ak_20081030

  IMPLICIT NONE
  !PUBLIC is already default
  PRIVATE :: dp
  ! op_bk_20130820+
#ifdef __ICON__
  PUBLIC
#endif
  ! op_bk_20130820-
  SAVE

! op_pj_20120120+ moved to this position from messy_main_tracer_bi.f90
  ! GLOBAL SWITCHES FOR TRACER SETS
  LOGICAL :: L_GP   = .FALSE. ! GLOBAL SWITCH   ! mz_ab_20100226
  LOGICAL :: L_LG   = .FALSE. ! GLOBAL SWITCH
  LOGICAL :: L_OM   = .FALSE. ! GLOBAL SWITCH   ! mz_ap_20071023
  LOGICAL :: L_CMAT = .FALSE. ! GLOBAL SWITCH   ! um_ak_20100128
  LOGICAL :: L_EC   = .FALSE. ! GLOBAL SWITCH   ! um_ak_20100128
! op_pj_20120120-

  ! TRACER SET: GRIDPOINT ##################################################
  !
  ! NAME OF TRACER SET
  CHARACTER(LEN=*), PARAMETER              :: GPTRSTR = 'gp'
  CHARACTER(LEN=*), PARAMETER              :: gp_channel = modstr//'_gp'
  !
  ! POINTER TO TRACER META INFORAMTION
  TYPE(t_trinfo_list),             POINTER :: trlist_gp => NULL()
  !
  ! ARRAY OF TRACER INFO STRUCTS
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_gp => NULL()
  !
  ! NUMBER OF TRACERS
  INTEGER                                  :: ntrac_gp = 0
  !
  ! POINTER TO TRACER FIELDS
  ! - FOR INTERNAL MEMORY MANAGEMENT
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxt   => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtte => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtm1 => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtf  => NULL()
  !
  ! - FOR OUTSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xt    => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtte  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtm1  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtf   => NULL()
  !
  ! - FOR INSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxt   => NULL()
  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtte => NULL()
  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtm1 => NULL()
  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtf  => NULL()
  !
  ! #######################################################################

  ! TRACER SET: LAGRANGE (ATTILA) #########################################
  INTEGER, PARAMETER :: NCELL  = 0
  INTEGER, PARAMETER :: NGCELL = 0
  !
  ! NAME OF TRACER SET
  CHARACTER(LEN=*), PARAMETER              :: LGTRSTR = 'lg'
  CHARACTER(LEN=*), PARAMETER              :: lg_channel = modstr//'_lg'
  !
  ! POINTER TO TRACER META INFORAMTION
  TYPE(t_trinfo_list),             POINTER :: trlist_lg => NULL()
  !
  ! ARRAY OF TRACER INFO STRUCTS
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_lg => NULL()
  !
  ! NUMBER OF TRACERS
  INTEGER                                  :: ntrac_lg = 0
  !
  ! POINTER TO TRACER FIELDS
  ! - FOR INTERNAL MEMORY MANAGEMENT
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxt_a   => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxtte_a => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxtm1_a => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxtf_a  => NULL()
  ! - FOR OUTSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xt_a    => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xtte_a  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xtm1_a  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xtf_a   => NULL()
  ! - FOR INSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxt_a   => NULL()
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxtte_a => NULL()
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxtm1_a => NULL()
  REAL(DP), DIMENSION(:,:),       POINTER   :: qxtf_a  => NULL()
  !
  ! FOR LG->GP CONVERSION
  CHARACTER(LEN=*), PARAMETER              :: LGGPTRSTR    = 'lggp'
  CHARACTER(LEN=*), PARAMETER              :: lggp_channel = modstr//'_lggp'
  TYPE(t_trinfo_list),             POINTER :: trlist_lggp => NULL()
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_lggp     => NULL()
  INTEGER                                  :: ntrac_lggp  = 0
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER  :: pxt_lggp    => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER  :: xt_lggp     => NULL()
  REAL(DP), DIMENSION(:,:,:),     POINTER  :: qxt_lggp    => NULL()
  !
  ! DECOMPOSITION
  ! -> NUMBER OF CELLS ON THIS PE
  INTEGER                                   :: NLCELL    = 0
  ! -> CELLS ON THIS PE ?
  LOGICAL                                   :: LG_ACTIVE = .false.
  !
  INTEGER :: number_mix = 0
  ! #######################################################################

! mz_ap_20071023+
  ! TRACER SET: OCEANIC (MPIOM) #########################################
  !
  ! NAME OF TRACER SET
  CHARACTER(LEN=*), PARAMETER              :: OMTRSTR = 'om'
  CHARACTER(LEN=*), PARAMETER              :: om_channel = modstr//'_om'
  !
  ! POINTER TO TRACER META INFORAMTION
  TYPE(t_trinfo_list),             POINTER :: trlist_om => NULL()
  !
  ! ARRAY OF TRACER INFO STRUCTS
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_om => NULL()
  !
  ! NUMBER OF TRACERS
  INTEGER                                  :: ntrac_om = 0
  !
  ! POINTER TO TRACER FIELDS
  ! - FOR INTERNAL MEMORY MANAGEMENT
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER   :: pxt_om   => NULL()
  ! - FOR USE IN SUBMODELS
  REAL(DP), DIMENSION(:,:,:,:),   POINTER   :: xt_om    => NULL()
  !
  ! #######################################################################
! mz_ap_20071023-

! um_ak_20080731+
  ! DEFINE BOUNDARY DATA VARIABLES FOR REGIONAL MODELS
  REAL(DP), DIMENSION(:,:,:,:,:),    POINTER :: xt_bd    => NULL()
  TYPE (PTR_4D_ARRAY), DIMENSION(:), POINTER :: xt_array => NULL()
! um_ak_20080731-

! mz_ab_20090622+
  ! TRACER SET: CMAT GRIDPOINT ############################################
  !
  ! NAME OF TRACER SET
  CHARACTER(LEN=*), PARAMETER              :: CMATTRSTR = 'cmat'
  CHARACTER(LEN=*), PARAMETER              :: cmat_channel = modstr//'_cm'
  !
  ! POINTER TO TRACER META INFORAMTION
  TYPE(t_trinfo_list),             POINTER :: trlist_cmat => NULL()
  !
  ! ARRAY OF TRACER INFO STRUCTS
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_cmat => NULL()
  !
  ! NUMBER OF TRACERS
  INTEGER                                  :: ntrac_cmat = 0
  !
  ! POINTER TO TRACER FIELDS
  ! - FOR INTERNAL MEMORY MANAGEMENT
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxt_cmat   => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtte_cmat => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtm1_cmat => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtm2_cmat  => NULL()
  !
  ! - FOR OUTSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xt_cmat    => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtte_cmat  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtm1_cmat  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtm2_cmat   => NULL()
  !
!!$  ! Probably not needed because CMAT runs only outside local loop
!!$  ! - FOR INSIDE 'LOCAL LOOP'
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxt_cmat   => NULL()
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtte_cmat => NULL()
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtm1_cmat => NULL()
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtf_cmat  => NULL()
  !
  ! ####################################################################### 

  ! TRACER SET: ECHAM+CMAT GRIDPOINT ######################################
  !
  ! NAME OF TRACER SET
  CHARACTER(LEN=*), PARAMETER              :: ECTRSTR = 'ec'
  CHARACTER(LEN=*), PARAMETER              :: ec_channel = modstr//'_ec'
  !
  ! POINTER TO TRACER META INFORAMTION
  TYPE(t_trinfo_list),             POINTER :: trlist_ec => NULL()
  !
  ! ARRAY OF TRACER INFO STRUCTS
  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti_ec => NULL()
  !
  ! NUMBER OF TRACERS
  INTEGER                                  :: ntrac_ec = 0
  !
  ! POINTER TO TRACER FIELDS
  ! - FOR INTERNAL MEMORY MANAGEMENT
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxt_ec   => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtte_ec => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtm1_ec => NULL()
  REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: pxtf_ec  => NULL()
  !
  ! - FOR OUTSIDE 'LOCAL LOOP'
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xt_ec    => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtte_ec  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtm1_ec  => NULL()
  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtf_ec   => NULL()
  !
!!$  ! Not Needed because CMAT accesses levels 1:ht_dim, i.e. upper part->no problem
!!$  ! - CMAT DOMAIN ACCESS ONLY: ecs=eCHAM cMAT subdomain
!!$  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xt_ecs   => NULL() 
!!$  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtte_ecs => NULL()
!!$  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtm1_ecs => NULL()
!!$  REAL(DP), DIMENSION(:,:,:,:),   POINTER :: xtf_ecs  => NULL()
!!$  !
!!$  ! Probably not needed because CMAT runs only outside local loop
!!$  ! - FOR INSIDE 'LOCAL LOOP'
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxt_ec   => NULL()
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtte_ec => NULL()
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtm1_ec => NULL()
!!$  REAL(DP), DIMENSION(:,:,:),     POINTER :: qxtf_ec  => NULL()
  !
  ! ####################################################################### 
! mz_ab_20090622-

! **********************************************************************+
END MODULE messy_main_tracer_mem_bi
! **********************************************************************+
