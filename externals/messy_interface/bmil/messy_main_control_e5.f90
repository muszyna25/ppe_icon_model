!==============================================================================
!
! This file controls the MESSy submodel CALLs
! Authors: Patrick Joeckel and Rolf Sander, MPICH, 2002-2007
!
!==============================================================================

SUBROUTINE messy_setup

  USE messy_main_switch   ! ONLY: USE_* 
  USE messy_main_switch_bi,  ONLY: main_switch_setup !initialize
  ! op_bk_20130828+
#ifdef __ICON__
  USE messy_main_mpi_bi,     ONLY: messy_mpi_initialize
#endif
  ! op_bk_20130828-
  ! op_bk_20130821
#ifndef __ICON__
  USE messy_main_qtimer_bi,  ONLY: main_qtimer_setup !initialize
! op_bk_20130821
#endif
  USE messy_main_timer_bi,   ONLY: main_timer_setup  !initialize
  USE messy_main_channel_bi, ONLY: messy_channel_init_restart

  IMPLICIT NONE

  ! op_bk_20130828+
#ifdef __ICON__
  CALL messy_mpi_initialize
#endif
  ! op_bk_20130828-
  ! TIMER SPLIT NEEDED HERE, TO ENABLE OTHER MODELS TO OVERWRITE NAMELIST
  ! ENTRY OF START DATE AND STOP DATE (e.g. MMD_CLIENT)
#ifdef MESSYTIMER
  CALL main_timer_setup(1)
!!$  CALL main_timer_setup(2) ! um_ak_20120307
  CALL messy_channel_init_restart
  CALL main_timer_setup(2)    ! um_ak_20120307
#endif
  ! do NOT change the following order of main_* subroutines
  ! op_bk_20130821
#ifndef __ICON__
  CALL main_qtimer_setup
! op_bk_20130821
#endif
  CALL main_switch_setup

  ! Decomposition setup is called directly from initialize.f90, 
  ! after the call to messy_setup. Reason: ngl and nlon are not yet
  ! known in messy_setup ...
  !!! CALL messy_main_decomp_setup 

END SUBROUTINE messy_setup

!==============================================================================

SUBROUTINE messy_initialize

  ! initialize submodels, read CTRL and CPL namelists (called by initialize)

  USE messy_main_switch   ! ONLY: USE_* 
  USE messy_main_channel_bi,ONLY: main_channel_initialize
  USE messy_main_timer_bi,  ONLY: main_timer_initialize
  USE messy_main_data_bi,   ONLY: main_data_initialize
  USE messy_main_tracer_bi, ONLY: main_tracer_initialize
  ! op_bk_20130820+
#ifndef __ICON__
  USE messy_main_tools_bi,  ONLY: main_tools_initialize
  USE messy_main_import_bi, ONLY: main_import_initialize
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi, ONLY: main_tendency_initialize
#endif

  USE messy_aeropt_e5,     ONLY:   aeropt_initialize
  USE messy_airsea_si,     ONLY:   airsea_initialize
  USE messy_bufly_e5,      ONLY:    bufly_initialize
  USE messy_cloud_e5,      ONLY:    cloud_initialize
  USE messy_convect_e5,    ONLY:  convect_initialize
  USE messy_cvtrans_si,    ONLY:  cvtrans_initialize
  USE messy_d14co_e5,      ONLY:    d14co_initialize
  USE messy_ddep_si,       ONLY:     ddep_initialize
  USE messy_dradon_si,     ONLY:   dradon_initialize
  USE messy_e4chem_si,     ONLY:   e4chem_initialize
  USE messy_gmxe_e5,       ONLY:     gmxe_initialize
  USE messy_gwave_si,      ONLY:    gwave_initialize
  USE messy_h2o_e5,        ONLY:      h2o_initialize
  USE messy_hetchem_e5,    ONLY:  hetchem_initialize 
  USE messy_jval_si,       ONLY:     jval_initialize
  USE messy_lnox_si,       ONLY:     lnox_initialize
  USE messy_m7_si,         ONLY:       m7_initialize
  USE messy_made_si,       ONLY:     made_initialize
  USE messy_mecca1_e5,     ONLY:   mecca1_initialize
  USE messy_mecca_si,      ONLY:    mecca_initialize
  USE messy_megan_si,      ONLY:    megan_initialize
  USE messy_mlocean_e5,    ONLY:  mlocean_initialize
  USE messy_mmforce_e5,    ONLY:  mmforce_initialize 
  USE messy_msbm_si,       ONLY:     msbm_initialize
  USE messy_o3orig_si,     ONLY:   o3orig_initialize
  USE messy_offemis_si,    ONLY:  offemis_initialize
  USE messy_onemis_si,     ONLY:   onemis_initialize
  USE messy_orbit_si,      ONLY:    orbit_initialize
  USE messy_photo_e5,      ONLY:    photo_initialize
  USE messy_plumegas_si,   ONLY: plumegas_initialize
  USE messy_psc_e5,        ONLY:      psc_initialize
  USE messy_ptrac_si,      ONLY:    ptrac_initialize
  USE messy_qbo_si,        ONLY:      qbo_initialize
  USE messy_rad4all_e5,    ONLY:  rad4all_initialize
  USE messy_scalc_si,      ONLY:    scalc_initialize
  USE messy_scav_si,       ONLY:     scav_initialize
  USE messy_satsims_e5,    ONLY:  satsims_initialize
  USE messy_scout_si,      ONLY:    scout_initialize
  USE messy_sedi_si,       ONLY:     sedi_initialize
  USE messy_s4d_si,        ONLY:      s4d_initialize
  USE messy_sorbit_si,     ONLY:   sorbit_initialize
  USE messy_spe_e5,        ONLY:      spe_initialize
  USE messy_spacenox_e5,   ONLY: spacenox_initialize
  USE messy_timepos_e5,    ONLY:  timepos_initialize
  USE messy_tnudge_si,     ONLY:   tnudge_initialize
  USE messy_trexp_si,      ONLY:    trexp_initialize
  USE messy_tropop_si,     ONLY:   tropop_initialize
  USE messy_vahr_e5,       ONLY:     vahr_initialize
  USE messy_viso_si,       ONLY:     viso_initialize
  ! op_bk_20130820+
#endif

  IMPLICIT NONE

  ! do NOT change the following order of main_* subroutines
#ifdef MESSYTIMER
  CALL main_timer_initialize
#endif
  CALL main_data_initialize
  CALL main_channel_initialize
  CALL main_tracer_initialize ! namelist(s)
  ! op_bk_20130820+
#ifndef __ICON__
  CALL main_tools_initialize
  CALL main_import_initialize
#ifdef MESSYTENDENCY
  CALL main_tendency_initialize
#endif  

  ! special requirements:
  !   (define number of lagrangian cells)
  ! - (the rest is in alphabetical order)
  IF (USE_AEROPT)   CALL   aeropt_initialize
  IF (USE_AIRSEA)   CALL   airsea_initialize
  IF (USE_BUFLY)    CALL    bufly_initialize
  IF (USE_CLOUD)    CALL    cloud_initialize
  IF (USE_CONVECT)  CALL  convect_initialize
  IF (USE_CVTRANS)  CALL  cvtrans_initialize
  IF (USE_D14CO)    CALL    d14co_initialize
  IF (USE_DDEP)     CALL     ddep_initialize
  IF (USE_DRADON)   CALL   dradon_initialize
  IF (USE_E4CHEM)   CALL   e4chem_initialize
  IF (USE_GMXE)     CALL     gmxe_initialize
  IF (USE_GWAVE)    CALL    gwave_initialize
  IF (USE_H2O)      CALL      h2o_initialize
  IF (USE_HETCHEM)  CALL  hetchem_initialize
  IF (USE_JVAL)     CALL     jval_initialize
  IF (USE_LNOX)     CALL     lnox_initialize
  IF (USE_M7)       CALL       m7_initialize
  IF (USE_MADE)     CALL     made_initialize
  IF (USE_MECCA1)   CALL   mecca1_initialize
  IF (USE_MECCA)    CALL    mecca_initialize
  IF (USE_MEGAN)    CALL    megan_initialize
  IF (USE_MLOCEAN)  CALL  mlocean_initialize
  IF (USE_MMFORCE)  CALL  mmforce_initialize
  IF (USE_MSBM)     CALL     msbm_initialize
  IF (USE_O3ORIG)   CALL   o3orig_initialize
  IF (USE_OFFEMIS)  CALL  offemis_initialize
  IF (USE_ONEMIS)   CALL   onemis_initialize
  IF (USE_ORBIT)    CALL    orbit_initialize
  IF (USE_PHOTO)    CALL    photo_initialize
  IF (USE_PSC)      CALL      psc_initialize
  IF (USE_PLUMEGAS) CALL plumegas_initialize
  IF (USE_PTRAC)    CALL    ptrac_initialize
  IF (USE_QBO)      CALL      qbo_initialize
  IF (USE_RAD4ALL)  CALL  rad4all_initialize
  IF (USE_SATSIMS)  CALL  satsims_initialize
  IF (USE_SCALC)    CALL    scalc_initialize
  IF (USE_SCAV)     CALL     scav_initialize
  IF (USE_SCOUT)    CALL    scout_initialize
  IF (USE_SEDI)     CALL     sedi_initialize
  IF (USE_S4D)      CALL      s4d_initialize
  IF (USE_SORBIT)   CALL   sorbit_initialize
  IF (USE_SPE)      CALL      spe_initialize
  IF (USE_SPACENOX) CALL spacenox_initialize
  IF (USE_TIMEPOS)  CALL  timepos_initialize
  IF (USE_TNUDGE)   CALL   tnudge_initialize
  IF (USE_TREXP)    CALL    trexp_initialize
  IF (USE_TROPOP)   CALL   tropop_initialize
  IF (USE_VAHR)     CALL     vahr_initialize
  IF (USE_VISO)     CALL     viso_initialize
  ! op_bk_20130820+
#endif

END SUBROUTINE messy_initialize

!==============================================================================

SUBROUTINE messy_new_tracer

  ! request tracers (called by initialize)

  USE messy_main_switch    ! ONLY: USE_*
  USE messy_main_tracer_bi,  ONLY: main_tracer_new_tracer

  ! op_bk_20130820+
#ifndef __ICON__
  USE messy_cloud_e5,      ONLY:   cloud_new_tracer
  USE messy_d14co_e5,      ONLY:   d14co_new_tracer
  USE messy_dradon_si,     ONLY:  dradon_new_tracer
  USE messy_e4chem_si,     ONLY:  e4chem_new_tracer
  USE messy_gmxe_e5,       ONLY:    gmxe_new_tracer
  USE messy_h2o_e5,        ONLY:     h2o_new_tracer
  USE messy_m7_si,         ONLY:      m7_new_tracer
  USE messy_made_si,       ONLY:    made_new_tracer
  USE messy_mecca1_e5,     ONLY:  mecca1_new_tracer
  USE messy_mecca_si,      ONLY:   mecca_new_tracer
  USE messy_msbm_si,       ONLY:    msbm_new_tracer
  USE messy_o3orig_si,     ONLY:  o3orig_new_tracer
  USE messy_plumegas_si,   ONLY: plumegas_new_tracer
  USE messy_psc_e5,        ONLY:     psc_new_tracer
  USE messy_ptrac_si,      ONLY:   ptrac_new_tracer
  USE messy_scav_si,       ONLY:    scav_new_tracer
  USE messy_trexp_si,      ONLY:   trexp_new_tracer
  ! op_bk_20130820+
#endif

  IMPLICIT NONE

  CALL main_tracer_new_tracer(1)  ! define tracer set

  ! op_bk_20130820+
#ifndef __ICON__
  ! special requirements:
  ! - h2o_new_tracer should be called before  gmXe_new_tracer
  ! - (the rest is in alphabetical order)
  IF (USE_MECCA1)  CALL  mecca1_new_tracer
  IF (USE_MECCA)   CALL   mecca_new_tracer
  IF (USE_CLOUD)   CALL   cloud_new_tracer
  IF (USE_D14CO)   CALL   d14co_new_tracer
  IF (USE_DRADON)  CALL  dradon_new_tracer
  IF (USE_E4CHEM)  CALL  e4chem_new_tracer
  IF (USE_H2O)     CALL     h2o_new_tracer
  IF (USE_GMXE)    CALL    gmxe_new_tracer
  IF (USE_M7)      CALL      m7_new_tracer
  IF (USE_MADE)    CALL    made_new_tracer
  IF (USE_MSBM)    CALL    msbm_new_tracer
  IF (USE_O3ORIG)  CALL  o3orig_new_tracer
  IF (USE_PLUMEGAS) CALL plumegas_new_tracer
  IF (USE_PSC)     CALL     psc_new_tracer
  IF (USE_PTRAC)   CALL   ptrac_new_tracer
  IF (USE_SCAV)    CALL    scav_new_tracer
  IF (USE_TREXP)   CALL   trexp_new_tracer
  ! op_bk_20130820+
#endif

  CALL main_tracer_new_tracer(2) ! define tracer families
  CALL main_tracer_new_tracer(3) ! diagnostic output

END SUBROUTINE messy_new_tracer

!==============================================================================

SUBROUTINE messy_init_memory

  ! allocate memory and define output channels
  ! (called by init_memory in mo_memory_streams)

  USE messy_main_switch    ! ONLY: USE_*
  ! op_bk_20130821
#ifndef __ICON__
  USE messy_main_qtimer_bi,  ONLY: main_qtimer_init_memory
! op_bk_20130821
#endif
  USE messy_main_data_bi,    ONLY: main_data_init_memory
  USE messy_main_tracer_bi,  ONLY: main_tracer_init_memory
  USE messy_main_channel_bi, ONLY: main_channel_init_memory
  ! op_bk_20130820+
#ifndef __ICON__
  USE messy_main_import_bi,  ONLY: main_import_init_memory
!!$#ifdef MESSYTENDENCY
!!$  USE messy_main_tendency_bi, ONLY: main_tendency_init_memory
!!$#endif

  USE messy_aeropt_e5,     ONLY:     aeropt_init_memory
  USE messy_airsea_si,     ONLY:     airsea_init_memory
  USE messy_cloud_e5,      ONLY:      cloud_init_memory
  USE messy_convect_e5,    ONLY:    convect_init_memory
  USE messy_cvtrans_si,    ONLY:    cvtrans_init_memory 
  USE messy_d14co_e5,      ONLY:      d14co_init_memory
  USE messy_ddep_si,       ONLY:       ddep_init_memory
  USE messy_dradon_si,     ONLY:     dradon_init_memory
  USE messy_e4chem_si,     ONLY:     e4chem_init_memory
  USE messy_gmxe_e5,       ONLY:       gmxe_init_memory
  USE messy_gwave_si,      ONLY:      gwave_init_memory
  USE messy_h2o_e5,        ONLY:        h2o_init_memory
  USE messy_hetchem_e5,    ONLY:    hetchem_init_memory 
  USE messy_jval_si,       ONLY:       jval_init_memory
  USE messy_lnox_si,       ONLY:       lnox_init_memory
  USE messy_m7_si,         ONLY:         m7_init_memory
  USE messy_made_si,       ONLY:       made_init_memory
  USE messy_mecca1_e5,     ONLY:     mecca1_init_memory
  USE messy_mecca_si,      ONLY:      mecca_init_memory
  USE messy_megan_si,      ONLY:      megan_init_memory
  USE messy_mlocean_e5,    ONLY:    mlocean_init_memory
  USE messy_mmforce_e5,    ONLY:    mmforce_init_memory
  USE messy_msbm_si,       ONLY:       msbm_init_memory
  USE messy_o3orig_si,     ONLY:     o3orig_init_memory
  USE messy_onemis_si,     ONLY:     onemis_init_memory
  USE messy_orbit_si,      ONLY:      orbit_init_memory
  USE messy_photo_e5,      ONLY:      photo_init_memory
  USE messy_plumegas_si,   ONLY:   plumegas_init_memory
  USE messy_psc_e5,        ONLY:        psc_init_memory
  USE messy_ptrac_si,      ONLY:      ptrac_init_memory
  USE messy_qbo_si,        ONLY:        qbo_init_memory
  USE messy_rad4all_e5,    ONLY:    rad4all_init_memory
  USE messy_s4d_si,        ONLY:        s4d_init_memory
  USE messy_satsims_e5,    ONLY:    satsims_init_memory
  USE messy_scav_si,       ONLY:       scav_init_memory
  USE messy_scout_si,      ONLY:      scout_init_memory
  USE messy_sedi_si,       ONLY:       sedi_init_memory
  USE messy_spe_e5,        ONLY:        spe_init_memory
  USE messy_spacenox_e5,   ONLY:   spacenox_init_memory
  USE messy_tnudge_si,     ONLY:     tnudge_init_memory
  USE messy_trexp_si,      ONLY:      trexp_init_memory
  USE messy_tropop_si,     ONLY:     tropop_init_memory
  ! op_bk_20130820+
#endif

  IMPLICIT NONE

  ! op_bk_20130821
#ifndef __ICON__
  CALL main_qtimer_init_memory
! op_bk_20130821
#endif
  ! associate ECHAM5 streams to MESSy data channels
  CALL main_channel_init_memory   ! setup tracer memory
  ! op_bk_20130903+
!   CALL main_tracer_init_memory(1) ! create data channel BML <-> SMIL
  ! op_bk_20130903-
  CALL main_data_init_memory
  ! op_bk_20130820+
#ifndef __ICON__
  CALL main_import_init_memory

  ! special requirements:
  ! - HAMOCC after MPIOM
  ! - (in alphabetical order)
  IF (USE_AEROPT)   CALL     aeropt_init_memory
  IF (USE_AIRSEA)   CALL     airsea_init_memory
  IF (USE_CLOUD)    CALL      cloud_init_memory
  IF (USE_CONVECT)  CALL    convect_init_memory
  IF (USE_CVTRANS)  CALL    cvtrans_init_memory
  IF (USE_D14CO)    CALL      d14co_init_memory
  IF (USE_DDEP)     CALL       ddep_init_memory
  IF (USE_DRADON)   CALL     dradon_init_memory
  IF (USE_E4CHEM)   CALL     e4chem_init_memory
  IF (USE_GMXE)     CALL       gmxe_init_memory
  IF (USE_GWAVE)    CALL      gwave_init_memory
  IF (USE_H2O)      CALL        h2o_init_memory
  IF (USE_HETCHEM)  CALL    hetchem_init_memory
  IF (USE_JVAL)     CALL       jval_init_memory
  IF (USE_LNOX)     CALL       lnox_init_memory
  IF (USE_M7)       CALL         m7_init_memory
  IF (USE_MADE)     CALL       made_init_memory
  IF (USE_MECCA1)   CALL     mecca1_init_memory
  IF (USE_MECCA)    CALL      mecca_init_memory
  IF (USE_MEGAN)    CALL      megan_init_memory
  IF (USE_MLOCEAN)  CALL    mlocean_init_memory
  IF (USE_MMFORCE)  CALL    mmforce_init_memory
  IF (USE_MSBM)     CALL       msbm_init_memory
  IF (USE_O3ORIG)   CALL     o3orig_init_memory
  IF (USE_ONEMIS)   CALL     onemis_init_memory
  IF (USE_ORBIT)    CALL      orbit_init_memory
  IF (USE_PHOTO)    CALL      photo_init_memory
  IF (USE_PLUMEGAS) CALL   plumegas_init_memory
  IF (USE_PSC)      CALL        psc_init_memory
  IF (USE_PTRAC)    CALL      ptrac_init_memory
  IF (USE_QBO)      CALL        qbo_init_memory
  IF (USE_RAD4ALL)  CALL    rad4all_init_memory
  IF (USE_S4D)      CALL        s4d_init_memory
  IF (USE_SATSIMS)  CALL    satsims_init_memory
  IF (USE_SCAV)     CALL       scav_init_memory
  IF (USE_SCOUT)    CALL      scout_init_memory
  IF (USE_SEDI)     CALL       sedi_init_memory
  IF (USE_SPE)      CALL        spe_init_memory
  IF (USE_SPACENOX) CALL   spacenox_init_memory
  IF (USE_TNUDGE)   CALL     tnudge_init_memory
  IF (USE_TREXP)    CALL      trexp_init_memory
  IF (USE_TROPOP)   CALL     tropop_init_memory
! op_bk_20130820+
#endif

  ! - associate tracer memory to MESSy channel(s)
  ! - setting meta information of family-members to fraction
  !   (for advection initialization)
  ! op_bk_20130903+
!   CALL main_tracer_init_memory(2)
  ! op_bk_20130903-
!!$#ifdef MESSYTENDENCY
!!$  ! - must be called after all submodels have placed their requests
!!$  !   to alter prognostic variables (and tracers!) and after the
!!$  !   number of tracers is fixated (i.e., after main_tracer_init_memory)
!!$  CALL main_tendency_init_memory
!!$#endif

END SUBROUTINE messy_init_memory

!==============================================================================

SUBROUTINE messy_init_coupling

  ! (will be called by control.f90)

  USE messy_main_switch ! ONLY: USE_*
  USE messy_main_channel_bi, ONLY: main_channel_init_coupling
  USE messy_main_tracer_bi,  ONLY: main_tracer_init_coupling
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi, ONLY: main_tendency_init_coupling
#endif

! op_bk_20130820+
#ifndef __ICON__
  USE messy_aeropt_e5,       ONLY:   aeropt_init_coupling
  USE messy_airsea_si,       ONLY:   airsea_init_coupling
  USE messy_bufly_e5,        ONLY:    bufly_init_coupling
  USE messy_cloud_e5,        ONLY:    cloud_init_coupling
  USE messy_convect_e5,      ONLY:  convect_init_coupling
  USE messy_cvtrans_si,      ONLY:  cvtrans_init_coupling
  USE messy_d14co_e5,        ONLY:    d14co_init_coupling
  USE messy_ddep_si,         ONLY:     ddep_init_coupling
  USE messy_e4chem_si,       ONLY:   e4chem_init_coupling
  USE messy_gmxe_e5,         ONLY:     gmxe_init_coupling
  USE messy_hetchem_e5,      ONLY:  hetchem_init_coupling
  USE messy_h2o_e5,          ONLY:      h2o_init_coupling
  USE messy_jval_si,         ONLY:     jval_init_coupling
  USE messy_lnox_si,         ONLY:     lnox_init_coupling
  USE messy_m7_si,           ONLY:       m7_init_coupling
  USE messy_made_si,         ONLY:     made_init_coupling
  USE messy_mecca1_e5,       ONLY:   mecca1_init_coupling
  USE messy_mecca_si,        ONLY:    mecca_init_coupling
  USE messy_megan_si,        ONLY:    megan_init_coupling
  USE messy_mlocean_e5,      ONLY:  mlocean_init_coupling
  USE messy_mmforce_e5,      ONLY:  mmforce_init_coupling
  USE messy_msbm_si,         ONLY:     msbm_init_coupling
  USE messy_o3orig_si,       ONLY:   o3orig_init_coupling
  USE messy_offemis_si,      ONLY:  offemis_init_coupling
  USE messy_onemis_si,       ONLY:   onemis_init_coupling
  USE messy_orbit_si,        ONLY:    orbit_init_coupling
  USE messy_photo_e5,        ONLY:    photo_init_coupling
  USE messy_plumegas_si,     ONLY: plumegas_init_coupling
  USE messy_psc_e5,          ONLY:      psc_init_coupling
  USE messy_qbo_si,          ONLY:      qbo_init_coupling
  USE messy_rad4all_e5,      ONLY:  rad4all_init_coupling
  USE messy_s4d_si,          ONLY:      s4d_init_coupling
  USE messy_sorbit_si,       ONLY:   sorbit_init_coupling
  USE messy_satsims_e5,      ONLY:  satsims_init_coupling
  USE messy_scalc_si,        ONLY:    scalc_init_coupling
  USE messy_scav_si,         ONLY:     scav_init_coupling
  USE messy_scout_si,        ONLY:    scout_init_coupling
  USE messy_sedi_si,         ONLY:     sedi_init_coupling
  USE messy_spe_e5,          ONLY:      spe_init_coupling
  USE messy_spacenox_e5,     ONLY: spacenox_init_coupling
  USE messy_timepos_e5,      ONLY:  timepos_init_coupling
  USE messy_tnudge_si,       ONLY:   tnudge_init_coupling
  USE messy_trexp_si,        ONLY:    trexp_init_coupling
  USE messy_vahr_e5,         ONLY:     vahr_init_coupling
  USE messy_viso_si,         ONLY:     viso_init_coupling
! op_bk_20130820
#endif

  IMPLICIT NONE

  ! - resetting meta information of family-members (to tracers)
  !   (after advection initialisation)
  CALL main_tracer_init_coupling

  ! special requirements:
  ! - Some SMs (as indicated) need to define new channel objects in
  !   init_coupling. Therefore, the diagnostic SMs S4D, SORBIT, SCOUT, SCALC,
  !   VISO, LGGP, and LGVFLUX 
  !   must be called at the end, since all other channel
  !   objects need to be present.
  ! - VISO must be called before LGVFLUX, because LGVFLUX requires the
  !   iso-surfaces defined in VISO
  ! - LGGP must be called before LGVFLUX
  ! - LGGP must be called twice (2nd time after  LGVFLUX) in order to be 
  !   able to convert Lagrangian information from LGVFLUX
  ! - GMXE must be called after offemis and onemis, since it uses objects
  !   from both submodels
  ! - ORBIT must be called before RADIATION, since it creates new
  !   objects for RADIATION
  ! - (the rest is in alphabetical order)

! op_bk_20130820
#ifndef __ICON__
  IF (USE_AEROPT)   CALL   aeropt_init_coupling
  IF (USE_AIRSEA)   CALL   airsea_init_coupling
  IF (USE_BUFLY)    CALL    bufly_init_coupling
  IF (USE_CLOUD)    CALL    cloud_init_coupling
  IF (USE_CONVECT)  CALL  convect_init_coupling
  IF (USE_CVTRANS)  CALL  cvtrans_init_coupling
  IF (USE_D14CO)    CALL    d14co_init_coupling
  IF (USE_DDEP)     CALL     ddep_init_coupling ! new channel objects
  IF (USE_E4CHEM)   CALL   e4chem_init_coupling
  IF (USE_H2O)      CALL      h2o_init_coupling
  IF (USE_HETCHEM)  CALL  hetchem_init_coupling
  IF (USE_JVAL)     CALL     jval_init_coupling
  IF (USE_LNOX)     CALL     lnox_init_coupling
  IF (USE_M7)       CALL       m7_init_coupling
  IF (USE_MADE)     CALL     made_init_coupling
  IF (USE_MECCA1)   CALL   mecca1_init_coupling
  IF (USE_MECCA)    CALL    mecca_init_coupling
  IF (USE_MEGAN)    CALL    megan_init_coupling
  IF (USE_MLOCEAN)  CALL  mlocean_init_coupling
  IF (USE_MMFORCE)  CALL  mmforce_init_coupling
  IF (USE_MSBM)     CALL     msbm_init_coupling
  IF (USE_O3ORIG)   CALL   o3orig_init_coupling
  IF (USE_OFFEMIS)  CALL  offemis_init_coupling ! new channel objects
  IF (USE_ONEMIS)   CALL   onemis_init_coupling
  IF (USE_ORBIT)    CALL    orbit_init_coupling ! new channel objects
  IF (USE_GMXE)     CALL     gmxe_init_coupling
  IF (USE_PHOTO)    CALL    photo_init_coupling
  IF (USE_PLUMEGAS) CALL plumegas_init_coupling
  IF (USE_PSC)      CALL      psc_init_coupling
  IF (USE_QBO)      CALL      qbo_init_coupling
  IF (USE_RAD4ALL)  CALL  rad4all_init_coupling
  IF (USE_SATSIMS)  CALL  satsims_init_coupling
  IF (USE_SCAV)     CALL     scav_init_coupling
  IF (USE_SEDI)     CALL     sedi_init_coupling ! new channel objects
  IF (USE_TIMEPOS)  CALL  timepos_init_coupling
  IF (USE_TNUDGE)   CALL   tnudge_init_coupling
  IF (USE_TREXP)    CALL    trexp_init_coupling
  IF (USE_SPE)      CALL      spe_init_coupling
  IF (USE_SPACENOX) CALL spacenox_init_coupling
  IF (USE_VAHR)     CALL     vahr_init_coupling
  IF (USE_S4D)      CALL      s4d_init_coupling ! new channel objects
  IF (USE_SCOUT)    CALL    scout_init_coupling ! new channel objects
  IF (USE_VISO)     CALL     viso_init_coupling ! new channel objects
  IF (USE_SORBIT)   CALL   sorbit_init_coupling ! new channel objects
  IF (USE_SCALC)    CALL    scalc_init_coupling ! new channel objects
! op_bk_20130820
#endif

#ifdef MESSYTENDENCY
  ! - must be called after all submodels have placed their requests
  !   to alter prognostic variables (and tracers!) and after the
  !   number of tracers is fixated (i.e., after main_tracer_init_memory)
  CALL main_tendency_init_coupling
#endif

  CALL main_channel_init_coupling

END SUBROUTINE messy_init_coupling

!==============================================================================

SUBROUTINE messy_init_tracer

  ! tracer initialization (called by xtini in mo_tracer.f90)

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_tracer_bi, ONLY: main_tracer_init_tracer

! op_bk_20130820
#ifndef __ICON__
!!$  USE messy_d14co_e5,      ONLY:   d14co_init_tracer
!!$  USE messy_dradon_si,     ONLY:  dradon_init_tracer
  USE messy_e4chem_si,     ONLY:  e4chem_init_tracer
  USE messy_h2o_e5,        ONLY:     h2o_init_tracer
  USE messy_mecca1_e5,     ONLY:  mecca1_init_tracer
  USE messy_mecca_si,      ONLY:   mecca_init_tracer
!!$  USE messy_msbm_si,       ONLY:    msbm_init_tracer
!!$  USE messy_psc_e5,        ONLY:     psc_init_tracer
!!$  USE messy_ptrac_si,      ONLY:   ptrac_init_tracer
! op_bk_20130820
#endif

  IMPLICIT NONE

  CALL main_tracer_init_tracer(1) ! check tracer init from rerun

  ! special requirements:
  !   main_tracer_init_tracer(4), to ensure that all fields for
  !   tracer initialisation (mass conserving transformation between gridpoint
  !   and Lagrangian space) are present

  CALL main_tracer_init_tracer(4) ! initialize via tracer_init (tracer.nml)

! op_bk_20130820
#ifndef __ICON__
  ! CHECK / MODIFY INDIVIDUALLY
  ! special requirements:
  ! - (the rest is in alphabetical order)
!!$  IF (USE_D14CO)   CALL   d14co_init_tracer    ! lagrange !
!!$  IF (USE_DRADON)  CALL  dradon_init_tracer    ! lagrange !
  IF (USE_E4CHEM)  CALL  e4chem_init_tracer
  IF (USE_H2O)     CALL     h2o_init_tracer    ! lagrange !
  IF (USE_MECCA1)  CALL  mecca1_init_tracer
  IF (USE_MECCA)   CALL   mecca_init_tracer    ! lagrange !
!!$  IF (USE_MSBM)    CALL    msbm_init_tracer
!!$  IF (USE_PSC)     CALL     psc_init_tracer
!!$  IF (USE_PTRAC)   CALL   ptrac_init_tracer    ! lagrange !
! op_bk_20130820
#endif
  CALL main_tracer_init_tracer(2) ! check tracer init by init_tracer
                                  ! set to constant if not yet initialized
  CALL main_tracer_init_tracer(3) ! diagnose tracer initialization

END SUBROUTINE messy_init_tracer

!==============================================================================

SUBROUTINE messy_global_start

  ! (called by stepon)

  USE messy_main_switch    ! ONLY: USE_*
  USE messy_main_timer_bi,   ONLY: main_timer_global_start
  USE messy_main_data_bi,    ONLY: main_data_global_start
  USE messy_main_tracer_bi,  ONLY: main_tracer_global_start
! op_bk_20130820
#ifndef __ICON__
  USE messy_main_import_bi,  ONLY: main_import_global_start
!!$  USE messy_main_channel_bi, ONLY: main_channel_global_start

  USE messy_bufly_e5,      ONLY:    bufly_global_start
  USE messy_d14co_e5,      ONLY:    d14co_global_start
  USE messy_ddep_si,       ONLY:     ddep_global_start
  USE messy_dradon_si,     ONLY:   dradon_global_start
  USE messy_gmxe_e5,       ONLY:     gmxe_global_start
  USE messy_h2o_e5,        ONLY:      h2o_global_start
  USE messy_hetchem_e5,    ONLY:  hetchem_global_start
  USE messy_jval_si,       ONLY:     jval_global_start
  USE messy_mlocean_e5,    ONLY:  mlocean_global_start
  USE messy_mmforce_e5,    ONLY:  mmforce_global_start
  USE messy_onemis_si,     ONLY:   onemis_global_start
  USE messy_orbit_si,      ONLY:    orbit_global_start
  USE messy_photo_e5,      ONLY:    photo_global_start
  USE messy_psc_e5,        ONLY:      psc_global_start
  USE messy_rad4all_e5,    ONLY:  rad4all_global_start
  USE messy_s4d_si,        ONLY:      s4d_global_start
  USE messy_sorbit_si,     ONLY:   sorbit_global_start
  USE messy_spe_e5,        ONLY:      spe_global_start
  USE messy_timepos_e5,    ONLY:  timepos_global_start
  USE messy_tnudge_si,     ONLY:   tnudge_global_start
  USE messy_trexp_si,      ONLY:    trexp_global_start
! op_bk_20130820
#endif

  IMPLICIT NONE

#ifdef MESSYTIMER
  CALL main_timer_global_start 
#endif
  CALL main_data_global_start
  CALL main_tracer_global_start
!!$  CALL main_channel_global_start
! op_bk_20130820
#ifndef __ICON__
  CALL main_import_global_start

  ! special requirements:
  ! - bufly_global_start must be called very first, because it modifies the
  !   initial condition
  !   mass conserving transformation between gridpoint and
  !   Lagrangian space are present
  ! - jval_global_start must be called after rad4all_global_start
  !   (cdisse must be up to date)
  ! - RAD4ALL must be called before ORBIT to provide the correct time offset
  !   Note: The new RADIATION scheme needs to encapsulate the call to ORBIT:
  !         1. RADIATON calculates trigger and time offset
  !         2. ORBIT calculates orbit parameters with time offset
  !         3. RADIATION calculates radiation
  !         Currently, not coupling from ORBIT back to RAD4ALL is implemented.
  ! - (the rest is in alphabetical order)
  IF (USE_BUFLY)    CALL    bufly_global_start
  IF (USE_D14CO)    CALL    d14co_global_start    ! lagrange !
  IF (USE_DDEP)     CALL     ddep_global_start    ! lagrange !
  IF (USE_DRADON)   CALL   dradon_global_start    ! lagrange !
  IF (USE_GMXE)     CALL     gmxe_global_start
  IF (USE_H2O)      CALL      h2o_global_start    ! lagrange !
  IF (USE_HETCHEM)  CALL  hetchem_global_start
  IF (USE_MLOCEAN)  CALL  mlocean_global_start
  IF (USE_MMFORCE)  CALL  mmforce_global_start
  IF (USE_ONEMIS)   CALL   onemis_global_start    ! lagrange !
  IF (USE_PHOTO)    CALL    photo_global_start
  IF (USE_PSC)      CALL      psc_global_start
  IF (USE_RAD4ALL)  CALL  rad4all_global_start
  IF (USE_ORBIT)    CALL    orbit_global_start
  IF (USE_JVAL)     CALL     jval_global_start    ! lagrange !
  IF (USE_S4D)      CALL      s4d_global_start
  IF (USE_SORBIT)   CALL   sorbit_global_start
  IF (USE_SPE)      CALL      spe_global_start
  IF (USE_TIMEPOS)  CALL  timepos_global_start
  IF (USE_TNUDGE)   CALL   tnudge_global_start
  IF (USE_TREXP)    CALL    trexp_global_start
! op_bk_20130820
#endif

  ! um_ak_20080709+
  ! renamed to messy_tracer_beforeadv
  ! for consistency reasons directly called from BML
!!$  CALL main_tracer_global_start(2) ! FAMILY: T2F, SUM
  ! um_ak_20080709-

END SUBROUTINE messy_global_start

!==============================================================================
! scan1.f90: ADEVECTION (scan1.f90)
! scan1.f90: CALL main_tracer_afteradv
!==============================================================================

SUBROUTINE messy_local_start

  ! (will be called by scan1)

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_data_bi,   ONLY: main_data_local_start
  USE messy_main_tracer_bi, ONLY: main_tracer_local_start

! op_bk_20130820
#ifndef __ICON__
  USE messy_msbm_si,        ONLY: msbm_local_start
  USE messy_psc_e5,         ONLY:  psc_local_start
  USE messy_e4chem_si,      ONLY: e4chem_local_start
! op_bk_20130820
#endif

  IMPLICIT NONE

  CALL main_data_local_start
  CALL main_tracer_local_start

! op_bk_20130820
#ifndef __ICON__
  ! special requirements:
  ! - call psc_local_start before modules which need flt_pscreg
  ! - call msbm_local_start before modules which need flt_pscreg
  IF (USE_MSBM)    CALL msbm_local_start
  IF (USE_PSC)     CALL  psc_local_start
  IF (USE_E4CHEM)  CALL e4chem_local_start
! op_bk_20130820
#endif

END SUBROUTINE messy_local_start

!==============================================================================

SUBROUTINE messy_radiation

  USE messy_main_switch   ! ONLY: USE_*

! op_bk_20130820
#ifndef __ICON__
  USE messy_aeropt_e5,     ONLY:  aeropt_radiation
  USE messy_cloud_e5,      ONLY:   cloud_radiation
  USE messy_gmxe_e5,       ONLY:    gmxe_radiation
  USE messy_rad4all_e5,    ONLY: rad4all_radiation
! op_bk_20130820
#endif

  IMPLICIT NONE

! op_bk_20130820
#ifndef __ICON__
  ! special requirements:
  ! - aeropt should be called directly before radiation
  !   such that no other process can influence the aerosol
  !   distribution after the optical properties have been determined
  IF (USE_CLOUD)   CALL   cloud_radiation
  IF (USE_GMXE)    CALL    gmxe_radiation
  IF (USE_AEROPT)  CALL  aeropt_radiation
  IF (USE_RAD4ALL) CALL rad4all_radiation
! op_bk_20130820
#endif

END SUBROUTINE messy_radiation

!==============================================================================

SUBROUTINE messy_vdiff

  ! (called by vdiff)

  USE messy_main_switch ! ONLY: USE_*

! op_bk_20130820
#ifndef __ICON__
  USE messy_airsea_si,         ONLY:     airsea_vdiff
  USE messy_ddep_si,           ONLY:       ddep_vdiff
  USE messy_dradon_si,         ONLY:     dradon_vdiff
  USE messy_gmxe_e5,           ONLY:       gmxe_vdiff
  USE messy_m7_si,             ONLY:         m7_vdiff
  USE messy_made_si,           ONLY:       made_vdiff
  USE messy_mecca1_e5,         ONLY:     mecca1_vdiff
  USE messy_mecca_si,          ONLY:      mecca_vdiff
  USE messy_megan_si,          ONLY:      megan_vdiff
  USE messy_offemis_si,        ONLY:    offemis_vdiff
  USE messy_onemis_si,         ONLY:     onemis_vdiff
  USE messy_tropop_si,         ONLY:     tropop_vdiff
  USE messy_viso_si,           ONLY:       viso_vdiff
! op_bk_20130820
#endif

  IMPLICIT NONE

! op_bk_20130820
#ifndef __ICON__
  ! special requirements:

  ! - call mecca_vdiff(1) after emissions (ONLEM, OFFLEM, OFFEMIS, AIRSEA,
  !   MEGAN)
  !   and before deposition (DRYDEP) to get the "real" new emitted aerosol
  ! - call mecca_vdiff(2) after DRYDEP
  ! - call m7_vdiff and gmxe_vdiff after emissions 
  !   (ONLEM, OFFLEM, OFFEMIS, AIRSEA, MEGAN) to get the
  !   current emissions if m7/gmxe calculates its own emissions
  !   to get the
  !   current emissions if cam does not calculate its own emissions
  ! - tropop_vdiff must be called before viso_vdiff
  ! - viso_vdiff and tropop_vdiff
  !   must be called first, since they calulate diagnostic quantities
  !   needed in other submodels

  IF (USE_TROPOP)  CALL  tropop_vdiff
  IF (USE_VISO)    CALL    viso_vdiff
  IF (USE_DRADON)  CALL  dradon_vdiff    ! lagrange !
  IF (USE_OFFEMIS) CALL offemis_vdiff  
  IF (USE_ONEMIS)  CALL  onemis_vdiff    ! lagrange !
  IF (USE_AIRSEA)  CALL  airsea_vdiff    ! lagrange !
  IF (USE_MEGAN)   CALL   megan_vdiff    
  IF (USE_M7)      CALL      m7_vdiff
  IF (USE_MADE)    CALL    made_vdiff
  IF (USE_GMXE)    CALL    gmxe_vdiff
  IF (USE_MECCA1)  CALL  mecca1_vdiff(1)
  IF (USE_MECCA)   CALL   mecca_vdiff(1)
  IF (USE_DDEP)    CALL    ddep_vdiff    ! lagrange !
  IF (USE_MECCA1)  CALL  mecca1_vdiff(2)
  IF (USE_MECCA)   CALL   mecca_vdiff(2)
! op_bk_20130820
#endif

END SUBROUTINE messy_vdiff

!==============================================================================

SUBROUTINE messy_radheat

  USE messy_main_switch   ! ONLY: USE_*

! op_bk_20130820
#ifndef __ICON__
  USE messy_rad4all_e5,    ONLY: rad4all_radheat
! op_bk_20130820
#endif

  IMPLICIT NONE

! op_bk_20130820
#ifndef __ICON__
  ! special requirements:
  IF (USE_RAD4ALL) CALL rad4all_radheat
! op_bk_20130820
#endif

END SUBROUTINE messy_radheat

!==============================================================================

SUBROUTINE messy_gwdrag

  USE messy_main_switch   ! ONLY: USE_*

! op_bk_20130820
#ifndef __ICON__
  USE messy_gwave_si,       ONLY: gwave_physc
! op_bk_20130820
#endif

  IMPLICIT NONE

! op_bk_20130820
#ifndef __ICON__
  ! special requirements:
   IF (USE_GWAVE) CALL gwave_physc
! op_bk_20130820
#endif

END SUBROUTINE messy_gwdrag

!==============================================================================

SUBROUTINE messy_convec

  ! (called by physc before convection subroutine cucall)

  USE messy_main_switch ! ONLY: USE_*

! op_bk_20130820
#ifndef __ICON__
  USE messy_convect_e5,     ONLY:  convect_convec
  USE messy_cvtrans_si,     ONLY:  cvtrans_convec
  USE messy_scav_si,        ONLY:     scav_convec
  USE messy_cloud_e5,       ONLY:    cloud_convec
! op_bk_20130820
#endif

  IMPLICIT NONE

! op_bk_20130820
#ifndef __ICON__
  ! special requirements:
  ! - convect_convec must be called first because it calculates
  !   the convection and all other submodels use convection data
  ! - cvtrans_convec(use_scav,1), scav_convec, and cvtrans_convec(use_scav,2)
  !   must be in this order and together
  IF (USE_CONVECT) CALL  convect_convec
  IF (USE_CVTRANS) CALL  cvtrans_convec(1)
  IF (USE_SCAV)    CALL     scav_convec
  IF (USE_CVTRANS) CALL  cvtrans_convec(2)
  IF (USE_CLOUD)   CALL    cloud_convec
! op_bk_20130820
#endif

END SUBROUTINE messy_convec

!==============================================================================

SUBROUTINE messy_mixlo

  USE messy_main_switch   ! ONLY: USE_*
! op_bk_20130820
#ifndef __ICON__
  USE messy_mlocean_e5,    ONLY: mlocean_mixlo

  ! special requirements:
  IF (USE_MLOCEAN) CALL mlocean_mixlo
! op_bk_20130820
#endif

END SUBROUTINE messy_mixlo

!==============================================================================

SUBROUTINE messy_physc

! op_bk_20130820
#ifndef __ICON__
  ! main submodel call (called by xtdriver in mo_tracer)

  ! ONLY FOR DEBUGGING +
  USE messy_main_data_bi,   ONLY: jrow, ngpblks, kproma 
#ifndef MESSYTIMER
  USE messy_main_data_bi,   ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
#else
  USE messy_main_timer,     ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
#endif
  USE messy_main_mpi_bi,               ONLY: p_pe
  ! ONLY FOR DEBUGGING -

  USE messy_main_switch ! ONLY: USE_*

  USE messy_bufly_e5,       ONLY:    bufly_physc
  USE messy_d14co_e5,       ONLY:    d14co_physc
  USE messy_dradon_si,      ONLY:   dradon_physc
  USE messy_e4chem_si,      ONLY:   e4chem_physc
  USE messy_gmxe_e5,        ONLY:     gmxe_physc
  USE messy_h2o_e5,         ONLY:      h2o_physc
  USE messy_hetchem_e5,     ONLY:  hetchem_physc
  USE messy_jval_si,        ONLY:     jval_physc
  USE messy_lnox_si,        ONLY:     lnox_physc
  USE messy_m7_si,          ONLY:       m7_physc
  USE messy_made_si,        ONLY:     made_physc
  USE messy_mecca1_e5,      ONLY:   mecca1_physc
  USE messy_mecca_si,       ONLY:    mecca_physc
  USE messy_mmforce_e5,     ONLY:  mmforce_physc
  USE messy_msbm_si,        ONLY:     msbm_physc
  USE messy_o3orig_si,      ONLY:   o3orig_physc
  USE messy_photo_e5,       ONLY:    photo_physc
  USE messy_plumegas_si,    ONLY: plumegas_physc
  USE messy_psc_e5,         ONLY:      psc_physc
  USE messy_qbo_si,         ONLY:      qbo_physc 
  USE messy_satsims_e5,     ONLY:  satsims_physc
  USE messy_scav_si,        ONLY:     scav_physc
  USE messy_sedi_si,        ONLY:     sedi_physc
  USE messy_spe_e5,         ONLY:      spe_physc
  USE messy_spacenox_e5,    ONLY: spacenox_physc
  USE messy_trexp_si,       ONLY:    trexp_physc
  USE messy_vahr_e5,        ONLY:     vahr_physc

  IMPLICIT NONE
  
  ! ONLY FOR DEBUGGING +
  IF (L_TIME_INFO) &
  WRITE(*, '(1x,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a3'//&
       &',a4,i4,a8,i5,a3,i5,a10,i5)') &
       YEAR,'-',MONTH,'-',DAY,' ',HOUR,':', MINUTE,':',SECOND,' # ' &
       ,' PE=',p_pe,' # JROW=',jrow,' / ',ngpblks,' # KPROMA=',kproma
  ! ONLY FOR DEBUGGING -

  ! special requirements: 
  ! - bufly_physc must be called first, because it modifies the initial cond.
  ! - jval_physc must be called before mecca_physc
  ! - photo_physc must be called before mecca_physc
  ! - psc_physc must be called before mecca_physc
  ! - msbm_physc must be called before mecca_physc
  ! - scav_physc(1) must be called before mecca_physc, m7_physc
  ! - scav_physc(2) must be called last, since species in the cloud phases are 
  !   released back into the gas phase
  ! - sedi_physc must be called before mecca_physc
  ! - m7_physc and gmxe_physc should be called before mecca_physc
  ! - hetchem_physc must be called:
  !   - after all aerosol modules (PSC, MSBM, M7, PTRAC)
  !   - directly before mecca_physc, to avoid problems due to the conversion
  !     of first order reaction rates in second order.
  ! - d14co_physc and trexp_physc must be called at the end, i.e., after all
  !   other chemistry
  ! - o3orig_physc must encapsulate mecca1_physc / e4chem_physc / mecca_physc

  IF (USE_BUFLY)    CALL    bufly_physc
  IF (USE_QBO)      CALL      qbo_physc
  IF (USE_SPE)      CALL      spe_physc
  IF (USE_SPACENOX) CALL spacenox_physc
  IF (USE_DRADON)   CALL   dradon_physc    ! lagrange !
  IF (USE_H2O)      CALL      h2o_physc    ! lagrange !
  IF (USE_LNOX)     CALL     lnox_physc    ! lagrange !
  IF (USE_PHOTO)    CALL    photo_physc
  IF (USE_JVAL)     CALL     jval_physc    ! lagrange !
  IF (USE_PSC)      CALL      psc_physc
  IF (USE_MSBM)     CALL     msbm_physc
  IF (USE_SCAV)     CALL     scav_physc(1)
  IF (USE_SEDI)     CALL     sedi_physc
  IF (USE_M7)       CALL       m7_physc
  IF (USE_MADE)     CALL     made_physc
  IF (USE_GMXE)     CALL     gmxe_physc
  IF (USE_HETCHEM)  CALL  hetchem_physc
  IF (USE_O3ORIG)   CALL   o3orig_physc(1)
  IF (USE_MECCA1)   CALL   mecca1_physc
  IF (USE_E4CHEM)   CALL   e4chem_physc
  IF (USE_MECCA)    CALL    mecca_physc    ! lagrange !
  IF (USE_O3ORIG)   CALL   o3orig_physc(2)
  IF (USE_MMFORCE)  CALL  mmforce_physc
  IF (USE_TREXP)    CALL    trexp_physc
  IF (USE_D14CO)    CALL    d14co_physc    ! lagrange !
  IF (USE_SATSIMS)  CALL  satsims_physc
  IF (USE_SCAV)     CALL     scav_physc(2)
  IF (USE_VAHR)     CALL     vahr_physc
  IF (USE_PLUMEGAS) CALL plumegas_physc
! op_bk_20130820
#endif

END SUBROUTINE messy_physc

!==============================================================================

SUBROUTINE messy_local_end

  ! diagnostics at the end of a time step (called by xtdiagn in mo_tracer)

  USE messy_main_switch   ! ONLY: USE_*

! op_bk_20130820
#ifndef __ICON__
  USE messy_h2o_e5,          ONLY:     h2o_local_end
  USE messy_hetchem_e5,      ONLY: hetchem_local_end
  USE messy_m7_si,           ONLY:      m7_local_end ! um_gg_20091030
  USE messy_msbm_si,         ONLY:    msbm_local_end
  USE messy_psc_e5,          ONLY:     psc_local_end
  USE messy_tnudge_si,       ONLY:  tnudge_local_end
  USE messy_viso_si,         ONLY:    viso_local_end
! op_bk_20130820
#endif

  IMPLICIT NONE

! op_bk_20130820
#ifndef __ICON__
  ! special requirements:
  ! - tnudge_local_end and viso_local_end must be called last, 
  !   since all channel objects must be up to date
  ! - (the rest is in alphabetical order)
  IF (USE_H2O)     CALL     h2o_local_end    ! lagrange !
  IF (USE_HETCHEM) CALL hetchem_local_end
  IF (USE_MSBM)    CALL    msbm_local_end
  IF (USE_M7)      CALL      m7_local_end
  IF (USE_PSC)     CALL     psc_local_end
  IF (USE_TNUDGE)  CALL  tnudge_local_end
  IF (USE_VISO)    CALL    viso_local_end
! op_bk_20130820
#endif

END SUBROUTINE messy_local_end

!==============================================================================

SUBROUTINE messy_global_end

  ! diagnostics at the end of a time step; outside the region loop
  ! (called by stepon)

  USE messy_main_switch    ! ONLY: USE_*
  ! op_bk_20130821
#ifndef __ICON__
  USE messy_main_qtimer_bi,  ONLY: main_qtimer_global_end
! op_bk_20130821
#endif
  USE messy_main_tracer_bi,  ONLY: main_tracer_global_end
!!$  USE messy_main_channel_bi, ONLY: main_channel_global_end
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,ONLY: main_tendency_global_end
#endif

! op_bk_20130820
#ifndef __ICON__
  USE messy_airsea_si,      ONLY:  airsea_global_end
  USE messy_cvtrans_si,     ONLY: cvtrans_global_end
  USE messy_ddep_si,        ONLY:    ddep_global_end
  USE messy_dradon_si,      ONLY:  dradon_global_end
  USE messy_d14co_e5,       ONLY:   d14co_global_end
  USE messy_gwave_si,       ONLY:   gwave_global_end
  USE messy_h2o_e5,         ONLY:     h2o_global_end
  USE messy_jval_si,        ONLY:    jval_global_end
  USE messy_lnox_si,        ONLY:    lnox_global_end
  USE messy_mecca_si,       ONLY:   mecca_global_end
  USE messy_megan_si,       ONLY:   megan_global_end
  USE messy_offemis_si,     ONLY: offemis_global_end
  USE messy_onemis_si,      ONLY:  onemis_global_end
  USE messy_scalc_si,       ONLY:   scalc_global_end
  USE messy_scav_si,        ONLY:    scav_global_end
  USE messy_sedi_si,        ONLY:    sedi_global_end
  USE messy_tnudge_si,      ONLY:  tnudge_global_end
  USE messy_trexp_si,       ONLY:   trexp_global_end
! op_bk_20130820
#endif

  IMPLICIT NONE

! op_bk_20130820
#ifndef __ICON__
  ! special requirements:
  !   since positions need to be updated
  ! - d14co_global_end and trexp_global_end must be called after
  !   all other lg-chemistry
  ! - jval_global_end must be called before mecca_global_end
  ! - lnox_global_end must be called before mecca_global_end
  ! - h2o_global_end must be called after mecca_global_end
  !   order) of all lg-submodels, since fields need to be up to date
  ! - mecca_global_end must be called after emissions (OFFLEM, OFFEMIS, ONLEM,
  !   AIRSEA, MEGAN) and before depositions (DRYDEP, DDEP) ; see vdiff for GP
  ! - scav_global_end(1) should be called before mecca_global_end to be 
  !   consistent with grid point order 
  ! - scav_global_end(2) should be called after all tracer modifying submodels
  ! - scalc_global_end must be called after depositions (DRYDEP, DDEP) and
  IF (USE_DRADON)  CALL  dradon_global_end    ! lagrange !
  IF (USE_OFFEMIS) CALL offemis_global_end    ! lagrange !
  IF (USE_ONEMIS)  CALL  onemis_global_end    ! lagrange !
  IF (USE_AIRSEA)  CALL  airsea_global_end    ! lagrange !
  IF (USE_MEGAN)   CALL   megan_global_end    
  IF (USE_GWAVE)   CALL   gwave_global_end
  IF (USE_CVTRANS) CALL cvtrans_global_end    ! lagrange !
  IF (USE_SCAV)    CALL    scav_global_end(1) ! lagrange !
  IF (USE_SEDI)    CALL    sedi_global_end    ! lagrange !
  IF (USE_JVAL)    CALL    jval_global_end    ! lagrange !
  IF (USE_LNOX)    CALL    lnox_global_end    ! lagrange !
  IF (USE_MECCA)   CALL   mecca_global_end    ! lagrange !
  IF (USE_DDEP)    CALL    ddep_global_end    ! lagrange !
  IF (USE_H2O)     CALL     h2o_global_end    ! lagrange !
  IF (USE_D14CO)   CALL   d14co_global_end    ! lagrange !
  IF (USE_SCAV)    CALL    scav_global_end(2) ! lagrange !
  IF (USE_SCALC)   CALL   scalc_global_end
  IF (USE_TNUDGE)  CALL  tnudge_global_end    ! lagrange !
  IF (USE_TREXP)   CALL   trexp_global_end    ! lagrange !
! op_bk_20130820
#endif

  CALL main_tracer_global_end
#ifdef MESSYTENDENCY
  CALL main_tendency_global_end
#endif
!!$  CALL main_channel_global_end ! op_pj_20110616
  ! op_bk_20130821
#ifndef __ICON__
  CALL main_qtimer_global_end
! op_bk_20130821
#endif

END SUBROUTINE messy_global_end

!==============================================================================

SUBROUTINE messy_write_output

  USE messy_main_switch       ! ONLY: USE_*
  USE messy_main_tracer_bi,   ONLY: main_tracer_write_output
  USE messy_main_channel_bi,  ONLY: messy_channel_write_output, IOMODE_OUT &
                                  , main_channel_write_output ! op_pj_20110616

! op_bk_20130820
#ifndef __ICON__
  USE messy_s4d_si,           ONLY:     s4d_write_output
  USE messy_scout_si,         ONLY:   scout_write_output
  USE messy_sorbit_si,        ONLY:  sorbit_write_output
  USE messy_timepos_e5,       ONLY: timepos_write_output
! op_bk_20130820
#endif

  ! op_bk_20130904+
!   CALL main_tracer_write_output(1)
  ! op_bk_20130904-

! op_bk_20130820
#ifndef __ICON__
  IF (USE_SCOUT)   CALL   scout_write_output
  IF (USE_S4D)     CALL     s4d_write_output    ! force output & new file
  IF (USE_SORBIT)  CALL  sorbit_write_output    ! force output
  IF (USE_TIMEPOS) CALL timepos_write_output
! op_bk_20130820
#endif

  CALL main_channel_write_output ! op_pj_20110616
  CALL messy_channel_write_output(IOMODE_OUT)

END SUBROUTINE messy_write_output

!==============================================================================

SUBROUTINE messy_free_memory

  ! free memory (called by free_memory in mo_memory_streams)

  USE messy_main_switch   ! ONLY: USE_*
  ! op_bk_20130821
#ifndef __ICON__
  USE messy_main_qtimer_bi,    ONLY: main_qtimer_free_memory
! op_bk_20130821
#endif
  ! op_bk_20130927+
!   USE messy_main_tracer_bi,    ONLY: main_tracer_free_memory
  ! op_bk_20130927-
  USE messy_main_channel_bi,   ONLY: main_channel_free_memory
  USE messy_main_data_bi,      ONLY: main_data_free_memory
! op_bk_20130820
#ifndef __ICON__
  USE messy_main_import_bi,    ONLY: main_import_free_memory
  USE messy_ncregrid_tools_e5, ONLY: messy_ncregrid_free_memory
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,  ONLY: main_tendency_free_memory
#endif  

  USE messy_aeropt_e5,      ONLY:     aeropt_free_memory
  USE messy_airsea_si,      ONLY:     airsea_free_memory
  USE messy_cloud_e5,       ONLY:      cloud_free_memory
  USE messy_convect_e5,     ONLY:    convect_free_memory
  USE messy_cvtrans_si,     ONLY:    cvtrans_free_memory
  USE messy_d14co_e5,       ONLY:      d14co_free_memory
  USE messy_ddep_si,        ONLY:       ddep_free_memory
  USE messy_e4chem_si,      ONLY:     e4chem_free_memory
  USE messy_gmxe_e5,        ONLY:       gmxe_free_memory
  USE messy_gwave_si,       ONLY:      gwave_free_memory
  USE messy_hetchem_e5,     ONLY:    hetchem_free_memory
  USE messy_jval_si,        ONLY:       jval_free_memory
  USE messy_lnox_si,        ONLY:       lnox_free_memory
  USE messy_m7_si,          ONLY:         m7_free_memory
  USE messy_made_si,        ONLY:       made_free_memory
  USE messy_mecca1_e5,      ONLY:     mecca1_free_memory
  USE messy_mecca_si,       ONLY:      mecca_free_memory
  USE messy_megan_si,       ONLY:      megan_free_memory
  USE messy_o3orig_si,      ONLY:     o3orig_free_memory
  USE messy_offemis_si,     ONLY:    offemis_free_memory
  USE messy_onemis_si,      ONLY:     onemis_free_memory
  USE messy_photo_e5,       ONLY:      photo_free_memory
  USE messy_plumegas_si,    ONLY:   plumegas_free_memory
  USE messy_qbo_si,         ONLY:        qbo_free_memory
  USE messy_rad4all_e5,     ONLY:    rad4all_free_memory
  USE messy_satsims_e5,     ONLY:    satsims_free_memory
  USE messy_scalc_si,       ONLY:      scalc_free_memory
  USE messy_scav_si,        ONLY:       scav_free_memory
  USE messy_scout_si,       ONLY:      scout_free_memory
  USE messy_sedi_si,        ONLY:       sedi_free_memory
  USE messy_s4d_si,         ONLY:        s4d_free_memory
  USE messy_sorbit_si,      ONLY:     sorbit_free_memory
  USE messy_spe_e5,         ONLY:        spe_free_memory
  USE messy_spacenox_e5,    ONLY:   spacenox_free_memory
  USE messy_timepos_e5,     ONLY:    timepos_free_memory
  USE messy_tnudge_si,      ONLY:     tnudge_free_memory
  USE messy_trexp_si,       ONLY:      trexp_free_memory
! op_bk_20130820
#endif

  IMPLICIT NONE

! op_bk_20130820
#ifndef __ICON__
  ! special requirements:
  ! - (none, all CALLs are in alphabetical order)
  IF (USE_AEROPT)   CALL     aeropt_free_memory
  IF (USE_AIRSEA)   CALL     airsea_free_memory
  IF (USE_CLOUD)    CALL      cloud_free_memory
  IF (USE_CONVECT)  CALL    convect_free_memory
  IF (USE_CVTRANS)  CALL    cvtrans_free_memory
  IF (USE_D14CO)    CALL      d14co_free_memory
  IF (USE_DDEP)     CALL       ddep_free_memory
  IF (USE_E4CHEM)   CALL     e4chem_free_memory
  IF (USE_GMXE)     CALL       gmxe_free_memory
  IF (USE_GWAVE)    CALL      gwave_free_memory
  IF (USE_HETCHEM)  CALL    hetchem_free_memory
  IF (USE_JVAL)     CALL       jval_free_memory
  IF (USE_LNOX)     CALL       lnox_free_memory
  IF (USE_M7)       CALL         m7_free_memory
  IF (USE_MADE)     CALL       made_free_memory
  IF (USE_MECCA1)   CALL     mecca1_free_memory
  IF (USE_MECCA)    CALL      mecca_free_memory
  IF (USE_MEGAN)    CALL      megan_free_memory
  IF (USE_O3ORIG)   CALL     o3orig_free_memory
  IF (USE_OFFEMIS)  CALL    offemis_free_memory
  IF (USE_ONEMIS)   CALL     onemis_free_memory
  IF (USE_PHOTO)    CALL      photo_free_memory
  IF (USE_PLUMEGAS) CALL   plumegas_free_memory
  IF (USE_QBO)      CALL        qbo_free_memory
  IF (USE_RAD4ALL)  CALL    rad4all_free_memory
  IF (USE_SATSIMS)  CALL    satsims_free_memory
  IF (USE_SCALC)    CALL      scalc_free_memory
  IF (USE_SCAV)     CALL       scav_free_memory
  IF (USE_SCOUT)    CALL      scout_free_memory
  IF (USE_SEDI)     CALL       sedi_free_memory
  IF (USE_S4D)      CALL        s4d_free_memory
  IF (USE_SORBIT)   CALL     sorbit_free_memory
  IF (USE_SPE)      CALL        spe_free_memory
  IF (USE_SPACENOX) CALL   spacenox_free_memory
  IF (USE_TIMEPOS)  CALL    timepos_free_memory
  IF (USE_TNUDGE)   CALL     tnudge_free_memory
  IF (USE_TREXP)    CALL      trexp_free_memory
! op_bk_20130820
#endif

! op_bk_20130820
#ifndef __ICON__
  CALL main_import_free_memory
! op_bk_20130820
#endif
  ! op_bk_20130927+
!   CALL main_tracer_free_memory
  ! op_bk_20130927-
  CALL main_channel_free_memory
  CALL main_data_free_memory
! op_bk_20130820
#ifndef __ICON__
  CALL messy_ncregrid_free_memory
! op_bk_20130820
#endif
  ! op_bk_20130821
#ifndef __ICON__
  CALL main_qtimer_free_memory
! op_bk_20130821
#endif
#ifdef MESSYTENDENCY
  CALL main_tendency_free_memory
#endif  


END SUBROUTINE messy_free_memory

!==============================================================================
