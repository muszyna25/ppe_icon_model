!*****************************************************************************
MODULE messy_main_data_bi
!*****************************************************************************
!*****************************************************************************
!                Time-stamp: <2013-09-04 15:46:12 kern_ba>
!*****************************************************************************
#if defined(ECHAM5)
! Authors: Rolf Sander,     MPICH,  2002-2004
!          Patrick Joeckel, MPICH,  2004
!
! The module messy_main_data_bi defines channels to transfer data
! between the ECHAM5 base model (physc.f90, vdiff.f90, etc.) and the MESSy
! submodels. Several quantities needed here are already
! defined in g2a or g3b. Submodels can and should still access them via
! messy_main_data_bi since they are linked here with
! new_channel_object_reference.

! WARNING: Do not use values in mo_memory_g1a or mo_memory_g1b (e.g. tf)
! because of problems with leap frog and time filter!

  ! NON-FIELD PARAMETER TRANSFER
  ! - TIME CONTROL AND TIME FILTER
#ifndef MESSYTIMER
  USE mo_time_control,  ONLY: time_step_len, delta_time   &  ! REAL
                            , lstart, lresume, lfirst_day &  ! LOGICAL
                            , L_TRIGGER_RESTART, lfirst_cycle ! LOGICAL
#else
  USE messy_main_timer, ONLY: time_step_len
#endif
  USE mo_semi_impl,     ONLY: eps                            ! REAL
  ! - GRID CONTROL
  USE mo_control,       ONLY: ngl, nlev, nlevp1, nlon     &  ! INTEGER
                            , nhgl                        &  ! INTEGER
                            , nsp                         &  ! INTEGER
                            , nmp1                        &  ! INTEGER
                            , nvclev, nn, nm, nk          &  ! INTEGER
                            , vct                         &  ! REAL(:)
#if defined(E5202A) || defined(E5301) 
                            , nrow                        &  ! INTEGER(:)
#endif
                            , lmidatm, lcolumn            &  ! LOGICAL
                            , lcouple                        ! LOGICAL

#if defined(E5302)
  USE mo_control,       ONLY: lipcc                          ! LOGICAL
#endif
  USE mo_gaussgrid,     ONLY: gl_gmu                      &  ! REAL(:)
                            , gl_gw, gridarea             &  ! REAL(:)
                            , philon, philat              &  ! REAL(:)
                            , coslon, sinlon              &  ! REAL(:)
                            , gl_twomu, gl_sqcst             ! REAL(:)
#if defined(E5202A) || defined(E5301) 
  USE mo_geoloc,        ONLY: gboxarea                    &  ! REAL(:)
                            , gboxarea_2d                 &  ! REAL(:,:)
                            , philat_2d, philon_2d        &  ! REAL(:,:)
                            , sqcst_2d                    &  ! REAL(:,:)
                            , ilat, ilon                  &  ! INTEGER(:,:)
                            , coriol_2d                   &
                            , twomu_2d                       ! REAL(:,:)
#endif
#if defined(E5302)
  USE mo_geoloc,        ONLY: gboxarea_2d                 &  ! REAL(:,:)
                            , philat_2d, philon_2d        &  ! REAL(:,:)
                            , sqcst_2d                    &  ! REAL(:,:)
                            , ilat, ilon                  &  ! INTEGER(:,:)
                            , coriol_2d                   &
                            , twomu_2d                       ! REAL(:,:)
#endif
  USE mo_hyb,           ONLY: apzero, apsurf              &  ! REAL
                            , NPLVP1, NLMSGL, nlevm1      &  ! INTEGER
                            , NPLVP2, NLMSLP              &  ! INTEGER
                            , ceta                        &  ! REAL(:)
                            , cetah                       &  ! REAL(:)
                            , nplev, delb, delpr

  ! - DECOMPOSITION -> nproma, npromz, ngpblks, nllev, nllevp1, lnsp, nlm, snsp
  USE messy_main_mpi_bi,    ONLY: p_pe, dcl

  ! -----------------------------------------------------------------
  ! - ONLY FOR USE WITHIN 'LOCAL LOOP'
  ! -----------------------------------------------------------------
#if defined(E5202A) || defined(E5301)
  ! ... at t
  USE mo_sc1,    ONLY: &
       alnpr    &
       , alpha  &
       , alps, alpste & ! SURFACE PRESSURE AND TENDENCY
       , d      &       ! DIVERGENCE
       , vo     &       ! VORTICITY
       , qte    &       ! SPECIFIC HUMIDITY TENDENCY
       , rh     &       ! SURFACE GEOPOTENTIAL
       , t, tte       & ! TEMPERATURE AND TENDENCY
       , u      &       ! U-WIND
       , dudl   &
       , v      &       ! V-WIND
       , dvdl   &
       , vervel &
       , vol    &
       , vom    &
       , xlte   &       ! LIQUID WATER TENDENCY
       , xite           ! ICE TENDENCY

  ! -----------------------------------------------------------------
  ! - FOR USE OUTSIDE 'LOCAL LOOP'
  ! -----------------------------------------------------------------
  ! ... at t
  USE mo_scan_buffer, ONLY:   & 
       alnpr_scb              &
       , alpha_scb            &
       , alps_scb, alpste_scb & ! SURFACE PRES. AND TEND.
       , d_scb                & ! DIVERGENCE  
       , vo_scb               & ! VORTICITY
       , qte_scb              & ! SPECIFIC HUMIDITY TENDENCY
       , rh_scb               & ! SURFACE GEOPOTENTIAL
       , t_scb, tte_scb       & ! TEMPERATURE AND TENDENCY
       , u_scb                & ! U-WIND
       , dudl_scb             &
       , v_scb                & ! V-WIND
       , dvdl_scb             &
       , vervel_scb           &
       , vol_scb              & 
       , vom_scb              &
       , xlte_scb             & ! LIQUID WATER TENDENCY
       , xite_scb             & ! ICE TENDENCY
       , dalpsl_scb, dalpsm_scb ! op_pj_20120315

#endif
#if defined(E5302)
  ! -----------------------------------------------------------------
  ! - FOR USE OUTSIDE 'LOCAL LOOP'
  ! -----------------------------------------------------------------
  ! ... at t
  USE mo_scan_buffer, ONLY:   & 
       alnpr_scb   => alnpr   &
       , alpha_scb => alpha   &
       , alps_scb  => alps    & ! SURFACE PRES.
       , alpste_scb=> alpste  & ! SURFACE PRES. TEND.
       , d_scb     => d       & ! DIVERGENCE  
       , vo_scb    => vo      & ! VORTICITY
       , qte_scb   => qte     & ! SPECIFIC HUMIDITY TENDENCY
       , rh_scb    => rh      & ! SURFACE GEOPOTENTIAL
       , t_scb     => t       & ! TEMPERATURE
       , tte_scb   => tte     & ! TEMPERATURE TENDENCY
       , u_scb     => u       & ! U-WIND
       , dudl_scb  => dudl    &
       , v_scb     => v       & ! V-WIND
       , dvdl_scb  => dvdl    &
       , vervel_scb=> vervel  &
       , vol_scb   => vol     & ! WIND TENDENCY
       , vom_scb   => vom     & ! WIND TENDENCY
       , xlte_scb  => xlte    & ! LIQUID WATER TENDENCY
       , xite_scb  => xite    & ! ICE TENDENCY
       , dalpsl_scb => dalpsl & ! op_pj_20120315
       , dalpsm_scb => dalpsm   ! op_pj_20120315
#endif

  USE mo_memory_gl,   ONLY:   &
       q                      & ! SPECIFIC HUMIDITY
       , xl                   & ! LIQUID WATER
       , xi                     ! ICE
  ! -----------------------------------------------------------------
  ! ... at t-1
  USE mo_memory_g1a,  ONLY:   &
       vom1                   & ! VORTICITY
       , dm1                  & ! DIVERGENCE
       , tm1                  &
       , alpsm1               & ! SURFACE PRESSURE
       , dalpslm1             &
       , dalpsmm1             &
       , qm1                  & ! SPECIFIC HUMIDITY
       , xlm1                 & ! LIQUID WATER
       , xim1                   ! ICE
  USE mo_memory_g2a,  ONLY:   &
       um1                    & ! U-WIND
       , vm1                  & ! V-WIND
       , dtlm1                &
       , dtmm1                &
       , dudlm1               &
       , dvdlm1
  ! -----------------------------------------------------------------
  ! .. at ???
  USE mo_memory_g3b, ONLY: geosp, slf, slm, seaice, forest, vgrat        &
  ! fb_mk_20100212+
  ! Variables added for sub-model mlocean.
                         , siced, fluxres, ahflw                         &
                         , trflw, trfli, soflw, sofli, amlcorr, amlcorac & 
                         , amlheatac, ahfres, ahfice, ahfcon             &
                         , qres, ahfli, friac                            &
  ! fb_mk_20100212-
                         , sni, tsi, alb                                 &
                         , vlt, albedo, aclc, u10, v10, tsoil, wsmx, ws  &
                         , az0, az0w, az0i, az0l, tslm1, aps, alake      &
                         , acdnc, aclcac, aclcov, relhum, rintop         &
                         , xvar, xskew, aprl, aprs, qvi, xlvi, topmax    &
                         , aprc, xtec, tsw, xivi                         &
                         , glac, tke &       ! mz_ht_20070629
                         , siced, ocu, ocv & ! mz_ap_20070913
                         , aprflux           ! mz_ab_20100829

  ! -----------------------------------------------------------------
  ! -- at ???
#if defined(E5202A) || defined(E5301)
  USE mo_tmp_buffer, ONLY: aphm1,aphp1,apm1,app1,loland,loglac,geom1 &
                         , rsfc, rsfl, ssfc, ssfl
#endif
  ! -----------------------------------------------------------------

  ! MESSy
  USE messy_main_constants_mem, ONLY: dp, DTR

  IMPLICIT NONE
  PUBLIC
  SAVE

  ! ECHAM5 PARAMETERS
  CHARACTER(LEN=*), PARAMETER :: modstr = 'ECHAM5'

#if defined(E5202A)
  CHARACTER(LEN=*), PARAMETER :: modver = '5.2.02a'
  CHARACTER(LEN=*), PARAMETER :: scb = '_scb'        ! op_pj_20120309
#endif
#if defined(E5301)
  CHARACTER(LEN=*), PARAMETER :: modver = '5.3.01'
  CHARACTER(LEN=*), PARAMETER :: scb = '_scb'        ! op_pj_20120309
#endif
#if defined(E5302)
  CHARACTER(LEN=*), PARAMETER :: modver = '5.3.02'
  CHARACTER(LEN=*), PARAMETER :: scb = ''            ! op_pj_20120309
  REAL(DP), POINTER :: vo(:,:)       => NULL()
  REAL(DP), POINTER :: d(:,:)        => NULL()
  REAL(DP), POINTER :: t(:,:)        => NULL()
  REAL(DP), POINTER :: alps(:)       => NULL()
  REAL(DP), POINTER :: u(:,:)        => NULL()
  REAL(DP), POINTER :: dudl(:,:)     => NULL()
  REAL(DP), POINTER :: v(:,:)        => NULL()
  REAL(DP), POINTER :: dvdl(:,:)     => NULL()
  REAL(DP), POINTER :: vol(:,:)      => NULL()
  REAL(DP), POINTER :: vom(:,:)      => NULL()
  REAL(DP), POINTER :: rh(:,:)       => NULL()
  REAL(DP), POINTER :: qte(:,:)      => NULL()
  REAL(DP), POINTER :: xlte(:,:)     => NULL()
  REAL(DP), POINTER :: xite(:,:)     => NULL()
  REAL(DP), POINTER :: xtte(:,:,:)   => NULL()
  REAL(DP), POINTER :: xtte_a(:,:,:) => NULL()
  REAL(DP), POINTER :: tte(:,:)      => NULL()
  REAL(DP), POINTER :: alpste(:)     => NULL()
  REAL(DP), POINTER :: alnpr(:,:)    => NULL()
  REAL(DP), POINTER :: alpha(:,:)    => NULL()
  REAL(DP), POINTER :: vervel(:,:)   => NULL()
  INTEGER, DIMENSION(3) :: nrow = (/0,0,0/)
  ! mo_geoloc
  REAL(DP), POINTER :: gboxarea(:)   => NULL()
  ! mo_tmp_buffer
  REAL(dp),POINTER, DIMENSION(:,:) :: aphm1  => NULL()
  REAL(dp),POINTER, DIMENSION(:,:) :: aphp1  => NULL()
  REAL(dp),POINTER, DIMENSION(:,:) :: apm1   => NULL()
  REAL(dp),POINTER, DIMENSION(:,:) :: app1   => NULL()
  LOGICAL, POINTER, DIMENSION(:)   :: loland => NULL()
  LOGICAL, POINTER, DIMENSION(:)   :: loglac => NULL()
  REAL(dp),POINTER, DIMENSION(:,:) :: geom1  => NULL()
  REAL(dp),POINTER, DIMENSION(:)   :: rsfc   => NULL()
  REAL(dp),POINTER, DIMENSION(:)   :: ssfc   => NULL()
  REAL(dp),POINTER, DIMENSION(:)   :: rsfl   => NULL()
  REAL(dp),POINTER, DIMENSION(:)   :: ssfl   => NULL()
#endif

  LOGICAL, PARAMETER :: l2tls = .FALSE.      ! mz_pj_20090414

  ! LOCALIZED PARAMETERS
  ! - GRID CONTROL
  INTEGER :: nllev
  INTEGER :: nllevp1
  INTEGER :: lnsp             ! LS:
  INTEGER :: snsp             ! SP:
  INTEGER :: nlm              ! FAS:
  INTEGER :: jrow             ! GP: CURRENT ROW
  INTEGER :: jglat            ! GP: global continuous latitude index
  INTEGER :: kproma           ! GP: VECTOR LENGTH OF CURRENT ROW
  INTEGER :: nproma           ! GP: VECTOR LENGTH
  INTEGER :: npromz           ! GP: VECTOR LENGTH OF LAST ROW
  INTEGER :: ngpblks          ! GP: NUMBER OF ROWS
  INTEGER :: nglon, nglat     ! GP: NUMBER OF LON AND LAT (1 if lcolumn=T)
#ifndef MESSYTIMER
  INTEGER :: YEAR,MONTH,DAY,HOUR,MINUTE,SECOND
  INTEGER :: YEAR_START,MONTH_START,DAY_START &
           , HOUR_START,MINUTE_START,SECOND_START
  INTEGER :: YEAR_NEXT,MONTH_NEXT,DAY_NEXT &
           , HOUR_NEXT,MINUTE_NEXT,SECOND_NEXT
  INTEGER :: DAYOFYEAR        ! day of year [day] ! mz_ab_20080309
  INTEGER :: current_time_step
#endif

  ! LOGICAL (needed for MMD COULPLING ECHAM5<->COSMO<->COSMO)
  LOGICAL :: L_IS_CLIENT = .FALSE.

  ! SCALAR
  ! NOTE: THESE ARE ONLY USED FOR THE CONNECTION BETWEEN RAD4ALL AND
  !       THE ORIGINAL ECHAM5-RADIATION ROUTINES!
  !       DO NOT USE THEM FOR MESSy-SUBMODELS. ACCESS INFORMATION
  !       FROM SUBMODEL ORBIT VIA CHANNEL-OBJECTS, INSTEAD.
  REAL(dp) :: cdisse   ! distance Sun - Earth in AU
  REAL(dp) :: cdissem  ! distance Sun - Earth in AU (at radiation time step)

  ! um_ak_20110627+
  ! 1-dimensional (vertical) fields
  REAL(dp), DIMENSION(:), POINTER :: hyam => NULL()
  REAL(dp), DIMENSION(:), POINTER :: hybm => NULL()
  ! um_ak_20110627-

  ! 2-dimensional fields
  REAL(dp), DIMENSION(:,:), POINTER :: &
    qflux       => NULL(), &
    s_heatflux  => NULL(), &
    ! mz_ht_20070421+    
    l_heatflux  => NULL(), &
    ! mz_ht_20070421-
    cdni        => NULL(), &
    cdnl        => NULL(), &
    cdnw        => NULL(), &
    cfmi        => NULL(), &
    cfml        => NULL(), &
    cfmw        => NULL(), &
    cfnci       => NULL(), &
    cfncl       => NULL(), &
    cfncw       => NULL(), &
    cossza_2d   => NULL(), &
    cvs         => NULL(), &
    cvsc        => NULL(), &
    cvw         => NULL(), &
    fws         => NULL(), &
    icecov      => NULL(), &
    seacov      => NULL(), &
    prc         => NULL(), &
    prl         => NULL(), &
    rco_leaf    => NULL(), &
    rh_2m       => NULL(), &
    rii         => NULL(), &
    ril         => NULL(), &
    riw         => NULL(), &
    srfl        => NULL(), &
    tsurf_2d    => NULL(), &
    tvi         => NULL(), &
    tvir        => NULL(), &
    tvl         => NULL(), &
    tvw         => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: &
    wind10_2d   => NULL(), &
    ! mz_ap_20070913+
    !for ocean model
    wind10w_2d  => NULL(), &
    awhea_2d    => NULL(), &
    aicon_2d   => NULL(), &
    aiqre_2d    => NULL(), &
    awfre_2d    => NULL(), &
    aifre_2d    => NULL(), &
    aiust_2d    => NULL(), &
    aivst_2d    => NULL(), &
    awust_2d    => NULL(), &
    awvst_2d    => NULL(), &
    ! for HD model
    aros_2d     => NULL(), &
    adrain_2d   => NULL(), &
    apmecal_2d  => NULL(), &
    disch_2d    => NULL(), &
    apmebco_2d  => NULL(), &  ! mz_ap_20090519
    ! mz_ap_20070913-
    zcdh_2d     => NULL(), &
    zlatkf_2d   => NULL(), &
    zsenkf_2d   => NULL(), &
    zust_2d     => NULL(), &
    tslnew      => NULL(), &
    decomp_gp_jp => NULL(), &
    decomp_gp_jr => NULL(), &
    decomp_gp_pe => NULL(), &
    ! mz_pj_20071130+
    ! op_pj_20101212+
    coslon_2d    => NULL(), &
    sinlon_2d    => NULL(), &
    ! op_pj_20101212-
    ! op_pj_20110714+
    coslat_2d    => NULL(), &
    sinlat_2d    => NULL()
    ! op_pj_20110714-


  REAL(dp), DIMENSION(:,:), POINTER :: &
    pmsl        => NULL(), & ! ub_ch_20100721
    zi0_2d      => NULL(), &
    srfl_2d     => NULL(), &
    ahfsw       => NULL(), &  ! fb_mk_20100514+
    ahfsi       => NULL(), &
    cvsi        => NULL(), &
    evapi       => NULL(), &
    fsnet       => NULL()     ! fb_mk_20100514-
  ! mz_pj_20071130-

  ! 3-dimensional fields
  REAL(dp), DIMENSION(:,:,:), POINTER :: &
    etadot_3d  => NULL(), &
    geopot_3d  => NULL(), &
    geopoti_3d => NULL(), &
    grmass     => NULL(), &
    grmassdry  => NULL(), &   ! op_pj_20100713
    grvol      => NULL(), &
    ilab       => NULL(), &
    press_3d   => NULL(), &
    pressi_3d  => NULL(), &
    qm1_3d     => NULL(), &
    qte_3d     => NULL(), &
    qtec       => NULL(), &
    relo3_3d   => NULL(), &
    rhum_3d    => NULL(), &
    rinum_3d   => NULL(), &
    tm1_3d     => NULL(), &
    tpot_3d    => NULL(), &
    tte_3d     => NULL(), &
    tvirt_3d   => NULL(), &
    xim1_3d    => NULL(), &
    xite_3d    => NULL(), &
    xlm1_3d    => NULL(), &
    xlte_3d    => NULL(), &
    flxs_3d    => NULL(), &     ! ut_kt_20041007
    flxt_3d    => NULL(), &     ! ut_kt_20041007
    vmixtau    => NULL(), &
    vdiffp     => NULL(), &
    gwdrag_u   => NULL(), &     ! mz_pj_20071126
    gwdrag_v   => NULL(), &     ! mz_pj_20071126
    gwflux_u   => NULL(), &     ! mz_pj_20071126
    gwflux_v   => NULL(), &     ! mz_pj_20071126
    xtecl      => NULL(), &     ! mz_ht_20071221
    xteci      => NULL() !!$, &     ! mz_ht_20071221
!!$    grmass0    => NULL()        ! um_ak_20100319

  ! 4-dimensional fields
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: &
    pxtems     => NULL()

  ! SUBROUTINES
  !PUBLIC :: main_data_initialize
  !PUBLIC :: main_data_init_memory
  !PUBLIC :: main_data_local_start
  !PUBLIC :: main_data_global_start
  !PUBLIC :: main_data_free_memory
  !
  PRIVATE :: channel_halt

CONTAINS

  !***************************************************************************
  SUBROUTINE main_data_initialize

#ifndef MESSYTIMER
    USE messy_main_bmluse_bi,    ONLY: get_date_components, start_date
#else
    USE messy_main_timer,        ONLY: timer_set_time_step_len     
#endif

    IMPLICIT NONE

    nllev   = dcl%nllev
    nllevp1 = dcl%nllevp1
    nlm     = dcl%nlm
    lnsp    = dcl%lnsp
    snsp    = dcl%snsp
    nproma  = dcl%nproma
    npromz  = dcl%npromz
    ngpblks = dcl%ngpblks
    nglon   = dcl%nglon
    nglat   = dcl%nglat

#ifndef MESSYTIMER
    ! SET START TIME
    CALL get_date_components(start_date     &
         ,YEAR_START,MONTH_START,DAY_START  &
         ,HOUR_START,MINUTE_START,SECOND_START)  
#endif    
    ! um_ak_20100419+
#ifdef MESSYTIMER
    ! set integration time step length
    CALL timer_set_time_step_len(l2tls)
#endif
    ! um_ak_20100419-

  END SUBROUTINE main_data_initialize
  !***************************************************************************

  !***************************************************************************
  SUBROUTINE main_data_init_memory

    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi,  ONLY: ntrac_gp
    USE messy_main_mpi_bi,         ONLY: dcl

    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_channel_object_reference    &
                                      , new_attribute, get_channel_object
    USE messy_main_channel_repr,  ONLY: get_representation_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_init_memory'
    INTEGER :: status
    INTEGER :: reprid    
    INTEGER :: jp, jr

    ! create new channel
    CALL new_channel (status, modstr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    ! um_ak_20110627+
    ! ########################################################################
    ! ---------------------------- GP_1D_LEV ---------------------------------
    ! ########################################################################
    ! ------------------------------------------------------------------------
    CALL get_representation_id(status, 'GP_1D_LEV', reprid)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'dhyam', &
         p1=hyam, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhyam', 'long_name' &
         , c='hybrid A coefficient at layer midpoints')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhyam', 'units', c='Pa')
    CALL channel_halt(substr, status)
    hyam(:) = (vct(1:nvclev-1) + vct(2:nvclev))/2
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'dhybm', &
         p1=hybm, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhybm', 'long_name' &
            , c='hybrid B coefficient at layer midpoints')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhybm', 'units', c='Pa')
    CALL channel_halt(substr, status)
    hybm(:) = (vct(nvclev+1:2*nvclev-1) + vct(nvclev+2:2*nvclev))/2
    ! ------------------------------------------------------------------------
    ! um_ak_20110627-

    ! ########################################################################
    ! ---------------------------- GP_2D_HORIZONTAL --------------------------
    ! ########################################################################
    CALL get_representation_id(status, 'GP_2D_HORIZONTAL', reprid)
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! --- CHANNEL OBJECT REFERENCES ------------------------------------------
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'u10', modstr, 'u10')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u10', 'long_name', c='10m u-velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u10', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'v10', modstr, 'v10')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v10', 'long_name', c='10m v-velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v10', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'tsoil', modstr, 'tsoil')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsoil', &
         'long_name', c='deep soil temperatures')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsoil', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'wsmx', modstr, 'wsmx')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wsmx', &
         'long_name', c='field capacity of soil')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wsmx', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'ws', modstr, 'ws')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ws', 'long_name', c='soil wetness')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ws', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'geosp', modstr, 'geosp')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geosp', &
         'long_name', c='surface geopotential')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geosp', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'slf', modstr, 'slf')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slf', 'long_name', c='land fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slf', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'alb', modstr, 'alb')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alb', 'long_name', &
         c='surface background albedo')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alb', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'tsi', modstr, 'tsi')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsi', 'long_name', &
         c='surface temperature of ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsi', 'units', c='K')
    CALL channel_halt(substr, status)
    ! fb_mk_20100215+
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'tsw', modstr, 'tsw')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsw', 'long_name', &
         c='surface temperature of water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsw', 'units', c='K')
    CALL channel_halt(substr, status)
    ! fb_mk_20100215-
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'sni', modstr, 'sni')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sni', 'long_name', &
         c='water equivalent of snow on ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sni', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! 0 = sea, 1 = land
    CALL new_channel_object_reference(status, 'g3b', 'slm', modstr, 'slm')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slm', 'long_name', c='land mask')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slm', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'glac',modstr, 'glac')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'glac', 'long_name', &
         c='fraction of land covered by glaciers')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'glac', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! 0 = lake, 1 = land
    CALL new_channel_object_reference(status, 'g3b', 'alake', modstr, 'alake')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alake', &
         'long_name', c='lake fraction of grid box')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alake', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! 0 = liquid ocean, 1 = ice
    ! This value is set to 0 over land, although it is actually undefined there
    CALL new_channel_object_reference(status, &
         'g3b', 'seaice', modstr, 'seaice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'seaice', &
         'long_name', c='seaice fraction rel to ocean')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'seaice', 'units', c='-')
    CALL channel_halt(substr, status)
    ! fb_mk_20100215+
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, &
         'g3b', 'siced', modstr, 'siced')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'siced', &
         'long_name', c='ice depth')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'siced', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, &
         'g3b', 'fluxres', modstr, 'fluxres')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fluxres', &
         'long_name', c='')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fluxres', 'units', c='')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object_reference(status, &
         'g3b', 'ahfli', modstr, 'ahfli')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfli', &
         'long_name', c='latent heat flux over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfli', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object_reference(status, &
         'g3b', 'ahflw', modstr, 'ahflw')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahflw', &
         'long_name', c='latent heat flux over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahflw', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object_reference(status, &
         'g3b', 'trflw', modstr, 'trflw')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'trflw', &
         'long_name', c='long wave flux over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'trflw', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object_reference(status, &
         'g3b', 'trfli', modstr, 'trfli')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'trfli', &
         'long_name', c='long wave flux over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'trfli', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object_reference(status, &
         'g3b', 'soflw', modstr, 'soflw')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'soflw', &
         'long_name', c='short wave flux over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'soflw', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object_reference(status, &
         'g3b', 'sofli', modstr, 'sofli')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sofli', &
         'long_name', c='short wave flux over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sofli', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, &
         'g3b', 'ahfice', modstr, 'ahfice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfice', &
         'long_name', c='conductive heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfice', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, &
         'g3b', 'amlcorr', modstr, 'amlcorr')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'amlcorr', &
         'long_name', c='mixed layer ocean flux correction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'amlcorr', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, &
         'g3b', 'qres', modstr, 'qres')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qres', &
         'long_name', c='residual heat flux for melting sea ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qres', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! fb_mk_20100215-
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, &
         'g3b', 'forest', modstr, 'forest')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'forest', &
         'long_name', c='forest fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'forest', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'vgrat', modstr, 'vgrat')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vgrat', &
         'long_name', c='vegetation fraction rel to land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vgrat', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'vlt', modstr, 'vlt')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vlt', &
         'long_name', c='leaf area index (LAI)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vlt', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, &
         'g3b', 'albedo', modstr, 'albedo')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'albedo', &
         'long_name', c='surface albedo')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'albedo', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'az0', modstr, 'az0')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0', &
         'long_name', c='roughness length')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'az0w', modstr, 'az0w')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0w', &
         'long_name', c='roughness length over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0w', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'az0i', modstr, 'az0i')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0i', &
         'long_name', c='roughness length over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0i', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'az0l', modstr, 'az0l')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0l', &
         'long_name', c='roughness length over land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0l', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'tslm1', modstr, 'tslm1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tslm1', &
         'long_name', c='surface temperature of land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tslm1', 'units', c='K')
    CALL channel_halt(substr, status)
!!$    ! ------------------------------------------------------------------------
!!$    CALL new_channel_object_reference(status, 'g3b', 'sn', modstr, 'sn')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr, 'sn', 'long_name', c='snow depth')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr, 'sn', 'units', c='m')
!!$    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    !
    ! ------------------------------------------------------------------------
    ! --- NEW CHANNEL OBJECTS ------------------------------------------------
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr, 'cdnl', p2=cdnl, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnl', &
         'long_name', c='neutral drag coeff., land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnl', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'cdnw', p2=cdnw, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnw', &
         'long_name', c='neutral drag coeff., water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnw', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'cdni', p2=cdni, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdni', &
         'long_name', c='neutral drag coeff., ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdni', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'cfml', p2=cfml, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfml', &
         'long_name', c='momentum drag coeff., land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfml', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'cfmw', p2=cfmw, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmw', &
         'long_name', c='momentum drag coeff., water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmw', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'cfmi', p2=cfmi, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmi', &
         'long_name', c='momentum drag coeff., ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmi', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    ! mz_lg_20041119+
    ! pressure thickness terms over land used to calculate momentum and
    ! heat exchange. See page 48, equations 3.3.3.2.b and 3.3.3.3.a in
    ! DKRZ Report No. 6: The ECHAM3 Atmospheric General Circulation Model,
    ! http://www.mpimet.mpg.de/en/extra/models/echam/echam3_DKRZ-ReportNo.6.pdf
    ! mz_lg_20041119-
    CALL new_channel_object(status, modstr,  'cfncl', p2=cfncl, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncl', &
         'long_name', c='exchange parameter, land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncl', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    ! mz_lg_20041119+
    ! over water
    ! mz_lg_20041119-
    CALL new_channel_object(status, modstr,  'cfncw', p2=cfncw, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncw', &
         'long_name', c='exchange parameter, water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncw', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    ! mz_lg_20041119+
    ! over ice
    ! mz_lg_20041119-
    CALL new_channel_object(status, modstr,  'cfnci', p2=cfnci, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfnci', &
         'long_name', c='exchange parameter, ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfnci', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'ril', p2=ril, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ril', &
         'long_name', c='Richardson number (land)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ril', 'units', c='1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'riw', p2=riw, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'riw', &
         'long_name', c='Richardson number (water)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'riw', 'units', c='1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'rii', p2=rii, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rii', &
         'long_name', c='Richardson number (ice)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rii', 'units', c='1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible, rank is lower than local
    CALL new_channel_object(status, modstr,  'tvir', p2=tvir, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvir', &
         'long_name', c='surface virtual temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvir', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'tvl', p2=tvl, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvl', &
         'long_name', c='surface virtual temperature (land)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvl', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'tvw', p2=tvw, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvw', &
         'long_name', c='surface virtual temperature (water)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvw', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'tvi', p2=tvi, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvi', &
         'long_name', c='surface virtual temperature (ice)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvi', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'srfl', p2=srfl, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'srfl', &
         'long_name', c='net surface radiative flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'srfl', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'fws', p2=fws, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fws', &
         'long_name', c='soil moisture stress function')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fws', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'rco_leaf', &
         p2=rco_leaf, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rco_leaf', &
         'long_name', c='leaf stomatal resistance')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rco_leaf', 'units', c='s m-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'rh_2m' &
         , p2=rh_2m, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rh_2m', &
         'long_name', c='relative humidity at 2m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rh_2m', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'zust' &
         , p2=zust_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zust', &
         'long_name', c='surface friction velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zust', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr, 'heat', &
         p2=zsenkf_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'heat', &
         'long_name', c='surface kinematic heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'heat', 'units', c='K m/s')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr, 'qflx', &
         p2=zlatkf_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qflx', &
         'long_name', c='surface kinematic moisture flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qflx', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'cdh', p2=zcdh_2d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdh', &
         'long_name', c='neutral drag for heat exchange ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdh', 'units', c='W m-2 s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'wind10', &
         p2=wind10_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10', &
         'long_name', c='10 m wind speed')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! mz_ap_20070913+
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'wind10w', &
         p2=wind10w_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10w', &
         'long_name', c='10 m wind speed over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10w', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'aros', &
         p2=aros_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aros', &
         'long_name', c='atmospheric runoff')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aros', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'adrain', &
         p2=adrain_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'adrain', &
         'long_name', c='atmospheric drainage')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'adrain', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'apmecal', &
         p2=apmecal_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'apmecal', &
         'long_name', c='(p - e) at glacier points')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'apmecal', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'disch', &
         p2=disch_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'disch', &
         'long_name', c='inflow on ocean grid (without ice melting [calvin model])')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'disch', 'units', c='m/s')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
! mz_ap_20090519+
    CALL new_channel_object(status, modstr,  'apmebco', &
         p2=apmebco_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'apmebco', &
         'long_name', c='vert.integr.tendencies of water for correction in hd')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'apmebco', 'units', c='m/s')
    CALL channel_halt(substr, status)
! mz_ap_20090519-
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'awhea', &
         p2=awhea_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awhea', &
         'long_name', c='net heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awhea', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'aicon', &
         p2=aicon_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aicon', &
         'long_name', c='surface downward heat flux on ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aicon', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'aiqre', &
         p2=aiqre_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aiqre', &
         'long_name', c='residual heat flux (sea-ice topmelt heat flux)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aiqre', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'awfre', &
         p2=awfre_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awfre', &
         'long_name', c='water flux into ocean')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awfre', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'aifre', &
         p2=aifre_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aifre', &
         'long_name', c='downward snow flux where sea-ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aifre', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'awust', &
         p2=awust_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awust', &
         'long_name', c='surface eastward stress where water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awust', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'awvst', &
         p2=awvst_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awvst', &
         'long_name', c='surface northward stress where water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'awvst', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'aiust', &
         p2=aiust_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aiust', &
         'long_name', c='surface eastward stress where sea ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aiust', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'aivst', &
         p2=aivst_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aivst', &
         'long_name', c='surface northward stress where sea ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aivst', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! mz_ap_20070913-
    ! ------------------------------------------------------------------------
    ! (from physc.f90); no local variable, pointer is used already
    CALL new_channel_object(status, modstr,  'qflux', &
         p2=qflux, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qflux', &
         'long_name', c='moisture flux at the surface')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qflux', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); no local variable, pointer is used already
    CALL new_channel_object(status, modstr,  's_heatflux', &
         p2=s_heatflux, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 's_heatflux', &
         'long_name', c='sensible heat flux at the surface')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 's_heatflux', 'units', c='W m-2')
    CALL channel_halt(substr, status)
! mz_ht_20070421+
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); no local variable, pointer is used already
    CALL new_channel_object(status, modstr,  'l_heatflux', &
         p2=l_heatflux, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'l_heatflux', &
         'long_name', c='latent heat flux at the surface')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'l_heatflux', 'units', c='W m-2')
    CALL channel_halt(substr, status)
! mz_ht_20070421-
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'tsurf', &
         p2=tsurf_2d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsurf', &
         'long_name', c='surface temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsurf', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'icecov', &
         p2=icecov, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'icecov', &
         'long_name', c='ice cover fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'icecov', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'seacov', &
         p2=seacov, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'seacov', &
         'long_name', c='sea cover fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'seacov', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); not convertible
    CALL new_channel_object(status, modstr,  'prc',  &
         p2=prc, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prc', &
         'long_name', c='convective precipitation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prc', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); not convertible
    CALL new_channel_object(status, modstr,  'prl', &
         p2=prl, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prl', &
         'long_name', c='large-scale precipitation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prl', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'cvs', p2=cvs, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvs', 'long_name', c='snow cover')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvs', 'units', c='fraction')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'cvsc', p2=cvsc, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvsc', 'long_name', &
         c='snow covered canopy')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvsc', 'units', c='fraction')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'cvw', p2=cvw, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvw', &
         'long_name', c='wet skin fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvw', 'units', c='fraction')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'tslnew', &
         p2=tslnew, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tslnew', &
         'long_name', c='land surface temperature for sensible heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tslnew', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! NOTE: THIS IS ONLY USED FOR THE CONNECTION BETWEEN RAD4ALL AND
    !       THE ORIGINAL ECHAM5-RADIATION ROUTINES!
    !       DO NOT USE IT FOR MESSy-SUBMODELS. ACCESS INFORMATION
    !       FROM SUBMODEL ORBIT VIA CHANNEL-OBJECT, INSTEAD.
    ! (from prerad.f90); no local variable, pointer already used
    CALL new_channel_object(status, modstr,  'cossza', &
         p2=cossza_2d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cossza', &
         'long_name', c='cos(solar zenith angle)')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! mz_pj_20071204+
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'zi0', &
         p2=zi0_2d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zi0', &
         'long_name', c='solar incidence')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'srfl_physc', &
         p2=srfl_2d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'srfl_physc', &
         'long_name', c='solar surface flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'srfl_physc', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! mz_pj_20071204-
    ! fb_mk_20100215+
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'ahfsi', &
         p2=ahfsi, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsi', &
         'long_name', c='sensible heat flux over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsi', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'ahfsw', &
         p2=ahfsw, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsw', &
         'long_name', c='sensible heat flux over water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ahfsw', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'evapi', &
         p2=evapi, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapi', &
         'long_name', c='evaporation over ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'evapi', 'units', c='')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'cvsi', &
         p2=cvsi, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvsi', &
         'long_name', c='snow cover ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvsi', 'units', c='fraction')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'fsnet', &
         p2=fsnet, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fsnet', &
         'long_name', c='net surface flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fsnet', 'units', c='W m-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! fb_mk_20100215-
    ! ------------------------------------------------------------------------
    ! (diagnostic; decomposition)
    CALL new_channel_object(status, modstr,  'decomp_gp_jp', &
         p2=decomp_gp_jp, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'decomp_gp_jp', &
         'long_name', c='vector index (1...nproma)')
    CALL channel_halt(substr, status)
    DO jp=1, nproma
       decomp_gp_jp(jp,:) = REAL(jp,DP)
    END DO
    ! ------------------------------------------------------------------------
    ! (diagnostic; decomposition)
    CALL new_channel_object(status, modstr,  'decomp_gp_jr', &
         p2=decomp_gp_jr, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'decomp_gp_jr', &
         'long_name', c='vector index (1...ngpblks)')
    CALL channel_halt(substr, status)
    DO jr=1, ngpblks
       decomp_gp_jr(:,jr) = REAL(jr,DP)
    END DO
    ! ------------------------------------------------------------------------
    ! (diagnostic; decomposition)
    CALL new_channel_object(status, modstr,  'decomp_gp_pe', &
         p2=decomp_gp_pe, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'decomp_gp_pe', &
         'long_name', c='processor number')
    CALL channel_halt(substr, status)
    decomp_gp_pe(:,:) = REAL(p_pe,DP)
    ! ------------------------------------------------------------------------
    ! ub_ch_20100721+
    ! pressure at sea level 
    CALL new_channel_object(status, modstr,  'pmsl', &
         p2=pmsl, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pmsl' &
         , 'long_name', c='sea-level pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pmsl', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ub_ch_20100721-
    ! ------------------------------------------------------------------------
    ! op_pj_20101212+
    ! ------------------------------------------------------------------------
    ! COS(longitude)
    CALL new_channel_object(status, modstr,  'coslon', &
         p2=coslon_2d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'coslon', &
         'long_name', c='cos(longitude)')
    CALL channel_halt(substr, status)
    coslon_2d(:,:) = COS(philon_2d(:,:)*DTR)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! SIN(longitude)
    CALL new_channel_object(status, modstr,  'sinlon', &
         p2=sinlon_2d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sinlon', &
         'long_name', c='sin(longitude)')
    CALL channel_halt(substr, status)
    sinlon_2d(:,:) = SIN(philon_2d(:,:)*DTR)
    ! ------------------------------------------------------------------------
    ! op_pj_20101212-
    ! op_pj_20110714+
    ! ------------------------------------------------------------------------
    ! COS(latitude)
    CALL new_channel_object(status, modstr,  'coslat', &
         p2=coslat_2d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'coslat', &
         'long_name', c='cos(latitude)')
    CALL channel_halt(substr, status)
    coslat_2d(:,:) = COS(philat_2d(:,:)*DTR)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! SIN(latitude)
    CALL new_channel_object(status, modstr,  'sinlat', &
         p2=sinlat_2d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sinlat', &
         'long_name', c='sin(latitude)')
    CALL channel_halt(substr, status)
    sinlat_2d(:,:) = SIN(philat_2d(:,:)*DTR)
    ! ------------------------------------------------------------------------
    ! op_pj_20110714-

    ! ########################################################################
    ! ----------------------------- GP_3D_MID --------------------------------
    ! ########################################################################
    CALL get_representation_id(status, 'GP_3D_MID', reprid)
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! --- CHANNEL OBJECT REFERENCES ------------------------------------------
    ! ------------------------------------------------------------------------
    ! qqq check if this is the correct wind data
    ! wind speed u-component
    CALL new_channel_object_reference(status, 'g2a', 'um1', modstr, 'um1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'um1', 'long_name', c= 'um1 * cos(lat)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'um1', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! mz_ab_20100906+
    CALL new_attribute(status, modstr, 'um1', 'positive', c='eastward')
    CALL channel_halt(substr, status)
    ! mz_ab_20100906-
    ! ------------------------------------------------------------------------
    ! qqq check if this is the correct wind data
    ! wind speed v-component
    CALL new_channel_object_reference(status, 'g2a', 'vm1', modstr, 'vm1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vm1', 'long_name', c= 'vm1 * cos(lat)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vm1', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! mz_ab_20100906+
    CALL new_attribute(status, modstr, 'vm1', 'positive', c='northward')
    CALL channel_halt(substr, status)
    ! mz_ab_20100906-
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'aclc', modstr, 'aclc')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aclc', &
         'long_name', c='large scale cloud cover')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aclc', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'g3b', 'acdnc', modstr, 'acdnc')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    !
    ! ------------------------------------------------------------------------
    ! --- NEW CHANNEL OBJECTS ------------------------------------------------
    ! ------------------------------------------------------------------------
    ! (from main_data_global_start)
    !                   pointer
    ! pressure at middle of box ("full level pressure")
    CALL new_channel_object(status, modstr,  'press', &
         p3=press_3d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press', 'long_name', c='pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'tpot', &
         p3=tpot_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tpot', &
         'long_name', c='potential temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tpot', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); not convertible
    CALL new_channel_object(status, modstr,  'rinum', &
         p3=rinum_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rinum', &
         'long_name', c='bulk Richardson number')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rinum', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'vdiffp', &
         p3=vdiffp, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vdiffp', &
         'long_name', c='rate of change of q due to vdiff scheme')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vdiffp', 'units', c='1/s')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'vmixtau', &
         p3=vmixtau, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vmixtau', &
         'long_name', c='inverse mixing timescale for vertical turbulence')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vmixtau', 'units', c='1/s')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); pointer already used
    CALL new_channel_object(status, modstr,  'ilab', p3=ilab, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ilab', &
         'long_name', c='index: below or in convective cloud')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ilab', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); pointer already used
    CALL new_channel_object(status, modstr,  'qtec', p3=qtec, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qtec', &
         'long_name', c='convective detrained humidity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qtec', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); pointer already used
    CALL new_channel_object(status, modstr,  'xtecl', p3=xtecl, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xtecl', &
         'long_name', c='convective detrained liquid')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xtecl', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); pointer already used
    CALL new_channel_object(status, modstr,  'xteci', p3=xteci, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xteci', &
         'long_name', c='convective detrained ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xteci', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
! op_pj_20120318+
!!$    ! (from physc.f90)
!!$    ! qqq should be reference to tm1 in g1a; unusable anyway
!!$    CALL new_channel_object(status, modstr,  'tm1', &
!!$         p3=tm1_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL new_channel_object_reference(status, 'g1a', 'tm1', modstr, 'tm1')
! op_pj_20120318-
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1', &
         'long_name', c='dry air temperature (tm1)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1', 'units', c='K')
    CALL channel_halt(substr, status)
! op_pj_20120318+
    CALL get_channel_object(status, modstr, 'tm1', p3=tm1_3d)
    CALL channel_halt(substr, status)
! op_pj_20120318-
    ! ------------------------------------------------------------------------
! op_pj_20120318+
!!$    ! (from physc.f90)
!!$    ! qqq mo_sc1.f90; this is not a channel object ; do not use from here !!!
!!$    !     -> no feedback possible; since only copy
!!$    CALL new_channel_object(status, modstr,  'tte', &
!!$         p3=tte_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL new_channel_object_reference(status, 'scnbuf', 'tte'//scb &
         , modstr, 'tte')
! op_pj_20120318-
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tte', &
         'long_name', c='dry air temperature tendency (tte)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tte', 'units', c='K s-1')
    CALL channel_halt(substr, status)
! op_pj_20120318+
    CALL get_channel_object(status, modstr, 'tte', p3=tte_3d)
    CALL channel_halt(substr, status)
! op_pj_20120318-
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable (from mo_tmp_buffer.f90)  converted to
    !                   pointer
    CALL new_channel_object(status, modstr,  'geopot', &
         p3=geopot_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopot', &
         'long_name', c='geopotential')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopot', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
! op_pj_20120309+
!!$    ! (from physc.f90)
!!$    ! qqq should be reference to qm1 in g1a; unusable anyway
!!$    CALL new_channel_object(status, modstr,  'qm1', p3=qm1_3d, reprid=reprid)
    CALL new_channel_object_reference(status, 'g1a', 'qm1', modstr, 'qm1')
! op_pj_20120309+
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qm1', &
         'long_name', c='specific humidity (qm1)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qm1', 'units', c='kg kg-1')
    CALL channel_halt(substr, status)
! op_pj_20120309+
    CALL get_channel_object(status, modstr, 'qm1', p3=qm1_3d)
    CALL channel_halt(substr, status)
! op_pj_20120309-
    ! ------------------------------------------------------------------------
! op_pj_20120309+
!!$    ! (from physc.f90)
!!$    ! qqq mo_sc1.f90; this is not a channel object ; do not use from here !!!
!!$    !     -> no feedback possible; since only copy
!!$    CALL new_channel_object(status, modstr,  'qte', p3=qte_3d, reprid=reprid)
    CALL new_channel_object_reference(status, 'scnbuf', 'qte'//scb &
         , modstr, 'qte')
! op_pj_20120309-
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qte', &
         'long_name', c='specific humidity tendency (qte)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qte', 'units', c='kg kg-1 s-1')
    CALL channel_halt(substr, status)
! op_pj_20120309+
    CALL get_channel_object(status, modstr, 'qte', p3=qte_3d)
    CALL channel_halt(substr, status)
! op_pj_20120309-
    ! ------------------------------------------------------------------------
! op_pj_20120309+
!!$    ! (from physc.f90)
!!$    ! qqq should be reference to xlm1 in g1a; unusable anyway
!!$    CALL new_channel_object(status, modstr,  'xlm1', p3=xlm1_3d, reprid=reprid)
    CALL new_channel_object_reference(status, 'g1a', 'xlm1', modstr, 'xlm1')
! op_pj_20120309-
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlm1', &
         'long_name', c='cloud water (xlm1)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlm1','units', c='kg kg-1')
    CALL channel_halt(substr, status)
! op_pj_20120309+
    CALL get_channel_object(status, modstr, 'xlm1', p3=xlm1_3d)
    CALL channel_halt(substr, status)
! op_pj_20120309-
    ! ------------------------------------------------------------------------
! op_pj_20120309+
!!$    ! (from physc.f90)
!!$    ! qqq mo_sc1.f90; this is not a channel object ; do not use from here !!!
!!$    !     -> no feedback possible; since only copy
!!$    CALL new_channel_object(status, modstr,  'xlte', p3=xlte_3d, reprid=reprid)
    CALL new_channel_object_reference(status, 'scnbuf', 'xlte'//scb &
         , modstr, 'xlte')
! op_pj_20120309-
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlte', &
         'long_name', c='cloud water tendency (xlte)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlte', 'units', c='kg kg-1 s-1')
    CALL channel_halt(substr, status)
! op_pj_20120309+
    CALL get_channel_object(status, modstr, 'xlte', p3=xlte_3d)
    CALL channel_halt(substr, status)
! op_pj_20120309-
    ! ------------------------------------------------------------------------
! op_pj_20120309+
!!$    ! (from physc.f90)
!!$    ! qqq should be reference to xim1 in g1a; unusable anyway
!!$    CALL new_channel_object(status, modstr,  'xim1', p3=xim1_3d, reprid=reprid)
    CALL new_channel_object_reference(status, 'g1a', 'xim1', modstr, 'xim1')
! op_pj_20120309-
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xim1', &
         'long_name', c='cloud ice (xim1)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xim1', 'units', c='kg kg-1')
    CALL channel_halt(substr, status)
! op_pj_20120309+
    CALL get_channel_object(status, modstr, 'xim1', p3=xim1_3d)
    CALL channel_halt(substr, status)
! op_pj_20120309-
    ! ------------------------------------------------------------------------
! op_pj_20120309+
!!$    ! (from physc.f90)
!!$    ! qqq mo_sc1.f90; this is not a channel object ; do not use from here !!!
!!$    !     -> no feedback possible; since only copy
!!$    CALL new_channel_object(status, modstr,  'xite', p3=xite_3d, reprid=reprid)
    CALL new_channel_object_reference(status, 'scnbuf', 'xite'//scb &
         , modstr, 'xite')
! op_pj_20120309-
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xite', &
         'long_name', c='cloud ice tendency (xite)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xite', 'units', c='kg kg-1 s-1')
! op_pj_20120309+
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, modstr, 'xite', p3=xite_3d)
! op_pj_20120309-
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); no equivalent local varibale available
    CALL new_channel_object(status, modstr,  'grmass', &
         p3=grmass, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grmass', 'long_name', c='grid mass')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grmass', 'units', c='kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! op_pj_20100713+
    CALL new_channel_object(status, modstr,  'grmassdry', &
         p3=grmassdry, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grmassdry', 'long_name' &
         , c='mass of dry air')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grmassdry', 'units', c='kg')
    CALL channel_halt(substr, status)
    ! op_pj_20100713-
    ! ------------------------------------------------------------------------
    ! (from physc.f90); no equivalent local varibale available
    CALL new_channel_object(status, modstr,  'grvol', p3=grvol, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grvol', 'long_name', c='grid volume')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grvol', 'units', c='m3')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! (from dyn.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr, 'tvirt', &
         p3=tvirt_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvirt', &
         'long_name', c='virtual temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvirt', 'units', c='degC')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from radiation.f90); not convertible
    CALL new_channel_object(status, modstr,  'relo3', &
         p3=relo3_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'relo3','long_name', c='ozone')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'relo3', 'units', c='mol mol-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from radiation.f90); not convertible
    ! qqq WARNING: only updated when l_trigrad=T, i.e. every 2 hours
    CALL new_channel_object(status, modstr, 'rhum', &
         p3=rhum_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rhum', &
         'long_name', c='relative humidity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rhum', 'units', c='%')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
!!$    ! um_ak_20100319+
!!$    ! (from physc.f90); no equivalent local varibale available
!!$    CALL new_channel_object(status, modstr,  'grmass0', &
!!$         p3=grmass0, reprid=reprid)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr, 'grmass0', 'long_name', c='grid mass')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr, 'grmass0', 'units', c='kg')
!!$    CALL channel_halt(substr, status)
!!$    ! um_ak_20100319-
    ! ------------------------------------------------------------------------

    ! ########################################################################
    ! ----------------------------- GP_3D_INT --------------------------------
    ! ########################################################################
    CALL get_representation_id(status, 'GP_3D_INT', reprid)
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! --- CHANNEL OBJECT REFERENCES ------------------------------------------
    ! ------------------------------------------------------------------------
    !
    ! ------------------------------------------------------------------------
    ! --- NEW CHANNEL OBJECTS ------------------------------------------------
    ! ------------------------------------------------------------------------
    ! (from main_data_global_start)
    ! pressure at level interfaces ("half level pressure")
    CALL new_channel_object(status, modstr,  'pressi', &
         p3=pressi_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressi', &
         'long_name', c='interface pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressi', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from physc.f90); local variable (from mo_tmp_buffer.f90)  converted to
    !                   pointer
    ! geopotential height at interface below current level
    CALL new_channel_object(status, modstr,  'geopoti', &
         p3=geopoti_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopoti', &
         'long_name', c='interface geopotential')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopoti', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from scan1.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr, 'etadot', &
         p3=etadot_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'etadot', &
         'long_name', c='vertical wind velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'etadot', 'units', c='1 s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ut_kt_20041009+
    ! (from radheat.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'flxs', &
         p3=flxs_3d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'flxs', 'long_name', c='shortwave flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'flxs', 'units', c='')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from radheat.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'flxt', &
         p3=flxt_3d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'flxt', 'long_name', c='thermal flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'flxt', 'units', c='')
    CALL channel_halt(substr, status)
    ! ut_kt_20041007-
    ! ------------------------------------------------------------------------
    ! mz_pj_20071126+
    ! (from mo_midatm.f90, gwdrag)
    IF (lmidatm) THEN
       CALL new_channel_object(status, modstr,  'gwdrag_u', &
            p3=gwdrag_u, reprid=reprid)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwdrag_u', 'long_name' &
            , c='zonal component of gravity wave drag')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwdrag_u', 'units', c='m/s^2')
       CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
       CALL new_channel_object(status, modstr,  'gwdrag_v', &
            p3=gwdrag_v, reprid=reprid)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwdrag_v', 'long_name' &
            , c='meridional component of gravity wave drag')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwdrag_v', 'units', c='m/s^2')
       CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
       CALL new_channel_object(status, modstr,  'gwflux_u', &
            p3=gwflux_u, reprid=reprid)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwflux_u', 'long_name' &
            , c='zonal component of vertical momentum flux')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwflux_u', 'units', c='Pa')
       CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
       CALL new_channel_object(status, modstr,  'gwflux_v', &
            p3=gwflux_v, reprid=reprid)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwflux_v', 'long_name' &
            , c='meridional component of vertical momentum flux')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'gwflux_v', 'units', c='Pa')
       CALL channel_halt(substr, status)
    END IF
    ! mz_pj_20071126-
    ! ------------------------------------------------------------------------

    ! ########################################################################
    ! -------------------- NON-CHANNEL-OBJECT MEMORY -------------------------
    ! ########################################################################
    ! ------------------------------------------------------------------------
    ! pxtems (from vdiff.f90); local variable converted to pointer
    ! ! long_name = vertical flux
    ! ! units     = (mol mol-1)*(kg m-2 s-1)
    ALLOCATE(pxtems(nproma, 1, ntrac_gp, ngpblks))
    ! ------------------------------------------------------------------------

  END SUBROUTINE main_data_init_memory
  !***************************************************************************

  !***************************************************************************
  SUBROUTINE main_data_local_start

    IMPLICIT NONE

    jrow = nrow(2)
    jglat = nrow(3)
    IF ( jrow == dcl%ngpblks ) THEN
       kproma = dcl%npromz
    ELSE
       kproma = dcl%nproma
    END IF

#if defined(E5302)
    CALL m_buftrow(jrow)
#endif

  END SUBROUTINE main_data_local_start
  !***************************************************************************

  !***************************************************************************
  SUBROUTINE main_data_global_start

    ! Author: Patrick Joeckel, MPICH, Nov 2004

    ! ECHAM5/MESSy
#ifndef MESSYTIMER
    USE messy_main_bmluse_bi,    ONLY: get_date_components &
                                     , get_time_step       &
                                     , next_date           &
                                     , current_date        &
                                     , get_year_day        ! mz_ab_20080309
#else
    USE messy_main_timer,        ONLY: timer_set_time_step_len &! mz_pj_20090414
                                     , lstart, lfirst_cycle
#endif 
    USE messy_main_blather_bi,   ONLY: error_bi
    ! MESSy
    USE messy_main_constants_mem, ONLY: g, M_air, R_gas, vtmpc1 &
                                      , rd ! ub_ch_20100721
    USE messy_main_tools,         ONLY: tlucua, jptlucu1, jptlucu2

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_global_start'
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: zaps
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: temp
    INTEGER                                 :: jk, jp, zjrow, it, jl
    INTEGER                                 :: zkproma
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zq, zqs
    LOGICAL                                 :: lookupoverflow = .FALSE.
!!$    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: pint   ! um_ak_20100319
!!$    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: psurf  ! um_ak_20100319
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: rho0_ll, dp0_ll ! ub_ch_20100721
    REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: hsurf           ! ub_ch_20100721

    ! mz_pj_20090414+
#ifdef MESSYTIMER
    ! set integration time step length
    CALL timer_set_time_step_len(l2tls)
#endif
    ! mz_pj_20090414-

#ifndef MESSYTIMER
    ! SET CURRENT DATE / TIME
    CALL get_date_components(current_date     &
         ,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)

    ! mz_ab_20080309+
    ! GET DAY OF YEAR (e.g. 1 Feb = 32)
    DAYOFYEAR=NINT(get_year_day(current_date))
    ! mz_ab_20080309-
#endif

#ifndef MESSYTIMER
    ! SET NEXT DATE / TIME
    CALL get_date_components(next_date     &
         ,YEAR_NEXT,MONTH_NEXT,DAY_NEXT,HOUR_NEXT,MINUTE_NEXT,SECOND_NEXT)

    ! SET CURRENT TIME STEP
    current_time_step = get_time_step()
#endif

    ! INIT
    ALLOCATE(zaps(dcl%nproma, dcl%ngpblks))
    ALLOCATE(temp(dcl%nproma, dcl%nlev, dcl%ngpblks))
    ALLOCATE(zq(dcl%nproma, dcl%nlev, dcl%ngpblks))
    ALLOCATE(zqs(dcl%nproma, dcl%nlev, dcl%ngpblks))
    
    ! CALCULATE PRESSURE
    ! WHAT ABOUT surface pressure tendency ???
    zaps(:,:)  = exp(alps_scb(:,:) + alpste_scb(:,:)*time_step_len) ! [Pa]
    ! PRESSURE AT LAYER INTERFACES
    DO jk=1, nlevp1
       pressi_3d(:,jk,:) = vct(jk) + vct(nvclev+jk) * zaps(:,:) ! [Pa]
    END DO
    ! PRESSURE AT LAYER MID
    DO jk=1, nlev
       press_3d(:,jk,:) = (pressi_3d(:,jk,:) + pressi_3d(:,jk + 1,:)) / 2.0_DP
    END DO

    ! AIR MASS IN GRID BOX
    grmass(:,:,:) = ((pressi_3d(:,2:nlevp1,:) - pressi_3d(:,1:nlev,:)) / g) &
         * SPREAD(gboxarea_2d,2,nlev)
    
!!$    ! um_ak_20100319+
!!$    ! AIR MASS FOR STANDARD ATMOSPHERE
!!$    IF (lfirst_cycle) THEN
!!$       ALLOCATE(psurf(dcl%nproma, dcl%ngpblks))
!!$       ALLOCATE(pint(dcl%nproma, dcl%nlev +1, dcl%ngpblks))
!!$       ! 1. Calculate pressure field at surface in [Pa]:
!!$       psurf(:,:) = 101325 *(1- (0.0065 * geosp(:,:)/(288.15* g)))**5.255
!!$       ! 2. Calculate interface pressure in Pa
!!$       DO jk=1, nlevp1
!!$          pint(:,jk,:) = vct(jk) + vct(nvclev+jk) * psurf(:,:)
!!$       END DO
!!$       ! 3. Calculate grid mass
!!$       grmass0(:,:,:) = ((pint(:,2:nlevp1,:) - pint(:,1:nlev,:)) / g) &
!!$            * SPREAD(gboxarea_2d,2,nlev)
!!$       DEALLOCATE(psurf)
!!$       DEALLOCATE(pint)
!!$    ENDIF
!!$    ! um_ak_20100319-

    ! METRIC VOLUME OF GRID BOX
    temp(:,:,:) = t_scb(:,:,:) + tte_scb(:,:,:)*time_step_len
    !
    DO zjrow=1, dcl%ngpblks
       IF ( zjrow == dcl%ngpblks ) THEN
          zkproma = dcl%npromz
       ELSE
          zkproma = dcl%nproma
       END IF
       grvol(1:zkproma,:,zjrow) = grmass(1:zkproma,:,zjrow) / &
            ( press_3d(1:zkproma,:,zjrow) * (1.0E-03 * M_air) &
            / (temp(1:zkproma,:,zjrow) * R_gas) )
    END DO

    ! CALCULATE RELATIVE HUMIDITY
    DO zjrow=1, dcl%ngpblks
       IF ( zjrow == dcl%ngpblks ) THEN
          zkproma = dcl%npromz
       ELSE
          zkproma = dcl%nproma
       END IF
       ! EPSILON(1.) serves to avoid water vapour content in a layer
       !          of less than EPSILON(1.).
       zq(1:zkproma,:,:)= &
            MAX( qm1(1:zkproma,:,:)+time_step_len*qte_scb(1:zkproma,:,:) &
            , EPSILON(1._dp) )
       !
       DO jk = 1, nlev
          DO jl = 1, zkproma
             it = INT(temp(jl,jk,zjrow)*1000._dp)
             IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
             it = MAX(MIN(it,jptlucu2),jptlucu1)
             zqs(jl,jk,zjrow) = tlucua(it)/press_3d(jl,jk,zjrow)
          END DO
       END DO
       IF (lookupoverflow) &
            CALL error_bi('lookuperror', substr)
       zqs(1:zkproma,:,zjrow)= MIN(zqs(1:zkproma,:,zjrow),0.5_dp)
       zqs(1:zkproma,:,zjrow)= zqs(1:zkproma,:,zjrow) &
            / (1._dp-vtmpc1*zqs(1:zkproma,:,zjrow))
       zqs(1:zkproma,:,zjrow)= MAX(2._dp*EPSILON(1._dp),zqs(1:zkproma,:,zjrow))
       rhum_3d(1:zkproma,:,zjrow) = 100._dp * zq(1:zkproma,:,zjrow) &
            / zqs(1:zkproma,:,zjrow)
    END DO

    ! op_pj_20100713+
    grmassdry(:,:,:) = grmass(:,:,:) * (1.0_dp - zq(:,:,:))
    ! op_pj_20100713-

    ! CLEAN UP
    DEALLOCATE(zaps)
    DEALLOCATE(temp)
    DEALLOCATE(zq)
    DEALLOCATE(zqs)

  END SUBROUTINE main_data_global_start
  !***************************************************************************

  !***************************************************************************
  SUBROUTINE main_data_free_memory

    IMPLICIT NONE

    ! CLEANUP NON-CHANNEL OBJECT MEMORY
    DEALLOCATE(pxtems)

  END SUBROUTINE main_data_free_memory
  !***************************************************************************

  ! -------------------------------------------------------------------
  ! PRIVATE HELPER ROUTINES 
  !  -> SAME AS IN MESSY_MAIN_CHANNEL_BI; DUE TO CIRCULAR DEPENDENCIES
  ! -------------------------------------------------------------------
  SUBROUTINE channel_halt(substr, status)

    ! MESSy
    USE messy_main_blather_bi,     ONLY: error_bi, info_bi
    USE messy_main_constants_mem,  ONLY: STRLEN_VLONG
    USE messy_main_channel_error,  ONLY: channel_error_str

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: substr
    INTEGER,          INTENT(IN)  :: status
    ! LOCAL
    CHARACTER(LEN=STRLEN_VLONG)   :: errstr

    IF (status == 0) RETURN

    errstr = channel_error_str(status)

    CALL error_bi(errstr, substr)

  END SUBROUTINE channel_halt
  ! -------------------------------------------------------------------

#if defined(E5302)
  ! -------------------------------------------------------------------
  SUBROUTINE m_buftrow(krow)

    USE messy_main_tracer_mem_bi, ONLY: xtte_scb=>xtte, xtte_a_scb=>xtte_a &
                                      , L_LG ! op_pj_20120120

    IMPLICIT NONE
    SAVE

    INTEGER, INTENT(IN) :: krow

    vo            => vo_scb(:,:,krow)
    d             => d_scb(:,:,krow)
    t             => t_scb(:,:,krow)
    alps          => alps_scb(:,krow)
    u             => u_scb(:,:,krow)
    dudl          => dudl_scb(:,:,krow)
    v             => v_scb(:,:,krow)
    dvdl          => dvdl_scb(:,:,krow)
    vol           => vol_scb(:,:,krow)
    vom           => vom_scb(:,:,krow)
    rh            => rh_scb(:,:,krow)
    qte           => qte_scb(:,:,krow)
    xlte          => xlte_scb(:,:,krow)
    xite          => xite_scb(:,:,krow)
    xtte          => xtte_scb(:,:,:,krow)
    ! mz_pj_20030930+
    ! KROW = 1 -> INIT AT START OF THE LOOP
    IF (L_LG) THEN  ! op_pj_20120120
       IF (krow == 1) THEN
          xtte_a     => xtte_a_scb(:,:,:,1)
       ENDIF
    END IF          ! op_pj_20120120
    ! mz_pj_20030930-
    tte           => tte_scb(:,:,krow)
    alpste        => alpste_scb(:,krow)
    alnpr         => alnpr_scb(:,:,krow)
    alpha         => alpha_scb(:,:,krow)
    vervel        => vervel_scb(:,:,krow)

    ! mo_geoloc
    gboxarea      => gboxarea_2d(:,krow)

  END SUBROUTINE m_buftrow
  ! -------------------------------------------------------------------
#endif

!!! END OF #if defined(ECHAM5)
#endif 

!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
  
#ifdef COSMO
  
  ! Author: Astrid Kerkweg,  IPA, Uni Mainz,  2008
  ! 
  ! the Module messy_main_data_bi 
  ! a) transfers data from the COSMO model
  !    to the messy submodel interface layer. I.e., especially it "renames"
  !    some fields in a way, that they are accessable under the same name 
  !    in the submodel interface layer
  ! b) serves as interface layer to provide information from the base model
  !    to the MESSy (e.g. the grid resolution etc.)
  ! Note: All variables associated with the time management are used by 
  !       messy_main_timer_bi directly  from the COSMO model.
  !       And all time information within MESSy is provided by MAIN_TIMER.
  !       Accordingly all variables providing information for the MPI 
  !       environment are use/provided by messy_main_mpi_bi

  ! COSMO
  USE data_runcontrol,   ONLY: lmulti_layer, psm0, dsem0, msem0, kem0, qcm0 &
                             , l2tls, nnow, nnew, nold, ntke, hnextrad
  USE data_modelconfig,  ONLY: ie, je, ke, ie_tot, je_tot, ke_tot &
       , ie_max, je_max, dlon, dlat, startlon_tot, startlat_tot   &
       , istartpar, iendpar, jstartpar, jendpar, ivctype ,dt      &
       , ke_soil, czmls, kcm, pollon, pollat, polgam , vcoord     &
       , sigmr,  hhlr, t0sl, p0sl, vcflat, dt0lp, delta_t, h_scal &
       , svc1, svc2, nfltvc, irefatm 
  USE data_io,           ONLY: root, pp_nl, var, ngribout
  USE data_fields,       ONLY:   &
       ! 1. constant fields for the reference atmosphere      (unit)
       ! -----------------------------------------------
       rho0       , & ! reference density at the full model levels     (kg/m3)
       dp0        , & ! pressure thickness of model layers             ( Pa  )
       p0         , & ! reference pressure at full levels              ( Pa  )
       p0hl       , & ! reference pressure at half levels              ( Pa  )
       dt0dz      , & ! temperature gradient of reference atmosphere   ( K/m )
       t0         , & ! reference temperature                          ( K   )
       hhl        , & ! geometical height of half levels               (  m  )
       sqrtg_r_s  , & ! reciprocal square root of G at skalar points   ( 1/m )
       sqrtg_r_u  , & ! reciprocal square root of G at u points        ( 1/m )
       sqrtg_r_v  , & ! reciprocal square root of G at v points        ( 1/m )
       sqrtg_r_w      ! reciprocal square root of G at w points        ( 1/m )

  USE data_fields     , ONLY :   &
       ! 2. external parameter fields                                  (unit)
       ! ----------------------------
       hsurf      , & ! geometical heigt of surface topography         (  m  )
       sso_stdh   , & ! standard deviation of sub-grid scale orography( m   )
       sso_gamma  , & ! anisotropy of sub-grid scale orography          --
       sso_theta  , & ! angle betw. principal axis of orography and E ( rad )
       sso_sigma  , & ! mean slope of sub-grid scale orography          --
       gz0        , & ! surface roughness  * g                         (m2/s2)
       fr_land    , & ! fraction of land in a grid element             (  -- )
       soiltyp    , & ! type of the soil (keys 0-9)                    (  -- )
       vio3       , & ! vertical integrated ozone contents             (pa O3)
       hmo3       , & ! ozone maximum                                  ( pa  )
       rlat       , & ! geographical latitude                          ( rad )
       rlon       , & ! geographical longitude                         ( rad )
       rlattot    , & ! geographical latitude                          ( rad )
       rlontot    , & ! geographical longitude                         ( rad )
       fc         , & ! coriolis-parameter                             ( 1/s )
       fccos      , & ! coriolis-parameter mit cos(phi)                ( 1/s )
       rmy        , & ! Davis-parameter for relaxation (mass, qv, qc)    --
       rmyq       , & ! Davis-parameter for relaxation (qr, qs, qg)      --
       hd_mask_dcoeff, & ! 3D-domain mask for horizontal diffusion * dcoeff --
       hd_mask    , & ! 3D-domain mask for horizontal diffusion * dcoeff --
       least_lbdz , & ! mask for eastern  lateral boundary zone
       lwest_lbdz , & ! mask for western  lateral boundary zone
       lnorth_lbdz, & ! mask for northern lateral boundary zone
       lsouth_lbdz, & ! mask for southern lateral boundary zone
       crlat      , & ! cosine of transformed latitude
       acrlat     , & ! 1 / ( crlat * radius of the earth )            ( 1/m )
       tgrlat     , & ! tangens of transformed latitude                  --
       aerlan     , & ! aerosol-distribution for rural areas             --
       aerurb     , & ! aerosol-distribution for urban areas             --
       aerdes     , & ! aerosol-distribution for desert areas            --
       aersea     , & ! aerosol-distribution for sea                     --
       plcov      , & ! fraction of plant cover                          --
       lai        , & ! leaf area index of plants                        --
       tai        , & ! transpiration area index                         --
       sai        , & ! surface area index                               --
       eai        , & ! (evaporative) earth area index                   --
       rootdp     , & ! depth of the roots                             (  m  )
       llandmask      ! landpoint mask
  USE data_fields     , ONLY :   &
       for_e      , & ! ground fraction covered by evergreen forest      --
       for_d      , & ! ground fraction covered by deciduous forest      --
       h_can      , & ! hight of the vertically resolved canopy        (  m  )
       d_pat      , & ! horizontal pattern length scale                (  m  )
       c_big      , & ! effective drag coefficient of canopy elements
                      ! larger than or equal to the tubulent lengthscale( 1/m )
       c_sml      , & ! effective drag coefficient of canopy elements
                      ! smaller than the tubulent length scale         ( 1/m )
       r_air      , & ! air containing fraction of a gridbox inside
                      ! the canopy                                     (  1  )
       fr_lake    , & ! lake fraction in a grid element [0,1]          (  -- )
       depth_lk   , & ! lake depth                                     (  m  )
       fetch_lk   , & ! wind fetch over lake                           (  m  )
       dp_bs_lk   , & ! depth of the thermally active layer
                      ! of bottom sediments                            (  m  )
       t_bs_lk    , & ! climatological temperature at the bottom of
                      ! the thermally active layer of sediments        (  K  )
       gamso_lk   , & ! attenuation coefficient for
                      ! solar radiation in lake water                  ( 1/m )

       ! 3. prognostic variables                                        (unit)
       ! -----------------------
       u          , & ! zonal wind speed                               ( m/s )
       v          , & ! meridional wind speed                          ( m/s )
       w          , & ! vertical wind speed (defined on half levels)   ( m/s )
       t          , & ! temperature                                    (  k  )
       qv         , & ! specific water vapor content                   (kg/kg)
       qc         , & ! specific cloud water content                   (kg/kg)
       qi         , & ! specific cloud ice   content                   (kg/kg)
       qr         , & ! specific rain content                          (kg/kg)
       qs         , & ! specific snow content                          (kg/kg)
       qg         , & ! specific graupel content                       (kg/kg)
       pp         , & ! deviation from the reference pressure          ( pa  )
       tke        , & ! SQRT(2 * turbulent kinetik energy)             ( m/s )
       ! (defined on half levels)
       edr,         & ! eddy dissipation rate of TKE (EDR)             (m2/s3)
       ! (defined on half levels)

       ! 4. tendency fields for the prognostic variables               (unit )
       ! -----------------------------------------------
       !    timely deviation  by diabatic and adiabatic processes 
       !    without sound-wave terms
       utens        ,  & ! u-tendency without sound-wave terms         ( m/s2)
       vtens        ,  & ! v-tendency without sound-wave terms         ( m/s2)
       wtens        ,  & ! w-tendency without sound-wave terms         ( m/s2)
                         ! (defined on half levels )
       ttens        ,  & ! t-tendency without sound-wave terms         ( m/s2)
       qvtens       ,  & ! qv-tendency                                 ( m/s2)
       qctens       ,  & ! qc-tendency                                 ( m/s2)
       qitens       ,  & ! qi-tendency                                 ( m/s2)
       pptens       ,  & ! pp-tendency without sound-wave terms        ( m/s2)
       tketens           !tke-tendency (defined on half-levels)        ( m/s )

  USE data_fields     , ONLY :&
       ! 5. fields for surface values and soil model variables         (unit )
       ! -----------------------------------------------------
       ps        ,  & ! surface pressure                               ( pa  )
       t_snow    ,  & ! temperature of the snow-surface                (  k  )
       t_s       ,  & ! temperature of the ground surface              (  k  )
       t_s_lake  ,  & ! temperature of the ground surface (lake)       (  k  )
       t_g       ,  & ! weighted surface temperature                   (  k  )
       qv_s      ,  & ! specific water vapor content on the surface    (kg/kg)
       t_m       ,  & ! temperature between upper and medium 
                      ! soil layer                                     (  k  )
       t_cl      ,  & ! temperature between medium and lower
                      ! soil layer                                     (  k  )
       t_so      ,  & ! multi-layer soil temperature                   (  k  )
       w_snow    ,  & ! water content of snow                          (m H2O)
       w_i       ,  & ! water content of interception water            (m H2O)
       w_g1      ,  & ! water content of the upper soil layer          (m H2O)
       w_g2      ,  & ! water content of the medium soil layer         (m H2O)
       w_g3      ,  & ! water content of the lower soil layer          (m H2O)
                      ! (if NLWB=3, unused otherwise)
       w_so      ,  & ! multi-layer soil moisture                      (m H2O)
       w_so_ice  ,  & ! multi-layer soil ice                           (m H2O)
       w_cl      ,  & ! climatological water content                   (m H2O) 
       freshsnow ,  & ! weighting function indicating 'freshness' of snow in
                      ! upper few centimeters of snow cover            ( -- )
       rho_snow  ,  & ! prognostic density of snow                     (kg/m3)
       h_snow    ,  & ! snow height                                    (  m  )
       t_e       ,  & ! surface temperature of the canopy elements     (  k  )
       qv_e      ,  & ! surface value of qd of the canopy elements     (Kg/Kg)
       t_ice     ,  & ! temperature at the snow-ice or
                      ! air-ice interface                              (  K  )
       t_mnw_lk  ,  & ! mean temperature of the water column           (  K  )
       t_wml_lk  ,  & ! mixed-layer temperature                        (  K  )
       t_bot_lk  ,  & ! temperature at the water-bottom sediment
                      ! interface                                      (  K  )
       t_b1_lk   ,  & ! temperature at the bottom of the upper layer
                      ! of the sediments                               (  K  )
       c_t_lk    ,  & ! shape factor with respect to the
                      ! temperature profile in lake thermocline        (  -  )
       h_ice     ,  & ! ice thickness                                  (  m  )
       h_ml_lk   ,  & ! thickness of the mixed-layer                   (  m  )
       h_b1_lk        ! thickness of the upper layer of bottom sediments(  m  )

  !--------------------------------------------------------------------------
  USE data_fields     , ONLY :&
       ! 6. fields that are computed in the parametrization and dynamics (unit)
       ! ---------------------------------------------------------------
       qvt_diff   , & ! humidity    tendency  due to diffusion         ( 1/s )
       tinc_lh    , & ! temperature increment due to latent heat       (  K  )
                      !   density of moist air 
       rho        , & ! total density of moist air                     (kg/m3)
       
       !   coefficients for turbulent diffusion in the atmosphere
       !   (defined on half levels)
       ! vertical   turbulent diffusion coefficients
       tkvm     ,   & ! ... for momentum                               (m2/s)
       tkvh     ,   & ! ... for heat and moisture                      (m2/s)
       ! horizontal turbulent diffusion coefficients
       tkhm     ,   & ! ... for momentum                               (m2/s)
       tkhh     ,   & ! ... for heat and moisture                      (m2/s)
       
       !   turbulence statistics in the atmosphere
       !   (defined on full levels)
       rcld       , & ! standard deviation of the saturation deficit      --
  
       !   turbulent coefficients at the surface 
       tcm      ,   & ! transfer coefficient for momentum               ( -- )
       tch      ,   & ! transfer coefficient for heat and moisture      ( -- )
       tfm      ,   & ! factor of laminar transfer of momentum            --
       tfh      ,   & ! factor of laminar transfer of scalars             --
       tfv      ,   & ! laminar reduction factor for evaporation        ( -- )
       
       !   fields from the radiation scheme
       sohr      ,  & ! rate of solar heating                          ( K/s )
       thhr      ,  & ! rate of thermal heating                        ( K/s )
       clc_sgs   ,  & ! subgrid-scale stratiform cloud cover              --
       alb_rad   ,  & ! albedo of the ground                              --
       sobs      ,  & ! solar radiation at the ground                  ( w/m2)
       thbs      ,  & ! thermal radiation at the ground                ( w/m2)
       pabs      ,  & ! photosynthetic active radiation at the ground  ( w/m2)
       sobt      ,  & ! solar radiation at the upper boundary          ( w/m2)
                      ! of the atmosphere
       thbt      ,  & ! thermal radiation at the upper boundary        ( w/m2)
                      ! of the atmosphere
       clch      ,  & ! cloud cover with high clouds                      --   
       clcm      ,  & ! cloud cover with medium clouds                    --   
       clcl      ,  & ! cloud cover with low clouds                       --   
       clct      ,  & ! total cloud cover                                 --   
       
       ! and used in the Climate-LM Version
       sodwddm   ,  & ! downward direct solar radiative flux / smu0    ( W/m2)
       qc_rad    ,  & ! subgrid-scale specific cloud water             (kg/kg)
       qi_rad         ! subgrid-scale specific ice water               (kg/kg)

  USE data_fields     , ONLY :&
       !   fields for the radiation correction scheme
       swdir_s  ,   & ! direct comp. of solar radiative flux at surface ( W/m2)
       swdifd_s ,   & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
       swdifu_s ,   & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
       lwd_s    ,   & !         downward comp. of long  wave rad. flux  ( W/m2)
       lwu_s    ,   & !         upward   comp. of long  wave rad. flux  ( W/m2)
       
       ! these are accumulated values
       aswdir_s ,   & ! direct comp. of solar radiative flux at surface ( W/m2)
       aswdifd_s,   & ! diffuse downward comp. of short wave rad. flux  ( W/m2)
       aswdifu_s,   & ! diffuse upward   comp. of short wave rad. flux  ( W/m2)
       alwd_s   ,   & !         downward comp. of long  wave rad. flux  ( W/m2)
       alwu_s   ,   & !         upward   comp. of long  wave rad. flux  ( W/m2)
       
       ! this is the essential correction factor
       swdir_cor,   & ! direct short wave radiation correction factor
                      ! actual value
       
       ! these are topographic parameters
       skyview  ,   & ! sky view
       slo_asp  ,   & ! slope aspect
       slo_ang  ,   & ! slope angle
       horizon           ! horizon

  !--------------------------------------------------------------------------
  USE data_fields     , ONLY :&
       !   fields from the convection scheme
       clc_con    , & ! cloud cover due to convection                     --
       clw_con    , & ! cloud liquid water due to convection              --
       prr_con    , & ! precipitation rate of rain, convective        (kg/m2s)
       prs_con    , & ! precipitation rate of snow, convective        (kg/m2s)
       prne_con   , & ! precipitation rate, no evaporat., convective  (kg/m2s)
       bas_con    , & ! level index of convective cloud base            --
       top_con    , & ! level index of convective cloud top             --
       tt_conv    , & ! temperature tendency due to convection        ( K/s  )
       qvt_conv   , & ! humidity    tendency due to convection        ( 1/s  )
       qct_conv   , & ! qc-tendency tendency due to convection        ( 1/s  )
       qit_conv   , & ! qi-tendency tendency due to convection        ( 1/s  )
       qrt_conv   , & ! qr-tendency tendency due to convection        ( 1/s  )
       qst_conv   , & ! qs-tendency tendency due to convection        ( 1/s  )
       ut_conv    , & ! u-tendency due to convection                  ( m/s^2)
       vt_conv    , & ! v-tendency due to convection                  ( m/s^2)
       mflx_con   , & ! convective massflux                           (kg/m2s)
       cape_con   , & ! convective available energy                   ( J/kg )
       tke_con    , & ! convective turbulent kinetic energy           ( J/kg )
       qcvg_con   , & ! moisture convergence for Kuo-type closure     ( 1/s  )
       w0avg      , &
       nca        , &
       
       !   fields of the precipitation
       qrs        , & ! precipitation water (water loading)           (kg/kg )
       prr_gsp    , & ! precipitation rate of rain, grid-scale        (kg/m2s)
       prs_gsp    , & ! precipitation rate of snow, grid-scale        (kg/m2s)
       prg_gsp    , & ! precipitation rate of graupel, grid-scale     (kg/m2s)
       
       !   fields that are computed in the dynamics
       dqvdt      , & ! threedimensional moisture convergence         ( 1/s  )
       qvsflx     , & ! surface flux of water vapour                  (1/m2s )
       dpsdt      , & ! tendency of the surface pressure              ( pa/s )
       umfl_s     , & ! u-momentum flux (surface)                     ( N/m2 )
       vmfl_s     , & ! v-momentum flux (surface)                     ( N/m2 )
       shfl_s     , & ! sensible heat flux (surface)                  ( W/m2 )
       lhfl_s     , & ! latent heat flux (surface)                    ( W/m2 )
       aumfl_s    , & ! average u-momentum flux (surface)             ( N/m2 )
       avmfl_s    , & ! average v-momentum flux (surface)             ( N/m2 )
       ashfl_s    , & ! average sensible heat flux (surface)          ( W/m2 )
       alhfl_s        ! average latent heat flux (surface)            ( W/m2 )

  !-------------------------------------------------------------------------

  USE data_fields     , ONLY :&
       ! 7. fields for model output and diagnostics                    (unit)
       ! ------------------------------------------
       t_2m       , & ! temperature in 2m                             (  K   )
       t_2m_av    , & ! time mean temperature in 2m                   (  K   )
       qv_2m      , & ! specific water vapor content in 2m            (kg/kg )
       td_2m      , & ! dew-point in 2m                               (  K   )
       td_2m_av   , & ! time mean dew-point in 2m                     (  K   )
       rh_2m      , & ! relative humidity in 2m                       (  %   )
       u_10m      , & ! zonal wind in 10m                             ( m/s  )
       u_10m_av   , & ! time mean zonal wind in 10m                   ( m/s  )
       v_10m      , & ! meridional wind in 10m                        ( m/s  )
       v_10m_av   , & ! time mean meridional wind in 10m              ( m/s  )
       tmin_2m    , & ! minimum temperature in 2m                     (  K   )
       tmax_2m    , & ! maximum temperature in 2m                     (  K   )
       vmax_10m   , & ! maximal windspeed in 10m                      ( m/s  )
       vgust_con  , & ! maximal convective wind gust in 10m           ( m/s )
       vgust_dyn  , & ! maximal dynamical wind gust in 10m            ( m/s )
       asob_s     , & ! average solar radiation budget (surface)      ( W/m2 )
       athb_s     , & ! average thermal radiation budget (surface)    ( W/m2 )
       apab_s     , & ! average photosynthetic active radiation (sfc) ( W/m2 )
       asob_t     , & ! average solar radiation budget (model top)    ( W/m2 )
       athb_t     , & ! average thermal radiation budget (model top)  ( W/m2 )
       sod_t      , & ! solar downward radiation at top of atmosphere (      )
       asod_t     , & ! averaged solar downward radiation at top      (      )
       dursun     , & ! sunshine duration                             (  s   )
       rain_gsp   , & ! amount of rain from grid-scale precip. (sum)  (kg/m2 )
       snow_gsp   , & ! amount of snow from grid-scale precip. (sum)  (kg/m2 )
       grau_gsp   , & ! amount of graupel from grid-scale prec. (sum) (kg/m2 )
       rain_con   , & ! amount of rain from convective precip. (sum)  (kg/m2 )
       snow_con   , & ! amount of snow from convective precip. (sum)  (kg/m2 )
       runoff_s   , & ! surface water runoff; sum over forecast       (kg/m2 )
       runoff_g   , & ! soil water runoff; sum over forecast          (kg/m2 )
       snow_melt ,  & ! amount of snow melt; sum over forecast        (kg/m2)
       tdiv_hum   , & ! vertical integral divergence of humidity      (kg/m2 )
       aevap_s        ! accumulated surface moisture flux             (kg/m2 )

  !------------------------------------------------------------------------------

  USE data_fields     , ONLY :&
       ma_usl_kmin, & ! level index of lower boundary of updraft source layer
       ma_usl_kmax, & ! level index of upper boundary of updraft source layer
       ma_lfs     , & !
       ma_etl_k   , & !
       ma_ml      , & !
       ma_ddt     , & !
       ma_usl_buoy, & ! buoyancy at start of ascent (T only)
       ma_nb_k    , & ! level index of neutral buoyancy below the LCL (T only)
       ma_nb_k_min, & ! minimum of ma_nb_k in last hour (T only)
       ma_lcl_k   , & ! level index of lifting condensation level (LCL)
       ma_lcl_t   , & ! parcel temperature at the LCL
       ma_lcl_dt  , & ! temperature perturbation at the LCL
       ma_lcl_tenv, & ! environment temperature at the LCL
       ma_trg     , & ! trigger criteria (degrees)
       ma_trg_max , & ! maximum of ma_trg in last hour
       ma_top_k   , & ! level index of cloud top
       ma_type    , & ! type of convection (1=deep, 2=shallow, 3=mid-level)
       ma_umf     , & ! updraft mass flux
       ma_udr     , & ! updraft detrainment rate
       ma_uer     , & ! updraft entrainment rate
       ma_urv     , & ! water vapour in updraft
       ma_urci    , & ! total condensat in updraft
       ma_ls_rad  , & ! updraft area
       ma_uw      , & ! velocity in updraft
       ma_wsub    , & ! compensating mass flux in environment
       ma_dmf     , & ! downdraft mass flux
       ma_der     , & ! downdraft entrainment rate
       ma_ddr     , & ! downdraft detrainment rate
       ma_drw     , & ! total downdraft water
       ma_prlflx  , & ! liquid precipitation flux
       ma_prsflx  , & ! solid precipitation flux
       ma_urr     , & ! liquid precipitation produced in model layer
       ma_urs         ! solid precipitation produced in model layer

  USE data_fields     , ONLY :&
       ! 8. fields for the boundary values                           (unit)
       ut_sso    , & ! u-tendency due to SSO                         ( m/s2)
       vt_sso    , & ! v-tendency due to SSO                         ( m/s2)
       tt_sso    , & ! temperature tendency due to SSO               ( K/s )
       ustr_sso  , & ! u-stress (surface momentum flux) due to SSO   ( N/m2)
       vstr_sso  , & ! v-stress (surface momentum flux) due to SSO   ( N/m2)
       vdis_sso  , & ! vert. int. dissipation of kin. en. due to SSO ( W/m2)
       austr_sso , & ! average of ustr_sso                           ( N/m2)
       avstr_sso , & ! average of vstr_sso                           ( N/m2)
       avdis_sso     ! average of vdis_sso                           ( W/m2)
  !---------------------------------------------------------------------------

  USE data_fields     , ONLY :&
       ! 8. fields for the boundary values                               (unit)
       ! ---------------------------------
       u_bd       , & ! boundary field for u                          ( m/s  )
       v_bd       , & ! boundary field for v                          ( m/s  )
       w_bd       , & ! boundary field for w                          ( m/s  )
       t_bd       , & ! boundary field for t                          (  k   )
       qv_bd      , & ! boundary field for qv                         (kg/kg )
       qc_bd      , & ! boundary field for qc                         (kg/kg )
       qi_bd      , & ! boundary field for qi                         (kg/kg )
       qr_bd      , & ! boundary field for qr                         (kg/kg )
       qs_bd      , & ! boundary field for qs                         (kg/kg )
       qg_bd      , & ! boundary field for qg                         (kg/kg )
       pp_bd      , & ! boundary field for pp                         (  pa  )
       qv_s_bd    , & ! boundary field for qv_s                       (kg/kg )
       t_snow_bd  , & ! boundary field for t_snow                     (  k   )
       t_s_bd     , & ! boundary field for t_s                        (  k   )
       t_m_bd     , & ! boundary field for t_m                        (  k   )
       w_snow_bd  , & ! boundary field for w_snow                     (m H2O )
       w_g1_bd    , & ! boundary field for w_g1                       (m H2O )
       w_g2_bd    , & ! boundary field for w_g2                       (m H2O )
       w_g3_bd    , & ! boundary field for w_g3                       (m H2O )
       hmo3_bd    , & ! boundary field for hmo3                       (m    )
       vio3_bd    , & ! boundary field for vio3                       (pa O3)
       w_cl_bd    , & ! boundary field for w_cl                       (m H2O)
       t_cl_bd    , & ! boundary field for t_cl                       (  K  )
       lai_bd     , & ! boundary field for lai                        ( --  )
       rootdp_bd  , & ! boundary field for rootdp                     (m    )
       plcov_bd   , & ! boundary field for plcov                      ( --  )
       
       ! 10. analysis increment fields
       ! -----------------------------
       ff_anai    , & ! wind velocity                                 ( m/s )
       dd_anai    , & ! wind direction                                ( rad )
       t_anai     , & ! temperature                                   (  k  )
       p_anai     , & ! deviation from the reference pressure         ( Pa  )
       qv_anai    , & ! specific water vapor content                  (kg/kg)
       qc_anai        ! specific cloud water content (via saturation adjustm)
  USE data_fields,  ONLY:  &
       synme7,      & ! Meteosat 7
       synmsg         ! Meteosat Second Generation

    USE data_lheat_nudge,  ONLY :   &
         llhn         ,& ! on/off switch for latent heat nudging (lhn)
         llhnverif    ,& ! on/off switch for latent heat nudging (lhn)
         lhn_qrs      ,& ! use integrated precipitation flux as reference
         tt_lheat     ,& ! profile of t-increments due to latent heating (K/s )
         qrsflux         ! total precipitation flux
    
    !   fields for lh-diagnostics
    USE data_lhn_diag, ONLY  :  &
         tinc_lhn_o  ,& ! temperature increments due to lhn  ( K/s  )
         tt_lheat_o  ,& ! array for cumulated latent heating
                        ! (grid scale + conv)( K )
         ttm_cv_o       ! array for test output of diverse 2D fields

    USE data_io,        ONLY :   &
         llb_qi,       & ! if .TRUE., take qi_bd-values from lateral boundaries
                         ! file else, qi_bd is set in the model
         llb_qr_qs,    & ! if .TRUE., take qr_bd- and qs_bd-values from lateral
                         ! bound. file else, qr_bd and qs_bd are set in model
         llb_qg,       & ! if .TRUE., take qg_bd-values from lateral boundaries
                         ! file else, qg_bd is set in the model
         lbdclim         ! boundary data in climate model

    USE data_runcontrol , ONLY :   &
         ! 3. controlling the physics
         ! --------------------------
         itype_gscp,   & ! type of grid-scale precipitation physics
         itype_turb,   & ! type of turbulent diffusion parametrization
         nlgw,         & ! number of prognostic soil water levels
         nlgw_ini,     & ! number of progn. soil water levels in initial data
         nlgw_bd,      & ! number of progn. soil water levels in boundary data
         ltur,         & ! forecast with vertical diffusion
         l3dturb,      & ! 3D-turbulence (additional horizontal diffusion)
         lprog_tke,    & ! prognostic treatment of TKE (for itype_turb=5/7)
         lphys,        & ! forecast with physical parametrizations
         itype_conv,   & ! type of convection parameterization
         llake,        & ! forecast with lake model
         lsso,         & ! forecast with sub-grid scale orography scheme
         lforest,      & ! if .true., run with forest (evergreen and deciduous)
         lprogprec,    & ! forecast with prognostic rain and snow (qr, qs)
         ! 4. controlling the dynamics
         ! ---------------------------
         lcori_deep,    & ! if =.TRUE.: account for cos(phi) coriolis terms
         ! 7. additional control variables
         ! -------------------------------
         lprog_qi,     & ! if .TRUE., running with cloud ice
         irunge_kutta, & ! type of Runge-Kutta scheme
         ldiabf_lh,    & ! include diabatic forcing due to latent heat in
                         ! RK-scheme
         lw_freeslip,  & ! if .TRUE.: with free slip lateral boundary condition
                         ! and if .FALSE. specified lateral bound. values for w
         lreproduce,   & ! the results are reproducible in parallel mode
         lradtopo,     & ! if .TRUE., calculate topographic correction of
                         ! radiation
         ! 8. diagnostic calculations
         ! --------------------------
         ldiagnos,     & ! perform diagnostic calculations
         ! 9. Other variables
         ! ------------------
         lout_anai,    & ! allocate fields to enable writing analysis increments
         luse_rttov      ! calculate satellite images
 ! COSMO/MESSy
  USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_abort , p_pe
  USE messy_main_blather_bi, ONLY: warning_bi

  ! MESSy
  USE messy_main_constants_mem, ONLY: dp

  !************************************************************************
  ! RENAMING USE STATEMENTS FOR USAGE IN MESSY-SUBMODELS
  !************************************************************************
  USE data_fields,  ONLY: slf       => fr_land  & ! sea-land fraction
                        , alake     => fr_lake  & ! lake fraction
                        , prl       => rain_gsp & ! rain from grid-scale precip.
                        , prc       => rain_con & ! rain from convective precip.
                        , vgrat     => plcov    & ! fraction of plant cover
                        , srfl      => sobs     & ! net surface radiative flux
                        , coriol_2d => fc       & ! coriolis coefficient
                        , albedo    => alb_rad  & ! surface albedo
                        , xite_3d   => qitens   & ! ice tendency
                        , xlte_3d   => qctens   & ! liquid water tendency
                        , tte_3d    => ttens    & ! temperature tendency
                        , qte_3d    => qvtens   &  ! water vapour tendency
                        , u10 => u_10m          &
                        , v10 => v_10m          !&
                        !!, seaice => fr_ice    not available at the time being

  IMPLICIT NONE
  PUBLIC
  SAVE

  INTRINSIC :: NULL

  ! COSMO PARAMETERS
  CHARACTER(LEN=*), PARAMETER :: modstr = 'COSMO'

  CHARACTER(LEN=*), PARAMETER :: modver = '4.8_clm12'

  REAL(dp) :: eps  = 0.1_dp

  INTEGER  :: ntime ! index for time level
                    ! l2tls = T => ntime = nnow
                    ! l2tls = F => ntime = nnold

  !****************************************************************************
  !****************************************************************************
  ! POINTER FOR CHANNEL OBJECTS
  !****************************************************************************
  !****************************************************************************
  REAL(dp), POINTER, DIMENSION(:,:,:) :: press_3d   => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: pressi_3d  => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: geopot_3d  => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: geopoti_3d => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: grmass     => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: grmassdry  => NULL() ! op_pj_20100713
  REAL(dp), POINTER, DIMENSION(:,:,:) :: grvol      => NULL() 

  REAL(dp), POINTER, DIMENSION(:,:)   :: crlat_2d    => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:)   :: crlati_2d   => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:)   :: gboxarea_2d => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:)   :: philat_2d   => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:)   :: philon_2d   => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:)   :: coslon_2d   => NULL() ! op_pj_20110714
  REAL(dp), POINTER, DIMENSION(:,:)   :: sinlon_2d   => NULL() ! op_pj_20110714
  REAL(dp), POINTER, DIMENSION(:,:)   :: coslat_2d   => NULL() ! op_pj_20110714
  REAL(dp), POINTER, DIMENSION(:,:)   :: sinlat_2d   => NULL() ! op_pj_20110714

  ! PROGNOSTIC VARIABLES
  REAL(dp), POINTER, DIMENSION(:,:,:) :: tm1_3d  => NULL() 
  !REAL(dp), POINTER, DIMENSION(:,:,:) :: tte_3d  => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qm1_3d  => NULL() 
  !REAL(dp), POINTER, DIMENSION(:,:,:) :: qte_3d  => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qc_3d   => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qi_3d   => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qr_3d   => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qs_3d   => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: qg_3d   => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: um1     => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: vm1     => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: wm1     => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: ppm1_3d => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: t_so_3d => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: t_so_ini => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: w_so_3d => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:,:) :: xim1_3d => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:) :: xlm1_3d => NULL()
  ! DIAGNOSTIC VARIABLES
  ! potential temperatur
  REAL(dp), POINTER, DIMENSION(:,:,:) :: tpot_3d => NULL()
  ! vorticity
  REAL(dp), POINTER, DIMENSION(:,:,:) :: vom1    => NULL()
  ! relative humidity
  REAL(dp), POINTER, DIMENSION(:,:,:) :: rhum_3d => NULL()
  
  ! 2D TLV FIELDS
  REAL(dp), POINTER, DIMENSION(:,:)   :: aps       => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:)   :: tsurf_2d  => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: t_snow_2d => NULL() ! bd
  REAL(dp), POINTER, DIMENSION(:,:)   :: w_snow_2d => NULL() ! bd
  REAL(dp), POINTER, DIMENSION(:,:)   :: w_i_2d    => NULL() 
  REAL(dp), POINTER, DIMENSION(:,:)   :: qv_s_2d   => NULL() ! bd
  ! DECOMPOSITION DIAGNOSTIC
  REAL(dp), POINTER, DIMENSION(:,:)   :: decomp_gp_ie   => NULL() ! bd
  REAL(dp), POINTER, DIMENSION(:,:)   :: decomp_gp_je   => NULL() ! bd
  REAL(dp), POINTER, DIMENSION(:,:)   :: decomp_gp_pe   => NULL() ! bd
  ! 2D FIELDS "RESCUED" FROM COSMO
  ! snow covered fraction
  REAL(dp), POINTER, DIMENSION(:,:)   :: cvs        => NULL()   
  ! water covered fraction
  REAL(dp), POINTER, DIMENSION(:,:)   :: cvw        => NULL()   
  ! leaf stomatal resistence
  REAL(dp), POINTER, DIMENSION(:,:)   :: rco_leaf   => NULL()   
  ! sea land mask sea=0 land=1
  REAL(dp), POINTER, DIMENSION(:,:)   :: slm => NULL()
  ! field capacity of soil
  REAL(dp), POINTER, DIMENSION(:,:)   :: wsmx => NULL()
  ! deep soil temperature
  REAL(dp), POINTER, DIMENSION(:,:,:)   :: tsoil    => NULL()
!!$  ! snow depth
!!$  REAL(dp), POINTER, DIMENSION(:,:)   :: sn => NULL()
  ! soil wetness
  REAL(dp), POINTER, DIMENSION(:,:)     :: ws => NULL()
  ! surface temperatur over land
  REAL(dp), POINTER, DIMENSION(:,:)   :: tslm1 => NULL()
  ! surface temperatur over water
  REAL(dp), POINTER, DIMENSION(:,:)   :: tsw => NULL()
  ! fraction of glacier covered land points
  REAL(dp), POINTER, DIMENSION(:,:)   :: glac => NULL()
  ! 10m wind speed
  REAL(dp), POINTER, DIMENSION(:,:)   :: wind10_2d => NULL()
!!$  ! 'cos(solar zenith angle)
!!$  REAL(dp), POINTER, DIMENSION(:,:)   :: cossza_2d => NULL()
  ! surface roughness
  REAL(dp), POINTER, DIMENSION(:,:)   :: az0 => NULL()
  ! soil moisture stress function
  REAL(dp), POINTER, DIMENSION(:,:)   :: fws => NULL()
  ! virtual temperature at surface
  REAL(dp), POINTER, DIMENSION(:,:)   :: tvir => NULL()
  ! surface virtual temperature(land)
  REAL(dp), POINTER, DIMENSION(:,:)   :: tvl => NULL()
  ! surface virtual temperature (water)
  REAL(dp), POINTER, DIMENSION(:,:)   :: tvw => NULL()
  ! surface virtual temperature (ice)
  REAL(dp), POINTER, DIMENSION(:,:)   :: tvi => NULL()
  ! neutral drag coeff., land
  REAL(dp), POINTER, DIMENSION(:,:)   :: cdnl => NULL()
  ! neutral drag coeff., water
  REAL(dp), POINTER, DIMENSION(:,:)   :: cdnw => NULL()
  ! neutral drag coeff., ice
  REAL(dp), POINTER, DIMENSION(:,:)   :: cdni => NULL()
  ! momentum drag coeff., land
  REAL(dp), POINTER, DIMENSION(:,:)   :: cfml => NULL()
  ! momentum drag coeff., water
  REAL(dp), POINTER, DIMENSION(:,:)   :: cfmw => NULL()
  ! momentum drag coeff., ice
  REAL(dp), POINTER, DIMENSION(:,:)   :: cfmi => NULL()
  ! exchange parameter, land
  REAL(dp), POINTER, DIMENSION(:,:)   :: cfncl => NULL()
  ! exchange parameter, water
  REAL(dp), POINTER, DIMENSION(:,:)   :: cfncw => NULL()
  ! exchange parameter, ice
  REAL(dp), POINTER, DIMENSION(:,:)   :: cfnci => NULL()
  ! Richardson number (land)
  REAL(dp), POINTER, DIMENSION(:,:)   :: ril => NULL()
  ! Richardson number (water)
  REAL(dp), POINTER, DIMENSION(:,:)   :: riw => NULL()
  ! Richardson number (ice)
  REAL(dp), POINTER, DIMENSION(:,:)   :: rii => NULL()
  ! Richardson number (ice)
  REAL(dp), POINTER, DIMENSION(:,:)   :: zust_2d => NULL()

  !### add POINTER

  ! 2D FIELDS NOT AVAILABLE IN COSMO
  ! sea ice fraction
  REAL(dp), POINTER, DIMENSION(:,:)   :: seaice => NULL()

  ! POINTER FOR GLOBAL FIELDS CONTAINING THE GEOGRAPHICAL LONG- AND LATITUDES
  REAL(dp), POINTER, DIMENSION(:,:)   :: lat_tot => NULL()
  REAL(dp), POINTER, DIMENSION(:,:)   :: lon_tot => NULL()
  
  ! POINTER FOR CONVECTION
  REAL  (dp), POINTER             ::           &
    massfu   (:,:,:)   , &  ! updraft mass flux          (kg/(m^2 s))
    u_detr   (:,:,:)   , &  ! updraft detrainment flux   (kg/(m^2 s))
    u_entr   (:,:,:)   , &  ! updraft entrainment flux   (kg/(m^2 s))
    massfd   (:,:,:)   , &  ! downdraft mass flux        (kg/(m^2 s))
    d_detr   (:,:,:)   , &  ! downdraft detrainment flux (kg/(m^2 s))
    d_entr   (:,:,:)   , &  ! downdraft entrainment flux (kg/(m^2 s))
    cv_precflx (:,:,:) , &  ! convective precipitation flux (3d) (kg/(m^2 s))
    cv_snowflx (:,:,:) , &  ! conv. snow precipitation flux (3d) (kg/(m^2 s))
    cv_precnew (:,:,:) , &  ! freshly formed conv. prec. flux (3d) (kg/(m^2 s))
    cv_snownew (:,:,:) , &  ! shly formed conv. snow flux(3d) (kg/(m^2 s))
    cv_lwc (:,:,:) ,     &  ! convective cloud water content(3d) (kg/kg)
    cv_iwc (:,:,:) ,     &  ! convective cloud ice content(3d)   (kg/kg)
    cv_rform(:,:,:),     &  ! convective precipitation formation (water) (kg/kg)
    cv_sform(:,:,:),     &  ! convective precipitation formation (snow)  (kg/kg)
    cu_top(:,:),         &  ! index of convective cloud top level (-)
    cu_bot(:,:),         &  ! index of convective cloud bottom level (-)
 !   cover_ls(:,:,:),     &  ! large scale cloud cover
    preccover_ls(:,:,:), &  ! large scale precipitation cloud cover
    precflx_ls(:,:,:),   &  ! large scale rain precipitation flux
    snowflx_ls(:,:,:),   &  ! large scale snow precipitation flux 
    precflxno_ls(:,:,:), &  ! large scale rain precipitation flux without
                            ! cloud production of new rain
    snowflxno_ls(:,:,:), &  ! large scale snow precipitation flux without 
                            ! cloud production of new snow
    rainform_bave(:,:,:),&  ! rain formation rate kg/kg, averaged over box
    snowform_bave(:,:,:),&  ! snow/ice formation rate kg/kg, averaged over box
    lwc_ls(:,:,:),       &  ! large scale liquid water content
    iwc_ls(:,:,:),       &  ! large scale snow/ice content
    precmelt_ls(:,:,:),  &  ! large scale frozen precipitation melting (kg/kg)
    sediice_ls(:,:,:)!,   &  ! large scale ice sedimentation (kg/kg)
  REAL(dp), POINTER, DIMENSION(:,:,:) :: aclc => NULL()

  !****************************************************************************
  ! POINTER FOR 1D FIELDS
  !****************************************************************************
  REAL(dp), DIMENSION(:), POINTER :: hyam => NULL()
  REAL(dp), DIMENSION(:), POINTER :: hyai => NULL()
  REAL(dp), DIMENSION(:), POINTER :: hybm => NULL()
  REAL(dp), DIMENSION(:), POINTER :: hybi => NULL()

  !****************************************************************************
  !****************************************************************************
  ! POINTER FOR REFERENCES
  !****************************************************************************
  !****************************************************************************
  ! 1D
  REAL(DP), DIMENSION(:),   POINTER :: gboxarea => NULL()
  ! 2D
  REAL(DP), DIMENSION(:,:), POINTER :: qte => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: tte => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: xlte => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: xite => NULL()

  ! 3D
  REAL(DP), DIMENSION(:,:,:), POINTER :: tm1 => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: qm1 => NULL()


  ! LOCALIZED PARAMETERS
  ! - GRID CONTROL
  INTEGER :: nlev             ! number vertical levels
  INTEGER :: nlevp1           ! nlev +1 (0 number interface levels)
  INTEGER :: nproma           ! length of x-direction in COSMO-grid (ie)
  INTEGER :: ngpblks          ! length of y-direction in COSMO-grid (je)
  INTEGER :: nlon, ngl
  INTEGER :: npromz           ! GP: longitud number of gridboxes on current PE
  INTEGER :: nglon, nglat     ! GP: NUMBER OF LON AND LAT (1 if lcolumn=T)
  INTEGER :: kepin, kezin     ! number of pressure and height levels for OUTPUT

  INTEGER :: kproma           ! defined in "local" loops (ie)
  INTEGER :: jrow             ! defined in "local" loops (je)

  ! imitate ECHAM5
  LOGICAL :: lmidatm = .FALSE.

  ! LOGICAL (needed for MMD COULPLING ECHAM5<->COSMO<->COSMO)
  LOGICAL :: L_IS_CLIENT = .FALSE.

  ! 4-dimensional fields ! DUMMY TO SIMPLIFY E5/C4 TREATMENT IN SUBMODELS
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: pxtems     => NULL()
  ! SUBROUTINES
  !
CONTAINS

  !***************************************************************************
  SUBROUTINE main_data_initialize

    USE messy_main_timer,    ONLY: timer_set_time_step_len

    IMPLICIT NONE

    nlev    = ke!_tot
    nlevp1  = ke+1
    nproma  = ie!_max
    ngpblks = je!_max
    ngl     = je_tot
    nlon    = ie_tot
    npromz  = nproma 
    nglon   = ie
    nglat   = je
    kepin   = root%kepin
    kezin   = root%kezin
    kproma  = nproma 

    ! set integration time step length
    CALL timer_set_time_step_len(l2tls)

  END SUBROUTINE main_data_initialize
  !----------------------------------------------------------------------------
  !***************************************************************************
  SUBROUTINE main_data_init_memory

    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_channel_object_reference    &
                                      , new_attribute, get_channel_object
    USE messy_main_channel_repr,  ONLY: get_representation_id
    USE messy_main_constants_mem, ONLY: r_earth => radius_earth, pi, g, DTR
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_init_memory'
    INTEGER :: status
    INTEGER :: repr_scalar ! TEST
    INTEGER :: reprid_2d  ! GP_2D_HORIZONTAL
    INTEGER :: repr_mid3d ! GP_3D_MID
    INTEGER :: repr_int3d ! GP_3D_INT
    INTEGER :: repr_soil1 ! GP_3D_SOIL1
    INTEGER :: repr_soil2 ! GP_3D_SOIL2
    INTEGER :: repr_lev
    INTEGER :: repr_ilev
    INTEGER :: i,j
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: p4 => NULL()

    ! create new channel
    CALL new_channel (status, modstr, lrestreq=.FALSE.)
    CALL channel_halt(substr, status)

    ! TEST
    CALL get_representation_id(status, 'SCALAR', repr_scalar)
    CALL channel_halt(substr, status)

    ! ########################################################################
    ! ---------------------------- GP_2D_HORIZONTAL --------------------------
    ! ########################################################################
    CALL get_representation_id(status, 'GP_2D_HORIZONTAL', reprid_2d)
    CALL channel_halt(substr, status)

    ! ------------------------------------------------------------------------
    ! create 2D_HORIZONTAL channel elements to output crlat 
    CALL new_channel_object(status, modstr,  'crlat_2d', &
         p2=crlat_2d, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'crlat_2d', 'long_name', c='crlat')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'crlat_2d', 'units', c=' ')
    CALL channel_halt(substr, status)

    DO i=1,ie 
       crlat_2d(i,:) = crlat(:,1)
    ENDDO
    ! ------------------------------------------------------------------------
    ! create 2D_HORIZONTAL channel elements to output crlat 
    CALL new_channel_object(status, modstr,  'crlati_2d', &
         p2=crlati_2d, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'crlati_2d' &
         , 'long_name', c='crlat interfaces')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'crlati_2d', 'units', c=' ')
    CALL channel_halt(substr, status)
    DO i=1,ie 
       crlati_2d(i,:) = crlat(:,2)
    ENDDO
    ! ------------------------------------------------------------------------
    ! create 2D_HORIZONTAL channel elements 2D grid box area
    CALL new_channel_object(status, modstr,  'gboxarea', &
         p2=gboxarea_2d, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gboxarea' &
         , 'long_name', c='gridbox area')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gboxarea', 'units', c='m2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    DO j = 1, je
       DO i = 1, ie
          gboxarea_2d(i,j) = crlat(j,1) * r_earth**2 * (pi/180.0_dp)**2 &
               * dlon * dlat
       ENDDO
    ENDDO
    ! ------------------------------------------------------------------------
    ! create 2D_HORIZONTAL channel elements 2D longitude
    CALL new_channel_object(status, modstr,  'philon_2d', &
         p2=philon_2d, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philon_2d' &
         , 'long_name', c='geographical longitude')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philon_2d', 'units', c='degree')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    DO j = 1, je
       DO i = 1, ie
          philon_2d(i,j) = rlon(i,j) * 180._dp/pi
       ENDDO
    ENDDO
    ! ------------------------------------------------------------------------
    ! create 2D_HORIZONTAL channel elements 2D longitude
    CALL new_channel_object(status, modstr,  'philat_2d', &
         p2=philat_2d, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philat_2d' &
         , 'long_name', c='geographical latitude')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philat_2d', 'units', c='degree')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    DO j = 1, je
       DO i = 1, ie
          philat_2d(i,j) = rlat(i,j) * 180._dp/pi
       ENDDO
    ENDDO
    ! op_pj_20110714+
    ! ------------------------------------------------------------------------
    ! COS(longitude)
    CALL new_channel_object(status, modstr,  'coslon', &
         p2=coslon_2d, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'coslon', &
         'long_name', c='cos(longitude)')
    CALL channel_halt(substr, status)
    coslon_2d(:,:) = COS(philon_2d(:,:)*DTR)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! SIN(longitude)
    CALL new_channel_object(status, modstr,  'sinlon', &
         p2=sinlon_2d, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sinlon', &
         'long_name', c='sin(longitude)')
    CALL channel_halt(substr, status)
    sinlon_2d(:,:) = SIN(philon_2d(:,:)*DTR)
    ! ------------------------------------------------------------------------
    ! COS(latitude)
    CALL new_channel_object(status, modstr,  'coslat', &
         p2=coslat_2d, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'coslat', &
         'long_name', c='cos(latitude)')
    CALL channel_halt(substr, status)
    coslat_2d(:,:) = COS(philat_2d(:,:)*DTR)
    ! ------------------------------------------------------------------------
    ! SIN(latitude)
    CALL new_channel_object(status, modstr,  'sinlat', &
         p2=sinlat_2d, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sinlat', &
         'long_name', c='sin(latitude)')
    CALL channel_halt(substr, status)
    sinlat_2d(:,:) = SIN(philat_2d(:,:)*DTR)
    ! ------------------------------------------------------------------------
    ! op_pj_20110714-

    ! (diagnostic; decomposition)
    CALL new_channel_object(status, modstr,  'decomp_gp_ie', &
         p2=decomp_gp_ie, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'decomp_gp_ie', &
         'long_name', c='x dimension index (1...ie)')
    CALL channel_halt(substr, status)
    DO i=1, ie
       decomp_gp_ie(i,:) = REAL(i,DP)
    END DO
    ! ------------------------------------------------------------------------
    ! (diagnostic; decomposition)
    CALL new_channel_object(status, modstr,  'decomp_gp_je', &
         p2=decomp_gp_je, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'decomp_gp_je', &
         'long_name', c='y dimension index (1...je)')
    CALL channel_halt(substr, status)
    DO j=1, je
       decomp_gp_je(:,j) = REAL(j,DP)
    END DO
    ! ------------------------------------------------------------------------
    ! (diagnostic; decomposition)
    CALL new_channel_object(status, modstr,  'decomp_gp_pe', &
         p2=decomp_gp_pe, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'decomp_gp_pe', &
         'long_name', c='processor number')
    CALL channel_halt(substr, status)
    decomp_gp_pe(:,:) = REAL(p_pe,DP)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! --- CHANNEL OBJECT REFERENCES ------------------------------------------
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------

    ! ########################################################################
    ! ----------------------------- GP_3D ------------------------------------
    ! ########################################################################
    CALL get_representation_id(status, 'GP_3D_MID', repr_mid3d)
    CALL channel_halt(substr, status)
    !
    CALL get_representation_id(status, 'GP_3D_INT', repr_int3d)
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! --- CHANNEL OBJECT REFERENCES ------------------------------------------
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'U_10M' &
         , modstr, 'u10')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u10', 'long_name', c='10m u-velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u10', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'V_10M' &
                                    , modstr, 'v10')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v10', 'long_name', c='10m v-velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v10', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    !
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'RAIN_GSP' &
                                    , modstr, 'prl')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prl' &
    , 'long_name', c='large-scale precipitation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prl', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'RAIN_CON' &
                                    , modstr, 'prc')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prc' &
    , 'long_name', c='convective precipitation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prc', 'units', c='m')
    CALL channel_halt(substr, status)
    !
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'FC' &
                                    , modstr, 'coriol')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'coriol' &
    , 'long_name', c='coriols parameter')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'coriol', 'units', c='s-1')
    CALL channel_halt(substr, status)    !
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'QVTENS' &
                                    , modstr, 'qte_3d')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qte_3d' &
    , 'long_name', c='specific humidity tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qte_3d', 'units', c='s-1')
    CALL channel_halt(substr, status)    !
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'TTENS' &
                                    , modstr, 'tte_3d')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tte_3d' &
    , 'long_name', c='temperature tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tte_3d', 'units', c='K s-1')
    CALL channel_halt(substr, status)    !
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'QITENS' &
                                    , modstr, 'xite_3d')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xite_3d' &
    , 'long_name', c='cloud ice tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xite_3d', 'units', c='s-1')
    CALL channel_halt(substr, status)    !
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'QCTENS' &
                                    , modstr, 'xlte_3d')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlte_3d' &
    , 'long_name', c='liquid water tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlte_3d', 'units', c='s-1')
    CALL channel_halt(substr, status)    !
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, &
         'COSMO_ORI', 'ALB_RAD', modstr, 'albedo')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'albedo', &
         'long_name', c='surface albedo')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'albedo', 'units', c='-')
    CALL channel_halt(substr, status)
    !
    ! ------------------------------------------------------------------------
    ! --- NEW CHANNEL OBJECTS ------------------------------------------------
    ! ------------------------------------------------------------------------
    ! (from main_data_global_start)
    !                   pointer
    ! pressure at middle of box ("full level pressure")
    CALL new_channel_object(status, modstr,  'press', &
         p3=press_3d, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press', 'long_name', c='pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! pressure at level interfaces ("half level pressure")
    CALL new_channel_object(status, modstr,  'pressi', &
         p3=pressi_3d, reprid=repr_int3d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressi', &
         'long_name', c='interface pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressi', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! geopotential at middle of box ("full level geopotential")
    CALL new_channel_object(status, modstr,  'geopot', &
         p3=geopot_3d, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopot' &
         , 'long_name', c='geopotential')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopot', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! geopotential height at interface below current level
    CALL new_channel_object(status, modstr,  'geopoti', &
         p3=geopoti_3d, reprid=repr_int3d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopoti', &
         'long_name', c='interface geopotential')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'geopoti', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    !  grid mass in kg
    CALL new_channel_object(status, modstr,  'grmass', &
         p3=grmass, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grmass', 'long_name', c='mass of gridbox')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grmass', 'units', c='kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! op_pj_20100713+
    !  mass of dry air in kg
    CALL new_channel_object(status, modstr,  'grmassdry', &
         p3=grmassdry, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grmassdry', 'long_name' &
         , c='mass of dry air')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grmassdry', 'units', c='kg')
    CALL channel_halt(substr, status)
    ! op_pj_20100713-
    ! ------------------------------------------------------------------------
    !  volume of grid box in m**3
    CALL new_channel_object(status, modstr,  'grvol', &
         p3=grvol, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grvol', 'long_name', c='volume of gridbox')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'grvol', 'units', c='m**3')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from vdiff.f90); local variable converted to pointer
    CALL new_channel_object(status, modstr,  'tpot', &
         p3=tpot_3d, reprid=repr_mid3d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tpot', &
         'long_name', c='potential temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tpot', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! vorticity
    CALL new_channel_object(status, modstr,'vom1', p3=vom1, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vom1'  &
         , 'long_name', c='vorticity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vom1', 'units', c='s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! DEFINE CHANNEL OBJECTS FOR PROGNOSTIC VARIABLES
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! 3 D
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! Temperature
    ! Temperature
    CALL new_channel_object(status, modstr,'tm1', p3=tm1_3d, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1'  &
         , 'long_name', c='temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#ifdef I2CINC
    ! Temperature BOUNDARIES
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'T_BD' &
         , modstr, 'tm1_BD')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1_BD' &
         , 'long_name', c='temperature boundaries')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1_BD', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#endif
    ! Temperature tendency
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'TTENS' &
         , modstr,'tte')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tte'  &
         , 'long_name', c='temperature tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tte', 'units', c='Ks-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    !pressure deviation
    CALL new_channel_object(status, modstr,'ppm1', p3=ppm1_3d, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ppm1'  &
         , 'long_name', c='pressure deviation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ppm1', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#ifdef I2CINC
    ! pressure variation BOUNDARIES
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'PP_BD' &
         , modstr, 'ppm1_BD')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ppm1_BD' &
         , 'long_name', c='pressure deviation boundaries')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ppm1_BD', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#endif
    ! ------------------------------------------------------------------------
    ! horizontal wind velocity
    CALL new_channel_object(status, modstr,'um1', p3=um1, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'um1'  &
         , 'long_name', c='horizontal wind velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'um1', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#ifdef I2CINC
    ! horizintal wind velocity BOUNDARIES
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'U_BD' &
         , modstr, 'um1_BD')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'um1_BD' &
         , 'long_name', c='horizontal wind boundaries')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'um1_BD', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#endif
    ! ------------------------------------------------------------------------
    ! horizontal wind velocity
    CALL new_channel_object(status, modstr,'vm1', p3=vm1, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vm1'  &
         , 'long_name', c='horizontal wind velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vm1', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#ifdef I2CINC
    ! horizintal wind velocity BOUNDARIES
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'V_BD' &
         , modstr, 'vm1_BD')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vm1_BD' &
         , 'long_name', c='horizontal wind boundaries')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vm1_BD', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#endif
    ! vertical wind velocity
    CALL new_channel_object(status, modstr,'wm1', p3=wm1, reprid=repr_int3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wm1'  &
         , 'long_name', c='vertical wind velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wm1', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ----------------------------------------------------------------------
#ifdef I2CINC
    IF ( .NOT. lw_freeslip ) THEN
       ! horizintal wind velocity BOUNDARIES
       CALL new_channel_object_reference(status, 'COSMO_ORI', 'W_BD' &
            , modstr, 'wm1_BD')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'wm1_BD' &
            , 'long_name', c='vertical wind boundaries')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'wm1_BD', 'units', c='m s-1')
       CALL channel_halt(substr, status)
       ! ----------------------------------------------------------------------
    ENDIF
#endif
    ! ------------------------------------------------------------------------
    ! water vapor
    CALL new_channel_object(status, modstr,'qv', p3=qm1_3d, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qv'  &
         , 'long_name', c='water vapor')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qv', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#ifdef I2CINC
    ! water vapor BOUNDARIES
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'QV_BD' &
         , modstr, 'qv_BD')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qv_BD' &
         , 'long_name', c='specific water vapor content')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qv_BD', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#endif
    ! ------------------------------------------------------------------------
    ! water vapor tendency
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'QVTENS' &
         , modstr,'qte')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qte'  &
         , 'long_name', c='specific humidity tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qte', 'units', c='kg/(kg s)')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    !----------------------------------------------------------------
    ! cloud water
    CALL new_channel_object(status, modstr,'qc', p3=qc_3d, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qc'  &
         , 'long_name', c='specific cloud water content')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qc', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#ifdef I2CINC
    ! cloud water content BOUNDARIES
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'QC_BD' &
         , modstr, 'qc_BD')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qc_BD' &
            , 'long_name', c='cloud water boundaries')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qc_BD', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#endif
    IF (lprog_qi) THEN
    ! ------------------------------------------------------------------------
    ! cloud Ice content
       CALL new_channel_object(status, modstr,'qi', p3=qi_3d, reprid=repr_mid3d)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qi'  &
            , 'long_name', c='specific cloud ice content')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qi', 'units', c='kg/kg')
       CALL channel_halt(substr, status)
       ! ----------------------------------------------------------------------
    ENDIF
#ifdef I2CINC
    IF (llb_qi) THEN
       ! Ice water  BOUNDARIES
       CALL new_channel_object_reference(status, 'COSMO_ORI', 'QI_BD' &
            , modstr, 'qi_BD')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qi_BD' &
            , 'long_name', c='cloud ice boundaries')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qi_BD', 'units', c='kg/kg')
       CALL channel_halt(substr, status)
       ! ---------------------------------------------------------------------
    ENDIF
#endif
    IF (lprogprec) THEN
       ! ---------------------------------------------------------------------
       ! specific rain content
       CALL new_channel_object(status, modstr,'qr', p3=qr_3d, reprid=repr_mid3d)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qr'  &
            , 'long_name', c='specific rain content')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qr', 'units', c='kg/kg')
       CALL channel_halt(substr, status)
       ! ----------------------------------------------------------------------
    ENDIF
#ifdef I2CINC
    IF (llb_qr_qs) THEN
       ! spec. rain content  BOUNDARIES
       CALL new_channel_object_reference(status, 'COSMO_ORI', 'QR_BD' &
            , modstr, 'qr_BD')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qr_BD' &
            , 'long_name', c='spec. rain content boundaries')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qr_BD', 'units', c='kg/kg')
       CALL channel_halt(substr, status)
       ! ---------------------------------------------------------------------
    ENDIF
#endif
    IF (lprogprec .AND. itype_gscp > 1) THEN
       ! -------------------------------------------------------------------
       ! specific snow content
       CALL new_channel_object(status, modstr,'qs' &
               , p3=qs_3d, reprid=repr_mid3d)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qs'  &
            , 'long_name', c='specific snow content')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qs', 'units', c='kg/kg')
       CALL channel_halt(substr, status)
       ! -------------------------------------------------------------------
    ENDIF
#ifdef I2CINC
    IF (llb_qr_qs) THEN
       ! spec. snow content  BOUNDARIES
       CALL new_channel_object_reference(status, 'COSMO_ORI', 'QS_BD' &
            , modstr, 'qs_BD')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qs_BD' &
            , 'long_name', c='spec. snow content boundaries')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qs_BD', 'units', c='kg/kg')
       CALL channel_halt(substr, status)
       ! ------------------------------------------------------------------
    ENDIF
#endif
    IF (itype_gscp ==4 .AND. lprogprec) THEN
       ! -------------------------------------------------------------------
       ! specific graupel content
       CALL new_channel_object(status, modstr,'qg' &
            , p3=qg_3d, reprid=repr_mid3d)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qg'  &
               , 'long_name', c='specific graupel content')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qg', 'units', c='kg/kg')
       CALL channel_halt(substr, status)
       ! -------------------------------------------------------------------
    ENDIF 
#ifdef I2CINC
    IF (llb_qg) THEN
       ! spec. snow content  BOUNDARIES
       CALL new_channel_object_reference(status, 'COSMO_ORI', 'QG_BD' &
            , modstr, 'qg_BD')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qg_BD' &
            , 'long_name', c='spec. graupel content boundaries')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'qg_BD', 'units', c='kg/kg')
       CALL channel_halt(substr, status)
       ! ------------------------------------------------------------------
    ENDIF 
#endif
    IF (lmulti_layer) THEN
       CALL get_representation_id(status, 'GP_3D_SOIL2', repr_soil2)
       CALL channel_halt(substr, status)
       ! -------------------------------------------------------------------
       ! multi-layer soil temperature
       CALL new_channel_object(status, modstr,'t_so' &
            , p3=t_so_3d, reprid=repr_soil2)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 't_so'  &
            , 'long_name', c='multi-layer soil temperature')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 't_so', 'units', c='K')
       CALL channel_halt(substr, status)
       ! -------------------------------------------------------------------
       ! multi-layer soil temperature
       CALL new_channel_object(status, modstr,'t_so_ini' &
            , p3=t_so_ini, reprid=repr_soil2)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 't_so_ini'  &
            , 'long_name', c='multi-layer soil temperature')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 't_so_ini', 'units', c='K')
       CALL channel_halt(substr, status)
       ! -------------------------------------------------------------------
       ! multi-layer soil moisture
       CALL get_representation_id(status, 'GP_3D_SOIL1', repr_soil1)
       CALL channel_halt(substr, status)
       
       CALL new_channel_object(status, modstr,'w_so' &
            , p3=w_so_3d, reprid=repr_soil1)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'w_so'  &
         , 'long_name', c='multi-layer soil moisture')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'w_so', 'units', c='m H2O')
       CALL channel_halt(substr, status)
       ! -------------------------------------------------------------------
    ENDIF
   ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'xim1' &
         , p3=xim1_3d, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xim1', &
         'long_name', c='cloud ice (xim1)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xim1', 'units', c='kg kg-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'xlm1' &
         , p3=xlm1_3d, reprid=repr_mid3d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlm1', &
         'long_name', c='cloud water (xlm1)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlm1','units', c='kg kg-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! 2 D
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    ! surface pressure
    CALL new_channel_object(status, modstr,'ps', p2=aps, reprid=reprid_2D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ps'  &
         , 'long_name', c='surface pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ps', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! temperature of the ground surface (soil) 
    CALL new_channel_object(status, modstr,'tsurf' &
         , p2=tsurf_2d, reprid=reprid_2D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsurf'  &
         , 'long_name', c='temperature of the ground surface (soil) ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsurf', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! temperature of the ground surface (water) 
    CALL new_channel_object(status, modstr,'tsw', p2=tsw, reprid=reprid_2D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsw'  &
         , 'long_name', c='temperature of the ground surface (water) ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsw', 'units', c='K')
    CALL channel_halt(substr, status)
     ! ------------------------------------------------------------------------
    ! temperature of the ground surface (soil) 
    CALL new_channel_object(status, modstr,'tslm1', p2=tslm1, reprid=reprid_2D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tslm1'  &
         , 'long_name', c='temperature of the ground surface (soil) ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tslm1', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr, 'wsmx', p2=wsmx, reprid=reprid_2D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wsmx', &
         'long_name', c='field capacity of soil')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wsmx', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
   ! ------------------------------------------------------------------------
#ifdef I2CINC
    ! temperature of the ground surface (soil)  BOUNDARIES
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'T_S_BD' &
         , modstr, 't_s_BD')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 't_s_BD' &
         , 'long_name', c='boundaries of t_s')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 't_s_BD', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#endif
    ! temperature of the snow-surface
    CALL new_channel_object(status, modstr,'t_snow' &
         , p2=t_snow_2d, reprid=reprid_2D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 't_snow'  &
         , 'long_name', c='temperature of the snow-surface')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 't_snow', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#ifdef I2CINC
    ! temperature of the snow-surface  BOUNDARIES
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'T_SNOW_BD' &
         , modstr, 't_snow_BD')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 't_snow_BD' &
         , 'long_name', c='boundaries of t_snow')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 't_snow_BD', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#endif
    ! water content of snow
    CALL new_channel_object(status, modstr,'w_snow'&
         , p2=w_snow_2d, reprid=reprid_2D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w_snow'  &
         , 'long_name', c='water content of snow')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w_snow', 'units', c='m H2O')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#ifdef I2CINC
    ! water content of snow BOUNDARIES
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'W_SNOW_BD' &
         , modstr, 'w_snow_BD')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w_snow_BD' &
         , 'long_name', c='boundaries of w_snow')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w_snow_BD', 'units', c='m H2O')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#endif
    ! water content of interception water 
    CALL new_channel_object(status, modstr,'w_i', p2=w_i_2d, reprid=reprid_2D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w_i'  &
         , 'long_name', c='water content of interception water ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w_i', 'units', c='m H2O')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! specific water vapor content at the surface
    CALL new_channel_object(status, modstr,'qv_s'&
         , p2=qv_s_2d, reprid=reprid_2D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qv_s'  &
         , 'long_name', c='specific water vapor content at the surface')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qv_s', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#ifdef I2CINC
    ! water content of snow BOUNDARIES
    CALL new_channel_object_reference(status,'COSMO_ORI', 'QV_S_BD' &
         , modstr, 'qv_s_BD')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qv_s_BD' &
         , 'long_name', c='boundaries of qv_s')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qv_s_BD', 'units', c='kg/kg')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
#endif
    ! DEFINED FOR CONVECTIVE TRANSPORT AND SCAVENGING

       ! CONVECTIVE FLUXES NEEDED BY CVTRANS 
       IF (itype_conv == 0) THEN ! sofar only implemented for TIEDKE SCHEME 
          ! updraft mass flux
          CALL new_channel_object(status, modstr, 'massfu' &
               , p3=massfu , reprid=repr_mid3d, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'massfu', &
               'long_name', c='updraft mass flux')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'massfu' &
               , 'units', c='kg /(m^2 s)')
          CALL channel_halt(substr, status)
          ! downward massflux
          CALL new_channel_object(status, modstr, 'massfd' &
               , p3=massfd , reprid=repr_mid3d, lrestreq=.TRUE.  )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'massfd', &
               'long_name', c='downward mass flux')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'massfd' &
               , 'units', c='kg /(m^2 s)')
          CALL channel_halt(substr, status)
          ! upward entraining mass flux
          CALL new_channel_object(status, modstr, 'u_entr'&
               , p3=u_entr , reprid=repr_mid3d, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'u_entr', &
               'long_name', c='upward entraining mass flux')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'u_entr' &
               , 'units', c='kg /(m^2 s)')
          CALL channel_halt(substr, status)
          ! downward entraining mass flux
          CALL new_channel_object(status, modstr, 'd_entr' &
               , p3=d_entr , reprid=repr_mid3d , lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'd_entr', &
               'long_name', c='downward entraining mass flux')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'd_entr' &
               , 'units', c='kg /(m^2 s)')
          CALL channel_halt(substr, status)
          ! upward detraining mass flux
          CALL new_channel_object(status, modstr, 'u_detr' &
               , p3=u_detr , reprid=repr_mid3d, lrestreq=.TRUE.  )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'u_detr', &
               'long_name', c='upward detraining mass flux')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'u_detr' &
               , 'units', c='kg /(m^2 s)')
          CALL channel_halt(substr, status)
          ! downward detraining mass flux
          CALL new_channel_object(status, modstr, 'd_detr' &
               , p3=d_detr  , reprid=repr_mid3d, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'd_detr', &
               'long_name', c='downward detraining mass flux')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'd_detr' &
               , 'units', c='kg /(m^2 s)')
          CALL channel_halt(substr, status)
          ! convective precipitation flux (3d)
          CALL new_channel_object(status, modstr, 'cv_precflx' &
               , p3=cv_precflx  , reprid=repr_mid3d, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_precflx', &
               'long_name', c='convective precipitation flux (3d)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_precflx' &
               , 'units', c='kg /(m^2 s)')
          CALL channel_halt(substr, status)
          ! convective snow precipitation flux (3d)
          CALL new_channel_object(status, modstr, 'cv_snowflx' &
               , p3=cv_snowflx  , reprid=repr_mid3d, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_snowflx', &
               'long_name', c='convective snow precipitation flux (3d)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_snowflx' &
               , 'units', c='kg /(m^2 s)')
          CALL channel_halt(substr, status)
          ! freshly formed convective precipitation flux (3d)
          CALL new_channel_object(status, modstr, 'cv_precnew' &
               , p3=cv_precnew  , reprid=repr_mid3d, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_precnew', &
               'long_name', c='convective precipitation flux (3d)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_precnew' &
               , 'units', c='kg /(m^2 s)')
          CALL channel_halt(substr, status)

          ! freshly formed convective snow precipitation flux (3d)
          CALL new_channel_object(status, modstr, 'cv_snownew' &
               , p3=cv_snownew  , reprid=repr_mid3d, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_snownew', &
               'long_name', c='convective snow precipitation flux (3d)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_snownew' &
               , 'units', c='kg /(m^2 s)')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'cv_lwc' &
               , p3=cv_lwc, reprid=repr_mid3d, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_lwc', &
               'long_name', c='convective cloud water content(3d)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_lwc', 'units', c='kg/kg')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'cv_iwc' &
               , p3=cv_iwc, reprid=repr_mid3d, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_iwc', &
               'long_name', c='convective cloud ice content(3d)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_iwc', 'units', c='kg/kg')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'cv_rform' &
               , p3=cv_rform, reprid=repr_mid3d, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_rform', &
               'long_name', c='convective precipitation formation (water)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_rform', 'units', c='kg/kg')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'cv_sform' &
               , p3=cv_sform, reprid=repr_mid3d, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_sform', &
               'long_name', c='convective precipitation formation (snow)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cv_sform', 'units', c='kg/kg')
          CALL channel_halt(substr, status)

          ! Define here  top and base level of convection
          ! for the COSMO inherent objects it depends on the namelist
          ! whether these are instantaneous values
          CALL new_channel_object(status, modstr, 'cu_top' &
               , p2=cu_top, reprid=reprid_2D, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cu_top', &
               'long_name', c='level index of convective cloud top')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cu_top', 'units', c='-')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'cu_bot' &
               , p2=cu_bot, reprid=reprid_2D, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cu_bot', &
               'long_name', c='level index of convective cloud bottom')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cu_bot', 'units', c='-')
          CALL channel_halt(substr, status)
       ELSE
          CALL warning_bi(substr, &
       'SCAV, LNOX and CVTRANS are only working with TIEDTKE CONVECTION SCHEME SO FAR')
       ENDIF
!!$       ! large scale cloud cover
!!$       CALL new_channel_object(status, modstr, 'cover' &
!!$            , p3=cover_ls, reprid=repr_mid3d, lrestreq=.TRUE. )
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr, 'cover', &
!!$            'long_name', c='large scale cloud cover')
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr, 'cover' &
!!$            , 'units', c='-')
!!$       CALL channel_halt(substr, status)

       ! large scale precipitation cloud cover
       CALL new_channel_object(status, modstr, 'prec_cover' &
            , p3=preccover_ls , reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'prec_cover', &
            'long_name', c='large scale precipitation cloud cover')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'prec_cover' &
            , 'units', c='-')
       CALL channel_halt(substr, status)

       ! large scale rain precipitation flux without
       ! cloud production of new rain
       CALL new_channel_object(status, modstr, 'precflx_no' &
            , p3=precflxno_ls , reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'precflx_no', &
            'long_name', c='large scale rain precipitation flux without'//&
         &' cloud production of new rain')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'precflx_no' &
            , 'units', c='kg /(m^2 s)')
       CALL channel_halt(substr, status)

       ! large scale snow precipitation flux without 
       ! cloud flxuction of new snow
       CALL new_channel_object(status, modstr, 'snowflx_no' &
            , p3=snowflxno_ls , reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'snowflx_no', &
            'long_name', c='large scale snow precipitation flux without'//&
         &' cloud production of new snow')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'snowflx_no' &
            , 'units', c='kg /(m^2 s)')
       CALL channel_halt(substr, status)
       !
      ! large scale rain precipitation flux without
       ! cloud production of new rain
       CALL new_channel_object(status, modstr, 'precflx' &
            , p3=precflx_ls , reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'precflx', &
            'long_name', c='large scale rain precipitation flux')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'precflx' &
            , 'units', c='kg /(m^2 s)')
       CALL channel_halt(substr, status)

       ! large scale snow precipitation flux
       CALL new_channel_object(status, modstr, 'snowflx' &
            , p3=snowflx_ls , reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'snowflx', &
            'long_name', c='large scale snow precipitation flux')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'snowflx' &
            , 'units', c='kg /(m^2 s)')
       CALL channel_halt(substr, status)
       !
       CALL new_channel_object(status, modstr, 'rain_form' &
            , p3=rainform_bave , reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'rain_form', &
            'long_name', c='large scale rain formation inside cloud')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'rain_form' &
            , 'units', c='kg/kg')
       CALL channel_halt(substr, status)
       !
       CALL new_channel_object(status, modstr, 'snow_form' &
            , p3=snowform_bave , reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'snow_form', &
            'long_name', c='large scale snow formation inside cloud')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'snow_form' &
            , 'units', c='kg/kg')
       CALL channel_halt(substr, status)
       !large scale cloud liquid water content
       CALL new_channel_object(status, modstr, 'lwc' &
            , p3=lwc_ls , reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'lwc', &
            'long_name', c='large scale cloud liquid water content')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'lwc' &
            , 'units', c='kg / kg')
       CALL channel_halt(substr, status)
       ! large scale cloud snow/ice content
       CALL new_channel_object(status, modstr, 'iwc' &
            , p3=iwc_ls, reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'iwc', &
            'long_name', c='large scale cloud snow/ice content')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'iwc' &
            , 'units', c='kg /kg')
       CALL channel_halt(substr, status)
       !large scale frozen precipitation melting
       CALL new_channel_object(status, modstr, 'prec_melt' &
            , p3=precmelt_ls , reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'prec_melt', &
            'long_name', c='large scale frozen precipitation melting')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'prec_melt' &
            , 'units', c='kg /(m^2 s)')
       CALL channel_halt(substr, status)
       !
       CALL new_channel_object(status, modstr, 'sedi_ice' &
            , p3=sediice_ls , reprid=repr_mid3d, lrestreq=.TRUE. )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'sedi_ice', &
            'long_name', c='large scale ice sedimentation')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'sedi_ice' &
            , 'units', c='kg /kg')
       CALL channel_halt(substr, status)
       ! ---------------------------------------------------------------------
       ! 
       CALL new_channel_object(status, modstr, 'rhum', &
            p3=rhum_3d, reprid=repr_mid3d, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'rhum', &
            'long_name', c='relative humidity')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'rhum', 'units', c='%')
       CALL channel_halt(substr, status)
       ! ----------------------------------------------------------------------
       CALL new_channel_object(status, modstr, 'aclc' &
            , p3=aclc, reprid=repr_mid3d)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'aclc', &
            'long_name', c='large scale cloud cover')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'aclc', 'units', c='-')
       CALL channel_halt(substr, status)
       ! ----------------------------------------------------------------------
!*************************************************************************
!*************************************************************************
! ADD CHANNEL REFERENCES OR CHANNEL OBJECTS FOR CONFORMATY WITH ECHAM HERE
!*************************************************************************
!*************************************************************************
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'FR_LAND' &
         , modstr, 'slf')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slf' &
         , 'long_name', c='sea land fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slf', 'units', c='0-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object_reference(status, 'COSMO_ORI', 'SOBS_RAD' &
         , modstr,  'srfl')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'srfl', &
         'long_name', c='net surface radiative flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'srfl', 'units', c='W m-2')
    CALL channel_halt(substr, status)
   ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr, 'cvs' &
         , p2=cvs, reprid=reprid_2D, lrestreq=.TRUE.) ! for dradon
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvs' &
         , 'long_name', c='snow covered fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvs', 'units', c='0-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr, 'cvw' &
         , p2=cvw, reprid=reprid_2D, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvw' &
         , 'long_name', c='water covered fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvw', 'units', c='0-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! stomatal resistance
    CALL new_channel_object(status, modstr,  'rco_leaf', &
         p2=rco_leaf, reprid=reprid_2d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rco_leaf', &
         'long_name', c='leaf stomatal resistance')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rco_leaf', 'units', c='s m-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
!!$    ! (from src_radiation.f90); local variable copied to pointer 
!!$    CALL new_channel_object(status, modstr,  'cossza', &
!!$         p2=cossza_2d, reprid=reprid_2d)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr, 'cossza', &
!!$         'long_name', c='cos(solar zenith angle)')
!!$    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! (from near_surface.f90); 
    CALL new_channel_object(status, modstr,  'wind10', &
         p2=wind10_2d, reprid=reprid_2d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10', &
         'long_name', c='10 m wind speed')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr, 'az0', p2=az0, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0', &
         'long_name', c='roughness length')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'az0', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! 0 = liquid ocean, 1 = ice
    ! This value is set to 0 over land, although it is actually undefined there
    CALL new_channel_object(status, modstr, 'seaice' &
         , p2=seaice, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'seaice', &
         'long_name', c='seaice fraction rel to ocean')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'seaice', 'units', c='-')
    CALL channel_halt(substr, status)
    ! INIT
    seaice(:,:) = 0._dp ! NO sea ice fraction currently available in COSMO
    ! ------------------------------------------------------------------------
    ! 0 = sea, 1 = land
    CALL new_channel_object(status, modstr, 'slm' &
         , p2=slm, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slm', 'long_name', c='land mask')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slm', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! fraction of land covered by glaciers 
    CALL new_channel_object(status, modstr, 'glac' &
         , p2=glac, reprid=reprid_2d )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'glac', 'long_name', &
         c='fraction of land covered by glaciers')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'glac', 'units', c='-')
    CALL channel_halt(substr, status)
    ! set to zero, as it is unknown in COSMO
    glac =0._dp
    ! ------------------------------------------------------------------------
    ! deep soil temperature
    IF (lmulti_layer) THEN
       CALL get_channel_object(status, 'COSMO', 't_so', p4=p4)
       CALL channel_halt(substr, status)
       tsoil => p4(:,:,:,1) 
    ELSE
       CALL get_channel_object(status, 'COSMO_ORI', 'T_CL', p4=p4)
       CALL channel_halt(substr, status)
       tsoil    => p4(:,:,1:1,1) 
    ENDIF

    ! ------------------------------------------------------------------------
    ! soil wetness
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status,  modstr, 'ws', p2= ws ,reprid=reprid_2D )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ws' &
         , 'long_name', c='frational water content (soil wetness)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ws', 'units', c='m')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in global start
    CALL new_channel_object(status, modstr,  'fws', p2=fws, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fws', &
         'long_name', c='soil moisture stress function')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fws', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in global start
    CALL new_channel_object(status, modstr,  'tvir', p2=tvir, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvir', &
         'long_name', c='surface virtual temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvir', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr,  'tvl', p2=tvl, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvl', &
         'long_name', c='surface virtual temperature(land)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvl', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr,  'tvw', p2=tvw, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvw', &
         'long_name', c='surface virtual temperature (water)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvw', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr,  'tvi', p2=tvi, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvi', &
         'long_name', c='surface virtual temperature (ice)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvi', 'units', c='K')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr,  'cdnl', p2=cdnl, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnl', &
         'long_name', c='neutral drag coeff., land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnl', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr,  'cdnw', p2=cdnw, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnw', &
         'long_name', c='neutral drag coeff., water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnw', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr,  'cdni', p2=cdni, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdni', &
         'long_name', c='neutral drag coeff., ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdni', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr,  'cfml', p2=cfml, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfml', &
         'long_name', c='momentum drag coeff., land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfml', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr,  'cfmw', p2=cfmw, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmw', &
         'long_name', c='momentum drag coeff., water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmw', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr, 'cfmi', p2=cfmi, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmi', &
         'long_name', c='momentum drag coeff., ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmi', 'units', c='-')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr,  'cfncl', p2=cfncl, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncl', &
         'long_name', c='exchange parameter, land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncl', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr, 'cfncw', p2=cfncw, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncw', &
         'long_name', c='exchange parameter, water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncw', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
     ! calculated in data_vdiff
    CALL new_channel_object(status, modstr, 'cfnci', p2=cfnci, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfnci', &
         'long_name', c='exchange parameter, ice ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfnci', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr,  'ril', p2=ril, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ril', &
         'long_name', c='Richardson number (land)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ril', 'units', c='1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
     ! calculated in data_vdiff
    CALL new_channel_object(status, modstr, 'riw', p2=riw, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'riw', &
         'long_name', c='Richardson number (water)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'riw', 'units', c='1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! calculated in data_vdiff
    CALL new_channel_object(status, modstr, 'rii', p2=rii, reprid=reprid_2d)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rii', &
         'long_name', c='Richardson number (ice)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rii', 'units', c='1')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'zust' &
         , p2=zust_2d, reprid=reprid_2d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zust', &
         'long_name', c='surface friction velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zust', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    zust_2d = 1.e-15
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------

    ! um_ak_20110627+
    ! 1-dimensional vertical fields
    CALL get_representation_id(status, 'GP_1D_LEV', repr_lev)
    CALL channel_halt(substr, status)

    CALL get_representation_id(status, 'GP_1D_ILEV', repr_ilev)
    CALL channel_halt(substr, status)

    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'dhyam', &
         p1=hyam, reprid=repr_lev)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhyam', 'long_name' &
         , c='hybrid A coefficient at layer midpoints')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhyam', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'dhybm', &
         p1=hybm, reprid=repr_lev)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhybm', 'long_name' &
            , c='hybrid B coefficient at layer midpoints')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhybm', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'dhyai', &
         p1=hyai, reprid=repr_ilev)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhyai', 'long_name' &
            , c='hybrid B coefficient at interface layer')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhyai', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    CALL new_channel_object(status, modstr,  'dhybi', &
         p1=hybi, reprid=repr_ilev)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhybi', 'long_name' &
            , c='hybrid B coefficient at interface layer')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dhybi', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! um_ak_20110627-

  END SUBROUTINE main_data_init_memory
  !-----------------------------------------------------------------------------
  SUBROUTINE main_data_init_coupling

    USE messy_main_blather_bi,    ONLY: error_bi

    IMPLICIT NONE

    INTEGER :: k
    CHARACTER(LEN=*), PARAMETER :: substr='main_data_init_memory'

    IF (ivctype /= 1) THEN
       write (0,*) 'ERROR  HYBRID COEFFICIENTS: ivctype==',ivctype
       write (0,*) 'ivctype is currently required to be 1 for COSMO/MESSy'

       CALL error_bi &
          ('SORRY! , COSMO/MESSy is currently only applicable with ivctype==1'&
            , substr)
    ENDIF
    DO k = 1, ke+1
      IF( vcoord(k) <= vcflat ) THEN
         hyai(k) = vcoord(k)*p0sl
         hybi(k) = 0.0_dp
      ELSE
         hyai(k) = vcflat*p0sl*(1.0_dp - vcoord(k))/(1.0_dp - vcflat)
         hybi(k) = (vcoord(k) - vcflat)/(1.0_dp - vcflat)
      ENDIF
   ENDDO
    DO k = 1, ke
       hyam(k) = (hyai(k)+hyai(k+1))/2._dp
       hybm(k) = (hybi(k)+hybi(k+1))/2._dp
    END DO
    ! ------------------------------------------------------------------------
  END SUBROUTINE main_data_init_coupling
  !***************************************************************************

  !***************************************************************************
  SUBROUTINE main_data_global_start

    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_constants_mem, ONLY: r_earth => radius_earth, pi, g &
                                      , c_vKar, vtmpc1
    USE messy_main_timer,         ONLY: timer_set_time_step_len, lstart &
                                      , YearDay
    USE messy_main_tools,         ONLY: tlucua, jptlucu1, jptlucu2
    USE data_soil,                ONLY: cfcap

    ! LOCAL
    CHARACTER (len=*), PARAMETER :: substr='main_data_global_start'
    REAL(DP),          PARAMETER :: dzeta = 1.0_dp
    REAL(DP),          PARAMETER :: plmin = 0.35_dp
    REAL(DP),          PARAMETER :: plmax = 0.75_dp
    REAL(DP),          PARAMETER :: pldiff= plmax-plmin
    INTEGER                      :: k,i,j 
!    INTEGER                      :: ntime ! index for time level
!                                          ! l2tls = T => ntime = nnow
!                                          ! l2tls = F => ntime = nnold
    INTEGER  :: idx(ie,je)
    INTEGER  :: it
    REAL(dp) :: pws(ie,je)
    REAL(dp) :: dx, dy
    REAL(dp) :: zq(ie,je,ke),zqs(ie,je,ke)
    LOGICAL  :: lookupoverflow = .FALSE.

    ! set integration time step length
    CALL timer_set_time_step_len(l2tls)

    IF (l2tls) THEN
       ntime = nnow
    ELSE
       ntime = nold
    ENDIF
    ! set prognostic variables:
    ! -- 3D 
    tm1    => T(1:ie,1:je,1:ke,ntime)
    tm1_3d(1:nproma,:,:)    = T(1:ie,:,1:ke,ntime)
!    tte_3d => ttens(1:ie,1:je,1:ke)
    qm1    => QV(1:ie,1:je,1:ke,ntime)
!    qte_3d => qvtens(1:ie,1:je,1:ke)

    um1(1:ie,1:je,1:ke)  = U(1:ie,1:je,1:ke,ntime)
    vm1(1:ie,1:je,1:ke)  = V(1:ie,1:je,1:ke,ntime)
    wm1(1:ie,1:je,1:ke+1) = W(1:ie,1:je,1:ke+1,ntime)
    ppm1_3d(1:nproma,:,:) = pp(1:ie,1:je,1:ke,ntime)

    qm1_3d(1:nproma,:,:)    = QV(1:ie,:,1:ke,ntime)
    IF (lprog_qi) THEN
       xim1_3d(1:nproma,:,:)   = QI(1:ie,:,1:ke,ntime)
    ELSE 
       xim1_3d(1:nproma,:,:)   =  0._dp
    ENDIF
    xlm1_3d(1:nproma,:,:)   = QC(1:ie,:,1:ke,ntime)
    IF (lprogprec) THEN
       xlm1_3d(1:nproma,:,:)   = xlm1_3d(1:nproma,:,:)+ QR(1:ie,:,1:ke,ntime)
       IF (itype_gscp > 1) THEN
          xim1_3d(1:nproma,:,:)   = xim1_3d(1:nproma,:,:) + QS(1:ie,:,1:ke,ntime)
          IF (itype_gscp == 4) xim1_3d(1:nproma,:,:) &
               = xim1_3d(1:nproma,:,:) + QG(1:ie,:,1:ke,ntime)
       ENDIF
    ENDIF

    qc_3d(1:ie,1:je,1:ke) = QC(1:ie,1:je,1:ke,ntime)
    IF (lprog_qi) &
         qi_3d(1:ie,1:je,1:ke) = QI(1:ie,1:je,1:ke,ntime)
    IF (lprogprec) THEN
       qr_3d(1:ie,1:je,1:ke) = QR(1:ie,1:je,1:ke,ntime)
       IF (itype_gscp > 1) &
            qs_3d(1:ie,1:je,1:ke) = QS(1:ie,1:je,1:ke,ntime)
       IF (itype_gscp ==4) &
            qg_3d(1:ie,1:je,1:ke) = QG(1:ie,1:je,1:ke,ntime)
    ENDIF
    IF (lmulti_layer) THEN
       t_so_3d(1:ie,1:je,1:ke_soil+2) = t_so(1:ie,1:je,0:ke_soil+1,ntime)
       w_so_3d(1:ie,1:je,1:ke_soil+1) = w_so(1:ie,1:je,1:ke_soil+1,ntime)
    ENDIF
    ! -- 2D 
    aps(1:ie,1:je)     = ps(1:ie,1:je,ntime)
    tsurf_2d(1:ie,1:je)  = t_s(1:ie,1:je,ntime)
    t_snow_2d(1:ie,1:je) = t_snow(1:ie,1:je,ntime)
    w_snow_2d(1:ie,1:je) = w_snow(1:ie,1:je,ntime)
    w_i_2d(1:ie,1:je)    = w_i(1:ie,1:je,ntime)
    qv_s_2d(1:ie,1:je)   = qv_s(1:ie,1:je,ntime)

    ! set pressure field
    press_3d(1:ie,1:je,1:ke)  = p0(1:ie,1:je,1:ke)+ pp(1:ie,1:je,1:ke,ntime)
    pressi_3d(1:ie,1:je,ke+1) = ps(1:ie,1:je,ntime)
    DO k=ke,1,-1
       DO i=1,ie
          DO j=1,je
             pressi_3d(i,j,k) = 2._dp*press_3d(i,j,k)-pressi_3d(i,j,k+1)
          ENDDO
       ENDDO
    ENDDO
    ! set geopotential fields
    !    define geopotential as zero at the surface
    !    => correct for height of surface
    geopoti_3d(:,:,ke+1) = 0._dp
    DO k=ke,1,-1
       geopoti_3d(:,:,k) = (hhl(:,:,k)-hhl(:,:,ke+1))*g
       geopot_3d(:,:,k) =  (hhl(:,:,k)+hhl(:,:,k+1)-2*hhl(:,:,ke+1))*g/2
    ENDDO
    ! calculate grid mass (grmass) in kg and grid volume (grvol) in m^3
    DO  k = 1, ke
       DO j = 1, je
          DO i = 1, ie
             sqrtg_r_s(i,j,k) = 1.0_dp / ( hhl(i,j,k)-hhl(i,j,k+1) )
             grvol(i,j,k)     = 1.0_dp / sqrtg_r_s(i,j,k) * gboxarea_2d(i,j)
!             grvol(i,j,k)     = 1.0_dp / sqrtg_r_s(i,j,k) *  crlat(j,1) &
!                  * r_earth**2 * (pi/180.0_dp)**2 * dlon * dlat * dzeta
             grmass(i,j,k) = rho(i,j,k) * grvol(i,j,k)
          ENDDO
       ENDDO
    ENDDO

    ! op_pj_20100713+
    grmassdry(:,:,:) = grmass(:,:,:) * (1.0_dp - qm1(:,:,:)) 
    ! op_pj_20100713-

    wind10_2d(:,:) = sqrt(v_10m(:,:)*v_10m(:,:) + u_10m(:,:) * u_10m(:,:))
    
    slm(:,:) = 0._dp           ! sea land mask
    tsw(:,:) = 273.15_dp
    tslm1(:,:) = t_s(:,:,ntime)
    WHERE (slf > 0.5_dp)
       slm = 1._dp
    ENDWHERE
    WHERE (slf < 0.5_dp .AND. t_g(:,:,ntime) > 0._dp )
       tsw = t_g(:,:,ntime)
    END WHERE
    ! surface roughness
    az0(:,:) = MAX(1.e-32_dp,gz0(:,:)/g)

    ! um_ak_20111011+
    ! ws  frational water content
    IF (lmulti_layer) THEN
       ws(:,:) = w_so(:,:,1,ntime) / (2._dp * czmls(1))
    ELSE
       ws(:,:) = w_g1(:,:,ntime) / (2._dp * czmls(1))
    ENDIF
    ! um_ak_20111011-

    ! set field capacity
    DO i=1,ie
       DO j=1,je
          wsmx(i,j) = cfcap(NINT(SOILTYP(i,j)))
       END DO
    END DO

    ! calculate soil moisture stress function (according to the
    ! calculation in ECHAM5 subroutine vdiff)
    pws(:,:) = ( ws(:,:)-plmin* wsmx(:,:) ) / (pldiff *  wsmx(:,:))
    fws(:,:) = MAX(0._dp, MIN(1._dp, pws(:,:)))

    ! calculate the virtual temperature (analogous to vdiff in ECHAM5)
    tvir(:,:) =  tm1_3d(:,:,ke) * (100000._dp/press_3d(:,:,ke))**c_vKar &
         * ( 1+ vtmpc1 * qm1_3d(:,:,ke) - (xlm1_3d(:,:,ke)+xim1_3d(:,:,ke) ))

    ! calculate relative vorticity
    dy = r_earth * (pi/180.0_dp) * dlon
    vom1(:,:,:) = 0._dp
    DO k = 1, ke
       DO j = 2, je-1 !jstart, jend
          dx = r_earth * (pi/180.0_dp) * dlon * crlat(j,1)
          DO i = 2, ie-1 !istart, iend
             vom1(i,j,k)= (( v(i+1,j,  k,ntime) + v(i+1,j-1,k,ntime) )        &
                  - (  v(i-1,j,  k,ntime) + v(i-1,j-1,k,ntime) ) )      &
                  * 0.5_dp / dx                               &
                  - (( u(i,  j+1,k,ntime) + u(i-1,j+1,k,ntime) )        &
                  - ( u(i,  j-1,k,ntime) + u(i-1,j-1,k,ntime) ) )       &
                  * 0.5_dp / dy
          END DO
        END DO
     END DO

!um_ak_20110704+
    ! CALCULATE RELATIVE HUMIDITY
    DO j=1, je
       ! EPSILON(1.) serves to avoid water vapour content in a layer
       !          of less than EPSILON(1.).
       zq(:,:,:)= MAX( qm1(:,:,:), EPSILON(1._dp) )
       !
       DO k = 1, ke
          DO i = 1, ie
             it = INT(tm1_3d(i,j,k)*1000._dp)
             IF (it<jptlucu1 .OR. it>jptlucu2) THEN 
                lookupoverflow = .TRUE.
                write (0,*) 'lookupoverflow ',i,j,k,tm1_3d(i,j,k)
             END IF
             it = MAX(MIN(it,jptlucu2),jptlucu1)
             zqs(i,j,k) = tlucua(it)/press_3d(i,j,k)
          END DO
       END DO
       IF (lookupoverflow) CALL error_bi('lookuperror', substr)
       zqs(1:ie,j,:)= MIN(zqs(1:ie,j,:),0.5_dp)
       zqs(1:ie,j,:)= zqs(1:ie,j,:) / (1._dp-vtmpc1*zqs(1:ie,j,:))
       zqs(1:ie,j,:)= MAX(2._dp*EPSILON(1._dp),zqs(1:ie,j,:))
       rhum_3d(1:ie,j,:) = 100._dp * zq(1:ie,j,:) / zqs(1:ie,j,:)
    END DO
!um_ak_20110704-
  END SUBROUTINE main_data_global_start
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE main_data_vdiff

    IMPLICIT NONE

    aclc(1:ie,1:je,1:ke) = clc_sgs(1:ie,1:je,1:ke)

    CALL calc_boundary_layer

  END SUBROUTINE main_data_vdiff
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE main_data_physc(jrow)

    INTEGER, INTENT(IN) :: jrow
    ! LOCAL
    INTEGER             :: ntime ! index for time level
                                 ! l2tls = T => ntime = nnow
                                 ! l2tls = F => ntime = nnold

    IF (l2tls) THEN
       ntime = nnow
    ELSE
       ntime = nold
    ENDIF
    ! If this model is MMD-Server some fields need to be
    ! set prior to calling MMDSERV
    qm1_3d(1:ie,jrow,1:ke)    = QV(1:ie,jrow,1:ke,ntime)
    tm1_3d(1:ie,jrow,1:ke)    = T(1:ie,jrow,1:ke,ntime)
    !qte_3d(1:ie,jrow,:)    = qvtens(1:ie,jrow,1:ke)
    !tte_3d(1:ie,jrow,:)    = ttens(1:ie,jrow,1:ke)
    IF (lprog_qi) THEN
       xim1_3d(1:ie,jrow,1:ke) = QI(1:ie,jrow,1:ke,ntime)
    ELSE 
       xim1_3d(1:ie,jrow,1:ke) = 0._dp
    ENDIF
    xlm1_3d(1:ie,jrow,1:ke) = QC(1:ie,jrow,1:ke,ntime)
    IF (lprogprec) THEN
       xlm1_3d(1:ie,jrow,1:ke)=xlm1_3d(1:ie,jrow,1:ke)+ QR(1:ie,jrow,1:ke,ntime)
       IF (itype_gscp > 1) THEN
          xim1_3d(1:ie,jrow,1:ke)= xim1_3d(1:ie,jrow,1:ke) &
               + QS(1:ie,jrow,1:ke,ntime)
          IF (itype_gscp == 4) xim1_3d(1:ie,jrow,1:ke) &
               = xim1_3d(1:ie,jrow,1:ke) + QG(1:ie,jrow,1:ke,ntime)
       ENDIF
    ENDIF
    aclc(1:ie,jrow,1:ke) = clc_sgs(1:ie,jrow,1:ke)

  END SUBROUTINE main_data_physc
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE main_data_2D_set_jrow(jrow)

    INTEGER, INTENT(IN) :: jrow
    ! LOCAL
    INTEGER             :: ntime ! index for time level
                                 ! l2tls = T => ntime = nnow
                                 ! l2tls = F => ntime = nnold
    IF (l2tls) THEN
       ntime = nnow
    ELSE
       ntime = nold
    ENDIF

    qte => qvtens(1:ie,jrow,:)
    xlte => qctens(1:ie,jrow,:)
    xite => qitens(1:ie,jrow,:)
    tte => ttens(1:ie,jrow,:)
    gboxarea => gboxarea_2d(1:ie,jrow)

  END SUBROUTINE main_data_2D_set_jrow
!-----------------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE HELPER ROUTINES 
  !  -> SAME AS IN MESSY_MAIN_CHANNEL_BI; DUE TO CIRCULAR DEPENDENCIES
  ! -------------------------------------------------------------------
  SUBROUTINE channel_halt(substr, status)

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,         ONLY: p_abort, p_pe
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

    CALL p_abort(substr, errstr)

  END SUBROUTINE channel_halt
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE calc_boundary_layer

    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_tools,          ONLY: jptlucu1, jptlucu2, tlucua
    USE messy_main_constants_mem,  ONLY: vtmpc1, vtmpc2, pi, cp_air, rd &
                                       , c_vkar, alv, als, tmelt, g
    USE messy_main_timer,          ONLY: time_step_len

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'calc_boundary_layer'

    INTEGER      :: it, it1        ! indices in lookup table
    REAL(dp)     :: es             ! saturation pressure
    REAL(dp)     :: esl(kproma)    ! saturation pressure over land
    REAL(dp)     :: qsl(kproma)    ! saturation ratio  over land
    REAL(dp)     :: pws(kproma)
    REAL(dp)     :: pwstop, pwslev 
    REAL(dp)     :: phum(kproma)
    REAL(dp)     :: zhsoil(kproma) ! soil humidity
    REAL(dp)     :: esw(kproma)    ! saturation pressure over water
    REAL(dp)     :: qsw(kproma)    ! saturation ratio  over water
    REAL(dp)     :: esi(kproma)    ! saturation pressure over ice
    REAL(dp)     :: qsi(kproma)    ! saturation ratio  over ice

    REAL(dp)     :: totclwater(kproma,nlev)  ! total cloud water
    REAL(dp)     :: zlteta1(kproma,nlev)     ! potential temperature ???
    REAL(dp)     :: faxe(kproma,nlev)        ! 
    REAL(dp)     :: zqss(kproma,nlev)        ! 
    REAL(dp)     :: zqddif
    REAL(dp)     :: zqlwi1, zqlwi2
    REAL(dp)     :: zteldif
    REAL(dp)     :: zdus1, zdus2
    REAL(dp)     :: zbuoy, zbet
    REAL(dp)     :: zfux, zfox, zusus1
    REAL(dp)     :: zmult1, zmult2, zmult3, zmult4, zmult5
    REAL(dp)     :: ztvir1(kproma,nlev)
    REAL(dp)     :: zfaxe(kproma,nlev)
    REAL(dp)     :: zdu2(kproma), zdu2oc(kproma)
    REAL(dp)     :: zucfl(kproma), zucfi(kproma), zucfw(kproma)
    REAL(dp)     :: zscfl(kproma), zscfi(kproma), zscfw(kproma)
    REAL(dp)     :: zvirmitte, ztemitte, zqmitte
    REAL(dp)     :: zcons, zalo
    REAL(dp)     :: zqsmit, ztmit, zqtmit

    ! fields here define locally (but global in ECHAM)
    REAL(dp)     :: ocu(kproma,je) 
    REAL(dp)     :: ocv(kproma,je) 
    REAL(dp)     :: tsi(kproma,je) 
    REAL(dp)     :: az0l(kproma,je)
    REAL(dp)     :: az0w(kproma,je)
    REAL(dp)     :: az0i(kproma,je)
    !REAL(dp)     :: aclc(kproma,je,ke)

    LOGICAL      :: lo

    INTEGER      :: jl, jk             ! loop index

    ! constant as set in echam5.3.xx 
    REAL(dp), PARAMETER :: kappa   = rd/cp_air
    REAL(dp), PARAMETER :: zrvrd   = 1._dp + vtmpc1  
    REAL(dp), PARAMETER :: zrdrv   = 1._dp/zrvrd
    REAL(dp)            :: const01 
    REAL(dp), PARAMETER :: const02 = 10._dp  ! from 0.2_dp * cb with cb = 5._dp 
    REAL(dp), PARAMETER :: const03 = 3._dp*5._dp*5._dp ! 3._dp * cb*cc

    !  set constants
    const01 = 1.5_dp * time_step_len * g / rd
    WHERE (az0 > 0._dp)
       az0l = az0
       az0w = az0
       az0i = az0
    ELSEWHERE
       ! set minimum value
       az0l = 1.e-10
       az0w = 1.e-10
       az0i = 1.e-10
    END WHERE
    !aclc(1:ie,1:je,1:ke) = clc_sgs(1:ie,1:je,1:ke)

    ocu(1:kproma,1:je) = 0._dp
    ocv(1:kproma,1:je) = 0._dp
    tsi(1:kproma,1:je) = t_s(1:ie,1:je,ntime)

    DO  jk=1,nlev
       DO  jl=1,kproma
          totclwater(jl,jk)=xlm1_3d(jl,jrow,jk)+xim1_3d(jl,jrow,jk) ! total cloud water
          tpot_3d(jl,jrow,jk)= &
               tm1(jl,jrow,jk)*(100000._dp/press_3d(jl,jrow,jk))**kappa
          ztvir1(jl,jk)=tpot_3d(jl,jrow,jk)* &
               (1._dp+vtmpc1*qm1(jl,jrow,jk)-totclwater(jl,jk))
          lo=tm1(jl,jrow,jk).GE.tmelt
          zfaxe(jl,jk)=MERGE(alv,als,lo)
          zbet=zfaxe(jl,jk)/cp_air
          zusus1=zbet*tpot_3d(jl,jrow,jk)/tm1(jl,jrow,jk)*totclwater(jl,jk)
          zlteta1(jl,jk)=tpot_3d(jl,jrow,jk)-zusus1
          it = NINT(tm1(jl,jrow,jk)*1000._dp)
          IF (it<jptlucu1 .OR. it>jptlucu2)  &
               CALL error_bi(substr, 'lookupoverflow 00')
          it = MAX(MIN(it,jptlucu2),jptlucu1)
          es=tlucua(it)/press_3d(jl,jrow,jk)
          es=MIN(es,0.5_dp)
          zqss(jl,jk)=es/(1._dp-vtmpc1*es)
       END DO
    END DO
    !
    !*      1.0   surface humidity and virtual temperature
    !                  for land, water and ice
    !
    DO jl=1,kproma
       !
       !    land ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !
       it = NINT(tslm1(jl,jrow)*1000._dp)
       IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          write (0,*) substr,jptlucu1,jptlucu2, 'it ' &
               ,it, jl, jrow, tslm1(jl,jrow)
          CALL error_bi(substr, 'lookupoverflow 01')
       END IF
       es=tlucua(it)/pressi_3d(jl,jrow,nlevp1)
       qsl(jl)=es/(1._dp-vtmpc1*es)

       it1= MAX(MIN(it+1,jptlucu2),jptlucu1)
       pws(jl)=MIN(ws(jl,jrow),wsmx(jl,jrow))
       pwstop=MIN(0.1_dp,wsmx(jl,jrow))
       pwslev=wsmx(jl,jrow)-pwstop
       IF(pws(jl).GT.pwslev.AND.pws(jl).GT. 0.35_dp*wsmx(jl,jrow)) THEN
          phum(jl)=0.5_dp*(1._dp-COS((pws(jl)-pwslev)*pi/pwstop))
       ELSE
          phum(jl)=0._dp
       END IF
       zhsoil(jl)= cvs(jl,jrow)+(1._dp-cvs(jl,jrow))                     &
            *(cvw(jl,jrow)+(1._dp-cvw(jl,jrow))*phum(jl))
       lo=qm1(jl,jrow,nlev).GT.qsl(jl)
       zhsoil(jl)=MERGE(1._dp,zhsoil(jl),lo)
       esl(jl)=tslm1(jl,jrow)*(1.e5_dp/pressi_3d(jl,jrow,nlevp1))**kappa
       tvl(jl,jrow)=esl(jl)*(1._dp+vtmpc1*zhsoil(jl)*qsl(jl))
       !
       !    water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !
       it = NINT(tsw(jl,jrow)*1000._dp)
       IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          write (0,*) substr, jptlucu1,jptlucu2, 'it ',it, tsw(jl,jrow)
          CALL error_bi(substr, 'lookupoverflow 02')
       END IF
       it = MAX(MIN(it,jptlucu2),jptlucu1)
       es=tlucua(it)/pressi_3d(jl,jrow,nlevp1)
       qsw(jl)=es/(1._dp-vtmpc1*es)
       esw(jl)=tsw(jl,jrow)*(1.e5_dp/pressi_3d(jl,jrow,nlevp1))**kappa
       tvw(jl,jrow)=esw(jl)*(1._dp+vtmpc1*qsw(jl))
       !
       !    ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !
       it = NINT(tsi(jl,jrow)*1000._dp)
       IF (it<jptlucu1 .OR. it>jptlucu2) THEN
          write (0,*) substr, jptlucu1,jptlucu2, 'it ',it, tsi(jl,jrow)
          CALL error_bi(substr, 'lookupoverflow 03')
       ENDIF
       it = MAX(MIN(it,jptlucu2),jptlucu1)
       es=tlucua(it)/pressi_3d(jl,jrow,nlevp1)
       qsi(jl)=es/(1._dp-vtmpc1*es)
       esi(jl)=tsi(jl,jrow)*(1.e5_dp/pressi_3d(jl,jrow,nlevp1))**kappa
       tvi(jl,jrow)=esi(jl)*(1._dp+vtmpc1*qsi(jl))
    END DO

    !     ------------------------------------------------------------------
    !
    !*         3.     COMPUTATION OF THE EXCHANGE COEFFICIENTS.
    !
    !        THE SURFACE LAYER IS NOW COMPUTED BEFORE THE OTHER LEVELS
    !
    !        3.1       COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
    !                  RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
    !                  AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
    !                  COMMON PART OF THE DRAG COEFFICIENTS.
    !
    DO jl=1,kproma
       !
       !     land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !
       zdu2(jl)=MAX(1.0_dp,um1(jl,jrow,nlev)**2+vm1(jl,jrow,nlev)**2)
       zqmitte=(qm1(jl,jrow,nlev)+qsl(jl)*zhsoil(jl))/2._dp
       zqtmit=totclwater(jl,nlev)*0.5_dp+zqmitte
       ztmit=(tm1(jl,jrow,nlev)+tslm1(jl,jrow))/2._dp
       zqsmit=(zqss(jl,nlev)+qsl(jl))/2._dp
       ztemitte=(tpot_3d(jl,jrow,nlev)+esl(jl))/2._dp
       zvirmitte=(ztvir1(jl,nlev)+tvl(jl,jrow))/2._dp
       zqlwi1=qm1(jl,jrow,nlev)+totclwater(jl,nlev)
       zqlwi2=qsl(jl)*zhsoil(jl)
       zfux=zfaxe(jl,nlev)/(cp_air*ztmit)
       zfox=zfaxe(jl,nlev)/(rd*ztmit)
       zmult1=1._dp+vtmpc1*zqtmit
       zmult2=zfux*zmult1-zrvrd
       zmult3=zrdrv*zfox*zqsmit/(1._dp+zrdrv*zfox*zfux*zqsmit)
       zmult5=zmult1-zmult2*zmult3
       zmult4=zfux*zmult5-1._dp
       zdus1=aclc(jl,jrow,nlev)*zmult5+(1._dp-aclc(jl,jrow,nlev))*zmult1
       zdus2=aclc(jl,jrow,nlev)*zmult4+(1._dp-aclc(jl,jrow,nlev))*vtmpc1
       zteldif=zlteta1(jl,nlev)-esl(jl)
       zqddif=zqlwi1-zqlwi2
       zbuoy=zdus1*zteldif+zdus2*ztemitte*zqddif
       ril(jl,jrow)=geopot_3d(jl,jrow,nlev)*zbuoy/(zvirmitte*zdu2(jl))
       cdnl(jl,jrow)=(c_vkar/LOG(1._dp+geopot_3d(jl,jrow,nlev)/ &
            (g*az0l(jl,jrow))))**2
       zucfl(jl)=1._dp/(1._dp+const03*cdnl(jl,jrow)*SQRT(ABS(ril(jl,jrow))      &
            *(1._dp +geopot_3d(jl,jrow,nlev)/(g*az0l(jl,jrow)))))
       zscfl(jl)=SQRT(1._dp+ABS(ril(jl,jrow)))
       zcons=const01*pressi_3d(jl,jrow,nlevp1)/                              &
            (tm1(jl,jrow,nlev)*(1._dp+vtmpc1*qm1(jl,jrow,nlev)-totclwater(jl,nlev)))
       cfncl(jl,jrow)=zcons*SQRT(zdu2(jl))*cdnl(jl,jrow)
       !
       !    water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !
       ! correction for water and ice points
       !
       zdu2oc(jl)=MAX(1.0_dp,(um1(jl,jrow,nlev)-ocu(jl,jrow))**2        &
            +(vm1(jl,jrow,nlev)-ocv(jl,jrow))**2)
       zqmitte=(qm1(jl,jrow,nlev)+qsw(jl))/2._dp
       zqtmit=totclwater(jl,nlev)*0.5_dp+zqmitte
       ztmit=(tm1(jl,jrow,nlev)+tsw(jl,jrow))/2._dp
       zqsmit=(zqss(jl,nlev)+qsw(jl))/2._dp
       ztemitte=(tpot_3d(jl,jrow,nlev)+esw(jl))/2._dp
       zvirmitte=(ztvir1(jl,nlev)+tvw(jl,jrow))/2._dp
       zqlwi1=qm1(jl,jrow,nlev)+totclwater(jl,nlev)
       zqlwi2=qsw(jl)
       zfux=zfaxe(jl,nlev)/(cp_air*ztmit)
       zfox=zfaxe(jl,nlev)/(rd*ztmit)
       zmult1=1._dp+vtmpc1*zqtmit
       zmult2=zfux*zmult1-zrvrd
       zmult3=zrdrv*zfox*zqsmit/(1._dp+zrdrv*zfox*zfux*zqsmit)
       zmult5=zmult1-zmult2*zmult3
       zmult4=zfux*zmult5-1._dp
       zdus1=aclc(jl,jrow,nlev)*zmult5+(1._dp-aclc(jl,jrow,nlev))*zmult1
       zdus2=aclc(jl,jrow,nlev)*zmult4+(1._dp-aclc(jl,jrow,nlev))*vtmpc1
       zteldif=zlteta1(jl,nlev)-esw(jl)
       zqddif=zqlwi1-zqlwi2
       zbuoy=zdus1*zteldif+zdus2*ztemitte*zqddif
       riw(jl,jrow)=geopot_3d(jl,jrow,nlev)*zbuoy/(zvirmitte*zdu2oc(jl))
       zalo=LOG(1._dp+geopot_3d(jl,jrow,nlev)/(g*az0w(jl,jrow)))
       cdnw(jl,jrow)=(c_vkar/zalo)**2
       zucfw(jl)=1._dp/(1._dp+const03*cdnw(jl,jrow)*SQRT(ABS(riw(jl,jrow))      &
            *(1._dp+geopot_3d(jl,jrow,nlev)/(g*az0w(jl,jrow)))))
       zscfw(jl)=SQRT(1._dp+ABS(riw(jl,jrow)))
       cfncw(jl,jrow)=zcons*SQRT(zdu2oc(jl))*cdnw(jl,jrow)
       !
       !     ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !
       zqmitte=(qm1(jl,jrow,nlev)+qsi(jl))/2._dp
       zqtmit=totclwater(jl,nlev)*0.5_dp+zqmitte
       ztmit=(tm1(jl,jrow,nlev)+tsi(jl,jrow))/2._dp
       zqsmit=(zqss(jl,nlev)+qsi(jl))/2._dp
       ztemitte=(tpot_3d(jl,jrow,nlev)+esi(jl))/2._dp
       zvirmitte=(ztvir1(jl,nlev)+tvi(jl,jrow))/2._dp
       zqlwi1=qm1(jl,jrow,nlev)+totclwater(jl,nlev)
       zqlwi2=qsi(jl)
       zfux=zfaxe(jl,nlev)/(cp_air*ztmit)
       zfox=zfaxe(jl,nlev)/(rd*ztmit)
       zmult1=1._dp+vtmpc1*zqtmit
       zmult2=zfux*zmult1-zrvrd
       zmult3=zrdrv*zfox*zqsmit/(1._dp+zrdrv*zfox*zfux*zqsmit)
       zmult5=zmult1-zmult2*zmult3
       zmult4=zfux*zmult5-1._dp
       zdus1=aclc(jl,jrow,nlev)*zmult5+(1._dp-aclc(jl,jrow,nlev))*zmult1
       zdus2=aclc(jl,jrow,nlev)*zmult4+(1._dp-aclc(jl,jrow,nlev))*vtmpc1
       zteldif=zlteta1(jl,nlev)-esi(jl)
       zqddif=zqlwi1-zqlwi2
       zbuoy=zdus1*zteldif+zdus2*ztemitte*zqddif
       rii(jl,jrow)=geopot_3d(jl,jrow,nlev)*zbuoy/(zvirmitte*zdu2oc(jl))
       zalo=LOG(1._dp+geopot_3d(jl,jrow,nlev)/(g*az0i(jl,jrow)))
       cdni(jl,jrow)=(c_vkar/zalo)**2
       zucfi(jl)=1._dp/(1._dp+const03*cdni(jl,jrow)*SQRT(ABS(rii(jl,jrow))     &
            *(1._dp+geopot_3d(jl,jrow,nlev)/(g*az0i(jl,jrow)))))
       zscfi(jl)=SQRT(1._dp+ABS(rii(jl,jrow)))
       cfnci(jl,jrow)=zcons*SQRT(zdu2oc(jl))*cdni(jl,jrow)
    END DO
    !
    !     3.2  DIMENSIONLESS HEAT TRANSFER COEFFICIENTS MULTIPLIED
    !          BY PRESSURE THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE
    !
    DO jl=1,kproma
       !
       !     land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !
       IF(ril(jl,jrow).GT.0._dp) THEN
          cfml(jl,jrow)=cfncl(jl,jrow)/(1._dp+const02*ril(jl,jrow)/zscfl(jl))
       ELSE
          cfml(jl,jrow)=cfncl(jl,jrow)*(1._dp-const02*ril(jl,jrow)*zucfl(jl))
       END IF
       !
       !     water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !
       IF(riw(jl,jrow).GT.0._dp) THEN
          cfmw(jl,jrow)=cfncw(jl,jrow)/(1._dp+const02*riw(jl,jrow)/zscfw(jl))
       ELSE
          cfmw(jl,jrow)=cfncw(jl,jrow)*(1._dp-const02*riw(jl,jrow)*zucfw(jl))
       END IF
       !
       !     ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       !
       IF(rii(jl,jrow).GT.0._dp) THEN
          cfmi(jl,jrow)=cfnci(jl,jrow)/(1._dp+const02*rii(jl,jrow)/zscfi(jl))
       ELSE
          cfmi(jl,jrow)=cfnci(jl,jrow)*(1._dp-const02*rii(jl,jrow)*zucfi(jl))
       END IF

    ENDDO
  END SUBROUTINE calc_boundary_layer
  ! -------------------------------------------------------------------
!!! END OF #ifdef COSMO
#endif

!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================

#ifdef BLANK

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather,       ONLY: start_message, end_message

  IMPLICIT NONE
  PUBLIC
  SAVE

  ! NAME AND VERSION OF THE BASEMODEL
  CHARACTER(LEN=*), PARAMETER :: modstr = 'BLANK'
  CHARACTER(LEN=*), PARAMETER :: modver = '1.0'

  LOGICAL, PARAMETER :: l2tls = .FALSE.      
  ! LOCALIZED PARAMETERS
  ! - GRID CONTROL HERE FIXED VALUES normally determined by grid definition
  INTEGER :: nlev = 2          ! number of vertical levels
  INTEGER :: nlon = 36
  INTEGER :: nlat = 18

  INTEGER :: nproma           ! vector length
  INTEGER :: ngpblks          ! number of vector rows
  INTEGER :: npromz           ! vector length of last row

  INTEGER :: kproma           ! vector length of current row
  INTEGER :: jrow             ! row loop index
  INTEGER :: je               !

  ! dummy
  LOGICAL, PARAMETER :: L_IS_CLIENT = .FALSE.

  ! defined constants
  REAL(dp) :: eps  = 0.1_dp

  ! exemplary variables for standard basemodel
  REAL(DP), POINTER :: yr,mo,dy,hr,mi,se,ms

CONTAINS

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_initialize

    IMPLICIT NONE

    nproma  = nlon
    ngpblks = nlat
    npromz  = nproma

  END SUBROUTINE main_data_initialize
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_init_memory

    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute
    USE messy_main_channel_repr,  ONLY: get_representation_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_init_memory'
    INTEGER :: status
    INTEGER :: reprid    

    ! create new channel
    CALL start_message(modstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

    CALL get_representation_id(status, 'SCALAR', reprid)
    CALL channel_halt(substr, status)

    CALL new_channel(status, modstr, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'channel_info' &
         , c = 'standard basemodel channel' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'yr', yr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'yr' &
         , 'long_name', c = 'year')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'mo', mo)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'mo' &
         , 'long_name', c = 'month')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'dy', dy)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dy' &
         , 'long_name', c = 'day')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'hr', hr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hr' &
         , 'long_name', c = 'hour')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'mi', mi)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'mi' &
         , 'long_name', c = 'minute')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'se', se)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'se' &
         , 'long_name', c = 'second')
    CALL channel_halt(substr, status)

    ! FORCE RESTART FILE CREATION
    CALL new_channel_object(status, modstr, 'ms', ms, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ms' &
         , 'long_name', c = 'millisecond')
    CALL channel_halt(substr, status)

    CALL end_message(modstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

  END SUBROUTINE main_data_init_memory
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_global_start

    ! BML/MESSy
    USE messy_main_timer,        ONLY: time_step_len

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, TRIM, ADJUSTL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_global_start'
    INTEGER :: status
    CHARACTER (8)   :: ydate
    CHARACTER (10)  :: ytime
    
    INTEGER :: iyr, imo, idy, ihr, imi, ise, ims

    CALL DATE_AND_TIME(ydate, ytime)

    READ(ydate,'(i4)') iyr
    READ(ydate,'(4x,i2)') imo
    READ(ydate,'(6x,i2)') idy

    READ(ytime,'(i2)') ihr
    READ(ytime,'(2x,i2)') imi
    READ(ytime,'(4x,i2)') ise
    READ(ytime,'(7x,i3)') ims

    yr = REAL(iyr,dp)
    mo = REAL(imo,dp)
    dy = REAL(idy,dp)
    hr = REAL(ihr,dp)
    mi = REAL(imi,dp)
    se = REAL(ise,dp)
    ms = REAL(ims,dp)

  END SUBROUTINE main_data_global_start
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_local_start

    IMPLICIT NONE

    ! jrow is set in main program loop (region loop)

    IF ( jrow == ngpblks ) THEN
       kproma = npromz
    ELSE
       kproma = nproma
    END IF

  END SUBROUTINE main_data_local_start
  !----------------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE HELPER ROUTINES 
  !  -> SAME AS IN MESSY_MAIN_CHANNEL_BI; DUE TO CIRCULAR DEPENDENCIES
  ! -------------------------------------------------------------------
  SUBROUTINE channel_halt(substr, status)

    ! MESSy
    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_constants_mem,  ONLY: STRLEN_VLONG
    USE messy_main_channel_error,  ONLY: channel_error_str

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: substr
    INTEGER,          INTENT(IN)  :: status
    ! LOCAL
    CHARACTER(LEN=STRLEN_VLONG)   :: errstr

    IF (status == 0) RETURN

    errstr = channel_error_str(status)

    CALL error_bi(errstr, substr)

  END SUBROUTINE channel_halt
  ! -------------------------------------------------------------------

!!! END OF #ifdef BLANK
#endif

!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================

! mz_ab_20100308+
#ifdef MBM_CMAT

  USE messy_main_constants_mem, ONLY: dp
  USE messy_cmat_grid,          ONLY: &
       nlon, nlat,                    &
       nlev => ht_dim,                &
       ht_dim, lon_dim, lat_dim,      &
       rlat, phi, pres
  USE messy_cmat_InputParamsPlanet, ONLY: cmat_idt_O3   => idt_O3  &
                                        , cmat_idt_H2O  => idt_H2O &
                                        , cmat_idt_H2O2 => idt_H2O2 
  USE messy_cmat_globedata, ONLY: R0_eff

  IMPLICIT NONE
  PUBLIC
  SAVE

  ! NAME AND VERSION OF THE BASEMODEL
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'CMAT2'
  CHARACTER(len=*), PARAMETER, PUBLIC :: modver = '1.4'

  ! RESTART CYCLE COUNTER
  CHARACTER(LEN=4) :: cystr = '0000'
  ! 
  LOGICAL, PARAMETER :: l2tls = .FALSE.      
  REAL(dp) :: eps  = 0.1_dp
!!$  ! LOCALIZED PARAMETERS
!!$
  INTEGER :: nproma           ! length of x-direction 
  INTEGER :: npromz  
  INTEGER :: ngpblks          ! length of y-direction 
!!$
  INTEGER :: kproma           ! defined in "local" loops (ie)
  INTEGER :: jrow             ! defined in "local" loops (je)
!!$  INTEGER :: je               ! defined in "local" loops (je)
  
  LOGICAL :: lmidatm = .TRUE.
  REAL(dp), ALLOCATABLE :: vct(:) ! vertical coefficients table
  INTEGER :: nlevp1 &
           , nvclev      ! number of levels with vertical coefficients
  REAL(dp) :: apzero = 1013.*1e2 ! surface pressure
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: philat, philon
  INTEGER :: ngl
  ! DUMMY VARIABLES
  INTEGER, PARAMETER :: nn = 42 ! max meridional wave number for m=0.
  LOGICAL, PARAMETER :: L_IS_CLIENT = .FALSE.

CONTAINS
  !***************************************************************************
  SUBROUTINE main_data_initialize

    USE messy_cmat_grid, ONLY: theta, phi, lat_dim, lon_dim, pres
    USE messy_cmat_globedata, ONLY: temp3d

    IMPLICIT NONE

    ! nlev    =       ! model number vertical level
    ! nlon    =       ! number longitude
    ! nlat    =       ! number latitudes
    ! nproma  = nlon
    ! ngpblks = nlat
    !npromz  = nproma 
    !nglon   = nlon
    !nglat   = nlat
    !kepin   = root%kepin
    !kezin   = root%kezin
    !kproma  = nproma 

    ALLOCATE(vct(nvclev*2))
    ALLOCATE(philat(lat_dim))
    ALLOCATE(philon(lon_dim))
    philat = theta ! deg latitude
    philon = phi   ! deg longitude
    nlevp1 = nlev+1
    nvclev = nlevp1     ! number of levels with vertical coefficients
    ngl = lat_dim

  END SUBROUTINE main_data_initialize
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  SUBROUTINE main_data_init_memory
    IMPLICIT NONE
  END SUBROUTINE main_data_init_memory
  SUBROUTINE main_data_global_start
    IMPLICIT NONE
  END SUBROUTINE main_data_global_start
#endif
! mz_ab_20100308-

! op_bk_20130820+
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================

#ifdef __ICON__

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather,       ONLY: start_message, end_message

  IMPLICIT NONE
  PUBLIC
  SAVE

  ! NAME AND VERSION OF THE BASEMODEL
  CHARACTER(LEN=*), PARAMETER :: modstr = 'ICON'
  CHARACTER(LEN=*), PARAMETER :: modver = '1.0.2'

  INTEGER :: nlev
  INTEGER :: nproma, kproma
  INTEGER :: npromz
  INTEGER :: ngpblks

  REAL(DP) :: eps
  INTEGER  :: jrow

  LOGICAL, PARAMETER :: l2tls = .TRUE.

CONTAINS

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_initialize

    USE messy_main_timer,    ONLY: timer_set_time_step_len

    IMPLICIT NONE

    CALL timer_set_time_step_len(.TRUE.)

  END SUBROUTINE main_data_initialize
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_init_memory

    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute
    USE messy_main_channel_repr,  ONLY: get_representation_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_init_memory'
    INTEGER :: status
    INTEGER :: reprid    

    ! create new channel
    CALL start_message(modstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

    CALL end_message(modstr,'CREATE STANDARD CHANNELS/OBJECTS',substr)

  END SUBROUTINE main_data_init_memory
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_global_start

    ! BML/MESSy
    USE messy_main_timer,        ONLY: timer_set_time_step_len

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, TRIM, ADJUSTL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_global_start'
    INTEGER :: status

    CALL timer_set_time_step_len(.TRUE.)

  END SUBROUTINE main_data_global_start
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE main_data_local_start

    IMPLICIT NONE

    ! jrow is set in main program loop (region loop)


  END SUBROUTINE main_data_local_start
  !----------------------------------------------------------------------------

  !***************************************************************************
  SUBROUTINE main_data_free_memory

    IMPLICIT NONE

    ! CLEANUP NON-CHANNEL OBJECT MEMORY

  END SUBROUTINE main_data_free_memory
  !***************************************************************************

  ! -------------------------------------------------------------------
  ! PRIVATE HELPER ROUTINES 
  !  -> SAME AS IN MESSY_MAIN_CHANNEL_BI; DUE TO CIRCULAR DEPENDENCIES
  ! -------------------------------------------------------------------
  SUBROUTINE channel_halt(substr, status)

    ! MESSy
    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_constants_mem,  ONLY: STRLEN_VLONG
    USE messy_main_channel_error,  ONLY: channel_error_str

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: substr
    INTEGER,          INTENT(IN)  :: status
    ! LOCAL
    CHARACTER(LEN=STRLEN_VLONG)   :: errstr

    IF (status == 0) RETURN

    errstr = channel_error_str(status)

    CALL error_bi(errstr, substr)

  END SUBROUTINE channel_halt
  ! -------------------------------------------------------------------

!!! END OF #ifdef __ICON__
#endif

!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================

! op_bk_20130820-

!*****************************************************************************
END MODULE messy_main_data_bi
!*****************************************************************************
