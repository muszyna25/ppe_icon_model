!>
!! cloud microphysics
!!
!! !------------------------------------------------------------------------------
!!
!! @par Description of *graupel*:
!!   This module procedure calculates the rates of change of temperature, cloud
!!   water, cloud ice, water vapor, rain, snow and graupel due to cloud microphysical
!!   processes related to the formation of grid scale precipitation. This
!!   includes the sedimentation of rain and snow. The precipitation fluxes at
!!   the surface are also calculated here.
!!
!! Method:
!!   Prognostic one-moment bulk microphysical parameterization.
!!   The sedimentation of rain and snow is computed implicitly.
!!
!! Current Code Owner: DWD, A. Seifert
!!    phone: +49-69-8062-2729,  fax:   +49-69-8062-3721
!!    email: Axel.Seifert@dwd.de
!!
!! @author Guenther Doms
!! @author Axel Seifert
!!
!! @par reference   This is an adaption of subroutine hydci_pp_gr in file src_gscp.f90
!!  of the COSMO-Model. Equation numbers refer to
!!  Doms, Foerstner, Heise, Herzog, Raschendorfer, Schrodin, Reinhardt, Vogel
!!    (September 2005): "A Description of the Nonhydrostatic Regional Model LM",
!!
!------------------------------------------------------------------------------
!!
!! $Id: n/a$
!!
!! @par Revision History
!! implemented into ICON by Felix Rieper (2012-06)
!! 
!! @par Revision History
!! Graupel scheme implemented into ICON by Carmen Koehler (2014-03)
!! Modifications from Felix Rieper for modifications to improve supercooled 
!! liquid water (SLW) prediction in cloudice
!!     - reduced number of ice crystal Ni(T), now according to Cooper (1986)
!!     - reduced rain freezing rate srfrzr, now according to Bigg (1953)
!!     - reduced depositional growth of ice and snow (zsidep, zssdep) 
!!       at cloud top, according to R. Forbes (2013) -> IFS model!  
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!

MODULE gscp_graupel

!------------------------------------------------------------------------------
!>
!! Description:
!!
!!   This subroutine in the module "gscp_graupel" calculates the rates of change of
!!   temperature, cloud condensate and water vapor due to cloud microphysical
!!   processes related to the formation of grid scale clouds and precipitation.
!!
!==============================================================================
!
! Declarations:
!
! Modules used:

!------------------------------------------------------------------------------
! Microphysical constants and variables
!------------------------------------------------------------------------------

#ifdef __COSMO__
USE kind_parameters, ONLY :   &
  wp ,       &    !! KIND-type parameter for real variables
  i4              !! KIND-type parameter for standard integer vaiables

USE data_constants  , ONLY :   &
!   pi,           & !!
!! 2. physical constants and related variables
!! -------------------------------------------
    t0=>t0_melt,  & !! melting temperature of ice
!   r_d,          & !! gas constant for dry air
    r_v,          & !! gas constant for water vapour
    rdv,          & !! r_d / r_v
    o_m_rdv,      & !! 1 - r_d/r_v
    rvd_m_o,      & !! r_v/r_d - 1
    cp_d,         & !! specific heat of dry air
    cpdr,         & !! 1 / cp_d
    lh_v,         & !! latent heat of vapourization
!   lh_f,         & !! latent heat of fusion
    lh_s,         & !! latent heat of sublimation
!   g,            & !! acceleration due to gravity
!   rho_w,        & !! density of liquid water (kg/m^3)
! 3. constants for parametrizations
!! ---------------------------------
    b1,           & !! variables for computing the saturation vapour pressure
    b2w,          & !! over water (w) and ice (i)
    b2i,          & !!               -- " --
    b3,           & !!               -- " --
    b4w,          & !!               -- " --
    b234w,        & !!               -- " --
    b4i,          & !!               -- " --
!> 4. tuning constants for radiation, cloud physics, turbulence
!! ------------------------------------------------------------
    qi0,          & !! cloud ice threshold for autoconversion
    qc0,          & !! cloud water threshold for autoconversion
!> 5. Precision-dependent security parameters (epsilons)
!! ------------------------------------------------------
    repsilon        !! precision of 1.0 in current floating point format

!! end of data_constants

USE data_runcontrol , ONLY :   &
    ldiabf_lh,      & ! include diabatic forcing due to latent heat in RK-scheme
    lsuper_coolw,   & ! switch for supercooled liquid water
    lsppt,          & ! switch, if .true., perturb the physical tendencies
    itype_qxpert_rn,& ! define which hum variables tend. are perturbed
    itype_qxlim_rn    ! type of reduction/removal of the perturbation 
                      ! in case of negative (qv, qc, qi) or 
                      ! supersaturated (qv) values

!------------------------------------------------------------------------------
!! COSMO environment modules
!------------------------------------------------------------------------------

USE utilities,                ONLY : message
USE meteo_utilities,          ONLY : satad       !! saturation adjustment
USE src_stoch_physics,        ONLY : apply_tqx_tend_adj
#endif
! of ifdef __COSMO__

!------------------------------------------------------------------------------

#ifdef NUDGING
USE data_lheat_nudge,           ONLY :  &
    llhn,         & ! main switch for latent heat nudging
    llhnverif,    & ! main switch for latent heat nudging
    lhn_qrs,      & ! use integrated precipitaion flux as reference
    tt_lheat,     & ! latent heat release of the model
    qrsflux         ! total precipitation flux
#endif

!------------------------------------------------------------------------------

#ifdef __ICON__
USE mo_kind,               ONLY: wp         , &
                                 i4
!USE mo_math_constants    , ONLY: pi
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                                 r_d   => rd    , & !! gas constant for dry air
                                 rvd_m_o=>vtmpc1, & !! r_v/r_d - 1
                                 o_m_rdv        , & !! 1 - r_d/r_v
                                 rdv            , & !! r_d / r_v
                                 lh_v  => alv   , & !! latent heat of vapourization
                                 lh_s  => als   , & !! latent heat of sublimation
!                                lh_f  => alf   , & !! latent heat of fusion
                                 cp_d  => cpd   , & !! specific heat of dry air at constant press
                                 cpdr  => rcpd  , & !! (spec. heat of dry air at constant press)^-1
                                 cvdr  => rcvd  , & !! (spec. heat of dry air at const vol)^-1
                                 b3    => tmelt , & !! melting temperature of ice/snow
!                                rho_w => rhoh2o, & !! density of liquid water (kg/m^3)
!                                g     => grav  , & !! acceleration due to gravity
                                 t0    => tmelt     !! melting temperature of ice/snow

USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b2i   => c3ies , & !!               -- " --
                                 b4w   => c4les , & !!               -- " --
                                 b4i   => c4ies , & !!               -- " --
                                 b234w => c5les     !!               -- " --
USE mo_satad,              ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
                                 sat_pres_ice!,   &  !! saturation vapor pressure w.r.t. ice
!                                 spec_humi          !! Specific humidity 
USE mo_exception,          ONLY: message, message_text
#endif

!------------------------------------------------------------------------------

! this can be used by ICON and COSMO
USE gscp_data, ONLY: &          ! all variables are used here

    ccsrim,    ccsagg,    ccsdep,    ccsvel,    ccsvxp,    ccslam,       &
    ccslxp,    ccsaxp,    ccsdxp,    ccshi1,    ccdvtp,    ccidep,       &
    ccswxp,    zconst,    zcev,      zbev,      zcevxp,    zbevxp,       &
    zvzxp,     zvz0r,                                                    &

    v0snow => v0snow_gr,  mu_rain,   rain_n0_factor,       cloud_num,    &

    x13o8,     x1o2,      x27o16,    x3o4,      x7o4,      x7o8,         &
    zbms,      zbvi,      zcac,      zccau,     zciau,     zcicri,       &
    zcrcri,    zcrfrz,    zcrfrz1,   zcrfrz2,   zeps,      zkcac,        &
    zkphi1,    zkphi2,    zkphi3,    zmi0,      zmimax,    zmsmin,       &
    zn0s0,     zn0s1,     zn0s2,     znimax_thom,          zqmin,        &
    zrho0,     zthet,     zthn,      ztmix,     ztrfrz,    zv1s,         &
    zvz0i,     x13o12,    x2o3,      x5o24,     zasmel,                  &
    zbsmel,    zcsmel,    x1o3,      zams => zams_gr,                    &
    iautocon,  isnow_n0temp, dist_cldtop_ref,   reduce_dep_ref,          &
    tmin_iceautoconv,     zceff_fac, zceff_min,                          &
    mma, mmb

#ifdef __ICON__
! this is (at the moment) an ICON part
USE gscp_data, ONLY: &          ! all variables are used here
    vtxexp,    & !
    kc_c1,     & !
    kc_c2,     & !
    kc_alpha,  & !
    kc_beta,   & !
    kc_gamma,  & !
    kc_sigma,  & !
    do_i,      & !
    co_i
#endif

!==============================================================================

IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: graupel

LOGICAL, PARAMETER :: &
  lsedi_ice    = .TRUE. , &  ! switch for sedimentation of cloud ice (Heymsfield & Donner 1990 *1/3)
  lstickeff    = .TRUE. , &  ! switch for sticking coeff. (work from Guenther Zaengl)
  lsuper_coolw = .TRUE. , &  ! switch for supercooled liquid water (work from Felix Rieper)
  lred_depgrow = .TRUE.      ! separate switch for reduced depositional growth near tops of water clouds
                             ! (now also used in ICON after correcting the cloud top diagnosis)
!------------------------------------------------------------------------------
!> Parameters and variables which are global in this module
!------------------------------------------------------------------------------

#ifdef __COSMO__
CHARACTER(132) :: message_text = ''
#endif

!==============================================================================

CONTAINS

!==============================================================================
!> Module procedure "graupel" in "gscp_graupel" for computing effects of grid
!! scale precipitation including cloud water, cloud ice, graupel, rain and snow
!------------------------------------------------------------------------------

SUBROUTINE graupel     (             &
  nvec,ke,                           & !> array dimensions
  ivstart,ivend, kstart,             & !! optional start/end indicies
  idbg,                              & !! optional debug level
  zdt, dz,                           & !! numerics parameters
  t,p,rho,qv,qc,qi,qr,qs,qg,qnc,     & !! prognostic variables
#ifdef __ICON__
  !xxx: this should become a module variable, e.g. in a new module mo_gscp_data.f90
  qi0,qc0,                           & !! cloud ice/water threshold for autoconversion
#endif
  prr_gsp,prs_gsp,prg_gsp,           & !! surface precipitation rates
#ifdef __COSMO__
  tinc_lh,                           & !  t-increment due to latent heat 
  pstoph,                            & !  stochastic multiplier of physics tendencies
#endif
#ifdef NUDGING
  tt_lheat,                          & !  t-increments due to latent heating (nud) 
  qrsflux,                           & !  total precipitation flux
#endif
  l_cv,                              &
  ddt_tend_t     , ddt_tend_qv     , &
  ddt_tend_qc    , ddt_tend_qi     , & !> ddt_tend_xx are tendencies
  ddt_tend_qr    , ddt_tend_qs     , & !!    necessary for dynamics
  ddt_diag_au    , ddt_diag_ac     , & !!
  ddt_diag_ev    , ddt_diag_nuc    , & !! ddt_diag_xxx are optional
  ddt_diag_idep  , ddt_diag_sdep   , & !!   diagnostic tendencies of all
  ddt_diag_agg   , ddt_diag_rim    , & !!   individual microphysical processes
  ddt_diag_rcri  , ddt_diag_icri   , & !!
  ddt_diag_dau   , ddt_diag_iau    , & !!
  ddt_diag_imelt , ddt_diag_smelt  , & !!
  ddt_diag_cfrz  , ddt_diag_rfrz   , & !!
  ddt_diag_shed  , ddt_tend_qg       ) !!

!------------------------------------------------------------------------------
! Description:
!   This module procedure calculates the rates of change of temperature, cloud
!   water, cloud ice, water vapor, rain, snow, and graupel due to cloud
!   microphysical processes related to the formation of grid scale
!   precipitation.
!   The variables are updated in this subroutine. Rain, snow, and graupel are
!   prognostic variables. The precipitation fluxes at the surface are stored
!   on the corresponding global fields.
!   The subroutine relies mostly on conversion terms used in cloudice
!   In contrast to cloudice, qc0 = 0.0002 (instead of 0.0) is used
!   This is version G29TK21.
!
! Method:
!   The sedimentation of rain, snow, and graupel is computed implicitly.
!
! Vectorization:
!   Most computations in this routine are grouped in IF-clauses. But the IFs
!   inside DO-loops often hinder or even disables vectorization. 
!   For the big IF-chunks, the condition is now checked at the beginning of 
!   the subroutine and the corresponding indices are stored in extra index
!   arrays. The following loops then are running only over these indices,
!   avoiding at least some IF-clauses inside the loops.
!
!------------------------------------------------------------------------------
!! Declarations:
!!
!------------------------------------------------------------------------------
!! Modules used: These are declared in the module declaration section
!! -------------

!! Subroutine arguments:
!! --------------------

  INTEGER, INTENT(IN) ::  &
    nvec          ,    & !> number of horizontal points
    ke                     !! number of grid points in vertical direction

  INTEGER, INTENT(IN), OPTIONAL ::  &
    ivstart   ,    & !> optional start index for horizontal direction
    ivend     ,    & !! optional end index   for horizontal direction
    kstart    ,    & !! optional start index for the vertical index
    idbg             !! optional debug level

  REAL(KIND=wp), INTENT(IN) :: &
    zdt                    !> time step for integration of microphysics     (  s  )

#ifdef __ICON__
  REAL(KIND=wp), INTENT(IN) :: &
    qi0,qc0          !> cloud ice/water threshold for autoconversion
#endif

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(in)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(IN) ::      &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
#endif
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv                   !! if true, cv is used instead of cp

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(INOUT) ::   &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
#endif
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    qr              ,    & !! specific rain content                         (kg/kg)
    qs              ,    & !! specific snow content                         (kg/kg)
    qg                     !! specific graupel content                      (kg/kg)

#ifdef __COSMO__    
  REAL(KIND=wp), INTENT(INOUT) :: &
       tinc_lh(:,:)    ! temperature increments due to heating             ( K/s )   

  REAL(KIND=wp), INTENT(IN)    :: &
       pstoph(:,:)     ! stochastic multiplier of physics tendencies
#endif

#ifdef NUDGING
  REAL(KIND=wp), INTENT(INOUT) :: &
       tt_lheat(:,:)  ,  & !  t-increments due to latent heating (nudg) ( K/s )
       qrsflux(:,:)       ! total precipitation flux (nudg)
#endif

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=wp), DIMENSION(nvec), INTENT(INOUT) ::   &
#else
  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
#endif
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    prs_gsp,             & !! precipitation rate of snow, grid-scale        (kg/(m2*s))
    prg_gsp,             & !! precipitation rate of graupel, grid-scale     (kg/(m2*s))
    qnc                    !! cloud number concentration

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(out)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(OUT), OPTIONAL ::   &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &     ! dim (ie,ke)
#endif
    ddt_tend_t      , & !> tendency T                                       ( 1/s )
    ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
    ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
    ddt_tend_qi     , & !! tendency qi                                      ( 1/s )
    ddt_tend_qr     , & !! tendency qr                                      ( 1/s )
    ddt_tend_qs     , & !! tendency qs                                      ( 1/s )
    ddt_tend_qg         !! tendency qg                                      ( 1/s )

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(out)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=wp), DIMENSION(nvec,ke), INTENT(OUT), OPTIONAL ::   &
#else
  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &   ! dim (ie,ke)
#endif
    ddt_diag_au     , & !> optional output autoconversion rate cloud to rain           ( 1/s )
    ddt_diag_ac     , & !! optional output accretion rate cloud to rain                ( 1/s )
    ddt_diag_ev     , & !! optional output evaporation of rain                         ( 1/s )
    ddt_diag_nuc    , & !! optional output mass nucleation of cloud ice                ( 1/s )
    ddt_diag_idep   , & !! optional output depositional growth of cloud ice            ( 1/s )
    ddt_diag_sdep   , & !! optional output depositional growth of snow                 ( 1/s )
    ddt_diag_agg    , & !! optional output aggregation snow collecting cloud ice       ( 1/s )
    ddt_diag_rim    , & !! optional output riming of snow by cloud water               ( 1/s )
    ddt_diag_rcri   , & !! optional output cloud ice + rain -> snow (rcri is sink qr)  ( 1/s )
    ddt_diag_icri   , & !! optional output cloud ice + rain -> snow (icri is sink qi)  ( 1/s )
    ddt_diag_dau    , & !! optional output depositional cloud ice autoconversion       ( 1/s )
    ddt_diag_iau    , & !! optional output aggregational cloud ice autoconversion      ( 1/s )
    ddt_diag_imelt  , & !! optional output melting of cloud ice                        ( 1/s )
    ddt_diag_smelt  , & !! optional output melting of snow                             ( 1/s )
    ddt_diag_cfrz   , & !! optional output freezing of cloud water                     ( 1/s )
    ddt_diag_rfrz   , & !! optional output rainwater freezing                          ( 1/s )
    ddt_diag_shed       !! optional output shedding                                    ( 1/s )


  !! Local parameters: None, parameters are in module header, gscp_data or data_constants
  !! ----------------
  
  !> Local scalars:
  !! -------------
  
  INTEGER (KIND=i4)        ::  &
    iv, k             !> loop indices

  REAL    (KIND=wp   ) :: nnr

  REAL    (KIND=wp   ) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)

  INTEGER ::  &
    iv_start     ,    & !> start index for horizontal direction
    iv_end       ,    & !! end index for horizontal direction
    k_start      ,    & !! model level where computations start
    izdebug             !! debug level

  REAL    (KIND=wp   ), PARAMETER ::  &
    zams  = 0.038_wp,     & ! Formfactor in the mass-size relation of snow particles
    zcsg=0.5_wp,          & !coefficient for snow-graupel conversion by riming
    zcrim_g=4.43_wp,      & !
    zrimexp_g=0.94878_wp, &
    zcagg_g = 2.46_wp ,   & !
    zasmel= 2.95E3_wp ,   & ! DIFF*lh_v*RHO/LHEAT
    zexpsedg=0.217_wp,    & ! exponent for graupel sedimentation
    zvz0g = 12.24_wp  ,   & ! coefficient of sedimentation velocity for graupel
    ztcrit=3339.5_wp        ! factor in calculation of critical temperature

 
  REAL    (KIND=wp   ) ::  &
    fpvsw, fpvsi, fqvs,& ! name of statement functions
    fxna ,             & ! statement function for ice crystal number
    fxna_cooper ,      & ! statement function for ice crystal number, Cooper(1986) 
    ztx  , zpv  , zpx ,& ! dummy arguments for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    znimix,            & ! number of ice crystals at ztmix -> threshold temp for mixed-phase clouds 
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
    zqvsw0diff,        & ! qv-zqvsw0  
    zdtr ,             & ! reciprocal of timestep for integration
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum,            & ! terms for limiting  total rain water depletion
    znin,              & ! number of cloud ice crystals at nucleation
    fnuc,              & !FR: coefficient needed for Forbes (2012) SLW layer parameterization 
    znid,              & ! number of cloud ice crystals for deposition
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,& ! terms for limiting total cloud ice depletion
    zqvsw,             & ! sat. specitic humidity at ice and water saturation
    zqvsidiff,         & ! qv-zqvsi
    ztfrzdiff,         & ! ztrfrz-t  
    zztau, zxfac, zx1,  ztt,  &   ! some help variables
    ztau, zphi, zhi, zdvtp, ztc, zeff, zlog_10

  REAL    (KIND=wp   ) ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice
    zqrt   ,& ! layer tendency of rain
    zqst   ,& ! layer tendency of snow
    zqgt      ! layer tendency of graupel

  REAL    (KIND=wp   ) ::  &
    zlnqrk,zlnqsk,zlnqik,     & !
    zlnlogmi,zlnqgk,ccswxp_ln1o2,zvzxp_ln1o2,zbvi_ln1o2,zexpsedg_ln1o2, &
    qcg,tg,qvg,qrg,qsg,qgg,qig,rhog,ppg,alf,bet,m2s,m3s,hlp,            &
    qcgk_1,maxevap,temp_c

  LOGICAL :: &
    llqs,llqc,llqi,llqg  !   switch for existence of qr, qs, qc, qi

  LOGICAL :: &
    llqr


  REAL(KIND=wp), DIMENSION(nvec,ke) ::   &
    t_in               ,    & !> temperature                                   (  K  )
    qv_in              ,    & !! specific water vapor content                  (kg/kg)
    qc_in              ,    & !! specific cloud water content                  (kg/kg)
    qi_in              ,    & !! specific cloud ice   content                  (kg/kg)
    qr_in              ,    & !! specific rain content                         (kg/kg)
    qs_in              ,    & !! specific snow content                         (kg/kg)
    qg_in                     !! specific graupel content                      (kg/kg)



!! Local (automatic) arrays:
!! -------------------------

  REAL    (KIND=wp   ) ::  &
    zqvsi             ,     & !> sat. specitic humidity at ice and water saturation
    zvzr        (nvec),     & !
    zvzs        (nvec),     & !
    zvzg        (nvec),     & ! 
    zvzi        (nvec),     & ! terminal fall velocity of ice
    zpkr        (nvec),     & !
    zpks        (nvec),     & !
    zpkg        (nvec),     & ! 
    zpki        (nvec),     & ! precipitation flux of ice
    zprvr       (nvec),     & !
    zprvs       (nvec),     & !
    zprvg       (nvec),     & !
    zprvi       (nvec),     & !
#ifdef __COSMO__
    zdummy      (nvec,8),   & !
#endif
    zcsdep            ,     & !
    zcidep            
    
 REAL    (KIND=wp   ) ::  &    
    zsrmax            ,     & !
    zssmax            ,     & !
    zsgmax            ,     & !
    zvz0s             ,     & !
    zcrim             ,     & !
    zcagg             ,     & !
    zbsdep            ,     & !
    zcslam            ,     & !
    zn0s              ,     & !
    zimr              ,     & !
    zims              ,     & !
    zimg              ,     & !
    zimi              ,     & !
    zzar              ,     & !
    zzas              ,     & !
    zzag              ,     & !
    zzai              ,     & !
    zqrk              ,     & !
    zqsk              ,     & !
    zqgk              ,     & !
    zqik              ,     & !
    zdtdh             ,     & !
    z1orhog           ,     & ! 1/rhog
    zrho1o2           ,     & ! (rho0/rhog)**1/2
    zrho1o3           ,     & ! (rho0/rhog)**1/3
    zeln7o8qrk        ,     & !
    zeln7o4qrk        ,     & ! FR new  
    zeln27o16qrk      ,     & !
    zeln13o8qrk       ,     & !
    zeln3o4qsk        ,     & ! 
    zeln6qgk          ,     & !
    zeln8qsk          ,     & !
    zelnrimexp_g      


  REAL    (KIND=wp   ) ::  &
    scau   , & ! transfer rate due to autoconversion of cloud water
    scac   , & ! transfer rate due to accretion of cloud water
    snuc   , & ! transfer rate due nucleation of cloud ice
    scfrz  , & ! transfer rate due homogeneous freezing of cloud water
    simelt , & ! transfer rate due melting of cloud ice
    sidep  , & ! transfer rate due depositional growth of cloud ice
    ssdep  , & ! transfer rate due depositional growth of snow
    sgdep  , & ! transfer rate due depositional growth of snow
    sdau   , & ! transfer rate due depositional cloud ice autoconversion
    srim   , & ! transfer rate due riming of snow
    srim2  , & ! transfer rate due riming of snow
    sconsg , & ! transfer rate due to conversion from snow to graupel by riming  
    sshed  , & ! transfer rate due shedding
    sicri  , & ! transfer rate due cloud ice collection by rain (sink qi)
    srcri  , & ! transfer rate due cloud ice collection by rain (sink qr)
    sagg   , & ! transfer rate due aggregation of snow and cloud ice
    sagg2  , & ! transfer rate due aggregation of snow and cloud ice
    siau   , & ! transfer rate due autoconversion of cloud ice
    ssmelt , & ! transfer rate due melting of snow
    sgmelt , & ! transfer rate due melting of snow
    sev    , & ! transfer rate due evaporation of rain
    sconr  , & ! transfer rate due to condensation on melting snow/graupel
    srfrz  , & ! transfer rate due to rainwater freezing
    reduce_dep,&!FR: coefficient: reduce deposition at cloud top (Forbes 2012)
    dist_cldtop(nvec) !FR: distance from cloud top layer 

  LOGICAL :: ldum

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine graupel
!------------------------------------------------------------------------------

!> Statement functions
! -------------------

! saturation vapour pressure over water (fpvsw), over ice (fpvsi)
! and specific humidity at vapour saturation (fqvs)
  fpvsw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) )
  fpvsi(ztx)     = b1*EXP( b2i*(ztx-b3)/(ztx-b4i) )
  fqvs (zpv,zpx) = rdv*zpv/( zpx - o_m_rdv*zpv )
  
! Number of activate ice crystals;  ztx is temperature
  fxna(ztx)   = 1.0E2_wp * EXP(0.2_wp * (t0 - ztx))
  fxna_cooper(ztx) = 5.0E+0_wp * EXP(0.304_wp * (t0 - ztx))   ! FR: Cooper (1986) used by Greg Thompson(2008)


! Define reciprocal of heat capacity of dry air (at constant pressure vs at constant volume)

#ifdef __COSMO__
  z_heat_cap_r = cpdr
#endif

#ifdef __ICON__
  IF (PRESENT(l_cv)) THEN
    IF (l_cv) THEN
      z_heat_cap_r = cvdr
    ELSE
      z_heat_cap_r = cpdr
    ENDIF
  ELSE
    z_heat_cap_r = cpdr
  ENDIF
#endif

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local and global variables
!------------------------------------------------------------------------------

! Some constant coefficients
  IF( lsuper_coolw) THEN
    znimax = znimax_Thom
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSE
    znimax = fxna(zthn) ! Maximum number of cloud ice crystals
    znimix = fxna(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  END IF

  zpvsw0 = fpvsw(t0)  ! sat. vap. pressure for t = t0
  zlog_10 = LOG(10._wp) ! logarithm of 10
  
  ! Precomputations for optimization
  ccswxp_ln1o2   = EXP (ccswxp * LOG (0.5_wp))
  zvzxp_ln1o2    = EXP (zvzxp * LOG (0.5_wp))
  zbvi_ln1o2     = EXP (zbvi * LOG (0.5_wp))
  zexpsedg_ln1o2 = EXP (zexpsedg * LOG (0.5_wp))

! Optional arguments

  IF (PRESENT(ddt_tend_t)) THEN
    ! save input arrays for final tendency calculation
    t_in  = t
    qv_in = qv
    qc_in = qc
    qi_in = qi
    qr_in = qr
    qs_in = qs
    qg_in = qg
  END IF
  IF (PRESENT(ivstart)) THEN
    iv_start = ivstart
  ELSE
    iv_start = 1
  END IF
  IF (PRESENT(ivend)) THEN
    iv_end = ivend
  ELSE
    iv_end = nvec
  END IF
  IF (PRESENT(kstart)) THEN
    k_start = kstart
  ELSE
    k_start = 1
  END IF
  IF (PRESENT(idbg)) THEN
    izdebug = idbg
  ELSE
    izdebug = 0
  END IF


! timestep for calculations
  zdtr  = 1.0_wp / zdt

! output for various debug levels
  IF (izdebug > 15) CALL message('gscp_graupel: ', 'Start of hydci_pp_gr')
  IF (izdebug > 20) THEN
    WRITE (message_text,*) '   nvec = ',nvec       ; CALL message('',message_text)
    WRITE (message_text,*) '   ke = ',ke           ; CALL message('',message_text)
    WRITE (message_text,*) '   ivstart = ',ivstart ; CALL message('',message_text)
    WRITE (message_text,*) '   ivend   = ',ivend   ; CALL message('',message_text)
  END IF
  IF (izdebug > 50) THEN
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN dz  = ',MAXVAL(dz),MINVAL(dz)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN T   = ',MAXVAL(t),MINVAL(t)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN p   = ',MAXVAL(p),MINVAL(p)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN rho = ',MAXVAL(rho),MINVAL(rho)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qv  = ',MAXVAL(qv),MINVAL(qv)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qc  = ',MAXVAL(qc),MINVAL(qc)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qr  = ',MAXVAL(qr),MINVAL(qr)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qi  = ',MAXVAL(qi),MINVAL(qi)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qs  = ',MAXVAL(qs),MINVAL(qs)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qg  = ',MAXVAL(qg),MINVAL(qg)
    CALL message('',message_text) 
  ENDIF

  ! Delete precipitation fluxes from previous timestep
  DO iv = iv_start, iv_end
    prr_gsp (iv) = 0.0_wp
    prs_gsp (iv) = 0.0_wp
    prg_gsp (iv) = 0.0_wp
    zpkr(iv)     = 0.0_wp
    zpks(iv)     = 0.0_wp
    zpkg(iv)     = 0.0_wp
    zpki(iv)     = 0.0_wp
    zprvr(iv)    = 0.0_wp
    zprvs(iv)    = 0.0_wp
    zprvg(iv)    = 0.0_wp
    zprvi(iv)    = 0.0_wp
    zvzr(iv)     = 0.0_wp
    zvzs(iv)     = 0.0_wp
    zvzg(iv)     = 0.0_wp
    zvzi(iv)     = 0.0_wp
    dist_cldtop(iv) = 0.0_wp
  END DO

! *********************************************************************
! Loop from the top of the model domain to the surface to calculate the
! transfer rates  and sedimentation terms
! *********************************************************************

  loop_over_levels: DO  k = k_start, ke

    DO iv = iv_start, iv_end  !loop over horizontal domain


#ifdef NUDGING
      ! add part of latent heating calculated in subroutine graupel to model latent
      ! heating field: subtract temperature from model latent heating field
      IF (llhn .OR. llhnverif) THEN
        IF (lhn_qrs) THEN
          qrsflux(iv,k) = 0.0_wp
        ENDIF
        ! replaces: CALL get_gs_lheating ('add',1,ke)
        tt_lheat(iv,k) = tt_lheat(iv,k) - t(iv,k)
      ENDIF
#endif

#ifdef __COSMO__
      IF ( ldiabf_lh ) THEN
        ! initialize temperature increment due to latent heat
        tinc_lh(iv,k) = tinc_lh(iv,k) - t(iv,k)
      ENDIF
#endif

      !----------------------------------------------------------------------------
      ! Section 2: Check for existence of rain and snow
      !            Initialize microphysics and sedimentation scheme
      !----------------------------------------------------------------------------

      zcrim  = 0.0_wp
      zcagg  = 0.0_wp
      zbsdep = 0.0_wp
      zvz0s  = 0.0_wp
      zn0s   = zn0s0
      reduce_dep = 1.0_wp  !FR: Reduction coeff. for dep. growth of rain and ice  

      !----------------------------------------------------------------------------
      ! 2.1: Preparations for computations and to check the different conditions
      !----------------------------------------------------------------------------

      qrg = qr(iv,k)
      qsg = qs(iv,k)
      qgg = qg(iv,k)
      qvg = qv(iv,k)
      qcg = qc(iv,k)
      qig = qi(iv,k)
      tg  = t(iv,k)
      ppg = p(iv,k)
      rhog = rho(iv,k)

      !..for density correction of fall speeds
      z1orhog = 1.0_wp/rhog
      hlp     = LOG(zrho0*z1orhog)
      zrho1o2 = EXP(hlp*x1o2)
      zrho1o3 = EXP(hlp*x1o3)

      zqrk = qrg * rhog
      zqsk = qsg * rhog
      zqgk = qgg * rhog
      zqik = qig * rhog

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqg = zqgk > zqmin
      llqi = zqik > zqmin

      zdtdh = 0.5_wp * zdt / dz(iv,k)

      zzar = zqrk/zdtdh + zprvr(iv) + zpkr(iv)
      zzas = zqsk/zdtdh + zprvs(iv) + zpks(iv)
      zzag = zqgk/zdtdh + zprvg(iv) + zpkg(iv)
      zzai = zqik/zdtdh + zprvi(iv) + zpki(iv)

      zpkr(iv) = 0.0_wp
      zpks(iv) = 0.0_wp
      zpkg(iv) = 0.0_wp
      zpki(iv) = 0.0_wp

      !-------------------------------------------------------------------------
      ! qs_prepare:
      !-------------------------------------------------------------------------
      IF (llqs) THEN
        IF (isnow_n0temp == 1) THEN
          ! Calculate n0s using the temperature-dependent
          ! formula of Field et al. (2005)
          ztc = tg - t0
          ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
          zn0s = zn0s1*EXP(zn0s2*ztc)
          zn0s = MIN(zn0s,1.0E9_wp)
          zn0s = MAX(zn0s,1.0E6_wp)
        ELSEIF (isnow_n0temp == 2) THEN
          ! Calculate n0s using the temperature-dependent moment
          ! relations of Field et al. (2005)
          ztc = tg - t0
          ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)

          nnr  = 3._wp
          hlp = mma(1) + mma(2)*ztc + mma(3)*nnr + mma(4)*ztc*nnr &
              + mma(5)*ztc**2 + mma(6)*nnr**2 + mma(7)*ztc**2*nnr &
              + mma(8)*ztc*nnr**2 + mma(9)*ztc**3 + mma(10)*nnr**3
          alf = EXP(hlp*zlog_10) ! 10.0_wp**hlp
          bet = mmb(1) + mmb(2)*ztc + mmb(3)*nnr + mmb(4)*ztc*nnr &
              + mmb(5)*ztc**2 + mmb(6)*nnr**2 + mmb(7)*ztc**2*nnr &
              + mmb(8)*ztc*nnr**2 + mmb(9)*ztc**3 + mmb(10)*nnr**3

          ! Here is the exponent bms=2.0 hardwired! not ideal! (Uli Blahak)
          m2s = qsg * rhog / zams   ! UB rho added as bugfix
          m3s = alf*EXP(bet*LOG(m2s))

          hlp  = zn0s1*EXP(zn0s2*ztc)
          zn0s = 13.50_wp * m2s * (m2s / m3s)**3
          zn0s = MAX(zn0s,0.5_wp*hlp)
          zn0s = MIN(zn0s,1.0E2_wp*hlp)
          zn0s = MIN(zn0s,1.0E9_wp)
          zn0s = MAX(zn0s,1.0E6_wp)
        ELSE
          ! Old constant n0s
          zn0s = 8.0E5_wp
        ENDIF
        zcrim  = ccsrim*zn0s
        zcagg  = ccsagg*zn0s
        zbsdep = ccsdep*SQRT(v0snow)
        zvz0s  = ccsvel*EXP(ccsvxp * LOG(zn0s))
        zlnqsk = zvz0s * EXP (ccswxp * LOG (zqsk)) * zrho1o2
        zpks(iv) = zqsk * zlnqsk
        IF (zvzs(iv) == 0.0_wp) THEN
          zvzs(iv) = zlnqsk * ccswxp_ln1o2
        ENDIF
      ENDIF ! qs_prepare
    
      ! sedimentation fluxes

      !-------------------------------------------------------------------------
      ! qr_sedi:
      !-------------------------------------------------------------------------

      IF (llqr) THEN
        zlnqrk = zvz0r * EXP (zvzxp * LOG (zqrk)) * zrho1o2
        zpkr(iv) = zqrk * zlnqrk
        IF (zvzr(iv) == 0.0_wp) THEN
          zvzr(iv) = zlnqrk * zvzxp_ln1o2
        ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      ! qg_sedi:
      !-------------------------------------------------------------------------

      IF (llqg) THEN
        zlnqgk = zvz0g * EXP (zexpsedg * LOG (zqgk)) * zrho1o2
        zpkg(iv) = zqgk * zlnqgk
        IF (zvzg(iv) == 0.0_wp) THEN
          zvzg(iv) = zlnqgk * zexpsedg_ln1o2
        ENDIF
      ENDIF ! qg_sedi

      !-------------------------------------------------------------------------
      ! qi_sedi:
      !-------------------------------------------------------------------------

      IF (llqi) THEN
        zlnqik = zvz0i * EXP (zbvi * LOG (zqik)) * zrho1o3
        zpki(iv) = zqik * zlnqik
        IF (zvzi(iv) == 0.0_wp) THEN
          zvzi(iv) = zlnqik * zbvi_ln1o2
        ENDIF
      ENDIF  ! qi_sedi

      !--------------------------------------------------------------------------
      ! 2.3: Second part of preparations
      !--------------------------------------------------------------------------

      zeln7o8qrk    = 0.0_wp
      zeln7o4qrk    = 0.0_wp
      zeln27o16qrk  = 0.0_wp
      zeln13o8qrk   = 0.0_wp
      zeln3o4qsk    = 0.0_wp
      zeln8qsk      = 0.0_wp
      zeln6qgk      = 0.0_wp
      zelnrimexp_g  = 0.0_wp
      zsrmax        = 0.0_wp
      zssmax        = 0.0_wp
      zsgmax        = 0.0_wp

      !FR old
      !   zcsdep    = 3.2E-2_wp
      zcsdep        = 3.367E-2_wp
      zcidep        = 1.3E-5_wp
      zcslam        = 1e10_wp

      scau          = 0.0_wp
      scac          = 0.0_wp
      snuc          = 0.0_wp
      scfrz         = 0.0_wp
      simelt        = 0.0_wp
      sidep         = 0.0_wp
      ssdep         = 0.0_wp
      sgdep         = 0.0_wp
      sdau          = 0.0_wp
      srim          = 0.0_wp
      srim2         = 0.0_wp
      sshed         = 0.0_wp
      sicri         = 0.0_wp
      srcri         = 0.0_wp
      sagg          = 0.0_wp
      sagg2         = 0.0_wp
      siau          = 0.0_wp
      ssmelt        = 0.0_wp
      sgmelt        = 0.0_wp
      sev           = 0.0_wp
      sconr         = 0.0_wp
      sconsg        = 0.0_wp
      srfrz         = 0.0_wp

      zpkr(iv)   = MIN( zpkr(iv) , zzar )
      zpks(iv)   = MIN( zpks(iv) , zzas )
      zpkg(iv)   = MIN( zpkg(iv) , MAX(0._wp,zzag) )
      zpki(iv)   = MIN( zpki(iv) , zzai )

      zzar   = zdtdh * (zzar-zpkr(iv))
      zzas   = zdtdh * (zzas-zpks(iv))
      zzag   = zdtdh * (zzag-zpkg(iv))
      zzai   = zdtdh * (zzai-zpki(iv))

      zimr   = 1.0_wp / (1.0_wp + zvzr(iv) * zdtdh)
      zims   = 1.0_wp / (1.0_wp + zvzs(iv) * zdtdh)
      zimg   = 1.0_wp / (1.0_wp + zvzg(iv) * zdtdh)
      zimi   = 1.0_wp / (1.0_wp + zvzi(iv) * zdtdh)

      zqrk   = zzar*zimr
      zqsk   = zzas*zims
      zqgk   = zzag*zimg
      zqik   = zzai*zimi

#ifdef __COSMO__
      zqvsi = fqvs( fpvsi(tg), ppg )
#endif

#ifdef __ICON__
      zqvsi = sat_pres_ice(tg)/(rhog * r_v * tg)
#endif

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqg = zqgk > zqmin
      llqc =  qcg > zqmin
      llqi =  qig > zqmin

      !!----------------------------------------------------------------------------
      !! 2.4: IF (llqr): ic1
      !!----------------------------------------------------------------------------

      IF (llqr) THEN
        zlnqrk   = LOG (zqrk)
        zsrmax   = zzar/rhog*zdtr  ! GZ: shifting this computation ahead of the IF condition changes results!
        IF ( qig+qcg > zqmin ) THEN
          zeln7o8qrk   = EXP (x7o8   * zlnqrk)
        ENDIF
        IF ( tg < ztrfrz ) THEN
          zeln7o4qrk   = EXP (x7o4   * zlnqrk) !FR new
          zeln27o16qrk = EXP (x27o16 * zlnqrk)
        ENDIF
        IF (llqi) THEN
          zeln13o8qrk  = EXP (x13o8  * zlnqrk)
        ENDIF
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.5: IF (llqs): ic2
      !!----------------------------------------------------------------------------

! ** GZ: the following computation differs substantially from the corresponding code in cloudice **
      IF (llqs) THEN
        zlnqsk   = LOG (zqsk)
        zssmax   = zzas / rhog*zdtr  ! GZ: shifting this computation ahead of the IF condition changes results!
        IF (qig+qcg > zqmin) THEN
          zeln3o4qsk = EXP (x3o4 *zlnqsk)
        ENDIF
        zeln8qsk = EXP (0.8_wp *zlnqsk)
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.6: IF (llqg): ic3
      !!----------------------------------------------------------------------------

      IF (zqgk > zqmin) THEN
        zlnqgk   = LOG (zqgk)
        zsgmax   = zzag / rhog*zdtr
        IF (qig+qcg > zqmin) THEN
          zelnrimexp_g = EXP (zrimexp_g * zlnqgk)
        ENDIF
        zeln6qgk = EXP (0.6_wp *zlnqgk)
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.7:  slope of snow PSD and coefficients for depositional growth (llqi,llqs)
      !!----------------------------------------------------------------------------    

      IF ((qig > zqmin) .OR. (zqsk > zqmin)) THEN
        zdvtp  = ccdvtp * EXP(1.94_wp * LOG(tg)) / ppg
        zhi    = ccshi1*zdvtp*rhog*zqvsi/(tg*tg)
        hlp    = zdvtp / (1.0_wp + zhi)
        zcidep = ccidep * hlp

        IF (llqs) THEN
          zcslam = EXP(ccslxp * LOG(ccslam * zn0s / zqsk ))
          zcslam = MIN(zcslam,1.0E15_wp)
          zcsdep = 4.0_wp * zn0s * hlp
        ENDIF
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.8: Deposition nucleation for low temperatures below a threshold (llqv)
      !!----------------------------------------------------------------------------    

      IF (( tg < zthet .AND. qvg >  8.E-6_wp &
                       .AND. qig <= 0.0_wp )) THEN
        IF( qvg > zqvsi ) THEN
          IF( lsuper_coolw) THEN
            znin  = MIN( fxna_cooper(tg), znimax )
          ELSE
            znin  = MIN( fxna(tg), znimax )
          END IF
          snuc = zmi0 * z1orhog * znin * zdtr
        ENDIF
      ENDIF

      !!--------------------------------------------------------------------------
      !! Section 3: Search for cloudy grid points with cloud water and
      !!            calculation of the conversion rates involving qc (ic6)
      !!--------------------------------------------------------------------------

      IF (qcg > zqmin) THEN
        llqs = zqsk > zqmin

        zscmax = qcg*zdtr
        IF( tg > zthn ) THEN
          IF (iautocon == 0) THEN
            ! Kessler (1969) autoconversion rate
            scau = zccau * MAX( qcg - qc0, 0.0_wp )
            scac = zcac  * qcg * zeln7o8qrk
          ELSEIF (iautocon == 1) THEN
            ! Seifert and Beheng (2001) autoconversion rate
            ! with constant cloud droplet number concentration qnc
            IF (qcg > 1.0E-6_wp) THEN
              ztau = MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp)
              ztau = MAX(ztau,1.E-30_wp)
              hlp  = EXP(zkphi2*LOG(ztau))
              zphi = zkphi1 * hlp * (1.0_wp - hlp)**3
              scau = zconst * qcg*qcg*qcg*qcg/(qnc(iv)*qnc(iv)) &
                   * (1.0_wp + zphi/(1.0_wp - ztau)**2)
              zphi = (ztau/(ztau+zkphi3))**4
              scac = zkcac * qcg * qrg * zphi
            ELSE
              scau = 0.0_wp
              scac = 0.0_wp
            ENDIF
          ENDIF
          IF (llqr) THEN
            ! Calculation of in-cloud rainwater freezing
            IF ( tg < ztrfrz .AND. qrg > 0.1_wp*qcg ) THEN
              IF (lsuper_coolw) THEN
                srfrz = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk
              ELSE
                ztfrzdiff=ztrfrz-tg
                srfrz = zcrfrz*ztfrzdiff*SQRT(ztfrzdiff)* zeln27o16qrk
              ENDIF
            ENDIF
          ENDIF
          IF (llqs) THEN
            srim = zcrim * qcg *  EXP(ccsaxp * LOG(zcslam))
          ENDIF
          srim2 = zcrim_g * qcg * zelnrimexp_g
          IF( tg >= t0 ) THEN
            sshed = srim+srim2
            srim  = 0.0_wp
            srim2 = 0.0_wp
          ELSE
            IF (qcg >= qc0) THEN
              sconsg = zcsg * qcg * zeln3o4qsk
            ENDIF
          ENDIF
          ! Check for maximum depletion of cloud water and adjust the
          ! transfer rates accordingly
          zscsum = scau + scac + srim + srim2 + sshed
          zcorr  = zscmax / MAX( zscmax, zscsum )
          scau   = zcorr*scau
          scac   = zcorr*scac
          srim   = zcorr*srim
          srim2  = zcorr*srim2
          sshed  = zcorr*sshed
          sconsg = MIN (sconsg, srim+zssmax)
        ELSE !tg >= tg: ! hom. freezing of cloud and rain water
          scfrz = zscmax
          srfrz = zsrmax
        ENDIF
        ! Calculation of heterogeneous nucleation of cloud ice.
        ! This is done in this section, because we require water saturation
        ! for this process (i.e. the existence of cloud water) to exist.
        ! Heterogeneous nucleation is assumed to occur only when no
        ! cloud ice is present and the temperature is below a nucleation
        ! threshold.

        IF( tg <= 267.15_wp .AND. .NOT.llqi ) THEN   
          IF (lsuper_coolw) THEN
            znin  = MIN( fxna_cooper(tg), znimax )
            snuc = zmi0 * z1orhog * znin * zdtr
          ELSE
            znin = MIN( fxna(tg), znimax )
            snuc = zmi0 / rhog * znin * zdtr
          END IF
        ENDIF
        ! Calculation of reduction of depositional growth at cloud top (Forbes 2012)
        IF( k>1 .AND. k<ke .AND. lred_depgrow ) THEN
          znin = MIN(fxna_cooper(tg), znimax )
          fnuc = MIN(znin/znimix, 1.0_wp)

          qcgk_1 = qc(iv,k-1) + qi(iv,k-1) + qs(iv,k-1)

          !! distance from cloud top
          IF( qcgk_1 .LT. zqmin ) THEN      ! upper cloud layer
            dist_cldtop(iv) = 0.0_wp    ! reset distance to upper cloud layer
          ELSE
            dist_cldtop(iv) = dist_cldtop(iv) + dz(iv,k)
          END IF

          ! with asymptotic behaviour dz -> 0 (xxx)
          !        reduce_dep = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
          !                             dist_cldtop(iv)/dist_cldtop_ref + &
          !                             (1.0_wp-reduce_dep_ref)*(zdh/dist_cldtop_ref)**4), 1.0_wp)

          ! without asymptotic behaviour dz -> 0
          reduce_dep = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
                        dist_cldtop(iv)/dist_cldtop_ref), 1.0_wp)

        END IF ! Reduction of dep. growth of snow/ice 

      ENDIF

      !------------------------------------------------------------------------
      ! Section 4: Search for cold grid points with cloud ice and/or snow and
      !            calculation of the conversion rates involving qi, qs and qg
      !------------------------------------------------------------------------

      IF ( (qig > zqmin) .OR. (zqsk > zqmin) .OR. zqgk > zqmin ) THEN
        llqs =  zqsk > zqmin
        llqi =   qig > zqmin

        IF (tg<=t0) THEN           ! cold case 

          zqvsidiff = qvg-zqvsi
          zsvmax    = zqvsidiff * zdtr
          IF (llqi) THEN

            IF( lsuper_coolw) THEN
              znin   = MIN( fxna_cooper(tg), znimax )
            ELSE
              znin   = MIN( fxna(tg), znimax )
            END IF
            ! Change in sticking efficiency needed in case of cloud ice sedimentation
            ! (based on Guenther Zaengls work)
            IF (lstickeff) THEN
              zeff     = MIN(EXP(0.09_wp*(tg-t0)),1.0_wp)
              zeff     = MAX(zeff, zceff_min, zceff_fac*(tg-tmin_iceautoconv)) 
            ELSE !original sticking efficiency of cloud ice
              zeff     = MIN(EXP(0.09_wp*(tg-t0)),1.0_wp)
              zeff     = MAX(zeff,0.2_wp)
            END IF
            sagg      = zeff * qig * zcagg * EXP(ccsaxp*LOG(zcslam))
            sagg2     = zeff * qig * zcagg_g * zelnrimexp_g
            siau      = zeff * zciau * MAX( qig - qi0, 0.0_wp )
            zmi       = MIN( rhog*qig/znin, zmimax )
            zmi       = MAX( zmi0, zmi )
            znid      = rhog * qig/zmi
            zlnlogmi  = LOG (zmi)
            sidep     = zcidep * znid * EXP(0.33_wp * zlnlogmi) * zqvsidiff
            zsvidep   = 0.0_wp
            zsvisub   = 0.0_wp
            zsimax    = qig*zdtr
            IF( sidep > 0.0_wp ) THEN
              IF (lred_depgrow ) THEN
                sidep = sidep * reduce_dep  !FR new: depositional growth reduction
              END IF
              zsvidep = MIN( sidep, zsvmax )
            ELSEIF ( sidep < 0.0_wp ) THEN
              zsvisub  =   MAX (   sidep,  zsvmax)
              zsvisub  = - MAX ( zsvisub, -zsimax)
            ENDIF
            zlnlogmi   = LOG  (zmsmin/zmi)
            zztau      = 1.5_wp*( EXP(0.66_wp*zlnlogmi) - 1.0_wp)
            sdau       = zsvidep/zztau
            sicri      = zcicri * qig * zeln7o8qrk
            IF (qsg > 1.e-7_wp) srcri = zcrcri * (qig/zmi) * zeln13o8qrk
          ELSE
            zsimax    =  0.0_wp
            zsvidep   =  0.0_wp
            zsvisub   =  0.0_wp
          ENDIF

          zxfac = 1.0_wp + zbsdep * EXP(ccsdxp*LOG(zcslam))
          ssdep = zcsdep * zxfac * zqvsidiff / (zcslam+zeps)**2
          !FR new: depositional growth reduction
          IF (lred_depgrow .AND. ssdep > 0.0_wp) THEN
            ssdep = ssdep*reduce_dep
          END IF
          ! GZ: This limitation, which was missing in the original graupel scheme,
          ! is crucial for numerical stability in the tropics!
          IF (ssdep > 0.0_wp) ssdep = MIN(ssdep, zsvmax-zsvidep)
          ! Suppress depositional growth of snow if the existing amount is too small for a
          ! a meaningful distiction between cloud ice and snow
          IF (qsg <= 1.e-7_wp) ssdep = MIN(ssdep, 0.0_wp)
! ** GZ: this numerical fit should be replaced with a physically more meaningful formulation **
          sgdep = (0.398561_wp-0.00152398_wp*tg                 &
                   + 2554.99_wp/ppg+ 2.6531E-7_wp*ppg) *        &
                 zqvsidiff * zeln6qgk
          ! Check for maximal depletion of cloud ice
          ! No check is done for depositional autoconversion because
          ! this is a always a fraction of the gain rate due to
          ! deposition (i.e the sum of this rates is always positive)
          zsisum = siau + sagg + sagg2 + sicri + zsvisub
          zcorr  = 0.0_wp
          IF( zsimax > 0.0_wp ) zcorr  = zsimax / MAX( zsimax, zsisum )
          sidep  = zsvidep - zcorr*zsvisub
          siau   = zcorr*siau
          sagg   = zcorr*sagg
          sagg2  = zcorr*sagg2
          sicri  = zcorr*sicri
          IF ( zqvsidiff < 0.0_wp ) THEN
            ssdep = MAX(ssdep, - zssmax)
            sgdep = MAX(sgdep, - zsgmax)
          ENDIF

        ELSE ! tg > 0 - warm case

          !------------------------------------------------------------------------
          ! Section 5: Search for warm grid points with cloud ice and/or snow and
          !            calculation of the melting rates of qi and ps
          !------------------------------------------------------------------------

          ! cloud ice melts instantaneously
          simelt = qig*zdtr

#ifdef __COSMO__
          zqvsw0     = fqvs( zpvsw0, ppg)
#endif

#ifdef __ICON__
          zqvsw0     = zpvsw0/(rhog * r_v *t0)
#endif
          zqvsw0diff = qvg-zqvsw0
! ** GZ: several numerical fits in this section should be replaced with physically more meaningful formulations **
          IF ( tg > (t0-ztcrit*zqvsw0diff) ) THEN
            !calculate melting rate
            zx1         = (tg - t0) + zasmel*zqvsw0diff
            ssmelt = (79.6863_wp/ppg+0.612654E-3_wp)* zx1 * zeln8qsk
            ssmelt = MIN (ssmelt,zssmax)
            sgmelt = (12.31698_wp/ppg+7.39441e-05_wp)* zx1 * zeln6qgk
            sgmelt = MIN (sgmelt, zsgmax)
            !deposition + melting, ice particle temperature: t0
            !calculation without howell-factor!
            ssdep  = (31282.3_wp/ppg+0.241897_wp)       &
                    * zqvsw0diff * zeln8qsk
            sgdep  = (0.153907_wp-ppg*7.86703e-07_wp)  &
                    * zqvsw0diff * zeln6qgk
            IF (zqvsw0diff < 0.0_wp) THEN
              !melting + evaporation of snow/graupel
              ssdep = MAX (-zssmax,ssdep)
              sgdep = MAX (-zsgmax,sgdep)
              !melt water evaporates
              ssmelt = ssmelt+ssdep
              sgmelt = sgmelt+sgdep
              ssmelt = MAX( ssmelt, 0.0_wp )
              sgmelt = MAX( sgmelt, 0.0_wp )
            ELSE
              !deposition on snow/graupel is interpreted as increase
              !in rain water ( qv --> qr, sconr)
              !therefore,  sconr=(zssdep+zsgdep)
              sconr=ssdep+sgdep
              ssdep=0.0_wp
              sgdep=0.0_wp
            ENDIF
          ELSE
            !if t<t_crit
            !no melting, only evaporation of snow/graupel
#ifdef __COSMO__
            zqvsw      = fqvs( fpvsw(tg), ppg )
#endif

#ifdef __ICON__
            zqvsw      = sat_pres_water(tg)/(rhog * r_v *tg)
#endif
            zqvsidiff  = qvg-zqvsw
            ssdep = (0.28003_wp-ppg*0.146293E-6_wp) &
                     * zqvsidiff * zeln8qsk
            sgdep = (0.0418521_wp-ppg*4.7524E-8_wp) &
                     * zqvsidiff *zeln6qgk
            ssdep = MAX(-zssmax ,ssdep )
            sgdep = MAX(-zsgmax ,sgdep )
          ENDIF !t_crit
        ENDIF !tg
      ENDIF

      !--------------------------------------------------------------------------
      ! Section 6: Search for grid points with rain in subsaturated areas
      !            and calculation of the evaporation rate of rain
      !--------------------------------------------------------------------------


#ifdef __COSMO__
        zqvsw    = fqvs( fpvsw(tg), ppg )
#endif

#ifdef __ICON__
        zqvsw    = sat_pres_water(tg)/(rhog * r_v *tg)
#endif

      IF( (llqr) .AND. (qvg+qcg <= zqvsw)) THEN

        zlnqrk   = LOG (zqrk)
        zx1      = 1.0_wp + zbev * EXP (zbevxp  * zlnqrk)
        !sev  = zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk)
        ! Limit evaporation rate in order to avoid overshoots towards supersaturation
        ! the pre-factor approximates (esat(T_wb)-e)/(esat(T)-e) at temperatures between 0 degC and 30 degC
        temp_c = tg - t0
        maxevap     = (0.61_wp-0.0163_wp*temp_c+1.111e-4_wp*temp_c**2)*(zqvsw-qvg)/zdt
        sev    = MIN(zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk), maxevap)

        IF( tg > zthn ) THEN
          ! Calculation of below-cloud rainwater freezing
          IF ( tg < ztrfrz ) THEN
            IF (lsuper_coolw) THEN
              !FR new: reduced rain freezing rate
              srfrz = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk
            ELSE
              srfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 ) * zeln27o16qrk
            ENDIF
          ENDIF
        ELSE ! Hom. freezing of rain water
          srfrz = zsrmax
        ENDIF
      ENDIF

      !--------------------------------------------------------------------------
      ! Section 7: Calculate the total tendencies of the prognostic variables.
      !            Update the prognostic variables in the interior domain.
      !--------------------------------------------------------------------------

      zsrsum = sev + srfrz + srcri
      zcorr  = 1.0_wp
      IF(zsrsum > 0) THEN
        zcorr  = zsrmax / MAX( zsrmax, zsrsum )
      ENDIF
      sev   = zcorr*sev
      srfrz = zcorr*srfrz
      srcri = zcorr*srcri
      
      zqvt =   sev    - sidep  - ssdep  - sgdep  - snuc   - sconr 
      zqct =   simelt - scau   - scfrz  - scac   - sshed  - srim   - srim2 
      zqit =   snuc   + scfrz  - simelt - sicri  + sidep  - sdau   - sagg   - sagg2  - siau
      zqrt =   scau   + sshed  + scac   + ssmelt + sgmelt - sev    - srcri  - srfrz  + sconr
      zqst =   siau   + sdau   - ssmelt + srim   + ssdep  + sagg   - sconsg
      zqgt =   sagg2  - sgmelt + sicri  + srcri  + sgdep  + srfrz  + srim2  + sconsg
      
      ztt = z_heat_cap_r*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst+zqgt) )

#ifdef __COSMO__
      IF (lsppt) THEN
!US how to check?        IF(ntstep>0.AND.i>=istart.and.i<=iend.and.j>=jstart.and.j<=jend) THEN
          CALL apply_tqx_tend_adj(itype_qxpert_rn, itype_qxlim_rn, ppg,     &
                                 tg,  qvg,  qcg,  qig,  qrg,  qsg,          &
                  pstoph(iv,k), ztt, zqvt, zqct, zqit, zqrt, zqst, ldum, qgg, zqgt)
!US        ENDIF
      ENDIF
#endif

      ! Update variables and add qi to qrs for water loading
      IF (lsedi_ice ) THEN
        qig = MAX ( 0.0_wp, (zzai*z1orhog + zqit*zdt)*zimi)
      ELSE
        qig = MAX ( 0.0_wp, qig + zqit*zdt)
      END IF
      qrg = MAX ( 0.0_wp, (zzar*z1orhog + zqrt*zdt)*zimr)
      qsg = MAX ( 0.0_wp, (zzas*z1orhog + zqst*zdt)*zims)
      qgg = MAX ( 0.0_wp, (zzag*z1orhog + zqgt*zdt)*zimg)

      
      !----------------------------------------------------------------------
      ! Section 10: Complete time step
      !----------------------------------------------------------------------

      IF ( k /= ke) THEN
        ! Store precipitation fluxes and sedimentation velocities 
        ! for the next level
        zprvr(iv) = qrg*rhog*zvzr(iv)
        zprvs(iv) = qsg*rhog*zvzs(iv)
        zprvg(iv) = qgg*rhog*zvzg(iv)
        zprvi(iv) = qig*rhog*zvzi(iv)
        IF (zprvr(iv) .LE. zqmin) zprvr(iv)=0.0_wp
        IF (zprvs(iv) .LE. zqmin) zprvs(iv)=0.0_wp
        IF (zprvg(iv) .LE. zqmin) zprvg(iv)=0.0_wp
        IF (zprvi(iv) .LE. zqmin) zprvi(iv)=0.0_wp


#ifdef NUDGING
        ! for the latent heat nudging
        IF ((llhn .OR. llhnverif) .AND. lhn_qrs ) THEN
          IF (lsedi_ice) THEN
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)+zprvg(iv)+zprvi(iv)
            qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv)+zpkg(iv)+zpki(iv))
          ELSE 
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)+zprvg(iv)
            qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv)+zpkg(iv))
          END IF
        ENDIF
#endif

        IF (qrg+qr(iv,k+1) <= zqmin) THEN
          zvzr(iv)= 0.0_wp
        ELSE
          zvzr(iv)= zvz0r * EXP(zvzxp*LOG((qrg+qr(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qsg+qs(iv,k+1) <= zqmin) THEN
          zvzs(iv)= 0.0_wp
        ELSE
          zvzs(iv)= zvz0s * EXP(ccswxp*LOG((qsg+qs(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qgg+qg(iv,k+1) <= zqmin ) THEN
          zvzg(iv)= 0.0_wp
        ELSE
          zvzg(iv)=zvz0g * EXP(zexpsedg*LOG((qgg+qg(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qig+qi(iv,k+1) <= zqmin ) THEN
          zvzi(iv)= 0.0_wp
        ELSE
          zvzi(iv)= zvz0i * EXP(zbvi*LOG((qig+qi(iv,k+1))*0.5_wp*rhog)) * zrho1o3
        ENDIF
          
      ELSE
        ! Precipitation fluxes at the ground
        prr_gsp(iv) = 0.5_wp * (qrg*rhog*zvzr(iv) + zpkr(iv))
        IF (lsedi_ice) THEN
          prs_gsp(iv) = 0.5_wp * (rhog*(qsg*zvzs(iv)+qig*zvzi(iv)) + zpks(iv)+zpki(iv))
        ELSE
          prs_gsp(iv) = 0.5_wp * (qsg*rhog*zvzs(iv) + zpks(iv))
        END IF
        prg_gsp(iv) = 0.5_wp * (qgg*rhog*zvzg(iv) + zpkg(iv))

#ifdef NUDGING
        ! for the latent heat nudging
        IF ((llhn .OR. llhnverif) .AND. lhn_qrs) THEN
          qrsflux(iv,k) = prr_gsp(iv)+prs_gsp(iv)+prg_gsp(iv)
        ENDIF
#endif

      ENDIF

      ! Update of prognostic variables or tendencies
      qr (iv,k) = MAX ( 0.0_wp, qrg )
      qs (iv,k) = MAX ( 0.0_wp, qsg )
      qi (iv,k) = MAX ( 0.0_wp, qig )
      qg (iv,k) = MAX ( 0.0_wp, qgg )
      t  (iv,k) = t (iv,k) + ztt*zdt 
      qv (iv,k) = MAX ( 0.0_wp, qv(iv,k) + zqvt*zdt )
      qc (iv,k) = MAX ( 0.0_wp, qc(iv,k) + zqct*zdt )

      IF (izdebug > 15) THEN
        ! Check for negative values
        IF (qr(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: graupel, negative value in qr'
          CALL message('',message_text)
        ENDIF
        IF (qc(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: graupel, negative value in qc'
          CALL message('',message_text)
        ENDIF
        IF (qi(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: graupel, negative value in qi'
          CALL message('',message_text)
        ENDIF
        IF (qs(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: graupel, negative value in qs'
          CALL message('',message_text)
        ENDIF
        IF (qv(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: graupel, negative value in qv'
          CALL message('',message_text)
        ENDIF
        IF (qg(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: graupel, negative value in qg'
          CALL message('',message_text)
        ENDIF
      ENDIF

    ENDDO  !loop over iv
  END DO loop_over_levels

#if defined (__COSMO__)
  !XL : currently satad is outside of the k loop as function call
  !     is not yet supported (without inlinining) on GPU
  !     Once it is supported satad will go back in the k,i loop
  ! Do a final saturation adjustment for new values of t, qv and qc
  DO  k = k_start, ke
    CALL satad ( 1, t(:,k), qv(:,k),                  &
               qc(:,k), t(:,k), p(:,k),               &
               zdummy(:,1),zdummy(:,2),zdummy(:,3),   &
               zdummy(:,4),zdummy(:,5),zdummy(:,6),   &
               zdummy(:,7),zdummy(:,8),               &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv, &
               rvd_m_o, lh_v, z_heat_cap_r, cp_d,     &
               nvec, 1, iv_start, iv_end, 1 , 1 )
  ENDDO

  DO  k = k_start, ke
    DO iv=iv_start,iv_end
      IF ( ldiabf_lh ) THEN
        ! compute temperature increment due to latent heat
        tinc_lh(iv,k) = tinc_lh(iv,k) + t(iv,k)
      ENDIF

#ifdef NUDGING
      ! add part of latent heating calculated in subroutine hydci to model latent
      ! heating field: add temperature to model latent heating field
      IF (llhn .OR. llhnverif) THEN
        !CALL get_gs_lheating ('inc',1,ke)  !XL: this should be called from within the block
        tt_lheat(iv,k) = tt_lheat(iv,k) + t(iv,k)
      ENDIF
#endif
    ENDDO
  ENDDO
#endif

!------------------------------------------------------------------------------
! final tendency calculation for ICON
!
! Note: as soon as we have a new satad subroutine in ICON, this tendency
! calculation will be done in the k-loop and the original 3D variables wont
! be used to store the new values. Then we wont need the _in variables anymore.
!------------------------------------------------------------------------------

  IF (PRESENT(ddt_tend_t)) THEN

    DO k=k_start,ke
       DO iv=iv_start,iv_end

          ! calculated pseudo-tendencies
          ddt_tend_t (iv,k) = (t (iv,k) - t_in (iv,k))*zdtr
          ddt_tend_qv(iv,k) = MAX(-qv_in(iv,k)*zdtr,(qv(iv,k) - qv_in(iv,k))*zdtr)
          ddt_tend_qc(iv,k) = MAX(-qc_in(iv,k)*zdtr,(qc(iv,k) - qc_in(iv,k))*zdtr)
          ddt_tend_qr(iv,k) = MAX(-qr_in(iv,k)*zdtr,(qr(iv,k) - qr_in(iv,k))*zdtr)
          ddt_tend_qs(iv,k) = MAX(-qs_in(iv,k)*zdtr,(qs(iv,k) - qs_in(iv,k))*zdtr)
          ddt_tend_qi(iv,k) = MAX(-qi_in(iv,k)*zdtr,(qi(iv,k) - qi_in(iv,k))*zdtr)
          ddt_tend_qg(iv,k) = MAX(-qg_in(iv,k)*zdtr,(qg(iv,k) - qg_in(iv,k))*zdtr)

      END DO
    END DO

  END IF

  IF (izdebug > 15) THEN
   CALL message('gscp_graupel', 'UPDATED VARIABLES')
   WRITE(message_text,'(A,2E20.9)') 'graupel  T= ',&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'graupel qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'graupel qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'graupel qi= ',&
    MAXVAL( qi(:,:)), MINVAL(qi(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'graupel qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'graupel qs= ',&
    MAXVAL( qs(:,:)), MINVAL(qs(:,:) )
   CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'graupel qg= ',&
    MAXVAL( qg(:,:)), MINVAL(qg(:,:) )
    CALL message('', TRIM(message_text))
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine graupel
!------------------------------------------------------------------------------

END SUBROUTINE graupel

!==============================================================================

END MODULE gscp_graupel
