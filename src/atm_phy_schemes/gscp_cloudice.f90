!>
!! cloud microphysics
!!
!! !------------------------------------------------------------------------------
!!
!! @par Description of *gscp_cloudice*:
!!   This module procedure calculates the rates of change of temperature, cloud
!!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
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
!! @par reference   This is an adaption of subroutine cloudice in file src_gscp.f90
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

MODULE gscp_cloudice

!------------------------------------------------------------------------------
!>
!! Description:
!!
!!   The subroutine in this module calculates the rates of change of
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
                                 g     => grav  , & !! acceleration due to gravity
                                 t0    => tmelt     !! melting temperature of ice/snow

USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b2i   => c3ies , & !!               -- " --
                                 b4w   => c4les , & !!               -- " --
                                 b4i   => c4ies , & !!               -- " --
                                 b234w => c5les     !!               -- " --
USE mo_satad,              ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
                                 sat_pres_ice!,   &  !! saturation vapor pressure w.r.t. ice
USE mo_exception,          ONLY: message, message_text
USE mo_run_config,         ONLY: ldass_lhn
#endif

!------------------------------------------------------------------------------

! this can be used by ICON and COSMO
USE gscp_data, ONLY: &          ! all variables are used here

    ccsrim,    ccsagg,    ccsdep,    ccsvel,    ccsvxp,    ccslam,       &
    ccslxp,    ccsaxp,    ccsdxp,    ccshi1,    ccdvtp,    ccidep,       &
    ccswxp,    zconst,    zcev,      zbev,      zcevxp,    zbevxp,       &
    zvzxp,     zvz0r,                                                    &

    v0snow,    mu_rain,   rain_n0_factor,       cloud_num,               &

    x13o8,     x1o2,      x27o16,    x3o4,      x7o4,      x7o8,         &
    zbms,      zbvi,      zcac,      zccau,     zciau,     zcicri,       &
    zcrcri,    zcrfrz,    zcrfrz1,   zcrfrz2,   zeps,      zkcac,        &
    zkphi1,    zkphi2,    zkphi3,    zmi0,      zmimax,    zmsmin,       &
    zn0s0,     zn0s1,     zn0s2,     znimax_thom,          zqmin,        &
    zrho0,     zthet,     zthn,      ztmix,     ztrfrz,    zv1s,         &
    zvz0i,     x13o12,    x2o3,      x5o24,     zams,      zasmel,       &
    zbsmel,    zcsmel,    icesedi_exp,                                   &
    iautocon,  isnow_n0temp, dist_cldtop_ref,   reduce_dep_ref,          &
    tmin_iceautoconv,     zceff_fac, zceff_min,                          &
    mma, mmb

#ifdef __ICON__
! this is (at the moment) an ICON part
USE gscp_data, ONLY: &          ! all variables are used here
    vtxexp,    & !  kc_c1,     & !
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

PUBLIC :: cloudice

LOGICAL, PARAMETER :: &
#ifdef __COSMO__
  lorig_icon   = .FALSE. , &  ! switch for original ICON setup (only for cloudice)
                              ! XL : should be false for COSMO ?
  lred_depgrowth = .TRUE. , & ! switch for reduced depositional growth near tops of stratus clouds
#else
  lorig_icon   = .TRUE.  , &  ! switch for original ICON setup (only for cloudice)
  lred_depgrowth = .FALSE., & ! switch for reduced depositional growth near tops of stratus clouds
                              ! (not used in ICON because it degrades scores at high latitudes)
#endif

  lsedi_ice    = .FALSE. , &  ! switch for sedimentation of cloud ice (Heymsfield & Donner 1990 *1/3)
  lstickeff    = .FALSE. , &  ! switch for sticking coeff. (work from Guenther Zaengl)
  lsuper_coolw = .TRUE.       ! switch for supercooled liquid water (work from Felix Rieper)

!------------------------------------------------------------------------------
!> Parameters and variables which are global in this module
!------------------------------------------------------------------------------

#ifdef __COSMO__
CHARACTER(132) :: message_text = ''
#endif

!==============================================================================

CONTAINS

!==============================================================================
!> Module procedure "cloudice" in "gscp_cloudice" for computing effects of 
!!  grid scale precipitation including cloud water, cloud ice, rain and snow
!------------------------------------------------------------------------------

SUBROUTINE cloudice (             &
  nvec,ke,                           & !> array dimensions
  ivstart,ivend, kstart,             & !! optional start/end indicies
  idbg,                              & !! optional debug level
  zdt, dz,                           & !! numerics parameters
  t,p,rho,qv,qc,qi,qr,qs,qnc,        & !! prognostic variables
#ifdef __ICON__
  !xxx: this should become a module variable, e.g. in a new module mo_gscp_data.f90
  qi0,qc0,                           & !! cloud ice/water threshold for autoconversion
#endif
  prr_gsp,prs_gsp,                   & !! surface precipitation rates
#ifdef __COSMO__
  tinc_lh,                           & !  t-increment due to latent heat 
  pstoph,                            & !  stochastic multiplier of physics tendencies
#endif
  qrsflux,                           & !  total precipitation flux
  l_cv,                              &
  ldiag_ttend,     ldiag_qtend     , &
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
  ddt_diag_shed                      ) !!

!------------------------------------------------------------------------------
!>
!! Description:
!!   This module procedure calculates the rates of change of temperature, cloud
!!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
!!   processes related to the formation of grid scale precipitation. The
!!   variables are updated in this subroutine. Rain and snow are prognostic
!!   variables. The precipitation fluxes at the surface are stored on the
!!   corresponding global fields.
!!
!! Method:
!!
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

  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv                   !! if true, cv is used instead of cp

  LOGICAL, INTENT(IN), OPTIONAL :: &
    ldiag_ttend,         & ! if true, temperature tendency shall be diagnosed
    ldiag_qtend            ! if true, moisture tendencies shall be diagnosed

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    qr              ,    & !! specific rain content                         (kg/kg)
    qs                     !! specific snow content                         (kg/kg)

#ifdef __COSMO__
  REAL(KIND=wp), INTENT(INOUT) :: &
       tinc_lh(:,:)    ! temperature increments due to heating             ( K/s ) 

  REAL(KIND=wp), INTENT(IN)    :: &
       pstoph(:,:)     ! stochastic multiplier of physics tendencies
#endif  

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
       qrsflux        ! total precipitation flux (nudg)

  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    prs_gsp,             & !! precipitation rate of snow, grid-scale        (kg/(m2*s))
    qnc                    !! cloud number concentration

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT), OPTIONAL ::   &     ! dim (ie,ke)
    ddt_tend_t      , & !> tendency T                                       ( 1/s )
    ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
    ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
    ddt_tend_qi     , & !! tendency qi                                      ( 1/s )
    ddt_tend_qr     , & !! tendency qr                                      ( 1/s )
    ddt_tend_qs         !! tendency qs                                      ( 1/s )

  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &   ! dim (ie,ke)
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
  
  INTEGER (KIND=i4) :: &
    iv, k             !> loop indices

  REAL    (KIND=wp   ) :: nnr

  REAL    (KIND=wp   ) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)

  INTEGER ::  &
    iv_start     ,    & !> start index for horizontal direction
    iv_end       ,    & !! end index for horizontal direction
    k_start      ,    & !! model level where computations start
    izdebug             !! debug level

  REAL    (KIND=wp   ) ::  &
    fpvsw, fpvsi, fqvs,& ! name of statement functions
    fxna ,             & ! statement function for ice crystal number
    fxna_cooper ,      & ! statement function for ice crystal number, Cooper(1986) 
    ztx  , zpv  , zpx ,& ! dummy arguments for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    znimix,            & ! number of ice crystals at ztmix -> threshold temp for mixed-phase clouds
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
    zdtr ,             & ! reciprocal of timestep for integration
    zscau , zscac  , zscrim , zscshe,         & ! local values of the
    zsiau , zsagg  , zsidep , zsicri, zsrcri, & ! transfer rates
    zsdau , zssdep , zssmelt, & ! defined below
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum, zsrmax,    & ! terms for limiting  total rain water depletion
    zssmax,            & ! term for limiting snow depletion
    znin,              & ! number of cloud ice crystals at nucleation
    fnuc,              & !FR: coefficient needed for Forbes (2012) SLW layer parameterization 
    znid,              & ! number of cloud ice crystals for deposition
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,   & ! terms for limiting total
    zztau, zxfac, zx1, zx2, ztt, &   ! some help variables
    ztau, zphi, zhi, zdvtp, ztc, zlog_10

  REAL    (KIND=wp   ) ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice
    zqrt   ,& ! layer tendency of rain
    zqst      ! layer tendency of snow

  REAL    (KIND=wp   ) ::  &
    zlnqrk,zlnqsk,zlnqik,     & !
    zlnlogmi,ccswxp_ln1o2,zvzxp_ln1o2,zbvi_ln1o2,  & !
    qcg,tg,qvg,qrg, qsg,qig,rhog,ppg, alf,bet,m2s,m3s,hlp,maxevap,temp_c,stickeff

  LOGICAL :: &
    llqr,llqs,llqc,llqi  !   switch for existence of qr, qs, qc, qi

  LOGICAL :: lldiag_ttend, lldiag_qtend

  REAL(KIND=wp), DIMENSION(nvec,ke) ::   &
    t_in               ,    & !> temperature                                   (  K  )
    qv_in              ,    & !! specific water vapor content                  (kg/kg)
    qc_in              ,    & !! specific cloud water content                  (kg/kg)
    qi_in              ,    & !! specific cloud ice   content                  (kg/kg)
    qr_in              ,    & !! specific rain content                         (kg/kg)
    qs_in                     !! specific snow content                         (kg/kg)


!! Local (automatic) arrays:
!! -------------------------

  REAL    (KIND=wp   ) ::  &
    zqvsi       (nvec),     & !> sat. specitic humidity at ice saturation
    zqvsw       (nvec),     & !> sat. specitic humidity at water saturation
    zqvsw_up    (nvec),     & !> sat. specitic humidity at water saturation in previous layer
    zvzr        (nvec),     & !
    zvzs        (nvec),     & !
    zvzi        (nvec),     & ! terminal fall velocity of ice
    zpkr        (nvec),     & !
    zpks        (nvec),     & !
    zpki        (nvec),     & ! precipitation flux of ice
    zprvr       (nvec),     & !
    zprvs       (nvec),     & !
    zprvi       (nvec),     & !
#ifdef __COSMO__
    zdummy      (nvec,8),   & !
#endif
    zcsdep      (nvec),     & !
    zcidep      (nvec),     & !
    zvz0s       (nvec),     & !
    zcrim       (nvec),     & !
    zcagg       (nvec),     & !
    zbsdep      (nvec),     & !
    zcslam      (nvec),     & !
    zn0s        (nvec),     & !
    zimr        (nvec),     & !
    zims        (nvec),     & !
    zimi        (nvec),     & !
    zzar        (nvec),     & !
    zzas        (nvec),     & !
    zzai        (nvec),     & !
    zqrk        (nvec),     & !
    zqsk        (nvec),     & !
    zqik        (nvec),     & !
    zdtdh       (nvec),     & !
    z1orhog     (nvec),     & ! 1/rhog
    zrho1o2     (nvec),     & ! (rho0/rhog)**1/2
    zrhofac_qi  (nvec),     & ! (rho0/rhog)**icesedi_exp
    zeln7o8qrk  (nvec),     & !
    zeln7o4qrk  (nvec),     & ! FR new     
    zeln27o16qrk(nvec),     & !
    zeln13o8qrk (nvec),     & !
!    zeln3o16qrk (nvec),     & !
    zeln5o24qsk (nvec),     & !
    zeln2o3qsk  (nvec)        !

  REAL    (KIND=wp   ) ::  &
    scau   (nvec), & ! transfer rate due to autoconversion of cloud water
    scac   (nvec), & ! transfer rate due to accretion of cloud water
    snuc   (nvec), & ! transfer rate due nucleation of cloud ice
    scfrz  (nvec), & ! transfer rate due homogeneous freezing of cloud water
    simelt (nvec), & ! transfer rate due melting of cloud ice
    sidep  (nvec), & ! transfer rate due depositional growth of cloud ice
    ssdep  (nvec), & ! transfer rate due depositional growth of snow
    sdau   (nvec), & ! transfer rate due depositional cloud ice autoconversion
    srim   (nvec), & ! transfer rate due riming of snow
    sshed  (nvec), & ! transfer rate due shedding
    sicri  (nvec), & ! transfer rate due cloud ice collection by rain (sink qi)
    srcri  (nvec), & ! transfer rate due cloud ice collection by rain (sink qr)
    sagg   (nvec), & ! transfer rate due aggregation of snow and cloud ice
    siau   (nvec), & ! transfer rate due autoconversion of cloud ice
    ssmelt (nvec), & ! transfer rate due melting of snow
    sev    (nvec), & ! transfer rate due evaporation of rain
    srfrz  (nvec), & ! transfer rate due to rainwater freezing
    reduce_dep(nvec),&!FR: coefficient: reduce deposition at cloud top (Forbes 2012)
    dist_cldtop(nvec) !FR: distance from cloud top layer

 
  ! Dimensions and loop counter for storing the indices
  INTEGER (KIND=i4)        ::  &
    ic1, ic2, ic3, ic4, ic5, ic6, i1d

  !> Integer arrays for a better vectorization
  INTEGER (KIND=i4)        ::  &
    ivdx1(nvec), & !!
    ivdx2(nvec), & !!
    ivdx3(nvec), & !!
    ivdx4(nvec), & !!
    ivdx5(nvec), & !!
    ivdx6(nvec)    !!

  LOGICAL :: ldum

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine cloudice
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
    znimax = znimax_Thom         !znimax_Thom = 250.E+3_wp,
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSEIF(lorig_icon) THEN
    znimax = 150.E+3_wp     ! from previous ICON code 
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSE
    znimax = fxna(zthn) ! Maximum number of cloud ice crystals
    znimix = fxna(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  END IF

  zpvsw0 = fpvsw(t0)  ! sat. vap. pressure for t = t0
  zlog_10 = LOG(10._wp) ! logarithm of 10

  ! Precomputations for optimization
  ccswxp_ln1o2 = EXP (ccswxp * LOG (0.5_wp))
  zvzxp_ln1o2  = EXP (zvzxp * LOG (0.5_wp))
  zbvi_ln1o2   = EXP (zbvi * LOG (0.5_wp))

! Delete precipitation fluxes from previous timestep
!CDIR BEGIN COLLAPSE
    prr_gsp (:) = 0.0_wp
    prs_gsp (:) = 0.0_wp
    zpkr    (:) = 0.0_wp
    zpks    (:) = 0.0_wp
    zpki    (:) = 0.0_wp
    zprvr   (:) = 0.0_wp
    zprvs   (:) = 0.0_wp
    zprvi   (:) = 0.0_wp
    zvzr    (:) = 0.0_wp
    zvzs    (:) = 0.0_wp
    zvzi    (:) = 0.0_wp
    zqvsw   (:) = 0.0_wp
    dist_cldtop(:) = 0.0_wp    
!CDIR END


! Optional arguments

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
  IF (PRESENT(ldiag_ttend)) THEN
    lldiag_ttend = ldiag_ttend
  ELSE
    lldiag_ttend = .FALSE.
  ENDIF
  IF (PRESENT(ldiag_qtend)) THEN
    lldiag_qtend = ldiag_qtend
  ELSE
    lldiag_qtend = .FALSE.
  ENDIF

  ! save input arrays for final tendency calculation
  IF (lldiag_ttend) THEN
    t_in  = t
  ENDIF
  IF (lldiag_qtend) THEN
    qv_in = qv
    qc_in = qc
    qi_in = qi
    qr_in = qr
    qs_in = qs
  END IF

! timestep for calculations
  zdtr  = 1.0_wp / zdt

  ! add part of latent heating calculated in subroutine cloudice to model latent
  ! heating field: subtract temperature from model latent heating field
  IF (ldass_lhn) THEN
!CDIR COLLAPSE
      qrsflux(:,:) = 0.0_wp
  ENDIF


! output for various debug levels
  IF (izdebug > 15) CALL message('gscp_cloudice: ','Start of cloudice')
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
  ENDIF

!CDIR COLLAPSE


! *********************************************************************
! Loop from the top of the model domain to the surface to calculate the
! transfer rates  and sedimentation terms
! *********************************************************************

  loop_over_levels: DO  k = k_start, ke


#ifdef __COSMO__
    IF ( ldiabf_lh ) THEN
      ! initialize temperature increment due to latent heat
      tinc_lh(:,k) = tinc_lh(:,k) - t(:,k)
    ENDIF
#endif

  !----------------------------------------------------------------------------
  ! Section 2: Check for existence of rain and snow
  !            Initialize microphysics and sedimentation scheme
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
    zcrim (:) = 0.0_wp
    zcagg (:) = 0.0_wp
    zbsdep(:) = 0.0_wp
    zvz0s (:) = 0.0_wp
    zn0s  (:) = zn0s0
    reduce_dep(:) = 1.0_wp  !FR: Reduction coeff. for dep. growth of rain and ice  
!CDIR END

    ic1 = 0
    ic2 = 0
    ic3 = 0

    DO iv = iv_start, iv_end

        qrg = qr(iv,k)
        qsg = qs(iv,k)
        qig = qi(iv,k)
        rhog = rho(iv,k)

        !..for density correction of fall speeds
        z1orhog(iv) = 1.0_wp/rhog
        hlp         = LOG(zrho0*z1orhog(iv))
        zrho1o2(iv) = EXP(hlp*x1o2) ! exponent 0.5 for rain and snow
        zrhofac_qi(iv) = EXP(hlp*icesedi_exp) ! user-defined exponent for cloud ice (default 0.33)

        zqrk(iv) = qrg * rhog
        zqsk(iv) = qsg * rhog
        zqik(iv) = qig * rhog

        zdtdh(iv) = 0.5_wp * zdt / dz(iv,k)

        zzar(iv)   = zqrk(iv)/zdtdh(iv) + zprvr(iv) + zpkr(iv)
        zzas(iv)   = zqsk(iv)/zdtdh(iv) + zprvs(iv) + zpks(iv)
        zzai(iv)   = zqik(iv)/zdtdh(iv) + zprvi(iv) + zpki(iv)

     ENDDO

    DO iv = iv_start, iv_end

        llqr = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin
        llqi = zqik(iv) > zqmin

        IF (llqs) THEN
          ic1 = ic1 + 1
          ivdx1(ic1) = iv
        ELSE
          zpks(iv) = 0._wp
        ENDIF
        IF (llqr) THEN
          ic2 = ic2 + 1
          ivdx2(ic2) = iv
        ELSE
          zpkr(iv) = 0._wp
        ENDIF
        IF (llqi) THEN
          ic3 = ic3 + 1
          ivdx3(ic3) = iv
        ELSE
          zpki(iv) = 0._wp
        ENDIF
    
     ENDDO

!DIR$ IVDEP
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs_prepare: DO i1d = 1, ic1
      iv = ivdx1(i1d)

      qsg = qs(iv,k)
      tg  = t(iv,k)

      IF (isnow_n0temp == 1) THEN
        ! Calculate n0s using the temperature-dependent
        ! formula of Field et al. (2005)
        ztc = tg - t0
        ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
        zn0s(iv) = zn0s1*EXP(zn0s2*ztc)
        zn0s(iv) = MIN(zn0s(iv),1e9_wp)
        zn0s(iv) = MAX(zn0s(iv),1e6_wp)
      ELSEIF (isnow_n0temp == 2) THEN
        ! Calculate n0s using the temperature-dependent moment
        ! relations of Field et al. (2005)
        ztc = tg - t0
        ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)

        nnr  = 3._wp
        hlp = mma(1) + mma(2)*ztc + mma(3)*nnr + mma(4)*ztc*nnr &
          & + mma(5)*ztc**2 + mma(6)*nnr**2 + mma(7)*ztc**2*nnr &
          & + mma(8)*ztc*nnr**2 + mma(9)*ztc**3 + mma(10)*nnr**3
        alf = EXP(hlp*zlog_10) ! 10.0_wp**hlp
        bet = mmb(1) + mmb(2)*ztc + mmb(3)*nnr + mmb(4)*ztc*nnr &
          & + mmb(5)*ztc**2 + mmb(6)*nnr**2 + mmb(7)*ztc**2*nnr &
          & + mmb(8)*ztc*nnr**2 + mmb(9)*ztc**3 + mmb(10)*nnr**3

        ! Here is the exponent bms=2.0 hardwired! not ideal! (Uli Blahak)
        m2s = qsg * rho(iv,k) / zams   ! UB rho added as bugfix
        m3s = alf*EXP(bet*LOG(m2s))

        hlp  = zn0s1*EXP(zn0s2*ztc)
        zn0s(iv) = 13.50_wp * m2s * (m2s / m3s) **3
        zn0s(iv) = MAX(zn0s(iv),0.5_wp*hlp)
        zn0s(iv) = MIN(zn0s(iv),1e2_wp*hlp)
        zn0s(iv) = MIN(zn0s(iv),1e9_wp)
        zn0s(iv) = MAX(zn0s(iv),1e6_wp)
      ELSE
        ! Old constant n0s
        zn0s(iv) = 8.0e5_wp
      ENDIF
      zcrim (iv) = ccsrim*zn0s(iv)
      zcagg (iv) = ccsagg*zn0s(iv)
      zbsdep(iv) = ccsdep*SQRT(v0snow)
      zvz0s (iv) = ccsvel*EXP(ccsvxp * LOG(zn0s(iv)))

      zlnqsk = zvz0s(iv) * EXP (ccswxp * LOG (zqsk(iv))) * zrho1o2(iv)
      zpks(iv) = zqsk (iv) * zlnqsk

      IF (zvzs(iv) == 0.0_wp) THEN
        zvzs(iv) = zlnqsk * ccswxp_ln1o2
      ENDIF
    ENDDO loop_over_qs_prepare

!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
    loop_over_qr_sedi: DO i1d = 1, ic2
      iv = ivdx2(i1d)

      zlnqrk = zvz0r * EXP (zvzxp * LOG (zqrk(iv))) * zrho1o2(iv)
      zpkr(iv) = zqrk(iv) * zlnqrk

      IF (zvzr(iv) == 0.0_wp) THEN
        zvzr(iv) = zlnqrk * zvzxp_ln1o2
      ENDIF
    ENDDO loop_over_qr_sedi

!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
    loop_over_qi_sedi: DO i1d = 1, ic3
      iv = ivdx3(i1d)
      
      zlnqik = zvz0i * EXP (zbvi * LOG (zqik(iv))) * zrhofac_qi(iv)
      zpki(iv) = zqik (iv) * zlnqik

      IF (zvzi(iv) == 0.0_wp) THEN
        zvzi(iv) = zlnqik * zbvi_ln1o2
      ENDIF
    ENDDO loop_over_qi_sedi


  !----------------------------------------------------------------------------
  ! Section 3:
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
    zeln7o8qrk   (:) = 0.0_wp
    zeln7o4qrk   (:) = 0.0_wp !FR
    zeln27o16qrk (:) = 0.0_wp
    zeln13o8qrk  (:) = 0.0_wp
    zeln5o24qsk  (:) = 0.0_wp
    zeln2o3qsk   (:) = 0.0_wp

!FR old  
!   zcsdep       (:) = 3.2E-2_wp
    zcsdep       (:) = 3.367E-2_wp    
    zcidep       (:) = 1.3E-5_wp
    zcslam       (:) = 1e10_wp

    scau         (:) = 0.0_wp
    scac         (:) = 0.0_wp
    snuc         (:) = 0.0_wp
    scfrz        (:) = 0.0_wp
    simelt       (:) = 0.0_wp
    sidep        (:) = 0.0_wp
    ssdep        (:) = 0.0_wp
    sdau         (:) = 0.0_wp
    srim         (:) = 0.0_wp
    sshed        (:) = 0.0_wp
    sicri        (:) = 0.0_wp
    srcri        (:) = 0.0_wp
    sagg         (:) = 0.0_wp
    siau         (:) = 0.0_wp
    ssmelt       (:) = 0.0_wp
    sev          (:) = 0.0_wp
    srfrz        (:) = 0.0_wp
!CDIR END

    ic1 = 0
    ic2 = 0
    ic3 = 0
    ic4 = 0
    ic5 = 0
    ic6 = 0

    DO iv = iv_start, iv_end

        tg   =  t(iv,k)
        rhog = rho(iv,k)

        zqvsw_up(iv) = zqvsw(iv)

#ifdef __COSMO__
        ppg  =  p(iv,k)
        zqvsi(iv) = fqvs( fpvsi(tg), ppg )
        zqvsw(iv) = fqvs( fpvsw(tg), ppg )
#endif

#ifdef __ICON__
        hlp = 1._wp/(rhog * r_v * tg)
        zqvsi(iv) = sat_pres_ice(tg) * hlp
        zqvsw(iv) = sat_pres_water(tg) * hlp
#endif


        zpkr(iv)   = MIN( zpkr(iv) , zzar(iv) )
        zpks(iv)   = MIN( zpks(iv) , zzas(iv) )
        zpki(iv)   = MIN( zpki(iv) , zzai(iv) )

        zzar(iv)   = zdtdh(iv) * (zzar(iv)-zpkr(iv))
        zzas(iv)   = zdtdh(iv) * (zzas(iv)-zpks(iv))
        zzai(iv)   = zdtdh(iv) * (zzai(iv)-zpki(iv))

        zimr(iv)   = 1.0_wp / (1.0_wp + zvzr(iv) * zdtdh(iv))
        zims(iv)   = 1.0_wp / (1.0_wp + zvzs(iv) * zdtdh(iv))
        zimi(iv)   = 1.0_wp / (1.0_wp + zvzi(iv) * zdtdh(iv))

        zqrk(iv)   = zzar(iv)*zimr(iv)
        zqsk(iv)   = zzas(iv)*zims(iv)
        zqik(iv)   = zzai(iv)*zimi(iv)

     ENDDO

    DO iv = iv_start, iv_end

        qvg  = qv(iv,k)
        qcg  = qc(iv,k)
        qig  = qi(iv,k)
        tg   =  t(iv,k)

        llqr = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin
        llqi = zqik(iv) > zqmin
        llqc =      qcg > zqmin

        IF (llqr) THEN
          ic1 = ic1 + 1
          ivdx1(ic1) = iv
        ENDIF
        IF (llqs) THEN
          ic2 = ic2 + 1
          ivdx2(ic2) = iv
        ENDIF
        IF (llqi .OR. llqs) THEN
          ic3 = ic3 + 1
          ivdx3(ic3) = iv
        ENDIF
        IF ( tg < zthet .AND. qvg >  8.E-6_wp &
                        .AND. qig <= 0.0_wp ) THEN
          ic4 = ic4 + 1
          ivdx4(ic4) = iv
        ENDIF
        IF (llqc) THEN
          ic5 = ic5 + 1
          ivdx5(ic5) = iv
        ENDIF
        ! Use qvg+qcg <= zqvsw(iv) instead of qcg == 0 in order to be independent of satad
        ! between turbulence and microphysics
        IF (llqr .AND. qvg+qcg <= zqvsw(iv)) THEN
          ic6 = ic6 + 1
          ivdx6(ic6) = iv
        ENDIF
    
     ENDDO


! ic1
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qr: DO i1d =1, ic1
      iv = ivdx1(i1d)

      qcg  = qc(iv,k)
      qig  = qi(iv,k)
      tg   =  t(iv,k)
      llqi =  qig > zqmin

      zlnqrk       = LOG (zqrk(iv))
      IF ( qig+qcg > zqmin ) THEN
        zeln7o8qrk(iv)   = EXP (x7o8   * zlnqrk)
      ENDIF
      IF ( tg < ztrfrz ) THEN
        zeln7o4qrk(iv)   = EXP (x7o4   * zlnqrk) !FR new
        zeln27o16qrk(iv) = EXP (x27o16 * zlnqrk)
      ENDIF
      IF (llqi) THEN
        zeln13o8qrk(iv)  = EXP (x13o8  * zlnqrk)
      ENDIF
!      IF (qcg <= 0.0_wp ) THEN
!        zeln3o16qrk(iv)  = EXP (x3o16  * zlnqrk)
!      ENDIF
    ENDDO loop_over_qr


! ic2
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
    loop_over_qs_coeffs: DO i1d =1, ic2
      iv = ivdx2(i1d)

      qcg = qc(iv,k)
      qig = qi(iv,k)

      zlnqsk       = LOG (zqsk(iv))
      zeln5o24qsk(iv)  = EXP (x5o24  * zlnqsk)
      zeln2o3qsk(iv)   = EXP (x2o3   * zlnqsk)
    ENDDO loop_over_qs_coeffs

! ic3
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qi_qs: DO i1d =1, ic3
      iv = ivdx3(i1d)

      tg   =   t(iv,k)
      ppg  =   p(iv,k)
      rhog = rho(iv,k)
      llqs = zqsk(iv) > zqmin
      zdvtp  = ccdvtp * EXP(1.94_wp * LOG(tg)) / ppg          
      zhi    = ccshi1*zdvtp*rhog*zqvsi(iv)/(tg*tg)
      hlp    = zdvtp / (1.0_wp + zhi)
      zcidep(iv) = ccidep * hlp
      IF (llqs) THEN
        zcslam(iv) = EXP(ccslxp * LOG(ccslam * zn0s(iv) / zqsk(iv) ))
        zcslam(iv) = MIN(zcslam(iv),1e15_wp)
        zcsdep(iv) = 4.0_wp * zn0s(iv) * hlp
      ENDIF


    ENDDO loop_over_qi_qs

    !--------------------------------------------------------------------------
    ! Section 4: Initialize the conversion rates with zeros in every layer
    !            Deposition nucleation for low temperatures
    !--------------------------------------------------------------------------

    ! Deposition nucleation for low temperatures below a threshold

! ic4
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_icenucleation: DO i1d = 1, ic4
      iv = ivdx4(i1d)

      qvg  =  qv(iv,k)
      tg   =   t(iv,k)
      rhog = rho(iv,k)

      IF( qvg > zqvsi(iv) ) THEN
        IF( lsuper_coolw .OR. lorig_icon) THEN
          znin  = MIN( fxna_cooper(tg), znimax )
        ELSE
          znin  = MIN( fxna(tg), znimax )
        END IF
        snuc(iv) = zmi0 * z1orhog(iv) * znin * zdtr
      ENDIF

    ENDDO loop_over_icenucleation

    !--------------------------------------------------------------------------
    ! Section 5: Search for cloudy grid points with cloud water and
    !            calculation of the conversion rates involving qc
    !--------------------------------------------------------------------------

! ic5
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qc: DO i1d =1, ic5
      iv = ivdx5(i1d)

      qrg  =   qr(iv,k)
      qvg  =   qv(iv,k)
      qcg  =   qc(iv,k)
      qig  =   qi(iv,k)
      tg   =    t(iv,k)
      ppg  =    p(iv,k)
      rhog =  rho(iv,k)
      llqs = zqsk(iv) > zqmin

      IF (iautocon == 0) THEN
        ! Kessler (1969) autoconversion rate
        zscau  = zccau * MAX( qcg - qc0, 0.0_wp )
        zscac  = zcac  * qcg * zeln7o8qrk(iv)
      ELSEIF (iautocon == 1) THEN
        ! Seifert and Beheng (2001) autoconversion rate
        ! with cloud droplet number concentration qnc depending on aerosol climatology
        IF (qcg > 1.e-6_wp) THEN
          ztau   = MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp)
          ztau   = MAX(ztau,1.E-30_wp)
          hlp    = EXP(zkphi2*LOG(ztau))
          zphi   = zkphi1 * hlp * (1.0_wp - hlp)**3
          zscau  = zconst * qcg*qcg*qcg*qcg/(qnc(iv)*qnc(iv)) &
            &    * (1.0_wp + zphi/(1.0_wp - ztau)**2)
          zphi   = (ztau/(ztau+zkphi3))**4
          zscac  = zkcac * qcg * qrg * zphi !* zrho1o2(iv)
        ELSE
          zscau  = 0.0_wp
          zscac  = 0.0_wp
        ENDIF
      ENDIF
      IF (llqs) THEN
        zscrim = zcrim(iv) * EXP(ccsaxp * LOG(zcslam(iv))) * qcg !* zrho1o2(iv)
      ELSE
        zscrim = 0.0_wp
      ENDIF

      zscshe = 0.0_wp
      IF( tg >= t0 ) THEN
        zscshe = zscrim
        zscrim = 0.0_wp
      ENDIF
      ! Check for maximum depletion of cloud water and adjust the
      ! transfer rates accordingly
      zscmax = qcg*zdtr 
      zscsum = zscau + zscac + zscrim + zscshe 
      zcorr  = zscmax / MAX( zscmax, zscsum )
      IF( tg <= zthn ) THEN
        scfrz(iv) = zscmax
      ELSE
        scau (iv) = zcorr*zscau
        scac (iv) = zcorr*zscac
        srim (iv) = zcorr*zscrim
        sshed(iv) = zcorr*zscshe 
      ENDIF

      ! Calculation of heterogeneous nucleation of cloud ice.
      ! This is done in this section, because we require water saturation
      ! for this process (i.e. the existence of cloud water) to exist.
      ! Hetrogeneous nucleation is assumed to occur only when no
      ! cloud ice is present and the temperature is below a nucleation
      ! threshold.
      IF( (tg <= 267.15_wp) .AND. (qig <= 0.0_wp)) THEN
        IF  (lsuper_coolw .OR. lorig_icon) THEN
          znin  = MIN( fxna_cooper(tg), znimax )
          snuc(iv) = zmi0 * z1orhog(iv) * znin * zdtr
        ELSE 
          znin  = MIN( fxna(tg), znimax )
          snuc(iv) = zmi0 * z1orhog(iv) * znin * zdtr
        END IF
      ENDIF
      ! Calculation of in-cloud rainwater freezing
      IF (tg < ztrfrz .AND. qrg > 0.1_wp*qcg)  THEN
        IF (lsuper_coolw) THEN
          srfrz(iv) = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk(iv)
        ELSE 
          srfrz(iv) = zcrfrz*SQRT( (ztrfrz-tg)**3 )* zeln27o16qrk(iv)
        ENDIF
      ENDIF
      
! Calculation of reduction of depositional growth at cloud top (Forbes 2012)
      IF( k>k_start .AND. k<ke .AND. lred_depgrowth ) THEN
        znin  = MIN( fxna_cooper(tg), znimax )
        fnuc = MIN(znin/znimix, 1.0_wp)
        
        !! distance from cloud top
        IF( qv(iv,k-1) + qc(iv,k-1) < zqvsw_up(iv) .AND. qi(iv,k-1) + qs(iv,k-1) < zqmin) THEN  ! cloud top layer
          dist_cldtop(iv) = 0.0_wp    ! reset distance to cloud top layer
        ELSE
          dist_cldtop(iv) = dist_cldtop(iv) + dz(iv,k)
        END IF
         
! with asymptotic behaviour dz -> 0 (xxx)
!        reduce_dep(iv) = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
!                             dist_cldtop(iv)/dist_cldtop_ref + &
!                             (1.0_wp-reduce_dep_ref)*(zdh/dist_cldtop_ref)**4), 1.0_wp)
        
! without asymptotic behaviour dz -> 0
        reduce_dep(iv) = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
                          dist_cldtop(iv)/dist_cldtop_ref), 1.0_wp)
       
      END IF ! Reduction of dep. growth of snow/ice 
      
    ENDDO loop_over_qc

      !------------------------------------------------------------------------
      ! Section 6: Search for cold grid points with cloud ice and/or snow and
      !            calculation of the conversion rates involving qi and qs
      !------------------------------------------------------------------------

! also ic3
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs: DO i1d =1, ic3
      iv = ivdx3(i1d)

      qvg  =  qv(iv,k)
      qig  =  qi(iv,k)
      qsg  =  qs(iv,k)
      tg   =   t(iv,k)
      ppg  =   p(iv,k)
      rhog = rho(iv,k)
      llqi =  qig > zqmin

      IF (tg<=t0) THEN
        IF( lsuper_coolw .OR. lorig_icon) THEN
          znin   = MIN( fxna_cooper(tg), znimax )
        ELSE
          znin   = MIN( fxna(tg), znimax )
        END IF
        IF (lstickeff .OR. lorig_icon) THEN
          stickeff = MIN(EXP(0.09_wp*(tg-t0)),1.0_wp)
          stickeff = MAX(stickeff, zceff_min, zceff_fac*(tg-tmin_iceautoconv))
        ELSE !original sticking efficiency of cloud ice
          stickeff = MIN(EXP(0.09_wp*(tg-t0)),1.0_wp)
          stickeff = MAX(stickeff,0.2_wp)
        END IF
        zmi      = MIN( rhog*qig/znin, zmimax )
        zmi      = MAX( zmi0, zmi )
        zsvmax   = (qvg - zqvsi(iv)) * zdtr
        zsagg    = zcagg(iv) * EXP(ccsaxp*LOG(zcslam(iv))) * qig
        zsagg    = MAX( zsagg, 0.0_wp ) * stickeff
        znid     = rhog * qig/zmi
        IF (llqi) THEN
          zlnlogmi = LOG (zmi)
          zsidep   = zcidep(iv) * znid * EXP(0.33_wp * zlnlogmi) * (qvg - zqvsi(iv))
        ELSE
          zsidep = 0.0_wp
        ENDIF
        zsvidep   = 0.0_wp
        zsvisub   = 0.0_wp
        ! for sedimenting quantities the maximum 
        ! allowed depletion is determined by the predictor value. 
        IF (lsedi_ice .OR. lorig_icon) THEN
          zsimax  = zzai(iv)*z1orhog(iv)*zdtr
        ELSE
          zsimax  = qig*zdtr
        ENDIF
        !
        IF( zsidep > 0.0_wp ) THEN
          IF (lred_depgrowth ) THEN
            zsidep = zsidep * reduce_dep(iv)  !FR new: SLW reduction
          END IF
          zsvidep = MIN( zsidep, zsvmax )
        ELSEIF (zsidep < 0.0_wp ) THEN
          zsvisub = - MAX(-zsimax, zsvmax )
        ENDIF
        zsiau = zciau * MAX( qig - qi0, 0.0_wp ) * stickeff
        IF (llqi) THEN
          zlnlogmi = LOG(zmsmin/zmi)
          zztau    = 1.5_wp*( EXP(0.66_wp*zlnlogmi) - 1.0_wp)
          zsdau    = zsvidep/MAX(zztau,zeps)
        ELSE
          zsdau    =  0.0_wp
        ENDIF
        zsicri  = zcicri * qig * zeln7o8qrk(iv)
        zsrcri  = zcrcri * (qig/zmi) * zeln13o8qrk(iv)
        zxfac   = 1.0_wp + zbsdep(iv) * EXP(ccsdxp*LOG(zcslam(iv)))
        zssdep  = zcsdep(iv) * zxfac * ( qvg - zqvsi(iv) ) / (zcslam(iv)+zeps)**2

        ! Check for maximal depletion of vapor by sdep
        IF (zssdep > 0.0_wp) THEN
          IF (lred_depgrowth ) THEN
            zssdep = zssdep*reduce_dep(iv)!FR new: SLW reduction
          END IF
          zssdep = MIN(zssdep, zsvmax-zsvidep)
        END IF
        
        zsisum = zsiau + zsdau + zsagg + zsicri + zsvisub
        zcorr  = 0.0_wp
        IF( zsimax > 0.0_wp ) zcorr  = zsimax / MAX( zsimax, zsisum )
        sidep(iv)  = zsvidep - zcorr*zsvisub
        sdau (iv)  = zcorr*zsdau
        siau (iv)  = zcorr*zsiau
        sagg (iv)  = zcorr*zsagg
        sicri(iv)  = zcorr*zsicri

        ! Allow growth of snow only if the existing amount of snow is sufficiently large 
        ! for a meaningful distiction between snow and cloud ice
        IF (qsg > 1.e-7_wp .OR. zssdep <= 0.0_wp) ssdep(iv) = zssdep
        IF (qsg > 1.e-7_wp)  srcri(iv) = zsrcri

      !------------------------------------------------------------------------
      ! Section 7: Search for warm grid points with cloud ice and/or snow and
      !            calculation of the melting rates of qi and ps
      !------------------------------------------------------------------------

      ELSE ! tg > 0
        IF (lsedi_ice .OR. lorig_icon) THEN
          simelt(iv) = zzai(iv)*z1orhog(iv)*zdtr
        ELSE
          simelt(iv) = qig*zdtr
        ENDIF
        zqvsw0      = zpvsw0 / (rhog * r_v * tg)
        zx1         = (tg - t0) + zasmel*(qvg - zqvsw0)
        zx2         = 1.0_wp + zbsmel * zeln5o24qsk(iv)
        zssmelt     = zcsmel * zx1 * zx2 * zeln2o3qsk(iv)
        ssmelt(iv) = MAX( zssmelt, 0.0_wp )
      ENDIF ! tg

    ENDDO loop_over_qs

    !--------------------------------------------------------------------------
    ! Section 8: Search for grid points with rain in subsaturated areas
    !            and calculation of the evaporation rate of rain
    !--------------------------------------------------------------------------

! ic6
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qr_nocloud: DO i1d =1, ic6
      iv = ivdx6(i1d)

      qvg = qv(iv,k)
      tg  =  t(iv,k)
      rhog = rho(iv,k)

      zlnqrk      = LOG (zqrk(iv))
      zx1         = 1.0_wp + zbev * EXP (zbevxp  * zlnqrk)
      ! Limit evaporation rate in order to avoid overshoots towards supersaturation
      ! the pre-factor approximates (esat(T_wb)-e)/(esat(T)-e) at temperatures between 0 degC and 30 degC
      temp_c = tg - t0
      maxevap     = (0.61_wp-0.0163_wp*temp_c+1.111e-4_wp*temp_c**2)*(zqvsw(iv)-qvg)/zdt
      sev(iv)    = MIN(zcev*zx1*(zqvsw(iv) - qvg) * EXP (zcevxp  * zlnqrk), maxevap)

      ! Calculation of below-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        IF ( lsuper_coolw ) THEN
          !FR new: reduced rain freezing rate
          srfrz(iv) = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk(iv)
        ELSE
          srfrz(iv)  = zcrfrz*SQRT( (ztrfrz-tg)**3 ) * zeln27o16qrk(iv)
        ENDIF
      ENDIF

    ENDDO loop_over_qr_nocloud

    !--------------------------------------------------------------------------
    ! Section 9: Calculate the total tendencies of the prognostic variables.
    !            Update the prognostic variables in the interior domain.
    !--------------------------------------------------------------------------

    loop_over_all_iv: DO iv = iv_start, iv_end

        qrg  = qr(iv,k)
        qsg  = qs(iv,k)
        qig  = qi(iv,k)
        rhog = rho(iv,k)

#ifdef __COSMO__
        qvg  = qv(iv,k)
        qcg  = qc(iv,k)
         tg  = t (iv,k)
        ppg  = p (iv,k)
#endif

        zsrmax = zzar(iv)*z1orhog(iv)*zdtr
        zssmax = zzas(iv)*z1orhog(iv)*zdtr
        zsrsum = sev(iv) + srfrz(iv) + srcri(iv)
        zcorr  = 1.0_wp
        IF(zsrsum > 0._wp) THEN
          zcorr  = zsrmax / MAX( zsrmax, zsrsum )
        ENDIF
        sev  (iv) = zcorr*sev(iv)
        srfrz(iv) = zcorr*srfrz(iv)
        srcri(iv) = zcorr*srcri(iv)
        ssmelt(iv) = MIN(ssmelt(iv), zssmax)
        IF (ssdep(iv) < 0.0_wp ) THEN
          ssdep(iv) = MAX(ssdep(iv), - zssmax)
        ENDIF
        zqvt = sev(iv)   - sidep(iv) - ssdep(iv)  - snuc(iv)
        zqct = simelt(iv)- scau(iv)  - scfrz(iv)  - scac(iv)   - sshed(iv) - srim(iv) 
        zqit = snuc(iv)  + scfrz(iv) - simelt(iv) - sicri(iv)  + sidep(iv) - sdau(iv)            &
                                                               - sagg(iv) - siau(iv)
        zqrt = scau(iv)  + sshed(iv) + scac(iv)   + ssmelt(iv) - sev(iv) - srcri(iv) - srfrz(iv) 
        zqst = siau(iv)  + sdau(iv)  + sagg(iv)   - ssmelt(iv) + sicri(iv) + srcri(iv)           &
                                                               + srim(iv) + ssdep(iv) + srfrz(iv)
        ztt = z_heat_cap_r*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst) )

#ifdef __COSMO__
        IF (lsppt) THEN
!US how to check?          IF(ntstep>0.AND.i>=istart.and.i<=iend.and.j>=jstart.and.j<=jend) THEN
            CALL apply_tqx_tend_adj(itype_qxpert_rn, itype_qxlim_rn, ppg,    &
                                    tg,  qvg,  qcg,  qig,  qrg,  qsg,        &
                    pstoph(iv,k),  ztt, zqvt, zqct, zqit, zqrt, zqst, ldum)
!US          ENDIF
        ENDIF
#endif

        ! Update variables and add qi to qrs for water loading 
        IF (lsedi_ice .OR. lorig_icon) THEN
          qig = MAX ( 0.0_wp, (zzai(iv)*z1orhog(iv) + zqit*zdt)*zimi(iv))
        ELSE
          qig = MAX ( 0.0_wp, qig + zqit*zdt)
        END IF
        qrg = MAX ( 0.0_wp, (zzar(iv)*z1orhog(iv) + zqrt*zdt)*zimr(iv))
        qsg = MAX ( 0.0_wp, (zzas(iv)*z1orhog(iv) + zqst*zdt)*zims(iv))

        !----------------------------------------------------------------------
        ! Section 10: Complete time step
        !----------------------------------------------------------------------

        IF ( k /= ke) THEN
          ! Store precipitation fluxes and sedimentation velocities 
          ! for the next level
          zprvr(iv) = qrg*rhog*zvzr(iv)
          zprvs(iv) = qsg*rhog*zvzs(iv)
          zprvi(iv) = qig*rhog*zvzi(iv)

          IF (zprvr(iv) <= zqmin) zprvr(iv)=0.0_wp
          IF (zprvs(iv) <= zqmin) zprvs(iv)=0.0_wp
          IF (zprvi(iv) <= zqmin) zprvi(iv)=0.0_wp

          ! for the latent heat nudging
          IF (ldass_lhn) THEN
            IF (lsedi_ice .OR. lorig_icon) THEN
              qrsflux(iv,k) = zprvr(iv)+zprvs(iv)+zprvi(iv)
              qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv)+zpki(iv)) 
            ELSE 
              qrsflux(iv,k) = zprvr(iv)+zprvs(iv)
              qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv))
            ENDIF
          ENDIF

          IF (qrg+qr(iv,k+1) <= zqmin) THEN
            zvzr(iv)= 0.0_wp
          ELSE
            zvzr(iv)= zvz0r * EXP(zvzxp*LOG((qrg+qr(iv,k+1))*0.5_wp*rhog)) * zrho1o2(iv)
          ENDIF
          IF (qsg+qs(iv,k+1) <= zqmin) THEN
            zvzs(iv)= 0.0_wp
          ELSE
            zvzs(iv)= zvz0s(iv) * EXP(ccswxp*LOG((qsg+qs(iv,k+1))*0.5_wp*rhog)) * zrho1o2(iv)
          ENDIF
          IF (qig+qi(iv,k+1) <= zqmin ) THEN
            zvzi(iv)= 0.0_wp
          ELSE
            zvzi(iv)= zvz0i * EXP(zbvi*LOG((qig+qi(iv,k+1))*0.5_wp*rhog)) * zrhofac_qi(iv)
          ENDIF

        ELSE
          ! Precipitation fluxes at the ground
          prr_gsp(iv) = 0.5_wp * (qrg*rhog*zvzr(iv) + zpkr(iv))
          IF (lsedi_ice .OR. lorig_icon) THEN
            prs_gsp(iv) = 0.5_wp * (rhog*(qsg*zvzs(iv)+qig*zvzi(iv)) + zpks(iv)+zpki(iv))
          ELSE
            prs_gsp(iv) = 0.5_wp * (qsg*rhog*zvzs(iv) + zpks(iv))
          END IF
          

          ! for the latent heat nudging
          IF (ldass_lhn) &
            qrsflux(iv,k) = prr_gsp(iv)+prs_gsp(iv)

        ENDIF

        ! Update of prognostic variables or tendencies
        qr (iv,k) = qrg
        qs (iv,k) = qsg
        qi (iv,k) = qig
!        qrs(iv,k   ) = qrg+qsg+qig       !qrs is now computed outside
        t  (iv,k) = t (iv,k) + ztt*zdt 
        qv (iv,k) = MAX ( 0.0_wp, qv(iv,k) + zqvt*zdt )
        qc (iv,k) = MAX ( 0.0_wp, qc(iv,k) + zqct*zdt )

      ENDDO loop_over_all_iv

  IF (izdebug > 15) THEN
    ! Check for negative values
     DO iv = iv_start, iv_end
        IF (qr(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(a)') ' WARNING: cloudice, negative value in qr'
          CALL message('',message_text)
        ENDIF
        IF (qc(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(a)') ' WARNING: cloudice, negative value in qc'
          CALL message('',message_text)
        ENDIF
        IF (qi(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(a)') ' WARNING: cloudice, negative value in qi'
          CALL message('',message_text)
        ENDIF
        IF (qs(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(a)') ' WARNING: cloudice, negative value in qs'
          CALL message('',message_text)
        ENDIF
        IF (qv(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(a)') ' WARNING: cloudice, negative value in qv'
          CALL message('',message_text)
        ENDIF
    ENDDO
  ENDIF

#if defined (__COSMO__)
  ! Do a final saturation adjustment for new values of t, qv and qc
    CALL satad ( 1, t(:,k), qv(:,k),                  &
               qc(:,k), t(:,k), p(:,k),               &
               zdummy(:,1),zdummy(:,2),zdummy(:,3),   &
               zdummy(:,4),zdummy(:,5),zdummy(:,6),   &
               zdummy(:,7),zdummy(:,8),               &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv, &
               rvd_m_o, lh_v, z_heat_cap_r, cp_d,     &
               nvec, 1, iv_start, iv_end, 1 , 1 )

  IF ( ldiabf_lh ) THEN
    ! compute temperature increment due to latent heat
!CDIR COLLAPSE
    tinc_lh(:,k) = tinc_lh(:,k) + t(:,k)
  ENDIF
#endif

ENDDO loop_over_levels

!------------------------------------------------------------------------------
! final tendency calculation for ICON
!
! Note: as soon as we have a new satad subroutine in ICON, this tendency
! calculation will be done in the k-loop and the original 3D variables wont
! be used to store the new values. Then we wont need the _in variables anymore.
!------------------------------------------------------------------------------

! calculated pseudo-tendencies

  IF ( lldiag_ttend ) THEN
    DO k=k_start,ke
      DO iv=iv_start,iv_end
        ddt_tend_t (iv,k) = (t (iv,k) - t_in (iv,k))*zdtr
      END DO
    END DO
  ENDIF

  IF ( lldiag_qtend ) THEN
    DO k=k_start,ke
      DO iv=iv_start,iv_end
        ddt_tend_qv(iv,k) = MAX(-qv_in(iv,k)*zdtr,(qv(iv,k) - qv_in(iv,k))*zdtr)
        ddt_tend_qc(iv,k) = MAX(-qc_in(iv,k)*zdtr,(qc(iv,k) - qc_in(iv,k))*zdtr)
        ddt_tend_qi(iv,k) = MAX(-qi_in(iv,k)*zdtr,(qi(iv,k) - qi_in(iv,k))*zdtr)
        ddt_tend_qr(iv,k) = MAX(-qr_in(iv,k)*zdtr,(qr(iv,k) - qr_in(iv,k))*zdtr)
        ddt_tend_qs(iv,k) = MAX(-qs_in(iv,k)*zdtr,(qs(iv,k) - qs_in(iv,k))*zdtr)
      END DO
    END DO
  ENDIF

  IF (izdebug > 15) THEN
    CALL message('gscp_cloudice', 'UPDATED VARIABLES')
   WRITE(message_text,'(a,2E20.9)') 'cloudice T= ',&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'cloudice qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'cloudice qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'cloudice qi= ',&
    MAXVAL( qi(:,:)), MINVAL(qi(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'cloudice qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'cloudice qs= ',&
    MAXVAL( qs(:,:)), MINVAL(qs(:,:) )
    CALL message('', TRIM(message_text))
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine cloudice
!------------------------------------------------------------------------------

END SUBROUTINE cloudice

!==============================================================================

END MODULE gscp_cloudice
