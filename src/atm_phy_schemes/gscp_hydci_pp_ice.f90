!>
!! cloud microphysics
!!
!! !------------------------------------------------------------------------------
!!
!! @par Description of *init_hydci_pp_ice*:
!!   Initializes some parameters for hydci_pp_ice
!! @par Description of *hydci_pp*:
!!   This module procedure calculates the rates of change of temperature, cloud
!!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
!!   processes related to the formation of grid scale precipitation. This
!!   includes the sedimentation of rain and snow. The precipitation fluxes at
!!   the surface are also calculated here.
!! @par Description of *satad*:
!!   Current COSMO-saturation adjustment, will later be accomplisehd by
!!   a mass conservative scheme for ICONAM
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
!! @par reference   This is an adaption of subroutine hydci_pp in file src_gscp.f90
!!  of the lm_f90 library (COSMO code). Equation numbers refer to
!!  Doms, Foerstner, Heise, Herzog, Raschendorfer, Schrodin, Reinhardt, Vogel
!!    (September 2005): "A Description of the Nonhydrostatic Regional Model LM",
!!
!------------------------------------------------------------------------------
!!
!!
!! @par Revision History
!! implemented into ICON by Felix Rieper (2012-06)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE gscp_hydci_pp_ice

!------------------------------------------------------------------------------
!>
!! Description:
!!
!!   The subroutines in the module "gscp" calculate the rates of change of
!!   temperature, cloud condensate and water vapor due to cloud microphysical
!!   processes related to the formation of grid scale clouds and precipitation.
!!   In the COSMO model the microphysical subroutines are either
!!   called from "organize_gscp" or from "organize_physics" itself.
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

USE data_parameters, ONLY :   &
  ireals,    & !! KIND-type parameter for real variables
  iintegers    !! KIND-type parameter for standard integer vaiables

USE data_constants  , ONLY :   &
    pi,           & !!
!! 2. physical constants and related variables
!! -------------------------------------------
    t0=>t0_melt,  & !! melting temperature of ice
    r_d,          & !! gas constant for dry air
    r_v,          & !! gas constant for water vapour
    rdv,          & !! r_d / r_v
    o_m_rdv,      & !! 1 - r_d/r_v
    rvd_m_o,      & !! r_v/r_d - 1
    cp_d,         & !! specific heat of dry air
    cpdr,         & !! 1 / cp_d
    con_m,        & !! kinematic viscosity (m2/s) !CK<    
    lh_v,         & !! latent heat of vapourization
    lh_f,         & !! latent heat of fusion
    lh_s,         & !! latent heat of sublimation
    g,            & !! acceleration due to gravity
    rho_w,        & !! density of liquid water (kg/m^3)
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
    qc0             !! cloud water threshold for autoconversion

!! end of data_constants

!!USE data_fields     , ONLY :   &
!!     tinc_lh    !! temperature increment due to latent heat      (  K  )

!------------------------------------------------------------------------------

USE data_gscp, ONLY: &          ! all variables are used here
!
!   variables for hydci_pp_ice
!  
    ccsrim,    & !
    ccsagg,    & !
    ccsdep,    & !
    ccsvel,    & !
    ccsvxp,    & !
    ccslam,    & !
    ccslxp,    & !
    ccsaxp,    & !
    ccsdxp,    & !
    ccshi1,    & !
    ccdvtp,    & !
    ccidep,    & !
    ccswxp,    & !
    zconst,    & !
    zcev,      & !
    zbev,      & !
    zcevxp,    & !
    zbevxp,    & !
    zvzxp,     & !
    zvz0r,     & !
!
!   variables for hydci_pp_ice and hydci_pp_gr
!
    v0snow,    & ! factor in the terminal velocity for snow
    mu_rain,   & ! 
    rain_n0_factor, &
    cloud_num, &
  !CK>
    vtxexp,    & !!
    kc_c1,     & !!
    kc_c2        !!


!------------------------------------------------------------------------------


USE data_runcontrol , ONLY :   &
     ldiabf_lh
!------------------------------------------------------------------------------
!! COSMO environment modules
!------------------------------------------------------------------------------

USE data_parallel,            ONLY : my_cart_id  !! rank of this subdomain in
                                                 !! the cartesian communicator
USE pp_utilities,             ONLY : gamma_fct   !! Gamma function
USE meteo_utilities,          ONLY : satad       !! saturation adjustment
USE environment,              ONLY : collapse, model_abort, get_free_unit, release_unit

#endif

#ifdef __GME__

USE prognostic_pp, ONLY : pi,gamma_fct, &
    t0=>T_melt,   & !! melting temperature of ice
    r_d,          & !! gas constant for dry air
    r_v,          & !! gas constant for water vapour
    rdv,          & !! r_d / r_v
    o_m_rdv,      & !! 1 - r_d/r_v
    rvd_m_o,      & !! r_v/r_d - 1
    cp_d=>cp,     & !! specific heat of dry air
    cpdr=>cpr,    & !! 1 / cp_d
    lh_v=>Hl_c,   & !! latent heat of vapourization
    lh_f=>Hl_f,   & !! latent heat of fusion
    lh_s=>Hl_s,   & !! latent heat of sublimation
    g,            & !! acceleration due to gravity
    b1,           & !! variables for computing the saturation vapour pressure
    b2w=>b2_w,    & !! over water (w) and ice (i)
    b2i=>b2_i,    & !!               -- " --
    b3,           & !!               -- " --
    b4w=>b4_w,    & !!               -- " --
    b234w=>b234_w,& !!               -- " --
    b4i=>b4_i,    & !!               -- " --
    qi0=>zqi0,    & !! cloud ice threshold for autoconversion
    qc0=>zqc0,    & !! cloud water threshold for autoconversion
    rho_w           !! density of liquid water (kg/m^3)

#endif

#ifdef __ICON__

USE mo_kind,               ONLY: ireals=>wp     , &
                                 iintegers=>i4
USE mo_math_constants    , ONLY: pi
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                                 r_d   => rd    , & !! gas constant for dry air
                                 rvd_m_o=>vtmpc1, & !! r_v/r_d - 1
                                 o_m_rdv        , & !! 1 - r_d/r_v
                                 rdv            , & !! r_d / r_v
                                 con_m          , & !! kinematic viscosity (m2/s) CK<
                                 lh_v  => alv   , & !! latent heat of vapourization
                                 lh_s  => als   , & !! latent heat of sublimation
                                 lh_f  => alf   , & !! latent heat of fusion
                                 cp_d  => cpd   , & !! specific heat of dry air at constant press
                                 cpdr  => rcpd  , & !! (spec. heat of dry air at constant press)^-1
                                 cvdr  => rcvd  , & !! (spec. heat of dry air at const vol)^-1
                                 b3    => tmelt , & !! melting temperature of ice/snow
                                 rho_w => rhoh2o, & !! density of liquid water (kg/m^3)
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
!                                 spec_humi          !! Specific humidity 
USE mo_exception,          ONLY: message, message_text
!xxx: common module COSMO/ICON, all variables are used here
USE gscp_data

USE data_hydci_pp_ice,     ONLY:    afrac_dust, &  !! look-up table of activated fraction of dust particles acting as ice nuclei
                                    afrac_soot, &  !! ... of soot particles
                                    afrac_orga     !! ... of organic material

USE mo_run_config,         ONLY: ldass_lhn

#endif

!==============================================================================

IMPLICIT NONE
PRIVATE


!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: hydci_pp_ice


!------------------------------------------------------------------------------
!> Parameters and variables which are global in this module
!------------------------------------------------------------------------------

#ifdef __COSMO__
CHARACTER(132) :: message_text = ''
#endif




!> Namelist Variables for hydci_pp_ice and hydci_pp_gr
! -------------------------------------------------
!    all moved to data_gscp module


CONTAINS

!==============================================================================

ELEMENTAL FUNCTION ice_nuclei_number(temp)
IMPLICIT NONE

REAL (KIND=ireals)              :: ice_nuclei_number
REAL (KIND=ireals), INTENT(IN)  :: temp

ice_nuclei_number = 1.0E2_ireals * EXP(0.2_ireals * (t0 - temp))

END FUNCTION ice_nuclei_number

!CK>
!==============================================================================
! dimensionless growth time scale (ECHAM Code)
! B.Kaercher and S. Solomon, JGR 104(D22), 27441-27459, 1999 
!==============================================================================

REAL(KIND=ireals) ELEMENTAL FUNCTION dep_growth_timescale(b,y,x_0) 
IMPLICIT NONE
REAL(KIND=ireals), INTENT(in) :: b,y,x_0

REAL(KIND=ireals), PARAMETER  :: &
  SQ31 = 0.577350269_ireals,       &
  SIX1 = 0.166666666_ireals

REAL(KIND=ireals) ::  X,F1,F10,F2,F20
    
IF (y.LE.x_0) THEN
  dep_growth_timescale = 0.
ELSE
  X     = MIN( y, 0.9999999_ireals )
  F1    = SIX1 * LOG( (1.+X +X**2)  / (1.-X)**2 ) 
  F10   = SIX1 * LOG( (1.+x_0+x_0**2) / (1.-x_0)**2 ) 
  F2    = SQ31 * ATAN( SQ31*(1.+2.*X) )
  F20   = SQ31 * ATAN( SQ31*(1.+2.*x_0) )      
  dep_growth_timescale = (b+1.)*(F1-F10) + (b-1.)*(F2-F20)            
END IF

END FUNCTION dep_growth_timescale
!CK<

#ifdef __COSMO__
!==============================================================================

ELEMENTAL FUNCTION sat_pres_water(temp)
IMPLICIT NONE

REAL (KIND=ireals)              :: sat_pres_water
REAL (KIND=ireals), INTENT(IN)  :: temp

sat_pres_water = b1*EXP( b2w*(temp-b3)/(temp-b4w) )

END FUNCTION sat_pres_water

!==============================================================================

ELEMENTAL FUNCTION sat_pres_ice(temp)
IMPLICIT NONE

REAL (KIND=ireals)              :: sat_pres_ice
REAL (KIND=ireals), INTENT(IN)  :: temp

sat_pres_ice = b1*EXP( b2i*(temp-b3)/(temp-b4i) )

END FUNCTION sat_pres_ice

!==============================================================================

ELEMENTAL FUNCTION spec_humi(pvap,pres)
IMPLICIT NONE

REAL (KIND=ireals)              :: spec_humi
REAL (KIND=ireals), INTENT(IN)  :: pres,pvap

spec_humi = rdv*pvap/( pres - o_m_rdv*pvap )

END FUNCTION spec_humi

!==============================================================================
SUBROUTINE message (name, text, all_print)
IMPLICIT NONE

CHARACTER (*) :: name, text
LOGICAL, INTENT(in), OPTIONAL :: all_print

LOGICAL :: lprint

IF (PRESENT(all_print)) THEN
  lprint = all_print
ELSE
  lprint = .FALSE.
ENDIF

IF (lprint .OR. my_cart_id==0) &
  WRITE (*,*) TRIM(name), TRIM(text)

END SUBROUTINE message

#endif


!==============================================================================
!> Module procedure "hydci_pp_ice" in "gscp" for computing effects of grid scale
!!  precipitation including cloud water, cloud ice, rain and snow in
!!  context with the Leapfrog and the Runge-Kutta time-integration
!------------------------------------------------------------------------------

SUBROUTINE hydci_pp_ice (             &
  nvec,ke,ke1,                        & !> array dimensions !CK> ke1
  ivstart,ivend,kstart,               & !! optional start/end indicies
  idbg,                               & !! optional debug level
  zdt,dz,                             & !! numerics parameters
  w,tke,                              & !! !CK> w,tke
  t,p,rho,qv,qc,qi,qni,qni_nuc,       & !! prognostic variables !CK> qni,qni_nuc
  qr,qs,                              & !! prognostic variables
#ifdef __ICON__
  !xxx: this should become a module variable, e.g. in a new module mo_data_gscp.f90
  qi0,qc0,                            & !! cloud ice/water threshold for autoconversion
#endif
  prr_gsp,prs_gsp,                    & !! surface precipitation rates
#ifdef __COSMO__
  tinc_lh,                            & !  t-increments due to latent heating 
#endif
  qrsflux,                            & !  total precipitation flux
  l_cv,                               &
  ldiag_ttend,     ldiag_qtend     , &
  ddt_tend_t     , ddt_tend_qv      , &
  ddt_tend_qc    , ddt_tend_qi      , & !> ddt_tend_xx are tendencies
  ddt_tend_qr    , ddt_tend_qs      , & !!    necessary for dynamics
  ddt_diag_au    , ddt_diag_ac      , & !!
  ddt_diag_ev    , ddt_diag_nuc     , & !! ddt_diag_xxx are optional
  ddt_diag_idep  , ddt_diag_sdep    , & !!   diagnostic tendencies of all
  ddt_diag_agg   , ddt_diag_rim     , & !!   individual microphysical processes
  ddt_diag_rcri  , ddt_diag_icri    , & !!
  ddt_diag_dau   , ddt_diag_iau     , & !!
  ddt_diag_imelt , ddt_diag_smelt   , & !!
  ddt_diag_cfrz  , ddt_diag_rfrz    , & !!
  ddt_diag_shed                       ) !!

!------------------------------------------------------------------------------
!>
!! Description:
!!   This module procedure calculates the rates of change of temperature, cloud
!!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
!!   processe related to the formation of grid scale precipitation. The
!!   variables are updated in this subroutine. Rain and snow are prognostic
!!   variables. The precipitation fluxes at the surface are stored on the
!!   corresponding global fields.
!!   The subroutine relies on conversion terms used in hydci.
!!
!! Method:
!!   The sedimentation of rain and snow is computed implicitly.
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
    nvec             ,    & !> number of horizontal points
    ke               ,    & !! number of grid points in vertical direction
!CK>
    ke1                     ! index of the lowest model half level (=ke+1)
!CK<
    
  INTEGER, INTENT(IN), OPTIONAL ::  &
    ivstart    ,    & !> optional start index for horizontal direction
    ivend      ,    & !! optional end index   for horizontal direction
    kstart    ,    & !! optional start index for the vertical index
    idbg             !! optional debug level

  REAL(KIND=ireals), INTENT(IN) :: &
    zdt                    !> time step for integration of microphysics     (  s  )

#ifdef __ICON__
  REAL(KIND=ireals), INTENT(IN) :: &
    qi0,qc0          !> cloud ice/water threshold for autoconversion
#endif

  REAL(KIND=ireals), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )
  
  REAL(KIND=ireals), DIMENSION(:,:), INTENT(IN) ::   &   ! dim (ie,ke1)
  w                 , & !! vertical wind speed (defined on half levels)  ( m/s )
  tke                   !! SQRT(2*TKE); TKE='turbul. kin. energy'        ( m/s )
                        !! (defined on half levels)  
  
  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv                   !! if true, cv is used instead of cp

  LOGICAL, INTENT(IN), OPTIONAL :: &
    ldiag_ttend,         & ! if true, temperature tendency shall be diagnosed
    ldiag_qtend            ! if true, moisture tendencies shall be diagnosed

  REAL(KIND=ireals), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    !CK>
    qni             ,    & !! specific cloud ice number                     ( 1/kg)   
    qni_nuc         ,    & !! activated ice nuclei                          ( 1/kg)
    !CK<    
    qr              ,    & !! specific rain content                         (kg/kg)
    qs                     !! specific snow content                         (kg/kg)

#ifdef __COSMO__
  REAL(KIND=ireals), INTENT(INOUT) :: &
       tinc_lh(:,:)   ,  & ! temperature increments due to heating             ( K/s )   
#endif

  REAL(KIND=ireals), INTENT(INOUT) :: &
       qrsflux(:,:)       ! total precipitation flux (nudg)

  REAL(KIND=ireals), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    prs_gsp                !! precipitation rate of snow, grid-scale        (kg/(m2*s))

  REAL(KIND=ireals), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &     ! dim (ie,ke)
    ddt_tend_t       , & !> tendency T                                       ( 1/s )
    ddt_tend_qv      , & !! tendency qv                                      ( 1/s )
    ddt_tend_qc      , & !! tendency qc                                      ( 1/s )
    ddt_tend_qi      , & !! tendency qi                                      ( 1/s )
    ddt_tend_qr      , & !! tendency qr                                      ( 1/s )
    ddt_tend_qs          !! tendency qs                                      ( 1/s )

  REAL(KIND=ireals), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &   ! dim (ie,ke)
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


  !! Local parameters: None, parameters are in module header, data_gscp or data_constants
  !! ----------------
  
  !> Local scalars:
  !! -------------
  
  INTEGER (KIND=iintegers) ::  &
    iv, k             !> loop indices

  REAL    (KIND=ireals   ) :: nnr

  REAL    (KIND=ireals   ) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)

  INTEGER ::  &
    iv_start     ,    & !> start index for horizontal direction
    iv_end       ,    & !! end index for horizontal direction
    k_start      ,    & !! model level where computations start
    izdebug             !! debug level

  REAL    (KIND=ireals   ) ::  &
    fpvsw, fpvsi, fqvs,& ! name of statement functions
    fxna ,             & ! statement function for ice crystal number
    ztx  , zpv  , zpx ,& ! dummy arguments for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    !CK>
    zximax,            & ! maximum mass of ice crystals
    zximin,            & ! minimum mass of ice crystals
    zxcloud,           & ! average mass of cloud droplets    
    !CK<    
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
!    zdt,               & ! timestep for integration (water / ice )
    zdtr ,             & ! reciprocal of timestep for integration
    zscau , zscac  , zscrim , zscshe,         & ! local values of the !CK removed zsnuc
    zsiau , zsagg  , zsidep , zsicri, zsrcri, & ! transfer rates
    zssdep , zssmelt,  & ! defined below        !CK removed zsdau
    zsrfrz,                                   & !
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum, zsrmax,    & ! terms for limiting  total rain water depletion
    zssmax,            & ! term for limiting snow depletion
    znid,              & ! number of cloud ice crystals for deposition !CK rm znin
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,   & ! terms for limiting total
    zqvsw,             & ! sat. specitic humidity at ice and water saturation
    zqvsidiff,         & !! qv-zqvsi !CK<
    zxfac, zx1, zx2, ztt, & ! some help variables !CK removed zztau
    ztau, zphi, zhi, zdvtp, ztc

  REAL    (KIND=ireals   ) ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice
    zqnit  ,& ! layer tendency of cloud ice !CK<    
    zqrt   ,& ! layer tendency of rain
    zqst      ! layer tendency of snow

  REAL (KIND=ireals)         ::       &
    mma(10), mmb(10)

  REAL    (KIND=ireals   ) ::         &
    zlnqrk,zlnqsk,                                       & !!
    qcg,tg,qvg,qrg,qsg,qig,rhog,ppg,alf,bet,m2s,m3s,hlp, & !!
!CK>
    qnig,qnig_nuc,asphere,ssi,wg,zri,xt,xs
!CK<
    
  LOGICAL :: &
    llqr,llqs,llqc,llqi  !   switch for existence of qr, qs, qc, qi

!CK>
  LOGICAL :: &
    lhom                       !! switch for homogeneous nucleation
  
  LOGICAL :: lldiag_ttend, lldiag_qtend

   REAL    (KIND=ireals   ), PARAMETER ::  &
    tau_mix    = 14000.0_ireals,  & !> mixing timescale for activated IN (2 hours = 7200 s)
    na_dust    = 162.e3_ireals,   & !! initial number density of dust [1/m3], Phillips08 value 162e3
    na_soot    =  15.e6_ireals,   & !! initial number density of soot [1/m3], Phillips08 value 15e6
    na_orga    = 177.e6_ireals,   & !! initial number density of organics [1/m3], Phillips08 value 177e6
    ni_het_max =  50.e3_ireals,   & !! max number of IN between 1-10 per liter, i.e. 1d3-10d3 [1/m3] 
    ni_hom_max = 5000e3_ireals      !! number of liquid aerosols between 100-5000 per liter   [1/m3]

  REAL    (KIND=ireals   ), PARAMETER ::  &
    rho0   = 1.225e+0_ireals, & ! reference air density
    ! Single bullets of Heymsfield and Iaquinta, see their Eq (A29) and text below
    zvz0i  = 3577_ireals,  & ! Factor in the terminal velocity for cloud ice
    zvzxi  = 1.31_ireals,  & ! Exponent in the terminal velocity for cloud ice
    zvimin = 0.01_ireals,  & ! Minimum terminal velocity for cloud ice
    zvimax = 1.00_ireals     ! Maximum terminal velocity for cloud ice

  ! some constants needed for Kaercher and Lohmann approach
  REAL(KIND=ireals), PARAMETER :: &
    rho_ice = 900.0             , &   !>..Materialdichte von Eis
    r_0     = 0.25e-6           , &   !! aerosol particle radius prior to freezing
    alpha_d = 0.5               , &   !! deposition coefficient (KL02; Spichtinger & Gierens 2009)
    k_b     = 1.38065e-23       , &   !! Boltzmann constant [J/K]
    M_w     = 18.01528e-3       , &   !! molecular mass of water [kg/mol]
    M_a     = 28.96e-3          , &   !! molecular mass of air [kg/mol]
    N_avo   = 6.022e23          , &   !! Avogadro number [1/mol]
    ma_w    = M_w / N_avo       , &   !! mass of water molecule [kg]
    svol    = ma_w / rho_ice          !! specific volume of a water molecule in ice
  
  !==============================================================================
  
  ! Look-up table for Phillips et al. 2008 nucleation
  INTEGER(KIND=iintegers), PARAMETER :: &
    ttmax  = 30,      &  !! sets limit for temperature in look-up table
    ssmax  = 60,      &  !! sets limit for ice supersaturation in look-up table
    ttstep = 2,       &  !! increment for temperature in look-up table
    ssstep = 1           !! increment for ice supersaturation in look-up table

!FR/CK: include file --> data_hydci_pp_ice.f90   
 
!CK<
  
  REAL(KIND=ireals), DIMENSION(nvec,ke) ::   &
    t_in               ,    & !> temperature                                   (  K  )
    qv_in              ,    & !! specific water vapor content                  (kg/kg)
    qc_in              ,    & !! specific cloud water content                  (kg/kg)
    qi_in              ,    & !! specific cloud ice   content                  (kg/kg)
    qr_in              ,    & !! specific rain content                         (kg/kg)
    qs_in                     !! specific snow content                         (kg/kg)


!! Local (automatic) arrays:
!! -------------------------
#ifdef __COSMO__
  REAL    (KIND=ireals   ) ::  &
    zpres       (nvec)        !! pressure
#endif

  REAL    (KIND=ireals   ) ::  &
    zqvsi       (nvec),     & !> sat. specitic humidity at ice and water saturation 
    zvzr        (nvec),     & !
    zvzs        (nvec),     & !
    zvzi        (nvec),     & ! !CK< ice sedimentation
    zpkr        (nvec),     & !
    zpks        (nvec),     & !
    zpki        (nvec),     & ! !CK< ice sedimentation
    zpkni       (nvec),     & ! !CK< ice sedimentation
    zprvr       (nvec),     & !
    zprvs       (nvec),     & !
    zprvi       (nvec),     & ! !CK< ice sedimentation
    zprvni      (nvec),     & ! !CK< ice sedimentation
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
    zimi        (nvec),     & ! !CK<
    zimni       (nvec),     & ! !CK<
    zzar        (nvec),     & !
    zzas        (nvec),     & !
    zzai        (nvec),     & ! !CK<
    zzani       (nvec),     & ! !CK< 
    zqrk        (nvec),     & !
    zqsk        (nvec),     & !
    !CK>
    zqik        (nvec),     & !
    zqnik       (nvec),     & !
    znidiag     (nvec),     & ! number of ice crystal from empirical diagnostic relation
    zDi         (nvec),     & ! mean diameter of cloud ice crystals
    zxi         (nvec),     & ! mean mass of cloud ice crystals
    vt_kc05     (nvec),     & ! terminal vel. ice
    zsi         (nvec),     & ! supersaturation over ice
    zpvsi       (nvec),     & ! saturation vapor pressure at over ice
    !CK<    
    zdtdh       (nvec),     & !
    z1orhog     (nvec),     & ! 1/rhog
    zrho1o2     (nvec),     & ! (rho0/rhog)**1/2
    zeln7o8qrk  (nvec),     & !
    zeln27o16qrk(nvec),     & !
    zeln13o8qrk (nvec),     & !
!    zeln3o16qrk (nvec),     & !
    zeln13o12qsk(nvec),     & !
    zeln5o24qsk (nvec),     & !
    zeln2o3qsk  (nvec)        !

  REAL    (KIND=ireals   ) ::  &
    scau   (nvec), & ! transfer rate due to autoconversion of cloud water
    scac   (nvec), & ! transfer rate due to accretion of cloud water
    !CK>
    snucn  (nvec), & !! transfer rate due to nucleation of cloud ice number
    snucq  (nvec), & !! transfer rate due to nucleation of cloud ice mass
    snucn2 (nvec), & !! transfer rate due to nucleation of cloud ice number
    snucq2 (nvec), & !! transfer rate due to nucleation of cloud ice mass
    !CK<    
    scfrz  (nvec), & ! transfer rate due homogeneous freezing of cloud water
    simelt (nvec), & ! transfer rate due melting of cloud ice
    sidep  (nvec), & ! transfer rate due depositional growth of cloud ice
    ssdep  (nvec), & ! transfer rate due depositional growth of snow
    srim   (nvec), & ! transfer rate due riming of snow
    sshed  (nvec), & ! transfer rate due shedding
    sicri  (nvec), & ! transfer rate due cloud ice collection by rain (sink qi)
    srcri  (nvec), & ! transfer rate due cloud ice collection by rain (sink qr)
    sagg   (nvec), & ! transfer rate due aggregation of snow and cloud ice
    siau   (nvec), & ! transfer rate due autoconversion of cloud ice
    ssmelt (nvec), & ! transfer rate due melting of snow
    sev    (nvec), & ! transfer rate due evaporation of rain
    srfrz  (nvec)    ! transfer rate due to rainwater freezing

  ! Dimensions and loop counter for storing the indices
  INTEGER (KIND=iintegers) ::  &
    ic1, ic2, ic3, ic4, ic5, ic6, ic7, i1d !CK> ic7

  !> Integer arrays for a better vectorization
  INTEGER (KIND=iintegers) ::  &
    ivdx1(nvec), & !!
    ivdx2(nvec), & !!
    ivdx3(nvec), & !!
    ivdx4(nvec), & !!
    ivdx5(nvec), & !!
    ivdx6(nvec), & !!
    ivdx7(nvec)    !! !CK< homogeneous nucleation

  !CK>
  INTEGER (KIND=iintegers) ::  &
    tt, ss, jj                 ! for look-up table 

  REAL(KIND=ireals)  :: infrac(3)
  REAL(KIND=ireals)  :: ni_hom,ri_hom,mi_hom
  REAL(KIND=ireals)  :: Xi_mck, Xi_mck_inv,tau_dep_ice_inv,tau_dep_snow_inv  
  REAL(KIND=ireals)  :: v_th,n_sat,flux,phi,cool,tau,delta,w_p,scr
  REAL(KIND=ireals)  :: ctau,acoeff(3),bcoeff(2),ri_dot
  REAL(KIND=ireals)  :: kappa,sqrtkap,ren,R_imfc,R_im,R_ik,ri_0
  REAL(KIND=ireals)  :: tgrow,ri_max,beta,xj,dxj,xmid,fmid,si_het,thom
  REAL(KIND=ireals)  :: Xbest , cbest, kc_b1, kc_b2, kc_a1, bracket
  !CK<

  
!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine hydci_pp_ice
!------------------------------------------------------------------------------



!> Statement functions
! -------------------

! saturation vapour pressure over water (fpvsw), over ice (fpvsi)
! and specific humidity at vapour saturation (fqvs)
  fpvsw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) )
  fpvsi(ztx)     = b1*EXP( b2i*(ztx-b3)/(ztx-b4i) )
  fqvs (zpv,zpx) = rdv*zpv/( zpx - o_m_rdv*zpv )
  
! Number of activate ice crystals;  ztx is temperature
  fxna(ztx)   = 1.0E2_ireals * EXP(0.2_ireals * (t0 - ztx))

!> Coeffs for moment relation based on 2nd moment (Field 2005)
  mma = (/5.065339_ireals, -0.062659_ireals, -3.032362_ireals, 0.029469_ireals, -0.000285_ireals, &
          0.312550_ireals,  0.000204_ireals,  0.003199_ireals, 0.000000_ireals, -0.015952_ireals /)
  mmb = (/0.476221_ireals, -0.015896_ireals,  0.165977_ireals, 0.007468_ireals, -0.000141_ireals, &
          0.060366_ireals,  0.000079_ireals,  0.000594_ireals, 0.000000_ireals, -0.003577_ireals /)

! Define reciprocal of heat capacity of dry air (at constant pressure vs at constant volume)

#ifdef __COSMO__
  z_heat_cap_r = cpdr
#endif

#ifdef __GME__
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
    znimax = fxna(zthn) ! Maximum number of cloud ice crystals
    zpvsw0 = fpvsw(t0)  ! sat. vap. pressure for t = t0

!CK>  
  asphere = 4._ireals/3._ireals * pi * rho_ice
  zxcloud = pi/6.0_ireals * rho_w * 50.e-6_ireals**3  ! Average size of cloud droplets
  zximin  = zami * 10.e-6_ireals**3      ! Minimum mean mass of cloud ice
  zximax  = zami * 250e-6_ireals**3      ! Maximum mean mass of cloud ice
!CK<
    
! Delete precipitation fluxes from previous timestep
!CDIR BEGIN COLLAPSE
    prr_gsp (:) = 0.0_ireals
    prs_gsp (:) = 0.0_ireals
    zpkr    (:) = 0.0_ireals
    zpks    (:) = 0.0_ireals
    zpki    (:) = 0.0_ireals !CK<
    zpkni   (:) = 0.0_ireals !CK<    
    zprvr   (:) = 0.0_ireals
    zprvs   (:) = 0.0_ireals
    zprvi   (:) = 0.0_ireals !CK<
    zprvni  (:) = 0.0_ireals !CK<    
    zvzr    (:) = 0.0_ireals
    zvzs    (:) = 0.0_ireals
    !CK>
    zzai    (:) = 0.0_ireals  
    zzani   (:) = 0.0_ireals
    zqik    (:) = 0.0_ireals  
    zqnik   (:) = 0.0_ireals  
    zimi    (:) = 0.0_ireals
    zimni   (:) = 0.0_ireals
    zvzr    (:) = 0.0_ireals
    zvzs    (:) = 0.0_ireals
    zvzi    (:) = 0.0_ireals
    vt_kc05 (:) = 0.0_ireals
    zDi     (:) = 0.0_ireals   
    zxi     (:) = zximin
    zsi     (:) = 0.0_ireals
    znidiag (:) = 0.0_ireals
    !CK<    
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
  zdtr  = 1.0_ireals / zdt

  ! add part of latent heating calculated in subroutine hydci_pp_ice to model latent
  ! heating field: subtract temperature from model latent heating field
  IF (ldass_lhn) THEN
      qrsflux(:,:) = 0.0_ireals
  ENDIF


! output for various debug levels
  IF (izdebug > 15) CALL message('','SRC_GSCP: Start of hydci_pp_ice')
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
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qni = ',MAXVAL(qni),MINVAL(qni) !CK<
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
    zcrim (:) = 0.0_ireals
    zcagg (:) = 0.0_ireals
    zbsdep(:) = 0.0_ireals
    zvz0s (:) = 0.0_ireals
    zn0s  (:) = zn0s0
!CDIR END

    ic1 = 0
    ic2 = 0
    ic4 = 0 !CK< cloud ice sedimentation
    
    DO iv = iv_start, iv_end

        qrg  = qr(iv,k)
        qsg  = qs(iv,k)
        qvg  = qv(iv,k)
        qcg  = qc(iv,k)
        qig  = qi(iv,k)
        qnig = qni(iv,k) !CK<       
        tg   = t(iv,k)
        ppg  = p(iv,k)
        rhog = rho(iv,k)

        !..for density correction of fall speeds
        z1orhog(iv) = 1.0_ireals/rhog
        zrho1o2(iv) = EXP(LOG(zrho0*z1orhog(iv))*x1o2)

        zqrk(iv)  = qrg  * rhog
        zqsk(iv)  = qsg  * rhog
        zqik(iv)  = qig  * rhog  !CK<
        zqnik(iv) = qnig * rhog  !CK<
      
        llqr = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin
        llqi = zqik(iv) > zqmin  !CK<
        
        zdtdh(iv) = 0.5_ireals * zdt / dz(iv,k)

        zzar(iv)   = zqrk(iv)/zdtdh(iv) + zprvr(iv) + zpkr(iv)
        zzas(iv)   = zqsk(iv)/zdtdh(iv) + zprvs(iv) + zpks(iv)
        !CK<
        zzai  (iv) = zqik(iv) / zdtdh(iv) + zprvi(iv)  + zpki(iv)
        zzani (iv) = zqnik(iv)/ zdtdh(iv) + zprvni(iv) + zpkni(iv)
        !CK<        

        IF (llqs) THEN
          ic1 = ic1 + 1
          ivdx1(ic1) = iv
        ENDIF
        IF (llqr) THEN
          ic2 = ic2 + 1
          ivdx2(ic2) = iv
        ENDIF
        !CK>
        IF (llqi) THEN
          ic4 = ic4 + 1
          ivdx4(ic4) = iv
        ENDIF
        !CK<  
     ENDDO

!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs_prepare: DO i1d = 1, ic1
      iv = ivdx1(i1d)

      qsg = qs(iv,k)
      tg  = t(iv,k)

      IF (isnow_n0temp == 1) THEN
        ! Calculate n0s using the temperature-dependent
        ! formula of Field et al. (2005)
        ztc = tg - t0
        ztc = MAX(MIN(ztc,0.0_ireals),-40.0_ireals)
        zn0s(iv) = zn0s1*EXP(zn0s2*ztc)
        zn0s(iv) = MIN(zn0s(iv),1e9_ireals)
        zn0s(iv) = MAX(zn0s(iv),1e6_ireals)
      ELSEIF (isnow_n0temp == 2) THEN
        ! Calculate n0s using the temperature-dependent moment
        ! relations of Field et al. (2005)
        ztc = tg - t0
        ztc = MAX(MIN(ztc,0.0_ireals),-40.0_ireals)

        nnr  = 3._ireals
        hlp = mma(1) + mma(2)*ztc + mma(3)*nnr + mma(4)*ztc*nnr &
          & + mma(5)*ztc**2 + mma(6)*nnr**2 + mma(7)*ztc**2*nnr &
          & + mma(8)*ztc*nnr**2 + mma(9)*ztc**3 + mma(10)*nnr**3
        alf = 10.0_ireals**hlp
        bet = mmb(1) + mmb(2)*ztc + mmb(3)*nnr + mmb(4)*ztc*nnr &
          & + mmb(5)*ztc**2 + mmb(6)*nnr**2 + mmb(7)*ztc**2*nnr &
          & + mmb(8)*ztc*nnr**2 + mmb(9)*ztc**3 + mmb(10)*nnr**3

        ! Here is the exponent bms=2.0 hardwired! not ideal! (Uli Blahak)
        m2s = qsg * rho(iv,k) / zams   ! UB rho added as bugfix
        m3s = alf*EXP(bet*LOG(m2s))

        hlp  = zn0s1*EXP(zn0s2*ztc)
        zn0s(iv) = 13.50_ireals * m2s**4 / m3s**3
        zn0s(iv) = MAX(zn0s(iv),0.5_ireals*hlp)
        zn0s(iv) = MIN(zn0s(iv),1e2_ireals*hlp)
        zn0s(iv) = MIN(zn0s(iv),1e9_ireals)
        zn0s(iv) = MAX(zn0s(iv),1e6_ireals)
      ELSE
        ! Old constant n0s
        zn0s(iv) = 8.0e5_ireals
      ENDIF
      zcrim (iv) = ccsrim*zn0s(iv)
      zcagg (iv) = ccsagg*zn0s(iv)
      zbsdep(iv) = ccsdep*SQRT(v0snow)
      zvz0s (iv) = ccsvel*EXP(ccsvxp * LOG(zn0s(iv)))

      IF (zvzs(iv) == 0.0_ireals) THEN
        zvzs(iv) = zvz0s(iv) * EXP (ccswxp * LOG (0.5_ireals*zqsk(iv))) * zrho1o2(iv)
      ENDIF
    ENDDO loop_over_qs_prepare

!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qr_sedi: DO i1d = 1, ic2
      iv = ivdx2(i1d)

      IF (zvzr(iv) == 0.0_ireals) THEN
!replaced: zvzr(iv) = zvz0r * EXP (x1o8  * LOG (0.5_ireals*zqrk(iv))) * zrho1o2(iv)
        zvzr(iv) = zvz0r * EXP (zvzxp  * LOG (0.5_ireals*zqrk(iv))) * zrho1o2(iv)

      ENDIF
    ENDDO loop_over_qr_sedi

!CK<    
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qi_sedi: DO i1d = 1, ic4
      iv = ivdx4(i1d)

      qig  = qi(iv,k)
      qnig = qni(iv,k)
      rhog = rho(iv,k)
      
      zxi(iv) = MAX(MIN(qig/(qnig+zeps),zximax),zximin)
      
      IF (zvzi(iv) == 0.0_ireals) THEN
        ! Cloud ice sedimentation based on Khvorostyanov and Curry 2005 
        cbest  = 2.0_ireals * kc_alpha * g / ( kc_gamma * rhog * con_m**2 )
        
        ! For consistency, hexagonal plates are assumed
        zDi(iv) = EXP((1./kc_beta)*LOG(zxi(iv)/kc_alpha))
        
        ! Best number, Eq. (8) of MH05
        Xbest = cbest * EXP(vtxexp*LOG(zDi(iv)))
        bracket= SQRT(1.0_ireals + kc_c1*SQRT(Xbest)) - 1.0_ireals
        
        ! Eq (2.8) of KC05
        kc_b1 = kc_c1*SQRT(Xbest) / ( 2.0*bracket * SQRT( 1.0 + kc_c1*SQRT(Xbest)) )
      
        ! Eq (2.7) of KC05
        kc_b2 =EXP( kc_b1*LOG(Xbest))
        kc_a1 = kc_c2 * bracket**2 / kc_b2
        
        ! Eq. (2.11)-(2.14) of KC05
        zvzi(iv) = kc_a1 * con_m **(1.0_ireals -2.0_ireals *kc_b1)    &
            * ( ( 2.0_ireals * kc_alpha * g ) / ( EXP( kc_b1 *LOG(rhog * kc_gamma) ) ))   &
            * EXP( (kc_b1*vtxexp - 1.0_ireals)*LOG (zDi(iv)))

!       zvzi(iv)  = MIN(MAX(vt_kc05(iv)),zvimin),zvimax)
      END IF
    END DO loop_over_qi_sedi
!CK<


  !----------------------------------------------------------------------------
  ! Section 3:
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
    zeln7o8qrk   (:) = 0.0_ireals
    zeln27o16qrk (:) = 0.0_ireals
    zeln13o8qrk  (:) = 0.0_ireals
!!    zeln3o16qrk  (:) = 0.0_ireals
    zeln13o12qsk (:) = 0.0_ireals
    zeln5o24qsk  (:) = 0.0_ireals
    zeln2o3qsk   (:) = 0.0_ireals
    zcsdep       (:) = 3.2E-2_ireals
    zcidep       (:) = 1.3E-5_ireals
    zcslam       (:) = 1e10_ireals

    scau         (:) = 0.0_ireals
    scac         (:) = 0.0_ireals
    !CK>
    snucn        (:) = 0.0_ireals
    snucq        (:) = 0.0_ireals
    snucn2       (:) = 0.0_ireals
    snucq2       (:) = 0.0_ireals
    !CK<    
    scfrz        (:) = 0.0_ireals
    simelt       (:) = 0.0_ireals
    sidep        (:) = 0.0_ireals
    ssdep        (:) = 0.0_ireals
    srim         (:) = 0.0_ireals
    sshed        (:) = 0.0_ireals
    sicri        (:) = 0.0_ireals
    srcri        (:) = 0.0_ireals
    sagg         (:) = 0.0_ireals
    siau         (:) = 0.0_ireals
    ssmelt       (:) = 0.0_ireals
    sev          (:) = 0.0_ireals
    srfrz        (:) = 0.0_ireals
!CDIR END

    ic1 = 0
    ic2 = 0
    ic3 = 0
    ic4 = 0
    ic5 = 0
    ic6 = 0
    ic7 = 0 !CK<   
    DO iv = iv_start, iv_end

        qrg  = qr(iv,k)
        qsg  = qs(iv,k)
        qvg  = qv(iv,k)
        qcg  = qc(iv,k)
        qig  = qi(iv,k)
        tg   =  t(iv,k)
        wg   =  w(iv,k) !CK<        
        ppg  =  p(iv,k)
        rhog = rho(iv,k)

        llqr = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin
        llqi = zqik(iv) > zqmin !CK<  

        !CK>
        ! saturation vapor pressure, saturation mixing ratio and saturation ratio over ice
        zpvsi(iv) = sat_pres_ice(tg)
        zqvsi(iv) = zpvsi(iv)/(rhog * r_v * tg) 
        zsi(iv)   = qvg * rhog * tg * r_v / zpvsi(iv)
        scr       = 2.349 - tg / 259.00 !critical ice supersaturation
        !CK<
        
        IF (llqr) THEN
          !US reported by Thorsten Reinhardt: Multiplication with zrho1o2 was missing
          zpkr(iv) = zqrk(iv) * zvz0r * EXP (zvzxp * LOG (zqrk(iv))) * zrho1o2(iv)
        ELSE
          zpkr(iv) = 0.0_ireals
        ENDIF

        IF (llqs) THEN
          !US reported by Thorsten Reinhardt: Multiplication with zrho1o2 was missing
          zpks(iv) = zqsk (iv) * zvz0s(iv) * EXP (ccswxp * LOG (zqsk(iv))) * zrho1o2(iv)
        ELSE
          zpks(iv) = 0.0_ireals
        ENDIF

        !CK>
        IF (llqi) THEN
          zpki(iv)  = zqik(iv)  * zvzi(iv)
          zpkni(iv) = zqnik(iv) * zvzi(iv)
        ELSE
          zpki(iv)  = 0.0_ireals
          zpkni(iv) = 0.0_ireals
        ENDIF
      !CK<
      
        zpkr(iv)   = MIN( zpkr(iv) , zzar(iv) )
        zpks(iv)   = MIN( zpks(iv) , zzas(iv) )
        zpki(iv)   = MIN( zpki(iv) , zzai(iv) ) 
        zpkni(iv)  = MIN( zpkni(iv), zzani(iv))

        zzar(iv)   = zdtdh(iv) * (zzar(iv)-zpkr(iv))
        zzas(iv)   = zdtdh(iv) * (zzas(iv)-zpks(iv))
        zzai(iv)   = zdtdh(iv) * (zzai(iv)-zpki(iv))
        zzani(iv)  = zdtdh(iv) * (zzani(iv)- zpkni(iv))
       
        zimr(iv)   = 1.0_ireals / (1.0_ireals + zvzr(iv) * zdtdh(iv))
        zims(iv)   = 1.0_ireals / (1.0_ireals + zvzs(iv) * zdtdh(iv))
        zimi(iv)   = 1.0_ireals / (1.0_ireals + zvzi(iv) * zdtdh(iv)) !
        zimni(iv)  = 1.0_ireals / (1.0_ireals + zvzi(iv) * zdtdh(iv)) !not a typo

        zqrk(iv)   = zzar(iv)*zimr(iv)
        zqsk(iv)   = zzas(iv)*zims(iv)
        zqik(iv)   = zzai(iv)*zimi(iv)  
        zqnik(iv)  = zzani(iv)*zimni(iv)

        llqr = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin
        llqc =       qcg > zqmin
        llqi =       qig > zqmin

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
        !CK>
        ! heterogeneous nucleation requires supersaturation
        IF ( (tg < t0 ) .AND. (zsi(iv) > 1.0_ireals) .AND. (zsi(iv) < scr) &
          .AND. (qni_nuc(iv,k) < ni_het_max) ) THEN
          ic4 = ic4 + 1
          ivdx4(ic4) = iv
        ENDIF
        ! homogeneous nucleation requires critical supersaturation
        IF ( (tg < t0 ) .AND. (zsi(iv) > scr) .AND. &
          (qni(iv,k) < ni_hom_max) .AND. (wg > 0._ireals) ) THEN
          ic7 = ic7 + 1
          ivdx7(ic7) = iv
        ENDIF
        !IF ( tg < zthet .AND. qvg >  8.E-6_ireals &
        !                .AND. qig <= 0.0_ireals ) THEN
        !  ic4 = ic4 + 1
        !  ivdx4(ic4) = iv
        !ENDIF
        !CK<
        IF (llqc) THEN
          ic5 = ic5 + 1
          ivdx5(ic5) = iv
        ENDIF
        IF (llqr .AND. qcg <= 0.0_ireals) THEN
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
        zeln27o16qrk(iv) = EXP (x27o16 * zlnqrk)
      ENDIF
      IF (llqi) THEN
        zeln13o8qrk(iv)  = EXP (x13o8  * zlnqrk)
      ENDIF
!      IF (qcg <= 0.0_ireals ) THEN
!        zeln3o16qrk(iv)  = EXP (x3o16  * zlnqrk)
!      ENDIF
    ENDDO loop_over_qr


! ic2
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs_coeffs: DO i1d =1, ic2
      iv = ivdx2(i1d)

      qcg = qc(iv,k)
      qig = qi(iv,k)

      zlnqsk       = LOG (zqsk(iv))
      IF (qig+qcg > zqmin) THEN
        zeln13o12qsk(iv) = EXP (x13o12 * zlnqsk)
      ENDIF
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
      zdvtp  = ccdvtp * EXP(1.94_ireals * LOG(tg)) / ppg          
      zhi    = ccshi1*zdvtp*rhog*zqvsi(iv)/(tg*tg)
      hlp    = zdvtp / (1.0_ireals + zhi)
      zcidep(iv) = ccidep * hlp
      IF (llqs) THEN
        zcslam(iv) = EXP(ccslxp * LOG(ccslam * zn0s(iv) / zqsk(iv) ))
        zcslam(iv) = MIN(zcslam(iv),1e15_ireals)
        zcsdep(iv) = 4.0_ireals * zn0s(iv) * hlp
      ENDIF

!>>>

    ENDDO loop_over_qi_qs

    !--------------------------------------------------------------------------
    ! Section 4: Initialize the conversion rates with zeros in every layer
    !            Deposition nucleation for low temperatures
    !--------------------------------------------------------------------------

    ! Deposition nucleation for low temperatures below a threshold

! ic4
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_heticenucleation: DO i1d = 1, ic4
      iv = ivdx4(i1d)

      ssi  = zsi(iv)       ! saturation ratio
      ppg  = p(iv,k)
      tg   = t(iv,k)
      qnig = qni(iv,k)
      qcg  = qc(iv,k)
      qvg  = qv(iv,k)

      qnig_nuc = qni_nuc(iv,k)

      llqc = qcg > zqmin

      ! Heterogeneous ice nucleation using Phillips et al (2008) scheme
      ! calculate indices used for look-up tables
      xt = (274._ireals - tg)  / ttstep
      xs = (100._ireals*(ssi-1.0_ireals)) / ssstep    
      xt = MIN(xt,ttmax-1.0_ireals)
      xs = MIN(xs,ssmax-1.0_ireals)
      tt = INT(xt)
      ss = INT(xs)
      
      IF (llqc) THEN
        ! immersion freezing at water saturation as function of temperature
        infrac(1) = (tt+1-xt) * afrac_dust(tt,99) + (xt-tt) * afrac_dust(tt+1,99)
        infrac(2) = (tt+1-xt) * afrac_soot(tt,99) + (xt-tt) * afrac_soot(tt+1,99)
        infrac(3) = (tt+1-xt) * afrac_orga(tt,99) + (xt-tt) * afrac_orga(tt+1,99)
      ELSE
        ! deposition freezing, bi-linear interpolation in look-up tables
        infrac(1) = (tt+1-xt) * (ss+1-xs) * afrac_dust(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_dust(tt+1,ss  ) &
                    + (tt+1-xt) * (xs-ss)   * afrac_dust(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_dust(tt+1,ss+1)
        infrac(2) = (tt+1-xt) * (ss+1-xs) * afrac_soot(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_soot(tt+1,ss  ) &
                    + (tt+1-xt) * (xs-ss)   * afrac_soot(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_soot(tt+1,ss+1)
        infrac(3) = (tt+1-xt) * (ss+1-xs) * afrac_orga(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_orga(tt+1,ss  ) &
                    + (tt+1-xt) * (xs-ss)   * afrac_orga(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_orga(tt+1,ss+1)
      END IF
 
      znidiag(iv) = na_dust * infrac(1) + na_soot * infrac(2) + na_orga * infrac(3)
      znidiag(iv) = MIN(znidiag(iv),ni_het_max)
      
      IF (znidiag(iv).LT.1e-9) znidiag(iv) = 0.0_ireals
       
      snucn(iv) = MAX( z1orhog(iv)*znidiag(iv) - qnig_nuc,0.d0) * zdtr
      snucq(iv) = zmi0 * snucn(iv)

     ! IF (znidiag(iv).GT.0) THEN
     !   WRITE(*,*)'ni_het,qnig_nuc,ssi',znidiag*1e-6, qnig_nuc*1e-6,ssi
     ! END IF
        
    ENDDO loop_over_heticenucleation

!!!!CDIR noNODEP,VOVERTAKE,VOB
!CDIR NOMOVE,NODEP
    loop_over_homicenucleation: DO i1d = 1, ic7
        iv = ivdx7(i1d)

        ! Homogeneous nucleation using KHL06 scheme
        ssi  = zsi(iv)       ! saturation ratio
        ppg  = p(iv,k)
        tg   = t(iv,k)
        wg   = w(iv,k) + 0.5 * tke(iv,k)
        qig  = qi(iv,k)
        qnig = qni(iv,k)
        qcg  = qc(iv,k)
        qvg  = qv(iv,k)
        rhog = rho(iv,k)

        zxi(iv) = MAX(MIN(qig/(qnig+zeps),zximax),zximin)
        
        ! critical supersaturation for homogenous nucleation, cf. Eq. (1) of KB08
        scr = 2.349_ireals - tg / 259.00_ireals
        
        ! pre-existing ice is spherical or not?
        zri = EXP( x1o3 * LOG(zxi(iv)/asphere) )

        zdvtp   = ccdvtp * EXP(1.94_ireals * LOG(tg)) / ppg          
        v_th    = SQRT( 8._ireals*k_b*tg/(pi*ma_w) ) 
        flux    = alpha_d * v_th/4._ireals
        n_sat   = zpvsi(iv) / (k_b*tg)  
        
        ! coeffs of supersaturation equation
        acoeff(1) = (lh_s * g) / (cp_d * r_v * tg**2) - g/(r_d * tg)
        acoeff(2) = 1.0_ireals / n_sat
        acoeff(3) = (lh_s**2 * M_w * ma_w)/(cp_d * ppg * tg * M_a) 
        
        ! coeffs of depositional growth equation
        bcoeff(1) = flux * svol * n_sat * (ssi - 1.0_ireals)
        bcoeff(2) = flux / zdvtp
        
        ! pre-existing ice crystals included as reduced updraft speed
        ri_dot = bcoeff(1) / (1.0_ireals + bcoeff(2) * zri)
        R_ik   = 4.0_ireals * pi / svol * rhog * qnig * zri**2 * ri_dot
        w_p    = (acoeff(2) + acoeff(3) * ssi)/(acoeff(1) * ssi) * R_ik  ! KHL06, Eq.(19)
        
        cool   = g / cp_d * wg
        ctau   = tg * ( 0.004_ireals *tg - 2._ireals ) + 304.4_ireals         
        tau    = 1.0_ireals / (ctau * cool)                    ! freezing timescale, Eq.(5)

        si_het =  ssi + acoeff(1) * ssi * (wg - w_p) * zdt     ! Si for next time step KHL06 Eq.(20)

        lhom = ( wg > w_p )

        IF (lhom .AND. (si_het > scr) ) THEN 
          ! do not trigger hom. nucleation if next time step is still within nucleation event          
          thom = (1.0_ireals/(acoeff(1)*wg) * LOG(scr/ssi) + tau)          
          lhom = ( thom < (0.5_ireals * zdt) ) ! nucleation event in less than half a time step
        ELSE
          lhom = .FALSE.
        END IF
 
        IF (lhom) THEN   ! homogeneous nucleation event
          
          ! timescales of freezing event (see KL02, RM05, KHL06)
          delta   = (bcoeff(2) * r_0)                         ! dimless aerosol radius, eq.(4)  
          phi     = acoeff(1)*ssi / ( acoeff(2) + acoeff(3)*ssi) * (wg - w_p) 
          
          ! monodisperse approximation following KHL06
          kappa   = 2._ireals * bcoeff(1) * bcoeff(2) * tau / (1._ireals+ delta)**2  ! kappa, Eq. 8 KHL06
          sqrtkap = SQRT(kappa)                                        ! root of kappa
          ren     = 3._ireals * sqrtkap / ( 2.0_ireals + SQRT(1._ireals+9._ireals*kappa/pi) ) ! analy. approx. of erfc by RM05

          R_imfc  = 4. * pi * bcoeff(1)/bcoeff(2)**2 / svol  
          R_im    = R_imfc / (1.0_ireals+ delta) * ( delta**2 - 1.0_ireals &
            & + (1.0_ireals+0.5_ireals*kappa*(1.0_ireals+ delta)**2) * ren/sqrtkap+ zeps)           ! Eq. 6 KHL06
          
          ! number concentration and radius of ice particles after nucleation
          ni_hom  = phi / R_im                                         ! Eq.9 KHL06
          ri_0    = 1.0_ireals + 0.5_ireals * sqrtkap * ren                           ! for Eq. 3 KHL06 
          ri_hom  = (ri_0 * (1.0_ireals + delta) - 1._ireals ) / bcoeff(2)            ! Eq. 3 KHL06 * REN = Eq.23 KHL06

          ! calculate change of ri by depositional growth
          ri_max = ( ri_hom**3 + 3._ireals * ma_w * n_sat * (scr - 1._ireals)  &     ! Eq. 12 KL02
                 / (asphere * ni_hom) )**x1o3     
         
          ! iteration 
          tgrow   = 0.75_ireals/(pi * ni_hom * zdv * ri_max)
          beta    = (4.0_ireals*zdv)/(v_th * alpha_d * ri_max) ! dimless, for Eq.15 KL02
          Xj      = ri_hom/ri_max                     ! ratio between initial particle size and final ri_max
          dXj     = 1.0_ireals - Xj
          
!CDIR EXPAND=8
          DO jj=1,8
            dXj     = dXj * 0.5_ireals
            Xmid    = Xj + dXj
            Fmid    = dep_growth_timescale(beta,Xmid,Xj) - zdt/tgrow
            IF (Fmid .LT. 0.) Xj = Xmid
          ENDDO

          ! Values after dt
          ri_hom = Xj * ri_max
          mi_hom = asphere * ri_hom**3

          ! ni_hom is a number density, need #/(kg*s)
          snucn2(iv) = MIN(MAX(z1orhog(iv)*ni_hom, 0.d0),ni_hom_max) * zdtr  
          snucq2(iv) = mi_hom * snucn2(iv)

         ! IF (ni_hom .GT. 0) THEN
         !   WRITE(*,*)'ni_hom,qig,ssi',ni_hom*1e-6, qig*1e6,ssi
         !   !snucq2(iv) = MIN(mi_hom * snucn2(iv),qvg - zqvsi)
         ! END IF
          
        ENDIF
 
      ENDDO loop_over_homicenucleation

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
      qnig =  qni(iv,k) !CK<    
      tg   =    t(iv,k)
      ppg  =    p(iv,k)
      rhog =  rho(iv,k)
      llqs = zqsk(iv) > zqmin

      IF (iautocon == 0) THEN
        ! Kessler (1969) autoconversion rate
        zscau  = zccau * MAX( qcg - qc0, 0.0_ireals )
        zscac  = zcac  * qcg * zeln7o8qrk(iv)
      ELSEIF (iautocon == 1) THEN
        ! Seifert and Beheng (2001) autoconversion rate
        ! with constant cloud droplet number concentration cloud_num
        IF (qcg > 1e-6) THEN
          ztau   = MIN(1.0_ireals-qcg/(qcg+qrg),0.9_ireals) 
          zphi   = zkphi1 * ztau**zkphi2 * (1.0_ireals - ztau**zkphi2)**3
          zscau  = zconst * qcg*qcg*qcg*qcg &
            &    * (1.0_ireals + zphi/(1.0_ireals - ztau)**2)
          zphi   = (ztau/(ztau+zkphi3))**4
          zscac  = zkcac * qcg * qrg * zphi !* zrho1o2(iv)
        ELSE
          zscau  = 0.0_ireals
          zscac  = 0.0_ireals
        ENDIF
      ENDIF
      IF (llqs) THEN
        zscrim = zcrim(iv) * EXP(ccsaxp * LOG(zcslam(iv))) * qcg !* zrho1o2(iv)
      ELSE
        zscrim = 0.0_ireals
      ENDIF

      zscshe = 0.0_ireals
      IF( tg >= t0 ) THEN
        zscshe = zscrim
        zscrim = 0.0_ireals
      ENDIF
      ! Check for maximum depletion of cloud water and adjust the
      ! transfer rates accordingly
      zscmax = qcg*zdtr
      !CK>
      ! heterogeneous freezing using Phillips et al. (2008) scheme
      scfrz(iv) = zxcloud * MAX( z1orhog(iv)*znidiag(iv) - qnig,0.d0) * zdtr
      zscsum = zscau + zscac + zscrim + zscshe + scfrz(iv)
      !CK<      
      !zscsum = zscau + zscac + zscrim + zscshe 
      zcorr  = zscmax / MAX( zscmax, zscsum )
      IF( tg <= zthn ) THEN
        scfrz(iv) = zscmax
      ELSE
        scau (iv) = zcorr*zscau
        scac (iv) = zcorr*zscac
        srim (iv) = zcorr*zscrim
        sshed(iv) = zcorr*zscshe 
      ENDIF

      !CK>
      ! Calculation of heterogeneous nucleation of cloud ice.
      ! This is done in this section, because we require water saturation
      ! for this process (i.e. the existence of cloud water) to exist.
      ! Hetrogeneous nucleation is assumed to occur only when no
      ! cloud ice is present and the temperature is below a nucleation
      ! threshold.
      ! IF( tg <= 267.15_ireals .AND. qig <= 0.0_ireals ) THEN
      !  znin  = MIN( fxna(tg), znimax )
      !  zsnuc = zmi0 * z1orhog(iv) * znin * zdtr
      !  snuc(iv) = zsnuc
      !ENDIF
      !CK<
      ! Calculation of in-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        zsrfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 )* zeln27o16qrk(iv)
        srfrz(iv) = zsrfrz
      ENDIF

    ENDDO loop_over_qc

      !------------------------------------------------------------------------
      ! Section 6: Search for cold grid points with cloud ice and/or snow and
      !            calculation of the conversion rates involving qi and ps
      !------------------------------------------------------------------------

! also ic3
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs: DO i1d =1, ic3
      iv = ivdx3(i1d)

      qvg  =  qv(iv,k)
      qig  =  qi(iv,k)
      qnig =  qni(iv,k)
      qsg  =  qs(iv,k)
      tg   =   t(iv,k)
      ppg  =   p(iv,k)
      rhog = rho(iv,k)
      llqi =  qig > zqmin
      !CK>
      zqvsidiff = qvg - zqvsi(iv)
      tau_dep_ice_inv  = 0.0_ireals
      tau_dep_snow_inv = 0.0_ireals
      !CK<
      
      IF (tg<=t0) THEN
        !CK>
        znid    = rhog * qnig 
        zmi     = MAX(MIN(qig/(qnig+zeps), zmimax ), zmi0)
        zsvmax  = (qvg - zqvsi(iv)) * zdtr
        zsagg   = zcagg(iv) * EXP(ccsaxp*LOG(zcslam(iv))) * qig
        zsagg   = MAX( zsagg, 0.0_ireals ) & !* zrho1o2(iv) &
          * MAX(0.2_ireals,MIN(EXP(0.09_ireals*(tg-t0)),1.0_ireals))
        
        IF (llqi) THEN
          ! depositional growth as relaxation for cloud ice
          tau_dep_ice_inv = zcidep(iv) * znid * EXP(0.33_ireals *  LOG (zmi))
!          WRITE (*,*) 'ni/cm3, qi mg/kg, 1/tau_i s ', znid*1e-6,qig*1e6,1./tau_dep_ice_inv
        END IF
          
        ! depositional growth as relaxation for snow
        zxfac      = 1.0_ireals + zbsdep(iv) * EXP(ccsdxp*LOG(zcslam(iv)))
        tau_dep_snow_inv = zcsdep(iv) * zxfac/ zcslam(iv)**2
        
        Xi_mck_inv = (tau_dep_ice_inv + tau_dep_snow_inv)!Eq. 7
        Xi_mck     = 1.0_ireals/ Xi_mck_inv

        IF (llqi) THEN
          zsidep  = tau_dep_ice_inv*Xi_mck*zqvsidiff * (1.0_ireals - EXP(-zdt*Xi_mck_inv)) * zdtr
          zsiau = zciau * MAX( qig - qi0, 0.0_ireals ) &
                            * MAX(0.2_ireals,MIN(EXP(0.09_ireals*(tg-t0)),1.0_ireals))

          zsicri    = zcicri * qig * zeln7o8qrk(iv)
          zsrcri    = zcrcri * (qig/zmi) * zeln13o8qrk(iv) 
        ELSE
          zsidep = 0.0_ireals
          zsiau = 0.0_ireals
          zsicri = 0.0_ireals
          zsrcri = 0.0_ireals
        ENDIF
          
        zssdep  = tau_dep_snow_inv*Xi_mck*zqvsidiff * (1.0_ireals - EXP(-zdt*Xi_mck_inv)) * zdtr
      
        zsvidep   = 0.0_ireals
        zsvisub   = 0.0_ireals
        zsimax    = qig*zdtr

        IF( zsidep > 0.0_ireals ) THEN
          zsvidep = MIN( zsidep, zsvmax )
        ELSEIF (zsidep < 0.0_ireals ) THEN
  !        zsvisub = - MAX(-zsimax, zsvmax ) ! this would imply sublimating everything within one time step
          zsvisub = - MAX(-zsidep, zsvmax )
        ENDIF

       
!        zxfac     = 1.0_ireals + zbsdep(iv) * EXP(ccsdxp*LOG(zcslam(iv)))
!        zssdep    = zcsdep(iv) * zxfac * ( qvg - zqvsi ) / (zcslam(iv)+zeps)**2

!        ! Check for maximal depletion of vapor by sdep
!        IF (zssdep > 0.0_ireals) zssdep = MIN(zssdep, zsvmax-zsvidep)
!        ! Check for maximal depletion of snow by sdep
!        IF (zssdep < 0.0_ireals) zssdep = MAX(zssdep, -qsg*zdtr)

        zsisum = zsiau + zsagg + zsicri + zsvisub
        zcorr  = 0.0_ireals
        IF( zsimax > 0.0_ireals ) zcorr  = zsimax / MAX( zsimax, zsisum )
        sidep(iv)  = zsvidep - zcorr*zsvisub
        siau (iv)  = zcorr*zsiau
        sagg (iv)  = zcorr*zsagg
        ssdep(iv)  = zssdep
        srcri(iv)  = zsrcri
        sicri(iv)  = zcorr*zsicri
!CK<        

      !------------------------------------------------------------------------
      ! Section 7: Search for warm grid points with cloud ice and/or snow and
      !            calculation of the melting rates of qi and ps
      !------------------------------------------------------------------------

      ELSE ! tg > 0
        simelt(iv)  = qig*zdtr
        zqvsw0      = zpvsw0 / (rhog * r_v * tg)
        zx1         = (tg - t0) + zasmel*(qvg - zqvsw0)
        zx2         = 1.0_ireals + zbsmel * zeln5o24qsk(iv)
        zssmelt     = zcsmel * zx1 * zx2 * zeln2o3qsk(iv)
        ssmelt(iv) = MAX( zssmelt, 0.0_ireals )
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
      zqvsw    = sat_pres_water(tg)/(rhog * r_v *tg)

      zlnqrk      = LOG (zqrk(iv))
      zx1         = 1.0_ireals + zbev * EXP (zbevxp  * zlnqrk)
      sev(iv)    = zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk)
!      zqvsw    = fqvs( fpvsw(tg), ppg )
!      zx1      = 1.0_ireals + zbev* zeln3o16qrk(iv)
!      zsev     = zcev*zx1*(zqvsw - qvg)*SQRT(zqrk(iv))
!      sev(iv) = MAX( zsev, 0.0_ireals )

      ! Calculation of below-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        zsrfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 ) * zeln27o16qrk(iv)
        srfrz(iv)  = zsrfrz
      ENDIF

    ENDDO loop_over_qr_nocloud

    !--------------------------------------------------------------------------
    ! Section 9: Calculate the total tendencies of the prognostic variables.
    !            Update the prognostic variables in the interior domain.
    !--------------------------------------------------------------------------

    loop_over_all_iv: DO iv = iv_start, iv_end

!        qvg = qv(iv,k)
!        qcg = qc(iv,k)
        qrg  = qr(iv,k)
        qsg  = qs(iv,k)
        qig  = qi(iv,k)
        qnig = qni(iv,k)         !CK<
        qnig_nuc = qni_nuc(iv,k) !CK<         
!        tg  = t (iv,k)
        rhog = rho(iv,k)

        zxi(iv) = MAX(MIN(qig/(qnig+zeps),zximax),zximin) !CK<
        
        zsrmax = zzar(iv)*z1orhog(iv)*zdtr
        zssmax = zzas(iv)*z1orhog(iv)*zdtr
        zsrsum = sev(iv) + srfrz(iv) + srcri(iv)
        zcorr  = 1.0_ireals
        IF(zsrsum > 0) THEN
          zcorr  = zsrmax / MAX( zsrmax, zsrsum )
        ENDIF
        sev  (iv) = zcorr*sev(iv)
        srfrz(iv) = zcorr*srfrz(iv)
        srcri(iv) = zcorr*srcri(iv)
        ssmelt(iv) = MIN(ssmelt(iv), zssmax)
        IF (ssdep(iv) < 0.0_ireals ) THEN
          ssdep(iv) = MAX(ssdep(iv), - zssmax)
        ENDIF
!CK>        
        zqvt = sev(iv)   - sidep(iv) - ssdep(iv)  - snucq(iv)  - snucq2(iv)
        zqct = simelt(iv)- scau(iv)  - scfrz(iv)  - scac(iv)   - sshed(iv) - srim(iv) 
        zqit = snucq(iv)  + snucq2(iv) + scfrz(iv) - simelt(iv) - sicri(iv)  + sidep(iv) &
                                                   - sagg(iv)   - siau(iv)
        zqnit = snucn (iv) + snucn2(iv) + scfrz (iv) / zxi(iv)  &
            - ( simelt(iv) + sicri (iv) + sagg  (iv) + siau(iv) ) / zxcloud      
!CK<        
        zqrt = scau(iv)  + sshed(iv) + scac(iv)   + ssmelt(iv) - sev(iv) - srcri(iv) - srfrz(iv) 
        zqst = siau(iv)  + sagg(iv)  - ssmelt(iv) + sicri(iv) + srcri(iv)            & 
                                                  + srim(iv) + ssdep(iv) + srfrz(iv)
        ztt = z_heat_cap_r*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst) )

        ! Update variables and add qi to qrs for water loading
!CK>
        qig  = MAX ( 0.0_ireals, (zzai (iv)*z1orhog(iv) + zqit  * zdt) * zimi(iv))
        qnig = MAX ( 0.0_ireals, (zzani(iv)*z1orhog(iv) + zqnit * zdt) * zimni(iv))
        !qig  = MAX ( 0.0_ireals, qig + zqit*zdt)
        
        ! tracking of activated heterogeneous IN
        qni_nuc(iv,k) = MAX(0.0_ireals, qnig_nuc + ( snucn(iv) - qnig_nuc/tau_mix ) * zdt)
!CK<
        
        qrg = MAX ( 0.0_ireals, (zzar(iv)*z1orhog(iv) + zqrt*zdt)*zimr(iv))
        qsg = MAX ( 0.0_ireals, (zzas(iv)*z1orhog(iv) + zqst*zdt)*zims(iv))

        !----------------------------------------------------------------------
        ! Section 10: Complete time step
        !----------------------------------------------------------------------

        IF ( k /= ke) THEN
          ! Store precipitation fluxes and sedimentation velocities 
          ! for the next level
          zprvr(iv) = qrg*rhog*zvzr(iv)
          zprvs(iv) = qsg*rhog*zvzs(iv)
          zprvi(iv) = qig*rhog*zvzi(iv)
          zprvni(iv) = qnig*rhog*zvzi(iv)

          IF (zprvr(iv) <= zqmin) zprvr(iv)=0.0_ireals
          IF (zprvs(iv) <= zqmin) zprvs(iv)=0.0_ireals
          !CK<
          IF (zprvi(iv) <= zqmin .or. zprvni(iv) <= zqmin) THEN
             zprvi(iv)  = 0.0_ireals
            zprvni(iv)  = 0.0_ireals
          END IF
          !CK>
        
          ! for the latent heat nudging
          IF (ldass_lhn) THEN
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)
            qrsflux(iv,k) = 0.5*(qrsflux(iv,k)+zpkr(iv)+zpks(iv))
          ENDIF

          IF (qrg+qr(iv,k+1) <= zqmin) THEN
            zvzr(iv)= 0.0_ireals
          ELSE
            zvzr(iv)= zvz0r * EXP(zvzxp*LOG((qrg+qr(iv,k+1))*0.5_ireals*rhog)) * zrho1o2(iv)
          ENDIF
          IF (qsg+qs(iv,k+1) <= zqmin) THEN
            zvzs(iv)= 0.0_ireals
          ELSE
            zvzs(iv)= zvz0s(iv) * EXP(zv1s/(zbms+1.0_ireals)*LOG((qsg+qs(iv,k+1))*0.5_ireals*rhog)) * zrho1o2(iv)
          ENDIF
!CK>
          IF (qig+qi(iv,k+1) <= zqmin .OR. qnig+qni(iv,k+1)<= zqmin) THEN
            zvzi(iv)= 0.0_ireals
          ELSE 
            cbest  = 2.0_ireals * kc_alpha * g / ( kc_gamma * rhog * con_m**2 )
            
            hlp =(qig+qi(iv,k+1))/(qnig+qni(iv,k+1)+zeps)
            zDi(iv) = EXP((1./kc_beta)*LOG(hlp/kc_alpha))
            
            ! Best number, Eq. (8) of MH05
            Xbest = cbest * EXP(vtxexp*LOG(zDi(iv)))
            
            bracket= SQRT(1.0_ireals + kc_c1*SQRT(Xbest)) - 1.0_ireals
            
            ! Eq (2.8) of KC05
            kc_b1 = kc_c1*SQRT(Xbest) / ( 2.0*bracket * SQRT( 1.0 + kc_c1*SQRT(Xbest)) )
            
            ! Eq (2.7) of KC05
            kc_b2 =EXP( kc_b1*LOG(Xbest))
            kc_a1 = kc_c2 * bracket**2 / kc_b2
            
            ! This is Eq. (2.11)-(2.14) of KC05
            zvzi(iv) = kc_a1 * con_m **(1.0_ireals -2.0_ireals *kc_b1)    &
              * ( ( 2.0_ireals * kc_alpha * g ) / ( EXP( kc_b1 *LOG(rhog * kc_gamma) ) ))   &
              * EXP( (kc_b1*vtxexp - 1.0_ireals)*LOG (zDi(iv)))
            
            ! zvzi(iv) = MIN(MAX(zvzi(iv),0.0_ireals),zvimax)
            
          ENDIF
          !CK<
          
        ELSE
          ! Precipitation fluxes at the ground 
          prr_gsp(iv) = 0.5_ireals * (qrg*rhog*zvzr(iv) + zpkr(iv))
!CK>
          prs_gsp(iv) = 0.5_ireals * (qsg*rhog*zvzs(iv) + zpks(iv)) + & 
                        0.5_ireals * (qig*rhog*zvzi(iv) + zpki(iv))     
!CK<        

          ! for the latent heat nudging
          IF (ldass_lhn) &
            qrsflux(iv,k) = prr_gsp(iv)+prs_gsp(iv)

        ENDIF

        ! Update of prognostic variables or tendencies
        qr (iv,k) = qrg
        qs (iv,k) = MAX ( 0.0_ireals, qsg )
        qi (iv,k) = qig
        qni(iv,k) = qnig !CK< 
!        qrs(iv,k   ) = qrg+qsg+qig
        t  (iv,k) = t (iv,k) + ztt*zdt 
        qv (iv,k) = MAX ( 0.0_ireals, qv(iv,k) + zqvt*zdt )
        qc (iv,k) = MAX ( 0.0_ireals, qc(iv,k) + zqct*zdt )

      ENDDO loop_over_all_iv

  IF (izdebug > 15) THEN
    ! Check for negative values
     DO iv = iv_start, iv_end
        IF (qr(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_ice, negative value in qr'
          CALL message('',message_text)
        ENDIF
        IF (qc(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_ice, negative value in qc'
          CALL message('',message_text)
        ENDIF
        IF (qi(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_ice, negative value in qi'
          CALL message('',message_text)
        ENDIF
        !CK>
        IF (qni(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_ice, negative value in qni'
          CALL message('',message_text)
        ENDIF
        !CK<
        IF (qs(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_ice, negative value in qs'
          CALL message('',message_text)
        ENDIF
        IF (qv(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_ice, negative value in qv'
          CALL message('',message_text)
        ENDIF
    ENDDO
  ENDIF

#if defined (__COSMO__)
  ! Do a final saturation adjustment for new values of t, qv and qc
!CDIR COLLAPSE
    zpres(:) = p(:,k)

    CALL satad ( 1, t(:,k), qv(:,k),              &
               qc(:,k), t(:,k), zpres,          &
               zdummy(:,1),zdummy(:,2),zdummy(:,3), &
               zdummy(:,4),zdummy(:,5),zdummy(:,6), &
               zdummy(:,7),zdummy(:,8),               &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,     &
               rvd_m_o, lh_v, z_heat_cap_r, cp_d,                 &
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

  IF (izdebug > 15) THEN ! for debugging
!  IF (izdebug > 1) THEN
    CALL message('gscp_hydci_pp_ice', 'UPDATED VARIABLES')
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_ice T= ',&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_ice qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_ice qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_ice qi= ',&
    MAXVAL( qi(:,:)), MINVAL(qi(:,:) )
   CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_ice qni= ',&
    MAXVAL( qni(:,:)), MINVAL(qni(:,:) )
   CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_ice qni_nuc= ',&
    MAXVAL( qni_nuc(:,:)), MINVAL(qni_nuc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_ice qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_ice qs= ',&
    MAXVAL( qs(:,:)), MINVAL(qs(:,:) )
    CALL message('', TRIM(message_text))
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine hydci_pp_ice
!------------------------------------------------------------------------------

END SUBROUTINE hydci_pp_ice



!==============================================================================

END MODULE gscp_hydci_pp_ice
