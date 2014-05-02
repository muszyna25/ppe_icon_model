!>
!! cloud microphysics
!!
!! !------------------------------------------------------------------------------
!!
!! @par Description of *init_hydci_pp*:
!!   Initializes some parameters for hydci_pp
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
!! $Id: n/a$
!!
!! @par Revision History
!! implemented into ICON by Felix Rieper (2012-06)
!! 
!! @par Revision History
!! Graupel scheme implemented into ICON by Carmen Köhler (2014-03)
!! Modifications from Felix Rieper for modifications to improve supercooled 
!! liquid water (SLW) prediction in hydci_pp
!!     - reduced number of ice crystal Ni(T), now according to Cooper (1986)
!!     - reduced rain freezing rate srfrzr, now according to Bigg (1953)
!!     - reduced depositional growth of ice and snow (zsidep, zssdep) 
!!       at cloud top, according to R. Forbes (2013) -> IFS model!
!! Unified Model Version for ICON and COSMO by Uli Schaettler and Xavier
!! Lapillone (Meteoswiss)
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

MODULE gscp_hydci_pp

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
!US for SP Version  ireals=>wp,    & !! KIND-type parameter for real variables
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
!   variables for hydci_pp
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
!   variables for hydci_pp and hydci_pp_gr
!
    v0snow,    & ! factor in the terminal velocity for snow
    mu_rain,   & ! 
    rain_n0_factor, &
    cloud_num
#ifdef __ICON__
!XL : there is no hydci_pp_ice in this file could we remove the variable ?
USE data_gscp, ONLY: &  
!
!   variables for hydci_pp_ice
!
    vtxexp,    & !
    kc_c1,     & !
    kc_c2      & !
    kc_alpha   & !
    kc_beta    & !
    kc_gamma   & !
    kc_sigma   & !
    do_i       & !
    co_i
#endif
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

#ifdef NUDGING
USE data_lheat_nudge,           ONLY :  &
    llhn,         & ! main switch for latent heat nudging
    llhnverif,    & ! main switch for latent heat nudging
    lhn_qrs,      & ! use integrated precipitaion flux as reference
    tt_lheat,     & ! latent heat release of the model
    qrsflux         ! total precipitation flux

!------------------------------------------------------------------------------

!USE src_lheating,             ONLY :  & 
!    get_gs_lheating            ! storage of grid scale latent heating for lhn
!this should not be called from block physics
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
USE mo_math_utilities    , ONLY: gamma_fct
USE mo_math_constants    , ONLY: pi
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                                 r_d   => rd    , & !! gas constant for dry air
                                 rvd_m_o=>vtmpc1, & !! r_v/r_d - 1
                                 o_m_rdv        , & !! 1 - r_d/r_v
                                 rdv            , & !! r_d / r_v
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

USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config

USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b2i   => c3ies , & !!               -- " --
                                 b4w   => c4les , & !!               -- " --
                                 b4i   => c4ies , & !!               -- " --
                                 b234w => c5les     !!               -- " --
USE mo_satad,              ONLY: satad_v_3d,     &  !! new saturation adjustment
                                 sat_pres_water, &  !! saturation vapor pressure w.r.t. water
                                 sat_pres_ice!,   &  !! saturation vapor pressure w.r.t. ice
!                                 spec_humi          !! Specific humidity 
USE mo_exception,          ONLY: message, message_text
USE data_gscp                 !xxx: common module COSMO/ICON, all variables are used here
#endif

!==============================================================================

IMPLICIT NONE
PRIVATE

CHARACTER(len=*), PARAMETER :: &
  &  version = '$Id$'

!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: hydci_pp, hydci_pp_gr, hydci_pp_init

LOGICAL, PARAMETER :: &
#ifdef __COSMO__
  lorig_icon   = .FALSE.  , &  ! switch for original ICON setup (only for hydci_pp)
                               ! XL : should be false for COSMO ?
#else
  lorig_icon   = .TRUE.  , &  ! switch for original ICON setup (only for hydci_pp)
#endif
  lsedi_ice    = .FALSE. , &  ! switch for sedimentation of cloud ice (Heymsfield & Donner 1990 *1/3)
  lstickeff    = .FALSE. , &  ! switch for sticking coeff. (work from Guenther Zaengl)
  lsuper_coolw = .FALSE.       ! switch for supercooled liquid water (work from Felix Rieper)
      
REAL(ireals), PARAMETER :: &
  zceff_fac        = 3.5E-3_ireals, & ! Scaling factor [1/K] for temperature-dependent cloud ice sticking efficiency
  tmin_iceautoconv = 188.15_ireals, & ! Temperature at which cloud ice autoconversion starts
  zceff_min        = 0.02_ireals      ! Minimum value for sticking efficiency

!------------------------------------------------------------------------------
!> Parameters and variables which are global in this module
!------------------------------------------------------------------------------

#ifdef __COSMO__
CHARACTER(132) :: message_text = ''
#endif


INTEGER (KIND=iintegers), PARAMETER ::  &
  iautocon       = 1,&
  isnow_n0temp   = 2

REAL    (KIND=ireals   ), PARAMETER ::  &
  zqmin = 1.0E-15_ireals,& ! threshold for computations
  zeps  = 1.0E-15_ireals,& ! small number
  
                                ! Parameters for autoconversion of cloud water and cloud ice 
  zccau  = 4.0E-4_ireals, & ! autoconversion coefficient (cloud water to rain)
  zciau  = 1.0E-3_ireals, & ! autoconversion coefficient (cloud ice   to snow)
  
  zkcau  = 9.44e+09_ireals, & ! kernel coeff for SB2001 autoconversion
  zkcac  = 5.25e+00_ireals, & ! kernel coeff for SB2001 accretion
  zcnue  = 2.00e+00_ireals, & ! gamma exponent for cloud distribution
  zxstar = 2.60e-10_ireals, & ! separating mass between cloud and rain
  zkphi1 = 6.00e+02_ireals, & ! constant in phi-function for autoconversion
  zkphi2 = 0.68e+00_ireals, & ! exponent in phi-function for autoconversion
  zkphi3 = 5.00e-05_ireals, & ! exponent in phi-function for accretion
  
  zhw   = 2.270603,       & ! Howell factor
  zecs  = 0.9_ireals,     & ! Collection efficiency for snow collecting cloud water
  
  zadi  = 0.217_ireals,   & ! Formfactor in the size-mass relation of ice particles
  zbdi  = 0.302_ireals,   & ! Exponent in the size-mass relation of ice particles
  zams  = 0.069_ireals,   & ! Formfactor in the mass-size relation of snow particles
  zbms  = 2.000_ireals,   & ! Exponent in the mass-size relation of snow particles
  
  zv1s  = 0.50_ireals,    & ! Exponent in the terminal velocity for snow
  
  zami  = 130.0_ireals,   & ! Formfactor in the mass-size relation of cloud ice
  zn0s0 = 8.0E5_ireals,   & ! 
  zn0s1 = 13.5_ireals * 5.65E5_ireals, & ! parameter in N0S(T)
  zn0s2 = -0.107_ireals , & ! parameter in N0S(T), Field et al
  zcac  = 1.72_ireals   , & ! (15/32)*(PI**0.5)*(ECR/RHOW)*V0R*AR**(1/8)
  zcicri= 1.72_ireals   , & ! (15/32)*(PI**0.5)*(EIR/RHOW)*V0R*AR**(1/8)
  zcrcri= 1.24E-3_ireals, & ! (PI/24)*EIR*V0R*Gamma(6.5)*AR**(-5/8)
  zcsmel= 1.48E-4_ireals, & ! 4*LHEAT*N0S*AS**(-2/3)/(RHO*lh_f)
  zbsmel= 20.32_ireals  , & !        0.26*sqrt(    RHO*v0s/eta)*Gamma(21/8)*AS**(-5/24)
!FR old    
! zasmel= 2.31E3_ireals , & ! DIFF*lh_v*RHO/LHEAT
  zasmel= 2.43E3_ireals , & ! DIFF*lh_v*RHO/LHEAT
 
  zcrfrz= 1.68_ireals    , & ! coefficient for raindrop freezing
  zcrfrz1= 9.95e-5_ireals, & !FR: 1. coefficient for immersion raindrop freezing: alpha_if
  zcrfrz2= 0.66_ireals   , & !FR: 2. coefficient for immersion raindrop freezing: a_if

  zrho0 = 1.225e+0_ireals, & ! reference air density
  zrhow = 1.000e+3_ireals, & ! density of liquid water
  
  zdv    = 2.22e-5_ireals, & ! molecular diffusion coefficient for water vapour
  zlheat = 2.40E-2_ireals, & ! thermal conductivity of dry air
  zeta   = 1.75e-5_ireals    ! kinematic viscosity of air 

REAL    (KIND=ireals   ), PARAMETER ::  &
  !! Additional parameters
  zthet = 248.15_ireals , & ! temperature for het. nuc. of cloud ice
  zthn  = 236.15_ireals , & ! temperature for hom. freezing of cloud water
  ztrfrz= 271.15_ireals , & ! threshold temperature for heterogeneous
!                           ! freezing of raindrops
  ztmix = 250.00_ireals , & ! threshold temperature for mixed-phase clouds (Forbes 2012)
!                           ! freezing of raindrops
  znimax_Thom = 250.E+3_ireals, & !FR: maximal number of ice crystals 
  zmi0   = 1.0E-12_ireals, & ! initial crystal mass for cloud ice nucleation
  zmimax = 1.0E-9_ireals , & ! maximum mass of cloud ice crystals   
  zmsmin = 3.0E-9_ireals , & ! initial mass of snow crystals        
  zvz0i  = 1.1_ireals    , & ! Terminal fall velocity of ice (Heymsfield+Donner 1990 multiplied by 1/3)
  zbvi   = 0.16_ireals   , & ! v = zvz0i*rhoqi^zbvi

  !! Constant exponents in the transfer rate equations
  x1o12  =  1.0_ireals/12.0_ireals  ,  x3o16  =  3.0_ireals/16.0_ireals, &
  x7o8   =  7.0_ireals/ 8.0_ireals  ,  x2o3   =  2.0_ireals/ 3.0_ireals, &
  x5o24  =  5.0_ireals/24.0_ireals  ,  x1o8   =  1.0_ireals/ 8.0_ireals, &
  x13o8  = 13.0_ireals/ 8.0_ireals  ,  x13o12 = 13.0_ireals/12.0_ireals, &
  x27o16 = 27.0_ireals/16.0_ireals  ,  x1o3   =  1.0_ireals/ 3.0_ireals, &
  x1o2   =  1.0_ireals/ 2.0_ireals  ,  x3o4   =  0.75_ireals           , &
  x7o4   =  7.0_ireals/ 4.0_ireals 

!------------------------------------------------------------------------------
! parameters relevant to support supercooled liquid water (SLW)
 REAL    (KIND=ireals   ), PARAMETER ::  &
  !FR: parameters for SLW improvement
   dist_cldtop_ref = 500.0_ireals, & ! Reference length for distance from cloud top (Forbes 2012)
   reduce_dep_ref  = 0.1_ireals      ! lower bound on snow/ice deposition reduction

!!==============================================================================

!> Global Variables
! -------------------

  REAL (KIND=ireals) ::          &
    zn0r,                        & ! N0_rain
    zar



!> Namelist Variables for hydci_pp and hydci_pp_gr
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
!>  Module procedure "hydci_pp_init" in "gscp" to initialize some
!!  coefficients which are used in "hydci_pp"
!------------------------------------------------------------------------------

SUBROUTINE hydci_pp_init(idbg)

!------------------------------------------------------------------------------
!> Description:
!!   Calculates some coefficients for hydci_pp. Usually called only once at
!!   model startup.
!------------------------------------------------------------------------------

  INTEGER, INTENT(IN), OPTIONAL::  &
    idbg              !! debug level

!------------------------------------------------------------------------------
!>  Initial setting of local and global variables
!------------------------------------------------------------------------------

!!  mu_rain = 0.5  ! is a namelist parameter
! FR new:
#ifdef __ICON__
  mu_rain = atm_phy_nwp_config(1)%mu_rain
!  mu_snow = atm_phy_nwp_config(1)%mu_snow
#endif


  zconst = zkcau / (20.0_ireals*zxstar*cloud_num*cloud_num) &
    * (zcnue+2.0_ireals)*(zcnue+4.0_ireals)/(zcnue+1.0_ireals)**2
  ccsrim = 0.25_ireals*pi*zecs*v0snow*gamma_fct(zv1s+3.0_ireals)
  ccsagg = 0.25_ireals*pi*v0snow*gamma_fct(zv1s+3.0_ireals)
  ccsdep = 0.26_ireals*gamma_fct((zv1s+5.0_ireals)/2.d0)*SQRT(1.0_ireals/zeta)
  ccsvxp = -(zv1s/(zbms+1.0_ireals)+1.0_ireals)
  ccsvel = zams*v0snow*gamma_fct(zbms+zv1s+1.0_ireals)&
    & *(zams*gamma_fct(zbms+1.0_ireals))**ccsvxp 
  ccsvxp = ccsvxp + 1.0_ireals
  ccslam = zams*gamma_fct(zbms+1.0_ireals)
  ccslxp = 1.0_ireals / (zbms+1.0_ireals)
  ccswxp = zv1s*ccslxp
  ccsaxp = -(zv1s+3.0_ireals)
  ccsdxp = -(zv1s+1.0_ireals)/2.0_ireals
  ccshi1 = lh_s*lh_s/(zlheat*r_v)
  ccdvtp = 2.22E-5 * t0**(-1.94) * 101325.0_ireals
  ccidep = 4.0_ireals * zami**(-x1o3)
  zn0r   = 8e6 * EXP(3.2*mu_rain) * (0.01)**(-mu_rain)  ! empirical relation adapted from Ulbrich (1983)
! to tune the zn0r variable
  zn0r   = zn0r * rain_n0_factor

  zar    = pi*zrhow/6.0 * zn0r * gamma_fct(mu_rain+4.0) ! pre-factor in lambda
  zcevxp = (mu_rain+2.)/(mu_rain+4.)
  zcev   = 2.0*pi*zdv/zhw*zn0r*zar**(-zcevxp) * gamma_fct(mu_rain+2.0)
  zbevxp = (2.*mu_rain+5.5_ireals)/(2.*mu_rain+8.)-zcevxp
  zbev   =  0.26 * SQRT(    zrho0*130.0_ireals/zeta)*zar**(-zbevxp) &
    &    * gamma_fct((2.0*mu_rain+5.5)/2.0) / gamma_fct(mu_rain+2.0)

  zvzxp  = 0.5/(mu_rain+4.0)
  zvz0r  = 130.0_ireals*gamma_fct(mu_rain+4.5)/gamma_fct(mu_rain+4.0)*zar**(-zvzxp)

#ifdef __ICON__  
  !CK> for cloud ice sedimentation based on KC05
  vtxexp = kc_beta + 2.0_ireals - kc_sigma 
  kc_c1  = 4.0_ireals / ( do_i**2 * SQRT(co_i) )
  kc_c2  = do_i **2 / 4.0_ireals
  !CK<
#endif    

  IF (PRESENT(idbg)) THEN
    IF (idbg > 10) THEN
      CALL message('gscp_hydci_pp','hydci_pp_init: Initialized coefficients for hydci_pp')
      WRITE (message_text,'(A,E10.3)') '      ccslam = ',ccslam ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsvel = ',ccsvel ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsrim = ',ccsrim ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsagg = ',ccsagg ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsdep = ',ccsdep ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccslxp = ',ccslxp ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccidep = ',ccidep ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      mu_r   = ',mu_rain; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zn0r   = ',zn0r   ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zbev   = ',zbev   ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zbevxp = ',zbevxp ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zcev   = ',zcev   ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zcevxp = ',zcevxp ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zvzxp  = ',zvzxp  ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      zvz0r  = ',zvz0r  ; CALL message('',message_text)
#ifdef __ICON__
      WRITE (message_text,'(A,E10.3)') '      vtxexp = ',vtxexp ; CALL message('',message_text) 
      WRITE (message_text,'(A,E10.3)') '       kc_c1 = ',kc_c1  ; CALL message('',message_text) 
      WRITE (message_text,'(A,E10.3)') '       kc_c2 = ',kc_c2  ; CALL message('',message_text) 
#endif
    ENDIF
  END IF

  CALL message('hydci_pp_init','microphysical values initialized')

END SUBROUTINE hydci_pp_init












!==============================================================================
!> Module procedure "hydci_pp" in "gscp" for computing effects of grid scale
!!  precipitation including cloud water, cloud ice, rain and snow in
!!  context with the Leapfrog and the Runge-Kutta time-integration
!------------------------------------------------------------------------------

SUBROUTINE hydci_pp (             &
  nvec,ke,                           & !> array dimensions
  ivstart,ivend, kstart,             & !! optional start/end indicies
  idbg,                              & !! optional debug level
  zdt, dz,                           & !! numerics parameters
  t,p,rho,qv,qc,qi,qr,qs,            & !! prognostic variables
#ifdef __ICON__
  !xxx: this should become a module variable, e.g. in a new module mo_data_gscp.f90
  qi0,qc0,                           & !! cloud ice/water threshold for autoconversion
#endif
  prr_gsp,prs_gsp,                   & !! surface precipitation rates
#ifdef NUDGING
  tinc_lh,                           & !  t-increment due to latent heat 
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
  ddt_diag_shed                      ) !!

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
!!   The
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

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(in)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=ireals), DIMENSION(nvec,ke), INTENT(IN) ::      &
#else
  REAL(KIND=ireals), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
#endif
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv                   !! if true, cv is used instead of cp

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=ireals), DIMENSION(nvec,ke), INTENT(INOUT) ::   &
#else
  REAL(KIND=ireals), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
#endif
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    qr              ,    & !! specific rain content                         (kg/kg)
    qs                     !! specific snow content                         (kg/kg)

#ifdef NUDGING
  REAL(KIND=ireals), INTENT(INOUT) :: &
       tinc_lh(:,:)   ,  & ! temperature increments due to heating             ( K/s )   
       tt_lheat(:,:)  ,  & !  t-increments due to latent heating (nudg) ( K/s )
       qrsflux(:,:)       ! total precipitation flux (nudg)
#endif

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=ireals), DIMENSION(nvec), INTENT(INOUT) ::   &
#else
  REAL(KIND=ireals), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
#endif
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    prs_gsp                !! precipitation rate of snow, grid-scale        (kg/(m2*s))

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(out)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=ireals), DIMENSION(nvec,ke), INTENT(OUT), OPTIONAL ::   &
#else
  REAL(KIND=ireals), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &     ! dim (ie,ke)
#endif
    ddt_tend_t      , & !> tendency T                                       ( 1/s )
    ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
    ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
    ddt_tend_qi     , & !! tendency qi                                      ( 1/s )
    ddt_tend_qr     , & !! tendency qr                                      ( 1/s )
    ddt_tend_qs         !! tendency qs                                      ( 1/s )

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(out)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=ireals), DIMENSION(nvec,ke), INTENT(OUT), OPTIONAL ::   &
#else
  REAL(KIND=ireals), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &   ! dim (ie,ke)
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
    zqvsw,             & ! sat. specitic humidity at ice and water saturation
    zztau, zxfac, zx1, zx2, ztt, &   ! some help variables
    ztau, zphi, zhi, zdvtp, ztc, zlog_10

  REAL    (KIND=ireals   ) ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice
    zqrt   ,& ! layer tendency of rain
    zqst      ! layer tendency of snow

  REAL (KIND=ireals)         ::       &
    mma(10), mmb(10)

  REAL    (KIND=ireals   ) ::  &
    zlnqrk,zlnqsk,     & !
    zlnlogmi,qcgk_1,               & !
    qcg,tg,qvg,qrg, qsg,qig,rhog,ppg, alf,bet,m2s,m3s,hlp,maxevap,temp_c,stickeff

  LOGICAL :: &
    llqr,llqs,llqc,llqi  !   switch for existence of qr, qs, qc, qi

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
    zeln7o8qrk  (nvec),     & !
    zeln7o4qrk  (nvec),     & ! FR new     
    zeln27o16qrk(nvec),     & !
    zeln13o8qrk (nvec),     & !
!    zeln3o16qrk (nvec),     & !
    zeln13o12qsk(nvec),     & !
    zeln5o24qsk (nvec),     & !
    zeln2o3qsk  (nvec)        !

  REAL    (KIND=ireals   ) ::  &
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
  INTEGER (KIND=iintegers) ::  &
    ic1, ic2, ic3, ic4, ic5, ic6, i1d

  !> Integer arrays for a better vectorization
  INTEGER (KIND=iintegers) ::  &
    ivdx1(nvec), & !!
    ivdx2(nvec), & !!
    ivdx3(nvec), & !!
    ivdx4(nvec), & !!
    ivdx5(nvec), & !!
    ivdx6(nvec)    !!


!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine hydci_pp
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
  fxna_cooper(ztx) = 5.0E+0_ireals * EXP(0.304_ireals * (t0 - ztx))   ! FR: Cooper (1986) used by Greg Thompson(2008)

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
  IF( lsuper_coolw) THEN
    znimax = znimax_Thom         !znimax_Thom = 250.E+3_ireals,
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSEIF(lorig_icon) THEN
    znimax = 150.E+3_ireals     ! from previous ICON code 
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSE
    znimax = fxna(zthn) ! Maximum number of cloud ice crystals
    znimix = fxna(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  END IF

  zpvsw0 = fpvsw(t0)  ! sat. vap. pressure for t = t0
  zlog_10 = LOG(10._ireals) ! logarithm of 10

! Delete precipitation fluxes from previous timestep
!CDIR BEGIN COLLAPSE
    prr_gsp (:) = 0.0_ireals
    prs_gsp (:) = 0.0_ireals
    zpkr    (:) = 0.0_ireals
    zpks    (:) = 0.0_ireals
    zpki    (:) = 0.0_ireals
    zprvr   (:) = 0.0_ireals
    zprvs   (:) = 0.0_ireals
    zprvi   (:) = 0.0_ireals
    zvzr    (:) = 0.0_ireals
    zvzs    (:) = 0.0_ireals
    zvzi    (:) = 0.0_ireals
    dist_cldtop(:) = 0.0_ireals    
!CDIR END


! Optional arguments

  IF (PRESENT(ddt_tend_t)) THEN
    ! save input arrays for final tendency calculation
    t_in  = t
    qv_in = qv
    qc_in = qc
    qi_in = qi
    qr_in = qr
    qs_in = qs
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
  zdtr  = 1.0_ireals / zdt

#ifdef NUDGING
  ! add part of latent heating calculated in subroutine hydci_pp to model latent
  ! heating field: subtract temperature from model latent heating field
  IF (llhn .OR. llhnverif) THEN
    IF (lhn_qrs) THEN
!CDIR COLLAPSE
      qrsflux(:,:) = 0.0_ireals
    ENDIF
    tt_lheat(:,:) = tt_lheat(:,:) - t(:,:)
!    CALL get_gs_lheating ('add',1,ke) !should not be called from block physics
  ENDIF
#endif


! output for various debug levels
  IF (izdebug > 15) CALL message('','SRC_GSCP: Start of hydci_pp')
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
    zcrim (:) = 0.0_ireals
    zcagg (:) = 0.0_ireals
    zbsdep(:) = 0.0_ireals
    zvz0s (:) = 0.0_ireals
    zn0s  (:) = zn0s0
    reduce_dep(:) = 1.0_ireals  !FR: Reduction coeff. for dep. growth of rain and ice  
!CDIR END

    ic1 = 0
    ic2 = 0
    ic3 = 0

    DO iv = iv_start, iv_end

        qrg = qr(iv,k)
        qsg = qs(iv,k)
        qvg = qv(iv,k)
        qcg = qc(iv,k)
        qig = qi(iv,k)
        tg  = t(iv,k)
        ppg = p(iv,k)
        rhog = rho(iv,k)

        !..for density correction of fall speeds
        z1orhog(iv) = 1.0_ireals/rhog
        zrho1o2(iv) = EXP(LOG(zrho0*z1orhog(iv))*x1o2)

        zqrk(iv) = qrg * rhog
        zqsk(iv) = qsg * rhog
        zqik(iv) = qig * rhog

        llqr = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin
        llqi = zqik(iv) > zqmin

        zdtdh(iv) = 0.5_ireals * zdt / dz(iv,k)

        zzar(iv)   = zqrk(iv)/zdtdh(iv) + zprvr(iv) + zpkr(iv)
        zzas(iv)   = zqsk(iv)/zdtdh(iv) + zprvs(iv) + zpks(iv)
        zzai(iv)   = zqik(iv)/zdtdh(iv) + zprvi(iv) + zpki(iv)

        IF (llqs) THEN
          ic1 = ic1 + 1
          ivdx1(ic1) = iv
        ENDIF
        IF (llqr) THEN
          ic2 = ic2 + 1
          ivdx2(ic2) = iv
       ENDIF
       IF (llqi) THEN
         ic3 = ic3 + 1
         ivdx3(ic3) = iv
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
        alf = EXP(hlp*zlog_10) ! 10.0_ireals**hlp
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

!      IF (zvzs(iv) == 0.0_ireals) THEN
        zvzs(iv) = zvz0s(iv) * EXP (ccswxp * LOG (0.5_ireals*zqsk(iv))) * zrho1o2(iv)
!      ENDIF
    ENDDO loop_over_qs_prepare

!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
    loop_over_qr_sedi: DO i1d = 1, ic2
      iv = ivdx2(i1d)

!      IF (zvzr(iv) == 0.0_ireals) THEN
!replaced: zvzr(iv) = zvz0r * EXP (x1o8  * LOG (0.5_ireals*zqrk(iv))) * zrho1o2(iv)
        zvzr(iv) = zvz0r * EXP (zvzxp  * LOG (0.5_ireals*zqrk(iv))) * zrho1o2(iv)

!      ENDIF
    ENDDO loop_over_qr_sedi

!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
    loop_over_qi_sedi: DO i1d = 1, ic3
      iv = ivdx3(i1d)
      
!      IF (zvzi(iv) == 0.0_ireals .AND. lnew_ice_sedi) THEN
        !! density correction not needed zrho1o2(iv)
        zvzi(iv) = zvz0i * EXP (zbvi  * LOG (0.5_ireals*zqik(iv)))
!      ENDIF
    ENDDO loop_over_qi_sedi


  !----------------------------------------------------------------------------
  ! Section 3:
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
    zeln7o8qrk   (:) = 0.0_ireals
    zeln7o4qrk   (:) = 0.0_ireals !FR
    zeln27o16qrk (:) = 0.0_ireals
    zeln13o8qrk  (:) = 0.0_ireals
!!    zeln3o16qrk  (:) = 0.0_ireals
    zeln13o12qsk (:) = 0.0_ireals
    zeln5o24qsk  (:) = 0.0_ireals
    zeln2o3qsk   (:) = 0.0_ireals

!FR old  
!   zcsdep       (:) = 3.2E-2_ireals
    zcsdep       (:) = 3.367E-2_ireals    
    zcidep       (:) = 1.3E-5_ireals
    zcslam       (:) = 1e10_ireals

    scau         (:) = 0.0_ireals
    scac         (:) = 0.0_ireals
    snuc         (:) = 0.0_ireals
    scfrz        (:) = 0.0_ireals
    simelt       (:) = 0.0_ireals
    sidep        (:) = 0.0_ireals
    ssdep        (:) = 0.0_ireals
    sdau         (:) = 0.0_ireals
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
    DO iv = iv_start, iv_end

        qrg  = qr(iv,k)
        qsg  = qs(iv,k)
        qvg  = qv(iv,k)
        qcg  = qc(iv,k)
        qig  = qi(iv,k)
        tg   =  t(iv,k)
        ppg  =  p(iv,k)
        rhog = rho(iv,k)

        zqvsi(iv) = sat_pres_ice(tg)/(rhog * r_v * tg)

        llqr = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin
        llqi = zqik(iv) > zqmin

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

        IF (llqi) THEN
          zpki(iv) = zqik (iv) * zvz0i * EXP (zbvi * LOG (zqik(iv))) * zrho1o2(iv)
        ELSE
          zpki(iv) = 0.0_ireals
        ENDIF

        zpkr(iv)   = MIN( zpkr(iv) , zzar(iv) )
        zpks(iv)   = MIN( zpks(iv) , zzas(iv) )
        zpki(iv)   = MIN( zpki(iv) , zzai(iv) )

        zzar(iv)   = zdtdh(iv) * (zzar(iv)-zpkr(iv))
        zzas(iv)   = zdtdh(iv) * (zzas(iv)-zpks(iv))
        zzai(iv)   = zdtdh(iv) * (zzai(iv)-zpki(iv))

        zimr(iv)   = 1.0_ireals / (1.0_ireals + zvzr(iv) * zdtdh(iv))
        zims(iv)   = 1.0_ireals / (1.0_ireals + zvzs(iv) * zdtdh(iv))
        zimi(iv)   = 1.0_ireals / (1.0_ireals + zvzi(iv) * zdtdh(iv))

        zqrk(iv)   = zzar(iv)*zimr(iv)
        zqsk(iv)   = zzas(iv)*zims(iv)
        zqik(iv)   = zzai(iv)*zimi(iv)

        llqr = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin
        llqi = zqik(iv) > zqmin
        llqc =       qcg > zqmin
      ! llqi =       qig > zqmin 

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
        IF ( tg < zthet .AND. qvg >  8.E-6_ireals &
                        .AND. qig <= 0.0_ireals ) THEN
          ic4 = ic4 + 1
          ivdx4(ic4) = iv
        ENDIF
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
        zeln7o4qrk(iv)   = EXP (x7o4   * zlnqrk) !FR new
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
!DIR$ IVDEP
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
      qcgk_1 = qc(iv,k-1) 
      qig  =   qi(iv,k)
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
          ztau   = MAX(ztau,1.E-30_ireals)
          hlp    = EXP(zkphi2*LOG(ztau))
          zphi   = zkphi1 * hlp * (1.0_ireals - hlp)**3
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
      IF( (tg <= 267.15_ireals) .AND. (qig <= 0.0_ireals)) THEN
        IF  (lsuper_coolw .OR. lorig_icon) THEN
          znin  = MIN( fxna_cooper(tg), znimax )
          snuc(iv) = zmi0 * z1orhog(iv) * znin * zdtr
        ELSE 
          znin  = MIN( fxna(tg), znimax )
          snuc(iv) = zmi0 * z1orhog(iv) * znin * zdtr
        END IF
      ENDIF
      ! Calculation of in-cloud rainwater freezing
      IF (tg < ztrfrz)  THEN
        IF (lsuper_coolw) THEN
          srfrz(iv) = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_ireals ) * zeln7o4qrk(iv)
        ELSE 
          srfrz(iv) = zcrfrz*SQRT( (ztrfrz-tg)**3 )* zeln27o16qrk(iv)
        ENDIF
      ENDIF
      
!FR>>> Calculation of reduction of depositional growth at cloud top (Forbes 2012)
      IF( k>1 .AND. k<ke .AND. lsuper_coolw ) THEN
        znin  = MIN( fxna_cooper(tg), znimax )
        fnuc = MIN(znin/znimix, 1.0_ireals)
        
        !! distance from cloud top
        IF( qcgk_1 .LT. zqmin ) THEN      ! upper cloud layer
          dist_cldtop(iv) = 0.0_ireals    ! reset distance to upper cloud layer
        ELSE
          dist_cldtop(iv) = dist_cldtop(iv) + dz(iv,k)
        END IF
         
! with asymptotic behaviour dz -> 0 (xxx)
!        reduce_dep(iv) = MIN(fnuc + (1.0_ireals-fnuc)*(reduce_dep_ref + &
!                             dist_cldtop(iv)/dist_cldtop_ref + &
!                             (1.0_ireals-reduce_dep_ref)*(zdh/dist_cldtop_ref)**4), 1.0_ireals)
        
! without asymptotic behaviour dz -> 0
        reduce_dep(iv) = MIN(fnuc + (1.0_ireals-fnuc)*(reduce_dep_ref + &
                          dist_cldtop(iv)/dist_cldtop_ref), 1.0_ireals)
       
      END IF ! Reduction of dep. growth of snow/ice 
!FR<<<
      
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
          stickeff = MIN(EXP(0.09_ireals*(tg-t0)),1.0_ireals)
          stickeff = MAX(stickeff, zceff_min, zceff_fac*(tg-tmin_iceautoconv))
        ELSE !original sticking efficiency of cloud ice
          stickeff = MIN(EXP(0.09_ireals*(tg-t0)),1.0_ireals)
          stickeff = MAX(stickeff,0.2_ireals)
        END IF
        zmi      = MIN( rhog*qig/znin, zmimax )
        zmi      = MAX( zmi0, zmi )
        zsvmax   = (qvg - zqvsi(iv)) * zdtr
        zsagg    = zcagg(iv) * EXP(ccsaxp*LOG(zcslam(iv))) * qig
        zsagg    = MAX( zsagg, 0.0_ireals ) * stickeff
        znid     = rhog * qig/zmi
        IF (llqi) THEN
          zlnlogmi = LOG (zmi)
          zsidep   = zcidep(iv) * znid * EXP(0.33_ireals * zlnlogmi) * (qvg - zqvsi(iv))
        ELSE
          zsidep = 0.0_ireals
        ENDIF
        zsvidep   = 0.0_ireals
        zsvisub   = 0.0_ireals
        zsimax    = qig*zdtr 
        IF( zsidep > 0.0_ireals ) THEN
          IF (lsuper_coolw ) THEN
            zsidep = zsidep * reduce_dep(iv)  !FR new: SLW reduction
          END IF
          zsvidep = MIN( zsidep, zsvmax )
        ELSEIF (zsidep < 0.0_ireals ) THEN
          zsvisub = - MAX(-zsimax, zsvmax )
        ENDIF
        zsiau = zciau * MAX( qig - qi0, 0.0_ireals ) * stickeff
        IF (llqi) THEN
          zlnlogmi = LOG(zmsmin/zmi)
          zztau    = 1.5_ireals*( EXP(0.66_ireals*zlnlogmi) - 1.0_ireals)
          zsdau    = zsvidep/MAX(zztau,zeps)
        ELSE
          zsdau    =  0.0_ireals
        ENDIF
        zsicri    = zcicri * qig * zeln7o8qrk(iv)
        zsrcri    = zcrcri * (qig/zmi) * zeln13o8qrk(iv)
        zxfac     = 1.0_ireals + zbsdep(iv) * EXP(ccsdxp*LOG(zcslam(iv)))
        zssdep    = zcsdep(iv) * zxfac * ( qvg - zqvsi(iv) ) / (zcslam(iv)+zeps)**2

        ! Check for maximal depletion of vapor by sdep
        IF (zssdep > 0.0_ireals) THEN
          IF (lsuper_coolw) THEN
            zssdep = zssdep*reduce_dep(iv)!FR new: SLW reduction
          END IF
          zssdep = MIN(zssdep, zsvmax-zsvidep)
        END IF
        
        ! Check for maximal depletion of snow by sdep
        IF (zssdep < 0.0_ireals) zssdep = MAX(zssdep, -qsg*zdtr)

        zsisum = zsiau + zsdau + zsagg + zsicri + zsvisub
        zcorr  = 0.0_ireals
        IF( zsimax > 0.0_ireals ) zcorr  = zsimax / MAX( zsimax, zsisum )
        sidep(iv)  = zsvidep - zcorr*zsvisub
        sdau (iv)  = zcorr*zsdau
        siau (iv)  = zcorr*zsiau
        sagg (iv)  = zcorr*zsagg
        ssdep(iv)  = zssdep
        srcri(iv)  = zsrcri
        sicri(iv)  = zcorr*zsicri


      !------------------------------------------------------------------------
      ! Section 7: Search for warm grid points with cloud ice and/or snow and
      !            calculation of the melting rates of qi and ps
      !------------------------------------------------------------------------

      ELSE ! tg > 0
        simelt(iv) = qig*zdtr
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
      ! Limit evaporation rate in order to avoid overshoots towards supersaturation
      ! the pre-factor approximates (esat(T_wb)-e)/(esat(T)-e) at temperatures between 0 degC and 30 degC
      temp_c = tg - t0
      maxevap     = (0.61_ireals-0.0163_ireals*temp_c+1.111e-4_ireals*temp_c**2)*(zqvsw-qvg)/zdt
      sev(iv)    = MIN(zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk), maxevap)
!      sev(iv)    = zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk)

!      zqvsw    = fqvs( fpvsw(tg), ppg )
!      zx1      = 1.0_ireals + zbev* zeln3o16qrk(iv)
!      zsev     = zcev*zx1*(zqvsw - qvg)*SQRT(zqrk(iv))
!      sev(iv) = MAX( zsev, 0.0_ireals )

      ! Calculation of below-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        IF ( lsuper_coolw ) THEN
          !FR new: reduced rain freezing rate
          srfrz(iv) = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_ireals ) * zeln7o4qrk(iv)
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

!        qvg = qv(iv,k)
!        qcg = qc(iv,k)
        qrg = qr(iv,k)
        qsg = qs(iv,k)
        qig = qi(iv,k)
!        tg  = t (iv,k)
        rhog = rho(iv,k)

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
        zqvt = sev(iv)   - sidep(iv) - ssdep(iv)  - snuc(iv)
        zqct = simelt(iv)- scau(iv)  - scfrz(iv)  - scac(iv)   - sshed(iv) - srim(iv) 
        zqit = snuc(iv)  + scfrz(iv) - simelt(iv) - sicri(iv)  + sidep(iv) - sdau(iv)            &
                                                               - sagg(iv) - siau(iv)
        zqrt = scau(iv)  + sshed(iv) + scac(iv)   + ssmelt(iv) - sev(iv) - srcri(iv) - srfrz(iv) 
        zqst = siau(iv)  + sdau(iv)  + sagg(iv)   - ssmelt(iv) + sicri(iv) + srcri(iv)           &
                                                               + srim(iv) + ssdep(iv) + srfrz(iv)
        ztt = z_heat_cap_r*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst) )

        ! Update variables and add qi to qrs for water loading 
        IF (lsedi_ice .OR. lorig_icon) THEN
          qig = MAX ( 0.0_ireals, (zzai(iv)*z1orhog(iv) + zqit*zdt)*zimi(iv))
        ELSE
          qig = MAX ( 0.0_ireals, qig + zqit*zdt)
        END IF
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

          IF (zprvr(iv) <= zqmin) zprvr(iv)=0.0_ireals
          IF (zprvs(iv) <= zqmin) zprvs(iv)=0.0_ireals
          IF (zprvi(iv) <= zqmin) zprvi(iv)=0.0_ireals

#ifdef NUDGING
          ! for the latent heat nudging
          IF ((llhn .OR. llhnverif) .AND. lhn_qrs) THEN
            IF (lsedi_ice .OR. lorig_icon) THEN
              qrsflux(iv,k) = zprvr(iv)+zprvs(iv)+zprvi(iv)
              qrsflux(iv,k) = 0.5*(qrsflux(iv,k)+zpkr(iv)+zpks(iv)+zpki(iv))
            ELSE 
              qrsflux(iv,k) = zprvr(iv)+zprvs(iv)                         
              qrsflux(iv,k) = 0.5*(qrsflux(iv,k)+zpkr(iv)+zpks(iv))       
            ENDIF
          ENDIF
#endif

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
          IF (qig+qi(iv,k+1) <= zqmin ) THEN
            zvzi(iv)= 0.0_ireals
          ELSE
            !! density correction not needed
            zvzi(iv)= zvz0i * EXP(zbvi*LOG((qig+qi(iv,k+1))*0.5_ireals*rhog))
          ENDIF

        ELSE
          ! Precipitation fluxes at the ground
          prr_gsp(iv) = 0.5_ireals * (qrg*rhog*zvzr(iv) + zpkr(iv))
          IF (lsedi_ice .OR. lorig_icon) THEN
            prs_gsp(iv) = 0.5_ireals * (rhog*(qsg*zvzs(iv)+qig*zvzi(iv)) + zpks(iv)+zpki(iv))
          ELSE
            prs_gsp(iv) = 0.5_ireals * (qsg*rhog*zvzs(iv) + zpks(iv))
          END IF
          

#ifdef NUDGING
          ! for the latent heat nudging
          IF ((llhn .OR. llhnverif) .AND. lhn_qrs)        &
            qrsflux(iv,k) = prr_gsp(iv)+prs_gsp(iv)
#endif

        ENDIF

        ! Update of prognostic variables or tendencies
        qr (iv,k) = qrg
        qs (iv,k) = MAX ( 0.0_ireals, qsg )
        qi (iv,k) = qig
!        qrs(iv,k   ) = qrg+qsg+qig       !qrs is now computed outside
        t  (iv,k) = t (iv,k) + ztt*zdt 
        qv (iv,k) = MAX ( 0.0_ireals, qv(iv,k) + zqvt*zdt )
        qc (iv,k) = MAX ( 0.0_ireals, qc(iv,k) + zqct*zdt )

      ENDDO loop_over_all_iv

  IF (izdebug > 15) THEN
    ! Check for negative values
     DO iv = iv_start, iv_end
        IF (qr(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qr'
          CALL message('',message_text)
        ENDIF
        IF (qc(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qc'
          CALL message('',message_text)
        ENDIF
        IF (qi(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qi'
          CALL message('',message_text)
        ENDIF
        IF (qs(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qs'
          CALL message('',message_text)
        ENDIF
        IF (qv(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qv'
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

#ifdef NUDGING
! add part of latent heating calculated in subroutine hydci to model latent
! heating field: add temperature to model latent heating field
IF (llhn .OR. llhnverif) &
! CALL get_gs_lheating ('inc',1,ke)  !this should be called from within the block
     tt_lheat(:,:) = tt_lheat(:,:) + t(:,:)
#endif

#ifdef __ICON__

CALL satad_v_3d (                             &
               & maxiter  = 10_iintegers ,& !> IN
               & tol      = 1.e-3_ireals ,& !> IN
               & te       = t            ,&
               & qve      = qv           ,&
               & qce      = qc           ,&
               & rhotot   = rho          ,&
               & idim     = nvec         ,&
               & kdim     = ke           ,&
               & ilo      = iv_start     ,&
               & iup      = iv_end       ,&
               & klo      = k_start      ,&
               & kup      = ke            &
               )

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

      END DO
    END DO

  END IF

  IF (izdebug > 15) THEN
    CALL message('gscp_hydci_pp', 'UPDATED VARIABLES')
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp T= ',&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp qi= ',&
    MAXVAL( qi(:,:)), MINVAL(qi(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp qs= ',&
    MAXVAL( qs(:,:)), MINVAL(qs(:,:) )
    CALL message('', TRIM(message_text))
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine hydci_pp
!------------------------------------------------------------------------------

END SUBROUTINE hydci_pp



!------------------------------------------------------------------------------
! Begin Graupel Scheme
!------------------------------------------------------------------------------

SUBROUTINE hydci_pp_gr (             &
  nvec,ke,                           & !> array dimensions
  ivstart,ivend, kstart,             & !! optional start/end indicies
  idbg,                              & !! optional debug level
  zdt, dz,                           & !! numerics parameters
  t,p,rho,qv,qc,qi,qr,qs,qg,         & !! prognostic variables
#ifdef __ICON__
  !xxx: this should become a module variable, e.g. in a new module mo_data_gscp.f90
  qi0,qc0,                           & !! cloud ice/water threshold for autoconversion
#endif
  prr_gsp,prs_gsp,prg_gsp,           & !! surface precipitation rates
#ifdef NUDGING
  tinc_lh,                           & !  t-increment due to latent heat 
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
!   The subroutine relies mostly on conversion terms used in hydci_pp.
!   In contrast to hydci_pp, qc0 = 0.0002 (instead of 0.0) is used!
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

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(in)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=ireals), DIMENSION(nvec,ke), INTENT(IN) ::      &
#else
  REAL(KIND=ireals), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
#endif
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv                   !! if true, cv is used instead of cp

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=ireals), DIMENSION(nvec,ke), INTENT(INOUT) ::   &
#else
  REAL(KIND=ireals), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
#endif
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    qr              ,    & !! specific rain content                         (kg/kg)
    qs              ,    & !! specific snow content                         (kg/kg)
    qg                     !! specific graupel content                      (kg/kg)
    
#ifdef NUDGING
  REAL(KIND=ireals), INTENT(INOUT) :: &
       tinc_lh(:,:)   ,  & ! temperature increments due to heating             ( K/s )   
       tt_lheat(:,:)  ,  & !  t-increments due to latent heating (nudg) ( K/s )
       qrsflux(:,:)       ! total precipitation flux (nudg)
#endif

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=ireals), DIMENSION(nvec), INTENT(INOUT) ::   &
#else
  REAL(KIND=ireals), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
#endif
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    prs_gsp,             & !! precipitation rate of snow, grid-scale        (kg/(m2*s))
    prg_gsp                !! precipitation rate of graupel, grid-scale     (kg/(m2*s))

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(out)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=ireals), DIMENSION(nvec,ke), INTENT(OUT), OPTIONAL ::   &
#else
  REAL(KIND=ireals), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &     ! dim (ie,ke)
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
  REAL(KIND=ireals), DIMENSION(nvec,ke), INTENT(OUT), OPTIONAL ::   &
#else
  REAL(KIND=ireals), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &   ! dim (ie,ke)
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

  REAL    (KIND=ireals   ), PARAMETER ::  &
    zams  = 0.038_ireals, & ! Formfactor in the mass-size relation of snow particles
    zcsg=0.5_ireals,        & !coefficient for snow-graupel conversion by riming
    zcrim_g=4.43_ireals,    & !
    zrimexp_g=0.94878_ireals, &
    zcagg_g = 2.46_ireals , & !
    zasmel= 2.95E3_ireals , & ! DIFF*lh_v*RHO/LHEAT
    zexpsedg=0.217_ireals,  & ! exponent for graupel sedimentation
    zvz0g = 12.24_ireals  , & ! coefficient of sedimentation velocity for graupel
    ztcrit=3339.5_ireals      ! factor in calculation of critical temperature

 
  REAL    (KIND=ireals   ) ::  &
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

  REAL    (KIND=ireals   ) ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice
    zqrt   ,& ! layer tendency of rain
    zqst   ,& ! layer tendency of snow
    zqgt      ! layer tendency of graupel

  REAL (KIND=ireals)         ::       &
    mma(10), mmb(10)

  REAL    (KIND=ireals   ) ::  &
    zlnqrk,zlnqsk,     & !
    zlnlogmi, zlnqgk, zln1o2,               & !
    qcg,tg,qvg,qrg,qsg,qgg,qig,rhog,ppg,alf,bet,m2s,m3s,hlp, &
    qcgk_1,maxevap,temp_c

  LOGICAL :: &
    llqs,llqc,llqi,llqg  !   switch for existence of qr, qs, qc, qi

  LOGICAL :: &
    llqr(nvec)


  REAL(KIND=ireals), DIMENSION(nvec,ke) ::   &
    t_in               ,    & !> temperature                                   (  K  )
    qv_in              ,    & !! specific water vapor content                  (kg/kg)
    qc_in              ,    & !! specific cloud water content                  (kg/kg)
    qi_in              ,    & !! specific cloud ice   content                  (kg/kg)
    qr_in              ,    & !! specific rain content                         (kg/kg)
    qs_in              ,    & !! specific snow content                         (kg/kg)
    qg_in                     !! specific graupel content                      (kg/kg)



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
    zcsdep      (nvec),     & !
    zcidep      (nvec)
    
 REAL    (KIND=ireals   ) ::  &    
    zsrmax      (nvec),     & !
    zssmax      (nvec),     & !
    zsgmax      (nvec),     & !
    zvz0s       (nvec),     & !
    zcrim       (nvec),     & !
    zcagg       (nvec),     & !
    zbsdep      (nvec),     & !
    zcslam      (nvec),     & !
    zn0s        (nvec),     & !
    zimr        (nvec),     & !
    zims        (nvec),     & !
    zimg        (nvec),     & !
    zimi        (nvec),     & !
    zzar        (nvec),     & !
    zzas        (nvec),     & !
    zzag        (nvec),     & !
    zzai        (nvec),     & !
    zqrk        (nvec),     & !
    zqsk        (nvec),     & !
    zqgk        (nvec),     & !
    zqik        (nvec),     & !
    zdtdh       (nvec),     & !
    z1orhog     (nvec),     & ! 1/rhog
    zrho1o2     (nvec),     & ! (rho0/rhog)**1/2
    zeln7o8qrk  (nvec),     & !
    zeln7o4qrk  (nvec),     & ! FR new  
    zeln27o16qrk(nvec),     & !
    zeln13o8qrk (nvec),     & !
    zeln3o4qsk  (nvec),     & ! 
    zeln6qgk    (nvec),     & !
    zeln8qsk    (nvec),     & !
    zelnrimexp_g(nvec)


  REAL    (KIND=ireals   ) ::  &
    scau   (nvec), & ! transfer rate due to autoconversion of cloud water
    scac   (nvec), & ! transfer rate due to accretion of cloud water
    snuc   (nvec), & ! transfer rate due nucleation of cloud ice
    scfrz  (nvec), & ! transfer rate due homogeneous freezing of cloud water
    simelt (nvec), & ! transfer rate due melting of cloud ice
    sidep  (nvec), & ! transfer rate due depositional growth of cloud ice
    ssdep  (nvec), & ! transfer rate due depositional growth of snow
    sgdep  (nvec), & ! transfer rate due depositional growth of snow
    sdau   (nvec), & ! transfer rate due depositional cloud ice autoconversion
    srim   (nvec), & ! transfer rate due riming of snow
    srim2  (nvec), & ! transfer rate due riming of snow
    sconsg (nvec), & ! transfer rate due to conversion from snow to graupel by riming  
    sshed  (nvec), & ! transfer rate due shedding
    sicri  (nvec), & ! transfer rate due cloud ice collection by rain (sink qi)
    srcri  (nvec), & ! transfer rate due cloud ice collection by rain (sink qr)
    sagg   (nvec), & ! transfer rate due aggregation of snow and cloud ice
    sagg2  (nvec), & ! transfer rate due aggregation of snow and cloud ice
    siau   (nvec), & ! transfer rate due autoconversion of cloud ice
    ssmelt (nvec), & ! transfer rate due melting of snow
    sgmelt (nvec), & ! transfer rate due melting of snow
    sev    (nvec), & ! transfer rate due evaporation of rain
    sconr  (nvec), & ! transfer rate due to condensation on melting snow/graupel
    srfrz  (nvec), & ! transfer rate due to rainwater freezing
    reduce_dep(nvec),&!FR: coefficient: reduce deposition at cloud top (Forbes 2012)
    dist_cldtop(nvec) !FR: distance from cloud top layer 

  ! Dimensions and loop counter for storing the indices
  INTEGER (KIND=iintegers) ::  &
    ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, i1d

  !> Integer arrays for a better vectorization
  INTEGER (KIND=iintegers) ::  &
    ivdx1(nvec), & !!
    ivdx2(nvec), & !!
    ivdx3(nvec), & !!
    ivdx4(nvec), & !!
    ivdx5(nvec), & !!
    ivdx6(nvec), & !!
    ivdx7(nvec), & !!
    ivdx8(nvec)    !!

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine hydci_pp_gr
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
  fxna_cooper(ztx) = 5.0E+0_ireals * EXP(0.304_ireals * (t0 - ztx))   ! FR: Cooper (1986) used by Greg Thompson(2008)

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
   IF( lsuper_coolw) THEN
    znimax = znimax_Thom
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSE
    znimax = fxna(zthn) ! Maximum number of cloud ice crystals
    znimix = fxna(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  END IF

  zpvsw0 = fpvsw(t0)  ! sat. vap. pressure for t = t0
  zlog_10 = LOG(10._ireals) ! logarithm of 10
  
! Delete precipitation fluxes from previous timestep
!CDIR BEGIN COLLAPSE
    prr_gsp (:) = 0.0_ireals
    prs_gsp (:) = 0.0_ireals
    prg_gsp (:) = 0.0_ireals
    zpkr    (:) = 0.0_ireals
    zpks    (:) = 0.0_ireals
    zpkg    (:) = 0.0_ireals
    zpki    (:) = 0.0_ireals
    zprvr   (:) = 0.0_ireals
    zprvs   (:) = 0.0_ireals
    zprvg   (:) = 0.0_ireals
    zprvi   (:) = 0.0_ireals
    zvzr    (:) = 0.0_ireals
    zvzs    (:) = 0.0_ireals
    zvzg    (:) = 0.0_ireals
    zvzi    (:) = 0.0_ireals
    dist_cldtop(:) = 0.0_ireals   
!CDIR END


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
  zdtr  = 1.0_ireals / zdt

#ifdef NUDGING
  ! add part of latent heating calculated in subroutine hydci_pp_gr to model latent
  ! heating field: subtract temperature from model latent heating field
  IF (llhn .OR. llhnverif) THEN
    IF (lhn_qrs) THEN
!CDIR COLLAPSE
      qrsflux(:,:) = 0.0_ireals
    ENDIF
!    CALL get_gs_lheating ('add',1,ke) !this should be called from within the block
    tt_lheat(:,:) = tt_lheat(:,:) - t(:,:)
  ENDIF
#endif


! output for various debug levels
  IF (izdebug > 15) CALL message('','SRC_GSCP: Start of hydci_pp_gr')
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
    reduce_dep(:) = 1.0_ireals  !FR: Reduction coeff. for dep. growth of rain and ice  
!CDIR END

  !----------------------------------------------------------------------------
  ! 2.1: Preparations for computations and to check the different conditions
  !----------------------------------------------------------------------------
    
  ! Nullify counters
  ic1 = 0
  ic2 = 0
  ic3 = 0
  ic4 = 0

    DO iv = iv_start, iv_end

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
        z1orhog(iv) = 1.0_ireals/rhog
        zrho1o2(iv) = EXP(LOG(zrho0*z1orhog(iv))*x1o2)

        zqrk(iv) = qrg * rhog
        zqsk(iv) = qsg * rhog
        zqgk(iv) = qgg * rhog
        zqik(iv) = qig * rhog
        
        llqr(iv) = zqrk(iv) > zqmin
        llqs     = zqsk(iv) > zqmin
        llqg     = zqgk(iv) > zqmin
        llqi     = zqik(iv) > zqmin
        
        zdtdh(iv) = 0.5_ireals * zdt / dz(iv,k)

        zzar(iv)   = zqrk(iv)/zdtdh(iv) + zprvr(iv) + zpkr(iv)
        zzas(iv)   = zqsk(iv)/zdtdh(iv) + zprvs(iv) + zpks(iv)
        zzag(iv)   = zqgk(iv)/zdtdh(iv) + zprvg(iv) + zpkg(iv)
        zzai(iv)   = zqik(iv)/zdtdh(iv) + zprvi(iv) + zpki(iv)
        
        zpkr(iv)   = MIN( zpkr(iv) , zzar(iv) )
        
        IF (llqs) THEN
          ic1 = ic1 + 1
          ivdx1(ic1) = iv
        ENDIF
        IF (llqr(iv)) THEN
          ic2 = ic2 + 1
          ivdx2(ic2) = iv
       ENDIF
       IF (llqg) THEN
         ic3 = ic3 + 1
         ivdx3(ic3) = iv
       ENDIF
       IF (llqi) THEN
         ic4 = ic4 + 1
         ivdx4(ic4) = iv
       ENDIF

    
     ENDDO

!CDIR BEGIN COLLAPSE
  zpkr  (:) = 0.0_ireals
  zpks  (:) = 0.0_ireals
  zpkg  (:) = 0.0_ireals
  zpki  (:) = 0.0_ireals
!CDIR END

  zln1o2=EXP (ccswxp * LOG (0.5_ireals))
     
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
        alf = EXP(hlp*zlog_10) ! 10.0_ireals**hlp
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
      zlnqsk = zvz0s(iv) * EXP (ccswxp * LOG (zqsk(iv))) * zrho1o2(iv)
      zpks  (iv) = zqsk(iv) * zlnqsk
      IF (zvzs(iv) == 0.0_ireals) THEN
        zvzs(iv) = zlnqsk * zln1o2
      ENDIF
    ENDDO loop_over_qs_prepare
    
! sedimentation fluxes
  zln1o2=EXP (zvzxp * LOG (0.5_ireals))    
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
    loop_over_qr_sedi: DO i1d = 1, ic2
      iv = ivdx2(i1d)

      zlnqrk = zvz0r * EXP (zvzxp * LOG (zqrk(iv))) * zrho1o2(iv)
      zpkr(iv) = zqrk(iv) * zlnqrk
      IF (zvzr(iv) == 0.0_ireals) THEN
        zvzr(iv) = zlnqrk * zln1o2
      ENDIF
    ENDDO loop_over_qr_sedi

  zln1o2=EXP (zexpsedg * LOG (0.5_ireals))   
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
    loop_over_qg_sedi: DO i1d = 1, ic3
      iv = ivdx3(i1d)
      zlnqgk = zvz0g * EXP (zexpsedg * LOG (zqgk(iv))) * zrho1o2(iv)
      zpkg(iv) = zqgk(iv) * zlnqgk
      IF (zvzg(iv) == 0.0_ireals) THEN
        zvzg(iv) = zlnqgk * zln1o2
      ENDIF
    ENDDO loop_over_qg_sedi

!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
    loop_over_qi_sedi: DO i1d = 1, ic4
      iv = ivdx4(i1d)
      
      IF (zvzi(iv) == 0.0_ireals) THEN
        !! density correction not needed zrho1o2(iv)
        zvzi(iv) = zvz0i * EXP (zbvi  * LOG (0.5_ireals*zqik(iv)))
      ENDIF
    ENDDO loop_over_qi_sedi

  !----------------------------------------------------------------------------
  ! 2.3: Second part of preparations
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
    zeln7o8qrk   (:) = 0.0_ireals
    zeln7o4qrk   (:) = 0.0_ireals 
    zeln27o16qrk (:) = 0.0_ireals
    zeln13o8qrk  (:) = 0.0_ireals
    zeln3o4qsk   (:) = 0.0_ireals
    zeln8qsk     (:) = 0.0_ireals
    zeln6qgk     (:) = 0.0_ireals
    zelnrimexp_g (:) = 0.0_ireals
    zsrmax       (:) = 0.0_ireals
    zssmax       (:) = 0.0_ireals
    zsgmax       (:) = 0.0_ireals
    
!FR old  
!   zcsdep       (:) = 3.2E-2_ireals
    zcsdep       (:) = 3.367E-2_ireals     
    zcidep       (:) = 1.3E-5_ireals
    zcslam       (:) = 1e10_ireals

    scau         (:) = 0.0_ireals
    scac         (:) = 0.0_ireals
    snuc         (:) = 0.0_ireals
    scfrz        (:) = 0.0_ireals
    simelt       (:) = 0.0_ireals
    sidep        (:) = 0.0_ireals
    ssdep        (:) = 0.0_ireals
    sgdep        (:) = 0.0_ireals
    sdau         (:) = 0.0_ireals
    srim         (:) = 0.0_ireals
    srim2        (:) = 0.0_ireals
    sshed        (:) = 0.0_ireals
    sicri        (:) = 0.0_ireals
    srcri        (:) = 0.0_ireals
    sagg         (:) = 0.0_ireals
    sagg2        (:) = 0.0_ireals
    siau         (:) = 0.0_ireals
    ssmelt       (:) = 0.0_ireals
    sgmelt       (:) = 0.0_ireals
    sev          (:) = 0.0_ireals
    sconr        (:) = 0.0_ireals
    sconsg       (:) = 0.0_ireals
    srfrz        (:) = 0.0_ireals
!CDIR END

    ! Nullify counters again
    ic1 = 0
    ic2 = 0
    ic3 = 0
    ic4 = 0
    ic5 = 0
    ic6 = 0
    ic7 = 0
    ic8 = 0
    
    DO iv = iv_start, iv_end

        qrg  = qr(iv,k)
        qsg  = qs(iv,k)
        qgg  = qg(iv,k)
        qvg  = qv(iv,k)
        qcg  = qc(iv,k)
        qig  = qi(iv,k)
        tg   =  t(iv,k)
        ppg  =  p(iv,k)
        rhog = rho(iv,k)
        
        zpkr(iv)   = MIN( zpkr(iv) , zzar(iv) )
        zpks(iv)   = MIN( zpks(iv) , zzas(iv) )
        zpkg(iv)   = MIN( zpkg(iv) , zzag(iv) )
        zpki(iv)   = MIN( zpki(iv) , zzai(iv) )

        zzar(iv)   = zdtdh(iv) * (zzar(iv)-zpkr(iv))
        zzas(iv)   = zdtdh(iv) * (zzas(iv)-zpks(iv))
        zzag(iv)   = zdtdh(iv) * (zzag(iv)-zpkg(iv))
        zzai(iv)   = zdtdh(iv) * (zzai(iv)-zpki(iv))

        zimr(iv)   = 1.0_ireals / (1.0_ireals + zvzr(iv) * zdtdh(iv))
        zims(iv)   = 1.0_ireals / (1.0_ireals + zvzs(iv) * zdtdh(iv))
        zimg(iv)   = 1.0_ireals / (1.0_ireals + zvzg(iv) * zdtdh(iv))
        zimi(iv)   = 1.0_ireals / (1.0_ireals + zvzi(iv) * zdtdh(iv))

        zqrk(iv)   = zzar(iv)*zimr(iv)
        zqsk(iv)   = zzas(iv)*zims(iv)
        zqgk(iv)   = zzag(iv)*zimg(iv)
        zqik(iv)   = zzai(iv)*zimi(iv)
        
        zqvsi(iv) = sat_pres_ice(tg)/(rhog * r_v * tg)
        
        llqr(iv) = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin
        llqg = zqgk(iv) > zqmin
       ! llqi = zqik(iv) > zqmin
        llqc =       qcg > zqmin
        llqi =       qig > zqmin

        IF (llqr(iv)) THEN
          ic1 = ic1 + 1
          ivdx1(ic1) = iv
        ENDIF
        IF (llqs) THEN
          ic2 = ic2 + 1
          ivdx2(ic2) = iv
        ENDIF
        IF (llqg) THEN
          ic3 = ic3 + 1
          ivdx3(ic3) = iv
        ENDIF
        IF (llqi .OR. llqs) THEN
          ic4 = ic4 + 1
          ivdx4(ic4) = iv
        ENDIF
        IF ( tg < zthet .AND. qvg >  8.E-6_ireals &
                        .AND. qig <= 0.0_ireals ) THEN
          ic5 = ic5 + 1
          ivdx5(ic5) = iv
        ENDIF
        IF (llqc) THEN
          ic6 = ic6 + 1
          ivdx6(ic6) = iv
        ENDIF
        IF (llqi .OR. llqs .OR. llqg ) THEN
          ic7 = ic7 + 1
          ivdx7(ic7) = iv
        ENDIF
        IF (llqr(iv) .AND. qcg <= 0.0_ireals) THEN
          ic8 = ic8 + 1
          ivdx8(ic8) = iv
        ENDIF
    
     ENDDO

  !!----------------------------------------------------------------------------
  !! 2.4: IF (llqr): ic1
  !!----------------------------------------------------------------------------
     
! ic1
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qr: DO i1d =1, ic1
      iv = ivdx1(i1d)

      rhog = rho(iv,k)
      qcg  = qc(iv,k)
      qig  = qi(iv,k)
      tg   =  t(iv,k)
      llqi = qig > zqmin

      zlnqrk       = LOG (zqrk(iv))
      zsrmax(iv)   = zzar(iv)/rhog*zdtr
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
    ENDDO loop_over_qr

  !!----------------------------------------------------------------------------
  !! 2.5: IF (llqs): ic2
  !!----------------------------------------------------------------------------
    
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
    loop_over_qs_coeffs: DO i1d =1, ic2
      iv = ivdx2(i1d)
      
      rhog = rho(iv,k)
      qcg  = qc(iv,k)
      qig  = qi(iv,k)
      
      zlnqsk       = LOG (zqsk(iv))
      zssmax(iv)   = zzas(iv) / rhog*zdtr
      IF (qig+qcg > zqmin) THEN
        zeln3o4qsk(iv) = EXP (x3o4 *zlnqsk)
      ENDIF
      zeln8qsk(iv) = EXP (0.8_ireals *zlnqsk)
    ENDDO loop_over_qs_coeffs

  !!----------------------------------------------------------------------------
  !! 2.6: IF (llqg): ic3
  !!----------------------------------------------------------------------------
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
    loop_over_qg_coeffs: DO i1d =1, ic3
      iv = ivdx3(i1d)
      
      rhog = rho(iv,k)
      qcg  = qc(iv,k)
      qig  = qi(iv,k)

      zlnqgk       = LOG (zqgk(iv))
      zsgmax(iv)   = zzag(iv) / rhog*zdtr
      IF (qig+qcg > zqmin) THEN
        zelnrimexp_g(iv) = EXP (zrimexp_g * zlnqgk)
      ENDIF
      zeln6qgk(iv) = EXP (0.6_ireals *zlnqgk)
    ENDDO loop_over_qg_coeffs

      
  !!----------------------------------------------------------------------------
  !! 2.7:  slope of snow PSD and coefficients for depositional growth (llqi,llqs)
  !!----------------------------------------------------------------------------    
    
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qi_qs: DO i1d =1, ic4
      iv = ivdx4(i1d)

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

    ENDDO loop_over_qi_qs

  
  !!----------------------------------------------------------------------------
  !! 2.8: Deposition nucleation for low temperatures below a threshold (llqv)
  !!----------------------------------------------------------------------------    

!CDIR NODEP,VOVERTAKE,VOB
    loop_over_icenucleation: DO i1d = 1, ic5
      iv = ivdx5(i1d)

      qvg  =  qv(iv,k)
      tg   =   t(iv,k)
      rhog = rho(iv,k)

      IF( qvg > zqvsi(iv) ) THEN
        IF( lsuper_coolw) THEN
          znin  = MIN( fxna_cooper(tg), znimax )
        ELSE
          znin  = MIN( fxna(tg), znimax )
        END IF
        snuc(iv) = zmi0 * z1orhog(iv) * znin * zdtr
      ENDIF

    ENDDO loop_over_icenucleation

  !!--------------------------------------------------------------------------
  !! Section 3: Search for cloudy grid points with cloud water and
  !!            calculation of the conversion rates involving qc (ic6)
  !!--------------------------------------------------------------------------

!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qc: DO i1d =1, ic6
      iv = ivdx6(i1d)

      qrg  =   qr(iv,k)
      qcg  =   qc(iv,k)
      qig  =   qi(iv,k)    
      qcgk_1 = qc(iv,k-1) 
      tg   =    t(iv,k)
      ppg  =    p(iv,k)
      rhog =  rho(iv,k)
      llqs = zqsk(iv) > zqmin
      llqi =  qig > zqmin 

      zscmax = qcg*zdtr
      IF( tg > zthn ) THEN
        IF (iautocon == 0) THEN
          ! Kessler (1969) autoconversion rate
          scau(iv) = zccau * MAX( qcg - qc0, 0.0_ireals )
          scac(iv) = zcac  * qcg * zeln7o8qrk(iv)
        ELSEIF (iautocon == 1) THEN
          ! Seifert and Beheng (2001) autoconversion rate
          ! with constant cloud droplet number concentration cloud_num
          IF (qcg > 1e-6) THEN
            ztau  = MIN(1.0_ireals-qcg/(qcg+qrg),0.9_ireals)
            zphi  = zkphi1 * ztau**zkphi2 * (1.0_ireals - ztau**zkphi2)**3
            scau(iv) = zconst * qcg*qcg*qcg*qcg &
                      * (1.0_ireals + zphi/(1.0_ireals - ztau)**2)
            zphi      = (ztau/(ztau+zkphi3))**4
            scac(iv) = zkcac * qcg * qrg * zphi
          ENDIF
        ENDIF
        IF (llqr(iv)) THEN
          ! Calculation of in-cloud rainwater freezing
          IF ( tg < ztrfrz ) THEN
            IF (lsuper_coolw) THEN
              srfrz(iv) = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_ireals ) * zeln7o4qrk(iv)
            ELSE
              ztfrzdiff=ztrfrz-tg
              srfrz(iv) = zcrfrz*ztfrzdiff*SQRT(ztfrzdiff)* zeln27o16qrk(iv)
            ENDIF
          ENDIF
        ENDIF
        IF (llqs) THEN
          srim(iv) = zcrim(iv) * qcg *  EXP(ccsaxp * LOG(zcslam(iv)))
        ENDIF
        srim2(iv)  = zcrim_g * qcg * zelnrimexp_g(iv)
        IF( tg >= t0 ) THEN
          sshed(iv) = srim(iv)+srim2(iv)
          srim (iv) = 0.0_ireals
          srim2(iv) = 0.0_ireals
        ELSE
          IF (qcg.GE.qc0) THEN
            sconsg(iv) = zcsg * qcg * zeln3o4qsk(iv)
          ENDIF
        ENDIF
        ! Check for maximum depletion of cloud water and adjust the
        ! transfer rates accordingly
        zscsum = scau(iv) + scac(iv) + srim(iv) + srim2(iv) + sshed(iv)
        zcorr  = zscmax / MAX( zscmax, zscsum )
        scau   (iv) = zcorr*scau(iv)
        scac   (iv) = zcorr*scac(iv)
        srim   (iv) = zcorr*srim(iv)
        srim2  (iv) = zcorr*srim2(iv)
        sshed  (iv) = zcorr*sshed(iv)
        sconsg (iv) = MIN (sconsg(iv), srim(iv)+zssmax(iv))
      ELSE !tg >= tg: ! hom. freezing of cloud and rain water
        scfrz(iv) = zscmax
        srfrz(iv) = zsrmax(iv)
      ENDIF
      ! Calculation of heterogeneous nucleation of cloud ice.
      ! This is done in this section, because we require water saturation
      ! for this process (i.e. the existence of cloud water) to exist.
      ! Heterogeneous nucleation is assumed to occur only when no
      ! cloud ice is present and the temperature is below a nucleation
      ! threshold.
      IF( tg <= 267.15_ireals .AND. .NOT.llqi ) THEN   
        IF (lsuper_coolw) THEN
          znin  = MIN( fxna_cooper(tg), znimax )
          snuc(iv) = zmi0 * z1orhog(iv) * znin * zdtr
        ELSE
          znin      = MIN( fxna(tg), znimax )
          snuc(iv) = zmi0 / rhog * znin * zdtr
        END IF
      ENDIF
!FR>>> Calculation of reduction of depositional growth at cloud top (Forbes 2012)
      IF( k>1 .AND. k<ke .AND. lsuper_coolw ) THEN
        znin = MIN( fxna_cooper(tg), znimax )
        fnuc = MIN(znin/znimix, 1.0_ireals)
        
        !! distance from cloud top
        IF( qcgk_1 .LT. zqmin ) THEN      ! upper cloud layer
          dist_cldtop(iv) = 0.0_ireals    ! reset distance to upper cloud layer
        ELSE
          dist_cldtop(iv) = dist_cldtop(iv) + dz(iv,k)
        END IF
         
! with asymptotic behaviour dz -> 0 (xxx)
!        reduce_dep(iv) = MIN(fnuc + (1.0_ireals-fnuc)*(reduce_dep_ref + &
!                             dist_cldtop(iv)/dist_cldtop_ref + &
!                             (1.0_ireals-reduce_dep_ref)*(zdh/dist_cldtop_ref)**4), 1.0_ireals)
        
! without asymptotic behaviour dz -> 0
        reduce_dep(iv) = MIN(fnuc + (1.0_ireals-fnuc)*(reduce_dep_ref + &
                          dist_cldtop(iv)/dist_cldtop_ref), 1.0_ireals)
       
      END IF ! Reduction of dep. growth of snow/ice 

      
    ENDDO loop_over_qc
    
      !------------------------------------------------------------------------
      ! Section 4: Search for cold grid points with cloud ice and/or snow and
      !            calculation of the conversion rates involving qi, qs and qg
      !------------------------------------------------------------------------

! ic7
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs_qi_qg: DO i1d =1, ic7
      iv = ivdx7(i1d)

      qvg  =  qv(iv,k)
      qig  =  qi(iv,k)
      tg   =   t(iv,k)
      ppg  =   p(iv,k)
      rhog = rho(iv,k)
      llqs =  zqsk(iv) > zqmin
      llqi =  qig > zqmin

      IF (tg<=t0) THEN           ! cold case 

        zqvsidiff = qvg-zqvsi(iv)
        zsvmax    = zqvsidiff * zdtr
        IF (llqi) THEN

          ! Change in sticking efficiency needed in case of cloud ice sedimentation
          ! (based on Guenther Zaengls work)
          IF (lstickeff) THEN
            zeff     = MIN(EXP(0.09_ireals*(tg-t0)),1.0_ireals)
            zeff     = MAX(zeff, zceff_min, zceff_fac*(tg-tmin_iceautoconv)) 
          ELSE !original sticking efficiency of cloud ice
            zeff     = MIN(EXP(0.09_ireals*(tg-t0)),1.0_ireals)
            zeff     = MAX(zeff,0.2_ireals)
          END IF
          sagg(iv)   = zeff * qig * zcagg(iv) * EXP(ccsaxp*LOG(zcslam(iv)))
          sagg2(iv)  = zeff * qig * zcagg_g * zelnrimexp_g(iv)
          siau(iv)   = zeff * zciau * MAX( qig - qi0, 0.0_ireals )
          znin       = MIN( fxna(tg), znimax )
          zmi        = MIN( rhog*qig/znin, zmimax )
          zmi        = MAX( zmi0, zmi )
          znid       = rhog * qig/zmi
          zlnlogmi   = LOG (zmi)
          sidep(iv)  = zcidep(iv) * znid * EXP(0.33_ireals * zlnlogmi) * zqvsidiff
          zsvidep    = 0.0_ireals
          zsvisub    = 0.0_ireals
          zsimax     = qig*zdtr
          IF( sidep(iv) > 0.0_ireals ) THEN
            IF (lsuper_coolw ) THEN
              sidep(iv) = sidep(iv) * reduce_dep(iv)  !FR new: SLW reduction
            END IF
            zsvidep = MIN( sidep(iv), zsvmax )
          ELSEIF ( sidep(iv) < 0.0_ireals ) THEN
            zsvisub  =   MAX ( sidep(iv),  zsvmax)
            zsvisub  = - MAX (    zsvisub, -zsimax)
          ENDIF
          zlnlogmi   = LOG  (zmsmin/zmi)
          zztau      = 1.5_ireals*( EXP(0.66_ireals*zlnlogmi) - 1.0_ireals)
          sdau(iv)  = zsvidep/zztau
          sicri(iv) = zcicri * qig * zeln7o8qrk(iv)
          srcri(iv) = zcrcri * (qig/zmi) * zeln13o8qrk(iv)
        ELSE
          zsimax    =  0.0_ireals
          zsvidep   =  0.0_ireals
          zsvisub   =  0.0_ireals
        ENDIF
        
        zxfac      = 1.0_ireals + zbsdep(iv) * EXP(ccsdxp*LOG(zcslam(iv)))
        ssdep(iv) = zcsdep(iv) * zxfac * zqvsidiff / (zcslam(iv)+zeps)**2
        !FR new: SLW reduction
        IF (lsuper_coolw .AND. ssdep(iv) > 0.0_ireals) THEN
          ssdep(iv) = ssdep(iv)*reduce_dep(iv)
        END IF
        sgdep(iv) = (0.398561_ireals-0.00152398_ireals*tg                 &
                     + 2554.99_ireals/ppg+ 2.6531E-7_ireals*ppg) *        &
                     zqvsidiff * zeln6qgk(iv)
        ! Check for maximal depletion of cloud ice
        ! No check is done for depositional autoconversion because
        ! this is a always a fraction of the gain rate due to
        ! deposition (i.e the sum of this rates is always positive)
        zsisum = siau(iv) + sagg(iv) + sagg2(iv) + sicri(iv) + zsvisub
        zcorr  = 0.0_ireals
        IF( zsimax > 0.0_ireals ) zcorr  = zsimax / MAX( zsimax, zsisum )
        sidep(iv)  = zsvidep - zcorr*zsvisub
        siau (iv)  = zcorr*siau(iv)
        sagg (iv)  = zcorr*sagg(iv)
        sagg2(iv)  = zcorr*sagg2(iv)
        sicri(iv)  = zcorr*sicri(iv)
        IF ( zqvsidiff < 0.0_ireals ) THEN
          ssdep(iv) = MAX(ssdep(iv), - zssmax(iv))
          sgdep(iv) = MAX(sgdep(iv), - zsgmax(iv))
        ENDIF
        
      !------------------------------------------------------------------------
      ! Section 5: Search for warm grid points with cloud ice and/or snow and
      !            calculation of the melting rates of qi and ps
      !------------------------------------------------------------------------

      ELSE ! tg > 0 - warm case
        
        ! cloud ice melts instantaneously
        simelt(iv) = qig*zdtr

        zqvsw0     = fqvs( zpvsw0, ppg)
        zqvsw0diff = qvg-zqvsw0
        
        IF ( tg > (t0-ztcrit*zqvsw0diff) ) THEN
          !calculate melting rate
          zx1         = (tg - t0) + zasmel*zqvsw0diff
          ssmelt(iv) = (79.6863_ireals/ppg+0.612654E-3_ireals)* zx1 * zeln8qsk(iv)
          ssmelt(iv) = MIN (ssmelt(iv),zssmax(iv))
          sgmelt(iv) = (12.31698_ireals/ppg+7.39441e-05_ireals)* zx1 * zeln6qgk(iv)
          sgmelt(iv) = MIN (sgmelt(iv), zsgmax(iv))
          !deposition + melting, ice particle temperature: t0
          !calculation without howell-factor!
          ssdep(iv)  = (31282.3_ireals/ppg+0.241897_ireals)       &
                      * zqvsw0diff * zeln8qsk(iv)
          sgdep(iv)  = (0.153907_ireals-ppg*7.86703e-07_ireals)  &
                      * zqvsw0diff * zeln6qgk(iv)
          IF (zqvsw0diff < 0.0_ireals) THEN
            !melting + evaporation of snow/graupel
            ssdep(iv) = MAX (-zssmax(iv),ssdep(iv))
            sgdep(iv) = MAX (-zsgmax(iv),sgdep(iv))
            !melt water evaporates
            ssmelt(iv) = ssmelt(iv)+ssdep(iv)
            sgmelt(iv) = sgmelt(iv)+sgdep(iv)
            ssmelt(iv) = MAX( ssmelt(iv), 0.0_ireals )
            sgmelt(iv) = MAX( sgmelt(iv), 0.0_ireals )
          ELSE
            !deposition on snow/graupel is interpreted as increase
            !in rain water ( qv --> qr, sconr)
            !therefore,  sconr=(zssdep+zsgdep)
            sconr(iv)=ssdep(iv)+sgdep(iv)
            ssdep(iv)=0.0_ireals
            sgdep(iv)=0.0_ireals
          ENDIF
        ELSE
          !if t<t_crit
          !no melting, only evaporation of snow/graupel
          zqvsw      = fqvs( fpvsw(tg), ppg )
          zqvsidiff  = qvg-zqvsw
          ssdep(iv) = (0.28003_ireals-ppg*0.146293E-6_ireals) &
                       * zqvsidiff * zeln8qsk(iv)
          sgdep(iv) = (0.0418521_ireals-ppg*4.7524E-8_ireals) &
                       * zqvsidiff *zeln6qgk(iv)
          ssdep(iv) = MAX(-zssmax(iv) ,ssdep(iv) )
          sgdep(iv) = MAX(-zsgmax(iv) ,sgdep(iv) )
        ENDIF !t_crit
      ENDIF !tg

    ENDDO loop_over_qs_qi_qg

    !--------------------------------------------------------------------------
    ! Section 6: Search for grid points with rain in subsaturated areas
    !            and calculation of the evaporation rate of rain
    !--------------------------------------------------------------------------

! ic8
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qr_nocloud: DO i1d =1, ic8
      iv = ivdx8(i1d)

      qvg = qv(iv,k)
      tg  =  t(iv,k)
      rhog = rho(iv,k)
      
      zqvsw    = sat_pres_water(tg)/(rhog * r_v *tg)
      zlnqrk   = LOG (zqrk(iv))
      zx1      = 1.0_ireals + zbev * EXP (zbevxp  * zlnqrk)
      !sev(iv)  = zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk)
      ! Limit evaporation rate in order to avoid overshoots towards supersaturation
      ! the pre-factor approximates (esat(T_wb)-e)/(esat(T)-e) at temperatures between 0 degC and 30 degC
      temp_c = tg - t0
      maxevap     = (0.61_ireals-0.0163_ireals*temp_c+1.111e-4_ireals*temp_c**2)*(zqvsw-qvg)/zdt
      sev(iv)    = MIN(zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk), maxevap)
      
      IF( tg > zthn ) THEN
        ! Calculation of below-cloud rainwater freezing
        IF ( tg < ztrfrz ) THEN
          IF (lsuper_coolw) THEN
            !FR new: reduced rain freezing rate
            srfrz(iv) = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_ireals ) * zeln7o4qrk(iv)
          ELSE
            srfrz(iv) = zcrfrz*SQRT( (ztrfrz-tg)**3 ) * zeln27o16qrk(iv)
          ENDIF
        ENDIF
      ELSE ! Hom. freezing of rain water
        srfrz(iv) = zsrmax(iv)
      ENDIF
        
    ENDDO loop_over_qr_nocloud

    !--------------------------------------------------------------------------
    ! Section 7: Calculate the total tendencies of the prognostic variables.
    !            Update the prognostic variables in the interior domain.
    !--------------------------------------------------------------------------

    loop_over_all_iv: DO iv = iv_start, iv_end
      
      qrg = qr(iv,k)
      qsg = qs(iv,k)
      qig = qi(iv,k)
      rhog = rho(iv,k)
      
      zsrsum = sev(iv) + srfrz(iv) + srcri(iv)
      zcorr  = 1.0_ireals
      IF(zsrsum > 0) THEN
        zcorr  = zsrmax(iv) / MAX( zsrmax(iv), zsrsum )
      ENDIF
      sev(iv)   = zcorr*sev(iv)
      srfrz(iv) = zcorr*srfrz(iv)
      srcri(iv) = zcorr*srcri(iv)
      
      zqvt =   sev   (iv) - sidep (iv) - ssdep (iv) - sgdep (iv)         &
             - snuc  (iv) - sconr (iv)
      zqct =   simelt(iv) - scau  (iv) - scfrz (iv) - scac  (iv)         &
             - sshed (iv) - srim  (iv) - srim2 (iv)
      zqit =   snuc  (iv) + scfrz (iv) - simelt(iv) - sicri (iv)         &
             + sidep (iv) - sdau  (iv) - sagg  (iv) - sagg2 (iv) - siau(iv)
      zqrt =   scau  (iv) + sshed (iv) + scac  (iv) + ssmelt(iv)         &
             + sgmelt(iv) - sev   (iv) - srcri (iv) - srfrz (iv) + sconr(iv)
      zqst =   siau  (iv) + sdau  (iv) - ssmelt(iv) + srim  (iv)         &
              + ssdep (iv) + sagg  (iv) - sconsg(iv)
      zqgt =   sagg2 (iv) - sgmelt(iv) + sicri (iv) + srcri (iv)         &
             + sgdep (iv) + srfrz (iv) + srim2 (iv) + sconsg(iv)
      
      ztt = z_heat_cap_r*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst+zqgt) )

      ! Update variables and add qi to qrs for water loading
      IF (lsedi_ice ) THEN
        qig = MAX ( 0.0_ireals, (zzai(iv)*z1orhog(iv) + zqit*zdt)*zimi(iv))
      ELSE
        qig = MAX ( 0.0_ireals, qig + zqit*zdt)
      END IF
      qrg = MAX ( 0.0_ireals, (zzar(iv)/rhog + zqrt*zdt)*zimr(iv))
      qsg = MAX ( 0.0_ireals, (zzas(iv)/rhog + zqst*zdt)*zims(iv))
      qgg = MAX ( 0.0_ireals, (zzag(iv)/rhog + zqgt*zdt)*zimg(iv))

      
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
        IF (zprvr(iv) .LE. zqmin) zprvr(iv)=0.0_ireals
        IF (zprvs(iv) .LE. zqmin) zprvs(iv)=0.0_ireals
        IF (zprvg(iv) .LE. zqmin) zprvg(iv)=0.0_ireals
        IF (zprvi(iv) .LE. zqmin) zprvi(iv)=0.0_ireals


#ifdef NUDGING
          ! for the latent heat nudging
        IF ((llhn .OR. llhnverif) .AND. lhn_qrs ) THEN
          IF (lsedi_ice) THEN
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)+zprvi(iv)+zprvg(iv)
            qrsflux(iv,k) = 0.5*(qrsflux(iv,k)+zpkr(iv)+zpks(iv)+zpkg(iv)+zpki(iv))
          ELSE 
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)+zprvg(iv)
            qrsflux(iv,k) = 0.5*(qrsflux(iv,k)+zpkr(iv)+zpks(iv)+zpkg(iv))
          END IF
        ENDIF
#endif

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
          IF (qgg+qg(iv,k+1) <= zqmin ) THEN
            zvzg(iv)= 0.0_ireals
          ELSE
            zvzg(iv)=zvz0g * EXP(zexpsedg*LOG((qgg+qg(iv,k+1))*0.5_ireals*rhog))* zrho1o2(iv)
          ENDIF
          IF (qig+qi(iv,k+1) <= zqmin ) THEN
            zvzi(iv)= 0.0_ireals
          ELSE
            !! density correction not needed
            zvzi(iv)= zvz0i * EXP(zbvi*LOG((qig+qi(iv,k+1))*0.5_ireals*rhog))
          ENDIF
          
        ELSE
          ! Precipitation fluxes at the ground
          prr_gsp(iv) = 0.5_ireals * (qrg*rhog*zvzr(iv) + zpkr(iv))
          IF (lsedi_ice) THEN
            prs_gsp(iv) = 0.5_ireals * (rhog*(qsg*zvzs(iv)+qig*zvzi(iv)) + zpks(iv)+zpki(iv))
          ELSE
            prs_gsp(iv) = 0.5_ireals * (qsg*rhog*zvzs(iv) + zpks(iv))
          END IF
          prg_gsp(iv) = 0.5_ireals * (qgg*rhog*zvzg(iv) + zpkg(iv))

          
#ifdef NUDGING
          ! for the latent heat nudging
          IF ((llhn .OR. llhnverif) .AND. lhn_qrs)        &
            qrsflux(iv,k) = prr_gsp(iv)+prs_gsp(iv)+prg_gsp(iv)
#endif

        ENDIF

        ! Update of prognostic variables or tendencies
        qr (iv,k) = MAX ( 0.0_ireals, qrg )
        qs (iv,k) = MAX ( 0.0_ireals, qsg )
        qi (iv,k) = MAX ( 0.0_ireals, qig )
        qg (iv,k) = MAX ( 0.0_ireals, qgg )
        t  (iv,k) = t (iv,k) + ztt*zdt 
        qv (iv,k) = MAX ( 0.0_ireals, qv(iv,k) + zqvt*zdt )
        qc (iv,k) = MAX ( 0.0_ireals, qc(iv,k) + zqct*zdt )

      ENDDO loop_over_all_iv

  IF (izdebug > 15) THEN
    ! Check for negative values
     DO iv = iv_start, iv_end
        IF (qr(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_gr, negative value in qr'
          CALL message('',message_text)
        ENDIF
        IF (qc(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_gr, negative value in qc'
          CALL message('',message_text)
        ENDIF
        IF (qi(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_gr, negative value in qi'
          CALL message('',message_text)
        ENDIF
        IF (qs(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_gr, negative value in qs'
          CALL message('',message_text)
        ENDIF
        IF (qv(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_gr, negative value in qv'
          CALL message('',message_text)
        ENDIF
        IF (qg(iv,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp_gr, negative value in qg'
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

#ifdef NUDGING
! add part of latent heating calculated in subroutine hydci to model latent
! heating field: add temperature to model latent heating field
IF (llhn .OR. llhnverif) &
 !CALL get_gs_lheating ('inc',1,ke)  !this should be called from within the block
 tt_lheat(:,:) = tt_lheat(:,:) + t(:,:)
#endif

#ifdef __ICON__

CALL satad_v_3d (                             &
               & maxiter  = 10_iintegers ,& !> IN
               & tol      = 1.e-3_ireals ,& !> IN
               & te       = t            ,&
               & qve      = qv           ,&
               & qce      = qc           ,&
               & rhotot   = rho          ,&
               & idim     = nvec         ,&
               & kdim     = ke           ,&
               & ilo      = iv_start     ,&
               & iup      = iv_end       ,&
               & klo      = k_start      ,&
               & kup      = ke            &
               )

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
    CALL message('gscp_hydci_pp_gr', 'UPDATED VARIABLES')
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_gr T= ',&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_gr qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_gr qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_gr qi= ',&
    MAXVAL( qi(:,:)), MINVAL(qi(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_gr qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_gr qs= ',&
    MAXVAL( qs(:,:)), MINVAL(qs(:,:) )
   CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp_gr qg= ',&
    MAXVAL( qg(:,:)), MINVAL(qg(:,:) )
    CALL message('', TRIM(message_text))
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine hydci_pp_gr
!------------------------------------------------------------------------------

  END SUBROUTINE hydci_pp_gr


END MODULE gscp_hydci_pp
