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

USE src_lheating,             ONLY :  &
    get_gs_lheating            ! storage of grid scale latent heating for lhn
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

!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: hydci_pp, hydci_pp_init


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
  
  !! Other process coefficients.
  !! These coefficients have been precalculated on the basis of the
  !! following basic parameters (original hydci):
  !! N0R = 8.0 E6 : Parameter in the size distrubution function for rain
  !! ECR = 0.8    : Collection efficiency for rain collecting cloud water
  !! ECS = 0.9    : Collection efficiency for snow collecting cloud water
  !! EIR = 0.8    : Collection efficiency for rain collecting cloud ice
  !! V0R = 130.0  : Factor in the terminal velocity for raindrops
  !!                VTR(D) = V0R*D**(1/2)
  !! AMS = 0.038  : Formfactor in the mass-size relation of snowparticles
  !!                m(DS) = AMS*DS**2
  !! AMI = 130.0  : Formfactor in the mass-size relation of cloud ice crystals
  !!                m(DI) = AMI*DI**3
  !! ETA  = 1.75 E-5 : Viscosity of air    
  !! DIFF = 2.22 E-5 : Molecular diffusion coefficient for water vapour
  !! LHEAT= 2.40 E-2 : Heat conductivity of air
  !! RHO  = 1.0      : Density of air 
  !! RHOW = 1.0 E3   : Density of water 
  !! lh_v            : latent heat of vapourization
  !! lh_f            : latent heat of fusion
  !! AR  = PI*RHOW*N0R
  !! AS  = 2*N0S*AMS
  !! HW              : Howell factor ( =(1/(1+H_w)) or =(1/(1+H_i))) 
  
  zhw   = 2.270603,     & ! Howell factor
  zecs  = 0.9_ireals,   & ! Collection efficiency for snow collecting cloud water
  
  zadi  = 0.217_ireals, & ! Formfactor in the size-mass relation of ice particles
  zbdi  = 0.302_ireals, & ! Exponent in the size-mass relation of ice particles
  zams  = 0.069_ireals, & ! Formfactor in the mass-size relation of snow particles
  zbms  = 2.000_ireals, & ! Exponent in the mass-size relation of snow particles
  
  zv1s  = 0.50_ireals,  & ! Exponent in the terminal velocity for snow
  
  zami  = 130.0_ireals, & ! Formfactor in the mass-size relation of cloud ice
  zn0s0 = 8.0E5_ireals, & ! 
  zn0s1 = 13.5_ireals * 5.65E5_ireals, & ! parameter in N0S(T)
  zn0s2 = -0.107_ireals , & ! parameter in N0S(T), Field et al
  zcac  = 1.72_ireals   , & ! (15/32)*(PI**0.5)*(ECR/RHOW)*V0R*AR**(1/8)
  zcicri= 1.72_ireals   , & ! (15/32)*(PI**0.5)*(EIR/RHOW)*V0R*AR**(1/8)
  zcrcri= 1.24E-3_ireals, & ! (PI/24)*EIR*V0R*Gamma(6.5)*AR**(-5/8)
! zcev  = 3.1E-3_ireals , & ! 2*PI*DIFF*HW*N0R*AR**(-1/2)
! zbev  = 9.0_ireals    , & ! 0.26*sqrt(0.5*RHO*v0r/eta)*Gamma(2.75)*AR**(-3/16)
  zcsmel= 1.48E-4_ireals, & ! 4*LHEAT*N0S*AS**(-2/3)/(RHO*lh_f)
  zbsmel= 14.37_ireals  , & ! 0.26*sqrt(0.5*RHO*v0s/eta)*Gamma(21/8)*AS**(-5/24)
  zasmel= 2.31E3_ireals , & ! DIFF*lh_v*RHO/LHEAT
  zcrfrz= 1.68_ireals   , & ! coefficient for raindrop freezing
  
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
  zmi0   = 1.0E-12_ireals, & ! initial crystal mass for cloud ice nucleation
  zmimax = 1.0E-9_ireals , & ! maximum mass of cloud ice crystals   
  zmsmin = 3.0E-9_ireals , & ! initial mass of snow crystals        
  
  !! Constant exponents in the transfer rate equations
  x1o12  =  1.0_ireals/12.0_ireals  ,  x3o16  =  3.0_ireals/16.0_ireals, &
  x7o8   =  7.0_ireals/ 8.0_ireals  ,  x2o3   =  2.0_ireals/ 3.0_ireals, &
  x5o24  =  5.0_ireals/24.0_ireals  ,  x1o8   =  1.0_ireals/ 8.0_ireals, &
  x13o8  = 13.0_ireals/ 8.0_ireals  ,  x13o12 = 13.0_ireals/12.0_ireals, &
  x27o16 = 27.0_ireals/16.0_ireals  ,  x1o3   =  1.0_ireals/ 3.0_ireals, &
  x1o2   =  1.0_ireals/ 2.0_ireals

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

  zconst = zkcau / (20.0_ireals*zxstar*cloud_num*cloud_num) &
    * (zcnue+2.0_ireals)*(zcnue+4.0_ireals)/(zcnue+1.0_ireals)**2
  ccsrim = 0.25_ireals*pi*zecs*v0snow*gamma_fct(zv1s+3.0_ireals)
  ccsagg = 0.25_ireals*pi*v0snow*gamma_fct(zv1s+3.0_ireals)
  ccsdep = 0.26_ireals*gamma_fct((zv1s+5.0_ireals)/2.d0)*SQRT(0.5/zeta)
  ccsvxp = -(zv1s/(zbms+1.0_ireals)+1.0_ireals)
  ccsvel = zams*v0snow*gamma_fct(zbms+zv1s+1.0_ireals)&
    & *(zams*gamma_fct(zbms+1.0_ireals))**ccsvxp 
  ccsvxp = ccsvxp + 1.0_ireals
  ccslam = zams*gamma_fct(zbms+1.0_ireals)
  ccslxp = 1.0_ireals / (zbms+1.0_ireals)
  ccswxp = zv1s*ccslxp
  ccsaxp = -(zv1s+3.0_ireals)
  ccsdxp = -(zbms+1.0_ireals)/2.0_ireals
  ccshi1 = lh_s*lh_s/(zlheat*r_v)
  ccdvtp = 2.11E-5 * t0**(-1.94) * 101325.0
  ccidep = 4.0_ireals * zami**(-x1o3)
  zn0r   = 8e6 * EXP(3.2*mu_rain) * (0.01)**(-mu_rain)  ! empirical relation adapted from Ulbrich (1983)
! to tune the zn0r variable
  zn0r   = zn0r * rain_n0_factor

  zar    = pi*zrhow/6.0 * zn0r * gamma_fct(mu_rain+4.0) ! pre-factor in lambda
  zcevxp = (mu_rain+2.)/(mu_rain+4.)
  zcev   = 2.0*pi*zdv/zhw*zn0r*zar**(-zcevxp) * gamma_fct(mu_rain+2.0)
  zbevxp = (2.*mu_rain+5.5_ireals)/(2.*mu_rain+8.)-zcevxp
  zbev   = 0.26 * SQRT(0.5*zrho0*130.0_ireals/zeta)*zar**(-zbevxp) &
    &    * gamma_fct((2.0*mu_rain+5.5)/2.0) / gamma_fct(mu_rain+2.0)
  zvzxp  = 0.5/(mu_rain+4.0)
  zvz0r  = 130.0_ireals*gamma_fct(mu_rain+4.5)/gamma_fct(mu_rain+4.0)*zar**(-zvzxp)


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
    ztx  , zpv  , zpx ,& ! dummy arguments for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
!    zdt,               & ! timestep for integration (water / ice )
    zdtr ,             & ! reciprocal of timestep for integration
    zscau , zscac  , zscrim , zscshe, zsnuc , & ! local values of the
    zsiau , zsagg  , zsidep , zsicri, zsrcri, & ! transfer rates
    zsdau , zssdep , zssmelt, & ! defined below
    zsrfrz,                                   & !
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum, zsrmax,    & ! terms for limiting  total rain water depletion
    zssmax,            & ! term for limiting snow depletion
    znin,              & ! number of cloud ice crystals at nucleation
    znid,              & ! number of cloud ice crystals for deposition
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,   & ! terms for limiting total
    zqvsi, zqvsw,      & ! sat. specitic humidity at ice and water saturation
    zztau, zxfac, zx1, zx2, ztt, &   ! some help variables
    ztau, zphi, zhi, zdvtp, ztc

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
    zlnlogmi,               & !
    qcg,tg,qvg,qrg, qsg,qig,rhog,ppg, alf,bet,m2s,m3s,hlp

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
    zvzr        (nvec),     & !
    zvzs        (nvec),     & !
    zpkr        (nvec),     & !
    zpks        (nvec),     & !
    zprvr       (nvec),     & !
    zprvs       (nvec),     & !
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
    zzar        (nvec),     & !
    zzas        (nvec),     & !
    zqrk        (nvec),     & !
    zqsk        (nvec),     & !
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
    srfrz  (nvec)    ! transfer rate due to rainwater freezing

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

! Delete precipitation fluxes from previous timestep
!CDIR BEGIN COLLAPSE
    prr_gsp (:) = 0.0_ireals
    prs_gsp (:) = 0.0_ireals
    zpkr    (:) = 0.0_ireals
    zpks    (:) = 0.0_ireals
    zprvr   (:) = 0.0_ireals
    zprvs   (:) = 0.0_ireals
    zvzr    (:) = 0.0_ireals
    zvzs    (:) = 0.0_ireals
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
    CALL get_gs_lheating ('add',1,ke)
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
!CDIR END

    ic1 = 0
    ic2 = 0

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

        llqr = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin

        zdtdh(iv) = 0.5_ireals * zdt / dz(iv,k)

        zzar(iv)   = zqrk(iv)/zdtdh(iv) + zprvr(iv) + zpkr(iv)
        zzas(iv)   = zqsk(iv)/zdtdh(iv) + zprvs(iv) + zpks(iv)

        IF (llqs) THEN
          ic1 = ic1 + 1
          ivdx1(ic1) = iv
        ENDIF
        IF (llqr) THEN
          ic2 = ic2 + 1
          ivdx2(ic2) = iv
       ENDIF
    
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

        llqr = zqrk(iv) > zqmin
        llqs = zqsk(iv) > zqmin

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

        zpkr(iv)   = MIN( zpkr(iv) , zzar(iv) )
        zpks(iv)   = MIN( zpks(iv) , zzas(iv) )

        zzar(iv)   = zdtdh(iv) * (zzar(iv)-zpkr(iv))
        zzas(iv)   = zdtdh(iv) * (zzas(iv)-zpks(iv))

        zimr(iv)   = 1.0_ireals / (1.0_ireals + zvzr(iv) * zdtdh(iv))
        zims(iv)   = 1.0_ireals / (1.0_ireals + zvzs(iv) * zdtdh(iv))

        zqrk(iv)   = zzar(iv)*zimr(iv)
        zqsk(iv)   = zzas(iv)*zims(iv)

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
      zqvsi  = fqvs( fpvsi(tg), ppg )
      zdvtp  = ccdvtp * EXP(1.94_ireals * LOG(tg)) / ppg          
      zhi    = ccshi1*zdvtp*rhog*zqvsi/(tg*tg)
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
      ppg  =   p(iv,k)
      rhog = rho(iv,k)

      zqvsi   = fqvs( fpvsi(tg), ppg )
      IF( qvg > zqvsi ) THEN
        znin  = MIN( fxna(tg), znimax )
        zsnuc = zmi0 * z1orhog(iv) * znin * zdtr
        snuc(iv) = zsnuc
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
      IF( tg <= 267.15_ireals .AND. qig <= 0.0_ireals ) THEN
        znin  = MIN( fxna(tg), znimax )
        zsnuc = zmi0 * z1orhog(iv) * znin * zdtr
        snuc(iv) = zsnuc
      ENDIF
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
      qsg  =  qs(iv,k)
      tg   =   t(iv,k)
      ppg  =   p(iv,k)
      rhog = rho(iv,k)
      llqi =  qig > zqmin

      IF (tg<=t0) THEN
        zqvsi   = fqvs( fpvsi(tg), ppg )
        znin    = MIN( fxna(tg), znimax )
        zmi     = MIN( rhog*qig/znin, zmimax )
        zmi     = MAX( zmi0, zmi )
        zsvmax  = (qvg - zqvsi) * zdtr
        zsagg   = zcagg(iv) * EXP(ccsaxp*LOG(zcslam(iv))) * qig
        zsagg   = MAX( zsagg, 0.0_ireals ) & !* zrho1o2(iv) &
          * MAX(0.2_ireals,MIN(EXP(0.09_ireals*(tg-t0)),1.0_ireals))
        znid      = rhog * qig/zmi
        IF (llqi) THEN
          zlnlogmi= LOG (zmi)
          zsidep    = zcidep(iv) * znid * EXP(0.33_ireals * zlnlogmi)   &
                        * ( qvg - zqvsi )
        ELSE
          zsidep = 0.0_ireals
        ENDIF
        zsvidep   = 0.0_ireals
        zsvisub   = 0.0_ireals
        zsimax    = qig*zdtr 
        IF( zsidep > 0.0_ireals ) THEN
          zsvidep = MIN( zsidep, zsvmax )
        ELSEIF (zsidep < 0.0_ireals ) THEN
          zsvisub = - MAX(-zsimax, zsvmax )
        ENDIF
        zsiau = zciau * MAX( qig - qi0, 0.0_ireals ) &
          * MAX(0.2_ireals,MIN(EXP(0.09_ireals*(tg-t0)),1.0_ireals))
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
        zssdep    = zcsdep(iv) * zxfac * ( qvg - zqvsi ) / (zcslam(iv)+zeps)**2

        ! Check for maximal depletion of vapor by sdep
        IF (zssdep > 0.0_ireals) zssdep = MIN(zssdep, zsvmax-zsvidep)
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
        zqvsw0      = fqvs( zpvsw0, ppg)
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
      ppg =  p(iv,k)

      zlnqrk      = LOG (zqrk(iv))
      zqvsw       = fqvs( fpvsw(tg), ppg )
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
        zqit = snuc(iv)  + scfrz(iv) - simelt(iv) - sicri(iv)  + sidep(iv) - sdau(iv)  - sagg(iv) &
          &  - siau(iv)
        zqrt = scau(iv)  + sshed(iv) + scac(iv)   + ssmelt(iv) - sev(iv)   - srcri(iv) - srfrz(iv) 
        zqst = siau(iv)  + sdau(iv)  + sagg(iv)   - ssmelt(iv) + sicri(iv) + srcri(iv) + srim(iv) &
								           + ssdep(iv) + srfrz(iv)
        ztt = cpdr*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst) )

        ! Update variables and add qi to qrs for water loading 
        qig = MAX ( 0.0_ireals, qig + zqit*zdt)
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

          IF (zprvr(iv) <= zqmin) zprvr(iv)=0.0_ireals
          IF (zprvs(iv) <= zqmin) zprvs(iv)=0.0_ireals

#ifdef NUDGING
          ! for the latent heat nudging
          IF ((llhn .OR. llhnverif) .AND. lhn_qrs) THEN
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)
            qrsflux(iv,k) = 0.5*(qrsflux(iv,k)+zpkr(iv)+zpks(iv))
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
        ELSE
          ! Precipitation fluxes at the ground
          prr_gsp(iv) = 0.5_ireals * (qrg*rhog*zvzr(iv) + zpkr(iv))
          prs_gsp(iv) = 0.5_ireals * (qsg*rhog*zvzs(iv) + zpks(iv))

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
!        qrs(iv,k   ) = qrg+qsg+qig
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
               rvd_m_o, lh_v, cpdr, cp_d,                 &
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
 CALL get_gs_lheating ('inc',1,ke)
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


!  CALL satad_v_3d (                             &
!                & maxiter  = 10_iintegers ,& !> IN
!                & tol      = 1.e-3_ireals ,& !> IN
!                & te       = t  (1,1,1)   ,&
!                & qve      = qv (1,1,1)   ,&
!                & qce      = qc (1,1,1)   ,&
!                & rhotot   = rho(1,1,1)   ,&
!                & idim     = nvec         ,&
!                & jdim     = 1            ,&
!                & kdim     = ke           ,&
!                & ilo      = iv_start     ,&
!                & iup      = iv_end       ,&
!                & jlo      = 1            ,&
!                & jup      = 1            ,&
!                & klo      = k_start      ,&
!                & kup      = ke            &
             !& count, errstat,
!                )
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

          ! restore input values
          t (iv,k) = t_in (iv,k)
          qv(iv,k) = qv_in(iv,k)
          qc(iv,k) = qc_in(iv,k)
          qi(iv,k) = qi_in(iv,k)
          qr(iv,k) = qr_in(iv,k)
          qs(iv,k) = qs_in(iv,k)

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


#ifdef __ICON__
SUBROUTINE SATAD ( kitera, te, qve, qce, tstart, phfe,                        &
  zdqd  , zqdwe, zh   , ztg0  , ztgn, zdqdt0, zgqd0, zphe ,  &
  b1, b2w, b3, b4w, b234w, rdrd, emrdrd, rddrm1, lh_v, cpdr, &
  cp_d, idim, jdim, ilo, iup, jlo, jup )

!-------------------------------------------------------------------------------
!
! Description:
!   This routine corrects the temperature (te), the specific humidity (qve)
!   and the cloud water content (qce) for condensation/evaporation.
!
! Method:
!   Saturation adjustment, i.e. reversible condensation/evaporation at
!   constant pressure by assuming chemical equilibrium of water and vapor
!
!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers), INTENT (IN)    ::  &
    kitera,              & !  Numver of iterations in the numerical scheme
    idim, jdim,          & !  Dimension of I/O-fields
    ilo, iup, jlo, jup     !  start- and end-indices for the computation

  REAL    (KIND=ireals),    INTENT (IN)    ::  &
    tstart  (idim,jdim), & ! Start temperature for iteration
    phfe    (idim,jdim)  ! Pressure (input)

  REAL    (KIND=ireals),    INTENT (INOUT) ::  &
    te      (idim,jdim), & ! Temperature on input/ouput
    qve     (idim,jdim), & ! Specific humidity on input/output
    qce     (idim,jdim), & ! Specific cloud water content on input/output
    zdqd    (idim,jdim), & !
    zqdwe   (idim,jdim), & !
    zh      (idim,jdim), & !
    ztg0    (idim,jdim), & !
    ztgn    (idim,jdim), & !
    zdqdt0  (idim,jdim), & !
    zgqd0   (idim,jdim), & !
    zphe    (idim,jdim)    !

  REAL    (KIND=ireals),    INTENT (IN)    ::  &
    b1, b2w, b3, b4w, b234w, rdrd, emrdrd, rddrm1, lh_v, cpdr, cp_d

! Local parameters: None
! ----------------
! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i, j,                & !  Loop indices
    nzit,                & !  Loop for iterations
    nsat,                & !  Number of saturated gridpoints
    iwrk(idim*jdim),     & !  i-index of saturated gridpoints
    jwrk(idim*jdim),     & !  j-index of saturated gridpoints
    indx                   !  loop index

  REAL    (KIND=ireals   ) ::  &
    zgeu  ,              & !
    zgqdu ,              & !
    zgew  ,              & !
    zqwmin,              & ! Minimum cloud water content for adjustment
    fgew  ,              & ! Name of satement function
    fgqd  ,              & ! ...
    fdqdt ,              & ! ...
    zt    ,              & ! Dummy argument for statement functions
    zge   ,              & ! ...
    zp    ,              & ! ...
    zgqd                   ! ...

  REAL    (KIND=ireals   ) ::  &
    minzdqd

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine satad
!-------------------------------------------------------------------------------

! STATEMENT FUNCTIONS

fgew(zt)       = b1*EXP( b2w*(zt-b3)/(zt-b4w) )
fgqd(zge,zp)   = rdrd*zge/( zp - emrdrd*zge )
fdqdt(zt,zgqd) = b234w*( 1.0_ireals + rddrm1*zgqd )*zgqd/( zt-b4w )**2

  zqwmin = 1.0E-20_ireals

  nsat = 0

  minzdqd= 1.0_ireals

  DO j = jlo , jup
    DO i = ilo , iup

      ! "save" the start values for the temperature
      ztg0 (i,j) = tstart(i,j)

      ! correction for negative values of qv and qc
      qve (i,j) = MAX( qve(i,j), 0.0_ireals )
      qce (i,j) = MAX( qce(i,j), 0.0_ireals )

      ! assume first subsaturation
      zqdwe(i,j)= qve(i,j) + qce(i,j)
      te (i,j)  = te(i,j) - lh_v*qce(i,j)*cpdr
      qve(i,j)  = zqdwe(i,j)
      qce(i,j)  = 0.0_ireals
      zgeu      = fgew(te(i,j))
      zgqdu     = fgqd(zgeu,phfe(i,j))
      zdqd(i,j) = zgqdu - zqdwe(i,j)
      minzdqd   = MIN(minzdqd,zdqd(i,j))

    ENDDO
  ENDDO

!NEC_CB if zdqd>=0, then for sure no points are found
  IF ( minzdqd >= 0.0_ireals ) RETURN

  DO j = jlo , jup
    DO i = ilo , iup

      IF (zdqd(i,j) < 0.0_ireals ) THEN
        nsat       = nsat+1
        iwrk(nsat) = i
        jwrk(nsat) = j
      ENDIF

    ENDDO
  ENDDO

  IF (nsat == 0) RETURN

! Do saturation adjustments for saturated gridpoints
! --------------------------------------------------

!cdir nodep
  DO indx = 1, nsat
     i = iwrk(indx)
     j = jwrk(indx)
     zh   (i,j) = cp_d*te(i,j) + lh_v*qve(i,j)
     zphe (i,j) = phfe(i,j)
     zgew       = fgew(ztg0(i,j))
     zgqd0(i,j) = fgqd(zgew,zphe(i,j))
  ENDDO

  IF ( kitera > 1 ) THEN
    DO  nzit  = 1 , kitera-1

!cdir nodep
      DO indx = 1, nsat
        i = iwrk(indx)
        j = jwrk(indx)
        zdqdt0(i,j) = fdqdt(ztg0(i,j),zgqd0(i,j))
        ztg0(i,j)   = (zh(i,j) - lh_v*(zgqd0(i,j)-zdqdt0(i,j)*ztg0(i,j)))/ &
                      ( cp_d + lh_v*zdqdt0(i,j) )
        zgew        = fgew(ztg0(i,j))
        zgqd0(i,j)  = fgqd(zgew,zphe(i,j))
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------

!cdir nodep
  DO indx = 1, nsat
      i = iwrk(indx)
      j = jwrk(indx)
      zdqdt0(i,j) = fdqdt(ztg0(i,j),zgqd0(i,j))
      ztgn(i,j)   = ( zh(i,j) - lh_v*(zgqd0(i,j)-zdqdt0(i,j)*ztg0(i,j)) ) / &
                    ( cp_d + lh_v*zdqdt0(i,j) )
      zgqd0(i,j)  = zgqd0(i,j) + zdqdt0(i,j)*( ztgn(i,j)-ztg0(i,j) )
  ENDDO

! Distribute the result on gridpoints
! -----------------------------------

!cdir nodep
  DO indx = 1, nsat
      i = iwrk(indx)
      j = jwrk(indx)
      te (i,j) =  ztgn(i,j)
      qve(i,j) = zgqd0(i,j)
      qce(i,j) = MAX( zqdwe(i,j) - zgqd0(i,j), zqwmin )
  ENDDO

! End of the subroutine

END SUBROUTINE satad
#endif

!==============================================================================

END MODULE gscp_hydci_pp
