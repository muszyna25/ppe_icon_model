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
!! implemented into ICOHAM by Kristina Froehlich and Axel Seifert (2010-06-10)
!! @par Revision History
!! warm-rain scheme (kessler_pp) implemented by Felix Rieper (2011-12-23)
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

MODULE mo_gscp_cosmo

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
    t0,           & !! melting temperature of ice
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

USE data_fields     , ONLY :   &
  tinc_lh    !! temperature increment due to latent heat      (  K  )

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

USE mo_kind              , ONLY: ireals=>wp     , &
                                 iintegers=>i4
USE mo_math_utilities    , ONLY: gamma_fct
USE mo_math_constants    , ONLY: pi
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                                 lh_v  => alv   , & !! latent heat of vapourization
                                 lh_s  => als   , & !! latent heat of sublimation
                                 cpdr  => rcpd  , & !! (spec. heat of dry air at constant press)^-1
                                 cvdr  => rcvd  , & !! (spec. heat of dry air at const vol)^-1
                                 t0    => tmelt     !! melting temperature of ice/snow
USE mo_satad             , ONLY: satad_v_3d     , & !! new saturation adjustment
                                 sat_pres_water , & !! saturation vapor pressure w.r.t. water
                                 sat_pres_ice!   , & !! saturation vapor pressure w.r.t. ice
!                                 spec_humi          !! Specific humidity 
USE mo_exception         , ONLY: message, message_text


#endif

!==============================================================================

IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: hydci_pp, hydci_pp_init, satad, kessler_pp

!------------------------------------------------------------------------------
!! Public variables
!------------------------------------------------------------------------------

PUBLIC :: cloud_num, mu_rain

!------------------------------------------------------------------------------
!> Parameters and variables which are global in this module
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers), PARAMETER ::  &
  iautocon       = 1,&       ! Autoconversion: 0 -> Kessler, 1 -> Seifert/Beheng
  isnow_n0temp   = 2

REAL    (KIND=ireals   ), PARAMETER ::  &
  zqmin = 1.0E-15_ireals,& !> threshold for computations
  zeps  = 1.0E-15_ireals,& !! small number

  !> Parameters for autoconversion of cloud water and cloud ice
  zccau  = 4.0E-4_ireals, & !! autoconversion coefficient (cloud water to rain)
  zciau  = 1.0E-3_ireals, & !! autoconversion coefficient (cloud ice   to snow)

  zkcau  = 9.44e+09_ireals, & !> kernel coeff for SB2001 autoconversion
  zkcac  = 5.25e+00_ireals, & !! kernel coeff for SB2001 accretion
  zcnue  = 2.00e+00_ireals, & !! gamma exponent for cloud distribution
  zxstar = 2.60e-10_ireals, & !! separating mass between cloud and rain
  zkphi1 = 6.00e+02_ireals, & !! constant in phi-function for autoconversion
  zkphi2 = 0.68e+00_ireals, & !! exponent in phi-function for autoconversion
  zkphi3 = 5.00e-05_ireals, & !! exponent in phi-function for accretion

  !> These coefficients have been precalculated on the basis of the
  !! following basic parameters (original hydci):
  !! N0R = 8.0 E6 : Parameter in the size distrubution function for rain
  !! ECR = 0.8    : Collection efficiency for rain collecting cloud water
  !! ECS = 0.9    : Collection efficiency for snow collecting cloud water
  !! EIR = 0.8    : Collection efficiency for rain collecting cloud ice
  !! V0R = 130.0  : Factor in the terminal velocity for raindrops
  !!                VTR(D) = V0R*D**(1/2)
  !! AMS = 0.069  : Formfactor in the mass-size relation of snowparticles
  !!                m(DS) = AMS*DS**2
  !! AMI = 130.0  : Formfactor in the mass-size relation of cloud ice crystals
  !!                m(DI) = AMI*DI**3
  !! ETA  = 1.75 E-5 : Viscosity of air
  !! DIFF = 2.22 E-5 : Molecular diffusion coefficient for water vapour
  !! LHEAT= 2.40 E-2 : Heat conductivity of air
  !! RHO  = 1.0      : Density of air
  !! RHOW = 1.0 E3   : Density of water

  zecs  = 0.9_ireals,   & !> Collection efficiency for snow collecting cloud water

  zadi  = 0.217_ireals, & !> Formfactor in the size-mass relation of ice particles
  zbdi  = 0.302_ireals, & !! Exponent in the size-mass relation of ice particles
  zams  = 0.069_ireals, & !! Formfactor in the mass-size relation of snow particles
  zbms  = 2.000_ireals, & !! Exponent in the mass-size relation of snow particles

  zv1s  = 0.50_ireals,  & !> Exponent in the terminal velocity for snow

  zami  = 130.0_ireals, & !> Formfactor in the mass-size relation of cloud ice
  zn0s0 = 8.0E5_ireals, & !!
  zn0s1 = 13.5_ireals * 5.65E5_ireals, & !! parameter in N0S(T)
  zn0s2 = -0.107_ireals , & !! parameter in N0S(T), Field et al
  zcac  = 1.72_ireals   , & !! (15/32)*(PI**0.5)*(ECR/RHOW)*V0R*AR**(1/8)
  zcicri= 1.72_ireals   , & !! (15/32)*(PI**0.5)*(EIR/RHOW)*V0R*AR**(1/8)
  zcrcri= 1.24E-3_ireals, & !! (PI/24)*EIR*V0R*Gamma(6.5)*AR**(-5/8)
  zcev  = 3.1E-3_ireals , & !! 2*PI*DIFF*HW*N0R*AR**(-1/2)
  zbev  = 9.0_ireals    , & !! 0.26*sqrt(0.5*RHO*v0r/eta)*Gamma(2.75)*AR**(-3/16)
  zcsmel= 1.48E-4_ireals, & !! 4*LHEAT*N0S*AS**(-2/3)/(RHO*lh_f)
  zbsmel= 14.37_ireals  , & !! 0.26*sqrt(0.5*RHO*v0s/eta)*Gamma(21/8)*AS**(-5/24)
  zasmel= 2.31E3_ireals , & !! DIFF*lh_v*RHO/LHEAT
  zcrfrz= 1.68_ireals   , & !! coefficient for raindrop freezing
  zrho0 = 1.225e+0_ireals, & !! reference air density
  zlheat = 2.40E-2_ireals, & !! thermal conductivity of dry air
  zeta   = 1.75e-5_ireals    !! kinematic viscosity of air

REAL    (KIND=ireals   ),  PARAMETER ::  &
  zthet  = 248.15_ireals  , & !> temperature for het. nuc. of cloud ice
  zthn   = 236.15_ireals  , & !! temperature for hom. freezing of cloud water
  ztrfrz = 271.15_ireals  , & !! threshold temperature for heterogeneous
  zmi0   = 1.0E-12_ireals, & !! initial crystal mass for cloud ice nucleation
  zmimax = 1.0E-9_ireals , & !! maximum mass of cloud ice crystals
  zmsmin = 3.0E-9_ireals , & !! initial mass of snow crystals
  zvz0r  = 12.63_ireals  , & !! coefficient of sedimentation velocity for rain

  !> Constant exponents in the transfer rate equations
  x1o12  =  1.0_ireals/12.0_ireals  ,  x3o16  =  3.0_ireals/16.0_ireals, &
  x7o8   =  7.0_ireals/ 8.0_ireals  ,  x2o3   =  2.0_ireals/ 3.0_ireals, &
  x5o24  =  5.0_ireals/24.0_ireals  ,  x1o8   =  1.0_ireals/ 8.0_ireals, &
  x13o8  = 13.0_ireals/ 8.0_ireals  ,  x13o12 = 13.0_ireals/12.0_ireals, &
  x27o16 = 27.0_ireals/16.0_ireals  ,  x1o3   =  1.0_ireals/ 3.0_ireals, &
  x1o2   = 0.5_ireals

!!==============================================================================

!> Global Variables for hydci_pp
! ----------------------

  REAL (KIND=ireals) ::           &
!!$    znimax,    & !> maximum number of cloud ice crystals
!!$    zpvsw0,    & !! sat.vap. pressure at melting temperature
    ccsrim,    & !!
    ccsagg,    & !!
    ccsdep,    & !!
    ccsvel,    & !!
    ccsvxp,    & !!
    ccslam,    & !!
    ccslxp,    & !!
    ccsaxp,    & !!
    ccsdxp,    & !!
    ccshi1,    & !!
    ccdvtp,    & !!
    ccidep,    & !!
    ccswxp,    & !!
    zconst!!$,    & !!
!!$    zxrmax,    & !! maximum mass of raindrops
!!$    zrbeva,    & !!
!!$    xpi6rhow, x1opi6rhow

!> Namelist Variables for hydci_pp and hydci_pp_gr
! --------------------------------------

  REAL (KIND=ireals) ::              &
    v0snow    = 20.0_ireals        , & !> factor in the terminal velocity for snow
    cloud_num = 200.00e+06_ireals  , & !! cloud droplet number concentration
    mu_rain   = 0.5_ireals             !! shape parameter of raindrop size distribution

CONTAINS

!==============================================================================

ELEMENTAL FUNCTION ice_nuclei_number(temp)
IMPLICIT NONE

REAL (KIND=ireals)              :: ice_nuclei_number
REAL (KIND=ireals), INTENT(IN)  :: temp

ice_nuclei_number = 1.0E2_ireals * EXP(0.2_ireals * (t0 - temp))

END FUNCTION ice_nuclei_number


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

  zconst = zkcau / (20.0_ireals*zxstar*cloud_num*cloud_num) &
         * (zcnue+2.0_ireals)*(zcnue+4.0_ireals)/(zcnue+1.0_ireals)**2
  ccsrim = 0.25_ireals*pi*zecs*v0snow*gamma_fct(zv1s+3.0_ireals)
  ccsagg = 0.25_ireals*pi*v0snow*gamma_fct(zv1s+3.0_ireals)
  ccsdep = 0.26_ireals*gamma_fct((zv1s+5.0_ireals)/2._ireals)*SQRT(0.5_ireals/zeta)
  ccsvxp = -(zv1s/(zbms+1.0_ireals)+1.0_ireals)
  ccsvel = zams*v0snow*gamma_fct(zbms+zv1s+1.0_ireals) &
         *(zams*gamma_fct(zbms+1.0_ireals))**ccsvxp
  ccsvxp = ccsvxp + 1.0_ireals
  ccslam = zams*gamma_fct(zbms+1.0_ireals)
  ccslxp = 1.0_ireals / (zbms+1.0_ireals)
  ccswxp = zv1s*ccslxp
  ccsaxp = -(zv1s+3.0_ireals)
  ccsdxp = -(zbms+1.0_ireals)/2.0_ireals
  ccshi1 = lh_s*lh_s/(zlheat*r_v)
  ccdvtp = 2.11E-5_ireals * t0**(-1.94_ireals) * 101325.0_ireals
  ccidep = 4.0_ireals * zami**(-x1o3)

  IF (PRESENT(idbg)) THEN
    IF (idbg > 12) THEN
      CALL message('mo_gscp_cosmo','hydci_pp_init: Initialized coefficients for hydci_pp')
      WRITE (message_text,'(A,E10.3)') '      ccslam = ',ccslam ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsvel = ',ccsvel ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsrim = ',ccsrim ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsagg = ',ccsagg ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccsdep = ',ccsdep ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccslxp = ',ccslxp ; CALL message('',message_text)
      WRITE (message_text,'(A,E10.3)') '      ccidep = ',ccidep ; CALL message('',message_text)
    ENDIF
  END IF

  CALL message('hydci_pp_init','microphysical values initialized')

END SUBROUTINE hydci_pp_init

!==============================================================================
!> Module procedure "hydci_pp" in "gscp" for computing effects of grid scale
!!  precipitation including cloud water, cloud ice, rain and snow in
!!  context with the Leapfrog and the Runge-Kutta time-integration
!------------------------------------------------------------------------------

SUBROUTINE hydci_pp(                 &
  ie,ke,                             & !> array dimensions
  istart,iend,kstart,                & !! optional start/end indicies
  idbg,                              & !! optional debug level
  l_cv,                              &
  dz,zdt,                            & !! numerics parameters
  t,p,rho,qv,qc,qi,qr,qs,            & !! prognostic variables
#ifdef __ICON__
  qi0,qc0,                           & !! cloud ice/water threshold for autoconversion
#endif
  prr_gsp,prs_gsp,                   & !! surface precipitation rates
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
    ie              ,    & !> number of grid points in first (zonal) direction
    ke                     !! number of grid points in vertical direction

  INTEGER, INTENT(IN), OPTIONAL ::  &
    istart    ,    & !> optional start index for computations in the parallel program
    iend      ,    & !! optional end index for computations in the parallel program
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
  REAL(KIND=ireals), DIMENSION(ie,ke), INTENT(INOUT) ::      &   ! (ie,ke)
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
  REAL(KIND=ireals), DIMENSION(ie,ke), INTENT(INOUT) ::   &   ! dim (ie,ke)
#else
  REAL(KIND=ireals), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
#endif
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    qr              ,    & !! specific rain content                         (kg/kg)
    qs                     !! specific snow content                         (kg/kg)

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  REAL(KIND=ireals), DIMENSION(ie), INTENT(INOUT) ::   &   ! dim (ie)
#else
  REAL(KIND=ireals), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
#endif
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    prs_gsp                !! precipitation rate of snow, grid-scale        (kg/(m2*s))

#ifdef __xlC__
  ! LL: xlc has trouble optimizing with the assumed shape, define the shape
  ! note: that these are actually intent(out)
  !       declared as intent(inout) to avoid copying
  REAL(KIND=ireals), DIMENSION(ie,ke), INTENT(INOUT), OPTIONAL ::   &     ! dim (ie,ke)
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
  REAL(KIND=ireals), DIMENSION(ie,ke), INTENT(INOUT), OPTIONAL ::   &   ! dim (ie,ke)
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


!! Local parameters: None, parameters are in data_gscp or data_constants
!! ----------------

!> Local scalars:
!! -------------

  INTEGER (KIND=iintegers) ::  &
    k, i             !> loop indees

  REAL    (KIND=ireals   ) :: nnr

  REAL    (KIND=ireals   ) ::  &
    znimax,            & !> maximum number of cloud ice crystals
    zpvsw0,            & !! sat.vap. pressure at melting temperature
    zqvsw0,            & !! sat.specific humidity at melting temperature
    zdtr ,             & !! reciprocal of timestep for integration
    zscau , zscac  , zscrim , zscshe, zsnuc , & !! local values of the
    zsiau , zsagg  , zsidep , zsicri, zsrcri, & !! transfer rates
    zsdau , zssdep , zssmelt, zsev,           & !! defined below (what about zsimelt?)
    zsrfrz,                                   & !!
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum, zsrmax,    & !! terms for limiting  total rain water depletion
    zssmax,            & !! term for limiting snow depletion
    znin,              & !! number of cloud ice crystals at nucleation
    znid,              & !! number of cloud ice crystals for deposition
    zmi ,              & !! mass of a cloud ice crystal
    zsvidep, zsvisub,  & !! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,   & !! terms for limiting total
    zqvsw,      & !! sat. specitic humidity at water saturation
    zztau, zxfac, zx1, zx2, ztt, &   !! some help variables
    ztau, zphi, zhi, zdvtp, ztc

  REAL    (KIND=ireals   ) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)

  INTEGER ::  &
    istartpar    ,    & !> start index for computations in the parallel program
    iendpar      ,    & !! end index for computations in the parallel program
    k_start      ,    & !! model level where computations start
    izdebug             !! debug level

  REAL(KIND=ireals), DIMENSION(ie,ke) ::   &
    t_in               ,    & !> temperature                                   (  K  )
    qv_in              ,    & !! specific water vapor content                  (kg/kg)
    qc_in              ,    & !! specific cloud water content                  (kg/kg)
    qi_in              ,    & !! specific cloud ice   content                  (kg/kg)
    qr_in              ,    & !! specific rain content                         (kg/kg)
    qs_in                     !! specific snow content                         (kg/kg)

  REAL    (KIND=ireals   ) ::  &
    zqct   ,& !> layer tendency of cloud water
    zqvt   ,& !! layer tendency of water vapour
    zqit   ,& !! layer tendency of cloud ice
    zqrt   ,& !! layer tendency of rain
    zqst      !! layer tendency of snow

  REAL (KIND=ireals)         ::       &
    mma(10), mmb(10)

  REAL    (KIND=ireals   ) ::  &
    zlnqrk,zlnqsk,     & !!
    zlnlogmi,                                                         & !!
    qcg,tg,qvg,qrg, qsg,qig,rhog,ppg, alf,bet,m2s,m3s,hlp

  LOGICAL :: &
    llqr,llqs,llqc,llqi  !!   switch for existence of qr, qs, qc, qi

!! Local (automatic) arrays:
!! -------------------------
#ifdef __COSMO__
  REAL    (KIND=ireals   ) ::  &
    zpres       (ie)        !! pressure
#endif

  REAL    (KIND=ireals   ) ::  &
    zqvsi       (ie),     & !> sat. specitic humidity at ice and water saturation
    zvzr        (ie),     & !!
    zvzs        (ie),     & !!
    zpkr        (ie),     & !!
    zpks        (ie),     & !!
    zprvr       (ie),     & !!
    zprvs       (ie),     & !!
    zdummy      (ie,8),   & !!
    zcsdep      (ie),     & !!
    zcidep      (ie),     & !!
    zvz0s       (ie),     & !!
    zcrim       (ie),     & !!
    zcagg       (ie),     & !!
    zbsdep      (ie),     & !!
    zcslam      (ie),     & !!
    zn0s        (ie),     & !!
    zimr        (ie),     & !!
    zims        (ie),     & !!
    zzar        (ie),     & !!
    zzas        (ie),     & !!
    zqrk        (ie),     & !!
    zqsk        (ie),     & !!
    zdtdh       (ie),     & !!
    z1orhog     (ie),     & !! 1/rhog
    zrho1o2     (ie),     & !! (rho0/rhog)**1/2
    zeln7o8qrk  (ie),     & !!
    zeln27o16qrk(ie),     & !!
    zeln13o8qrk (ie),     & !!
    zeln3o16qrk (ie),     & !!
    zeln13o12qsk(ie),     & !!
    zeln5o24qsk (ie),     & !!
    zeln2o3qsk  (ie)        !!

  REAL    (KIND=ireals   ) ::  &
    scau   (ie), & !> transfer rate due to autoconversion of cloud water
    scac   (ie), & !! transfer rate due to accretion of cloud water
    snuc   (ie), & !! transfer rate due nucleation of cloud ice
    scfrz  (ie), & !! transfer rate due homogeneous freezing of cloud water
    simelt (ie), & !! transfer rate due melting of cloud ice
    sidep  (ie), & !! transfer rate due depositional growth of cloud ice
    ssdep  (ie), & !! transfer rate due depositional growth of snow
    sdau   (ie), & !! transfer rate due depositional cloud ice autoconversion
    srim   (ie), & !! transfer rate due riming of snow
    sshed  (ie), & !! transfer rate due shedding
    sicri  (ie), & !! transfer rate due cloud ice collection by rain (sink qi)
    srcri  (ie), & !! transfer rate due cloud ice collection by rain (sink qr)
    sagg   (ie), & !! transfer rate due aggregation of snow and cloud ice
    siau   (ie), & !! transfer rate due autoconversion of cloud ice
    ssmelt (ie), & !! transfer rate due melting of snow
    sev    (ie), & !! transfer rate due evaporation of rain
    srfrz  (ie)    !! transfer rate due to rainwater freezing


  !> Integer arrays for a better vectorization
  INTEGER (KIND=iintegers) ::  &
    idx1(ie),   & !!
    idx2(ie),   & !!
    idx3(ie),   & !!
    idx4(ie),   & !!
    idx5(ie),   & !!
    idx6(ie)

  !> Dimensions and loop counter for storing the indices
  INTEGER (KIND=iintegers) ::  &
    ic1, ic2, ic3, ic4, ic5, ic6, i1d

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine hydci_pp
!------------------------------------------------------------------------------

!> Statement functions
! -------------------

!> Coeffs for moment relation based on 2nd moment (Field 2005)
  mma = (/5.065339_ireals, -0.062659_ireals, -3.032362_ireals, 0.029469_ireals, -0.000285_ireals, &
          0.312550_ireals,  0.000204_ireals,  0.003199_ireals, 0.000000_ireals, -0.015952_ireals /)
  mmb = (/0.476221_ireals, -0.015896_ireals,  0.165977_ireals, 0.007468_ireals, -0.000141_ireals, &
          0.060366_ireals,  0.000079_ireals,  0.000594_ireals, 0.000000_ireals, -0.003577_ireals /)

!> Some constant coefficients
  znimax = ice_nuclei_number(zthn) !! Maximum number of cloud ice crystals
  zpvsw0 = sat_pres_water(t0)      !! sat. vap. pressure for t = t0


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

!! Delete precipitation fluxes from previous timestep
  prr_gsp (:) = 0.0_ireals
  prs_gsp (:) = 0.0_ireals
  zpkr    (:) = 0.0_ireals
  zpks    (:) = 0.0_ireals
  zprvr   (:) = 0.0_ireals
  zprvs   (:) = 0.0_ireals
  zvzr    (:) = 0.0_ireals
  zvzs    (:) = 0.0_ireals
!CDIR COLLAPSE
  zdummy(:,:) = 0.0_ireals

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
  IF (PRESENT(istart)) THEN
    istartpar = istart
  ELSE
    istartpar = 1
  END IF
  IF (PRESENT(iend)) THEN
    iendpar = iend
  ELSE
    iendpar = ie
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

! output for various debug levels
  IF (izdebug > 25) CALL message('','SRC_GSCP: Start of hydci_pp')

  IF (izdebug > 25) THEN
    WRITE (message_text,'(A,E10.3)') '      ccslam = ',ccslam ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccsvel = ',ccsvel ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccsrim = ',ccsrim ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccsagg = ',ccsagg ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccsdep = ',ccsdep ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccslxp = ',ccslxp ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccidep = ',ccidep ; CALL message('',message_text)
  ENDIF
  IF (izdebug > 30) THEN
    WRITE (message_text,*) '   ie = ',ie ; CALL message('',message_text)
    WRITE (message_text,*) '   ke = ',ke ; CALL message('',message_text)
    WRITE (message_text,*) '   istartpar = ',istartpar ; CALL message('',message_text)
    WRITE (message_text,*) '   iendpar   = ',iendpar   ; CALL message('',message_text)
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

    zcrim (:) = 0.0_ireals
    zcagg (:) = 0.0_ireals
    zbsdep(:) = 0.0_ireals
    zvz0s (:) = 0.0_ireals
    zn0s  (:) = zn0s0

    ic1 = 0
    ic2 = 0

    DO i = istartpar, iendpar

      qrg = qr(i,k)
      qsg = qs(i,k)
      qvg = qv(i,k)
      qcg = qc(i,k)
      qig = qi(i,k)
      tg  = t(i,k)
      rhog = rho(i,k)

      !..for density correction of fall speeds
      z1orhog(i) = 1.0_ireals/rhog
      zrho1o2(i) = EXP(LOG(zrho0*z1orhog(i))*x1o2)

      zqrk(i) = qrg * rhog
      zqsk(i) = qsg * rhog

      llqr = zqrk(i) > zqmin
      llqs = zqsk(i) > zqmin

      zdtdh(i) = 0.5_ireals * zdt / dz(i,k)

      zzar(i)   = zqrk(i)/zdtdh(i) + zprvr(i) + zpkr(i)
      zzas(i)   = zqsk(i)/zdtdh(i) + zprvs(i) + zpks(i)

      IF (llqs) THEN
        ic1 = ic1 + 1
        idx1(ic1) = i
      ENDIF
      IF (llqr) THEN
        ic2 = ic2 + 1
        idx2(ic2) = i
      ENDIF
    ENDDO


!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs_prepare: DO i1d = 1, ic1
      i = idx1(i1d)

      qsg = qs(i,k)
      tg  = t(i,k)

      IF (isnow_n0temp == 1) THEN
        ! Calculate n0s using the temperature-dependent
        ! formula of Field et al. (2005)
        ztc = tg - t0
        ztc = MAX(MIN(ztc,0.0_ireals),-40.0_ireals)
        zn0s(i) = zn0s1*EXP(zn0s2*ztc)
        zn0s(i) = MIN(zn0s(i),1e9_ireals)
        zn0s(i) = MAX(zn0s(i),1e6_ireals)
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
        m2s = qsg * rho(i,k) / zams
        m3s = alf*EXP(bet*LOG(m2s))

        hlp  = zn0s1*EXP(zn0s2*ztc)
        zn0s(i) = 13.50_ireals * m2s**4 / m3s**3
        zn0s(i) = MAX(zn0s(i),0.5_ireals*hlp)
        zn0s(i) = MIN(zn0s(i),1e2_ireals*hlp)
        zn0s(i) = MIN(zn0s(i),1e9_ireals)
        zn0s(i) = MAX(zn0s(i),1e6_ireals)
      ELSE
        ! Old constant n0s
        zn0s(i) = 8.0e5_ireals
      ENDIF
      zcrim (i) = ccsrim*zn0s(i)
      zcagg (i) = ccsagg*zn0s(i)
      zbsdep(i) = ccsdep*SQRT(v0snow)
      zvz0s (i) = ccsvel*EXP(ccsvxp * LOG(zn0s(i)))

      IF (zvzs(i) == 0.0_ireals) THEN
        zvzs(i) = zvz0s(i) * EXP (ccswxp * LOG (0.5_ireals*zqsk(i))) * zrho1o2(i)
      ENDIF
    ENDDO loop_over_qs_prepare

!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qr_sedi: DO i1d = 1, ic2
      i = idx2(i1d)

      IF (zvzr(i) == 0.0_ireals) THEN
        zvzr(i) = zvz0r * EXP (x1o8  * LOG (0.5_ireals*zqrk(i))) * zrho1o2(i)
      ENDIF
    ENDDO loop_over_qr_sedi

  !----------------------------------------------------------------------------
  ! Section 3:
  !----------------------------------------------------------------------------

    zeln7o8qrk   (:) = 0.0_ireals
    zeln27o16qrk (:) = 0.0_ireals
    zeln13o8qrk  (:) = 0.0_ireals
    zeln3o16qrk  (:) = 0.0_ireals
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

    ic1 = 0
    ic2 = 0
    ic3 = 0
    ic4 = 0
    ic5 = 0
    ic6 = 0
    DO i = istartpar, iendpar

      qrg  = qr(i,k)
      qsg  = qs(i,k)
      qvg  = qv(i,k)
      qcg  = qc(i,k)
      qig  = qi(i,k)
      tg   =  t(i,k)
      rhog = rho(i,k)

      zqvsi(i) = sat_pres_ice(tg)/(rhog * r_v * tg)

      llqr = zqrk(i) > zqmin
      llqs = zqsk(i) > zqmin

      IF (llqr) THEN
        zpkr(i) = zqrk(i) * zvz0r * EXP (x1o8  * LOG (zqrk(i))) * zrho1o2(i)
      ELSE
        zpkr(i) = 0.0_ireals
      ENDIF

      IF (llqs) THEN
        zpks(i) = zqsk (i) * zvz0s(i) * EXP (ccswxp * LOG (zqsk(i))) * zrho1o2(i)
      ELSE
        zpks(i) = 0.0_ireals
      ENDIF

      zpkr(i)   = MIN( zpkr(i) , zzar(i) )
      zpks(i)   = MIN( zpks(i) , zzas(i) )

      zzar(i)   = zdtdh(i) * (zzar(i)-zpkr(i))
      zzas(i)   = zdtdh(i) * (zzas(i)-zpks(i))

      zimr(i)   = 1.0_ireals / (1.0_ireals + zvzr(i) * zdtdh(i))
      zims(i)   = 1.0_ireals / (1.0_ireals + zvzs(i) * zdtdh(i))

      zqrk(i)   = zzar(i)*zimr(i)
      zqsk(i)   = zzas(i)*zims(i)

      llqr = zqrk(i) > zqmin
      llqs = zqsk(i) > zqmin
      llqc =       qcg > zqmin
      llqi =       qig > zqmin

      IF (llqr) THEN
        ic1 = ic1 + 1
        idx1(ic1) = i
      ENDIF
      IF (llqs) THEN
        ic2 = ic2 + 1
        idx2(ic2) = i
      ENDIF
      IF (llqi .OR. llqs) THEN
        ic3 = ic3 + 1
        idx3(ic3) = i
      ENDIF
      IF ( tg < zthet .AND. qvg >  8.E-6_ireals &
                      .AND. qig <= 0.0_ireals ) THEN
        ic4 = ic4 + 1
        idx4(ic4) = i
      ENDIF
      IF (llqc) THEN
        ic5 = ic5 + 1
        idx5(ic5) = i
      ENDIF
      IF (llqr .AND. qcg <= 0.0_ireals) THEN
        ic6 = ic6 + 1
        idx6(ic6) = i
      ENDIF
    ENDDO


! ic1
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qr: DO i1d =1, ic1
      i = idx1(i1d)

      qcg  = qc(i,k)
      qig  = qi(i,k)
      tg   =  t(i,k)
      llqi =  qig > zqmin

      zlnqrk       = LOG (zqrk(i))
      IF ( qig+qcg > zqmin ) THEN
        zeln7o8qrk(i)   = EXP (x7o8   * zlnqrk)
      ENDIF
      IF ( tg < ztrfrz ) THEN
        zeln27o16qrk(i) = EXP (x27o16 * zlnqrk)
      ENDIF
      IF (llqi) THEN
        zeln13o8qrk(i)  = EXP (x13o8  * zlnqrk)
      ENDIF
      IF (qcg <= 0.0_ireals ) THEN
        zeln3o16qrk(i)  = EXP (x3o16  * zlnqrk)
      ENDIF
    ENDDO loop_over_qr

! ic2
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs_coeffs: DO i1d =1, ic2
      i = idx2(i1d)

      qcg = qc(i,k)
      qig = qi(i,k)

      zlnqsk       = LOG (zqsk(i))
      IF (qig+qcg > zqmin) THEN
        zeln13o12qsk(i) = EXP (x13o12 * zlnqsk)
      ENDIF
      zeln5o24qsk(i)  = EXP (x5o24  * zlnqsk)
      zeln2o3qsk(i)   = EXP (x2o3   * zlnqsk)
    ENDDO loop_over_qs_coeffs

! ic3
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qi_qs: DO i1d =1, ic3
      i = idx3(i1d)

      tg   =   t(i,k)
      ppg  =   p(i,k)
      rhog = rho(i,k)
      llqs = zqsk(i) > zqmin

      zdvtp  = ccdvtp * EXP(1.94_ireals * LOG(tg)) / ppg
      zhi    = ccshi1*zdvtp*rhog*zqvsi(i)/(tg*tg)
      hlp    = zdvtp / (1.0_ireals + zhi)
      zcidep(i) = ccidep * hlp
      IF (llqs) THEN
        zcslam(i) = EXP(ccslxp * LOG(ccslam * zn0s(i) / zqsk(i) ))
        zcslam(i) = MIN(zcslam(i),1e15_ireals)
        zcsdep(i) = 4.0_ireals * zn0s(i) * hlp
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
      i = idx4(i1d)

      qvg  =  qv(i,k)
      tg   =   t(i,k)
      rhog = rho(i,k)

      IF( qvg > zqvsi(i) ) THEN
        znin  = MIN( ice_nuclei_number(tg), znimax )
        zsnuc = zmi0 / rhog * znin * zdtr
        snuc(i) = zsnuc
      ENDIF
    ENDDO loop_over_icenucleation

    !--------------------------------------------------------------------------
    ! Section 5: Search for cloudy grid points with cloud water and
    !            calculation of the conversion rates involving qc
    !--------------------------------------------------------------------------

! ic5
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qc: DO i1d =1, ic5
      i = idx5(i1d)

      qrg  =   qr(i,k)
      qvg  =   qv(i,k)
      qcg  =   qc(i,k)
      qig  =   qi(i,k)
      tg   =    t(i,k)
      llqs = zqsk(i) > zqmin

      IF (iautocon == 0) THEN
        ! Kessler (1969) autoconversion rate
        zscau  = zccau * MAX( qcg - qc0, 0.0_ireals )
        zscac  = zcac  * qcg * zeln7o8qrk(i)
      ELSEIF (iautocon == 1) THEN
        ! Seifert and Beheng (2001) autoconversion rate
        ! with constant cloud droplet number concentration cloud_num
        IF (qcg > 1.0e-6_ireals) THEN
          ztau   = MIN(1.0_ireals-qcg/(qcg+qrg),0.9_ireals)
          zphi   = zkphi1 * ztau**zkphi2 * (1.0_ireals - ztau**zkphi2)**3
          zscau  = zconst * qcg*qcg*qcg*qcg &
            &    * (1.0_ireals + zphi/(1.0_ireals - ztau)**2)
          zphi   = (ztau/(ztau+zkphi3))**4
          zscac  = zkcac * qcg * qrg * zphi !* zrho1o2(i)
        ELSE
          zscau  = 0.0_ireals
          zscac  = 0.0_ireals
        ENDIF
      ENDIF
      IF (llqs) THEN
        zscrim = zcrim(i) * EXP(ccsaxp * LOG(zcslam(i))) * qcg !* zrho1o2(i)
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
        scfrz(i) = zscmax
      ELSE
        scau (i) = zcorr*zscau
        scac (i) = zcorr*zscac
        srim (i) = zcorr*zscrim
        sshed(i) = zcorr*zscshe
      ENDIF

      ! Calculation of heterogeneous nucleation of cloud ice.
      ! This is done in this section, because we require water saturation
      ! for this process (i.e. the existence of cloud water) to exist.
      ! Hetrogeneous nucleation is assumed to occur only when no
      ! cloud ice is present and the temperature is below a nucleation
      ! threshold.
      IF( tg <= 267.15_ireals .AND. qig <= 0.0_ireals ) THEN
        znin  = MIN( ice_nuclei_number(tg), znimax )
        zsnuc = zmi0 * z1orhog(i) * znin * zdtr
        snuc(i) = zsnuc
      ENDIF
      ! Calculation of in-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        zsrfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 )* zeln27o16qrk(i)
        srfrz(i) = zsrfrz
      ENDIF
    ENDDO loop_over_qc

      !------------------------------------------------------------------------
      ! Section 6: Search for cold grid points with cloud ice and/or snow and
      !            calculation of the conversion rates involving qi and ps
      !------------------------------------------------------------------------

! also ic3
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs: DO i1d =1, ic3
      i = idx3(i1d)

      qvg  =  qv(i,k)
      qig  =  qi(i,k)
      qsg  =  qs(i,k)
      tg   =   t(i,k)
      rhog = rho(i,k)
      llqi =  qig > zqmin

      IF (tg<=t0) THEN
        znin    = MIN( ice_nuclei_number(tg), znimax )
        zmi     = MIN( rhog*qig/znin, zmimax )
        zmi     = MAX( zmi0, zmi )
        zsvmax  = (qvg - zqvsi(i)) * zdtr
        zsagg   = zcagg(i) * EXP(ccsaxp*LOG(zcslam(i))) * qig
        zsagg   = MAX( zsagg, 0.0_ireals ) & !* zrho1o2(i) &
          * MAX(0.2_ireals,MIN(EXP(0.09_ireals*(tg-t0)),1.0_ireals))
        znid      = rhog * qig/zmi
        IF (llqi) THEN
          zlnlogmi= LOG (zmi)
          zsidep    = zcidep(i) * znid * EXP(0.33_ireals * zlnlogmi)   &
                        * ( qvg - zqvsi(i) )
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
        zsicri    = zcicri * qig * zeln7o8qrk(i)
        zsrcri    = zcrcri * (qig/zmi) * zeln13o8qrk(i)
        zxfac     = 1.0_ireals + zbsdep(i) * EXP(ccsdxp*LOG(zcslam(i)))
        zssdep    = zcsdep(i) * zxfac * ( qvg - zqvsi(i) ) / (zcslam(i)+zeps)**2

        ! Check for maximal depletion of vapor by sdep
        IF (zssdep > 0.0_ireals) zssdep = MIN(zssdep, zsvmax-zsvidep)
        ! Check for maximal depletion of snow by sdep
        IF (zssdep < 0.0_ireals) zssdep = MAX(zssdep, -qsg*zdtr)

        zsisum = zsiau + zsdau + zsagg + zsicri + zsvisub
        zcorr  = 0.0_ireals
        IF( zsimax > 0.0_ireals ) zcorr  = zsimax / MAX( zsimax, zsisum )
        sidep(i)  = zsvidep - zcorr*zsvisub
        sdau (i)  = zcorr*zsdau
        siau (i)  = zcorr*zsiau
        sagg (i)  = zcorr*zsagg
        ssdep(i)  = zssdep
        srcri(i)  = zsrcri
        sicri(i)  = zcorr*zsicri

      !------------------------------------------------------------------------
      ! Section 7: Search for warm grid points with cloud ice and/or snow and
      !            calculation of the melting rates of qi and ps
      !------------------------------------------------------------------------

      ELSE ! tg > 0
        simelt(i) = qig*zdtr
        zqvsw0      = zpvsw0 / (rhog * r_v * tg)
        zx1         = (tg - t0) + zasmel*(qvg - zqvsw0)
        zx2         = 1.0_ireals + zbsmel * zeln5o24qsk(i)
        zssmelt     = zcsmel * zx1 * zx2 * zeln2o3qsk(i)
        ssmelt(i) = MAX( zssmelt, 0.0_ireals )
      ENDIF ! tg
    ENDDO loop_over_qs

    !--------------------------------------------------------------------------
    ! Section 8: Search for grid points with rain in subsaturated areas
    !            and calculation of the evaporation rate of rain
    !--------------------------------------------------------------------------

! ic6
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qr_nocloud: DO i1d =1, ic6
      i = idx6(i1d)

      qvg = qv(i,k)
      tg  =  t(i,k)
      rhog = rho(i,k)

      zqvsw    = sat_pres_water(tg)/(rhog * r_v *tg)
      zx1      = 1.0_ireals + zbev* zeln3o16qrk(i)
      zsev     = zcev*zx1*(zqvsw - qvg)*SQRT(zqrk(i))
      sev(i) = MAX( zsev, 0.0_ireals )
      ! Calculation of below-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        zsrfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 ) * zeln27o16qrk(i)
        srfrz(i)  = zsrfrz
      ENDIF
    ENDDO loop_over_qr_nocloud

    !--------------------------------------------------------------------------
    ! Section 9: Calculate the total tendencies of the prognostic variables.
    !            Update the prognostic variables in the interior domain.
    !--------------------------------------------------------------------------

    DO i = istartpar, iendpar

      qvg = qv(i,k)
      qcg = qc(i,k)
      qrg = qr(i,k)
      qsg = qs(i,k)
      qig = qi(i,k)
      tg  = t (i,k)
      rhog = rho(i,k)

      zsrmax = zzar(i)*z1orhog(i)*zdtr
      zssmax = zzas(i)*z1orhog(i)*zdtr
      zsrsum = sev(i) + srfrz(i) + srcri(i)
      zcorr  = 1.0_ireals
      IF(zsrsum > 0.0_ireals) THEN
        zcorr  = zsrmax / MAX( zsrmax, zsrsum )
      ENDIF
      sev  (i) = zcorr*sev(i)
      srfrz(i) = zcorr*srfrz(i)
      srcri(i) = zcorr*srcri(i)

      ssmelt(i) = MIN(ssmelt(i), zssmax)
      IF (ssdep(i) < 0.0_ireals ) THEN
        ssdep(i) = MAX(ssdep(i), - zssmax)
      ENDIF
      zqvt = sev(i)   - sidep(i) - ssdep(i)  - snuc(i)
      zqct = simelt(i)- scau(i)  - scfrz(i)  - scac(i)   - sshed(i) - srim(i)
      zqit = snuc(i)  + scfrz(i) - simelt(i) - sicri(i)  + sidep(i) - sdau(i)  &
           & - sagg(i) - siau(i)
      zqrt = scau(i)  + sshed(i) + scac(i)   + ssmelt(i) - sev(i)   - srcri(i) &
           & - srfrz(i)
      zqst = siau(i)  + sdau(i)  + sagg(i)   - ssmelt(i) + sicri(i) + srcri(i) &
           & + srim(i) + ssdep(i) + srfrz(i)
      ztt = z_heat_cap_r*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst) )

      ! Update local variables
      qig = MAX ( 0.0_ireals, qig + zqit*zdt)
      qrg = MAX ( 0.0_ireals, (zzar(i)/rhog + zqrt*zdt)*zimr(i))
      qsg = MAX ( 0.0_ireals, (zzas(i)/rhog + zqst*zdt)*zims(i))
      qvg = MAX ( 0.0_ireals, qvg + zqvt*zdt )
      qcg = MAX ( 0.0_ireals, qcg + zqct*zdt )
      tg  = tg + ztt*zdt

      !----------------------------------------------------------------------
      ! Section 10: Complete time step
      !----------------------------------------------------------------------

      IF ( k /= ke) THEN
        ! Store precipitation fluxes and sedimentation velocities
        ! for the next level
        zprvr(i) = qrg*rhog*zvzr(i)
        zprvs(i) = qsg*rhog*zvzs(i)
        IF (zprvr(i) <= zqmin) zprvr(i)=0.0_ireals
        IF (zprvs(i) <= zqmin) zprvs(i)=0.0_ireals

        IF (qrg+qr(i,k+1) <= zqmin) THEN
          zvzr(i)= 0.0_ireals
        ELSE
          zvzr(i) = zvz0r                                               &
               &      * EXP(x1o8 *LOG((qrg+qr(i,k+1))*0.5_ireals*rhog)) &
               &      * zrho1o2(i)
        ENDIF
        IF (qsg+qs(i,k+1) <= zqmin) THEN
          zvzs(i)= 0.0_ireals
        ELSE
          zvzs(i) = zvz0s(i)                                                           &
               &      * EXP(zv1s/(zbms+1.0_ireals)*LOG((qsg+qs(i,k+1))*0.5_ireals*rhog)) &
               &      * zrho1o2(i)
        ENDIF
      ELSE
        ! Precipitation fluxes at the ground
        prr_gsp(i) = 0.5_ireals * (qrg*rhog*zvzr(i) + zpkr(i))
        prs_gsp(i) = 0.5_ireals * (qsg*rhog*zvzs(i) + zpks(i))

      ENDIF

      ! Update of prognostic variables or tendencies
      qr (i,k) = qrg
      qs (i,k) = qsg
      qi (i,k) = qig
      t  (i,k) = tg
      qv (i,k) = qvg
      qc (i,k) = qcg

    ENDDO

#ifndef __ICON__
    ! Store optional microphysical rates for diagnostics
    IF (PRESENT(ddt_diag_au   )) ddt_diag_au   (:,k) = scau  (:)
    IF (PRESENT(ddt_diag_ac   )) ddt_diag_ac   (:,k) = scac  (:)
    IF (PRESENT(ddt_diag_ev   )) ddt_diag_ev   (:,k) = sev   (:)
    IF (PRESENT(ddt_diag_nuc  )) ddt_diag_nuc  (:,k) = snuc  (:)
    IF (PRESENT(ddt_diag_idep )) ddt_diag_idep (:,k) = sidep (:)
    IF (PRESENT(ddt_diag_sdep )) ddt_diag_sdep (:,k) = ssdep (:)
    IF (PRESENT(ddt_diag_agg  )) ddt_diag_agg  (:,k) = sagg  (:)
    IF (PRESENT(ddt_diag_rim  )) ddt_diag_rim  (:,k) = srim  (:)
    IF (PRESENT(ddt_diag_rcri )) ddt_diag_rcri (:,k) = srcri (:)
    IF (PRESENT(ddt_diag_icri )) ddt_diag_icri (:,k) = sicri (:)
    IF (PRESENT(ddt_diag_dau  )) ddt_diag_dau  (:,k) = sdau  (:)
    IF (PRESENT(ddt_diag_iau  )) ddt_diag_iau  (:,k) = siau  (:)
    IF (PRESENT(ddt_diag_imelt)) ddt_diag_imelt(:,k) = simelt(:)
    IF (PRESENT(ddt_diag_smelt)) ddt_diag_smelt(:,k) = ssmelt(:)
    IF (PRESENT(ddt_diag_cfrz )) ddt_diag_cfrz (:,k) = scfrz (:)
    IF (PRESENT(ddt_diag_rfrz )) ddt_diag_rfrz (:,k) = srfrz (:)
    IF (PRESENT(ddt_diag_shed )) ddt_diag_shed (:,k) = sshed (:)
#endif

  IF (izdebug > 25) THEN
    ! Check for negative values
    DO i = istartpar, iendpar
      IF (qr(i,k) < 0.0_ireals) THEN
        WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qr'
        CALL message('',message_text)
      ENDIF
      IF (qc(i,k) < 0.0_ireals) THEN
        WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qc'
        CALL message('',message_text)
      ENDIF
      IF (qi(i,k) < 0.0_ireals) THEN
        WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qi'
        CALL message('',message_text)
      ENDIF
      IF (qs(i,k) < 0.0_ireals) THEN
        WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qs'
        CALL message('',message_text)
      ENDIF
      IF (qv(i,k) < 0.0_ireals) THEN
        WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qv'
        CALL message('',message_text)
      ENDIF
    ENDDO
  ENDIF

#if defined (__COSMO__)
  ! Do a final saturation adjustment for new values of t, qv and qc
!CDIR COLLAPSE
    zpres(:) = p(:,k)

    CALL satad ( 1, t(:,k), qv(:,k),                   &
               qc(:,k), t(:,k), zpres,                 &
               zdummy(:,1),zdummy(:,2),zdummy(:,3),    &
               zdummy(:,4),zdummy(:,5),zdummy(:,6),    &
               zdummy(:,7),zdummy(:,8),                &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,  &
               rvd_m_o, lh_v, cpdr, cp_d,              &
               ie, istartpar, iendpar )

  IF ( ldiabf_lh ) THEN
    ! compute temperature increment due to latent heat
!CDIR COLLAPSE
    tinc_lh(:,k) = tinc_lh(:,k) + t(:,k)
  ENDIF
#endif

ENDDO loop_over_levels

#ifdef __ICON__

 CALL satad_v_3d (                             &
               & maxiter  = 10_iintegers ,& !> IN
               & tol      = 1.e-3_ireals ,& !> IN
               & te       = t            ,&
               & qve      = qv           ,&
               & qce      = qc           ,&
               & rhotot   = rho          ,&
               & idim     = ie           ,&
               & kdim     = ke           ,&
               & ilo      = istartpar    ,&
               & iup      = iendpar      ,&
               & klo      = k_start      ,&
               & kup      = ke            &
!              !& count, errstat,
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
      DO i=istartpar,iendpar

        ! calculated pseudo-tendencies
        ddt_tend_t (i,k) = (t (i,k) - t_in (i,k))*zdtr
        ddt_tend_qv(i,k) = MAX(-qv_in(i,k)*zdtr,(qv(i,k) - qv_in(i,k))*zdtr)
        ddt_tend_qc(i,k) = MAX(-qc_in(i,k)*zdtr,(qc(i,k) - qc_in(i,k))*zdtr)
        ddt_tend_qr(i,k) = MAX(-qr_in(i,k)*zdtr,(qr(i,k) - qr_in(i,k))*zdtr)
        ddt_tend_qs(i,k) = MAX(-qs_in(i,k)*zdtr,(qs(i,k) - qs_in(i,k))*zdtr)
        ddt_tend_qi(i,k) = MAX(-qi_in(i,k)*zdtr,(qi(i,k) - qi_in(i,k))*zdtr)

        ! restore input values
        t (i,k) = t_in (i,k)
        qv(i,k) = qv_in(i,k)
        qc(i,k) = qc_in(i,k)
        qi(i,k) = qi_in(i,k)
        qr(i,k) = qr_in(i,k)
        qs(i,k) = qs_in(i,k)

      END DO
    END DO

  END IF

  IF (izdebug > 25) THEN
    CALL message('mo_gscp', 'UPDATED VARIABLES')
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





SUBROUTINE kessler_pp(               &
  ie,ke,                             & ! array dimensions
  istart,iend,kstart,                & ! optional start/end indicies
  idbg,                              & ! optional debug level
  l_cv,                              & ! switch for cv, cp 
  dz,zdt,                            & ! numerics parameters
  t,p,rho,qv,qc,qr,                  & ! prognostic variables
#ifdef __ICON__
  qc0,                               & ! cloud water threshold for autoconversion
#endif
  prr_gsp,                           & ! surface precipitation rates
  ddt_tend_t     , ddt_tend_qv     , &
  ddt_tend_qc    ,                   & ! ddt_tend_xx are tendencies
  ddt_tend_qr    ,                   & !    necessary for dynamics
  ddt_diag_au    , ddt_diag_ac     , & !
  ddt_diag_ev                        ) ! ddt_diag_xxx are optional
                                       !   diagnostic tendencies of all
                                       !   individual microphysical processes

!------------------------------------------------------------------------------
!>
!! Description:
!!   This module procedure calculates the rates of change of temperature, 
!!   cloud water, water vapor, and rain due to cloud microphysical processes 
!!   related to the formation of grid scale precipitation. The variables will
!!   be updated in this subroutine.
!!   The precipitation fluxes at the 
!!   surface are stored on the corresponding global fields.
!!   This subroutine relies on conversion rates used in the subroutine 
!!   hydorprog.
!!
!! Method:
!!   The tendencies involving cloud water are calculated using a full implicit 
!!   scheme whereas evaporation is computed explicitly.
!!   Rain is a prognostic variable and the sedimentation term is 
!!   computed implicitly.
!
!------------------------------------------------------------------------------
!
! Declarations:

! Subroutine arguments: 
! --------------------

  INTEGER, INTENT(IN) ::  &
    ie               ,    & !< number of grid points in first (zonal) direction
    ke                      !< number of grid points in vertical direction

  INTEGER, INTENT(IN), OPTIONAL ::  &
    istart    ,    & !< optional start index for computations in the parallel program
    iend      ,    & !< optional end index for computations in the parallel program
    kstart    ,    & !< optional start index for the vertical index
    idbg             !< optional debug level

  REAL(KIND=ireals), INTENT(IN) :: &
    zdt              !< time step for integration of microphysics     (  s  )
#ifdef __ICON__
  REAL(KIND=ireals), INTENT(IN) :: &
    qc0              !< cloud ice/water threshold for autoconversion
#endif

  REAL(KIND=ireals), DIMENSION(ie,ke), INTENT(IN) ::      &
    dz              ,    & !< layer thickness of full levels                (  m  )
    rho             ,    & !< density of moist air                          (kg/m3)
    p                      !< pressure                                      ( Pa  )

  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv                   !< if true, cv is used instead of cp

  REAL(KIND=ireals), DIMENSION(ie,ke), INTENT(INOUT) ::   &
    t               ,    & !< temperature                                   (  K  )
    qv              ,    & !< specific water vapor content                  (kg/kg)
    qc              ,    & !< specific cloud water content                  (kg/kg)
    qr                     !< specific rain content                         (kg/kg)

  REAL(KIND=ireals), DIMENSION(ie), INTENT(INOUT) ::   &
    prr_gsp                !<  precipitation rate of rain, grid-scale        (kg/(m2*s))

  REAL(KIND=ireals), DIMENSION(ie,ke), INTENT(OUT), OPTIONAL ::   &
    ddt_tend_t      , & !< tendency T                                       ( 1/s )
    ddt_tend_qv     , & !< tendency qv                                      ( 1/s )
    ddt_tend_qc     , & !< tendency qc                                      ( 1/s )
    ddt_tend_qr         !< tendency qr                                      ( 1/s )

  REAL(KIND=ireals), DIMENSION(ie,ke), INTENT(OUT), OPTIONAL ::   &
    ddt_diag_au     , & !< optional output autoconversion rate cloud to rain           ( 1/s )     
    ddt_diag_ac     , & !< optional output accretion rate cloud to rain                ( 1/s )
    ddt_diag_ev         !< optional output evaporation of rain                         ( 1/s )


! Local parameters:
! ----------------
! TODO[FR] : put parameters to module header

  REAL    (KIND=ireals   ), PARAMETER ::  &
    ! basic constants of the parameterization scheme

    zaau   = 1.0_ireals/1000.0_ireals,         & ! coef. for autoconversion  
    zaac   = 1.72_ireals,                      & ! coef. for accretion (neu)
!    zbev   = 9.1_ireals,                      & ! coef. for drop ventilation,
    !                                              already module var

    ! constants for the process rates
    z1d8   = 1.0_ireals/8.0_ireals,               & !
    z7d8   = 7.0_ireals/8.0_ireals,               & !
    z3d16  = 3.0_ireals/16.0_ireals,              & !

    ! constants for sedimentation
!    zvz0r = 12.63_ireals,                        & ! coefficient of sedimentation velocity
    !                                                 for rain, alread module variable
 
    ! to avoid illegal operations in power expressions
    znull = 1.E-20_ireals

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    k, i                 ! loop indices
  
  

  REAL    (KIND=ireals   ) ::  &
    ztx   ,            & ! 
!xxx    zpx   ,            & ! dummy arguments for statement function
!xxx    zgex  ,            & ! -""-
!xxx    fspw  ,            & ! function for equilibrium vapour pressure over water
!xxx    fsqv  ,            & ! function for specific humidity at 
    fsa3  ,            & !
!xxx    zspw  ,            & ! equilibrium vapour pressure over water
    zsa3  ,            & !
    zc3   ,            & !
    zc1c  ,            & !
    zc1   ,            & !
    zx    ,            & !
    zsrmax

  REAL    (KIND=ireals   ) ::  & 
!    zphf  ,            & !  pressure in a k-layer
    zsqvw ,            & !  specific humidity at water saturation
    zqvts ,            & !  qv-tendency in a layer
    zqcts ,            & !  qc-tendency in a layer
    zqrts ,            & !  qr-tendency in a layer
    ztts  ,            & !  t -tendency in a layer
    ztau  ,            & !  some auxiliary variables
    zphi  
  

  REAL    (KIND=ireals   ) :: z_heat_cap_r ! reciprocal of cpdr or cvdr (depending on l_cv)
  
  INTEGER ::  &
    istartpar    ,    & ! start index for computations in the parallel program
    iendpar      ,    & ! end index for computations in the parallel program
    k_start      ,    & ! model level where computations start
    izdebug             ! debug level

  REAL(KIND=ireals), DIMENSION(ie,ke) ::   &
    t_in         ,    & ! temperature                                   (  K  )
    qv_in        ,    & ! specific water vapor content                  (kg/kg)
    qc_in        ,    & ! specific cloud water content                  (kg/kg)
    qr_in                     ! specific rain content                         (kg/kg)

  REAL    (KIND=ireals   ) ::  & 
    zdtr         ,     & ! 
    zimr         ,     & !
    qvg          ,     & !  specific water vapor content: local grid cell value
    qcg          ,     & !  specfic cloud water content:  - "" -
    tg           ,     & !  temperature:                  - "" -
    qrg          ,     & !  specific rain content:        - "" -
    rhog         ,     & !  density                       - "" -
    rhogr        ,     & !  reciprocal density            - "" -
    ppg                  !  full level pressure           - "" -
  

  LOGICAL :: &
    llqr, llqc           !  switch for existence of qr, qc

 
  !! Local (automatic) arrays:
  !! -------------------------
#ifdef __COSMO__
  REAL    (KIND=ireals   ) :: &
    zpres       (ie)     !  pressure
#endif
  
  REAL    (KIND=ireals   ) ::  &
    lnzqrk (ie),       & !
    zzar   (ie),       & !
    zdtdh  (ie),       & !
    zvzr   (ie),       & !
    zprvr  (ie),       & !
    zpkr   (ie),       & !
    zqrk   (ie),       & !
    zpkm1r (ie),       & !
    zdummy (ie,8),     & !
    zsrd   (ie),       & !  local evaporation rate S_ev
    zswra  (ie),       & !  local autoconversion rate S_au
    zswrk  (ie)          !  local accretion rate S_ac

    
  ! Debugging / auxiliary variables
  INTEGER, PARAMETER ::  iautocon1 = 0     ! 0 -> Kessler 1969, 1-> Seifert, Beheng 2001
  


    
!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
! Begin Subroutine kessler_pp
!------------------------------------------------------------------------------

! Statement functions
! -------------------
 
!xxx fspw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) ) ! sat.vap. pressure over water
!xxx fsqv(zgex,zpx) = rdv*zgex/( zpx - o_m_rdv*zgex ) ! specific humidity at sat.
  
! coefficients for conversion rates
fsa3(ztx) = 3.86E-3_ireals - 9.41E-5_ireals*(ztx-t0) 

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local variables
!
!------------------------------------------------------------------------------

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

  ! Delete precipitation fluxes from previous timestep
  prr_gsp (:) = 0.0_ireals
  zpkr    (:) = 0.0_ireals
  zprvr   (:) = 0.0_ireals
  zvzr    (:) = 0.0_ireals
  ! CDIR COLLAPSE
  zdummy(:,:)=0.0_ireals
  
  ! Optional arguments
  
  IF (PRESENT(ddt_tend_t)) THEN
    ! save input arrays for final tendency calculation
    t_in(:,:)  = t(:,:)
    qv_in(:,:) = qv(:,:)
    qc_in(:,:) = qc(:,:)
    qr_in(:,:) = qr(:,:)
  END IF

  IF (PRESENT(istart)) THEN
    istartpar = istart
  ELSE
    istartpar = 1
  END IF

  IF (PRESENT(iend)) THEN
    iendpar = iend
  ELSE
    iendpar = ie
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
  
  !------------------------------------------------------------------------------
  !>  Initial setting of local and global variables
  !------------------------------------------------------------------------------
  
  zconst = zkcau / (20.0_ireals*zxstar*cloud_num*cloud_num) &
    * (zcnue+2.0_ireals)*(zcnue+4.0_ireals)/(zcnue+1.0_ireals)**2
  
  ! output for various debug levels
  IF (izdebug > 25) CALL message('','SRC_GSCP: Start of keppler_pp')

!!$  IF (izdebug > 15) THEN
!!$    WRITE (message_text,'(A,E10.3)') '      ccslam = ',ccslam ; CALL message('',message_text)
!!$    WRITE (message_text,'(A,E10.3)') '      ccsvel = ',ccsvel ; CALL message('',message_text)
!!$    WRITE (message_text,'(A,E10.3)') '      ccsrim = ',ccsrim ; CALL message('',message_text)
!!$    WRITE (message_text,'(A,E10.3)') '      ccsagg = ',ccsagg ; CALL message('',message_text)
!!$    WRITE (message_text,'(A,E10.3)') '      ccsdep = ',ccsdep ; CALL message('',message_text)
!!$    WRITE (message_text,'(A,E10.3)') '      ccslxp = ',ccslxp ; CALL message('',message_text)
!!$    WRITE (message_text,'(A,E10.3)') '      ccidep = ',ccidep ; CALL message('',message_text)
!!$  ENDIF

  IF (izdebug > 30) THEN
    WRITE (message_text,*) '   ie = ',ie ; CALL message('',message_text)
    WRITE (message_text,*) '   ke = ',ke ; CALL message('',message_text)
    WRITE (message_text,*) '   istartpar = ',istartpar ; CALL message('',message_text)
    WRITE (message_text,*) '   iendpar   = ',iendpar   ; CALL message('',message_text)
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
  ENDIF


!!$  ! select timelevel and timestep for calculations
!!$  nu    = nnew
!!$  IF ( l2tls ) THEN
!!$    zdt   = dt
!!$  ELSE
!!$    zdt   = dt2
!!$  ENDIF


! ----------------------------------------------------------------------
! Loop from the top of the model domain to the surface to calculate the
! transfer rates  and sedimentation terms
! ----------------------------------------------------------------------


!!$#ifdef NUDGING
!!$  ! add part of latent heating calculated in subroutine kessler_pp to model latent
!!$  ! heating field: subtract temperature from model latent heating field
!!$    IF (llhn) &
!!$       CALL get_gs_lheating ('add',1,ke)
!!$#endif

  loop_over_levels: DO k = k_start, ke

#ifdef __COSMO__
    IF ( ldiabf_lh ) THEN
      ! initialize temperature increment due to latent heat
      tinc_lh(:,k) = tinc_lh(:,k) - t(:,k)
    ENDIF
#endif
    
    
  !----------------------------------------------------------------------------
  !  Section 2: Test for clouds and precipitation in present layer.
  !             If no cloud or precipitation points have been found
  !             go to the next layer.
  !----------------------------------------------------------------------------

    DO i = istartpar, iendpar 
      IF(qr (i,k) < 1.0e-15_ireals) qr (i,k) = 0.0_ireals
      IF(qc (i,k) < 1.0e-15_ireals) qc (i,k) = 0.0_ireals
    END DO
    
    loop_over_i: DO i = istartpar, iendpar 
      
      qcg  = qc(i,k)
      qvg  = qv(i,k)
      qrg  = qr(i,k)
      tg   = t(i,k)
      ppg  = p(i,k)
      rhog = rho(i,k)
      rhogr= 1.0_ireals / rhog
!      zphf = p0(i,k) + pp(i,k)  !xxx

      zqrk(i) = qrg * rhog

      ! check existence of rain and cloud
      llqr = zqrk(i) > zqmin    ! zqmin: module var
      llqc = qcg     > zqmin


      zdtdh(i)  = 0.5_ireals * zdt / dz(i,k)

      zpkm1r(i) = zpkr(i)  
      zzar(i)   = zqrk(i)/zdtdh(i) + zprvr(i) + zpkm1r(i) 

      IF (llqr) THEN
        zpkr(i)   = zqrk(i) * zvz0r * EXP(z1d8*LOG(zqrk(i)))
      ELSE
        zpkr(i)   = 0.0_ireals
      ENDIF  
      zpkr(i)   = MIN( zpkr(i), zzar(i) )
      zzar(i)   = zdtdh(i) * (zzar(i)-zpkr(i))

      IF( zvzr(i) == 0.0_ireals ) THEN 
        ! xxx: why different from hydci_pp_ICON?
        zvzr(i) = zvz0r * EXP(z1d8 * LOG(MAX(0.5_ireals*zqrk(i),znull)))
      END IF
 
      zimr        = 1.0_ireals / (1.0_ireals + zvzr(i) * zdtdh(i))    ! implicit factor
      zqrk(i)     = zzar(i)*zimr

!      lnzqrk(i)      = LOG(MAX(zqrk(i),znull))  ! Optimization
      IF (llqr) THEN
        lnzqrk(i)      = LOG(zqrk(i))
      ELSE
        lnzqrk(i)      = 0.0_ireals
      ENDIF
        
      ztts        = 0.0_ireals
      zqvts       = 0.0_ireals
      zqcts       = 0.0_ireals
      
      !------------------------------------------------------------------------
      !  Section 5: Calculation of cloud microphysics for cloud case
      !             ( qc > 0)
      !------------------------------------------------------------------------

      IF (llqc) THEN
! Calculate conversion rates
        
! Coefficients
        
        IF (iautocon1 == 0) THEN
! Kessler (1969) autoconversion rate 
          zc1c  = zaac *  EXP(lnzqrk(i)*z7d8)
          zc1   = zaau + zc1c
          zx    = qcg / (1.0_ireals + zc1*zdt)
          
! Conversion rates
          zswra(i)  = zaau * zx    ! S_au
          zswrk(i)  = zc1c * zx    ! S_ac

        ELSE
! Seifert and Beheng (2001) autoconversion rate
! with constant cloud droplet number concentration cloud_num

          IF (qcg > 1.0e-6_ireals) THEN
            ztau      = MIN(1.0_ireals-qcg/(qcg+qrg),0.9_ireals)
            zphi      = zkphi1 * ztau**zkphi2 * (1.0_ireals - ztau**zkphi2)**3
            zswra(i)  = zconst * qcg*qcg*qcg*qcg &                      ! S_au
              &    * (1.0_ireals + zphi/(1.0_ireals - ztau)**2)
            zphi      = (ztau/(ztau+zkphi3))**4
            zswrk(i)  = zkcac * qcg * qrg * zphi !* zrho1o2(i)          ! S_ac
          ELSE
            zswra(i)  = 0.0_ireals  ! S_au                    
            zswrk(i)  = 0.0_ireals  ! S_ac
          ENDIF
          
        END IF
 
        ! Store tendencies
        zqcts  = - zswra(i) - zswrk(i)
        zqrts  =   zswra(i) + zswrk(i)

!>FR: test Kessler scheme with losing water load
!        zqrts = 0.0_ireals
!<FR
 
        ! Update values

        qrg = MAX(0.0_ireals,(zzar(i)*rhogr + zqrts*zdt)*zimr)

      !------------------------------------------------------------------------
      !  Section 7: Calculation of cloud microphysics for
      !             precipitation case without cloud ( qc = 0 )
      !             -> Evaporation
      !------------------------------------------------------------------------

      ELSEIF ( llqr .AND. .NOT.llqc )    THEN

!xxx        zspw  = fspw(tg)
!xxx        zspw  = sat_pres_water(tg)
!xxx        zsqvw = fsqv(zspw, ppg)   
        zsqvw = sat_pres_water(tg)/(rhog * r_v *tg)

        ! Coefficients
        zsrmax = zzar(i)*rhogr * zdtr
        zsa3   = fsa3(tg)
        zc3    = zsa3 * SQRT(zqrk(i)) * (1.0_ireals + zbev * EXP(lnzqrk(i)*z3d16))
 
        ! Conversion rate: S_ev
        zsrd(i)   = -zc3 * (qvg - zsqvw)
        zsrd(i)   = MIN (zsrmax, zsrd(i))         

 
        ! Store tendencies
        zqvts =   zsrd(i)
        ztts  = - lh_v*zsrd(i)*z_heat_cap_r
        zqrts = - zsrd(i)

        ! Update values

        qrg = MAX(0.0_ireals,(zzar(i)*rhogr + zqrts*zdt)*zimr)
 
      ENDIF         

      !------------------------------------------------------------------------
      !  Section 8: Complete time step
      !------------------------------------------------------------------------

      IF ( k /= ke ) THEN
        ! Store precipitation fluxes and sedimentation velocities for the next level
        zprvr(i) = qrg*rhog*zvzr(i)
        !zvzr(i)  = zvz0r * EXP(z1d8 * LOG(MAX((qrg+qr(i,k+1,nu))*0.5_ireals*rhog,znull)))
        zvzr(i)  = zvz0r * EXP(z1d8 * LOG(MAX((qrg+qr(i,k+1))*0.5_ireals*rhog,znull)))
      ELSE
        ! Precipitation flux at the ground
        prr_gsp(i) = 0.5_ireals * (qrg*rhog*zvzr(i) + zpkr(i))
      ENDIF

      ! Update of prognostic variables or tendencies
      qr (i,k) = qrg
      t  (i,k) = t (i,k) + ztts*zdt 
      qv (i,k) = MAX ( 0.0_ireals, qv(i,k) + zqvts*zdt )
      qc (i,k) = MAX ( 0.0_ireals, qc(i,k) + zqcts*zdt )

    ENDDO loop_over_i

#ifndef __ICON__
    ! Store optional microphysical rates for diagnostics
    IF (PRESENT(ddt_diag_au   )) ddt_diag_au(:,k) = zswra(:)
    IF (PRESENT(ddt_diag_ac   )) ddt_diag_ac(:,k) = zswrk(:)
    IF (PRESENT(ddt_diag_ev   )) ddt_diag_ev(:,k) = zsrd(:)
#endif

    IF (izdebug > 25) THEN
      ! Check for negative values
      DO i = istartpar, iendpar
        IF (qr(i,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: kessler_pp, negative value in qr'
          CALL message('',message_text)
        ENDIF
        IF (qc(i,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: kessler_pp, negative value in qc'
          CALL message('',message_text)
        ENDIF
        IF (qv(i,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: kessler_pp, negative value in qv'
          CALL message('',message_text)
        ENDIF
      ENDDO
    ENDIF
    
#if defined (__COSMO__)
    !Do a final saturation adjustment for new values of t, qv and qc
    zpres(:) = p(:,k)

    CALL satad ( 1, t(:,k), qv(:,k),          &
      qc(:,k), t(:,k), zpres,                 &
      zdummy(:,1),zdummy(:,2),zdummy(:,3),    &
      zdummy(:,4),zdummy(:,5),zdummy(:,6),    &
      zdummy(:,7),zdummy(:,8),                &
      b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,  &
      rvd_m_o, lh_v, cpdr, cp_d,              &
      ie, istartpar, iendpar )

    IF ( ldiabf_lh ) THEN
      ! compute temperature increment due to latent heat
      tinc_lh(:,k) = tinc_lh(:,k) + t(:,k)
    ENDIF
#endif
    
ENDDO loop_over_levels


#ifdef __ICON__

 CALL satad_v_3d (                             &
               & maxiter  = 10_iintegers ,& !> IN
               & tol      = 1.e-3_ireals ,& !> IN
               & te       = t            ,&
               & qve      = qv           ,&
               & qce      = qc           ,&
               & rhotot   = rho          ,&
               & idim     = ie           ,&
               & kdim     = ke           ,&
               & ilo      = istartpar    ,&
               & iup      = iendpar      ,&
               & klo      = k_start      ,&
               & kup      = ke            &
!              !& count, errstat,
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
      DO i=istartpar,iendpar

        ! calculated pseudo-tendencies
        ddt_tend_t (i,k) = (t (i,k) - t_in (i,k))*zdtr
        ddt_tend_qv(i,k) = MAX(-qv_in(i,k)*zdtr,(qv(i,k) - qv_in(i,k))*zdtr)
        ddt_tend_qc(i,k) = MAX(-qc_in(i,k)*zdtr,(qc(i,k) - qc_in(i,k))*zdtr)
        ddt_tend_qr(i,k) = MAX(-qr_in(i,k)*zdtr,(qr(i,k) - qr_in(i,k))*zdtr)

        ! restore input values
        t (i,k) = t_in (i,k)
        qv(i,k) = qv_in(i,k)
        qc(i,k) = qc_in(i,k)
        qr(i,k) = qr_in(i,k)

      END DO
    END DO

  END IF

  IF (izdebug > 25) THEN
    CALL message('mo_gscp', 'UPDATED VARIABLES')
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp T= ' ,&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
  ENDIF

!!$#ifdef NUDGING
!!$! add part of latent heating calculated in subroutine kessler_pp to model latent
!!$! heating field: add temperature to model latent heating field
!!$IF (llhn) &
!!$ CALL get_gs_lheating ('inc',1,ke)
!!$#endif

!------------------------------------------------------------------------------
! End of module  kessler_pp
!------------------------------------------------------------------------------

END SUBROUTINE kessler_pp




SUBROUTINE SATAD ( kitera, te, qve, qce, tstart, phfe,                        &
  zdqd  , zqdwe, zh   , ztg0  , ztgn, zdqdt0, zgqd0, zphe ,  &
  b1, b2w, b3, b4w, b234w, rdrd, emrdrd, rddrm1, lh_v, cpdr, &
  cp_d, idim, ilo, iup )

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
    idim,                & !  Dimension of I/O-fields
    ilo, iup               !  start- and end-indices for the computation

  REAL    (KIND=ireals),    INTENT (IN)    ::  &
    tstart  (idim), & ! Start temperature for iteration
    phfe    (idim)  ! Pressure (input)

  REAL    (KIND=ireals),    INTENT (INOUT) ::  &
    te      (idim), & ! Temperature on input/ouput
    qve     (idim), & ! Specific humidity on input/output
    qce     (idim), & ! Specific cloud water content on input/output
    zdqd    (idim), & !
    zqdwe   (idim), & !
    zh      (idim), & !
    ztg0    (idim), & !
    ztgn    (idim), & !
    zdqdt0  (idim), & !
    zgqd0   (idim), & !
    zphe    (idim)    !

  REAL    (KIND=ireals),    INTENT (IN)    ::  &
    b1, b2w, b3, b4w, b234w, rdrd, emrdrd, rddrm1, lh_v, cpdr, cp_d

! Local parameters: None
! ----------------
! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    i,                   & !  Loop indices
    nzit,                & !  Loop for iterations
    nsat,                & !  Number of saturated gridpoints
    iwrk(idim),          & !  i-index of saturated gridpoints
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

  DO i = ilo , iup

    ! "save" the start values for the temperature
    ztg0 (i) = tstart(i)

    ! correction for negative values of qv and qc
    qve (i) = MAX( qve(i), 0.0_ireals )
    qce (i) = MAX( qce(i), 0.0_ireals )

    ! assume first subsaturation
    zqdwe(i)= qve(i) + qce(i)
    te (i)  = te(i) - lh_v*qce(i)*cpdr
    qve(i)  = zqdwe(i)
    qce(i)  = 0.0_ireals
    zgeu      = fgew(te(i))
    zgqdu     = fgqd(zgeu,phfe(i))
    zdqd(i) = zgqdu - zqdwe(i)
    minzdqd   = MIN(minzdqd,zdqd(i))

  ENDDO

!NEC_CB if zdqd>=0, then for sure no points are found
  IF ( minzdqd >= 0.0_ireals ) RETURN

  DO i = ilo , iup

    IF (zdqd(i) < 0.0_ireals ) THEN
      nsat       = nsat+1
      iwrk(nsat) = i
    ENDIF

  ENDDO

  IF (nsat == 0) RETURN

! Do saturation adjustments for saturated gridpoints
! --------------------------------------------------

!cdir nodep
  DO indx = 1, nsat
     i = iwrk(indx)
     zh   (i) = cp_d*te(i) + lh_v*qve(i)
     zphe (i) = phfe(i)
     zgew       = fgew(ztg0(i))
     zgqd0(i) = fgqd(zgew,zphe(i))
  ENDDO

  IF ( kitera > 1 ) THEN
    DO  nzit  = 1 , kitera-1

!cdir nodep
      DO indx = 1, nsat
        i = iwrk(indx)
        zdqdt0(i) = fdqdt(ztg0(i),zgqd0(i))
        ztg0(i)   = (zh(i) - lh_v*(zgqd0(i)-zdqdt0(i)*ztg0(i)))/ &
                      ( cp_d + lh_v*zdqdt0(i) )
        zgew        = fgew(ztg0(i))
        zgqd0(i)  = fgqd(zgew,zphe(i))
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------

!cdir nodep
  DO indx = 1, nsat
      i = iwrk(indx)
      zdqdt0(i) = fdqdt(ztg0(i),zgqd0(i))
      ztgn(i)   = ( zh(i) - lh_v*(zgqd0(i)-zdqdt0(i)*ztg0(i)) ) / &
                    ( cp_d + lh_v*zdqdt0(i) )
      zgqd0(i)  = zgqd0(i) + zdqdt0(i)*( ztgn(i)-ztg0(i) )
  ENDDO

! Distribute the result on gridpoints
! -----------------------------------

!cdir nodep
  DO indx = 1, nsat
      i = iwrk(indx)
      te (i) =  ztgn(i)
      qve(i) = zgqd0(i)
      qce(i) = MAX( zqdwe(i) - zgqd0(i), zqwmin )
  ENDDO

! End of the subroutine

END SUBROUTINE satad


!==============================================================================

END MODULE mo_gscp_cosmo
