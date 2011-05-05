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
USE mo_atm_phy_nwp_nml,    ONLY: qi0,qc0 !!,satad
#endif

!==============================================================================

IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: hydci_pp, hydci_pp_init,satad

!------------------------------------------------------------------------------
!! Public variables
!------------------------------------------------------------------------------

PUBLIC :: cloud_num, mu_rain

!------------------------------------------------------------------------------
!> Parameters and variables which are global in this module
!------------------------------------------------------------------------------

INTEGER (KIND=iintegers), PARAMETER ::  &
  iautocon       = 1,&
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
    IF (idbg > 10) THEN
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
  ie,je,ke,                          & !> array dimensions
  istart,iend,jstart,jend,kstart,    & !! optional start/end indicies
  idbg,                              & !! optional debug level
  l_cv,                              &
  dz,zdt,                            & !! numerics parameters
  t,p,rho,qv,qc,qi,qr,qs,            & !! prognostic variables
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
    je              ,    & !! number of grid points in second (meridional) direction
    ke                     !! number of grid points in vertical direction

  INTEGER, INTENT(IN), OPTIONAL ::  &
    istart    ,    & !> optional start index for computations in the parallel program
    iend      ,    & !! optional end index for computations in the parallel program
    jstart    ,    & !! optional start index for computations in the parallel program
    jend      ,    & !! optional end index for computations in the parallel program
    kstart    ,    & !! optional start index for the vertical index
    idbg             !! optional debug level

  REAL(KIND=ireals), INTENT(IN) :: &
    zdt                    !> time step for integration of microphysics     (  s  )

  REAL(KIND=ireals), DIMENSION(ie,je,ke), INTENT(IN) ::      &
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv                   !! if true, cv is used instead of cp

  REAL(KIND=ireals), DIMENSION(ie,je,ke), INTENT(INOUT) ::   &
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    qr              ,    & !! specific rain content                         (kg/kg)
    qs                     !! specific snow content                         (kg/kg)

  REAL(KIND=ireals), DIMENSION(ie,je), INTENT(INOUT) ::   &
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    prs_gsp                !! precipitation rate of snow, grid-scale        (kg/(m2*s))

  REAL(KIND=ireals), DIMENSION(ie,je,ke), INTENT(OUT), OPTIONAL ::   &
    ddt_tend_t      , & !> tendency T                                       ( 1/s )
    ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
    ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
    ddt_tend_qi     , & !! tendency qi                                      ( 1/s )
    ddt_tend_qr     , & !! tendency qr                                      ( 1/s )
    ddt_tend_qs         !! tendency qs                                      ( 1/s )

  REAL(KIND=ireals), DIMENSION(ie,je,ke), INTENT(OUT), OPTIONAL ::   &
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
    k, i, j             !> loop indees

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
    jstartpar    ,    & !! start index for computations in the parallel program
    jendpar      ,    & !! end index for computations in the parallel program
    k_start      ,    & !! model level where computations start
    izdebug             !! debug level

  REAL(KIND=ireals), DIMENSION(ie,je,ke) ::   &
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
    zpres       (ie,je)        !! pressure
#endif

  REAL    (KIND=ireals   ) ::  &
    zqvsi       (ie,je),     & !> sat. specitic humidity at ice and water saturation
    zvzr        (ie,je),     & !!
    zvzs        (ie,je),     & !!
    zpkr        (ie,je),     & !!
    zpks        (ie,je),     & !!
    zprvr       (ie,je),     & !!
    zprvs       (ie,je),     & !!
    zdummy      (ie,je,8),   & !!
    zcsdep      (ie,je),     & !!
    zcidep      (ie,je),     & !!
    zvz0s       (ie,je),     & !!
    zcrim       (ie,je),     & !!
    zcagg       (ie,je),     & !!
    zbsdep      (ie,je),     & !!
    zcslam      (ie,je),     & !!
    zn0s        (ie,je),     & !!
    zimr        (ie,je),     & !!
    zims        (ie,je),     & !!
    zzar        (ie,je),     & !!
    zzas        (ie,je),     & !!
    zqrk        (ie,je),     & !!
    zqsk        (ie,je),     & !!
    zdtdh       (ie,je),     & !!
    z1orhog     (ie,je),     & !! 1/rhog
    zrho1o2     (ie,je),     & !! (rho0/rhog)**1/2
    zeln7o8qrk  (ie,je),     & !!
    zeln27o16qrk(ie,je),     & !!
    zeln13o8qrk (ie,je),     & !!
    zeln3o16qrk (ie,je),     & !!
    zeln13o12qsk(ie,je),     & !!
    zeln5o24qsk (ie,je),     & !!
    zeln2o3qsk  (ie,je)        !!

  REAL    (KIND=ireals   ) ::  &
    scau   (ie,je), & !> transfer rate due to autoconversion of cloud water
    scac   (ie,je), & !! transfer rate due to accretion of cloud water
    snuc   (ie,je), & !! transfer rate due nucleation of cloud ice
    scfrz  (ie,je), & !! transfer rate due homogeneous freezing of cloud water
    simelt (ie,je), & !! transfer rate due melting of cloud ice
    sidep  (ie,je), & !! transfer rate due depositional growth of cloud ice
    ssdep  (ie,je), & !! transfer rate due depositional growth of snow
    sdau   (ie,je), & !! transfer rate due depositional cloud ice autoconversion
    srim   (ie,je), & !! transfer rate due riming of snow
    sshed  (ie,je), & !! transfer rate due shedding
    sicri  (ie,je), & !! transfer rate due cloud ice collection by rain (sink qi)
    srcri  (ie,je), & !! transfer rate due cloud ice collection by rain (sink qr)
    sagg   (ie,je), & !! transfer rate due aggregation of snow and cloud ice
    siau   (ie,je), & !! transfer rate due autoconversion of cloud ice
    ssmelt (ie,je), & !! transfer rate due melting of snow
    sev    (ie,je), & !! transfer rate due evaporation of rain
    srfrz  (ie,je)    !! transfer rate due to rainwater freezing


  !> Integer arrays for a better vectorization
  INTEGER (KIND=iintegers) ::  &
    idx1(ie*je),       jdx1(ie*je),   & !!
    idx2(ie*je),       jdx2(ie*je),   & !!
    idx3(ie*je),       jdx3(ie*je),   & !!
    idx4(ie*je),       jdx4(ie*je),   & !!
    idx5(ie*je),       jdx5(ie*je),   & !!
    idx6(ie*je),       jdx6(ie*je)      !!

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
!CDIR BEGIN COLLAPSE
  prr_gsp (:,:) = 0.0_ireals
  prs_gsp (:,:) = 0.0_ireals
  zpkr    (:,:) = 0.0_ireals
  zpks    (:,:) = 0.0_ireals
  zprvr   (:,:) = 0.0_ireals
  zprvs   (:,:) = 0.0_ireals
  zvzr    (:,:) = 0.0_ireals
  zvzs    (:,:) = 0.0_ireals
  zdummy(:,:,:) = 0.0_ireals
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
  IF (PRESENT(jstart)) THEN
    jstartpar = jstart
  ELSE
    jstartpar = 1
  END IF
  IF (PRESENT(jend)) THEN
    jendpar = jend
  ELSE
    jendpar = je
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
  IF (izdebug > 15) CALL message('','SRC_GSCP: Start of hydci_pp')

  IF (izdebug > 15) THEN
    WRITE (message_text,'(A,E10.3)') '      ccslam = ',ccslam ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccsvel = ',ccsvel ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccsrim = ',ccsrim ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccsagg = ',ccsagg ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccsdep = ',ccsdep ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccslxp = ',ccslxp ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccidep = ',ccidep ; CALL message('',message_text)
  ENDIF
  IF (izdebug > 20) THEN
    WRITE (message_text,*) '   ie = ',ie ; CALL message('',message_text)
    WRITE (message_text,*) '   je = ',je ; CALL message('',message_text)
    WRITE (message_text,*) '   ke = ',ke ; CALL message('',message_text)
    WRITE (message_text,*) '   istartpar = ',istartpar ; CALL message('',message_text)
    WRITE (message_text,*) '   iendpar   = ',iendpar   ; CALL message('',message_text)
    WRITE (message_text,*) '   jstartpar = ',jstartpar ; CALL message('',message_text)
    WRITE (message_text,*) '   jendpar   = ',jendpar   ; CALL message('',message_text)
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
      tinc_lh(:,:,k) = tinc_lh(:,:,k) - t(:,:,k)
    ENDIF
#endif

  !----------------------------------------------------------------------------
  ! Section 2: Check for existence of rain and snow
  !            Initialize microphysics and sedimentation scheme
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
    zcrim (:,:) = 0.0_ireals
    zcagg (:,:) = 0.0_ireals
    zbsdep(:,:) = 0.0_ireals
    zvz0s (:,:) = 0.0_ireals
    zn0s  (:,:) = zn0s0
!CDIR END

    ic1 = 0
    ic2 = 0

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar

        qrg = qr(i,j,k)
        qsg = qs(i,j,k)
        qvg = qv(i,j,k)
        qcg = qc(i,j,k)
        qig = qi(i,j,k)
        tg  = t(i,j,k)
        ppg = p(i,j,k)
        rhog = rho(i,j,k)

        !..for density correction of fall speeds
        z1orhog(i,j) = 1.0_ireals/rhog
        zrho1o2(i,j) = EXP(LOG(zrho0*z1orhog(i,j))*x1o2)

        zqrk(i,j) = qrg * rhog
        zqsk(i,j) = qsg * rhog

        llqr = zqrk(i,j) > zqmin
        llqs = zqsk(i,j) > zqmin

        zdtdh(i,j) = 0.5_ireals * zdt / dz(i,j,k)

        zzar(i,j)   = zqrk(i,j)/zdtdh(i,j) + zprvr(i,j) + zpkr(i,j)
        zzas(i,j)   = zqsk(i,j)/zdtdh(i,j) + zprvs(i,j) + zpks(i,j)

        IF (llqs) THEN
          ic1 = ic1 + 1
          idx1(ic1) = i
          jdx1(ic1) = j
        ENDIF
        IF (llqr) THEN
          ic2 = ic2 + 1
          idx2(ic2) = i
          jdx2(ic2) = j
        ENDIF
      ENDDO
    ENDDO

#ifdef __ICON__
    ! j is always 1 in ICON
    j = 1
#endif

!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs_prepare: DO i1d = 1, ic1
      i = idx1(i1d)
#ifndef __ICON__
      j = jdx1(i1d)
#endif
      qsg = qs(i,j,k)
      tg  = t(i,j,k)

      IF (isnow_n0temp == 1) THEN
        ! Calculate n0s using the temperature-dependent
        ! formula of Field et al. (2005)
        ztc = tg - t0
        ztc = MAX(MIN(ztc,0.0_ireals),-40.0_ireals)
        zn0s(i,j) = zn0s1*EXP(zn0s2*ztc)
        zn0s(i,j) = MIN(zn0s(i,j),1e9_ireals)
        zn0s(i,j) = MAX(zn0s(i,j),1e6_ireals)
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
        m2s = qsg / zams
        m3s = alf*EXP(bet*LOG(m2s))

        hlp  = zn0s1*EXP(zn0s2*ztc)
        zn0s(i,j) = 13.50_ireals * m2s**4 / m3s**3
        zn0s(i,j) = MAX(zn0s(i,j),0.5_ireals*hlp)
        zn0s(i,j) = MIN(zn0s(i,j),1e2_ireals*hlp)
        zn0s(i,j) = MIN(zn0s(i,j),1e9_ireals)
        zn0s(i,j) = MAX(zn0s(i,j),1e6_ireals)
      ELSE
        ! Old constant n0s
        zn0s(i,j) = 8.0e5_ireals
      ENDIF
      zcrim (i,j) = ccsrim*zn0s(i,j)
      zcagg (i,j) = ccsagg*zn0s(i,j)
      zbsdep(i,j) = ccsdep*SQRT(v0snow)
      zvz0s (i,j) = ccsvel*EXP(ccsvxp * LOG(zn0s(i,j)))

      IF (zvzs(i,j) == 0.0_ireals) THEN
        zvzs(i,j) = zvz0s(i,j) * EXP (ccswxp * LOG (0.5_ireals*zqsk(i,j))) * zrho1o2(i,j)
      ENDIF
    ENDDO loop_over_qs_prepare

!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qr_sedi: DO i1d = 1, ic2
      i = idx2(i1d)
#ifndef __ICON__
      j = jdx2(i1d)
#endif
      IF (zvzr(i,j) == 0.0_ireals) THEN
        zvzr(i,j) = zvz0r * EXP (x1o8  * LOG (0.5_ireals*zqrk(i,j))) * zrho1o2(i,j)
      ENDIF
    ENDDO loop_over_qr_sedi

  !----------------------------------------------------------------------------
  ! Section 3:
  !----------------------------------------------------------------------------

!CDIR BEGIN COLLAPSE
    zeln7o8qrk   (:,:) = 0.0_ireals
    zeln27o16qrk (:,:) = 0.0_ireals
    zeln13o8qrk  (:,:) = 0.0_ireals
    zeln3o16qrk  (:,:) = 0.0_ireals
    zeln13o12qsk (:,:) = 0.0_ireals
    zeln5o24qsk  (:,:) = 0.0_ireals
    zeln2o3qsk   (:,:) = 0.0_ireals
    zcsdep       (:,:) = 3.2E-2_ireals
    zcidep       (:,:) = 1.3E-5_ireals
    zcslam       (:,:) = 1e10_ireals

    scau         (:,:) = 0.0_ireals
    scac         (:,:) = 0.0_ireals
    snuc         (:,:) = 0.0_ireals
    scfrz        (:,:) = 0.0_ireals
    simelt       (:,:) = 0.0_ireals
    sidep        (:,:) = 0.0_ireals
    ssdep        (:,:) = 0.0_ireals
    sdau         (:,:) = 0.0_ireals
    srim         (:,:) = 0.0_ireals
    sshed        (:,:) = 0.0_ireals
    sicri        (:,:) = 0.0_ireals
    srcri        (:,:) = 0.0_ireals
    sagg         (:,:) = 0.0_ireals
    siau         (:,:) = 0.0_ireals
    ssmelt       (:,:) = 0.0_ireals
    sev          (:,:) = 0.0_ireals
    srfrz        (:,:) = 0.0_ireals
!CDIR END

    ic1 = 0
    ic2 = 0
    ic3 = 0
    ic4 = 0
    ic5 = 0
    ic6 = 0
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar

        qrg  = qr(i,j,k)
        qsg  = qs(i,j,k)
        qvg  = qv(i,j,k)
        qcg  = qc(i,j,k)
        qig  = qi(i,j,k)
        tg   =  t(i,j,k)
        ppg  =  p(i,j,k)
        rhog = rho(i,j,k)

!        zqvsi(i,j) = spec_humi( sat_pres_ice(tg), ppg )
        zqvsi(i,j) = sat_pres_ice(tg)/(rhog * r_v * tg)

        llqr = zqrk(i,j) > zqmin
        llqs = zqsk(i,j) > zqmin

        IF (llqr) THEN
          zpkr(i,j) = zqrk(i,j) * zvz0r * EXP (x1o8  * LOG (zqrk(i,j))) * zrho1o2(i,j)
        ELSE
          zpkr(i,j) = 0.0_ireals
        ENDIF

        IF (llqs) THEN
          zpks(i,j) = zqsk (i,j) * zvz0s(i,j) * EXP (ccswxp * LOG (zqsk(i,j))) * zrho1o2(i,j)
        ELSE
          zpks(i,j) = 0.0_ireals
        ENDIF

        zpkr(i,j)   = MIN( zpkr(i,j) , zzar(i,j) )
        zpks(i,j)   = MIN( zpks(i,j) , zzas(i,j) )

        zzar(i,j)   = zdtdh(i,j) * (zzar(i,j)-zpkr(i,j))
        zzas(i,j)   = zdtdh(i,j) * (zzas(i,j)-zpks(i,j))

        zimr(i,j)   = 1.0_ireals / (1.0_ireals + zvzr(i,j) * zdtdh(i,j))
        zims(i,j)   = 1.0_ireals / (1.0_ireals + zvzs(i,j) * zdtdh(i,j))

        zqrk(i,j)   = zzar(i,j)*zimr(i,j)
        zqsk(i,j)   = zzas(i,j)*zims(i,j)

        llqr = zqrk(i,j) > zqmin
        llqs = zqsk(i,j) > zqmin
        llqc =       qcg > zqmin
        llqi =       qig > zqmin

        IF (llqr) THEN
          ic1 = ic1 + 1
          idx1(ic1) = i
          jdx1(ic1) = j
        ENDIF
        IF (llqs) THEN
          ic2 = ic2 + 1
          idx2(ic2) = i
          jdx2(ic2) = j
        ENDIF
        IF (llqi .OR. llqs) THEN
          ic3 = ic3 + 1
          idx3(ic3) = i
          jdx3(ic3) = j
        ENDIF
        IF ( tg < zthet .AND. qvg >  8.E-6_ireals &
                        .AND. qig <= 0.0_ireals ) THEN
          ic4 = ic4 + 1
          idx4(ic4) = i
          jdx4(ic4) = j
        ENDIF
        IF (llqc) THEN
          ic5 = ic5 + 1
          idx5(ic5) = i
          jdx5(ic5) = j
        ENDIF
        IF (llqr .AND. qcg <= 0.0_ireals) THEN
          ic6 = ic6 + 1
          idx6(ic6) = i
          jdx6(ic6) = j
        ENDIF
      ENDDO
    ENDDO

#ifdef __ICON__
    ! j is always 1 in ICON
    j = 1
#endif

! ic1
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qr: DO i1d =1, ic1
      i = idx1(i1d)
#ifndef __ICON__
      j = jdx1(i1d)
#endif
      qcg  = qc(i,j,k)
      qig  = qi(i,j,k)
      tg   =  t(i,j,k)
      llqi =  qig > zqmin

      zlnqrk       = LOG (zqrk(i,j))
      IF ( qig+qcg > zqmin ) THEN
        zeln7o8qrk(i,j)   = EXP (x7o8   * zlnqrk)
      ENDIF
      IF ( tg < ztrfrz ) THEN
        zeln27o16qrk(i,j) = EXP (x27o16 * zlnqrk)
      ENDIF
      IF (llqi) THEN
        zeln13o8qrk(i,j)  = EXP (x13o8  * zlnqrk)
      ENDIF
      IF (qcg <= 0.0_ireals ) THEN
        zeln3o16qrk(i,j)  = EXP (x3o16  * zlnqrk)
      ENDIF
    ENDDO loop_over_qr

! ic2
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qs_coeffs: DO i1d =1, ic2
      i = idx2(i1d)
#ifndef __ICON__
      j = jdx2(i1d)
#endif
      qcg = qc(i,j,k)
      qig = qi(i,j,k)

      zlnqsk       = LOG (zqsk(i,j))
      IF (qig+qcg > zqmin) THEN
        zeln13o12qsk(i,j) = EXP (x13o12 * zlnqsk)
      ENDIF
      zeln5o24qsk(i,j)  = EXP (x5o24  * zlnqsk)
      zeln2o3qsk(i,j)   = EXP (x2o3   * zlnqsk)
    ENDDO loop_over_qs_coeffs

! ic3
!CDIR NODEP,VOVERTAKE,VOB
    loop_over_qi_qs: DO i1d =1, ic3
      i = idx3(i1d)
#ifndef __ICON__
      j = jdx3(i1d)
#endif
      tg   =   t(i,j,k)
      ppg  =   p(i,j,k)
      rhog = rho(i,j,k)
      llqs = zqsk(i,j) > zqmin

      zdvtp  = ccdvtp * EXP(1.94_ireals * LOG(tg)) / ppg
      zhi    = ccshi1*zdvtp*rhog*zqvsi(i,j)/(tg*tg)
      hlp    = zdvtp / (1.0_ireals + zhi)
      zcidep(i,j) = ccidep * hlp
      IF (llqs) THEN
        zcslam(i,j) = EXP(ccslxp * LOG(ccslam * zn0s(i,j) / zqsk(i,j) ))
        zcslam(i,j) = MIN(zcslam(i,j),1e15_ireals)
        zcsdep(i,j) = 4.0_ireals * zn0s(i,j) * hlp
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
#ifndef __ICON__
      j = jdx4(i1d)
#endif
      qvg  =  qv(i,j,k)
      tg   =   t(i,j,k)
      ppg  =   p(i,j,k)
      rhog = rho(i,j,k)

      IF( qvg > zqvsi(i,j) ) THEN
        znin  = MIN( ice_nuclei_number(tg), znimax )
        zsnuc = zmi0 / rhog * znin * zdtr
        snuc(i,j) = zsnuc
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
#ifndef __ICON__
      j = jdx5(i1d)
#endif
      qrg  =   qr(i,j,k)
      qvg  =   qv(i,j,k)
      qcg  =   qc(i,j,k)
      qig  =   qi(i,j,k)
      tg   =    t(i,j,k)
      ppg  =    p(i,j,k)
      rhog =  rho(i,j,k)
      llqs = zqsk(i,j) > zqmin

      if (iautocon == 0) THEN
        ! Kessler (1969) autoconversion rate
        zscau  = zccau * MAX( qcg - qc0, 0.0_ireals )
        zscac  = zcac  * qcg * zeln7o8qrk(i,j)
      ELSEIF (iautocon == 1) THEN
        ! Seifert and Beheng (2001) autoconversion rate
        ! with constant cloud droplet number concentration cloud_num
        IF (qcg > 1.0e-6_ireals) THEN
          ztau   = MIN(1.0_ireals-qcg/(qcg+qrg),0.9_ireals)
          zphi   = zkphi1 * ztau**zkphi2 * (1.0_ireals - ztau**zkphi2)**3
          zscau  = zconst * qcg*qcg*qcg*qcg &
            &    * (1.0_ireals + zphi/(1.0_ireals - ztau)**2)
          zphi   = (ztau/(ztau+zkphi3))**4
          zscac  = zkcac * qcg * qrg * zphi !* zrho1o2(i,j)
        ELSE
          zscau  = 0.0_ireals
          zscac  = 0.0_ireals
        ENDIF
      ENDIF
      IF (llqs) THEN
        zscrim = zcrim(i,j) * EXP(ccsaxp * LOG(zcslam(i,j))) * qcg !* zrho1o2(i,j)
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
        scfrz(i,j) = zscmax
      ELSE
        scau (i,j) = zcorr*zscau
        scac (i,j) = zcorr*zscac
        srim (i,j) = zcorr*zscrim
        sshed(i,j) = zcorr*zscshe
      ENDIF

      ! Calculation of heterogeneous nucleation of cloud ice.
      ! This is done in this section, because we require water saturation
      ! for this process (i.e. the existence of cloud water) to exist.
      ! Hetrogeneous nucleation is assumed to occur only when no
      ! cloud ice is present and the temperature is below a nucleation
      ! threshold.
      IF( tg <= 267.15_ireals .AND. qig <= 0.0_ireals ) THEN
        znin  = MIN( ice_nuclei_number(tg), znimax )
        zsnuc = zmi0 * z1orhog(i,j) * znin * zdtr
        snuc(i,j) = zsnuc
      ENDIF
      ! Calculation of in-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        zsrfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 )* zeln27o16qrk(i,j)
        srfrz(i,j) = zsrfrz
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
#ifndef __ICON__
      j = jdx3(i1d)
#endif
      qvg  =  qv(i,j,k)
      qig  =  qi(i,j,k)
      qsg  =  qs(i,j,k)
      tg   =   t(i,j,k)
      ppg  =   p(i,j,k)
      rhog = rho(i,j,k)
      llqi =  qig > zqmin

      IF (tg<=t0) THEN
        znin    = MIN( ice_nuclei_number(tg), znimax )
        zmi     = MIN( rhog*qig/znin, zmimax )
        zmi     = MAX( zmi0, zmi )
        zsvmax  = (qvg - zqvsi(i,j)) * zdtr
        zsagg   = zcagg(i,j) * EXP(ccsaxp*LOG(zcslam(i,j))) * qig
        zsagg   = MAX( zsagg, 0.0_ireals ) & !* zrho1o2(i,j) &
          * MAX(0.2_ireals,MIN(EXP(0.09_ireals*(tg-t0)),1.0_ireals))
        znid      = rhog * qig/zmi
        IF (llqi) THEN
          zlnlogmi= LOG (zmi)
          zsidep    = zcidep(i,j) * znid * EXP(0.33_ireals * zlnlogmi)   &
                        * ( qvg - zqvsi(i,j) )
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
        zsicri    = zcicri * qig * zeln7o8qrk(i,j)
        zsrcri    = zcrcri * (qig/zmi) * zeln13o8qrk(i,j)
        zxfac     = 1.0_ireals + zbsdep(i,j) * EXP(ccsdxp*LOG(zcslam(i,j)))
        zssdep    = zcsdep(i,j) * zxfac * ( qvg - zqvsi(i,j) ) / (zcslam(i,j)+zeps)**2

        ! Check for maximal depletion of vapor by sdep
        IF (zssdep > 0.0_ireals) zssdep = MIN(zssdep, zsvmax-zsvidep)
        ! Check for maximal depletion of snow by sdep
        IF (zssdep < 0.0_ireals) zssdep = MAX(zssdep, -qsg*zdtr)

        zsisum = zsiau + zsdau + zsagg + zsicri + zsvisub
        zcorr  = 0.0_ireals
        IF( zsimax > 0.0_ireals ) zcorr  = zsimax / MAX( zsimax, zsisum )
        sidep(i,j)  = zsvidep - zcorr*zsvisub
        sdau (i,j)  = zcorr*zsdau
        siau (i,j)  = zcorr*zsiau
        sagg (i,j)  = zcorr*zsagg
        ssdep(i,j)  = zssdep
        srcri(i,j)  = zsrcri
        sicri(i,j)  = zcorr*zsicri

      !------------------------------------------------------------------------
      ! Section 7: Search for warm grid points with cloud ice and/or snow and
      !            calculation of the melting rates of qi and ps
      !------------------------------------------------------------------------

      ELSE ! tg > 0
        simelt(i,j) = qig*zdtr
!        zqvsw0      = spec_humi( zpvsw0, ppg)
        zqvsw0      = zpvsw0 / (rhog * r_v * tg)
        zx1         = (tg - t0) + zasmel*(qvg - zqvsw0)
        zx2         = 1.0_ireals + zbsmel * zeln5o24qsk(i,j)
        zssmelt     = zcsmel * zx1 * zx2 * zeln2o3qsk(i,j)
        ssmelt(i,j) = MAX( zssmelt, 0.0_ireals )
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
#ifndef __ICON__
      j = jdx6(i1d)
#endif
      qvg = qv(i,j,k)
      tg  =  t(i,j,k)
      ppg =  p(i,j,k)
      rhog = rho(i,j,k)

!      zqvsw    = spec_humi( sat_pres_water(tg), ppg )
      zqvsw    = sat_pres_water(tg)/(rhog * r_v *tg)
      zx1      = 1.0_ireals + zbev* zeln3o16qrk(i,j)
      zsev     = zcev*zx1*(zqvsw - qvg)*SQRT(zqrk(i,j))
      sev(i,j) = MAX( zsev, 0.0_ireals )
      ! Calculation of below-cloud rainwater freezing
      IF ( tg < ztrfrz ) THEN
        zsrfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 ) * zeln27o16qrk(i,j)
        srfrz(i,j)  = zsrfrz
      ENDIF
    ENDDO loop_over_qr_nocloud

    !--------------------------------------------------------------------------
    ! Section 9: Calculate the total tendencies of the prognostic variables.
    !            Update the prognostic variables in the interior domain.
    !--------------------------------------------------------------------------

    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar

        qvg = qv(i,j,k)
        qcg = qc(i,j,k)
        qrg = qr(i,j,k)
        qsg = qs(i,j,k)
        qig = qi(i,j,k)
        tg  = t (i,j,k)
        rhog = rho(i,j,k)

        zsrmax = zzar(i,j)*z1orhog(i,j)*zdtr
        zssmax = zzas(i,j)*z1orhog(i,j)*zdtr
        zsrsum = sev(i,j) + srfrz(i,j) + srcri(i,j)
        zcorr  = 1.0_ireals
        IF(zsrsum > 0.0_ireals) THEN
          zcorr  = zsrmax / MAX( zsrmax, zsrsum )
        ENDIF
        sev  (i,j) = zcorr*sev(i,j)
        srfrz(i,j) = zcorr*srfrz(i,j)
        srcri(i,j) = zcorr*srcri(i,j)

        ssmelt(i,j) = MIN(ssmelt(i,j), zssmax)
        IF (ssdep(i,j) < 0.0_ireals ) THEN
          ssdep(i,j) = MAX(ssdep(i,j), - zssmax)
        ENDIF
        zqvt = sev(i,j)   - sidep(i,j) - ssdep(i,j)  - snuc(i,j)
        zqct = simelt(i,j)- scau(i,j)  - scfrz(i,j)  - scac(i,j)   - sshed(i,j) - srim(i,j)
        zqit = snuc(i,j)  + scfrz(i,j) - simelt(i,j) - sicri(i,j)  + sidep(i,j) - sdau(i,j)  &
             & - sagg(i,j) - siau(i,j)
        zqrt = scau(i,j)  + sshed(i,j) + scac(i,j)   + ssmelt(i,j) - sev(i,j)   - srcri(i,j) &
             & - srfrz(i,j)
        zqst = siau(i,j)  + sdau(i,j)  + sagg(i,j)   - ssmelt(i,j) + sicri(i,j) + srcri(i,j) &
             & + srim(i,j) + ssdep(i,j) + srfrz(i,j)
!        ztt = cpdr*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst) )
        ztt = z_heat_cap_r*( lh_v*(zqct+zqrt) + lh_s*(zqit+zqst) )

        ! Update local variables
        qig = MAX ( 0.0_ireals, qig + zqit*zdt)
        qrg = MAX ( 0.0_ireals, (zzar(i,j)/rhog + zqrt*zdt)*zimr(i,j))
        qsg = MAX ( 0.0_ireals, (zzas(i,j)/rhog + zqst*zdt)*zims(i,j))
        qvg = MAX ( 0.0_ireals, qvg + zqvt*zdt )
        qcg = MAX ( 0.0_ireals, qcg + zqct*zdt )
        tg  = tg + ztt*zdt

        !----------------------------------------------------------------------
        ! Section 10: Complete time step
        !----------------------------------------------------------------------

        IF ( k /= ke) THEN
          ! Store precipitation fluxes and sedimentation velocities
          ! for the next level
          zprvr(i,j) = qrg*rhog*zvzr(i,j)
          zprvs(i,j) = qsg*rhog*zvzs(i,j)
          IF (zprvr(i,j) <= zqmin) zprvr(i,j)=0.0_ireals
          IF (zprvs(i,j) <= zqmin) zprvs(i,j)=0.0_ireals

          IF (qrg+qr(i,j,k+1) <= zqmin) THEN
            zvzr(i,j)= 0.0_ireals
          ELSE
            zvzr(i,j) = zvz0r                                               &
                 &      * EXP(x1o8 *LOG((qrg+qr(i,j,k+1))*0.5_ireals*rhog)) &
                 &      * zrho1o2(i,j)
          ENDIF
          IF (qsg+qs(i,j,k+1) <= zqmin) THEN
            zvzs(i,j)= 0.0_ireals
          ELSE
            zvzs(i,j) = zvz0s(i,j)                                                           &
                 &      * EXP(zv1s/(zbms+1.0_ireals)*LOG((qsg+qs(i,j,k+1))*0.5_ireals*rhog)) &
                 &      * zrho1o2(i,j)
          ENDIF
        ELSE
          ! Precipitation fluxes at the ground
          prr_gsp(i,j) = 0.5_ireals * (qrg*rhog*zvzr(i,j) + zpkr(i,j))
          prs_gsp(i,j) = 0.5_ireals * (qsg*rhog*zvzs(i,j) + zpks(i,j))

        ENDIF

        ! Update of prognostic variables or tendencies
        qr (i,j,k) = qrg
        qs (i,j,k) = qsg
        qi (i,j,k) = qig
        t  (i,j,k) = tg
        qv (i,j,k) = qvg
        qc (i,j,k) = qcg

        ! Store optional microphysical rates for diagnostics
        IF (PRESENT(ddt_diag_au   )) ddt_diag_au   (i,j,k) = scau  (i,j)
        IF (PRESENT(ddt_diag_ac   )) ddt_diag_ac   (i,j,k) = scac  (i,j)
        IF (PRESENT(ddt_diag_ev   )) ddt_diag_ev   (i,j,k) = sev   (i,j)
        IF (PRESENT(ddt_diag_nuc  )) ddt_diag_nuc  (i,j,k) = snuc  (i,j)
        IF (PRESENT(ddt_diag_idep )) ddt_diag_idep (i,j,k) = sidep (i,j)
        IF (PRESENT(ddt_diag_sdep )) ddt_diag_sdep (i,j,k) = ssdep (i,j)
        IF (PRESENT(ddt_diag_agg  )) ddt_diag_agg  (i,j,k) = sagg  (i,j)
        IF (PRESENT(ddt_diag_rim  )) ddt_diag_rim  (i,j,k) = srim  (i,j)
        IF (PRESENT(ddt_diag_rcri )) ddt_diag_rcri (i,j,k) = srcri (i,j)
        IF (PRESENT(ddt_diag_icri )) ddt_diag_icri (i,j,k) = sicri (i,j)
        IF (PRESENT(ddt_diag_dau  )) ddt_diag_dau  (i,j,k) = sdau  (i,j)
        IF (PRESENT(ddt_diag_iau  )) ddt_diag_iau  (i,j,k) = siau  (i,j)
        IF (PRESENT(ddt_diag_imelt)) ddt_diag_imelt(i,j,k) = simelt(i,j)
        IF (PRESENT(ddt_diag_smelt)) ddt_diag_smelt(i,j,k) = ssmelt(i,j)
        IF (PRESENT(ddt_diag_cfrz )) ddt_diag_cfrz (i,j,k) = scfrz (i,j)
        IF (PRESENT(ddt_diag_rfrz )) ddt_diag_rfrz (i,j,k) = srfrz (i,j)
        IF (PRESENT(ddt_diag_shed )) ddt_diag_shed (i,j,k) = sshed (i,j)

      ENDDO
    ENDDO

  IF (izdebug > 15) THEN
    ! Check for negative values
    DO j = jstartpar, jendpar
      DO i = istartpar, iendpar
        IF (qr(i,j,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qr'
          CALL message('',message_text)
        ENDIF
        IF (qc(i,j,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qc'
          CALL message('',message_text)
        ENDIF
        IF (qi(i,j,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qi'
          CALL message('',message_text)
        ENDIF
        IF (qs(i,j,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qs'
          CALL message('',message_text)
        ENDIF
        IF (qv(i,j,k) < 0.0_ireals) THEN
          WRITE(message_text,'(a)') ' WARNING: hydci_pp, negative value in qv'
          CALL message('',message_text)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

#if defined (__COSMO__)
  ! Do a final saturation adjustment for new values of t, qv and qc
!CDIR COLLAPSE
    zpres(:,:) = p(:,:,k)

    CALL satad ( 1, t(:,:,k), qv(:,:,k),              &
               qc(:,:,k), t(:,:,k), zpres,          &
               zdummy(:,:,1),zdummy(:,:,2),zdummy(:,:,3), &
               zdummy(:,:,4),zdummy(:,:,5),zdummy(:,:,6), &
               zdummy(:,:,7),zdummy(:,:,8),               &
               b1, b2w, b3, b4w, b234w, rdv, o_m_rdv,     &
               rvd_m_o, lh_v, cpdr, cp_d,                 &
               ie, je, istartpar, iendpar, jstartpar, jendpar )

  IF ( ldiabf_lh ) THEN
    ! compute temperature increment due to latent heat
!CDIR COLLAPSE
    tinc_lh(:,:,k) = tinc_lh(:,:,k) + t(:,:,k)
  ENDIF
#endif

ENDDO loop_over_levels

#ifdef __ICON__

 CALL satad_v_3d (                             &
               & maxiter  = 10_iintegers ,& !> IN
               & tol      = 1.e-3_ireals ,& !> IN
               & te       = t  (1,1,1)   ,&
               & qve      = qv (1,1,1)   ,&
               & qce      = qc (1,1,1)   ,&
               & rhotot   = rho(1,1,1)   ,&
               & idim     = ie           ,&
               & jdim     = je           ,&
               & kdim     = ke           ,&
               & ilo      = istartpar    ,&
               & iup      = iendpar      ,&
               & jlo      = jstartpar    ,&
               & jup      = jendpar      ,&
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
      DO j=jstartpar,jendpar
        DO i=istartpar,iendpar

          ! calculated pseudo-tendencies
          ddt_tend_t (i,j,k) = (t (i,j,k) - t_in (i,j,k))*zdtr
          ddt_tend_qv(i,j,k) = MAX(-qv_in(i,j,k)*zdtr,(qv(i,j,k) - qv_in(i,j,k))*zdtr)
          ddt_tend_qc(i,j,k) = MAX(-qc_in(i,j,k)*zdtr,(qc(i,j,k) - qc_in(i,j,k))*zdtr)
          ddt_tend_qr(i,j,k) = MAX(-qr_in(i,j,k)*zdtr,(qr(i,j,k) - qr_in(i,j,k))*zdtr)
          ddt_tend_qs(i,j,k) = MAX(-qs_in(i,j,k)*zdtr,(qs(i,j,k) - qs_in(i,j,k))*zdtr)
          ddt_tend_qi(i,j,k) = MAX(-qi_in(i,j,k)*zdtr,(qi(i,j,k) - qi_in(i,j,k))*zdtr)

          ! restore input values
          t (i,j,k) = t_in (i,j,k)
          qv(i,j,k) = qv_in(i,j,k)
          qc(i,j,k) = qc_in(i,j,k)
          qi(i,j,k) = qi_in(i,j,k)
          qr(i,j,k) = qr_in(i,j,k)
          qs(i,j,k) = qs_in(i,j,k)

        END DO
      END DO
    END DO

  END IF

  IF (izdebug > 15) THEN
    CALL message('mo_gscp', 'UPDATED VARIABLES')
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp T= ',&
    MAXVAL( t(:,:,:)), MINVAL(t(:,:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp qv= ',&
    MAXVAL( qv(:,:,:)), MINVAL(qv(:,:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp qc= ',&
    MAXVAL( qc(:,:,:)), MINVAL(qc(:,:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp qi= ',&
    MAXVAL( qi(:,:,:)), MINVAL(qi(:,:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp qr= ',&
    MAXVAL( qr(:,:,:)), MINVAL(qr(:,:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'hydci_pp qs= ',&
    MAXVAL( qs(:,:,:)), MINVAL(qs(:,:,:) )
    CALL message('', TRIM(message_text))
  ENDIF


!------------------------------------------------------------------------------
! End of subroutine hydci_pp
!------------------------------------------------------------------------------

END SUBROUTINE hydci_pp


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


!==============================================================================

END MODULE mo_gscp_cosmo
