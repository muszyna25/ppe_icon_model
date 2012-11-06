!>
!! This module provides physical constants for the ICON general circulation models.
!!
!! Physical constants are grouped as follows:
!! - Natural constants
!! - Molar weights
!! - Earth and Earth orbit constants
!! - Thermodynamic constants for the dry and moist atmosphere
!! - Constants used for the computation of lookup tables of the saturation
!!    mixing ratio over liquid water (*c_les*) or ice(*c_ies*)
!!    (to be shifted to the module that computes the lookup tables)

!!
!! @par Revision History
!!  Developed  by Luis Kornblueh and Luca Bonaventura (2002-03)
!!  Modified to ProTeX-style by  Luca Bonaventura and Thomas Heinze (2004).
!!  Modified according to style guide by Thomas Heinze (2005-06-24):
!!   - module renamed from mo_constants to mo_physical_constants
!!   - eps moved to mo_math_constants
!!   - su0 renamed to u0 (as in Williamson et al. (1992) paper)
!!  Adding units to comments by Thomas Heinze (2005-07-26)
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Hui Wan (2007-01-12):
!!  - added more constants from ECHAM5
!!  Modification by Hui Wan (2007-01-16):
!!  - parameter u0 moved to <i>mo_sw_testcases</i>
!!  Modification by Kristina Froehlich (2010-05-07):
!!  - start to introduce consolidated physical constants
!!  Modification by Marco Giorgetta (2010-07-16):
!!  - improve comments and structure
!!  - add "day" for length of day in [s]
!!  - move constants "*_bg" for the background structure of the atmosphere of the nh-model
!!    to mo_vertical_grid, which is the only module using these constants
!!  - move "grav_o_cpd" to mo_divergent_modes.f90, which is th only module using this
!!    derived constant
!!  Modification by Felix Rieper (2012-02)
!!  - added some constants needed in cloud physics scheme
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_physical_constants

  USE mo_kind,            ONLY: wp

  IMPLICIT NONE

  PUBLIC

  ! Natural constants
  ! -----------------
  !
!!$  ! ECHAM values
!!$  REAL(wp), PARAMETER :: avo   = 6.022045e23_wp   ! [1/mo]    Avogadro constant
!!$  REAL(wp), PARAMETER :: ak    = 1.380662e-23_wp  ! [J/K]     Boltzmann constant
!!$  REAL(wp), PARAMETER :: argas = 8.314409_wp      ! [J/K/mol] molar/universal/ideal gas constant
!!$  REAL(wp), PARAMETER :: stbo  = 5.67E-8_wp       ! [W/m2/K4] Stephan-Boltzmann constant
  ! WMO/SI values
  REAL(wp), PARAMETER :: avo   = 6.02214179e23_wp !> [1/mo]    Avogadro constant
  REAL(wp), PARAMETER :: ak    = 1.3806504e-23_wp !! [J/K]     Boltzmann constant
  REAL(wp), PARAMETER :: argas = 8.314472_wp      !! [J/K/mol] molar/universal/ideal gas constant
  REAL(wp), PARAMETER :: stbo  = 5.6704E-8_wp     !! [W/m2/K4] Stephan-Boltzmann constant


  !> Molar weights
  !! -------------
  !!
  !! Pure species
  REAL(wp), PARAMETER :: amco2 = 44.011_wp        !>[g/mol] CO2
  REAL(wp), PARAMETER :: amch4 = 16.043_wp        !! [g/mol] CH4
  REAL(wp), PARAMETER :: amo3  = 47.9982_wp       !! [g/mol] O3
  REAL(wp), PARAMETER :: amo2  = 31.9988_wp       !! [g/mol] O2
  REAL(wp), PARAMETER :: amn2o = 44.013_wp        !! [g/mol] N2O
  REAL(wp), PARAMETER :: amc11 =137.3686_wp       !! [g/mol] CFC11
  REAL(wp), PARAMETER :: amc12 =120.9140_wp       !! [g/mol] CFC12
  REAL(wp), PARAMETER :: amw   = 18.0154_wp       !! [g/mol] H2O
  !
  !> Mixed species
  REAL(wp), PARAMETER :: amd   = 28.970_wp        !> [g/mol] dry air
  !
  !> Auxiliary constants
  ! ppmv2gg converts ozone from volume mixing ratio in ppmv
  ! to mass mixing ratio in g/g
  REAL(wp), PARAMETER :: ppmv2gg=1.e-6_wp*amo3/amd

  !> Earth and Earth orbit constants
  !! -------------------------------
  !!
!   REAL(wp), PARAMETER :: re    = 6.371229e6_wp    !! [m]    average radius
!   REAL(wp), PARAMETER :: rre   = 1._wp/ph_re         !! [1/m]
  REAL(wp), PARAMETER :: earth_radius           = 6.371229e6_wp    !! [m]    average radius
  REAL(wp), PARAMETER :: inverse_earth_radius   = 1._wp/earth_radius         !! [1/m]
  REAL(wp), PARAMETER :: earth_angular_velocity = 7.29212e-5_wp    !! [1/s]  angular velocity
  !
!!$  ! ECHAM values
!!  REAL(wp), PARAMETER :: grav  = 9.80616_wp       ! [m/s2] av. gravitational acceleration
  ! WMO/SI value
  REAL(wp), PARAMETER :: grav  = 9.80665_wp       !> [m/s2] av. gravitational acceleration
  REAL(wp), PARAMETER :: rgrav = 1._wp/grav       !! [s2/m]
  !
  REAL(wp), PARAMETER :: rae   = 0.1277E-2_wp     !> [m/m]  ratio of atm. scale height
  !                                               !!        to Earth radius


  ! Thermodynamic constants for the dry and moist atmosphere
  ! --------------------------------------------------------
  !
  !> Dry air
  REAL(wp), PARAMETER :: rd    = 287.04_wp        !> [J/K/kg] gas constant
  REAL(wp), PARAMETER :: cpd   = 1004.64_wp       !! [J/K/kg] specific heat at constant pressure
  REAL(wp), PARAMETER :: cvd   = cpd-rd           !! [J/K/kg] specific heat at constant volume
  REAL(wp), PARAMETER :: con_m = 1.50E-5_wp       !! [m^2/s]  kinematic viscosity of dry air
  REAL(wp), PARAMETER :: con_h = 2.20E-5_wp       !! [m^2/s]  scalar conductivity of dry air
  REAL(wp), PARAMETER :: con0_h= 2.40e-2_wp       !! [J/m/s/K]thermal conductivity of dry air
  REAL(wp), PARAMETER :: eta0d = 1.717e-5_wp      !! [N*s/m2] dyn viscosity of dry air at tmelt
  !
  !> H2O
  !! - gas
  REAL(wp), PARAMETER :: rv    = 461.51_wp        !> [J/K/kg] gas constant for water vapor
  REAL(wp), PARAMETER :: cpv   = 1869.46_wp       !! [J/K/kg] specific heat at constant pressure
  REAL(wp), PARAMETER :: cvv   = cpv-rv           !! [J/K/kg] specific heat at constant volume
  REAL(wp), PARAMETER :: dv0   = 2.22e-5_wp       !! [m^2/s]  diff coeff of H2O vapor in dry air at tmelt
  !> - liquid / water
  REAL(wp), PARAMETER :: rhoh2o= 1000._wp         !> [kg/m3]  density of liquid water

 !REAL(wp), PARAMETER :: clw   = 4186.84_wp       !! [J/K/kg] specific heat of water
                                                  !!  see below 
  REAL(wp), PARAMETER ::  cv_i =  2000.0_wp
  !> - phase changes
  REAL(wp), PARAMETER :: alv   = 2.5008e6_wp      !> [J/kg]   latent heat for vaporisation
  REAL(wp), PARAMETER :: als   = 2.8345e6_wp      !! [J/kg]   latent heat for sublimation
  REAL(wp), PARAMETER :: alf   = als-alv          !! [J/kg]   latent heat for fusion
  REAL(wp), PARAMETER :: tmelt = 273.15_wp        !! [K]      melting temperature of ice/snow
  !
  !> Auxiliary constants
  REAL(wp), PARAMETER :: rdv   = rd/rv            !> [ ]
  REAL(wp), PARAMETER :: vtmpc1= rv/rd-1._wp      !! [ ]
  REAL(wp), PARAMETER :: vtmpc2= cpv/cpd-1._wp    !! [ ]
  REAL(wp), PARAMETER :: rcpv  = cpd/cpv-1._wp    !! [ ]
  REAL(wp), PARAMETER :: alvdcp= alv/cpd          !! [K]
  REAL(wp), PARAMETER :: alsdcp= als/cpd          !! [K]
  REAL(wp), PARAMETER :: rcpd  = 1._wp/cpd        !! [K*kg/J]
  REAL(wp), PARAMETER :: rcvd  = 1._wp/cvd        !! [K*kg/J]
  REAL(wp), PARAMETER :: rcpl  =  3.1733_wp       !! cp_d / cp_l - 1
  !
  REAL(wp), PARAMETER :: clw   = (rcpl + 1.0_wp) * cpd !> specific heat capacity of liquid water
  REAL(wp), PARAMETER :: cv_v  = (rcpv + 1.0_wp) * cpd - rv
  !
  REAL(wp), PARAMETER :: O_m_rdv  = 1._wp-rd/rv   !> [ ]
  REAL(wp), PARAMETER :: rd_o_cpd = rd/cpd        !! [ ]
  REAL(wp), PARAMETER :: cvd_o_rd = cvd/rd        !! [ ]
  !
  REAL(wp), PARAMETER :: p0ref     = 100000.0_wp   !> [Pa]  reference pressure for Exner function

  !> Variables for computing cloud cover in RH scheme
  REAL(wp), PARAMETER ::  uc1  = 0.8_wp
  REAL(wp), PARAMETER ::  ucl  = 1.00_wp

  !> vertical profile parameters (vpp) of CH4 and N2O
  REAL(wp), PARAMETER :: vpp_ch4(3) = (/1.25e-01_wp,  683.0_wp, -1.43_wp/)
  REAL(wp), PARAMETER :: vpp_n2o(3) = (/1.20e-02_wp, 1395.0_wp, -1.43_wp/)

  !> constants for radiation module
  REAL(wp), PARAMETER :: zemiss_def = 0.996_wp  !> lw sfc default emissivity factor


!------------below are parameters for ocean model---------------
  ! coefficients in linear EOS
  REAL(wp), PARAMETER :: a_T = 2.55E-04_wp,     & ! thermal expansion coefficient (kg/m3/K)
    &                    b_S = 7.64E-01_wp         ! haline contraction coefficient (kg/m3/psu)

  ! density reference values, to be constant in Boussinesq ocean models
  REAL(wp), PARAMETER :: rho_ref = 1025.022_wp         ! reference density [kg/m^3]
  REAL(wp), PARAMETER :: rho_inv = 0.0009755881663_wp  ! inverse reference density [m^3/kg]
  REAL(wp), PARAMETER :: sal_ref = 35.0_wp             ! reference salinity [psu]

  REAL(wp), PARAMETER :: SItodBar = 1.0e-4_wp          !Conversion from pressure [p] to pressure [bar]
                                                       !used in ocean thermodynamics
  REAL(wp),PARAMETER :: sfc_press_pascal = 101300.0_wp
  REAL(wp),PARAMETER :: sfc_press_bar    = 101300.0_wp*SItodBar

  REAL(wp), PARAMETER :: p0sl_bg         = 101325._wp  ! [Pa]     sea level pressure

!----------below are parameters for sea-ice model---------------
  REAL(wp), PARAMETER ::           &
    ks           = 0.31_wp,         & ! heat conductivity snow     [J  / (m s K)]
    ki           = 2.1656_wp,       & ! heat conductivity ice      [J  / (m s K)]   
    rhoi         = 917.0_wp,        & ! density of sea ice         [kg / m�]
    rhos         = 300.0_wp,        & ! density of snow            [kg / m�]
    ci           = 2106.0_wp,       & ! Heat capacity of ice       [J / (kg K)]
    Tf           = -1.80_wp,        & ! Temperature ice bottom     [C]
    Sice         = 5.0_wp,          & ! Sea-ice bulk salinity      [ppt]
    mu           = 0.054_wp,        & ! Constant in linear freezing-
                                      ! point relationship         [�C/ppt]
    muS          = mu*Sice,         & ! = - (sea-ice liquidus 
                                      ! (aka melting) temperature) [�C]
!   muS          = -(-0.0575 + 1.710523E-3*Sqrt(Sice) - 2.154996E-4*Sice) * Sice
    albs         = 0.75_wp,         & ! Albedo of snow (not melting)
    albsm        = 0.65_wp,         & ! Albedo of snow (melting)    
    albi         = 0.66_wp,         & ! Albedo of ice (not melting)
    albim        = 0.64_wp,         & ! Albedo of ice (melting)    
    albedoW      = 0.07_wp,         & ! albedo of the ocean 
    !I_0          = 0.3             ! Ice-surface penetrating shortwave fraction
    I_0          = 0.17_wp            ! Ice-surface penetrating shortwave fraction

!------------------------------------------------------------


END MODULE mo_physical_constants
