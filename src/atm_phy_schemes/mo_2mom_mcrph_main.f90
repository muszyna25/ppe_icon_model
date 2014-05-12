!===============================================================================!
!
! Two-moment bulk microphysics by Axel Seifert, Klaus Beheng and Uli Blahak
!
! Description:
! Provides various modules and subroutines for two-moment bulk microphysics
!
! Current Code Owner: Axel Seifert, DWD
!                     axel.seifert@dwd.de
!
! Language: Fortran 2003
!
! Some code standards or recommendations, at least:
!
! - All changes that potentially change the results need to 
!   be approved by AS and UB
! - All new variables/subroutines should be named in English
! - In the future also comments should be written in English, 
!   but temporary use of some German is OK, too.
! - Length of names of subroutines should be <= 20
! - Length of names of variables should be <= 15
! - Length of lines has to be < 120 including comments,
!   recommended is <= 100 for code without comments. 
! - Temporary modifications for experiments should be marked by
!
!     AS_YYYYMMDD>
!         ! Change for / Bugfix ....
!     <AS_YYYYMMDD
!
!   until they are validated as improvements and made official
!   (with AS, or whatever, being the initials of the guy doing the stuff
!   and YYYYMMDD=year/month/day).
!
!===============================================================================!
! Re-write for ICON 04/2014 by AS:
! Some general notes:
! - This version may need an up-to-date compiler due to some Fortran2003 features
! - Adapted physical constants to ICON
! - Atlas-type fall speed of rain has been changed to SBB2014, GMD
! Some notes on optimization (tests on thunder Thunder):
! - Small penalty for the particle%meanmass etc. functions, but compensated
!   by the exp(b*log(x)) instead of the original power law.
! - Increased q_crit from 1e-9 to 1e-7 for efficiency. Looks ok for WK-test,
!   but has to be tested in a real-case setup with stratiform and cirrus clouds.
! - Replaced some more power laws by exp(a*log(x)), e.g.,
!   in graupel_hail_conv_wet_gamlook()
! - Replaced almost all ()**0.5 by sqrt(), including the **m_f in ventilation 
!   coefficients
! - Clipping in rain_freeze is necessary, but removed everywhere else
!===============================================================================!
! To Do:
! - Check conservation of water mass
! - Further optimization might be possible in rain_freeze.
!===============================================================================!
! Further plans (physics) including HDCP2 project:
! - Implement prognostic cloud droplet number
! - Implement budget equations for CCN and IN
! - Implement new collision rate parameterizations of SBB2014
! - Implement improved height dependency of terminal fall velocity similar as 
!   used in COSMO two-moment code (but may be quite expensive).
! - Write a version with three or four different ice particle species 
!   (hom, het, frz, and splinters from ice multiplication)
! - New melting with prognostic melt water
!===============================================================================!
! Small stuff:
! - Increase alpha_spacefilling?
! - Are the minor differences in the sticking efficiencies important?
!===============================================================================!
! Further plans (restructuring and numerics):
! - Better understand performance issues of semi-implicit solver
! - Introduce logicals llqi_crit=(qi>q_crit), and llqi_zero = (qi>0.0), etc. 
!   which are calculated once in the driver
! - Include pointers for q and n in particle type and pass by argument
!   instead of global pointers in module (q_cloud, n_cloud, etc.)
!===============================================================================!

MODULE mo_2mom_mcrph_main

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message, txt => message_text
  USE mo_math_constants,     ONLY: pi, pi4 => pi_4
  USE mo_physical_constants, ONLY: &
       & R_l   => rd,     & ! gas constant of dry air (luft)
       & R_d   => rv,     & ! gas constant of water vapor (dampf)
       & cp    => cpd,    & ! specific heat capacity of air at constant pressure
       & c_w   => clw,    & ! specific heat capacity of water
       & L_wd  => alv,    & ! specific heat of vaporization (wd: wasser->dampf)
       & L_ed  => als,    & ! specific heat of sublimation (ed: eis->dampf)
       & L_ew  => alf,    & ! specific heat of fusion (ew: eis->wasser)
       & T_3   => tmelt,  & ! melting temperature of ice
       & rho_w => rhoh2o, & ! density of liquid water
       & nu_l  => con_m,  & ! kinematic viscosity of air
       & D_v   => dv0,    & ! diffusivity of water vapor in air at 0 C
       & K_t   => con0_h, & ! heat conductivity of air
       & N_avo => avo,    & ! Avogadro number [1/mol]
       & k_b   => ak,     & ! Boltzmann constant [J/K]
       & grav               ! acceleration due to Earth's gravity
  USE mo_satad, ONLY:     &
       & e_ws  => sat_pres_water,  & ! saturation pressure over liquid water
       & e_es  => sat_pres_ice       ! saturation pressure over ice

  USE mo_2mom_mcrph_util, ONLY: &
       & gfct,                       &  ! Gamma function (becomes intrinsic in Fortran2008)
       & gamlookuptable,             &  ! For look-up table of incomplete Gamma function
       & nlookup, nlookuphr_dummy,   &  !   array size of table
       & incgfct_lower_lookupcreate, &  !   create table
       & incgfct_lower_lookup,       &  !   interpolation in table, lower incomplete Gamma function
       & incgfct_upper_lookup,       &  !   interpolation in talbe, upper incomplete Gamma function
       & dmin_wg_gr_ltab_equi,       &  ! For look-up table of wet growth diameter
       & ltabdminwgg, lookupt_4D        !

  IMPLICIT NONE 

  PUBLIC

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_main'

#ifndef __SX__

  ! start and end indices for 2D slices
  INTEGER :: istart, iend, kstart, kend

  ! time step within two-moment scheme
  REAL(wp) :: dt

  ! switches for ice scheme, ice nucleation, drop activation and autoconversion 
  INTEGER  :: ice_typ, nuc_i_typ, nuc_c_typ, auto_typ

  ! Pointers for arrays which are global in this module, but allocated
  ! on the stack. Pointers are set in mo_2mom_mcrph_driver.
  real(wp), pointer, dimension(:,:) ::  &
       & w_p, p_p, t_p, rho_p, qv,                          &
       & rrho_04, rrho_c,                                   &
       & q_cloud, q_ice, q_rain, q_graupel, q_snow, q_hail, &
       & n_cloud, n_ice, n_rain, n_graupel, n_snow, n_hail
!$omp threadprivate (w_p)
!$omp threadprivate (p_p)
!$omp threadprivate (t_p)
!$omp threadprivate (rho_p)
!$omp threadprivate (qv)
!$omp threadprivate (q_cloud)
!$omp threadprivate (n_cloud)
!$omp threadprivate (q_rain)
!$omp threadprivate (n_rain)
!$omp threadprivate (q_ice)
!$omp threadprivate (n_ice)
!$omp threadprivate (q_snow)
!$omp threadprivate (n_snow)
!$omp threadprivate (q_graupel)
!$omp threadprivate (n_graupel)
!$omp threadprivate (q_hail)
!$omp threadprivate (n_hail)

  real(wp) :: qnc_const

  ! Derived type for hydrometeor species
  TYPE PARTICLE
    CHARACTER(20) :: name       !..name of particle class
    REAL(wp)      :: nu         !..first shape parameter of size distribution
    REAL(wp)      :: mu         !..2nd shape parameter
    REAL(wp)      :: x_max      !..max mean particle mass
    REAL(wp)      :: x_min      !..min mean particle mass
    REAL(wp)      :: a_geo      !..pre-factor in diameter-mass relation
    REAL(wp)      :: b_geo      !..exponent in diameter-mass relation
    REAL(wp)      :: a_vel      !..pre-factor in power law fall speed   (for now all particles have a power law fall
    REAL(wp)      :: b_vel      !..exponent in power law fall speed      speed, some have also an Atlas-type relation)
    REAL(wp)      :: a_ven      !..first parameter in ventilation coefficient
    REAL(wp)      :: b_ven      !..2nd parameter in ventilation coefficient 
    REAL(wp)      :: cap        !..coefficient for capacity of particle
    REAL(wp)      :: vsedi_max  !..max bulk sedimentation velocity
    REAL(wp)      :: vsedi_min  !..min bulk sedimentation velocity
  CONTAINS
    PROCEDURE     :: meanmass => particle_meanmass
    PROCEDURE     :: diameter => particle_diameter
    PROCEDURE     :: velocity => particle_velocity
  END TYPE PARTICLE

  ! .. for spherical particles we need to store the coefficients for the
  !    power law bulk sedimentation velocity
  TYPE particle_sphere
    REAL(wp)      :: coeff_alfa_n 
    REAL(wp)      :: coeff_alfa_q
    REAL(wp)      :: coeff_lambda
  END TYPE particle_sphere

  ! .. non-spherical particles have an Atlas-type terminal fall velocity relation
  TYPE particle_nonsphere
    REAL(wp)      :: alfa   !..1st parameter in Atlas-type fall speed
    REAL(wp)      :: beta   !..2nd parameter in Atlas-type fall speed
    REAL(wp)      :: gama   !..3rd parameter in Atlas-type fall speed
  END TYPE particle_nonsphere

  ! .. raindrops have an Atlas-type terminal fall velocity relation
  !    and a mu-D-relation which is used in sedimentation and evaporation
  !    (see Seifert 2008, J. Atmos. Sci.)
  TYPE particle_rain
    REAL(wp)      :: alfa   !..1st parameter in Atlas-type fall speed
    REAL(wp)      :: beta   !..2nd parameter in Atlas-type fall speed
    REAL(wp)      :: gama   !..3rd parameter in Atlas-type fall speed
    REAL(wp)      :: cmu0   !..Parameters for mu-D-relation of rain
    REAL(wp)      :: cmu1   !     max of left branch
    REAL(wp)      :: cmu2   !     max of right branch
    REAL(wp)      :: cmu3   !     min value of relation
    REAL(wp)      :: cmu4   !     location of min value = breakup equilibrium diameter
    INTEGER       :: cmu5   !     exponent
  CONTAINS
    PROCEDURE :: mue_Dm_relation => rain_mue_Dm_relation
  END TYPE particle_rain

  ! These are the fundamental particle variables for the two-moment scheme
  ! These pointers will be specified from the list of pre-defined particles 
  TYPE(particle), POINTER        :: cloud, ice, snow, graupel, hail, rain
  TYPE(particle_sphere)          :: ice_coeffs, snow_coeffs, graupel_coeffs, hail_coeffs
  TYPE(particle_rain),   POINTER :: rain_coeffs

  ! Physical parameters and coefficients which occur only in the two-moment scheme

  ! ... some physical parameters not found in ICON
  REAL(wp), PARAMETER :: rho0    = 1.225_wp     !..Norm-Luftdichte
  REAL(wp), PARAMETER :: T_f     = 233.0_wp     !..Bei T < T_f kein Fl.wasser
  REAL(wp), PARAMETER :: rho_ice = 916.7_wp     !..Materialdichte von Eis

  ! .. some cloud physics parameters
  REAL(wp), PARAMETER :: N_sc = 0.710_wp        !..Schmidt-Zahl (PK, S.541)
  REAL(wp), PARAMETER :: n_f  = 0.333_wp        !..Exponent von N_sc im Vent-koeff. (PK, S.541)
  REAL(wp), PARAMETER :: m_f  = 0.500_wp        !..Exponent von N_re im Vent-koeff. (PK, S.541)
                                                !  (replaced by sqrt almost everywhere)

  ! .. for old saturation pressure relations (keep this for some time for testing)
  REAL(wp), PARAMETER :: A_e  = 2.18745584e1_wp !..Konst. Saettigungsdamppfdruck - Eis
  REAL(wp), PARAMETER :: A_w  = 1.72693882e1_wp !..Konst. Saettigungsdamppfdruck - Wasser
  REAL(wp), PARAMETER :: B_e  = 7.66000000e0_wp !..Konst. Saettigungsdamppfdruck - Eis
  REAL(wp), PARAMETER :: B_w  = 3.58600000e1_wp !..Konst. Saettigungsdamppfdruck - Wasser
  REAL(wp), PARAMETER :: e_3  = 6.10780000e2_wp !..Saettigungsdamppfdruck bei T = T_3

  ! .. Hallet-Mossop ice multiplication
  REAL(wp), PARAMETER ::           &
       &    C_mult     = 3.5e8_wp, &    !..Koeff. fuer Splintering
       &    T_mult_min = 265.0_wp, &    !..Minimale Temp. Splintering
       &    T_mult_max = 270.0_wp, &    !..Maximale Temp. Splintering
       &    T_mult_opt = 268.0_wp       !..Optimale Temp. Splintering

  ! .. exponents for simple height dependency of terminal fall velocity
  REAL(wp), PARAMETER :: rho_vel    = 0.4e0_wp    !..exponent for density correction 
  REAL(wp), PARAMETER :: rho_vel_c  = 1.0e0_wp    !..for cloud droplets

  ! .. Phillips et al. ice nucleation scheme, see ice_nucleation_homhet() for more details 
  REAL(wp) ::                         & 
       &    na_dust    = 162.e3_wp,   & ! initial number density of dust [1/m³], Phillips08 value 162e3
       &    na_soot    =  15.e6_wp,   & ! initial number density of soot [1/m³], Phillips08 value 15e6
       &    na_orga    = 177.e6_wp,   & ! initial number density of organics [1/m3], Phillips08 value 177e6
       &    ni_het_max = 100.0e3_wp,  & ! max number of IN between 1-10 per liter, i.e. 1d3-10d3 
       &    ni_hom_max = 5000.0e3_wp    ! number of liquid aerosols between 100-5000 per liter

  INTEGER, PARAMETER ::               & ! Look-up table for Phillips et al. nucleation
       &    ttmax  = 30,              & ! sets limit for temperature in look-up table
       &    ssmax  = 60,              & ! sets limit for ice supersaturation in look-up table
       &    ttstep = 2,               & ! increment for temperature in look-up table
       &    ssstep = 1                  ! increment for ice supersaturation in look-up table

  REAL(wp), DIMENSION(0:100,0:100), SAVE :: &
       &    afrac_dust, &  ! look-up table of activated fraction of dust particles acting as ice nuclei
       &    afrac_soot, &  ! ... of soot particles
       &    afrac_orga     ! ... of organic material

  INCLUDE 'phillips_nucleation_2010.incf'

  ! Size thresholds for partioning of freezing rain in the hail scheme:
  ! Raindrops smaller than D_rainfrz_ig freeze into cloud ice, 
  ! drops between D_rainfrz_ig and D_rainfrz_gh freeze to graupel, and the
  ! largest raindrop freeze directly to hail.
  REAL(wp), PARAMETER ::               &
       &    D_rainfrz_ig = 0.50e-3_wp, & ! rain --> ice oder graupel
       &    D_rainfrz_gh = 1.25e-3_wp    ! rain --> graupel oder hail

  ! Various parameters for collision and conversion rates
  REAL(wp), PARAMETER ::                  &
       &    ecoll_ic     = 0.80_wp,       &  !..max. eff. for ice_cloud_riming
       &    ecoll_sc     = 0.80_wp,       &  !..max. eff. for snow_cloud_riming
       &    ecoll_gc     = 1.00_wp,       &  !..max. eff. for graupel_cloud_riming
       &    ecoll_hc     = 1.00_wp,       &  !..max. eff. for hail_cloud_riming
       &    ecoll_min    = 0.01_wp,       &  !..min. eff. for graupel_cloud, ice_cloud and snow_cloud
       &    ecoll_gg     = 0.10_wp,       &  !..collision efficiency for graupel selfcollection
       &    ecoll_gg_wet = 0.40_wp,       &  !    in case of wet graupel
       &    alpha_spacefilling = 0.01_wp, &  !..Raumerfuellungskoeff (max. 0.68)
       &    ice_s_vel  = 0.00_wp,         &  !..dispersion of fall velocity for collection kernel
       &    snow_s_vel = 0.25_wp,         &  !..dispersion of fall velocity (see SB2006, Eqs 60-63)
       &    r_shedding = 500.0e-6_wp,     &  !..mean radius of shedding droplets
       &    T_shed     = 263.2_wp            !..crit temperature for sheddding

  ! Even more parameters for collision and conversion rates
  REAL(wp), PARAMETER :: &
       &    q_crit_ii = 1.000e-6_wp, & ! q-threshold for ice_selfcollection 
       &    D_crit_ii = 100.0e-6_wp, & ! D-threshold for ice_selfcollection
       &    D_conv_ii = 75.00e-6_wp, & ! D-threshold for conversion in ice_selfcollection
       &    q_crit_ic = 1.000e-5_wp, & ! q-threshold for ice_cloud_riming
       &    D_crit_ic = 150.0e-6_wp, & ! D-threshold for ice_cloud_riming
       &    q_crit_ir = 1.000e-5_wp, & ! q-threshold for ice_rain_riming
       &    D_crit_ir = 100.0e-6_wp, & ! D-threshold for ice_rain_riming
       &    q_crit_sc = 1.000e-5_wp, & ! q-threshold for snow_cloud_riming
       &    D_crit_sc = 150.0e-6_wp, & ! D-threshold for snow_cloud_riming
       &    q_crit_sr = 1.000e-5_wp, & ! q-threshold for snow_rain_riming
       &    D_crit_sr = 100.0e-6_wp, & ! D-threshold for snow_rain_riming
       &    q_crit_gc = 1.000e-6_wp, & ! q-threshold for graupel_cloud_riming
       &    D_crit_gc = 100.0e-6_wp, & ! D-threshold for graupel_cloud_riming
       &    q_crit_hc = 1.000e-6_wp, & ! q-threshold for hail_cloud_riming
       &    D_crit_hc = 100.0e-6_wp, & ! D-threshold for hail_cloud_riming
       &    q_crit_fr = 1.000e-6_wp, & ! q-threshold for rain_freeze
       &    q_crit_c  = 1.000e-6_wp, & ! q-threshold for cloud water
       &    q_crit    = 1.000e-7_wp, & ! q-threshold elsewhere 1e-4 g/m3 = 0.1 mg/m3
       &    D_conv_sg = 200.0e-6_wp, & ! D-threshold for conversion of snow to graupel
       &    D_conv_ig = 200.0e-6_wp, & ! D-threshold for conversion of ice to graupel
       &    x_conv    = 0.100e-9_wp, & ! minimum mass of conversion due to riming
       &    D_shed_g  = 3.000e-3_wp, & ! D-threshold for graupel shedding
       &    D_shed_h  = 5.000e-3_wp, & ! D-threshold for hail shedding
       &    D_crit_c  = 10.00e-6_wp, & ! D-threshold for cloud drop collection efficiency
       &    D_coll_c  = 40.00e-6_wp    ! upper bound for diameter in collision efficiency

  REAL(wp), PARAMETER ::           &
       &    T_nuc     = 268.15_wp, & ! lower temperature threshold for ice nucleation, -5 C
       &    T_freeze  = 273.15_wp    ! lower temperature threshold for raindrop freezing

  ! choice of mu-D relation for rain, default is mu_Dm_rain_typ = 1
  INTEGER, PARAMETER     :: mu_Dm_rain_typ = 1     ! see init_twomoment() for possible choices

  ! Parameter for evaporation of rain, determines change of n_rain during evaporation
  REAL(wp) :: rain_gfak   ! this is set in init_twomoment depending on mu_Dm_rain_typ

  ! debug switch
  LOGICAL, PARAMETER     :: isdebug = .false.   ! use only when really desperate
                                     
  ! some cloud microphysical switches
  LOGICAL, PARAMETER     :: ice_multiplication = .TRUE.  ! default is .true.
  LOGICAL, PARAMETER     :: enhanced_melting   = .TRUE.  ! default is .true.

  ! switches for shedding
  ! UB: Das angefrorene Wasser wird bei T > T_shed, 
  !     ggf. nach enhanced-melting, wieder abgeworfen und zu Regen.
  !     Default is .false. for both parameters, because this simple bulk approach
  !     would limit the growth of graupel and hail too much.
  LOGICAL, PARAMETER     :: graupel_shedding = .FALSE.
  LOGICAL, PARAMETER     :: hail_shedding    = .FALSE.

  !..Jason's mu-Dm-relation for snow (see also Milbrandt & Yau 2005), 
  !  currently not used, because it shows no significant improvement
  LOGICAL,  PARAMETER :: use_mu_Dm_snow_sedi = .FALSE.
  REAL(wp), PARAMETER :: snow_cmu1 = 4.5_wp      ! a 
  REAL(wp), PARAMETER :: snow_cmu2 = 0.5e+3_wp   ! b 
  REAL(wp), PARAMETER :: snow_cmu3 = 5.0e-3_wp   ! Dnue
  REAL(wp), PARAMETER :: snow_cmu4 = 5.5_wp      ! a
  INTEGER,  PARAMETER :: snow_cmu5 = 1           ! exponent

  !..Pre-defined particle types
  TYPE(particle), TARGET :: graupelhail_cosmo5 = particle( & ! graupelhail2test4
       &                                 'graupelhail_cosmo5' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 5.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 1.00d-09, & !.x_min..minimale Teilchenmasse 
       &                                 1.42d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.314000, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 86.89371, & !.a_vel..Koeff. Fallgesetz
       &                                 0.268325, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.00,     & !.cap....Koeff. Kapazitaet
       &                                 30.0,     & !.vsedi_max
       &                                 0.1 )       !.vsedi_min
 
  TYPE(particle), TARGET :: hail_cosmo5 = particle( & ! hailULItest
       &                                 'hail_cosmo5' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 5.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-9, & !.x_min..minimale Teilchenmasse 
       &                                 0.1366 , & !.a_geo..Koeff. Geometrie
       &                                 0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                                 39.3    , & !.a_vel..Koeff. Fallgesetz
       &                                 0.166667, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.00,     & !.cap....Koeff. Kapazitaet
       &                                 30.0,     & !.vsedi_max
       &                                 0.1 )       !.vsedi_min

  TYPE(particle), TARGET :: cloud_cosmo5 = PARTICLE( &
       &                               'cloud_cosmo5',  & !.name...Bezeichnung der Partikelklasse
       &                                 0.0,      & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 2.60d-10, & !.x_max..maximale Teilchenmasse D=80e-6m
       &                                 4.20d-15, & !.x_min..minimale Teilchenmasse D=2.e-6m
       &                                 1.24d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                                 3.75d+05, & !.a_vel..Koeff. Fallgesetz
       &                                 0.666667, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.00,     & !.cap....Koeff. Kapazitaet
       &                                 1.0,      & !.vsedi_max
       &                                 0.0)        !.vsedi_min

  TYPE(particle), TARGET :: cloud_nue1mue1 = PARTICLE( &
       &                               'cloud_nue1mue1',  & !.name...Bezeichnung der Partikelklasse
       &                               1.000000, & !.nu.....Breiteparameter der Verteil.
       &                               1.000000, & !.mu.....Exp.-parameter der Verteil.
       &                               2.60d-10, & !.x_max..maximale Teilchenmasse D=80e-6m
       &                               4.20d-15, & !.x_min..minimale Teilchenmasse D=2.e-6m
       &                               1.24d-01, & !.a_geo..Koeff. Geometrie
       &                               0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                               3.75d+05, & !.a_vel..Koeff. Fallgesetz
       &                               0.666667, & !.b_vel..Koeff. Fallgesetz
       &                               0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                               0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.00,     & !.cap....Koeff. Kapazitaet
       &                                 1.0,      & !.vsedi_max
       &                                 0.0)        !.vsedi_min

  TYPE(particle), TARGET :: ice_cosmo5 =  particle( & ! iceCRY2test
       &                              'ice_cosmo5', & !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu...e..Breiteparameter der Verteil.
       &                              0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                              1.00d-06, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.00d-12, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              0.835000, & !.a_geo..Koeff. Geometrie
       &                              0.390000, & !.b_geo..Koeff. Geometrie
       &                              2.77d+01, & !.a_vel..Koeff. Fallgesetz 
       &                              0.215790, & !.b_vel..Koeff. Fallgesetz = 0.41/1.9
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0,      & !.cap....Koeff. Kapazitaet
       &                              3.0,      & !.vsedi_max
       &                              0.0 )       !.vsedi_min

  TYPE(particle), TARGET :: snow_cosmo5 = particle( & ! nach Andy Heymsfield (CRYSTAL-FACE)
       &                              'snow_cosmo5', & !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu.....Breiteparameter der Verteil.
       &                              0.500000, & !.mu.....Exp.-parameter der Verteil.
       &                              2.00d-05, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.00d-10, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              2.400000, & !.a_geo..Koeff. Geometrie
       &                              0.455000, & !.b_geo..Koeff. Geometrie
       &                              8.800000, & !.a_vel..Koeff. Fallgesetz 
       &                              0.150000, & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.00,     & !.cap....Koeff. Kapazitaet
       &                              3.0,      & !.vsedi_max
       &                              0.1 )       !.vsedi_min

  TYPE(particle), TARGET :: snowSBB =  particle( & ! 
       &                              'snowSBB',& !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu.....Breiteparameter der Verteil.
       &                              0.500000, & !.mu.....Exp.-parameter der Verteil.
       &                              2.00d-05, & !.x_max..maximale Teilchenmasse
       &                              1.00d-10, & !.x_min..minimale Teilchenmasse
       &                              5.130000, & !.a_geo..Koeff. Geometrie, x = 0.038*D**2
       &                              0.500000, & !.b_geo..Koeff. Geometrie = 1/2
       &                              8.294000, & !.a_vel..Koeff. Fallgesetz 
       &                              0.125000, & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.00,     & !.cap....Koeff. Kapazitaet
       &                              3.0,      & !.vsedi_max
       &                              0.1 )       !.vsedi_min

  TYPE(particle), TARGET :: rainULI = particle( & ! Blahak, v=v(x) gefittet 6.9.2005
       &                              'rainULI', & !..name
       &                              0.000000,  & !..nu
       &                              0.333333,  & !..mu
       &                              3.00d-06,  & !..x_max
       &                              2.60d-10,  & !..x_min
       &                              1.24d-01,  & !..a_geo
       &                              0.333333,  & !..b_geo
       &                              114.0137,  & !..a_vel
       &                              0.234370,  & !..b_vel
       &                              0.780000,  & !..a_ven
       &                              0.308000,  & !..b_ven
       &                              2.00,      & !..cap....Koeff. Kapazitaet
       &                              20.0,      & !..vsedi_max
       &                              0.1)         !..vsedi_min

  TYPE(particle), TARGET :: rainSBB = particle( & 
       &                              'rainSBB', & !..name
       &                              0.000000,  & !..nu
       &                              0.333333,  & !..mu
       &                              3.00d-06,  & !..x_max
       &                              2.60d-10,  & !..x_min
       &                              1.24d-01,  & !..a_geo
       &                              0.333333,  & !..b_geo
       &                              114.0137,  & !..a_vel
       &                              0.234370,  & !..b_vel
       &                              0.780000,  & !..a_ven
       &                              0.308000,  & !..b_ven
       &                              2.000000,  & !..cap
       &                              2.000d+1,  & !..vsedi_max
       &                              0.100000 )   !..vsedi_min

  TYPE(particle_rain), TARGET :: rainSBBcoeffs = particle_rain( & 
       &                              9.292000,  & !..alfa
       &                              9.623000,  & !..beta
       &                              6.222d+2,  & !..gama
       &                              6.000000,  & !..cmu0
       &                              3.000d+1,  & !..cmu1
       &                              1.000d+3,  & !..cmu2
       &                              1.100d-3,  & !..cmu3 = D_br
       &                              1.000000,  & !..cmu4
       &                              2 )          !..cmu5

  REAL(wp), PARAMETER :: pi6 = pi/6.0_wp, pi8 = pi/8.0_wp ! more pieces of pi

CONTAINS

  !*******************************************************************************
  ! This subroutine has to be called once at the start of the model run by
  ! the main program. It properly sets the parameters for the different hydrometeor
  ! classes according to predefined parameter sets (see above).
  !*******************************************************************************

  SUBROUTINE init_2mom_scheme( cloud_type )
    IMPLICIT NONE
    INTEGER, INTENT(in), OPTIONAL :: cloud_type

    INTEGER :: wolke_typ = 2403 ! default cloud type

    REAL(wp), DIMENSION(1:1) :: x_r,q_c,vn_rain_min, vq_rain_min, vn_rain_max, vq_rain_max

    IF (PRESENT(cloud_type)) THEN
       wolke_typ = cloud_type
    END IF

    ice_typ   = wolke_typ/1000           ! (0) no ice, (1) no hail (2) with hail
    nuc_i_typ = MOD(wolke_typ/100,10)    ! choice of ice nucleation, see ice_nucleation_homhet()
    nuc_c_typ = MOD(wolke_typ/10,10)     ! choice of CCN assumptions, see cloud_nucleation()
    auto_typ  = MOD(wolke_typ,10)        ! choice of warm rain scheme, see clouds_twomoment()

    IF (cloud_type < 2000) THEN
      ! Without hail class:
      cloud   => cloud_nue1mue1
      rain    => rainSBB
      ice     => ice_cosmo5
      snow    => snowSBB
      graupel => graupelhail_cosmo5
      ! Dummy value for hail:
      hail    => hail_cosmo5
    ELSE
      ! Including hail class:
      cloud   => cloud_nue1mue1
      rain    => rainSBB
      ice     => ice_cosmo5
      snow    => snowSBB
      graupel => graupelhail_cosmo5
      hail    => hail_cosmo5
    END IF

    rain_coeffs => rainSBBcoeffs

    ! initialize bulk sedimentation velocities
    ! calculates coeff_alfa_n, coeff_alfa_q, and coeff_lambda
    call init_sedi_vel(ice,ice_coeffs)
    call init_sedi_vel(snow,snow_coeffs)
    call init_sedi_vel(graupel,graupel_coeffs)
    call init_sedi_vel(hail,hail_coeffs)

    ! other options for mue-D-relation of raindrops (for sensitivity studies)
    IF (mu_Dm_rain_typ.EQ.0) THEN
       !..constant mue value
       rain_coeffs%cmu0 = 0.0
       rain_coeffs%cmu1 = 0.0
       rain_coeffs%cmu2 = 1.0
       rain_coeffs%cmu3 = 1.0
       rain_coeffs%cmu4 = (rain%nu+1.0_wp)/rain%b_geo - 1.0_wp ! <-- this is the (constant) mue value 
       rain_coeffs%cmu5 = 1
       rain_gfak = -1.0  ! In this case gamma = 1 in rain_evaporation    
    ELSEIF (mu_Dm_rain_typ.EQ.1) THEN
       !..Axel's mu-Dm-relation for raindrops based on 1d-bin model
       !  (this is the default and the cmus are set in the particle constructor)
       rain_gfak = 1.0      
    ELSEIF (mu_Dm_rain_typ.EQ.2) THEN
       !..Modifikation of mu-Dm-relation for experiments with increased evaporation
       rain_coeffs%cmu0 = 11.0            ! instead of 6.0      
       rain_coeffs%cmu1 = 30.0            ! 
       rain_coeffs%cmu2 = 1.00d+3         ! 
       rain_coeffs%cmu3 = 1.10d-3         ! 
       rain_coeffs%cmu4 = 4.0             ! instead of 1.0  
       rain_coeffs%cmu5 = 2               ! 
       rain_gfak = 0.5             ! instead of 1.0      
    ELSEIF (mu_Dm_rain_typ.EQ.3) THEN
      !..Jason Milbrandts mu-Dm-relation'
       rain_coeffs%cmu0 = 19.0           !
       rain_coeffs%cmu1 = 19.0           ! Jason Milbrandt's mu-Dm-relation for rain_coeffs
       rain_coeffs%cmu2 = 0.60d+3        ! (Milbrandt&Yau 2005, JAS, Table 1)
       rain_coeffs%cmu3 = 1.80d-3        !
       rain_coeffs%cmu4 = 17.0           !
       rain_coeffs%cmu5 = 1              !
       rain_gfak = -1.0  ! In this case gamma = 1 in rain_evaporation    
    ENDIF

    CALL message(TRIM(routine), "init_2mom_scheme: rain coeffs and sedi vel")
    WRITE (txt,'(2A)') "    name  = ",rain%name ; CALL message(routine,TRIM(txt))        
    WRITE(txt,'(A,D10.3)') "     alfa  = ",rain_coeffs%alfa ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     beta  = ",rain_coeffs%beta ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     gama  = ",rain_coeffs%gama ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     cmu0  = ",rain_coeffs%cmu0 ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     cmu1  = ",rain_coeffs%cmu1 ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     cmu2  = ",rain_coeffs%cmu2 ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     cmu3  = ",rain_coeffs%cmu3 ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     cmu4  = ",rain_coeffs%cmu4 ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,I10)')   "     cmu5  = ",rain_coeffs%cmu5 ; CALL message(routine,TRIM(txt))
    x_r = rain%x_min ; CALL sedi_vel_rain(rain,rain_coeffs,x_r,vn_rain_min,vq_rain_min,1,1)
    x_r = rain%x_max ; CALL sedi_vel_rain(rain,rain_coeffs,x_r,vn_rain_max,vq_rain_max,1,1)
    WRITE(txt,'(A)')       "    out-of-cloud: " ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vn_rain_min  = ",vn_rain_min ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vn_rain_max  = ",vn_rain_max ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vq_rain_min  = ",vq_rain_min ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vq_rain_max  = ",vq_rain_max ; CALL message(routine,TRIM(txt))
    q_c = 1e-3_wp
    x_r = rain%x_min ; CALL sedi_vel_rain(rain,rain_coeffs,x_r,vn_rain_min,vq_rain_min,1,1,q_c)
    x_r = rain%x_max ; CALL sedi_vel_rain(rain,rain_coeffs,x_r,vn_rain_max,vq_rain_max,1,1,q_c)
    WRITE(txt,'(A)')       "    in-cloud: " ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vn_rain_min  = ",vn_rain_min ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vn_rain_max  = ",vn_rain_max ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vq_rain_min  = ",vq_rain_min ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vq_rain_max  = ",vq_rain_max ; CALL message(routine,TRIM(txt))

  END SUBROUTINE init_2mom_scheme

  !*******************************************************************************
  ! Main subroutine of the two-moment microphysics                               
  !
  ! All individual processes are called in sequence and do their own time 
  ! integration for the mass and number densities, i.e., Marchuk-type operator 
  ! splitting. Temperature is only updated in the driver.
  !*******************************************************************************
  
  SUBROUTINE clouds_twomoment (nsize, nlev)

    INTEGER, INTENT (IN) :: nsize, nlev

    REAL(wp), DIMENSION(nsize,nlev) :: dep_rate_ice, dep_rate_snow

    INTEGER :: k

    IF (isdebug) CALL message(TRIM(routine),"clouds_twomoment start")

    dep_rate_ice(:,:)  = 0.0_wp
    dep_rate_snow(:,:) = 0.0_wp

!    IF (isdebug) WRITE (txt,*) 'CLOUDS: cloud_nucleation'
      IF (nuc_c_typ .EQ. 0) THEN
        IF (isdebug) WRITE (txt,*) '  ... force constant cloud droplet number conc.'
        DO k=kstart,kend
           n_cloud(:,k) = qnc_const
        END DO
     END IF
!      ELSEIF (nuc_c_typ < 6) THEN
!        IF (isdebug) WRITE (txt,*) '  ... according to SB2006'
!        CALL cloud_nucleation()
!      ELSE
!        IF (isdebug) WRITE (txt,*) '  ... look-up tables according to Segal& Khain'
!        CALL cloud_nucleation_SK()
!      END IF
!    END IF

    IF (nuc_c_typ.ne.0) THEN
       DO k=kstart,kend
          n_cloud(:,k) = MAX(n_cloud(:,k), q_cloud(:,k) / cloud%x_max)
          n_cloud(:,k) = MIN(n_cloud(:,k), q_cloud(:,k) / cloud%x_min)
       END DO
    END IF

    IF (ice_typ > 0) THEN

       ! homogeneous and heterogeneous ice nucleation
       CALL ice_nucleation_homhet()
       DO k=kstart,kend
          n_ice(:,k) = MIN(n_ice(:,k), q_ice(:,k)/ice%x_min)
          n_ice(:,k) = MAX(n_ice(:,k), q_ice(:,k)/ice%x_max)
       END DO

       ! homogeneous freezing of cloud droplets
       CALL cloud_freeze ()
       
       ! depositional growth of all ice particles
       ! ( store deposition rate of ice and snow for conversion calculation in 
       !   ice_riming and snow_riming )
       CALL vapor_dep_relaxation (nsize, nlev, dt, dep_rate_ice, dep_rate_snow)

       ! ice-ice collisions
       CALL ice_selfcollection ()
       CALL snow_selfcollection ()
       CALL snow_ice_collection ()
       CALL graupel_selfcollection ()
       CALL graupel_ice_collection ()
       CALL graupel_snow_collection ()
       
       IF (ice_typ > 1) THEN
          ! conversion of graupel to hail in wet growth regime
          CALL graupel_hail_conv_wet_gamlook ()
          ! hail collisions
          CALL hail_ice_collection ()    ! Important?
          CALL hail_snow_collection ()   ! Important?
          ! AS Question: If yes, then implement a hail_graupel_collection with low efficiency?
       END IF

       ! riming of ice with cloud droplets and rain drops, and conversion to graupel
       CALL ice_riming (nsize,nlev,dep_rate_ice)

       ! riming of snow with cloud droplets and rain drops, and conversion to graupel
       CALL snow_riming (nsize,nlev,dep_rate_snow)
       
       ! more riming
       IF (ice_typ > 1) THEN
          CALL hail_cloud_riming ()
          CALL hail_rain_riming ()
       END IF
       CALL graupel_cloud_riming ()
       CALL graupel_rain_riming ()

       ! freezing of rain and conversion to ice/graupel/hail
       CALL rain_freeze_gamlook ()

       ! melting
       CALL ice_melting ()
       CALL snow_melting ()
       CALL graupel_melting ()
       IF (ice_typ > 1) CALL hail_melting ()

       ! evaporation from melting ice particles
       CALL snow_evaporation ()
       CALL graupel_evaporation ()
       IF (ice_typ > 1) CALL hail_evaporation ()

    ENDIF

    ! warm rain processes 
    ! (using something other than SB is somewhat inconsistent and not recommended)
    IF (auto_typ == 1) THEN
       CALL autoconversionKB ()   ! Beheng (1994)
       CALL accretionKB ()
       CALL rain_selfcollectionSB ()
    ELSE IF (auto_typ == 2) THEN
       ! Khairoutdinov and Kogan (2000)
       ! (KK2000 originally assume a 25 micron size threshold)
       CALL autoconversionKK ()   
       CALL accretionKK ()
       CALL rain_selfcollectionSB ()
    ELSE IF (auto_typ == 3) THEN
       CALL autoconversionSB ()   ! Seifert and Beheng (2001)
       CALL accretionSB ()
       CALL rain_selfcollectionSB ()
    ENDIF

    ! evaporation of rain following Seifert (2008)
    CALL rain_evaporation ()

    ! size limits for all hydrometeors
    IF (nuc_c_typ > 0) THEN
       DO k=kstart,kend
          n_cloud(:,k) = MIN(n_cloud(:,k), q_cloud(:,k)/cloud%x_min)
          n_cloud(:,k) = MAX(n_cloud(:,k), q_cloud(:,k)/cloud%x_max)
          ! Hard upper limit for cloud number conc.
          n_cloud(:,k) = MIN(n_cloud(:,k), 5000d6)
       END DO
    end if
       DO k=kstart,kend
          n_rain(:,k) = MIN(n_rain(:,k), q_rain(:,k)/rain%x_min)
          n_rain(:,k) = MAX(n_rain(:,k), q_rain(:,k)/rain%x_max)
       END DO
    IF (ice_typ > 0) THEN
       DO k=kstart,kend
          n_ice(:,k) = MIN(n_ice(:,k), q_ice(:,k)/ice%x_min)
          n_ice(:,k) = MAX(n_ice(:,k), q_ice(:,k)/ice%x_max)
          n_snow(:,k) = MIN(n_snow(:,k), q_snow(:,k)/snow%x_min)
          n_snow(:,k) = MAX(n_snow(:,k), q_snow(:,k)/snow%x_max)
          n_graupel(:,k) = MIN(n_graupel(:,k), q_graupel(:,k)/graupel%x_min)
          n_graupel(:,k) = MAX(n_graupel(:,k), q_graupel(:,k)/graupel%x_max)
       END DO
    END IF
    IF (ice_typ > 1) THEN
       DO k=kstart,kend
          n_hail(:,k) = MIN(n_hail(:,k), q_hail(:,k)/hail%x_min)
          n_hail(:,k) = MAX(n_hail(:,k), q_hail(:,k)/hail%x_max)
       END DO
    END IF

    IF (isdebug) CALL message(TRIM(routine),"clouds_twomoment end")

  END SUBROUTINE clouds_twomoment

  !*******************************************************************************
  ! Functions and subroutines working on particle class
  !*******************************************************************************
  ! (1) CLASS procedures for particle class
  !*******************************************************************************

  ! mean mass with limiters, Eq. (94) of SB2006
  PURE FUNCTION particle_meanmass(this,q,n) RESULT(xmean)
    CLASS(particle), INTENT(in) :: this
    REAL(wp),        INTENT(in) :: q,n
    REAL(wp)                    :: xmean
    REAL(wp), PARAMETER         :: eps = 1e-20_wp

    xmean = MIN(MAX(q/(n+eps),this%x_min),this%x_max) 
    RETURN
  END FUNCTION particle_meanmass

  ! mass-diameter relation, power law, Eq. (32) of SB2006 
  PURE FUNCTION particle_diameter(this,x) RESULT(D)
    CLASS(particle), INTENT(in) :: this
    REAL(wp),        INTENT(in) :: x
    REAL(wp)                    :: D
    
    D = this%a_geo * exp(this%b_geo*log(x))    ! v = a_geo * x**b_geo
    RETURN
  END FUNCTION particle_diameter

  ! terminal fall velocity of particles, cf. Eq. (33) of SB2006
  PURE FUNCTION particle_velocity(this,x) RESULT(v)
    CLASS(particle), INTENT(in) :: this
    REAL(wp),        INTENT(in) :: x
    REAL(wp)                    :: v

    v = this%a_vel * exp(this%b_vel * log(x))  ! v = a_vel * x**b_vel 
    RETURN
  END FUNCTION particle_velocity

  ! mue-Dm relation of raindrops
  PURE FUNCTION rain_mue_dm_relation(this,D_m) result(mue)
    CLASS(particle_rain), intent(in) :: this
    real(wp),             intent(in) :: D_m
    real(wp)                         :: mue, delta

    delta = this%cmu2*(D_m-this%cmu3)
    if (D_m.le.this%cmu3) then    
       mue = this%cmu0*tanh((4.*delta)**2) + this%cmu4
    else
       mue = this%cmu1*tanh(delta**2) + this%cmu4
    endif
    return
  END FUNCTION rain_mue_dm_relation

  !*******************************************************************************
  ! (2) More functions working on particle class, these are not CLASS procedures)
  !*******************************************************************************

  ! bulk ventilation coefficient, Eq. (88) of SB2006
  REAL(wp) FUNCTION vent_coeff_a(parti,n) 
    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: n
    TYPE(particle), INTENT(IN) :: parti

    vent_coeff_a = parti%a_ven * gfct((parti%nu+n+parti%b_geo)/parti%mu)              &
         &                     / gfct((parti%nu+1.0)/parti%mu)                        & 
         &                   * ( gfct((parti%nu+1.0)/parti%mu)                        & 
         &                     / gfct((parti%nu+2.0)/parti%mu) )**(parti%b_geo+n-1.0) 
  END FUNCTION vent_coeff_a

  ! bulk ventilation coefficient, Eq. (89) of SB2006
  REAL(wp) FUNCTION vent_coeff_b(parti,n)
    IMPLICIT NONE
    INTEGER, INTENT(in)         :: n
    TYPE(particle), INTENT(in) :: parti

    REAL(wp), PARAMETER :: m_f = 0.500 ! see PK, S.541. Do not change.

    vent_coeff_b = parti%b_ven                                                  & 
         & * gfct((parti%nu+n+(m_f+1.0)*parti%b_geo+m_f*parti%b_vel)/parti%mu)  &
         &             / gfct((parti%nu+1.0)/parti%mu)                          & 
         &           * ( gfct((parti%nu+1.0)/parti%mu)                          &
         &             / gfct((parti%nu+2.0)/parti%mu)                          &
         &             )**((m_f+1.0)*parti%b_geo+m_f*parti%b_vel+n-1.0)
  END FUNCTION vent_coeff_b

  ! complete mass moment of particle size distribution, Eq (82) of SB2006
  REAL(wp) FUNCTION moment_gamma(p,n) 
    IMPLICIT NONE
    INTEGER, INTENT(in)           :: n
    TYPE(particle), INTENT(in)   :: p

    moment_gamma  = gfct((n+p%nu+1.0)/p%mu) / gfct((p%nu+1.0)/p%mu)        &
         &      * ( gfct((  p%nu+1.0)/p%mu) / gfct((p%nu+2.0)/p%mu) )**n
  END FUNCTION moment_gamma

  ! fractional mass moment of particle size distribution, i.e., this is 
  ! Eq (82) of SB2006 with a non-integer exponent fexp
  REAL(wp) FUNCTION fracmoment_gamma(p,fexp) 
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: fexp
    TYPE(particle), INTENT(in)   :: p

    fracmoment_gamma  = gfct((fexp+p%nu+1.0)/p%mu) / gfct((p%nu+1.0)/p%mu)        &
         &          * ( gfct((    p%nu+1.0)/p%mu) / gfct((p%nu+2.0)/p%mu) )**fexp
  END FUNCTION fracmoment_gamma

  ! coefficient for slope of PSD, i.e., for lambda in Eq. (80) of SB2006
  REAL(wp) FUNCTION lambda_gamma(p,x) 
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: x
    TYPE(particle), INTENT(in)   :: p

    lambda_gamma  = ( gfct((p%nu+1.0)/p%mu) / gfct((p%nu+2.0)/p%mu) * x)**(-p%mu)
  END FUNCTION lambda_gamma

  ! coefficient for general collision integral, Eq. (90) of SB2006
  REAL(wp) FUNCTION coll_delta(p1,n)
    IMPLICIT NONE
    TYPE(particle), INTENT(in) :: p1
    INTEGER, INTENT(in)         :: n

    coll_delta = gfct((2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)         &
         &                  / gfct((p1%nu+1.0  )/p1%mu)         &
         &        * gfct((p1%nu+1.0)/p1%mu)**(2.0*p1%b_geo+n)   &
         &        / gfct((p1%nu+2.0)/p1%mu)**(2.0*p1%b_geo+n)
    RETURN
  END FUNCTION coll_delta

  ! wrapper for coll_delta (unnecessary and unused argument p2, but do not remove this)
  REAL(wp) FUNCTION coll_delta_11(p1,p2,n)
    TYPE(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n
    coll_delta_11 = coll_delta(p1,n)
    RETURN
  END FUNCTION coll_delta_11

  ! wrapper for coll_delta (unnecessary and unused argument p2, but do not remove this)
  REAL(wp) FUNCTION coll_delta_22(p1,p2,n)
    TYPE(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n
    coll_delta_22 = coll_delta(p2,n)
    RETURN
  END FUNCTION coll_delta_22

  ! coefficient for general collision integral, Eq. (91) of SB2006
  REAL(wp) FUNCTION coll_delta_12(p1,p2,n)
    TYPE(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n

    coll_delta_12 = 2.0 * gfct((p1%b_geo+p1%nu+1.0)/p1%mu)               &
         &                       / gfct((p1%nu+1.0)/p1%mu)               &
         &                * gfct((p1%nu+1.0)/p1%mu)**(p1%b_geo)          &
         &                / gfct((p1%nu+2.0)/p1%mu)**(p1%b_geo)          &
         &              * gfct((p2%b_geo+p2%nu+1.0+n)/p2%mu)             &
         &                        /gfct((p2%nu+1.0  )/p2%mu)             &
         &                * gfct((p2%nu+1.0)/p2%mu)**(p2%b_geo+n)        &
         &                / gfct((p2%nu+2.0)/p2%mu)**(p2%b_geo+n)
    RETURN
  END FUNCTION coll_delta_12

  ! coefficient for general collision integral, Eq. (92) of SB2006
  REAL(wp) FUNCTION coll_theta(p1,n)
    TYPE(particle), INTENT(in) :: p1
    INTEGER, INTENT(in)         :: n

    coll_theta = gfct((2.0*p1%b_vel+2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)    &
         &                  / gfct((2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)    &
         &             * gfct((p1%nu+1.0)/p1%mu)**(2.0*p1%b_vel)        &
         &             / gfct((p1%nu+2.0)/p1%mu)**(2.0*p1%b_vel)
    RETURN
  END FUNCTION coll_theta

  ! wrapper for coll_theta (unnecessary and unused argument p2, but do not remove this)
  REAL(wp) FUNCTION coll_theta_11(p1,p2,n)
    TYPE(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n

    coll_theta_11 = coll_theta(p1,n)
    RETURN
  END FUNCTION coll_theta_11

  ! wrapper for coll_theta (unnecessary and unused argument p2, but do not remove this)
  REAL(wp) FUNCTION coll_theta_22(p1,p2,n)
    TYPE(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n

    coll_theta_22 = coll_theta(p2,n)
    RETURN
  END FUNCTION coll_theta_22

  ! coefficient for general collision integral, Eq. (93) of SB2006
  REAL(wp) FUNCTION coll_theta_12(p1,p2,n)
    TYPE(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n

    coll_theta_12 = 2.0 * gfct((p1%b_vel+2.0*p1%b_geo+p1%nu+1.0)/p1%mu)       &  
         &                       / gfct((2.0*p1%b_geo+p1%nu+1.0)/p1%mu)       &
         &                 * gfct((p1%nu+1.0)/p1%mu)**(p1%b_vel)              &
         &                 / gfct((p1%nu+2.0)/p1%mu)**(p1%b_vel)              &
         &              * gfct((p2%b_vel+2.0*p2%b_geo+p2%nu+1.0+n)/p2%mu)     &  
         &                       / gfct((2.0*p2%b_geo+p2%nu+1.0+n)/p2%mu)     &
         &                 * gfct((p2%nu+1.0)/p2%mu)**(p2%b_vel)              &
         &                 / gfct((p2%nu+2.0)/p2%mu)**(p2%b_vel)
    RETURN
  END FUNCTION coll_theta_12

  ! bulk sedimentation velocities
  subroutine sedi_vel_rain(this,thisCoeffs,x,vn,vq,its,ite,qc)
    TYPE(particle), intent(in)      :: this
    TYPE(particle_rain), intent(in) :: thisCoeffs
    integer,  intent(in)  :: its,ite
    real(wp), intent(in)  :: x(its:ite)
    real(wp), intent(in), optional  :: qc(its:ite)
    real(wp), intent(inout) :: vn(its:ite), vq(its:ite)

    integer  :: i
    real(wp) :: D_m,mue,D_p

    do i=its,ite
       D_m = this%diameter(x(i))
       IF (PRESENT(qc)) THEN
          if (qc(i) >= q_crit) THEN 
             mue = (this%nu+1.0_wp)/this%b_geo - 1.0_wp
          else
             mue = thisCoeffs%mue_Dm_relation(D_m)
          end if
       ELSE
          mue = thisCoeffs%mue_Dm_relation(D_m)
       END IF
       D_p = D_m * exp((-1./3.)*log((mue+3.)*(mue+2.)*(mue+1.)))          
       vn(i) = thisCoeffs%alfa - thisCoeffs%beta * exp(-(mue+1.)*log(1.0 + thisCoeffs%gama*D_p))
       vq(i) = thisCoeffs%alfa - thisCoeffs%beta * exp(-(mue+4.)*log(1.0 + thisCoeffs%gama*D_p))
    end do
  end subroutine sedi_vel_rain

  ! bulk sedimentation velocities
  subroutine sedi_vel_sphere(this,thisCoeffs,x,vn,vq,its,ite)
    TYPE(particle), intent(in)        :: this
    TYPE(particle_sphere), intent(in) :: thisCoeffs
    integer,  intent(in)  :: its,ite
    real(wp), intent(in)  :: x(its:ite)
    real(wp), intent(out) :: vn(its:ite), vq(its:ite)

    integer  :: i
    real(wp) :: lam,v_n,v_q

    do i=its,ite
       lam = exp(this%b_vel* log(thisCoeffs%coeff_lambda*x(i)))
       v_n = thisCoeffs%coeff_alfa_n * lam
       v_q = thisCoeffs%coeff_alfa_q * lam
       v_n = MAX(v_n,this%vsedi_min)
       v_q = MAX(v_q,this%vsedi_min)
       v_n = MIN(v_n,this%vsedi_max)
       v_q = MIN(v_q,this%vsedi_max)
       vn(i) = v_n
       vq(i) = v_q
    end do
  end subroutine sedi_vel_sphere

  ! initialize coefficients for bulk sedimentation velocity
  subroutine init_sedi_vel(this,thisCoeffs)
    TYPE(particle)        :: this
    TYPE(particle_sphere) :: thisCoeffs

    thisCoeffs%coeff_alfa_n = this%a_vel * gfct((this%nu+this%b_vel+1.0)/this%mu) / gfct((this%nu+1.0)/this%mu)
    thisCoeffs%coeff_alfa_q = this%a_vel * gfct((this%nu+this%b_vel+2.0)/this%mu) / gfct((this%nu+2.0)/this%mu)
    thisCoeffs%coeff_lambda = gfct((this%nu+1.0)/this%mu)/gfct((this%nu+2.0)/this%mu)

    CALL message(TRIM(routine), "init_sedi_vel: ")
    WRITE (txt,'(2A)') "    name  = ",this%name ; CALL message(routine,TRIM(txt))        
    WRITE (txt,'(A,D10.3)') "    c_lam = ",thisCoeffs%coeff_lambda ; CALL message(routine,TRIM(txt))        
    WRITE (txt,'(A,D10.3)') "    alf_n = ",thisCoeffs%coeff_alfa_n ; CALL message(routine,TRIM(txt))
    WRITE (txt,'(A,D10.3)') "    alf_q = ",thisCoeffs%coeff_alfa_q ; CALL message(routine,TRIM(txt))
    
  end subroutine init_sedi_vel

  ! currently not used
  FUNCTION D_average_factor (parti)
    ! UB: Faktor zur Berechnung des mittleren Durchmessers von verallg. gammaverteilten Hydrometeoren:
    !     gueltig fuer D = a_geo * x^b_geo
    !     Berechnung des mittleren Durchmessers: D_average = parti%b_geo * D_average_factor * (q/qn)**parti%b_geo
    REAL(wp) :: D_average_factor
    TYPE(particle), INTENT(in) :: parti
    
    D_average_factor = &
         ( gfct( (parti%b_geo+parti%nu+1.0)/parti%mu ) / &
         gfct( (parti%nu+1.0)/parti%mu ) ) * &
         (gfct( (parti%nu+1.0)/parti%mu ) / gfct( (parti%nu+2.0)/parti%mu ) ) ** parti%b_geo
  END FUNCTION D_average_factor

  !*******************************************************************************
  ! Fundamental physical relations
  !*******************************************************************************

  ! Molecular diffusivity of water vapor
  ELEMENTAL FUNCTION diffusivity(T,p) result(D_v)
    REAL(wp), INTENT(IN) :: T,p 
    REAL(wp) :: D_v
    ! This is D_v = 8.7602e-5 * T_a**(1.81) / p_a
    D_v = 8.7602e-5 * exp(1.81*log(T)) / p
    RETURN
  END FUNCTION diffusivity

  !*******************************************************************************
  ! saturation pressure over ice and liquid water                                *
  !*******************************************************************************

  !  ELEMENTAL REAL(wp) FUNCTION e_es (ta)
  !    REAL(wp), INTENT(IN) :: ta
  !    e_es  = e_3 * EXP (A_e * (ta - T_3) / (ta - B_e))
  !  END FUNCTION e_es_old

  !  ELEMENTAL REAL(wp) FUNCTION e_ws (ta)
  !    REAL(wp), INTENT (IN) :: ta
  !    e_ws  = e_3 * EXP (A_w * (ta - T_3) / (ta - B_w))
  !  END FUNCTION e_ws_old

  ! FUNCTION e_ws_vec (ta,idim,jdim)
  !   INTEGER :: idim, jdim
  !   REAL(wp)               :: e_ws_vec(idim,jdim)
  !   REAL(wp), INTENT (IN)  :: ta(idim,jdim)
  !   e_ws_vec  = e_3 * EXP (A_w * (ta - T_3) / (ta - B_w))
  ! END FUNCTION e_ws_vec

  ! FUNCTION e_es_vec (ta,idim,jdim)
  !   INTEGER :: idim, jdim
  !   REAL(wp)               :: e_es_vec(idim,jdim)
  !   REAL(wp), INTENT (IN)  :: ta(idim,jdim)
  !   e_es_vec  = e_3 * EXP (A_e * (ta - T_3) / (ta - B_e))
  ! END FUNCTION e_es_vec

  SUBROUTINE autoconversionSB ()
    !*******************************************************************************
    ! Autoconversion of Seifert and Beheng (2001, Atmos. Res.)                     *
    !*******************************************************************************
    INTEGER         :: i,k
    INTEGER,  SAVE  :: firstcall
    REAL(wp), SAVE  :: k_au,k_sc
    REAL(wp)        :: q_c, q_r, n_c, x_c, nu, mu, tau, phi, x_s, au, sc

    REAL(wp), PARAMETER :: k_c  = 9.44e+9_wp   !..Long-Kernel
    REAL(wp), PARAMETER :: k_1  = 6.00e+2_wp   !..Parameter for Phi
    REAL(wp), PARAMETER :: k_2  = 0.68e+0_wp   !..Parameter fof Phi
    REAL(wp), PARAMETER :: eps  = 1.00e-25_wp

    IF (isdebug) THEN
      WRITE(txt,*) "autoconversionSB" ; CALL message(routine,TRIM(txt)) 
    END IF

    x_s = cloud%x_max     

    IF (firstcall.NE.1) THEN
       nu = cloud%nu
       mu = cloud%mu
       IF (mu == 1.0) THEN 
          !.. see SB2001
          k_au  = k_c / (20.0_wp*x_s) * (nu+2.0)*(nu+4.0)/(nu+1.0)**2
          k_sc  = k_c * (nu+2.0)/(nu+1.0)
       ELSE 
          !.. see Eq. (3.44) of Seifert (2002)
          k_au = k_c / (20.0_wp*x_s)                                       &
               & * ( 2.0_wp * gfct((nu+4.0)/mu)**1                          &
               &            * gfct((nu+2.0)/mu)**1 * gfct((nu+1.0)/mu)**2   &
               &   - 1.0_wp * gfct((nu+3.0)/mu)**2 * gfct((nu+1.0)/mu)**2 ) &
               &   / gfct((nu+2.0)/mu)**4
          k_sc = k_c * moment_gamma(cloud,2)
       ENDIF
    END IF

    DO k = kstart,kend
       DO i = istart,iend

          q_c = q_cloud(i,k) 
          IF (q_c > q_crit) THEN

            n_c = n_cloud(i,k)
            q_r = q_rain(i,k) 
            x_c = cloud%meanmass(q_c,n_c)

            au  = k_au * q_c**2 * x_c**2 * dt * rrho_c(i,k)
            tau = MIN(MAX(1.0-q_c/(q_c+q_r+eps),eps),0.9_wp)
            phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
            au  = au * (1.0 + phi/(1.0 - tau)**2)

            au  = MAX(MIN(q_c,au),0.0_wp)

            sc  = k_sc * q_c**2 * dt * rrho_c(i,k)

            n_rain(i,k)  = n_rain(i,k)  + au / x_s
            q_rain(i,k)  = q_rain(i,k)  + au
            n_cloud(i,k) = n_cloud(i,k) - MIN(n_c,sc) 
            q_cloud(i,k) = q_cloud(i,k) - au             
            
          ENDIF
       END DO
    END DO

  END SUBROUTINE autoconversionSB

  SUBROUTINE accretionSB ()
    !*******************************************************************************
    ! Accretion of Seifert and Beheng (2001, Atmos. Res.)                          *
    !*******************************************************************************
    INTEGER     :: i, k
    REAL(wp)    :: ac
    REAL(wp)    :: q_c, q_r, tau, phi, n_c, x_c

    REAL(wp), PARAMETER :: k_r = 5.78_wp       ! kernel
    REAL(wp), PARAMETER :: k_1 = 5.00e-04_wp   ! Phi function
    REAL(wp), PARAMETER :: eps = 1.00e-25_wp

    IF (isdebug) THEN
      WRITE(txt,*) "accretionSB" ; CALL message(routine,TRIM(txt)) 
    END IF

    DO k = kstart,kend
       DO i = istart,iend

          n_c = n_cloud(i,k)       
          q_c = q_cloud(i,k) 
          q_r = q_rain(i,k)  
          IF (q_c > 0.0_wp.AND.q_r > 0.0_wp) THEN
             
             !..accretion rate of SB2001
             tau = MIN(MAX(1.0-q_c/(q_c+q_r+eps),eps),1.0_wp)
             phi = (tau/(tau+k_1))**4
             ac  = k_r *  q_c * q_r * phi * dt

             ac = MIN(q_c,ac)

             x_c = cloud%meanmass(q_c,n_c)

             q_rain(i,k)  = q_rain(i,k)  + ac
             q_cloud(i,k) = q_cloud(i,k) - ac
             n_cloud(i,k) = n_cloud(i,k) - MIN(n_c,ac/x_c)
          ENDIF
       END DO
    END DO

  END SUBROUTINE accretionSB

  SUBROUTINE rain_selfcollectionSB ()
    !*******************************************************************************
    ! Selfcollection of Seifert and Beheng (2001, Atmos. Res.)                     *
    !*******************************************************************************
    INTEGER     :: i, k
    REAL(wp)    :: sc, br
    REAL(wp)    :: q_r, n_r, x_r, d_r

    !..Parameters based on Seifert (2008, JAS)
    REAL(wp), PARAMETER :: D_br = 1.10e-3_wp
    REAL(wp), PARAMETER :: k_rr = 4.33e+0_wp
    REAL(wp), PARAMETER :: k_br = 1.00e+3_wp

    IF (isdebug) THEN
      WRITE(txt,*) "rain_selfcollectionSB" ; CALL message(routine,TRIM(txt)) 
    END IF

    DO k = kstart,kend
       DO i = istart,iend
          
          n_r = n_rain(i,k)
          q_r = q_rain(i,k)          
          IF (q_r > 0.0_wp) THEN
            x_r = rain%meanmass(q_r,n_r)
            D_r = rain%diameter(x_r)

            !..Selfcollection as in SB2001
            sc = k_rr *  n_r * q_r * rrho_04(i,k) * dt
            
            !..Breakup as in Seifert (2008, JAS), Eq. (A13)
            br = 0.0_wp
            IF (D_r.GT.0.30e-3_wp) THEN
               br = (k_br * (D_r - D_br) + 1.0)  * sc                   
            ENDIF
           
            n_rain(i,k) = n_r - MIN(n_r,sc-br)
         ENDIF
      END DO
   END DO

  END SUBROUTINE rain_selfcollectionSB

  SUBROUTINE autoconversionKB ()

    INTEGER  :: i,k
    REAL(wp) :: q_c, x_c, nu_c, n_c, k_a, x_s, au

    nu_c = 9.59
    x_s  = cloud%x_max                  
    k_a  = 6.0d+25 * nu_c**(-1.7)

    !..Parameterization of Beheng (1994)
    DO k = kstart,kend
       DO i = istart,iend
          
          q_c = q_cloud(i,k)
          n_c = n_cloud(i,k)
          IF (q_c > q_crit) THEN
             
             x_c = cloud%meanmass(q_c,n_c)
             
             !..Berechnung der Autokonversionsrate nach Beheng (1994)
             au = k_a * (x_c*1e3)**(3.3) * (q_c*1e-3)**(1.4) * dt * 1e3
             au = MIN(q_c,au)
             
             n_rain(i,k)  = n_rain(i,k)  + au / x_s
             q_rain(i,k)  = q_rain(i,k)  + au
             n_cloud(i,k) = n_cloud(i,k) - au / x_s * 2.0
             q_cloud(i,k) = q_cloud(i,k) - au             
             
          ENDIF
       END DO
    END DO

  END SUBROUTINE autoconversionKB

  SUBROUTINE accretionKB ()

    !..Parameter of Beheng (1994)
    REAL(wp), PARAMETER :: k_r = 6.00d+00

    INTEGER      :: i,k
    REAL(wp)     :: ac
    REAL(wp)     :: q_c, q_r

    DO k = kstart,kend
       DO i = istart,iend

          q_c = q_cloud(i,k)
          q_r = q_rain(i,k) 
          
          IF (q_c > q_crit .and. q_r > q_crit) THEN

             ac = k_r *  q_c * q_r * dt
             ac = MIN(q_c,ac)

             q_rain(i,k)  = q_rain(i,k)  + ac
             q_cloud(i,k) = q_cloud(i,k) - ac             

          ENDIF
       END DO
    END DO

  END SUBROUTINE accretionKB

  SUBROUTINE autoconversionKK ()

    INTEGER             :: i,k
    REAL(wp)            :: q_c, x_s, au, n_c
    REAL(wp), PARAMETER :: k_a  = 1350.0_wp

    x_s  = cloud%x_max   

    !..Parameterization of Khairoutdinov and Kogan (2000), MWR 128, 229-243
    DO k = kstart,kend
       DO i = istart,iend

          q_c = q_cloud(i,k) 
          IF (q_c > q_crit) THEN

             n_c = n_cloud(i,k) * 1e6  ! in 1/cm3

             !..autoconversion rate of KK2000
             au = k_a * q_c**(2.47) * n_c**(-1.79) * dt
             au = MIN(q_c,au)

             n_rain(i,k)  = n_rain(i,k)  + au / x_s
             q_rain(i,k)  = q_rain(i,k)  + au
             n_cloud(i,k) = n_cloud(i,k) - au / x_s * 2.0
             q_cloud(i,k) = q_cloud(i,k) - au             

          ENDIF
       END DO
    END DO

  END SUBROUTINE autoconversionKK

  SUBROUTINE accretionKK ()

    INTEGER             :: i,k
    REAL(wp)            :: ac, q_c, q_r
    REAL(wp), PARAMETER :: k_a = 6.70d+01

    !..Parameterization of Khairoutdinov and Kogan (2000), MWR 128, 229-243
    DO k = kstart,kend
       DO i = istart,iend

          q_c = q_cloud(i,k)
          q_r = q_rain(i,k) 

          IF (q_c > q_crit .AND. q_r > q_crit) THEN

             ac = k_a *  (q_c * q_r)**1.15 * dt * 1e3
             ac = MIN(q_c,ac)
             
             q_rain(i,k)  = q_rain(i,k)  + ac
             q_cloud(i,k) = q_cloud(i,k) - ac             
          ENDIF
       END DO
    END DO

  END SUBROUTINE accretionKK

  SUBROUTINE rain_evaporation ()
    !*******************************************************************************
    ! Evaporation of rain based on Seifert (2008, J. Atmos. Sci.)                  *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a,p_a,e_sw,s_sw,g_d,eva_q,eva_n
    REAL(wp)            :: q_r,n_r,x_r,e_d,f_v
    REAL(wp)            :: mue,d_m,gamma_eva,lam,d_vtp,gfak
    REAL(wp)            :: aa,bb,cc,mm
    REAL(wp), SAVE      :: a_q,b_q 

    REAL(wp), PARAMETER :: D_br = 1.1e-3_wp

    aa = rain_coeffs%alfa
    bb = rain_coeffs%beta
    cc = rain_coeffs%gama

    IF (firstcall.NE.1) THEN
      a_q = vent_coeff_a(rain,1)
      b_q = vent_coeff_b(rain,1)
      firstcall = 1
      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "rain_evaporation:" ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     nu    = ",rain%nu ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     mu    = ",rain%mu ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     a_geo = ",rain%a_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     b_geo = ",rain%b_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     a_q   = ",a_q ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     b_q   = ",b_q ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     cmu0  = ",rain_coeffs%cmu0 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     cmu1  = ",rain_coeffs%cmu1 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     cmu2  = ",rain_coeffs%cmu2 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     cmu3  = ",rain_coeffs%cmu3 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     cmu4  = ",rain_coeffs%cmu4 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     g_fak = ",rain_gfak ; CALL message(routine,TRIM(txt))
      END IF
    ELSEIF (isdebug) THEN
      WRITE(txt,*) "rain_evaporation " ; CALL message(routine,TRIM(txt)) 
    END IF

    DO k = kstart,kend
       DO i = istart,iend
          
          q_r  = q_rain(i,k)
          n_r  = n_rain(i,k)
          p_a  = p_p(i,k)
          T_a  = T_p(i,k)

          e_d  = qv(i,k) * R_d * T_a
          e_sw = e_ws(T_a)
          s_sw = e_d / e_sw - 1.0 

          IF (s_sw < 0.0 .AND. q_r > 0.0_wp .AND. q_cloud(i,k) < q_crit) THEN

             D_vtp = diffusivity(T_a,p_a)

             ! note that 2*pi is correct, because c_r = 1/2 is assumed
             g_d = 2.0*pi / ( L_wd**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_sw) )

             x_r = rain%meanmass(q_r,n_r)
             D_m = rain%diameter(x_r)

             ! Eq. (20) of Seifert (2008)
             IF (D_m.LE.rain_coeffs%cmu3) THEN 
                mue = rain_coeffs%cmu0*TANH((4.*rain_coeffs%cmu2*(D_m-rain_coeffs%cmu3))**rain_coeffs%cmu5) &
                     & + rain_coeffs%cmu4
             ELSE
                mue = rain_coeffs%cmu1*TANH((1.*rain_coeffs%cmu2*(D_m-rain_coeffs%cmu3))**rain_coeffs%cmu5) &
                     & + rain_coeffs%cmu4
             ENDIF

             ! Eq. (A8)
             lam = exp(1./3.*log(pi6*rho_w*(mue+3.0)*(mue+2.0)*(mue+1.0)/x_r))

              ! chebyshev approximation of Gamma(mue+5/2)/Gamma(mue+2)
             gfak =  0.1357940435E+01 &
                  &  + mue * ( +0.3033273220E+00  &
                  &  + mue * ( -0.1299313363E-01  &
                  &  + mue * ( +0.4002257774E-03  &
                  &  - mue * 0.4856703981E-05 ) ) )

             mm = mue+5.0/2.0
             
             ! Eq. (A7) rewritten with (A5) and (A9)
             f_v  = rain%a_ven + rain%b_ven * N_sc**n_f * gfak                     &
                  &            * sqrt(aa/nu_l*rrho_04(i,k) / lam)                  &
                  &    * (1.0 - 1./2.  * (bb/aa)**1 * exp(mm*log(lam/(1.*cc+lam))) &
                  &           - 1./8.  * (bb/aa)**2 * exp(mm*log(lam/(2.*cc+lam))) &
                  &           - 1./16. * (bb/aa)**3 * exp(mm*log(lam/(3.*cc+lam))) &
                  &           - 5./127.* (bb/aa)**4 * exp(mm*log(lam/(4.*cc+lam))) )
         
             IF (rain_gfak.GT.0) THEN
                ! Eq. (23)
                gamma_eva = rain_gfak * (D_br/D_m) * exp(-0.2*mue)
             ELSE
                gamma_eva = 1.0
             END IF

             ! Eq (A5) with (A9) and Gamma(mue+2) from (A7) 
             eva_q = g_d * n_r * (mue+1.0) / lam * f_v * s_sw * dt
             eva_n = gamma_eva * eva_q / x_r
             
             eva_q = MAX(-eva_q,0.0_wp) 
             eva_n = MAX(-eva_n,0.0_wp) 
             
             eva_q = MIN(eva_q,q_r) 
             eva_n = MIN(eva_n,n_r) 

             qv(i,k)     = qv(i,k)     + eva_q
             q_rain(i,k) = q_rain(i,k) - eva_q
             n_rain(i,k) = n_rain(i,k) - eva_n
          END IF
       END DO
    END DO

  END SUBROUTINE rain_evaporation

  SUBROUTINE graupel_evaporation ()
    !*******************************************************************************
    ! Evaporation from melting graupel, see SB2006                                 *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a,e_sw,s_sw,g_d,eva
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g,f_v,e_d
    REAL(wp), SAVE      :: c_g       
    REAL(wp), SAVE      :: a_f,b_f   

    IF (firstcall.NE.1) THEN
      a_f = vent_coeff_a(graupel,1)
      b_f = vent_coeff_b(graupel,1) * N_sc**n_f / sqrt(nu_l)
      c_g = 1.0 / graupel%cap
      firstcall = 1
      IF (isdebug) THEN
        WRITE(txt,'(A)') "graupel_evaporation:"  ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     a_f = ",a_f ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     b_f = ",b_f ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     c_g = ",c_g ; CALL message(routine,TRIM(txt))
      END IF
    ELSEIF (isdebug) THEN
      WRITE(txt,*) "graupel_evaporation " ; CALL message(routine,TRIM(txt)) 
    END IF

    DO k = kstart,kend
       DO i = istart,iend

          q_g = q_graupel(i,k)
          T_a = T_p(i,k)

          IF (q_g > 0.0_wp .AND. T_a > T_3) THEN

            e_d  = qv(i,k) * R_d * T_a 
            e_sw = e_ws(T_a)
            s_sw = e_d / e_sw - 1.0       

            !..Eq. (37) of SB2006, note that 4*pi is correct because c_g is used below
            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_g = n_graupel(i,k)          
            x_g = graupel%meanmass(q_g,n_g)
            d_g = graupel%diameter(x_g)
            v_g = graupel%velocity(x_g)

            !..note that a_f includes more than just ventilation, do never ever set f_v=1
            f_v = a_f + b_f * sqrt(v_g*d_g)   

            eva = g_d * n_g * c_g * d_g * f_v * s_sw * dt

            eva = MAX(-eva,0.0_wp) 
            eva = MIN(eva,q_g) 

            qv(i,k)        = qv(i,k)        + eva
            q_graupel(i,k) = q_graupel(i,k) - eva

          END IF
       END DO
    END DO

  END SUBROUTINE graupel_evaporation

  SUBROUTINE hail_evaporation ()
    !*******************************************************************************
    ! Evaporation of melting graupel, see SB2006                                   *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a,e_sw,s_sw,g_d,eva
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h,f_v,e_d
    REAL(wp), SAVE      :: c_h        
    REAL(wp), SAVE      :: a_f,b_f    

    IF (firstcall.NE.1) THEN
      a_f = vent_coeff_a(hail,1)
      b_f = vent_coeff_b(hail,1) * N_sc**n_f / sqrt(nu_l)
      c_h = 1.0 / hail%cap
      firstcall = 1
      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "hail_evaporation:" ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     a_f = ",a_f   ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     b_f = ",b_f   ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     c_h = ",c_h   ; CALL message(routine,TRIM(txt))
      END IF
    ELSEIF (isdebug) THEN
      WRITE(txt,*) "hail_evaporation " ; CALL message(routine,TRIM(txt)) 
    END IF

    DO k = kstart,kend
       DO i = istart,iend

          q_h = q_hail(i,k)
          T_a = T_p(i,k)

          IF (q_h > 0.0_wp .AND. T_a > T_3) THEN

            e_d  = qv(i,k) * R_d * T_a 
            e_sw = e_ws(T_a)
            s_sw = e_d / e_sw - 1.0   

            !..Eq. (37) of SB2006, note that 4*pi is correct because c_h is used below
            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_h = n_hail(i,k)            
            x_h = hail%meanmass(q_h,n_h)
            d_h = hail%diameter(x_h)
            v_h = hail%velocity(x_h)

            f_v  = a_f + b_f * sqrt(v_h*D_h)

            eva = g_d * n_h * c_h * d_h * f_v * s_sw * dt

            eva = MAX(-eva,0.0_wp) 
            eva = MIN(eva,q_h) 

            qv(i,k)     = qv(i,k)     + eva
            q_hail(i,k) = q_hail(i,k) - eva

          END IF
       END DO
    END DO
  END SUBROUTINE hail_evaporation

  SUBROUTINE snow_evaporation ()
    !*******************************************************************************
    ! Evaporation of melting snow, see SB2006                                      *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a,e_sw,s_sw,g_d,eva
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s,f_v,e_d
    REAL(wp), SAVE      :: c_s       
    REAL(wp), SAVE      :: a_f,b_f   

    IF (firstcall.NE.1) THEN
      a_f = vent_coeff_a(snow,1)
      b_f = vent_coeff_b(snow,1) * N_sc**n_f / sqrt(nu_l)
      c_s = 1.0/snow%cap
      firstcall = 1
      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "snow_evaporation:" ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     a_f = ",a_f ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     b_f = ",b_f ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     c_s = ",c_s ; CALL message(routine,TRIM(txt))
      END IF
    ELSEIF (isdebug) THEN
      WRITE(txt,*) "snow_evaporation " ; CALL message(routine,TRIM(txt)) 
    END IF

    DO k = kstart,kend
       DO i = istart,iend

          q_s = q_snow(i,k)
          T_a = T_p(i,k) 

          IF (q_s > 0.0_wp .AND. T_a > T_3) THEN

            e_d  = qv(i,k) * R_d * T_a 
            e_sw = e_ws(T_a)
            s_sw = e_d / e_sw - 1.0  

            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_s = n_snow(i,k)   

            x_s = snow%meanmass(q_s,n_s)
            D_s = snow%diameter(x_s)
            v_s = snow%velocity(x_s) * rrho_04(i,k)

            f_v  = a_f + b_f * sqrt(v_s*D_s)

            eva = g_d * n_s * c_s * d_s * f_v * s_sw * dt

            eva = MAX(-eva,0.0_wp) 
            eva = MIN(eva,q_s) 

            qv(i,k)     = qv(i,k)     + eva
            q_snow(i,k) = q_snow(i,k) - eva

          END IF
       END DO
    END DO
  END SUBROUTINE snow_evaporation

  SUBROUTINE cloud_freeze ()
    !*******************************************************************************
    ! This is only the homogeneous freezing of liquid water droplets.              *
    ! Immersion freezing and homogeneous freezing of liquid aerosols are           *
    ! treated in the subroutine ice_nucleation_homhet()                            *
    !*******************************************************************************
    INTEGER             :: i, k
    REAL(wp)            :: fr_q,fr_n,T_a,q_c,x_c,n_c,j_hom,T_c
    REAL(wp), SAVE      :: coeff_z
    INTEGER,  SAVE      :: firstcall

    IF (firstcall.NE.1) THEN
       firstcall = 1
       coeff_z   = moment_gamma(cloud,2)  
    END IF

    DO k = kstart,kend
       DO i = istart,iend

          T_a = T_p(i,k)
          IF (T_a < T_3) THEN

             T_c = T_a - T_3
             q_c = q_cloud(i,k)
             n_c = n_cloud(i,k)
             IF (q_c > 0.0_wp) THEN
                IF (T_c < -50.0_wp) THEN            
                   fr_q = q_c             !..instantaneous freezing
                   fr_n = n_c             !..below -50 C
                ELSE          
                   x_c = cloud%meanmass(q_c,n_c)                            
                
                   !..Hom. freezing based on Jeffrey und Austin (1997), see also Cotton und Field (2001)
                   IF (T_c > -30.0_wp) THEN            
                      j_hom = 1.0e6_wp/rho_w &
                           &  * EXP(-7.63-2.996*(T_c+30.0))           !..J in 1/(m3 s) 
                   ELSE
                      j_hom = 1.0e6_wp/rho_w &
                           &  * EXP(-243.4-14.75*T_c-0.307*T_c**2-0.00287*T_c**3-0.0000102*T_c**4)
                   ENDIF

                   fr_n  = j_hom * q_c *  dt
                   fr_q  = j_hom * q_c * x_c * dt * coeff_z
                   fr_q  = MIN(fr_q,q_c)
                   fr_n  = MIN(fr_n,n_c)
                END IF

                q_cloud(i,k) = q_cloud(i,k) - fr_q
                n_cloud(i,k) = n_cloud(i,k) - fr_n

                fr_n  = MAX(fr_n,fr_q/cloud%x_max)
              
                !..special treatment for constant drop number
                IF (nuc_c_typ .EQ. 0) THEN
                   ! ... force upper bound in cloud_freeze'
                   fr_n = MAX(MIN(fr_n,qnc_const-n_ice(i,k)),0.0_wp)
                ENDIF
              
                q_ice(i,k)   = q_ice(i,k) + fr_q
                n_ice(i,k)   = n_ice(i,k) + fr_n
             ENDIF
          END IF
       END DO
    END DO
  END SUBROUTINE cloud_freeze

  SUBROUTINE ice_nucleation_homhet()
  !*******************************************************************************
  !                                                                              *
  ! Homogeneous and heterogeneous ice nucleation                                 *
  !                                                                              *
  ! Nucleation scheme is based on the papers:                                    *
  !                                                                              *
  ! "A parametrization of cirrus cloud formation: Homogenous                     *
  ! freezing of supercooled aerosols" by B. Kaercher and                         *
  ! U. Lohmann 2002 (KL02 hereafter)                                             *
  !                                                                              *
  ! "Physically based parameterization of cirrus cloud formation                 *
  ! for use in global atmospheric models" by B. Kaercher, J. Hendricks           *
  ! and U. Lohmann 2006 (KHL06 hereafter)                                        *
  !                                                                              *
  ! and Phillips et al. (2008) with extensions                                   *
  !                                                                              *
  ! implementation by Carmen Koehler and AS                                      *
  !*******************************************************************************
  INTEGER             :: i,k,nuc_typ
  REAL(wp)            :: nuc_n, nuc_q
  REAL(wp)            :: T_a,p_a,ssi
  REAL(wp)            :: q_i,n_i,x_i,r_i
  REAL(wp)            :: ndiag

  ! switch for version of Phillips et al. scheme 
  ! (but make sure you have the correct INCLUDE file)
  INTEGER             :: iphillips = 2010

  ! some more constants needed for homogeneous nucleation scheme
  REAL(wp), PARAMETER ::            &
    r_0     = 0.25e-6_wp          , &    ! aerosol particle radius prior to freezing
    alpha_d = 0.5_wp              , &    ! deposition coefficient (KL02; Spichtinger & Gierens 2009)
    M_w     = 18.01528e-3_wp      , &    ! molecular mass of water [kg/mol]
    M_a     = 28.96e-3_wp         , &    ! molecular mass of air [kg/mol]
    ma_w    = M_w / N_avo         , &    ! mass of water molecule [kg]
    svol    = ma_w / rho_ice             ! specific volume of a water molecule in ice
  
  REAL(wp)  :: e_si
  REAL(wp)  :: ni_hom,ri_hom,mi_hom
  REAL(wp)  :: v_th,n_sat,flux,phi,cool,tau,delta,w_pre,scr
  REAL(wp)  :: ctau, tau_g,acoeff(3),bcoeff(2), ri_dot
  REAL(wp)  :: kappa,sqrtkap,ren,R_imfc,R_im,R_ik,ri_0

  ! variables for interpolation in look-up table (real is good enough here)
  REAL      :: xt,xs 
  INTEGER   :: ss,tt

  LOGICAL   :: use_homnuc = .true.

  REAL(wp), DIMENSION(3) :: infrac

  nuc_typ = nuc_i_typ
    
  SELECT CASE (nuc_typ)
  CASE(0)
    ! Heterogeneous nucleation ONLY 
    use_homnuc = .FALSE.
  CASE(1:9)
    ! Homog. and het. nucleation"
    use_homnuc = .TRUE.
  END SELECT

  IF (isdebug) THEN
    IF (.NOT.use_homnuc) THEN
       WRITE(txt,*) "ice_nucleation_homhet: Heterogeneous nucleation only"
    ELSE 
       WRITE(txt,*) "ice_nucleation_homhet: Homogeneous and heterogeneous nucleation"
    END IF
    CALL message(routine,TRIM(txt))
  END IF

  ! Heterogeneous nucleation using Phillips et al. scheme

  ! possible pre-defined choices
  IF (iphillips == 2010) THEN
     IF (nuc_typ.EQ.4) THEN  ! with no organics and rather high soot, coming close to Meyers formula at -20 C
        na_dust  = 160.e4_wp    ! initial number density of dust [1/m3]
        na_soot  =  30.e6_wp    ! initial number density of soot [1/m3]
        na_orga  =   0.e0_wp    ! initial number density of organics [1/m3]
     END IF
     IF (nuc_typ.EQ.5) THEN     ! with some organics and rather high soot, 
        na_dust  = 160.e4_wp    !          coming close to Meyers formula at -20 C
        na_soot  =  25.e6_wp 
        na_orga  =  30.e6_wp 
     END IF
     IF (nuc_typ.EQ.5) THEN     ! no organics, no soot, coming close to DeMott et al. 2010 at -20 C
        na_dust  =  70.e4_wp    ! i.e. roughly one order in magnitude lower than Meyers
        na_soot  =   0.e6_wp 
        na_orga  =   0.e6_wp 
     END IF
  END IF
      
  DO k = kstart,kend
     DO i = istart,iend
        p_a  = p_p(i,k)
        T_a  = T_p(i,k)
        e_si = e_es(T_a)
        ssi  = qv(i,k) * R_d * T_a / e_si
        
        IF (T_p(i,k) < T_nuc .AND. ssi > 0.0_wp  &
             & .AND. ( n_ice(i,k)+n_snow(i,k) < ni_het_max ) ) THEN

           IF (q_cloud(i,k) > 0.0_wp) THEN

              ! immersion freezing at water saturation
              xt = (274.- real(T_p(i,k)))  / ttstep
              xt = MIN(xt,real(ttmax-1))      
              tt = INT(xt)
              infrac(1) = (tt+1-xt) * afrac_dust(tt,99) + (xt-tt) * afrac_dust(tt+1,99) 
              infrac(2) = (tt+1-xt) * afrac_soot(tt,99) + (xt-tt) * afrac_soot(tt+1,99) 
              infrac(3) = (tt+1-xt) * afrac_orga(tt,99) + (xt-tt) * afrac_orga(tt+1,99)                
           ELSE
              ! deposition nucleation below water saturation

              ! calculate indices used for 2D look-up tables
              xt = (274.- real(T_p(i,k)))  / ttstep
              xs = 100. * real(ssi) / ssstep    
              xt = MIN(xt,real(ttmax-1))
              xs = MIN(xs,real(ssmax-1))          
              tt = INT(xt)
              ss = INT(xs)
                
              ! bi-linear interpolation in look-up tables
              infrac(1) = (tt+1-xt)*(ss+1-xs) * afrac_dust(tt,ss  ) + (xt-tt)*(ss+1-xs) * afrac_dust(tt+1,ss  ) &
                   &    + (tt+1-xt)*(xs-ss)   * afrac_dust(tt,ss+1) + (xt-tt)*(xs-ss)   * afrac_dust(tt+1,ss+1)
              infrac(2) = (tt+1-xt)*(ss+1-xs) * afrac_soot(tt,ss  ) + (xt-tt)*(ss+1-xs) * afrac_soot(tt+1,ss  ) &
                   &    + (tt+1-xt)*(xs-ss)   * afrac_soot(tt,ss+1) + (xt-tt)*(xs-ss)   * afrac_soot(tt+1,ss+1)
              infrac(3) = (tt+1-xt)*(ss+1-xs) * afrac_orga(tt,ss  ) + (xt-tt)*(ss+1-xs) * afrac_orga(tt+1,ss  ) &
                   &    + (tt+1-xt)*(xs-ss)   * afrac_orga(tt,ss+1) + (xt-tt)*(xs-ss)   * afrac_orga(tt+1,ss+1)
           ENDIF

           ndiag = na_dust * infrac(1) + na_soot * infrac(2) + na_orga * infrac(3)
           ndiag = MIN(ndiag,ni_het_max)
           
           nuc_n = MAX(ndiag - (n_ice(i,k) + n_snow(i,k)),0.0_wp)
           nuc_q = MIN(nuc_n * ice%x_min, qv(i,k))
           
           n_ice(i,k) = n_ice(i,k) + nuc_n
           q_ice(i,k) = q_ice(i,k) + nuc_q
           qv(i,k)    = qv(i,k)    - nuc_q

        ENDIF

     END DO
  END DO

  ! Homogeneous nucleation using KHL06 approach
  IF (use_homnuc) THEN
     DO k = kstart,kend
        DO i = istart,iend
          p_a  = p_p(i,k)
          T_a  = T_p(i,k)
          e_si = e_es(T_a)
          ssi  = qv(i,k) * R_d * T_a / e_si

          ! critical supersaturation for homogeneous nucleation
          scr  = 2.349 - T_a / 259.00
          
          IF (ssi > scr .AND. n_ice(i,k) < ni_hom_max ) THEN

            n_i = n_ice(i,k)
            q_i = q_ice(i,k)
            x_i = ice%meanmass(q_i,n_i) 
            r_i = (x_i/(4./3.*pi*rho_ice))**(1./3.)

            v_th  = SQRT( 8.*k_b*T_a/(pi*ma_w) ) 
            flux  = alpha_d * v_th/4.
            n_sat = e_si / (k_b*T_a)  

            ! coeffs of supersaturation equation
            acoeff(1) = (L_ed * grav) / (cp * R_d * T_a**2) - grav/(R_l * T_a)
            acoeff(2) = 1.0/n_sat
            acoeff(3) = (L_ed**2 * M_w * ma_w)/(cp * p_a * T_a * M_a) 
          
            ! coeffs of depositional growth equation
            bcoeff(1) = flux * svol * n_sat * (ssi - 1.)
            bcoeff(2) = flux / diffusivity(T_a,p_a)

            ! pre-existing ice crystals included as reduced updraft speed
            ri_dot = bcoeff(1) / (1. + bcoeff(2) * r_i)
            R_ik   = (4 * pi) / svol * n_i * r_i**2 * ri_dot
            w_pre  = (acoeff(2) + acoeff(3) * ssi)/(acoeff(1) * ssi) * R_ik  ! KHL06 Eq. 19
            w_pre  = MAX(w_pre,0.0_wp)

            IF (w_p(i,k) > w_pre) THEN   ! homogenous nucleation event

              ! timescales of freezing event (see KL02, RM05, KHL06)
              cool    = grav / cp * w_p(i,k)
              ctau    = T_a * ( 0.004*T_a - 2. ) + 304.4         
              tau     = 1.0 / (ctau * cool)                       ! freezing timescale, eq. (5)
              delta   = (bcoeff(2) * r_0)                         ! dimless aerosol radius, eq.(4)  
              tau_g   = (bcoeff(1) / r_0) / (1 + delta)           ! timescale for initial growth, eq.(4)
              phi     = acoeff(1)*ssi / ( acoeff(2) + acoeff(3)*ssi) * (w_p(i,k) - w_pre) 
     
              ! monodisperse approximation following KHL06
              kappa   = 2. * bcoeff(1) * bcoeff(2) * tau / (1.+ delta)**2  ! kappa, Eq. 8 KHL06
              sqrtkap = SQRT(kappa)                                        ! root of kappa
              ren     = 3. * sqrtkap / ( 2. + SQRT(1.+9.*kappa/pi) )       ! analy. approx. of erfc by RM05
              R_imfc  = 4. * pi * bcoeff(1)/bcoeff(2)**2 / svol  
              R_im    = R_imfc / (1.+ delta) * ( delta**2 - 1. &
                & + (1.+0.5*kappa*(1.+ delta)**2) * ren/sqrtkap)           ! RIM Eq. 6 KHL06
              
              ! number concentration and radius of ice particles
              ni_hom  = phi / R_im                                         ! ni Eq.9 KHL06
              ri_0    = 1. + 0.5 * sqrtkap * ren                           ! for Eq. 3 KHL06
              ri_hom  = (ri_0 * (1. + delta) - 1. ) / bcoeff(2)            ! Eq. 3 KHL06 * REN = Eq.23 KHL06
              mi_hom  = (4./3. * pi * rho_ice) * ni_hom * ri_hom**3
              mi_hom  = MAX(mi_hom,ice%x_min)

              nuc_n = MAX(MIN(ni_hom, ni_hom_max), 0.0_wp)
              nuc_q = MIN(nuc_n * mi_hom, qv(i,k))
              
              n_ice(i,k) = n_ice(i,k) + nuc_n
              q_ice(i,k) = q_ice(i,k) + nuc_q
              qv(i,k)  = qv(i,k)  - nuc_q
              
            END IF
          END IF
       ENDDO
    ENDDO
  END IF

  END SUBROUTINE ice_nucleation_homhet

  SUBROUTINE vapor_dep_relaxation(nsize, nlev, dt_local, dep_rate_ice, dep_rate_snow)
    !*******************************************************************************
    ! Deposition and sublimation                                                   *
    !*******************************************************************************
    INTEGER,  INTENT(IN) :: nsize, nlev
    REAL(wp), INTENT(IN) :: dt_local
    REAL(wp), INTENT(INOUT), DIMENSION(nsize,nlev) :: dep_rate_ice, dep_rate_snow

    REAL(wp), DIMENSION(nsize,nlev) :: &
         & s_si,g_i,dep_ice,dep_snow,dep_graupel,dep_hail

    INTEGER             :: i,k
    REAL(wp)            :: D_vtp
    REAL(wp)            :: zdt,qvsidiff,Xi_i,Xfac
    REAL(wp)            :: tau_i_i,tau_s_i,tau_g_i,tau_h_i
    REAL(wp), PARAMETER :: eps  = 1.d-20
    REAL(wp)            :: T_a         
    REAL(wp)            :: e_si            !..saturation water pressure over ice
    REAL(wp)            :: e_sw            !..saturation water pressure over liquid
    REAL(wp)            :: e_d,p_a,dep_sum

    IF (isdebug) THEN
       WRITE(txt,*) "vapor_deposition_growth " ; CALL message(routine,TRIM(txt))
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend
          p_a  = p_p(i,k)
          T_a  = T_p(i,k)
          e_d  = qv(i,k) * R_d * T_a
          e_si = e_es(T_a)
          e_sw = e_ws(T_a)
          s_si(i,k) = e_d / e_si - 1.0    !..supersaturation over ice
          D_vtp = diffusivity(T_a,p_a)    !  D_v = 8.7602e-5 * T_a**(1.81) / p_a
          IF (T_a < T_3) THEN
             g_i(i,k) = 4.0*pi / ( L_ed**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_si) )
          ELSE
             g_i(i,k)  = 0.0
             s_si(i,k) = 0.0
          ENDIF
       ENDDO
    ENDDO

    DO k = kstart,kend
       dep_ice(:,k)     = 0.0
       dep_snow(:,k)    = 0.0
       dep_graupel(:,k) = 0.0
       dep_hail(:,k)    = 0.0
    END DO

    CALL vapor_deposition_ice()
    CALL vapor_deposition_snow()
    CALL vapor_deposition_graupel()
    IF (ice_typ > 1) CALL vapor_deposition_hail()

    zdt = 1.0/dt_local

    DO k = kstart,kend
       DO i = istart,iend

          T_a  = T_p(i,k)

          ! Deposition only below T_3, evaporation of melting particles at warmer T is treated elsewhere
          IF (T_a < T_3) THEN

             ! Depositional growth with relaxation time-scale approach based on:
             ! "A New Double-Moment Microphysics Parameterization for Application in Cloud and
             ! Climate Models. Part 1: Description" by H. Morrison, J.A.Curry, V.I. Khvorostyanov
             
             qvsidiff  = qv(i,k) - e_es(T_a)/(R_d*T_a)

             if (abs(qvsidiff).lt.eps) then                
                dep_ice(i,k)     = 0.0
                dep_snow(i,k)    = 0.0
                dep_graupel(i,k) = 0.0               
                dep_hail(i,k)    = 0.0
                dep_sum          = 0.0
             else
             
                ! deposition rates are already multiplied with dt_local, therefore divide them here
                tau_i_i  = zdt/qvsidiff*dep_ice(i,k)
                tau_s_i  = zdt/qvsidiff*dep_snow(i,k)
                tau_g_i  = zdt/qvsidiff*dep_graupel(i,k)
                tau_h_i  = zdt/qvsidiff*dep_hail(i,k)
                
                Xi_i = ( tau_i_i + tau_s_i + tau_g_i + tau_h_i ) 
                
                if (Xi_i.lt.eps) then
                   Xfac = 0.0_wp
                else
                   Xfac =  qvsidiff / Xi_i * (1.0 - EXP(- dt_local*Xi_i))
                end if
                
                dep_ice(i,k)     = Xfac * tau_i_i
                dep_snow(i,k)    = Xfac * tau_s_i
                dep_graupel(i,k) = Xfac * tau_g_i
                dep_hail(i,k)    = Xfac * tau_h_i
             
                IF (qvsidiff < 0.0) THEN
                   dep_ice(i,k)     = MAX(dep_ice(i,k),    -q_ice(i,k))
                   dep_snow(i,k)    = MAX(dep_snow(i,k),   -q_snow(i,k))
                   dep_graupel(i,k) = MAX(dep_graupel(i,k),-q_graupel(i,k)) 
                   dep_hail(i,k)    = MAX(dep_hail(i,k),   -q_hail(i,k))  
                END IF

                dep_sum = dep_ice(i,k) + dep_graupel(i,k) + dep_snow(i,k) + dep_hail(i,k)
             END IF

             q_graupel(i,k) = q_graupel(i,k) + dep_graupel(i,k)
             
             q_ice(i,k)     = q_ice(i,k)     + dep_ice(i,k)
             q_snow(i,k)    = q_snow(i,k)    + dep_snow(i,k)
             IF (ice_typ > 1) THEN
                q_hail(i,k)    = q_hail(i,k)    + dep_hail(i,k)
             END IF

             qv(i,k) = qv(i,k) - dep_sum
             
             dep_rate_ice(i,k)  = dep_rate_ice(i,k)  + dep_ice(i,k)
             dep_rate_snow(i,k) = dep_rate_snow(i,k) + dep_snow(i,k)
           
          ENDIF
       ENDDO
    ENDDO

  CONTAINS

    SUBROUTINE vapor_deposition_ice()
      INTEGER             :: i,k
      INTEGER, SAVE       :: firstcall      
      REAL(wp)            :: q_i,n_i,x_i,d_i,v_i,f_v
      REAL(wp), SAVE      :: c_i        ! coeff for capacity
      REAL(wp), SAVE      :: a_f,b_f    ! coeffs for ventilation

      IF (firstcall.NE.1) THEN
         c_i = 1.0 / ice%cap
         a_f = vent_coeff_a(ice,1)
         b_f = vent_coeff_b(ice,1) * N_sc**n_f / sqrt(nu_l)
         IF (isdebug) THEN
            CALL message(routine,"vapor_deposition_ice:")
            WRITE (txt,'(A,D10.3)') "    a_geo   = ",ice%a_geo ; CALL message(routine,TRIM(txt))        
            WRITE (txt,'(A,D10.3)') "    b_geo   = ",ice%b_geo ; CALL message(routine,TRIM(txt))        
            WRITE (txt,'(A,D10.3)') "    a_vel   = ",ice%a_vel ; CALL message(routine,TRIM(txt))        
            WRITE (txt,'(A,D10.3)') "    b_vel   = ",ice%b_vel ; CALL message(routine,TRIM(txt))        
            WRITE (txt,'(A,D10.3)') "    c_i     = ",c_i ; CALL message(routine,TRIM(txt))        
            WRITE (txt,'(A,D10.3)') "    a_f     = ",a_f ; CALL message(routine,TRIM(txt))        
            WRITE (txt,'(A,D10.3)') "    b_f     = ",b_f ; CALL message(routine,TRIM(txt))        
            firstcall = 1
         ELSEIF (isdebug) THEN
            CALL message(routine, "vapor_deposition_ice")
         ENDIF
      END IF
      
      DO k = kstart,kend
         DO i = istart,iend
            
            IF (q_ice(i,k) == 0.0_wp) THEN
               dep_ice(i,k) = 0.0_wp
            ELSE  
              n_i = n_ice(i,k)                                 
              q_i = q_ice(i,k)                                               
              x_i = ice%meanmass(q_i,n_i)
              D_i = ice%diameter(x_i)
              v_i = ice%velocity(x_i) * rrho_04(i,k)

              !..note that a_f includes more than just ventilation, do never ever set f_v=1
              f_v  = a_f + b_f * sqrt(v_i*d_i)
              f_v  = MAX(f_v,a_f/ice%a_ven) 
              
              dep_ice(i,k) = g_i(i,k) * n_i * c_i * d_i * f_v * s_si(i,k) * dt_local
           ENDIF
        ENDDO
     ENDDO     
   END SUBROUTINE vapor_deposition_ice
   
   SUBROUTINE vapor_deposition_graupel()
     INTEGER             :: i,k
     INTEGER, SAVE       :: firstcall
     REAL(wp)            :: q_g,n_g,x_g,d_g,v_g,f_v
     REAL(wp), SAVE      :: c_g           
     REAL(wp), SAVE      :: a_f,b_f
     
     IF (firstcall.NE.1) THEN
        c_g = 1.0 / graupel%cap
        a_f = vent_coeff_a(graupel,1)
        b_f = vent_coeff_b(graupel,1) * N_sc**n_f / sqrt(nu_l)
        IF (isdebug) THEN
          WRITE(txt,*) "  vapor_deposition_graupel: " ; CALL message(routine,TRIM(txt)) 
          WRITE(txt,'(A,D10.3)') "    a_geo = ",graupel%a_geo ; CALL message(routine,TRIM(txt))   
          WRITE(txt,'(A,D10.3)') "    b_geo = ",graupel%b_geo ; CALL message(routine,TRIM(txt)) 
          WRITE(txt,'(A,D10.3)') "    a_vel = ",graupel%a_vel ; CALL message(routine,TRIM(txt))  
          WRITE(txt,'(A,D10.3)') "    b_vel = ",graupel%b_vel ; CALL message(routine,TRIM(txt))  
          WRITE(txt,'(A,D10.3)') "    c_g   = ",c_g ; CALL message(routine,TRIM(txt))  
          WRITE(txt,'(A,D10.3)') "    a_f   = ",a_f ; CALL message(routine,TRIM(txt))  
          WRITE(txt,'(A,D10.3)') "    b_f   = ",b_f ; CALL message(routine,TRIM(txt))  
        END IF
        firstcall = 1
      ELSEIF (isdebug) THEN
        WRITE(txt,*) "  vapor_deposition_graupel" ; CALL message(routine,TRIM(txt))  
      ENDIF

      DO k = kstart,kend
         DO i = istart,iend

            IF (q_graupel(i,k) == 0.0_wp) THEN
               dep_graupel(i,k) = 0.0_wp
            ELSE
               n_g = n_graupel(i,k)  
               q_g = q_graupel(i,k)  
               x_g = graupel%meanmass(q_g,n_g)
               d_g = graupel%diameter(x_g)
               v_g = graupel%velocity(x_g)
               
               f_v  = a_f + b_f * sqrt(v_g*d_g)
               f_v  = MAX(f_v,a_f/graupel%a_ven)  

               dep_graupel(i,k) = g_i(i,k) * n_g * c_g * d_g * f_v * s_si(i,k) * dt_local
            ENDIF
         ENDDO
      ENDDO
    END SUBROUTINE vapor_deposition_graupel

    SUBROUTINE vapor_deposition_hail()
      INTEGER             :: i,k
      INTEGER, SAVE       :: firstcall
      REAL(wp)            :: q_h,n_h,x_h,d_h,v_h,f_v
      REAL(wp), SAVE      :: c_h               
      REAL(wp), SAVE      :: a_f,b_f 

      IF (firstcall.NE.1) THEN
         c_h = 1.0 / hail%cap
         a_f = vent_coeff_a(hail,1)
         b_f = vent_coeff_b(hail,1) * N_sc**n_f / sqrt(nu_l)
         IF (isdebug) THEN
            WRITE(txt,*) "  vapor_deposition_hail: " 
            WRITE(txt,'(A,D10.3)') "    a_geo = ",hail%a_geo ; CALL message(routine,TRIM(txt))  
            WRITE(txt,'(A,D10.3)') "    b_geo = ",hail%b_geo ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    a_vel = ",hail%a_vel ; CALL message(routine,TRIM(txt))  
            WRITE(txt,'(A,D10.3)') "    b_vel = ",hail%b_vel ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "    c_h   = ",c_h ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "    a_f   = ",a_f ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "    b_f   = ",b_f ; CALL message(routine,TRIM(txt)) 
         END IF
         firstcall = 1
      ELSEIF (isdebug) THEN
         WRITE(txt,*) "  vapor_deposition_hail" ; CALL message(routine,TRIM(txt)) 
      ENDIF
      
      DO k = kstart,kend
         DO i = istart,iend

            IF (q_hail(i,k) == 0.0_wp) THEN
               dep_hail(i,k)   = 0.0_wp
            ELSE
               n_h = n_hail(i,k)                                 
               q_h = q_hail(i,k)            
               x_h = hail%meanmass(q_h,n_h)
               D_h = hail%diameter(x_h)
               v_h = hail%velocity(x_h) * rrho_04(i,k)
                                    
               f_v  = a_f + b_f * sqrt(v_h*d_h)
               f_v  = MAX(f_v,a_f/hail%a_ven)  
               
               dep_hail(i,k) = g_i(i,k) * n_h * c_h * d_h * f_v * s_si(i,k) * dt_local
            ENDIF
         ENDDO
      ENDDO
    END SUBROUTINE vapor_deposition_hail

    SUBROUTINE vapor_deposition_snow()
      REAL(wp)            :: q_s,n_s,x_s,d_s,v_s,f_v
      REAL(wp), SAVE      :: c_s     
      REAL(wp), SAVE      :: a_f,b_f 
      INTEGER             :: i,k
      INTEGER, SAVE       :: firstcall

      IF (firstcall.NE.1) THEN
        c_s = 1.0 / snow%cap
        a_f = vent_coeff_a(snow,1)
        b_f = vent_coeff_b(snow,1) * N_sc**n_f / sqrt(nu_l)
        IF (isdebug) THEN
          WRITE(txt,*) "  vapor_deposition_snow: " 
          WRITE(txt,'(A,D10.3)') "    a_geo = ",snow%a_geo ; CALL message(routine,TRIM(txt))  
          WRITE(txt,'(A,D10.3)') "    b_geo = ",snow%b_geo ; CALL message(routine,TRIM(txt))
          WRITE(txt,'(A,D10.3)') "    a_vel = ",snow%a_vel ; CALL message(routine,TRIM(txt)) 
          WRITE(txt,'(A,D10.3)') "    b_vel = ",snow%b_vel ; CALL message(routine,TRIM(txt)) 
          WRITE(txt,'(A,D10.3)') "    c_s   = ",c_s ; CALL message(routine,TRIM(txt)) 
          WRITE(txt,'(A,D10.3)') "    a_f   = ",a_f ; CALL message(routine,TRIM(txt)) 
          WRITE(txt,'(A,D10.3)') "     b_f  = ",b_f ; CALL message(routine,TRIM(txt)) 
        END IF
        firstcall = 1
      ELSEIF (isdebug) THEN
        WRITE(txt,*) "  vapor_deposition_snow " ; CALL message(routine,TRIM(txt)) 
      ENDIF

      DO k = kstart,kend
         DO i = istart,iend

            IF (q_snow(i,k) == 0.0_wp) THEN
               dep_snow(i,k) = 0.0_wp
            ELSE
               n_s = n_snow(i,k)       
               q_s = q_snow(i,k)       

               x_s = snow%meanmass(q_s,n_s)
               D_s = snow%diameter(x_s)
               v_s = snow%velocity(x_s) * rrho_04(i,k)
               
               f_v  = a_f + b_f * sqrt(D_s*v_s)  
               f_v  = MAX(f_v,a_f/snow%a_ven)   

               dep_snow(i,k) = g_i(i,k) * n_s * c_s * d_s * f_v * s_si(i,k) * dt_local
            ENDIF

         ENDDO
      ENDDO
    END SUBROUTINE vapor_deposition_snow

  END SUBROUTINE vapor_dep_relaxation

  SUBROUTINE rain_freeze_gamlook ()
    !*******************************************************************************
    ! Freezing of raindrops                                                        *
    ! by Uli Blahak                                                                *
    !                                                                              *
    ! incomplete gamma functions are implemented as look-up tables                 *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: fr_q,fr_n,T_a,q_r,x_r,n_r,j_het, &
         &  fr_q_i,fr_n_i,fr_q_g,fr_n_g,fr_q_h,fr_n_h,n_0,lam,xmax_ice,xmax_gr,fr_q_tmp,fr_n_tmp
    REAL(wp), PARAMETER :: a_HET = 6.5d-1 ! Data of Barklie and Gokhale (PK S.350)
    REAL(wp), PARAMETER :: b_HET = 2.0d+2 !         Barklie and Gokhale (PK S.350)

    REAL(wp), SAVE             :: coeff_z
    REAL(wp), SAVE             :: nm1, nm2, nm3, g1, g2
    TYPE(gamlookuptable), SAVE :: ltable1, ltable2, ltable3

    LOGICAL, PARAMETER         :: lclipping = .true.   
    REAL(wp), PARAMETER        :: eps = 1e-15_wp      ! for clipping

    IF (firstcall.NE.1) THEN
      firstcall = 1
      coeff_z = moment_gamma(rain,2)  ! coeff for 2nd moment 
      nm1 = (rain%nu+1.0)/rain%mu
      nm2 = (rain%nu+2.0)/rain%mu
      nm3 = (rain%nu+3.0)/rain%mu
      CALL incgfct_lower_lookupcreate(nm1, ltable1, nlookup, nlookuphr_dummy)
      CALL incgfct_lower_lookupcreate(nm2, ltable2, nlookup, nlookuphr_dummy)
      CALL incgfct_lower_lookupcreate(nm3, ltable3, nlookup, nlookuphr_dummy)      
      g1 = ltable1%igf(ltable1%n) ! ordinary gamma function of nm1 is the last value in the lookup table 1:      
      g2 = ltable2%igf(ltable2%n) ! ordinary gamma function of nm2 is the last value in the lookup table 2:
      IF (isdebug) THEN
        WRITE(txt,*) "rain_freeze_gamlook:" ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "    coeff_z= ",coeff_z ; CALL message(routine,TRIM(txt))
      ENDIF
    ELSE IF (isdebug) THEN
      WRITE(txt,*) "rain_freeze_gamlook" ; CALL message(routine,TRIM(txt))
    ENDIF

    xmax_ice = ( (D_rainfrz_ig/rain%a_geo)**(1.0_wp/rain%b_geo) )**rain%mu
    xmax_gr  = ( (D_rainfrz_gh/rain%a_geo)**(1.0_wp/rain%b_geo) )**rain%mu

    DO k = kstart,kend
       DO i = istart,iend

          T_a = T_p(i,k)
          q_r = q_rain(i,k)
          n_r = n_rain(i,k)

          IF (T_a < T_freeze) THEN
             IF (q_r <= q_crit_fr) THEN
                IF (T_a < T_f) THEN  ! instantaneous freezing below T_f or -40 C
                   fr_q = q_r
                   fr_n = n_r
                   fr_n_i= n_r
                   fr_q_i= q_r
                   fr_n_g= 0.0 ; fr_q_g= 0.0 ; fr_n_h= 0.0 ; fr_q_h= 0.0
                   fr_n_tmp = 1.0 ; fr_q_tmp = 1.0
                ELSE
                   fr_q = 0.0   ; fr_n = 0.0
                   fr_n_i = 0.0 ; fr_q_i = 0.0
                   fr_n_g = 0.0 ; fr_q_g = 0.0
                   fr_n_h = 0.0 ; fr_q_h = 0.0
                   fr_n_tmp = 0.0 ; fr_q_tmp = 0.0
                END IF
             ELSE
                x_r = rain%meanmass(q_r,n_r)
                n_r = q_r / x_r
                IF (T_a < T_f) THEN
                   ! Diesen Zweig koennte man auch weglassen. ist zudem zwar quantitativ richtig, 
                   ! aber nicht konsistent zum
                   ! Grenzfall fuer komplettes Gefrieren der Rechnung im T_a >= T_f - Zweig weiter unten
                   fr_q = q_r                  !  Ausfrieren unterhalb T_f \approx -40 C
                   fr_n = n_r
                   ! Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                   ! oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                   ! bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                   ! und von xmax_gr bis unendlich (--> Hagel).
                   lam = exp( log( g1/g2*x_r) * (-rain%mu) )
                   n_0 = rain%mu * n_r * lam**(nm1) / g1            
                   fr_n_i = n_0/(rain%mu*lam**(nm1)) * incgfct_lower_lookup(lam*xmax_ice,ltable1)
                   fr_q_i = n_0/(rain%mu*lam**(nm2)) * incgfct_lower_lookup(lam*xmax_ice,ltable2)
                   fr_n_g = n_0/(rain%mu*lam**(nm1)) * incgfct_lower_lookup(lam*xmax_gr, ltable1)
                   fr_q_g = n_0/(rain%mu*lam**(nm2)) * incgfct_lower_lookup(lam*xmax_gr, ltable2)
                   
                   fr_n_h = fr_n - fr_n_g
                   fr_q_h = fr_q - fr_q_g
                   fr_n_g = fr_n_g - fr_n_i
                   fr_q_g = fr_q_g - fr_q_i
                   fr_n_tmp = n_r/MAX(fr_n,n_r)
                   fr_q_tmp = q_r/MAX(fr_q,q_r)
                ELSE                           
                   !..heterogeneous freezing
                   j_het = MAX(b_HET * ( EXP( a_HET * (T_3 - T_a)) - 1.0 ),0.0_wp) / rho_w * dt

                   ! Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                   ! oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                   ! bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                   ! und von xmax_gr bis unendlich (--> Hagel).                   
                   IF (j_het >= 1d-20) THEN
                      fr_n  = j_het * q_r
                      fr_q  = j_het * q_r * x_r * coeff_z
                      
                      lam = ( g1 / g2 * x_r)**(-rain%mu)
                      n_0 = rain%mu * n_r * lam**(nm1) / g1
                      fr_n_i = j_het * n_0/(rain%mu*lam**(nm2)) * incgfct_lower_lookup(lam*xmax_ice,ltable2)
                      fr_q_i = j_het * n_0/(rain%mu*lam**(nm3)) * incgfct_lower_lookup(lam*xmax_ice,ltable3)
                      fr_n_g = j_het * n_0/(rain%mu*lam**(nm2)) * incgfct_lower_lookup(lam*xmax_gr, ltable2)
                      fr_q_g = j_het * n_0/(rain%mu*lam**(nm3)) * incgfct_lower_lookup(lam*xmax_gr, ltable3)
                      
                      fr_n_h = fr_n - fr_n_g
                      fr_q_h = fr_q - fr_q_g
                      fr_n_g = fr_n_g - fr_n_i
                      fr_q_g = fr_q_g - fr_q_i
                      fr_n_tmp = n_r/MAX(fr_n,n_r)
                      fr_q_tmp = q_r/MAX(fr_q,q_r)
                   ELSE
                      fr_n= 0.0
                      fr_q= 0.0
                      fr_n_i= 0.0
                      fr_q_i= 0.0
                      fr_n_g= 0.0
                      fr_q_g= 0.0
                      fr_n_h= 0.0
                      fr_q_h= 0.0
                      fr_n_tmp = 0.0
                      fr_q_tmp = 0.0
                   END IF

                END IF
                
                fr_n = fr_n * fr_n_tmp
                fr_q = fr_q * fr_q_tmp
                fr_n_i = fr_n_i * fr_n_tmp
                fr_n_g = fr_n_g * fr_n_tmp
                fr_n_h = fr_n_h * fr_n_tmp
                fr_q_i = fr_q_i * fr_q_tmp
                fr_q_g = fr_q_g * fr_q_tmp
                fr_q_h = fr_q_h * fr_q_tmp                
             END IF

             q_rain(i,k) = q_rain(i,k) - fr_q
             n_rain(i,k) = n_r - fr_n

             IF (ice_typ < 2) THEN 
                ! ohne Hagelklasse,  gefrierender Regen wird Eis oder Graupel
                q_ice(i,k) = q_ice(i,k)  + fr_q_i
                n_ice(i,k) = n_ice(i,k)  + fr_n_i
                q_graupel(i,k) = q_graupel(i,k)  + fr_q_h + fr_q_g
                n_graupel(i,k) = n_graupel(i,k)  + fr_n_h + fr_n_g
             ELSE
                ! mit Hagelklasse, gefrierender Regen wird Eis, Graupel oder Hagel
                q_ice(i,k) = q_ice(i,k)  + fr_q_i
                n_ice(i,k) = n_ice(i,k)  + fr_n_i
                q_graupel(i,k) = q_graupel(i,k)  + fr_q_g
                n_graupel(i,k) = n_graupel(i,k)  + fr_n_g
                q_hail(i,k) = q_hail(i,k)  + fr_q_h
                n_hail(i,k) = n_hail(i,k)  + fr_n_h
             ENDIF

             ! clipping of small negatives is necessary here
             if (lclipping) then
                IF (q_rain(i,k) < 0.0 .and. abs(q_rain(i,k)) < eps) q_rain(i,k) = 0.0_wp
                IF (n_rain(i,k) < 0.0 .and. abs(n_rain(i,k)) < eps) n_rain(i,k) = 0.0_wp
                IF (q_graupel(i,k) < 0.0 .and. abs(q_graupel(i,k)) < eps) q_graupel(i,k) = 0.0_wp
                IF (n_graupel(i,k) < 0.0 .and. abs(q_graupel(i,k)) < eps) n_graupel(i,k) = 0.0_wp
                IF (q_hail(i,k) < 0.0 .and. abs(q_hail(i,k)) < eps) q_hail(i,k) = 0.0_wp
                IF (n_hail(i,k) < 0.0 .and. abs(n_hail(i,k)) < eps) n_hail(i,k) = 0.0_wp
             end if

          END IF
       END DO
    END DO
  END SUBROUTINE rain_freeze_gamlook

  SUBROUTINE ice_selfcollection()
    !*******************************************************************************
    ! selfcollection of ice crystals, see SB2006 or Seifert (2002)                 *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a       
    REAL(wp)            :: q_i,n_i,x_i,d_i,v_i,e_coll,x_conv_ii
    REAL(wp)            :: self_n,self_q
    REAL(wp)            :: delta_n_11,delta_n_12,delta_n_22
    REAL(wp)            :: delta_q_11,delta_q_12,delta_q_22
    REAL(wp)            :: theta_n_11,theta_n_12,theta_n_22
    REAL(wp)            :: theta_q_11,theta_q_12,theta_q_22
    REAL(wp),SAVE       :: delta_n,delta_q
    REAL(wp),SAVE       :: theta_n,theta_q

    IF (isdebug) THEN
      WRITE(txt,*) "ice_selfcollection " ; CALL message(routine,TRIM(txt)) 
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_11 = coll_delta_11(ice,ice,0)
      delta_n_12 = coll_delta_12(ice,ice,0)
      delta_n_22 = coll_delta_22(ice,ice,0)
      delta_q_11 = coll_delta_11(ice,ice,0) 
      delta_q_12 = coll_delta_12(ice,ice,1)
      delta_q_22 = coll_delta_22(ice,ice,1)

      theta_n_11 = coll_theta_11(ice,ice,0)
      theta_n_12 = coll_theta_12(ice,ice,0)
      theta_n_22 = coll_theta_22(ice,ice,0)
      theta_q_11 = coll_theta_11(ice,ice,0)
      theta_q_12 = coll_theta_12(ice,ice,1)
      theta_q_22 = coll_theta_22(ice,ice,1)

      delta_n = delta_n_11 + delta_n_12 + delta_n_22
      delta_q = delta_q_11 + delta_q_12 + delta_q_22
      theta_n = theta_n_11 - theta_n_12 + theta_n_22
      theta_q = theta_q_11 - theta_q_12 + theta_q_22

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "    a_ice      = ",ice%a_geo ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "    b_ice      = ",ice%b_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    alf_ice    = ",ice%a_vel ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "    bet_ice    = ",ice%b_vel ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "    delta_n_11 = ",delta_n_11 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_12 = ",delta_n_12 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_22 = ",delta_n_22 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n    = ",delta_n    ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "    theta_n_11 = ",theta_n_11 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_12 = ",theta_n_12 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_22 = ",theta_n_22 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n    = ",theta_n    ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_11 = ",delta_q_11 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_12 = ",delta_q_12 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_22 = ",delta_q_22 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q    = ",delta_q    ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_11 = ",theta_q_11 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_12 = ",theta_q_12 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_22 = ",theta_q_22 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q    = ",theta_q    ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ENDIF

    x_conv_ii = (D_conv_ii/snow%a_geo)**(1./snow%b_geo)

    DO k = kstart,kend
       DO i = istart,iend
          
          q_i = q_ice(i,k)                                  
          n_i = n_ice(i,k)           

          x_i = ice%meanmass(q_i,n_i)
          D_i = ice%diameter(x_i)                      

          IF ( n_i > 0.0_wp .AND. q_i > q_crit_ii .AND. D_i > D_crit_ii ) THEN

             T_a = T_p(i,k)

             !.. Temperaturabhaengige Efficiency nach Cotton et al. (1986) 
             !   (siehe auch Straka, 1989; S. 53)
             e_coll = MIN(10**(0.035*(T_a-T_3)-0.7),0.2_wp)
             
             v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,k)  
             
             self_n = pi4 * e_coll * delta_n * n_i * n_i * D_i * D_i &
                  & * sqrt( theta_n * v_i * v_i + 2.0*ice_s_vel**2 ) * dt
             
             self_q = pi4 * e_coll * delta_q * n_i * q_i * D_i * D_i &
                  & * sqrt( theta_q * v_i * v_i + 2.0*ice_s_vel**2 ) * dt
             
             self_q = MIN(self_q,q_i)
             self_n = MIN(MIN(self_n,self_q/x_conv_ii),n_i)
             
             q_ice(i,k)  = q_ice(i,k)  - self_q
             q_snow(i,k) = q_snow(i,k) + self_q
             
             n_ice(i,k)  = n_ice(i,k)  - self_n
             n_snow(i,k) = n_snow(i,k) + self_n / 2.0
             
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE ice_selfcollection
 
  SUBROUTINE snow_selfcollection()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a             
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s,e_coll
    REAL(wp)            :: self_n
    REAL(wp)            :: delta_n_11,delta_n_12
    REAL(wp)            :: theta_n_11,theta_n_12
    REAL(wp),SAVE       :: delta_n
    REAL(wp),SAVE       :: theta_n

    IF (isdebug) THEN
      WRITE(txt,*) "snow_selfcollection " ; CALL message(routine,TRIM(txt)) 
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_11 = coll_delta_11(snow,snow,0)
      delta_n_12 = coll_delta_12(snow,snow,0)
      theta_n_11 = coll_theta_11(snow,snow,0)
      theta_n_12 = coll_theta_12(snow,snow,0)

      delta_n = (2.0*delta_n_11 + delta_n_12)
      theta_n = (2.0*theta_n_11 - theta_n_12)

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "    a_snow     = ",snow%a_geo ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "    b_snow     = ",snow%b_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    alf_snow   = ",snow%a_vel ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    bet_snow   = ",snow%b_vel ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "    delta_n_11 = ",delta_n_11  ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_12 = ",delta_n_12  ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n    = ",delta_n     ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_11 = ",theta_n_11  ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_12 = ",theta_n_12  ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n    = ",theta_n     ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          q_s = q_snow(i,k)   
          T_a = T_p(i,k)
 
          IF ( q_s > q_crit ) THEN

             !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
             e_coll = MAX(0.1_wp,MIN(EXP(0.09*(T_a-T_3)),1.0_wp))

             n_s = n_snow(i,k)     
             x_s = snow%meanmass(q_s,n_s)
             D_s = snow%diameter(x_s)
             v_s = snow%velocity(x_s) * rrho_04(i,k)
                
             self_n = pi8 * e_coll * n_s * n_s * delta_n * D_s * D_s * &
                  &          sqrt( theta_n * v_s * v_s + 2.0 * snow_s_vel**2 ) * dt

             self_n = MIN(self_n,n_s)

             n_snow(i,k) = n_snow(i,k) - self_n

         ENDIF
      ENDDO
   ENDDO

  END SUBROUTINE snow_selfcollection

  SUBROUTINE snow_melting()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s,T_a,e_a
    REAL(wp)            :: melt,melt_v,melt_h,melt_n,melt_q
    REAL(wp)            :: fh_q,fv_q
    REAL(wp), SAVE      :: a_vent,b_vent

    IF (isdebug) THEN
      WRITE(txt,*) "snow_melting " ; CALL message(routine,TRIM(txt))
    END IF

    IF (firstcall.NE.1) THEN
      a_vent = vent_coeff_a(snow,1)
      b_vent = vent_coeff_b(snow,1) * N_sc**n_f / sqrt(nu_l)

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "    a_geo  = ",snow%a_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    b_geo  = ",snow%b_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    a_vel  = ",snow%a_vel ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    b_vel  = ",snow%b_vel ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    a_ven  = ",snow%a_ven ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    b_ven  = ",snow%b_ven ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    a_vent = ",a_vent ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    b_vent = ",b_vent ; CALL message(routine,TRIM(txt))
      ENDIF
      firstcall = 1
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend
          
          T_a = T_p(i,k)
          q_s = q_snow(i,k)                   
          IF (T_a > T_3 .AND. q_s > 0.0_wp) THEN
            e_a = e_ws(T_a)            ! saturation pressure
            n_s = n_snow(i,k)                               

            x_s = snow%meanmass(q_s,n_s)
            D_s = snow%diameter(x_s)
            v_s = snow%velocity(x_s) * rrho_04(i,k)

            fv_q = a_vent + b_vent * sqrt(v_s*D_s)

            ! UB: Based on Rasmussen and Heymsfield (1987) the ratio fh_q / fv_q is approx. 1.05
            !     for a wide temperature- and pressure range:
            fh_q = 1.05 * fv_q

            melt   = 2.0*pi/L_ew * D_s * n_s * dt

            melt_h = melt * K_T * (T_a - T_3)
            melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)
            melt_q = (melt_h * fh_q + melt_v * fv_q)

            ! UB: for melt_n we assume that x_s is constant during melting
            melt_n = MIN(MAX( (melt_q - q_s) / x_s + n_s, 0.0_wp), n_s)
 
            melt_q = MIN(q_s,MAX(melt_q,0.0_wp))
            melt_n = MIN(n_s,MAX(melt_n,0.0_wp))

            ! UB: snow melts instantaneously at 10 C
            IF (T_a - T_3 > 10.0_wp) THEN
              melt_q = q_s  
              melt_n = n_s
            ENDIF

            q_snow(i,k) = q_snow(i,k) - melt_q
            q_rain(i,k) = q_rain(i,k) + melt_q

            n_snow(i,k) = n_snow(i,k) - melt_n
            n_rain(i,k) = n_rain(i,k) + melt_n

            n_snow(i,k) = MAX(n_snow(i,k), q_snow(i,k)/snow%x_max)

         ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE snow_melting

  SUBROUTINE graupel_snow_collection()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************

    ! Locale Variablen 
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s
    REAL(wp)            :: coll_n,coll_q,e_coll
    REAL(wp), SAVE      :: delta_n_gg,delta_n_gs,delta_n_ss
    REAL(wp), SAVE      :: delta_q_gg,delta_q_gs,delta_q_ss
    REAL(wp), SAVE      :: theta_n_gg,theta_n_gs,theta_n_ss
    REAL(wp), SAVE      :: theta_q_gg,theta_q_gs,theta_q_ss

    IF (isdebug) THEN
      WRITE(txt,*) "graupel_snow_collection" ; CALL message(routine,TRIM(txt))
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,snow,0)
      delta_n_gs = coll_delta_12(graupel,snow,0)
      delta_n_ss = coll_delta_22(graupel,snow,0)
      delta_q_gg = coll_delta_11(graupel,snow,0) 
      delta_q_gs = coll_delta_12(graupel,snow,1)
      delta_q_ss = coll_delta_22(graupel,snow,1)

      theta_n_gg = coll_theta_11(graupel,snow,0)
      theta_n_gs = coll_theta_12(graupel,snow,0)
      theta_n_ss = coll_theta_22(graupel,snow,0)
      theta_q_gg = coll_theta_11(graupel,snow,0)
      theta_q_gs = coll_theta_12(graupel,snow,1)
      theta_q_ss = coll_theta_22(graupel,snow,1)

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "    delta_n_gg = ",delta_n_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_gs = ",delta_n_gs ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_ss = ",delta_n_ss ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_gg = ",theta_n_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_gs = ",theta_n_gs ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_ss = ",theta_n_ss ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_gg = ",delta_q_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_gs = ",delta_q_gs ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_ss = ",delta_q_ss ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_gg = ",theta_q_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_gs = ",theta_q_gs ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_ss = ",theta_q_ss ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          q_s = q_snow(i,k)    
          q_g = q_graupel(i,k) 

          IF (q_s > q_crit .AND. q_g > q_crit) THEN
            T_a = T_p(i,k) 

            !.. sticking efficiency of Lin (1983)
            e_coll = min(exp(0.09*(T_a-T_3)),1.0_wp)

            n_s = n_snow(i,k)                                      
            n_g = n_graupel(i,k)                                   

            x_g = graupel%meanmass(q_g,n_g)
            d_g = graupel%diameter(x_g)
            v_g = graupel%velocity(x_g) * rrho_04(i,k)

            x_s = snow%meanmass(q_s,n_s)
            d_s = snow%diameter(x_s)
            v_s = snow%velocity(x_s) * rrho_04(i,k)

            coll_n = pi4 * n_g * n_s * e_coll * dt & 
                 & *     (delta_n_gg * D_g**2 + delta_n_gs * D_g*D_s + delta_n_ss * D_s**2) &
                 & * sqrt(theta_n_gg * v_g**2 - theta_n_gs * v_g*v_s + theta_n_ss * v_s**2  &
                 &       +snow_s_vel**2)

            coll_q = pi4 * n_g * q_s * e_coll * dt & 
                 & *     (delta_q_gg * D_g**2 + delta_q_gs * D_g*D_s + delta_q_ss * D_s**2) &
                 & * sqrt(theta_q_gg * v_g**2 - theta_q_gs * v_g*v_s + theta_q_ss * v_s**2  &
                 &       +snow_s_vel**2)

            coll_n = MIN(n_s,coll_n)
            coll_q = MIN(q_s,coll_q)

            q_graupel(i,k) = q_graupel(i,k) + coll_q
            q_snow(i,k)    = q_snow(i,k)    - coll_q
            n_snow(i,k)    = n_snow(i,k)    - coll_n

         ENDIF
      ENDDO
   ENDDO

  END SUBROUTINE graupel_snow_collection

  SUBROUTINE hail_snow_collection()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s
    REAL(wp)            :: coll_n,coll_q,e_coll
    REAL(wp), SAVE      :: delta_n_hh,delta_n_hs,delta_n_ss
    REAL(wp), SAVE      :: delta_q_hh,delta_q_hs,delta_q_ss
    REAL(wp), SAVE      :: theta_n_hh,theta_n_hs,theta_n_ss
    REAL(wp), SAVE      :: theta_q_hh,theta_q_hs,theta_q_ss

    IF (isdebug) THEN
      WRITE(txt,*) "hail_snow_collection" ; CALL message(routine,TRIM(txt))
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,snow,0)
      delta_n_hs = coll_delta_12(hail,snow,0)
      delta_n_ss = coll_delta_22(hail,snow,0)
      delta_q_hh = coll_delta_11(hail,snow,0) 
      delta_q_hs = coll_delta_12(hail,snow,1)
      delta_q_ss = coll_delta_22(hail,snow,1)

      theta_n_hh = coll_theta_11(hail,snow,0)
      theta_n_hs = coll_theta_12(hail,snow,0)
      theta_n_ss = coll_theta_22(hail,snow,0)
      theta_q_hh = coll_theta_11(hail,snow,0)
      theta_q_hs = coll_theta_12(hail,snow,1)
      theta_q_ss = coll_theta_22(hail,snow,1)

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "    delta_n_hh = ",delta_n_hh ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_hs = ",delta_n_hs ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_ss = ",delta_n_ss ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_hh = ",theta_n_hh ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_hs = ",theta_n_hs ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_ss = ",theta_n_ss ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_hh = ",delta_q_hh ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_hs = ",delta_q_hs ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_ss = ",delta_q_ss ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_hh = ",theta_q_hh ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_hs = ",theta_q_hs ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_ss = ",theta_q_ss ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          q_s = q_snow(i,k) 
          q_h = q_hail(i,k) 

          IF (q_s > q_crit .AND. q_h > q_crit) THEN
            T_a = T_p(i,k)

            !.. sticking efficiency of Lin (1983)
            e_coll = min(exp(0.09*(T_a-T_3)),1.0_wp)

            n_s = n_snow(i,k)                                  
            n_h = n_hail(i,k)        

            x_s = snow%meanmass(q_s,n_s)
            D_s = snow%diameter(x_s)
            v_s = snow%velocity(x_s) * rrho_04(i,k)

            x_h = hail%meanmass(q_h,n_h)
            D_h = hail%diameter(x_h)
            v_h = hail%velocity(x_h) * rrho_04(i,k)
                          
            coll_n = pi4 * n_h * n_s * e_coll * dt & 
                 & *     (delta_n_hh * D_h**2 + delta_n_hs * D_h*D_s + delta_n_ss * D_s**2) &
                 & * sqrt(theta_n_hh * v_h**2 - theta_n_hs * v_h*v_s + theta_n_ss * v_s**2  &
                 &       +snow_s_vel**2)

            coll_q = pi4 * n_h * q_s * e_coll * dt & 
                 & *     (delta_q_hh * D_h**2 + delta_q_hs * D_h*D_s + delta_q_ss * D_s**2) &
                 & * sqrt(theta_q_hh * v_h**2 - theta_q_hs * v_h*v_s + theta_q_ss * v_s**2  &
                 &       +snow_s_vel**2)

            coll_n = MIN(n_s,coll_n)
            coll_q = MIN(q_s,coll_q)

            q_hail(i,k) = q_hail(i,k) + coll_q
            q_snow(i,k) = q_snow(i,k) - coll_q
            n_snow(i,k) = n_snow(i,k) - coll_n

          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE hail_snow_collection

  SUBROUTINE graupel_ice_collection()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
    REAL(wp)            :: coll_n,coll_q,e_coll
    REAL(wp), SAVE      :: delta_n_gg,delta_n_gi,delta_n_ii
    REAL(wp), SAVE      :: delta_q_gg,delta_q_gi,delta_q_ii
    REAL(wp), SAVE      :: theta_n_gg,theta_n_gi,theta_n_ii
    REAL(wp), SAVE      :: theta_q_gg,theta_q_gi,theta_q_ii

    IF (isdebug) THEN
      WRITE(txt,*) " graupel_ice_collection" ; CALL message(routine,TRIM(txt))
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,ice,0)
      delta_n_gi = coll_delta_12(graupel,ice,0)
      delta_n_ii = coll_delta_22(graupel,ice,0)
      delta_q_gg = coll_delta_11(graupel,ice,0) 
      delta_q_gi = coll_delta_12(graupel,ice,1)
      delta_q_ii = coll_delta_22(graupel,ice,1)

      theta_n_gg = coll_theta_11(graupel,ice,0)
      theta_n_gi = coll_theta_12(graupel,ice,0)
      theta_n_ii = coll_theta_22(graupel,ice,0)
      theta_q_gg = coll_theta_11(graupel,ice,0)
      theta_q_gi = coll_theta_12(graupel,ice,1)
      theta_q_ii = coll_theta_22(graupel,ice,1)

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "   a_graupel   = ",graupel%a_geo ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "   b_graupel   = ",graupel%b_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   alf_graupel = ",graupel%a_vel ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   bet_graupel = ",graupel%b_vel ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   a_ice       = ",ice%a_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   b_ice       = ",ice%b_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   alf_ice     = ",ice%a_vel ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   bet_ice     = ",ice%b_vel ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_n_gg  = ",delta_n_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_n_gi  = ",delta_n_gi ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_n_ii  = ",delta_n_ii ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_n_gg  = ",theta_n_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_n_gi  = ",theta_n_gi ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_n_ii  = ",theta_n_ii ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_q_gg  = ",delta_q_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_q_gi  = ",delta_q_gi ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_q_ii  = ",delta_q_ii ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_q_gg  = ",theta_q_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_q_gi  = ",theta_q_gi ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_q_ii  = ",theta_q_ii ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          q_i = q_ice(i,k)        
          q_g = q_graupel(i,k)    
          IF (q_i > q_crit .AND. q_g > q_crit) THEN
            T_a = T_p(i,k)

            !.. sticking efficiency of Lin (1983)
            e_coll = min(exp(0.09*(T_a-T_3)),1.0_wp)

            n_i = n_ice(i,k)                                        
            n_g = n_graupel(i,k)      

            x_g = graupel%meanmass(q_g,n_g)
            d_g = graupel%diameter(x_g)
            v_g = graupel%velocity(x_g) * rrho_04(i,k)

            x_i = ice%meanmass(q_i,n_i)
            d_i = ice%diameter(x_i)
            v_i = ice%velocity(x_i) * rrho_04(i,k)
                              
            coll_n = pi4 * n_g * n_i * e_coll * dt & 
                 & *     (delta_n_gg * D_g**2 + delta_n_gi * D_g*D_i + delta_n_ii * D_i**2) &
                 & * sqrt(theta_n_gg * v_g**2 - theta_n_gi * v_g*v_i + theta_n_ii * v_i**2)

            coll_q = pi4 * n_g * q_i * e_coll * dt & 
                 & *     (delta_q_gg * D_g**2 + delta_q_gi * D_g*D_i + delta_q_ii * D_i**2) &
                 & * sqrt(theta_q_gg * v_g**2 - theta_q_gi * v_g*v_i + theta_q_ii * v_i**2)

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            q_graupel(i,k) = q_graupel(i,k) + coll_q
            q_ice(i,k)     = q_ice(i,k)     - coll_q
            n_ice(i,k)     = n_ice(i,k)     - coll_n

          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE graupel_ice_collection

  SUBROUTINE hail_ice_collection()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h
    REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
    REAL(wp)            :: coll_n,coll_q,e_coll
    REAL(wp), SAVE      :: delta_n_hh,delta_n_hi,delta_n_ii
    REAL(wp), SAVE      :: delta_q_hh,delta_q_hi,delta_q_ii
    REAL(wp), SAVE      :: theta_n_hh,theta_n_hi,theta_n_ii
    REAL(wp), SAVE      :: theta_q_hh,theta_q_hi,theta_q_ii

    IF (isdebug) THEN
      WRITE(txt,*) "hail_ice_collection" ; CALL message(routine,TRIM(txt))
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,ice,0)
      delta_n_hi = coll_delta_12(hail,ice,0)
      delta_n_ii = coll_delta_22(hail,ice,0)
      delta_q_hh = coll_delta_11(hail,ice,0) 
      delta_q_hi = coll_delta_12(hail,ice,1)
      delta_q_ii = coll_delta_22(hail,ice,1)

      theta_n_hh = coll_theta_11(hail,ice,0)
      theta_n_hi = coll_theta_12(hail,ice,0)
      theta_n_ii = coll_theta_22(hail,ice,0)
      theta_q_hh = coll_theta_11(hail,ice,0)
      theta_q_hi = coll_theta_12(hail,ice,1)
      theta_q_ii = coll_theta_22(hail,ice,1)

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "    delta_n_hh = ",delta_n_hh ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_hi = ",delta_n_hi ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_ii = ",delta_n_ii ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_hh = ",theta_n_hh ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_hi = ",theta_n_hi ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_ii = ",theta_n_ii ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_hh = ",delta_q_hh ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_hi = ",delta_q_hi ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_ii = ",delta_q_ii ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_hh = ",theta_q_hh ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_hi = ",theta_q_hi ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_ii = ",theta_q_ii ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          q_i = q_ice(i,k)    
          q_h = q_hail(i,k)   
          IF (q_i > q_crit .AND. q_h > q_crit) THEN
            T_a = T_p(i,k)

            !.. sticking efficiency of Lin (1983)
            e_coll = min(exp(0.09*(T_a-T_3)),1.0_wp)

            n_i = n_ice(i,k)        
            n_h = n_hail(i,k)       

            x_h = hail%meanmass(q_h,n_h)
            D_h = hail%diameter(x_h)
            v_h = hail%velocity(x_h) * rrho_04(i,k)

            x_i = ice%meanmass(q_i,n_i)
            D_i = ice%diameter(x_i)
            v_i = ice%velocity(x_i) * rrho_04(i,k)

            coll_n = pi4 * n_h * n_i * e_coll * dt & 
                 & *     (delta_n_hh * D_h**2 + delta_n_hi * D_h*D_i + delta_n_ii * D_i**2) &
                 & * sqrt(theta_n_hh * v_h**2 - theta_n_hi * v_h*v_i + theta_n_ii * v_i**2)

            coll_q = pi4 * n_h * q_i * e_coll * dt & 
                 & *     (delta_q_hh * D_h**2 + delta_q_hi * D_h*D_i + delta_q_ii * D_i**2) &
                 & * sqrt(theta_q_hh * v_h**2 - theta_q_hi * v_h*v_i + theta_q_ii * v_i**2)

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            q_hail(i,k) = q_hail(i,k) + coll_q
            q_ice(i,k)  = q_ice(i,k)  - coll_q
            n_ice(i,k)  = n_ice(i,k)  - coll_n

          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE hail_ice_collection

  SUBROUTINE snow_ice_collection()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s
    REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
    REAL(wp)            :: coll_n,coll_q,e_coll
    REAL(wp), SAVE      :: delta_n_ss,delta_n_si,delta_n_ii
    REAL(wp), SAVE      :: delta_q_ss,delta_q_si,delta_q_ii
    REAL(wp), SAVE      :: theta_n_ss,theta_n_si,theta_n_ii
    REAL(wp), SAVE      :: theta_q_ss,theta_q_si,theta_q_ii

    IF (isdebug) THEN
      WRITE(txt,*) "snow_ice_collection" ; CALL message(routine,TRIM(txt))
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_ss = coll_delta_11(snow,ice,0)
      delta_n_si = coll_delta_12(snow,ice,0)
      delta_n_ii = coll_delta_22(snow,ice,0)
      delta_q_ss = coll_delta_11(snow,ice,0) 
      delta_q_si = coll_delta_12(snow,ice,1)
      delta_q_ii = coll_delta_22(snow,ice,1)

      theta_n_ss = coll_theta_11(snow,ice,0)
      theta_n_si = coll_theta_12(snow,ice,0)
      theta_n_ii = coll_theta_22(snow,ice,0)
      theta_q_ss = coll_theta_11(snow,ice,0)
      theta_q_si = coll_theta_12(snow,ice,1)
      theta_q_ii = coll_theta_22(snow,ice,1)

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "   delta_n_ss = ",delta_n_ss ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_n_si = ",delta_n_si ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_n_ii = ",delta_n_ii ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_n_ss = ",theta_n_ss ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_n_si = ",theta_n_si ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_n_ii = ",theta_n_ii ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_q_ss = ",delta_q_ss ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_q_si = ",delta_q_si ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   delta_q_ii = ",delta_q_ii ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_q_ss = ",theta_q_ss ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_q_si = ",theta_q_si ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   theta_q_ii = ",theta_q_ii ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          q_i = q_ice(i,k)
          q_s = q_snow(i,k)
          IF (q_i > q_crit .AND. q_s > q_crit) THEN
            T_a = T_p(i,k) 

            !.. sticking efficiency of Lin (1983)
            e_coll = max(0.1_wp,min(exp(0.09*(T_a-T_3)),1.0_wp))

            n_i = n_ice(i,k)   
            n_s = n_snow(i,k)  

            x_i = ice%meanmass(q_i,n_i)
            d_i = ice%diameter(x_i)
            v_i = ice%velocity(x_i) * rrho_04(i,k)

            x_s = snow%meanmass(q_s,n_s)
            d_s = snow%diameter(x_s)
            v_s = snow%velocity(x_s) * rrho_04(i,k)

            coll_n = pi4 * n_s * n_i * e_coll * dt & 
                 & *     (delta_n_ss * D_s**2 + delta_n_si * D_s*D_i + delta_n_ii * D_i**2) &
                 & * sqrt(theta_n_ss * v_s**2 - theta_n_si * v_s*v_i + theta_n_ii * v_i**2  &
                 &       +snow_s_vel**2 + ice_s_vel**2)

            coll_q = pi4 * n_s * q_i * e_coll * dt & 
                 & *     (delta_q_ss * D_s**2 + delta_q_si * D_s*D_i + delta_q_ii * D_i**2) &
                 & * sqrt(theta_q_ss * v_s**2 - theta_q_si * v_s*v_i + theta_q_ii * v_i**2  &
                 &       +snow_s_vel**2 + ice_s_vel**2)

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            q_snow(i,k) = q_snow(i,k) + coll_q
            q_ice(i,k)  = q_ice(i,k)  - coll_q
            n_ice(i,k)  = n_ice(i,k)  - coll_n
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE snow_ice_collection

  SUBROUTINE graupel_selfcollection()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: self_n
    REAL(wp)            :: delta_n_11,delta_n_12
    REAL(wp)            :: theta_n_11,theta_n_12
    REAL(wp)            :: delta_n, theta_n
    REAL(wp),SAVE       :: coll_n

    IF (isdebug) THEN
      WRITE(txt,*) "graupel_selfcollection" ; CALL message(routine,TRIM(txt)) 
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_11 = coll_delta_11(graupel,graupel,0)
      delta_n_12 = coll_delta_12(graupel,graupel,0)
      theta_n_11 = coll_theta_11(graupel,graupel,0)
      theta_n_12 = coll_theta_12(graupel,graupel,0)

      delta_n = (2.0*delta_n_11 + delta_n_12)
      theta_n = (2.0*theta_n_11 - theta_n_12)**0.5
      coll_n  = pi8 * delta_n * theta_n

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "    delta_n_11 = ",delta_n_11 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_12 = ",delta_n_12 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n    = ",delta_n    ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_11 = ",theta_n_11 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_12 = ",theta_n_12 ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n    = ",theta_n    ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    coll_n     = ",coll_n     ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          q_g = q_graupel(i,k)
          IF ( q_g > q_crit ) THEN

             n_g = n_graupel(i,k)     
             x_g = graupel%meanmass(q_g,n_g)
             d_g = graupel%diameter(x_g)
             v_g = graupel%velocity(x_g) * rrho_04(i,k)
             
             self_n = coll_n * n_g**2 * D_g**2 * v_g * dt

             ! sticking efficiency does only distinguish dry and wet based on T_3
             IF (T_p(i,k) > T_3) THEN
                self_n = self_n * ecoll_gg_wet
             ELSE
                self_n = self_n * ecoll_gg
             END IF

             self_n = MIN(self_n,n_g)

             n_graupel(i,k) = n_graupel(i,k) - self_n
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE graupel_selfcollection

  SUBROUTINE ice_melting ()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    INTEGER             :: i,k
    REAL(wp)            :: q_i,x_i,n_i
    REAL(wp)            :: melt_q,melt_n

    IF (isdebug) THEN
      WRITE(txt,*) "ice_melting" ; CALL message(routine,TRIM(txt)) 
    END IF

    DO k = kstart,kend
       DO i = istart,iend
 
          q_i = q_ice(i,k)

          IF (T_p(i,k) > T_3 .AND. q_i > 0.0) THEN    

            n_i = n_ice(i,k)
            x_i = ice%meanmass(q_i,n_i)

            ! complete melting within this time step
            melt_q = q_i         
            melt_n = n_i            
            q_ice(i,k) = 0.0_wp
            n_ice(i,k) = 0.0_wp

            ! ice either melts into cloud droplets or rain depending on x_i
            IF (x_i > cloud%x_max) THEN
               q_rain(i,k)  = q_rain(i,k)  + melt_q
               n_rain(i,k)  = n_rain(i,k)  + melt_n
            ELSE
               q_cloud(i,k) = q_cloud(i,k) + melt_q
               n_cloud(i,k) = n_cloud(i,k) + melt_n
            ENDIF

          END IF
       END DO
    END DO
  END SUBROUTINE ice_melting

  SUBROUTINE graupel_cloud_riming()
    !*******************************************************************************
    ! see SB2006                                                                   *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: q_c,n_c,x_c,d_c,v_c,x_coll_c
    REAL(wp)            :: rime_n,rime_q,e_coll_n
    REAL(wp)            :: melt_n,melt_q,e_coll_q
    REAL(wp)            :: shed_n,shed_q,x_shed
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp), SAVE      :: delta_n_gg,delta_n_gc,delta_n_cc
    REAL(wp), SAVE      :: delta_q_gg,delta_q_gc,delta_q_cc
    REAL(wp), SAVE      :: theta_n_gg,theta_n_gc,theta_n_cc
    REAL(wp), SAVE      :: theta_q_gg,theta_q_gc,theta_q_cc
    REAL(wp)            :: const1,const2,const3,const4

    IF (isdebug) THEN
       WRITE(txt,*) "graupel_cloud_riming" ; CALL message(routine,TRIM(txt)) 
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,cloud,0)
      delta_n_gc = coll_delta_12(graupel,cloud,0)
      delta_n_cc = coll_delta_22(graupel,cloud,0)
      delta_q_gg = coll_delta_11(graupel,cloud,0) 
      delta_q_gc = coll_delta_12(graupel,cloud,1)
      delta_q_cc = coll_delta_22(graupel,cloud,1)

      theta_n_gg = coll_theta_11(graupel,cloud,0)
      theta_n_gc = coll_theta_12(graupel,cloud,0)
      theta_n_cc = coll_theta_22(graupel,cloud,0)
      theta_q_gg = coll_theta_11(graupel,cloud,0)
      theta_q_gc = coll_theta_12(graupel,cloud,1)
      theta_q_cc = coll_theta_22(graupel,cloud,1)

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "    delta_n_gg  = ",delta_n_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_gc  = ",delta_n_gc ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_n_cc  = ",delta_n_cc ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_gg  = ",theta_n_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_gc  = ",theta_n_gc ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_n_cc  = ",theta_n_cc ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_gg  = ",delta_q_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_gc  = ",delta_q_gc ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    delta_q_cc  = ",delta_q_cc ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_gg  = ",theta_q_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_gc  = ",theta_q_gc ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    theta_q_cc  = ",theta_q_cc ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ENDIF

    x_shed = 4./3.*pi * rho_w * r_shedding**3     !..mean mass of shedding drops

    x_coll_c = (D_coll_c/cloud%a_geo)**3          !..lower threshold for collection

    const1 = ecoll_gc/(D_coll_c - D_crit_c)    
    const2 = 1.0/(T_mult_opt - T_mult_min)
    const3 = 1.0/(T_mult_opt - T_mult_max)
    const4 = c_w / L_ew

    DO k = kstart,kend
       DO i = istart,iend

          q_c = q_cloud(i,k)                                      
          q_g = q_graupel(i,k)                                    
          n_c = n_cloud(i,k)                                      
          n_g = n_graupel(i,k)       
                             
          x_g = graupel%meanmass(q_g,n_g)
          D_g = graupel%diameter(x_g)
          x_c = cloud%meanmass(q_c,n_c)
          D_c = cloud%diameter(x_c)

          T_a = T_p(i,k)
          IF (q_c > q_crit_c .AND. q_g > q_crit_gc .AND. D_g > D_crit_gc .AND. &
               & D_c > D_crit_c) THEN

             v_g = graupel%velocity(x_g) * rrho_04(i,k)
             v_c = cloud%velocity(x_c)   * rrho_04(i,k)

             e_coll_n = MIN(ecoll_gc, MAX(const1*(D_c - D_crit_c),ecoll_min))
             e_coll_q = e_coll_n

             rime_n = pi4 * e_coll_n * n_g * n_c * dt & 
                  & *     (delta_n_gg * D_g*D_g + delta_n_gc * D_g*D_c + delta_n_cc * D_c*D_c) &
                  & * sqrt(theta_n_gg * v_g*v_g - theta_n_gc * v_g*v_c + theta_n_cc * v_c*v_c)

             rime_q = pi4 * e_coll_q * n_g * q_c * dt & 
                  & *     (delta_q_gg * D_g*D_g + delta_q_gc * D_g*D_c + delta_q_cc * D_c*D_c) &
                  & * sqrt(theta_q_gg * v_g*v_g - theta_q_gc * v_g*v_c + theta_q_cc * v_c*v_c)

             rime_q = MIN(q_c,rime_q)
             rime_n = MIN(n_c,rime_n)

             q_graupel(i,k) = q_graupel(i,k) + rime_q
             q_cloud(i,k)   = q_cloud(i,k)   - rime_q
             n_cloud(i,k)   = n_cloud(i,k)   - rime_n

             ! ice multiplication based on Hallet and Mossop
             mult_q = 0.0_wp
             IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = const2*(T_a - T_mult_min) 
                mult_2 = const3*(T_a - T_mult_max) 
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,k)     = n_ice(i,k)     + mult_n
                q_ice(i,k)     = q_ice(i,k)     + mult_q
                q_graupel(i,k) = q_graupel(i,k) - mult_q
             ENDIF

             ! enhancement of melting of graupel
             IF (T_a > T_3 .AND. enhanced_melting) THEN
                melt_q = const4 * (T_a - T_3) * rime_q
                melt_n = melt_q / x_g

                melt_q = MIN(q_graupel(i,k),melt_q)
                melt_n = MIN(n_graupel(i,k),melt_n)
                
                q_graupel(i,k) = q_graupel(i,k) - melt_q
                q_rain(i,k)    = q_rain(i,k)    + melt_q
                n_graupel(i,k) = n_graupel(i,k) - melt_n
                n_rain(i,k)    = n_rain(i,k)    + melt_n
             ELSE
                melt_q = 0.0_wp
             ENDIF

             ! Shedding
             IF ((graupel_shedding .AND. D_g > D_shed_g .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
                q_g = q_graupel(i,k)
                n_g = n_graupel(i,k)
                x_g = graupel%meanmass(q_g,n_g)
                
                shed_q = MIN(q_g,rime_q)
                shed_n = shed_q / MIN(x_shed,x_g)
                
                q_graupel(i,k) = q_graupel(i,k) - shed_q
                q_rain(i,k)    = q_rain(i,k)    + shed_q                   
                n_rain(i,k)    = n_rain(i,k)    + shed_n
             ELSE
                shed_q = 0.0
             ENDIF
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE graupel_cloud_riming

  SUBROUTINE hail_cloud_riming()
    !*******************************************************************************
    ! see SB2006 for the equations                                                 *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h
    REAL(wp)            :: q_c,n_c,x_c,d_c,v_c,x_coll_c
    REAL(wp)            :: rime_n,rime_q,e_coll_n
    REAL(wp)            :: melt_n,melt_q,e_coll_q
    REAL(wp)            :: shed_n,shed_q,x_shed
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp), SAVE      :: delta_n_hh,delta_n_hc,delta_n_cc
    REAL(wp), SAVE      :: delta_q_hh,delta_q_hc,delta_q_cc
    REAL(wp), SAVE      :: theta_n_hh,theta_n_hc,theta_n_cc
    REAL(wp), SAVE      :: theta_q_hh,theta_q_hc,theta_q_cc
    REAL(wp)            :: const1,const2,const3,const4

    IF (isdebug) THEN
       WRITE(txt,*) " hail_cloud_riming" ; CALL message(routine,TRIM(txt)) 
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,cloud,0)
      delta_n_hc = coll_delta_12(hail,cloud,0)
      delta_n_cc = coll_delta_22(hail,cloud,0)
      delta_q_hh = coll_delta_11(hail,cloud,0) 
      delta_q_hc = coll_delta_12(hail,cloud,1)
      delta_q_cc = coll_delta_22(hail,cloud,1)

      theta_n_hh = coll_theta_11(hail,cloud,0)
      theta_n_hc = coll_theta_12(hail,cloud,0)
      theta_n_cc = coll_theta_22(hail,cloud,0)
      theta_q_hh = coll_theta_11(hail,cloud,0)
      theta_q_hc = coll_theta_12(hail,cloud,1)
      theta_q_cc = coll_theta_22(hail,cloud,1)

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "     delta_n_hh  = ",delta_n_hh ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     delta_n_hc  = ",delta_n_hc ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     delta_n_cc  = ",delta_n_cc ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     theta_n_hh  = ",theta_n_hh ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     theta_n_hc  = ",theta_n_hc ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     theta_n_cc  = ",theta_n_cc ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     delta_q_hh  = ",delta_q_hh ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     delta_q_hc  = ",delta_q_hc ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     delta_q_cc  = ",delta_q_cc ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     theta_q_hh  = ",theta_q_hh ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     theta_q_hc  = ",theta_q_hc ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "     theta_q_cc  = ",theta_q_cc ; CALL message(routine,TRIM(txt)) 
      END IF
      firstcall = 1
    ENDIF

    x_shed = 4./3.*pi * rho_w * r_shedding**3     !..mean mass of shedding droplets

    x_coll_c = (D_coll_c/cloud%a_geo)**3          !..lower mass for collection

    const1 = ecoll_hc/(D_coll_c - D_crit_c)    
    const2 = 1.0/(T_mult_opt - T_mult_min)
    const3 = 1.0/(T_mult_opt - T_mult_max)
    const4 = c_w / L_ew

    DO k = kstart,kend
       DO i = istart,iend

          T_a = T_p(i,k)
          q_c = q_cloud(i,k)                                 
          q_h = q_hail(i,k)                                  
          n_c = n_cloud(i,k)                                 
          n_h = n_hail(i,k)         

          x_h = hail%meanmass(q_h,n_h)
          D_h = hail%diameter(x_h)                  
          x_c = cloud%meanmass(q_c,n_c)
          D_c = cloud%diameter(x_c)
       
          IF (q_c > q_crit_c .AND. q_h > q_crit_hc .AND. D_h > D_crit_hc .AND. &
               & D_c > D_crit_c) THEN

            v_h = hail%velocity(x_h)  * rrho_04(i,k)
            v_c = cloud%velocity(x_c) * rrho_04(i,k)

            e_coll_n = MIN(ecoll_hc, MAX(const1*(D_c - D_crit_c),ecoll_min))
            e_coll_q = e_coll_n

            rime_n = pi4 * e_coll_n * n_h * n_c * dt & 
                 &   * (delta_n_hh * D_h*D_h + delta_n_hc * D_h*D_c + delta_n_cc * D_c*D_c) &
                 &   * SQRT(theta_n_hh * v_h*v_h - theta_n_hc * v_h*v_c + theta_n_cc * v_c*v_c)

            rime_q = pi4 * e_coll_q * n_h * q_c * dt & 
                 &   * (delta_q_hh * D_h*D_h + delta_q_hc * D_h*D_c + delta_q_cc * D_c*D_c) &
                 &   * SQRT(theta_q_hh * v_h*v_h - theta_q_hc * v_h*v_c + theta_q_cc * v_c*v_c)

            rime_q = MIN(q_c,rime_q)
            rime_n = MIN(n_c,rime_n)

            q_hail(i,k)  = q_hail(i,k) + rime_q
            q_cloud(i,k) = q_cloud(i,k) - rime_q
            n_cloud(i,k) = n_cloud(i,k) - rime_n

            ! ice multiplication
            mult_q = 0.0_wp
            IF (T_a < T_3 .AND. ice_multiplication) THEN
              mult_1 = const2*(T_a - T_mult_min) 
              mult_2 = const3*(T_a - T_mult_max) 
              mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
              mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
              mult_n = C_mult * mult_1 * mult_2 * rime_q
              mult_q = mult_n * ice%x_min
              mult_q = MIN(rime_q,mult_q)

              n_ice(i,k)  = n_ice(i,k)  + mult_n
              q_ice(i,k)  = q_ice(i,k)  + mult_q
              q_hail(i,k) = q_hail(i,k) - mult_q
            ENDIF

            ! enhancement of melting of hail
            IF (T_a > T_3 .AND. enhanced_melting) THEN
              melt_q = const4 * (T_a - T_3) * rime_q
              melt_n = melt_q / x_h

              melt_q = MIN(q_hail(i,k),melt_q)
              melt_n = MIN(n_hail(i,k),melt_n)

              q_hail(i,k) = q_hail(i,k) - melt_q
              q_rain(i,k) = q_rain(i,k) + melt_q
              n_hail(i,k) = n_hail(i,k) - melt_n
              n_rain(i,k) = n_rain(i,k) + melt_n
            ELSE
              melt_q = 0.0
            ENDIF

            ! Shedding
            IF ((D_h > D_shed_h .AND. T_a > T_shed .AND. hail_shedding) .OR. T_a > T_3 ) THEN
              q_h = q_hail(i,k)
              n_h = n_hail(i,k)
              x_h = hail%meanmass(q_h,n_h)

              shed_q = MIN(q_h,rime_q)
              shed_n = shed_q / MIN(x_shed,x_h)

              q_hail(i,k) = q_hail(i,k) - shed_q
              q_rain(i,k) = q_rain(i,k) + shed_q                   
              n_rain(i,k) = n_rain(i,k) + shed_n
            ELSE
              shed_q = 0.0
            ENDIF

          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE hail_cloud_riming

  SUBROUTINE graupel_rain_riming()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: q_r,n_r,x_r,d_r,v_r,x_shed
    REAL(wp)            :: rime_n,rime_q,melt_n,melt_q,shed_n,shed_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp), SAVE      :: delta_n_gg,delta_n_gr,delta_n_rr
    REAL(wp), SAVE      :: delta_q_gg,delta_q_gr,delta_q_rg,delta_q_rr
    REAL(wp), SAVE      :: theta_n_gg,theta_n_gr,theta_n_rr
    REAL(wp), SAVE      :: theta_q_gg,theta_q_gr,theta_q_rg,theta_q_rr
    REAL(wp)            :: const3,const4

    IF (isdebug) THEN
       WRITE(txt,*) " graupel_rain_riming" ; CALL message(routine,TRIM(txt)) 
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,rain,0)
      delta_n_gr = coll_delta_12(graupel,rain,0)
      delta_n_rr = coll_delta_22(graupel,rain,0)
      delta_q_gg = coll_delta_11(graupel,rain,0)
      delta_q_gr = coll_delta_12(graupel,rain,1)
      delta_q_rg = coll_delta_12(rain,graupel,1)
      delta_q_rr = coll_delta_22(graupel,rain,1)

      theta_n_gg = coll_theta_11(graupel,rain,0)
      theta_n_gr = coll_theta_12(graupel,rain,0)
      theta_n_rr = coll_theta_22(graupel,rain,0)
      theta_q_gg = coll_theta_11(graupel,rain,0)
      theta_q_gr = coll_theta_12(graupel,rain,1)
      theta_q_rg = coll_theta_12(rain,graupel,1)
      theta_q_rr = coll_theta_22(graupel,rain,1)

      IF (isdebug) THEN
        WRITE(txt,'(A,D10.3)') "     delta_n_gg = ",delta_n_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     delta_n_gr = ",delta_n_gr ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     delta_n_rr = ",delta_n_rr ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     theta_n_gg = ",theta_n_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     theta_n_gr = ",theta_n_gr ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     theta_n_rr = ",theta_n_rr ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     delta_q_gg = ",delta_q_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     delta_q_gr = ",delta_q_gr ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     delta_q_rg = ",delta_q_rg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     delta_q_rr = ",delta_q_rr ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     theta_q_gg = ",theta_q_gg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     theta_q_gr = ",theta_q_gr ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     theta_q_rg = ",theta_q_rg ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "     theta_q_rr = ",theta_q_rr ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ENDIF

    x_shed = 4./3. * pi * rho_w * r_shedding**3

    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)

    DO k = kstart,kend
       DO i = istart,iend

          q_r = q_rain(i,k)      
          q_g = q_graupel(i,k)   
          T_a = T_p(i,k)

          IF (q_r > q_crit .AND. q_g > q_crit) THEN

            n_r = n_rain(i,k)                                 
            n_g = n_graupel(i,k)                              

            x_g = graupel%meanmass(q_g,n_g)
            d_g = graupel%diameter(x_g)
            v_g = graupel%velocity(x_g) * rrho_04(i,k)
            x_r = rain%meanmass(q_r,n_r)
            d_r = rain%diameter(x_r)
            v_r = rain%velocity(x_r) * rrho_04(i,k)

            rime_n = pi4 * n_g * n_r * dt & 
                 & *     (delta_n_gg * D_g**2 + delta_n_gr * D_g*D_r + delta_n_rr * D_r**2) &
                 & * sqrt(theta_n_gg * v_g**2 - theta_n_gr * v_g*v_r + theta_n_rr * v_r**2)

            rime_q = pi4 * n_g * q_r * dt & 
                 & *     (delta_n_gg * D_g**2 + delta_q_gr * D_g*D_r + delta_q_rr * D_r**2) &
                 & * sqrt(theta_n_gg * v_g**2 - theta_q_gr * v_g*v_r + theta_q_rr * v_r**2)

            rime_n = MIN(n_r,rime_n)
            rime_q = MIN(q_r,rime_q)

            q_graupel(i,k) = q_graupel(i,k) + rime_q
            q_rain(i,k)    = q_rain(i,k)    - rime_q
            n_rain(i,k)    = n_rain(i,k)    - rime_n

            ! ice multiplication based on Hallet and Mossop
            mult_q = 0.0_wp
            IF (T_a < T_3 .AND. ice_multiplication) THEN
               mult_1 = (T_a - T_mult_min) * const3
               mult_2 = (T_a - T_mult_max) * const4
               mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
               mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
               mult_n = C_mult * mult_1 * mult_2 * rime_q
               mult_q = mult_n * ice%x_min
               mult_q = MIN(rime_q,mult_q)

               n_ice(i,k)     = n_ice(i,k)     + mult_n
               q_ice(i,k)     = q_ice(i,k)     + mult_q
               q_graupel(i,k) = q_graupel(i,k) - mult_q
            ENDIF

            ! enhancement of melting of graupel
            IF (T_a > T_3 .AND. enhanced_melting) THEN
               melt_q = c_w / L_ew * (T_a - T_3) * rime_q
               melt_n = melt_q / x_g
               
               melt_q = MIN(q_graupel(i,k),melt_q)
               melt_n = MIN(n_graupel(i,k),melt_n)
               
               q_graupel(i,k) = q_graupel(i,k) - melt_q
               q_rain(i,k)    = q_rain(i,k)    + melt_q

               n_graupel(i,k) = n_graupel(i,k) - melt_n
               n_rain(i,k)    = n_rain(i,k)    + melt_n
            ELSE
               melt_q = 0.0_wp
            ENDIF
            
            ! shedding
            IF ((graupel_shedding .AND. D_g > D_shed_g .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
               q_g = q_graupel(i,k)
               n_g = n_graupel(i,k)
               x_g = graupel%meanmass(q_g,n_g)
               
               shed_q = MIN(q_g,rime_q)
               IF (T_a <= T_3) THEN
                  shed_n = shed_q / MIN(x_shed,x_g)
               ELSE
                  shed_n = shed_q / MAX(x_r,x_g)
               ENDIF

               q_graupel(i,k) = q_graupel(i,k) - shed_q
               q_rain(i,k)    = q_rain(i,k)    + shed_q                   
               n_rain(i,k)    = n_rain(i,k)    + shed_n
            ELSE
               shed_q = 0.0
            ENDIF

         ENDIF
      ENDDO
   ENDDO

  END SUBROUTINE graupel_rain_riming

  SUBROUTINE hail_rain_riming()
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************

    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h
    REAL(wp)            :: q_r,n_r,x_r,d_r,v_r,x_shed
    REAL(wp)            :: rime_n,rime_q,melt_n,melt_q,shed_n,shed_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp), SAVE      :: delta_n_hh,delta_n_hr,delta_n_rr
    REAL(wp), SAVE      :: delta_q_hh,delta_q_hr,delta_q_rr
    REAL(wp), SAVE      :: theta_n_hh,theta_n_hr,theta_n_rr
    REAL(wp), SAVE      :: theta_q_hh,theta_q_hr,theta_q_rr
    REAL(wp)            :: const3,const4

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,rain,0)
      delta_n_hr = coll_delta_12(hail,rain,0)
      delta_n_rr = coll_delta_22(hail,rain,0)
      delta_q_hh = coll_delta_11(hail,rain,0) 
      delta_q_hr = coll_delta_12(hail,rain,1)
      delta_q_rr = coll_delta_22(hail,rain,1)

      theta_n_hh = coll_theta_11(hail,rain,0)
      theta_n_hr = coll_theta_12(hail,rain,0)
      theta_n_rr = coll_theta_22(hail,rain,0)
      theta_q_hh = coll_theta_11(hail,rain,0)
      theta_q_hr = coll_theta_12(hail,rain,1)
      theta_q_rr = coll_theta_22(hail,rain,1)

      IF (isdebug) THEN
         WRITE(txt,*) " hail_rain_riming:" ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     a_hail     = ",hail%a_geo ; CALL message(routine,TRIM(txt)) 
         WRITE(txt,'(A,D10.3)') "     b_hail     = ",hail%b_geo ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     alf_hail   = ",hail%a_vel ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     bet_hail   = ",hail%b_vel ; CALL message(routine,TRIM(txt)) 
         WRITE(txt,'(A,D10.3)') "     a_rain     = ",rain%a_geo ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     b_rain     = ",rain%b_geo ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     alf_rain   = ",rain%a_vel ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     bet_rain   = ",rain%b_vel ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     delta_n_hh = ",delta_n_hh ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     delta_n_hr = ",delta_n_hr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     delta_n_rr = ",delta_n_rr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     theta_n_hh = ",theta_n_hh ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     theta_n_hr = ",theta_n_hr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     theta_n_rr = ",theta_n_rr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     delta_q_hh = ",delta_q_hh ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     delta_q_hr = ",delta_q_hr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     delta_q_rr = ",delta_q_rr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     theta_q_hh = ",theta_q_hh ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     theta_q_hr = ",theta_q_hr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "     theta_q_rr = ",theta_q_rr ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
    ELSEIF (isdebug) THEN
       WRITE(txt,*) " hail_rain_riming" ; CALL message(routine,TRIM(txt))
    ENDIF

    x_shed = 4./3. * pi * rho_w * r_shedding**3

    const3 = 1.0/(T_mult_opt - T_mult_min)
    const4 = 1.0/(T_mult_opt - T_mult_max)

    DO k = kstart,kend
       DO i = istart,iend

          q_r = q_rain(i,k) 
          q_h = q_hail(i,k) 
          T_a = T_p(i,k) 

          IF (q_r > q_crit .AND. q_h > q_crit) THEN
            n_r = n_rain(i,k)                              
            n_h = n_hail(i,k)    

            x_h = hail%meanmass(q_h,n_h)
            D_h = hail%diameter(x_h)
            v_h = hail%velocity(x_h) * rrho_04(i,k)
                          
            x_r = rain%meanmass(q_r,n_r)
            D_r = rain%diameter(x_r)
            v_r = rain%velocity(x_r) * rrho_04(i,k)

            rime_n = pi4 * n_h * n_r * dt & 
                 & *     (delta_n_hh * D_h**2 + delta_n_hr * D_h*D_r + delta_n_rr * D_r**2) &
                 & * sqrt(theta_n_hh * v_h**2 - theta_n_hr * v_h*v_r + theta_n_rr * v_r**2)

            rime_q = pi4 * n_h * q_r * dt & 
                 & *     (delta_q_hh * D_h**2 + delta_q_hr * D_h*D_r + delta_q_rr * D_r**2) &
                 & * sqrt(theta_q_hh * v_h**2 - theta_q_hr * v_h*v_r + theta_q_rr * v_r**2)

            rime_n = MIN(n_r,rime_n)
            rime_q = MIN(q_r,rime_q)

            q_hail(i,k) = q_hail(i,k) + rime_q
            q_rain(i,k) = q_rain(i,k) - rime_q
            n_rain(i,k) = n_rain(i,k) - rime_n

            ! ice multiplication based on Hallet and Mossop
            mult_q = 0.0
            IF (T_a < T_3 .AND. ice_multiplication) THEN
              mult_1 = (T_a - T_mult_min) * const3
              mult_2 = (T_a - T_mult_max) * const4
              mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
              mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
              mult_n = C_mult * mult_1 * mult_2 * rime_q
              mult_q = mult_n * ice%x_min
              mult_q = MIN(rime_q,mult_q)

              n_ice(i,k)  = n_ice(i,k)  + mult_n
              q_ice(i,k)  = q_ice(i,k)  + mult_q
              q_hail(i,k) = q_hail(i,k) - mult_q
            ENDIF

            ! enhancement of melting of hail
            IF (T_a > T_3 .AND. enhanced_melting) THEN
              melt_q = c_w / L_ew * (T_a - T_3) * rime_q
              melt_n = melt_q / x_h

              melt_q = MIN(q_hail(i,k),melt_q)
              melt_n = MIN(n_hail(i,k),melt_n)

              q_hail(i,k) = q_hail(i,k) - melt_q
              q_rain(i,k) = q_rain(i,k) + melt_q

              n_hail(i,k) = n_hail(i,k) - melt_n
              n_rain(i,k) = n_rain(i,k) + melt_n
            ELSE
              melt_q = 0.0
            ENDIF

            ! shedding
            IF ((hail_shedding .AND. D_h > D_shed_h .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
              q_h = q_hail(i,k)
              n_h = n_hail(i,k)
              x_h = hail%meanmass(q_h,n_h)

              shed_q = MIN(q_h,rime_q)
              IF (T_a <= T_3) THEN
                shed_n = shed_q / MIN(x_shed,x_h)
              ELSE
                !  shed_n = shed_q / x_shed
                shed_n = shed_q / MAX(x_r,x_h)
              ENDIF

              q_hail(i,k) = q_hail(i,k) - shed_q
              q_rain(i,k) = q_rain(i,k)    + shed_q                   
              n_rain(i,k) = n_rain(i,k)    + shed_n
            ELSE
              shed_q = 0.0
            ENDIF

          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE hail_rain_riming

  SUBROUTINE graupel_melting()
    !*******************************************************************************
    ! Melting of graupel                                                           *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g,T_a,e_a
    REAL(wp)            :: melt,melt_v,melt_h,melt_n,melt_q
    REAL(wp)            :: fh_q,fv_q
    REAL(wp), SAVE      :: a_vent,b_vent

    IF (firstcall.NE.1) THEN
      a_vent = vent_coeff_a(graupel,1)
      b_vent = vent_coeff_b(graupel,1) * N_sc**n_f / sqrt(nu_l)
      firstcall = 1
      IF (isdebug) THEN
        WRITE(txt,*) " graupel_melting " 
        WRITE(txt,'(A,D10.3)') "   a_geo    = ",graupel%a_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   b_geo    = ",graupel%b_geo ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   a_vel    = ",graupel%a_vel ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "   b_vel    = ",graupel%b_vel ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "   a_ven    = ",graupel%a_ven ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "   b_ven    = ",graupel%b_ven ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "   a_vent = ",a_vent ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "   b_vent = ",b_vent ; CALL message(routine,TRIM(txt))
      ENDIF
    ELSEIF (isdebug) THEN
      WRITE(txt,*) " graupel_melting " ; CALL message(routine,TRIM(txt)) 
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          T_a = T_p(i,k)
          q_g = q_graupel(i,k)                  

          IF (T_a > T_3 .AND. q_g > 0.0) THEN
             e_a = e_ws(T_a)                     
             n_g = n_graupel(i,k)                
             
             x_g = graupel%meanmass(q_g,n_g)
             D_g = graupel%diameter(x_g)
             v_g = graupel%velocity(x_g) * rrho_04(i,k)

             fv_q = a_vent + b_vent * sqrt(v_g*D_g)
             fh_q = 1.05 * fv_q
             
             melt   = 2.0*pi / L_ew * D_g * n_g * dt

             melt_h = melt * K_T * (T_a - T_3)
             melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)
             
             melt_q = (melt_h * fh_q + melt_v * fv_q)
            
             ! UB: assume that x_g is constant during melting
             melt_n = MIN(MAX( (melt_q - q_g) / x_g + n_g, 0.0_wp), n_g)
             
             melt_q = MIN(q_g,melt_q)
             melt_n = MIN(n_g,melt_n)
             melt_q = MAX(0.0_wp,melt_q)
             melt_n = MAX(0.0_wp,melt_n)

             q_graupel(i,k) = q_graupel(i,k) - melt_q
             q_rain(i,k)    = q_rain(i,k)    + melt_q

             n_graupel(i,k) = n_graupel(i,k) - melt_n
             n_rain(i,k)    = n_rain(i,k)    + melt_n
             
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE graupel_melting

  SUBROUTINE hail_melting()
    !*******************************************************************************
    ! Melting of hail                                                              *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h,T_a,e_a
    REAL(wp)            :: melt,melt_v,melt_h,melt_n,melt_q
    REAL(wp)            :: fh_q,fv_q
    REAL(wp), SAVE      :: a_vent,b_vent

    IF (firstcall.NE.1) THEN
      a_vent = vent_coeff_a(hail,1)
      b_vent = vent_coeff_b(hail,1) * N_sc**n_f / sqrt(nu_l)
      firstcall = 1
      IF (isdebug) THEN
        WRITE(txt,*) " hail_melting: " ; CALL message(routine,TRIM(txt)) 
        WRITE(txt,'(A,D10.3)') "    a_vent = ",a_vent ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,D10.3)') "    b_vent = ",b_vent ; CALL message(routine,TRIM(txt))
      ENDIF
    ELSEIF (isdebug) THEN
      WRITE(txt,*) " hail_melting " ; CALL message(routine,TRIM(txt)) 
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          T_a = T_p(i,k)
          q_h = q_hail(i,k)    

          IF (T_a > T_3 .AND. q_h > 0.0_wp) THEN
            e_a = e_ws(T_a)    
            n_h = n_hail(i,k)  

            x_h = hail%meanmass(q_h,n_h)
            D_h = hail%diameter(x_h)
            v_h = hail%velocity(x_h) * rrho_04(i,k)

            fv_q = a_vent + b_vent * sqrt(v_h*D_h)
            fh_q = 1.05 * fv_q                            ! UB: based on Rasmussen and Heymsfield

            melt   = 2.0*pi / L_ew * D_h * n_h * dt

            melt_h = melt * K_T * (T_a - T_3)
            melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)

            melt_q = (melt_h * fh_q + melt_v * fv_q)

            ! UB: assume that x_h is constant during melting
            melt_n = MIN(MAX( (melt_q - q_h) / x_h + n_h, 0.0_wp), n_h)

            melt_q = MIN(q_h,melt_q)
            melt_n = MIN(n_h,melt_n)
            melt_q = MAX(0.0_wp,melt_q)
            melt_n = MAX(0.0_wp,melt_n)

            q_hail(i,k) = q_hail(i,k) - melt_q
            q_rain(i,k) = q_rain(i,k) + melt_q

            n_hail(i,k) = n_hail(i,k) - melt_n
            n_rain(i,k) = n_rain(i,k) + melt_n

         ENDIF
      ENDDO
   ENDDO

  END SUBROUTINE hail_melting

  SUBROUTINE graupel_hail_conv_wet_gamlook()
    !*******************************************************************************
    !  Wet growth and conversion of graupel to hail                                *
    !  (uses look-up table for incomplete gamma functions)                         *
    !  by Uli Blahak                                                               *
    !*******************************************************************************
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a, p_a, d_trenn, qw_a, qi_a, N_0, lam, xmin
    REAL(wp)            :: q_g,n_g,x_g,d_g
    REAL(wp)            :: q_c,q_r
    REAL(wp)            :: conv_n,conv_q

    TYPE(gamlookuptable), SAVE :: ltable1, ltable2
    REAL(wp), SAVE             :: nm1, nm2, g1, g2

    IF (firstcall.NE.1) THEN
      firstcall = 1
      nm1 = (graupel%nu+1.0)/graupel%mu
      nm2 = (graupel%nu+2.0)/graupel%mu
      CALL incgfct_lower_lookupcreate(nm1, ltable1, nlookup, nlookuphr_dummy)
      CALL incgfct_lower_lookupcreate(nm2, ltable2, nlookup, nlookuphr_dummy)     
      g1 = ltable1%igf(ltable1%n) ! ordinary gamma function of nm1 is the last value in lookup table 1      
      g2 = ltable2%igf(ltable2%n) ! ordinary gamma function of nm2 is the last value in lookup table 2
      IF (isdebug) THEN
        WRITE(txt,*) "graupel_hail_conv_wet_gamlook" ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,ES12.5)') "(nug+1.0)/mug :    ", nm1  ; CALL message(routine,TRIM(txt))
        WRITE(txt,'(A,ES12.5)') "(nug+2.0)/mug :    ", nm2  ; CALL message(routine,TRIM(txt))
      ENDIF
    ELSE IF (isdebug) THEN
      WRITE(txt,*) "graupel_hail_conv_wet_gamlook" ; CALL message(routine,TRIM(txt)) 
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          q_c = q_cloud(i,k)                       
          q_r = q_rain(i,k)                        
          q_g = q_graupel(i,k)                     
          n_g = n_graupel(i,k)            

          x_g = graupel%meanmass(q_g,n_g)
          D_g = graupel%diameter(x_g)
          n_g = q_g / x_g  ! for consistency for limiters, n_g is used explicitly below

          T_a = T_p(i,k)
          p_a = p_p(i,k)

          !..supercooled liquid water in the cloud environment = sum of rain and cloud water
          qw_a = q_r + q_c

          IF (T_a < T_3 .AND. q_g > q_crit_gc .AND. qw_a > 1e-3) THEN

            !.. Umgebungsgehalt Eispartikel (vernachl. werden Graupel und Hagel wg. geringer Kollisionseff.)
            !.. koennte problematisch sein, weil in konvekt. Wolken viel mehr Graupel und Hagel enthalten ist!!!
            qi_a = q_ice(i,k) + q_snow(i,k) 
            d_trenn = dmin_wg_gr_ltab_equi(p_a,T_a,qw_a,qi_a,ltabdminwgg)

            IF (d_trenn > 0.0_wp .AND. d_trenn < 10.0_wp * D_g) THEN 

               xmin = exp(log(d_trenn/graupel%a_geo)*(1.0_wp/graupel%b_geo))
               
               lam  = exp(log(g2/(g1*x_g))*(graupel%mu))
               xmin = exp(log(xmin)*graupel%mu)
               n_0  = graupel%mu * n_g * exp(log(lam)*nm1) / g1

               conv_n = n_0 / (graupel%mu * exp(log(lam)*nm1)) * incgfct_upper_lookup(lam*xmin, ltable1)
               conv_q = n_0 / (graupel%mu * exp(log(lam)*nm2)) * incgfct_upper_lookup(lam*xmin, ltable2)
               
               conv_n = MIN(conv_n,n_g)
               conv_q = MIN(conv_q,q_g)
               
               q_graupel(i,k) = q_g - conv_q
               n_graupel(i,k) = n_g - conv_n
               
               q_hail(i,k) = q_hail(i,k) + conv_q
               n_hail(i,k) = n_hail(i,k) + conv_n

            END IF
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE graupel_hail_conv_wet_gamlook

  SUBROUTINE ice_riming (nsize, nlev, dep_rate_ice)
    !*******************************************************************************
    !  Riming of ice with cloud droplet and rain drops. First the process rates    *
    !  are calculated in                                                           *
    !      snow_cloud_riming ()                                                    *
    !      snow_rain_riming ()                                                     *
    !  using those rates and the previously calculated and stored deposition       *
    !  rate the conversion of snow to graupel and rain is done.                    *
    !*******************************************************************************

    INTEGER,  INTENT (IN) :: nsize, nlev
    REAL(wp), INTENT (IN), DIMENSION(nsize,nlev) :: dep_rate_ice

    REAL(wp), DIMENSION(nsize,nlev) :: rime_rate_qc, rime_rate_nc, &
         &               rime_rate_qi, rime_rate_qr, rime_rate_nr

    INTEGER             :: i,k
    REAL(wp)            :: q_i,n_i,x_i,d_i
    REAL(wp)            :: T_a,x_coll_c,x_r
    REAL(wp)            :: rime_n,rime_q,rime_qr,rime_qi
    REAL(wp)            :: conv_n,conv_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp)            :: const2,const3,const4,const5

    IF (isdebug) THEN
       WRITE(txt,*) " ice riming" ; CALL message(routine,TRIM(txt))
    END IF

    rime_rate_qc(:,:) = 0.0_wp
    rime_rate_qr(:,:) = 0.0_wp
    rime_rate_qi(:,:) = 0.0_wp
    rime_rate_nc(:,:) = 0.0_wp
    rime_rate_qr(:,:) = 0.0_wp

    CALL ice_cloud_riming()
    CALL ice_rain_riming()

    x_coll_c = (D_coll_c/cloud%a_geo)**3   !..lower mean mass for collection

    const2 = 1.0/x_coll_c
    const3 = 1.0/(T_mult_opt - T_mult_min)
    const4 = 1.0/(T_mult_opt - T_mult_max)
    const5 = alpha_spacefilling * rho_w/rho_ice
    !
    ! Complete ice-cloud and ice-rain riming

    DO k = kstart,kend
       DO i = istart,iend

          T_a = T_p(i,k)

          IF (dep_rate_ice(i,k) > 0.0_wp &
               & .and. dep_rate_ice(i,k) .ge. rime_rate_qc(i,k)+rime_rate_qr(i,k)) THEN

            !
            ! 1) Depositional growth is stronger than riming growth, therefore ice stays ice
            !

            !.. ice_cloud_riming

            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(q_cloud(i,k),rime_q)
              rime_n = MIN(n_cloud(i,k),rime_n)

              q_ice(i,k)   = q_ice(i,k)  + rime_q
              q_cloud(i,k) = q_cloud(i,k) - rime_q
              n_cloud(i,k) = n_cloud(i,k) - rime_n

              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4 
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                n_ice(i,k) = n_ice(i,k)  + mult_n
              ENDIF
            END IF

            !.. ice_rain_riming
            IF (rime_rate_qr(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qr(i,k)
              rime_n = rime_rate_nr(i,k)
              rime_q = MIN(q_rain(i,k),rime_q)
              rime_n = MIN(n_rain(i,k),rime_n)
              
              q_ice(i,k)  = q_ice(i,k)  + rime_q
              q_rain(i,k) = q_rain(i,k) - rime_q
              n_rain(i,k) = n_rain(i,k) - rime_n

              !..ice multiplication               
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                
                n_ice(i,k) = n_ice(i,k)  + mult_n
              ENDIF
              
            END IF

         ELSE

            !
            !.. 2) Depositional growth negative or smaller than riming growth, therefore ice is 
            !      allowed to convert to graupel and / or hail
            !

            !.. ice_cloud_riming
            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              n_i = n_ice(i,k)
              q_i = q_ice(i,k)
              x_i = ice%meanmass(q_i,n_i)
              D_i = ice%diameter(x_i)

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(q_cloud(i,k),rime_q)
              rime_n = MIN(n_cloud(i,k),rime_n)

              q_ice(i,k)   = q_ice(i,k)   + rime_q
              q_cloud(i,k) = q_cloud(i,k) - rime_q
              n_cloud(i,k) = n_cloud(i,k) - rime_n

              ! ice multiplication
              mult_q = 0.0_wp
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4 
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                n_ice(i,k) = n_ice(i,k)  + mult_n              
              ENDIF

              ! conversion ice -> graupel (depends on alpha_spacefilling)
              IF (D_i > D_conv_ig) THEN
                 q_i = q_ice(i,k)
                 conv_q = (rime_q - mult_q) / ( const5 * (pi6*rho_ice*d_i**3/x_i - 1.0) )
                 conv_q = MIN(q_i,conv_q)
                 x_i    = ice%meanmass(q_i,n_i)
                 conv_n = conv_q / MAX(x_i,x_conv) 
                 conv_n = MIN(n_ice(i,k),conv_n)
              ELSE
                 conv_q = 0.0_wp
                 conv_n = 0.0_wp
              ENDIF

              q_ice(i,k)     = q_ice(i,k)     - conv_q
              q_graupel(i,k) = q_graupel(i,k) + conv_q
              n_ice(i,k)     = n_ice(i,k)     - conv_n
              n_graupel(i,k) = n_graupel(i,k) + conv_n
            END IF

            !.. ice_rain_riming
            IF (rime_rate_qi(i,k) > 0.0_wp) THEN

              rime_qi = rime_rate_qi(i,k)
              rime_qr = rime_rate_qr(i,k)
              rime_n  = rime_rate_nr(i,k)
              rime_n  = MIN(MIN(n_rain(i,k),n_ice(i,k)),rime_n)
              rime_qr = MIN(q_rain(i,k),rime_qr)
              rime_qi = MIN(q_ice(i,k),rime_qi)

              n_ice(i,k)  = n_ice(i,k)  - rime_n
              n_rain(i,k) = n_rain(i,k) - rime_n
              q_ice(i,k)  = q_ice(i,k)  - rime_qi
              q_rain(i,k) = q_rain(i,k) - rime_qr

              ! ice multiplication
              mult_q = 0.0
              mult_n = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)
              ENDIF

              IF (T_a >= T_3) THEN
                 ! shedding of rain at warm temperatures
                 ! i.e. undo time integration, but with modified n_rain
                 x_r = rain%meanmass(q_rain(i,k),n_rain(i,k))
                 n_ice(i,k)  = n_ice(i,k)  + rime_n
                 n_rain(i,k) = n_rain(i,k) + rime_qr / x_r
                 q_ice(i,k)  = q_ice(i,k)  + rime_qi     
                 q_rain(i,k) = q_rain(i,k) + rime_qr
              ELSE
                 ! new ice particles from multiplication
                 n_ice(i,k) = n_ice(i,k) + mult_n 
                 q_ice(i,k) = q_ice(i,k) + mult_q 
                 ! riming to graupel
                 n_graupel(i,k) = n_graupel(i,k) + rime_n
                 q_graupel(i,k) = q_graupel(i,k) + rime_qi + rime_qr - mult_q
              END IF
            END IF
          END IF
       END DO
    END DO

  CONTAINS

    SUBROUTINE ice_cloud_riming()
      !*******************************************************************************
      !  Riming rate of ice collecting cloud droplets                                *
      !*******************************************************************************
      INTEGER             :: i,k
      INTEGER, SAVE       :: firstcall
      REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
      REAL(wp)            :: q_c,n_c,x_c,d_c,v_c,e_coll,x_coll_c
      REAL(wp)            :: rime_n,rime_q
      REAL(wp), SAVE      :: delta_n_ii,delta_n_ic,delta_n_cc
      REAL(wp), SAVE      :: delta_q_ii,delta_q_ic,delta_q_cc
      REAL(wp), SAVE      :: theta_n_ii,theta_n_ic,theta_n_cc
      REAL(wp), SAVE      :: theta_q_ii,theta_q_ic,theta_q_cc
      REAL(wp)            :: const1

      IF (isdebug) THEN
         WRITE(txt,*) "ice_cloud_riming" ; CALL message(routine,TRIM(txt)) 
      END IF

      IF (firstcall.NE.1) THEN
         delta_n_ii = coll_delta_11(ice,cloud,0)
         delta_n_ic = coll_delta_12(ice,cloud,0)
         delta_n_cc = coll_delta_22(ice,cloud,0)
         delta_q_ii = coll_delta_11(ice,cloud,0) 
         delta_q_ic = coll_delta_12(ice,cloud,1)
         delta_q_cc = coll_delta_22(ice,cloud,1)

         theta_n_ii = coll_theta_11(ice,cloud,0)
         theta_n_ic = coll_theta_12(ice,cloud,0)
         theta_n_cc = coll_theta_22(ice,cloud,0)
         theta_q_ii = coll_theta_11(ice,cloud,0)
         theta_q_ic = coll_theta_12(ice,cloud,1)
         theta_q_cc = coll_theta_22(ice,cloud,1)
         
         IF (isdebug) THEN
            WRITE(txt,'(A,D10.3)') "    a_ice      = ",ice%a_geo ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "    b_ice      = ",ice%b_geo ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    alf_ice    = ",ice%a_vel ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "    bet_ice    = ",ice%b_vel ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "    a_cloud    = ",cloud%a_geo ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    b_cloud    = ",cloud%b_geo ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    alf_cloud  = ",cloud%a_vel ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "    bet_cloud  = ",cloud%b_vel ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "    delta_n_ii = ",delta_n_ii ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    delta_n_ic = ",delta_n_ic ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    delta_n_cc = ",delta_n_cc ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    theta_n_ii = ",theta_n_ii ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    theta_n_ic = ",theta_n_ic ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    theta_n_cc = ",theta_n_cc ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    delta_q_ii = ",delta_q_ii ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    delta_q_ic = ",delta_q_ic ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    delta_q_cc = ",delta_q_cc ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    theta_q_ii = ",theta_q_ii ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    theta_q_ic = ",theta_q_ic ; CALL message(routine,TRIM(txt))
            WRITE(txt,'(A,D10.3)') "    theta_q_cc = ",theta_q_cc ; CALL message(routine,TRIM(txt))
         END IF
         firstcall = 1
      ENDIF
      
      const1   = ecoll_ic/(D_coll_c - D_crit_c)    !..collision efficiency coeff
      x_coll_c = (D_coll_c/cloud%a_geo)**3         !..lower mass threshold for collection
    
      DO k = kstart,kend
         DO i = istart,iend
            
            n_c = n_cloud(i,k)      
            q_c = q_cloud(i,k)      
            n_i = n_ice(i,k)        
            q_i = q_ice(i,k)        
            
            x_c = cloud%meanmass(q_c,n_c)
            D_c = cloud%diameter(x_c)
            x_i = ice%meanmass(q_i,n_i)
            D_i = ice%diameter(x_i)

            IF (q_c > q_crit_c .AND. q_i > q_crit_ic .AND. D_i > D_crit_ic .AND. D_c > D_crit_c) THEN

               v_c = cloud%velocity(x_c) * rrho_04(i,k)
               v_i = ice%velocity(x_i)   * rrho_04(i,k)
               
               e_coll = MIN(ecoll_ic, MAX(const1*(D_c - D_crit_c), ecoll_min))

               rime_n = pi4 * e_coll * n_i * n_c * dt & 
                    &   *     (delta_n_ii * D_i*D_i + delta_n_ic * D_i*D_c + delta_n_cc * D_c*D_c) &
                    &   * SQRT(theta_n_ii * v_i*v_i - theta_n_ic * v_i*v_c + theta_n_cc * v_c*v_c  &
                    &         +ice_s_vel**2)

               rime_q = pi4 * e_coll * n_i * q_c * dt & 
                    &   *     (delta_q_ii * D_i*D_i + delta_q_ic * D_i*D_c + delta_q_cc * D_c*D_c) &
                    &   * SQRT(theta_q_ii * v_i*v_i - theta_q_ic * v_i*v_c + theta_q_cc * v_c*v_c  &
                    &          +ice_s_vel**2)

               rime_rate_qc(i,k) = rime_q
               rime_rate_nc(i,k) = rime_n
               
            ENDIF
         ENDDO
      ENDDO
    END SUBROUTINE ice_cloud_riming

    SUBROUTINE ice_rain_riming()
      !*******************************************************************************
      !  Riming rate of ice collecting rain drop, or rain collecting ice             *
      !*******************************************************************************      
      INTEGER             :: i,k
      INTEGER, SAVE       :: firstcall
      REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
      REAL(wp)            :: q_r,n_r,x_r,d_r,v_r
      REAL(wp)            :: rime_n,rime_qi,rime_qr
      REAL(wp), SAVE      :: delta_n_ii,delta_n_ir,           delta_n_rr
      REAL(wp), SAVE      :: delta_q_ii,delta_q_ir,delta_q_ri,delta_q_rr
      REAL(wp), SAVE      :: theta_n_ii,theta_n_ir,           theta_n_rr
      REAL(wp), SAVE      :: theta_q_ii,theta_q_ir,theta_q_ri,theta_q_rr
      
      IF (isdebug) THEN
         WRITE(txt,*) "ice_rain_riming " ; CALL message(routine,TRIM(txt)) 
      END IF

      IF (firstcall.NE.1) THEN
         delta_n_ii = coll_delta_11(ice,rain,0)
         delta_n_ir = coll_delta_12(ice,rain,0)
         delta_n_rr = coll_delta_22(ice,rain,0)
         delta_q_ii = coll_delta_11(ice,rain,0) 
         delta_q_ir = coll_delta_12(ice,rain,1)
         delta_q_ri = coll_delta_12(rain,ice,1)
         delta_q_rr = coll_delta_22(ice,rain,1)

         theta_n_ii = coll_theta_11(ice,rain,0)
         theta_n_ir = coll_theta_12(ice,rain,0)
         theta_n_rr = coll_theta_22(ice,rain,0)
         theta_q_ii = coll_theta_11(ice,rain,0)
         theta_q_ir = coll_theta_12(ice,rain,1)
         theta_q_ri = coll_theta_12(rain,ice,1)
         theta_q_rr = coll_theta_22(ice,rain,1)
         
         IF (isdebug) THEN
            WRITE(txt,'(A,D10.3)') "     a_ice      = ",ice%a_geo ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     b_ice      = ",ice%b_geo ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     alf_ice    = ",ice%a_vel ; CALL message(routine,TRIM(txt))  
            WRITE(txt,'(A,D10.3)') "     bet_ice    = ",ice%b_vel ; CALL message(routine,TRIM(txt))  
            WRITE(txt,'(A,D10.3)') "     a_rain    = ",rain%a_geo ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     b_rain    = ",rain%b_geo ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     alf_rain  = ",rain%a_vel ; CALL message(routine,TRIM(txt))  
            WRITE(txt,'(A,D10.3)') "     bet_rain  = ",rain%b_vel ; CALL message(routine,TRIM(txt))  
            WRITE(txt,'(A,D10.3)') "     delta_n_ii = ",delta_n_ii ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     delta_n_ir = ",delta_n_ir ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     delta_n_rr = ",delta_n_rr ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     theta_n_ii = ",theta_n_ii ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     theta_n_ir = ",theta_n_ir ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     theta_n_rr = ",theta_n_rr ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     delta_q_ii = ",delta_q_ii ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     delta_q_ir = ",delta_q_ir ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     delta_q_ri = ",delta_q_ri ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     delta_q_rr = ",delta_q_rr ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     theta_q_ii = ",theta_q_ii ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     theta_q_ir = ",theta_q_ir ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     theta_q_ri = ",theta_q_ri ; CALL message(routine,TRIM(txt)) 
            WRITE(txt,'(A,D10.3)') "     theta_q_rr = ",theta_q_rr ; CALL message(routine,TRIM(txt)) 
         END IF
         firstcall = 1
      ENDIF

      DO k = kstart,kend
         DO i = istart,iend
            
            q_r = q_rain(i,k)              
            q_i = q_ice(i,k)               
            n_r = n_rain(i,k)              
            n_i = n_ice(i,k)               
            
            x_i = ice%meanmass(q_i,n_i)
            D_i = ice%diameter(x_i)
            
            IF (q_r > q_crit .AND. q_i > q_crit_ir .AND. D_i > D_crit_ir) THEN

               x_r = rain%meanmass(q_r,n_r)
               D_r = rain%diameter(x_r)
               v_r = rain%velocity(x_r) * rrho_04(i,k)
               v_i = ice%velocity(x_i) * rrho_04(i,k)

               rime_n  = pi4 * n_i * n_r * dt & 
                    &  *     (delta_n_ii * D_i*D_i + delta_n_ir * D_i*D_r + delta_n_rr * D_r*D_r) &
                    &  * sqrt(theta_n_ii * v_i*v_i - theta_n_ir * v_i*v_r + theta_n_rr * v_r*v_r  &
                    &         +ice_s_vel**2)

               rime_qr = pi4 * n_i * q_r * dt & 
                    &  *     (delta_n_ii * D_i*D_i + delta_q_ir * D_i*D_r + delta_q_rr * D_r*D_r) &
                    &  * sqrt(theta_n_ii * v_i*v_i - theta_q_ir * v_i*v_r + theta_q_rr * v_r*v_r  &
                    &        +ice_s_vel**2)
               
               rime_qi = pi4 * n_r * q_i * dt & 
                    &  *     (delta_q_ii * D_i*D_i + delta_q_ri * D_i*D_r + delta_n_rr * D_r*D_r) &
                    &  * sqrt(theta_q_ii * v_i*v_i - theta_q_ri * v_i*v_r + theta_n_rr * v_r*v_r  &
                    &        +ice_s_vel**2)

               rime_rate_nr(i,k) = rime_n
               rime_rate_qi(i,k) = rime_qi
               rime_rate_qr(i,k) = rime_qr
               
            ENDIF
         ENDDO
      ENDDO
   
    END SUBROUTINE ice_rain_riming

  END SUBROUTINE ice_riming

  SUBROUTINE snow_riming (nsize, nlev, dep_rate_snow)
    !*******************************************************************************
    !  Riming of snow with cloud droplet and rain drops. First the process rates   *
    !  are calculated in                                                           *
    !      snow_cloud_riming ()                                                    *
    !      snow_rain_riming ()                                                     *
    !  using those rates and the previously calculated and stored deposition       *
    !  rate the conversion of snow to graupel and rain is done.                    *
    !*******************************************************************************

    INTEGER, INTENT (IN) :: nsize, nlev
    REAL(wp), INTENT(IN), DIMENSION(nsize,nlev) :: dep_rate_snow

    REAL(wp), DIMENSION(nsize,nlev) ::              &
         & rime_rate_qc, rime_rate_nc,    &
         & rime_rate_qs, rime_rate_qr, rime_rate_nr

    INTEGER             :: i,k
    REAL(wp)            :: T_a  
    REAL(wp)            :: q_s,n_s,x_s,d_s
    REAL(wp)            :: x_coll_c,x_r
    REAL(wp)            :: rime_n,rime_q,rime_qr,rime_qs
    REAL(wp)            :: conv_n,conv_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp)            :: const2,const3,const4,const5

    IF (isdebug) THEN
       WRITE(txt,*) "snow_riming" ; CALL message(routine,TRIM(txt))
    END IF

    rime_rate_qc(:,:) = 0.0_wp
    rime_rate_qr(:,:) = 0.0_wp
    rime_rate_qs(:,:) = 0.0_wp
    rime_rate_nc(:,:) = 0.0_wp
    rime_rate_qr(:,:) = 0.0_wp

    CALL snow_cloud_riming()
    CALL snow_rain_riming()

    x_coll_c = (D_coll_c/cloud%a_geo)**3   !..lower mean mass for collection

    const2 = 1.0/x_coll_c
    const3 = 1.0/(T_mult_opt - T_mult_min)
    const4 = 1.0/(T_mult_opt - T_mult_max)
    const5 = alpha_spacefilling * rho_w/rho_ice

    DO k = kstart,kend
       DO i = istart,iend

          T_a = T_p(i,k)

          IF (dep_rate_snow(i,k) > 0.0_wp &
               & .AND. dep_rate_snow(i,k) .ge. rime_rate_qc(i,k) + rime_rate_qr(i,k)) THEN
             
            !
            ! 1) Depositional growth is stronger than riming growth, therefore snow stays snow:
            !

            !.. time integration of snow_cloud_riming

            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(q_cloud(i,k),rime_q)
              rime_n = MIN(n_cloud(i,k),rime_n)

              q_snow(i,k)  = q_snow(i,k)  + rime_q
              q_cloud(i,k) = q_cloud(i,k) - rime_q
              n_cloud(i,k) = n_cloud(i,k) - rime_n

              ! ice multiplication
              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,k)  = n_ice(i,k)  + mult_n
                q_ice(i,k)  = q_ice(i,k)  + mult_q
                q_snow(i,k) = q_snow(i,k) - mult_q
              ENDIF

            END IF

            !.. time integration snow_rain_riming

            IF (rime_rate_qr(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qr(i,k)
              rime_n = rime_rate_nr(i,k)
              rime_q = MIN(q_rain(i,k),rime_q)
              rime_n = MIN(n_rain(i,k),rime_n)

              q_snow(i,k) = q_snow(i,k) + rime_q
              q_rain(i,k) = q_rain(i,k) - rime_q
              n_rain(i,k) = n_rain(i,k) - rime_n

              ! ice multiplication
              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,k)  = n_ice(i,k)  + mult_n
                q_ice(i,k)  = q_ice(i,k)  + mult_q
                q_snow(i,k) = q_snow(i,k) - mult_q
              ENDIF
            END IF

          ELSE

            !
            !.. 2) Depositional growth is negative or smaller than riming growth, therefore snow is 
            !      allowed to convert to graupel and / or hail:
            !

            !.. time integration of snow_cloud_riming

            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              n_s = n_snow(i,k)
              q_s = q_snow(i,k)
              x_s = snow%meanmass(q_s,n_s)
              D_s = snow%diameter(x_s)

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(q_cloud(i,k),rime_q)
              rime_n = MIN(n_cloud(i,k),rime_n)

              q_snow(i,k)  = q_snow(i,k)  + rime_q
              q_cloud(i,k) = q_cloud(i,k) - rime_q
              n_cloud(i,k) = n_cloud(i,k) - rime_n

              ! ice multiplication
              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,k)  = n_ice(i,k)  + mult_n
                q_ice(i,k)  = q_ice(i,k)  + mult_q
                q_snow(i,k) = q_snow(i,k) - mult_q
              ENDIF

              !.. conversion of snow to graupel, depends on alpha_spacefilling

              IF (D_s > D_conv_sg) THEN
                 conv_q = (rime_q - mult_q) / ( const5 * (pi6*rho_ice*d_s**3/x_s - 1.0) )
                 conv_q = MIN(q_s,conv_q)
                 x_s    = snow%meanmass(q_s,n_s)
                 conv_n = conv_q / MAX(x_s,x_conv) 
                 conv_n = MIN(n_snow(i,k),conv_n)
              ELSE
                 conv_q = 0.0_wp
                 conv_n = 0.0_wp
              ENDIF

              q_snow(i,k)    = q_snow(i,k)    - conv_q
              q_graupel(i,k) = q_graupel(i,k) + conv_q

              n_snow(i,k)    = n_snow(i,k)    - conv_n
              n_graupel(i,k) = n_graupel(i,k) + conv_n

            END IF

            !.. time integration of snow_rain_riming

            IF (rime_rate_qs(i,k) > 0.0_wp) THEN

              rime_qs = rime_rate_qs(i,k)
              rime_qr = rime_rate_qr(i,k)
              rime_n  = rime_rate_nr(i,k)
              rime_qr = MIN(q_rain(i,k),rime_qr)
              rime_qs = MIN(q_snow(i,k),rime_qs)
              rime_n  = MIN(n_rain(i,k),rime_n)
              rime_n  = MIN(n_snow(i,k),rime_n)

              n_snow(i,k) = n_snow(i,k) - rime_n
              n_rain(i,k) = n_rain(i,k) - rime_n
              q_snow(i,k) = q_snow(i,k) - rime_qs
              q_rain(i,k) = q_rain(i,k) - rime_qr

              ! ice multiplication
              mult_q = 0.0_wp
              mult_n = 0.0_wp
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)
              ENDIF

              IF (T_a >= T_3) THEN
                 ! shedding of rain at warm temperatures
                 ! i.e. undo time integration, but with modified n_rain
                 x_r = rain%meanmass(q_rain(i,k),n_rain(i,k))
                 n_snow(i,k) = n_snow(i,k) + rime_n
                 n_rain(i,k) = n_rain(i,k) + rime_qr / x_r
                 q_snow(i,k) = q_snow(i,k) + rime_qs
                 q_rain(i,k) = q_rain(i,k) + rime_qr
              ELSE
                 ! new ice particles from multiplication
                 n_ice(i,k)  = n_ice(i,k)  + mult_n
                 q_ice(i,k)  = q_ice(i,k)  + mult_q
                 ! riming to graupel
                 n_graupel(i,k) = n_graupel(i,k) + rime_n
                 q_graupel(i,k) = q_graupel(i,k) + rime_qr + rime_qs - mult_q
              END IF              

           END IF           
        END IF

      END DO
    END DO

 CONTAINS

   SUBROUTINE snow_rain_riming()
     !*******************************************************************************
     ! Riming of snow with raindrops                                                *
     !*******************************************************************************     
     INTEGER             :: i,k
     INTEGER, SAVE       :: firstcall
     REAL(wp)            :: q_s,n_s,x_s,d_s,v_s
     REAL(wp)            :: q_r,n_r,x_r,d_r,v_r
     REAL(wp)            :: rime_n,rime_qr,rime_qs
     REAL(wp), SAVE      :: delta_n_ss,delta_n_sr,delta_n_rr
     REAL(wp), SAVE      :: delta_q_ss,delta_q_sr,delta_q_rs,delta_q_rr
     REAL(wp), SAVE      :: theta_n_ss,theta_n_sr,theta_n_rr
     REAL(wp), SAVE      :: theta_q_ss,theta_q_sr,theta_q_rs,theta_q_rr
     
     IF (isdebug) THEN
        WRITE(txt,*) "snow_rain_riming" ; CALL message(routine,TRIM(txt))
     END IF

     IF (firstcall.NE.1) THEN
        delta_n_ss = coll_delta_11(snow,rain,0)
        delta_n_sr = coll_delta_12(snow,rain,0)
        delta_n_rr = coll_delta_22(snow,rain,0)
        delta_q_ss = coll_delta_11(snow,rain,0) 
        delta_q_sr = coll_delta_12(snow,rain,1)
        delta_q_rs = coll_delta_12(rain,snow,1)
        delta_q_rr = coll_delta_22(snow,rain,1)
      
        theta_n_ss = coll_theta_11(snow,rain,0)
        theta_n_sr = coll_theta_12(snow,rain,0)
        theta_n_rr = coll_theta_22(snow,rain,0)
        theta_q_ss = coll_theta_11(snow,rain,0)
        theta_q_sr = coll_theta_12(snow,rain,1)
        theta_q_rs = coll_theta_12(rain,snow,1)
        theta_q_rr = coll_theta_22(snow,rain,1)
        
        IF (isdebug) THEN
         WRITE(txt,'(A,D10.3)') "    a_snow     = ",snow%a_geo ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    b_snow     = ",snow%b_geo ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    alf_snow   = ",snow%a_vel ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    bet_snow   = ",snow%b_vel ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    a_rain     = ",rain%a_geo ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    b_rain     = ",rain%b_geo ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    alf_rain   = ",rain%a_vel ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    bet_rain   = ",rain%b_vel ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    delta_n_ss = ",delta_n_ss ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    delta_n_sr = ",delta_n_sr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    delta_n_rr = ",delta_n_rr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    theta_n_ss = ",theta_n_ss ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    theta_n_sr = ",theta_n_sr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    theta_n_rr = ",theta_n_rr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    delta_q_ss = ",delta_q_ss ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    delta_q_sr = ",delta_q_sr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    delta_q_rs = ",delta_q_rs ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    delta_q_rr = ",delta_q_rr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    theta_q_ss = ",theta_q_ss ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    theta_q_sr = ",theta_q_sr ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    theta_q_rs = ",theta_q_rs ; CALL message(routine,TRIM(txt))
         WRITE(txt,'(A,D10.3)') "    theta_q_rr = ",theta_q_rr ; CALL message(routine,TRIM(txt))
      END IF
      firstcall = 1
   ENDIF

   DO k = kstart,kend
      DO i = istart,iend

         q_r = q_rain(i,k)            
         q_s = q_snow(i,k)            
         n_r = n_rain(i,k)            
         n_s = n_snow(i,k)            
         
         x_s = snow%meanmass(q_s,n_s)
         D_s = snow%diameter(x_s)
                  
         IF (q_r > q_crit .AND. q_s > q_crit_sr .AND. D_s > D_crit_sr) THEN
                        
            x_r = rain%meanmass(q_r,n_r)
            D_r = rain%diameter(x_r)
            v_r = rain%velocity(x_r) * rrho_04(i,k)
            v_s = snow%velocity(x_s) * rrho_04(i,k)
            
            rime_n = pi4 * n_s * n_r * dt & 
                 & *     (delta_n_ss * D_s**2 + delta_n_sr * D_s*D_r + delta_n_rr * D_r**2) &
                 & * sqrt(theta_n_ss * v_s**2 - theta_n_sr * v_s*v_r + theta_n_rr * v_r**2  &
                 &       +snow_s_vel**2)

            rime_qr = pi4 * n_s * q_r * dt & 
                 &  *     (delta_n_ss * D_s**2 + delta_q_sr * D_s*D_r + delta_q_rr * D_r**2) &
                 &  * sqrt(theta_n_ss * v_s**2 - theta_q_sr * v_s*v_r + theta_q_rr * v_r**2  &
                 &        +snow_s_vel**2)

            rime_qs = pi4 * n_r * q_s * dt & 
                 &  *     (delta_q_ss * D_s**2 + delta_q_rs * D_s*D_r + delta_n_rr * D_r**2) &
                 &  * sqrt(theta_q_ss * v_s**2 - theta_q_rs * v_s*v_r + theta_n_rr * v_r**2  &
                 &        +snow_s_vel**2)

            rime_rate_nr(i,k)  = rime_n
            rime_rate_qs(i,k)  = rime_qs
            rime_rate_qr(i,k)  = rime_qr
            
         ENDIF
       ENDDO
     ENDDO    
   END SUBROUTINE snow_rain_riming

   SUBROUTINE snow_cloud_riming()
     !*******************************************************************************
     ! Riming of snow with cloud droplets                                           *
     !*******************************************************************************     
     INTEGER             :: i,k
     INTEGER, SAVE       :: firstcall
     REAL(wp)            :: q_s,n_s,x_s,d_s,v_s
     REAL(wp)            :: q_c,n_c,x_c,d_c,v_c,e_coll
     REAL(wp)            :: rime_n,rime_q
     REAL(wp)            :: const1
     REAL(wp), SAVE      :: delta_n_ss,delta_n_sc,delta_n_cc
     REAL(wp), SAVE      :: delta_q_ss,delta_q_sc,delta_q_cc
     REAL(wp), SAVE      :: theta_n_ss,theta_n_sc,theta_n_cc
     REAL(wp), SAVE      :: theta_q_ss,theta_q_sc,theta_q_cc
     
     IF (isdebug) THEN
        WRITE(txt,*) "snow_cloud_riming " ; CALL message(routine,TRIM(txt))
     END IF
     
     IF (firstcall.NE.1) THEN
        delta_n_ss = coll_delta_11(snow,cloud,0)
        delta_n_sc = coll_delta_12(snow,cloud,0)
        delta_n_cc = coll_delta_22(snow,cloud,0)
        delta_q_ss = coll_delta_11(snow,cloud,0) 
        delta_q_sc = coll_delta_12(snow,cloud,1)
        delta_q_cc = coll_delta_22(snow,cloud,1)
        
        theta_n_ss = coll_theta_11(snow,cloud,0)
        theta_n_sc = coll_theta_12(snow,cloud,0)
        theta_n_cc = coll_theta_22(snow,cloud,0)
        theta_q_ss = coll_theta_11(snow,cloud,0)
        theta_q_sc = coll_theta_12(snow,cloud,1)
        theta_q_cc = coll_theta_22(snow,cloud,1)
        
        IF (isdebug) THEN
           WRITE(txt,'(A,D10.3)') "   a_snow      = ",snow%a_geo ; CALL message(routine,TRIM(txt)) 
           WRITE(txt,'(A,D10.3)') "   b_snow      = ",snow%b_geo ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   alf_snow    = ",snow%a_vel ; CALL message(routine,TRIM(txt)) 
           WRITE(txt,'(A,D10.3)') "   bet_snow    = ",snow%b_vel ; CALL message(routine,TRIM(txt)) 
           WRITE(txt,'(A,D10.3)') "   a_cloud    = ",cloud%a_geo ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   b_cloud    = ",cloud%b_geo ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   alf_cloud  = ",cloud%a_vel ; CALL message(routine,TRIM(txt)) 
           WRITE(txt,'(A,D10.3)') "   bet_cloud  = ",cloud%b_vel ; CALL message(routine,TRIM(txt)) 
           WRITE(txt,'(A,D10.3)') "   delta_n_ss = ",delta_n_ss ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   delta_n_sc = ",delta_n_sc ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   delta_n_cc = ",delta_n_cc ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   theta_n_ss = ",theta_n_ss ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   theta_n_sc = ",theta_n_sc ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   theta_n_cc = ",theta_n_cc ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   delta_q_ss = ",delta_q_ss ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   delta_q_sc = ",delta_q_sc ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   delta_q_cc = ",delta_q_cc ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   theta_q_ss = ",theta_q_ss ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   theta_q_sc = ",theta_q_sc ; CALL message(routine,TRIM(txt))
           WRITE(txt,'(A,D10.3)') "   theta_q_cc = ",theta_q_cc ; CALL message(routine,TRIM(txt))
        END IF
        firstcall = 1
     ENDIF
     
     const1 = ecoll_ic/(D_coll_c - D_crit_c)

     DO k = kstart,kend
        DO i = istart,iend

           n_c = n_cloud(i,k)  
           q_c = q_cloud(i,k)  
           n_s = n_snow(i,k)   
           q_s = q_snow(i,k)   
           
           x_c = cloud%meanmass(q_c,n_c)
           D_c = cloud%diameter(x_c)           
           x_s = snow%meanmass(q_s,n_s)
           D_s = snow%diameter(x_s)

           IF (q_c > q_crit_c .AND. q_s > q_crit_sc .AND. D_s > D_crit_sc .AND. D_c > D_crit_c) THEN
              
              v_c = cloud%velocity(x_c) * rrho_04(i,k)
              v_s = snow%velocity(x_s)  * rrho_04(i,k)
              
              e_coll = MIN(ecoll_sc, MAX(const1*(D_c - D_crit_c), ecoll_min))
              
              rime_n = pi4 * e_coll * n_s * n_c * dt & 
                   & *     (delta_n_ss * D_s**2 + delta_n_sc * D_s*D_c + delta_n_cc * D_c**2) &
                   & * sqrt(theta_n_ss * v_s**2 - theta_n_sc * v_s*v_c + theta_n_cc * v_c**2  &
                   &       +snow_s_vel**2)
              
              rime_q = pi4 * e_coll * n_s * q_c * dt & 
                   & *     (delta_q_ss * D_s**2 + delta_q_sc * D_s*D_c + delta_q_cc * D_c**2) &
                   & * sqrt(theta_q_ss * v_s**2 - theta_q_sc * v_s*v_c + theta_q_cc * v_c**2  &
                   &       +snow_s_vel**2)
              
              rime_rate_qc(i,k) = rime_q
              rime_rate_nc(i,k) = rime_n

           ENDIF
        ENDDO
     ENDDO

   END SUBROUTINE snow_cloud_riming

 END SUBROUTINE snow_riming

  !*******************************************************************************
  ! Sedimentation subroutines for ICON
  !*******************************************************************************

  SUBROUTINE sedi_icon_rain (qp,np,precrate,qc,rhocorr,adz,dt, &
      &                  its,ite,kts,kte,cmax) !

    INTEGER, INTENT(IN) :: its,ite,kts,kte
    REAL(wp), DIMENSION(its:ite,kts:kte), INTENT(INOUT) :: qp,np,qc
    REAL(wp), DIMENSION(its:ite,kts:kte), INTENT(IN)    :: adz,rhocorr
    REAL(wp), DIMENSION(its:ite),         INTENT(OUT)   :: precrate
    REAL(wp) :: dt
    REAL(wp), INTENT(INOUT), OPTIONAL :: cmax

    INTEGER  :: i, k, kk
    REAL(wp) :: x_p,D_m,D_p,mue,v_n,v_q

    REAL(wp), DIMENSION(its:ite,kts-1:kte+1) :: q_fluss, n_fluss
    REAL(wp), DIMENSION(its:ite,kts-1:kte+1) :: v_n_sedi,v_q_sedi
    REAL(wp), DIMENSION(its:ite) :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    LOGICAL,  DIMENSION(its:ite) :: cflag

    v_n_sedi = 0.0_wp
    v_q_sedi = 0.0_wp
    q_fluss  = 0.0_wp
    n_fluss  = 0.0_wp

    DO k = kts,kte
       DO i = its,ite

          IF (qp(i,k) > q_crit) THEN

             x_p = rain%meanmass(qp(i,k) ,np(i,k))
             D_m = rain%diameter(x_p)
             IF (qc(i,k) >= q_crit) THEN 
                mue = (rain%nu+1.0_wp)/rain%b_geo - 1.0_wp
             ELSE
                mue = rain_coeffs%mue_Dm_relation(D_m)
             END IF
             D_p = D_m * exp((-1./3.)*log((mue+3.)*(mue+2.)*(mue+1.)))          

             v_n = rain_coeffs%alfa - rain_coeffs%beta * exp(-(mue+1.)*log(1.0 + rain_coeffs%gama*D_p))
             v_q = rain_coeffs%alfa - rain_coeffs%beta * exp(-(mue+4.)*log(1.0 + rain_coeffs%gama*D_p))
             v_n = v_n * rhocorr(i,k)
             v_q = v_q * rhocorr(i,k)
             
             v_n_sedi(i,k) = - v_n
             v_q_sedi(i,k) = - v_q
          ENDIF
       END DO
    END DO

    v_n_sedi(:,kte-1) = v_n_sedi(:,kte)   ! untere Randbedingung fuer Fallgeschw.
    v_q_sedi(:,kte-1) = v_q_sedi(:,kte)   ! lower BC for the terminal veloc.

    DO k = kts+1,kte

        DO i = its,ite
          v_nv(i) = 0.5 * (v_n_sedi(i,k-1)+v_n_sedi(i,k))   ! Das sollte wohl k-1 sein.
          v_qv(i) = 0.5 * (v_q_sedi(i,k-1)+v_q_sedi(i,k))
          ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
          c_nv(i) = -v_nv(i) * adz(i,k) * dt 
          c_qv(i) = -v_qv(i) * adz(i,k) * dt
        END DO
        IF (PRESENT(cmax)) cmax = MAX(cmax,MAXVAL(c_qv))

        kk = k
        s_nv = 0.0
        DO i = its,ite
          IF (c_nv(i) <= 1) THEN
            s_nv(i) = v_nv(i) * np(i,k)
          END IF
        END DO
        IF (ANY(c_nv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_nv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_nv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_nv(i) = s_nv(i) + np(i,kk)/adz(i,kk)
                c_nv(i) = (c_nv(i) - 1) * adz(i,kk-1)/adz(i,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_nv(i) = s_nv(i) + np(i,kk)/adz(i,kk)*MIN(c_nv(i),1.0_wp)
              s_nv(i) = -s_nv(i) / dt
            END IF
          END DO
        ENDIF

        kk = k
        s_qv = 0.0
        DO i = its,ite
          IF (c_qv(i) <= 1) THEN
            s_qv(i) = v_qv(i) * qp(i,k)
          END IF
        END DO
        IF (ANY(c_qv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_qv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_qv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_qv(i) = s_qv(i) + qp(i,kk)/adz(i,kk)
                c_qv(i) = (c_qv(i) - 1) * adz(i,kk-1)/adz(i,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_qv(i) = s_qv(i) + qp(i,kk)/adz(i,kk)*MIN(c_qv(i),1.0_wp)
              s_qv(i) = -s_qv(i) / dt
            END IF
          END DO
        ENDIF

        ! Flux-limiter to avoid negative values
        DO i = its,ite
          n_fluss(i,k) = MAX(s_nv(i),n_fluss(i,k-1)-np(i,k)/(adz(i,k)*dt))
          q_fluss(i,k) = MAX(s_qv(i),q_fluss(i,k-1)-qp(i,k)/(adz(i,k)*dt))
        END DO

    END DO

    n_fluss(:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,kts-1) = 0.0 ! upper BC

    DO k = kts,kte
        DO i = its,ite
          np(i,k) = np(i,k) + ( n_fluss(i,k) - n_fluss(i,k-1) )*adz(i,k)*dt
          qp(i,k) = qp(i,k) + ( q_fluss(i,k) - q_fluss(i,k-1) )*adz(i,k)*dt
        ENDDO
    ENDDO
    precrate(its:ite) = - q_fluss(its:ite,kte) ! Regenrate

    RETURN

  END SUBROUTINE sedi_icon_rain


  SUBROUTINE sedi_icon_sphere (ptype,pcoeffs,qp,np,precrate,rhocorr,adz,dt, &
      &                  its,ite,kts,kte,cmax) !

    TYPE(particle), INTENT(in) :: ptype
    TYPE(particle_sphere), INTENT(in) :: pcoeffs
    INTEGER, INTENT(IN) :: its,ite,kts,kte
    REAL(wp), DIMENSION(its:ite,kts:kte), INTENT(INOUT) :: qp,np
    REAL(wp), DIMENSION(its:ite,kts:kte), INTENT(IN)    :: adz,rhocorr
    REAL(wp), DIMENSION(its:ite),         INTENT(OUT)   :: precrate
    REAL(wp) :: dt
    REAL(wp), INTENT(INOUT), OPTIONAL :: cmax

    INTEGER  :: i, k, kk
    REAL(wp) :: x_p,v_n,v_q,lam

    REAL(wp), DIMENSION(its:ite,kts-1:kte+1) :: q_fluss, n_fluss
    REAL(wp), DIMENSION(its:ite,kts-1:kte+1) :: v_n_sedi,v_q_sedi
    REAL(wp), DIMENSION(its:ite) :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    LOGICAL,  DIMENSION(its:ite) :: cflag

    v_n_sedi = 0.0_wp
    v_q_sedi = 0.0_wp
    q_fluss  = 0.0_wp
    n_fluss  = 0.0_wp

    DO k = kts,kte
       DO i = its,ite
          IF (qp(i,k) > q_crit) THEN
             
             x_p = ptype%meanmass(qp(i,k),np(i,k))
             lam = exp(ptype%b_vel* log(pcoeffs%coeff_lambda*x_p))
             
             v_n = pcoeffs%coeff_alfa_n * lam
             v_q = pcoeffs%coeff_alfa_q * lam
             v_n = MAX(v_n,ptype%vsedi_min)
             v_q = MAX(v_q,ptype%vsedi_min)
             v_n = MIN(v_n,ptype%vsedi_max)
             v_q = MIN(v_q,ptype%vsedi_max)             
             v_n = v_n * rhocorr(i,k)
             v_q = v_q * rhocorr(i,k)
             
             v_n_sedi(i,k) = -v_n
             v_q_sedi(i,k) = -v_q
          END IF
       END DO
    END DO
    
    v_n_sedi(:,kte-1) = v_n_sedi(:,kte)   ! untere Randbedingung fuer Fallgeschw.
    v_q_sedi(:,kte-1) = v_q_sedi(:,kte)   ! lower BC for the terminal veloc.

    DO k = kts+1,kte

        DO i = its,ite
          v_nv(i) = 0.5 * (v_n_sedi(i,k-1)+v_n_sedi(i,k))   ! Das sollte wohl k-1 sein.
          v_qv(i) = 0.5 * (v_q_sedi(i,k-1)+v_q_sedi(i,k))
          ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
          c_nv(i) = -v_nv(i) * adz(i,k) * dt 
          c_qv(i) = -v_qv(i) * adz(i,k) * dt
        END DO
        IF (PRESENT(cmax)) cmax = MAX(cmax,MAXVAL(c_qv))

        kk = k
        s_nv = 0.0
        DO i = its,ite
          IF (c_nv(i) <= 1) THEN
            s_nv(i) = v_nv(i) * np(i,k)
          END IF
        END DO
        IF (ANY(c_nv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_nv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_nv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_nv(i) = s_nv(i) + np(i,kk)/adz(i,kk)
                c_nv(i) = (c_nv(i) - 1) * adz(i,kk-1)/adz(i,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_nv(i) = s_nv(i) + np(i,kk)/adz(i,kk)*MIN(c_nv(i),1.0_wp)
              s_nv(i) = -s_nv(i) / dt
            END IF
          END DO
        ENDIF

        kk = k
        s_qv = 0.0
        DO i = its,ite
          IF (c_qv(i) <= 1) THEN
            s_qv(i) = v_qv(i) * qp(i,k)
          END IF
        END DO
        IF (ANY(c_qv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_qv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_qv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_qv(i) = s_qv(i) + qp(i,kk)/adz(i,kk)
                c_qv(i) = (c_qv(i) - 1) * adz(i,kk-1)/adz(i,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_qv(i) = s_qv(i) + qp(i,kk)/adz(i,kk)*MIN(c_qv(i),1.0_wp)
              s_qv(i) = -s_qv(i) / dt
            END IF
          END DO
        ENDIF

        ! Flux-limiter to avoid negative values
        DO i = its,ite
          n_fluss(i,k) = MAX(s_nv(i),n_fluss(i,k-1)-np(i,k)/(adz(i,k)*dt))
          q_fluss(i,k) = MAX(s_qv(i),q_fluss(i,k-1)-qp(i,k)/(adz(i,k)*dt))
        END DO

    END DO

    n_fluss(:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,kts-1) = 0.0 ! upper BC

    DO k = kts,kte
        DO i = its,ite
          np(i,k) = np(i,k) + ( n_fluss(i,k) - n_fluss(i,k-1) )*adz(i,k)*dt
          qp(i,k) = qp(i,k) + ( q_fluss(i,k) - q_fluss(i,k-1) )*adz(i,k)*dt
        ENDDO
    ENDDO
    precrate(its:ite) = - q_fluss(its:ite,kte) ! Regenrate

    RETURN

  END SUBROUTINE sedi_icon_sphere

#endif
  
END MODULE mo_2mom_mcrph_main
