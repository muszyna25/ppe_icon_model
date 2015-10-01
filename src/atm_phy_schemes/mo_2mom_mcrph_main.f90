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
! Version of May 2015 by AS:
! - New IN and CCN routines implemented based on Hande et al. (HDCP2-M3)
! - gscp=4 has now prognostic QNC and IN depletion (n_inact)
! - gscp=5 has additional budget equations for IN and CCN
!===============================================================================!
! To Do:
! - Check conservation of water mass
! - Further optimization might be possible in rain_freeze.
!===============================================================================!
! Further plans (physics) including HDCP2 project:
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
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
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

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_main'


  ! switches for ice scheme, ice nucleation, drop activation and autoconversion
  INTEGER  :: ice_typ, nuc_i_typ, nuc_c_typ, auto_typ

  ! for constant droplet number runs
  real(wp) :: qnc_const

  ! Derived type for atmospheric variables
  TYPE ATMOSPHERE
     REAL(wp), pointer, dimension(:,:) :: w, p, t, rho, qv
  END TYPE ATMOSPHERE

  ! Derived type for hydrometeor species including pointers to data
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
    REAL(wp), pointer, dimension(:,:) :: n     !..number density
    REAL(wp), pointer, dimension(:,:) :: q     !..mass density
    REAL(wp), pointer, dimension(:,:) :: rho_v !..density correction of terminal fall velocity
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
  TYPE, EXTENDS(particle_nonsphere) :: particle_rain
    REAL(wp)      :: cmu0   !..Parameters for mu-D-relation of rain
    REAL(wp)      :: cmu1   !     max of left branch
    REAL(wp)      :: cmu2   !     max of right branch
    REAL(wp)      :: cmu3   !     min value of relation
    REAL(wp)      :: cmu4   !     location of min value = breakup equilibrium diameter
    INTEGER       :: cmu5   !     exponent
  CONTAINS
    PROCEDURE :: mue_Dm_relation => rain_mue_Dm_relation
  END TYPE particle_rain

  TYPE aerosol_ccn
     REAL(wp)      :: Ncn0      ! CN concentration at ground
     REAL(wp)      :: Nmin      ! minimum value for CCN
     REAL(wp)      :: lsigs     ! log(sigma_s)
     REAL(wp)      :: R2        ! in mum
     REAL(wp)      :: etas      ! soluble fraction
     REAL(wp)      :: z0        ! parameter for height-dependency, constant up to z0_nccn
     REAL(wp)      :: z1e       ! 1/e scale height of N_ccn profile
  END TYPE aerosol_ccn

  TYPE aerosol_in
     REAL(wp)      :: N0     ! CN concentration at ground
     REAL(wp)      :: z0     ! parameter for height-dependency, constant up to z0_nccn
     REAL(wp)      :: z1e    ! 1/e scale height of N_ccn profile
  END TYPE aerosol_in

  ! These are the fundamental particle variables for the two-moment scheme
  ! These pointers will be specified from the list of pre-defined particles
  TYPE(particle_sphere)          :: ice_coeffs, snow_coeffs, graupel_coeffs, hail_coeffs
  TYPE(particle_rain) :: rain_coeffs
  TYPE(aerosol_ccn)              :: ccn_coeffs
  TYPE(aerosol_in)               :: in_coeffs

  ! Physical parameters and coefficients which occur only in the two-moment scheme

  ! ... some physical parameters not found in ICON
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

  ! .. Phillips et al. ice nucleation scheme, see ice_nucleation_homhet() for more details
  REAL(wp) ::                         &
       &    na_dust    = 162.e3_wp,   & ! initial number density of dust [1/m³], Phillips08 value 162e3
       &    na_soot    =  15.e6_wp,   & ! initial number density of soot [1/m³], Phillips08 value 15e6
       &    na_orga    = 177.e6_wp,   & ! initial number density of organics [1/m3], Phillips08 value 177e6
       &    ni_het_max = 50.0e3_wp,   & ! max number of IN between 1-10 per liter, i.e. 1d3-10d3
       &    ni_hom_max = 5000.0e3_wp    ! number of liquid aerosols between 100-5000 per liter

  INTEGER, PARAMETER ::               & ! Look-up table for Phillips et al. nucleation
       &    ttmax  = 30,              & ! sets limit for temperature in look-up table
       &    ssmax  = 60,              & ! sets limit for ice supersaturation in look-up table
       &    ttstep = 2,               & ! increment for temperature in look-up table
       &    ssstep = 1                  ! increment for ice supersaturation in look-up table

  REAL(wp), DIMENSION(0:100,0:100)  :: &
       &    afrac_dust, &  ! look-up table of activated fraction of dust particles acting as ice nuclei
       &    afrac_soot, &  ! ... of soot particles
       &    afrac_orga     ! ... of organic material

  INCLUDE 'phillips_nucleation_2010.incf'

  ! Look up tables for graupel_hail_conv_wet_gamlook and rain_freeze_gamlook

  TYPE(gamlookuptable) :: graupel_ltable1, graupel_ltable2
  TYPE(gamlookuptable) :: rain_ltable1, rain_ltable2, rain_ltable3
  REAL(wp)             :: rain_nm1, rain_nm2, rain_nm3, rain_g1, rain_g2
  REAL(wp)             :: graupel_nm1, graupel_nm2, graupel_g1, graupel_g2

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

  ! debug switches
  LOGICAL, PARAMETER     :: isdebug = .false.   ! use only when really desperate
  LOGICAL, PARAMETER     :: ischeck = .true.   ! frequently check for positive definite q's

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


    TYPE(particle_rain), PARAMETER :: rainSBBcoeffs = particle_rain( &
         &     9.292000,  & !..alfa
         &     9.623000,  & !..beta
         &     6.222d+2,  & !..gama
         &     6.000000,  & !..cmu0
         &     3.000d+1,  & !..cmu1
         &     1.000d+3,  & !..cmu2
         &     1.100d-3,  & !..cmu3 = D_br
         &     1.000000,  & !..cmu4
         &     2 )          !..cmu5

  REAL(wp), PARAMETER :: pi6 = pi/6.0_wp, pi8 = pi/8.0_wp ! more pieces of pi

  PUBLIC :: atmosphere, particle

  PUBLIC :: init_2mom_scheme, init_2mom_scheme_once, clouds_twomoment

  PUBLIC :: sedi_vel_rain, sedi_vel_sphere, sedi_icon_rain,sedi_icon_sphere

  PUBLIC :: q_crit
  ! these should be fixed to private for reasons of encapsulation
  PUBLIC :: rain_coeffs, ice_coeffs, snow_coeffs, graupel_coeffs, hail_coeffs, &
       ccn_coeffs, in_coeffs

  PUBLIC :: qnc_const

  TYPE sym_riming_params
     REAL(wp) :: delta_n_aa, delta_n_ab, delta_n_bb, &
          &      delta_q_aa, delta_q_ab, delta_q_bb, &
          &      theta_n_aa, theta_n_ab, theta_n_bb, &
          &      theta_q_aa, theta_q_ab, theta_q_bb
  END TYPE sym_riming_params

  TYPE asym_riming_params
    REAL(wp) :: delta_n_aa,delta_n_ab,delta_n_bb, &
         &      delta_q_aa,delta_q_ab,delta_q_ba,delta_q_bb, &
         &      theta_n_aa,theta_n_ab,theta_n_bb, &
         &      theta_q_aa,theta_q_ab,theta_q_ba,theta_q_bb
  END TYPE asym_riming_params

  !> run-time- and location-invariant snow cloud riming parameters
  TYPE(sym_riming_params), SAVE :: scr_params
  !> run-time- and location-invariant snow rain riming parameters
  TYPE(asym_riming_params), SAVE :: srr_params
  !> run-time- and location-invariant ice rain riming parameters
  TYPE(asym_riming_params), SAVE :: irr_params
  !> run-time- and location-invariant ice cloud riming parameters
  TYPE(sym_riming_params), SAVE :: icr_params
  !> run-time- and location-invariant hail rain riming parameters
  TYPE(sym_riming_params), SAVE :: hrr_params
  !> run-time- and location-invariant graupel rain riming parameters
  TYPE(asym_riming_params), SAVE :: grr_params
  !> run-time- and location-invariant hail cloud riming parameters
  TYPE(sym_riming_params), SAVE :: hcr_params
  !> run-time- and location-invariant graupel cloud riming parameters
  TYPE(sym_riming_params), SAVE :: gcr_params
  !> run-time- and location-invariant snow ice collection parameters
  TYPE(sym_riming_params), SAVE :: sic_params
  !> run-time- and location-invariant hail ice collection parameters
  TYPE(sym_riming_params), SAVE :: hic_params
  !> run-time- and location-invariant graupel ice collection parameters
  TYPE(sym_riming_params), SAVE :: gic_params
  !> run-time- and location-invariant hail snow collection parameters
  TYPE(sym_riming_params), SAVE :: hsc_params
  !> run-time- and location-invariant graupel snow collection parameters
  TYPE(sym_riming_params), SAVE :: gsc_params

  TYPE evaporation_deposition_params
    REAL(wp) :: c          !< coeff for capacity
    REAL(wp) :: a_f,b_f    !< coeffs for ventilation
  END TYPE evaporation_deposition_params

  !> run-time- and location-invariant vapor ice deposition parameters
  TYPE(evaporation_deposition_params), SAVE :: vid_params
  !> run-time- and location-invariant vapor graupel deposition parameters
  TYPE(evaporation_deposition_params), SAVE :: vgd_params
  !> run-time- and location-invariant vapor hail deposition parameters
  TYPE(evaporation_deposition_params), SAVE :: vhd_params
  !> run-time- and location-invariant vapor snow deposition parameters
  TYPE(evaporation_deposition_params), SAVE :: vsd_params

  !> run-time- and location-invariant graupel evaporation parameters
  TYPE(evaporation_deposition_params), SAVE :: ge_params
  !> run-time- and location-invariant hail evaporation parameters
  TYPE(evaporation_deposition_params), SAVE :: he_params
  !> run-time- and location-invariant snow evaporation parameters
  TYPE(evaporation_deposition_params), SAVE :: se_params

  TYPE melt_params
    REAL(wp) :: a_vent, b_vent
  END TYPE melt_params

  !> run-time- and location-invariant graupel melting parameters
  TYPE(melt_params), SAVE :: gm_params
  !> run-time- and location-invariant hail melting parameters
  TYPE(melt_params), SAVE :: hm_params
  !> run-time- and location-invariant snow melting parameters
  TYPE(melt_params), SAVE :: sm_params

  !> run-time- and location invariant for graupel selfcollection
  REAL(wp), SAVE :: graupel_sc_coll_n

  !> run-time- and location invariant for snow selfcollection
  REAL(wp), SAVE :: snow_sc_delta_n, snow_sc_theta_n

CONTAINS

  !*******************************************************************************
  ! This subroutine has to be called once at the start of the model run by
  ! the main program. It properly sets the parameters for the different hydrometeor
  ! classes according to predefined parameter sets (see above).
  !*******************************************************************************

  SUBROUTINE init_2mom_scheme(cloud,rain,ice,snow,graupel,hail)
    TYPE(particle), INTENT(out) :: cloud, rain, ice, snow, graupel, hail

    !..Pre-defined particle types
    TYPE(particle), PARAMETER :: graupelhail_cosmo5 = particle( & ! graupelhail2test4
         &        'graupelhail_cosmo5' ,& !.name...Bezeichnung
         &        1.000000, & !..nu.....Breiteparameter der Verteil.
         &        0.333333, & !..mu.....Exp.-parameter der Verteil.
         &        5.00d-04, & !..x_max..maximale Teilchenmasse
         &        1.00d-09, & !..x_min..minimale Teilchenmasse
         &        1.42d-01, & !..a_geo..Koeff. Geometrie
         &        0.314000, & !..b_geo..Koeff. Geometrie = 1/3.10
         &        86.89371, & !..a_vel..Koeff. Fallgesetz
         &        0.268325, & !..b_vel..Koeff. Fallgesetz
         &        0.780000, & !..a_ven..Koeff. Ventilation (PK, S.541)
         &        0.308000, & !..b_ven..Koeff. Ventilation (PK, S.541)
         &        2.00,     & !..cap....Koeff. Kapazitaet
         &        30.0,     & !..vsedi_max
         &        0.10,     & !..vsedi_min
         &        null(),   & !..n pointer
         &        null(),   & !..q pointer
         &        null() )    !..rho_v pointer
    TYPE(particle), PARAMETER :: hail_cosmo5 = particle( & ! hailULItest
         &        'hail_cosmo5' ,& !.name...Bezeichnung
         &        1.000000, & !..nu.....Breiteparameter der Verteil.
         &        0.333333, & !..mu.....Exp.-parameter der Verteil.
         &        5.00d-04, & !..x_max..maximale Teilchenmasse
         &        2.60d-9,  & !..x_min..minimale Teilchenmasse
         &        0.1366 ,  & !..a_geo..Koeff. Geometrie
         &        0.333333, & !..b_geo..Koeff. Geometrie = 1/3
         &        39.3    , & !..a_vel..Koeff. Fallgesetz
         &        0.166667, & !..b_vel..Koeff. Fallgesetz
         &        0.780000, & !..a_ven..Koeff. Ventilation (PK, S.541)
         &        0.308000, & !..b_ven..Koeff. Ventilation (PK, S.541)
         &        2.00,     & !..cap....Koeff. Kapazitaet
         &        30.0,     & !..vsedi_max
         &        0.1,      & !..vsedi_min
         &        null(),   & !..n pointer
         &        null(),   & !..q pointer
         &        null() )    !..rho_v pointer

    TYPE(particle), PARAMETER :: cloud_cosmo5 = PARTICLE( &
         &      'cloud_cosmo5',  & !.name...Bezeichnung der Partikelklasse
         &        0.0,      & !..nu.....Breiteparameter der Verteil.
         &        0.333333, & !..mu.....Exp.-parameter der Verteil.
         &        2.60d-10, & !..x_max..maximale Teilchenmasse D=80e-6m
         &        4.20d-15, & !..x_min..minimale Teilchenmasse D=2.e-6m
         &        1.24d-01, & !..a_geo..Koeff. Geometrie
         &        0.333333, & !..b_geo..Koeff. Geometrie = 1/3
         &        3.75d+05, & !..a_vel..Koeff. Fallgesetz
         &        0.666667, & !..b_vel..Koeff. Fallgesetz
         &        0.780000, & !..a_ven..Koeff. Ventilation (PK, S.541)
         &        0.308000, & !..b_ven..Koeff. Ventilation (PK, S.541)
         &        2.00,     & !..cap....Koeff. Kapazitaet
         &        1.0,      & !..vsedi_max
         &        0.0,      & !..vsedi_min
         &        null(),   & !..n pointer
         &        null(),   & !..q pointer
         &        null() )    !..rho_v pointer

    TYPE(particle), PARAMETER :: cloud_nue1mue1 = PARTICLE( &
         &        'cloud_nue1mue1',  & !.name...Bezeichnung der Partikelklasse
         &        1.000000, & !..nu.....Breiteparameter der Verteil.
         &        1.000000, & !..mu.....Exp.-parameter der Verteil.
         &        2.60d-10, & !..x_max..maximale Teilchenmasse D=80e-6m
         &        4.20d-15, & !..x_min..minimale Teilchenmasse D=2.e-6m
         &        1.24d-01, & !..a_geo..Koeff. Geometrie
         &        0.333333, & !..b_geo..Koeff. Geometrie = 1/3
         &        3.75d+05, & !..a_vel..Koeff. Fallgesetz
         &        0.666667, & !..b_vel..Koeff. Fallgesetz
         &        0.780000, & !..a_ven..Koeff. Ventilation (PK, S.541)
         &        0.308000, & !..b_ven..Koeff. Ventilation (PK, S.541)
         &        2.00,     & !..cap....Koeff. Kapazitaet
         &        1.0,      & !..vsedi_max
         &        0.0,      & !..vsedi_min
         &        null(),   & !..n pointer
         &        null(),   & !..q pointer
         &        null() )    !..rho_v pointer

    TYPE(particle), PARAMETER :: ice_cosmo5 =  particle( & ! iceCRY2test
         &        'ice_cosmo5', & !.name...Bezeichnung der Partikelklasse
         &        0.000000, & !..nu...e..Breiteparameter der Verteil.
         &        0.333333, & !..mu.....Exp.-parameter der Verteil.
         &        1.00d-05, & !..x_max..maximale Teilchenmasse D=???e-2m
         &        1.00d-12, & !..x_min..minimale Teilchenmasse D=200e-6m
         &        0.835000, & !..a_geo..Koeff. Geometrie
         &        0.390000, & !..b_geo..Koeff. Geometrie
         &        2.77d+01, & !..a_vel..Koeff. Fallgesetz
         &        0.215790, & !..b_vel..Koeff. Fallgesetz = 0.41/1.9
         &        0.780000, & !..a_ven..Koeff. Ventilation (PK, S.541)
         &        0.308000, & !..b_ven..Koeff. Ventilation (PK, S.541)
         &        2.0,      & !..cap....Koeff. Kapazitaet
         &        3.0,      & !..vsedi_max
         &        0.0,      & !..vsedi_min
         &        null(),   & !..n pointer
         &        null(),   & !..q pointer
         &        null() )    !..rho_v pointer

    TYPE(particle), PARAMETER :: snow_cosmo5 = particle( & ! nach Andy Heymsfield (CRYSTAL-FACE)
         &        'snow_cosmo5', & !.name...Bezeichnung der Partikelklasse
         &        0.000000, & !..nu.....Breiteparameter der Verteil.
         &        0.500000, & !..mu.....Exp.-parameter der Verteil.
         &        2.00d-05, & !..x_max..maximale Teilchenmasse D=???e-2m
         &        1.00d-10, & !..x_min..minimale Teilchenmasse D=200e-6m
         &        2.400000, & !..a_geo..Koeff. Geometrie
         &        0.455000, & !..b_geo..Koeff. Geometrie
         &        8.800000, & !..a_vel..Koeff. Fallgesetz
         &        0.150000, & !..b_vel..Koeff. Fallgesetz
         &        0.780000, & !..a_ven..Koeff. Ventilation (PK, S.541)
         &        0.308000, & !..b_ven..Koeff. Ventilation (PK, S.541)
         &        2.00,     & !..cap....Koeff. Kapazitaet
         &        3.0,      & !..vsedi_max
         &        0.1,      & !..vsedi_min
         &        null(),   & !..n pointer
         &        null(),   & !..q pointer
         &        null() )    !..rho_v pointer

    TYPE(particle), PARAMETER :: snowSBB =  particle(   & !
         &        'snowSBB',& !..name...Bezeichnung der Partikelklasse
         &        0.000000, & !..nu.....Breiteparameter der Verteil.
         &        0.500000, & !..mu.....Exp.-parameter der Verteil.
         &        2.00d-05, & !..x_max..maximale Teilchenmasse
         &        1.00d-10, & !..x_min..minimale Teilchenmasse
         &        5.130000, & !..a_geo..Koeff. Geometrie, x = 0.038*D**2
         &        0.500000, & !..b_geo..Koeff. Geometrie = 1/2
         &        8.294000, & !..a_vel..Koeff. Fallgesetz
         &        0.125000, & !..b_vel..Koeff. Fallgesetz
         &        0.780000, & !..a_ven..Koeff. Ventilation (PK, S.541)
         &        0.308000, & !..b_ven..Koeff. Ventilation (PK, S.541)
         &        2.00,     & !..cap....Koeff. Kapazitaet
         &        3.0,      & !..vsedi_max
         &        0.1,      & !..vsedi_min
         &        null(),   & !..n pointer
         &        null(),   & !..q pointer
         &        null() )    !..rho_v pointer

    TYPE(particle), PARAMETER :: rainULI = particle( & ! Blahak, v=v(x) gefittet 6.9.2005
         &        'rainULI', & !..name
         &        0.000000,  & !..nu
         &        0.333333,  & !..mu
         &        3.00d-06,  & !..x_max
         &        2.60d-10,  & !..x_min
         &        1.24d-01,  & !..a_geo
         &        0.333333,  & !..b_geo
         &        114.0137,  & !..a_vel
         &        0.234370,  & !..b_vel
         &        0.780000,  & !..a_ven
         &        0.308000,  & !..b_ven
         &        2.00,      & !..cap
         &        20.0,      & !..vsedi_max
         &        0.1,       & !..vsedi_min
         &        null(),    & !..n pointer
         &        null(),    & !..q pointer
         &        null() )     !..rho_v pointer

    TYPE(particle), PARAMETER :: rainSBB = particle( &
         &        'rainSBB', & !..name
         &        0.000000,  & !..nu
         &        0.333333,  & !..mu
         &        3.00d-06,  & !..x_max
         &        2.60d-10,  & !..x_min
         &        1.24d-01,  & !..a_geo
         &        0.333333,  & !..b_geo
         &        114.0137,  & !..a_vel
         &        0.234370,  & !..b_vel
         &        0.780000,  & !..a_ven
         &        0.308000,  & !..b_ven
         &        2.000000,  & !..cap
         &        2.000d+1,  & !..vsedi_max
         &        0.1,       & !..vsedi_min
         &        null(),    & !..n pointer
         &        null(),    & !..q pointer
         &        null() )     !..rho_v pointer

    cloud   = cloud_nue1mue1
    rain    = rainSBB
    ice     = ice_cosmo5
    snow    = snowSBB
    graupel = graupelhail_cosmo5
    hail    = hail_cosmo5

  END SUBROUTINE init_2mom_scheme

  SUBROUTINE init_2mom_scheme_once(cloud_type)
    INTEGER, INTENT(in)  :: cloud_type

    CHARACTER(len=*), PARAMETER :: routine = 'init_2mom_scheme_once'
    REAL(wp), DIMENSION(1:1) :: q_r,x_r,q_c,vn_rain_min, vq_rain_min, vn_rain_max, vq_rain_max
    TYPE(particle) :: cloud, rain, ice, snow, graupel, hail

    CALL init_2mom_scheme(cloud,rain,ice,snow,graupel,hail)

    ice_typ   = cloud_type/1000           ! (0) no ice, (1) no hail (2) with hail
    nuc_i_typ = MOD(cloud_type/100,10)    ! choice of ice nucleation, see ice_nucleation_homhet()
    nuc_c_typ = MOD(cloud_type/10,10)     ! choice of CCN assumptions, see cloud_nucleation()
    auto_typ  = MOD(cloud_type,10)        ! choice of warm rain scheme, see clouds_twomoment()

    rain_coeffs = rainSBBcoeffs

    CALL message(TRIM(routine), "calculate run-time coefficients")
    WRITE(txt,'(A,I10)')   "  cloud_type = ",cloud_type ; CALL message(routine,TRIM(txt))

    ! initialize bulk sedimentation velocities
    ! calculates coeff_alfa_n, coeff_alfa_q, and coeff_lambda
    call init_sedi_vel(ice,ice_coeffs)
    call init_sedi_vel(snow,snow_coeffs)
    call init_sedi_vel(graupel,graupel_coeffs)
    call init_sedi_vel(hail,hail_coeffs)

    ! look-up table and parameters for rain_freeze_gamlook
    rain_nm1 = (rain%nu+1.0)/rain%mu
    rain_nm2 = (rain%nu+2.0)/rain%mu
    rain_nm3 = (rain%nu+3.0)/rain%mu
    CALL incgfct_lower_lookupcreate(rain_nm1, rain_ltable1, nlookup, nlookuphr_dummy)
    CALL incgfct_lower_lookupcreate(rain_nm2, rain_ltable2, nlookup, nlookuphr_dummy)
    CALL incgfct_lower_lookupcreate(rain_nm3, rain_ltable3, nlookup, nlookuphr_dummy)
    rain_g1 = rain_ltable1%igf(rain_ltable1%n) ! ordinary gamma function of nm1 is the last value in table 1
    rain_g2 = rain_ltable2%igf(rain_ltable2%n) ! ordinary gamma function of nm2 is the last value in table 2

    ! table and parameters for graupel_hail_conv_wet_gamlook
    graupel_nm1 = (graupel%nu+1.0)/graupel%mu
    graupel_nm2 = (graupel%nu+2.0)/graupel%mu
    CALL incgfct_lower_lookupcreate(graupel_nm1, graupel_ltable1, nlookup, nlookuphr_dummy)
    CALL incgfct_lower_lookupcreate(graupel_nm2, graupel_ltable2, nlookup, nlookuphr_dummy)
    graupel_g1 = graupel_ltable1%igf(graupel_ltable1%n) ! ordinary gamma function of nm1 is the last value in table 1
    graupel_g2 = graupel_ltable2%igf(graupel_ltable2%n) ! ordinary gamma function of nm2 is the last value in table 2

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
    q_r = 1.0e-3_wp
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
    x_r = rain%x_min ; CALL sedi_vel_rain(rain,rain_coeffs,q_r,x_r,vn_rain_min,vq_rain_min,1,1)
    x_r = rain%x_max ; CALL sedi_vel_rain(rain,rain_coeffs,q_r,x_r,vn_rain_max,vq_rain_max,1,1)
    WRITE(txt,'(A)')       "    out-of-cloud: " ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vn_rain_min  = ",vn_rain_min ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vn_rain_max  = ",vn_rain_max ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vq_rain_min  = ",vq_rain_min ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vq_rain_max  = ",vq_rain_max ; CALL message(routine,TRIM(txt))
    q_c = 1e-3_wp
    x_r = rain%x_min ; CALL sedi_vel_rain(rain,rain_coeffs,q_r,x_r,vn_rain_min,vq_rain_min,1,1,q_c)
    x_r = rain%x_max ; CALL sedi_vel_rain(rain,rain_coeffs,q_r,x_r,vn_rain_max,vq_rain_max,1,1,q_c)
    WRITE(txt,'(A)')       "    in-cloud: " ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vn_rain_min  = ",vn_rain_min ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vn_rain_max  = ",vn_rain_max ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vq_rain_min  = ",vq_rain_min ; CALL message(routine,TRIM(txt))
    WRITE(txt,'(A,D10.3)') "     vq_rain_max  = ",vq_rain_max ; CALL message(routine,TRIM(txt))

    ! initialization for snow_cloud_riming
    scr_params%delta_n_aa = coll_delta_11(snow,cloud,0)
    scr_params%delta_n_ab = coll_delta_12(snow,cloud,0)
    scr_params%delta_n_bb = coll_delta_22(snow,cloud,0)
    scr_params%delta_q_aa = coll_delta_11(snow,cloud,0)
    scr_params%delta_q_ab = coll_delta_12(snow,cloud,1)
    scr_params%delta_q_bb = coll_delta_22(snow,cloud,1)

    scr_params%theta_n_aa = coll_theta_11(snow,cloud,0)
    scr_params%theta_n_ab = coll_theta_12(snow,cloud,0)
    scr_params%theta_n_bb = coll_theta_22(snow,cloud,0)
    scr_params%theta_q_aa = coll_theta_11(snow,cloud,0)
    scr_params%theta_q_ab = coll_theta_12(snow,cloud,1)
    scr_params%theta_q_bb = coll_theta_22(snow,cloud,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "   a_snow      = ",snow%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   b_snow      = ",snow%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   alf_snow    = ",snow%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   bet_snow    = ",snow%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   a_cloud    = ",cloud%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   b_cloud    = ",cloud%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   alf_cloud  = ",cloud%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   bet_cloud  = ",cloud%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_n_ss = ",scr_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_n_sc = ",scr_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_n_cc = ",scr_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_n_ss = ",scr_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_n_sc = ",scr_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_n_cc = ",scr_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_q_ss = ",scr_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_q_sc = ",scr_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_q_cc = ",scr_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_q_ss = ",scr_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_q_sc = ",scr_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_q_cc = ",scr_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! coefficients for snow_rain_riming
    srr_params%delta_n_aa = coll_delta_11(snow,rain,0)
    srr_params%delta_n_ab = coll_delta_12(snow,rain,0)
    srr_params%delta_n_bb = coll_delta_22(snow,rain,0)
    srr_params%delta_q_aa = coll_delta_11(snow,rain,1) ! mass weighted
    srr_params%delta_q_ab = coll_delta_12(snow,rain,1)
    srr_params%delta_q_ba = coll_delta_12(rain,snow,1)
    srr_params%delta_q_bb = coll_delta_22(snow,rain,1)

    srr_params%theta_n_aa = coll_theta_11(snow,rain,0)
    srr_params%theta_n_ab = coll_theta_12(snow,rain,0)
    srr_params%theta_n_bb = coll_theta_22(snow,rain,0)
    srr_params%theta_q_aa = coll_theta_11(snow,rain,1) ! mass weighted
    srr_params%theta_q_ab = coll_theta_12(snow,rain,1)
    srr_params%theta_q_ba = coll_theta_12(rain,snow,1)
    srr_params%theta_q_bb = coll_theta_22(snow,rain,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "    a_snow     = ",snow%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_snow     = ",snow%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    alf_snow   = ",snow%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    bet_snow   = ",snow%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_rain     = ",rain%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_rain     = ",rain%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    alf_rain   = ",rain%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    bet_rain   = ",rain%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_ss = ",srr_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_sr = ",srr_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_rr = ",srr_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_ss = ",srr_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_sr = ",srr_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_rr = ",srr_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_ss = ",srr_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_sr = ",srr_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_rs = ",srr_params%delta_q_ba ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_rr = ",srr_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_ss = ",srr_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_sr = ",srr_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_rs = ",srr_params%theta_q_ba ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_rr = ",srr_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! ice rain riming parameters
    irr_params%delta_n_aa = coll_delta_11(ice,rain,0)
    irr_params%delta_n_ab = coll_delta_12(ice,rain,0)
    irr_params%delta_n_bb = coll_delta_22(ice,rain,0)
    irr_params%delta_q_aa = coll_delta_11(ice,rain,1) ! here mass weighted
    irr_params%delta_q_ab = coll_delta_12(ice,rain,1)
    irr_params%delta_q_ba = coll_delta_12(rain,ice,1)
    irr_params%delta_q_bb = coll_delta_22(ice,rain,1)

    irr_params%theta_n_aa = coll_theta_11(ice,rain,0)
    irr_params%theta_n_ab = coll_theta_12(ice,rain,0)
    irr_params%theta_n_bb = coll_theta_22(ice,rain,0)
    irr_params%theta_q_aa = coll_theta_11(ice,rain,1) ! here mass weighted
    irr_params%theta_q_ab = coll_theta_12(ice,rain,1)
    irr_params%theta_q_ba = coll_theta_12(rain,ice,1)
    irr_params%theta_q_bb = coll_theta_22(ice,rain,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "     a_ice      = ",ice%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     b_ice      = ",ice%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     alf_ice    = ",ice%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     bet_ice    = ",ice%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     a_rain    = ",rain%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     b_rain    = ",rain%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     alf_rain  = ",rain%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     bet_rain  = ",rain%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_n_ii = ", irr_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_n_ir = ", irr_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_n_rr = ", irr_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_ii = ", irr_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_ir = ", irr_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_rr = ", irr_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_ii = ", irr_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_ir = ", irr_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_ri = ", irr_params%delta_q_ba ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_rr = ", irr_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_ii = ", irr_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_ir = ", irr_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_ri = ", irr_params%theta_q_ba ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_rr = ", irr_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! ice cloud riming parameter setup
    icr_params%delta_n_aa = coll_delta_11(ice,cloud,0)
    icr_params%delta_n_ab = coll_delta_12(ice,cloud,0)
    icr_params%delta_n_bb = coll_delta_22(ice,cloud,0)
    icr_params%delta_q_aa = coll_delta_11(ice,cloud,0)
    icr_params%delta_q_ab = coll_delta_12(ice,cloud,1)
    icr_params%delta_q_bb = coll_delta_22(ice,cloud,1)

    icr_params%theta_n_aa = coll_theta_11(ice,cloud,0)
    icr_params%theta_n_ab = coll_theta_12(ice,cloud,0)
    icr_params%theta_n_bb = coll_theta_22(ice,cloud,0)
    icr_params%theta_q_aa = coll_theta_11(ice,cloud,0)
    icr_params%theta_q_ab = coll_theta_12(ice,cloud,1)
    icr_params%theta_q_bb = coll_theta_22(ice,cloud,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "    a_ice      = ",ice%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_ice      = ",ice%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    alf_ice    = ",ice%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    bet_ice    = ",ice%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_cloud    = ",cloud%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_cloud    = ",cloud%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    alf_cloud  = ",cloud%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    bet_cloud  = ",cloud%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_ii = ", icr_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_ic = ", icr_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_cc = ", icr_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_ii = ", icr_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_ic = ", icr_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_cc = ", icr_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_ii = ", icr_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_ic = ", icr_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_cc = ", icr_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_ii = ", icr_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_ic = ", icr_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_cc = ", icr_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! hail rain riming
    hrr_params%delta_n_aa = coll_delta_11(hail,rain,0)
    hrr_params%delta_n_ab = coll_delta_12(hail,rain,0)
    hrr_params%delta_n_bb = coll_delta_22(hail,rain,0)
    hrr_params%delta_q_aa = coll_delta_11(hail,rain,0)
    hrr_params%delta_q_ab = coll_delta_12(hail,rain,1)
    hrr_params%delta_q_bb = coll_delta_22(hail,rain,1)

    hrr_params%theta_n_aa = coll_theta_11(hail,rain,0)
    hrr_params%theta_n_ab = coll_theta_12(hail,rain,0)
    hrr_params%theta_n_bb = coll_theta_22(hail,rain,0)
    hrr_params%theta_q_aa = coll_theta_11(hail,rain,0)
    hrr_params%theta_q_ab = coll_theta_12(hail,rain,1)
    hrr_params%theta_q_bb = coll_theta_22(hail,rain,1)

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
      WRITE(txt,'(A,D10.3)') "     delta_n_hh = ",hrr_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_n_hr = ",hrr_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_n_rr = ",hrr_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_hh = ",hrr_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_hr = ",hrr_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_rr = ",hrr_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_hh = ",hrr_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_hr = ",hrr_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_rr = ",hrr_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_hh = ",hrr_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_hr = ",hrr_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_rr = ",hrr_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! graupel rain riming parameter setup
    grr_params%delta_n_aa = coll_delta_11(graupel,rain,0)
    grr_params%delta_n_ab = coll_delta_12(graupel,rain,0)
    grr_params%delta_n_bb = coll_delta_22(graupel,rain,0)
    grr_params%delta_q_aa = coll_delta_11(graupel,rain,0)
    grr_params%delta_q_ab = coll_delta_12(graupel,rain,1)
    grr_params%delta_q_ba = coll_delta_12(rain,graupel,1)
    grr_params%delta_q_bb = coll_delta_22(graupel,rain,1)

    grr_params%theta_n_aa = coll_theta_11(graupel,rain,0)
    grr_params%theta_n_ab = coll_theta_12(graupel,rain,0)
    grr_params%theta_n_bb = coll_theta_22(graupel,rain,0)
    grr_params%theta_q_aa = coll_theta_11(graupel,rain,0)
    grr_params%theta_q_ab = coll_theta_12(graupel,rain,1)
    grr_params%theta_q_ba = coll_theta_12(rain,graupel,1)
    grr_params%theta_q_bb = coll_theta_22(graupel,rain,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "     delta_n_gg = ",grr_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_n_gr = ",grr_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_n_rr = ",grr_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_gg = ",grr_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_gr = ",grr_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_rr = ",grr_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_gg = ",grr_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_gr = ",grr_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_rg = ",grr_params%delta_q_ba ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_rr = ",grr_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_gg = ",grr_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_gr = ",grr_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_rg = ",grr_params%theta_q_ba ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_rr = ",grr_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! hail cloud riming parameter setup
    hcr_params%delta_n_aa = coll_delta_11(hail,cloud,0)
    hcr_params%delta_n_ab = coll_delta_12(hail,cloud,0)
    hcr_params%delta_n_bb = coll_delta_22(hail,cloud,0)
    hcr_params%delta_q_aa = coll_delta_11(hail,cloud,0)
    hcr_params%delta_q_ab = coll_delta_12(hail,cloud,1)
    hcr_params%delta_q_bb = coll_delta_22(hail,cloud,1)

    hcr_params%theta_n_aa = coll_theta_11(hail,cloud,0)
    hcr_params%theta_n_ab = coll_theta_12(hail,cloud,0)
    hcr_params%theta_n_bb = coll_theta_22(hail,cloud,0)
    hcr_params%theta_q_aa = coll_theta_11(hail,cloud,0)
    hcr_params%theta_q_ab = coll_theta_12(hail,cloud,1)
    hcr_params%theta_q_bb = coll_theta_22(hail,cloud,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "     delta_n_hh  = ", hcr_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_n_hc  = ", hcr_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_n_cc  = ", hcr_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_hh  = ", hcr_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_hc  = ", hcr_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_n_cc  = ", hcr_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_hh  = ", hcr_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_hc  = ", hcr_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     delta_q_cc  = ", hcr_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_hh  = ", hcr_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_hc  = ", hcr_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_cc  = ", hcr_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! graupel cloud riming parameters
    gcr_params%delta_n_aa = coll_delta_11(graupel,cloud,0)
    gcr_params%delta_n_ab = coll_delta_12(graupel,cloud,0)
    gcr_params%delta_n_bb = coll_delta_22(graupel,cloud,0)
    gcr_params%delta_q_aa = coll_delta_11(graupel,cloud,0)
    gcr_params%delta_q_ab = coll_delta_12(graupel,cloud,1)
    gcr_params%delta_q_bb = coll_delta_22(graupel,cloud,1)

    gcr_params%theta_n_aa = coll_theta_11(graupel,cloud,0)
    gcr_params%theta_n_ab = coll_theta_12(graupel,cloud,0)
    gcr_params%theta_n_bb = coll_theta_22(graupel,cloud,0)
    gcr_params%theta_q_aa = coll_theta_11(graupel,cloud,0)
    gcr_params%theta_q_ab = coll_theta_12(graupel,cloud,1)
    gcr_params%theta_q_bb = coll_theta_22(graupel,cloud,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "    delta_n_gg  = ", gcr_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_gc  = ", gcr_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_cc  = ", gcr_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_gg  = ", gcr_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_gc  = ", gcr_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_cc  = ", gcr_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_gg  = ", gcr_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_gc  = ", gcr_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_cc  = ", gcr_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_gg  = ", gcr_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_gc  = ", gcr_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_cc  = ", gcr_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! snow ice collection parameters setup
    sic_params%delta_n_aa = coll_delta_11(snow,ice,0)
    sic_params%delta_n_ab = coll_delta_12(snow,ice,0)
    sic_params%delta_n_bb = coll_delta_22(snow,ice,0)
    sic_params%delta_q_aa = coll_delta_11(snow,ice,0)
    sic_params%delta_q_ab = coll_delta_12(snow,ice,1)
    sic_params%delta_q_bb = coll_delta_22(snow,ice,1)

    sic_params%theta_n_aa = coll_theta_11(snow,ice,0)
    sic_params%theta_n_ab = coll_theta_12(snow,ice,0)
    sic_params%theta_n_bb = coll_theta_22(snow,ice,0)
    sic_params%theta_q_aa = coll_theta_11(snow,ice,0)
    sic_params%theta_q_ab = coll_theta_12(snow,ice,1)
    sic_params%theta_q_bb = coll_theta_22(snow,ice,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "   delta_n_ss = ", sic_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_n_si = ", sic_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_n_ii = ", sic_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_n_ss = ", sic_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_n_si = ", sic_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_n_ii = ", sic_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_q_ss = ", sic_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_q_si = ", sic_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_q_ii = ", sic_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_q_ss = ", sic_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_q_si = ", sic_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_q_ii = ", sic_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! hail ice collection parameter setup
    hic_params%delta_n_aa = coll_delta_11(hail,ice,0)
    hic_params%delta_n_ab = coll_delta_12(hail,ice,0)
    hic_params%delta_n_bb = coll_delta_22(hail,ice,0)
    hic_params%delta_q_aa = coll_delta_11(hail,ice,0)
    hic_params%delta_q_ab = coll_delta_12(hail,ice,1)
    hic_params%delta_q_bb = coll_delta_22(hail,ice,1)

    hic_params%theta_n_aa = coll_theta_11(hail,ice,0)
    hic_params%theta_n_ab = coll_theta_12(hail,ice,0)
    hic_params%theta_n_bb = coll_theta_22(hail,ice,0)
    hic_params%theta_q_aa = coll_theta_11(hail,ice,0)
    hic_params%theta_q_ab = coll_theta_12(hail,ice,1)
    hic_params%theta_q_bb = coll_theta_22(hail,ice,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "    delta_n_hh = ", hic_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_hi = ", hic_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_ii = ", hic_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_hh = ", hic_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_hi = ", hic_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_ii = ", hic_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_hh = ", hic_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_hi = ", hic_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_ii = ", hic_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_hh = ", hic_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_hi = ", hic_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_ii = ", hic_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! graupel ice collection parameter setup
    gic_params%delta_n_aa = coll_delta_11(graupel,ice,0)
    gic_params%delta_n_ab = coll_delta_12(graupel,ice,0)
    gic_params%delta_n_bb = coll_delta_22(graupel,ice,0)
    gic_params%delta_q_aa = coll_delta_11(graupel,ice,0)
    gic_params%delta_q_ab = coll_delta_12(graupel,ice,1)
    gic_params%delta_q_bb = coll_delta_22(graupel,ice,1)

    gic_params%theta_n_aa = coll_theta_11(graupel,ice,0)
    gic_params%theta_n_ab = coll_theta_12(graupel,ice,0)
    gic_params%theta_n_bb = coll_theta_22(graupel,ice,0)
    gic_params%theta_q_aa = coll_theta_11(graupel,ice,0)
    gic_params%theta_q_ab = coll_theta_12(graupel,ice,1)
    gic_params%theta_q_bb = coll_theta_22(graupel,ice,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "   a_graupel   = ",graupel%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   b_graupel   = ",graupel%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   alf_graupel = ",graupel%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   bet_graupel = ",graupel%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   a_ice       = ",ice%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   b_ice       = ",ice%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   alf_ice     = ",ice%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   bet_ice     = ",ice%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_n_gg  = ", gic_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_n_gi  = ", gic_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_n_ii  = ", gic_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_n_gg  = ", gic_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_n_gi  = ", gic_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_n_ii  = ", gic_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_q_gg  = ", gic_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_q_gi  = ", gic_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   delta_q_ii  = ", gic_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_q_gg  = ", gic_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_q_gi  = ", gic_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   theta_q_ii  = ", gic_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! hail snow collection parameter setup
    hsc_params%delta_n_aa = coll_delta_11(hail,snow,0)
    hsc_params%delta_n_ab = coll_delta_12(hail,snow,0)
    hsc_params%delta_n_bb = coll_delta_22(hail,snow,0)
    hsc_params%delta_q_aa = coll_delta_11(hail,snow,0)
    hsc_params%delta_q_ab = coll_delta_12(hail,snow,1)
    hsc_params%delta_q_bb = coll_delta_22(hail,snow,1)

    hsc_params%theta_n_aa = coll_theta_11(hail,snow,0)
    hsc_params%theta_n_ab = coll_theta_12(hail,snow,0)
    hsc_params%theta_n_bb = coll_theta_22(hail,snow,0)
    hsc_params%theta_q_aa = coll_theta_11(hail,snow,0)
    hsc_params%theta_q_ab = coll_theta_12(hail,snow,1)
    hsc_params%theta_q_bb = coll_theta_22(hail,snow,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "    delta_n_hh = ", hsc_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_hs = ", hsc_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_ss = ", hsc_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_hh = ", hsc_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_hs = ", hsc_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_ss = ", hsc_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_hh = ", hsc_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_hs = ", hsc_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_ss = ", hsc_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_hh = ", hsc_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_hs = ", hsc_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_ss = ", hsc_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! graupel snow collection parameter setup
    gsc_params%delta_n_aa = coll_delta_11(graupel,snow,0)
    gsc_params%delta_n_ab = coll_delta_12(graupel,snow,0)
    gsc_params%delta_n_bb = coll_delta_22(graupel,snow,0)
    gsc_params%delta_q_aa = coll_delta_11(graupel,snow,0)
    gsc_params%delta_q_ab = coll_delta_12(graupel,snow,1)
    gsc_params%delta_q_bb = coll_delta_22(graupel,snow,1)

    gsc_params%theta_n_aa = coll_theta_11(graupel,snow,0)
    gsc_params%theta_n_ab = coll_theta_12(graupel,snow,0)
    gsc_params%theta_n_bb = coll_theta_22(graupel,snow,0)
    gsc_params%theta_q_aa = coll_theta_11(graupel,snow,0)
    gsc_params%theta_q_ab = coll_theta_12(graupel,snow,1)
    gsc_params%theta_q_bb = coll_theta_22(graupel,snow,1)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "    delta_n_gg = ", gsc_params%delta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_gs = ", gsc_params%delta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_ss = ", gsc_params%delta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_gg = ", gsc_params%theta_n_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_gs = ", gsc_params%theta_n_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_ss = ", gsc_params%theta_n_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_gg = ", gsc_params%delta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_gs = ", gsc_params%delta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_q_ss = ", gsc_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_gg = ", gsc_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_gs = ", gsc_params%theta_q_ab ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_q_ss = ", gsc_params%theta_q_bb ; CALL message(routine,TRIM(txt))
    END IF

    ! vapor ice deposition params
    vid_params%c = 1.0 / ice%cap
    vid_params%a_f = vent_coeff_a(ice,1)
    vid_params%b_f = vent_coeff_b(ice,1) * N_sc**n_f / SQRT(nu_l)
    IF (isdebug) THEN
      WRITE (txt,'(A,D10.3)') "    a_geo   = ",ice%a_geo ; CALL message(routine,TRIM(txt))
      WRITE (txt,'(A,D10.3)') "    b_geo   = ",ice%b_geo ; CALL message(routine,TRIM(txt))
      WRITE (txt,'(A,D10.3)') "    a_vel   = ",ice%a_vel ; CALL message(routine,TRIM(txt))
      WRITE (txt,'(A,D10.3)') "    b_vel   = ",ice%b_vel ; CALL message(routine,TRIM(txt))
      WRITE (txt,'(A,D10.3)') "    c_i     = ", vid_params%c ; CALL message(routine,TRIM(txt))
      WRITE (txt,'(A,D10.3)') "    a_f     = ", vid_params%a_f ; CALL message(routine,TRIM(txt))
      WRITE (txt,'(A,D10.3)') "    b_f     = ", vid_params%b_f ; CALL message(routine,TRIM(txt))
    END IF

    ! vapor graupel deposition parameter setup
    vgd_params%c = 1.0 / graupel%cap
    vgd_params%a_f = vent_coeff_a(graupel,1)
    vgd_params%b_f = vent_coeff_b(graupel,1) * N_sc**n_f / SQRT(nu_l)
    IF (isdebug) THEN
      WRITE(txt,*) "  vapor_deposition_graupel: " ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_geo = ",graupel%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_geo = ",graupel%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_vel = ",graupel%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_vel = ",graupel%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    c_g   = ", vgd_params%c ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_f   = ", vgd_params%a_f ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_f   = ", vgd_params%b_f ; CALL message(routine,TRIM(txt))
    END IF

    ! vapor hail deposition parameter setup
    vhd_params%c = 1.0 / hail%cap
    vhd_params%a_f = vent_coeff_a(hail,1)
    vhd_params%b_f = vent_coeff_b(hail,1) * N_sc**n_f / sqrt(nu_l)
    IF (isdebug) THEN
      WRITE(txt,*) "  vapor_deposition_hail: "
      WRITE(txt,'(A,D10.3)') "    a_geo = ",hail%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_geo = ",hail%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_vel = ",hail%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_vel = ",hail%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    c_h   = ",vhd_params%c ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_f   = ",vhd_params%a_f ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_f   = ",vhd_params%b_f ; CALL message(routine,TRIM(txt))
    END IF

    ! vapor snow deposition parameter setup
    vsd_params%c = 1.0 / snow%cap
    vsd_params%a_f = vent_coeff_a(snow,1)
    vsd_params%b_f = vent_coeff_b(snow,1) * N_sc**n_f / sqrt(nu_l)
    IF (isdebug) THEN
      WRITE(txt,*) "  vapor_depositiosnow%n: "
      WRITE(txt,'(A,D10.3)') "    a_geo = ",snow%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_geo = ",snow%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_vel = ",snow%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_vel = ",snow%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    c_s   = ",vsd_params%c ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_f   = ",vsd_params%a_f ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     b_f  = ",vsd_params%b_f ; CALL message(routine,TRIM(txt))
    END IF

    ! graupel evaporation parameters
    ge_params%a_f = vent_coeff_a(graupel,1)
    ge_params%b_f = vent_coeff_b(graupel,1) * N_sc**n_f / SQRT(nu_l)
    ge_params%c = 1.0_wp / graupel%cap
    IF (isdebug) THEN
      WRITE(txt,'(A)') "graupel_evaporation:"  ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     a_f = ",ge_params%a_f ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     b_f = ",ge_params%b_f ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     c_g = ",ge_params%c ; CALL message(routine,TRIM(txt))
    END IF

    ! hail evaporation parameter setup
    he_params%a_f = vent_coeff_a(hail,1)
    he_params%b_f = vent_coeff_b(hail,1) * N_sc**n_f / sqrt(nu_l)
    he_params%c = 1.0_wp / hail%cap
    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "hail_evaporation:" ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     a_f = ",he_params%a_f   ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     b_f = ",he_params%b_f   ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     c_h = ",he_params%c   ; CALL message(routine,TRIM(txt))
    END IF

    ! snow evaporation parameter setup
    se_params%a_f = vent_coeff_a(snow,1)
    se_params%b_f = vent_coeff_b(snow,1) * N_sc**n_f / sqrt(nu_l)
    se_params%c = 1.0_wp / snow%cap
    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "snow_evaporation:" ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     a_f = ",se_params%a_f ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     b_f = ",se_params%b_f ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     c_s = ",se_params%c ; CALL message(routine,TRIM(txt))
    END IF

    ! graupel melting parameter setup
    gm_params%a_vent = vent_coeff_a(graupel,1)
    gm_params%b_vent = vent_coeff_b(graupel,1) * N_sc**n_f / SQRT(nu_l)
    IF (isdebug) THEN
      WRITE(txt,*) " graupel_melting "
      WRITE(txt,'(A,D10.3)') "   a_geo    = ",graupel%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   b_geo    = ",graupel%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   a_vel    = ",graupel%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   b_vel    = ",graupel%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   a_ven    = ",graupel%a_ven ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   b_ven    = ",graupel%b_ven ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   a_vent = ",gm_params%a_vent ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "   b_vent = ",gm_params%b_vent ; CALL message(routine,TRIM(txt))
    ENDIF

    ! hail melting parameters
    hm_params%a_vent = vent_coeff_a(hail,1)
    hm_params%b_vent = vent_coeff_b(hail,1) * N_sc**n_f / SQRT(nu_l)
    IF (isdebug) THEN
      WRITE(txt,*) " hail_melting: " ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_vent = ",hm_params%a_vent ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_vent = ",hm_params%b_vent ; CALL message(routine,TRIM(txt))
    ENDIF

    ! snow melting parameter setup
    sm_params%a_vent = vent_coeff_a(snow,1)
    sm_params%b_vent = vent_coeff_b(snow,1) * N_sc**n_f / SQRT(nu_l)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "    a_geo  = ",snow%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_geo  = ",snow%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_vel  = ",snow%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_vel  = ",snow%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_ven  = ",snow%a_ven ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_ven  = ",snow%b_ven ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    a_vent = ",sm_params%a_vent ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_vent = ",sm_params%b_vent ; CALL message(routine,TRIM(txt))
    ENDIF

    CALL setup_graupel_selfcollection(graupel)
    CALL setup_snow_selfcollection(snow)

  END SUBROUTINE init_2mom_scheme_once

  !*******************************************************************************
  ! Main subroutine of the two-moment microphysics
  !
  ! All individual processes are called in sequence and do their own time
  ! integration for the mass and number densities, i.e., Marchuk-type operator
  ! splitting. Temperature is only updated in the driver.
  !*******************************************************************************

  SUBROUTINE clouds_twomoment(ik_slice, dt, use_prog_in, atmo, &
       cloud, rain, ice, snow, graupel, hail, n_inact, n_cn, n_inpot)

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    LOGICAL, INTENT(in) :: use_prog_in
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud, rain, ice, snow, graupel, hail

    REAL(wp), DIMENSION(:,:) :: n_inact

    ! optional arguments for version with prognostic CCN and IN
    ! (for nuc_c_typ > 0 and  use_prog_in=true)
    REAL(wp), DIMENSION(:,:) ,OPTIONAL :: n_inpot, n_cn

    REAL(wp), DIMENSION(size(cloud%n,1),size(cloud%n,2)) :: dep_rate_ice, dep_rate_snow

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER :: k, i

    IF (isdebug) CALL message(TRIM(routine),"clouds_twomoment start")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    dep_rate_ice(:,:)  = 0.0_wp
    dep_rate_snow(:,:) = 0.0_wp

    IF (isdebug) CALL message(TRIM(routine),'cloud_nucleation')
    IF (nuc_c_typ .EQ. 0) THEN

       IF (isdebug) CALL message(TRIM(routine),'  ... force constant cloud droplet number')
       cloud%n(:,:) = qnc_const

    ELSEIF (nuc_c_typ < 6) THEN
       IF (isdebug) CALL message(TRIM(routine),'  ... Hande et al CCN activation')
       IF (PRESENT(n_cn)) THEN
          CALL finish(TRIM(routine),&
               & 'Error in two_moment_mcrph: Hande et al activation not supported for progn. aerosol')
       ELSE
          CALL ccn_activation_hdcp2(ik_slice, atmo, cloud)
       END IF
    ELSE
       IF (isdebug) CALL message(TRIM(routine), &
            & '  ... CCN activation using look-up tables according to Segal& Khain')
       IF (PRESENT(n_cn)) THEN
          CALL ccn_activation_sk(ik_slice, atmo, cloud, n_cn)
       ELSE
          CALL finish(TRIM(routine),&
               & 'Error in two_moment_mcrph: Segal and Khain activation only supported for progn. aerosol')
       END IF
    END IF

    IF (ischeck) CALL check(ik_slice, 'start',cloud,rain,ice,snow,graupel,hail)

    IF (nuc_c_typ.ne.0) THEN
      DO k=kstart,kend
        DO i=istart,iend
          cloud%n(i,k) = MAX(cloud%n(i,k), cloud%q(i,k) / cloud%x_max)
          cloud%n(i,k) = MIN(cloud%n(i,k), cloud%q(i,k) / cloud%x_min)
        END DO
      END DO
    END IF

    IF (ice_typ > 0) THEN

       ! homogeneous and heterogeneous ice nucleation
       ! FIXME: optional arguments passed conditionally
       IF (present(n_inpot)) THEN
          CALL ice_nucleation_homhet(ik_slice, use_prog_in, atmo, cloud, ice, snow, n_inact, n_inpot)
       ELSE
          CALL ice_nucleation_homhet(ik_slice, use_prog_in, atmo, cloud, ice, snow, n_inact)
       END IF

       ! homogeneous freezing of cloud droplets
       CALL cloud_freeze(ik_slice, dt, atmo, cloud, ice)
       IF (ischeck) CALL check(ik_slice, 'cloud_freeze', cloud, rain, ice, snow, graupel,hail)

       DO k=kstart,kend
        DO i=istart,iend
          ice%n(i,k) = MIN(ice%n(i,k), ice%q(i,k)/ice%x_min)
          ice%n(i,k) = MAX(ice%n(i,k), ice%q(i,k)/ice%x_max)
        END DO
       END DO
       IF (ischeck) CALL check(ik_slice, 'ice nucleation',cloud,rain,ice,snow,graupel,hail)

       ! depositional growth of all ice particles
       ! ( store deposition rate of ice and snow for conversion calculation in
       !   ice_riming and snow_riming )
       CALL vapor_dep_relaxation(ik_slice, dt,atmo,ice,snow,graupel,hail,dep_rate_ice,dep_rate_snow)
       IF (ischeck) CALL check(ik_slice, 'vapor_dep_relaxation',cloud,rain,ice,snow,graupel,hail)

       if (.true.) then

       ! ice-ice collisions
       CALL ice_selfcollection(ik_slice, dt, atmo, ice, snow)
       CALL snow_selfcollection(ik_slice, dt, atmo, snow)
       CALL snow_ice_collection(ik_slice, dt, atmo, ice, snow)
       IF (ischeck) CALL check(ik_slice, 'ice and snow collection',cloud,rain,ice,snow,graupel,hail)

       CALL graupel_selfcollection(ik_slice, dt, atmo, graupel)
       CALL graupel_ice_collection(ik_slice, dt, atmo, ice, graupel)
       CALL graupel_snow_collection(ik_slice, dt, atmo, snow, graupel)
       IF (ischeck) CALL check(ik_slice, 'graupel collection',cloud,rain,ice,snow,graupel,hail)

       IF (ice_typ > 1) THEN

          ! conversion of graupel to hail in wet growth regime
          CALL graupel_hail_conv_wet_gamlook(ik_slice, atmo, graupel, cloud, &
               rain, ice, snow, hail)
          IF (ischeck) CALL check(ik_slice, 'graupel_hail_conv_wet_gamlook',cloud,rain,ice,snow,graupel,hail)

          ! hail collisions
          CALL hail_ice_collection(ik_slice, dt, atmo, ice, hail)    ! Important?
          CALL hail_snow_collection(ik_slice, dt, atmo, snow, hail)
          IF (ischeck) CALL check(ik_slice, 'hail collection',cloud,rain,ice,snow,graupel,hail)
       END IF

       ! riming of ice with cloud droplets and rain drops, and conversion to graupel
       CALL ice_riming(ik_slice, dt, atmo, ice, cloud, rain, graupel, &
            dep_rate_ice)
       IF (ischeck) CALL check(ik_slice, 'ice_riming',cloud,rain,ice,snow,graupel,hail)

       ! riming of snow with cloud droplets and rain drops, and conversion to graupel
       CALL snow_riming(ik_slice, dt, atmo, &
            snow, cloud, rain, ice, graupel, dep_rate_snow)
       IF (ischeck) CALL check(ik_slice, 'snow_riming',cloud,rain,ice,snow,graupel,hail)

       ! more riming
       IF (ice_typ > 1) THEN
          CALL hail_cloud_riming(ik_slice, dt, atmo, hail, cloud, rain, ice)
          CALL hail_rain_riming(ik_slice, dt, atmo, hail, rain, ice)
      IF (ischeck) CALL check(ik_slice, 'hail riming',cloud,rain,ice,snow,graupel,hail)
       END IF
       CALL graupel_cloud_riming(ik_slice, dt, atmo, graupel, cloud, rain, ice)
       CALL graupel_rain_riming(ik_slice, dt, atmo, graupel, rain, ice)
      IF (ischeck) CALL check(ik_slice, 'graupel riming',cloud,rain,ice,snow,graupel,hail)

       ! freezing of rain and conversion to ice/graupel/hail
       CALL rain_freeze_gamlook(ik_slice, dt, atmo,rain,ice,snow,graupel,hail)
       IF (ischeck) CALL check(ik_slice, 'rain_freeze_gamlook',cloud,rain,ice,snow,graupel,hail)

       ! melting
       CALL ice_melting(ik_slice, atmo, ice, cloud, rain)
       CALL snow_melting(ik_slice, dt, atmo,snow,rain)
       CALL graupel_melting(ik_slice, dt, atmo, graupel, rain)
       IF (ice_typ > 1) CALL hail_melting(ik_slice, dt, atmo, hail, rain)
       IF (ischeck) CALL check(ik_slice, 'melting',cloud,rain,ice,snow,graupel,hail)

       ! evaporation from melting ice particles
       CALL snow_evaporation(ik_slice, dt, atmo, snow)
       CALL graupel_evaporation(ik_slice, dt, atmo, graupel)
       IF (ice_typ > 1) CALL hail_evaporation(ik_slice, dt, atmo, hail)
       IF (ischeck) CALL check(ik_slice, 'evaporation of ice',cloud,rain,ice,snow,graupel,hail)

       end if
    ENDIF

    ! warm rain processes
    ! (using something other than SB is somewhat inconsistent and not recommended)
    IF (auto_typ == 1) THEN
       CALL autoconversionKB(ik_slice, dt, cloud, rain)   ! Beheng (1994)
       CALL accretionKB(ik_slice, dt, cloud, rain)
       CALL rain_selfcollectionSB(ik_slice, dt, rain)
    ELSE IF (auto_typ == 2) THEN
       ! Khairoutdinov and Kogan (2000)
       ! (KK2000 originally assume a 25 micron size threshold)
       CALL autoconversionKK(ik_slice, dt, cloud, rain)
       CALL accretionKK(ik_slice, dt, cloud, rain)
       CALL rain_selfcollectionSB(ik_slice, dt, rain)
    ELSE IF (auto_typ == 3) THEN
       CALL autoconversionSB(ik_slice, dt, cloud, rain)  ! Seifert and Beheng (2001)
       CALL accretionSB(ik_slice, dt, cloud, rain)
       CALL rain_selfcollectionSB(ik_slice, dt, rain)
    ENDIF
    IF (ischeck) CALL check(ik_slice, 'warm rain',cloud,rain,ice,snow,graupel,hail)

    ! evaporation of rain following Seifert (2008)
    CALL rain_evaporation(ik_slice, dt, atmo, cloud, rain)

    ! size limits for all hydrometeors
    IF (nuc_c_typ > 0) THEN
       DO k=kstart,kend
        DO i=istart,iend
          cloud%n(i,k) = MIN(cloud%n(i,k), cloud%q(i,k)/cloud%x_min)
          cloud%n(i,k) = MAX(cloud%n(i,k), cloud%q(i,k)/cloud%x_max)
          ! Hard upper limit for cloud number conc.
          cloud%n(i,k) = MIN(cloud%n(i,k), 5000d6)
        END DO
       END DO
    END IF
    DO k=kstart,kend
       DO i=istart,iend
          rain%n(i,k) = MIN(rain%n(i,k), rain%q(i,k)/rain%x_min)
          rain%n(i,k) = MAX(rain%n(i,k), rain%q(i,k)/rain%x_max)
       END DO
    END DO
    IF (ice_typ > 0) THEN
       DO k=kstart,kend
          DO i=istart,iend
             ice%n(i,k) = MIN(ice%n(i,k), ice%q(i,k)/ice%x_min)
             ice%n(i,k) = MAX(ice%n(i,k), ice%q(i,k)/ice%x_max)
             snow%n(i,k) = MIN(snow%n(i,k), snow%q(i,k)/snow%x_min)
             snow%n(i,k) = MAX(snow%n(i,k), snow%q(i,k)/snow%x_max)
             graupel%n(i,k) = MIN(graupel%n(i,k), graupel%q(i,k)/graupel%x_min)
             graupel%n(i,k) = MAX(graupel%n(i,k), graupel%q(i,k)/graupel%x_max)
          END DO
       END DO
    END IF
    IF (ice_typ > 1) THEN
       DO k=kstart,kend
          DO i=istart,iend
             hail%n(i,k) = MIN(hail%n(i,k), hail%q(i,k)/hail%x_min)
             hail%n(i,k) = MAX(hail%n(i,k), hail%q(i,k)/hail%x_max)
          END DO
       END DO
    END IF

    IF (ischeck) CALL check(ik_slice, 'clouds_twomoment end',cloud,rain,ice,snow,graupel,hail)
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
  subroutine sedi_vel_rain(this,thisCoeffs,q,x,vn,vq,its,ite,qc)
    TYPE(particle), intent(in)      :: this
    TYPE(particle_rain), intent(in) :: thisCoeffs
    integer,  intent(in)  :: its,ite
    real(wp), intent(in)  :: q(:),x(:)
    real(wp), intent(in), optional  :: qc(:)
    real(wp), intent(inout) :: vn(:), vq(:)

    integer  :: i
    real(wp) :: D_m,mue,D_p

    do i=its,ite
       if (q(i).gt.q_crit) then
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
       else
          vn(i) = 0.0_wp
          vq(i) = 0.0_wp
       end if
    end do
  end subroutine sedi_vel_rain

  ! bulk sedimentation velocities
  subroutine sedi_vel_sphere(this,thisCoeffs,q,x,vn,vq,its,ite)
    TYPE(particle), intent(in)        :: this
    TYPE(particle_sphere), intent(in) :: thisCoeffs
    integer,  intent(in)  :: its,ite
    real(wp), intent(in)  :: q(:),x(:)
    real(wp), intent(out) :: vn(:), vq(:)

    integer  :: i
    real(wp) :: lam,v_n,v_q

    do i=its,ite
       if (q(i).gt.q_crit) then
          lam = exp(this%b_vel* log(thisCoeffs%coeff_lambda*x(i)))
          v_n = thisCoeffs%coeff_alfa_n * lam
          v_q = thisCoeffs%coeff_alfa_q * lam
          v_n = MAX(v_n,this%vsedi_min)
          v_q = MAX(v_q,this%vsedi_min)
          v_n = MIN(v_n,this%vsedi_max)
          v_q = MIN(v_q,this%vsedi_max)
          vn(i) = v_n
          vq(i) = v_q
       else
          vn(i) = 0.0_wp
          vq(i) = 0.0_wp
       end if
    end do
  end subroutine sedi_vel_sphere

  ! initialize coefficients for bulk sedimentation velocity
  subroutine init_sedi_vel(this,thisCoeffs)
    TYPE(particle), INTENT(in) :: this
    TYPE(particle_sphere), INTENT(out) :: thisCoeffs

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

  SUBROUTINE autoconversionSB(ik_slice, dt, cloud,rain)
    !*******************************************************************************
    ! Autoconversion of Seifert and Beheng (2001, Atmos. Res.)                     *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(particle), INTENT(inout) :: cloud, rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER          :: i,k
    INTEGER,  SAVE   :: firstcall
    REAL(wp), SAVE   :: k_au,k_sc
    REAL(wp)         :: q_c, q_r, n_c, x_c, nu, mu, tau, phi, x_s_i, au, sc

    REAL(wp), PARAMETER :: k_c  = 9.44e+9_wp   !..Long-Kernel
    REAL(wp), PARAMETER :: k_1  = 6.00e+2_wp   !..Parameter for Phi
    REAL(wp), PARAMETER :: k_2  = 0.68e+0_wp   !..Parameter fof Phi
    REAL(wp), PARAMETER :: eps  = 1.00e-25_wp
!$omp threadprivate (firstcall)
!$omp threadprivate (k_au)
!$omp threadprivate (k_sc)

    IF (isdebug) CALL message(routine, "autoconversionSB")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    x_s_i = 1.0_wp / cloud%x_max

    IF (firstcall.NE.1) THEN
       nu = cloud%nu
       mu = cloud%mu
       IF (mu == 1.0) THEN
          !.. see SB2001
          k_au  = k_c * x_s_i * (1.0_wp / 20.0_wp) * (nu+2.0)*(nu+4.0)/(nu+1.0)**2
          k_sc  = k_c * (nu+2.0)/(nu+1.0)
       ELSE
          !.. see Eq. (3.44) of Seifert (2002)
         k_au = k_c * x_s_i * (1.0_wp / 20.0_wp)                            &
               & * ( 2.0_wp * gfct((nu+4.0)/mu)**1                          &
               &            * gfct((nu+2.0)/mu)**1 * gfct((nu+1.0)/mu)**2   &
               &   - 1.0_wp * gfct((nu+3.0)/mu)**2 * gfct((nu+1.0)/mu)**2 ) &
               &   / gfct((nu+2.0)/mu)**4
          k_sc = k_c * moment_gamma(cloud,2)
       ENDIF
    END IF

    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          IF (q_c > q_crit) THEN

            n_c = cloud%n(i,k)
            q_r = rain%q(i,k)
            x_c = cloud%meanmass(q_c,n_c)

            au  = k_au * q_c**2 * x_c**2 * dt * cloud%rho_v(i,k)
            tau = MIN(MAX(1.0-q_c/(q_c+q_r+eps),eps),0.9_wp)
            phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
            au  = au * (1.0 + phi/(1.0 - tau)**2)

            au  = MAX(MIN(q_c,au),0.0_wp)

            sc  = k_sc * q_c**2 * dt * cloud%rho_v(i,k)

            rain%n(i,k)  = rain%n(i,k)  + au * x_s_i
            rain%q(i,k)  = rain%q(i,k)  + au
            cloud%n(i,k) = cloud%n(i,k) - MIN(n_c,sc)
            cloud%q(i,k) = cloud%q(i,k) - au

          ENDIF
       END DO
    END DO

  END SUBROUTINE autoconversionSB

  SUBROUTINE accretionSB(ik_slice, dt, cloud,rain)
    !*******************************************************************************
    ! Accretion of Seifert and Beheng (2001, Atmos. Res.)                          *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(particle), INTENT(inout)  :: cloud, rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER     :: i, k
    REAL(wp)    :: ac
    REAL(wp)    :: q_c, q_r, tau, phi, n_c, x_c

    REAL(wp), PARAMETER :: k_r = 5.78_wp       ! kernel
    REAL(wp), PARAMETER :: k_1 = 5.00e-04_wp   ! Phi function
    REAL(wp), PARAMETER :: eps = 1.00e-25_wp

    IF (isdebug) CALL message(routine, "accretionSB")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          n_c = cloud%n(i,k)
          q_c = cloud%q(i,k)
          q_r = rain%q(i,k)
          IF (q_c > 0.0_wp.AND.q_r > 0.0_wp) THEN

             !..accretion rate of SB2001
             tau = MIN(MAX(1.0-q_c/(q_c+q_r+eps),eps),1.0_wp)
             phi = (tau/(tau+k_1))**4
             ac  = k_r *  q_c * q_r * phi * dt

             ac = MIN(q_c,ac)

             x_c = cloud%meanmass(q_c,n_c)

             rain%q(i,k)  = rain%q(i,k)  + ac
             cloud%q(i,k) = cloud%q(i,k) - ac
             cloud%n(i,k) = cloud%n(i,k) - MIN(n_c,ac/x_c)
          ENDIF
       END DO
    END DO

  END SUBROUTINE accretionSB

  SUBROUTINE rain_selfcollectionSB(ik_slice, dt, rain)
    !*******************************************************************************
    ! Selfcollection of Seifert and Beheng (2001, Atmos. Res.)                     *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(particle) , INTENT(inout) :: rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER        :: i, k
    REAL(wp)       :: sc, br
    REAL(wp)       :: q_r, n_r, x_r, d_r

    !..Parameters based on Seifert (2008, JAS)
    REAL(wp), PARAMETER :: D_br = 1.10e-3_wp
    REAL(wp), PARAMETER :: k_rr = 4.33e+0_wp
    REAL(wp), PARAMETER :: k_br = 1.00e+3_wp

    IF (isdebug) CALL message(routine, "rain_selfcollectionSB")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          n_r = rain%n(i,k)
          q_r = rain%q(i,k)
          IF (q_r > 0.0_wp) THEN
            x_r = rain%meanmass(q_r,n_r)
            D_r = rain%diameter(x_r)

            !..Selfcollection as in SB2001
            sc = k_rr *  n_r * q_r * rain%rho_v(i,k) * dt

            !..Breakup as in Seifert (2008, JAS), Eq. (A13)
            br = 0.0_wp
            IF (D_r.GT.0.30e-3_wp) THEN
               br = (k_br * (D_r - D_br) + 1.0)  * sc
            ENDIF

            rain%n(i,k) = n_r - MIN(n_r,sc-br)
         ENDIF
      END DO
   END DO

  END SUBROUTINE rain_selfcollectionSB

  SUBROUTINE autoconversionKB(ik_slice, dt, cloud, rain)
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(particle), INTENT(inout) :: cloud, rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER          :: i,k
    REAL(wp)         :: q_c, x_c, nu_c, n_c, k_a, x_s_i, au

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    nu_c = 9.59
    x_s_i = 1.0_wp / cloud%x_max
    k_a  = 6.0d+25 * nu_c**(-1.7)

    !..Parameterization of Beheng (1994)
    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          n_c = cloud%n(i,k)
          IF (q_c > q_crit) THEN

             x_c = cloud%meanmass(q_c,n_c)

             !..Berechnung der Autokonversionsrate nach Beheng (1994)
             au = k_a * (x_c*1e3)**(3.3) * (q_c*1e-3)**(1.4) * dt * 1e3
             au = MIN(q_c,au)

             rain%n(i,k)  = rain%n(i,k)  + au * x_s_i
             rain%q(i,k)  = rain%q(i,k)  + au
             cloud%n(i,k) = cloud%n(i,k) - au * x_s_i * 2.0
             cloud%q(i,k) = cloud%q(i,k) - au

          ENDIF
       END DO
    END DO

  END SUBROUTINE autoconversionKB

  SUBROUTINE accretionKB(ik_slice, dt, cloud, rain)
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(particle), INTENT(inout)   :: cloud, rain

    !..Parameter of Beheng (1994)
    REAL(wp), PARAMETER :: k_r = 6.00d+00
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER      :: i,k
    REAL(wp)     :: ac
    REAL(wp)     :: q_c, q_r

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          q_r = rain%q(i,k)

          IF (q_c > q_crit .and. q_r > q_crit) THEN

             ac = k_r *  q_c * q_r * dt
             ac = MIN(q_c,ac)

             rain%q(i,k)  = rain%q(i,k)  + ac
             cloud%q(i,k) = cloud%q(i,k) - ac

          ENDIF
       END DO
    END DO

  END SUBROUTINE accretionKB

  SUBROUTINE autoconversionKK(ik_slice, dt, cloud,rain)

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(particle), INTENT(inout) :: cloud, rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_c, x_s_i, au, n_c
    REAL(wp), PARAMETER :: k_a  = 1350.0_wp

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    x_s_i  = 1.0_wp / cloud%x_max

    !..Parameterization of Khairoutdinov and Kogan (2000), MWR 128, 229-243
    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          IF (q_c > q_crit) THEN

             n_c = cloud%n(i,k) * 1e6  ! in 1/cm3

             !..autoconversion rate of KK2000
             au = k_a * q_c**(2.47) * n_c**(-1.79) * dt
             au = MIN(q_c,au)

             rain%n(i,k)  = rain%n(i,k)  + au * x_s_i
             rain%q(i,k)  = rain%q(i,k)  + au
             cloud%n(i,k) = cloud%n(i,k) - au * x_s_i * 2.0
             cloud%q(i,k) = cloud%q(i,k) - au

          ENDIF
       END DO
    END DO

  END SUBROUTINE autoconversionKK

  SUBROUTINE accretionKK(ik_slice, dt, cloud, rain)
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(particle), INTENT(inout) :: cloud, rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: ac, q_c, q_r
    REAL(wp), PARAMETER :: k_a = 6.70d+01

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    !..Parameterization of Khairoutdinov and Kogan (2000), MWR 128, 229-243
    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          q_r = rain%q(i,k)

          IF (q_c > q_crit .AND. q_r > q_crit) THEN

             ac = k_a *  (q_c * q_r)**1.15 * dt * 1e3
             ac = MIN(q_c,ac)

             rain%q(i,k)  = rain%q(i,k)  + ac
             cloud%q(i,k) = cloud%q(i,k) - ac
          ENDIF
       END DO
    END DO

  END SUBROUTINE accretionKK

  SUBROUTINE rain_evaporation(ik_slice, dt, atmo,cloud,rain)
    !*******************************************************************************
    ! Evaporation of rain based on Seifert (2008, J. Atmos. Sci.)                  *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle)      :: cloud,rain
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    INTEGER, SAVE       :: firstcall
    REAL(wp)            :: T_a,p_a,e_sw,s_sw,g_d,eva_q,eva_n
    REAL(wp)            :: q_r,n_r,x_r,e_d,f_v
    REAL(wp)            :: mue,d_m,gamma_eva,lam,d_vtp,gfak
    REAL(wp)            :: aa,bb,cc,mm
    REAL(wp), SAVE      :: a_q,b_q
    REAL(wp), PARAMETER :: D_br = 1.1e-3_wp
!$omp threadprivate (firstcall)
!$omp threadprivate (a_q)
!$omp threadprivate (b_q)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

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
      CALL message(routine, "rain_evaporation")
    END IF

    DO k = kstart,kend
       DO i = istart,iend

          q_r  = rain%q(i,k)
          n_r  = rain%n(i,k)
          p_a  = atmo%p(i,k)
          T_a  = atmo%T(i,k)

          e_d  = atmo%qv(i,k) * R_d * T_a
          e_sw = e_ws(T_a)
          s_sw = e_d / e_sw - 1.0

          IF (s_sw < 0.0 .AND. q_r > 0.0_wp .AND. cloud%q(i,k) < q_crit) THEN

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
                  &            * sqrt(aa/nu_l * rain%rho_v(i,k) / lam)                  &
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

             atmo%qv(i,k)     = atmo%qv(i,k)     + eva_q
             rain%q(i,k) = rain%q(i,k) - eva_q
             rain%n(i,k) = rain%n(i,k) - eva_n
          END IF
       END DO
    END DO

  END SUBROUTINE rain_evaporation

  SUBROUTINE graupel_evaporation(ik_slice, dt, atmo, graupel)
    !*******************************************************************************
    ! Evaporation from melting graupel, see SB2006                                 *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle)      :: graupel
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a,e_sw,s_sw,g_d,eva
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g,f_v,e_d

    IF (isdebug) CALL message(routine, "graupel_evaporation")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_g = graupel%q(i,k)
          T_a = atmo%T(i,k)

          IF (q_g > 0.0_wp .AND. T_a > T_3) THEN

            e_d  = atmo%qv(i,k) * R_d * T_a
            e_sw = e_ws(T_a)
            s_sw = e_d / e_sw - 1.0

            !..Eq. (37) of SB2006, note that 4*pi is correct because c_g is used below
            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_g = graupel%n(i,k)
            x_g = graupel%meanmass(q_g,n_g)
            d_g = graupel%diameter(x_g)
            v_g = graupel%velocity(x_g) * graupel%rho_v(i,k)

            !..note that a_f includes more than just ventilation, do never ever set f_v=1
            f_v = ge_params%a_f + ge_params%b_f * sqrt(v_g*d_g)

            eva = g_d * n_g * ge_params%c * d_g * f_v * s_sw * dt

            eva = MAX(-eva,0.0_wp)
            eva = MIN(eva,q_g)

            atmo%qv(i,k)   = atmo%qv(i,k)   + eva
            graupel%q(i,k) = graupel%q(i,k) - eva

          END IF
       END DO
    END DO

  END SUBROUTINE graupel_evaporation

  SUBROUTINE hail_evaporation(ik_slice, dt, atmo, hail)
    !*******************************************************************************
    ! Evaporation of melting graupel, see SB2006                                   *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: hail
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a,e_sw,s_sw,g_d,eva
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h,f_v,e_d

    IF (isdebug) CALL message(routine, "hail_evaporation")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_h = hail%q(i,k)
          T_a = atmo%T(i,k)

          IF (q_h > 0.0_wp .AND. T_a > T_3) THEN

            e_d  = atmo%qv(i,k) * R_d * T_a
            e_sw = e_ws(T_a)
            s_sw = e_d / e_sw - 1.0

            !..Eq. (37) of SB2006, note that 4*pi is correct because c_h is used below
            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_h = hail%n(i,k)
            x_h = hail%meanmass(q_h,n_h)
            d_h = hail%diameter(x_h)
            v_h = hail%velocity(x_h) * hail%rho_v(i,k)

            f_v  = he_params%a_f + he_params%b_f * sqrt(v_h*D_h)

            eva = g_d * n_h * he_params%c * d_h * f_v * s_sw * dt

            eva = MAX(-eva,0.0_wp)
            eva = MIN(eva,q_h)

            atmo%qv(i,k) = atmo%qv(i,k) + eva
            hail%q(i,k)  = hail%q(i,k)  - eva

          END IF
       END DO
    END DO
  END SUBROUTINE hail_evaporation

  SUBROUTINE snow_evaporation(ik_slice, dt, atmo, snow)
    !*******************************************************************************
    ! Evaporation of melting snow, see SB2006                                      *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: snow

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a,e_sw,s_sw,g_d,eva
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s,f_v,e_d

    IF (isdebug) CALL message(routine, "snow_evaporation")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_s = snow%q(i,k)
          T_a = atmo%T(i,k)

          IF (q_s > 0.0_wp .AND. T_a > T_3) THEN

            e_d  = atmo%qv(i,k) * R_d * T_a
            e_sw = e_ws(T_a)
            s_sw = e_d / e_sw - 1.0

            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_s = snow%n(i,k)

            x_s = snow%meanmass(q_s,n_s)
            D_s = snow%diameter(x_s)
            v_s = snow%velocity(x_s) * snow%rho_v(i,k)

            f_v  = se_params%a_f + se_params%b_f * sqrt(v_s*D_s)

            eva = g_d * n_s * se_params%c * d_s * f_v * s_sw * dt

            eva = MAX(-eva,0.0_wp)
            eva = MIN(eva,q_s)

            atmo%qv(i,k) = atmo%qv(i,k) + eva
            snow%q(i,k)  = snow%q(i,k)  - eva

          END IF
       END DO
    END DO
  END SUBROUTINE snow_evaporation

  SUBROUTINE cloud_freeze(ik_slice, dt, atmo, cloud, ice)
    !*******************************************************************************
    ! This is only the homogeneous freezing of liquid water droplets.              *
    ! Immersion freezing and homogeneous freezing of liquid aerosols are           *
    ! treated in the subroutine ice_nucleation_homhet()                            *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout)  :: cloud, ice
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER            :: i, k
    REAL(wp)           :: fr_q,fr_n,T_a,q_c,x_c,n_c,j_hom,T_c
    REAL(wp), SAVE     :: coeff_z
    INTEGER,  SAVE     :: firstcall
!$omp threadprivate (firstcall)
!$omp threadprivate (coeff_z)

    IF (firstcall.NE.1) THEN
       firstcall = 1
       coeff_z   = moment_gamma(cloud,2)
    END IF

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)
          IF (T_a < T_3) THEN

             T_c = T_a - T_3
             q_c = cloud%q(i,k)
             n_c = cloud%n(i,k)
             IF (q_c > 0.0_wp .and. T_c < -30.0_wp) THEN
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

                cloud%q(i,k) = cloud%q(i,k) - fr_q
                cloud%n(i,k) = cloud%n(i,k) - fr_n

                fr_n  = MAX(fr_n,fr_q/cloud%x_max)

                !..special treatment for constant drop number
                IF (nuc_c_typ .EQ. 0) THEN
                   ! ... force upper bound in cloud_freeze'
                   fr_n = MAX(MIN(fr_n,qnc_const-ice%n(i,k)),0.0_wp)
                ENDIF

                ice%q(i,k)   = ice%q(i,k) + fr_q
                ice%n(i,k)   = ice%n(i,k) + fr_n
             ENDIF
          END IF
       END DO
    END DO
  END SUBROUTINE cloud_freeze

  SUBROUTINE ice_nucleation_homhet(ik_slice, use_prog_in, &
       atmo, cloud, ice, snow, n_inact, n_inpot)
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
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    LOGICAL, INTENT(in) :: use_prog_in

    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud, ice, snow
    REAL(wp), DIMENSION(:,:) :: n_inact
    REAL(wp), DIMENSION(:,:), OPTIONAL :: n_inpot

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER              :: i,k,nuc_typ
    REAL(wp)             :: nuc_n, nuc_q
    REAL(wp)             :: T_a, p_a, ssi
    REAL(wp)             :: q_i,n_i,x_i,r_i
    REAL(wp)             :: ndiag, ndiag_dust, ndiag_all

    ! switch for version of Phillips et al. scheme
    ! (but make sure you have the correct INCLUDE file)
    INTEGER, PARAMETER :: iphillips = 2010

    ! switch for Hande et al. ice nucleation, if .true. this turns off Phillips scheme
    LOGICAL              :: use_hdcp2_het = .false.

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

    ! variables for Hande et al. nucleation parameterization for HDCP2 simulations
    REAL(wp)  :: alf_dep, bet_dep, nin_dep
    REAL(wp)  :: alf_imm, bet_imm, nin_imm

    ! parameters for deposition formula, Eq (2) of Hande et al.
    REAL(wp), PARAMETER :: a_dep =  2.7626_wp
    REAL(wp), PARAMETER :: b_dep =  6.2100_wp  ! 0.0621*100, because we use ssi instead of RHi
    REAL(wp), PARAMETER :: c_dep = -1.3107_wp
    REAL(wp), PARAMETER :: d_dep =  2.6789_wp

    REAL(wp), PARAMETER :: eps = 1.e-20_wp

    ! variables for interpolation in look-up table (real is good enough here)
    REAL      :: xt,xs,ssr
    INTEGER   :: ss,tt

    LOGICAL   :: use_homnuc, lwrite_n_inpot

    REAL(wp), DIMENSION(3) :: infrac

    LOGICAL :: ndiag_mask(ik_slice(1):ik_slice(2), ik_slice(3):ik_slice(4))
    REAL(wp) :: &
         nuc_n_a(ik_slice(1):ik_slice(2), ik_slice(3):ik_slice(4))



    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

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

    ! switch for Hande et al. ice nucleation, if .true. this turns off Phillips scheme
    IF (nuc_typ.LE.5) use_hdcp2_het = .TRUE.

    ! Heterogeneous nucleation using Hande et al. scheme
    IF (use_hdcp2_het) THEN
      IF (nuc_typ.EQ.1) THEN
        ! Spring of Table 1
        nin_imm = 1.5684e5
        alf_imm = 0.2466
        bet_imm = 1.2293
        nin_dep = 1.7836e5
        alf_dep = 0.0075
        bet_dep = 2.0341
      ELSEIF (nuc_typ.EQ.2) THEN
        ! Summer
        nin_imm = 2.9694e4
        alf_imm = 0.2813
        bet_imm = 1.1778
        nin_dep = 2.6543e4
        alf_dep = 0.0020
        bet_dep = 2.5128
      ELSEIF (nuc_typ.EQ.3) THEN
        ! Autumn
        nin_imm = 4.9920e4
        alf_imm = 0.2622
        bet_imm = 1.2044
        nin_dep = 7.7167e4
        alf_dep = 0.0406
        bet_dep = 1.4705
      ELSEIF (nuc_typ.EQ.4) THEN
        ! Winter
        nin_imm = 1.0259e5
        alf_imm = 0.2073
        bet_imm = 1.2873
        nin_dep = 1.1663e4
        alf_dep = 0.0194
        bet_dep = 1.6943
      ELSEIF (nuc_typ.EQ.5) THEN
        ! Spring with 95th percentile scaling factor
        nin_imm = 1.5684e5 * 17.82
        alf_imm = 0.2466
        bet_imm = 1.2293
        nin_dep = 1.7836e5 * 5.87
        alf_dep = 0.0075
        bet_dep = 2.0341
      ELSE
        CALL finish(TRIM(routine),&
             & 'Error in two_moment_mcrph: Invalid value nuc_typ in case of use_hdcp2_het=.true.')
      END IF

    ELSE
      ! Heterogeneous nucleation using Phillips et al. scheme
      IF (iphillips == 2010) THEN
        ! possible pre-defined choices
        IF (nuc_typ.EQ.6) THEN  ! with no organics and rather high soot, coming close to Meyers formula at -20 C
          na_dust  = 160.e4_wp    ! initial number density of dust [1/m3]
          na_soot  =  30.e6_wp    ! initial number density of soot [1/m3]
          na_orga  =   0.e0_wp    ! initial number density of organics [1/m3]
        ELSEIF (nuc_typ.EQ.7) THEN     ! with some organics and rather high soot,
          na_dust  = 160.e4_wp    !          coming close to Meyers formula at -20 C
          na_soot  =  25.e6_wp
          na_orga  =  30.e6_wp
        ELSEIF (nuc_typ.EQ.8) THEN     ! no organics, no soot, coming close to DeMott et al. 2010 at -20 C
          na_dust  =  70.e4_wp    ! i.e. roughly one order in magnitude lower than Meyers
          na_soot  =   0.e6_wp
          na_orga  =   0.e6_wp
        ELSE
          CALL finish(TRIM(routine),&
               & 'Error in two_moment_mcrph: Invalid value nuc_typ in case of use_hdcp2_het=.false.')
        END IF
      END IF
    END IF

    DO k = kstart,kend
      DO i = istart,iend

        p_a  = atmo%p(i,k)
        T_a  = atmo%T(i,k)
        e_si = e_es(T_a)
        ssi  = atmo%qv(i,k) * R_d * T_a / e_si

        IF (T_a < T_nuc .AND. T_a > 180.0_wp .AND. ssi > 1.0_wp  &
             & .AND. ( ice%n(i,k)+snow%n(i,k) < ni_het_max ) ) THEN

          xt = (274.- REAL(atmo%T(i,k)))  / ttstep
          xt = MIN(xt,REAL(ttmax-1))
          tt = INT(xt)

           IF (cloud%q(i,k) > 0.0_wp) THEN

              ! immersion freezing at water saturation

              IF (use_hdcp2_het) THEN
                 ! Hande et al. scheme, Eq. (1)
                 if (T_a.lt.261.15_wp) then
                    ndiag = nin_imm * exp( - alf_imm * exp(bet_imm*log(MAX(eps,T_a - 237.15_wp))) )
                 else
                    ndiag = 0.0_wp
                 end if
              ELSE
                 ! Phillips scheme
                 ! immersion freezing at water saturation
                 infrac(1) = (AINT(xt)+1.0-xt) * afrac_dust(tt,99) &
                      &        + (xt-AINT(xt)) * afrac_dust(tt+1,99)
                 infrac(2) = (AINT(xt)+1.0-xt) * afrac_soot(tt,99) &
                      &        + (xt-AINT(xt)) * afrac_soot(tt+1,99)
                 infrac(3) = (AINT(xt)+1.0-xt) * afrac_orga(tt,99) &
                      &        + (xt-AINT(xt)) * afrac_orga(tt+1,99)
              END IF
            ELSE
              ! deposition nucleation below water saturation

              IF (use_hdcp2_het) THEN
                 ! Hande et al. scheme, Eq. (3) with (2) and (1)
                 if (T_a.lt.253.0_wp) then
                    ndiag = nin_dep * exp( - alf_dep * exp(bet_dep*log(MAX(eps,T_a - 220.0_wp))) )
                    ndiag = ndiag * (a_dep * atan(b_dep*(ssi-1.0_wp)+c_dep) + d_dep)
                 else
                    ndiag = 0.0_wp
                 end if
              ELSE
                ! deposition nucleation below water saturation
                ! calculate indices used for 2D look-up tables
                xs = 100. * REAL(ssi-1.0_wp) / ssstep
                xs = MIN(xs,REAL(ssmax-1))
                ss = MAX(1,INT(xs))
                ssr = MAX(1.0, AINT(xs))
                ! bi-linear interpolation in look-up tables
                infrac(1) =   (AINT(xt) + 1.0 - xt) * (ssr + 1.0 - xs) &
                     &        * afrac_dust(tt, ss) &
                     &      + (xt - AINT(xt)) * (ssr + 1.0 - xs) &
                     &        * afrac_dust(tt+1, ss) &
                     &      + (AINT(xt) + 1.0 - xt) * (xs - ssr) &
                     &        * afrac_dust(tt, ss+1) &
                     &      + (xt - AINT(xt)) * (xs - ssr) &
                     &        * afrac_dust(tt+1, ss+1)
                infrac(2) =   (AINT(xt) + 1.0 - xt) * (ssr + 1.0 - xs) &
                     &        * afrac_soot(tt, ss) &
                     &      + (xt - AINT(xt)) * (ssr + 1.0 - xs) &
                     &        * afrac_soot(tt+1, ss  ) &
                     &      + (AINT(xt) + 1.0 - xt) * (xs - ssr) &
                     &        * afrac_soot(tt, ss + 1) &
                     &      + (xt - AINT(xt)) * (xs - ssr) &
                     &        * afrac_soot(tt+1, ss+1)
                infrac(3) = (AINT(xt) + 1.0 - xt) * (ssr + 1.0 - xs) &
                     &        * afrac_orga(tt,ss) &
                     &      + (xt - AINT(xt)) * (ssr + 1.0 - xs) &
                     &        * afrac_orga(tt+1, ss) &
                     &      + (AINT(xt) + 1.0 - xt) * (xs - ssr) &
                     &        * afrac_orga(tt, ss+1) &
                     &      + (xt - AINT(xt)) * (xs - ssr) &
                     &        * afrac_orga(tt+1, ss+1)
              END IF
            END IF

          IF (.NOT. use_hdcp2_het) THEN
            ! only for Phillips scheme we have to sum up the three modes
            ndiag = na_soot * infrac(2) + na_orga * infrac(3)
            IF (use_prog_in) THEN
              ! n_inpot replaces na_dust, na_soot and na_orga are assumed to constant
              ndiag  = n_inpot(i,k) * infrac(1) + ndiag
              ndiag_dust = n_inpot(i,k) * infrac(1)
              ndiag_all = ndiag
            ELSE
              ! all aerosol species are diagnostic
              ndiag = na_dust * infrac(1) + ndiag
            END IF
            ndiag = MIN(ndiag,ni_het_max)
          END IF

          nuc_n = MAX(ndiag - n_inact(i,k),0.0_wp)
          nuc_q = MIN(nuc_n * ice%x_min, atmo%qv(i,k))
          nuc_n = nuc_q / ice%x_min

          ice%n(i,k)   = ice%n(i,k)   + nuc_n
          ice%q(i,k)   = ice%q(i,k)   + nuc_q
          atmo%qv(i,k) = atmo%qv(i,k) - nuc_q
          n_inact(i,k) = n_inact(i,k) + nuc_n

          lwrite_n_inpot = use_prog_in .AND. ndiag .GT. 1d-12
          ndiag_mask(i, k) = lwrite_n_inpot

          IF (lwrite_n_inpot .AND. .NOT. use_hdcp2_het) THEN
            nuc_n = nuc_n * ndiag_dust / ndiag_all
          END IF
          nuc_n_a(i, k) = nuc_n
        ELSE
          nuc_n_a(i, k) = 0.0_wp
          ndiag_mask(i, k) = .FALSE.
        ENDIF

      END DO
    END DO

    IF (use_prog_in) THEN
      DO k = kstart, kend
        DO i = istart, iend
          n_inpot(i,k) = MERGE(MAX(n_inpot(i,k) - nuc_n_a(i, k), 0.0_wp), &
               &               n_inpot(i, k), &
               &               ndiag_mask(i, k))
        END DO
      END DO
    END IF

    ! Homogeneous nucleation using KHL06 approach
    IF (use_homnuc) THEN
      DO k = kstart,kend
        DO i = istart,iend
          p_a  = atmo%p(i,k)
          T_a  = atmo%T(i,k)
          e_si = e_es(T_a)
          ssi  = atmo%qv(i,k) * R_d * T_a / e_si

          ! critical supersaturation for homogeneous nucleation
          scr  = 2.349 - T_a * (1.0_wp/ 259.00_wp)

          IF (ssi > scr .AND. T_a < 235.0 .AND. ice%n(i,k) < ni_hom_max ) THEN

            n_i = ice%n(i,k)
            q_i = ice%q(i,k)
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

            IF (atmo%w(i,k) > w_pre) THEN   ! homogenous nucleation event

              ! timescales of freezing event (see KL02, RM05, KHL06)
              cool    = grav / cp * atmo%w(i,k)
              ctau    = T_a * ( 0.004*T_a - 2. ) + 304.4
              tau     = 1.0 / (ctau * cool)                       ! freezing timescale, eq. (5)
              delta   = (bcoeff(2) * r_0)                         ! dimless aerosol radius, eq.(4)
              tau_g   = (bcoeff(1) / r_0) / (1 + delta)           ! timescale for initial growth, eq.(4)
              phi     = acoeff(1)*ssi / ( acoeff(2) + acoeff(3)*ssi) * (atmo%w(i,k) - w_pre)

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
              nuc_q = MIN(nuc_n * mi_hom, atmo%qv(i,k))

              ice%n(i,k) = ice%n(i,k) + nuc_n
              ice%q(i,k) = ice%q(i,k) + nuc_q
              atmo%qv(i,k)  = atmo%qv(i,k)  - nuc_q

            END IF
          END IF
        ENDDO
      ENDDO
    END IF

  END SUBROUTINE ice_nucleation_homhet

  SUBROUTINE vapor_dep_relaxation(ik_slice, dt_local, &
       &               atmo, ice, snow, graupel, hail, dep_rate_ice, dep_rate_snow)
    !*******************************************************************************
    ! Deposition and sublimation                                                   *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    TYPE(atmosphere)         :: atmo
    TYPE(particle)           :: ice, snow, graupel, hail
    REAL(wp), INTENT(IN)     :: dt_local
    REAL(wp), INTENT(INOUT), DIMENSION(:,:) :: dep_rate_ice, dep_rate_snow

    REAL(wp), DIMENSION(size(dep_rate_ice,1),size(dep_rate_ice,2)) :: &
                                    & s_si,g_i,dep_ice,dep_snow,dep_graupel,dep_hail

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: D_vtp
    REAL(wp)            :: zdt,qvsidiff,Xi_i,Xfac
    REAL(wp)            :: tau_i_i,tau_s_i,tau_g_i,tau_h_i
    REAL(wp), PARAMETER :: eps  = 1.d-20
    REAL(wp)            :: T_a
    REAL(wp)            :: e_si            !..saturation water pressure over ice
    REAL(wp)            :: e_sw            !..saturation water pressure over liquid
    REAL(wp)            :: e_d,p_a,dep_sum !,weight

    IF (isdebug) CALL message(routine, "vapor_deposition_growth")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend
          p_a  = atmo%p(i,k)
          T_a  = atmo%T(i,k)
          IF (T_a < T_3) THEN
             e_d  = atmo%qv(i,k) * R_d * T_a
             e_si = e_es(T_a)
             e_sw = e_ws(T_a)
             s_si(i,k) = e_d / e_si - 1.0    !..supersaturation over ice
             D_vtp = diffusivity(T_a,p_a)    !  D_v = 8.7602e-5 * T_a**(1.81) / p_a
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

    CALL vapor_deposition_generic(ik_slice, ice, vid_params, g_i, s_si, &
         dt_local, dep_ice)
    CALL vapor_deposition_generic(ik_slice, snow, vsd_params, g_i, s_si, &
         dt_local, dep_snow)
    CALL vapor_deposition_graupel()
    IF (ice_typ > 1) CALL vapor_deposition_hail()

    zdt = 1.0/dt_local

    DO k = kstart,kend
       DO i = istart,iend

          T_a  = atmo%T(i,k)

          ! Deposition only below T_3, evaporation of melting particles at warmer T is treated elsewhere
          IF (T_a < T_3) THEN

             ! Depositional growth with relaxation time-scale approach based on:
             ! "A New Double-Moment Microphysics Parameterization for Application in Cloud and
             ! Climate Models. Part 1: Description" by H. Morrison, J.A.Curry, V.I. Khvorostyanov

             qvsidiff  = atmo%qv(i,k) - e_es(T_a)/(R_d*T_a)

             if (abs(qvsidiff).gt.eps) then

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

                ! this limiter should not be necessary
                IF (qvsidiff < 0.0_wp) THEN
                   dep_ice(i,k)     = MAX(dep_ice(i,k),    -ice%q(i,k))
                   dep_snow(i,k)    = MAX(dep_snow(i,k),   -snow%q(i,k))
                   dep_graupel(i,k) = MAX(dep_graupel(i,k),-graupel%q(i,k))
                   dep_hail(i,k)    = MAX(dep_hail(i,k),   -hail%q(i,k))
                END IF

                dep_sum = dep_ice(i,k) + dep_graupel(i,k) + dep_snow(i,k) + dep_hail(i,k)

                ! this limiter should not be necessary
                !IF (qvsidiff > 0.0_wp .and. dep_sum > qvsidiff) then
                !   weight = qvsidiff / dep_sum
                !   dep_sum          = weight * dep_sum
                !   dep_ice(i,k)     = weight * dep_ice(i,k)
                !   dep_snow(i,k)    = weight * dep_snow(i,k)
                !   dep_graupel(i,k) = weight * dep_graupel(i,k)
                !   dep_hail(i,k)    = weight * dep_hail(i,k)
                !END IF

                ice%q(i,k)     = ice%q(i,k)     + dep_ice(i,k)
                snow%q(i,k)    = snow%q(i,k)    + dep_snow(i,k)
                graupel%q(i,k) = graupel%q(i,k) + dep_graupel(i,k)
                IF (ice_typ > 1) THEN
                   hail%q(i,k)    = hail%q(i,k)    + dep_hail(i,k)
                END IF

                atmo%qv(i,k) = atmo%qv(i,k) - dep_sum

                dep_rate_ice(i,k)  = dep_rate_ice(i,k)  + dep_ice(i,k)
                dep_rate_snow(i,k) = dep_rate_snow(i,k) + dep_snow(i,k)

             END IF

          ENDIF
       ENDDO
    ENDDO

  CONTAINS

    SUBROUTINE vapor_deposition_ice()
      INTEGER             :: i,k
      REAL(wp)            :: q_i,n_i,x_i,d_i,v_i,f_v

      DO k = kstart,kend
         DO i = istart,iend

            IF (ice%q(i,k) == 0.0_wp) THEN
               dep_ice(i,k) = 0.0_wp
            ELSE
              n_i = ice%n(i,k)
              q_i = ice%q(i,k)
              x_i = ice%meanmass(q_i,n_i)
              D_i = ice%diameter(x_i)
              v_i = ice%velocity(x_i) * ice%rho_v(i,k)

              !..note that a_f includes more than just ventilation, do never ever set f_v=1
              f_v  = vid_params%a_f + vid_params%b_f * sqrt(v_i*d_i)
              f_v  = MAX(f_v,vid_params%a_f/ice%a_ven)

              dep_ice(i,k) = g_i(i,k) * n_i * vid_params%c * d_i * f_v * s_si(i,k) * dt_local
           ENDIF
        ENDDO
     ENDDO
   END SUBROUTINE vapor_deposition_ice

   SUBROUTINE vapor_deposition_graupel()
     INTEGER             :: i,k
     REAL(wp)            :: q_g,n_g,x_g,d_g,v_g,f_v

     IF (isdebug) CALL message(routine, "vapor_deposition_graupel")

      DO k = kstart,kend
         DO i = istart,iend

            IF (graupel%q(i,k) == 0.0_wp) THEN
               dep_graupel(i,k) = 0.0_wp
            ELSE
               n_g = graupel%n(i,k)
               q_g = graupel%q(i,k)
               x_g = graupel%meanmass(q_g,n_g)
               d_g = graupel%diameter(x_g)
               v_g = graupel%velocity(x_g) * graupel%rho_v(i,k)

               f_v  = vgd_params%a_f + vgd_params%b_f * sqrt(v_g*d_g)
               f_v  = MAX(f_v,vgd_params%a_f/graupel%a_ven)

               dep_graupel(i,k) = g_i(i,k) * n_g * vgd_params%c * d_g * f_v * s_si(i,k) * dt_local
            ENDIF
         ENDDO
      ENDDO
    END SUBROUTINE vapor_deposition_graupel

    SUBROUTINE vapor_deposition_hail()
      INTEGER             :: i,k
      REAL(wp)            :: q_h,n_h,x_h,d_h,v_h,f_v

      IF (isdebug) CALL message(routine, "vapor_deposition_hail")

      DO k = kstart,kend
         DO i = istart,iend

            IF (hail%q(i,k) == 0.0_wp) THEN
               dep_hail(i,k)   = 0.0_wp
            ELSE
               n_h = hail%n(i,k)
               q_h = hail%q(i,k)
               x_h = hail%meanmass(q_h,n_h)
               D_h = hail%diameter(x_h)
               v_h = hail%velocity(x_h) * hail%rho_v(i,k)

               f_v  = vhd_params%a_f + vhd_params%b_f * sqrt(v_h*d_h)
               f_v  = MAX(f_v,vhd_params%a_f/hail%a_ven)

               dep_hail(i,k) = g_i(i,k) * n_h * vhd_params%c * d_h * f_v * s_si(i,k) * dt_local
            ENDIF
         ENDDO
      ENDDO
    END SUBROUTINE vapor_deposition_hail

    SUBROUTINE vapor_deposition_snow()
      REAL(wp)            :: q_s,n_s,x_s,d_s,v_s,f_v
      INTEGER             :: i,k

      IF (isdebug) CALL message(routine, "vapor_depositiosnow%n")

      DO k = kstart,kend
         DO i = istart,iend

            IF (snow%q(i,k) == 0.0_wp) THEN
               dep_snow(i,k) = 0.0_wp
            ELSE
               n_s = snow%n(i,k)
               q_s = snow%q(i,k)

               x_s = snow%meanmass(q_s,n_s)
               D_s = snow%diameter(x_s)
               v_s = snow%velocity(x_s) * snow%rho_v(i,k)

               f_v  = vsd_params%a_f + vsd_params%b_f * sqrt(D_s*v_s)
               f_v  = MAX(f_v,vsd_params%a_f/snow%a_ven)

               dep_snow(i,k) = g_i(i,k) * n_s * vsd_params%c * d_s * f_v * s_si(i,k) * dt_local
            ENDIF

         ENDDO
      ENDDO
    END SUBROUTINE vapor_deposition_snow

  END SUBROUTINE vapor_dep_relaxation

  SUBROUTINE vapor_deposition_generic(ik_slice, prtcl, params, g_i, s_si, &
       dt, dep)
    INTEGER, INTENT(in) :: ik_slice(4)
    TYPE(particle), INTENT(in) :: prtcl
    TYPE(evaporation_deposition_params), INTENT(in) :: params
    REAL(wp), INTENT(in) :: g_i(:, :), s_si(:, :)
    REAL(wp), INTENT(in) :: dt
    REAL(wp), INTENT(out) :: dep(:, :)
    REAL(wp)            :: q,n,x,d,v,f_v
    INTEGER             :: i,k
    INTEGER :: istart, iend, kstart, kend

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
      DO i = istart,iend
        IF (prtcl%q(i,k) == 0.0_wp) THEN
          dep(i,k) = 0.0_wp
        ELSE
          n = prtcl%n(i,k)
          q = prtcl%q(i,k)

          x = prtcl%meanmass(q,n)
          D = prtcl%diameter(x)
          v = particle_velocity(prtcl, x) * prtcl%rho_v(i,k)

          f_v  = params%a_f + params%b_f * SQRT(D*v)
          f_v  = MAX(f_v,params%a_f/prtcl%a_ven)

          dep(i,k) = g_i(i,k) * n * params%c * d * f_v * s_si(i,k) * dt
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE vapor_deposition_generic

  SUBROUTINE rain_freeze_gamlook(ik_slice, dt, atmo,rain,ice,snow,graupel,hail)
    !*******************************************************************************
    ! Freezing of raindrops                                                        *
    ! by Uli Blahak                                                                *
    !                                                                              *
    ! incomplete gamma functions are implemented as look-up tables                 *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: rain, ice, snow, graupel, hail

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: fr_q,fr_n,T_a,q_r,x_r,n_r,j_het,               &
         &                 fr_q_i,fr_n_i,fr_q_g,fr_n_g,fr_q_h,fr_n_h,n_0, &
         &                 lam,xmax_ice,xmax_gr,fr_q_tmp,fr_n_tmp

    REAL(wp), PARAMETER :: a_HET = 6.5d-1      ! Data of Barklie and Gokhale (PK S.350)
    REAL(wp), PARAMETER :: b_HET = 2.0d+2      !         Barklie and Gokhale (PK S.350)
    REAL(wp), PARAMETER :: eps = 1e-15_wp      ! for clipping
    LOGICAL,  PARAMETER :: lclipping = .true.

    INTEGER, SAVE       :: firstcall
    REAL(wp), SAVE      :: coeff_z
!$omp threadprivate (firstcall)
!$omp threadprivate (coeff_z)

    IF (firstcall.NE.1) THEN
      firstcall = 1
      coeff_z = moment_gamma(rain,2)  ! coeff for 2nd moment
      IF (isdebug) THEN
        CALL message(routine, "rain_freeze_gamlook:")
        WRITE(txt,'(A,D10.3)') "    coeff_z= ",coeff_z ; CALL message(routine,TRIM(txt))
      ENDIF
    ELSE IF (isdebug) THEN
      CALL message(routine, "rain_freeze_gamlook")
    ENDIF

    xmax_ice = ( (D_rainfrz_ig/rain%a_geo)**(1.0_wp/rain%b_geo) )**rain%mu
    xmax_gr  = ( (D_rainfrz_gh/rain%a_geo)**(1.0_wp/rain%b_geo) )**rain%mu

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)
          q_r = rain%q(i,k)
          n_r = rain%n(i,k)

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
                   lam = exp( log( rain_g1/rain_g2*x_r) * (-rain%mu) )
                   n_0 = rain%mu * n_r * lam**(rain_nm1) / rain_g1
                   fr_n_i = n_0/(rain%mu*lam**(rain_nm1)) * incgfct_lower_lookup(lam*xmax_ice,rain_ltable1)
                   fr_q_i = n_0/(rain%mu*lam**(rain_nm2)) * incgfct_lower_lookup(lam*xmax_ice,rain_ltable2)
                   fr_n_g = n_0/(rain%mu*lam**(rain_nm1)) * incgfct_lower_lookup(lam*xmax_gr, rain_ltable1)
                   fr_q_g = n_0/(rain%mu*lam**(rain_nm2)) * incgfct_lower_lookup(lam*xmax_gr, rain_ltable2)

                   fr_n_h = fr_n - fr_n_g
                   fr_q_h = fr_q - fr_q_g
                   fr_n_g = fr_n_g - fr_n_i
                   fr_q_g = fr_q_g - fr_q_i
                   fr_n_tmp = n_r/MAX(fr_n,n_r)
                   fr_q_tmp = q_r/MAX(fr_q,q_r)
                ELSE
                   !..heterogeneous freezing
                   j_het = MAX(b_HET * ( EXP( a_HET * (T_3 - T_a)) - 1.0 ),0.0_wp) / rho_w * dt
!!                   if (use_prog_in) j_het = MIN(j_het, n_inact(i,k)/q_r)

                   ! Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                   ! oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                   ! bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                   ! und von xmax_gr bis unendlich (--> Hagel).
                   IF (j_het >= 1d-20) THEN
                      fr_n  = j_het * q_r
                      fr_q  = j_het * q_r * x_r * coeff_z

                      lam = ( rain_g1 / rain_g2 * x_r)**(-rain%mu)
                      n_0 = rain%mu * n_r * lam**(rain_nm1) / rain_g1
                      fr_n_i = j_het * n_0/(rain%mu*lam**(rain_nm2)) * incgfct_lower_lookup(lam*xmax_ice,rain_ltable2)
                      fr_q_i = j_het * n_0/(rain%mu*lam**(rain_nm3)) * incgfct_lower_lookup(lam*xmax_ice,rain_ltable3)
                      fr_n_g = j_het * n_0/(rain%mu*lam**(rain_nm2)) * incgfct_lower_lookup(lam*xmax_gr, rain_ltable2)
                      fr_q_g = j_het * n_0/(rain%mu*lam**(rain_nm3)) * incgfct_lower_lookup(lam*xmax_gr, rain_ltable3)

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

             rain%q(i,k) = rain%q(i,k) - fr_q
             rain%n(i,k) = n_r - fr_n
             !if (use_prog_in) then
             !   n_inact(i,k) = n_inact(i,k) + fr_n
             !end if

             IF (ice_typ < 2) THEN
                ! ohne Hagelklasse,  gefrierender Regen wird Eis oder Graupel
                ice%q(i,k) = ice%q(i,k)  + fr_q_i
                ice%n(i,k) = ice%n(i,k)  + fr_n_i
                graupel%q(i,k) = graupel%q(i,k)  + fr_q_h + fr_q_g
                graupel%n(i,k) = graupel%n(i,k)  + fr_n_h + fr_n_g
             ELSE
                ! mit Hagelklasse, gefrierender Regen wird Eis, Graupel oder Hagel
                snow%q(i,k) = snow%q(i,k)  + fr_q_i
                snow%n(i,k) = snow%n(i,k)  + fr_n_i   ! put this into snow
                !ice%q(i,k) = ice%q(i,k)  + fr_q_i    ! ... or into ice?
                !ice%n(i,k) = ice%n(i,k)  + fr_n_i
                graupel%q(i,k) = graupel%q(i,k)  + fr_q_g
                graupel%n(i,k) = graupel%n(i,k)  + fr_n_g
                hail%q(i,k) = hail%q(i,k)  + fr_q_h
                hail%n(i,k) = hail%n(i,k)  + fr_n_h
             ENDIF

             ! clipping of small negatives is necessary here
             if (lclipping) then
                IF (rain%q(i,k) < 0.0 .and. abs(rain%q(i,k)) < eps) rain%q(i,k) = 0.0_wp
                IF (rain%n(i,k) < 0.0 .and. abs(rain%n(i,k)) < eps) rain%n(i,k) = 0.0_wp
                IF (graupel%q(i,k) < 0.0 .and. abs(graupel%q(i,k)) < eps) graupel%q(i,k) = 0.0_wp
                IF (graupel%n(i,k) < 0.0 .and. abs(graupel%q(i,k)) < eps) graupel%n(i,k) = 0.0_wp
                IF (hail%q(i,k) < 0.0 .and. abs(hail%q(i,k)) < eps) hail%q(i,k) = 0.0_wp
                IF (hail%n(i,k) < 0.0 .and. abs(hail%n(i,k)) < eps) hail%n(i,k) = 0.0_wp
             end if

          END IF
       END DO
    END DO
  END SUBROUTINE rain_freeze_gamlook

  SUBROUTINE ice_selfcollection(ik_slice, dt, atmo, ice, snow)
    !*******************************************************************************
    ! selfcollection of ice crystals, see SB2006 or Seifert (2002)                 *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: ice, snow
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
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
!$omp threadprivate (firstcall)
!$omp threadprivate (delta_n)
!$omp threadprivate (delta_q)
!$omp threadprivate (theta_n)
!$omp threadprivate (theta_q)

    IF (isdebug) CALL message(routine, "ice_selfcollection")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

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

          q_i = ice%q(i,k)
          n_i = ice%n(i,k)

          x_i = ice%meanmass(q_i,n_i)
          D_i = ice%diameter(x_i)

          IF ( n_i > 0.0_wp .AND. q_i > q_crit_ii .AND. D_i > D_crit_ii ) THEN

             T_a = atmo%T(i,k)

             !.. Temperaturabhaengige Efficiency nach Cotton et al. (1986)
             !   (siehe auch Straka, 1989; S. 53)
             e_coll = MIN(10**(0.035*(T_a-T_3)-0.7),0.2_wp)

             v_i = ice%a_vel * x_i**ice%b_vel * ice%rho_v(i,k)

             self_n = pi4 * e_coll * delta_n * n_i * n_i * D_i * D_i &
                  & * sqrt( theta_n * v_i * v_i + 2.0*ice_s_vel**2 ) * dt

             self_q = pi4 * e_coll * delta_q * n_i * q_i * D_i * D_i &
                  & * sqrt( theta_q * v_i * v_i + 2.0*ice_s_vel**2 ) * dt

             self_q = MIN(self_q,q_i)
             self_n = MIN(MIN(self_n,self_q/x_conv_ii),n_i)

             ice%q(i,k)  = ice%q(i,k)  - self_q
             snow%q(i,k) = snow%q(i,k) + self_q

             ice%n(i,k)  = ice%n(i,k)  - self_n
             snow%n(i,k) = snow%n(i,k) + self_n / 2.0

          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE ice_selfcollection

  SUBROUTINE setup_snow_selfcollection(snow)
    TYPE(particle), INTENT(in) :: snow
    REAL(wp) :: delta_n_11,delta_n_12
    REAL(wp) :: theta_n_11,theta_n_12

    delta_n_11 = coll_delta_11(snow,snow,0)
    delta_n_12 = coll_delta_12(snow,snow,0)
    theta_n_11 = coll_theta_11(snow,snow,0)
    theta_n_12 = coll_theta_12(snow,snow,0)

    snow_sc_delta_n = (2.0*delta_n_11 + delta_n_12)
    snow_sc_theta_n = (2.0*theta_n_11 - theta_n_12)

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "    a_snow     = ",snow%a_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    b_snow     = ",snow%b_geo ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    alf_snow   = ",snow%a_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    bet_snow   = ",snow%b_vel ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_11 = ",delta_n_11  ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_12 = ",delta_n_12  ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n    = ",snow_sc_delta_n     ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_11 = ",theta_n_11  ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_12 = ",theta_n_12  ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n    = ",snow_sc_theta_n     ; CALL message(routine,TRIM(txt))
    END IF
  END SUBROUTINE setup_snow_selfcollection

  SUBROUTINE snow_selfcollection(ik_slice, dt, atmo, snow)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: snow
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s,e_coll
    REAL(wp)            :: self_n

    IF (isdebug) CALL message(routine, "snow_selfcollection")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_s = snow%q(i,k)
          T_a = atmo%T(i,k)

          IF ( q_s > q_crit ) THEN

             !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
             e_coll = MAX(0.1_wp,MIN(EXP(0.09*(T_a-T_3)),1.0_wp))

             n_s = snow%n(i,k)
             x_s = snow%meanmass(q_s,n_s)
             D_s = snow%diameter(x_s)
             v_s = snow%velocity(x_s) * snow%rho_v(i,k)

             self_n = pi8 * e_coll * n_s * n_s * snow_sc_delta_n * D_s * D_s * &
                  &          sqrt(  snow_sc_theta_n * v_s * v_s                &
                  &               + 2.0 * snow_s_vel**2 ) * dt

             self_n = MIN(self_n,n_s)

             snow%n(i,k) = snow%n(i,k) - self_n

         ENDIF
      ENDDO
   ENDDO

  END SUBROUTINE snow_selfcollection

  SUBROUTINE snow_melting(ik_slice, dt, atmo,snow,rain)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: rain, snow
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s,T_a,e_a
    REAL(wp)            :: melt,melt_v,melt_h,melt_n,melt_q
    REAL(wp)            :: fh_q,fv_q

    IF (isdebug) CALL message(routine, "snow_melting")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)
          q_s = snow%q(i,k)
          IF (T_a > T_3 .AND. q_s > 0.0_wp) THEN
            e_a = e_ws(T_a)            ! saturation pressure
            n_s = snow%n(i,k)

            x_s = snow%meanmass(q_s,n_s)
            D_s = snow%diameter(x_s)
            v_s = snow%velocity(x_s) * snow%rho_v(i,k)

            fv_q = sm_params%a_vent + sm_params%b_vent * sqrt(v_s*D_s)

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

            snow%q(i,k) = snow%q(i,k) - melt_q
            rain%q(i,k) = rain%q(i,k) + melt_q

            snow%n(i,k) = snow%n(i,k) - melt_n
            rain%n(i,k) = rain%n(i,k) + melt_n

            snow%n(i,k) = MAX(snow%n(i,k), snow%q(i,k)/snow%x_max)

         ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE snow_melting

  SUBROUTINE graupel_snow_collection(ik_slice, dt, atmo, snow, graupel)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: snow, graupel

    ! Locale Variablen
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s
    REAL(wp)            :: coll_n,coll_q,e_coll

    IF (isdebug) CALL message(routine, "graupel_snow_collection")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_s = snow%q(i,k)
          q_g = graupel%q(i,k)

          IF (q_s > q_crit .AND. q_g > q_crit) THEN
            T_a = atmo%T(i,k)

            !.. sticking efficiency of Lin (1983)
            e_coll = min(exp(0.09*(T_a-T_3)),1.0_wp)

            n_s = snow%n(i,k)
            n_g = graupel%n(i,k)

            x_g = graupel%meanmass(q_g,n_g)
            d_g = graupel%diameter(x_g)
            v_g = graupel%velocity(x_g) * graupel%rho_v(i,k)

            x_s = snow%meanmass(q_s,n_s)
            d_s = snow%diameter(x_s)
            v_s = snow%velocity(x_s) * snow%rho_v(i,k)

            coll_n = pi4 * n_g * n_s * e_coll * dt &
                 & *     (gsc_params%delta_n_aa * D_g**2 + gsc_params%delta_n_ab * D_g*D_s + gsc_params%delta_n_bb * D_s**2) &
                 & * sqrt(gsc_params%theta_n_aa * v_g**2 - gsc_params%theta_n_ab * v_g*v_s + gsc_params%theta_n_bb * v_s**2  &
                 &       +snow_s_vel**2)

            coll_q = pi4 * n_g * q_s * e_coll * dt &
                 & *     (gsc_params%delta_q_aa * D_g**2 + gsc_params%delta_q_ab * D_g*D_s + gsc_params%delta_q_bb * D_s**2) &
                 & * sqrt(gsc_params%theta_q_aa * v_g**2 - gsc_params%theta_q_ab * v_g*v_s + gsc_params%theta_q_bb * v_s**2  &
                 &       +snow_s_vel**2)

            coll_n = MIN(n_s,coll_n)
            coll_q = MIN(q_s,coll_q)

            graupel%q(i,k) = graupel%q(i,k) + coll_q
            snow%q(i,k)    = snow%q(i,k)    - coll_q
            snow%n(i,k)    = snow%n(i,k)    - coll_n

         ENDIF
      ENDDO
   ENDDO

  END SUBROUTINE graupel_snow_collection

  SUBROUTINE hail_snow_collection(ik_slice, dt, atmo, snow, hail)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: snow, hail
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s
    REAL(wp)            :: coll_n,coll_q,e_coll

    IF (isdebug) CALL message(routine, "hail_snow_collection")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_s = snow%q(i,k)
          q_h = hail%q(i,k)

          IF (q_s > q_crit .AND. q_h > q_crit) THEN
            T_a = atmo%T(i,k)

            !.. sticking efficiency of Lin (1983)
            e_coll = min(exp(0.09*(T_a-T_3)),1.0_wp)

            n_s = snow%n(i,k)
            n_h = hail%n(i,k)

            x_s = snow%meanmass(q_s,n_s)
            D_s = snow%diameter(x_s)
            v_s = snow%velocity(x_s) * snow%rho_v(i,k)

            x_h = hail%meanmass(q_h,n_h)
            D_h = hail%diameter(x_h)
            v_h = hail%velocity(x_h) * hail%rho_v(i,k)

            coll_n = pi4 * n_h * n_s * e_coll * dt &
                 & *     (hsc_params%delta_n_aa * D_h**2 + hsc_params%delta_n_ab * D_h*D_s + hsc_params%delta_n_bb * D_s**2) &
                 & * sqrt(hsc_params%theta_n_aa * v_h**2 - hsc_params%theta_n_ab * v_h*v_s + hsc_params%theta_n_bb * v_s**2  &
                 &       +snow_s_vel**2)

            coll_q = pi4 * n_h * q_s * e_coll * dt &
                 & *     (hsc_params%delta_q_aa * D_h**2 + hsc_params%delta_q_ab * D_h*D_s + hsc_params%delta_q_bb * D_s**2) &
                 & * sqrt(hsc_params%theta_q_aa * v_h**2 - hsc_params%theta_q_ab * v_h*v_s + hsc_params%theta_q_bb * v_s**2  &
                 &       +snow_s_vel**2)

            coll_n = MIN(n_s,coll_n)
            coll_q = MIN(q_s,coll_q)

            hail%q(i,k) = hail%q(i,k) + coll_q
            snow%q(i,k) = snow%q(i,k) - coll_q
            snow%n(i,k) = snow%n(i,k) - coll_n

          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE hail_snow_collection

  SUBROUTINE graupel_ice_collection(ik_slice, dt, atmo, ice, graupel)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: ice, graupel
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
    REAL(wp)            :: coll_n,coll_q,e_coll

    IF (isdebug) CALL message(routine, "graupel_ice_collection")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_i = ice%q(i,k)
          q_g = graupel%q(i,k)
          IF (q_i > q_crit .AND. q_g > q_crit) THEN
            T_a = atmo%T(i,k)

            !.. sticking efficiency of Lin (1983)
            e_coll = min(exp(0.09*(T_a-T_3)),1.0_wp)

            n_i = ice%n(i,k)
            n_g = graupel%n(i,k)

            x_g = graupel%meanmass(q_g,n_g)
            d_g = graupel%diameter(x_g)
            v_g = graupel%velocity(x_g) * graupel%rho_v(i,k)

            x_i = ice%meanmass(q_i,n_i)
            d_i = ice%diameter(x_i)
            v_i = ice%velocity(x_i) * ice%rho_v(i,k)

            coll_n = pi4 * n_g * n_i * e_coll * dt &
                 & *     (gic_params%delta_n_aa * D_g**2 + gic_params%delta_n_ab * D_g*D_i + gic_params%delta_n_bb * D_i**2) &
                 & * sqrt(gic_params%theta_n_aa * v_g**2 - gic_params%theta_n_ab * v_g*v_i + gic_params%theta_n_bb * v_i**2)

            coll_q = pi4 * n_g * q_i * e_coll * dt &
                 & *     (gic_params%delta_q_aa * D_g**2 + gic_params%delta_q_ab * D_g*D_i + gic_params%delta_q_bb * D_i**2) &
                 & * sqrt(gic_params%theta_q_aa * v_g**2 - gic_params%theta_q_ab * v_g*v_i + gic_params%theta_q_bb * v_i**2)

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            graupel%q(i,k) = graupel%q(i,k) + coll_q
            ice%q(i,k)     = ice%q(i,k)     - coll_q
            ice%n(i,k)     = ice%n(i,k)     - coll_n

          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE graupel_ice_collection

  SUBROUTINE hail_ice_collection(ik_slice, dt, atmo, ice, hail)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: ice, hail
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h
    REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
    REAL(wp)            :: coll_n,coll_q,e_coll

    IF (isdebug) CALL message(routine, "hail_ice_collection")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_i = ice%q(i,k)
          q_h = hail%q(i,k)
          IF (q_i > q_crit .AND. q_h > q_crit) THEN
            T_a = atmo%T(i,k)

            !.. sticking efficiency of Lin (1983)
            e_coll = min(exp(0.09*(T_a-T_3)),1.0_wp)

            n_i = ice%n(i,k)
            n_h = hail%n(i,k)

            x_h = hail%meanmass(q_h,n_h)
            D_h = hail%diameter(x_h)
            v_h = hail%velocity(x_h) * hail%rho_v(i,k)

            x_i = ice%meanmass(q_i,n_i)
            D_i = ice%diameter(x_i)
            v_i = ice%velocity(x_i) * ice%rho_v(i,k)

            coll_n = pi4 * n_h * n_i * e_coll * dt &
                 & *     (hic_params%delta_n_aa * D_h**2 + hic_params%delta_n_ab * D_h*D_i + hic_params%delta_n_bb * D_i**2) &
                 & * sqrt(hic_params%theta_n_aa * v_h**2 - hic_params%theta_n_ab * v_h*v_i + hic_params%theta_n_bb * v_i**2)

            coll_q = pi4 * n_h * q_i * e_coll * dt &
                 & *     (hic_params%delta_q_aa * D_h**2 + hic_params%delta_q_ab * D_h*D_i + hic_params%delta_q_bb * D_i**2) &
                 & * sqrt(hic_params%theta_q_aa * v_h**2 - hic_params%theta_q_ab * v_h*v_i + hic_params%theta_q_bb * v_i**2)

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            hail%q(i,k) = hail%q(i,k) + coll_q
            ice%q(i,k)  = ice%q(i,k)  - coll_q
            ice%n(i,k)  = ice%n(i,k)  - coll_n

          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE hail_ice_collection

  SUBROUTINE snow_ice_collection(ik_slice, dt, atmo,ice,snow)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: ice, snow
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_s,n_s,x_s,d_s,v_s
    REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
    REAL(wp)            :: coll_n,coll_q,e_coll

    IF (isdebug) CALL message(routine, "snow_ice_collection")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_i = ice%q(i,k)
          q_s = snow%q(i,k)
          IF (q_i > q_crit .AND. q_s > q_crit) THEN
            T_a = atmo%T(i,k)

            !.. sticking efficiency of Lin (1983)
            e_coll = max(0.1_wp,min(exp(0.09*(T_a-T_3)),1.0_wp))

            n_i = ice%n(i,k)
            n_s = snow%n(i,k)

            x_i = ice%meanmass(q_i,n_i)
            d_i = ice%diameter(x_i)
            v_i = ice%velocity(x_i) * ice%rho_v(i,k)

            x_s = snow%meanmass(q_s,n_s)
            d_s = snow%diameter(x_s)
            v_s = snow%velocity(x_s) * snow%rho_v(i,k)

            coll_n = pi4 * n_s * n_i * e_coll * dt &
                 & *     (sic_params%delta_n_aa * D_s**2 + sic_params%delta_n_ab * D_s*D_i + sic_params%delta_n_bb * D_i**2) &
                 & * sqrt(sic_params%theta_n_aa * v_s**2 - sic_params%theta_n_ab * v_s*v_i + sic_params%theta_n_bb * v_i**2  &
                 &       +snow_s_vel**2 + ice_s_vel**2)

            coll_q = pi4 * n_s * q_i * e_coll * dt &
                 & *     (sic_params%delta_q_aa * D_s**2 + sic_params%delta_q_ab * D_s*D_i + sic_params%delta_q_bb * D_i**2) &
                 & * sqrt(sic_params%theta_q_aa * v_s**2 - sic_params%theta_q_ab * v_s*v_i + sic_params%theta_q_bb * v_i**2  &
                 &       +snow_s_vel**2 + ice_s_vel**2)

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            snow%q(i,k) = snow%q(i,k) + coll_q
            ice%q(i,k)  = ice%q(i,k)  - coll_q
            ice%n(i,k)  = ice%n(i,k)  - coll_n
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE snow_ice_collection

  SUBROUTINE setup_graupel_selfcollection(graupel)
    TYPE(particle), INTENT(in) :: graupel
    REAL(wp) :: delta_n_11,delta_n_12
    REAL(wp) :: theta_n_11,theta_n_12
    REAL(wp) :: delta_n, theta_n

    delta_n_11 = coll_delta_11(graupel,graupel,0)
    delta_n_12 = coll_delta_12(graupel,graupel,0)
    theta_n_11 = coll_theta_11(graupel,graupel,0)
    theta_n_12 = coll_theta_12(graupel,graupel,0)

    delta_n = (2.0*delta_n_11 + delta_n_12)
    theta_n = (2.0*theta_n_11 - theta_n_12)**0.5
    graupel_sc_coll_n  = pi8 * delta_n * theta_n

    IF (isdebug) THEN
      WRITE(txt,'(A,D10.3)') "    delta_n_11 = ",delta_n_11 ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n_12 = ",delta_n_12 ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    delta_n    = ",delta_n    ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_11 = ",theta_n_11 ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n_12 = ",theta_n_12 ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    theta_n    = ",theta_n    ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "    coll_n     = ",graupel_sc_coll_n     ; CALL message(routine,TRIM(txt))
    END IF
  END SUBROUTINE setup_graupel_selfcollection


  SUBROUTINE graupel_selfcollection(ik_slice, dt, atmo, graupel)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: graupel
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: self_n

    IF (isdebug) CALL message(routine, "graupel_selfcollection")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_g = graupel%q(i,k)
          IF ( q_g > q_crit ) THEN

             n_g = graupel%n(i,k)
             x_g = graupel%meanmass(q_g,n_g)
             d_g = graupel%diameter(x_g)
             v_g = graupel%velocity(x_g) * graupel%rho_v(i,k)

             self_n = graupel_sc_coll_n * n_g**2 * D_g**2 * v_g * dt

             ! sticking efficiency does only distinguish dry and wet based on T_3
             IF (atmo%T(i,k) > T_3) THEN
                self_n = self_n * ecoll_gg_wet
             ELSE
                self_n = self_n * ecoll_gg
             END IF

             self_n = MIN(self_n,n_g)

             graupel%n(i,k) = graupel%n(i,k) - self_n
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE graupel_selfcollection

  SUBROUTINE ice_melting(ik_slice, atmo, ice, cloud, rain)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud, rain, ice
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_i,x_i,n_i
    REAL(wp)            :: melt_q,melt_n

    IF (isdebug) CALL message(routine, "ice_melting")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_i = ice%q(i,k)

          IF (atmo%T(i,k) > T_3 .AND. q_i > 0.0) THEN

            n_i = ice%n(i,k)
            x_i = ice%meanmass(q_i,n_i)

            ! complete melting within this time step
            melt_q = q_i
            melt_n = n_i
            ice%q(i,k) = 0.0_wp
            ice%n(i,k) = 0.0_wp

            ! ice either melts into cloud droplets or rain depending on x_i
            IF (x_i > cloud%x_max) THEN
               rain%q(i,k)  = rain%q(i,k)  + melt_q
               rain%n(i,k)  = rain%n(i,k)  + melt_n
            ELSE
               cloud%q(i,k) = cloud%q(i,k) + melt_q
               cloud%n(i,k) = cloud%n(i,k) + melt_n
            ENDIF

          END IF
       END DO
    END DO
  END SUBROUTINE ice_melting

  SUBROUTINE graupel_cloud_riming(ik_slice, dt, atmo, graupel, cloud, rain, ice)
    !*******************************************************************************
    ! see SB2006                                                                   *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud, rain, ice, graupel
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: q_c,n_c,x_c,d_c,v_c,x_coll_c
    REAL(wp)            :: rime_n,rime_q,e_coll_n
    REAL(wp)            :: melt_n,melt_q,e_coll_q
    REAL(wp)            :: shed_n,shed_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp)            :: const1
    REAL(wp), PARAMETER :: const2 = 1.0/(T_mult_opt - T_mult_min), &
         const3 = 1.0/(T_mult_opt - T_mult_max), &
         const4 = c_w / L_ew, &
         !..mean mass of shedding drops
         x_shed = 4./3.*pi * rho_w * r_shedding**3

    IF (isdebug) CALL message(routine, "graupel_cloud_riming")

    x_coll_c = (D_coll_c/cloud%a_geo)**3          !..lower threshold for collection
    const1 = ecoll_gc/(D_coll_c - D_crit_c)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          q_g = graupel%q(i,k)
          n_c = cloud%n(i,k)
          n_g = graupel%n(i,k)

          x_g = graupel%meanmass(q_g,n_g)
          D_g = graupel%diameter(x_g)
          x_c = cloud%meanmass(q_c,n_c)
          D_c = cloud%diameter(x_c)

          T_a = atmo%T(i,k)
          IF (q_c > q_crit_c .AND. q_g > q_crit_gc .AND. D_g > D_crit_gc .AND. &
               & D_c > D_crit_c) THEN

             v_g = graupel%velocity(x_g) * graupel%rho_v(i,k)
             v_c = cloud%velocity(x_c)   * cloud%rho_v(i,k)

             e_coll_n = MIN(ecoll_gc, MAX(const1*(D_c - D_crit_c),ecoll_min))
             e_coll_q = e_coll_n

             rime_n = pi4 * e_coll_n * n_g * n_c * dt &
                  & *     (gcr_params%delta_n_aa * D_g*D_g + gcr_params%delta_n_ab * D_g*D_c + gcr_params%delta_n_bb * D_c*D_c) &
                  & * sqrt(gcr_params%theta_n_aa * v_g*v_g - gcr_params%theta_n_ab * v_g*v_c + gcr_params%theta_n_bb * v_c*v_c)

             rime_q = pi4 * e_coll_q * n_g * q_c * dt &
                  & *     (gcr_params%delta_q_aa * D_g*D_g + gcr_params%delta_q_ab * D_g*D_c + gcr_params%delta_q_bb * D_c*D_c) &
                  & * sqrt(gcr_params%theta_q_aa * v_g*v_g - gcr_params%theta_q_ab * v_g*v_c + gcr_params%theta_q_bb * v_c*v_c)

             rime_q = MIN(q_c,rime_q)
             rime_n = MIN(n_c,rime_n)

             graupel%q(i,k) = graupel%q(i,k) + rime_q
             cloud%q(i,k)   = cloud%q(i,k)   - rime_q
             cloud%n(i,k)   = cloud%n(i,k)   - rime_n

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

                ice%n(i,k)     = ice%n(i,k)     + mult_n
                ice%q(i,k)     = ice%q(i,k)     + mult_q
                graupel%q(i,k) = graupel%q(i,k) - mult_q
             ENDIF

             ! enhancement of melting of graupel
             IF (T_a > T_3 .AND. enhanced_melting) THEN
                melt_q = const4 * (T_a - T_3) * rime_q
                melt_n = melt_q / x_g

                melt_q = MIN(graupel%q(i,k),melt_q)
                melt_n = MIN(graupel%n(i,k),melt_n)

                graupel%q(i,k) = graupel%q(i,k) - melt_q
                rain%q(i,k)    = rain%q(i,k)    + melt_q
                graupel%n(i,k) = graupel%n(i,k) - melt_n
                rain%n(i,k)    = rain%n(i,k)    + melt_n
             ELSE
                melt_q = 0.0_wp
             ENDIF

             ! Shedding
             IF ((graupel_shedding .AND. D_g > D_shed_g .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
                q_g = graupel%q(i,k)
                n_g = graupel%n(i,k)
                x_g = graupel%meanmass(q_g,n_g)

                shed_q = MIN(q_g,rime_q)
                shed_n = shed_q / MIN(x_shed,x_g)

                graupel%q(i,k) = graupel%q(i,k) - shed_q
                rain%q(i,k)    = rain%q(i,k)    + shed_q
                rain%n(i,k)    = rain%n(i,k)    + shed_n
             ELSE
                shed_q = 0.0
             ENDIF
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE graupel_cloud_riming

  SUBROUTINE hail_cloud_riming(ik_slice, dt, atmo, hail, cloud, rain, ice)
    !*******************************************************************************
    ! see SB2006 for the equations                                                 *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud, rain, ice, hail
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h
    REAL(wp)            :: q_c,n_c,x_c,d_c,v_c,x_coll_c
    REAL(wp)            :: rime_n,rime_q,e_coll_n
    REAL(wp)            :: melt_n,melt_q,e_coll_q
    REAL(wp)            :: shed_n,shed_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp), PARAMETER :: &
         const1 = ecoll_hc/(D_coll_c - D_crit_c), &
         const2 = 1.0/(T_mult_opt - T_mult_min), &
         const3 = 1.0/(T_mult_opt - T_mult_max), &
         const4 = c_w / L_ew, &
         !..mean mass of shedding droplets
         x_shed = 4./3.*pi * rho_w * r_shedding**3


    IF (isdebug) CALL message(routine, "hail_cloud_riming")

    x_coll_c = (D_coll_c/cloud%a_geo)**3          !..lower mass for collection


    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)
          q_c = cloud%q(i,k)
          q_h = hail%q(i,k)
          n_c = cloud%n(i,k)
          n_h = hail%n(i,k)

          x_h = hail%meanmass(q_h,n_h)
          D_h = hail%diameter(x_h)
          x_c = cloud%meanmass(q_c,n_c)
          D_c = cloud%diameter(x_c)

          IF (q_c > q_crit_c .AND. q_h > q_crit_hc .AND. D_h > D_crit_hc .AND. &
               & D_c > D_crit_c) THEN

            v_h = hail%velocity(x_h)  * hail%rho_v(i,k)
            v_c = cloud%velocity(x_c) * cloud%rho_v(i,k)

            e_coll_n = MIN(ecoll_hc, MAX(const1*(D_c - D_crit_c),ecoll_min))
            e_coll_q = e_coll_n

            rime_n = pi4 * e_coll_n * n_h * n_c * dt &
                 &   * (hcr_params%delta_n_aa * D_h*D_h + hcr_params%delta_n_ab * D_h*D_c + hcr_params%delta_n_bb * D_c*D_c) &
                 &   * SQRT(hcr_params%theta_n_aa * v_h*v_h - hcr_params%theta_n_ab * v_h*v_c + hcr_params%theta_n_bb * v_c*v_c)

            rime_q = pi4 * e_coll_q * n_h * q_c * dt &
                 &   * (hcr_params%delta_q_aa * D_h*D_h + hcr_params%delta_q_ab * D_h*D_c + hcr_params%delta_q_bb * D_c*D_c) &
                 &   * SQRT(hcr_params%theta_q_aa * v_h*v_h - hcr_params%theta_q_ab * v_h*v_c + hcr_params%theta_q_bb * v_c*v_c)

            rime_q = MIN(q_c,rime_q)
            rime_n = MIN(n_c,rime_n)

            hail%q(i,k)  = hail%q(i,k) + rime_q
            cloud%q(i,k) = cloud%q(i,k) - rime_q
            cloud%n(i,k) = cloud%n(i,k) - rime_n

            ! ice multiplication
            IF (T_a < T_3 .AND. ice_multiplication) THEN
              mult_1 = const2*(T_a - T_mult_min)
              mult_2 = const3*(T_a - T_mult_max)
              mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
              mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
              mult_n = C_mult * mult_1 * mult_2 * rime_q
              mult_q = mult_n * ice%x_min
              mult_q = MIN(rime_q,mult_q)

              ice%n(i,k)  = ice%n(i,k)  + mult_n
              ice%q(i,k)  = ice%q(i,k)  + mult_q
              hail%q(i,k) = hail%q(i,k) - mult_q
            ENDIF

            ! enhancement of melting of hail
            IF (T_a > T_3 .AND. enhanced_melting) THEN
              melt_q = const4 * (T_a - T_3) * rime_q
              melt_n = melt_q / x_h

              melt_q = MIN(hail%q(i,k),melt_q)
              melt_n = MIN(hail%n(i,k),melt_n)

              hail%q(i,k) = hail%q(i,k) - melt_q
              rain%q(i,k) = rain%q(i,k) + melt_q
              hail%n(i,k) = hail%n(i,k) - melt_n
              rain%n(i,k) = rain%n(i,k) + melt_n
            ENDIF

            ! Shedding
            IF ((D_h > D_shed_h .AND. T_a > T_shed .AND. hail_shedding) .OR. T_a > T_3 ) THEN
              q_h = hail%q(i,k)
              n_h = hail%n(i,k)
              x_h = hail%meanmass(q_h,n_h)

              shed_q = MIN(q_h,rime_q)
              shed_n = shed_q / MIN(x_shed,x_h)

              hail%q(i,k) = hail%q(i,k) - shed_q
              rain%q(i,k) = rain%q(i,k) + shed_q
              rain%n(i,k) = rain%n(i,k) + shed_n
            ELSE
              shed_q = 0.0
            ENDIF

          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE hail_cloud_riming

  SUBROUTINE graupel_rain_riming(ik_slice, dt, atmo,graupel,rain,ice)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: rain, ice, graupel
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g
    REAL(wp)            :: q_r,n_r,x_r,d_r,v_r
    REAL(wp)            :: rime_n,rime_q,melt_n,melt_q,shed_n,shed_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp), PARAMETER :: &
         x_shed = 4./3. * pi * rho_w * r_shedding**3, &
         const3 = 1/(T_mult_opt - T_mult_min), &
         const4 = 1/(T_mult_opt - T_mult_max)


    IF (isdebug) CALL message(routine, "graupel_rain_riming")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_r = rain%q(i,k)
          q_g = graupel%q(i,k)
          T_a = atmo%T(i,k)

          IF (q_r > q_crit .AND. q_g > q_crit) THEN

            n_r = rain%n(i,k)
            n_g = graupel%n(i,k)

            x_g = graupel%meanmass(q_g,n_g)
            d_g = graupel%diameter(x_g)
            v_g = graupel%velocity(x_g) * graupel%rho_v(i,k)
            x_r = rain%meanmass(q_r,n_r)
            d_r = rain%diameter(x_r)
            v_r = rain%velocity(x_r) * rain%rho_v(i,k)

            rime_n = pi4 * n_g * n_r * dt &
                 & *     (grr_params%delta_n_aa * D_g**2 + grr_params%delta_n_ab * D_g*D_r + grr_params%delta_n_bb * D_r**2) &
                 & * sqrt(grr_params%theta_n_aa * v_g**2 - grr_params%theta_n_ab * v_g*v_r + grr_params%theta_n_bb * v_r**2)

            rime_q = pi4 * n_g * q_r * dt &
                 & *     (grr_params%delta_n_aa * D_g**2 + grr_params%delta_q_ab * D_g*D_r + grr_params%delta_q_bb * D_r**2) &
                 & * sqrt(grr_params%theta_n_aa * v_g**2 - grr_params%theta_q_ab * v_g*v_r + grr_params%theta_q_bb * v_r**2)

            rime_n = MIN(n_r,rime_n)
            rime_q = MIN(q_r,rime_q)

            graupel%q(i,k) = graupel%q(i,k) + rime_q
            rain%q(i,k)    = rain%q(i,k)    - rime_q
            rain%n(i,k)    = rain%n(i,k)    - rime_n

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

               ice%n(i,k)     = ice%n(i,k)     + mult_n
               ice%q(i,k)     = ice%q(i,k)     + mult_q
               graupel%q(i,k) = graupel%q(i,k) - mult_q
            ENDIF

            ! enhancement of melting of graupel
            IF (T_a > T_3 .AND. enhanced_melting) THEN
               melt_q = c_w / L_ew * (T_a - T_3) * rime_q
               melt_n = melt_q / x_g

               melt_q = MIN(graupel%q(i,k),melt_q)
               melt_n = MIN(graupel%n(i,k),melt_n)

               graupel%q(i,k) = graupel%q(i,k) - melt_q
               rain%q(i,k)    = rain%q(i,k)    + melt_q

               graupel%n(i,k) = graupel%n(i,k) - melt_n
               rain%n(i,k)    = rain%n(i,k)    + melt_n
            ELSE
               melt_q = 0.0_wp
            ENDIF

            ! shedding
            IF ((graupel_shedding .AND. D_g > D_shed_g .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
               q_g = graupel%q(i,k)
               n_g = graupel%n(i,k)
               x_g = graupel%meanmass(q_g,n_g)

               shed_q = MIN(q_g,rime_q)
               IF (T_a <= T_3) THEN
                  shed_n = shed_q / MIN(x_shed,x_g)
               ELSE
                  shed_n = shed_q / MAX(x_r,x_g)
               ENDIF

               graupel%q(i,k) = graupel%q(i,k) - shed_q
               rain%q(i,k)    = rain%q(i,k)    + shed_q
               rain%n(i,k)    = rain%n(i,k)    + shed_n
            ELSE
               shed_q = 0.0
            ENDIF

         ENDIF
      ENDDO
   ENDDO

  END SUBROUTINE graupel_rain_riming

  SUBROUTINE hail_rain_riming(ik_slice, dt, atmo, hail, rain, ice)
    !*******************************************************************************
    !                                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: rain, ice, hail

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h
    REAL(wp)            :: q_r,n_r,x_r,d_r,v_r
    REAL(wp)            :: rime_n,rime_q,melt_n,melt_q,shed_n,shed_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp), PARAMETER :: &
         x_shed = 4./3. * pi * rho_w * r_shedding**3, &
         const3 = 1.0/(T_mult_opt - T_mult_min), &
         const4 = 1.0/(T_mult_opt - T_mult_max)


    IF (isdebug) CALL message(routine, "hail_rain_riming")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_r = rain%q(i,k)
          q_h = hail%q(i,k)
          T_a = atmo%T(i,k)

          IF (q_r > q_crit .AND. q_h > q_crit) THEN
            n_r = rain%n(i,k)
            n_h = hail%n(i,k)

            x_h = hail%meanmass(q_h,n_h)
            D_h = hail%diameter(x_h)
            v_h = hail%velocity(x_h) * hail%rho_v(i,k)

            x_r = rain%meanmass(q_r,n_r)
            D_r = rain%diameter(x_r)
            v_r = rain%velocity(x_r) * rain%rho_v(i,k)

            rime_n = pi4 * n_h * n_r * dt &
                 & *     (hrr_params%delta_n_aa * D_h**2 + hrr_params%delta_n_ab * D_h*D_r + hrr_params%delta_n_bb * D_r**2) &
                 & * sqrt(hrr_params%theta_n_aa * v_h**2 - hrr_params%theta_n_ab * v_h*v_r + hrr_params%theta_n_bb * v_r**2)

            rime_q = pi4 * n_h * q_r * dt &
                 & *     (hrr_params%delta_q_aa * D_h**2 + hrr_params%delta_q_ab * D_h*D_r + hrr_params%delta_q_bb * D_r**2) &
                 & * sqrt(hrr_params%theta_q_aa * v_h**2 - hrr_params%theta_q_ab * v_h*v_r + hrr_params%theta_q_bb * v_r**2)

            rime_n = MIN(n_r,rime_n)
            rime_q = MIN(q_r,rime_q)

            hail%q(i,k) = hail%q(i,k) + rime_q
            rain%q(i,k) = rain%q(i,k) - rime_q
            rain%n(i,k) = rain%n(i,k) - rime_n

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

              ice%n(i,k)  = ice%n(i,k)  + mult_n
              ice%q(i,k)  = ice%q(i,k)  + mult_q
              hail%q(i,k) = hail%q(i,k) - mult_q
            ENDIF

            ! enhancement of melting of hail
            IF (T_a > T_3 .AND. enhanced_melting) THEN
              melt_q = c_w / L_ew * (T_a - T_3) * rime_q
              melt_n = melt_q / x_h

              melt_q = MIN(hail%q(i,k),melt_q)
              melt_n = MIN(hail%n(i,k),melt_n)

              hail%q(i,k) = hail%q(i,k) - melt_q
              rain%q(i,k) = rain%q(i,k) + melt_q

              hail%n(i,k) = hail%n(i,k) - melt_n
              rain%n(i,k) = rain%n(i,k) + melt_n
            ELSE
              melt_q = 0.0
            ENDIF

            ! shedding
            IF ((hail_shedding .AND. D_h > D_shed_h .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
              q_h = hail%q(i,k)
              n_h = hail%n(i,k)
              x_h = hail%meanmass(q_h,n_h)

              shed_q = MIN(q_h,rime_q)
              IF (T_a <= T_3) THEN
                shed_n = shed_q / MIN(x_shed,x_h)
              ELSE
                !  shed_n = shed_q / x_shed
                shed_n = shed_q / MAX(x_r,x_h)
              ENDIF

              hail%q(i,k) = hail%q(i,k) - shed_q
              rain%q(i,k) = rain%q(i,k)    + shed_q
              rain%n(i,k) = rain%n(i,k)    + shed_n
            ELSE
              shed_q = 0.0
            ENDIF

          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE hail_rain_riming

  SUBROUTINE graupel_melting(ik_slice, dt, atmo, graupel, rain)
    !*******************************************************************************
    ! Melting of graupel                                                           *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: rain, graupel
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_g,n_g,x_g,d_g,v_g,T_a,e_a
    REAL(wp)            :: melt,melt_v,melt_h,melt_n,melt_q
    REAL(wp)            :: fh_q,fv_q

    IF (isdebug) CALL message(routine, "graupel_melting")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)
          q_g = graupel%q(i,k)

          IF (T_a > T_3 .AND. q_g > 0.0) THEN
             e_a = e_ws(T_a)
             n_g = graupel%n(i,k)

             x_g = graupel%meanmass(q_g,n_g)
             D_g = graupel%diameter(x_g)
             v_g = graupel%velocity(x_g) * graupel%rho_v(i,k)

             fv_q = gm_params%a_vent + gm_params%b_vent * SQRT(v_g*D_g)
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

             graupel%q(i,k) = graupel%q(i,k) - melt_q
             rain%q(i,k)    = rain%q(i,k)    + melt_q

             graupel%n(i,k) = graupel%n(i,k) - melt_n
             rain%n(i,k)    = rain%n(i,k)    + melt_n

          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE graupel_melting

  SUBROUTINE hail_melting(ik_slice, dt, atmo, hail, rain)
    !*******************************************************************************
    ! Melting of hail                                                              *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: rain, hail
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_h,n_h,x_h,d_h,v_h,T_a,e_a
    REAL(wp)            :: melt,melt_v,melt_h,melt_n,melt_q
    REAL(wp)            :: fh_q,fv_q

    IF (isdebug) CALL message(routine, "hail_melting")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)
          q_h = hail%q(i,k)

          IF (T_a > T_3 .AND. q_h > 0.0_wp) THEN
            e_a = e_ws(T_a)
            n_h = hail%n(i,k)

            x_h = hail%meanmass(q_h,n_h)
            D_h = hail%diameter(x_h)
            v_h = hail%velocity(x_h) * hail%rho_v(i,k)

            fv_q = hm_params%a_vent + hm_params%b_vent * sqrt(v_h*D_h)
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

            hail%q(i,k) = hail%q(i,k) - melt_q
            rain%q(i,k) = rain%q(i,k) + melt_q

            hail%n(i,k) = hail%n(i,k) - melt_n
            rain%n(i,k) = rain%n(i,k) + melt_n

         ENDIF
      ENDDO
   ENDDO

  END SUBROUTINE hail_melting

  SUBROUTINE graupel_hail_conv_wet_gamlook(ik_slice, atmo, graupel, cloud, &
       rain, ice, snow, hail)
    !*******************************************************************************
    !  Wet growth and conversion of graupel to hail                                *
    !  (uses look-up table for incomplete gamma functions)                         *
    !  by Uli Blahak                                                               *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud, rain, ice, graupel, snow, hail
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a, p_a, d_trenn, qw_a, qi_a, N_0, lam, xmin
    REAL(wp)            :: q_g,n_g,x_g,d_g
    REAL(wp)            :: q_c,q_r
    REAL(wp)            :: conv_n,conv_q

    IF (isdebug) CALL message(routine, "graupel_hail_conv_wet_gamlook")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          q_c = cloud%q(i,k)
          q_r = rain%q(i,k)
          q_g = graupel%q(i,k)
          n_g = graupel%n(i,k)

          x_g = graupel%meanmass(q_g,n_g)
          D_g = graupel%diameter(x_g)
          n_g = q_g / x_g  ! for consistency for limiters, n_g is used explicitly below

          T_a = atmo%T(i,k)
          p_a = atmo%p(i,k)

          !..supercooled liquid water in the cloud environment = sum of rain and cloud water
          qw_a = q_r + q_c

          IF (T_a < T_3 .AND. q_g > q_crit_gc .AND. qw_a > 1e-3) THEN

            !.. Umgebungsgehalt Eispartikel (vernachl. werden Graupel und Hagel wg. geringer Kollisionseff.)
            !.. koennte problematisch sein, weil in konvekt. Wolken viel mehr Graupel und Hagel enthalten ist!!!
            qi_a = ice%q(i,k) + snow%q(i,k)
            d_trenn = dmin_wg_gr_ltab_equi(p_a,T_a,qw_a,qi_a,ltabdminwgg)

            IF (d_trenn > 0.0_wp .AND. d_trenn < 10.0_wp * D_g) THEN

               xmin = exp(log(d_trenn/graupel%a_geo)*(1.0_wp/graupel%b_geo))

               lam  = exp(log(graupel_g2/(graupel_g1*x_g))*(graupel%mu))
               xmin = exp(log(xmin)*graupel%mu)
               n_0  = graupel%mu * n_g * exp(log(lam)*graupel_nm1) / graupel_g1

               conv_n = n_0 / (graupel%mu * exp(log(lam)*graupel_nm1)) * incgfct_upper_lookup(lam*xmin,graupel_ltable1)
               conv_q = n_0 / (graupel%mu * exp(log(lam)*graupel_nm2)) * incgfct_upper_lookup(lam*xmin,graupel_ltable2)

               conv_n = MIN(conv_n,n_g)
               conv_q = MIN(conv_q,q_g)

               graupel%q(i,k) = q_g - conv_q
               graupel%n(i,k) = n_g - conv_n

               hail%q(i,k) = hail%q(i,k) + conv_q
               hail%n(i,k) = hail%n(i,k) + conv_n

            END IF
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE graupel_hail_conv_wet_gamlook

  SUBROUTINE ice_riming(ik_slice, dt, atmo, ice, cloud, rain, graupel, &
       dep_rate_ice)
    !*******************************************************************************
    !  Riming of ice with cloud droplet and rain drops. First the process rates    *
    !  are calculated in                                                           *
    !      snow_cloud_riming ()                                                    *
    !      snow_rain_riming ()                                                     *
    !  using those rates and the previously calculated and stored deposition       *
    !  rate the conversion of snow to graupel and rain is done.                    *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud, ice, rain, graupel
    REAL(wp), INTENT (IN), DIMENSION(:,:) :: dep_rate_ice

    REAL(wp), DIMENSION(size(dep_rate_ice,1),size(dep_rate_ice,2)) ::       &
         &               rime_rate_qc, rime_rate_nc,                        &
         &               rime_rate_qi, rime_rate_qr, rime_rate_nr

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: q_i,n_i,x_i,d_i
    REAL(wp)            :: T_a,x_coll_c,x_r
    REAL(wp)            :: rime_n,rime_q,rime_qr,rime_qi
    REAL(wp)            :: conv_n,conv_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp)            :: const2
    REAL(wp), PARAMETER :: &
         const3 = 1.0_wp/(T_mult_opt - T_mult_min), &
         const4 = 1.0_wp/(T_mult_opt - T_mult_max), &
         const5 = alpha_spacefilling * rho_w/rho_ice


    IF (isdebug) CALL message(routine, "ice riming")

    rime_rate_qc(:,:) = 0.0_wp
    rime_rate_qr(:,:) = 0.0_wp
    rime_rate_qi(:,:) = 0.0_wp
    rime_rate_nc(:,:) = 0.0_wp
    rime_rate_qr(:,:) = 0.0_wp

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    CALL ice_cloud_riming()
    CALL ice_rain_riming()

    x_coll_c = (D_coll_c/cloud%a_geo)**3   !..lower mean mass for collection

    const2 = 1.0/x_coll_c
    !
    ! Complete ice-cloud and ice-rain riming

    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)

          IF (dep_rate_ice(i,k) > 0.0_wp &
               & .and. dep_rate_ice(i,k) .ge. rime_rate_qc(i,k)+rime_rate_qr(i,k)) THEN

            !
            ! 1) Depositional growth is stronger than riming growth, therefore ice stays ice
            !

            !.. ice_cloud_riming

            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(cloud%q(i,k),rime_q)
              rime_n = MIN(cloud%n(i,k),rime_n)

              ice%q(i,k)   = ice%q(i,k)  + rime_q
              cloud%q(i,k) = cloud%q(i,k) - rime_q
              cloud%n(i,k) = cloud%n(i,k) - rime_n

              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                ice%n(i,k) = ice%n(i,k)  + mult_n
              ENDIF
            END IF

            !.. ice_rain_riming
            IF (rime_rate_qr(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qr(i,k)
              rime_n = rime_rate_nr(i,k)
              rime_q = MIN(rain%q(i,k),rime_q)
              rime_n = MIN(rain%n(i,k),rime_n)

              ice%q(i,k)  = ice%q(i,k)  + rime_q
              rain%q(i,k) = rain%q(i,k) - rime_q
              rain%n(i,k) = rain%n(i,k) - rime_n

              !..ice multiplication
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                ice%n(i,k) = ice%n(i,k)  + mult_n
              ENDIF

            END IF

         ELSE

            !
            !.. 2) Depositional growth negative or smaller than riming growth, therefore ice is
            !      allowed to convert to graupel and / or hail
            !

            !.. ice_cloud_riming
            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              n_i = ice%n(i,k)
              q_i = ice%q(i,k)
              x_i = ice%meanmass(q_i,n_i)
              D_i = ice%diameter(x_i)

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(cloud%q(i,k),rime_q)
              rime_n = MIN(cloud%n(i,k),rime_n)

              ice%q(i,k)   = ice%q(i,k)   + rime_q
              cloud%q(i,k) = cloud%q(i,k) - rime_q
              cloud%n(i,k) = cloud%n(i,k) - rime_n

              ! ice multiplication
              mult_q = 0.0_wp
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.0_wp,MIN(mult_1,1.0_wp))
                mult_2 = MAX(0.0_wp,MIN(mult_2,1.0_wp))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                ice%n(i,k) = ice%n(i,k)  + mult_n
              ENDIF

              ! conversion ice -> graupel (depends on alpha_spacefilling)
              IF (D_i > D_conv_ig) THEN
                 q_i = ice%q(i,k)
                 conv_q = (rime_q - mult_q) / ( const5 * (pi6*rho_ice*d_i**3/x_i - 1.0) )
                 conv_q = MIN(q_i,conv_q)
                 x_i    = ice%meanmass(q_i,n_i)
                 conv_n = conv_q / MAX(x_i,x_conv)
                 conv_n = MIN(ice%n(i,k),conv_n)
              ELSE
                 conv_q = 0.0_wp
                 conv_n = 0.0_wp
              ENDIF

              ice%q(i,k)     = ice%q(i,k)     - conv_q
              graupel%q(i,k) = graupel%q(i,k) + conv_q
              ice%n(i,k)     = ice%n(i,k)     - conv_n
              graupel%n(i,k) = graupel%n(i,k) + conv_n
            END IF

            !.. ice_rain_riming
            IF (rime_rate_qi(i,k) > 0.0_wp) THEN

              rime_qi = rime_rate_qi(i,k)
              rime_qr = rime_rate_qr(i,k)
              rime_n  = rime_rate_nr(i,k)
              rime_n  = MIN(MIN(rain%n(i,k),ice%n(i,k)),rime_n)
              rime_qr = MIN(rain%q(i,k),rime_qr)
              rime_qi = MIN(ice%q(i,k),rime_qi)

              ice%n(i,k)  = ice%n(i,k)  - rime_n
              rain%n(i,k) = rain%n(i,k) - rime_n
              ice%q(i,k)  = ice%q(i,k)  - rime_qi
              rain%q(i,k) = rain%q(i,k) - rime_qr

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
                 ! i.e. undo time integration, but with modified rain%n
                 x_r = rain%meanmass(rain%q(i,k),rain%n(i,k))
                 ice%n(i,k)  = ice%n(i,k)  + rime_n
                 rain%n(i,k) = rain%n(i,k) + rime_qr / x_r
                 ice%q(i,k)  = ice%q(i,k)  + rime_qi
                 rain%q(i,k) = rain%q(i,k) + rime_qr
              ELSE
                 ! new ice particles from multiplication
                 ice%n(i,k) = ice%n(i,k) + mult_n
                 ice%q(i,k) = ice%q(i,k) + mult_q
                 ! riming to graupel
                 graupel%n(i,k) = graupel%n(i,k) + rime_n
                 graupel%q(i,k) = graupel%q(i,k) + rime_qi + rime_qr - mult_q
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
      REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
      REAL(wp)            :: q_c,n_c,x_c,d_c,v_c,e_coll,x_coll_c
      REAL(wp)            :: rime_n,rime_q
      REAL(wp), PARAMETER :: &
           !..collision efficiency coeff
           const1   = ecoll_ic/(D_coll_c - D_crit_c)

      IF (isdebug) CALL message(routine, "ice_cloud_riming")

      x_coll_c = (D_coll_c/cloud%a_geo)**3         !..lower mass threshold for collection

      DO k = kstart,kend
         DO i = istart,iend

            n_c = cloud%n(i,k)
            q_c = cloud%q(i,k)
            n_i = ice%n(i,k)
            q_i = ice%q(i,k)

            x_c = cloud%meanmass(q_c,n_c)
            D_c = cloud%diameter(x_c)
            x_i = ice%meanmass(q_i,n_i)
            D_i = ice%diameter(x_i)

            IF (q_c > q_crit_c .AND. q_i > q_crit_ic .AND. D_i > D_crit_ic .AND. D_c > D_crit_c) THEN

               v_c = cloud%velocity(x_c) * cloud%rho_v(i,k)
               v_i = ice%velocity(x_i)   * ice%rho_v(i,k)

               e_coll = MIN(ecoll_ic, MAX(const1*(D_c - D_crit_c), ecoll_min))

               rime_n = pi4 * e_coll * n_i * n_c * dt &
                    &   *     (icr_params%delta_n_aa * D_i * D_i &
                    &          + icr_params%delta_n_ab * D_i * D_c &
                    &          + icr_params%delta_n_bb * D_c * D_c) &
                    &   * SQRT(icr_params%theta_n_aa * v_i * v_i &
                    &          - icr_params%theta_n_ab * v_i * v_c &
                    &          + icr_params%theta_n_bb * v_c * v_c &
                    &          + ice_s_vel**2)

               rime_q = pi4 * e_coll * n_i * q_c * dt &
                    &   *     (icr_params%delta_q_aa * D_i * D_i &
                    &          + icr_params%delta_q_ab * D_i * D_c &
                    &          + icr_params%delta_q_bb * D_c * D_c) &
                    &   * SQRT(icr_params%theta_q_aa * v_i * v_i &
                    &          - icr_params%theta_q_ab * v_i * v_c &
                    &          + icr_params%theta_q_bb * v_c * v_c  &
                    &          + ice_s_vel**2)

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
      REAL(wp)            :: q_i,n_i,x_i,d_i,v_i
      REAL(wp)            :: q_r,n_r,x_r,d_r,v_r
      REAL(wp)            :: rime_n,rime_qi,rime_qr

      IF (isdebug) CALL message(routine, "ice_rain_riming")

      DO k = kstart,kend
         DO i = istart,iend

            q_r = rain%q(i,k)
            q_i = ice%q(i,k)
            n_r = rain%n(i,k)
            n_i = ice%n(i,k)

            x_i = ice%meanmass(q_i,n_i)
            D_i = ice%diameter(x_i)

            IF (q_r > q_crit .AND. q_i > q_crit_ir .AND. D_i > D_crit_ir) THEN

               x_r = rain%meanmass(q_r,n_r)
               D_r = rain%diameter(x_r)
               v_r = rain%velocity(x_r) * rain%rho_v(i,k)
               v_i = ice%velocity(x_i) * ice%rho_v(i,k)

               rime_n  = pi4 * n_i * n_r * dt &
                    &  *     (irr_params%delta_n_aa * D_i * D_i &
                    &         + irr_params%delta_n_ab * D_i * D_r &
                    &         + irr_params%delta_n_bb * D_r * D_r) &
                    &  * SQRT(irr_params%theta_n_aa * v_i * v_i &
                    &         - irr_params%theta_n_ab * v_i * v_r &
                    &         + irr_params%theta_n_bb * v_r * v_r &
                    &         + ice_s_vel**2)

               rime_qr = pi4 * n_i * q_r * dt &
                    &  *     (irr_params%delta_n_aa * D_i * D_i &
                    &         + irr_params%delta_q_ab * D_i * D_r &
                    &         + irr_params%delta_q_bb * D_r * D_r) &
                    &  * SQRT(irr_params%theta_n_aa * v_i * v_i &
                    &         - irr_params%theta_q_ab * v_i * v_r &
                    &         + irr_params%theta_q_bb * v_r * v_r &
                    &         + ice_s_vel**2)

               rime_qi = pi4 * n_r * q_i * dt &
                    &  *     (irr_params%delta_q_aa * D_i * D_i &
                    &         + irr_params%delta_q_ba * D_i * D_r &
                    &         + irr_params%delta_n_bb * D_r * D_r) &
                    &  * SQRT(irr_params%theta_q_aa * v_i * v_i &
                    &         - irr_params%theta_q_ba * v_i * v_r &
                    &         + irr_params%theta_n_bb * v_r * v_r &
                    &         + ice_s_vel**2)

               rime_rate_nr(i,k) = rime_n
               rime_rate_qi(i,k) = rime_qi
               rime_rate_qr(i,k) = rime_qr

            ENDIF
         ENDDO
      ENDDO

    END SUBROUTINE ice_rain_riming

  END SUBROUTINE ice_riming

  SUBROUTINE snow_riming(ik_slice, dt, atmo, &
       snow, cloud, rain, ice, graupel, dep_rate_snow)
    !*******************************************************************************
    !  Riming of snow with cloud droplet and rain drops. First the process rates   *
    !  are calculated in                                                           *
    !      snow_cloud_riming ()                                                    *
    !      snow_rain_riming ()                                                     *
    !  using those rates and the previously calculated and stored deposition       *
    !  rate the conversion of snow to graupel and rain is done.                    *
    !*******************************************************************************
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    ! time step within two-moment scheme
    REAL(wp), INTENT(in) :: dt
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud, rain, snow, ice, graupel
    REAL(wp), INTENT(IN), DIMENSION(:,:) :: dep_rate_snow

    REAL(wp), DIMENSION(size(dep_rate_snow,1),size(dep_rate_snow,2)) ::       &
         & rime_rate_qc, rime_rate_nc,                                        &
         & rime_rate_qs, rime_rate_qr, rime_rate_nr

    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER             :: i,k
    REAL(wp)            :: T_a
    REAL(wp)            :: q_s,n_s,x_s,d_s
    REAL(wp)            :: x_coll_c,x_r
    REAL(wp)            :: rime_n,rime_q,rime_qr,rime_qs
    REAL(wp)            :: conv_n,conv_q
    REAL(wp)            :: mult_n,mult_q,mult_1,mult_2
    REAL(wp)            :: const2
    REAL(wp), PARAMETER :: &
         const3 = 1.0/(T_mult_opt - T_mult_min), &
         const4 = 1.0/(T_mult_opt - T_mult_max), &
         const5 = alpha_spacefilling * rho_w/rho_ice

    IF (isdebug) CALL message(routine, "snow_riming")

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    rime_rate_qc(:,:) = 0.0_wp
    rime_rate_qr(:,:) = 0.0_wp
    rime_rate_qs(:,:) = 0.0_wp
    rime_rate_nc(:,:) = 0.0_wp
    rime_rate_qr(:,:) = 0.0_wp

    CALL snow_cloud_riming()
    CALL snow_rain_riming()

    x_coll_c = (D_coll_c/cloud%a_geo)**3   !..lower mean mass for collection

    const2 = 1.0/x_coll_c

    DO k = kstart,kend
       DO i = istart,iend

          T_a = atmo%T(i,k)

          IF (dep_rate_snow(i,k) > 0.0_wp &
               & .AND. dep_rate_snow(i,k) .ge. rime_rate_qc(i,k) + rime_rate_qr(i,k)) THEN

            !
            ! 1) Depositional growth is stronger than riming growth, therefore snow stays snow:
            !

            !.. time integration of snow_cloud_riming

            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(cloud%q(i,k),rime_q)
              rime_n = MIN(cloud%n(i,k),rime_n)

              snow%q(i,k)  = snow%q(i,k)  + rime_q
              cloud%q(i,k) = cloud%q(i,k) - rime_q
              cloud%n(i,k) = cloud%n(i,k) - rime_n

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

                ice%n(i,k)  = ice%n(i,k)  + mult_n
                ice%q(i,k)  = ice%q(i,k)  + mult_q
                snow%q(i,k) = snow%q(i,k) - mult_q
              ENDIF

            END IF

            !.. time integration snow_rain_riming

            IF (rime_rate_qr(i,k) > 0.0_wp) THEN

              rime_q = rime_rate_qr(i,k)
              rime_n = rime_rate_nr(i,k)
              rime_q = MIN(rain%q(i,k),rime_q)
              rime_n = MIN(rain%n(i,k),rime_n)

              snow%q(i,k) = snow%q(i,k) + rime_q
              rain%q(i,k) = rain%q(i,k) - rime_q
              rain%n(i,k) = rain%n(i,k) - rime_n

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

                ice%n(i,k)  = ice%n(i,k)  + mult_n
                ice%q(i,k)  = ice%q(i,k)  + mult_q
                snow%q(i,k) = snow%q(i,k) - mult_q
              ENDIF
            END IF

          ELSE

            !
            !.. 2) Depositional growth is negative or smaller than riming growth, therefore snow is
            !      allowed to convert to graupel and / or hail:
            !

            !.. time integration of snow_cloud_riming

            IF (rime_rate_qc(i,k) > 0.0_wp) THEN

              n_s = snow%n(i,k)
              q_s = snow%q(i,k)
              x_s = snow%meanmass(q_s,n_s)
              D_s = snow%diameter(x_s)

              rime_q = rime_rate_qc(i,k)
              rime_n = rime_rate_nc(i,k)
              rime_q = MIN(cloud%q(i,k),rime_q)
              rime_n = MIN(cloud%n(i,k),rime_n)

              snow%q(i,k)  = snow%q(i,k)  + rime_q
              cloud%q(i,k) = cloud%q(i,k) - rime_q
              cloud%n(i,k) = cloud%n(i,k) - rime_n

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

                ice%n(i,k)  = ice%n(i,k)  + mult_n
                ice%q(i,k)  = ice%q(i,k)  + mult_q
                snow%q(i,k) = snow%q(i,k) - mult_q
              ENDIF

              !.. conversion of snow to graupel, depends on alpha_spacefilling

              IF (D_s > D_conv_sg) THEN
                 conv_q = (rime_q - mult_q) / ( const5 * (pi6*rho_ice*d_s**3/x_s - 1.0) )
                 conv_q = MIN(q_s,conv_q)
                 x_s    = snow%meanmass(q_s,n_s)
                 conv_n = conv_q / MAX(x_s,x_conv)
                 conv_n = MIN(snow%n(i,k),conv_n)
              ELSE
                 conv_q = 0.0_wp
                 conv_n = 0.0_wp
              ENDIF

              snow%q(i,k)    = snow%q(i,k)    - conv_q
              graupel%q(i,k) = graupel%q(i,k) + conv_q

              snow%n(i,k)    = snow%n(i,k)    - conv_n
              graupel%n(i,k) = graupel%n(i,k) + conv_n

            END IF

            !.. time integration of snow_rain_riming

            IF (rime_rate_qs(i,k) > 0.0_wp) THEN

              rime_qs = rime_rate_qs(i,k)
              rime_qr = rime_rate_qr(i,k)
              rime_n  = rime_rate_nr(i,k)
              rime_qr = MIN(rain%q(i,k),rime_qr)
              rime_qs = MIN(snow%q(i,k),rime_qs)
              rime_n  = MIN(rain%n(i,k),rime_n)
              rime_n  = MIN(snow%n(i,k),rime_n)

              snow%n(i,k) = snow%n(i,k) - rime_n
              rain%n(i,k) = rain%n(i,k) - rime_n
              snow%q(i,k) = snow%q(i,k) - rime_qs
              rain%q(i,k) = rain%q(i,k) - rime_qr

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
                 ! i.e. undo time integration, but with modified rain%n
                 x_r = rain%meanmass(rain%q(i,k),rain%n(i,k))
                 snow%n(i,k) = snow%n(i,k) + rime_n
                 rain%n(i,k) = rain%n(i,k) + rime_qr / x_r
                 snow%q(i,k) = snow%q(i,k) + rime_qs
                 rain%q(i,k) = rain%q(i,k) + rime_qr
              ELSE
                 ! new ice particles from multiplication
                 ice%n(i,k)  = ice%n(i,k)  + mult_n
                 ice%q(i,k)  = ice%q(i,k)  + mult_q
                 ! riming to graupel
                 graupel%n(i,k) = graupel%n(i,k) + rime_n
                 graupel%q(i,k) = graupel%q(i,k) + rime_qr + rime_qs - mult_q
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
     REAL(wp)            :: q_s,n_s,x_s,d_s,v_s
     REAL(wp)            :: q_r,n_r,x_r,d_r,v_r
     REAL(wp)            :: rime_n,rime_qr,rime_qs

     IF (isdebug) CALL message(routine, "snow_rain_riming")

     DO k = kstart,kend
       DO i = istart,iend

         q_r = rain%q(i,k)
         q_s = snow%q(i,k)
         n_r = rain%n(i,k)
         n_s = snow%n(i,k)

         x_s = snow%meanmass(q_s,n_s)
         D_s = snow%diameter(x_s)

         IF (q_r > q_crit .AND. q_s > q_crit_sr .AND. D_s > D_crit_sr) THEN

            x_r = rain%meanmass(q_r,n_r)
            D_r = rain%diameter(x_r)
            v_r = rain%velocity(x_r) * rain%rho_v(i,k)
            v_s = snow%velocity(x_s) * snow%rho_v(i,k)

            rime_n = pi4 * n_s * n_r * dt &
                 & *     (srr_params%delta_n_aa * D_s**2 + srr_params%delta_n_ab * D_s*D_r + srr_params%delta_n_bb * D_r**2) &
                 & * sqrt(srr_params%theta_n_aa * v_s**2 - srr_params%theta_n_ab * v_s*v_r + srr_params%theta_n_bb * v_r**2  &
                 &       +snow_s_vel**2)

            rime_qr = pi4 * n_s * q_r * dt &
                 &  *     (srr_params%delta_n_aa * D_s**2 + srr_params%delta_q_ab * D_s*D_r + srr_params%delta_q_bb * D_r**2) &
                 &  * sqrt(srr_params%theta_n_aa * v_s**2 - srr_params%theta_q_ab * v_s*v_r + srr_params%theta_q_bb * v_r**2  &
                 &        +snow_s_vel**2)

            rime_qs = pi4 * n_r * q_s * dt &
                 &  *     (srr_params%delta_q_aa * D_s**2 + srr_params%delta_q_ba * D_s*D_r + srr_params%delta_n_bb * D_r**2) &
                 &  * sqrt(srr_params%theta_q_aa * v_s**2 - srr_params%theta_q_ba * v_s*v_r + srr_params%theta_n_bb * v_r**2  &
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
     REAL(wp)            :: q_s,n_s,x_s,d_s,v_s
     REAL(wp)            :: q_c,n_c,x_c,d_c,v_c,e_coll
     REAL(wp)            :: rime_n,rime_q
     REAL(wp), PARAMETER :: const1 = ecoll_sc/(D_coll_c - D_crit_c)

     IF (isdebug) CALL message(routine, "snow_cloud_riming")

     DO k = kstart,kend
        DO i = istart,iend

           n_c = cloud%n(i,k)
           q_c = cloud%q(i,k)
           n_s = snow%n(i,k)
           q_s = snow%q(i,k)

           x_c = cloud%meanmass(q_c,n_c)
           D_c = cloud%diameter(x_c)
           x_s = snow%meanmass(q_s,n_s)
           D_s = snow%diameter(x_s)

           IF (q_c > q_crit_c .AND. q_s > q_crit_sc .AND. D_s > D_crit_sc .AND. D_c > D_crit_c) THEN

              v_c = cloud%velocity(x_c) * cloud%rho_v(i,k)
              v_s = snow%velocity(x_s)  * snow%rho_v(i,k)

              e_coll = MIN(ecoll_sc, MAX(const1*(D_c - D_crit_c), ecoll_min))

              rime_n = pi4 * e_coll * n_s * n_c * dt &
                   & *     (scr_params%delta_n_aa * D_s**2 + scr_params%delta_n_ab * D_s*D_c + scr_params%delta_n_bb * D_c**2) &
                   & * sqrt(scr_params%theta_n_aa * v_s**2 - scr_params%theta_n_ab * v_s*v_c + scr_params%theta_n_bb * v_c**2  &
                   &       +snow_s_vel**2)

              rime_q = pi4 * e_coll * n_s * q_c * dt &
                   & *     (scr_params%delta_q_aa * D_s**2 + scr_params%delta_q_ab * D_s*D_c + scr_params%delta_q_bb * D_c**2) &
                   & * sqrt(scr_params%theta_q_aa * v_s**2 - scr_params%theta_q_ab * v_s*v_c + scr_params%theta_q_bb * v_c**2  &
                   &       +snow_s_vel**2)

              rime_rate_qc(i,k) = rime_q
              rime_rate_nc(i,k) = rime_n

           ENDIF
        ENDDO
     ENDDO

   END SUBROUTINE snow_cloud_riming

 END SUBROUTINE snow_riming

 SUBROUTINE ccn_activation_sk(ik_slice, atmo,cloud,n_cn)
   !*******************************************************************************
   !       Calculation of ccn activation                                          *
   !       using the look-up tables by Segal and Khain 2006 (JGR, vol.11)         *
   !       (implemented by Heike Noppel)                                          *
   !*******************************************************************************

    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud
    REAL(wp), DIMENSION(:,:) :: n_cn

    ! Locale Variablen
    INTEGER, PARAMETER :: n_ncn=8, n_r2=3, n_lsigs=5, n_wcb=4
    INTEGER            :: i_lsigs, i_R2
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER            :: i,k,nuc_typ
    REAL(wp)           :: n_c,q_c,rho_a
    REAL(wp)           :: nuc_n, nuc_q
    REAL(wp)           :: Ncn, wcb, &
         tab_Ndrop(n_wcb,n_ncn), & ! number of cloud droplets in look_up-Table
         tab_Ndrop_i(n_wcb)        ! number of cloud droplets in look_up-Table,
                                   ! interpolated with respect to Ncn

    REAL(wp), PARAMETER :: &
          ! look-up-table for Ncn
         tab_Ncn(n_ncn) = (/   50.d06,  100.d06,  200.d06,  400.d06, &
         &                    800.d06, 1600.d06, 3200.d06, 6400.d06 /), &
         ! look-up_tbale for R2, in 10^(-6) m
         tab_R2(n_r2) = (/0.02d0, 0.03d0, 0.04d0/), &
         ! look-up-table for log(sigma_s)
         tab_lsigs(n_lsigs)  = (/0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0/), &
         ! look-up-table for w at cloud base
         tab_wcb(n_wcb) = (/0.5d0, 1.0d0, 2.5d0, 5.0d0/)

    nuc_typ = nuc_c_typ

    ! ATTENTION: At the moment only the values given above can be chosen for R2,
    ! and lsigs (see below).
    ! Only wcb and N_cn0 can be interpolated. For the others this possibility is still missing.
    ! For high N_cn0 and small values of wcb a kind of "saturation" for cloud droplet
    ! nucleation was assumed, because the original look-up tables don't give these values.
    ! This assumption might be wrong!!! (comment by Heike Noppel)

    IF(isdebug) THEN
       WRITE(txt,*) "cloud_activation_SK: nuc_typ = ",nuc_typ
       CALL message(routine, TRIM(txt))
    ENDIF

    ! this could be done in an IF(firstcall.eq.false) block and i-value stored in ccn_coeffs
    i_lsigs = 0
    i_R2   = 0
    DO k=1,n_lsigs
      IF (ccn_coeffs%lsigs==tab_lsigs(k)) THEN
         i_lsigs=k
         EXIT
      ENDIF
    END DO
    DO k=1,n_r2
      IF (ccn_coeffs%R2==tab_R2(k)) THEN
         i_R2=k
         EXIT
      ENDIF
    END DO
    IF (i_lsigs==0) THEN
       CALL finish(TRIM(routine),'Error in two_moment_mcrph: Invalid value for LSIGS in ccn_activation_sk')
    END IF
    IF (i_R2==0) THEN
       CALL finish(TRIM(routine),'Error in two_moment_mcrph: Invalid value for R2 in ccn_activation_sk')
    END IF

    CALL lookuptable(tab_Ndrop,i_lsigs,i_R2)

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       DO i = istart,iend

          nuc_q = 0.0d0
          nuc_n = 0.d0
          n_c   = cloud%n(i,k)
          q_c   = cloud%q(i,k)
          rho_a = atmo%rho(i,k)
          wcb   = atmo%w(i,k)

          if (q_c > 0.0_wp .and. wcb > 0.0_wp) then

             ! hard upper limit for number conc that
             ! eliminates also unrealistic high value
             ! that would come from the dynamical core

             cloud%n(i,k) = MIN(cloud%n(i,k),ccn_coeffs%Ncn0)

             Ncn = n_cn(i,k) ! number of CN from prognostic variable

             ! min value for vertical velocity (instead of Nmin of older code)
             wcb = max(wcb,0.2_wp)

             ! Interpolation of the look-up tables with respect to Ncn

             ! If Ncn is outside the range of the lookup table values, resulting
             ! NCCN are clipped to the margin values. For the case of these margin values
             ! beeing larger than Ncn, limit NCCN by Ncn:
             tab_Ndrop_i = MIN(ip_ndrop_ncn(tab_Ndrop,tab_Ncn,Ncn), Ncn)

             ! interpol. with respect to wcb = ip_ndrop_wcb
             nuc_n = MAX(ccn_coeffs%etas * ip_ndrop_wcb(tab_Ndrop_i,tab_wcb,wcb),ccn_coeffs%Nmin) - n_c

             nuc_n = MAX(nuc_n,0.0d0)

             nuc_q = MIN(nuc_n * cloud%x_min, atmo%qv(i,k))
             nuc_n = nuc_q / cloud%x_min

             n_cn(i,k)    = n_cn(i,k)    - MIN(Ncn,nuc_n)
             cloud%n(i,k) = cloud%n(i,k) + nuc_n
             cloud%q(i,k) = cloud%q(i,k) + nuc_q
             atmo%qv(i,k) = atmo%qv(i,k) - nuc_q

          END IF

       END DO
    END DO

  CONTAINS

    SUBROUTINE lookuptable(tab_ndrop,i_lsigs,i_R2)
      REAL(wp), INTENT(out) :: tab_ndrop(n_wcb,n_ncn)
      INTEGER, INTENT(in) :: i_lsigs, i_R2
      INTEGER :: i
      REAL(wp), PARAMETER :: ndrop(n_ncn,5,n_wcb,3) &
           ! look up tables
           ! Ncn              50       100       200       400       800       1600      3200      6400
           ! table4a (R2=0.02mum, wcb=0.5m/s) (for Ncn=3200  and Ncn=6400 "extrapolated")
           = RESHAPE((/ &
           !ndrop1_11 =
           42.2d06,  70.2d06, 112.2d06, 173.1d06, 263.7d06, 397.5d06, 397.5d06, 397.5d06, &
           !ndrop1_12 =
           35.5d06,  60.1d06, 100.0d06, 163.9d06, 264.5d06, 418.4d06, 418.4d06, 418.4d06, &
           !ndrop1_13
           32.6d06,  56.3d06,  96.7d06, 163.9d06, 272.0d06, 438.5d06, 438.5d06, 438.5d06, &
           !ndrop1_14
           30.9d06,  54.4d06,  94.6d06, 162.4d06, 271.9d06, 433.5d06, 433.5d06, 433.5d06, &
           !ndrop1_15
           29.4d06,  51.9d06,  89.9d06, 150.6d06, 236.5d06, 364.4d06, 364.4d06, 364.4d06, &
           ! table4b (R2=0.02mum, wcb=1.0m/s) (for Ncn=50 "interpolted" and Ncn=6400 extrapolated)
           !ndrop1_21
           45.3d06,  91.5d06, 158.7d06, 264.4d06, 423.1d06, 672.5d06, 397.5d06, 397.5d06, &
           !ndrop1_22
           38.5d06,  77.1d06, 133.0d06, 224.9d06, 376.5d06, 615.7d06, 418.4d06, 418.4d06, &
           !ndrop1_23
           35.0d06,  70.0d06, 122.5d06, 212.0d06, 362.1d06, 605.3d06, 438.5d06, 438.5d06, &
           !ndrop1_24
           32.4d06,  65.8d06, 116.4d06, 204.0d06, 350.6d06, 584.4d06, 433.5d06, 433.5d06, &
           !ndrop1_25
           31.2d06,  62.3d06, 110.1d06, 191.3d06, 320.6d06, 501.3d06, 364.4d06, 364.4d06, &
           ! table4c (R2=0.02mum, wcb=2.5m/s) (for Ncn=50 and Ncn=100 "interpolated")
           !ndrop1_31
           50.3d06, 100.5d06, 201.1d06, 373.1d06, 664.7d06,1132.8d06,1876.8d06,2973.7d06, &
           !ndrop1_32
           44.1d06,  88.1d06, 176.2d06, 314.0d06, 546.9d06, 941.4d06,1579.2d06,2542.2d06, &
           !ndrop1_33
           39.7d06,  79.5d06, 158.9d06, 283.4d06, 498.9d06, 865.9d06,1462.6d06,2355.8d06, &
           !ndrop1_34
           37.0d06,  74.0d06, 148.0d06, 264.6d06, 468.3d06, 813.3d06,1371.3d06,2137.2d06, &
           !ndrop1_35
           34.7d06,  69.4d06, 138.8d06, 246.9d06, 432.9d06, 737.8d06,1176.7d06,1733.0d06, &
           ! table4d (R2=0.02mum, wcb=5.0m/s) (for Ncn=50,100,200 "interpolated")
           !ndrop1_41
           51.5d06, 103.1d06, 206.1d06, 412.2d06, 788.1d06,1453.1d06,2585.1d06,4382.5d06, &
           !ndrop1_42
           46.6d06,  93.2d06, 186.3d06, 372.6d06, 657.2d06,1202.8d06,2098.0d06,3556.9d06, &
           !ndrop1_43
           70.0d06,  70.0d06, 168.8d06, 337.6d06, 606.7d06,1078.5d06,1889.0d06,3206.9d06, &
           !ndrop1_44
           42.2d06,  84.4d06, 166.4d06, 312.7d06, 562.2d06,1000.3d06,1741.1d06,2910.1d06, &
           !ndrop1_45
           36.5d06,  72.9d06, 145.8d06, 291.6d06, 521.0d06, 961.1d06,1551.1d06,2444.6d06, &
           ! table5a (R2=0.03mum, wcb=0.5m/s)
           !ndrop2_11
           50.0d06,  95.8d06, 176.2d06, 321.6d06, 562.3d06, 835.5d06, 835.5d06, 835.5d06, &
           !ndrop2_12
           44.7d06,  81.4d06, 144.5d06, 251.5d06, 422.7d06, 677.8d06, 677.8d06, 677.8d06, &
           !ndrop2_13
           40.2d06,  72.8d06, 129.3d06, 225.9d06, 379.9d06, 606.5d06, 606.5d06, 606.5d06, &
           !ndrop2_14
           37.2d06,  67.1d06, 119.5d06, 206.7d06, 340.5d06, 549.4d06, 549.4d06, 549.4d06, &
           !ndrop2_15
           33.6d06,  59.0d06,  99.4d06, 150.3d06, 251.8d06, 466.0d06, 466.0d06, 466.0d06, &
           ! table5b (R2=0.03mum, wcb=1.0m/s) (Ncn=50 "interpolated", Ncn=6400 "extrapolated)
           !ndrop2_21
           50.7d06, 101.4d06, 197.6d06, 357.2d06, 686.6d06,1186.4d06,1892.2d06,1892.2d06, &
           !ndrop2_22
           46.6d06,  93.3d06, 172.2d06, 312.1d06, 550.7d06, 931.6d06,1476.6d06,1476.6d06, &
           !ndrop2_23
           42.2d06,  84.4d06, 154.0d06, 276.3d06, 485.6d06, 811.2d06,1271.7d06,1271.7d06, &
           !ndrop2_24
           39.0d06,  77.9d06, 141.2d06, 251.8d06, 436.7d06, 708.7d06,1117.7d06,1117.7d06, &
           !ndrop2_25
           35.0d06,  70.1d06, 123.9d06, 210.2d06, 329.9d06, 511.9d06, 933.4d06, 933.4d06, &
           ! table5c (R2=0.03mum, wcb=2.5m/s)
           !ndrop2_31
           51.5d06, 103.0d06, 205.9d06, 406.3d06, 796.4d06,1524.0d06,2781.4d06,4609.3d06, &
           !ndrop2_32
           49.6d06,  99.1d06, 198.2d06, 375.5d06, 698.3d06,1264.1d06,2202.8d06,3503.6d06, &
           !ndrop2_33
           45.8d06,  91.6d06, 183.2d06, 339.5d06, 618.9d06,1105.2d06,1881.8d06,2930.9d06, &
           !ndrop2_34
           42.3d06,  84.7d06, 169.3d06, 310.3d06, 559.5d06, 981.7d06,1611.6d06,2455.6d06, &
           !ndrop2_35
           38.2d06,  76.4d06, 152.8d06, 237.3d06, 473.3d06, 773.1d06,1167.9d06,1935.0d06, &
           ! table5d (R2=0.03mum, wcb=5.0m/s)
           !ndrop2_41
           51.9d06, 103.8d06, 207.6d06, 415.1d06, 819.6d06,1616.4d06,3148.2d06,5787.9d06, &
           !ndrop2_42
           50.7d06, 101.5d06, 203.0d06, 405.9d06, 777.0d06,1463.8d06,2682.6d06,4683.0d06, &
           !ndrop2_43
           47.4d06,  94.9d06, 189.7d06, 379.4d06, 708.7d06,1301.3d06,2334.3d06,3951.8d06, &
           !ndrop2_44
           44.0d06,  88.1d06, 176.2d06, 352.3d06, 647.8d06,1173.0d06,2049.7d06,3315.6d06, &
           !ndrop2_45
           39.7d06,  79.4d06, 158.8d06, 317.6d06, 569.5d06, 988.5d06,1615.6d06,2430.3d06, &
           ! table6a (R2=0.04mum, wcb=0.5m/s)
           !ndrop3_11
           50.6d06, 100.3d06, 196.5d06, 374.7d06, 677.3d06,1138.9d06,1138.9d06,1138.9d06, &
           !ndrop3_12
           48.4d06,  91.9d06, 170.6d06, 306.9d06, 529.2d06, 862.4d06, 862.4d06, 862.4d06, &
           !ndrop3_13
           44.4d06,  82.5d06, 150.3d06, 266.4d06, 448.0d06, 740.7d06, 740.7d06, 740.7d06, &
           !ndrop3_14
           40.9d06,  75.0d06, 134.7d06, 231.9d06, 382.1d06, 657.6d06, 657.6d06, 657.6d06, &
           !ndrop3_15
           34.7d06,  59.3d06,  93.5d06, 156.8d06, 301.9d06, 603.8d06, 603.8d06, 603.8d06, &
           ! table6b (R2=0.04mum, wcb=1.0m/s)
           !ndrop3_21
           50.9d06, 101.7d06, 201.8d06, 398.8d06, 773.7d06,1420.8d06,2411.8d06,2411.8d06, &
           !ndrop3_22
           49.4d06,  98.9d06, 189.7d06, 356.2d06, 649.5d06,1117.9d06,1805.2d06,1805.2d06, &
           !ndrop3_23
           45.6d06,  91.8d06, 171.5d06, 214.9d06, 559.0d06, 932.8d06,1501.6d06,1501.6d06, &
           !ndrop3_24
           42.4d06,  84.7d06, 155.8d06, 280.5d06, 481.9d06, 779.0d06,1321.9d06,1321.9d06, &
           !ndrop3_25
           36.1d06,  72.1d06, 124.4d06, 198.4d06, 319.1d06, 603.8d06,1207.6d06,1207.6d06, &
           ! table6c (R2=0.04mum, wcb=2.5m/s)
           !ndrop3_31
           51.4d06, 102.8d06, 205.7d06, 406.9d06, 807.6d06,1597.5d06,3072.2d06,5393.9d06, &
           !ndrop3_32
           50.8d06, 101.8d06, 203.6d06, 396.0d06, 760.4d06,1422.1d06,2517.4d06,4062.8d06, &
           !ndrop3_33
           48.2d06,  96.4d06, 193.8d06, 367.3d06, 684.0d06,1238.3d06,2087.3d06,3287.1d06, &
           !ndrop3_34
           45.2d06,  90.4d06, 180.8d06, 335.7d06, 611.2d06,1066.3d06,1713.4d06,2780.3d06, &
           !ndrop3_35
           38.9d06,  77.8d06, 155.5d06, 273.7d06, 455.2d06, 702.2d06,1230.7d06,2453.7d06, &
           ! table6d (R2=0.04mum, wcb=5.0m/s)
           !ndrop3_41
           53.1d06, 106.2d06, 212.3d06, 414.6d06, 818.3d06,1622.2d06,3216.8d06,6243.9d06, &
           !ndrop3_42
           51.6d06, 103.2d06, 206.3d06, 412.5d06, 805.3d06,1557.4d06,2940.4d06,5210.1d06, &
           !ndrop3_43
           49.6d06,  99.2d06, 198.4d06, 396.7d06, 755.5d06,1414.5d06,2565.3d06,4288.1d06, &
           !ndrop3_44
           46.5d06,  93.0d06, 186.0d06, 371.9d06, 692.9d06,1262.0d06,2188.3d06,3461.2d06, &
           !ndrop3_45
           39.9d06,  79.9d06, 159.7d06, 319.4d06, 561.7d06, 953.9d06,1493.9d06,2464.7d06 /), (/ n_ncn, 5, n_wcb, 3 /) )
      IF (i_lsigs < 1 .OR. i_lsigs > 5 .OR. i_r2 < 1 .OR. i_r2 > 3) THEN
        WRITE(0, '(a)') "!!!! wrong value for lsigs or R2 in cloud_nucleation_SK !!!!!"
      END IF
      DO i = 1, 4
        tab_ndrop(i, :) = ndrop(:, i_lsigs, i, i_r2)
      END DO
    END SUBROUTINE lookuptable

    !--------------

    FUNCTION ip_ndrop_ncn(tab_ndrop,tab_ncn,Ncn)

      ! Interpolation of the look-up table with respect to aerosol concentration Ncn
      REAL(wp) :: ip_ndrop_ncn(n_wcb)
      REAL(wp), INTENT(in) :: Ncn, tab_ndrop(n_wcb,n_ncn), tab_ncn(n_ncn)
      INTEGER            ::  ki
      LOGICAL :: found

      ! Interpolation of Ndrop from the values of Ncn in the look-up-table to given Ncn
      found = .FALSE.
      IF (Ncn <= tab_ncn(1)) THEN
        ip_ndrop_ncn = tab_ndrop(:,1)
        found = .TRUE.
      ELSE IF (Ncn >= tab_ncn(n_ncn)) then
        ip_ndrop_ncn = tab_ndrop(:,n_ncn)
        found = .TRUE.
      ELSE
        DO ki = 1,n_ncn-1
          IF (Ncn >= tab_ncn(ki) .AND.  Ncn <= tab_ncn(ki+1)) THEN
            ip_ndrop_ncn(:) = tab_ndrop(:,ki) + &
                 (tab_ndrop(:,ki+1)-tab_ndrop(:,ki)) / (tab_ncn(ki+1)-tab_ncn(ki)) * (Ncn-tab_ncn(ki))
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF (.NOT.found) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, ip_ndrop_ncn: lookup table interpolation failed! ')
      END IF

      RETURN

    END FUNCTION  ip_ndrop_ncn

    !------------------
    ! Interpolation of the interpolated look-up table with respect to w_cb
    FUNCTION ip_ndrop_wcb(tab_Ndrop_i,tab_wcb,wcb)
      REAL(wp) :: ip_ndrop_wcb
      REAL(wp), INTENT(in) :: tab_wcb(1:n_wcb), tab_Ndrop_i(1:n_wcb), wcb
      INTEGER            :: ki
      LOGICAL :: found

      !... Interpolation for Ndrop from the values of wcb in the look-up-tabel to detected wcb

      found = .FALSE.
      ip_ndrop_wcb = 0
      IF (wcb <= tab_wcb(1)) THEN
        ip_ndrop_wcb = tab_ndrop_i(1) / tab_wcb(1) * wcb
        found = .TRUE.
      ELSE IF (wcb  >= tab_wcb(n_wcb)) then
        ip_ndrop_wcb = tab_ndrop_i(n_wcb)
        found = .TRUE.
      ELSE
        DO ki = 1,n_wcb-1
          IF (wcb >= tab_wcb(ki) .AND.  wcb <= tab_wcb(ki+1)) THEN
            ip_ndrop_wcb = tab_ndrop_i(ki) + &
              (tab_ndrop_i(ki+1)-tab_ndrop_i(ki)) / (tab_wcb(ki+1)-tab_wcb(ki)) * (wcb-tab_wcb(ki))
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF (.NOT.found) THEN
         CALL finish(TRIM(routine),'Error in two_moment_mcrph, ip_ndrop_wcb: lookup table interpolation failed! ')
      END IF

      RETURN

    END FUNCTION  ip_ndrop_wcb

  END SUBROUTINE ccn_activation_sk

 SUBROUTINE ccn_activation_hdcp2(ik_slice, atmo, cloud)
   !*******************************************************************************
   !       Calculation of ccn activation                                          *
   !       using the approach of Hande et al 2015                                 *
   !*******************************************************************************
    IMPLICIT NONE
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)

    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout)  :: cloud

    ! Locale Variablen
    ! start and end indices for 2D slices
    INTEGER :: istart, iend, kstart, kend
    INTEGER            :: i,k,nuc_typ
    REAL(wp)           :: n_c,q_c
    REAL(wp)           :: nuc_n,nuc_q
    REAL(wp)           :: wcb,pres

    REAL(wp)           :: acoeff,bcoeff,ccoeff,dcoeff
    REAL(wp)           :: a_ccn(4),b_ccn(4),c_ccn(4),d_ccn(4)

    ! Data from HDCP2_CCN_params.txt for 20130417
    DATA a_ccn(1),b_ccn(1),c_ccn(1),d_ccn(1) /183230691.161_wp,0.0001984051994_wp,16.2420263911_wp,287736034.13_wp/
    DATA a_ccn(2),b_ccn(2),c_ccn(2),d_ccn(2) /0.10147358938_wp,4.473190485e-05_wp,3.22011836758_wp,0.6258809883_wp/
    DATA a_ccn(3),b_ccn(3),c_ccn(3),d_ccn(3) /-0.2922395814_wp,0.0001843225275_wp,13.8499423719_wp,0.8907491812_wp/
    DATA a_ccn(4),b_ccn(4),c_ccn(4),d_ccn(4) /229189886.226_wp,0.0001986158191_wp,16.2461600644_wp,360848977.55_wp/

    istart = ik_slice(1)
    iend   = ik_slice(2)
    kstart = ik_slice(3)
    kend   = ik_slice(4)

    nuc_typ = nuc_c_typ

    IF(isdebug) THEN
       WRITE(txt,*) "cloud_activation_hdcp2: nuc_typ = ",nuc_typ ; CALL message(routine,TRIM(txt))
    ENDIF

    DO k = kstart,kend
       DO i = istart,iend

          nuc_q = 0.0d0
          nuc_n = 0.d0
          n_c   = cloud%n(i,k)
          q_c   = cloud%q(i,k)
          pres  = atmo%p(i,k)
          wcb   = atmo%w(i,k)

          if (q_c > 0.0_wp .and. wcb > 0.0_wp) then

             ! Based on write-up of Luke Hande of 6 May 2015

             acoeff = a_ccn(1) * atan(b_ccn(1) * pres - c_ccn(1)) + d_ccn(1)
             bcoeff = a_ccn(2) * atan(b_ccn(2) * pres - c_ccn(2)) + d_ccn(2)
             ccoeff = a_ccn(3) * atan(b_ccn(3) * pres - c_ccn(3)) + d_ccn(3)
             dcoeff = a_ccn(4) * atan(b_ccn(4) * pres - c_ccn(4)) + d_ccn(4)

             nuc_n = acoeff * atan(bcoeff * log(wcb) + ccoeff) + dcoeff

             nuc_n = MAX(MAX(nuc_n,10.0e-6_wp) - n_c,0.0_wp)

             nuc_q = MIN(nuc_n * cloud%x_min, atmo%qv(i,k))
             nuc_n = nuc_q / cloud%x_min

             cloud%n(i,k) = cloud%n(i,k) + nuc_n
             cloud%q(i,k) = cloud%q(i,k) + nuc_q
             atmo%qv(i,k) = atmo%qv(i,k) - nuc_q

          END IF

       END DO
    END DO

  END SUBROUTINE ccn_activation_hdcp2

  !*******************************************************************************
  ! Sedimentation subroutines for ICON
  !*******************************************************************************

  SUBROUTINE sedi_icon_rain (rain,qp,np,precrate,qc,rhocorr,adz,dt, &
      &                  its,ite,kts,kte,cmax) !

    TYPE(particle), INTENT(in)             :: rain
    INTEGER,  INTENT(IN)                   :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(INOUT):: qp,np,qc
    REAL(wp), DIMENSION(:,:), INTENT(IN)   :: adz,rhocorr
    REAL(wp), DIMENSION(:),   INTENT(OUT)  :: precrate
    REAL(wp), INTENT(IN)                   :: dt
    REAL(wp), INTENT(INOUT), OPTIONAL      :: cmax

    INTEGER  :: i, k, kk
    REAL(wp) :: x_p,D_m,D_p,mue,v_n,v_q, cmax_temp

    REAL(wp), DIMENSION(its:ite,kts-1:kte) :: q_fluss, n_fluss
    REAL(wp), DIMENSION(its:ite,kts-1:kte) :: v_n_sedi,v_q_sedi
    REAL(wp), DIMENSION(its:ite) :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    LOGICAL,  DIMENSION(its:ite) :: cflag

    v_n_sedi(:,kts-1) = 0.0_wp
    v_q_sedi(:,kts-1) = 0.0_wp
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
          D_p = D_m * EXP((-1./3.)*LOG((mue+3.)*(mue+2.)*(mue+1.)))

          v_n = rain_coeffs%alfa - rain_coeffs%beta * EXP(-(mue+1.)*LOG(1.0 + rain_coeffs%gama*D_p))
          v_q = rain_coeffs%alfa - rain_coeffs%beta * EXP(-(mue+4.)*LOG(1.0 + rain_coeffs%gama*D_p))
          v_n = v_n * rhocorr(i,k)
          v_q = v_q * rhocorr(i,k)

          v_n_sedi(i,k) = - v_n
          v_q_sedi(i,k) = - v_q
        ELSE
          v_n_sedi(i,k) = 0.0_wp
          v_q_sedi(i,k) = 0.0_wp
        ENDIF
      END DO
    END DO


    IF (PRESENT(cmax)) THEN
      cmax_temp = cmax
    ELSE
      cmax_temp = 0.0_wp
    END IF
    DO k = kts,kte

        DO i = its,ite
          v_nv(i) = 0.5 * (v_n_sedi(i,k-1)+v_n_sedi(i,k))   ! Das sollte wohl k-1 sein.
          v_qv(i) = 0.5 * (v_q_sedi(i,k-1)+v_q_sedi(i,k))
          ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
          c_nv(i) = -v_nv(i) * adz(i,k) * dt
          c_qv(i) = -v_qv(i) * adz(i,k) * dt
          cmax_temp = MAX(cmax_temp, c_qv(i))
        END DO

        kk = k
        s_nv = 0._wp
        DO i = its,ite
          IF (c_nv(i) <= 1._wp) THEN
            s_nv(i) = v_nv(i) * np(i,k)
          END IF
        END DO
        IF (ANY(c_nv > 1._wp)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_nv > 1._wp) .AND. kk > 2)
            DO i = its,ite
              IF (c_nv(i) > 1._wp) THEN
                cflag(i) = .TRUE.
                s_nv(i) = s_nv(i) + np(i,kk)/adz(i,kk)
                c_nv(i) = (c_nv(i) - 1._wp) * adz(i,kk-1)/adz(i,kk)
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
          IF (c_qv(i) <= 1._wp) THEN
            s_qv(i) = v_qv(i) * qp(i,k)
          END IF
        END DO
        IF (ANY(c_qv > 1._wp)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_qv > 1._wp) .AND. kk > 2)
            DO i = its,ite
              IF (c_qv(i) > 1._wp) THEN
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
    IF (PRESENT(cmax)) cmax = cmax_temp

    n_fluss(:,kts-1) = 0._wp ! obere Randbedingung
    q_fluss(:,kts-1) = 0._wp ! upper BC

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

    TYPE(particle), INTENT(in)              :: ptype
    TYPE(particle_sphere), INTENT(in)       :: pcoeffs
    INTEGER, INTENT(IN)                     :: its,ite,kts,kte
    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: qp,np
    REAL(wp), DIMENSION(:,:), INTENT(IN)    :: adz,rhocorr
    REAL(wp), DIMENSION(:),   INTENT(OUT)   :: precrate
    REAL(wp), INTENT(IN)                    :: dt
    REAL(wp), INTENT(INOUT), OPTIONAL       :: cmax

    INTEGER  :: i, k, kk
    REAL(wp) :: x_p,v_n,v_q,lam,cmax_temp

    REAL(wp), DIMENSION(its:ite,kts-1:kte+1) :: q_fluss, n_fluss
    REAL(wp), DIMENSION(its:ite,kts-1:kte+1) :: v_n_sedi,v_q_sedi
    REAL(wp), DIMENSION(its:ite) :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    LOGICAL,  DIMENSION(its:ite) :: cflag

    v_n_sedi(:, kts-1) = 0.0_wp
    v_q_sedi(:, kts-1) = 0.0_wp
    q_fluss  = 0.0_wp
    n_fluss  = 0.0_wp

    DO k = kts,kte
      DO i = its,ite
        IF (qp(i,k) > q_crit) THEN

          x_p = ptype%meanmass(qp(i,k),np(i,k))
          lam = EXP(ptype%b_vel* LOG(pcoeffs%coeff_lambda*x_p))

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
        ELSE
          v_n_sedi(i,k) = 0.0_wp
          v_q_sedi(i,k) = 0.0_wp
        END IF
      END DO
    END DO

    IF (PRESENT(cmax)) THEN
      cmax_temp = cmax
    ELSE
      cmax_temp = 0.0_wp
    END IF
    DO k = kts,kte

        DO i = its,ite
          v_nv(i) = 0.5 * (v_n_sedi(i,k-1)+v_n_sedi(i,k))   ! Das sollte wohl k-1 sein.
          v_qv(i) = 0.5 * (v_q_sedi(i,k-1)+v_q_sedi(i,k))
          ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
          c_nv(i) = -v_nv(i) * adz(i,k) * dt
          c_qv(i) = -v_qv(i) * adz(i,k) * dt
          cmax_temp = MAX(cmax_temp, c_qv(i))
        END DO

        kk = k
        s_nv = 0.0
        DO i = its,ite
          IF (c_nv(i) <= 1._wp) THEN
            s_nv(i) = v_nv(i) * np(i,k)
          END IF
        END DO
        IF (ANY(c_nv > 1._wp)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_nv > 1._wp) .AND. kk > 2)
            DO i = its,ite
              IF (c_nv(i) > 1._wp) THEN
                cflag(i) = .TRUE.
                s_nv(i) = s_nv(i) + np(i,kk)/adz(i,kk)
                c_nv(i) = (c_nv(i) - 1._wp) * adz(i,kk-1)/adz(i,kk)
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
        IF (ANY(c_qv > 1._wp)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_qv > 1._wp) .AND. kk > 2)
            DO i = its,ite
              IF (c_qv(i) > 1._wp) THEN
                cflag(i) = .TRUE.
                s_qv(i) = s_qv(i) + qp(i,kk)/adz(i,kk)
                c_qv(i) = (c_qv(i) - 1._wp) * adz(i,kk-1)/adz(i,kk)
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

    IF (PRESENT(cmax)) cmax = cmax_temp

    n_fluss(:,kts-1) = 0._wp ! obere Randbedingung
    q_fluss(:,kts-1) = 0._wp ! upper BC

    DO k = kts,kte
        DO i = its,ite
          np(i,k) = np(i,k) + ( n_fluss(i,k) - n_fluss(i,k-1) )*adz(i,k)*dt
          qp(i,k) = qp(i,k) + ( q_fluss(i,k) - q_fluss(i,k-1) )*adz(i,k)*dt
        ENDDO
    ENDDO
    precrate(its:ite) = - q_fluss(its:ite,kte) ! Regenrate

    DO k=kts,kte
      DO i = its,ite
        np(i,k) = MIN(np(i,k), qp(i,k)/ptype%x_min)
        np(i,k) = MAX(np(i,k), qp(i,k)/ptype%x_max)
      END DO
    END DO

    RETURN

  END SUBROUTINE sedi_icon_sphere

  SUBROUTINE check(ik_slice, mtxt,cloud,rain,ice,snow,graupel,hail)
    ! start and end indices for 2D slices
    ! istart = slice(1), iend = slice(2), kstart = slice(3), kend = slice(4)
    INTEGER, INTENT(in) :: ik_slice(4)
    CHARACTER(len=*), INTENT(in) :: mtxt
    TYPE(particle), INTENT(in) :: cloud, rain, ice, snow, graupel, hail
    INTEGER :: k, kstart, kend
    REAL(wp), PARAMETER  :: meps = -1e-12

    kstart = ik_slice(3)
    kend   = ik_slice(4)

    DO k = kstart,kend
       IF (MINVAL(cloud%q(:,k)) < meps) THEN
          WRITE (txt,'(1X,A,I4,A)') '  qc < 0 at k = ',k,' after '//TRIM(mtxt)
          CALL message(routine,TRIM(txt))
          CALL finish(TRIM(routine),txt)
       ENDIF
       IF (MINVAL(rain%q(:,k)) < meps) THEN
          WRITE (txt,'(1X,A,I4,A)') '  qr < 0 at k = ',k,' after '//TRIM(mtxt)
          CALL message(routine,TRIM(txt))
          CALL finish(TRIM(routine),txt)
       ENDIF
       IF (MINVAL(ice%q(:,k)) < meps) THEN
          WRITE (txt,'(1X,A,I4,A)') '  qi < 0 at k = ',k,' after '//TRIM(mtxt)
          CALL message(routine,TRIM(txt))
          CALL finish(TRIM(routine),txt)
       ENDIF
       IF (MINVAL(snow%q(:,k)) < meps) THEN
          WRITE (txt,'(1X,A,I4,A)') '  qs < 0 at k = ',k,' after '//TRIM(mtxt)
          CALL message(routine,TRIM(txt))
          CALL finish(TRIM(routine),txt)
       ENDIF
       IF (MINVAL(graupel%q(:,k)) < meps) THEN
          WRITE (txt,'(1X,A,I4,A)') '  qg < 0 at k = ',k,' after '//TRIM(mtxt)
          CALL message(routine,TRIM(txt))
          CALL finish(TRIM(routine),txt)
       ENDIF
       IF (MINVAL(hail%q(:,k)) < meps) THEN
          WRITE (txt,'(1X,A,I4,A)') '  qh < 0 at k = ',k,' after '//TRIM(mtxt)
          CALL message(routine,TRIM(txt))
          CALL finish(TRIM(routine),txt)
       ENDIF
    END DO

  END SUBROUTINE check

END MODULE mo_2mom_mcrph_main
