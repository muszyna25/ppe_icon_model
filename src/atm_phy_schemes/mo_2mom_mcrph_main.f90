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
! Version of November, 2015, by D. Rieger:
! - Created an own module for process subroutines (mo_2mom_mcrph_processes)
! - Extended the argument lists of several process subroutines as the declaration
!   of the according variables happens still in mo_2mom_mcrph_main
! - Usage of a separate routine for type declaration as developed by AS
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

  USE mo_kind,               ONLY: sp, wp
  USE mo_exception,          ONLY: finish, message, txt => message_text
  USE mo_math_constants,     ONLY: pi
  USE mo_physical_constants, ONLY:    &
       & R_l   => rd,                 & ! gas constant of dry air (luft)
       & R_d   => rv,                 & ! gas constant of water vapor (dampf)
       & nu_l  => con_m                 ! kinematic viscosity of air

  USE mo_satad, ONLY:                 &
       & e_ws  => sat_pres_water,     & ! saturation pressure over liquid water
       & e_es  => sat_pres_ice          ! saturation pressure over ice

  USE mo_2mom_mcrph_util, ONLY: &
       & gfct,                        & ! Gamma function (becomes intrinsic in Fortran2008)
       & gamlookuptable,              & ! For look-up table of incomplete Gamma function
       & nlookup, nlookuphr_dummy,    & !   array size of table
       & incgfct_lower_lookupcreate,  & !   create table
       & incgfct_lower_lookup,        & !   interpolation in table, lower incomplete Gamma function
       & incgfct_upper_lookup,        & !   interpolation in talbe, upper incomplete Gamma function
       & dmin_wg_gr_ltab_equi,        & ! For look-up table of wet growth diameter
       & ltabdminwgg

  USE mo_2mom_mcrph_processes, ONLY:                                         &
       &  coll_delta_11, coll_delta_12, coll_delta_22,                       &
       &  coll_theta_11, coll_theta_12, coll_theta_22,                       &
       &  vent_coeff_a, vent_coeff_b, moment_gamma,                          &
       &  sedi_vel_rain, init_sedi_vel, autoconversionSB,                    &
       &  accretionSB, rain_selfcollectionSB, autoconversionKB, accretionKB, &
       &  autoconversionKK, accretionKK, rain_evaporation, evaporation,      &
       &  cloud_freeze, ice_nucleation_homhet,                               &
       &  vapor_dep_relaxation, rain_freeze_gamlook,                         &
       &  setup_ice_selfcollection, ice_selfcollection,                      &
       &  setup_snow_selfcollection, snow_selfcollection,                    &
       &  snow_melting, graupel_snow_collection, hail_snow_collection,       &
       &  graupel_ice_collection, hail_ice_collection, snow_ice_collection,  &
       &  setup_graupel_selfcollection, graupel_selfcollection, ice_melting, &
       &  prtcl_cloud_riming, prtcl_rain_riming, graupel_melting,            &
       &  hail_melting, graupel_hail_conv_wet_gamlook, ice_riming,           &
       &  snow_riming, ccn_activation_sk, ccn_activation_hdcp2
  ! And some parameters declared in the process module that are required here. 
  USE mo_2mom_mcrph_processes, ONLY:                                         &
       &  N_sc, n_f, ecoll_gc, ecoll_hc,                                     &
       &  q_crit_gc, D_crit_gc, q_crit_hc, D_crit_hc,                        &
       &  D_shed_g, D_shed_h, D_crit_c, D_coll_c, rainSBBcoeffs
  ! And switches...
  USE mo_2mom_mcrph_processes, ONLY:                                         &
       &  ice_typ, nuc_i_typ, nuc_c_typ, auto_typ, isdebug, ischeck

  USE mo_2mom_mcrph_types,     ONLY:             &
       &  ATMOSPHERE, PARTICLE, particle_sphere, &
       &  particle_rain_coeffs,                  &
       &  aerosol_ccn, aerosol_in,               &
       &  sym_riming_params, asym_riming_params, &
       &  evaporation_deposition_params,         &
       &  melt_params

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_main'

  ! for constant droplet number runs
  real(wp) :: qnc_const = 200.0e6_wp

  ! These are the fundamental particle variables for the two-moment scheme
  ! These pointers will be specified from the list of pre-defined particles
  TYPE(particle_sphere)          :: ice_coeffs, snow_coeffs, graupel_coeffs, hail_coeffs
  TYPE(particle_rain_coeffs)     :: rain_coeffs
  TYPE(aerosol_ccn)              :: ccn_coeffs
  TYPE(aerosol_in)               :: in_coeffs

  ! Look up tables for graupel_hail_conv_wet_gamlook and rain_freeze_gamlook

  TYPE(gamlookuptable) :: graupel_ltable1, graupel_ltable2
  TYPE(gamlookuptable) :: rain_ltable1, rain_ltable2, rain_ltable3
  REAL(wp)             :: rain_nm1, rain_nm2, rain_nm3, rain_g1, rain_g2
  REAL(wp)             :: graupel_nm1, graupel_nm2, graupel_g1, graupel_g2

  ! choice of mu-D relation for rain, default is mu_Dm_rain_typ = 1
  INTEGER, PARAMETER     :: mu_Dm_rain_typ = 1     ! see init_twomoment() for possible choices

  ! Parameter for evaporation of rain, determines change of n_rain during evaporation
  REAL(wp) :: rain_gfak   ! this is set in init_twomoment depending on mu_Dm_rain_typ

  ! switches for shedding
  ! UB: Das angefrorene Wasser wird bei T > T_shed,
  !     ggf. nach enhanced-melting, wieder abgeworfen und zu Regen.
  !     Default is .false. for both parameters, because this simple bulk approach
  !     would limit the growth of graupel and hail too much.
  LOGICAL, PARAMETER     :: graupel_shedding = .FALSE.
  LOGICAL, PARAMETER     :: hail_shedding    = .FALSE.

  PUBLIC :: init_2mom_scheme, init_2mom_scheme_once, clouds_twomoment

  ! these should be fixed to private for reasons of encapsulation
  PUBLIC :: rain_coeffs, ice_coeffs, snow_coeffs, graupel_coeffs, hail_coeffs, &
       ccn_coeffs, in_coeffs

  PUBLIC :: qnc_const

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
  TYPE(sym_riming_params), SAVE :: grr_params
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

  !> run-time- and location invariant for cloud freeze
  REAL(wp), SAVE :: cloud_freeze_coeff_z

  !> run-time- and location invariants for ice selfcollection
  REAL(wp), SAVE :: ice_sc_delta_n, ice_sc_delta_q, &
       ice_sc_theta_n, ice_sc_theta_q

  !> run-time- and location invariants for rain_freeze_gamlook
  REAL(wp), SAVE :: rain_freeze_coeff_z

  !> run-time- and location invariants for autoconversionSB
  REAL(wp), SAVE :: autoconversion_sb_k_au, autoconversion_sb_k_sc

  
  ! DR: The following block is necessarily public as these parameters/coefficients
  !     are required by the ART code
  PUBLIC :: scr_params, srr_params, irr_params, icr_params
  PUBLIC :: hrr_params, grr_params, hcr_params, gcr_params
  PUBLIC :: sic_params, hic_params, gic_params, hsc_params, gsc_params
  PUBLIC :: vid_params, vgd_params, vhd_params, vsd_params
  PUBLIC :: ge_params, se_params, gm_params, hm_params, sm_params
  PUBLIC :: graupel_sc_coll_n, snow_sc_delta_n, snow_sc_theta_n
  PUBLIC :: cloud_freeze_coeff_z, ice_sc_delta_n, ice_sc_delta_q
  PUBLIC :: ice_sc_theta_n, ice_sc_theta_q
  PUBLIC :: rain_freeze_coeff_z, autoconversion_sb_k_au, autoconversion_sb_k_sc
  PUBLIC :: graupel_ltable1, graupel_ltable2
  PUBLIC :: rain_ltable1, rain_ltable2, rain_ltable3
  PUBLIC :: rain_nm1, rain_nm2, rain_nm3, rain_g1, rain_g2
  PUBLIC :: graupel_nm1, graupel_nm2, graupel_g1, graupel_g2
  ! END DR

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
    REAL(wp) :: nu, mu, x_s_i
    REAL(wp), PARAMETER :: k_c  = 9.44e+9_wp   !..Long-Kernel


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
    grr_params%delta_q_bb = coll_delta_22(graupel,rain,1)

    grr_params%theta_n_aa = coll_theta_11(graupel,rain,0)
    grr_params%theta_n_ab = coll_theta_12(graupel,rain,0)
    grr_params%theta_n_bb = coll_theta_22(graupel,rain,0)
    grr_params%theta_q_aa = coll_theta_11(graupel,rain,0)
    grr_params%theta_q_ab = coll_theta_12(graupel,rain,1)
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
      WRITE(txt,'(A,D10.3)') "     delta_q_rr = ",grr_params%delta_q_bb ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_gg = ",grr_params%theta_q_aa ; CALL message(routine,TRIM(txt))
      WRITE(txt,'(A,D10.3)') "     theta_q_gr = ",grr_params%theta_q_ab ; CALL message(routine,TRIM(txt))
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

    CALL setup_graupel_selfcollection(graupel,graupel_sc_coll_n)
    CALL setup_snow_selfcollection(snow, snow_sc_delta_n, snow_sc_theta_n)
    CALL setup_ice_selfcollection(ice, ice_sc_delta_n, ice_sc_delta_q, ice_sc_theta_n, ice_sc_theta_q)

    ! setup coefficient for cloud_freeze
    cloud_freeze_coeff_z = moment_gamma(cloud,2)

    rain_freeze_coeff_z = moment_gamma(rain,2)  ! coeff for 2nd moment
    IF (isdebug) THEN
      CALL message(routine, "rain_freeze_gamlook:")
      WRITE(txt,'(A,D10.3)') "    coeff_z= ",rain_freeze_coeff_z
      CALL message(routine,TRIM(txt))
    ENDIF

    nu = cloud%nu
    mu = cloud%mu
    x_s_i = 1.0_wp / cloud%x_max

    IF (mu == 1.0) THEN
      !.. see SB2001
      autoconversion_sb_k_au  = k_c * x_s_i * (1.0_wp / 20.0_wp) * (nu+2.0)*(nu+4.0)/(nu+1.0)**2
      autoconversion_sb_k_sc  = k_c * (nu+2.0)/(nu+1.0)
    ELSE
      !.. see Eq. (3.44) of Seifert (2002)
      autoconversion_sb_k_au = k_c * x_s_i * (1.0_wp / 20.0_wp)                            &
           & * ( 2.0_wp * gfct((nu+4.0)/mu)**1                          &
           &            * gfct((nu+2.0)/mu)**1 * gfct((nu+1.0)/mu)**2   &
           &   - 1.0_wp * gfct((nu+3.0)/mu)**2 * gfct((nu+1.0)/mu)**2 ) &
           &   / gfct((nu+2.0)/mu)**4
      autoconversion_sb_k_sc = k_c * moment_gamma(cloud,2)
    ENDIF


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
          CALL ccn_activation_sk(ik_slice, ccn_coeffs, atmo, cloud, n_cn)
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
      CALL ice_nucleation_homhet(ik_slice, use_prog_in, atmo, cloud, ice, snow, n_inact, n_inpot)

       ! homogeneous freezing of cloud droplets
       CALL cloud_freeze(ik_slice, dt, cloud_freeze_coeff_z, qnc_const, atmo, cloud, ice)
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
       CALL vapor_dep_relaxation(ik_slice, vid_params, vgd_params, vhd_params, vsd_params, dt, &
         &                       atmo,ice,snow,graupel,hail,dep_rate_ice,dep_rate_snow)
       IF (ischeck) CALL check(ik_slice, 'vapor_dep_relaxation',cloud,rain,ice,snow,graupel,hail)

       if (.true.) then

       ! ice-ice collisions
       CALL ice_selfcollection(ik_slice, dt, ice_sc_delta_n, ice_sc_delta_q, &
         &                     ice_sc_theta_n, ice_sc_theta_q, atmo, ice, snow)
       CALL snow_selfcollection(ik_slice, dt, snow_sc_delta_n, snow_sc_theta_n, atmo, snow)
       CALL snow_ice_collection(ik_slice, dt, sic_params, atmo, ice, snow)
       IF (ischeck) CALL check(ik_slice, 'ice and snow collection',cloud,rain,ice,snow,graupel,hail)

       CALL graupel_selfcollection(ik_slice, dt, graupel_sc_coll_n, atmo, graupel)
       CALL graupel_ice_collection(ik_slice, dt, gic_params, atmo, ice, graupel)
       CALL graupel_snow_collection(ik_slice, dt, gsc_params, atmo, snow, graupel)
       IF (ischeck) CALL check(ik_slice, 'graupel collection',cloud,rain,ice,snow,graupel,hail)

       IF (ice_typ > 1) THEN

          ! conversion of graupel to hail in wet growth regime
          CALL graupel_hail_conv_wet_gamlook(ik_slice, graupel_ltable1, graupel_ltable2,       &
                                             graupel_nm1, graupel_nm2, graupel_g1, graupel_g2, &
                                             atmo, graupel, cloud, rain, ice, snow, hail)
          IF (ischeck) CALL check(ik_slice, 'graupel_hail_conv_wet_gamlook',cloud,rain,ice,snow,graupel,hail)

          ! hail collisions
          CALL hail_ice_collection(ik_slice, dt, hic_params, atmo, ice, hail)    ! Important?
          CALL hail_snow_collection(ik_slice, dt, hsc_params, atmo, snow, hail)
          IF (ischeck) CALL check(ik_slice, 'hail collection',cloud,rain,ice,snow,graupel,hail)
       END IF

       ! riming of ice with cloud droplets and rain drops, and conversion to graupel
       CALL ice_riming(ik_slice, dt, icr_params, irr_params, atmo, ice, cloud, rain, graupel, &
            dep_rate_ice)
       IF (ischeck) CALL check(ik_slice, 'ice_riming',cloud,rain,ice,snow,graupel,hail)

       ! riming of snow with cloud droplets and rain drops, and conversion to graupel
       CALL snow_riming(ik_slice, dt, scr_params, srr_params, atmo, &
            snow, cloud, rain, ice, graupel, dep_rate_snow)
       IF (ischeck) CALL check(ik_slice, 'snow_riming',cloud,rain,ice,snow,graupel,hail)

       ! more riming
       IF (ice_typ > 1) THEN
          CALL prtcl_cloud_riming(ik_slice, dt, atmo, hail, &
               ecoll_hc/(D_coll_c - D_crit_c), ecoll_hc, d_shed_h, &
               q_crit_hc, d_crit_hc, hail_shedding, hcr_params, &
               cloud, rain, ice)
          CALL prtcl_rain_riming(ik_slice, dt, atmo, hail, hrr_params, &
               d_shed_h, hail_shedding, rain, ice)
         IF (ischeck) CALL check(ik_slice, 'hail riming',cloud,rain,ice,snow,graupel,hail)
       END IF
       CALL prtcl_cloud_riming(ik_slice, dt, atmo, graupel, &
            ecoll_gc/(D_coll_c - D_crit_c), ecoll_gc, d_shed_g, &
            q_crit_gc, d_crit_gc, graupel_shedding, gcr_params, &
            cloud, rain, ice)
       CALL prtcl_rain_riming(ik_slice, dt, atmo, graupel, grr_params, &
            d_shed_g, graupel_shedding, rain, ice)
       IF (ischeck) CALL check(ik_slice, 'graupel riming',cloud,rain,ice,snow,graupel,hail)

       ! freezing of rain and conversion to ice/graupel/hail
       CALL rain_freeze_gamlook(ik_slice, dt, rain_ltable1, rain_ltable2, rain_ltable3, &
                                rain_nm1, rain_nm2, rain_nm3, rain_g1, rain_g2,         &
                                rain_freeze_coeff_z, atmo,rain,ice,snow,graupel,hail)
       IF (ischeck) CALL check(ik_slice, 'rain_freeze_gamlook',cloud,rain,ice,snow,graupel,hail)

       ! melting
       CALL ice_melting(ik_slice, atmo, ice, cloud, rain)
       CALL snow_melting(ik_slice, dt, sm_params, atmo,snow,rain)
       CALL graupel_melting(ik_slice, dt, gm_params, atmo, graupel, rain)
       IF (ice_typ > 1) CALL hail_melting(ik_slice, dt, hm_params, atmo, hail, rain)
       IF (ischeck) CALL check(ik_slice, 'melting',cloud,rain,ice,snow,graupel,hail)

       ! evaporation from melting ice particles
       CALL evaporation(ik_slice, dt, atmo, snow, se_params)
       CALL evaporation(ik_slice, dt, atmo, graupel, ge_params)
       IF (ice_typ > 1) CALL evaporation(ik_slice, dt, atmo, hail, he_params)
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
       CALL autoconversionSB(ik_slice, dt, autoconversion_sb_k_au, autoconversion_sb_k_sc, &
         &                   cloud, rain)  ! Seifert and Beheng (2001)
       CALL accretionSB(ik_slice, dt, cloud, rain)
       CALL rain_selfcollectionSB(ik_slice, dt, rain)
    ENDIF
    IF (ischeck) CALL check(ik_slice, 'warm rain',cloud,rain,ice,snow,graupel,hail)

    ! evaporation of rain following Seifert (2008)
    CALL rain_evaporation(ik_slice, dt, rain_coeffs, rain_gfak, atmo, cloud, rain)

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
