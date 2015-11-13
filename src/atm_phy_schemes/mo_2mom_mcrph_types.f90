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

MODULE mo_2mom_mcrph_types

  USE mo_kind,               ONLY: sp, wp
  USE mo_exception,          ONLY: finish, message, txt => message_text

  IMPLICIT NONE

  PRIVATE

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
    REAL(wp)      :: a_vel      !..pre-factor in power law fall speed (all particles have a power law fall speed,
    REAL(wp)      :: b_vel      !..exponent in power law fall speed    some have also an Atlas-type relation)
    REAL(wp)      :: a_ven      !..first parameter in ventilation coefficient
    REAL(wp)      :: b_ven      !..2nd parameter in ventilation coefficient
    REAL(wp)      :: cap        !..coefficient for capacity of particle
    REAL(wp)      :: vsedi_max  !..max bulk sedimentation velocity
    REAL(wp)      :: vsedi_min  !..min bulk sedimentation velocity
    REAL(wp), pointer, dimension(:,:) :: n     !..number density
    REAL(wp), pointer, dimension(:,:) :: q     !..mass density
    REAL(wp), pointer, dimension(:,:) :: rho_v !..density correction of terminal fall velocity
  END TYPE PARTICLE

  TYPE, EXTENDS(particle) :: particle_frozen
    REAL(wp)      :: ecoll_c    !..maximum collision efficiency with cloud droplets
    REAL(wp)      :: D_crit_c   !..D-threshold for cloud riming
    REAL(wp)      :: q_crit_c   !..q-threshold for cloud riming
    REAL(wp)      :: s_vel      !..dispersion of fall velocity for collection kernel (see SB2006, Eqs 60-63)
  END TYPE particle_frozen

  TYPE, EXTENDS(particle_frozen) :: particle_lwf
    REAL(wp)      :: lwf_cnorm1  !..1st parameter for normalized diameter
    REAL(wp)      :: lwf_cnorm2  !..2nd parameter for normalized diameter
    REAL(wp)      :: lwf_cnorm3  !..3rd parameter for normalized diameter
    REAL(wp)      :: lwf_cmelt1  !..1st parameter for melting integral
    REAL(wp)      :: lwf_cmelt2  !..2nd parameter for melting integral
    REAL(wp), pointer, dimension(:,:) :: l  !..mass density of liquid water on ice (per unit volume of air)
  END TYPE particle_lwf

  ! .. Because of OpenMP we have to separate the data pointers from the run-time-invariant coefficients.
  !    Therefore we carry 2 data structures for each particle species, e.g. graupel and graupel_coeff.
  !    The following derived types are for the run-time coefficients
  
  TYPE particle_params
    REAL(wp)      :: a_f  ! ventilation coefficient, vent_coeff_a(particle,1)
    REAL(wp)      :: b_f  ! ventilation coefficient, vent_coeff_b(particle,1) * N_sc**n_f / SQRT(nu_l)
    REAL(wp)      :: c_i  ! 1.0/particle%cap
  END type particle_params
  
  ! .. for spherical particles we need to store the coefficients for the
  !    power law bulk sedimentation velocity
  TYPE, EXTENDS(particle_params) :: particle_sphere
    REAL(wp)      :: coeff_alfa_n
    REAL(wp)      :: coeff_alfa_q
    REAL(wp)      :: coeff_lambda
  END TYPE particle_sphere

  ! .. non-spherical particles have an Atlas-type terminal fall velocity relation
  TYPE, EXTENDS(particle_params) :: particle_nonsphere
    REAL(wp)      :: alfa   !..1st parameter in Atlas-type fall speed
    REAL(wp)      :: beta   !..2nd parameter in Atlas-type fall speed
    REAL(wp)      :: gama   !..3rd parameter in Atlas-type fall speed
  END TYPE particle_nonsphere

  ! .. raindrops have an Atlas-type terminal fall velocity relation
  !    and a mu-D-relation which is used in sedimentation and evaporation
  !    (see Seifert 2008, J. Atmos. Sci.)
  TYPE, EXTENDS(particle_nonsphere) :: particle_rain_coeffs
    REAL(wp)      :: cmu0   !..Parameters for mu-D-relation of rain
    REAL(wp)      :: cmu1   !     max of left branch
    REAL(wp)      :: cmu2   !     max of right branch
    REAL(wp)      :: cmu3   !     min value of relation
    REAL(wp)      :: cmu4   !     location of min value = breakup equilibrium diameter
    INTEGER       :: cmu5   !     exponent
  END TYPE particle_rain_coeffs

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

    !..these are coefficients for collection processes of the type a+b->a
  TYPE collection_params
     REAL(wp) :: delta_n_aa, delta_n_ab, delta_n_bb, &
          &      delta_q_aa, delta_q_ab, delta_q_bb, &
          &      theta_n_aa, theta_n_ab, theta_n_bb, &
          &      theta_q_aa, theta_q_ab, theta_q_bb
  END TYPE collection_params

  !..these are coefficients for collection processes of the type a+b->c
  TYPE rain_riming_params
    REAL(wp) :: delta_n_aa,delta_n_ab,delta_n_bb, &
         &      delta_q_aa,delta_q_ab,delta_q_ba,delta_q_bb, &
         &      theta_n_aa,theta_n_ab,theta_n_bb, &
         &      theta_q_aa,theta_q_ab,theta_q_ba,theta_q_bb
  END TYPE rain_riming_params

  TYPE dep_imm_params
    REAL(wp) :: alf_dep, bet_dep, nin_dep, &
                alf_imm, bet_imm, nin_imm
  END TYPE dep_imm_params

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

  TYPE evaporation_deposition_params
    REAL(wp) :: c          !< coeff for capacity
    REAL(wp) :: a_f,b_f    !< coeffs for ventilation
  END TYPE evaporation_deposition_params

  TYPE melt_params
    REAL(wp) :: a_vent, b_vent
  END TYPE melt_params

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_types'

  PUBLIC :: ATMOSPHERE
  PUBLIC :: PARTICLE
  PUBLIC :: particle_frozen, particle_lwf, particle_params
  PUBLIC :: particle_sphere, particle_nonsphere
  PUBLIC :: particle_rain_coeffs
  PUBLIC :: aerosol_ccn, aerosol_in
  PUBLIC :: sym_riming_params, asym_riming_params
  PUBLIC :: collection_params, rain_riming_params
  PUBLIC :: dep_imm_params, evaporation_deposition_params, melt_params
  
END MODULE mo_2mom_mcrph_types
