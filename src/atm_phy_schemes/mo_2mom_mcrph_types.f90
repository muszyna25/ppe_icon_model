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

  PUBLIC

  ! Derived type for atmospheric variables
  TYPE ATMOSPHERE
     REAL(wp), pointer, dimension(:,:) :: w, p, t, rho, qv, zh
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
  
  TYPE particle_coeffs
    REAL(wp)      :: a_f  ! ventilation coefficient, vent_coeff_a(particle,1)
    REAL(wp)      :: b_f  ! ventilation coefficient, vent_coeff_b(particle,1) * N_sc**n_f / SQRT(nu_l)
    REAL(wp)      :: c_i  ! 1.0/particle%cap
    REAL(wp)      :: c_z  ! coefficient for 2nd mass moment
  END type particle_coeffs
  
  ! .. for spherical particles we need to store the coefficients for the
  !    power law bulk sedimentation velocity
  TYPE, EXTENDS(particle_coeffs) :: particle_sphere
    REAL(wp)      :: coeff_alfa_n
    REAL(wp)      :: coeff_alfa_q
    REAL(wp)      :: coeff_lambda
  END TYPE particle_sphere

  ! .. non-spherical particles have an Atlas-type terminal fall velocity relation
  TYPE, EXTENDS(particle_coeffs) :: particle_nonsphere
    REAL(wp)      :: alfa   !..1st parameter in Atlas-type fall speed
    REAL(wp)      :: beta   !..2nd parameter in Atlas-type fall speed
    REAL(wp)      :: gama   !..3rd parameter in Atlas-type fall speed
  END TYPE particle_nonsphere

  ! .. raindrops have an Atlas-type terminal fall velocity relation
  !    and a mu-D-relation which is used in sedimentation and evaporation
  !    (see Seifert 2008, J. Atmos. Sci.)
  TYPE, EXTENDS(particle_nonsphere) :: particle_rain_coeffs
    REAL(wp)      :: cmu0   !..Parameters for mu-D-relation of rain: max of left branch
    REAL(wp)      :: cmu1   !     max of right branch
    REAL(wp)      :: cmu2   !     scaling factor
    REAL(wp)      :: cmu3   !     location of min value = breakup equilibrium diameter
    REAL(wp)      :: cmu4   !     min value of relation
    INTEGER       :: cmu5   !     exponent
  END TYPE particle_rain_coeffs

  TYPE, EXTENDS(particle_coeffs) :: particle_cloud_coeffs
    REAL(wp)      :: k_au   !..Parameters for autoconversion
    REAL(wp)      :: k_sc   !    and selfcollection
  END TYPE particle_cloud_coeffs

  TYPE, EXTENDS(particle_sphere) :: particle_graupel_coeffs
    REAL(wp)      :: sc_coll_n  !..Parameters for self-collection
  END TYPE particle_graupel_coeffs

  TYPE, EXTENDS(particle_sphere) :: particle_snow_coeffs
    REAL(wp)      :: sc_delta_n !..Parameters for self-collection
    REAL(wp)      :: sc_theta_n !   of snow
  END TYPE particle_snow_coeffs

  TYPE, EXTENDS(particle_sphere) :: particle_ice_coeffs
    REAL(wp)      :: sc_delta_n !..Parameters for self-collection
    REAL(wp)      :: sc_delta_q !   of cloud ice
    REAL(wp)      :: sc_theta_n
    REAL(wp)      :: sc_theta_q
  END TYPE particle_ice_coeffs

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
  TYPE collection_coeffs
     REAL(wp) :: delta_n_aa, delta_n_ab, delta_n_bb, &
          &      delta_q_aa, delta_q_ab, delta_q_bb, &
          &      theta_n_aa, theta_n_ab, theta_n_bb, &
          &      theta_q_aa, theta_q_ab, theta_q_bb
  END TYPE collection_coeffs

  !..these are coefficients for collection processes of the type a+b->c
  TYPE rain_riming_coeffs
    REAL(wp) :: delta_n_aa,delta_n_ab,delta_n_bb, &
         &      delta_q_aa,delta_q_ab,delta_q_ba,delta_q_bb, &
         &      theta_n_aa,theta_n_ab,theta_n_bb, &
         &      theta_q_aa,theta_q_ab,theta_q_ba,theta_q_bb
  END TYPE rain_riming_coeffs

  TYPE dep_imm_coeffs
    REAL(wp) :: alf_dep, bet_dep, nin_dep, &
                alf_imm, bet_imm, nin_imm
  END TYPE dep_imm_coeffs

  ! Type declaration for a general 4D equidistant lookup table:
  TYPE lookupt_4D
    INTEGER :: n1  ! number of grid points in x1-direction
    INTEGER :: n2  ! number of grid points in x2-direction
    INTEGER :: n3  ! number of grid points in x3-direction
    INTEGER :: n4  ! number of grid points in x4-direction
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x1  ! grid vector in x1-direction
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x2  ! grid vector in x1-direction
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x3  ! grid vector in x1-direction
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x4  ! grid vector in x1-direction
    DOUBLE PRECISION                     :: dx1          ! dx1   (grid distance w.r.t. x1)
    DOUBLE PRECISION                     :: dx2          ! dx2   (grid distance w.r.t. x2)
    DOUBLE PRECISION                     :: dx3          ! dx3   (grid distance w.r.t. x3)
    DOUBLE PRECISION                     :: dx4          ! dx4   (grid distance w.r.t. x4)
    DOUBLE PRECISION                     :: odx1         ! one over dx 1
    DOUBLE PRECISION                     :: odx2         ! one over dx 2
    DOUBLE PRECISION                     :: odx3         ! one over dx 3
    DOUBLE PRECISION                     :: odx4         ! one over dx 4
    DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER :: ltable
  END TYPE lookupt_4D

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_types'

  PRIVATE :: routine

END MODULE mo_2mom_mcrph_types
