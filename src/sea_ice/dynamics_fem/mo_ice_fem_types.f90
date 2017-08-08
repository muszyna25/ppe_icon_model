!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_ice_fem_types
!
! Arrays used to store ice variables and organize coupling
!
  USE mo_kind,    ONLY: wp

  IMPLICIT NONE

  PUBLIC
  SAVE

  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: u_ice, v_ice, m_ice, a_ice  
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: rhs_u, rhs_v, m_snow
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: rhs_m, rhs_a, rhs_mis
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: u_w, v_w
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: elevation
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: mass_matrix
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: lmass_matrix
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: sigma11, sigma12, sigma22  
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: stress_atmice_x         
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: stress_atmice_y
  ! auxiliary arrays required for fct
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: m_icel, a_icel, m_snowl
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: dm_ice, da_ice, dm_snow
  REAL(wp), ALLOCATABLE, DIMENSION(:,:)       :: icefluxes
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: icepplus, icepminus

! There are, in principle, 3 matrices which share the same structure stored in icestiff.
! We use icestiff%values for both implicit advection schemes and VP algorithms.
! Since they are re-assembled at each time step, there is no problem with sharing it.
! Additionally, there is mass matrix with the same structure, which is used in FCT advection scheme.
! So we use mass_matrix of the size of icestiff%values only to keep  the values.

  TYPE sparse_matrix
     INTEGER :: nza
     INTEGER :: dim
     REAL(wp),          POINTER, DIMENSION(:)   :: values
     INTEGER(KIND=4),   POINTER, DIMENSION(:)   :: colind
     INTEGER(KIND=4),   POINTER, DIMENSION(:)   :: rowptr
  END TYPE sparse_matrix

  TYPE(sparse_matrix) :: icestiff

END MODULE mo_ice_fem_types
