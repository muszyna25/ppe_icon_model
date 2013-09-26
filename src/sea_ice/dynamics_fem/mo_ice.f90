
MODULE mo_ice
!
! Arrays used to store ice variables and organize coupling
!

USE mo_ice_data_types
USE mo_kind,    ONLY: wp
implicit none
PUBLIC
save
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: u_ice, v_ice, m_ice, a_ice  
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: rhs_u, rhs_v, m_snow
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: rhs_m, rhs_a, rhs_ms
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: u_w, v_w
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: elevation
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: mass_matrix
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: lmass_matrix
  TYPE(sparse_matrix)                         :: icestiff
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: sigma11, sigma12, sigma22  
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: fresh_wa_flux
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: net_heat_flux
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: S_oc_array, T_oc_array
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: stress_iceoce_x         
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: stress_iceoce_y
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: stress_atmice_x         
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: stress_atmice_y
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: stress_atmoce_x         
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: stress_atmoce_y
  ! auxiliary arrays required for fct
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: m_icel, a_icel, m_snowl
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: dm_ice, da_ice, dm_snow
  REAL(wp), ALLOCATABLE, DIMENSION(:,:)       :: icefluxes
  REAL(wp), ALLOCATABLE, DIMENSION(:)         :: icepplus, icepminus
END MODULE mo_ice

!=====================================================================
