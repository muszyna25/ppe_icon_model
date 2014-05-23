!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!----------------------------------------------------------------------------
!
module mo_ice_mesh
  !
  ! Arrays defined here are used to keep mesh information
  !
  use mo_ice_data_types
  USE mo_kind,    ONLY: wp
  implicit none
  PUBLIC
  save
  !
  integer                                      :: nod2D
  REAL(wp), allocatable, dimension(:,:)    :: coord_nod2D
  integer, allocatable, dimension(:)           :: index_nod2D
  REAL(wp), allocatable, dimension(:)      :: cos_elem2D
  REAL(wp), allocatable, dimension(:)      :: sin_elem2D
  integer, allocatable, dimension(:)           :: col_pos
  !
  type(addresstype), allocatable, dimension(:) :: nod_in_elem2D
  type(addresstype), allocatable, dimension(:) :: nghbr_nod2D
  !
end module mo_ice_mesh
!
!----------------------------------------------------------------------------
