!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
module mo_ice_elements
  !
  ! Arrays defined here are used to keep mesh information
  ! Separation between MESH and ELEMENTS is for
  ! historical reasons (and is inherited from FESOM)
  !
  USE mo_kind,    ONLY: wp
  implicit none
  PUBLIC
  save
  integer                                      :: elem2D
  integer(KIND=4), allocatable, dimension(:,:) :: elem2D_nodes
  integer(KIND=4), allocatable, dimension(:,:) :: nod2D_elems
  integer(KIND=4), allocatable, dimension(:,:) :: elem2D_nghbrs
  REAL(wp)                                 :: Vol2D
  REAL(wp), allocatable, dimension(:,:)    :: derivative_stdbf
  REAL(wp), allocatable, dimension(:,:)    :: bafux, bafuy
  REAL(wp), allocatable, dimension(:)      :: voltriangle
end module mo_ice_elements
