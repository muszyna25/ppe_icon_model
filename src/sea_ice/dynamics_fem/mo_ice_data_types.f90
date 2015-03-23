!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!----------------------------------------------------------------------------
!
module mo_ice_data_types
USE mo_kind,    ONLY: wp
  !
  ! Defines structure to keep sparse matrices (VP solver)
  ! + neighborhood information
  !
  implicit none
  PUBLIC :: sparse_matrix, addresstype
  save
  type sparse_matrix
     integer :: nza
     integer :: dim
     REAL(wp), pointer, dimension(:)      :: values
     integer(KIND=4), pointer,   dimension(:) :: colind
     integer(KIND=4), pointer,   dimension(:) :: rowptr
  end type sparse_matrix
  type addresstype
     integer                                  :: nmb
     integer(KIND=4), dimension(:), pointer   :: addresses
  end type addresstype
end module mo_ice_data_types
!
!----------------------------------------------------------------------------
