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
