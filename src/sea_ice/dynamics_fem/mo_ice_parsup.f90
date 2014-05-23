!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!----------------------------------------------------------------------------!
module mo_ice_parsup
  !
  ! Arrays used to organize the parallel support
  !
  USE mo_kind,    ONLY: wp
  implicit none
  PUBLIC
  save

! Einar: No parallelization yet
!#ifdef PETSC
!#include "finclude/petsc.h"
!#else
!  include 'mpif.h'
!#endif
  integer      :: maxPEnum
  type com_struct
     integer    :: rPEnum
     integer, dimension(:), pointer :: rPE
     integer, dimension(:), pointer :: rptr
     integer, dimension(:), pointer :: rlist
     integer    :: sPEnum
     integer, dimension(:), pointer :: sPE
     integer, dimension(:), pointer :: sptr
     integer, dimension(:), pointer :: slist
  end type com_struct

  type(com_struct)   :: com_nod2D, com_nod3D
  ! Buffer arrays to store information to be communicated
  type com_array
     REAL(wp), dimension(:), pointer :: array
  end  type com_array
  type(com_array), allocatable             :: s_buff_ssh(:), r_buff_ssh(:)
  type(com_array), allocatable             :: s_buff_ts(:), r_buff_ts(:)

  ! general MPI part
  integer            :: MPIERR
  integer            :: npes
  integer            :: mype
  integer, allocatable, dimension(:)  :: part, part3D

  ! Mesh partition
  integer                             :: myDim_nod2D, eDim_nod2D
  integer, allocatable, dimension(:)  :: myList_nod2D
  integer                             :: myDim_nod3D, eDim_nod3D
  integer, allocatable, dimension(:)  :: myList_nod3D
  integer                             :: myDim_elem2D
  integer, allocatable, dimension(:)  :: myList_elem2D
  integer                             :: myDim_elem3Dr  !(regular)
  integer                             :: myDim_elem3D   !(total)
  integer, allocatable, dimension(:)  :: myList_elem3D
end module mo_ice_parsup
!
!===========================================================================
