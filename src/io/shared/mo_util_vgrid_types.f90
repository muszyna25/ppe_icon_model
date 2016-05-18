!>
!! Module containing types for the reading / on-the-fly generation and the writing of the vertical grid.
!!
!! @author L. Kornblueh, F. Prill
!!
!! @par Revision History
!! Initial implementation by L. Kornblueh (MPI-M), F. Prill (DWD) : 2014-01-28
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_util_vgrid_types

  USE mo_kind, ONLY: wp
  USE mo_util_uuid, ONLY: t_uuid
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_vgrid_buffer
  PUBLIC :: vgrid_buffer

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_vgrid_types'


  TYPE t_vgrid_buffer
    ! Geom. height at vertical interface of cells and vertices (nproma,nlevp1,nblks_c/nblks_v): 
    !
    ! Note: This array is deallocated after being used in
    !       mo_vertical_grid::set_nh_metrics
    !
    REAL(wp), POINTER      &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS         &
#endif
      ::                    &
      z_ifc(:,:,:)

    ! UUID of vertical grid
    TYPE(t_uuid) :: uuid
  END TYPE t_vgrid_buffer

  ! module variable: temporary buffer for coordinate arrays
  TYPE (t_vgrid_buffer), ALLOCATABLE :: vgrid_buffer(:)


END MODULE mo_util_vgrid_types
