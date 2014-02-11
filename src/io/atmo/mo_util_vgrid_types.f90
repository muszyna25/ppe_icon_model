!>
!! Module containing types for the reading / on-the-fly generation and the writing of the vertical grid.
!!
!! @author L. Kornblueh, F. Prill
!!
!! @par Revision History
!! Initial implementation by L. Kornblueh (MPI-M), F. Prill (DWD) : 2014-01-28
!!
MODULE mo_util_vgrid_types

  USE mo_kind,                              ONLY: wp
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
#ifdef _CRAYFTN
      , CONTIGUOUS         &
#endif
      ::                    &
      z_ifc(:,:,:)

    ! UUID of vertical grid
    CHARACTER(len=1) :: uuid(16)
  END TYPE t_vgrid_buffer

  ! module variable: temporary buffer for coordinate arrays
  TYPE (t_vgrid_buffer), ALLOCATABLE :: vgrid_buffer(:)


END MODULE mo_util_vgrid_types
