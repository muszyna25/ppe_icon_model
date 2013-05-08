!-------------------------------------------------------------------------------------
!>
!! Set of methods for simple statistics
!!
!! @author Leonidas Linardakis, MPI-M
!!
!! @par Revision History
!!   First implementation by Leonidas Linardakis, MPI-M, 2012-01-19
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!-------------------------------------------------------------------------------------
MODULE mo_statistics_utils
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, message, warning, finish
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range
  USE mo_mpi,                ONLY: process_mpi_stdio_id, get_my_mpi_work_communicator, p_max, p_min, &
    & my_process_is_mpi_parallel
!   USE mo_io_units,           ONLY: nnml, filename_max
!   USE mo_namelist,           ONLY: position_nml, open_nml, positioned

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: global_minmax
  !-------------------------------------------------------------------------

  INTERFACE global_minmax
    MODULE PROCEDURE globalspace_2D_minmax
    MODULE PROCEDURE globalspace_3D_minmax
  END INTERFACE
  

CONTAINS

  !-----------------------------------------------------------------------
  !>
  FUNCTION globalspace_2D_minmax(values, subset) result(minmax)
    REAL(wp), INTENT(in) :: values(:,:)
    TYPE(t_subset_range), TARGET, OPTIONAL :: subset
    REAL(wp) :: minmax(2)

    REAL(wp) :: min_in_block, max_in_block, min_value, max_value
    INTEGER :: jb, startidx, endidx
    INTEGER :: communicator


    IF (PRESENT(subset)) THEN
      ! init the min, max values
      jb=subset%start_block
      CALL get_index_range(subset, jb, startidx, endidx)
      min_value = values(startidx, jb)
      max_value = values(startidx, jb)
!ICON_OMP_PARALLEL_DO PRIVATE(jb, startidx, endidx, min_in_block, max_in_block) reduction(min:min_value) reduction(max:max_value)
      DO jb = subset%start_block, subset%end_block
        CALL get_index_range(subset, jb, startidx, endidx)
        min_in_block = MINVAL(values(startidx:endidx, jb))
        max_in_block = MAXVAL(values(startidx:endidx, jb))
        min_value    = MIN(min_value, min_in_block)
        max_value    = MAX(max_value, max_in_block)
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    ELSE
      min_value      = values(1,1)
      max_value      = values(1,1)
!ICON_OMP_PARALLEL_DO PRIVATE(jb, min_in_block, max_in_block) reduction(min:min_value) reduction(max:max_value)
      DO jb = 1,  SIZE(values, 2)
        min_in_block = MINVAL(values(:, jb))
        max_in_block = MINVAL(values(:, jb))
        min_value      = MIN(min_value, min_in_block)
        max_value      = MAX(max_value, max_in_block)
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    ENDIF
    
    ! the global min, max is avaliable only to stdio process
    ! the global min, max is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN
      communicator = get_my_mpi_work_communicator()
      minmax(1) = p_min( min_value,  comm=communicator ) ! only mpi_all_reduce is avaliable
      minmax(2) = p_max( max_value,  comm=communicator, root=process_mpi_stdio_id )
    ELSE
      minmax(1) = min_value
      minmax(2) = max_value
    ENDIF

  END FUNCTION globalspace_2D_minmax
  !-----------------------------------------------------------------------
  

  !-----------------------------------------------------------------------
  !>
  FUNCTION globalspace_3D_minmax(values, subset, start_level, end_level) result(minmax)
    REAL(wp), INTENT(in) :: values(:,:,:)
    TYPE(t_subset_range), TARGET, OPTIONAL :: subset
    INTEGER, OPTIONAL :: start_level, end_level
    REAL(wp) :: minmax(2)

    REAL(wp) :: min_in_block, max_in_block, min_value, max_value
    INTEGER :: jb, level, startidx, endidx, start_vertical, end_vertical
    INTEGER :: communicator
    INTEGER :: idx
    CHARACTER(LEN=*), PARAMETER :: method_name='mo_statistics_utils:arc_length'


    IF (PRESENT(start_level)) THEN
      start_vertical = start_level
    ELSE
      start_vertical = 1
    ENDIF
    IF (PRESENT(end_level)) THEN
      end_vertical = start_level
    ELSE
      end_vertical = SIZE(values, 2)
    ENDIF

    IF (start_vertical > end_vertical) &
      CALL finish(method_name, "start_vertical > end_vertical")

    IF (PRESENT(subset)) THEN
      ! init the min, max values
      jb=subset%start_block
      CALL get_index_range(subset, jb, startidx, endidx)
      min_value = values(startidx, start_vertical, jb)
      max_value = values(startidx, start_vertical, jb)
!ICON_OMP_PARALLEL_DO PRIVATE(jb, level, startidx, endidx, min_in_block, max_in_block) reduction(min:min_value) reduction(max:max_value)
      DO jb = subset%start_block, subset%end_block
        CALL get_index_range(subset, jb, startidx, endidx)
        DO level = start_vertical, end_vertical
          DO idx = startidx, endidx
            WRITE(0,*) "checking ", jb, level, idx, values(idx, level, jb)
          ENDDO
          min_in_block = MINVAL(values(startidx:endidx, level, jb))
          max_in_block = MAXVAL(values(startidx:endidx, level, jb))
          min_value    = MIN(min_value, min_in_block)
          max_value    = MAX(max_value, max_in_block)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    ELSE
      min_value      = values(1, start_vertical, 1)
      max_value      = values(1, start_vertical, 1)
!ICON_OMP_PARALLEL_DO PRIVATE(jb, level, min_in_block, max_in_block) reduction(min:min_value) reduction(max:max_value)
      DO jb = 1,  SIZE(values, 3)
        DO level = start_vertical, end_vertical
          min_in_block = MINVAL(values(:, level, jb))
          max_in_block = MINVAL(values(:, level, jb))
          min_value    = MIN(min_value, min_in_block)
          max_value    = MAX(max_value, max_in_block)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    ENDIF

    ! the global min, max is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN
      communicator = get_my_mpi_work_communicator()
      minmax(1) = p_min( min_value,  comm=communicator ) ! only mpi_all_reduce is avaliable
      minmax(2) = p_max( max_value,  comm=communicator, root=process_mpi_stdio_id )
    ELSE
      minmax(1) = min_value
      minmax(2) = max_value
    ENDIF

  END FUNCTION globalspace_3D_minmax
  !-----------------------------------------------------------------------



END MODULE mo_statistics_utils

