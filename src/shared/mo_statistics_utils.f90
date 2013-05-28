!-------------------------------------------------------------------------------------
!>
!! Set of methods for simple statistics
!! NOTE: in order to get correct results make sure you provide the proper subset!
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
    & my_process_is_mpi_parallel, p_sum
!   USE mo_io_units,           ONLY: nnml, filename_max
!   USE mo_namelist,           ONLY: position_nml, open_nml, positioned

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC :: global_minmaxmean
  !-------------------------------------------------------------------------

  ! NOTE: in order to get correct results make sure you provide the proper subset (ie, owned)!
  INTERFACE global_minmaxmean
    MODULE PROCEDURE globalspace_2D_minmaxmean
    MODULE PROCEDURE globalspace_3D_minmaxmean
  END INTERFACE

CONTAINS

  !-----------------------------------------------------------------------
  !>
  FUNCTION globalspace_2D_minmaxmean(values, subset) result(minmaxmean)
    REAL(wp), INTENT(in) :: values(:,:)
    TYPE(t_subset_range), TARGET, OPTIONAL :: subset
    REAL(wp) :: minmaxmean(3)

    REAL(wp) :: min_in_block, max_in_block, min_value, max_value, sum_value, global_number_of_values
    INTEGER :: jb, startidx, endidx, number_of_values
    INTEGER :: communicator

    IF (PRESENT(subset)) THEN
      ! init the min, max values
      jb = subset%start_block
      CALL get_index_range(subset, jb, startidx, endidx)
      min_value = values(startidx, jb)
      max_value = values(startidx, jb)
      sum_value = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb, startidx, endidx) FIRSTPRIVATE(min_in_block, max_in_block, sum_value) &
!ICON_OMP  reduction(min:min_value) reduction(max:max_value) reduction(sum:sum_value)
      DO jb = subset%start_block, subset%end_block
        CALL get_index_range(subset, jb, startidx, endidx)
        min_in_block = MINVAL(values(startidx:endidx, jb))
        max_in_block = MAXVAL(values(startidx:endidx, jb))
        min_value    = MIN(min_value, min_in_block)
        max_value    = MAX(max_value, max_in_block)
        sum_value    = sum_value + SUM(values(startidx:endidx, jb))
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      ! compute the total number of values
      IF ((subset%end_block - subset%start_block) > 1) THEN
        number_of_values = (subset%end_block - subset%start_block -1) * subset%block_size
      ELSE
        number_of_values = 0
      ENDIF
      number_of_values = number_of_values + subset%end_index + (subset%block_size - subset%start_index + 1)

    ELSE
      min_value = values(1,1)
      max_value = values(1,1)
      sum_value = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb, startidx, endidx) FIRSTPRIVATE(min_in_block, max_in_block, sum_value) &
!ICON_OMP  reduction(min:min_value) reduction(max:max_value) reduction(sum:sum_value)
      DO jb = 1,  SIZE(values, 2)
        min_in_block = MINVAL(values(:, jb))
        max_in_block = MAXVAL(values(:, jb))
        min_value    = MIN(min_value, min_in_block)
        max_value    = MAX(max_value, max_in_block)
        sum_value    = sum_value + SUM(values(:, jb))
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      number_of_values = SIZE(values)
    ENDIF
    

    ! the global min, max, mean, is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN
      communicator = get_my_mpi_work_communicator()
      minmaxmean(1) = p_min( min_value,  comm=communicator ) ! only mpi_all_reduce is avaliable
      minmaxmean(2) = p_max( max_value,  comm=communicator, root=process_mpi_stdio_id )

      ! these are avaliable to all processes
      global_number_of_values = p_sum( REAL(number_of_values,wp),  comm=communicator)
      minmaxmean(3) = p_sum( sum_value,  comm=communicator) / global_number_of_values

 !     write(0,*) "globalspace_2D_minmaxmean parallel:",  sum_value, number_of_values

    ELSE

      minmaxmean(1) = min_value
      minmaxmean(2) = max_value
      minmaxmean(3) = sum_value / REAL(number_of_values, wp)
  !    write(0,*) "globalspace_2D_minmaxmean seq:", min_value, max_value, sum_value, number_of_values

    ENDIF

  END FUNCTION globalspace_2D_minmaxmean
  !-----------------------------------------------------------------------
  

  !-----------------------------------------------------------------------
  !>
  FUNCTION globalspace_3D_minmaxmean(values, subset, start_level, end_level) result(minmaxmean)
    REAL(wp), INTENT(in) :: values(:,:,:)
    TYPE(t_subset_range), TARGET, OPTIONAL :: subset
    INTEGER, OPTIONAL :: start_level, end_level
    REAL(wp) :: minmaxmean(3)

    REAL(wp) :: min_in_block, max_in_block, min_value, max_value, sum_value, global_number_of_values
    INTEGER :: jb, level, startidx, endidx, start_vertical, end_vertical, number_of_values
    INTEGER :: communicator
!    INTEGER :: idx
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
      sum_value = 0

!ICON_OMP_PARALLEL_DO PRIVATE(jb, startidx, endidx) FIRSTPRIVATE(min_in_block, max_in_block, sum_value) &
!ICON_OMP  reduction(min:min_value) reduction(max:max_value) reduction(sum:sum_value)
      DO jb = subset%start_block, subset%end_block
        CALL get_index_range(subset, jb, startidx, endidx)
        DO level = start_vertical, end_vertical
!          DO idx = startidx, endidx
!            WRITE(0,*) "checking ", jb, level, idx
!            WRITE(0,*) "                  =   ", values(idx, level, jb)
!          ENDDO
          min_in_block = MINVAL(values(startidx:endidx, level, jb))
          max_in_block = MAXVAL(values(startidx:endidx, level, jb))
          min_value    = MIN(min_value, min_in_block)
          max_value    = MAX(max_value, max_in_block)
          sum_value    = sum_value + SUM(values(startidx:endidx, level, jb))
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      IF ((subset%end_block - subset%start_block) > 1) THEN
        number_of_values = (subset%end_block - subset%start_block -1) * subset%block_size
      ELSE
        number_of_values = 0
      ENDIF
      number_of_values = (number_of_values + subset%end_index + (subset%block_size - subset%start_index + 1)) * &
        & (end_vertical - start_vertical)

    ELSE

      min_value      = values(1, start_vertical, 1)
      max_value      = values(1, start_vertical, 1)
      sum_value = 0

!ICON_OMP_PARALLEL_DO PRIVATE(jb, startidx, endidx) FIRSTPRIVATE(min_in_block, max_in_block, sum_value) &
!ICON_OMP  reduction(min:min_value) reduction(max:max_value) reduction(sum:sum_value)
      DO jb = 1,  SIZE(values, 3)
        DO level = start_vertical, end_vertical
          min_in_block = MINVAL(values(:, level, jb))
          max_in_block = MAXVAL(values(:, level, jb))
          min_value    = MIN(min_value, min_in_block)
          max_value    = MAX(max_value, max_in_block)
          sum_value    = sum_value + SUM(values(:, level, jb))
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      number_of_values = SIZE(values(:,start_vertical:end_vertical,:))
    ENDIF

    ! the global min, max is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN
      communicator = get_my_mpi_work_communicator()
      minmaxmean(1) = p_min( min_value,  comm=communicator ) ! only mpi_all_reduce is avaliable
      minmaxmean(2) = p_max( max_value,  comm=communicator, root=process_mpi_stdio_id )

      ! these are avaliable to all processes
      global_number_of_values = p_sum( REAL(number_of_values,wp),  comm=communicator)
      minmaxmean(3) = p_sum( sum_value,  comm=communicator) / global_number_of_values

    ELSE

      minmaxmean(1) = min_value
      minmaxmean(2) = max_value
      minmaxmean(3) = sum_value / REAL(number_of_values, wp)

    ENDIF

  END FUNCTION globalspace_3D_minmaxmean
  !-----------------------------------------------------------------------



END MODULE mo_statistics_utils

