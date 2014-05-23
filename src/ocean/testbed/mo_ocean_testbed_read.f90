!>
!! @brief tests the mo_read_interface methods
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ocean_testbed_read

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_timer,               ONLY: init_timer, ltimer, new_timer, timer_start, timer_stop, &
    & activate_sync_timers, timers_level, timer_barrier, timer_radiaton_recv
  USE mo_parallel_config,     ONLY: nproma, icon_comm_method
  USE mo_ocean_nml
  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no
  USE mo_io_units,            ONLY: filename_max

  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_read_interface
  USE mo_netcdf_read
!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: ocean_test_read


CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !! Test reading T
  SUBROUTINE ocean_test_read(namelist_filename,shr_namelist_filename, patch_3d)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename
    TYPE(t_patch_3d ),TARGET, INTENT(in)    :: patch_3d

    REAL(wp), POINTER :: T(:,:,:,:), T_check(:,:,:,:)   ! is (nproma, levels, blocks, time )
    INTEGER :: levels, lnwl_size, return_status, stream_id
    TYPE(t_patch),POINTER            :: patch_2d
    CHARACTER(filename_max) :: OutputFileName   !< file name for reading in

    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_testbed_read:ocean_test_read"

    CALL message(method_name,   initialState_InputFileName)
    patch_2d => patch_3d%p_patch_2d(1)
    !---------------------------------------------------------------------
    CALL read_onCells_3D_time(                &
      & filename=initialState_InputFileName,  &
      & variable_name="T",                    &
      & return_pointer=T,                     &
      & patch=patch_2d,                       &
      & return_status=return_status )

    !---------------------------------------------------------------------
    ! check
    OutputFileName="testOut.nc"
    CALL message(method_name,   OutputFileName)
    return_status = netcdf_write_oncells_3d_time(  &
      & filename=OutputFileName,                   &
      & variable_name="T",                         &
      & write_array=T,                             &
      & patch=patch_2d)
    CALL read_onCells_3D_time(                &
      & filename=OutputFileName,              &
      & variable_name="T",                    &
      & return_pointer=T_check,               &
      & patch=patch_2d,                       &
      & return_status=return_status )
    IF ( MAXVAL(ABS(T - T_check )) > 0.0_wp ) &
      CALL finish(method_name, "Check failed")
    !---------------------------------------------------------------------
    DEALLOCATE(T)
    DEALLOCATE(T_check)

  END SUBROUTINE ocean_test_read
  !-------------------------------------------------------------------------


END MODULE mo_ocean_testbed_read

