!>
!! @brief tests the mo_netcdf_read methods
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
!!
MODULE mo_test_netcdf_read

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_timer,               ONLY: init_timer, ltimer, new_timer, timer_start, timer_stop, &
    & print_timer, activate_sync_timers, timers_level, timer_barrier, timer_radiaton_recv
  USE mo_parallel_config,     ONLY: nproma, icon_comm_method

  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no
  USE mo_icon_testbed_config, ONLY: testbed_model, testfile_3D_time, testfile_2D_time

  USE mo_model_domain,        ONLY: p_patch
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model
  USE mo_atmo_hydrostatic,    ONLY: construct_atmo_hydrostatic, destruct_atmo_hydrostatic
  USE mo_netcdf_read

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: test_netcdf_read


CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_netcdf_read(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename
    
    REAL(wp), POINTER :: fill_3D_time_array(:,:,:,:)   ! is (nproma, vertical_levels, blocks, time (months))
    REAL(wp), POINTER :: fill_2D_time_array(:,:,:)    ! is (nproma, blocks, time )
    INTEGER :: return_status

    CHARACTER(*), PARAMETER :: method_name = "test_netcdf_read"

    !---------------------------------------------------------------------
    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    CALL construct_atmo_hydrostatic()
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! test sst
    ! read sst
    CALL message(method_name,   testfile_2D_time(1))
    NULLIFY(fill_2D_time_array)
    return_status = netcdf_read_oncells_2D_time(testfile_2D_time(1), TRIM(testfile_2D_time(2)), fill_2D_time_array, p_patch(1), &
      & start_timestep=2, end_timestep=2)
    ! write(0,*) "shape of ", TRIM(testfile_2D_time(2)), ":", SHAPE(fill_2D_time_array)
    !---------------------------------------------------------------------
    ! write the sst array to output for comparing
    return_status = netcdf_write_oncells_2d_time("testfile_2D_time.nc", TRIM(testfile_2D_time(2)), fill_2D_time_array, p_patch(1))
    DEALLOCATE(fill_2D_time_array)
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! test O3
    ! read o3
    CALL message(method_name,   testfile_3D_time(1))
    NULLIFY(fill_3D_time_array)
    return_status = netcdf_read_oncells_3D_time(testfile_3D_time(1), testfile_3D_time(2), fill_3D_time_array, p_patch(1))
    ! write(0,*) "shape of ", TRIM(testfile_3D_time(2)), ":", SHAPE(fill_3D_time_array)
    !---------------------------------------------------------------------
    ! write the o3 array to output for comparing
    return_status = netcdf_write_oncells_3d_time("testfile_3D_time.nc", TRIM(testfile_3D_time(2)), fill_3D_time_array, p_patch(1))
    DEALLOCATE(fill_3D_time_array)
    !---------------------------------------------------------------------


    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    CALL destruct_atmo_hydrostatic()
    CALL destruct_atmo_model()
     

  END SUBROUTINE test_netcdf_read
  !-------------------------------------------------------------------------
        

END MODULE mo_test_netcdf_read

