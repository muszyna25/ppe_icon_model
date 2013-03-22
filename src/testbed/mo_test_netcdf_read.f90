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

  USE mo_model_domain,        ONLY: t_patch, p_patch
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

    TYPE(t_patch), POINTER   :: patch
    REAL(wp), POINTER :: fill_3D_time_array(:,:,:,:)   ! is (nproma, vertical_levels, blocks, time (months))
    REAL(wp), POINTER :: fill_2D_time_array(:,:,:),  point_2D_time(:,:,:)   ! is (nproma, blocks, time )
    INTEGER :: levels, return_status

    CHARACTER(*), PARAMETER :: method_name = "test_netcdf_read"

    !---------------------------------------------------------------------
    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    CALL construct_atmo_hydrostatic()
    patch => p_patch(1)
    levels  = patch%nlev
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! test sst
    CALL message(method_name,   testfile_2D_time(1))

    !---------------------------------------------------------------------
    ! will read all timesteps in the file
    fill_2D_time_array => netcdf_read_oncells_2D_time(   &
      & filename=testfile_2D_time(1),              &
      & variable_name=TRIM(testfile_2D_time(2)),   &
      & patch=patch)

    ! here the fill_2D_time_array is allocated and filled
    ! write the sst array to output ONLY for comparing
    return_status = netcdf_write_oncells_2d_time("testfile_2D_time_alltimes.nc", TRIM(testfile_2D_time(2)), &
      & fill_2D_time_array, patch)
    DEALLOCATE(fill_2D_time_array)
    !---------------------------------------------------------------------

    ! Example of reading selective timesteps with explicit shape
    ! say 13 months
    ALLOCATE(fill_2D_time_array(nproma, patch%nblks_c, 0:13 ))
    point_2D_time => netcdf_read_oncells_2D_time(   &
      & filename=testfile_2D_time(1),               &
      & variable_name=TRIM(testfile_2D_time(2)),    &
      & fill_array=fill_2D_time_array(:,:,0:),      &
      & patch=patch,                                &
      & start_timestep=12, end_timestep=12)
    point_2D_time => netcdf_read_oncells_2D_time(   &
      & filename=testfile_2D_time(1),               &
      & variable_name=TRIM(testfile_2D_time(2)),    &
      & fill_array=fill_2D_time_array(:,:,1:),      &
      & patch=patch,                                &
      & start_timestep=1,  end_timestep=12)
    point_2D_time => netcdf_read_oncells_2D_time(   &
      & filename=testfile_2D_time(1),               &
      & variable_name=TRIM(testfile_2D_time(2)),    &
      & fill_array=fill_2D_time_array(:,:,13:),     &
      & patch=patch,                                &
      & start_timestep=1, end_timestep=1)
    ! write(0,*) "shape of ", TRIM(testfile_2D_time(2)), ":", SHAPE(fill_2D_time_array)
    !---------------------------------------------------------------------
    ! write the sst array to output for comparing
    return_status = netcdf_write_oncells_2d_time("testfile_2D_14times.nc", TRIM(testfile_2D_time(2)), &
      & fill_2D_time_array, patch)
    DEALLOCATE(fill_2D_time_array)
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! test O3
    CALL message(method_name,   testfile_3D_time(1))
    fill_3D_time_array => netcdf_read_oncells_3D_time(   &
      & filename=testfile_3D_time(1),              &
      & variable_name=TRIM(testfile_3D_time(2)),   &
      & patch=patch)
    !---------------------------------------------------------------------
    ! write the o3 array to output for comparing
    CALL message(method_name,   "write testfile_3D_time(1)")
    return_status = netcdf_write_oncells_3d_time(  &
      & filename=testfile_3D_time(1),              &
      & variable_name=TRIM(testfile_3D_time(2)),   &
      & write_array=fill_3D_time_array,            &
      & patch=patch)
    DEALLOCATE(fill_3D_time_array)
    !---------------------------------------------------------------------
!    NULLIFY(fill_3D_time_array)
!    ALLOCATE(fill_3D_time_array(nproma, levels,  patch%nblks_c, 0:0 ))


    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    CALL destruct_atmo_hydrostatic()
    CALL destruct_atmo_model()
     

  END SUBROUTINE test_netcdf_read
  !-------------------------------------------------------------------------
        

END MODULE mo_test_netcdf_read

