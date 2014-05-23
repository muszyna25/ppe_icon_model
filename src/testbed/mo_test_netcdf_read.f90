!>
!! @brief tests the mo_netcdf_read methods
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
MODULE mo_test_netcdf_read

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_timer,               ONLY: ltimer, activate_sync_timers, timers_level
  USE mo_parallel_config,     ONLY: nproma

  USE mo_icon_testbed_config, ONLY: testfile_3D_time, testfile_2D_time

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
  !! Test reading amip aerosol
  !! tesetd with /pool/data/ICON/import/incoming/bc_aeropt_kinne/bc_aeropt_kinne_sw_b14_fin_2000.nc
  SUBROUTINE test_netcdf_read(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_patch), POINTER   :: patch
    REAL(wp), POINTER :: lnwl_array(:), levels_array(:), times_array(:)
    REAL(wp), POINTER :: aod(:,:,:,:), asy(:,:,:,:)   ! is (nproma, lnwl, blocks, time (months) )
    REAL(wp), POINTER :: z_aer_fine_mo(:,:,:,:), return_pointer(:,:,:,:)
    INTEGER :: levels, lnwl_size, stream_id

    CHARACTER(*), PARAMETER :: method_name = "mo_test_netcdf_read:test_netcdf_read"

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
    CALL message(method_name,   testfile_3D_time(1))
    stream_id = netcdf_open_input(filename = testfile_3D_time(1))

    !---------------------------------------------------------------------
    lnwl_array => netcdf_read_1D(     &
      & file_id       = stream_id,    &
      & variable_name = "lnwl")

    levels_array => netcdf_read_1D(     &
      & file_id       = stream_id,    &
      & variable_name = "lev")

    times_array => netcdf_read_1D(     &
      & file_id       = stream_id,    &
      & variable_name = "time")

    !-----------------------------------------------------
    ! example with non-allocated arrays
    ! mote that allocation time dim will be 1:12, as in the file
    aod => netcdf_read_oncells_3D_time( &
      & file_id       = stream_id,      &
      & variable_name = "aod",          &
      & levelsdim_name = "lnwl",        & ! optional, just for checking
      & patch         = patch)

    !-----------------------------------------------------
    ! example with allocated arrays, time dim is 0:13
    lnwl_size = SIZE(lnwl_array, 1)
    ALLOCATE(asy(nproma, lnwl_size, patch%nblks_c, 0:13 ))
    ! read previous month, should be filled form the previous read though
    return_pointer => netcdf_read_oncells_3D_time( &
      & filename       = testfile_3D_time(1),      &
      & variable_name  = "asy",                    &
      & fill_array     = asy(:,:,:,0:0),           &
      & patch          = patch,                    &
      & levelsdim_name = "lnwl",                   & ! optional, just for checking
      & start_timestep = 12,                       &
      & end_timestep   = 12)
    ! read current year
    return_pointer => netcdf_read_oncells_3D_time( &
      & file_id        = stream_id,                &
      & variable_name  = "asy",                    &
      & fill_array     = asy(:,:,:,1:12),          &
      & patch          = patch,                    &
      & levelsdim_name = "lnwl")                    ! optional, just for checking
    ! read next month
    return_pointer => netcdf_read_oncells_3D_time( &
      & filename       = testfile_3D_time(1),      &
      & variable_name  = "asy",                    &
      & fill_array     = asy(:,:,:,13:13),         &
      & patch          = patch,                    &
      & start_timestep = 1,                        &
      & end_timestep   = 1)

    !-----------------------------------------------------
    ! check Z-axis vertical levels
    IF (SIZE(levels_array, 1) /= levels) THEN
      write(0,*) "patch levels=", levels, " SIZE(levels_array)=", SIZE(levels_array, 1)
      CALL finish(method_name, "SIZE(levels_array, 1) /= levels")
    ENDIF
    ALLOCATE(z_aer_fine_mo(nproma, levels, patch%nblks_c, 0:13 ))
    ! read previous month, should be filled form the previous read though
    ! read current year
    return_pointer => netcdf_read_oncells_3D_time(  &
      & file_id        = stream_id,                 &
      & variable_name  = "z_aer_fine_mo",           &
      & fill_array     = z_aer_fine_mo(:,:,:,1:12), &
      & patch          = patch,                     &
      & levelsdim_name = "lev")                    ! optional, just for checking

    !---------------------------------------------------------------------
    DEALLOCATE(lnwl_array, levels_array, times_array, aod, asy, z_aer_fine_mo)

    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    CALL destruct_atmo_hydrostatic()
    CALL destruct_atmo_model()


  END SUBROUTINE test_netcdf_read
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_netcdf_read_demo_1(namelist_filename,shr_namelist_filename)

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
     

  END SUBROUTINE test_netcdf_read_demo_1
  !-------------------------------------------------------------------------
        

END MODULE mo_test_netcdf_read

