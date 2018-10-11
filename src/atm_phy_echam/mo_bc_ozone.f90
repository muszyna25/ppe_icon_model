!>
!! @brief Subroutine read_bc_ozone reads monthly ozone
!! concentrations from yearly files. The ozone concentrations are
!! used as boundary conditions for the radiative forcing of the
!! atmosphereAMIP. The routine is called from mo_echam_phy_interface.f90.
!!
!! @author Sebastian Rast, MPI-M
!!
!! @par Revision History
!!  Original version March 2013
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_bc_ozone

  USE mo_kind,                     ONLY: wp, i8
  USE mo_exception,                ONLY: message, message_text
  USE mo_model_domain,             ONLY: t_patch
  USE mo_parallel_config,          ONLY: p_test_run
  USE mo_io_config,                ONLY: default_read_method
  USE mo_read_interface,           ONLY: nf, openInputFile, closeFile, &
  &                                      read_3D_time, t_stream_id, on_cells
  USE mo_mpi,                      ONLY: my_process_is_stdio, p_bcast, &
                                  &      p_comm_work_test, p_comm_work, p_io
  USE mo_physical_constants,       ONLY: amo3, amd
  USE mo_impl_constants,           ONLY: max_dom
  USE mo_grid_config,              ONLY: n_dom
  USE mo_echam_rad_config,         ONLY: echam_rad_config

  IMPLICIT NONE
  PRIVATE
  REAL(wp), PARAMETER               :: vmr2mmr_o3=amo3/amd  ! Volume mixing ratio to mass mixing ratio
  INTEGER(i8), SAVE                 :: pre_year(max_dom)=-999999 ! Variable to check if it is time to read

  PUBLIC                            :: ext_ozone
  PUBLIC                            :: read_bc_ozone

  TYPE t_ext_ozone
    REAL(wp), ALLOCATABLE           :: o3_plev(:,:,:,:)     ! Monthly ozone mass mixing ratio at pressure levels
    REAL(wp), ALLOCATABLE           :: plev_full_o3(:)      ! Full pressure levels
    REAL(wp), ALLOCATABLE           :: plev_half_o3(:)      ! Half pressure levels, derived from plev_full_o3
    INTEGER                         :: nplev_o3             ! Number of full pressure levels
  END TYPE t_ext_ozone

  TYPE(t_ext_ozone), ALLOCATABLE, TARGET :: ext_ozone(:)

  INCLUDE 'netcdf.inc'
  
CONTAINS

  SUBROUTINE read_bc_ozone(year, p_patch)

    INTEGER(i8)  , INTENT(in)         :: year
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch

    CHARACTER(len=512)                :: fname
    TYPE(t_stream_id)                 :: stream_id
    CHARACTER(len=4)                  :: cyear
    CHARACTER(len=2)                  :: cjg
    CHARACTER(len=*), PARAMETER       :: subprog_name &
         = 'mo_bc_ozone:read_bc_ozone'
    INTEGER                           :: ncid, varid, mpi_comm, jg
    INTEGER                           :: nplev_o3
    REAL(wp), POINTER                 :: zo3_plev(:,:,:,:)           ! (nproma, levels, blocks, time)

    jg    = p_patch%id
    WRITE(cjg,'(i2.2)') jg

    ! allocate once only structure for all grids
    IF (.NOT. ALLOCATED(ext_ozone)) ALLOCATE(ext_ozone(n_dom))

    IF (year > pre_year(jg)) THEN
      !
      ! If year = pre_year, then the external monthly ozone data are already stored.
      ! Nothing needs to be done.
      !
      IF (ALLOCATED(ext_ozone(jg)% o3_plev)) THEN
        !
        ! If the memory for o3_plev is already allocated, then this is a new
        ! year of a simulation that started before this year, meaning that
        ! a year of external monthly ozone data is already stored.
        !
        IF (echam_rad_config(jg)% irad_o3==8) THEN
          !
          ! If irad_o3=8, ozone is transient and the external monthly ozone data
          ! of this year and January of the next year  must be read from file.
          !
          ! For other irad_o3 cases no new data must be read.
          ! Nothing needs to be done.
          !
          WRITE(message_text,'(2a)') 'Copy ozone for months 12:13 to months 0:1'
          CALL message('read_bc_ozone', message_text)
          ext_ozone(jg)% o3_plev(:,:,:,0:1) = ext_ozone(jg)% o3_plev(:,:,:,12:13)

          WRITE(cyear,'(i4)') year
          IF (n_dom > 1) THEN
            fname = 'bc_ozone_DOM'//TRIM(cjg)//'_'//TRIM(cyear)//'.nc'
          ELSE
            fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
          END IF

          WRITE(message_text,'(2a)') 'Read ozone for months  2:12 from file: ',TRIM(fname)
          CALL message('read_bc_ozone', message_text)
          stream_id = openInputFile(fname, p_patch, default_read_method)
          CALL read_3D_time(stream_id=stream_id, location=on_cells,      &
            &               variable_name='O3', return_pointer=zo3_plev, &
            &               start_timestep=2,end_timestep=12)
          CALL closeFile(stream_id)
          ext_ozone(jg)% o3_plev(:,:,:,2:12) = vmr2mmr_o3*zo3_plev(:,:,:,1:11)

          WRITE(cyear,'(i4)') year+1
          IF (n_dom > 1) THEN
            fname = 'bc_ozone_DOM'//TRIM(cjg)//'_'//TRIM(cyear)//'.nc'
          ELSE
            fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
          END IF

          WRITE(message_text,'(2a)') 'Read ozone for month     13 from file: ',TRIM(fname)
          CALL message('read_bc_ozone', message_text)
          stream_id = openInputFile(fname, p_patch, default_read_method)
          CALL read_3D_time(stream_id=stream_id, location=on_cells,      &
            &               variable_name='O3', return_pointer=zo3_plev, &
            &               start_timestep=1,end_timestep=1)
          CALL closeFile(stream_id)
          ext_ozone(jg)% o3_plev(:,:,:,13)=vmr2mmr_o3*zo3_plev(:,:,:,1)

        END IF
        !
      ELSE
        !
        ! If the memory of o3_plev is not allocated, then this is
        ! a new simulation or a simulation that starts from a restart file
        ! and ozone data must be read.
        !
        ! Depending on the irad_o3 case different amounts of data must be read.
        !
        SELECT CASE (echam_rad_config(jg)% irad_o3)
          !
        CASE (2) ! Ozone has a climatological annual cycle defined by monthly data in an annual file
          !
          IF (n_dom > 1) THEN
            fname = 'bc_ozone_DOM'//TRIM(cjg)//'.nc'
          ELSE
            fname = 'bc_ozone'//'.nc'
          END IF
          !
          WRITE(message_text,'(2a)') 'Read clim. annual cycle of ozone from file: ',TRIM(fname)
          CALL message('read_bc_ozone', message_text)
          !
          stream_id = openInputFile(fname, p_patch, default_read_method)
          CALL read_3D_time(stream_id=stream_id, location=on_cells,      &
            &               variable_name='O3', return_pointer=zo3_plev, &
            &               start_timestep=1,end_timestep=12)
          CALL closeFile(stream_id)
          !
          ! Now the spatial dimensions are known --> allocate memory for months 0:13
          ALLOCATE(ext_ozone(jg)% o3_plev(SIZE(zo3_plev,1),SIZE(zo3_plev,2),SIZE(zo3_plev,3),0:13))
          !
          ext_ozone(jg)% o3_plev(:,:,:,   0)=vmr2mmr_o3*zo3_plev(:,:,:,  12)
          ext_ozone(jg)% o3_plev(:,:,:,1:12)=vmr2mmr_o3*zo3_plev(:,:,:,1:12)
          ext_ozone(jg)% o3_plev(:,:,:,  13)=vmr2mmr_o3*zo3_plev(:,:,:,   1)
          !
          !
        CASE (4) ! Ozone is constant in time
          !
          IF (n_dom > 1) THEN
            fname = 'bc_ozone_DOM'//TRIM(cjg)//'.nc'
          ELSE
            fname = 'bc_ozone'//'.nc'
          END IF
          !
          WRITE(message_text,'(2a)') 'Read constant-in-time ozone from file: ',TRIM(fname)
          CALL message('read_bc_ozone', message_text)
          !
          stream_id = openInputFile(fname, p_patch, default_read_method)
          CALL read_3D_time(stream_id=stream_id, location=on_cells,      &
            &               variable_name='O3', return_pointer=zo3_plev, &
            &               start_timestep=1,end_timestep=1)
          CALL closeFile(stream_id)
          !
          ! Now the spatial dimensions are known --> allocate memory for one time slice
          ALLOCATE(ext_ozone(jg)% o3_plev(SIZE(zo3_plev,1),SIZE(zo3_plev,2),SIZE(zo3_plev,3),1))
          !
          ext_ozone(jg)% o3_plev(:,:,:,1)=vmr2mmr_o3*zo3_plev(:,:,:,1)
          !
          !
        CASE (8) ! Ozone is transient and defined by monthly data in annual files
           !
           ! 1. Read December of the previous year
           !
           WRITE(cyear,'(i4)') year-1
           IF (n_dom > 1) THEN
             fname = 'bc_ozone_DOM'//TRIM(cjg)//'_'//TRIM(cyear)//'.nc'
           ELSE
             fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
           END IF
           !
           WRITE(message_text,'(2a)') 'Read ozone for month      0 from file: ',TRIM(fname)
           CALL message('read_bc_ozone', message_text)
           stream_id = openInputFile(fname, p_patch, default_read_method)
           CALL read_3D_time(stream_id=stream_id, location=on_cells,         &
                &               variable_name='O3', return_pointer=zo3_plev, &
                &               start_timestep=12,end_timestep=12)
           CALL closeFile(stream_id)
           !
           ! Now the spatial dimensions are known --> allocate memory for months 0:13
           ALLOCATE(ext_ozone(jg)% o3_plev(SIZE(zo3_plev,1),SIZE(zo3_plev,2),SIZE(zo3_plev,3),0:13))
           !
           ext_ozone(jg)% o3_plev(:,:,:,0)=vmr2mmr_o3*zo3_plev(:,:,:,1)
           !
           ! 2. Read January-December of this year
           !
           WRITE(cyear,'(i4)') year
           IF (n_dom > 1) THEN
             fname = 'bc_ozone_DOM'//TRIM(cjg)//'_'//TRIM(cyear)//'.nc'
           ELSE
             fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
           END IF
           !
           WRITE(message_text,'(2a)') 'Read ozone for months  1:12 from file: ',TRIM(fname)
           CALL message('read_bc_ozone', message_text)
           stream_id = openInputFile(fname, p_patch, default_read_method)
           CALL read_3D_time(stream_id=stream_id, location=on_cells,         &
                &               variable_name='O3', return_pointer=zo3_plev, &
                &               start_timestep=1,end_timestep=12)
           CALL closeFile(stream_id)
           !
           ext_ozone(jg)% o3_plev(:,:,:,1:12)=vmr2mmr_o3*zo3_plev(:,:,:,1:12)
           !
           ! 3. Read January of the next year
           !
           WRITE(cyear,'(i4)') year+1
           IF (n_dom > 1) THEN
             fname = 'bc_ozone_DOM'//TRIM(cjg)//'_'//TRIM(cyear)//'.nc'
           ELSE
             fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
           END IF
           !
           WRITE(message_text,'(2a)') 'Read ozone for month     13 from file: ',TRIM(fname)
           CALL message('read_bc_ozone', message_text)
           stream_id = openInputFile(fname, p_patch, default_read_method)
           CALL read_3D_time(stream_id=stream_id, location=on_cells,         &
                &               variable_name='O3', return_pointer=zo3_plev, &
                &               start_timestep=1,end_timestep=1)
           CALL closeFile(stream_id)
           !
           ext_ozone(jg)% o3_plev(:,:,:,13)=vmr2mmr_o3*zo3_plev(:,:,:,1)
        CASE(10)
                       !
           ! 1. Read December of the previous year
           !
           WRITE(cyear,'(i4)') year-1
           IF (jg > 1) THEN
             fname = 'bc_ozone_'//TRIM(cyear)//'_DOM'//TRIM(cjg)//'.nc'
           ELSE
             fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
           END IF
           !
           WRITE(message_text,'(2a)') 'Read ozone for month      0 from file: ',TRIM(fname)
           CALL message('read_bc_ozone', message_text)
           stream_id = openInputFile(fname, p_patch, default_read_method)
           CALL read_3D_time(stream_id=stream_id, location=on_cells,         &
                &               variable_name='O3', return_pointer=zo3_plev, &
                &               start_timestep=12,end_timestep=12)
           CALL closeFile(stream_id)
           !
           ! Now the spatial dimensions are known --> allocate memory for months 0:13
           ALLOCATE(ext_ozone(jg)% o3_plev(SIZE(zo3_plev,1),SIZE(zo3_plev,2),SIZE(zo3_plev,3),0:13))
           !
           ext_ozone(jg)% o3_plev(:,:,:,0)=vmr2mmr_o3*zo3_plev(:,:,:,1)
           !
           ! 2. Read January-December of this year
           !
           WRITE(cyear,'(i4)') year
           IF (jg > 1) THEN
             fname = 'bc_ozone_'//TRIM(cyear)//'_DOM'//TRIM(cjg)//'.nc'
           ELSE
             fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
           END IF
           !
           WRITE(message_text,'(2a)') 'Read ozone for months  1:12 from file: ',TRIM(fname)
           CALL message('read_bc_ozone', message_text)
           stream_id = openInputFile(fname, p_patch, default_read_method)
           CALL read_3D_time(stream_id=stream_id, location=on_cells,         &
                &               variable_name='O3', return_pointer=zo3_plev, &
                &               start_timestep=1,end_timestep=12)
           CALL closeFile(stream_id)
           !
           ext_ozone(jg)% o3_plev(:,:,:,1:12)=vmr2mmr_o3*zo3_plev(:,:,:,1:12)
           !
           ! 3. Read January of the next year
           !
           WRITE(cyear,'(i4)') year+1
           IF (jg > 1) THEN
             fname = 'bc_ozone_'//TRIM(cyear)//'_DOM'//TRIM(cjg)//'.nc'
           ELSE
             fname = 'bc_ozone_'//TRIM(cyear)//'.nc'
           END IF
           !
           WRITE(message_text,'(2a)') 'Read ozone for month     13 from file: ',TRIM(fname)
           CALL message('read_bc_ozone', message_text)
           stream_id = openInputFile(fname, p_patch, default_read_method)
           CALL read_3D_time(stream_id=stream_id, location=on_cells,         &
                &               variable_name='O3', return_pointer=zo3_plev, &
                &               start_timestep=1,end_timestep=1)
           CALL closeFile(stream_id)
           !
           ext_ozone(jg)% o3_plev(:,:,:,13)=vmr2mmr_o3*zo3_plev(:,:,:,1)
           !
        END SELECT

        ! Read pressure level grid on which the external ozone data are valid
        ! This needs to be done equally for all irad_o3 cases.

        nplev_o3 = SIZE(ext_ozone(jg)% o3_plev,2)

        ext_ozone(jg)% nplev_o3 = nplev_o3

        IF(ALLOCATED(ext_ozone(jg)% plev_full_o3)) DEALLOCATE(ext_ozone(jg)% plev_full_o3)
        IF(ALLOCATED(ext_ozone(jg)% plev_half_o3)) DEALLOCATE(ext_ozone(jg)% plev_half_o3)
        ALLOCATE(ext_ozone(jg)% plev_full_o3(nplev_o3  ))
        ALLOCATE(ext_ozone(jg)% plev_half_o3(nplev_o3+1))

        mpi_comm = MERGE(p_comm_work_test, p_comm_work, p_test_run)

        IF(my_process_is_stdio()) THEN
          CALL nf(nf_open(TRIM(fname), NF_NOWRITE, ncid), subprog_name)
          CALL nf(nf_inq_varid(ncid, 'plev', varid), subprog_name)
          CALL nf(nf_get_var_double(ncid, varid, ext_ozone(jg)% plev_full_o3), subprog_name)
          CALL nf(nf_close(ncid), subprog_name)
        END IF
        CALL p_bcast(ext_ozone(jg)% plev_full_o3(:), p_io, mpi_comm)

        ! define half levels of ozone pressure grid
        ! upper boundary: ph =      0.Pa -> extrapolation of uppermost value
        ! lower boundary: ph = 125000.Pa -> extrapolation of lowermost value
        ext_ozone(jg)% plev_half_o3(1)          = 0._wp
        ext_ozone(jg)% plev_half_o3(2:nplev_o3) = 0.5_wp*( ext_ozone(jg)% plev_full_o3(1:nplev_o3-1) &
          &                                               +ext_ozone(jg)% plev_full_o3(2:nplev_o3  ) )
        ext_ozone(jg)% plev_half_o3(nplev_o3+1) = 125000._wp

      END IF

      pre_year(jg) = year

    END IF

  END SUBROUTINE read_bc_ozone

END MODULE mo_bc_ozone
