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
  USE mo_model_domain,             ONLY: t_patch
  USE mo_parallel_config,          ONLY: p_test_run
  USE mo_io_config,                ONLY: default_read_method
  USE mo_read_interface,           ONLY: nf, openInputFile, closeFile, &
  &                                      read_3D_time, t_stream_id, on_cells
  USE mo_mpi,                      ONLY: my_process_is_stdio, p_bcast, &
                                  &      p_comm_work_test, p_comm_work, p_io
  USE mo_physical_constants,       ONLY: amo3, amd
  USE mo_echam_phy_config,         ONLY: echam_phy_config

  IMPLICIT NONE
  PRIVATE
  REAL(wp), PARAMETER               :: vmr2mmr_o3=amo3/amd               ! Volume mixing ratio to mass mixing ratio
  INTEGER(i8), SAVE                 :: pre_year=-999999                  ! Variable to check if it is time to read

  PUBLIC                            :: o3_plev, plev_half_o3, plev_full_o3, nplev_o3
  PUBLIC                            :: read_bc_ozone

  REAL(wp), ALLOCATABLE             :: o3_plev(:,:,:,:)                  ! Ozonemass mixing ratio at pressure levels
  REAL(wp), ALLOCATABLE             :: plev_half_o3(:), plev_full_o3(:)  ! Pressure levels in ozone file
  INTEGER                           :: nplev_o3

  INCLUDE 'netcdf.inc'
  
CONTAINS

  SUBROUTINE read_bc_ozone(year, p_patch)

    INTEGER(i8)  , INTENT(in)         :: year
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch

    CHARACTER(len=16)                 :: fname
    TYPE(t_stream_id)                 :: stream_id
    CHARACTER(len=4)                  :: cyear
    CHARACTER(len=*), PARAMETER       :: subprog_name &
         = 'mo_bc_ozone:read_bc_ozone'

    INTEGER                           :: ncid, varid, mpi_comm
    REAL(wp), POINTER                 :: zo3_plev(:,:,:,:)

    IF (year > pre_year) THEN

      IF (ALLOCATED(o3_plev)) THEN
        o3_plev(:,:,:,0:1)=o3_plev(:,:,:,12:13)

        IF ( echam_phy_config(p_patch%id)%lamip ) THEN
          WRITE(cyear,'(i4)') year
          fname='bc_ozone_'//TRIM(cyear)//'.nc'
        ELSE
          fname='bc_ozone'//'.nc'
        ENDIF

        write(0,*) 'Read ozone from file: ',fname 
        stream_id = openInputFile(fname, p_patch, default_read_method)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, &
          &               variable_name='O3', return_pointer=zo3_plev, &
          &               start_timestep=2,end_timestep=12)
        CALL closeFile(stream_id)
        o3_plev(:,:,:,2:12)=vmr2mmr_o3*zo3_plev(:,:,:,1:11)
        write(cyear,'(i4)') year+1

        IF ( echam_phy_config(p_patch%id)%lamip ) THEN
          fname='bc_ozone_'//TRIM(cyear)//'.nc'
        ELSE
          fname='bc_ozone'//'.nc'
        ENDIF

        stream_id = openInputFile(fname, p_patch, default_read_method)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, &
          &               variable_name='O3', return_pointer=zo3_plev, &
          &               start_timestep=1,end_timestep=1)
        CALL closeFile(stream_id)
        o3_plev(:,:,:,13)=vmr2mmr_o3*zo3_plev(:,:,:,1)
      ELSE

        IF ( echam_phy_config(p_patch%id)%lamip ) THEN
          WRITE(cyear,'(i4)') year
          fname='bc_ozone_'//TRIM(cyear)//'.nc'
        ELSE
          fname='bc_ozone'//'.nc'
        ENDIF

        write(0,*) 'Read ozone from file: ',fname 
        stream_id = openInputFile(fname, p_patch, default_read_method)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, &
          &               variable_name='O3', return_pointer=zo3_plev)
        CALL closeFile(stream_id)
        ALLOCATE(o3_plev(SIZE(zo3_plev,1),SIZE(zo3_plev,2),SIZE(zo3_plev,3),0:13))
        o3_plev(:,:,:,1:12)=vmr2mmr_o3*zo3_plev

        IF ( echam_phy_config(p_patch%id)%lamip ) THEN
          WRITE(cyear,'(i4)') year-1
          fname='bc_ozone_'//TRIM(cyear)//'.nc'
        ELSE
          fname='bc_ozone'//'.nc'
        ENDIF

        stream_id = openInputFile(fname, p_patch, default_read_method)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, &
          &               variable_name='O3', return_pointer=zo3_plev, &
          &               start_timestep=12,end_timestep=12)
        CALL closeFile(stream_id)
        o3_plev(:,:,:,0)=vmr2mmr_o3*zo3_plev(:,:,:,1)

        IF ( echam_phy_config(p_patch%id)%lamip ) THEN
          WRITE(cyear,'(i4)') year+1
          fname='bc_ozone_'//TRIM(cyear)//'.nc'
        ELSE
          fname='bc_ozone'//'.nc'
        ENDIF

        stream_id = openInputFile(fname, p_patch, default_read_method)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, &
          &               variable_name='O3', return_pointer=zo3_plev, &
          &               start_timestep=1,end_timestep=1)
        CALL closeFile(stream_id)
        o3_plev(:,:,:,13)=vmr2mmr_o3*zo3_plev(:,:,:,1)       
      END IF

      nplev_o3=SIZE(o3_plev,2)

      IF(ALLOCATED(plev_full_o3)) DEALLOCATE(plev_full_o3)
      IF(ALLOCATED(plev_half_o3)) DEALLOCATE(plev_half_o3)
      ALLOCATE(plev_full_o3(nplev_o3))
      ALLOCATE(plev_half_o3(nplev_o3+1))

      mpi_comm = MERGE(p_comm_work_test, p_comm_work, p_test_run)

      IF(my_process_is_stdio()) THEN
        CALL nf(nf_open(TRIM(fname), NF_NOWRITE, ncid), subprog_name)
        CALL nf(nf_inq_varid(ncid, 'plev', varid), subprog_name)
        CALL nf(nf_get_var_double(ncid, varid, plev_full_o3), subprog_name)
        CALL nf(nf_close(ncid), subprog_name)
      END IF
      CALL p_bcast(plev_full_o3(:), p_io, mpi_comm)

      ! define half levels of ozone pressure grid
      ! upper boundary: ph =      0.Pa -> extrapolation of uppermost value
      ! lower boundary: ph = 125000.Pa -> extrapolation of lowermost value
      plev_half_o3(1)=0._wp
      plev_half_o3(2:nplev_o3)=0.5_wp*(plev_full_o3(1:nplev_o3-1)+plev_full_o3(2:nplev_o3))
      plev_half_o3(nplev_o3+1)=125000._wp

      pre_year = year

    END IF

  END SUBROUTINE read_bc_ozone

  !>
  !! Subroutine to calculate interpolation indices and weights
  !! Given: Monthly mean values valid at the middle of a month for months 0 to 13
  !! The December of the previous year has index 0, January of next year, index 13
  !! inm1, inm2 are the indices of the months between the middles of which the event_days
  !! is located. Linear interpolation between these values with weights wgt1 and wgt2 resp.
  !! 
  !! @par Revision History
  !! Rewritten after echam6 by J.S. Rast, MPI-M Hamburg (2013-07-31)
  !!
END MODULE mo_bc_ozone
