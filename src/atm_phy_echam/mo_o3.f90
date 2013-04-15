!>
!! @brief Subroutine read_amip_o3 reads ozone from files for the AMIP 
!! experiment. The routine is called from mo_echam_phy_interface.f90.
!!
!! @author Sebastian Rast, MPI-M
!!
!! @par Revision History
!!  Original version March 2013
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
!!
MODULE mo_o3

  USE mo_kind,                     ONLY: wp
  USE mo_model_domain,             ONLY: t_patch
  USE mo_parallel_config,          ONLY: nproma, p_test_run
  USE mo_run_config,               ONLY: nlev
  USE mo_netcdf_read,              ONLY: netcdf_read_oncells_3D_time, nf
  USE mo_mpi,                      ONLY: my_process_is_stdio, p_bcast, &
                                  &      p_comm_work_test, p_comm_work, p_io
  USE mo_physical_constants,       ONLY: amo3, amd

  IMPLICIT NONE
  PRIVATE
  REAL(wp), PARAMETER               :: vmr2mmr_o3=amo3/amd               ! Volume mixing ratio to mass mixing ratio
  INTEGER, SAVE                     :: pre_year=-999999                  ! Variable to check if it is time to read

  PUBLIC                            :: o3_plev, plev_half_o3, plev_full_o3, nplev_o3
  PUBLIC                            :: read_amip_o3

  REAL(wp), POINTER                 :: o3_plev(:,:,:,:)                  ! Ozone at pressure levels
  REAL(wp), ALLOCATABLE             :: plev_half_o3(:), plev_full_o3(:)  ! Pressure levels in ozone file
  INTEGER                           :: nplev_o3

  include 'netcdf.inc'

CONTAINS

  SUBROUTINE read_amip_o3(year, p_patch)

    INTEGER, INTENT(in)               :: year
    TYPE(t_patch), INTENT(in)         :: p_patch
    CHARACTER(len=12)                 :: fname
    CHARACTER(len=4)                  :: cyear
    CHARACTER(len=18)                 :: subprog_name
    INTEGER                           :: ncid, varid, mpi_comm,i
  
    IF (year > pre_year) THEN

      write(cyear,'(i4)') year
      fname='ozone'//TRIM(cyear)//'.nc'
      write(0,*) 'Read ozone from file: ',fname 
      o3_plev=>netcdf_read_oncells_3D_time(filename=fname,variable_name='O3',patch=p_patch)
      o3_plev=vmr2mmr_o3*o3_plev

      subprog_name='mo_o3:read_amip_o3'
      nplev_o3=SIZE(o3_plev,2)

      ALLOCATE(plev_full_o3(nplev_o3))
      ALLOCATE(plev_half_o3(nplev_o3+1))

      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF

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

      pre_year=year

    END IF

  END SUBROUTINE read_amip_o3
END MODULE mo_o3
