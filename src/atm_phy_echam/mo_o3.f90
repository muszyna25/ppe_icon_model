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
  USE mo_datetime,                 ONLY: t_datetime, aux_datetime, print_datetime_all

  IMPLICIT NONE
  PRIVATE
  REAL(wp), PARAMETER               :: vmr2mmr_o3=amo3/amd               ! Volume mixing ratio to mass mixing ratio
  INTEGER, SAVE                     :: pre_year=-999999                  ! Variable to check if it is time to read

  PUBLIC                            :: o3_plev, plev_half_o3, plev_full_o3, nplev_o3
  PUBLIC                            :: read_amip_o3, time_weights

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

      IF(ALLOCATED(plev_full_o3)) DEALLOCATE(plev_full_o3)
      IF(ALLOCATED(plev_half_o3)) DEALLOCATE(plev_half_o3)
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
  SUBROUTINE time_weights(event_date, wgt1, wgt2, inm1, inm2)
    !
    TYPE(t_datetime), INTENT(in) :: event_date !any event date
    REAL(wp), INTENT(out)        :: wgt1, wgt2 !interpolation weights for months with
    INTEGER, INTENT(out)         :: inm1, inm2 !inm1 and inm2 as indices [0,..,13] respectively
    
    REAL(wp)                     :: zcmonfrc ! current month fraction
    REAL(wp)                     :: zevent_tim !time in month
    REAL(wp)                     :: zcmlen2, znmlen2 !half of current/nearest month length
    TYPE(t_datetime)             :: znevent_date

    zcmonfrc=event_date%monfrc
    zevent_tim=event_date%montim ! event time since first of month in frac. days
    zcmlen2=event_date%monlen*0.5_wp
    IF (zcmonfrc<=0.5_wp) THEN
      inm1=event_date%month-1  !interpolate between value of previous and current month, 
                               !"nearest" is previous month
      inm2=event_date%month
      znevent_date=event_date
      znevent_date%day=1
      IF (inm1 == 0) THEN
        znevent_date%month=12
        znevent_date%year=event_date%year-1
      END IF
      CALL aux_datetime(znevent_date)
      znmlen2=znevent_date%monlen*0.5_wp
      wgt1=(zcmlen2-zevent_tim)/(zcmlen2+znmlen2)
      wgt2=1._wp-wgt1
    ELSE
      inm1=event_date%month
      inm2=event_date%month+1
      znevent_date=event_date
      znevent_date%day=1
      IF (inm2 == 13) THEN
        znevent_date%month=1
        znevent_date%year=znevent_date%year+1
      END IF
      CALL aux_datetime(znevent_date)
      znmlen2=znevent_date%monlen*0.5_wp
      wgt2=(zevent_tim-zcmlen2)/(zcmlen2+znmlen2)
      wgt1=1._wp-wgt2
    END IF
  write(0,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  CALL print_datetime_all(event_date)
  write(0,*) 'inm1,inm2,wgt1,wgt2= ',inm1, inm2, wgt1, wgt2
  write(0,*) '=============================================================='
  END SUBROUTINE time_weights 
END MODULE mo_o3
