!>
!! Preliminary read and time interpolation of greenhouse gases data
!!
!! This is  a clone of the respective ECHAM routine
!!
!! Time series of various greenhouse gases are read from
!! file bc_greenhouse_gases.nc (CO2, CH4, N2O, and CFC's).
!! Provides interpolation in time and conversion from volume mixing ratio 
!! to mass mixing ratio of CO2, CH4, and N2O - not for CFC's!
!!
!! U. Schlese, DKRZ, June 1995, original source
!! L. Kornblueh, MPI, November 2001, changed to read netCDF input,
!!             packed in a module, f90 rewrite, and parallelization        
!! M. Esch, MPI, May 2004, modified for scenarios
!! M. Esch, MPI, December 2009, modified for CMIP5
!! J. S. Rast, MPI, August 2010, modified interpolation to time step of radiation
!! R. Schnur,  MPI, November 2010, for current time step and CO2 only
!! L. Kornblueh, MPI, March 2013, adapted as temporary reader in ICON
!! 
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_bc_greenhouse_gases

  USE mo_kind,               ONLY: wp, dp, i8
  USE mo_exception,          ONLY: finish, message, message_text, warning
  USE mo_physical_constants, ONLY: amd, amco2, amch4, amn2o, amc11, amc12
  USE mo_netcdf_parallel,    ONLY: p_nf_open, p_nf_inq_dimid, p_nf_inq_dimlen, &
       &                           p_nf_inq_varid, p_nf_get_var_double, p_nf_close, &
       &                           nf_read, nf_noerr, nf_strerror
  USE mtime,                 ONLY: datetime, no_of_sec_in_a_day, &
       &                           getNoOfDaysInYearDateTime, &
       &                           getdayofyearfromdatetime,  &
       &                           getnoofsecondselapsedindaydatetime

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_bc_greenhouse_gases
  PUBLIC :: bc_greenhouse_gases_time_interpolation
  PUBLIC :: cleanup_greenhouse_gases

  PUBLIC :: ghg_no_cfc

  PUBLIC :: bc_greenhouse_gases_file_read
  PUBLIC :: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcmmr

  INTEGER, PARAMETER :: ghg_no_cfc = 2
  CHARACTER(len=*), PARAMETER :: ghg_cfc_names(ghg_no_cfc) = (/ "CFC_11", "CFC_12" /)

  REAL(wp) :: ghg_base_year

  INTEGER :: ghg_no_years

  REAL(wp), ALLOCATABLE :: ghg_years(:)
  REAL(wp), ALLOCATABLE :: ghg_co2(:)
  REAL(wp), ALLOCATABLE :: ghg_ch4(:)
  REAL(wp), ALLOCATABLE :: ghg_n2o(:)
  REAL(wp), ALLOCATABLE :: ghg_cfc(:,:)

  REAL(wp) :: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr
  REAL(wp) :: ghg_cfcmmr(ghg_no_cfc)

  LOGICAL, SAVE :: bc_greenhouse_gases_file_read = .FALSE.

CONTAINS

  SUBROUTINE read_bc_greenhouse_gases

    INTEGER :: ncid, ndimid, nvarid
    INTEGER :: i

    IF (bc_greenhouse_gases_file_read) THEN
      CALL message('','Greenhouse gases already read ...')
      RETURN
    ENDIF

    CALL message('','Use transient, annually resolved greenhouse gases secenario based on CMIP5')
    CALL nf_check(p_nf_open('bc_greenhouse_gases.nc', nf_read, ncid))
    CALL nf_check(p_nf_inq_dimid (ncid, 'time', ndimid)) 
    CALL nf_check(p_nf_inq_dimlen (ncid, ndimid, ghg_no_years))

    ALLOCATE (ghg_years(ghg_no_years))
    ALLOCATE (ghg_co2(ghg_no_years))
    ALLOCATE (ghg_ch4(ghg_no_years))
    ALLOCATE (ghg_n2o(ghg_no_years))
    ALLOCATE (ghg_cfc(ghg_no_years,ghg_no_cfc))
    !$ACC ENTER DATA PCREATE( ghg_years, ghg_co2, ghg_ch4, ghg_n2o, ghg_cfc, ghg_cfcmmr )
    
    CALL nf_check(p_nf_inq_varid(ncid, 'time', nvarid))
    CALL nf_check(p_nf_get_var_double (ncid, nvarid, ghg_years))
      
    CALL nf_check(p_nf_inq_varid (ncid, 'CO2', nvarid))
    CALL nf_check(p_nf_get_var_double (ncid, nvarid, ghg_co2))
      
    CALL nf_check(p_nf_inq_varid (ncid, 'CH4', nvarid))
    CALL nf_check(p_nf_get_var_double (ncid, nvarid, ghg_ch4))
      
    CALL nf_check(p_nf_inq_varid (ncid, 'N2O', nvarid))
    CALL nf_check(p_nf_get_var_double (ncid, nvarid, ghg_n2o))
      
    DO i = 1, ghg_no_cfc
      CALL nf_check(p_nf_inq_varid (ncid, TRIM(ghg_cfc_names(i)), nvarid))
      CALL nf_check(p_nf_get_var_double (ncid, nvarid, ghg_cfc(:,i)))
    ENDDO
      
    bc_greenhouse_gases_file_read = .TRUE.

    CALL nf_check(p_nf_close(ncid))

    ghg_base_year = ghg_years(1)

#ifdef _OPENACC
    CALL warning("GPU:read_bc_greenhouse_gases", "GPU device synchronization")
#endif
    !$ACC UPDATE DEVICE( ghg_years, ghg_co2, ghg_ch4, ghg_n2o, ghg_cfc )
    
  END SUBROUTINE read_bc_greenhouse_gases
  
  SUBROUTINE bc_greenhouse_gases_time_interpolation(radiation_date)

    TYPE(datetime), POINTER, INTENT(in) :: radiation_date 

    REAL(dp) :: zsecref, zsecnow
    REAL(dp) :: zw1, zw2
    REAL(wp) :: zco2int, zch4int, zn2oint
    REAL(wp) :: zcfc(ghg_no_cfc)
    INTEGER(i8) :: yearlen, yearday
    INTEGER :: iyear, iyearm, iyearp

    ! interpolation in time

    yearlen = getNoOfDaysInYearDateTime(radiation_date)*no_of_sec_in_a_day
    yearday = (getdayofyearfromdatetime(radiation_date)-1)*no_of_sec_in_a_day &
         &   +getnoofsecondselapsedindaydatetime(radiation_date)    
    zsecref = REAL(yearlen, dp)
    zsecnow = REAL(yearday, dp)

    iyear  = radiation_date%date%year - INT(ghg_base_year) + 1   ! set right index to access in ghg fields
    iyearm = iyear - 1
    iyearp = iyear + 1

    ! Data are allocated from 1 to ghg_no_years, thus
    ! iyear, iyearm and iyearp shall stay within this
    ! range

    IF (radiation_date%date%month <= 6) THEN     ! first half of year

      IF ( iyear  < 1 .OR. iyear  > ghg_no_years .OR. &
   &       iyearm < 1 .OR. iyearm > ghg_no_years ) THEN

        WRITE (message_text,'(a,i8,a,i8,a,i8)') 'iyear ', iyear,    &
   &                                         ' or iyearm ', iyearm, &
   &                                         ' are out of range 1 - ', ghg_no_years
        CALL finish('mo_bc_greenhouse_gases', message_text)
 
      ENDIF

      zw1 = zsecnow/zsecref + 0.5_dp
      zw2 = 1.0_dp - zw1

      zco2int   = 1.0e-06_wp * ( zw1*ghg_co2(iyear)   + zw2*ghg_co2(iyearm)   )
      zch4int   = 1.0e-09_wp * ( zw1*ghg_ch4(iyear)   + zw2*ghg_ch4(iyearm)   )
      zn2oint   = 1.0e-09_wp * ( zw1*ghg_n2o(iyear)   + zw2*ghg_n2o(iyearm)   )
      zcfc(:)   = 1.0e-12_wp * ( zw1*ghg_cfc(iyear,:) + zw2*ghg_cfc(iyearm,:) )

    ELSE                                         ! second half of year

      IF ( iyear  < 1 .OR. iyear  > ghg_no_years .OR. &
   &       iyearp < 1 .OR. iyearp > ghg_no_years ) THEN

        WRITE (message_text,'(a,i8,a,i8,a,i8)') 'iyear ', iyear,    &
   &                                         ' or iyearp ', iyearp, &
   &                                         ' are out of range 1 - ', ghg_no_years
        CALL finish('mo_bc_greenhouse_gases', message_text)

      ENDIF

      zw2= zsecnow/zsecref - 0.5_dp
      zw1= 1.0_dp - zw2

      zco2int   = 1.0e-06_wp * ( zw1*ghg_co2(iyear)   + zw2*ghg_co2(iyearp)   )
      zch4int   = 1.0e-09_wp * ( zw1*ghg_ch4(iyear)   + zw2*ghg_ch4(iyearp)   )
      zn2oint   = 1.0e-09_wp * ( zw1*ghg_n2o(iyear)   + zw2*ghg_n2o(iyearp)   )
      zcfc(:)   = 1.0e-12_wp * ( zw1*ghg_cfc(iyear,:) + zw2*ghg_cfc(iyearp,:) )
    END IF

    ! convert from volume to mass mixing ratio

    ghg_co2mmr    = zco2int*amco2/amd 
    ghg_ch4mmr    = zch4int*amch4/amd
    ghg_n2ommr    = zn2oint*amn2o/amd

    ghg_cfcmmr(1) = zcfc(1)*amc11/amd
    ghg_cfcmmr(2) = zcfc(2)*amc12/amd
#ifdef _OPENACC
    CALL warning("GPU:bc_greenhouse_gases_time_interpolation", "GPU device synchronization")
#endif
    !$ACC UPDATE DEVICE( ghg_cfcmmr )

  END SUBROUTINE bc_greenhouse_gases_time_interpolation

  SUBROUTINE cleanup_greenhouse_gases
    !$ACC EXIT DATA DELETE( ghg_years ) IF( ALLOCATED(ghg_years) )
    !$ACC EXIT DATA DELETE( ghg_co2 ) IF( ALLOCATED(ghg_co2) )
    !$ACC EXIT DATA DELETE( ghg_ch4 ) IF( ALLOCATED(ghg_ch4) )
    !$ACC EXIT DATA DELETE( ghg_n2o ) IF( ALLOCATED(ghg_n2o) )
    !$ACC EXIT DATA DELETE( ghg_cfc ) IF( ALLOCATED(ghg_cfc) )
    !$ACC EXIT DATA DELETE( ghg_cfcmmr )
    IF (ALLOCATED(ghg_years)) DEALLOCATE(ghg_years)
    IF (ALLOCATED(ghg_co2))   DEALLOCATE(ghg_co2)
    IF (ALLOCATED(ghg_ch4))   DEALLOCATE(ghg_ch4)
    IF (ALLOCATED(ghg_n2o))   DEALLOCATE(ghg_n2o)
    IF (ALLOCATED(ghg_cfc))   DEALLOCATE(ghg_cfc)
  END SUBROUTINE cleanup_greenhouse_gases

  SUBROUTINE nf_check(iret)
    INTEGER, INTENT(in) :: iret
    IF (iret /= nf_noerr) THEN
      CALL finish('mo_bc_greenhouse_gases', nf_strerror(iret))
    ENDIF
  END SUBROUTINE nf_check

END MODULE mo_bc_greenhouse_gases
