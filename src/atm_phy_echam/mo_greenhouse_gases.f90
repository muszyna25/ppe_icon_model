!>
!! Preliminary read and time interpolation of greenhouse gases data
!!
!! This is  a clone of the respective ECHAM routine
!!
!! Time series of various greenhouse gases are read from
!! file greenhouse_gases.nc (CO2, CH4, N2O, and CFC's).
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
MODULE mo_greenhouse_gases

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_physical_constants,   ONLY: amd, amco2, amch4, amn2o
  USE mo_datetime,             ONLY: t_datetime, idaylen
  USE mo_read_netcdf_parallel, ONLY: p_nf_open, p_nf_inq_dimid, p_nf_inq_dimlen, &
       &                             p_nf_inq_varid, p_nf_get_var_double, p_nf_close, &
       &                             nf_read, nf_noerr, nf_strerror

  ! scaling of co2 concentration not implemented yet
  ! USE mo_radiation_parameters, ONLY: fco2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_ghg_bc
  PUBLIC :: ghg_time_interpolation, co2_time_interpolation
  PUBLIC :: cleanup_greenhouse_gases

  PUBLIC :: ghg_no_cfc
  PUBLIC :: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcvmr

  PUBLIC :: ghg_file_read

  INTEGER, PARAMETER :: ghg_no_cfc = 2
  CHARACTER(len=*), PARAMETER :: ghg_cfc_names(ghg_no_cfc) = (/ "CFC_11", "CFC_12" /)

  REAL(wp) :: ghg_base_year
  REAL(wp) :: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcvmr(ghg_no_cfc)

  INTEGER :: ghg_no_years

  REAL(wp), ALLOCATABLE :: ghg_years(:)
  REAL(wp), ALLOCATABLE :: ghg_co2(:)
  REAL(wp), ALLOCATABLE :: ghg_ch4(:)
  REAL(wp), ALLOCATABLE :: ghg_n2o(:)
  REAL(wp), ALLOCATABLE :: ghg_cfc(:,:)

  LOGICAL, SAVE :: ghg_file_read = .FALSE.

CONTAINS

  SUBROUTINE read_ghg_bc(ighg)

    INTEGER, INTENT(in) :: ighg

    INTEGER :: ncid, ndimid, nvarid
    INTEGER :: i

    IF (ghg_file_read) THEN
      CALL message('','Greenhouse gases already read ...')
      RETURN
    ENDIF

    SELECT CASE(ighg) ! select scenario
    CASE (0) 
      CALL message('','Use predefined greenhouse gases secenario from 1990 based on CMIP5')
      RETURN  
    CASE (1)
      CALL nf_check(p_nf_open('bc_greenhouse_gases.nc', nf_read, ncid))
      CALL nf_check(p_nf_inq_dimid (ncid, 'time', ndimid)) 
      CALL nf_check(p_nf_inq_dimlen (ncid, ndimid, ghg_no_years))
    CASE DEFAULT
      CALL finish('','Greenhouse gases scenario not available ...')
    END SELECT

    ALLOCATE (ghg_years(ghg_no_years))
    ALLOCATE (ghg_co2(ghg_no_years))
    ALLOCATE (ghg_ch4(ghg_no_years))
    ALLOCATE (ghg_n2o(ghg_no_years))
    ALLOCATE (ghg_cfc(ghg_no_years,ghg_no_cfc))
    
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
      
    ghg_file_read = .TRUE.

    CALL nf_check(p_nf_close(ncid))

    ghg_base_year = ghg_years(1)
    
  END SUBROUTINE read_ghg_bc
  
  SUBROUTINE ghg_time_interpolation(radiation_date)

    TYPE(t_datetime), INTENT(in) :: radiation_date 

    REAL(wp) :: zsecref, zsecnow
    REAL(wp) :: zw1, zw2
    REAL(wp) :: zco2int, zch4int, zn2oint
    REAL(wp) :: zcfc(ghg_no_cfc)

    INTEGER :: iyear, iyearm, iyearp

    CHARACTER(len=32)  :: cdate, cformat
    CHARACTER(len=256) :: ccfc

    ! interpolation in time

    zsecref = REAL(radiation_date%yealen * idaylen, wp)
    zsecnow = REAL(((radiation_date%yeaday - 1) * idaylen), wp) + radiation_date%daysec

    iyear =  radiation_date%year - INT(ghg_base_year) + 1   ! set right index to access in ghg fields
    iyearm = iyear - 1
    iyearp = iyear + 1

    IF (radiation_date%month <= 6) THEN     ! first half of year

      zw1 = zsecnow/zsecref + 0.5_wp
      zw2 = 1.0_wp - zw1

      zco2int = zw1*ghg_co2(iyear)+zw2*ghg_co2(iyearm)
      zch4int = zw1*ghg_ch4(iyear)+zw2*ghg_ch4(iyearm)
      zn2oint = zw1*ghg_n2o(iyear)+zw2*ghg_n2o(iyearm)
      zcfc(:) = zw1*ghg_cfc(iyear,:)+zw2*ghg_cfc(iyearm,:)

    ELSE                                    ! second half of year

      zw2= zsecnow/zsecref - 0.5_wp
      zw1= 1.0_wp - zw2

      zco2int = zw1*ghg_co2(iyear)+zw2*ghg_co2(iyearp)
      zch4int = zw1*ghg_ch4(iyear)+zw2*ghg_ch4(iyearp)
      zn2oint = zw1*ghg_n2o(iyear)+zw2*ghg_n2o(iyearp)
      zcfc(:) = zw1*ghg_cfc(iyear,:)+zw2*ghg_cfc(iyearp,:)

    END IF

    ! IF (ABS(fco2-1.0_wp) > EPSILON(1.0_wp)) zco2int = fco2 * zco2int

    WRITE(cdate,'(i6,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,f9.6,a)') &
         &      radiation_date%year, '-', radiation_date%month, '-', radiation_date%day, 'T', &
         &      radiation_date%hour, ':', radiation_date%minute, ':', radiation_date%second,'Z'
    WRITE(cformat,'(a,i0,a)') '(a,', ghg_no_cfc, 'f7.2)'
    WRITE(ccfc,cformat) ' CFC = ', zcfc(1:ghg_no_cfc)
    WRITE (message_text,'(a,a,a,f9.4,a,f8.2,a,f8.3,a)') &
         'Greenhouse gas concentrations ', TRIM(cdate), ' CO2 = ', zco2int, &
         ' CH4 = ', zch4int,' N2O = ', zn2oint, TRIM(ccfc)
    CALL message('', TRIM(message_text))

    ! convert to mass mixing ratio

    ghg_co2mmr = zco2int * 1.0e-06_wp * amco2/amd 
    ghg_ch4mmr = zch4int * 1.0e-09_wp * amch4/amd
    ghg_n2ommr = zn2oint * 1.0e-09_wp * amn2o/amd

    ! scale CFCs only, keep the volume mixing ratio 

    ghg_cfcvmr(:) = zcfc(:) * 1.0e-12_wp

  END SUBROUTINE ghg_time_interpolation

  FUNCTION co2_time_interpolation(current_date) RESULT(co2mmr)
    REAL(wp) :: co2mmr
    TYPE(t_datetime), INTENT(in) :: current_date 

    REAL(wp) :: zsecref, zsecnow
    REAL(wp) :: zw1, zw2
    REAL(wp) :: zco2int

    INTEGER :: iyear, iyearm, iyearp
    
    ! interpolation in time

    zsecref = REAL(current_date%yealen * idaylen, wp)
    zsecnow = REAL(((current_date%yeaday - 1) * idaylen), wp) + current_date%second

    iyear =  current_date%year - INT(ghg_base_year) + 1   ! set right index to access in ghg fields
    iyearm = iyear - 1
    iyearp = iyear + 1

    IF (current_date%month <= 6) THEN      ! first half of year

      zw1 = zsecnow/zsecref + 0.5_wp
      zw2 = 1.0_wp - zw1

      zco2int = zw1*ghg_co2(iyear)+zw2*ghg_co2(iyearm)

    ELSE                                   ! second half of year

      zw2= zsecnow/zsecref - 0.5_wp
      zw1= 1.0_wp - zw2

      zco2int = zw1*ghg_co2(iyear)+zw2*ghg_co2(iyearp)

    END IF

    ! convert to mass mixing ratio

    co2mmr = zco2int * 1.0e-06_wp * amco2/amd 

  END FUNCTION co2_time_interpolation

  SUBROUTINE cleanup_greenhouse_gases
    IF (ALLOCATED(ghg_years)) DEALLOCATE(ghg_years)
    IF (ALLOCATED(ghg_co2))   DEALLOCATE(ghg_co2)
    IF (ALLOCATED(ghg_ch4))   DEALLOCATE(ghg_ch4)
    IF (ALLOCATED(ghg_n2o))   DEALLOCATE(ghg_n2o)
    IF (ALLOCATED(ghg_cfc))   DEALLOCATE(ghg_cfc)
  END SUBROUTINE cleanup_greenhouse_gases

  SUBROUTINE nf_check(iret)
    INTEGER, INTENT(in) :: iret
    IF (iret /= nf_noerr) THEN
      CALL finish('mo_greenhouse_gases', nf_strerror(iret))
    ENDIF
  END SUBROUTINE nf_check

END MODULE mo_greenhouse_gases
