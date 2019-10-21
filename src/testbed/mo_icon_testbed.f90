!>
!! @brief Main program for the ICON atmospheric model
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!  Hui Wan             (MPI-M)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_icon_testbed

  USE mo_exception,           ONLY: finish
  USE mo_master_control,      ONLY: get_my_process_name

  USE mo_icon_testbed_config, ONLY: testbed_model, null_model, test_coupler_model, &
    & test_jitter_model, test_halo_communication, test_netcdf_read_model,          &
    & test_gather_communication, test_exchange_communication,                      &
    & test_bench_exchange_data_mult
  USE mo_icon_testbed_nml,    ONLY: read_icon_testbed_namelist

#ifndef __NO_ICON_ATMO__
  USE mo_test_communication,  ONLY: test_communication
  USE mo_test_jitter,         ONLY: test_jitter
  USE mo_test_netcdf_read,    ONLY: test_netcdf_read
#endif
!-------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: icon_testbed

CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE icon_testbed(testbed_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: testbed_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_icon_testbed:icon_testbed"

    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
    
    CALL read_icon_testbed_namelist(testbed_namelist_filename)

    SELECT CASE(testbed_model)
    
    CASE(null_model)
      ! do nothing
      RETURN


#ifndef __NO_ICON_ATMO__

    CASE(test_halo_communication, test_gather_communication, &
         test_exchange_communication, test_bench_exchange_data_mult)
      CALL test_communication(testbed_namelist_filename,shr_namelist_filename)

    CASE(test_jitter_model)
      CALL test_jitter(testbed_namelist_filename,shr_namelist_filename)

    CASE(test_netcdf_read_model)
      CALL test_netcdf_read(testbed_namelist_filename,shr_namelist_filename)
#endif

    CASE default
      CALL finish(method_name, "Unrecognized testbed_model")

    END SELECT    
   

  END SUBROUTINE icon_testbed
  !-------------------------------------------------------------------------

END MODULE mo_icon_testbed

