!>
!! @brief Main program for the ICON atmospheric model
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!  Hui Wan             (MPI-M)
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
MODULE mo_icon_testbed

  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: global_mpi_barrier
  USE mo_master_control,      ONLY: get_my_process_name

  USE mo_icon_testbed_config, ONLY: testbed_model, null_model, test_coupler_model, &
    & test_jitter_model, test_halo_communication,test_radiation_communication,     &
    & test_netcdf_read_model
  USE mo_icon_testbed_nml,    ONLY: read_icon_testbed_namelist

  USE mo_test_coupler,        ONLY: test_coupler
  USE mo_test_communication,  ONLY: test_communication
  USE mo_test_jitter,         ONLY: test_jitter
  USE mo_test_netcdf_read,    ONLY: test_netcdf_read

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

    CASE(test_coupler_model)
      CALL test_coupler(testbed_namelist_filename,shr_namelist_filename)

    CASE(test_halo_communication,test_radiation_communication)
      CALL test_communication(testbed_namelist_filename,shr_namelist_filename)

    CASE(test_jitter_model)
      CALL test_jitter(testbed_namelist_filename,shr_namelist_filename)

    CASE(test_netcdf_read_model)
      CALL test_netcdf_read(testbed_namelist_filename,shr_namelist_filename)

    CASE default
      CALL finish(method_name, "Unrecognized testbed_model")

    END SELECT    
   

  END SUBROUTINE icon_testbed
  !-------------------------------------------------------------------------

END MODULE mo_icon_testbed

