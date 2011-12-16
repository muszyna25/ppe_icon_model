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
  USE mo_master_control,      ONLY: is_restart_run, get_my_process_name, &
                                    get_my_model_no

  USE mo_icon_testbed_config, ONLY: testbed_model, null_model, test_coupler
  USE mo_icon_testbed_nml,    ONLY: read_icon_testbed_namelist

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

    CALL read_icon_testbed_namelist(testbed_namelist_filename)
    
    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
   

  END SUBROUTINE icon_testbed
  !-------------------------------------------------------------------------

END MODULE mo_icon_testbed

