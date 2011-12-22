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
MODULE mo_test_communication

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: global_mpi_barrier, p_pe_work
  USE mo_timer,               ONLY: init_timer
  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no

  USE mo_model_domain,        ONLY:  p_patch
  
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: test_communication

CONTAINS
!>
!!
  SUBROUTINE test_communication(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_communication"


    !---------------------------------------------------------------------

    CALL global_mpi_barrier()
    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
    
    !---------------------------------------------------------------------

    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    
    CALL global_mpi_barrier()
    

    !---------------------------------------------------------------------
    ! Call cmmunication methods
    !---------------------------------------------------------------------


    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_model()
  

    CALL message(TRIM(method_name),'clean-up finished')

  END SUBROUTINE test_communication
  !-------------------------------------------------------------------------



END MODULE mo_test_communication

