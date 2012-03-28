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
MODULE mo_test_nh_communication

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_timer,               ONLY: init_timer, ltimer, new_timer, timer_start, timer_stop, &
    & print_timer

  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no
  USE mo_icon_testbed_config, ONLY: testbed_iterations, calculate_iterations

  USE mo_model_domain,        ONLY:  p_patch  
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model
  
  USE mo_parallel_config,    ONLY: itype_comm, iorder_sendrecv
  USE mo_sync,               ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_patch_array_mult
  USE mo_icon_comm_interface,ONLY: construct_icoham_communication, destruct_icoham_communication
  USE mo_icon_comm_lib
  USE mo_atmo_nonhydrostatic, ONLY: construct_atmo_nonhydrostatic, destruct_atmo_nonhydrostatic

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: test_nh_communication

CONTAINS
  

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_nh_communication(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    ! 2D variables    

    
    INTEGER :: patch_no, i

    
    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_nh_communication"

    !---------------------------------------------------------------------
    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
    
    !---------------------------------------------------------------------

    ltimer = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    CALL construct_atmo_nonhydrostatic()
        
    CALL work_mpi_barrier()
    !---------------------------------------------------------------------
         
    !---------------------------------------------------------------------
    

    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_nonhydrostatic()
    CALL destruct_atmo_model()
    CALL destruct_icoham_communication()
    CALL message(TRIM(method_name),'clean-up finished')

    !---------------------------------------------------------------------
    ! print the timers
!    IF (my_process_is_stdio()) THEN
      CALL message("===================", "=======================")
      WRITE(message_text,*) "Communication Iterations=", testbed_iterations
      CALL message(method_name, TRIM(message_text))
      CALL print_timer()
!    ENDIF
    !---------------------------------------------------------------------
     

  END SUBROUTINE test_nh_communication
  !-------------------------------------------------------------------------

END MODULE mo_test_nh_communication

