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
  USE mo_timer,               ONLY: init_timer, ltimer, new_timer, timer_start, timer_stop, &
    & print_timer

  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no
  USE mo_icon_testbed_config, ONLY: testbed_iterations

  USE mo_model_domain,        ONLY:  p_patch  
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model
  
  USE mo_sync,               ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_patch_array_mult

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

    INTEGER :: timer_sync_2D_cells, timer_sync_2D_edges, timer_sync_2D_verts, timer_sync_2D_all
    INTEGER :: patch_no, i
    
    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_communication"


    !---------------------------------------------------------------------

    CALL global_mpi_barrier()
    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
    
    !---------------------------------------------------------------------

    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    
    CALL global_mpi_barrier()
    
    patch_no=1

    !---------------------------------------------------------------------
    ! Call cmmunication methods
    !---------------------------------------------------------------------
    ! test the 2D sync
    timer_sync_2D_cells  = new_timer  ("sync_2D_cells")
    timer_sync_2D_edges  = new_timer  ("sync_2D_edges")
    timer_sync_2D_verts  = new_timer  ("sync_2D_verts")
    timer_sync_2D_all  = new_timer  ("sync_2D_all")

    ! test the 2D sync on cells
    CALL timer_start(timer_sync_2D_cells)
    DO i=1,testbed_iterations
      CALL sync_patch_array( SYNC_C, p_patch(patch_no), p_patch(patch_no)%cells%area(:,:) )
    ENDDO    
    CALL timer_stop(timer_sync_2D_cells)
    
    ! test the 2D sync on edges
    CALL timer_start(timer_sync_2D_edges)
    DO i=1,testbed_iterations
      CALL sync_patch_array( SYNC_E, p_patch(patch_no), &
        & p_patch(patch_no)%edges%primal_edge_length(:,:) )
    ENDDO    
    CALL timer_stop(timer_sync_2D_edges)
    
    ! test the 2D sync on verts
    CALL timer_start(timer_sync_2D_verts)
    DO i=1,testbed_iterations
      CALL sync_patch_array( SYNC_V, p_patch(patch_no), p_patch(patch_no)%verts%dual_area(:,:) )
    ENDDO    
    CALL timer_stop(timer_sync_2D_verts)

    ! test the 2D sync on cells and edges and verts
    CALL timer_start(timer_sync_2D_all)
    DO i=1,testbed_iterations
      CALL sync_patch_array( SYNC_C, p_patch(patch_no), p_patch(patch_no)%cells%area(:,:) )
      CALL sync_patch_array( SYNC_E, p_patch(patch_no), &
        & p_patch(patch_no)%edges%primal_edge_length(:,:) )
      CALL sync_patch_array( SYNC_V, p_patch(patch_no), p_patch(patch_no)%verts%dual_area(:,:) )
    ENDDO    
    CALL timer_stop(timer_sync_2D_all)
    !---------------------------------------------------------------------
    ! test the 3D sync
    

    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_model()  
    CALL message(TRIM(method_name),'clean-up finished')

    !---------------------------------------------------------------------
    ! print the timers
    CALL message("===================", "=======================")
    WRITE(message_text,*) "Communication Iterations=", testbed_iterations
    CALL message(method_name, TRIM(message_text))
    CALL print_timer()
    !---------------------------------------------------------------------
     

  END SUBROUTINE test_communication
  !-------------------------------------------------------------------------



END MODULE mo_test_communication

