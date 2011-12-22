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
  USE mo_ha_stepping,         ONLY: prepare_ha_dyn!, initcond_ha_dyn
  USE mo_icoham_dyn_memory,   ONLY: p_hydro_state, destruct_icoham_dyn_state
  
  USE mo_parallel_config,    ONLY: itype_comm, iorder_sendrecv
  USE mo_sync,               ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_patch_array_mult
  USE mo_icon_comm_interface,ONLY: construct_icoham_communication, destruct_icoham_communication
  USE mo_icon_comm_lib

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

    ! 2D variables
    REAL(wp), POINTER :: pnt_2D_cells(:,:), pnt_2D_edges(:,:), pnt_2D_verts(:,:)
    
    INTEGER :: timer_sync_1_1_2D_cells, timer_sync_1_1_2D_edges, timer_sync_1_1_2D_verts, &
      & timer_sync_1_1_2D_all
    INTEGER :: timer_sync_1_2_2D_cells, timer_sync_1_2_2D_edges, timer_sync_1_2_2D_verts, &
      & timer_sync_1_2_2D_all
    INTEGER :: timer_sync_2_1_2D_cells, timer_sync_2_1_2D_edges, timer_sync_2_1_2D_verts, &
      & timer_sync_2_1_2D_all
    INTEGER :: timer_sync_2_2_2D_cells, timer_sync_2_2_2D_edges, timer_sync_2_2_2D_verts, &
      & timer_sync_2_2_2D_all
    INTEGER :: timer_sync_3_1_2D_cells, timer_sync_3_1_2D_edges, timer_sync_3_1_2D_verts, &
      & timer_sync_3_1_2D_all
    INTEGER :: timer_sync_3_2_2D_cells, timer_sync_3_2_2D_edges, timer_sync_3_2_2D_verts, &
      & timer_sync_3_2_2D_all
    
    INTEGER :: timer_iconcom_2D_cells, timer_iconcom_2D_edges, timer_iconcom_2D_verts, &
      & timer_iconcom_2D_all, timer_iconcom_2D_comb, timer_iconcom_2D_keep

    ! 3D variables
    REAL(wp), POINTER :: pnt_3D_cells(:,:,:), pnt_3D_edges(:,:,:), pnt_3D_verts(:,:,:)
    
    INTEGER :: timer_sync_1_1_3D_cells, timer_sync_1_1_3D_edges, timer_sync_1_1_3D_verts, &
      & timer_sync_1_1_3D_all
    INTEGER :: timer_sync_1_2_3D_cells, timer_sync_1_2_3D_edges, timer_sync_1_2_3D_verts, &
      & timer_sync_1_2_3D_all
    INTEGER :: timer_sync_2_1_3D_cells, timer_sync_2_1_3D_edges, timer_sync_2_1_3D_verts, &
      & timer_sync_2_1_3D_all
    INTEGER :: timer_sync_2_2_3D_cells, timer_sync_2_2_3D_edges, timer_sync_2_2_3D_verts, &
      & timer_sync_2_2_3D_all
    INTEGER :: timer_sync_3_1_3D_cells, timer_sync_3_1_3D_edges, timer_sync_3_1_3D_verts, &
      & timer_sync_3_1_3D_all
    INTEGER :: timer_sync_3_2_3D_cells, timer_sync_3_2_3D_edges, timer_sync_3_2_3D_verts, &
      & timer_sync_3_2_3D_all
    
    INTEGER :: timer_iconcom_3D_cells, timer_iconcom_3D_edges, timer_iconcom_3D_verts, &
      & timer_iconcom_3D_all, timer_iconcom_3D_comb, timer_iconcom_3D_keep


    ! other
    INTEGER :: comm_1, comm_2, comm_3
    
    INTEGER :: patch_no, i

    
    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_communication"


    !---------------------------------------------------------------------

    CALL global_mpi_barrier()
    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
    
    !---------------------------------------------------------------------

    ltimer = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    CALL prepare_ha_dyn( p_patch(1:) )
    CALL construct_icoham_communication()
    
    CALL global_mpi_barrier()
    
    !---------------------------------------------------------------------
    patch_no=1
    
    pnt_2D_cells => p_patch(patch_no)%cells%area(:,:)
    pnt_2D_edges => p_patch(patch_no)%edges%primal_edge_length(:,:)
    pnt_2D_verts => p_patch(patch_no)%verts%dual_area(:,:)
    
    pnt_3D_cells => p_hydro_state(patch_no)%prog(1)%temp(:,:,:)
    pnt_3D_edges => p_hydro_state(patch_no)%prog(1)%vn(:,:,:)
    pnt_3D_verts => p_hydro_state(patch_no)%diag%rel_vort(:,:,:)
    
    !---------------------------------------------------------------------
    ! Call cmmunication methods
  ! INTEGER :: itype_comm = 1
  ! 1 = synchronous communication with local memory for exchange buffers
  ! 2 = synchronous communication with global memory for exchange buffers
  ! 3 = asynchronous communication within dynamical core with global memory 
  !     for exchange buffers (not yet implemented)

  ! Order of send/receive sequence in exchange routines
 ! INTEGER :: iorder_sendrecv = 1
  ! 1 = irecv, send
  ! 2 = isend, recv
    
    !---------------------------------------------------------------------
    ! test the 2D sync
    !---------------------------------------------------------------------
    itype_comm = 1
    iorder_sendrecv = 1
    timer_sync_1_1_2D_cells  = new_timer  ("sync_1_1_2D_cells")
    timer_sync_1_1_2D_edges  = new_timer  ("sync_1_1_2D_edges")
    timer_sync_1_1_2D_verts  = new_timer  ("sync_1_1_2D_verts")
    timer_sync_1_1_2D_all    = new_timer  ("sync_1_1_2D_all")

    CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_1_1_2D_cells)
    CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_1_1_2D_edges)
    CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_1_1_2D_verts)
    CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_1_1_2D_all)
        
    !---------------------------------------------------------------------
    itype_comm = 1
    iorder_sendrecv = 2
    timer_sync_1_2_2D_cells  = new_timer  ("sync_1_2_2D_cells")
    timer_sync_1_2_2D_edges  = new_timer  ("sync_1_2_2D_edges")
    timer_sync_1_2_2D_verts  = new_timer  ("sync_1_2_2D_verts")
    timer_sync_1_2_2D_all    = new_timer  ("sync_1_2_2D_all")

    CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_1_2_2D_cells)
    CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_1_2_2D_edges)
    CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_1_2_2D_verts)
    CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_1_2_2D_all)

    !---------------------------------------------------------------------
    itype_comm = 2
    iorder_sendrecv = 1
    timer_sync_2_1_2D_cells  = new_timer  ("sync_2_1_2D_cells")
    timer_sync_2_1_2D_edges  = new_timer  ("sync_2_1_2D_edges")
    timer_sync_2_1_2D_verts  = new_timer  ("sync_2_1_2D_verts")
    timer_sync_2_1_2D_all    = new_timer  ("sync_2_1_2D_all")

    CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_2_1_2D_cells)
    CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_2_1_2D_edges)
    CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_2_1_2D_verts)
    CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_2_1_2D_all)
        
    !---------------------------------------------------------------------
    itype_comm = 2
    iorder_sendrecv = 2
    timer_sync_2_2_2D_cells  = new_timer  ("sync_2_2_2D_cells")
    timer_sync_2_2_2D_edges  = new_timer  ("sync_2_2_2D_edges")
    timer_sync_2_2_2D_verts  = new_timer  ("sync_2_2_2D_verts")
    timer_sync_2_2_2D_all    = new_timer  ("sync_2_2_2D_all")

    CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_2_2_2D_cells)
    CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_2_2_2D_edges)
    CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_2_2_2D_verts)
    CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_2_2_2D_all)
        
    !---------------------------------------------------------------------
    itype_comm = 3
    iorder_sendrecv = 1
    timer_sync_3_1_2D_cells  = new_timer  ("sync_3_1_2D_cells")
    timer_sync_3_1_2D_edges  = new_timer  ("sync_3_1_2D_edges")
    timer_sync_3_1_2D_verts  = new_timer  ("sync_3_1_2D_verts")
    timer_sync_3_1_2D_all    = new_timer  ("sync_3_1_2D_all")

    CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_3_1_2D_cells)
    CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_3_1_2D_edges)
    CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_3_1_2D_verts)
    CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_3_1_2D_all)
        
    !---------------------------------------------------------------------
    itype_comm = 3
    iorder_sendrecv = 2
    timer_sync_3_2_2D_cells  = new_timer  ("sync_3_2_2D_cells")
    timer_sync_3_2_2D_edges  = new_timer  ("sync_3_2_2D_edges")
    timer_sync_3_2_2D_verts  = new_timer  ("sync_3_2_2D_verts")
    timer_sync_3_2_2D_all    = new_timer  ("sync_3_2_2D_all")

    CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_3_2_2D_cells)
    CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_3_2_2D_edges)
    CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_3_2_2D_verts)
    CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_3_2_2D_all)
    

    !---------------------------------------------------------------------
    ! test the 2D icon_comm_lib
    !---------------------------------------------------------------------
    timer_iconcom_2D_cells  = new_timer  ("iconcom_2D_cells")
    timer_iconcom_2D_edges  = new_timer  ("iconcom_2D_edges")
    timer_iconcom_2D_verts  = new_timer  ("iconcom_2D_verts")
    timer_iconcom_2D_all    = new_timer  ("iconcom_2D_all")
    timer_iconcom_2D_comb   = new_timer  ("iconcom_2D_comb")
    timer_iconcom_2D_keep   = new_timer  ("iconcom_2D_keep")
    
    ! test the 2D iconcom on cells
    
    CALL timer_start(timer_iconcom_2D_cells)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_cells, on_cells, p_patch(patch_no))
    ENDDO
    CALL timer_stop(timer_iconcom_2D_cells)

    ! test the 2D iconcom on edges
    CALL timer_start(timer_iconcom_2D_edges)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_edges, on_edges, p_patch(patch_no))
    ENDDO
    CALL timer_stop(timer_iconcom_2D_edges)
    
    ! test the 2D iconcom on verts
    CALL timer_start(timer_iconcom_2D_verts)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_verts, on_verts, p_patch(patch_no))
    ENDDO
    CALL timer_stop(timer_iconcom_2D_verts)
    
    ! test the 2D iconcom on all
    CALL timer_start(timer_iconcom_2D_all)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_cells, on_cells, p_patch(patch_no))
       CALL icon_comm_sync(pnt_2D_edges, on_edges, p_patch(patch_no))
       CALL icon_comm_sync(pnt_2D_verts, on_verts, p_patch(patch_no))
    ENDDO
    CALL timer_stop(timer_iconcom_2D_all)
    
    ! test the 2D iconcom on combined
    CALL timer_start(timer_iconcom_2D_comb)
    DO i=1,testbed_iterations
    
       comm_1 = new_icon_comm_variable(pnt_2D_cells, on_cells, &
         & p_patch(patch_no), status=is_ready, scope=until_sync)
         
       comm_2 = new_icon_comm_variable(pnt_2D_edges, &
         & on_edges, p_patch(patch_no), status=is_ready, scope=until_sync)
       
       comm_3 = new_icon_comm_variable(pnt_2D_verts, &
         & on_verts, p_patch(patch_no), status=is_ready, scope=until_sync)

       CALL icon_comm_sync_all()

    ENDDO
    CALL timer_stop(timer_iconcom_2D_comb)
    
    ! test the 2D iconcom on keeping the communicators
    comm_1 = new_icon_comm_variable(pnt_2D_cells, on_cells, &
      & p_patch(patch_no))
    comm_2 = new_icon_comm_variable(pnt_2D_edges, &
      & on_edges, p_patch(patch_no))
    comm_3 = new_icon_comm_variable(pnt_2D_verts, &
      & on_verts, p_patch(patch_no))
    
    CALL timer_start(timer_iconcom_2D_keep)
    DO i=1,testbed_iterations
    
       CALL icon_comm_var_is_ready(comm_1)
       CALL icon_comm_var_is_ready(comm_2)
       CALL icon_comm_var_is_ready(comm_3)
       CALL icon_comm_sync_all()
       
    ENDDO
    CALL timer_stop(timer_iconcom_2D_keep)
    CALL delete_icon_comm_variable(comm_1)
    CALL delete_icon_comm_variable(comm_2)
    CALL delete_icon_comm_variable(comm_3)
    !---------------------------------------------------------------------
    

    !---------------------------------------------------------------------
    ! test the 3D sync
    !---------------------------------------------------------------------
    itype_comm = 1
    iorder_sendrecv = 1
    timer_sync_1_1_3D_cells  = new_timer  ("sync_1_1_3D_cells")
    timer_sync_1_1_3D_edges  = new_timer  ("sync_1_1_3D_edges")
    timer_sync_1_1_3D_verts  = new_timer  ("sync_1_1_3D_verts")
    timer_sync_1_1_3D_all    = new_timer  ("sync_1_1_3D_all")

    CALL test_sync_3D( SYNC_C, pnt_3D_cells, timer_sync_1_1_3D_cells)
    CALL test_sync_3D( SYNC_E, pnt_3D_edges, timer_sync_1_1_3D_edges)
    CALL test_sync_3D( SYNC_V, pnt_3D_verts, timer_sync_1_1_3D_verts)
    CALL test_sync_3D_all(pnt_3D_cells, pnt_3D_edges, pnt_3D_verts, timer_sync_1_1_3D_all)
        
    !---------------------------------------------------------------------
    itype_comm = 1
    iorder_sendrecv = 2
    timer_sync_1_2_3D_cells  = new_timer  ("sync_1_2_3D_cells")
    timer_sync_1_2_3D_edges  = new_timer  ("sync_1_2_3D_edges")
    timer_sync_1_2_3D_verts  = new_timer  ("sync_1_2_3D_verts")
    timer_sync_1_2_3D_all    = new_timer  ("sync_1_2_3D_all")

    CALL test_sync_3D( SYNC_C, pnt_3D_cells, timer_sync_1_2_3D_cells)
    CALL test_sync_3D( SYNC_E, pnt_3D_edges, timer_sync_1_2_3D_edges)
    CALL test_sync_3D( SYNC_V, pnt_3D_verts, timer_sync_1_2_3D_verts)
    CALL test_sync_3D_all(pnt_3D_cells, pnt_3D_edges, pnt_3D_verts, timer_sync_1_2_3D_all)

    !---------------------------------------------------------------------
    itype_comm = 2
    iorder_sendrecv = 1
    timer_sync_2_1_3D_cells  = new_timer  ("sync_2_1_3D_cells")
    timer_sync_2_1_3D_edges  = new_timer  ("sync_2_1_3D_edges")
    timer_sync_2_1_3D_verts  = new_timer  ("sync_2_1_3D_verts")
    timer_sync_2_1_3D_all    = new_timer  ("sync_2_1_3D_all")

    CALL test_sync_3D( SYNC_C, pnt_3D_cells, timer_sync_2_1_3D_cells)
    CALL test_sync_3D( SYNC_E, pnt_3D_edges, timer_sync_2_1_3D_edges)
    CALL test_sync_3D( SYNC_V, pnt_3D_verts, timer_sync_2_1_3D_verts)
    CALL test_sync_3D_all(pnt_3D_cells, pnt_3D_edges, pnt_3D_verts, timer_sync_2_1_3D_all)
        
    !---------------------------------------------------------------------
    itype_comm = 2
    iorder_sendrecv = 2
    timer_sync_2_2_3D_cells  = new_timer  ("sync_2_2_3D_cells")
    timer_sync_2_2_3D_edges  = new_timer  ("sync_2_2_3D_edges")
    timer_sync_2_2_3D_verts  = new_timer  ("sync_2_2_3D_verts")
    timer_sync_2_2_3D_all    = new_timer  ("sync_2_2_3D_all")

    CALL test_sync_3D( SYNC_C, pnt_3D_cells, timer_sync_2_2_3D_cells)
    CALL test_sync_3D( SYNC_E, pnt_3D_edges, timer_sync_2_2_3D_edges)
    CALL test_sync_3D( SYNC_V, pnt_3D_verts, timer_sync_2_2_3D_verts)
    CALL test_sync_3D_all(pnt_3D_cells, pnt_3D_edges, pnt_3D_verts, timer_sync_2_2_3D_all)
        
    !---------------------------------------------------------------------
    itype_comm = 3
    iorder_sendrecv = 1
    timer_sync_3_1_3D_cells  = new_timer  ("sync_3_1_3D_cells")
    timer_sync_3_1_3D_edges  = new_timer  ("sync_3_1_3D_edges")
    timer_sync_3_1_3D_verts  = new_timer  ("sync_3_1_3D_verts")
    timer_sync_3_1_3D_all    = new_timer  ("sync_3_1_3D_all")

    CALL test_sync_3D( SYNC_C, pnt_3D_cells, timer_sync_3_1_3D_cells)
    CALL test_sync_3D( SYNC_E, pnt_3D_edges, timer_sync_3_1_3D_edges)
    CALL test_sync_3D( SYNC_V, pnt_3D_verts, timer_sync_3_1_3D_verts)
    CALL test_sync_3D_all(pnt_3D_cells, pnt_3D_edges, pnt_3D_verts, timer_sync_3_1_3D_all)
        
    !---------------------------------------------------------------------
    itype_comm = 3
    iorder_sendrecv = 2
    timer_sync_3_2_3D_cells  = new_timer  ("sync_3_2_3D_cells")
    timer_sync_3_2_3D_edges  = new_timer  ("sync_3_2_3D_edges")
    timer_sync_3_2_3D_verts  = new_timer  ("sync_3_2_3D_verts")
    timer_sync_3_2_3D_all    = new_timer  ("sync_3_2_3D_all")

    CALL test_sync_3D( SYNC_C, pnt_3D_cells, timer_sync_3_2_3D_cells)
    CALL test_sync_3D( SYNC_E, pnt_3D_edges, timer_sync_3_2_3D_edges)
    CALL test_sync_3D( SYNC_V, pnt_3D_verts, timer_sync_3_2_3D_verts)
    CALL test_sync_3D_all(pnt_3D_cells, pnt_3D_edges, pnt_3D_verts, timer_sync_3_2_3D_all)
    
    

    !---------------------------------------------------------------------
    ! test the 3D icon_comm_lib
    !---------------------------------------------------------------------
    timer_iconcom_3D_cells  = new_timer  ("iconcom_3D_cells")
    timer_iconcom_3D_edges  = new_timer  ("iconcom_3D_edges")
    timer_iconcom_3D_verts  = new_timer  ("iconcom_3D_verts")
    timer_iconcom_3D_all    = new_timer  ("iconcom_3D_all")
    timer_iconcom_3D_comb   = new_timer  ("iconcom_3D_comb")
    timer_iconcom_3D_keep   = new_timer  ("iconcom_3D_keep")
    
    ! test the 3D iconcom on cells
    CALL timer_start(timer_iconcom_3D_cells)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_3D_cells, on_cells, p_patch(patch_no))
    ENDDO
    CALL timer_stop(timer_iconcom_3D_cells)

    ! test the 3D iconcom on edges
    CALL timer_start(timer_iconcom_3D_edges)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_3D_edges, on_edges, p_patch(patch_no))
    ENDDO
    CALL timer_stop(timer_iconcom_3D_edges)
    
    ! test the 3D iconcom on verts
    CALL timer_start(timer_iconcom_3D_verts)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_3D_verts, on_verts, p_patch(patch_no))
    ENDDO
    CALL timer_stop(timer_iconcom_3D_verts)
    
    ! test the 3D iconcom on all
    CALL timer_start(timer_iconcom_3D_all)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_3D_cells, on_cells, p_patch(patch_no))
       CALL icon_comm_sync(pnt_3D_edges, on_edges, p_patch(patch_no))
       CALL icon_comm_sync(pnt_3D_verts, on_verts, p_patch(patch_no))
    ENDDO
    CALL timer_stop(timer_iconcom_3D_all)
    
    ! test the 3D iconcom on combined
    CALL timer_start(timer_iconcom_3D_comb)
    DO i=1,testbed_iterations
    
       comm_1 = new_icon_comm_variable(pnt_3D_cells, on_cells, &
         & p_patch(patch_no), status=is_ready, scope=until_sync)
         
       comm_2 = new_icon_comm_variable(pnt_3D_edges, &
         & on_edges, p_patch(patch_no), status=is_ready, scope=until_sync)
       
       comm_3 = new_icon_comm_variable(pnt_3D_verts, &
         & on_verts, p_patch(patch_no), status=is_ready, scope=until_sync)

       CALL icon_comm_sync_all()

    ENDDO
    CALL timer_stop(timer_iconcom_3D_comb)
    
    ! test the 3D iconcom on keeping the communicators
    comm_1 = new_icon_comm_variable(pnt_3D_cells, on_cells, &
      & p_patch(patch_no))
    comm_2 = new_icon_comm_variable(pnt_3D_edges, &
      & on_edges, p_patch(patch_no))
    comm_3 = new_icon_comm_variable(pnt_3D_verts, &
      & on_verts, p_patch(patch_no))
    
    CALL timer_start(timer_iconcom_3D_keep)
    DO i=1,testbed_iterations
    
       CALL icon_comm_var_is_ready(comm_1)
       CALL icon_comm_var_is_ready(comm_2)
       CALL icon_comm_var_is_ready(comm_3)
       CALL icon_comm_sync_all()
       
    ENDDO
    CALL timer_stop(timer_iconcom_3D_keep)
    CALL delete_icon_comm_variable(comm_1)
    CALL delete_icon_comm_variable(comm_2)
    CALL delete_icon_comm_variable(comm_3)
    !---------------------------------------------------------------------
    
    

    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_icoham_dyn_state()
    CALL destruct_atmo_model()
    CALL destruct_icoham_communication()
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

  !-------------------------------------------------------------------------
  SUBROUTINE test_sync_2D(comm_pattern, var, timer)
    INTEGER :: comm_pattern, timer
    REAL(wp) , POINTER:: var(:,:)

    INTEGER :: i, patch_no
    
    patch_no = 1
    CALL timer_start(timer)
    DO i=1,testbed_iterations
      CALL sync_patch_array( comm_pattern, p_patch(patch_no), var )
    ENDDO    
    CALL timer_stop(timer)
  END SUBROUTINE test_sync_2D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_sync_2D_all(cell_var, edge_var, vert_var, timer)
    INTEGER :: comm_pattern, timer
    REAL(wp) , POINTER:: cell_var(:,:), edge_var(:,:), vert_var(:,:)

    INTEGER :: i, patch_no
    
    patch_no = 1
    CALL timer_start(timer)
    DO i=1,testbed_iterations
      CALL sync_patch_array(SYNC_C , p_patch(patch_no), cell_var )
      CALL sync_patch_array(SYNC_E , p_patch(patch_no), edge_var )
      CALL sync_patch_array(SYNC_V , p_patch(patch_no), vert_var )
    ENDDO    
    CALL timer_stop(timer)
  END SUBROUTINE test_sync_2D_all
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE test_sync_3D(comm_pattern, var, timer)
    INTEGER :: comm_pattern, timer
    REAL(wp) , POINTER:: var(:,:,:)

    INTEGER :: i, patch_no
    
    patch_no = 1
    CALL timer_start(timer)
    DO i=1,testbed_iterations
      CALL sync_patch_array( comm_pattern, p_patch(patch_no), var )
    ENDDO    
    CALL timer_stop(timer)
  END SUBROUTINE test_sync_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_sync_3D_all(cell_var, edge_var, vert_var, timer)
    INTEGER :: comm_pattern, timer
    REAL(wp) , POINTER:: cell_var(:,:,:), edge_var(:,:,:), vert_var(:,:,:)

    INTEGER :: i, patch_no
    
    patch_no = 1
    CALL timer_start(timer)
    DO i=1,testbed_iterations
      CALL sync_patch_array(SYNC_C , p_patch(patch_no), cell_var )
      CALL sync_patch_array(SYNC_E , p_patch(patch_no), edge_var )
      CALL sync_patch_array(SYNC_V , p_patch(patch_no), vert_var )
    ENDDO    
    CALL timer_stop(timer)
  END SUBROUTINE test_sync_3D_all
  !-------------------------------------------------------------------------

END MODULE mo_test_communication

