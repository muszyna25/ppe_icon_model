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
MODULE mo_test_communication

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, p_n_work, p_pe_work, &
    &                               p_comm_work, p_comm_rank, p_comm_size, &
    &                               p_barrier
  USE mo_timer,               ONLY: ltimer, new_timer, timer_start, timer_stop, &
    & print_timer, activate_sync_timers, timers_level, timer_barrier, timer_radiaton_recv
  USE mo_parallel_config,     ONLY: nproma, icon_comm_method
  USE mo_communication,       ONLY: t_comm_gather_pattern, &
    &                               setup_comm_gather_pattern, &
    &                               delete_comm_gather_pattern, exchange_data, &
    &                               t_comm_allgather_pattern, &
    &                               setup_comm_allgather_pattern, &
    &                               delete_comm_allgather_pattern

  USE mo_master_control,      ONLY: get_my_process_name

  USE mo_model_domain,        ONLY: p_patch
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model
  USE mo_icoham_dyn_memory,   ONLY: p_hydro_state
  USE mo_math_gradients,      ONLY: grad_fd_norm

  !nh utils
  USE mo_nonhydro_state,      ONLY: p_nh_state
  USE mo_atmo_nonhydrostatic, ONLY: construct_atmo_nonhydrostatic, destruct_atmo_nonhydrostatic

  USE mo_parallel_config,    ONLY: itype_comm, iorder_sendrecv
  USE mo_sync,               ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_patch_array_mult
  USE mo_icon_comm_lib

  USE mo_icon_testbed_config, ONLY: testbed_model, test_halo_communication, &
    & test_radiation_communication, testbed_iterations, calculate_iterations, &
    & test_gather_communication


!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: test_communication

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_communication(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_communication"


    !---------------------------------------------------------------------

    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name

    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    CALL construct_atmo_nonhydrostatic()
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! DO the tests
!     ltimer = .true.
!     activate_sync_timers = .true.
    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    SELECT CASE(testbed_model)

    CASE(test_halo_communication)
      CALL halo_communication_3D_testbed()

    CASE(test_radiation_communication)
      CALL radiation_communication_testbed()

    CASE(test_gather_communication)
      CALL gather_communication_testbed()

    CASE default
      CALL finish(method_name, "Unrecognized testbed_model")

    END SELECT
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! print the timers
!    IF (my_process_is_stdio()) THEN
    CALL message("===================", "=======================")
    WRITE(message_text,*) "Communication Iterations=", testbed_iterations
    CALL message(method_name, TRIM(message_text))
    CALL print_timer()
!    ENDIF
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_nonhydrostatic()
    CALL destruct_atmo_model()
!     CALL destruct_icon_communication()
    CALL message(TRIM(method_name),'clean-up finished')



  END SUBROUTINE test_communication
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE halo_communication_3D_testbed()


    !---------------------------------------------------------------------
    ! test the 3D sync
    !---------------------------------------------------------------------
    itype_comm = 1
    iorder_sendrecv = 1
    CALL test_oldsync_cells_3D("sync_1_1")
    CALL test_oldsync_edges_3D("sync_1_1")

    itype_comm = 1
    iorder_sendrecv = 3
    CALL test_oldsync_cells_3D("sync_1_3")
    CALL test_oldsync_edges_3D("sync_1_3")

    icon_comm_method = 1
    CALL test_iconcom_cells_3D("comm_1")
    CALL test_iconcom_edges_3D("comm_1")

    icon_comm_method = 2
    CALL test_iconcom_cells_3D("com_2")
    CALL test_iconcom_edges_3D("com_2")

    icon_comm_method = 3
    CALL test_iconcom_cells_3D("com_3")
    CALL test_iconcom_edges_3D("com_3")

    icon_comm_method = 102
    CALL test_iconcom_cells_3D("com_102")
    CALL test_iconcom_edges_3D("com_102")

    icon_comm_method = 103
    CALL test_iconcom_cells_3D("com_103")
    CALL test_iconcom_edges_3D("com_103")

    RETURN
  END SUBROUTINE halo_communication_3D_testbed
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_iconcom_cells_3D(timer_descr)
    CHARACTER(len=*), INTENT(in) :: timer_descr

    ! 3D variables
    REAL(wp), POINTER :: pnt_3D_cells_1(:,:,:), pnt_3D_cells_2(:,:,:), pnt_3D_cells_3(:,:,:), &
      & pnt_3D_cells_4(:,:,:)

    INTEGER :: timer_3D_cells_1, timer_3D_cells_2, timer_3D_cells_3, timer_3D_cells_4

    INTEGER :: patch_no

    patch_no=1

    pnt_3D_cells_1 => p_nh_state(patch_no)%prog(1)%w(:,:,:)
    pnt_3D_cells_2 => p_nh_state(patch_no)%prog(1)%rho(:,:,:)
    pnt_3D_cells_3 => p_nh_state(patch_no)%prog(1)%exner(:,:,:)
    pnt_3D_cells_4 => p_nh_state(patch_no)%prog(1)%theta_v(:,:,:)

    pnt_3D_cells_1(:,:,:) = 0.0_wp
    pnt_3D_cells_2(:,:,:) = 0.0_wp
    pnt_3D_cells_3(:,:,:) = 0.0_wp
    pnt_3D_cells_4(:,:,:) = 0.0_wp

    timer_3D_cells_1  = new_timer  (timer_descr//"_3dcells_1")
    CALL test_iconcom_3D(p_patch(patch_no)%sync_cells_not_in_domain, &
      & var1=pnt_3D_cells_1, &
      & timer_id=timer_3D_cells_1)

    timer_3D_cells_2  = new_timer  (timer_descr//"_3dcells_2")
    CALL test_iconcom_3D(p_patch(patch_no)%sync_cells_not_in_domain, &
      & var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & timer_id=timer_3D_cells_2)

    timer_3D_cells_3  = new_timer  (timer_descr//"_3dcells_3")
    CALL test_iconcom_3D(p_patch(patch_no)%sync_cells_not_in_domain, &
      & var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & var3=pnt_3D_cells_3, &
      & timer_id=timer_3D_cells_3)

    timer_3D_cells_4  = new_timer  (timer_descr//"_3dcells_4")
    CALL test_iconcom_3D(p_patch(patch_no)%sync_cells_not_in_domain, &
      & var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & var3=pnt_3D_cells_3,var4=pnt_3D_cells_4, &
      & timer_id=timer_3D_cells_4)

  END SUBROUTINE test_iconcom_cells_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_iconcom_edges_3D(timer_descr)
    CHARACTER(len=*), INTENT(in) :: timer_descr

    ! 3D variables
    REAL(wp), POINTER :: pnt_3D_edges_1(:,:,:), pnt_3D_edges_2(:,:,:), pnt_3D_edges_3(:,:,:)

    INTEGER :: timer_3D_edges_1, timer_3D_edges_2, timer_3D_edges_3

    INTEGER :: patch_no

    patch_no=1

    pnt_3D_edges_1 => p_nh_state(patch_no)%diag%ddt_vn_phy(:,:,:)
    pnt_3D_edges_2 => p_nh_state(patch_no)%diag%mass_fl_e(:,:,:)
    pnt_3D_edges_3 => p_nh_state(patch_no)%diag%vt(:,:,:)
!     pnt_3D_edges_4 => p_nh_state(patch_no)%diag%hfl_tracer(:,:,:,1)

    pnt_3D_edges_1(:,:,:) = 0.0_wp
    pnt_3D_edges_2(:,:,:) = 0.0_wp
    pnt_3D_edges_3(:,:,:) = 0.0_wp
!     pnt_3D_edges_4(:,:,:) = 0.0_wp

    timer_3D_edges_1  = new_timer  (timer_descr//"_3dedges_1")
    CALL test_iconcom_3D(p_patch(patch_no)%sync_cells_not_in_domain, &
      & var1=pnt_3D_edges_1, timer_id=timer_3D_edges_1)

    timer_3D_edges_2  = new_timer  (timer_descr//"_3dedges_2")
    CALL test_iconcom_3D(p_patch(patch_no)%sync_cells_not_in_domain, &
      & var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
      & timer_id=timer_3D_edges_2)

    timer_3D_edges_3  = new_timer  (timer_descr//"_3dedges_3")
    CALL test_iconcom_3D(p_patch(patch_no)%sync_cells_not_in_domain, &
      & var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
      & var3=pnt_3D_edges_3, &
      & timer_id=timer_3D_edges_3)

!     timer_3D_edges_4  = new_timer  (timer_descr//"_3dedges_4")
!     CALL test_iconcom_3D(p_patch(patch_no)%sync_cells_not_in_domain, &
!       & var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
!       & var3=pnt_3D_edges_3,var4=pnt_3D_edges_4, &
!       & timer_id=timer_3D_edges_4)

  END SUBROUTINE test_iconcom_edges_3D
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_oldsync_cells_3D(timer_descr)
    CHARACTER(len=*), INTENT(in) :: timer_descr

    ! 3D variables
    REAL(wp), POINTER :: pnt_3D_cells_1(:,:,:), pnt_3D_cells_2(:,:,:), pnt_3D_cells_3(:,:,:), &
      & pnt_3D_cells_4(:,:,:)

    INTEGER :: timer_3D_cells_1, timer_3D_cells_2, timer_3D_cells_3, timer_3D_cells_4

    INTEGER :: timer_sync_1_2_3D_cells


    INTEGER :: patch_no

    patch_no=1

    pnt_3D_cells_1 => p_nh_state(patch_no)%prog(1)%w(:,:,:)
    pnt_3D_cells_2 => p_nh_state(patch_no)%prog(1)%rho(:,:,:)
    pnt_3D_cells_3 => p_nh_state(patch_no)%prog(1)%exner(:,:,:)
    pnt_3D_cells_4 => p_nh_state(patch_no)%prog(1)%theta_v(:,:,:)
    pnt_3D_cells_1(:,:,:) = 0.0_wp
    pnt_3D_cells_2(:,:,:) = 0.0_wp
    pnt_3D_cells_3(:,:,:) = 0.0_wp
    pnt_3D_cells_4(:,:,:) = 0.0_wp

    timer_3D_cells_1  = new_timer  (timer_descr//"_3dcells_1")
    CALL test_oldsync_3D(SYNC_C, var1=pnt_3D_cells_1, &
      & timer_id=timer_3D_cells_1, patch_no=patch_no)

    timer_3D_cells_2  = new_timer  (timer_descr//"_3dcells_2")
    CALL test_oldsync_3D(SYNC_C, var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & timer_id=timer_3D_cells_2, patch_no=patch_no)

    timer_3D_cells_3  = new_timer  (timer_descr//"_3dcells_3")
    CALL test_oldsync_3D(SYNC_C, var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & var3=pnt_3D_cells_3, &
      & timer_id=timer_3D_cells_3, patch_no=patch_no)

    timer_3D_cells_4  = new_timer  (timer_descr//"_3dcells_4")
    CALL test_oldsync_3D(SYNC_C, var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & var3=pnt_3D_cells_3,var4=pnt_3D_cells_4, &
      & timer_id=timer_3D_cells_4, patch_no=patch_no)

  END SUBROUTINE test_oldsync_cells_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_oldsync_edges_3D(timer_descr)
    CHARACTER(len=*), INTENT(in) :: timer_descr

    ! 3D variables
    REAL(wp), POINTER :: pnt_3D_edges_1(:,:,:), pnt_3D_edges_2(:,:,:), pnt_3D_edges_3(:,:,:)

    INTEGER :: timer_3D_edges_1, timer_3D_edges_2, timer_3D_edges_3

    INTEGER :: patch_no

    patch_no=1

    pnt_3D_edges_1 => p_nh_state(patch_no)%diag%ddt_vn_phy(:,:,:)
    pnt_3D_edges_2 => p_nh_state(patch_no)%diag%mass_fl_e(:,:,:)
    pnt_3D_edges_3 => p_nh_state(patch_no)%diag%vt(:,:,:)
!     pnt_3D_edges_4 => p_nh_state(patch_no)%diag%hfl_tracer(:,:,:,1)
    pnt_3D_edges_1(:,:,:) = 0.0_wp
    pnt_3D_edges_2(:,:,:) = 0.0_wp
    pnt_3D_edges_3(:,:,:) = 0.0_wp
!     pnt_3D_edges_4(:,:,:) = 0.0_wp

    timer_3D_edges_1  = new_timer  (timer_descr//"_3dedges_1")
    CALL test_oldsync_3D(SYNC_E, var1=pnt_3D_edges_1, &
      & timer_id=timer_3D_edges_1, patch_no=patch_no)

    timer_3D_edges_2  = new_timer  (timer_descr//"_3dedges_2")
    CALL test_oldsync_3D(SYNC_E, var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
      & timer_id=timer_3D_edges_2, patch_no=patch_no)

    timer_3D_edges_3  = new_timer  (timer_descr//"_3dedges_3")
    CALL test_oldsync_3D(SYNC_E, var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
      & var3=pnt_3D_edges_3, &
      & timer_id=timer_3D_edges_3, patch_no=patch_no)

!     timer_3D_edges_4  = new_timer  (timer_descr//"_3dedges_4")
!     CALL test_oldsync_3D(SYNC_E, var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
!       & var3=pnt_3D_edges_3,var4=pnt_3D_edges_4, &
!       & timer_id=timer_3D_edges_4, patch_no=patch_no)

  END SUBROUTINE test_oldsync_edges_3D
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE test_communication_2D

    ! 2D variables
    REAL(wp), POINTER :: pnt_2D_cells(:,:), pnt_2D_edges(:,:), pnt_2D_verts(:,:)

    INTEGER :: timer_sync_1_1_2D_cells, timer_sync_1_1_2D_edges, timer_sync_1_1_2D_verts, &
      & timer_sync_1_1_2D_all
    INTEGER :: timer_sync_1_2_2D_cells, timer_sync_1_2_2D_edges, timer_sync_1_2_2D_verts, &
      & timer_sync_1_2_2D_all
!     INTEGER :: timer_sync_2_1_2D_cells, timer_sync_2_1_2D_edges, timer_sync_2_1_2D_verts, &
!       & timer_sync_2_1_2D_all
!     INTEGER :: timer_sync_2_2_2D_cells, timer_sync_2_2_2D_edges, timer_sync_2_2_2D_verts, &
!       & timer_sync_2_2_2D_all
!     INTEGER :: timer_sync_3_1_2D_cells, timer_sync_3_1_2D_edges, timer_sync_3_1_2D_verts, &
!       & timer_sync_3_1_2D_all
!     INTEGER :: timer_sync_3_2_2D_cells, timer_sync_3_2_2D_edges, timer_sync_3_2_2D_verts, &
!       & timer_sync_3_2_2D_all

    INTEGER :: timer_iconcom_2D_cells, timer_iconcom_2D_edges, timer_iconcom_2D_verts, &
      & timer_iconcom_2D_all, timer_iconcom_2D_comb, timer_iconcom_2D_keep


    ! other
    INTEGER :: comm_1, comm_2, comm_3

    INTEGER :: patch_no, i


    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_communication"

    !---------------------------------------------------------------------
    patch_no=1

    pnt_2D_cells => p_patch(patch_no)%cells%area(:,:)
    pnt_2D_edges => p_patch(patch_no)%edges%primal_edge_length(:,:)
    pnt_2D_verts => p_patch(patch_no)%verts%dual_area(:,:)
    !---------------------------------------------------------------------
    ! test the 2D sync
    !---------------------------------------------------------------------
    itype_comm = 1
    iorder_sendrecv = 1
    timer_sync_1_1_2D_cells  = new_timer  ("sync_1_1_2D_cells")
    timer_sync_1_1_2D_edges  = new_timer  ("sync_1_1_2D_edges")
    timer_sync_1_1_2D_verts  = new_timer  ("sync_1_1_2D_verts")
    timer_sync_1_1_2D_all    = new_timer  ("sync_1_1_2D_all")

    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_1_1_2D_cells)
    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_1_1_2D_edges)
    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_1_1_2D_verts)
    CALL work_mpi_barrier()
    CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_1_1_2D_all)

    !---------------------------------------------------------------------
    itype_comm = 1
    iorder_sendrecv = 2
    timer_sync_1_2_2D_cells  = new_timer  ("sync_1_2_2D_cells")
    timer_sync_1_2_2D_edges  = new_timer  ("sync_1_2_2D_edges")
    timer_sync_1_2_2D_verts  = new_timer  ("sync_1_2_2D_verts")
    timer_sync_1_2_2D_all    = new_timer  ("sync_1_2_2D_all")

    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_1_2_2D_cells)
    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_1_2_2D_edges)
    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_1_2_2D_verts)
    CALL work_mpi_barrier()
    CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_1_2_2D_all)

    !---------------------------------------------------------------------
!     itype_comm = 2
!     iorder_sendrecv = 1
!     timer_sync_2_1_2D_cells  = new_timer  ("sync_2_1_2D_cells")
!     timer_sync_2_1_2D_edges  = new_timer  ("sync_2_1_2D_edges")
!     timer_sync_2_1_2D_verts  = new_timer  ("sync_2_1_2D_verts")
!     timer_sync_2_1_2D_all    = new_timer  ("sync_2_1_2D_all")
!
!     CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_2_1_2D_cells)
!     CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_2_1_2D_edges)
!     CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_2_1_2D_verts)
!     CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_2_1_2D_all)
!
!     !---------------------------------------------------------------------
!     itype_comm = 2
!     iorder_sendrecv = 2
!     timer_sync_2_2_2D_cells  = new_timer  ("sync_2_2_2D_cells")
!     timer_sync_2_2_2D_edges  = new_timer  ("sync_2_2_2D_edges")
!     timer_sync_2_2_2D_verts  = new_timer  ("sync_2_2_2D_verts")
!     timer_sync_2_2_2D_all    = new_timer  ("sync_2_2_2D_all")
!
!     CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_2_2_2D_cells)
!     CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_2_2_2D_edges)
!     CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_2_2_2D_verts)
!     CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_2_2_2D_all)
!
!     !---------------------------------------------------------------------
!     itype_comm = 3
!     iorder_sendrecv = 1
!     timer_sync_3_1_2D_cells  = new_timer  ("sync_3_1_2D_cells")
!     timer_sync_3_1_2D_edges  = new_timer  ("sync_3_1_2D_edges")
!     timer_sync_3_1_2D_verts  = new_timer  ("sync_3_1_2D_verts")
!     timer_sync_3_1_2D_all    = new_timer  ("sync_3_1_2D_all")
!
!     CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_3_1_2D_cells)
!     CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_3_1_2D_edges)
!     CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_3_1_2D_verts)
!     CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_3_1_2D_all)
!
!     !---------------------------------------------------------------------
!     itype_comm = 3
!     iorder_sendrecv = 2
!     timer_sync_3_2_2D_cells  = new_timer  ("sync_3_2_2D_cells")
!     timer_sync_3_2_2D_edges  = new_timer  ("sync_3_2_2D_edges")
!     timer_sync_3_2_2D_verts  = new_timer  ("sync_3_2_2D_verts")
!     timer_sync_3_2_2D_all    = new_timer  ("sync_3_2_2D_all")
!
!     CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_3_2_2D_cells)
!     CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_3_2_2D_edges)
!     CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_3_2_2D_verts)
!     CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_3_2_2D_all)


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

    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_cells)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_cells, p_patch(patch_no)%sync_cells_not_in_domain)
    ENDDO
    CALL timer_stop(timer_iconcom_2D_cells)

    ! test the 2D iconcom on edges
    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_edges)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_edges, p_patch(patch_no)%sync_edges_not_owned)
    ENDDO
    CALL timer_stop(timer_iconcom_2D_edges)

    ! test the 2D iconcom on verts
    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_verts)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_verts, p_patch(patch_no)%sync_verts_not_owned)
    ENDDO
    CALL timer_stop(timer_iconcom_2D_verts)

    ! test the 2D iconcom on all
    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_all)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_cells, p_patch(patch_no)%sync_cells_not_in_domain)
       CALL icon_comm_sync(pnt_2D_edges, p_patch(patch_no)%sync_edges_not_owned)
       CALL icon_comm_sync(pnt_2D_verts, p_patch(patch_no)%sync_verts_not_owned)
    ENDDO
    CALL timer_stop(timer_iconcom_2D_all)

    ! test the 2D iconcom on combined
    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_comb)
    DO i=1,testbed_iterations

       comm_1 = new_icon_comm_variable(pnt_2D_cells, &
         &  p_patch(patch_no)%sync_cells_not_in_domain, &
         &  status=is_ready, scope=until_sync)

       comm_2 = new_icon_comm_variable(pnt_2D_edges, &
         & p_patch(patch_no)%sync_edges_not_owned, status=is_ready, &
         & scope=until_sync)

       comm_3 = new_icon_comm_variable(pnt_2D_verts, &
         & p_patch(patch_no)%sync_verts_not_owned, status=is_ready, &
         & scope=until_sync)

       CALL icon_comm_sync_all()

    ENDDO
    CALL timer_stop(timer_iconcom_2D_comb)

    ! test the 2D iconcom on keeping the communicators
    comm_1 = new_icon_comm_variable(pnt_2D_cells, p_patch(patch_no)%sync_cells_not_in_domain)
    comm_2 = new_icon_comm_variable(pnt_2D_edges, &
      & p_patch(patch_no)%sync_edges_not_owned)
    comm_3 = new_icon_comm_variable(pnt_2D_verts, &
      & p_patch(patch_no)%sync_verts_not_owned)

    CALL work_mpi_barrier()
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
  END SUBROUTINE test_communication_2D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_iconcom_3D(comm_pattern, var1, var2, var3, var4, timer_id)
    INTEGER :: comm_pattern, timer_id, patch_no
    REAL(wp) , POINTER:: var1(:,:,:)
    REAL(wp) , POINTER, OPTIONAL :: var2(:,:,:), var3(:,:,:), var4(:,:,:)

    INTEGER :: i

    CALL work_mpi_barrier()

    IF (.NOT. PRESENT(var2)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL icon_comm_sync(var1, comm_pattern)
        CALL timer_stop(timer_id)
      ENDDO
      RETURN
    ENDIF

    IF (PRESENT(var4)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL icon_comm_sync(var1, var2, var3, var4, comm_pattern)
        CALL timer_stop(timer_id)
      ENDDO

   ELSEIF (PRESENT(var3)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL icon_comm_sync(var1, var2, var3,  comm_pattern)
        CALL timer_stop(timer_id)
      ENDDO
   ELSE
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL icon_comm_sync(var1, var2,  comm_pattern)
        CALL timer_stop(timer_id)
      ENDDO
   ENDIF

   CALL work_mpi_barrier()

  END SUBROUTINE test_iconcom_3D
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE test_oldsync_3D(comm_pattern, var1, var2, var3, var4, timer_id, patch_no)
    INTEGER :: comm_pattern, timer_id, patch_no
    REAL(wp) , POINTER:: var1(:,:,:)
    REAL(wp) , POINTER, OPTIONAL :: var2(:,:,:), var3(:,:,:), var4(:,:,:)

    INTEGER :: i

    CALL work_mpi_barrier()

    IF (.NOT. PRESENT(var2)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL sync_patch_array( comm_pattern, p_patch(patch_no), var1 )
        CALL timer_stop(timer_id)
      ENDDO
      RETURN
    ENDIF

    IF (PRESENT(var4)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL sync_patch_array_mult(comm_pattern, p_patch(patch_no), 4, var1, var2, var3, var4)
        CALL timer_stop(timer_id)
      ENDDO

   ELSEIF (PRESENT(var3)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL sync_patch_array_mult(comm_pattern, p_patch(patch_no), 3, var1, var2, var3)
        CALL timer_stop(timer_id)
      ENDDO
   ELSE
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL sync_patch_array_mult(comm_pattern, p_patch(patch_no), 2, var1, var2)
        CALL timer_stop(timer_id)
      ENDDO
   ENDIF

   CALL work_mpi_barrier()

  END SUBROUTINE test_oldsync_3D
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
      CALL do_calculations()
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
      CALL do_calculations()
    ENDDO
    CALL timer_stop(timer)
  END SUBROUTINE test_sync_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_sync_3D_all(cell_var, edge_var, vert_var, timer)

    REAL(wp) , POINTER:: cell_var(:,:,:), edge_var(:,:,:), vert_var(:,:,:)

    INTEGER :: i, patch_no, timer

    patch_no = 1
    CALL timer_start(timer)
    DO i=1,testbed_iterations
      CALL sync_patch_array(SYNC_C , p_patch(patch_no), cell_var )
      CALL do_calculations()
      CALL sync_patch_array(SYNC_E , p_patch(patch_no), edge_var )
      CALL do_calculations()
      CALL sync_patch_array(SYNC_V , p_patch(patch_no), vert_var )
      CALL do_calculations()
    ENDDO
    CALL timer_stop(timer)
  END SUBROUTINE test_sync_3D_all
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE do_calculations()

    INTEGER :: i, patch_no

    patch_no = 1
    DO i=1,calculate_iterations
      CALL grad_fd_norm (p_hydro_state(patch_no)%prog(1)%temp(:,:,:), &
        & p_patch(patch_no), p_hydro_state(patch_no)%prog(1)%vn(:,:,:))
    ENDDO

  END SUBROUTINE do_calculations
  !-------------------------------------------------------------------------
!     pnt_3D_verts_1 => p_hydro_state(patch_no)%diag%rel_vort(:,:,:)
!
!
!     !---------------------------------------------------------------------
!     ! Call cmmunication methods
!   ! INTEGER :: itype_comm = 1
!   ! 1 = synchronous communication with local memory for exchange buffers
!   ! 2 = synchronous communication with global memory for exchange buffers
!   ! 3 = asynchronous communication within dynamical core with global memory
!   !     for exchange buffers (not yet implemented)
!
!   ! Order of send/receive sequence in exchange routines
!  ! INTEGER :: iorder_sendrecv = 1
!   ! 1 = irecv, send
!   ! 2 = isend, recv
!
!
!
!     timer_sync_1_1_3D_edges  = new_timer  ("sync_1_1_3D_edges")
!     timer_sync_1_1_3D_verts  = new_timer  ("sync_1_1_3D_verts")
!     timer_sync_1_1_3D_all    = new_timer  ("sync_1_1_3D_all")
!
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_1_1_3D_cells)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_1_1_3D_edges)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_1_1_3D_verts)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_1_1_3D_all)
!
!     !---------------------------------------------------------------------
!     itype_comm = 1
!     iorder_sendrecv = 2
!     timer_sync_1_2_3D_cells  = new_timer  ("sync_1_2_3D_cells")
!     timer_sync_1_2_3D_edges  = new_timer  ("sync_1_2_3D_edges")
!     timer_sync_1_2_3D_verts  = new_timer  ("sync_1_2_3D_verts")
!     timer_sync_1_2_3D_all    = new_timer  ("sync_1_2_3D_all")
!
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_1_2_3D_cells)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_1_2_3D_edges)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_1_2_3D_verts)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_1_2_3D_all)
!
!     !---------------------------------------------------------------------
! !     itype_comm = 2
! !     iorder_sendrecv = 1
! !     timer_sync_2_1_3D_cells  = new_timer  ("sync_2_1_3D_cells")
! !     timer_sync_2_1_3D_edges  = new_timer  ("sync_2_1_3D_edges")
! !     timer_sync_2_1_3D_verts  = new_timer  ("sync_2_1_3D_verts")
! !     timer_sync_2_1_3D_all    = new_timer  ("sync_2_1_3D_all")
! !
! !     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_2_1_3D_cells)
! !     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_2_1_3D_edges)
! !     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_2_1_3D_verts)
! !     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_2_1_3D_all)
! !
! !     !---------------------------------------------------------------------
! !     itype_comm = 2
! !     iorder_sendrecv = 2
! !     timer_sync_2_2_3D_cells  = new_timer  ("sync_2_2_3D_cells")
! !     timer_sync_2_2_3D_edges  = new_timer  ("sync_2_2_3D_edges")
! !     timer_sync_2_2_3D_verts  = new_timer  ("sync_2_2_3D_verts")
! !     timer_sync_2_2_3D_all    = new_timer  ("sync_2_2_3D_all")
! !
! !     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_2_2_3D_cells)
! !     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_2_2_3D_edges)
! !     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_2_2_3D_verts)
! !     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_2_2_3D_all)
! !
! !     !---------------------------------------------------------------------
! !     itype_comm = 3
! !     iorder_sendrecv = 1
! !     timer_sync_3_1_3D_cells  = new_timer  ("sync_3_1_3D_cells")
! !     timer_sync_3_1_3D_edges  = new_timer  ("sync_3_1_3D_edges")
! !     timer_sync_3_1_3D_verts  = new_timer  ("sync_3_1_3D_verts")
! !     timer_sync_3_1_3D_all    = new_timer  ("sync_3_1_3D_all")
! !
! !     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_3_1_3D_cells)
! !     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_3_1_3D_edges)
! !     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_3_1_3D_verts)
! !     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_3_1_3D_all)
! !
! !     !---------------------------------------------------------------------
! !     itype_comm = 3
! !     iorder_sendrecv = 2
! !     timer_sync_3_2_3D_cells  = new_timer  ("sync_3_2_3D_cells")
! !     timer_sync_3_2_3D_edges  = new_timer  ("sync_3_2_3D_edges")
! !     timer_sync_3_2_3D_verts  = new_timer  ("sync_3_2_3D_verts")
! !     timer_sync_3_2_3D_all    = new_timer  ("sync_3_2_3D_all")
! !
! !     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_3_2_3D_cells)
! !     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_3_2_3D_edges)
! !     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_3_2_3D_verts)
! !     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_3_2_3D_all)
!
!
!
!     !---------------------------------------------------------------------
!     ! test the 3D icon_comm_lib
!     !---------------------------------------------------------------------
!     timer_iconcom_3D_cells  = new_timer  ("iconcom_3D_cells")
!     timer_iconcom_3D_edges  = new_timer  ("iconcom_3D_edges")
!     timer_iconcom_3D_verts  = new_timer  ("iconcom_3D_verts")
!     timer_iconcom_3D_all    = new_timer  ("iconcom_3D_all")
!     timer_iconcom_3D_comb   = new_timer  ("iconcom_3D_comb")
!     timer_iconcom_3D_keep   = new_timer  ("iconcom_3D_keep")
!
!     ! test the 3D iconcom on cells
!     CALL work_mpi_barrier()
!     CALL timer_start(timer_iconcom_3D_cells)
!     DO i=1,testbed_iterations
!       CALL icon_comm_sync(pnt_3D_cells_1, cells_not_in_domain)
!       CALL do_calculations()
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_cells)
!
!     CALL work_mpi_barrier()
!
!     ! test the 3D iconcom on edges
!     CALL work_mpi_barrier()
!     CALL timer_start(timer_iconcom_3D_edges)
!     DO i=1,testbed_iterations
!       CALL icon_comm_sync(pnt_3D_edges_1, edges_not_owned)
!       CALL do_calculations()
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_edges)
!
!     ! test the 3D iconcom on verts
!     CALL timer_start(timer_iconcom_3D_verts)
!     DO i=1,testbed_iterations
!       CALL icon_comm_sync(pnt_3D_verts_1, verts_not_owned)
!       CALL do_calculations()
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_verts)
!
!     ! test the 3D iconcom on all
!     CALL work_mpi_barrier()
!     CALL timer_start(timer_iconcom_3D_all)
!     DO i=1,testbed_iterations
!        CALL icon_comm_sync(pnt_3D_cells_1, cells_not_in_domain)
!        CALL icon_comm_sync(pnt_3D_edges_1, edges_not_owned)
!        CALL icon_comm_sync(pnt_3D_verts_1, verts_not_owned)
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_all)
!
!     ! test the 3D iconcom on combined
!     CALL work_mpi_barrier()
!     CALL timer_start(timer_iconcom_3D_comb)
!     DO i=1,testbed_iterations
!
!       comm_1 = new_icon_comm_variable(pnt_3D_cells_1, cells_not_in_domain, &
!         & p_patch(patch_no), status=is_ready, scope=until_sync)
!
!       comm_2 = new_icon_comm_variable(pnt_3D_edges_1, &
!         & edges_not_owned, p_patch(patch_no), status=is_ready, scope=until_sync)
!
!       comm_3 = new_icon_comm_variable(pnt_3D_verts_1, &
!         & verts_not_owned, p_patch(patch_no), status=is_ready, scope=until_sync)
!
!       CALL icon_comm_sync_all()
!
!       CALL do_calculations()
!       CALL do_calculations()
!       CALL do_calculations()
!
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_comb)
!
!     ! test the 3D iconcom on keeping the communicators
!     comm_1 = new_icon_comm_variable(pnt_3D_cells_1, cells_not_in_domain, &
!       & p_patch(patch_no))
!     comm_2 = new_icon_comm_variable(pnt_3D_edges_1, &
!       & edges_not_owned)
!     comm_3 = new_icon_comm_variable(pnt_3D_verts_1, &
!       & verts_not_owned)
!
!     CALL work_mpi_barrier()
!     CALL timer_start(timer_iconcom_3D_keep)
!     DO i=1,testbed_iterations
!
!       CALL icon_comm_var_is_ready(comm_1)
!       CALL icon_comm_var_is_ready(comm_2)
!       CALL icon_comm_var_is_ready(comm_3)
!       CALL icon_comm_sync_all()
!
!       CALL do_calculations()
!       CALL do_calculations()
!       CALL do_calculations()
!
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_keep)
!     CALL delete_icon_comm_variable(comm_3)
!     CALL delete_icon_comm_variable(comm_2)
!     CALL delete_icon_comm_variable(comm_1)
    !---------------------------------------------------------------------


  SUBROUTINE gather_communication_testbed()

    INTEGER :: i, j

    INTEGER :: local_size, p_n_intercomm_remote, global_size
    INTEGER, ALLOCATABLE :: owner_local(:), glb_index(:)
    TYPE(t_comm_gather_pattern) :: gather_pattern
    TYPE(t_comm_allgather_pattern) :: allgather_pattern
    LOGICAL :: disable_consistency_check

    INTEGER :: nlev
    REAL(wp), ALLOCATABLE :: in_array_r_1d(:,:), in_array_r_2d(:,:,:)
    INTEGER, ALLOCATABLE :: in_array_i_1d(:,:), in_array_i_2d(:,:,:)
    REAL(wp), ALLOCATABLE :: out_array_r_1d(:), out_array_r_2d(:,:)
    REAL(wp), ALLOCATABLE :: ref_out_array_r_1d(:), ref_out_array_r_2d(:,:)
    INTEGER, ALLOCATABLE :: out_array_i_1d(:), out_array_i_2d(:,:)
    INTEGER, ALLOCATABLE :: ref_out_array_i_1d(:), ref_out_array_i_2d(:,:)
    REAL(wp) :: fill_value
    INTEGER :: p_comm_work_backup, p_pe_work_backup, p_n_work_backup
    INTEGER :: intercomm, p_comm_work_new, ierr
    LOGICAL :: is_active
    

    CHARACTER(*), PARAMETER :: method_name = &
      "mo_test_communication:gather_communication_testbed"

    !---------------------------------------------------------------------------
    ! simple test in which each process has its own local contiguous part of the
    ! global array
    !---------------------------------------------------------------------------
    ! generate gather pattern
    local_size = 10 * nproma
    global_size = p_n_work * local_size
    ALLOCATE(owner_local(local_size), glb_index(local_size))
    DO i = 1, local_size
      owner_local(i) = p_pe_work
      glb_index(i) = local_size * p_pe_work + i
    END DO
    disable_consistency_check = .FALSE.
    CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
      &                            gather_pattern)

    ! initialise in- and reference out data
    nlev = 5
    fill_value = -1
    ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
      &      in_array_r_2d(nproma, nlev, local_size / nproma), &
      &      in_array_i_1d(nproma, local_size / nproma), &
      &      in_array_i_2d(nproma, nlev, local_size / nproma))
    DO i = 0, local_size-1
      in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      DO j = 1, nlev
        in_array_r_2d(MOD(i,nproma)+1, j, i/nproma+1) = (j - 1) * global_size + &
          &                                             p_pe_work * local_size + i
        in_array_i_2d(MOD(i,nproma)+1, j, i/nproma+1) = (j - 1) * global_size + &
          &                                             p_pe_work * local_size + i
      END DO
    END DO
    IF (p_pe_work == 0) THEN
      ALLOCATE(out_array_r_1d(global_size), &
        &      out_array_r_2d(global_size, nlev), &
        &      out_array_i_1d(global_size), &
        &      out_array_i_2d(global_size, nlev), &
        &      ref_out_array_r_1d(global_size), &
        &      ref_out_array_r_2d(global_size, nlev), &
        &      ref_out_array_i_1d(global_size), &
        &      ref_out_array_i_2d(global_size, nlev))
      DO i = 0, global_size - 1
        ref_out_array_r_1d(i+1) = i
        ref_out_array_i_1d(i+1) = i
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = (j - 1) * global_size + i
          ref_out_array_i_2d(i+1, j) = (j - 1) * global_size + i
        END DO
      END DO
    ELSE
      ALLOCATE(out_array_r_1d(0), &
        &      out_array_r_2d(0, 0), &
        &      out_array_i_1d(0), &
        &      out_array_i_2d(0, 0), &
        &      ref_out_array_r_1d(0), &
        &      ref_out_array_r_2d(0, 0), &
        &      ref_out_array_i_1d(0), &
        &      ref_out_array_i_2d(0, 0))
    END IF

    ! check gather pattern
    CALL check_exchange(in_array_r_1d, in_array_r_2d, &
      &                 in_array_i_1d, in_array_i_2d, &
      &                 out_array_r_1d, out_array_r_2d, &
      &                 ref_out_array_r_1d, ref_out_array_r_2d, &
      &                 out_array_i_1d, out_array_i_2d, &
      &                 ref_out_array_i_1d, ref_out_array_i_2d, &
      &                 gather_pattern, __LINE__)
    CALL check_exchange(in_array_r_1d, in_array_r_2d, &
      &                 in_array_i_1d, in_array_i_2d, &
      &                 out_array_r_1d, out_array_r_2d, &
      &                 ref_out_array_r_1d, ref_out_array_r_2d, &
      &                 out_array_i_1d, out_array_i_2d, &
      &                 ref_out_array_i_1d, ref_out_array_i_2d, &
      &                 gather_pattern, __LINE__, fill_value)

    ! delete gather pattern and other arrays
    DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
    DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
      &        ref_out_array_r_1d, ref_out_array_r_2d, &
      &        ref_out_array_i_1d, ref_out_array_i_2d)
    CALL delete_comm_gather_pattern(gather_pattern)
    DEALLOCATE(owner_local, glb_index)

    !---------------------------------------------------------------------------
    ! simple test in which each process has its own local contiguous part of the
    ! global array plus some overlap with neighbouring processes
    !---------------------------------------------------------------------------
    ! generate gather pattern
    local_size = 12 * nproma
    global_size = p_n_work * 10 * nproma
    ALLOCATE(owner_local(local_size), glb_index(local_size))
    DO i = 1, local_size
      owner_local(i) = p_pe_work
      glb_index(i) =  &
        MOD(global_size + (p_pe_work * 10 - 1) * nproma + i - 1, global_size) + 1
    END DO
    DO i = 1, nproma
      owner_local(i) = MOD(p_n_work + p_pe_work - 1, p_n_work)
      owner_local(11 * nproma + i) = MOD(p_pe_work + 1, p_n_work)
    END DO
    disable_consistency_check = .FALSE.
    CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
      &                            gather_pattern)
    ! initialise in- and reference out data
    nlev = 5
    fill_value = -1
    ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
      &      in_array_r_2d(nproma, nlev, local_size / nproma), &
      &      in_array_i_1d(nproma, local_size / nproma), &
      &      in_array_i_2d(nproma, nlev, local_size / nproma))

    DO i = 0, local_size-1
      in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = &
        MOD(global_size + (p_pe_work * 10 - 1) * nproma + i, global_size)
      in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = &
        MOD(global_size + (p_pe_work * 10 - 1) * nproma + i, global_size)
      DO j = 1, nlev
        in_array_r_2d(MOD(i,nproma)+1, j, i/nproma+1) = &
          in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) + (j - 1) * global_size
        in_array_i_2d(MOD(i,nproma)+1, j, i/nproma+1) = &
          in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) + (j - 1) * global_size
      END DO
    END DO
    IF (p_pe_work == 0) THEN
      ALLOCATE(out_array_r_1d(global_size), &
        &      out_array_r_2d(global_size, nlev), &
        &      out_array_i_1d(global_size), &
        &      out_array_i_2d(global_size, nlev), &
        &      ref_out_array_r_1d(global_size), &
        &      ref_out_array_r_2d(global_size, nlev), &
        &      ref_out_array_i_1d(global_size), &
        &      ref_out_array_i_2d(global_size, nlev))
      DO i = 0, global_size - 1
        ref_out_array_r_1d(i+1) = i
        ref_out_array_i_1d(i+1) = i
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = (j - 1) * global_size + i
          ref_out_array_i_2d(i+1, j) = (j - 1) * global_size + i
        END DO
      END DO
    ELSE
      ALLOCATE(out_array_r_1d(0), &
        &      out_array_r_2d(0, 0), &
        &      out_array_i_1d(0), &
        &      out_array_i_2d(0, 0), &
        &      ref_out_array_r_1d(0), &
        &      ref_out_array_r_2d(0, 0), &
        &      ref_out_array_i_1d(0), &
        &      ref_out_array_i_2d(0, 0))
    END IF

    ! check gather pattern
    CALL check_exchange(in_array_r_1d, in_array_r_2d, &
      &                 in_array_i_1d, in_array_i_2d, &
      &                 out_array_r_1d, out_array_r_2d, &
      &                 ref_out_array_r_1d, ref_out_array_r_2d, &
      &                 out_array_i_1d, out_array_i_2d, &
      &                 ref_out_array_i_1d, ref_out_array_i_2d, &
      &                 gather_pattern, __LINE__)
    CALL check_exchange(in_array_r_1d, in_array_r_2d, &
      &                 in_array_i_1d, in_array_i_2d, &
      &                 out_array_r_1d, out_array_r_2d, &
      &                 ref_out_array_r_1d, ref_out_array_r_2d, &
      &                 out_array_i_1d, out_array_i_2d, &
      &                 ref_out_array_i_1d, ref_out_array_i_2d, &
      &                 gather_pattern, __LINE__, fill_value)

    ! delete gather pattern and other arrays
    DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
    DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
      &        ref_out_array_r_1d, ref_out_array_r_2d, &
      &        ref_out_array_i_1d, ref_out_array_i_2d)
    CALL delete_comm_gather_pattern(gather_pattern)
    DEALLOCATE(owner_local, glb_index)

    !---------------------------------------------------------------------------
    ! test in which only odd numbered global indices are owned by processes
    !---------------------------------------------------------------------------
    ! generate gather pattern
    local_size = 10 * nproma
    global_size = p_n_work * local_size
    ALLOCATE(owner_local(local_size), glb_index(local_size))
    DO i = 1, local_size, 2
      owner_local(i) = p_pe_work
    END DO
    DO i = 2, local_size, 2
      owner_local(i) = -1
    END DO
    DO i = 1, local_size
      glb_index(i) = local_size * p_pe_work + i
    END DO
    disable_consistency_check = .TRUE.
    CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
      &                            gather_pattern, disable_consistency_check)

    ! initialise in- and reference out data
    nlev = 5
    fill_value = -1
    ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
      &      in_array_r_2d(nproma, nlev, local_size / nproma), &
      &      in_array_i_1d(nproma, local_size / nproma), &
      &      in_array_i_2d(nproma, nlev, local_size / nproma))
    DO i = 0, local_size-1
      in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      DO j = 1, nlev
        in_array_r_2d(MOD(i,nproma)+1, j, i/nproma+1) = (j - 1) * global_size + &
          &                                             p_pe_work * local_size + i
        in_array_i_2d(MOD(i,nproma)+1, j, i/nproma+1) = (j - 1) * global_size + &
          &                                             p_pe_work * local_size + i
      END DO
    END DO
    IF (p_pe_work == 0) THEN
      ALLOCATE(out_array_r_1d(global_size), &
        &      out_array_r_2d(global_size, nlev), &
        &      out_array_i_1d(global_size), &
        &      out_array_i_2d(global_size, nlev), &
        &      ref_out_array_r_1d(global_size), &
        &      ref_out_array_r_2d(global_size, nlev), &
        &      ref_out_array_i_1d(global_size), &
        &      ref_out_array_i_2d(global_size, nlev))
      DO i = 0, global_size - 1, 2
        ref_out_array_r_1d(i+1) = i
        ref_out_array_i_1d(i+1) = i
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = (j - 1) * global_size + i
          ref_out_array_i_2d(i+1, j) = (j - 1) * global_size + i
        END DO
      END DO
      DO i = 1, global_size - 1, 2
        ref_out_array_r_1d(i+1) = fill_value
        ref_out_array_i_1d(i+1) = fill_value
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = fill_value
          ref_out_array_i_2d(i+1, j) = fill_value
        END DO
      END DO
    ELSE
      ALLOCATE(out_array_r_1d(0), &
        &      out_array_r_2d(0, 0), &
        &      out_array_i_1d(0), &
        &      out_array_i_2d(0, 0), &
        &      ref_out_array_r_1d(0), &
        &      ref_out_array_r_2d(0, 0), &
        &      ref_out_array_i_1d(0), &
        &      ref_out_array_i_2d(0, 0))
    END IF

    ! check gather pattern
    CALL check_exchange(in_array_r_1d, in_array_r_2d, &
      &                 in_array_i_1d, in_array_i_2d, &
      &                 out_array_r_1d, out_array_r_2d, &
      &                 ref_out_array_r_1d, ref_out_array_r_2d, &
      &                 out_array_i_1d, out_array_i_2d, &
      &                 ref_out_array_i_1d, ref_out_array_i_2d, &
      &                 gather_pattern, __LINE__, fill_value)

    ! delete gather pattern and other arrays
    DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
    DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
      &        ref_out_array_r_1d, ref_out_array_r_2d, &
      &        ref_out_array_i_1d, ref_out_array_i_2d)
    CALL delete_comm_gather_pattern(gather_pattern)
    DEALLOCATE(owner_local, glb_index)

    !---------------------------------------------------------------------------
    ! simple test in which each process has its own local contiguous part of the
    ! global array plus some overlap with neighbouring processes
    ! only odd numbered global indices are owned by processes
    !---------------------------------------------------------------------------
    ! generate gather pattern
    local_size = 12 * nproma
    global_size = p_n_work * 10 * nproma
    ALLOCATE(owner_local(local_size), glb_index(local_size))
    DO i = 1, local_size
      owner_local(i) = p_pe_work
      glb_index(i) =  &
        MOD(global_size + (p_pe_work * 10 - 1) * nproma + i - 1, global_size) + 1
    END DO
    DO i = 1, nproma
      owner_local(i) = MOD(p_n_work + p_pe_work - 1, p_n_work)
      owner_local(11 * nproma + i) = MOD(p_pe_work + 1, p_n_work)
    END DO
    DO i = 2, local_size, 2
      owner_local(i) = -1
    END DO
    disable_consistency_check = .TRUE.
    CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
      &                            gather_pattern, disable_consistency_check)

    ! initialise in- and reference out data
    nlev = 5
    fill_value = -1
    ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
      &      in_array_r_2d(nproma, nlev, local_size / nproma), &
      &      in_array_i_1d(nproma, local_size / nproma), &
      &      in_array_i_2d(nproma, nlev, local_size / nproma))

    DO i = 0, local_size-1
      in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = &
        MOD(global_size + (p_pe_work * 10 - 1) * nproma + i, global_size)
      in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = &
        MOD(global_size + (p_pe_work * 10 - 1) * nproma + i, global_size)
      DO j = 1, nlev
        in_array_r_2d(MOD(i,nproma)+1, j, i/nproma+1) = &
          in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) + (j - 1) * global_size
        in_array_i_2d(MOD(i,nproma)+1, j, i/nproma+1) = &
          in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) + (j - 1) * global_size
      END DO
    END DO
    IF (p_pe_work == 0) THEN
      ALLOCATE(out_array_r_1d(global_size), &
        &      out_array_r_2d(global_size, nlev), &
        &      out_array_i_1d(global_size), &
        &      out_array_i_2d(global_size, nlev), &
        &      ref_out_array_r_1d(global_size), &
        &      ref_out_array_r_2d(global_size, nlev), &
        &      ref_out_array_i_1d(global_size), &
        &      ref_out_array_i_2d(global_size, nlev))
      DO i = 0, global_size - 1, 2
        ref_out_array_r_1d(i+1) = i
        ref_out_array_i_1d(i+1) = i
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = (j - 1) * global_size + i
          ref_out_array_i_2d(i+1, j) = (j - 1) * global_size + i
        END DO
      END DO
      DO i = 1, global_size - 1, 2
        ref_out_array_r_1d(i+1) = fill_value
        ref_out_array_i_1d(i+1) = fill_value
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = fill_value
          ref_out_array_i_2d(i+1, j) = fill_value
        END DO
      END DO
    ELSE
      ALLOCATE(out_array_r_1d(0), &
        &      out_array_r_2d(0, 0), &
        &      out_array_i_1d(0), &
        &      out_array_i_2d(0, 0), &
        &      ref_out_array_r_1d(0), &
        &      ref_out_array_r_2d(0, 0), &
        &      ref_out_array_i_1d(0), &
        &      ref_out_array_i_2d(0, 0))
    END IF

    ! check gather pattern
    CALL check_exchange(in_array_r_1d, in_array_r_2d, &
      &                 in_array_i_1d, in_array_i_2d, &
      &                 out_array_r_1d, out_array_r_2d, &
      &                 ref_out_array_r_1d, ref_out_array_r_2d, &
      &                 out_array_i_1d, out_array_i_2d, &
      &                 ref_out_array_i_1d, ref_out_array_i_2d, &
      &                 gather_pattern, __LINE__, fill_value)

    ! initialise reference out data (for no fill_value case)
    DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
      &        ref_out_array_r_1d, ref_out_array_r_2d, &
      &        ref_out_array_i_1d, ref_out_array_i_2d)
    IF (p_pe_work == 0) THEN
      ALLOCATE(out_array_r_1d(global_size/2), &
        &      out_array_r_2d(global_size/2, nlev), &
        &      out_array_i_1d(global_size/2), &
        &      out_array_i_2d(global_size/2, nlev), &
        &      ref_out_array_r_1d(global_size/2), &
        &      ref_out_array_r_2d(global_size/2, nlev), &
        &      ref_out_array_i_1d(global_size/2), &
        &      ref_out_array_i_2d(global_size/2, nlev))
      DO i = 1, global_size/2
        ref_out_array_r_1d(i) = 2 * (i - 1)
        ref_out_array_i_1d(i) = 2 * (i - 1)
        DO j = 1, nlev
          ref_out_array_r_2d(i, j) = (j - 1) * global_size + 2 * (i - 1)
          ref_out_array_i_2d(i, j) = (j - 1) * global_size + 2 * (i - 1)
        END DO
      END DO
    ELSE
      ALLOCATE(out_array_r_1d(0), &
        &      out_array_r_2d(0, 0), &
        &      out_array_i_1d(0), &
        &      out_array_i_2d(0, 0), &
        &      ref_out_array_r_1d(0), &
        &      ref_out_array_r_2d(0, 0), &
        &      ref_out_array_i_1d(0), &
        &      ref_out_array_i_2d(0, 0))
    END IF

    ! check gather pattern
    CALL check_exchange(in_array_r_1d, in_array_r_2d, &
      &                 in_array_i_1d, in_array_i_2d, &
      &                 out_array_r_1d, out_array_r_2d, &
      &                 ref_out_array_r_1d, ref_out_array_r_2d, &
      &                 out_array_i_1d, out_array_i_2d, &
      &                 ref_out_array_i_1d, ref_out_array_i_2d, &
      &                 gather_pattern, __LINE__)

    ! delete gather pattern and other arrays
    DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
    DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
      &        ref_out_array_r_1d, ref_out_array_r_2d, &
      &        ref_out_array_i_1d, ref_out_array_i_2d)
    CALL delete_comm_gather_pattern(gather_pattern)
    DEALLOCATE(owner_local, glb_index)

#ifndef NOMPI
    !---------------------------------------------------------------------------
    ! initial setup for allgather_intercomm tests
    !---------------------------------------------------------------------------
    p_comm_work_backup = p_comm_work
    p_pe_work_backup = p_pe_work
    p_n_work_backup = p_n_work
    CALL MPI_Comm_split(p_comm_work, p_pe_work/((p_n_work*2)/3), p_pe_work, &
      &                 p_comm_work_new, ierr)
    CALL MPI_Intercomm_create(p_comm_work_new, 0, p_comm_work, &
      MERGE((p_n_work*2)/3, 0, p_pe_work < (p_n_work*2)/3), 2, intercomm, ierr)
    p_comm_work = p_comm_work_new
    p_pe_work = p_comm_rank(p_comm_work_new)
    p_n_work = p_comm_size(p_comm_work_new)
    CALL MPI_Comm_remote_size(intercomm, p_n_intercomm_remote, ierr)
    IF (p_n_intercomm_remote /= p_n_work_backup - p_n_work) &
      CALL finish(method_name, "problem with intercomm")

    !---------------------------------------------------------------------------
    ! simple allgather_intercomm test in which each process has its own local
    ! contiguous part of the global array
    !---------------------------------------------------------------------------
    ! generate gather pattern
    local_size = 10 * nproma
    global_size = p_n_work * local_size
    ALLOCATE(owner_local(local_size), glb_index(local_size))
    DO i = 1, local_size
      owner_local(i) = p_pe_work
      glb_index(i) = local_size * p_pe_work + i
    END DO
    disable_consistency_check = .FALSE.
    CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
      &                            gather_pattern)
    CALL setup_comm_allgather_pattern(gather_pattern, intercomm, &
      &                               allgather_pattern)

    ! initialise in- and reference out data
    nlev = 5
    fill_value = -1
    ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
      &      in_array_i_1d(nproma, local_size / nproma))
    DO i = 0, local_size-1
      in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
    END DO
    global_size = p_n_intercomm_remote * local_size
    ALLOCATE(out_array_r_1d(global_size), &
      &      out_array_i_1d(global_size), &
      &      ref_out_array_r_1d(global_size), &
      &      ref_out_array_i_1d(global_size))
    DO i = 0, global_size - 1
      ref_out_array_r_1d(i+1) = i
      ref_out_array_i_1d(i+1) = i
    END DO

    ! check gather pattern
    CALL check_exchange_allgather(in_array_r_1d, in_array_i_1d, &
      &                           out_array_r_1d, ref_out_array_r_1d, &
      &                           out_array_i_1d, ref_out_array_i_1d, &
      &                           allgather_pattern, __LINE__)
    CALL check_exchange_allgather(in_array_r_1d, in_array_i_1d, &
      &                           out_array_r_1d, ref_out_array_r_1d, &
      &                           out_array_i_1d, ref_out_array_i_1d, &
      &                           allgather_pattern, __LINE__, fill_value)
!
    ! delete gather pattern and other arrays
    DEALLOCATE(in_array_r_1d, in_array_i_1d, &
      &        out_array_r_1d, out_array_i_1d, &
      &        ref_out_array_r_1d, ref_out_array_i_1d)
    CALL delete_comm_allgather_pattern(allgather_pattern)
    CALL delete_comm_gather_pattern(gather_pattern)
    DEALLOCATE(owner_local, glb_index)

    !---------------------------------------------------------------------------
    ! simple allgather_intercomm test in which each process has its own local
    ! contiguous part of the global array (only one side of the intercomm has
    ! data)
    !---------------------------------------------------------------------------
    do j = 0, 1

      is_active = (p_n_work < (p_n_work_backup / 2)) .EQV. (j == 0)

      ! generate gather pattern
      local_size = MERGE(10 * nproma, 0, is_active)
      global_size = p_n_work * local_size
      ALLOCATE(owner_local(local_size), glb_index(local_size))
      DO i = 1, local_size
        owner_local(i) = p_pe_work
        glb_index(i) = local_size * p_pe_work + i
      END DO
      disable_consistency_check = .FALSE.
      CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
        &                            gather_pattern)
      CALL setup_comm_allgather_pattern(gather_pattern, intercomm, &
        &                               allgather_pattern)

      ! initialise in- and reference out data
      nlev = 5
      fill_value = -1
      ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
        &      in_array_i_1d(nproma, local_size / nproma))
      DO i = 0, local_size-1
        in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
        in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      END DO
      global_size = p_n_intercomm_remote * MERGE(0, 10 * nproma, is_active)
      ALLOCATE(out_array_r_1d(global_size), &
        &      out_array_i_1d(global_size), &
        &      ref_out_array_r_1d(global_size), &
        &      ref_out_array_i_1d(global_size))
      DO i = 0, global_size - 1
        ref_out_array_r_1d(i+1) = i
        ref_out_array_i_1d(i+1) = i
      END DO

      ! check gather pattern
      CALL check_exchange_allgather(in_array_r_1d, in_array_i_1d, &
        &                           out_array_r_1d, ref_out_array_r_1d, &
        &                           out_array_i_1d, ref_out_array_i_1d, &
        &                           allgather_pattern, __LINE__)
      CALL check_exchange_allgather(in_array_r_1d, in_array_i_1d, &
        &                           out_array_r_1d, ref_out_array_r_1d, &
        &                           out_array_i_1d, ref_out_array_i_1d, &
        &                           allgather_pattern, __LINE__, fill_value)
  !
      ! delete gather pattern and other arrays
      DEALLOCATE(in_array_r_1d, in_array_i_1d, &
        &        out_array_r_1d, out_array_i_1d, &
        &        ref_out_array_r_1d, ref_out_array_i_1d)
      CALL delete_comm_allgather_pattern(allgather_pattern)
      CALL delete_comm_gather_pattern(gather_pattern)
      DEALLOCATE(owner_local, glb_index)
    END DO

    !---------------------------------------------------------------------------
    ! clean up allgather_intercomm stuff
    !---------------------------------------------------------------------------

    CALL MPI_Comm_free(intercomm, ierr)
    CALL MPI_Comm_free(p_comm_work, ierr)

    p_comm_work = p_comm_work_backup
    p_pe_work = p_pe_work_backup
    p_n_work = p_n_work_backup
#endif

  CONTAINS

    SUBROUTINE check_exchange(in_array_r_1d, in_array_r_2d, &
      &                       in_array_i_1d, in_array_i_2d, &
      &                       out_array_r_1d, out_array_r_2d, &
      &                       ref_out_array_r_1d, ref_out_array_r_2d, &
      &                       out_array_i_1d, out_array_i_2d, &
      &                       ref_out_array_i_1d, ref_out_array_i_2d, &
      &                       gather_pattern, line, fill_value)

      REAL(wp), INTENT(IN) :: in_array_r_1d(:,:), in_array_r_2d(:,:,:)
      INTEGER, INTENT(IN) :: in_array_i_1d(:,:), in_array_i_2d(:,:,:)
      REAL(wp), INTENT(OUT) :: out_array_r_1d(:), out_array_r_2d(:,:)
      REAL(wp), INTENT(IN) :: ref_out_array_r_1d(:), ref_out_array_r_2d(:,:)
      INTEGER, INTENT(OUT) :: out_array_i_1d(:), out_array_i_2d(:,:)
      INTEGER, INTENT(IN) :: ref_out_array_i_1d(:), ref_out_array_i_2d(:,:)
      TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
      INTEGER, INTENT(in) :: line
      REAL(wp), OPTIONAL, INTENT(IN) :: fill_value
      INTEGER :: i

      CALL exchange_data(in_array=in_array_r_1d, out_array=out_array_r_1d, &
        &                gather_pattern=gather_pattern, fill_value=fill_value)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_r_1d /= ref_out_array_r_1d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result r_1d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
      CALL exchange_data(in_array=in_array_r_2d, out_array=out_array_r_2d, &
        &                gather_pattern=gather_pattern, fill_value=fill_value)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_r_2d /= ref_out_array_r_2d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result r_2d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
      CALL exchange_data(in_array=in_array_i_1d, out_array=out_array_i_1d, &
        &                gather_pattern=gather_pattern, fill_value=fill_value)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_i_1d /= ref_out_array_i_1d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result i_1d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
      CALL exchange_data(in_array=in_array_i_2d, out_array=out_array_i_2d, &
        &                gather_pattern=gather_pattern, fill_value=fill_value)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_i_2d /= ref_out_array_i_2d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result i_2d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
    END SUBROUTINE

    SUBROUTINE check_exchange_allgather(in_array_r_1d, &
      &                                 in_array_i_1d, &
      &                                 out_array_r_1d, &
      &                                 ref_out_array_r_1d, &
      &                                 out_array_i_1d, &
      &                                 ref_out_array_i_1d, &
      &                                 allgather_pattern, line, &
      &                                 fill_value)

      REAL(wp), INTENT(IN) :: in_array_r_1d(:,:)
      INTEGER, INTENT(IN) :: in_array_i_1d(:,:)
      REAL(wp), INTENT(OUT) :: out_array_r_1d(:)
      REAL(wp), INTENT(IN) :: ref_out_array_r_1d(:)
      INTEGER, INTENT(OUT) :: out_array_i_1d(:)
      INTEGER, INTENT(IN) :: ref_out_array_i_1d(:)
      TYPE(t_comm_allgather_pattern), INTENT(IN) :: allgather_pattern
      INTEGER, INTENT(in) :: line
      REAL(wp), OPTIONAL, INTENT(IN) :: fill_value
      INTEGER :: i

      CALL exchange_data(in_array=in_array_r_1d, out_array=out_array_r_1d, &
        &                allgather_pattern=allgather_pattern, &
        &                fill_value=fill_value)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_r_1d /= ref_out_array_r_1d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result r_1d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
      CALL exchange_data(in_array=in_array_i_1d, out_array=out_array_i_1d, &
        &                allgather_pattern=allgather_pattern, &
        &                fill_value=fill_value)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_i_1d /= ref_out_array_i_1d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result i_1d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
    END SUBROUTINE

  END SUBROUTINE

END MODULE mo_test_communication

