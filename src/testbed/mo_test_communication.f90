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
  USE mo_communication,       ONLY: t_comm_gather_pattern, t_comm_pattern, &
    &                               setup_comm_pattern, delete_comm_pattern, &
    &                               setup_comm_gather_pattern, &
    &                               delete_comm_gather_pattern, exchange_data, &
    &                               exchange_data_4de1, &
    &                               exchange_data_mult, &
    &                               t_comm_allgather_pattern, &
    &                               setup_comm_allgather_pattern, &
    &                               delete_comm_allgather_pattern, &
    &                               exchange_data_grf, &
    &                               t_comm_pattern_collection, &
    &                               setup_comm_pattern_collection, &
    &                               delete_comm_pattern_collection
  USE mo_decomposition_tools, ONLY: t_glb2loc_index_lookup, &
    &                               init_glb2loc_index_lookup, &
    &                               set_inner_glb_index, &
    &                               deallocate_glb2loc_index_lookup

  USE mo_master_control,      ONLY: get_my_process_name

  USE mo_model_domain,        ONLY: p_patch
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model
  USE mo_icoham_dyn_memory,   ONLY: p_hydro_state
  USE mo_math_gradients,      ONLY: grad_fd_norm

  !nh utils
  USE mo_nonhydro_state,      ONLY: p_nh_state
  USE mo_atmo_nonhydrostatic, ONLY: construct_atmo_nonhydrostatic, destruct_atmo_nonhydrostatic
  USE mo_async_latbc_types,   ONLY: t_latbc_data

  USE mo_parallel_config,    ONLY: itype_comm, iorder_sendrecv
  USE mo_sync,               ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, &
    &                              sync_patch_array_mult, sync_patch_array_4de1
  USE mo_icon_comm_lib

  USE mo_rrtm_data_interface, ONLY: t_rrtm_data, recv_rrtm_input, &
    & init_rrtm_model_repart

  USE mo_icon_testbed_config, ONLY: testbed_model, test_halo_communication, &
    & test_radiation_communication, testbed_iterations, calculate_iterations, &
    & test_gather_communication, test_exchange_communication


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

    TYPE(t_latbc_data) :: latbc !< data structure for async latbc prefetching

    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_communication"


    !---------------------------------------------------------------------

    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name

    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    CALL construct_atmo_nonhydrostatic(latbc)
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

    CASE(test_exchange_communication)
      CALL exchange_communication_testbed()
      CALL exchange_communication_grf_testbed()

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
    CALL destruct_atmo_nonhydrostatic(latbc)
    CALL destruct_atmo_model()
!     CALL destruct_icon_communication()
    CALL message(TRIM(method_name),'clean-up finished')



  END SUBROUTINE test_communication
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE radiation_communication_testbed()

    TYPE(t_rrtm_data), TARGET :: rrtm_local_data
    TYPE(t_rrtm_data), POINTER :: rrtm_rad_data
    INTEGER ::  timer_barrier_only, i

    CALL init_rrtm_model_repart()
!     ! now allocate the data for the radiation interface
!     CALL init_rrtm_data( &
!       & rrtm_data   = rrtm_local_data   , &
!       & no_of_cells = p_patch(1)%n_patch_cells, &
!       & full_levels = p_patch(1)%nlev,         &
!       & half_levels = p_patch(1)%nlevp1,       &
!       & block_size  = nproma)

    !---------------------------------------------------------------------
    ! test the barrier
    !---------------------------------------------------------------------
    timer_barrier_only  = new_timer("mpi_barrier_only")
    CALL work_mpi_barrier()
    DO i=1,testbed_iterations
      CALL timer_start(timer_barrier_only)
      CALL work_mpi_barrier()
      CALL timer_stop(timer_barrier_only)
    ENDDO
    CALL work_mpi_barrier()

    DO i=1,testbed_iterations
      CALL timer_start(timer_radiaton_recv)
      CALL recv_rrtm_input( &
      ktype       = rrtm_local_data%convection_type,&!< in     type of convection
      zland       = rrtm_local_data%fr_land_smt    ,&  !< in     land fraction
      zglac       = rrtm_local_data%fr_glac_smt    ,&  !< in     land glacier fraction
      cos_mu0     = rrtm_local_data%cosmu0         ,&  !< in  cos of zenith angle mu0
 !     alb_vis_dir = rrtm_local_data%albedo_vis_dir ,&  !< in surface albedo for visible range, direct
 !     alb_nir_dir = rrtm_local_data%albedo_nir_dir ,&  !< in surface albedo for near IR range, direct
      alb_vis_dif = rrtm_local_data%albedo_vis_dif ,&  !< in surface albedo for visible range, diffuse
 !     alb_nir_dif = rrtm_local_data%albedo_nir_dif ,&  !< in surface albedo for near IR range, diffuse
      emis_rad    = rrtm_local_data%emis_rad       ,&  !< in longwave surface emissivity
      tk_sfc      = rrtm_local_data%tsfctrad       ,&  !< in surface temperature
      pp_hl       = rrtm_local_data%pres_ifc       ,&  !< in  pres at half levels at t-dt [Pa]
      pp_fl       = rrtm_local_data%pres           ,&  !< in  pres at full levels at t-dt [Pa]
      tk_fl       = rrtm_local_data%temp           ,&  !< in  temperature at full level at t-dt
      qm_vap      = rrtm_local_data%qm_vapor       ,&  !< in  water vapor mass mix ratio at t-dt
      qm_liq      = rrtm_local_data%qm_liquid      ,&  !< in cloud water mass mix ratio at t-dt
      qm_ice      = rrtm_local_data%qm_ice         ,&  !< in cloud ice mass mixing ratio at t-dt
      qm_o3       = rrtm_local_data%qm_o3          ,&  !< in o3 mass mixing ratio at t-dt
      cdnc        = rrtm_local_data%acdnc          ,&  !< in  cloud droplet numb conc. [1/m**3]
      cld_frc     = rrtm_local_data%cld_frc        ,&  !< in  cloud fraction [m2/m2]
      zaeq1       = rrtm_local_data%zaeq1          ,&  !< in aerosol continental
      zaeq2       = rrtm_local_data%zaeq2          ,&  !< in aerosol maritime
      zaeq3       = rrtm_local_data%zaeq3          ,&  !< in aerosol urban
      zaeq4       = rrtm_local_data%zaeq4          ,&  !< in aerosol volcano ashes
      zaeq5       = rrtm_local_data%zaeq5          ,&  !< in aerosol stratospheric background
      patch       = p_patch(1)     ,&  !< in
      rrtm_data   = rrtm_rad_data )  !< pointer out
      CALL timer_stop(timer_radiaton_recv)
      CALL timer_start(timer_barrier)
      CALL work_mpi_barrier()
      CALL timer_stop(timer_barrier)
    ENDDO


  END SUBROUTINE radiation_communication_testbed
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
    REAL(wp), ALLOCATABLE :: pnt_4D_cells(:,:,:,:)

    INTEGER :: timer_3D_cells_1, timer_3D_cells_2, timer_3D_cells_3, &
      &        timer_3D_cells_4, timer_4DE1_cells

    INTEGER :: timer_sync_1_2_3D_cells

    INTEGER :: i, j, k, l, m, n
    CHARACTER(len=128) :: str_i


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

    DO i = 1, 16
      ALLOCATE(pnt_4D_cells(i, SIZE(pnt_3D_cells_1, 1), &
        &                      SIZE(pnt_3D_cells_1, 2), &
        &                      SIZE(pnt_3D_cells_1, 3)))

      n = 1
      DO j = 1, SIZE(pnt_4D_cells, 1)
        DO k = 1, SIZE(pnt_4D_cells, 2)
          DO l = 1, SIZE(pnt_4D_cells, 3)
            DO m = 1, SIZE(pnt_4D_cells, 4)
              pnt_4D_cells(j, k, l, m) = n
              n = n + 1
            END DO
          END DO
        END DO
      END DO
      write (str_i, '(I4)') i
      timer_4DE1_cells = new_timer(timer_descr//"_4de1cells_"//ADJUSTL(TRIM(str_i)))
      CALL test_oldsync_4DE1(SYNC_C, var=pnt_4D_cells, nfields=i, &
        &                    timer_id=timer_4DE1_cells, patch_no=patch_no)
      DEALLOCATE(pnt_4D_cells)
    END DO

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
  SUBROUTINE test_oldsync_4DE1(comm_pattern, var, nfields, timer_id, patch_no)
    INTEGER :: comm_pattern, nfields, timer_id, patch_no
    REAL(wp) :: var(:,:,:,:)

    INTEGER :: i

    CALL work_mpi_barrier()
    DO i=1,testbed_iterations
      CALL timer_start(timer_id)
      CALL sync_patch_array_4de1(comm_pattern, p_patch(patch_no), nfields, var)
      CALL timer_stop(timer_id)
    ENDDO

   CALL work_mpi_barrier()

  END SUBROUTINE test_oldsync_4DE1
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


  SUBROUTINE exchange_communication_testbed()

    INTEGER, ALLOCATABLE :: owner_local_src(:), owner_local_dst(:), &
      &                     glb_index_src(:), glb_index_dst(:)
    INTEGER :: i, j, n, local_size_src, local_size_dst, global_size
    TYPE(t_glb2loc_index_lookup) :: send_glb2loc_index
    TYPE(t_comm_pattern) :: comm_pattern

    INTEGER :: nlev
    REAL(wp), ALLOCATABLE :: in_array_r_2d(:,:), in_array_r_3d(:,:,:), &
      &                      in_array_r_4d(:,:,:,:)
    INTEGER, ALLOCATABLE ::  in_array_i_2d(:,:), in_array_i_3d(:,:,:)
    LOGICAL, ALLOCATABLE ::  in_array_l_2d(:,:), in_array_l_3d(:,:,:)
    REAL(wp), ALLOCATABLE :: out_array_r_2d(:,:), out_array_r_3d(:,:,:), &
      &                      out_array_r_4d(:,:,:,:)
    INTEGER, ALLOCATABLE ::  out_array_i_2d(:,:), out_array_i_3d(:,:,:)
    LOGICAL, ALLOCATABLE ::  out_array_l_2d(:,:), out_array_l_3d(:,:,:)
    REAL(wp), ALLOCATABLE :: add_array_r_2d(:,:), add_array_r_3d(:,:,:), &
      &                      add_array_r_4d(:,:,:,:)
    INTEGER, ALLOCATABLE ::  add_array_i_2d(:,:), add_array_i_3d(:,:,:)
    REAL(wp), ALLOCATABLE :: ref_out_array_r_2d(:,:), ref_out_array_r_3d(:,:,:), &
      &                      ref_out_array_r_4d(:,:,:,:)
    INTEGER, ALLOCATABLE ::  ref_out_array_i_2d(:,:), ref_out_array_i_3d(:,:,:)
    LOGICAL, ALLOCATABLE ::  ref_out_array_l_2d(:,:), ref_out_array_l_3d(:,:,:)

    CHARACTER(*), PARAMETER :: method_name = &
      "mo_test_communication:exchange_communication_testbed"

    ! generate communication pattern
    local_size_src = 10 * nproma
    local_size_dst = 16 * nproma
    global_size = local_size_src * p_n_work
    ALLOCATE(owner_local_src(local_size_src), glb_index_src(local_size_src), &
      &      owner_local_dst(local_size_dst), glb_index_dst(local_size_dst))

    owner_local_src(1:local_size_src:2) = -1
    owner_local_src(2:local_size_src:2) = p_pe_work

    owner_local_dst(1:local_size_dst:2) = -1
    owner_local_dst(2:3*nproma:2) = MOD(p_pe_work + p_n_work - 1, p_n_work)
    owner_local_dst(3*nproma+2:13*nproma:2) = p_pe_work
    owner_local_dst(13*nproma+2:16*nproma:2) = MOD(p_pe_work + 1, p_n_work)
    DO i = 1, local_size_src
      glb_index_src(i) = local_size_src * p_pe_work + i
    END DO
    DO i = 1, local_size_dst
      glb_index_dst(i) = MOD(local_size_src * p_pe_work - 3 * nproma - 1 + i + &
      &                      global_size, global_size) + 1
    END DO
    CALL init_glb2loc_index_lookup(send_glb2loc_index, global_size)

    CALL set_inner_glb_index(send_glb2loc_index, glb_index_src, &
      &                      (/(i, i=1, local_size_src)/))

    CALL setup_comm_pattern(local_size_dst, owner_local_dst, glb_index_dst, &
      &                     send_glb2loc_index, local_size_src, &
      &                     owner_local_src, glb_index_src, comm_pattern)

    nlev = 7

    ALLOCATE(in_array_r_2d(nproma, 10), in_array_r_3d(nproma, nlev, 10), &
      &      in_array_i_2d(nproma, 10), in_array_i_3d(nproma, nlev, 10), &
      &      in_array_l_2d(nproma, 10), in_array_l_3d(nproma, nlev, 10), &
      &      out_array_r_2d(nproma, 16), out_array_r_3d(nproma, nlev, 16), &
      &      out_array_i_2d(nproma, 16), out_array_i_3d(nproma, nlev, 16), &
      &      out_array_l_2d(nproma, 16), out_array_l_3d(nproma, nlev, 16), &
      &      add_array_r_2d(nproma, 16), add_array_r_3d(nproma, nlev, 16), &
      &      add_array_i_2d(nproma, 16), add_array_i_3d(nproma, nlev, 16), &
      &      ref_out_array_r_2d(nproma, 16), &
      &      ref_out_array_r_3d(nproma, nlev, 16), &
      &      ref_out_array_i_2d(nproma, 16), &
      &      ref_out_array_i_3d(nproma, nlev, 16), &
      &      ref_out_array_l_2d(nproma, 16), &
      &      ref_out_array_l_3d(nproma, nlev, 16))

    in_array_r_2d = RESHAPE(glb_index_src, (/nproma, 10/))
    in_array_i_2d = RESHAPE(glb_index_src, (/nproma, 10/))
    in_array_l_2d = .TRUE.
    DO i = 1, nlev
      in_array_r_3d(:,i,:) = in_array_r_2d + (i - 1) * global_size
      in_array_i_3d(:,i,:) = in_array_i_2d + (i - 1) * global_size
    END DO
    in_array_l_3d = .TRUE.

    out_array_r_2d = -1
    out_array_i_2d = -1
    out_array_l_2d = .FALSE.
    out_array_r_3d = -1
    out_array_i_3d = -1
    out_array_l_3d = .FALSE.

    ref_out_array_r_2d = RESHAPE(MERGE(-1, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 16/))
    ref_out_array_i_2d = RESHAPE(MERGE(-1, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 16/))
    ref_out_array_l_2d = RESHAPE(MERGE(.FALSE., .TRUE., owner_local_dst == -1), &
      &                          (/nproma, 16/))
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = MERGE(-1._wp, ref_out_array_r_2d + &
        &                               (i - 1) * global_size, &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                                             (/nproma, 16/)))
      ref_out_array_i_3d(:,i,:) = MERGE(-1, ref_out_array_i_2d + &
        &                               (i - 1) * global_size, &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                                             (/nproma, 16/)))
      ref_out_array_l_3d(:,i,:) = ref_out_array_l_2d
    END DO

    CALL check_exchange(in_array_r_2d = in_array_r_2d, &
      &                 in_array_r_3d = in_array_r_3d, &
      &                 in_array_i_2d = in_array_i_2d, &
      &                 in_array_i_3d = in_array_i_3d, &
      &                 in_array_l_2d = in_array_l_2d, &
      &                 in_array_l_3d = in_array_l_3d, &
      &                 out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    add_array_r_2d = MERGE(-1._wp, ref_out_array_r_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
    add_array_i_2d = MERGE(-1, ref_out_array_i_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
    DO i = 1, nlev
      add_array_r_3d(:,i,:) = MERGE(-1._wp, ref_out_array_r_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
      add_array_i_3d(:,i,:) = MERGE(-1, ref_out_array_i_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
    END DO

    ref_out_array_r_2d = MERGE(-1._wp, 2._wp * ref_out_array_r_2d, &
      &                        -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
    ref_out_array_i_2d = MERGE(-1, 2 * ref_out_array_i_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = MERGE(-1._wp, 2._wp * ref_out_array_r_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 16/)))
      ref_out_array_i_3d(:,i,:) = MERGE(-1, 2 * ref_out_array_i_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 16/)))
    END DO

    CALL check_exchange(in_array_r_2d = in_array_r_2d, &
      &                 in_array_r_3d = in_array_r_3d, &
      &                 in_array_i_2d = in_array_i_2d, &
      &                 in_array_i_3d = in_array_i_3d, &
      &                 in_array_l_2d = in_array_l_2d, &
      &                 in_array_l_3d = in_array_l_3d, &
      &                 out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 add_array_r_2d = add_array_r_2d, &
      &                 add_array_r_3d = add_array_r_3d, &
      &                 add_array_i_2d = add_array_i_2d, &
      &                 add_array_i_3d = add_array_i_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    CALL delete_comm_pattern(comm_pattern)

    ! generate communication pattern
    local_size_dst = 10 * nproma
    DEALLOCATE(owner_local_dst, glb_index_dst)
    ALLOCATE(owner_local_dst(local_size_dst), glb_index_dst(local_size_dst))

    owner_local_dst(1:local_size_dst:2) = -1
    owner_local_dst(2:7*nproma:2) = p_pe_work
    owner_local_dst(7*nproma+2:local_size_dst:2) = &
      MOD(p_pe_work + p_n_work + 1, p_n_work)
    DO i = 1, local_size_dst
      glb_index_dst(i) = MOD(local_size_src * p_pe_work + 3 * nproma - 1 + i + &
      &                      global_size, global_size) + 1
    END DO

    CALL setup_comm_pattern(local_size_dst, owner_local_dst, glb_index_dst, &
      &                     send_glb2loc_index, local_size_src, &
      &                     owner_local_src, glb_index_src, comm_pattern)

    DEALLOCATE(out_array_r_2d, out_array_r_3d, out_array_i_2d, out_array_i_3d, &
      &        out_array_l_2d, out_array_l_3d, add_array_r_2d, add_array_r_3d, &
      &        add_array_i_2d, add_array_i_3d, ref_out_array_r_2d, &
      &        ref_out_array_r_3d, ref_out_array_i_2d, ref_out_array_i_3d, &
      &        ref_out_array_l_2d, ref_out_array_l_3d)
    ALLOCATE(out_array_r_2d(nproma, 10), out_array_r_3d(nproma, nlev, 10), &
      &      out_array_i_2d(nproma, 10), out_array_i_3d(nproma, nlev, 10), &
      &      out_array_l_2d(nproma, 10), out_array_l_3d(nproma, nlev, 10), &
      &      add_array_r_2d(nproma, 10), add_array_r_3d(nproma, nlev, 10), &
      &      add_array_i_2d(nproma, 10), add_array_i_3d(nproma, nlev, 10), &
      &      ref_out_array_r_2d(nproma, 10), &
      &      ref_out_array_r_3d(nproma, nlev, 10), &
      &      ref_out_array_i_2d(nproma, 10), &
      &      ref_out_array_i_3d(nproma, nlev, 10), &
      &      ref_out_array_l_2d(nproma, 10), &
      &      ref_out_array_l_3d(nproma, nlev, 10))

    out_array_r_2d = -1
    out_array_i_2d = -1
    out_array_l_2d = .FALSE.
    out_array_r_3d = -1
    out_array_i_3d = -1
    out_array_l_3d = .FALSE.

    ref_out_array_r_2d = RESHAPE(MERGE(-1, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    ref_out_array_i_2d = RESHAPE(MERGE(-1, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    ref_out_array_l_2d = RESHAPE(MERGE(.FALSE., .TRUE., owner_local_dst == -1), &
      &                          (/nproma, 16/))
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = MERGE(-1._wp, ref_out_array_r_2d + &
        &                               (i - 1) * global_size, &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                                             (/nproma, 10/)))
      ref_out_array_i_3d(:,i,:) = MERGE(-1, ref_out_array_i_2d + &
        &                               (i - 1) * global_size, &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                                             (/nproma, 10/)))
      ref_out_array_l_3d(:,i,:) = ref_out_array_l_2d
    END DO

    CALL check_exchange(in_array_r_2d = in_array_r_2d, &
      &                 in_array_r_3d = in_array_r_3d, &
      &                 in_array_i_2d = in_array_i_2d, &
      &                 in_array_i_3d = in_array_i_3d, &
      &                 in_array_l_2d = in_array_l_2d, &
      &                 in_array_l_3d = in_array_l_3d, &
      &                 out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    add_array_r_2d = MERGE(-1._wp, ref_out_array_r_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    add_array_i_2d = MERGE(-1, ref_out_array_i_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    DO i = 1, nlev
      add_array_r_3d(:,i,:) = MERGE(-1._wp, ref_out_array_r_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
      add_array_i_3d(:,i,:) = MERGE(-1, ref_out_array_i_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    END DO

    ref_out_array_r_2d = MERGE(-1._wp, 2._wp * ref_out_array_r_2d, &
      &                        -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    ref_out_array_i_2d = MERGE(-1, 2 * ref_out_array_i_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = MERGE(-1._wp, 2._wp * ref_out_array_r_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 10/)))
      ref_out_array_i_3d(:,i,:) = MERGE(-1, 2 * ref_out_array_i_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 10/)))
    END DO

    CALL check_exchange(in_array_r_2d = in_array_r_2d, &
      &                 in_array_r_3d = in_array_r_3d, &
      &                 in_array_i_2d = in_array_i_2d, &
      &                 in_array_i_3d = in_array_i_3d, &
      &                 in_array_l_2d = in_array_l_2d, &
      &                 in_array_l_3d = in_array_l_3d, &
      &                 out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 add_array_r_2d = add_array_r_2d, &
      &                 add_array_r_3d = add_array_r_3d, &
      &                 add_array_i_2d = add_array_i_2d, &
      &                 add_array_i_3d = add_array_i_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    out_array_r_2d = in_array_r_2d
    out_array_i_2d = in_array_i_2d
    out_array_l_2d = in_array_l_2d
    out_array_r_3d = in_array_r_3d
    out_array_i_3d = in_array_i_3d
    out_array_l_3d = in_array_l_3d

    ref_out_array_r_2d = RESHAPE(MERGE(glb_index_src, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    ref_out_array_i_2d = RESHAPE(MERGE(glb_index_src, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    ref_out_array_l_2d = .TRUE.
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = ref_out_array_r_2d + (i - 1) * global_size
      ref_out_array_i_3d(:,i,:) = ref_out_array_i_2d + (i - 1) * global_size
    END DO
    ref_out_array_l_3d = .TRUE.

    CALL check_exchange(out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    out_array_r_2d = in_array_r_2d
    out_array_i_2d = in_array_i_2d
    out_array_l_2d = in_array_l_2d
    out_array_r_3d = in_array_r_3d
    out_array_i_3d = in_array_i_3d
    out_array_l_3d = in_array_l_3d

    add_array_r_2d = MERGE(-1._wp, ref_out_array_r_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    add_array_i_2d = MERGE(-1, ref_out_array_i_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    DO i = 1, nlev
      add_array_r_3d(:,i,:) = MERGE(-1._wp, ref_out_array_r_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
      add_array_i_3d(:,i,:) = MERGE(-1, ref_out_array_i_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    END DO

    ref_out_array_r_2d = RESHAPE(MERGE(glb_index_src, 2 * glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    ref_out_array_i_2d = RESHAPE(MERGE(glb_index_src, 2 * glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = MERGE(ref_out_array_r_3d(:,i,:), &
        &                               2._wp * ref_out_array_r_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 10/)))
      ref_out_array_i_3d(:,i,:) = MERGE(ref_out_array_i_3d(:,i,:), &
        &                               2 * ref_out_array_i_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 10/)))
    END DO

    CALL check_exchange(out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 add_array_r_2d = add_array_r_2d, &
      &                 add_array_r_3d = add_array_r_3d, &
      &                 add_array_i_2d = add_array_i_2d, &
      &                 add_array_i_3d = add_array_i_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    DO n = 1, 16

      ALLOCATE(in_array_r_4d(n, nproma, nlev, 10), &
        &      out_array_r_4d(n, nproma, nlev, 10), &
        &      ref_out_array_r_4d(n, nproma, nlev, 10))

      DO i = 1, nlev
        DO j = 1, n
          in_array_r_4d(j,:,i,:) = RESHAPE(glb_index_src, (/nproma, 10/)) + &
            &                              (i - 1) * global_size + 0.1_wp * j
        END DO
      END DO
      out_array_r_4d = -1
      DO i = 1, nlev
        DO j = 1, n
          ref_out_array_r_4d(j,:,i,:) = &
            RESHAPE(MERGE(-1._wp, glb_index_dst + (i - 1) * global_size + &
            &             0.1_wp * j, owner_local_dst == -1), (/nproma, 10/))
        END DO
      END DO

      CALL check_exchange_4de1(in_array=in_array_r_4d, &
        &                      out_array=out_array_r_4d, &
        &                      ref_out_array=ref_out_array_r_4d, &
        &                      comm_pattern=comm_pattern)

      out_array_r_4d = in_array_r_4d
      DO i = 1, nlev
        DO j = 1, n
          ref_out_array_r_4d(j,:,i,:) = &
            RESHAPE(MERGE(glb_index_src, glb_index_dst, &
            &             owner_local_dst == -1), (/nproma, 10/)) + &
            &             (i - 1) * global_size + 0.1_wp * j
        END DO
      END DO

      CALL check_exchange_4de1(out_array=out_array_r_4d, &
        &                      ref_out_array=ref_out_array_r_4d, &
        &                      comm_pattern=comm_pattern)

      DEALLOCATE(in_array_r_4d, out_array_r_4d, ref_out_array_r_4d)
    END DO

    ALLOCATE(in_array_r_4d(nproma,nlev,10,20), &
      &      add_array_r_4d(nproma,nlev,10,20), &
      &      out_array_r_4d(nproma,nlev,10,20), &
      &      ref_out_array_r_4d(nproma,nlev,10,20))

    DO n = 1, 20
      DO i = 1, nlev
        in_array_r_4d(:,i,:,n) = &
          RESHAPE(MERGE(-1._wp, glb_index_src + (i - 1) * global_size + &
          &             n * 0.1_wp, -1 == owner_local_src), (/nproma, 10/))
        ref_out_array_r_4d(:,i,:,n) = &
          RESHAPE(MERGE(-1._wp, glb_index_dst + (i - 1) * global_size + &
          &             n * 0.1_wp, -1 == owner_local_dst), (/nproma, 10/))
      END DO
    END DO

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      in_array4 = in_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      in_array5 = in_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array6 = out_array_r_4d(:,:,:,6), &
      in_array6 = in_array_r_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_r_4d(:,:,:,6), &
      out_array7 = out_array_r_4d(:,:,:,7), &
      in_array7 = in_array_r_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_r_4d(:,:,:,7), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      in_array4d = in_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      in_array4 = in_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      in_array5 = in_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array6 = out_array_r_4d(:,:,:,6), &
      in_array6 = in_array_r_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_r_4d(:,:,:,6), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      in_array4d = in_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      in_array4 = in_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      in_array5 = in_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      in_array4d = in_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      in_array4 = in_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      in_array4d = in_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      in_array4d = in_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      in_array4d = in_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      in_array4d = in_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      in_array4d = in_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      in_array4 = in_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      in_array5 = in_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array6 = out_array_r_4d(:,:,:,6), &
      in_array6 = in_array_r_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_r_4d(:,:,:,6), &
      out_array7 = out_array_r_4d(:,:,:,7), &
      in_array7 = in_array_r_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_r_4d(:,:,:,7), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      in_array4 = in_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      in_array5 = in_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array6 = out_array_r_4d(:,:,:,6), &
      in_array6 = in_array_r_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_r_4d(:,:,:,6), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      in_array4 = in_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      in_array5 = in_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      in_array4 = in_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      comm_pattern = comm_pattern)

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array6 = out_array_r_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_r_4d(:,:,:,6), &
      out_array7 = out_array_r_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_r_4d(:,:,:,7), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array6 = out_array_r_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_r_4d(:,:,:,6), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array6 = out_array_r_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_r_4d(:,:,:,6), &
      out_array7 = out_array_r_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_r_4d(:,:,:,7), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array6 = out_array_r_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_r_4d(:,:,:,6), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      comm_pattern = comm_pattern)

    out_array_r_4d = in_array_r_4d
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      comm_pattern = comm_pattern)

    DO n = 1, 20
      DO i = 1, nlev
        in_array_r_4d(:,i,:,n) = &
          RESHAPE(MERGE(-1._wp, glb_index_src + (i - 1) * global_size + &
          &             n * 0.1_wp, -1 == owner_local_src), (/nproma, 10/))
        ref_out_array_r_4d(:,i,:,n) = &
          RESHAPE(MERGE(-1._wp, glb_index_dst + (i - 1) * global_size + &
          &             n * 0.1_wp, -1 == owner_local_dst), (/nproma, 10/))
      END DO
    END DO

    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,:,:,1), &
      in_array1 = in_array_r_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,:,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      in_array4 = in_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      in_array5 = in_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array6 = out_array_r_4d(:,:,:,6), &
      in_array6 = in_array_r_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_r_4d(:,:,:,6), &
      out_array7 = out_array_r_4d(:,:,:,7), &
      in_array7 = in_array_r_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_r_4d(:,:,:,7), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      in_array4d = in_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      nshift = 0, comm_pattern = comm_pattern)

    ref_out_array_r_4d(:,2:,:,1) = -1._wp
    ref_out_array_r_4d(:,:5,:,2:) = -1._wp
    out_array_r_4d = -1
    CALL check_exchange_mult( &
      out_array1 = out_array_r_4d(:,1:1,:,1), &
      in_array1 = in_array_r_4d(:,1:1,:,1), &
      ref_out_array1 = ref_out_array_r_4d(:,1:1,:,1), &
      out_array2 = out_array_r_4d(:,:,:,2), &
      in_array2 = in_array_r_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_r_4d(:,:,:,2), &
      out_array3 = out_array_r_4d(:,:,:,3), &
      in_array3 = in_array_r_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_r_4d(:,:,:,3), &
      out_array4 = out_array_r_4d(:,:,:,4), &
      in_array4 = in_array_r_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_r_4d(:,:,:,4), &
      out_array5 = out_array_r_4d(:,:,:,5), &
      in_array5 = in_array_r_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_r_4d(:,:,:,5), &
      out_array6 = out_array_r_4d(:,:,:,6), &
      in_array6 = in_array_r_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_r_4d(:,:,:,6), &
      out_array7 = out_array_r_4d(:,:,:,7), &
      in_array7 = in_array_r_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_r_4d(:,:,:,7), &
      out_array4d = out_array_r_4d(:,:,:,8:), &
      in_array4d = in_array_r_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_r_4d(:,:,:,8:), &
      nshift = 5, comm_pattern = comm_pattern)

    DEALLOCATE(in_array_r_4d, add_array_r_4d, out_array_r_4d, &
      &        ref_out_array_r_4d)

    CALL delete_comm_pattern(comm_pattern)
    CALL deallocate_glb2loc_index_lookup(send_glb2loc_index)

  CONTAINS

    SUBROUTINE check_exchange(in_array_r_2d, in_array_r_3d, &
      &                       in_array_i_2d, in_array_i_3d, &
      &                       in_array_l_2d, in_array_l_3d, &
      &                       out_array_r_2d, out_array_r_3d, &
      &                       out_array_i_2d, out_array_i_3d, &
      &                       out_array_l_2d, out_array_l_3d, &
      &                       add_array_r_2d, add_array_r_3d, &
      &                       add_array_i_2d, add_array_i_3d, &
      &                       ref_out_array_r_2d, ref_out_array_r_3d, &
      &                       ref_out_array_i_2d, ref_out_array_i_3d, &
      &                       ref_out_array_l_2d, ref_out_array_l_3d, &
      &                       comm_pattern)

      REAL(wp), OPTIONAL, INTENT(IN) :: in_array_r_2d(:,:), in_array_r_3d(:,:,:)
      INTEGER, OPTIONAL, INTENT(IN) ::  in_array_i_2d(:,:), in_array_i_3d(:,:,:)
      LOGICAL, OPTIONAL, INTENT(IN) ::  in_array_l_2d(:,:), in_array_l_3d(:,:,:)
      REAL(wp), INTENT(INOUT) :: out_array_r_2d(:,:), out_array_r_3d(:,:,:)
      INTEGER, INTENT(INOUT) ::  out_array_i_2d(:,:), out_array_i_3d(:,:,:)
      LOGICAL, INTENT(INOUT) ::  out_array_l_2d(:,:), out_array_l_3d(:,:,:)
      REAL(wp), OPTIONAL, INTENT(IN) :: add_array_r_2d(:,:), add_array_r_3d(:,:,:)
      INTEGER, OPTIONAL, INTENT(IN) ::  add_array_i_2d(:,:), add_array_i_3d(:,:,:)
      REAL(wp), INTENT(IN) :: ref_out_array_r_2d(:,:), ref_out_array_r_3d(:,:,:)
      INTEGER, INTENT(IN) ::  ref_out_array_i_2d(:,:), ref_out_array_i_3d(:,:,:)
      LOGICAL, INTENT(IN) ::  ref_out_array_l_2d(:,:), ref_out_array_l_3d(:,:,:)
      TYPE(t_comm_pattern), INTENT(INOUT) :: comm_pattern

      CALL exchange_data(p_pat=comm_pattern, recv=out_array_r_2d, &
        &                send=in_array_r_2d, add=add_array_r_2d, &
        &                l_recv_exists=.TRUE.)
      IF (ANY(out_array_r_2d /= ref_out_array_r_2d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result r_2d"
        CALL finish(method_name, message_text)
      END IF
      CALL exchange_data(p_pat=comm_pattern, recv=out_array_r_3d, &
        &                send=in_array_r_3d, add=add_array_r_3d)
      IF (ANY(out_array_r_3d /= ref_out_array_r_3d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result r_3d"
        CALL finish(method_name, message_text)
      END IF

      CALL exchange_data(p_pat=comm_pattern, recv=out_array_i_2d, &
        &                send=in_array_i_2d, add=add_array_i_2d, &
        &                l_recv_exists=.TRUE.)
      IF (ANY(out_array_i_2d /= ref_out_array_i_2d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result i_2d"
        CALL finish(method_name, message_text)
      END IF
      CALL exchange_data(p_pat=comm_pattern, recv=out_array_i_3d, &
        &                send=in_array_i_3d, add=add_array_i_3d)
      IF (ANY(out_array_i_3d /= ref_out_array_i_3d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result i_3d"
        CALL finish(method_name, message_text)
      END IF

      CALL exchange_data(p_pat=comm_pattern, recv=out_array_l_2d, &
        &                send=in_array_l_2d, l_recv_exists=.TRUE.)
      IF (ANY(out_array_l_2d /= ref_out_array_l_2d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result l_2d"
        CALL finish(method_name, message_text)
      END IF
      CALL exchange_data(p_pat=comm_pattern, recv=out_array_l_3d, &
        &                send=in_array_l_3d)
      IF (ANY(out_array_l_3d /= ref_out_array_l_3d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result l_3d"
        CALL finish(method_name, message_text)
      END IF
    END SUBROUTINE check_exchange

    SUBROUTINE check_exchange_4de1(in_array, out_array, ref_out_array, &
      &                            comm_pattern)

      REAL(wp), OPTIONAL, INTENT(IN) :: in_array(:,:,:,:)
      REAL(wp), INTENT(INOUT) :: out_array(:,:,:,:)
      REAL(wp), INTENT(IN) :: ref_out_array(:,:,:,:)
      TYPE(t_comm_pattern), INTENT(INOUT) :: comm_pattern

      INTEGER :: nfields, ndim2tot

      nfields = SIZE(out_array, 1)
      ndim2tot = nfields * SIZE(out_array, 3)

      CALL exchange_data_4de1(comm_pattern, nfields, ndim2tot, out_array, &
        &                     in_array)
      IF (ANY(out_array /= ref_out_array)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result 4de1"
        CALL finish(method_name, message_text)
      END IF

    END SUBROUTINE check_exchange_4de1

    SUBROUTINE check_exchange_mult( &
      & out_array1, in_array1, ref_out_array1, &
      & out_array2, in_array2, ref_out_array2, &
      & out_array3, in_array3, ref_out_array3, &
      & out_array4, in_array4, ref_out_array4, &
      & out_array5, in_array5, ref_out_array5, &
      & out_array6, in_array6, ref_out_array6, &
      & out_array7, in_array7, ref_out_array7, &
      & out_array4d, in_array4d, ref_out_array4d, &
      & nshift, comm_pattern)

      REAL(wp), INTENT(INOUT), OPTIONAL ::  &
        out_array1(:,:,:), out_array2(:,:,:), out_array3(:,:,:), &
        out_array4(:,:,:), out_array5(:,:,:), out_array6(:,:,:), &
        out_array7(:,:,:)
      REAL(wp), INTENT(IN), OPTIONAL ::  &
        in_array1(:,:,:), in_array2(:,:,:), in_array3(:,:,:), &
        in_array4(:,:,:), in_array5(:,:,:), in_array6(:,:,:), &
        in_array7(:,:,:)
      REAL(wp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1(:,:,:), ref_out_array2(:,:,:), ref_out_array3(:,:,:), &
        ref_out_array4(:,:,:), ref_out_array5(:,:,:), ref_out_array6(:,:,:), &
        ref_out_array7(:,:,:)
      REAL(wp), INTENT(INOUT), OPTIONAL :: out_array4d(:,:,:,:)
      REAL(wp), INTENT(IN   ), OPTIONAL :: in_array4d(:,:,:,:)
      REAL(wp), INTENT(INOUT), OPTIONAL :: ref_out_array4d(:,:,:,:)
      INTEGER, OPTIONAL, INTENT(IN) :: nshift
      TYPE(t_comm_pattern), INTENT(INOUT) :: comm_pattern

      INTEGER :: nfields, ndim2tot, ndim2, kshift

      kshift = 0
      IF (PRESENT(nshift)) kshift = nshift

      nfields = 0
      ndim2tot = 0
      IF (PRESENT(out_array4d)) THEN
        nfields = nfields + SIZE(out_array4d, 4)
        ndim2 = SIZE(out_array4d, 2)
        ndim2tot = ndim2tot + &
          &        MERGE(1, ndim2 - kshift, ndim2 == 1) * SIZE(out_array4d, 4)
      END IF
      IF (PRESENT(out_array1)) THEN
        nfields = nfields + 1
        ndim2 = SIZE(out_array1, 2)
        ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
        IF (PRESENT(out_array2)) THEN
          nfields = nfields + 1
          ndim2 = SIZE(out_array2, 2)
          ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
          IF (PRESENT(out_array3)) THEN
            nfields = nfields + 1
            ndim2 = SIZE(out_array3, 2)
            ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
            IF (PRESENT(out_array4)) THEN
              nfields = nfields + 1
              ndim2 = SIZE(out_array4, 2)
              ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
              IF (PRESENT(out_array5)) THEN
                nfields = nfields + 1
                ndim2 = SIZE(out_array5, 2)
                ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
                IF (PRESENT(out_array6)) THEN
                  nfields = nfields + 1
                  ndim2 = SIZE(out_array6, 2)
                  ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
                  IF (PRESENT(out_array7)) THEN
                    nfields = nfields + 1
                    ndim2 = SIZE(out_array7, 2)
                    ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF

      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=ndim2tot, &
        recv1=out_array1, send1=in_array1, &
        recv2=out_array2, send2=in_array2, &
        recv3=out_array3, send3=in_array3, &
        recv4=out_array4, send4=in_array4, &
        recv5=out_array5, send5=in_array5, &
        recv6=out_array6, send6=in_array6, &
        recv7=out_array7, send7=in_array7, &
        recv4d=out_array4d, send4d=in_array4d, &
        nshift=kshift)

      IF (PRESENT(out_array1)) THEN
        IF (ANY(out_array1 /= ref_out_array1)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result r_3d_1"
          CALL finish(method_name, message_text)
        END IF
      END IF
      IF (PRESENT(out_array2)) THEN
        IF (ANY(out_array2 /= ref_out_array2)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result r_3d_2"
          CALL finish(method_name, message_text)
        END IF
      END IF
      IF (PRESENT(out_array3)) THEN
        IF (ANY(out_array3 /= ref_out_array3)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result r_3d_3"
          CALL finish(method_name, message_text)
        END IF
      END IF
      IF (PRESENT(out_array4)) THEN
        IF (ANY(out_array4 /= ref_out_array4)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result r_3d_4"
          CALL finish(method_name, message_text)
        END IF
      END IF
      IF (PRESENT(out_array5)) THEN
        IF (ANY(out_array5 /= ref_out_array5)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result r_3d_5"
          CALL finish(method_name, message_text)
        END IF
      END IF
      IF (PRESENT(out_array6)) THEN
        IF (ANY(out_array6 /= ref_out_array6)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result r_3d_6"
          CALL finish(method_name, message_text)
        END IF
      END IF
      IF (PRESENT(out_array7)) THEN
        IF (ANY(out_array7 /= ref_out_array7)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result r_3d_7"
          CALL finish(method_name, message_text)
        END IF
      END IF
      IF (PRESENT(out_array4d)) THEN
        IF (ANY(out_array4d /= ref_out_array4d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result r_4d"
          CALL finish(method_name, message_text)
        END IF
      END IF
    END SUBROUTINE check_exchange_mult

  END SUBROUTINE exchange_communication_testbed

  SUBROUTINE exchange_communication_grf_testbed()

    CHARACTER(*), PARAMETER :: method_name = &
      "mo_test_communication:exchange_communication_grf_testbed"

    INTEGER, ALLOCATABLE :: owner_local_src(:), owner_local_dst(:,:), &
      &                     glb_index_src(:), glb_index_dst(:)
    REAL(wp), ALLOCATABLE :: ref_recv_add(:,:)
    INTEGER :: i, j, n, local_size_src, local_size_dst, global_size, count
    TYPE(t_glb2loc_index_lookup) :: send_glb2loc_index
    TYPE(t_comm_pattern) :: comm_pattern(4)
    TYPE(t_comm_pattern_collection) :: comm_pattern_collection

    REAL(wp) :: recv1(nproma,2,18), recv2(nproma,4,18), &
      &         recv3(nproma,6,18), recv4(nproma,8,18), &
      &         recv5(nproma,10,18), recv6(nproma,12,18), &
      &         recv4d1(nproma,4,18,6), recv4d2(nproma,8,18,6)
    REAL(wp) :: send1(2,12*nproma,4), send2(4,12*nproma,4), &
      &         send3(6,12*nproma,4), send4(8,12*nproma,4), &
      &         send5(10,12*nproma,4), send6(12,12*nproma,4), &
      &         send4d1(4,12*nproma,4,6), send4d2(8,12*nproma,4,6)
    REAL(wp) :: ref_recv1(nproma,2,18), ref_recv2(nproma,4,18), &
      &         ref_recv3(nproma,6,18), ref_recv4(nproma,8,18), &
      &         ref_recv5(nproma,10,18), ref_recv6(nproma,12,18), &
      &         ref_recv4d1(nproma,4,18,6), ref_recv4d2(nproma,8,18,6)

    !generate communication patterns
    local_size_src = 12 * nproma
    local_size_dst = 18 * nproma
    global_size = local_size_src * p_n_work
    ALLOCATE(owner_local_src(local_size_src), glb_index_src(local_size_src), &
      &      owner_local_dst(local_size_dst,4), glb_index_dst(local_size_dst), &
      &      ref_recv_add(nproma,18))

    owner_local_src = p_pe_work
    owner_local_dst = -1
    DO i = 1, 4
      owner_local_dst(i:3*nproma:4,i) = MOD(p_pe_work + p_n_work - 1, p_n_work)
      owner_local_dst(3*nproma+i:15*nproma:4,i) = p_pe_work
      owner_local_dst(15*nproma+i:18*nproma:4,i) = MOD(p_pe_work + 1, p_n_work)
    END DO
    DO i = 1, local_size_src
      glb_index_src(i) = local_size_src * p_pe_work + i
    END DO
    DO i = 1, local_size_dst
      glb_index_dst(i) = MOD(local_size_src * p_pe_work - 3 * nproma - 1 + i + &
      &                      global_size, global_size) + 1
    END DO
    CALL init_glb2loc_index_lookup(send_glb2loc_index, global_size)
    CALL set_inner_glb_index(send_glb2loc_index, glb_index_src, &
      &                      (/(i, i=1, local_size_src)/))

    DO i = 1, 4
      CALL setup_comm_pattern(local_size_dst, owner_local_dst(:,i), &
        &                     glb_index_dst, send_glb2loc_index, &
        &                     local_size_src, owner_local_src, glb_index_src, &
        &                     comm_pattern(i))
    END DO

    CALL setup_comm_pattern_collection(comm_pattern, comm_pattern_collection)

    recv1 = -1
    recv2 = -1
    recv3 = -1
    recv4 = -1
    recv5 = -1
    recv6 = -1
    recv4d1 = -1
    recv4d2 = -1

    count = 0
    DO n = 1, SIZE(send1, 1)
      DO i = 1, local_size_src
        send1(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(send2, 1)
      DO i = 1, local_size_src
        send2(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(send3, 1)
      DO i = 1, local_size_src
        send3(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(send4, 1)
      DO i = 1, local_size_src
        send4(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(send5, 1)
      DO i = 1, local_size_src
        send5(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(send6, 1)
      DO i = 1, local_size_src
        send6(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO

    count = 0
    DO j = 1, SIZE(send4d1, 4)
      DO n = 1, SIZE(send4d1, 1)
        DO i = 1, local_size_src
          send4d1(n,i,:, j) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
            &                 p_pe_work * local_size_src
        END DO
        count = count + local_size_src
      END DO
    END DO
    DO j = 1, SIZE(send4d2, 4)
      DO n = 1, SIZE(send4d2, 1)
        DO i = 1, local_size_src
          send4d2(n,i,:, j) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &                   p_pe_work * local_size_src
        END DO
        count = count + local_size_src
      END DO
    END DO

    ref_recv_add(1::4,:) = 0.1_wp
    ref_recv_add(2::4,:) = 0.2_wp
    ref_recv_add(3::4,:) = 0.3_wp
    ref_recv_add(4::4,:) = 0.4_wp

    count = 0
    DO n = 1, SIZE(ref_recv1, 2)
      ref_recv1(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(ref_recv2, 2)
      ref_recv2(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(ref_recv3, 2)
      ref_recv3(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(ref_recv4, 2)
      ref_recv4(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(ref_recv5, 2)
      ref_recv5(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(ref_recv6, 2)
      ref_recv6(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO

    count = 0
    DO i = 1, SIZE(ref_recv4d1, 4)
      DO n = 1, SIZE(ref_recv4d1, 2)
        ref_recv4d1(:,n,:,i) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
        count = count + local_size_src
      END DO
    END DO
    DO i = 1, SIZE(ref_recv4d2, 4)
      DO n = 1, SIZE(ref_recv4d2, 2)
        ref_recv4d2(:,n,:,i) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
        count = count + local_size_src
      END DO
    END DO

    CALL check_exchange_grf(comm_pattern_collection, recv1, send1, ref_recv1, &
      &                     recv2, send2, ref_recv2, recv3, send3, ref_recv3, &
      &                     recv4, send4, ref_recv4, recv5, send5, ref_recv5, &
      &                     recv6, send6, ref_recv6, recv4d1, send4d1, &
      &                     ref_recv4d1, recv4d2, send4d2, ref_recv4d2)

    CALL delete_comm_pattern_collection(comm_pattern_collection)

    CONTAINS

    SUBROUTINE check_exchange_grf( &
      p_pat_coll, recv1, send1, ref_recv1, recv2, send2, ref_recv2, recv3, &
      send3, ref_recv3, recv4, send4, ref_recv4, recv5, send5, ref_recv5, &
      recv6, send6, ref_recv6, recv4d1, send4d1, ref_recv4d1, recv4d2, &
      send4d2, ref_recv4d2)

      TYPE(t_comm_pattern_collection), INTENT(INOUT), TARGET :: p_pat_coll

      ! recv3d (nproma,nlev,blk)
      ! recv4d (nproma,nlev,blk,nfield)
      REAL(wp), INTENT(INOUT) :: recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), &
        &                        recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), &
        &                        recv4d1(:,:,:,:), recv4d2(:,:,:,:)
      ! send3d (nlev,i,npat)
      ! send4d (nlev,i,npat,nfield)
      REAL(wp), INTENT(IN) ::  send1(:,:,:), send2(:,:,:), send3(:,:,:), &
        &                      send4(:,:,:), send5(:,:,:), send6(:,:,:), &
        &                      send4d1(:,:,:,:), send4d2(:,:,:,:)

      REAL(wp), INTENT(IN) :: ref_recv1(:,:,:), ref_recv2(:,:,:), &
        &                     ref_recv3(:,:,:), ref_recv4(:,:,:), &
        &                     ref_recv5(:,:,:), ref_recv6(:,:,:), &
        &                     ref_recv4d1(:,:,:,:), ref_recv4d2(:,:,:,:)

      REAL(wp), ALLOCATABLE :: tmp_recv1(:,:,:), tmp_recv2(:,:,:), &
        &                      tmp_recv3(:,:,:), tmp_recv4(:,:,:), &
        &                      tmp_recv5(:,:,:), tmp_recv6(:,:,:), &
        &                      tmp_recv4d1(:,:,:,:), tmp_recv4d2(:,:,:,:)

      INTEGER :: nfields, ndim2tot, i

      ALLOCATE(tmp_recv1(SIZE(recv1,1),SIZE(recv1,2),SIZE(recv1,3)), &
        &      tmp_recv2(SIZE(recv2,1),SIZE(recv2,2),SIZE(recv2,3)), &
        &      tmp_recv3(SIZE(recv3,1),SIZE(recv3,2),SIZE(recv3,3)), &
        &      tmp_recv4(SIZE(recv4,1),SIZE(recv4,2),SIZE(recv4,3)), &
        &      tmp_recv5(SIZE(recv5,1),SIZE(recv5,2),SIZE(recv5,3)), &
        &      tmp_recv6(SIZE(recv6,1),SIZE(recv6,2),SIZE(recv6,3)), &
        &      tmp_recv4d1(SIZE(recv4d1,1),SIZE(recv4d1,2),&
        &                  SIZE(recv4d1,3),SIZE(recv4d1,4)), &
        &      tmp_recv4d2(SIZE(recv4d2,1),SIZE(recv4d2,2),&
        &                  SIZE(recv4d2,3),SIZE(recv4d2,4)))

      tmp_recv1 = recv1
      tmp_recv2 = recv2
      tmp_recv3 = recv3
      tmp_recv4 = recv4
      tmp_recv5 = recv5
      tmp_recv6 = recv6
      tmp_recv4d1 = recv4d1
      tmp_recv4d2 = recv4d2

      i = 1
      nfields = 1
      ndim2tot = SIZE(recv1, 2)
      recv1 = tmp_recv1
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1)
      IF (ANY(recv1 /= ref_recv1)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 2
      ndim2tot = SIZE(recv1, 2) + SIZE(recv2, 2)
      recv1 = tmp_recv1
      recv2 = tmp_recv2
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1, &
        &                    recv2=recv2, send2=send2)
      IF (ANY(recv1 /= ref_recv1) .OR. ANY(recv2 /= ref_recv2)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 3
      ndim2tot = SIZE(recv1, 2) + SIZE(recv2, 2) + SIZE(recv3, 2)
      recv1 = tmp_recv1
      recv2 = tmp_recv2
      recv3 = tmp_recv3
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1, &
        &                    recv2=recv2, send2=send2, recv3=recv3, &
        &                    send3=send3)
      IF (ANY(recv1 /= ref_recv1) .OR. ANY(recv2 /= ref_recv2) .OR. &
        & ANY(recv3 /= ref_recv3)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 4
      ndim2tot = SIZE(recv1, 2) + SIZE(recv2, 2) + SIZE(recv3, 2) + &
        &        SIZE(recv4, 2)
      recv1 = tmp_recv1
      recv2 = tmp_recv2
      recv3 = tmp_recv3
      recv4 = tmp_recv4
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1, &
        &                    recv2=recv2, send2=send2, recv3=recv3, &
        &                    send3=send3, recv4=recv4, send4=send4)
      IF (ANY(recv1 /= ref_recv1) .OR. ANY(recv2 /= ref_recv2) .OR. &
        & ANY(recv3 /= ref_recv3) .OR. ANY(recv4 /= ref_recv4)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 5
      ndim2tot = SIZE(recv1, 2) + SIZE(recv2, 2) + SIZE(recv3, 2) + &
        &        SIZE(recv4, 2) + SIZE(recv5, 2)
      recv1 = tmp_recv1
      recv2 = tmp_recv2
      recv3 = tmp_recv3
      recv4 = tmp_recv4
      recv5 = tmp_recv5
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1, &
        &                    recv2=recv2, send2=send2, recv3=recv3, &
        &                    send3=send3, recv4=recv4, send4=send4, &
        &                    recv5=recv5, send5=send5)
      IF (ANY(recv1 /= ref_recv1) .OR. ANY(recv2 /= ref_recv2) .OR. &
        & ANY(recv3 /= ref_recv3) .OR. ANY(recv4 /= ref_recv4) .OR. &
        & ANY(recv5 /= ref_recv5)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 6
      ndim2tot = SIZE(recv1, 2) + SIZE(recv2, 2) + SIZE(recv3, 2) + &
        &        SIZE(recv4, 2) + SIZE(recv5, 2) + SIZE(recv6, 2)
      recv1 = tmp_recv1
      recv2 = tmp_recv2
      recv3 = tmp_recv3
      recv4 = tmp_recv4
      recv5 = tmp_recv5
      recv6 = tmp_recv6
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1, &
        &                    recv2=recv2, send2=send2, recv3=recv3, &
        &                    send3=send3, recv4=recv4, send4=send4, &
        &                    recv5=recv5, send5=send5, recv6=recv6, &
        &                    send6=send6)
      IF (ANY(recv1 /= ref_recv1) .OR. ANY(recv2 /= ref_recv2) .OR. &
        & ANY(recv3 /= ref_recv3) .OR. ANY(recv4 /= ref_recv4) .OR. &
        & ANY(recv5 /= ref_recv5) .OR. ANY(recv6 /= ref_recv6)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = SIZE(recv4d1, 4)
      ndim2tot = SIZE(recv4d1, 4) * SIZE(recv4d1, 2)
      recv4d1 = tmp_recv4d1
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv4d1=recv4d1, &
        &                    send4d1=send4d1)
      IF (ANY(recv4d1 /= ref_recv4d1)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 2 * SIZE(recv4d1, 4)
      ndim2tot = SIZE(recv4d1, 4) * SIZE(recv4d1, 2) + &
        &        SIZE(recv4d2, 4) * SIZE(recv4d2, 2)
      recv4d1 = tmp_recv4d1
      recv4d2 = tmp_recv4d2
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv4d1=recv4d1, &
        &                    send4d1=send4d1, recv4d2=recv4d2, send4d2=send4d2)
      IF (ANY(recv4d1 /= ref_recv4d1) .OR. ANY(recv4d2 /= ref_recv4d2)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF
    END SUBROUTINE check_exchange_grf

  END SUBROUTINE exchange_communication_grf_testbed

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
    INTEGER :: intercomm, intercomm_key, p_comm_work_new, ierr
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
    CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
      &                        in_array_i_1d, in_array_i_2d, &
      &                        out_array_r_1d, out_array_r_2d, &
      &                        ref_out_array_r_1d, ref_out_array_r_2d, &
      &                        out_array_i_1d, out_array_i_2d, &
      &                        ref_out_array_i_1d, ref_out_array_i_2d, &
      &                        gather_pattern, __LINE__)
    CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
      &                        in_array_i_1d, in_array_i_2d, &
      &                        out_array_r_1d, out_array_r_2d, &
      &                        ref_out_array_r_1d, ref_out_array_r_2d, &
      &                        out_array_i_1d, out_array_i_2d, &
      &                        ref_out_array_i_1d, ref_out_array_i_2d, &
      &                        gather_pattern, __LINE__, fill_value)

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
    IF (p_n_work > 1) THEN
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
      CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
        &                        in_array_i_1d, in_array_i_2d, &
        &                        out_array_r_1d, out_array_r_2d, &
        &                        ref_out_array_r_1d, ref_out_array_r_2d, &
        &                        out_array_i_1d, out_array_i_2d, &
        &                        ref_out_array_i_1d, ref_out_array_i_2d, &
        &                        gather_pattern, __LINE__)
      CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
        &                        in_array_i_1d, in_array_i_2d, &
        &                        out_array_r_1d, out_array_r_2d, &
        &                        ref_out_array_r_1d, ref_out_array_r_2d, &
        &                        out_array_i_1d, out_array_i_2d, &
        &                        ref_out_array_i_1d, ref_out_array_i_2d, &
        &                        gather_pattern, __LINE__, fill_value)

      ! delete gather pattern and other arrays
      DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
      DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
        &        ref_out_array_r_1d, ref_out_array_r_2d, &
        &        ref_out_array_i_1d, ref_out_array_i_2d)
      CALL delete_comm_gather_pattern(gather_pattern)
      DEALLOCATE(owner_local, glb_index)
    END IF ! (p_n_work > 1)

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
    CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
      &                        in_array_i_1d, in_array_i_2d, &
      &                        out_array_r_1d, out_array_r_2d, &
      &                        ref_out_array_r_1d, ref_out_array_r_2d, &
      &                        out_array_i_1d, out_array_i_2d, &
      &                        ref_out_array_i_1d, ref_out_array_i_2d, &
      &                        gather_pattern, __LINE__, fill_value)

    ! delete gather pattern and other arrays
    DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
    DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
      &        ref_out_array_r_1d, ref_out_array_r_2d, &
      &        ref_out_array_i_1d, ref_out_array_i_2d)
    CALL delete_comm_gather_pattern(gather_pattern)
    DEALLOCATE(owner_local, glb_index)
    !---------------------------------------------------------------------------
    ! test in which only even numbered global indices are owned by processes
    !---------------------------------------------------------------------------
    ! generate gather pattern
    local_size = 10 * nproma
    global_size = p_n_work * local_size
    ALLOCATE(owner_local(local_size), glb_index(local_size))
    DO i = 2, local_size, 2
      owner_local(i) = p_pe_work
    END DO
    DO i = 1, local_size, 2
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
      DO i = 1, global_size - 1, 2
        ref_out_array_r_1d(i+1) = i
        ref_out_array_i_1d(i+1) = i
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = (j - 1) * global_size + i
          ref_out_array_i_2d(i+1, j) = (j - 1) * global_size + i
        END DO
      END DO
      DO i = 0, global_size - 1, 2
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
    CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
      &                        in_array_i_1d, in_array_i_2d, &
      &                        out_array_r_1d, out_array_r_2d, &
      &                        ref_out_array_r_1d, ref_out_array_r_2d, &
      &                        out_array_i_1d, out_array_i_2d, &
      &                        ref_out_array_i_1d, ref_out_array_i_2d, &
      &                        gather_pattern, __LINE__, fill_value)

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
    IF (p_n_work > 1) THEN
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
      CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
        &                        in_array_i_1d, in_array_i_2d, &
        &                        out_array_r_1d, out_array_r_2d, &
        &                        ref_out_array_r_1d, ref_out_array_r_2d, &
        &                        out_array_i_1d, out_array_i_2d, &
        &                        ref_out_array_i_1d, ref_out_array_i_2d, &
        &                        gather_pattern, __LINE__, fill_value)

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
      CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
        &                        in_array_i_1d, in_array_i_2d, &
        &                        out_array_r_1d, out_array_r_2d, &
        &                        ref_out_array_r_1d, ref_out_array_r_2d, &
        &                        out_array_i_1d, out_array_i_2d, &
        &                        ref_out_array_i_1d, ref_out_array_i_2d, &
        &                        gather_pattern, __LINE__)

      ! delete gather pattern and other arrays
      DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
      DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
        &        ref_out_array_r_1d, ref_out_array_r_2d, &
        &        ref_out_array_i_1d, ref_out_array_i_2d)
      CALL delete_comm_gather_pattern(gather_pattern)
      DEALLOCATE(owner_local, glb_index)
    END IF ! (p_n_work > 1)

#ifndef NOMPI
    IF (p_n_work > 1) THEN
      !---------------------------------------------------------------------------
      ! initial setup for allgather_intercomm tests
      !---------------------------------------------------------------------------
      p_comm_work_backup = p_comm_work
      p_pe_work_backup = p_pe_work
      p_n_work_backup = p_n_work
      intercomm_key = MERGE(0, 1, p_pe_work < (p_n_work * 2)/3)
      CALL MPI_Comm_split(p_comm_work, intercomm_key, p_pe_work, &
        &                 p_comm_work_new, ierr)
      CALL MPI_Intercomm_create(p_comm_work_new, 0, p_comm_work, &
        MERGE((p_n_work*2)/3, 0, intercomm_key == 0), 2, intercomm, ierr)
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

        is_active = (intercomm_key == 0) .EQV. (j == 0)

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
    END IF ! (p_n_work > 1)
#endif

  CONTAINS

    SUBROUTINE check_exchange_gather(in_array_r_1d, in_array_r_2d, &
      &                              in_array_i_1d, in_array_i_2d, &
      &                              out_array_r_1d, out_array_r_2d, &
      &                              ref_out_array_r_1d, ref_out_array_r_2d, &
      &                              out_array_i_1d, out_array_i_2d, &
      &                              ref_out_array_i_1d, ref_out_array_i_2d, &
      &                              gather_pattern, line, fill_value)

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

