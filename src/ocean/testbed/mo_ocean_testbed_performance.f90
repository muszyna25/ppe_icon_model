!>
!! @brief Testbed for ocean
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_testbed_ocean_performance

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_timer,               ONLY: init_timer, ltimer, new_timer, timer_start, timer_stop, &
    & print_timer, activate_sync_timers, timers_level, timer_barrier, timer_radiaton_recv
  USE mo_parallel_config,     ONLY: nproma, icon_comm_method
  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no

  ! USE mo_icon_testbed_config, ONLY: testbed_iterations

  USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_impl_constants,            ONLY: sea_boundary, sea
  USE mo_math_constants,            ONLY: pi
  USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer,                                                  &
    &                                     threshold_min_T, threshold_max_T, threshold_min_S, threshold_max_S, &
    &                                     type_3dimRelax_Temp, para_3dimRelax_Temp,                           &
    &                                     type_3dimRelax_Salt, para_3dimRelax_Salt,                           &
    &                                     iswm_oce, l_edge_based,                                             &
    &                                     FLUX_CALCULATION_HORZ, FLUX_CALCULATION_VERT,                       &
    &                                     forcing_enable_freshwater
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer, test_mode
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state
  USE mo_model_domain,              ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_ext_data_types,            ONLY: t_external_data
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state
  USE mo_ocean_physics_types,       ONLY: t_ho_params
  USE mo_sea_ice_types,             ONLY: t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean, t_sea_ice
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff

  USE mo_grid_config,         ONLY: n_dom

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: ocean_test_performance

INTEGER :: testbed_iterations = 100

CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE ocean_test_performance( namelist_filename, shr_namelist_filename, &
    & patch_3d, ocean_state, external_data,                     &
    & surface_fluxes, physics_parameters,             &
    & oceans_atmosphere, oceans_atmosphere_fluxes, ocean_ice, operators_coefficients)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_external_data), TARGET, INTENT(in)        :: external_data(n_dom)
    TYPE(t_sfc_flx)                                  :: surface_fluxes
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: oceans_atmosphere
    TYPE(t_atmos_fluxes ),    INTENT(inout)          :: oceans_atmosphere_fluxes
    TYPE (t_sea_ice),         INTENT(inout)          :: ocean_ice
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    
    CHARACTER(*), PARAMETER :: method_name = "ocean_test_performance"

    SELECT CASE (test_mode)  !  1000 - 1100
      CASE (1000)
        CALL test_prepare_tracer_transport( patch_3d, ocean_state, operators_coefficients)

      CASE DEFAULT
        CALL finish(method_name, "Unknown test_mode")

    END SELECT



  END SUBROUTINE ocean_test_performance

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_prepare_tracer_transport( patch_3d, ocean_state, operators_coefficients)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients

    INTEGER :: timer_prep_trace_trans_0, timer_prep_trace_trans_1, timer_prep_trace_trans_2, &
      & timer_prep_trace_trans_3
    INTEGER :: iter
    CHARACTER(*), PARAMETER :: method_name = "mo_testbed_ocean_performance:ocean_test_performance"

    !---------------------------------------------------------------------
    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name
    
    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    !---------------------------------------------------------------------
        
    !---------------------------------------------------------------------
    ! DO the tests
    timer_prep_trace_trans_0  = new_timer("prep_trace_trans_0")
    timer_prep_trace_trans_1  = new_timer("prep_trace_trans_1")
    timer_prep_trace_trans_2  = new_timer("prep_trace_trans_2")
    timer_prep_trace_trans_3  = new_timer("prep_trace_trans_3")

    !---------------------------------------------------------------------
    ! measure the original implementation
    CALL work_mpi_barrier()
    DO iter=1, testbed_iterations
      CALL timer_start(timer_prep_trace_trans_0)
      CALL prepare_tracer_transport_0( patch_3D, ocean_state(1), operators_coefficients)
      CALL timer_stop(timer_prep_trace_trans_0)
    ENDDO
    !---------------------------------------------------------------------
    ! measure 1st implementation
    CALL work_mpi_barrier()
    DO iter=1, testbed_iterations
      CALL timer_start(timer_prep_trace_trans_1)
      CALL prepare_tracer_transport_1( patch_3D, ocean_state(1), operators_coefficients)
      CALL timer_stop(timer_prep_trace_trans_1)
    ENDDO
    !---------------------------------------------------------------------
    ! measure 2nd implementation
    CALL work_mpi_barrier()
    DO iter=1, testbed_iterations
      CALL timer_start(timer_prep_trace_trans_2)
      CALL prepare_tracer_transport_2( patch_3D, ocean_state(1), operators_coefficients)
      CALL timer_stop(timer_prep_trace_trans_2)
    ENDDO
    !---------------------------------------------------------------------
    ! measure 3d implementation
    CALL work_mpi_barrier()
    DO iter=1, testbed_iterations
      CALL timer_start(timer_prep_trace_trans_3)
      CALL prepare_tracer_transport_3( patch_3D, ocean_state(1), operators_coefficients)
      CALL timer_stop(timer_prep_trace_trans_3)
    ENDDO
    !---------------------------------------------------------------------


    !---------------------------------------------------------------------
    ! print the timers
    CALL work_mpi_barrier()
    CALL message("===================", "=======================")
    CALL message(method_name, TRIM(message_text))
    CALL print_timer()
    !---------------------------------------------------------------------
    

  END SUBROUTINE test_prepare_tracer_transport
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !!  0 version (original), without the syncs
  !!
  SUBROUTINE prepare_tracer_transport_0(patch_3D, p_os, p_op_coeff)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff),INTENT(INOUT) :: p_op_coeff
    !
    !Local variables
    INTEGER  :: slev, elev
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: i_startidx_e, i_endidx_e
    INTEGER  :: je, jk, jb,jc         !< index of edge, vert level, block
    INTEGER  :: il_c1, il_c2, ib_c1, ib_c2
    INTEGER  :: il_c, ib_c
    REAL(wp) :: delta_z
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc
    TYPE(t_cartesian_coordinates):: flux_sum
    !-------------------------------------------------------------------------------
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    !-------------------------------------------------------------------------------
    patch_2d         => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2d%cells%in_domain
    edges_in_domain => patch_2d%edges%in_domain

    slev = 1
    elev = n_zlev

    ! This should be moved to the vertical advection
    ! just moving data around should not take place
!    p_os%p_diag%w_time_weighted(1:nproma, 1:n_zlev+1, 1:patch_2d%nblks_c) = &
!      &  p_os%p_diag%w(1:nproma, 1:n_zlev+1, 1:patch_2d%nblks_c)

    DO jk = slev, elev
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e

          IF (patch_3D%lsm_e(je,jk,jb) == sea) THEN

            !Get indices of two adjacent cells
            il_c1 = patch_2d%edges%cell_idx(je,jb,1)
            ib_c1 = patch_2d%edges%cell_blk(je,jb,1)
            il_c2 = patch_2d%edges%cell_idx(je,jb,2)
            ib_c2 = patch_2d%edges%cell_blk(je,jb,2)

            p_os%p_diag%p_vn_mean(je,jk,jb)%x = 0.5_wp*&
              &(p_os%p_diag%p_vn(il_c1,jk,ib_c1)%x+p_os%p_diag%p_vn(il_c2,jk,ib_c2)%x)

            p_op_coeff%moved_edge_position_cc(je,jk,jb)%x = &
              & p_op_coeff%edge_position_cc(je,jk,jb)%x     &
              &  - 0.5_wp*dtime*p_os%p_diag%p_vn_mean(je,jk,jb)%x

            IF ( p_os%p_diag%vn_time_weighted(je,jk,jb) > 0.0_wp ) THEN
              il_c = patch_2d%edges%cell_idx(je,jb,1)
              ib_c = patch_2d%edges%cell_blk(je,jb,1)
            ELSE  ! p_os%p_diag%vn_time_weighted <= 0.0
              il_c = patch_2d%edges%cell_idx(je,jb,2)
              ib_c = patch_2d%edges%cell_blk(je,jb,2)
            ENDIF

            p_op_coeff%upwind_cell_idx(je,jk,jb) = il_c
            p_op_coeff%upwind_cell_blk(je,jk,jb) = ib_c

            p_op_coeff%upwind_cell_position_cc(je,jk,jb)%x = &
              & patch_2d%cells%cartesian_center(il_c,ib_c)%x

          ENDIF

        END DO
      END DO
    END DO

  END SUBROUTINE prepare_tracer_transport_0
  !-------------------------------------------------------------------------
    
  !-------------------------------------------------------------------------
  !>
  !!  1 version
  !!
  SUBROUTINE prepare_tracer_transport_1(patch_3D, p_os, p_op_coeff)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff),INTENT(INOUT) :: p_op_coeff
    !
    !Local variables
    INTEGER  :: slev, elev
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: i_startidx_e, i_endidx_e
    INTEGER  :: je, jk, jb,jc         !< index of edge, vert level, block
    INTEGER  :: il_c1, il_c2, ib_c1, ib_c2
    INTEGER  :: il_c, ib_c
    REAL(wp) :: delta_z
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc
    TYPE(t_cartesian_coordinates):: flux_sum
    !-------------------------------------------------------------------------------
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    !-------------------------------------------------------------------------------
    patch_2d         => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2d%cells%in_domain
    edges_in_domain => patch_2d%edges%in_domain

    slev = 1
    elev = n_zlev

    ! This should be moved to the vertical advection
    ! just moving data around should not take place
!    p_os%p_diag%w_time_weighted(1:nproma, 1:n_zlev+1, 1:patch_2d%nblks_c) = &
!      &  p_os%p_diag%w(1:nproma, 1:n_zlev+1, 1:patch_2d%nblks_c)

    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e
        !Get indices of two adjacent cells
        il_c1 = patch_2d%edges%cell_idx(je,jb,1)
        ib_c1 = patch_2d%edges%cell_blk(je,jb,1)
        il_c2 = patch_2d%edges%cell_idx(je,jb,2)
        ib_c2 = patch_2d%edges%cell_blk(je,jb,2)

        DO jk = slev, elev

          IF (patch_3D%lsm_e(je,jk,jb) == sea) THEN


            p_os%p_diag%p_vn_mean(je,jk,jb)%x = 0.5_wp*&
              &(p_os%p_diag%p_vn(il_c1,jk,ib_c1)%x+p_os%p_diag%p_vn(il_c2,jk,ib_c2)%x)

            p_op_coeff%moved_edge_position_cc(je,jk,jb)%x = &
              & p_op_coeff%edge_position_cc(je,jk,jb)%x     &
              &  - 0.5_wp*dtime*p_os%p_diag%p_vn_mean(je,jk,jb)%x

            IF ( p_os%p_diag%vn_time_weighted(je,jk,jb) > 0.0_wp ) THEN
              il_c = patch_2d%edges%cell_idx(je,jb,1)
              ib_c = patch_2d%edges%cell_blk(je,jb,1)
            ELSE  ! p_os%p_diag%vn_time_weighted <= 0.0
              il_c = patch_2d%edges%cell_idx(je,jb,2)
              ib_c = patch_2d%edges%cell_blk(je,jb,2)
            ENDIF

            p_op_coeff%upwind_cell_idx(je,jk,jb) = il_c
            p_op_coeff%upwind_cell_blk(je,jk,jb) = ib_c

            p_op_coeff%upwind_cell_position_cc(je,jk,jb)%x = &
              & patch_2d%cells%cartesian_center(il_c,ib_c)%x

          ENDIF

        END DO
      END DO
    END DO

  END SUBROUTINE prepare_tracer_transport_1
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!  2 version
  !!
  SUBROUTINE prepare_tracer_transport_2(patch_3D, p_os, p_op_coeff)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff),INTENT(INOUT) :: p_op_coeff
    !
    !Local variables
    INTEGER  :: slev, elev
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: i_startidx_e, i_endidx_e
    INTEGER  :: je, jk, jb,jc         !< index of edge, vert level, block
    INTEGER  :: edge_cell_index(2), edge_cell_block(2)
    INTEGER  :: upwind_index
    REAL(wp) :: delta_z, half_time
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc
    TYPE(t_cartesian_coordinates):: flux_sum
    !-------------------------------------------------------------------------------
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    !-------------------------------------------------------------------------------
    patch_2d         => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2d%cells%in_domain
    edges_in_domain => patch_2d%edges%in_domain

    slev = 1
    half_time = 0.5_wp * dtime

    ! This should be moved to the vertical advection
    ! just moving data around should not take place
!    p_os%p_diag%w_time_weighted(1:nproma, 1:n_zlev+1, 1:patch_2d%nblks_c) = &
!      &  p_os%p_diag%w(1:nproma, 1:n_zlev+1, 1:patch_2d%nblks_c)

    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e
        !Get indices of two adjacent cells
        edge_cell_index(1) = patch_2d%edges%cell_idx(je,jb,1)
        edge_cell_block(2) = patch_2d%edges%cell_blk(je,jb,1)
        edge_cell_index(2) = patch_2d%edges%cell_idx(je,jb,2)
        edge_cell_block(2) = patch_2d%edges%cell_blk(je,jb,2)
        elev  = patch_3D%p_patch_1D(1)%dolic_e(je,jb)

        DO jk = slev, elev

          upwind_index = MERGE(1, 2, p_os%p_diag%vn_time_weighted(je,jk,jb) > 0.0_wp)

          p_os%p_diag%p_vn_mean(je,jk,jb)%x = 0.5_wp *                          &
            & (p_os%p_diag%p_vn(edge_cell_index(1), jk, edge_cell_block(1))%x + &
            &  p_os%p_diag%p_vn(edge_cell_index(2), jk, edge_cell_block(2))%x)

          p_op_coeff%upwind_cell_idx(je,jk,jb) = edge_cell_index(upwind_index)
          p_op_coeff%upwind_cell_blk(je,jk,jb) = edge_cell_block(upwind_index)

          p_op_coeff%upwind_cell_position_cc(je,jk,jb)%x = &
            & patch_2d%cells%cartesian_center(edge_cell_index(upwind_index), edge_cell_block(upwind_index))%x
          ! & p_op_coeff%cell_position_cc(edge_cell_index(upwind_index), jk, edge_cell_block(upwind_index))%x

          p_op_coeff%moved_edge_position_cc(je,jk,jb)%x =   &
            &  p_op_coeff%edge_position_cc(je,jk,jb)%x      &
            &  - half_time * p_os%p_diag%p_vn_mean(je,jk,jb)%x

        END DO
      END DO
    END DO

  END SUBROUTINE prepare_tracer_transport_2
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!  3 version
  !!
  SUBROUTINE prepare_tracer_transport_3(patch_3D, p_os, p_op_coeff)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET    :: p_os
    TYPE(t_operator_coeff),INTENT(INOUT) :: p_op_coeff
    !
    !Local variables
    INTEGER  :: slev, elev
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: i_startidx_e, i_endidx_e
    INTEGER  :: je, jk, jb,jc         !< index of edge, vert level, block
    INTEGER  :: edge_cell_index(2), edge_cell_block(2)
    INTEGER  :: upwind_index
    REAL(wp) :: delta_z, half_time
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc
    TYPE(t_cartesian_coordinates):: flux_sum
    !-------------------------------------------------------------------------------
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
    !-------------------------------------------------------------------------------
    patch_2d        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2d%cells%in_domain
    edges_in_domain => patch_2d%edges%in_domain

    slev = 1
    half_time = 0.5_wp * dtime

    ! This should be moved to the vertical advection
    ! just moving data around should not take place
!    p_os%p_diag%w_time_weighted(1:nproma, 1:n_zlev+1, 1:patch_2d%nblks_c) = &
!      &  p_os%p_diag%w(1:nproma, 1:n_zlev+1, 1:patch_2d%nblks_c)

    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e
        !Get indices of two adjacent cells
        edge_cell_index(1) = patch_2d%edges%cell_idx(je,jb,1)
        edge_cell_block(1) = patch_2d%edges%cell_blk(je,jb,1)
        edge_cell_index(2) = patch_2d%edges%cell_idx(je,jb,2)
        edge_cell_block(2) = patch_2d%edges%cell_blk(je,jb,2)
        elev  = patch_3D%p_patch_1D(1)%dolic_e(je,jb)

        DO jk = slev, elev

          p_os%p_diag%p_vn_mean(je,jk,jb)%x = 0.5_wp *                          &
            & (p_os%p_diag%p_vn(edge_cell_index(1), jk, edge_cell_block(1))%x + &
            &  p_os%p_diag%p_vn(edge_cell_index(2), jk, edge_cell_block(2))%x)

          p_op_coeff%moved_edge_position_cc(je,jk,jb)%x =   &
            &  p_op_coeff%edge_position_cc(je,jk,jb)%x      &
            &  - half_time * p_os%p_diag%p_vn_mean(je,jk,jb)%x

        END DO

        DO jk = slev, elev
          upwind_index = MERGE(1, 2, p_os%p_diag%vn_time_weighted(je,jk,jb) > 0.0_wp)


          p_op_coeff%upwind_cell_idx(je,jk,jb) = edge_cell_index(upwind_index)
          p_op_coeff%upwind_cell_blk(je,jk,jb) = edge_cell_block(upwind_index)

          p_op_coeff%upwind_cell_position_cc(je,jk,jb)%x = &
            & patch_2d%cells%cartesian_center(edge_cell_index(upwind_index), edge_cell_block(upwind_index))%x
          ! & p_op_coeff%cell_position_cc(edge_cell_index(upwind_index), jk, edge_cell_block(upwind_index))%x
        END DO

      END DO
    END DO

  END SUBROUTINE prepare_tracer_transport_3
  !-------------------------------------------------------------------------


END MODULE mo_testbed_ocean_performance

