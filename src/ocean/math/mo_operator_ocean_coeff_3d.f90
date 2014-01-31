!>
!! Contains the definition of coefficients used for div-grad-curl and reconstruction/scalar product.
!! All coefficients are three-dimensional arrays including the number of vertical levels. This is necessary
!! if one has coefficients that vary within the vertical level but not in time such that one can precompute the
!! coefficients. This is in the ocean model where the land-sea mask is different at each level and therefore
!! the expansion coefficients associated with land and boundary vary with the vertical level but are constant in time.
!!
!!
!! @par Revision History
!! Developed  by Peter Korn (2012)
!!
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_operator_ocean_coeff_3d
  !-------------------------------------------------------------------------

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: success,&
    &                               max_char_length, beta_plane_coriolis,full_coriolis, &
    &                               SEA_BOUNDARY, BOUNDARY, SEA, min_dolic
  USE mo_math_constants,      ONLY: deg2rad, pi!, rad2deg
  USE mo_physical_constants,  ONLY: earth_radius
  USE mo_math_utilities,      ONLY: gc2cc, cc2gc, t_cartesian_coordinates,      &
    &                               t_geographical_coordinates, vector_product, &
    &                               arc_length
  USE mo_ocean_nml,           ONLY: n_zlev, no_tracer, &
    &                               coriolis_type, basin_center_lat, basin_height_deg
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array!, sync_idx, global_max
  !USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_oce_types,           ONLY: t_hydro_ocean_state, t_ptr3d, t_operator_coeff
  USE mo_oce_physics,         ONLY: t_ho_params
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_grid_config,         ONLY: grid_sphere_radius, grid_angular_velocity
  USE mo_run_config,          ONLY: dtime
  USE mo_var_list,            ONLY: add_var, add_ref
  USE mo_var_metadata,        ONLY: groups
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var
  USE mo_cdi_constants

  IMPLICIT NONE


  PRIVATE

  PUBLIC  :: t_operator_coeff
  PUBLIC  :: construct_operators_coefficients
  PUBLIC  :: destruct_operators_coefficients
  PUBLIC  :: update_diffusion_matrices


  !these two parameters are set below in sbr "allocate_operators_coefficients"
  !according to MAXVAL(patch_2D%cells%num_edges) and MAXVAL(patch_2D%verts%num_edges)
  INTEGER,PUBLIC :: no_dual_edges
  INTEGER,PUBLIC :: no_primal_edges 

  ! flags for computing ocean coefficients
  LOGICAL, PARAMETER :: MID_POINT_DUAL_EDGE = .TRUE. !Please do not change this unless you are sure, you know what you do.
  LOGICAL, PARAMETER :: LARC_LENGTH = .FALSE.

CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE construct_operators_coefficients( patch_3D, operators_coefficients, var_list)
    TYPE(t_patch_3D),TARGET,INTENT(inout) :: patch_3D
    TYPE(t_operator_coeff), INTENT(inout) :: operators_coefficients
    TYPE(t_var_list)                      :: var_list

    CALL allocate_operators_coefficients( patch_3d%p_patch_2d(1), operators_coefficients, var_list)
    CALL par_init_operator_coeff( patch_3d, operators_coefficients)

  END SUBROUTINE construct_operators_coefficients
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE destruct_operators_coefficients( operators_coefficients )
    TYPE(t_operator_coeff), INTENT(inout) :: operators_coefficients

    CALL deallocate_operators_coefficients( operators_coefficients )

  END SUBROUTINE destruct_operators_coefficients
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !> Allocation of operators coefficients.
  !
  ! @par Revision History
  ! Peter Korn (2012-2)
  !
  SUBROUTINE allocate_operators_coefficients( patch_2D, operators_coefficients, var_list)
    !
    TYPE(t_patch),TARGET,INTENT(in)       :: patch_2D
    TYPE(t_operator_coeff), INTENT(inout) :: operators_coefficients
    TYPE(t_var_list)                      :: var_list

    INTEGER :: nblks_c, nblks_e, nblks_v, nz_lev
    INTEGER :: ist,ie,i
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: jc,je,jk,jb
    CHARACTER(len=max_char_length) :: var_suffix

    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells

    !-----------------------------------------------------------------------
    !
    ! determine size of arrays, i.e.
    ! values for the blocking
    !
    nblks_c  = patch_2D%alloc_cell_blocks
    nblks_e  = patch_2D%nblks_e
    nblks_v  = patch_2D%nblks_v
    nz_lev   = n_zlev

    no_primal_edges = MAXVAL(patch_2D%cells%num_edges)
    no_dual_edges   = MAXVAL(patch_2D%verts%num_edges)

    ALLOCATE(operators_coefficients%div_coeff(nproma,n_zlev,nblks_c,no_primal_edges),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for geofac_div failed')
    ENDIF

    ALLOCATE(operators_coefficients%grad_coeff(nproma,n_zlev,nblks_e),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for geofac_grad failed')
    ENDIF

    ALLOCATE(operators_coefficients%rot_coeff(nproma,n_zlev,nblks_v,no_dual_edges),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d', &
        & 'allocation for geofac_rot failed')
    ENDIF

    ALLOCATE(operators_coefficients%n2s_coeff(nproma,n_zlev,nblks_c,no_primal_edges+1),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                       &
        & 'allocation for geofac_n2s failed')
    ENDIF
    ALLOCATE(operators_coefficients%n2v_coeff(nproma,n_zlev,nblks_e),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                       &
        & 'allocation for geofac_n2v failed')
    ENDIF
    !
    ALLOCATE(operators_coefficients%dist_cell2edge(nproma,n_zlev,nblks_e,2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating dist_cell2edge failed')
    ENDIF

    ALLOCATE(operators_coefficients%bnd_edge_idx(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating bnd_edge_idx failed')
    ENDIF
    ALLOCATE(operators_coefficients%bnd_edge_blk(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating bnd_edge_blk failed')
    ENDIF
    ALLOCATE(operators_coefficients%edge_idx(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge_idx failed')
    ENDIF
    ALLOCATE(operators_coefficients%orientation(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating orientation failed')
    ENDIF
    ALLOCATE(operators_coefficients%bnd_edges_per_vertex(nproma,n_zlev,nblks_v),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating bnd_edges_per_vertex failed')
    ENDIF
    ALLOCATE(operators_coefficients%upwind_cell_idx(nproma,n_zlev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind_cell_idx failed')
    ENDIF
    ALLOCATE(operators_coefficients%upwind_cell_blk(nproma,n_zlev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind_cell_blk failed')
    ENDIF
    !
    ! arrays that are required for setting up the scalar product
    !
    !coefficients for edge to cell mapping, one half of the scalar product.
    !Dimension: nproma,nblks_c encode number of cells, 1:3 corresponds to number
    !of edges per cell, 1:2 is for u and v component of cell vector
    !     ALLOCATE(operators_coefficients%edge2cell_coeff(nproma,nblks_c,1:3, 1:2),STAT=ist)
    !     IF (ist /= SUCCESS) THEN
    !       CALL finish ('allocating edge2cell_coeff failed')
    !     ENDIF
    ALLOCATE(operators_coefficients%edge2edge_viacell_coeff(nproma,nz_lev,nblks_e,1:2*no_primal_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2edge_viacell_coeff failed')
    ENDIF
    ALLOCATE(operators_coefficients%edge2edge_viacell_coeff_top(1:2*no_primal_edges, nproma, nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2edge_viacell_coeff_top failed')
    ENDIF
    ALLOCATE(operators_coefficients%edge2edge_viacell_coeff_integrated(1:2*no_primal_edges, nproma, nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2edge_viacell_coeff_integrated failed')
    ENDIF
    ALLOCATE(operators_coefficients%edge2edge_viacell_coeff_all(1:2*no_primal_edges, nproma, nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2edge_viacell_coeff_all failed')
    ENDIF

    ALLOCATE(operators_coefficients%edge2cell_coeff_cc(nproma,nz_lev,nblks_c,1:no_primal_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc failed')
    ENDIF

    ALLOCATE(operators_coefficients%edge2cell_coeff_cc_dyn(nproma,1,nblks_c,1:no_primal_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc_dyn failed')
    ENDIF
    !ALLOCATE(operators_coefficients%edge2vert_coeff_cc_dyn(nproma,1,nblks_v,1:no_dual_edges),stat=ist)
    !IF (ist /= success) THEN
    !  CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff_cc_dyn failed')
    !ENDIF

    !coefficients for transposed of edge to cell mapping, second half of the scalar product.
    !Dimension: nproma,nblks_e encode number of edges, 1:2 is for cell neighbors of an edge
    ALLOCATE(operators_coefficients%edge2cell_coeff_cc_t(nproma,nz_lev,nblks_e,1:2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating transposed edge2cell_coeff failed')
    ENDIF

    !
    !coefficients for edge to vertex mapping.
    !
    !Dimension: nproma,nblks_v encode number of vertices,
    !1:6 is number of edges of a vertex,
    !1:2 is for u and v component of vertex vector
    ALLOCATE(operators_coefficients%edge2vert_coeff_cc(nproma,nz_lev,nblks_v,1:no_dual_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff failed')
    ENDIF

    ALLOCATE(operators_coefficients%edge2vert_coeff_cc_t(nproma,nz_lev,nblks_e,1:2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff failed')
    ENDIF
    ALLOCATE(operators_coefficients%edge2vert_vector_cc(nproma,nz_lev,nblks_v,1:no_dual_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_vector failed')
    ENDIF

   ALLOCATE(operators_coefficients%edge2edge_viavert_coeff(nproma,nz_lev,nblks_e,1:2*no_dual_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc failed')
    ENDIF

    ALLOCATE(operators_coefficients%upwind_cell_position_cc(nproma,nz_lev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind cell failed')
    ENDIF
    ALLOCATE(operators_coefficients%moved_edge_position_cc(nproma,nz_lev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge failed')
    ENDIF
    ALLOCATE(operators_coefficients%edge_position_cc(nproma,nz_lev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge failed')
    ENDIF
!    ALLOCATE(operators_coefficients%cell_position_cc(nproma,nz_lev,nblks_c),stat=ist)
!    IF (ist /= success) THEN
!      CALL finish ('mo_operator_ocean_coeff_3d:allocating cell failed')
!    ENDIF
    !
    !normalizing factors for edge to cell mapping.
    !
    !Either by fixed volume or by variable one taking the surface elevation
    !into account. The later one depends on time and space.
    ALLOCATE(operators_coefficients%fixed_vol_norm(nproma,nz_lev,nblks_c),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating fixed_vol_norm failed')
    ENDIF
    ALLOCATE(operators_coefficients%variable_vol_norm(nproma,nz_lev,nblks_c,1:no_primal_edges),stat=ist)
    IF (ist /= success) THEN
    ENDIF

    ALLOCATE(operators_coefficients%variable_dual_vol_norm(nproma,nz_lev,nblks_v,1:no_dual_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating variable_dual_vol_norm failed')
    ENDIF

    !The last index "3" comes from the fact that we use a tridiagonal matrix
!   ALLOCATE(operators_coefficients%matrix_vert_diff_c(nproma,n_zlev,nblks_c, 3),&
!     & stat=ist)
!   IF (ist /= success) THEN
!     CALL finish ('mo_operator_ocean_coeff_3d',                 &
!       & 'allocation for matrix_vert_diff_c failed')
!   ENDIF
!   !The last index "3" comes from the fact that we use a tridiagonal matrix
!   ALLOCATE(operators_coefficients%matrix_vert_diff_e(nproma,n_zlev,nblks_e, 3),&
!     & stat=ist)
!   IF (ist /= success) THEN
!     CALL finish ('mo_operator_ocean_coeff_3d',                 &
!       & 'allocation for matrix_vert_diff_e failed')
!   ENDIF
    CALL add_var(var_list, 'matrix_vert_diff_c', operators_coefficients%matrix_vert_diff_c, &
      &          GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_SEA, &
      &          t_cf_var('matrix_vert_diff_c','','for each edge',DATATYPE_FLT64),&
      &          t_grib2_var(255,255,255,DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,n_zlev,nblks_c,no_primal_edges/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
    ALLOCATE(operators_coefficients%matrix_vert_diff_c_ptr(no_primal_edges))
    DO i=1,no_primal_edges
      WRITE(var_suffix,'(a,i1.1)') '_',i
      CALL add_ref( var_list, 'matrix_vert_diff_c', &
        &           'matrix_vert_diff_c'//TRIM(var_suffix), &
        &           operators_coefficients%matrix_vert_diff_c_ptr(i)%p,    &
        &           GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_SEA,&
        &           t_cf_var('matrix_vert_diff_c'//TRIM(var_suffix),'','', DATATYPE_FLT64), &
        &           t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
        &           ldims=(/nproma,n_zlev,nblks_c/),in_group=groups("oce_coeffs"))
    ENDDO

    CALL add_var(var_list, 'matrix_vert_diff_e', operators_coefficients%matrix_vert_diff_e, &
      &          GRID_UNSTRUCTURED_EDGE, ZA_DEPTH_BELOW_SEA, &
      &          t_cf_var('matrix_vert_diff_e','','for each cell',DATATYPE_FLT64),&
      &          t_grib2_var(255,255,255,DATATYPE_PACK16, GRID_REFERENCE, GRID_EDGE),&
      &          ldims=(/nproma,n_zlev,nblks_e,no_primal_edges/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    !
    ! initialize all components
    !
    DO ie = 1,3
      operators_coefficients%edge2cell_coeff_cc%x(ie)     = 0._wp
      operators_coefficients%edge2cell_coeff_cc_t%x(ie)   = 0._wp
      operators_coefficients%edge2vert_coeff_cc%x(ie)     = 0._wp
      operators_coefficients%edge2vert_coeff_cc_t%x(ie)   = 0._wp
      operators_coefficients%edge2vert_vector_cc%x(ie)    = 0._wp
      operators_coefficients%edge2cell_coeff_cc_dyn%x(ie) = 0._wp
      !operators_coefficients%edge2vert_coeff_cc_dyn%x(ie) = 0._wp
    END DO

    all_cells => patch_2D%cells%all
    all_edges => patch_2D%edges%all
    !all_verts => patch_2D%verts%all

    DO jk = 1, nz_lev
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO je =  i_startidx_e, i_endidx_e
!           operators_coefficients%edge_position_cc(je,jk,jb)             = gc2cc(patch_2D%edges%center(je,jb))
          operators_coefficients%edge_position_cc(je,jk,jb)             = patch_2D%edges%cartesian_center(je,jb)
          operators_coefficients%moved_edge_position_cc(je,jk,jb)%x(:)  = 0._wp
          operators_coefficients%upwind_cell_position_cc(je,jk,jb)%x(:) = 0._wp
        END DO
      END DO
    END DO

!    DO jk = 1, nz_lev
!      DO jb = all_cells%start_block, all_cells%end_block
!        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!        DO jc = i_startidx_c, i_endidx_c
!          operators_coefficients%cell_position_cc(jc,jk,jb) &
!            & = patch_2D%cells%cartesian_center(jc,jb)
! !             & = gc2cc(patch_2D%cells%center(jc,jb))
!
!        END DO
!      END DO
!    END DO

    operators_coefficients%edge2edge_viacell_coeff= 0._wp
    operators_coefficients%edge2edge_viavert_coeff= 0._wp

    operators_coefficients%fixed_vol_norm         = 0._wp
    operators_coefficients%variable_vol_norm      = 0._wp
    operators_coefficients%variable_dual_vol_norm = 0._wp

    operators_coefficients%dist_cell2edge = 0._wp

    operators_coefficients%div_coeff  = 0._wp
    operators_coefficients%rot_coeff  = 0._wp
    operators_coefficients%grad_coeff = 0._wp
    !operators_coefficients%n2s_coeff  = 0._wp
    !operators_coefficients%n2v_coeff  = 0._wp

    operators_coefficients%bnd_edge_idx = 0
    operators_coefficients%bnd_edge_blk = 0
    operators_coefficients%edge_idx     = 0
    operators_coefficients%orientation  = 0.0_wp
    operators_coefficients%bnd_edges_per_vertex= 0

    operators_coefficients%upwind_cell_idx = 1
    operators_coefficients%upwind_cell_blk = 1

    operators_coefficients%matrix_vert_diff_c= 0.0_wp
    operators_coefficients%matrix_vert_diff_e= 0.0_wp

    CALL message ('mo_operator_ocean_coeff_3d:allocate_operators_coefficients',&
      & 'memory allocation finished')

  END SUBROUTINE allocate_operators_coefficients
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Deallocation of operators coefficients.
  SUBROUTINE deallocate_operators_coefficients( operators_coefficients )
    ! !
    TYPE(t_operator_coeff), INTENT(inout) :: operators_coefficients

    !-----------------------------------------------------------------------
    DEALLOCATE(operators_coefficients%div_coeff)
    DEALLOCATE(operators_coefficients%grad_coeff)

    DEALLOCATE(operators_coefficients%rot_coeff)

    DEALLOCATE(operators_coefficients%n2s_coeff)
    DEALLOCATE(operators_coefficients%n2v_coeff)
    !
    DEALLOCATE(operators_coefficients%dist_cell2edge)

    DEALLOCATE(operators_coefficients%bnd_edge_idx)
    DEALLOCATE(operators_coefficients%bnd_edge_blk)
    DEALLOCATE(operators_coefficients%edge_idx)
    DEALLOCATE(operators_coefficients%orientation)
    DEALLOCATE(operators_coefficients%bnd_edges_per_vertex)
    DEALLOCATE(operators_coefficients%upwind_cell_idx)
    DEALLOCATE(operators_coefficients%upwind_cell_blk)

    DEALLOCATE(operators_coefficients%edge2edge_viacell_coeff)
    DEALLOCATE(operators_coefficients%edge2edge_viacell_coeff_top)
    DEALLOCATE(operators_coefficients%edge2edge_viacell_coeff_integrated)
    DEALLOCATE(operators_coefficients%edge2edge_viacell_coeff_all)

    DEALLOCATE(operators_coefficients%edge2cell_coeff_cc)

    DEALLOCATE(operators_coefficients%edge2cell_coeff_cc_dyn)

    DEALLOCATE(operators_coefficients%edge2cell_coeff_cc_t)

    DEALLOCATE(operators_coefficients%edge2vert_coeff_cc)

    DEALLOCATE(operators_coefficients%edge2vert_coeff_cc_t)
    DEALLOCATE(operators_coefficients%edge2vert_vector_cc)

    DEALLOCATE(operators_coefficients%edge2edge_viavert_coeff)

    DEALLOCATE(operators_coefficients%upwind_cell_position_cc)
    DEALLOCATE(operators_coefficients%moved_edge_position_cc)
    DEALLOCATE(operators_coefficients%edge_position_cc)

    DEALLOCATE(operators_coefficients%fixed_vol_norm)
    DEALLOCATE(operators_coefficients%variable_vol_norm)
    DEALLOCATE(operators_coefficients%variable_dual_vol_norm)

  END SUBROUTINE deallocate_operators_coefficients
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !> Initialize expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !!
  SUBROUTINE par_init_operator_coeff( patch_3D, operators_coefficients)
    !
    TYPE(t_patch_3D ),TARGET, INTENT(INOUT) :: patch_3D
    TYPE(t_operator_coeff),   INTENT(inout) :: operators_coefficients
   !
   !Local variables:
   TYPE(t_patch),POINTER :: patch_2D
    !-----------------------------------------------------------------------
    !TYPE(t_cartesian_coordinates) :: check_v(nproma, n_zlev, patch_2D%nblks_v, 6)
    !REAL(wp) :: check_r(nproma, n_zlev, patch_2D%nblks_c, 3)
    !REAL(wp) :: max_diff, max_val
    !-----------------------------------------------------------------------
    patch_2D => patch_3D%p_patch_2D(1)

    CALL init_operator_coeffs( patch_2D, operators_coefficients)

    CALL init_diff_operator_coeff_3D ( patch_2D, operators_coefficients )

    CALL apply_boundary2coeffs(patch_3D, operators_coefficients)

  END SUBROUTINE par_init_operator_coeff
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !> Initialize 3D expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  SUBROUTINE update_diffusion_matrices(   patch_3D,          &
                                         & p_phys_param,       &
                                         & matrix_vert_diff_e, &
                                         & matrix_vert_diff_c)
 
     TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
     TYPE (t_ho_params),       INTENT(IN)    :: p_phys_param
     REAL(wp), INTENT(INOUT) :: matrix_vert_diff_e(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e,1:3)
     REAL(wp), INTENT(INOUT) :: matrix_vert_diff_c(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks,1:3)
    !
    !
    ! matrix_vert_diff_c(:,:,:,1) : lower diagonal term
    ! matrix_vert_diff_c(:,:,:,2) : diagonal term
    ! matrix_vert_diff_c(:,:,:,3) : upper diagonal term
    !
    !Local variables
    !
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp), POINTER :: A_v(:,:,:)
    REAL(wp) :: inv_zinv_i(1:n_zlev)
    REAL(wp) :: inv_zinv_m(1:n_zlev)
    REAL(wp) :: dt_inv,dz_m(1:n_zlev),dz_i(1:n_zlev)
    INTEGER  :: je,jc,jb,jk,i_no_t
    INTEGER  :: slev,z_dolic
    INTEGER  :: i_startidx_e, i_endidx_e,  i_startidx_c, i_endidx_c
    TYPE(t_patch), POINTER :: patch_2D
    !---------------------------------------------------------
    patch_2D     => patch_3D%p_patch_2D(1)
    all_cells => patch_2D%cells%all
    all_edges => patch_2D%edges%all
    !---------------------------------------------------------
    slev   = 1
    dt_inv = 1.0_wp/dtime

    !The vertical diffusion matrices for the tracers
    DO i_no_t = 1,no_tracer

      A_v => p_phys_param%A_tracer_v(:,:,:, i_no_t)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)

          IF ( patch_3D%lsm_c(jc,1,jb) <= SEA_BOUNDARY ) THEN
            IF ( z_dolic >=MIN_DOLIC ) THEN

              inv_zinv_i(:) = patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,:,jb)
              inv_zinv_m(:) = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,:,jb)
                    dz_i(:) = patch_3D%p_patch_1D(1)%prism_center_dist_c(jc,:,jb)
                    dz_m(:) = patch_3D%p_patch_1D(1)%prism_thick_c(jc,:,jb)

              !first level
              matrix_vert_diff_c(jc,slev,jb,1) = 0.0_wp           
              matrix_vert_diff_c(jc,slev,jb,3) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)
              matrix_vert_diff_c(jc,slev,jb,2) = dt_inv- matrix_vert_diff_c(jc,slev,jb,3) 

              !Fill triangular matrix
              !b is diagonal a and c are upper and lower band
              !This corresponds to the 4th indices: "2", "1" and "3"
              DO jk = slev+1, z_dolic-1
                matrix_vert_diff_c(jc,jk,jb,1) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)
                matrix_vert_diff_c(jc,jk,jb,3) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)
                matrix_vert_diff_c(jc,jk,jb,2)  = dt_inv&
                                               & -matrix_vert_diff_c(jc,jk,jb,1)&
                                               & -matrix_vert_diff_c(jc,jk,jb,3)
              END DO
  
              !bottom
              matrix_vert_diff_c(jc,z_dolic,jb,1) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)
              matrix_vert_diff_c(jc,z_dolic,jb,3) = 0.0_wp
              matrix_vert_diff_c(jc,z_dolic,jb,2) = dt_inv - matrix_vert_diff_c(jc,z_dolic,jb,1)
            ENDIF
          ENDIF
        END DO
      END DO
    END DO 
   !---------------------------------------------------------


  !Now the velocity matrix
  A_v => p_phys_param%A_veloc_v

  DO jb = all_edges%start_block, all_edges%end_block
    CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
    DO je = i_startidx_e, i_endidx_e
      z_dolic = patch_3D%p_patch_1D(1)%dolic_e(je,jb)

      IF (  patch_3D%lsm_e(je,1,jb) <= SEA_BOUNDARY ) THEN
        IF ( z_dolic >= MIN_DOLIC ) THEN

          !inv_zinv_i(:)=1.0_wp/v_base%del_zlev_i(:)
          !inv_zinv_m(:)=1.0_wp/v_base%del_zlev_m(:)
          !inv_zinv_i(:) = p_os%p_diag%inv_prism_center_dist_e(je,:,jb)
          !inv_zinv_m(:) = p_os%p_diag%inv_prism_thick_e(je,:,jb)
          inv_zinv_i(:) = patch_3D%p_patch_1D(1)%inv_prism_center_dist_e(je,:,jb)
          inv_zinv_m(:) = patch_3D%p_patch_1D(1)%inv_prism_thick_e(je,:,jb)


          !Fill triangular matrix
          !b is diagonal a and c are upper and lower band
          DO jk = slev+1, z_dolic-1
            !a(jk) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)
            !c(jk) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)
            !b(jk) = dt_inv-a(jk)-c(jk)
            matrix_vert_diff_e(je,jk,jb,1) = -A_v(je,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)
            matrix_vert_diff_e(je,jk,jb,3) = -A_v(je,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)
            matrix_vert_diff_e(je,jk,jb,2)  = dt_inv&
                                            & -matrix_vert_diff_e(je,jk,jb,1)&
                                            & -matrix_vert_diff_e(je,jk,jb,3)

          END DO

          ! The first row
          !c(slev) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)!*dtime
          !a(slev) = 0.0_wp           
          !b(slev) = dt_inv- c(slev) !- a(slev) 
          matrix_vert_diff_e(je,slev,jb,3) = -A_v(je,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)
          matrix_vert_diff_e(je,slev,jb,1) = 0.0_wp           
          matrix_vert_diff_e(je,slev,jb,2) = dt_inv- matrix_vert_diff_e(je,slev,jb,3) 

          ! The last row
          !a(z_dolic) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)!*dtime
          !c(z_dolic) = 0.0_wp
          !b(z_dolic) = dt_inv - a(z_dolic)! - c(z_dolic)
          matrix_vert_diff_e(je,z_dolic,jb,1) = -A_v(je,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)
          matrix_vert_diff_e(je,z_dolic,jb,3) = 0.0_wp
          matrix_vert_diff_e(je,z_dolic,jb,2) = dt_inv - matrix_vert_diff_e(je,z_dolic,jb,1)
        ENDIF
      ENDIF
    END DO
  END DO
   END SUBROUTINE update_diffusion_matrices
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes the coefficients that determine the scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The coefficients are applied in module "mo_scalar_product".
  !! The following components of the data type "ocean_patch" are filled:
  !!   edge2cell_coeff  : coefficients for edge to cell mapping
  !!   edge2cell_coeff_t: coefficients for transposed of edge to cell mappings
  !!   edge2vert_coeff  : coefficients for edge to vertex mapping
  !!   edge2vert_coeff_t: coefficients for transposed of edge to vertex mappings
  !!   fixed_vol_norm   : summed volume weight of moved cell
  !!   variable_vol_norm: volume weight at the edges of moved cell
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M  2010-09
  !!  Modification by Stephan Lorenz, 2010-11
  !!
  SUBROUTINE init_operator_coeffs( patch_2D, operators_coefficients)
    TYPE(t_patch)    , TARGET, INTENT(INOUT)     :: patch_2D
    TYPE(t_operator_coeff),    INTENT(inout)     :: operators_coefficients

!Local variables
!
    REAL(wp)                      :: prime_edge_length      (1:nproma,patch_2D%nblks_e)
    REAL(wp)                      :: dual_edge_length       (1:nproma,patch_2D%nblks_e)
    REAL(wp)                      :: dist_cell2edge         (1:nproma,1:patch_2D%nblks_e,1:2)

    REAL(wp)                      :: div_coeff              (1:nproma,1:patch_2D%alloc_cell_blocks,1:no_primal_edges)
    REAL(wp)                      :: rot_coeff              (1:nproma,1:patch_2D%nblks_v,1:no_dual_edges)
    REAL(wp)                      :: grad_coeff             (1:nproma,1:patch_2D%nblks_e)

    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc     (1:nproma,1:patch_2D%alloc_cell_blocks,1:no_primal_edges)


    TYPE(t_subset_range), POINTER :: owned_edges         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_cells         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_verts         ! these are the owned entities
    TYPE(t_cartesian_coordinates) ::  edge_center,vertex_position!, vertex_center
    TYPE(t_cartesian_coordinates) :: dist_vector!, dist_vector_basic
    TYPE(t_cartesian_coordinates), POINTER :: dual_edge_middle(:,:)
    TYPE(t_cartesian_coordinates) :: coriolis_cartesian_coordinates
    TYPE(t_geographical_coordinates) :: coriolis_geo_coordinates, geo_coordinates
    REAL(wp) :: basin_center_lat_rad, basin_height_rad
    REAL(wp) :: length
    REAL(wp) :: inverse_sphere_radius
    !REAL(wp) :: dist_edge_cell, dist_edge_cell_basic
    !INTEGER :: edge_block_cell, edge_index_cell, ictr
    !INTEGER :: cell_edge, vert_edge
    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: cell_1_index, cell_1_block, cell_2_index, cell_2_block
    INTEGER :: vertex_1_index, vertex_1_block, vertex_2_index, vertex_2_block
    INTEGER :: level
    !-----------------------------------------------------------------------
    inverse_sphere_radius = 1.0_wp / grid_sphere_radius

    owned_edges => patch_2D%edges%owned
    owned_cells => patch_2D%cells%owned
    owned_verts => patch_2D%verts%owned

    !edge2vert_coeff_cc(:,:,:)%x(1) = 0.0_wp
    !edge2vert_coeff_cc(:,:,:)%x(2) = 0.0_wp
    !edge2vert_coeff_cc(:,:,:)%x(3) = 0.0_wp

    !edge2vert_coeff_cc_t(:,:,:)%x(1) = 0.0_wp
    !edge2vert_coeff_cc_t(:,:,:)%x(2) = 0.0_wp
    !edge2vert_coeff_cc_t(:,:,:)%x(3) = 0.0_wp

    !edge2edge_viacell_coeff(:,:,:) = 0.0_wp
    !edge2edge_viavert_coeff(:,:,:) = 0.0_wp
    rot_coeff(:,:,:)        = 0.0_wp
    div_coeff(:,:,:)        = 0.0_wp
    grad_coeff(:,:)         = 0.0_wp
    !-------------------------------------------
    ! compute some basic distances
    ! this is required if the cartesian distance is used
    ! instead of the spherical
    !
    ! computes_dist_cell2edge( patch_2D, intp_2D_coeff)
    !
    IF ( MID_POINT_DUAL_EDGE ) THEN
      dual_edge_middle => patch_2D%edges%cartesian_dual_middle
    ELSE
      dual_edge_middle => patch_2D%edges%cartesian_center
    ENDIF

    ! 1) calcultate prima and dual length as cartesian distance
    IF (LARC_LENGTH) THEN

      ! 1a) we just need to get them from the grid
      ! NOTE:  these are earth's distances, translate on a unit sphere
      dist_cell2edge(:,:,:) = &
        & patch_2D%edges%edge_cell_length(:,:,:) * inverse_sphere_radius
      prime_edge_length(:,:) = &
        & patch_2D%edges%primal_edge_length(:,:) * inverse_sphere_radius
      dual_edge_length(:,:) = &
        & patch_2D%edges%dual_edge_length(:,:) * inverse_sphere_radius

    ELSE

      !1b) calcultate prima and dual length as cartesian distance
      prime_edge_length(:,:) = 0.0_wp
      dual_edge_length (:,:) = 0.0_wp
      dist_cell2edge (:,:,:) =  0.0_wp

      DO edge_block = owned_edges%start_block, owned_edges%end_block
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index = start_index, end_index

          !----------------------------------------
          ! calculate the cartesian edge length
          vertex_1_index = patch_2D%edges%vertex_idx(edge_index, edge_block, 1)
          vertex_1_block = patch_2D%edges%vertex_blk(edge_index, edge_block, 1)
          vertex_2_index = patch_2D%edges%vertex_idx(edge_index, edge_block, 2)
          vertex_2_block = patch_2D%edges%vertex_blk(edge_index, edge_block, 2)

          dist_vector%x = &
            & patch_2D%verts%cartesian(vertex_1_index, vertex_1_block)%x - &
            & patch_2D%verts%cartesian(vertex_2_index, vertex_2_block)%x

            prime_edge_length(edge_index,edge_block) = &
              & SQRT(SUM((  dist_vector%x *  dist_vector%x)))
          !----------------------------------------

          !----------------------------------------
          ! calculate the cartesian distance of the edge center to the cell center
          DO neigbor = 1,2

            dist_cell2edge(edge_index,edge_block,neigbor) = 0.0_wp

            cell_index = patch_2D%edges%cell_idx(edge_index,edge_block,neigbor)
            cell_block = patch_2D%edges%cell_blk(edge_index,edge_block,neigbor)

            IF (cell_index > 0) THEN
              dist_vector%x = &
                & patch_2D%edges%cartesian_center(edge_index,edge_block)%x - &
                & patch_2D%cells%cartesian_center(cell_index,cell_block)%x

              dist_cell2edge(edge_index,edge_block,neigbor) = &
                & SQRT(SUM((  dist_vector%x *  dist_vector%x)))
            ENDIF

          ENDDO ! neigbor = 1,2
          !----------------------------------------

          !----------------------------------------
          ! calculate the cartesian dual edge length
          cell_1_index = patch_2D%edges%cell_idx(edge_index, edge_block, 1)
          cell_1_block = patch_2D%edges%cell_blk(edge_index, edge_block, 1)
          cell_2_index = patch_2D%edges%cell_idx(edge_index, edge_block, 2)
          cell_2_block = patch_2D%edges%cell_blk(edge_index, edge_block, 2)

          IF (cell_1_index > 0 .AND. cell_2_index > 0) THEN
            dist_vector%x = &
              & patch_2D%cells%cartesian_center(cell_1_index, cell_1_block)%x - &
              & patch_2D%cells%cartesian_center(cell_2_index, cell_2_block)%x

              dual_edge_length(edge_index,edge_block) = &
                & SQRT(SUM((  dist_vector%x *  dist_vector%x)))
          ELSE
              dual_edge_length(edge_index,edge_block) =              &
                & dist_cell2edge(edge_index,edge_block,1) + &
                & dist_cell2edge(edge_index,edge_block,2)
           ENDIF
          !----------------------------------------

        ENDDO ! edge_index=start_index,end_index
      ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block

      ! synchronize the edge distances
      CALL sync_patch_array(SYNC_E, patch_2D, dist_cell2edge(:,:,1))
      CALL sync_patch_array(SYNC_E, patch_2D, dist_cell2edge(:,:,2))
      CALL sync_patch_array(SYNC_E, patch_2D, prime_edge_length(:,:))
      CALL sync_patch_array(SYNC_E, patch_2D, dual_edge_length(:,:))
    ENDIF
    ! primal end dual edge lenght have been computed
    !-------------------------------------------



    !-------------------------------------------
    !2) calculate coefficients for difference operators
    !
    !2a) divergence
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        DO neigbor=1, patch_2D%cells%num_edges(cell_index,cell_block)!no_primal_edges

          edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x = 0.0_wp

          edge_index = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
          edge_block = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

          div_coeff(cell_index,cell_block,neigbor) =                           &
              & patch_2D%edges%primal_edge_length(edge_index,edge_block) *        &
              & patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)  / &
              & patch_2D%cells%area(cell_index,cell_block)

        ENDDO !neigbor=1,patch_2D%cell_type
      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
 

   !2b) gradient
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        grad_coeff(edge_index,edge_block)&
        & =1.0_wp/ dual_edge_length(edge_index,edge_block)!patch_2D%edges%inv_dual_edge_length(edge_index, edge_block)

      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
    

   !2c) curl coefficients
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index

        vertex_position%x = patch_2D%verts%cartesian(vertex_index, vertex_block)%x

        DO neigbor=1, patch_2D%verts%num_edges(vertex_index,vertex_block)!no_dual_edges

          edge_index = patch_2D%verts%edge_idx(vertex_index, vertex_block, neigbor)
          edge_block = patch_2D%verts%edge_blk(vertex_index, vertex_block, neigbor)

          IF (edge_block > 0) THEN
            rot_coeff(vertex_index,vertex_block,neigbor)           &
            &= patch_2D%edges%dual_edge_length(edge_index,edge_block) &
            &* patch_2D%verts%edge_orientation(vertex_index,vertex_block,neigbor)
          ENDIF
        ENDDO !neigbor=1,6
      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
   
    !Copy coefficients to 3D
    DO level=1,n_zlev
    operators_coefficients%div_coeff(:,level,:,:) = div_coeff(:,:,:)
    operators_coefficients%rot_coeff(:,level,:,:) = rot_coeff(:,:,:)
    operators_coefficients%grad_coeff(:,level,:)  = grad_coeff(:,:)
    END DO
    !-------------------
    ! sync the results
    CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%grad_coeff(:,:,:))
    DO neigbor=1,no_primal_edges
      CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%div_coeff(:,:,:,neigbor))
    END DO
    DO neigbor=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%rot_coeff(:,:,:,neigbor))
    END DO

    !----------------------------------------------------
    CALL init_operator_coeffs_cell( patch_2D, operators_coefficients,prime_edge_length, dual_edge_length )
    CALL init_operator_coeffs_vertex( patch_2D, operators_coefficients, prime_edge_length, dual_edge_length)

    !----------------------------------------------------
    ! 9) recalculate the coriolis coefficient
    ! It is required if we use the middle of the dual_edge_length
    IF (MID_POINT_DUAL_EDGE) THEN

      IF (CORIOLIS_TYPE == full_coriolis) THEN

        DO edge_block = owned_edges%start_block, owned_edges%end_block
          CALL get_index_range(owned_edges, edge_block, start_index, end_index)
          DO edge_index = start_index, end_index

             coriolis_geo_coordinates = cc2gc(dual_edge_middle(edge_index,edge_block))
             patch_2D%edges%f_e(edge_index,edge_block) = &
               & 2._wp * grid_angular_velocity * SIN(coriolis_geo_coordinates%lat)

          ENDDO
        ENDDO

      ELSEIF (CORIOLIS_TYPE == BETA_PLANE_CORIOLIS) THEN

        basin_center_lat_rad = basin_center_lat * deg2rad
        basin_height_rad     = basin_height_deg * deg2rad
        coriolis_geo_coordinates%lat = basin_center_lat_rad - 0.5_wp * basin_height_rad
        coriolis_geo_coordinates%lon = 0.0_wp
        coriolis_cartesian_coordinates  = gc2cc(coriolis_geo_coordinates)

        DO edge_block = owned_edges%start_block, owned_edges%end_block
          CALL get_index_range(owned_edges, edge_block, start_index, end_index)
          DO edge_index = start_index, end_index

          geo_coordinates     = cc2gc(dual_edge_middle(edge_index,edge_block))
          geo_coordinates%lon = 0.0_wp
          edge_center         = gc2cc(geo_coordinates)
          length              = grid_sphere_radius * &
            & arc_length(edge_center, coriolis_cartesian_coordinates)

          patch_2D%edges%f_e(edge_index,edge_block) =  2.0_wp * grid_angular_velocity * &
            & ( sin(basin_center_lat_rad) + (cos(basin_center_lat_rad) / &
            &   grid_sphere_radius) * length)

          ENDDO
        ENDDO

      ENDIF !(CORIOLIS_TYPE==full_coriolis)
    ENDIF ! (MID_POINT_DUAL_EDGE)
    !-------------------
    ! sync patch_2D%edges%f_e
    CALL sync_patch_array(SYNC_E, patch_2D, patch_2D%edges%f_e)
    !---------------------------------------------------------

  END SUBROUTINE init_operator_coeffs
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes the coefficients that determine the scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The coefficients are applied in module "mo_scalar_product".
  !! The following components of the data type "ocean_patch" are filled:
  !!   edge2cell_coeff  : coefficients for edge to cell mapping
  !!   edge2cell_coeff_t: coefficients for transposed of edge to cell mappings
  !!   edge2vert_coeff  : coefficients for edge to vertex mapping
  !!   edge2vert_coeff_t: coefficients for transposed of edge to vertex mappings
  !!   fixed_vol_norm   : summed volume weight of moved cell
  !!   variable_vol_norm: volume weight at the edges of moved cell
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M  2010-09
  !!  Modification by Stephan Lorenz, 2010-11
  SUBROUTINE init_operator_coeffs_cell( patch_2D, operators_coefficients,prime_edge_length, dual_edge_length )
    TYPE(t_patch), TARGET,  INTENT(INOUT)  :: patch_2D
    TYPE(t_operator_coeff), INTENT(inout)  :: operators_coefficients
    REAL(wp),               INTENT(IN)     :: prime_edge_length(1:nproma,1:patch_2D%nblks_e)
    REAL(wp),               INTENT(IN)     :: dual_edge_length (1:nproma,1:patch_2D%nblks_e)

!Local variables
!
    REAL(wp)                      :: edge2edge_viacell_coeff_2D(1:nproma,1:patch_2D%nblks_e,1:2*no_primal_edges)
    !REAL(wp)                      :: dist_cell2edge         (1:nproma,1:patch_2D%nblks_e,1:2)
    REAL(wp)                      :: fixed_vol_norm         (1:nproma,patch_2D%alloc_cell_blocks)
    REAL(wp)                      :: variable_vol_norm      (1:nproma,1:patch_2D%alloc_cell_blocks,1:no_primal_edges)
    REAL(wp)                      :: norm, orientation
    REAL(wp)                      :: dist_edge_cell, dist_edge_cell_basic

    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc     (1:nproma,1:patch_2D%alloc_cell_blocks,1:no_primal_edges)
    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc_t   (1:nproma,1:patch_2D%nblks_e,1:2)
    TYPE(t_cartesian_coordinates) :: cell_center, edge_center
    TYPE(t_cartesian_coordinates) :: dist_vector, dist_vector_basic

    INTEGER :: edge_block_cell, edge_index_cell, ictr
    INTEGER :: cell_edge
    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: level

    TYPE(t_subset_range), POINTER :: owned_edges, all_edges         
    TYPE(t_subset_range), POINTER :: owned_cells, all_cells        
    !-----------------------------------------------------------------------
    owned_edges => patch_2D%edges%owned
    all_edges   => patch_2D%edges%all
    owned_cells => patch_2D%cells%owned
    all_cells => patch_2D%cells%all
    !-------------------------------------------
    ! 3) compute:
    !   edge2cell_coeff_cc
    !   fixed_vol_norm
    !   variable_vol_norm
    edge2cell_coeff_cc(:,:,:)%x(1) = 0.0_wp
    edge2cell_coeff_cc(:,:,:)%x(2) = 0.0_wp
    edge2cell_coeff_cc(:,:,:)%x(3) = 0.0_wp

    edge2cell_coeff_cc_t(:,:,:)%x(1) = 0.0_wp
    edge2cell_coeff_cc_t(:,:,:)%x(2) = 0.0_wp
    edge2cell_coeff_cc_t(:,:,:)%x(3) = 0.0_wp

    fixed_vol_norm(:,:)       = 0.0_wp
    variable_vol_norm(:,:,:)  = 0.0_wp
    edge2edge_viacell_coeff_2D(:,:,:) = 0.0_wp

    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        cell_center%x = patch_2D%cells%cartesian_center(cell_index, cell_block)%x
        fixed_vol_norm(cell_index,cell_block) = 0.0_wp

        !-------------------------------
        DO neigbor=1, patch_2D%cells%num_edges(cell_index,cell_block)!no_primal_edges

          edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x = 0.0_wp
          variable_vol_norm(cell_index, cell_block, neigbor) =  0.0_wp

          edge_index = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
          edge_block = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

          IF (edge_block > 0 ) THEN
            ! we have an edge
            dist_vector%x = &
              & patch_2D%edges%cartesian_center(edge_index,edge_block)%x - &
              & cell_center%x

            norm  = SQRT(SUM( dist_vector%x * dist_vector%x))
            ! compute edge2cell_coeff_cc
            edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x =  &
              & dist_vector%x *                                             &
              & prime_edge_length(edge_index,edge_block) *                  &
              & patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)

            fixed_vol_norm(cell_index,cell_block) = &
              & fixed_vol_norm(cell_index,cell_block) + &
              & 0.5_wp * norm * prime_edge_length(edge_index,edge_block)

            variable_vol_norm(cell_index, cell_block, neigbor) = &
              & 0.5_wp * norm * prime_edge_length(edge_index,edge_block)

          ENDIF !(edge_block > 0 )
        ENDDO !neigbor=1,patch_2D%cell_type
        !-------------------------------
      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block


    !-------------------
    ! sync the results
    CALL sync_patch_array(SYNC_C, patch_2D, fixed_vol_norm(:,:))
    DO neigbor=1,patch_2D%cell_type
      CALL sync_patch_array(SYNC_C, patch_2D, edge2cell_coeff_cc(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_C, patch_2D, edge2cell_coeff_cc(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_C, patch_2D, edge2cell_coeff_cc(:,:,neigbor)%x(3))
      CALL sync_patch_array(SYNC_C, patch_2D, variable_vol_norm(:,:,neigbor))
    ENDDO
    !-------------------

   !copy calculated 2D arrays to 3D structure
    DO cell_block = all_cells%start_block, all_cells%end_block
      DO level = 1, n_zlev

       operators_coefficients%fixed_vol_norm(:,level,cell_block) = fixed_vol_norm(:,cell_block)

       DO neigbor=1,patch_2D%cell_type

         operators_coefficients%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(1)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(1)

         operators_coefficients%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(2)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(2)

         operators_coefficients%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(3)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(3)

         operators_coefficients%variable_vol_norm(:,level,cell_block,neigbor)  &
           &= variable_vol_norm(:,cell_block,neigbor)

        ENDDO ! neigbor=1,patch_2D%cell_type
      ENDDO  !  level = 1, n_zlev
    ENDDO ! cell_block
! no need for sync
!     CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%fixed_vol_norm(:,:,:))
!     DO neigbor=1,no_primal_edges
!       CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%edge2cell_coeff_cc(:,:,:,neigbor)%x(1))
!       CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%edge2cell_coeff_cc(:,:,:,neigbor)%x(2))
!       CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%edge2cell_coeff_cc(:,:,:,neigbor)%x(3))
!       CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%variable_vol_norm(:,:,:,neigbor))
!     ENDDO

    !-------------------------------------------
    ! 4) compute:
    !   edge2cell_coeff_cc_t

    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

!        edge2cell_coeff_cc_t(edge_index, edge_block, 2)%x = 0.0_wp
        edge_center%x = patch_2D%edges%cartesian_center(edge_index, edge_block)%x

        DO neigbor=1,2

 !         edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x = 0.0_wp
          cell_index = patch_2D%edges%cell_idx(edge_index, edge_block, neigbor)
          cell_block = patch_2D%edges%cell_blk(edge_index, edge_block, neigbor)

          IF (cell_index > 0) THEN

            dist_vector%x =  edge_center%x -                             &
              patch_2D%cells%cartesian_center(cell_index, cell_block)%x

            orientation = DOT_PRODUCT(dist_vector%x, &
              & patch_2D%edges%primal_cart_normal(edge_index, edge_block)%x)
            IF (orientation < 0.0_wp) dist_vector%x = - dist_vector%x

            edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x = &
              dist_vector%x / dual_edge_length(edge_index, edge_block)

          ENDIF ! (cell_block > 0)
        ENDDO ! neigbor=1,2
      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
    !-------------------
    ! sync the results
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch_2D, edge2cell_coeff_cc_t(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2cell_coeff_cc_t(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2cell_coeff_cc_t(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2
    !   edge2cell_coeff_cc_t is computed

    !copy 2D to 3D structure
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index =  start_index, end_index
        DO level = 1, n_zlev

          operators_coefficients%edge2cell_coeff_cc_t(edge_index,level,edge_block,1)%x &
          &= edge2cell_coeff_cc_t(edge_index,edge_block,1)%x

          operators_coefficients%edge2cell_coeff_cc_t(edge_index,level,edge_block,2)%x &
          &= edge2cell_coeff_cc_t(edge_index,edge_block,2)%x

        ENDDO
      ENDDO
    ENDDO
    ! sync the results
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2
    !-------------------------------------------


    !-------------------------------------------
    ! 5) compute
    !   calculate edge2edge_viacell_coeff_2D
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge_center%x = patch_2D%edges%cartesian_center(edge_index, edge_block)%x

        !ictr=0
        DO neigbor=1,2

          cell_index    = patch_2D%edges%cell_idx(edge_index, edge_block, neigbor)
          cell_block    = patch_2D%edges%cell_blk(edge_index, edge_block, neigbor)

          IF (cell_index <= 0) CYCLE

          IF(neigbor==1) ictr = 0
          IF(neigbor==2) ictr = no_primal_edges

          cell_center%x = patch_2D%cells%cartesian_center(cell_index, cell_block)%x

          !dist_vector_basic%x = edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x
          dist_vector_basic%x = edge_center%x - cell_center%x

          dist_edge_cell_basic  = SQRT(SUM( dist_vector_basic%x * dist_vector_basic%x))
          dist_vector_basic%x = dist_vector_basic%x/dist_edge_cell_basic

          orientation = DOT_PRODUCT(dist_vector_basic%x, &
          & patch_2D%edges%primal_cart_normal(edge_index,edge_block)%x)
          IF (orientation < 0.0_wp) dist_vector_basic%x = - dist_vector_basic%x

          !loop over the edges of neighbor 1 and 2
          DO cell_edge=1,patch_2D%cells%num_edges(cell_index,cell_block)!no_primal_edges!patch_2D%cell_type

            ictr=ictr+1
            !actual edge
            edge_index_cell = patch_2D%cells%edge_idx(cell_index, cell_block, cell_edge)
            edge_block_cell = patch_2D%cells%edge_blk(cell_index, cell_block, cell_edge)

            !dist_vector%x = edge2cell_coeff_cc(cell_index,cell_block,cell_edge)%x
            dist_vector%x =  patch_2D%edges%cartesian_center(edge_index_cell, edge_block_cell)%x  &
            & -cell_center%x

            dist_edge_cell  = SQRT(SUM( dist_vector%x * dist_vector%x))
            dist_vector%x = dist_vector%x/dist_edge_cell
            dist_vector%x = dist_vector%x*patch_2D%cells%edge_orientation(cell_index,cell_block,cell_edge)

            !This is the cosine of the angle between vectors from cell center
            !to cell edges
            edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr)&
            & =DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)

            !IF(abs(edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr)-1.0_wp)<1.0E-6_wp)THEN
            !  write(*,*)'ran into'
              !edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr)=1.0_wp
            !ENDIF

            !multiply the cosine by length and orientation and divide by
            !dual length
            edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr)=        &
              &edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr)        &
              &*prime_edge_length(edge_index_cell,edge_block_cell)        &
              &* dist_edge_cell *dist_edge_cell_basic                     &
              &/dual_edge_length(edge_index, edge_block)

! IF(edge_index==1.and.edge_block==1)THEN
! write(123,*)'actual angle',neigbor, edge_index_cell, edge_block_cell,ictr,&
! & edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr),&
! &DOT_PRODUCT(edge2cell_coeff_cc(cell_index,cell_block,cell_edge)%x,edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x),&
! &acos(edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr))*rad2deg
! !IF(edge_index_cell==edge_index.and.edge_block_cell==edge_block)THEN
! !write(123,*)'vecs',neigbor,dist_vector_basic%x,dist_vector%x
! !ENDIF
! ENDIF
          END DO
        ENDDO ! neigbor=1,2
      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block

    ! these coeffecients will not be used for-non owened edges,
    ! sync only for safety
    DO ictr=1, 2*no_primal_edges
      CALL sync_patch_array(SYNC_E, patch_2D, edge2edge_viacell_coeff_2D(:,:,ictr))
    ENDDO

    !copy 2D to 3D structure
    DO edge_block = all_edges%start_block, all_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(all_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          operators_coefficients%edge2edge_viacell_coeff(edge_index,level,edge_block,1:2*no_primal_edges) = &
            & edge2edge_viacell_coeff_2D(edge_index,edge_block,1:2*no_primal_edges)

        ENDDO
      ENDDO
    ENDDO    
    !DO neigbor=1, 2*no_primal_edges
    !  CALL sync_patch_array(SYNC_E, patch_2D,  operators_coefficients%edge2edge_viacell_coeff_2D(:,:,:,neigbor))
    !ENDDO
   !-------------------------------------------
!Do ictr=1,12
! write(*,*)'max coeff',&
! & maxval(edge2edge_viacell_coeff(:,:,ictr))
! write(*,*)'min coeff',&
! & minval(edge2edge_viacell_coeff(:,:,ictr))
! write(*,*)'max angle',&
! & maxval(acos(edge2edge_viacell_coeff(:,:,ictr))*rad2deg)
! write(*,*)'min angle',&
! & minval(acos(edge2edge_viacell_coeff(:,:,ictr))*rad2deg)
!END DO
!-------------------------------------------
  END SUBROUTINE init_operator_coeffs_cell
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes the coefficients that determine the scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The coefficients are applied in module "mo_scalar_product".
  !! The following components of the data type "ocean_patch" are filled:
  !!   edge2cell_coeff  : coefficients for edge to cell mapping
  !!   edge2cell_coeff_t: coefficients for transposed of edge to cell mappings
  !!   edge2vert_coeff  : coefficients for edge to vertex mapping
  !!   edge2vert_coeff_t: coefficients for transposed of edge to vertex mappings
  !!   fixed_vol_norm   : summed volume weight of moved cell
  !!   variable_vol_norm: volume weight at the edges of moved cell
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M  2010-09
  !!  Modification by Stephan Lorenz, 2010-11
  !!
  SUBROUTINE init_operator_coeffs_vertex( patch_2D, operators_coefficients, prime_edge_length, dual_edge_length)
    TYPE(t_patch), TARGET, INTENT(INOUT) :: patch_2D
    TYPE(t_operator_coeff),INTENT(INOUT) :: operators_coefficients
    REAL(wp), INTENT(IN)                 :: prime_edge_length(1:nproma,1:patch_2D%nblks_e)
    REAL(wp), INTENT(IN)                 :: dual_edge_length (1:nproma,1:patch_2D%nblks_e)

!Local variables
!
    TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc     (1:nproma,1:patch_2D%nblks_v,1:no_dual_edges)
    TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc_t   (1:nproma,1:patch_2D%nblks_e,1:2)
    TYPE(t_cartesian_coordinates) :: vertex_position, edge_center, vertex_center
    TYPE(t_cartesian_coordinates) :: dist_vector, dist_vector_basic
    TYPE(t_cartesian_coordinates), POINTER :: dual_edge_middle(:,:)

    REAL(wp)                      :: edge2edge_viavert_coeff(1:nproma,1:patch_2D%nblks_e,1:2*no_dual_edges )
    !REAL(wp)                      :: variable_dual_vol_norm (1:nproma,1:patch_2D%nblks_e,1:no_dual_edges)
    REAL(wp)                      :: norm, orientation, length

    INTEGER :: ictr,edge_block_cell, edge_index_cell
    INTEGER :: vert_edge
    INTEGER :: edge_block, edge_index
    !INTEGER :: cell_index, cell_block
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: level

    TYPE(t_subset_range), POINTER :: owned_edges, all_edges        
    TYPE(t_subset_range), POINTER :: owned_verts         
    !-----------------------------------------------------------------------
    owned_edges => patch_2D%edges%owned
    all_edges   => patch_2D%edges%all
    owned_verts => patch_2D%verts%owned

    IF ( MID_POINT_DUAL_EDGE ) THEN
      dual_edge_middle => patch_2D%edges%cartesian_dual_middle
    ELSE
      dual_edge_middle => patch_2D%edges%cartesian_center
    ENDIF

    edge2vert_coeff_cc(:,:,:)%x(1) = 0.0_wp
    edge2vert_coeff_cc(:,:,:)%x(2) = 0.0_wp
    edge2vert_coeff_cc(:,:,:)%x(3) = 0.0_wp

    edge2vert_coeff_cc_t(:,:,:)%x(1) = 0.0_wp
    edge2vert_coeff_cc_t(:,:,:)%x(2) = 0.0_wp
    edge2vert_coeff_cc_t(:,:,:)%x(3) = 0.0_wp
    edge2edge_viavert_coeff(:,:,:)   = 0.0_wp
   !-------------------------------------------
    ! 1) compute:
    !   edge2vert_coeff_cc
    !   variable_dual_vol_norm will be handled in boundary_coeff-sbr
    !variable_dual_vol_norm(:,:,:)  = 0.0_wp

    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index

        vertex_position%x = patch_2D%verts%cartesian(vertex_index, vertex_block)%x

        DO neigbor=1, patch_2D%verts%num_edges(vertex_index,vertex_block)  !no_dual_edges

          edge_index = patch_2D%verts%edge_idx(vertex_index, vertex_block, neigbor)
          edge_block = patch_2D%verts%edge_blk(vertex_index, vertex_block, neigbor)

          IF (edge_block > 0) THEN
            ! we got an adjacent edge
            dist_vector%x = &
              dual_edge_middle(edge_index, edge_block)%x - &
              vertex_position%x

            ! the dist_vector has cartesian length
            ! if we use spherical distance we need to recalculate
            ! its length
            IF (LARC_LENGTH) THEN
              length = arc_length(vertex_position, dual_edge_middle(edge_index, edge_block))
              norm = SQRT(SUM( dist_vector%x * dist_vector%x ))
              dist_vector%x = dist_vector%x * length / norm
            ELSE
              length = SQRT(SUM( dist_vector%x * dist_vector%x ))
            ENDIF

            dist_vector = vector_product(dist_vector, dual_edge_middle(edge_index, edge_block))
            orientation = DOT_PRODUCT( dist_vector%x,                         &
               & patch_2D%edges%primal_cart_normal(edge_index, edge_block)%x)
            IF (orientation < 0.0_wp) dist_vector%x = - dist_vector%x

              edge2vert_coeff_cc(vertex_index, vertex_block, neigbor)%x = &
              & dist_vector%x                                *                    &
              & dual_edge_length(edge_index, edge_block) 
          ENDIF !(edge_block > 0) THEN

          !rot_coeff(vertex_index,vertex_block,neigbor)     &
          !    &= patch_2D%edges%dual_edge_length(edge_index,edge_block) * &
          !    & patch_2D%verts%edge_orientation(vertex_index,vertex_block,neigbor)

        ENDDO !neigbor=1,6
      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
    !-------------------
    ! sync the results
    DO neigbor=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch_2D, edge2vert_coeff_cc(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_V, patch_2D, edge2vert_coeff_cc(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_V, patch_2D, edge2vert_coeff_cc(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,6
    ! edge2vert_coeff_cc is computed

    !copy 2D to 3D structure
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      DO level = 1, n_zlev
        DO neigbor=1,no_dual_edges

          operators_coefficients%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(1)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(1)

          operators_coefficients%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(2)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(2)

          operators_coefficients%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(3)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(3)

!          operators_coefficients%variable_dual_vol_norm(:,level,vertex_block,neigbor)&
!          &=variable_dual_vol_norm(:,vertex_block,neigbor)
        ENDDO ! neigbor=1,patch_2D%cell_type
      ENDDO  !  level = 1, n_zlev
    ENDDO ! vertex_block
    DO neigbor=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,6
    !----------------------------------------------------

    !----------------------------------------------------
    ! 7) compute:
    !   edge2vert_coeff_cc_t
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge_center%x = dual_edge_middle(edge_index, edge_block)%x

        DO neigbor=1,2

          edge2vert_coeff_cc_t(edge_index, edge_block, neigbor)%x = 0.0_wp

          vertex_index = patch_2D%edges%vertex_idx(edge_index, edge_block, neigbor)
          vertex_block = patch_2D%edges%vertex_blk(edge_index, edge_block, neigbor)

          edge2vert_coeff_cc_t(edge_index, edge_block, neigbor)%x =              &
            & (edge_center%x - patch_2D%verts%cartesian(vertex_index, vertex_block)%x) * &
            & patch_2D%edges%system_orientation(edge_index, edge_block)                / &
            & prime_edge_length(edge_index, edge_block)

        ENDDO !neigbor=1,2

      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
    !-------------------
    ! sync the results
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2
    ! edge2vert_coeff_cc_t is computed


    DO edge_block = owned_edges%start_block, owned_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          operators_coefficients%edge2vert_coeff_cc_t(edge_index,level,edge_block,1)%x &
          &=  edge2vert_coeff_cc_t(edge_index,edge_block,1)%x

          operators_coefficients%edge2vert_coeff_cc_t(edge_index,level,edge_block,2)%x &
          &=  edge2vert_coeff_cc_t(edge_index,edge_block,2)%x

        ENDDO
      ENDDO
    ENDDO
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2

    !----------------------------------------------------
    ! 8) compute
    !   edge2edge_viavert calculation
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge_center%x = dual_edge_middle(edge_index, edge_block)%x

        DO neigbor=1,2

          vertex_index   = patch_2D%edges%vertex_idx(edge_index, edge_block, neigbor)
          vertex_block   = patch_2D%edges%vertex_blk(edge_index, edge_block, neigbor)
          vertex_center%x= patch_2D%verts%cartesian(vertex_index, vertex_block)%x

          dist_vector_basic%x = (edge_center%x - vertex_center%x)

          !IF(neigbor==1)ictr = 0
          !IF(neigbor==2)ictr = no_dual_edges

          ictr = (neigbor - 1)*no_dual_edges

          DO vert_edge=1,patch_2D%verts%num_edges(vertex_index,vertex_block)!no_dual_edges
            ictr=ictr+1
            !actual edge
            edge_index_cell = patch_2D%verts%edge_idx(vertex_index, vertex_block, vert_edge)
            edge_block_cell = patch_2D%verts%edge_blk(vertex_index, vertex_block, vert_edge)
            dist_vector%x  =  dual_edge_middle(edge_index_cell, edge_block_cell)%x - vertex_center%x

            dist_vector = vector_product(dist_vector, dual_edge_middle(edge_index_cell, edge_block_cell))
            orientation = DOT_PRODUCT( dist_vector%x,                         &
               & patch_2D%edges%primal_cart_normal(edge_index_cell, edge_block_cell)%x)
            ! orientation should not be 0, since this would mean that the prime and dual are parallel
            ! overall this calculation should be derived from the verts%edge_orientation
            ! orientation will recieve a value -1, or 1 based on the previous,
            ! then multuplied by -1 if neigbor=2, otherwise unchanged
            orientation = SIGN(1.0_wp,orientation) * (3.0_wp - 2.0_wp * REAL(neigbor,wp))

            !The dot product is the cosine of the angle between vectors from dual cell centers
            !to dual cell edges 
            edge2edge_viavert_coeff(edge_index,edge_block,ictr)         &
              & = orientation                                           &
              & * DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)        &
              & * patch_2D%edges%system_orientation(edge_index, edge_block)&
              & * (dual_edge_length(edge_index_cell, edge_block_cell)   &
              &    / prime_edge_length(edge_index, edge_block))

          END DO
        ENDDO !neigbor=1,2
      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block

    !-------------------
    DO ictr=1, 2*no_dual_edges
!      write(0,*)'ictr:',ictr
      CALL sync_patch_array(SYNC_E, patch_2D, edge2edge_viavert_coeff(:,:,ictr))
    ENDDO

    DO edge_block = all_edges%start_block, all_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(all_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          operators_coefficients%edge2edge_viavert_coeff(edge_index,level,edge_block,1:2*no_dual_edges) = &
            & edge2edge_viavert_coeff(edge_index,edge_block,1:2*no_dual_edges)

        ENDDO
      ENDDO
    ENDDO
!    DO neigbor=1,2*no_dual_edges
!      CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2edge_viavert_coeff(:,:,:,neigbor))
!    END DO
   !-------------------------------------------
  END SUBROUTINE init_operator_coeffs_vertex
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Modify coefficients at the boundary
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !!
  SUBROUTINE apply_boundary2coeffs( patch_3D, operators_coefficients)
    ! !
    TYPE(t_patch_3D ),TARGET, INTENT(INOUT) :: patch_3D
    TYPE(t_operator_coeff), INTENT(INOUT)   :: operators_coefficients

    !Local variables
    INTEGER :: jk, jc, jb, je, ibe, ile, jev, jv,ie
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_v, i_endidx_v

    INTEGER :: sea_edges_per_vertex(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: ibnd_edge_idx(4), ibnd_edge_blk(4)  !maximal 4 boundary edges in a dual loop.
    INTEGER :: i_edge_idx(4)
    !REAL(wp) :: z_orientation(4)!,z_area_scaled 
    REAL(wp) :: zarea_fraction(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: icell_idx_1, icell_blk_1
    INTEGER :: icell_idx_2, icell_blk_2
    INTEGER :: boundary_counter
    !INTEGER :: cell_index, cell_block
    !INTEGER :: edge_index_cell, edge_block_cell
    INTEGER :: neigbor, k
    INTEGER :: cell_index, cell_block, edge_index_of_cell, edge_block_of_cell, k_coeff

    !INTEGER :: vertex_edge
    TYPE(t_cartesian_coordinates) :: cell1_cc, cell2_cc, vertex_cc

    TYPE(t_subset_range), POINTER :: all_edges, owned_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: owned_verts!, in_domain_verts
    TYPE(t_patch), POINTER        :: patch_2D
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:apply_boundary2coeffs')
    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')
    patch_2D    => patch_3D%p_patch_2D(1)
    all_cells   => patch_2D%cells%all
    all_edges   => patch_2D%edges%all
    owned_edges => patch_2D%edges%owned
    owned_verts => patch_2D%verts%owned
    !in_domain_verts  => patch_3D%p_patch_2D(1)%verts%in_domain

    sea_edges_per_vertex(:,:,:)                       = 0
    zarea_fraction(1:nproma,1:n_zlev,1:patch_2D%nblks_v) = 0.0_wp

    !-------------------------------------------------------------
    !0. check the coefficients for edges, these are:
    !     grad_coeff
    !     edge2cell_coeff_cc_t
    !     edge2vert_coeff_cc_t
    !
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO jk = 1, n_zlev
        DO je = i_startidx_e, i_endidx_e

          IF ( patch_3D%lsm_e(je,jk,jb) /= sea ) THEN
            operators_coefficients%grad_coeff          (je,jk,jb) = 0.0_wp
            operators_coefficients%edge2cell_coeff_cc_t(je,jk,jb,1)%x(1:3) = 0.0_wp
            operators_coefficients%edge2cell_coeff_cc_t(je,jk,jb,2)%x(1:3) = 0.0_wp
            operators_coefficients%edge2vert_coeff_cc_t(je,jk,jb,1)%x(1:3) = 0.0_wp
            operators_coefficients%edge2vert_coeff_cc_t(je,jk,jb,2)%x(1:3) = 0.0_wp
          ENDIF
        ENDDO
      END DO
    END DO
    !-------------------------------------------------------------
    !All the coefficients "edge2edge_viacell_coeff" are set to zero
    !if the actual edge under consideration is an land edge.
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO jk = 1, n_zlev
        DO je = i_startidx_e, i_endidx_e
          IF ( patch_3D%lsm_e(je,jk,jb) /= sea ) THEN
            operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,:) = 0.0_wp
            operators_coefficients%edge2edge_viavert_coeff(je,jk,jb,:) = 0.0_wp
          ENDIF
        END DO
      END DO
    END DO
    !-------------------------------------------------------------
    !The coefficients "edge2edge_viacell_coeff" are set to zero
    !for which the stencil contains is an land edge.
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO jk = 1, n_zlev
        DO je = i_startidx_e, i_endidx_e

          !Handle neighbour cell 1
          ictr  = 0
          il_c  = patch_2D%edges%cell_idx(je,jb,1)
          ib_c  = patch_2D%edges%cell_blk(je,jb,1)
          IF (il_c > 0) THEN
            DO ie=1, no_primal_edges
              ictr =ictr+1
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

              IF ( patch_3D%lsm_e(il_e,jk,ib_e) /= sea ) THEN
                operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr)=0.0_wp
              ENDIF
            END DO
          ENDIF

          !Handle neighbour cell 2
          ictr  = no_primal_edges
          il_c  = patch_2D%edges%cell_idx(je,jb,2)
          ib_c  = patch_2D%edges%cell_blk(je,jb,2)
          IF (il_c > 0) THEN
            DO ie=1, no_primal_edges
              ictr =ictr+1
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

              IF ( patch_3D%lsm_e(il_e,jk,ib_e) /= sea ) THEN
                operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr)=0.0_wp
              ENDIF
            END DO
          ENDIF
        END DO
      END DO
    END DO 

    !-------------------------------------------------------------
    ! Normalize "edge2edge_viacell_coeff" are set to zero
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO jk = 1, n_zlev
        DO je = i_startidx_e, i_endidx_e
          IF ( patch_3D%lsm_e(je,jk,jb) == sea ) THEN
            icell_idx_1 = patch_2D%edges%cell_idx(je,jb,1)
            icell_blk_1 = patch_2D%edges%cell_blk(je,jb,1)

            icell_idx_2 = patch_2D%edges%cell_idx(je,jb,2)
            icell_blk_2 = patch_2D%edges%cell_blk(je,jb,2)

            operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,1:no_primal_edges) &
            &= operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,1:no_primal_edges)&
            &/operators_coefficients%fixed_vol_norm(icell_idx_1,jk,icell_blk_1)

            operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,no_primal_edges+1:2*no_primal_edges) &
            &= operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,no_primal_edges+1:2*no_primal_edges)&
            &/operators_coefficients%fixed_vol_norm(icell_idx_2,jk,icell_blk_2)
          ENDIF

        END DO
      END DO
    END DO 

    !-------------------------------------------------------------
    !Fill edge2edge_viacell_coeff_top and edge2edge_viacell_coeff_integrated from edge2edge_viacell_coeff
    !
    ! the top is just the first level rearranged:
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e
        DO k=1, 2 * no_primal_edges
          operators_coefficients%edge2edge_viacell_coeff_top(k, je, jb) = &
             operators_coefficients%edge2edge_viacell_coeff(je, 1, jb, k)
        ENDDO
      ENDDO
    ENDDO

    ! the integrated are the rest of levels > 1, weighted by prism_thick_e and integrated:
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e

        ! zero for two cells/ six edges
        DO k_coeff = 1, 2 * no_primal_edges
           operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, jb) = 0.0_wp
        ENDDO

        ! the first cell
        cell_index   = patch_2D%edges%cell_idx(je,jb,1)
        cell_block = patch_2D%edges%cell_blk(je,jb,1)
        IF (cell_index > 0) THEN ! in case it's a lateral boundary edge
          DO k=1, no_primal_edges
            edge_index_of_cell = patch_2D%cells%edge_idx(cell_index, cell_block, k)
            edge_block_of_cell = patch_2D%cells%edge_blk(cell_index, cell_block, k)
            k_coeff = k
            DO jk=2, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
              operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, jb) = &
                 operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, jb) + &
                 operators_coefficients%edge2edge_viacell_coeff(je, jk, jb, k_coeff) * &
                 patch_3D%p_patch_1D(1)%prism_thick_e(edge_index_of_cell, jk, edge_block_of_cell)

           !   write(0,*) jb, je, jk, k_coeff, operators_coefficients%edge2edge_viacell_coeff(je, jk, jb, k_coeff), &
           !     & patch_3D%p_patch_1D(1)%prism_thick_e(edge_index_of_cell, jk, edge_block_of_cell), &
           !    &  operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, jb)
            ENDDO
          ENDDO
        ENDIF

        ! the second cell
        cell_index   = patch_2D%edges%cell_idx(je,jb,2)
        cell_block = patch_2D%edges%cell_blk(je,jb,2)
        IF (cell_index > 0) THEN ! in case it's a lateral boundary edge
          DO k=1, no_primal_edges
            edge_index_of_cell = patch_2D%cells%edge_idx(cell_index, cell_block, k)
            edge_block_of_cell = patch_2D%cells%edge_blk(cell_index, cell_block, k)
            k_coeff = no_primal_edges + k
            DO jk=2, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
              operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, jb) = &
                 operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, jb) + &
                 operators_coefficients%edge2edge_viacell_coeff(je, jk, jb, k_coeff) * &
                 patch_3D%p_patch_1D(1)%prism_thick_e(edge_index_of_cell, jk, edge_block_of_cell)

            !  write(0,*) jb, je, jk, k_coeff, operators_coefficients%edge2edge_viacell_coeff(je, jk, jb, k_coeff), &
            !    & patch_3D%p_patch_1D(1)%prism_thick_e(edge_index_of_cell, jk, edge_block_of_cell), &
            !    &  operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, jb)
            ENDDO
          ENDDO
        ENDIF

      ENDDO
    ENDDO

    ! zeros edge2edge_viacell_coeff_all
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e
        DO k=1, 2 * no_primal_edges
          operators_coefficients%edge2edge_viacell_coeff_all(k, je, jb) = 0.0_wp
        ENDDO
      ENDDO
    ENDDO
    !-------------------------------------------------------------
    !1) Set coefficients for div and grad to zero at boundary edges
    ! Also for operators_coefficients%edge2cell_coeff_cc
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jk=1,n_zlev
        DO jc = i_startidx_c, i_endidx_c
          DO je = 1, patch_2D%cells%num_edges(jc,jb)!no_primal_edges!

            ile = patch_2D%cells%edge_idx(jc,jb,je)
            ibe = patch_2D%cells%edge_blk(jc,jb,je)

            IF ( patch_3D%lsm_e(ile,jk,ibe) /= sea ) THEN
              operators_coefficients%div_coeff(jc,jk,jb,je) = 0.0_wp
              operators_coefficients%edge2cell_coeff_cc(jc,jk,jb,je)%x(1:3) = 0.0_wp
            ENDIF
          ENDDO ! je = 1, patch_2D%cells%num_edges(jc,jb)
        ENDDO ! jc = i_startidx_c, i_endidx_c
      END DO ! jk=1,n_zlev
    END DO ! jb = all_cells%start_block, all_cells%end_block
    !-------------------------------------------------------------
    !2) prepare coefficients for rot at boundary edges
    ! this is done on owned vertices, so sync is required
    DO jb = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, jb, i_startidx_v, i_endidx_v)
      DO jk = 1, n_zlev
        DO jv = i_startidx_v, i_endidx_v

          ibnd_edge_idx(1:4)      = 0
          ibnd_edge_blk(1:4)      = 0
          !i_edge_idx(1:4)         = 0
          !z_orientation(1:4)      = 0.0_wp
          boundary_counter        = 0

          DO jev = 1, patch_2D%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = patch_2D%verts%edge_idx(jv,jb,jev)
            ibe = patch_2D%verts%edge_blk(jv,jb,jev)
            !Check, if edge is sea or boundary edge and take care of dummy edge
            ! edge with indices ile, ibe is sea edge
            ! edge with indices ile, ibe is boundary edge

            IF ( patch_3D%lsm_e(ile,jk,ibe) == SEA ) THEN
              sea_edges_per_vertex(jv,jk,jb) = sea_edges_per_vertex(jv,jk,jb) + 1
            ELSEIF ( patch_3D%lsm_e(ile,jk,ibe) == BOUNDARY ) THEN

              !increase boundary edge counter
              boundary_counter = boundary_counter + 1

              operators_coefficients%bnd_edges_per_vertex(jv,jk,jb) &
                & = operators_coefficients%bnd_edges_per_vertex(jv,jk,jb) +1

              IF (boundary_counter > 4) THEN
                !maximal 4 boundary edges per dual loop are allowed: somethings wrong with the grid
                CALL message (TRIM('sbr nonlinear Coriolis'), &
                  & 'more than 4 boundary edges per dual loop: something is wrong with the grid')
                CALL finish (routine,'Grid-boundary error !!')
              ENDIF
              ibnd_edge_idx(boundary_counter) = ile
              ibnd_edge_blk(boundary_counter) = ibe
              !z_orientation(boundary_counter) = patch_2D%verts%edge_orientation(jv,jb,jev)
              i_edge_idx(boundary_counter)    = jev

              operators_coefficients%bnd_edge_idx(jv,jk,jb,boundary_counter)= ile
              operators_coefficients%bnd_edge_blk(jv,jk,jb,boundary_counter)= ibe
              operators_coefficients%orientation(jv,jk,jb,boundary_counter) = &
                & patch_2D%verts%edge_orientation(jv,jb,jev)
              operators_coefficients%edge_idx(jv,jk,jb,boundary_counter)    = jev

            END IF
          END DO ! jev = 1, patch_2D%verts%num_edges(jv,jb)

          IF( MOD(boundary_counter,2) /= 0 ) THEN
            CALL finish (routine,'MOD(boundary_counter,2) /= 0 !!')
          ENDIF

          !---------------------------------------------------------------------------------
          !Modified area calculation
          vertex_cc = patch_2D%verts%cartesian(jv,jb)
          DO jev = 1, patch_2D%verts%num_edges(jv,jb)
            ! get line and block indices of edge jev around vertex jv
            ile = patch_2D%verts%edge_idx(jv,jb,jev)
            ibe = patch_2D%verts%edge_blk(jv,jb,jev)
            !get neighbor cells
            icell_idx_1 = patch_2D%edges%cell_idx(ile,ibe,1)
            icell_idx_2 = patch_2D%edges%cell_idx(ile,ibe,2)
            icell_blk_1 = patch_2D%edges%cell_blk(ile,ibe,1)
            icell_blk_2 = patch_2D%edges%cell_blk(ile,ibe,2)

            IF ( patch_3D%lsm_e(ile,jk,ibe) <= sea_boundary ) THEN
              cell1_cc%x  = patch_2D%cells%cartesian_center(icell_idx_1,icell_blk_1)%x
              cell2_cc%x  = patch_2D%cells%cartesian_center(icell_idx_2,icell_blk_2)%x

              !Check, if edge is sea or boundary edge and take care of dummy edge
              !edge with indices ile, ibe is sea edge
              !Add up for wet dual area.
              !IF ( v_base%lsm_e(ile,jk,ibe) <= sea_boundary ) THEN
              operators_coefficients%variable_dual_vol_norm(jv,jk,jb,jev)= triangle_area(cell1_cc, vertex_cc, cell2_cc)
              ! edge with indices ile, ibe is boundary edge
            ELSE IF ( patch_3D%lsm_e(ile,jk,ibe) == boundary ) THEN
              operators_coefficients%variable_dual_vol_norm(jv,jk,jb,jev)=0.0_wp!0.5_wp*triangle_area(cell1_cc, vertex_cc, cell2_cc)
            END IF
          END DO

          !---------------------------------------------------------------------------------------------
          DO je = 1, boundary_counter
            ibnd_edge_idx(je) = operators_coefficients%bnd_edge_idx(jv,jk,jb,je)
            ibnd_edge_blk(je) = operators_coefficients%bnd_edge_blk(jv,jk,jb,je)

            operators_coefficients%rot_coeff(jv,jk,jb,i_edge_idx(je) )=&
              & 0.5_wp*patch_2D%edges%system_orientation(ibnd_edge_idx(je),ibnd_edge_blk(je)) * &
              & patch_2D%edges%primal_edge_length(ibnd_edge_idx(je),ibnd_edge_blk(je))

          ENDDO
        END DO ! jv = i_startidx_v, i_endidx_v

        !!$OMP END PARALLEL DO

      END DO ! jk = 1, n_zlev
    END DO ! jb = owned_verts%start_block, owned_verts%end_block
    ! sync rot_coeff is done after area normalization at the end
    !-------------------------------------------------------------


    !-------------------------------------------------------------
    !The dynamical changing coefficient for the surface layer
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c

        DO je = 1, patch_2D%cells%num_edges(jc,jb)
          operators_coefficients%edge2cell_coeff_cc_dyn(jc,1,jb,je)%x = &
            operators_coefficients%edge2cell_coeff_cc(jc,1,jb,je)%x
        ENDDO 
      END DO ! jc = i_startidx_c, i_endidx_c
    END DO ! jb = all_cells%start_block, all_cells%end_block
    !-------------------------------------------------------------

    !-------------------------------------------------------------
    !3.3) Edge to vert coefficient
    DO jb = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, jb, i_startidx_v, i_endidx_v)
      DO jk = 1, n_zlev
        DO jv = i_startidx_v, i_endidx_v

          DO je = 1, patch_2D%verts%num_edges(jv,jb)
            ile = patch_2D%verts%edge_idx(jv,jb,je)
            ibe = patch_2D%verts%edge_blk(jv,jb,je)

             IF ( patch_3D%lsm_e(ile,jk,ibe) /= sea) THEN
              operators_coefficients%edge2vert_coeff_cc(jv,jk,jb,je)%x(1:3) = 0.0_wp
              operators_coefficients%variable_dual_vol_norm(jv,jk,jb,je)    = 0.0_wp
            ENDIF
          ENDDO ! je = 1, patch_2D%verts%num_edges(jv,jb)
        ENDDO ! jv = i_startidx_v, i_endidx_v
      END DO ! jk = 1, n_zlev
    END DO ! jb = owned_verts%start_block, owned_verts%end_block
    ! sync the result
    DO je=1,no_dual_edges
    ! these will be synced in the next loop
!       CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,:,:, je)%x(1))
!       CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,:,:, je)%x(2))
!       CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,:,:, je)%x(3))
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%variable_dual_vol_norm(:,:,:, je))
    ENDDO
    !-------------------------------------------------------------

    !-------------------------------------------------------------
    !Merge dual area calculation with coefficients
    ! note: sea_edges_per_vertex has been calculated on the owned_verts
    !       it does not need to be synced
    DO jb = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, jb, i_startidx_v, i_endidx_v)
      DO jk = 1, n_zlev
!CDIR nextscalar
        DO jv = i_startidx_v, i_endidx_v

          IF ( sea_edges_per_vertex(jv,jk,jb) == no_dual_edges ) THEN ! we have to count for lateral boundaries at the top
            zarea_fraction(jv,jk,jb)= patch_2D%verts%dual_area(jv,jb)/(earth_radius*earth_radius)
            !zarea_fraction(jv,jk,jb)=SUM(operators_coefficients%variable_dual_vol_norm(jv,jk,jb,:))

            !ELSEIF(operators_coefficients%bnd_edges_per_vertex(jv,jk,jb)/=0)THEN!boundary edges are involved
          ELSEIF ( sea_edges_per_vertex(jv,jk,jb) /= 0 ) THEN

            !Modified area calculation
            vertex_cc = patch_2D%verts%cartesian(jv,jb)
            DO jev = 1, patch_2D%verts%num_edges(jv,jb)
              ! get line and block indices of edge jev around vertex jv
              ile = patch_2D%verts%edge_idx(jv,jb,jev)
              ibe = patch_2D%verts%edge_blk(jv,jb,jev)
              !get neighbor cells
              icell_idx_1 = patch_2D%edges%cell_idx(ile,ibe,1)
              icell_idx_2 = patch_2D%edges%cell_idx(ile,ibe,2)
              icell_blk_1 = patch_2D%edges%cell_blk(ile,ibe,1)
              icell_blk_2 = patch_2D%edges%cell_blk(ile,ibe,2)

              !Check, if edge is sea or boundary edge and take care of dummy edge
              !edge with indices ile, ibe is sea edge
              !Add up for wet dual area.
              ! Note that this should be modified.
              !   sea_boundary means an open boundary
              !   boundary means that only the sea cell are should be added
              IF ( patch_3D%lsm_e(ile,jk,ibe) <= sea_boundary ) THEN
                cell1_cc%x  = patch_2D%cells%cartesian_center(icell_idx_1,icell_blk_1)%x
                cell2_cc%x  = patch_2D%cells%cartesian_center(icell_idx_2,icell_blk_2)%x
                zarea_fraction(jv,jk,jb) = zarea_fraction(jv,jk,jb)  &
                  & + triangle_area(cell1_cc, vertex_cc, cell2_cc)
                ! edge with indices ile, ibe is boundary edge                
              ELSE IF ( patch_3D%lsm_e(ile,jk,ibe) == boundary ) THEN
                ! at least one of the two cells exists and is sea cell
                IF (icell_idx_2 <= 0) THEN
                  cell1_cc%x  = patch_2D%cells%cartesian_center(icell_idx_1,icell_blk_1)%x
                ELSE IF (icell_idx_1 <= 0) THEN
                  cell1_cc%x  = patch_2D%cells%cartesian_center(icell_idx_2,icell_blk_2)%x
                ELSE IF (patch_3D%lsm_c(icell_idx_1,jk,icell_blk_1) <= sea_boundary) THEN
                  cell1_cc%x  = patch_2D%cells%cartesian_center(icell_idx_1,icell_blk_1)%x
                ELSE
                  cell1_cc%x  = patch_2D%cells%cartesian_center(icell_idx_2,icell_blk_2)%x
                ENDIF
!                 zarea_fraction(jv,jk,jb) = zarea_fraction(jv,jk,jb)  &
!                   & + 0.5_wp*triangle_area(cell1_cc, vertex_cc, cell2_cc)
                ! add only the sea dual area, ie the triagle area between
                ! the vertex, edge centre, and sea cell center
                zarea_fraction(jv,jk,jb) = zarea_fraction(jv,jk,jb)  &
                  & + triangle_area(cell1_cc, vertex_cc, patch_2D%edges%cartesian_center(ile,ibe))
              ENDIF
              
            END DO ! jev = 1, patch_2D%verts%num_edges(jv,jb)
            
          ENDIF !( sea_edges_per_vertex(jv,jk,jb) == patch_2D%verts%num_edges(jv,jb) )
          !The two quantities: 
          !zarea_fraction(jv,jk,jb)*(earth_radius*earth_radius) 
          !and 
          !patch_2D%verts%dual_area(jv,jb)
          !are identical
!CDIR nextscalar

          !Final coefficient calculation
          IF(zarea_fraction(jv,jk,jb)/=0.0_wp)THEN

            operators_coefficients%rot_coeff(jv,jk,jb,:)&
            &=operators_coefficients%rot_coeff(jv,jk,jb,:)/(zarea_fraction(jv,jk,jb)*(earth_radius*earth_radius))
            
            DO jev = 1, patch_2D%verts%num_edges(jv,jb)
              operators_coefficients%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)&
                & =operators_coefficients%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)/zarea_fraction(jv,jk,jb)
                !SUM(operators_coefficients%variable_dual_vol_norm(jv,jk,jb,:))!
            END DO

          ELSE
            DO jev = 1, patch_2D%verts%num_edges(jv,jb)
                operators_coefficients%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)=0.0_wp
            END DO
            operators_coefficients%rot_coeff(jv,jk,jb,:)=0.0_wp
          ENDIF
         !!ENDIF !( sea_edges_per_vertex(jv,jk,jb) == patch_2D%verts%num_edges(jv,jb) )

        ENDDO!jv = i_startidx_v, i_endidx_v
      END DO!jk = 1, n_zlev
    END DO!jb = owned_verts%start_block, owned_verts%end_block
    ! sync the result
    DO jev=1,no_dual_edges
      DO jk = 1, n_zlev
        CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,jk,:, jev)%x(1))
        CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,jk,:, jev)%x(2))
        CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,jk,:, jev)%x(3))
      ENDDO
    ENDDO
    DO je=1,patch_2D%cell_type
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%rot_coeff(:,:,:, je))
    ENDDO

    DO jk = 1, n_zlev
      CALL sync_patch_array(SYNC_V, patch_2D, zarea_fraction(:,jk,:))
    ENDDO
    DO jev=1,2*no_dual_edges
      DO jk = 1, n_zlev
        CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2edge_viavert_coeff(:,jk,:, jev))
      ENDDO
    ENDDO
    
    DO jb = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e
        DO jk = 1, n_zlev
          
          IF ( patch_3D%lsm_e(je,jk,jb) <= sea_boundary ) THEN

            DO neigbor=1,2

              jv  = patch_2D%edges%vertex_idx(je, jb, neigbor)
              jev = patch_2D%edges%vertex_blk(je, jb, neigbor)

              IF(neigbor==1)THEN
                operators_coefficients%edge2edge_viavert_coeff(je,jk,jb,1:no_dual_edges)&
                &=operators_coefficients%edge2edge_viavert_coeff(je,jk,jb,1:no_dual_edges)&
                &/zarea_fraction(jv,jk,jev)!SUM(operators_coefficients%variable_dual_vol_norm(jv,jk,jev,:))
              ELSEIF(neigbor==2)THEN
                operators_coefficients%edge2edge_viavert_coeff(je,jk,jb,no_dual_edges+1:2*no_dual_edges)&
                &=operators_coefficients%edge2edge_viavert_coeff(je,jk,jb,no_dual_edges+1:2*no_dual_edges)&
                &/zarea_fraction(jv,jk,jev)!SUM(operators_coefficients%variable_dual_vol_norm(jv,jk,jev,:))
              ENDIF
            
            END DO !neigbor=1,2

          ENDIF !  patch_3D%lsm_e(je,jk,jb) <= sea_boundary
          
        END DO
      END DO
    END DO

    DO jev=1,2*no_dual_edges
      DO jk = 1, n_zlev
        CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2edge_viavert_coeff(:,jk,:, jev))
      ENDDO
    ENDDO

    !-------------------------------------------------------------

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE apply_boundary2coeffs
!-------------------------------------------------------------------------
  !>
  !! Precomputes the geometrical factors used in the divergence, rotation.
  !!
  !! @par Revision History
  !!  developed by Guenther Zaengl, 2009-03-17
  !!  Modification by Almut Gassmann, 2009-12-19
  !!  - Vorticity is computed on quads in case of the hexagonal grid
  !!  Modification by Almut Gassmann, 2010-02-05
  !!  - Added feature for poor men's 3rd order advection, where a directional
  !!    laplace is needed at the edges.
  !!  Modification by Stephan Lorenz, 2010-06-02
  !!  - Storage moved from int_state into patch_oce since it is static
  !!    geometric information used in the ocean model
  !!  Modification by Peter Korn, 2010-11
  !!  - Calculation of cell area changed to achieve compatibility with
  !!    sw-model (cell area and consequently divergence different)
  !!  Modification by Stephan Lorenz, 2011-07
  !!   - 3-dim structures moved from patch_oce to hydro_ocean_base for parallelization
  !!
  SUBROUTINE init_diff_operator_coeff_3D( patch_2D, operators_coefficients )
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: patch_2D
    TYPE(t_operator_coeff),     INTENT(inout) :: operators_coefficients
    !

    INTEGER :: ie
    !INTEGER :: rl_start, rl_end
    INTEGER :: i_nchdom!,i_startblk, i_endblk, i_startidx, i_endidx

    !INTEGER :: ile, ibe!, ilc1, ibc1, ilc2, ibc2, ifac, ic, ilnc, ibnc
    !INTEGER :: ile1, ibe1,ile2,ibe2,ile3,ibe3
    !INTEGER, PARAMETER :: i_cell_type = 3
    !TYPE(cartesian_coordinates)::z_pn_k,z_pn_j
    !REAL(wp) :: z_lon, z_lat, z_nu, z_nv, z_proj
    !REAL(wp) :: cell_area
  
    TYPE(t_subset_range), POINTER :: all_edges, owned_cells, owned_verts
    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor, level

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:init_diff_operator_coeff_3D')
    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')
    i_nchdom   = MAX(1,patch_2D%n_childdom)
    all_edges => patch_2D%edges%all
    owned_cells => patch_2D%cells%owned
    owned_verts => patch_2D%verts%owned

!!$OMP PARALLEL
    ! 1) coefficients for divergence
!!$OMP DO PRIVATE(cell_block, start_index, end_index, cell_index, neigbor, &
!!$OMP edge_index, edge_block, level) ICON_OMP_DEFAULT_SCHEDULE
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index
        DO neigbor=1, patch_2D%cells%num_edges(cell_index, cell_block)

          edge_index = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
          edge_block = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

          operators_coefficients%div_coeff(cell_index, 1, cell_block, neigbor)               = &
            & patch_2D%edges%primal_edge_length(edge_index, edge_block)        * &
            & patch_2D%cells%edge_orientation(cell_index, cell_block, neigbor) / &
            & patch_2D%cells%area(cell_index, cell_block)

          DO level=2, n_zlev
            operators_coefficients%div_coeff(cell_index, level, cell_block, neigbor) = &
               operators_coefficients%div_coeff(cell_index, 1, cell_block, neigbor)
          ENDDO !levels

        ENDDO !neigbor
      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block
!!$OMP ENDDO NOWAIT

    ! 2) coefficients for curl
!!$OMP DO PRIVATE(vertex_block, start_index, end_index, vertex_index, neigbor, &
!!$OMP edge_index, edge_block, level) ICON_OMP_DEFAULT_SCHEDULE
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index
        DO neigbor=1, patch_2D%verts%num_edges(vertex_index, vertex_block)
          edge_index = patch_2D%verts%edge_idx(vertex_index, vertex_block, neigbor)
          edge_block = patch_2D%verts%edge_blk(vertex_index, vertex_block, neigbor)

          operators_coefficients%rot_coeff(vertex_index,1,vertex_block,neigbor)  =     &
            & patch_2D%edges%dual_edge_length(edge_index,edge_block) *     &
            & patch_2D%verts%edge_orientation(vertex_index,vertex_block,neigbor)

          DO level=2, n_zlev
            operators_coefficients%rot_coeff(vertex_index,level,vertex_block,neigbor) = &
               operators_coefficients%rot_coeff(vertex_index,1,vertex_block,neigbor)
          ENDDO !levels
        
        ENDDO !neigbor
      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
!!$OMP END DO NOWAIT


    ! 4) coefficients for gradient
!!$OMP DO PRIVATE(edge_block, start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index
        ! this should be calculated in the patch_2D setup
        patch_2D%edges%inv_dual_edge_length(edge_index,edge_block) = &
          & 1._wp / patch_2D%edges%dual_edge_length(edge_index,edge_block)

        DO level=1, n_zlev
          operators_coefficients%grad_coeff(edge_index,level, edge_block)   &
            & = patch_2D%edges%inv_dual_edge_length(edge_index,edge_block)
        ENDDO !levels

      ENDDO ! edge_index = start_index, end_index
    ENDDO  ! edge_block
!!$OMP ENDDO NOWAIT
!!$OMP END PARALLEL

    ! no need to synchronize all elements of operators_coefficients%grad_coeff
!     CALL sync_patch_array(sync_e, patch_2D, operators_coefficients%grad_coeff)

    DO ie = 1, no_primal_edges
      CALL sync_patch_array(sync_c, patch_2D, operators_coefficients%div_coeff(:,:,:,ie))
    END DO

    DO ie = 1, no_dual_edges!9-i_cell_type
      CALL sync_patch_array(sync_v, patch_2D, operators_coefficients%rot_coeff(:,:,:,ie))
    END DO

    CALL message (TRIM(routine), 'end')

    ! 3) coefficients for nabla2_scalar
    !rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    !rl_end = min_rlcell_int
    ! values for the blocking
    !i_startblk = patch_2D%cells%start_blk(rl_start,1)
    !i_endblk   = patch_2D%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch_2D cells (and blocks)
    !
    !!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
    !!$OMP    ilc2,ibc2,ilnc,ibnc)
    ! !       DO jk=1,n_zlev
    ! !       DO jb = i_startblk, i_endblk
    ! !
    ! !         CALL get_indices_c(patch_2D, jb, i_startblk, i_endblk, &
    ! !                            i_startidx, i_endidx, rl_start, rl_end)
    ! !
    ! !         DO je = 1, i_cell_type
    ! !           DO jc = i_startidx, i_endidx
    ! !
    ! !             ile = patch_2D%cells%edge_idx(jc,jb,je)
    ! !             ibe = patch_2D%cells%edge_blk(jc,jb,je)
    ! !
    ! !             ilc1 = patch_2D%edges%cell_idx(ile,ibe,1)
    ! !             ibc1 = patch_2D%edges%cell_blk(ile,ibe,1)
    ! !             ilc2 = patch_2D%edges%cell_idx(ile,ibe,2)
    ! !             ibc2 = patch_2D%edges%cell_blk(ile,ibe,2)
    ! !
    ! !             IF (jc == ilc1 .AND. jb == ibc1) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 operators_coefficients%geofac_n2s(jc,1,jb)     =  &
    ! !                 &  operators_coefficients%geofac_n2s(jc,1,jb)  -  &
    ! !                 &  operators_coefficients%geofac_div(jc,je,jb) /  &
    ! !                 &  patch_2D%edges%dual_edge_length(ile,ibe)
    ! !             ENDIF
    ! !           ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
    ! !             IF (i_cell_type == 3) THEN
    ! !               operators_coefficients%geofac_n2s(jc,1,jb)       =  &
    ! !                 &  operators_coefficients%geofac_n2s(jc,1,jb)  +  &
    ! !                 &  operators_coefficients%geofac_div(jc,je,jb) /  &
    ! !                 &  patch_2D%edges%dual_edge_length(ile,ibe)
    ! !             ENDIF
    ! !           ENDIF
    ! !           DO ic = 1, i_cell_type
    ! !             ilnc = patch_2D%cells%neighbor_idx(jc,jb,ic)
    ! !             ibnc = patch_2D%cells%neighbor_blk(jc,jb,ic)
    ! !             IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 operators_coefficients%geofac_n2s(jc,ic+1,jb)     = &
    ! !                   &  operators_coefficients%geofac_n2s(jc,ic+1,jb)- &
    ! !                   &  operators_coefficients%geofac_div(jc,je,jb)  / &
    ! !                   &  patch_2D%edges%dual_edge_length(ile,ibe)
    ! !               ENDIF
    ! !             ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 operators_coefficients%geofac_n2s(jc,ic+1,jb)     = &
    ! !                   &  operators_coefficients%geofac_n2s(jc,ic+1,jb)+ &
    ! !                   &  operators_coefficients%geofac_div(jc,je,jb)  / &
    ! !                   &  patch_2D%edges%dual_edge_length(ile,ibe)
    ! !               ENDIF
    ! !             ENDIF
    ! !           ENDDO
    ! !           ! To ensure that dummy edges have a factor of 0:
    ! !           IF (je > patch_2D%cells%num_edges(jc,jb)) THEN
    ! !             operators_coefficients%geofac_n2s(jc,je+1,jb) = 0._wp
    ! !           ENDIF
    ! !
    ! !         ENDDO !cell loop
    ! !       ENDDO
    ! !     END DO !block loop
    ! !     END DO
    !!$OMP END DO

  END SUBROUTINE init_diff_operator_coeff_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes area of triangular cell.
  !!
  !! @par Revision History
  !! Developed  by Luis Kornblueh  (2004).
  ELEMENTAL FUNCTION triangle_area (x0, x1, x2) result(area)

    TYPE(t_cartesian_coordinates), INTENT(in) :: x0, x1, x2

    REAL(wp) :: area
    REAL(wp) :: z_s12, z_s23, z_s31, z_ca1, z_ca2, z_ca3, z_a1, z_a2, z_a3

    TYPE(t_cartesian_coordinates) :: u12, u23, u31

    ! This variant to calculate the area of a spherical triangle
    ! is more precise.

    !  Compute cross products Uij = Vi x Vj.
    u12 = vector_product (x0, x1)
    u23 = vector_product (x1, x2)
    u31 = vector_product (x2, x0)

    !  Normalize Uij to unit vectors.
    z_s12 = DOT_PRODUCT ( u12%x(1:3), u12%x(1:3) )
    z_s23 = DOT_PRODUCT ( u23%x(1:3), u23%x(1:3) )
    z_s31 = DOT_PRODUCT ( u31%x(1:3), u31%x(1:3) )

    !  Test for a degenerate triangle associated with collinear vertices.
    IF ( z_s12 == 0.0_wp .or. z_s23 == 0.0_wp  .or. z_s31 == 0.0_wp ) THEN
      area = 0.0_wp
      RETURN
    END IF

    z_s12 = SQRT(z_s12)
    z_s23 = SQRT(z_s23)
    z_s31 = SQRT(z_s31)

    u12%x(1:3) = u12%x(1:3)/z_s12
    u23%x(1:3) = u23%x(1:3)/z_s23
    u31%x(1:3) = u31%x(1:3)/z_s31

    !  Compute interior angles Ai as the dihedral angles between planes:
    !  CA1 = cos(A1) = -<U12,U31>
    !  CA2 = cos(A2) = -<U23,U12>
    !  CA3 = cos(A3) = -<U31,U23>
    z_ca1 = -u12%x(1)*u31%x(1)-u12%x(2)*u31%x(2)-u12%x(3)*u31%x(3)
    z_ca2 = -u23%x(1)*u12%x(1)-u23%x(2)*u12%x(2)-u23%x(3)*u12%x(3)
    z_ca3 = -u31%x(1)*u23%x(1)-u31%x(2)*u23%x(2)-u31%x(3)*u23%x(3)

    IF (z_ca1 < -1.0_wp) z_ca1 = -1.0_wp
    IF (z_ca1 >  1.0_wp) z_ca1 =  1.0_wp
    IF (z_ca2 < -1.0_wp) z_ca2 = -1.0_wp
    IF (z_ca2 >  1.0_wp) z_ca2 =  1.0_wp
    IF (z_ca3 < -1.0_wp) z_ca3 = -1.0_wp
    IF (z_ca3 >  1.0_wp) z_ca3 =  1.0_wp

    z_a1 = ACOS(z_ca1)
    z_a2 = ACOS(z_ca2)
    z_a3 = ACOS(z_ca3)

    !  Compute areas = z_a1 + z_a2 + z_a3 - pi.
    area = z_a1+z_a2+z_a3-pi

    IF ( area < 0.0_wp ) area = 0.0_wp

  END FUNCTION triangle_area
 !-------------------------------------------------------------------------
!   !-------------------------------------------------------------------------
! 
!   !> Initialize 3D expansion coefficients.
!   !!
!   !! @par Revision History
!   !! Peter Korn (2012-2)
!   !! Parellelized by Leonidas Linardakis 2012-3
! 
! 
!   SUBROUTINE copy_2D_to_3D_coeff( patch_2D,                 &
!                                & operators_coefficients,            &
!                                & edge2edge_viacell_coeff,&
!                                & edge2edge_viavert_coeff,&
!                                & edge2cell_coeff_cc,     &
!                                & edge2cell_coeff_cc_t,   &
!                                & edge2vert_coeff_cc,     &
!                                & edge2vert_coeff_cc_t,   &
!                                & dist_cell2edge,         &
!                                & fixed_vol_norm,         &
!                                & variable_vol_norm,      &
!                                & variable_dual_vol_norm)
!     TYPE(t_patch)    , TARGET, INTENT(INOUT)     :: patch_2D
!     TYPE(t_operator_coeff), INTENT(inout) :: operators_coefficients
!     REAL(wp), INTENT(IN)                      :: edge2edge_viacell_coeff(1:nproma,1:patch_2D%nblks_e,1:2*no_primal_edges)
!     REAL(wp), INTENT(IN)                      :: edge2edge_viavert_coeff(1:nproma,1:patch_2D%nblks_e,1:2*no_dual_edges)
!     TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2cell_coeff_cc     (1:nproma,1:patch_2D%nblks_c,1:patch_2D%cell_type)
!     TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2cell_coeff_cc_t   (1:nproma,1:patch_2D%nblks_e,1:2)
!     TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2vert_coeff_cc     (1:nproma,1:patch_2D%nblks_v,1:no_dual_edges)
!     TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2vert_coeff_cc_t   (1:nproma,1:patch_2D%nblks_e,1:2)
!     REAL(wp), INTENT(IN) :: dist_cell2edge        (1:nproma,1:patch_2D%nblks_e,1:2)
!     REAL(wp), INTENT(IN) :: fixed_vol_norm        (1:nproma,patch_2D%nblks_c)
!     REAL(wp), INTENT(IN) :: variable_vol_norm     (1:nproma,1:patch_2D%nblks_c,1:no_primal_edges)
!     REAL(wp), INTENT(IN) :: variable_dual_vol_norm(1:nproma,1:patch_2D%nblks_v,1:no_dual_edges)
!     !
!     !Local variables
!     !
!     TYPE(t_subset_range), POINTER :: all_edges
!     TYPE(t_subset_range), POINTER :: all_cells
!     TYPE(t_subset_range), POINTER :: all_verts
! 
!     INTEGER :: edge_block, cell_block, vertex_block, level, neigbor
!     INTEGER :: je,jk, i_startidx_e, i_endidx_e
! 
!     all_cells => patch_2D%cells%all
!     all_edges => patch_2D%edges%all
!     all_verts => patch_2D%verts%all
! 
!     !---------------------------------------------------------
!     ! the following coefficients will be copied:
!     !
!     ! operators_coefficients%edge_position_cc(:,:,:)            on edges
!     ! operators_coefficients%dist_cell2edge(:,:,:,1-2)          on edges
!     !
!     ! operators_coefficients%edge2cell_coeff_cc(:,:,:,1-3)%x    on cells
!     ! operators_coefficients%fixed_vol_norm(:,:,:)              on cells
!     ! operators_coefficients%variable_vol_norm(:,:,:,1-3)       on cells
!     !
!     ! operators_coefficients%edge2cell_coeff_cc_t(:,:,:,1-2)%x  on edges
!     ! operators_coefficients%edge2vert_coeff_cc(:,:,:,1-6)%x   on verts
!     ! operators_coefficients%edge2vert_coeff_cc_t(:,:,:,1-2)%x  on edges
!     !
!     ! patch_2D%edges%f_e(:, :) is already calculated in par_init_scalar_product_oce    !
!     !---------------------------------------------------------
! 
! 
!     !---------------------------------------------------------
!     ! calculate operators_coefficients%edge_position_cc(:,:,:) on edges
!     ! this is the same as the 2D, just copy it
!     ! it does not change, maybe turn it to 2D?
!     DO edge_block = all_edges%start_block, all_edges%end_block
!       DO level = 1, n_zlev
!         operators_coefficients%edge_position_cc(:,level,edge_block) = &
!           patch_2D%edges%cartesian_center(:,edge_block)
!       ENDDO
!     ENDDO
!     !---------------------------------------------------------
! 
!     !---------------------------------------------------------
!     ! calculate operators_coefficients%dist_cell2edge(:,:,:,1-2) on edges
!     ! this is the same as the 2D, just copy it
!     ! it does not change, maybe turn it to 2D?
!     DO edge_block = all_edges%start_block, all_edges%end_block
!       DO level = 1, n_zlev
!         operators_coefficients%dist_cell2edge(:,level,edge_block,1) = &
!           dist_cell2edge(:,edge_block,1)
!         operators_coefficients%dist_cell2edge(:,level,edge_block,2) = &
!           dist_cell2edge(:,edge_block,2)
!       ENDDO
!     ENDDO
!     !---------------------------------------------------------
! 
!     !---------------------------------------------------------
!     ! For the rest of coefficients,
!     ! copy the 2D and then apply_boundary2coeffs
!     !
!     ! on edges
!     DO edge_block = all_edges%start_block, all_edges%end_block
!       DO level = 1, n_zlev
!         CALL get_index_range(all_edges, edge_block, i_startidx_e, i_endidx_e)
!         DO je =  i_startidx_e, i_endidx_e
! 
!           operators_coefficients%edge2cell_coeff_cc_t(je,level,edge_block,1)%x = &
!             edge2cell_coeff_cc_t(je,edge_block,1)%x
!           operators_coefficients%edge2cell_coeff_cc_t(je,level,edge_block,2)%x = &
!             edge2cell_coeff_cc_t(je,edge_block,2)%x
! 
!           operators_coefficients%edge2vert_coeff_cc_t(je,level,edge_block,1)%x = &
!             edge2vert_coeff_cc_t(je,edge_block,1)%x
!           operators_coefficients%edge2vert_coeff_cc_t(je,level,edge_block,2)%x = &
!             edge2vert_coeff_cc_t(je,edge_block,2)%x
! 
!           !operators_coefficients%edge2edge_viacell_coeff(je,level,edge_block,1:6)=&
!           !&edge2edge_viacell_coeff(je,edge_block,1:6)
! 
!           !operators_coefficients%edge2edge_viavert_coeff(je,level,edge_block,1:12)=&
!           !&edge2edge_viavert_coeff(je,edge_block,1:12)
! 
!         ENDDO
!       ENDDO
!     ENDDO
! 
!     ! on cells
!     DO cell_block = all_cells%start_block, all_cells%end_block
!       DO level = 1, n_zlev
! 
!        operators_coefficients%fixed_vol_norm(:,level,cell_block) = &
!          fixed_vol_norm(:,cell_block)
! 
!        DO neigbor=1,no_primal_edges!patch_2D%cell_type
! 
!          operators_coefficients%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(1) = &
!            & edge2cell_coeff_cc(:,cell_block,neigbor)%x(1)
!          operators_coefficients%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(2) = &
!            & edge2cell_coeff_cc(:,cell_block,neigbor)%x(2)
!          operators_coefficients%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(3) = &
!            & edge2cell_coeff_cc(:,cell_block,neigbor)%x(3)
! 
!          operators_coefficients%variable_vol_norm(:,level,cell_block,neigbor) = &
!            & variable_vol_norm(:,cell_block,neigbor)
! 
!         ENDDO ! neigbor=1,patch_2D%cell_type
! 
!       ENDDO  !  level = 1, n_zlev
!     ENDDO ! cell_block
! 
!     ! on verts
!     DO vertex_block = all_verts%start_block, all_verts%end_block
!       DO level = 1, n_zlev
!         DO neigbor=1,no_dual_edges
! 
!           operators_coefficients%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(1) &
!              &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(1)
!           operators_coefficients%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(2) &
!             &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(2)
!           operators_coefficients%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(3) &
!             &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(3)
! 
!          operators_coefficients%variable_dual_vol_norm(:,level,vertex_block,neigbor)&
!          &=variable_dual_vol_norm(:,vertex_block,neigbor)
!         ENDDO ! neigbor=1,patch_2D%cell_type
!       ENDDO  !  level = 1, n_zlev
!     ENDDO ! vertex_block
! 
!    !DO neigbor=1,2*no_dual_edges
!    !  CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2edge_viavert_coeff(:,:,:,neigbor))
!    !END DO
! 
!   END SUBROUTINE copy_2D_to_3D_coeff
!   !-------------------------------------------------------------------------
!     !-------------------------------------------------------------------------
!   !>
!   !! Computes the coefficients that determine the scalar product on the primal grid. This
!   !! scalar product depends on the grid geometry only and  is used to formulate the primitive
!   !! equations in weak form. The coefficients are applied in module "mo_scalar_product".
!   !! The following components of the data type "ocean_patch" are filled:
!   !!   edge2cell_coeff  : coefficients for edge to cell mapping
!   !!   edge2cell_coeff_t: coefficients for transposed of edge to cell mappings
!   !!   edge2vert_coeff  : coefficients for edge to vertex mapping
!   !!   edge2vert_coeff_t: coefficients for transposed of edge to vertex mappings
!   !!   fixed_vol_norm   : summed volume weight of moved cell
!   !!   variable_vol_norm: volume weight at the edges of moved cell
!   !!
!   !! @par Revision History
!   !!  developed by Peter Korn, MPI-M  2010-09
!   !!  Modification by Stephan Lorenz, 2010-11
!   !!
!   SUBROUTINE par_init_coeff_2D( patch_2D,               &
!                                         & edge2cell_coeff_cc,   &
!                                         & edge2cell_coeff_cc_t, &
!                                         & edge2vert_coeff_cc,   &
!                                         & edge2vert_coeff_cc_t, &
!                                         & dist_cell2edge,       &
!                                         & fixed_vol_norm,       &
!                                         & variable_vol_norm,    &
!                                         & variable_dual_vol_norm)
!     TYPE(t_patch)    , TARGET, INTENT(INOUT)     :: patch_2D
!     TYPE(t_cartesian_coordinates), INTENT(INOUT) :: edge2cell_coeff_cc(1:nproma,1:patch_2D%nblks_c,1:no_primal_edges)
!     TYPE(t_cartesian_coordinates), INTENT(INOUT) :: edge2cell_coeff_cc_t(1:nproma,1:patch_2D%nblks_e,1:2)
!     TYPE(t_cartesian_coordinates), INTENT(INOUT) :: edge2vert_coeff_cc(1:nproma,1:patch_2D%nblks_v,1:no_dual_edges)
!     TYPE(t_cartesian_coordinates), INTENT(INOUT) :: edge2vert_coeff_cc_t(1:nproma,1:patch_2D%nblks_e,1:2)
!     REAL(wp), INTENT(INOUT) :: dist_cell2edge(1:nproma,1:patch_2D%nblks_e,1:2)
!     REAL(wp), INTENT(INOUT) :: fixed_vol_norm(1:nproma,patch_2D%nblks_c)
!     REAL(wp), INTENT(INOUT) :: variable_vol_norm(1:nproma,1:patch_2D%nblks_c,1:no_primal_edges)
!     REAL(wp), INTENT(INOUT) :: variable_dual_vol_norm(1:nproma,1:patch_2D%nblks_e,1:no_dual_edges)
! !
! !Local variables
! !
!     REAL(wp), ALLOCATABLE :: prime_edge_length( :, : )
!     REAL(wp), ALLOCATABLE :: dual_edge_length ( :, : )
! !     REAL(wp), ALLOCATABLE :: cell_area( :, : )
! !     REAL(wp), ALLOCATABLE :: dual_cell_area ( :, : )
! 
!     TYPE(t_subset_range), POINTER :: owned_edges         ! these are the owned entities
!     TYPE(t_subset_range), POINTER :: owned_cells         ! these are the owned entities
!     TYPE(t_subset_range), POINTER :: owned_verts         ! these are the owned entities
!     TYPE(t_cartesian_coordinates) :: vertex_position, cell_center, edge_center
!     TYPE(t_cartesian_coordinates) :: dist_vector
!     TYPE(t_cartesian_coordinates), POINTER :: dual_edge_middle(:,:)
! 
!     TYPE(t_cartesian_coordinates) :: coriolis_cartesian_coordinates
!     TYPE(t_geographical_coordinates) :: coriolis_geo_coordinates, geo_coordinates
!     REAL(wp) :: basin_center_lat_rad, basin_height_rad
! 
!     REAL(wp) :: norm, orientation, length
!     REAL(wp) :: inverse_sphere_radius
! 
!     INTEGER :: edge_block, edge_index
!     INTEGER :: cell_index, cell_block
!     INTEGER :: vertex_index, vertex_block
!     INTEGER :: start_index, end_index, neigbor
!     INTEGER :: cell_1_index, cell_1_block, cell_2_index, cell_2_block
!     INTEGER :: vertex_1_index, vertex_1_block, vertex_2_index, vertex_2_block
!     !-----------------------------------------------------------------------
! !     REAL(wp) :: dist_cell2edge(nproma, patch_2D%nblks_e,2)
! !     TYPE(t_cartesian_coordinates) :: check_v1(nproma, patch_2D%nblks_v, 6)
! !     TYPE(t_cartesian_coordinates) :: check_v2(nproma, patch_2D%nblks_e, 2)
! !     REAL(wp) :: max_diff, max_val
!     !-----------------------------------------------------------------------
!     inverse_sphere_radius = 1.0_wp / grid_sphere_radius
! 
!     owned_edges => patch_2D%edges%owned
!     owned_cells => patch_2D%cells%owned
!     owned_verts => patch_2D%verts%owned
! 
!     edge2vert_coeff_cc(:,:,:)%x(1) = 0.0_wp
!     edge2vert_coeff_cc(:,:,:)%x(2) = 0.0_wp
!     edge2vert_coeff_cc(:,:,:)%x(3) = 0.0_wp
! 
!     edge2vert_coeff_cc_t(:,:,:)%x(1) = 0.0_wp
!     edge2vert_coeff_cc_t(:,:,:)%x(2) = 0.0_wp
!     edge2vert_coeff_cc_t(:,:,:)%x(3) = 0.0_wp
! 
! 
!     !-------------------------------------------
!     ! compute some basic distances
!     ! this is required if the cartesian distance is used
!     ! instead of the spherical
!     !
!     ! computes_dist_cell2edge( patch_2D, intp_2D_coeff)
!     !
!     ALLOCATE( prime_edge_length( nproma, patch_2D%nblks_e))
!     ALLOCATE( dual_edge_length ( nproma, patch_2D%nblks_e))
! !     ALLOCATE( cell_area        ( nproma, patch_2D%nblks_c))
! !     ALLOCATE( dual_cell_area   ( nproma, patch_2D%nblks_v))
! 
!     IF ( MID_POINT_DUAL_EDGE ) THEN
!       dual_edge_middle => patch_2D%edges%cartesian_dual_middle
!     ELSE
!       dual_edge_middle => patch_2D%edges%cartesian_center
!     ENDIF
! 
!     ! get the areas on a unit sphere
! !     cell_area(:,:)      = patch_2D%cells%area(:,:)      * inverse_earth_radius * inverse_earth_radius
! !     dual_cell_area(:,:) = patch_2D%verts%dual_area(:,:) * inverse_earth_radius * inverse_earth_radius
! 
!     IF (LARC_LENGTH) THEN
! 
!       ! we just need to get them from the grid
!       ! NOTE:  these are earth's distances, translate on a unit sphere
!       dist_cell2edge(:,:,:) = &
!         & patch_2D%edges%edge_cell_length(:,:,:) * inverse_sphere_radius
!       prime_edge_length(:,:) = &
!         & patch_2D%edges%primal_edge_length(:,:) * inverse_sphere_radius
!       dual_edge_length(:,:) = &
!         & patch_2D%edges%dual_edge_length(:,:) * inverse_sphere_radius
! 
!     ELSE
! 
!       ! calcultate cartesian distance
!       prime_edge_length(:,:) = 0.0_wp
!       dual_edge_length(:,:) = 0.0_wp
!       dist_cell2edge(:,:,:) =  0.0_wp
! 
!       DO edge_block = owned_edges%start_block, owned_edges%end_block
!         CALL get_index_range(owned_edges, edge_block, start_index, end_index)
!         DO edge_index = start_index, end_index
! 
!           !----------------------------------------
!           ! calculate the cartesian edge length
!           vertex_1_index = patch_2D%edges%vertex_idx(edge_index, edge_block, 1)
!           vertex_1_block = patch_2D%edges%vertex_blk(edge_index, edge_block, 1)
!           vertex_2_index = patch_2D%edges%vertex_idx(edge_index, edge_block, 2)
!           vertex_2_block = patch_2D%edges%vertex_blk(edge_index, edge_block, 2)
! 
!           dist_vector%x = &
!             & patch_2D%verts%cartesian(vertex_1_index, vertex_1_block)%x - &
!             & patch_2D%verts%cartesian(vertex_2_index, vertex_2_block)%x
! 
!             prime_edge_length(edge_index,edge_block) = &
!               & SQRT(SUM((  dist_vector%x *  dist_vector%x)))
!           !----------------------------------------
! 
!           !----------------------------------------
!           ! calculate the cartesian distance of the edge center to the cell center
!           DO neigbor = 1,2
! 
!             dist_cell2edge(edge_index,edge_block,neigbor) = 0.0_wp
! 
!             cell_index = patch_2D%edges%cell_idx(edge_index,edge_block,neigbor)
!             cell_block = patch_2D%edges%cell_blk(edge_index,edge_block,neigbor)
! 
!             IF (cell_block > 0) THEN
!               dist_vector%x = &
!                 & patch_2D%edges%cartesian_center(edge_index,edge_block)%x - &
!                 & patch_2D%cells%cartesian_center(cell_index,cell_block)%x
! 
!               dist_cell2edge(edge_index,edge_block,neigbor) = &
!                 & SQRT(SUM((  dist_vector%x *  dist_vector%x)))
!             ENDIF
! 
!           ENDDO ! neigbor = 1,2
!           !----------------------------------------
! 
!           !----------------------------------------
!           ! calculate the cartesian dual edge length
!           cell_1_index = patch_2D%edges%cell_idx(edge_index, edge_block, 1)
!           cell_1_block = patch_2D%edges%cell_blk(edge_index, edge_block, 1)
!           cell_2_index = patch_2D%edges%cell_idx(edge_index, edge_block, 2)
!           cell_2_block = patch_2D%edges%cell_blk(edge_index, edge_block, 2)
! 
!           IF (cell_1_block > 0 .AND. cell_2_block > 0) THEN
!             dist_vector%x = &
!               & patch_2D%cells%cartesian_center(cell_1_index, cell_1_block)%x - &
!               & patch_2D%cells%cartesian_center(cell_2_index, cell_2_block)%x
! 
!               dual_edge_length(edge_index,edge_block) = &
!                 & SQRT(SUM((  dist_vector%x *  dist_vector%x)))
!           ELSE
!               dual_edge_length(edge_index,edge_block) =              &
!                 & dist_cell2edge(edge_index,edge_block,1) + &
!                 & dist_cell2edge(edge_index,edge_block,2)
!            ENDIF
!           !----------------------------------------
! 
!         ENDDO ! edge_index=start_index,end_index
!       ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
! 
!       ! synchronize the edge distances
!       CALL sync_patch_array(SYNC_E, patch_2D, dist_cell2edge(:,:,1))
!       CALL sync_patch_array(SYNC_E, patch_2D, dist_cell2edge(:,:,2))
!       CALL sync_patch_array(SYNC_E, patch_2D, prime_edge_length(:,:))
!       CALL sync_patch_array(SYNC_E, patch_2D, dual_edge_length(:,:))
!     ENDIF
!     ! distances have been computed
!     !-------------------------------------------
! 
!     !-------------------------------------------
!     ! compute:
!     !   edge2cell_coeff_cc
!     !   fixed_vol_norm
!     !   variable_vol_norm
!     edge2cell_coeff_cc(:,:,:)%x(1) = 0.0_wp
!     edge2cell_coeff_cc(:,:,:)%x(2) = 0.0_wp
!     edge2cell_coeff_cc(:,:,:)%x(3) = 0.0_wp
! 
!     edge2cell_coeff_cc_t(:,:,:)%x(1) = 0.0_wp
!     edge2cell_coeff_cc_t(:,:,:)%x(2) = 0.0_wp
!     edge2cell_coeff_cc_t(:,:,:)%x(3) = 0.0_wp
! 
!     fixed_vol_norm(:,:)       = 0.0_wp
!     variable_vol_norm(:,:,:)  = 0.0_wp
!     DO cell_block = owned_cells%start_block, owned_cells%end_block
!       CALL get_index_range(owned_cells, cell_block, start_index, end_index)
!       DO cell_index = start_index, end_index
! 
!         cell_center%x = patch_2D%cells%cartesian_center(cell_index, cell_block)%x
!         fixed_vol_norm(cell_index,cell_block) = 0.0_wp
! 
!         !-------------------------------
!         DO neigbor=1,patch_2D%cells%num_edges(cell_index,cell_block)  !no_primal_edges!patch_2D%cell_type
! 
!           edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x = 0.0_wp
!           variable_vol_norm(cell_index, cell_block, neigbor) =  0.0_wp
! 
!           edge_index = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
!           edge_block = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)
! 
!           IF (edge_block > 0 ) THEN ! this if should not be necessary but for safety let it be
!             ! we have an edge
!             dist_vector%x = &
!               & patch_2D%edges%cartesian_center(edge_index,edge_block)%x - &
!               & cell_center%x
! 
!             norm  = SQRT(SUM( dist_vector%x * dist_vector%x))
! 
!             ! compute edge2cell_coeff_cc
!             edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x =  &
!               & dist_vector%x *                                             &
!               & prime_edge_length(edge_index,edge_block) *                  &
!               & patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)! / &
!               ! & cell_area(cell_index, cell_block)
!               ! Note: here we do not divide by the cell area !
! 
!             fixed_vol_norm(cell_index,cell_block) = &
!               & fixed_vol_norm(cell_index,cell_block) + &
!               & 0.5_wp * norm * prime_edge_length(edge_index,edge_block)
! 
!             variable_vol_norm(cell_index, cell_block, neigbor) = &
!               & 0.5_wp * norm * prime_edge_length(edge_index,edge_block)
! 
!           ENDIF !(edge_block > 0 )
! 
!         ENDDO !neigbor=1,patch_2D%cell_type
!         !-------------------------------
! 
!       ENDDO ! cell_index = start_index, end_index
!     ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
!     !-------------------
!     ! sync the results
!     CALL sync_patch_array(SYNC_C, patch_2D, fixed_vol_norm(:,:))
!     DO neigbor=1,patch_2D%cell_type
!       CALL sync_patch_array(SYNC_C, patch_2D, edge2cell_coeff_cc(:,:,neigbor)%x(1))
!       CALL sync_patch_array(SYNC_C, patch_2D, edge2cell_coeff_cc(:,:,neigbor)%x(2))
!       CALL sync_patch_array(SYNC_C, patch_2D, edge2cell_coeff_cc(:,:,neigbor)%x(3))
!       CALL sync_patch_array(SYNC_C, patch_2D, variable_vol_norm(:,:,neigbor))
!     ENDDO
!     !-------------------
! 
!     !-------------------------------------------
!     ! compute:
!     !   edge2cell_coeff_cc_t
! 
!     DO edge_block = owned_edges%start_block, owned_edges%end_block
!       CALL get_index_range(owned_edges, edge_block, start_index, end_index)
!       DO edge_index = start_index, end_index
! 
!         edge2cell_coeff_cc_t(edge_index, edge_block, 2)%x = 0.0_wp
!         edge_center%x = patch_2D%edges%cartesian_center(edge_index, edge_block)%x
! 
!         DO neigbor=1,2
! 
!           edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x = 0.0_wp
!           cell_index = patch_2D%edges%cell_idx(edge_index, edge_block, neigbor)
!           cell_block = patch_2D%edges%cell_blk(edge_index, edge_block, neigbor)
! 
!           IF (cell_block > 0) THEN
! 
!             dist_vector%x =  edge_center%x -                             &
!               patch_2D%cells%cartesian_center(cell_index, cell_block)%x
! 
!             orientation = DOT_PRODUCT(dist_vector%x, &
!               & patch_2D%edges%primal_cart_normal(edge_index, edge_block)%x)
!             IF (orientation < 0.0_wp) dist_vector%x = - dist_vector%x
! 
!             edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x = &
!               dist_vector%x / dual_edge_length(edge_index, edge_block)
! 
!           ENDIF ! (cell_block > 0)
! 
!         ENDDO ! neigbor=1,2
! 
!       ENDDO ! edge_index = start_index, end_index
!     ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
!     !-------------------
!     ! sync the results
!     DO neigbor=1,2
!       CALL sync_patch_array(SYNC_E, patch_2D, edge2cell_coeff_cc_t(:,:,neigbor)%x(1))
!       CALL sync_patch_array(SYNC_E, patch_2D, edge2cell_coeff_cc_t(:,:,neigbor)%x(2))
!       CALL sync_patch_array(SYNC_E, patch_2D, edge2cell_coeff_cc_t(:,:,neigbor)%x(3))
!     ENDDO ! neigbor=1,2
!     !   edge2cell_coeff_cc_t is computed
!     !-------------------------------------------
! 
!     !-------------------------------------------
!     ! compute:
!     !   edge2vert_coeff_cc
!     !   variable_dual_vol_norm
!     variable_dual_vol_norm(:,:,:)  = 0.0_wp
!     DO vertex_block = owned_verts%start_block, owned_verts%end_block
!       CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
!       DO vertex_index = start_index, end_index
! 
!         vertex_position%x = patch_2D%verts%cartesian(vertex_index, vertex_block)%x
! 
!         DO neigbor=1, patch_2D%verts%num_edges(vertex_index,vertex_block) !no_dual_edges
!! we have to change this to accomodate the dual grid
! 
!           variable_dual_vol_norm(vertex_index, vertex_block, neigbor) = 0.0_wp
! 
!           edge_index = patch_2D%verts%edge_idx(vertex_index, vertex_block, neigbor)
!           edge_block = patch_2D%verts%edge_blk(vertex_index, vertex_block, neigbor)
! 
!           IF (edge_block > 0) THEN
!             ! we got an adjacent edge
!             dist_vector%x = &
!               dual_edge_middle(edge_index, edge_block)%x - &
!               vertex_position%x
! 
! 
!             ! the dist_vector has cartesian length
!             ! if we use spherical distance we need to recalculate
!             ! its length
!             IF (LARC_LENGTH) THEN
!               length = arc_length(vertex_position, dual_edge_middle(edge_index, edge_block))
!               norm = SQRT(SUM( dist_vector%x * dist_vector%x ))
!               dist_vector%x = dist_vector%x * length / norm
!             ELSE
!               length = SQRT(SUM( dist_vector%x * dist_vector%x ))
!             ENDIF
! 
!             dist_vector = vector_product(dist_vector, dual_edge_middle(edge_index, edge_block))
!             orientation = DOT_PRODUCT( dist_vector%x,                         &
!                & patch_2D%edges%primal_cart_normal(edge_index, edge_block)%x)
!             IF (orientation < 0.0_wp) dist_vector%x = - dist_vector%x
! 
!               edge2vert_coeff_cc(vertex_index, vertex_block, neigbor)%x = &
!               & dist_vector%x                                *                    &
!               & dual_edge_length(edge_index, edge_block) !    /                    &
!               !& dual_cell_area(vertex_index, vertex_block)
! 
!               variable_dual_vol_norm(vertex_index, vertex_block, neigbor) = &
!               & 0.5_wp * dual_edge_length(edge_index, edge_block) * length
! 
!           ENDIF !(edge_block > 0) THEN
! 
!         ENDDO !neigbor=1,6
! 
!       ENDDO ! vertex_index = start_index, end_index
!     ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
!     !-------------------
!     ! sync the results
!     DO neigbor=1,no_dual_edges
!       CALL sync_patch_array(SYNC_V, patch_2D, edge2vert_coeff_cc(:,:,neigbor)%x(1))
!       CALL sync_patch_array(SYNC_V, patch_2D, edge2vert_coeff_cc(:,:,neigbor)%x(2))
!       CALL sync_patch_array(SYNC_V, patch_2D, edge2vert_coeff_cc(:,:,neigbor)%x(3))
!       CALL sync_patch_array(SYNC_V, patch_2D, variable_dual_vol_norm(:,:, neigbor))
!     ENDDO ! neigbor=1,6
!     ! edge2vert_coeff_cc
!     ! variable_dual_vol_norm
!     !   are computed
!     !----------------------------------------------------
! 
!     !----------------------------------------------------
!     ! compute:
!     !   edge2vert_coeff_cc_t
!     DO edge_block = owned_edges%start_block, owned_edges%end_block
!       CALL get_index_range(owned_edges, edge_block, start_index, end_index)
!       DO edge_index = start_index, end_index
! 
!         edge_center%x = dual_edge_middle(edge_index, edge_block)%x
! 
!         DO neigbor=1,2
! 
!           edge2vert_coeff_cc_t(edge_index, edge_block, neigbor)%x = 0.0_wp
! 
!           vertex_index = patch_2D%edges%vertex_idx(edge_index, edge_block, neigbor)
!           vertex_block = patch_2D%edges%vertex_blk(edge_index, edge_block, neigbor)
! 
!           edge2vert_coeff_cc_t(edge_index, edge_block, neigbor)%x =              &
!             & (edge_center%x - patch_2D%verts%cartesian(vertex_index, vertex_block)%x) * &
!             & patch_2D%edges%system_orientation(edge_index, edge_block)                / &
!             & prime_edge_length(edge_index, edge_block)
! 
!         ENDDO !neigbor=1,2
! 
!       ENDDO ! edge_index = start_index, end_index
!     ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
!     !-------------------
!     ! sync the results
!     DO neigbor=1,2
!       CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(1))
!       CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(2))
!       CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(3))
!     ENDDO ! neigbor=1,2
!     ! edge2vert_coeff_cc_t is computed
!     !----------------------------------------------------
! 
!     !----------------------------------------------------
!     ! recalculate the coriolis coefficient
!     ! It is required if we use the middle of the dual_edge_length
!     IF (MID_POINT_DUAL_EDGE) THEN
! 
!       IF (CORIOLIS_TYPE == full_coriolis) THEN
! 
!         DO edge_block = owned_edges%start_block, owned_edges%end_block
!           CALL get_index_range(owned_edges, edge_block, start_index, end_index)
!           DO edge_index = start_index, end_index
! 
!              coriolis_geo_coordinates = cc2gc(dual_edge_middle(edge_index,edge_block))
!              patch_2D%edges%f_e(edge_index,edge_block) = &
!                & 2._wp * grid_angular_velocity * SIN(coriolis_geo_coordinates%lat)
! 
!           ENDDO
!         ENDDO
! 
!       ELSEIF (CORIOLIS_TYPE == BETA_PLANE_CORIOLIS) THEN
! 
!         basin_center_lat_rad = basin_center_lat * deg2rad
!         basin_height_rad     = basin_height_deg * deg2rad
!         coriolis_geo_coordinates%lat = basin_center_lat_rad - 0.5_wp * basin_height_rad
!         coriolis_geo_coordinates%lon = 0.0_wp
!         coriolis_cartesian_coordinates  = gc2cc(coriolis_geo_coordinates)
! 
!         DO edge_block = owned_edges%start_block, owned_edges%end_block
!           CALL get_index_range(owned_edges, edge_block, start_index, end_index)
!           DO edge_index = start_index, end_index
! 
!           geo_coordinates     = cc2gc(dual_edge_middle(edge_index,edge_block))
!           geo_coordinates%lon = 0.0_wp
!           edge_center         = gc2cc(geo_coordinates)
!           length              = grid_sphere_radius * &
!             & arc_length(edge_center, coriolis_cartesian_coordinates)
! 
!           patch_2D%edges%f_e(edge_index,edge_block) =  2.0_wp * grid_angular_velocity * &
!             & ( sin(basin_center_lat_rad) + (cos(basin_center_lat_rad) / &
!             &   grid_sphere_radius) * length)
! 
!           ENDDO
!         ENDDO
! 
!       ENDIF !(CORIOLIS_TYPE==full_coriolis)
!     ENDIF ! (MID_POINT_DUAL_EDGE)
!     !-------------------
!     ! sync patch_2D%edges%f_e
!     CALL sync_patch_array(SYNC_E, patch_2D, patch_2D%edges%f_e)
! 
! 
!     !----------------------------------------------------
! !    DEALLOCATE( prime_edge_length)
! 
!     DEALLOCATE( dual_edge_length )
! !     DEALLOCATE( cell_area )
! !     DEALLOCATE( dual_cell_area )
!     !---------------------------------------------------------
! !     RETURN
!     !---------------------------------------------------------
!     ! checks
! !
! !      !---------------------------------------------------------
! !      check_v1 = intp_2D_coeff%edge2vert_coeff_cc
! !      check_v2 = intp_2D_coeff%edge2vert_coeff_cc_t
! !      !---------------------------------------------------------
! !
! !      CALL init_scalar_product_oce( patch_2D, intp_2D_coeff )
! !
! !      !---------------------------------------------------------
! !      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc(:,:,:)%x(1) - &
! !        &  check_v1(:,:,:)%x(1) ))
! !      max_val  =  MAXVAL(ABS( check_v1(:,:,:)%x(1)))
! !      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(1)=", max_diff, max_val
! !
! !      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc(:,:,:)%x(2) - &
! !        & check_v1(:,:,:)%x(2) ))
! !      max_val  =  MAXVAL(ABS( check_v1(:,:,:)%x(2)))
! !      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(2)=", max_diff, max_val
! !
! !      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc(:,:,:)%x(3) - &
! !        & check_v1(:,:,:)%x(3) ))
! !      max_val  =  MAXVAL(ABS( check_v1(:,:,:)%x(3)))
! !      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(3)=", max_diff, max_val
! !      !---------------------------------------------------------
! !
! !      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc_t(:,:,:)%x(1) - &
! !        &  check_v2(:,:,:)%x(1) ))
! !      max_val  =  MAXVAL(ABS( check_v2(:,:,:)%x(1)))
! !      Write(0,*) "max diff of edge2vert_coeff_cc_t(:,:,:,:)%x(1)=", max_diff, max_val
! !
! !      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc_t(:,:,:)%x(2) - &
! !        & check_v2(:,:,:)%x(2) ))
! !      max_val  =  MAXVAL(ABS( check_v2(:,:,:)%x(2)))
! !      Write(0,*) "max diff of edge2vert_coeff_cc_t(:,:,:,:)%x(2)=", max_diff, max_val
! !
! !      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc_t(:,:,:)%x(3) - &
! !        & check_v2(:,:,:)%x(3) ))
! !      max_val  =  MAXVAL(ABS( check_v2(:,:,:)%x(3)))
! !      Write(0,*) "max diff of edge2vert_coeff_cc_t(:,:,:,:)%x(3)=", max_diff, max_val
! !
! 
!   END SUBROUTINE par_init_coeff_2D
! ! !-----------------------------------------------------------------------------------------
!   !-------------------------------------------------------------------------
!   !> Initialize expansion coefficients.
!   !!
!   !! @par Revision History
!   !! Peter Korn (2012-2)
!   !!
!   SUBROUTINE par_init_operator_coeff2( patch_3D, p_os, p_phys_param, operators_coefficients)
!     !
!     TYPE(t_patch_3D ),TARGET, INTENT(INOUT) :: patch_3D
!     TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os
!     TYPE (t_ho_params),       INTENT(IN)    :: p_phys_param
!     TYPE(t_operator_coeff),   INTENT(inout) :: operators_coefficients
! 
!    !
!    !Local variables: strcutures for 2D coefficients
!    !
!     TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc(1:nproma,1:patch_3D%p_patch_2D(1)%nblks_c,&
!                                      &1:patch_3D%p_patch_2D(1)%cell_type)
!     TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc_t(1:nproma,1:patch_3D%p_patch_2D(1)%nblks_e,1:2)
! 
!     TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc(1:nproma,1:patch_3D%p_patch_2D(1)%nblks_v,1:6)
!     TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc_t(1:nproma,1:patch_3D%p_patch_2D(1)%nblks_e,1:2)
!     TYPE(t_patch), POINTER :: patch_2D
! 
!     REAL(wp) :: dist_cell2edge(1:nproma,1:patch_3D%p_patch_2D(1)%nblks_e,1:2)
!     REAL(wp) :: fixed_vol_norm(1:nproma,patch_3D%p_patch_2D(1)%nblks_c)
!     REAL(wp) :: variable_vol_norm(1:nproma,1:patch_3D%p_patch_2D(1)%nblks_c,1:3)
!     REAL(wp) :: variable_dual_vol_norm(1:nproma,1:patch_3D%p_patch_2D(1)%nblks_e,1:6)
!     REAL(wp) :: edge2edge_viacell_coeff(1:nproma,1:patch_3D%p_patch_2D(1)%nblks_e,1:6)
!     REAL(wp) :: edge2edge_viavert_coeff(1:nproma,1:patch_3D%p_patch_2D(1)%nblks_e,1:12)
!     !-----------------------------------------------------------------------
!    patch_2D => patch_3D%p_patch_2D(1)
!    ! CALL init_operator_coeffs( patch_2D, operators_coefficients)
! 
! !--Old initialization of operator coefficients
!       CALL  par_init_coeff_2D( patch_2D,               &
!                           & edge2cell_coeff_cc,   &
!                           & edge2cell_coeff_cc_t, &
!                           & edge2vert_coeff_cc,   &
!                           & edge2vert_coeff_cc_t, &
!                           & dist_cell2edge,       &
!                           & fixed_vol_norm,       &
!                           & variable_vol_norm,    &
!                           & variable_dual_vol_norm)
! 
!     CALL copy_2D_to_3D_coeff( patch_2D,                  &
!                             & operators_coefficients,            &
!                             & edge2edge_viacell_coeff,&
!                             & edge2edge_viavert_coeff,&
!                             & edge2cell_coeff_cc,     &
!                             & edge2cell_coeff_cc_t,   &
!                             & edge2vert_coeff_cc,     &
!                             & edge2vert_coeff_cc_t,   &
!                             & dist_cell2edge,         &
!                             & fixed_vol_norm,         &
!                             & variable_vol_norm,      &
!                             & variable_dual_vol_norm)
! 
!     CALL init_diff_operator_coeff_3D ( patch_2D, operators_coefficients )
!     
!     CALL apply_boundary2coeffs(patch_3D, operators_coefficients)
! 
! 
!      CALL update_diffusion_matrices(   patch_3D,                    &
!                                      & p_os,                          &
!                                      & p_phys_param,                  &
!                                      & operators_coefficients%matrix_vert_diff_e,&
!                                      & operators_coefficients%matrix_vert_diff_c)
!   END SUBROUTINE par_init_operator_coeff2
!   !-------------------------------------------------------------------------


END MODULE mo_operator_ocean_coeff_3d

