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
MODULE mo_operator_ocean_coeff_3d
  !-------------------------------------------------------------------------

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert, max_dom,success,&
    & max_char_length, beta_plane_coriolis,full_coriolis, &
    & min_rledge_int,min_rlcell_int,min_rlvert_int,&
    & sea_boundary, boundary, sea
  USE mo_math_constants,      ONLY: deg2rad, pi
  USE mo_physical_constants,  ONLY: re,omega
  USE mo_math_utilities,      ONLY: gc2cc, cc2gc, t_cartesian_coordinates,      &
    & t_geographical_coordinates, vector_product, &
    & arc_length
  USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m, no_tracer, t_ref, s_ref,          &
    & coriolis_type, basin_center_lat, basin_height_deg
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, sync_idx, global_max
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_oce_state,           ONLY: v_base
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_coeffs,         ONLY: par_init_scalar_product_oce
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  IMPLICIT NONE


  PRIVATE

  PUBLIC :: allocate_exp_coeff
  PUBLIC :: init_operator_coeff, par_init_operator_coeff
  PUBLIC :: t_operator_coeff
  PUBLIC :: apply_boundary2coeffs
  PRIVATE :: init_scalar_product_oce_3d
  PRIVATE :: init_geo_factors_oce_3d

  INTEGER,PARAMETER :: no_dual_edges   = 6
  INTEGER,PARAMETER :: no_primal_edges = 3


  TYPE t_operator_coeff

    ! 1) precomputed 3D-factors for mathematical operators (for efficiency).
    !------------------------------------------------------------------------------
    REAL(wp), ALLOCATABLE :: div_coeff(:,:,:,:)    ! factor for divergence (nproma,nlev,nblks_c,no_primal_edges)
    REAL(wp), ALLOCATABLE :: rot_coeff(:,:,:,:)    ! factor for divergence (nproma,nlev,nblks_v,no_dual_edges)
    REAL(wp), ALLOCATABLE :: grad_coeff(:,:,:)     ! factor for nabla2-scalar (nproma,nlev,nblks_e)
    REAL(wp), ALLOCATABLE :: n2s_coeff(:,:,:,:)    ! factor for nabla2-scalar (nproma,nlev,nblks_c)
    REAL(wp), ALLOCATABLE :: n2v_coeff(:,:,:)      ! factor for nabla2-vector (nproma,nlev,nblks_e)


    !2) Required for description of boundary around a vertex
    !------------------------------------------------------------------------------
    INTEGER, ALLOCATABLE :: bnd_edges_per_vertex(:,:,:)
    INTEGER, POINTER :: bnd_edge_idx(:,:,:,:)!(nproma,nlev,nblks_v,1:NO_DUAL_EDGES-2)
    INTEGER, POINTER :: bnd_edge_blk(:,:,:,:)!(nproma,nlev,nblks_v,1:NO_DUAL_EDGES-2)
    INTEGER, POINTER :: edge_idx(:,:,:,:)   !(nproma,nlev,nblks_v,1:NO_DUAL_EDGES-2)
    REAL(wp),POINTER :: orientation(:,:,:,:)!(nproma,nlev,nblks_v,1:NO_DUAL_EDGES-2)

    INTEGER, ALLOCATABLE :: upwind_cell_idx(:,:,:)
    INTEGER, ALLOCATABLE :: upwind_cell_blk(:,:,:)
    !3) Scalarproduct: The following arrays are required for the reconstruction process.
    !------------------------------------------------------------------------------
    !
    ! Vector pointing from cell circumcenter to edge midpoint. In the associated
    ! cell2edge_weight-array the cell2edge_vec is multiplied by some other geometric
    ! quantities (edge-length, cell area). The weight is used in the reconstruction
    ! the vector is used in the transposed reconstruction.
    ! index=1,nproma, index2=1,nblks_c, index3=1,3
    ! other choice would be index2=1,nblks_e, index3=1,2
    ! Eventually switch to other second indexing if this is more appropriate

    ! Vector pointing from vertex (dual center) to midpoint of dual edge
    ! (/= midpoint of primal edge).
    ! In the associated vertex2dualedge_mid_weight-array the vertex2dualedge_mid_vec
    ! is multiplied by some other geometric quantities (dual edge-length, dual cell
    ! area). The weight is used in the reconstruction the vector is used in the
    ! transposed reconstruction.
    ! index=1,nproma, index2=1,nblks_v, index3=1,6
    ! other choice index2=1,nblks_e, index3=1,2
    ! Eventually switch to other second indexing if this is more appropriate
    ! new constructs for mimetic core:
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2cell_coeff_cc(:,:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2cell_coeff_cc_t(:,:,:,:)

    !coefficient for surface layer, changes in time, in contrast to other coefficients
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2cell_coeff_cc_dyn(:,:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc_dyn(:,:,:,:)

    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc(:,:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc_t(:,:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_vector_cc(:,:,:,:)

    REAL(wp), ALLOCATABLE :: fixed_vol_norm(:,:,:)
    REAL(wp), ALLOCATABLE :: variable_vol_norm(:,:,:,:)
    REAL(wp), ALLOCATABLE :: variable_dual_vol_norm(:,:,:,:)

    !!$    TYPE(t_geographical_coordinates), ALLOCATABLE :: mid_dual_edge(:,:)
    ! Cartesian distance from vertex1 to vertex2 via dual edge midpoint
    REAL(wp), ALLOCATABLE :: dist_cell2edge(:,:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: cell_position_cc(:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge_position_cc(:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: moved_edge_position_cc(:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: upwind_cell_position_cc(:,:,:)

  END TYPE t_operator_coeff

CONTAINS

  ! !-------------------------------------------------------------------------
  ! !
  ! !
  ! !> Allocation of expansion coefficients.
  ! !!
  ! !! @par Revision History
  ! !! Peter Korn (2012-2)
  ! !!
  SUBROUTINE allocate_exp_coeff( ptr_patch, ptr_coeff)
    ! !
    TYPE(t_patch),TARGET,INTENT(in)       :: ptr_patch
    TYPE(t_operator_coeff), INTENT(inout) :: ptr_coeff

    INTEGER :: nblks_c, nblks_e, nblks_v, nz_lev
    INTEGER :: ist,ie
    INTEGER :: rl_start,rl_end
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: jc,je,jk,jb

    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: all_verts

    INTEGER :: edge_block, cell_block, vertex_block, level, neigbor
    INTEGER :: i_startidx_e, i_endidx_e

    !-----------------------------------------------------------------------

    !
    ! determine size of arrays, i.e.
    ! values for the blocking
    !
    nblks_c  = ptr_patch%nblks_c
    nblks_e  = ptr_patch%nblks_e
    nblks_v  = ptr_patch%nblks_v
    nz_lev   = n_zlev

    ALLOCATE(ptr_coeff%div_coeff(nproma,n_zlev,nblks_c,no_primal_edges),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for geofac_div failed')
    ENDIF

    ALLOCATE(ptr_coeff%grad_coeff(nproma,n_zlev,nblks_e),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for geofac_grad failed')
    ENDIF

    ALLOCATE(ptr_coeff%rot_coeff(nproma,n_zlev,nblks_v,no_dual_edges),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d', &
        & 'allocation for geofac_rot failed')
    ENDIF

    ALLOCATE(ptr_coeff%n2s_coeff(nproma,n_zlev,nblks_c,no_primal_edges+1),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                       &
        & 'allocation for geofac_n2s failed')
    ENDIF
    ALLOCATE(ptr_coeff%n2v_coeff(nproma,n_zlev,nblks_e),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                       &
        & 'allocation for geofac_n2v failed')
    ENDIF
    !
    ALLOCATE(ptr_coeff%dist_cell2edge(nproma,n_zlev,nblks_e,2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating dist_cell2edge failed')
    ENDIF

    ALLOCATE(ptr_coeff%bnd_edge_idx(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating bnd_edge_idx failed')
    ENDIF
    ALLOCATE(ptr_coeff%bnd_edge_blk(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating bnd_edge_blk failed')
    ENDIF
    ALLOCATE(ptr_coeff%edge_idx(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge_idx failed')
    ENDIF
    ALLOCATE(ptr_coeff%orientation(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating orientation failed')
    ENDIF
    ALLOCATE(ptr_coeff%bnd_edges_per_vertex(nproma,n_zlev,nblks_v),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating bnd_edges_per_vertex failed')
    ENDIF
    ALLOCATE(ptr_coeff%upwind_cell_idx(nproma,n_zlev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind_cell_idx failed')
    ENDIF
    ALLOCATE(ptr_coeff%upwind_cell_blk(nproma,n_zlev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind_cell_blk failed')
    ENDIF
    !
    ! arrays that are required for setting up the scalar product
    !
    !coefficients for edge to cell mapping, one half of the scalar product.
    !Dimension: nproma,nblks_c encode number of cells, 1:3 corresponds to number
    !of edges per cell, 1:2 is for u and v component of cell vector
    !     ALLOCATE(ptr_coeff%edge2cell_coeff(nproma,nblks_c,1:3, 1:2),STAT=ist)
    !     IF (ist /= SUCCESS) THEN
    !       CALL finish ('allocating edge2cell_coeff failed')
    !     ENDIF
    ALLOCATE(ptr_coeff%edge2cell_coeff_cc(nproma,nz_lev,nblks_c,1:3),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc failed')
    ENDIF

    ALLOCATE(ptr_coeff%edge2cell_coeff_cc_dyn(nproma,1,nblks_c,1:3),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc_dyn failed')
    ENDIF
    ALLOCATE(ptr_coeff%edge2vert_coeff_cc_dyn(nproma,1,nblks_c,1:3),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff_cc_dyn failed')
    ENDIF

    !coefficients for transposed of edge to cell mapping, second half of the scalar product.
    !Dimension: nproma,nblks_e encode number of edges, 1:2 is for cell neighbors of an edge
    ALLOCATE(ptr_coeff%edge2cell_coeff_cc_t(nproma,nz_lev,nblks_e,1:2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating transposed edge2cell_coeff failed')
    ENDIF

    !
    !coefficients for edge to vertex mapping.
    !
    !Dimension: nproma,nblks_v encode number of vertices,
    !1:6 is number of edges of a vertex,
    !1:2 is for u and v component of vertex vector
    ALLOCATE(ptr_coeff%edge2vert_coeff_cc(nproma,nz_lev,nblks_v,1:6),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff failed')
    ENDIF

    ALLOCATE(ptr_coeff%edge2vert_coeff_cc_t(nproma,nz_lev,nblks_e,1:2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff failed')
    ENDIF
    ALLOCATE(ptr_coeff%edge2vert_vector_cc(nproma,nz_lev,nblks_v,1:6),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_vector failed')
    ENDIF

    ALLOCATE(ptr_coeff%upwind_cell_position_cc(nproma,nz_lev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind cell failed')
    ENDIF
    ALLOCATE(ptr_coeff%moved_edge_position_cc(nproma,nz_lev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge failed')
    ENDIF
    ALLOCATE(ptr_coeff%edge_position_cc(nproma,nz_lev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge failed')
    ENDIF
    ALLOCATE(ptr_coeff%cell_position_cc(nproma,nz_lev,nblks_c),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating cell failed')
    ENDIF
    !
    !normalizing factors for edge to cell mapping.
    !
    !Either by fixed volume or by variable one taking the surface elevation
    !into account. The later one depends on time and space.
    ALLOCATE(ptr_coeff%fixed_vol_norm(nproma,nz_lev,nblks_c),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating fixed_vol_norm failed')
    ENDIF
    ALLOCATE(ptr_coeff%variable_vol_norm(nproma,nz_lev,nblks_c,1:3),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating variable_vol_norm failed')
    ENDIF

    ALLOCATE(ptr_coeff%variable_dual_vol_norm(nproma,nz_lev,nblks_v,1:6),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating variable_dual_vol_norm failed')
    ENDIF
    !   ALLOCATE (ptr_coeff%n2s_coeff(nproma, nz_lev, nblks_c,3+1), STAT=ist )
    !   IF (ist /= SUCCESS) THEN
    !     CALL finish ('mo_operator_ocean_coeff_3d',                       &
    !       &             'allocation for geofac_n2s failed')
    !   ENDIF
    !   ALLOCATE (ptr_coeff%n2v_coeff(nproma, nz_lev, nblks_e), STAT=ist )
    !   IF (ist /= SUCCESS) THEN
    ! write(*,*)'data',nproma, nz_lev, nblks_e, ist
    !     CALL finish ('mo_operator_ocean_coeff_3d',                       &
    !       &             'allocation for geofac_n2v failed')
    !   ENDIF
    !
    ! initialize all components
    !
    DO ie = 1,3
      ptr_coeff%edge2cell_coeff_cc%x(ie)     = 0._wp
      ptr_coeff%edge2cell_coeff_cc_t%x(ie)   = 0._wp
      ptr_coeff%edge2vert_coeff_cc%x(ie)     = 0._wp
      ptr_coeff%edge2vert_coeff_cc_t%x(ie)   = 0._wp
      ptr_coeff%edge2vert_vector_cc%x(ie)    = 0._wp
      ptr_coeff%edge2cell_coeff_cc_dyn%x(ie) = 0._wp
      ptr_coeff%edge2vert_coeff_cc_dyn%x(ie) = 0._wp
    END DO

    all_cells => ptr_patch%cells%all
    all_edges => ptr_patch%edges%all
    all_verts => ptr_patch%verts%all

    DO jk = 1, nz_lev
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO je =  i_startidx_e, i_endidx_e
          ptr_coeff%edge_position_cc(je,jk,jb)             = gc2cc(ptr_patch%edges%center(je,jb))
          ptr_coeff%moved_edge_position_cc(je,jk,jb)%x(:)  = 0._wp
          ptr_coeff%upwind_cell_position_cc(je,jk,jb)%x(:) = 0._wp
        END DO
      END DO
    END DO

    DO jk = 1, nz_lev
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          ptr_coeff%cell_position_cc(jc,jk,jb)&
            & = gc2cc(ptr_patch%cells%center(jc,jb))

        END DO
      END DO
    END DO

    ptr_coeff%fixed_vol_norm         = 0._wp
    ptr_coeff%variable_vol_norm      = 0._wp
    ptr_coeff%variable_dual_vol_norm = 0._wp

    ptr_coeff%dist_cell2edge = 0._wp

    ptr_coeff%div_coeff  = 0._wp
    ptr_coeff%rot_coeff  = 0._wp
    ptr_coeff%grad_coeff = 0._wp
    !ptr_coeff%n2s_coeff  = 0._wp
    !ptr_coeff%n2v_coeff  = 0._wp

    ptr_coeff%bnd_edge_idx = 0
    ptr_coeff%bnd_edge_blk = 0
    ptr_coeff%edge_idx     = 0
    ptr_coeff%orientation  = 0.0_wp
    ptr_coeff%bnd_edges_per_vertex= 0

    ptr_coeff%upwind_cell_idx = 1
    ptr_coeff%upwind_cell_blk = 1

    CALL message ('mo_operator_ocean_coeff_3d:allocate_exp_coeff',&
      & 'memory allocation finished')

  END SUBROUTINE allocate_exp_coeff
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !> Initialize expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !!
  SUBROUTINE par_init_operator_coeff( patch, ocean_coeff, intp_2D_coeff)
    !
    TYPE(t_patch),          INTENT(inout) :: patch
    TYPE(t_operator_coeff), INTENT(inout) :: ocean_coeff
    TYPE(t_int_state),      INTENT(inout) :: intp_2D_coeff

    INTEGER :: rl_start_e,rl_end_e,rl_start_c,rl_end_c
    INTEGER :: i_startblk_e, i_endblk_e,i_startidx_e,i_endidx_e
    INTEGER :: i_startblk_c, i_endblk_c,i_startidx_c,i_endidx_c
    INTEGER :: jk,je,jb,jc
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
!     TYPE(t_cartesian_coordinates) :: check_v(nproma, n_zlev, patch%nblks_v, 6)
!     REAL(wp) :: check_r(nproma, n_zlev, patch%nblks_c, 3)
!     REAL(wp) :: max_diff, max_val
    !-----------------------------------------------------------------------
    ocean_coeff%dist_cell2edge(:,:,:,:) = 0.0_wp
    CALL par_init_scalar_product_oce(patch, intp_2D_coeff)
    CALL copy_2D_to_3D_intp_coeff( patch, ocean_coeff, intp_2D_coeff)
    CALL init_geo_factors_oce_3d ( patch, ocean_coeff )
    !---------------------------------------------------------
    CALL par_apply_boundary2coeffs(patch, ocean_coeff)

    RETURN
    !---------------------------------------------------------
    ! checks
!     check_v = ocean_coeff%edge2vert_coeff_cc
!     check_r = ocean_coeff%variable_vol_norm
!
!     !---------------------------------------------------------
!      CALL init_operator_coeff( patch, ocean_coeff )
!
!      max_diff = MAXVAL(ABS(ocean_coeff%edge2vert_coeff_cc(:,:,:,:)%x(1) - &
!        &  check_v(:,:,:,:)%x(1) ))
!      max_val  =  MAXVAL(ABS( check_v(:,:,:,:)%x(1)))
!      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(1)=", max_diff, max_val
!
!      max_diff = MAXVAL(ABS(ocean_coeff%edge2vert_coeff_cc(:,:,:,:)%x(2) - &
!        & check_v(:,:,:,:)%x(2) ))
!      max_val  =  MAXVAL(ABS( check_v(:,:,:,:)%x(2)))
!      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(2)=", max_diff, max_val
!
!      max_diff = MAXVAL(ABS(ocean_coeff%edge2vert_coeff_cc(:,:,:,:)%x(3) - &
!        & check_v(:,:,:,:)%x(3) ))
!      max_val  =  MAXVAL(ABS( check_v(:,:,:,:)%x(3)))
!      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(3)=", max_diff, max_val
!
!      max_diff = MAXVAL(ABS(ocean_coeff%variable_vol_norm - check_r ))
!      max_val  =  MAXVAL(ABS( check_r ))
!      Write(0,*) "max diff of variable_vol_norm=", max_diff, max_val
!
!      STOP

  END SUBROUTINE par_init_operator_coeff
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Initialize 3D expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !! Parellelized by Leonidas Linardakis 2012-3
  SUBROUTINE copy_2D_to_3D_intp_coeff( patch, ocean_coeff, intp_2D_coeff)
    !
    TYPE(t_patch),  TARGET, INTENT(inout) :: patch
    TYPE(t_operator_coeff), INTENT(inout) :: ocean_coeff
    TYPE(t_int_state),      INTENT(inout) :: intp_2D_coeff

    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: all_verts

    INTEGER :: edge_block, cell_block, vertex_block, level, neigbor
    INTEGER :: je,jk, i_startidx_e, i_endidx_e

    all_cells => patch%cells%all
    all_edges => patch%edges%all
    all_verts => patch%verts%all

    !---------------------------------------------------------
    ! the following coefficients will be copied:
    !
    ! ocean_coeff%edge_position_cc(:,:,:)            on edges
    ! ocean_coeff%dist_cell2edge(:,:,:,1-2)          on edges
    !
    ! ocean_coeff%edge2cell_coeff_cc(:,:,:,1-3)%x    on cells
    ! ocean_coeff%fixed_vol_norm(:,:,:)              on cells
    ! ocean_coeff%variable_vol_norm(:,:,:,1-3)       on cells
    !
    ! ocean_coeff%edge2cell_coeff_cc_t(:,:,:,1-2)%x  on edges
    ! ocean_coeff%edge2vert_coeff_cc(:,:,:,1-6)%x   on verts
    ! ocean_coeff%edge2vert_coeff_cc_t(:,:,:,1-2)%x  on edges
    !
    ! ptr_patch%edges%f_e(:, :) is already calculated in par_init_scalar_product_oce    !
    !---------------------------------------------------------


    !---------------------------------------------------------
    ! calculate ocean_coeff%edge_position_cc(:,:,:) on edges
    ! this is the same as the 2D, just copy it
    ! it does not change, maybe turn it to 2D?
    DO edge_block = all_edges%start_block, all_edges%end_block
      DO level = 1, n_zlev
        ocean_coeff%edge_position_cc(:,level,edge_block) = &
          patch%edges%cartesian_center(:,edge_block)
      ENDDO
    ENDDO
    !---------------------------------------------------------

    !---------------------------------------------------------
    ! calculate ocean_coeff%dist_cell2edge(:,:,:,1-2) on edges
    ! this is the same as the 2D, just copy it
    ! it does not change, maybe turn it to 2D?
    DO edge_block = all_edges%start_block, all_edges%end_block
      DO level = 1, n_zlev
        ocean_coeff%dist_cell2edge(:,level,edge_block,1) = &
          intp_2D_coeff%dist_cell2edge(:,edge_block,1)
        ocean_coeff%dist_cell2edge(:,level,edge_block,2) = &
          intp_2D_coeff%dist_cell2edge(:,edge_block,2)
      ENDDO
    ENDDO
    !---------------------------------------------------------

    !---------------------------------------------------------
    ! For the rest of coefficients,
    ! copy the 2D and then apply_boundary2coeffs
    !
    ! on edges
    DO edge_block = all_edges%start_block, all_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(all_edges, edge_block, i_startidx_e, i_endidx_e)
        DO je =  i_startidx_e, i_endidx_e

          ocean_coeff%edge2cell_coeff_cc_t(je,level,edge_block,1)%x = &
            intp_2D_coeff%edge2cell_coeff_cc_t(je,edge_block,1)%x
          ocean_coeff%edge2cell_coeff_cc_t(je,level,edge_block,2)%x = &
            intp_2D_coeff%edge2cell_coeff_cc_t(je,edge_block,2)%x

          ocean_coeff%edge2vert_coeff_cc_t(je,level,edge_block,1)%x = &
            intp_2D_coeff%edge2vert_coeff_cc_t(je,edge_block,1)%x
          ocean_coeff%edge2vert_coeff_cc_t(je,level,edge_block,2)%x = &
            intp_2D_coeff%edge2vert_coeff_cc_t(je,edge_block,2)%x

        ENDDO
      ENDDO
    ENDDO

    ! on cells
    DO cell_block = all_cells%start_block, all_cells%end_block
      DO level = 1, n_zlev

       ocean_coeff%fixed_vol_norm(:,level,cell_block) = &
         intp_2D_coeff%fixed_vol_norm(:,cell_block)

       DO neigbor=1,patch%cell_type

         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(1) = &
           & intp_2D_coeff%edge2cell_coeff_cc(:,cell_block,neigbor)%x(1)
         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(2) = &
           & intp_2D_coeff%edge2cell_coeff_cc(:,cell_block,neigbor)%x(2)
         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(3) = &
           & intp_2D_coeff%edge2cell_coeff_cc(:,cell_block,neigbor)%x(3)

         ocean_coeff%variable_vol_norm(:,level,cell_block,neigbor) = &
           & intp_2D_coeff%variable_vol_norm(:,cell_block,neigbor)

        ENDDO ! neigbor=1,patch%cell_type

      ENDDO  !  level = 1, n_zlev
    ENDDO ! cell_block

    ! on verts
    DO vertex_block = all_verts%start_block, all_verts%end_block
      DO level = 1, n_zlev
        DO neigbor=1,6

          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(1) = &
            & intp_2D_coeff%edge2vert_coeff_cc(:,vertex_block,neigbor)%x(1)
          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(2) = &
            & intp_2D_coeff%edge2vert_coeff_cc(:,vertex_block,neigbor)%x(2)
          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(3) = &
            & intp_2D_coeff%edge2vert_coeff_cc(:,vertex_block,neigbor)%x(3)

        ENDDO ! neigbor=1,patch%cell_type
      ENDDO  !  level = 1, n_zlev
    ENDDO ! vertex_block
    !---------------------------------------------------------

  END SUBROUTINE copy_2D_to_3D_intp_coeff
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Initialize expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !!
  SUBROUTINE par_apply_boundary2coeffs( patch, ocean_coeff)
    ! !
    TYPE(t_patch), TARGET,  INTENT(inout) :: patch
    TYPE(t_operator_coeff), INTENT(inout) :: ocean_coeff

    !Local variables
    INTEGER :: jk, jc, jb, je, ibe, ile, jev, jv
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_v, i_endidx_v

    INTEGER :: i_v_ctr(nproma,n_zlev,patch%nblks_v)
    INTEGER :: ibnd_edge_idx(4), ibnd_edge_blk(4)  !maximal 4 boundary edges in a dual loop.
    INTEGER :: i_edge_idx(4)
    REAL(wp) :: z_orientation(4),z_area_scaled, zarea_fraction
    INTEGER :: icell_idx_1, icell_blk_1
    INTEGER :: icell_idx_2, icell_blk_2
    INTEGER :: boundary_counter

    TYPE(t_cartesian_coordinates) :: cell1_cc, cell2_cc, vertex_cc

    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: owned_verts

 !   REAL(wp) :: comm_edges(nproma,n_zlev,patch%nblks_e)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:apply_boundary2coeffs')

    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')

    all_cells   => patch%cells%all
    all_edges   => patch%edges%all
    owned_verts => patch%verts%owned

    i_v_ctr(:,:,:)          = 0
    ! this should not be done here
!    comm_edges = REAL(v_base%lsm_oce_e,wp)
!    CALL sync_patch_array(SYNC_E, patch,  comm_edges )
!    v_base%lsm_oce_e = INT(comm_edges)
!      DO je=1,2
!        DO jk=1,3
!          CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2cell_coeff_cc_t(:,:,:,je)%x(jk))
!          CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2vert_coeff_cc_t(:,:,:,je)%x(jk))
!        ENDDO
!      ENDDO

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

          IF(v_base%lsm_oce_e(je,jk,jb) /= sea) THEN
            ocean_coeff%grad_coeff          (je,jk,jb) = 0.0_wp
            ocean_coeff%edge2cell_coeff_cc_t(je,jk,jb,1)%x(1:3) = 0.0_wp
            ocean_coeff%edge2cell_coeff_cc_t(je,jk,jb,2)%x(1:3) = 0.0_wp
            ocean_coeff%edge2vert_coeff_cc_t(je,jk,jb,1)%x(1:3) = 0.0_wp
            ocean_coeff%edge2vert_coeff_cc_t(je,jk,jb,2)%x(1:3) = 0.0_wp
          ENDIF

        ENDDO
      END DO
    END DO
    ! these are computed an all edges,  thus no sync is required
    ! we sync only in p_test_run for checking
    IF (p_test_run) THEN
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%grad_coeff(:,:,:))
      DO jk = 1, n_zlev
        DO je=1,2
          DO jb=1,3
            CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2cell_coeff_cc_t(:,jk,:,je)%x(jb))
            CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2vert_coeff_cc_t(:,jk,:,je)%x(jb))
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !-------------------------------------------------------------

    !-------------------------------------------------------------
    !1) Set coefficients for div and grad to zero at boundary edges
    ! Aslo for ocean_coeff%edge2cell_coeff_cc
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jk=1,n_zlev
        DO jc = i_startidx_c, i_endidx_c

          DO je = 1, patch%cells%num_edges(jc,jb)

            ile = patch%cells%edge_idx(jc,jb,je)
            ibe = patch%cells%edge_blk(jc,jb,je)
            IF ( v_base%lsm_oce_e(ile,jk,ibe) /= sea) THEN

              ocean_coeff%div_coeff(jc,jk,jb,je) = 0.0_wp
              ocean_coeff%edge2cell_coeff_cc(jc,jk,jb,je)%x(1:3) = 0.0_wp

            ENDIF

          ENDDO ! je = 1, patch%cells%num_edges(jc,jb)
          !write(1234,*)'div coeff 3D',jk,jc,jb,ptr_intp%div_coeff(jc,jk,jb,:)
        ENDDO ! jc = i_startidx_c, i_endidx_c
      END DO ! jk=1,n_zlev
    END DO ! jb = all_cells%start_block, all_cells%end_block

    ! these are computed an all cells,  thus no sync is required
    ! we sync only in p_test_run for checking
!    IF (p_test_run) THEN
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%grad_coeff(:,:,:))

      DO je=1,patch%cell_type
        CALL sync_patch_array(SYNC_C, patch, ocean_coeff%div_coeff(:,:,:, je))

        DO jk = 1, n_zlev
          CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc(:,jk,:,je)%x(1))
          CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc(:,jk,:,je)%x(2))
          CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc(:,jk,:,je)%x(3))
        ENDDO

      ENDDO
!    ENDIF
    !-------------------------------------------------------------

    !-------------------------------------------------------------
    !2) prepare coefficients for rot at boundary edges
    ! this is done on owned vertices, so sync is required
    DO jb = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, jb, i_startidx_v, i_endidx_v)
      DO jk = 1, n_zlev
        DO jv = i_startidx_v, i_endidx_v

          ibnd_edge_idx(1:4)      = 0
          ibnd_edge_blk(1:4)      = 0
          i_edge_idx(1:4)         = 0
          z_orientation(1:4)      = 0.0_wp
          boundary_counter        = 0

          DO jev = 1, patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = patch%verts%edge_idx(jv,jb,jev)
            ibe = patch%verts%edge_blk(jv,jb,jev)
            !Check, if edge is sea or boundary edge and take care of dummy edge
            ! edge with indices ile, ibe is sea edge
            ! edge with indices ile, ibe is boundary edge

            IF ( v_base%lsm_oce_e(ile,jk,ibe) == sea) THEN
              i_v_ctr(jv,jk,jb)=i_v_ctr(jv,jk,jb)+1
            ELSEIF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN

              !increase boundary edge counter
!               i_v_bnd_edge_ctr(jv,jk,jb)=i_v_bnd_edge_ctr(jv,jk,jb)+1
              boundary_counter = boundary_counter + 1

              ocean_coeff%bnd_edges_per_vertex(jv,jk,jb) &
                & = ocean_coeff%bnd_edges_per_vertex(jv,jk,jb) +1

              IF (boundary_counter > 4) THEN
                !maximal 4 boundary edges per dual loop are allowed: somethings wrong with the grid
                CALL message (TRIM('sbr nonlinear Coriolis'), &
                  & 'more than 4 boundary edges per dual loop: something is wrong with the grid')
                CALL finish (routine,'Grid-boundary error !!')
              ENDIF

              ibnd_edge_idx(boundary_counter) = ile
              ibnd_edge_blk(boundary_counter) = ibe
              z_orientation(boundary_counter) = patch%verts%edge_orientation(jv,jb,jev)
              i_edge_idx(boundary_counter)    = jev

              ocean_coeff%bnd_edge_idx(jv,jk,jb,boundary_counter)= ile
              ocean_coeff%bnd_edge_blk(jv,jk,jb,boundary_counter)= ibe
              ocean_coeff%orientation(jv,jk,jb,boundary_counter) = &
                & patch%verts%edge_orientation(jv,jb,jev)
              ocean_coeff%edge_idx(jv,jk,jb,boundary_counter)    = jev

!               write(0,*) "boundary_counter, jev=", boundary_counter, jev

               !LL: this should be equivelant to the above
!               IF(i_v_bnd_edge_ctr(jv,jk,jb)==1)THEN
!                 ibnd_edge_idx(1) = ile
!                 ibnd_edge_blk(1) = ibe
!                 z_orientation(1) = patch%verts%edge_orientation(jv,jb,jev)
!                 i_edge_idx(1)    = jev
!
!                 ocean_coeff%bnd_edge_idx(jv,jk,jb,1)= ile
!                 ocean_coeff%bnd_edge_blk(jv,jk,jb,1)= ibe
!                 ocean_coeff%orientation(jv,jk,jb,1) = patch%verts%edge_orientation(jv,jb,jev)
!                 ocean_coeff%edge_idx(jv,jk,jb,1)    = jev
!
!               ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==2)THEN
!                 ibnd_edge_idx(2) = ile
!                 ibnd_edge_blk(2) = ibe
!                 z_orientation(2) = patch%verts%edge_orientation(jv,jb,jev)
!                 i_edge_idx(2)    = jev
!
!                 ocean_coeff%bnd_edge_idx(jv,jk,jb,2)= ile
!                 ocean_coeff%bnd_edge_blk(jv,jk,jb,2)= ibe
!                 ocean_coeff%orientation(jv,jk,jb,2) = patch%verts%edge_orientation(jv,jb,jev)
!                 ocean_coeff%edge_idx(jv,jk,jb,2)    = jev
!
!               ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==3)THEN
!                 ibnd_edge_idx(3) = ile
!                 ibnd_edge_blk(3) = ibe
!                 z_orientation(3) = patch%verts%edge_orientation(jv,jb,jev)
!                 i_edge_idx(3)    = jev
!
!                 ocean_coeff%bnd_edge_idx(jv,jk,jb,3)= ile
!                 ocean_coeff%bnd_edge_blk(jv,jk,jb,3)= ibe
!                 ocean_coeff%orientation(jv,jk,jb,3) = patch%verts%edge_orientation(jv,jb,jev)
!                 ocean_coeff%edge_idx(jv,jk,jb,3)    = jev
!
!               ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==4)THEN
!                 ibnd_edge_idx(4) = ile
!                 ibnd_edge_blk(4) = ibe
!                 z_orientation(4) = patch%verts%edge_orientation(jv,jb,jev)
!                 i_edge_idx(4)    = jev
!
!                 ocean_coeff%bnd_edge_idx(jv,jk,jb,4)= ile
!                 ocean_coeff%bnd_edge_blk(jv,jk,jb,4)= ibe
!                 ocean_coeff%orientation(jv,jk,jb,4) = patch%verts%edge_orientation(jv,jb,jev)
!                 ocean_coeff%edge_idx(jv,jk,jb,4)    = jev
!
!               ELSE
!                 !maximal 4 boundary edges per dual loop are allowed: somethings wrong withe the grid
!                 CALL message (TRIM('sbr nonlinear Coriolis'), &
!                   & 'more than 4 boundary edges per dual loop: something is wrong with the grid')
!                 CALL finish ('TRIM(sbr nonlinear Coriolis)','Grid-boundary error !!')
!               ENDIF
            END IF
          END DO ! jev = 1, patch%verts%num_edges(jv,jb)

          IF( MOD(boundary_counter,2) /= 0 ) THEN
            CALL finish (routine,'MOD(boundary_counter,2) /= 0 !!')
          ENDIF

          DO je = 1, boundary_counter
!-original
!             ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(je) )=&
!               & 0.5_wp*z_orientation(je) * &
!               & patch%edges%primal_edge_length(ibnd_edge_idx(je),ibnd_edge_blk(je))

            ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(je) )=&
              & 0.5_wp*patch%edges%system_orientation(ibnd_edge_idx(je),ibnd_edge_blk(je)) * &
              & patch%edges%primal_edge_length(ibnd_edge_idx(je),ibnd_edge_blk(je))

          ENDDO

          ! LL: the above loop should be equivelant to the following
!           IF(boundary_counter == 2)THEN
!             ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(1) )=&
!               & 0.5_wp*z_orientation(1)&
!               & *patch%edges%primal_edge_length(ibnd_edge_idx(1),ibnd_edge_blk(1))
!
!             ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(2))=&
!               & 0.5_wp*z_orientation(2)&
!               & *patch%edges%primal_edge_length(ibnd_edge_idx(2),ibnd_edge_blk(2))
!
!           ELSEIF(boundary_counter == 4)THEN
!
!             ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(1))=&
!               & 0.5_wp*z_orientation(1)&
!               & *patch%edges%primal_edge_length(ibnd_edge_idx(1),ibnd_edge_blk(1))
!
!             ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(2))=&
!               & 0.5_wp*z_orientation(2)&
!               & *patch%edges%primal_edge_length(ibnd_edge_idx(2),ibnd_edge_blk(2))
!
!             ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(3))=&
!               & 0.5_wp*z_orientation(3)&
!               & *patch%edges%primal_edge_length(ibnd_edge_idx(3),ibnd_edge_blk(3))
!
!             ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(4))=&
!               & 0.5_wp*z_orientation(4)&
!               & *patch%edges%primal_edge_length(ibnd_edge_idx(4),ibnd_edge_blk(4))
!           ENDIF
          !ocean_coeff%rot_coeff(jv,jk,jb,:)=ocean_coeff%rot_coeff(jv,jk,jb,:)/patch%verts%dual_area(jv,jb)

        END DO ! jv = i_startidx_v, i_endidx_v
        !!$OMP END PARALLEL DO

      END DO ! jk = 1, n_zlev
    END DO ! jb = owned_verts%start_block, owned_verts%end_block
    ! sync the results, can be done at the end
    DO je=1,patch%cell_type
      CALL sync_patch_array(SYNC_V, patch, ocean_coeff%rot_coeff(:,:,:, je))
    ENDDO
    !-------------------------------------------------------------


    !-------------------------------------------------------------
    !3) Handle scalar product coefficients
    !3.1) Edge to cell coefficient
    ! This is done in step 1
!     DO jk=1,n_zlev
!       DO jb = i_startblk_c, i_endblk_c
!         CALL get_indices_c(patch, jb, i_startblk_c, i_endblk_c,      &
!           & i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
!         DO jc = i_startidx_c, i_endidx_c
!           DO je = 1, patch%cells%num_edges(jc,jb)
!
!             ile = patch%cells%edge_idx(jc,jb,je)
!             ibe = patch%cells%edge_blk(jc,jb,je)
!             IF ( v_base%lsm_oce_e(ile,jk,ibe) /= sea) THEN
!               ocean_coeff%edge2cell_coeff_cc(jc,jk,jb,je)%x(1:3) = 0.0_wp
!             ENDIF
!           ENDDO
!         ENDDO
!       END DO
!     END DO

    !-------------------------------------------------------------
    !The dynamical changing coefficient for the surface layer
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c

        DO je = 1, patch%cells%num_edges(jc,jb)
          ocean_coeff%edge2cell_coeff_cc_dyn(jc,1,jb,je)%x = &
            ocean_coeff%edge2cell_coeff_cc(jc,1,jb,je)%x 
        ENDDO ! je = 1, patch%cells%num_edges(jc,jb)

      END DO ! jc = i_startidx_c, i_endidx_c
    END DO ! jb = all_cells%start_block, all_cells%end_block
    ! these are computed an all cells, thus no sync is required
    ! we sync only in p_test_run for checking
!    IF (p_test_run) THEN
      DO je=1,patch%cell_type
          CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc_dyn(:,1,:,je)%x(1))
          CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc_dyn(:,1,:,je)%x(2))
          CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc_dyn(:,1,:,je)%x(3))
      ENDDO
!    ENDIF
    !-------------------------------------------------------------

    !-------------------------------------------------------------
    !3.3) Edge to vert coefficient
    DO jb = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, jb, i_startidx_v, i_endidx_v)
      DO jk = 1, n_zlev
        DO jv = i_startidx_v, i_endidx_v

          DO je = 1, patch%verts%num_edges(jv,jb)
            ile = patch%verts%edge_idx(jv,jb,je)
            ibe = patch%verts%edge_blk(jv,jb,je)
            IF ( v_base%lsm_oce_e(ile,jk,ibe) /= sea) THEN
              ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,je)%x(1:3) = 0.0_wp
              ocean_coeff%variable_dual_vol_norm(jv,jk,jb,je)=0.0_wp
            ENDIF
          ENDDO ! je = 1, patch%verts%num_edges(jv,jb)

        ENDDO ! jv = i_startidx_v, i_endidx_v
      END DO ! jk = 1, n_zlev
    END DO ! jb = owned_verts%start_block, owned_verts%end_block
    ! sync the result
    DO je=1,6
    ! these will be synced in the next loop
!       CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,:,:, je)%x(1))
!       CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,:,:, je)%x(2))
!       CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,:,:, je)%x(3))
      CALL sync_patch_array(SYNC_V, patch, ocean_coeff%variable_dual_vol_norm(:,:,:, je))
    ENDDO
    !-------------------------------------------------------------

    !-------------------------------------------------------------
    !Merge dual area calculation with coefficients
    ! note: i_v_ctr has been calculated on the owned_verts
    !       it does not nedd to be synced
    DO jb = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, jb, i_startidx_v, i_endidx_v)
      DO jk = 1, n_zlev
        DO jv = i_startidx_v, i_endidx_v

          zarea_fraction   = 0.0_wp
          z_area_scaled    = 0.0_wp
          !IF ( ocean_coeff%bnd_edges_per_vertex(jv,jk,jb) == 0 ) THEN
          IF ( i_v_ctr(jv,jk,jb) == patch%verts%num_edges(jv,jb) ) THEN
            z_area_scaled = patch%verts%dual_area(jv,jb)/(re*re)!SUM(ocean_coeff%variable_dual_vol_norm(jv,jk,jb,:))

            !Final coefficient calculation
            DO jev = 1, patch%verts%num_edges(jv,jb)

              IF(z_area_scaled/=0.0_wp)THEN
                ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)&
                  & =ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)!/z_area_scaled
              ELSE
                ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)=0.0_wp
              ENDIF
            END DO

            !ELSEIF(ocean_coeff%bnd_edges_per_vertex(jv,jk,jb)/=0)THEN!boundary edges are involved
          ELSEIF ( i_v_ctr(jv,jk,jb) /= 0 ) THEN

            !Modified area calculation
            vertex_cc = patch%verts%cartesian(jv,jb)
            DO jev = 1, patch%verts%num_edges(jv,jb)
              ! get line and block indices of edge jev around vertex jv
              ile = patch%verts%edge_idx(jv,jb,jev)
              ibe = patch%verts%edge_blk(jv,jb,jev)
              !get neighbor cells
              icell_idx_1 = patch%edges%cell_idx(ile,ibe,1)
              icell_idx_2 = patch%edges%cell_idx(ile,ibe,2)
              icell_blk_1 = patch%edges%cell_blk(ile,ibe,1)
              icell_blk_2 = patch%edges%cell_blk(ile,ibe,2)
              cell1_cc    = gc2cc(patch%cells%center(icell_idx_1,icell_blk_1))
              cell2_cc    = gc2cc(patch%cells%center(icell_idx_2,icell_blk_2))

              !Check, if edge is sea or boundary edge and take care of dummy edge
              !edge with indices ile, ibe is sea edge
              !Add up for wet dual area.
              IF ( v_base%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN
                zarea_fraction = zarea_fraction  &
                  & + triangle_area(cell1_cc, vertex_cc, cell2_cc)
                ! edge with indices ile, ibe is boundary edge
              ELSE IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
                zarea_fraction = zarea_fraction  &
                  & + 0.5_wp*triangle_area(cell1_cc, vertex_cc, cell2_cc)
              END IF
            END DO
            ! no division by zero
            !IF (zarea_fraction /= 0.0_wp) THEN
            z_area_scaled   = patch%verts%dual_area(jv,jb)/(re*re)!zarea_fraction !SUM(ocean_coeff%variable_dual_vol_norm(jv,jk,jb,:))
            !ENDIF
            !Final coefficient calculation
            DO jev = 1, patch%verts%num_edges(jv,jb)
              IF(z_area_scaled/=0.0_wp)THEN
                ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)&
                  & =ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)!/z_area_scaled
              ELSE
                ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)=0.0_wp
              ENDIF
            END DO

          ENDIF !( i_v_ctr(jv,jk,jb) == patch%verts%num_edges(jv,jb) )

        ENDDO
      END DO
    END DO
    ! sync the result
    DO jev=1,6
      DO jk = 1, n_zlev
        CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,jk,:, jev)%x(1))
        CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,jk,:, jev)%x(2))
        CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,jk,:, jev)%x(3))
      ENDDO
    ENDDO
    !-------------------------------------------------------------

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE par_apply_boundary2coeffs
  !--------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Initialize expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !!
  SUBROUTINE init_operator_coeff( ptr_patch, ptr_coeff)
    !
    TYPE(t_patch),      INTENT(inout)     :: ptr_patch
    TYPE(t_operator_coeff), INTENT(inout) :: ptr_coeff
    !-----------------------------------------------------------------------

    CALL init_scalar_product_oce_3d( ptr_patch, ptr_coeff)
    CALL init_geo_factors_oce_3d( ptr_patch, ptr_coeff )

    CALL apply_boundary2coeffs(ptr_patch, ptr_coeff)

  END SUBROUTINE init_operator_coeff
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Initialize expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !!
  SUBROUTINE apply_boundary2coeffs( ptr_patch, ptr_coeff)
    ! !
    TYPE(t_patch),      INTENT(inout)     :: ptr_patch
    TYPE(t_operator_coeff), INTENT(inout) :: ptr_coeff

    !Local variables
    INTEGER :: jk, jc, jb, je, ibe, ile, jev, jv
    INTEGER :: rl_start_e, rl_end_e
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: rl_start_c, rl_end_c
    INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
    INTEGER :: rl_start_v != 2
    INTEGER :: rl_end_v
    INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v

    INTEGER :: i_v_ctr(nproma,n_zlev,ptr_patch%nblks_v)
    INTEGER :: i_v_bnd_edge_ctr(nproma,n_zlev,ptr_patch%nblks_v)
    INTEGER :: ibnd_edge_idx(4), ibnd_edge_blk(4)  !maximal 4 boundary edges in a dual loop.
    INTEGER :: i_edge_idx(4)
    REAL(wp) :: z_orientation(4),z_area_scaled, zarea_fraction
    INTEGER :: icell_idx_1, icell_blk_1
    INTEGER :: icell_idx_2, icell_blk_2
    TYPE(t_cartesian_coordinates) :: cell1_cc, cell2_cc, vertex_cc

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:apply_boundary2coeffs')

    !-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'start')

    rl_start_c   = 1
    rl_end_c     = min_rlcell_int
    i_startblk_c = ptr_patch%cells%start_blk(rl_start_c,1)
    i_endblk_c   = ptr_patch%cells%end_blk(rl_end_c,1)
    rl_start_e   = 1
    rl_end_e     = min_rledge ! Loop over the whole local domain
    i_startblk_e = ptr_patch%edges%start_blk(rl_start_e,1)
    i_endblk_e   = ptr_patch%edges%end_blk(rl_end_e,1)
    rl_start_v   = 1
    rl_end_v     = min_rlvert
    i_startblk_v = ptr_patch%verts%start_blk(rl_start_v,1)
    i_endblk_v   = ptr_patch%verts%end_blk(rl_end_v,1)

    i_v_ctr(:,:,:)          = 0
    i_v_bnd_edge_ctr(:,:,:) = 0

    !1) Set coefficients for div and grad to zero at boundary edges
    DO jk=1,n_zlev
      DO jb = i_startblk_c, i_endblk_c
        CALL get_indices_c(ptr_patch, jb, i_startblk_c, i_endblk_c,      &
          & i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
          DO je = 1, ptr_patch%cells%num_edges(jc,jb)

            ile = ptr_patch%cells%edge_idx(jc,jb,je)
            ibe = ptr_patch%cells%edge_blk(jc,jb,je)

            IF ( v_base%lsm_oce_e(ile,jk,ibe) /= sea) THEN
              ptr_coeff%div_coeff(jc,jk,jb,je) = 0.0_wp
              ptr_coeff%grad_coeff(ile,jk,ibe) = 0.0_wp
            ENDIF
          ENDDO
          !write(1234,*)'div coeff 3D',jk,jc,jb,ptr_intp%div_coeff(jc,jk,jb,:)
        ENDDO
      END DO
    END DO

    !2) prepare coefficients for rot at boundary edges
    DO jb = i_startblk_v, i_endblk_v

      CALL get_indices_v(ptr_patch, jb, i_startblk_v, i_endblk_v, &
        & i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
      DO jk = 1, n_zlev
        !!$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(PRIVATE)  &
        !!$OMP   SHARED(u_vec_e,v_vec_e,ptr_patch,rot_vec_v,jb) FIRSTPRIVATE(jk)
        DO jv = i_startidx_v, i_endidx_v

          ibnd_edge_idx(1:4)      = 0
          ibnd_edge_blk(1:4)      = 0
          i_edge_idx(1:4)         = 0
          z_orientation(1:4)      = 0.0_wp

          DO jev = 1, ptr_patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = ptr_patch%verts%edge_idx(jv,jb,jev)
            ibe = ptr_patch%verts%edge_blk(jv,jb,jev)
            !Check, if edge is sea or boundary edge and take care of dummy edge
            ! edge with indices ile, ibe is sea edge
            ! edge with indices ile, ibe is boundary edge

            IF ( v_base%lsm_oce_e(ile,jk,ibe) == sea) THEN
              i_v_ctr(jv,jk,jb)=i_v_ctr(jv,jk,jb)+1
            ELSEIF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN

              !increase boundary edge counter
              i_v_bnd_edge_ctr(jv,jk,jb)=i_v_bnd_edge_ctr(jv,jk,jb)+1

              ptr_coeff%bnd_edges_per_vertex(jv,jk,jb) &
                & = ptr_coeff%bnd_edges_per_vertex(jv,jk,jb) +1

              !Store actual boundary edge indices
              IF(i_v_bnd_edge_ctr(jv,jk,jb)==1)THEN
                ibnd_edge_idx(1) = ile
                ibnd_edge_blk(1) = ibe
                z_orientation(1) = ptr_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(1)    = jev

                ptr_coeff%bnd_edge_idx(jv,jk,jb,1)= ile
                ptr_coeff%bnd_edge_blk(jv,jk,jb,1)= ibe
                ptr_coeff%orientation(jv,jk,jb,1) = ptr_patch%verts%edge_orientation(jv,jb,jev)
                ptr_coeff%edge_idx(jv,jk,jb,1)    = jev

              ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==2)THEN
                ibnd_edge_idx(2) = ile
                ibnd_edge_blk(2) = ibe
                z_orientation(2) = ptr_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(2)    = jev

                ptr_coeff%bnd_edge_idx(jv,jk,jb,2)= ile
                ptr_coeff%bnd_edge_blk(jv,jk,jb,2)= ibe
                ptr_coeff%orientation(jv,jk,jb,2) = ptr_patch%verts%edge_orientation(jv,jb,jev)
                ptr_coeff%edge_idx(jv,jk,jb,2)    = jev

              ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==3)THEN
                ibnd_edge_idx(3) = ile
                ibnd_edge_blk(3) = ibe
                z_orientation(3) = ptr_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(3)    = jev

                ptr_coeff%bnd_edge_idx(jv,jk,jb,3)= ile
                ptr_coeff%bnd_edge_blk(jv,jk,jb,3)= ibe
                ptr_coeff%orientation(jv,jk,jb,3) = ptr_patch%verts%edge_orientation(jv,jb,jev)
                ptr_coeff%edge_idx(jv,jk,jb,3)    = jev

              ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==4)THEN
                ibnd_edge_idx(4) = ile
                ibnd_edge_blk(4) = ibe
                z_orientation(4) = ptr_patch%verts%edge_orientation(jv,jb,jev)
                i_edge_idx(4)    = jev

                ptr_coeff%bnd_edge_idx(jv,jk,jb,4)= ile
                ptr_coeff%bnd_edge_blk(jv,jk,jb,4)= ibe
                ptr_coeff%orientation(jv,jk,jb,4) = ptr_patch%verts%edge_orientation(jv,jb,jev)
                ptr_coeff%edge_idx(jv,jk,jb,4)    = jev

              ELSE
                !maximal 4 boundary edges per dual loop are allowed: somethings wrong withe the grid
                CALL message (TRIM('sbr nonlinear Coriolis'), &
                  & 'more than 4 boundary edges per dual loop: something is wrong with the grid')
                CALL finish ('TRIM(sbr nonlinear Coriolis)','Grid-boundary error !!')
              ENDIF
            END IF
          END DO

          IF(i_v_bnd_edge_ctr(jv,jk,jb)==2)THEN
            ptr_coeff%rot_coeff(jv,jk,jb,i_edge_idx(1) )=&
              & 0.5_wp*z_orientation(1)&
              & *ptr_patch%edges%primal_edge_length(ibnd_edge_idx(1),ibnd_edge_blk(1))

            ptr_coeff%rot_coeff(jv,jk,jb,i_edge_idx(2))=&
              & 0.5_wp*z_orientation(2)&
              & *ptr_patch%edges%primal_edge_length(ibnd_edge_idx(2),ibnd_edge_blk(2))

          ELSEIF(i_v_bnd_edge_ctr(jv,jk,jb)==4)THEN

            ptr_coeff%rot_coeff(jv,jk,jb,i_edge_idx(1))=&
              & 0.5_wp*z_orientation(1)&
              & *ptr_patch%edges%primal_edge_length(ibnd_edge_idx(1),ibnd_edge_blk(1))

            ptr_coeff%rot_coeff(jv,jk,jb,i_edge_idx(2))=&
              & 0.5_wp*z_orientation(2)&
              & *ptr_patch%edges%primal_edge_length(ibnd_edge_idx(2),ibnd_edge_blk(2))

            ptr_coeff%rot_coeff(jv,jk,jb,i_edge_idx(3))=&
              & 0.5_wp*z_orientation(3)&
              & *ptr_patch%edges%primal_edge_length(ibnd_edge_idx(3),ibnd_edge_blk(3))

            ptr_coeff%rot_coeff(jv,jk,jb,i_edge_idx(4))=&
              & 0.5_wp*z_orientation(4)&
              & *ptr_patch%edges%primal_edge_length(ibnd_edge_idx(4),ibnd_edge_blk(4))
          ENDIF
          !ptr_coeff%rot_coeff(jv,jk,jb,:)=ptr_coeff%rot_coeff(jv,jk,jb,:)/ptr_patch%verts%dual_area(jv,jb)
        END DO
        !!$OMP END PARALLEL DO
      END DO
    END DO
    ! DO jb = i_startblk_v, i_endblk_v
    !   CALL get_indices_v(ptr_patch, jb, i_startblk_v, i_endblk_v, &
    !                      i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
    !   DO jk = 1, n_zlev
    !     DO jv = i_startidx_v, i_endidx_v
    ! !IF(rot_vec_v(jv,jk,jb)/=0.0_wp)THEN
    ! write(1234567,*)'rot 3D COEFF:',jk,jv,jb,i_v_bnd_edge_ctr(jv,jk,jb),i_v_ctr(jv,jk,jb)
    ! !ENDIF
    !     END DO
    !   END DO
    ! END DO


    !3) Handle scalar product coefficients
    !3.1) Edge to cell coefficient
    DO jk=1,n_zlev
      DO jb = i_startblk_c, i_endblk_c
        CALL get_indices_c(ptr_patch, jb, i_startblk_c, i_endblk_c,      &
          & i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
          DO je = 1, ptr_patch%cells%num_edges(jc,jb)

            ile = ptr_patch%cells%edge_idx(jc,jb,je)
            ibe = ptr_patch%cells%edge_blk(jc,jb,je)
            IF ( v_base%lsm_oce_e(ile,jk,ibe) /= sea) THEN
              ptr_coeff%edge2cell_coeff_cc(jc,jk,jb,je)%x(1:3) = 0.0_wp
            ENDIF
          ENDDO
        ENDDO
      END DO
    END DO
    !The dynamical changing coefficient for the surface layer
    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c(ptr_patch, jb, i_startblk_c, i_endblk_c,      &
        & i_startidx_c, i_endidx_c, rl_start_c, rl_end_c)
      DO jc = i_startidx_c, i_endidx_c
        DO je = 1, ptr_patch%cells%num_edges(jc,jb)

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)
          IF ( v_base%lsm_oce_e(ile,1,ibe) /= sea) THEN
            ptr_coeff%edge2cell_coeff_cc_dyn(jc,1,jb,je)%x(1:3) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
    END DO

    DO jk=1,n_zlev
      DO jb = i_startblk_e, i_endblk_e
        CALL get_indices_e(ptr_patch, jb, i_startblk_e, i_endblk_e,      &
          & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
        DO je = i_startidx_e, i_endidx_e

          IF(v_base%lsm_oce_e(je,jk,jb) /= sea) THEN
            ptr_coeff%edge2cell_coeff_cc_t(je,jk,jb,1)%x(1:3) = 0.0_wp
            ptr_coeff%edge2cell_coeff_cc_t(je,jk,jb,2)%x(1:3) = 0.0_wp
          ENDIF
        ENDDO
      END DO
    END DO

    !3.3) Edge to vert coefficient
    DO jk=1,n_zlev
      DO jb = i_startblk_v, i_endblk_v
        CALL get_indices_v(ptr_patch, jb, i_startblk_v, i_endblk_v,      &
          & i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
        DO jv = i_startidx_v, i_endidx_v
          DO je = 1, ptr_patch%verts%num_edges(jv,jb)
            ile = ptr_patch%verts%edge_idx(jv,jb,je)
            ibe = ptr_patch%verts%edge_blk(jv,jb,je)
            IF ( v_base%lsm_oce_e(ile,jk,ibe) /= sea) THEN
              ptr_coeff%edge2vert_coeff_cc(jv,jk,jb,je)%x(1:3) = 0.0_wp
              ptr_coeff%variable_dual_vol_norm(jv,jk,jb,je)=0.0_wp
            ENDIF
          ENDDO
        ENDDO
      END DO
    END DO

    !Merge dual area calculation with coefficients
    DO jk=1,n_zlev
      DO jb = i_startblk_v, i_endblk_v
        CALL get_indices_v(ptr_patch, jb, i_startblk_v, i_endblk_v,      &
          & i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
        DO jv = i_startidx_v, i_endidx_v

          zarea_fraction   = 0.0_wp
          z_area_scaled    = 0.0_wp
          !IF ( ptr_coeff%bnd_edges_per_vertex(jv,jk,jb) == 0 ) THEN
          IF ( i_v_ctr(jv,jk,jb) == ptr_patch%verts%num_edges(jv,jb) ) THEN
            z_area_scaled = ptr_patch%verts%dual_area(jv,jb)/(re*re)!SUM(ptr_coeff%variable_dual_vol_norm(jv,jk,jb,:))

            !Final coefficient calculation
            DO jev = 1, ptr_patch%verts%num_edges(jv,jb)

              IF(z_area_scaled/=0.0_wp)THEN
                ptr_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)&
                  & =ptr_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)!/z_area_scaled
              ELSE
                ptr_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)=0.0_wp
              ENDIF
            END DO

            !ELSEIF(ptr_coeff%bnd_edges_per_vertex(jv,jk,jb)/=0)THEN!boundary edges are involved
          ELSEIF ( i_v_ctr(jv,jk,jb) /= 0 ) THEN

            !Modified area calculation
            vertex_cc = gc2cc(ptr_patch%verts%vertex(jv,jb))
            DO jev = 1, ptr_patch%verts%num_edges(jv,jb)
              ! get line and block indices of edge jev around vertex jv
              ile = ptr_patch%verts%edge_idx(jv,jb,jev)
              ibe = ptr_patch%verts%edge_blk(jv,jb,jev)
              !get neighbor cells
              icell_idx_1 = ptr_patch%edges%cell_idx(ile,ibe,1)
              icell_idx_2 = ptr_patch%edges%cell_idx(ile,ibe,2)
              icell_blk_1 = ptr_patch%edges%cell_blk(ile,ibe,1)
              icell_blk_2 = ptr_patch%edges%cell_blk(ile,ibe,2)
              cell1_cc    = gc2cc(ptr_patch%cells%center(icell_idx_1,icell_blk_1))
              cell2_cc    = gc2cc(ptr_patch%cells%center(icell_idx_2,icell_blk_2))

              !Check, if edge is sea or boundary edge and take care of dummy edge
              !edge with indices ile, ibe is sea edge
              !Add up for wet dual area.
              IF ( v_base%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN
                zarea_fraction = zarea_fraction  &
                  & + triangle_area(cell1_cc, vertex_cc, cell2_cc)
                ! edge with indices ile, ibe is boundary edge
              ELSE IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
                zarea_fraction = zarea_fraction  &
                  & + 0.5_wp*triangle_area(cell1_cc, vertex_cc, cell2_cc)
              END IF
            END DO
            ! no division by zero
            !IF (zarea_fraction /= 0.0_wp) THEN
            z_area_scaled   = ptr_patch%verts%dual_area(jv,jb)/(re*re)!zarea_fraction !SUM(ptr_coeff%variable_dual_vol_norm(jv,jk,jb,:))
            !ENDIF
            !Final coefficient calculation
            DO jev = 1, ptr_patch%verts%num_edges(jv,jb)
              IF(z_area_scaled/=0.0_wp)THEN
                ptr_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)&
                  & =ptr_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)!/z_area_scaled
              ELSE
                ptr_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)=0.0_wp
              ENDIF
            END DO
          ENDIF
        ENDDO
      END DO
    END DO


    !3.4) Vert to edge coefficient
    DO jk=1,n_zlev
      DO jb = i_startblk_e, i_endblk_e
        CALL get_indices_e(ptr_patch, jb, i_startblk_e, i_endblk_e,      &
          & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
        DO je = i_startidx_e, i_endidx_e

          IF(v_base%lsm_oce_e(je,jk,jb) /= sea) THEN
            ptr_coeff%edge2vert_coeff_cc_t(je,jk,jb,1)%x(1:3) = 0.0_wp
            ptr_coeff%edge2vert_coeff_cc_t(je,jk,jb,2)%x(1:3) = 0.0_wp
          ENDIF
        ENDDO
      END DO
    END DO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE apply_boundary2coeffs
  !--------------------------------------------------------------------------------------

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
  SUBROUTINE init_scalar_product_oce_3d( ptr_patch, ptr_intp)

    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch
    TYPE(t_operator_coeff),INTENT(inout) :: ptr_intp

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:init_scalar_product_oce_3d')

    INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER, PARAMETER :: no_vert_edges = 6
    INTEGER :: jb, je, jv, ie, ie_1, ie_2, icc,jk, jc
    INTEGER :: il_e,ib_e,k
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: il_v1, il_v2, ib_v1, ib_v2

    INTEGER :: iil_c1(no_cell_edges), iil_c2(no_cell_edges)
    INTEGER :: iib_c1(no_cell_edges), iib_c2(no_cell_edges)

    INTEGER :: jil_c1, jib_c1,jil_c2, jib_c2
    INTEGER :: rl_start_e, rl_end_e
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

    REAL(wp) :: z_tmp
    REAL(wp) :: norm_c1_c2, norm_v1_v2, norm
    REAL(wp) :: dual_edge_length(no_vert_edges)
    REAL(wp) :: vert_edge_dist(no_vert_edges,2)
    REAL(wp) :: vert_dual_mid_dist(no_vert_edges,2)
    REAL(wp) :: vert_edge_distance, vert_dual_mid_distance

    TYPE(t_geographical_coordinates) :: gc_mid_dual_edge(no_vert_edges)
    TYPE(t_geographical_coordinates) :: gc1,gc2

    TYPE(t_cartesian_coordinates)    :: cc_dual_edge(no_vert_edges), cc_edge(no_cell_edges)
    TYPE(t_cartesian_coordinates)    :: xx1,xx2
    TYPE(t_cartesian_coordinates)    :: vert1_midedge_cc(nproma,ptr_patch%nblks_v,no_vert_edges)
    TYPE(t_cartesian_coordinates)    :: vert2_midedge_cc(nproma,ptr_patch%nblks_v,no_vert_edges)
    TYPE(t_cartesian_coordinates)    :: cell2cell_cc
    TYPE(t_cartesian_coordinates)    :: cc_e0, cc_c1,cc_c2,cc_v0
    TYPE(t_cartesian_coordinates)    :: cv_c1_e0, cv_c2_e0, cv_c1_c2
    TYPE(t_cartesian_coordinates)    :: cc_mid_dual_edge(no_vert_edges)
    TYPE(t_cartesian_coordinates)    :: recon_vec_cc
    TYPE(t_cartesian_coordinates)    :: z_vec_c1(no_cell_edges),z_vec_c2(no_cell_edges)
    TYPE(t_cartesian_coordinates)    :: recon_vec_cc_v1(no_vert_edges)
    TYPE(t_cartesian_coordinates)    :: recon_vec_cc_v2(no_vert_edges)

    REAL(wp) :: z_edge_length(no_cell_edges)
    REAL(wp) :: z_cell_edge_dist_c1(no_cell_edges,2),z_cell_edge_dist_c2(no_cell_edges,2)
    REAL(wp) :: z_y

    REAL(wp) :: z_sync_c(nproma,n_zlev, ptr_patch%nblks_c)
    REAL(wp) :: z_sync_e(nproma,n_zlev, ptr_patch%nblks_e)
    REAL(wp) :: z_sync_v(nproma,n_zlev, ptr_patch%nblks_v)

    LOGICAL, PARAMETER :: mid_point_dual_edge = .TRUE. !Please do not change this unless
    !you are sure, you know what you do.
    LOGICAL, PARAMETER :: larc_length = .FALSE.
    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')

    rl_start     = 1
    rl_end       = min_rlcell_int
    i_startblk   = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk     = ptr_patch%cells%end_blk(rl_end,1)

    rl_start_e   = 1
    rl_end_e     = min_rledge ! Loop over the whole local domain
    i_startblk_e = ptr_patch%edges%start_blk(rl_start_e,1)
    i_endblk_e   = ptr_patch%edges%end_blk(rl_end_e,1)


    ptr_intp%fixed_vol_norm(:,:,:) = 0._wp

    !-----------------------------------------------------------------------
    !STEP 1: edge2cell and cell2edge coefficients
    DO jk=1,n_zlev
      edge_blk_loop_primal: DO jb = i_startblk_e, i_endblk_e

        CALL get_indices_e(ptr_patch, jb,&
          & i_startblk_e, i_endblk_e,&
          & i_startidx_e, i_endidx_e,&
          & rl_start_e, rl_end_e)

        edge_idx_loop_primal: DO je =  i_startidx_e, i_endidx_e

          !Get indices of two adjacent triangles
          il_c1 = ptr_patch%edges%cell_idx(je,jb,1)
          ib_c1 = ptr_patch%edges%cell_blk(je,jb,1)
          il_c2 = ptr_patch%edges%cell_idx(je,jb,2)
          ib_c2 = ptr_patch%edges%cell_blk(je,jb,2)

          ! Go only over edges where at least one neighboring triangle is not in the halo
          IF(il_c1<=0 .OR. il_c2<=0) CYCLE
          IF(.NOT.ptr_patch%cells%owner_mask(il_c1,ib_c1) .AND. &
            & .NOT.ptr_patch%cells%owner_mask(il_c2,ib_c2) ) CYCLE

          cc_c1 = gc2cc(ptr_patch%cells%center(il_c1, ib_c1))
          cc_c2 = gc2cc(ptr_patch%cells%center(il_c2, ib_c2))

          z_cell_edge_dist_c1 = 0.0_wp
          z_cell_edge_dist_c2 = 0.0_wp

          !normals in cell 1
          DO ie = 1, no_cell_edges

            !actual edges of cell c1
            iil_c1(ie)  = ptr_patch%cells%edge_idx(il_c1,ib_c1,ie)
            iib_c1(ie)  = ptr_patch%cells%edge_blk(il_c1,ib_c1,ie)

            cc_edge(ie) = gc2cc(ptr_patch%edges%center(iil_c1(ie),iib_c1(ie)))

            !calculate edge length
            !get vertex indices adjacent to actual edge
            il_v1       = ptr_patch%edges%vertex_idx(iil_c1(ie),iib_c1(ie),1)
            ib_v1       = ptr_patch%edges%vertex_blk(iil_c1(ie),iib_c1(ie),1)
            il_v2       = ptr_patch%edges%vertex_idx(iil_c1(ie),iib_c1(ie),2)
            ib_v2       = ptr_patch%edges%vertex_blk(iil_c1(ie),iib_c1(ie),2)

            !get vertex positions
            xx1         = gc2cc(ptr_patch%verts%vertex(il_v1,ib_v1))
            xx2         = gc2cc(ptr_patch%verts%vertex(il_v2,ib_v2))

            IF(larc_length)THEN
              norm              = SQRT(SUM(xx1%x*xx1%x))
              xx1%x             = xx1%x/norm
              norm              = SQRT(SUM(xx2%x*xx2%x))
              xx2%x             = xx2%x/norm
              z_edge_length(ie) = arc_length(xx2,xx1)
            ELSE
              z_edge_length(ie) = SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
            ENDIF

            !calculate cell-edge distance as half of cell-cell distance
            !get cell indices adjacent to actual edge
            jil_c1 = ptr_patch%edges%cell_idx(iil_c1(ie),iib_c1(ie),1)
            jib_c1 = ptr_patch%edges%cell_blk(iil_c1(ie),iib_c1(ie),1)
            jil_c2 = ptr_patch%edges%cell_idx(iil_c1(ie),iib_c1(ie),2)
            jib_c2 = ptr_patch%edges%cell_blk(iil_c1(ie),iib_c1(ie),2)

            !get cell positions
            xx1    = gc2cc(ptr_patch%cells%center(jil_c1,jib_c1))
            xx2    = gc2cc(ptr_patch%cells%center(jil_c2,jib_c2))

            IF(jil_c1==il_c1.AND.jib_c1==ib_c1)THEN
              k=1
            ELSEIF(jil_c2==il_c1.AND.jib_c2==ib_c1)THEN
              k=2
            ENDIF

            IF(larc_length)THEN
              norm                      = SQRT(SUM(xx1%x*xx1%x))
              xx1%x                     = xx1%x/norm
              norm                      = SQRT(SUM(xx2%x*xx2%x))
              xx2%x                     = xx2%x/norm
              norm                      = SQRT(SUM(cc_edge(ie)%x*cc_edge(ie)%x))
              cc_edge(ie)%x             = cc_edge(ie)%x/norm
              z_cell_edge_dist_c1(ie,1) = arc_length(cc_edge(ie),xx1)
              z_cell_edge_dist_c1(ie,2) = arc_length(cc_edge(ie),xx2)
            ELSE
              z_cell_edge_dist_c1(ie,1) = SQRT(SUM((cc_edge(ie)%x-xx1%x)*(cc_edge(ie)%x-xx1%x)))
              z_cell_edge_dist_c1(ie,2) = SQRT(SUM((cc_edge(ie)%x-xx2%x)*(cc_edge(ie)%x-xx2%x)))
            ENDIF
            ptr_intp%dist_cell2edge(iil_c1(ie),jk,iib_c1(ie),1) = z_cell_edge_dist_c1(ie,1)
            ptr_intp%dist_cell2edge(iil_c1(ie),jk,iib_c1(ie),2) = z_cell_edge_dist_c1(ie,2)

            z_vec_c1(ie)%x = cc_edge(ie)%x - cc_c1%x     !ptr_patch%edges%primal_cart_normal(iil_c1(ie),iib_c1(ie))
            norm           = SQRT(SUM( z_vec_c1(ie)%x* z_vec_c1(ie)%x))
            !write(*,*)'NORM:',norm !TODOram

            ptr_intp%edge2cell_coeff_cc(il_c1,jk,ib_c1,ie)%x = &
              & z_vec_c1(ie)%x*ptr_patch%cells%edge_orientation(il_c1,ib_c1,ie)*z_edge_length(ie)

            ptr_intp%fixed_vol_norm(il_c1,jk,ib_c1) = ptr_intp%fixed_vol_norm(il_c1,jk,ib_c1)&
              & +  0.5_wp*norm*z_edge_length(ie)
            ptr_intp%variable_vol_norm(il_c1,jk,ib_c1,ie) = 0.5_wp*norm*z_edge_length(ie)

            !write(*,*)'edge length   :',z_edge_length(ie),ptr_patch%edges%primal_edge_length(iil_c1(ie),iib_c1(ie))/re
            !write(*,*)'cell-edge dist:', z_cell_edge_dist_c1(ie,k),ptr_patch%edges%edge_cell_length(iil_c1(ie),iib_c1(ie),k)/re
          END DO

          !normals in cell 2
          DO ie = 1, no_cell_edges

            !actual edges of cell c2
            iil_c2(ie)  = ptr_patch%cells%edge_idx(il_c2,ib_c2,ie)
            iib_c2(ie)  = ptr_patch%cells%edge_blk(il_c2,ib_c2,ie)


            cc_edge(ie) = gc2cc(ptr_patch%edges%center(iil_c2(ie),iib_c2(ie)))

            !calculate edge length
            !get vertex indices adjacent to actual edge
            il_v1       = ptr_patch%edges%vertex_idx(iil_c2(ie),iib_c2(ie),1)
            ib_v1       = ptr_patch%edges%vertex_blk(iil_c2(ie),iib_c2(ie),1)
            il_v2       = ptr_patch%edges%vertex_idx(iil_c2(ie),iib_c2(ie),2)
            ib_v2       = ptr_patch%edges%vertex_blk(iil_c2(ie),iib_c2(ie),2)

            !get vertex positions
            xx1         = gc2cc(ptr_patch%verts%vertex(il_v1,ib_v1))
            xx2         = gc2cc(ptr_patch%verts%vertex(il_v2,ib_v2))

            IF(larc_length)THEN
              norm              = SQRT(SUM(xx1%x*xx1%x))
              xx1%x             = xx1%x/norm

              norm              = SQRT(SUM(xx2%x*xx2%x))
              xx2%x             = xx2%x/norm

              z_edge_length(ie) = arc_length(xx2,xx1)
              !z_edge_length(ie) = ptr_patch%edges%primal_edge_length(iil_c2(ie),iib_c2(ie))/re
              !write(*,*)'arc length',arc_length(xx2,xx1),z_edge_length(ie),SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
            ELSE
              z_edge_length(ie) = SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
            ENDIF
            !calculate cell-edge distance as half of cell-cell distance
            !get cell indices adjacent to actual edge
            jil_c1 = ptr_patch%edges%cell_idx(iil_c2(ie),iib_c2(ie),1)
            jib_c1 = ptr_patch%edges%cell_blk(iil_c2(ie),iib_c2(ie),1)
            jil_c2 = ptr_patch%edges%cell_idx(iil_c2(ie),iib_c2(ie),2)
            jib_c2 = ptr_patch%edges%cell_blk(iil_c2(ie),iib_c2(ie),2)

            IF (jil_c2 < 0) THEN
              WRITE(0,*)'ptr_patch%edges%cell_idx:',ptr_patch%edges%cell_idx
            ENDIF
            !get cell positions
            xx1 = gc2cc(ptr_patch%cells%center(jil_c1,jib_c1))
            xx2 = gc2cc(ptr_patch%cells%center(jil_c2,jib_c2))

            IF(jil_c1==il_c2.AND.jib_c1==ib_c2)THEN
              k=1
            ELSEIF(jil_c2==il_c2.AND.jib_c2==ib_c2)THEN
              k=2
            ENDIF

            IF(larc_length)THEN
              norm                      = SQRT(SUM(xx1%x*xx1%x))
              xx1%x                     = xx1%x/norm
              norm                      = SQRT(SUM(xx2%x*xx2%x))
              xx2%x                     = xx2%x/norm
              norm                      = SQRT(SUM(cc_edge(ie)%x*cc_edge(ie)%x))
              cc_edge(ie)%x             = cc_edge(ie)%x/norm
              z_cell_edge_dist_c2(ie,1) = arc_length(cc_edge(ie),xx1)
              z_cell_edge_dist_c2(ie,2) = arc_length(cc_edge(ie),xx2)
            ELSE
              z_cell_edge_dist_c2(ie,1) = SQRT(SUM((cc_edge(ie)%x-xx1%x)*(cc_edge(ie)%x-xx1%x)))
              z_cell_edge_dist_c2(ie,2) = SQRT(SUM((cc_edge(ie)%x-xx2%x)*(cc_edge(ie)%x-xx2%x)))
            ENDIF
            ptr_intp%dist_cell2edge(iil_c2(ie),jk,iib_c2(ie),1) = z_cell_edge_dist_c2(ie,1)
            ptr_intp%dist_cell2edge(iil_c2(ie),jk,iib_c2(ie),2) = z_cell_edge_dist_c2(ie,2)

            z_vec_c2(ie)%x = cc_edge(ie)%x - cc_c2%x  !ptr_patch%edges%primal_cart_normal(iil_c2(ie),iib_c2(ie))
            norm           = SQRT(SUM( z_vec_c2(ie)%x* z_vec_c2(ie)%x))

            ptr_intp%edge2cell_coeff_cc(il_c2,jk,ib_c2,ie)%x = &
              &  z_vec_c2(ie)%x * ptr_patch%cells%edge_orientation(il_c2,ib_c2,ie) * &
              &  z_edge_length(ie)

            ptr_intp%fixed_vol_norm(il_c2,jk,ib_c2) = ptr_intp%fixed_vol_norm(il_c2,jk,ib_c2)&
              & + 0.5_wp*norm*z_edge_length(ie)
            ptr_intp%variable_vol_norm(il_c2,jk,ib_c2,ie) = 0.5_wp*norm*z_edge_length(ie)

          END DO
        END DO edge_idx_loop_primal
      END DO edge_blk_loop_primal
    END DO
    !In the edge loop above each triangle is visisted three times. Since the "fixed_vol_norm" is
    !accumulated we correct its value here:
    ptr_intp%fixed_vol_norm = ptr_intp%fixed_vol_norm/3.0_wp


    !Assign values to dynamical coefficients for surface layer
    edge_blk_loop_dyn: DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb,&
        & i_startblk, i_endblk,&
        & i_startidx, i_endidx,&
        & rl_start, rl_end)

      edge_idx_loop_dyn: DO jc =  i_startidx, i_endidx
        DO ie = 1, no_cell_edges
          ptr_intp%edge2cell_coeff_cc_dyn(jc,1,jb,ie)%x=ptr_intp%edge2cell_coeff_cc(jc,1,jb,ie)%x
        END DO
      END DO edge_idx_loop_dyn
    END DO edge_blk_loop_dyn


    !commented put for testing--------------------------------

    !    !merge fixed volume and edge2cell coeff
    !     DO jk=1,n_zlev
    !       DO jb = i_startblk, i_endblk
    !
    !       CALL get_indices_c(ptr_patch, jb,&
    !                        & i_startblk, i_endblk,&
    !                        & i_startidx, i_endidx,&
    !                        & rl_start, rl_end)
    !
    !         DO jc =  i_startidx, i_endidx
    !           DO ie=1,no_cell_edges
    !             ptr_intp%edge2cell_coeff_cc(jc,jk,jb,ie)%x&
    !             & = ptr_intp%edge2cell_coeff_cc(jc,jk,jb,ie)%x/ptr_intp%fixed_vol_norm(jc,jk,jb)
    !           END DO
    !         END DO
    !       END DO
    !     END DO

    rl_start   = 1
    rl_end     = min_rledge_int
    i_startblk = ptr_patch%edges%start_blk(rl_start,1)
    i_endblk   = ptr_patch%edges%end_blk(rl_end,1)

    DO jk=1,n_zlev
      edge_blk_loop_secondary: DO jb = i_startblk, i_endblk
        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk,&
          & i_startidx, i_endidx, rl_start, rl_end)
        edge_idx_loop_secondary: DO je =  i_startidx, i_endidx

          !Get indices of two adjacent triangles
          il_c1      = ptr_patch%edges%cell_idx(je,jb,1)
          ib_c1      = ptr_patch%edges%cell_blk(je,jb,1)
          il_c2      = ptr_patch%edges%cell_idx(je,jb,2)
          ib_c2      = ptr_patch%edges%cell_blk(je,jb,2)

          !cartesian coordinates of edge and neighbor cells on 1-sphere
          cc_e0      = gc2cc(ptr_patch%edges%center(je,jb))
          cc_c1      = gc2cc(ptr_patch%cells%center(il_c1,ib_c1))
          cc_c2      = gc2cc(ptr_patch%cells%center(il_c2,ib_c2))

          !cartesian vectors from:
          !cell 2 to cell 1, cell 1 to edge je and cell 2 to edge je
          cv_c1_c2%x = cc_c1%x - cc_c2%x
          cv_c1_e0%x = cc_e0%x - cc_c1%x
          cv_c2_e0%x = cc_e0%x - cc_c2%x

          IF(larc_length)THEN
            norm        = SQRT(SUM(cc_e0%x*cc_e0%x))
            cc_e0%x     = cc_e0%x/norm
            norm        = SQRT(SUM(cc_c1%x*cc_c1%x))
            cc_c1%x     = cc_c1%x/norm
            norm        = SQRT(SUM(cc_c2%x*cc_c2%x))
            cc_c2%x     = cc_c2%x/norm
            norm_c1_c2  = arc_length(cc_c1, cc_c2)
          ELSE
            norm_c1_c2  = SQRT(SUM(cv_c1_e0%x*cv_c1_e0%x))+SQRT(SUM(cv_c2_e0%x*cv_c2_e0%x)) !SQRT(SUM(cv_c1_c2%x*cv_c1_c2%x))!
          ENDIF

          !Determine which edge of both of the two adjacent cells corresponds to the
          !actual edge "je". This information is used below for the edge-orientation.
          DO ie = 1, no_cell_edges
            IF (ptr_patch%cells%edge_idx(il_c1,ib_c1,ie) == je.AND.&
              & ptr_patch%cells%edge_blk(il_c1,ib_c1,ie) == jb) THEN
              ie_1 = ie
            END IF
            IF (ptr_patch%cells%edge_idx(il_c2,ib_c2,ie) == je.AND.&
              & ptr_patch%cells%edge_blk(il_c2,ib_c2,ie) == jb) THEN
              ie_2 = ie
            END IF
          END DO

          ptr_intp%edge2cell_coeff_cc_t(je,jk,jb,1)%x&
            & = cv_c1_e0%x * ptr_patch%cells%edge_orientation(il_c1,ib_c1,ie_1)/norm_c1_c2

          ptr_intp%edge2cell_coeff_cc_t(je,jk,jb,2)%x&
            & = cv_c2_e0%x * ptr_patch%cells%edge_orientation(il_c2,ib_c2,ie_2)/norm_c1_c2

        END DO edge_idx_loop_secondary
      END DO edge_blk_loop_secondary
    END DO
    !------------------------------------------------------------------------------
    !STEP 2: edge2vert coefficients for dual grid
    !------------------------------------------------------------------------------

    rl_start = 1
    rl_end   = min_rlvert ! Loop over the whole local domain

    i_startblk = ptr_patch%verts%start_blk(rl_start,1)
    i_endblk   = ptr_patch%verts%end_blk(rl_end,1)

    DO jk=1,n_zlev
      vert_blk_loop: DO jb = i_startblk, i_endblk
        CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk,&
          & i_startidx, i_endidx, rl_start, rl_end)
        vert_idx_loop: DO jv =  i_startidx, i_endidx
          ! current number of edges around vertex (5 or 6)
          cc_v0        = gc2cc(ptr_patch%verts%vertex(jv,jb))

          ! Go only over vertices where all edges are in the local domain (but maybe in the halo)
          IF(ANY(ptr_patch%verts%edge_idx(jv,jb,:)<0)) CYCLE

          DO ie = 1, no_vert_edges  ! #slo# it_vertedges ??

            il_e             = ptr_patch%verts%edge_idx(jv,jb,ie)
            ib_e             = ptr_patch%verts%edge_blk(jv,jb,ie)

            ! #slo# - I assume all geographical coordinates are already synchronized

            cc_dual_edge(ie) = gc2cc(ptr_patch%edges%center(il_e,ib_e))
            !Parts of this code parallels the implementation in the grid-generator
            !module "mo_geometry".
            !
            !1) determine normal vector from adjacent cell to adjacent cell
            !   in cartesian coordinate for moved dual cell
            !Get indices of two adjacent triangles
            il_c1            = ptr_patch%edges%cell_idx(il_e,ib_e,1)
            ib_c1            = ptr_patch%edges%cell_blk(il_e,ib_e,1)
            il_c2            = ptr_patch%edges%cell_idx(il_e,ib_e,2)
            ib_c2            = ptr_patch%edges%cell_blk(il_e,ib_e,2)

            xx1              = gc2cc(ptr_patch%cells%center(il_c1,ib_c1))
            norm             = SQRT(SUM(xx1%x*xx1%x))
            xx1%x            = xx1%x/norm

            xx2              = gc2cc(ptr_patch%cells%center(il_c2,ib_c2))
            norm             = SQRT(SUM(xx2%x*xx2%x))
            xx2%x            = xx2%x/norm

            cell2cell_cc%x   = xx2%x - xx1%x
            IF(larc_length)THEN
              norm_c1_c2 = arc_length(xx1,xx2)
            ELSE
              norm_c1_c2 = SQRT(SUM(cell2cell_cc%x*cell2cell_cc%x))
            ENDIF
            dual_edge_length(ie) = norm_c1_c2
            !          cell2cell_cc%x       = cell2cell_cc%x/norm_c1_c2

            IF(mid_point_dual_edge)THEN
              cc_mid_dual_edge(ie)%x = 0.5_wp*(xx2%x+xx1%x)
              gc_mid_dual_edge(ie)   = cc2gc(cc_mid_dual_edge(ie))

              IF(coriolis_type==full_coriolis)THEN
                ptr_patch%edges%f_e(il_e, ib_e) = 2._wp*omega*SIN(gc_mid_dual_edge(ie)%lat)
              ELSEIF(coriolis_type==beta_plane_coriolis)THEN
                gc1%lat = basin_center_lat* deg2rad - 0.5_wp*basin_height_deg*deg2rad
                gc1%lon = 0.0_wp
                xx1     = gc2cc(gc1)

                gc2%lat = gc_mid_dual_edge(ie)%lat!*deg2rad
                gc2%lon = 0.0_wp
                xx2     = gc2cc(gc2)
                z_y     = re*arc_length(xx2,xx1)

                !z_y = ptr_patch%edges%center(je,jb)%lat - z_lat_basin_center
                ptr_patch%edges%f_e(il_e, ib_e) = &
                  & 2.0_wp*omega*( SIN(basin_center_lat * deg2rad) + &
                  & (COS(basin_center_lat * deg2rad)/re)*z_y)
              ENDIF
            ELSE
              cc_mid_dual_edge(ie)%x = cc_dual_edge(ie)%x
              gc_mid_dual_edge(ie)   = cc2gc(cc_mid_dual_edge(ie))
            ENDIF

            !2) determine vector from adjacent vertex to adjacent vertex
            !   in cartesian coordinate for moved dual cell
            !Get indices of two adjacent vertices
            il_v1 = ptr_patch%edges%vertex_idx(il_e,ib_e,1)
            ib_v1 = ptr_patch%edges%vertex_blk(il_e,ib_e,1)
            il_v2 = ptr_patch%edges%vertex_idx(il_e,ib_e,2)
            ib_v2 = ptr_patch%edges%vertex_blk(il_e,ib_e,2)

            xx1   = gc2cc(ptr_patch%verts%vertex(il_v1,ib_v1))
            norm  = SQRT(SUM(xx1%x*xx1%x))
            xx1%x = xx1%x/norm

            xx2   = gc2cc(ptr_patch%verts%vertex(il_v2,ib_v2))
            norm  = SQRT(SUM(xx2%x*xx2%x))
            xx2%x = xx2%x/norm

            vert1_midedge_cc(jv, jb, ie)%x = cc_mid_dual_edge(ie)%x - xx1%x
            vert2_midedge_cc(jv, jb, ie)%x = cc_mid_dual_edge(ie)%x - xx2%x
            norm = SQRT(SUM(vert1_midedge_cc(jv,jb,ie)%x*vert1_midedge_cc(jv,jb,ie)%x))
            vert1_midedge_cc(jv, jb, ie)%x = vert1_midedge_cc(jv, jb, ie)%x/norm

            norm = SQRT(SUM(vert2_midedge_cc(jv,jb,ie)%x*vert2_midedge_cc(jv,jb,ie)%x))
            vert2_midedge_cc(jv, jb, ie)%x = vert2_midedge_cc(jv, jb, ie)%x/norm


            !calculate vertex edge distance
            IF(larc_length)THEN
              vert_edge_dist(ie,1)     = arc_length (cc_dual_edge(ie), xx1)
              vert_edge_dist(ie,2)     = arc_length (cc_dual_edge(ie), xx2)
              vert_dual_mid_dist(ie,1) = arc_length (cc_mid_dual_edge(ie), xx1)
              vert_dual_mid_dist(ie,2) = arc_length (cc_mid_dual_edge(ie), xx2)
            ELSE
              vert_edge_dist(ie,1)&
                & = SQRT(SUM((cc_dual_edge(ie)%x - xx1%x)*(cc_dual_edge(ie)%x - xx1%x)))
              vert_edge_dist(ie,2)&
                & = SQRT(SUM((cc_dual_edge(ie)%x - xx2%x)*(cc_dual_edge(ie)%x - xx2%x)))
              vert_dual_mid_dist(ie,1)&
                & = SQRT(SUM((cc_mid_dual_edge(ie)%x - xx1%x)*(cc_mid_dual_edge(ie)%x - xx1%x)))
              vert_dual_mid_dist(ie,2)&
                & = SQRT(SUM((cc_mid_dual_edge(ie)%x - xx2%x)*(cc_mid_dual_edge(ie)%x - xx2%x)))
            ENDIF

            !calculate normal vector that is perpendicular to vertex-vertex- and edge position vector
            !If one uses the edge position vector this results in the moved primal normal.
            recon_vec_cc_v1(ie)   = vector_product(vert1_midedge_cc(jv, jb, ie),&
              & cc_mid_dual_edge(ie))
            norm                  = SQRT(SUM(recon_vec_cc_v1(ie)%x*recon_vec_cc_v1(ie)%x))
            recon_vec_cc_v1(ie)%x = recon_vec_cc_v1(ie)%x/norm

            recon_vec_cc_v2(ie)   = vector_product(vert2_midedge_cc(jv,jb,ie),&
              & cc_mid_dual_edge(ie))
            norm                  = SQRT(SUM(recon_vec_cc_v2(ie)%x*recon_vec_cc_v2(ie)%x))
            recon_vec_cc_v2(ie)%x = recon_vec_cc_v2(ie)%x/norm

            !Fix orientation
            z_tmp = DOT_PRODUCT(recon_vec_cc_v1(ie)%x,&
              & ptr_patch%edges%primal_cart_normal(il_e,ib_e)%x)
            IF (z_tmp <0._wp) recon_vec_cc_v1(ie)%x = -1._wp * recon_vec_cc_v1(ie)%x

            z_tmp = DOT_PRODUCT(recon_vec_cc_v2(ie)%x,&
              & ptr_patch%edges%primal_cart_normal(il_e,ib_e)%x)
            IF (z_tmp <0._wp) recon_vec_cc_v2(ie)%x = -1._wp * recon_vec_cc_v2(ie)%x


            IF ( (ptr_patch%edges%vertex_idx(il_e,ib_e,1) == jv) .AND. &
              & (ptr_patch%edges%vertex_blk(il_e,ib_e,1) == jb)         ) THEN

              vert_edge_distance     = vert_edge_dist(ie,1)
              vert_dual_mid_distance = vert_dual_mid_dist(ie,1)
              recon_vec_cc           = recon_vec_cc_v1(ie)
              !PK: not used
              !              ptr_intp%edge2vert_vector_cc(jv,jb,ie)=&
              !              &vector_product(vert1_midedge_cc(jv, jb, ie), cc_mid_dual_edge(ie))
              ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x=vert1_midedge_cc(jv, jb, ie)%x

              z_tmp = DOT_PRODUCT(ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x,&
                & ptr_patch%edges%primal_cart_normal(il_e,ib_e)%x)

              IF (z_tmp <0._wp) ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x&
                & = -1._wp * ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x


              ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x = vert_edge_distance&
                & *ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x

            ELSE IF ( (ptr_patch%edges%vertex_idx(il_e,ib_e,2) == jv) .AND. &
              & (ptr_patch%edges%vertex_blk(il_e,ib_e,2) == jb) ) THEN

              vert_edge_distance     = vert_edge_dist(ie,2)
              vert_dual_mid_distance = vert_dual_mid_dist(ie,2)
              recon_vec_cc           = recon_vec_cc_v2(ie)

              !              ptr_intp%edge2vert_vector_cc(jv,jb,ie)=&
              !              &vector_product(vert2_midedge_cc(jv, jb, ie), cc_mid_dual_edge(ie))

              ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x=vert2_midedge_cc(jv,jb,ie)%x

              z_tmp = DOT_PRODUCT(ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x,&
                & ptr_patch%edges%primal_cart_normal(il_e,ib_e)%x)

              IF (z_tmp <0._wp) ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x =&
                & -1._wp * ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x

              ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x = vert_edge_distance&
                & *ptr_intp%edge2vert_vector_cc(jv,jk,jb,ie)%x


            ELSE
              CALL message (TRIM(routine), 'WARNING - vert_edge_distance not found')
              WRITE(0,'(a,7i5)') &
                & 'jv, jb, edge, ptr_patch%edges%vertex_idx/blk(il_e,ib_e,1-2)=', &
                & jv, jb, ie, &
                & ptr_patch%edges%vertex_idx(il_e,ib_e,1), &
                & ptr_patch%edges%vertex_blk(il_e,ib_e,1), &
                & ptr_patch%edges%vertex_idx(il_e,ib_e,2), &
                & ptr_patch%edges%vertex_blk(il_e,ib_e,2)

            END IF

            ptr_intp%variable_dual_vol_norm(jv,jk,jb,ie) = &
              & 0.5_wp*dual_edge_length(ie)*vert_dual_mid_distance
            !vert_edge_distance*dual_edge_length(ie)!

            ptr_intp%edge2vert_coeff_cc(jv,jk,jb,ie)%x   = &
              & recon_vec_cc%x*dual_edge_length(ie)*vert_dual_mid_distance

            norm_v1_v2 = &
              & SQRT(SUM(vert1_midedge_cc(jv, jb, ie)%x*vert1_midedge_cc(jv, jb, ie)%x)) + &
              & SQRT(SUM(vert2_midedge_cc(jv, jb, ie)%x*vert2_midedge_cc(jv, jb, ie)%x))

            ptr_intp%edge2vert_coeff_cc_t(il_e,jk,ib_e,1)%x = vert1_midedge_cc(jv, jb, ie)%x * &
              & ( ptr_patch%edges%system_orientation(il_e,ib_e)/norm_v1_v2 )

            ptr_intp%edge2vert_coeff_cc_t(il_e,jk,ib_e,2)%x = vert2_midedge_cc(jv, jb, ie)%x * &
              & ( ptr_patch%edges%system_orientation(il_e,ib_e)/norm_v1_v2 )

          END DO
        ENDDO vert_idx_loop
      END DO vert_blk_loop
    END DO

    !     !--------------------------------------------------------------------------
    !     ! SYNCHRONIZE ALL ELEMENTS OF V_BASE:
    !     ! synchronize elements on cells
    !     DO ie = 1, no_cell_edges
    !       DO icc = 1, 3
    !         z_sync_c(:,:,:) =  ptr_intp%edge2cell_coeff_cc(:,:,:,ie)%x(icc)
    !         CALL sync_patch_array(SYNC_C, ptr_patch, z_sync_c(:,:,:))
    !         ptr_intp%edge2cell_coeff_cc(:,:,:,ie)%x(icc) = z_sync_c(:,:,:)
    !       END DO
    !       z_sync_c(:,:,:) = ptr_intp%variable_vol_norm(:,:,:,ie)
    !       CALL sync_patch_array(SYNC_C, ptr_patch, z_sync_c(:,:,:))
    !       ptr_intp%variable_vol_norm(:,:,:,ie) = z_sync_c(:,:,:)
    !     END DO
    !     CALL sync_patch_array(SYNC_C, ptr_patch,ptr_intp%fixed_vol_norm)
    !
    !     ! synchronize elements on edges
    !     DO ie = 1, 2
    !       DO icc = 1, 3
    !         z_sync_e(:,:,:) =  ptr_intp%edge2vert_coeff_cc_t(:,:,:,ie)%x(icc)
    !         CALL sync_patch_array(SYNC_E, ptr_patch, z_sync_e(:,:,:))
    !         ptr_intp%edge2vert_coeff_cc_t(:,:,:,ie)%x(icc) = z_sync_e(:,:,:)
    !
    !         z_sync_e(:,:,:) =  ptr_intp%edge2cell_coeff_cc_t(:,:,:,ie)%x(icc)
    !         CALL sync_patch_array(SYNC_E, ptr_patch, z_sync_e(:,:,:))
    !         ptr_intp%edge2cell_coeff_cc_t(:,:,:,ie)%x(icc) = z_sync_e(:,:,:)
    !       END DO
    !     END DO
    !
    !     ! synchronize cartesian coordinates on vertices:
    !     DO ie = 1, no_vert_edges
    !       DO icc = 1, 3
    !         z_sync_v(:,:,:) =  ptr_intp%edge2vert_vector_cc(:,:,:,ie)%x(icc)
    !         CALL sync_patch_array(SYNC_V, ptr_patch, z_sync_v(:,:,:))
    !         ptr_intp%edge2vert_vector_cc(:,:,:,ie)%x(icc) = z_sync_v(:,:,:)
    !
    !         z_sync_v(:,:,:) = ptr_intp%edge2vert_coeff_cc(:,:,:,ie)%x(icc)
    !         CALL sync_patch_array(SYNC_V, ptr_patch, z_sync_v(:,:,:))
    !         ptr_intp%edge2vert_coeff_cc(:,:,:,ie)%x(icc) = z_sync_v(:,:,:)
    !       END DO
    !       z_sync_v(:,:,:) = ptr_intp%variable_dual_vol_norm(:,:,:,ie)
    !       CALL sync_patch_array(SYNC_V, ptr_patch, z_sync_v(:,:,:))
    !       ptr_intp%variable_dual_vol_norm(:,:,:,ie) = z_sync_v(:,:,:)
    !     END DO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_scalar_product_oce_3d
  !-------------------------------------------------------------------------


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
  SUBROUTINE init_geo_factors_oce_3d( ptr_patch, ptr_intp )
    !
    IMPLICIT NONE
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch
    TYPE(t_operator_coeff),     INTENT(inout) :: ptr_intp
    !

    INTEGER :: jc, jb, je, jv, ie,jk
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER :: ile, ibe!, ilc1, ibc1, ilc2, ibc2, ifac, ic, ilnc, ibnc
    !INTEGER :: ile1, ibe1,ile2,ibe2,ile3,ibe3
    INTEGER, PARAMETER :: i_cell_type = 3
    !TYPE(cartesian_coordinates)::z_pn_k,z_pn_j
    !REAL(wp) :: z_lon, z_lat, z_nu, z_nv, z_proj
    !REAL(wp) :: cell_area

    REAL(wp) :: z_sync_c(nproma,n_zlev,ptr_patch%nblks_c)
    !REAL(wp) :: z_sync_e(nproma,n_zlev,ptr_patch%nblks_e)
    REAL(wp) :: z_sync_v(nproma,n_zlev,ptr_patch%nblks_v)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:init_geo_factors_oce_3d')
    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')
    i_nchdom   = MAX(1,ptr_patch%n_childdom)

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

    ! 1) coefficients for divergence
    rl_start = 1
    rl_end = min_rlcell

    ! values for the blocking
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jc,i_startidx,i_endidx,ile,ibe)
    DO jk=1,n_zlev
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,      &
          & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          !         ile1 = ptr_patch%cells%edge_idx(jc,jb,1)
          !         ibe1 = ptr_patch%cells%edge_blk(jc,jb,1)
          !         ile2 = ptr_patch%cells%edge_idx(jc,jb,2)
          !         ibe2 = ptr_patch%cells%edge_blk(jc,jb,2)
          !         ile3 = ptr_patch%cells%edge_idx(jc,jb,3)
          !         ibe3 = ptr_patch%cells%edge_blk(jc,jb,3)
          !         cell_area =  0.25_wp&
          !         & *( ptr_patch%edges%primal_edge_length(ile1,ibe1)*ptr_patch%edges%dual_edge_length(ile1,ibe1)&
          !         &   +ptr_patch%edges%primal_edge_length(ile2,ibe2)*ptr_patch%edges%dual_edge_length(ile2,ibe2)&
          !         &   +ptr_patch%edges%primal_edge_length(ile3,ibe3)*ptr_patch%edges%dual_edge_length(ile3,ibe3))

          DO je = 1, i_cell_type

            ile = ptr_patch%cells%edge_idx(jc,jb,je)
            ibe = ptr_patch%cells%edge_blk(jc,jb,je)

            ptr_intp%div_coeff(jc,jk,jb,je) =                &
              & ptr_patch%edges%primal_edge_length(ile,ibe) * &
              & ptr_patch%cells%edge_orientation(jc,jb,je)  / &
              & ptr_patch%cells%area(jc,jb)
          ENDDO !edge loop
          !write(1234,*)'div coeff 3D',jk,jc,jb,ptr_intp%div_coeff(jc,jk,jb,:)
        ENDDO !idx loop
      END DO !block loop
    END DO
!$OMP END DO

    ! 2) coefficients for curl
    rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    rl_end = min_rlvert_int

    ! values for the blocking
    i_startblk = ptr_patch%verts%start_blk(rl_start,1)
    i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jv,i_startidx,i_endidx,ile,ibe)
    DO jk=1,n_zlev
      DO jb = i_startblk, i_endblk

        CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO jv = i_startidx, i_endidx
          DO je = 1, ptr_patch%verts%num_edges(jv,jb)

            ile = ptr_patch%verts%edge_idx(jv,jb,je)
            ibe = ptr_patch%verts%edge_blk(jv,jb,je)

            ptr_intp%rot_coeff(jv,jk,jb,je) =                &
              & ptr_patch%edges%dual_edge_length(ile,ibe) * &
              & ptr_patch%verts%edge_orientation(jv,jb,je)!/ &
            !&    ptr_patch%verts%dual_area(jv,jb)
          ENDDO
        ENDDO
      END DO
    END DO
!$OMP END DO

    ! 3) coefficients for nabla2_scalar
    rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    rl_end = min_rlcell_int

    ! values for the blocking
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
    !!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
    !!$OMP    ilc2,ibc2,ilnc,ibnc)
    ! !       DO jk=1,n_zlev
    ! !       DO jb = i_startblk, i_endblk
    ! !
    ! !         CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
    ! !                            i_startidx, i_endidx, rl_start, rl_end)
    ! !
    ! !         DO je = 1, i_cell_type
    ! !           DO jc = i_startidx, i_endidx
    ! !
    ! !             ile = ptr_patch%cells%edge_idx(jc,jb,je)
    ! !             ibe = ptr_patch%cells%edge_blk(jc,jb,je)
    ! !
    ! !             ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
    ! !             ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
    ! !             ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
    ! !             ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)
    ! !
    ! !             IF (jc == ilc1 .AND. jb == ibc1) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 ptr_intp%geofac_n2s(jc,1,jb)     =  &
    ! !                 &  ptr_intp%geofac_n2s(jc,1,jb)  -  &
    ! !                 &  ptr_intp%geofac_div(jc,je,jb) /  &
    ! !                 &  ptr_patch%edges%dual_edge_length(ile,ibe)
    ! !             ENDIF
    ! !           ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
    ! !             IF (i_cell_type == 3) THEN
    ! !               ptr_intp%geofac_n2s(jc,1,jb)       =  &
    ! !                 &  ptr_intp%geofac_n2s(jc,1,jb)  +  &
    ! !                 &  ptr_intp%geofac_div(jc,je,jb) /  &
    ! !                 &  ptr_patch%edges%dual_edge_length(ile,ibe)
    ! !             ENDIF
    ! !           ENDIF
    ! !           DO ic = 1, i_cell_type
    ! !             ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
    ! !             ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
    ! !             IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 ptr_intp%geofac_n2s(jc,ic+1,jb)     = &
    ! !                   &  ptr_intp%geofac_n2s(jc,ic+1,jb)- &
    ! !                   &  ptr_intp%geofac_div(jc,je,jb)  / &
    ! !                   &  ptr_patch%edges%dual_edge_length(ile,ibe)
    ! !               ENDIF
    ! !             ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 ptr_intp%geofac_n2s(jc,ic+1,jb)     = &
    ! !                   &  ptr_intp%geofac_n2s(jc,ic+1,jb)+ &
    ! !                   &  ptr_intp%geofac_div(jc,je,jb)  / &
    ! !                   &  ptr_patch%edges%dual_edge_length(ile,ibe)
    ! !               ENDIF
    ! !             ENDIF
    ! !           ENDDO
    ! !
    ! !           ! To ensure that dummy edges have a factor of 0:
    ! !           IF (je > ptr_patch%cells%num_edges(jc,jb)) THEN
    ! !             ptr_intp%geofac_n2s(jc,je+1,jb) = 0._wp
    ! !           ENDIF
    ! !
    ! !         ENDDO !cell loop
    ! !       ENDDO
    ! !     END DO !block loop
    ! !     END DO
    !!$OMP END DO


    ! 4) coefficients for gradient

    rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    rl_end = min_rledge_int

    ! Second step: computed projected orientation vectors and related information
    i_startblk = ptr_patch%edges%start_blk(rl_start,1)
    i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch edges
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je)
    DO jk=1,n_zlev
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO je =  i_startidx, i_endidx
          ! compute inverse dual edge length (undefined for refin_ctrl=1)
          ptr_patch%edges%inv_dual_edge_length(je,jb) = &
            & 1._wp/ptr_patch%edges%dual_edge_length(je,jb)

          ptr_intp%grad_coeff(je,jk,jb)&
            & = ptr_patch%edges%inv_dual_edge_length(je,jb)

        ENDDO
      END DO !block loop
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! synchronize all elements of ptr_intp:
    DO ie = 1, i_cell_type

      z_sync_c(:,:,:) = ptr_intp%div_coeff(:,:,:,ie)
      CALL sync_patch_array(sync_c, ptr_patch, z_sync_c(:,:,:))
      ptr_intp%div_coeff(:,:,:,ie) = z_sync_c(:,:,:)
      !         z_sync_c(:,:,:) = ptr_intp%n2s_coeff(:,:,:,ie)
      !         CALL sync_patch_array(SYNC_C, ptr_patch, z_sync_c(:,:,:))
      !         ptr_intp%n2s_coeff(:,:,:,ie) = z_sync_c(:,:,:)
    END DO

    CALL sync_patch_array(sync_e, ptr_patch, ptr_patch%edges%inv_dual_edge_length(:,:))

    DO ie = 1, 9-i_cell_type
      z_sync_v(:,:,:) = ptr_intp%rot_coeff(:,:,:,ie)
      CALL sync_patch_array(sync_v, ptr_patch, z_sync_v(:,:,:))
      ptr_intp%rot_coeff(:,:,:,ie) = z_sync_v(:,:,:)
    END DO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_geo_factors_oce_3d
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

END MODULE mo_operator_ocean_coeff_3d

