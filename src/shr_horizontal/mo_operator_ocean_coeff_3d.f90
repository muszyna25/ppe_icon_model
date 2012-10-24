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
    &                               max_char_length, beta_plane_coriolis,full_coriolis, &
    &                               min_rledge_int,min_rlcell_int,min_rlvert_int,&
    &                               sea_boundary, boundary, sea, min_dolic
  USE mo_math_constants,      ONLY: deg2rad, pi
  USE mo_physical_constants,  ONLY: earth_radius
  USE mo_math_utilities,      ONLY: gc2cc, cc2gc, t_cartesian_coordinates,      &
    &                               t_geographical_coordinates, vector_product, &
    &                               arc_length
  USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m, no_tracer, t_ref, s_ref, &
    &                               coriolis_type, basin_center_lat, basin_height_deg
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D_oce
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, sync_idx, global_max
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_oce_state,           ONLY: t_hydro_ocean_state
  USE mo_oce_physics,         ONLY: t_ho_params
  !USE mo_intp_data_strc,      ONLY: t_int_state
  !USE mo_intp_coeffs,         ONLY: par_init_scalar_product_oce
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_grid_config,         ONLY: grid_sphere_radius, grid_angular_velocity
  USE mo_run_config,          ONLY: dtime

  IMPLICIT NONE


  PRIVATE

  PUBLIC  :: t_operator_coeff
  PUBLIC  :: allocate_exp_coeff
  PUBLIC  :: par_init_operator_coeff2

  PRIVATE :: par_init_coeff_2D
  PRIVATE :: copy_2D_to_3D_coeff
  PRIVATE :: par_apply_boundary2coeffs
  PRIVATE :: init_geo_factors_oce_3d
  PUBLIC  :: update_diffusion_matrices
  !PRIVATE  :: par_init_operator_coeff
  !PRIVATE :: copy_2D_to_3D_intp_coeff
  !PRIVATE :: init_operator_coeff 
  !PRIVATE  :: apply_boundary2coeffs
  !PRIVATE :: init_scalar_product_oce_3d


  INTEGER,PARAMETER :: no_dual_edges   = 6
  !INTEGER,PARAMETER :: no_primal_edges = 3

  ! flags for computing ocean coefficients
  LOGICAL, PARAMETER :: MID_POINT_DUAL_EDGE = .TRUE. !Please do not change this unless you are sure, you know what you do.
  LOGICAL, PARAMETER :: LARC_LENGTH = .FALSE.



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

    REAL(wp), ALLOCATABLE :: matrix_vert_diff_c(:,:,:,:)
    REAL(wp), ALLOCATABLE :: matrix_vert_diff_e(:,:,:,:)


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
  SUBROUTINE allocate_exp_coeff( p_patch, p_coeff)
    ! !
    TYPE(t_patch),TARGET,INTENT(in)       :: p_patch
    TYPE(t_operator_coeff), INTENT(inout) :: p_coeff

    INTEGER :: nblks_c, nblks_e, nblks_v, nz_lev
    INTEGER :: ist,ie
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: jc,je,jk,jb

    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: all_verts

    !INTEGER :: edge_block, cell_block, vertex_block, level, neigbor
    INTEGER :: i_startidx_e, i_endidx_e

    !-----------------------------------------------------------------------
    !
    ! determine size of arrays, i.e.
    ! values for the blocking
    !
    nblks_c  = p_patch%nblks_c
    nblks_e  = p_patch%nblks_e
    nblks_v  = p_patch%nblks_v
    nz_lev   = n_zlev

    ALLOCATE(p_coeff%div_coeff(nproma,n_zlev,nblks_c,p_patch%cell_type),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for geofac_div failed')
    ENDIF

    ALLOCATE(p_coeff%grad_coeff(nproma,n_zlev,nblks_e),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for geofac_grad failed')
    ENDIF

    ALLOCATE(p_coeff%rot_coeff(nproma,n_zlev,nblks_v,no_dual_edges),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d', &
        & 'allocation for geofac_rot failed')
    ENDIF

    ALLOCATE(p_coeff%n2s_coeff(nproma,n_zlev,nblks_c,p_patch%cell_type+1),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                       &
        & 'allocation for geofac_n2s failed')
    ENDIF
    ALLOCATE(p_coeff%n2v_coeff(nproma,n_zlev,nblks_e),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                       &
        & 'allocation for geofac_n2v failed')
    ENDIF
    !
    ALLOCATE(p_coeff%dist_cell2edge(nproma,n_zlev,nblks_e,2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating dist_cell2edge failed')
    ENDIF

    ALLOCATE(p_coeff%bnd_edge_idx(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating bnd_edge_idx failed')
    ENDIF
    ALLOCATE(p_coeff%bnd_edge_blk(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating bnd_edge_blk failed')
    ENDIF
    ALLOCATE(p_coeff%edge_idx(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge_idx failed')
    ENDIF
    ALLOCATE(p_coeff%orientation(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating orientation failed')
    ENDIF
    ALLOCATE(p_coeff%bnd_edges_per_vertex(nproma,n_zlev,nblks_v),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating bnd_edges_per_vertex failed')
    ENDIF
    ALLOCATE(p_coeff%upwind_cell_idx(nproma,n_zlev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind_cell_idx failed')
    ENDIF
    ALLOCATE(p_coeff%upwind_cell_blk(nproma,n_zlev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind_cell_blk failed')
    ENDIF
    !
    ! arrays that are required for setting up the scalar product
    !
    !coefficients for edge to cell mapping, one half of the scalar product.
    !Dimension: nproma,nblks_c encode number of cells, 1:3 corresponds to number
    !of edges per cell, 1:2 is for u and v component of cell vector
    !     ALLOCATE(p_coeff%edge2cell_coeff(nproma,nblks_c,1:3, 1:2),STAT=ist)
    !     IF (ist /= SUCCESS) THEN
    !       CALL finish ('allocating edge2cell_coeff failed')
    !     ENDIF
    ALLOCATE(p_coeff%edge2cell_coeff_cc(nproma,nz_lev,nblks_c,1:3),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc failed')
    ENDIF

    ALLOCATE(p_coeff%edge2cell_coeff_cc_dyn(nproma,1,nblks_c,1:3),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc_dyn failed')
    ENDIF
    ALLOCATE(p_coeff%edge2vert_coeff_cc_dyn(nproma,1,nblks_c,1:3),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff_cc_dyn failed')
    ENDIF

    !coefficients for transposed of edge to cell mapping, second half of the scalar product.
    !Dimension: nproma,nblks_e encode number of edges, 1:2 is for cell neighbors of an edge
    ALLOCATE(p_coeff%edge2cell_coeff_cc_t(nproma,nz_lev,nblks_e,1:2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating transposed edge2cell_coeff failed')
    ENDIF

    !
    !coefficients for edge to vertex mapping.
    !
    !Dimension: nproma,nblks_v encode number of vertices,
    !1:6 is number of edges of a vertex,
    !1:2 is for u and v component of vertex vector
    ALLOCATE(p_coeff%edge2vert_coeff_cc(nproma,nz_lev,nblks_v,1:6),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff failed')
    ENDIF

    ALLOCATE(p_coeff%edge2vert_coeff_cc_t(nproma,nz_lev,nblks_e,1:2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff failed')
    ENDIF
    ALLOCATE(p_coeff%edge2vert_vector_cc(nproma,nz_lev,nblks_v,1:6),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_vector failed')
    ENDIF

    ALLOCATE(p_coeff%upwind_cell_position_cc(nproma,nz_lev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind cell failed')
    ENDIF
    ALLOCATE(p_coeff%moved_edge_position_cc(nproma,nz_lev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge failed')
    ENDIF
    ALLOCATE(p_coeff%edge_position_cc(nproma,nz_lev,nblks_e),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge failed')
    ENDIF
    ALLOCATE(p_coeff%cell_position_cc(nproma,nz_lev,nblks_c),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating cell failed')
    ENDIF
    !
    !normalizing factors for edge to cell mapping.
    !
    !Either by fixed volume or by variable one taking the surface elevation
    !into account. The later one depends on time and space.
    ALLOCATE(p_coeff%fixed_vol_norm(nproma,nz_lev,nblks_c),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating fixed_vol_norm failed')
    ENDIF
    ALLOCATE(p_coeff%variable_vol_norm(nproma,nz_lev,nblks_c,1:3),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating variable_vol_norm failed')
    ENDIF

    ALLOCATE(p_coeff%variable_dual_vol_norm(nproma,nz_lev,nblks_v,1:6),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating variable_dual_vol_norm failed')
    ENDIF

    !The last index "3" comes from the fact that we use a tridiagonal matrix
    ALLOCATE(p_coeff%matrix_vert_diff_c(nproma,n_zlev,nblks_c, 3),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for matrix_vert_diff_c failed')
    ENDIF
    !The last index "3" comes from the fact that we use a tridiagonal matrix
    ALLOCATE(p_coeff%matrix_vert_diff_e(nproma,n_zlev,nblks_e, 3),&
      & stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for matrix_vert_diff_e failed')
    ENDIF

    !
    ! initialize all components
    !
    DO ie = 1,3
      p_coeff%edge2cell_coeff_cc%x(ie)     = 0._wp
      p_coeff%edge2cell_coeff_cc_t%x(ie)   = 0._wp
      p_coeff%edge2vert_coeff_cc%x(ie)     = 0._wp
      p_coeff%edge2vert_coeff_cc_t%x(ie)   = 0._wp
      p_coeff%edge2vert_vector_cc%x(ie)    = 0._wp
      p_coeff%edge2cell_coeff_cc_dyn%x(ie) = 0._wp
      p_coeff%edge2vert_coeff_cc_dyn%x(ie) = 0._wp
    END DO

    all_cells => p_patch%cells%all
    all_edges => p_patch%edges%all
    all_verts => p_patch%verts%all

    DO jk = 1, nz_lev
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO je =  i_startidx_e, i_endidx_e
          p_coeff%edge_position_cc(je,jk,jb)             = gc2cc(p_patch%edges%center(je,jb))
          p_coeff%moved_edge_position_cc(je,jk,jb)%x(:)  = 0._wp
          p_coeff%upwind_cell_position_cc(je,jk,jb)%x(:) = 0._wp
        END DO
      END DO
    END DO

    DO jk = 1, nz_lev
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          p_coeff%cell_position_cc(jc,jk,jb)&
            & = gc2cc(p_patch%cells%center(jc,jb))

        END DO
      END DO
    END DO

    p_coeff%fixed_vol_norm         = 0._wp
    p_coeff%variable_vol_norm      = 0._wp
    p_coeff%variable_dual_vol_norm = 0._wp

    p_coeff%dist_cell2edge = 0._wp

    p_coeff%div_coeff  = 0._wp
    p_coeff%rot_coeff  = 0._wp
    p_coeff%grad_coeff = 0._wp
    !p_coeff%n2s_coeff  = 0._wp
    !p_coeff%n2v_coeff  = 0._wp

    p_coeff%bnd_edge_idx = 0
    p_coeff%bnd_edge_blk = 0
    p_coeff%edge_idx     = 0
    p_coeff%orientation  = 0.0_wp
    p_coeff%bnd_edges_per_vertex= 0

    p_coeff%upwind_cell_idx = 1
    p_coeff%upwind_cell_blk = 1

    p_coeff%matrix_vert_diff_c= 0.0_wp
    p_coeff%matrix_vert_diff_e= 0.0_wp

    CALL message ('mo_operator_ocean_coeff_3d:allocate_exp_coeff',&
      & 'memory allocation finished')

  END SUBROUTINE allocate_exp_coeff
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
!  

  !-------------------------------------------------------------------------
  !> Initialize expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !!
  SUBROUTINE par_init_operator_coeff2( patch, p_patch_3D, p_phys_param, ocean_coeff)
    !
    TYPE(t_patch),            INTENT(inout) :: patch
    TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
    TYPE (t_ho_params),       INTENT(IN)    :: p_phys_param
    TYPE(t_operator_coeff),   INTENT(inout) :: ocean_coeff

   !
   !Local variables: strcutures for 2D coefficients
   !
    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc(1:nproma,1:patch%nblks_c,1:patch%cell_type)
    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc_t(1:nproma,1:patch%nblks_e,1:2)

    TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc(1:nproma,1:patch%nblks_v,1:6)
    TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc_t(1:nproma,1:patch%nblks_e,1:2)

    REAL(wp) :: dist_cell2edge(1:nproma,1:patch%nblks_e,1:2)
    REAL(wp) :: fixed_vol_norm(1:nproma,patch%nblks_c)
    REAL(wp) :: variable_vol_norm(1:nproma,1:patch%nblks_c,1:3)
    REAL(wp) :: variable_dual_vol_norm(1:nproma,1:patch%nblks_e,1:6)
    !-----------------------------------------------------------------------
!     TYPE(t_cartesian_coordinates) :: check_v(nproma, n_zlev, patch%nblks_v, 6)
!     REAL(wp) :: check_r(nproma, n_zlev, patch%nblks_c, 3)
!     REAL(wp) :: max_diff, max_val
    !-----------------------------------------------------------------------
    !ocean_coeff%dist_cell2edge(:,:,:,:) = 0.0_wp
!write(*,*)'before par_init_coeff_2D'
    CALL par_init_coeff_2D( patch,                &
                          & edge2cell_coeff_cc,   &
                          & edge2cell_coeff_cc_t, &
                          & edge2vert_coeff_cc,   &
                          & edge2vert_coeff_cc_t, &
                          & dist_cell2edge,       &
                          & fixed_vol_norm,       &
                          & variable_vol_norm,    &
                          & variable_dual_vol_norm)
!write(*,*)'beforecopy_2D_to_3D_coeff'
    CALL copy_2D_to_3D_coeff( patch,                &
                            & ocean_coeff,          &
                            & edge2cell_coeff_cc,   &
                            & edge2cell_coeff_cc_t, &
                            & edge2vert_coeff_cc,   &
                            & edge2vert_coeff_cc_t, &
                            & dist_cell2edge,       &
                            & fixed_vol_norm,       &
                            & variable_vol_norm)!,    &
                            !& variable_dual_vol_norm)
!write(*,*)'before init_geo_factors_oce_3d'
    CALL init_geo_factors_oce_3d ( patch, ocean_coeff )
    !---------------------------------------------------------
    CALL par_apply_boundary2coeffs(patch, p_patch_3D, ocean_coeff)


     CALL update_diffusion_matrices( patch, p_patch_3D,               &
                                     & p_phys_param,                  &
                                     & ocean_coeff%matrix_vert_diff_e,&
                                     & ocean_coeff%matrix_vert_diff_c)
    RETURN
    !---------------------------------------------------------
    ! checks
!     check_v = ocean_coeff%edge2vert_coeff_cc
!     check_r = ocean_coeff%variable_vol_norm
!     !---------------------------------------------------------
!      CALL init_operator_coeff( patch, ocean_coeff )
!
!      max_diff = MAXVAL(ABS(ocean_coeff%edge2vert_coeff_cc(:,:,:,:)%x(1) - &
!        &  check_v(:,:,:,:)%x(1) ))
!      max_val  =  MAXVAL(ABS( check_v(:,:,:,:)%x(1)))
!      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(1)=", max_diff, max_val
!      max_diff = MAXVAL(ABS(ocean_coeff%edge2vert_coeff_cc(:,:,:,:)%x(2) - &
!        & check_v(:,:,:,:)%x(2) ))
!      max_val  =  MAXVAL(ABS( check_v(:,:,:,:)%x(2)))
!      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(2)=", max_diff, max_val
!      max_diff = MAXVAL(ABS(ocean_coeff%edge2vert_coeff_cc(:,:,:,:)%x(3) - &
!        & check_v(:,:,:,:)%x(3) ))
!      max_val  =  MAXVAL(ABS( check_v(:,:,:,:)%x(3)))
!      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(3)=", max_diff, max_val
!      max_diff = MAXVAL(ABS(ocean_coeff%variable_vol_norm - check_r ))
!      max_val  =  MAXVAL(ABS( check_r ))
!      Write(0,*) "max diff of variable_vol_norm=", max_diff, max_val
!
!      STOP

  END SUBROUTINE par_init_operator_coeff2
  !-------------------------------------------------------------------------
  !> Initialize 3D expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !! Parellelized by Leonidas Linardakis 2012-3
 
 
   SUBROUTINE update_diffusion_matrices( patch, p_patch_3D,    &
                                         & p_phys_param,       &
                                         & matrix_vert_diff_e, &
                                         & matrix_vert_diff_c)
 
     TYPE(t_patch), TARGET, INTENT(INOUT)  :: patch
     TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
     TYPE (t_ho_params),       INTENT(IN) :: p_phys_param
     REAL(wp), INTENT(INOUT) :: matrix_vert_diff_e(1:nproma,1:n_zlev,1:patch%nblks_e,1:3)
     REAL(wp), INTENT(INOUT) :: matrix_vert_diff_c(1:nproma,1:n_zlev,1:patch%nblks_c,1:3)
    !
    !Local variables
    !
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp), POINTER :: A_v(:,:,:)
    REAL(wp) :: inv_zinv_i(1:n_zlev)
    REAL(wp) :: inv_zinv_m(1:n_zlev)
    REAL(wp) :: dt_inv
    INTEGER  :: je,jc,jb,jk, i_no_t
    INTEGER  :: slev,z_dolic
    INTEGER  :: i_startidx_e, i_endidx_e,  i_startidx_c, i_endidx_c
    !---------------------------------------------------------
    all_cells => patch%cells%all
    all_edges => patch%edges%all
    !---------------------------------------------------------
    slev   = 1
    dt_inv = 1.0_wp/dtime

    !The vertical diffusion matrices for the tracers
    DO i_no_t = 1,no_tracer

      A_v => p_phys_param%A_tracer_v(:,:,:, i_no_t)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)

          !IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN 
          IF ( p_patch_3D%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN 
            IF ( z_dolic >=MIN_DOLIC ) THEN

              !inv_zinv_i(:) = 1.0_wp/v_base%del_zlev_i(:)
              !inv_zinv_m(:) = 1.0_wp/v_base%del_zlev_m(:)
               !inv_zinv_i(:) = p_os%p_diag%inv_prism_center_dist_c(jc,:,jb)
               !inv_zinv_m(:) = p_os%p_diag%inv_prism_thick_c(jc,:,jb)
              inv_zinv_i(:) = p_patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,:,jb)
              inv_zinv_m(:) = p_patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,:,jb)


              !Fill triangular matrix
              !b is diagonal a and c are upper and lower band
              !This corresponds to the 4th indices: "2", "1" and "3"
              DO jk = slev+1, z_dolic-1
                !a(jk) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)!*dtime
                !c(jk) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)!*dtime
                !b(jk) = dt_inv-a(jk)-c(jk)
                matrix_vert_diff_c(jc,jk,jb,1) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)
                matrix_vert_diff_c(jc,jk,jb,3) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)
                matrix_vert_diff_c(jc,jk,jb,2)  = dt_inv&
                                               & -matrix_vert_diff_c(jc,jk,jb,1)&
                                               & -matrix_vert_diff_c(jc,jk,jb,3)
              END DO
!write(*,*)'coeff 521',matrix_vert_diff_c(5,2,1,2),matrix_vert_diff_c(5,2,1,1),matrix_vert_diff_c(5,2,1,3),&
!&-A_v(5,2,1),inv_zinv_m(2) ,inv_zinv_i(2),A_v(5,3,1),inv_zinv_m(2),inv_zinv_i(3)
              ! The first row
              !c(slev) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)
              !a(slev) = 0.0_wp           
              !b(slev) = dt_inv- c(slev)!! - a(slev) 
               
              matrix_vert_diff_c(jc,slev,jb,3) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)
              matrix_vert_diff_c(jc,slev,jb,1) = 0.0_wp           
              matrix_vert_diff_c(jc,slev,jb,2) = dt_inv- matrix_vert_diff_c(jc,slev,jb,3) 
!write(*,*)'coeff 511',matrix_vert_diff_c(5,1,1,2),matrix_vert_diff_c(5,1,1,1),matrix_vert_diff_c(5,1,1,3),&
!&-A_v(jc,slev+1,jb),inv_zinv_m(slev),inv_zinv_i(slev+1)
              ! The last row
              !a(z_dolic) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)
              !c(z_dolic) = 0.0_wp
              !b(z_dolic) = dt_inv - a(z_dolic)!! - c(z_dolic)
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
      z_dolic = p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)!v_base%dolic_e(je,jb)

      !IF ( v_base%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN
      IF (  p_patch_3D%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN
        IF ( z_dolic >= MIN_DOLIC ) THEN

          !inv_zinv_i(:)=1.0_wp/v_base%del_zlev_i(:)
          !inv_zinv_m(:)=1.0_wp/v_base%del_zlev_m(:)
          !inv_zinv_i(:) = p_os%p_diag%inv_prism_center_dist_e(je,:,jb)
          !inv_zinv_m(:) = p_os%p_diag%inv_prism_thick_e(je,:,jb)
          inv_zinv_i(:) = p_patch_3D%p_patch_1D(1)%inv_prism_center_dist_e(je,:,jb)
          inv_zinv_m(:) = p_patch_3D%p_patch_1D(1)%inv_prism_thick_e(je,:,jb)


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
  !!  Parallelized by Leonidas Linardakis, 2012-3
  SUBROUTINE par_init_coeff_2D( patch,               &
                                        & edge2cell_coeff_cc,   &
                                        & edge2cell_coeff_cc_t, &
                                        & edge2vert_coeff_cc,   &
                                        & edge2vert_coeff_cc_t, &
                                        & dist_cell2edge,       &
                                        & fixed_vol_norm,       &
                                        & variable_vol_norm,    &
                                        & variable_dual_vol_norm)
    TYPE(t_patch)    , TARGET, INTENT(INOUT)     :: patch
    TYPE(t_cartesian_coordinates), INTENT(INOUT) :: edge2cell_coeff_cc(1:nproma,1:patch%nblks_c,1:patch%cell_type)
    TYPE(t_cartesian_coordinates), INTENT(INOUT) :: edge2cell_coeff_cc_t(1:nproma,1:patch%nblks_e,1:2)
    TYPE(t_cartesian_coordinates), INTENT(INOUT) :: edge2vert_coeff_cc(1:nproma,1:patch%nblks_v,1:6)
    TYPE(t_cartesian_coordinates), INTENT(INOUT) :: edge2vert_coeff_cc_t(1:nproma,1:patch%nblks_e,1:2)
    REAL(wp), INTENT(INOUT) :: dist_cell2edge(1:nproma,1:patch%nblks_e,1:2)
    REAL(wp), INTENT(INOUT) :: fixed_vol_norm(1:nproma,patch%nblks_c)
    REAL(wp), INTENT(INOUT) :: variable_vol_norm(1:nproma,1:patch%nblks_c,1:3)
    REAL(wp), INTENT(INOUT) :: variable_dual_vol_norm(1:nproma,1:patch%nblks_e,1:6)
!
!Local variables
!
    REAL(wp), ALLOCATABLE :: prime_edge_length( :, : )
    REAL(wp), ALLOCATABLE :: dual_edge_length ( :, : )
!     REAL(wp), ALLOCATABLE :: cell_area( :, : )
!     REAL(wp), ALLOCATABLE :: dual_cell_area ( :, : )

    TYPE(t_subset_range), POINTER :: owned_edges         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_cells         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_verts         ! these are the owned entities
    TYPE(t_cartesian_coordinates) :: vertex_position, cell_center, edge_center
    TYPE(t_cartesian_coordinates) :: dist_vector
    TYPE(t_cartesian_coordinates), POINTER :: dual_edge_middle(:,:)

    TYPE(t_cartesian_coordinates) :: coriolis_cartesian_coordinates
    TYPE(t_geographical_coordinates) :: coriolis_geo_coordinates, geo_coordinates
    REAL(wp) :: basin_center_lat_rad, basin_height_rad

    REAL(wp) :: norm, orientation, length
    REAL(wp) :: inverse_sphere_radius

    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: cell_1_index, cell_1_block, cell_2_index, cell_2_block
    INTEGER :: vertex_1_index, vertex_1_block, vertex_2_index, vertex_2_block
    !-----------------------------------------------------------------------
!     REAL(wp) :: dist_cell2edge(nproma, patch%nblks_e,2)
!     TYPE(t_cartesian_coordinates) :: check_v1(nproma, patch%nblks_v, 6)
!     TYPE(t_cartesian_coordinates) :: check_v2(nproma, patch%nblks_e, 2)
!     REAL(wp) :: max_diff, max_val
    !-----------------------------------------------------------------------
    inverse_sphere_radius = 1.0_wp / grid_sphere_radius

    owned_edges => patch%edges%owned
    owned_cells => patch%cells%owned
    owned_verts => patch%verts%owned

    edge2vert_coeff_cc(:,:,:)%x(1) = 0.0_wp
    edge2vert_coeff_cc(:,:,:)%x(2) = 0.0_wp
    edge2vert_coeff_cc(:,:,:)%x(3) = 0.0_wp

    edge2vert_coeff_cc_t(:,:,:)%x(1) = 0.0_wp
    edge2vert_coeff_cc_t(:,:,:)%x(2) = 0.0_wp
    edge2vert_coeff_cc_t(:,:,:)%x(3) = 0.0_wp


    !-------------------------------------------
    ! compute some basic distances
    ! this is required if the cartesian distance is used
    ! instead of the spherical
    !
    ! computes_dist_cell2edge( patch, intp_2D_coeff)
    !
    ALLOCATE( prime_edge_length( nproma, patch%nblks_e))
    ALLOCATE( dual_edge_length ( nproma, patch%nblks_e))
!     ALLOCATE( cell_area        ( nproma, patch%nblks_c))
!     ALLOCATE( dual_cell_area   ( nproma, patch%nblks_v))

    IF ( MID_POINT_DUAL_EDGE ) THEN
      dual_edge_middle => patch%edges%cartesian_dual_middle
    ELSE
      dual_edge_middle => patch%edges%cartesian_center
    ENDIF

    ! get the areas on a unit sphere
!     cell_area(:,:)      = patch%cells%area(:,:)      * inverse_earth_radius * inverse_earth_radius
!     dual_cell_area(:,:) = patch%verts%dual_area(:,:) * inverse_earth_radius * inverse_earth_radius

    IF (LARC_LENGTH) THEN

      ! we just need to get them from the grid
      ! NOTE:  these are earth's distances, translate on a unit sphere
      dist_cell2edge(:,:,:) = &
        & patch%edges%edge_cell_length(:,:,:) * inverse_sphere_radius
      prime_edge_length(:,:) = &
        & patch%edges%primal_edge_length(:,:) * inverse_sphere_radius
      dual_edge_length(:,:) = &
        & patch%edges%dual_edge_length(:,:) * inverse_sphere_radius

    ELSE

      ! calcultate cartesian distance
      prime_edge_length(:,:) = 0.0_wp
      dual_edge_length(:,:) = 0.0_wp
      dist_cell2edge(:,:,:) =  0.0_wp

      DO edge_block = owned_edges%start_block, owned_edges%end_block
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index = start_index, end_index

          !----------------------------------------
          ! calculate the cartesian edge length
          vertex_1_index = patch%edges%vertex_idx(edge_index, edge_block, 1)
          vertex_1_block = patch%edges%vertex_blk(edge_index, edge_block, 1)
          vertex_2_index = patch%edges%vertex_idx(edge_index, edge_block, 2)
          vertex_2_block = patch%edges%vertex_blk(edge_index, edge_block, 2)

          dist_vector%x = &
            & patch%verts%cartesian(vertex_1_index, vertex_1_block)%x - &
            & patch%verts%cartesian(vertex_2_index, vertex_2_block)%x

            prime_edge_length(edge_index,edge_block) = &
              & SQRT(SUM((  dist_vector%x *  dist_vector%x)))
          !----------------------------------------

          !----------------------------------------
          ! calculate the cartesian distance of the edge center to the cell center
          DO neigbor = 1,2

            dist_cell2edge(edge_index,edge_block,neigbor) = 0.0_wp

            cell_index = patch%edges%cell_idx(edge_index,edge_block,neigbor)
            cell_block = patch%edges%cell_blk(edge_index,edge_block,neigbor)

            IF (cell_block > 0) THEN
              dist_vector%x = &
                & patch%edges%cartesian_center(edge_index,edge_block)%x - &
                & patch%cells%cartesian_center(cell_index,cell_block)%x

              dist_cell2edge(edge_index,edge_block,neigbor) = &
                & SQRT(SUM((  dist_vector%x *  dist_vector%x)))
            ENDIF

          ENDDO ! neigbor = 1,2
          !----------------------------------------

          !----------------------------------------
          ! calculate the cartesian dual edge length
          cell_1_index = patch%edges%cell_idx(edge_index, edge_block, 1)
          cell_1_block = patch%edges%cell_blk(edge_index, edge_block, 1)
          cell_2_index = patch%edges%cell_idx(edge_index, edge_block, 2)
          cell_2_block = patch%edges%cell_blk(edge_index, edge_block, 2)

          IF (cell_1_block > 0 .AND. cell_2_block > 0) THEN
            dist_vector%x = &
              & patch%cells%cartesian_center(cell_1_index, cell_1_block)%x - &
              & patch%cells%cartesian_center(cell_2_index, cell_2_block)%x

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
      CALL sync_patch_array(SYNC_E, patch, dist_cell2edge(:,:,1))
      CALL sync_patch_array(SYNC_E, patch, dist_cell2edge(:,:,2))
      CALL sync_patch_array(SYNC_E, patch, prime_edge_length(:,:))
      CALL sync_patch_array(SYNC_E, patch, dual_edge_length(:,:))
    ENDIF
    ! distances have been computed
    !-------------------------------------------

    !-------------------------------------------
    ! compute:
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
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        cell_center%x = patch%cells%cartesian_center(cell_index, cell_block)%x
        fixed_vol_norm(cell_index,cell_block) = 0.0_wp

        !-------------------------------
        DO neigbor=1,patch%cell_type

          edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x = 0.0_wp
          variable_vol_norm(cell_index, cell_block, neigbor) =  0.0_wp

          edge_index = patch%cells%edge_idx(cell_index, cell_block, neigbor)
          edge_block = patch%cells%edge_blk(cell_index, cell_block, neigbor)

          IF (edge_block > 0 ) THEN
            ! we have an edge
            dist_vector%x = &
              & patch%edges%cartesian_center(edge_index,edge_block)%x - &
              & cell_center%x

            norm  = SQRT(SUM( dist_vector%x * dist_vector%x))

            ! compute edge2cell_coeff_cc
            edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x =  &
              & dist_vector%x *                                             &
              & prime_edge_length(edge_index,edge_block) *                  &
              & patch%cells%edge_orientation(cell_index,cell_block,neigbor)! / &
              ! & cell_area(cell_index, cell_block)
              ! Note: here we do not divide by the cell area !

            fixed_vol_norm(cell_index,cell_block) = &
              & fixed_vol_norm(cell_index,cell_block) + &
              & 0.5_wp * norm * prime_edge_length(edge_index,edge_block)

            variable_vol_norm(cell_index, cell_block, neigbor) = &
              & 0.5_wp * norm * prime_edge_length(edge_index,edge_block)

          ENDIF !(edge_block > 0 )

        ENDDO !neigbor=1,patch%cell_type
        !-------------------------------

      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
    !-------------------
    ! sync the results
    CALL sync_patch_array(SYNC_C, patch, fixed_vol_norm(:,:))
    DO neigbor=1,patch%cell_type
      CALL sync_patch_array(SYNC_C, patch, edge2cell_coeff_cc(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_C, patch, edge2cell_coeff_cc(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_C, patch, edge2cell_coeff_cc(:,:,neigbor)%x(3))
      CALL sync_patch_array(SYNC_C, patch, variable_vol_norm(:,:,neigbor))
    ENDDO
    !-------------------

    !-------------------------------------------
    ! compute:
    !   edge2cell_coeff_cc_t

    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge2cell_coeff_cc_t(edge_index, edge_block, 2)%x = 0.0_wp
        edge_center%x = patch%edges%cartesian_center(edge_index, edge_block)%x

        DO neigbor=1,2

          edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x = 0.0_wp
          cell_index = patch%edges%cell_idx(edge_index, edge_block, neigbor)
          cell_block = patch%edges%cell_blk(edge_index, edge_block, neigbor)

          IF (cell_block > 0) THEN

            dist_vector%x =  edge_center%x -                             &
              patch%cells%cartesian_center(cell_index, cell_block)%x

            orientation = DOT_PRODUCT(dist_vector%x, &
              & patch%edges%primal_cart_normal(edge_index, edge_block)%x)
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
      CALL sync_patch_array(SYNC_E, patch, edge2cell_coeff_cc_t(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch, edge2cell_coeff_cc_t(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch, edge2cell_coeff_cc_t(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2
    !   edge2cell_coeff_cc_t is computed
    !-------------------------------------------

    !-------------------------------------------
    ! compute:
    !   edge2vert_coeff_cc
    !   variable_dual_vol_norm
    variable_dual_vol_norm(:,:,:)  = 0.0_wp
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index

        vertex_position%x = patch%verts%cartesian(vertex_index, vertex_block)%x

        DO neigbor=1, 6 ! we have to change this to accomodate the dual grid

          variable_dual_vol_norm(vertex_index, vertex_block, neigbor) = 0.0_wp

          edge_index = patch%verts%edge_idx(vertex_index, vertex_block, neigbor)
          edge_block = patch%verts%edge_blk(vertex_index, vertex_block, neigbor)

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
               & patch%edges%primal_cart_normal(edge_index, edge_block)%x)
            IF (orientation < 0.0_wp) dist_vector%x = - dist_vector%x

              edge2vert_coeff_cc(vertex_index, vertex_block, neigbor)%x = &
              & dist_vector%x                                *                    &
              & dual_edge_length(edge_index, edge_block) !    /                    &
              !& dual_cell_area(vertex_index, vertex_block)

              variable_dual_vol_norm(vertex_index, vertex_block, neigbor) = &
              & 0.5_wp * dual_edge_length(edge_index, edge_block) * length

          ENDIF !(edge_block > 0) THEN

        ENDDO !neigbor=1,6

      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
    !-------------------
    ! sync the results
    DO neigbor=1,6
      CALL sync_patch_array(SYNC_V, patch, edge2vert_coeff_cc(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_V, patch, edge2vert_coeff_cc(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_V, patch, edge2vert_coeff_cc(:,:,neigbor)%x(3))
      CALL sync_patch_array(SYNC_V, patch, variable_dual_vol_norm(:,:, neigbor))
    ENDDO ! neigbor=1,6
    ! edge2vert_coeff_cc
    ! variable_dual_vol_norm
    !   are computed
    !----------------------------------------------------

    !----------------------------------------------------
    ! compute:
    !   edge2vert_coeff_cc_t
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge_center%x = dual_edge_middle(edge_index, edge_block)%x

        DO neigbor=1,2

          edge2vert_coeff_cc_t(edge_index, edge_block, neigbor)%x = 0.0_wp

          vertex_index = patch%edges%vertex_idx(edge_index, edge_block, neigbor)
          vertex_block = patch%edges%vertex_blk(edge_index, edge_block, neigbor)

          edge2vert_coeff_cc_t(edge_index, edge_block, neigbor)%x =              &
            & (edge_center%x - patch%verts%cartesian(vertex_index, vertex_block)%x) * &
            & patch%edges%system_orientation(edge_index, edge_block)                / &
            & prime_edge_length(edge_index, edge_block)

        ENDDO !neigbor=1,2

      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
    !-------------------
    ! sync the results
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch, edge2vert_coeff_cc_t(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch, edge2vert_coeff_cc_t(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch, edge2vert_coeff_cc_t(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2
    ! edge2vert_coeff_cc_t is computed
    !----------------------------------------------------

    !----------------------------------------------------
    ! recalculate the coriolis coefficient
    ! It is required if we use the middle of the dual_edge_length
    IF (MID_POINT_DUAL_EDGE) THEN

      IF (CORIOLIS_TYPE == full_coriolis) THEN

        DO edge_block = owned_edges%start_block, owned_edges%end_block
          CALL get_index_range(owned_edges, edge_block, start_index, end_index)
          DO edge_index = start_index, end_index

             coriolis_geo_coordinates = cc2gc(dual_edge_middle(edge_index,edge_block))
             patch%edges%f_e(edge_index,edge_block) = &
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

          patch%edges%f_e(edge_index,edge_block) =  2.0_wp * grid_angular_velocity * &
            & ( sin(basin_center_lat_rad) + (cos(basin_center_lat_rad) / &
            &   grid_sphere_radius) * length)

          ENDDO
        ENDDO

      ENDIF !(CORIOLIS_TYPE==full_coriolis)
    ENDIF ! (MID_POINT_DUAL_EDGE)
    !-------------------
    ! sync patch%edges%f_e
    CALL sync_patch_array(SYNC_E, patch, patch%edges%f_e)


    !----------------------------------------------------
!    DEALLOCATE( prime_edge_length)

    DEALLOCATE( dual_edge_length )
!     DEALLOCATE( cell_area )
!     DEALLOCATE( dual_cell_area )
    !---------------------------------------------------------
!     RETURN
    !---------------------------------------------------------
    ! checks
!
!      !---------------------------------------------------------
!      check_v1 = intp_2D_coeff%edge2vert_coeff_cc
!      check_v2 = intp_2D_coeff%edge2vert_coeff_cc_t
!      !---------------------------------------------------------
!
!      CALL init_scalar_product_oce( patch, intp_2D_coeff )
!
!      !---------------------------------------------------------
!      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc(:,:,:)%x(1) - &
!        &  check_v1(:,:,:)%x(1) ))
!      max_val  =  MAXVAL(ABS( check_v1(:,:,:)%x(1)))
!      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(1)=", max_diff, max_val
!
!      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc(:,:,:)%x(2) - &
!        & check_v1(:,:,:)%x(2) ))
!      max_val  =  MAXVAL(ABS( check_v1(:,:,:)%x(2)))
!      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(2)=", max_diff, max_val
!
!      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc(:,:,:)%x(3) - &
!        & check_v1(:,:,:)%x(3) ))
!      max_val  =  MAXVAL(ABS( check_v1(:,:,:)%x(3)))
!      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(3)=", max_diff, max_val
!      !---------------------------------------------------------
!
!      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc_t(:,:,:)%x(1) - &
!        &  check_v2(:,:,:)%x(1) ))
!      max_val  =  MAXVAL(ABS( check_v2(:,:,:)%x(1)))
!      Write(0,*) "max diff of edge2vert_coeff_cc_t(:,:,:,:)%x(1)=", max_diff, max_val
!
!      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc_t(:,:,:)%x(2) - &
!        & check_v2(:,:,:)%x(2) ))
!      max_val  =  MAXVAL(ABS( check_v2(:,:,:)%x(2)))
!      Write(0,*) "max diff of edge2vert_coeff_cc_t(:,:,:,:)%x(2)=", max_diff, max_val
!
!      max_diff = MAXVAL(ABS(intp_2D_coeff%edge2vert_coeff_cc_t(:,:,:)%x(3) - &
!        & check_v2(:,:,:)%x(3) ))
!      max_val  =  MAXVAL(ABS( check_v2(:,:,:)%x(3)))
!      Write(0,*) "max diff of edge2vert_coeff_cc_t(:,:,:,:)%x(3)=", max_diff, max_val
!

  END SUBROUTINE par_init_coeff_2D
  !-------------------------------------------------------------------------

  !> Initialize 3D expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !! Parellelized by Leonidas Linardakis 2012-3


  SUBROUTINE copy_2D_to_3D_coeff( patch,               &
                               & ocean_coeff,          &
                               & edge2cell_coeff_cc,   &
                               & edge2cell_coeff_cc_t, &
                               & edge2vert_coeff_cc,   &
                               & edge2vert_coeff_cc_t, &
                               & dist_cell2edge,       &
                               & fixed_vol_norm,       &
                               & variable_vol_norm)!,    &
                               !& variable_dual_vol_norm)
    TYPE(t_patch)    , TARGET, INTENT(INOUT)     :: patch
    TYPE(t_operator_coeff), INTENT(inout) :: ocean_coeff
    TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2cell_coeff_cc(1:nproma,1:patch%nblks_c,1:patch%cell_type)
    TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2cell_coeff_cc_t(1:nproma,1:patch%nblks_e,1:2)
    TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2vert_coeff_cc(1:nproma,1:patch%nblks_v,1:6)
    TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2vert_coeff_cc_t(1:nproma,1:patch%nblks_e,1:2)
    REAL(wp), INTENT(IN) :: dist_cell2edge(1:nproma,1:patch%nblks_e,1:2)
    REAL(wp), INTENT(IN) :: fixed_vol_norm(1:nproma,patch%nblks_c)
    REAL(wp), INTENT(IN) :: variable_vol_norm(1:nproma,1:patch%nblks_c,1:3)
    !REAL(wp), INTENT(IN) :: variable_dual_vol_norm(1:nproma,1:patch%nblks_v,1:6)
    !
    !Local variables
    !
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: all_verts

    INTEGER :: edge_block, cell_block, vertex_block, level, neigbor
    INTEGER :: je,i_startidx_e, i_endidx_e

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
    ! p_patch%edges%f_e(:, :) is already calculated in par_init_scalar_product_oce    !
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
          dist_cell2edge(:,edge_block,1)
        ocean_coeff%dist_cell2edge(:,level,edge_block,2) = &
          dist_cell2edge(:,edge_block,2)
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
            edge2cell_coeff_cc_t(je,edge_block,1)%x
          ocean_coeff%edge2cell_coeff_cc_t(je,level,edge_block,2)%x = &
            edge2cell_coeff_cc_t(je,edge_block,2)%x

          ocean_coeff%edge2vert_coeff_cc_t(je,level,edge_block,1)%x = &
            edge2vert_coeff_cc_t(je,edge_block,1)%x
          ocean_coeff%edge2vert_coeff_cc_t(je,level,edge_block,2)%x = &
            edge2vert_coeff_cc_t(je,edge_block,2)%x

        ENDDO
      ENDDO
    ENDDO

    ! on cells
    DO cell_block = all_cells%start_block, all_cells%end_block
      DO level = 1, n_zlev

       ocean_coeff%fixed_vol_norm(:,level,cell_block) = &
         fixed_vol_norm(:,cell_block)

       DO neigbor=1,patch%cell_type

         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(1) = &
           & edge2cell_coeff_cc(:,cell_block,neigbor)%x(1)
         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(2) = &
           & edge2cell_coeff_cc(:,cell_block,neigbor)%x(2)
         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(3) = &
           & edge2cell_coeff_cc(:,cell_block,neigbor)%x(3)

         ocean_coeff%variable_vol_norm(:,level,cell_block,neigbor) = &
           & variable_vol_norm(:,cell_block,neigbor)

        ENDDO ! neigbor=1,patch%cell_type

      ENDDO  !  level = 1, n_zlev
    ENDDO ! cell_block

    ! on verts
    DO vertex_block = all_verts%start_block, all_verts%end_block
      DO level = 1, n_zlev
        DO neigbor=1,6

          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(1) = &
            & edge2vert_coeff_cc(:,vertex_block,neigbor)%x(1)
          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(2) = &
            & edge2vert_coeff_cc(:,vertex_block,neigbor)%x(2)
          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(3) = &
            & edge2vert_coeff_cc(:,vertex_block,neigbor)%x(3)

        ENDDO ! neigbor=1,patch%cell_type
      ENDDO  !  level = 1, n_zlev
    ENDDO ! vertex_block
    !---------------------------------------------------------

  END SUBROUTINE copy_2D_to_3D_coeff
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !> Initialize expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !!
  SUBROUTINE par_apply_boundary2coeffs( patch, p_patch_3D, ocean_coeff)
    ! !
    TYPE(t_patch), TARGET,  INTENT(inout) :: patch
    TYPE(t_patch_3D_oce ),TARGET, INTENT(INOUT)   :: p_patch_3D
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

    !all_cells   => patch%cells%all
    !all_edges   => patch%edges%all
    !owned_verts => patch%verts%owned

    all_cells   => p_patch_3D%p_patch_2D(1)%cells%all
    all_edges   => p_patch_3D%p_patch_2D(1)%edges%all
    owned_verts => p_patch_3D%p_patch_2D(1)%verts%owned

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

          !IF(v_base%lsm_oce_e(je,jk,jb) /= sea) THEN
          IF ( p_patch_3D%lsm_oce_e(je,jk,jb) /= sea ) THEN
            ocean_coeff%grad_coeff          (je,jk,jb) = 0.0_wp
            ocean_coeff%edge2cell_coeff_cc_t(je,jk,jb,1)%x(1:3) = 0.0_wp
            ocean_coeff%edge2cell_coeff_cc_t(je,jk,jb,2)%x(1:3) = 0.0_wp
            ocean_coeff%edge2vert_coeff_cc_t(je,jk,jb,1)%x(1:3) = 0.0_wp
            ocean_coeff%edge2vert_coeff_cc_t(je,jk,jb,2)%x(1:3) = 0.0_wp
          ENDIF

        ENDDO
      END DO
    END DO
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
            !IF ( v_base%lsm_oce_e(ile,jk,ibe) /= sea) THEN
            IF ( p_patch_3D%lsm_oce_e(ile,jk,ibe) /= sea ) THEN
              ocean_coeff%div_coeff(jc,jk,jb,je) = 0.0_wp
              ocean_coeff%edge2cell_coeff_cc(jc,jk,jb,je)%x(1:3) = 0.0_wp
            ENDIF
          ENDDO ! je = 1, patch%cells%num_edges(jc,jb)
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

            !IF ( v_base%lsm_oce_e(ile,jk,ibe) == sea) THEN
            IF ( p_patch_3D%lsm_oce_e(ile,jk,ibe) == sea ) THEN
              i_v_ctr(jv,jk,jb)=i_v_ctr(jv,jk,jb)+1
            ELSEIF ( p_patch_3D%lsm_oce_e(ile,jk,ibe) == boundary ) THEN

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
!---------------------------------------------------------------------------------
            !Modified area calculation
            vertex_cc = patch%verts%cartesian(jv,jb)
            zarea_fraction=0.0_wp
            DO jev = 1, patch%verts%num_edges(jv,jb)
              ! get line and block indices of edge jev around vertex jv
              ile = patch%verts%edge_idx(jv,jb,jev)
              ibe = patch%verts%edge_blk(jv,jb,jev)
              !get neighbor cells
              icell_idx_1 = patch%edges%cell_idx(ile,ibe,1)
              icell_idx_2 = patch%edges%cell_idx(ile,ibe,2)
              icell_blk_1 = patch%edges%cell_blk(ile,ibe,1)
              icell_blk_2 = patch%edges%cell_blk(ile,ibe,2)
              ! cell1_cc    = gc2cc(patch%cells%center(icell_idx_1,icell_blk_1))
              ! cell2_cc    = gc2cc(patch%cells%center(icell_idx_2,icell_blk_2))
              cell1_cc%x  = patch%cells%cartesian_center(icell_idx_1,icell_blk_1)%x
              cell2_cc%x  = patch%cells%cartesian_center(icell_idx_2,icell_blk_2)%x

              !Check, if edge is sea or boundary edge and take care of dummy edge
              !edge with indices ile, ibe is sea edge
              !Add up for wet dual area.
              !IF ( v_base%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN
              IF ( p_patch_3D%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN
                zarea_fraction = zarea_fraction  &
                  & + triangle_area(cell1_cc, vertex_cc, cell2_cc)
                ! edge with indices ile, ibe is boundary edge
              ELSE IF ( p_patch_3D%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
                zarea_fraction = zarea_fraction  &
                  & + 0.5_wp*triangle_area(cell1_cc, vertex_cc, cell2_cc)
              END IF
            END DO
!IF(zarea_fraction*(re*re)/= patch%verts%dual_area(jv,jb))&
!&write(*,*)'areas',zarea_fraction*(re*re), patch%verts%dual_area(jv,jb)
            ! no division by zero
            !IF (zarea_fraction /= 0.0_wp) THEN
            z_area_scaled   =  zarea_fraction*(earth_radius*earth_radius)
            !patch%verts%dual_area(jv,jb)! !SUM(ocean_coeff%variable_dual_vol_norm(jv,jk,jb,:))
            !ENDIF
!---------------------------------------------------------------------------------------------
          DO je = 1, boundary_counter
!-original
!             ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(je) )=&
!               & 0.5_wp*z_orientation(je) * &
!               & patch%edges%primal_edge_length(ibnd_edge_idx(je),ibnd_edge_blk(je))

            ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(je) )=&
              & 0.5_wp*patch%edges%system_orientation(ibnd_edge_idx(je),ibnd_edge_blk(je)) * &
              & patch%edges%primal_edge_length(ibnd_edge_idx(je),ibnd_edge_blk(je))

          ENDDO

          IF(z_area_scaled/=0.0_wp)THEN
            ocean_coeff%rot_coeff(jv,jk,jb,:)&
            &=ocean_coeff%rot_coeff(jv,jk,jb,:)/z_area_scaled
          ELSE
            ocean_coeff%rot_coeff(jv,jk,jb,:)=0.0_wp
          ENDIF
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
            !IF ( v_base%lsm_oce_e(ile,jk,ibe) /= sea) THEN
             IF ( p_patch_3D%lsm_oce_e(ile,jk,ibe) /= sea) THEN
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
    !       it does not need to be synced
    DO jb = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, jb, i_startidx_v, i_endidx_v)
      DO jk = 1, n_zlev
!CDIR nextscalar
        DO jv = i_startidx_v, i_endidx_v

          zarea_fraction   = 0.0_wp
          z_area_scaled    = 0.0_wp
          !IF ( ocean_coeff%bnd_edges_per_vertex(jv,jk,jb) == 0 ) THEN
          IF ( i_v_ctr(jv,jk,jb) == patch%verts%num_edges(jv,jb) ) THEN
            z_area_scaled = patch%verts%dual_area(jv,jb)/(earth_radius*earth_radius)
            !SUM(ocean_coeff%variable_dual_vol_norm(jv,jk,jb,:))

!TODO ram   !Final coefficient calculation
!            DO jev = 1, patch%verts%num_edges(jv,jb)
!
!              IF(z_area_scaled/=0.0_wp)THEN
!                ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)&
!                  & =ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)!/z_area_scaled
!              ELSE
!                ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)=0.0_wp
!              ENDIF
!            END DO

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
!               cell1_cc    = gc2cc(patch%cells%center(icell_idx_1,icell_blk_1))
!               cell2_cc    = gc2cc(patch%cells%center(icell_idx_2,icell_blk_2))
              cell1_cc%x  = patch%cells%cartesian_center(icell_idx_1,icell_blk_1)%x
              cell2_cc%x  = patch%cells%cartesian_center(icell_idx_2,icell_blk_2)%x

              !Check, if edge is sea or boundary edge and take care of dummy edge
              !edge with indices ile, ibe is sea edge
              !Add up for wet dual area.
              ! Note that this should be modified.
              !   sea_boundary means an open boundary
              !   boundary means that only the sea cell are should be added
              !IF ( v_base%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN
              IF ( p_patch_3D%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN
                zarea_fraction = zarea_fraction  &
                  & + triangle_area(cell1_cc, vertex_cc, cell2_cc)
                ! edge with indices ile, ibe is boundary edge                
              !ELSE IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
              ELSE IF ( p_patch_3D%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
                zarea_fraction = zarea_fraction  &
                  & + 0.5_wp*triangle_area(cell1_cc, vertex_cc, cell2_cc)
              END IF
            END DO
            ! no division by zero
            !IF (zarea_fraction /= 0.0_wp) THEN
            z_area_scaled = zarea_fraction !patch%verts%dual_area(jv,jb)/(re*re)!!SUM(ocean_coeff%variable_dual_vol_norm(jv,jk,jb,:))
            !ENDIF
          ENDIF !( i_v_ctr(jv,jk,jb) == patch%verts%num_edges(jv,jb) )

          !Final coefficient calculation

! !CDIR nextscalar
          IF(z_area_scaled/=0.0_wp)THEN
            DO jev = 1, patch%verts%num_edges(jv,jb)
              ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)&
                & =ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)/z_area_scaled
            END DO
          ELSE
            DO jev = 1, patch%verts%num_edges(jv,jb)
                ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)=0.0_wp
            END DO
          ENDIF
!
!          ENDIF !( i_v_ctr(jv,jk,jb) == patch%verts%num_edges(jv,jb) )

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
  SUBROUTINE init_geo_factors_oce_3d( p_patch, p_coeff )
    !
    IMPLICIT NONE
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch
    TYPE(t_operator_coeff),     INTENT(inout) :: p_coeff
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

    REAL(wp) :: z_sync_c(nproma,n_zlev,p_patch%nblks_c)
    !REAL(wp) :: z_sync_e(nproma,n_zlev,p_patch%nblks_e)
    REAL(wp) :: z_sync_v(nproma,n_zlev,p_patch%nblks_v)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:init_geo_factors_oce_3d')
    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')
    i_nchdom   = MAX(1,p_patch%n_childdom)

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

    ! 1) coefficients for divergence
    rl_start = 1
    rl_end = min_rlcell

    ! values for the blocking
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jc,i_startidx,i_endidx,ile,ibe)
    DO jk=1,n_zlev
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
          & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          !         ile1 = p_patch%cells%edge_idx(jc,jb,1)
          !         ibe1 = p_patch%cells%edge_blk(jc,jb,1)
          !         ile2 = p_patch%cells%edge_idx(jc,jb,2)
          !         ibe2 = p_patch%cells%edge_blk(jc,jb,2)
          !         ile3 = p_patch%cells%edge_idx(jc,jb,3)
          !         ibe3 = p_patch%cells%edge_blk(jc,jb,3)
          !         cell_area =  0.25_wp&
          !         & *( p_patch%edges%primal_edge_length(ile1,ibe1)*p_patch%edges%dual_edge_length(ile1,ibe1)&
          !         &   +p_patch%edges%primal_edge_length(ile2,ibe2)*p_patch%edges%dual_edge_length(ile2,ibe2)&
          !         &   +p_patch%edges%primal_edge_length(ile3,ibe3)*p_patch%edges%dual_edge_length(ile3,ibe3))

          DO je = 1, i_cell_type

            ile = p_patch%cells%edge_idx(jc,jb,je)
            ibe = p_patch%cells%edge_blk(jc,jb,je)

            p_coeff%div_coeff(jc,jk,jb,je) =                &
              & p_patch%edges%primal_edge_length(ile,ibe) * &
              & p_patch%cells%edge_orientation(jc,jb,je)  / &
              & p_patch%cells%area(jc,jb)
          ENDDO !edge loop
          !write(1234,*)'div coeff 3D',jk,jc,jb,p_coeff%div_coeff(jc,jk,jb,:)
        ENDDO !idx loop
      END DO !block loop
    END DO
!$OMP END DO

    ! 2) coefficients for curl
    rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    rl_end = min_rlvert_int

    ! values for the blocking
    i_startblk = p_patch%verts%start_blk(rl_start,1)
    i_endblk   = p_patch%verts%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jv,i_startidx,i_endidx,ile,ibe)
    DO jk=1,n_zlev
      DO jb = i_startblk, i_endblk

        CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO jv = i_startidx, i_endidx
          DO je = 1, p_patch%verts%num_edges(jv,jb)

            ile = p_patch%verts%edge_idx(jv,jb,je)
            ibe = p_patch%verts%edge_blk(jv,jb,je)

            p_coeff%rot_coeff(jv,jk,jb,je) =                &
              & p_patch%edges%dual_edge_length(ile,ibe) * &
              & p_patch%verts%edge_orientation(jv,jb,je)!/ &
            !&    p_patch%verts%dual_area(jv,jb)
          ENDDO
        ENDDO
      END DO
    END DO
!$OMP END DO

    ! 3) coefficients for nabla2_scalar
    rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    rl_end = min_rlcell_int

    ! values for the blocking
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
    !!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
    !!$OMP    ilc2,ibc2,ilnc,ibnc)
    ! !       DO jk=1,n_zlev
    ! !       DO jb = i_startblk, i_endblk
    ! !
    ! !         CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
    ! !                            i_startidx, i_endidx, rl_start, rl_end)
    ! !
    ! !         DO je = 1, i_cell_type
    ! !           DO jc = i_startidx, i_endidx
    ! !
    ! !             ile = p_patch%cells%edge_idx(jc,jb,je)
    ! !             ibe = p_patch%cells%edge_blk(jc,jb,je)
    ! !
    ! !             ilc1 = p_patch%edges%cell_idx(ile,ibe,1)
    ! !             ibc1 = p_patch%edges%cell_blk(ile,ibe,1)
    ! !             ilc2 = p_patch%edges%cell_idx(ile,ibe,2)
    ! !             ibc2 = p_patch%edges%cell_blk(ile,ibe,2)
    ! !
    ! !             IF (jc == ilc1 .AND. jb == ibc1) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 p_coeff%geofac_n2s(jc,1,jb)     =  &
    ! !                 &  p_coeff%geofac_n2s(jc,1,jb)  -  &
    ! !                 &  p_coeff%geofac_div(jc,je,jb) /  &
    ! !                 &  p_patch%edges%dual_edge_length(ile,ibe)
    ! !             ENDIF
    ! !           ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
    ! !             IF (i_cell_type == 3) THEN
    ! !               p_coeff%geofac_n2s(jc,1,jb)       =  &
    ! !                 &  p_coeff%geofac_n2s(jc,1,jb)  +  &
    ! !                 &  p_coeff%geofac_div(jc,je,jb) /  &
    ! !                 &  p_patch%edges%dual_edge_length(ile,ibe)
    ! !             ENDIF
    ! !           ENDIF
    ! !           DO ic = 1, i_cell_type
    ! !             ilnc = p_patch%cells%neighbor_idx(jc,jb,ic)
    ! !             ibnc = p_patch%cells%neighbor_blk(jc,jb,ic)
    ! !             IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 p_coeff%geofac_n2s(jc,ic+1,jb)     = &
    ! !                   &  p_coeff%geofac_n2s(jc,ic+1,jb)- &
    ! !                   &  p_coeff%geofac_div(jc,je,jb)  / &
    ! !                   &  p_patch%edges%dual_edge_length(ile,ibe)
    ! !               ENDIF
    ! !             ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 p_coeff%geofac_n2s(jc,ic+1,jb)     = &
    ! !                   &  p_coeff%geofac_n2s(jc,ic+1,jb)+ &
    ! !                   &  p_coeff%geofac_div(jc,je,jb)  / &
    ! !                   &  p_patch%edges%dual_edge_length(ile,ibe)
    ! !               ENDIF
    ! !             ENDIF
    ! !           ENDDO
    ! !
    ! !           ! To ensure that dummy edges have a factor of 0:
    ! !           IF (je > p_patch%cells%num_edges(jc,jb)) THEN
    ! !             p_coeff%geofac_n2s(jc,je+1,jb) = 0._wp
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
    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch edges
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je)
    DO jk=1,n_zlev
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO je =  i_startidx, i_endidx
          ! compute inverse dual edge length (undefined for refin_ctrl=1)
          p_patch%edges%inv_dual_edge_length(je,jb) = &
            & 1._wp/p_patch%edges%dual_edge_length(je,jb)

          p_coeff%grad_coeff(je,jk,jb)&
            & = p_patch%edges%inv_dual_edge_length(je,jb)

        ENDDO
      END DO !block loop
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! synchronize all elements of p_coeff:
    CALL sync_patch_array(sync_e, p_patch, p_coeff%grad_coeff)
    DO ie = 1, i_cell_type

      z_sync_c(:,:,:) = p_coeff%div_coeff(:,:,:,ie)
      CALL sync_patch_array(sync_c, p_patch, z_sync_c(:,:,:))
      p_coeff%div_coeff(:,:,:,ie) = z_sync_c(:,:,:)
      !         z_sync_c(:,:,:) = p_coeff%n2s_coeff(:,:,:,ie)
      !         CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:,:))
      !         p_coeff%n2s_coeff(:,:,:,ie) = z_sync_c(:,:,:)
    END DO

    CALL sync_patch_array(sync_e, p_patch, p_patch%edges%inv_dual_edge_length(:,:))

    DO ie = 1, 9-i_cell_type
      z_sync_v(:,:,:) = p_coeff%rot_coeff(:,:,:,ie)
      CALL sync_patch_array(sync_v, p_patch, z_sync_v(:,:,:))
      p_coeff%rot_coeff(:,:,:,ie) = z_sync_v(:,:,:)
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
! !---------------------------------------------------------------------------------------
! !> Initialize expansion coefficients.
!   !!
!   !! @par Revision History
!   !! Peter Korn (2012-2)
!   !!
!   SUBROUTINE par_init_operator_coeff( patch, p_os, p_phys_param,ocean_coeff, intp_2D_coeff)
!     !
!     TYPE(t_patch),          INTENT(inout) :: patch
!     TYPE(t_hydro_ocean_state),INTENT(IN)  :: p_os
!     TYPE (t_ho_params),       INTENT(IN)  :: p_phys_param
!     TYPE(t_operator_coeff), INTENT(inout) :: ocean_coeff
!     TYPE(t_int_state),      INTENT(inout) :: intp_2D_coeff
! 
!     !INTEGER :: rl_start_e,rl_end_e,rl_start_c,rl_end_c
!     !INTEGER :: i_startblk_e, i_endblk_e,i_startidx_e,i_endidx_e
!     !INTEGER :: i_startblk_c, i_endblk_c,i_startidx_c,i_endidx_c
!     !INTEGER :: jk,je,jb,jc
!     !-----------------------------------------------------------------------
!     !-----------------------------------------------------------------------
! !     TYPE(t_cartesian_coordinates) :: check_v(nproma, n_zlev, patch%nblks_v, 6)
! !     REAL(wp) :: check_r(nproma, n_zlev, patch%nblks_c, 3)
! !     REAL(wp) :: max_diff, max_val
!     !-----------------------------------------------------------------------
!     ocean_coeff%dist_cell2edge(:,:,:,:) = 0.0_wp
! 
!     CALL par_init_scalar_product_oce(patch, intp_2D_coeff)
!     CALL copy_2D_to_3D_intp_coeff( patch, ocean_coeff, intp_2D_coeff)
!     CALL init_geo_factors_oce_3d ( patch, ocean_coeff )
!     !---------------------------------------------------------
!     CALL par_apply_boundary2coeffs(patch, ocean_coeff)
! 
! !     CALL update_diffusion_matrices( patch,                          &
! !                                      & p_os,                          &
! !                                      & p_phys_param,                  &
! !                                      & ocean_coeff%matrix_vert_diff_e,&
! !                                      & ocean_coeff%matrix_vert_diff_c)
! 
!     RETURN
!     !---------------------------------------------------------
!     ! checks
! !     check_v = ocean_coeff%edge2vert_coeff_cc
! !     check_r = ocean_coeff%variable_vol_norm
! !
! !     !---------------------------------------------------------
! !      CALL init_operator_coeff( patch, ocean_coeff )
! !
! !      max_diff = MAXVAL(ABS(ocean_coeff%edge2vert_coeff_cc(:,:,:,:)%x(1) - &
! !        &  check_v(:,:,:,:)%x(1) ))
! !      max_val  =  MAXVAL(ABS( check_v(:,:,:,:)%x(1)))
! !      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(1)=", max_diff, max_val
! !
! !      max_diff = MAXVAL(ABS(ocean_coeff%edge2vert_coeff_cc(:,:,:,:)%x(2) - &
! !        & check_v(:,:,:,:)%x(2) ))
! !      max_val  =  MAXVAL(ABS( check_v(:,:,:,:)%x(2)))
! !      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(2)=", max_diff, max_val
! !
! !      max_diff = MAXVAL(ABS(ocean_coeff%edge2vert_coeff_cc(:,:,:,:)%x(3) - &
! !        & check_v(:,:,:,:)%x(3) ))
! !      max_val  =  MAXVAL(ABS( check_v(:,:,:,:)%x(3)))
! !      Write(0,*) "max diff of edge2vert_coeff_cc(:,:,:,:)%x(3)=", max_diff, max_val
! !
! !      max_diff = MAXVAL(ABS(ocean_coeff%variable_vol_norm - check_r ))
! !      max_val  =  MAXVAL(ABS( check_r ))
! !      Write(0,*) "max diff of variable_vol_norm=", max_diff, max_val
! !
! !      STOP
! 
!   END SUBROUTINE par_init_operator_coeff
  !-------------------------------------------------------------------------
 
  !-------------------------------------------------------------------------
! !   !> Initialize 3D expansion coefficients.
! !   !!
! !   !! @par Revision History
! !   !! Peter Korn (2012-2)
! !   !! Parellelized by Leonidas Linardakis 2012-3
! !   SUBROUTINE copy_2D_to_3D_intp_coeff( patch, ocean_coeff, intp_2D_coeff)
! !     !
! !     TYPE(t_patch),  TARGET, INTENT(inout) :: patch
! !     TYPE(t_operator_coeff), INTENT(inout) :: ocean_coeff
! !     TYPE(t_int_state),      INTENT(inout) :: intp_2D_coeff
! ! 
! !     TYPE(t_subset_range), POINTER :: all_edges
! !     TYPE(t_subset_range), POINTER :: all_cells
! !     TYPE(t_subset_range), POINTER :: all_verts
! ! 
! !     INTEGER :: edge_block, cell_block, vertex_block, level, neigbor
! !     INTEGER :: je,jk, i_startidx_e, i_endidx_e
! ! 
! !     all_cells => patch%cells%all
! !     all_edges => patch%edges%all
! !     all_verts => patch%verts%all
! ! 
! !     !---------------------------------------------------------
! !     ! the following coefficients will be copied:
! !     !
! !     ! ocean_coeff%edge_position_cc(:,:,:)            on edges
! !     ! ocean_coeff%dist_cell2edge(:,:,:,1-2)          on edges
! !     !
! !     ! ocean_coeff%edge2cell_coeff_cc(:,:,:,1-3)%x    on cells
! !     ! ocean_coeff%fixed_vol_norm(:,:,:)              on cells
! !     ! ocean_coeff%variable_vol_norm(:,:,:,1-3)       on cells
! !     !
! !     ! ocean_coeff%edge2cell_coeff_cc_t(:,:,:,1-2)%x  on edges
! !     ! ocean_coeff%edge2vert_coeff_cc(:,:,:,1-6)%x   on verts
! !     ! ocean_coeff%edge2vert_coeff_cc_t(:,:,:,1-2)%x  on edges
! !     !
! !     ! p_patch%edges%f_e(:, :) is already calculated in par_init_scalar_product_oce    !
! !     !---------------------------------------------------------
! ! 
! ! 
! !     !---------------------------------------------------------
! !     ! calculate ocean_coeff%edge_position_cc(:,:,:) on edges
! !     ! this is the same as the 2D, just copy it
! !     ! it does not change, maybe turn it to 2D?
! !     DO edge_block = all_edges%start_block, all_edges%end_block
! !       DO level = 1, n_zlev
! !         ocean_coeff%edge_position_cc(:,level,edge_block) = &
! !           patch%edges%cartesian_center(:,edge_block)
! !       ENDDO
! !     ENDDO
! !     !---------------------------------------------------------
! ! 
! !     !---------------------------------------------------------
! !     ! calculate ocean_coeff%dist_cell2edge(:,:,:,1-2) on edges
! !     ! this is the same as the 2D, just copy it
! !     ! it does not change, maybe turn it to 2D?
! !     DO edge_block = all_edges%start_block, all_edges%end_block
! !       DO level = 1, n_zlev
! !         ocean_coeff%dist_cell2edge(:,level,edge_block,1) = &
! !           intp_2D_coeff%dist_cell2edge(:,edge_block,1)
! !         ocean_coeff%dist_cell2edge(:,level,edge_block,2) = &
! !           intp_2D_coeff%dist_cell2edge(:,edge_block,2)
! !       ENDDO
! !     ENDDO
! !     !---------------------------------------------------------
! ! 
! !     !---------------------------------------------------------
! !     ! For the rest of coefficients,
! !     ! copy the 2D and then apply_boundary2coeffs
! !     !
! !     ! on edges
! !     DO edge_block = all_edges%start_block, all_edges%end_block
! !       DO level = 1, n_zlev
! !         CALL get_index_range(all_edges, edge_block, i_startidx_e, i_endidx_e)
! !         DO je =  i_startidx_e, i_endidx_e
! ! 
! !           ocean_coeff%edge2cell_coeff_cc_t(je,level,edge_block,1)%x = &
! !             intp_2D_coeff%edge2cell_coeff_cc_t(je,edge_block,1)%x
! !           ocean_coeff%edge2cell_coeff_cc_t(je,level,edge_block,2)%x = &
! !             intp_2D_coeff%edge2cell_coeff_cc_t(je,edge_block,2)%x
! ! 
! !           ocean_coeff%edge2vert_coeff_cc_t(je,level,edge_block,1)%x = &
! !             intp_2D_coeff%edge2vert_coeff_cc_t(je,edge_block,1)%x
! !           ocean_coeff%edge2vert_coeff_cc_t(je,level,edge_block,2)%x = &
! !             intp_2D_coeff%edge2vert_coeff_cc_t(je,edge_block,2)%x
! ! 
! !         ENDDO
! !       ENDDO
! !     ENDDO
! ! 
! !     ! on cells
! !     DO cell_block = all_cells%start_block, all_cells%end_block
! !       DO level = 1, n_zlev
! ! 
! !        ocean_coeff%fixed_vol_norm(:,level,cell_block) = &
! !          intp_2D_coeff%fixed_vol_norm(:,cell_block)
! ! 
! !        DO neigbor=1,patch%cell_type
! ! 
! !          ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(1) = &
! !            & intp_2D_coeff%edge2cell_coeff_cc(:,cell_block,neigbor)%x(1)
! !          ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(2) = &
! !            & intp_2D_coeff%edge2cell_coeff_cc(:,cell_block,neigbor)%x(2)
! !          ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(3) = &
! !            & intp_2D_coeff%edge2cell_coeff_cc(:,cell_block,neigbor)%x(3)
! ! 
! !          ocean_coeff%variable_vol_norm(:,level,cell_block,neigbor) = &
! !            & intp_2D_coeff%variable_vol_norm(:,cell_block,neigbor)
! ! 
! !         ENDDO ! neigbor=1,patch%cell_type
! ! 
! !       ENDDO  !  level = 1, n_zlev
! !     ENDDO ! cell_block
! ! 
! !     ! on verts
! !     DO vertex_block = all_verts%start_block, all_verts%end_block
! !       DO level = 1, n_zlev
! !         DO neigbor=1,6
! ! 
! !           ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(1) = &
! !             & intp_2D_coeff%edge2vert_coeff_cc(:,vertex_block,neigbor)%x(1)
! !           ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(2) = &
! !             & intp_2D_coeff%edge2vert_coeff_cc(:,vertex_block,neigbor)%x(2)
! !           ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(3) = &
! !             & intp_2D_coeff%edge2vert_coeff_cc(:,vertex_block,neigbor)%x(3)
! ! 
! !         ENDDO ! neigbor=1,patch%cell_type
! !       ENDDO  !  level = 1, n_zlev
! !     ENDDO ! vertex_block
! !     !---------------------------------------------------------
! ! 
! !   END SUBROUTINE copy_2D_to_3D_intp_coeff
! !   !-------------------------------------------------------------------------


 
! !   !-------------------------------------------------------------------------
! !   !>
! !   !! Computes the coefficients that determine the scalar product on the primal grid. This
! !   !! scalar product depends on the grid geometry only and  is used to formulate the primitive
! !   !! equations in weak form. The coefficients are applied in module "mo_scalar_product".
! !   !! The following components of the data type "ocean_patch" are filled:
! !   !!   edge2cell_coeff  : coefficients for edge to cell mapping
! !   !!   edge2cell_coeff_t: coefficients for transposed of edge to cell mappings
! !   !!   edge2vert_coeff  : coefficients for edge to vertex mapping
! !   !!   edge2vert_coeff_t: coefficients for transposed of edge to vertex mappings
! !   !!   fixed_vol_norm   : summed volume weight of moved cell
! !   !!   variable_vol_norm: volume weight at the edges of moved cell
! !   !!
! !   !! @par Revision History
! !   !!  developed by Peter Korn, MPI-M  2010-09
! !   !!  Modification by Stephan Lorenz, 2010-11
! !   !!
! !   SUBROUTINE init_scalar_product_oce_3d( p_patch, p_coeff)
! ! 
! !     !  patch on which computation is performed
! !     !
! !     TYPE(t_patch), TARGET, INTENT(inout) :: p_patch
! !     TYPE(t_operator_coeff),INTENT(inout) :: p_coeff
! ! 
! !     CHARACTER(LEN=max_char_length), PARAMETER :: &
! !       & routine = ('mo_operator_ocean_coeff_3d:init_scalar_product_oce_3d')
! ! 
! !     INTEGER, PARAMETER :: no_cell_edges = 3
! !     INTEGER, PARAMETER :: no_vert_edges = 6
! !     INTEGER :: jb, je, jv, ie, ie_1, ie_2, icc,jk, jc
! !     INTEGER :: il_e,ib_e,k
! !     INTEGER :: il_c1, ib_c1, il_c2, ib_c2
! !     INTEGER :: il_v1, il_v2, ib_v1, ib_v2
! ! 
! !     INTEGER :: iil_c1(no_cell_edges), iil_c2(no_cell_edges)
! !     INTEGER :: iib_c1(no_cell_edges), iib_c2(no_cell_edges)
! ! 
! !     INTEGER :: jil_c1, jib_c1,jil_c2, jib_c2
! !     INTEGER :: rl_start_e, rl_end_e
! !     INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
! !     INTEGER :: rl_start, rl_end
! !     INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! ! 
! !     REAL(wp) :: z_tmp
! !     REAL(wp) :: norm_c1_c2, norm_v1_v2, norm
! !     REAL(wp) :: dual_edge_length(no_vert_edges)
! !     REAL(wp) :: vert_edge_dist(no_vert_edges,2)
! !     REAL(wp) :: vert_dual_mid_dist(no_vert_edges,2)
! !     REAL(wp) :: vert_edge_distance, vert_dual_mid_distance
! ! 
! !     TYPE(t_geographical_coordinates) :: gc_mid_dual_edge(no_vert_edges)
! !     TYPE(t_geographical_coordinates) :: gc1,gc2
! ! 
! !     TYPE(t_cartesian_coordinates)    :: cc_dual_edge(no_vert_edges), cc_edge(no_cell_edges)
! !     TYPE(t_cartesian_coordinates)    :: xx1,xx2
! !     TYPE(t_cartesian_coordinates)    :: vert1_midedge_cc(nproma,p_patch%nblks_v,no_vert_edges)
! !     TYPE(t_cartesian_coordinates)    :: vert2_midedge_cc(nproma,p_patch%nblks_v,no_vert_edges)
! !     TYPE(t_cartesian_coordinates)    :: cell2cell_cc
! !     TYPE(t_cartesian_coordinates)    :: cc_e0, cc_c1,cc_c2,cc_v0
! !     TYPE(t_cartesian_coordinates)    :: cv_c1_e0, cv_c2_e0, cv_c1_c2
! !     TYPE(t_cartesian_coordinates)    :: cc_mid_dual_edge(no_vert_edges)
! !     TYPE(t_cartesian_coordinates)    :: recon_vec_cc
! !     TYPE(t_cartesian_coordinates)    :: z_vec_c1(no_cell_edges),z_vec_c2(no_cell_edges)
! !     TYPE(t_cartesian_coordinates)    :: recon_vec_cc_v1(no_vert_edges)
! !     TYPE(t_cartesian_coordinates)    :: recon_vec_cc_v2(no_vert_edges)
! ! 
! !     REAL(wp) :: z_edge_length(no_cell_edges)
! !     REAL(wp) :: z_cell_edge_dist_c1(no_cell_edges,2),z_cell_edge_dist_c2(no_cell_edges,2)
! !     REAL(wp) :: z_y
! ! 
! !     REAL(wp) :: z_sync_c(nproma,n_zlev, p_patch%nblks_c)
! !     REAL(wp) :: z_sync_e(nproma,n_zlev, p_patch%nblks_e)
! !     REAL(wp) :: z_sync_v(nproma,n_zlev, p_patch%nblks_v)
! ! 
! !     REAL(wp) :: omega
! ! 
! !     LOGICAL, PARAMETER :: mid_point_dual_edge = .TRUE. !Please do not change this unless
! !     !you are sure, you know what you do.
! !     LOGICAL, PARAMETER :: larc_length = .FALSE.
! !     !-----------------------------------------------------------------------
! !     CALL message (TRIM(routine), 'start')
! ! 
! !     omega = grid_angular_velocity
! ! 
! !     rl_start     = 1
! !     rl_end       = min_rlcell_int
! !     i_startblk   = p_patch%cells%start_blk(rl_start,1)
! !     i_endblk     = p_patch%cells%end_blk(rl_end,1)
! ! 
! !     rl_start_e   = 1
! !     rl_end_e     = min_rledge ! Loop over the whole local domain
! !     i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
! !     i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
! ! 
! ! 
! !     p_coeff%fixed_vol_norm(:,:,:) = 0._wp
! ! 
! !     !-----------------------------------------------------------------------
! !     !STEP 1: edge2cell and cell2edge coefficients
! !     DO jk=1,n_zlev
! !       edge_blk_loop_primal: DO jb = i_startblk_e, i_endblk_e
! ! 
! !         CALL get_indices_e(p_patch, jb,&
! !           & i_startblk_e, i_endblk_e,&
! !           & i_startidx_e, i_endidx_e,&
! !           & rl_start_e, rl_end_e)
! ! 
! !         edge_idx_loop_primal: DO je =  i_startidx_e, i_endidx_e
! ! 
! !           !Get indices of two adjacent triangles
! !           il_c1 = p_patch%edges%cell_idx(je,jb,1)
! !           ib_c1 = p_patch%edges%cell_blk(je,jb,1)
! !           il_c2 = p_patch%edges%cell_idx(je,jb,2)
! !           ib_c2 = p_patch%edges%cell_blk(je,jb,2)
! ! 
! !           ! Go only over edges where at least one neighboring triangle is not in the halo
! !           IF(il_c1<=0 .OR. il_c2<=0) CYCLE
! !           IF(.NOT.p_patch%cells%owner_mask(il_c1,ib_c1) .AND. &
! !             & .NOT.p_patch%cells%owner_mask(il_c2,ib_c2) ) CYCLE
! ! 
! !           cc_c1 = gc2cc(p_patch%cells%center(il_c1, ib_c1))
! !           cc_c2 = gc2cc(p_patch%cells%center(il_c2, ib_c2))
! ! 
! !           z_cell_edge_dist_c1 = 0.0_wp
! !           z_cell_edge_dist_c2 = 0.0_wp
! ! 
! !           !normals in cell 1
! !           DO ie = 1, no_cell_edges
! ! 
! !             !actual edges of cell c1
! !             iil_c1(ie)  = p_patch%cells%edge_idx(il_c1,ib_c1,ie)
! !             iib_c1(ie)  = p_patch%cells%edge_blk(il_c1,ib_c1,ie)
! ! 
! !             cc_edge(ie) = gc2cc(p_patch%edges%center(iil_c1(ie),iib_c1(ie)))
! ! 
! !             !calculate edge length
! !             !get vertex indices adjacent to actual edge
! !             il_v1       = p_patch%edges%vertex_idx(iil_c1(ie),iib_c1(ie),1)
! !             ib_v1       = p_patch%edges%vertex_blk(iil_c1(ie),iib_c1(ie),1)
! !             il_v2       = p_patch%edges%vertex_idx(iil_c1(ie),iib_c1(ie),2)
! !             ib_v2       = p_patch%edges%vertex_blk(iil_c1(ie),iib_c1(ie),2)
! ! 
! !             !get vertex positions
! !             xx1         = gc2cc(p_patch%verts%vertex(il_v1,ib_v1))
! !             xx2         = gc2cc(p_patch%verts%vertex(il_v2,ib_v2))
! ! 
! !             IF(larc_length)THEN
! !               norm              = SQRT(SUM(xx1%x*xx1%x))
! !               xx1%x             = xx1%x/norm
! !               norm              = SQRT(SUM(xx2%x*xx2%x))
! !               xx2%x             = xx2%x/norm
! !               z_edge_length(ie) = arc_length(xx2,xx1)
! !             ELSE
! !               z_edge_length(ie) = SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
! !             ENDIF
! ! 
! !             !calculate cell-edge distance as half of cell-cell distance
! !             !get cell indices adjacent to actual edge
! !             jil_c1 = p_patch%edges%cell_idx(iil_c1(ie),iib_c1(ie),1)
! !             jib_c1 = p_patch%edges%cell_blk(iil_c1(ie),iib_c1(ie),1)
! !             jil_c2 = p_patch%edges%cell_idx(iil_c1(ie),iib_c1(ie),2)
! !             jib_c2 = p_patch%edges%cell_blk(iil_c1(ie),iib_c1(ie),2)
! ! 
! !             !get cell positions
! !             xx1    = gc2cc(p_patch%cells%center(jil_c1,jib_c1))
! !             xx2    = gc2cc(p_patch%cells%center(jil_c2,jib_c2))
! ! 
! !             IF(jil_c1==il_c1.AND.jib_c1==ib_c1)THEN
! !               k=1
! !             ELSEIF(jil_c2==il_c1.AND.jib_c2==ib_c1)THEN
! !               k=2
! !             ENDIF
! ! 
! !             IF(larc_length)THEN
! !               norm                      = SQRT(SUM(xx1%x*xx1%x))
! !               xx1%x                     = xx1%x/norm
! !               norm                      = SQRT(SUM(xx2%x*xx2%x))
! !               xx2%x                     = xx2%x/norm
! !               norm                      = SQRT(SUM(cc_edge(ie)%x*cc_edge(ie)%x))
! !               cc_edge(ie)%x             = cc_edge(ie)%x/norm
! !               z_cell_edge_dist_c1(ie,1) = arc_length(cc_edge(ie),xx1)
! !               z_cell_edge_dist_c1(ie,2) = arc_length(cc_edge(ie),xx2)
! !             ELSE
! !               z_cell_edge_dist_c1(ie,1) = SQRT(SUM((cc_edge(ie)%x-xx1%x)*(cc_edge(ie)%x-xx1%x)))
! !               z_cell_edge_dist_c1(ie,2) = SQRT(SUM((cc_edge(ie)%x-xx2%x)*(cc_edge(ie)%x-xx2%x)))
! !             ENDIF
! !             p_coeff%dist_cell2edge(iil_c1(ie),jk,iib_c1(ie),1) = z_cell_edge_dist_c1(ie,1)
! !             p_coeff%dist_cell2edge(iil_c1(ie),jk,iib_c1(ie),2) = z_cell_edge_dist_c1(ie,2)
! ! 
! !             z_vec_c1(ie)%x = cc_edge(ie)%x - cc_c1%x     !p_patch%edges%primal_cart_normal(iil_c1(ie),iib_c1(ie))
! !             norm           = SQRT(SUM( z_vec_c1(ie)%x* z_vec_c1(ie)%x))
! !             !write(*,*)'NORM:',norm !TODOram
! ! 
! !             p_coeff%edge2cell_coeff_cc(il_c1,jk,ib_c1,ie)%x = &
! !               & z_vec_c1(ie)%x*p_patch%cells%edge_orientation(il_c1,ib_c1,ie)*z_edge_length(ie)
! ! 
! !             p_coeff%fixed_vol_norm(il_c1,jk,ib_c1) = p_coeff%fixed_vol_norm(il_c1,jk,ib_c1)&
! !               & +  0.5_wp*norm*z_edge_length(ie)
! !             p_coeff%variable_vol_norm(il_c1,jk,ib_c1,ie) = 0.5_wp*norm*z_edge_length(ie)
! ! 
! !             !write(*,*)'edge length   :',z_edge_length(ie),p_patch%edges%primal_edge_length(iil_c1(ie),iib_c1(ie))/re
! !             !write(*,*)'cell-edge dist:', z_cell_edge_dist_c1(ie,k),p_patch%edges%edge_cell_length(iil_c1(ie),iib_c1(ie),k)/re
! !           END DO
! ! 
! !           !normals in cell 2
! !           DO ie = 1, no_cell_edges
! ! 
! !             !actual edges of cell c2
! !             iil_c2(ie)  = p_patch%cells%edge_idx(il_c2,ib_c2,ie)
! !             iib_c2(ie)  = p_patch%cells%edge_blk(il_c2,ib_c2,ie)
! ! 
! ! 
! !             cc_edge(ie) = gc2cc(p_patch%edges%center(iil_c2(ie),iib_c2(ie)))
! ! 
! !             !calculate edge length
! !             !get vertex indices adjacent to actual edge
! !             il_v1       = p_patch%edges%vertex_idx(iil_c2(ie),iib_c2(ie),1)
! !             ib_v1       = p_patch%edges%vertex_blk(iil_c2(ie),iib_c2(ie),1)
! !             il_v2       = p_patch%edges%vertex_idx(iil_c2(ie),iib_c2(ie),2)
! !             ib_v2       = p_patch%edges%vertex_blk(iil_c2(ie),iib_c2(ie),2)
! ! 
! !             !get vertex positions
! !             xx1         = gc2cc(p_patch%verts%vertex(il_v1,ib_v1))
! !             xx2         = gc2cc(p_patch%verts%vertex(il_v2,ib_v2))
! ! 
! !             IF(larc_length)THEN
! !               norm              = SQRT(SUM(xx1%x*xx1%x))
! !               xx1%x             = xx1%x/norm
! ! 
! !               norm              = SQRT(SUM(xx2%x*xx2%x))
! !               xx2%x             = xx2%x/norm
! ! 
! !               z_edge_length(ie) = arc_length(xx2,xx1)
! !               !z_edge_length(ie) = p_patch%edges%primal_edge_length(iil_c2(ie),iib_c2(ie))/re
! !               !write(*,*)'arc length',arc_length(xx2,xx1),z_edge_length(ie),SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
! !             ELSE
! !               z_edge_length(ie) = SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
! !             ENDIF
! !             !calculate cell-edge distance as half of cell-cell distance
! !             !get cell indices adjacent to actual edge
! !             jil_c1 = p_patch%edges%cell_idx(iil_c2(ie),iib_c2(ie),1)
! !             jib_c1 = p_patch%edges%cell_blk(iil_c2(ie),iib_c2(ie),1)
! !             jil_c2 = p_patch%edges%cell_idx(iil_c2(ie),iib_c2(ie),2)
! !             jib_c2 = p_patch%edges%cell_blk(iil_c2(ie),iib_c2(ie),2)
! ! 
! !             IF (jil_c2 < 0) THEN
! !               WRITE(0,*)'p_patch%edges%cell_idx:',p_patch%edges%cell_idx
! !             ENDIF
! !             !get cell positions
! !             xx1 = gc2cc(p_patch%cells%center(jil_c1,jib_c1))
! !             xx2 = gc2cc(p_patch%cells%center(jil_c2,jib_c2))
! ! 
! !             IF(jil_c1==il_c2.AND.jib_c1==ib_c2)THEN
! !               k=1
! !             ELSEIF(jil_c2==il_c2.AND.jib_c2==ib_c2)THEN
! !               k=2
! !             ENDIF
! ! 
! !             IF(larc_length)THEN
! !               norm                      = SQRT(SUM(xx1%x*xx1%x))
! !               xx1%x                     = xx1%x/norm
! !               norm                      = SQRT(SUM(xx2%x*xx2%x))
! !               xx2%x                     = xx2%x/norm
! !               norm                      = SQRT(SUM(cc_edge(ie)%x*cc_edge(ie)%x))
! !               cc_edge(ie)%x             = cc_edge(ie)%x/norm
! !               z_cell_edge_dist_c2(ie,1) = arc_length(cc_edge(ie),xx1)
! !               z_cell_edge_dist_c2(ie,2) = arc_length(cc_edge(ie),xx2)
! !             ELSE
! !               z_cell_edge_dist_c2(ie,1) = SQRT(SUM((cc_edge(ie)%x-xx1%x)*(cc_edge(ie)%x-xx1%x)))
! !               z_cell_edge_dist_c2(ie,2) = SQRT(SUM((cc_edge(ie)%x-xx2%x)*(cc_edge(ie)%x-xx2%x)))
! !             ENDIF
! !             p_coeff%dist_cell2edge(iil_c2(ie),jk,iib_c2(ie),1) = z_cell_edge_dist_c2(ie,1)
! !             p_coeff%dist_cell2edge(iil_c2(ie),jk,iib_c2(ie),2) = z_cell_edge_dist_c2(ie,2)
! ! 
! !             z_vec_c2(ie)%x = cc_edge(ie)%x - cc_c2%x  !p_patch%edges%primal_cart_normal(iil_c2(ie),iib_c2(ie))
! !             norm           = SQRT(SUM( z_vec_c2(ie)%x* z_vec_c2(ie)%x))
! ! 
! !             p_coeff%edge2cell_coeff_cc(il_c2,jk,ib_c2,ie)%x = &
! !               &  z_vec_c2(ie)%x * p_patch%cells%edge_orientation(il_c2,ib_c2,ie) * &
! !               &  z_edge_length(ie)
! ! 
! !             p_coeff%fixed_vol_norm(il_c2,jk,ib_c2) = p_coeff%fixed_vol_norm(il_c2,jk,ib_c2)&
! !               & + 0.5_wp*norm*z_edge_length(ie)
! !             p_coeff%variable_vol_norm(il_c2,jk,ib_c2,ie) = 0.5_wp*norm*z_edge_length(ie)
! ! 
! !           END DO
! !         END DO edge_idx_loop_primal
! !       END DO edge_blk_loop_primal
! !     END DO
! !     !In the edge loop above each triangle is visisted three times. Since the "fixed_vol_norm" is
! !     !accumulated we correct its value here:
! !     p_coeff%fixed_vol_norm = p_coeff%fixed_vol_norm/3.0_wp
! ! 
! ! 
! !     !Assign values to dynamical coefficients for surface layer
! !     edge_blk_loop_dyn: DO jb = i_startblk, i_endblk
! ! 
! !       CALL get_indices_c(p_patch, jb,&
! !         & i_startblk, i_endblk,&
! !         & i_startidx, i_endidx,&
! !         & rl_start, rl_end)
! ! 
! !       edge_idx_loop_dyn: DO jc =  i_startidx, i_endidx
! !         DO ie = 1, no_cell_edges
! !           p_coeff%edge2cell_coeff_cc_dyn(jc,1,jb,ie)%x=p_coeff%edge2cell_coeff_cc(jc,1,jb,ie)%x
! !         END DO
! !       END DO edge_idx_loop_dyn
! !     END DO edge_blk_loop_dyn
! ! 
! ! 
! !     !commented put for testing--------------------------------
! ! 
! !     !    !merge fixed volume and edge2cell coeff
! !     !     DO jk=1,n_zlev
! !     !       DO jb = i_startblk, i_endblk
! !     !
! !     !       CALL get_indices_c(p_patch, jb,&
! !     !                        & i_startblk, i_endblk,&
! !     !                        & i_startidx, i_endidx,&
! !     !                        & rl_start, rl_end)
! !     !
! !     !         DO jc =  i_startidx, i_endidx
! !     !           DO ie=1,no_cell_edges
! !     !             p_coeff%edge2cell_coeff_cc(jc,jk,jb,ie)%x&
! !     !             & = p_coeff%edge2cell_coeff_cc(jc,jk,jb,ie)%x/p_coeff%fixed_vol_norm(jc,jk,jb)
! !     !           END DO
! !     !         END DO
! !     !       END DO
! !     !     END DO
! ! 
! !     rl_start   = 1
! !     rl_end     = min_rledge_int
! !     i_startblk = p_patch%edges%start_blk(rl_start,1)
! !     i_endblk   = p_patch%edges%end_blk(rl_end,1)
! ! 
! !     DO jk=1,n_zlev
! !       edge_blk_loop_secondary: DO jb = i_startblk, i_endblk
! !         CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,&
! !           & i_startidx, i_endidx, rl_start, rl_end)
! !         edge_idx_loop_secondary: DO je =  i_startidx, i_endidx
! ! 
! !           !Get indices of two adjacent triangles
! !           il_c1      = p_patch%edges%cell_idx(je,jb,1)
! !           ib_c1      = p_patch%edges%cell_blk(je,jb,1)
! !           il_c2      = p_patch%edges%cell_idx(je,jb,2)
! !           ib_c2      = p_patch%edges%cell_blk(je,jb,2)
! ! 
! !           !cartesian coordinates of edge and neighbor cells on 1-sphere
! !           cc_e0      = gc2cc(p_patch%edges%center(je,jb))
! !           cc_c1      = gc2cc(p_patch%cells%center(il_c1,ib_c1))
! !           cc_c2      = gc2cc(p_patch%cells%center(il_c2,ib_c2))
! ! 
! !           !cartesian vectors from:
! !           !cell 2 to cell 1, cell 1 to edge je and cell 2 to edge je
! !           cv_c1_c2%x = cc_c1%x - cc_c2%x
! !           cv_c1_e0%x = cc_e0%x - cc_c1%x
! !           cv_c2_e0%x = cc_e0%x - cc_c2%x
! ! 
! !           IF(larc_length)THEN
! !             norm        = SQRT(SUM(cc_e0%x*cc_e0%x))
! !             cc_e0%x     = cc_e0%x/norm
! !             norm        = SQRT(SUM(cc_c1%x*cc_c1%x))
! !             cc_c1%x     = cc_c1%x/norm
! !             norm        = SQRT(SUM(cc_c2%x*cc_c2%x))
! !             cc_c2%x     = cc_c2%x/norm
! !             norm_c1_c2  = arc_length(cc_c1, cc_c2)
! !           ELSE
! !             norm_c1_c2  = SQRT(SUM(cv_c1_e0%x*cv_c1_e0%x))+SQRT(SUM(cv_c2_e0%x*cv_c2_e0%x)) !SQRT(SUM(cv_c1_c2%x*cv_c1_c2%x))!
! !           ENDIF
! ! 
! !           !Determine which edge of both of the two adjacent cells corresponds to the
! !           !actual edge "je". This information is used below for the edge-orientation.
! !           DO ie = 1, no_cell_edges
! !             IF (p_patch%cells%edge_idx(il_c1,ib_c1,ie) == je.AND.&
! !               & p_patch%cells%edge_blk(il_c1,ib_c1,ie) == jb) THEN
! !               ie_1 = ie
! !             END IF
! !             IF (p_patch%cells%edge_idx(il_c2,ib_c2,ie) == je.AND.&
! !               & p_patch%cells%edge_blk(il_c2,ib_c2,ie) == jb) THEN
! !               ie_2 = ie
! !             END IF
! !           END DO
! ! 
! !           p_coeff%edge2cell_coeff_cc_t(je,jk,jb,1)%x&
! !             & = cv_c1_e0%x * p_patch%cells%edge_orientation(il_c1,ib_c1,ie_1)/norm_c1_c2
! ! 
! !           p_coeff%edge2cell_coeff_cc_t(je,jk,jb,2)%x&
! !             & = cv_c2_e0%x * p_patch%cells%edge_orientation(il_c2,ib_c2,ie_2)/norm_c1_c2
! ! 
! !         END DO edge_idx_loop_secondary
! !       END DO edge_blk_loop_secondary
! !     END DO
! !     !------------------------------------------------------------------------------
! !     !STEP 2: edge2vert coefficients for dual grid
! !     !------------------------------------------------------------------------------
! ! 
! !     rl_start = 1
! !     rl_end   = min_rlvert ! Loop over the whole local domain
! ! 
! !     i_startblk = p_patch%verts%start_blk(rl_start,1)
! !     i_endblk   = p_patch%verts%end_blk(rl_end,1)
! ! 
! !     DO jk=1,n_zlev
! !       vert_blk_loop: DO jb = i_startblk, i_endblk
! !         CALL get_indices_v(p_patch, jb, i_startblk, i_endblk,&
! !           & i_startidx, i_endidx, rl_start, rl_end)
! !         vert_idx_loop: DO jv =  i_startidx, i_endidx
! !           ! current number of edges around vertex (5 or 6)
! !           cc_v0        = gc2cc(p_patch%verts%vertex(jv,jb))
! ! 
! !           ! Go only over vertices where all edges are in the local domain (but maybe in the halo)
! !           IF(ANY(p_patch%verts%edge_idx(jv,jb,:)<0)) CYCLE
! ! 
! !           DO ie = 1, no_vert_edges  ! #slo# it_vertedges ??
! ! 
! !             il_e             = p_patch%verts%edge_idx(jv,jb,ie)
! !             ib_e             = p_patch%verts%edge_blk(jv,jb,ie)
! ! 
! !             ! #slo# - I assume all geographical coordinates are already synchronized
! ! 
! !             cc_dual_edge(ie) = gc2cc(p_patch%edges%center(il_e,ib_e))
! !             !Parts of this code parallels the implementation in the grid-generator
! !             !module "mo_geometry".
! !             !
! !             !1) determine normal vector from adjacent cell to adjacent cell
! !             !   in cartesian coordinate for moved dual cell
! !             !Get indices of two adjacent triangles
! !             il_c1            = p_patch%edges%cell_idx(il_e,ib_e,1)
! !             ib_c1            = p_patch%edges%cell_blk(il_e,ib_e,1)
! !             il_c2            = p_patch%edges%cell_idx(il_e,ib_e,2)
! !             ib_c2            = p_patch%edges%cell_blk(il_e,ib_e,2)
! ! 
! !             xx1              = gc2cc(p_patch%cells%center(il_c1,ib_c1))
! !             norm             = SQRT(SUM(xx1%x*xx1%x))
! !             xx1%x            = xx1%x/norm
! ! 
! !             xx2              = gc2cc(p_patch%cells%center(il_c2,ib_c2))
! !             norm             = SQRT(SUM(xx2%x*xx2%x))
! !             xx2%x            = xx2%x/norm
! ! 
! !             cell2cell_cc%x   = xx2%x - xx1%x
! !             IF(larc_length)THEN
! !               norm_c1_c2 = arc_length(xx1,xx2)
! !             ELSE
! !               norm_c1_c2 = SQRT(SUM(cell2cell_cc%x*cell2cell_cc%x))
! !             ENDIF
! !             dual_edge_length(ie) = norm_c1_c2
! !             !          cell2cell_cc%x       = cell2cell_cc%x/norm_c1_c2
! ! 
! !             IF(mid_point_dual_edge)THEN
! !               cc_mid_dual_edge(ie)%x = 0.5_wp*(xx2%x+xx1%x)
! !               gc_mid_dual_edge(ie)   = cc2gc(cc_mid_dual_edge(ie))
! ! 
! !               IF(coriolis_type==full_coriolis)THEN
! !                 p_patch%edges%f_e(il_e, ib_e) = 2._wp*omega*SIN(gc_mid_dual_edge(ie)%lat)
! !               ELSEIF(coriolis_type==beta_plane_coriolis)THEN
! !                 gc1%lat = basin_center_lat* deg2rad - 0.5_wp*basin_height_deg*deg2rad
! !                 gc1%lon = 0.0_wp
! !                 xx1     = gc2cc(gc1)
! ! 
! !                 gc2%lat = gc_mid_dual_edge(ie)%lat!*deg2rad
! !                 gc2%lon = 0.0_wp
! !                 xx2     = gc2cc(gc2)
! !                 z_y     = earth_radius*arc_length(xx2,xx1)
! ! 
! !                 !z_y = p_patch%edges%center(je,jb)%lat - z_lat_basin_center
! !                 p_patch%edges%f_e(il_e, ib_e) = &
! !                   & 2.0_wp*omega*( SIN(basin_center_lat * deg2rad) + &
! !                   & (COS(basin_center_lat * deg2rad)/earth_radius)*z_y)
! !               ENDIF
! !             ELSE
! !               cc_mid_dual_edge(ie)%x = cc_dual_edge(ie)%x
! !               gc_mid_dual_edge(ie)   = cc2gc(cc_mid_dual_edge(ie))
! !             ENDIF
! ! 
! !             !2) determine vector from adjacent vertex to adjacent vertex
! !             !   in cartesian coordinate for moved dual cell
! !             !Get indices of two adjacent vertices
! !             il_v1 = p_patch%edges%vertex_idx(il_e,ib_e,1)
! !             ib_v1 = p_patch%edges%vertex_blk(il_e,ib_e,1)
! !             il_v2 = p_patch%edges%vertex_idx(il_e,ib_e,2)
! !             ib_v2 = p_patch%edges%vertex_blk(il_e,ib_e,2)
! ! 
! !             xx1   = gc2cc(p_patch%verts%vertex(il_v1,ib_v1))
! !             norm  = SQRT(SUM(xx1%x*xx1%x))
! !             xx1%x = xx1%x/norm
! ! 
! !             xx2   = gc2cc(p_patch%verts%vertex(il_v2,ib_v2))
! !             norm  = SQRT(SUM(xx2%x*xx2%x))
! !             xx2%x = xx2%x/norm
! ! 
! !             vert1_midedge_cc(jv, jb, ie)%x = cc_mid_dual_edge(ie)%x - xx1%x
! !             vert2_midedge_cc(jv, jb, ie)%x = cc_mid_dual_edge(ie)%x - xx2%x
! !             norm = SQRT(SUM(vert1_midedge_cc(jv,jb,ie)%x*vert1_midedge_cc(jv,jb,ie)%x))
! !             vert1_midedge_cc(jv, jb, ie)%x = vert1_midedge_cc(jv, jb, ie)%x/norm
! ! 
! !             norm = SQRT(SUM(vert2_midedge_cc(jv,jb,ie)%x*vert2_midedge_cc(jv,jb,ie)%x))
! !             vert2_midedge_cc(jv, jb, ie)%x = vert2_midedge_cc(jv, jb, ie)%x/norm
! ! 
! ! 
! !             !calculate vertex edge distance
! !             IF(larc_length)THEN
! !               vert_edge_dist(ie,1)     = arc_length (cc_dual_edge(ie), xx1)
! !               vert_edge_dist(ie,2)     = arc_length (cc_dual_edge(ie), xx2)
! !               vert_dual_mid_dist(ie,1) = arc_length (cc_mid_dual_edge(ie), xx1)
! !               vert_dual_mid_dist(ie,2) = arc_length (cc_mid_dual_edge(ie), xx2)
! !             ELSE
! !               vert_edge_dist(ie,1)&
! !                 & = SQRT(SUM((cc_dual_edge(ie)%x - xx1%x)*(cc_dual_edge(ie)%x - xx1%x)))
! !               vert_edge_dist(ie,2)&
! !                 & = SQRT(SUM((cc_dual_edge(ie)%x - xx2%x)*(cc_dual_edge(ie)%x - xx2%x)))
! !               vert_dual_mid_dist(ie,1)&
! !                 & = SQRT(SUM((cc_mid_dual_edge(ie)%x - xx1%x)*(cc_mid_dual_edge(ie)%x - xx1%x)))
! !               vert_dual_mid_dist(ie,2)&
! !                 & = SQRT(SUM((cc_mid_dual_edge(ie)%x - xx2%x)*(cc_mid_dual_edge(ie)%x - xx2%x)))
! !             ENDIF
! ! 
! !             !calculate normal vector that is perpendicular to vertex-vertex- and edge position vector
! !             !If one uses the edge position vector this results in the moved primal normal.
! !             recon_vec_cc_v1(ie)   = vector_product(vert1_midedge_cc(jv, jb, ie),&
! !               & cc_mid_dual_edge(ie))
! !             norm                  = SQRT(SUM(recon_vec_cc_v1(ie)%x*recon_vec_cc_v1(ie)%x))
! !             recon_vec_cc_v1(ie)%x = recon_vec_cc_v1(ie)%x/norm
! ! 
! !             recon_vec_cc_v2(ie)   = vector_product(vert2_midedge_cc(jv,jb,ie),&
! !               & cc_mid_dual_edge(ie))
! !             norm                  = SQRT(SUM(recon_vec_cc_v2(ie)%x*recon_vec_cc_v2(ie)%x))
! !             recon_vec_cc_v2(ie)%x = recon_vec_cc_v2(ie)%x/norm
! ! 
! !             !Fix orientation
! !             z_tmp = DOT_PRODUCT(recon_vec_cc_v1(ie)%x,&
! !               & p_patch%edges%primal_cart_normal(il_e,ib_e)%x)
! !             IF (z_tmp <0._wp) recon_vec_cc_v1(ie)%x = -1._wp * recon_vec_cc_v1(ie)%x
! ! 
! !             z_tmp = DOT_PRODUCT(recon_vec_cc_v2(ie)%x,&
! !               & p_patch%edges%primal_cart_normal(il_e,ib_e)%x)
! !             IF (z_tmp <0._wp) recon_vec_cc_v2(ie)%x = -1._wp * recon_vec_cc_v2(ie)%x
! ! 
! ! 
! !             IF ( (p_patch%edges%vertex_idx(il_e,ib_e,1) == jv) .AND. &
! !               & (p_patch%edges%vertex_blk(il_e,ib_e,1) == jb)         ) THEN
! ! 
! !               vert_edge_distance     = vert_edge_dist(ie,1)
! !               vert_dual_mid_distance = vert_dual_mid_dist(ie,1)
! !               recon_vec_cc           = recon_vec_cc_v1(ie)
! !               !PK: not used
! !               !              p_coeff%edge2vert_vector_cc(jv,jb,ie)=&
! !               !              &vector_product(vert1_midedge_cc(jv, jb, ie), cc_mid_dual_edge(ie))
! !               p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x=vert1_midedge_cc(jv, jb, ie)%x
! ! 
! !               z_tmp = DOT_PRODUCT(p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x,&
! !                 & p_patch%edges%primal_cart_normal(il_e,ib_e)%x)
! ! 
! !               IF (z_tmp <0._wp) p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x&
! !                 & = -1._wp * p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x
! ! 
! ! 
! !               p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x = vert_edge_distance&
! !                 & *p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x
! ! 
! !             ELSE IF ( (p_patch%edges%vertex_idx(il_e,ib_e,2) == jv) .AND. &
! !               & (p_patch%edges%vertex_blk(il_e,ib_e,2) == jb) ) THEN
! ! 
! !               vert_edge_distance     = vert_edge_dist(ie,2)
! !               vert_dual_mid_distance = vert_dual_mid_dist(ie,2)
! !               recon_vec_cc           = recon_vec_cc_v2(ie)
! ! 
! !               !              p_coeff%edge2vert_vector_cc(jv,jb,ie)=&
! !               !              &vector_product(vert2_midedge_cc(jv, jb, ie), cc_mid_dual_edge(ie))
! ! 
! !               p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x=vert2_midedge_cc(jv,jb,ie)%x
! ! 
! !               z_tmp = DOT_PRODUCT(p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x,&
! !                 & p_patch%edges%primal_cart_normal(il_e,ib_e)%x)
! ! 
! !               IF (z_tmp <0._wp) p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x =&
! !                 & -1._wp * p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x
! ! 
! !               p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x = vert_edge_distance&
! !                 & *p_coeff%edge2vert_vector_cc(jv,jk,jb,ie)%x
! ! 
! ! 
! !             ELSE
! !               CALL message (TRIM(routine), 'WARNING - vert_edge_distance not found')
! !               WRITE(0,'(a,7i5)') &
! !                 & 'jv, jb, edge, p_patch%edges%vertex_idx/blk(il_e,ib_e,1-2)=', &
! !                 & jv, jb, ie, &
! !                 & p_patch%edges%vertex_idx(il_e,ib_e,1), &
! !                 & p_patch%edges%vertex_blk(il_e,ib_e,1), &
! !                 & p_patch%edges%vertex_idx(il_e,ib_e,2), &
! !                 & p_patch%edges%vertex_blk(il_e,ib_e,2)
! ! 
! !             END IF
! ! 
! !             p_coeff%variable_dual_vol_norm(jv,jk,jb,ie) = &
! !               & 0.5_wp*dual_edge_length(ie)*vert_dual_mid_distance
! !             !vert_edge_distance*dual_edge_length(ie)!
! ! 
! !             p_coeff%edge2vert_coeff_cc(jv,jk,jb,ie)%x   = &
! !               & recon_vec_cc%x*dual_edge_length(ie)*vert_dual_mid_distance
! ! 
! !             norm_v1_v2 = &
! !               & SQRT(SUM(vert1_midedge_cc(jv, jb, ie)%x*vert1_midedge_cc(jv, jb, ie)%x)) + &
! !               & SQRT(SUM(vert2_midedge_cc(jv, jb, ie)%x*vert2_midedge_cc(jv, jb, ie)%x))
! ! 
! !             p_coeff%edge2vert_coeff_cc_t(il_e,jk,ib_e,1)%x = vert1_midedge_cc(jv, jb, ie)%x * &
! !               & ( p_patch%edges%system_orientation(il_e,ib_e)/norm_v1_v2 )
! ! 
! !             p_coeff%edge2vert_coeff_cc_t(il_e,jk,ib_e,2)%x = vert2_midedge_cc(jv, jb, ie)%x * &
! !               & ( p_patch%edges%system_orientation(il_e,ib_e)/norm_v1_v2 )
! ! 
! !           END DO
! !         ENDDO vert_idx_loop
! !       END DO vert_blk_loop
! !     END DO
! ! 
! !     !     !--------------------------------------------------------------------------
! !     !     ! SYNCHRONIZE ALL ELEMENTS OF V_BASE:
! !     !     ! synchronize elements on cells
! !     !     DO ie = 1, no_cell_edges
! !     !       DO icc = 1, 3
! !     !         z_sync_c(:,:,:) =  p_coeff%edge2cell_coeff_cc(:,:,:,ie)%x(icc)
! !     !         CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:,:))
! !     !         p_coeff%edge2cell_coeff_cc(:,:,:,ie)%x(icc) = z_sync_c(:,:,:)
! !     !       END DO
! !     !       z_sync_c(:,:,:) = p_coeff%variable_vol_norm(:,:,:,ie)
! !     !       CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:,:))
! !     !       p_coeff%variable_vol_norm(:,:,:,ie) = z_sync_c(:,:,:)
! !     !     END DO
! !     !     CALL sync_patch_array(SYNC_C, p_patch,p_coeff%fixed_vol_norm)
! !     !
! !     !     ! synchronize elements on edges
! !     !     DO ie = 1, 2
! !     !       DO icc = 1, 3
! !     !         z_sync_e(:,:,:) =  p_coeff%edge2vert_coeff_cc_t(:,:,:,ie)%x(icc)
! !     !         CALL sync_patch_array(SYNC_E, p_patch, z_sync_e(:,:,:))
! !     !         p_coeff%edge2vert_coeff_cc_t(:,:,:,ie)%x(icc) = z_sync_e(:,:,:)
! !     !
! !     !         z_sync_e(:,:,:) =  p_coeff%edge2cell_coeff_cc_t(:,:,:,ie)%x(icc)
! !     !         CALL sync_patch_array(SYNC_E, p_patch, z_sync_e(:,:,:))
! !     !         p_coeff%edge2cell_coeff_cc_t(:,:,:,ie)%x(icc) = z_sync_e(:,:,:)
! !     !       END DO
! !     !     END DO
! !     !
! !     !     ! synchronize cartesian coordinates on vertices:
! !     !     DO ie = 1, no_vert_edges
! !     !       DO icc = 1, 3
! !     !         z_sync_v(:,:,:) =  p_coeff%edge2vert_vector_cc(:,:,:,ie)%x(icc)
! !     !         CALL sync_patch_array(SYNC_V, p_patch, z_sync_v(:,:,:))
! !     !         p_coeff%edge2vert_vector_cc(:,:,:,ie)%x(icc) = z_sync_v(:,:,:)
! !     !
! !     !         z_sync_v(:,:,:) = p_coeff%edge2vert_coeff_cc(:,:,:,ie)%x(icc)
! !     !         CALL sync_patch_array(SYNC_V, p_patch, z_sync_v(:,:,:))
! !     !         p_coeff%edge2vert_coeff_cc(:,:,:,ie)%x(icc) = z_sync_v(:,:,:)
! !     !       END DO
! !     !       z_sync_v(:,:,:) = p_coeff%variable_dual_vol_norm(:,:,:,ie)
! !     !       CALL sync_patch_array(SYNC_V, p_patch, z_sync_v(:,:,:))
! !     !       p_coeff%variable_dual_vol_norm(:,:,:,ie) = z_sync_v(:,:,:)
! !     !     END DO
! ! 
! !     CALL message (TRIM(routine), 'end')
! ! 
! !   END SUBROUTINE init_scalar_product_oce_3d
! !   !-------------------------------------------------------------------------


  

END MODULE mo_operator_ocean_coeff_3d

