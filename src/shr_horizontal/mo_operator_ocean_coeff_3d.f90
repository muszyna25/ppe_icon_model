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
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert, max_dom,success,&
    &                               max_char_length, beta_plane_coriolis,full_coriolis, &
    &                               min_rledge_int,min_rlcell_int,min_rlvert_int,&
    &                               SEA_BOUNDARY, BOUNDARY, SEA, min_dolic
  USE mo_math_constants,      ONLY: deg2rad, pi, rad2deg
  USE mo_physical_constants,  ONLY: earth_radius
  USE mo_math_utilities,      ONLY: gc2cc, cc2gc, t_cartesian_coordinates,      &
    &                               t_geographical_coordinates, vector_product, &
    &                               arc_length
  USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m, no_tracer, t_ref, s_ref, &
    &                               coriolis_type, basin_center_lat, basin_height_deg
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
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
  PUBLIC  :: init_operator_coeffs

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
  INTEGER,PARAMETER :: no_primal_edges = 3

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
    REAL(wp), ALLOCATABLE                      :: edge2edge_viacell_coeff(:,:,:,:)

    !coefficient for surface layer, changes in time, in contrast to other coefficients
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2cell_coeff_cc_dyn(:,:,:,:)
    !TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc_dyn(:,:,:,:)

    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc(:,:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc_t(:,:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_vector_cc(:,:,:,:)
    REAL(wp), ALLOCATABLE                      :: edge2edge_viavert_coeff(:,:,:,:)

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
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: jc,je,jk,jb

    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: all_verts
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
    ALLOCATE(p_coeff%edge2edge_viacell_coeff(nproma,nz_lev,nblks_e,1:2*no_primal_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc failed')
    ENDIF

    ALLOCATE(p_coeff%edge2cell_coeff_cc(nproma,nz_lev,nblks_c,1:no_primal_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc failed')
    ENDIF

    ALLOCATE(p_coeff%edge2cell_coeff_cc_dyn(nproma,1,nblks_c,1:no_primal_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc_dyn failed')
    ENDIF
    !ALLOCATE(p_coeff%edge2vert_coeff_cc_dyn(nproma,1,nblks_v,1:no_dual_edges),stat=ist)
    !IF (ist /= success) THEN
    !  CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff_cc_dyn failed')
    !ENDIF

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
    ALLOCATE(p_coeff%edge2vert_coeff_cc(nproma,nz_lev,nblks_v,1:no_dual_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff failed')
    ENDIF

    ALLOCATE(p_coeff%edge2vert_coeff_cc_t(nproma,nz_lev,nblks_e,1:2),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff failed')
    ENDIF
    ALLOCATE(p_coeff%edge2vert_vector_cc(nproma,nz_lev,nblks_v,1:no_dual_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_vector failed')
    ENDIF

   ALLOCATE(p_coeff%edge2edge_viavert_coeff(nproma,nz_lev,nblks_e,1:2*no_dual_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc failed')
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
    ALLOCATE(p_coeff%variable_vol_norm(nproma,nz_lev,nblks_c,1:no_primal_edges),stat=ist)
    IF (ist /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating variable_vol_norm failed')
    ENDIF

    ALLOCATE(p_coeff%variable_dual_vol_norm(nproma,nz_lev,nblks_v,1:no_dual_edges),stat=ist)
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
      !p_coeff%edge2vert_coeff_cc_dyn%x(ie) = 0._wp
    END DO

    all_cells => p_patch%cells%all
    all_edges => p_patch%edges%all
    all_verts => p_patch%verts%all

    DO jk = 1, nz_lev
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO je =  i_startidx_e, i_endidx_e
!           p_coeff%edge_position_cc(je,jk,jb)             = gc2cc(p_patch%edges%center(je,jb))
          p_coeff%edge_position_cc(je,jk,jb)             = p_patch%edges%cartesian_center(je,jb)
          p_coeff%moved_edge_position_cc(je,jk,jb)%x(:)  = 0._wp
          p_coeff%upwind_cell_position_cc(je,jk,jb)%x(:) = 0._wp
        END DO
      END DO
    END DO

    DO jk = 1, nz_lev
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          p_coeff%cell_position_cc(jc,jk,jb) &
            & = p_patch%cells%cartesian_center(jc,jb)
!             & = gc2cc(p_patch%cells%center(jc,jb))

        END DO
      END DO
    END DO

    p_coeff%edge2edge_viacell_coeff= 0._wp
    p_coeff%edge2edge_viavert_coeff= 0._wp

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
  SUBROUTINE par_init_operator_coeff2( patch, p_patch_3D, p_os, p_phys_param, ocean_coeff)
    !
    TYPE(t_patch),            INTENT(inout)     :: patch
    TYPE(t_patch_3D ),TARGET, INTENT(INOUT) :: p_patch_3D
    TYPE(t_hydro_ocean_state),INTENT(IN)        :: p_os
    TYPE (t_ho_params),       INTENT(IN)        :: p_phys_param
    TYPE(t_operator_coeff),   INTENT(inout)     :: ocean_coeff

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
    REAL(wp) :: edge2edge_viacell_coeff(1:nproma,1:patch%nblks_e,1:6)
    REAL(wp) :: edge2edge_viavert_coeff(1:nproma,1:patch%nblks_e,1:12)
    !-----------------------------------------------------------------------
!     TYPE(t_cartesian_coordinates) :: check_v(nproma, n_zlev, patch%nblks_v, 6)
!     REAL(wp) :: check_r(nproma, n_zlev, patch%nblks_c, 3)
!     REAL(wp) :: max_diff, max_val
    !-----------------------------------------------------------------------
    !ocean_coeff%dist_cell2edge(:,:,:,:) = 0.0_wp

    !CALL init_operator_coeffs( patch, ocean_coeff)

      CALL  par_init_coeff_2D( patch,               &
                          & edge2cell_coeff_cc,   &
                          & edge2cell_coeff_cc_t, &
                          & edge2vert_coeff_cc,   &
                          & edge2vert_coeff_cc_t, &
                          & dist_cell2edge,       &
                          & fixed_vol_norm,       &
                          & variable_vol_norm,    &
                          & variable_dual_vol_norm)

    CALL copy_2D_to_3D_coeff( patch,                  & 
                            & ocean_coeff,            &
                            & edge2edge_viacell_coeff,&
                            & edge2edge_viavert_coeff,&
                            & edge2cell_coeff_cc,     &
                            & edge2cell_coeff_cc_t,   &
                            & edge2vert_coeff_cc,     &
                            & edge2vert_coeff_cc_t,   &
                            & dist_cell2edge,         &
                            & fixed_vol_norm,         &
                            & variable_vol_norm,      &
                            & variable_dual_vol_norm)

     CALL init_geo_factors_oce_3d ( patch, ocean_coeff )
    !---------------------------------------------------------
    CALL par_apply_boundary2coeffs(patch, p_patch_3D, ocean_coeff)


     CALL update_diffusion_matrices( patch, p_patch_3D,               &
                                     & p_os,                          &
                                     & p_phys_param,                  &
                                     & ocean_coeff%matrix_vert_diff_e,&
                                     & ocean_coeff%matrix_vert_diff_c)
  END SUBROUTINE par_init_operator_coeff2
  !-------------------------------------------------------------------------
  !> Initialize 3D expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !! Parellelized by Leonidas Linardakis 2012-3
 
 
   SUBROUTINE update_diffusion_matrices( patch, p_patch_3D,    &
                                         & p_os,               &  
                                         & p_phys_param,       &
                                         & matrix_vert_diff_e, &
                                         & matrix_vert_diff_c)
 
     TYPE(t_patch), TARGET, INTENT(INOUT)  :: patch
     TYPE(t_patch_3D ),TARGET, INTENT(INOUT)   :: p_patch_3D
     TYPE(t_hydro_ocean_state),INTENT(IN) :: p_os
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

          !IF ( v_base%lsm_c(jc,1,jb) <= sea_boundary ) THEN 
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= SEA_BOUNDARY ) THEN 
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

      !IF ( v_base%lsm_e(je,1,jb) <= sea_boundary ) THEN
      IF (  p_patch_3D%lsm_e(je,1,jb) <= SEA_BOUNDARY ) THEN
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

        DO neigbor=1, patch%verts%num_edges(vertex_index, vertex_block) 

          variable_dual_vol_norm(vertex_index, vertex_block, neigbor) = 0.0_wp

          edge_index = patch%verts%edge_idx(vertex_index, vertex_block, neigbor)
          edge_block = patch%verts%edge_blk(vertex_index, vertex_block, neigbor)

          IF (edge_index > 0) THEN
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
  SUBROUTINE init_operator_coeffs( patch, ocean_coeff)
    TYPE(t_patch)    , TARGET, INTENT(INOUT)     :: patch
    TYPE(t_operator_coeff),    INTENT(inout)     :: ocean_coeff

!Local variables
!
    REAL(wp)                      :: prime_edge_length      (1:nproma,patch%nblks_e)
    REAL(wp)                      :: dual_edge_length       (1:nproma,patch%nblks_e)
    REAL(wp)                      :: edge2edge_viacell_coeff(1:nproma,1:patch%nblks_e,1:2*no_primal_edges)
    REAL(wp)                      :: edge2edge_viavert_coeff(1:nproma,1:patch%nblks_e,1:2*no_dual_edges )
    REAL(wp)                      :: dist_cell2edge         (1:nproma,1:patch%nblks_e,1:2)
    REAL(wp)                      :: fixed_vol_norm         (1:nproma,patch%nblks_c)
    REAL(wp)                      :: variable_vol_norm      (1:nproma,1:patch%nblks_c,1:no_primal_edges)
    REAL(wp)                      :: variable_dual_vol_norm (1:nproma,1:patch%nblks_e,1:no_dual_edges)

    REAL(wp)                      :: div_coeff              (1:nproma,1:patch%nblks_c,1:no_primal_edges)
    REAL(wp)                      :: rot_coeff              (1:nproma,1:patch%nblks_v,1:no_dual_edges)
    REAL(wp)                      :: grad_coeff             (1:nproma,1:patch%nblks_e)

    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc     (1:nproma,1:patch%nblks_c,1:no_primal_edges)
    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc_t   (1:nproma,1:patch%nblks_e,1:2)
    TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc     (1:nproma,1:patch%nblks_v,1:no_dual_edges)
    TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc_t   (1:nproma,1:patch%nblks_e,1:2)


    TYPE(t_subset_range), POINTER :: owned_edges         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_cells         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_verts         ! these are the owned entities
    TYPE(t_cartesian_coordinates) :: vertex_position, cell_center, edge_center, vertex_center
    TYPE(t_cartesian_coordinates) :: dist_vector, dist_vector_basic
    TYPE(t_cartesian_coordinates), POINTER :: dual_edge_middle(:,:)
    TYPE(t_cartesian_coordinates) :: coriolis_cartesian_coordinates
    TYPE(t_geographical_coordinates) :: coriolis_geo_coordinates, geo_coordinates
    REAL(wp) :: basin_center_lat_rad, basin_height_rad
    REAL(wp) :: norm, orientation, length
    REAL(wp) :: inverse_sphere_radius
    REAL(wp) :: dist_edge_cell, dist_edge_cell_basic
    INTEGER :: edge_block_cell, edge_index_cell, ictr
    INTEGER :: cell_edge, vert_edge
    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: cell_1_index, cell_1_block, cell_2_index, cell_2_block
    INTEGER :: vertex_1_index, vertex_1_block, vertex_2_index, vertex_2_block
    INTEGER :: level
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
    IF ( MID_POINT_DUAL_EDGE ) THEN
      dual_edge_middle => patch%edges%cartesian_dual_middle
    ELSE
      dual_edge_middle => patch%edges%cartesian_center
    ENDIF

    ! 1) calcultate prima and dual length as cartesian distance
    IF (LARC_LENGTH) THEN

      ! 1a) we just need to get them from the grid
      ! NOTE:  these are earth's distances, translate on a unit sphere
      dist_cell2edge(:,:,:) = &
        & patch%edges%edge_cell_length(:,:,:) * inverse_sphere_radius
      prime_edge_length(:,:) = &
        & patch%edges%primal_edge_length(:,:) * inverse_sphere_radius
      dual_edge_length(:,:) = &
        & patch%edges%dual_edge_length(:,:) * inverse_sphere_radius

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
    ! primal end dual edge lenght have been computed
    !-------------------------------------------



    !-------------------------------------------
    !2) calculate coefficients for difference operators
    !
    !2a) divergence
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        DO neigbor=1, no_primal_edges

          edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x = 0.0_wp
          variable_vol_norm(cell_index, cell_block, neigbor) =  0.0_wp

          edge_index = patch%cells%edge_idx(cell_index, cell_block, neigbor)
          edge_block = patch%cells%edge_blk(cell_index, cell_block, neigbor)

          div_coeff(cell_index,cell_block,neigbor) =                           &
              & patch%edges%primal_edge_length(edge_index,edge_block) *        &
              & patch%cells%edge_orientation(cell_index,cell_block,neigbor)  / &
              & patch%cells%area(cell_index,cell_block)

        ENDDO !neigbor=1,patch%cell_type
      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
 

   !2b) gradient
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        grad_coeff(edge_index,edge_block)&
        & =1.0_wp/ dual_edge_length(edge_index,edge_block)!patch%edges%inv_dual_edge_length(edge_index, edge_block)

      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
    

   !2c) curl coefficients
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index

        vertex_position%x = patch%verts%cartesian(vertex_index, vertex_block)%x

        DO neigbor=1, patch%verts%num_edges(vertex_index, vertex_block)  ! no_dual_edges

          edge_index = patch%verts%edge_idx(vertex_index, vertex_block, neigbor)
          edge_block = patch%verts%edge_blk(vertex_index, vertex_block, neigbor)

!           IF (edge_index > 0) THEN
          rot_coeff(vertex_index,vertex_block,neigbor)           &
            &= patch%edges%dual_edge_length(edge_index,edge_block) &
            &* patch%verts%edge_orientation(vertex_index,vertex_block,neigbor)
!           ENDIF
        ENDDO !neigbor=1,6
      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
   
  !Copy coefficients to 3D 
   DO level=1,n_zlev
     ocean_coeff%div_coeff(:,level,:,:) = div_coeff(:,:,:)
     ocean_coeff%rot_coeff(:,level,:,:) = rot_coeff(:,:,:)
     ocean_coeff%grad_coeff(:,level,:)  = grad_coeff(:,:)
   END DO
   !-------------------
   ! sync the results
   CALL sync_patch_array(SYNC_E, patch, ocean_coeff%grad_coeff(:,:,:))
   DO neigbor=1,no_primal_edges
     CALL sync_patch_array(SYNC_C, patch, ocean_coeff%div_coeff(:,:,:,neigbor))
   END DO
   DO neigbor=1,no_dual_edges
     CALL sync_patch_array(SYNC_V, patch, ocean_coeff%rot_coeff(:,:,:,neigbor))
   END DO

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
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        cell_center%x = patch%cells%cartesian_center(cell_index, cell_block)%x
        fixed_vol_norm(cell_index,cell_block) = 0.0_wp

        !-------------------------------
        DO neigbor=1, no_primal_edges

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
              & patch%cells%edge_orientation(cell_index,cell_block,neigbor)

            fixed_vol_norm(cell_index,cell_block) = &
              & fixed_vol_norm(cell_index,cell_block) + &
              & 0.5_wp * norm * prime_edge_length(edge_index,edge_block)

            variable_vol_norm(cell_index, cell_block, neigbor) = &
              & 0.5_wp * norm * prime_edge_length(edge_index,edge_block)

          ENDIF !(edge_block > 0 )
          !div_coeff(cell_index,cell_block,neigbor) =                           &
          !    & patch%edges%primal_edge_length(edge_index,edge_block) *        &
          !    & patch%cells%edge_orientation(cell_index,cell_block,neigbor)  / &
          !    & patch%cells%area(cell_index,cell_block)

        ENDDO !neigbor=1,patch%cell_type
        !-------------------------------
      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block


    !-------------------
    ! sync the results
    CALL sync_patch_array(SYNC_C, patch, fixed_vol_norm(:,:))
    DO neigbor=1,no_primal_edges
      CALL sync_patch_array(SYNC_C, patch, edge2cell_coeff_cc(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_C, patch, edge2cell_coeff_cc(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_C, patch, edge2cell_coeff_cc(:,:,neigbor)%x(3))
      CALL sync_patch_array(SYNC_C, patch, variable_vol_norm(:,:,neigbor))
    ENDDO
    !-------------------

   !copy calculated 2D arrays to 3D structure
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      DO level = 1, n_zlev

       ocean_coeff%fixed_vol_norm(:,level,cell_block) = fixed_vol_norm(:,cell_block)

       DO neigbor=1,no_primal_edges

         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(1)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(1)

         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(2)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(2)

         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(3)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(3)

         ocean_coeff%variable_vol_norm(:,level,cell_block,neigbor)  &
           &= variable_vol_norm(:,cell_block,neigbor)

        ENDDO ! neigbor=1,patch%cell_type
      ENDDO  !  level = 1, n_zlev
    ENDDO ! cell_block

    CALL sync_patch_array(SYNC_C, patch, ocean_coeff%fixed_vol_norm(:,:,:))
    DO neigbor=1,no_primal_edges
      CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc(:,:,:,neigbor)%x(3))
      CALL sync_patch_array(SYNC_C, patch, ocean_coeff%variable_vol_norm(:,:,:,neigbor))
    ENDDO




    !-------------------------------------------
    ! 4) compute:
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

        !grad_coeff(edge_index,edge_block)&
        !& =1.0_wp/ dual_edge_length(edge_index,edge_block)!patch%edges%inv_dual_edge_length(edge_index, edge_block)

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

    !copy 2D to 3D structure
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          ocean_coeff%edge2cell_coeff_cc_t(edge_index,level,edge_block,1)%x &
          &= edge2cell_coeff_cc_t(edge_index,edge_block,1)%x

          ocean_coeff%edge2cell_coeff_cc_t(edge_index,level,edge_block,2)%x &
          &= edge2cell_coeff_cc_t(edge_index,edge_block,2)%x

        ENDDO
      ENDDO
    ENDDO
    ! sync the results
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2
    !-------------------------------------------


    !-------------------------------------------
    ! 5) compute
    !   calculate edge2edge_viacell_coeff
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

       edge_center%x = patch%edges%cartesian_center(edge_index, edge_block)%x

        ictr=0
        DO neigbor=1,2

          cell_index    = patch%edges%cell_idx(edge_index, edge_block, neigbor)
          cell_block    = patch%edges%cell_blk(edge_index, edge_block, neigbor)         
          cell_center%x = patch%cells%cartesian_center(cell_index, cell_block)%x

          !dist_vector_basic%x = edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x
          dist_vector_basic%x = edge_center%x - cell_center%x

          dist_edge_cell_basic  = SQRT(SUM( dist_vector_basic%x * dist_vector_basic%x))
          dist_vector_basic%x = dist_vector_basic%x/dist_edge_cell_basic

          orientation = DOT_PRODUCT(dist_vector_basic%x, &
          & patch%edges%primal_cart_normal(edge_index,edge_block)%x)
          IF (orientation < 0.0_wp) dist_vector_basic%x = - dist_vector_basic%x

          IF (cell_block > 0) THEN

            !loop over the edges of neighbor 1 and 2
            DO cell_edge=1,patch%cell_type

              ictr=ictr+1
              !actual edge
              edge_index_cell = patch%cells%edge_idx(cell_index, cell_block, cell_edge)
              edge_block_cell = patch%cells%edge_blk(cell_index, cell_block, cell_edge) 

              !dist_vector%x = edge2cell_coeff_cc(cell_index,cell_block,cell_edge)%x
              dist_vector%x =  patch%edges%cartesian_center(edge_index_cell, edge_block_cell)%x  &
              & -cell_center%x

              dist_edge_cell  = SQRT(SUM( dist_vector%x * dist_vector%x))
              dist_vector%x = dist_vector%x/dist_edge_cell
              dist_vector%x = dist_vector%x*patch%cells%edge_orientation(cell_index,cell_block,cell_edge)

              !This is the cosine of the angle between vectors from cell center
              !to cell edges 
              edge2edge_viacell_coeff(edge_index,edge_block,ictr)&
              & =DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)

              !IF(abs(edge2edge_viacell_coeff(edge_index,edge_block,ictr)-1.0_wp)<1.0E-6_wp)THEN
              !  write(*,*)'ran into'
                !edge2edge_viacell_coeff(edge_index,edge_block,ictr)=1.0_wp
              !ENDIF

              !multiply the cosine by length and orientation and divide by
              !dual length
                edge2edge_viacell_coeff(edge_index,edge_block,ictr)=        &
                &edge2edge_viacell_coeff(edge_index,edge_block,ictr)        &
                &*prime_edge_length(edge_index_cell,edge_block_cell)        &
                &* dist_edge_cell *dist_edge_cell_basic                     &
                &/dual_edge_length(edge_index, edge_block)                 

! IF(edge_index==1.and.edge_block==1)THEN
! write(123,*)'actual angle',neigbor, edge_index_cell, edge_block_cell,ictr,&
! & edge2edge_viacell_coeff(edge_index,edge_block,ictr),&
! &DOT_PRODUCT(edge2cell_coeff_cc(cell_index,cell_block,cell_edge)%x,edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x),&
! &acos(edge2edge_viacell_coeff(edge_index,edge_block,ictr))*rad2deg
! !IF(edge_index_cell==edge_index.and.edge_block_cell==edge_block)THEN
! !write(123,*)'vecs',neigbor,dist_vector_basic%x,dist_vector%x
! !ENDIF
! ENDIF
            END DO
          ENDIF ! (cell_block > 0)
        ENDDO ! neigbor=1,2
      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block

    DO neigbor=1, 2*no_primal_edges
      CALL sync_patch_array(SYNC_E, patch, edge2edge_viacell_coeff(:,:,neigbor))
    ENDDO
    !copy 2D to 3D structure
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          ocean_coeff%edge2edge_viacell_coeff(edge_index,level,edge_block,1:2*no_primal_edges) &
          &= edge2edge_viacell_coeff(edge_index,edge_block,1:2*no_primal_edges)

        ENDDO
      ENDDO
    ENDDO
    DO neigbor=1, 2*no_primal_edges
      CALL sync_patch_array(SYNC_E, patch,  ocean_coeff%edge2edge_viacell_coeff(:,:,:,neigbor))
    ENDDO
   !-------------------------------------------
!Do ictr=1,6
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


    !-------------------------------------------
    ! 6) compute:
    !   edge2vert_coeff_cc
    !   variable_dual_vol_norm
    variable_dual_vol_norm(:,:,:)  = 0.0_wp
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index

        vertex_position%x = patch%verts%cartesian(vertex_index, vertex_block)%x

        DO neigbor=1, no_dual_edges

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
              & dual_edge_length(edge_index, edge_block) 
          ENDIF !(edge_block > 0) THEN

          !rot_coeff(vertex_index,vertex_block,neigbor)     &
          !    &= patch%edges%dual_edge_length(edge_index,edge_block) * &
          !    & patch%verts%edge_orientation(vertex_index,vertex_block,neigbor)

        ENDDO !neigbor=1,6
      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
    !-------------------
    ! sync the results
    DO neigbor=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch, edge2vert_coeff_cc(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_V, patch, edge2vert_coeff_cc(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_V, patch, edge2vert_coeff_cc(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,6
    ! edge2vert_coeff_cc is computed

    !copy 2D to 3D structure
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      DO level = 1, n_zlev
        DO neigbor=1,no_dual_edges

          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(1)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(1)

          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(2)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(2)

          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(3)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(3)

!          ocean_coeff%variable_dual_vol_norm(:,level,vertex_block,neigbor)&
!          &=variable_dual_vol_norm(:,vertex_block,neigbor)
        ENDDO ! neigbor=1,patch%cell_type
      ENDDO  !  level = 1, n_zlev
    ENDDO ! vertex_block
    ! sync the results
    DO neigbor=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,:,:,neigbor)%x(3))
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


    DO edge_block = owned_edges%start_block, owned_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          ocean_coeff%edge2vert_coeff_cc_t(edge_index,level,edge_block,1)%x &
          &=  edge2vert_coeff_cc_t(edge_index,edge_block,1)%x

          ocean_coeff%edge2vert_coeff_cc_t(edge_index,level,edge_block,2)%x &
          &=  edge2vert_coeff_cc_t(edge_index,edge_block,2)%x

        ENDDO
      ENDDO
    ENDDO
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2vert_coeff_cc_t(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2vert_coeff_cc_t(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2vert_coeff_cc_t(:,:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2

    !----------------------------------------------------

   !----------------------------------------------------
   ! 8) compute
   !   edge2edge_viavert calculation
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge_center%x = dual_edge_middle(edge_index, edge_block)%x

        ictr= 0
        DO neigbor=1,2

          vertex_index   = patch%edges%vertex_idx(edge_index, edge_block, neigbor)
          vertex_block   = patch%edges%vertex_blk(edge_index, edge_block, neigbor)
          vertex_center%x= patch%verts%cartesian(vertex_index, vertex_block)%x

          dist_vector_basic%x = edge_center%x - vertex_center%x

          dist_vector_basic%x = dist_vector_basic%x * patch%edges%system_orientation(edge_index, edge_block)

          DO vert_edge=1,no_dual_edges
            ictr=ictr+1
            !actual edge
            edge_index_cell = patch%verts%edge_idx(vertex_index, vertex_block, vert_edge)
            edge_block_cell = patch%verts%edge_blk(vertex_index, vertex_block, vert_edge) 
            dist_vector%x  =  dual_edge_middle(edge_index_cell, edge_block_cell)%x &
                           & -vertex_center%x

            dist_vector = vector_product(dist_vector, dual_edge_middle(edge_index_cell, edge_block_cell))
            orientation = DOT_PRODUCT( dist_vector%x,                         &
               & patch%edges%primal_cart_normal(edge_index_cell, edge_block_cell)%x)
            IF (orientation < 0.0_wp) dist_vector%x = - dist_vector%x

            !This is the cosine of the angle between vectors from dual cell centers
            !to dual cell edges 
            IF(neigbor==1)THEN
            edge2edge_viavert_coeff(edge_index,edge_block,ictr)&
              & =DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)

            ELSEIF(neigbor==2)THEN
            edge2edge_viavert_coeff(edge_index,edge_block,ictr)&
              & =-DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)
            ENDIF
!edge2edge_viavert_coeff(edge_index,edge_block,ictr)&
!&=DOT_PRODUCT(edge2vert_coeff_cc(vertex_index, vertex_block, vert_edge)%x,&
!&edge2vert_coeff_cc_t(edge_index, edge_block, neigbor)%x)
! IF(edge_index==1 .and. edge_block==1)THEN
! write(*,*)'oce coeff',ictr,edge2edge_viavert_coeff(edge_index,edge_block,ictr),&
! &acos(edge2edge_viavert_coeff(edge_index,edge_block,ictr))*rad2deg, edge_index,edge_block,edge_index_cell,edge_block_cell
! ENDIF
            edge2edge_viavert_coeff(edge_index,edge_block,ictr)&
            &=edge2edge_viavert_coeff(edge_index,edge_block,ictr)&
            &*(dual_edge_length(edge_index_cell, edge_block_cell) &
            &/prime_edge_length(edge_index, edge_block))

          END DO
        ENDDO !neigbor=1,2
      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
    !-------------------
   DO neigbor=1,2*no_dual_edges
   CALL sync_patch_array(SYNC_V, patch, edge2edge_viavert_coeff(:,:,neigbor))
   END DO

    DO edge_block = owned_edges%start_block, owned_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          ocean_coeff%edge2edge_viavert_coeff(edge_index,level,edge_block,1:2*no_dual_edges)=&
          &edge2edge_viavert_coeff(edge_index,edge_block,1:2*no_dual_edges)

        ENDDO
      ENDDO
    ENDDO
   DO neigbor=1,2*no_dual_edges
   CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2edge_viavert_coeff(:,:,:,neigbor))
   END DO
   !-------------------------------------------


!Do ictr=1,12
 ! write(*,*)'max coeff',&
! & maxval(edge2edge_viacell_coeff(:,:,1)),&
! & maxval(edge2edge_viacell_coeff(:,:,2)),&
! & maxval(edge2edge_viacell_coeff(:,:,3)),&
! & maxval(edge2edge_viacell_coeff(:,:,4)),&
! & maxval(edge2edge_viacell_coeff(:,:,5)),&
! & maxval(edge2edge_viacell_coeff(:,:,6))
! 
! write(*,*)'min coeff',&
! & minval(edge2edge_viacell_coeff(:,:,1)),&
! & minval(edge2edge_viacell_coeff(:,:,2)),&
! & minval(edge2edge_viacell_coeff(:,:,3)),&
! & minval(edge2edge_viacell_coeff(:,:,4)),&
! & minval(edge2edge_viacell_coeff(:,:,5)),&
! & minval(edge2edge_viacell_coeff(:,:,6))
!  
! write(*,*)'max coeff',&
!  & maxval(edge2edge_viavert_coeff(:,:,ictr)),&
!  & minval(edge2edge_viavert_coeff(:,:,ictr))
! write(*,*)'max angle',ictr,&
! & maxval(acos(edge2edge_viavert_coeff(:,:,ictr))*rad2deg)
! 
! write(*,*)'min angle',&
! & minval(acos(edge2edge_viavert_coeff(:,:,ictr))*rad2deg)
! END DO
!-------------------------------------------



    !----------------------------------------------------
    ! 9) recalculate the coriolis coefficient
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
  !!  Parallelized by Leonidas Linardakis, 2012-3
  SUBROUTINE init_operator_coeffs_cell( patch, ocean_coeff,prime_edge_length, dual_edge_length )
    TYPE(t_patch), TARGET,  INTENT(INOUT)  :: patch
    TYPE(t_operator_coeff), INTENT(inout)  :: ocean_coeff
    REAL(wp),               INTENT(IN)     :: prime_edge_length(1:nproma,1:patch%nblks_e)
    REAL(wp),               INTENT(IN)     :: dual_edge_length (1:nproma,1:patch%nblks_e)

!Local variables
!
    REAL(wp)                      :: edge2edge_viacell_coeff(1:nproma,1:patch%nblks_e,1:2*no_primal_edges)
    REAL(wp)                      :: dist_cell2edge         (1:nproma,1:patch%nblks_e,1:2)
    REAL(wp)                      :: fixed_vol_norm         (1:nproma,patch%nblks_c)
    REAL(wp)                      :: variable_vol_norm      (1:nproma,1:patch%nblks_c,1:no_primal_edges)
    REAL(wp)                      :: norm, orientation, length
    REAL(wp)                      :: dist_edge_cell, dist_edge_cell_basic

    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc     (1:nproma,1:patch%nblks_c,1:no_primal_edges)
    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc_t   (1:nproma,1:patch%nblks_e,1:2)
    TYPE(t_cartesian_coordinates) :: cell_center, edge_center
    TYPE(t_cartesian_coordinates) :: dist_vector, dist_vector_basic

    INTEGER :: edge_block_cell, edge_index_cell, ictr
    INTEGER :: cell_edge
    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: level

    TYPE(t_subset_range), POINTER :: owned_edges         
    TYPE(t_subset_range), POINTER :: owned_cells        
    !-----------------------------------------------------------------------
    owned_edges => patch%edges%owned
    owned_cells => patch%cells%owned
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

    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        cell_center%x = patch%cells%cartesian_center(cell_index, cell_block)%x
        fixed_vol_norm(cell_index,cell_block) = 0.0_wp

        !-------------------------------
        DO neigbor=1, no_primal_edges

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
              & patch%cells%edge_orientation(cell_index,cell_block,neigbor)

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

   !copy calculated 2D arrays to 3D structure
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      DO level = 1, n_zlev

       ocean_coeff%fixed_vol_norm(:,level,cell_block) = fixed_vol_norm(:,cell_block)

       DO neigbor=1,no_primal_edges

         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(1)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(1)

         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(2)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(2)

         ocean_coeff%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(3)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(3)

         ocean_coeff%variable_vol_norm(:,level,cell_block,neigbor)  &
           &= variable_vol_norm(:,cell_block,neigbor)

        ENDDO ! neigbor=1,patch%cell_type
      ENDDO  !  level = 1, n_zlev
    ENDDO ! cell_block
    CALL sync_patch_array(SYNC_C, patch, ocean_coeff%fixed_vol_norm(:,:,:))
    DO neigbor=1,no_primal_edges
      CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_C, patch, ocean_coeff%edge2cell_coeff_cc(:,:,:,neigbor)%x(3))
      CALL sync_patch_array(SYNC_C, patch, ocean_coeff%variable_vol_norm(:,:,:,neigbor))
    ENDDO

    !-------------------------------------------
    ! 4) compute:
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

    !copy 2D to 3D structure
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          ocean_coeff%edge2cell_coeff_cc_t(edge_index,level,edge_block,1)%x &
          &= edge2cell_coeff_cc_t(edge_index,edge_block,1)%x

          ocean_coeff%edge2cell_coeff_cc_t(edge_index,level,edge_block,2)%x &
          &= edge2cell_coeff_cc_t(edge_index,edge_block,2)%x

        ENDDO
      ENDDO
    ENDDO
    ! sync the results
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2
    !-------------------------------------------


    !-------------------------------------------
    ! 5) compute
    !   calculate edge2edge_viacell_coeff
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

       edge_center%x = patch%edges%cartesian_center(edge_index, edge_block)%x

        ictr=0
        DO neigbor=1,2

          cell_index    = patch%edges%cell_idx(edge_index, edge_block, neigbor)
          cell_block    = patch%edges%cell_blk(edge_index, edge_block, neigbor)         
          cell_center%x = patch%cells%cartesian_center(cell_index, cell_block)%x

          !dist_vector_basic%x = edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x
          dist_vector_basic%x = edge_center%x - cell_center%x

          dist_edge_cell_basic  = SQRT(SUM( dist_vector_basic%x * dist_vector_basic%x))
          dist_vector_basic%x = dist_vector_basic%x/dist_edge_cell_basic

          orientation = DOT_PRODUCT(dist_vector_basic%x, &
          & patch%edges%primal_cart_normal(edge_index,edge_block)%x)
          IF (orientation < 0.0_wp) dist_vector_basic%x = - dist_vector_basic%x

          IF (cell_block > 0) THEN

            !loop over the edges of neighbor 1 and 2
            DO cell_edge=1,patch%cell_type

              ictr=ictr+1
              !actual edge
              edge_index_cell = patch%cells%edge_idx(cell_index, cell_block, cell_edge)
              edge_block_cell = patch%cells%edge_blk(cell_index, cell_block, cell_edge) 

              !dist_vector%x = edge2cell_coeff_cc(cell_index,cell_block,cell_edge)%x
              dist_vector%x =  patch%edges%cartesian_center(edge_index_cell, edge_block_cell)%x  &
              & -cell_center%x

              dist_edge_cell  = SQRT(SUM( dist_vector%x * dist_vector%x))
              dist_vector%x = dist_vector%x/dist_edge_cell
              dist_vector%x = dist_vector%x*patch%cells%edge_orientation(cell_index,cell_block,cell_edge)

              !This is the cosine of the angle between vectors from cell center
              !to cell edges 
              edge2edge_viacell_coeff(edge_index,edge_block,ictr)&
              & =DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)

              !IF(abs(edge2edge_viacell_coeff(edge_index,edge_block,ictr)-1.0_wp)<1.0E-6_wp)THEN
              !  write(*,*)'ran into'
                !edge2edge_viacell_coeff(edge_index,edge_block,ictr)=1.0_wp
              !ENDIF

              !multiply the cosine by length and orientation and divide by
              !dual length
                edge2edge_viacell_coeff(edge_index,edge_block,ictr)=        &
                &edge2edge_viacell_coeff(edge_index,edge_block,ictr)        &
                &*prime_edge_length(edge_index_cell,edge_block_cell)        &
                &* dist_edge_cell *dist_edge_cell_basic                     &
                &/dual_edge_length(edge_index, edge_block)                 

! IF(edge_index==1.and.edge_block==1)THEN
! write(123,*)'actual angle',neigbor, edge_index_cell, edge_block_cell,ictr,&
! & edge2edge_viacell_coeff(edge_index,edge_block,ictr),&
! &DOT_PRODUCT(edge2cell_coeff_cc(cell_index,cell_block,cell_edge)%x,edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x),&
! &acos(edge2edge_viacell_coeff(edge_index,edge_block,ictr))*rad2deg
! !IF(edge_index_cell==edge_index.and.edge_block_cell==edge_block)THEN
! !write(123,*)'vecs',neigbor,dist_vector_basic%x,dist_vector%x
! !ENDIF
! ENDIF
            END DO
          ENDIF ! (cell_block > 0)
        ENDDO ! neigbor=1,2
      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block

    CALL sync_patch_array(SYNC_E, patch, edge2edge_viacell_coeff)

    !copy 2D to 3D structure
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          ocean_coeff%edge2edge_viacell_coeff(edge_index,level,edge_block,1:2*no_primal_edges) &
          &= edge2edge_viacell_coeff(edge_index,edge_block,1:2*no_primal_edges)

        ENDDO
      ENDDO
    ENDDO    
    DO neigbor=1, 2*no_primal_edges
      CALL sync_patch_array(SYNC_E, patch,  ocean_coeff%edge2edge_viacell_coeff(:,:,:,neigbor))
    ENDDO
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
  !!  Parallelized by Leonidas Linardakis, 2012-3
  SUBROUTINE init_operator_coeffs_vertex( patch, ocean_coeff, prime_edge_length, dual_edge_length)
    TYPE(t_patch), TARGET, INTENT(INOUT) :: patch
    TYPE(t_operator_coeff),INTENT(INOUT) :: ocean_coeff
    REAL(wp), INTENT(IN)                 :: prime_edge_length(1:nproma,1:patch%nblks_e)
    REAL(wp), INTENT(IN)                 :: dual_edge_length (1:nproma,1:patch%nblks_e)

!Local variables
!
    TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc     (1:nproma,1:patch%nblks_v,1:no_dual_edges)
    TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc_t   (1:nproma,1:patch%nblks_e,1:2)
    TYPE(t_cartesian_coordinates) :: vertex_position, cell_center, edge_center, vertex_center
    TYPE(t_cartesian_coordinates) :: dist_vector, dist_vector_basic
    TYPE(t_cartesian_coordinates), POINTER :: dual_edge_middle(:,:)

    REAL(wp)                      :: edge2edge_viavert_coeff(1:nproma,1:patch%nblks_e,1:2*no_dual_edges )
    !REAL(wp)                      :: variable_dual_vol_norm (1:nproma,1:patch%nblks_e,1:no_dual_edges)
    REAL(wp)                      :: norm, orientation, length

    INTEGER :: ictr,edge_block_cell, edge_index_cell
    INTEGER :: vert_edge
    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: level

    TYPE(t_subset_range), POINTER :: owned_edges        
    TYPE(t_subset_range), POINTER :: owned_verts         
    !-----------------------------------------------------------------------
    owned_edges => patch%edges%owned
    owned_verts => patch%verts%owned

    edge2vert_coeff_cc(:,:,:)%x(1) = 0.0_wp
    edge2vert_coeff_cc(:,:,:)%x(2) = 0.0_wp
    edge2vert_coeff_cc(:,:,:)%x(3) = 0.0_wp

    edge2vert_coeff_cc_t(:,:,:)%x(1) = 0.0_wp
    edge2vert_coeff_cc_t(:,:,:)%x(2) = 0.0_wp
    edge2vert_coeff_cc_t(:,:,:)%x(3) = 0.0_wp

   !-------------------------------------------
    ! 6) compute:
    !   edge2vert_coeff_cc
    !   variable_dual_vol_norm will be handled in boundary_coeff-sbr
    !variable_dual_vol_norm(:,:,:)  = 0.0_wp

    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index

        vertex_position%x = patch%verts%cartesian(vertex_index, vertex_block)%x

        DO neigbor=1, no_dual_edges

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
              & dual_edge_length(edge_index, edge_block) 
          ENDIF !(edge_block > 0) THEN

          !rot_coeff(vertex_index,vertex_block,neigbor)     &
          !    &= patch%edges%dual_edge_length(edge_index,edge_block) * &
          !    & patch%verts%edge_orientation(vertex_index,vertex_block,neigbor)

        ENDDO !neigbor=1,6
      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
    !-------------------
    ! sync the results
    DO neigbor=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch, edge2vert_coeff_cc(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_V, patch, edge2vert_coeff_cc(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_V, patch, edge2vert_coeff_cc(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,6
    ! edge2vert_coeff_cc is computed

    !copy 2D to 3D structure
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      DO level = 1, n_zlev
        DO neigbor=1,no_dual_edges

          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(1)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(1)

          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(2)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(2)

          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(3)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(3)

!          ocean_coeff%variable_dual_vol_norm(:,level,vertex_block,neigbor)&
!          &=variable_dual_vol_norm(:,vertex_block,neigbor)
        ENDDO ! neigbor=1,patch%cell_type
      ENDDO  !  level = 1, n_zlev
    ENDDO ! vertex_block
    DO neigbor=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,:,:,neigbor)%x(3))
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


    DO edge_block = owned_edges%start_block, owned_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          ocean_coeff%edge2vert_coeff_cc_t(edge_index,level,edge_block,1)%x &
          &=  edge2vert_coeff_cc_t(edge_index,edge_block,1)%x

          ocean_coeff%edge2vert_coeff_cc_t(edge_index,level,edge_block,2)%x &
          &=  edge2vert_coeff_cc_t(edge_index,edge_block,2)%x

        ENDDO
      ENDDO
    ENDDO
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch, edge2vert_coeff_cc_t(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch, edge2vert_coeff_cc_t(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch, edge2vert_coeff_cc_t(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2

    !----------------------------------------------------

   !----------------------------------------------------
   ! 8) compute
   !   edge2edge_viavert calculation
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge_center%x = dual_edge_middle(edge_index, edge_block)%x

        ictr= 0
        DO neigbor=1,2

          vertex_index   = patch%edges%vertex_idx(edge_index, edge_block, neigbor)
          vertex_block   = patch%edges%vertex_blk(edge_index, edge_block, neigbor)
          vertex_center%x= patch%verts%cartesian(vertex_index, vertex_block)%x

          dist_vector_basic%x = edge_center%x - vertex_center%x

          dist_vector_basic%x = dist_vector_basic%x * patch%edges%system_orientation(edge_index, edge_block)

          DO vert_edge=1,no_dual_edges
            ictr=ictr+1
            !actual edge
            edge_index_cell = patch%verts%edge_idx(vertex_index, vertex_block, vert_edge)
            edge_block_cell = patch%verts%edge_blk(vertex_index, vertex_block, vert_edge) 
            dist_vector%x  =  dual_edge_middle(edge_index_cell, edge_block_cell)%x &
                           & -vertex_center%x

            dist_vector = vector_product(dist_vector, dual_edge_middle(edge_index_cell, edge_block_cell))
            orientation = DOT_PRODUCT( dist_vector%x,                         &
               & patch%edges%primal_cart_normal(edge_index_cell, edge_block_cell)%x)
            IF (orientation < 0.0_wp) dist_vector%x = - dist_vector%x

            !This is the cosine of the angle between vectors from dual cell centers
            !to dual cell edges 
            IF(neigbor==1)THEN
            edge2edge_viavert_coeff(edge_index,edge_block,ictr)&
              & =DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)

            ELSEIF(neigbor==2)THEN
            edge2edge_viavert_coeff(edge_index,edge_block,ictr)&
              & =-DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)
            ENDIF
!edge2edge_viavert_coeff(edge_index,edge_block,ictr)&
!&=DOT_PRODUCT(edge2vert_coeff_cc(vertex_index, vertex_block, vert_edge)%x,&
!&edge2vert_coeff_cc_t(edge_index, edge_block, neigbor)%x)
! IF(edge_index==1 .and. edge_block==1)THEN
! write(*,*)'oce coeff',ictr,edge2edge_viavert_coeff(edge_index,edge_block,ictr),&
! &acos(edge2edge_viavert_coeff(edge_index,edge_block,ictr))*rad2deg, edge_index,edge_block,edge_index_cell,edge_block_cell
! ENDIF
            edge2edge_viavert_coeff(edge_index,edge_block,ictr)&
            &=edge2edge_viavert_coeff(edge_index,edge_block,ictr)&
            &*(dual_edge_length(edge_index_cell, edge_block_cell) &
            &/prime_edge_length(edge_index, edge_block))

          END DO
        ENDDO !neigbor=1,2
      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
    !-------------------
   CALL sync_patch_array(SYNC_V, patch, edge2edge_viavert_coeff)

    DO edge_block = owned_edges%start_block, owned_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          ocean_coeff%edge2edge_viavert_coeff(edge_index,level,edge_block,1:2*no_dual_edges)&
          &=edge2edge_viavert_coeff(edge_index,edge_block,1:2*no_dual_edges)

        ENDDO
      ENDDO
    ENDDO
   DO neigbor=1,2*no_dual_edges
     CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2edge_viavert_coeff(:,:,:,neigbor))
   END DO
   !-------------------------------------------

!Do ictr=1,12
 ! write(*,*)'max coeff',&
! & maxval(edge2edge_viacell_coeff(:,:,1)),&
! & maxval(edge2edge_viacell_coeff(:,:,2)),&
! & maxval(edge2edge_viacell_coeff(:,:,3)),&
! & maxval(edge2edge_viacell_coeff(:,:,4)),&
! & maxval(edge2edge_viacell_coeff(:,:,5)),&
! & maxval(edge2edge_viacell_coeff(:,:,6))
! 
! write(*,*)'min coeff',&
! & minval(edge2edge_viacell_coeff(:,:,1)),&
! & minval(edge2edge_viacell_coeff(:,:,2)),&
! & minval(edge2edge_viacell_coeff(:,:,3)),&
! & minval(edge2edge_viacell_coeff(:,:,4)),&
! & minval(edge2edge_viacell_coeff(:,:,5)),&
! & minval(edge2edge_viacell_coeff(:,:,6))
!  
! write(*,*)'max coeff',&
!  & maxval(edge2edge_viavert_coeff(:,:,ictr)),&
!  & minval(edge2edge_viavert_coeff(:,:,ictr))
! write(*,*)'max angle',ictr,&
! & maxval(acos(edge2edge_viavert_coeff(:,:,ictr))*rad2deg)
! 
! write(*,*)'min angle',&
! & minval(acos(edge2edge_viavert_coeff(:,:,ictr))*rad2deg)
! END DO
!-------------------------------------------



  END SUBROUTINE init_operator_coeffs_vertex
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------

  !> Initialize 3D expansion coefficients.
  !!
  !! @par Revision History
  !! Peter Korn (2012-2)
  !! Parellelized by Leonidas Linardakis 2012-3


  SUBROUTINE copy_2D_to_3D_coeff( patch,                 &
                               & ocean_coeff,            &
                               & edge2edge_viacell_coeff,&
                               & edge2edge_viavert_coeff,&
                               & edge2cell_coeff_cc,     &
                               & edge2cell_coeff_cc_t,   &
                               & edge2vert_coeff_cc,     &
                               & edge2vert_coeff_cc_t,   &
                               & dist_cell2edge,         &
                               & fixed_vol_norm,         &
                               & variable_vol_norm,      &
                               & variable_dual_vol_norm)
    TYPE(t_patch)    , TARGET, INTENT(INOUT)     :: patch
    TYPE(t_operator_coeff), INTENT(inout) :: ocean_coeff
    REAL(wp), INTENT(IN)                      :: edge2edge_viacell_coeff(1:nproma,1:patch%nblks_e,1:2*no_primal_edges)
    REAL(wp), INTENT(IN)                      :: edge2edge_viavert_coeff(1:nproma,1:patch%nblks_e,1:2*no_dual_edges)
    TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2cell_coeff_cc     (1:nproma,1:patch%nblks_c,1:patch%cell_type)
    TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2cell_coeff_cc_t   (1:nproma,1:patch%nblks_e,1:2)
    TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2vert_coeff_cc     (1:nproma,1:patch%nblks_v,1:no_dual_edges)
    TYPE(t_cartesian_coordinates), INTENT(IN) :: edge2vert_coeff_cc_t   (1:nproma,1:patch%nblks_e,1:2)
    REAL(wp), INTENT(IN) :: dist_cell2edge        (1:nproma,1:patch%nblks_e,1:2)
    REAL(wp), INTENT(IN) :: fixed_vol_norm        (1:nproma,patch%nblks_c)
    REAL(wp), INTENT(IN) :: variable_vol_norm     (1:nproma,1:patch%nblks_c,1:no_primal_edges)
    REAL(wp), INTENT(IN) :: variable_dual_vol_norm(1:nproma,1:patch%nblks_v,1:no_dual_edges)
    !
    !Local variables
    !
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

          !ocean_coeff%edge2edge_viacell_coeff(je,level,edge_block,1:6)=&
          !&edge2edge_viacell_coeff(je,edge_block,1:6)

          !ocean_coeff%edge2edge_viavert_coeff(je,level,edge_block,1:12)=&
          !&edge2edge_viavert_coeff(je,edge_block,1:12)

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
        DO neigbor=1,no_dual_edges

          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(1) &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(1)
          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(2) &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(2)
          ocean_coeff%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(3) &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(3)

         ocean_coeff%variable_dual_vol_norm(:,level,vertex_block,neigbor)&
         &=variable_dual_vol_norm(:,vertex_block,neigbor)
        ENDDO ! neigbor=1,patch%cell_type
      ENDDO  !  level = 1, n_zlev
    ENDDO ! vertex_block

   !DO neigbor=1,2*no_dual_edges
   !  CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2edge_viavert_coeff(:,:,:,neigbor))
   !END DO

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
    TYPE(t_patch), TARGET,  INTENT(INOUT)       :: patch
    TYPE(t_patch_3D ),TARGET, INTENT(INOUT) :: p_patch_3D
    TYPE(t_operator_coeff), INTENT(INOUT)       :: ocean_coeff

    !Local variables
    INTEGER :: jk, jc, jb, je, ibe, ile, jev, jv, iiv,iib
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_v, i_endidx_v

    INTEGER :: i_v_ctr(nproma,n_zlev,patch%nblks_v)
    INTEGER :: ibnd_edge_idx(4), ibnd_edge_blk(4)  !maximal 4 boundary edges in a dual loop.
    INTEGER :: i_edge_idx(4)
    !REAL(wp) :: z_orientation(4)!,z_area_scaled 
    REAL(wp) :: zarea_fraction(nproma,n_zlev,patch%nblks_v)
    INTEGER :: icell_idx_1, icell_blk_1
    INTEGER :: icell_idx_2, icell_blk_2
    INTEGER :: boundary_counter
    INTEGER :: cell_index, cell_block
    INTEGER :: edge_index_cell, edge_block_cell
    INTEGER :: ictr, cell_edge,neigbor
    TYPE(t_cartesian_coordinates) :: cell1_cc, cell2_cc, vertex_cc

    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: owned_verts

 !   REAL(wp) :: comm_edges(nproma,n_zlev,patch%nblks_e)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:apply_boundary2coeffs')
    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')

    all_cells   => p_patch_3D%p_patch_2D(1)%cells%all
    all_edges   => p_patch_3D%p_patch_2D(1)%edges%all
    owned_verts => p_patch_3D%p_patch_2D(1)%verts%owned

    i_v_ctr(:,:,:)          = 0
    zarea_fraction(1:nproma,1:n_zlev,1:patch%nblks_v) = 0.0_wp

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

          !IF(v_base%lsm_e(je,jk,jb) /= sea) THEN
          IF ( p_patch_3D%lsm_e(je,jk,jb) /= sea ) THEN
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
!     DO jb = all_edges%start_block, all_edges%end_block
!       CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
!       DO jk = 1, n_zlev
!         DO je = i_startidx_e, i_endidx_e
!           ictr=0
!           DO neigbor=1,2
!             cell_index    = patch%edges%cell_idx(je, jb, neigbor)
!             cell_block    = patch%edges%cell_blk(je, jb, neigbor)   
!       
!             DO cell_edge=1,no_primal_edges!patch%cell_type
!               ictr=ictr+1
!               edge_index_cell = patch%cells%edge_idx(cell_index,cell_block,cell_edge)
!               edge_block_cell = patch%cells%edge_idx(cell_index,cell_block,cell_edge)
! 
!                IF ( p_patch_3D%lsm_e(edge_index_cell,jk,edge_block_cell) /= sea ) THEN
!                  ocean_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr)=0.0_wp
!                ENDIF
!             ENDDO
!           END DO
!         END DO
!       END DO
!     END DO
 
    !-------------------------------------------------------------
    !1) Set coefficients for div and grad to zero at boundary edges
    ! Also for ocean_coeff%edge2cell_coeff_cc
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jk=1,n_zlev
        DO jc = i_startidx_c, i_endidx_c

          DO je = 1, patch%cells%num_edges(jc,jb)

            ile = patch%cells%edge_idx(jc,jb,je)
            ibe = patch%cells%edge_blk(jc,jb,je)
            !IF ( v_base%lsm_e(ile,jk,ibe) /= sea) THEN
            IF ( p_patch_3D%lsm_e(ile,jk,ibe) /= sea ) THEN
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
          !i_edge_idx(1:4)         = 0
          !z_orientation(1:4)      = 0.0_wp
          boundary_counter        = 0

!           write(0,*) "-----------------------------"
          DO jev = 1, patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = patch%verts%edge_idx(jv,jb,jev)
            ibe = patch%verts%edge_blk(jv,jb,jev)

!             write(0,*) jb, jv, patch%verts%num_edges(jv,jb), jk, ":", ile, ibe, p_patch_3D%lsm_e(ile,jk,ibe)
            
            !Check, if edge is sea or boundary edge and take care of dummy edge
            ! edge with indices ile, ibe is sea edge
            ! edge with indices ile, ibe is boundary edge

            IF ( p_patch_3D%lsm_e(ile,jk,ibe) == SEA ) THEN
              i_v_ctr(jv,jk,jb) = i_v_ctr(jv,jk,jb)+1
            ELSEIF ( p_patch_3D%lsm_e(ile,jk,ibe) == BOUNDARY ) THEN

              !increase boundary edge counter
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
              !z_orientation(boundary_counter) = patch%verts%edge_orientation(jv,jb,jev)
              i_edge_idx(boundary_counter)    = jev

              ocean_coeff%bnd_edge_idx(jv,jk,jb,boundary_counter)= ile
              ocean_coeff%bnd_edge_blk(jv,jk,jb,boundary_counter)= ibe
              ocean_coeff%orientation(jv,jk,jb,boundary_counter) = &
                & patch%verts%edge_orientation(jv,jb,jev)
              ocean_coeff%edge_idx(jv,jk,jb,boundary_counter)    = jev

            END IF
          END DO ! jev = 1, patch%verts%num_edges(jv,jb)

          IF( MOD(boundary_counter,2) /= 0 ) THEN
            CALL finish (routine,'MOD(boundary_counter,2) /= 0 !!')
          ENDIF
!---------------------------------------------------------------------------------
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

              cell1_cc%x  = patch%cells%cartesian_center(icell_idx_1,icell_blk_1)%x
              cell2_cc%x  = patch%cells%cartesian_center(icell_idx_2,icell_blk_2)%x

              !Check, if edge is sea or boundary edge and take care of dummy edge
              !edge with indices ile, ibe is sea edge
              !Add up for wet dual area.
              !IF ( v_base%lsm_e(ile,jk,ibe) <= sea_boundary ) THEN
              IF ( p_patch_3D%lsm_e(ile,jk,ibe) <= sea_boundary ) THEN
                ocean_coeff%variable_dual_vol_norm(jv,jk,jb,jev)= triangle_area(cell1_cc, vertex_cc, cell2_cc)
                ! edge with indices ile, ibe is boundary edge
              ELSE IF ( p_patch_3D%lsm_e(ile,jk,ibe) == boundary ) THEN
                ocean_coeff%variable_dual_vol_norm(jv,jk,jb,jev)=0.0_wp!0.5_wp*triangle_area(cell1_cc, vertex_cc, cell2_cc)
              END IF
            END DO

         !---------------------------------------------------------------------------------------------
          DO je = 1, boundary_counter
            ibnd_edge_idx(je) = ocean_coeff%bnd_edge_idx(jv,jk,jb,je)
            ibnd_edge_blk(je) = ocean_coeff%bnd_edge_blk(jv,jk,jb,je)

            ocean_coeff%rot_coeff(jv,jk,jb,i_edge_idx(je) )=&
              & 0.5_wp*patch%edges%system_orientation(ibnd_edge_idx(je),ibnd_edge_blk(je)) * &
              & patch%edges%primal_edge_length(ibnd_edge_idx(je),ibnd_edge_blk(je))

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

        DO je = 1, patch%cells%num_edges(jc,jb)
          ocean_coeff%edge2cell_coeff_cc_dyn(jc,1,jb,je)%x = &
            ocean_coeff%edge2cell_coeff_cc(jc,1,jb,je)%x 
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

          DO je = 1, patch%verts%num_edges(jv,jb)
            ile = patch%verts%edge_idx(jv,jb,je)
            ibe = patch%verts%edge_blk(jv,jb,je)

             IF ( p_patch_3D%lsm_e(ile,jk,ibe) /= sea) THEN
              ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,je)%x(1:3) = 0.0_wp
              ocean_coeff%variable_dual_vol_norm(jv,jk,jb,je)    = 0.0_wp
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

          IF ( i_v_ctr(jv,jk,jb) == patch%verts%num_edges(jv,jb) ) THEN
            zarea_fraction(jv,jk,jb)= patch%verts%dual_area(jv,jb)/(earth_radius*earth_radius)
            !zarea_fraction(jv,jk,jb)=SUM(ocean_coeff%variable_dual_vol_norm(jv,jk,jb,:))

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

              cell1_cc%x  = patch%cells%cartesian_center(icell_idx_1,icell_blk_1)%x
              cell2_cc%x  = patch%cells%cartesian_center(icell_idx_2,icell_blk_2)%x

              !Check, if edge is sea or boundary edge and take care of dummy edge
              !edge with indices ile, ibe is sea edge
              !Add up for wet dual area.
              ! Note that this should be modified.
              !   sea_boundary means an open boundary
              !   boundary means that only the sea cell are should be added
              IF ( p_patch_3D%lsm_e(ile,jk,ibe) <= sea_boundary ) THEN
                zarea_fraction(jv,jk,jb) = zarea_fraction(jv,jk,jb)  &
                  & + triangle_area(cell1_cc, vertex_cc, cell2_cc)
                ! edge with indices ile, ibe is boundary edge                
              ELSE IF ( p_patch_3D%lsm_e(ile,jk,ibe) == boundary ) THEN
!                 zarea_fraction(jv,jk,jb) = zarea_fraction(jv,jk,jb)  &
!                   & + 0.5_wp*triangle_area(cell1_cc, vertex_cc, cell2_cc)
                ! add only the sea dual area, ie the trinalge are between
                ! the vertex, edge centre, and sea cell center
                ! use cell1_cc%x for the center of the sea cell
                IF (p_patch_3D%lsm_c(icell_idx_1,jk,icell_blk_1) <= sea_boundary) THEN
                  ! just a consistence check
                  IF (p_patch_3D%lsm_c(icell_idx_2,jk,icell_blk_2) <= sea_boundary) &
                     CALL finish(routine, "Incosistent cell and edge lsm")
                ELSE
                  ! the idx_2 cell is the sea cell
                  cell1_cc%x  = patch%cells%cartesian_center(icell_idx_2,icell_blk_2)%x
                  IF (p_patch_3D%lsm_c(icell_idx_1,jk,icell_blk_1) <= sea_boundary) &
                     CALL finish(routine, "Incosistent cell and edge lsm")
                ENDIF  
                  
                zarea_fraction(jv,jk,jb) = zarea_fraction(jv,jk,jb)  &
                  & + triangle_area(cell1_cc, vertex_cc, patch%edges%cartesian_center(ile,ibe))
              
              END IF
            END DO
          ENDIF !( i_v_ctr(jv,jk,jb) == patch%verts%num_edges(jv,jb) )

          !Final coefficient calculation
!CDIR nextscalar
          IF(zarea_fraction(jv,jk,jb)/=0.0_wp)THEN

            ocean_coeff%rot_coeff(jv,jk,jb,:)&
            &=ocean_coeff%rot_coeff(jv,jk,jb,:)/(zarea_fraction(jv,jk,jb)*(earth_radius*earth_radius))

            DO jev = 1, patch%verts%num_edges(jv,jb)
              ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)&
                & =ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)/zarea_fraction(jv,jk,jb)
            END DO

          ELSE
            DO jev = 1, patch%verts%num_edges(jv,jb)
                ocean_coeff%edge2vert_coeff_cc(jv,jk,jb,jev)%x(1:3)=0.0_wp
            END DO
            ocean_coeff%rot_coeff(jv,jk,jb,:)=0.0_wp
          ENDIF
         !!ENDIF !( i_v_ctr(jv,jk,jb) == patch%verts%num_edges(jv,jb) )

        ENDDO!jv = i_startidx_v, i_endidx_v
      END DO!jk = 1, n_zlev
    END DO!jb = owned_verts%start_block, owned_verts%end_block
    ! sync the result
    DO jev=1,6
      DO jk = 1, n_zlev
        CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,jk,:, jev)%x(1))
        CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,jk,:, jev)%x(2))
        CALL sync_patch_array(SYNC_V, patch, ocean_coeff%edge2vert_coeff_cc(:,jk,:, jev)%x(3))
      ENDDO
    ENDDO
    DO je=1,patch%cell_type
      CALL sync_patch_array(SYNC_V, patch, ocean_coeff%rot_coeff(:,:,:, je))
    ENDDO


! !     DO jk = 1, n_zlev
! !       CALL sync_patch_array(SYNC_V, patch, zarea_fraction(:,jk,:))
! !     ENDDO
! !     DO jev=1,2*no_dual_edges
! !       DO jk = 1, n_zlev
! !         CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2edge_viavert_coeff(:,jk,:, jev))
! !       ENDDO
! !     ENDDO
! !     DO jb = all_edges%start_block, all_edges%end_block
! !       CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
! !       DO jk = 1, n_zlev
! !         DO je = i_startidx_e, i_endidx_e
! ! 
! !           DO neigbor=1,2
! !             jv  = patch%edges%vertex_idx(je, jb, neigbor)
! !             jev = patch%edges%vertex_blk(je, jb, neigbor)   
! ! 
! !             IF(zarea_fraction(jv,jk,jev)/=0.0_wp)THEN
! ! 
! !               IF(neigbor==1)THEN
! !                 ocean_coeff%edge2edge_viavert_coeff(je,jk,jb,1:6)&
! !                 &=ocean_coeff%edge2edge_viavert_coeff(je,jk,jb,1:6)/zarea_fraction(jv,jk,jev)
! !               ELSEIF(neigbor==2)THEN
! !                 ocean_coeff%edge2edge_viavert_coeff(je,jk,jb,7:12)&
! !                 &=ocean_coeff%edge2edge_viavert_coeff(je,jk,jb,7:12)/zarea_fraction(jv,jk,jev)
! !               ENDIF
! !             ENDIF
! !           END DO
! !         END DO
! !       END DO
! !     END DO
! ! 
! !     DO jev=1,2*no_dual_edges
! !       DO jk = 1, n_zlev
! !         CALL sync_patch_array(SYNC_E, patch, ocean_coeff%edge2edge_viavert_coeff(:,jk,:, jev))
! !       ENDDO
! !     ENDDO

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
  SUBROUTINE init_geo_factors_oce_3d( patch, p_coeff )
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: patch
    TYPE(t_operator_coeff),INTENT(inout) :: p_coeff
    !

    INTEGER :: jc, jb, je, jv, ie,jk
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER :: ile, ibe!, ilc1, ibc1, ilc2, ibc2, ifac, ic, ilnc, ibnc
    !INTEGER :: ile1, ibe1,ile2,ibe2,ile3,ibe3
    !INTEGER, PARAMETER :: i_cell_type = 3
    !TYPE(cartesian_coordinates)::z_pn_k,z_pn_j
    !REAL(wp) :: z_lon, z_lat, z_nu, z_nv, z_proj
    !REAL(wp) :: cell_area
    TYPE(t_subset_range), POINTER :: owned_edges         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_cells         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_verts         ! these are the owned entities
    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor, level
    

    REAL(wp) :: z_sync_c(nproma,n_zlev,patch%nblks_c)
    !REAL(wp) :: z_sync_e(nproma,n_zlev,patch%nblks_e)
    REAL(wp) :: z_sync_v(nproma,n_zlev,patch%nblks_v)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:init_geo_factors_oce_3d')
    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')
    owned_edges => patch%edges%owned
    owned_cells => patch%cells%owned
    owned_verts => patch%verts%owned
    
    i_nchdom   = MAX(1,patch%n_childdom)

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

!     ! 1) coefficients for divergence
!     rl_start = 1
!     rl_end = min_rlcell
! 
!     ! values for the blocking
!     i_startblk = patch%cells%start_blk(rl_start,1)
!     i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)
!     !
!     ! loop through all patch cells 
!     !
! !$OMP DO PRIVATE(jb,je,jc,i_startidx,i_endidx,ile,ibe)
!     DO jk=1,n_zlev
!       DO jb = i_startblk, i_endblk
!         CALL get_indices_c(patch, jb, i_startblk, i_endblk,      &
!           & i_startidx, i_endidx, rl_start, rl_end)
!         DO jc = i_startidx, i_endidx
!           !         ile1 = patch%cells%edge_idx(jc,jb,1)
!           !         ibe1 = patch%cells%edge_blk(jc,jb,1)
!           !         ile2 = patch%cells%edge_idx(jc,jb,2)
!           !         ibe2 = patch%cells%edge_blk(jc,jb,2)
!           !         ile3 = patch%cells%edge_idx(jc,jb,3)
!           !         ibe3 = patch%cells%edge_blk(jc,jb,3)
!           !         cell_area =  0.25_wp&
!           !         & *( patch%edges%primal_edge_length(ile1,ibe1)*patch%edges%dual_edge_length(ile1,ibe1)&
!           !         &   +patch%edges%primal_edge_length(ile2,ibe2)*patch%edges%dual_edge_length(ile2,ibe2)&
!           !         &   +patch%edges%primal_edge_length(ile3,ibe3)*patch%edges%dual_edge_length(ile3,ibe3))
! 
!           DO je = 1, no_primal_edges
! 
!             ile = patch%cells%edge_idx(jc,jb,je)
!             ibe = patch%cells%edge_blk(jc,jb,je)
! 
!             p_coeff%div_coeff(jc,jk,jb,je) =                &
!               & patch%edges%primal_edge_length(ile,ibe) * &
!               & patch%cells%edge_orientation(jc,jb,je)  / &
!               & patch%cells%area(jc,jb)
!           ENDDO !edge loop
!           !write(1234,*)'div coeff 3D',jk,jc,jb,p_coeff%div_coeff(jc,jk,jb,:)
!         ENDDO !idx loop
!       END DO !block loop
!     END DO
! !$OMP END DO


!$OMP DO PRIVATE(cell_block, start_index, end_index, cell_index, neigbor, &
!$OMP edge_index, edge_block, level) ICON_OMP_DEFAULT_SCHEDULE
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        DO neigbor=1, patch%cells%num_edges(cell_index, cell_block)  ! no_dual_edges

          edge_index = patch%cells%edge_idx(cell_index, cell_block, neigbor)
          edge_block = patch%cells%edge_blk(cell_index, cell_block, neigbor)

          p_coeff%div_coeff(cell_index, 1, cell_block, neigbor) =                &
            & patch%edges%primal_edge_length(edge_index, edge_block) * &
            & patch%cells%edge_orientation(cell_index, cell_block, neigbor)  / &
            & patch%cells%area(cell_index, cell_block)

          DO level=2, n_zlev
            p_coeff%div_coeff(cell_index, level, cell_block, neigbor) = &
               p_coeff%div_coeff(cell_index, 1, cell_block, neigbor)      
          ENDDO !levels
        
        ENDDO !neigbor
        
      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block 
!$OMP END DO NOWAIT


!     ! 2) coefficients for curl
!     rl_start = 1  ! #slo# changed to 1 - 2010-12-07
!     rl_end = min_rlvert_int
! 
!     ! values for the blocking
!     i_startblk = patch%verts%start_blk(rl_start,1)
!     i_endblk   = patch%verts%end_blk(rl_end,i_nchdom)
!     !
!     ! loop through all patch cells (and blocks)
!     !
! !$OMP DO PRIVATE(jb,je,jv,i_startidx,i_endidx,ile,ibe)
!     DO jk=1,n_zlev
!       DO jb = i_startblk, i_endblk
! 
!         CALL get_indices_v(patch, jb, i_startblk, i_endblk, &
!           & i_startidx, i_endidx, rl_start, rl_end)
! 
!         DO jv = i_startidx, i_endidx
!           DO je = 1, no_dual_edges!patch%verts%num_edges(jv,jb)
! 
!             ile = patch%verts%edge_idx(jv,jb,je)
!             ibe = patch%verts%edge_blk(jv,jb,je)
! 
!             p_coeff%rot_coeff(jv,jk,jb,je) =              &
!               & patch%edges%dual_edge_length(ile,ibe) * &
!               & patch%verts%edge_orientation(jv,jb,je)
!           ENDDO
!         ENDDO
!       END DO
!     END DO
! !$OMP END DO

   
!$OMP DO PRIVATE(vertex_block, start_index, end_index, vertex_index, neigbor, &
!$OMP edge_index, edge_block, level) ICON_OMP_DEFAULT_SCHEDULE
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index

        DO neigbor=1, patch%verts%num_edges(vertex_index, vertex_block)  ! no_dual_edges

          edge_index = patch%verts%edge_idx(vertex_index, vertex_block, neigbor)
          edge_block = patch%verts%edge_blk(vertex_index, vertex_block, neigbor)

          p_coeff%rot_coeff(vertex_index,1,vertex_block,neigbor)  =     &
            & patch%edges%dual_edge_length(edge_index,edge_block) *     &
            & patch%verts%edge_orientation(vertex_index,vertex_block,neigbor)
          DO level=2, n_zlev
            p_coeff%rot_coeff(vertex_index,level,vertex_block,neigbor) = &
               p_coeff%rot_coeff(vertex_index,1,vertex_block,neigbor)      
          ENDDO !levels

        ENDDO !neigbor
      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
!$OMP END DO NOWAIT



    ! 4) coefficients for gradient

    rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    rl_end = min_rledge_int

    ! Second step: computed projected orientation vectors and related information
    i_startblk = patch%edges%start_blk(rl_start,1)
    i_endblk   = patch%edges%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch edges
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je)
    DO jk=1,n_zlev
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO je =  i_startidx, i_endidx
          ! compute inverse dual edge length (undefined for refin_ctrl=1)
          patch%edges%inv_dual_edge_length(je,jb) = &
            & 1._wp/patch%edges%dual_edge_length(je,jb)

          p_coeff%grad_coeff(je,jk,jb)&
            & = patch%edges%inv_dual_edge_length(je,jb)

        ENDDO
      END DO !block loop
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! synchronize all elements of p_coeff:
    CALL sync_patch_array(sync_e, patch, p_coeff%grad_coeff)

    DO ie = 1, no_primal_edges

      z_sync_c(:,:,:) = p_coeff%div_coeff(:,:,:,ie)
      CALL sync_patch_array(sync_c, patch, z_sync_c(:,:,:))
      p_coeff%div_coeff(:,:,:,ie) = z_sync_c(:,:,:)
      !         z_sync_c(:,:,:) = p_coeff%n2s_coeff(:,:,:,ie)
      !         CALL sync_patch_array(SYNC_C, patch, z_sync_c(:,:,:))
      !         p_coeff%n2s_coeff(:,:,:,ie) = z_sync_c(:,:,:)
    END DO

    CALL sync_patch_array(sync_e, patch, patch%edges%inv_dual_edge_length(:,:))

    DO ie = 1, no_dual_edges!9-i_cell_type
      z_sync_v(:,:,:) = p_coeff%rot_coeff(:,:,:,ie)
      CALL sync_patch_array(sync_v, patch, z_sync_v(:,:,:))
      p_coeff%rot_coeff(:,:,:,ie) = z_sync_v(:,:,:)
    END DO

    CALL message (TRIM(routine), 'end')

    ! 3) coefficients for nabla2_scalar
    !rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    !rl_end = min_rlcell_int
    ! values for the blocking
    !i_startblk = patch%cells%start_blk(rl_start,1)
    !i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
    !!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
    !!$OMP    ilc2,ibc2,ilnc,ibnc)
    ! !       DO jk=1,n_zlev
    ! !       DO jb = i_startblk, i_endblk
    ! !
    ! !         CALL get_indices_c(patch, jb, i_startblk, i_endblk, &
    ! !                            i_startidx, i_endidx, rl_start, rl_end)
    ! !
    ! !         DO je = 1, i_cell_type
    ! !           DO jc = i_startidx, i_endidx
    ! !
    ! !             ile = patch%cells%edge_idx(jc,jb,je)
    ! !             ibe = patch%cells%edge_blk(jc,jb,je)
    ! !
    ! !             ilc1 = patch%edges%cell_idx(ile,ibe,1)
    ! !             ibc1 = patch%edges%cell_blk(ile,ibe,1)
    ! !             ilc2 = patch%edges%cell_idx(ile,ibe,2)
    ! !             ibc2 = patch%edges%cell_blk(ile,ibe,2)
    ! !
    ! !             IF (jc == ilc1 .AND. jb == ibc1) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 p_coeff%geofac_n2s(jc,1,jb)     =  &
    ! !                 &  p_coeff%geofac_n2s(jc,1,jb)  -  &
    ! !                 &  p_coeff%geofac_div(jc,je,jb) /  &
    ! !                 &  patch%edges%dual_edge_length(ile,ibe)
    ! !             ENDIF
    ! !           ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
    ! !             IF (i_cell_type == 3) THEN
    ! !               p_coeff%geofac_n2s(jc,1,jb)       =  &
    ! !                 &  p_coeff%geofac_n2s(jc,1,jb)  +  &
    ! !                 &  p_coeff%geofac_div(jc,je,jb) /  &
    ! !                 &  patch%edges%dual_edge_length(ile,ibe)
    ! !             ENDIF
    ! !           ENDIF
    ! !           DO ic = 1, i_cell_type
    ! !             ilnc = patch%cells%neighbor_idx(jc,jb,ic)
    ! !             ibnc = patch%cells%neighbor_blk(jc,jb,ic)
    ! !             IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 p_coeff%geofac_n2s(jc,ic+1,jb)     = &
    ! !                   &  p_coeff%geofac_n2s(jc,ic+1,jb)- &
    ! !                   &  p_coeff%geofac_div(jc,je,jb)  / &
    ! !                   &  patch%edges%dual_edge_length(ile,ibe)
    ! !               ENDIF
    ! !             ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
    ! !               IF (i_cell_type == 3) THEN
    ! !                 p_coeff%geofac_n2s(jc,ic+1,jb)     = &
    ! !                   &  p_coeff%geofac_n2s(jc,ic+1,jb)+ &
    ! !                   &  p_coeff%geofac_div(jc,je,jb)  / &
    ! !                   &  patch%edges%dual_edge_length(ile,ibe)
    ! !               ENDIF
    ! !             ENDIF
    ! !           ENDDO
    ! !           ! To ensure that dummy edges have a factor of 0:
    ! !           IF (je > patch%cells%num_edges(jc,jb)) THEN
    ! !             p_coeff%geofac_n2s(jc,je+1,jb) = 0._wp
    ! !           ENDIF
    ! !
    ! !         ENDDO !cell loop
    ! !       ENDDO
    ! !     END DO !block loop
    ! !     END DO
    !!$OMP END DO

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

