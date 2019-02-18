!>
!! Contains the implementation of velocity and tracer diffusion for the ICON ocean model.
!!
!!
!! @par Revision History
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_tracer_diffusion
  
  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: boundary, sea_boundary, min_dolic ! ,max_char_length
  USE mo_parallel_config,     ONLY: nproma
  USE mo_ocean_nml,           ONLY: n_zlev, iswm_oce, VelocityDiffusion_order, laplacian_form
  USE mo_run_config,          ONLY: dtime
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_ocean_types,         ONLY: t_hydro_ocean_diag, t_hydro_ocean_aux
  USE mo_ocean_tracer_transport_types,  ONLY: t_ocean_tracer
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_ocean_physics_types, ONLY: t_ho_params
  USE mo_scalar_product,      ONLY: map_cell2edges_3d, map_edges2edges_viacell_3d_const_z
  USE mo_ocean_math_operators,ONLY: div_oce_3d, rot_vertex_ocean_3d,&
    & map_edges2vert_3d, grad_fd_norm_oce_3D, grad_vector, div_vector_onTriangle
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, sync_patch_array_mult
  USE mo_exception,           ONLY: finish !, message_text, message
 
  IMPLICIT NONE
  
  PRIVATE
  
  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(LEN=12)           :: str_module    = 'oceDiffusion'  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  !
  ! PUBLIC INTERFACE
  !
  INTEGER, PARAMETER :: top=1
  PUBLIC :: tracer_diffusion_horz
  PUBLIC :: tracer_diffusion_vert_explicit
  PUBLIC :: tracer_diffusion_vertical_implicit
  
CONTAINS
    
  !-------------------------------------------------------------------------
  !Subroutine computes the horizontal diffusive flux of an arbitrary tracer.
   SUBROUTINE tracer_diffusion_horz(patch_3D, trac_in, diff_flx, k_t, subset_range)
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)              :: trac_in(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)           :: diff_flx(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp), OPTIONAL                :: k_t(:,:,:) !mixing coefficient for tracer    
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    !
    !Local variables
    INTEGER :: level, blockNo, edge_index
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_ocediffusion:tracer_diffusion_horz')
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-------------------------------------------------------------------------------
    
  IF(PRESENT(k_t))THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      diff_flx(:,:,blockNo) = 0.0_wp
      DO edge_index = start_edge_index, end_edge_index
        !Get indices of two adjacent triangles
        il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
        ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
        il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
        ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)

        DO level=1,  patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
          
          diff_flx(edge_index,level,blockNo) = &
            &   k_t(edge_index,level,blockNo) &
            & * patch_3D%p_patch_1d(1)%prism_thick_e(edge_index,level,blockNo)  &
            & * (trac_in(il_c2,level,ib_c2) - trac_in(il_c1,level,ib_c1))       &
            & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
            
        ENDDO
        
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2D, diff_flx)
    ENDIF    
    
  ELSEIF(.NOT.PRESENT(k_t))THEN  
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      diff_flx(:,:,blockNo) = 0.0_wp
      DO edge_index = start_edge_index, end_edge_index
        !Get indices of two adjacent triangles
        il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
        ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
        il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
        ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)

        DO level=1,  patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
          
          diff_flx(edge_index,level,blockNo) = &
            & patch_3D%p_patch_1d(1)%prism_thick_e(edge_index,level,blockNo)  &
            & * (trac_in(il_c2,level,ib_c2) - trac_in(il_c1,level,ib_c1))       &
            & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
            
        ENDDO
        
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2D, diff_flx)
    ENDIF 
  ENDIF  
    
  END SUBROUTINE tracer_diffusion_horz
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!Subroutine computes the vertical diffusive flux of an arbitrary tracer.
  !>
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  SUBROUTINE tracer_diffusion_vert_explicit(patch_3D,        &
    & trac_c,          &
    & top_bc_tracer,   &
    & a_v,             &
    & div_diff_flx)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3D
    REAL(wp), INTENT(in)              :: trac_c       (nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)              :: top_bc_tracer(nproma, patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)              :: a_v(:,:,:)
    REAL(wp), INTENT(inout)             :: div_diff_flx(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    !
    !Local variables
    INTEGER :: start_level, end_level
    INTEGER :: cell_index, level, blockNo
    INTEGER :: start_index, end_index
    INTEGER :: z_dolic
    ! vertical diffusive tracer flux
    REAL(wp)                      :: z_diff_flx(nproma, n_zlev+1,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_ocean_tracer_diffusion:tracer_diffusion_vert')
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2d(1)
    all_cells => patch_2D%cells%ALL
    !-----------------------------------------------------------------------
    start_level              = 1
    end_level              = n_zlev
    z_diff_flx(1:nproma, 1:n_zlev+1,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks) = 0.0_wp
    
    !First vertical flux then vertical divergence
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO cell_index = start_index, end_index
        z_dolic  = patch_3D%p_patch_1d(1)%dolic_c(cell_index,blockNo)!v_base%dolic_c(cell_index,blockNo)
        IF ( z_dolic >=min_dolic ) THEN
          !1a) 0cean surface
          z_diff_flx(cell_index,start_level,blockNo) =a_v(cell_index,start_level,blockNo)*top_bc_tracer(cell_index,blockNo)
          
          !1b) ocean interior
          DO level = start_level+1, z_dolic

            z_diff_flx(cell_index,level,blockNo)&
              & = a_v(cell_index,level,blockNo) &
              & * (trac_c(cell_index,level-1,blockNo)-trac_c(cell_index,level,blockNo))/ patch_3D%p_patch_1d(1)%del_zlev_i(level)
          END DO
          
          
          DO level = 1, z_dolic
            ! positive vertical divergence in direction of w (upward positive)
            div_diff_flx(cell_index,level,blockNo) = &
              & (z_diff_flx(cell_index,level,blockNo) - z_diff_flx(cell_index,level+1,blockNo))&
              & /patch_3D%p_patch_1d(1)%del_zlev_m(level)
          END DO
          !1c) ocean bottom zero bottom boundary condition
          !diff_flx(cell_index,z_dolic+1,blockNo) = bot_bc_tracer(cell_index,level)
        ENDIF
      END DO
    END DO
    
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=5  ! output print level (1-5, fix)
    CALL dbg_print('TrcDiffExpl: z_diff_flx'   ,z_diff_flx               ,str_module,idt_src)
    CALL dbg_print('TrcDiffExpl: div_diff_flx' ,div_diff_flx             ,str_module,idt_src)
    !---------------------------------------------------------------------
    
  END SUBROUTINE tracer_diffusion_vert_explicit
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE tracer_diffusion_vertical_implicit( &
    & patch_3d,                  &
    & ocean_tracer,              &
    & a_v,                       &
    & operators_coefficients) !,  &
    ! & diff_column)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: ocean_tracer
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)
    TYPE(t_operator_coeff),TARGET :: operators_coefficients
    !
    INTEGER :: cell_block, start_index, end_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d

    !-----------------------------------------------------------------------
    cells_in_domain       =>  patch_3d%p_patch_2d(1)%cells%in_domain
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO cell_block = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, cell_block, start_index, end_index)

      CALL tracer_diffusion_vertical_implicit_onBlock( &
        & patch_3d,                  &
        & ocean_tracer,              &
        & a_v,                       &
        & operators_coefficients,    &
        & cell_block, start_index, end_index)

    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE tracer_diffusion_vertical_implicit
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for scalar fields.
  !>
  !! sbr identical to sbr above but now with homogeneous boundary conditions
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !! Preconditioning by Leonidas Linardakis, MPI-M (2014)
  !!
  !! The result ocean_tracer%concetration is calculated on domain_cells
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE tracer_diffusion_vertical_implicit_onBlock( &
    & patch_3d,                &
    & ocean_tracer,            &
    & a_v,                     &
    & operators_coefficients,  &
    & blockNo, start_index, end_index) !,  &
    ! & diff_column)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: ocean_tracer
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)
    TYPE(t_operator_coeff),TARGET :: operators_coefficients
    INTEGER, INTENT(in) :: blockNo, start_index, end_index
    !
    !
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)! , nb(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    REAL(wp) :: column_tracer(1:n_zlev)
    REAL(wp) :: dt_inv, diagonal_product
    REAL(wp), POINTER :: field_column(:,:,:)
    INTEGER :: bottom_level
    INTEGER :: cell_index, level
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d

!     REAL(wp) :: tmp, tmp_add
    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2d%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    dt_inv = 1.0_wp/dtime
    
    DO cell_index = start_index, end_index
      bottom_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)

      IF (bottom_level < 2 ) CYCLE ! nothing to diffuse

      DO level=1,bottom_level
        inv_prism_thickness(level)        = patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_index,level,blockNo)
        inv_prisms_center_distance(level) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(cell_index,level,blockNo)
      ENDDO

      !------------------------------------
      ! Fill triangular matrix
      ! b is diagonal, a is the upper diagonal, c is the lower
      !   top level
      a(1) = 0.0_wp
      c(1) = -a_v(cell_index,2,blockNo) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
      b(1) = dt_inv - c(1)
      DO level = 2, bottom_level-1
        a(level) = - a_v(cell_index,level,blockNo)   * inv_prism_thickness(level) * inv_prisms_center_distance(level)
        c(level) = - a_v(cell_index,level+1,blockNo) * inv_prism_thickness(level) * inv_prisms_center_distance(level+1)
        b(level) = dt_inv - a(level) - c(level)
      END DO
      ! bottom
      a(bottom_level) = -a_v(cell_index,bottom_level,blockNo) * &
        & inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)
      b(bottom_level) = dt_inv - a(bottom_level)

      ! precondition: set diagonal equal to diagonal_product
      diagonal_product = PRODUCT(b(1:bottom_level))
!       tmp_add = 0.0_wp

      DO level = 1, bottom_level
        fact(level) = diagonal_product / b(level)
        a(level)  = a(level)  * fact(level)
        b(level)  = diagonal_product
        c(level)  = dt_inv * fact(level) - a(level) - b(level)
 
!         tmp = field_column(cell_index,level,blockNo)
!         column_tracer(level) = tmp
!         tmp_add = tmp_add + column_tracer(level)
!         tmp = dt_inv
!         column_tracer(level) = tmp
!         tmp_add = tmp_add + column_tracer(level)
!         tmp = diagonal_product
!         column_tracer(level) = tmp
!         tmp_add = tmp_add + column_tracer(level)
!         tmp = b(level)
!         column_tracer(level) = tmp
!         tmp_add = tmp_add + column_tracer(level)
!         tmp = fact(level)
!         column_tracer(level) = tmp
!         tmp_add = tmp_add + column_tracer(level)

        column_tracer(level) = field_column(cell_index,level,blockNo) * dt_inv * fact(level)

      ENDDO
      c(bottom_level) = 0.0_wp
!       write(0,*) tmp_add

      !------------------------------------
      ! solver from lapack
      !
      ! eliminate lower diagonal
      DO level=bottom_level-1, 1, -1
        fact(level+1)  = c( level ) / b( level+1 )
        b( level ) = b( level ) - fact(level+1) * a( level +1 )
        column_tracer( level ) = column_tracer( level ) - fact(level+1) * column_tracer( level+1 )
      ENDDO

      !     Back solve with the matrix U from the factorization.
      column_tracer( 1 ) = column_tracer( 1 ) / b( 1 )
      DO level =  2, bottom_level
        column_tracer( level ) = ( column_tracer( level ) - a( level ) * column_tracer( level-1 ) ) / b( level )
      ENDDO

      DO level = 1, bottom_level
        ocean_tracer%concentration(cell_index,level,blockNo) = column_tracer(level)
      ENDDO
    
    ENDDO ! cell_index
    
  END SUBROUTINE tracer_diffusion_vertical_implicit_onBlock
  !------------------------------------------------------------------------
    
END MODULE mo_ocean_tracer_diffusion
