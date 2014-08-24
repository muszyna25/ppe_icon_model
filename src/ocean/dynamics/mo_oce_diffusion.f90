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
!! This code is subedge_indexct to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_oce_diffusion
  
  USE mo_kind,                ONLY: wp
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates, gvec2cvec !, gc2cc
  USE mo_impl_constants,      ONLY: boundary, sea_boundary, min_dolic ! ,max_char_length
  USE mo_parallel_config,     ONLY: nproma
  USE mo_ocean_nml,           ONLY: n_zlev, iswm_oce, veloc_diffusion_order, veloc_diffusion_form
  USE mo_run_config,          ONLY: dtime
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_oce_types,           ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, t_ocean_tracer, t_hydro_ocean_aux
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_oce_physics,         ONLY: t_ho_params
  USE mo_scalar_product,      ONLY: map_cell2edges_3d
  USE mo_oce_math_operators,  ONLY: div_oce_3d, rot_vertex_ocean_3d,&
    & map_edges2vert_3d
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
  PUBLIC :: velocity_diffusion
  PUBLIC :: veloc_diff_harmonic_div_grad, veloc_diff_biharmonic_div_grad
  PUBLIC :: veloc_diff_harmonic_curl_curl, veloc_diff_biharmonic_curl_curl
  ! PUBLIC :: velocity_diffusion_vert_explicit
  PUBLIC :: velocity_diffusion_vertical_implicit
  PUBLIC :: velocity_diffusion_vertical_implicit_onBlock
  PUBLIC :: tracer_diffusion_horz
  PUBLIC :: tracer_diffusion_vert_explicit
  PUBLIC :: tracer_diffusion_vertical_implicit
  
CONTAINS
  
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! !  SUBROUTINE calculates horizontal diffusion of edge velocity via laplacian diffusion
  !!    implemented as P^T div( K_H grad P v).
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  SUBROUTINE velocity_diffusion( patch_3D, vn_in, p_param, p_diag,p_op_coeff, laplacian_vn_out)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp)                               :: vn_in(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    TYPE(t_ho_params), INTENT(in)       :: p_param !mixing parameters
    TYPE(t_hydro_ocean_diag),INTENT(in) :: p_diag
    TYPE(t_operator_coeff),INTENT(in)   :: p_op_coeff
    REAL(wp), INTENT(inout)             :: laplacian_vn_out(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    
    !Local variables
    !REAL(wp) :: z_lapl(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    !INTEGER  :: level
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_oce_diffusion:velocity_diffusion_horz')
    !-------------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    
    IF(veloc_diffusion_order==1)THEN
      
      !divgrad laplacian is chosen
      IF(veloc_diffusion_form==2)THEN
        
        CALL veloc_diff_harmonic_div_grad( patch_3D,&
          & p_param,   &
          & p_diag,    &
          & p_op_coeff,&
          & laplacian_vn_out)
        
      ELSEIF(veloc_diffusion_form==1)THEN
        ! inUse
        CALL veloc_diff_harmonic_curl_curl( vn_in, p_diag%vort,&
          & patch_3D,           &
          & p_op_coeff,        &
          & p_diag%p_vn_dual,  &
          & laplacian_vn_out,  &
          & p_param%k_veloc_h)
        
      ENDIF
      
    ELSEIF(veloc_diffusion_order==2)THEN
      
      IF(veloc_diffusion_form==2)THEN
        
        CALL veloc_diff_biharmonic_div_grad( patch_3D,   &
          & p_param,      &
          & p_diag,       &
          & p_op_coeff,   &
          & laplacian_vn_out)
      ELSEIF(veloc_diffusion_form==1)THEN
        
        CALL veloc_diff_biharmonic_curl_curl( patch_3D,        &
          & vn_in,             &
          & p_diag%vort,       &
          & p_param%k_veloc_h, &
          & p_op_coeff,        &
          & p_diag%p_vn_dual,  &
          & laplacian_vn_out)
      ENDIF
    ENDIF
    
    ! CALL sync_patch_array(sync_e, patch_3D%p_patch_2d(1), laplacian_vn_out)
    
  END SUBROUTINE velocity_diffusion
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates horizontal diffusion of edge velocity via laplacian diffusion
  !!    implemented as P^T div( K_H grad P v).
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  SUBROUTINE veloc_diff_harmonic_div_grad( patch_3D, p_param, p_diag,&
    & p_op_coeff, laplacian_vn_out)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    TYPE(t_ho_params), INTENT(in)     :: p_param !mixing parameters
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)           :: laplacian_vn_out(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    
    !Local variables
    INTEGER :: start_level, elev
    INTEGER :: level, blockNo, edge_index,cell_index
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: start_index, end_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: idx_cartesian
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_cartesian_coordinates) :: z_grad_u    (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: z_div_grad_u(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells, all_edges, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-------------------------------------------------------------------------------
    patch_2D         => patch_3D%p_patch_2d(1)
    all_cells       => patch_2D%cells%ALL
    all_edges       => patch_2D%edges%ALL
    edges_in_domain => patch_2D%edges%in_domain
    !-------------------------------------------------------------------------------
    start_level = 1
    elev = n_zlev
    
    laplacian_vn_out(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e) = 0.0_wp
    
    ! loop over cells in local domain + halo
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO level = start_level, elev
        DO cell_index = start_index, end_index
          z_div_grad_u(cell_index,level,blockNo)%x =  0.0_wp
        END DO
      END DO
    END DO
    
    ! loop over edges in local domain + halo
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
      DO level = start_level, elev
        DO edge_index = start_edge_index, end_edge_index
          z_grad_u(edge_index,level,blockNo)%x = 0.0_wp
        ENDDO
      END DO
    END DO
    
    !-------------------------------------------------------------------------------------------------------
    !Step 1: Calculate gradient of cell velocity vector.
    !Result is a gradient vector, located at edges
    !Step 2: Multiply each component of gradient vector with mixing coefficients
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      DO level = start_level, elev
        DO edge_index = start_edge_index, end_edge_index
          !IF ( v_base%lsm_e(edge_index,level,blockNo) <= sea_boundary ) THEN
          IF (patch_3D%lsm_e(edge_index,level,blockNo) <= sea_boundary) THEN
            !Get indices of two adjacent triangles
            il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
            ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
            il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
            ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)
            
            z_grad_u(edge_index,level,blockNo)%x = p_param%k_veloc_h(edge_index,level,blockNo)*  &
              & ( p_diag%p_vn(il_c2,level,ib_c2)%x &
              & - p_diag%p_vn(il_c1,level,ib_c1)%x)&
              & / patch_2D%edges%dual_edge_length(edge_index,blockNo)
            
          ENDIF
        ENDDO
      END DO
    END DO
    DO idx_cartesian = 1,3
      CALL sync_patch_array(sync_e, patch_2D,z_grad_u(:,:,:)%x(idx_cartesian) )
    END DO
    
    !Step 2: Apply divergence to each component of mixing times gradient vector
    iidx => patch_2D%cells%edge_idx
    iblk => patch_2D%cells%edge_blk
    
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      
#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO level = start_level, elev
        DO cell_index = start_index, end_index
          
          IF (patch_3D%lsm_c(cell_index,level,blockNo) >= boundary) THEN
            z_div_grad_u(cell_index,level,blockNo)%x = 0.0_wp
          ELSE
            z_div_grad_u(cell_index,level,blockNo)%x =  &
              & z_grad_u(iidx(cell_index,blockNo,1),level,iblk(cell_index,blockNo,1))%x&
              & * p_op_coeff%div_coeff(cell_index,level,blockNo,1)+&
              & z_grad_u(iidx(cell_index,blockNo,2),level,iblk(cell_index,blockNo,2))%x&
              & * p_op_coeff%div_coeff(cell_index,level,blockNo,2)+&
              & z_grad_u(iidx(cell_index,blockNo,3),level,iblk(cell_index,blockNo,3))%x&
              & * p_op_coeff%div_coeff(cell_index,level,blockNo,3)
            
          ENDIF
        END DO
      END DO
    END DO
    DO idx_cartesian = 1,3
      CALL sync_patch_array(sync_c, patch_2D,z_div_grad_u(:,:,:)%x(idx_cartesian) )
    END DO
    
    !Step 3: Map divergence back to edges
    CALL map_cell2edges_3d( patch_3D, z_div_grad_u, laplacian_vn_out, p_op_coeff)
    
    ! write(*,*)'lapla',maxval(p_param%K_veloc_h(:,1,:)*patch_2D%edges%primal_edge_length),&
    ! &minval(p_param%K_veloc_h(:,1,:)*patch_2D%edges%primal_edge_length)
    ! !  CALL sync_patch_array(SYNC_E, patch_2D, laplacian_vn_out)
  END SUBROUTINE veloc_diff_harmonic_div_grad
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates horizontal diffusion of edge velocity via bilaplacian diffusion
  !!    implemented as P^T (div K grad(divgrad P v)). Due to the position of the mixing matrix
  !!    which is following MPI-OM, the operator can not be written as simple iteration of
  !!    the laplacian in divgrad form.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  !<Optimize:inUse>
  SUBROUTINE veloc_diff_biharmonic_div_grad( patch_3D, p_param, p_diag,&
    & p_op_coeff, laplacian_vn_out)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    !REAL(wp), INTENT(in)              :: vn_in(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_ho_params), INTENT(in)     :: p_param !mixing parameters
    TYPE(t_hydro_ocean_diag)          :: p_diag
    TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)           :: laplacian_vn_out(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    
    !Local variables
    INTEGER :: start_level, elev
    INTEGER :: level, blockNo, edge_index,cell_index
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: start_index, end_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: idx_cartesian
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_cartesian_coordinates) :: z_grad_u          (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: z_div_grad_u      (nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, cells_oneEdgeInDomain
    TYPE(t_subset_range), POINTER :: all_edges, edges_in_domain, edges_gradIsCalculable
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_oce_diffusion:velocity_diffusion_horz')
    !-------------------------------------------------------------------------------
    patch_2D         => patch_3D%p_patch_2d(1)
    all_cells       => patch_2D%cells%ALL
    cells_in_domain       => patch_2D%cells%in_domain
    cells_oneEdgeInDomain => patch_2D%cells%one_edge_in_domain
    all_edges       => patch_2D%edges%ALL
    edges_in_domain => patch_2D%edges%in_domain
    edges_gradIsCalculable => patch_2D%edges%gradIsCalculable
    !-------------------------------------------------------------------------------
    start_level = 1
    elev = n_zlev
    
!     laplacian_vn_out  (1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e) = 0.0_wp
    
    ! loop over cells in local domain + halo
!     DO blockNo = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, blockNo, start_index, end_index)
!       DO level = start_level, elev
!         DO cell_index = start_index, end_index
!           z_div_grad_u(cell_index,level,blockNo)%x =  0.0_wp
!         END DO
!       END DO
!     END DO
    
!ICON_OMP_PARALLEL PRIVATE(iidx, iblk )
    !-------------------------------------------------------------------------------------------------------
    !Step 1: Calculate gradient of cell velocity vector.
    !Result is a gradient vector, located at edges
!ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2 ) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_gradIsCalculable%start_block, edges_gradIsCalculable%end_block
   ! DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(edges_gradIsCalculable, blockNo, start_edge_index, end_edge_index)

      DO edge_index = start_edge_index, end_edge_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)

          !Get indices of two adjacent triangles
          il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
          ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
          il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
          ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)

          z_grad_u(edge_index,level,blockNo)%x = & !p_param%k_veloc_h(edge_index,level,blockNo)*&
            & (p_diag%p_vn(il_c2,level,ib_c2)%x - p_diag%p_vn(il_c1,level,ib_c1)%x)&
            & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
        ENDDO
        ! zero the land levels
        DO level = patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)+1, elev
          z_grad_u(edge_index,level,blockNo)%x = 0.0_wp
        ENDDO
      END DO
    END DO
!ICON_OMP_END_DO
    
!     DO idx_cartesian = 1,3
!       CALL sync_patch_array(sync_e, patch_2D,z_grad_u(:,:,:)%x(idx_cartesian) )
!     END DO
    
    !Step 2: Apply divergence to each component of gradient vector
    iidx => patch_2D%cells%edge_idx
    iblk => patch_2D%cells%edge_blk
    
!ICON_OMP_DO PRIVATE(start_index, end_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_oneEdgeInDomain%start_block, cells_oneEdgeInDomain%end_block
    ! DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(cells_oneEdgeInDomain, blockNo, start_index, end_index)
      DO cell_index = start_index, end_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(cell_index, blockNo)
        
          z_div_grad_u(cell_index,level,blockNo)%x =  &
            & z_grad_u(iidx(cell_index,blockNo,1),level,iblk(cell_index,blockNo,1))%x &
            & * p_op_coeff%div_coeff(cell_index,level,blockNo,1) + &
            & z_grad_u(iidx(cell_index,blockNo,2),level,iblk(cell_index,blockNo,2))%x &
            & * p_op_coeff%div_coeff(cell_index,level,blockNo,2) + &
            & z_grad_u(iidx(cell_index,blockNo,3),level,iblk(cell_index,blockNo,3))%x &
            & * p_op_coeff%div_coeff(cell_index,level,blockNo,3)
            
        END DO
        
        DO level = patch_3D%p_patch_1d(1)%dolic_c(cell_index, blockNo)+1, elev
          z_div_grad_u(cell_index,level,blockNo)%x = 0.0_wp
        ENDDO
        
      END DO
    END DO
!ICON_OMP_END_DO
    
!     DO idx_cartesian = 1,3
!       CALL sync_patch_array(sync_c, patch_2D,z_div_grad_u(:,:,:)%x(idx_cartesian) )
!     END DO
    
    !Step 4: Repeat the application of div and grad and take the mixing coefficients into account
    !First the grad of previous result
    !now times the mixiing/friction coefficient
!ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2 ) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO edge_index = start_edge_index, end_edge_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
          
            !Get indices of two adjacent triangles
            il_c1 = patch_2D%edges%cell_idx(edge_index,blockNo,1)
            ib_c1 = patch_2D%edges%cell_blk(edge_index,blockNo,1)
            il_c2 = patch_2D%edges%cell_idx(edge_index,blockNo,2)
            ib_c2 = patch_2D%edges%cell_blk(edge_index,blockNo,2)
            
            z_grad_u(edge_index,level,blockNo)%x = p_param%k_veloc_h(edge_index,level,blockNo) * &
              & ( z_div_grad_u(il_c2,level,ib_c2)%x      &
              &   - z_div_grad_u(il_c1,level,ib_c1)%x)   &
              & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
              
        ENDDO
      END DO
    END DO
!ICON_OMP_END_DO
!     DO idx_cartesian = 1,3
!       CALL sync_patch_array(sync_e, patch_2D,z_grad_u(:,:,:)%x(idx_cartesian) )
!     END DO
    
    !Step 5: Apply divergence to each component of gradient vector
    iidx => patch_2D%cells%edge_idx
    iblk => patch_2D%cells%edge_blk
    
!ICON_OMP_DO PRIVATE(start_index, end_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      DO cell_index = start_index, end_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(cell_index, blockNo)
          
           z_div_grad_u(cell_index,level,blockNo)%x =  - ( &        ! take the negative div in order to avoid the negation of the laplacian
              & z_grad_u(iidx(cell_index,blockNo,1),level,iblk(cell_index,blockNo,1))%x &
              & * p_op_coeff%div_coeff(cell_index,level,blockNo,1)+ &
              & z_grad_u(iidx(cell_index,blockNo,2),level,iblk(cell_index,blockNo,2))%x &
              & * p_op_coeff%div_coeff(cell_index,level,blockNo,2)+ &
              & z_grad_u(iidx(cell_index,blockNo,3),level,iblk(cell_index,blockNo,3))%x &
              & * p_op_coeff%div_coeff(cell_index,level,blockNo,3))
              
        END DO
      END DO
    END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL


    CALL sync_patch_array_mult(sync_c, patch_2D, 3,  &
     & z_div_grad_u(:,:,:)%x(1),  &
     & z_div_grad_u(:,:,:)%x(2),  &
     & z_div_grad_u(:,:,:)%x(3)   )
    
    
    !Step 6: Map divergence back to edges
    CALL map_cell2edges_3d( patch_3D, z_div_grad_u,laplacian_vn_out,p_op_coeff)! requires cells_oneEdgeInDomain 

    ! this is not needed since we take the negative div
!     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!       CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
!       DO level = start_level, elev
!         DO edge_index = start_edge_index, end_edge_index
!           
!           laplacian_vn_out(edge_index,level,blockNo)&
!           &=-laplacian_vn_out(edge_index,level,blockNo)
!  
!         ENDDO
!       END DO
!     END DO
    
    
!        DO level=1,n_zlev
!         write(*,*)'Biharmonic divgrad',level,maxval(laplacian_vn_out(:,level,:)),&
!         &minval(laplacian_vn_out(:,level,:))
!        END DO
    ! !  CALL sync_patch_array(SYNC_E, patch_2D, laplacian_vn_out)
  END SUBROUTINE veloc_diff_biharmonic_div_grad
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!  Computes  laplacian of a vector field in curl curl form.
  !!
  !! input:  lives on edges (velocity points)
  !! output: lives on edges
  !!
  !! @par Revision History
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
  !<Optimize:inUse>
  SUBROUTINE veloc_diff_harmonic_curl_curl( u_vec_e, vort, patch_3D, p_op_coeff,&
    & p_vn_dual, nabla2_vec_e,k_h )
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)      :: patch_3D
    REAL(wp), INTENT(in)                      :: u_vec_e(nproma, n_zlev, patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)                      :: vort   (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff), INTENT(in)        :: p_op_coeff
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_dual    (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout)                   :: nabla2_vec_e(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp),OPTIONAL, INTENT(in)             :: k_h(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    !
    !Local variables
    INTEGER :: start_level, elev     ! vertical start and end level
    INTEGER :: edge_index, level, blockNo
    INTEGER :: start_index, end_index
    REAL(wp) ::  z_div_c(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)!, &
    !REAL(wp) ::  z_vn_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    !REAL(wp) ::  z_rot_v(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk, ividx, ivblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    ! note that this will go through the lateral boundaries
    patch_2D         => patch_3D%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-----------------------------------------------------------------------
    
    start_level = 1
    
    icidx => patch_2D%edges%cell_idx
    icblk => patch_2D%edges%cell_blk
    ividx => patch_2D%edges%vertex_idx
    ivblk => patch_2D%edges%vertex_blk
    
    ! compute divergence of vector field
    ! z_div_c(:,:,patch_2D%alloc_cell_blocks) = 0.0_wp
    
#ifdef NAGFOR
    z_div_c = 0.0_wp
#endif

    CALL div_oce_3d( u_vec_e, patch_3D, p_op_coeff%div_coeff, z_div_c)
    CALL sync_patch_array(sync_c,patch_2D,z_div_c)
    
    ! compute rotation of vector field for the ocean
    !CALL rot_vertex_ocean_3D( patch_2D, u_vec_e, p_vn_dual, p_op_coeff, z_rot_v)!
    !CALL sync_patch_array(SYNC_V,patch_2D,z_rot_v)
    !z_rot_v=vort
    
    IF(PRESENT(k_h))THEN
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        nabla2_vec_e(:,:,blockNo) = 0.0_wp
        DO edge_index = start_index, end_index
          DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
            
            !              write(0, *) "==============================="
            !              write(0, *) "0",  edge_index,level,blockNo
            !              write(0, *) "1",  p_patch_3D%wet_e(edge_index,level,blockNo)
            !              write(0,*)  "2",  k_h(edge_index,level,blockNo)
            !              write(0,*)  "3",  patch_2D%edges%system_orientation(edge_index,blockNo)
            !              write(0,*)  "4",  vort(ividx(edge_index,blockNo,2),level,ivblk(edge_index,blockNo,2))
            !              write(0,*)  "5",  vort(ividx(edge_index,blockNo,1),level,ivblk(edge_index,blockNo,1))
            !              write(0,*)  "6",  patch_2D%edges%inv_primal_edge_length(edge_index,blockNo)
            !              write(0,*)  "7",  z_div_c(icidx(edge_index,blockNo,2),level,icblk(edge_index,blockNo,2))
            !              write(0,*)  "8",  z_div_c(icidx(edge_index,blockNo,1),level,icblk(edge_index,blockNo,1))
            !              write(0,*)  "9",  patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
            !IF(v_base%lsm_e(edge_index,level,blockNo) < land_boundary)THEN
            
            nabla2_vec_e(edge_index,level,blockNo) = patch_3D%wet_e(edge_index,level,blockNo)* &    
              & k_h(edge_index,level,blockNo) * (                                &
              & patch_2D%edges%system_orientation(edge_index,blockNo) *        &
              & ( vort(ividx(edge_index,blockNo,2),level,ivblk(edge_index,blockNo,2))         &
              & - vort(ividx(edge_index,blockNo,1),level,ivblk(edge_index,blockNo,1)) )       &
              & * patch_2D%edges%inv_primal_edge_length(edge_index,blockNo)    &
              & +                                                &
              & ( z_div_c(icidx(edge_index,blockNo,2),level,icblk(edge_index,blockNo,2))      &
              & - z_div_c(icidx(edge_index,blockNo,1),level,icblk(edge_index,blockNo,1)) )    &
              & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo))
          END DO
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ELSEIF(.NOT.PRESENT(k_h))THEN
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
        nabla2_vec_e(:,:,blockNo) = 0.0_wp
        DO edge_index = start_index, end_index
          DO level = start_level, patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
          
            nabla2_vec_e(edge_index,level,blockNo) = patch_3D%wet_e(edge_index,level,blockNo)*&
              & (   &
              & patch_2D%edges%system_orientation(edge_index,blockNo) *     &
              & ( vort(ividx(edge_index,blockNo,2),level,ivblk(edge_index,blockNo,2))     &
              & - vort(ividx(edge_index,blockNo,1),level,ivblk(edge_index,blockNo,1)) )   &
              & * patch_2D%edges%inv_primal_edge_length(edge_index,blockNo) &
              & + &
              & ( z_div_c(icidx(edge_index,blockNo,2),level,icblk(edge_index,blockNo,2))     &
              & - z_div_c(icidx(edge_index,blockNo,1),level,icblk(edge_index,blockNo,1)) )   &
              & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo))
          END DO
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO
      
    ENDIF
   
  END SUBROUTINE veloc_diff_harmonic_curl_curl
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!  Computes  Bilaplacian of a vector field.The placement of the mixing tensor
  !! follows the practice in NEMO. There is no corresponding operator in MPI-OM
  !!
  !! input:  lives on edges (velocity points)
  !! output: lives on edges
  !!
  !! @par Revision History
  !! Developed by P. Korn, MPI-M(2012)
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
  SUBROUTINE veloc_diff_biharmonic_curl_curl( patch_3D,u_vec_e,vort, k_h, p_op_coeff,&
    & p_vn_dual, nabla4_vec_e)
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch_3d ),TARGET, INTENT(in)  :: patch_3D
    REAL(wp), INTENT(in)                      :: u_vec_e(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)                      :: vort   (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(in)                      :: k_h(:,:,:)
    TYPE(t_operator_coeff), INTENT(in)        :: p_op_coeff
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_dual   (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout)                   :: nabla4_vec_e(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    
    !
    !Local variables
    INTEGER :: start_level, elev     ! vertical start and end level
    INTEGER :: edge_index, level, blockNo
    INTEGER :: start_index, end_index
    REAL(wp) ::  z_div_c   (nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) ::  z_rot_v   (nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_v)
    REAL(wp) ::  z_nabla2_e(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: h_e        (nproma,patch_3D%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates)  :: p_nabla2_dual(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_v)
    INTEGER,  DIMENSION(:,:,:), POINTER :: icidx, icblk, ividx, ivblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    ! note that this will go through the lateral boundaries
    patch_2D         => patch_3D%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-----------------------------------------------------------------------
    
    start_level = 1
    elev = n_zlev
    
    icidx => patch_2D%edges%cell_idx
    icblk => patch_2D%edges%cell_blk
    ividx => patch_2D%edges%vertex_idx
    ivblk => patch_2D%edges%vertex_blk
    
    z_nabla2_e(1:nproma,1:n_zlev,1:patch_2D%nblks_e) = 0.0_wp
    p_nabla2_dual(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v)%x(1)=0.0_wp
    p_nabla2_dual(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v)%x(2)=0.0_wp
    p_nabla2_dual(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v)%x(3)=0.0_wp
    z_div_c(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%alloc_cell_blocks) =0.0_wp   
    z_rot_v(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_v) =0.0_wp
    
    CALL veloc_diff_harmonic_curl_curl( u_vec_e, vort, patch_3D, p_op_coeff,&
      & p_vn_dual, z_nabla2_e)!, k_h)     
      
    ! compute divergence of vector field
    !     CALL div_oce_3d( u_vec_e, patch_2D, p_op_coeff%div_coeff, z_div_c)
    !     ! DO level = start_level, elev
    !     ! write(*,*)'vort1:',level,maxval(vort(:,level,:)),minval(vort(:,level,:))
    !     ! END DO
    !     ! compute rotation of vector field for the ocean
    !     !CALL rot_vertex_ocean_3D( patch_2D, u_vec_e, p_vn_dual, p_op_coeff, z_rot_v)!
    !     !z_rot_v=vort
    !     !
    !     !  loop through all patch edges (and blocks)
    !     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    !       CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
    !       DO edge_index = start_index, end_index
    !         DO level = start_level, elev
    !           !DO edge_index = start_index, end_index
    !           !IF(v_base%lsm_e(edge_index,level,blockNo) < land_boundary)THEN
    !           z_nabla2_e(edge_index,level,blockNo) =  &
    !             & v_base%wet_e(edge_index,level,blockNo)*     &
    !             & (patch_2D%edges%system_orientation(edge_index,blockNo) *  &
    !             & ( vort(ividx(edge_index,blockNo,2),level,ivblk(edge_index,blockNo,2))  &
    !             & - vort(ividx(edge_index,blockNo,1),level,ivblk(edge_index,blockNo,1)) )  &
    !             & * patch_2D%edges%inv_primal_edge_length(edge_index,blockNo))  &
    !             & +v_base%wet_e(edge_index,level,blockNo)*&
    !             & (( z_div_c(icidx(edge_index,blockNo,2),level,icblk(edge_index,blockNo,2))    &
    !             & - z_div_c(icidx(edge_index,blockNo,1),level,icblk(edge_index,blockNo,1)) )  &
    !             & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo))
    !            z_nabla2_e(edge_index,level,blockNo)=sqrt(k_h(edge_index,level,blockNo))*z_nabla2_e(edge_index,level,blockNo)
    !         END DO
    !       END DO
    !     END DO
    
    ! compute divergence of vector field
    CALL div_oce_3d( z_nabla2_e, patch_3D, p_op_coeff%div_coeff, z_div_c)
    
    ! compute rotation of vector field for the ocean
    CALL map_edges2vert_3d( patch_2D, &
      & z_nabla2_e,&
      & p_op_coeff%edge2vert_coeff_cc,&
      & p_nabla2_dual)
      
    CALL rot_vertex_ocean_3d( patch_3D, z_nabla2_e, p_nabla2_dual, p_op_coeff, z_rot_v)
    
    !combine divergence and vorticity
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_index, end_index)
      DO edge_index = start_index, end_index
        DO level = start_level, elev
          
          nabla4_vec_e(edge_index,level,blockNo) =  &
            & patch_3D%wet_e(edge_index,level,blockNo)&
            &* k_h(edge_index,level,blockNo)  &
            & *(patch_2D%edges%system_orientation(edge_index,blockNo) *  &
            & ( z_rot_v(ividx(edge_index,blockNo,2),level,ivblk(edge_index,blockNo,2))  &
            & - z_rot_v(ividx(edge_index,blockNo,1),level,ivblk(edge_index,blockNo,1)) )  &
            & * patch_2D%edges%inv_primal_edge_length(edge_index,blockNo))   &
            & +patch_3D%wet_e(edge_index,level,blockNo)                      &
            &* k_h(edge_index,level,blockNo) * &            
            & (( z_div_c(icidx(edge_index,blockNo,2),level,icblk(edge_index,blockNo,2))    &
            & - z_div_c(icidx(edge_index,blockNo,1),level,icblk(edge_index,blockNo,1)) )  &
            & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo))

           nabla4_vec_e(edge_index,level,blockNo)=- nabla4_vec_e(edge_index,level,blockNo)

        END DO
      END DO
    END DO

      DO level=1,n_zlev
       write(*,*)'Biharmonic curlcurls',level,maxval(nabla4_vec_e(:,level,:)),&
       &minval(nabla4_vec_e(:,level,:))
      END DO


  END SUBROUTINE veloc_diff_biharmonic_curl_curl
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !Subroutine computes the horizontal diffusive flux of an arbitrary tracer.
   SUBROUTINE tracer_diffusion_horz(patch_3D, trac_in, p_os, diff_flx, k_t, subset_range)
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)              :: trac_in(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_hydro_ocean_state), TARGET :: p_os
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
    INTEGER :: start_level, elev
    INTEGER :: cell_index, level, blockNo
    INTEGER :: start_index, end_index
    INTEGER :: z_dolic
    ! vertical diffusive tracer flux
    REAL(wp)                      :: z_diff_flx(nproma, n_zlev+1,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_oce_diffusion:tracer_diffusion_vert')
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2d(1)
    all_cells => patch_2D%cells%ALL
    !-----------------------------------------------------------------------
    start_level              = 1
    elev              = n_zlev
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
      ! b is diagonal, a is the lower diagonal, c is the upper
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

      DO level = 1, bottom_level
        fact(level) = diagonal_product / b(level)
        a(level)  = a(level)  * fact(level)
        b(level)  = diagonal_product
        c(level)  = dt_inv * fact(level) - a(level) - b(level)
        column_tracer(level) = field_column(cell_index,level,blockNo) * dt_inv * fact(level)
      ENDDO
      c(bottom_level) = 0.0_wp

      !------------------------------------
      ! solver from lapack
      !
      ! eliminate upper diagonal
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
  
  
  
  !-------------------------------------------------------------------------
  !
  !!Subroutine implements implicit vertical diffusion for horizontal velocity fields
  !!by inverting a scalar field..
  !>
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !! Optimized ofr round-off erros, Leonidas Linardakis, MPI-M (2014)
  !------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE velocity_diffusion_vertical_implicit( patch_3d,           &
    & velocity, a_v,                                                   &
    & operators_coefficients ) !,  &
    ! & diff_column)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)              :: velocity(:,:,:)   ! on edges
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)
    TYPE(t_operator_coeff),TARGET :: operators_coefficients
    !
    INTEGER :: start_index, end_index, edge_block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    
    !-----------------------------------------------------------------------
    edges_in_domain       =>  patch_3d%p_patch_2d(1)%edges%in_domain
    !-----------------------------------------------------------------------
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO edge_block = edges_in_domain%start_block, edges_in_domain%end_block    
      CALL get_index_range(edges_in_domain, edge_block, start_index, end_index)
      
      CALL velocity_diffusion_vertical_implicit_onBlock(  &
        & patch_3d,                                       &
        & velocity(:,:,edge_block),                       &
        & a_v(:,:,edge_block),                            &
        & operators_coefficients,                         &
        & start_index, end_index, edge_block)
        
    END DO
!ICON_OMP_END_PARALLEL_DO 
    
  END SUBROUTINE velocity_diffusion_vertical_implicit
  !------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !
  !!Subroutine implements implicit vertical diffusion for horizontal velocity fields
  !!by inverting a scalar field..
  !>
  !! sbr identical to previous one, except for homogeneous boundary conditions
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !! Optimized ofr round-off erros, Leonidas Linardakis, MPI-M (2014)
  !------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE velocity_diffusion_vertical_implicit_onBlock( &
    & patch_3d,                            &
    & velocity,                            &
    & a_v,                                 &
    & operators_coefficients,              &
    & start_index, end_index, edge_block)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)              :: velocity(:,:)   ! on edges, (nproma, levels)
    REAL(wp), INTENT(inout)              :: a_v(:,:)      ! on edges, (nproma, levels)
    TYPE(t_operator_coeff),TARGET :: operators_coefficients
    INTEGER , INTENT(in):: start_index, end_index, edge_block
    !
    REAL(wp) :: dt_inv
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev), diagonal_product
    REAL(wp) :: column_velocity(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)

    INTEGER :: bottom_level
    INTEGER :: edge_index, level
    TYPE(t_patch), POINTER :: patch_2d

    !-----------------------------------------------------------------------
    dt_inv = 1.0_wp/dtime

    DO edge_index = start_index, end_index
      bottom_level = patch_3d%p_patch_1d(1)%dolic_e(edge_index,edge_block)

      IF (bottom_level < 2 ) CYCLE ! nothing to diffuse

      ! Note : the inv_prism_thick_e, inv_prism_center_dist_e should be updated in calculate_thickness
      DO level=1, bottom_level
        inv_prism_thickness(level)        = patch_3d%p_patch_1d(1)%inv_prism_thick_e(edge_index,level,edge_block)
        inv_prisms_center_distance(level) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_e(edge_index,level,edge_block)
      ENDDO

      !------------------------------------
      ! Fill triangular matrix
      ! b is diagonal, a is the lower diagonal, c is the upper
      !   top level
      a(1) = 0.0_wp
      c(1) = -a_v(edge_index,2) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
      b(1) = dt_inv - c(1)
      DO level = 2, bottom_level-1
        a(level) = - a_v(edge_index,level)   * inv_prism_thickness(level) * inv_prisms_center_distance(level)
        c(level) = - a_v(edge_index,level+1) * inv_prism_thickness(level) * inv_prisms_center_distance(level+1)
        b(level) = dt_inv - a(level) - c(level)
      END DO
      ! bottom
      a(bottom_level) = -a_v(edge_index,bottom_level) * inv_prism_thickness(bottom_level) * &
        & inv_prisms_center_distance(bottom_level)
      b(bottom_level) = dt_inv - a(bottom_level)

      ! precondition: set diagonal equal to diagonal_product
      diagonal_product = PRODUCT(b(1:bottom_level))

      DO level = 1, bottom_level
        fact(level) = diagonal_product / b(level)
        a(level)  = a(level)  * fact(level)
        b(level)  = diagonal_product
        c(level)  = dt_inv * fact(level) - a(level) - b(level)
        column_velocity(level) = velocity(edge_index,level) * dt_inv * fact(level)
      ENDDO
      c(bottom_level) = 0.0_wp

      !------------------------------------
      ! solver from lapack
      !
      ! eliminate upper diagonal
      DO level = bottom_level-1, 1, -1
        fact(level+1)  = c( level ) / b( level+1 )
        b( level ) = b( level ) - fact(level+1) * a( level +1 )
        column_velocity( level ) = column_velocity( level ) - fact(level+1) * column_velocity( level+1 )
      ENDDO

      !     Back solve with the matrix U from the factorization.
      column_velocity( 1 ) = column_velocity( 1 ) / b( 1 )
      DO level = 2, bottom_level
        column_velocity( level ) = ( column_velocity( level ) - a( level ) * column_velocity( level-1 ) ) / b( level )
      ENDDO

      DO level = 1, bottom_level
        velocity(edge_index,level) = column_velocity(level)
      ENDDO

    END DO ! edge_index = start_index, end_index

  END SUBROUTINE velocity_diffusion_vertical_implicit_onBlock
  !------------------------------------------------------------------------

!   !-------------------------------------------------------------------------
!   !!Subroutine implements implicit vertical diffusion for scalar fields.
!   !>
!   !! sbr identical to sbr above but now with homogeneous boundary conditions
!   !!
!   !! @par Revision History
!   !! Developed  by  Peter Korn, MPI-M (2011).
!   !! Preconditioning by Leonidas Linardakis, MPI-M (2014)
!   !!
!   !! The result ocean_tracer%concetration is calculated on domain_cells
!   !-------------------------------------------------------------------------
!   SUBROUTINE tracer_diffusion_vertical_implicit( patch_3D,               &
!                                            & ocean_tracer, A_v,   &
!                                            & operators_coefficients ) !,  &
!                                           ! & diff_column)
! 
!     TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
!     TYPE(t_ocean_tracer), TARGET         :: ocean_tracer
!     REAL(wp), INTENT(inout)              :: A_v(:,:,:)
!     TYPE(t_operator_coeff),TARGET        :: operators_coefficients
!     !
!     !
!     REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
!     REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)! , nb(1:n_zlev)
!     REAL(wp) :: fact(1:n_zlev)
!     REAL(wp) :: column_tracer(1:n_zlev)
!     REAL(wp) :: dt_inv, diagonal_product
!     REAL(wp), POINTER   :: field_column(:,:,:)
!     INTEGER  :: bottom_level
!     INTEGER :: jc, jk, jb
!     INTEGER :: start_index, end_index
!     TYPE(t_subset_range), POINTER :: cells_in_domain
!     TYPE(t_patch), POINTER         :: patch_2D
!     !-----------------------------------------------------------------------
!     patch_2D        => patch_3D%p_patch_2D(1)
!     cells_in_domain => patch_2D%cells%in_domain
!     field_column    => ocean_tracer%concentration
!     !-----------------------------------------------------------------------
!     dt_inv = 1.0_wp/dtime
! 
!     DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!       CALL get_index_range(cells_in_domain, jb, start_index, end_index)
!       DO jc = start_index, end_index
!         bottom_level = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
! 
!         IF (bottom_level < 1 ) CYCLE
! 
!         DO jk=1,bottom_level
!           inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
!           inv_prisms_center_distance(jk) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,jb)
!         ENDDO
! 
!         !------------------------------------
!         ! Fill triangular matrix
!         ! b is diagonal, a is the lower diagonal, c is the upper
!         !   top level
!         a(1) = 0.0_wp
!         c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
!         b(1) = dt_inv - c(1)
!         DO jk = 2, bottom_level-1
!           a(jk) = - A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
!           c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
!           b(jk) = dt_inv - a(jk) - c(jk)
!         END DO
!         ! bottom
!         a(bottom_level) = -A_v(jc,bottom_level,jb) * inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)
!         b(bottom_level) = dt_inv - a(bottom_level)
! 
!         ! precondition: set diagonal equal to diagonal_product
!         diagonal_product = PRODUCT(b(1:bottom_level))
! 
!         DO jk = 1, bottom_level
!            fact(jk) = diagonal_product / b(jk)
!            a(jk)  = a(jk)  * fact(jk)
!            b(jk)  = diagonal_product
!            c(jk)  = dt_inv * fact(jk) - a(jk) - b(jk)
!            column_tracer(jk) = field_column(jc,jk,jb) * dt_inv * fact(jk)
!         ENDDO
!         c(bottom_level) = 0.0_wp
! 
!         !------------------------------------
!         ! solver from lapack
!         !
!         ! eliminate upper diagonal
!         DO jk=bottom_level-1, 1, -1
!           fact(jk+1)  = c( jk ) / b( jk+1 )
!           b( jk ) = b( jk ) - fact(jk+1) * a( jk +1 )
!           column_tracer( jk ) = column_tracer( jk ) - fact(jk+1) * column_tracer( jk+1 )
!         ENDDO
! 
!         !     Back solve with the matrix U from the factorization.
!         column_tracer( 1 ) = column_tracer( 1 ) / b( 1 )
!         DO jk =  2, bottom_level
!           column_tracer( jk ) = ( column_tracer( jk ) - a( jk ) * column_tracer( jk-1 ) ) / b( jk )
!         ENDDO
! 
!         DO jk = 1, bottom_level
!           ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk)
!         ENDDO
! 
!       END DO ! jc = start_index, end_index
!     END DO ! jb = cells_in_domain%start_block, cells_in_domain%end_block
! 
!   END SUBROUTINE tracer_diffusion_vertical_implicit
!   !------------------------------------------------------------------------

  
END MODULE mo_oce_diffusion
