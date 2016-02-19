!>
!!               The module <i>mo_scalar_product</i>
!! implements discrete scalar products which depend on the grid geometry only
!! are used to formulate the primitive equations in weak form.
!! The coefficients used for projections from triangular edges to cell centers
!! and vice versa are contained in the ocean part of the model domain and
!! are calculated in <i>mo_ocean_topo</i>.
!!
!! @par Revision History
!! Initial version  by Peter Korn and Stephan Lorenz,  MPI-M, Hamburg, October 2010
!! Modification by Stephan Lorenz, MPI-M, (2010-11-02)
!! - initial primal_flip_flop divided into basic parts of primal_map_e2c and c2e
!! Modification by Stephan Lorenz, MPI-M, (2010-11-16)
!! - implementation as primal_map_e2c, primal_map_e2c_no_edge_height, primal_map_c2e
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
#include "icon_definitions.inc"
!----------------------------
MODULE mo_scalar_product
  !-------------------------------------------------------------------------
  !
  USE mo_kind,               ONLY: wp, sp
  USE mo_exception,          ONLY: finish
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: sea_boundary, sea, min_dolic
  USE mo_model_domain,       ONLY: t_patch, t_patch_3d
  USE mo_ocean_types,        ONLY: t_hydro_ocean_diag, t_solvercoeff_singleprecision
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce, fast_performance_level,l_ANTICIPATED_VORTICITY,&
    & KineticEnergy_type, KineticEnergy_onPrimalGrid, KineticEnergy_onDualGrid
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff, no_primal_edges, no_dual_edges, &
    & Get3DVectorTo2DLocal_array3D
  USE mo_ocean_math_operators,ONLY: grad_fd_norm_oce_3D, rot_vertex_ocean_3d, map_edges2vert_3d
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_e, sync_v,sync_patch_array,sync_c, sync_patch_array_mult!,  & sync_idx, global_max
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_timer,               ONLY: timers_level, timer_start, timer_stop, timer_extra25, timer_extra26, timer_extra27, &
    & timer_extra28, timer_extra29


  
  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: calc_scalar_product_veloc_3d
  PUBLIC :: nonlinear_coriolis_3d
  PUBLIC :: map_edges2vert_3d
  PUBLIC :: map_edges2cell_3d
  PUBLIC :: map_cell2edges_3d
  !PRIVATE :: map_edges2edges_viacell_3d
  PUBLIC :: map_edges2edges_viavert_3d
  PUBLIC :: map_edges2edges_viacell_3d_const_z
  PUBLIC :: map_edges2edges_viacell_2d_constZ_onTriangles_sp
  PUBLIC :: map_vec_prismtop2center_on_block
  PUBLIC :: map_scalar_prismtop2center
  PUBLIC :: map_scalar_center2prismtop
  PUBLIC :: map_edges2edges_viacell_2D_per_level
  
!   INTERFACE map_edges2edges_viacell_3d
!     MODULE PROCEDURE map_edges2edges_viacell_3d_1lev
!     MODULE PROCEDURE map_edges2edges_viacell_3d_mlev
!   END INTERFACE
  
  INTERFACE map_edges2edges_viacell_3d_const_z
    MODULE PROCEDURE map_edges2edges_viacell_2D_constZ
!     MODULE PROCEDURE map_edges2edges_viacell_2D_constZs
    MODULE PROCEDURE map_edges2edges_viacell_3d_mlev_const_z
    MODULE PROCEDURE map_edges2edges_viacell_3d_mlev_constZs
  END INTERFACE
  
  
  INTERFACE map_cell2edges_3d
    MODULE PROCEDURE map_cell2edges_3d_1level
    MODULE PROCEDURE map_cell2edges_3d_mlevels
  END INTERFACE
  
  INTERFACE map_edges2cell_3d    
    !MODULE PROCEDURE map_edges2cell_with_height_3d
    MODULE PROCEDURE map_edges2cell_no_height_3d    
  END INTERFACE
  CHARACTER(LEN=16)          :: module_name    = 'mo_scalarprod'  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug  
CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !<Optimize:inUse>
  SUBROUTINE calc_scalar_product_veloc_3d( patch_3d, vn_e,&
    & p_diag, operators_coefficients)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)      :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_hydro_ocean_diag)  :: p_diag
    TYPE(t_operator_coeff)    :: operators_coefficients
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: cell_index, blockNo, level
    INTEGER :: edge_idx, edge_blk, n
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk
    REAL(wp) :: w1, w2, w3, w4, w
    TYPE(t_subset_range), POINTER :: all_cells ! , cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    all_cells       => patch_2d%cells%all
    
    edge_of_cell_idx  => patch_2d%cells%edge_idx
    edge_of_cell_blk  => patch_2d%cells%edge_blk    
!     cells_in_domain => patch_2d%cells%in_domain
    !-----------------------------------------------------------------------
    startLevel = 1
    endLevel = n_zlev

    start_detail_timer(timer_extra25,5)
    CALL map_edges2vert_3d(patch_3d%p_patch_2d(1), vn_e, operators_coefficients%edge2vert_coeff_cc, &
      & p_diag%p_vn_dual)
    stop_detail_timer(timer_extra25,5)

    start_detail_timer(timer_extra26,5)
    !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
    CALL map_edges2cell_3d(patch_3d, vn_e, operators_coefficients, p_diag%p_vn) !, subset_range=cells_in_domain)
    stop_detail_timer(timer_extra26,5)

    start_detail_timer(timer_extra27,5)
    CALL map_edges2edges_viacell_3d_const_z( patch_3D,       &
                                    & vn_e,                  &
                                    & operators_coefficients,&
                                    & p_diag%ptp_vn)    
    stop_detail_timer(timer_extra27,5)
    ! CALL sync_patch_array(sync_e, patch_2d, p_diag%ptp_vn)
   
    start_detail_timer(timer_extra28,5)
    !--------------------------------------------------------------
    !calculate kinetic energy
    ! First option is our default configuration
    SELECT CASE(KineticEnergy_type)
    CASE (KineticEnergy_onDualGrid)
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)

!         p_diag%kin(:,:,blockNo)=0.0_wp        
        DO cell_index =  start_cell_index, end_cell_index
          DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)      

            p_diag%kin(cell_index,level,blockNo) = 0.5_wp * &
              & DOT_PRODUCT(p_diag%p_vn(cell_index,level,blockNo)%x, p_diag%p_vn(cell_index,level,blockNo)%x)
              
          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    CASE (KineticEnergy_onPrimalGrid)
      IF (no_primal_edges /= 3) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level, n, &
!ICON_OMP edge_idx, edge_blk, w1, w2, w3, w4, w) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)

          p_diag%kin(:,:,blockNo)=0.0_wp
          DO cell_index =  start_cell_index, end_cell_index
            w = 0.0_wp
            DO n=1,patch_2d%cells%max_connectivity
              edge_idx = edge_of_cell_idx(cell_index,blockNo,n)
              edge_blk = edge_of_cell_blk(cell_index,blockNo,n)
              IF (edge_idx > 0) THEN
                w1 = patch_2d%edges%area_edge(edge_idx, edge_blk)
                w  = w + w1
                DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
                  p_diag%kin(cell_index,level,blockNo) = p_diag%kin(cell_index,level,blockNo) + &
                    & vn_e(edge_idx, level, edge_blk)**2 * w1
                ENDDO
              ENDIF
            ENDDO

            DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
              p_diag%kin(cell_index,level,blockNo) = 0.5_wp * &
                (p_diag%kin(cell_index,level,blockNo) / w)! -   &
               ! ((p_diag%w(cell_index,level,blockNo)+p_diag%w(cell_index,level+1,blockNo))*0.5_wp)**2)
            END DO

          END DO
        END DO ! block
!ICON_OMP_END_PARALLEL_DO
!         CALL finish("calc_scalar_product_veloc_3d", "NONLINEAR_CORIOLIS_PRIMAL_GRID works only for triangles")
      ELSE
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level, w1, w2, w3, w4, w) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)

  !         p_diag%kin(:,:,blockNo)=0.0_wp

          DO cell_index =  start_cell_index, end_cell_index
            w1 =  patch_2d%edges%area_edge      &
              &      (edge_of_cell_idx(cell_index,blockNo,1),edge_of_cell_blk(cell_index,blockNo,1))
            w2 =  patch_2d%edges%area_edge      &
              &      (edge_of_cell_idx(cell_index,blockNo,2),edge_of_cell_blk(cell_index,blockNo,2))
            w3 =  patch_2d%edges%area_edge      &
              &      (edge_of_cell_idx(cell_index,blockNo,3),edge_of_cell_blk(cell_index,blockNo,3))
  !           w4 = 2.0_wp * patch_2d%cells%area(cell_index,blockNo)
            w = w1 + w2 + w3
            DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)

              p_diag%kin(cell_index,level,blockNo)= (&!p_diag%kin(cell_index,level,blockNo)&
                & vn_e(edge_of_cell_idx(cell_index,blockNo,1),level,edge_of_cell_blk(cell_index,blockNo,1))**2 &
                !*p_diag%ptp_vn                         (edge_of_cell_idx(cell_index,blockNo,1),level,edge_of_cell_blk(cell_index,blockNo,1))&
                & *  w1 &
                &+&
                & vn_e(edge_of_cell_idx(cell_index,blockNo,2),level,edge_of_cell_blk(cell_index,blockNo,2))**2 &
                !&*p_diag%ptp_vn                         (edge_of_cell_idx(cell_index,blockNo,2),level,edge_of_cell_blk(cell_index,blockNo,2))&
                & * w2 &
                &+&
                & vn_e(edge_of_cell_idx(cell_index,blockNo,3),level,edge_of_cell_blk(cell_index,blockNo,3))**2 &
                ! &*p_diag%ptp_vn                         (edge_of_cell_idx(cell_index,blockNo,3),level,edge_of_cell_blk(cell_index,blockNo,3))&
                & * w3 &
                & ) / w

              p_diag%kin(cell_index,level,blockNo) = 0.5_wp * &
                & (p_diag%kin(cell_index,level,blockNo) - &
                & ((p_diag%w(cell_index,level,blockNo)+p_diag%w(cell_index,level+1,blockNo))*0.5_wp)**2)

            END DO
          END DO
        END DO ! block
!ICON_OMP_END_PARALLEL_DO
      ENDIF

    END SELECT
          
    stop_detail_timer(timer_extra28,5)
!       CALL sync_patch_array(sync_c, patch_2d,p_diag%kin)
      
      !convert cartesian velocity vector p_diag%p_vn(cell_index,level,blockNo)%x to geographical coordinate system
      !for output, sea-ice and coupling
    start_detail_timer(timer_extra29,5)
    CALL Get3DVectorTo2DLocal_array3D(vector=p_diag%p_vn, &
      & position_local=patch_2d%cells%center, &
      & levels=patch_3d%p_patch_1d(1)%dolic_c, &
      & subset=all_cells,                      &
      & geometry_info=patch_2d%geometry_info, &
      & x=p_diag%u, y=p_diag%v)
    stop_detail_timer(timer_extra29,5)
 
!       DO cell_index =  start_cell_index, end_cell_index
!         DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
!           CALL Get3DVectorToPlanarLocal(p_diag%p_vn(cell_index,level,blockNo), &
!             & patch_2d%cells%center(cell_index,blockNo), patch_2d%geometry_info, &
!             & p_diag%u(cell_index,level,blockNo), p_diag%v(cell_index,level,blockNo) )
!         END DO
!       END DO
      
    !--------------------------------------------------------------    
    
    CALL dbg_print('veloc_3d: vn_e',                   vn_e,      module_name,1, &
          patch_2D%edges%owned )
    CALL dbg_print('veloc_3d: kin energy',             p_diag%kin,module_name,1, &
          patch_2D%cells%owned )
   CALL dbg_print('veloc_3d: East-West  :u',p_diag%u,  module_name,1, &
         patch_2D%cells%owned )
    CALL dbg_print('veloc_3d: North-South :v',p_diag%v, module_name,1, &
         patch_2D%cells%owned )

          
  END SUBROUTINE calc_scalar_product_veloc_3d
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Optimized version of the nonlinear_coriolis_3d
  !! Note: intel -o3 will produce very different results when using this version,
  !!  the model exhibits great sensitivity using the AtlanticBoxACC setup
  !! Note: vn should be zero on land edges, otherwise it will give wrong result.
  !!
  !!Option for Anticipated Vorticity added. Code for this option is mostly identical to
  !! code without this option, but optimization suggests a bit of code doubling. 
  !<Optimize:inUse>
  SUBROUTINE nonlinear_coriolis_3d_fast(patch_3d, vn, p_vn_dual, vort_v, &
    & operators_coefficients, vort_flux)
    TYPE(t_patch_3d ),TARGET,INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(inout)  :: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout)                    :: vort_v   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in)          :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    
    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: startLevel! , endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, vertex_edge
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk
    INTEGER ::  edgeOfVertex_index, edgeOfVertex_block
    REAL(wp):: this_vort_flux(n_zlev, 2) ! for each of the two vertices
    REAL(wp):: thick_edge(n_zlev,2), thick_vert(n_zlev,2)
    REAL(wp):: numOfEdges(n_zlev,2)
    
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp) :: vort_flux_old(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    startLevel    = 1
    ! endLevel    = n_zlev
    
    CALL rot_vertex_ocean_3d( patch_3d, vn, p_vn_dual, operators_coefficients, vort_v)
    ! this is not needed, since vort_v is on vertices in domain
    ! CALL sync_patch_array(SYNC_V, patch_2D, vort_v)

    IF(.NOT.l_ANTICIPATED_VORTICITY)THEN
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,level,je,start_edge_index,end_edge_index, this_vort_flux, &
!ICON_OMP  vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk,  &
!ICON_OMP vertex_edge, ictr, thick_vert, numOfEdges, edgeOfVertex_index, edgeOfVertex_block, &
!ICON_OMP thick_edge) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      ! vort_flux(:,:,blockNo) = 0.0_wp
      
      DO je =  start_edge_index, end_edge_index
        
        this_vort_flux(:,:) = 0.0_wp
        
        vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
        vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
        vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
        vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)

        ! vertex 1
        ictr = 0
        
        thick_vert(1:n_zlev,1)=0.0_wp
        thick_vert(1:n_zlev,2)=0.0_wp
        numOfEdges(1:n_zlev,1)=0.0_wp
        numOfEdges(1:n_zlev,2)=0.0_wp
        
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex1_idx,vertex1_blk)

          ictr =ictr+1
          edgeOfVertex_index = patch_2d%verts%edge_idx(vertex1_idx,vertex1_blk,vertex_edge)
          edgeOfVertex_block = patch_2d%verts%edge_blk(vertex1_idx,vertex1_blk,vertex_edge)

          DO level = startLevel, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,blockNo), &
            &                        patch_3d%p_patch_1d(1)%dolic_e(edgeOfVertex_index, edgeOfVertex_block))

            numOfEdges(level,1) = numOfEdges(level,1) + 1.0_wp
            thick_edge(level,1) = patch_3D%p_patch_1d(1)%prism_thick_e(                  &
              & edgeOfVertex_index, level, edgeOfVertex_block)

            ! This should be calculated in the calculate_thickness routine
            thick_vert(level,1) = thick_vert(level,1)+thick_edge(level,1)
            
            this_vort_flux(level, 1) =  this_vort_flux(level, 1) +                        &
              & vn( edgeOfVertex_index, level, edgeOfVertex_block)                        &
              & * operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)   &
              & * thick_edge(level,1)

          ENDDO

        END DO ! edges of this vertex
        
        ! vertex 2
        ictr = no_dual_edges
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex2_idx,vertex2_blk)!no_dual_cell_edges

          ictr =ictr+1
          edgeOfVertex_index = patch_2d%verts%edge_idx(vertex2_idx,vertex2_blk,vertex_edge)
          edgeOfVertex_block = patch_2d%verts%edge_blk(vertex2_idx,vertex2_blk,vertex_edge)

          DO level = startLevel, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,blockNo), &
            &                        patch_3d%p_patch_1d(1)%dolic_e(edgeOfVertex_index, edgeOfVertex_block))

            numOfEdges(level,2) = numOfEdges(level,2) + 1.0_wp
            thick_edge(level,2) = patch_3D%p_patch_1d(1)%prism_thick_e(                  &
              & edgeOfVertex_index, level, edgeOfVertex_block)
            
              
            ! This should be calculated in the calculate_thickness routine
            thick_vert(level,2) = thick_vert(level,2)+thick_edge(level,2)
          
            this_vort_flux(level, 2) =  this_vort_flux(level, 2) +                        &
              & vn( edgeOfVertex_index, level, edgeOfVertex_block)                        &
              & * operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)   &
              & * thick_edge(level,2)

          ENDDO
        END DO ! edges of this vertex
        
        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          
          vort_flux(je,level,blockNo) =  &
            & this_vort_flux(level,1) &
            ! divide by the vertex thickness
            &  * numOfEdges(level,1) / thick_vert(level,1)  &
            &  * (vort_v(vertex1_idx, level, vertex1_blk) + patch_2d%verts%f_v(vertex1_idx, vertex1_blk))  &
            & + this_vort_flux(level,2)  &
            ! divide by the vertex thickness
            &  * numOfEdges(level,2) / thick_vert(level,2)  &
            &  * (vort_v(vertex2_idx, level, vertex2_blk) + patch_2d%verts%f_v(vertex2_idx, vertex2_blk))

        ENDDO
        
      END DO
      
    END DO ! blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!ICON_OMP_END_PARALLEL_DO

  ELSEIF(l_ANTICIPATED_VORTICITY)THEN

    !This is only for debugging needed
    vort_flux_old(:,:,:)=0.0_wp

!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,level,je,start_edge_index,end_edge_index, this_vort_flux, &
!ICON_OMP  vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk,  &
!ICON_OMP vertex_edge, ictr, thick_vert, numOfEdges, edgeOfVertex_index, edgeOfVertex_block, &
!ICON_OMP thick_edge) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      ! vort_flux(:,:,blockNo) = 0.0_wp
      
      DO je =  start_edge_index, end_edge_index
        
        this_vort_flux(:,:) = 0.0_wp
        
        vertex1_idx = patch_2d%edges%vertex_idx(je,blockNo,1)
        vertex1_blk = patch_2d%edges%vertex_blk(je,blockNo,1)
        vertex2_idx = patch_2d%edges%vertex_idx(je,blockNo,2)
        vertex2_blk = patch_2d%edges%vertex_blk(je,blockNo,2)

        ! vertex 1
        ictr = 0
        
        thick_vert(1:n_zlev,1)=0.0_wp
        thick_vert(1:n_zlev,2)=0.0_wp
        numOfEdges(1:n_zlev,1)=0.0_wp
        numOfEdges(1:n_zlev,2)=0.0_wp
        
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex1_idx,vertex1_blk)

          ictr =ictr+1
          edgeOfVertex_index = patch_2d%verts%edge_idx(vertex1_idx,vertex1_blk,vertex_edge)
          edgeOfVertex_block = patch_2d%verts%edge_blk(vertex1_idx,vertex1_blk,vertex_edge)

          DO level = startLevel, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,blockNo), &
            &                        patch_3d%p_patch_1d(1)%dolic_e(edgeOfVertex_index, edgeOfVertex_block))

            numOfEdges(level,1) = numOfEdges(level,1) + 1.0_wp
            thick_edge(level,1) = patch_3D%p_patch_1d(1)%prism_thick_e(                  &
              & edgeOfVertex_index, level, edgeOfVertex_block)

            ! This should be calculated in the calculate_thickness routine
            thick_vert(level,1) = thick_vert(level,1)+thick_edge(level,1)
            
            this_vort_flux(level, 1) =  this_vort_flux(level, 1) +                        &
              & vn( edgeOfVertex_index, level, edgeOfVertex_block)                        &
              & * operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)   &
              & * thick_edge(level,1)

          ENDDO

        END DO ! edges of this vertex
        
        ! vertex 2
        ictr = no_dual_edges
        DO vertex_edge=1, patch_2d%verts%num_edges(vertex2_idx,vertex2_blk)!no_dual_cell_edges

          ictr =ictr+1
          edgeOfVertex_index = patch_2d%verts%edge_idx(vertex2_idx,vertex2_blk,vertex_edge)
          edgeOfVertex_block = patch_2d%verts%edge_blk(vertex2_idx,vertex2_blk,vertex_edge)

          DO level = startLevel, MIN(patch_3d%p_patch_1d(1)%dolic_e(je,blockNo), &
            &                        patch_3d%p_patch_1d(1)%dolic_e(edgeOfVertex_index, edgeOfVertex_block))

            numOfEdges(level,2) = numOfEdges(level,2) + 1.0_wp
            thick_edge(level,2) = patch_3D%p_patch_1d(1)%prism_thick_e(                  &
              & edgeOfVertex_index, level, edgeOfVertex_block)
            
              
            ! This should be calculated in the calculate_thickness routine
            thick_vert(level,2) = thick_vert(level,2)+thick_edge(level,2)
          
            this_vort_flux(level, 2) =  this_vort_flux(level, 2) +                        &
              & vn( edgeOfVertex_index, level, edgeOfVertex_block)                        &
              & * operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)   &
              & * thick_edge(level,2)

          ENDDO
        END DO ! edges of this vertex
        
        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          
          vort_flux(je,level,blockNo) =  &
            & this_vort_flux(level,1) &
            ! divide by the vertex thickness
            &  * numOfEdges(level,1) / thick_vert(level,1)  &
            &  * (vort_v(vertex1_idx, level, vertex1_blk) + patch_2d%verts%f_v(vertex1_idx, vertex1_blk))  &
            & + this_vort_flux(level,2)  &
            ! divide by the vertex thickness
            &  * numOfEdges(level,2) / thick_vert(level,2)  &
            &  * (vort_v(vertex2_idx, level, vertex2_blk) + patch_2d%verts%f_v(vertex2_idx, vertex2_blk))

          vort_flux_old(je,level,blockNo)=vort_flux(je,level,blockNo)

          vort_flux(je,level,blockNo) =  vort_flux(je,level,blockNo)-&
               &   (this_vort_flux(level,1) * numOfEdges(level,1) / thick_vert(level,1)  &
               & + this_vort_flux(level,2)  * numOfEdges(level,2) / thick_vert(level,2)) &
               & *0.5_wp*(vort_v(vertex2_idx, level, vertex2_blk)-vort_v(vertex1_idx, level, vertex1_blk))&
                &/patch_2d%edges%primal_edge_length(je,blockNo)
  
!write(123,*)'details',0.5_wp*(vort_v(vertex2_idx, level, vertex2_blk)-vort_v(vertex1_idx, level, vertex1_blk))&
!&/patch_2d%edges%primal_edge_length(je,blockNo),&
!&this_vort_flux(level,1) * numOfEdges(level,1) / thick_vert(level,1)  &
!& + this_vort_flux(level,2)  * numOfEdges(level,2) / thick_vert(level,2),&
!&(this_vort_flux(level,1) * numOfEdges(level,1) / thick_vert(level,1)  &
!              & + this_vort_flux(level,2)  * numOfEdges(level,2) / thick_vert(level,2)) &
!              & *0.5_wp*(vort_v(vertex2_idx, level, vertex2_blk)-vort_v(vertex1_idx, level, vertex1_blk))&
!              &/patch_2d%edges%primal_edge_length(je,blockNo)         
        ENDDO
        
      END DO
      
    END DO ! blockNo = edges_inDomain%start_block, edges_inDomain%end_block	
!ICON_OMP_END_PARALLEL_DO

  ENDIF


  END SUBROUTINE nonlinear_coriolis_3d_fast
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Note: vn must habve been synced before this routine
  !! the resulting vort_v is synced,
  !! vort_flux id calculated on edges in_domain
!<Optimize:inUse>
  SUBROUTINE nonlinear_coriolis_3d(patch_3d, vn, p_vn_dual, vort_v, &
    & operators_coefficients, vort_flux)
    TYPE(t_patch_3d ),TARGET,INTENT(in) :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(inout)  :: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    REAL(wp), INTENT(inout)                    :: vort_v   (nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in)          :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    
    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: startLevel, endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    REAL(wp) :: vort_global, thick_edge, thick_vert
    REAL(wp) :: vort_flux_budget_perlevel(n_zlev), volume(n_zlev)
    REAL(wp) :: vort_flux_budget_perlevel_pos(n_zlev),vort_flux_budget_perlevel_neg(n_zlev)
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    
    !-----------------------------------------------------------------------
      IF (fast_performance_level > 10) THEN
        CALL nonlinear_coriolis_3d_fast(patch_3d, vn, p_vn_dual, vort_v, &
                  & operators_coefficients, vort_flux)
        RETURN
      ENDIF
     
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    startLevel    = 1
    endLevel    = n_zlev
    
    !vn=0.0_wp
    !vn(3,1:n_zlev,5)=1.0E14
    
    CALL rot_vertex_ocean_3d( patch_3d, vn, p_vn_dual, operators_coefficients, vort_v)
    ! this is not needed
    ! CALL sync_patch_array(SYNC_V, patch_2D, vort_v)
    

    ! !$OMP PARALLEL
    ! !$OMP DO PRIVATE(blockNo,level,je,start_edge_index,end_edge_index)
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      level_loop: DO level = startLevel, endLevel
        
        edge_idx_loop: DO je =  start_edge_index, end_edge_index
          
          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN
            
            vort_flux(je,level,blockNo) = 0.0_wp
            
            DO neighbor=1,2
              IF(neighbor==1) ictr = 0
              IF(neighbor==2) ictr = no_dual_edges
              
              il_v = patch_2d%edges%vertex_idx(je,blockNo,neighbor)
              ib_v = patch_2d%edges%vertex_blk(je,blockNo,neighbor)
              
              vort_global = (vort_v(il_v,level,ib_v) + patch_2d%verts%f_v(il_v,ib_v))
              DO vertex_edge=1, patch_2d%verts%num_edges(il_v,ib_v)
                
                ictr =ictr+1
                
                il_e = patch_2d%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = patch_2d%verts%edge_blk(il_v,ib_v,vertex_edge)
                
!                 thick_edge = patch_3D%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
                vort_flux(je,level,blockNo) =  vort_flux(je,level,blockNo)+vn(il_e,level,ib_e)*vort_global&
                  & *operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)!*(thick_edge/thick_vert)
              END DO
            END DO
          ELSE
            vort_flux(je,level,blockNo)= 0.0_wp
          ENDIF
        END DO edge_idx_loop
      END DO level_loop
    END DO ! blockNo = edges_inDomain%start_blockold, edges_inDomain%end_block
    ! !$OMP END DO NOWAIT
    ! !$OMP END PARALLEL
      
!     !Diagnostic of energetic neutrality
!     vort_flux_budget_perlevel(1:n_zlev)=0.0_wp
!     vort_flux_budget_perlevel_pos(1:n_zlev)=0.0_wp
!     vort_flux_budget_perlevel_neg(1:n_zlev)=0.0_wp    
!     ! !$OMP PARALLEL
!     ! !$OMP DO PRIVATE(blockNo,level,je,start_edge_index,end_edge_index)
!     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!       CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
!       DO level = startLevel, endLevel
!       volume(level)=0.0_wp
!         
!         DO je =  start_edge_index, end_edge_index
!           
!           IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN  
!           vort_flux_budget_perlevel(level)=vort_flux_budget_perlevel(level)&
!           &+vort_flux(je,level,blockNo)*vn(je,level,blockNo)&
!           &*patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)!&
!           !&*patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
!           
!           IF(vort_flux(je,level,blockNo)>0.0_wp)THEN
!            vort_flux_budget_perlevel_pos(1:n_zlev)= vort_flux_budget_perlevel_pos(1:n_zlev)&
!            &+vort_flux(je,level,blockNo)*vn(je,level,blockNo)&
!           &*patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
!           &*patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
!           ELSEIF(vort_flux(je,level,blockNo)<0.0_wp)THEN
!            vort_flux_budget_perlevel_neg(1:n_zlev)= vort_flux_budget_perlevel_neg(1:n_zlev)&
!            &+vort_flux(je,level,blockNo)*vn(je,level,blockNo)&
!           &*patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
!           &*patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
! 
!           ENDIF
!           volume(level)=volume(level)+patch_2D%edges%primal_edge_length(je,blockNo)*patch_2D%edges%dual_edge_length(je,blockNo)&
!           &*patch_3D%p_patch_1d(1)%prism_thick_e(je,level,blockNo)
!           ENDIF
!        END DO
!       END DO 
!     END DO ! blockNo = edges_inDomain%start_blockold, edges_inDomain%end_block
!     ! !$OMP END DO NOWAIT
!     ! !$OMP END PARALLEL    
!     
!     DO level = startLevel, endLevel
!     write(0,*)'global vorticity budget', level,&
!     &vort_flux_budget_perlevel(level)/volume(level),&
!     &vort_flux_budget_perlevel_pos(level)/volume(level),&
!     &vort_flux_budget_perlevel_neg(level)/volume(level)
!     write(0,*)'global vorticity budget pos+neg', level,&
!     &vort_flux_budget_perlevel_pos(level)/volume(level)+vort_flux_budget_perlevel_neg(level)/volume(level)
!     END DO
!     
!     stop
  END SUBROUTINE nonlinear_coriolis_3d
  !-------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------
  ! map_edges2cell_without_height: This sbr is used for kinetic energy calculation
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !<Optimize:inUse>
  SUBROUTINE map_edges2cell_no_height_3d( patch_3d, vn_e, operators_coefficients, p_vn_c, opt_startLevel, &
    & opt_endLevel, subset_range)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    TYPE(t_cartesian_coordinates)              :: p_vn_c(:,:,:)  ! output (nproma,n_zlev,alloc_cell_blocks)
    ! intent(inout) for nag compiler
    INTEGER, INTENT(in), OPTIONAL :: opt_startLevel,  opt_endLevel      ! optional vertical start level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: il_e, ib_e
    INTEGER :: cell_index, blockNo, level, ie
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
     IF (no_primal_edges == 3) THEN
       CALL map_edges2cell_no_height_3d_onTriangles( patch_3d, vn_e, operators_coefficients, p_vn_c, opt_startLevel, &
         & opt_endLevel, subset_range)
       RETURN
     ENDIF
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2d%cells%ALL
    ENDIF
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_startLevel) ) THEN
      startLevel = opt_startLevel
    ELSE
      startLevel = 1
    END IF
    IF ( PRESENT(opt_endLevel) ) THEN
      endLevel = opt_endLevel
    ELSE
      endLevel = n_zlev
    END IF
       
    !Calculation of Pv in cartesian coordinates
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      
      cell_idx_loop: DO cell_index =  start_cell_index, end_cell_index
 
        endLevel=patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
        level_loop: DO level = startLevel, endLevel
           !calculate velocity reconstruction at cell center
          p_vn_c(cell_index,level,blockNo)%x = 0.0_wp
          
          DO ie=1, no_primal_edges
            il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
            ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)
            IF (il_e > 0) THEN
              p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x&
                & + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,ie)%x&
                & * vn_e(il_e,level,ib_e)*patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
            ENDIF
          END DO
          
          p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x &
              & / (operators_coefficients%fixed_vol_norm(cell_index,level,blockNo) &
              &   *patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo))

!           IF (abs(p_vn_c(cell_index,level,blockNo)%x(3)) > 0.5) THEN
!             write(0,*) "---------------------------"
!             write(0,*) blockNo, cell_index, " velocity at cell:", p_vn_c(cell_index,level,blockNo)%x
!             p_vn_c(cell_index,level,blockNo)%x = 0.0_wp
!             DO ie=1, no_primal_edges
!               il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
!               ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)
!               IF (il_e > 0) THEN
!                 write(0,*) ie, ":"
!                 write(0,*) " ", patch_2d%edges%cartesian_center(il_e,ib_e)%x
!                 write(0,*) " ", patch_2d%verts%cartesian( &
!                   patch_2d%edges%vertex_idx(il_e,ib_e,1), patch_2d%edges%vertex_blk(il_e,ib_e,1)) %x
!                 write(0,*) " ", patch_2d%verts%cartesian( &
!                   patch_2d%edges%vertex_idx(il_e,ib_e,2), patch_2d%edges%vertex_blk(il_e,ib_e,2)) %x
!                 write(0,*) " ", (patch_2d%verts%cartesian( &
!                   patch_2d%edges%vertex_idx(il_e,ib_e,2), patch_2d%edges%vertex_blk(il_e,ib_e,2))%x + &
!                   patch_2d%verts%cartesian( &
!                   patch_2d%edges%vertex_idx(il_e,ib_e,2), patch_2d%edges%vertex_blk(il_e,ib_e,2))%x) * 0.5_wp
! 
!                 write(0,*) " "
!                 write(0,*) " ", patch_2d%cells%cartesian_center(cell_index,blockNo)%x
!                 write(0,*) " ", patch_2d%edges%cartesian_center(il_e,ib_e)%x - &
!                   & patch_2d%cells%cartesian_center(cell_index,blockNo)%x
! 
!                 write(0,*) " "
!                 write(0,*) "   ", vn_e(il_e,level,ib_e)
! 
!                 p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x&
!                   & + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,ie)%x&
!                   & * vn_e(il_e,level,ib_e)*patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
!                 write(0,*) "   velocity:", p_vn_c(cell_index,level,blockNo)%x
!               ENDIF
!             END DO
!           ENDIF

        END DO level_loop
      END DO  cell_idx_loop
      
    END DO ! blockNo = all_cells%start_block, all_cells%end_block

!     CALL dbg_print('x(1)', p_vn_c(:,:,:)%x(1), &
!       & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     CALL dbg_print('x(2)', p_vn_c(:,:,:)%x(2), &
!       & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    CALL dbg_print('x(3)', p_vn_c(:,:,:)%x(3), &
      & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)

  END SUBROUTINE map_edges2cell_no_height_3d
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !: This sbr is used for kinetic energy calculation
  SUBROUTINE map_edges2cell_no_height_3d_onTriangles( patch_3d, vn_e, operators_coefficients, p_vn_c, &
    & opt_startLevel, opt_endLevel, subset_range)

    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    TYPE(t_cartesian_coordinates)              :: p_vn_c(:,:,:)  ! output (nproma,n_zlev,alloc_cell_blocks)
    ! intent(inout) for nag compiler
    INTEGER, INTENT(in), OPTIONAL :: opt_startLevel,  opt_endLevel      ! optional vertical start level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: edge_1_index, edge_1_block, edge_2_index, edge_2_block, edge_3_index, edge_3_block
    INTEGER :: cell_index, blockNo, level
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2d%cells%ALL
    ENDIF
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_startLevel) ) THEN
      startLevel = opt_startLevel
    ELSE
      startLevel = 1
    END IF
    IF ( PRESENT(opt_endLevel) ) THEN
      endLevel = opt_endLevel
    ELSE
      endLevel = n_zlev
    END IF

    !Calculation of Pv in cartesian coordinates
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, &
!ICON_OMP edge_1_index, edge_1_block, edge_2_index, edge_2_block, edge_3_index, edge_3_block,  &
!ICON_OMP level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      DO cell_index = start_cell_index, end_cell_index

        edge_1_index = patch_2d%cells%edge_idx(cell_index,blockNo,1)
        edge_1_block = patch_2d%cells%edge_blk(cell_index,blockNo,1)
        edge_2_index = patch_2d%cells%edge_idx(cell_index,blockNo,2)
        edge_2_block = patch_2d%cells%edge_blk(cell_index,blockNo,2)
        edge_3_index = patch_2d%cells%edge_idx(cell_index,blockNo,3)
        edge_3_block = patch_2d%cells%edge_blk(cell_index,blockNo,3)


        DO level = startLevel, MIN(patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo), endLevel)
          p_vn_c(cell_index,level,blockNo)%x =                                            &
            & (  operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,1)%x  &
            &      * vn_e(edge_1_index,level,edge_1_block)                                &
            &*patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_index,level,edge_1_block)       &
            &  + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,2)%x  &
            &      * vn_e(edge_2_index,level,edge_2_block)                                &
            &*patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_index,level,edge_2_block)       &            
            &  + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,3)%x  &
            &       * vn_e(edge_3_index,level,edge_3_block)                               &
            &*patch_3d%p_patch_1d(1)%prism_thick_e(edge_3_index,level,edge_3_block))      &
            & / (operators_coefficients%fixed_vol_norm(cell_index,level,blockNo)          &
            & *patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo))
        END DO
          
      END DO

    END DO ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2cell_no_height_3d_onTriangles
  !-----------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  !This sbr is used in continuity equation for the calculation of vertical velocity, and
  !in setting up he right hand side for the free-surface equation.
  !With an additional (tracer)-scalar it is used optionally in tracer transport.
  !It has a 2d-version which is used in left-hand-side of free-surface equation.
  SUBROUTINE map_edges2edges_viacell_3d_mlev_const_z( patch_3d, vn_e, operators_coefficients, out_vn_e)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, blockNo, level, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    IF ( patch_2d%cells%max_connectivity == 3 .and. fast_performance_level > 10) THEN
      CALL map_edges2edges_viacell_3d_mlev_constZ_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e)
      RETURN
    ENDIF
    !-----------------------------------------------------------------------
    edges_inDomain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    
    DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
      CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
      
     out_vn_e(:,:,blockNo) = 0.0_wp
      
      level_loop_e: DO level = startLevel, endLevel
        edge_idx_loop: DO je =  start_edge_index, end_edge_index
          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN
            out_vn_e(je,level,blockNo)= 0.0_wp
            
            ictr = 0
            il_c = patch_2d%edges%cell_idx(je,blockNo,1)
            ib_c = patch_2d%edges%cell_blk(je,blockNo,1)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,level,ib_c)
            
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
              IF (il_e > 0) THEN
                thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)

                out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
                  & +vn_e(il_e,level,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                  & *thick_edge
              ENDIF
              
            END DO
            
            ictr = no_primal_edges
            il_c = patch_2d%edges%cell_idx(je,blockNo,2)
            ib_c = patch_2d%edges%cell_blk(je,blockNo,2)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,level,ib_c)
            
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
              IF (il_e > 0) THEN
                thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)

                out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
                  & +vn_e(il_e,level,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                  & *thick_edge
              ENDIF
              
            END DO
          ENDIF
        END DO edge_idx_loop
      END DO level_loop_e
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    
  END SUBROUTINE map_edges2edges_viacell_3d_mlev_const_z
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! the map_edges2edges_viacell_3d_mlev_constZ optimized for triangles
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_3d_mlev_constZ_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, blockNo, level

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
    INTEGER :: edge_11_block, edge_12_block, edge_13_block
    INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
    INTEGER :: edge_21_block, edge_22_block, edge_23_block

    REAL(wp), POINTER :: coeffs(:,:,:,:)
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')

    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    coeffs => operators_coefficients%edge2edge_viacell_coeff
    !-----------------------------------------------------------------------

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_edge_index, end_edge_index, je, cell_1_index, cell_1_block, &
!ICON_OMP   cell_2_index, cell_2_block, edge_11_index, edge_12_index, edge_13_index, &
!ICON_OMP  edge_11_block, edge_12_block, edge_13_block, edge_21_index, edge_22_index, &
!ICON_OMP  edge_23_index, edge_21_block, edge_22_block, edge_23_block, level)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      out_vn_e(:, :, blockNo) = 0.0_wp
      DO je =  start_edge_index, end_edge_index

        IF (patch_3d%p_patch_1d(1)%dolic_e(je,blockNo) < 1) CYCLE

        cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
        cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)

        edge_11_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
        edge_12_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
        edge_13_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
        edge_11_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
        edge_12_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
        edge_13_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)

        edge_21_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
        edge_22_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
        edge_23_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
        edge_21_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
        edge_22_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
        edge_23_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)

        ! levels
        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)

          out_vn_e(je, level, blockNo) =  &
            & (  vn_e(edge_11_index, level, edge_11_block) * coeffs(je, level, blockNo, 1)               &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_11_index, level, edge_11_block)  +  &
            & vn_e(edge_12_index, level, edge_12_block) * coeffs(je, level, blockNo, 2)               &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_12_index, level, edge_12_block)  +  &
            & vn_e(edge_13_index, level, edge_13_block) * coeffs(je, level, blockNo, 3)               &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_13_index, level, edge_13_block)     &
            & )                               &
            & + &
            & (  vn_e(edge_21_index, level, edge_21_block) * coeffs(je, level, blockNo, 4)               &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_21_index, level, edge_21_block)  +  &
            & vn_e(edge_22_index, level, edge_22_block) * coeffs(je, level, blockNo, 5)               &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_22_index, level, edge_22_block)  +  &
            & vn_e(edge_23_index, level, edge_23_block) * coeffs(je, level, blockNo, 6)               &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_23_index, level, edge_23_block)     &
            & ) 

        END DO ! levels

      END DO

    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

  END SUBROUTINE map_edges2edges_viacell_3d_mlev_constZ_onTriangles
  !-----------------------------------------------------------------------------
 
 
  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  !Same sbr as above, now with an additional scalar (e.g. tracer). It is used optionally in tracer transport.
  SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs( patch_3d, vn_e, operators_coefficients, out_vn_e, scalar)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)       :: scalar(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, blockNo, level, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER :: patch_2d
    
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    IF ( patch_2d%cells%max_connectivity == 3) THEN
      IF (fast_performance_level > 10 ) THEN
        CALL map_edges2edges_viacell_3d_mlev_constZs_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e, scalar)
        RETURN
      ENDIF
    ENDIF
    
    !-----------------------------------------------------------------------
    edges_inDomain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    
    DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
      CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
      
      out_vn_e(:,:,blockNo) = 0.0_wp
      
      level_loop_e2: DO level = startLevel, endLevel
        edge_idx_loop2: DO je =  start_edge_index, end_edge_index
          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN
            out_vn_e(je,level,blockNo)= 0.0_wp
            
            ictr = 0
            il_c        = patch_2d%edges%cell_idx(je,blockNo,1)
            ib_c        = patch_2d%edges%cell_blk(je,blockNo,1)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,level,ib_c)
            scalar_cell = scalar(il_c,level,ib_c)
            
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
              
              out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
                & +vn_e(il_e,level,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                & *thick_edge*scalar_cell
              
            END DO
            
            ictr        = no_primal_edges
            il_c        = patch_2d%edges%cell_idx(je,blockNo,2)
            ib_c        = patch_2d%edges%cell_blk(je,blockNo,2)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,level,ib_c)
            scalar_cell = scalar(il_c,level,ib_c)
            
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
              
              out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
                & +vn_e(il_e,level,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                & *thick_edge*scalar_cell
              
            END DO
          ENDIF
        END DO edge_idx_loop2
      END DO level_loop_e2
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    
  END SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! the map_edges2edges_viacell_3d_mlev_constZs optimized for triangles
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e,scalar)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    REAL(wp), INTENT(in)                 :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)   :: operators_coefficients
    REAL(wp), INTENT(inout)              :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp), INTENT(in)                 :: scalar(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !Local variables
    INTEGER :: startLevel, endLevel
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp), POINTER :: coeffs(:,:,:,:)

    ! omp private
    ! defined in a struct for omp clearness
    TYPE omp_local_private
    
      INTEGER :: start_edge_index, end_edge_index
      INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
      INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
      INTEGER :: edge_11_block, edge_12_block, edge_13_block
      INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
      INTEGER :: edge_21_block, edge_22_block, edge_23_block
      
    END TYPE omp_local_private
    
!     TYPE(omp_local_private) :: omp_this 
    INTEGER :: je, blockNo, level
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
    INTEGER :: edge_11_block, edge_12_block, edge_13_block
    INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
    INTEGER :: edge_21_block, edge_22_block, edge_23_block
   !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_in_domain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = n_zlev
    coeffs => operators_coefficients%edge2edge_viacell_coeff
    !-----------------------------------------------------------------------
    
!ICON_OMP_PARALLEL_DO PRIVATE( je, level, start_edge_index, end_edge_index, &
!ICON_OMP cell_1_index, cell_2_index, cell_1_block, cell_2_block, &
!ICON_OMP edge_11_index, edge_12_index, edge_13_index, &
!ICON_OMP edge_21_index, edge_22_index, edge_23_index, &
!ICON_OMP edge_21_block, edge_22_block, edge_23_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      DO je =  start_edge_index, end_edge_index
        
        cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
        cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)
        
        edge_11_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
        edge_12_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
        edge_13_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
        edge_11_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
        edge_12_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
        edge_13_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
        
        edge_21_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
        edge_22_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
        edge_23_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
        edge_21_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
        edge_22_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
        edge_23_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)
                
        ! levels
        DO level = startLevel, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          
          out_vn_e(je, level, blockNo) =  &
            & (  vn_e(edge_11_index, level, edge_11_block) * coeffs(je, level, blockNo, 1)      &
            &    * patch_3d%p_patch_1d(1)%prism_thick_e(edge_11_index, level, edge_11_block)    &
            &  + vn_e(edge_12_index, level, edge_12_block) * coeffs(je, level, blockNo, 2)      &
            &    * patch_3d%p_patch_1d(1)%prism_thick_e(edge_12_index, level, edge_12_block)    &
            &  + vn_e(edge_13_index, level, edge_13_block) * coeffs(je, level, blockNo, 3)      &
            &    * patch_3d%p_patch_1d(1)%prism_thick_e(edge_13_index, level, edge_13_block)    &
            & ) * scalar(cell_1_index, level, cell_1_block)                                     &
            & + &
            & (  vn_e(edge_21_index, level, edge_21_block) * coeffs(je, level, blockNo, 4)           &
            &   * patch_3d%p_patch_1d(1)%prism_thick_e(edge_21_index, level, edge_21_block)  +  &
            &  vn_e(edge_22_index, level, edge_22_block) * coeffs(je, level, blockNo, 5)             &
            &  * patch_3d%p_patch_1d(1)%prism_thick_e(edge_22_index, level, edge_22_block)  +   &
            &  vn_e(edge_23_index, level, edge_23_block) * coeffs(je, level, blockNo, 6)             &
            & * patch_3d%p_patch_1d(1)%prism_thick_e(edge_23_index, level, edge_23_block)       &
            & ) * scalar(cell_2_index, level, cell_2_block)
          
        END DO !levels
        
      END DO
      
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO
     
  END SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs_onTriangles
  !-----------------------------------------------------------------------------
  


  
  !-----------------------------------------------------------------------------
!<Optimize:inUse>
  !This sbr is used in left-hand-side of free-surface equation for certain solver configurations.
  SUBROUTINE map_edges2edges_viacell_2D_constZ( patch_3d, vn_e, operators_coefficients, out_vn_e)!&   subset_range)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, blockNo, level, ie
    !REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell!, thick_frac
    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER :: patch_2d
    
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    IF ( patch_2d%cells%max_connectivity == 3 .and. fast_performance_level > 10 ) THEN
      CALL map_edges2edges_viacell_2D_constZ_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e )
      RETURN
    ENDIF
    !-----------------------------------------------------------------------
    
    edges_inDomain => patch_2d%edges%in_domain
    startLevel = 1
    endLevel = startLevel!n_zlev

    DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
      CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
      
      out_vn_e(:,blockNo) = 0.0_wp
!       
      edge_idx_loop: DO je =  start_edge_index, end_edge_index
        endLevel=patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        level_loop_e: DO level = startLevel, endLevel
        
          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN
            
            ictr = 0
            il_c = patch_2d%edges%cell_idx(je,blockNo,1)
            ib_c = patch_2d%edges%cell_blk(je,blockNo,1)
            !thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
            
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
              IF (il_e > 0) THEN
                thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)

                out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
                  & +vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                  & *thick_edge
              ENDIF
              
            END DO
            
            ictr = no_primal_edges
            il_c = patch_2d%edges%cell_idx(je,blockNo,2)
            ib_c = patch_2d%edges%cell_blk(je,blockNo,2)
            !thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
            
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
              IF (il_e > 0) THEN
                thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)

                out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
                  & +vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                  & *thick_edge
              ENDIF
              
            END DO
          ENDIF
        END DO  level_loop_e
      END DO edge_idx_loop
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    
  END SUBROUTINE map_edges2edges_viacell_2D_constZ
  !-------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_2D_constZ_onTriangles( patch_3d, in_vn_e, operators_coefficients, out_vn_e )!&   subset_range)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                       :: in_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    
    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, blockNo, start_edge_index, end_edge_index, level
    
    REAL(wp), POINTER :: all_coeffs(:,:,:)
    
    TYPE(t_subset_range), POINTER :: edges_indomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
!     IF (no_primal_edges /= 3) &
!       & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_indomain    => patch_2d%edges%in_domain
    all_coeffs        => operators_coefficients%edge2edge_viacell_coeff_all
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, level, cell_1_index, cell_1_block, &
!ICON_OMP  cell_2_index, cell_2_block, edge_1_1_index, edge_1_2_index, edge_1_3_index, &
!ICON_OMP  edge_1_1_block, edge_1_2_block, edge_1_3_block, edge_2_1_index, edge_2_2_index, &
!ICON_OMP  edge_2_3_index, edge_2_1_block, edge_2_2_block, edge_2_3_block)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_indomain%start_block, edges_indomain%end_block
      CALL get_index_range(edges_indomain, blockNo, start_edge_index, end_edge_index)
      
      DO je = start_edge_index, end_edge_index
        
        out_vn_e(je,blockNo) = 0.0_wp
        
        DO level = 1, MIN(1,patch_3d%p_patch_1d(1)%dolic_e(je,blockNo))
!         IF (patch_3d%lsm_e(je,1,blockNo) == sea) THEN ! this if should be removed
          
          ! get the two cells of the edge
          cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
          cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
          cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
          cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)
          
          ! get the six edges of the two cells
          edge_1_1_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
          edge_1_2_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
          edge_1_3_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
          edge_2_1_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
          edge_2_2_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
          edge_2_3_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
          edge_1_1_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
          edge_1_2_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
          edge_1_3_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
          edge_2_1_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
          edge_2_2_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
          edge_2_3_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)
          
          out_vn_e(je,blockNo) = &
            & in_vn_e(edge_1_1_index, edge_1_1_block) * all_coeffs(1, je, blockNo) + &
            & in_vn_e(edge_1_2_index, edge_1_2_block) * all_coeffs(2, je, blockNo) + &
            & in_vn_e(edge_1_3_index, edge_1_3_block) * all_coeffs(3, je, blockNo) + &
            & in_vn_e(edge_2_1_index, edge_2_1_block) * all_coeffs(4, je, blockNo) + &
            & in_vn_e(edge_2_2_index, edge_2_2_block) * all_coeffs(5, je, blockNo) + &
            & in_vn_e(edge_2_3_index, edge_2_3_block) * all_coeffs(6, je, blockNo)
          
        END DO
      END DO
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2edges_viacell_2D_constZ_onTriangles
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_2d_constZ_onTriangles_sp( patch_3d, in_vn_e, operators_coefficients, out_vn_e )!&   subset_range)

    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(sp), INTENT(in)                       :: in_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_solvercoeff_singleprecision), INTENT(in)  :: operators_coefficients
    REAL(sp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, blockNo, start_edge_index, end_edge_index, level

    REAL(sp), POINTER :: all_coeffs(:,:,:)

    TYPE(t_subset_range), POINTER :: edges_indomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
!     IF (no_primal_edges /= 3) &
!       & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_indomain    => patch_2d%edges%in_domain
    all_coeffs        => operators_coefficients%edge2edge_viacell_coeff_all
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, level, cell_1_index, cell_1_block, &
!ICON_OMP  cell_2_index, cell_2_block, edge_1_1_index, edge_1_2_index, edge_1_3_index, &
!ICON_OMP  edge_1_1_block, edge_1_2_block, edge_1_3_block, edge_2_1_index, edge_2_2_index, &
!ICON_OMP  edge_2_3_index, edge_2_1_block, edge_2_2_block, edge_2_3_block)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_indomain%start_block, edges_indomain%end_block
      CALL get_index_range(edges_indomain, blockNo, start_edge_index, end_edge_index)

      DO je = start_edge_index, end_edge_index

        out_vn_e(je,blockNo) = 0.0_wp

        DO level = 1, MIN(1,patch_3d%p_patch_1d(1)%dolic_e(je,blockNo))

          ! get the two cells of the edge
          cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
          cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
          cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
          cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)

          ! get the six edges of the two cells
          edge_1_1_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
          edge_1_2_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
          edge_1_3_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
          edge_2_1_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
          edge_2_2_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
          edge_2_3_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
          edge_1_1_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
          edge_1_2_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
          edge_1_3_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
          edge_2_1_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
          edge_2_2_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
          edge_2_3_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)

          out_vn_e(je,blockNo) = &
            & in_vn_e(edge_1_1_index, edge_1_1_block) * all_coeffs(1, je, blockNo) + &
            & in_vn_e(edge_1_2_index, edge_1_2_block) * all_coeffs(2, je, blockNo) + &
            & in_vn_e(edge_1_3_index, edge_1_3_block) * all_coeffs(3, je, blockNo) + &
            & in_vn_e(edge_2_1_index, edge_2_1_block) * all_coeffs(4, je, blockNo) + &
            & in_vn_e(edge_2_2_index, edge_2_2_block) * all_coeffs(5, je, blockNo) + &
            & in_vn_e(edge_2_3_index, edge_2_3_block) * all_coeffs(6, je, blockNo)

        ENDDO
      END DO
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2edges_viacell_2d_constZ_onTriangles_sp
  !----------------------------------------------------------------------------- 
  
  

  !-----------------------------------------------------------------------------
!<Optimize:inUse>
  !This sbr is used in left-hand-side of free-surface equation for certain solver configurations.
  SUBROUTINE map_edges2edges_viacell_2D_per_level( patch_3d, vn_e, operators_coefficients, out_vn_e, level)!&   subset_range)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                       :: vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER, INTENT(IN)                        :: level
    !Local variables
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, blockNo, ie
    REAL(wp) :: thick_edge, thick_cell
    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER :: patch_2d
    
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
!     IF ( patch_2d%cells%max_connectivity == 3 .and. fast_performance_level > 10 ) THEN
!       CALL map_edges2edges_viacell_2D_per_level_onTriangles( patch_3d, vn_e, operators_coefficients, out_vn_e, level)!&   subset_range)
!       RETURN
!     ENDIF
    !-----------------------------------------------------------------------
    
    edges_inDomain => patch_2d%edges%in_domain
  
    DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
      CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
      
      out_vn_e(:,blockNo) = 0.0_wp

      edge_idx_loop: DO je =  start_edge_index, end_edge_index
        
          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN
            
            ictr = 0
            il_c = patch_2d%edges%cell_idx(je,blockNo,1)
            ib_c = patch_2d%edges%cell_blk(je,blockNo,1)
            thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
            
            
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
              !IF(il_e==je.and.ib_e==blockNo)THEN
              out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
                & +vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                & *(thick_edge/thick_cell)
              !ENDIF
            END DO
            
            ictr = no_primal_edges
            il_c = patch_2d%edges%cell_idx(je,blockNo,2)
            ib_c = patch_2d%edges%cell_blk(je,blockNo,2)
            thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)            
            
            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
              !IF(il_e==je.and.ib_e==blockNo)THEN              
              out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
                & +vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr)&
                & *(thick_edge/thick_cell)
              !ENDIF
            END DO
          ENDIF
      END DO edge_idx_loop
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
    
  END SUBROUTINE map_edges2edges_viacell_2D_per_level
  !-------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_2D_per_level_onTriangles( patch_3d, in_vn_e, operators_coefficients, out_vn_e, level)!&   subset_range)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(wp), INTENT(in)                       :: in_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER, INTENT(IN)                        :: level
    
    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, blockNo, start_edge_index, end_edge_index
    
    REAL(wp), POINTER :: all_coeffs(:,:,:)
    TYPE(t_subset_range), POINTER :: edges_indomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_indomain    => patch_2d%edges%in_domain
    all_coeffs        => operators_coefficients%edge2edge_viacell_coeff_all
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, cell_1_index, cell_1_block, &
!ICON_OMP   cell_2_index, cell_2_block, edge_1_1_index, edge_1_2_index, edge_1_3_index, &
!ICON_OMP  edge_1_1_block, edge_1_2_block, edge_1_3_block, edge_2_1_index, edge_2_2_index, &
!ICON_OMP  edge_2_3_index, edge_2_1_block, edge_2_2_block, edge_2_3_block)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_indomain%start_block, edges_indomain%end_block
      CALL get_index_range(edges_indomain, blockNo, start_edge_index, end_edge_index)
      
      DO je = start_edge_index, end_edge_index
        
        out_vn_e(je,blockNo) = 0.0_wp
        
        IF (patch_3d%lsm_e(je,1,blockNo) == sea) THEN ! this if should be removed
          
          ! get the two cells of the edge
          cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
          cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
          cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
          cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)
          
          ! get the six edges of the two cells
          edge_1_1_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
          edge_1_2_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
          edge_1_3_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
          edge_2_1_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
          edge_2_2_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
          edge_2_3_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
          edge_1_1_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
          edge_1_2_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
          edge_1_3_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
          edge_2_1_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
          edge_2_2_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
          edge_2_3_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)
          !not that all_coeffs already includes the edge thickness
          out_vn_e(je,blockNo) = &
            & in_vn_e(edge_1_1_index, edge_1_1_block) * all_coeffs(1, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_1_index,level, edge_1_1_block)& 
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block)     &
            &+ &
            & in_vn_e(edge_1_2_index, edge_1_2_block) * all_coeffs(2, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_2_index,level, edge_1_2_block)&
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block)     &
            &+ &
            & in_vn_e(edge_1_3_index, edge_1_3_block) * all_coeffs(3, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_3_index,level, edge_1_3_block)&
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block)     &
            &+ &
            & in_vn_e(edge_2_1_index, edge_2_1_block) * all_coeffs(4, je, blockNo)       &
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block)     &
            &+ &
            & in_vn_e(edge_2_2_index, edge_2_2_block) * all_coeffs(5, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_2_index,level, edge_2_2_block)&
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block)     &
            &+ &
            & in_vn_e(edge_2_3_index, edge_2_3_block) * all_coeffs(6, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_3_index,level, edge_2_3_block)&
            & * patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block)     
          
        ENDIF
      END DO
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2edges_viacell_2D_per_level_onTriangles
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE map_edges2edges_viacell_2D_per_level_onTriangles_sp( patch_3d, in_vn_e, operators_coefficients, out_vn_e, level)!&   subset_range)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
    REAL(sp), INTENT(in)                       :: in_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e) 
    TYPE(t_solvercoeff_singleprecision), INTENT(in)  :: operators_coefficients
    REAL(sp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
    INTEGER, INTENT(IN)                        :: level
    
    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, blockNo, start_edge_index, end_edge_index
    
    REAL(sp), POINTER :: all_coeffs(:,:,:)
    !REAL(wp) :: thick_edge, thick_cell
    TYPE(t_subset_range), POINTER :: edges_indomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      & CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_indomain    => patch_2d%edges%in_domain
    all_coeffs        => operators_coefficients%edge2edge_viacell_coeff_all
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, cell_1_index, cell_1_block, &
!ICON_OMP   cell_2_index, cell_2_block, edge_1_1_index, edge_1_2_index, edge_1_3_index, &
!ICON_OMP  edge_1_1_block, edge_1_2_block, edge_1_3_block, edge_2_1_index, edge_2_2_index, &
!ICON_OMP  edge_2_3_index, edge_2_1_block, edge_2_2_block, edge_2_3_block)  ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_indomain%start_block, edges_indomain%end_block
      CALL get_index_range(edges_indomain, blockNo, start_edge_index, end_edge_index)
      
      DO je = start_edge_index, end_edge_index
        
        out_vn_e(je,blockNo) = 0.0_wp
        
        IF (patch_3d%lsm_e(je,1,blockNo) == sea) THEN ! this if should be removed
          
          ! get the two cells of the edge
          cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
          cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
          cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
          cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)
          
          ! get the six edges of the two cells
          edge_1_1_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 1)
          edge_1_2_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 2)
          edge_1_3_index = patch_2d%cells%edge_idx(cell_1_index, cell_1_block, 3)
          edge_2_1_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 1)
          edge_2_2_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 2)
          edge_2_3_index = patch_2d%cells%edge_idx(cell_2_index, cell_2_block, 3)
          edge_1_1_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 1)
          edge_1_2_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 2)
          edge_1_3_block = patch_2d%cells%edge_blk(cell_1_index, cell_1_block, 3)
          edge_2_1_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 1)
          edge_2_2_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 2)
          edge_2_3_block = patch_2d%cells%edge_blk(cell_2_index, cell_2_block, 3)
          !not that all_coeffs already includes the edge thickness
          out_vn_e(je,blockNo) = &
            & in_vn_e(edge_1_1_index, edge_1_1_block) * all_coeffs(1, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_1_index,level, edge_1_1_block)& 
            &* (patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block))     &
            &+ &
            & in_vn_e(edge_1_2_index, edge_1_2_block) * all_coeffs(2, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_2_index,level, edge_1_2_block)&
            & * (patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block))     &
            &+ &
            & in_vn_e(edge_1_3_index, edge_1_3_block) * all_coeffs(3, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_1_3_index,level, edge_1_3_block)&
            & * (patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_1_index,level,cell_1_block))     &
            &+ &
            & in_vn_e(edge_2_1_index, edge_2_1_block) * all_coeffs(4, je, blockNo)       &
            & * (patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block))     &
            &+ &
            & in_vn_e(edge_2_2_index, edge_2_2_block) * all_coeffs(5, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_2_index,level, edge_2_2_block)&
            &/(patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block))     &
            &+ &
            & in_vn_e(edge_2_3_index, edge_2_3_block) * all_coeffs(6, je, blockNo)       &
            !&*(patch_3d%p_patch_1d(1)%prism_thick_e(edge_2_3_index,level, edge_2_3_block)&
            & * (patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_2_index,level,cell_2_block))     
          
        ENDIF
      END DO
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_edges2edges_viacell_2D_per_level_onTriangles_sp
  !-----------------------------------------------------------------------------
  
  
  
  !-----------------------------------------------------------------------------
  SUBROUTINE map_edges2edges_viavert_3d(patch_3d, vn, p_vn_dual,operators_coefficients, vort_flux)
    TYPE(t_patch_3d ),TARGET,INTENT(in)        :: patch_3d
    REAL(wp), INTENT(inout)                    :: vn(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(inout)  :: p_vn_dual(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in)          :: operators_coefficients
    REAL(wp), INTENT(inout)                    :: vort_flux(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
    
    !Local variables
    INTEGER :: startLevel, endLevel     ! vertical start and end level
    INTEGER :: je, level, blockNo
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_inDomain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    startLevel    = 1
    endLevel    = n_zlev
    
    DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
      CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
      
      level_loop: DO level = startLevel, endLevel

        vort_flux(:,:,blockNo) = 0.0_wp
        edge_idx_loop: DO je =  start_edge_index, end_edge_index
          
          IF (patch_3d%lsm_e(je,level,blockNo) == sea) THEN
            
            DO neighbor=1,2
              IF(neighbor==1) ictr = 0
              IF(neighbor==2) ictr = no_dual_edges
              
              il_v = patch_2d%edges%vertex_idx(je,blockNo,neighbor)
              ib_v = patch_2d%edges%vertex_blk(je,blockNo,neighbor)
              
              DO vertex_edge=1, patch_2d%verts%num_edges(il_v,ib_v)!no_dual_cell_edges
                
                ictr =ictr+1
                
                il_e = patch_2d%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = patch_2d%verts%edge_blk(il_v,ib_v,vertex_edge)
                
                vort_flux(je,level,blockNo) =  vort_flux(je,level,blockNo)+vn(il_e,level,ib_e)&
                  & *operators_coefficients%edge2edge_viavert_coeff(je,level,blockNo,ictr)
              END DO
            END DO
            
          ELSE
            vort_flux(je,level,blockNo)= 0.0_wp
          ENDIF ! (v_base%lsm_e(je,level,blockNo) <= sea_boundary)
        END DO edge_idx_loop
      END DO level_loop
    END DO ! blockNo = edges_inDomain%start_block, edges_inDomain%end_block
    
  END SUBROUTINE map_edges2edges_viavert_3d
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !<Optimize:inUse>
  SUBROUTINE map_cell2edges_3d_mlevels( patch_3d, p_vn_c, ptp_vn, operators_coefficients,&
    & opt_startLevel, opt_endLevel, subset_range )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:,:)    ! input vector (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)                      :: ptp_vn(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff)                     :: operators_coefficients
    INTEGER, INTENT(in), OPTIONAL :: opt_startLevel        ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_endLevel        ! optional vertical end level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    
    !Local variables
    INTEGER :: startLevel, endLevel
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, blockNo, level
    INTEGER :: cell_1_index,cell_1_block, cell_2_index,cell_2_block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    edges_in_domain   => patch_2d%edges%in_domain
    
    
    ! check optional arguments
    IF ( PRESENT(opt_startLevel) ) THEN
      startLevel = opt_startLevel
    ELSE
      startLevel = 1
    END IF
    IF ( PRESENT(opt_endLevel) ) THEN
      endLevel = opt_endLevel
    ELSE
      endLevel = n_zlev
    END IF
    
    ! calculation of transposed P^TPv from Pv (incart coord)
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index,je, level,  &
!ICON_OMP cell_1_index, cell_1_block, cell_2_index, cell_2_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      ptp_vn(:,:,blockNo) = 0.0_wp
      
      DO je = start_edge_index, end_edge_index
        !Get indices of two adjacent triangles
        cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
        cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
        cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
        cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)
        DO level = 1, patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
          
          ptp_vn(je,level,blockNo) =&
            & DOT_PRODUCT(p_vn_c(cell_1_index,level,cell_1_block)%x,&
            & operators_coefficients%edge2cell_coeff_cc_t(je,level,blockNo,1)%x)&
            & +DOT_PRODUCT(p_vn_c(cell_2_index,level,cell_2_block)%x,&
            & operators_coefficients%edge2cell_coeff_cc_t(je,level,blockNo,2)%x)
          
        END DO 
      END DO
      
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO    

    ! sync the result if necessary
    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2d, ptp_vn)
    ENDIF
    
  END SUBROUTINE map_cell2edges_3d_mlevels
  !-----------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !<Optimize:inUse:done>
  SUBROUTINE map_cell2edges_3d_1level( patch_3d, p_vn_c, ptp_vn,operators_coefficients, level, subset_range )
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:)    ! input vector (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)                      :: ptp_vn(:,:)    ! output vector (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff)                     :: operators_coefficients
    INTEGER, INTENT(in) :: level          ! vertical level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    
    !Local variables
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, blockNo
    INTEGER :: cell_1_index,cell_1_block, cell_2_index,cell_2_block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    !-----------------------------------------------------------------------
    edges_in_domain   => patch_2d%edges%in_domain
    
    ! calculation of transposed P^TPv from Pv (incart coord)
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index, end_edge_index, je, &
!ICON_OMP  cell_1_index,cell_1_block, cell_2_index,cell_2_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      ptp_vn(:,blockNo) = 0.0_wp
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      DO je =  start_edge_index, end_edge_index
        
        IF (patch_3d%p_patch_1d(1)%dolic_e(je, blockNo) > 0) THEN
          !Get indices of two adjacent triangles
          cell_1_index = patch_2d%edges%cell_idx(je,blockNo,1)
          cell_1_block = patch_2d%edges%cell_blk(je,blockNo,1)
          cell_2_index = patch_2d%edges%cell_idx(je,blockNo,2)
          cell_2_block = patch_2d%edges%cell_blk(je,blockNo,2)
          ptp_vn(je,blockNo) = &
            & DOT_PRODUCT(p_vn_c(cell_1_index,cell_1_block)%x,               &
            & operators_coefficients%edge2cell_coeff_cc_t(je,level,blockNo,1)%x)  &
            & +DOT_PRODUCT(p_vn_c(cell_2_index,cell_2_block)%x,              &
            & operators_coefficients%edge2cell_coeff_cc_t(je,level,blockNo,2)%x)

        ENDIF
        
      END DO
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO    
    ! sync the result if necessary
    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2d, ptp_vn)
    ENDIF
    
  END SUBROUTINE map_cell2edges_3d_1level
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE maps for a fluid column a scalar value from the top/bottom of a 3D prism to the central level of the prism.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE map_vec_prismtop2center_on_block(patch_3d, vec_top, p_op_coeff, vec_center, &
    & blockNo, start_cell_index, end_cell_index)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_cartesian_coordinates), INTENT(in)        :: vec_top(nproma, n_zlev)
    TYPE(t_operator_coeff), INTENT(in)               :: p_op_coeff
    INTEGER, INTENT(in)                              :: blockNo, start_cell_index, end_cell_index
    TYPE(t_cartesian_coordinates), INTENT(inout)     :: vec_center(nproma, n_zlev,blockNo) ! out
    
    !Local variables
    INTEGER :: level, jc!,jb
    INTEGER :: start_level, end_level
    REAL(wp), POINTER :: prism_center_distance(:,:), prism_thick(:,:)
    !-------------------------------------------------------------------------------
    start_level = 1
!     vec_center(1:nproma,1:n_zlev,blockNo)%x(1)=0.0_wp
!     vec_center(1:nproma,1:n_zlev,blockNo)%x(2)=0.0_wp
!     vec_center(1:nproma,1:n_zlev,blockNo)%x(3)=0.0_wp
   !-------------------------------------------------------------------------------    
    ! these do not include the height
    prism_center_distance => patch_3D%p_patch_1D(1)%constantPrismCenters_Zdistance  (:,:,blockNo)
    prism_thick           => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,blockNo)

    DO jc = start_cell_index, end_cell_index
      end_level  = patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)
!       IF ( end_level >=min_dolic ) THEN
        DO level = start_level, end_level-1
          vec_center(jc,level,blockNo)%x &
          & = (prism_center_distance(jc,level)   * vec_top(jc,level)%x    &
          & +  prism_center_distance(jc,level+1) * vec_top(jc,level+1)%x) &
          & / (2.0_wp*prism_thick(jc,level))
              
        END DO          
!       ENDIF
    END DO
!     CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vec_center(:,:,:)%x(1))
!     CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vec_center(:,:,:)%x(2))
!     CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vec_center(:,:,:)%x(3))
  END SUBROUTINE map_vec_prismtop2center_on_block
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE maps for a fluid column a scalar value from the top/bottom of a 3D prism to the central level of the prism.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE map_scalar_prismtop2center(patch_3d, scalar_top, p_op_coeff,scalar_center)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    REAL(wp)                                         :: scalar_top(nproma, n_zlev+1,patch_3D%p_patch_2d(1)%nblks_c)   
    TYPE(t_operator_coeff),INTENT(in)                :: p_op_coeff
    REAL(wp)                                         :: scalar_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_c)       
    
    !Local variables
    INTEGER :: level, blockNo, jc!,jb
    INTEGER :: start_cell_index, end_cell_index!, cell_index
    !INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level!, level
    TYPE(t_subset_range), POINTER :: cells_in_domain!, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D 
    REAL(wp), POINTER ::  prism_center_distance(:,:),prism_thick(:,:)
    ! INTEGER :: dolic
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    !edges_in_domain => patch_2D%edges%in_domain
    start_level = 1

   !-------------------------------------------------------------------------------  
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, level,&
!ICON_OMP prism_center_distance,prism_thick, jc) ICON_OMP_DEFAULT_SCHEDULE  
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

      ! these do not include the height
      ! if they are used for the GM-Redi then we should re-consider if height is needed
      prism_center_distance => patch_3D%p_patch_1D(1)%constantPrismCenters_Zdistance  (:,:,blockNo)
      prism_thick => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,blockNo)

      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
!        dolic  = patch_3D%p_patch_1d(1)%dolic_c(cell_index,blockNo)
!        IF ( dolic >=min_dolic ) THEN
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)-1
          scalar_center(jc,level,blockNo) &
          & = (prism_center_distance(jc,level)   * scalar_top(jc,level,blockNo)    &
          & +  prism_center_distance(jc,level+1) * scalar_top(jc,level+1,blockNo)) &
          & / (2.0_wp*prism_thick(jc,level))
              
        END DO
!        ENDIF
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO
   !CALL sync_patch_array(sync_c, patch_2D, scalar_center)
  END SUBROUTINE map_scalar_prismtop2center
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE maps within a fluid column a scalar value from the the central level of the prism to top/bottom of a 3D prism.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE map_scalar_center2prismtop(patch_3d, scalar_center, p_op_coeff,scalar_top)
    TYPE(t_patch_3d ),TARGET, INTENT(in   )          :: patch_3d
    REAL(wp), INTENT(in)                             :: scalar_center(:,:,:)   
    TYPE(t_operator_coeff),INTENT(in)                :: p_op_coeff
    REAL(wp), INTENT(out)                            :: scalar_top(:,:,:)       
    
    !Local variables
    INTEGER :: blockNo
    INTEGER :: start_cell_index, end_cell_index, cell_index
    !INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level
    TYPE(t_subset_range), POINTER :: cells_in_domain!, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D 
    REAL(wp), POINTER ::  prism_center_distance(:,:),prism_thick(:,:)
!     INTEGER :: dolic
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    !edges_in_domain => patch_2D%edges%in_domain
    start_level = 1

   !-------------------------------------------------------------------------------
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, level,&
!ICON_OMP prism_center_distance,prism_thick) ICON_OMP_DEFAULT_SCHEDULE     
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
    
      scalar_top(:,:,blockNo) = 0.0_wp
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      DO cell_index = start_cell_index, end_cell_index
        DO level = start_level, patch_3D%p_patch_1d(1)%dolic_c(cell_index,blockNo)-1
          scalar_top(cell_index,level+1,blockNo) &
          & = 0.5_wp*( scalar_center(cell_index,level,blockNo)    &
          & +          scalar_center(cell_index,level+1,blockNo))              
        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE map_scalar_center2prismtop
  !-------------------------------------------------------------------------

  ! !   !-----------------------------------------------------------------------------
!   !-----------------------------------------------------------------------------
!   SUBROUTINE map_edges2edges_viacell_3d_mlev( patch_3d, vn_e, operators_coefficients, out_vn_e,scalar, &
!     & opt_startLevel, opt_endLevel, subset_range)
!     
!     TYPE(t_patch_3d ),TARGET, INTENT(in)       :: patch_3d
!     REAL(wp), INTENT(in)                       :: vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
!     TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
!     REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,n_zlev,patch_3d%p_patch_2d(1)%nblks_e)
!     REAL(wp), INTENT(in), OPTIONAL :: scalar(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     INTEGER, INTENT(in), OPTIONAL :: opt_startLevel       ! optional vertical start level
!     INTEGER, INTENT(in), OPTIONAL :: opt_endLevel       ! optional vertical end level
!     TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
!     !Local variables
!     INTEGER :: startLevel, endLevel
!     INTEGER :: start_edge_index, end_edge_index
!     INTEGER :: il_e, ib_e, il_c, ib_c, ictr!, neighbor
!     INTEGER :: je, blockNo, level, ie
!     REAL(wp) :: scalar_cell
!     REAL(wp) :: thick_edge, thick_cell, thick_frac
!     TYPE(t_subset_range), POINTER :: edges_inDomain
!     TYPE(t_patch), POINTER :: patch_2d
!     !-----------------------------------------------------------------------
!     patch_2d   => patch_3d%p_patch_2d(1)
!     !-----------------------------------------------------------------------
!     IF ( PRESENT(subset_range) ) THEN
!       edges_inDomain => subset_range
!     ELSE
!       edges_inDomain => patch_2d%edges%in_domain
!     ENDIF
!     !-----------------------------------------------------------------------
!     IF ( PRESENT(opt_startLevel) ) THEN
!       startLevel = opt_startLevel
!     ELSE
!       startLevel = 1
!     END IF
!     IF ( PRESENT(opt_endLevel) ) THEN
!       endLevel = opt_endLevel
!     ELSE
!       endLevel = n_zlev
!     END IF
!     
!     IF(.NOT.PRESENT(scalar))THEN
!       
!       DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!         CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
!         out_vn_e(:,:,blockNo)= 0.0_wp
!         
!         level_loop_e: DO level = startLevel, endLevel
!           edge_idx_loop: DO je =  start_edge_index, end_edge_index
!             IF (patch_3d%lsm_e(je,level,blockNo) <= sea_boundary) THEN
!               ictr            = 0
!               !IF(patch_3D%lsm_e(je,level,blockNo) == sea)THEN
!               il_c = patch_2d%edges%cell_idx(je,blockNo,1)
!               ib_c = patch_2d%edges%cell_blk(je,blockNo,1)
!               thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
!               
!               DO ie=1, no_primal_edges
!                 ictr =ictr+1
!                 il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!                 ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!                 
!                 thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
!                 thick_frac = thick_edge/thick_cell
!                 out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
!                   & +vn_e(il_e,level,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr))*thick_frac
!                 
!               END DO
!               
!               ictr = no_primal_edges
!               il_c = patch_2d%edges%cell_idx(je,blockNo,2)
!               ib_c = patch_2d%edges%cell_blk(je,blockNo,2)
!               thick_cell = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
!               
!               DO ie=1, no_primal_edges
!                 ictr =ictr+1
!                 il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!                 ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!                 
!                 thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
!                 thick_frac = thick_edge/thick_cell
!                 out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
!                   & +vn_e(il_e,level,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr))*thick_frac
!                 
!               END DO
!             ENDIF
!           END DO edge_idx_loop
!         END DO level_loop_e
!       END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!       
!     ELSEIF(PRESENT(scalar))THEN
!       
!       DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!         CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
!         
!         level_loop_e2: DO level = startLevel, endLevel
!           edge_idx_loop2: DO je =  start_edge_index, end_edge_index
!             IF (patch_3d%lsm_e(je,level,blockNo) <= sea_boundary) THEN
!               ictr            = 0
!               out_vn_e(je,level,blockNo)= 0.0_wp
!               il_c        = patch_2d%edges%cell_idx(je,blockNo,1)
!               ib_c        = patch_2d%edges%cell_blk(je,blockNo,1)
!               scalar_cell = scalar(il_c,level,ib_c)
!               thick_cell  = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
!               
!               DO ie=1, no_primal_edges
!                 ictr =ictr+1
!                 il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!                 ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!                 thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
!                 thick_frac = thick_edge/thick_cell
!                 
!                 out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
!                   & +vn_e(il_e,level,ib_e)*scalar_cell   &
!                   & *(operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr))*thick_frac
!               END DO
!               
!               
!               ictr        = no_primal_edges
!               il_c        = patch_2d%edges%cell_idx(je,blockNo,2)
!               ib_c        = patch_2d%edges%cell_blk(je,blockNo,2)
!               scalar_cell = scalar(il_c,level,ib_c)
!               thick_cell  = patch_3d%p_patch_1d(1)%prism_thick_c(il_c,level,ib_c)
!               
!               DO ie=1, no_primal_edges
!                 ictr =ictr+1
!                 il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!                 ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!                 thick_edge = patch_3d%p_patch_1d(1)%prism_thick_e(il_e,level,ib_e)
!                 thick_frac = thick_edge/thick_cell
!                 
!                 out_vn_e(je,level,blockNo) = out_vn_e(je,level,blockNo) &
!                   & +vn_e(il_e,level,ib_e)*scalar_cell   &
!                   & *(operators_coefficients%edge2edge_viacell_coeff(je,level,blockNo,ictr))*thick_frac
!                 
!               END DO
!               
!             ENDIF
!           END DO edge_idx_loop2
!         END DO level_loop_e2
!       END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!     ENDIF
!     
!   END SUBROUTINE map_edges2edges_viacell_3d_mlev
!   !-----------------------------------------------------------------------------
!   
!   !-----------------------------------------------------------------------------
!   SUBROUTINE map_edges2edges_viacell_3d_1lev( patch_3d, vn_e, operators_coefficients, out_vn_e,scalar,scalar_e, level, &
!     & subset_range)
!     
!     TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
!     REAL(wp), INTENT(in)                       :: vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
!     TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
!     REAL(wp), INTENT(inout)                    :: out_vn_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
!     REAL(wp), INTENT(in), OPTIONAL :: scalar(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     REAL(wp), INTENT(in), OPTIONAL :: scalar_e(nproma,patch_3d%p_patch_2d(1)%nblks_e)
!     INTEGER, INTENT(in), OPTIONAL :: level       ! optional vertical start level
!     TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
!     !Local variables
!     INTEGER :: startLevel, endLevel
!     INTEGER :: start_edge_index, end_edge_index
!     INTEGER :: il_e, ib_e, il_c, ib_c, ictr
!     INTEGER :: je, blockNo, ie, lev
!     REAL(wp) :: scalar_cell
!     TYPE(t_subset_range), POINTER :: edges_inDomain
!     TYPE(t_patch), POINTER :: patch_2d
!     !-----------------------------------------------------------------------
!     patch_2d   => patch_3d%p_patch_2d(1)
!     !-----------------------------------------------------------------------
!     IF ( PRESENT(subset_range) ) THEN
!       edges_inDomain => subset_range
!     ELSE
!       edges_inDomain => patch_2d%edges%in_domain
!     ENDIF
!     !--------------------------
!     IF ( PRESENT(level) ) THEN
!       lev = level
!     ELSE
!       lev = 1
!     END IF
!     
!     IF(.NOT.PRESENT(scalar))THEN
!       
!       DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!         CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
!         
!         out_vn_e(:,blockNo) = 0.0_wp
!         edge_idx_loop: DO je =  start_edge_index, end_edge_index
!           IF (patch_3d%lsm_e(je,lev,blockNo) <= sea_boundary) THEN
!             ictr          = 0
!             il_c = patch_2d%edges%cell_idx(je,blockNo,1)
!             ib_c = patch_2d%edges%cell_blk(je,blockNo,1)
!             
!             !thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)
!             DO ie=1, no_primal_edges
!               ictr =ictr+1
!               il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!               ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!               !thick_edge=patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
!               !thick_frac=thick_edge/thick_cell
!               out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
!                 & +vn_e(il_e,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,lev,blockNo,ictr))!*thick_frac
!               
!             END DO
!             
!             ictr = no_primal_edges
!             il_c = patch_2d%edges%cell_idx(je,blockNo,2)
!             ib_c = patch_2d%edges%cell_blk(je,blockNo,2)
!             
!             !thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)
!             DO ie=1, no_primal_edges
!               ictr =ictr+1
!               il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!               ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!               !thick_edge=patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
!               !thick_frac=thick_edge/thick_cell
!               out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
!                 & +vn_e(il_e,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,lev,blockNo,ictr))!*thick_frac
!               
!             END DO
!             
!           ENDIF
!         END DO edge_idx_loop
!       END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!       
!     ELSEIF(PRESENT(scalar))THEN
!       
!       DO blockNo = edges_inDomain%start_block, edges_inDomain%end_block
!         CALL get_index_range(edges_inDomain, blockNo, start_edge_index, end_edge_index)
!         
!         edge_idx_loop2: DO je =  start_edge_index, end_edge_index
!           IF (patch_3d%lsm_e(je,lev,blockNo) <= sea_boundary) THEN
!             ictr = 0
!             out_vn_e(je,blockNo) = 0.0_wp
!             il_c        = patch_2d%edges%cell_idx(je,blockNo,1)
!             ib_c        = patch_2d%edges%cell_blk(je,blockNo,1)
!             scalar_cell = scalar(il_c,ib_c)
!             
!             DO ie=1, no_primal_edges
!               ictr =ictr+1
!               il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!               ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!               
!               out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
!                 & +vn_e(il_e,ib_e)*scalar_cell &
!                 & *(operators_coefficients%edge2edge_viacell_coeff(je,lev,blockNo,ictr))!*thick_frac
!               
!             END DO
!             ictr        = no_primal_edges
!             il_c        = patch_2d%edges%cell_idx(je,blockNo,2)
!             ib_c        = patch_2d%edges%cell_blk(je,blockNo,2)
!             scalar_cell = scalar(il_c,ib_c)
!             !thick_cell  = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)
!             
!             DO ie=1, no_primal_edges
!               ictr =ictr+1
!               il_e = patch_2d%cells%edge_idx(il_c,ib_c,ie)
!               ib_e = patch_2d%cells%edge_blk(il_c,ib_c,ie)
!               
!               !thick_edge=patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
!               !thick_frac=thick_edge/thick_cell
!               out_vn_e(je,blockNo) = out_vn_e(je,blockNo) &
!                 & +vn_e(il_e,ib_e)*scalar_cell&
!                 & *(operators_coefficients%edge2edge_viacell_coeff(je,lev,blockNo,ictr))!*thick_frac
!               
!             END DO
!           ENDIF
!         END DO edge_idx_loop2
!       END DO
!     ENDIF
!   END SUBROUTINE map_edges2edges_viacell_3d_1lev
!   !-----------------------------------------------------------------------------
!   !-----------------------------------------------------------------------------
!   !>
!   !!
!   !! @par Revision History
!   !!  developed by Peter Korn, MPI-M (2010-11)
!   SUBROUTINE map_edges2cell_with_height_3d( patch_3d, vn_e, operators_coefficients, p_vn_c, h_e,&
!     & opt_startLevel, opt_endLevel, subset_range)
!     
!     TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
!     REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
!     ! 3D case: h_e is surface endLevelation at edges
!     TYPE(t_cartesian_coordinates),INTENT(inout):: p_vn_c(:,:,:)  ! outputput (nproma,n_zlev,alloc_cell_blocks)
!     REAL(wp), INTENT(in)                       :: h_e(:,:)       ! SW-case: h_e is thickness at edges
!     TYPE(t_operator_coeff)                     :: operators_coefficients
!     INTEGER, INTENT(in), OPTIONAL :: opt_startLevel       ! optional vertical start level
!     INTEGER, INTENT(in), OPTIONAL :: opt_endLevel       ! optional vertical end level
!     TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
!     
!     !Local variables
!     !INTEGER, PARAMETER :: no_primal_edges = 3
!     INTEGER :: startLevel, endLevel
!     INTEGER :: start_cell_index, end_cell_index
!     INTEGER :: il_e, ib_e
!     INTEGER :: cell_index, blockNo, level, ie!,je
!     REAL(wp) :: z_weight
!     REAL(wp) :: z_thick_e
!     
!     TYPE(t_subset_range), POINTER :: all_cells
!     !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
!     !  & routine = ('mo_scalar_product:primal_map_e2c')
!     TYPE(t_patch), POINTER :: patch_2d
!     !-----------------------------------------------------------------------
!     patch_2d   => patch_3d%p_patch_2d(1)
!     !-----------------------------------------------------------------------
!     IF ( PRESENT(subset_range) ) THEN
!       all_cells => subset_range
!     ELSE
!       all_cells => patch_2d%cells%ALL
!     ENDIF
!     !-----------------------------------------------------------------------
!     IF ( PRESENT(opt_startLevel) ) THEN
!       startLevel = opt_startLevel
!     ELSE
!       startLevel = 1
!     END IF
!     IF ( PRESENT(opt_endLevel) ) THEN
!       endLevel = opt_endLevel
!     ELSE
!       endLevel = n_zlev
!     END IF
!     
!     IF ( iswm_oce == 1 ) THEN
!       
!       !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
!       DO blockNo = all_cells%start_block, all_cells%end_block
!         CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
!         level_loop_swm: DO level = startLevel, endLevel
!           cell_idx_loop_swm: DO cell_index =  start_cell_index, end_cell_index
!             !calculate velocity reconstruction at cell center
!             z_weight           = 0.0_wp
!             p_vn_c(cell_index,level,blockNo)%x = 0.0_wp
!             DO ie=1, no_primal_edges
!               
!               il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
!               ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)
!               
!               z_thick_e =patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(il_e,level,ib_e)&
!                 & + h_e(il_e,ib_e)
!               z_weight = z_weight + operators_coefficients%variable_vol_norm(cell_index,level,blockNo,ie) * z_thick_e
!               
!               p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x&
!                 & + operators_coefficients%edge2cell_coeff_cc_dyn(cell_index,level,blockNo,ie)%x&
!                 & * vn_e(il_e,level,ib_e)* z_thick_e
!             END DO
!             IF( z_weight/=0.0_wp)THEN
!               p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x / z_weight
!             ELSE
!               p_vn_c(cell_index,level,blockNo)%x=0.0_wp
!             ENDIF
!           END DO cell_idx_loop_swm
!         END DO level_loop_swm
!       END DO ! blockNo = all_cells%start_block, all_cells%end_block
!       
!     ELSEIF( iswm_oce /= 1 ) THEN
!       
!       !Step 1: Calculation of Pv in cartesian coordinates
!       DO blockNo = all_cells%start_block, all_cells%end_block
!         CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
!         
!         !We are dealing with the surface layer first
!         cell_idx_loop_top: DO cell_index =  start_cell_index, end_cell_index
!           z_weight             = 0.0_wp
!           p_vn_c(cell_index,startLevel,blockNo)%x = 0.0_wp
!           DO ie=1, no_primal_edges
!             
!             il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
!             ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)
!             
!             z_thick_e = patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(il_e,startLevel,ib_e)&
!               & + h_e(il_e,ib_e)
!             z_weight = z_weight + operators_coefficients%variable_vol_norm(cell_index,startLevel,blockNo,ie) * z_thick_e
!             
!             p_vn_c(cell_index,startLevel,blockNo)%x = p_vn_c(cell_index,startLevel,blockNo)%x&
!               & + operators_coefficients%edge2cell_coeff_cc_dyn(cell_index,1,blockNo,ie)%x&
!               & * vn_e(il_e,startLevel,ib_e) * z_thick_e
!             
!           END DO
!           
!           IF(z_weight/=0.0_wp)THEN
!             p_vn_c(cell_index,startLevel,blockNo)%x = p_vn_c(cell_index,startLevel,blockNo)%x / z_weight
!           ELSE
!             p_vn_c(cell_index,startLevel,blockNo)%x=0.0_wp
!           ENDIF
!         END DO cell_idx_loop_top
!         
!         !Now we calculate at the levels below the surface
!         level_loop: DO level = startLevel+1, endLevel
!           cell_idx_loop: DO cell_index =  start_cell_index, end_cell_index
!             p_vn_c(cell_index,level,blockNo)%x = 0.0_wp
!             !z_weight = 0.0_wp
!             DO ie=1, no_primal_edges
!               
!               il_e = patch_2d%cells%edge_idx(cell_index,blockNo,ie)
!               ib_e = patch_2d%cells%edge_blk(cell_index,blockNo,ie)
!               p_vn_c(cell_index,level,blockNo)%x = p_vn_c(cell_index,level,blockNo)%x&
!                 & + operators_coefficients%edge2cell_coeff_cc(cell_index,level,blockNo,ie)%x&
!                 & * vn_e(il_e,level,ib_e)
!             END DO
!           END DO cell_idx_loop
!         END DO level_loop
!       END DO ! blockNo = all_cells%start_block, all_cells%end_block
!     ENDIF
!     
!     ! LL no sync required
!     
!   END SUBROUTINE map_edges2cell_with_height_3d
!   !----------------------------------------------------------------
  
  ! !   !-----------------------------------------------------------------------------
  ! !   !>
  ! !   !! Discrete mapping of cell-based vectors to edges on the primal grid.
  ! !   !!
  ! !   !!
  ! !   !! @par Revision History
  ! !   !!  developed by Peter Korn, MPI-M (2010-11)
  ! !   !!  mpi parallelized LL, result not synced
  ! !   SUBROUTINE map_cell2edges_2d( patch_3D, p_vn_c, ptp_vn, operators_coefficients)
  ! !
  ! !     TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3D
  ! !     TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:)    ! input vector (nproma,n_zlev,alloc_cell_blocks)
  ! !     REAL(wp), INTENT(inout)                    :: ptp_vn(:,:)    ! output vector (nproma,n_zlev,nblks_e)
  ! !     TYPE(t_operator_coeff)                     :: operators_coefficients
  ! !
  ! !     !Local variables
  ! !     INTEGER :: start_edge_index, end_edge_index
  ! !     INTEGER :: cell_1_index, cell_1_block, cell_2_index, cell_2_block
  ! !     INTEGER :: je, blockNo
  ! !     TYPE(t_subset_range), POINTER :: edges_in_domain
  ! !     TYPE(t_patch), POINTER        :: patch_2D
  ! !     !-----------------------------------------------------------------------
  ! !     patch_2D   => patch_3D%p_patch_2D(1)
  ! !     !-----------------------------------------------------------------------
  ! !     !CALL message (TRIM(routine), 'start')
  ! !     edges_in_domain   => patch_2D%edges%in_domain
  ! !
  ! !     ! calculation of transposed P^TPv from Pv (incart coord)
  ! !     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
  ! !       CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
  ! !
  ! !       edge_idx_loop: DO je =  start_edge_index, end_edge_index
  ! !
  ! !         IF(patch_3D%lsm_e(je,1,blockNo) <= sea_boundary)THEN
  ! !           !Get indices of two adjacent triangles
  ! !           cell_1_index = patch_2D%edges%cell_idx(je,blockNo,1)
  ! !           cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
  ! !           cell_2_index = patch_2D%edges%cell_idx(je,blockNo,2)
  ! !           cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)
  ! !
  ! !           ptp_vn(je,blockNo) = &
  ! !             &  DOT_PRODUCT(p_vn_c(cell_1_index,cell_1_block)%x, operators_coefficients%edge2cell_coeff_cc_t(je,1,blockNo,1)%x)&
  ! !             & +DOT_PRODUCT(p_vn_c(cell_2_index,cell_2_block)%x, operators_coefficients%edge2cell_coeff_cc_t(je,1,blockNo,2)%x)
  ! !         ELSE
  ! !           ptp_vn(je,blockNo) = 0.0_wp
  ! !         ENDIF
  ! !
  ! !       END DO edge_idx_loop
  ! !     END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
  ! !
  ! ! !     CALL sync_patch_array(SYNC_E, patch_2D, ptp_vn)
  ! !
  ! !   END SUBROUTINE map_cell2edges_2d
  !-----------------------------------------------------------------------------
  
  
  
END MODULE mo_scalar_product

