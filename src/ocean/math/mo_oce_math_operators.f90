!>
!!   Contains the implementation of the mathematical operators for the ocean.
!!
!!   Contains the implementation of the mathematical operators for the ocean.
!!
!! @par Revision History
!!  Developed  by Peter Korn and Stephan Lorenz 2010-04
!!  Modified by Stephan Lorenz                  2011-02
!!    correct implementation of ocean boundaries
!!
!! @par To Do
!! Boundary exchange, nblks in presence of halos and dummy edge
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_oce_math_operators
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp, sp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: finish,message
  USE mo_run_config,         ONLY: ltimer, dtime
  USE mo_math_constants
  USE mo_physical_constants
  USE mo_impl_constants,     ONLY: boundary, sea_boundary, min_dolic !,sea,land, land_boundary, sea, max_char_length, &
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce, &
    & select_solver, select_restart_mixedprecision_gmres
  
  USE mo_dynamics_config,    ONLY: nold
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_div, timer_grad
  USE mo_oce_types,          ONLY: t_hydro_ocean_state, t_solvercoeff_singleprecision, &
    & t_verticaladvection_ppm_coefficients, t_operator_coeff
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, vector_product !, gc2cc
!   USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array
  USE mo_grid_config,         ONLY: n_dom
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=12) :: str_module    = 'oceMathOps  '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  
  PUBLIC :: grad_fd_norm_oce_3D
  PUBLIC :: grad_fd_norm_oce_3D_onblock
  PUBLIC :: div_oce_3D, div_oce_2D_sp
  PUBLIC :: div_oce_2D_onTriangles_onBlock, div_oce_2D_onTriangles_onBlock_sp, div_oce_3D_onTriangles_onBlock
  PUBLIC :: rot_vertex_ocean_3D
  PUBLIC :: grad_fd_norm_oce_2D_3D, grad_fd_norm_oce_2D_3D_sp
  PUBLIC :: grad_fd_norm_oce_2D_onBlock
  PUBLIC :: verticalDeriv_vec_midlevel_on_block
  PUBLIC :: verticalDeriv_scalar_midlevel_on_block
  PUBLIC :: verticalDiv_scalar_midlevel
  PUBLIC :: verticalDiv_scalar_midlevel_on_block
  PUBLIC :: calculate_thickness
  PUBLIC :: map_edges2vert_3D
  PUBLIC :: check_cfl_horizontal, check_cfl_vertical
  PUBLIC :: smooth_onCells
  
  
  INTERFACE div_oce_3D
    MODULE PROCEDURE div_oce_3D_mlevels
    MODULE PROCEDURE div_oce_3D_1level
  END INTERFACE
  
  INTERFACE smooth_onCells
    MODULE PROCEDURE   smooth_onCells_3D
    MODULE PROCEDURE   smooth_onCells_2D
  END INTERFACE
  
CONTAINS
  
  
  !-------------------------------------------------------------------------
  !>
!<Optimize:inUse>
  SUBROUTINE map_edges2vert_3D(patch_2D, vn, edge2vert_coeff_cc, vn_dual)
    
    TYPE(t_patch), TARGET, INTENT(in)       :: patch_2D
    REAL(wp), INTENT(in)                    :: vn(:,:,:)
    TYPE(t_cartesian_coordinates),INTENT(in):: edge2vert_coeff_cc(:,:,:,:)
    TYPE(t_cartesian_coordinates)           :: vn_dual(nproma,n_zlev,patch_2D%nblks_v)
    
    INTEGER :: start_level, end_level
    INTEGER :: vertexIndex, level, blockNo,vertexConnect
    INTEGER :: edgeOfVertex_index, edgeOfVertex_block
    INTEGER :: start_index_v, end_index_v
    TYPE(t_subset_range), POINTER :: verts_in_domain
    !-----------------------------------------------------------------------
    verts_in_domain => patch_2D%verts%in_domain
    
    !i_v_ctr(:,:,:) = 0
    start_level         = 1
    end_level         = n_zlev
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_index_v,end_index_v, vertexIndex, vertexConnect, &
!ICON_OMP edgeOfVertex_index, edgeOfVertex_block, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, blockNo, start_index_v, end_index_v)
      vn_dual(:,:,blockNo)%x(1) = 0.0_wp
      vn_dual(:,:,blockNo)%x(2) = 0.0_wp
      vn_dual(:,:,blockNo)%x(3) = 0.0_wp
      DO vertexIndex = start_index_v, end_index_v
          
        DO vertexConnect = 1, patch_2D%verts%num_edges(vertexIndex,blockNo)

          edgeOfVertex_index = patch_2D%verts%edge_idx(vertexIndex,blockNo,vertexConnect)
          edgeOfVertex_block = patch_2D%verts%edge_blk(vertexIndex,blockNo,vertexConnect)
            
          DO level = start_level, end_level
            vn_dual(vertexIndex,level,blockNo)%x = vn_dual(vertexIndex,level,blockNo)%x        &
              & + edge2vert_coeff_cc(vertexIndex,level,blockNo,vertexConnect)%x &
              & * vn(edgeOfVertex_index,level,edgeOfVertex_block)
            
          END DO
        END DO ! vertexIndex = start_index_v, end_index_v
      END DO ! level = start_level, end_level
    END DO ! blockNo = verts_in_domain%start_block, verts_in_domain%end_block
!ICON_OMP_END_PARALLEL_DO    
    
  END SUBROUTINE map_edges2vert_3D
  !-------------------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE grad_fd_norm_oce_3D( psi_c, patch_3D, grad_coeff, grad_norm_psi_e)
    
    TYPE(t_patch_3D ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)                   :: grad_coeff(:,:,:)!(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(in)                   :: psi_c          (nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: grad_norm_psi_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    
    !
    INTEGER :: start_edge_index, end_edge_index, blockNo
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => patch_3D%p_patch_2D(1)%edges%in_domain
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      CALL grad_fd_norm_oce_3D_onblock( psi_c, patch_3D, &
        & grad_coeff(:,:,blockNo),      &
        & grad_norm_psi_e(:,:,blockNo), &
        & start_edge_index, end_edge_index, blockNo)
    END DO
!ICON_OMP_END_PARALLEL_DO
    
  END SUBROUTINE grad_fd_norm_oce_3D
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!  Computes directional derivative of a cell centered variable in presence of lateral boundaries
  !!  as in the ocean model setting.
  !!
  !!  Computes directional  derivative of a cell centered variable
  !!  with respect to direction normal to triangle edge.
  !! input: lives on centres of triangles
  !! output:  lives on edges (velocity points)
  !!
  !! @par Revision History
  !! Developed  by  Luca Bonaventura, MPI-M (2002-5).
  !! Adapted to new data structure by Peter Korn
  !! and Luca Bonaventura, MPI-M (2005).
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! -Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch_2D
  !!Boundary handling for triangles by P. Korn (2009)
  !!
  !!  mpi note: the result is on edges_in_domain.
!<Optimize:inUse>
  SUBROUTINE grad_fd_norm_oce_3D_onblock( psi_c, patch_3D, grad_coeff, grad_norm_psi_e, &
    & start_edge_index, end_edge_index, blockNo)
    
    TYPE(t_patch_3D ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)                   :: grad_coeff(:,:)!(nproma,n_zlev)
    REAL(wp), INTENT(in)                   :: psi_c          (nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: grad_norm_psi_e(nproma,n_zlev)
    INTEGER, INTENT(in)                    :: start_edge_index, end_edge_index, blockNo
    
    INTEGER :: je, level
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
    !-----------------------------------------------------------------------
    
    iidx => patch_3D%p_patch_2D(1)%edges%cell_idx
    iblk => patch_3D%p_patch_2D(1)%edges%cell_blk
!     grad_norm_psi_e(:,:) = 0.0_wp
    
    DO je = start_edge_index, end_edge_index
      DO level = 1, patch_3D%p_patch_1d(1)%dolic_e(je,blockNo)
        grad_norm_psi_e(je,level) =                                        &
          & grad_coeff(je,level) *                                         &
          & ( psi_c(iidx(je,blockNo,2),level,iblk(je,blockNo,2)) -       &
          & psi_c(iidx(je,blockNo,1),level,iblk(je,blockNo,1)) )
      ENDDO
    END DO
    
  END SUBROUTINE grad_fd_norm_oce_3D_onblock
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Computes discrete divergence of a vector field in presence of lateral boundaries as in ocean setting.
  !!
  !! Computes discrete divergence of a vector field
  !! given by its components in the directions normal to triangle edges.
  !! The midpoint rule is used for quadrature.
  !! input:  lives on edges (velocity points)
  !! output: lives on centers of triangles
  !!
  !! @par Revision History
  !! Developed  by  Luca Bonaventura, MPI-M (2002-5).
  !! Changes according to programming guide by Thomas Heinze, DWD (2006-08-18).
  !! Modification by Thomas Heinze, DWD (2006-09-11):
  !! - loop only over the inner cells of a patch_2D, not any more over halo cells
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! - Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch_2D
  !! Modification by Guenther Zaengl, DWD (2009-03-17)
  !! - vector optimization
  !! Modification by Peter Korn, MPI-M    (2009)
  !! - Boundary treatment for the ocean
  !! Modification by Stephan Lorenz, MPI-M (2010-08-05)
  !! - New boundary definition with inner and boundary points on land/sea
  !!
!<Optimize:inUse>
  SUBROUTINE div_oce_3D_mlevels_onTriangles( vec_e, patch_3D, div_coeff, div_vec_c, opt_start_level, opt_end_level, &
    & subset_range)
    
    TYPE(t_patch_3D ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)          :: vec_e(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)          :: div_coeff(:,:,:,:)
    REAL(wp), INTENT(inout)       :: div_vec_c(:,:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level       ! optional vertical end level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range

    INTEGER :: start_level, end_level
    INTEGER :: blockNo
    INTEGER ::start_index, end_index
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: cells_subset
    !-----------------------------------------------------------------------
    IF (PRESENT(subset_range)) THEN
      cells_subset => subset_range
    ELSE
      cells_subset => patch_3D%p_patch_2D(1)%cells%in_domain
    ENDIF
    
    IF ( PRESENT(opt_start_level) ) THEN
      start_level = opt_start_level
    ELSE
      start_level = 1
    END IF
    IF ( PRESENT(opt_end_level) ) THEN
      end_level = opt_end_level
    ELSE
      end_level = n_zlev
    END IF
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_subset%start_block, cells_subset%end_block
      CALL get_index_range(cells_subset, blockNo, start_index, end_index)

      CALL div_oce_3D_onTriangles_onBlock( vec_e, patch_3D, div_coeff, div_vec_c(:,:,blockNo), &
        & blockNo, start_index, end_index, start_level, end_level)
    END DO
!ICON_OMP_END_PARALLEL_DO
    
  END SUBROUTINE div_oce_3D_mlevels_onTriangles
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  ! compute the discrete divergence for cell jc by finite volume
  ! approximation (see Bonaventura and Ringler MWR 2005);
  ! multiplication of the normal vector component vec_e at the edges
  ! by the appropriate cell based edge_orientation is required to
  ! obtain the correct value for the application of Gauss theorem
  ! (which requires the scalar product of the vector field with the
  ! OUTWARD pointing unit vector with respect to cell jc; since the
  ! positive direction for the vector components is not necessarily
  ! the outward pointing one with respect to cell jc, a correction
  ! coefficient (equal to +-1) is necessary, given by
  ! patch_2D%grid%cells%edge_orientation)
  !
  ! Distinghuish: case of a land cell (put div to zero), and
  ! cases where one of the edges are boundary or land
  ! (put corresponding velocity to zero).
  ! sea, sea_boundary, boundary (edges only), land_boundary, land =
  !  -2,      -1,         0,                  1,             2
  !This information is stored inside the divergence coefficients.
!<Optimize:inUse>
  SUBROUTINE div_oce_3D_onTriangles_onBlock( vec_e, patch_3D, div_coeff, div_vec_c, &
    & blockNo, start_index, end_index, start_level, end_level)

    TYPE(t_patch_3D ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)          :: vec_e(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)          :: div_coeff(:,:,:,:)
    REAL(wp), INTENT(inout)       :: div_vec_c(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
    INTEGER, INTENT(in)           :: blockNo, start_index, end_index
    INTEGER, INTENT(in) :: start_level, end_level     ! vertical start and end level

    INTEGER :: jc, level
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: cells_subset
    !-----------------------------------------------------------------------

    iidx => patch_3D%p_patch_2D(1)%cells%edge_idx
    iblk => patch_3D%p_patch_2D(1)%cells%edge_blk

    div_vec_c(:,:) = 0.0_wp
    DO jc = start_index, end_index
      DO level = start_level, MIN(end_level, patch_3D%p_patch_1d(1)%dolic_c(jc, blockNo))
        div_vec_c(jc,level) =  &
          & vec_e(iidx(jc,blockNo,1),level,iblk(jc,blockNo,1)) * div_coeff(jc,level,blockNo,1) + &
          & vec_e(iidx(jc,blockNo,2),level,iblk(jc,blockNo,2)) * div_coeff(jc,level,blockNo,2) + &
          & vec_e(iidx(jc,blockNo,3),level,iblk(jc,blockNo,3)) * div_coeff(jc,level,blockNo,3)
      END DO
    END DO

  END SUBROUTINE div_oce_3D_onTriangles_onBlock
  !-------------------------------------------------------------------------

 
  !-------------------------------------------------------------------------
  !>
  !! Computes discrete divergence of a vector field in presence of lateral boundaries as in ocean setting.
  !!
  !! Computes discrete divergence of a vector field
  !! given by its components in the directions normal to triangle edges.
  !! The midpoint rule is used for quadrature.
  !! input:  lives on edges (velocity points)
  !! output: lives on centers of triangles
  !!
  !! @par Revision History
  !! Developed  by  Luca Bonaventura, MPI-M (2002-5).
  !! Changes according to programming guide by Thomas Heinze, DWD (2006-08-18).
  !! Modification by Thomas Heinze, DWD (2006-09-11):
  !! - loop only over the inner cells of a patch_2D, not any more over halo cells
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! - Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch_2D
  !! Modification by Guenther Zaengl, DWD (2009-03-17)
  !! - vector optimization
  !! Modification by Peter Korn, MPI-M    (2009)
  !! - Boundary treatment for the ocean
  !! Modification by Stephan Lorenz, MPI-M (2010-08-05)
  !! - New boundary definition with inner and boundary points on land/sea
  !!
!<Optimize:inUse>
  SUBROUTINE div_oce_3D_mlevels( vec_e, patch_3D, div_coeff, div_vec_c, opt_start_level, opt_end_level, &
    & subset_range)
    
    TYPE(t_patch_3D ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)          :: vec_e(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)          :: div_coeff(:,:,:,:)
    REAL(wp), INTENT(inout)       :: div_vec_c(:,:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level       ! optional vertical end level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    
    INTEGER :: start_level, end_level
    INTEGER :: jc, level, blockNo, max_connectivity, edgeofcell
    INTEGER ::start_index, end_index
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: cells_subset
    !-----------------------------------------------------------------------
    IF ( patch_3D%p_patch_2D(1)%cells%max_connectivity == 3) THEN
      CALL div_oce_3D_mlevels_onTriangles(vec_e, patch_3D, div_coeff, div_vec_c, &
        & opt_start_level, opt_end_level, subset_range)
      RETURN
    !-----------------------------------------------------------------------
    ENDIF
    
    IF (PRESENT(subset_range)) THEN
      cells_subset => subset_range
    ELSE
      cells_subset => patch_3D%p_patch_2D(1)%cells%in_domain
    ENDIF
    
    IF ( PRESENT(opt_start_level) ) THEN
      start_level = opt_start_level
    ELSE
      start_level = 1
    END IF
    IF ( PRESENT(opt_end_level) ) THEN
      end_level = opt_end_level
    ELSE
      end_level = n_zlev
    END IF
    
!ICON_OMP_PARALLEL PRIVATE(iidx, iblk, max_connectivity)
    iidx => patch_3D%p_patch_2D(1)%cells%edge_idx
    iblk => patch_3D%p_patch_2D(1)%cells%edge_blk
    max_connectivity = patch_3D%p_patch_2D(1)%cells%max_connectivity
    
!ICON_OMP_DO PRIVATE(start_index,end_index, jc, level, edgeOfCell) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_subset%start_block, cells_subset%end_block
      CALL get_index_range(cells_subset, blockNo, start_index, end_index)
      div_vec_c(:,:,blockNo) = 0.0_wp
      DO jc = start_index, end_index
        DO level = start_level, MIN(end_level, patch_3D%p_patch_1d(1)%dolic_c(jc, blockNo))
          ! compute the discrete divergence for cell jc by finite volume
          ! approximation (see Bonaventura and Ringler MWR 2005);
          ! multiplication of the normal vector component vec_e at the edges
          ! by the appropriate cell based edge_orientation is required to
          ! obtain the correct value for the application of Gauss theorem
          ! (which requires the scalar product of the vector field with the
          ! OUTWARD pointing unit vector with respect to cell jc; since the
          ! positive direction for the vector components is not necessarily
          ! the outward pointing one with respect to cell jc, a correction
          ! coefficient (equal to +-1) is necessary, given by
          ! patch_2D%grid%cells%edge_orientation)
          !
          ! Distinghuish: case of a land cell (put div to zero), and
          ! cases where one of the edges are boundary or land
          ! (put corresponding velocity to zero).
          ! sea, sea_boundary, boundary (edges only), land_boundary, land =
          !  -2,      -1,         0,                  1,             2
          !This information is stored inside the divergence coefficients.
          div_vec_c(jc,level,blockNo) = 0.0_wp
          
          DO edgeofcell = 1, max_connectivity
            div_vec_c(jc,level,blockNo) = div_vec_c(jc,level,blockNo) + &
              & vec_e(iidx(jc,blockNo,edgeofcell),level,iblk(jc,blockNo,edgeofcell)) * &
              & div_coeff(jc,level,blockNo,edgeofcell)
          END DO
          
        END DO
      END DO
    END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE div_oce_3D_mlevels
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Computes discrete divergence of a vector field in presence of lateral boundaries as in ocean setting.
  !!
  !! Computes discrete divergence of a vector field
  !! given by its components in the directions normal to triangle edges.
  !! The midpoint rule is used for quadrature.
  !! input:  lives on edges (velocity points)
  !! output: lives on centers of triangles
  !!
  !! @par Revision History
  !! Developed  by  Luca Bonaventura, MPI-M (2002-5).
  !! Changes according to programming guide by Thomas Heinze, DWD (2006-08-18).
  !! Modification by Thomas Heinze, DWD (2006-09-11):
  !! - loop only over the inner cells of a patch_2D, not any more over halo cells
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! - Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch_2D
  !! Modification by Guenther Zaengl, DWD (2009-03-17)
  !! - vector optimization
  !! Modification by Peter Korn, MPI-M    (2009)
  !! - Boundary treatment for the ocean
  !! Modification by Stephan Lorenz, MPI-M (2010-08-05)
  !! - New boundary definition with inner and boundary points on land/sea
  !!
  SUBROUTINE div_oce_3D_1level( vec_e, patch_2D, div_coeff, div_vec_c,  &
    & level, subset_range)
    
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    !
    ! edge based variable of which divergence
    ! is computed
    !
    REAL(wp), INTENT(inout)       :: vec_e(:,:) ! dim: (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)          :: div_coeff(:,:,:,:)
    REAL(wp), INTENT(inout)       :: div_vec_c(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
    INTEGER,  INTENT(in)          :: level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    
    INTEGER :: jc, blockNo
    INTEGER :: start_index, end_index
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    IF (PRESENT(subset_range)) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2D%cells%ALL
    ENDIF
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      CALL div_oce_2D_onTriangles_onBlock( vec_e, patch_2D, div_coeff, div_vec_c(:,blockNo),  &
         & level, blockNo, start_index, end_index)
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE div_oce_3D_1level
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  ! compute the discrete divergence for cell jc by finite volume
  ! approximation (see Bonaventura and Ringler MWR 2005);
  ! multiplication of the normal vector component vec_e at the edges
  ! by the appropriate cell based edge_orientation is required to
  ! obtain the correct value for the application of Gauss theorem
  ! (which requires the scalar product of the vector field with the
  ! OUTWARD pointing unit vector with respect to cell jc; since the
  ! positive direction for the vector components is not necessarily
  ! the outward pointing one with respect to cell jc, a correction
  ! coefficient (equal to +-1) is necessary, given by
  ! patch_2D%grid%cells%edge_orientation)
  !
  ! Distinghuish: case of a land cell (put div to zero), and
  ! cases where one of the edges are boundary or land
  ! (put corresponding velocity to zero).
  ! sea, sea_boundary, boundary (edges only), land_boundary, land =
  !  -2,      -1,         0,                  1,             2
  !This information is stored inside the divergence coefficients.
!<Optimize:inUse>
  SUBROUTINE div_oce_2D_onTriangles_onBlock( vec_e, patch_2D, div_coeff, div_vec_c,  &
    & level, blockNo, start_index, end_index)

    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    !
    ! edge based variable of which divergence
    ! is computed
    !
    REAL(wp), INTENT(inout)       :: vec_e(:,:) ! dim: (nproma,nblks_e)
    REAL(wp), INTENT(in)          :: div_coeff(:,:,:,:)
    REAL(wp), INTENT(inout)       :: div_vec_c(:) ! dim: (nproma)
    INTEGER,  INTENT(in)          :: level
    INTEGER,  INTENT(in) :: blockNo, start_index, end_index

    INTEGER :: jc
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    !-----------------------------------------------------------------------

    iidx => patch_2D%cells%edge_idx
    iblk => patch_2D%cells%edge_blk

    DO jc = start_index, end_index

      div_vec_c(jc) =  &
        & vec_e(iidx(jc,blockNo,1),iblk(jc,blockNo,1)) * div_coeff(jc,level,blockNo,1) + &
        & vec_e(iidx(jc,blockNo,2),iblk(jc,blockNo,2)) * div_coeff(jc,level,blockNo,2) + &
        & vec_e(iidx(jc,blockNo,3),iblk(jc,blockNo,3)) * div_coeff(jc,level,blockNo,3)
    END DO
    
  END SUBROUTINE div_oce_2D_onTriangles_onBlock
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE div_oce_2D_onTriangles_onBlock_sp( vec_e, patch_2D, div_coeff, div_vec_c,  &
    &  blockNo, start_index, end_index)

    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    !
    ! edge based variable of which divergence
    ! is computed
    !
    REAL(sp), INTENT(inout)       :: vec_e(:,:) ! dim: (nproma,nblks_e)
    REAL(sp), INTENT(in)          :: div_coeff(:,:,:)
    REAL(sp), INTENT(inout)       :: div_vec_c(:) ! dim: (nproma)
    INTEGER,  INTENT(in) :: blockNo, start_index, end_index

    INTEGER :: jc
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    !-----------------------------------------------------------------------

    iidx => patch_2D%cells%edge_idx
    iblk => patch_2D%cells%edge_blk

    DO jc = start_index, end_index

      div_vec_c(jc) =  &
        & vec_e(iidx(jc,blockNo,1),iblk(jc,blockNo,1)) * div_coeff(jc,blockNo,1) + &
        & vec_e(iidx(jc,blockNo,2),iblk(jc,blockNo,2)) * div_coeff(jc,blockNo,2) + &
        & vec_e(iidx(jc,blockNo,3),iblk(jc,blockNo,3)) * div_coeff(jc,blockNo,3)
    END DO

  END SUBROUTINE div_oce_2D_onTriangles_onBlock_sp
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! as div_oce_3D_1level in single precisison and 2D
  SUBROUTINE div_oce_2D_sp( vec_e, patch_2D, div_coeff, div_vec_c,  &
    & subset_range)
    
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    !
    ! edge based variable of which divergence is computed
    REAL(sp), INTENT(inout)       :: vec_e(:,:) ! dim: (nproma,n_zlev,nblks_e)
    REAL(sp), INTENT(in)          :: div_coeff(:,:,:)
    REAL(sp), INTENT(inout)       :: div_vec_c(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    
    INTEGER :: jc, blockNo
    INTEGER :: start_index, end_index
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    IF (PRESENT(subset_range)) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2D%cells%ALL
    ENDIF
    
    IF (ltimer) CALL timer_start(timer_div)
    ! !$OMP PARALLEL
    
    iidx => patch_2D%cells%edge_idx
    iblk => patch_2D%cells%edge_blk
    
    ! !$OMP DO PRIVATE(blockNo,start_index,end_index,jc)
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index, end_index)
      DO jc = start_index, end_index
        
        div_vec_c(jc,blockNo) =  &
          & vec_e(iidx(jc,blockNo,1),iblk(jc,blockNo,1)) * div_coeff(jc,blockNo,1) + &
          & vec_e(iidx(jc,blockNo,2),iblk(jc,blockNo,2)) * div_coeff(jc,blockNo,2) + &
          & vec_e(iidx(jc,blockNo,3),iblk(jc,blockNo,3)) * div_coeff(jc,blockNo,3)
      END DO
    END DO
    ! !$OMP END DO
    
    ! !$OMP END PARALLEL
    IF (ltimer) CALL timer_stop(timer_div)
  END SUBROUTINE div_oce_2D_sp
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE grad_fd_norm_oce_2D_3D( psi_c, patch_2D, grad_coeff, grad_norm_psi_e, subset_range)
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    REAL(wp), INTENT(in)    :: psi_c(:,:)             ! dim: (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(in)    :: grad_coeff(:,:)
    REAL(wp), INTENT(inout) ::  grad_norm_psi_e(:,:)  ! dim: (nproma,nblks_e)
    TYPE(t_subset_range), TARGET, OPTIONAL :: subset_range

    INTEGER :: blockNo
    INTEGER :: start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: edges_in_domain

   !-----------------------------------------------------------------------
    IF (PRESENT(subset_range)) THEN
      edges_in_domain => subset_range
    ELSE
      edges_in_domain => patch_2D%edges%in_domain
    ENDIF    
    !-----------------------------------------------------------------------
    
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_edge_index,end_edge_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      CALL grad_fd_norm_oce_2D_onBlock(psi_c,  patch_2D, grad_coeff(:,blockNo), grad_norm_psi_e(:,blockNo), &
        & start_edge_index, end_edge_index, blockNo)
      
    END DO
!ICON_OMP_END_PARALLEL_DO
    
  END SUBROUTINE grad_fd_norm_oce_2D_3D
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!  Computes directional derivative of a cell centered variable in presence of lateral boundaries
  !!  as in the ocean model setting.
  !!
  !!  Computes directional  derivative of a cell centered variable
  !!  with respect to direction normal to triangle edge.
  !! input: lives on centres of triangles
  !! output:  lives on edges (velocity points)
  !!
  !! @par Revision History
  !! Developed  by  Luca Bonaventura, MPI-M (2002-5).
  !! Adapted to new data structure by Peter Korn
  !! and Luca Bonaventura, MPI-M (2005).
  !! Modifications by P. Korn, MPI-M(2007-2)
  !! -Switch fom array arguments to pointers
  !! Modification by Almut Gassmann, MPI-M (2007-04-20)
  !! - abandon grid for the sake of patch_2D
  !! Boundary handling for triangles by P. Korn (2009)
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  !!
!<Optimize:inUse>
  SUBROUTINE grad_fd_norm_oce_2D_onBlock(psi_c,  patch_2D, grad_coeff, grad_norm_psi_e, start_index, end_index, blockNo)
    !
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    REAL(wp), INTENT(in)    ::  psi_c(:,:)               ! dim: (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(in)    ::  grad_coeff(:)
    REAL(wp), INTENT(inout) ::  grad_norm_psi_e(:)   ! dim: (nproma)
    INTEGER, INTENT(in)     :: start_index, end_index, blockNo
    
    INTEGER :: je
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    !-----------------------------------------------------------------------
    
    iidx => patch_2D%edges%cell_idx
    iblk => patch_2D%edges%cell_blk
    
    DO je = start_index, end_index
      ! compute the normal derivative
      ! by the finite difference approximation
      ! (see Bonaventura and Ringler MWR 2005)
!       IF (iidx(je,blockNo,1) < 1 .or. iidx(je,blockNo,2) < 1) THEN
!         WRITE(0,*) "je=", je, " blockNo=", blockNo, &
!           & " iidx(je,blockNo,1)=", iidx(je,blockNo,1), " iidx(je,blockNo,2)=", iidx(je,blockNo,2)
!         CALL finish("invalid connectivity", "")
!       ENDIF
      
      grad_norm_psi_e(je) =  &
        & (psi_c(iidx(je,blockNo,2),iblk(je,blockNo,2))-psi_c(iidx(je,blockNo,1),iblk(je,blockNo,1)))&
        & * grad_coeff(je)
      
    END DO
    
  END SUBROUTINE grad_fd_norm_oce_2D_onBlock
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  ! the same as grad_fd_norm_oce_2D_3D_sp in single precisison
!<Optimize:inUse>
  SUBROUTINE grad_fd_norm_oce_2D_3D_sp( psi_c, patch_2D, grad_coeff, grad_norm_psi_e, subset_range)
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    REAL(sp), INTENT(in)    :: psi_c(:,:)             ! dim: (nproma,alloc_cell_blocks)
    REAL(sp), INTENT(in)    :: grad_coeff(:,:)
    REAL(sp), INTENT(inout) ::  grad_norm_psi_e(:,:)  ! dim: (nproma,nblks_e)
    TYPE(t_subset_range), TARGET, OPTIONAL :: subset_range
    
    INTEGER :: je, blockNo
    INTEGER :: start_edge_index, end_edge_index
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    IF (PRESENT(subset_range)) THEN
      edges_in_domain => subset_range
    ELSE
      edges_in_domain => patch_2D%edges%in_domain
    ENDIF
    iidx => patch_2D%edges%cell_idx
    iblk => patch_2D%edges%cell_blk
    
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_edge_index,end_edge_index,je) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      DO je = start_edge_index, end_edge_index
        ! compute the normal derivative
        ! by the finite difference approximation
        ! (see Bonaventura and Ringler MWR 2005)
        grad_norm_psi_e(je,blockNo) =  &
          & (psi_c(iidx(je,blockNo,2),iblk(je,blockNo,2))-psi_c(iidx(je,blockNo,1),iblk(je,blockNo,1)))&
          & * grad_coeff(je,blockNo)
        
      END DO
      
    END DO
!ICON_OMP_END_PARALLEL_DO
    
  END SUBROUTINE grad_fd_norm_oce_2D_3D_sp
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !! Computes the discrete rotation at vertices in presence of boundaries as in the ocean setting.
  !! Computes in presence of boundaries the discrete rotation at vertices
  !! of triangle cells (centers of dual grid cells) from a vector field
  !! given by its components in the directions normal to triangle edges and
  !! takes the presence of boundaries into account.
  !!
  !! This sbr calculates the vorticity for mimetic discretization. A second one for the RBF-branch
  !! can be found below. The two sbr use the same approach to calculate the curl, but they differ
  !! in the calculation of the tangential velocity, which is only need at lateral boundaries. Mimetic
  !! does the tangential velocity calculate from velocity vector at vertices (vn_dual), whedgeOfVertex_index RBF uses
  !! a specific routine for that purpose.
  !!
  !! mpi note: the results is not synced. should be done by the calling method if necessary
  !!     vn, vn_dual must have been synced on level 2 (in_domain + 1)
!<Optimize:inUse>
  SUBROUTINE rot_vertex_ocean_3D( patch_3D, vn, vn_dual, p_op_coeff, rot_vec_v)
    !>
    !!
    TYPE(t_patch_3D ),TARGET, INTENT(in)      :: patch_3D
    REAL(wp), INTENT(in)                      :: vn(:,:,:)
    TYPE(t_cartesian_coordinates), INTENT(in) :: vn_dual(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),TARGET, INTENT(in) :: p_op_coeff
    REAL(wp), INTENT(inout)                   :: rot_vec_v(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)

    !Local variables
    !
    REAL(wp) :: z_vort_internal(n_zlev)
    REAL(wp) :: z_vort_boundary(n_zlev)
    REAL(wp) :: z_vt(4)  ! max boundary edges on a on a boundary vertex
    INTEGER :: start_level, end_level
    INTEGER :: vertexIndex, level, blockNo, vertexConnect
    INTEGER :: edge_index, edge_block, boundaryEdge_index, boundaryEdge_block, boundaryEdge_inVertex
    INTEGER :: il_v1, il_v2,ib_v1, ib_v2
    INTEGER :: start_index_v, end_index_v

    INTEGER, POINTER :: vertex_boundaryEdgeIndex(:,:,:,:), vertex_boundaryEdgeBlock(:,:,:,:), coeffs_VertexEdgeIndex(:,:,:,:)
    !REAL(wp), POINTER :: z_orientation(:,:,:,:)

    TYPE(t_subset_range), POINTER :: verts_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D          => patch_3D%p_patch_2D(1)
    verts_in_domain   => patch_2D%verts%in_domain
    start_level       = 1
    !set pointer that carry edge information
    vertex_boundaryEdgeIndex    => p_op_coeff%bnd_edge_idx
    vertex_boundaryEdgeBlock    => p_op_coeff%bnd_edge_blk
    coeffs_VertexEdgeIndex      => p_op_coeff%edge_idx
    !z_orientation    => p_op_coeff%orientation

    !In this loop vorticity at vertices is calculated
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_index_v,end_index_v,vertexIndex,end_level,vertexConnect,edge_index,edge_block,    &
!ICON_OMP z_vort_internal, level, z_vort_boundary, z_vt, boundaryEdge_inVertex, boundaryEdge_index, boundaryEdge_block, &
!ICON_OMP  il_v1,ib_v1,il_v2,ib_v2) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, blockNo, start_index_v, end_index_v)
      rot_vec_v(:,:,blockNo) = 0.0_wp
      DO vertexIndex = start_index_v, end_index_v
        end_level = patch_3D%p_patch_1d(1)%vertex_bottomLevel(vertexIndex, blockNo)
        z_vort_internal(:) = 0.0_wp

        DO vertexConnect = 1, patch_2D%verts%num_edges(vertexIndex,blockNo)

          ! get line and block indices of edge vertexConnect around vertex vertexIndex
          edge_index = patch_2D%verts%edge_idx(vertexIndex,blockNo,vertexConnect)
          edge_block = patch_2D%verts%edge_blk(vertexIndex,blockNo,vertexConnect)
          DO level = start_level, end_level

            !add contribution of normal velocity at edge (edgeOfVertex_index,edgeOfVertex_block) to rotation
            !IF ( v_base%lsm_e(edgeOfVertex_index,level,edgeOfVertex_block) == sea) THEN
            ! sea, sea_boundary, boundary (edges only), land_boundary, land =
            !  -2,      -1,         0,                  1,             2
            !Distinction between sea-lean-boundary is taken into account by coeffcients.
            !It is assumed here that vn is already zero at boundary edges.
            z_vort_internal(level) = z_vort_internal(level) + vn(edge_index,level,edge_block) * &
              & p_op_coeff%rot_coeff(vertexIndex,level,blockNo,vertexConnect)

          END DO ! level
        ENDDO ! verts%num_edges

        !Finalize vorticity calculation by closing the dual loop along boundary edges
        z_vort_boundary(start_level:end_level) = 0.0_wp
        z_vt(:) = 0.0_wp
        DO level = start_level, end_level
!           IF ( .NOT. (p_op_coeff%bnd_edges_per_vertex(vertexIndex,level,blockNo) == 0 .or. &
!                       p_op_coeff%bnd_edges_per_vertex(vertexIndex,level,blockNo) == 2 .or. &
!                       p_op_coeff%bnd_edges_per_vertex(vertexIndex,level,blockNo) == 4)) &
!             CALL finish("rot_vertex_ocean_3D", "wrong bnd_edges_per_vertex")
          DO boundaryEdge_inVertex = 1, p_op_coeff%bnd_edges_per_vertex(vertexIndex,level,blockNo)
            boundaryEdge_index = vertex_boundaryEdgeIndex(vertexIndex,level,blockNo,boundaryEdge_inVertex)
            boundaryEdge_block = vertex_boundaryEdgeBlock(vertexIndex,level,blockNo,boundaryEdge_inVertex)
            !calculate tangential velocity
            il_v1 = patch_2D%edges%vertex_idx(boundaryEdge_index,boundaryEdge_block,1)
            ib_v1 = patch_2D%edges%vertex_blk(boundaryEdge_index,boundaryEdge_block,1)
            il_v2 = patch_2D%edges%vertex_idx(boundaryEdge_index,boundaryEdge_block,2)
            ib_v2 = patch_2D%edges%vertex_blk(boundaryEdge_index,boundaryEdge_block,2)

            z_vt(boundaryEdge_inVertex)=   &
              & - DOT_PRODUCT(vn_dual(il_v1,level,ib_v1)%x,                                               &
              &     p_op_coeff%edge2vert_coeff_cc_t(boundaryEdge_index,level,boundaryEdge_block,1)%x)     &
              & + DOT_PRODUCT(vn_dual(il_v2,level,ib_v2)%x,                                               &
              &     p_op_coeff%edge2vert_coeff_cc_t(boundaryEdge_index,level,boundaryEdge_block,2)%x)

          ENDDO
          DO boundaryEdge_inVertex = 1, p_op_coeff%bnd_edges_per_vertex(vertexIndex,level,blockNo)

            z_vort_boundary(level) = z_vort_boundary(level) + &
              & z_vt(boundaryEdge_inVertex) * &
              &    p_op_coeff%rot_coeff(vertexIndex,level,blockNo, &
              &       coeffs_VertexEdgeIndex(vertexIndex,level,blockNo,boundaryEdge_inVertex))

          ENDDO ! boundaryEdge_inVertex
        ENDDO ! levels

        DO level = start_level, end_level
          !Final vorticity calculation
          !TODO ram
          !       rot_vec_v(vertexIndex,level,blockNo) = (z_vort_internal + z_vort_boundary(vertexIndex,level,blockNo)) / &
          !         & patch_2D%verts%dual_area(vertexIndex,blockNo)
          rot_vec_v(vertexIndex,level,blockNo) = z_vort_internal(level) + z_vort_boundary(level)

        END DO ! levels

      END DO ! vertexIndex
    END DO ! vertexBlock
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE rot_vertex_ocean_3D
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates vertical derivative for a vector that is located at cell center and at midelevel, i.e. at the center of a 3D prism.
  !!    start level has to be specifed, at end level value zero is assigned to vert. derivative   
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
!<Optimize:inUse>
  SUBROUTINE verticalDeriv_vec_midlevel_on_block(patch_3d, vec_in, vertDeriv_vec,start_level, &
    & blockNo, start_index, end_index)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_cartesian_coordinates), INTENT(in)        :: vec_in(nproma, n_zlev)
    INTEGER, INTENT(in)                              :: start_level
    INTEGER, INTENT(in)                              :: blockNo, start_index, end_index
    TYPE(t_cartesian_coordinates), INTENT(inout)       :: vertDeriv_vec(nproma, n_zlev)    ! out
    
    !Local variables
    INTEGER :: jk, jc!,jb
    ! REAL(wp), POINTER ::  prism_center_distance(:,:)
!     INTEGER :: end_level
    !-------------------------------------------------------------------------------
    ! prism_center_distance => patch_3D%p_patch_1D(1)%prism_center_dist_c  (:,:,blockNo)

    DO jc = start_index, end_index
!       end_level  = patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)
!      IF ( end_level >=min_dolic ) THEN
        DO jk = start_level,patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo) - 1
          vertDeriv_vec(jc,jk)%x &
          & = (vec_in(jc,jk-1)%x - vec_in(jc,jk)%x)  & !/ prism_center_distance(jc,jk)
              & * patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,jk,blockNo)
              
        END DO    
        ! vertDeriv_vec(jc,end_level)%x = 0.0_wp ! this is not needed 
!      ENDIF
    END DO
    
  END SUBROUTINE verticalDeriv_vec_midlevel_on_block
  !------------------------------------------------------------------------- 
  
  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates vertical derivative for a scalar that is located at cell center and at midelevel, i.e. at the center of a 3D prism.
  !!    start level has to be specifed, at end level value zero is assigned to vert. derivative 
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
!<Optimize:inUse>
  SUBROUTINE verticalDeriv_scalar_midlevel_on_block(patch_3d, scalar_in, vertDeriv_scalar, start_level, &
    & blockNo, start_index, end_index)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    REAL(wp), INTENT(in)                             :: scalar_in(nproma, n_zlev)
    INTEGER, INTENT(in)                              :: start_level
    INTEGER, INTENT(in)                              :: blockNo, start_index, end_index
    REAL(wp), INTENT(inout)                          :: vertDeriv_scalar(nproma, n_zlev)    ! out

    !Local variables
    INTEGER :: jk, jc!,jb
    ! REAL(wp), POINTER ::  prism_center_distance(:,:)
!     INTEGER :: end_level
    !-------------------------------------------------------------------------------
    ! prism_center_distance => patch_3D%p_patch_1D(1)%prism_center_dist_c  (:,:,blockNo)

    DO jc = start_index, end_index
!       end_level  = patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)
!      IF ( end_level >=min_dolic ) THEN
        DO jk = start_level,patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo) - 1
          vertDeriv_scalar(jc,jk) &
          & = (scalar_in(jc,jk) - scalar_in(jc,jk+1))  & 
              & * patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,jk,blockNo)

        END DO
        ! vertDeriv_vec(jc,end_level)%x = 0.0_wp ! this is not needed
!      ENDIF
    END DO
   !CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vertDeriv_scalar(:,:))
  END SUBROUTINE verticalDeriv_scalar_midlevel_on_block
  !-------------------------------------------------------------------------
  
   !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE verticalDiv_scalar_midlevel( patch_3d, scalar_in, vertDiv_scalar, subset_range)
    TYPE(t_patch_3d), TARGET, INTENT(in) :: patch_3D
    REAL(wp), INTENT(in)                 :: scalar_in(nproma, n_zlev, patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)              :: vertDiv_scalar(nproma, n_zlev, patch_3D%p_patch_2D(1)%alloc_cell_blocks)    ! out
    TYPE(t_subset_range), TARGET, OPTIONAL :: subset_range

    INTEGER :: blockNo,start_level,jk
    INTEGER :: start_cell_index, end_cell_index
    TYPE(t_subset_range), POINTER :: cells_in_domain

   !-----------------------------------------------------------------------
    IF (PRESENT(subset_range)) THEN
      cells_in_domain => subset_range
    ELSE
      cells_in_domain => patch_3D%p_patch_2D(1)%cells%in_domain
    ENDIF    
    start_level=1
    !-----------------------------------------------------------------------
    
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_cell_index,end_cell_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
    
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      vertDiv_scalar(:,:,blockNo)=0.0_wp
      CALL verticalDiv_scalar_midlevel_on_block(patch_3d, scalar_in(:,:,blockNo), vertDiv_scalar(:,:,blockNo), start_level, &
        & blockNo, start_cell_index, end_cell_index)
      
    END DO
!ICON_OMP_END_PARALLEL_DO

  !CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vertDiv_scalar)
  END SUBROUTINE verticalDiv_scalar_midlevel
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates vertical derivative for a scalar that is located at cell center and at midelevel, i.e. at the center of a 3D prism.
  !!    start level has to be specifed, at end level value zero is assigned to vert. derivative 
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
!<Optimize:inUse>
  SUBROUTINE verticalDiv_scalar_midlevel_on_block(patch_3d, scalar_in, vertDiv_scalar, start_level, &
    & blockNo, start_index, end_index)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    REAL(wp), INTENT(in)                             :: scalar_in(nproma, n_zlev)
    INTEGER, INTENT(in)                              :: start_level
    INTEGER, INTENT(in)                              :: blockNo, start_index, end_index
    REAL(wp), INTENT(inout)                          :: vertDiv_scalar(nproma, n_zlev)    ! out

    !Local variables
    INTEGER :: jk, jc!,jb
    ! REAL(wp), POINTER ::  prism_center_distance(:,:)
!     INTEGER :: end_level
    !-------------------------------------------------------------------------------
    ! prism_center_distance => patch_3D%p_patch_1D(1)%prism_center_dist_c  (:,:,blockNo)

    DO jc = start_index, end_index
!       end_level  = patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)
!      IF ( end_level >=min_dolic ) THEN
        DO jk = start_level,patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo) - 1
          vertDiv_scalar(jc,jk) &
          & = (scalar_in(jc,jk) - scalar_in(jc,jk+1))  & !/ prism_center_distance(jc,jk)
              & * patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,blockNo)

        END DO
        ! vertDeriv_vec(jc,end_level)%x = 0.0_wp ! this is not needed
!      ENDIF
    END DO
     !CALL sync_patch_array(sync_c, patch_3D%p_patch_2D(1), vertDiv_scalar)
  END SUBROUTINE verticalDiv_scalar_midlevel_on_block
  !-------------------------------------------------------------------------
 
 
  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE smooth_onCells_3D( patch_3D, in_value, out_value, smooth_weights)

    TYPE(t_patch_3D ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)          :: in_value(:,:,:)  ! dim: (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)       :: out_value(:,:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(in)          :: smooth_weights(1:2) ! 1st=weight for this cell, 2nd=weight for the some of the neigbors


    INTEGER :: max_connectivity, blockNo, start_index,end_index, jc, level, neigbor, neigbor_index,neigbor_block 
    REAL(wp) :: numberOfNeigbors, neigbors_weight
    TYPE(t_subset_range), POINTER :: cells_inDomain
    !-----------------------------------------------------------------------
    cells_inDomain => patch_3D%p_patch_2D(1)%cells%in_domain
    max_connectivity = patch_3D%p_patch_2D(1)%cells%max_connectivity
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, level, neigbor, neigbor_index,neigbor_block, &
!ICON_OMP numberOfNeigbors, neigbors_weight) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_inDomain%start_block, cells_inDomain%end_block
      CALL get_index_range(cells_inDomain, blockNo, start_index, end_index)
      out_value(:,:,blockNo) = 0.0_wp
      
      DO jc = start_index, end_index
        DO level = 1, patch_3D%p_patch_1d(1)%dolic_c(jc, blockNo)

          ! calculate how many sea neigbors we have 
          numberOfNeigbors = 0.0_wp            
          ! now compute out_value, out_value at this point is zeroe
          DO neigbor = 1, max_connectivity
            neigbor_index = patch_3D%p_patch_2D(1)%cells%neighbor_idx(jc,blockNo,neigbor)
            neigbor_block = patch_3D%p_patch_2D(1)%cells%neighbor_blk(jc,blockNo,neigbor)

            IF (patch_3D%p_patch_1d(1)%dolic_c(neigbor_index, neigbor_block) >= level) THEN
              out_value(jc,level,blockNo) = out_value(jc,level,blockNo) + &
                & in_value(neigbor_index,level,neigbor_block)
                ! & * patch_3D%p_patch_2D(1)%cells%area(neigbor_index,neigbor_block) 
              numberOfNeigbors = numberOfNeigbors + 1.0_wp
            ENDIF
          ENDDO
          
          neigbors_weight = smooth_weights(2) / numberOfNeigbors          
          out_value(jc,level,blockNo) = &
            &  out_value(jc,level,blockNo) * neigbors_weight + &
            &  in_value(jc,level,blockNo) * smooth_weights(1)

        END DO ! level
      END DO 
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE smooth_onCells_3D
  !-------------------------------------------------------------------------
 
  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE smooth_onCells_2D( patch_3D, in_value, out_value, smooth_weights)

    TYPE(t_patch_3D ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)          :: in_value(:,:)  ! dim: (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)       :: out_value(:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(in)          :: smooth_weights(1:2) ! 1st=weight for this cell, 2nd=weight for the some of the neigbors


    INTEGER :: max_connectivity, blockNo, start_index,end_index, jc, level, neigbor, neigbor_index,neigbor_block
    REAL(wp) :: numberOfNeigbors, neigbors_weight
    TYPE(t_subset_range), POINTER :: cells_inDomain
    !-----------------------------------------------------------------------
    cells_inDomain => patch_3D%p_patch_2D(1)%cells%in_domain
    max_connectivity = patch_3D%p_patch_2D(1)%cells%max_connectivity

!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, jc, level, neigbor, neigbor_index,neigbor_block, &
!ICON_OMP numberOfNeigbors, neigbors_weight) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_inDomain%start_block, cells_inDomain%end_block
      CALL get_index_range(cells_inDomain, blockNo, start_index, end_index)
      out_value(:,blockNo) = 0.0_wp

      DO jc = start_index, end_index
        DO level = 1, MAX(patch_3D%p_patch_1d(1)%dolic_c(jc, blockNo),1)

          ! calculate how many sea negbors we have 
          numberOfNeigbors = 0.0_wp

          ! now compute out_value, out_value at this point is zeroe
          DO neigbor = 1, max_connectivity
            neigbor_index = patch_3D%p_patch_2D(1)%cells%neighbor_idx(jc,blockNo,neigbor)
            neigbor_block = patch_3D%p_patch_2D(1)%cells%neighbor_blk(jc,blockNo,neigbor)

            IF (patch_3D%p_patch_1d(1)%dolic_c(neigbor_index, neigbor_block) >= level) THEN
              out_value(jc,blockNo) = out_value(jc,blockNo) + &
                & in_value(neigbor_index,neigbor_block) 
                ! & * patch_3D%p_patch_2D(1)%cells%area(neigbor_index,neigbor_block)
              numberOfNeigbors = numberOfNeigbors + 1.0_wp
            ENDIF
          ENDDO
          
          neigbors_weight = smooth_weights(2) / numberOfNeigbors

          out_value(jc,blockNo) = &
            & out_value(jc,blockNo) * neigbors_weight  + &
            & in_value(jc,blockNo)  * smooth_weights(1)

        END DO ! 1 level
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE smooth_onCells_2D
  !-------------------------------------------------------------------------

  !---------------------------------------------------------------------------------
  !>
  !!
  !!  Calculation of total fluid thickness at cell centers and surface elevation at
  !!  cell edges from prognostic surface height at cell centers. We use height at
  !!  old timelevel "n"
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
!<Optimize:inUse>
  SUBROUTINE calculate_thickness( patch_3D, ocean_state, p_ext_data, operators_coefficients, solvercoeff_sp, inTopCellThickness)
    !SUBROUTINE calculate_thickness( p_patch_3D, ocean_state, p_ext_data, ice_hi)
    !
    ! patch_2D on which computation is performed
    TYPE(t_patch_3D ),TARGET, INTENT(in)   :: patch_3D
    !
    ! Type containing ocean state
    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
    !
    ! Type containing external data
    TYPE(t_external_data), TARGET, INTENT(in) :: p_ext_data
    !REAL(wp), INTENT(IN)                      :: ice_hi(nproma,1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_operator_coeff), INTENT(in)      :: operators_coefficients
    TYPE(t_solvercoeff_singleprecision), INTENT(inout) :: solvercoeff_sp
    REAL(wp), OPTIONAL :: inTopCellThickness(:,:)
    
    !  local variables
    INTEGER :: cell_StartIndex, cell_EndIndex
    INTEGER :: edge_StartIndex, edge_EndIndex
    INTEGER :: jc, blockNo, je, level
    INTEGER :: thislevel, levelabove, levelbelow, level2below, cell_levels
    
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    REAL(wp)           :: z_dist_e_c1, z_dist_e_c2
    TYPE(t_subset_range), POINTER :: all_cells, all_edges, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    
    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    
    REAL(wp) :: top_vn_1, top_vn_2, integrated_vn
    REAL(wp), POINTER :: top_coeffs(:,:,:), integrated_coeffs(:,:,:), sum_to_2D_coeffs(:,:,:)
    REAL(wp), POINTER :: cell_thickness(:,:,:), edge_thickness(:,:,:)
    REAL(wp), POINTER :: inv_cell_thickness(:,:,:), inv_edge_thickness(:,:,:)
    REAL(wp), POINTER :: inv_prisms_center_distance(:,:,:), inv_edgefaces_middle_distance(:,:,:)
    REAL(wp)  :: cell_thickness_1, cell_thickness_2
    !-------------------------------------------------------------------------------
    ! pointers for the ppm vertical transport
    TYPE(t_verticaladvection_ppm_coefficients), POINTER :: vertadvppm
    !-------------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    patch_2D            => patch_3D%p_patch_2D(1)
    all_cells           => patch_2D%cells%ALL
    all_edges           => patch_2D%edges%ALL
    edges_in_domain     => patch_2D%edges%in_domain
    cell_thickness     => patch_3D%p_patch_1d(1)%prism_thick_c
    inv_cell_thickness => patch_3D%p_patch_1d(1)%inv_prism_thick_c
    inv_prisms_center_distance => patch_3D%p_patch_1d(1)%inv_prism_center_dist_c
    edge_thickness     => patch_3D%p_patch_1d(1)%prism_thick_e
    inv_edge_thickness => patch_3D%p_patch_1d(1)%inv_prism_thick_e
    inv_edgefaces_middle_distance => patch_3D%p_patch_1d(1)%inv_prism_center_dist_e
    
    
    ! already done after update fluxes
    !    CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_prog(nold(1))%h)
    
    !Step 1: calculate cell-located variables for 2D and 3D case
    !For 3D and for SWE thick_c contains thickness of fluid column


    !Update prism thickness. The prism-thickness below the surface is
    !not updated it is initialized in construct_hydro_ocean_diag
    !with z-coordinate-thickness.
    !1) Thickness at cells
!ICON_OMP_PARALLEL PRIVATE(top_coeffs, integrated_coeffs, sum_to_2D_coeffs)
    IF (PRESENT(inTopCellThickness)) THEN
!ICON_OMP_DO PRIVATE(cell_StartIndex, cell_EndIndex, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, cell_StartIndex, cell_EndIndex)
        DO jc = cell_StartIndex, cell_EndIndex
          IF ( patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo) > 0 ) THEN
            cell_thickness(jc,1,blockNo) = inTopCellThickness(jc,blockNo)
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO

    ELSE
!ICON_OMP_DO PRIVATE(cell_StartIndex, cell_EndIndex, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, cell_StartIndex, cell_EndIndex)
        DO jc = cell_StartIndex, cell_EndIndex
          IF ( patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo) > 0 ) THEN
            cell_thickness(jc,1,blockNo) = &
              & patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(jc,1,blockNo) + &
              & ocean_state%p_prog(nold(1))%h(jc,blockNo)
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO

    ENDIF
    IF ( iswm_oce /= 1 ) THEN  !  3D case      
!ICON_OMP_DO PRIVATE(cell_StartIndex, cell_EndIndex, jc, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, cell_StartIndex, cell_EndIndex)
      DO jc = cell_StartIndex, cell_EndIndex
        IF ( patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo) > 0 ) THEN
          
          patch_3D%p_patch_1d(1)%prism_center_dist_c(jc,2,blockNo) = 0.5_wp * &
            & (cell_thickness(jc,1,blockNo) + cell_thickness(jc,2,blockNo))
          
          patch_3D%p_patch_1d(1)%prism_volume(jc,1,blockNo) = cell_thickness(jc,1,blockNo) * &
            & patch_2D%cells%area(jc,blockNo)
          
          inv_cell_thickness(jc,1,blockNo) = 1.0_wp / cell_thickness(jc,1,blockNo)
          
          inv_prisms_center_distance(jc,2,blockNo) = &
            & 1.0_wp / patch_3D%p_patch_1d(1)%prism_center_dist_c(jc,2,blockNo)
          
          ocean_state%p_diag%thick_c(jc,blockNo) = ocean_state%p_prog(nold(1))%h(jc,blockNo) + patch_3D%column_thick_c(jc,blockNo)
          
          patch_3D%p_patch_1d(1)%depth_cellmiddle(jc,1,blockNo) = cell_thickness(jc,1,blockNo) * 0.5_wp
          patch_3D%p_patch_1d(1)%depth_cellinterface(jc,2,blockNo) = cell_thickness(jc,1,blockNo)
          
        ENDIF
        
        DO level=2, patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)
          patch_3D%p_patch_1d(1)%depth_cellmiddle(jc,level,blockNo) = &
            & patch_3D%p_patch_1d(1)%depth_cellinterface(jc,level,blockNo) + cell_thickness(jc,level,blockNo) * 0.5_wp
          patch_3D%p_patch_1d(1)%depth_cellinterface(jc,level+1,blockNo) = &
            & patch_3D%p_patch_1d(1)%depth_cellinterface(jc,level,blockNo) + cell_thickness(jc,level,blockNo)
        ENDDO
        
      END DO
    END DO
!ICON_OMP_END_DO
    ENDIF
    !----------------------------------------------------------------------------------------
    IF ( iswm_oce == 1 ) THEN  !  SWM
!ICON_OMP_DO PRIVATE(cell_StartIndex, cell_EndIndex, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, cell_StartIndex, cell_EndIndex)
        !calculate for each fluid colum the total depth, i.e.
        !from bottom boundary to surface height, i.e. using individual bathymetry for SWM
        DO jc = cell_StartIndex, cell_EndIndex
          IF(patch_3D%lsm_c(jc,1,blockNo) <= sea_boundary)THEN
            
            ocean_state%p_diag%thick_c(jc,blockNo) = ocean_state%p_prog(nold(1))%h(jc,blockNo)&
              & - p_ext_data%oce%bathymetry_c(jc,blockNo)
            !        &  - ice_hi(jc,1,blockNo) 
            cell_thickness(jc,1,blockNo)=ocean_state%p_diag%thick_c(jc,blockNo)
          ELSE
            ocean_state%p_diag%thick_c(jc,blockNo) = 0.0_wp
          ENDIF
        END DO
      END DO!write(*,*)'bathymetry',maxval(p_ext_data%oce%bathymetry_c),minval(p_ext_data%oce%bathymetry_c)
!ICON_OMP_END_DO
      write(*,*)'bathymetry cell',&
      &maxval(p_ext_data%oce%bathymetry_c),minval(p_ext_data%oce%bathymetry_c),&
      &maxval(ocean_state%p_diag%thick_c),minval(ocean_state%p_diag%thick_c),&
      &maxval(ocean_state%p_prog(nold(1))%h),minval(ocean_state%p_prog(nold(1))%h)
      
      
      !Step 2: calculate edge-located variables for 2D and 3D case from respective cell variables
      !For SWE : thick_e = thickness of fluid column at edges
      !         h_e     = surface elevation at edges, without depth of first layer
!ICON_OMP_DO PRIVATE(edge_StartIndex, edge_EndIndex, je, il_c1, ib_c1, il_c2, ib_c2, &
!ICON_OMP z_dist_e_c1, z_dist_e_c2) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, edge_StartIndex, edge_EndIndex)
        DO je = edge_StartIndex, edge_EndIndex
          
          il_c1 = patch_2D%edges%cell_idx(je,blockNo,1)
          ib_c1 = patch_2D%edges%cell_blk(je,blockNo,1)
          il_c2 = patch_2D%edges%cell_idx(je,blockNo,2)
          ib_c2 = patch_2D%edges%cell_blk(je,blockNo,2)
          
          z_dist_e_c1 = 0.5_wp!z_dist_e_c1=p_patch%edges%edge_cell_length(je,blockNo,1)
          z_dist_e_c2 = 0.5_wp!z_dist_e_c2=p_patch%edges%edge_cell_length(je,blockNo,2)
          
          IF(patch_3D%lsm_e(je,1,blockNo) <= sea_boundary)THEN
         
            ocean_state%p_diag%thick_e(je,blockNo) = ( z_dist_e_c1*ocean_state%p_diag%thick_c(il_c1,ib_c1)&
              & +   z_dist_e_c2*ocean_state%p_diag%thick_c(il_c2,ib_c2) )&
              & /(z_dist_e_c1+z_dist_e_c2)
            
            ocean_state%p_diag%h_e(je,blockNo) = ( z_dist_e_c1*ocean_state%p_prog(nold(1))%h(il_c1,ib_c1)&
              & +   z_dist_e_c2*ocean_state%p_prog(nold(1))%h(il_c2,ib_c2) )&
              & /(z_dist_e_c1+z_dist_e_c2)
              
            patch_3d%p_patch_1d(1)%prism_thick_e(je,1,blockNo)=ocean_state%p_diag%thick_e(je,blockNo)
          ELSE
            ocean_state%p_diag%h_e(je,blockNo) = 0.0_wp
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO
      write(*,*)'bathymetry edge',&
      &maxval(ocean_state%p_diag%thick_e),minval(ocean_state%p_diag%thick_e),&
      &maxval(ocean_state%p_diag%h_e),minval(ocean_state%p_diag%h_e)      

!ICON_OMP_MASTER
      CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%thick_e)
      CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%h_e)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER

      !2) Thickness at edges
!ICON_OMP_DO PRIVATE(edge_StartIndex, edge_EndIndex, je) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, blockNo, edge_StartIndex, edge_EndIndex)
        DO je = edge_StartIndex, edge_EndIndex
          IF ( patch_3D%p_patch_1d(1)%dolic_e(je,blockNo) > 0 ) THEN
            
            edge_thickness(je,1,blockNo) = ocean_state%p_diag%thick_e(je,blockNo)
            !patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_e(je,1,blockNo) + ocean_state%p_diag%h_e(je,blockNo)
            
            inv_edge_thickness(je,1,blockNo)= 1.0_wp / edge_thickness(je,1,blockNo)
            
            !Not possible for SWE
            !inv_edgefaces_middle_distance(je,2,blockNo) = 2.0_wp / &
            !  & (edge_thickness(je,1,blockNo) + edge_thickness(je,2,blockNo))
            
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO
      
    !----------------------------------------------------------------------------------------
    ELSE!IF 3D model
      
      !Step 2: calculate edge-located variables for 2D and 3D case from respective cell variables
      !For 3D: thick_e = thickness of fluid column at edges
      !         h_e     = surface elevation at edges, without depth of first layer
!ICON_OMP_DO PRIVATE(edge_StartIndex, edge_EndIndex, je, il_c1, ib_c1, il_c2, ib_c2, &
!ICON_OMP z_dist_e_c1, z_dist_e_c2) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockNo, edge_StartIndex, edge_EndIndex)
        DO je = edge_StartIndex, edge_EndIndex
          
          il_c1 = patch_2D%edges%cell_idx(je,blockNo,1)
          ib_c1 = patch_2D%edges%cell_blk(je,blockNo,1)
          il_c2 = patch_2D%edges%cell_idx(je,blockNo,2)
          ib_c2 = patch_2D%edges%cell_blk(je,blockNo,2)
          
          !z_dist_e_c1 = 0.5_wp!z_dist_e_c1=p_patch%edges%edge_cell_length(je,blockNo,1)
          !z_dist_e_c2 = 0.5_wp!z_dist_e_c2=p_patch%edges%edge_cell_length(je,blockNo,2)
          
          IF ( patch_3D%p_patch_1d(1)%dolic_e(je,blockNo) > 0 ) THEN
            
!             ocean_state%p_diag%h_e(je,blockNo) = ( z_dist_e_c1 * ocean_state%p_prog(nold(1))%h(il_c1,ib_c1)   &
!               & +   z_dist_e_c2 * ocean_state%p_prog(nold(1))%h(il_c2,ib_c2) )                                &
!               & /(z_dist_e_c1 + z_dist_e_c2)              
            ocean_state%p_diag%h_e(je,blockNo) = ( 0.5_wp * ocean_state%p_prog(nold(1))%h(il_c1,ib_c1)   &
              & +   0.5_wp * ocean_state%p_prog(nold(1))%h(il_c2,ib_c2) )
            
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO

!ICON_OMP_MASTER      
      CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%h_e)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER

      !2) Thickness at edges
!ICON_OMP_DO PRIVATE(edge_StartIndex, edge_EndIndex, je) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, blockNo, edge_StartIndex, edge_EndIndex)
        DO je = edge_StartIndex, edge_EndIndex
          IF ( patch_3D%p_patch_1d(1)%dolic_e(je,blockNo) > 0 ) THEN
            
            edge_thickness(je,1,blockNo)&
              & = patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_e(je,1,blockNo) + ocean_state%p_diag%h_e(je,blockNo)
            
            inv_edge_thickness(je,1,blockNo)= 1.0_wp / edge_thickness(je,1,blockNo)
            
            inv_edgefaces_middle_distance(je,2,blockNo) = 2.0_wp / &
              & (edge_thickness(je,1,blockNo) + edge_thickness(je,2,blockNo))
            
            ocean_state%p_diag%thick_e(je,blockNo) = ocean_state%p_diag%h_e(je,blockNo) &
              & + patch_3D%column_thick_e(je,blockNo)
            
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO
      !---------------------------------------------------------------------
    ENDIF  ! shallow water model
    
    !---------------------------------------------------------------------
    ! update the coefficients for the edge2edge_viacell_1D fast operator
    top_coeffs        => operators_coefficients%edge2edge_viacell_coeff_top
    integrated_coeffs => operators_coefficients%edge2edge_viacell_coeff_integrated
    sum_to_2D_coeffs  => operators_coefficients%edge2edge_viacell_coeff_all
!ICON_OMP_DO PRIVATE(edge_StartIndex, edge_EndIndex, je, cell_1_index, cell_1_block,  &
!ICON_OMP cell_2_index, cell_2_block, edge_1_1_index, edge_1_2_index, edge_1_3_index, &
!ICON_OMP edge_2_1_index, edge_2_2_index, edge_2_3_index, edge_1_1_block, edge_1_2_block, &
!ICON_OMP edge_1_3_block, edge_2_1_block, edge_2_2_block, edge_2_3_block) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, edge_StartIndex, edge_EndIndex)
      DO je = edge_StartIndex, edge_EndIndex

        IF ( patch_3D%p_patch_1d(1)%dolic_e(je,blockNo) > 0 ) THEN

          ! get the two cells of the edge
          cell_1_index = patch_2D%edges%cell_idx(je,blockNo,1)
          cell_1_block = patch_2D%edges%cell_blk(je,blockNo,1)
          cell_2_index = patch_2D%edges%cell_idx(je,blockNo,2)
          cell_2_block = patch_2D%edges%cell_blk(je,blockNo,2)
!           cell_thickness_1 = cell_thickness(cell_1_index, 1, cell_1_block)
!           cell_thickness_2 = cell_thickness(cell_2_index, 1, cell_2_block)

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

!           sum_to_2D_coeffs(1, je, blockNo) = top_coeffs(1, je, blockNo) * cell_thickness_1 + integrated_coeffs(1, je, blockNo)
!           sum_to_2D_coeffs(2, je, blockNo) = top_coeffs(2, je, blockNo) * cell_thickness_1 + integrated_coeffs(2, je, blockNo)
!           sum_to_2D_coeffs(3, je, blockNo) = top_coeffs(3, je, blockNo) * cell_thickness_1 + integrated_coeffs(3, je, blockNo)
!
!           sum_to_2D_coeffs(4, je, blockNo) = top_coeffs(4, je, blockNo) * cell_thickness_2 + integrated_coeffs(4, je, blockNo)
!           sum_to_2D_coeffs(5, je, blockNo) = top_coeffs(5, je, blockNo) * cell_thickness_2 + integrated_coeffs(5, je, blockNo)
!           sum_to_2D_coeffs(6, je, blockNo) = top_coeffs(6, je, blockNo) * cell_thickness_2 + integrated_coeffs(6, je, blockNo)

          ! now we ues the edge thickeness instead of the cell thickeness for the top layer

          sum_to_2D_coeffs(1, je, blockNo) = &
            & top_coeffs(1, je, blockNo) * edge_thickness(edge_1_1_index, 1, edge_1_1_block) &
            & + integrated_coeffs(1, je, blockNo)

          sum_to_2D_coeffs(2, je, blockNo) = &
            & top_coeffs(2, je, blockNo) * edge_thickness(edge_1_2_index, 1, edge_1_2_block) &
            & + integrated_coeffs(2, je, blockNo)

          sum_to_2D_coeffs(3, je, blockNo) = &
            & top_coeffs(3, je, blockNo) * edge_thickness(edge_1_3_index, 1, edge_1_3_block) &
            & + integrated_coeffs(3, je, blockNo)

           sum_to_2D_coeffs(4, je, blockNo) = &
            & top_coeffs(4, je, blockNo) * edge_thickness(edge_2_1_index, 1, edge_2_1_block) &
            & + integrated_coeffs(4, je, blockNo)

          sum_to_2D_coeffs(5, je, blockNo) = &
            & top_coeffs(5, je, blockNo) * edge_thickness(edge_2_2_index, 1, edge_2_2_block) &
            & + integrated_coeffs(5, je, blockNo)

          sum_to_2D_coeffs(6, je, blockNo) = &
            & top_coeffs(6, je, blockNo) * edge_thickness(edge_2_3_index, 1, edge_2_3_block) &
            & + integrated_coeffs(6, je, blockNo)

        ENDIF
      END DO
    END DO ! blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!ICON_OMP_END_DO

    
    IF (select_solver == select_restart_mixedprecision_gmres) THEN
!ICON_OMP WORKSHARE
      solvercoeff_sp%edge_thickness(:,:)  = REAL(ocean_state%p_diag%thick_e(:,:), sp)
      solvercoeff_sp%cell_thickness(:,:)  = REAL(ocean_state%p_diag%thick_c(:,:), sp)
!ICON_OMP_END_WORKSHARE      
    ENDIF
    !-------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------
    ! update the coefficients for the upwind_vflux_ppm_fast vertical advection
!ICON_OMP_DO PRIVATE(cell_StartIndex, cell_EndIndex, vertAdvPPM, jc, cell_levels, thisLevel, levelAbove, &
!ICON_OMP  levelBelow, level2Below   ) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, cell_StartIndex, cell_EndIndex)
      vertadvppm => operators_coefficients%verticaladvectionppmcoeffs(blockNo)
      DO jc = cell_StartIndex, cell_EndIndex
        
        cell_levels = patch_3D%p_patch_1d(1)%dolic_c(jc,blockNo)
        
        thislevel  = 1
        levelbelow = 2
        IF ( cell_levels >= levelbelow ) THEN
          
          vertadvppm%cellheightratio_this_tobelow(jc, thislevel) = &
            & cell_thickness(jc, thislevel, blockNo) / cell_thickness(jc, levelbelow, blockNo)
          
          vertadvppm%cellheightratio_this_tothisbelow(jc, thislevel) = &
            & cell_thickness(jc, thislevel, blockNo) / &
            & (cell_thickness(jc, thislevel, blockNo) + cell_thickness(jc, levelbelow, blockNo))
          
          vertadvppm%cellheight_2xbelow_x_ratiothis_tothisbelow(jc,thislevel) = &
            & 2._wp * cell_thickness(jc,levelbelow, blockNo) * &
            & vertadvppm%cellheightratio_this_tothisbelow(jc, thislevel)
          
        ENDIF
        
        thislevel  = 2
        levelabove = 1
        levelbelow = 3
        level2below = 4
        IF ( cell_levels >= levelbelow ) THEN
          
          
          vertadvppm%cellheightratio_this_tothisabovebelow(jc,thislevel) = &
            & cell_thickness(jc, thislevel ,blockNo) / &
            & (cell_thickness(jc,levelabove,blockNo) + cell_thickness(jc,thislevel,blockNo)    &
            & + cell_thickness(jc,levelbelow,blockNo))
          
          vertadvppm%cellheightratio_2xaboveplusthis_tothisbelow(jc,thislevel) = &
            & (2._wp * cell_thickness(jc,levelabove,blockNo) + cell_thickness(jc,thislevel,blockNo))     &
            & / (cell_thickness(jc,levelbelow,blockNo) + cell_thickness(jc,thislevel,blockNo))
          
          vertadvppm%cellheightratio_2xbelowplusthis_tothisabove(jc,thislevel) = &
            & + (cell_thickness(jc,thislevel,blockNo) + 2._wp * cell_thickness(jc,levelbelow,blockNo))   &
            & / (cell_thickness(jc,levelabove,blockNo) + cell_thickness(jc,thislevel,blockNo))
          
          vertadvppm%cellheightratio_thisabove_to2xthisplusbelow(jc,thislevel) =                         &
            & (cell_thickness(jc,levelabove,blockNo) + cell_thickness(jc,thislevel,blockNo))            &
            & / (2._wp*cell_thickness(jc,thislevel,blockNo) + cell_thickness(jc,levelbelow,blockNo))
          
          vertadvppm%cellheightratio_thisbelow_to2xthisplusabove(jc,thislevel) =                 &
            & (cell_thickness(jc,levelbelow,blockNo) + cell_thickness(jc,thislevel,blockNo))                  &
            & / (2._wp*cell_thickness(jc,thislevel,blockNo) + cell_thickness(jc,levelabove,blockNo))
          ! = 1 / cellHeightRatio_2xBelowplusThis_toThisAbove(levelBelow)
          
        ENDIF
        
        IF ( cell_levels >= level2below ) THEN
          vertadvppm%cellheight_inv_thisabovebelow2below(jc,thislevel) =                                  &
            & 1._wp / (cell_thickness(jc,levelabove,blockNo) + cell_thickness(jc,thislevel,blockNo)       &
            & + cell_thickness(jc,levelbelow,blockNo) + cell_thickness(jc,level2below,blockNo))
        ENDIF
        
      END DO
    END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL    
    !-------------------------------------------------------------------------
    
    
    !---------Debug Diagnostics-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('heightRelQuant: h_e'    ,ocean_state%p_diag%h_e        ,str_module,idt_src, &
      & in_subset=patch_2D%edges%owned)
    idt_src=3
    CALL dbg_print('heightRelQuant: h_c'    ,ocean_state%p_prog(nold(1))%h ,str_module,idt_src, &
      & in_subset=patch_2D%cells%owned)
    CALL dbg_print('heightRelQuant: thick_c',ocean_state%p_diag%thick_c    ,str_module,idt_src, &
      & in_subset=patch_2D%cells%owned)
    CALL dbg_print('heightRelQuant: thick_e',ocean_state%p_diag%thick_e    ,str_module,idt_src, &
      & in_subset=patch_2D%edges%owned)
    CALL dbg_print('depth_CellMiddle', &
      & patch_3D%p_patch_1d(1)%depth_cellmiddle   ,str_module,idt_src, &
      & in_subset=patch_2D%cells%owned)
    CALL dbg_print('depth_CellInterface', &
      & patch_3D%p_patch_1d(1)%depth_cellinterface   ,str_module,idt_src, &
      & in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------
  END SUBROUTINE calculate_thickness
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE check_cfl_horizontal(normal_velocity,inv_dual_edge_length,timestep,edges,threshold, &
    & cfl_diag, stop_on_violation, output)
    REAL(wp),POINTER :: normal_velocity(:,:,:)
    TYPE(t_subset_range)  :: edges
    REAL(wp), INTENT(in)  :: inv_dual_edge_length(:,:), threshold, timestep
    REAL(wp), POINTER :: cfl_diag(:,:,:)
    LOGICAL, INTENT(in)   :: stop_on_violation, output
    REAL(wp), POINTER :: cfl(:,:,:)
    
    INTEGER :: je,level,blockNo,start_index,end_index
    
    ALLOCATE(cfl(LBOUND(normal_velocity,1):UBOUND(normal_velocity,1),&
      & LBOUND(normal_velocity,2):UBOUND(normal_velocity,2),&
      & LBOUND(normal_velocity,3):UBOUND(normal_velocity,3)))
    cfl = 0.0_wp
    
    DO blockNo = edges%start_block, edges%end_block
      CALL get_index_range(edges,blockNo,start_index,end_index)
      DO je = start_index,end_index
        DO level = 1, n_zlev
          cfl(je,level,blockNo) = ABS(dtime*normal_velocity(je,level,blockNo)*inv_dual_edge_length(je,blockNo))
        END DO
      END DO
    END DO
    IF (output) THEN
      IF (ASSOCIATED(cfl_diag)) THEN
        cfl_diag(:,:,:) = cfl(:,:,:)
      ELSE
        CALL finish('check_cfl_vertical','cfl_diag pointer for output NOT ASSOCIATED')
      ENDIF
    ENDIF
    
    CALL dbg_print('check horiz. CFL',cfl ,str_module,3,in_subset=edges)
    
    CALL check_cfl_threshold(MAXVAL(cfl),threshold,'horz',stop_on_violation)
    
    DEALLOCATE(cfl)
  END SUBROUTINE check_cfl_horizontal
  
  SUBROUTINE check_cfl_vertical(vertical_velocity, thicknesses, timestep, cells, threshold, &
    & cfl_diag, stop_on_violation, output)
    REAL(wp),POINTER :: vertical_velocity(:,:,:), thicknesses(:,:,:)
    REAL(wp), INTENT(in) :: timestep, threshold
    TYPE(t_subset_range) :: cells
    REAL(wp), POINTER :: cfl_diag(:,:,:)
    LOGICAL, INTENT(in)  :: stop_on_violation,output
    REAL(wp), POINTER :: cfl(:,:,:)
    
    INTEGER :: jc, level, blockNo, cell_StartIndex, cell_EndIndex
    
    ALLOCATE(cfl(LBOUND(vertical_velocity,1):UBOUND(vertical_velocity,1),&
      & LBOUND(vertical_velocity,2):UBOUND(vertical_velocity,2),&
      & LBOUND(vertical_velocity,3):UBOUND(vertical_velocity,3)))
    cfl = 0.0_wp
    
    DO blockNo = cells%start_block, cells%end_block
      CALL get_index_range(cells, blockNo, cell_StartIndex, cell_EndIndex)
      DO jc = cell_StartIndex, cell_EndIndex
        DO level=1, cells%vertical_levels(jc,blockNo)
          cfl(jc,level,blockNo)=ABS(dtime*vertical_velocity(jc,level,blockNo)/thicknesses(jc,level,blockNo))
        END DO
      END DO
    END DO
    
    IF (output) THEN
      IF (ASSOCIATED(cfl_diag)) THEN
        cfl_diag(:,:,:) = cfl(:,:,:)
      ELSE
        CALL finish('check_cfl_vertical','cfl_diag pointer for output NOT ASSOCIATED')
      ENDIF
    ENDIF
    
    CALL dbg_print('check vert.  CFL',cfl ,str_module,3,in_subset=cells)
    
    CALL check_cfl_threshold(MAXVAL(cfl),threshold,'vert',stop_on_violation)
    
    DEALLOCATE(cfl)
  END SUBROUTINE check_cfl_vertical
  
  SUBROUTINE check_cfl_threshold(maxcfl,threshold,orientation, stop_on_violation)
    REAL(wp),INTENT(in)          :: maxcfl, threshold
    CHARACTER(LEN=4), INTENT(in) :: orientation
    LOGICAL, INTENT(in)          :: stop_on_violation
    
    IF (threshold < maxcfl) THEN
      IF (stop_on_violation) THEN
        ! location lookup
        ! location print
        ! throw error
        CALL finish('check_cfl','Found violation of CFL ('//TRIM(orientation)//') criterion')
      ELSE
        CALL message('check_cfl','Found violation of CFL ('//TRIM(orientation)//') criterion')
      END IF
    END IF
  END SUBROUTINE check_cfl_threshold
  
END MODULE mo_oce_math_operators
