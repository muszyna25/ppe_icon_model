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
  USE mo_impl_constants,     ONLY: boundary, sea_boundary !,sea,land, land_boundary, sea, max_char_length, &
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce, &
    & select_solver, select_restart_mixedprecision_gmres
  
  USE mo_dynamics_config,    ONLY: nold
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_div, timer_grad
  USE mo_oce_types,          ONLY: t_hydro_ocean_state, t_solvercoeff_singleprecision, &
    & t_verticaladvection_ppm_coefficients
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, vector_product !, gc2cc
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array
  USE mo_grid_config,         ONLY: n_dom
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=12)           :: str_module    = 'oceMathOps  '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  
  PUBLIC :: grad_fd_norm_oce_3D
  PUBLIC :: grad_fd_norm_oce_3D_onblock
  PUBLIC :: div_oce_3D, div_oce_2D_sp
  PUBLIC :: rot_vertex_ocean_3D
  PUBLIC :: grad_fd_norm_oce_2D_3D, grad_fd_norm_oce_2D_3D_sp
  PUBLIC :: grad_fd_norm_oce_2D_onblock
  PUBLIC :: calculate_thickness
  PUBLIC :: map_edges2vert_3D
  PUBLIC :: check_cfl_horizontal, check_cfl_vertical
  
  
  INTERFACE div_oce_3D
    MODULE PROCEDURE div_oce_3D_mlevels
    MODULE PROCEDURE div_oce_3D_1level
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
    INTEGER :: jv, jk, blockno,jev
    INTEGER :: ile, ibe
    INTEGER :: start_index_v, end_index_v
    TYPE(t_subset_range), POINTER :: verts_in_domain
    !-----------------------------------------------------------------------
    verts_in_domain => patch_2D%verts%in_domain
    
    !i_v_ctr(:,:,:) = 0
    start_level         = 1
    end_level         = n_zlev
    
    DO blockno = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, blockno, start_index_v, end_index_v)
      DO jk = start_level, end_level
        DO jv = start_index_v, end_index_v
          
          vn_dual(jv,jk,blockno)%x = 0.0_wp
          DO jev = 1, patch_2D%verts%num_edges(jv,blockno)
            
            ! get line and block indices of edge jev around vertex jv
            ile = patch_2D%verts%edge_idx(jv,blockno,jev)
            ibe = patch_2D%verts%edge_blk(jv,blockno,jev)
            
            vn_dual(jv,jk,blockno)%x = vn_dual(jv,jk,blockno)%x        &
              & +edge2vert_coeff_cc(jv,jk,blockno,jev)%x &
              & *vn(ile,jk,ibe)
            
          END DO
        END DO ! jv = start_index_v, end_index_v
      END DO ! jk = start_level, end_level
    END DO ! blockNo = verts_in_domain%start_block, verts_in_domain%end_block
    
    ! sync the result
    !     CALL sync_patch_array(SYNC_V, patch_2D, vn_dual(:,:,:)%x(1))
    !     CALL sync_patch_array(SYNC_V, patch_2D, vn_dual(:,:,:)%x(2))
    !     CALL sync_patch_array(SYNC_V, patch_2D, vn_dual(:,:,:)%x(3))
    
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
    INTEGER :: start_edge_index, end_edge_index, blockno
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => patch_3D%p_patch_2D(1)%edges%in_domain
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockno = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockno, start_edge_index, end_edge_index)
      CALL grad_fd_norm_oce_3D_onblock( psi_c, patch_3D, &
        & grad_coeff(:,:,blockno),      &
        & grad_norm_psi_e(:,:,blockno), &
        & start_edge_index, end_edge_index, blockno)
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
    & start_edge_index, end_edge_index, blockno)
    
    TYPE(t_patch_3D ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)                   :: grad_coeff(:,:)!(nproma,n_zlev)
    REAL(wp), INTENT(in)                   :: psi_c          (nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)                :: grad_norm_psi_e(nproma,n_zlev)
    INTEGER, INTENT(in)                    :: start_edge_index, end_edge_index, blockno
    
    INTEGER :: je, jk
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
    !-----------------------------------------------------------------------
    
    iidx => patch_3D%p_patch_2D(1)%edges%cell_idx
    iblk => patch_3D%p_patch_2D(1)%edges%cell_blk
    
    DO je = start_edge_index, end_edge_index
      DO jk = 1, patch_3D%p_patch_1d(1)%dolic_e(je,blockno)
        grad_norm_psi_e(je,jk) =                                        &
          & grad_coeff(je,jk) *                                         &
          & ( psi_c(iidx(je,blockno,2),jk,iblk(je,blockno,2)) -       &
          & psi_c(iidx(je,blockno,1),jk,iblk(je,blockno,1)) )
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
  SUBROUTINE div_oce_3D_mlevels_ontriangles( vec_e, patch_3D, div_coeff, div_vec_c, opt_start_level, opt_end_level, &
    & subset_range)
    
    TYPE(t_patch_3D ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)          :: vec_e(:,:,:) ! dim: (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)          :: div_coeff(:,:,:,:)
    REAL(wp), INTENT(inout)       :: div_vec_c(:,:,:) ! dim: (nproma,n_zlev,alloc_cell_blocks)
    INTEGER, INTENT(in), OPTIONAL :: opt_start_level       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_end_level       ! optional vertical end level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    
    INTEGER :: start_level, end_level     ! vertical start and end level
    INTEGER :: jc, jk, blockno
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
    
!ICON_OMP_PARALLEL PRIVATE(iidx, iblk)
    iidx => patch_3D%p_patch_2D(1)%cells%edge_idx
    iblk => patch_3D%p_patch_2D(1)%cells%edge_blk
    
!ICON_OMP_DO PRIVATE(start_index,end_index, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockno = cells_subset%start_block, cells_subset%end_block
      CALL get_index_range(cells_subset, blockno, start_index, end_index)
      div_vec_c(:,:,blockno) = 0.0_wp
      DO jc = start_index, end_index
        DO jk = start_level, MIN(end_level, patch_3D%p_patch_1d(1)%dolic_c(jc, blockno))
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
          div_vec_c(jc,jk,blockno) =  &
            & vec_e(iidx(jc,blockno,1),jk,iblk(jc,blockno,1)) * div_coeff(jc,jk,blockno,1) + &
            & vec_e(iidx(jc,blockno,2),jk,iblk(jc,blockno,2)) * div_coeff(jc,jk,blockno,2) + &
            & vec_e(iidx(jc,blockno,3),jk,iblk(jc,blockno,3)) * div_coeff(jc,jk,blockno,3)
        END DO
      END DO
    END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE div_oce_3D_mlevels_ontriangles
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
    INTEGER :: jc, jk, blockno, max_connectivity, edgeofcell
    INTEGER ::start_index, end_index
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: cells_subset
    !-----------------------------------------------------------------------
    IF ( patch_3D%p_patch_2D(1)%cells%max_connectivity == 3) THEN
      CALL div_oce_3D_mlevels_ontriangles(vec_e, patch_3D, div_coeff, div_vec_c, &
        & opt_start_level, opt_end_level, subset_range)
      RETURN
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
    
!ICON_OMP_DO PRIVATE(start_index,end_index, jc, jk, edgeOfCell) ICON_OMP_DEFAULT_SCHEDULE
    DO blockno = cells_subset%start_block, cells_subset%end_block
      CALL get_index_range(cells_subset, blockno, start_index, end_index)
      div_vec_c(:,:,blockno) = 0.0_wp
      DO jc = start_index, end_index
        DO jk = start_level, MIN(end_level, patch_3D%p_patch_1d(1)%dolic_c(jc, blockno))
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
          div_vec_c(jc,jk,blockno) = 0.0_wp
          
          DO edgeofcell = 1, max_connectivity
            div_vec_c(jc,jk,blockno) = div_vec_c(jc,jk,blockno) + &
              & vec_e(iidx(jc,blockno,edgeofcell),jk,iblk(jc,blockno,edgeofcell)) * &
              & div_coeff(jc,jk,blockno,edgeofcell)
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
  !<Optimize:inUse>
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
    
    INTEGER :: jc, blockno
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
    DO blockno = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockno, start_index, end_index)
      DO jc = start_index, end_index
        
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
        div_vec_c(jc,blockno) =  &
          & vec_e(iidx(jc,blockno,1),iblk(jc,blockno,1)) * div_coeff(jc,level,blockno,1) + &
          & vec_e(iidx(jc,blockno,2),iblk(jc,blockno,2)) * div_coeff(jc,level,blockno,2) + &
          & vec_e(iidx(jc,blockno,3),iblk(jc,blockno,3)) * div_coeff(jc,level,blockno,3)
      END DO
    END DO
    ! !$OMP END DO
    
    ! !$OMP END PARALLEL
    IF (ltimer) CALL timer_stop(timer_div)
  END SUBROUTINE div_oce_3D_1level
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
    
    INTEGER :: jc, blockno
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
    DO blockno = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockno, start_index, end_index)
      DO jc = start_index, end_index
        
        div_vec_c(jc,blockno) =  &
          & vec_e(iidx(jc,blockno,1),iblk(jc,blockno,1)) * div_coeff(jc,blockno,1) + &
          & vec_e(iidx(jc,blockno,2),iblk(jc,blockno,2)) * div_coeff(jc,blockno,2) + &
          & vec_e(iidx(jc,blockno,3),iblk(jc,blockno,3)) * div_coeff(jc,blockno,3)
      END DO
    END DO
    ! !$OMP END DO
    
    ! !$OMP END PARALLEL
    IF (ltimer) CALL timer_stop(timer_div)
  END SUBROUTINE div_oce_2D_sp
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE grad_fd_norm_oce_2D_3D( psi_c, patch_2D, grad_coeff, grad_norm_psi_e)
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    REAL(wp), INTENT(in)    :: psi_c(:,:)             ! dim: (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(in)    :: grad_coeff(:,:)
    REAL(wp), INTENT(inout) ::  grad_norm_psi_e(:,:)  ! dim: (nproma,nblks_e)
    
    INTEGER :: blockno
    INTEGER :: start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => patch_2D%edges%in_domain
    
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_edge_index,end_edge_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockno = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockno, start_edge_index, end_edge_index)
      
      CALL grad_fd_norm_oce_2D_onblock(psi_c,  patch_2D, grad_coeff(:,blockno), grad_norm_psi_e(:,blockno), &
        & start_edge_index, end_edge_index, blockno)
      
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
  SUBROUTINE grad_fd_norm_oce_2D_onblock(psi_c,  patch_2D, grad_coeff, grad_norm_psi_e, start_index, end_index, blockno)
    !
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    REAL(wp), INTENT(in)    ::  psi_c(:,:)               ! dim: (nproma,alloc_cell_blocks)
    REAL(wp), INTENT(in)    ::  grad_coeff(:)
    REAL(wp), INTENT(inout) ::  grad_norm_psi_e(:)   ! dim: (nproma)
    INTEGER, INTENT(in)     :: start_index, end_index, blockno
    
    INTEGER :: je
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    !-----------------------------------------------------------------------
    
    iidx => patch_2D%edges%cell_idx
    iblk => patch_2D%edges%cell_blk
    
    DO je = start_index, end_index
      ! compute the normal derivative
      ! by the finite difference approximation
      ! (see Bonaventura and Ringler MWR 2005)
!       IF (iidx(je,blockno,1) < 1 .or. iidx(je,blockno,2) < 1) THEN
!         WRITE(0,*) "je=", je, " blockno=", blockno, &
!           & " iidx(je,blockno,1)=", iidx(je,blockno,1), " iidx(je,blockno,2)=", iidx(je,blockno,2)
!         CALL finish("invalid connectivity", "")
!       ENDIF
      
      grad_norm_psi_e(je) =  &
        & (psi_c(iidx(je,blockno,2),iblk(je,blockno,2))-psi_c(iidx(je,blockno,1),iblk(je,blockno,1)))&
        & * grad_coeff(je)
      
    END DO
    
  END SUBROUTINE grad_fd_norm_oce_2D_onblock
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  ! the same as grad_fd_norm_oce_2D_3D_sp in single precisison
  SUBROUTINE grad_fd_norm_oce_2D_3D_sp( psi_c, patch_2D, grad_coeff, grad_norm_psi_e)
    TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
    REAL(sp), INTENT(in)    :: psi_c(:,:)             ! dim: (nproma,alloc_cell_blocks)
    REAL(sp), INTENT(in)    :: grad_coeff(:,:)
    REAL(sp), INTENT(inout) ::  grad_norm_psi_e(:,:)  ! dim: (nproma,nblks_e)
    
    INTEGER :: je, blockno
    INTEGER :: start_edge_index, end_edge_index
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => patch_2D%edges%in_domain
    
    iidx => patch_2D%edges%cell_idx
    iblk => patch_2D%edges%cell_blk
    
!ICON_OMP_PARALLEL_DO PRIVATE(blockNo,start_edge_index,end_edge_index,je) ICON_OMP_DEFAULT_SCHEDULE
    DO blockno = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockno, start_edge_index, end_edge_index)
      
      DO je = start_edge_index, end_edge_index
        ! compute the normal derivative
        ! by the finite difference approximation
        ! (see Bonaventura and Ringler MWR 2005)
        grad_norm_psi_e(je,blockno) =  &
          & (psi_c(iidx(je,blockno,2),iblk(je,blockno,2))-psi_c(iidx(je,blockno,1),iblk(je,blockno,1)))&
          & * grad_coeff(je,blockno)
        
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
  !! does the tangential velocity calculate from velocity vector at vertices (vn_dual), while RBF uses
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
    REAL(wp) :: z_vort_int
    REAL(wp) :: z_vort_boundary(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: start_level, end_level
    INTEGER :: jv, jk, blockno, jev
    INTEGER :: ile, ibe
    INTEGER :: il_v1, il_v2,ib_v1, ib_v2
    INTEGER :: start_index_v, end_index_v
    
    REAL(wp) :: z_vt(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    INTEGER, POINTER :: ibnd_edge_idx(:,:,:,:), ibnd_edge_blk(:,:,:,:), i_edge_idx(:,:,:,:)
    !REAL(wp), POINTER :: z_orientation(:,:,:,:)
    
    TYPE(t_subset_range), POINTER :: verts_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D           => patch_3D%p_patch_2D(1)
    verts_in_domain => patch_2D%verts%in_domain
    start_level         = 1
    end_level         = n_zlev
    
    
    z_vort_boundary  = 0.0_wp
    rot_vec_v(:,:,:) = 0.0_wp
    z_vt(:,:,:)      = 0.0_wp
    
    !set pointer that carry edge information
    ibnd_edge_idx    => p_op_coeff%bnd_edge_idx
    ibnd_edge_blk    => p_op_coeff%bnd_edge_blk
    i_edge_idx       => p_op_coeff%edge_idx
    !z_orientation    => p_op_coeff%orientation
    
    
    !In this loop tangential velocity is calulated at boundaries
    DO blockno = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, blockno, start_index_v, end_index_v)
      DO jk = start_level, end_level
        !!$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(PRIVATE)  &
        !!$OMP   SHARED(u_vec_e,v_vec_e,patch_2D,rot_vec_v,blockNo) FIRSTPRIVATE(jk)
        DO jv = start_index_v, end_index_v
          
          DO jev = 1, patch_2D%verts%num_edges(jv,blockno)
            
            ! get line and block indices of edge jev around vertex jv
            ile = patch_2D%verts%edge_idx(jv,blockno,jev)
            ibe = patch_2D%verts%edge_blk(jv,blockno,jev)
            
            !IF ( v_base%lsm_e(ile,jk,ibe) == boundary ) THEN
            IF(patch_3D%lsm_e(ile,jk,ibe) == boundary)THEN
              !calculate tangential velocity
              il_v1 = patch_2D%edges%vertex_idx(ile,ibe,1)
              ib_v1 = patch_2D%edges%vertex_blk(ile,ibe,1)
              il_v2 = patch_2D%edges%vertex_idx(ile,ibe,2)
              ib_v2 = patch_2D%edges%vertex_blk(ile,ibe,2)
              
              z_vt(ile,jk,ibe)= &
                & - DOT_PRODUCT(vn_dual(il_v1,jk,ib_v1)%x,&
                & p_op_coeff%edge2vert_coeff_cc_t(ile,jk,ibe,1)%x)&
                & + DOT_PRODUCT(vn_dual(il_v2,jk,ib_v2)%x,&
                & p_op_coeff%edge2vert_coeff_cc_t(ile,jk,ibe,2)%x)
            ENDIF
          END DO
        END DO
        !!$OMP END PARALLEL DO
      END DO
    END DO
    
    !In this loop vorticity at vertices is calculated
    DO blockno = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, blockno, start_index_v, end_index_v)
      DO jk = start_level, end_level
        !!$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(PRIVATE)  &
        !!$OMP   SHARED(u_vec_e,v_vec_e,patch_2D,rot_vec_v,blockNo) FIRSTPRIVATE(jk)
        DO jv = start_index_v, end_index_v
          
          z_vort_int = 0.0_wp
          
          DO jev = 1, patch_2D%verts%num_edges(jv,blockno)
            
            ! get line and block indices of edge jev around vertex jv
            ile = patch_2D%verts%edge_idx(jv,blockno,jev)
            ibe = patch_2D%verts%edge_blk(jv,blockno,jev)
            !add contribution of normal velocity at edge (ile,ibe) to rotation
            !IF ( v_base%lsm_e(ile,jk,ibe) == sea) THEN
            ! sea, sea_boundary, boundary (edges only), land_boundary, land =
            !  -2,      -1,         0,                  1,             2
            !Distinction between sea-lean-boundary is taken into account by coeffcients.
            !It is assumed here that vn is already zero at boundary edges.
            
            z_vort_int = z_vort_int + vn(ile,jk,ibe)*p_op_coeff%rot_coeff(jv,jk,blockno,jev)
            
          END DO
          
          !Finalize vorticity calculation by closing the dual loop along boundary edges
          !IF(i_v_bnd_edge_ctr(jv,jk,blockNo)==2)THEN
          IF(p_op_coeff%bnd_edges_per_vertex(jv,jk,blockno)==2)THEN
            
            z_vort_boundary(jv,jk,blockno)=&
              & z_vt(ibnd_edge_idx(jv,jk,blockno,1),jk,ibnd_edge_blk(jv,jk,blockno,1))&
              & *p_op_coeff%rot_coeff(jv,jk,blockno,i_edge_idx(jv,jk,blockno,1))&
              & +&
              & z_vt(ibnd_edge_idx(jv,jk,blockno,2),jk,ibnd_edge_blk(jv,jk,blockno,2))&
              & *p_op_coeff%rot_coeff(jv,jk,blockno,i_edge_idx(jv,jk,blockno,2))
            
            !ELSEIF(i_v_bnd_edge_ctr(jv,jk,blockNo)==4)THEN
          ELSEIF(p_op_coeff%bnd_edges_per_vertex(jv,jk,blockno)==4)THEN
            !In case of 4 boundary edges within a dual loop, we have 2 land triangles
            !around the vertex. these two land triangles have one vertex in common and are
            !seperated by two wet triangles
            
            z_vort_boundary(jv,jk,blockno)=&
              & z_vt(ibnd_edge_idx(jv,jk,blockno,1),jk,ibnd_edge_blk(jv,jk,blockno,1))&
              & *p_op_coeff%rot_coeff(jv,jk,blockno,i_edge_idx(jv,jk,blockno,1))&
              & +&
              & z_vt(ibnd_edge_idx(jv,jk,blockno,2),jk,ibnd_edge_blk(jv,jk,blockno,2))&
              & *p_op_coeff%rot_coeff(jv,jk,blockno,i_edge_idx(jv,jk,blockno,2))&
              & +&
              & z_vt(ibnd_edge_idx(jv,jk,blockno,3),jk,ibnd_edge_blk(jv,jk,blockno,3))&
              & *p_op_coeff%rot_coeff(jv,jk,blockno,i_edge_idx(jv,jk,blockno,3))&
              & +&
              & z_vt(ibnd_edge_idx(jv,jk,blockno,4),jk,ibnd_edge_blk(jv,jk,blockno,4))&
              & *p_op_coeff%rot_coeff(jv,jk,blockno,i_edge_idx(jv,jk,blockno,4))
          ENDIF
          
          !Final vorticity calculation
          !TODO ram
          !       rot_vec_v(jv,jk,blockNo) = (z_vort_int + z_vort_boundary(jv,jk,blockNo)) / &
          !         & patch_2D%verts%dual_area(jv,blockNo)
          rot_vec_v(jv,jk,blockno) = z_vort_int + z_vort_boundary(jv,jk,blockno)
          
        END DO
        !!$OMP END PARALLEL DO
      END DO
    END DO
    
    ! DO blockNo = i_startblk_v, i_endblk_v
    !   CALL get_indices_v(patch_2D, blockNo, i_startblk_v, i_endblk_v, &
    !                      start_index_v, end_index_v, rl_start_v, rl_end_v)
    !   DO jk = start_level, end_level
    !     DO jv = start_index_v, end_index_v
    ! IF(rot_vec_v(jv,jk,blockNo)/=0.0_wp)THEN
    ! write(123456,*)'rot 3D:',jk,jv,blockNo,rot_vec_v(jv,jk,blockNo),&
    ! &z_vort_boundary(jv,jk,blockNo), p_op_coeff%bnd_edges_per_vertex(jv,jk,blockNo)
    ! ENDIF
    !     END DO
    !   END DO
    ! END DO
  END SUBROUTINE rot_vertex_ocean_3D
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
    !REAL(wp), INTENT(IN)                      :: ice_hi(nproma,1,p_patch_3D%p_patch_2D(1)%nblks_c)
    TYPE(t_operator_coeff), INTENT(in)      :: operators_coefficients
    TYPE(t_solvercoeff_singleprecision), INTENT(inout) :: solvercoeff_sp
    REAL(wp), OPTIONAL :: inTopCellThickness(:,:)
    
    !  local variables
    INTEGER :: cell_StartIndex, cell_EndIndex
    INTEGER :: edge_StartIndex, edge_EndIndex
    INTEGER :: jc, blockno, je, jk
    INTEGER :: thislevel, levelabove, levelbelow, level2below, cell_levels
    
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    REAL(wp)           :: z_dist_e_c1, z_dist_e_c2
    TYPE(t_subset_range), POINTER :: all_cells, all_edges, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    
    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    REAL(wp) :: top_vn_1, top_vn_2, integrated_vn
    REAL(wp), POINTER :: top_coeffs(:,:,:), integrated_coeffs(:,:,:), sum_to_2D_coeffs(:,:,:)
    REAL(wp), POINTER :: cell_thickeness(:,:,:), edge_thickeness(:,:,:)
    REAL(wp), POINTER :: inv_cell_thickeness(:,:,:), inv_edge_thickeness(:,:,:)
    REAL(wp), POINTER :: inv_prisms_center_distance(:,:,:), inv_edgefaces_middle_distance(:,:,:)
    REAL(wp)  :: cell_thickeness_1, cell_thickeness_2
    !-------------------------------------------------------------------------------
    ! pointers for the ppm vertical transport
    TYPE(t_verticaladvection_ppm_coefficients), POINTER :: vertadvppm
    !-------------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    patch_2D            => patch_3D%p_patch_2D(1)
    all_cells           => patch_2D%cells%ALL
    all_edges           => patch_2D%edges%ALL
    edges_in_domain     => patch_2D%edges%in_domain
    cell_thickeness     => patch_3D%p_patch_1d(1)%prism_thick_c
    inv_cell_thickeness => patch_3D%p_patch_1d(1)%inv_prism_thick_c
    inv_prisms_center_distance => patch_3D%p_patch_1d(1)%inv_prism_center_dist_c
    edge_thickeness     => patch_3D%p_patch_1d(1)%prism_thick_e
    inv_edge_thickeness => patch_3D%p_patch_1d(1)%inv_prism_thick_e
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
      DO blockno = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockno, cell_StartIndex, cell_EndIndex)
        DO jc = cell_StartIndex, cell_EndIndex
          IF ( patch_3D%p_patch_1d(1)%dolic_c(jc,blockno) > 0 ) THEN
            cell_thickeness(jc,1,blockno) = inTopCellThickness(jc,blockno)
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO

    ELSE
!ICON_OMP_DO PRIVATE(cell_StartIndex, cell_EndIndex, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO blockno = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockno, cell_StartIndex, cell_EndIndex)
        DO jc = cell_StartIndex, cell_EndIndex
          IF ( patch_3D%p_patch_1d(1)%dolic_c(jc,blockno) > 0 ) THEN
            cell_thickeness(jc,1,blockno) = &
              & patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(jc,1,blockno) + &
              & ocean_state%p_prog(nold(1))%h(jc,blockno)
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO

    ENDIF
    
!ICON_OMP_DO PRIVATE(cell_StartIndex, cell_EndIndex, jc, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockno = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockno, cell_StartIndex, cell_EndIndex)
      DO jc = cell_StartIndex, cell_EndIndex
        IF ( patch_3D%p_patch_1d(1)%dolic_c(jc,blockno) > 0 ) THEN
          
          patch_3D%p_patch_1d(1)%prism_center_dist_c(jc,2,blockno) = 0.5_wp * &
            & (cell_thickeness(jc,1,blockno) + cell_thickeness(jc,2,blockno))
          
          patch_3D%p_patch_1d(1)%prism_volume(jc,1,blockno) = cell_thickeness(jc,1,blockno) * &
            & patch_2D%cells%area(jc,blockno)
          
          inv_cell_thickeness(jc,1,blockno) = 1.0_wp / cell_thickeness(jc,1,blockno)
          
          inv_prisms_center_distance(jc,2,blockno) = &
            & 1.0_wp / patch_3D%p_patch_1d(1)%prism_center_dist_c(jc,2,blockno)
          
          ocean_state%p_diag%thick_c(jc,blockno) = ocean_state%p_prog(nold(1))%h(jc,blockno) + patch_3D%column_thick_c(jc,blockno)
          
          patch_3D%p_patch_1d(1)%depth_cellmiddle(jc,1,blockno) = cell_thickeness(jc,1,blockno) * 0.5_wp
          patch_3D%p_patch_1d(1)%depth_cellinterface(jc,2,blockno) = cell_thickeness(jc,1,blockno)
          
        ENDIF
        
        DO jk=2, patch_3D%p_patch_1d(1)%dolic_c(jc,blockno)
          patch_3D%p_patch_1d(1)%depth_cellmiddle(jc,jk,blockno) = &
            & patch_3D%p_patch_1d(1)%depth_cellinterface(jc,jk,blockno) + cell_thickeness(jc,jk,blockno) * 0.5_wp
          patch_3D%p_patch_1d(1)%depth_cellinterface(jc,jk+1,blockno) = &
            & patch_3D%p_patch_1d(1)%depth_cellinterface(jc,jk,blockno) + cell_thickeness(jc,jk,blockno)
        ENDDO
        
      END DO
    END DO
!ICON_OMP_END_DO

    !----------------------------------------------------------------------------------------
    IF ( iswm_oce == 1 ) THEN  !  SWM
!ICON_OMP_DO PRIVATE(cell_StartIndex, cell_EndIndex, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO blockno = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockno, cell_StartIndex, cell_EndIndex)
        !calculate for each fluid colum the total depth, i.e.
        !from bottom boundary to surface height, i.e. using individual bathymetry for SWM
        DO jc = cell_StartIndex, cell_EndIndex
          IF(patch_3D%lsm_c(jc,1,blockno) <= sea_boundary)THEN
            
            ocean_state%p_diag%thick_c(jc,blockno) = ocean_state%p_prog(nold(1))%h(jc,blockno)&
              & - p_ext_data%oce%bathymetry_c(jc,blockno)
            !        &  - ice_hi(jc,1,blockNo)
          ELSE
            ocean_state%p_diag%thick_c(jc,blockno) = 0.0_wp
          ENDIF
        END DO
      END DO!write(*,*)'bathymetry',maxval(p_ext_data%oce%bathymetry_c),minval(p_ext_data%oce%bathymetry_c)
!ICON_OMP_END_DO
      
      !Step 2: calculate edge-located variables for 2D and 3D case from respective cell variables
      !For SWE : thick_e = thickness of fluid column at edges
      !         h_e     = surface elevation at edges, without depth of first layer
!ICON_OMP_DO PRIVATE(edge_StartIndex, edge_EndIndex, je, il_c1, ib_c1, il_c2, ib_c2, &
!ICON_OMP z_dist_e_c1, z_dist_e_c2) ICON_OMP_DEFAULT_SCHEDULE
      DO blockno = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockno, edge_StartIndex, edge_EndIndex)
        DO je = edge_StartIndex, edge_EndIndex
          
          il_c1 = patch_2D%edges%cell_idx(je,blockno,1)
          ib_c1 = patch_2D%edges%cell_blk(je,blockno,1)
          il_c2 = patch_2D%edges%cell_idx(je,blockno,2)
          ib_c2 = patch_2D%edges%cell_blk(je,blockno,2)
          
          z_dist_e_c1 = 0.5_wp!z_dist_e_c1=p_patch%edges%edge_cell_length(je,blockNo,1)
          z_dist_e_c2 = 0.5_wp!z_dist_e_c2=p_patch%edges%edge_cell_length(je,blockNo,2)
          
          IF(patch_3D%lsm_e(je,1,blockno) <= sea_boundary)THEN
            ocean_state%p_diag%thick_e(je,blockno) = ( z_dist_e_c1*ocean_state%p_diag%thick_c(il_c1,ib_c1)&
              & +   z_dist_e_c2*ocean_state%p_diag%thick_c(il_c2,ib_c2) )&
              & /(z_dist_e_c1+z_dist_e_c2)
            
            ocean_state%p_diag%h_e(je,blockno) = ( z_dist_e_c1*ocean_state%p_prog(nold(1))%h(il_c1,ib_c1)&
              & +   z_dist_e_c2*ocean_state%p_prog(nold(1))%h(il_c2,ib_c2) )&
              & /(z_dist_e_c1+z_dist_e_c2)
          ELSE
            ocean_state%p_diag%h_e(je,blockno) = 0.0_wp
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO

!ICON_OMP_MASTER
      CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%thick_e)
      CALL sync_patch_array(sync_e, patch_2D, ocean_state%p_diag%h_e)
!ICON_OMP_END_MASTER
!ICON_OMP_BARRIER

      !2) Thickness at edges
!ICON_OMP_DO PRIVATE(edge_StartIndex, edge_EndIndex, je) ICON_OMP_DEFAULT_SCHEDULE
      DO blockno = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, blockno, edge_StartIndex, edge_EndIndex)
        DO je = edge_StartIndex, edge_EndIndex
          IF ( patch_3D%p_patch_1d(1)%dolic_e(je,blockno) > 0 ) THEN
            
            edge_thickeness(je,1,blockno)&
              & = patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_e(je,1,blockno) + ocean_state%p_diag%h_e(je,blockno)
            
            inv_edge_thickeness(je,1,blockno)= 1.0_wp / edge_thickeness(je,1,blockno)
            
            inv_edgefaces_middle_distance(je,2,blockno) = 2.0_wp / &
              & (edge_thickeness(je,1,blockno) + edge_thickeness(je,2,blockno))
            
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
      DO blockno = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, blockno, edge_StartIndex, edge_EndIndex)
        DO je = edge_StartIndex, edge_EndIndex
          
          il_c1 = patch_2D%edges%cell_idx(je,blockno,1)
          ib_c1 = patch_2D%edges%cell_blk(je,blockno,1)
          il_c2 = patch_2D%edges%cell_idx(je,blockno,2)
          ib_c2 = patch_2D%edges%cell_blk(je,blockno,2)
          
          z_dist_e_c1 = 0.5_wp!z_dist_e_c1=p_patch%edges%edge_cell_length(je,blockNo,1)
          z_dist_e_c2 = 0.5_wp!z_dist_e_c2=p_patch%edges%edge_cell_length(je,blockNo,2)
          
          IF ( patch_3D%p_patch_1d(1)%dolic_e(je,blockno) > 0 ) THEN
            
            ocean_state%p_diag%h_e(je,blockno) = ( z_dist_e_c1*ocean_state%p_prog(nold(1))%h(il_c1,ib_c1)&
              & +   z_dist_e_c2*ocean_state%p_prog(nold(1))%h(il_c2,ib_c2) )&
              & /(z_dist_e_c1+z_dist_e_c2)
            
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
      DO blockno = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, blockno, edge_StartIndex, edge_EndIndex)
        DO je = edge_StartIndex, edge_EndIndex
          IF ( patch_3D%p_patch_1d(1)%dolic_e(je,blockno) > 0 ) THEN
            
            edge_thickeness(je,1,blockno)&
              & = patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_e(je,1,blockno) + ocean_state%p_diag%h_e(je,blockno)
            
            inv_edge_thickeness(je,1,blockno)= 1.0_wp / edge_thickeness(je,1,blockno)
            
            inv_edgefaces_middle_distance(je,2,blockno) = 2.0_wp / &
              & (edge_thickeness(je,1,blockno) + edge_thickeness(je,2,blockno))
            
            ocean_state%p_diag%thick_e(je,blockno) = ocean_state%p_diag%h_e(je,blockno) &
              & + patch_3D%column_thick_e(je,blockno)
            
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
!ICON_OMP_DO PRIVATE(edge_StartIndex, edge_EndIndex, je, cell_1_index, cell_1_block, &
!ICON_OMP cell_2_index, cell_2_block, cell_thickeness_1, cell_thickeness_2) ICON_OMP_DEFAULT_SCHEDULE
    DO blockno = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockno, edge_StartIndex, edge_EndIndex)
      DO je = edge_StartIndex, edge_EndIndex
        
        IF ( patch_3D%p_patch_1d(1)%dolic_e(je,blockno) > 0 ) THEN
          
          ! get the two cells of the edge
          cell_1_index = patch_2D%edges%cell_idx(je,blockno,1)
          cell_1_block = patch_2D%edges%cell_blk(je,blockno,1)
          cell_2_index = patch_2D%edges%cell_idx(je,blockno,2)
          cell_2_block = patch_2D%edges%cell_blk(je,blockno,2)
          cell_thickeness_1 = cell_thickeness(cell_1_index, 1, cell_1_block)
          cell_thickeness_2 = cell_thickeness(cell_2_index, 1, cell_2_block)
          
          sum_to_2D_coeffs(1, je, blockno) = top_coeffs(1, je, blockno) * cell_thickeness_1 + integrated_coeffs(1, je, blockno)
          sum_to_2D_coeffs(2, je, blockno) = top_coeffs(2, je, blockno) * cell_thickeness_1 + integrated_coeffs(2, je, blockno)
          sum_to_2D_coeffs(3, je, blockno) = top_coeffs(3, je, blockno) * cell_thickeness_1 + integrated_coeffs(3, je, blockno)
          sum_to_2D_coeffs(4, je, blockno) = top_coeffs(4, je, blockno) * cell_thickeness_2 + integrated_coeffs(4, je, blockno)
          sum_to_2D_coeffs(5, je, blockno) = top_coeffs(5, je, blockno) * cell_thickeness_2 + integrated_coeffs(5, je, blockno)
          sum_to_2D_coeffs(6, je, blockno) = top_coeffs(6, je, blockno) * cell_thickeness_2 + integrated_coeffs(6, je, blockno)
          
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
    DO blockno = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockno, cell_StartIndex, cell_EndIndex)
      vertadvppm => operators_coefficients%verticaladvectionppmcoeffs(blockno)
      DO jc = cell_StartIndex, cell_EndIndex
        
        cell_levels = patch_3D%p_patch_1d(1)%dolic_c(jc,blockno)
        
        thislevel  = 1
        levelbelow = 2
        IF ( cell_levels >= levelbelow ) THEN
          
          vertadvppm%cellheightratio_this_tobelow(jc, thislevel) = &
            & cell_thickeness(jc, thislevel, blockno) / cell_thickeness(jc, levelbelow, blockno)
          
          vertadvppm%cellheightratio_this_tothisbelow(jc, thislevel) = &
            & cell_thickeness(jc, thislevel, blockno) / &
            & (cell_thickeness(jc, thislevel, blockno) + cell_thickeness(jc, levelbelow, blockno))
          
          vertadvppm%cellheight_2xbelow_x_ratiothis_tothisbelow(jc,thislevel) = &
            & 2._wp * cell_thickeness(jc,levelbelow, blockno) * &
            & vertadvppm%cellheightratio_this_tothisbelow(jc, thislevel)
          
        ENDIF
        
        thislevel  = 2
        levelabove = 1
        levelbelow = 3
        level2below = 4
        IF ( cell_levels >= levelbelow ) THEN
          
          
          vertadvppm%cellheightratio_this_tothisabovebelow(jc,thislevel) = &
            & cell_thickeness(jc, thislevel ,blockno) / &
            & (cell_thickeness(jc,levelabove,blockno) + cell_thickeness(jc,thislevel,blockno)    &
            & + cell_thickeness(jc,levelbelow,blockno))
          
          vertadvppm%cellheightratio_2xaboveplusthis_tothisbelow(jc,thislevel) = &
            & (2._wp * cell_thickeness(jc,levelabove,blockno) + cell_thickeness(jc,thislevel,blockno))     &
            & / (cell_thickeness(jc,levelbelow,blockno) + cell_thickeness(jc,thislevel,blockno))
          
          vertadvppm%cellheightratio_2xbelowplusthis_tothisabove(jc,thislevel) = &
            & + (cell_thickeness(jc,thislevel,blockno) + 2._wp * cell_thickeness(jc,levelbelow,blockno))   &
            & / (cell_thickeness(jc,levelabove,blockno) + cell_thickeness(jc,thislevel,blockno))
          
          vertadvppm%cellheightratio_thisabove_to2xthisplusbelow(jc,thislevel) =                         &
            & (cell_thickeness(jc,levelabove,blockno) + cell_thickeness(jc,thislevel,blockno))            &
            & / (2._wp*cell_thickeness(jc,thislevel,blockno) + cell_thickeness(jc,levelbelow,blockno))
          
          vertadvppm%cellheightratio_thisbelow_to2xthisplusabove(jc,thislevel) =                 &
            & (cell_thickeness(jc,levelbelow,blockno) + cell_thickeness(jc,thislevel,blockno))                  &
            & / (2._wp*cell_thickeness(jc,thislevel,blockno) + cell_thickeness(jc,levelabove,blockno))
          ! = 1 / cellHeightRatio_2xBelowplusThis_toThisAbove(levelBelow)
          
        ENDIF
        
        IF ( cell_levels >= level2below ) THEN
          vertadvppm%cellheight_inv_thisabovebelow2below(jc,thislevel) =                                  &
            & 1._wp / (cell_thickeness(jc,levelabove,blockno) + cell_thickeness(jc,thislevel,blockno)       &
            & + cell_thickeness(jc,levelbelow,blockno) + cell_thickeness(jc,level2below,blockno))
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
  !<Optimize:inUse>
  SUBROUTINE check_cfl_horizontal(normal_velocity,inv_dual_edge_length,timestep,edges,threshold, &
    & cfl_diag, stop_on_violation, output)
    REAL(wp),POINTER :: normal_velocity(:,:,:)
    TYPE(t_subset_range)  :: edges
    REAL(wp), INTENT(in)  :: inv_dual_edge_length(:,:), threshold, timestep
    REAL(wp), POINTER :: cfl_diag(:,:,:)
    LOGICAL, INTENT(in)   :: stop_on_violation, output
    REAL(wp), POINTER :: cfl(:,:,:)
    
    INTEGER :: je,jk,blockno,start_index,end_index
    
    ALLOCATE(cfl(LBOUND(normal_velocity,1):UBOUND(normal_velocity,1),&
      & LBOUND(normal_velocity,2):UBOUND(normal_velocity,2),&
      & LBOUND(normal_velocity,3):UBOUND(normal_velocity,3)))
    cfl = 0.0_wp
    
    DO blockno = edges%start_block, edges%end_block
      CALL get_index_range(edges,blockno,start_index,end_index)
      DO je = start_index,end_index
        DO jk = 1, n_zlev
          cfl(je,jk,blockno) = ABS(dtime*normal_velocity(je,jk,blockno)*inv_dual_edge_length(je,blockno))
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
  
  !<Optimize:inUse>
  SUBROUTINE check_cfl_vertical(vertical_velocity, thicknesses, timestep, cells, threshold, &
    & cfl_diag, stop_on_violation, output)
    REAL(wp),POINTER :: vertical_velocity(:,:,:), thicknesses(:,:,:)
    REAL(wp), INTENT(in) :: timestep, threshold
    TYPE(t_subset_range) :: cells
    REAL(wp), POINTER :: cfl_diag(:,:,:)
    LOGICAL, INTENT(in)  :: stop_on_violation,output
    REAL(wp), POINTER :: cfl(:,:,:)
    
    INTEGER :: jc, jk, blockno, cell_StartIndex, cell_EndIndex
    
    ALLOCATE(cfl(LBOUND(vertical_velocity,1):UBOUND(vertical_velocity,1),&
      & LBOUND(vertical_velocity,2):UBOUND(vertical_velocity,2),&
      & LBOUND(vertical_velocity,3):UBOUND(vertical_velocity,3)))
    cfl = 0.0_wp
    
    DO blockno = cells%start_block, cells%end_block
      CALL get_index_range(cells, blockno, cell_StartIndex, cell_EndIndex)
      DO jc = cell_StartIndex, cell_EndIndex
        DO jk=1, cells%vertical_levels(jc,blockno)
          cfl(jc,jk,blockno)=ABS(dtime*vertical_velocity(jc,jk,blockno)/thicknesses(jc,jk,blockno))
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
  
  !<Optimize:inUse>
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
