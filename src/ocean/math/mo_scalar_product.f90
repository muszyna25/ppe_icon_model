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
MODULE mo_scalar_product
  !-------------------------------------------------------------------------
  !
  USE mo_kind,               ONLY: wp, sp
  USE mo_exception,          ONLY: finish
  USE mo_parallel_config,    ONLY: nproma, p_test_run
  USE mo_impl_constants,     ONLY: sea_boundary, sea
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_oce_types,          ONLY: t_hydro_ocean_diag, t_solverCoeff_singlePrecision
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce, fast_performance_level
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates,cvec2gvec!, gc2cc, vector_product
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff, no_primal_edges, no_dual_edges
  USE mo_oce_math_operators,  ONLY: rot_vertex_ocean_3d, map_edges2vert_3d
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_e, sync_v,sync_patch_array!,sync_c,  & sync_idx, global_max
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  
  IMPLICIT NONE
  
  PRIVATE
  
 
  PUBLIC :: calc_scalar_product_veloc_3d
  PUBLIC :: nonlinear_coriolis_3d
  PUBLIC :: nonlinear_coriolis_3d_old
  PUBLIC :: map_edges2vert_3d
  PUBLIC :: map_edges2cell_3d
  PUBLIC :: map_cell2edges_3D
  PUBLIC :: map_edges2edges_viacell_3d
  PUBLIC :: map_edges2edges_viavert_3D
  PUBLIC :: map_edges2edges_viacell_3D_const_z, map_edges2edges_viacell_2d_constZ_sp
  ! PUBLIC :: map_edges2edges_viacell_2d_1lev_const_z_fast0
  !PRIVATE :: map_cell2edges_2d


  INTERFACE map_edges2edges_viacell_3d 
    MODULE PROCEDURE map_edges2edges_viacell_3d_1lev
    MODULE PROCEDURE map_edges2edges_viacell_3d_mlev
  END INTERFACE

  INTERFACE map_edges2edges_viacell_3d_const_z 
    MODULE PROCEDURE map_edges2edges_viacell_3D_1lev_constZ
    MODULE PROCEDURE map_edges2edges_viacell_3d_1lev_constZs
    MODULE PROCEDURE map_edges2edges_viacell_3d_mlev_const_z
    MODULE PROCEDURE map_edges2edges_viacell_3d_mlev_constZs
  END INTERFACE



  INTERFACE map_cell2edges_3D
    MODULE PROCEDURE map_cell2edges_3D_1level
    MODULE PROCEDURE map_cell2edges_3D_mlevels
  END INTERFACE

  INTERFACE map_edges2cell_3d
    
    MODULE PROCEDURE map_edges2cell_with_height_3d
    MODULE PROCEDURE map_edges2cell_no_height_3d
    
  END INTERFACE

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
    !!  mpi parallelized by LL
!<Optimize_Used>
  SUBROUTINE calc_scalar_product_veloc_3D( patch_3D, vn_e_old, vn_e_new,&
    & p_diag, operators_coefficients)
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    REAL(wp), INTENT(in)      :: vn_e_old(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(in)      :: vn_e_new(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_hydro_ocean_diag)  :: p_diag
    TYPE(t_operator_coeff)    :: operators_coefficients
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: jc, jb, jk

    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_cells => patch_2D%cells%all
    !-----------------------------------------------------------------------
    slev = 1
    elev = n_zlev

    CALL map_edges2vert_3d(patch_3D%p_patch_2D(1), vn_e_old, operators_coefficients%edge2vert_coeff_cc, &
      & p_diag%p_vn_dual)

    !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
    CALL map_edges2cell_3d(patch_3D, vn_e_old, operators_coefficients, p_diag%p_vn)

    CALL map_cell2edges_3D( patch_3D, p_diag%p_vn, p_diag%ptp_vn, operators_coefficients)

!      CALL map_edges2edges_viacell_3D( patch_3D,    &
!                                     & vn_e_old,      &
!                                     & operators_coefficients,    &
!                                     & p_diag%ptp_vn)

   CALL sync_patch_array(SYNC_E, patch_2D, p_diag%ptp_vn)

    !--------------------------------------------------------------    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      !calculate kinetic energy
      DO jk = slev, elev
        DO jc =  i_startidx_c, i_endidx_c

          !IF ( v_base%lsm_c(jc,jk,jb) > sea_boundary ) THEN
          IF(patch_3D%lsm_c(jc,jk,jb) > sea_boundary)THEN
            p_diag%kin(jc,jk,jb) = 0.0_wp
          ELSE
            p_diag%kin(jc,jk,jb) = 0.5_wp * &
              & DOT_PRODUCT(p_diag%p_vn(jc,jk,jb)%x, p_diag%p_vn(jc,jk,jb)%x)
            !p_diag%kin(jc,jk,jb) = 0.5_wp*DOT_PRODUCT(z_pv_cc(jc,jk,jb)%x,z_pv_cc(jc,jk,jb)%x)
            !IF(p_diag%kin(jc,jk,jb)/=0.0_wp)&
            !&write(*,*)'Pv',jc,jb,p_diag%p_vn(jc,jk,jb)%x, z_pv_cc(jc,jk,jb)%x
          ENDIF
        END DO
      END DO
    END DO
    ! LL: no sync is required
    !--------------------------------------------------------------    

    !--------------------------------------------------------------    
    ! DO jk = slev, elev
    !   write(*,*)'max/min kin energy:',maxval(p_diag%kin(:,jk,:)), minval(p_diag%kin(:,jk,:))!,&
    ! !&maxval(z_kin(:,1,:)), minval(z_kin(:,1,:))
    ! END DO
    !convert cartesian velocity vector p_diag%p_vn(jc,jk,jb)%x to geographical coordinate system
    !for output
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jk = slev, elev
        DO jc =  i_startidx_c, i_endidx_c
          CALL cvec2gvec ( p_diag%p_vn(jc,jk,jb)%x(1),     &
            & p_diag%p_vn(jc,jk,jb)%x(2),     &
            & p_diag%p_vn(jc,jk,jb)%x(3),     &
            & patch_2D%cells%center(jc,jb)%lon,&
            & patch_2D%cells%center(jc,jb)%lat,&
            & p_diag%u(jc,jk,jb), p_diag%v(jc,jk,jb) )
        END DO
      END DO
    END DO

  END SUBROUTINE calc_scalar_product_veloc_3d
  !-------------------------------------------------------------------------
  
 
  !-------------------------------------------------------------------------
  !>
  !! Optimized version of the nonlinear_coriolis_3d
  !! Note: intel -o# will produce very different results when using this version,
  !!  not clear if the model exhibits great sensitivity or there's another problem
  !<Optimize_Used.done>
  SUBROUTINE linear_coriolis_3d_fast(patch_3D, vn, p_vn_dual, vort_v, &
    & operators_coefficients, vort_flux)
    TYPE(t_patch_3D ),TARGET,INTENT(IN) :: patch_3D
    REAL(wp), INTENT(INOUT)                    :: vn(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(INout)  :: p_vn_dual(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp), INTENT(INOUT)                    :: vort_v   (nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(IN)          :: operators_coefficients
    REAL(wp), INTENT(INOUT)                    :: vort_flux(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: slev! , elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk
    REAL(wp) :: this_vort_flux(n_zlev, 2) ! for each of the two vertices

    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-----------------------------------------------------------------------
    slev    = 1
    ! elev    = n_zlev

    CALL rot_vertex_ocean_3d( patch_3D, vn, p_vn_dual, operators_coefficients, vort_v)
    ! this is not needed, since vort_v is on vertices in domain
    ! CALL sync_patch_array(SYNC_V, patch_2D, vort_v)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,je,start_edge_index,end_edge_index, this_vort_flux, &
!ICON_OMP  vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk, neighbor, il_v, ib_v, &
!ICON_OMP vertex_edge, ictr, il_e, ib_e) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_edge_index, end_edge_index)

      vort_flux(:,:,jb) = 0.0_wp

      edge_idx_loop: DO je =  start_edge_index, end_edge_index

        this_vort_flux(:,:) = 0.0_wp

        vertex1_idx = patch_2D%edges%vertex_idx(je,jb,1)
        vertex1_blk = patch_2D%edges%vertex_blk(je,jb,1)
        vertex2_idx = patch_2D%edges%vertex_idx(je,jb,2)
        vertex2_blk = patch_2D%edges%vertex_blk(je,jb,2)

        DO neighbor=1,2
          IF(neighbor==1) ictr = 0
          IF(neighbor==2) ictr = no_dual_edges

          il_v = patch_2D%edges%vertex_idx(je,jb,neighbor)
          ib_v = patch_2D%edges%vertex_blk(je,jb,neighbor)

          DO vertex_edge=1, patch_2D%verts%num_edges(il_v,ib_v)!no_dual_cell_edges

            ictr =ictr+1

            il_e = patch_2D%verts%edge_idx(il_v,ib_v,vertex_edge)
            ib_e = patch_2D%verts%edge_blk(il_v,ib_v,vertex_edge)

            DO jk = slev, patch_3d%p_patch_1d(1)%dolic_e(je,jb)

              this_vort_flux(jk, neighbor) =  this_vort_flux(jk, neighbor) + &
                & vn(il_e, jk, ib_e) * operators_coefficients%edge2edge_viavert_coeff(je,jk,jb,ictr)

            ENDDO

          END DO ! edges of this vertex

        END DO ! neighbor=1,2

        DO jk = slev, patch_3d%p_patch_1d(1)%dolic_e(je,jb)

          vort_flux(je,jk,jb) =  &
            &   this_vort_flux(jk,1) * patch_2D%verts%f_v(vertex1_idx, vertex1_blk)  &
            & + this_vort_flux(jk,2) * patch_2D%verts%f_v(vertex2_idx, vertex2_blk)

        ENDDO

      END DO edge_idx_loop

    END DO ! jb = all_edges%start_block, all_edges%end_block
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

  END SUBROUTINE linear_coriolis_3d_fast
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Optimized version of the nonlinear_coriolis_3d
  !! Note: intel -o3 will produce very different results when using this version,
  !!  the model exhibits great sensitivity using the AtlanticBoxACC setup
  !<Optimize_Used.done>
  SUBROUTINE nonlinear_coriolis_3d_fast(patch_3D, vn, p_vn_dual, vort_v, &
    & operators_coefficients, vort_flux)
    TYPE(t_patch_3D ),TARGET,INTENT(IN) :: patch_3D
    REAL(wp), INTENT(INOUT)                    :: vn(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(INout)  :: p_vn_dual(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp), INTENT(INOUT)                    :: vort_v   (nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(IN)          :: operators_coefficients
    REAL(wp), INTENT(INOUT)                    :: vort_flux(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: slev! , elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    INTEGER :: vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk
    REAL(wp) :: this_vort_flux(n_zlev, 2) ! for each of the two vertices

    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-----------------------------------------------------------------------
    slev    = 1
    ! elev    = n_zlev

    CALL rot_vertex_ocean_3d( patch_3D, vn, p_vn_dual, operators_coefficients, vort_v)
    ! this is not needed, since vort_v is on vertices in domain
    ! CALL sync_patch_array(SYNC_V, patch_2D, vort_v)

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb,jk,je,start_edge_index,end_edge_index, this_vort_flux, &
!ICON_OMP  vertex1_idx, vertex1_blk, vertex2_idx, vertex2_blk, neighbor, il_v, ib_v, &
!ICON_OMP vertex_edge, ictr, il_e, ib_e) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_edge_index, end_edge_index)

      vort_flux(:,:,jb) = 0.0_wp

      edge_idx_loop: DO je =  start_edge_index, end_edge_index

        this_vort_flux(:,:) = 0.0_wp

        vertex1_idx = patch_2D%edges%vertex_idx(je,jb,1)
        vertex1_blk = patch_2D%edges%vertex_blk(je,jb,1)
        vertex2_idx = patch_2D%edges%vertex_idx(je,jb,2)
        vertex2_blk = patch_2D%edges%vertex_blk(je,jb,2)

        DO neighbor=1,2
          IF(neighbor==1) ictr = 0
          IF(neighbor==2) ictr = no_dual_edges

          il_v = patch_2D%edges%vertex_idx(je,jb,neighbor)
          ib_v = patch_2D%edges%vertex_blk(je,jb,neighbor)

          DO vertex_edge=1, patch_2D%verts%num_edges(il_v,ib_v)!no_dual_cell_edges

            ictr =ictr+1

            il_e = patch_2D%verts%edge_idx(il_v,ib_v,vertex_edge)
            ib_e = patch_2D%verts%edge_blk(il_v,ib_v,vertex_edge)

            DO jk = slev, patch_3d%p_patch_1d(1)%dolic_e(je,jb)

              this_vort_flux(jk, neighbor) =  this_vort_flux(jk, neighbor) + &
                & vn(il_e, jk, ib_e) * operators_coefficients%edge2edge_viavert_coeff(je,jk,jb,ictr)

            ENDDO

          END DO ! edges of this vertex

        END DO ! neighbor=1,2

        DO jk = slev, patch_3d%p_patch_1d(1)%dolic_e(je,jb)

          vort_flux(je,jk,jb) =  &
            &   this_vort_flux(jk,1) * &
            &     (vort_v(vertex1_idx, jk, vertex1_blk) + patch_2D%verts%f_v(vertex1_idx, vertex1_blk))  &
            & + this_vort_flux(jk,2) * &
            &     (vort_v(vertex2_idx, jk, vertex2_blk) + patch_2D%verts%f_v(vertex2_idx, vertex2_blk))

        ENDDO

      END DO edge_idx_loop

    END DO ! jb = all_edges%start_block, all_edges%end_block
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL

  END SUBROUTINE nonlinear_coriolis_3d_fast
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Note: vn must habve been synced before this routine
  !! the resulting vort_v is synced,
  !! vort_flux id calculated on edges in_domain
  SUBROUTINE nonlinear_coriolis_3d_fast0(patch_3D, vn, p_vn_dual, vort_v, &
    & operators_coefficients, vort_flux)
    TYPE(t_patch_3D ),TARGET,INTENT(IN) :: patch_3D
    REAL(wp), INTENT(INOUT)                    :: vn(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(INout)  :: p_vn_dual(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp), INTENT(INOUT)                    :: vort_v   (nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(IN)          :: operators_coefficients
    REAL(wp), INTENT(INOUT)                    :: vort_flux(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    REAL(wp) :: vort_global(n_zlev)
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: patch_2D

    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-----------------------------------------------------------------------
    slev    = 1
    elev    = n_zlev

    CALL rot_vertex_ocean_3d( patch_3D, vn, p_vn_dual, operators_coefficients, vort_v)
    ! this is not needed
    ! CALL sync_patch_array(SYNC_V, patch_2D, vort_v)


! !$OMP PARALLEL
! !$OMP DO PRIVATE(jb,jk,je,start_edge_index,end_edge_index)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_edge_index, end_edge_index)


       vort_flux(:,:,jb) = 0.0_wp

       edge_idx_loop: DO je =  start_edge_index, end_edge_index

          DO neighbor=1,2
            IF(neighbor==1) ictr = 0
            IF(neighbor==2) ictr = no_dual_edges

            il_v = patch_2D%edges%vertex_idx(je,jb,neighbor)
            ib_v = patch_2D%edges%vertex_blk(je,jb,neighbor)

            DO jk = slev, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
              vort_global(jk) = (vort_v(il_v,jk,ib_v) + patch_2D%verts%f_v(il_v,ib_v))
            ENDDO

            DO vertex_edge=1, patch_2D%verts%num_edges(il_v,ib_v)!no_dual_cell_edges

              ictr =ictr+1

              il_e = patch_2D%verts%edge_idx(il_v,ib_v,vertex_edge)
              ib_e = patch_2D%verts%edge_blk(il_v,ib_v,vertex_edge)

              DO jk = slev, patch_3d%p_patch_1d(1)%dolic_e(je,jb)
                    vort_flux(je,jk,jb) =  vort_flux(je,jk,jb)+vn(il_e,jk,ib_e)*vort_global(jk) &
                         & * operators_coefficients%edge2edge_viavert_coeff(je,jk,jb,ictr)
              END DO

            END DO ! vertex_edge=1,
          END DO ! neighbor

        END DO edge_idx_loop
    END DO ! jb = all_edges%start_blockold, all_edges%end_block
! !$OMP END DO NOWAIT
! !$OMP END PARALLEL

  END SUBROUTINE nonlinear_coriolis_3d_fast0
  !-------------------------------------------------------------------------

   !-------------------------------------------------------------------------
  !>
  !! Note: vn must habve been synced before this routine
  !! the resulting vort_v is synced,
  !! vort_flux id calculated on edges in_domain
  !<Optimize_Used.done>
  SUBROUTINE nonlinear_coriolis_3d(patch_3D, vn, p_vn_dual, vort_v, &
    & operators_coefficients, vort_flux)
    TYPE(t_patch_3D ),TARGET,INTENT(IN) :: patch_3D
    REAL(wp), INTENT(INOUT)                    :: vn(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(INout)  :: p_vn_dual(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp), INTENT(INOUT)                    :: vort_v   (nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(IN)          :: operators_coefficients
    REAL(wp), INTENT(INOUT)                    :: vort_flux(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)

    !Local variables
    !TYPE(t_patch), POINTER         :: patch_2D
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    REAL(wp) :: vort_global
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: patch_2D

    !-----------------------------------------------------------------------
    IF (fast_performance_level > 20) THEN
      CALL nonlinear_coriolis_3d_fast(patch_3D, vn, p_vn_dual, vort_v, &
        & operators_coefficients, vort_flux)
      RETURN
    ELSEIF  (fast_performance_level > 12) THEN
      CALL nonlinear_coriolis_3d_fast0(patch_3D, vn, p_vn_dual, vort_v, &
        & operators_coefficients, vort_flux)
      RETURN
    ENDIF

    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-----------------------------------------------------------------------
    slev    = 1
    elev    = n_zlev

    CALL rot_vertex_ocean_3d( patch_3D, vn, p_vn_dual, operators_coefficients, vort_v)
    ! this is not needed
    ! CALL sync_patch_array(SYNC_V, patch_2D, vort_v)


! !$OMP PARALLEL
! !$OMP DO PRIVATE(jb,jk,je,start_edge_index,end_edge_index)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_edge_index, end_edge_index)

      level_loop: DO jk = slev, elev

        edge_idx_loop: DO je =  start_edge_index, end_edge_index

          IF (patch_3D%lsm_e(je,jk,jb) == sea) THEN

            vort_flux(je,jk,jb) = 0.0_wp 

            DO neighbor=1,2   
              IF(neighbor==1) ictr = 0
              IF(neighbor==2) ictr = no_dual_edges

              il_v = patch_2D%edges%vertex_idx(je,jb,neighbor)
              ib_v = patch_2D%edges%vertex_blk(je,jb,neighbor)

              vort_global = (vort_v(il_v,jk,ib_v) + patch_2D%verts%f_v(il_v,ib_v))

              DO vertex_edge=1, patch_2D%verts%num_edges(il_v,ib_v)!no_dual_cell_edges

                ictr =ictr+1 

                il_e = patch_2D%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = patch_2D%verts%edge_blk(il_v,ib_v,vertex_edge)

                vort_flux(je,jk,jb) =  vort_flux(je,jk,jb)+vn(il_e,jk,ib_e)*vort_global&
                                    &*operators_coefficients%edge2edge_viavert_coeff(je,jk,jb,ictr)
              END DO
            END DO
            ELSE
              vort_flux(je,jk,jb)= 0.0_wp
          ENDIF ! (v_base%lsm_e(je,jk,jb) <= sea_boundary)
        END DO edge_idx_loop
      END DO level_loop
    END DO ! jb = all_edges%start_blockold, all_edges%end_block
! !$OMP END DO NOWAIT
! !$OMP END PARALLEL

  END SUBROUTINE nonlinear_coriolis_3d
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  SUBROUTINE map_edges2cell_with_height_3d( patch_3D, vn_e, operators_coefficients, p_vn_c, h_e,&
    & opt_slev, opt_elev, subset_range)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    ! 3D case: h_e is surface elevation at edges
    TYPE(t_cartesian_coordinates),INTENT(inout):: p_vn_c(:,:,:)  ! outputput (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(in)                       :: h_e(:,:)       ! SW-case: h_e is thickness at edges
    TYPE(t_operator_coeff)                     :: operators_coefficients
    INTEGER, INTENT(in), OPTIONAL :: opt_slev       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_elev       ! optional vertical end level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
 
    !Local variables
    !INTEGER, PARAMETER :: no_primal_edges = 3
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: il_e, ib_e
    INTEGER :: jc, jb, jk, ie!,je
    REAL(wp) :: z_weight
    REAL(wp) :: z_thick_e

    TYPE(t_subset_range), POINTER :: all_cells 
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
   TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2D%cells%all
    ENDIF
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF

    IF ( iswm_oce == 1 ) THEN

      !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
        level_loop_swm: DO jk = slev, elev
          cell_idx_loop_swm: DO jc =  i_startidx_c, i_endidx_c
            !calculate velocity reconstruction at cell center
            z_weight           = 0.0_wp
            p_vn_c(jc,jk,jb)%x = 0.0_wp
            DO ie=1, no_primal_edges

              il_e = patch_2D%cells%edge_idx(jc,jb,ie)
              ib_e = patch_2D%cells%edge_blk(jc,jb,ie)

              z_thick_e =patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_e(il_e,jk,ib_e)&
              & + h_e(il_e,ib_e)  
              z_weight = z_weight + operators_coefficients%variable_vol_norm(jc,jk,jb,ie) * z_thick_e

              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
                & + operators_coefficients%edge2cell_coeff_cc_dyn(jc,jk,jb,ie)%x&
                & * vn_e(il_e,jk,ib_e)* z_thick_e
            END DO
            IF( z_weight/=0.0_wp)THEN
              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x / z_weight
            ELSE
              p_vn_c(jc,jk,jb)%x=0.0_wp
            ENDIF
          END DO cell_idx_loop_swm
        END DO level_loop_swm
      END DO ! jb = all_cells%start_block, all_cells%end_block

    ELSEIF( iswm_oce /= 1 ) THEN

      !Step 1: Calculation of Pv in cartesian coordinates
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

        !We are dealing with the surface layer first
        cell_idx_loop_top: DO jc =  i_startidx_c, i_endidx_c
          z_weight             = 0.0_wp
          p_vn_c(jc,slev,jb)%x = 0.0_wp
          DO ie=1, no_primal_edges

            il_e = patch_2D%cells%edge_idx(jc,jb,ie)
            ib_e = patch_2D%cells%edge_blk(jc,jb,ie)

            z_thick_e = patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_e(il_e,slev,ib_e)&
            & + h_e(il_e,ib_e) 
            z_weight = z_weight + operators_coefficients%variable_vol_norm(jc,slev,jb,ie) * z_thick_e

            p_vn_c(jc,slev,jb)%x = p_vn_c(jc,slev,jb)%x&
              & + operators_coefficients%edge2cell_coeff_cc_dyn(jc,1,jb,ie)%x&
              & * vn_e(il_e,slev,ib_e) * z_thick_e

          END DO

          IF(z_weight/=0.0_wp)THEN
            p_vn_c(jc,slev,jb)%x = p_vn_c(jc,slev,jb)%x / z_weight
          ELSE
            p_vn_c(jc,slev,jb)%x=0.0_wp
          ENDIF
        END DO cell_idx_loop_top

        !Now we calculate at the levels below the surface
        level_loop: DO jk = slev+1, elev
          cell_idx_loop: DO jc =  i_startidx_c, i_endidx_c
            p_vn_c(jc,jk,jb)%x = 0.0_wp
            !z_weight = 0.0_wp
            DO ie=1, no_primal_edges

              il_e = patch_2D%cells%edge_idx(jc,jb,ie)
              ib_e = patch_2D%cells%edge_blk(jc,jb,ie)
              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
                & + operators_coefficients%edge2cell_coeff_cc(jc,jk,jb,ie)%x&
                & * vn_e(il_e,jk,ib_e)
            END DO
          END DO cell_idx_loop
        END DO level_loop
      END DO ! jb = all_cells%start_block, all_cells%end_block
    ENDIF

    ! LL no sync required    

  END SUBROUTINE map_edges2cell_with_height_3d
  !----------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! map_edges2cell_without_height
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
!<Optimize_Used>
  SUBROUTINE map_edges2cell_no_height_3d( patch_3D, vn_e, operators_coefficients, p_vn_c, opt_slev, opt_elev, &
    &                                     subset_range)
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    TYPE(t_cartesian_coordinates)              :: p_vn_c(:,:,:)  ! output (nproma,n_zlev,alloc_cell_blocks)
                                                                 ! intent(inout) for nag compiler
    INTEGER, INTENT(in), OPTIONAL :: opt_slev       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_elev       ! optional vertical end level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: il_e, ib_e
    INTEGER :: jc, jb, jk, ie
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => patch_2D%cells%all
    ENDIF
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF

    !Calculation of Pv in cartesian coordinates
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      level_loop: DO jk = slev, elev
        cell_idx_loop: DO jc =  i_startidx_c, i_endidx_c
          !calculate velocity reconstruction at cell center
          p_vn_c(jc,jk,jb)%x = 0.0_wp

          DO ie=1, no_primal_edges
            il_e = patch_2D%cells%edge_idx(jc,jb,ie)
            ib_e = patch_2D%cells%edge_blk(jc,jb,ie)

            p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
              & + operators_coefficients%edge2cell_coeff_cc(jc,jk,jb,ie)%x&
              & * vn_e(il_e,jk,ib_e)
          END DO
          IF(operators_coefficients%fixed_vol_norm(jc,jk,jb)/=0.0_wp)THEN
            p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x/operators_coefficients%fixed_vol_norm(jc,jk,jb)
          ENDIF
        END DO cell_idx_loop
      END DO level_loop

    END DO ! jb = all_cells%start_block, all_cells%end_block
    ! LL no sync required

  END SUBROUTINE map_edges2cell_no_height_3d
  !-----------------------------------------------------------------------------
  SUBROUTINE map_edges2edges_viacell_3d_mlev( patch_3D, vn_e, operators_coefficients, out_vn_e,scalar, opt_slev, opt_elev, &
    &                                     subset_range)
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)       :: patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(INOUT)                    :: out_vn_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(IN), OPTIONAL             :: scalar(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    INTEGER, INTENT(in), OPTIONAL              :: opt_slev       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL              :: opt_elev       ! optional vertical end level
    TYPE(t_subset_range), TARGET,  OPTIONAL    :: subset_range
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr!, neighbor
    INTEGER :: je, jb, jk, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
     TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER         :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_edges => subset_range
    ELSE
      all_edges => patch_2D%edges%all
    ENDIF
    !-----------------------------------------------------------------------
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF

    IF(.NOT.PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)
 
        level_loop_e: DO jk = slev, elev
          edge_idx_loop: DO je =  start_edge_index, end_edge_index
          IF (patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
            ictr            = 0
            out_vn_e(je,jk,jb)= 0.0_wp
            !IF(patch_3D%lsm_e(je,jk,jb) == sea)THEN
            il_c = patch_2D%edges%cell_idx(je,jb,1)
            ib_c = patch_2D%edges%cell_blk(je,jb,1)
            thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

              thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)
              thick_frac = thick_edge/thick_cell
              out_vn_e(je,jk,jb) = out_vn_e(je,jk,jb) &
              &+vn_e(il_e,jk,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr))*thick_frac

            END DO

            ictr = no_primal_edges
            il_c = patch_2D%edges%cell_idx(je,jb,2)
            ib_c = patch_2D%edges%cell_blk(je,jb,2)
            thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

              thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)
              thick_frac = thick_edge/thick_cell
              out_vn_e(je,jk,jb) = out_vn_e(je,jk,jb) &
              &+vn_e(il_e,jk,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr))*thick_frac

            END DO
            ENDIF
          END DO edge_idx_loop
        END DO level_loop_e
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

    ELSEIF(PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)

        level_loop_e2: DO jk = slev, elev
          edge_idx_loop2: DO je =  start_edge_index, end_edge_index
          IF (patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
            ictr            = 0
            out_vn_e(je,jk,jb)= 0.0_wp
            il_c        = patch_2D%edges%cell_idx(je,jb,1)
            ib_c        = patch_2D%edges%cell_blk(je,jb,1)
            scalar_cell = scalar(il_c,jk,ib_c)
            thick_cell  = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

              DO ie=1, no_primal_edges
                ictr =ictr+1 
                il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
                ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
                thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)
                thick_frac = thick_edge/thick_cell

                out_vn_e(je,jk,jb) = out_vn_e(je,jk,jb) &
                &+vn_e(il_e,jk,ib_e)*scalar_cell   &
                &*(operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr))*thick_frac
              END DO


            ictr        = no_primal_edges
            il_c        = patch_2D%edges%cell_idx(je,jb,2)
            ib_c        = patch_2D%edges%cell_blk(je,jb,2)
            scalar_cell = scalar(il_c,jk,ib_c)
            thick_cell  = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

              DO ie=1, no_primal_edges
                ictr =ictr+1 
                il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
                ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
                thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)
                thick_frac = thick_edge/thick_cell

                out_vn_e(je,jk,jb) = out_vn_e(je,jk,jb) &
                &+vn_e(il_e,jk,ib_e)*scalar_cell   &
                &*(operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr))*thick_frac

              END DO

            ENDIF
          END DO edge_idx_loop2
        END DO level_loop_e2
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block
    ENDIF

  END SUBROUTINE map_edges2edges_viacell_3D_mlev
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE map_edges2edges_viacell_3d_1lev( patch_3D, vn_e, operators_coefficients, out_vn_e,scalar,scalar_e, level, &
    &                                     subset_range)
   
    TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(INOUT)                    :: out_vn_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(IN), OPTIONAL             :: scalar(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(IN), OPTIONAL             :: scalar_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    INTEGER, INTENT(in), OPTIONAL              :: level       ! optional vertical start level
    TYPE(t_subset_range), TARGET,  OPTIONAL    :: subset_range
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, jb, ie, lev
    REAL(wp) :: scalar_cell
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_edges => subset_range
    ELSE
      all_edges => patch_2D%edges%all
    ENDIF
    !--------------------------
    IF ( PRESENT(level) ) THEN
      lev = level
    ELSE
      lev = 1
    END IF

    IF(.NOT.PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)

          edge_idx_loop: DO je =  start_edge_index, end_edge_index
          IF (patch_3D%lsm_e(je,lev,jb) <= sea_boundary) THEN
            ictr          = 0
            out_vn_e(je,jb) = 0.0_wp
            il_c = patch_2D%edges%cell_idx(je,jb,1)
            ib_c = patch_2D%edges%cell_blk(je,jb,1)

            !thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)
            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
              !thick_edge=patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
              !thick_frac=thick_edge/thick_cell
              out_vn_e(je,jb) = out_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,lev,jb,ictr))!*thick_frac

            END DO

            ictr = no_primal_edges
            il_c = patch_2D%edges%cell_idx(je,jb,2)
            ib_c = patch_2D%edges%cell_blk(je,jb,2)

            !thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)
            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
              !thick_edge=patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
              !thick_frac=thick_edge/thick_cell
              out_vn_e(je,jb) = out_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*(operators_coefficients%edge2edge_viacell_coeff(je,lev,jb,ictr))!*thick_frac

            END DO

          ENDIF
          END DO edge_idx_loop
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

    ELSEIF(PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)

          edge_idx_loop2: DO je =  start_edge_index, end_edge_index
          IF (patch_3D%lsm_e(je,lev,jb) <= sea_boundary) THEN
            ictr = 0
            out_vn_e(je,jb) = 0.0_wp
            il_c        = patch_2D%edges%cell_idx(je,jb,1)
            ib_c        = patch_2D%edges%cell_blk(je,jb,1)
            scalar_cell = scalar(il_c,ib_c)

            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

              out_vn_e(je,jb) = out_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*scalar_cell &
              &  *(operators_coefficients%edge2edge_viacell_coeff(je,lev,jb,ictr))!*thick_frac

            END DO
            ictr        = no_primal_edges
            il_c        = patch_2D%edges%cell_idx(je,jb,2)
            ib_c        = patch_2D%edges%cell_blk(je,jb,2)
            scalar_cell = scalar(il_c,ib_c)
            !thick_cell  = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)

            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

              !thick_edge=patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
              !thick_frac=thick_edge/thick_cell
              out_vn_e(je,jb) = out_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*scalar_cell&
              &  *(operators_coefficients%edge2edge_viacell_coeff(je,lev,jb,ictr))!*thick_frac

            END DO
           ENDIF
           END DO edge_idx_loop2
      END DO
    ENDIF
  END SUBROUTINE map_edges2edges_viacell_3D_1lev
  !----------------------------------------------------------------------------- 

  !-----------------------------------------------------------------------------
!<Optimize_Used>
  SUBROUTINE map_edges2edges_viacell_3d_mlev_const_z( patch_3D, vn_e, operators_coefficients, out_vn_e)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    REAL(wp), INTENT(IN)                 :: vn_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(IN)   :: operators_coefficients
    REAL(wp), INTENT(INOUT)              :: out_vn_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, jb, jk, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    all_edges => patch_2D%edges%all
    slev = 1
    elev = n_zlev

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)

          edge_idx_loop_sfc: DO je =  start_edge_index, end_edge_index
          IF (patch_3D%lsm_e(je,slev,jb) == sea) THEN

            out_vn_e(je,slev,jb) = 0.0_wp

            ictr = 0
            il_c = patch_2D%edges%cell_idx(je,jb,1)
            ib_c = patch_2D%edges%cell_blk(je,jb,1)
            thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

              out_vn_e(je,slev,jb) = out_vn_e(je,slev,jb) &
              &+vn_e(il_e,slev,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,slev,jb,ictr)&
              &*thick_cell

            END DO

            ictr = no_primal_edges
            il_c = patch_2D%edges%cell_idx(je,jb,2)
            ib_c = patch_2D%edges%cell_blk(je,jb,2)
            thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

              out_vn_e(je,slev,jb) = out_vn_e(je,slev,jb) &
              &+vn_e(il_e,slev,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,slev,jb,ictr)&
              &*thick_cell

            END DO
            ENDIF
          END DO edge_idx_loop_sfc

        level_loop_e: DO jk = slev+1, elev
          edge_idx_loop: DO je =  start_edge_index, end_edge_index
          IF (patch_3D%lsm_e(je,jk,jb) == sea) THEN
            out_vn_e(je,jk,jb)= 0.0_wp

            ictr = 0
            il_c = patch_2D%edges%cell_idx(je,jb,1)
            ib_c = patch_2D%edges%cell_blk(je,jb,1)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

              out_vn_e(je,jk,jb) = out_vn_e(je,jk,jb) &
              &+vn_e(il_e,jk,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr)&
              &*thick_edge

            END DO

            ictr = no_primal_edges
            il_c = patch_2D%edges%cell_idx(je,jb,2)
            ib_c = patch_2D%edges%cell_blk(je,jb,2)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

              out_vn_e(je,jk,jb) = out_vn_e(je,jk,jb) &
              &+vn_e(il_e,jk,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr)&
              &*thick_edge

            END DO
            ENDIF
          END DO edge_idx_loop
        END DO level_loop_e
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

  END SUBROUTINE map_edges2edges_viacell_3D_mlev_const_z
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
!<Optimize_Used>
  SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs( patch_3D, vn_e, operators_coefficients, out_vn_e, scalar)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    REAL(wp), INTENT(IN)                 :: vn_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(IN)   :: operators_coefficients
    REAL(wp), INTENT(INOUT)              :: out_vn_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(IN)       :: scalar(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, jb, jk, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: patch_2D

    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    IF ( patch_2D%cells%max_connectivity == 3) THEN
      IF (fast_performance_level > 10 ) THEN
        CALL map_edges2edges_viacell_3d_mlev_constZs_fast( patch_3D, vn_e, operators_coefficients, out_vn_e, scalar)
        RETURN
      ENDIF
    ENDIF

    !-----------------------------------------------------------------------
    all_edges => patch_2D%edges%all
    slev = 1
    elev = n_zlev

    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)

        edge_idx_loop_sfc2: DO je =  start_edge_index, end_edge_index
        IF (patch_3D%lsm_e(je,slev,jb) == sea) THEN

          out_vn_e(je,slev,jb) = 0.0_wp

          ictr = 0
          il_c       = patch_2D%edges%cell_idx(je,jb,1)
          ib_c       = patch_2D%edges%cell_blk(je,jb,1)
          thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)
          scalar_cell= scalar(il_c,slev,ib_c)

          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

            out_vn_e(je,slev,jb) = out_vn_e(je,slev,jb)&
            &+vn_e(il_e,slev,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,slev,jb,ictr)&
            &*thick_cell *scalar_cell

          END DO

          ictr = no_primal_edges
          il_c = patch_2D%edges%cell_idx(je,jb,2)
          ib_c = patch_2D%edges%cell_blk(je,jb,2)
          thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)
          scalar_cell = scalar(il_c,slev,ib_c)

          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

            out_vn_e(je,slev,jb) = out_vn_e(je,slev,jb) &
            &+vn_e(il_e,slev,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,slev,jb,ictr)&
            &*thick_cell *scalar_cell

          END DO
          ENDIF
        END DO edge_idx_loop_sfc2

      level_loop_e2: DO jk = slev+1, elev
        edge_idx_loop2: DO je =  start_edge_index, end_edge_index
        IF (patch_3D%lsm_e(je,jk,jb) == sea) THEN
          out_vn_e(je,jk,jb)= 0.0_wp

          ictr = 0
          il_c        = patch_2D%edges%cell_idx(je,jb,1)
          ib_c        = patch_2D%edges%cell_blk(je,jb,1)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)
          scalar_cell = scalar(il_c,jk,ib_c)

          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
            thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

            out_vn_e(je,jk,jb) = out_vn_e(je,jk,jb) &
            &+vn_e(il_e,jk,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr)&
            &*thick_edge*scalar_cell

          END DO

          ictr        = no_primal_edges
          il_c        = patch_2D%edges%cell_idx(je,jb,2)
          ib_c        = patch_2D%edges%cell_blk(je,jb,2)!thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)
          scalar_cell = scalar(il_c,jk,ib_c)

          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
            thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

            out_vn_e(je,jk,jb) = out_vn_e(je,jk,jb) &
            &+vn_e(il_e,jk,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr)&
            &*thick_edge*scalar_cell

          END DO
          ENDIF
        END DO edge_idx_loop2
      END DO level_loop_e2
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

  END SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! the map_edges2edges_viacell_3d_mlev_constZs optimized for triangles
!<Optimize_Used>
  SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs_fast( patch_3D, vn_e, operators_coefficients, out_vn_e,scalar)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    REAL(wp), INTENT(IN)                 :: vn_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(IN)   :: operators_coefficients
    REAL(wp), INTENT(INOUT)              :: out_vn_e(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(IN)                 :: scalar(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, jb, level

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_11_index, edge_12_index, edge_13_index ! edges of cell_1
    INTEGER :: edge_11_block, edge_12_block, edge_13_block
    INTEGER :: edge_21_index, edge_22_index, edge_23_index ! edges of cell_2
    INTEGER :: edge_21_block, edge_22_block, edge_23_block

    REAL(wp) :: scalar_cell
    REAL(wp), POINTER :: coeffs(:,:,:,:)
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')

    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_edges => patch_2D%edges%all
    slev = 1
    elev = n_zlev
    coeffs => operators_coefficients%edge2edge_viacell_coeff
    !-----------------------------------------------------------------------

    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)

      DO je =  start_edge_index, end_edge_index

        IF (patch_3D%surface_edge_sea_land_mask(je,jb) /= sea) CYCLE

        cell_1_index = patch_2D%edges%cell_idx(je,jb,1)
        cell_1_block = patch_2D%edges%cell_blk(je,jb,1)
        cell_2_index = patch_2D%edges%cell_idx(je,jb,2)
        cell_2_block = patch_2D%edges%cell_blk(je,jb,2)

        edge_11_index = patch_2D%cells%edge_idx(cell_1_index, cell_1_block, 1)
        edge_12_index = patch_2D%cells%edge_idx(cell_1_index, cell_1_block, 2)
        edge_13_index = patch_2D%cells%edge_idx(cell_1_index, cell_1_block, 3)
        edge_11_block = patch_2D%cells%edge_blk(cell_1_index, cell_1_block, 1)
        edge_12_block = patch_2D%cells%edge_blk(cell_1_index, cell_1_block, 2)
        edge_13_block = patch_2D%cells%edge_blk(cell_1_index, cell_1_block, 3)

        edge_21_index = patch_2D%cells%edge_idx(cell_2_index, cell_2_block, 1)
        edge_22_index = patch_2D%cells%edge_idx(cell_2_index, cell_2_block, 2)
        edge_23_index = patch_2D%cells%edge_idx(cell_2_index, cell_2_block, 3)
        edge_21_block = patch_2D%cells%edge_blk(cell_2_index, cell_2_block, 1)
        edge_22_block = patch_2D%cells%edge_blk(cell_2_index, cell_2_block, 2)
        edge_23_block = patch_2D%cells%edge_blk(cell_2_index, cell_2_block, 3)

        ! top level
        out_vn_e(je, slev, jb) =  &
          &  ( vn_e(edge_11_index, slev, edge_11_block) * coeffs(je, slev, jb, 1) +      &
          &    vn_e(edge_12_index, slev, edge_12_block) * coeffs(je, slev, jb, 2) +      &
          &    vn_e(edge_13_index, slev, edge_13_block) * coeffs(je, slev, jb, 3)        &
          &  )  * patch_3D%p_patch_1D(1)%prism_thick_c(cell_1_index, slev, cell_1_block) &
          &     *                               scalar(cell_1_index, slev, cell_1_block) &
          & + &
          &  ( vn_e(edge_21_index, slev, edge_21_block) * coeffs(je, slev, jb, 4) +      &
          &    vn_e(edge_22_index, slev, edge_22_block) * coeffs(je, slev, jb, 5) +      &
          &    vn_e(edge_23_index, slev, edge_23_block) * coeffs(je, slev, jb, 6)        &
          &  )  * patch_3D%p_patch_1D(1)%prism_thick_c(cell_2_index, slev, cell_2_block) &
          &     *                               scalar(cell_2_index, slev, cell_2_block)

        ! next levels
        DO level = slev+1, patch_3d%p_patch_1d(1)%dolic_e(je,jb)

          out_vn_e(je, level, jb) =  &
            & (  vn_e(edge_11_index, level, edge_11_block) * coeffs(je, level, jb, 1)               &
            &       * patch_3D%p_patch_1D(1)%prism_thick_e(edge_11_index, level, edge_11_block)  +  &
            &    vn_e(edge_12_index, level, edge_12_block) * coeffs(je, level, jb, 2)               &
            &       * patch_3D%p_patch_1D(1)%prism_thick_e(edge_12_index, level, edge_12_block)  +  &
            &    vn_e(edge_13_index, level, edge_13_block) * coeffs(je, level, jb, 3)               &
            &       * patch_3D%p_patch_1D(1)%prism_thick_e(edge_13_index, level, edge_13_block)     &
            & ) * scalar(cell_1_index, level, cell_1_block)                                         &
            & + &
            & (  vn_e(edge_21_index, level, edge_21_block) * coeffs(je, level, jb, 4)               &
            &       * patch_3D%p_patch_1D(1)%prism_thick_e(edge_21_index, level, edge_21_block)  +  &
            &    vn_e(edge_22_index, level, edge_22_block) * coeffs(je, level, jb, 5)               &
            &       * patch_3D%p_patch_1D(1)%prism_thick_e(edge_22_index, level, edge_22_block)  +  &
            &    vn_e(edge_23_index, level, edge_23_block) * coeffs(je, level, jb, 6)               &
            &       * patch_3D%p_patch_1D(1)%prism_thick_e(edge_23_index, level, edge_23_block)     &
            & ) * scalar(cell_2_index, level, cell_2_block)

        END DO ! next levels

      END DO

    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

  END SUBROUTINE map_edges2edges_viacell_3d_mlev_constZs_fast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE map_edges2edges_viacell_3D_1lev_constZ( patch_3D, vn_e, operators_coefficients, out_vn_e)!&   subset_range)
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)       :: patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(INOUT)                    :: out_vn_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, jb, jk, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: patch_2D

    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    IF ( patch_2D%cells%max_connectivity == 3) THEN
      IF (fast_performance_level > 10 ) THEN
        CALL map_edges2edges_viacell_2d_1lev_constZ_fast( patch_3D, vn_e, operators_coefficients, out_vn_e )
        RETURN
      ENDIF
    ENDIF
    !-----------------------------------------------------------------------

    all_edges => patch_2D%edges%all
    slev = 1
    elev = n_zlev

    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)

        edge_idx_loop_sfc: DO je =  start_edge_index, end_edge_index
        IF (patch_3D%lsm_e(je,slev,jb) == sea) THEN

          out_vn_e(je,jb) = 0.0_wp

          ictr = 0
          il_c = patch_2D%edges%cell_idx(je,jb,1)
          ib_c = patch_2D%edges%cell_blk(je,jb,1)
          thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)

          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

            out_vn_e(je,jb) = out_vn_e(je,jb) &
            &+vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,slev,jb,ictr)&
            &*thick_cell

          END DO

          ictr = no_primal_edges
          il_c = patch_2D%edges%cell_idx(je,jb,2)
          ib_c = patch_2D%edges%cell_blk(je,jb,2)
          thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)

          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

            out_vn_e(je,jb) = out_vn_e(je,jb) &
            &+vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,slev,jb,ictr)&
            &*thick_cell

          END DO
          ENDIF
        END DO edge_idx_loop_sfc

      level_loop_e: DO jk = slev+1, elev
        edge_idx_loop: DO je =  start_edge_index, end_edge_index
        IF (patch_3D%lsm_e(je,jk,jb) == sea) THEN

          ictr = 0
          il_c = patch_2D%edges%cell_idx(je,jb,1)
          ib_c = patch_2D%edges%cell_blk(je,jb,1)

          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
            thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

            out_vn_e(je,jb) = out_vn_e(je,jb) &
            &+vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr)&
            &*thick_edge

          END DO

          ictr = no_primal_edges
          il_c = patch_2D%edges%cell_idx(je,jb,2)
          ib_c = patch_2D%edges%cell_blk(je,jb,2)

          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
            thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

            out_vn_e(je,jb) = out_vn_e(je,jb) &
            &+vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr)&
            &*thick_edge

          END DO
          ENDIF
        END DO edge_idx_loop
      END DO level_loop_e
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

  END SUBROUTINE map_edges2edges_viacell_3D_1lev_constZ
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  SUBROUTINE map_edges2edges_viacell_3d_1lev_constZs( patch_3D, vn_e, operators_coefficients, out_vn_e,scalar)!&   subset_range)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)       :: patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(INOUT)                    :: out_vn_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(IN)                       :: scalar(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, jb, jk, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: patch_2D

    !-----------------------------------------------------------------------
    patch_2D  => patch_3D%p_patch_2D(1)
    all_edges => patch_2D%edges%all
    slev = 1
    elev = n_zlev

    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)

        edge_idx_loop_sfc2: DO je =  start_edge_index, end_edge_index
        IF (patch_3D%lsm_e(je,slev,jb) == sea) THEN

          out_vn_e(je,jb) = 0.0_wp

          ictr = 0
          il_c = patch_2D%edges%cell_idx(je,jb,1)
          ib_c = patch_2D%edges%cell_blk(je,jb,1)
          thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)
          scalar_cell= scalar(il_c,ib_c)

          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

            out_vn_e(je,jb) = out_vn_e(je,jb) &
            &+vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,slev,jb,ictr)&
            &*thick_cell*scalar_cell

          END DO

          ictr = no_primal_edges
          il_c = patch_2D%edges%cell_idx(je,jb,2)
          ib_c = patch_2D%edges%cell_blk(je,jb,2)
          thick_cell = patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)
          scalar_cell= scalar(il_c,ib_c)
          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

            out_vn_e(je,jb) = out_vn_e(je,jb) &
            &+vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,slev,jb,ictr)&
            &*thick_cell*scalar_cell

          END DO
          ENDIF
        END DO edge_idx_loop_sfc2

      level_loop_e2: DO jk = slev+1, elev
        edge_idx_loop2: DO je =  start_edge_index, end_edge_index
        IF (patch_3D%lsm_e(je,jk,jb) == sea) THEN

          ictr = 0
          il_c = patch_2D%edges%cell_idx(je,jb,1)
          ib_c = patch_2D%edges%cell_blk(je,jb,1)
          scalar_cell= scalar(il_c,ib_c)

          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
            thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

            out_vn_e(je,jb) = out_vn_e(je,jb) &
            &+vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr)&
            &*thick_edge*scalar_cell

          END DO

          ictr = no_primal_edges
          il_c = patch_2D%edges%cell_idx(je,jb,2)
          ib_c = patch_2D%edges%cell_blk(je,jb,2)
          scalar_cell= scalar(il_c,ib_c)
          DO ie=1, no_primal_edges
            ictr = ictr+1
            il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
            ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
            thick_edge = patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

            out_vn_e(je,jb) = out_vn_e(je,jb) &
            &+vn_e(il_e,ib_e)*operators_coefficients%edge2edge_viacell_coeff(je,jk,jb,ictr)&
            &*thick_edge*scalar_cell

          END DO
          ENDIF
        END DO edge_idx_loop2
      END DO level_loop_e2
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

  END SUBROUTINE map_edges2edges_viacell_3d_1lev_constZs
  !-------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE map_edges2edges_viacell_2d_1lev_constZ_fast( patch_3D, in_vn_e, operators_coefficients, out_vn_e )!&   subset_range)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)       :: patch_3D
    REAL(wp), INTENT(in)                       :: in_vn_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(INOUT)                    :: out_vn_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, jb, start_edge_index, end_edge_index

    REAL(wp), POINTER :: all_coeffs(:,:,:)

    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    IF (no_primal_edges /= 3) &
      CALL finish ('map_edges2edges_viacell triangle version', 'no_primal_edges /= 3')
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    edges_inDomain    => patch_2D%edges%in_domain
    all_coeffs        => operators_coefficients%edge2edge_viacell_coeff_all

    DO jb = edges_inDomain%start_block, edges_inDomain%end_block
      CALL get_index_range(edges_inDomain, jb, start_edge_index, end_edge_index)

      DO je = start_edge_index, end_edge_index

        out_vn_e(je,jb) = 0.0_wp

        IF (patch_3D%lsm_e(je,1,jb) == sea) THEN ! this if should be removed

          ! get the two cells of the edge
          cell_1_index = patch_2D%edges%cell_idx(je,jb,1)
          cell_1_block = patch_2D%edges%cell_blk(je,jb,1)
          cell_2_index = patch_2D%edges%cell_idx(je,jb,2)
          cell_2_block = patch_2D%edges%cell_blk(je,jb,2)

          ! get the six edges of the two cells
          edge_1_1_index = patch_2D%cells%edge_idx(cell_1_index, cell_1_block, 1)
          edge_1_2_index = patch_2D%cells%edge_idx(cell_1_index, cell_1_block, 2)
          edge_1_3_index = patch_2D%cells%edge_idx(cell_1_index, cell_1_block, 3)
          edge_2_1_index = patch_2D%cells%edge_idx(cell_2_index, cell_2_block, 1)
          edge_2_2_index = patch_2D%cells%edge_idx(cell_2_index, cell_2_block, 2)
          edge_2_3_index = patch_2D%cells%edge_idx(cell_2_index, cell_2_block, 3)
          edge_1_1_block = patch_2D%cells%edge_blk(cell_1_index, cell_1_block, 1)
          edge_1_2_block = patch_2D%cells%edge_blk(cell_1_index, cell_1_block, 2)
          edge_1_3_block = patch_2D%cells%edge_blk(cell_1_index, cell_1_block, 3)
          edge_2_1_block = patch_2D%cells%edge_blk(cell_2_index, cell_2_block, 1)
          edge_2_2_block = patch_2D%cells%edge_blk(cell_2_index, cell_2_block, 2)
          edge_2_3_block = patch_2D%cells%edge_blk(cell_2_index, cell_2_block, 3)

          out_vn_e(je,jb) = &
                          in_vn_e(edge_1_1_index, edge_1_1_block) * all_coeffs(1, je, jb) + &
                          in_vn_e(edge_1_2_index, edge_1_2_block) * all_coeffs(2, je, jb) + &
                          in_vn_e(edge_1_3_index, edge_1_3_block) * all_coeffs(3, je, jb) + &
                          in_vn_e(edge_2_1_index, edge_2_1_block) * all_coeffs(4, je, jb) + &
                          in_vn_e(edge_2_2_index, edge_2_2_block) * all_coeffs(5, je, jb) + &
                          in_vn_e(edge_2_3_index, edge_2_3_block) * all_coeffs(6, je, jb)

        ENDIF
      END DO
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

  END SUBROUTINE map_edges2edges_viacell_2d_1lev_constZ_fast
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! single precision version of the map_edges2edges_viacell_2d_1lev_constZ_fast
  SUBROUTINE map_edges2edges_viacell_2d_constZ_sp( patch_3D, in_vn_e, operators_coefficients, out_vn_e )!&   subset_range)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)       :: patch_3D
    REAL(sp), INTENT(in)                       :: in_vn_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_solverCoeff_singlePrecision), INTENT(in)         :: operators_coefficients
    REAL(sp), INTENT(INOUT)                    :: out_vn_e(nproma,patch_3D%p_patch_2D(1)%nblks_e)

    INTEGER :: cell_1_index, cell_2_index, cell_1_block, cell_2_block
    INTEGER :: edge_1_1_index, edge_1_2_index, edge_1_3_index
    INTEGER :: edge_2_1_index, edge_2_2_index, edge_2_3_index
    INTEGER :: edge_1_1_block, edge_1_2_block, edge_1_3_block
    INTEGER :: edge_2_1_block, edge_2_2_block, edge_2_3_block
    INTEGER :: je, jb, start_edge_index, end_edge_index

    REAL(sp), POINTER :: all_coeffs(:,:,:)

    TYPE(t_subset_range), POINTER :: edges_inDomain
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)

    edges_inDomain    => patch_2D%edges%in_domain
    all_coeffs        => operators_coefficients%edge2edge_viacell_coeff_all

    DO jb = edges_inDomain%start_block, edges_inDomain%end_block
      CALL get_index_range(edges_inDomain, jb, start_edge_index, end_edge_index)

      DO je = start_edge_index, end_edge_index

        out_vn_e(je,jb) = 0.0_wp

        IF (patch_3D%lsm_e(je,1,jb) == sea) THEN ! this if should be removed

          ! get the two cells of the edge
          cell_1_index = patch_2D%edges%cell_idx(je,jb,1)
          cell_1_block = patch_2D%edges%cell_blk(je,jb,1)
          cell_2_index = patch_2D%edges%cell_idx(je,jb,2)
          cell_2_block = patch_2D%edges%cell_blk(je,jb,2)

          ! get the six edges of the two cells
          edge_1_1_index = patch_2D%cells%edge_idx(cell_1_index, cell_1_block, 1)
          edge_1_2_index = patch_2D%cells%edge_idx(cell_1_index, cell_1_block, 2)
          edge_1_3_index = patch_2D%cells%edge_idx(cell_1_index, cell_1_block, 3)
          edge_2_1_index = patch_2D%cells%edge_idx(cell_2_index, cell_2_block, 1)
          edge_2_2_index = patch_2D%cells%edge_idx(cell_2_index, cell_2_block, 2)
          edge_2_3_index = patch_2D%cells%edge_idx(cell_2_index, cell_2_block, 3)
          edge_1_1_block = patch_2D%cells%edge_blk(cell_1_index, cell_1_block, 1)
          edge_1_2_block = patch_2D%cells%edge_blk(cell_1_index, cell_1_block, 2)
          edge_1_3_block = patch_2D%cells%edge_blk(cell_1_index, cell_1_block, 3)
          edge_2_1_block = patch_2D%cells%edge_blk(cell_2_index, cell_2_block, 1)
          edge_2_2_block = patch_2D%cells%edge_blk(cell_2_index, cell_2_block, 2)
          edge_2_3_block = patch_2D%cells%edge_blk(cell_2_index, cell_2_block, 3)

          out_vn_e(je,jb) = &
                          in_vn_e(edge_1_1_index, edge_1_1_block) * all_coeffs(1, je, jb) + &
                          in_vn_e(edge_1_2_index, edge_1_2_block) * all_coeffs(2, je, jb) + &
                          in_vn_e(edge_1_3_index, edge_1_3_block) * all_coeffs(3, je, jb) + &
                          in_vn_e(edge_2_1_index, edge_2_1_block) * all_coeffs(4, je, jb) + &
                          in_vn_e(edge_2_2_index, edge_2_2_block) * all_coeffs(5, je, jb) + &
                          in_vn_e(edge_2_3_index, edge_2_3_block) * all_coeffs(6, je, jb)

        ENDIF
      END DO
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

  END SUBROUTINE map_edges2edges_viacell_2d_constZ_sp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  SUBROUTINE map_edges2edges_viavert_3D(patch_3D, vn, p_vn_dual,operators_coefficients, vort_flux)
    TYPE(t_patch_3D ),TARGET,INTENT(IN)        :: patch_3D
    REAL(wp), INTENT(INOUT)                    :: vn(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(INout)  :: p_vn_dual(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(IN)          :: operators_coefficients
    REAL(wp), INTENT(INOUT)                    :: vort_flux(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)

    !Local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: il_e, ib_e
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_edges => patch_2D%edges%all
    !-----------------------------------------------------------------------
    slev    = 1
    elev    = n_zlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,start_edge_index,end_edge_index)
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)

      level_loop: DO jk = slev, elev

        edge_idx_loop: DO je =  start_edge_index, end_edge_index

          IF (patch_3D%lsm_e(je,jk,jb) == sea) THEN

            vort_flux(je,jk,jb) = 0.0_wp 

            DO neighbor=1,2   
              IF(neighbor==1) ictr = 0
              IF(neighbor==2) ictr = no_dual_edges

              il_v = patch_2D%edges%vertex_idx(je,jb,neighbor)
              ib_v = patch_2D%edges%vertex_blk(je,jb,neighbor)

              DO vertex_edge=1, patch_2D%verts%num_edges(il_v,ib_v)!no_dual_cell_edges

                ictr =ictr+1 

                il_e = patch_2D%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = patch_2D%verts%edge_blk(il_v,ib_v,vertex_edge)

                vort_flux(je,jk,jb) =  vort_flux(je,jk,jb)+vn(il_e,jk,ib_e)&
                                    &*operators_coefficients%edge2edge_viavert_coeff(je,jk,jb,ictr)
              END DO
            END DO

            ELSE
              vort_flux(je,jk,jb)= 0.0_wp
          ENDIF ! (v_base%lsm_e(je,jk,jb) <= sea_boundary)
        END DO edge_idx_loop
      END DO level_loop
    END DO ! jb = all_edges%start_block, all_edges%end_block
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ! LL: no sync required

  END SUBROUTINE map_edges2edges_viavert_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized by LL, result not synced
!<Optimize_Used>
  SUBROUTINE map_cell2edges_3D_mlevels( patch_3D, p_vn_c, ptp_vn, operators_coefficients,&
                                   & opt_slev, opt_elev, subset_range )
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:,:)    ! input vector (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)                      :: ptp_vn(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff)                     :: operators_coefficients
    INTEGER, INTENT(in), OPTIONAL              :: opt_slev        ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL              :: opt_elev        ! optional vertical end level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range

    !Local variables
    INTEGER :: slev, elev
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, jb, jk
    INTEGER :: il_c1,ib_c1, il_c2,ib_c2
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    edges_in_domain   => patch_2D%edges%in_domain

    ptp_vn(:,:,:) = 0.0_wp

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF

    ! calculation of transposed P^TPv from Pv (incart coord)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_edge_index, end_edge_index)

      level_loop_e: DO jk = slev, elev
        edge_idx_loop: DO je =  start_edge_index, end_edge_index

          IF(patch_3D%lsm_e(je,jk,jb) <= sea_boundary)THEN

            !Get indices of two adjacent triangles
            il_c1 = patch_2D%edges%cell_idx(je,jb,1)
            ib_c1 = patch_2D%edges%cell_blk(je,jb,1)
            il_c2 = patch_2D%edges%cell_idx(je,jb,2)
            ib_c2 = patch_2D%edges%cell_blk(je,jb,2)
            ptp_vn(je,jk,jb) =&
              & DOT_PRODUCT(p_vn_c(il_c1,jk,ib_c1)%x,&
              & operators_coefficients%edge2cell_coeff_cc_t(je,jk,jb,1)%x)&
              & +DOT_PRODUCT(p_vn_c(il_c2,jk,ib_c2)%x,&
              & operators_coefficients%edge2cell_coeff_cc_t(je,jk,jb,2)%x)

          ELSE
            ptp_vn(je,jk,jb) = 0.0_wp
          ENDIF

        END DO edge_idx_loop
      END DO level_loop_e
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

    ! sync the result if necessary
    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
       & CALL sync_patch_array(SYNC_E, patch_2D, ptp_vn)
    ENDIF
       
    !stop
  END SUBROUTINE map_cell2edges_3D_mlevels
  !-----------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized by LL, result not synced
!<Optimize_Used>
  SUBROUTINE map_cell2edges_3D_1level( patch_3D, p_vn_c, ptp_vn,operators_coefficients, level, subset_range )
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:)    ! input vector (nproma,n_zlev,alloc_cell_blocks)
    REAL(wp), INTENT(inout)                      :: ptp_vn(:,:)    ! output vector (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff)                     :: operators_coefficients
    INTEGER, INTENT(in) :: level          ! vertical level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range

    !Local variables
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: je, jb
    INTEGER :: il_c1,ib_c1, il_c2,ib_c2
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    edges_in_domain   => patch_2D%edges%in_domain
    ptp_vn(:,:) = 0.0_wp

    ! calculation of transposed P^TPv from Pv (incart coord)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, start_edge_index, end_edge_index)
      edge_idx_loop: DO je =  start_edge_index, end_edge_index

        IF(patch_3D%lsm_e(je,level,jb) <= sea_boundary)THEN
          !Get indices of two adjacent triangles
          il_c1 = patch_2D%edges%cell_idx(je,jb,1)
          ib_c1 = patch_2D%edges%cell_blk(je,jb,1)
          il_c2 = patch_2D%edges%cell_idx(je,jb,2)
          ib_c2 = patch_2D%edges%cell_blk(je,jb,2)
          ptp_vn(je,jb) =&
            & DOT_PRODUCT(p_vn_c(il_c1,ib_c1)%x,&
            & operators_coefficients%edge2cell_coeff_cc_t(je,level,jb,1)%x)&
            & +DOT_PRODUCT(p_vn_c(il_c2,ib_c2)%x,&
            & operators_coefficients%edge2cell_coeff_cc_t(je,level,jb,2)%x)
       ELSE
          ptp_vn(je,jb) = 0.0_wp
        ENDIF

      END DO edge_idx_loop
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

    ! sync the result if necessary
    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
       & CALL sync_patch_array(SYNC_E, patch_2D, ptp_vn)
    ENDIF
       
    !stop
  END SUBROUTINE map_cell2edges_3D_1level
  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE nonlinear_coriolis_3d_old(patch_3D, vn, p_vn_dual, vort_v, &
    & operators_coefficients, vort_flux)

    TYPE(t_patch_3D ),TARGET,INTENT(IN):: patch_3D
    REAL(wp), INTENT(inout)                   :: vn       (nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(inout) :: p_vn_dual(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp), INTENT(inout)                   :: vort_v   (nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in)         :: operators_coefficients
    REAL(wp), INTENT(inout)                   :: vort_flux(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)

    !Local variables
    !
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jk, jb, jev,je
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: il_v1, il_v2, ib_v1, ib_v2
    TYPE(t_cartesian_coordinates) :: u_v1_cc, u_v2_cc
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_edges => patch_2D%edges%all
    !-----------------------------------------------------------------------
    slev         = 1
    elev         = n_zlev
    !CALL map_edges2vert_3d(patch_3D%p_patch_2D(1), vn, operators_coefficients%edge2vert_coeff_cc, &
    !   & p_vn_dual)

    CALL rot_vertex_ocean_3d( patch_3D, vn, p_vn_dual, operators_coefficients, vort_v)
    CALL sync_patch_array(SYNC_V, patch_2D, vort_v)

! !$OMP PARALLEL
! !$OMP DO PRIVATE(jb,jk,je,start_edge_index,end_edge_index,il_v1,ib_v1,il_v2,ib_v2, u_v1_cc, u_v2_cc)
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)

      level_loop: DO jk = slev, elev

        edge_idx_loop: DO je =  start_edge_index, end_edge_index

          IF (patch_3D%lsm_e(je,jk,jb) == sea) THEN
            !Get indices of two adjacent vertices
            il_v1 = patch_2D%edges%vertex_idx(je,jb,1)
            ib_v1 = patch_2D%edges%vertex_blk(je,jb,1)
            il_v2 = patch_2D%edges%vertex_idx(je,jb,2)
            ib_v2 = patch_2D%edges%vertex_blk(je,jb,2)

            !Multiply velocity reconstruction at vertex by vorticity
            u_v1_cc%x = p_vn_dual(il_v1,jk,ib_v1)%x &
            & *(vort_v(il_v1,jk,ib_v1) + patch_2D%verts%f_v(il_v1,ib_v1))
            u_v2_cc%x = p_vn_dual(il_v2,jk,ib_v2)%x &
            & *(vort_v(il_v2,jk,ib_v2) + patch_2D%verts%f_v(il_v2,ib_v2))

            !calculate finall vortex flux by mapping the vorticity-velocity-product
            !from vertices back to edges
            vort_flux(je,jk,jb) =( &
               & - DOT_PRODUCT(u_v2_cc%x, operators_coefficients%edge2vert_coeff_cc_t(je,jk,jb,2)%x)&
               & + DOT_PRODUCT(u_v1_cc%x, operators_coefficients%edge2vert_coeff_cc_t(je,jk,jb,1)%x))
          ELSE
            vort_flux(je,jk,jb)= 0.0_wp
          ENDIF ! (v_base%lsm_e(je,jk,jb) <= sea_boundary)
        END DO edge_idx_loop
      END DO level_loop
    END DO ! jb = all_edges%start_block, all_edges%end_block
! !$OMP END DO NOWAIT
! !$OMP END PARALLEL
  ! LL: no sync required

  END SUBROUTINE nonlinear_coriolis_3d_old
  !-------------------------------------------------------------------------

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
! !     INTEGER :: il_c1, ib_c1, il_c2, ib_c2
! !     INTEGER :: je, jb
! !     TYPE(t_subset_range), POINTER :: edges_in_domain
! !     TYPE(t_patch), POINTER        :: patch_2D
! !     !-----------------------------------------------------------------------
! !     patch_2D   => patch_3D%p_patch_2D(1)
! !     !-----------------------------------------------------------------------
! !     !CALL message (TRIM(routine), 'start')
! !     edges_in_domain   => patch_2D%edges%in_domain
! ! 
! !     ! calculation of transposed P^TPv from Pv (incart coord)
! !     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, jb, start_edge_index, end_edge_index)
! ! 
! !       edge_idx_loop: DO je =  start_edge_index, end_edge_index
! ! 
! !         IF(patch_3D%lsm_e(je,1,jb) <= sea_boundary)THEN
! !           !Get indices of two adjacent triangles
! !           il_c1 = patch_2D%edges%cell_idx(je,jb,1)
! !           ib_c1 = patch_2D%edges%cell_blk(je,jb,1)
! !           il_c2 = patch_2D%edges%cell_idx(je,jb,2)
! !           ib_c2 = patch_2D%edges%cell_blk(je,jb,2)
! ! 
! !           ptp_vn(je,jb) = &
! !             &  DOT_PRODUCT(p_vn_c(il_c1,ib_c1)%x, operators_coefficients%edge2cell_coeff_cc_t(je,1,jb,1)%x)&
! !             & +DOT_PRODUCT(p_vn_c(il_c2,ib_c2)%x, operators_coefficients%edge2cell_coeff_cc_t(je,1,jb,2)%x)
! !         ELSE
! !           ptp_vn(je,jb) = 0.0_wp
! !         ENDIF
! ! 
! !       END DO edge_idx_loop
! !     END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block
! ! 
! ! !     CALL sync_patch_array(SYNC_E, patch_2D, ptp_vn)
! !     
! !   END SUBROUTINE map_cell2edges_2d
  !-----------------------------------------------------------------------------


 
END MODULE mo_scalar_product

