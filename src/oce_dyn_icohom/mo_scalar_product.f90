!>
!!               The module <i>mo_scalar_product</i>.
!!
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
!!  methods used are mpi parallelized, LL
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_scalar_product
  !-------------------------------------------------------------------------
  !
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma, p_test_run
  !USE mo_run_config,         ONLY: dtime
  USE mo_intp_data_strc,     ONLY: p_int_state
  USE mo_impl_constants,     ONLY: sea_boundary, boundary, sea,   &!max_char_length, &
  !  &                             land, land_boundary, boundary  &
    & min_rlcell, min_rledge , min_rlvert
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_model_domain,       ONLY: t_patch
  USE mo_oce_state,          ONLY: t_hydro_ocean_diag, v_base!, t_hydro_ocean_state
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce !, ab_gam
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, gc2cc, vector_product,&
    & gvec2cvec, cvec2gvec
  !USE mo_oce_index,          ONLY: ne_b, ne_i, nv_b, nv_i, form4ar, ldbg, c_k!, c_b, c_i
  USE mo_exception,                 ONLY: message, finish
  USE mo_physical_constants, ONLY: earth_radius
  USE mo_math_constants,     ONLY: pi
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_oce_math_operators,  ONLY: rot_vertex_ocean_3d, map_edges2vert_3d
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, &
    & sync_idx, global_max
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  PUBLIC :: calc_scalar_product_veloc_3d
  PUBLIC :: nonlinear_coriolis_3d
  PUBLIC :: map_edges2vert_3d
  PUBLIC :: map_edges2cell_3d
  PRIVATE :: map_edges2cell_with_height_3d
  PRIVATE :: map_edges2cell_no_height_3d
  
  PUBLIC :: map_cell2edges
  PUBLIC :: map_cell2edges_2d
  PUBLIC :: map_edges2cell
  
  PUBLIC :: primal_map_c2e
  PUBLIC :: dual_flip_flop
  
  PUBLIC :: map_edges2edges
  
  
  INTERFACE map_cell2edges
    MODULE PROCEDURE map_cell2edges_1level
    MODULE PROCEDURE map_cell2edges_mlevels
  END INTERFACE
    
  
  INTERFACE map_edges2cell
    
    MODULE PROCEDURE map_edges2cell_with_height
    MODULE PROCEDURE map_edges2cell_without_height
    MODULE PROCEDURE map_edges2cell_with_h_1vl      ! 1 vertical level
    MODULE PROCEDURE map_edges2cell_without_h_1vl   ! 1 vertical level
    
  END INTERFACE
  
  
  INTERFACE map_edges2cell_3d
    
    MODULE PROCEDURE map_edges2cell_with_height_3d
    MODULE PROCEDURE map_edges2cell_no_height_3d
    MODULE PROCEDURE map_edges2cell_no_h_3d_1vl
    
  END INTERFACE
  
  INTERFACE map_edges2edges
    
    MODULE PROCEDURE map_edges2edges_with_height
    MODULE PROCEDURE map_edges2edges_without_height
    MODULE PROCEDURE map_edges2edges_with_h_1vl      ! 1 vertical level
    MODULE PROCEDURE map_edges2edges_without_h_1vl   ! 1 vertical level
    
  END INTERFACE
  
CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
    !!  mpi parallelized by LL
  SUBROUTINE calc_scalar_product_veloc_3d( p_patch, vn_e_old, vn_e_new,&
    & h_e, p_diag, p_op_coeff)
    
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch            ! patch on which computation is performed
    REAL(wp), INTENT(in)      :: vn_e_old(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)      :: vn_e_new(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)      :: h_e(:,:)           ! SW-case: h_e is thicknerss at edges ! 3D case: h_e is surface elevation at edges
    TYPE(t_hydro_ocean_diag)  :: p_diag
    TYPE(t_operator_coeff)    :: p_op_coeff
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: jc, jb, jk

    TYPE(t_subset_range), POINTER :: all_cells
    
    !TYPE(t_cartesian_coordinates)    :: z_pv_cc(nproma,n_zlev,p_patch%nblks_c)
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    all_cells => p_patch%cells%all
    
    slev = 1
    elev = n_zlev
        
    CALL map_edges2vert_3d(p_patch, vn_e_old, h_e, p_op_coeff%edge2vert_coeff_cc, &
      & p_diag%p_vn_dual)
    !CALL map_edges2vert(p_patch, vn_e_old, h_e, p_diag%p_vn_dual)
    !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
    
    !CALL map_edges2cell(p_patch, vn_e_old, p_diag%p_vn)
    CALL map_edges2cell_3d(p_patch, vn_e_old, p_op_coeff, p_diag%p_vn)
    !CALL map_edges2cell_3D(p_patch, vn_e_old, p_diag%p_vn, p_op_coeff)
        
    !--------------------------------------------------------------    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      !calculate kinetic energy
      DO jk = slev, elev
        DO jc =  i_startidx_c, i_endidx_c
          
          IF ( v_base%lsm_oce_c(jc,jk,jb) > sea_boundary ) THEN
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
            & p_patch%cells%center(jc,jb)%lon,&
            & p_patch%cells%center(jc,jb)%lat,&
            & p_diag%u(jc,jk,jb), p_diag%v(jc,jk,jb) )
        END DO
      END DO
    END DO
    ! LL: no sync is required
    !--------------------------------------------------------------    
    
    CALL map_cell2edges( p_patch, p_diag%p_vn, p_diag%ptp_vn)

    CALL sync_patch_array(SYNC_E, p_patch, p_diag%ptp_vn)    
    
  END SUBROUTINE calc_scalar_product_veloc_3d
  !-------------------------------------------------------------------------
  
 
  !-------------------------------------------------------------------------
  !>
  !!
  !! mpi parallelized by LL, openmp corrected
  SUBROUTINE nonlinear_coriolis_3d(p_patch, vn, p_vn, p_vn_dual, h_e, vort_v, &
    & p_op_coeff, vort_flux)
    
    TYPE(t_patch), TARGET, INTENT(in)         :: p_patch
    REAL(wp), INTENT(inout)                   :: vn(:,:,:)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn(nproma,n_zlev,p_patch%nblks_c)
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_dual(nproma,n_zlev,p_patch%nblks_v)
    REAL(wp), INTENT(in)                      :: h_e(:,:)
    REAL(wp), INTENT(inout)                   :: vort_v(:,:,:)
    TYPE(t_operator_coeff),INTENT(in)         :: p_op_coeff
    REAL(wp), INTENT(inout)                   :: vort_flux(:,:,:)
    
    !Local variables
    !
    REAL(wp) :: z_vort_tmp, z_vort_tmp_boundary
    !REAL(wp) :: z_weight(nproma,n_zlev,p_patch%nblks_v)
    REAL(wp) :: zarea_fraction
    REAL(wp) :: z_area_scaled
    
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jv, jk, jb, jev,je
    INTEGER :: ile, ibe!, il, ib
    INTEGER :: i_startidx_e, i_endidx_e
    
    INTEGER :: icell_idx_1, icell_blk_1
    INTEGER :: icell_idx_2, icell_blk_2
    INTEGER :: il_v1, il_v2, ib_v1, ib_v2
    TYPE(t_cartesian_coordinates) :: u_v1_cc, u_v2_cc

    TYPE(t_subset_range), POINTER :: all_edges            
    !-----------------------------------------------------------------------
    all_edges => p_patch%edges%all
    
    slev         = 1
    elev         = n_zlev
    
    CALL rot_vertex_ocean_3d( p_patch, vn, p_vn_dual, p_op_coeff, vort_v)
    CALL sync_patch_array(SYNC_V, p_patch, vort_v)
    
    
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx_e,i_endidx_e,il_v1,ib_v1,il_v2,ib_v2, u_v1_cc, u_v2_cc)
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
      
      level_loop: DO jk = slev, elev
        
        edge_idx_loop: DO je =  i_startidx_e, i_endidx_e
          
          IF (v_base%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN
            !Get indices of two adjacent vertices
            il_v1 = p_patch%edges%vertex_idx(je,jb,1)
            ib_v1 = p_patch%edges%vertex_blk(je,jb,1)
            il_v2 = p_patch%edges%vertex_idx(je,jb,2)
            ib_v2 = p_patch%edges%vertex_blk(je,jb,2)
            
            !Multiply velocity reconstruction at vertex by vorticity
            u_v1_cc%x = p_vn_dual(il_v1,jk,ib_v1)%x &
              & *(vort_v(il_v1,jk,ib_v1) + p_patch%verts%f_v(il_v1,ib_v1))
            u_v2_cc%x = p_vn_dual(il_v2,jk,ib_v2)%x &
              & *(vort_v(il_v2,jk,ib_v2) + p_patch%verts%f_v(il_v2,ib_v2))
            
            !calculate finall vortex flux by mapping the vorticity-velocity-product
            !from vertices back to edges
            vort_flux(je,jk,jb) = &
              & - DOT_PRODUCT(u_v2_cc%x, p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
              & + DOT_PRODUCT(u_v1_cc%x, p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            
            !          IF(   i_v_ctr(il_v1,jk,ib_v1)==p_patch%verts%num_edges(il_v1,ib_v1)&
            !          &.AND.i_v_ctr(il_v2,jk,ib_v2)==p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !            vort_flux(je,jk,jb) = &
            !            &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !            &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !           ELSEIF( i_v_ctr(il_v1,jk,ib_v1)==p_patch%verts%num_edges(il_v1,ib_v1)&
            !           &.AND.  i_v_ctr(il_v2,jk,ib_v2)<p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !             IF(i_v_ctr(il_v2,jk,ib_v2)>=p_patch%verts%num_edges(il_v2,ib_v2)-3)THEN
            !              vort_flux(je,jk,jb) = &
            !               &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !               &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !             ELSE
            !               vort_flux(je,jk,jb) = &
            !               & 2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !             ENDIF
            !           ELSEIF( i_v_ctr(il_v1,jk,ib_v1)<p_patch%verts%num_edges(il_v1,ib_v1)&
            !             &.AND.i_v_ctr(il_v2,jk,ib_v2)==p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !             IF(i_v_ctr(il_v1,jk,ib_v1)>=p_patch%verts%num_edges(il_v1,ib_v1)-3)THEN
            !               vort_flux(je,jk,jb) = &
            !               &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !               &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !             ELSE
            !               vort_flux(je,jk,jb) = &
            !               &- 2.0_wp*DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)
            !             ENDIF
            !           ELSE
            !             vort_flux(je,jk,jb) = &
            !             &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !             &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !         ENDIF
            
          ELSE
            vort_flux(je,jk,jb)= 0.0_wp
          ENDIF ! (v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)
          
        END DO edge_idx_loop
      END DO level_loop
    END DO ! jb = all_edges%start_block, all_edges%end_block
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ! LL: no sync required
    
  END SUBROUTINE nonlinear_coriolis_3d
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  !!  mpi note: in and out is not synced. Should be done in the calling method if required
  SUBROUTINE map_edges2edges_with_height( p_patch, vn_in, vn_out,&
    & h_e, opt_slev,opt_elev )
    
    TYPE(t_patch), INTENT(in) :: p_patch            ! patch on which computation is performed
    REAL(wp), INTENT(inout)   :: vn_in(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(inout)   :: vn_out(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)      :: h_e(:,:)           ! SW-case: h_e is thicknerss at edges ! 3D case: h_e is surface elevation at edges
    INTEGER, INTENT(in), OPTIONAL ::  opt_slev  ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL ::  opt_elev  ! optional vertical end level
    
    !Local variables
    TYPE(t_cartesian_coordinates)    :: z_p_vn_cc(nproma,n_zlev,p_patch%nblks_c)
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    
    CALL map_edges2cell( p_patch, vn_in, z_p_vn_cc, h_e)
    CALL map_cell2edges( p_patch, z_p_vn_cc, vn_out )
    
  END SUBROUTINE map_edges2edges_with_height
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  !!  mpi note: in and out is not synced. Should be done in the calling method if required
  SUBROUTINE map_edges2edges_with_h_1vl( p_patch, vn_in, vn_out,&
    & h_e, level  )
    
    TYPE(t_patch), INTENT(in) :: p_patch            ! patch on which computation is performed
    REAL(wp), INTENT(inout)   :: vn_in(:,:)    ! input vector (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(inout)   :: vn_out(:,:)    ! input vector (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(in)      :: h_e(:,:)           ! SW-case: h_e is thicknerss at edges ! 3D case: h_e is surface elevation at edges
    !TYPE(t_operator_coeff)    :: p_op_coeff
    INTEGER, INTENT(in)       ::  level  ! vertical  level
    
    !Local variables
    !   INTEGER, PARAMETER :: no_cell_edges = 3
    !   INTEGER :: slev, elev
    !   INTEGER :: rl_start_c, rl_end_c, rl_start_e, rl_end_e
    !   INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
    !   INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    !   INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    !   INTEGER :: il_e, ib_e
    !   INTEGER :: jc, jb, jk, ie,je
    TYPE(t_cartesian_coordinates)    :: z_p_vn_cc(nproma,p_patch%nblks_c)
    !   REAL(wp) :: vort_flux_e(nproma,1,p_patch%nblks_e)
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    
    CALL map_edges2cell( p_patch, vn_in, z_p_vn_cc, h_e, level)
    CALL map_cell2edges( p_patch, z_p_vn_cc, vn_out, level )
    
  END SUBROUTINE map_edges2edges_with_h_1vl
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  !!  mpi note: in and out is not synced. Should be done in the calling method if required
  SUBROUTINE map_edges2edges_without_height( p_patch, vn_in, vn_out, opt_slev, opt_elev)
    
    TYPE(t_patch), INTENT(in) :: p_patch            ! patch on which computation is performed
    REAL(wp), INTENT(inout)   :: vn_in(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(inout)   :: vn_out(:,:,:)    ! input vector (nproma,n_zlev,nblks_e)
    !TYPE(t_operator_coeff)    :: p_op_coeff
    INTEGER, INTENT(in), OPTIONAL ::  opt_slev  ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL ::  opt_elev  ! optional vertical end level
    
    !Local variables
    !   INTEGER, PARAMETER :: no_cell_edges = 3
    !   INTEGER :: slev, elev
    !   INTEGER :: rl_start_c, rl_end_c, rl_start_e, rl_end_e
    !   INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
    !   INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    !   INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    !   INTEGER :: il_e, ib_e
    !   INTEGER :: jc, jb, jk, ie,je
    TYPE(t_cartesian_coordinates)    :: z_p_vn_cc(nproma,n_zlev,p_patch%nblks_c)
    !   REAL(wp) :: vort_flux_e(nproma,1,p_patch%nblks_e)
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
        
    CALL map_edges2cell( p_patch, vn_in, z_p_vn_cc )
    CALL map_cell2edges( p_patch, z_p_vn_cc, vn_out)
    
  END SUBROUTINE map_edges2edges_without_height
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!
  !!  mpi note: in and out is not synced. Should be done in the calling method if required
  SUBROUTINE map_edges2edges_without_h_1vl( p_patch, vn_in, vn_out, level)
    
    TYPE(t_patch), INTENT(in) :: p_patch            ! patch on which computation is performed
    REAL(wp), INTENT(inout)   :: vn_in(:,:)    ! input vector (nproma,n_zlev,nblks_e)
    REAL(wp), INTENT(inout)   :: vn_out(:,:)    ! input vector (nproma,n_zlev,nblks_e)
    !TYPE(t_operator_coeff)    :: p_op_coeff
    INTEGER, INTENT(in)  ::  level  !  vertical  level
    
    !Local variables
    !   INTEGER, PARAMETER :: no_cell_edges = 3
    !   INTEGER :: slev, elev
    !   INTEGER :: rl_start_c, rl_end_c, rl_start_e, rl_end_e
    !   INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
    !   INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    !   INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    !   INTEGER :: il_e, ib_e
    !   INTEGER :: jc, jb, jk, ie,je
    TYPE(t_cartesian_coordinates)    :: z_p_vn_cc(nproma,p_patch%nblks_c)
    !   REAL(wp) :: vort_flux_e(nproma,1,p_patch%nblks_e)
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
        
    CALL map_edges2cell( p_patch, vn_in, z_p_vn_cc, level )
    CALL map_cell2edges( p_patch, z_p_vn_cc, vn_out, level)
    
  END SUBROUTINE map_edges2edges_without_h_1vl
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized by LL, result not synced
  SUBROUTINE map_cell2edges_mlevels( p_patch, p_vn_c, ptp_vn, &
                                   & opt_slev, opt_elev, subset_range )
    
    TYPE(t_patch), TARGET,  INTENT(in)        :: p_patch          ! patch on which computation is performed
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_c(:,:,:)    ! input vector (nproma,n_zlev,nblks_c)
    REAL(wp), INTENT(inout)                   :: ptp_vn(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
    !TYPE(t_operator_coeff)                    :: p_op_coeff
    INTEGER, INTENT(in), OPTIONAL :: opt_slev        ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_elev        ! optional vertical end level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    
    
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: i_startidx_e, i_endidx_e
    !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: je, jb, jk!, i, ie,je
    INTEGER :: il_c1,ib_c1, il_c2,ib_c2!, ile_c1 , ibe_c1, ile_c2 , ibe_c2
    
    TYPE(t_subset_range), POINTER :: edges_in_domain        

    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    
    edges_in_domain   => p_patch%edges%in_domain
    IF (p_test_run) ptp_vn(:,:,:) = 0.0_wp
    
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
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
          
      level_loop_e: DO jk = slev, elev
        edge_idx_loop: DO je =  i_startidx_e, i_endidx_e
          
          IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
            
            !Get indices of two adjacent triangles
            il_c1 = p_patch%edges%cell_idx(je,jb,1)
            ib_c1 = p_patch%edges%cell_blk(je,jb,1)
            il_c2 = p_patch%edges%cell_idx(je,jb,2)
            ib_c2 = p_patch%edges%cell_blk(je,jb,2)
            ptp_vn(je,jk,jb) =&
              & DOT_PRODUCT(p_vn_c(il_c1,jk,ib_c1)%x,&
              & p_int_state(1)%edge2cell_coeff_cc_t(je,jb,1)%x)&
              & +DOT_PRODUCT(p_vn_c(il_c2,jk,ib_c2)%x,&
              & p_int_state(1)%edge2cell_coeff_cc_t(je,jb,2)%x)
          ELSE
            ptp_vn(je,jk,jb) = 0.0_wp
          ENDIF
          
        END DO edge_idx_loop
      END DO level_loop_e
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

    ! sync the result if necessary
    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
       & CALL sync_patch_array(SYNC_E, p_patch, ptp_vn)
    ENDIF
       
    !stop
  END SUBROUTINE map_cell2edges_mlevels
  !-----------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized by LL, result not synced
  SUBROUTINE map_cell2edges_1level( p_patch, p_vn_c, ptp_vn, level, subset_range )
    
    TYPE(t_patch), TARGET,  INTENT(in)        :: p_patch          ! patch on which computation is performed
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_c(:,:)    ! input vector (nproma,n_zlev,nblks_c)
    REAL(wp), INTENT(out)                   :: ptp_vn(:,:)    ! output vector (nproma,n_zlev,nblks_e)
    !TYPE(t_operator_coeff)                    :: p_op_coeff
    INTEGER, INTENT(in) :: level          ! vertical level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    
    !Local variables
    INTEGER :: i_startidx_e, i_endidx_e
    !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: je, jb !, i, ie,je
    INTEGER :: il_c1,ib_c1, il_c2,ib_c2!, ile_c1 , ibe_c1, ile_c2 , ibe_c2
    
    TYPE(t_subset_range), POINTER :: edges_in_domain        

    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    
    edges_in_domain   => p_patch%edges%in_domain
    IF (p_test_run) ptp_vn(:,:) = 0.0_wp
    
    
    ! calculation of transposed P^TPv from Pv (incart coord)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)          
      edge_idx_loop: DO je =  i_startidx_e, i_endidx_e

        IF(v_base%lsm_oce_e(je,level,jb) <= sea_boundary)THEN

          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(je,jb,1)
          ib_c1 = p_patch%edges%cell_blk(je,jb,1)
          il_c2 = p_patch%edges%cell_idx(je,jb,2)
          ib_c2 = p_patch%edges%cell_blk(je,jb,2)
          ptp_vn(je,jb) =&
            & DOT_PRODUCT(p_vn_c(il_c1,ib_c1)%x,&
            & p_int_state(1)%edge2cell_coeff_cc_t(je,jb,1)%x)&
            & +DOT_PRODUCT(p_vn_c(il_c2,ib_c2)%x,&
            & p_int_state(1)%edge2cell_coeff_cc_t(je,jb,2)%x)
        ELSE
          ptp_vn(je,jb) = 0.0_wp
        ENDIF

      END DO edge_idx_loop
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

    ! sync the result if necessary
    IF (PRESENT(subset_range)) THEN
      IF (.NOT. subset_range%is_in_domain) &
       & CALL sync_patch_array(SYNC_E, p_patch, ptp_vn)
    ENDIF
       
    !stop
  END SUBROUTINE map_cell2edges_1level
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized LL, result not synced
  SUBROUTINE map_cell2edges_2d( p_patch, p_vn_c, ptp_vn)
    
    TYPE(t_patch), TARGET, INTENT(in)         :: p_patch        ! patch on which computation is performed
    TYPE(t_cartesian_coordinates), INTENT(in) :: p_vn_c(:,:)    ! input vector (nproma,n_zlev,nblks_c)
    REAL(wp), INTENT(inout)                   :: ptp_vn(:,:)    ! output vector (nproma,n_zlev,nblks_e)
        
    !Local variables
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: je, jb

    TYPE(t_subset_range), POINTER :: edges_in_domain        
    
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    edges_in_domain   => p_patch%edges%in_domain
        
    ! calculation of transposed P^TPv from Pv (incart coord)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      
      edge_idx_loop: DO je =  i_startidx_e, i_endidx_e
        
        IF(v_base%lsm_oce_e(je,1,jb) < sea_boundary)THEN
          
          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(je,jb,1)
          ib_c1 = p_patch%edges%cell_blk(je,jb,1)
          il_c2 = p_patch%edges%cell_idx(je,jb,2)
          ib_c2 = p_patch%edges%cell_blk(je,jb,2)
          
          ptp_vn(je,jb) = &
            &  DOT_PRODUCT(p_vn_c(il_c1,ib_c1)%x, p_int_state(1)%edge2cell_coeff_cc_t(je,jb,1)%x)&
            & +DOT_PRODUCT(p_vn_c(il_c2,ib_c2)%x, p_int_state(1)%edge2cell_coeff_cc_t(je,jb,2)%x)
        ELSE
          ptp_vn(je,jb) = 0.0_wp
        ENDIF
        
      END DO edge_idx_loop
    END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block
    
!     CALL sync_patch_array(SYNC_E, p_patch, ptp_vn)
    
  END SUBROUTINE map_cell2edges_2d
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !>
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized LL
  SUBROUTINE map_edges2cell_with_height( p_patch, vn_e, p_vn_c, h_e, opt_slev, opt_elev, &
    & subset_range)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: p_patch        ! patch on which computation is performed
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    ! 3D case: h_e is surface elevation at edges
    TYPE(t_cartesian_coordinates),INTENT(inout):: p_vn_c(:,:,:)  ! outputput (nproma,n_zlev,nblks_c)
    REAL(wp), INTENT(in)                       :: h_e(:,:)       ! SW-case: h_e is thickness at edges
    !TYPE(t_operator_coeff)                     :: p_op_coeff
    INTEGER, INTENT(in), OPTIONAL :: opt_slev       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_elev       ! optional vertical end level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    
    !Local variables
    INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    !INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: il_e, ib_e
    INTEGER :: jc, jb, jk, ie!,je
    REAL(wp) :: z_weight
    REAL(wp) :: z_thick_e

    TYPE(t_subset_range), POINTER :: all_cells

    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => p_patch%cells%all
    ENDIF
    
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
            DO ie=1, no_cell_edges
              
              il_e = p_patch%cells%edge_idx(jc,jb,ie)
              ib_e = p_patch%cells%edge_blk(jc,jb,ie)
              !IF(v_base%lsm_oce_e(il_e,jk,ib_e) <= sea_boundary)THEN
              z_thick_e = v_base%del_zlev_m(slev) + h_e(il_e,ib_e)
              z_weight = z_weight + p_int_state(1)%variable_vol_norm(jc,jb,ie) * z_thick_e
              !ENDIF
              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x             &
                & + p_int_state(1)%edge2cell_coeff_cc(jc,jb,ie)%x &
                & * vn_e(il_e,jk,ib_e)* z_thick_e
            END DO
            IF ( z_weight/=0.0_wp) THEN!v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x / z_weight
            ELSE
              p_vn_c(jc,jk,jb)%x=0.0_wp
            ENDIF
          END DO cell_idx_loop_swm
        END DO level_loop_swm
      END DO  ! jb = all_cells%start_block, all_cells%end_block
      
    ELSEIF( iswm_oce /= 1 ) THEN
      
      !Step 1: Calculation of Pv in cartesian coordinates
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !We are dealing with the surface layer first
        cell_idx_loop_top: DO jc =  i_startidx_c, i_endidx_c
          z_weight             = 0.0_wp
          p_vn_c(jc,slev,jb)%x = 0.0_wp
          DO ie=1, no_cell_edges
            
            il_e = p_patch%cells%edge_idx(jc,jb,ie)
            ib_e = p_patch%cells%edge_blk(jc,jb,ie)
            !IF(v_base%lsm_oce_e(il_e,jk,ib_e) <= sea_boundary)THEN
            z_thick_e = v_base%del_zlev_m(slev) + h_e(il_e,ib_e)
            z_weight = z_weight + p_int_state(1)%variable_vol_norm(jc,jb,ie) * z_thick_e
            !ENDIF
            p_vn_c(jc,slev,jb)%x = p_vn_c(jc,slev,jb)%x&
              & + p_int_state(1)%edge2cell_coeff_cc(jc,jb,ie)%x&
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
            z_weight = 0.0_wp
            DO ie = 1, no_cell_edges
              
              il_e = p_patch%cells%edge_idx(jc,jb,ie)
              ib_e = p_patch%cells%edge_blk(jc,jb,ie)
              
              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x             &
                & + p_int_state(1)%edge2cell_coeff_cc(jc,jb,ie)%x &
                & * vn_e(il_e,jk,ib_e)
            END DO
            
            p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x/p_int_state(1)%fixed_vol_norm(jc,jb)!z_weight!
            
          END DO cell_idx_loop
        END DO level_loop
      END DO ! jb = all_cells%start_block, all_cells%end_block
    ENDIF
    ! LL no sync required
    
  END SUBROUTINE map_edges2cell_with_height
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !>
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized LL
  SUBROUTINE map_edges2cell_with_h_1vl( p_patch, vn_e, p_vn_c, h_e, level, &
    & subset_range)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: p_patch        ! patch on which computation is performed
    REAL(wp), INTENT(in)                       :: vn_e(:,:)    ! input (nproma,n_zlev,nblks_e)
    ! 3D case: h_e is surface elevation at edges
    TYPE(t_cartesian_coordinates),INTENT(inout):: p_vn_c(:,:)  ! outputput (nproma,n_zlev,nblks_c)
    REAL(wp), INTENT(in)                       :: h_e(:,:)       ! SW-case: h_e is thickness at edges
    !TYPE(t_operator_coeff)                     :: p_op_coeff
    INTEGER, INTENT(in)                        :: level       !  vertical level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    
    !Local variables
    INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    !INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: il_e, ib_e
    INTEGER :: jc, jb, ie!,je
    REAL(wp) :: z_weight
    REAL(wp) :: z_thick_e

    TYPE(t_subset_range), POINTER :: all_cells

    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => p_patch%cells%all
    ENDIF
    
        
    IF ( iswm_oce == 1 ) THEN
      
      !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
        cell_idx_loop_swm: DO jc =  i_startidx_c, i_endidx_c
          !calculate velocity reconstruction at cell center
          z_weight           = 0.0_wp
          p_vn_c(jc,jb)%x = 0.0_wp

          DO ie=1, no_cell_edges

            il_e = p_patch%cells%edge_idx(jc,jb,ie)
            ib_e = p_patch%cells%edge_blk(jc,jb,ie)
            !IF(v_base%lsm_oce_e(il_e,jk,ib_e) <= sea_boundary)THEN
            z_thick_e = v_base%del_zlev_m(level) + h_e(il_e,ib_e)
            z_weight = z_weight + p_int_state(1)%variable_vol_norm(jc,jb,ie) * z_thick_e
            !ENDIF
            p_vn_c(jc,jb)%x = p_vn_c(jc,jb)%x             &
              & + p_int_state(1)%edge2cell_coeff_cc(jc,jb,ie)%x &
              & * vn_e(il_e,ib_e)* z_thick_e

          END DO ! ie=1, no_cell_edges

          IF ( z_weight/=0.0_wp) THEN!v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary)THEN
            p_vn_c(jc,jb)%x = p_vn_c(jc,jb)%x / z_weight
          ELSE
            p_vn_c(jc,jb)%x=0.0_wp
          ENDIF

        END DO cell_idx_loop_swm
      END DO  ! jb = all_cells%start_block, all_cells%end_block
      
    ELSEIF( iswm_oce /= 1 ) THEN
      
      !Step 1: Calculation of Pv in cartesian coordinates
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !We are dealing with the surface layer first
        cell_idx_loop_top: DO jc =  i_startidx_c, i_endidx_c
          z_weight             = 0.0_wp
          p_vn_c(jc,jb)%x = 0.0_wp
          DO ie=1, no_cell_edges
            
            il_e = p_patch%cells%edge_idx(jc,jb,ie)
            ib_e = p_patch%cells%edge_blk(jc,jb,ie)
            !IF(v_base%lsm_oce_e(il_e,jk,ib_e) <= sea_boundary)THEN
            z_thick_e = v_base%del_zlev_m(level) + h_e(il_e,ib_e)
            z_weight = z_weight + p_int_state(1)%variable_vol_norm(jc,jb,ie) * z_thick_e
            !ENDIF
            p_vn_c(jc,jb)%x = p_vn_c(jc,jb)%x&
              & + p_int_state(1)%edge2cell_coeff_cc(jc,jb,ie)%x&
              & * vn_e(il_e,ib_e) * z_thick_e
            
          END DO
          
          IF(z_weight/=0.0_wp)THEN
            p_vn_c(jc,jb)%x = p_vn_c(jc,jb)%x / z_weight
          ELSE
            p_vn_c(jc,jb)%x=0.0_wp
          ENDIF

          
        END DO cell_idx_loop_top       
      END DO ! jb = all_cells%start_block, all_cells%end_block
    ENDIF
    ! LL no sync required
    
  END SUBROUTINE map_edges2cell_with_h_1vl
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized LL
  SUBROUTINE map_edges2cell_with_height_3d( p_patch, vn_e, p_op_coeff, p_vn_c, h_e,&
    & opt_slev, opt_elev, subset_range)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: p_patch        ! patch on which computation is performed
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    ! 3D case: h_e is surface elevation at edges
    TYPE(t_cartesian_coordinates),INTENT(inout):: p_vn_c(:,:,:)  ! outputput (nproma,n_zlev,nblks_c)
    REAL(wp), INTENT(in)                       :: h_e(:,:)       ! SW-case: h_e is thickness at edges
    TYPE(t_operator_coeff)                     :: p_op_coeff
    INTEGER, INTENT(in), OPTIONAL :: opt_slev       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_elev       ! optional vertical end level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    
    !Local variables
    INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    !INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: il_e, ib_e
    INTEGER :: jc, jb, jk, ie!,je
    REAL(wp) :: z_weight
    REAL(wp) :: z_thick_e

    TYPE(t_subset_range), POINTER :: all_cells
    
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => p_patch%cells%all
    ENDIF
    
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
            DO ie=1, no_cell_edges
              
              il_e = p_patch%cells%edge_idx(jc,jb,ie)
              ib_e = p_patch%cells%edge_blk(jc,jb,ie)
              
              z_thick_e = v_base%del_zlev_m(slev) + h_e(il_e,ib_e)
              z_weight = z_weight + p_int_state(1)%variable_vol_norm(jc,jb,ie) * z_thick_e
              
              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
                & + p_op_coeff%edge2cell_coeff_cc_dyn(jc,jk,jb,ie)%x&
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
          DO ie=1, no_cell_edges
            
            il_e = p_patch%cells%edge_idx(jc,jb,ie)
            ib_e = p_patch%cells%edge_blk(jc,jb,ie)
            
            z_thick_e = v_base%del_zlev_m(slev) + h_e(il_e,ib_e)
            z_weight = z_weight + p_int_state(1)%variable_vol_norm(jc,jb,ie) * z_thick_e
            
            p_vn_c(jc,slev,jb)%x = p_vn_c(jc,slev,jb)%x&
              & + p_op_coeff%edge2cell_coeff_cc_dyn(jc,1,jb,ie)%x&
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
            z_weight = 0.0_wp
            DO ie=1, no_cell_edges
              
              il_e = p_patch%cells%edge_idx(jc,jb,ie)
              ib_e = p_patch%cells%edge_blk(jc,jb,ie)
              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
                & + p_op_coeff%edge2cell_coeff_cc(jc,jk,jb,ie)%x&
                & * vn_e(il_e,jk,ib_e)
              !         IF(jk==2)THEN
              !         IF(p_vn_c(jc,jk,jb)%x(1)/=0.0_wp)&
              !          &write(54321,*)'p_vn 3D',jc,jk,jb,ie, p_vn_c(jc,jk,jb)%x, vn_e(il_e,jk,ib_e),p_op_coeff%edge2cell_coeff_cc(jc,jk,jb,ie)%x
              !         ENDIF
            END DO
          END DO cell_idx_loop
        END DO level_loop
      END DO ! jb = all_cells%start_block, all_cells%end_block
    ENDIF

    ! LL no sync required    
    
  END SUBROUTINE map_edges2cell_with_height_3d
  !----------------------------------------------------------------
    
  !----------------------------------------------------------------
  !  map_edges2cell_without_height
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized LL
  SUBROUTINE map_edges2cell_without_height( p_patch, vn_e, p_vn_c,opt_slev, opt_elev, &
    & subset_range )
    
    TYPE(t_patch), TARGET, INTENT(in)          :: p_patch        ! patch on which computation is performed
    REAL(wp)                                   :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    ! 3D case: h_e is surface elevation at edges
    TYPE(t_cartesian_coordinates),INTENT(inout):: p_vn_c(:,:,:)  ! outputput (nproma,n_zlev,nblks_c)
    INTEGER, INTENT(in), OPTIONAL :: opt_slev       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_elev       ! optional vertical end level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    !TYPE(t_operator_coeff)    :: p_op_coeff
    !Local variables
    INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    !INTEGER :: i_startblk_e, i_endblk_e!, i_startidx_e, i_endidx_e
    !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: il_e, ib_e
    INTEGER :: jc, jb, jk, ie
    REAL(wp) :: z_weight

    TYPE(t_subset_range), POINTER :: all_cells
    
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => p_patch%cells%all
    ENDIF
    
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
          !z_weight           = 0.0_wp
          DO ie=1, no_cell_edges
            il_e = p_patch%cells%edge_idx(jc,jb,ie)
            ib_e = p_patch%cells%edge_blk(jc,jb,ie)
            
            p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
              & + p_int_state(1)%edge2cell_coeff_cc(jc,jb,ie)%x&
              & * vn_e(il_e,jk,ib_e)
          END DO
          !IF(z_weight/=0.0_wp)THEN
          p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x /p_int_state(1)%fixed_vol_norm(jc,jb)
          !ELSE
          !  p_vn_c(jc,jk,jb)%x=0.0_wp
          !ENDIF
        END DO cell_idx_loop
      END DO level_loop
    END DO ! jb = all_cells%start_block, all_cells%end_block
    ! LL no sync required
    
  END SUBROUTINE map_edges2cell_without_height
  !-----------------------------------------------------------------------
  
  !----------------------------------------------------------------
  !  map_edges2cell_without_height
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized LL
  SUBROUTINE map_edges2cell_without_h_1vl( p_patch, vn_e, p_vn_c, level, &
    & subset_range )
    
    TYPE(t_patch), TARGET, INTENT(in)          :: p_patch        ! patch on which computation is performed
    REAL(wp)                                   :: vn_e(:,:)    ! input (nproma,n_zlev,nblks_e)
    ! 3D case: h_e is surface elevation at edges
    TYPE(t_cartesian_coordinates),INTENT(inout):: p_vn_c(:,:)  ! outputput (nproma,n_zlev,nblks_c)
    INTEGER, INTENT(in)            :: level          ! vertical level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    !TYPE(t_operator_coeff)    :: p_op_coeff
    !Local variables
    INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    !INTEGER :: i_startblk_e, i_endblk_e!, i_startidx_e, i_endidx_e
    !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: il_e, ib_e
    INTEGER :: jc, jb, ie
    REAL(wp) :: z_weight

    TYPE(t_subset_range), POINTER :: all_cells
    
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => p_patch%cells%all
    ENDIF
        
    !Calculation of Pv in cartesian coordinates
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
        cell_idx_loop: DO jc =  i_startidx_c, i_endidx_c
          !calculate velocity reconstruction at cell center
          p_vn_c(jc,jb)%x = 0.0_wp
          !z_weight           = 0.0_wp
          DO ie=1, no_cell_edges
            il_e = p_patch%cells%edge_idx(jc,jb,ie)
            ib_e = p_patch%cells%edge_blk(jc,jb,ie)
            
            p_vn_c(jc,jb)%x = p_vn_c(jc,jb)%x&
              & + p_int_state(1)%edge2cell_coeff_cc(jc,jb,ie)%x&
              & * vn_e(il_e,ib_e)
          END DO
          !IF(z_weight/=0.0_wp)THEN
          p_vn_c(jc,jb)%x = p_vn_c(jc,jb)%x /p_int_state(1)%fixed_vol_norm(jc,jb)
          !ELSE
          !  p_vn_c(jc,jk,jb)%x=0.0_wp
          !ENDIF
          
        END DO cell_idx_loop
    END DO ! jb = all_cells%start_block, all_cells%end_block
    ! LL no sync required
    
  END SUBROUTINE map_edges2cell_without_h_1vl
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  ! map_edges2cell_without_height
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized LL
  SUBROUTINE map_edges2cell_no_height_3d( p_patch, vn_e, p_op_coeff, p_vn_c, opt_slev, opt_elev, &
    &                                     subset_range)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: p_patch        ! patch on which computation is performed
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: p_op_coeff
    TYPE(t_cartesian_coordinates)              :: p_vn_c(:,:,:)  ! output (nproma,n_zlev,nblks_c)
                                                                 ! intent(inout) for nag compiler
    INTEGER, INTENT(in), OPTIONAL :: opt_slev       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL :: opt_elev       ! optional vertical end level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    !Local variables
    INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: il_e, ib_e
    INTEGER :: jc, jb, jk, ie

     TYPE(t_subset_range), POINTER :: all_cells
   
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => p_patch%cells%all
    ENDIF
    
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
          
          DO ie=1, no_cell_edges
            il_e = p_patch%cells%edge_idx(jc,jb,ie)
            ib_e = p_patch%cells%edge_blk(jc,jb,ie)
            
            p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
              & + p_op_coeff%edge2cell_coeff_cc(jc,jk,jb,ie)%x&
              & * vn_e(il_e,jk,ib_e)
          END DO
          IF(p_op_coeff%fixed_vol_norm(jc,jk,jb)/=0.0_wp)THEN
            p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x/p_op_coeff%fixed_vol_norm(jc,jk,jb)
          ENDIF
        END DO cell_idx_loop
      END DO level_loop
      
    END DO ! jb = all_cells%start_block, all_cells%end_block
    ! LL no sync required
    
  END SUBROUTINE map_edges2cell_no_height_3d
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  ! map_edges2cell_without_height
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized LL
  SUBROUTINE map_edges2cell_no_h_3d_1vl( p_patch, vn_e,p_op_coeff, p_vn_c, level, &
    &  subset_range)
    
    TYPE(t_patch), TARGET, INTENT(in)          :: p_patch        ! patch on which computation is performed
    REAL(wp), INTENT(in)                       :: vn_e(:,:)    ! input (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: p_op_coeff
    TYPE(t_cartesian_coordinates),INTENT(out)  :: p_vn_c(:,:)  ! outputput (nproma,n_zlev,nblks_c)
    INTEGER, INTENT(in)           :: level       ! vertical level
    TYPE(t_subset_range), TARGET,  OPTIONAL :: subset_range
    !Local variables
    INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: il_e, ib_e
    INTEGER :: jc, jb, ie

     TYPE(t_subset_range), POINTER :: all_cells
   
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_e2c')
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => p_patch%cells%all
    ENDIF
    
    
    !Calculation of Pv in cartesian coordinates
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      cell_idx_loop: DO jc =  i_startidx_c, i_endidx_c
        !calculate velocity reconstruction at cell center
        p_vn_c(jc,jb)%x = 0.0_wp

        DO ie=1, no_cell_edges
          il_e = p_patch%cells%edge_idx(jc,jb,ie)
          ib_e = p_patch%cells%edge_blk(jc,jb,ie)

          p_vn_c(jc,jb)%x = p_vn_c(jc,jb)%x&
            & + p_op_coeff%edge2cell_coeff_cc(jc,level,jb,ie)%x&
            & * vn_e(il_e,ib_e)
        END DO
        
        IF(p_op_coeff%fixed_vol_norm(jc,level,jb)/=0.0_wp)THEN
          p_vn_c(jc,jb)%x = p_vn_c(jc,jb)%x/p_op_coeff%fixed_vol_norm(jc,level,jb)
        ENDIF
        
      END DO cell_idx_loop
      
    END DO ! jb = all_cells%start_block, all_cells%end_block
    ! LL no sync required
    
  END SUBROUTINE map_edges2cell_no_h_3d_1vl
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !>
  !! Discrete mapping of cell-based into edge-based vectors on the primal grid.
  !!
  !! Discrete mapping of cell-based into edge-based vectors on the primal grid.
  !! This mapping depends on the grid geometry only and is used to formulate
  !! the primitive equations in weak form.
  !! The coefficients are calculated in sbrt init_scalar_product and stored in
  !! the ocean part of the patch.
  !! Input lives on cells (2-dim, 3-dim), output lives on edges (3-dim)
  !!
  !! @par Revision History
  !!  developed by Stephan Lorenz, MPI-M (2010-11)
  !!  mpi parallelized LL
  !!
  SUBROUTINE primal_map_c2e( p_patch, u_c, v_c, vn_e, h_c, opt_slev, opt_elev )
    
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch ! patch on which computation is performed
    REAL(wp),    INTENT(in)   :: u_c(:,:,:)     ! zonal input component (nproma,n_zlev,nblks_c)
    REAL(wp),    INTENT(in)   :: v_c(:,:,:)     ! meridional input component (nproma,n_zlev,nblks_c)
    REAL(wp),    INTENT(out)  :: vn_e(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
    REAL(wp),INTENT(in), OPTIONAL :: h_c(:,:)    ! SWE-case: thickness at cell centers
    INTEGER, INTENT(in), OPTIONAL ::  opt_slev  ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL ::  opt_elev  ! optional vertical end level
    
    !Local variables
    TYPE(t_cartesian_coordinates)   :: zu_cc(nproma,p_patch%nblks_c)
    TYPE(t_cartesian_coordinates)   :: zu_cc2(nproma,p_patch%nblks_c)
    INTEGER :: slev, elev
    INTEGER :: i_startidx, i_endidx
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: jb, jc, je, jk, ie
    
    TYPE(t_subset_range), POINTER :: all_cells, edges_in_domain 
    
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:primal_map_c2e')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    all_cells       => p_patch%cells%all
    edges_in_domain => p_patch%edges%in_domain
    
    DO ie = 1,3
      zu_cc(:,:)%x(ie) = 0.0_wp
    END DO
    
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
    
    ! #slo# due to nag -nan compiler-option - remove for production
    vn_e(:,slev:elev,:) = 0.0_wp
    
    level_loop: DO jk = slev, elev
      
      ! cells      
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
        
        cell_idx_loop: DO jc =  i_startidx, i_endidx
          
          CALL gvec2cvec( u_c(jc,jk,jb), v_c(jc,jk,jb),   &
            & p_patch%cells%center(jc,jb)%lon, &
            & p_patch%cells%center(jc,jb)%lat, &
            & zu_cc(jc,jb)%x(1),zu_cc(jc,jb)%x(2), zu_cc(jc,jb)%x(3))
          
          ! IF(jk==1)THEN
          ! IF(jc==1.and.jb<=10)THEN
          !  write(*,*)'P-cart-old:',jc,jb, zu_cc(jc,jb)%x
          ! ENDIF
          ! ENDIF
        END DO cell_idx_loop
      END DO  ! jb = all_cells%start_block, all_cells%end_block
      
      ! edges
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
        
        edge_idx_loop: DO je =  i_startidx, i_endidx
          
          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(je,jb,1)
          ib_c1 = p_patch%edges%cell_blk(je,jb,1)
          il_c2 = p_patch%edges%cell_idx(je,jb,2)
          ib_c2 = p_patch%edges%cell_blk(je,jb,2)
          
          IF (jk == 1) THEN
            IF ( PRESENT(h_c) ) THEN
              !IF(jb==17)then
              !write(*,*)'c2e:vec',jb,je,zu_cc(il_c1,ib_c1)%x, zu_cc(il_c2,ib_c2)%x,&
              !&h_c(il_c1,ib_c1),h_c(il_c2,ib_c2)
              !ENDIF
              !Multiply cartesian vector with height vaiable
              zu_cc2(il_c1,ib_c1)%x = zu_cc(il_c1,ib_c1)%x * h_c(il_c1,ib_c1)
              zu_cc2(il_c2,ib_c2)%x = zu_cc(il_c2,ib_c2)%x * h_c(il_c2,ib_c2)
              
              IF (v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
                vn_e(je,jk,jb) = &
                  &   DOT_PRODUCT(zu_cc2(il_c1,ib_c1)%x, &
                  &               p_int_state(1)%edge2cell_coeff_cc_t(je,jb,1)%x) &
                  & + DOT_PRODUCT(zu_cc2(il_c2,ib_c2)%x, &
                  &               p_int_state(1)%edge2cell_coeff_cc_t(je,jb,2)%x)
              ELSE
                vn_e(je,jk,jb) = 0.0_wp
              ENDIF
              
              ! !IF(je=12.and.jb==17)then
              ! IF(jb==17)then
              ! write(*,*)'c2e:vec*h',jb,je,zu_cc2(il_c1,ib_c1)%x, zu_cc2(il_c2,ib_c2)%x
              ! write(*,*)'c2e:',je,jb,vn_e(je,jk,jb),&
              ! &DOT_PRODUCT(zu_cc2(il_c1,ib_c1)%x, p_int_state(1)%edge2cell_coeff_t(je,jb,1)%x),&
              ! &DOT_PRODUCT(zu_cc2(il_c2,ib_c2)%x, p_int_state(1)%edge2cell_coeff_t(je,jb,2)%x)
              ! ENDIF
            ELSEIF (.NOT.PRESENT(h_c)) THEN
              
              IF (v_base%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN
                vn_e(je,jk,jb) = &
                  & DOT_PRODUCT(zu_cc(il_c1,ib_c1)%x,&
                  & p_int_state(1)%edge2cell_coeff_cc_t(je,jb,1)%x)&
                  & + DOT_PRODUCT(zu_cc(il_c2,ib_c2)%x,&
                  & p_int_state(1)%edge2cell_coeff_cc_t(je,jb,2)%x)
              ELSE
                vn_e(je,jk,jb) = 0.0_wp
              ENDIF
              ! IF(jb==17)then
              ! write(*,*)'c2e:vec',je,jb,zu_cc(il_c1,ib_c1)%x, zu_cc(il_c2,ib_c2)%x
              !  write(*,*)'c2e:',je,jb,il_c1,ib_c1,il_c2,ib_c2, vn_e(je,jk,jb),&
              !  &DOT_PRODUCT(zu_cc(il_c1,ib_c1)%x, p_int_state(1)%edge2cell_coeff_t(je,jb,1)%x),&
              !  &DOT_PRODUCT(zu_cc(il_c2,ib_c2)%x, p_int_state(1)%edge2cell_coeff_t(je,jb,2)%x)
              ! ENDIF
            ENDIF
          ELSEIF ( jk > 1 ) THEN
            IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
              vn_e(je,jk,jb) =&
                &   DOT_PRODUCT(zu_cc(il_c1,ib_c1)%x, &
                &               p_int_state(1)%edge2cell_coeff_cc_t(je,jb,1)%x) &
                & + DOT_PRODUCT(zu_cc(il_c2,ib_c2)%x, &
                &               p_int_state(1)%edge2cell_coeff_cc_t(je,jb,2)%x)
            ELSE
              vn_e(je,jk,jb) = 0.0_wp
            ENDIF
          ENDIF!jk-condition
          
        END DO edge_idx_loop
      END DO !jb = edges_in_domain%start_block, edges_in_domain%end_block
      
    END DO level_loop

   ! sync the result
    CALL sync_patch_array(SYNC_E, p_patch, vn_e(:,:,:))
  
  END SUBROUTINE primal_map_c2e
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! !   !>
  ! !   !! Discrete mapping of cell-based vectors to edges on the primal grid.
  ! !   !!
  ! !   !!
  ! !   !! @par Revision History
  ! !   !!  developed by Peter Korn, MPI-M (2010-11)
  ! !   !!
  ! !   SUBROUTINE map_cell2edges_upwind( p_patch, p_vec_c, p_vn_e, p_os, opt_slev, opt_elev )
  ! !
  ! !   TYPE(t_patch), INTENT(IN)                 :: p_patch          ! patch on which computation is performed
  ! !   TYPE(t_cartesian_coordinates), INTENT(IN) :: p_vec_c(:,:,:)    ! input vector (nproma,n_zlev,nblks_c)
  ! !   REAL(wp), INTENT(OUT)                     :: p_vn_e(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
  ! !   TYPE(t_hydro_ocean_state)                 :: p_os
  ! !   INTEGER, INTENT(IN), OPTIONAL             :: opt_slev        ! optional vertical start level
  ! !   INTEGER, INTENT(IN), OPTIONAL             :: opt_elev        ! optional vertical end level
  ! !
  ! !
  ! !   !Local variables
  ! !   INTEGER :: slev, elev
  ! !   INTEGER :: rl_start_e, rl_end_e
  ! !   INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  ! !   INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  ! !   INTEGER :: je, jb, jk!, ie,je
  ! !   REAL(wp) :: z_dot_c1, z_dot_c2, norm_c1_c2, norm
  ! !   TYPE(t_cartesian_coordinates)    :: cc_e0, cc_c1,cc_c2
  ! !   TYPE(t_cartesian_coordinates)    :: cv_c1_e0, cv_c2_e0, cv_c1_c2
  ! !   !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  ! !   !  & routine = ('mo_scalar_product:primal_map_e2c')
  ! !   !-----------------------------------------------------------------------
  ! !   !CALL message (TRIM(routine), 'start')
  ! !
  ! !   ! check optional arguments
  ! !   IF ( PRESENT(opt_slev) ) THEN
  ! !     slev = opt_slev
  ! !   ELSE
  ! !     slev = 1
  ! !   END IF
  ! !   IF ( PRESENT(opt_elev) ) THEN
  ! !     elev = opt_elev
  ! !   ELSE
  ! !     elev = n_zlev
  ! !   END IF
  ! !
  ! ! rl_start_e = 1
  ! ! rl_end_e  = min_rledge
  ! !
  ! ! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
  ! ! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
  ! !
  ! !
  ! ! ! calculation of transposed P^TPv from Pv (incart coord)
  ! ! EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e
  ! !
  ! !   CALL get_indices_e(p_patch, jb,&
  ! !                    & i_startblk_e, i_endblk_e,&
  ! !                    & i_startidx_e, i_endidx_e,&
  ! !                    & rl_start_e, rl_end_e)
  ! !
  ! !   LEVEL_LOOP_E: DO jk = slev, elev
  ! !     EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
  ! !
  ! !       IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
  ! !
  ! !         !Get indices of two adjacent triangles
  ! !         il_c1 = p_patch%edges%cell_idx(je,jb,1)
  ! !         ib_c1 = p_patch%edges%cell_blk(je,jb,1)
  ! !         il_c2 = p_patch%edges%cell_idx(je,jb,2)
  ! !         ib_c2 = p_patch%edges%cell_blk(je,jb,2)
  ! !
  ! !         cc_e0 = gc2cc(p_patch%edges%center(je,jb))
  ! !         cc_c1 = gc2cc(p_patch%cells%center(il_c1,ib_c1))
  ! !         cc_c2 = gc2cc(p_patch%cells%center(il_c2,ib_c2))
  ! !
  ! !         cv_c1_c2%x = cc_c1%x - cc_c2%x
  ! !         cv_c1_e0%x = cc_e0%x - cc_c1%x
  ! !         cv_c2_e0%x = cc_e0%x - cc_c2%x
  ! !
  ! !         norm_c1_c2 = SQRT(SUM(cv_c1_e0%x*cv_c1_e0%x))+SQRT(SUM(cv_c2_e0%x*cv_c2_e0%x))
  ! !
  ! !         IF(p_os%p_diag%ptp_vn(je,jk,jb)>=0.0_wp)THEN !velocity is pointing
  ! !           norm=SQRT(SUM(cv_c1_e0%x*cv_c1_e0%x))
  ! !           p_vn_e(je,jk,jb) =&
  ! !           & DOT_PRODUCT(p_vec_c(il_c1,jk,ib_c1)%x,&
  ! !           &(p_int_state(1)%edge2cell_coeff_cc_t(je,jb,1)%x)*(norm_c1_c2/norm))!&
  ! !           !&+&
  ! !           !& DOT_PRODUCT(p_vec_c(il_c2,jk,ib_c2)%x,&
  ! !           !&(p_int_state(1)%edge2cell_coeff_cc_t(je,jb,2)%x)*2.0_wp)!(norm_c1_c2/norm))
  ! !         !&+DOT_PRODUCT(p_vn_c(il_c1,jk,ib_c1)%x,p_int_state(1)%edge2cell_coeff_cc_t(je,jb,1)%x)
  ! !
  ! !         ELSEIF(p_os%p_diag%ptp_vn(je,jk,jb)<0.0_wp)THEN
  ! !           norm=SQRT(SUM(cv_c2_e0%x*cv_c2_e0%x))
  ! !           p_vn_e(je,jk,jb) =&
  ! !           & DOT_PRODUCT(p_vec_c(il_c2,jk,ib_c2)%x,&
  ! !           & (p_int_state(1)%edge2cell_coeff_cc_t(je,jb,2)%x)*(norm_c1_c2/norm))!&
  ! !           !&+&
  ! !           !& DOT_PRODUCT(p_vec_c(il_c1,jk,ib_c1)%x,&
  ! !           !&(p_int_state(1)%edge2cell_coeff_cc_t(je,jb,1)%x)*2.0_wp)!(norm_c1_c2/norm))
  ! !           !write(*,*)'frac',(norm_c1_c2/norm)
  ! !         ENDIF
  ! !
  ! !        ELSE
  ! !          p_vn_e(je,jk,jb) = 0.0_wp
  ! !        ENDIF
  ! !     END DO EDGE_IDX_LOOP
  ! !   END DO LEVEL_LOOP_E
  ! ! END DO EDGE_BLK_LOOP
  ! !
  ! !   END SUBROUTINE map_cell2edges_upwind
  ! ! !-----------------------------------------------------------------------------
  ! !   !>
  ! !   !! Discrete mapping of cell-based vectors to edges on the primal grid.
  ! !   !!
  ! !   !!
  ! !   !! @par Revision History
  ! !   !!  developed by Peter Korn, MPI-M (2010-11)
  ! !   !!
  ! !   SUBROUTINE map_cell2edges_upwind2( p_patch, p_vec_c, p_vn_e, p_os, opt_slev, opt_elev )
  ! !
  ! !   TYPE(t_patch), INTENT(IN)                 :: p_patch          ! patch on which computation is performed
  ! !   TYPE(t_cartesian_coordinates), INTENT(IN) :: p_vec_c(:,:,:)    ! input vector (nproma,n_zlev,nblks_c)
  ! !   REAL(wp), INTENT(OUT)                     :: p_vn_e(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
  ! !   TYPE(t_hydro_ocean_state)                 :: p_os
  ! !   INTEGER, INTENT(IN), OPTIONAL             :: opt_slev        ! optional vertical start level
  ! !   INTEGER, INTENT(IN), OPTIONAL             :: opt_elev        ! optional vertical end level
  ! !
  ! !
  ! !   !Local variables
  ! !   INTEGER :: slev, elev
  ! !   INTEGER :: rl_start_e, rl_end_e
  ! !   INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  ! !   INTEGER :: il_c1, ib_c1, il_c2, ib_c2
  ! !   INTEGER :: je, jb, jk!, ie,je
  ! !   REAL(wp) :: z_dot_c1, z_dot_c2, norm_c1_c2, norm
  ! !   !TYPE(t_cartesian_coordinates)    :: cc_e0, cc_c1,cc_c2
  ! !   TYPE(t_cartesian_coordinates)    :: cv_c1_e0, cv_c2_e0, cv_c1_c2
  ! !   !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  ! !   !  & routine = ('mo_scalar_product:primal_map_e2c')
  ! !   !-----------------------------------------------------------------------
  ! !   !CALL message (TRIM(routine), 'start')
  ! !
  ! !   ! check optional arguments
  ! !   IF ( PRESENT(opt_slev) ) THEN
  ! !     slev = opt_slev
  ! !   ELSE
  ! !     slev = 1
  ! !   END IF
  ! !   IF ( PRESENT(opt_elev) ) THEN
  ! !     elev = opt_elev
  ! !   ELSE
  ! !     elev = n_zlev
  ! !   END IF
  ! !
  ! ! rl_start_e = 1
  ! ! rl_end_e  = min_rledge
  ! !
  ! ! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
  ! ! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
  ! !
  ! !
  ! ! ! calculation of transposed P^TPv from Pv (incart coord)
  ! ! EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e
  ! !
  ! !   CALL get_indices_e(p_patch, jb,&
  ! !                    & i_startblk_e, i_endblk_e,&
  ! !                    & i_startidx_e, i_endidx_e,&
  ! !                    & rl_start_e, rl_end_e)
  ! !
  ! !   LEVEL_LOOP_E: DO jk = slev, elev
  ! !     EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
  ! !
  ! !       IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
  ! !
  ! !         !Get indices of two adjacent triangles
  ! !         il_c1 = p_patch%edges%cell_idx(je,jb,1)
  ! !         ib_c1 = p_patch%edges%cell_blk(je,jb,1)
  ! !         il_c2 = p_patch%edges%cell_idx(je,jb,2)
  ! !         ib_c2 = p_patch%edges%cell_blk(je,jb,2)
  ! !
  ! !           p_vn_e(je,jk,jb) =&
  ! !           & DOT_PRODUCT(p_vec_c(il_c1,jk,ib_c1)%x,&
  ! !           &(p_int_state(1)%edge2cell_coeff_cc_t(je,jb,1)%x))&
  ! !           &+&
  ! !           & DOT_PRODUCT(p_vec_c(il_c2,jk,ib_c2)%x,&
  ! !           &(p_int_state(1)%edge2cell_coeff_cc_t(je,jb,2)%x))&
  ! !           &-&
  ! !           p_os%p_diag%ptp_vn(je,jk,jb)*dtime*&
  ! !           &SQRT(SUM( (p_vec_c(il_c2,jk,ib_c2)%x-p_vec_c(il_c1,jk,ib_c1)%x)&
  ! !           &*p_patch%edges%inv_dual_edge_length(je,jb)&
  ! !           &*(p_vec_c(il_c2,jk,ib_c2)%x- p_vec_c(il_c1,jk,ib_c1)%x)&
  ! !           &*p_patch%edges%inv_dual_edge_length(je,jb)))
  ! !
  ! !        ELSE
  ! !          p_vn_e(je,jk,jb) = 0.0_wp
  ! !        ENDIF
  ! !     END DO EDGE_IDX_LOOP
  ! !   END DO LEVEL_LOOP_E
  ! ! END DO EDGE_BLK_LOOP
  ! !
  ! !   END SUBROUTINE map_cell2edges_upwind2
  !-----------------------------------------------------------------------------
  !   !
  !   !>
  !   !! Function implements discrete scalar product on the primal grid. This
  !   !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !   !! equations in weak form. The functions uses coefficients that are calculated in sbr
  !   !! "init_scalar_product".
  !   !!
  !   !! @par Revision History
  !   !!  developed by Peter Korn, MPI-M (2010-9)
  !   !!
  !   FUNCTION dual_flip_flop(p_patch, vn_old_e, vn_new_e, vort_v,&
  !                         & h_e, opt_slev, opt_elev) RESULT(vn_out_e)
  !
  !   TYPE(t_patch), INTENT(IN) :: p_patch
  !   REAL(wp), INTENT(inout)      :: vn_old_e(:,:,:)
  !   REAL(wp), INTENT(inout)      :: vn_new_e(:,:,:)
  !   REAL(wp), INTENT(in)      :: vort_v(:,:,:)
  !   REAL(wp), INTENT(IN)      :: h_e(:,:)
  !   INTEGER, INTENT(IN), OPTIONAL ::  opt_slev  ! optional vertical start level
  !   INTEGER, INTENT(IN), OPTIONAL ::  opt_elev  ! optional vertical end level
  !   REAL(wp)                  :: vn_out_e(SIZE(vn_old_e,1), SIZE(vn_old_e,2), SIZE(vn_old_e,3))
  !
  !   !Local variables
  !   INTEGER, PARAMETER :: no_vert_edges = 6
  !   INTEGER :: slev, elev
  !   INTEGER :: rl_start_v, rl_end_v
  !   INTEGER :: rl_start_e, rl_end_e
  !   INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v
  !   INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  !   INTEGER :: il_v1, il_v2, ib_v1, ib_v2, il_e, ib_e
  !   INTEGER :: je, jb, jv, ie, jk
  !
  !   REAL(wp) :: z_thick, z_weight
  !   REAL(wp) :: z_vn
  !   REAL(wp) :: z_vn_e(nproma,n_zlev,p_patch%nblks_e)
  !   TYPE(t_cartesian_coordinates) :: u_v_cc(SIZE(vort_v,1),SIZE(vort_v,3))
  !   !TYPE(t_cartesian_coordinates) ::  u_tmp(SIZE(vort_v,1),SIZE(vort_v,3))
  !   REAL(wp) :: u_v(SIZE(vort_v,1),SIZE(vort_v,3))
  !   REAL(wp) :: v_v(SIZE(vort_v,1),SIZE(vort_v,3))
  !   INTEGER  :: i_ctr
  !
  !   !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !   !  & routine = ('mo_scalar_product:dual_flip_flop')
  !   !-----------------------------------------------------------------------
  !   !CALL message (TRIM(routine), 'start')
  !    vn_out_e(:,:,:) = 0.0_wp
  !
  !   ! check optional arguments
  !   IF ( PRESENT(opt_slev) ) THEN
  !     slev = opt_slev
  !   ELSE
  !     slev = 1
  !   END IF
  !   IF ( PRESENT(opt_elev) ) THEN
  !     elev = opt_elev
  !   ELSE
  !     elev = n_zlev
  !   END IF
  !
  !  rl_start_v = 1
  !  rl_end_v = min_rlvert
  !  i_startblk_v = p_patch%verts%start_blk(rl_start_v,1)
  !  i_endblk_v   = p_patch%verts%end_blk(rl_end_v,1)
  !
  !  rl_start_e = 1
  !  rl_end_e = min_rledge
  !  i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
  !  i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
  !
  !  z_vn_e = ab_gam*vn_new_e + (1.0_wp-ab_gam)*vn_old_e
  !
  !   LEVEL_LOOP: DO jk = slev, elev
  !
  !   u_v(:,:) = 0.0_wp
  !   v_v(:,:) = 0.0_wp
  !
  !  !Step 1: edge to vertex mapping
  !
  ! !$OMP PARALLEL
  ! !$OMP DO PRIVATE(jb,jv,ie,i_startidx_v,i_endidx_v,il_e,ib_e)
  !     VERT_BLK_LOOP: DO jb = i_startblk_v, i_endblk_v
  !
  !       CALL get_indices_v(p_patch, jb, i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v, &
  !         &                rl_start_v, rl_end_v)
  !
  !       VERT_IDX_LOOP: DO jv = i_startidx_v, i_endidx_v
  !         u_v(jv,jb) = 0.0_wp
  !         v_v(jv,jb) = 0.0_wp
  !         u_v_cc(jv,jb)%x= 0.0_wp
  !         z_weight   = 0.0_wp
  !         i_ctr = 0
  !         DO ie=1, p_patch%verts%num_edges(jv,jb) !no_vert_edges
  !
  !           il_e = p_patch%verts%edge_idx(jv,jb,ie)
  !           ib_e = p_patch%verts%edge_blk(jv,jb,ie)
  !
  !           z_thick = v_base%del_zlev_m(jk)
  !            IF ( iswm_oce == 1 ) THEN
  !              z_thick = h_e(il_e,ib_e)
  !            ELSEIF( iswm_oce /= 1 ) THEN
  !              IF (jk == 1 )THEN
  !                z_thick = v_base%del_zlev_m(jk) + h_e(il_e,ib_e)
  !              ENDIF
  !            ENDIF
  !
  !            IF ( v_base%lsm_oce_e(il_e,jk,ib_e) <= sea_boundary ) THEN
  !              z_weight = z_weight + p_int_state(1)%variable_dual_vol_norm(jv,jb,ie)*z_thick
  !              u_v_cc(jv,jb)%x = u_v_cc(jv,jb)%x +                                  &
  !              &                 p_int_state(1)%edge2vert_coeff_cc(jv,jb,ie)%x * &
  !              &                       z_vn_e(il_e,jk,ib_e)*z_thick
  !              i_ctr=i_ctr+1
  !            ENDIF
  ! ! write(*,*)'dual coeff',jv,jb,ie,&
  ! ! & p_int_state(1)%edge2vert_coeff(jv,jb,ie,1),&
  ! ! & p_int_state(1)%edge2vert_coeff(jv,jb,ie,2), z_thick
  !
  !         END DO
  !         IF(z_weight/=0.0_wp)THEN
  !           IF(i_ctr<=1)THEN
  !             u_v_cc(jv,jb)%x=0.0_wp
  !           ELSE
  !             u_v_cc(jv,jb)%x = u_v_cc(jv,jb)%x/z_weight
  !           ENDIF
  !         ELSE
  !           u_v_cc(jv,jb)%x=0.0_wp
  !         ENDIF
  !
  !         !Multiply vector with global vorticity
  !         u_v_cc(jv,jb)%x = u_v_cc(jv,jb)%x * vort_v(jv,jk,jb)
  ! ! IF(z_weight/=0.0_wp)THEN
  ! ! write(*,*)'compare',jv,jb,u_v_cc(jv,jb)%x
  ! ! write(*,*)'compare',jv,jb,u_tmp(jv,jb)%x
  ! !ENDIF
  !       END DO VERT_IDX_LOOP
  !     END DO VERT_BLK_LOOP
  ! !$OMP END DO
  !
  !
  ! !$OMP DO PRIVATE(jb,i_startidx_e,i_endidx_e,il_v1,ib_v1,il_v2,ib_v2)
  !     EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e
  !
  !       CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,&
  !                        & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
  !
  !       EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
  !         !Get indices of two adjacent vertices
  !         il_v1 = p_patch%edges%vertex_idx(je,jb,1)
  !         ib_v1 = p_patch%edges%vertex_blk(je,jb,1)
  !         il_v2 = p_patch%edges%vertex_idx(je,jb,2)
  !         ib_v2 = p_patch%edges%vertex_blk(je,jb,2)
  !
  ! !        IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
  !           vn_out_e(je,jk,jb) = &
  !    &- DOT_PRODUCT(u_v_cc(il_v2,ib_v2)%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
  !    &+ DOT_PRODUCT(u_v_cc(il_v1,ib_v1)%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
  ! !        ELSE
  ! !          vn_out_e(je,jk,jb)= 0.0_wp
  ! !        ENDIF
  ! ! write(*,*)'vertex vect 1:',il_v1,ib_v1,(jb-1)*nproma+je,jb,je, u_v_cc(il_v1,ib_v1)%x
  ! ! write(*,*)'vertex vect 2:',il_v2,ib_v2,(jb-1)*nproma+je,jb,je, u_v_cc(il_v2,ib_v2)%x
  !
  ! !write(*,*)'P^T coeff 1:',je,jb,1, p_int_state(1)%edge2vert_coeff_t(je,jb,1)%x
  ! !write(*,*)'P^T coeff 2:',je,jb,2, p_int_state(1)%edge2vert_coeff_t(je,jb,2)%x
  !       END DO EDGE_IDX_LOOP
  !     END DO EDGE_BLK_LOOP
  ! !$OMP END DO
  ! !$OMP END PARALLEL
  !   END DO LEVEL_LOOP
  !   END FUNCTION dual_flip_flop
  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !
  !>
  !! Function implements discrete scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The functions uses coefficients that are calculated in sbr
  !! "init_scalar_product".
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-9)
  !!  no-mpi parallelized
  !!
  FUNCTION dual_flip_flop2(p_patch, vn_old_e, vn_new_e, vort_v,&
    & h_e,p_vn, opt_slev, opt_elev) result(vn_out_e)
    
    TYPE(t_patch), INTENT(in)    :: p_patch
    REAL(wp), INTENT(inout)      :: vn_old_e(:,:,:)
    REAL(wp), INTENT(inout)      :: vn_new_e(:,:,:)
    REAL(wp), INTENT(inout)         :: vort_v(:,:,:)
    TYPE(t_cartesian_coordinates):: p_vn(nproma,n_zlev,p_patch%nblks_c)
    REAL(wp), INTENT(in)      :: h_e(:,:)
    INTEGER, INTENT(in), OPTIONAL ::  opt_slev  ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL ::  opt_elev  ! optional vertical end level
    REAL(wp)                  :: vn_out_e(SIZE(vn_old_e,1), SIZE(vn_old_e,2), SIZE(vn_old_e,3))
    
    !Local variables
    !INTEGER, PARAMETER :: no_vert_edges = 6
    INTEGER :: slev, elev
    !INTEGER :: rl_start_v, rl_end
    INTEGER :: rl_start_e, rl_end_e
    !INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: il_v1, il_v2, ib_v1, ib_v2, il_e, ib_e
    INTEGER :: je, jb, ie, jk!, iie_v1, iie_v2!, jv
    !INTEGER :: ic1, ib1, ic2, ib2
    !INTEGER :: il_star1,ib_star1,il_star2,ib_star2, ie_star1, ie_star2
    
    REAL(wp) :: z_thick, z_weight!, z_tmp1, z_tmp2, z_sum1, z_sum2
    !REAL(wp) :: vn_tmp!z_vn
    REAL(wp) :: z_vn_e(nproma,n_zlev,p_patch%nblks_e)
    TYPE(t_cartesian_coordinates) :: u_v1_cc!(SIZE(vort_v,1),SIZE(vort_v,3))
    TYPE(t_cartesian_coordinates) :: u_v2_cc!(SIZE(vort_v,1),SIZE(vort_v,3))
    INTEGER :: i_v1_ctr, i_v2_ctr
    
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:dual_flip_flop')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    vn_out_e(:,:,:) = 0.0_wp
    
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
    
    rl_start_e = 1
    rl_end_e     = min_rledge
    i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
    i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
    
    z_vn_e = vn_old_e!ab_gam*vn_new_e + (1.0_wp-ab_gam)*vn_old_e
    
    level_loop: DO jk = slev, elev
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx_e,i_endidx_e,il_v1,ib_v1,il_v2,ib_v2)
      edge_blk_loop: DO jb = i_startblk_e, i_endblk_e
        
        CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,&
          & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
        
        edge_idx_loop: DO je =  i_startidx_e, i_endidx_e
          
          IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
            
            !Get indices of two adjacent vertices
            il_v1 = p_patch%edges%vertex_idx(je,jb,1)
            ib_v1 = p_patch%edges%vertex_blk(je,jb,1)
            il_v2 = p_patch%edges%vertex_idx(je,jb,2)
            ib_v2 = p_patch%edges%vertex_blk(je,jb,2)
            
            u_v1_cc%x  = 0.0_wp
            z_weight   = 0.0_wp
            i_v1_ctr = 0
            
            DO ie=1, p_patch%verts%num_edges(il_v1,ib_v1) !no_vert_edges
              
              il_e = p_patch%verts%edge_idx(il_v1,ib_v1,ie)
              ib_e = p_patch%verts%edge_blk(il_v1,ib_v1,ie)
              
              z_thick = v_base%del_zlev_m(jk)
              IF ( iswm_oce == 1 ) THEN
                z_thick    = v_base%del_zlev_m(jk) +h_e(il_e,ib_e)
              ELSEIF( iswm_oce /= 1 ) THEN
                IF (jk == 1 )THEN
                  z_thick = v_base%del_zlev_m(jk) !+ h_e(il_e,ib_e)
                ENDIF
              ENDIF
              
              IF ( v_base%lsm_oce_e(il_e,jk,ib_e) < sea_boundary ) THEN
                i_v1_ctr=i_v1_ctr+1
                !              z_weight = z_weight + p_int_state(1)%variable_dual_vol_norm(il_v1,ib_v1,ie)!*z_thick
                !            ELSE
                !              z_weight = z_weight + 0.5_wp*p_int_state(1)%variable_dual_vol_norm(il_v1,ib_v1,ie)!*z_thick
              ENDIF
              z_weight = z_weight + p_int_state(1)%variable_dual_vol_norm(il_v1,ib_v1,ie)!*z_thick
              
              u_v1_cc%x = u_v1_cc%x +                      &
                & p_int_state(1)%edge2vert_coeff_cc(il_v1,ib_v1,ie)%x * &
                & z_vn_e(il_e,jk,ib_e)!*z_thick
            END DO
            
            !IF(z_weight/=0.0_wp)THEN
            u_v1_cc%x = u_v1_cc%x/(z_weight)!p_patch%verts%dual_area(il_v1,ib_v1)!z_weight
            !ENDIF
            
            u_v2_cc%x   = 0.0_wp
            z_weight    = 0.0_wp
            i_v2_ctr    = 0
            
            DO ie=1, p_patch%verts%num_edges(il_v2,ib_v2) !no_vert_edges
              
              il_e = p_patch%verts%edge_idx(il_v2,ib_v2,ie)
              ib_e = p_patch%verts%edge_blk(il_v2,ib_v2,ie)
              
              z_thick = v_base%del_zlev_m(jk)
              IF ( iswm_oce == 1 ) THEN
                z_thick    = v_base%del_zlev_m(jk) + h_e(il_e,ib_e)
              ELSEIF( iswm_oce /= 1 ) THEN
                IF (jk == 1 )THEN
                  z_thick = v_base%del_zlev_m(jk) !+ h_e(il_e,ib_e)
                ENDIF
              ENDIF
              
              IF ( v_base%lsm_oce_e(il_e,jk,ib_e)  < sea_boundary ) THEN
                i_v2_ctr=i_v2_ctr+1
                !              z_weight = z_weight + p_int_state(1)%variable_dual_vol_norm(il_v2,ib_v2,ie)!*z_thick
                !            ELSE
                !              z_weight = z_weight + 0.5_wp*p_int_state(1)%variable_dual_vol_norm(il_v2,ib_v2,ie)!*z_thick
              ENDIF
              z_weight = z_weight + p_int_state(1)%variable_dual_vol_norm(il_v2,ib_v2,ie)!*z_thick
              
              u_v2_cc%x = u_v2_cc%x +                      &
                & p_int_state(1)%edge2vert_coeff_cc(il_v2,ib_v2,ie)%x * &
                & z_vn_e(il_e,jk,ib_e)!*z_thick
              
            END DO
            !IF(z_weight/=0.0_wp)THEN
            u_v2_cc%x = u_v2_cc%x/(z_weight)!p_patch%verts%dual_area(il_v2,ib_v2)!            /z_weight
            !ENDIF
            
            
            IF(   i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
              & .AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
              
              u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
              u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
              
              vn_out_e(je,jk,jb) = &
                & - DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
                & + DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
              
              
            ELSEIF( i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
              & .AND.i_v2_ctr<p_patch%verts%num_edges(il_v2,ib_v2))THEN
              
              IF(i_v2_ctr>=p_patch%verts%num_edges(il_v2,ib_v2)-3)THEN
                u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
                u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
                
                vn_out_e(je,jk,jb) = &
                  & - DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
                  & + DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
              ELSE
                u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
                
                vn_out_e(je,jk,jb) = &
                  & 2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
              ENDIF
              
            ELSEIF( i_v1_ctr<p_patch%verts%num_edges(il_v1,ib_v1)&
              & .AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
              
              
              
              IF(i_v1_ctr>=p_patch%verts%num_edges(il_v1,ib_v1)-3)THEN
                u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
                u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
                
                vn_out_e(je,jk,jb) = &
                  & - DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
                  & + DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
              ELSE
                u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
                
                vn_out_e(je,jk,jb) = - 2.0_wp * &
                  & DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)
              ENDIF
              
            ELSE
              
              u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
              u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
              
              vn_out_e(je,jk,jb) = &
                & - DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
                & + DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
              
              
              
              !              ic1 = p_patch%edges%cell_idx(je,jb,1)
              !              ib1 = p_patch%edges%cell_blk(je,jb,1)
              !              ic2 = p_patch%edges%cell_idx(je,jb,2)
              !              ib2 = p_patch%edges%cell_blk(je,jb,2)
              !
              !              vn_tmp = 0.5_wp*&
              !              &(DOT_PRODUCT(p_vn(ic1,jk,ib1)%x,p_patch%edges%dual_cart_normal(je,jb)%x)&
              !              &+DOT_PRODUCT(p_vn(ic2,jk,ib2)%x,p_patch%edges%dual_cart_normal(je,jb)%x))
              !
              !
              !               vn_out_e(je,jk,jb) =&
              !               & vn_tmp*p_patch%edges%f_e(je,jb)&
              !               &*0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2))
              
            ENDIF
            ! !--------------------------------------------------
            !               IF(   i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
            !               &.AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                 u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
            !                 u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
            !
            !               ELSEIF(   i_v1_ctr/=p_patch%verts%num_edges(il_v1,ib_v1)&
            !               &.AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                 u_v1_cc%x=u_v1_cc%x*(p_patch%verts%f_v(il_v1,ib_v1))
            !                 u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
            !
            !               ELSEIF(   i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
            !               &.AND.i_v2_ctr/=p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                 u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
            !                 u_v2_cc%x=u_v2_cc%x*(p_patch%verts%f_v(il_v2,ib_v2))
            !               ELSEIF(   i_v1_ctr/=p_patch%verts%num_edges(il_v1,ib_v1)&
            !               &.AND.i_v2_ctr/=p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !                 u_v1_cc%x=u_v1_cc%x*(p_patch%verts%f_v(il_v1,ib_v1))
            !                 u_v2_cc%x=u_v2_cc%x*(p_patch%verts%f_v(il_v2,ib_v2))
            !               ENDIF
            
            !                  u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
            !                  u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
            !
            !                  vn_out_e(je,jk,jb) = &
            !                 &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !                 &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !----------------------------------------------------------------------------------------------
            !               vn_out_e(je,jk,jb) = &
            !              &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !              &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !
            ! !             z_thick = v_base%del_zlev_m(jk)
            ! !             IF ( iswm_oce == 1 ) THEN
            ! !               z_thick = h_e(il_e,ib_e)
            ! !             ELSEIF( iswm_oce /= 1 ) THEN
            ! !               IF (jk == 1 )THEN
            ! !                 z_thick = v_base%del_zlev_m(jk) + h_e(il_e,ib_e)
            ! !               ENDIF
            ! !             ENDIF
            !
            ! !             vn_out_e(je,jk,jb)= &!"z_thick*&
            ! !             &vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            ! !             &+(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            ! !            &/(i_v1_ctr+i_v2_ctr))
            !
            !           vn_out_e(je,jk,jb)= &!"z_thick*&
            !             &vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            !             +0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2)))
            ! !
            ! ! !            vn_out_e(je,jk,jb)= &!"z_thick*&
            ! ! !            &vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)
            !-------------------------------------------------------------
            
            
            !              IF(   i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
            !              &.AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                vn_out_e(je,jk,jb) = &
            !            &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !            &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !
            !                !vn_out_e(je,jk,jb)= &!"z_thick*&
            !                !&vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            !                !&+0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2)))
            !                vn_out_e(je,jk,jb) =&
            !                & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            !                &*(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            !                &/(i_v1_ctr+i_v2_ctr)
            ! !               ic1 = p_patch%edges%cell_idx(je,jb,1)
            ! !               ib1 = p_patch%edges%cell_blk(je,jb,1)
            ! !               ic2 = p_patch%edges%cell_idx(je,jb,2)
            ! !               ib2 = p_patch%edges%cell_blk(je,jb,2)
            ! !              vn_tmp = 0.5_wp*&
            ! !              &(DOT_PRODUCT(p_vn(ic1,jk,ib1)%x,p_patch%edges%dual_cart_normal(je,jb)%x)&
            ! !              &+DOT_PRODUCT(p_vn(ic2,jk,ib2)%x,p_patch%edges%dual_cart_normal(je,jb)%x))
            ! !              vn_tmp = 2.0_wp*&
            ! !              &(DOT_PRODUCT(p_vn(ic1,jk,ib1)%x,p_patch%edges%dual_cart_normal(je,jb)%x)&
            ! !              &+DOT_PRODUCT(p_vn(ic2,jk,ib2)%x,p_patch%edges%dual_cart_normal(je,jb)%x))
            !
            !               vn_out_e(je,jk,jb) =&
            !              & vn_old_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            !              &*0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2))
            ! !
            ! ! !write(*,*)'compare', je,jb,vn_out_e(je,jk,jb), vn_tmp
            ! ! !          ENDIF
            ! ! !   vn_out_e(je,jk,jb)  = vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            ! ! !   &+0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2)))
            ! ! ! write(*,*)'comparison',je,jb,&
            ! ! ! & DOT_PRODUCT(u_v1_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x),&
            ! ! ! &  DOT_PRODUCT(u_v2_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x),&
            ! ! ! &- DOT_PRODUCT(u_v2_p_int_state(1)_base%edge2vert_coeff_cc_t(je,jb,2)%x),&
            ! ! ! &DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x),vn_out_e(je,jk,jb)
            ! ! ! ! IF(jk==1)THEN
            ! ! ! ! write(*,*)'dual-flip-flop: ',je,jb,vn_out_e(je,jk,jb),&
            ! ! ! ! &- 2.0_wp*DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x),&
            ! ! ! ! &  2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! ! ! ! write(*,*)'dual-flip-flop2:',je,jb,&
            ! ! ! ! &2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_vector_cc(il_v1,ib_v1,iie_v1)%x),&
            ! ! ! ! &-2.0_wp*DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_vector_cc(il_v2,ib_v2,iie_v2)%x )
            ! ! ! ! ENDIF
            !
            !
            !              IF( i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
            !              &.AND.  i_v2_ctr< p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                IF(i_v2_ctr<=2)THEN
            !                  !vort_v(il_v2,jk,ib_v2) = 0.0_wp
            !                  !u_v2_cc%x              = 0.0_wp
            !                  vn_out_e(je,jk,jb) = &
            !           & 2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !
            !                  vn_out_e(je,jk,jb) =&
            !                  & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)*vort_v(il_v1,jk,ib_v1)
            !                 ENDIF
            ! !                 ELSEIF(i_v2_ctr>=3)THEN
            ! !                   vn_out_e(je,jk,jb) = &
            ! !                   &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            ! !                   &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! !                   vn_out_e(je,jk,jb) =&
            ! !                   & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            ! !                   &*(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            ! !                   &/(i_v1_ctr+i_v2_ctr)
            ! !                 ENDIF
            !              ENDIF
            ! !            vn_out_e(je,jk,jb) = &
            ! !             &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            ! !             &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! !
            ! !             vn_out_e(je,jk,jb)= &!"z_thick*&
            ! !             &vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            ! !             &+(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            ! !             &/(i_v1_ctr+i_v2_ctr))
            ! ! !             z_mean = (REAL(i_v1_ctr,wp)*vort_v(il_v1,jk,ib_v1)&
            ! ! !                     &+REAL(i_v2_ctr,wp)*vort_v(il_v2,jk,ib_v2))&
            ! ! !                     &/REAL(i_v1_ctr+i_v2_ctr,wp)
            ! ! !             u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
            ! ! !             u_v2_cc%x=u_v2_cc%x*(z_mean+p_patch%verts%f_v(il_v2,ib_v2))
            ! ! !             !vn_out_e(je,jk,jb) = &
            ! ! !             !& 2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! ! !              vn_out_e(je,jk,jb) = DOT_PRODUCT(u_v1_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x)
            ! ! !
            !               IF( i_v1_ctr<p_patch%verts%num_edges(il_v1,ib_v1)&
            !               &.AND.  i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                IF(i_v1_ctr<=2)THEN
            !                  !vort_v(il_v1,jk,ib_v1) = 0.0_wp
            !                  !u_v1_cc%x              = 0.0_wp
            !                  vn_out_e(je,jk,jb) = &
            !                  & -2.0_wp*DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)
            !
            !                  vn_out_e(je,jk,jb) =&
            !                  & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)*vort_v(il_v1,jk,ib_v1)
            !                 ENDIF
            !                ELSEIF(i_v1_ctr>=3)THEN
            !                  vn_out_e(je,jk,jb) = &
            !           &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !           &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !
            !                  vn_out_e(je,jk,jb) =&
            !                  & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            !                  &*(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            !                  &/(i_v1_ctr+i_v2_ctr)
            !                ENDIF
            !              ENDIF
            ! !            vn_out_e(je,jk,jb) = &
            ! !             &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            ! !             &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! !
            ! !             vn_out_e(je,jk,jb)= &!"z_thick*&
            ! !             &vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            ! !             &+(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            ! !             &/(i_v1_ctr+i_v2_ctr))
            ! ! !
            ! ! !              z_mean = (REAL(i_v1_ctr,wp)*vort_v(il_v1,jk,ib_v1)&
            ! ! !                      &+REAL(i_v2_ctr,wp)*vort_v(il_v2,jk,ib_v2))&
            ! ! !                      &/REAL(i_v1_ctr+i_v2_ctr,wp)
            ! ! !
            ! ! !              u_v1_cc%x=u_v1_cc%x*(z_mean+p_patch%verts%f_v(il_v1,ib_v1))
            ! ! !              u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
            ! ! ! !  !           vn_out_e(je,jk,jb) = &
            ! ! ! !  !          & -2.0_wp*DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)
            ! ! ! !
            ! ! !              vn_out_e(je,jk,jb) = DOT_PRODUCT(u_v2_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x)
            ! ! !
            !             ELSEIF( i_v1_ctr<p_patch%verts%num_edges(il_v1,ib_v1)&
            !              &.AND.  i_v2_ctr<p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            ! write(*,*)'incomplete dual cells: indices',il_v1,ib_v1,il_v2,ib_v2, je,jb
            ! ! ! !
            ! !ELSE
            ! ! ! !    z_mean = (REAL(i_v1_ctr,wp)*vort_v(il_v1,jk,ib_v1)+REAL(i_v2_ctr,wp)*vort_v(il_v2,jk,ib_v2))&
            ! ! ! !    &/REAL(i_v1_ctr+i_v2_ctr,wp)
            ! ! ! !    u_v1_cc%x=u_v1_cc%x*(z_mean+p_patch%verts%f_v(il_v1,ib_v1))
            ! ! ! !    u_v2_cc%x=u_v2_cc%x*(z_mean+p_patch%verts%f_v(il_v2,ib_v2))
            ! ! ! !            vn_out_e(je,jk,jb) = &
            ! ! ! !           &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            ! ! ! !           &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! ! ! !
            ! !               ic1 = p_patch%edges%cell_idx(je,jb,1)
            ! !               ib1 = p_patch%edges%cell_blk(je,jb,1)
            ! !               ic2 = p_patch%edges%cell_idx(je,jb,2)
            ! !               ib2 = p_patch%edges%cell_blk(je,jb,2)
            ! !
            ! !               vn_out_e(je,jk,jb) = 0.5_wp*&
            ! !               &(DOT_PRODUCT(p_vn(ic1,jk,ib1)%x,p_patch%edges%dual_cart_normal(je,jb)%x)&
            ! !               &+DOT_PRODUCT(p_vn(ic2,jk,ib2)%x,p_patch%edges%dual_cart_normal(je,jb)%x))
            ! !
            ! !               IF(i_v1_ctr<i_v2_ctr)THEN
            ! !                 vn_out_e(je,jk,jb) =&
            ! !                 & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            ! !                 & *vort_v(il_v2,jk,ib_v2)
            ! !               ELSEIF(i_v1_ctr>=i_v2_ctr)THEN
            ! !                 vn_out_e(je,jk,jb) =&
            ! !                 & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            ! !                 & *vort_v(il_v1,jk,ib_v1)
            ! !               ENDIF
            ! ! ! !             vn_out_e(je,jk,jb) =&
            ! ! ! !             & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            ! ! ! !             &*(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            ! ! ! !             &/(i_v1_ctr+i_v2_ctr)
            ! ! !             vn_out_e(je,jk,jb) =&
            ! ! !             & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            ! ! !             &*0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2))
            !             ENDIF
          ELSE
            vn_out_e(je,jk,jb)= 0.0_wp
          ENDIF
        END DO edge_idx_loop
      END DO edge_blk_loop
      
!$OMP END DO
!$OMP END PARALLEL
    END DO level_loop
  END FUNCTION dual_flip_flop2
  !  !-------------------------------------------------------------------------
  
  
  
  
  !-------------------------------------------------------------------------
  !>
  !! Function implements discrete scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The functions uses coefficients that are calculated in sbr
  !! "init_scalar_product".
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-9)
  !!  no-mpi parallelized
  !!
  FUNCTION dual_flip_flop(p_patch, vn_old_e, vn_new_e, vort_v,&
    & h_e,p_vn, opt_slev, opt_elev) result(vn_out_e)
    
    TYPE(t_patch), INTENT(in)    :: p_patch
    REAL(wp), INTENT(inout)      :: vn_old_e(:,:,:)
    REAL(wp), INTENT(inout)      :: vn_new_e(:,:,:)
    REAL(wp), INTENT(inout)         :: vort_v(:,:,:)
    TYPE(t_cartesian_coordinates):: p_vn(nproma,n_zlev,p_patch%nblks_c)
    REAL(wp), INTENT(in)      :: h_e(:,:)
    INTEGER, INTENT(in), OPTIONAL ::  opt_slev  ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL ::  opt_elev  ! optional vertical end level
    REAL(wp)                  :: vn_out_e(SIZE(vn_old_e,1), SIZE(vn_old_e,2), SIZE(vn_old_e,3))
    
    !Local variables
    !INTEGER, PARAMETER :: no_vert_edges = 6
    INTEGER :: slev, elev
    !INTEGER :: rl_start_v, rl_end
    INTEGER :: rl_start_e, rl_end_e
    !INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: il_v1, il_v2, ib_v1, ib_v2, il_e, ib_e
    INTEGER :: je, jb, ie, jk!, iie_v1, iie_v2!, jv
    !INTEGER :: ic1, ib1, ic2,ib2
    !INTEGER :: il_star1,ib_star1,il_star2,ib_star2, ie_star1, ie_star2
    
    REAL(wp) :: z_thick, z_weight!, z_tmp1, z_tmp2, z_sum1, z_sum2
    !REAL(wp) :: vn_tmp!z_vn
    REAL(wp) :: z_vn_e(nproma,n_zlev,p_patch%nblks_e)
    TYPE(t_cartesian_coordinates) :: u_v1_cc!(SIZE(vort_v,1),SIZE(vort_v,3))
    TYPE(t_cartesian_coordinates) :: u_v2_cc!(SIZE(vort_v,1),SIZE(vort_v,3))
    INTEGER :: i_v1_ctr, i_v2_ctr
    
    !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    !  & routine = ('mo_scalar_product:dual_flip_flop')
    !-----------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    vn_out_e(:,:,:) = 0.0_wp
    
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
    
    rl_start_e = 1
    rl_end_e     = min_rledge
    i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
    i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
    
    z_vn_e = vn_old_e!ab_gam*vn_new_e + (1.0_wp-ab_gam)*vn_old_e
    
    level_loop: DO jk = slev, elev
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx_e,i_endidx_e,il_v1,ib_v1,il_v2,ib_v2)
      edge_blk_loop: DO jb = i_startblk_e, i_endblk_e
        
        CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,&
          & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
        
        edge_idx_loop: DO je =  i_startidx_e, i_endidx_e
          
          IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
            
            !Get indices of two adjacent vertices
            il_v1 = p_patch%edges%vertex_idx(je,jb,1)
            ib_v1 = p_patch%edges%vertex_blk(je,jb,1)
            il_v2 = p_patch%edges%vertex_idx(je,jb,2)
            ib_v2 = p_patch%edges%vertex_blk(je,jb,2)
            
            u_v1_cc%x  = 0.0_wp
            z_weight   = 0.0_wp
            i_v1_ctr = 0
            
            DO ie=1, p_patch%verts%num_edges(il_v1,ib_v1) !no_vert_edges
              
              il_e = p_patch%verts%edge_idx(il_v1,ib_v1,ie)
              ib_e = p_patch%verts%edge_blk(il_v1,ib_v1,ie)
              
              z_thick = v_base%del_zlev_m(jk)
              IF ( iswm_oce == 1 ) THEN
                z_thick    = 1.0_wp!v_base%del_zlev_m(jk) +h_e(il_e,ib_e)
              ELSEIF( iswm_oce /= 1 ) THEN
                !IF (jk == 1 )THEN
                z_thick = v_base%del_zlev_m(jk) !+ h_e(il_e,ib_e)
                !ENDIF
              ENDIF
              
              IF ( v_base%lsm_oce_e(il_e,jk,ib_e) < sea_boundary ) THEN
                i_v1_ctr=i_v1_ctr+1
                !              z_weight = z_weight + p_int_state(1)%variable_dual_vol_norm(il_v1,ib_v1,ie)!*z_thick
                !            ELSE
                !              z_weight = z_weight + 0.5_wp*p_int_state(1)%variable_dual_vol_norm(il_v1,ib_v1,ie)!*z_thick
              ENDIF
              z_weight = z_weight + p_int_state(1)%variable_dual_vol_norm(il_v1,ib_v1,ie)!*z_thick
              
              u_v1_cc%x = u_v1_cc%x +                      &
                & p_int_state(1)%edge2vert_coeff_cc(il_v1,ib_v1,ie)%x * &
                & z_vn_e(il_e,jk,ib_e)!*z_thick
              
              ! write(1234,*)'dual-coeff',il_v1,ib_v1,ie,&
              ! &p_int_state(1)%edge2vert_coeff_cc(il_v1,ib_v1,ie)%x,&
              ! &z_weight
              !ENDIF
            END DO
            
            !           IF(z_weight/=0.0_wp)THEN
            u_v1_cc%x = u_v1_cc%x/(z_weight)!p_patch%verts%dual_area(il_v1,ib_v1)!z_weight
            !           ENDIF
            
            u_v2_cc%x   = 0.0_wp
            z_weight    = 0.0_wp
            i_v2_ctr    = 0
            
            DO ie=1, p_patch%verts%num_edges(il_v2,ib_v2) !no_vert_edges
              
              il_e = p_patch%verts%edge_idx(il_v2,ib_v2,ie)
              ib_e = p_patch%verts%edge_blk(il_v2,ib_v2,ie)
              
              z_thick = v_base%del_zlev_m(jk)
              IF ( iswm_oce == 1 ) THEN
                z_thick    = 1.0_wp!v_base%del_zlev_m(jk) + h_e(il_e,ib_e)
              ELSEIF( iswm_oce /= 1 ) THEN
                !IF (jk == 1 )THEN
                z_thick = v_base%del_zlev_m(jk) !+ h_e(il_e,ib_e)
                !ENDIF
              ENDIF
              
              IF ( v_base%lsm_oce_e(il_e,jk,ib_e)  < sea_boundary ) THEN
                i_v2_ctr=i_v2_ctr+1
                !              z_weight = z_weight + p_int_state(1)%variable_dual_vol_norm(il_v2,ib_v2,ie)!*z_thick
                !            ELSE
                !              z_weight = z_weight + 0.5_wp*p_int_state(1)%variable_dual_vol_norm(il_v2,ib_v2,ie)!*z_thick
              ENDIF
              z_weight = z_weight + p_int_state(1)%variable_dual_vol_norm(il_v2,ib_v2,ie)!*z_thick
              
              
              u_v2_cc%x = u_v2_cc%x +                      &
                & p_int_state(1)%edge2vert_coeff_cc(il_v2,ib_v2,ie)%x * &
                & z_vn_e(il_e,jk,ib_e)!*z_thick
              
              !ENDIF
            END DO
            !           IF(z_weight/=0.0_wp)THEN
            u_v2_cc%x = u_v2_cc%x/(z_weight)!p_patch%verts%dual_area(il_v2,ib_v2)!            /z_weight
            !           ENDIF
            
            
            
            
            IF(   i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
              & .AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
              
              u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
              u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
              
              vn_out_e(je,jk,jb) = &
                & - DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
                & + DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
              
              
            ELSEIF( i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
              & .AND.i_v2_ctr<p_patch%verts%num_edges(il_v2,ib_v2))THEN
              
              IF(i_v2_ctr>=p_patch%verts%num_edges(il_v2,ib_v2)-3)THEN
                u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
                u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
                
                vn_out_e(je,jk,jb) = &
                  & - DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
                  & + DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
              ELSE
                u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
                
                vn_out_e(je,jk,jb) = &
                  & 2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
              ENDIF
              
            ELSEIF( i_v1_ctr<p_patch%verts%num_edges(il_v1,ib_v1)&
              & .AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
              
              
              
              IF(i_v1_ctr>=p_patch%verts%num_edges(il_v1,ib_v1)-3)THEN
                u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
                u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
                
                vn_out_e(je,jk,jb) = &
                  & - DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
                  & + DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
              ELSE
                u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
                
                vn_out_e(je,jk,jb) = - 2.0_wp * &
                  & DOT_PRODUCT(u_v2_cc%x, p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)
              ENDIF
              
            ELSE
              
              u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
              u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
              
              vn_out_e(je,jk,jb) = &
                & - DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
                & + DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
              
              
              
              !              ic1 = p_patch%edges%cell_idx(je,jb,1)
              !              ib1 = p_patch%edges%cell_blk(je,jb,1)
              !              ic2 = p_patch%edges%cell_idx(je,jb,2)
              !              ib2 = p_patch%edges%cell_blk(je,jb,2)
              !
              !              vn_tmp = 0.5_wp*&
              !              &(DOT_PRODUCT(p_vn(ic1,jk,ib1)%x,p_patch%edges%dual_cart_normal(je,jb)%x)&
              !              &+DOT_PRODUCT(p_vn(ic2,jk,ib2)%x,p_patch%edges%dual_cart_normal(je,jb)%x))
              !
              !
              !               vn_out_e(je,jk,jb) =&
              !               & vn_tmp*p_patch%edges%f_e(je,jb)&
              !               &*0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2))
              
            ENDIF
            ! !--------------------------------------------------
            !               IF(   i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
            !               &.AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                 u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
            !                 u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
            !
            !               ELSEIF(   i_v1_ctr/=p_patch%verts%num_edges(il_v1,ib_v1)&
            !               &.AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                 u_v1_cc%x=u_v1_cc%x*(p_patch%verts%f_v(il_v1,ib_v1))
            !                 u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
            !
            !               ELSEIF(   i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
            !               &.AND.i_v2_ctr/=p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                 u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
            !                 u_v2_cc%x=u_v2_cc%x*(p_patch%verts%f_v(il_v2,ib_v2))
            !               ELSEIF(   i_v1_ctr/=p_patch%verts%num_edges(il_v1,ib_v1)&
            !               &.AND.i_v2_ctr/=p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !                 u_v1_cc%x=u_v1_cc%x*(p_patch%verts%f_v(il_v1,ib_v1))
            !                 u_v2_cc%x=u_v2_cc%x*(p_patch%verts%f_v(il_v2,ib_v2))
            !               ENDIF
            
            !                  u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
            !                  u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
            !
            !                  vn_out_e(je,jk,jb) = &
            !                 &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !                 &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !----------------------------------------------------------------------------------------------
            !               vn_out_e(je,jk,jb) = &
            !              &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !              &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !
            ! !             z_thick = v_base%del_zlev_m(jk)
            ! !             IF ( iswm_oce == 1 ) THEN
            ! !               z_thick = h_e(il_e,ib_e)
            ! !             ELSEIF( iswm_oce /= 1 ) THEN
            ! !               IF (jk == 1 )THEN
            ! !                 z_thick = v_base%del_zlev_m(jk) + h_e(il_e,ib_e)
            ! !               ENDIF
            ! !             ENDIF
            !
            ! !             vn_out_e(je,jk,jb)= &!"z_thick*&
            ! !             &vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            ! !             &+(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            ! !            &/(i_v1_ctr+i_v2_ctr))
            !
            !           vn_out_e(je,jk,jb)= &!"z_thick*&
            !             &vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            !             +0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2)))
            ! !
            ! ! !            vn_out_e(je,jk,jb)= &!"z_thick*&
            ! ! !            &vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)
            !-------------------------------------------------------------
            
            
            !              IF(   i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
            !              &.AND.i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                vn_out_e(je,jk,jb) = &
            !            &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !            &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !
            !                !vn_out_e(je,jk,jb)= &!"z_thick*&
            !                !&vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            !                !&+0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2)))
            !                vn_out_e(je,jk,jb) =&
            !                & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            !                &*(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            !                &/(i_v1_ctr+i_v2_ctr)
            ! !               ic1 = p_patch%edges%cell_idx(je,jb,1)
            ! !               ib1 = p_patch%edges%cell_blk(je,jb,1)
            ! !               ic2 = p_patch%edges%cell_idx(je,jb,2)
            ! !               ib2 = p_patch%edges%cell_blk(je,jb,2)
            ! !              vn_tmp = 0.5_wp*&
            ! !              &(DOT_PRODUCT(p_vn(ic1,jk,ib1)%x,p_patch%edges%dual_cart_normal(je,jb)%x)&
            ! !              &+DOT_PRODUCT(p_vn(ic2,jk,ib2)%x,p_patch%edges%dual_cart_normal(je,jb)%x))
            ! !              vn_tmp = 2.0_wp*&
            ! !              &(DOT_PRODUCT(p_vn(ic1,jk,ib1)%x,p_patch%edges%dual_cart_normal(je,jb)%x)&
            ! !              &+DOT_PRODUCT(p_vn(ic2,jk,ib2)%x,p_patch%edges%dual_cart_normal(je,jb)%x))
            !
            !               vn_out_e(je,jk,jb) =&
            !              & vn_old_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            !              &*0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2))
            ! !
            ! ! !write(*,*)'compare', je,jb,vn_out_e(je,jk,jb), vn_tmp
            ! ! !          ENDIF
            ! ! !   vn_out_e(je,jk,jb)  = vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            ! ! !   &+0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2)))
            ! ! ! write(*,*)'comparison',je,jb,&
            ! ! ! & DOT_PRODUCT(u_v1_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x),&
            ! ! ! &  DOT_PRODUCT(u_v2_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x),&
            ! ! ! &- DOT_PRODUCT(u_v2_p_int_state(1)_base%edge2vert_coeff_cc_t(je,jb,2)%x),&
            ! ! ! &DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x),vn_out_e(je,jk,jb)
            ! ! ! ! IF(jk==1)THEN
            ! ! ! ! write(*,*)'dual-flip-flop: ',je,jb,vn_out_e(je,jk,jb),&
            ! ! ! ! &- 2.0_wp*DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x),&
            ! ! ! ! &  2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! ! ! ! write(*,*)'dual-flip-flop2:',je,jb,&
            ! ! ! ! &2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_vector_cc(il_v1,ib_v1,iie_v1)%x),&
            ! ! ! ! &-2.0_wp*DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_vector_cc(il_v2,ib_v2,iie_v2)%x )
            ! ! ! ! ENDIF
            !
            !
            !              IF( i_v1_ctr==p_patch%verts%num_edges(il_v1,ib_v1)&
            !              &.AND.  i_v2_ctr< p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                IF(i_v2_ctr<=2)THEN
            !                  !vort_v(il_v2,jk,ib_v2) = 0.0_wp
            !                  !u_v2_cc%x              = 0.0_wp
            !                  vn_out_e(je,jk,jb) = &
            !           & 2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !
            !                  vn_out_e(je,jk,jb) =&
            !                  & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)*vort_v(il_v1,jk,ib_v1)
            !                 ENDIF
            ! !                 ELSEIF(i_v2_ctr>=3)THEN
            ! !                   vn_out_e(je,jk,jb) = &
            ! !                   &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            ! !                   &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! !                   vn_out_e(je,jk,jb) =&
            ! !                   & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            ! !                   &*(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            ! !                   &/(i_v1_ctr+i_v2_ctr)
            ! !                 ENDIF
            !              ENDIF
            ! !            vn_out_e(je,jk,jb) = &
            ! !             &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            ! !             &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! !
            ! !             vn_out_e(je,jk,jb)= &!"z_thick*&
            ! !             &vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            ! !             &+(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            ! !             &/(i_v1_ctr+i_v2_ctr))
            ! ! !             z_mean = (REAL(i_v1_ctr,wp)*vort_v(il_v1,jk,ib_v1)&
            ! ! !                     &+REAL(i_v2_ctr,wp)*vort_v(il_v2,jk,ib_v2))&
            ! ! !                     &/REAL(i_v1_ctr+i_v2_ctr,wp)
            ! ! !             u_v1_cc%x=u_v1_cc%x*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
            ! ! !             u_v2_cc%x=u_v2_cc%x*(z_mean+p_patch%verts%f_v(il_v2,ib_v2))
            ! ! !             !vn_out_e(je,jk,jb) = &
            ! ! !             !& 2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! ! !              vn_out_e(je,jk,jb) = DOT_PRODUCT(u_v1_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x)
            ! ! !
            !               IF( i_v1_ctr<p_patch%verts%num_edges(il_v1,ib_v1)&
            !               &.AND.  i_v2_ctr==p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            !                IF(i_v1_ctr<=2)THEN
            !                  !vort_v(il_v1,jk,ib_v1) = 0.0_wp
            !                  !u_v1_cc%x              = 0.0_wp
            !                  vn_out_e(je,jk,jb) = &
            !                  & -2.0_wp*DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)
            !
            !                  vn_out_e(je,jk,jb) =&
            !                  & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)*vort_v(il_v1,jk,ib_v1)
            !                 ENDIF
            !                ELSEIF(i_v1_ctr>=3)THEN
            !                  vn_out_e(je,jk,jb) = &
            !           &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            !           &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            !
            !                  vn_out_e(je,jk,jb) =&
            !                  & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            !                  &*(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            !                  &/(i_v1_ctr+i_v2_ctr)
            !                ENDIF
            !              ENDIF
            ! !            vn_out_e(je,jk,jb) = &
            ! !             &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            ! !             &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! !
            ! !             vn_out_e(je,jk,jb)= &!"z_thick*&
            ! !             &vn_out_e(je,jk,jb)*(p_patch%edges%f_e(je,jb)&
            ! !             &+(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            ! !             &/(i_v1_ctr+i_v2_ctr))
            ! ! !
            ! ! !              z_mean = (REAL(i_v1_ctr,wp)*vort_v(il_v1,jk,ib_v1)&
            ! ! !                      &+REAL(i_v2_ctr,wp)*vort_v(il_v2,jk,ib_v2))&
            ! ! !                      &/REAL(i_v1_ctr+i_v2_ctr,wp)
            ! ! !
            ! ! !              u_v1_cc%x=u_v1_cc%x*(z_mean+p_patch%verts%f_v(il_v1,ib_v1))
            ! ! !              u_v2_cc%x=u_v2_cc%x*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
            ! ! ! !  !           vn_out_e(je,jk,jb) = &
            ! ! ! !  !          & -2.0_wp*DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)
            ! ! ! !
            ! ! !              vn_out_e(je,jk,jb) = DOT_PRODUCT(u_v2_cc%x,p_patch%edges%dual_cart_normal(je,jb)%x)
            ! ! !
            !             ELSEIF( i_v1_ctr<p_patch%verts%num_edges(il_v1,ib_v1)&
            !              &.AND.  i_v2_ctr<p_patch%verts%num_edges(il_v2,ib_v2))THEN
            !
            ! write(*,*)'incomplete dual cells: indices',il_v1,ib_v1,il_v2,ib_v2, je,jb
            ! ! ! !
            ! !ELSE
            ! ! ! !    z_mean = (REAL(i_v1_ctr,wp)*vort_v(il_v1,jk,ib_v1)+REAL(i_v2_ctr,wp)*vort_v(il_v2,jk,ib_v2))&
            ! ! ! !    &/REAL(i_v1_ctr+i_v2_ctr,wp)
            ! ! ! !    u_v1_cc%x=u_v1_cc%x*(z_mean+p_patch%verts%f_v(il_v1,ib_v1))
            ! ! ! !    u_v2_cc%x=u_v2_cc%x*(z_mean+p_patch%verts%f_v(il_v2,ib_v2))
            ! ! ! !            vn_out_e(je,jk,jb) = &
            ! ! ! !           &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
            ! ! ! !           &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
            ! ! ! !
            ! !               ic1 = p_patch%edges%cell_idx(je,jb,1)
            ! !               ib1 = p_patch%edges%cell_blk(je,jb,1)
            ! !               ic2 = p_patch%edges%cell_idx(je,jb,2)
            ! !               ib2 = p_patch%edges%cell_blk(je,jb,2)
            ! !
            ! !               vn_out_e(je,jk,jb) = 0.5_wp*&
            ! !               &(DOT_PRODUCT(p_vn(ic1,jk,ib1)%x,p_patch%edges%dual_cart_normal(je,jb)%x)&
            ! !               &+DOT_PRODUCT(p_vn(ic2,jk,ib2)%x,p_patch%edges%dual_cart_normal(je,jb)%x))
            ! !
            ! !               IF(i_v1_ctr<i_v2_ctr)THEN
            ! !                 vn_out_e(je,jk,jb) =&
            ! !                 & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            ! !                 & *vort_v(il_v2,jk,ib_v2)
            ! !               ELSEIF(i_v1_ctr>=i_v2_ctr)THEN
            ! !                 vn_out_e(je,jk,jb) =&
            ! !                 & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            ! !                 & *vort_v(il_v1,jk,ib_v1)
            ! !               ENDIF
            ! ! ! !             vn_out_e(je,jk,jb) =&
            ! ! ! !             & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            ! ! ! !             &*(i_v1_ctr*vort_v(il_v1,jk,ib_v1)+i_v2_ctr*vort_v(il_v2,jk,ib_v2))&
            ! ! ! !             &/(i_v1_ctr+i_v2_ctr)
            ! ! !             vn_out_e(je,jk,jb) =&
            ! ! !             & vn_out_e(je,jk,jb)*p_patch%edges%f_e(je,jb)&
            ! ! !             &*0.5_wp*(vort_v(il_v1,jk,ib_v1)+vort_v(il_v2,jk,ib_v2))
            !             ENDIF
          ELSE
            vn_out_e(je,jk,jb)= 0.0_wp
          ENDIF
        END DO edge_idx_loop
      END DO edge_blk_loop
      
!$OMP END DO
!$OMP END PARALLEL
      
    END DO level_loop
  END FUNCTION dual_flip_flop
  
  
  !-------------------------------------------------------------------------
  !
  !>
  !!
  ! ! SUBROUTINE nonlinear_Coriolis(p_patch, vn, p_vn, p_vn_dual,h_e,              &
  ! !                         & vort_v, vort_flux)
  ! !
  ! !   TYPE(t_patch), INTENT(IN)      :: p_patch
  ! !   REAL(wp), INTENT(in)           :: vn(:,:,:)
  ! !   TYPE(t_cartesian_coordinates)  :: p_vn(nproma,n_zlev,p_patch%nblks_c)
  ! !   TYPE(t_cartesian_coordinates)  :: p_vn_dual(nproma,n_zlev,p_patch%nblks_v)
  ! !   REAL(wp), INTENT(IN)           :: h_e(:,:)
  ! !   REAL(wp), INTENT(inout)        :: vort_v(:,:,:)
  ! !   REAL(wp), INTENT(inout)        :: vort_flux(:,:,:)
  ! !
  ! ! !Local variables
  ! ! !
  ! ! REAL(wp) :: z_vort_tmp
  ! ! REAL(wp) :: z_weight(nproma,n_zlev,p_patch%nblks_v)
  ! ! REAL(wp) :: zarea_fraction
  ! ! REAL(wp) :: z_area_scaled
  ! !
  ! ! INTEGER :: slev, elev     ! vertical start and end level
  ! ! INTEGER :: jv, jk, jb, jev,je
  ! ! INTEGER :: ile, ibe, il, ib, ill
  ! ! INTEGER :: ik, ikk
  ! ! INTEGER :: rl_start_e, rl_end_e
  ! ! INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v, i_nchdom
  ! ! INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
  ! !
  ! ! INTEGER :: i_bdr_ctr
  ! ! INTEGER :: iboundary_edge_idx_1, iboundary_edge_blk_1
  ! ! INTEGER :: iboundary_edge_idx_2, iboundary_edge_blk_2
  ! ! INTEGER :: icell_idx_1, icell_blk_1
  ! ! INTEGER :: icell_idx_2, icell_blk_2
  ! ! INTEGER :: il_v1, il_v2,ib_v1, ib_v2
  ! ! INTEGER  :: i_v_ctr(nproma,n_zlev,p_patch%nblks_v)
  ! ! TYPE(t_cartesian_coordinates) :: cell1_cc, cell2_cc, vertex_cc
  ! ! TYPE(t_cartesian_coordinates) :: u_v_cc(nproma,n_zlev,p_patch%nblks_v)
  ! ! TYPE(t_cartesian_coordinates) :: u_v1_cc, u_v2_cc
  ! ! INTEGER,PARAMETER :: ino_dual_edges = 6
  ! ! INTEGER,PARAMETER :: rl_start_v = 2
  ! ! INTEGER,PARAMETER :: rl_end_v   = min_rlvert
  ! ! !-----------------------------------------------------------------------
  ! ! slev         = 1
  ! ! elev         = n_zlev
  ! ! rl_start_e   = 1
  ! ! rl_end_e     = min_rledge
  ! !
  ! ! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
  ! ! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
  ! !
  ! ! i_startblk_v = p_patch%verts%start_blk(rl_start_v,1)
  ! ! i_endblk_v   = p_patch%verts%end_blk(rl_end_v,1)
  ! !
  ! ! i_v_ctr = 0
  ! ! ! #slo# due to nag -nan compiler-option
  ! ! vort_v(:,:,:) = 0.0_wp
  ! !
  ! ! !In this loop vorticity and velocity reconstruction at vertices are calculated
  ! ! DO jb = i_startblk_v, i_endblk_v
  ! !
  ! !   CALL get_indices_v(p_patch, jb, i_startblk_v, i_endblk_v, &
  ! !                      i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
  ! !   DO jk = slev, elev
  ! ! !$OMP PARALLEL DO SCHEDULE(runtime) DEFAULT(PRIVATE)  &
  ! ! !$OMP   SHARED(u_vec_e,v_vec_e,ptr_patch,rot_vec_v,jb) FIRSTPRIVATE(jk)
  ! !     DO jv = i_startidx_v, i_endidx_v
  ! !
  ! !       z_vort_tmp          = 0.0_wp
  ! !       zarea_fraction      = 0.0_wp
  ! !       !z_weight(jv,jk,jb) = 0.0_wp
  ! !       !u_v_cc(jv,jk,jb)%x = 0.0_wp
  ! !
  ! !       vertex_cc = gc2cc(p_patch%verts%vertex(jv,jb))
  ! !
  ! !       DO jev = 1, p_patch%verts%num_edges(jv,jb)
  ! !
  ! !         ! get line and block indices of edge jev around vertex jv
  ! !         ile = p_patch%verts%edge_idx(jv,jb,jev)
  ! !         ibe = p_patch%verts%edge_blk(jv,jb,jev)
  ! !
  ! !         !Check, if edge is sea or boundary edge and take care of dummy edge
  ! !         ! edge with indices ile, ibe is sea edge
  ! !         IF ( v_base%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN
  ! !
  ! !           !Distinguish the following cases
  ! !           ! edge ie_k is
  ! !           !a) ocean edge: compute as usual,
  ! !           !b) land edge: do not consider it
  ! !           !c) boundary edge take:
  ! !           !  no-slip boundary condition:  normal and tangential velocity at boundary are zero
  ! !           ! sea, sea_boundary, boundary (edges only), land_boundary, land =
  ! !           !  -2,      -1,         0,                  1,             2
  ! !           !add contribution of normal velocity at edge (ile,ibe) to rotation
  ! !           z_vort_tmp = z_vort_tmp + vn(ile,jk,ibe)                   &
  ! !             & * p_patch%edges%dual_edge_length(ile,ibe)  &
  ! !             & * p_patch%verts%edge_orientation(jv,jb,jev)
  ! !
  ! !           u_v_cc(jv,jk,jb)%x = u_v_cc(jv,jk,jb)%x        &
  ! !              &           +p_int_state(1)%edge2vert_coeff_cc(jv,jb,jev)%x &
  ! !              &           *vn(ile,jk,ibe)!*z_thick
  ! !
  ! !           !z_weight might be an alternative to dual_area and can include
  ! !           !varying height in top layer. Differences have to be explored.
  ! !           !z_weight(jv,jk,jb) = z_weight(jv,jk,jb) &
  ! !           !&+ p_int_state(1)%variable_dual_vol_norm(jv,jb,jev)!*z_thick
  ! !
  ! !
  ! !           !increase wet cell ctr
  ! !           i_v_ctr(jv,jk,jb)=i_v_ctr(jv,jk,jb)+1
  ! !
  ! ! !         ! edge with indices ile, ibe is boundary edge
  ! !          ELSE IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
  ! !
  ! !            !increase boundary edge counter
  ! !            i_bdr_ctr=i_bdr_ctr+1
  ! !
  ! !            !Store actual boundary edge indices
  ! !            IF(i_bdr_ctr==1)THEN
  ! !              iboundary_edge_idx_1 = ile
  ! !              iboundary_edge_blk_1 = ibe
  ! !            ELSEIF(i_bdr_ctr==2)THEN
  ! !              iboundary_edge_idx_2 = ile
  ! !              iboundary_edge_blk_2 = ibe
  ! !            ELSE
  ! !            !only 2 boundray edges per dual loop are allowed: somethings wrong withe the grid
  ! !              CALL message (TRIM('sbr nonlinear Coriolis'), &
  ! !              &'more than 2 boundary edges per dual loop: something is wrong with the grid')
  ! !              CALL finish ('TRIM(sbr nonlinear Coriolis)','Grid-boundary error !!')
  ! !            ENDIF
  ! !
  ! ! !            z_vort_tmp = z_vort_tmp + vn(ile,jk,ibe)                          &
  ! ! !              &  * 0.5_wp * p_patch%edges%dual_edge_length(ile,ibe)  &
  ! ! !              &           * p_patch%verts%edge_orientation(jv,jb,jev)
  ! !          END IF
  ! !       END DO
  ! !       !
  ! !       !divide by hex/pentagon area, if all dual cells are in the ocean interior
  ! !       !divide by apropriate fraction if boundaries are involved
  ! !       IF ( i_v_ctr(jv,jk,jb) == p_patch%verts%num_edges(jv,jb) ) THEN
  ! !
  ! !         !vort_v(jv,jk,jb) = z_vort_tmp /p_patch%verts%dual_area(jv,jb)! (re*re*z_weight(jv,jk,jb))!
  ! !
  ! !         z_area_scaled       = p_patch%verts%dual_area(jv,jb)/(re*re)
  ! !         u_v_cc(jv,jk,jb)%x  = u_v_cc(jv,jk,jb)%x/z_area_scaled
  ! !
  ! !         !u_v_cc(jv,jk,jb)%x  = u_v_cc(jv,jk,jb)%x/z_weight(jv,jk,jb)
  ! !
  ! !       ELSE!boundary edges are involved
  ! !         DO jev = 1, p_patch%verts%num_edges(jv,jb)
  ! !           !
  ! !           ! get line and block indices of edge jev around vertex jv
  ! !           ile = p_patch%verts%edge_idx(jv,jb,jev)
  ! !           ibe = p_patch%verts%edge_blk(jv,jb,jev)
  ! !
  ! !           !get neighbor cells
  ! !           icell_idx_1 = p_patch%edges%cell_idx(ile,ibe,1)
  ! !           icell_idx_2 = p_patch%edges%cell_idx(ile,ibe,2)
  ! !           icell_blk_1 = p_patch%edges%cell_blk(ile,ibe,1)
  ! !           icell_blk_2 = p_patch%edges%cell_blk(ile,ibe,2)
  ! !
  ! !           cell1_cc = gc2cc(p_patch%cells%center(icell_idx_1,icell_blk_1))
  ! !           cell2_cc = gc2cc(p_patch%cells%center(icell_idx_2,icell_blk_2))
  ! !
  ! !           !Check, if edge is sea or boundary edge and take care of dummy edge
  ! !           ! edge with indices ile, ibe is sea edge
  ! !           !Add up for wet dual area.
  ! !           IF ( v_base%lsm_oce_e(ile,jk,ibe) <= sea_boundary ) THEN
  ! !
  ! !            zarea_fraction = zarea_fraction  &
  ! !               &     + re*re*triangle_area(cell1_cc, vertex_cc, cell2_cc)
  ! !
  ! !           ! edge with indices ile, ibe is boundary edge
  ! !           ELSE IF ( v_base%lsm_oce_e(ile,jk,ibe) == boundary ) THEN
  ! !
  ! !            zarea_fraction = zarea_fraction  &
  ! !               &  + 0.5_wp*re*re*triangle_area(cell1_cc, vertex_cc, cell2_cc)
  ! !
  ! !           END IF
  ! !
  ! !         END DO
  ! !
  ! !         ! no division by zero in case of zero-test (#slo# 2010-06-09)
  ! !         IF (zarea_fraction == 0.0_wp) THEN
  ! !           !vort_v(jv,jk,jb) = 0.0_wp
  ! !           u_v_cc(jv,jk,jb)%x=0.0_wp
  ! !         ELSE
  ! !           !vort_v(jv,jk,jb)    = z_vort_tmp /zarea_fraction!(re*re*z_weight(jv,jk,jb))!
  ! !           z_area_scaled       = zarea_fraction/(re*re)
  ! !           u_v_cc(jv,jk,jb)%x  = u_v_cc(jv,jk,jb)%x/z_area_scaled!z_weight(jv,jk,jb)!
  ! !         ENDIF
  ! !       ENDIF
  ! !     END DO
  ! ! !$OMP END PARALLEL DO
  ! !   END DO
  ! ! END DO
  ! !
  ! !   LEVEL_LOOP: DO jk = slev, elev
  ! ! !$OMP PARALLEL
  ! ! !$OMP DO PRIVATE(jb,i_startidx_e,i_endidx_e,il_v1,ib_v1,il_v2,ib_v2)
  ! !     EDGE_BLK_LOOP: DO jb = i_startblk_e, i_endblk_e
  ! !
  ! !       CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,&
  ! !                        & i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
  ! !
  ! !       EDGE_IDX_LOOP: DO je =  i_startidx_e, i_endidx_e
  ! !
  ! !         IF(v_base%lsm_oce_e(je,jk,jb) <= sea_boundary)THEN
  ! !
  ! !           !Get indices of two adjacent vertices
  ! !           il_v1 = p_patch%edges%vertex_idx(je,jb,1)
  ! !           ib_v1 = p_patch%edges%vertex_blk(je,jb,1)
  ! !           il_v2 = p_patch%edges%vertex_idx(je,jb,2)
  ! !           ib_v2 = p_patch%edges%vertex_blk(je,jb,2)
  ! !
  ! !           u_v1_cc%x=u_v_cc(il_v1,jk,ib_v1)%x&
  ! !             &*(vort_v(il_v1,jk,ib_v1)+p_patch%verts%f_v(il_v1,ib_v1))
  ! !           u_v2_cc%x=u_v_cc(il_v2,jk,ib_v2)%x&
  ! !             &*(vort_v(il_v2,jk,ib_v2)+p_patch%verts%f_v(il_v2,ib_v2))
  ! !
  ! ! !          IF(   i_v_ctr(il_v1,jk,ib_v1)==p_patch%verts%num_edges(il_v1,ib_v1)&
  ! ! !          &.AND.i_v_ctr(il_v2,jk,ib_v2)==p_patch%verts%num_edges(il_v2,ib_v2))THEN
  ! !
  ! !             vort_flux(je,jk,jb) = &
  ! !             &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
  ! !             &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
  ! !
  ! ! !           ELSEIF( i_v_ctr(il_v1,jk,ib_v1)==p_patch%verts%num_edges(il_v1,ib_v1)&
  ! ! !           &.AND.  i_v_ctr(il_v2,jk,ib_v2)<p_patch%verts%num_edges(il_v2,ib_v2))THEN
  ! ! !
  ! ! !             IF(i_v_ctr(il_v2,jk,ib_v2)>=p_patch%verts%num_edges(il_v2,ib_v2)-3)THEN
  ! ! !
  ! ! !               vort_flux(je,jk,jb) = &
  ! ! !               &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
  ! ! !               &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
  ! ! !
  ! ! !             ELSE
  ! ! !               vort_flux(je,jk,jb) = &
  ! ! !               & 2.0_wp*DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
  ! ! !             ENDIF
  ! ! !
  ! ! !           ELSEIF( i_v_ctr(il_v1,jk,ib_v1)<p_patch%verts%num_edges(il_v1,ib_v1)&
  ! ! !             &.AND.i_v_ctr(il_v2,jk,ib_v2)==p_patch%verts%num_edges(il_v2,ib_v2))THEN
  ! ! !
  ! ! !             IF(i_v_ctr(il_v1,jk,ib_v1)>=p_patch%verts%num_edges(il_v1,ib_v1)-3)THEN
  ! ! !
  ! ! !               vort_flux(je,jk,jb) = &
  ! ! !               &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
  ! ! !               &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
  ! ! !
  ! ! !             ELSE
  ! ! !               vort_flux(je,jk,jb) = &
  ! ! !               &- 2.0_wp*DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)
  ! ! !             ENDIF
  ! ! !           ELSE
  ! ! !             vort_flux(je,jk,jb) = &
  ! ! !             &- DOT_PRODUCT(u_v2_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,2)%x)&
  ! ! !             &+ DOT_PRODUCT(u_v1_cc%x,p_int_state(1)%edge2vert_coeff_cc_t(je,jb,1)%x)
  ! ! !         ENDIF
  ! !
  ! !        ELSE
  ! !          vort_flux(je,jk,jb)= 0.0_wp
  ! !        ENDIF
  ! !       END DO EDGE_IDX_LOOP
  ! !     END DO EDGE_BLK_LOOP
  ! ! !$OMP END DO
  ! ! !$OMP END PARALLEL
  ! !   END DO LEVEL_LOOP
  ! !
  ! ! END SUBROUTINE nonlinear_Coriolis
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
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
    IF ( z_s12 == 0.0_wp .OR. z_s23 == 0.0_wp  .OR. z_s31 == 0.0_wp ) THEN
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
  
  
  
END MODULE mo_scalar_product

