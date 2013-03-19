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
  USE mo_impl_constants,     ONLY: sea_boundary, sea
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_oce_state,          ONLY: t_hydro_ocean_diag!, t_hydro_ocean_state
  USE mo_ocean_nml,          ONLY: n_zlev, iswm_oce 
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates,cvec2gvec!, gc2cc, vector_product
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff, no_primal_edges, no_dual_edges
  USE mo_oce_math_operators,  ONLY: rot_vertex_ocean_3d, map_edges2vert_3d
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_e, sync_v,sync_patch_array!,sync_c,  & sync_idx, global_max
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  PUBLIC :: calc_scalar_product_veloc_3d
  PUBLIC :: nonlinear_coriolis_3d
  PUBLIC :: nonlinear_coriolis_3d_old
  PUBLIC :: map_edges2vert_3d
  PUBLIC :: map_edges2cell_3d
  PUBLIC :: map_cell2edges_3D
  PUBLIC :: map_edges2edges_viacell_3d
  PUBLIC :: map_edges2edges_viavert_3D
  PUBLIC :: map_edges2edges_viacell_3D_const_z
  !PRIVATE :: map_cell2edges_2d


  INTERFACE map_edges2edges_viacell_3d 
    MODULE PROCEDURE map_edges2edges_viacell_3d_1lev
    MODULE PROCEDURE map_edges2edges_viacell_3d_mlev
  END INTERFACE

  INTERFACE map_edges2edges_viacell_3d_const_z 
    MODULE PROCEDURE map_edges2edges_viacell_3D_1lev_const_z
    MODULE PROCEDURE map_edges2edges_viacell_3D_mlev_const_z
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
  SUBROUTINE calc_scalar_product_veloc_3D( p_patch_3D, vn_e_old, vn_e_new,&
    & p_diag, p_op_coeff)
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(in)      :: vn_e_old(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(in)      :: vn_e_new(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_hydro_ocean_diag)  :: p_diag
    TYPE(t_operator_coeff)    :: p_op_coeff
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: jc, jb, jk

    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    !-----------------------------------------------------------------------
    slev = 1
    elev = n_zlev

    CALL map_edges2vert_3d(p_patch_3D%p_patch_2D(1), vn_e_old, p_op_coeff%edge2vert_coeff_cc, &
      & p_diag%p_vn_dual)

    !Step 1: Calculation of Pv in cartesian coordinates and of kinetic energy
    CALL map_edges2cell_3d(p_patch_3D, vn_e_old, p_op_coeff, p_diag%p_vn)

    CALL map_cell2edges_3D( p_patch_3D, p_diag%p_vn, p_diag%ptp_vn, p_op_coeff)

!      CALL map_edges2edges_viacell_3D( p_patch_3D,    &
!                                     & vn_e_old,      &
!                                     & p_op_coeff,    &
!                                     & p_diag%ptp_vn)

   CALL sync_patch_array(SYNC_E, p_patch, p_diag%ptp_vn)

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
          IF(p_patch_3D%lsm_c(jc,jk,jb) > sea_boundary)THEN
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

  END SUBROUTINE calc_scalar_product_veloc_3d
  !-------------------------------------------------------------------------
  
 
  !-------------------------------------------------------------------------
  !>
  !!
  !! mpi parallelized by LL, openmp corrected
  SUBROUTINE nonlinear_coriolis_3d(p_patch_3D, vn, p_vn_dual,h_e, vort_v, &
    & p_op_coeff, vort_flux)
    TYPE(t_patch_3D ),TARGET,INTENT(IN) :: p_patch_3D
    REAL(wp), INTENT(INOUT)                    :: vn(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(INout)  :: p_vn_dual(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp), INTENT(IN)                       :: h_e      (nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(INOUT)                    :: vort_v   (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(IN)          :: p_op_coeff
    REAL(wp), INTENT(INOUT)                    :: vort_flux(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

    !Local variables
    !TYPE(t_patch), POINTER         :: p_patch 
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: il_e, ib_e
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    REAL(wp) :: vort_global
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_edges => p_patch%edges%all
    !-----------------------------------------------------------------------
    slev    = 1
    elev    = n_zlev

    CALL rot_vertex_ocean_3d( p_patch_3D, vn, p_vn_dual, p_op_coeff, vort_v)
    CALL sync_patch_array(SYNC_V, p_patch, vort_v)


! !$OMP PARALLEL
! !$OMP DO PRIVATE(jb,jk,je,i_startidx_e,i_endidx_e)
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)

      level_loop: DO jk = slev, elev

        edge_idx_loop: DO je =  i_startidx_e, i_endidx_e

          IF (p_patch_3D%lsm_e(je,jk,jb) == sea) THEN

            vort_flux(je,jk,jb) = 0.0_wp 

            DO neighbor=1,2   
              IF(neighbor==1) ictr = 0
              IF(neighbor==2) ictr = no_dual_edges

              il_v = p_patch%edges%vertex_idx(je,jb,neighbor)
              ib_v = p_patch%edges%vertex_blk(je,jb,neighbor) 

              vort_global = (vort_v(il_v,jk,ib_v) + p_patch%verts%f_v(il_v,ib_v))

              DO vertex_edge=1, p_patch%verts%num_edges(il_v,ib_v)!no_dual_cell_edges

                ictr =ictr+1 

                il_e = p_patch%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = p_patch%verts%edge_blk(il_v,ib_v,vertex_edge)

                vort_flux(je,jk,jb) =  vort_flux(je,jk,jb)+vn(il_e,jk,ib_e)*vort_global&
                                    &*p_op_coeff%edge2edge_viavert_coeff(je,jk,jb,ictr)
              END DO
            END DO
            ELSE
              vort_flux(je,jk,jb)= 0.0_wp
          ENDIF ! (v_base%lsm_e(je,jk,jb) <= sea_boundary)
        END DO edge_idx_loop
      END DO level_loop
    END DO ! jb = all_edges%start_block, all_edges%end_block
! !$OMP END DO NOWAIT
! !$OMP END PARALLEL
  ! LL: no sync required

  END SUBROUTINE nonlinear_coriolis_3d
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !>
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized LL
  SUBROUTINE map_edges2cell_with_height_3d( p_patch_3D, vn_e, p_op_coeff, p_vn_c, h_e,&
    & opt_slev, opt_elev, subset_range)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    ! 3D case: h_e is surface elevation at edges
    TYPE(t_cartesian_coordinates),INTENT(inout):: p_vn_c(:,:,:)  ! outputput (nproma,n_zlev,nblks_c)
    REAL(wp), INTENT(in)                       :: h_e(:,:)       ! SW-case: h_e is thickness at edges
    TYPE(t_operator_coeff)                     :: p_op_coeff
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
   TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => p_patch%cells%all
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

              il_e = p_patch%cells%edge_idx(jc,jb,ie)
              ib_e = p_patch%cells%edge_blk(jc,jb,ie)

              z_thick_e =p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_e(il_e,jk,ib_e)&
              & + h_e(il_e,ib_e)  
              z_weight = z_weight + p_op_coeff%variable_vol_norm(jc,jk,jb,ie) * z_thick_e

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
          DO ie=1, no_primal_edges

            il_e = p_patch%cells%edge_idx(jc,jb,ie)
            ib_e = p_patch%cells%edge_blk(jc,jb,ie)

            z_thick_e = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_e(il_e,slev,ib_e)&
            & + h_e(il_e,ib_e) 
            z_weight = z_weight + p_op_coeff%variable_vol_norm(jc,slev,jb,ie) * z_thick_e

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
            !z_weight = 0.0_wp
            DO ie=1, no_primal_edges

              il_e = p_patch%cells%edge_idx(jc,jb,ie)
              ib_e = p_patch%cells%edge_blk(jc,jb,ie)
              p_vn_c(jc,jk,jb)%x = p_vn_c(jc,jk,jb)%x&
                & + p_op_coeff%edge2cell_coeff_cc(jc,jk,jb,ie)%x&
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
  !!  mpi parallelized LL
  SUBROUTINE map_edges2cell_no_height_3d( p_patch_3D, vn_e, p_op_coeff, p_vn_c, opt_slev, opt_elev, &
    &                                     subset_range)
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(:,:,:)    ! input (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: p_op_coeff
    TYPE(t_cartesian_coordinates)              :: p_vn_c(:,:,:)  ! output (nproma,n_zlev,nblks_c)
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
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_cells => subset_range
    ELSE
      all_cells => p_patch%cells%all
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
  SUBROUTINE map_edges2edges_viacell_3d_mlev( p_patch_3D, vn_e, p_op_coeff, p_vn_e,scalar, opt_slev, opt_elev, &
    &                                     subset_range)
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)       :: p_patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: p_op_coeff
    REAL(wp), INTENT(INOUT)                    :: p_vn_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(IN), OPTIONAL             :: scalar(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)  
    INTEGER, INTENT(in), OPTIONAL              :: opt_slev       ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL              :: opt_elev       ! optional vertical end level
    TYPE(t_subset_range), TARGET,  OPTIONAL    :: subset_range
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr!, neighbor
    INTEGER :: je, jb, jk, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
     TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER         :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_edges => subset_range
    ELSE
      all_edges => p_patch%edges%all
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
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
 
        level_loop_e: DO jk = slev, elev
          edge_idx_loop: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
            ictr            = 0
            p_vn_e(je,jk,jb)= 0.0_wp
            !IF(p_patch_3D%lsm_e(je,jk,jb) == sea)THEN
            il_c = p_patch%edges%cell_idx(je,jb,1)
            ib_c = p_patch%edges%cell_blk(je,jb,1) 
            thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)
              thick_frac = thick_edge/thick_cell
              p_vn_e(je,jk,jb) = p_vn_e(je,jk,jb) &
              &+vn_e(il_e,jk,ib_e)*(p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr))*thick_frac

            END DO

            ictr = no_primal_edges
            il_c = p_patch%edges%cell_idx(je,jb,2)
            ib_c = p_patch%edges%cell_blk(je,jb,2)
            thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)
              thick_frac = thick_edge/thick_cell
              p_vn_e(je,jk,jb) = p_vn_e(je,jk,jb) &
              &+vn_e(il_e,jk,ib_e)*(p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr))*thick_frac

            END DO
            ENDIF
          END DO edge_idx_loop
        END DO level_loop_e
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

    ELSEIF(PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)

        level_loop_e2: DO jk = slev, elev
          edge_idx_loop2: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary) THEN
            ictr            = 0
            p_vn_e(je,jk,jb)= 0.0_wp
            il_c        = p_patch%edges%cell_idx(je,jb,1)
            ib_c        = p_patch%edges%cell_blk(je,jb,1)
            scalar_cell = scalar(il_c,jk,ib_c)
            thick_cell  = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

              DO ie=1, no_primal_edges
                ictr =ictr+1 
                il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
                ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
                thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)
                thick_frac = thick_edge/thick_cell

                p_vn_e(je,jk,jb) = p_vn_e(je,jk,jb) &
                &+vn_e(il_e,jk,ib_e)*scalar_cell   &
                &*(p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr))*thick_frac
              END DO


            ictr        = no_primal_edges
            il_c        = p_patch%edges%cell_idx(je,jb,2)
            ib_c        = p_patch%edges%cell_blk(je,jb,2)
            scalar_cell = scalar(il_c,jk,ib_c)
            thick_cell  = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

              DO ie=1, no_primal_edges
                ictr =ictr+1 
                il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
                ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
                thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)
                thick_frac = thick_edge/thick_cell

                p_vn_e(je,jk,jb) = p_vn_e(je,jk,jb) &
                &+vn_e(il_e,jk,ib_e)*scalar_cell   &
                &*(p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr))*thick_frac

              END DO

            ENDIF
          END DO edge_idx_loop2
        END DO level_loop_e2
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block
    ENDIF

  END SUBROUTINE map_edges2edges_viacell_3D_mlev
  !-----------------------------------------------------------------------------

SUBROUTINE map_edges2edges_viacell_3d_1lev( p_patch_3D, vn_e, p_op_coeff, p_vn_e,scalar,scalar_e, level, &
    &                                     subset_range)
   
    TYPE(t_patch_3D ),TARGET, INTENT(IN):: p_patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: p_op_coeff
    REAL(wp), INTENT(INOUT)                    :: p_vn_e(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(IN), OPTIONAL             :: scalar(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)  
    REAL(wp), INTENT(IN), OPTIONAL             :: scalar_e(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)  
    INTEGER, INTENT(in), OPTIONAL              :: level       ! optional vertical start level
    TYPE(t_subset_range), TARGET,  OPTIONAL    :: subset_range
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, jb, ie, lev
    REAL(wp) :: scalar_cell
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    IF ( PRESENT(subset_range) ) THEN
      all_edges => subset_range
    ELSE
      all_edges => p_patch%edges%all
    ENDIF
    !--------------------------
    IF ( PRESENT(level) ) THEN
      lev = level
    ELSE
      lev = 1
    END IF

    IF(.NOT.PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)

          edge_idx_loop: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,lev,jb) <= sea_boundary) THEN
            ictr          = 0
            p_vn_e(je,jb) = 0.0_wp
            il_c = p_patch%edges%cell_idx(je,jb,1)
            ib_c = p_patch%edges%cell_blk(je,jb,1)

            !thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)
            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
              !thick_edge=p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
              !thick_frac=thick_edge/thick_cell
              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*(p_op_coeff%edge2edge_viacell_coeff(je,lev,jb,ictr))!*thick_frac

            END DO

            ictr = no_primal_edges
            il_c = p_patch%edges%cell_idx(je,jb,2)
            ib_c = p_patch%edges%cell_blk(je,jb,2)

            !thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)
            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
              !thick_edge=p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
              !thick_frac=thick_edge/thick_cell
              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*(p_op_coeff%edge2edge_viacell_coeff(je,lev,jb,ictr))!*thick_frac

            END DO

          ENDIF
          END DO edge_idx_loop
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

    ELSEIF(PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)

          edge_idx_loop2: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,lev,jb) <= sea_boundary) THEN
            ictr = 0
            p_vn_e(je,jb) = 0.0_wp
            il_c        = p_patch%edges%cell_idx(je,jb,1)
            ib_c        = p_patch%edges%cell_blk(je,jb,1)
            scalar_cell = scalar(il_c,ib_c)

            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*scalar_cell &
              &  *(p_op_coeff%edge2edge_viacell_coeff(je,lev,jb,ictr))!*thick_frac

            END DO
            ictr        = no_primal_edges
            il_c        = p_patch%edges%cell_idx(je,jb,2)
            ib_c        = p_patch%edges%cell_blk(je,jb,2)
            scalar_cell = scalar(il_c,ib_c)
            !thick_cell  = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,1,ib_c)

            DO ie=1, no_primal_edges
              ictr =ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              !thick_edge=p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,1,ib_e)
              !thick_frac=thick_edge/thick_cell
              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*scalar_cell&
              &  *(p_op_coeff%edge2edge_viacell_coeff(je,lev,jb,ictr))!*thick_frac

            END DO
           ENDIF
           END DO edge_idx_loop2
      END DO
    ENDIF
  END SUBROUTINE map_edges2edges_viacell_3D_1lev
  !----------------------------------------------------------------------------- 
  SUBROUTINE map_edges2edges_viacell_3d_mlev_const_z( p_patch_3D, vn_e, p_op_coeff, p_vn_e,scalar)
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: p_patch_3D
    REAL(wp), INTENT(IN)                 :: vn_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(IN)   :: p_op_coeff
    REAL(wp), INTENT(INOUT)              :: p_vn_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(IN), OPTIONAL       :: scalar(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)  
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, jb, jk, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    all_edges => p_patch%edges%all
    slev = 1
    elev = n_zlev

    IF(.NOT.PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)

          edge_idx_loop_sfc: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,slev,jb) == sea) THEN

            p_vn_e(je,slev,jb) = 0.0_wp

            ictr = 0
            il_c = p_patch%edges%cell_idx(je,jb,1)
            ib_c = p_patch%edges%cell_blk(je,jb,1) 
            thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              p_vn_e(je,slev,jb) = p_vn_e(je,slev,jb) &
              &+vn_e(il_e,slev,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,slev,jb,ictr)&
              &*thick_cell

            END DO

            ictr = no_primal_edges
            il_c = p_patch%edges%cell_idx(je,jb,2)
            ib_c = p_patch%edges%cell_blk(je,jb,2)
            thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              p_vn_e(je,slev,jb) = p_vn_e(je,slev,jb) &
              &+vn_e(il_e,slev,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,slev,jb,ictr)&
              &*thick_cell

            END DO
            ENDIF
          END DO edge_idx_loop_sfc
 
        level_loop_e: DO jk = slev+1, elev
          edge_idx_loop: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,jk,jb) == sea) THEN
            p_vn_e(je,jk,jb)= 0.0_wp

            ictr = 0
            il_c = p_patch%edges%cell_idx(je,jb,1)
            ib_c = p_patch%edges%cell_blk(je,jb,1)!thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

              p_vn_e(je,jk,jb) = p_vn_e(je,jk,jb) &
              &+vn_e(il_e,jk,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr)&
              &*thick_edge

            END DO

            ictr = no_primal_edges
            il_c = p_patch%edges%cell_idx(je,jb,2)
            ib_c = p_patch%edges%cell_blk(je,jb,2)!thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

              p_vn_e(je,jk,jb) = p_vn_e(je,jk,jb) &
              &+vn_e(il_e,jk,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr)&
              &*thick_edge

            END DO
            ENDIF
          END DO edge_idx_loop
        END DO level_loop_e
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

    ELSEIF(PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)

          edge_idx_loop_sfc2: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,slev,jb) == sea) THEN

            p_vn_e(je,slev,jb) = 0.0_wp

            ictr = 0
            il_c       = p_patch%edges%cell_idx(je,jb,1)
            ib_c       = p_patch%edges%cell_blk(je,jb,1) 
            thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)
            scalar_cell= scalar(il_c,slev,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              p_vn_e(je,slev,jb) = p_vn_e(je,slev,jb)&
              &+vn_e(il_e,slev,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,slev,jb,ictr)&
              &*thick_cell *scalar_cell

            END DO

            ictr = no_primal_edges
            il_c = p_patch%edges%cell_idx(je,jb,2)
            ib_c = p_patch%edges%cell_blk(je,jb,2)
            thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)
            scalar_cell = scalar(il_c,slev,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              p_vn_e(je,slev,jb) = p_vn_e(je,slev,jb) &
              &+vn_e(il_e,slev,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,slev,jb,ictr)&
              &*thick_cell *scalar_cell

            END DO
            ENDIF
          END DO edge_idx_loop_sfc2
 
        level_loop_e2: DO jk = slev+1, elev
          edge_idx_loop2: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,jk,jb) == sea) THEN
            p_vn_e(je,jk,jb)= 0.0_wp

            ictr = 0
            il_c        = p_patch%edges%cell_idx(je,jb,1)
            ib_c        = p_patch%edges%cell_blk(je,jb,1)!thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)
            scalar_cell = scalar(il_c,jk,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

              p_vn_e(je,jk,jb) = p_vn_e(je,jk,jb) &
              &+vn_e(il_e,jk,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr)&
              &*thick_edge*scalar_cell

            END DO

            ictr        = no_primal_edges
            il_c        = p_patch%edges%cell_idx(je,jb,2)
            ib_c        = p_patch%edges%cell_blk(je,jb,2)!thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,jk,ib_c)
            scalar_cell = scalar(il_c,jk,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

              p_vn_e(je,jk,jb) = p_vn_e(je,jk,jb) &
              &+vn_e(il_e,jk,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr)&
              &*thick_edge*scalar_cell

            END DO
            ENDIF
          END DO edge_idx_loop2
        END DO level_loop_e2
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block
    ENDIF
  END SUBROUTINE map_edges2edges_viacell_3D_mlev_const_z
  !----------------------------------------------------------------------------- 

  SUBROUTINE map_edges2edges_viacell_3d_1lev_const_z( p_patch_3D, vn_e, p_op_coeff, p_vn_e,scalar)!&   subset_range)
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)       :: p_patch_3D
    REAL(wp), INTENT(in)                       :: vn_e(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_operator_coeff), INTENT(in)         :: p_op_coeff
    REAL(wp), INTENT(INOUT)                    :: p_vn_e(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(IN),OPTIONAL              :: scalar(nproma,p_patch_3D%p_patch_2D(1)%nblks_c)  
    !Local variables
    INTEGER :: slev, elev
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: je, jb, jk, ie
    REAL(wp) :: scalar_cell
    REAL(wp) :: thick_edge, thick_cell, thick_frac
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------

    all_edges => p_patch%edges%all
    slev = 1
    elev = n_zlev

    IF(.NOT.PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)

          edge_idx_loop_sfc: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,slev,jb) == sea) THEN

            p_vn_e(je,jb) = 0.0_wp

            ictr = 0
            il_c = p_patch%edges%cell_idx(je,jb,1)
            ib_c = p_patch%edges%cell_blk(je,jb,1) 
            thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,slev,jb,ictr)&
              &*thick_cell

            END DO

            ictr = no_primal_edges
            il_c = p_patch%edges%cell_idx(je,jb,2)
            ib_c = p_patch%edges%cell_blk(je,jb,2)
            thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,slev,jb,ictr)&
              &*thick_cell

            END DO
            ENDIF
          END DO edge_idx_loop_sfc
 
        level_loop_e: DO jk = slev+1, elev
          edge_idx_loop: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,jk,jb) == sea) THEN

            ictr = 0
            il_c = p_patch%edges%cell_idx(je,jb,1)
            ib_c = p_patch%edges%cell_blk(je,jb,1)            

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr)&
              &*thick_edge

            END DO

            ictr = no_primal_edges
            il_c = p_patch%edges%cell_idx(je,jb,2)
            ib_c = p_patch%edges%cell_blk(je,jb,2)            

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr)&
              &*thick_edge

            END DO
            ENDIF
          END DO edge_idx_loop
        END DO level_loop_e
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block

    ELSEIF(PRESENT(scalar))THEN

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)

          edge_idx_loop_sfc2: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,slev,jb) == sea) THEN

            p_vn_e(je,jb) = 0.0_wp

            ictr = 0
            il_c = p_patch%edges%cell_idx(je,jb,1)
            ib_c = p_patch%edges%cell_blk(je,jb,1) 
            thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)
            scalar_cell= scalar(il_c,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,slev,jb,ictr)&
              &*thick_cell*scalar_cell

            END DO

            ictr = no_primal_edges
            il_c = p_patch%edges%cell_idx(je,jb,2)
            ib_c = p_patch%edges%cell_blk(je,jb,2)
            thick_cell = p_patch_3D%p_patch_1D(1)%prism_thick_c(il_c,slev,ib_c)
            scalar_cell= scalar(il_c,ib_c)
            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)

              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,slev,jb,ictr)&
              &*thick_cell*scalar_cell

            END DO
            ENDIF
          END DO edge_idx_loop_sfc2
 
        level_loop_e2: DO jk = slev+1, elev
          edge_idx_loop2: DO je =  i_startidx_e, i_endidx_e
          IF (p_patch_3D%lsm_e(je,jk,jb) == sea) THEN

            ictr = 0
            il_c = p_patch%edges%cell_idx(je,jb,1)
            ib_c = p_patch%edges%cell_blk(je,jb,1)            
            scalar_cell= scalar(il_c,ib_c)

            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr)&
              &*thick_edge*scalar_cell

            END DO

            ictr = no_primal_edges
            il_c = p_patch%edges%cell_idx(je,jb,2)
            ib_c = p_patch%edges%cell_blk(je,jb,2)            
            scalar_cell= scalar(il_c,ib_c)
            DO ie=1, no_primal_edges
              ictr = ictr+1 
              il_e = p_patch%cells%edge_idx(il_c,ib_c,ie)
              ib_e = p_patch%cells%edge_blk(il_c,ib_c,ie)
              thick_edge = p_patch_3D%p_patch_1D(1)%prism_thick_e(il_e,jk,ib_e)

              p_vn_e(je,jb) = p_vn_e(je,jb) &
              &+vn_e(il_e,ib_e)*p_op_coeff%edge2edge_viacell_coeff(je,jk,jb,ictr)&
              &*thick_edge*scalar_cell

            END DO
            ENDIF
          END DO edge_idx_loop2
        END DO level_loop_e2
      END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block
    ENDIF
  END SUBROUTINE map_edges2edges_viacell_3D_1lev_const_z
  !-------------------------------------------------------------------------

  !>
  !!
  !! mpi parallelized by LL, openmp corrected
  SUBROUTINE map_edges2edges_viavert_3D(p_patch_3D, vn, p_vn_dual,p_op_coeff, vort_flux)
    TYPE(t_patch_3D ),TARGET,INTENT(IN)        :: p_patch_3D
    REAL(wp), INTENT(INOUT)                    :: vn(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(INout)  :: p_vn_dual(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(IN)          :: p_op_coeff
    REAL(wp), INTENT(INOUT)                    :: vort_flux(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

    !Local variables
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: je, jk, jb
    INTEGER :: il_e, ib_e
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: ictr, neighbor, vertex_edge
    INTEGER :: il_v, ib_v
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_edges => p_patch%edges%all
    !-----------------------------------------------------------------------
    slev    = 1
    elev    = n_zlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx_e,i_endidx_e)
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)

      level_loop: DO jk = slev, elev

        edge_idx_loop: DO je =  i_startidx_e, i_endidx_e

          IF (p_patch_3D%lsm_e(je,jk,jb) == sea) THEN

            vort_flux(je,jk,jb) = 0.0_wp 

            DO neighbor=1,2   
              IF(neighbor==1) ictr = 0
              IF(neighbor==2) ictr = no_dual_edges

              il_v = p_patch%edges%vertex_idx(je,jb,neighbor)
              ib_v = p_patch%edges%vertex_blk(je,jb,neighbor) 

              DO vertex_edge=1, p_patch%verts%num_edges(il_v,ib_v)!no_dual_cell_edges

                ictr =ictr+1 

                il_e = p_patch%verts%edge_idx(il_v,ib_v,vertex_edge)
                ib_e = p_patch%verts%edge_blk(il_v,ib_v,vertex_edge)

                vort_flux(je,jk,jb) =  vort_flux(je,jk,jb)+vn(il_e,jk,ib_e)&
                                    &*p_op_coeff%edge2edge_viavert_coeff(je,jk,jb,ictr)
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
  !>
  !! Discrete mapping of cell-based vectors to edges on the primal grid.
  !!
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M (2010-11)
  !!  mpi parallelized by LL, result not synced
  SUBROUTINE map_cell2edges_3D_mlevels( p_patch_3D, p_vn_c, ptp_vn, p_op_coeff,&
                                   & opt_slev, opt_elev, subset_range )
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:,:)    ! input vector (nproma,n_zlev,nblks_c)
    REAL(wp), INTENT(out)                      :: ptp_vn(:,:,:)    ! output vector (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff)                     :: p_op_coeff
    INTEGER, INTENT(in), OPTIONAL              :: opt_slev        ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL              :: opt_elev        ! optional vertical end level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range

    !Local variables
    INTEGER :: slev, elev
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: je, jb, jk
    INTEGER :: il_c1,ib_c1, il_c2,ib_c2
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    edges_in_domain   => p_patch%edges%in_domain

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
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)

      level_loop_e: DO jk = slev, elev
        edge_idx_loop: DO je =  i_startidx_e, i_endidx_e

          IF(p_patch_3D%lsm_e(je,jk,jb) <= sea_boundary)THEN

            !Get indices of two adjacent triangles
            il_c1 = p_patch%edges%cell_idx(je,jb,1)
            ib_c1 = p_patch%edges%cell_blk(je,jb,1)
            il_c2 = p_patch%edges%cell_idx(je,jb,2)
            ib_c2 = p_patch%edges%cell_blk(je,jb,2)
            ptp_vn(je,jk,jb) =&
              & DOT_PRODUCT(p_vn_c(il_c1,jk,ib_c1)%x,&
              & p_op_coeff%edge2cell_coeff_cc_t(je,jk,jb,1)%x)&
              & +DOT_PRODUCT(p_vn_c(il_c2,jk,ib_c2)%x,&
              & p_op_coeff%edge2cell_coeff_cc_t(je,jk,jb,2)%x)

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
  SUBROUTINE map_cell2edges_3D_1level( p_patch_3D, p_vn_c, ptp_vn,p_op_coeff, level, subset_range )
    
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:)    ! input vector (nproma,n_zlev,nblks_c)
    REAL(wp), INTENT(out)                      :: ptp_vn(:,:)    ! output vector (nproma,n_zlev,nblks_e)
    TYPE(t_operator_coeff)                     :: p_op_coeff
    INTEGER, INTENT(in) :: level          ! vertical level
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range

    !Local variables
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: je, jb
    INTEGER :: il_c1,ib_c1, il_c2,ib_c2
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    edges_in_domain   => p_patch%edges%in_domain
    ptp_vn(:,:) = 0.0_wp

    ! calculation of transposed P^TPv from Pv (incart coord)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)          
      edge_idx_loop: DO je =  i_startidx_e, i_endidx_e

        IF(p_patch_3D%lsm_e(je,level,jb) <= sea_boundary)THEN
          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(je,jb,1)
          ib_c1 = p_patch%edges%cell_blk(je,jb,1)
          il_c2 = p_patch%edges%cell_idx(je,jb,2)
          ib_c2 = p_patch%edges%cell_blk(je,jb,2)
          ptp_vn(je,jb) =&
            & DOT_PRODUCT(p_vn_c(il_c1,ib_c1)%x,&
            & p_op_coeff%edge2cell_coeff_cc_t(je,level,jb,1)%x)&
            & +DOT_PRODUCT(p_vn_c(il_c2,ib_c2)%x,&
            & p_op_coeff%edge2cell_coeff_cc_t(je,level,jb,2)%x)
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
  END SUBROUTINE map_cell2edges_3D_1level
  !-------------------------------------------------------------------------
  !>
  !!
  !! mpi parallelized by LL, openmp corrected
  SUBROUTINE nonlinear_coriolis_3d_old(p_patch_3D, vn, p_vn_dual, h_e, vort_v, &
    & p_op_coeff, vort_flux)

    TYPE(t_patch_3D ),TARGET,INTENT(IN):: p_patch_3D
    REAL(wp), INTENT(inout)                   :: vn       (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates), INTENT(inout) :: p_vn_dual(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp), INTENT(in)                      :: h_e      (nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(inout)                   :: vort_v   (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_operator_coeff),INTENT(in)         :: p_op_coeff
    REAL(wp), INTENT(inout)                   :: vort_flux(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)

    !Local variables
    !
    INTEGER :: slev, elev     ! vertical start and end level
    INTEGER :: jk, jb, jev,je
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: il_v1, il_v2, ib_v1, ib_v2
    TYPE(t_cartesian_coordinates) :: u_v1_cc, u_v2_cc
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_edges => p_patch%edges%all
    !-----------------------------------------------------------------------
    slev         = 1
    elev         = n_zlev
    !CALL map_edges2vert_3d(p_patch_3D%p_patch_2D(1), vn, p_op_coeff%edge2vert_coeff_cc, &
    !   & p_vn_dual)

    CALL rot_vertex_ocean_3d( p_patch_3D, vn, p_vn_dual, p_op_coeff, vort_v)
    CALL sync_patch_array(SYNC_V, p_patch, vort_v)

! !$OMP PARALLEL
! !$OMP DO PRIVATE(jb,jk,je,i_startidx_e,i_endidx_e,il_v1,ib_v1,il_v2,ib_v2, u_v1_cc, u_v2_cc)
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)

      level_loop: DO jk = slev, elev

        edge_idx_loop: DO je =  i_startidx_e, i_endidx_e

          IF (p_patch_3D%lsm_e(je,jk,jb) == sea) THEN
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
            vort_flux(je,jk,jb) =( &
               & - DOT_PRODUCT(u_v2_cc%x, p_op_coeff%edge2vert_coeff_cc_t(je,jk,jb,2)%x)&
               & + DOT_PRODUCT(u_v1_cc%x, p_op_coeff%edge2vert_coeff_cc_t(je,jk,jb,1)%x))
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
! !   SUBROUTINE map_cell2edges_2d( p_patch_3D, p_vn_c, ptp_vn, p_op_coeff)
! !     
! !     TYPE(t_patch_3D ),TARGET, INTENT(IN):: p_patch_3D
! !     TYPE(t_cartesian_coordinates), INTENT(in)  :: p_vn_c(:,:)    ! input vector (nproma,n_zlev,nblks_c) 
! !     REAL(wp), INTENT(inout)                    :: ptp_vn(:,:)    ! output vector (nproma,n_zlev,nblks_e)
! !     TYPE(t_operator_coeff)                     :: p_op_coeff
! ! 
! !     !Local variables
! !     INTEGER :: i_startidx_e, i_endidx_e
! !     INTEGER :: il_c1, ib_c1, il_c2, ib_c2
! !     INTEGER :: je, jb
! !     TYPE(t_subset_range), POINTER :: edges_in_domain
! !     TYPE(t_patch), POINTER        :: p_patch 
! !     !-----------------------------------------------------------------------
! !     p_patch   => p_patch_3D%p_patch_2D(1)
! !     !-----------------------------------------------------------------------
! !     !CALL message (TRIM(routine), 'start')
! !     edges_in_domain   => p_patch%edges%in_domain
! ! 
! !     ! calculation of transposed P^TPv from Pv (incart coord)
! !     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
! ! 
! !       edge_idx_loop: DO je =  i_startidx_e, i_endidx_e
! ! 
! !         IF(p_patch_3D%lsm_e(je,1,jb) <= sea_boundary)THEN          
! !           !Get indices of two adjacent triangles
! !           il_c1 = p_patch%edges%cell_idx(je,jb,1)
! !           ib_c1 = p_patch%edges%cell_blk(je,jb,1)
! !           il_c2 = p_patch%edges%cell_idx(je,jb,2)
! !           ib_c2 = p_patch%edges%cell_blk(je,jb,2)
! ! 
! !           ptp_vn(je,jb) = &
! !             &  DOT_PRODUCT(p_vn_c(il_c1,ib_c1)%x, p_op_coeff%edge2cell_coeff_cc_t(je,1,jb,1)%x)&
! !             & +DOT_PRODUCT(p_vn_c(il_c2,ib_c2)%x, p_op_coeff%edge2cell_coeff_cc_t(je,1,jb,2)%x)
! !         ELSE
! !           ptp_vn(je,jb) = 0.0_wp
! !         ENDIF
! ! 
! !       END DO edge_idx_loop
! !     END DO ! jb = edges_in_domain%start_block, edges_in_domain%end_block
! ! 
! ! !     CALL sync_patch_array(SYNC_E, p_patch, ptp_vn)
! !     
! !   END SUBROUTINE map_cell2edges_2d
  !-----------------------------------------------------------------------------


 
END MODULE mo_scalar_product

