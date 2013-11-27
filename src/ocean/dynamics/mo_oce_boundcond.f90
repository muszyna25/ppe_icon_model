!>
!! Contains the implementation of the top and bottom ocean boundary conditions
!!
!! @author Peter Korn, MPI
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2010-04)
!!  Modified by Stephan Lorenz,     MPI-M (2010-07)
!!  methods used are mpi parallelized, LL
!!
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_oce_boundcond
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_physical_constants, ONLY: rho_ref
  USE mo_impl_constants,     ONLY: max_char_length, sea_boundary, sea, min_rlcell, min_dolic
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ocean_nml,          ONLY: iswm_oce, i_bc_veloc_top, i_bc_veloc_bot
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_run_config,         ONLY: dtime
  USE mo_exception,          ONLY: message, finish
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_oce_state,          ONLY: t_hydro_ocean_state
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_scalar_product,     ONLY: map_cell2edges_3D
  USE mo_sea_ice_types,      ONLY: t_sfc_flx
  USE mo_oce_physics,        ONLY: t_ho_params
  USE mo_oce_math_operators, ONLY: grad_fd_norm_oce_2d_3d, div_oce_3D
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, gvec2cvec
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range
  USE mo_sync,               ONLY: SYNC_E, sync_patch_array
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  
  PUBLIC :: bot_bound_cond_horz_veloc
  PUBLIC :: top_bound_cond_horz_veloc
  !PUBLIC :: bot_bound_cond_vert_veloc
  !PUBLIC :: top_bound_cond_vert_veloc
  
  PUBLIC :: top_bound_cond_tracer
  PUBLIC :: bot_bound_cond_tracer

  INTEGER, PARAMETER :: top=1
  CHARACTER(len=12)  :: str_module = 'oceBoundCond'  ! Output of module for 1 line debug
  INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug

CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !! Computes top boundary condition for horizontal velocity. This information
  !! is required for cell velocity, i.e. we have to prescribe values for
  !! du_c/dz and dv_c/dz.
  !!
  !! The forcing fluxes are provided by mo_ho_forcing
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn,         MPI-M (2010)
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !!  mpi parallelized LL
  !!
  SUBROUTINE top_bound_cond_horz_veloc( patch_3D, p_os, p_op_coeff, p_sfc_flx, &
    & top_bc_u_c, top_bc_v_c, top_bc_u_cc )
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3D
    TYPE(t_hydro_ocean_state), INTENT(inout)   :: p_os            ! ocean state variable
    TYPE(t_operator_coeff), INTENT(IN)         :: p_op_coeff
    TYPE(t_sfc_flx)                            :: p_sfc_flx       ! external data
    REAL(wp)                                   :: top_bc_u_c(:,:) ! Top boundary condition
    REAL(wp)                                   :: top_bc_v_c(:,:) ! dim: (nproma,alloc_cell_blocks)
    TYPE(t_cartesian_coordinates), INTENT(inout) :: top_bc_u_cc(:,:)
    
    !Local variables
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp):: z_scale(nproma,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_subset_range), POINTER :: all_cells   
    TYPE(t_patch), POINTER        :: patch_2D
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !& routine = ('mo_oce_boundcond:top_bound_cond_veloc')
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    all_cells => patch_2D%cells%all
    !-----------------------------------------------------------------------
    
    ! Modification of surface wind forcing according to surface boundary condition
    IF(iswm_oce == 1)THEN
      !z_scale(:,:) = v_base%del_zlev_m(1)*rho_ref
      !z_scale(:,:) = rho_ref*patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,1,:)
       z_scale(:,:) = rho_ref*p_os%p_diag%thick_c(:,:)
    ELSEIF(iswm_oce /= 1)THEN
      z_scale = rho_ref
    ENDIF

    ! set to zero (NAG)
    top_bc_u_c (:,:)      = 0.0_wp
    top_bc_v_c (:,:)      = 0.0_wp
    top_bc_u_cc(:,:)%x(1) = 0.0_wp
    top_bc_u_cc(:,:)%x(2) = 0.0_wp
    top_bc_u_cc(:,:)%x(3) = 0.0_wp

    SELECT CASE (i_bc_veloc_top)

   CASE (0)

   !  ! CALL message (TRIM(routine),'ZERO top velocity boundary conditions chosen')
     DO jb = all_cells%start_block, all_cells%end_block
       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
       DO jc = i_startidx_c, i_endidx_c
         top_bc_u_c(jc,jb)    =0.0_wp
         top_bc_v_c(jc,jb)    =0.0_wp
         top_bc_u_cc(jc,jb)%x =0.0_wp
       END DO
     END DO

    CASE (1) ! Forced by wind stress stored in p_sfc_flx

      ! CALL message (TRIM(routine),'(1) top velocity boundary condition: use surface wind stress')
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF(patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
            top_bc_u_c(jc,jb)    = p_sfc_flx%forc_wind_u(jc,jb)/z_scale(jc,jb)
            top_bc_v_c(jc,jb)    = p_sfc_flx%forc_wind_v(jc,jb)/z_scale(jc,jb)
            top_bc_u_cc(jc,jb)%x = p_sfc_flx%forc_wind_cc(jc,jb)%x/z_scale(jc,jb)
          ENDIF
        END DO
      END DO

    CASE (2) ! Forced by difference between wind velocity stored in p_sfc_flx and ocean velocity at top layer

      ! CALL message (TRIM(routine),'(2) top velocity boundary condition: use forc-u minus U(1) ')
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF(patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
          top_bc_u_c(jc,jb)    = ( p_sfc_flx%forc_wind_u(jc,jb)   &
            & - p_os%p_diag%u(jc,1,jb) ) / z_scale(jc,jb)
          top_bc_v_c(jc,jb)    = ( p_sfc_flx%forc_wind_u(jc,jb)   &
            & - p_os%p_diag%u(jc,1,jb) ) / z_scale(jc,jb)
          top_bc_u_cc(jc,jb)%x = ( p_sfc_flx%forc_wind_cc(jc,jb)%x &
            & - p_os%p_diag%p_vn(jc,1,jb)%x)/z_scale(jc,jb)
          ENDIF
        END DO
      END DO
    END SELECT

    CALL map_cell2edges_3D( patch_3D, top_bc_u_cc,p_os%p_aux%bc_top_vn,p_op_coeff,  level=1)
    CALL sync_patch_array(SYNC_E, patch_3D%p_patch_2D(1), p_os%p_aux%bc_top_vn)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('top bound.cond. u_c'         ,top_bc_u_c               ,str_module,idt_src, &
      in_subset=patch_2D%cells%owned)
    CALL dbg_print('top bound.cond. v_c'         ,top_bc_v_c               ,str_module,idt_src, &
      in_subset=patch_2D%cells%owned)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('top bound.cond. vn'          ,p_os%p_aux%bc_top_vn     ,str_module,idt_src, &
      in_subset=patch_2D%edges%owned)
    !---------------------------------------------------------------------
    
  END SUBROUTINE top_bound_cond_horz_veloc
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Computes bottom boundary condition for horizontal velocity.
  !!
  !! Computes bottom boundary condition for horizontal velocity. This information
  !! is required for the cell velocity vector, i.e. we have to prescribe values for
  !! du_c/dz and dv_c/dz.
  !!
  !!
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !! Modified by Stephan Lorenz,     MPI-M (2010-07)
  !!  mpi parallelized LL
  !!
  SUBROUTINE bot_bound_cond_horz_veloc( patch_3D, p_os, p_phys_param, p_op_coeff)
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3D
    TYPE(t_hydro_ocean_state), INTENT(inout) :: p_os            ! ocean state variable
    TYPE(t_ho_params), INTENT(in)            :: p_phys_param    ! physical parameters
    TYPE(t_operator_coeff), INTENT(IN)       :: p_op_coeff
    !REAL(wp), INTENT(in)                     :: div_coeff(:,:,:,:)
    
    ! Local variables
    INTEGER :: jb, je,jc
    !INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: i_startidx_c, i_endidx_c
    !INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: z_dolic!, z_dolic_c1,z_dolic_c2
    REAL(wp) :: z_norm
    TYPE(t_subset_range), POINTER :: all_cells, edges_in_domain
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_oce_boundcond:bot_bound_cond_veloc')
    TYPE(t_patch), POINTER        :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3D%p_patch_2D(1)
    !-----------------------------------------------------------------------
    all_cells       => patch_2D%cells%all
    edges_in_domain => patch_2D%edges%in_domain

    SELECT CASE (i_bc_veloc_bot)

    CASE(0)
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          p_os%p_aux%bc_bot_veloc_cc(jc,jb)%x = 0.0_wp
        END DO
      END DO

    CASE(1)!Bottom friction
 
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
          IF ( z_dolic > min_dolic ) THEN  ! wet points only

            z_norm  = SQRT(2.0_wp * p_os%p_diag%kin(jc,z_dolic,jb))

            p_os%p_aux%bc_bot_veloc_cc(jc,jb)%x = &
              & p_phys_param%bottom_drag_coeff * z_norm * p_os%p_diag%p_vn(jc,z_dolic,jb)%x

            !p_os%p_aux%bc_bot_u(jc,jb) = &
            !  & p_phys_param%bottom_drag_coeff * z_norm * p_os%p_diag%u(jc,z_dolic,jb)
            !p_os%p_aux%bc_bot_v(jc,jb) = &
            !  & p_phys_param%bottom_drag_coeff * z_norm * p_os%p_diag%v(jc,z_dolic,jb)

          END IF
        END DO
      END DO
      CALL map_cell2edges_3D( patch_3D, p_os%p_aux%bc_bot_veloc_cc,p_os%p_aux%bc_bot_vn,p_op_coeff,level=1)
      !CALL map_cell2edges_2d( patch_3D, p_os%p_aux%bc_bot_veloc_cc, p_os%p_aux%bc_bot_vn,p_op_coeff)
      CALL sync_patch_array(SYNC_E, patch_2D, p_os%p_aux%bc_bot_vn)

    CASE(2) !Bottom friction and topographic slope
      CALL message (TRIM(routine), &
        & 'TOPOGRAPHY_SLOPE bottom velocity boundary conditions not implemented yet')
      CALL finish (TRIM(routine), 'TOPOGRAPHY_SLOPE bottom velocity boundary conditions not implemented yet') 
! !       DO jb = all_cells%start_block, all_cells%end_block
! !         CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
! !         DO jc = i_startidx_c, i_endidx_c
! !           
! !           !z_dolic = v_base%dolic_c(jc,jb)
! !           z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
! !           IF ( z_dolic >= min_dolic ) THEN  ! wet points only
! !             
! !             z_norm  = SQRT(2.0_wp*p_os%p_diag%kin(jc,z_dolic,jb))
! !             
! !             p_os%p_aux%bc_bot_veloc_cc(jc,jb)%x =&
! !               & p_phys_param%bottom_drag_coeff*z_norm*p_os%p_diag%p_vn(jc,z_dolic,jb)%x
! !             
! !             !Only for RBF relevant: there should be an if, PK
! !             p_os%p_aux%bc_bot_u(jc,jb)=&
! !               & p_phys_param%bottom_drag_coeff*z_norm*p_os%p_diag%u(jc,z_dolic,jb)
! !             p_os%p_aux%bc_bot_v(jc,jb)=&
! !               & p_phys_param%bottom_drag_coeff*z_norm*p_os%p_diag%v(jc,z_dolic,jb)
! !             
! !           END IF
! !         END DO
! !       END DO
! !       
! !       !z_depth(:,1,:)=p_os%p_diag%thick_e
! !       z_depth(:,1,:)=patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_e(:,1,:)
! !       CALL div_oce_3d( z_depth, patch_2D, p_op_coeff%div_coeff, z_div_depth, opt_slev=1,opt_elev=1 )
! ! 
! !       
! !       ! LL: the whole loop seems to be doing nothing,
! !       !     no sync
! !       DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !         CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
! !         DO je = i_startidx_e, i_endidx_e
! !           
! !           !Get indices of two adjacent triangles
! !           il_c1 = patch_2D%edges%cell_idx(je,jb,1)
! !           ib_c1 = patch_2D%edges%cell_blk(je,jb,1)
! !           il_c2 = patch_2D%edges%cell_idx(je,jb,2)
! !           ib_c2 = patch_2D%edges%cell_blk(je,jb,2)
! !           
! !           !z_dolic_c1 = v_base%dolic_c(il_c1,ib_c1)
! !           !z_dolic_c2 = v_base%dolic_c(il_c2,ib_c2)
! ! 
! !           z_dolic_c1 = patch_3D%p_patch_1D(1)%dolic_c(il_c1,ib_c1)
! !           z_dolic_c2 = patch_3D%p_patch_1D(1)%dolic_c(il_c2,ib_c2)
! !           ! LL: this is not used
! !           IF(z_dolic_c1 >= min_dolic .AND. z_dolic_c2 >= min_dolic) THEN
! !             z_grad_u(je,1,jb)%x = (p_os%p_diag%p_vn(il_c2,z_dolic_c2,ib_c2)%x &
! !               & - p_os%p_diag%p_vn(il_c1,z_dolic_c1,ib_c1)%x) &
! !               & / patch_2D%edges%dual_edge_length(je,jb)
! !           ELSE
! !             z_grad_u(je,1,jb)%x = 0.0_wp
! !           ENDIF
! ! 
! !           ! LL: this is not used
! !           z_e(je,1,jb) = &
! !             & DOT_PRODUCT(patch_2D%edges%primal_cart_normal(je,jb)%x,z_grad_u(je,1,jb)%x)
! !           
! !         END DO
! !       END DO
! !       
! !       CALL map_cell2edges_2d(patch_3D, p_os%p_aux%bc_bot_veloc_cc, p_os%p_aux%bc_bot_vn,p_op_coeff)
! !       CALL sync_patch_array(SYNC_E, patch_2D, p_os%p_aux%bc_bot_v)
! !       
! !       !p_os%p_aux%bc_bot_vn(:,:) = p_os%p_aux%bc_bot_vn(:,:) - z_e(:,1,:)
! !       
    CASE default
      CALL message (TRIM(routine),'choosen wrong bottom velocity boundary conditions')
    END SELECT

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    !CALL dbg_print('bot bound.cond. u_c'         ,p_os%p_aux%bc_bot_u      ,str_module,idt_src)
    !CALL dbg_print('bot bound.cond. v_c'         ,p_os%p_aux%bc_bot_v      ,str_module,idt_src)
    CALL dbg_print('bot bound.cond. vn'          ,p_os%p_aux%bc_bot_vn     ,str_module,idt_src, &
      in_subset=patch_2D%edges%owned)
    !---------------------------------------------------------------------
    
  END SUBROUTINE bot_bound_cond_horz_veloc
  !-------------------------------------------------------------------------
  
! !   !-------------------------------------------------------------------------
! !   !>
! !   !! Computes bottom boundary condition for vertical velocity.
! !   !! sbr calulates  Pu dot P (nabla H), this corresponds to
! !   !! continuous top boundary conditiopn u dot nabla H
! !   !!
! !   !! @par Revision History
! !   !! Developed  by  Peter Korn, MPI-M (2010).
! !   !!  mpi parallelized LL
! !   !!
! !   SUBROUTINE bot_bound_cond_vert_veloc( patch_2D, patch_3D, p_os, bot_bc_w )
! !     !
! !     TYPE(t_patch), TARGET, INTENT(in) :: patch_2D     !  patch on which computation is performed
! !     TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3D
! !     !
! !     ! Normal verlocity at edges
! !     TYPE(t_hydro_ocean_state), TARGET :: p_os
! !     !DR TYPE(external_data), INTENT(in) :: p_ext_data  !< external data
! !     !
! !     ! Bottom boundary condition at cells
! !     REAL(wp), INTENT(inout)           :: bot_bc_w(:,:) ! dim: (nproma,alloc_cell_blocks)
! !     !
! !     ! Local variables
! !     INTEGER :: jb, jc, je, i_dolic
! !     INTEGER :: i_startidx_c, i_endidx_c
! !     INTEGER :: i_startidx_e, i_endidx_e
! !     REAL(wp) :: z_grad_h(nproma,1,patch_2D%nblks_e)
! !     INTEGER, DIMENSION(:,:,:),POINTER :: iidx, iblk
! !     INTEGER, DIMENSION(:,:),  POINTER :: p_dolic
! !     REAL(wp), DIMENSION(:,:,:),   POINTER :: p_bathy
! !     TYPE(t_cartesian_coordinates) :: z_grad_h_cc(nproma,1,patch_2D%alloc_cell_blocks)
! ! 
! !     TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells    
! !     !-----------------------------------------------------------------------
! !     edges_in_domain => patch_2D%edges%in_domain
! !     all_cells => patch_2D%cells%all
! !         
! !     bot_bc_w(:,:) = 0.0_wp
! !     z_grad_h_cc(nproma,1,patch_2D%alloc_cell_blocks)%x = 0.0_wp
! !     
! !     iidx      => patch_2D%edges%cell_idx
! !     iblk      => patch_2D%edges%cell_blk
! !     !p_bathy   => v_base%zlev_m
! !     !p_dolic   => v_base%dolic_c
! !   
! !     p_bathy   => patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c
! !     p_dolic   => patch_3D%p_patch_1D(1)%dolic_c
! ! 
! !     !----------------------------------------
! !     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
! !       DO je = i_startidx_e, i_endidx_e
! !       
! !         !i_dolic = v_base%dolic_e(je,jb)
! !         i_dolic =patch_3D%p_patch_1D(1)%dolic_e(je,jb)
! !         !IF ( v_base%lsm_e(je,i_dolic,jb) <= sea ) THEN
! !         IF(patch_3D%lsm_e(je,i_dolic,jb) <= sea_boundary)THEN
! !           z_grad_h(je,1,jb) =  &
! !             & ( p_bathy(iidx(je,jb,2),i_dolic,iblk(je,jb,2)) &
! !             & -  p_bathy(iidx(je,jb,1),i_dolic,iblk(je,jb,1))) &
! !             & * patch_2D%edges%inv_dual_edge_length(je,jb)
! !         ELSE
! !           z_grad_h(je,1,jb) =  0.0_wp
! !         ENDIF
! !         
! !       ENDDO
! !     END DO
! !     CALL sync_patch_array(SYNC_E, patch_2D, z_grad_h(:,:,:))
! !     !----------------------------------------
! !     
! !     !----------------------------------------
! ! !     CALL map_edges2cell( patch_2D, &
! ! !       & z_grad_h,&
! ! !       & z_grad_h_cc,&
! ! !       & opt_slev=1, opt_elev=1)
! ! !     CALL map_edges2cell_3D( patch_2D, &
! ! !       & z_grad_h,&
! ! !       & z_grad_h_cc,&
! ! !       & opt_slev=1, opt_elev=1)
! !     !----------------------------------------
! !     
! !     DO jb = all_cells%start_block, all_cells%end_block
! !       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
! !       DO jc = i_startidx_c, i_endidx_c
! !         !calulate  Pu dot P (nabla H), this corresponds to continuous top boundary condition u dot nabla H
! !         bot_bc_w(jc,jb) = -DOT_PRODUCT(z_grad_h_cc(jc,1,jb)%x,&
! !           & p_os%p_diag%p_vn(jc,1,jb)%x)
! !       END DO
! !     END DO
! !     
! !   END SUBROUTINE bot_bound_cond_vert_veloc
! !   !-------------------------------------------------------------------------
! !   
! !   !-------------------------------------------------------------------------
! !   !>
! !   !! Computes top boundary condition for vertical velocity.
! !   !! sbr calulates (h^(n+1)-h^n)/dt + Pu dot P (nabla h), this corresponds to
! !   !! continuous top boundary condition d_t h +u dot nabla h
! !   !!
! !   !! @par Revision History
! !   !! Developed  by  Peter Korn, MPI-M (2010).
! !   !!   no-mpi parallelized
! !   !!
! !   SUBROUTINE top_bound_cond_vert_veloc( patch_2D, p_os, top_bc_w, timestep)!, p_int )
! !     !
! !     TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
! !     TYPE(t_hydro_ocean_state), TARGET :: p_os
! !     REAL(wp),POINTER                  :: grad_coeff(:,:,:)
! !     REAL(wp), INTENT(inout)           :: top_bc_w(nproma,patch_2D%alloc_cell_blocks)
! !     INTEGER                           :: timestep
! !     !TYPE(t_int_state),TARGET,INTENT(in), OPTIONAL :: p_int
! !     
! !     ! Local variables
! !     INTEGER :: jb, jc
! !     INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
! !     INTEGER :: rl_start, rl_end
! !     REAL(wp) :: z_grad_h(nproma,1,patch_2D%nblks_e)
! !     REAL(wp) :: z_u_times_gradh_c
! !     TYPE(t_cartesian_coordinates) :: z_grad_h_cc_vec(1:nproma,1,1:patch_2D%alloc_cell_blocks)
! !     REAL(wp) :: grad_h_u(1:nproma,1,1:patch_2D%alloc_cell_blocks)
! !     REAL(wp) :: grad_h_v(1:nproma,1,1:patch_2D%alloc_cell_blocks)
! !     ! CHARACTER(len=max_char_length), PARAMETER :: &
! !     !          & routine = ('mo_oce_boundcond:bot_bound_cond_veloc')
! !     !-----------------------------------------------------------------------
! !     rl_start = 1
! !     rl_end = min_rlcell
! !     i_startblk = patch_2D%cells%start_blk(rl_start,1)
! !     i_endblk   = patch_2D%cells%end_blk(rl_end,1)
! !     
! !     top_bc_w(:,:) = 0.0_wp
! !  !  z_grad_h_cc_vec(nproma,1,patch_2D%alloc_cell_blocks)%x(:) = 0.0_wp
! !     
! !     DO jb = i_startblk, i_endblk
! !       CALL get_indices_c(patch_2D, jb, i_startblk, i_endblk,&
! !         & i_startidx, i_endidx, rl_start, rl_end)
! !       DO jc = i_startidx, i_endidx
! !         z_grad_h_cc_vec(nproma,1,patch_2D%alloc_cell_blocks)%x(:) = 0.0_wp
! !       END DO
! !     END DO
! !     
! !     !calculate normal derivative of new height field
! !     !CALL grad_fd_norm_oce_2d(p_os%p_prog(nnew(1))%h, &
! !     !  & patch_2D,                 &
! !     !  & z_grad_h(:,1,:))
! !      CALL grad_fd_norm_oce_2d_3d( p_os%p_prog(nnew(1))%h,&
! !                                 & patch_2D,               &
! !                                 & grad_coeff(:,1,:),     &
! !                                 & z_grad_h(:,1,:))
! ! ! CALL grad_fd_norm_oce_2d_3D( p_os%p_prog(nold(1))%h,     &
! ! !          &                  p_patch_horz,                  &
! ! !          &                  p_op_coeff%grad_coeff(:,1,:),  &
! ! !          &                  z_gradh_e(:,1,:))
! !     CALL sync_patch_array(SYNC_E, patch_2D, z_grad_h(:,1,:))
! !     
! !     IF(idisc_scheme==1)THEN
! ! !       CALL map_edges2cell( patch_2D,        &
! ! !         & z_grad_h,       &
! ! !         & z_grad_h_cc_vec,&
! ! !       !                         & p_os%p_diag%h_e,&
! ! !         & opt_slev=1,opt_elev=1 )
! !       
! ! ! !     ELSEIF(idisc_scheme==2)THEN
! ! ! !       
! ! ! !       CALL rbf_vec_interpol_cell( z_grad_h,&
! ! ! !         & patch_2D,    &
! ! ! !         & p_int,      &
! ! ! !         & grad_h_u,   &
! ! ! !         & grad_h_v,   &
! ! ! !         & opt_slev=1, opt_elev=1)
! ! ! !       
! ! ! !       DO jb = i_startblk, i_endblk
! ! ! !         CALL get_indices_c(patch_2D, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
! ! ! !           & rl_start, rl_end)
! ! ! !         DO jc = i_startidx, i_endidx
! ! ! !           CALL gvec2cvec(grad_h_u(jc,1,jb),        &
! ! ! !             & grad_h_v(jc,1,jb),        &
! ! ! !             & patch_2D%cells%center(jc,jb)%lon,&
! ! ! !             & patch_2D%cells%center(jc,jb)%lat,&
! ! ! !             & z_grad_h_cc_vec(jc,1,jb)%x(1),&
! ! ! !             & z_grad_h_cc_vec(jc,1,jb)%x(2),&
! ! ! !             & z_grad_h_cc_vec(jc,1,jb)%x(3) )
! ! ! !           ! if(jb==900)then
! ! ! !           ! write(*,*)'top w',grad_h_u(jc,1,jb),grad_h_v(jc,1,jb),&
! ! ! !           ! &z_grad_h_cc_vec(jc,1,jb)%x
! ! ! !           ! endif
! ! ! !         END DO
! ! ! !       END DO
! !     ENDIF
! !     
! !     !CALL message (TRIM(routine),'ZERO bottom velocity boundary conditions chosen')
! !     IF(timestep>1)THEN
! !       DO jb = i_startblk, i_endblk
! !         CALL get_indices_c(patch_2D, jb, i_startblk, i_endblk,&
! !           & i_startidx, i_endidx, rl_start, rl_end)
! !         DO jc = i_startidx, i_endidx
! !           !calulate  Pu dot P (nabla h), this corresponds to continuous top boundary condition u dot nabla h
! !           !z_h_u(jc,1,jb)*z_u_v(jc,1,jb) + z_h_v(jc,1,jb)*z_v_v(jc,1,jb)
! !           z_u_times_gradh_c = DOT_PRODUCT(z_grad_h_cc_vec(jc,1,jb)%x,p_os%p_diag%p_vn(jc,1,jb)%x)
! !           
! !           top_bc_w(jc,jb) = (p_os%p_prog(nnew(1))%h(jc,jb) - p_os%p_prog(nold(1))%h(jc,jb))/dtime&
! !             & + z_u_times_gradh_c
! !           !write(*,*)'top bc W:',jc,jb,p_os%p_diag%p_vn(jc,1,jb)%x
! !           !p_os%p_prog(nnew(1))%h(jc,jb), p_os%p_prog(nold(1))%h(jc,jb),&
! !           !          & z_u_times_gradh_c p_diag%p_vn(jc,1,jb)%x
! !         END DO
! !       END DO
! !     ELSE
! !       DO jb = i_startblk, i_endblk
! !         CALL get_indices_c(patch_2D, jb, i_startblk, i_endblk,&
! !           & i_startidx, i_endidx, rl_start, rl_end)
! !         DO jc = i_startidx, i_endidx
! !           !calulate  Pu dot P (nabla h), this corresponds to continuous top boundary condition u dot nabla h
! !           !z_h_u(jc,1,jb)*z_u_v(jc,1,jb) + z_h_v(jc,1,jb)*z_v_v(jc,1,jb)
! !           
! !           top_bc_w(jc,jb) = DOT_PRODUCT(z_grad_h_cc_vec(jc,1,jb)%x,p_os%p_diag%p_vn(jc,1,jb)%x)
! !           
! !           !write(*,*)'top bc W:',jc,jb,top_bc_w(jc,jb)!p_os%p_diag%p_vn(jc,1,jb)%x
! !           !p_os%p_prog(nnew(1))%h(jc,jb), p_os%p_prog(nold(1))%h(jc,jb),&
! !           !          & z_u_times_gradh_c p_diag%p_vn(jc,1,jb)%x
! !         END DO
! !       END DO
! !       
! !     ENDIF
! !     !write(*,*)'MAX/MIN top boundary cond: w:', maxval(top_bc_w(1:nproma,1:patch_2D%nblks_c))!,&
! !     !& minval(top_bc_w(1:nproma,1:patch_2D%nblks_c))
! !   END SUBROUTINE top_bound_cond_vert_veloc
! !   !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Computes top boundary condition for tracer specified by tracer_id.
  !! d C/dz
  !!
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!  mpi parallelized LL (no sync required)
  SUBROUTINE top_bound_cond_tracer( patch_2D, pstate_oce, tracer_id, p_sfc_flx, top_bc_tracer)
    
    TYPE(t_patch)    , TARGET, INTENT(in) :: patch_2D             ! patch on which computation is performed
    TYPE(t_hydro_ocean_state), INTENT(in) :: pstate_oce          ! ocean state variable
    INTEGER, INTENT(in)                   :: tracer_id
    TYPE(t_sfc_flx), INTENT(in)           :: p_sfc_flx
    REAL(wp), INTENT(inout)               :: top_bc_tracer(:,:,:) !Top boundary condition at cells for all tracers
    !
    !Local variables
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp):: z_c(nproma,patch_2D%alloc_cell_blocks)

    TYPE(t_subset_range), POINTER :: all_cells
    
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_oce_boundcond:top_bound_cond_tracer')
    !-----------------------------------------------------------------------
    all_cells => patch_2D%cells%all
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        top_bc_tracer(jc,jb, tracer_id) = p_sfc_flx%forc_tracer(jc,jb, tracer_id)
      END DO
    END DO

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    z_c(:,:)=top_bc_tracer(:,:,tracer_id)
    CALL dbg_print('top bound.cond.tracer' ,z_c, str_module, idt_src, in_subset=patch_2D%cells%owned)
    !---------------------------------------------------------------------
    
  END SUBROUTINE top_bound_cond_tracer
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Computes bottom boundary condition for tracer specified by tracer_id.
  !!
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!  mpi parallelized LL (no sync required)
  !!
  SUBROUTINE bot_bound_cond_tracer( patch_2D, pstate_oce, tracer_id, bot_bc_tracer)
    
    TYPE(t_patch)    , TARGET, INTENT(in) :: patch_2D              ! patch on which computation is performed
    TYPE(t_hydro_ocean_state), INTENT(in) :: pstate_oce           ! ocean state variable
    INTEGER, INTENT(in)                   :: tracer_id
    REAL(wp), INTENT(inout)                 :: bot_bc_tracer(:,:,:) !Bottom boundary condition at cells for all tracers
    
    !Local variables
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    all_cells => patch_2D%cells%all
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        bot_bc_tracer(jc,jb, tracer_id) = 0.0_wp
      END DO
    END DO
  END SUBROUTINE bot_bound_cond_tracer
  !-------------------------------------------------------------------------
  
END MODULE mo_oce_boundcond
