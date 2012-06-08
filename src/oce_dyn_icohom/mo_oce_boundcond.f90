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
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ocean_nml,          ONLY: idisc_scheme, iswm_oce, i_bc_veloc_top, i_bc_veloc_bot
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_run_config,         ONLY: dtime
  USE mo_exception,          ONLY: message
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_util_dbg_prnt,      ONLY: dbg_print
  USE mo_oce_state,          ONLY: t_hydro_ocean_state, v_base
  USE mo_scalar_product,     ONLY: map_edges2cell, map_cell2edges_2d,&
                                 & map_cell2edges
  USE mo_sea_ice_types,      ONLY: t_sfc_flx
  USE mo_oce_physics,        ONLY: t_ho_params
  USE mo_oce_math_operators, ONLY: grad_fd_norm_oce_2d, div_oce_3D
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates, gvec2cvec!, cvec2gvec
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_intp_rbf,           ONLY: rbf_vec_interpol_cell
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range
  USE mo_sync,               ONLY: SYNC_E, sync_patch_array
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  
  
  PUBLIC :: bot_bound_cond_horz_veloc
  PUBLIC :: top_bound_cond_horz_veloc
  PUBLIC :: bot_bound_cond_vert_veloc
  !PUBLIC :: top_bound_cond_vert_veloc
  
  PUBLIC :: top_bound_cond_tracer
  PUBLIC :: bot_bound_cond_tracer
  !PUBLIC :: update_ocean_forcing_CORE
  !PUBLIC :: update_ocean_surface_fluxes
  
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
  SUBROUTINE top_bound_cond_horz_veloc( p_patch, p_os, p_sfc_flx, &
    & top_bc_u_c, top_bc_v_c, top_bc_u_cc )
    !
    TYPE(t_patch), TARGET                      :: p_patch
    TYPE(t_hydro_ocean_state), INTENT(inout)   :: p_os            ! ocean state variable
    TYPE(t_sfc_flx)                            :: p_sfc_flx       ! external data
    REAL(wp)                                   :: top_bc_u_c(:,:) ! Top boundary condition
    REAL(wp)                                   :: top_bc_v_c(:,:) ! dim: (nproma,nblks_c)
    TYPE(t_cartesian_coordinates), INTENT(out) :: top_bc_u_cc(:,:)
    
    !Local variables
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp):: z_scale
    REAL(wp):: z_e(nproma,1,p_patch%nblks_e)
    
    TYPE(t_subset_range), POINTER :: all_cells
    
    !CHARACTER(len=max_char_length), PARAMETER :: &
    !& routine = ('mo_oce_boundcond:top_bound_cond_veloc')
    !-----------------------------------------------------------------------

    z_e(:,1,:) = 0.0_wp

    all_cells => p_patch%cells%all
    
    ! Modification of surface wind forcing according to surface boundary condition
    IF(iswm_oce == 1)THEN
      z_scale = v_base%del_zlev_m(1)*rho_ref
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

   !CASE (0)

   !  ! CALL message (TRIM(routine),'ZERO top velocity boundary conditions chosen')
   !  DO jb = all_cells%start_block, all_cells%end_block
   !    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
   !    DO jc = i_startidx_c, i_endidx_c
   !      top_bc_u_c(jc,jb)    =0.0_wp
   !      top_bc_v_c(jc,jb)    =0.0_wp
   !      top_bc_u_cc(jc,jb)%x =0.0_wp
   !    END DO
   !  END DO

    CASE (1) ! Forced by wind stress stored in p_sfc_flx

      ! CALL message (TRIM(routine),'(1) top velocity boundary condition: use surface wind stress')
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary)THEN
            top_bc_u_c(jc,jb)    = p_sfc_flx%forc_wind_u(jc,jb)/z_scale
            top_bc_v_c(jc,jb)    = p_sfc_flx%forc_wind_v(jc,jb)/z_scale
            top_bc_u_cc(jc,jb)%x = p_sfc_flx%forc_wind_cc(jc,jb)%x/z_scale
          !ELSE
          !  top_bc_u_c(jc,jb)    =0.0_wp
          !  top_bc_v_c(jc,jb)    =0.0_wp
          !  top_bc_u_cc(jc,jb)%x =0.0_wp
          ENDIF
        END DO
      END DO

    CASE (2) ! Forced by difference between wind velocity stored in p_sfc_flx and ocean velocity at top layer

      ! CALL message (TRIM(routine),'(2) top velocity boundary condition: use forc-u minus U(1) ')
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          !IF(v_base%lsm_oce_c(jc,1,jb) <= sea_boundary)THEN

          top_bc_u_c(jc,jb)    = ( p_sfc_flx%forc_wind_u(jc,jb)   &
            & - p_os%p_diag%u(jc,1,jb) ) / z_scale
          top_bc_v_c(jc,jb)    = ( p_sfc_flx%forc_wind_u(jc,jb)   &
            & - p_os%p_diag%u(jc,1,jb) ) / z_scale
          top_bc_u_cc(jc,jb)%x = ( p_sfc_flx%forc_wind_cc(jc,jb)%x &
            & - p_os%p_diag%p_vn(jc,1,jb)%x)/z_scale
          !ELSE
          ! top_bc_u_c(jc,jb)    =0.0_wp
          ! top_bc_v_c(jc,jb)    =0.0_wp
          ! top_bc_u_cc(jc,jb)%x =0.0_wp
          !ENDIF
        END DO
      END DO
    END SELECT
    ! LL: no sync rquired

    CALL map_cell2edges( p_patch, top_bc_u_cc, p_os%p_aux%bc_top_vn, level=1)
    CALL sync_patch_array(SYNC_E, p_patch, p_os%p_aux%bc_top_vn)

    !---------Debug Diagnostics-------------------------------------------
    idt_src=2  ! output print level (1-5, fix)
    CALL dbg_print('top bound.cond. u_c'         ,top_bc_u_c               ,str_module,idt_src)
    CALL dbg_print('top bound.cond. v_c'         ,top_bc_v_c               ,str_module,idt_src)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('top bound.cond. vn'          ,p_os%p_aux%bc_top_vn     ,str_module,idt_src)
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
  SUBROUTINE bot_bound_cond_horz_veloc( p_patch, p_os, p_phys_param, div_coeff)
    !
    TYPE(t_patch), TARGET, INTENT(in)        :: p_patch         ! patch on which computation is performed
    TYPE(t_hydro_ocean_state), INTENT(inout) :: p_os            ! ocean state variable
    TYPE(t_ho_params), INTENT(in)            :: p_phys_param    ! physical parameters
    REAL(wp), INTENT(in)                     :: div_coeff(:,:,:,:)
    
    ! Local variables
    INTEGER :: jb, je,jc
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: z_dolic, z_dolic_c1,z_dolic_c2
    REAL(wp) :: z_norm
    REAL(wp) :: z_e(nproma,1,p_patch%nblks_e)
    REAL(wp) :: z_depth(nproma,1,p_patch%nblks_e)
    REAL(wp) :: z_div_depth(nproma,1,p_patch%nblks_c)
    TYPE(t_cartesian_coordinates) :: z_grad_u(nproma,1,p_patch%nblks_e)
    
    TYPE(t_subset_range), POINTER :: all_cells, edges_in_domain
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_oce_boundcond:bot_bound_cond_veloc')
      
    !-----------------------------------------------------------------------
    all_cells       => p_patch%cells%all
    edges_in_domain => p_patch%edges%in_domain
        
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
          
          z_dolic = v_base%dolic_c(jc,jb)
          IF ( z_dolic > min_dolic ) THEN  ! wet points only
            
            z_norm  = SQRT(2.0_wp * p_os%p_diag%kin(jc,z_dolic,jb))
            
            p_os%p_aux%bc_bot_veloc_cc(jc,jb)%x = &
              & p_phys_param%bottom_drag_coeff * z_norm * p_os%p_diag%p_vn(jc,z_dolic,jb)%x
            
            p_os%p_aux%bc_bot_u(jc,jb) = &
              & p_phys_param%bottom_drag_coeff * z_norm * p_os%p_diag%u(jc,z_dolic,jb)
            p_os%p_aux%bc_bot_v(jc,jb) = &
              & p_phys_param%bottom_drag_coeff * z_norm * p_os%p_diag%v(jc,z_dolic,jb)
            
          END IF
        END DO
      END DO
      
      CALL map_cell2edges_2d( p_patch, p_os%p_aux%bc_bot_veloc_cc, p_os%p_aux%bc_bot_vn)
      CALL sync_patch_array(SYNC_E, p_patch, p_os%p_aux%bc_bot_v)
      
    CASE(2) !Bottom friction and topographic slope
      CALL message (TRIM(routine), &
        & 'TOPOGRAPHY_SLOPE bottom velocity boundary conditions not implemented yet')
      
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          
          z_dolic = v_base%dolic_c(jc,jb)
          IF ( z_dolic >= min_dolic ) THEN  ! wet points only
            
            z_norm  = SQRT(2.0_wp*p_os%p_diag%kin(jc,z_dolic,jb))
            
            p_os%p_aux%bc_bot_veloc_cc(jc,jb)%x =&
              & p_phys_param%bottom_drag_coeff*z_norm*p_os%p_diag%p_vn(jc,z_dolic,jb)%x
            
            !Only for RBF relevant: there should be an if, PK
            p_os%p_aux%bc_bot_u(jc,jb)=&
              & p_phys_param%bottom_drag_coeff*z_norm*p_os%p_diag%u(jc,z_dolic,jb)
            p_os%p_aux%bc_bot_v(jc,jb)=&
              & p_phys_param%bottom_drag_coeff*z_norm*p_os%p_diag%v(jc,z_dolic,jb)
            
          END IF
        END DO
      END DO
      
      z_depth(:,1,:)=p_os%p_diag%thick_e
      CALL div_oce_3d( z_depth, p_patch, div_coeff, z_div_depth, opt_slev=1,opt_elev=1 )

      
      ! LL: the whole loop seems to be doing nothing,
      !     no sync
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e
          
          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(je,jb,1)
          ib_c1 = p_patch%edges%cell_blk(je,jb,1)
          il_c2 = p_patch%edges%cell_idx(je,jb,2)
          ib_c2 = p_patch%edges%cell_blk(je,jb,2)
          
          z_dolic_c1 = v_base%dolic_c(il_c1,ib_c1)
          z_dolic_c2 = v_base%dolic_c(il_c2,ib_c2)

          ! LL: this is not used
          IF(z_dolic_c1 >= min_dolic .AND. z_dolic_c2 >= min_dolic) THEN
            z_grad_u(je,1,jb)%x = (p_os%p_diag%p_vn(il_c2,z_dolic_c2,ib_c2)%x &
              & - p_os%p_diag%p_vn(il_c1,z_dolic_c1,ib_c1)%x) &
              & / p_patch%edges%dual_edge_length(je,jb)
          ELSE
            z_grad_u(je,1,jb)%x = 0.0_wp
          ENDIF

          ! LL: this is not used
          z_e(je,1,jb) = &
            & DOT_PRODUCT(p_patch%edges%primal_cart_normal(je,jb)%x,z_grad_u(je,1,jb)%x)
          
        END DO
      END DO
      
      CALL map_cell2edges_2d( p_patch, p_os%p_aux%bc_bot_veloc_cc, p_os%p_aux%bc_bot_vn)
      CALL sync_patch_array(SYNC_E, p_patch, p_os%p_aux%bc_bot_v)
      
      !p_os%p_aux%bc_bot_vn(:,:) = p_os%p_aux%bc_bot_vn(:,:) - z_e(:,1,:)
      
    CASE default
      CALL message (TRIM(routine),'choosen wrong bottom velocity boundary conditions')
    END SELECT

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('bot bound.cond. u_c'         ,p_os%p_aux%bc_bot_u      ,str_module,idt_src)
    CALL dbg_print('bot bound.cond. v_c'         ,p_os%p_aux%bc_bot_v      ,str_module,idt_src)
    CALL dbg_print('bot bound.cond. vn'          ,p_os%p_aux%bc_bot_vn     ,str_module,idt_src)
    !---------------------------------------------------------------------
    
  END SUBROUTINE bot_bound_cond_horz_veloc
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Computes bottom boundary condition for vertical velocity.
  !! sbr calulates  Pu dot P (nabla H), this corresponds to
  !! continuous top boundary conditiopn u dot nabla H
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!  mpi parallelized LL
  !!
  SUBROUTINE bot_bound_cond_vert_veloc( p_patch, p_os, bot_bc_w )
    !
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch     !  patch on which computation is performed
    !
    ! Normal verlocity at edges
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    !DR TYPE(external_data), INTENT(in) :: p_ext_data  !< external data
    !
    ! Bottom boundary condition at cells
    REAL(wp), INTENT(inout)           :: bot_bc_w(:,:) ! dim: (nproma,nblks_c)
    !
    ! Local variables
    INTEGER :: jb, jc, je, i_dolic
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_e, i_endidx_e
    REAL(wp) :: z_grad_h(nproma,1,p_patch%nblks_e)
    INTEGER, DIMENSION(:,:,:),POINTER :: iidx, iblk
    INTEGER, DIMENSION(:,:),  POINTER :: p_dolic
    REAL(wp), DIMENSION(:),   POINTER :: p_bathy
    TYPE(t_cartesian_coordinates) :: z_grad_h_cc(nproma,1,p_patch%nblks_c)

    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells    
    !-----------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain
    all_cells => p_patch%cells%all
        
    bot_bc_w(:,:) = 0.0_wp
    z_grad_h_cc(nproma,1,p_patch%nblks_c)%x = 0.0_wp
    
    iidx      => p_patch%edges%cell_idx
    iblk      => p_patch%edges%cell_blk
    p_bathy   => v_base%zlev_m
    p_dolic   => v_base%dolic_c
    
    !----------------------------------------
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      DO je = i_startidx_e, i_endidx_e
      
        i_dolic = v_base%dolic_c(je,jb)
        IF ( v_base%lsm_oce_e(je,i_dolic,jb) <= sea ) THEN
          
          z_grad_h(je,1,jb) =  &
            & ( p_bathy(p_dolic(iidx(je,jb,2),iblk(je,jb,2))) &
            & -  p_bathy(p_dolic(iidx(je,jb,1),iblk(je,jb,1)))) &
            & * p_patch%edges%inv_dual_edge_length(je,jb)
        ELSE
          z_grad_h(je,1,jb) =  0.0_wp
        ENDIF
        
      ENDDO
    END DO
    CALL sync_patch_array(SYNC_E, p_patch, z_grad_h(:,:,:))
    !----------------------------------------
    
    !----------------------------------------
    CALL map_edges2cell( p_patch, &
      & z_grad_h,&
      & z_grad_h_cc,&
      & opt_slev=1, opt_elev=1)
    !----------------------------------------
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        !calulate  Pu dot P (nabla H), this corresponds to continuous top boundary condition u dot nabla H
        bot_bc_w(jc,jb) = -DOT_PRODUCT(z_grad_h_cc(jc,1,jb)%x,&
          & p_os%p_diag%p_vn(jc,1,jb)%x)
      END DO
    END DO
    
  END SUBROUTINE bot_bound_cond_vert_veloc
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Computes top boundary condition for vertical velocity.
  !! sbr calulates (h^(n+1)-h^n)/dt + Pu dot P (nabla h), this corresponds to
  !! continuous top boundary condition d_t h +u dot nabla h
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!   no-mpi parallelized
  !!
  SUBROUTINE top_bound_cond_vert_veloc( p_patch, p_os, top_bc_w, timestep, p_int )
    !
    TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    REAL(wp), INTENT(inout)                       :: top_bc_w(nproma,p_patch%nblks_c)
    INTEGER :: timestep
    TYPE(t_int_state),TARGET,INTENT(in), OPTIONAL :: p_int
    
    ! Local variables
    INTEGER :: jb, jc
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    REAL(wp) :: z_grad_h(nproma,1,p_patch%nblks_e)
    REAL(wp) :: z_u_times_gradh_c
    TYPE(t_cartesian_coordinates) :: z_grad_h_cc_vec(1:nproma,1,1:p_patch%nblks_c)
    REAL(wp) :: grad_h_u(1:nproma,1,1:p_patch%nblks_c)
    REAL(wp) :: grad_h_v(1:nproma,1,1:p_patch%nblks_c)
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !          & routine = ('mo_oce_boundcond:bot_bound_cond_veloc')
    !-----------------------------------------------------------------------
    rl_start = 1
    rl_end = min_rlcell
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,1)
    
    top_bc_w(:,:) = 0.0_wp
 !  z_grad_h_cc_vec(nproma,1,p_patch%nblks_c)%x(:) = 0.0_wp
    
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,&
        & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        z_grad_h_cc_vec(nproma,1,p_patch%nblks_c)%x(:) = 0.0_wp
      END DO
    END DO
    
    !calculate normal derivative of new height field
    CALL grad_fd_norm_oce_2d(p_os%p_prog(nnew(1))%h, &
      & p_patch,                 &
      & z_grad_h(:,1,:))
    CALL sync_patch_array(SYNC_E, p_patch, z_grad_h(:,1,:))        
    
    IF(idisc_scheme==1)THEN
      CALL map_edges2cell( p_patch,        &
        & z_grad_h,       &
        & z_grad_h_cc_vec,&
      !                         & p_os%p_diag%h_e,&
        & opt_slev=1,opt_elev=1 )
      
    ELSEIF(idisc_scheme==2)THEN
      
      CALL rbf_vec_interpol_cell( z_grad_h,&
        & p_patch,    &
        & p_int,      &
        & grad_h_u,   &
        & grad_h_v,   &
        & opt_slev=1, opt_elev=1)
      
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
          & rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          CALL gvec2cvec(grad_h_u(jc,1,jb),        &
            & grad_h_v(jc,1,jb),        &
            & p_patch%cells%center(jc,jb)%lon,&
            & p_patch%cells%center(jc,jb)%lat,&
            & z_grad_h_cc_vec(jc,1,jb)%x(1),&
            & z_grad_h_cc_vec(jc,1,jb)%x(2),&
            & z_grad_h_cc_vec(jc,1,jb)%x(3) )
          ! if(jb==900)then
          ! write(*,*)'top w',grad_h_u(jc,1,jb),grad_h_v(jc,1,jb),&
          ! &z_grad_h_cc_vec(jc,1,jb)%x
          ! endif
        END DO
      END DO
    ENDIF
    
    !CALL message (TRIM(routine),'ZERO bottom velocity boundary conditions chosen')
    IF(timestep>1)THEN
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,&
          & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          !calulate  Pu dot P (nabla h), this corresponds to continuous top boundary condition u dot nabla h
          !z_h_u(jc,1,jb)*z_u_v(jc,1,jb) + z_h_v(jc,1,jb)*z_v_v(jc,1,jb)
          z_u_times_gradh_c = DOT_PRODUCT(z_grad_h_cc_vec(jc,1,jb)%x,p_os%p_diag%p_vn(jc,1,jb)%x)
          
          top_bc_w(jc,jb) = (p_os%p_prog(nnew(1))%h(jc,jb) - p_os%p_prog(nold(1))%h(jc,jb))/dtime&
            & + z_u_times_gradh_c
          !write(*,*)'top bc W:',jc,jb,p_os%p_diag%p_vn(jc,1,jb)%x
          !p_os%p_prog(nnew(1))%h(jc,jb), p_os%p_prog(nold(1))%h(jc,jb),&
          !          & z_u_times_gradh_c p_diag%p_vn(jc,1,jb)%x
        END DO
      END DO
    ELSE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,&
          & i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          !calulate  Pu dot P (nabla h), this corresponds to continuous top boundary condition u dot nabla h
          !z_h_u(jc,1,jb)*z_u_v(jc,1,jb) + z_h_v(jc,1,jb)*z_v_v(jc,1,jb)
          
          top_bc_w(jc,jb) = DOT_PRODUCT(z_grad_h_cc_vec(jc,1,jb)%x,p_os%p_diag%p_vn(jc,1,jb)%x)
          
          !write(*,*)'top bc W:',jc,jb,top_bc_w(jc,jb)!p_os%p_diag%p_vn(jc,1,jb)%x
          !p_os%p_prog(nnew(1))%h(jc,jb), p_os%p_prog(nold(1))%h(jc,jb),&
          !          & z_u_times_gradh_c p_diag%p_vn(jc,1,jb)%x
        END DO
      END DO
      
    ENDIF
    !write(*,*)'MAX/MIN top boundary cond: w:', maxval(top_bc_w(1:nproma,1:p_patch%nblks_c))!,&
    !& minval(top_bc_w(1:nproma,1:p_patch%nblks_c))
  END SUBROUTINE top_bound_cond_vert_veloc
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Computes top boundary condition for tracer specified by tracer_id.
  !! d C/dz
  !!
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!  mpi parallelized LL (no sync required)
  SUBROUTINE top_bound_cond_tracer( p_patch, pstate_oce, tracer_id, p_sfc_flx, top_bc_tracer)
    
    TYPE(t_patch)    , TARGET, INTENT(in) :: p_patch             ! patch on which computation is performed
    TYPE(t_hydro_ocean_state), INTENT(in) :: pstate_oce          ! ocean state variable
    INTEGER, INTENT(in)                   :: tracer_id
    TYPE(t_sfc_flx), INTENT(in)           :: p_sfc_flx
    REAL(wp), INTENT(inout)               :: top_bc_tracer(:,:,:) !Top boundary condition at cells for all tracers
    !
    !Local variables
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp):: z_c(nproma,p_patch%nblks_c)

    TYPE(t_subset_range), POINTER :: all_cells
    
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_oce_boundcond:top_bound_cond_tracer')
    !-----------------------------------------------------------------------
    all_cells => p_patch%cells%all
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        top_bc_tracer(jc,jb, tracer_id) = p_sfc_flx%forc_tracer(jc,jb, tracer_id)
      END DO
    END DO

    !---------Debug Diagnostics-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    z_c(:,:)=top_bc_tracer(:,:,tracer_id)
    CALL dbg_print('top bound.cond.tracer'       ,z_c                      ,str_module,idt_src)
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
  SUBROUTINE bot_bound_cond_tracer( p_patch, pstate_oce, tracer_id, bot_bc_tracer)
    
    TYPE(t_patch)    , TARGET, INTENT(in) :: p_patch              ! patch on which computation is performed
    TYPE(t_hydro_ocean_state), INTENT(in) :: pstate_oce           ! ocean state variable
    INTEGER, INTENT(in)                   :: tracer_id
    REAL(wp), INTENT(out)                 :: bot_bc_tracer(:,:,:) !Bottom boundary condition at cells for all tracers
    
    !Local variables
    INTEGER :: jc, jb
    INTEGER :: i_startidx_c, i_endidx_c
    TYPE(t_subset_range), POINTER :: all_cells
    !-----------------------------------------------------------------------
    all_cells => p_patch%cells%all
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        bot_bc_tracer(jc,jb, tracer_id) = 0.0_wp
      END DO
    END DO
  END SUBROUTINE bot_bound_cond_tracer
  !-------------------------------------------------------------------------
  
END MODULE mo_oce_boundcond
