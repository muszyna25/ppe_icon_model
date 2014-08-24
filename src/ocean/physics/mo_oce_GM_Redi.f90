!>
!! Contains the implementation of the tracer transport routines for the ICON ocean model.
!! This comprises advection and diffusion in horizontal and vertical direction.
!!
!!
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
! #include "omp_definitions.inc"
!----------------------------
MODULE mo_oce_GM_Redi
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_impl_constants,            ONLY: sea_boundary, sea, min_dolic
  USE mo_math_constants,            ONLY: pi, dbl_eps
  USE mo_physical_constants,  ONLY: grav, rho_ref, sal_ref, rho_inv, a_t, b_s, &
    & sitodbar, sfc_press_bar
  USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer,                    &
    &                                     implicit_diffusion,GMRedi_configuration,&
    &                                     GMRedi_combined, GM_only, Redi_only,Cartesian_Mixing,&
    &                                     k_tracer_dianeutral, k_tracer_isoneutral, k_tracer_GM_kappa,&
    &                                     implicit_diffusion,explicit_diffusion,vertical_tracer_diffusion_type,&
    &                                     GMRedi_configuration,GMRedi_combined, GM_only,Redi_only,Cartesian_Mixing, &
    &                                     tapering_scheme,tapering_DanaMcWilliams,tapering_Large,tapering_Griffies, &
    &                                     S_max, S_d, c_speed
    

  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
  USE mo_oce_types,                 ONLY: t_hydro_ocean_state, t_ocean_tracer !, v_base
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_oce_boundcond,             ONLY: top_bound_cond_tracer
  USE mo_oce_physics
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_c, sync_e, sync_patch_array
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_dif_vert
  USE mo_statistics,                ONLY: global_minmaxmean
  USE mo_mpi,                       ONLY: my_process_is_stdio !global_mpi_barrier
 
  USE mo_oce_math_operators,  ONLY: grad_fd_norm_oce_3d_onBlock, verticalDeriv_scalar_midlevel_on_block
  USE mo_scalar_product,            ONLY: map_cell2edges_3d,map_edges2cell_3d, &
    & map_edges2edges_viacell_3d_const_z, map_scalar_center2prismtop, map_scalar_prismtop2center
  IMPLICIT NONE
  
  PRIVATE
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'GM_Redi'
  CHARACTER(LEN=16)           :: str_module = 'GM_Redi'  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug

  PUBLIC  :: prepare_ocean_physics
  PUBLIC  :: calc_ocean_physics
  PRIVATE :: calc_flux_neutral_diffusion
  PRIVATE :: calc_flux_GentMcWilliams  
  PRIVATE :: calc_neutral_slopes
  PRIVATE :: calc_neutralslope_coeff
  PRIVATE :: calc_neutralslope_coeff_func 
  PRIVATE :: apply_tapering_function
  PRIVATE :: calc_tapering_function
  
CONTAINS


  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE prepare_ocean_physics(patch_3d, ocean_state, param, op_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(inout) :: op_coeff
   !-------------------------------------------------------------------------------

    CALL calc_neutral_slopes(patch_3d, ocean_state, param, op_coeff)

    CALL calc_tapering_function(patch_3d, ocean_state)

    CALL apply_tapering_function(patch_3d, ocean_state, param)

  END SUBROUTINE prepare_ocean_physics
  !-------------------------------------------------------------------------




  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE calc_ocean_physics(patch_3d, ocean_state, param, op_coeff, tracer_index)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params), INTENT(inout)         :: param
    TYPE(t_operator_coeff),INTENT(in)        :: op_coeff
    INTEGER, INTENT(IN)                      :: tracer_index
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, jc, jb
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain   
    
    
   SELECT CASE(GMRedi_configuration)!GMRedi_configuration==Cartesian_Mixing)RETURN

    CASE(Redi_only)
      CALL calc_flux_neutral_diffusion( patch_3d,    &
                                      & ocean_state, &
                                      & param,     &
                                      & op_coeff,  &
                                      & ocean_state%p_diag%GMRedi_flux_horz(:,:,:,tracer_index),&
                                      & ocean_state%p_diag%GMRedi_flux_vert(:,:,:,tracer_index),&
                                      & tracer_index)
    

    CASE(GM_only)
      CALL calc_flux_GentMcWilliams( patch_3d,    &
                                   & ocean_state, &
                                   & param,     &
                                   & op_coeff,  &
                                   & ocean_state%p_diag%GMRedi_flux_horz(:,:,:,tracer_index),&
                                   & ocean_state%p_diag%GMRedi_flux_vert(:,:,:,tracer_index),&
                                   & tracer_index)

    CASE(GMRedi_combined)
        CALL calc_combined_GentMcWilliamsRedi_flux( patch_3d,   &
                                                 & ocean_state,&
                                                 & param,      &
                                                 & op_coeff,   &
                                                 & ocean_state%p_diag%GMRedi_flux_horz(:,:,:,tracer_index),&
                                                 & ocean_state%p_diag%GMRedi_flux_vert(:,:,:,tracer_index),&
                                                 & tracer_index)
    CASE DEFAULT
    CALL finish(TRIM('mo_oce_GM_Redi'), 'This GMRedi_configuration is not supported')
    
    END SELECT

                                   

  END SUBROUTINE calc_ocean_physics
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE calc_flux_neutral_diffusion(patch_3d, ocean_state, param, op_coeff, redi_flux_horz, redi_flux_vert, tracer_index)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params),      INTENT(inout)    :: param
    TYPE(t_operator_coeff), INTENT(in)       :: op_coeff
    REAL(wp), INTENT(inout)                  :: redi_flux_horz(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(inout)                  :: redi_flux_vert(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_c)
    INTEGER, INTENT(IN)                      :: tracer_index  
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, jc, jb,jk,start_level,end_level, blockNo
    INTEGER :: start_edge_index, end_edge_index, je     
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain,edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: redi_flux_vert_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_c)
    TYPE(t_cartesian_coordinates) :: redi_flux_vec_horz_center(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
   
    TYPE(t_cartesian_coordinates),POINTER :: tracer_gradient_horz_vec_center(:,:,:), slopes(:,:,:)
    REAL(wp),POINTER :: tracer_gradient_vert_center(:,:,:)
    REAL(wp),POINTER :: K_I(:,:,:),K_D(:,:,:), slopes_squared(:,:,:)
    REAL(wp) tracer_gradient_vert_center2top(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    edges_in_domain => patch_2D%edges%in_domain 
    slopes          => ocean_state%p_aux%slopes 
    
    K_I           => param%k_tracer_isoneutral
    K_D           => param%k_tracer_dianeutral
    slopes_squared=> ocean_state%p_aux%slopes_squared

    start_level=1
    redi_flux_horz        (1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)=0.0_wp
    redi_flux_vert_center (1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_c)            =0.0_wp    
    !redi_flux_vert_top             (1:nproma, 1:n_zlev+1,1:patch_3D%p_patch_2d(1)%nblks_c)=0.0_wp   
    tracer_gradient_vert_center2top(1:nproma, 1:n_zlev+1,1:patch_3D%p_patch_2d(1)%nblks_c)=0.0_wp 

    IF(no_tracer<=2)THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index,jc, jk, end_level, &
!ICON_OMP z_adv_u_i) ICON_OMP_DEFAULT_SCHEDULE

      IF(tracer_index==1)THEN
        tracer_gradient_horz_vec_center =>ocean_state%p_aux%PgradTemperature_horz_center
        tracer_gradient_vert_center     => ocean_state%p_aux%DerivTemperature_vert_center
      ELSEIF(tracer_index==2)THEN
        tracer_gradient_horz_vec_center =>ocean_state%p_aux%PgradSalinity_horz_center
        tracer_gradient_vert_center     => ocean_state%p_aux%DerivSalinity_vert_center
      ENDIF
        
 Do jk=1,5
 write(*,*)'K_I',jk,maxval(K_I(:,jk,:)),minval(K_I(:,jk,:))
 write(*,*)'K_D',jk,maxval(K_I(:,jk,:)),minval(K_D(:,jk,:))
 END DO
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

        redi_flux_vert_center     (1:nproma,1:n_zlev,blockNo)     =0.0_wp
        redi_flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(1)=0.0_wp
        redi_flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(2)=0.0_wp
        redi_flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(3)=0.0_wp        
        
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
        DO jc = start_cell_index, end_cell_index
        
          DO jk = start_level, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)
          
            !horizontal Redi Flux
            redi_flux_vec_horz_center(jc,jk,blockNo)%x &
            &=K_I(jc,jk,blockNo)&          
            &*( tracer_gradient_horz_vec_center(jc,jk,blockNo)%x&
            &+tracer_gradient_vert_center(jc,jk,blockNo)*slopes(jc,jk,blockNo)%x)
          
!write(123,*)'GM h',jk,&
!&  tracer_gradient_horz_vec_center(jc,jk,blockNo)%x(1),tracer_gradient_vert_center(jc,jk,blockNo)*slopes(jc,jk,blockNo)%x(1)            
            !vertical Redi Flux
             redi_flux_vert_center(jc,jk,blockNo) &
             &= K_D(jc,jk,blockNo)&!+K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo))&
             &*tracer_gradient_vert_center(jc,jk,blockNo)&
             &+&              
             &K_I(jc,jk,blockNo)&               
             &*Dot_Product( tracer_gradient_horz_vec_center(jc,jk,blockNo)%x,slopes(jc,jk,blockNo)%x)
!             redi_flux_vert_center(jc,jk,blockNo) &
!             &= &              
!             &K_I(jc,jk,blockNo)&               
!             &*Dot_Product( tracer_gradient_horz_vec_center(jc,jk,blockNo)%x,slopes(jc,jk,blockNo)%x)
!write(123,*)'GM v',jk,&
!&  K_D(jc,jk,blockNo)*tracer_gradient_vert_center(jc,jk,blockNo),K_I(jc,jk,blockNo)&               
!             &*Dot_Product( tracer_gradient_horz_vec_center(jc,jk,blockNo)%x,slopes(jc,jk,blockNo)%x)            

          END DO                  
        END DO                
      END DO
 
! Do jk=1,5
! write(*,*)'vert deriv & slope',jk,&
! &maxval(tracer_gradient_vert_center(:,jk,:)),minval(tracer_gradient_vert_center(:,jk,:)),&
! &maxval(slopes(:,jk,:)%x(1)),minval(slopes(:,jk,:)%x(1)),&
! &maxval(tracer_gradient_vert_center(:,jk,:)*slopes(:,jk,:)%x(1)),&
! &minval(tracer_gradient_vert_center(:,jk,:)*slopes(:,jk,:)%x(1))
! !write(*,*)'K_D',jk,maxval(K_I(:,jk,:)),minval(K_D(:,jk,:))
! END DO

 
      CALL sync_patch_array(sync_c, patch_2D, redi_flux_vec_horz_center(:,:,:)%x(1))
      CALL sync_patch_array(sync_c, patch_2D, redi_flux_vec_horz_center(:,:,:)%x(2))
      CALL sync_patch_array(sync_c, patch_2D, redi_flux_vec_horz_center(:,:,:)%x(3))
      CALL sync_patch_array(sync_c, patch_2D, redi_flux_vert_center)

      IF(vertical_tracer_diffusion_type == explicit_diffusion)THEN
      
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
          DO jc = start_cell_index, end_cell_index
            DO jk = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)
         
              !vertical GM-Redi Flux
               redi_flux_vert_center(jc,jk,blockNo)   &
               &=redi_flux_vert_center(jc,jk,blockNo) &
               &+(K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo))*tracer_gradient_vert_center(jc,jk,blockNo)
!               redi_flux_vert_center(jc,jk,blockNo)   &
!               &=redi_flux_vert_center(jc,jk,blockNo) &
!               &+(K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo))*tracer_gradient_vert_center(jc,jk,blockNo)&
!               &+K_D(jc,jk,blockNo)&!+K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo))&
!               &*tracer_gradient_vert_center(jc,jk,blockNo)                       
            END DO                  
          END DO                
        END DO
        CALL sync_patch_array(sync_c, patch_2D, redi_flux_vert_center)          
        
      ELSEIF(vertical_tracer_diffusion_type == implicit_diffusion)THEN
      
        CALL map_scalar_center2prismtop( patch_3d, &
          &                              K_I(:,1:n_zlev,:)*slopes_squared(:,1:n_zlev,:),&
          &                              op_coeff,           &
          &                              param%a_tracer_v(:,:,:, tracer_index))  
     Do jk=1,n_zlev
     write(0,*)'New vert coeff',tracer_index,jk,maxval(param%a_tracer_v(:,jk,:, tracer_index)),&
     &minval(param%a_tracer_v(:,jk,:, tracer_index))
     END DO
!         DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
!       
!           CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
!       
!           DO jc = start_cell_index, end_cell_index
!             DO jk = start_level, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)-1
!               param%a_tracer_v(jc,jk,blockNo, tracer_index)&
!               &=0.5*(K_D(jc,jk,blockNo)*tracer_gradient_vert_center(jc,jk,blockNo)&
!               &    +K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo)&
!               &+K_D(jc,jk+1,blockNo)*tracer_gradient_vert_center(jc,jk+1,blockNo)&
!               &    +K_I(jc,jk+1,blockNo)*slopes_squared(jc,jk+1,blockNo))
!             
! !               write(1234,*)'vert coeff',jc,jk,blockNo, &
! !               &param%a_tracer_v(jc,jk,blockNo, tracer_index),K_I(jc,jk,blockNo),slopes_squared(jc,jk,blockNo),&
! !               &K_D(jc,jk,blockNo),tracer_gradient_vert_center(jc,jk,blockNo), &
! !               &K_D(jc,jk,blockNo)*tracer_gradient_vert_center(jc,jk,blockNo),K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo)        
!             END DO                  
!           END DO                
!         END DO          
      ENDIF
                
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('InRedi: vert_center',redi_flux_vert_center(:,2,:),&
      & str_module, idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------
       
      !map quantities to cell boundary
      CALL map_scalar_center2prismtop(patch_3d, redi_flux_vert_center, op_coeff,redi_flux_vert)
      CALL map_cell2edges_3D         ( patch_3D,redi_flux_vec_horz_center, redi_flux_horz(:,:,:),op_coeff)
        
!     Do jk=1,n_zlev
!     write(0,*)'mapping:Redicenter',tracer_index,jk,maxval(redi_flux_vert_center(:,jk,:)),minval(redi_flux_vert_center(:,jk,:)),&
!     &maxval(redi_flux_vert(:,jk,:)),minval(redi_flux_vert(:,jk,:))!,tracer_index))
! ! !   write(0,*)'mapping:gradient_center2top',jk,maxval(tracer_gradient_vert_center2top(:,jk,:)),&
!     END DO
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('InRedi: Redi_vert',redi_flux_vert(:,:,:),&!,tracer_index),&
      & str_module, idt_src, in_subset=cells_in_domain)
      CALL dbg_print('InRedi: Redi_horz',redi_flux_horz(:,:,:),&!,tracer_index),&
      & str_module, idt_src, in_subset=edges_in_domain)
      !---------------------------------------------------------------------

       
    !ICON_OMP_END_PARALLEL_DO    
     CALL sync_patch_array(sync_e, patch_2D, redi_flux_horz(:,:,:))
     CALL sync_patch_array(sync_c, patch_2D, redi_flux_vert(:,:,:))
     
  ELSEIF( no_tracer>2)THEN
    CALL finish(TRIM('calc_flux_neutral_diffusion'),&
    & 'calc_flux_neutral_diffusion beyond temperature and salinity is not impemented yet')
  ENDIF
  

  END SUBROUTINE calc_flux_neutral_diffusion
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the Gent-McWilliams eddy parametrization.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE calc_flux_GentMcWilliams(patch_3d, ocean_state, param, op_coeff, GM_flux_horz, GM_flux_vert, tracer_index)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params),      INTENT(inout)    :: param
    TYPE(t_operator_coeff), INTENT(in)       :: op_coeff
    REAL(wp), INTENT(inout)                  :: GM_flux_horz(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(inout)                  :: GM_flux_vert(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_c)
    INTEGER, INTENT(IN)                      :: tracer_index  
    
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index,start_edge_index, end_edge_index,start_level,end_level, blockNo
    INTEGER :: jc, je,jb,jk
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: GM_flux_vert_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_c)
    TYPE(t_cartesian_coordinates) :: GM_flux_vec_horz_center(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    
    TYPE(t_cartesian_coordinates),POINTER :: tracer_gradient_vec_horz(:,:,:), slopes(:,:,:)
    REAL(wp),POINTER :: tracer_gradient_vert(:,:,:)
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    edges_in_domain => patch_2D%edges%in_domain
    slopes          => ocean_state%p_aux%slopes  
    
    start_level=1
    
    IF(no_tracer<=2)THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index,jc, jk, end_level, &
!ICON_OMP z_adv_u_i) ICON_OMP_DEFAULT_SCHEDULE
      
      IF(tracer_index==1)THEN
        tracer_gradient_vec_horz => ocean_state%p_aux%PgradTemperature_horz_center
        tracer_gradient_vert     => ocean_state%p_aux%DerivTemperature_vert_center
      ELSEIF(tracer_index==2)THEN
        tracer_gradient_vec_horz => ocean_state%p_aux%PgradSalinity_horz_center
        tracer_gradient_vert     => ocean_state%p_aux%DerivSalinity_vert_center
      ENDIF
        
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
        
        GM_flux_vert_center(1:nproma,1:n_zlev,blockNo)         =0.0_wp
        GM_flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(1)=0.0_wp
        GM_flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(2)=0.0_wp
        GM_flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(3)=0.0_wp
        
        DO jc = start_cell_index, end_cell_index

          GM_flux_vert_center(jc,start_level,blockNo) = &
            &-param%k_tracer_GM_kappa(jc,start_level,blockNo)&
            &*Dot_Product(tracer_gradient_vec_horz(jc,start_level,blockNo)%x,slopes(jc,start_level,blockNo)%x)

        
          DO jk = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)

            !horizontal GM Flux
            GM_flux_vec_horz_center(jc,jk,blockNo)%x &
            &=param%k_tracer_GM_kappa(jc,jk,blockNo)* &
            & tracer_gradient_vert(jc,jk,blockNo)*slopes(jc,jk,blockNo)%x

              
            !vertical GM Flux
            GM_flux_vert_center(jc,jk,blockNo) = &
            &-param%k_tracer_GM_kappa(jc,jk,blockNo)&
            &*Dot_Product(tracer_gradient_vec_horz(jc,jk,blockNo)%x,slopes(jc,jk,blockNo)%x)
               
          END DO                  
        END DO          
      END DO
    !ICON_OMP_END_PARALLEL_DO

      CALL sync_patch_array(sync_c, patch_2D, GM_flux_vec_horz_center(:,:,:)%x(1))
      CALL sync_patch_array(sync_c, patch_2D, GM_flux_vec_horz_center(:,:,:)%x(2))
      CALL sync_patch_array(sync_c, patch_2D, GM_flux_vec_horz_center(:,:,:)%x(3))
      CALL sync_patch_array(sync_c, patch_2D, GM_flux_vert_center)

      !map quantities to cell boundary 
      CALL map_cell2edges_3D         ( patch_3d, GM_flux_vec_horz_center, GM_flux_horz, op_coeff)
      CALL map_scalar_center2prismtop( patch_3D, GM_flux_vert_center,op_coeff, GM_flux_vert)        

!        DO blockNo = all_edges%start_block, all_edges%end_block
!           CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
!           DO je = start_edge_index, end_edge_index
!             DO jk = start_level+1, patch_3D%p_patch_1D(1)%dolic_e(je,blockNo)
!               GM_flux_horz(je,jk,blockNo)&
!               & =param%k_tracer_GM_kappa_e(je,jk,blockNo)*GM_flux_horz(je,jk,blockNo)
!             END DO                  
!           END DO          
!         END DO        
!         DO blockNo = all_cells%start_block, all_cells%end_block
!           CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
!       
!           DO jc = start_cell_index, end_cell_index
!               DO jk = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)
!                 GM_flux_vert(jc,jk,blockNo)&
!                 & = -param%k_tracer_GM_kappa_e(je,jk,blockNo)*GM_flux_vert(jc,jk,blockNo)
!               END DO                  
!           END DO          
!         END DO
    
     Do jk=1,n_zlev
     write(0,*)'GM',tracer_index,jk,&
     &maxval(GM_flux_vert_center(:,jk,:)),minval(GM_flux_vert_center(:,jk,:)),&
     &maxval(GM_flux_vert(:,jk,:)),minval(GM_flux_vert(:,jk,:)),&
     &maxval(GM_flux_horz(:,jk,:)),minval(GM_flux_horz(:,jk,:))
     END DO
      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=3  ! output print level (1-5, fix)
      CALL dbg_print('InGM: GM_vert',GM_flux_vert(:,:,:),&!,tracer_index),&
      & str_module, idt_src, in_subset=cells_in_domain)
      CALL dbg_print('InGM: GM_horz',GM_flux_horz(:,:,:),&!,tracer_index),&
      & str_module, idt_src, in_subset=edges_in_domain)
      !---------------------------------------------------------------------
    ELSEIF( no_tracer>2)THEN
      CALL finish(TRIM('calc_flux_neutral_diffusion'),&
      & 'calc_flux_neutral_diffusion beyond temperature and salinityis not impemented yet')
  ENDIF
        
  END SUBROUTINE calc_flux_GentMcWilliams
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
   SUBROUTINE calc_combined_GentMcWilliamsRedi_flux(patch_3d, ocean_state, param, op_coeff, GMredi_flux_horz, GMredi_flux_vert,&
                                                  & tracer_index)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params),      INTENT(inout)    :: param
    TYPE(t_operator_coeff), INTENT(in)       :: op_coeff
    REAL(wp), INTENT(inout)                  :: GMredi_flux_horz(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp), INTENT(inout)                  :: GMredi_flux_vert(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_c) 
    INTEGER, INTENT(IN)                      :: tracer_index 
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, jc, jb,jk,start_level,end_level,blockNo
    INTEGER :: start_edge_index, end_edge_index, je     
    TYPE(t_subset_range), POINTER :: cells_in_domain,edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: flux_vert_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_c)
    TYPE(t_cartesian_coordinates) :: flux_vec_horz_center(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates),POINTER :: tracer_gradient_horz_vec_center(:,:,:), slopes(:,:,:)
    REAL(wp), POINTER :: tracer_gradient_vert_center(:,:,:),K_I(:,:,:), K_D(:,:,:), kappa(:,:,:), slopes_squared(:,:,:)  
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    edges_in_domain => patch_2D%edges%in_domain 
    slopes          => ocean_state%p_aux%slopes 

    K_I           => param%k_tracer_isoneutral
    K_D           => param%k_tracer_dianeutral
    kappa         => param%k_tracer_GM_kappa
    slopes_squared=> ocean_state%p_aux%slopes_squared
    
    start_level=1
    flux_vert_center    (1:nproma,1:n_zlev,  1:patch_3D%p_patch_2d(1)%nblks_c)=0.0_wp    
    flux_vec_horz_center(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_c)%x(1)=0.0_wp
    flux_vec_horz_center(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_c)%x(2)=0.0_wp
    flux_vec_horz_center(1:nproma,1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_c)%x(3)=0.0_wp        

  
    IF(no_tracer<=2)THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index,jc, jk, end_level, &
!ICON_OMP ) ICON_OMP_DEFAULT_SCHEDULE

      IF(tracer_index==1)THEN
        tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradTemperature_horz_center
        tracer_gradient_vert_center     => ocean_state%p_aux%DerivTemperature_vert_center
      ELSEIF(tracer_index==2)THEN
        tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradSalinity_horz_center
        tracer_gradient_vert_center     => ocean_state%p_aux%DerivSalinity_vert_center
      ENDIF
 Do jk=1,5
 write(*,*)'K_I',jk,maxval(K_I(:,jk,:)),minval(K_I(:,jk,:))
 write(*,*)'K_D',jk,maxval(K_I(:,jk,:)),minval(K_D(:,jk,:))
 END DO
        
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

        flux_vert_center    (1:nproma,1:n_zlev,blockNo)     =0.0_wp
        flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(1)=0.0_wp
        flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(2)=0.0_wp
        flux_vec_horz_center(1:nproma,1:n_zlev,blockNo)%x(3)=0.0_wp        
        
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
        DO jc = start_cell_index, end_cell_index
        
          !horizontal GMRedi Flux at top layer
          !flux_vec_horz_center(jc,start_level,blockNo)%x &
          !&=K_I(jc,start_level,blockNo) &
          !&*tracer_gradient_horz_vec_center(jc,start_level,blockNo)%x&
          !&+(K_I(jc,start_level,blockNo)-kappa(jc,start_level,blockNo))&
          !&*tracer_gradient_vert_center(jc,start_level,blockNo)*slopes(jc,start_level,blockNo)%x

          DO jk = start_level, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)
          
            !horizontal GM-Redi Flux
            flux_vec_horz_center(jc,jk,blockNo)%x &
            &=K_I(jc,jk,blockNo)* tracer_gradient_horz_vec_center(jc,jk,blockNo)%x!&
!            &+(K_I(jc,jk,blockNo)-kappa(jc,jk,blockNo))&
!            &*tracer_gradient_vert_center(jc,jk,blockNo)*slopes(jc,jk,blockNo)%x
              
            !vertical GM-Redi Flux
            flux_vert_center(jc,jk,blockNo) &
            &=&
            & K_D(jc,jk,blockNo)& !+K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo))&
            &*tracer_gradient_vert_center(jc,jk,blockNo)&
            &+&              
            &(K_I(jc,jk,blockNo)+kappa(jc,jk,blockNo))&               
            &*Dot_Product(tracer_gradient_horz_vec_center(jc,jk,blockNo)%x,slopes(jc,jk,blockNo)%x)
                
          END DO                  
        END DO                
      END DO
      CALL sync_patch_array(sync_c, patch_2D, flux_vec_horz_center(:,:,:)%x(1))
      CALL sync_patch_array(sync_c, patch_2D, flux_vec_horz_center(:,:,:)%x(2))
      CALL sync_patch_array(sync_c, patch_2D, flux_vec_horz_center(:,:,:)%x(3))
      CALL sync_patch_array(sync_c, patch_2D, flux_vert_center)

IF(tracer_index==1)THEN
   Do jk=1,n_zlev
   write(0,*)'GMR',jk,&
   &maxval(tracer_gradient_vert_center(:,jk,:)),minval(tracer_gradient_vert_center(:,jk,:)),&
   &maxval(slopes_squared(:,jk,:)),minval(slopes_squared(:,jk,:))
   END DO
 ENDIF


      IF(vertical_tracer_diffusion_type == explicit_diffusion)THEN
        
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
          DO jc = start_cell_index, end_cell_index
            DO jk = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)
              !vertical GM-Redi Flux
               flux_vert_center(jc,jk,blockNo)   &
               &=flux_vert_center(jc,jk,blockNo) &
               &+(K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo))*tracer_gradient_vert_center(jc,jk,blockNo)
               
!               flux_vert_center(jc,jk,blockNo)   &
!               &=flux_vert_center(jc,jk,blockNo) &
!               &+(K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo))*tracer_gradient_vert_center(jc,jk,blockNo)&
!               &+K_D(jc,jk,blockNo)&!+K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo))&
!               &*tracer_gradient_vert_center(jc,jk,blockNo)                       

            END DO                  
          END DO                
        END DO       
        CALL sync_patch_array(sync_c, patch_2D, flux_vert_center)          
      ELSEIF(vertical_tracer_diffusion_type == implicit_diffusion)THEN
      
        CALL map_scalar_center2prismtop( patch_3d, &
          &                              K_I(:,1:n_zlev,:)*slopes_squared(:,1:n_zlev,:),&
          &                              op_coeff,           &
          &                              param%a_tracer_v(:,:,:, tracer_index))  
     Do jk=1,n_zlev
     write(0,*)'New vert coeff',tracer_index,jk,maxval(param%a_tracer_v(:,jk,:, tracer_index)),&
     &minval(param%a_tracer_v(:,jk,:, tracer_index))
     END DO
! !         DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
! !       
! !           CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
! !       
! !           DO jc = start_cell_index, end_cell_index
! !             DO jk = start_level, patch_3D%p_patch_1D(1)%dolic_c(jc,blockNo)-1
! !               param%a_tracer_v(jc,jk,blockNo, tracer_index)&
! !               &=0.5*(K_D(jc,jk,blockNo)*tracer_gradient_vert_center(jc,jk,blockNo)&
! !               &    +K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo)&
! !               &+K_D(jc,jk+1,blockNo)*tracer_gradient_vert_center(jc,jk+1,blockNo)&
! !               &    +K_I(jc,jk+1,blockNo)*slopes_squared(jc,jk+1,blockNo))
! !             
! ! !               write(1234,*)'vert coeff',jc,jk,blockNo, &
! ! !               &param%a_tracer_v(jc,jk,blockNo, tracer_index),K_I(jc,jk,blockNo),slopes_squared(jc,jk,blockNo),&
! ! !               &K_D(jc,jk,blockNo),tracer_gradient_vert_center(jc,jk,blockNo), &
! ! !               &K_D(jc,jk,blockNo)*tracer_gradient_vert_center(jc,jk,blockNo),K_I(jc,jk,blockNo)*slopes_squared(jc,jk,blockNo)        
! !             END DO                  
! !           END DO                
! !         END DO

          
      ENDIF

        
   
      !map quantities to cell boundary
      CALL map_scalar_center2prismtop(patch_3d, flux_vert_center, op_coeff,GMredi_flux_vert)
      CALL map_cell2edges_3D         ( patch_3D,flux_vec_horz_center, GMredi_flux_horz, op_coeff)
       !---------DEBUG DIAGNOSTICS-------------------------------------------
       idt_src=3  ! output print level (1-5, fix)
       CALL dbg_print('InGMRedi: GMRedi_vert',GMredi_flux_vert(:,:,:),&
       & str_module, idt_src, in_subset=cells_in_domain)
       CALL dbg_print('InGMRedi: GMRedi_horz',GMredi_flux_horz(:,:,:),&
       & str_module, idt_src, in_subset=cells_in_domain)

       !---------------------------------------------------------------------
    !ICON_OMP_END_PARALLEL_DO    
    
     CALL sync_patch_array(sync_e, patch_2D, GMredi_flux_horz(:,:,:))
     CALL sync_patch_array(sync_c, patch_2D, GMredi_flux_vert(:,:,:))
     
    ELSEIF( no_tracer>2)THEN
      CALL finish(TRIM('calc_GMRediflux'),&
      & 'calc_flux_neutral_diffusion beyond temperature and salinity is not impemented yet')
  ENDIF
  
 
  Do jk=1,n_zlev
  write(0,*)'Horz/vert GMREDI flux',tracer_index,jk,&
  &maxval(GMredi_flux_horz(:,jk,:)),minval(GMredi_flux_horz(:,jk,:)),&
  &maxval(GMredi_flux_vert(:,jk,:)),minval(GMredi_flux_vert(:,jk,:))  
  END DO
  
!   Do jk=1,n_zlev
!   write(0,*)'VERT GMREDI flux',tracer_index,jk,maxval(GMredi_flux_vert(:,jk,:)),minval(GMredi_flux_vert(:,jk,:))
!   END DO


  END SUBROUTINE calc_combined_GentMcWilliamsRedi_flux
  !-------------------------------------------------------------------------
! 



  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the slopes required for isopycnial diffusion and eddy parametrization.
  !!
  !!    !Note that in case of 1-component fluid we use the density for the slope calculation,
  !!    !all relevant imformation is stored in the tracer%temperature structure
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE calc_neutral_slopes(patch_3d, ocean_state, param, op_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(inout) :: op_coeff
    
    !Local variables
    REAL(wp) :: grad_T_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: grad_S_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: grad_T_vert(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: grad_S_vert(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    
    TYPE(t_cartesian_coordinates),POINTER :: grad_T_vec(:,:,:)    
    TYPE(t_cartesian_coordinates),POINTER :: grad_S_vec(:,:,:)        
    REAL(wp),POINTER :: grad_T_vert_center(:,:,:)
    REAL(wp),POINTER :: grad_S_vert_center(:,:,:)
    
    REAL(wp) :: neutral_alpha(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: neutral_beta (nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    
    INTEGER :: jk, blockNo, je, jc,jb
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    INTEGER :: start_level, level
    REAL(wp) :: alpha, beta
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), POINTER :: pot_temp(:,:,:), salinity(:,:,:)
    INTEGER :: end_level
    REAL(wp):: size_grad_T_horz_vec(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: size_grad_S_horz_vec(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1) 
    all_cells       => patch_2D%cells%all
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain


    pot_temp          => ocean_state%p_prog(nold(1))%ocean_tracers(1)%concentration
    grad_T_vec        => ocean_state%p_aux%PgradTemperature_horz_center
    grad_T_vert_center=> ocean_state%p_aux%DerivTemperature_vert_center

    IF(no_tracer>=2)THEN
    
      pot_temp          => ocean_state%p_prog(nold(1))%ocean_tracers(1)%concentration
      grad_T_vec        => ocean_state%p_aux%PgradTemperature_horz_center
      grad_T_vert_center=> ocean_state%p_aux%DerivTemperature_vert_center
    
      salinity          => ocean_state%p_prog(nold(1))%ocean_tracers(2)%concentration
      grad_S_vec        => ocean_state%p_aux%PgradSalinity_horz_center
      grad_S_vert_center=> ocean_state%p_aux%DerivSalinity_vert_center
      
    !Note that in case of 1-component fluid we use the density for the slope calculation,
    !all relevant imformation is stored in the tracer%temperature structure
    ELSEIF(no_tracer==1)THEN      
      pot_temp          => ocean_state%p_diag%rho
      grad_T_vec        => ocean_state%p_aux%PgradTemperature_horz_center
      grad_T_vert_center=> ocean_state%p_aux%DerivTemperature_vert_center
    
    ENDIF

    
    start_level = 1
    
    grad_T_horz(1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)=0.0_wp
    grad_S_horz(1:nproma, 1:n_zlev,1:patch_3D%p_patch_2d(1)%nblks_e)=0.0_wp
    grad_T_vert(1:nproma, 1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    grad_S_vert(1:nproma, 1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    
    size_grad_T_horz_vec(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp
    size_grad_S_horz_vec(1:nproma,1:n_zlev,1:patch_3d%p_patch_2d(1)%alloc_cell_blocks)=0.0_wp    
   !-------------------------------------------------------------------------------    

    !1) calculation of horizontal and vertical gradient for potential temperature and salinity
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      
      !1a) calculate horizontal gradient of temperature
      CALL grad_fd_norm_oce_3d_onBlock ( &
        & pot_temp, &
        & patch_3D,                           &
        & op_coeff%grad_coeff(:,:,blockNo), &
        & grad_T_horz(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo)

     !1b)calculate horizontal  gradient of salinity
     IF(no_tracer>=2)THEN
      CALL grad_fd_norm_oce_3d_onBlock ( &
        & salinity, &
        & patch_3D,                    &
        & op_coeff%grad_coeff(:,:,blockNo), &
        & grad_S_horz(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo)
     ENDIF 
    END DO ! blocks
   CALL sync_patch_array(sync_e, patch_2D, grad_T_horz)
   IF(no_tracer>=2)   CALL sync_patch_array(sync_e, patch_2D, grad_S_horz)
   
   
   !---------DEBUG DIAGNOSTICS-------------------------------------------
   idt_src=3  ! output print level (1-5, fix)
   CALL dbg_print('calc_slopes: grad_T_horz',grad_T_horz,&
   &str_module,idt_src, in_subset=edges_in_domain)      
   CALL dbg_print('calc_slopes: grad_S_horz',grad_S_horz,&
   &str_module,idt_src, in_subset=edges_in_domain)      
  !---------------------------------------------------------------------  
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

      !1c) calculation of vertical derivative for temperature and salinity
      CALL verticalDeriv_scalar_midlevel_on_block( patch_3d,                &
                                                 & pot_temp(:,1:n_zlev,blockNo),   &
                                                 & grad_T_vert(:,1:n_zlev,blockNo),&
                                                 & start_level+1,             &
                                                 & blockNo,                 &
                                                 & start_cell_index,        &
                                                 & end_cell_index)

      IF(no_tracer>=2)THEN                                           
      CALL verticalDeriv_scalar_midlevel_on_block( patch_3d,                &
                                                 & salinity(:,:,blockNo),   &
                                                 & grad_S_vert(:,:,blockNo),&
                                                 & start_level+1,             &
                                                 & blockNo,                 &
                                                 & start_cell_index,        &
                                                 & end_cell_index)

      ENDIF
    END DO ! blocks
!ICON_OMP_END_PARALLEL_DO    

   CALL sync_patch_array(sync_c, patch_2D, grad_T_vert)       
   IF(no_tracer>=2)   CALL sync_patch_array(sync_c, patch_2D, grad_S_vert)   

   !---------DEBUG DIAGNOSTICS-------------------------------------------
!    idt_src=3  ! output print level (1-5, fix)
!    CALL dbg_print('neutral_slopes: grad_T_vert',grad_T_vert,&
!    &str_module,idt_src, in_subset=cells_in_domain)      
!    CALL dbg_print('neutral_slopes: grad_S_vert',grad_S_vert,&
!    &str_module,idt_src, in_subset=cells_in_domain)      
  !---------------------------------------------------------------------   
   
   
   
    !2) map horizontal and vertial derivative to cell centered vector
    CALL map_edges2cell_3d(patch_3D, &
        & grad_T_horz,               &
        & op_coeff,                &
        & ocean_state%p_aux%PgradTemperature_horz_center)
    CALL map_scalar_prismtop2center(patch_3d,&
        & grad_T_vert,                &
        & op_coeff,                 &
        &  ocean_state%p_aux%DerivTemperature_vert_center)

    CALL sync_patch_array(sync_c, patch_2D, grad_T_vec(:,:,:)%x(1))
    CALL sync_patch_array(sync_c, patch_2D, grad_T_vec(:,:,:)%x(2))
    CALL sync_patch_array(sync_c, patch_2D, grad_T_vec(:,:,:)%x(3))
       
    IF(no_tracer>=2)THEN        
      CALL map_edges2cell_3d(patch_3D, &
          & grad_S_horz,               &
          & op_coeff,                &
          & ocean_state%p_aux%PgradSalinity_horz_center)
      CALL map_scalar_prismtop2center(patch_3d,&
          & grad_S_vert,                &
          & op_coeff,                 &
          & ocean_state%p_aux%DerivSalinity_vert_center)
          
      CALL sync_patch_array(sync_c, patch_2D, grad_S_vec(:,:,:)%x(1))
      CALL sync_patch_array(sync_c, patch_2D, grad_S_vec(:,:,:)%x(2))
      CALL sync_patch_array(sync_c, patch_2D, grad_S_vec(:,:,:)%x(3))
          
    ENDIF

    CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_aux%DerivTemperature_vert_center)
    IF(no_tracer>=2)   CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_aux%DerivSalinity_vert_center)
    
    !4) calculate slope coefficients as thermal expanso and saline contraction coefficients
    CALL calc_neutralslope_coeff( patch_3d, &
                                & ocean_state%p_prog(nold(1))%tracer, &
                                & ocean_state%p_prog(nold(1))%h,      &
                                & neutral_alpha, neutral_beta)
                                
   CALL sync_patch_array(sync_c, patch_2D, neutral_alpha)
   CALL sync_patch_array(sync_c, patch_2D, neutral_beta)
   
  !---------DEBUG DIAGNOSTICS-------------------------------------------
   idt_src=3  ! output print level (1-5, fix)
   CALL dbg_print('calc_slopes: gradT_vert',grad_T_vert_center,&
   &str_module,idt_src, in_subset=cells_in_domain)      
   IF(no_tracer>=2)CALL dbg_print('calc_slopes: gradS_vert',grad_S_vert_center,&
   &str_module,idt_src, in_subset=cells_in_domain)      
  !---------------------------------------------------------------------   
    
    !5) calculate slope as cell centered vector 
    IF(no_tracer>=2)THEN    
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      
        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          IF(end_level >= min_dolic) THEN
        
            DO level = start_level+1, end_level-1
              size_grad_T_horz_vec(cell_index,level,blockNo)&
              &=SQRT(DOT_PRODUCT(grad_T_vec(cell_index,level,blockNo)%x,grad_T_vec(cell_index,level,blockNo)%x))
              size_grad_S_horz_vec(cell_index,level,blockNo)&
              &=SQRT(DOT_PRODUCT(grad_S_vec(cell_index,level,blockNo)%x,grad_S_vec(cell_index,level,blockNo)%x))
              
              alpha = -neutral_alpha(cell_index,level,blockNo)
              beta  =  neutral_beta(cell_index,level,blockNo)
            
               ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
               & = -(alpha*grad_T_vec(cell_index,level,blockNo)%x + beta*grad_S_vec(cell_index,level,blockNo)%x)&
               &/(alpha*grad_T_vert_center(cell_index,level,blockNo)+beta*grad_S_vert_center(cell_index,level,blockNo)-dbl_eps)

               ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)=&
               &DOT_PRODUCT(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,&
                           &ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)
!write(123,*)'slope squared',level,cell_index,blockNo,ocean_state%p_aux%slopes_squared(cell_index,level,blockNo),&
!&sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))
            END DO
          ENDIF
        END DO
      END DO
    ELSEIF(no_tracer==1)THEN
      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      
        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          IF(end_level >= min_dolic) THEN
        
            DO level = start_level+1, end_level-1
          
              alpha = -neutral_alpha(cell_index,level,blockNo)
              beta  = neutral_beta(cell_index,level,blockNo)
            
              ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
              & = -(alpha*grad_T_vec(cell_index,level,blockNo)%x)&
              &/(alpha*grad_T_vert_center(cell_index,level,blockNo)-dbl_eps)
            
             ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)=&
             &DOT_PRODUCT(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,&
                         &ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)

            END DO
          ENDIF
        END DO
      END DO
    
    ENDIF
!write(123,*)'--------------------------------------------------------------'    
    CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_aux%slopes(:,:,:)%x(1))
    CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_aux%slopes(:,:,:)%x(2))
    CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_aux%slopes(:,:,:)%x(3))
    CALL sync_patch_array(sync_c, patch_2D, ocean_state%p_aux%slopes_squared)


 !---------DEBUG DIAGNOSTICS-------------------------------------------
   idt_src=3  ! output print level (1-5, fix)
   CALL dbg_print('calc_slopes: sizegradT',size_grad_T_horz_vec,&
   &str_module,idt_src, in_subset=cells_in_domain)      
   IF(no_tracer>=2)CALL dbg_print('calc:slopes: sizegradS',size_grad_S_horz_vec,&
   &str_module,idt_src, in_subset=cells_in_domain)      
  !---------------------------------------------------------------------   

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=5  ! output print level (1-5, fix)
  DO level=1,n_zlev
  CALL dbg_print('calc_slopes: squared',(ocean_state%p_aux%slopes_squared(:,level,:)),&
  &str_module,idt_src, in_subset=cells_in_domain)
  END DO
  !---------------------------------------------------------------------

!  DO level= start_level+1, end_level-1  
!  write(*,*)'max-min SLOPE squared',level,maxval(ocean_state%p_aux%slopes_squared(:,level,:)),&
!  & minval(ocean_state%p_aux%slopes_squared(:,level,:))
!  END DO

  END SUBROUTINE calc_neutral_slopes
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the Gent-McWilliams eddy parametrization.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE calc_tapering_function(patch_3d, ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    !TYPE(t_ho_params),                 INTENT(inout) :: param
    !TYPE(t_operator_coeff),            INTENT(in)    :: op_coeff
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level, blockNo
    INTEGER :: start_edge_index, end_edge_index, je     
    TYPE(t_subset_range), POINTER :: cells_in_domain!all_cells, ,edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: lambda
    REAL(wp) :: depth_scale, depth
    REAL(wp) :: Coriolis_abs
    
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    !edges_in_domain => patch_2D%edges%in_domain 
    start_level=1
    !-------------------------------------------------------------------------------    
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          IF(end_level >= min_dolic) THEN
        
            DO level = start_level, end_level
              ocean_state%p_aux%taper_function_1(cell_index,level,blockNo) =&
              & 0.5_wp*(1.0_wp + tanh((S_max - sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)))/S_d))
              

            END DO
          ENDIF
        END DO
      END DO
      CALL sync_patch_array(sync_c, patch_2D,ocean_state%p_aux%taper_function_1)
! Do level=1,n_zlev      
! write(0,*)'max-min taper 1',maxval( ocean_state%p_aux%taper_function_1(:,level,:)),&
! &minval( ocean_state%p_aux%taper_function_1(:,level,:))     
! End do

      !tapering schemes other than Danabasoglu-McWilliams require a second
      !tapering function
      IF(tapering_scheme/=tapering_DanaMcWilliams)THEN
      
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

            IF(end_level >= min_dolic) THEN
        
              DO level = start_level, end_level-1
                Coriolis_abs=abs(patch_2d%cells%f_c(cell_index,blockNo))
                lambda=min(max(15000._WP, c_speed/Coriolis_abs),100000.0_WP)
                depth_scale = lambda*sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))
                depth=patch_3d%p_patch_1d(1)%depth_CellMiddle(cell_index,level,blockNo)

                IF(depth<=depth_scale)THEN
                  ocean_state%p_aux%taper_function_2(cell_index,level,blockNo) &
                  &= 0.5_wp*(1.0_wp+sin(pi*depth/depth_scale-pi/2.0_wp))
                ELSE
                  ocean_state%p_aux%taper_function_2(cell_index,level,blockNo) =1.0_wp
                ENDIF
! write(123,*)'tapering', level,cell_index,blockNo,ocean_state%p_aux%taper_function_1(cell_index,level,blockNo),&
! & ocean_state%p_aux%taper_function_2(cell_index,level,blockNo),&
! & ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)            
              END DO
            ENDIF
          END DO
        END DO
        CALL sync_patch_array(sync_c, patch_2D,ocean_state%p_aux%taper_function_2)
      ENDIF
!  Do level=1,n_zlev      
!  write(0,*)'max-min taper1/2',level,&
!   &maxval( ocean_state%p_aux%taper_function_1(:,level,:)),&
!   &minval( ocean_state%p_aux%taper_function_1(:,level,:)),& 
!  &maxval( ocean_state%p_aux%taper_function_2(:,level,:)),&
!  &minval( ocean_state%p_aux%taper_function_2(:,level,:)),& 
!  &maxval( ocean_state%p_aux%slopes_squared(:,level,:)),&
! &minval( ocean_state%p_aux%slopes_squared(:,level,:)) 
!  End do        
!  stop
     CALL dbg_print('calc_tapering_function: function 1', ocean_state%p_aux%taper_function_1 ,&
     & this_mod_name, 3, patch_2D%cells%in_domain)
     IF(tapering_scheme/=tapering_DanaMcWilliams)THEN
     CALL dbg_print('calc_tapering_function: function 2', ocean_state%p_aux%taper_function_2 ,&
     & this_mod_name, 3, patch_2D%cells%in_domain)
     ENDIF
     
   
   
    
  END SUBROUTINE calc_tapering_function
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE applies the tapering functions. Several options are available
  !! !
  !!         !A) Danabasoglu,G. and J. C.McWilliams, 1995:
  !!         ! Sensitivity of the global ocean circulation to 
  !!         ! parameterizations of mesoscale tracer transports
  !!         ! Journal of Climate, 8, 2967-2987
  !!         ! For steep slope regions:
  !!         !    Exponential taper applied to both the neutral and GM operator
  !!         !    parts.
  !!
  !!         !B) Large, W.G. et al. 1997
  !!         ! Sensitivity to surface forcing and boundary layer mixing in a global 
  !!         !  ocean model: annual-mean climatology
  !!         ! JPO, 27, 2418-2447
  !!         ! For steep slope regions:
  !!         !    Exponential taper applied to both neutral operator and GM operator
  !!         ! For near surface part:
  !!         !    Sine taper also applied to both neutral operator and GM operator. 
  !!         !C) Grffies, Fundamentals of Ocean Climte Models, 2004
  !!         ! For steep slope region:
  !!         ! a) no taper applied to diagonal piece of horizontal neutral operator
  !!         ! b) hyperbolic tangent(exponential) taper applied to off-diagonal piece of
  !!         !    horizontal operator and to diagonal and off-diagonal piece of vertical
  !!         !    neutral diffusion operator. a)+b) means we transfer the tracer diffusion
  !!         !    to a horizontal-vertical manner in regions of steep neutral slopes.
  !!         ! c) Exponential taper applied to GM operator.
  !!         ! For surface layer with small slope:
  !!         ! a) sine taper applied to both neutral operator and GM operator, except the
  !!         !    diagonal piece of the horizontal diffusion.
  
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  SUBROUTINE apply_tapering_function(patch_3d, ocean_state, param)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    !TYPE(t_operator_coeff),            INTENT(in)    :: op_coeff
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level, blockNo
    INTEGER :: start_edge_index, end_edge_index, je     
    TYPE(t_subset_range), POINTER :: cells_in_domain!,edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D    
    REAL(wp), POINTER :: K_I(:,:,:), K_D(:,:,:), kappa(:,:,:)
    REAL(wp)           :: geometric_scale, geometric_scale_factor_GMR
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    !edges_in_domain => patch_2D%edges%in_domain 
    start_level=1
    !-------------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------------
    K_I           => param%k_tracer_isoneutral
    K_D           => param%k_tracer_dianeutral
    kappa         => param%k_tracer_GM_kappa
    !-------------------------------------------------------------------------------
    
    
      SELECT CASE(tapering_scheme)
      
      CASE(tapering_DanaMcWilliams)
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

            IF(end_level >= min_dolic) THEN
              geometric_scale=1.0_wp!sqrt(patch_2D%cells%area(cell_index,blockNo))*geometric_scale_factor_GMR
              DO level = start_level, end_level
                 K_I  (cell_index,level,blockNo)  &
               &=k_tracer_isoneutral&
               &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)*geometric_scale
                
                 K_D  (cell_index,level,blockNo)  &
               &=k_tracer_dianeutral&
               &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)*geometric_scale
                
                kappa(cell_index,level,blockNo)   &
               &= k_tracer_GM_kappa*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)*geometric_scale
              END DO
            ENDIF
          END DO
        END DO
      CASE(tapering_Large)
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

            IF(end_level >= min_dolic) THEN
              geometric_scale=1.0_wp!sqrt(patch_2D%cells%area(cell_index,blockNo))*geometric_scale_factor_GMR        
              DO level = start_level, end_level
                 K_I  (cell_index,level,blockNo)  &
               &=k_tracer_isoneutral&
               &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
               &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)*geometric_scale
                
                 K_D  (cell_index,level,blockNo)  &
               &=k_tracer_dianeutral&
               &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
               &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)*geometric_scale 
               
                kappa(cell_index,level,blockNo)   &
               &= k_tracer_GM_kappa &
               &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
               &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)*geometric_scale               
              END DO
            ENDIF
          END DO
        END DO
      
      CASE(tapering_Griffies)
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

            IF(end_level >= min_dolic) THEN
              geometric_scale=1.0_wp!sqrt(patch_2D%cells%area(cell_index,blockNo))*geometric_scale_factor_GMR                
              DO level = start_level, end_level
                 K_I  (cell_index,level,blockNo)  &
               &=k_tracer_isoneutral&
               &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
               &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)*geometric_scale                               
               
                 K_D  (cell_index,level,blockNo)  &
               &=k_tracer_isoneutral*geometric_scale!&
               !&ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)
                
                kappa(cell_index,level,blockNo)   &
               &= k_tracer_GM_kappa &
               &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
               &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)*geometric_scale                              
              END DO
            ENDIF
          END DO
        END DO
      
      END SELECT
      CALL sync_patch_array(sync_c, patch_2D,K_I)
      CALL sync_patch_array(sync_c, patch_2D,K_D)
      CALL sync_patch_array(sync_c, patch_2D,kappa)
      
! write(0,*)'geometric factor',&
! & maxval(sqrt(patch_2D%cells%area)*geometric_scale_factor_GMR ),&
! & minval(sqrt(patch_2D%cells%area)*geometric_scale_factor_GMR )
     CALL dbg_print('apply_tapering: K_I', K_I , this_mod_name, 3, patch_2D%cells%in_domain)
     CALL dbg_print('apply_tapering: K_D', K_D , this_mod_name, 3, patch_2D%cells%in_domain)
     CALL dbg_print('apply_tapering: Kappa', kappa , this_mod_name, 3, patch_2D%cells%in_domain)
   
  END SUBROUTINE apply_tapering_function
  !-------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! Calculates polynomial coefficients for thermal expansion and saline contraction
  !! matching the equation of state as in Gill, Atmosphere-Ocean Dynamics, Appendix 3
  !!
  !! @par Revision History
  !! Initial version by Stephan Lorenz, MPI-M (2014)
  !!
  SUBROUTINE calc_neutralslope_coeff(patch_3d, tracer, surface_elevation, neutral_alph, neutral_beta)
    !
    !-----------------------------------------------------------------
    ! REFERENCE:
    !    McDougall, T.J. 1987.  Neutral Surfaces
    !    Journal of Physical Oceanography, vol 17, 1950-1964,
    !-----------------------------------------------------------------

    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3d
    REAL(wp), INTENT(in)                   :: tracer(:,:,:,:)         !  tracer(1): temperature, tracer(2): salinity
    REAL(wp), INTENT(in)                   :: surface_elevation(:,:)  !  surface elevation due to height equation
    REAL(wp), INTENT(inout)                :: neutral_alph(:,:,:)     !  thermal expansion coefficient [1/C]
    REAL(wp), INTENT(inout)                :: neutral_beta(:,:,:)     !  saline contraction coefficient [1/psu]
    
    ! !LOCAL VARIABLES:
    ! loop indices
    REAL(wp):: pressure, neutral_coeff(2)
    INTEGER :: jc, jk, jb
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, start_index, end_index
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D   => patch_3d%p_patch_2d(1)
    !-------------------------------------------------------------------------
    all_cells => patch_2D%cells%ALL

    !  tracer 1: potential temperature
    !  tracer 2: salinity
    IF(no_tracer==2)THEN
      
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk, pressure, neutral_coeff) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
          DO jk=1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb) ! operate on wet ocean points only
            ! compute pressure in dezi-bar, i.e. depth of water column in vertical centre (meter)
            !  - account for individual layer depth at bottom for use of partial cells (prism_thick_flat_sfc_c)
            !  - add elevation by passing old, new, or intermediate value of surface elevation (e.g. p_prog(nold(1)%h)
            pressure = patch_3d%p_patch_1d(1)%zlev_i(jk) &
              &      + patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(jc,jk,jb)*0.5_wp &
              &      + surface_elevation(jc,jb)
            neutral_coeff = calc_neutralslope_coeff_func( tracer(jc,jk,jb,1), tracer(jc,jk,jb,2), pressure)
            neutral_alph(jc,jk,jb) = neutral_coeff(1)
            neutral_beta(jc,jk,jb) = neutral_coeff(2)
          END DO
        END DO
      END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
      
    ELSEIF(no_tracer==1)THEN
      
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_index, end_index, jc, jk, pressure, neutral_coeff) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc = start_index, end_index
          DO jk=1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb) ! operate on wet ocean points only
            pressure = patch_3d%p_patch_1d(1)%zlev_i(jk) &
              &      + patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(jc,jk,jb)*0.5_wp &
              &      + surface_elevation(jc,jb)
            neutral_coeff = calc_neutralslope_coeff_func( tracer(jc,jk,jb,1), sal_ref, pressure)
            neutral_alph(jc,jk,jb) = neutral_coeff(1)
            neutral_beta(jc,jk,jb) = neutral_coeff(2)
          END DO
        END DO
      END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
      
    ENDIF

    CALL dbg_print('calc_neutral_coeff: alpha', neutral_alph , this_mod_name, 3, patch_2D%cells%in_domain)
    CALL dbg_print('calc_neutral_coeff: beta ', neutral_beta , this_mod_name, 3, patch_2D%cells%in_domain)

  END SUBROUTINE calc_neutralslope_coeff
 
 !-------------------------------------------------------------------------
  !>
  !! Calculates polynomial coefficients for thermal expansion and saline contraction
  !! matching the equation of state as described in (UNESCO)
  !!   Fofonoff and Millard, 1984, UNESCO, Paris, Tech. Pap. Mar. Sci., 44, 53pp
  !! This method is using the older !! IPTS (International Practical Temperature Scale) of 1968.
  !! The code below is adopted from FESOM (Quiang Wang, Sergey Danilov)
  !!
  !! @par Revision History
  !! Initial version by Stephan Lorenz, MPI-M (2014)
  !!
  FUNCTION calc_neutralslope_coeff_func(t,s,p) result(coeff)
    !
    !-----------------------------------------------------------------
    ! REFERENCES:
    !    McDougall, T.J. 1987.  Neutral Surfaces
    !    Journal of Physical Oceanography, Vol 17, 1950-1964,
    !-----------------------------------------------------------------
    ! CHECK VALUE:
    !    sw_beta=0.72088e-3 psu^-1 @ S=40.0psu, ptmp=10.0C (ITS-90), p=4000db
    !    a_over_b=0.34765 psu*C^-1 @ S=40.0psu, ptmp=10.0C, p=4000db
    ! Valid Range:
    !    S=25 to 40psu, p=0 to 4000db (ptmp=10C)
    !                   p=0 to 1000db (ptmp=20-40C)
    !-----------------------------------------------------------------
    !
    REAL(wp), INTENT(in)  :: t        !  potential temperature (in ITS-90) [C]
    REAL(wp), INTENT(in)  :: s        !  salinity (in PSS-78) [psu]
    REAL(wp), INTENT(in)  :: p        !  pressure (in dezi-bar) [db]
    REAL(wp)              :: coeff(2) !  thermal expansion [1/C] and saline contraction [1/psu] coefficients

    ! local variables, following the naming of the FESOM implementation
    REAL(wp):: aob, t1, t2, t3, t4, s35, s35sq, s1, s2, s3, p1, p2, p3
  
    !  polynomial parameter for calculation of saline contraction coeff beta
    REAL(wp), PARAMETER :: &
      & bet_t0   = 0.785567e-3_wp,  &
      & bet_t1   = 0.301985e-5_wp,  &
      & bet_t2   = 0.555579e-7_wp,  &
      & bet_t3   = 0.415613e-9_wp,  &
      & bet_st0  = 0.356603e-6_wp,  &
      & bet_st1  = 0.788212e-8_wp,  &
      & bet_sp1  = 0.408195e-10_wp, &
      & bet_sp2  = 0.602281e-15_wp, &
      & bet_s2   = 0.515032e-8_wp,  &
      & bet_p1t0 = 0.121555e-7_wp,  &
      & bet_p1t1 = 0.192867e-9_wp,  &
      & bet_p1t2 = 0.213127e-11_wp, &
      & bet_p2t0 = 0.176621e-12_wp, &
      & bet_p2t1 = 0.175379e-14_wp, &
      & bet_p3   = 0.121551e-17_wp
  
    !  polynomial parameter for calculation of thermal expansion coefficient alpha
    !  via fraction alpha over beta (aob)
    REAL(wp), PARAMETER :: &
      & aob_t0   = 0.665157e-1_wp,  &
      & aob_t1   = 0.170907e-1_wp,  &
      & aob_t2   = 0.203814e-3_wp,  &
      & aob_t3   = 0.298357e-5_wp,  &
      & aob_t4   = 0.255019e-7_wp,  &
      & aob_st0  = 0.378110e-2_wp,  &
      & aob_st1  = 0.846960e-4_wp,  &
      & aob_sp1  = 0.164759e-6_wp,  &
      & aob_sp2  = 0.251520e-11_wp, &
      & aob_s2   = 0.678662e-5_wp,  &
      & aob_p1t0 = 0.380374e-4_wp,  &
      & aob_p1t1 = 0.933746e-6_wp,  &
      & aob_p1t2 = 0.791325e-8_wp,  &
      & aob_p2t2 = 0.512857e-12_wp, &
      & aob_p3   = 0.302285e-13_wp

     t1 = t
     s1 = s
     p1 = p

   ! correction factor for conversion of 1990 to 1968 temperature standard (IPTS-68 to IPTS-90)
   ! the correction is less than 0.01 K in ocean water temperature range
   !  - T68 = 1.00024*T90
   !  - above mentioned CHECK VALUES of the paper are better met by this correction
     t1 = t*1.00024_wp
     
     t2    = t1*t1
     t3    = t2*t1
     t4    = t3*t1
     p2    = p1*p1
     p3    = p2*p1
     s35   = s-35.0_wp
     s35sq = s35*s35

     ! calculate beta, saline contraction
     coeff(2) = bet_t0 - bet_t1*t1                            &
       &         + bet_t2*t2 - bet_t3*t3                      &
       &         + s35*(-bet_st0    + bet_st1*t1              &
       &         +       bet_sp1*p1 - bet_sp2*p2)             &
       &         + s35sq*bet_s2                               & 
       &         + p1*(-bet_p1t0 + bet_p1t1*t1 - bet_p1t2*t2) &
       &         + p2*( bet_p2t0 - bet_p2t1*t1)               &
       &         + p3*bet_p3

     ! calculate alpha/beta
     aob      = aob_t0 + aob_t1*t1                            &
       &         - aob_t2*t2 + aob_t3*t3                      &
       &         - aob_t4*t4                                  &
       &         + s35*(+aob_st0    - aob_st1*t1              &
       &                -aob_sp1*p1 - aob_sp2*p2)             &
       &         - s35sq*aob_s2                               &
       &         + p1*(+aob_p1t0 - aob_p1t1*t1 + aob_p1t2*t2) &
       &         + p2*t2*aob_p2t2                             &
       &         - p3*aob_p3

     ! calculate alpha, thermal expansion
     coeff(1) = aob*coeff(2)
    
  END FUNCTION calc_neutralslope_coeff_func



END MODULE mo_oce_GM_Redi


