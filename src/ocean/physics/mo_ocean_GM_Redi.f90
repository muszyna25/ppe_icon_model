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
!=============================================================================================
#include "omp_definitions.inc"
#include "iconfor_dsl_definitions.inc"
!=============================================================================================
MODULE mo_ocean_GM_Redi
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  USE mo_impl_constants,            ONLY: sea_boundary, sea, min_dolic
  USE mo_math_constants,            ONLY: pi, dbl_eps
  USE mo_physical_constants,        ONLY: grav, sal_ref, rho_inv, a_t, b_s, &
    & sitodbar, sfc_press_bar
  USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer,                    &
    & GMRedi_configuration,&
    & GMRedi_combined, GM_only, Redi_only,Cartesian_Mixing,&
    & k_tracer_dianeutral_parameter, k_tracer_isoneutral_parameter, &
    & k_tracer_GM_kappa_parameter,EOS_TYPE,&
    & GMRedi_configuration,GMRedi_combined, GM_only,Redi_only,Cartesian_Mixing, &
    & tapering_scheme,tapering_DanaMcWilliams,tapering_Large,tapering_Griffies, &
    & S_max, S_d, S_critical, c_speed, GMRedi_usesRelativeMaxSlopes,            &
    & RossbyRadius_max, RossbyRadius_min,switch_off_diagonal_vert_expl,         &
    & GMREDI_COMBINED_DIAGNOSTIC,GM_INDIVIDUAL_DIAGNOSTIC,REDI_INDIVIDUAL_DIAGNOSTIC,&
    &TEST_MODE_REDI_ONLY,TEST_MODE_GM_ONLY,LinearThermoExpansionCoefficient, LinearHalineContractionCoefficient,&
    & SWITCH_OFF_TAPERING,SWITCH_ON_REDI_BALANCE_DIAGONSTIC, SWITCH_ON_TAPERING_HORIZONTAL_DIFFUSION,&
    &SLOPE_CALC_VIA_TEMPERTURE_SALINITY,BOLUS_VELOCITY_DIAGNOSTIC,REVERT_VERTICAL_RECON_AND_TRANSPOSED,&
    & INCLUDE_SLOPE_SQUARED_IMPLICIT

  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state, t_ocean_tracer !, v_base
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish, message !, message_text, message
  USE mo_ocean_boundcond,           ONLY: top_bound_cond_tracer
  USE mo_ocean_physics_types,       ONLY: t_ho_params 
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff, Get3DVectorTo2DLocal_array3D
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_patch_array_mult, sync_c, sync_e!, sync_patch_array
  USE mo_timer,                     ONLY: timer_start, timer_stop, timer_dif_vert
  USE mo_statistics,                ONLY: global_minmaxmean
  USE mo_mpi,                       ONLY: my_process_is_stdio !global_mpi_barrier
 
  USE mo_ocean_math_operators,      ONLY: grad_fd_norm_oce_3d_onBlock, verticalDeriv_scalar_onHalfLevels_on_block,&
    & verticalDeriv_vec_midlevel_on_block,div_oce_3d_ontriangles_onblock
  USE mo_scalar_product,            ONLY: map_cell2edges_3d,map_edges2cell_3d, &
    & map_scalar_center2prismtop, map_scalar_prismtop2center,map_edges2cell_with_height_3d,&
    & map_vec_prismtop2center_on_block,&
    & map_scalar_center2prismtop_GM, map_scalar_prismtop2center_GM,&
    & map_vec_prismtop2center_on_block_GM
    

  IMPLICIT NONE
  
  PRIVATE
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'GM_Redi'
  CHARACTER(LEN=16)           :: str_module = 'GM_Redi'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src    = 1          ! Level of detail for 1 line debug

  PUBLIC  :: prepare_ocean_physics
  PUBLIC  :: calc_ocean_physics
  PUBLIC  :: calc_neutralslope_coeff_func_onColumn
  
  PRIVATE :: calc_combined_GentMcWilliamsRedi_flux
  PRIVATE :: calc_neutral_slopes
  PRIVATE :: calc_tapering_function
  PRIVATE :: calc_entries_mixing_tensor
  PUBLIC  :: diagnose_Redi_flux_balance
  PRIVATE :: calc_bolus_velocity
  
CONTAINS


  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
!<Optimize:inUse:done>
  SUBROUTINE prepare_ocean_physics(patch_3d, ocean_state, param, op_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(inout) :: op_coeff
        
   !-------------------------------------------------------------------------------
   ! IF(tracer_index>2)THEN
   !   CALL message (TRIM('calc_ocean_physics'), 'GM-Redi parametrization &
   !   &requires only temperature and salinity, i.e.. tracers 1 and 2')
   !   return
   ! ENDIF


    CALL calc_neutral_slopes(patch_3d, ocean_state, param, op_coeff)

    CALL calc_tapering_function(patch_3d, param, ocean_state)
    
     
  END SUBROUTINE prepare_ocean_physics
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  !<Optimize:inUse:done>
  SUBROUTINE calc_ocean_physics(patch_3d, ocean_state, param, op_coeff, tracer_index)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params), INTENT(inout)         :: param
    TYPE(t_operator_coeff),INTENT(in)        :: op_coeff
    INTEGER, INTENT(IN)                      :: tracer_index
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index, blockNo
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    
    onEdges                :: GMredi_flux_horz
    onCells_HalfLevels     :: GMredi_flux_vert
    !-------------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain  => patch_2D%cells%in_domain   
    GMredi_flux_horz => ocean_state%p_diag%GMRedi_flux_horz(:,:,:,tracer_index)
    GMredi_flux_vert => ocean_state%p_diag%GMRedi_flux_vert(:,:,:,tracer_index)
        
    SELECT CASE(GMRedi_configuration)!GMRedi_configuration==Cartesian_Mixing)RETURN
 
    CASE(GMRedi_combined)
    
    
        CALL calc_combined_GentMcWilliamsRedi_flux( patch_3d,        &
            & ocean_state,      &
            & param,            &
            & op_coeff,         &
            & GMRedi_flux_horz, &
            & GMRedi_flux_vert, &
            & tracer_index)
    CASE DEFAULT
      CALL finish(TRIM('mo_ocean_GM_Redi'), 'This GMRedi_configuration is not supported')
    
    END SELECT
    
    IF(SWITCH_ON_REDI_BALANCE_DIAGONSTIC.AND.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN

     CALL diagnose_Redi_flux_balance(patch_3d, ocean_state, param, op_coeff,&
     & tracer_index)
    ENDIF           
    

  END SUBROUTINE calc_ocean_physics
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi and the eddy flux of GM.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  !<Optimize:inUse>
  SUBROUTINE calc_combined_GentMcWilliamsRedi_flux(patch_3d, ocean_state, param, op_coeff,&
    &GMredi_flux_horz, GMredi_flux_vert, tracer_index)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params),      INTENT(inout)    :: param
    TYPE(t_operator_coeff), INTENT(in)       :: op_coeff
    onEdges, INTENT(inout)                   :: GMredi_flux_horz
    onCells_HalfLevels, INTENT(inout)        :: GMredi_flux_vert
    INTEGER, INTENT(IN)                      :: tracer_index 
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level,blockNo
    INTEGER :: start_edge_index, end_edge_index, je     
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: flux_vert_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: flux_vec_horz_center(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates),POINTER :: tracer_gradient_horz_vec_center(:,:,:), slopes(:,:,:)
    REAL(wp), POINTER :: tracer_gradient_vert_center(:,:,:)

    TYPE(t_cartesian_coordinates)            :: taper_off_diagonal_vert(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)    
    TYPE(t_cartesian_coordinates)            :: taper_off_diagonal_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                                 :: taper_diagonal_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                                 :: taper_diagonal_vert_expl(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)    
    REAL(wp)                                 :: taper_diagonal_vert_impl(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)            
    REAL(wp) :: mapped_vertical_diagonal_impl(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)      
    
    REAL(wp) :: nabla_T_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: nabla_S_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)   
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    all_cells       => patch_2D%cells%all
    edges_in_domain => patch_2D%edges%in_domain 
    slopes          => ocean_state%p_aux%slopes 


    start_level=1
    IF(TEST_MODE_REDI_ONLY.OR.TEST_MODE_GM_ONLY)THEN
      CALL selective_GMRedi(patch_3d, ocean_state, param, &
                      & taper_diagonal_horz,           &
                      & taper_diagonal_vert_expl,      &
                      & taper_diagonal_vert_impl,      &
                      & taper_off_diagonal_horz,       &
                      & taper_off_diagonal_vert,       &
                      & tracer_index )
    ELSE
      CALL calc_entries_mixing_tensor(patch_3d, ocean_state, param, &
                      & taper_diagonal_horz,           &
                      & taper_diagonal_vert_expl,      &
                      & taper_diagonal_vert_impl,      &
                      & taper_off_diagonal_horz,       &
                      & taper_off_diagonal_vert,       &
                      & tracer_index )
    ENDIF    
     
    !Set pointers for tracer gradients, according to actual tracer index
    IF(tracer_index==1)THEN
      !write(0,*) "DerivTemperature_vec_center"
      tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradTemperature_horz_center
      tracer_gradient_vert_center     => ocean_state%p_aux%DerivTemperature_vert_center
    ELSEIF(tracer_index==2)THEN
      !write(0,*) "DerivSalinity_vec_center"
      tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradSalinity_horz_center
      tracer_gradient_vert_center     => ocean_state%p_aux%DerivSalinity_vert_center
    ELSEIF(tracer_index>2)THEN
    
      !Here we have to provide a sbr that calculates derivatives below    
      CALL calc_tracer_derivatives( patch_3d,&
                                  & ocean_state%p_prog(nold(1))%ocean_tracers(tracer_index)%concentration,&
                                  & ocean_state, &
                                  & op_coeff, &
                                  & tracer_index)

      tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradTracer_horz_center
      tracer_gradient_vert_center     => ocean_state%p_aux%DerivTracer_vert_center
       
    ENDIF

    DO blockNo = all_cells%start_block, all_cells%end_block        
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)      
      DO cell_index = start_cell_index, end_cell_index
        DO level = 1, n_zlev
          flux_vec_horz_center(cell_index,level,blockNo)%x = 0.0_wp
        ENDDO
      ENDDO
    ENDDO
      
    !The code below has to be executed for temperature, salinity and HAMMOC tracers, i.e. for all tracers.
    !
!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block        
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)      
      DO cell_index = start_cell_index, end_cell_index
 
        !GMRedi Flux at top layer: with tapering this is the just horizontal diffusion
        DO level = start_level, MIN(patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo),start_level)
          flux_vec_horz_center(cell_index,start_level,blockNo)%x &
          &=taper_diagonal_horz(cell_index,start_level,blockNo)  &
          &*tracer_gradient_horz_vec_center(cell_index,start_level,blockNo)%x

          !flux_vert_center(cell_index,start_level,blockNo) &
          !&=taper_diagonal_vert_expl(cell_index,start_level,blockNo)&
          !&*tracer_gradient_vert_center(cell_index,start_level,blockNo)

          ! the top level flux_vert_center will be filled from the second level, if it exists
          flux_vert_center(cell_index,start_level,blockNo) = 0.0_wp
          
        ENDDO

        DO level = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
          
          !horizontal GM-Redi Flux
          flux_vec_horz_center(cell_index,level,blockNo)%x &
          &=taper_diagonal_horz(cell_index,level,blockNo)*  &
          &tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x&
          !the second term vanishes if GM-Kappa=isoneutral diffusion
          &+taper_off_diagonal_horz(cell_index,level,blockNo)%x&
          &*tracer_gradient_vert_center(cell_index,level,blockNo)
              
          !vertical GM-Redi Flux: this is the part that is explicit in time
          !If namelist option "switch_off_diagonal_vert_expl=TRUE" the
          !first (diagonal) term is set to zero. Default option "switch_off_diagonal_vert_expl=FALSE"
          !Second term contains isoneutral and GM components
          flux_vert_center(cell_index,level,blockNo)= &
          &taper_diagonal_vert_expl(cell_index,level,blockNo)&
          &*tracer_gradient_vert_center(cell_index,level,blockNo)&
          &+&
          &Dot_Product(tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x,&
          &            taper_off_diagonal_vert(cell_index,level,blockNo)%x)
          
          !flux_vert_center(cell_index,level,blockNo)&
          !&= flux_vert_center(cell_index,level,blockNo)&
          !&*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(cell_index,level,blockNo)
                
        END DO

      END DO                
    END DO
!ICON_OMP_END_DO_PARALLEL


    !Map the explicit horizontal tracer flux from cell centers to edges (where the horizontal divergence is calculated)
    !
    ! use a vector communicator
    CALL sync_patch_array_mult(sync_c, patch_2D, 3, &
        & flux_vec_horz_center(:,:,:)%x(1), flux_vec_horz_center(:,:,:)%x(2), flux_vec_horz_center(:,:,:)%x(3))
    !    
    CALL map_cell2edges_3D( patch_3D,flux_vec_horz_center, GMredi_flux_horz, op_coeff)


    !Map the (explicit) vertical tracer flux to the prsim top (where the vertical divergence is calculated later)
    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_center2prismtop(patch_3d, flux_vert_center, op_coeff,GMredi_flux_vert)
    ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_center2prismtop_GM(patch_3d, flux_vert_center, op_coeff,GMredi_flux_vert)    
    ENDIF
    
    IF(INCLUDE_SLOPE_SQUARED_IMPLICIT)THEN
      ! Now we treat the vertical isoneutral coefficient that is discretized implicitely in time.  
      ! This is only neccessary once for temperature and salinity, the HAMOCC tracers use these value  
      IF(tracer_index<=2)THEN
        ! 
        !1.) Interpolate the tapered coefficient for the vertical tracer flux from prism center to prism top: 
        !this is the diagonal part that is handled implicitely in time
      
        IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN        
          CALL map_scalar_center2prismtop( patch_3d, &
          &                              taper_diagonal_vert_impl,&
          &                              op_coeff,                &
          &                              ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit)        

        ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
          CALL map_scalar_center2prismtop_GM( patch_3d, &
          &                              taper_diagonal_vert_impl,&
          &                              op_coeff,                &
          &                              ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit)        
    
        ENDIF      
        IF(tracer_index==1)THEN
          Do level=1,n_zlev
            !CALL dbg_print('Old vert coeff: A_v', param%a_tracer_v(:,level,:, tracer_index), this_mod_name, 4, patch_2D%cells%in_domain)
            CALL dbg_print('Old vert coeff: A_v', ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit,&
            & this_mod_name, 4, patch_2D%cells%in_domain)
          End do
        ENDIF
        !
        !2.) Here we combine the vertical GMRedicoefficient that is treated implicitely (mapped_vertical_diagonal_impl, this
        !term involves the slopes-squared times the isopycnal mixing coefficient and is potentially large, therefore
        !it is discretized implicitely) with the vertical mixing coefficient from the PP-scheme. 
        !We follow the approach in POP, where these two contributions are added
        !(see Reference manual POP, sect 5.1.3, in particular p. 41, after eq (150)).
        !
!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block     
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)      
          DO cell_index = start_cell_index, end_cell_index
            DO level = start_level, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
              param%a_tracer_v(cell_index,level,blockNo, tracer_index) = &!0.0_wp!  &
              !& param%a_tracer_v(cell_index,level,blockNo, tracer_index) + &
              & ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit(cell_index,level,blockNo)
            END DO                  
          END DO                
        END DO
!ICON_OMP_END_DO_PARALLEL
        IF(tracer_index==1)THEN
          Do level=1,n_zlev
            CALL dbg_print('New vert coeff: A_v', param%a_tracer_v(:,level,:, tracer_index),&
            & this_mod_name, 4, patch_2D%cells%in_domain)
            !CALL dbg_print('New vert coeff: A_v', &
            !&ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit(:,level,:),&
            !& this_mod_name, 4, patch_2D%cells%in_domain)
          END DO 
        ENDIF 
      ENDIF!IF(tracer_index<=1)THEN 
    ENDIF  
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    Do level=1,n_zlev
    idt_src=1  ! output print level (1-5, fix)
    CALL dbg_print('InGMRedi: GMRedi_vert',GMredi_flux_vert(:,level,:),&
    & str_module, idt_src, in_subset=cells_in_domain)
    END DO
    Do level=1,n_zlev
    CALL dbg_print('InGMRedi: GMRedi_horz',GMredi_flux_horz(:,level,:),&
    & str_module, idt_src, in_subset=edges_in_domain)
    END DO
    !---------------------------------------------------------------------
    
     


  END SUBROUTINE calc_combined_GentMcWilliamsRedi_flux
  !-------------------------------------------------------------------------
 

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the slopes required for isopycnial diffusion and eddy parametrization.
  !!
  !!    !Note that in case of 1-component fluid we use the density for the slope calculation,
  !!    !all relevant imformation is stored in the tracer%temperature structure
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_neutral_slopes(patch_3d, ocean_state, param, op_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(inout) :: op_coeff
    
    !Local variables
    REAL(wp) :: grad_T_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: grad_S_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: grad_T_vert(nproma, n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: grad_S_vert(nproma, n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: grad_rho_GM_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: grad_rho_GM_vert(nproma, n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
   
    TYPE(t_cartesian_coordinates),POINTER :: grad_T_vec(:,:,:)    
    TYPE(t_cartesian_coordinates),POINTER :: grad_S_vec(:,:,:)        
    TYPE(t_cartesian_coordinates),POINTER :: grad_rho_GM_vec(:,:,:)            
    REAL(wp),POINTER :: grad_T_vert_center(:,:,:)
    REAL(wp),POINTER :: grad_S_vert_center(:,:,:)
    REAL(wp),POINTER :: grad_rho_GM_vert_center(:,:,:)
       
    INTEGER :: level, blockNo, je, start_level, end_level
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    REAL(wp) :: alpha, beta
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), POINTER :: pot_temp(:,:,:), salinity(:,:,:),rho_GM(:,:,:)
!     REAL(wp):: size_grad_T_horz_vec(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     REAL(wp):: size_grad_S_horz_vec(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    REAL(wp), POINTER :: depth_cellinterface(:,:,:)
    REAL(wp):: neutral_coeff(1:n_zlev, 2), salinityColumn(1:n_zlev)
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    all_cells       => patch_2D%cells%all
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain
    depth_cellinterface => patch_3D%p_patch_1d(1)%depth_cellinterface
    !-------------------------------------------------------------------------

    IF(no_tracer>=2)THEN
      IF(SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
        pot_temp          => ocean_state%p_prog(nold(1))%ocean_tracers(1)%concentration
        grad_T_vec        => ocean_state%p_aux%PgradTemperature_horz_center
        grad_T_vert_center=> ocean_state%p_aux%DerivTemperature_vert_center
    
        salinity          => ocean_state%p_prog(nold(1))%ocean_tracers(2)%concentration
        grad_S_vec        => ocean_state%p_aux%PgradSalinity_horz_center
        grad_S_vert_center=> ocean_state%p_aux%DerivSalinity_vert_center
      
      ELSE
        rho_GM                 => ocean_state%p_diag%rho_GM   
        grad_rho_GM_vec        => ocean_state%p_aux%PgradDensity_horz_center
        grad_rho_GM_vert_center=> ocean_state%p_aux%DerivDensity_vert_center

        pot_temp          => ocean_state%p_prog(nold(1))%ocean_tracers(1)%concentration
        grad_T_vec        => ocean_state%p_aux%PgradTemperature_horz_center
        grad_T_vert_center=> ocean_state%p_aux%DerivTemperature_vert_center
        
        salinity          => ocean_state%p_prog(nold(1))%ocean_tracers(2)%concentration
        grad_S_vec        => ocean_state%p_aux%PgradSalinity_horz_center
        grad_S_vert_center=> ocean_state%p_aux%DerivSalinity_vert_center
      
      ENDIF
    !Note that in case of 1-component fluid we use the density for the slope calculation,
    !all relevant imformation is stored in the tracer%temperature structure
    ELSEIF(no_tracer==1)THEN      
      pot_temp          => ocean_state%p_diag%rho	
      grad_T_vec        => ocean_state%p_aux%PgradTemperature_horz_center
      grad_T_vert_center=> ocean_state%p_aux%DerivTemperature_vert_center
    
    ENDIF
    
    !grad_T_vec(:,:,:)%x(1)=0.0_wp
    !grad_T_vec(:,:,:)%x(2)=0.0_wp
    !grad_T_vec(:,:,:)%x(3)=0.0_wp
    !grad_T_vert_center(:,:,:)=0.0_wp

    start_level = 1

    !-------------------------------------------------------------------------------
    !1) calculation of horizontal and vertical gradient for potential temperature and salinity
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_edge_index,end_edge_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      grad_T_horz(:,:,blockNo) = 0.0_wp
      !1a) calculate horizontal gradient of temperature
      CALL grad_fd_norm_oce_3d_onBlock ( &
        & pot_temp, &
        & patch_3D,                           &
        & op_coeff%grad_coeff(:,:,blockNo), &
        & grad_T_horz(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo)

     !1b)calculate horizontal  gradient of salinity
      IF(no_tracer>=2)THEN
        grad_S_horz(:,:,blockNo) = 0.0_wp
        CALL grad_fd_norm_oce_3d_onBlock ( &
          & salinity, &
          & patch_3D,                    &
          & op_coeff%grad_coeff(:,:,blockNo), &
          & grad_S_horz(:,:,blockNo),           &
          & start_edge_index, end_edge_index, blockNo)
     ENDIF
     IF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
        grad_rho_GM_horz(:,:,blockNo) = 0.0_wp
        CALL grad_fd_norm_oce_3d_onBlock ( &
          & rho_GM, &
          & patch_3D,                    &
          & op_coeff%grad_coeff(:,:,blockNo), &
          & grad_rho_GM_horz(:,:,blockNo),           &
          & start_edge_index, end_edge_index, blockNo)
     
     ENDIF
    END DO ! blocks
!ICON_OMP_END_DO


    !---------------------------------------------------------------------
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index) ICON_OMP_DEFAULT_SCHEDULE
 
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

      grad_T_vert(:,:,blockNo) = 0.0_wp ! this is only for the top level
      !1c) calculation of vertical derivative for temperature and salinity
      CALL verticalDeriv_scalar_onHalfLevels_on_block( patch_3d,                &
                                                  & pot_temp(:,:,blockNo),   &
                                                  & grad_T_vert(:,:,blockNo),&
                                                  & start_level+1,             &
                                                  & blockNo,                 &
                                                  & start_cell_index,        &
                                                  & end_cell_index)
      IF((no_tracer>=2))THEN
        grad_S_vert(:,:,blockNo) = 0.0_wp! this is only for the top level
        CALL verticalDeriv_scalar_onHalfLevels_on_block( patch_3d,                &
                                                    & salinity(:,:,blockNo),   &
                                                    & grad_S_vert(:,:,blockNo),&
                                                    & start_level+1,             &
                                                    & blockNo,                 &
                                                    & start_cell_index,        &
                                                    & end_cell_index)
      ENDIF
    IF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
        grad_rho_GM_vert(:,:,blockNo) = 0.0_wp! this is only for the top level
!        CALL verticalDeriv_scalar_onHalfLevels_on_block( patch_3d,                &
!                                                    & rho_GM(:,:,blockNo),   &
!                                                    & grad_rho_GM_vert(:,:,blockNo),&
!                                                    & start_level+1,             &
!                                                    & blockNo,                 &
!                                                    & start_cell_index,        &
!                                                    & end_cell_index)
        grad_rho_GM_vert(:,:,blockNo) = ocean_state%p_diag%grad_rho_PP_vert(:,:,blockNo)
    ENDIF
    END DO ! blocks
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('calc_slopes: grad_T_horz',grad_T_horz,&
      & str_module,idt_src, in_subset=edges_in_domain)
    IF(no_tracer==2)THEN  
    CALL dbg_print('calc_slopes: grad_S_horz',grad_S_horz,&
      & str_module,idt_src, in_subset=edges_in_domain)
    ENDIF
    IF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
    CALL dbg_print('calc_slopes: grad_rho_GM_horz',grad_rho_GM_horz,&
      & str_module,idt_src, in_subset=edges_in_domain)
        
    ENDIF    
!    CALL sync_patch_array(sync_c, patch_2D, grad_T_vert)       
!    IF(no_tracer>=2)   CALL sync_patch_array(sync_c, patch_2D, grad_S_vert)   

   !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('neutral_slopes: grad_T_vert',grad_T_vert,&
      & str_module,idt_src, in_subset=cells_in_domain)
    IF(no_tracer==2)THEN        
    CALL dbg_print('neutral_slopes: grad_S_vert',grad_S_vert,&
      & str_module,idt_src, in_subset=cells_in_domain) 
    ENDIF
    IF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN  
    CALL dbg_print('neutral_slopes: grad_rho_GM_vert',grad_rho_GM_vert,&
      & str_module,idt_src, in_subset=cells_in_domain) 
    
    ENDIF   
  !---------------------------------------------------------------------   
   
    !2) map horizontal and vertial derivative to cell centered vector
    CALL map_edges2cell_3d(patch_3D,  &
        & grad_T_horz,                &
        & op_coeff,                   &
        & grad_T_vec,                 &
        & subset_range=cells_in_domain)

!     CALL sync_patch_array(sync_c, patch_2D, grad_T_vert)
    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
    CALL map_scalar_prismtop2center(patch_3d,&
        & grad_T_vert,                &
        & op_coeff,                 &
        & grad_T_vert_center)
    ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
    CALL map_scalar_prismtop2center_GM(patch_3d,&
        & grad_T_vert,                &
        & op_coeff,                 &
        & grad_T_vert_center)    
    
    ENDIF    
        
       
    IF((no_tracer>=2))THEN        
      
      CALL map_edges2cell_3d(patch_3D,  &      
          & grad_S_horz,                &
          & op_coeff,                   &
          & grad_S_vec,                 &
          & subset_range=cells_in_domain)
    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN          
      CALL map_scalar_prismtop2center(patch_3d,&
          & grad_S_vert,                       &
          & op_coeff,                          &
          & grad_S_vert_center)
    ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN      
      CALL map_scalar_prismtop2center_GM(patch_3d,&
          & grad_S_vert,                       &
          & op_coeff,                          &
          & grad_S_vert_center)  
    ENDIF              
    ENDIF
    IF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN

      CALL map_edges2cell_3d(patch_3D,  &      
          & grad_rho_GM_horz,           &
          & op_coeff,                   &
          & grad_rho_GM_vec,            &
          & subset_range=cells_in_domain)
    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN                    
      CALL map_scalar_prismtop2center(patch_3d,&
          & grad_rho_GM_vert,                  &
          & op_coeff,                          &
          & grad_rho_GM_vert_center)
     ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN  
      CALL map_scalar_prismtop2center_GM(patch_3d,&
          & grad_rho_GM_vert,                  &
          & op_coeff,                          &
          & grad_rho_GM_vert_center)     
     
     ENDIF         
    
    ENDIF    
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------    
    SELECT CASE (no_tracer)

      CASE(2)!Two tracer
        IF(SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
!ICON_OMP_PARALLEL PRIVATE(salinityColumn)                             
    salinityColumn(1:n_zlev) = sal_ref  ! in case of absent salinty tracer!
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,neutral_coeff, &
!ICON_OMP  level) ICON_OMP_DEFAULT_SCHEDULE
          DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
            CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

            DO cell_index = start_cell_index, end_cell_index
              end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
              IF(end_level <= min_dolic) CYCLE
                salinityColumn(1:end_level) = salinity(cell_index,1:end_level,blockNo)

                !4.1) calculate slope coefficients as thermal expansion and saline contraction coefficients
                !
                !Nonlinear EOS, slope coefficients are calculated via the McDougall-method
                IF(EOS_TYPE/=1)THEN
          
                  neutral_coeff = calc_neutralslope_coeff_func_onColumn(         &
                  & pot_temp(cell_index,1:end_level,blockNo), salinityColumn(1:end_level),  &
                  & depth_cellinterface(cell_index,2:end_level+1,blockNo), end_level)
                !Linear EOS: slope coefficients are equal to EOS-coefficients
                ELSEIF(EOS_TYPE==1)THEN          
                  neutral_coeff(:,1) = LinearThermoExpansionCoefficient 
                  neutral_coeff(:,2) = LinearHalineContractionCoefficient
                ENDIF
        
                DO level = start_level+1, end_level-1

           ocean_state%p_aux%slopes(cell_index,level,blockNo)%x &
              & = -neutral_coeff(level,1) * grad_T_vec(cell_index,level,blockNo)%x + &
              &      neutral_coeff(level,2) * grad_S_vec(cell_index,level,blockNo)%x

           ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo)=&
             & DOT_PRODUCT(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,&
                           &ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)
           ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo)=&
             & sqrt(ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo))

           ocean_state%p_aux%slopes_drdz(cell_index,level,blockNo)=&
                 & -(-neutral_coeff(level,1) * grad_T_vert_center(cell_index,level,blockNo)+ &
                 &      neutral_coeff(level,2) * grad_S_vert_center(cell_index,level,blockNo))
              
                  ocean_state%p_aux%slopes(cell_index,level,blockNo)%x &
                    & = -(-neutral_coeff(level,1) * grad_T_vec(cell_index,level,blockNo)%x + &
                    &      neutral_coeff(level,2) * grad_S_vec(cell_index,level,blockNo)%x)  &           
                    &   /(-neutral_coeff(level,1) * grad_T_vert_center(cell_index,level,blockNo)+ &
                    &      neutral_coeff(level,2) * grad_S_vert_center(cell_index,level,blockNo)-dbl_eps)

                  ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)=&
                    & DOT_PRODUCT(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,&
                                 &ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)
!write(123,*)'slope squared',level,ocean_state%p_aux%slopes_squared(cell_index,level,blockNo),&
!&sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))
                                  
                END DO
         !Perform nearest neighbor interpolation at level where slopes are not well-defined
!           ocean_state%p_aux%slopes(cell_index,start_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,start_level+1,blockNo)%x   
!           ocean_state%p_aux%slopes(cell_index,end_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,end_level-1,blockNo)%x        
            END DO ! cell_index = start_cell_index, end_cell_index
          END DO  ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
    
        ELSEIF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level, &
!ICON_OMP  level) ICON_OMP_DEFAULT_SCHEDULE
          DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
            CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

            DO cell_index = start_cell_index, end_cell_index
              end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
              IF(end_level <= min_dolic) CYCLE
              !salinityColumn(1:end_level) = salinity(cell_index,1:end_level,blockNo)
                
              DO level = start_level+1, end_level-1

                ocean_state%p_aux%slopes(cell_index,level,blockNo)%x                      &
                & =   grad_rho_GM_vec(cell_index,level,blockNo)%x  &
                &    /(grad_rho_GM_vert_center(cell_index,level,blockNo)-dbl_eps)

                ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)=&
                & DOT_PRODUCT(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,&
                             &ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)

                ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo)=&
               & DOT_PRODUCT(grad_rho_GM_vec(cell_index,level,blockNo)%x,&
               & grad_rho_GM_vec(cell_index,level,blockNo)%x)

               ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo)=&
               & sqrt(ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo))

               ocean_state%p_aux%slopes_drdz(cell_index,level,blockNo)=&
               &         grad_rho_GM_vert_center(cell_index,level,blockNo)            

              END DO                         
!Perform nearest neighbor interpolation at level where slopes are not well-defined
!           ocean_state%p_aux%slopes(cell_index,start_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,start_level+1,blockNo)%x   
!           ocean_state%p_aux%slopes(cell_index,end_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,end_level-1,blockNo)%x                    
            END DO ! cell_index = start_cell_index, end_cell_index
          END DO  ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
    
        ENDIF
        
    CASE(1)!Case of single tracer
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level, &
!ICON_OMP  level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
          IF(end_level <= min_dolic) CYCLE

          DO level = start_level+1, end_level-1
          
 
            ocean_state%p_aux%slopes(cell_index,level,blockNo)%x                      &
             & = - grad_T_vec(cell_index,level,blockNo)%x  &
             &    /(grad_T_vert_center(cell_index,level,blockNo)-dbl_eps)
            
            ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)=&
            & DOT_PRODUCT(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,&
                         &ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)
          END DO
          !Perform nearest neighbor interpolation at level where slopes are not well-defined
!           ocean_state%p_aux%slopes(cell_index,start_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,start_level+1,blockNo)%x   
!           ocean_state%p_aux%slopes(cell_index,end_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,end_level-1,blockNo)%x    
        END DO ! cell_index = start_cell_index, end_cell_index
      END DO  ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL  
    END SELECT


 !---------DEBUG DIAGNOSTICS-------------------------------------------
!     idt_src=3  ! output print level (1-5, fix)
!     CALL dbg_print('calc_slopes: sizegradT1',grad_T_vec(:,:,:)%x(1), &
!       &  str_module,idt_src, in_subset=cells_in_domain)      
!     CALL dbg_print('calc_slopes: sizegradT2',grad_T_vec(:,:,:)%x(2), &
!       &  str_module,idt_src, in_subset=cells_in_domain)      
!     CALL dbg_print('calc_slopes: sizegradT3',grad_T_vec(:,:,:)%x(3), &
!       &  str_module,idt_src, in_subset=cells_in_domain)      
!      
!      
!    IF(no_tracer>=2)  &
!      & CALL dbg_print('calc:slopes: sizegradS',size_grad_S_horz_vec,&
!      &                str_module,idt_src, in_subset=cells_in_domain)
  !---------------------------------------------------------------------   

  !---------------------------------------------------------------------
  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=1  ! output print level (1-5, fix)
  CALL dbg_print('calc_slopes: squared',(ocean_state%p_aux%slopes_squared(:,:,:)),&
    & str_module,idt_src, in_subset=cells_in_domain)
   DO level=1,n_zlev
     CALL dbg_print('calc_slopes: slobe abs',sqrt(ocean_state%p_aux%slopes_squared(:,level,:)),&
       & str_module,idt_src, in_subset=cells_in_domain)
       
  END DO

   DO level=1,n_zlev       
       CALL dbg_print('calc_slopes: slopes dz',ocean_state%p_aux%slopes_drdz(:,level,:),&
       & str_module,idt_src, in_subset=cells_in_domain)
  END DO


   DO level=1,n_zlev       
       CALL dbg_print('calc_slopes: slopes dx',ocean_state%p_aux%slopes_drdx(:,level,:),&
       & str_module,idt_src, in_subset=cells_in_domain)
  END DO
  
   !---------------------------------------------------------------------

 !  DO level= 1, n_zlev  
 !  write(0,*)'max-min vert deriv',level,maxval(grad_T_vert_center(:,level,:)),&
 !  & minval(grad_T_vert_center(:,level,:))!,maxval(grad_S_vert_center(:,level,:)),&
 !  !& minval(grad_S_vert_center(:,level,:))
 !  END DO
 
 
  !--------------------------------------------------------------------- 
  IF(BOLUS_VELOCITY_DIAGNOSTIC)THEN
    CALL calc_bolus_velocity(patch_3d, ocean_state, param, op_coeff)
  ENDIF
  !---------------------------------------------------------------------  
  
  
  END SUBROUTINE calc_neutral_slopes
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates horizontal and vertical derivatives of tracer and the corresponding
  !!    !reconstructions on cell centers.
  !!
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2016).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_tracer_derivatives(patch_3d, tracer, ocean_state, op_coeff, tracer_index)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    REAL(wp)                                         :: tracer(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_operator_coeff),            INTENT(in)    :: op_coeff
    INTEGER, INTENT(in)                              :: tracer_index
    
    !Local variables
    TYPE(t_cartesian_coordinates),POINTER :: grad_horz_center_vec(:,:,:)    
    REAL(wp),POINTER :: grad_vert_center(:,:,:)
    REAL(wp),POINTER :: grad_vert(:,:,:)    
    REAL(wp),POINTER :: grad_horz(:,:,:)
               
    INTEGER :: level, blockNo, je, start_level, end_level
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    all_cells       => patch_2D%cells%all
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain
    !-------------------------------------------------------------------------
    IF(tracer_index==1)THEN
    
      grad_horz            => ocean_state%p_aux%temperature_grad_horz
      grad_vert            => ocean_state%p_aux%temperature_deriv_vert        

      grad_horz_center_vec => ocean_state%p_aux%PgradTemperature_horz_center
      grad_vert_center     => ocean_state%p_aux%DerivTemperature_vert_center
          
    ELSEIF(tracer_index==2)THEN      
    
      grad_horz             => ocean_state%p_aux%salinity_grad_horz
      grad_vert             => ocean_state%p_aux%salinity_deriv_vert  
          
      grad_horz_center_vec  => ocean_state%p_aux%PgradSalinity_horz_center
      grad_vert_center      => ocean_state%p_aux%DerivSalinity_vert_center
  
    ELSEIF(tracer_index>2)THEN           

      grad_horz            => ocean_state%p_aux%tracer_grad_horz
      grad_vert            => ocean_state%p_aux%tracer_deriv_vert
      
      grad_horz_center_vec => ocean_state%p_aux%PgradTracer_horz_center
      grad_vert_center     => ocean_state%p_aux%DerivTracer_vert_center
    
    ENDIF


    !-------------------------------------------------------------------------------
    !1) calculation of horizontal and vertical gradient for potential temperature and salinity
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_edge_index,end_edge_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      grad_horz(:,:,blockNo) = 0.0_wp
      !1a) calculate horizontal gradient
      CALL grad_fd_norm_oce_3d_onBlock ( &
        & tracer, &
        & patch_3D,                           &
        & op_coeff%grad_coeff(:,:,blockNo), &
        & grad_horz(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo)

    END DO ! blocks
!ICON_OMP_END_DO


    !---------------------------------------------------------------------
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index) ICON_OMP_DEFAULT_SCHEDULE
 
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

      grad_vert(:,:,blockNo) = 0.0_wp ! this is only for the top level
      !1c) calculation of vertical derivative
      CALL verticalDeriv_scalar_onHalfLevels_on_block( patch_3d,                &
                                                  & tracer(:,:,blockNo),   &
                                                  & grad_vert(:,:,blockNo),&
                                                  & start_level+1,             &
                                                  & blockNo,                 &
                                                  & start_cell_index,        &
                                                  & end_cell_index)
    END DO ! blocks
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('calc_derivatives: grad_horz',grad_horz,&
      & str_module,idt_src, in_subset=edges_in_domain)

    CALL dbg_print('calc_derivatives: grad_vert',grad_vert,&
      & str_module,idt_src, in_subset=cells_in_domain)
  !---------------------------------------------------------------------   
   
    !2) map horizontal and vertial derivative to cell centered vector
    CALL map_edges2cell_3d(patch_3D,  &
        & grad_horz,                  &
        & op_coeff,                   &
        & grad_horz_center_vec,       &
        & subset_range=cells_in_domain)


    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
    CALL map_scalar_prismtop2center(patch_3d,&
        & grad_vert,                &
        & op_coeff,                 &
        & grad_vert_center)
    ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
    CALL map_scalar_prismtop2center_GM(patch_3d,&
        & grad_vert,                &
        & op_coeff,                 &
        & grad_vert_center)    
    
    ENDIF    
    !------------------------------------------------------------------------------


  END SUBROUTINE calc_tracer_derivatives
  !-------------------------------------------------------------------------



 !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the bolus velocity from the isopycnal slopes.
  !! !  This is currently used for diagnostic purposes, the GM-Redi paramerization
  !! !  uses the skew-flux approach.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2016).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_bolus_velocity(patch_3d, ocean_state, param, op_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(inout) :: op_coeff
    
    !Local variables
    REAL(wp) :: div_bolus_vn(nproma, n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: cell_max_slope, inv_cell_characteristic_length,slope_abs
    
    TYPE(t_cartesian_coordinates) :: grad_slope_vec(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_cartesian_coordinates) :: kappa_times_slope(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: slope_deriv_atcenter(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
           
    INTEGER :: level, blockNo, start_level, end_level
    INTEGER :: start_cell_index, end_cell_index, cell_index
    TYPE(t_subset_range), POINTER :: cells_in_domain, all_cells!, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    all_cells       => patch_2D%cells%all
    cells_in_domain => patch_2D%cells%in_domain
    !edges_in_domain => patch_2D%edges%in_domain
    !-------------------------------------------------------------------------

    start_level = 1
    
    
    !Step 1: Calculate kappa times slopes
    !-------------------------------------------------------------------------------
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level, &
!ICON_OMP  level) ICON_OMP_DEFAULT_SCHEDULE
          DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
            CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

            DO cell_index = start_cell_index, end_cell_index
            
              inv_cell_characteristic_length = 1.0_wp / SQRT(patch_2D%cells%area(cell_index,blockNo))
              end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
              IF(end_level <= min_dolic) CYCLE
                              
              DO level = start_level+1, end_level-1

                cell_max_slope      = S_max  &
                & * patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo) &
                & * inv_cell_characteristic_length
              
                slope_abs = sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))

                IF(slope_abs <= cell_max_slope)THEN
                  kappa_times_slope(cell_index,level,blockNo)%x=        &
                  &-param%k_tracer_GM_kappa(cell_index,level,blockNo)   &
                  &+ocean_state%p_aux%slopes(cell_index,level,blockNo)%x
                ELSE
                  kappa_times_slope(cell_index,level,blockNo)%x=0.0_wp                      
                ENDIF
              END DO                         
            END DO ! cell_index = start_cell_index, end_cell_index
          END DO  ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
    !-------------------------------------------------------------------------------

    !Step 2: Calculate vertical derivative of kappa times slopes and map it back to prism center(midlevel)
    !-------------------------------------------------------------------------------           
!ICON_OMP_PARALLEL   
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index) ICON_OMP_DEFAULT_SCHEDULE
 
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
       
      CALL verticalDeriv_vec_midlevel_on_block( patch_3d,                      &
                                              & kappa_times_slope(:,:,blockNo),&
                                              & grad_slope_vec(:,:,blockNo),   &
                                              & start_level+1,                 &
                                              & blockNo,                       &
                                              & start_cell_index,              &
                                              & end_cell_index)
       IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN                                       
       CALL map_vec_prismtop2center_on_block( patch_3d,                   &
                                            & grad_slope_vec(:,:,blockNo),&
                                            & op_coeff,                   &
                                            & slope_deriv_atcenter,       &
                                            & blockNo, start_cell_index, end_cell_index) 
       ELSE
       CALL map_vec_prismtop2center_on_block_GM( patch_3d,                   &
                                            & grad_slope_vec(:,:,blockNo),&
                                            & op_coeff,                   &
                                            & slope_deriv_atcenter,       &
                                            & blockNo, start_cell_index, end_cell_index) 
              
       
       ENDIF                                                                           
                                              
    END DO ! blocks
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL        
    !------------------------------------------------------------------------------

    !Step 3: Map result back to edges to obtain the horizontal bolus velocity 
    !------------------------------------------------------------------------------
    CALL  map_cell2edges_3d( patch_3d,                   &
                           & slope_deriv_atcenter,       &
                           & ocean_state%p_diag%vn_bolus,&
                           & op_coeff,                   &
                           & start_level+1)
    !------------------------------------------------------------------------------


    !Step 4: Calculate vertical bolus as horizontal divergence of horizontal bolus velocity
    !------------------------------------------------------------------------------
!ICON_OMP_PARALLEL    
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index) ICON_OMP_DEFAULT_SCHEDULE
          DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
            CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
        
            CALL div_oce_3D_onTriangles_onBlock( &
              & ocean_state%p_diag%vn_bolus,     &
              & patch_3D,op_coeff%div_coeff,     &
              & div_bolus_vn(:,:,blockNo),       &
              & blockNo,start_cell_index, end_cell_index,start_level+1,n_zlev)!,      &
              !& start_level=1, end_level=n_zlev)
        
            DO cell_index = start_cell_index, end_cell_index          
              !use bottom boundary condition for vertical velocity at bottom of prism
              DO level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo), 1, -1
                ocean_state%p_diag%w_bolus(cell_index,level,blockNo) &
                  &= ocean_state%p_diag%w_bolus(cell_index,level+1,blockNo)&
                  & - div_bolus_vn(cell_index,level,blockNo)
              END DO
            END DO
        
          END DO ! blockNo
!ICON_OMP_END_PARALLEL_DO 
!ICON_OMP_END_PARALLEL
    !------------------------------------------------------------------------------
!    CALL map_edges2cell_3d(patch_3d, ocean_state%p_diag%vn_bolus, op_coeff,ocean_state%p_diag%p_vn)
!    CALL Get3DVectorTo2DLocal_array3D(vector=ocean_state%p_diag%p_vn, &
!      & position_local=patch_2d%cells%center, &
!      & levels=patch_3d%p_patch_1d(1)%dolic_c, &
!      & subset=all_cells,                      &
!      & geometry_info=patch_2d%geometry_info, &
!      & x=ocean_state%p_diag%u, y=ocean_state%p_diag%v)

  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=1  ! output print level (1-5, fix)
  CALL dbg_print('calc_bolus:',(ocean_state%p_diag%w_bolus(:,:,:)),&
    & str_module,idt_src, in_subset=cells_in_domain)
   DO level=1,n_zlev
     CALL dbg_print('calc_bolus:',ocean_state%p_diag%w_bolus(:,level,:),&
       & str_module,idt_src, in_subset=cells_in_domain)
  END DO
!stop
  !---------------------------------------------------------------------

  END SUBROUTINE calc_bolus_velocity
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the tapering functions used in Gent-McWilliams-Redi parametrization.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
!<Optimize:inUse>
  SUBROUTINE calc_tapering_function(patch_3d, param, ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    !TYPE(t_operator_coeff),            INTENT(in)    :: op_coeff
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level, blockNo
    INTEGER :: start_edge_index, end_edge_index, je     
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: lambda
    REAL(wp) :: depth_scale, depth
    REAL(wp) :: Coriolis_abs
    REAL(wp) :: inv_S_d, slope_abs, inv_cell_characteristic_length, cell_max_slope, cell_critical_slope
    REAL(wp) :: cell_characteristic_length
    
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    start_level=1
    inv_S_d = 1.0_wp / S_d
    
    
    !-------------------------------------------------------------------------------
!ICON_OMP_PARALLEL   
    IF (GMRedi_usesRelativeMaxSlopes) THEN
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, slope_abs, inv_cell_characteristic_length, &
!ICON_OMP cell_max_slope, cell_critical_slope) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
          inv_cell_characteristic_length = 1.0_wp / SQRT(patch_2D%cells%area(cell_index,blockNo))
          cell_characteristic_length     = SQRT(patch_2D%cells%area(cell_index,blockNo))
          
          IF(end_level >= min_dolic) THEN

            DO level = start_level, end_level

              cell_critical_slope      = S_critical  &
                & * patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo) &
                & * inv_cell_characteristic_length
              !cell_max_slope      = S_max*SQRT(patch_2D%cells%area(cell_index,blockNo))&
              !& / patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo)
!              cell_critical_slope = S_critical &
!                & * patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo) &
!                & * inv_cell_characteristic_length

              cell_max_slope = S_max &
                & * patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo) &
                & * cell_characteristic_length&
                &/(2.0_wp*param%k_tracer_isoneutral(cell_index,level,blockNo)*dtime)
                
              slope_abs    = sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))           

              IF(slope_abs <= cell_max_slope)THEN

                ocean_state%p_aux%taper_function_1(cell_index,level,blockNo) &
                  &= 0.5_wp*(1.0_wp + tanh((cell_critical_slope - slope_abs)*inv_S_d))

              ELSE
                ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)=0.0_wp
              ENDIF

            END DO
          ENDIF
        END DO
      END DO

!ICON_OMP_END_DO
    ELSE! not GMRedi_usesRelativeMaxSlopes
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, slope_abs, inv_cell_characteristic_length, &
!ICON_OMP cell_max_slope, cell_critical_slope) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        cell_max_slope      = S_max
        cell_critical_slope = S_critical
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
          inv_cell_characteristic_length = 1.0_wp / SQRT(patch_2D%cells%area(cell_index,blockNo))

          IF(end_level >= min_dolic) THEN

            DO level = start_level, end_level
            
             slope_abs    = sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))           
              
              !IF(slope_abs <= cell_max_slope)THEN
                ocean_state%p_aux%taper_function_1(cell_index,level,blockNo) &
                  &= 0.5_wp*(1.0_wp + tanh((cell_max_slope - slope_abs)*inv_S_d))
              !ELSE
              !  ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)=0.0_wp
              !ENDIF

            END DO
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO
    ENDIF !GMRedi_usesRelativeMaxSlopes

    !tapering schemes other than Danabasoglu-McWilliams require a second
    !tapering function
    IF(tapering_scheme/=tapering_DanaMcWilliams)THEN

!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, &
!ICON_OMP Coriolis_abs,lambda,depth_scale,depth) ICON_OMP_DEFAULT_SCHEDULE

      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          IF(end_level >= min_dolic) THEN

            DO level = start_level, end_level-1
              Coriolis_abs=abs(patch_2d%cells%f_c(cell_index,blockNo))
              lambda      = min(max(15000._WP, c_speed/Coriolis_abs),100000.0_WP)
              depth_scale = lambda*sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))
              depth       = patch_3d%p_patch_1d(1)%depth_CellMiddle(cell_index,level,blockNo)

              IF(depth<=depth_scale)THEN
                ocean_state%p_aux%taper_function_2(cell_index,level,blockNo) &
                  & = 0.5_wp*(1.0_wp+sin(pi*depth/depth_scale-pi/2.0_wp))
              ELSE
                ocean_state%p_aux%taper_function_2(cell_index,level,blockNo) =1.0_wp
              ENDIF

            END DO
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO
!       CALL sync_patch_array(sync_c, patch_2D,ocean_state%p_aux%taper_function_2)
    ENDIF
!ICON_OMP_END_PARALLEL

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
     CALL dbg_print('calc_tapering_fct1:', ocean_state%p_aux%taper_function_1 ,&
       & this_mod_name, 3, patch_2D%cells%in_domain)
     IF(tapering_scheme/=tapering_DanaMcWilliams) THEN
       CALL dbg_print('calc_tapering_fct2:', ocean_state%p_aux%taper_function_2 ,&
       & this_mod_name, 3, patch_2D%cells%in_domain)
     ENDIF
   
  END SUBROUTINE calc_tapering_function
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the entries of he Gent-McWilliams mixing tensor, 
  !! !             this includes the tapering for which different options are available
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
!<Optimize:inUse>
  SUBROUTINE calc_entries_mixing_tensor(patch_3d, ocean_state, param,taper_diagonal_horz,taper_diagonal_vert_expl,&
    & taper_diagonal_vert_impl, taper_off_diagonal_horz, taper_off_diagonal_vert,tracer_index )
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    REAL(wp), INTENT(inout)                          :: taper_diagonal_horz(:,:,:)
    REAL(wp), INTENT(inout)                          :: taper_diagonal_vert_expl(:,:,:)        
    REAL(wp), INTENT(inout)                          :: taper_diagonal_vert_impl(:,:,:)    
    TYPE(t_cartesian_coordinates), INTENT(inout)     :: taper_off_diagonal_horz(:,:,:)
    TYPE(t_cartesian_coordinates), INTENT(inout)     :: taper_off_diagonal_vert(:,:,:) 
    INTEGER                                          :: tracer_index   
    
    
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
    start_level=1
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !
    !The dianeutral diffusivity is the number determined by the PP-scheme
    K_I           => param%k_tracer_isoneutral
    K_D           => param%k_tracer_dianeutral !param%a_tracer_v(:,:,:,tracer_index) !  
    kappa         => param%k_tracer_GM_kappa
    
    IF(SWITCH_OFF_TAPERING)THEN
     ocean_state%p_aux%taper_function_1(:,:,:)=1.0_wp
     ocean_state%p_aux%taper_function_2(:,:,:)=1.0_wp
    ENDIF
    taper_diagonal_horz(:,:,:)     =0.0_wp
    taper_diagonal_vert_expl(:,:,:)=0.0_wp        
    taper_diagonal_vert_impl(:,:,:)=0.0_wp    
    taper_off_diagonal_horz(:,:,:)%x(1)=0.0_wp
    taper_off_diagonal_horz(:,:,:)%x(2)=0.0_wp
    taper_off_diagonal_horz(:,:,:)%x(3)=0.0_wp    
    taper_off_diagonal_vert(:,:,:)%x(1)=0.0_wp
    taper_off_diagonal_vert(:,:,:)%x(2)=0.0_wp
    taper_off_diagonal_vert(:,:,:)%x(3)=0.0_wp
    !-------------------------------------------------------------------------------
    
    SELECT CASE(tapering_scheme)

    CASE(tapering_DanaMcWilliams)
      IF(SWITCH_ON_TAPERING_HORIZONTAL_DIFFUSION)THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE

      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
          
            DO level = start_level, end_level
            
              !Following recommendations in Griffies Ocean Climate Model book (sect. 15.3.4.2)
              !Danabasoglou-McWilliams tapering is for horizontal flux only applied to
              !off-diagonal terms
              
              !coefficients for horizontal fluxes
              taper_diagonal_horz(cell_index,level,blockNo)  &
              & = ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)*&
              &   K_I(cell_index,level,blockNo)

              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
              & = (K_I(cell_index,level,blockNo)-kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)

              !coefficients for vertical fluxes
              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
              &= (K_I(cell_index,level,blockNo)+kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&                
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x

              
              taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
              &=K_I(cell_index,level,blockNo)&
              &*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)
            END DO
        END DO
      END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL      
      
      ELSEIF(.NOT.SWITCH_ON_TAPERING_HORIZONTAL_DIFFUSION)THEN

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE

      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

            DO level = start_level, end_level
            
              !Following recommendations in Griffies Ocean Climate Model book (sect. 15.3.4.2)
              !Danabasoglou-McWilliams tapering is for horizontal flux only applied to
              !off-diagonal terms
              
              !coefficients for horizontal fluxes
              taper_diagonal_horz(cell_index,level,blockNo)  &
              & = K_I(cell_index,level,blockNo)

              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
              & = (K_I(cell_index,level,blockNo)-kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)

              !coefficients for vertical fluxes
              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
              &= (K_I(cell_index,level,blockNo)+kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&                
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x

              
              taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
              &=K_I(cell_index,level,blockNo)&
              &*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)

            END DO

        END DO
      END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL      
      ENDIF
      IF(switch_off_diagonal_vert_expl)THEN
!ICON_OMP_PARALLEL      
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

!             IF(end_level >= min_dolic) THEN
          
              DO level = start_level, end_level
                taper_diagonal_vert_expl(cell_index,level,blockNo)=0.0_wp
              END DO
!             ENDIF
          END DO
        END DO      
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
      ELSEIF(.NOT.switch_off_diagonal_vert_expl)THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE
       DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

!             IF(end_level >= min_dolic) THEN
          
              DO level = start_level, end_level
      
                taper_diagonal_vert_expl(cell_index,level,blockNo)  &
                &=K_D(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)
               END DO
!               ENDIF
            END DO
          END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
      ENDIF            

    CASE(tapering_Large)
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, geometric_scale) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          IF(end_level >= min_dolic) THEN
             !geometric_scale=1.0_wp!sqrt(patch_2D%cells%area(cell_index,blockNo))*geometric_scale_factor_GMR
            DO level = start_level, end_level
            
              !coefficients for horizontal fluxes           
              taper_diagonal_horz(cell_index,level,blockNo)  &
              &=K_I(cell_index,level,blockNo)

              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
              & = ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
              & * ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)&
              & *(K_I(cell_index,level,blockNo)-kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x

              !coefficients for vertical fluxes
              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
              &= (K_I(cell_index,level,blockNo)+kappa(cell_index,level,blockNo))&
              &+ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)&
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x

              
              taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
              &=K_I(cell_index,level,blockNo)&
              &*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)
 
            END DO
          ENDIF
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO
! !ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE     
      IF(switch_off_diagonal_vert_expl)THEN
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

            IF(end_level >= min_dolic) THEN
          
              DO level = start_level, end_level
                taper_diagonal_vert_expl(cell_index,level,blockNo)=0.0_wp
              END DO
            ENDIF
          END DO
        END DO      
      ELSEIF(.NOT.switch_off_diagonal_vert_expl)THEN
       DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

            IF(end_level >= min_dolic) THEN
          
              DO level = start_level, end_level
      
                taper_diagonal_vert_expl(cell_index,level,blockNo)  &
                &=K_D(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)
               END DO
              ENDIF
            END DO
          END DO
        
      ENDIF      
! !ICON_OMP_END_PARALLEL_DO
               

    CASE(tapering_Griffies)
    CALL finish(TRIM('mo_ocean_GM_Redi'), 'The tapering option = tapering_Griffies is not supported yet')    
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, geometric_scale) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          IF(end_level >= min_dolic) THEN
            geometric_scale=1.0_wp!sqrt(patch_2D%cells%area(cell_index,blockNo))*geometric_scale_factor_GMR
            DO level = start_level, end_level
            
              taper_diagonal_horz(cell_index,level,blockNo)  &
                &=K_I(cell_index,level,blockNo)

              taper_diagonal_vert_expl(cell_index,level,blockNo)  &
                &=K_D(cell_index,level,blockNo)              

               taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
                 &=(K_I(cell_index,level,blockNo)*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
                 &+K_D(cell_index,level,blockNo))*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
                 &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)              
              
              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
                = (K_I(cell_index,level,blockNo)-kappa(cell_index,level,blockNo))&
                &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)

              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
                = (K_I(cell_index,level,blockNo)+kappa(cell_index,level,blockNo))&
                &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)
            
            END DO
          ENDIF
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    END SELECT

Do level=1,n_zlev    
    CALL dbg_print('calc_mixing_tensor: horz diag', taper_diagonal_horz(:,level,:),&
    & this_mod_name, 3, patch_2D%cells%in_domain)
END DO   
Do level=1,n_zlev
    CALL dbg_print('calc_mixing_tensor: vert diag expl', taper_diagonal_vert_expl(:,level,:),&
    & this_mod_name, 3, patch_2D%cells%in_domain)
END DO    
Do level=1,n_zlev
    CALL dbg_print('calc_mixing_tensor: vert diag impl', taper_diagonal_vert_impl(:,level,:),&
    & this_mod_name, 3, patch_2D%cells%in_domain)
END DO
  END SUBROUTINE calc_entries_mixing_tensor
  !-------------------------------------------------------------------------
  
  
   !-------------------------------------------------------------------------
  !
  !>
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
!<Optimize:inUse>
  SUBROUTINE selective_GMRedi(patch_3d, ocean_state, param,taper_diagonal_horz,taper_diagonal_vert_expl,&
    & taper_diagonal_vert_impl, taper_off_diagonal_horz, taper_off_diagonal_vert,tracer_index )
    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    REAL(wp), INTENT(inout)                          :: taper_diagonal_horz(:,:,:)
    REAL(wp), INTENT(inout)                          :: taper_diagonal_vert_expl(:,:,:)        
    REAL(wp), INTENT(inout)                          :: taper_diagonal_vert_impl(:,:,:)    
    TYPE(t_cartesian_coordinates), INTENT(inout)     :: taper_off_diagonal_horz(:,:,:)
    TYPE(t_cartesian_coordinates), INTENT(inout)     :: taper_off_diagonal_vert(:,:,:) 
    INTEGER                                          :: tracer_index   
    
    
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
    start_level=1
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !
    !The dianeutral diffusivity is the number determined by the PP-scheme
    K_I           => param%k_tracer_isoneutral
    K_D           => param%k_tracer_dianeutral !param%a_tracer_v(:,:,:,tracer_index) !  
    kappa         => param%k_tracer_GM_kappa
    !-------------------------------------------------------------------------------
    
    SELECT CASE(tapering_scheme)

    CASE(tapering_DanaMcWilliams)
    
    
    IF(TEST_MODE_GM_ONLY)THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

!           IF(end_level >= min_dolic) THEN
          
            DO level = start_level, end_level
            
              !Following recommendations in Griffies Ocean Climate Model book (sect. 15.3.4.2)
              !Danabasoglou-McWilliams tapering is for horizontal flux only applied to
              !off-diagonal terms
              
              !coefficients for horizontal fluxes: horizontal diffusion is retained here !
              taper_diagonal_horz(cell_index,level,blockNo)  &
              & =&! ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)*&
              &   K_I(cell_index,level,blockNo)

              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
              & =(-kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)

              !coefficients for vertical fluxes
              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
              &= (kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&                
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x

              !implicit part of vertical diffusion due to GM is here set to zero.
              !The implicit part due to PP-scheme is retained.
              taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
              &=0.0_wp!
              
              taper_diagonal_vert_expl(cell_index,level,blockNo)=0.0_wp
            END DO
!           ENDIF
        END DO
      END DO
!ICON_OMP_END_DO

!ICON_OMP_END_PARALLEL
     ELSEIF(TEST_MODE_REDI_ONLY)THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

!           IF(end_level >= min_dolic) THEN
          
            DO level = start_level, end_level
            
              !Following recommendations in Griffies Ocean Climate Model book (sect. 15.3.4.2)
              !Danabasoglou-McWilliams tapering is for horizontal flux only applied to
              !off-diagonal terms
              
              !coefficients for horizontal fluxes: horizontal diffusion is untapered.
              taper_diagonal_horz(cell_index,level,blockNo)  &
              & =&! ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)*&
              &   K_I(cell_index,level,blockNo)

              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
              & = (K_I(cell_index,level,blockNo))&
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)

              !coefficients for vertical fluxes
              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
              &= (K_I(cell_index,level,blockNo))&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&                
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x

              
              taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
              &=K_I(cell_index,level,blockNo)&
              &*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)
              
            END DO
!           ENDIF
        END DO
      END DO
!ICON_OMP_END_DO

      IF(switch_off_diagonal_vert_expl)THEN
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

!             IF(end_level >= min_dolic) THEN
          
              DO level = start_level, end_level
                taper_diagonal_vert_expl(cell_index,level,blockNo)=0.0_wp
              END DO
!             ENDIF
          END DO
        END DO      
!ICON_OMP_END_DO
      ELSEIF(.NOT.switch_off_diagonal_vert_expl)THEN
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE
       DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

!             IF(end_level >= min_dolic) THEN
          
              DO level = start_level, end_level
      
                taper_diagonal_vert_expl(cell_index,level,blockNo)  &
                &=K_D(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)
               END DO
!               ENDIF
            END DO
          END DO
!ICON_OMP_END_DO
        
      ENDIF      
!ICON_OMP_END_PARALLEL
     
     
        
     ENDIF
    CASE(tapering_Large)
       CALL finish(TRIM('mo_ocean_GM_Redi'), 'The testing mode for LARGE-tapering option is not supported yet')        

    CASE(tapering_Griffies)
    CALL finish(TRIM('mo_ocean_GM_Redi'), 'The tapering option = tapering_Griffies is not supported yet')    
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, geometric_scale) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          IF(end_level >= min_dolic) THEN
            geometric_scale=1.0_wp!sqrt(patch_2D%cells%area(cell_index,blockNo))*geometric_scale_factor_GMR
            DO level = start_level, end_level
            
              taper_diagonal_horz(cell_index,level,blockNo)  &
                &=K_I(cell_index,level,blockNo)

              taper_diagonal_vert_expl(cell_index,level,blockNo)  &
                &=K_D(cell_index,level,blockNo)              

               taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
                 &=(K_I(cell_index,level,blockNo)*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
                 &+K_D(cell_index,level,blockNo))*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
                 &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)              
              
              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
                = (K_I(cell_index,level,blockNo)-kappa(cell_index,level,blockNo))&
                &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)

              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
                = (K_I(cell_index,level,blockNo)+kappa(cell_index,level,blockNo))&
                &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)
            
! write(123,*)'data',taper_diagonal_vert_impl(cell_index,level,blockNo),&
! &ocean_state%p_aux%slopes_squared(cell_index,level,blockNo),&
! &ocean_state%p_aux%taper_function_1(cell_index,level,blockNo),&
! &(K_I(cell_index,level,blockNo)*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
! &+K_D(cell_index,level,blockNo)),&
! Dot_Product(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)
            
            END DO
          ENDIF
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    END SELECT
! Do level=start_level,4!end_level
! write(*,*)'max/min',level,&
! & maxval(taper_diagonal_horz(:,level,:)),minval(taper_diagonal_horz(:,level,:)),&
! & maxval(taper_diagonal_vert_expl(:,level,:)),minval(taper_diagonal_vert_expl(:,level,:)),&
! & maxval(taper_diagonal_vert_impl(:,level,:)),minval(taper_diagonal_vert_impl(:,level,:)) 
! END DO
    
!     CALL sync_patch_array(sync_c, patch_2D,K_I)
!     CALL sync_patch_array(sync_c, patch_2D,K_D)
!     CALL sync_patch_array(sync_c, patch_2D,kappa)
Do level=1,n_zlev
    CALL dbg_print('selective_GMRedi: vert diag expl', taper_diagonal_vert_expl(:,level,:),&
    & this_mod_name, 1, patch_2D%cells%in_domain)
END DO    
Do level=1,n_zlev
    CALL dbg_print('selective_GMRedi: vert off-diag impl', taper_diagonal_vert_impl(:,level,:),&
    & this_mod_name, 1, patch_2D%cells%in_domain)
END DO
Do level=1,n_zlev    
    CALL dbg_print('selective_GMRedi: horz diag', taper_diagonal_horz(:,level,:),&
    & this_mod_name, 1, patch_2D%cells%in_domain)
END DO   
  END SUBROUTINE selective_GMRedi
  !-------------------------------------------------------------------------
  
   
  

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
  !<Optimize:inUse>
  FUNCTION calc_neutralslope_coeff_func_onColumn(t,s,p,levels) result(coeff)
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
    INTEGER, INTENT(in)   :: levels
    REAL(wp), INTENT(in)  :: t(:)        !  potential temperature (in ITS-90) [C]
    REAL(wp), INTENT(in)  :: s(:)        !  salinity (in PSS-78) [psu]
    REAL(wp), INTENT(in)  :: p(:)       !  pressure (in dezi-bar) [db]
    REAL(wp)              :: coeff(1:n_zlev,2) !  thermal expansion [1/C] and saline contraction [1/psu] coefficients

    ! local variables, following the naming of the FESOM implementation
    REAL(wp):: aob(1:levels), t1(1:levels), t2(1:levels), t3(1:levels), t4(1:levels), &
      & s35(1:levels), s35sq(1:levels), s1(1:levels), s2(1:levels), s3(1:levels), p1(1:levels), p2(1:levels), p3(1:levels)

    INTEGER :: level

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

    ! t1 = t
    s1(1:levels) = s(1:levels)
    p1(1:levels) = p(1:levels)

   ! correction factor for conversion of 1990 to 1968 temperature standard (IPTS-68 to IPTS-90)
   ! the correction is less than 0.01 K in ocean water temperature range
   !  - T68 = 1.00024*T90
   !  - above mentioned CHECK VALUES of the paper are better met by this correction
    t1(1:levels) = t(1:levels) * 1.00024_wp

    t2(1:levels)    = t1(1:levels) * t1(1:levels)
    t3(1:levels)    = t2(1:levels) * t1(1:levels)
    t4(1:levels)    = t3(1:levels) * t1(1:levels)
    p2(1:levels)    = p1(1:levels )* p1(1:levels)
    p3(1:levels)    = p2(1:levels )* p1(1:levels)
    s35(1:levels)   = s(1:levels) - 35.0_wp
    s35sq(1:levels) = s35(1:levels) * s35(1:levels)

    DO level=1,levels
      ! calculate beta, saline contraction
      coeff(level,2) = bet_t0 - bet_t1*t1(level)                             &
        &         + bet_t2*t2(level) - bet_t3*t3(level)                      &
        &         + s35(level)*(-bet_st0    + bet_st1*t1(level)              &
        &         +       bet_sp1*p1(level) - bet_sp2*p2(level))             &
        &         + s35sq(level)*bet_s2                                      &
        &         + p1(level)*(-bet_p1t0 + bet_p1t1*t1(level) - bet_p1t2*t2(level)) &
        &         + p2(level)*( bet_p2t0 - bet_p2t1*t1(level))               &
        &         + p3(level)*bet_p3

      ! calculate alpha/beta
      aob(level) = aob_t0 + aob_t1*t1 (level)                                &
        &         - aob_t2*t2(level) + aob_t3*t3(level)                      &
        &         - aob_t4*t4(level)                                         &
        &         + s35(level)*(+aob_st0    - aob_st1*t1(level)              &
        &                -aob_sp1*p1(level) - aob_sp2*p2(level))             &
        &         - s35sq(level)*aob_s2                                      &
        &         + p1(level)*(+aob_p1t0 - aob_p1t1*t1(level) + aob_p1t2*t2(level)) &
        &         + p2(level)*t2(level)*aob_p2t2                             &
        &         - p3(level)*aob_p3

      ! calculate alpha, thermal expansion
      coeff(level,1) = aob(level)* coeff(level, 2)
    ENDDO

  END FUNCTION calc_neutralslope_coeff_func_onColumn

  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi and the eddy flux of GM.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2014).
  !!
  !<Optimize:inUse>
  SUBROUTINE diagnose_Redi_flux_balance(patch_3d, ocean_state, param, op_coeff,&
    &tracer_index)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)  :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params),      INTENT(inout)    :: param
    TYPE(t_operator_coeff), INTENT(in)       :: op_coeff
    INTEGER, INTENT(IN)                      :: tracer_index 
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level,blockNo
    INTEGER :: start_edge_index, end_edge_index, je     
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: flux_vert_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: flux_vec_horz_center(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates),POINTER :: tracer_gradient_horz_vec_center(:,:,:), slopes(:,:,:)
    REAL(wp), POINTER :: tracer_gradient_vert_center(:,:,:)

    TYPE(t_cartesian_coordinates)            :: taper_off_diagonal_vert(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)    
    TYPE(t_cartesian_coordinates)            :: taper_off_diagonal_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                                 :: taper_diagonal_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                                 :: taper_diagonal_vert_expl(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)    
    REAL(wp)                                 :: taper_diagonal_vert_impl(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)            
    REAL(wp) :: mapped_vertical_diagonal_impl(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)     
    TYPE(t_cartesian_coordinates) :: flux_sum_horz(nproma,n_zlev, patch_3d%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                      :: flux_sum_vert(nproma,n_zlev, patch_3d%p_patch_2D(1)%alloc_cell_blocks)                      
    REAL(wp) :: flux_vert_temp(nproma,n_zlev, patch_3d%p_patch_2D(1)%alloc_cell_blocks)           
    REAL(wp) :: flux_vert_sal(nproma,n_zlev, patch_3d%p_patch_2D(1)%alloc_cell_blocks)           
    REAL(wp):: neutral_coeff(1:n_zlev, 2), salinityColumn(1:n_zlev)
    REAL(wp), POINTER :: pot_temp(:,:,:), salinity(:,:,:)
    REAL(wp), POINTER :: depth_cellinterface(:,:,:)
    
    REAL(wp)        :: GMredi_flux_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp)        :: GMredi_flux_vert(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    all_cells       => patch_2D%cells%all
    edges_in_domain => patch_2D%edges%in_domain 
    slopes          => ocean_state%p_aux%slopes 
    depth_cellinterface => patch_3D%p_patch_1d(1)%depth_cellinterface
    
    
    GMredi_flux_horz(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)=0.0_wp
    GMredi_flux_vert(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)=0.0_wp
    
   IF(no_tracer>=2)THEN
      pot_temp          => ocean_state%p_prog(nold(1))%ocean_tracers(1)%concentration
      salinity          => ocean_state%p_prog(nold(1))%ocean_tracers(2)%concentration
    ELSEIF(no_tracer==1)THEN      
      pot_temp          => ocean_state%p_diag%rho
    ENDIF
    
    start_level=1

    CALL calc_entries_mixing_tensor(patch_3d, ocean_state, param, &
                    & taper_diagonal_horz,           &
                    & taper_diagonal_vert_expl,      &
                    & taper_diagonal_vert_impl,      &
                    & taper_off_diagonal_horz,       &
                    & taper_off_diagonal_vert,       &
                    & tracer_index )
!    ENDIF   
             
    IF(no_tracer<=2)THEN

      flux_sum_horz(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2D(1)%alloc_cell_blocks)%x(1)=0.0_wp 
      flux_sum_horz(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2D(1)%alloc_cell_blocks)%x(2)=0.0_wp         
      flux_sum_horz(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2D(1)%alloc_cell_blocks)%x(3)=0.0_wp
      
      flux_sum_vert(:,:,:) =0.0_wp                
      flux_vert_temp(:,:,:)=0.0_wp  


      IF(tracer_index==1)THEN
        
!       write(0,*) "DerivTemperature_vert_center"
        tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradTemperature_horz_center
        tracer_gradient_vert_center     => ocean_state%p_aux%DerivTemperature_vert_center
      ELSEIF(tracer_index==2)THEN
        
        flux_vert_sal(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2D(1)%alloc_cell_blocks)=0.0_wp        

!       write(0,*) "DerivSalinity_vert_center"
        tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradSalinity_horz_center
        tracer_gradient_vert_center     => ocean_state%p_aux%DerivSalinity_vert_center
      ENDIF

      DO blockNo = all_cells%start_block, all_cells%end_block        
        CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)      
        DO cell_index = start_cell_index, end_cell_index
          DO level = 1, n_zlev
            flux_vec_horz_center(cell_index,level,blockNo)%x = 0.0_wp
          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block        
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)      
        DO cell_index = start_cell_index, end_cell_index
 
          !horizontal GMRedi Flux at top layer: with the tapering this is the just horizontal diffusion
          DO level = start_level, MIN(patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo),start_level)
            flux_vec_horz_center(cell_index,start_level,blockNo)%x &
              &=taper_diagonal_horz(cell_index,start_level,blockNo) &
              &*tracer_gradient_horz_vec_center(cell_index,start_level,blockNo)%x

!             flux_vert_center(cell_index,start_level,blockNo) &
!               &=taper_diagonal_vert_expl(cell_index,start_level,blockNo)&
!               &*tracer_gradient_vert_center(cell_index,start_level,blockNo)
!            ! the top level flux_vert_center will be filled from the second level, if it exists
             flux_vert_center(cell_index,start_level,blockNo) = 0.0_wp
          ENDDO

          DO level = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
          
            !horizontal GM-Redi Flux
            flux_vec_horz_center(cell_index,level,blockNo)%x &
              &=taper_diagonal_horz(cell_index,level,blockNo)*  &
              &tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x&
              !the second term vanishes if GM-Kappa=isoneutral diffusion
              &+taper_off_diagonal_horz(cell_index,level,blockNo)%x&
              &*tracer_gradient_vert_center(cell_index,level,blockNo)
              
            !vertical GM-Redi Flux: this is the part that is explicit in time
            !If namelist option "switch_off_diagonal_vert_expl=TRUE" the
            !first (diagonal) term is set to zero. Default option "switch_off_diagonal_vert_expl=FALSE"
            !Second term contains isoneutral and GM components
            !
            !For diagnostic purposes everything is handled exlicitely.
            flux_vert_center(cell_index,level,blockNo)= &
              &(taper_diagonal_vert_expl(cell_index,level,blockNo)&
              &+taper_diagonal_vert_impl(cell_index,level,blockNo))&
              &*tracer_gradient_vert_center(cell_index,level,blockNo)&
              &+&
              &Dot_Product(tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x,&
              &            taper_off_diagonal_vert(cell_index,level,blockNo)%x)
                
          END DO
          ! fill the top level flux_vert_center from the second level, if it exists
!           DO level = start_level+1, MIN(patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo),start_level+1)
!             flux_vert_center(cell_index,start_level,blockNo) = flux_vert_center(cell_index,start_level+1,blockNo)
!           END DO
        END DO                
      END DO
!ICON_OMP_END_DO_PARALLEL

!    !Map the (explicit) vertical tracer flux to the prism top (where the vertical divergence is calculated later)
!    CALL map_scalar_center2prismtop(patch_3d, flux_vert_center, op_coeff,GMredi_flux_vert)

      ! now we treat the vertical isoneutral flux that is discretized implicitely in time.
      ! 
      !1.) Interpolate the tapered coefficient for the vertical tracer flux from prism center to prism top: 
      !this is the diagonal part that is handled implicitely in time  
      !CALL map_scalar_center2prismtop( patch_3d, &
      !  &                              taper_diagonal_vert_impl,&
      !  &                              op_coeff,                &
      !  &                              ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit)        
      !  &                              mapped_vertical_diagonal_impl)!param%a_tracer_v(:,:,:, tracer_index))
      !
      !Do level=1,n_zlev
      !CALL dbg_print('Old vert coeff: A_v', param%a_tracer_v(:,level,:, tracer_index), this_mod_name, 4, patch_2D%cells%in_domain)
      !CALL dbg_print('Old vert coeff: A_v', ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit,&
      !& this_mod_name, 4, patch_2D%cells%in_domain)
      !End do
      !
      !2.) Here we combine the vertical GMRedicoefficient that is treated implicitely (mapped_vertical_diagonal_impl, this
      !term involves the slopes-squared times the isopycnal mixing coefficient and is potentially large, therefore
      !it is discretized implicitely) with the vertical mixing coefficient from the PP-scheme. 
      !We follow the approach in POP, where these two contributions are added
      !(see Reference manual POP, sect 5.1.3, in particular p. 41, after eq (150)).
      !
!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
      !DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block     
      !  CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)      
      !  DO cell_index = start_cell_index, end_cell_index
      !    DO level = start_level, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
      !      param%a_tracer_v(cell_index,level,blockNo, tracer_index) =                &
      !       & param%a_tracer_v(cell_index,level,blockNo, tracer_index) + &
      !       & ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit(cell_index,level,blockNo)
      !    END DO                  
      !  END DO                
      !END DO
!ICON_OMP_END_DO_PARALLEL
     !Do level=1,n_zlev
     ! CALL dbg_print('New vert coeff: A_v', param%a_tracer_v(:,level,:, tracer_index), this_mod_name, 4, patch_2D%cells%in_domain)
      !CALL dbg_print('New vert coeff: A_v', &
      !&ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit(:,level,:),&
      !& this_mod_name, 4, patch_2D%cells%in_domain)
     !END DO 
!    !Map the (explicit) vertical tracer flux to the prism top (where the vertical divergence is calculated later)
!    CALL map_scalar_center2prismtop(patch_3d, flux_vert_center, op_coeff,GMredi_flux_vert)

      ! now we treat the vertical isoneutral flux that is discretized implicitely in time.
      ! 
      !1.) Interpolate the tapered coefficient for the vertical tracer flux from prism center to prism top: 
      !this is the diagonal part that is handled implicitely in time  
      !CALL map_scalar_center2prismtop( patch_3d, &
      !  &                              taper_diagonal_vert_impl,&
      !  &                              op_coeff,                &
      !  &                              ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit)        
      !  &                              mapped_vertical_diagonal_impl)!param%a_tracer_v(:,:,:, tracer_index))
      !
      !Do level=1,n_zlev
      !CALL dbg_print('Old vert coeff: A_v', param%a_tracer_v(:,level,:, tracer_index), this_mod_name, 4, patch_2D%cells%in_domain)
      !CALL dbg_print('Old vert coeff: A_v', ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit,&
      !& this_mod_name, 4, patch_2D%cells%in_domain)
      !End do
      !
      !2.) Here we combine the vertical GMRedicoefficient that is treated implicitely (mapped_vertical_diagonal_impl, this
      !term involves the slopes-squared times the isopycnal mixing coefficient and is potentially large, therefore
      !it is discretized implicitely) with the vertical mixing coefficient from the PP-scheme. 
      !We follow the approach in POP, where these two contributions are added
      !(see Reference manual POP, sect 5.1.3, in particular p. 41, after eq (150)).
      !
! !ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
      !DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block     
      !  CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)      
      !  DO cell_index = start_cell_index, end_cell_index
      !    DO level = start_level, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
      !      param%a_tracer_v(cell_index,level,blockNo, tracer_index) =                &
      !       & param%a_tracer_v(cell_index,level,blockNo, tracer_index) + &
      !       & ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit(cell_index,level,blockNo)
      !    END DO                  
      !  END DO                
      !END DO
! !ICON_OMP_END_DO_PARALLEL
     !Do level=1,n_zlev
     ! CALL dbg_print('New vert coeff: A_v', param%a_tracer_v(:,level,:, tracer_index), this_mod_name, 4, patch_2D%cells%in_domain)
      !CALL dbg_print('New vert coeff: A_v', &
      !&ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit(:,level,:),&
      !& this_mod_name, 4, patch_2D%cells%in_domain)
     !END DO 
        
    ! use a vector communicator
    CALL sync_patch_array_mult(sync_c, patch_2D, 3, &
        & flux_vec_horz_center(:,:,:)%x(1), flux_vec_horz_center(:,:,:)%x(2), flux_vec_horz_center(:,:,:)%x(3))

! !ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,neutral_coeff, &
! !ICON_OMP  level) ICON_OMP_DEFAULT_SCHEDULE

!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,neutral_coeff,level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

      DO cell_index = start_cell_index, end_cell_index
        end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
        IF(end_level <= min_dolic) CYCLE

        
         IF (no_tracer>=2)THEN
         
           salinityColumn(1:end_level) = salinity(cell_index,1:end_level,blockNo)
           
           !Nonlinear EOS 
           IF(EOS_TYPE/=1)THEN
             neutral_coeff = calc_neutralslope_coeff_func_onColumn(         &
             & pot_temp(cell_index,1:end_level,blockNo), salinityColumn(1:end_level),  &
             & depth_cellinterface(cell_index,2:end_level+1,blockNo), end_level)
             
           !Linear EOS: slope coefficients are equal to EOS-coefficients
           ELSEIF(EOS_TYPE==1)THEN
              neutral_coeff(:,1) = LinearThermoExpansionCoefficient 
              neutral_coeff(:,2) = LinearHalineContractionCoefficient       
            ENDIF
          
         ELSEIF(no_tracer==1)THEN         
           neutral_coeff(:,1)=1.0         
         ENDIF 


         IF(tracer_index==1)THEN
           DO level = start_level+1, end_level-1
                       
            ocean_state%p_aux%diagnose_Redi_flux_temp(cell_index,level,blockNo)%x&
            &=-neutral_coeff(level,1)*flux_vec_horz_center(cell_index,level,blockNo)%x
            
            flux_sum_horz(cell_index,level,blockNo)%x&
            &= ocean_state%p_aux%diagnose_Redi_flux_temp(cell_index,level,blockNo)%x


            flux_vert_temp(cell_index,level,blockNo)&
            &=-neutral_coeff(level,1)*flux_vert_center(cell_index,level,blockNo)
            
            ocean_state%p_aux%diagnose_Redi_flux_vert(cell_index,level,blockNo)  = flux_vert_temp(cell_index,level,blockNo)
            
           END DO
         ELSEIF(tracer_index==2)THEN
          DO level = start_level+1, end_level-1         
           
           ocean_state%p_aux%diagnose_Redi_flux_sal(cell_index,level,blockNo)%x&
           &=neutral_coeff(level,2)*flux_vec_horz_center(cell_index,level,blockNo)%x  
                       
            flux_sum_horz(cell_index,level,blockNo)%x&
            &=ocean_state%p_aux%diagnose_Redi_flux_temp(cell_index,level,blockNo)%x&
            &+ ocean_state%p_aux%diagnose_Redi_flux_sal(cell_index,level,blockNo)%x
            !flux_sum_horz(cell_index,level,blockNo)%x&
            !&=-neutral_coeff(level,1)*ocean_state%p_aux%PgradTemperature_horz_center(cell_index,level,blockNo)%x&
            !&+neutral_coeff(level,2)*ocean_state%p_aux%PgradSalinity_horz_center(cell_index,level,blockNo)%x


            flux_vert_sal(cell_index,level,blockNo)&
            &=neutral_coeff(level,2)*flux_vert_center(cell_index,level,blockNo)

            flux_sum_vert(cell_index,level,blockNo)  = ocean_state%p_aux%diagnose_Redi_flux_vert(cell_index,level,blockNo)&
            &+flux_vert_sal(cell_index,level,blockNo)
            !flux_sum_vert(cell_index,level,blockNo)&
            !&=-neutral_coeff(level,1)*ocean_state%p_aux%DerivTemperature_vert_center(cell_index,level,blockNo)&
            !&+neutral_coeff(level,2)*ocean_state%p_aux%DerivSalinity_vert_center(cell_index,level,blockNo)
          END DO                        
         ENDIF
          
      END DO ! cell_index = start_cell_index, end_cell_index
    END DO  ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_DO_PARALLEL


   !Map the explicit horizontal tracer flux from cell centers to edges (where the horizontal divergence is calculated)
    !
    ! use a vector communicator
    CALL sync_patch_array_mult(sync_c, patch_2D, 3, &
        & flux_sum_horz(:,:,:)%x(1), flux_sum_horz(:,:,:)%x(2), flux_sum_horz(:,:,:)%x(3))
    !    
!    CALL map_cell2edges_3D( patch_3D,flux_sum_horz, GMredi_flux_horz, op_coeff)

    !Map the (explicit) vertical tracer flux to the prism top (where the vertical divergence is calculated later)
    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_center2prismtop(patch_3d, flux_sum_vert, op_coeff,GMredi_flux_vert)
    ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_center2prismtop_GM(patch_3d, flux_sum_vert, op_coeff,GMredi_flux_vert)    
    ENDIF
!    !---------DEBUG DIAGNOSTICS-------------------------------------------edges_in_domain
    IF(tracer_index==no_tracer)THEN



    CALL map_cell2edges_3D( patch_3D,flux_sum_horz, GMredi_flux_horz, op_coeff)
    Do level=start_level+1,n_zlev-1
    idt_src=1  ! output print level (1-5, fix)
    write(0,*)'level',level
    CALL dbg_print('diagnose_Redi:horz',GMredi_flux_horz(:,level,:),&
    & str_module, idt_src, in_subset=edges_in_domain)
    END DO
    
    Do level=start_level+1,n_zlev-1
    idt_src=1  ! output print level (1-5, fix)
    write(0,*)'level',level,tracer_index
    CALL dbg_print('diagnose_Redi:vert',GMredi_flux_vert(:,level,:),&
    & str_module, idt_src, in_subset=cells_in_domain)
    END DO
    ENDIF
    
    !---------------------------------------------------------------------
    
!      CALL sync_patch_array(sync_e, patch_2D, GMredi_flux_horz(:,:,:))
!      CALL sync_patch_array(sync_c, patch_2D, GMredi_flux_vert(:,:,:))
     
  ELSEIF( no_tracer>2)THEN
    CALL finish(TRIM('calc_GMRediflux'),&
    & 'diagnose_Redi_flux_balance beyond temperature and salinity is not impemented yet')
  ENDIF

   ! for debugging
   ! GMredi_flux_vert(:,:,:) = 0.0_wp

!IF(tracer_index==no_tracer)stop

  END SUBROUTINE diagnose_Redi_flux_balance
  !-------------------------------------------------------------------------

!---------------------------------------------------------------------------

END MODULE mo_ocean_GM_Redi


