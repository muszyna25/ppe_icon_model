#ifndef __NO_ICON_OCEAN__

!>
!! @file bgc.f90
!! @brief Main biogeochemical subroutine, called at each time step
!!
!! This subroutine calls all routines that calculate changes of pelagic biogeochemical 
!! tracers due to local processes (like photosythesis, heterotrophic
!! processes, N-fixation, and denitrification), the air-sea gas
!! exchange of carbon dioxide, oxygen, dinitrogen, and the
!! benthic processes. It further calls the computation of the vertical displacement of
!! particles
!!
!! called by mo_hydro_ocean_run:perform_ho_stepping
!!
!!
#include "icon_definitions.inc"

SUBROUTINE BGC_ICON(p_patch_3D, hamocc_ocean_state)

  USE mo_model_domain,        ONLY: t_patch,t_patch_3D

  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  USE mo_bgc_icon_comm,       ONLY: update_icon, update_bgc, hamocc_state, &
       &                            set_bgc_tendencies_output
  USE mo_dynamics_config,     ONLY: nold 
  USE mo_hamocc_nml,          ONLY: i_settling, l_cyadyn,l_bgc_check,io_stdo_bgc,l_implsed, &
       &                            l_dynamic_pi, l_pdm_settling 
  USE mo_control_bgc,         ONLY: ndtdaybgc,  &
       &                        ldtrunbgc, bgc_nproma

  USE mo_cyano,               ONLY: cyano, cyadyn
  USE mo_bgc_surface,         ONLY: gasex, update_weathering, dust_deposition, &
&                                   nitrogen_deposition, update_linage
  USE mo_bgc_bcond,           ONLY: ext_data_bgc
  USE mo_hamocc_diagnostics,  ONLY: get_inventories, get_omz
  USE mo_exception, ONLY: message
  USE mo_carchm,              ONLY: calc_dissol 
  USE mo_powach,              ONLY: powach, powach_impl
  USE mo_sedmnt, ONLY         : ini_bottom
  USE mo_timer, ONLY          : timer_bgc_up_bgc, timer_bgc_swr, timer_bgc_wea,timer_bgc_depo, &
    &                           timer_bgc_chemcon, timer_bgc_ocprod, timer_bgc_sett,timer_bgc_cya,&
    &                           timer_bgc_gx, timer_bgc_calc, timer_bgc_powach, timer_bgc_up_ic, &
    &                           timer_bgc_tend, timer_start, timer_stop, timers_level
  USE mo_settling,   ONLY    : settling, settling_pdm
  USE mo_ocean_hamocc_couple_state, ONLY: t_ocean_to_hamocc_state, t_hamocc_to_ocean_state, &
    & t_hamocc_ocean_state
!   USE mo_util_dbg_prnt,          ONLY: dbg_print

  IMPLICIT NONE

  TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d
  TYPE(t_hamocc_ocean_state), TARGET     :: hamocc_ocean_state

  ! Local variables
  INTEGER ::  jb
  INTEGER :: start_index, end_index
  INTEGER :: itrig_chemcon
  !INTEGER, POINTER :: levels(:)
  INTEGER :: levels(bgc_nproma)

  TYPE(t_subset_range), POINTER :: all_cells
  INTEGER :: alloc_cell_blocks
  TYPE(t_patch),POINTER    :: p_patch

  TYPE(t_ocean_to_hamocc_state), POINTER           :: ocean_to_hamocc_state
  TYPE(t_hamocc_to_ocean_state), POINTER           :: hamocc_to_ocean_state

  INTEGER :: i
  
!   CHARACTER(LEN=*), PARAMETER  :: str_module = 'BGC_ICON'  ! Output of module for 1 line debug
!   INTEGER :: idt_src
!   TYPE(t_subset_range), POINTER :: owned_cells

  ocean_to_hamocc_state => hamocc_ocean_state%ocean_to_hamocc_state
  hamocc_to_ocean_state => hamocc_ocean_state%hamocc_to_ocean_state


  !-----------------------------------------------------------------------
  p_patch   => p_patch_3D%p_patch_2d(1)
  alloc_cell_blocks =  p_patch%alloc_cell_blocks
  !--------------------------------------------
  !
  !----------------------------------------------------------------------
  !
  all_cells => p_patch%cells%ALL
!   owned_cells => p_patch%cells%owned

  
  ! trigger chemcon at depth only once per day
  itrig_chemcon=mod(ldtrunbgc,ndtdaybgc)+1
  !
  !
  !----------------------------------------------------------------------

  !---------DEBUG DIAGNOSTICS-------------------------------------------
!   idt_src=1  ! output print level (1-5, fix)
!   CALL dbg_print('h_old'           ,ocean_to_hamocc_state%h_old, str_module,idt_src, owned_cells)
!   CALL dbg_print('co2_mixing'      ,ocean_to_hamocc_state%co2_mixing_ratio, str_module,idt_src, owned_cells)
!   CALL dbg_print('short_wave_flux' ,ocean_to_hamocc_state%short_wave_flux, str_module,idt_src, owned_cells) 
!   CALL dbg_print('concSum'         ,ocean_to_hamocc_state%ice_concentration_sum, str_module,idt_src, owned_cells)
!   CALL dbg_print('salinity'        ,ocean_to_hamocc_state%salinity, str_module,idt_src, owned_cells)
!   CALL dbg_print('temperature'     ,ocean_to_hamocc_state%temperature, str_module,idt_src, owned_cells)
!   CALL dbg_print('wind10m'         ,ocean_to_hamocc_state%wind10m, str_module,idt_src, owned_cells)
 !----------------------------------------------------------------------
 
IF(l_bgc_check)THEN
 call message('before loop','inventories',io_stdo_bgc)
 call get_inventories(hamocc_state, ocean_to_hamocc_state%h_old, hamocc_state%p_prog(nold(1))%tracer, p_patch_3d) 
ENDIF

!DIR$ INLINE
  DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        !  tracer 1: potential temperature
        !  tracer 2: salinity
        levels(start_index:end_index) = p_patch_3D%p_patch_1d(1)%dolic_c(start_index:end_index,jb)

        start_detail_timer(timer_bgc_up_bgc,5)
        CALL update_bgc(start_index,end_index,levels,&
             & p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&  ! cell thickness
             &jb, hamocc_state%p_prog(nold(1))%tracer(:,:,jb,:), &
             & ocean_to_hamocc_state%co2_mixing_ratio(:,jb)                        & ! co2mixing ratio
             & ,hamocc_state%p_diag,hamocc_state%p_sed, hamocc_state%p_tend)
        stop_detail_timer(timer_bgc_up_bgc,5)

        CALL ini_bottom(start_index,end_index,levels,p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb))

        start_detail_timer(timer_bgc_swr,5)
       ! Net solar radiation update and swr_frac
        CALL swr_absorption(start_index,end_index,levels,                   &
 &                          ocean_to_hamocc_state%short_wave_flux(:,jb),                                 & ! SW radiation
 &                          ocean_to_hamocc_state%ice_concentration_sum(:,jb),                           & ! sea ice concentration
 &                          p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb))   ! level thickness

        stop_detail_timer(timer_bgc_swr,5)

       ! Linear age
        CALL update_linage(levels, start_index, end_index,  & ! index range, levels, salinity
   &                 p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb))! cell thickness (check for z0)

       ! Biogeochemistry

        start_detail_timer(timer_bgc_wea,5)
       ! Weathering fluxes 
        CALL update_weathering(start_index, end_index,  & ! index range, levels, salinity
   &                 p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&! cell thickness (check for z0)
   &                 ocean_to_hamocc_state%h_old(:,jb)) ! surface_height

        stop_detail_timer(timer_bgc_wea,5)

        start_detail_timer(timer_bgc_depo,5)
      ! Dust deposition
        CALL dust_deposition(start_index, end_index,  & ! index range, levels,
   &                 p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb), &! cell thickness (check for z0)
   &                 ocean_to_hamocc_state%h_old(:,jb),& ! surface_height
   &                 ext_data_bgc%dusty(:,jb))       ! dust input
      ! Nitrogen deposition
        CALL nitrogen_deposition(start_index, end_index,  & ! index range, levels
   &                 p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb), &! cell thickness (check for z0)
   &                 ocean_to_hamocc_state%h_old(:,jb),& ! surface_height
   &                 ext_data_bgc%nitro(:,jb))      ! nitrogen input

        stop_detail_timer(timer_bgc_depo,5)
       !----------------------------------------------------------------------
       ! Calculate chemical properties 

        start_detail_timer(timer_bgc_chemcon,5)
        CALL chemcon(start_index, end_index,levels,  ocean_to_hamocc_state%salinity(:,:,jb), & ! index range, levels, salinity
   &                 ocean_to_hamocc_state%temperature(:,:,jb),                              & ! pot. temperature
   &                 p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),                    & ! cell thickness
   &                 p_patch_3d%p_patch_1d(1)%depth_CellInterface(:,:,jb) ,  &           ! depths at interface  
   &                 itrig_chemcon)           
        stop_detail_timer(timer_bgc_chemcon,5)
       !----------------------------------------------------------------------
       ! Calculate plankton dynamics and particle settling 

       IF(i_settling.ne.2)then !2==agg

        start_detail_timer(timer_bgc_ocprod,5)
         ! plankton dynamics and remineralization  
         CALL ocprod(levels, start_index,end_index, ocean_to_hamocc_state%temperature(:,:,jb),&
   &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb), & ! cell thickness
   &               ocean_to_hamocc_state%h_old(:,jb),& ! surface height
   &               p_patch_3d%p_patch_1d(1)%depth_CellInterface(:,:,jb),& ! depths at interface  
   &               l_dynamic_pi ) ! depths at interface  
        stop_detail_timer(timer_bgc_ocprod,5)

        start_detail_timer(timer_bgc_sett,5)
         ! particle settling 
        IF (l_pdm_settling)then
         CALL settling_pdm(levels,start_index, end_index, &
   &                   p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb)) ! cell thickness
 
        ELSE
         CALL settling(levels,start_index, end_index, &
   &                   p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),& ! cell thickness
   &                   ocean_to_hamocc_state%h_old(:,jb))                   ! surface height
       
        ENDIF 
       endif

       stop_detail_timer(timer_bgc_sett,5)


       !----------------------------------------------------------------------
       ! Calculate N2 fixation 

       start_detail_timer(timer_bgc_cya,5)
       IF (l_cyadyn) THEN 
        ! dynamic cyanobacteria
        CALL cyadyn(levels, start_index,end_index, &  ! vertical range, cell range,
     &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),& ! cell thickness
     &               ocean_to_hamocc_state%h_old(:,jb), &                 ! surface height
     &               ocean_to_hamocc_state%temperature(:,:,jb), &        ! pot. temperature 
     &               p_patch_3d%p_patch_1d(1)%depth_CellInterface(:,:,jb),& ! depths at interface  
     &               l_dynamic_pi ) ! depths at interface  
       ELSE
        ! diagnostic N2 fixation
        CALL cyano (start_index, end_index,p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&
     &               ocean_to_hamocc_state%h_old(:,jb))                 ! surface height    
       endif
       stop_detail_timer(timer_bgc_cya,5)


       !----------------------------------------------------------------------
       ! Calculate gas exchange

        start_detail_timer(timer_bgc_gx,5)
        CALL gasex( start_index, end_index,    & 
  &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&  ! cell thickness
  &               ocean_to_hamocc_state%h_old(:,jb), &                   ! surface height
  &               ocean_to_hamocc_state%temperature(:,:,jb), &          ! pot. temperature 
  &               ocean_to_hamocc_state%salinity(:,:,jb), &          ! salinity
  &               ocean_to_hamocc_state%wind10m(:,jb)            , &          ! 10m wind speed 
  &               ocean_to_hamocc_state%ice_concentration_sum(:,jb))                              ! sea ice concentration

       stop_detail_timer(timer_bgc_gx,5)
        !----------------------------------------------------------------------
        ! Calculate carbonate dissolution
 
        start_detail_timer(timer_bgc_calc,5)
         CALL calc_dissol( start_index, end_index, levels,   & 
   &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),& ! cell thickness
   &               ocean_to_hamocc_state%salinity(:,:,jb))           ! salinity
 
        stop_detail_timer(timer_bgc_calc,5)
       !----------------------------------------------------------------------
        ! Calculate sediment dynamics
        start_detail_timer(timer_bgc_powach,5)
        if(l_implsed)then 
         CALL powach_impl( start_index, end_index,    & 
   &               ocean_to_hamocc_state%salinity(:,:,jb))          ! salinity
         else
 
        CALL powach( start_index, end_index, &    
   &               ocean_to_hamocc_state%salinity(:,:,jb),          &! salinity
   &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb))  ! cell thickness
         endif
        stop_detail_timer(timer_bgc_powach,5)

        if(mod(ldtrunbgc,ndtdaybgc).eq.0) CALL sedshi(start_index,end_index)
 
        start_detail_timer(timer_bgc_up_ic,5)
        CALL update_icon(start_index,end_index,levels,&
  &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&  ! cell thickness
  &               hamocc_state%p_prog(nold(1))%tracer(:,:,jb,:),            &
  &               hamocc_to_ocean_state%co2_flux(:,jb)                    )          ! co2flux for coupling
        stop_detail_timer(timer_bgc_up_ic,5)

        start_detail_timer(timer_bgc_tend,5)
        CALL set_bgc_tendencies_output(start_index,end_index,levels, &
  &               p_patch_3D%p_patch_1d(1)%prism_thick_c(:,:,jb),&  ! cell thickness
  &                                   jb, &
  &                                   hamocc_state%p_tend,            &
  &                                   hamocc_state%p_diag,            &
  &                                   hamocc_state%p_sed)

        stop_detail_timer(timer_bgc_tend,5)
 ENDDO

! O2 min depth & value diagnostics
CALL get_omz(hamocc_state,ocean_to_hamocc_state%h_old,p_patch_3d)

  ! Increment bgc time step counter of run (initialized in INI_BGC).
  !
  ldtrunbgc = ldtrunbgc + 1

IF(l_bgc_check)THEN
 call message('after loop','inventories',io_stdo_bgc)
 call get_inventories(hamocc_state, ocean_to_hamocc_state%h_old, hamocc_state%p_prog(nold(1))%tracer, p_patch_3d) 
ENDIF
  

END SUBROUTINE 

#endif
