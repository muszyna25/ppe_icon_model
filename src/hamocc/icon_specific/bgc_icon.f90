!>
!! @file bgc.f90
!! @brief Main biogeochemical subroutine, called at each time step
!!
!! This subroutine computes all changes of pelagic biogeochemical 
!! tracers due to local processes (like photosythesis, heterotrophic
!! processes, N-fixation, and denitrification), the air-sea gas
!! exchange of carbon dioxide, oxygen, dinitrogen, and the
!! benthic processes. It further computes the vertical displacement of
!! particles
!!
!! called by mo_hydro_ocean_run:perform_ho_stepping
!!
!! @author Irene Stemmler, MPI-Met, HH
!!
!! @par Revision History
!!
SUBROUTINE BGC_ICON(p_patch_3D, p_os, p_as, p_ice)

  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: n_zlev, nbgctra, no_tracer, l_partial_cells
  USE mo_model_domain,        ONLY: t_patch,t_patch_3D

  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  USE mo_bgc_icon_comm,       ONLY: update_icon, update_bgc, hamocc_state, &
       &                            set_bgc_tendencies_output
  USE mo_dynamics_config,     ONLY: nold, nnew
  USE mo_sea_ice_types,       ONLY: t_atmos_for_ocean, t_sea_ice ! for now, use different later
  USE mo_hamocc_nml,          ONLY: i_settling, l_cyadyn,l_bgc_check,io_stdo_bgc,l_implsed 
  USE mo_control_bgc,         ONLY: dtb, dtbgc, inv_dtbgc, ndtdaybgc, icyclibgc,  &
       &                        ndtrunbgc, ldtrunbgc, bgc_zlevs, bgc_nproma

  USE mo_cyano,               ONLY: cyano, cyadyn
  USE mo_bgc_surface,         ONLY: gasex, update_weathering, dust_deposition
  USE mo_bgc_bcond,           ONLY: ext_data_bgc
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_hamocc_diagnostics,  ONLY: get_monitoring, get_inventories
  USE mo_exception, ONLY: message
  USE mo_carchm,              ONLY: calc_dissol 
  USE mo_sedmnt, ONLY         : ini_bottom

  IMPLICIT NONE

  TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d
  TYPE(t_hydro_ocean_state)              :: p_os
  TYPE(t_atmos_for_ocean)                :: p_as
  TYPE(t_sea_ice)                        :: p_ice

  ! Local variables
  INTEGER :: jc, jk, jb, itest
  INTEGER :: start_index, end_index
  INTEGER :: itrig_chemcon
  !INTEGER, POINTER :: levels(:)
  INTEGER :: levels(bgc_nproma)

  TYPE(t_subset_range), POINTER :: all_cells
  INTEGER :: alloc_cell_blocks
  TYPE(t_patch),POINTER    :: p_patch



  !-----------------------------------------------------------------------
  p_patch   => p_patch_3D%p_patch_2d(1)
  alloc_cell_blocks =  p_patch%alloc_cell_blocks
  !--------------------------------------------
  !
  !----------------------------------------------------------------------
  !
  all_cells => p_patch%cells%ALL

  
  ! trigger chemcon at depth only once per day
  itrig_chemcon=mod(ldtrunbgc,ndtdaybgc)+1
  !
  !
  !----------------------------------------------------------------------

 
IF(l_bgc_check)THEN
 call message('before loop','inventories',io_stdo_bgc)
! call get_inventories(hamocc_state, p_os, p_patch_3d)
ENDIF
  !CALL dbg_print('h new   ',p_os%p_prog(nnew(1))%h,       'bgc_icon', 1, in_subset=p_patch%cells%owned)
 ! CALL dbg_print('thick_c_flat   ', p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,1,:),       'bgc_icon', 1, in_subset=p_patch%cells%owned)


  DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        !  tracer 1: potential temperature
        !  tracer 2: salinity
        levels(start_index:end_index) = p_patch_3D%p_patch_1d(1)%dolic_c(start_index:end_index,jb)


        CALL ini_bottom(start_index,end_index,levels,p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb))

         CALL update_bgc(start_index,end_index,levels,&
             & p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&  ! cell thickness
             &jb, p_os%p_prog(nold(1))%tracer(:,:,jb,:)&
             & ,hamocc_state%p_diag,hamocc_state%p_sed, hamocc_state%p_tend)


       !
       ! Net solar radiation update and swr_frac
        CALL swr_absorption(start_index,end_index,levels,                   &
 &                          p_as%fswr(:,jb),                                & ! SW radiation
 &                          p_ice%concSum(:,jb),                            & ! sea ice concentration
 &                          p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb))   ! level thickness


       ! Biogeochemistry

       ! Weathering fluxes 
        CALL update_weathering(start_index, end_index,  & ! index range, levels, salinity
   &                 p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&! cell thickness (check for z0)
   &                 p_os%p_prog(nold(1))%h(:,jb)) ! surface_height

      ! Dust deposition
        CALL dust_deposition(start_index, end_index,  & ! index range, levels, salinity
   &                 p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb), &! cell thickness (check for z0)
   &                 p_os%p_prog(nold(1))%h(:,jb),& ! surface_height
   &                 ext_data_bgc%dusty(:,jb))       ! dust input

       !----------------------------------------------------------------------
       ! Calculate chemical properties 

        CALL chemcon(start_index, end_index,levels,  p_os%p_prog(nold(1))%tracer(:,:,jb,2), & ! index range, levels, salinity
   &                 p_os%p_prog(nold(1))%tracer(:,:,jb,1),                              & ! pot. temperature
   &                 p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),                    & ! cell thickness
   &                 p_patch_3d%p_patch_1d(1)%depth_CellInterface(:,:,jb) ,  &           ! depths at interface  
   &                 itrig_chemcon)           
       !----------------------------------------------------------------------
       ! Calculate plankton dynamics and particle settling 

       IF(i_settling.ne.2)then !2==agg

         ! plankton dynamics and remineralization  
         CALL ocprod(levels, start_index,end_index, p_os%p_prog(nold(1))%tracer(:,:,jb,1),&
   &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb), & ! cell thickness
   &               p_os%p_prog(nold(1))%h(:,jb)) ! surface height

         ! particle settling
         CALL settling(levels,start_index, end_index, &
   &                   p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),& ! cell thickness
   &                   p_os%p_prog(nold(1))%h(:,jb))                   ! surface height

        endif


       !----------------------------------------------------------------------
       ! Calculate N2 fixation 

       IF (l_cyadyn) THEN 
        ! dynamic cyanobacteria
        CALL cyadyn(levels, start_index,end_index, &  ! vertical range, cell range,
     &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),& ! cell thickness
     &               p_os%p_prog(nold(1))%h(:,jb), &                 ! surface height
     &               p_os%p_prog(nold(1))%tracer(:,:,jb,1) )          ! pot. temperature 
       ELSE
        ! diagnostic N2 fixation
        CALL cyano (start_index, end_index,p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&
     &               p_os%p_prog(nold(1))%h(:,jb))                 ! surface height    
       endif


       !----------------------------------------------------------------------
       ! Calculate gas exchange

        CALL gasex( start_index, end_index,    & 
  &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&  ! cell thickness
  &               p_os%p_prog(nold(1))%h(:,jb), &                   ! surface height
  &               p_os%p_prog(nold(1))%tracer(:,:,jb,1), &          ! pot. temperature 
  &               p_os%p_prog(nold(1))%tracer(:,:,jb,2), &          ! salinity
  &               p_as%fu10(:,jb)                      , &          ! 10m wind speed 
  &               p_ice%concSum(:,jb))                              ! sea ice concentration

        !----------------------------------------------------------------------
        ! Calculate carbonate dissolution
 
         CALL calc_dissol( start_index, end_index, levels,   & 
   &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),& ! cell thickness
   &               p_os%p_prog(nold(1))%tracer(:,:,jb,2))           ! salinity
 
       !----------------------------------------------------------------------
        ! Calculate sediment dynamics
        if(l_implsed)then 
         CALL powach_impl( start_index, end_index, levels,   & 
   &               p_os%p_prog(nold(1))%tracer(:,:,jb,2))          ! salinity
         else
 
        CALL powach( start_index, end_index, &    
   &               p_os%p_prog(nold(1))%tracer(:,:,jb,2),          &! salinity
   &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb))  ! cell thickness
         endif

        if(mod(ldtrunbgc,ndtdaybgc).eq.0) CALL sedshi(start_index,end_index)
 
        CALL update_icon(start_index,end_index,levels,&
  &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&  ! cell thickness
  &               p_os%p_prog(nold(1))%tracer(:,:,jb,:))
 
        CALL set_bgc_tendencies_output(start_index,end_index,levels, &
  &               p_patch_3D%p_patch_1d(1)%prism_thick_c(:,:,jb),&  ! cell thickness
  &                                   jb, &
  &                                   hamocc_state%p_tend,            &
  &                                   hamocc_state%p_diag,            &
  &                                   hamocc_state%p_sed)

 ENDDO
  ! Increment bgc time step counter of run (initialized in INI_BGC).
  !
  ldtrunbgc = ldtrunbgc + 1

IF(l_bgc_check)THEN
 call message('after loop','inventories',io_stdo_bgc)
! call get_inventories(hamocc_state, p_os, p_patch_3d)
ENDIF
  

END SUBROUTINE 
