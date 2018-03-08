#ifndef __NO_ICON_OCEAN__

SUBROUTINE INI_BGC_ICON(p_patch_3D, p_os,p_as,l_is_restart)

  USE mo_kind, ONLY           : wp, dp
  USE mo_hamocc_nml, ONLY     : l_cpl_co2, io_stdo_bgc
  USE mo_exception, ONLY      : message
  USE mo_control_bgc, ONLY    : dtb, dtbgc, inv_dtbgc, ndtdaybgc, icyclibgc,  &
       &                        ndtrunbgc, ldtrunbgc, bgc_zlevs, bgc_nproma, &
       &                        inv_dtb
!#ifdef AVFLUX
!  USE mo_avflux, ONLY         : avflux_ini
!#endif


  USE mo_memory_bgc, ONLY      : alloc_mem_biomod, n90depth,n1000depth,n2000depth,&
  &                              alloc_mem_carbch,totalarea
  USE mo_bgc_icon_comm, ONLY  : ini_bgc_regions, initial_update_icon, hamocc_state, &
      &                         print_bgc_parameters,print_wpoc, update_bgc

  USE mo_sedmnt, ONLY         : alloc_mem_sedmnt, ini_bottom, sediment_bottom
  USE mo_ini_bgc, ONLY        : ini_aquatic_tracers,            &
       &                        ini_pore_water_tracers,         &
       &                        ini_atmospheric_concentrations, &
       &                        set_parameters_bgc, &
       &                        ini_continental_carbon_input, &
       &                        ini_wpoc, bgc_param_conv_unit

  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_run_config,          ONLY: dtime
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_nml,           ONLY: n_zlev, nbgctra, no_tracer,dzlev_m
  USE mo_dynamics_config,     ONLY: nold,nnew
  USE mo_parallel_config,     ONLY: nproma
  USE mo_sync,                ONLY: global_sum_array
  USE mo_math_utilities,      ONLY: set_zlev
  USE mo_ocean_surface_types, ONLY: t_atmos_for_ocean 
 ! USE mo_hamocc_diagnostics,  ONLY: get_inventories

  IMPLICIT NONE

  !! Arguments
  TYPE(t_patch_3D),             TARGET,INTENT(IN)    :: p_patch_3D
  TYPE(t_hydro_ocean_state)                   :: p_os
  TYPE(t_atmos_for_ocean)                     :: p_as
  LOGICAL, INTENT(in):: l_is_restart

  INTEGER :: alloc_cell_blocks

  !! Local variables

  CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'ini_bgc_icon'
  TYPE(t_patch),POINTER    :: p_patch
  TYPE(t_subset_range), POINTER :: all_cells

  INTEGER :: jc, jk, jb
  INTEGER :: start_index, end_index
  INTEGER  :: levels(nproma)

  INTEGER, POINTER              :: regions(:,:)

  REAL(wp),ALLOCATABLE:: dlevels_m(:),dlevels_i(:)
 
  CALL message(TRIM(routine), 'start')

  !
  !----------------------------------------------------------------------
  !
  ! Set control constants ( mo_control_bgc )
  !
  dtbgc = dtime                         !  time step length [sec]
  inv_dtbgc = 1.0_wp / dtbgc
  ndtdaybgc = NINT(86400._wp / dtbgc) !  time steps per day [no.]
  dtb = 1._wp / REAL(ndtdaybgc, wp)   !  time step length [days]
  inv_dtb = 1._wp/ dtb


   ! determine size of arrays 
   p_patch           => p_patch_3D%p_patch_2D(1)
   alloc_cell_blocks =  p_patch%alloc_cell_blocks

   bgc_zlevs = n_zlev
   bgc_nproma = nproma

   all_cells => p_patch%cells%ALL

   regions => p_patch_3D%regio_c

  !
  ! Initialize time step counter of run.
  !
  ldtrunbgc = 0

  
  CALL message(TRIM(routine), 'set_parameters_bgc' )

  CALL set_parameters_bgc


  CALL print_bgc_parameters 

  CALL message(TRIM(routine), 'bgc_param_conv_unit')


  ! region indices
  CALL message(TRIM(routine), 'ini_bgc_regions' )

  CALL INI_BGC_REGIONS


  !
  ! Allocate memory : biology
  !
  CALL message(TRIM(routine), 'alloc_mem_biomod')
  CALL ALLOC_MEM_BIOMOD

  !
  ! Allocate memory : sediment
  !
  CALL message(TRIM(routine), 'alloc_mem_sedmnt' )
  CALL ALLOC_MEM_SEDMNT

  !
  ! Allocate memory : inorganic carbon cycle
  !
  CALL message(TRIM(routine), 'alloc_mem_carbch')
  CALL ALLOC_MEM_CARBCH

  !

  !
  ! Initialize sediment layering
  !
  CALL message(TRIM(routine), 'sediment_bottom')
  CALL sediment_bottom
  
  ! convert 1/d to 1/ts
  CALL bgc_param_conv_unit


  CALL message(TRIM(routine), 'ini weathering fluxes')
  totalarea = 0._wp
  DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        DO jc=start_index, end_index
           totalarea = totalarea + p_patch%cells%area(jc,jb) 
        ENDDO
  ENDDO
  totalarea     = global_sum_array(totalarea)
  CALL ini_continental_carbon_input(totalarea)

! set level for 90m, 1000m, 2000m for diagnostic output
  ALLOCATE(dlevels_m(n_zlev))
  ALLOCATE(dlevels_i(n_zlev+1))

  CALL set_zlev(dlevels_i, dlevels_m, n_zlev, dzlev_m)

  n90depth   = maxloc(dlevels_i,DIM=1,MASK=(dlevels_i < 90._wp))
  n1000depth = maxloc(dlevels_i,DIM=1,MASK=(dlevels_i < 1000._wp))
  n2000depth = maxloc(dlevels_i,DIM=1,MASK=(dlevels_i < 2000._wp))

! Initialize POC sinking speed
  CALL ini_wpoc(dlevels_i)
  CALL print_wpoc

  DEALLOCATE(dlevels_m)
  DEALLOCATE(dlevels_i)
   
!DIR$ INLINE
 DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
        !  tracer 1: potential temperature
        !  tracer 2: salinity
        
        levels(start_index:end_index) = p_patch_3d%p_patch_1d(1)%dolic_c(start_index:end_index,jb)



        IF(l_is_restart)then

          CALL update_bgc(start_index,end_index,levels,&
             & p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&  ! cell thickness
             &jb, p_os%p_prog(nold(1))%tracer(:,:,jb,:)&
             &,p_as%co2(:,jb)&
             & ,hamocc_state%p_diag,hamocc_state%p_sed, hamocc_state%p_tend)

          CALL ini_bottom(start_index,end_index,levels,p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb))

         ELSE

          CALL ini_aquatic_tracers(start_index,end_index,levels,regions(:,jb))

          CALL ini_bottom(start_index,end_index,levels,p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb))

          CALL ini_pore_water_tracers(start_index,end_index)

         ENDIF 

        CALL ini_atmospheric_concentrations

        CALL initial_update_icon(start_index,end_index,levels, &
              &               p_patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,jb),&  ! cell thickness
              &              jb, &
              &              p_os%p_prog(nold(1))%tracer(:,:,jb,:),&
              &              hamocc_state%p_sed, hamocc_state%p_diag,&
              &              p_as%co2flx(:,jb))
     

  ENDDO

!  IF (l_cpl_co2) THEN
!     !     initializes fields used for redistribution of co2 fluxes
!     CALL avflux_ini
!  END IF

 CALL message(TRIM(routine), 'end ini bgc')
END SUBROUTINE 

#endif
