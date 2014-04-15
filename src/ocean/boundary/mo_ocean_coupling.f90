!>
!! @par Revision History
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_ocean_coupling
  
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: message, finish
  USE mo_impl_constants,      ONLY: success, max_char_length
  USE mo_physical_constants,  ONLY: tmelt, rho_ref
  USE mo_master_control,      ONLY: is_restart_run, get_my_process_name, get_my_model_no
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs , num_restart_procs
  USE mo_mpi,                 ONLY: my_process_is_io,set_mpi_work_communicators,p_pe_work, process_mpi_io_size
  USE mo_grid_config,         ONLY: n_dom
  USE mo_datetime,            ONLY: t_datetime
  USE mo_time_config,         ONLY: time_config
  USE mo_run_config,          ONLY: ltimer
  USE mo_dynamics_config,     ONLY: nold
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling
  USE mo_sync,                ONLY: sync_c, sync_patch_array
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d, p_patch_local_parent
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  
  USE mo_oce_types
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean
  
  !-------------------------------------------------------------
  ! For the coupling
#ifndef __NO_ICON_ATMO__
# ifdef YAC_coupling
  USE mo_parallel_config,     ONLY: nproma
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE finterface_description, ONLY: yac_finit, yac_fdef_comp,                    &
    & yac_fdef_subdomain, yac_fconnect_subdomains, &
    & yac_fdef_elements, yac_fdef_points,          &
    & yac_fdef_mask, yac_fdef_field, yac_fsearch,  &
    & yac_ffinalize, yac_fput, yac_fget, yac_fget_nbr_fields, yac_fget_field_ids
  USE mo_coupling_config,     ONLY: is_coupled_run
# else
  USE mo_icon_cpl_init,       ONLY: icon_cpl_init
  USE mo_icon_cpl_init_comp,  ONLY: icon_cpl_init_comp
  USE mo_coupling_config,     ONLY: is_coupled_run, config_debug_coupler_level
  USE mo_icon_cpl_def_grid,   ONLY: icon_cpl_def_grid, icon_cpl_def_location
  USE mo_icon_cpl_def_field,  ONLY: icon_cpl_def_field, icon_cpl_get_nbr_fields, icon_cpl_get_field_ids
  USE mo_icon_cpl_search,     ONLY: icon_cpl_search
  USE mo_icon_cpl_finalize,   ONLY: icon_cpl_finalize
  USE mo_icon_cpl_restart,    ONLY: icon_cpl_write_restart
  USE mo_icon_cpl_exchg,      ONLY: icon_cpl_put, icon_cpl_get
# endif
#endif
  !-------------------------------------------------------------
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: construct_ocean_coupling, destruct_ocean_coupling
  PUBLIC :: couple_ocean_toatmo_fluxes
  
  CHARACTER(LEN=12)  :: module_name    = 'ocean_coupli'
  
CONTAINS
  
#ifdef __NO_ICON_ATMO__
  ! ---------------------------------------------------
  ! Dummy routines for compiling without atmo
  SUBROUTINE construct_ocean_coupling(patch_3d)
    TYPE(t_patch_3d ), TARGET, INTENT(in)    :: patch_3d
    RETURN
  END SUBROUTINE construct_ocean_coupling
  
  SUBROUTINE destruct_ocean_coupling()
    RETURN
  END SUBROUTINE destruct_ocean_coupling
  
  SUBROUTINE couple_ocean_toatmo_fluxes(patch_3d, ocean_state, atmos_for_ocean, ice, atmos_fluxes, &
    & surface_fluxes, jstep, datetime)
    TYPE(t_patch_3d ),TARGET, INTENT(in)        :: patch_3d
    TYPE(t_hydro_ocean_state)                   :: ocean_state
    TYPE(t_atmos_for_ocean)                     :: atmos_for_ocean
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_sea_ice)                             :: ice
    TYPE(t_sfc_flx)                             :: surface_fluxes
    INTEGER, INTENT(in)                         :: jstep
    TYPE(t_datetime), INTENT(inout)             :: datetime
  END SUBROUTINE couple_ocean_toatmo_fluxes
  
  !  SUBROUTINE init_coupled_ocean(patch_2d, ocean_state)
  !    TYPE(t_patch)                     :: patch_2d
  !    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
  !    RETURN
  !  END SUBROUTINE init_coupled_ocean
  ! ---------------------------------------------------
#else
  
  !--------------------------------------------------------------------------
  ! Prepare the coupling
  !
  ! For the time being this could all go into a subroutine which is
  ! common to atmo and ocean. Does this make sense if the setup deviates
  ! too much in future.
  !------------------------------------------------------------------
  SUBROUTINE construct_ocean_coupling(patch_3d)
    TYPE(t_patch_3d ), TARGET, INTENT(in)    :: patch_3d
    
    INTEGER, PARAMETER :: no_of_fields = 10
    CHARACTER(LEN=max_char_length) ::  field_name(no_of_fields)
    INTEGER :: field_id(no_of_fields)
    INTEGER :: grid_id
    INTEGER :: grid_shape(2)
    INTEGER :: field_shape(3)
    INTEGER :: i, error_status
    
    INTEGER :: patch_no
    
# ifdef YAC_coupling
    INTEGER, PARAMETER :: nbr_vertices_per_cell = 3 ! Triangle
    
    INTEGER, PARAMETER :: nbr_subdomain_ids = 1
    
    INTEGER :: comp_id
    INTEGER :: cell_point_id
    INTEGER :: edge_point_id
    INTEGER :: mask_id
    INTEGER :: subdomain_id
    INTEGER :: subdomain_ids(nbr_subdomain_ids)
    
    INTEGER, PARAMETER :: cell     = 0 ! one point per cell
    INTEGER, PARAMETER :: corner   = 1 ! one point per vertex
    INTEGER, PARAMETER :: edge     = 2 ! one point per edge
    ! (see definition of enum location in points.h)
    INTEGER :: jb, jc, je, INDEX
    INTEGER :: cell_start_idx, cell_end_idx
    INTEGER :: edge_start_idx, edge_end_idx
    
    REAL(wp), ALLOCATABLE :: buffer_x(:)
    REAL(wp), ALLOCATABLE :: buffer_y(:)
    REAL(wp), ALLOCATABLE :: buffer_c(:)
    
    TYPE(t_subset_range), POINTER :: all_cells, all_edges
    TYPE(t_patch), POINTER :: patch_horz
    
    IF (.NOT. is_coupled_run()) RETURN
    
    patch_no = 1
    
    ! Initialise the coupler
    CALL yac_finit ( "couling.xml", "coupling.xsd" )
    
    ! Inform the coupler about what we are
    CALL yac_fdef_comp ( "ICON_ocean", comp_id )
    
    ! Announce one subdomain (patch) to the coupler
    CALL yac_fdef_subdomain ( comp_id, "ICON_ocean", subdomain_id )
    
    patch_horz => patch_3d%p_patch_2d(patch_no)
    all_cells  => patch_horz%cells%ALL
    all_edges  => patch_horz%edges%ALL
    
    ! Extract cell information
    !
    ! cartesian coordinates of cell vertices are stored in
    ! patch_horz%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.
    
    ALLOCATE(buffer_x(nproma*(all_cells%end_block-all_cells%start_block+1)  ))
    ALLOCATE(buffer_y(nproma*(all_cells%end_block-all_cells%start_block+1)  ))
    ALLOCATE(buffer_c(nproma*(all_cells%end_block-all_cells%start_block+1)*3))
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, cell_start_idx, cell_end_idx)
      DO jc = cell_start_idx, cell_end_idx
        INDEX = (jb-1)*nproma+jc
        buffer_x(INDEX) = patch_horz%verts(jc,jb)%vertex%lon
        buffer_y(INDEX) = patch_horz%verts(jc,jb)%vertex%lat
        buffer_c((INDEX-1)*3+1) = patch_horz%cells%vertex_idx(jc,jb,1)
        buffer_c((INDEX-1)*3+2) = patch_horz%cells%vertex_idx(jc,jb,2)
        buffer_c((INDEX-1)*3+3) = patch_horz%cells%vertex_idx(jc,jb,3)
      ENDDO
    ENDDO
    
    ! Description of elements, here as unstructured grid
    CALL yac_fdef_elements ( subdomain_id,              &
      & patch_horz%n_patch_verts,  &
      & patch_horz%n_patch_cells,  &
      & nbr_vertices_per_cell,     &
      & buffer_x,                  &
      & buffer_y,                  &
      & buffer_c )
    
    ! Can we have two fdef_point calls for the same subdomain, i.e.
    ! one single set of cells?
    !
    ! Define cell center points (location = 0)
    !
    ! cartesian coordinates of cell centers are stored in
    ! patch_horz%cells%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.
    
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, cell_start_idx, cell_end_idx)
      DO jc = cell_start_idx, cell_end_idx
        INDEX = (jb-1)*nproma+jc
        buffer_x(INDEX) = patch_horz%cells%center(jc,jb)%lon
        buffer_x(INDEX) = patch_horz%cells%center(jc,jb)%lat
      ENDDO
    ENDDO
    
    CALL yac_fdef_points ( subdomain_id,            &
      & patch_horz%n_patch_cells,   &
      & cell,                    &
      & buffer_x,                &
      & buffer_y,                &
      & cell_point_id )
    
    ! Define edge center points (location = 2)
    !
    ! cartesian coordinates of cell centers are stored in
    ! patch_horz%edges%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.
    
    DEALLOCATE (buffer_x, buffer_y)
    
    ALLOCATE(buffer_x(nproma*(all_edges%end_block-all_edges%start_block+1)))
    ALLOCATE(buffer_y(nproma*(all_edges%end_block-all_edges%start_block+1)))
    
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, edge_start_idx, edge_end_idx)
      DO je = edge_start_idx, edge_end_idx
        INDEX = (jb-1)*nproma+je
        buffer_x(INDEX) = patch_horz%edges%center(jc,jb)%lon
        buffer_x(INDEX) = patch_horz%edges%center(jc,jb)%lat
      ENDDO
    ENDDO
    
    CALL yac_fdef_points ( subdomain_id,             &
      & patch_horz%n_patch_cells, &
      & edge,                     &
      & buffer_x,                 &
      & buffer_y,                 &
      & edge_point_id )
    
    ! Connect subdomains
    CALL yac_fconnect_subdomains ( comp_id,           &
      & nbr_subdomain_ids, &
      & subdomain_ids,     &
      & domain_id )
    !
    ! mask generation : ... not yet defined ...
    !
    ! We could use the patch_horz%cells%decomp_info%owner_local information
    ! e.g. to mask out halo points. We do we get the info about what is local and what
    ! is remote.
    !
    ! The land-sea mask for the ocean is available in patch_3D%surface_cell_sea_land_mask(:,:)
    !
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! CALL yac_fdef_mask ( mask_size,     &
    !                      imask,         &
    !                      cell_point_id, &
    !                      mask_id )
    
    CALL yac_fdef_mask ( mask_size,  &  !rr TODO
      & imask,      &  !rr TODO
      & points_id,  &
      & mask_id )
    
    DEALLOCATE (buffer_x, buffer_y, buffer_c)
    
# else
    IF (.NOT. is_coupled_run()) RETURN
    
    !------------------------------------------------------------
    CALL icon_cpl_init(debug_level=config_debug_coupler_level)
    ! Inform the coupler about what we are
    CALL icon_cpl_init_comp ( get_my_process_name(), get_my_model_no(), error_status )
    ! split the global_mpi_communicator into the components
    !------------------------------------------------------------
    patch_no      = 1
    
    grid_shape(1) = 1
    grid_shape(2) = patch_3d%p_patch_2d(patch_no)%n_patch_cells
    
    CALL icon_cpl_def_grid ( &
      & grid_shape, patch_3d%p_patch_2d(patch_no)%cells%decomp_info%glb_index, & ! input
      & grid_id, error_status )                          ! output
    
    ! Marker for internal and halo points, a list which contains the
    ! rank where the native cells are located.
    CALL icon_cpl_def_location ( &
      & grid_id, grid_shape, patch_3d%p_patch_2d(patch_no)%cells%decomp_info%owner_local, & ! input
      & p_pe_work,  & ! this owner id
      & error_status )                                            ! output
    
# endif
    
    field_name(1) = "TAUX"   ! bundled field containing two components
    field_name(2) = "TAUY"   ! bundled field containing two components
    field_name(3) = "SFWFLX" ! bundled field containing three components
    field_name(4) = "SFTEMP"
    field_name(5) = "THFLX"  ! bundled field containing four components
    field_name(6) = "ICEATM" ! bundled field containing four components
    field_name(7) = "SST"
    field_name(8) = "OCEANU"
    field_name(9) = "OCEANV"
    field_name(10) = "ICEOCE" ! bundled field containing four components
    
# ifdef YAC_coupling
    DO i = 1, no_of_fields
      CALL yac_fdef_field ( field_name(i),            &
        & comp_id,                  &
        & domain_id,                &
        & point_id,                 &
        & mask_id,                  &
        & patch_horz%n_patch_cells, &
        & field_id(i) )
    ENDDO
    
    CALL yac_fsearch ( nbr_components, comp_id, no_of_fields, field_id, error_status )
# else
    
    field_shape(1:2) = grid_shape(1:2)
    
    DO i = 1, no_of_fields
      IF ( i == 1 .OR. i == 2 ) THEN
        field_shape(3) = 2
      ELSE IF ( i == 3 ) THEN
        field_shape(3) = 3
      ELSE IF ( i == 6 .OR. i == 5 ) THEN
        field_shape(3) = 4
      ELSE IF ( i == 10 ) THEN
        field_shape(3) = 5
      ELSE
        field_shape(3) = 1
      ENDIF
      CALL icon_cpl_def_field ( field_name(i), grid_id, field_id(i), &
        & field_shape, error_status )
    ENDDO
    
    CALL icon_cpl_search
#endif
  END SUBROUTINE construct_ocean_coupling
  !--------------------------------------------------------------------------
  
  
  !--------------------------------------------------------------------------
  SUBROUTINE destruct_ocean_coupling()
# ifdef YAC_coupling
    IF ( is_coupled_run() ) CALL yac_ffinalize
# else
    IF ( is_coupled_run() ) CALL icon_cpl_finalize ()
# endif
  END SUBROUTINE destruct_ocean_coupling
  !--------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------
  SUBROUTINE couple_ocean_toatmo_fluxes(patch_3d, ocean_state, atmos_for_ocean, ice, atmos_fluxes, &
    & surface_fluxes, jstep, datetime)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in)        :: patch_3d
    TYPE(t_hydro_ocean_state)                   :: ocean_state
    TYPE(t_atmos_for_ocean)                     :: atmos_for_ocean
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_sea_ice)                             :: ice
    TYPE(t_sfc_flx)                             :: surface_fluxes
    INTEGER, INTENT(in)                         :: jstep
    TYPE(t_datetime), INTENT(inout)             :: datetime
    !
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = 'couple_ocean_toatmo_fluxes'
    INTEGER :: jmon, jdmon, jmon1, jmon2, ylen, yday
    INTEGER :: iniyear, curyear, offset
    INTEGER :: jc, jb, no_set
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp) :: z_tmin, z_relax, rday1, rday2, dtm1, dsec, z_smax, z_forc_tracer_old
    REAL(wp) ::  z_c2(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) ::   tfw(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), POINTER :: t_top(:,:), s_top(:,:)
    
    ! Local declarations for coupling:
    LOGICAL :: write_coupler_restart
    INTEGER :: info, ierror   !< return values form cpl_put/get calls
    INTEGER :: nbr_hor_points ! = inner and halo points
    INTEGER :: nbr_points     ! = nproma * nblks
    INTEGER :: nbr_fields
    INTEGER, ALLOCATABLE :: field_id(:)
    INTEGER :: field_shape(3)
    REAL(wp), ALLOCATABLE :: buffer(:,:)
    REAL(wp), PARAMETER :: seconds_per_month = 2.592e6_wp !TODO: use real month lenght
    TYPE(t_patch), POINTER:: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain
    INTEGER :: idt_src  ! Level of detail for 1 line debug
    
    
    IF (.NOT. is_coupled_run() ) RETURN
    
    IF (ltimer) CALL timer_start(timer_coupling)
    
    patch_2d   => patch_3D%p_patch_2D(1)
    time_config%cur_datetime = datetime
    
    nbr_hor_points = patch_2d%n_patch_cells
    nbr_points     = nproma * patch_2d%nblks_c
    ALLOCATE(buffer(nbr_points,5))
    buffer(:,:) = 0.0_wp
    
    !
    !  see drivers/mo_ocean_model.f90:
    !
    !   field_id(1) represents "TAUX"   wind stress component
    !   field_id(2) represents "TAUY"   wind stress component
    !   field_id(3) represents "SFWFLX" surface fresh water flux
    !   field_id(4) represents "SFTEMP" surface temperature
    !   field_id(5) represents "THFLX"  total heat flux
    !   field_id(6) represents "ICEATM" ice temperatures and melt potential
    !
    !   field_id(7) represents "SST"    sea surface temperature
    !   field_id(8) represents "OCEANU" u component of ocean surface current
    !   field_id(9) represents "OCEANV" v component of ocean surface current
    !   field_id(10)represents "ICEOCE" ice thickness, concentration and temperatures
    !
    !
#ifdef YAC_Coupling
    CALL yac_fget_nbr_fields ( nbr_fields )
    ALLOCATE(field_id(nbr_fields))
    CALL yac_fget_field_ids ( nbr_fields, field_id )
#else
    CALL icon_cpl_get_nbr_fields ( nbr_fields )
    ALLOCATE(field_id(nbr_fields))
    CALL icon_cpl_get_field_ids ( nbr_fields, field_id )
#endif
    !
    field_shape(1) = 1
    field_shape(2) = patch_2d%n_patch_cells
    field_shape(3) = 1
    
    !
    ! buffer is allocated over nproma only
    
    !
    ! Send fields from ocean to atmosphere
    ! ------------------------------------
    !
    write_coupler_restart = .FALSE.
    !
    ! SST
    buffer(:,1) = RESHAPE(ocean_state%p_prog(nold(1))%tracer(:,1,:,1), (/nbr_points /) ) &
      & + tmelt
    
#ifdef YAC_coupling
    CALL yac_fput ( field_id(7), nbr_hor_points, 1, 1, 1, buffer, ierror )
#else
    CALL icon_cpl_put ( field_id(7), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
    !
    ! zonal velocity
    buffer(:,1) = RESHAPE(ocean_state%p_diag%u(:,1,:), (/nbr_points /) )
#ifdef YAC_coupling
    CALL yac_fput ( field_id(8), nbr_hor_points, 1, 1, 1, buffer, ierror )
#else
    CALL icon_cpl_put ( field_id(8), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
    !
    ! meridional velocity
    buffer(:,1) = RESHAPE(ocean_state%p_diag%v(:,1,:), (/nbr_points /) )
#ifdef YAC_coupling
    CALL yac_fput ( field_id(9), nbr_hor_points, 1, 1, 1, buffer, ierror )
#else
    CALL icon_cpl_put ( field_id(9), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
    !
    ! Ice thickness, concentration, T1 and T2
    buffer(:,1) = RESHAPE(ice%hi  (:,1,:), (/nbr_points /) )
    buffer(:,2) = RESHAPE(ice%hs  (:,1,:), (/nbr_points /) )
    buffer(:,3) = RESHAPE(ice%conc(:,1,:), (/nbr_points /) )
    buffer(:,4) = RESHAPE(ice%t1  (:,1,:), (/nbr_points /) )
    buffer(:,5) = RESHAPE(ice%t2  (:,1,:), (/nbr_points /) )
    field_shape(3) = 5
#ifdef YAC_coupling
    CALL yac_fput ( field_id(10), nbr_hor_points, 4, 1, 1, buffer, ierror )
#else
    CALL icon_cpl_put ( field_id(10), field_shape, buffer(1:nbr_hor_points,1:5), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
    
    IF ( write_coupler_restart ) CALL icon_cpl_write_restart ( 4, field_id(7:10), ierror )
    !
    ! Receive fields from atmosphere
    ! ------------------------------
    
    !
    ! Apply wind stress - records 0 and 1 of field_id
    
    ! zonal wind stress
    field_shape(3) = 2
#ifdef YAC_coupling
    CALL yac_fget ( field_id(1), nbr_hor_points, 1, 1, 1, buffer, info, ierror )
#else
    CALL icon_cpl_get ( field_id(1), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:field_shape(3)) = 0.0_wp
      atmos_fluxes%stress_xw(:,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )
      atmos_fluxes%stress_x (:,:) = RESHAPE(buffer(:,2),(/ nproma, patch_2d%nblks_c /) )
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%stress_xw(:,:))
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%stress_x (:,:))
    ENDIF
    !
    ! meridional wind stress
#ifdef YAC_coupling
    CALL yac_fget ( field_id(2), nbr_hor_points, 1, 1, 1, buffer, info, ierror )
#else
    CALL icon_cpl_get ( field_id(2), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:field_shape(3)) = 0.0_wp
      atmos_fluxes%stress_yw(:,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )
      atmos_fluxes%stress_y (:,:) = RESHAPE(buffer(:,2),(/ nproma, patch_2d%nblks_c /) )
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%stress_yw(:,:))
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%stress_y (:,:))
    ENDIF
    !
    ! Apply freshwater flux - 2 parts, precipitation and evaporation - record 3
    !  - here freshwater can be bracketed by forcing_enable_freshwater, i.e. it must not be passed through coupler if not used
    ! IF (forcing_enable_freshwater) THEN
    field_shape(3) = 3
#ifdef YAC_coupling
    CALL yac_fget ( field_id(3), nbr_hor_points, 3, 1, 1, buffer, info, ierror )
#else
    CALL icon_cpl_get ( field_id(3), field_shape, buffer(1:nbr_hor_points,1:3), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:3) = 0.0_wp
      surface_fluxes%FrshFlux_Precipitation(:,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )
      surface_fluxes%FrshFlux_SnowFall  (:,:) = RESHAPE(buffer(:,2),(/ nproma, patch_2d%nblks_c /) )
      surface_fluxes%FrshFlux_Evaporation  (:,:) = RESHAPE(buffer(:,3),(/ nproma, patch_2d%nblks_c /) )
      
      surface_fluxes%FrshFlux_Precipitation(:,:) = surface_fluxes%FrshFlux_Precipitation(:,:)/rho_ref
      surface_fluxes%FrshFlux_SnowFall  (:,:) = surface_fluxes%FrshFlux_SnowFall(:,:)/rho_ref
      surface_fluxes%FrshFlux_Evaporation  (:,:) = surface_fluxes%FrshFlux_Evaporation(:,:)/rho_ref
      
      CALL sync_patch_array(sync_c, patch_2d, surface_fluxes%FrshFlux_Precipitation(:,:))
      CALL sync_patch_array(sync_c, patch_2d, surface_fluxes%FrshFlux_SnowFall(:,:))
      CALL sync_patch_array(sync_c, patch_2d, surface_fluxes%FrshFlux_Evaporation(:,:))
    END IF
    ! ENDIF ! forcing_enable_freshwater
    !
    ! Apply surface air temperature
    !  - it can be used for relaxing SST to T_a with type_surfRelax_Temp=1
    !  - set to 0 to omit relaxation to T_a=data_surfRelax_Temp(:,:)
    ! IF (type_surfRelax_Temp >=1) THEN
    field_shape(3) = 1
#ifdef YAC_coupling
    CALL yac_fget ( field_id(4), nbr_hor_points, 1, 1, 1, buffer, info, ierror )
#else
    CALL icon_cpl_get ( field_id(4), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
      surface_fluxes%data_surfRelax_Temp(:,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )
      !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
      surface_fluxes%data_surfRelax_Temp(:,:) = surface_fluxes%data_surfRelax_Temp(:,:) - tmelt
    END IF
    ! ENDIF  ! type_surfRelax_Temp >=1
    !
    ! Apply total heat flux - 4 parts - record 5
    ! surface_fluxes%swflx(:,:)  ocean short wave heat flux                              [W/m2]
    ! surface_fluxes%lwflx(:,:)  ocean long  wave heat fluxe                             [W/m2]
    ! surface_fluxes%ssflx(:,:)  ocean sensible heat fluxes                              [W/m2]
    ! surface_fluxes%slflx(:,:)  ocean latent heat fluxes                                [W/m2]
    field_shape(3) = 4
#ifdef YAC_coupling
    CALL yac_fget ( field_id(5), nbr_hor_points, 4, 1, 1, buffer, info, ierror )
#else
    CALL icon_cpl_get ( field_id(5), field_shape, buffer(1:nbr_hor_points,1:4), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:4) = 0.0_wp
      surface_fluxes%HeatFlux_ShortWave(:,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )
      surface_fluxes%HeatFlux_LongWave (:,:) = RESHAPE(buffer(:,2),(/ nproma, patch_2d%nblks_c /) )
      surface_fluxes%HeatFlux_Sensible (:,:) = RESHAPE(buffer(:,3),(/ nproma, patch_2d%nblks_c /) )
      surface_fluxes%HeatFlux_Latent   (:,:) = RESHAPE(buffer(:,4),(/ nproma, patch_2d%nblks_c /) )
      CALL sync_patch_array(sync_c, patch_2d, surface_fluxes%HeatFlux_ShortWave(:,:))
      CALL sync_patch_array(sync_c, patch_2d, surface_fluxes%HeatFlux_LongWave (:,:))
      CALL sync_patch_array(sync_c, patch_2d, surface_fluxes%HeatFlux_Sensible (:,:))
      CALL sync_patch_array(sync_c, patch_2d, surface_fluxes%HeatFlux_Latent   (:,:))
      ! sum of fluxes for ocean boundary condition
      surface_fluxes%HeatFlux_Total(:,:) = surface_fluxes%HeatFlux_ShortWave(:,:) + surface_fluxes%HeatFlux_LongWave(:,:) &
          & + surface_fluxes%HeatFlux_Sensible(:,:) + surface_fluxes%HeatFlux_Latent(:,:)
    ENDIF
    ! ice%Qtop(:,:)         Surface melt potential of ice                           [W/m2]
    ! ice%Qbot(:,:)         Bottom melt potential of ice                            [W/m2]
    ! ice%T1  (:,:)         Temperature of the upper ice layer                      [degC]
    ! ice%T2  (:,:)         Temperature of the lower ice layer                      [degC]
    field_shape(3) = 4
#ifdef YAC_coupling
    CALL yac_fget ( field_id(6), nbr_hor_points, 4, 1, 1, buffer, info, ierror )
#else
    CALL icon_cpl_get ( field_id(6), field_shape, buffer(1:nbr_hor_points,1:4), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:4) = 0.0_wp
      ice%qtop(:,1,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )
      ice%qbot(:,1,:) = RESHAPE(buffer(:,2),(/ nproma, patch_2d%nblks_c /) )
      ice%t1  (:,1,:) = RESHAPE(buffer(:,3),(/ nproma, patch_2d%nblks_c /) )
      ice%t2  (:,1,:) = RESHAPE(buffer(:,4),(/ nproma, patch_2d%nblks_c /) )
      CALL sync_patch_array(sync_c, patch_2d, ice%qtop(:,1,:))
      CALL sync_patch_array(sync_c, patch_2d, ice%qbot(:,1,:))
      CALL sync_patch_array(sync_c, patch_2d, ice%t1  (:,1,:))
      CALL sync_patch_array(sync_c, patch_2d, ice%t2  (:,1,:))
    END IF
    
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print(' CPL: Total  HF'     ,surface_fluxes%HeatFlux_Total    ,module_name,1,in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: SW-flux'       ,surface_fluxes%HeatFlux_ShortWave,module_name,2,in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: non-solar flux',surface_fluxes%HeatFlux_LongWave ,module_name,2,in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: Melt-pot. top' ,ice%qtop                         ,module_name,1,in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: Melt-pot. bot' ,ice%qbot                         ,module_name,1,in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: Precip.'       ,surface_fluxes%FrshFlux_Precipitation       ,module_name,1,in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: Evaporation'   ,surface_fluxes%FrshFlux_Evaporation         ,module_name,1,in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: Freshw. Flux'  ,surface_fluxes%FrshFlux_TotalSalt        ,module_name,1,in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------
    
    DEALLOCATE(buffer)
    DEALLOCATE(field_id)
    
    IF (ltimer) CALL timer_stop(timer_coupling)
    
  END SUBROUTINE couple_ocean_toatmo_fluxes
  !--------------------------------------------------------------------------
  
  
  !-------------------------------------------------------------------------
  !>
  !! Send data from atmosphere to ocean after initialization of ocean state
  !
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M, 2011-09
  !! LL: NOT Used
  !  SUBROUTINE init_coupled_ocean(patch_2d, ocean_state)
  !    TYPE(t_patch)                     :: patch_2d
  !    TYPE(t_hydro_ocean_state), TARGET :: ocean_state
  !
  !    sphere_radius = grid_sphere_radius
  !    u0 =(2.0_wp*pi*sphere_radius)/(12.0_wp*24.0_wp*3600.0_wp)
  
  !rr  ! Local declarations for coupling:
  !rr  INTEGER               :: ierror         !   return values form cpl_put/get calls
  !rr  INTEGER               :: nbr_hor_points ! = inner and halo points
  !rr  INTEGER               :: nbr_points     ! = nproma * nblks
  !rr  INTEGER               :: nbr_fields
  !rr  INTEGER, ALLOCATABLE  :: field_id(:)
  !rr  INTEGER               :: field_shape(3)
  !rr  REAL(wp), ALLOCATABLE :: buffer(:,:)
  !rr
  !rr  !-------------------------------------------------------------------------
  !rr
  !rr  IF ( is_coupled_run() ) THEN
  !rr
  !rr     nbr_hor_points = patch_2D%n_patch_cells
  !rr     nbr_points     = nproma * patch_2D%nblks_c
  !rr     ALLOCATE(buffer(nbr_points,1))
  !rr     !
  !rr     !  see drivers/mo_atmo_model.f90:
  !rr     !
  !rr     !   field_id(1) represents "TAUX"   wind stress component
  !rr     !   field_id(2) represents "TAUY"   wind stress component
  !rr     !   field_id(3) represents "SFWFLX" surface fresh water flux
  !rr     !   field_id(4) represents "SHFLX"  sensible heat flux
  !rr     !   field_id(5) represents "LHFLX"  latent heat flux
  !rr     !
  !rr     !   field_id(6) represents "SST"    sea surface temperature
  !rr     !   field_id(7) represents "OCEANU" u component of ocean surface current
  !rr     !   field_id(8) represents "OCEANV" v component of ocean surface current
  !rr     !
  !rr     CALL ICON_cpl_get_nbr_fields ( nbr_fields )
  !rr     ALLOCATE(field_id(nbr_fields))
  !rr     CALL ICON_cpl_get_field_ids ( nbr_fields, field_id )
  !rr     !
  !rr     field_shape(1) = 1
  !rr     field_shape(2) = patch_2D%n_patch_cells
  !rr     field_shape(3) = 1
  !rr
  !rr     !
  !rr     ! buffer is allocated over nproma only
  !rr
  !rr     !
  !rr     ! Send fields from ocean to atmosphere
  !rr     ! ------------------------------------
  !rr     !
  !rr     ! SST:
  !rr     buffer(:,1) = RESHAPE(ocean_state%p_prog(nold(1))%tracer(:,1,:,1), (/nbr_points /) ) + tmelt
  !rr     CALL ICON_cpl_put_init ( field_id(6), field_shape, buffer, ierror )
  !rr     !
  !rr     ! zonal velocity
  !rr     buffer(:,1) = RESHAPE(ocean_state%p_diag%u(:,1,:), (/nbr_points /) )
  !rr     CALL ICON_cpl_put_init ( field_id(7), field_shape, buffer, ierror )
  !rr     !
  !rr     ! meridional velocity
  !rr     buffer(:,1) = RESHAPE(ocean_state%p_diag%v(:,1,:), (/nbr_points /) )
  !rr     CALL ICON_cpl_put_init ( field_id(8), field_shape, buffer, ierror )
  !rr
  !rr     DeALLOCATE(field_id)
  !rr     DEALLOCATE(buffer)
  !rr
  !rr  END IF
  
  ! END SUBROUTINE init_coupled_ocean
  !-------------------------------------------------------------------------
  
#endif
  
END MODULE mo_ocean_coupling

