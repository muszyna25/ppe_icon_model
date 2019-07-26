!>
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

#ifdef YAC_coupling

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ocean_coupling

  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: warning, message
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_physical_constants,  ONLY: tmelt, rhoh2o
  USE mo_mpi,                 ONLY: p_pe_work
  USE mo_run_config,          ONLY: ltimer
  USE mo_dynamics_config,     ONLY: nold, nnew
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling, &
       &                            timer_coupling_put, timer_coupling_get,  &
       &                            timer_coupling_1stget, timer_coupling_init
  USE mo_sync,                ONLY: sync_c, sync_patch_array
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d

  USE mo_ocean_types
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_ocean_surface_types, ONLY: t_atmos_for_ocean
  USE mo_time_config,         ONLY: set_tc_current_date
  USE mtime,                  ONLY: datetime, datetimeToString, &
    &                               MAX_DATETIME_STR_LEN

  !-------------------------------------------------------------
  ! For the coupling
  !
  USE mo_math_constants,      ONLY: pi
  USE mo_parallel_config,     ONLY: nproma
  USE mo_yac_finterface,      ONLY: yac_finit, yac_fdef_comp, yac_fget_version,  &
    &                               yac_fdef_datetime,                           &
    &                               yac_fdef_subdomain, yac_fconnect_subdomains, &
    &                               yac_fdef_elements, yac_fdef_points,          &
    &                               yac_fdef_mask, yac_fdef_field, yac_fsearch,  &
    &                               yac_ffinalize, yac_fput, yac_fget,           &
    &                               YAC_LOCATION_CELL, COUPLING, OUT_OF_BOUND
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_time_config,         ONLY: time_config 
  USE mo_hamocc_nml,          ONLY: l_cpl_co2

  !-------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_ocean_coupling, destruct_ocean_coupling
  PUBLIC :: couple_ocean_toatmo_fluxes

  INTEGER, PARAMETER    :: no_of_fields = 14
  INTEGER               :: field_id(no_of_fields)

  REAL(wp), ALLOCATABLE :: buffer(:,:)
  INTEGER, SAVE         :: nbr_inner_cells
  LOGICAL, SAVE         :: lyac_very_1st_get

  CHARACTER(len=12)     :: str_module    = 'oceanCouplng'  ! Output of module for 1 line debug
  INTEGER               :: idt_src       = 1               ! Level of detail for 1 line debug

CONTAINS

  !--------------------------------------------------------------------------
  ! Prepare the coupling
  !
  ! For the time being this could all go into a subroutine which is
  ! common to atmo and ocean. Does this make sense if the setup deviates
  ! too much in future.
  !------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE construct_ocean_coupling(patch_3d)
    TYPE(t_patch_3d ), TARGET, INTENT(in)    :: patch_3d

    CHARACTER(LEN=max_char_length) ::  field_name(no_of_fields)
    INTEGER :: error_status

    INTEGER                :: patch_no
    TYPE(t_patch), POINTER :: patch_horz

    INTEGER, PARAMETER :: nbr_subdomain_ids = 1

    REAL(wp), PARAMETER :: deg = 180.0_wp / pi

    CHARACTER(LEN=max_char_length) :: xml_filename
    CHARACTER(LEN=max_char_length) :: xsd_filename
    CHARACTER(LEN=max_char_length) :: grid_name
    CHARACTER(LEN=max_char_length) :: comp_name

    INTEGER :: comp_id
    INTEGER :: comp_ids(1)
    INTEGER :: cell_point_ids(1)
    INTEGER :: cell_mask_ids(2)
    INTEGER :: domain_id
    INTEGER :: subdomain_id
    INTEGER :: subdomain_ids(nbr_subdomain_ids)
    INTEGER :: nbr_vertices_per_cell

    INTEGER :: mask_checksum
    INTEGER :: nblks
    INTEGER :: BLOCK, idx, INDEX

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER, ALLOCATABLE  :: buffer_c(:,:)
    INTEGER, ALLOCATABLE  :: ibuffer(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring

    ! Skip time measurement of the very first yac_fget
    ! as this will measure mainly the wait time caused
    ! by the initialisation of the model components
    ! and does not tell us much about the load balancing
    ! in subsequent calls.

    lyac_very_1st_get = .TRUE.

    IF (.NOT. is_coupled_run()) RETURN

    IF (ltimer) CALL timer_start(timer_coupling_init)

    comp_name = TRIM(get_my_process_name())

    patch_no = 1
    patch_horz => patch_3d%p_patch_2d(patch_no)

    ! Initialise the coupler
    xml_filename = "coupling.xml"
    xsd_filename = "coupling.xsd"
    CALL yac_finit ( TRIM(xml_filename), TRIM(xsd_filename) )

    ! Inform the coupler about what we are
    CALL yac_fdef_comp ( TRIM(comp_name), comp_id )
    comp_ids(1) = comp_id

    ! Print the YAC version
    CALL message('Running ICON ocean in coupled mode with YAC version ', TRIM(yac_fget_version()) )

    ! Overwrite job start and end date with component data
    CALL datetimeToString(time_config%tc_startdate, startdatestring)
    CALL datetimeToString(time_config%tc_stopdate, stopdatestring)

    CALL yac_fdef_datetime ( start_datetime = TRIM(startdatestring), &
         &                   end_datetime   = TRIM(stopdatestring)   )

    ! Announce one subdomain (patch) to the coupler
    grid_name = "grid1"
    CALL yac_fdef_subdomain ( comp_id, TRIM(grid_name), subdomain_id )

    subdomain_ids(1) = subdomain_id

    ! Extract cell information
    !
    ! cartesian coordinates of cell vertices are stored in
    ! patch_horz%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

    nblks = max(patch_horz%nblks_c,patch_horz%nblks_v)

    ALLOCATE(buffer_lon(nproma*nblks))
    ALLOCATE(buffer_lat(nproma*nblks))
    ALLOCATE(buffer_c(3,nproma*nblks))

    nbr_vertices_per_cell = 3

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_DEFAULT_SCHEDULE
    DO BLOCK = 1, patch_horz%nblks_v
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer_lon(INDEX) = patch_horz%verts%vertex(idx,BLOCK)%lon * deg
        buffer_lat(INDEX) = patch_horz%verts%vertex(idx,BLOCK)%lat * deg
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT

!ICON_OMP_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_DEFAULT_SCHEDULE
    DO BLOCK = 1, patch_horz%nblks_c
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer_c(1,INDEX) = (patch_horz%cells%vertex_blk(idx,BLOCK,1)-1)*nproma + &
          &                  patch_horz%cells%vertex_idx(idx,BLOCK,1)
        buffer_c(2,INDEX) = (patch_horz%cells%vertex_blk(idx,BLOCK,2)-1)*nproma + &
          &                  patch_horz%cells%vertex_idx(idx,BLOCK,2)
        buffer_c(3,INDEX) = (patch_horz%cells%vertex_blk(idx,BLOCK,3)-1)*nproma + &
                             patch_horz%cells%vertex_idx(idx,BLOCK,3)
      ENDDO
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    ! Description of elements, here as unstructured grid
    CALL yac_fdef_elements (      &
      & subdomain_id,             &
      & patch_horz%n_patch_verts, &
      & patch_horz%n_patch_cells, &
      & nbr_vertices_per_cell,    &
      & buffer_lon,               &
      & buffer_lat,               &
      & buffer_c )

    ! Can we have two fdef_point calls for the same subdomain, i.e.
    ! one single set of cells?
    !
    ! Define cell center points (location = 0)
    !
    ! cartesian coordinates of cell centers are stored in
    ! patch_horz%cells%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

!ICON_OMP_PARALLEL_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_DEFAULT_SCHEDULE
    DO BLOCK = 1, patch_horz%nblks_c
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer_lon(INDEX) = patch_horz%cells%center(idx,BLOCK)%lon * deg
        buffer_lat(INDEX) = patch_horz%cells%center(idx,BLOCK)%lat * deg
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    ! center points in cells (needed e.g. for patch recovery and nearest neighbour)
    CALL yac_fdef_points (        &
      & subdomain_id,             &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & buffer_lon,               &
      & buffer_lat,               &
      & cell_point_ids(1) )

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c)

    ALLOCATE(ibuffer(nproma*patch_horz%nblks_c))

    nbr_inner_cells = 0
!ICON_OMP_PARALLEL DO PRIVATE(idx) REDUCTION(+:nbr_inner_cells) ICON_OMP_DEFAULT_SCHEDULE
    DO idx = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(idx) ) THEN
         ibuffer(idx) = -1
         nbr_inner_cells = nbr_inner_cells + 1
       ELSE
         ibuffer(idx) = patch_horz%cells%decomp_info%owner_local(idx)
       ENDIF
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    ! decomposition information
    CALL yac_fdef_index_location (              &
      & subdomain_id,                           &
      & patch_horz%n_patch_cells,               &
      & YAC_LOCATION_CELL,                      &
      & patch_horz%cells%decomp_info%glb_index, &
      & ibuffer )

    ! Connect subdomains
    CALL yac_fconnect_subdomains ( &
      & comp_id,                   &
      & nbr_subdomain_ids,         &
      & subdomain_ids,             &
      & domain_id )

    !
    ! mask generation : ... not yet defined ...
    !
    ! We could use the patch_horz%cells%decomp_info%owner_local information
    ! e.g. to mask out halo points. We do we get the info about what is local and what
    ! is remote.
    !
    ! The integer land-sea mask:
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! This integer mask for the ocean is available in patch_3D%surface_cell_sea_land_mask(:,:)
    ! The logical mask for the coupler is set to .FALSE. for land points to exclude them from mapping by yac.
    ! These points are not touched by yac.

    mask_checksum = 0
!ICON_OMP_PARALLEL_DO PRIVATE(BLOCK,idx) REDUCTION(+:mask_checksum) ICON_OMP_DEFAULT_SCHEDULE
    DO BLOCK = 1, patch_horz%nblks_c
      DO idx = 1, nproma
        mask_checksum = mask_checksum + ABS(patch_3d%surface_cell_sea_land_mask(idx, BLOCK))
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF ( mask_checksum > 0 ) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_DEFAULT_SCHEDULE
      DO BLOCK = 1, patch_horz%nblks_c
        DO idx = 1, nproma
          IF ( patch_3d%surface_cell_sea_land_mask(idx, BLOCK) < 0 ) THEN
            ! ocean and ocean-coast is valid (-2, -1)
            ibuffer((BLOCK-1)*nproma+idx) = 0
          ELSE
            ! land is undef (1, 2)
            ibuffer((BLOCK-1)*nproma+idx) = 1
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ELSE

!ICON_OMP_PARALLEL_DO PRIVATE(idx) ICON_OMP_DEFAULT_SCHEDULE
      DO idx = 1, patch_horz%nblks_c * nproma
        ibuffer(idx) = 0
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ENDIF

    CALL yac_fdef_mask (          &
      & patch_horz%n_patch_cells, &
      & ibuffer,                  &
      & cell_point_ids(1),        &
      & cell_mask_ids(1) )

    field_name(1) = "surface_downward_eastward_stress"   ! bundled field containing two components
    field_name(2) = "surface_downward_northward_stress"  ! bundled field containing two components
    field_name(3) = "surface_fresh_water_flux"           ! bundled field containing three components
    field_name(4) = "total_heat_flux"                    ! bundled field containing four components
    field_name(5) = "atmosphere_sea_ice_bundle"          ! bundled field containing two components
    field_name(6) = "sea_surface_temperature"
    field_name(7) = "eastward_sea_water_velocity"
    field_name(8) = "northward_sea_water_velocity"
    field_name(9) = "ocean_sea_ice_bundle"               ! bundled field containing three components
    field_name(10) = "10m_wind_speed"
    field_name(11) = "river_runoff"
    field_name(12) = "co2_mixing_ratio"
    field_name(13) = "co2_flux"
    field_name(14) = "sea_level_pressure"

    ! Define the mask for all fields but the runoff

    DO idx = 1, no_of_fields 
      if(field_name(idx).ne. "river_runoff")then
      CALL yac_fdef_field (      &
        & TRIM(field_name(idx)), &
        & comp_id,               &
        & domain_id,             &
        & cell_point_ids,        &
        & cell_mask_ids(1),      &
        & 1,                     &
        & field_id(idx) )
      endif
    ENDDO

    ! Define cell_mask_ids(2) for runoff: all ocean points are valid.
    !!slo! Define cell_mask_ids(2) for runoff: ocean coastal points only are valid.
    !!slo!  - todo: use same mask as for other ones: all points, better wet points only

    IF ( mask_checksum > 0 ) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_DEFAULT_SCHEDULE
      DO BLOCK = 1, patch_horz%nblks_c
        DO idx = 1, nproma
          ! ocean coast (-1) is valid
!         IF ( patch_3d%surface_cell_sea_land_mask(idx, BLOCK) == -1 ) THEN
          ! all ocean points (-1, -2) are valid
          IF ( patch_3d%surface_cell_sea_land_mask(idx, BLOCK) <= -1 ) THEN
            ibuffer((BLOCK-1)*nproma+idx) = 0
          ELSE
            ! elsewhere (land or open ocean 1, 2, -2) is undef
            ibuffer((BLOCK-1)*nproma+idx) = 1
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ELSE

!ICON_OMP_PARALLEL_DO PRIVATE(idx) ICON_OMP_DEFAULT_SCHEDULE
      DO idx = 1, patch_horz%nblks_c * nproma
        ibuffer(idx) = 0
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ENDIF

    CALL yac_fdef_mask (          &
      & patch_horz%n_patch_cells, &
      & ibuffer,                  &
      & cell_point_ids(1),        &
      & cell_mask_ids(2) )


    DEALLOCATE(ibuffer)

    ! Define the mask for runoff

    ! Utilize mask field for runoff
    !  - cell_mask_ids(1) is whole ocean for nearest neighbor interpolation (pre03)
    !  - cell_mask_ids(2) is ocean coast points only for source point mapping (pre04, source_to_target_map)

    CALL yac_fdef_field (               &
      & TRIM("river_runoff"), &
      & comp_id,                        &
      & domain_id,                      &
      & cell_point_ids,                 &
      & cell_mask_ids(2),               &
      & 1,                              &
      & field_id(11) )

    CALL yac_fsearch ( 1, comp_ids, no_of_fields, field_id, error_status )

    ALLOCATE(buffer(nproma * patch_horz%nblks_c,5))

!ICON_OMP_PARALLEL_DO PRIVATE(BLOCK, INDEX, idx) ICON_OMP_DEFAULT_SCHEDULE
    DO BLOCK = 1, patch_horz%nblks_c
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer(INDEX,1) = 0.0_wp
        buffer(INDEX,2) = 0.0_wp
        buffer(INDEX,3) = 0.0_wp
        buffer(INDEX,4) = 0.0_wp
        buffer(INDEX,5) = 0.0_wp
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    IF (ltimer) CALL timer_stop(timer_coupling_init)

  END SUBROUTINE construct_ocean_coupling

  !--------------------------------------------------------------------------

!<Optimize:inUse>
  SUBROUTINE destruct_ocean_coupling()

    IF (.NOT. is_coupled_run()) RETURN

    DEALLOCATE(buffer)

    CALL yac_ffinalize

  END SUBROUTINE destruct_ocean_coupling

  !--------------------------------------------------------------------------

  SUBROUTINE couple_ocean_toatmo_fluxes(patch_3d, ocean_state, ice, atmos_fluxes, atmos_forcing, this_datetime)

    TYPE(t_patch_3d ),TARGET, INTENT(in)        :: patch_3d
    TYPE(t_hydro_ocean_state)                   :: ocean_state
    TYPE(t_sea_ice)                             :: ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_atmos_for_ocean)                     :: atmos_forcing
    TYPE(datetime), POINTER                     :: this_datetime

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = 'couple_ocean_toatmo_fluxes'
    
    ! Local declarations for coupling:
    LOGICAL :: write_coupler_restart
    INTEGER :: nbr_hor_cells  ! = inner and halo points
    INTEGER :: n              ! nproma loop count
    INTEGER :: nn             ! block offset
    INTEGER :: i_blk          ! block loop count
    INTEGER :: nlen           ! nproma/npromz
    INTEGER :: no_arr         !  no of arrays in bundle for put/get calls
    TYPE(t_patch), POINTER:: patch_horz

    INTEGER                             :: info, ierror   !< return values from cpl_put/get calls
    REAL(wp), PARAMETER                 :: dummy = 0.0_wp
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: datestring
    REAL(wp)                            :: total_rain

    IF (.NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start(timer_coupling)

    patch_horz   => patch_3D%p_patch_2D(1)

!    CALL datetimeToString(this_datetime, datestring)
!    CALL set_tc_current_date(TRIM(datestring))

    nbr_hor_cells = patch_horz%n_patch_cells

    !
    !
    !  Receive fields from atmosphere
    !   field_id(1) represents "surface_downward_eastward_stress" bundle  - zonal wind stress component over ice and water
    !   field_id(2) represents "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !   field_id(3) represents "surface_fresh_water_flux" bundle          - liquid rain, snowfall, evaporation
    !   field_id(4) represents "total heat flux" bundle                   - short wave, long wave, sensible, latent heat flux
    !   field_id(5) represents "atmosphere_sea_ice_bundle"                - sea ice surface and bottom melt potentials
    !   field_id(10) represents "10m_wind_speed"                          - atmospheric wind speed
    !   field_id(12) represents "co2_mixing_ratio"                        - co2 mixing ratio
    !
    !  Receive field from HD-model:
    !   field_id(11) represents "river_runoff"                            - river discharge into the ocean
    !
    !  Send fields to atmosphere:
    !   field_id(6) represents "sea_surface_temperature"                  - SST
    !   field_id(7) represents "eastward_sea_water_velocity"              - zonal velocity, u component of ocean surface current
    !   field_id(8) represents "northward_sea_water_velocity"             - meridional velocity, v component of ocean surface current
    !   field_id(9) represents "ocean_sea_ice_bundle"                     - ice thickness, snow thickness, ice concentration
    !

    !  add another comment here - 


    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from ocean to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    write_coupler_restart = .FALSE.

    !
    ! ------------------------------
    !  Send SST
    !   field_id(6) represents "sea_surface_temperature" - SST
    !
    buffer(:,:) = 0.0_wp  ! temporarily
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO i_blk = 1, patch_horz%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO n = 1, nlen
        ! use SST updated by ocean dynamics - coupling is done at end of timestep
        buffer(nn+n,1) = ocean_state%p_prog(nnew(1))%tracer(n,1,i_blk,1) + tmelt
        !buffer(nn+n,1) = ocean_state%p_prog(nold(1))%tracer(n,1,i_blk,1) + tmelt
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !    
    IF (ltimer) CALL timer_start(timer_coupling_put)

    CALL yac_fput ( field_id(6), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fput called after end of run - id=6, SST')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send zonal velocity
    !   field_id(7) represents "eastward_sea_water_velocity" - zonal velocity, u component of ocean surface current
    !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO i_blk = 1, patch_horz%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO n = 1, nlen
        buffer(nn+n,1) = ocean_state%p_diag%u(n,1,i_blk)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    CALL yac_fput ( field_id(7), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fput called after end of run - id=7, u velocity')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send meridional velocity
    !   field_id(8) represents "northward_sea_water_velocity" - meridional velocity, v component of ocean surface current
    !
    buffer(:,:) = 0.0_wp  ! temporarily
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO i_blk = 1, patch_horz%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO n = 1, nlen
        buffer(nn+n,1) = ocean_state%p_diag%v(n,1,i_blk)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    CALL yac_fput ( field_id(8), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fput called after end of run - id=8, v velocity')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send sea ice bundle
    !   field_id(9) represents "ocean_sea_ice_bundle" - ice thickness, snow thickness, ice concentration
    !
    buffer(:,:) = 0.0_wp  ! temporarily
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO i_blk = 1, patch_horz%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO n = 1, nlen
        buffer(nn+n,1) = ice%hi  (n,1,i_blk)
        buffer(nn+n,2) = ice%hs  (n,1,i_blk)
        buffer(nn+n,3) = ice%conc(n,1,i_blk)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 3
    CALL yac_fput ( field_id(9), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) write_coupler_restart = .TRUE.
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fput called after end of run - id=8, sea ice bundle')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( write_coupler_restart ) THEN
       CALL message('couple_ocean_toatmo_fluxes', 'YAC says it is put for restart - ids = 6 to 9, ocean fields')
    ENDIF


    IF (l_cpl_co2) THEN
      !
      ! ------------------------------
      !  Send CO2 flux
      !   field_id(13) represents "co2_flux" - co2flux
      !
      buffer(:,:) = 0.0_wp  ! temporarily
      !ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          buffer(nn+n,1) = atmos_forcing%co2flx(n,i_blk) 
        ENDDO
      ENDDO
      !ICON_OMP_END_PARALLEL_DO
      !    
      IF (ltimer) CALL timer_start(timer_coupling_put)

      CALL yac_fput ( field_id(13), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
      IF ( info > COUPLING .AND. info < OUT_OF_BOUND ) THEN
        CALL message('couple_ocean_toatmo_fluxes', 'YAC says it is put for restart - id=13, CO2 flux')
      ENDIF
      IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                     'YAC says fput called after end of run - id=13, CO2 flux')

      IF (ltimer) CALL timer_stop(timer_coupling_put)

    ENDIF

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    !

    !
    ! ------------------------------
    !  Receive zonal wind stress bundle
    !   field_id(1) represents "surface_downward_eastward_stress" bundle - zonal wind stress component over ice and water
    !
    IF ( .NOT. lyac_very_1st_get ) THEN
      IF (ltimer) CALL timer_start(timer_coupling_1stget)
    ENDIF

    buffer(:,:) = 0.0_wp  ! temporarily
    no_arr = 2
    CALL yac_fget ( field_id(1), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND )                       &
         &                      CALL message('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says it is get for restart - id=1, u-stress')
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fget called after end of run - id=1, u-stress')

    IF ( .NOT. lyac_very_1st_get ) THEN
       IF (ltimer) CALL timer_stop(timer_coupling_1stget)
    ENDIF

    lyac_very_1st_get = .FALSE.

    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            atmos_fluxes%stress_xw(n,i_blk) = dummy
            atmos_fluxes%stress_x (n,i_blk) = dummy
          ELSE
            atmos_fluxes%stress_xw(n,i_blk) = buffer(nn+n,1)
            atmos_fluxes%stress_x (n,i_blk) = buffer(nn+n,2)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%stress_xw(:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%stress_x (:,:))
    ENDIF

    !
    ! ------------------------------
    !  Receive meridional wind stress bundle
    !   field_id(2) represents "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp  ! temporarily
    no_arr = 2
    CALL yac_fget ( field_id(2), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND )                       &
         &                      CALL message('couple_ocean_toatmo_fluxes', &
         &                                  'YAC says it is get for restart - id=2, v-stress')
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fget called after end of run - id=2, v-stress')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            atmos_fluxes%stress_yw(n,i_blk) = dummy
            atmos_fluxes%stress_y (n,i_blk) = dummy
          ELSE
            atmos_fluxes%stress_yw(n,i_blk) = buffer(nn+n,1)
            atmos_fluxes%stress_y (n,i_blk) = buffer(nn+n,2)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%stress_yw(:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%stress_y (:,:))
    ENDIF

    !
    ! ------------------------------
    !  Receive surface fresh water flux bundle
    !   field_id(3) represents "surface_fresh_water_flux" bundle - liquid rain, snowfall, evaporation

    ! Note: freshwater fluxes are received in kg/m^2/s and are converted to m/s by division by rhoh2o below.
    ! Note: precipitation is the sum of rain and snowfall
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp  ! temporarily
    no_arr = 3
    CALL yac_fget ( field_id(3), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND )                       &
         &                      CALL message('couple_ocean_toatmo_fluxes', &
         &                                      'YAC says it is get for restart - id=3, surface fresh water flux')
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fget called after end of run - id=3, surface fresh water flux')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen, total_rain) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            total_rain                                   = dummy
            atmos_fluxes%FrshFlux_SnowFall     (n,i_blk) = dummy
            atmos_fluxes%FrshFlux_Evaporation  (n,i_blk) = dummy
            atmos_fluxes%FrshFlux_Precipitation(n,i_blk) = dummy
          ELSE
            total_rain                                   = buffer(nn+n,1) / rhoh2o
            atmos_fluxes%FrshFlux_SnowFall     (n,i_blk) = buffer(nn+n,2) / rhoh2o
            atmos_fluxes%FrshFlux_Evaporation  (n,i_blk) = buffer(nn+n,3) / rhoh2o
            atmos_fluxes%FrshFlux_Precipitation(n,i_blk) = total_rain + &
              &  atmos_fluxes%FrshFlux_SnowFall(n,i_blk)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%FrshFlux_Precipitation(:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%FrshFlux_SnowFall     (:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%FrshFlux_Evaporation  (:,:))
    END IF

    !
    ! ------------------------------
    !  Receive total heat flux bundle
    !   field_id(4) represents "total heat flux" bundle - short wave, long wave, sensible, latent heat flux
    !
    ! atmos_fluxes%swflx(:,:)  ocean short wave heat flux                              [W/m2]
    ! atmos_fluxes%lwflx(:,:)  ocean long  wave heat fluxe                             [W/m2]
    ! atmos_fluxes%ssflx(:,:)  ocean sensible heat fluxes                              [W/m2]
    ! atmos_fluxes%slflx(:,:)  ocean latent heat fluxes                                [W/m2]
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp  ! temporarily
    no_arr = 4
    CALL yac_fget ( field_id(4), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND )                       &
         &                      CALL message('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says it is get for restart - id=4, heat flux')
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fget called after end of run - id=4, heat flux')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            atmos_fluxes%HeatFlux_ShortWave(n,i_blk) = dummy
            atmos_fluxes%HeatFlux_LongWave (n,i_blk) = dummy
            atmos_fluxes%HeatFlux_Sensible (n,i_blk) = dummy
            atmos_fluxes%HeatFlux_Latent   (n,i_blk) = dummy
          ELSE
            atmos_fluxes%HeatFlux_ShortWave(n,i_blk) = buffer(nn+n,1)
            atmos_fluxes%HeatFlux_LongWave (n,i_blk) = buffer(nn+n,2)
            atmos_fluxes%HeatFlux_Sensible (n,i_blk) = buffer(nn+n,3)
            atmos_fluxes%HeatFlux_Latent   (n,i_blk) = buffer(nn+n,4)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%HeatFlux_ShortWave(:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%HeatFlux_LongWave (:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%HeatFlux_Sensible (:,:))
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%HeatFlux_Latent   (:,:))

      ! sum of fluxes for ocean boundary condition
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          atmos_fluxes%HeatFlux_Total(n,i_blk) = atmos_fluxes%HeatFlux_ShortWave(n,i_blk) &
        &                                      + atmos_fluxes%HeatFlux_LongWave (n,i_blk) &
        &                                      + atmos_fluxes%HeatFlux_Sensible (n,i_blk) &
        &                                      + atmos_fluxes%HeatFlux_Latent   (n,i_blk)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    ENDIF

    !
    ! ------------------------------
    !  Receive sea ice flux bundle
    !   field_id(5) represents "atmosphere_sea_ice_bundle" - sea ice surface and bottom melt potentials
    !
    ! ice%Qtop(:,:)         Surface melt potential of ice                           [W/m2]
    ! ice%Qbot(:,:)         Bottom melt potential of ice                            [W/m2]
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp  ! temporarily
    no_arr = 2
    CALL yac_fget ( field_id(5), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND )                       &
         &                      CALL message('couple_ocean_toatmo_fluxes', &
         &                                      'YAC says it is get for restart - id=5, sea ice')
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fget called after end of run - id=5, sea ice')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            ice%qtop(n,1,i_blk) = dummy
            ice%qbot(n,1,i_blk) = dummy
          ELSE
            ice%qtop(n,1,i_blk) = buffer(nn+n,1)
            ice%qbot(n,1,i_blk) = buffer(nn+n,2)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, ice%qtop(:,1,:))
      CALL sync_patch_array(sync_c, patch_horz, ice%qbot(:,1,:))
    END IF

    !
    ! ------------------------------
    !  Receive 10m wind speed
    !   field_id(10) represents "10m_wind_speed" - atmospheric wind speed
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp  ! temporarily
    no_arr = 1
    CALL yac_fget ( field_id(10), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND )                       &
         &                      CALL message('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says it is get for restart - id=10, wind speed')
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fget called after end of run - id=10, wind speed')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            atmos_forcing%fu10(n,i_blk) = dummy
          ELSE
            atmos_forcing%fu10(n,i_blk) = buffer(nn+n,1)
            IF ( atmos_forcing%fu10(n,i_blk) < 0.0_wp ) atmos_forcing%fu10(n,i_blk) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
!!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_forcing%fu10(:,:))
    END IF

    !  Receive slp
    !   field_id(14) represents atmospheric sea level pressure
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp  ! temporarily
    no_arr = 1
    CALL yac_fget ( field_id(14), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND )                       &
         &                      CALL message('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says it is get for restart - id=14, sea level pressure')
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fget called after end of run - id=14, sea level pressure')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            atmos_forcing%pao(n,i_blk) = dummy
          ELSE
            atmos_forcing%pao(n,i_blk) = buffer(nn+n,1)
          ENDIF
        ENDDO
      ENDDO
!!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_forcing%pao(:,:))
    END IF
    IF(l_cpl_co2)then
    !
    ! ------------------------------
    !  Receive co2 mixing ratio
    !   field_id(12) represents "co2 mixing ratio" - atmospheric co2 mixing
    !   ratio
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp  ! temporarily
    no_arr = 1
    CALL yac_fget ( field_id(12), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND )                       &
         &                      CALL message('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says it is get for restart - id=12, co2 mr')
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fget called after end of run - id=12, co2 mixing ratio')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            atmos_forcing%co2(n,i_blk) = dummy
          ELSE
            atmos_forcing%co2(n,i_blk) = buffer(nn+n,1)
            IF ( atmos_forcing%co2(n,i_blk) < 0.0_wp ) atmos_forcing%co2(n,i_blk) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
!!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_forcing%co2(:,:))
    END IF
    END IF !l_cpl_co2

    !
    ! ------------------------------
    !  Receive river runoff
    !   field_id(11) represents "river_runoff" - river discharge into the ocean
    !
    ! Note: river runoff fluxes are received in m^3/s and are converted to m/s by division by whole grid area
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp  ! temporarily
    CALL yac_fget ( field_id(11), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > COUPLING .AND. info < OUT_OF_BOUND )                       &
         &                      CALL message('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says it is get for restart - id=11, runoff')
    IF ( info == OUT_OF_BOUND ) CALL warning('couple_ocean_toatmo_fluxes', &
         &                                   'YAC says fget called after end of run - id=11, runoff')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF (info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO i_blk = 1, patch_horz%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= patch_horz%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = patch_horz%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            atmos_fluxes%FrshFlux_Runoff(n,i_blk) = dummy
          ELSE
    ! !!! Note: freshwater fluxes are received in kg/m^2/s and are converted to m/s by division by rhoh2o below.
    ! !!!   atmos_fluxes%FrshFlux_Runoff(n,i_blk) = buffer(nn+n,1) / rhoh2o
    ! discharge_ocean is in m3/s
            atmos_fluxes%FrshFlux_Runoff(n,i_blk) = buffer(nn+n,1) / patch_horz%cells%area(n,i_blk)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%FrshFlux_Runoff(:,:))
    END IF

    !---------DEBUG DIAGNOSTICS-------------------------------------------

    CALL dbg_print('toatmo: AtmFluxStress_x  ',atmos_fluxes%stress_x              ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: AtmFluxStress_xw ',atmos_fluxes%stress_xw             ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: AtmFluxStress_y  ',atmos_fluxes%stress_y              ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: AtmFluxStress_yw ',atmos_fluxes%stress_yw             ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxPrecip   ',atmos_fluxes%FrshFlux_Precipitation,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxEvapo    ',atmos_fluxes%FrshFlux_Evaporation  ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxSnowFall ',atmos_fluxes%FrshFlux_SnowFall     ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxTotal    ',atmos_fluxes%HeatFlux_Total        ,str_module,2,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxShortwave',atmos_fluxes%HeatFlux_ShortWave    ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxLongwave ',atmos_fluxes%HeatFlux_Longwave     ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxSensible ',atmos_fluxes%HeatFlux_Sensible     ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxLatent   ',atmos_fluxes%HeatFlux_Latent       ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: ice%Qtop         ',ice%qtop                           ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: ice%Qbot         ',ice%qbot                           ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxRunoff   ',atmos_fluxes%FrshFlux_Runoff       ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: 10m_wind_speed   ',atmos_forcing%fu10                 ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: tracer(1): SST   ', ocean_state%p_prog(nold(1))%tracer(:,1,:,1) + tmelt &
      &                                                                           ,str_module,2,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: p_diag%u         ',ocean_state%p_diag%u(:,1,:)        ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: p_diag%v         ',ocean_state%p_diag%v(:,1,:)        ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: ice%hi           ',ice%hi(:,1,:)                      ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: ice%hs           ',ice%hs(:,1,:)                      ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: ice%conc         ',ice%conc(:,1,:)                    ,str_module,4,in_subset=patch_horz%cells%owned)
    !---------------------------------------------------------------------

    IF (ltimer) CALL timer_stop(timer_coupling)

  END SUBROUTINE couple_ocean_toatmo_fluxes
  !--------------------------------------------------------------------------

END MODULE mo_ocean_coupling

#else

MODULE mo_ocean_coupling

  USE mo_model_domain,        ONLY: t_patch_3d
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_ocean_surface_types, ONLY: t_atmos_for_ocean
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_exception,           ONLY: finish
  USE mtime,                  ONLY: datetime

  PUBLIC :: construct_ocean_coupling, destruct_ocean_coupling
  PUBLIC :: couple_ocean_toatmo_fluxes

CONTAINS

  SUBROUTINE construct_ocean_coupling(patch_3d)

    TYPE(t_patch_3d ), TARGET, INTENT(in)    :: patch_3d

    IF ( is_coupled_run() ) THEN
       CALL finish('construct_ocean_coupling: unintentionally called. Check your source code and configure.')
    ELSE
       RETURN
    ENDIF

  END SUBROUTINE construct_ocean_coupling

  SUBROUTINE couple_ocean_toatmo_fluxes(patch_3d, ocean_state, ice, atmos_fluxes, atmos_forcing, this_datetime)

    TYPE(t_patch_3d ),TARGET, INTENT(in)        :: patch_3d
    TYPE(t_hydro_ocean_state)                   :: ocean_state
    TYPE(t_sea_ice)                             :: ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_atmos_for_ocean)                     :: atmos_forcing
    TYPE(datetime), POINTER                     :: this_datetime

    IF ( is_coupled_run() ) THEN
       CALL finish('couple_ocean_toatmo_fluxes: unintentionally called. Check your source code and configure.')
    ELSE
       RETURN
    ENDIF

  END SUBROUTINE couple_ocean_toatmo_fluxes

  SUBROUTINE destruct_ocean_coupling()

    IF ( is_coupled_run() ) THEN
       CALL finish('destruct_ocean_coupling: unintentionally called. Check your source code and configure.')
    ELSE
       RETURN
    ENDIF

  END SUBROUTINE destruct_ocean_coupling

END MODULE mo_ocean_coupling

#endif
