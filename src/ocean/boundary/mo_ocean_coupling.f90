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
  USE mo_exception,           ONLY: warning
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_physical_constants,  ONLY: tmelt, rhoh2o
  USE mo_mpi,                 ONLY: p_pe_work
  USE mo_datetime,            ONLY: t_datetime
  USE mo_time_config,         ONLY: time_config
  USE mo_run_config,          ONLY: ltimer
  USE mo_dynamics_config,     ONLY: nold
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling, &
       &                            timer_coupling_put, timer_coupling_get,  &
       &                            timer_coupling_1stget, timer_coupling_init
  USE mo_sync,                ONLY: sync_c, sync_patch_array
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d

  USE mo_ocean_types
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes

  !-------------------------------------------------------------
  ! For the coupling
  !
  USE mo_math_constants,      ONLY: pi
  USE mo_parallel_config,     ONLY: nproma
  USE mo_yac_finterface,      ONLY: yac_finit, yac_fdef_comp,                    &
    &                               yac_fdef_datetime,                           &
    &                               yac_fdef_subdomain, yac_fconnect_subdomains, &
    &                               yac_fdef_elements, yac_fdef_points,          &
    &                               yac_fdef_mask, yac_fdef_field, yac_fsearch,  &
    &                               yac_ffinalize, yac_fput, yac_fget
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_mtime_extensions,    ONLY: get_datetime_string
  USE mo_output_event_types,  ONLY: t_sim_step_info

  !-------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_ocean_coupling, destruct_ocean_coupling
  PUBLIC :: couple_ocean_toatmo_fluxes

! CHARACTER(LEN=12)     :: module_name    = 'ocean_coupli'

  INTEGER, PARAMETER    :: no_of_fields = 11
  INTEGER               :: field_id(no_of_fields)

  REAL(wp), ALLOCATABLE :: buffer(:,:)
  INTEGER, SAVE         :: nbr_inner_cells

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
    INTEGER, PARAMETER :: CELL = 0 ! one point per cell
    ! (see definition of enum location in points.h)

    REAL(wp), PARAMETER :: deg = 180.0_wp / pi

    CHARACTER(LEN=max_char_length) :: xml_filename
    CHARACTER(LEN=max_char_length) :: xsd_filename
    CHARACTER(LEN=max_char_length) :: grid_name
    CHARACTER(LEN=max_char_length) :: comp_name

    INTEGER :: comp_id
    INTEGER :: comp_ids(1)
    INTEGER :: cell_point_ids(1)
    INTEGER :: cell_mask_ids(1)
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

    TYPE(t_sim_step_info) :: sim_step_info

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

    ! Overwrite job start and end date with component data
    CALL get_datetime_string(sim_step_info%run_start,    time_config%cur_datetime)
    CALL get_datetime_string(sim_step_info%restart_time, time_config%cur_datetime, &
      & INT(time_config%dt_restart))

    CALL yac_fdef_datetime ( start_datetime = TRIM(sim_step_info%run_start), &
      &                      end_datetime   = TRIM(sim_step_info%restart_time)   )

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
      & CELL,                     &
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
      & CELL,                                   &
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
    ! The land-sea mask for the ocean is available in patch_3D%surface_cell_sea_land_mask(:,:)
    !
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !

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
            ! water (-2, -1)
            ibuffer((BLOCK-1)*nproma+idx) = 0
          ELSE
            ! land or boundary
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

    DEALLOCATE (ibuffer)

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

    DO idx = 1, no_of_fields 
      CALL yac_fdef_field (      &
        & TRIM(field_name(idx)), &
        & comp_id,               &
        & domain_id,             &
        & cell_point_ids,        &
        & cell_mask_ids,         &
        & 1,                     &
        & field_id(idx) )
    ENDDO

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

  SUBROUTINE couple_ocean_toatmo_fluxes(patch_3d, ocean_state, ice, atmos_fluxes, atm_wind_speed, datetime)

    TYPE(t_patch_3d ),TARGET, INTENT(in)        :: patch_3d
    TYPE(t_hydro_ocean_state)                   :: ocean_state
    TYPE(t_sea_ice)                             :: ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes !atmos_fluxes
    TYPE(t_datetime), INTENT(inout)             :: datetime
    ! currently, wind speed comes from t_atmos_for_ocean%fu10
    REAL(wp),         INTENT(inout)             :: atm_wind_speed(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)

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

    INTEGER :: info, ierror   !< return values from cpl_put/get calls

    REAL(wp) :: total_rain
    REAL(wp), PARAMETER :: dummy = 0.0_wp

    IF (.NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start(timer_coupling)

    patch_horz   => patch_3D%p_patch_2D(1)
    time_config%cur_datetime = datetime

    nbr_hor_cells = patch_horz%n_patch_cells

    !
    !  Receive fields from atmosphere
    !   field_id(1) represents "surface_downward_eastward_stress" bundle  - zonal wind stress component over ice and water
    !   field_id(2) represents "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !   field_id(3) represents "surface_fresh_water_flux" bundle          - liquid rain, snowfall, evaporation
    !   field_id(4) represents "total heat flux" bundle                   - short wave, long wave, sensible, latent heat flux
    !   field_id(5) represents "atmosphere_sea_ice_bundle"                - sea ice surface and bottom melt potentials
    !   field_id(10) represents "10m_wind_speed"                          - atmospheric wind speed
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
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO i_blk = 1, patch_horz%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= patch_horz%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = patch_horz%npromz_c
      END IF
      DO n = 1, nlen
        buffer(nn+n,1) = ocean_state%p_prog(nold(1))%tracer(n,1,i_blk,1) + tmelt
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !    
    IF (ltimer) CALL timer_start(timer_coupling_put)

    CALL yac_fput ( field_id(6), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fput called after end of run - id=6, SST')

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
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fput called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send meridional velocity
    !   field_id(8) represents "northward_sea_water_velocity" - meridional velocity, v component of ocean surface current
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
        buffer(nn+n,1) = ocean_state%p_diag%v(n,1,i_blk)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    CALL yac_fput ( field_id(8), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fput called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send sea ice bundle
    !   field_id(9) represents "ocean_sea_ice_bundle" - ice thickness, snow thickness, ice concentration
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
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fput called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( write_coupler_restart ) THEN
       CALL warning('couple_ocean_toatmo_fluxes', 'YAC says it is put for restart - id=9, ocean sea ice bundle')
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
    IF (ltimer) CALL timer_start(timer_coupling_1stget)

    no_arr = 2
    CALL yac_fget ( field_id(1), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says it is get for restart - id=1, u-stress')
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fget called after end of run - id=1, u-stress')

    IF (ltimer) CALL timer_stop(timer_coupling_1stget)
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

    no_arr = 2
    CALL yac_fget ( field_id(2), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says it is get for restart')
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fget called after end of run')

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

    no_arr = 3
    CALL yac_fget ( field_id(3), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says it is get for restart')
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fget called after end of run')

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

    no_arr = 4
    CALL yac_fget ( field_id(4), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says it is get for restart')
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fget called after end of run')

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

    no_arr = 2
    CALL yac_fget ( field_id(5), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says it is get for restart')
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fget called after end of run')

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

    no_arr = 1
    CALL yac_fget ( field_id(11), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says it is get for restart')
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fget called after end of run')

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
            atm_wind_speed(n,i_blk) = dummy
          ELSE
            atm_wind_speed(n,i_blk) = buffer(nn+n,1)
          ENDIF
        ENDDO
      ENDDO
!!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atm_wind_speed(:,:))
    END IF

    !
    ! ------------------------------
    !  Receive river runoff
    !   field_id(10) represents "river_runoff" - river discharge into the ocean
    !
    ! Note: freshwater fluxes are received in kg/m^2/s and are converted to m/s by division by rhoh2o below.
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    CALL yac_fget ( field_id(10), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > 1 .AND. info < 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says it is get for restart')
    IF ( info == 7 ) CALL warning('couple_ocean_toatmo_fluxes', 'YAC says fget called after end of run')

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
            atmos_fluxes%FrshFlux_Runoff(n,i_blk) = buffer(nn+n,1) / rhoh2o
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, patch_horz, atmos_fluxes%FrshFlux_Runoff(:,:))
    END IF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('toatmo: AtmFluxStress_x  ',atmos_fluxes%stress_x             ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: AtmFluxStress_xw ',atmos_fluxes%stress_xw            ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: AtmFluxStress_y  ',atmos_fluxes%stress_y             ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: AtmFluxStress_yw ',atmos_fluxes%stress_yw            ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxPrecip  ',atmos_fluxes%FrshFlux_Precipitation,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxEvapo    ',atmos_fluxes%FrshFlux_Evaporation ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxSnowFall ',atmos_fluxes%FrshFlux_SnowFall    ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxTotal    ',atmos_fluxes%HeatFlux_Total       ,str_module,2,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxShortwave',atmos_fluxes%HeatFlux_ShortWave   ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxLongwave ',atmos_fluxes%HeatFlux_Longwave    ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxSensible ',atmos_fluxes%HeatFlux_Sensible    ,str_module,4,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: HeatFluxLatent   ',atmos_fluxes%HeatFlux_Latent      ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: FrshFluxRunoff   ',atmos_fluxes%FrshFlux_Runoff      ,str_module,3,in_subset=patch_horz%cells%owned)
    CALL dbg_print('toatmo: 10m_wind_speed   ',atm_wind_speed                    ,str_module,3,in_subset=patch_horz%cells%owned)
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
  USE mo_datetime,            ONLY: t_datetime
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_exception,           ONLY: finish

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

  SUBROUTINE couple_ocean_toatmo_fluxes(patch_3d, ocean_state, ice, atmos_fluxes, datetime)

    TYPE(t_patch_3d ),TARGET, INTENT(in)        :: patch_3d
    TYPE(t_hydro_ocean_state)                   :: ocean_state
    TYPE(t_sea_ice)                             :: ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes !atmos_fluxes
    TYPE(t_datetime), INTENT(inout)             :: datetime

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
