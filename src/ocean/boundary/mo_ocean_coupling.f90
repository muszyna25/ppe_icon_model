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
MODULE mo_ocean_coupling

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: max_char_length
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

  USE mo_ocean_types
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean

  !-------------------------------------------------------------
  ! For the coupling
#ifndef __NO_ICON_ATMO__
# ifdef YAC_coupling
  USE mo_mpi,                 ONLY: p_n_work
  USE mo_math_constants,      ONLY: pi
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_yac_finterface,      ONLY: yac_finit, yac_fdef_comp,                    &
    &                               yac_fdef_datetime,                           &
    &                               yac_fdef_subdomain, yac_fconnect_subdomains, &
    &                               yac_fdef_elements, yac_fdef_points,          &
    &                               yac_fdef_mask, yac_fdef_field, yac_fsearch,  &
    &                               yac_ffinalize, yac_fput, yac_fget,           &
    &                               yac_fget_nbr_fields, yac_fget_field_ids,     &
    &                               yac_redirstdout
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_mtime_extensions,    ONLY: get_datetime_string
  USE mo_output_event_types,  ONLY: t_sim_step_info
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
!<Optimize:inUse>
  SUBROUTINE construct_ocean_coupling(patch_3d)
    TYPE(t_patch_3d ), TARGET, INTENT(in)    :: patch_3d
    RETURN
  END SUBROUTINE construct_ocean_coupling

!<Optimize:inUse>
  SUBROUTINE destruct_ocean_coupling()
    RETURN
  END SUBROUTINE destruct_ocean_coupling

  SUBROUTINE couple_ocean_toatmo_fluxes(patch_3d, ocean_state, ice, atmos_fluxes, jstep, datetime)
    TYPE(t_patch_3d ),TARGET, INTENT(in)        :: patch_3d
    TYPE(t_hydro_ocean_state)                   :: ocean_state
    TYPE(t_sea_ice)                             :: ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
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
!<Optimize:inUse>
  SUBROUTINE construct_ocean_coupling(patch_3d)
    TYPE(t_patch_3d ), TARGET, INTENT(in)    :: patch_3d

    INTEGER, PARAMETER :: no_of_fields = 10
    CHARACTER(LEN=max_char_length) ::  field_name(no_of_fields)
    INTEGER :: field_id(no_of_fields)
    INTEGER :: i, error_status

    INTEGER :: patch_no

# ifdef YAC_coupling

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
    INTEGER :: i_startidx, i_endidx
    INTEGER :: nblks
    INTEGER :: BLOCK, idx, INDEX

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER, ALLOCATABLE  :: buffer_c(:,:)
    INTEGER, ALLOCATABLE  :: ibuffer(:)

    TYPE(t_sim_step_info) :: sim_step_info
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_horz

    IF (.NOT. is_coupled_run()) RETURN

    comp_name = TRIM(get_my_process_name())

    patch_no = 1

    i = LEN_TRIM(comp_name)
    CALL yac_redirstdout ( TRIM(comp_name), i, 1, p_pe_work, p_n_work, error_status )

    ! Initialise the coupler
    xml_filename = "coupling.xml"
    xsd_filename = "coupling.xsd"
    CALL yac_finit ( xml_filename, xsd_filename )

    ! Inform the coupler about what we are
    CALL yac_fdef_comp ( comp_name, comp_id )
    comp_ids(1) = comp_id

    ! Overwrite job start and end date with component data
    CALL get_datetime_string(sim_step_info%restart_time,  time_config%cur_datetime, &
      &                      INT(time_config%dt_restart))
    CALL get_datetime_string(sim_step_info%run_start, time_config%cur_datetime)

    CALL yac_fdef_datetime ( start_datetime = TRIM(sim_step_info%sim_start), &
      &                      end_datetime   = TRIM(sim_step_info%sim_end)   )

    ! Announce one subdomain (patch) to the coupler
    grid_name = "grid1"
    CALL yac_fdef_subdomain ( comp_id, grid_name, subdomain_id )

    subdomain_ids(1) = subdomain_id

    patch_horz => patch_3d%p_patch_2d(patch_no)
    cells_in_domain  => patch_horz%cells%in_domain

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

    DO BLOCK = 1, patch_horz%nblks_v
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer_lon(INDEX) = patch_horz%verts%vertex(idx,BLOCK)%lon * deg
        buffer_lat(INDEX) = patch_horz%verts%vertex(idx,BLOCK)%lat * deg
      ENDDO
    ENDDO

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

    DO BLOCK = 1, patch_horz%nblks_c
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer_lon(INDEX) = patch_horz%cells%center(idx,BLOCK)%lon * deg
        buffer_lat(INDEX) = patch_horz%cells%center(idx,BLOCK)%lat * deg
      ENDDO
    ENDDO

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

    DO idx = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(idx) ) THEN
         ibuffer(idx) = -1
       ELSE
         ibuffer(idx) = patch_horz%cells%decomp_info%owner_local(idx)
       ENDIF
    ENDDO

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

!rr    DO BLOCK = 1, patch_horz%nblks_c
!rr       CALL get_indices_c ( patch_horz, BLOCK, 1, patch_horz%nblks_c,  &
!rr                               i_startidx, i_endidx, 2 )
!rr       DO idx = i_startidx, i_endidx
!rr          mask_checksum = mask_checksum + ABS(patch_3d%surface_cell_sea_land_mask(idx, BLOCK))
!rr       ENDDO
!rr    ENDDO

    DO i = 1, patch_horz%n_patch_cells
       ibuffer(i) = 0
    ENDDO

    IF ( mask_checksum > 0 ) THEN
       DO BLOCK = 1, patch_horz%nblks_c
          CALL get_indices_c ( patch_horz, BLOCK, 1, patch_horz%nblks_c,  &
                               i_startidx, i_endidx, 2 )
          DO idx = i_startidx, i_endidx
             IF ( patch_3d%surface_cell_sea_land_mask(idx, BLOCK) < 0 ) THEN
                ibuffer((BLOCK-1)*nproma+idx) = 0
             ELSE
                ibuffer((BLOCK-1)*nproma+idx) = 1
             ENDIF
          ENDDO
       ENDDO
    ELSE
       DO i = 1, patch_horz%n_patch_cells
          ibuffer(i) = 0
       ENDDO
    ENDIF

    CALL yac_fdef_mask (          &
      & patch_horz%n_patch_cells, &
      & ibuffer,                &
      & cell_point_ids(1),        &
      & cell_mask_ids(1) )

    DEALLOCATE (ibuffer)

# else

    INTEGER :: grid_id
    INTEGER :: grid_shape(2)
    INTEGER :: field_shape(3)

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
      & grid_id, error_status )                                                  ! output

    ! Marker for internal and halo points, a list which contains the
    ! rank where the native cells are located.
    CALL icon_cpl_def_location ( &
      & grid_id, grid_shape, patch_3d%p_patch_2d(patch_no)%cells%decomp_info%owner_local, & ! input
      & p_pe_work,  & ! this owner id
      & error_status )                                            ! output

# endif

    field_name(1) =  "TAUX"   ! bundled field containing two components
    field_name(2) =  "TAUY"   ! bundled field containing two components
    field_name(3) =  "SFWFLX" ! bundled field containing three components
    field_name(4) =  "SFTEMP"
    field_name(5) =  "THFLX"  ! bundled field containing four components
    field_name(6) =  "ICEATM" ! bundled field containing four components
    field_name(7) =  "SST"
    field_name(8) =  "OCEANU"
    field_name(9) =  "OCEANV"
    field_name(10) = "ICEOCE" ! bundled field containing four components

# ifdef YAC_coupling

    DO i = 1, no_of_fields
      CALL yac_fdef_field (    &
        & field_name(i),       &
        & comp_id,             &
        & domain_id,           &
        & cell_point_ids,      &
        & cell_mask_ids,       &
        & 1,                   &
        & field_id(i) )
    ENDDO

    CALL yac_fsearch ( 1, comp_ids, no_of_fields, field_id, error_status )

# else

    field_shape(1:2) = grid_shape(1:2)

    ! see equivalent ocean counterpart in drivers/mo_atmo_model.f90
    ! routine construct_atmo_coupler 

    DO i = 1, no_of_fields

      IF ( i == 1 .OR. i == 2 ) THEN
        field_shape(3) = 2
      ELSE IF ( i == 3 ) THEN
        field_shape(3) = 3
      ELSE IF ( i == 5 .OR. i == 6 ) THEN
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
!<Optimize:inUse>
  SUBROUTINE destruct_ocean_coupling()
# ifdef YAC_coupling
    IF ( is_coupled_run() ) CALL yac_ffinalize
# else
    IF ( is_coupled_run() ) CALL icon_cpl_finalize ()
# endif
  END SUBROUTINE destruct_ocean_coupling
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE couple_ocean_toatmo_fluxes(patch_3d, ocean_state, ice, atmos_fluxes, jstep, datetime)

    TYPE(t_patch_3d ),TARGET, INTENT(in)        :: patch_3d
    TYPE(t_hydro_ocean_state)                   :: ocean_state
    TYPE(t_sea_ice)                             :: ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes !atmos_fluxes
    INTEGER, INTENT(in)                         :: jstep
    TYPE(t_datetime), INTENT(inout)             :: datetime
    !
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = 'couple_ocean_toatmo_fluxes'
    INTEGER :: jmon, jdmon, jmon1, jmon2, ylen, yday
    INTEGER :: iniyear, curyear, offset
    INTEGER :: jc, jb, no_set
    INTEGER :: i_startidx_c, i_endidx_c
    REAL(wp) :: z_tmin, z_relax, rday1, rday2, dtm1, dsec, z_smax
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
    TYPE(t_patch), POINTER:: patch_2d
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
#ifdef YAC_coupling
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
    field_shape(2) = nbr_hor_points
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
    buffer(:,1) = RESHAPE(ocean_state%p_prog(nold(1))%tracer(:,1,:,1), (/nbr_points /) ) + tmelt
    
#ifdef YAC_coupling
    CALL yac_fput ( field_id(7), nbr_hor_points, 1, 1, 1, buffer(1:nbr_hor_points,1:1), info, ierror )
    IF ( info > 0 ) WRITE ( 6 , * ) "ocean CPL SST", minval(buffer(1:nbr_hor_points,1:1)), maxval(buffer(1:nbr_hor_points,1:1))
#else
    CALL icon_cpl_put ( field_id(7), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
    !
    ! zonal velocity
    buffer(:,1) = RESHAPE(ocean_state%p_diag%u(:,1,:), (/nbr_points /) )
#ifdef YAC_coupling
    CALL yac_fput ( field_id(8), nbr_hor_points, 1, 1, 1, buffer(1:nbr_hor_points,1:1), info, ierror )
    IF ( info > 0 ) WRITE ( 6 , * ) "ocean CPL U", minval(buffer(1:nbr_hor_points,1:1)), maxval(buffer(1:nbr_hor_points,1:1))
#else
    CALL icon_cpl_put ( field_id(8), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.
    !
    ! meridional velocity
    buffer(:,1) = RESHAPE(ocean_state%p_diag%v(:,1,:), (/nbr_points /) )
#ifdef YAC_coupling
    CALL yac_fput ( field_id(9), nbr_hor_points, 1, 1, 1, buffer(1:nbr_hor_points,1:1), info, ierror )
    IF ( info > 0 ) WRITE ( 6 , * ) "ocean CPL V", minval(buffer(1:nbr_hor_points,1:1)), maxval(buffer(1:nbr_hor_points,1:1))
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
    CALL yac_fput ( field_id(10), nbr_hor_points, 5, 1, 1, buffer(1:nbr_hor_points,1:5), info, ierror )
    IF ( info > 0 ) THEN
        WRITE ( 6 , * ) "ocean CPL ice 1", minval(buffer(1:nbr_hor_points,1:1)), maxval(buffer(1:nbr_hor_points,1:1))
        WRITE ( 6 , * ) "ocean CPL ice 2", minval(buffer(1:nbr_hor_points,2:2)), maxval(buffer(1:nbr_hor_points,2:2))
        WRITE ( 6 , * ) "ocean CPL ice 3", minval(buffer(1:nbr_hor_points,3:3)), maxval(buffer(1:nbr_hor_points,3:3))
        WRITE ( 6 , * ) "ocean CPL ice 4", minval(buffer(1:nbr_hor_points,4:4)), maxval(buffer(1:nbr_hor_points,4:4))
        WRITE ( 6 , * ) "ocean CPL ice 5", minval(buffer(1:nbr_hor_points,5:5)), maxval(buffer(1:nbr_hor_points,5:5))
    ENDIF
#else
    CALL icon_cpl_put ( field_id(10), field_shape, buffer(1:nbr_hor_points,1:5), info, ierror )
#endif
    IF ( info == 2 ) write_coupler_restart = .TRUE.

#ifdef YAC_coupling
    IF ( write_coupler_restart ) &
      & CALL finish ( 'couple_ocean_toatmo_fluxes', 'ERROR: restart not supported yet for YAC')
#else
    IF ( write_coupler_restart ) CALL icon_cpl_write_restart ( 4, field_id(7:10), ierror )
#endif
    !
    ! Receive fields from atmosphere
    ! ------------------------------

    !
    ! Apply wind stress - records 0 and 1 of field_id

    ! zonal wind stress
    field_shape(3) = 2
#ifdef YAC_coupling
    CALL yac_fget ( field_id(1), nbr_hor_points, 2, 1, 1, buffer(1:nbr_hor_points,1:2), info, ierror )
#else
    CALL icon_cpl_get ( field_id(1), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:field_shape(3)) = 0.0_wp
      atmos_fluxes%stress_xw(:,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )
      atmos_fluxes%stress_x (:,:) = RESHAPE(buffer(:,2),(/ nproma, patch_2d%nblks_c /) ) !TODO + 100.0_wp
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%stress_xw(:,:))
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%stress_x (:,:))
    ENDIF
    !
    ! meridional wind stress
#ifdef YAC_coupling
    CALL yac_fget ( field_id(2), nbr_hor_points, 2, 1, 1, buffer(1:nbr_hor_points,1:2), info, ierror )
#else
    CALL icon_cpl_get ( field_id(2), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:field_shape(3)) = 0.0_wp
      atmos_fluxes%stress_yw(:,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )
      atmos_fluxes%stress_y (:,:) = RESHAPE(buffer(:,2),(/ nproma, patch_2d%nblks_c /) )  !TODO+ 100.0_wp
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%stress_yw(:,:))
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%stress_y (:,:))
    ENDIF
    !
    ! Apply freshwater flux - 2 parts, precipitation and evaporation - record 3
    !  - here freshwater can be bracketed by forcing_enable_freshwater, i.e. it must not be passed through coupler if not used
    ! IF (forcing_enable_freshwater) THEN
    field_shape(3) = 3
#ifdef YAC_coupling
    CALL yac_fget ( field_id(3), nbr_hor_points, 3, 1, 1, buffer(1:nbr_hor_points,1:3), info, ierror )
#else
    CALL icon_cpl_get ( field_id(3), field_shape, buffer(1:nbr_hor_points,1:3), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:3) = 0.0_wp
      atmos_fluxes%FrshFlux_Precipitation(:,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )
      atmos_fluxes%FrshFlux_SnowFall  (:,:) = RESHAPE(buffer(:,2),(/ nproma, patch_2d%nblks_c /) )
      atmos_fluxes%FrshFlux_Evaporation  (:,:) = RESHAPE(buffer(:,3),(/ nproma, patch_2d%nblks_c /) )

      atmos_fluxes%FrshFlux_Precipitation(:,:) = atmos_fluxes%FrshFlux_Precipitation(:,:)/rho_ref
      atmos_fluxes%FrshFlux_SnowFall  (:,:) = atmos_fluxes%FrshFlux_SnowFall(:,:)/rho_ref
      atmos_fluxes%FrshFlux_Evaporation  (:,:) = atmos_fluxes%FrshFlux_Evaporation(:,:)/rho_ref

      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%FrshFlux_Precipitation(:,:))
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%FrshFlux_SnowFall(:,:))
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%FrshFlux_Evaporation(:,:))
    END IF
    ! ENDIF ! forcing_enable_freshwater
    !
    ! Apply surface air temperature
    !  - it can be used for relaxing SST to T_a with type_surfRelax_Temp=1
    !  - set to 0 to omit relaxation to T_a=data_surfRelax_Temp(:,:)
    ! IF (type_surfRelax_Temp >=1) THEN
    field_shape(3) = 1
#ifdef YAC_coupling
    CALL yac_fget ( field_id(4), nbr_hor_points, 1, 1, 1, buffer(1:nbr_hor_points,1:1), info, ierror )
#else
    CALL icon_cpl_get ( field_id(4), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
      atmos_fluxes%data_surfRelax_Temp(:,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )
      !  - change units to deg C, subtract tmelt (0 deg C, 273.15)
      atmos_fluxes%data_surfRelax_Temp(:,:) = atmos_fluxes%data_surfRelax_Temp(:,:) - tmelt
    END IF
    ! ENDIF  ! type_surfRelax_Temp >=1
    !
    ! Apply total heat flux - 4 parts - record 5
    ! atmos_fluxes%swflx(:,:)  ocean short wave heat flux                              [W/m2]
    ! atmos_fluxes%lwflx(:,:)  ocean long  wave heat fluxe                             [W/m2]
    ! atmos_fluxes%ssflx(:,:)  ocean sensible heat fluxes                              [W/m2]
    ! atmos_fluxes%slflx(:,:)  ocean latent heat fluxes                                [W/m2]
    field_shape(3) = 4
#ifdef YAC_coupling
    CALL yac_fget ( field_id(5), nbr_hor_points, 4, 1, 1, buffer(1:nbr_hor_points,1:4), info, ierror )
#else
    CALL icon_cpl_get ( field_id(5), field_shape, buffer(1:nbr_hor_points,1:4), info, ierror )
#endif
    IF (info > 0 ) THEN
      buffer(nbr_hor_points+1:nbr_points,1:4) = 0.0_wp
      atmos_fluxes%HeatFlux_ShortWave(:,:) = RESHAPE(buffer(:,1),(/ nproma, patch_2d%nblks_c /) )  !TODO+ 300.0_wp
      atmos_fluxes%HeatFlux_LongWave (:,:) = RESHAPE(buffer(:,2),(/ nproma, patch_2d%nblks_c /) )  !TODO+ 300.0_wp
      atmos_fluxes%HeatFlux_Sensible (:,:) = RESHAPE(buffer(:,3),(/ nproma, patch_2d%nblks_c /) )  !TODO+ 300.0_wp
      atmos_fluxes%HeatFlux_Latent   (:,:) = RESHAPE(buffer(:,4),(/ nproma, patch_2d%nblks_c /) )  !TODO+ 300.0_wp
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%HeatFlux_ShortWave(:,:))
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%HeatFlux_LongWave (:,:))
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%HeatFlux_Sensible (:,:))
      CALL sync_patch_array(sync_c, patch_2d, atmos_fluxes%HeatFlux_Latent   (:,:))
      ! sum of fluxes for ocean boundary condition
      atmos_fluxes%HeatFlux_Total(:,:) = atmos_fluxes%HeatFlux_ShortWave(:,:) &
        &                              + atmos_fluxes%HeatFlux_LongWave(:,:) &
        &                              + atmos_fluxes%HeatFlux_Sensible(:,:) &
        &                              + atmos_fluxes%HeatFlux_Latent(:,:)
    ENDIF
    ! ice%Qtop(:,:)         Surface melt potential of ice                           [W/m2]
    ! ice%Qbot(:,:)         Bottom melt potential of ice                            [W/m2]
    ! ice%T1  (:,:)         Temperature of the upper ice layer                      [degC]
    ! ice%T2  (:,:)         Temperature of the lower ice layer                      [degC]
    field_shape(3) = 4
#ifdef YAC_coupling
    CALL yac_fget ( field_id(6), nbr_hor_points, 4, 1, 1, buffer(1:nbr_hor_points,1:4), info, ierror )
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
    CALL dbg_print(' CPL: Melt-pot. top' , ice%qtop                   , module_name, 1, in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: Melt-pot. bot' , ice%qbot                   , module_name, 1, in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: Total  HF'     , atmos_fluxes%HeatFlux_Total        , module_name, 1, in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: SW-flux'       , atmos_fluxes%HeatFlux_ShortWave    , module_name, 2, in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: non-solar flux', atmos_fluxes%HeatFlux_LongWave     , module_name, 2, in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: Precip.'       , atmos_fluxes%FrshFlux_Precipitation, module_name, 1, in_subset=patch_2d%cells%owned)
    CALL dbg_print(' CPL: Evaporation'   , atmos_fluxes%FrshFlux_Evaporation  , module_name, 1, in_subset=patch_2d%cells%owned)
    !---------------------------------------------------------------------

    DEALLOCATE(buffer)
    DEALLOCATE(field_id)

    IF (ltimer) CALL timer_stop(timer_coupling)

  END SUBROUTINE couple_ocean_toatmo_fluxes
  !--------------------------------------------------------------------------

#endif

END MODULE mo_ocean_coupling
