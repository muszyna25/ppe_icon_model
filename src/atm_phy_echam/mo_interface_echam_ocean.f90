!>
!! @brief Interface between ECHAM physics and the ocean, through a coupler
!!
!! @author Marco Giorgetta (MPI-M)
!!
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_interface_echam_ocean

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_echam_phy_memory    ,ONLY: prm_field
                                
  USE mo_parallel_config     ,ONLY: nproma
  
  USE mo_run_config          ,ONLY: nlev, ltimer
  USE mo_timer,               ONLY: timer_start, timer_stop,                &
       &                            timer_coupling_put, timer_coupling_get, &
       &                            timer_coupling_1stget, timer_coupling_init
  USE mo_echam_sfc_indices   ,ONLY: iwtr, iice

  USE mo_sync                ,ONLY: SYNC_C, sync_patch_array
  USE mo_impl_constants      ,ONLY: MAX_CHAR_LENGTH

  USE mo_ext_data_state      ,ONLY: ext_data
  USE mo_time_config         ,ONLY: time_config      ! variable


#ifdef YAC_coupling
  USE mo_master_control      ,ONLY: get_my_process_name

  USE mo_mpi                 ,ONLY: p_pe_work
  USE mo_math_constants      ,ONLY: pi
  USE mo_parallel_config     ,ONLY: nproma

  USE mo_coupling_config     ,ONLY: is_coupled_run
  USE mo_model_domain        ,ONLY: t_patch

  USE mo_exception           ,ONLY: warning
  USE mo_output_event_types  ,ONLY: t_sim_step_info
  USE mo_mtime_extensions    ,ONLY: get_datetime_string

  USE mo_yac_finterface      ,ONLY: yac_fput, yac_fget,                          &
    &                               yac_finit, yac_fdef_comp,                    &
    &                               yac_fdef_datetime,                           &
    &                               yac_fdef_subdomain, yac_fconnect_subdomains, &
    &                               yac_fdef_elements, yac_fdef_points,          &
    &                               yac_fdef_mask, yac_fdef_field, yac_fsearch,  &
    &                               yac_ffinalize

#else
  USE mo_master_control      ,ONLY: get_my_process_name, get_my_model_no

  USE mo_mpi                 ,ONLY: p_pe_work
  USE mo_icon_cpl_exchg      ,ONLY: ICON_cpl_put, ICON_cpl_get
  USE mo_icon_cpl_def_field  ,ONLY: ICON_cpl_get_nbr_fields, ICON_cpl_get_field_ids
  USE mo_icon_cpl_init       ,ONLY: icon_cpl_init
  USE mo_icon_cpl_init_comp  ,ONLY: icon_cpl_init_comp
  USE mo_icon_cpl_def_grid   ,ONLY: icon_cpl_def_grid, icon_cpl_def_location
  USE mo_icon_cpl_def_field  ,ONLY: icon_cpl_def_field
  USE mo_icon_cpl_search     ,ONLY: icon_cpl_search
  USE mo_icon_cpl_finalize   ,ONLY: icon_cpl_finalize
  USE mo_coupling_config     ,ONLY: is_coupled_run, config_debug_coupler_level

#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interface_echam_ocean
  PUBLIC :: construct_atmo_coupler, destruct_atmo_coupler

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_interface_echam_ocean'

  INTEGER, PARAMETER    :: no_of_fields = 10
  INTEGER               :: field_id(no_of_fields)

  REAL(wp), ALLOCATABLE :: buffer(:,:)
  INTEGER, SAVE         :: nbr_inner_cells

CONTAINS

  !>
  !! SUBROUTINE construct_atmo_coupler -- the initialisation for the coupling
  !! of ECHAM physics and the ocean, through a coupler

  SUBROUTINE construct_atmo_coupler (p_patch)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)

    CHARACTER(LEN=MAX_CHAR_LENGTH)    ::  field_name(no_of_fields)

    INTEGER :: i
    INTEGER :: error_status

    INTEGER                :: patch_no
    TYPE(t_patch), POINTER :: patch_horz

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------

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

    INTEGER :: mask_checksum
    INTEGER :: nblks
    INTEGER :: BLOCK, idx, INDEX
    INTEGER :: nbr_vertices_per_cell

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER,  ALLOCATABLE :: buffer_c(:,:)
    INTEGER,  ALLOCATABLE :: ibuffer(:)

    TYPE(t_sim_step_info)   :: sim_step_info  

    IF ( .NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start (timer_coupling_init)

    comp_name = TRIM(get_my_process_name())

    patch_no = 1
    patch_horz => p_patch(patch_no)

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

    CALL yac_fdef_datetime ( start_datetime = TRIM(sim_step_info%run_start),  &
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

    nbr_inner_cells = 0

    DO idx = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(idx) ) THEN
         ibuffer(idx) = -1
         nbr_inner_cells = nbr_inner_cells + 1
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
    ! The land-sea mask for the ocean is available in patch_3D%surface_cell_sea_land_mask(:,:)
    !
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! The mask which is used in the ECHAM physics is prm_field(1)%lsmask(:,:).
    ! This locial mask is set to .TRUE. for land points.
    ! We can get access to via "USE mo_echam_phy_memory,ONLY: prm_field"
    !
    ! Here we use a mask which is hopefully identical to the one used by the
    ! ocean, and which works independent of the physics chosen. 
    !
    mask_checksum = 0

    DO BLOCK = 1, patch_horz%nblks_c
      DO idx = 1, nproma
        mask_checksum = mask_checksum + ABS(ext_data(1)%atm%lsm_ctr_c(idx, BLOCK))
      ENDDO
    ENDDO

    IF ( mask_checksum > 0 ) THEN
       DO BLOCK = 1, patch_horz%nblks_c
          DO idx = 1, nproma
             IF ( ext_data(1)%atm%lsm_ctr_c(idx, BLOCK) < 0 ) THEN
               ! Ocean point
               ibuffer((BLOCK-1)*nproma+idx) = 0
             ELSE
               ! Land point
               ibuffer((BLOCK-1)*nproma+idx) = 1
             ENDIF
          ENDDO
       ENDDO
    ELSE
       DO i = 1, patch_horz%n_patch_cells
          ibuffer(i) = 0
       ENDDO
    ENDIF

    CALL yac_fdef_mask (           &
      & patch_horz%n_patch_cells,  &
      & ibuffer,                   &
      & cell_point_ids(1),         &
      & cell_mask_ids(1) )

    DEALLOCATE (ibuffer)

    field_name(1) = "surface_downward_eastward_stress"   ! bundled field containing two components
    field_name(2) = "surface_downward_northward_stress"  ! bundled field containing two components
    field_name(3) = "surface_fresh_water_flux"           ! bundled field containing three components
    field_name(4) = "surface_temperature"
    field_name(5) = "total_heat_flux"                    ! bundled field containing four components
    field_name(6) = "atmosphere_sea_ice_bundle"          ! bundled field containing four components
    field_name(7) = "sea_surface_temperature"
    field_name(8) = "eastward_sea_water_velocity"
    field_name(9) = "northward_sea_water_velocity"
    field_name(10) = "ocean_sea_ice_bundle"              ! bundled field containing five components

    DO i = 1, no_of_fields
      CALL yac_fdef_field (    &
        & TRIM(field_name(i)), &
        & comp_id,             &
        & domain_id,           &
        & cell_point_ids,      &
        & cell_mask_ids,       &
        & 1,                   &
        & field_id(i) )
    ENDDO

    CALL yac_fsearch ( 1, comp_ids, no_of_fields, field_id, error_status )

# else

    INTEGER :: grid_id
    INTEGER :: grid_shape(2)
    INTEGER :: field_shape(3)

    IF ( .NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start (timer_coupling_init)

    !------------------------------------------------------------
    CALL icon_cpl_init(debug_level=config_debug_coupler_level)
    ! Inform the coupler about what we are
    CALL icon_cpl_init_comp ( get_my_process_name(), get_my_model_no(), error_status )
    ! split the global_mpi_communicator into the components
    !------------------------------------------------------------

    patch_no = 1
    patch_horz => p_patch(patch_no)

    grid_shape(1)=1
    grid_shape(2)=patch_horz%n_patch_cells

    ! CALL get_patch_global_indexes ( patch_no, CELLS, no_of_entities, grid_glob_index )
    ! should grid_glob_index become a pointer in icon_cpl_def_grid as well?
    CALL icon_cpl_def_grid ( &
      & grid_shape, patch_horz%cells%decomp_info%glb_index, & ! input
      & grid_id, error_status )                               ! output

    ! Marker for internal and halo points, a list which contains the
    ! rank where the native cells are located.
    CALL icon_cpl_def_location ( &
      & grid_id, grid_shape, patch_horz%cells%decomp_info%owner_local, & ! input
      & p_pe_work,  &                                                    ! this owner id
      & error_status )                                                   ! output

    field_name(1) = "TAUX"   ! bundled field containing two components
    field_name(2) = "TAUY"   ! bundled field containing two components
    field_name(3) = "SFWFLX" ! bundled field containing three components
    field_name(4) = "SFTEMP"
    field_name(5) = "THFLX"  ! bundled field containing four components
    field_name(6) = "ICEATM" ! bundled field containing four components
    field_name(7) = "SST"
    field_name(8) = "OCEANU"
    field_name(9) = "OCEANV"
    field_name(10) = "ICEOCE" ! bundled field containing five components

    field_shape(1:2) = grid_shape(1:2)

    nbr_inner_cells = 0

    DO i = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(i) ) &
      &   nbr_inner_cells = nbr_inner_cells + 1
    ENDDO

    ! see equivalent atmosphere counterpart in ocean/boundary/mo_ocean_coupling.f90
    ! routine construct_ocean_coupling

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

       CALL icon_cpl_def_field ( &
         & field_name(i), grid_id, field_id(i), &
         & field_shape, error_status )

    ENDDO

    CALL icon_cpl_search

#endif

    ALLOCATE(buffer(nproma*patch_horz%nblks_c,5))
    buffer(:,:) = 0.0_wp

    IF (ltimer) CALL timer_stop(timer_coupling_init)

  END SUBROUTINE construct_atmo_coupler

  !>
  !! SUBROUTINE interface_icoham_echam -- the interface between
  !! ECHAM physics and the ocean, through a coupler
  !!
  !! This subroutine is called in the time loop of the ICONAM model.
  !! It takes the following as input:
  !! <ol>
  !! <li> prognostic and diagnostic variables of the dynamical core;
  !! <li> tendency of the prognostic varibles induced by adiabatic dynamics;
  !! <li> time step;
  !! <li> information about the dynamics grid;
  !! <li> interplation coefficients.
  !! </ol>
  !!
  !! The output includes tendencies of the prognostic variables caused by
  !! the parameterisations.
  !!
  !! Note that each call of this subroutine deals with a single grid level
  !! rather than the entire grid tree.

  SUBROUTINE interface_echam_ocean( jg      ,&! in
    &                               p_patch ) ! in

    ! Arguments

    INTEGER,               INTENT(IN)    :: jg            !< grid level/domain index
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: nbr_hor_cells  ! = inner and halo points
    INTEGER               :: nbr_cells      ! = nproma * nblks
    INTEGER               :: n              ! nproma loop count
    INTEGER               :: nn             ! block offset
    INTEGER               :: i_blk          ! block loop count
#ifndef YAC_coupling
    INTEGER               :: field_shape(3)
#endif
    INTEGER               :: info, ierror !< return values from cpl_put/get calls

    IF ( .NOT. is_coupled_run() ) RETURN

    !-------------------------------------------------------------------------
    ! If running in atm-oce coupled mode, exchange information 
    !-------------------------------------------------------------------------

    ! Possible fields that contain information to be sent to the ocean include
    !
    ! 1. prm_field(jg)% u_stress_tile(:,:,iwtr)  and 
    !    prm_field(jg)% v_stress_tile(:,:,iwtr)  which are the wind stress components;
    !
    ! 2. prm_field(jg)% evap_tile(:,:,iwtr) evaporation rate
    !
    ! 3. prm_field(jg)%rsfl + prm_field(jg)%rsfc + prm_field(jg)%ssfl + prm_field(jg)%ssfc
    !    which gives the precipitation rate;
    !
    ! 4. prm_field(jg)% temp(:,nlev,:)  temperature at the lowest model level, or
    !    prm_field(jg)% temp_2m(:,:)    2-m temperature, not available yet, or
    !    prm_field(jg)% shflx_tile(:,:,iwtr) sensible heat flux
    !
    ! 5  prm_field(jg)% lhflx_tile(:,:,iwtr) latent heat flux
    ! 6. shortwave radiation flux at the surface
    !
    ! Possible fields to receive from the ocean include
    !
    ! 1. prm_field(jg)% tsfc_tile(:,:,iwtr)   SST
    ! 2. prm_field(jg)% ocu(:,:) and ocv(:,:) ocean surface current
    ! 

    nbr_hor_cells = p_patch%n_patch_cells
    nbr_cells     = nproma * p_patch%nblks_c

    !
    !  see drivers/mo_atmo_model.f90:
    !
    !   field_id(1) represents "TAUX"   wind stress component
    !   field_id(2) represents "TAUY"   wind stress component
    !   field_id(3) represents "SFWFLX" surface fresh water flux
    !   field_id(4) represents "SFTEMP" surface temperature
    !   field_id(5) represents "THFLX"  total heat flux
    !   field_id(6) represents "ICEATM" ice temperatures and melt potential
    !
    !   field_id(7) represents "SST"    sea surface temperature
    !   field_id(9) represents "OCEANU" u component of ocean surface current
    !   field_id(9) represents "OCEANV" v component of ocean surface current
    !   field_id(10)represents "ICEOCE" ice thickness, concentration and temperatures
    !
#ifndef YAC_coupling
    field_shape(1) = 1
    field_shape(2) = nbr_hor_cells
    field_shape(3) = 1
#endif
    !
    ! Send fields away
    ! ----------------
    !
    write_coupler_restart = .FALSE.
    !
    ! TAUX
    !
    ! buffer(:,1)     = RESHAPE ( prm_field(jg)%u_stress_tile(:,:,iwtr), (/ nbr_cells /) )
    ! buffer(:,2)     = RESHAPE ( prm_field(jg)%u_stress_tile(:,:,iice), (/ nbr_cells /) )
    !
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      DO n = 1, nproma
         buffer(nn+n,1) = prm_field(jg)%u_stress_tile(n,i_blk,iwtr)
         buffer(nn+n,2) = prm_field(jg)%u_stress_tile(n,i_blk,iice)
      ENDDO
    ENDDO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
#ifdef YAC_coupling
    CALL yac_fput ( field_id(1), nbr_hor_cells, 2, 1, 1, buffer(1:nbr_hor_cells,1:2), info, ierror )
    IF ( info > 1 ) write_coupler_restart = .TRUE.
#else
    field_shape(3) = 2
    CALL ICON_cpl_put ( field_id(1), field_shape, buffer(1:nbr_hor_cells,1:2), info, ierror )
    IF ( info == 2 ) write_coupler_restart = .TRUE.
#endif
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    !
    ! TAUY
    !
    ! buffer(:,1)     = RESHAPE ( prm_field(jg)%v_stress_tile(:,:,iwtr), (/ nbr_cells /) )
    ! buffer(:,2)     = RESHAPE ( prm_field(jg)%v_stress_tile(:,:,iice), (/ nbr_cells /) )
    !
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      DO n = 1, nproma
         buffer(nn+n,1) = prm_field(jg)%v_stress_tile(n,i_blk,iwtr)
         buffer(nn+n,2) = prm_field(jg)%v_stress_tile(n,i_blk,iice)
      ENDDO
    ENDDO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
#ifdef YAC_coupling
    CALL yac_fput ( field_id(2), nbr_hor_cells, 2, 1, 1, buffer(1:nbr_hor_cells,1:2), info, ierror )
    IF ( info > 1 ) write_coupler_restart = .TRUE.
#else
    CALL ICON_cpl_put ( field_id(2), field_shape, buffer(1:nbr_hor_cells,1:2), info, ierror )
    IF ( info == 2 ) write_coupler_restart = .TRUE.
#endif
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    !
    ! SFWFLX Note: the evap_tile should be properly updated and added
    !
    !       write(0,*)  prm_field(jg)%rsfl(:,:)
    !       write(0,*)  prm_field(jg)%rsfc(:,:)
    !       write(0,*)  prm_field(jg)%ssfl(:,:)
    !       write(0,*)  prm_field(jg)%ssfc(:,:)
    !       write(0,*)  prm_field(jg)%evap_tile(:,:,iwtr)
    !
    ! buffer(:,1)     = RESHAPE ( prm_field(jg)%rsfl(:,:), (/ nbr_cells /) ) + &
    !   &               RESHAPE ( prm_field(jg)%rsfc(:,:), (/ nbr_cells /) ) ! total rain
    ! buffer(:,2)     = RESHAPE ( prm_field(jg)%ssfl(:,:), (/ nbr_cells /) ) + &
    !   &               RESHAPE ( prm_field(jg)%ssfc(:,:), (/ nbr_cells /) ) ! total snow
    ! buffer(:,3)     = RESHAPE ( prm_field(jg)%evap_tile(:,:,iwtr), (/ nbr_cells /) )
    !
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      DO n = 1, nproma
        buffer(nn+n,1) = prm_field(jg)%rsfl(n,i_blk) + prm_field(jg)%rsfc(n,i_blk) ! total rain
        buffer(nn+n,2) = prm_field(jg)%ssfl(n,i_blk) + prm_field(jg)%ssfc(n,i_blk) ! total snow
        buffer(nn+n,3) = prm_field(jg)%evap_tile(n,i_blk,iwtr)
      ENDDO
    ENDDO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
#ifdef YAC_coupling
    CALL yac_fput ( field_id(3), nbr_hor_cells, 3, 1, 1, buffer(1:nbr_hor_cells,1:3), info, ierror )
    IF ( info > 1 ) write_coupler_restart = .TRUE.
#else
    field_shape(3)  = 3
    CALL ICON_cpl_put ( field_id(3), field_shape, buffer(1:nbr_hor_cells,1:3), info, ierror )
    IF ( info == 2 ) write_coupler_restart = .TRUE.
#endif
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    !
    ! SFTEMP
    !
    ! buffer(:,1) =  RESHAPE ( prm_field(jg)%temp(:,nlev,:), (/ nbr_cells /) )
    !
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      DO n = 1, nproma
        buffer(nn+n,1) = prm_field(jg)%temp(n,nlev,i_blk)
      ENDDO
    ENDDO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
#ifdef YAC_coupling
    CALL yac_fput ( field_id(4), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info > 1 ) write_coupler_restart = .TRUE.
#else
    field_shape(3) = 1
    CALL ICON_cpl_put ( field_id(4), field_shape, buffer(1:nbr_hor_cells,1:1), info, ierror )
    IF ( info == 2 ) write_coupler_restart = .TRUE.
#endif
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    !
    ! THFLX, total heat flux
    !
    ! buffer(:,1)     =  RESHAPE ( prm_field(jg)%swflxsfc_tile(:,:,iwtr), (/ nbr_cells /) ) !net shortwave flux for ocean
    ! buffer(:,2)     =  RESHAPE ( prm_field(jg)%lwflxsfc_tile(:,:,iwtr), (/ nbr_cells /) ) !net longwave flux
    ! buffer(:,3)     =  RESHAPE ( prm_field(jg)%shflx_tile(:,:,iwtr),    (/ nbr_cells /) ) !sensible heat flux
    ! buffer(:,4)     =  RESHAPE ( prm_field(jg)%lhflx_tile(:,:,iwtr),    (/ nbr_cells /) ) !latent heat flux for ocean
    !
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      DO n = 1, nproma
        buffer(nn+n,1) = prm_field(jg)%swflxsfc_tile(n,i_blk,iwtr)
        buffer(nn+n,2) = prm_field(jg)%lwflxsfc_tile(n,i_blk,iwtr)
        buffer(nn+n,3) = prm_field(jg)%shflx_tile   (n,i_blk,iwtr)
        buffer(nn+n,4) = prm_field(jg)%lhflx_tile   (n,i_blk,iwtr)
      ENDDO
    ENDDO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
#ifdef YAC_coupling
    CALL yac_fput ( field_id(5), nbr_hor_cells, 4, 1, 1, buffer(1:nbr_hor_cells,1:4), info, ierror )
    IF ( info > 1 ) write_coupler_restart = .TRUE.
#else
    field_shape(3)  = 4
    CALL ICON_cpl_put ( field_id(5), field_shape, buffer(1:nbr_hor_cells,1:4), info, ierror )
    IF ( info == 2 ) write_coupler_restart = .TRUE.
#endif
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    !
    ! ICEATM, Ice state determined by atmosphere
    !
    ! buffer(:,1)     =  RESHAPE ( prm_field(jg)%Qtop(:,1,:), (/ nbr_cells /) ) !Melt-potential for ice - top
    ! buffer(:,2)     =  RESHAPE ( prm_field(jg)%Qbot(:,1,:), (/ nbr_cells /) ) !Melt-potential for ice - bottom
    ! buffer(:,3)     =  RESHAPE ( prm_field(jg)%T1  (:,1,:), (/ nbr_cells /) ) !Temperature of upper ice layer
    ! buffer(:,4)     =  RESHAPE ( prm_field(jg)%T2  (:,1,:), (/ nbr_cells /) ) !Temperature of lower ice layer
    !
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      DO n = 1, nproma
        buffer(nn+n,1) = prm_field(jg)%Qtop(n,1,i_blk)
        buffer(nn+n,2) = prm_field(jg)%Qbot(n,1,i_blk)
        buffer(nn+n,3) = prm_field(jg)%T1  (n,1,i_blk)
        buffer(nn+n,4) = prm_field(jg)%T2  (n,1,i_blk)
      ENDDO
    ENDDO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)
#ifdef YAC_coupling
    CALL yac_fput ( field_id(6), nbr_hor_cells, 4, 1, 1, buffer(1:nbr_hor_cells,1:4), info, ierror )
    IF ( info > 1 ) write_coupler_restart = .TRUE.
#else
    field_shape(3)  = 4
    CALL ICON_cpl_put ( field_id(6), field_shape, buffer(1:nbr_hor_cells,1:4), info, ierror )
    IF ( info == 2 ) write_coupler_restart = .TRUE.
#endif
    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( write_coupler_restart ) THEN
#ifdef YAC_coupling
       CALL warning('interface_echam_ocean', 'YAC says it is put for restart')
#else
       WRITE ( 6 , * ) "interface_echam_ocean: cpl layer says it is put for restart"
#endif
    ENDIF
    !
    ! Receive fields, only assign values if something was received ( info > 0 )
    ! -------------------------------------------------------------------------
    !
    buffer(nbr_hor_cells+1:nbr_cells,1:5) = 0.0_wp
    !
    ! SST
    !
    IF (ltimer) CALL timer_start(timer_coupling_1stget)
#ifdef YAC_coupling
    CALL yac_fget ( field_id(7), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    if ( info > 1 ) CALL warning('interface_echam_ocean', 'YAC says it is get for restart')
#else
    field_shape(3) = 1
    CALL ICON_cpl_get ( field_id(7), field_shape, buffer(1:nbr_hor_cells,1:1), info, ierror )
#endif
    IF (ltimer) CALL timer_stop(timer_coupling_1stget)
    !
    IF ( info > 0 ) THEN
      !
      ! prm_field(jg)%tsfc_tile(:,:,iwtr) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
      !
      buffer(nbr_inner_cells+1:nbr_cells,1) = 0.0_wp
      !
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        DO n = 1, nproma
          prm_field(jg)%tsfc_tile(n,i_blk,iwtr) = buffer(nn+n,1)
        ENDDO
      ENDDO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%tsfc_tile(:,:,iwtr))
    END IF
    !
    ! OCEANU
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)
#ifdef YAC_coupling
    CALL yac_fget ( field_id(8), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    if ( info > 1 ) CALL warning('interface_echam_ocean', 'YAC says it is get for restart')
#else
    CALL ICON_cpl_get ( field_id(8), field_shape, buffer(1:nbr_hor_cells,1:1), info, ierror )
#endif
    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > 0 ) THEN
      !
      ! prm_field(jg)%ocu(:,:) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
      !
      buffer(nbr_inner_cells+1:nbr_cells,1) = 0.0_wp
      !
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        DO n = 1, nproma
          prm_field(jg)%ocu(n,i_blk) = buffer(nn+n,1)
        ENDDO
      ENDDO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocu(:,:))
    END IF
    !
    ! OCEANV
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)
#ifdef YAC_coupling
    CALL yac_fget ( field_id(9), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    if ( info > 1 ) CALL warning('interface_echam_ocean', 'YAC says it is get for restart')
#else
    CALL ICON_cpl_get ( field_id(9), field_shape, buffer(1:nbr_hor_cells,1:1), info, ierror )
#endif
    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > 0 ) THEN
      !
      ! prm_field(jg)%ocv(:,:) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
      !
      buffer(nbr_inner_cells+1:nbr_cells,1) = 0.0_wp
      !
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        DO n = 1, nproma
          prm_field(jg)%ocv(n,i_blk) = buffer(nn+n,1)
        ENDDO
      ENDDO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocv(:,:))
    END IF
    !
    ! ICEOCE
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)
#ifdef YAC_coupling
    CALL yac_fget ( field_id(10), nbr_hor_cells, 5, 1, 1, buffer(1:nbr_hor_cells,1:5), info, ierror )
    if ( info > 1 ) CALL warning('interface_echam_ocean', 'YAC says it is get for restart')
#else
    field_shape(3) = 5
    CALL ICON_cpl_get ( field_id(10), field_shape, buffer(1:nbr_hor_cells,1:5), info, ierror )
#endif
    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > 0 ) THEN
      !
      ! prm_field(jg)%hi  (:,1,:) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
      ! prm_field(jg)%hs  (:,1,:) = RESHAPE (buffer(:,2), (/ nproma, p_patch%nblks_c /) )
      ! prm_field(jg)%conc(:,1,:) = RESHAPE (buffer(:,3), (/ nproma, p_patch%nblks_c /) )
      ! prm_field(jg)%T1  (:,1,:) = RESHAPE (buffer(:,4), (/ nproma, p_patch%nblks_c /) )
      ! prm_field(jg)%T2  (:,1,:) = RESHAPE (buffer(:,5), (/ nproma, p_patch%nblks_c /) )
      !
      buffer(nbr_inner_cells+1:nbr_cells,1:5) = 0.0_wp
      !
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        DO n = 1, nproma
          prm_field(jg)%hi  (n,1,i_blk) = buffer(nn+n,1)
          prm_field(jg)%hs  (n,1,i_blk) = buffer(nn+n,2)
          prm_field(jg)%conc(n,1,i_blk) = buffer(nn+n,3)
          prm_field(jg)%T1  (n,1,i_blk) = buffer(nn+n,4)
          prm_field(jg)%T2  (n,1,i_blk) = buffer(nn+n,5)
        ENDDO
      ENDDO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hi  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hs  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%conc(:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%T1  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%T2  (:,1,:))
      prm_field(jg)%seaice(:,:) = prm_field(jg)%conc(:,1,:)

    END IF

  END SUBROUTINE interface_echam_ocean

  !>
  !! SUBROUTINE destruct_echam_ocean_coupling -- terminates the coupling
  !! between ECHAM physics and the ocean.
  !!
  !! This subroutine is called at the end of the time loop of the ICONAM model.

  SUBROUTINE destruct_atmo_coupler

    IF ( .NOT. is_coupled_run() ) RETURN

    DEALLOCATE(buffer)

#ifdef YAC_coupling
    CALL yac_ffinalize
#else
    CALL icon_cpl_finalize ()
#endif

  END SUBROUTINE destruct_atmo_coupler
  
END MODULE mo_interface_echam_ocean
