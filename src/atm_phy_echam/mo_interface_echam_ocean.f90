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

#ifdef YAC_coupling

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_interface_echam_ocean

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_echam_phy_memory    ,ONLY: prm_field
  USE mo_echam_phy_config    ,ONLY: echam_phy_config
                                
  USE mo_parallel_config     ,ONLY: nproma
  
  USE mo_run_config          ,ONLY: ltimer, ico2, nlev
  USE mo_timer,               ONLY: timer_start, timer_stop,                &
       &                            timer_coupling_put, timer_coupling_get, &
       &                            timer_coupling_1stget, timer_coupling_init
  USE mo_echam_sfc_indices   ,ONLY: iwtr, iice, ilnd, nsfc_type

  USE mo_sync                ,ONLY: SYNC_C, sync_patch_array
  USE mo_impl_constants      ,ONLY: MAX_CHAR_LENGTH

  USE mo_ext_data_state      ,ONLY: ext_data

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)
  USE mo_interface_hd_ocean  ,ONLY: jsb_fdef_hd_fields
#endif

  USE mo_master_control      ,ONLY: get_my_process_name

  USE mo_mpi                 ,ONLY: p_pe_work
  USE mo_math_constants      ,ONLY: pi
  USE mo_parallel_config     ,ONLY: nproma

  USE mo_coupling_config     ,ONLY: is_coupled_run
  USE mo_time_config         ,ONLY: time_config
  
  USE mo_model_domain        ,ONLY: t_patch

  USE mo_exception           ,ONLY: warning, finish, message

  USE mo_yac_finterface      ,ONLY: yac_fput, yac_fget, yac_fget_version,        &
    &                               yac_fget_nbr_fields, yac_fget_field_ids,     &
    &                               yac_finit, yac_fdef_comp,                    &
    &                               yac_fdef_datetime,                           &
    &                               yac_fdef_subdomain, yac_fconnect_subdomains, &
    &                               yac_fdef_elements, yac_fdef_points,          &
    &                               yac_fdef_mask, yac_fdef_field, yac_fsearch,  &
    &                               yac_ffinalize, YAC_LOCATION_CELL

  USE mtime                  ,ONLY: datetimeToString, MAX_DATETIME_STR_LEN
  
  USE mo_util_dbg_prnt       ,ONLY: dbg_print

  USE mo_physical_constants  ,ONLY: amd, amco2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interface_echam_ocean
  PUBLIC :: construct_atmo_coupler, destruct_atmo_coupler

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_interface_echam_ocean'

  INTEGER, PARAMETER    :: no_of_fields = 12
  INTEGER               :: field_id(no_of_fields)

  REAL(wp), ALLOCATABLE :: buffer(:,:)
  INTEGER, SAVE         :: nbr_inner_cells

  CHARACTER(len=12)     :: str_module    = 'InterFaceOce'  ! Output of module for 1 line debug

CONTAINS

  !>
  !! SUBROUTINE construct_atmo_coupler -- the initialisation for the coupling
  !! of ECHAM physics and the ocean, through a coupler

  SUBROUTINE construct_atmo_coupler (p_patch)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)

    CHARACTER(LEN=MAX_CHAR_LENGTH)    ::  field_name(no_of_fields)

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
    INTEGER :: no_of_fields_total

    INTEGER :: mask_checksum
    INTEGER :: nblks
    INTEGER :: BLOCK, idx, INDEX
    INTEGER :: nbr_vertices_per_cell

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER,  ALLOCATABLE :: buffer_c(:,:)
    INTEGER,  ALLOCATABLE :: ibuffer(:)
    INTEGER,  ALLOCATABLE :: field_ids_total(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring

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

    ! Print the YAC version
    CALL message('Running ICON atmosphere in coupled mode with YAC version ', TRIM(yac_fget_version()) )

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
!ICON_OMP_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_RUNTIME_SCHEDULE
    DO BLOCK = 1, patch_horz%nblks_v
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer_lon(INDEX) = patch_horz%verts%vertex(idx,BLOCK)%lon * deg
        buffer_lat(INDEX) = patch_horz%verts%vertex(idx,BLOCK)%lat * deg
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT

!ICON_OMP_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_RUNTIME_SCHEDULE
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

!ICON_OMP_PARALLEL_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_RUNTIME_SCHEDULE
    DO BLOCK = 1, patch_horz%nblks_c
      DO idx = 1, nproma
        INDEX = (BLOCK-1)*nproma+idx
        buffer_lon(INDEX) = patch_horz%cells%center(idx,BLOCK)%lon * deg
        buffer_lat(INDEX) = patch_horz%cells%center(idx,BLOCK)%lat * deg
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    ! center points in cells (needed e.g. for patch recovery and nearest neighbour interpolation)
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
!ICON_OMP_PARALLEL_DO PRIVATE(idx) REDUCTION(+:nbr_inner_cells) ICON_OMP_RUNTIME_SCHEDULE
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
    ! The integer land-sea mask:
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! This integer mask for the atmosphere is available in ext_data(1)%atm%lsm_ctr_c(:,:).
    ! The (fractional) mask which is used in the ECHAM physics is prm_field(1)%lsmask(:,:).
    !
    ! The logical mask for the coupler must be generated from the fractional mask by setting
    !   only those gridpoints to land that have no ocean part at all (lsf<1 is ocean).
    ! The logical mask is then set to .FALSE. for land points to exclude them from mapping by yac.
    ! These points are not touched by yac.
    !

    mask_checksum = 0
!ICON_OMP_PARALLEL_DO PRIVATE(BLOCK,idx) REDUCTION(+:mask_checksum) ICON_OMP_RUNTIME_SCHEDULE
    DO BLOCK = 1, patch_horz%nblks_c
      DO idx = 1, nproma
        mask_checksum = mask_checksum + ABS(ext_data(1)%atm%lsm_ctr_c(idx, BLOCK))
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    !
    ! Define cell_mask_ids(1): all ocean and coastal points are valid
    !   This is the standard for the coupling fields listed below
    !
    IF ( mask_checksum > 0 ) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_RUNTIME_SCHEDULE
       DO BLOCK = 1, patch_horz%nblks_c
          DO idx = 1, nproma
             IF ( ext_data(1)%atm%lsm_ctr_c(idx, BLOCK) < 0 ) THEN
               ! Ocean point (lsm_ctr_c = -1 or -2) is valid
               ibuffer((BLOCK-1)*nproma+idx) = 0
             ELSE
               ! Land point (lsm_ctr_c = 1 or 2) is undef
               ibuffer((BLOCK-1)*nproma+idx) = 1
             ENDIF
          ENDDO
       ENDDO
!ICON_OMP_END_PARALLEL_DO
    ELSE
!ICON_OMP_PARALLEL_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_RUNTIME_SCHEDULE
       DO idx = 1,patch_horz%nblks_c * nproma
          ibuffer(idx) = 0
       ENDDO
!ICON_OMP_END_PARALLEL_DO
    ENDIF

    CALL yac_fdef_mask (           &
      & patch_horz%n_patch_cells,  &
      & ibuffer,                   &
      & cell_point_ids(1),         &
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
    field_name(11) = "co2_mixing_ratio"
    field_name(12) = "co2_flux"

    DO idx = 1, no_of_fields
      CALL yac_fdef_field (      &
        & TRIM(field_name(idx)), &
        & comp_id,               &
        & domain_id,             &
        & cell_point_ids,        &
        & cell_mask_ids(1),      &
        & 1,                     &
        & field_id(idx) )
    ENDDO

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)
    !
    ! Define cell_mask_ids(2) for runoff:
    !   Ocean coastal points with respect to HDmodel mask only are valid.
    !   The integer mask for the HDmodel is ext_data(1)%atm%lsm_hd_c(:,:).
    !   Caution: jg=1 is only valid for coupling to ocean
    !
    IF ( mask_checksum > 0 ) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(BLOCK, idx, INDEX) ICON_OMP_RUNTIME_SCHEDULE
        DO BLOCK = 1, patch_horz%nblks_c
          DO idx = 1, nproma
             IF ( ext_data(1)%atm%lsm_hd_c(idx, BLOCK) == -1 ) THEN
!            write(0,'(a,3i10)') 'BLOCK,IDX,SLM:', block,idx,ocean_coast(idx,block)
!            ibuffer((BLOCK-1)*nproma+idx) = ocean_coast(idx,BLOCK)
!            IF ( prm_field(1)%hdmask(idx, BLOCK) < -0.9_wp .AND. &
!              &  prm_field(1)%hdmask(idx, BLOCK) > -1.1_wp) THEN
                ! Ocean point at coast is valid
                ibuffer((BLOCK-1)*nproma+idx) = 0
             ELSE
                ! Land point or ocean point without coast is undef
                ibuffer((BLOCK-1)*nproma+idx) = 1
             ENDIF
          ENDDO
        ENDDO
!ICON_OMP_END_PARALLEL_DO

        CALL yac_fdef_mask (           &
          & patch_horz%n_patch_cells,  &
          & ibuffer,                   &
          & cell_point_ids(1),         &
          & cell_mask_ids(2) )

      ENDIF

      ! Define additional coupling field(s) for JSBACH/HD
      ! Utilize mask field for runoff
      !  - cell_mask_ids(2:2) is ocean coast points only for source point mapping (source_to_target_map)
      CALL jsb_fdef_hd_fields(comp_id, domain_id, cell_point_ids, cell_mask_ids(2:2))

#endif

    DEALLOCATE (ibuffer)

    ! End definition of coupling fields and search
    CALL yac_fget_nbr_fields(no_of_fields_total)
    ALLOCATE(field_ids_total(no_of_fields_total))
    CALL yac_fget_field_ids(no_of_fields_total, field_ids_total)
    CALL yac_fsearch ( 1, comp_ids, no_of_fields_total, field_ids_total, error_status )

    ALLOCATE(buffer(nproma*patch_horz%nblks_c,5))

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

  SUBROUTINE interface_echam_ocean( p_patch ) ! in

    ! Arguments

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: nbr_hor_cells  ! = inner and halo points
    INTEGER               :: jg             ! grid index
    INTEGER               :: n              ! nproma loop count
    INTEGER               :: nn             ! block offset
    INTEGER               :: i_blk          ! block loop count
    INTEGER               :: nlen           ! nproma/npromz
    INTEGER               :: info, ierror   !< return values from cpl_put/get calls
    INTEGER               :: no_arr         !  no of arrays in bundle for put/get calls

    REAL(wp), PARAMETER   :: dummy = 0.0_wp
    REAL(wp)              :: scr(nproma,p_patch%alloc_cell_blocks)
    REAL(wp)              :: frac_oce, fwf_fac

    IF ( .NOT. is_coupled_run() ) RETURN

    jg = p_patch%id

    !-------------------------------------------------------------------------
    ! If running in atm-oce coupled mode, exchange information 
    !-------------------------------------------------------------------------

    ! Possible fields that contain information to be sent to the ocean include
    !
    ! 1. prm_field(jg)% u_stress_tile(:,:,iwtr/iice)  and 
    !    prm_field(jg)% v_stress_tile(:,:,iwtr/iice)  which are the wind stress components over water and ice respectively
    !
    ! 2. prm_field(jg)% evap_tile(:,:,iwtr/iice)  evaporation rate over ice-covered and open ocean/lakes, no land;
    !
    ! 3. prm_field(jg)%rsfl + prm_field(jg)%rsfc + prm_field(jg)%ssfl + prm_field(jg)%ssfc
    !    which gives the precipitation rate;
    !
    ! 4. prm_field(jg)% ta(:,nlev,:)  temperature at the lowest model level, or
    !    prm_field(jg)% tas(:,:)      2-m temperature, not available yet, or
    !    prm_field(jg)% shflx_tile(:,:,iwtr) sensible heat flux
    !    ... tbc
    !
    ! 5  prm_field(jg)% lhflx_tile(:,:,iwtr) latent heat flux
    ! 6. shortwave radiation flux at the surface
    !
    ! Possible fields to receive from the ocean include
    !
    ! 1. prm_field(jg)% ts_tile(:,:,iwtr)   SST
    ! 2. prm_field(jg)% ocu(:,:) and ocv(:,:) ocean surface current
    ! 3. ... tbc
    ! 

    nbr_hor_cells = p_patch%n_patch_cells

    !
    !  Send fields to ocean:
    !   field_id(1) represents "surface_downward_eastward_stress" bundle  - zonal wind stress component over ice and water
    !   field_id(2) represents "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !   field_id(3) represents "surface_fresh_water_flux" bundle          - liquid rain, snowfall, evaporation
    !   field_id(4) represents "total heat flux" bundle                   - short wave, long wave, sensible, latent heat flux
    !   field_id(5) represents "atmosphere_sea_ice_bundle"                - sea ice surface and bottom melt potentials
    !   field_id(10) represents "10m_wind_speed"                          - atmospheric wind speed
    !   field_id(11) represents "qtrc(nlev,co2)"                          - co2 mixing ratio
    !
    !  Receive fields from ocean:
    !   field_id(6) represents "sea_surface_temperature"                  - SST
    !   field_id(7) represents "eastward_sea_water_velocity"              - zonal velocity, u component of ocean surface current
    !   field_id(8) represents "northward_sea_water_velocity"             - meridional velocity, v component of ocean surface current
    !   field_id(9) represents "ocean_sea_ice_bundle"                     - ice thickness, snow thickness, ice concentration
    !   field_id(12) represents "co2_flux"                                - ocean co2 flux
    !

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from atmosphere to ocean
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    write_coupler_restart = .FALSE.

    !
    ! ------------------------------
    !  Send zonal wind stress bundle
    !   field_id(1) represents "surface_downward_eastward_stress" bundle - zonal wind stress component over ice and water
    !
    buffer(:,:) = 0.0_wp  ! temporarily
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
         buffer(nn+n,1) = prm_field(jg)%u_stress_tile(n,i_blk,iwtr)
         buffer(nn+n,2) = prm_field(jg)%u_stress_tile(n,i_blk,iice)
      ENDDO
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 2
    CALL yac_fput ( field_id(1), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fput called after end of run - id=1, u-stress')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send meridional wind stress bundle
    !   field_id(2) represents "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !
    buffer(:,:) = 0.0_wp  ! temporarily
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
         buffer(nn+n,1) = prm_field(jg)%v_stress_tile(n,i_blk,iwtr)
         buffer(nn+n,2) = prm_field(jg)%v_stress_tile(n,i_blk,iice)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 2
    CALL yac_fput ( field_id(2), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fput called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send surface fresh water flux bundle
    !   field_id(3) represents "surface_fresh_water_flux" bundle - liquid rain, snowfall, evaporation
    !
    !   Note: the evap_tile should be properly updated and added;
    !         as long as evaporation over sea-ice is not used in ocean thermodynamics, the evaporation over the
    !         whole ocean part of grid-cell is passed to the ocean
    !         for pre04 a preliminary solution for evaporation in ocean model is to exclude the land fraction
    !         evap.oce = (evap.wtr*frac.wtr + evap.ice*frac.ice)/(1-frac.lnd)
    !
    buffer(:,:) = 0.0_wp  ! temporarily
    !
    ! Preliminary: hard-coded correction factor for freshwater imbalance stemming from the atmosphere
    ! Precipitation is reduced by Factor fwf_fac
    ! factor calculated from run slo1014, used in run slo1016:
    ! Global imbalance D=1.7 mm/y; Precip P=1070 mm/y; D/P~0.0016; 1-D/P= 0.9984
    ! fwf_fac = 0.9984_wp   
    fwf_fac = 1.0_wp      ! neutral factor

    ! Aquaplanet coupling: surface types ocean and ice only
    IF (nsfc_type == 2) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
     
          ! total rates of rain and snow over whole cell
          !buffer(nn+n,1) = prm_field(jg)%rsfl(n,i_blk) + prm_field(jg)%rsfc(n,i_blk)
          !buffer(nn+n,2) = prm_field(jg)%ssfl(n,i_blk) + prm_field(jg)%ssfc(n,i_blk)
          buffer(nn+n,1) = (prm_field(jg)%rsfl(n,i_blk) + prm_field(jg)%rsfc(n,i_blk))*fwf_fac
          buffer(nn+n,2) = (prm_field(jg)%ssfl(n,i_blk) + prm_field(jg)%ssfc(n,i_blk))*fwf_fac
     
          ! evaporation over ice-free and ice-covered water fraction - of whole ocean part
          frac_oce = prm_field(jg)%frac_tile(n,i_blk,iwtr) + prm_field(jg)%frac_tile(n,i_blk,iice) ! 1.0?
          buffer(nn+n,3) = prm_field(jg)%evap_tile(n,i_blk,iwtr)*prm_field(jg)%frac_tile(n,i_blk,iwtr) + &
            &              prm_field(jg)%evap_tile(n,i_blk,iice)*prm_field(jg)%frac_tile(n,i_blk,iice)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ! Full coupling including jsbach: surface types ocean, ice, land
    ELSE IF (nsfc_type == 3) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
    
          ! total rates of rain and snow over whole cell
          !buffer(nn+n,1) = prm_field(jg)%rsfl(n,i_blk) + prm_field(jg)%rsfc(n,i_blk)
          !buffer(nn+n,2) = prm_field(jg)%ssfl(n,i_blk) + prm_field(jg)%ssfc(n,i_blk)
          buffer(nn+n,1) = (prm_field(jg)%rsfl(n,i_blk) + prm_field(jg)%rsfc(n,i_blk))*fwf_fac
          buffer(nn+n,2) = (prm_field(jg)%ssfl(n,i_blk) + prm_field(jg)%ssfc(n,i_blk))*fwf_fac
    
          ! evaporation over ice-free and ice-covered water fraction, of whole ocean part, without land part
          frac_oce=1.0_wp-prm_field(jg)%frac_tile(n,i_blk,ilnd)
          IF (frac_oce <= 0.0_wp) THEN
            ! land part is zero
            buffer(nn+n,3) = 0.0_wp
          ELSE
            buffer(nn+n,3) = (prm_field(jg)%evap_tile(n,i_blk,iwtr)*prm_field(jg)%frac_tile(n,i_blk,iwtr) + &
              &               prm_field(jg)%evap_tile(n,i_blk,iice)*prm_field(jg)%frac_tile(n,i_blk,iice))/frac_oce
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    ELSE
      CALL finish('interface_echam_ocean: coupling only for nsfc_type equals 2 or 3. Check your code/configuration!')
    ENDIF  !  nsfc_type
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 3
    CALL yac_fput ( field_id(3), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fput called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send total heat flux bundle
    !   field_id(4) represents "total heat flux" bundle - short wave, long wave, sensible, latent heat flux
    !
    buffer(:,:) = 0.0_wp  ! temporarily
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
        buffer(nn+n,1) = prm_field(jg)%swflxsfc_tile(n,i_blk,iwtr)
        buffer(nn+n,2) = prm_field(jg)%lwflxsfc_tile(n,i_blk,iwtr)
        buffer(nn+n,3) = prm_field(jg)%shflx_tile   (n,i_blk,iwtr)
        buffer(nn+n,4) = prm_field(jg)%lhflx_tile   (n,i_blk,iwtr)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 4
    CALL yac_fput ( field_id(4), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fput called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    !
    ! ------------------------------
    !  Send sea ice flux bundle
    !   field_id(5) represents "atmosphere_sea_ice_bundle" - sea ice surface and bottom melt potentials Qtop, Qbot
    !
    buffer(:,:) = 0.0_wp  ! temporarily
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
        buffer(nn+n,1) = prm_field(jg)%Qtop(n,1,i_blk)
        buffer(nn+n,2) = prm_field(jg)%Qbot(n,1,i_blk)
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 2
    CALL yac_fput ( field_id(5), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fput called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( write_coupler_restart ) THEN
       CALL warning('interface_echam_ocean', 'YAC says it is put for restart - id=5, atmos sea ice bundle')
    ENDIF

    !
    ! ------------------------------
    !  Send 10m wind speed
    !   field_id(10) represents "10m_wind_speed" - atmospheric wind speed
    !
    buffer(:,:) = 0.0_wp
!!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
        ! as far as no tiles (pre04) are correctly implemented, use the grid-point mean of 10m wind for coupling
        buffer(nn+n,1) = prm_field(jg)%sfcWind(n,i_blk)
      ENDDO
    ENDDO
!!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 1
    CALL yac_fput ( field_id(10), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fput called after end of run - id=10, wind speed')

    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( write_coupler_restart ) THEN
       CALL warning('interface_echam_ocean', 'YAC says it is put for restart - id=10, wind speed')
    ENDIF

#ifndef __NO_ICON_OCEAN__
    IF(ANY(echam_phy_config(:)%lcpl_co2_atmoce))then
    !
    ! ------------------------------
    !  Send co2 mixing ratio
    !   field_id(11) represents "co2_mixing_ratio" - CO2 mixing ratio
    !
    buffer(:,:) = 0.0_wp
!!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      nn = (i_blk-1)*nproma
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      DO n = 1, nlen
        buffer(nn+n,1) = prm_field(jg)%qtrc(n,nlev,i_blk,ico2) * amd/amco2 * 1.0e6_wp
      ENDDO
    ENDDO
!!ICON_OMP_END_PARALLEL_DO
    !
    IF (ltimer) CALL timer_start(timer_coupling_put)

    no_arr = 1
    CALL yac_fput ( field_id(11), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.
    IF ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fput called after end of run - id=11, co2mmr')

    IF (ltimer) CALL timer_stop(timer_coupling_put)
    !
    IF ( write_coupler_restart ) THEN
       CALL warning('interface_echam_ocean', 'YAC says it is put for restart - id=11, co2mmr')
    ENDIF
    ENDIF
#endif



    !

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from ocean to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    !  Receive fields, only assign values if something was received ( info > 0 )
    !   - ocean fields have undefined values on land, which are not sent to the atmosphere,
    !     therefore buffer is set to zero to avoid unintended usage of ocean values over land
    !

    !
    ! ------------------------------
    !  Receive SST
    !   field_id(6) represents "sea_surface_temperature" - SST
    !
    IF (ltimer) CALL timer_start(timer_coupling_1stget)

    !buffer(:,:) = 0.0_wp
    ! buffer for tsfc in Kelvin
    buffer(:,:) = 199.99_wp

    CALL yac_fget ( field_id(6), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    if ( info > 1 .AND. info < 7 ) CALL warning('interface_echam_ocean', 'YAC says it is get for restart')
    if ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fget called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_1stget)
    !
    IF ( info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%ts_tile(n,i_blk,iwtr) = dummy
          ELSE
            prm_field(jg)%ts_tile(n,i_blk,iwtr) = buffer(nn+n,1)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ts_tile(:,:,iwtr))
    END IF
    !
    ! ------------------------------
    !  Receive zonal velocity
    !   field_id(7) represents "eastward_sea_water_velocity" - zonal velocity, u component of ocean surface current
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp
    CALL yac_fget ( field_id(7), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    if ( info > 1 .AND. info < 7 ) CALL warning('interface_echam_ocean', 'YAC says it is get for restart')
    if ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fget called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%ocu(n,i_blk) = dummy
          ELSE
            prm_field(jg)%ocu(n,i_blk) = buffer(nn+n,1)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocu(:,:))
    END IF
    !
    !
    ! ------------------------------
    !  Receive meridional velocity
    !   field_id(8) represents "northward_sea_water_velocity" - meridional velocity, v component of ocean surface current
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp
    CALL yac_fget ( field_id(8), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    if ( info > 1 .AND. info < 7 ) CALL warning('interface_echam_ocean', 'YAC says it is get for restart')
    if ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fget called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%ocv(n,i_blk) = dummy
          ELSE
            prm_field(jg)%ocv(n,i_blk) = buffer(nn+n,1)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocv(:,:))
    END IF
    !
    ! ------------------------------
    !  Receive sea ice bundle
    !   field_id(9) represents "ocean_sea_ice_bundle" - ice thickness, snow thickness, ice concentration
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp
    no_arr = 3
    CALL yac_fget ( field_id(9), nbr_hor_cells, no_arr, 1, 1, buffer(1:nbr_hor_cells,1:no_arr), info, ierror )
    if ( info > 1 .AND. info < 7 ) CALL warning('interface_echam_ocean', 'YAC says it is get for restart')
    if ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fget called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%hi  (n,1,i_blk) = dummy
            prm_field(jg)%hs  (n,1,i_blk) = dummy
            prm_field(jg)%conc(n,1,i_blk) = dummy
          ELSE
            prm_field(jg)%hi  (n,1,i_blk) = buffer(nn+n,1)
            prm_field(jg)%hs  (n,1,i_blk) = buffer(nn+n,2)
            prm_field(jg)%conc(n,1,i_blk) = buffer(nn+n,3)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hi  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hs  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%conc(:,1,:))
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          prm_field(jg)%seaice(n,i_blk) = prm_field(jg)%conc(n,1,i_blk)
          prm_field(jg)%siced(n,i_blk)  = prm_field(jg)%hi(n,1,i_blk)
!  snow thickness not yet connected
!!$          prm_field(jg)%...(n,i_blk)    = prm_field(jg)%hs(n,1,i_blk)
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      
    END IF
  !    !
    IF(ANY(echam_phy_config(:)%lcpl_co2_atmoce))then
    !
    ! ------------------------------
    !  Receive co2 flux
    !   field_id(12) represents "co2_flux" - ocean co2 flux
    !
    IF (ltimer) CALL timer_start(timer_coupling_get)

    buffer(:,:) = 0.0_wp
    CALL yac_fget ( field_id(12), nbr_hor_cells, 1, 1, 1, buffer(1:nbr_hor_cells,1:1), info, ierror )
    if ( info > 1 .AND. info < 7 ) CALL warning('interface_echam_ocean', 'YAC says it is get for restart')
    if ( info == 7 ) CALL warning('interface_echam_ocean', 'YAC says fget called after end of run')

    IF (ltimer) CALL timer_stop(timer_coupling_get)
    !
    IF ( info > 0 .AND. info < 7 ) THEN
      !
!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%co2_flux_tile(n,i_blk,iwtr) = dummy
          ELSE
            prm_field(jg)%co2_flux_tile(n,i_blk,iwtr) = buffer(nn+n,1)
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      !
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%co2_flux_tile(:,:,iwtr))
    ENDIF ! lcpl_co2_atmoce

    END IF    !---------DEBUG DIAGNOSTICS-------------------------------------------

    ! u/v-stress on ice and water sent
    scr(:,:) = prm_field(jg)%u_stress_tile(:,:,iwtr)
    CALL dbg_print('EchOce: u_stress.wtr',scr,str_module,3,in_subset=p_patch%cells%owned)
    scr(:,:) = prm_field(jg)%u_stress_tile(:,:,iice)
    CALL dbg_print('EchOce: u_stress.ice',scr,str_module,3,in_subset=p_patch%cells%owned)
    scr(:,:) = prm_field(jg)%v_stress_tile(:,:,iwtr)
    CALL dbg_print('EchOce: v_stress.wtr',scr,str_module,4,in_subset=p_patch%cells%owned)
    scr(:,:) = prm_field(jg)%v_stress_tile(:,:,iice)
    CALL dbg_print('EchOce: v_stress.ice',scr,str_module,4,in_subset=p_patch%cells%owned)

    ! rain, snow, evaporation
    scr(:,:) = prm_field(jg)%rsfl(:,:) + prm_field(jg)%rsfc(:,:)
    CALL dbg_print('EchOce: total rain  ',scr,str_module,3,in_subset=p_patch%cells%owned)
    scr(:,:) = prm_field(jg)%ssfl(:,:) + prm_field(jg)%ssfc(:,:)
    CALL dbg_print('EchOce: total snow  ',scr,str_module,4,in_subset=p_patch%cells%owned)
    CALL dbg_print('EchOce: evaporation ',prm_field(jg)%evap   ,str_module,4,in_subset=p_patch%cells%owned)

    ! short wave, long wave, sensible, latent heat flux sent
    scr(:,:) = prm_field(jg)%swflxsfc_tile(:,:,iwtr)
    CALL dbg_print('EchOce: swflxsfc.wtr',scr,str_module,3,in_subset=p_patch%cells%owned)
    scr(:,:) = prm_field(jg)%lwflxsfc_tile(:,:,iwtr)
    CALL dbg_print('EchOce: lwflxsfc.wtr',scr,str_module,3,in_subset=p_patch%cells%owned)
    scr(:,:) = prm_field(jg)%shflx_tile(:,:,iwtr)
    CALL dbg_print('EchOce: shflx.wtr   ',scr,str_module,3,in_subset=p_patch%cells%owned)
    scr(:,:) = prm_field(jg)%lhflx_tile(:,:,iwtr)
    CALL dbg_print('EchOce: lhflx.wtr   ',scr,str_module,3,in_subset=p_patch%cells%owned)

    ! Qtop and Qbot, windspeed sent
    !scr(:,:) = prm_field(jg)%Qtop(:,1,:)
    !CALL dbg_print('EchOce: u_stress.wtr',scr,str_module,3,in_subset=p_patch%cells%owned)
    CALL dbg_print('EchOce: ice-Qtop    ',prm_field(jg)%Qtop   ,str_module,4,in_subset=p_patch%cells%owned)
    CALL dbg_print('EchOce: ice-Qbot    ',prm_field(jg)%Qbot   ,str_module,3,in_subset=p_patch%cells%owned)
    CALL dbg_print('EchOce: sfcWind     ',prm_field(jg)%sfcWind,str_module,3,in_subset=p_patch%cells%owned)

    ! SST, sea ice, ocean velocity received
    scr(:,:) = prm_field(jg)%ts_tile(:,:,iwtr)
    CALL dbg_print('EchOce: ts_tile.iwtr',scr                  ,str_module,2,in_subset=p_patch%cells%owned)
    CALL dbg_print('EchOce: siced       ',prm_field(jg)%siced  ,str_module,3,in_subset=p_patch%cells%owned)
    CALL dbg_print('EchOce: seaice      ',prm_field(jg)%seaice ,str_module,4,in_subset=p_patch%cells%owned)
    scr(:,:) = prm_field(jg)%ocu(:,:)
    CALL dbg_print('EchOce: ocu         ',prm_field(jg)%ocu    ,str_module,4,in_subset=p_patch%cells%owned)
    CALL dbg_print('EchOce: ocv         ',prm_field(jg)%ocv    ,str_module,4,in_subset=p_patch%cells%owned)

    ! Fraction of tiles:
    scr(:,:) = prm_field(jg)%frac_tile(:,:,iwtr)
    CALL dbg_print('EchOce: frac_tile.wtr',scr,str_module,2,in_subset=p_patch%cells%owned)
    scr(:,:) = prm_field(jg)%frac_tile(:,:,iice)
    CALL dbg_print('EchOce: frac_tile.ice',scr,str_module,3,in_subset=p_patch%cells%owned)
    if (nsfc_type == 3) THEN
      scr(:,:) = prm_field(jg)%frac_tile(:,:,ilnd)
      CALL dbg_print('EchOce: frac_tile.lnd',scr,str_module,4,in_subset=p_patch%cells%owned)
    ENDIF

    !---------------------------------------------------------------------


  END SUBROUTINE interface_echam_ocean

  !>
  !! SUBROUTINE destruct_echam_ocean_coupling -- terminates the coupling
  !! between ECHAM physics and the ocean.
  !!
  !! This subroutine is called at the end of the time loop of the ICONAM model.

  SUBROUTINE destruct_atmo_coupler

    IF ( .NOT. is_coupled_run() ) RETURN

    DEALLOCATE(buffer)

    CALL yac_ffinalize

  END SUBROUTINE destruct_atmo_coupler
  
END MODULE mo_interface_echam_ocean

#else

MODULE mo_interface_echam_ocean

  USE mo_model_domain,    ONLY: t_patch
  USE mo_exception,       ONLY: finish
  USE mo_coupling_config, ONLY: is_coupled_run

  PUBLIC :: interface_echam_ocean
  PUBLIC :: construct_atmo_coupler, destruct_atmo_coupler

CONTAINS

  SUBROUTINE construct_atmo_coupler (p_patch)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)

    IF ( is_coupled_run() ) THEN
       CALL finish('construct_atmo_coupler: unintentionally called. Check your source code and configure.')
    ELSE
       RETURN
    ENDIF

  END SUBROUTINE construct_atmo_coupler

  SUBROUTINE interface_echam_ocean ( p_patch )

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch

    IF ( is_coupled_run() ) THEN
       CALL finish('interface_echam_ocean: unintentionally called. Check your source code and configure.')
    ELSE
       RETURN
    ENDIF

  END SUBROUTINE interface_echam_ocean

  SUBROUTINE destruct_atmo_coupler

    IF ( is_coupled_run() ) THEN
       CALL finish('destruct_atmo_coupler: unintentionally called. Check your source code and configure.')
    ELSE
       RETURN
    ENDIF

  END SUBROUTINE destruct_atmo_coupler

END MODULE mo_interface_echam_ocean

#endif
