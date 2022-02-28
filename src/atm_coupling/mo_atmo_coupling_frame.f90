!>
!! @brief Initialisation of atmosphere coupling
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

MODULE mo_atmo_coupling_frame

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_nonhydro_types      ,ONLY: t_nh_diag
  USE mo_ext_data_state      ,ONLY: ext_data
#ifndef __NO_ECHAM__
  USE mo_echam_phy_memory    ,ONLY: prm_field
#endif
  USE mo_lnd_nwp_config      ,ONLY: isub_water

  USE mo_parallel_config     ,ONLY: nproma

  USE mo_run_config          ,ONLY: iforcing, ltimer
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling_init

  USE mo_impl_constants      ,ONLY: MAX_CHAR_LENGTH, inwp, iecham

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)
  USE mo_interface_hd_ocean  ,ONLY: jsb_fdef_hd_fields
#endif

  USE mo_master_control      ,ONLY: get_my_process_name

  USE mo_mpi                 ,ONLY: p_pe_work, p_comm_work, p_sum
  USE mo_math_constants      ,ONLY: pi
  USE mo_parallel_config     ,ONLY: nproma

  USE mo_coupling_config     ,ONLY: is_coupled_run
  USE mo_time_config         ,ONLY: time_config

  USE mo_exception           ,ONLY: finish, message

  USE mo_yac_finterface      ,ONLY: yac_fget_version,                            &
    &                               yac_finit, yac_fdef_comp,                    &
    &                               yac_fdef_datetime,                           &
    &                               yac_fdef_grid, yac_fdef_points,              &
    &                               yac_fset_global_index, yac_fset_core_mask,   &
    &                               yac_fdef_mask, yac_fdef_field_mask,          &
    &                               yac_fsearch, yac_ffinalize,                  &
    &                               YAC_LOCATION_CELL, COUPLING, OUT_OF_BOUND

  USE mtime                  ,ONLY: datetimeToString, MAX_DATETIME_STR_LEN

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_atmo_coupling, destruct_atmo_coupling
  PUBLIC :: lyac_very_1st_get, nbr_inner_cells, mask_checksum, field_id

  INTEGER, PARAMETER    :: no_of_fields = 13
  INTEGER               :: field_id(no_of_fields)

  INTEGER, SAVE         :: nbr_inner_cells
  INTEGER, SAVE         :: mask_checksum
  LOGICAL, SAVE         :: lyac_very_1st_get

CONTAINS

  !>
  !! SUBROUTINE construct_atmo_coupling -- the initialisation for the coupling
  !! of atmosphere and the ocean, through a coupler

  SUBROUTINE construct_atmo_coupling (p_patch)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)
    CHARACTER(LEN=MAX_CHAR_LENGTH)    :: field_name(no_of_fields)

    INTEGER :: error_status

    TYPE(t_patch), POINTER :: patch_horz

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------

    CHARACTER(LEN=max_char_length) :: xml_filename
    CHARACTER(LEN=max_char_length) :: xsd_filename
    CHARACTER(LEN=max_char_length) :: grid_name
    CHARACTER(LEN=max_char_length) :: comp_name

    INTEGER :: comp_id
    INTEGER :: comp_ids(1)
    INTEGER :: cell_point_ids(1)
    INTEGER :: cell_mask_ids(1)
    INTEGER :: grid_id

    INTEGER :: jg
    INTEGER :: nblks
    INTEGER :: jb, jc, nn
    INTEGER :: nbr_vertices_per_cell

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER,  ALLOCATABLE :: buffer_c(:,:)
    LOGICAL,  ALLOCATABLE :: is_valid(:)

    REAL(wp), ALLOCATABLE :: lsmnolake(:,:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring

    ! Skip time measurement of the very first yac_fget
    ! as this will measure mainly the wait time caused
    ! by the initialisation of the model components
    ! and does not tell us much about the load balancing
    ! in subsequent calls.

    lyac_very_1st_get = .TRUE.

    IF ( .NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start (timer_coupling_init)

    comp_name = TRIM(get_my_process_name())

    jg = 1
    patch_horz => p_patch(jg)

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
 
    ! Announce one grid (patch) to the coupler
    grid_name = "icon_atmos_grid"

    ! Extract cell information
    !
    ! cartesian coordinates of cell vertices are stored in
    ! patch_horz%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes in rad.

    nblks = max(patch_horz%nblks_c,patch_horz%nblks_v)

    ALLOCATE(buffer_lon(nproma*nblks))
    ALLOCATE(buffer_lat(nproma*nblks))
    ALLOCATE(buffer_c(3,nproma*nblks))

    ALLOCATE(lsmnolake(nproma,nblks))

    nbr_vertices_per_cell = 3

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_v
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_lon(nn) = patch_horz%verts%vertex(jc,jb)%lon
        buffer_lat(nn) = patch_horz%verts%vertex(jc,jb)%lat
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT

!ICON_OMP_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_c
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_c(1,nn) = (patch_horz%cells%vertex_blk(jc,jb,1)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(jc,jb,1)
        buffer_c(2,nn) = (patch_horz%cells%vertex_blk(jc,jb,2)-1)*nproma + &
          &               patch_horz%cells%vertex_idx(jc,jb,2)
        buffer_c(3,nn) = (patch_horz%cells%vertex_blk(jc,jb,3)-1)*nproma + &
                          patch_horz%cells%vertex_idx(jc,jb,3)
      ENDDO
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    ! Description of elements, here as unstructured grid
    CALL yac_fdef_grid(           &
      & TRIM(grid_name),          &
      & patch_horz%n_patch_verts, &
      & patch_horz%n_patch_cells, &
      & nbr_vertices_per_cell,    &
      & buffer_lon,               &
      & buffer_lat,               &
      & buffer_c,                 &
      & grid_id)

    !
    ! Define cell center points (location = 0)
    !
    ! cartesian coordinates of cell centers are stored in
    ! patch_horz%cells%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_c
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_lon(nn) = patch_horz%cells%center(jc,jb)%lon
        buffer_lat(nn) = patch_horz%cells%center(jc,jb)%lat
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    ! center points in cells (needed e.g. for patch recovery and nearest neighbour interpolation)
    CALL yac_fdef_points (        &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & buffer_lon,               &
      & buffer_lat,               &
      & cell_point_ids(1) )

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c)

    CALL yac_fset_global_index (                &
      & patch_horz%cells%decomp_info%glb_index, &
      & YAC_LOCATION_CELL,                      &
      & grid_id )

    ALLOCATE(is_valid(nproma*patch_horz%nblks_c))

    nbr_inner_cells = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jc) REDUCTION(+:nbr_inner_cells) ICON_OMP_RUNTIME_SCHEDULE
    DO jc = 1, patch_horz%n_patch_cells
       IF ( p_pe_work == patch_horz%cells%decomp_info%owner_local(jc) ) THEN
         is_valid(jc) = .TRUE.
         nbr_inner_cells = nbr_inner_cells + 1
       ELSE
         is_valid(jc) = .FALSE.
       ENDIF
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    CALL yac_fset_core_mask ( &
      & is_valid,             &
      & YAC_LOCATION_CELL,    &
      & grid_id )

    !
    ! The integer land-sea mask:
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! The (fractional) mask which is used in the ECHAM physics is prm_field(1)%lsmask(:,:).
    !
    ! The logical mask for the coupler must be generated from the fractional mask by setting
    !   only those gridpoints to land that have no ocean part at all (lsf<1 is ocean).
    ! The logical mask is then set to .FALSE. for land points to exclude them from mapping by yac.
    ! These points are not touched by yac.
    !

    SELECT CASE( iforcing )

       CASE ( inwp )

 !ICON_OMP_PARALLEL_DO PRIVATE(jb,jc) ICON_OMP_RUNTIME_SCHEDULE
          DO jb = 1, patch_horz%nblks_c
             DO jc = 1, nproma
                !RR check 
                lsmnolake(jc, jb) = ext_data(jg)%atm%frac_t(jc,jb,isub_water) + ext_data(jg)%atm%fr_lake(jc,jb)
             ENDDO
          ENDDO
!ICON_OMP_END_PARALLEL_DO

       CASE ( iecham )
#ifdef __NO_ECHAM__
          CALL finish ('mo_atmo_coupling_frame:construct_atmo_coupling' &
            & //'coupled model needs echam; remove --disable-echam and reconfigure')
#else
!ICON_OMP_PARALLEL_DO PRIVATE(jb,jc) ICON_OMP_RUNTIME_SCHEDULE
          DO jb = 1, patch_horz%nblks_c
             DO jc = 1, nproma
                !  slo: caution - lsmask includes alake, must be added to refetch pure lsm:
                lsmnolake(jc, jb) = prm_field(1)%lsmask(jc,jb) + prm_field(1)%alake(jc,jb)
             ENDDO
          ENDDO
!ICON_OMP_END_PARALLEL_DO
#endif
       CASE DEFAULT

          CALL finish ('Please mask handling for new forcing in ' &
            & //'src/coupling/mo_atmo_coupling_frame: construct_atmo_coupling. Thank you!')

    END SELECT
    
    mask_checksum = 0
!ICON_OMP_PARALLEL_DO PRIVATE(jb,jc) REDUCTION(+:mask_checksum) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, patch_horz%nblks_c
      DO jc = 1, nproma
        mask_checksum = mask_checksum + ABS( lsmnolake(jc,jb))
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    mask_checksum = p_sum(mask_checksum, comm=p_comm_work)
    !
    ! Define cell_mask_ids(1): all ocean and coastal points are valid
    !   This is the standard for the coupling of atmospheric fields listed below
    !
    IF ( mask_checksum > 0 ) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
       DO jb = 1, patch_horz%nblks_c
          DO jc = 1, nproma

             IF ( lsmnolake(jc, jb) .LT. 1.0_wp ) THEN
               ! ocean point (fraction of ocean is >0., lsmnolake .lt. 1.) is valid
               is_valid((jb-1)*nproma+jc) = .TRUE.
             ELSE
               ! land point (fraction of land is one, no sea water, lsmnolake=1.) is undef
               is_valid((jb-1)*nproma+jc) = .FALSE.
             ENDIF

          ENDDO
       ENDDO
!ICON_OMP_END_PARALLEL_DO
    ELSE
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
       DO jc = 1, patch_horz%nblks_c * nproma
          is_valid(jc) = .TRUE.
       ENDDO
!ICON_OMP_END_PARALLEL_DO
    ENDIF

    CALL yac_fdef_mask (          &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & is_valid,                 &
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
    field_name(13) = "sea_level_pressure"

    DO jc = 1, no_of_fields
      CALL yac_fdef_field_mask ( &
        & TRIM(field_name(jc)),  &
        & comp_id,               &
        & cell_point_ids,        &
        & cell_mask_ids(1),      &
        & 1,                     &
        & field_id(jc) )
    ENDDO

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)
    !
    ! ! Define a new cell_mask_ids(1) for runoff:
    ! !slo old!   Ocean coastal points with respect to HDmodel mask only are valid.
    ! !slo old!   The integer mask for the HDmodel is ext_data(1)%atm%lsm_hd_c(:,:).
    ! !slo old!   Caution: jg=1 is only valid for coupling to ocean
    ! !
    ! Define cell_mask_ids(1) for runoff - same as above, ocean wet points are valid
    IF ( mask_checksum > 0 ) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
        DO jb = 1, patch_horz%nblks_c
          DO jc = 1, nproma

             IF ( lsmnolake(jc, jb) .LT. 1.0_wp ) THEN
               ! ocean point (fraction of ocean is >0., lsmnolake .lt. 1.) is valid
               is_valid((jb-1)*nproma+jc) = .TRUE.
             ELSE
               ! land point (fraction of land is one, lsmnolake=1.) is undef
               is_valid((jb-1)*nproma+jc) = .FALSE.
             ENDIF

          ENDDO
        ENDDO
!ICON_OMP_END_PARALLEL_DO
    ELSE
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
       DO jc = 1,patch_horz%nblks_c * nproma
          is_valid(jc) = .TRUE.
       ENDDO
!ICON_OMP_END_PARALLEL_DO

    ENDIF

    CALL yac_fdef_mask (          &
      & grid_id,                  &
      & patch_horz%n_patch_cells, &
      & YAC_LOCATION_CELL,        &
      & is_valid,                 &
      & cell_mask_ids(1) )

    ! Define additional coupling field(s) for JSBACH/HD
    ! Utilize mask field for runoff
    ! cell_mask_ids(1) shall contain ocean coast points only for source point mapping (source_to_target_map)
    ! Currently it is the same mask as for the rest. 

    CALL jsb_fdef_hd_fields(comp_id, cell_point_ids, cell_mask_ids )

#endif

    DEALLOCATE (is_valid)

    DEALLOCATE (lsmnolake)

    ! End definition of coupling fields and search

    CALL yac_fsearch ( error_status )

    IF (ltimer) CALL timer_stop(timer_coupling_init)

  END SUBROUTINE construct_atmo_coupling


  !>
  !! SUBROUTINE destruct_echam_ocean_coupling -- terminates the coupling
  !! between ECHAM physics and the ocean.
  !!
  !! This subroutine is called at the end of the time loop of the ICONAM model.

  SUBROUTINE destruct_atmo_coupling

    IF ( .NOT. is_coupled_run() ) RETURN

    CALL yac_ffinalize

  END SUBROUTINE destruct_atmo_coupling
  
END MODULE mo_atmo_coupling_frame

