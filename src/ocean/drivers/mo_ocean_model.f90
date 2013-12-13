!>
!! @page pagecontrolmodelf901 Main program for the ICON ocean model
!!
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
MODULE mo_ocean_model
  
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_master_control,      ONLY: is_restart_run, get_my_process_name, get_my_model_no
  USE mo_parallel_config,     ONLY: p_test_run, l_test_openmp, num_io_procs , num_restart_procs
  USE mo_mpi,                 ONLY: my_process_is_io,set_mpi_work_communicators,p_pe_work, process_mpi_io_size
  USE mo_timer,               ONLY: init_timer, timer_start, timer_stop, print_timer, timer_model_init
  USE mo_datetime,            ONLY: t_datetime
  USE mo_name_list_output_init, ONLY: init_name_list_output, parse_variable_groups
  USE mo_name_list_output,    ONLY: close_name_list_output
  USE mo_name_list_output_config,  ONLY: use_async_name_list_io
  USE mo_grid_config,         ONLY: n_dom
  USE mo_dynamics_config,     ONLY: iequations
  
  !  USE mo_advection_config,    ONLY: configure_advection
  USE mo_dynamics_config,     ONLY: configure_dynamics  ! subroutine
  USE mo_run_config,          ONLY: configure_run, output_mode
  USE mo_gribout_config,      ONLY: configure_gribout
  
  ! Control parameters: run control, dynamics, i/o
  !
  USE mo_run_config,          ONLY: &
    & dtime,                  & !    :
    & nsteps,                 & !    :
    & ltimer,                 & !    :
    & iforcing,               & !
    & num_lev, num_levp1,     &
    & nlev, nlevp1,           &
    & iqc, iqi, iqr, iqs,     &
    & iqni, iqni_nuc, iqg,    &
    & nshift, ntracer,        &
    & grid_generatingcenter,  & ! grid generating center
    & grid_generatingsubcenter  ! grid generating subcenter
  
  USE mo_ocean_nml_crosscheck,   ONLY: oce_crosscheck
  
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d, p_patch_local_parent
  
  ! Horizontal grid
  !
  USE mo_grid_config,         ONLY: n_dom
  
  USE mo_oce_state,           ONLY: t_hydro_ocean_state, setup_ocean_namelists, destruct_patch_3D
  USE mo_build_decomposition, ONLY: build_decomposition
  USE mo_complete_subdivision,ONLY: setup_phys_patches
  
  USE mo_impl_constants,      ONLY: success !, ihs_ocean
  
  ! External data
  USE mo_ocean_ext_data,      ONLY: ext_data, construct_ocean_ext_data, destruct_ocean_ext_data
  
  USE mo_hydro_ocean_run,     ONLY: perform_ho_stepping,&
    & prepare_ho_stepping, &
    & construct_ocean_states, &
    & finalise_ho_integration
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_oce_physics,         ONLY: v_params!, t_ho_params, t_ho_physics
  USE mo_sea_ice_types,       ONLY: t_atmos_fluxes, t_atmos_for_ocean, &
    & v_sfc_flx, v_sea_ice!, t_sfc_flx
  USE mo_impl_constants,      ONLY: max_char_length

  USE mo_alloc_patches,       ONLY: destruct_patches
  USE mo_ocean_read_namelists, ONLY: read_ocean_namelists
  USE mo_io_restart,          ONLY: read_restart_info_file, read_restart_files
  USE mo_io_restart_namelist, ONLY: read_restart_namelists
  USE mo_io_restart_attributes,ONLY: read_restart_attributes, get_restart_attribute
  USE mo_oce_patch_setup,     ONLY: complete_ocean_patch
  USE mo_time_config,         ONLY: time_config
  USE mo_icon_comm_interface, ONLY: construct_icon_communication, destruct_icon_communication
  USE mo_mtime_extensions,    ONLY: get_datetime_string
  USE mo_output_event_types,  ONLY: t_sim_step_info
  USE mtime,                  ONLY: setCalendar, PROLEPTIC_GREGORIAN
  
  !-------------------------------------------------------------
  ! For the coupling
#ifndef __NO_ICON_ATMO__
# ifdef YAC_coupling
  USE mo_parallel_config,     ONLY: nproma
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE finterface_description, ONLY: yac_finit, yac_fdef_comp,                    &
    &                               yac_fdef_subdomain, yac_fconnect_subdomains, &
    &                               yac_fdef_elements, yac_fdef_points,          &
    &                               yac_fdef_mask, yac_fdef_field, yac_fsearch,  &
    &                               yac_ffinalize
  USE mo_coupling_config,     ONLY: is_coupled_run
# else
  USE mo_icon_cpl_init,       ONLY: icon_cpl_init
  USE mo_icon_cpl_init_comp,  ONLY: icon_cpl_init_comp
  USE mo_coupling_config,     ONLY: is_coupled_run, config_debug_coupler_level
  USE mo_icon_cpl_def_grid,   ONLY: icon_cpl_def_grid, icon_cpl_def_location
  USE mo_icon_cpl_def_field,  ONLY: icon_cpl_def_field
  USE mo_icon_cpl_search,     ONLY: icon_cpl_search
  USE mo_icon_cpl_finalize,   ONLY: icon_cpl_finalize
# endif
#endif
  !-------------------------------------------------------------
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: ocean_model
  PUBLIC :: construct_ocean_model, destruct_ocean_model
  PUBLIC :: ocean_patch_3d, ocean_state, operators_coefficients

  TYPE(t_patch_3d), POINTER :: ocean_patch_3d
  TYPE(t_atmos_for_ocean)                         :: p_as
  TYPE(t_atmos_fluxes)                            :: p_atm_f
  TYPE(t_operator_coeff)                          :: operators_coefficients
  TYPE(t_hydro_ocean_state), ALLOCATABLE, TARGET  :: ocean_state(:)
  TYPE(t_datetime)                                :: start_datetime

  
CONTAINS


  !--------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE ocean_model(oce_namelist_filename,shr_namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: oce_namelist_filename,shr_namelist_filename
    
    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:ocean_model"
    
    INTEGER :: jg, ist
    INTEGER :: error_status
    TYPE(t_sim_step_info) :: sim_step_info  
    INTEGER :: jstep0

    !-------------------------------------------------------------------
    IF (is_restart_run()) THEN
      CALL read_restart_ocean_namelists('restart_oce.nc')
    END IF

    !-------------------------------------------------------------------
    CALL construct_ocean_model(oce_namelist_filename,shr_namelist_filename)
    
    !-------------------------------------------------------------------
    IF (is_restart_run()) THEN
      ! This is an resumed integration. Read model state from restart file(s).
#ifdef NOMPI
      CALL read_restart_files
#else
      jg = 1 !no nesting
      !DO jg = ,n_dom
      CALL read_restart_files( ocean_patch_3d%p_patch_2d(jg) )
      !END DO
#endif
      CALL message(TRIM(routine),'normal exit from read_restart_files')
      !ELSE
      !  Prepare the initial conditions:
      !  forcing is part of the restart file
    END IF ! is_restart_run()
    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------

    IF (output_mode%l_nml) THEN
      WRITE(0,*)'process_mpi_io_size:',process_mpi_io_size
      IF (process_mpi_io_size > 0) use_async_name_list_io = .TRUE.
      CALL parse_variable_groups()
      CALL setCalendar(PROLEPTIC_GREGORIAN)
      ! compute sim_start, sim_end
      CALL get_datetime_string(sim_step_info%sim_start, time_config%ini_datetime)
      CALL get_datetime_string(sim_step_info%sim_end,   time_config%end_datetime)
      sim_step_info%dtime      = dtime
      sim_step_info%iadv_rcf   = 1
      sim_step_info%dt_restart = time_config%dt_restart
      jstep0 = 0
      IF (is_restart_run()) THEN
        ! get start counter for time loop from restart file:
        CALL get_restart_attribute("jstep", jstep0)
      END IF
      sim_step_info%jstep0    = jstep0
      CALL init_name_list_output(sim_step_info, opt_lprintlist=.TRUE.,opt_l_is_ocean=.TRUE.)
    ENDIF

    CALL prepare_ho_stepping(ocean_patch_3d,operators_coefficients,ocean_state(1),is_restart_run())

    !------------------------------------------------------------------
    CALL perform_ho_stepping( ocean_patch_3d, ocean_state,                    &
      & ext_data, start_datetime,                     &
      & (nsteps == INT(time_config%dt_restart/dtime)),&
      & v_sfc_flx,                                    &
      & v_params, p_as, p_atm_f,v_sea_ice,operators_coefficients)
    !------------------------------------------------------------------
    
    CALL print_timer()
    
    !------------------------------------------------------------------
    !  cleaning up process
    CALL destruct_ocean_model()

  END SUBROUTINE ocean_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE destruct_ocean_model()

    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:destruct_ocean_model"

    INTEGER :: error_status

    !------------------------------------------------------------------
    !  cleaning up process
    !------------------------------------------------------------------
    CALL message(TRIM(routine),'start to clean up')

    CALL finalise_ho_integration(ocean_state, v_params, &
      & p_as, p_atm_f, v_sea_ice, v_sfc_flx)

    !---------------------------------------------------------------------
    ! 13. Integration finished. Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    ! Destruct external data state
    CALL destruct_ocean_ext_data

    ! deallocate ext_data array
    DEALLOCATE(ext_data, stat=error_status)
    IF (error_status/=success) THEN
      CALL finish(TRIM(routine), 'deallocation of ext_data')
    ENDIF

    !The 3D-ocean version of previous calls
    CALL destruct_patches( ocean_patch_3d%p_patch_2d )
    CALL destruct_patches( p_patch_local_parent )
    NULLIFY( ocean_patch_3d%p_patch_2d )
    CALL destruct_patch_3D( ocean_patch_3d )


    ! Delete variable lists

    IF (output_mode%l_nml) THEN
      CALL close_name_list_output
    ENDIF

    CALL destruct_icon_communication()

#ifndef __NO_ICON_ATMO__
# ifdef YAC_coupling
    IF ( is_coupled_run() ) CALL yac_ffinalize
# else
    IF ( is_coupled_run() ) CALL icon_cpl_finalize ()
# endif
#endif

    CALL message(TRIM(routine),'clean-up finished')

  END SUBROUTINE destruct_ocean_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !!
  !! It does not include the restart processes, these are called from the calling routine ocean_model
  !!
  SUBROUTINE construct_ocean_model(oce_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: oce_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:construct_ocean_model"

    INTEGER :: ist

    INTEGER :: error_status
    !-------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL read_ocean_namelists(oce_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------
    CALL oce_crosscheck()

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs, num_restart_procs)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    CALL init_timer

    IF (ltimer) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! Initialize date and time
    !-------------------------------------------------------------------
    start_datetime = time_config%cur_datetime

    !-------------------------------------------------------------------
    ! 4. Import patches
    !-------------------------------------------------------------------
    CALL build_decomposition(num_lev,num_levp1,nshift, is_ocean_decomposition =.TRUE., &
      &                      patch_3d=ocean_patch_3d)
    CALL construct_icon_communication(ocean_patch_3d%p_patch_2d(:), n_dom=1)
    CALL complete_ocean_patch(ocean_patch_3d%p_patch_2d(1))

    !--------------------------------------------
    ! Setup the information for the physical patches
    CALL setup_phys_patches

    CALL setup_ocean_namelists(ocean_patch_3d%p_patch_2d(1))
    !------------------------------------------------------------------
    ! step 5b: allocate state variables
    !------------------------------------------------------------------
    ALLOCATE (ocean_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for ocean_state failed')
    ENDIF

    !---------------------------------------------------------------------
    ! 9. Horizontal and vertical grid(s) are now defined.
    !    Assign values to derived variables in the configuration states
    !---------------------------------------------------------------------

    CALL configure_dynamics ( n_dom )

    CALL configure_gribout(grid_generatingcenter, grid_generatingsubcenter, n_dom)

    !    DO jg =1,n_dom
    !      !The 3D-ocean version of previous calls
    !      CALL configure_advection( jg, ocean_patch_3d%p_patch_2D(jg)%nlev, ocean_patch_3d%p_patch_2D(1)%nlev, &
    !        &                      iequations, iforcing, iqc, iqi, iqr, iqs, iqni, iqni_nuc, iqg, &
    !        &                      0, 1, .false., .true., ntracer )
    !    ENDDO

    !------------------------------------------------------------------
    ! 10. Create and optionally read external data fields
    !------------------------------------------------------------------
    ALLOCATE (ext_data(n_dom), stat=error_status)
    IF (error_status /= success) THEN
      CALL finish(TRIM(routine),'allocation for ext_data failed')
    ENDIF

    ! allocate memory for oceanic external data and
    ! optionally read those data from netCDF file.
    CALL construct_ocean_ext_data(ocean_patch_3d%p_patch_2d(1:), ext_data)

#ifndef __NO_ICON_ATMO__
    IF ( is_coupled_run() ) THEN
      CALL construct_ocean_coupling()
    ENDIF
#endif

    ! Prepare time integration
    CALL construct_ocean_states(ocean_patch_3d, ocean_state, ext_data, v_sfc_flx, &
      & v_params, p_as, p_atm_f, v_sea_ice,operators_coefficients)!,p_int_state(1:))


  END SUBROUTINE construct_ocean_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
#ifndef __NO_ICON_ATMO__
  !------------------------------------------------------------------
  ! Prepare the coupling
  !
  ! For the time being this could all go into a subroutine which is
  ! common to atmo and ocean. Does this make sense if the setup deviates
  ! too much in future.
  !------------------------------------------------------------------
  SUBROUTINE construct_ocean_coupling()

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

    INTEGER            :: comp_id
    INTEGER            :: cell_point_id
    INTEGER            :: edge_point_id
    INTEGER            :: mask_id
    INTEGER            :: subdomain_id
    INTEGER            :: subdomain_ids(nbr_subdomain_ids)

    INTEGER, PARAMETER :: CELL     = 0 ! one point per cell
    INTEGER, PARAMETER :: CORNER   = 1 ! one point per vertex
    INTEGER, PARAMETER :: EDGE     = 2 ! one point per edge
                                       ! (see definition of enum location in points.h)
    INTEGER            :: jb, jc, je, index
    INTEGER            :: cell_start_idx, cell_end_idx
    INTEGER            :: edge_start_idx, edge_end_idx

    REAL(wp), ALLOCATABLE :: buffer_x(:)
    REAL(wp), ALLOCATABLE :: buffer_y(:)
    REAL(wp), ALLOCATABLE :: buffer_c(:)

    TYPE(t_subset_range), POINTER :: all_cells, all_edges
    TYPE(t_patch), POINTER        :: patch_horz

    patch_no = 1

    ! Initialise the coupler
    CALL yac_finit ( "couling.xml", "coupling.xsd" )

    ! Inform the coupler about what we are
    CALL yac_fdef_comp ( "ICON_ocean", comp_id )

    ! Announce one subdomain (patch) to the coupler
    CALL yac_fdef_subdomain ( comp_id, "ICON_ocean", subdomain_id )

    patch_horz => patch_3D%p_patch_2D(patch_no)
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
          index = (jb-1)*nproma+jc
          buffer_x(index) = patch_horz%verts(jc,jb)%vertex%lon
          buffer_y(index) = patch_horz%verts(jc,jb)%vertex%lat
          buffer_c((index-1)*3+1) = patch_horz%cells%vertex_idx(jc,jb,1)
          buffer_c((index-1)*3+2) = patch_horz%cells%vertex_idx(jc,jb,2)
          buffer_c((index-1)*3+3) = patch_horz%cells%vertex_idx(jc,jb,3)
       ENDDO
    ENDDO

    ! Description of elements, here as unstructured grid
    CALL yac_fdef_elements ( subdomain_id,              &
                             patch_horz%n_patch_verts,  &
                             patch_horz%n_patch_cells,  &
                             nbr_vertices_per_cell,     &
                             buffer_x,                  &
                             buffer_y,                  &
                             buffer_c )

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
          index = (jb-1)*nproma+jc
          buffer_x(index) = patch_horz%cells%center(jc,jb)%lon
          buffer_x(index) = patch_horz%cells%center(jc,jb)%lat
       ENDDO
    ENDDO

    CALL yac_fdef_points ( subdomain_id,            &
                           patch_horz%n_patch_cells,   &
                           CELL,                    &
                           buffer_x,                &
                           buffer_y,                &
                           cell_point_id )

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
          index = (jb-1)*nproma+je
          buffer_x(index) = patch_horz%edges%center(jc,jb)%lon
          buffer_x(index) = patch_horz%edges%center(jc,jb)%lat
       ENDDO
    ENDDO

    CALL yac_fdef_points ( subdomain_id,             &
                           patch_horz%n_patch_cells, &
                           EDGE,                     &
                           buffer_x,                 &
                           buffer_y,                 &
                           edge_point_id )

    ! Connect subdomains
    CALL yac_fconnect_subdomains ( comp_id,           &
                                   nbr_subdomain_ids, &
                                   subdomain_ids,     &
                                   domain_id )
    !
    ! mask generation : ... not yet defined ...
    !
    ! We could use the patch_horz%cells%decomp_info%owner_local information
    ! e.g. to mask out halo points. We do we get the info about what is local and what
    ! is remote.
    !
    ! The land-sea mask for the ocean is available in p_patch_3D%surface_cell_sea_land_mask(:,:)
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
                         imask,      &  !rr TODO
                         points_id,  &
                         mask_id )

    DEALLOCATE (buffer_x, buffer_y, buffer_c)

# else

    !------------------------------------------------------------
    CALL icon_cpl_init(debug_level=config_debug_coupler_level)
    ! Inform the coupler about what we are
    CALL icon_cpl_init_comp ( get_my_process_name(), get_my_model_no(), error_status )
    ! split the global_mpi_communicator into the components
    !------------------------------------------------------------
    patch_no      = 1

    grid_shape(1) = 1
    grid_shape(2) = ocean_patch_3d%p_patch_2d(patch_no)%n_patch_cells

    CALL icon_cpl_def_grid ( &
      & grid_shape, ocean_patch_3d%p_patch_2d(patch_no)%cells%decomp_info%glb_index, & ! input
      & grid_id, error_status )                          ! output

    ! Marker for internal and halo points, a list which contains the
    ! rank where the native cells are located.
    CALL icon_cpl_def_location ( &
      & grid_id, grid_shape, ocean_patch_3d%p_patch_2d(patch_no)%cells%decomp_info%owner_local, & ! input
      & p_pe_work,  & ! this owner id
      & error_status )                                            ! output

# endif

    field_name(1) = "TAUX"   ! bundled field containing two components
    field_name(2) = "TAUY"   ! bundled field containing two components
    field_name(3) = "SFWFLX" ! bundled field containing two components
    field_name(4) = "SFTEMP"
    field_name(5) = "THFLX"  ! bundled field containing two components
    field_name(6) = "ICEATM" ! bundled field containing four components
    field_name(7) = "SST"
    field_name(8) = "OCEANU"
    field_name(9) = "OCEANV"
    field_name(10) = "ICEOCE" ! bundled field containing four components

# ifdef YAC_coupling
       DO i = 1, no_of_fields
          CALL yac_fdef_field ( field_name(i),            &
                                comp_id,                  &
                                domain_id,                &
                                point_id,                 &
                                mask_id,                  &
                                patch_horz%n_patch_cells, &
                                field_id(i) )
       ENDDO

       CALL yac_fsearch ( nbr_components, comp_id, no_of_fields, field_id, error_status )
# else

    field_shape(1:2) = grid_shape(1:2)

    DO i = 1, no_of_fields
      IF ( i == 1 .OR. i == 2 .OR. i == 3 .OR. i == 5 ) THEN
        field_shape(3) = 2
      ELSE IF ( i == 6 ) THEN
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
#endif
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE read_restart_ocean_namelists(restart_filename)

    CHARACTER(LEN=*), INTENT(in) :: restart_filename

    CHARACTER(*), PARAMETER :: routine = "mo_ocean_model:read_restart_namelists"

    CHARACTER(LEN=max_char_length) :: grid_file_name
    LOGICAL :: lsuccess
    !-------------------------------------------------------------------

    ! First read the restart master file (ASCII format) to find out
    ! which NetCDF files the model should read.
    ! Comment by Hui Wan:
    ! The namelist variable atmo_restart_info_filename should be
    ! an input argument of the subroutine read_restart_info_file.
    CALL read_restart_info_file(grid_file_name, lsuccess) ! out, out

    IF (lsuccess) THEN
      CALL message( TRIM(routine),                          &
        & 'Running model in restart mode. '       &
        & //'Horizontal grid should be read from '&
        & //TRIM(grid_file_name) )
    ELSE
      CALL finish(TRIM(routine),'Failed to read restart info file')
    END IF

    ! Read all namelists used in the previous run
    ! and store them in a buffer. These values will overwrite the
    ! model default, and will later be overwritten if the user has
    ! specified something different for this integraion.
    CALL read_restart_namelists(restart_filename)
    CALL message(TRIM(routine), 'read namelists from ocean restart file')

    ! Read all global attributs in the restart file and store them in a buffer.

    CALL read_restart_attributes(restart_filename)
    CALL message(TRIM(routine), 'read global attributes from ocean restart file')

  END SUBROUTINE read_restart_ocean_namelists
  !--------------------------------------------------------------------------

  
END MODULE mo_ocean_model

