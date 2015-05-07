!>
!! Module handling the initialization of synchronous and asynchronous output.
!!
!! @author R. Johanni
!!
!! @par Revision History
!! Initial implementation  by  R. Johanni  (2011)
!! Major changes: F. Prill, DWD (2012-2013)
!! A-priori calculation of output step events: F. Prill, DWD (10/2013)
!!
!! Define USE_CRAY_POINTER for platforms having problems with ISO_C_BINDING
!! BUT understand CRAY pointers
!!
!!   #define USE_CRAY_POINTER
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_name_list_output_init

#ifndef USE_CRAY_POINTER
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_intptr_t, c_f_pointer, c_int64_t
#endif
! USE_CRAY_POINTER

  ! constants and global settings
  USE mo_cdi_constants          ! We need all
  USE mo_kind,                              ONLY: wp, i8, dp, sp
  USE mo_impl_constants,                    ONLY: max_phys_dom, max_dom, SUCCESS,                 &
    &                                             max_var_ml, max_var_pl, max_var_hl, max_var_il, &
    &                                             MAX_TIME_LEVELS, vname_len,                     &
    &                                             MAX_CHAR_LENGTH, MAX_NUM_IO_PROCS,              &
    &                                             MAX_TIME_INTERVALS, ihs_ocean, MAX_NPLEVS,      &
    &                                             MAX_NZLEVS, MAX_NILEVS
  USE mo_io_units,                          ONLY: filename_max, nnml, nnml_output
  USE mo_master_nml,                        ONLY: model_base_dir
  USE mo_master_control,                    ONLY: is_restart_run, my_process_is_ocean
  ! basic utility modules
  USE mo_exception,                         ONLY: finish, message, message_text
  USE mo_dictionary,                        ONLY: t_dictionary, dict_init,                        &
    &                                             dict_loadfile, dict_get, DICT_MAX_STRLEN
  USE mo_fortran_tools,                     ONLY: assign_if_present
  USE mo_util_cdi,                          ONLY: set_additional_GRIB2_keys
  USE mo_util_uuid,                         ONLY: uuid2char
  USE mo_io_util,                           ONLY: get_file_extension
  USE mo_util_string,                       ONLY: t_keyword_list, associate_keyword,              &
    &                                             with_keywords, insert_group,                    &
    &                                             tolower, int2string, difference,                &
    &                                             sort_and_compress_list, one_of
  USE mo_datetime,                          ONLY: t_datetime
  USE mo_cf_convention,                     ONLY: t_cf_var
  USE mo_io_restart_attributes,             ONLY: get_restart_attribute
  USE mo_model_domain,                      ONLY: p_patch, p_phys_patch
  USE mo_mtime_extensions,                  ONLY: get_datetime_string, get_duration_string, &
                                                  get_duration_string_real
  USE mo_math_utilities,                    ONLY: merge_values_into_set
  ! config modules
  USE mo_parallel_config,                   ONLY: nproma, p_test_run, use_dp_mpi2io

  USE mo_run_config,                        ONLY: dtime, msg_level, output_mode,                  &
    &                                             number_of_grid_used
  USE mo_grid_config,                       ONLY: n_dom, n_phys_dom, start_time, end_time,        &
    &                                             DEFAULT_ENDTIME
  USE mo_io_config,                         ONLY: netcdf_dict, output_nml_dict
  USE mo_name_list_output_config,           ONLY: use_async_name_list_io,                         &
    &                                             first_output_name_list,                         &
    &                                             add_var_desc
  USE mo_time_config,                       ONLY: time_config
  USE mo_gribout_config,                    ONLY: gribout_config
  USE mo_dynamics_config,                   ONLY: iequations

#ifndef __NO_ICON_ATMO__
  USE mo_nh_pzlev_config,                   ONLY: nh_pzlev_config
#endif

  ! MPI Communication routines
  USE mo_mpi,                               ONLY: p_bcast, get_my_mpi_work_id, p_max,             &
    &                                             get_my_mpi_work_communicator,                   &
    &                                             p_comm_work, p_comm_work_2_io,                  &
    &                                             p_comm_io, p_comm_work_io,                      &
    &                                             p_int, p_real_dp, p_real_sp,                    &
    &                                             my_process_is_stdio, my_process_is_mpi_test,    &
    &                                             my_process_is_mpi_workroot,                     &
    &                                             my_process_is_mpi_seq, my_process_is_io,        &
    &                                             my_process_is_mpi_ioroot,                       &
    &                                             process_mpi_stdio_id, process_work_io0,         &
    &                                             process_mpi_io_size, num_work_procs, p_n_work,  &
    &                                             p_pe_work, p_io_pe0, p_pe
  USE mo_communication,                     ONLY: idx_no, blk_no
  ! namelist handling
  USE mo_namelist,                          ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_nml_annotate,                      ONLY: temp_defaults, temp_settings
  ! variable lists
  USE mo_var_metadata_types,                ONLY: t_var_metadata, VARNAME_LEN
  USE mo_linked_list,                       ONLY: t_var_list, t_list_element
  USE mo_var_list,                          ONLY: nvar_lists, max_var_lists, var_lists,           &
    &                                             new_var_list,                                   &
    &                                             total_number_of_variables, collect_group,       &
    &                                             get_var_timelevel, get_var_name,                &
    &                                             get_var_tileidx
  USE mo_var_list_element,                  ONLY: level_type_ml, level_type_pl, level_type_hl,    &
    &                                             level_type_il
  ! lon-lat interpolation
  USE mo_lonlat_grid,                       ONLY: t_lon_lat_grid, compute_lonlat_blocking,        &
    &                                             compute_lonlat_specs, threshold_delta_or_intvls
  USE mo_intp_data_strc,                    ONLY: t_lon_lat_intp,                                 &
    &                                             t_lon_lat_data, get_free_lonlat_grid,           &
    &                                             lonlat_grid_list, n_lonlat_grids,               &
    &                                             get_lonlat_grid_ID
  ! output events
  USE mtime,                                ONLY: MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN,    &
    &                                             timedelta, newTimedelta, deallocateTimedelta,   &
    &                                             OPERATOR(<), newDatetime, deallocateDatetime,   &
    &                                             getTotalMilliSecondsTimeDelta, datetime,        &
    &                                             OPERATOR(+), datetimeToString, OPERATOR(>)
  USE mo_mtime_extensions,                  ONLY: get_datetime_string, get_duration_string
  USE mo_output_event_types,                ONLY: t_sim_step_info, MAX_EVENT_NAME_STR_LEN,        &
    &                                             DEFAULT_EVENT_NAME, t_par_output_event
  USE mo_output_event_control,              ONLY: compute_matching_sim_steps,                     &
    &                                             generate_output_filenames
  USE mo_output_event_handler,              ONLY: new_parallel_output_event,                      &
    &                                             complete_event_setup, union_of_all_events,      &
    &                                             print_output_event, trigger_output_step_irecv,  &
    &                                             set_event_to_simstep
  ! name list output
  USE mo_name_list_output_types,            ONLY: l_output_phys_patch, t_output_name_list,        &
    &                                             t_output_file, t_var_desc,                      &
    &                                             t_patch_info, t_reorder_info,                   &
    &                                             REMAP_NONE, REMAP_REGULAR_LATLON,               &
    &                                             sfs_name_list, ffs_name_list, second_tos,       &
    &                                             first_tos, GRP_PREFIX, t_fname_metadata,        &
    &                                             all_events, t_patch_info_ll,GRB2_GRID_INFO,     &
    &                                             is_grid_info_var, GRB2_GRID_INFO_NAME
  USE mo_name_list_output_gridinfo,         ONLY: set_grid_info_grb2,                             &
    &                                             set_grid_info_netcdf, collect_all_grid_info,    &
    &                                             copy_grid_info, bcast_grid_info,                &
    &                                             deallocate_all_grid_info,                       &
    &                                             GRID_INFO_NONE, GRID_INFO_FILE, GRID_INFO_BCAST
  USE mo_name_list_output_metadata,         ONLY: metainfo_allocate_memory_window
  USE mo_name_list_output_zaxes,            ONLY: setup_ml_axes_atmo, setup_pl_axis_atmo,         &
    &                                             setup_hl_axis_atmo, setup_il_axis_atmo,         &
    &                                             setup_zaxes_oce
  USE mo_util_vgrid_types,                  ONLY: vgrid_buffer

#ifndef __NO_ICON_ATMO__
  USE mo_vertical_coord_table,              ONLY: vct
#endif

#if !defined (__NO_ICON_ATMO__) && !defined (__NO_ICON_OCEAN__)
  USE mo_coupling_config,                   ONLY: is_coupled_run
  USE mo_master_control,                    ONLY: get_my_process_name
#endif

  IMPLICIT NONE

  PRIVATE

  ! variables and data types
  PUBLIC :: out_varnames_dict
  PUBLIC :: varnames_dict
  PUBLIC :: output_file
  PUBLIC :: patch_info
  PUBLIC :: lonlat_info
  ! subroutines
  PUBLIC :: read_name_list_output_namelists
  PUBLIC :: parse_variable_groups
  PUBLIC :: init_name_list_output
  PUBLIC :: setup_output_vlist
  PUBLIC :: collect_requested_ipz_levels
#ifdef USE_CRAY_POINTER
  PUBLIC :: set_mem_ptr_sp
  PUBLIC :: set_mem_ptr_dp
#endif


  !------------------------------------------------------------------------------------------------

  TYPE(t_output_file),   ALLOCATABLE, TARGET :: output_file(:)
  TYPE(t_patch_info),    ALLOCATABLE, TARGET :: patch_info (:)
  TYPE(t_patch_info_ll), ALLOCATABLE, TARGET :: lonlat_info(:,:)

  ! Number of output domains. This depends on l_output_phys_patch and is either the number
  ! of physical or the number of logical domains.
  INTEGER :: n_dom_out

  ! Broadcast root for intercommunicator broadcasts form compute PEs to IO PEs using p_comm_work_2_io
  INTEGER :: bcast_root

  !------------------------------------------------------------------------------------------------
  ! dictionaries for variable names:
  TYPE (t_dictionary) ::     &
    &   varnames_dict      , & !< maps variable names onto the internal ICON names.
    &   out_varnames_dict      !< maps internal variable names onto names in output file (NetCDF only).
  !------------------------------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_init'


CONTAINS

  !------------------------------------------------------------------------------------------------
  !> Read configuration for namelist controlled output module.
  !
  !  The name is a bit strange, but if follows the convention to read
  !  namelists with a routine called read_XXX_namelist (plural is
  !  used here since several namelists are read).
  !
  !  Please note the difference between name_list and namelist!
  !
  SUBROUTINE read_name_list_output_namelists( filename )
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    ! local variables
    CHARACTER(LEN=*), PARAMETER       :: routine = 'read_name_list_output_namelists'

    INTEGER                               :: istat, i, j
    TYPE(t_output_name_list), POINTER     :: p_onl
    INTEGER                               :: nnamelists
    LOGICAL                               :: lrewind

    ! Local variables corresponding to members of output_name_list
    INTEGER                               :: filetype
    INTEGER                               :: mode
    INTEGER                               :: taxis_tunit
    INTEGER                               :: dom(max_phys_dom)
    INTEGER                               :: steps_per_file
    LOGICAL                               :: steps_per_file_inclfirst         !< Flag. Do not count first step in files count
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)  :: file_interval                    !< length of a file (ISO8601 duration)
    LOGICAL                               :: include_last
    LOGICAL                               :: output_grid
    CHARACTER(LEN=filename_max)           :: output_filename
    CHARACTER(LEN=filename_max)           :: filename_format
    CHARACTER(LEN=filename_max)           :: filename_extn                    !< user-specified filename extension (or "default")
    CHARACTER(LEN=vname_len)              :: ml_varlist(max_var_ml)
    CHARACTER(LEN=vname_len)              :: pl_varlist(max_var_pl)
    CHARACTER(LEN=vname_len)              :: hl_varlist(max_var_hl)
    CHARACTER(LEN=vname_len)              :: il_varlist(max_var_il)
    CHARACTER(len=MAX_CHAR_LENGTH)        :: m_levels                         !< level selection: model levels
    REAL(wp)                              :: p_levels(MAX_NPLEVS)             !< pressure levels
    REAL(wp)                              :: h_levels(MAX_NZLEVS)             !< height levels
    REAL(wp)                              :: i_levels(MAX_NILEVS)             !< isentropic levels
    INTEGER                               :: remap
    REAL(wp)                              :: reg_lon_def(3)
    REAL(wp)                              :: reg_lat_def(3)
    INTEGER                               :: reg_def_mode
    REAL(wp)                              :: north_pole(2)

    REAL(wp)                              :: output_bounds(3*MAX_TIME_INTERVALS)
    INTEGER                               :: output_time_unit
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: output_start(MAX_TIME_INTERVALS), &
      &                                      output_end(MAX_TIME_INTERVALS),   &
      &                                      output_interval(MAX_TIME_INTERVALS)
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN) :: ready_file  !< ready filename prefix (=output event name)

    TYPE(t_lon_lat_data),  POINTER        :: lonlat
    TYPE (t_keyword_list), POINTER        :: keywords => NULL()
    CHARACTER(len=MAX_CHAR_LENGTH)        :: cfilename
    INTEGER                               :: iunit, lonlat_id, jg
    TYPE (t_lon_lat_grid)                 :: new_grid

    !> "stream_partitions": Split one namelist into concurrent,
    !> alternating files:
    INTEGER                               :: stream_partitions_ml, &
      &                                      stream_partitions_pl, &
      &                                      stream_partitions_hl, &
      &                                      stream_partitions_il

    !> MPI ranks which were explicitly specified by the user:
    INTEGER                               :: pe_placement_ml(MAX_NUM_IO_PROCS), &
      &                                      pe_placement_pl(MAX_NUM_IO_PROCS), &
      &                                      pe_placement_hl(MAX_NUM_IO_PROCS), &
      &                                      pe_placement_il(MAX_NUM_IO_PROCS)

    !> RBF shape parameter.
    REAL(wp)                              :: rbf_scale

    ! The namelist containing all variables above
    NAMELIST /output_nml/ &
      mode, taxis_tunit, dom,                                &
      filetype, filename_format, output_filename,            &
      steps_per_file, steps_per_file_inclfirst,              &
      include_last, output_grid,                             &
      remap, reg_lon_def, reg_lat_def, north_pole,           &
      ml_varlist, pl_varlist, hl_varlist, il_varlist,        &
      m_levels, p_levels, h_levels, i_levels,                &
      output_start, output_end, output_interval,             &
      output_bounds, output_time_unit,                       &
      ready_file, file_interval, reg_def_mode,               &
      stream_partitions_ml, stream_partitions_pl,            &
      stream_partitions_hl, stream_partitions_il,            &
      pe_placement_ml, pe_placement_pl,                      &
      pe_placement_hl, pe_placement_il,                      &
      filename_extn, rbf_scale

    ! -- preliminary checks:
    !
    ! We need dtime
    IF(dtime<=0._wp) CALL finish(routine, 'dtime must be set before reading output namelists')

    ! Before we start: prepare the levels set objects for the vertical
    ! interpolation.
    !
    ! level ordering: zlevels, plevels, ilevels must be ordered from
    ! TOA to bottom:
    !
#ifndef __NO_ICON_ATMO__
    DO jg = 1, max_dom
      ! pressure levels:
      nh_pzlev_config(jg)%plevels%sort_smallest_first = .TRUE.
      ! height levels
      nh_pzlev_config(jg)%zlevels%sort_smallest_first = .FALSE.
      ! isentropic levels
      nh_pzlev_config(jg)%ilevels%sort_smallest_first = .FALSE.
    END DO
#endif

    ! -- Open input file and position to first namelist 'output_nml'

    CALL open_nml(TRIM(filename))

    ! As in COSMO, there may exist several output_nml namelists in the input file
    ! Loop until EOF is reached

    p_onl                  => NULL()
    first_output_name_list => NULL()
    nnamelists             =  0
    lrewind                = .TRUE.

    IF (.NOT. output_mode%l_nml) RETURN ! do not read output namelists if main switch is set to false

    DO
      CALL position_nml ('output_nml', lrewind=lrewind, status=istat)
      IF(istat /= POSITIONED) THEN

        ! if no "output_nml" has been found at all, we disable this
        ! mode (i.e. the user's namelist settings were inconsistent).
        IF (.NOT.ASSOCIATED(first_output_name_list))  output_mode%l_nml = .FALSE.

        CALL close_nml
        RETURN
      ENDIF
      lrewind = .FALSE.

      ! -- Set all variables in output_nml to their default values

      filetype                 = FILETYPE_NC2 ! NetCDF
      mode                     = 2
      taxis_tunit              = TUNIT_MINUTE
      dom(:)                   = -1
      steps_per_file           = -1
      steps_per_file_inclfirst = .TRUE.
      file_interval            = ' '
      include_last             = .TRUE.
      output_grid              = .FALSE.
      output_filename          = ' '
      filename_format          = "<output_filename>_DOM<physdom>_<levtype>_<jfile>"
      filename_extn            = "default"
      ml_varlist(:)            = ' '
      pl_varlist(:)            = ' '
      hl_varlist(:)            = ' '
      il_varlist(:)            = ' '
      m_levels                 = " "
      p_levels(:)              = -1._wp
      h_levels(:)              = -1._wp
      i_levels(:)              = -1._wp
      remap                    = REMAP_NONE
      reg_lon_def(:)           = 0._wp
      reg_lat_def(:)           = 0._wp
      reg_def_mode             = 0
      north_pole(:)            = (/ 0._wp, 90._wp /)
      output_start(:)          = ' '
      output_end(:)            = ' '
      output_interval(:)       = ' '
      output_bounds(:)         = -1._wp
      output_time_unit         = 1
      ready_file               = DEFAULT_EVENT_NAME
      lonlat_id                = 0
      stream_partitions_ml     = 1
      stream_partitions_pl     = 1
      stream_partitions_hl     = 1
      stream_partitions_il     = 1
      pe_placement_ml(:)       = -1 !< i.e. MPI rank undefined (round-robin placement)
      pe_placement_pl(:)       = -1 !< i.e. MPI rank undefined (round-robin placement)
      pe_placement_hl(:)       = -1 !< i.e. MPI rank undefined (round-robin placement)
      pe_placement_il(:)       = -1 !< i.e. MPI rank undefined (round-robin placement)
      rbf_scale                = -1._wp

      ! -- Read output_nml

      IF (my_process_is_stdio())  THEN
        iunit = temp_defaults()
        WRITE(iunit, output_nml)                                     ! write defaults to temporary text file
      END IF
      READ (nnml, output_nml)                          ! overwrite default settings
      IF (my_process_is_stdio())  THEN
        iunit = temp_settings()
        WRITE(iunit, output_nml)                                     ! write settings to temporary text file
      END IF

      nnamelists = nnamelists+1

      ! -- Consistency checks:

      IF ((steps_per_file == -1) .AND. (TRIM(file_interval) == "")) THEN
        CALL finish(routine, "Please specify either <steps_per_file> or <file_interval>!")
      END IF
      IF ((steps_per_file /= -1) .AND. (TRIM(file_interval) /= "")) THEN
        CALL finish(routine, "User has specified conflicting parameters <steps_per_file>, <file_interval>!")
      END IF
      IF(remap/=REMAP_NONE .AND. remap/=REMAP_REGULAR_LATLON) THEN
        CALL finish(routine,'Unsupported value for remap')
      END IF
      IF ((remap==REMAP_REGULAR_LATLON) .AND. ALL(reg_lon_def(:) == 0.) .AND. ALL(reg_lat_def(:) == 0.)) THEN
        CALL finish(routine,'Lon-lat output: Grid not specified in namelist!')
      END IF
      IF ((reg_lon_def(3) >  reg_lon_def(1)) .AND. &
        & (reg_lon_def(2) <= 0._wp)) THEN
        CALL finish(routine,'Illegal LON increment')
      END IF
      IF ((reg_lat_def(3) /= reg_lat_def(1)) .AND. &
        & (reg_lat_def(2) == 0._wp)) THEN
        CALL finish(routine,'Illegal LAT increment')
      END IF
      IF(reg_lon_def(3)<reg_lon_def(1)) CALL finish(routine,'end lon < start lon')

      ! -- Scale output bounds

      IF (output_bounds(1) > -1._wp) THEN
        ! output_bounds is always in seconds - the question is what to do with months or years
        ! output_time_unit: 1 = second, 2=minute, 3=hour, 4=day, 5=month, 6=year
        SELECT CASE(output_time_unit)
        CASE(1); output_bounds(:) = output_bounds(:)
        CASE(2); output_bounds(:) = output_bounds(:)*60._wp
        CASE(3); output_bounds(:) = output_bounds(:)*3600._wp
        CASE(4); output_bounds(:) = output_bounds(:)*86400._wp
        CASE(5); output_bounds(:) = output_bounds(:)*86400._wp*30._wp  ! Not a real calender month
        CASE(6); output_bounds(:) = output_bounds(:)*86400._wp*365._wp ! Not a real calender year
        CASE DEFAULT
          CALL finish(routine,'Illegal output_time_unit')
        END SELECT
      END IF

      ! -- Read the map files into dictionary data structures

      CALL dict_init(varnames_dict,     lcase_sensitive=.FALSE.)
      CALL dict_init(out_varnames_dict, lcase_sensitive=.FALSE.)

      CALL associate_keyword("<path>", TRIM(model_base_dir), keywords)
      IF(output_nml_dict     /= ' ') THEN
        cfilename = TRIM(with_keywords(keywords, output_nml_dict))
        CALL message(routine, "load dictionary file.")
        CALL dict_loadfile(varnames_dict, cfilename)
      END IF
      IF(netcdf_dict /= ' ') THEN
        cfilename = TRIM(with_keywords(keywords, netcdf_dict))
        CALL message(routine, "load dictionary file (output names).")
        CALL dict_loadfile(out_varnames_dict, cfilename, linverse=.TRUE.)
      END IF

      ! -- If "remap=1": lon-lat interpolation requested

#ifndef __NO_ICON_ATMO__
      IF (remap == REMAP_REGULAR_LATLON) THEN
        ! define lon-lat grid specification
        new_grid%reg_lon_def(:) = reg_lon_def(:)
        new_grid%reg_lat_def(:) = reg_lat_def(:)
        new_grid%north_pole(:)  = north_pole(:)
        new_grid%reg_def_mode   = reg_def_mode
        ! compute some additional entries of lon-lat grid specification:
        CALL compute_lonlat_specs(new_grid)
        CALL compute_lonlat_blocking(new_grid, nproma)
        ! If the user has explicitly specified an interpolation
        ! parameter, then we always register this as a new lon-lat
        ! grid. Otherwise we might share the lon-lat coefficients with
        ! other output namelists.
        IF (rbf_scale > 0._wp) THEN
          lonlat_id             =  get_free_lonlat_grid()
          lonlat                => lonlat_grid_list(lonlat_id)
          lonlat%grid           =  new_grid
          lonlat%intp%rbf_scale =  rbf_scale
        ELSE
          ! check, if lon-lat grids has already been registered
          lonlat_id = get_lonlat_grid_ID(new_grid)
          IF (lonlat_id == -1) THEN
            ! Register a lon-lat grid data structure in global list
            lonlat_id             =  get_free_lonlat_grid()
            lonlat                => lonlat_grid_list(lonlat_id)
            lonlat%grid           =  new_grid
            lonlat%intp%rbf_scale =  rbf_scale
          ELSE
            lonlat => lonlat_grid_list(lonlat_id)
          END IF
        END IF

        ! Flag those domains, which are used for this lon-lat grid:
        !     If dom(:) was not specified in namelist input, it is set
        !     completely to -1.  In this case all domains are wanted in
        !     the output
        IF (dom(1) < 0)  THEN
          lonlat%l_dom(:) = .TRUE.
        ELSE
          DOM_LOOP : DO i = 1, max_dom
            IF (dom(i) < 0) exit DOM_LOOP
            lonlat%l_dom( dom(i) ) = .TRUE.
          ENDDO DOM_LOOP
        END IF
      ENDIF
#endif
! #ifndef __NO_ICON_ATMO__

      ! -- Allocate next output_name_list

      IF(.NOT.ASSOCIATED(first_output_name_list)) THEN
        ! Allocate first name_list
        ALLOCATE(first_output_name_list)
        p_onl => first_output_name_list
      ELSE
        ! This is not the first one, p_onl points to the last one which was created
        ALLOCATE(p_onl%next)
        p_onl => p_onl%next
      ENDIF

      ! -- Set next output_name_list from values read

      p_onl%filetype                 = filetype
      p_onl%mode                     = mode
      p_onl%taxis_tunit              = taxis_tunit
      p_onl%dom(:)                   = dom(:)
      p_onl%steps_per_file           = steps_per_file
      ! conditional default: steps_per_file_inclfirst=.FALSE. for GRIB output
      p_onl%steps_per_file_inclfirst = steps_per_file_inclfirst .AND. (filetype /= FILETYPE_GRB2)
      p_onl%file_interval            = file_interval
      p_onl%include_last             = include_last
      p_onl%output_grid              = output_grid
      p_onl%output_filename          = output_filename
      p_onl%filename_format          = filename_format
      p_onl%filename_extn            = filename_extn
      p_onl%ml_varlist(:)            = ml_varlist(:)
      p_onl%pl_varlist(:)            = pl_varlist(:)
      p_onl%hl_varlist(:)            = hl_varlist(:)
      p_onl%il_varlist(:)            = il_varlist(:)
      p_onl%m_levels                 = m_levels
      p_onl%p_levels                 = p_levels
      p_onl%z_levels                 = h_levels
      p_onl%i_levels                 = i_levels
      p_onl%remap                    = remap
      p_onl%lonlat_id                = -1
      p_onl%output_start(:)          = output_start(:)
      p_onl%output_end(:)            = output_end
      p_onl%output_interval(:)       = output_interval
      p_onl%additional_days(:)       = 0
      p_onl%output_bounds(:)         = output_bounds(:)
      p_onl%ready_file               = ready_file
      p_onl%lonlat_id                = lonlat_id
      p_onl%stream_partitions_ml     = stream_partitions_ml
      p_onl%stream_partitions_pl     = stream_partitions_pl
      p_onl%stream_partitions_hl     = stream_partitions_hl
      p_onl%stream_partitions_il     = stream_partitions_il
      p_onl%pe_placement_ml(:)       = pe_placement_ml(:)
      p_onl%pe_placement_pl(:)       = pe_placement_pl(:)
      p_onl%pe_placement_hl(:)       = pe_placement_hl(:)
      p_onl%pe_placement_il(:)       = pe_placement_il(:)

      ! -- translate variables names according to variable name
      !    dictionary:
      DO i=1,max_var_ml
        p_onl%ml_varlist(i) = dict_get(varnames_dict, p_onl%ml_varlist(i), &
          &                            default=p_onl%ml_varlist(i))
      END DO
      DO i=1,max_var_pl
        p_onl%pl_varlist(i) = dict_get(varnames_dict, p_onl%pl_varlist(i), &
          &                            default=p_onl%pl_varlist(i))
      END DO
      DO i=1,max_var_hl
        p_onl%hl_varlist(i) = dict_get(varnames_dict, p_onl%hl_varlist(i), &
          &                            default=p_onl%hl_varlist(i))
      END DO
      DO i=1,max_var_il
        p_onl%il_varlist(i) = dict_get(varnames_dict, p_onl%il_varlist(i), &
          &                            default=p_onl%il_varlist(i))
      END DO

      ! allow case-insensitive variable names:
      DO i=1,max_var_ml
        p_onl%ml_varlist(i) = tolower(p_onl%ml_varlist(i))
      END DO
      DO i=1,max_var_pl
        p_onl%pl_varlist(i) = tolower(p_onl%pl_varlist(i))
      END DO
      DO i=1,max_var_hl
        p_onl%hl_varlist(i) = tolower(p_onl%hl_varlist(i))
      END DO
      DO i=1,max_var_il
        p_onl%il_varlist(i) = tolower(p_onl%il_varlist(i))
      END DO
      
      p_onl%next => NULL()

      ! -- if the namelist switch "output_grid" has been enabled: add
      !    "clon, "clat", "elon", "elat", etc. to the list of
      !    variables:
      !
      IF (p_onl%output_grid) THEN
        ! model levels
        IF (TRIM(p_onl%ml_varlist(1)) /=  "") THEN
          SELECT CASE(p_onl%remap)
          CASE (REMAP_NONE)
            DO i=1,3
              DO j=1,2
                CALL append_varname(p_onl%ml_varlist, GRB2_GRID_INFO_NAME(i,j))
              END DO
            END DO
          CASE (REMAP_REGULAR_LATLON)
            DO j=1,2
              CALL append_varname(p_onl%ml_varlist, GRB2_GRID_INFO_NAME(0,j))
            END DO
          END SELECT
        END IF
        ! pressure levels
        IF (TRIM(p_onl%pl_varlist(1)) /=  "") THEN
          SELECT CASE(p_onl%remap)
          CASE (REMAP_NONE)
            DO i=1,3
              DO j=1,2
                CALL append_varname(p_onl%pl_varlist, GRB2_GRID_INFO_NAME(i,j))
              END DO
            END DO
          CASE (REMAP_REGULAR_LATLON)
            DO j=1,2
              CALL append_varname(p_onl%pl_varlist, GRB2_GRID_INFO_NAME(0,j))
            END DO
          END SELECT
        END IF
        ! height levels
        IF (TRIM(p_onl%hl_varlist(1)) /=  "") THEN
          SELECT CASE(p_onl%remap)
          CASE (REMAP_NONE)
            DO i=1,3
              DO j=1,2
                CALL append_varname(p_onl%hl_varlist, GRB2_GRID_INFO_NAME(i,j))
              END DO
            END DO
          CASE (REMAP_REGULAR_LATLON)
            DO j=1,2
              CALL append_varname(p_onl%hl_varlist, GRB2_GRID_INFO_NAME(0,j))
            END DO
          END SELECT
        END IF
        ! isentropic levels
        IF (TRIM(p_onl%il_varlist(1)) /=  "") THEN
          SELECT CASE(p_onl%remap)
          CASE (REMAP_NONE)
            DO i=1,3
              DO j=1,2
                CALL append_varname(p_onl%il_varlist, GRB2_GRID_INFO_NAME(i,j))
              END DO
            END DO
          CASE (REMAP_REGULAR_LATLON)
            DO j=1,2
              CALL append_varname(p_onl%il_varlist, GRB2_GRID_INFO_NAME(0,j))
            END DO
          END SELECT
        END IF
      END IF

      ! -- write the contents of the namelist to an ASCII file

      IF(my_process_is_stdio()) WRITE(nnml_output,nml=output_nml)

    ENDDO

    CALL close_nml

  END SUBROUTINE read_name_list_output_namelists


  !------------------------------------------------------------------------------------------------
  !> Utility routine: searches for the end of a list of variable name
  !  and appends another entry.
  SUBROUTINE append_varname(p_varlist, new_varname)
    CHARACTER(LEN=vname_len), INTENT(INOUT) :: p_varlist(:)
    CHARACTER(len=*),         INTENT(IN)    :: new_varname
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::append_varname"
    INTEGER :: ivar

    ! Get the number of variables in varlist
    DO ivar = 1, SIZE(p_varlist)
      IF (p_varlist(ivar) == ' ') EXIT ! Last one reached
    ENDDO
    IF (ivar > SIZE(p_varlist)) THEN
      CALL finish(routine, "Insufficient array size!")
    END IF
    p_varlist(ivar) = TRIM(tolower(TRIM(new_varname)))
  END SUBROUTINE append_varname


  !------------------------------------------------------------------------------------------------
  !> Appends the chosen p-levels, z-levels, i-levels to the levels
  !  sets for the corresponding domains (note that we do the vertical
  !  interpolation for the union of all chosen levels and only do a
  !  selection for each output namelist):
  !
  SUBROUTINE collect_requested_ipz_levels()
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::collect_requested_ipz_levels"
    !
    TYPE(t_output_name_list), POINTER     :: p_onl
    INTEGER :: n_dom_out, jp, log_patch_id, nlevs

#ifndef __NO_ICON_ATMO__
    ! Loop over the output namelists and create a union set of
    ! all requested vertical levels (per domain):
    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT

      n_dom_out = n_dom
      IF(l_output_phys_patch)  n_dom_out = n_phys_dom

      DO jp = 1, n_dom_out
        ! append the chosen p-levels, z-levels, i-levels

        log_patch_id = jp
        IF(l_output_phys_patch)  log_patch_id = p_phys_patch(jp)%logical_id

        IF ((p_onl%dom(1) <= 0) .OR. (ANY(p_onl%dom(:) == log_patch_id))) THEN

          ! pressure levels
          !
          ! count the no. of levels
          DO nlevs=1,SIZE(p_onl%p_levels)
            IF (p_onl%p_levels(nlevs) < 0._wp) EXIT
          END DO
          nlevs = nlevs - 1
          IF ((nlevs == 0) .AND. (p_onl%pl_varlist(1) /= ' ')) THEN
            CALL finish(routine, "Input error (output_nml): User has requested pressure interpolation without specifying levels!")
          END IF
          ! append nlevs pressure levels (domain log_patch_id)
          IF (nlevs > 0)  CALL merge_values_into_set(nlevs, p_onl%p_levels, &
            &                                        nh_pzlev_config(log_patch_id)%plevels)

          ! height levels
          !
          ! count the no. of levels
          DO nlevs=1,SIZE(p_onl%z_levels)
            IF (p_onl%z_levels(nlevs) < 0._wp) EXIT
          END DO
          nlevs = nlevs - 1
          IF ((nlevs == 0) .AND. (p_onl%hl_varlist(1) /= ' ')) THEN
            CALL finish(routine, "Input error (output_nml): User has requested height interpolation without specifying levels!")
          END IF
          ! append nlevs height levels
          IF (nlevs > 0)  CALL merge_values_into_set(nlevs, p_onl%z_levels, &
            &                                        nh_pzlev_config(log_patch_id)%zlevels)

          ! isentropic levels
          !
          ! count the no. of levels
          DO nlevs=1,SIZE(p_onl%i_levels)
            IF (p_onl%i_levels(nlevs) < 0._wp) EXIT
          END DO
          nlevs = nlevs - 1
          IF ((nlevs == 0) .AND. (p_onl%il_varlist(1) /= ' ')) THEN
            CALL finish(routine, "Input error (output_nml): User has requested isentropic interpolation without specifying levels!")
          END IF
          ! append nlevs isentropic levels
          IF (nlevs > 0)  CALL merge_values_into_set(nlevs, p_onl%i_levels, &
            &                                        nh_pzlev_config(log_patch_id)%ilevels)
        END IF
      END DO
      p_onl => p_onl%next

    END DO ! p_onl
#endif

  END SUBROUTINE collect_requested_ipz_levels


  !------------------------------------------------------------------------------------------------
  !> Looks for variable groups ("group:xyz") and replaces them
  !
  !  @note This subroutine cannot be called directly from
  !         "read_name_list_output_namelists" because the latter
  !         routine is processed earlier, before all variable lists
  !         have been registered through "add_vars".
  !
  !  @note In more detail, this subroutine looks for variable groups
  !        ("group:xyz") and replaces them by all variables belonging
  !        to the group. Afterwards, variables can be REMOVED from
  !        this union set with the syntax "-varname". Note that typos
  !        are not detected but that the corresponding variable is
  !        simply not removed!
  !
  SUBROUTINE parse_variable_groups()
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::parse_variable_groups"
    !
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: varlist(:), grp_vars(:), new_varlist(:)
    CHARACTER(LEN=VARNAME_LEN)              :: vname, grp_name
    INTEGER                                 :: nvars, ngrp_vars, i_typ, ierrstat, &
      &                                        ivar, ntotal_vars, jvar, i,        &
      &                                        nsubtract_vars
    CHARACTER(LEN=vname_len),  POINTER      :: in_varlist(:)
    TYPE (t_output_name_list), POINTER      :: p_onl

    ntotal_vars = total_number_of_variables()
    ! temporary variables needed for variable group parsing
    ALLOCATE(varlist(ntotal_vars), grp_vars(ntotal_vars), &
      &      new_varlist(ntotal_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! -- loop over all output namelists
    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT

      ! process i_typ=ml_varlist, pl_varlist, hl_varlist, il_varlist:
      DO i_typ = 1, 4

        IF (i_typ == level_type_ml)  in_varlist => p_onl%ml_varlist
        IF (i_typ == level_type_pl)  in_varlist => p_onl%pl_varlist
        IF (i_typ == level_type_hl)  in_varlist => p_onl%hl_varlist
        IF (i_typ == level_type_il)  in_varlist => p_onl%il_varlist

        ! Get the number of variables in varlist
        nvars = 1
        DO
          IF (nvars>SIZE(in_varlist))   EXIT
          IF (in_varlist(nvars) == ' ') EXIT
          nvars = nvars + 1
        END DO
        nvars = nvars - 1

        IF (nvars>ntotal_vars)  CALL finish(routine, "Internal error: nvars > ntotal_vars")

        if (nvars > 0)  varlist(1:nvars) = in_varlist(1:nvars)
        varlist((nvars+1):ntotal_vars) = " "
        ! look for variable groups ("group:xyz") and replace them:
        DO ivar = 1, nvars
          vname = in_varlist(ivar)
          IF (INDEX(vname, GRP_PREFIX) > 0) THEN
            ! this is a group identifier
            grp_name = vname((LEN(TRIM(GRP_PREFIX))+1) : LEN(vname))
            ! loop over all variables and collects the variables names
            ! corresponding to the group "grp_name"
            CALL collect_group(grp_name, grp_vars, ngrp_vars, &
              &               loutputvars_only=.TRUE.,        &
              &               lremap_lonlat=(p_onl%remap == REMAP_REGULAR_LATLON), &
              &               opt_vlevel_type=i_typ)
            DO i=1,ngrp_vars
              grp_vars(i) = tolower(grp_vars(i))
            END DO
            ! generate varlist where "grp_name" has been replaced;
            ! duplicates are removed
            CALL insert_group(varlist, VARNAME_LEN, ntotal_vars, &
              &               TRIM(GRP_PREFIX)//TRIM(grp_name),  &
              &               grp_vars(1:ngrp_vars), new_varlist)
            varlist(:) = new_varlist(:)

            ! status output
            IF (msg_level >= 12) THEN
              CALL message(routine, "Activating group of variables: "//TRIM(grp_name))
              DO jvar=1,ngrp_vars
                CALL message(routine, "   "//TRIM(grp_vars(jvar)))
              END DO
            END IF
          END IF
        END DO

        ! Again, count the number of variables in varlist
        nvars = 1
        DO
          IF (nvars>SIZE(varlist))   EXIT
          IF (varlist(nvars) == ' ') EXIT
          nvars = nvars + 1
        END DO
        nvars = nvars - 1

        IF (i_typ == level_type_ml)  p_onl%ml_varlist(1:nvars) = varlist(1:nvars)
        IF (i_typ == level_type_pl)  p_onl%pl_varlist(1:nvars) = varlist(1:nvars)
        IF (i_typ == level_type_hl)  p_onl%hl_varlist(1:nvars) = varlist(1:nvars)
        IF (i_typ == level_type_il)  p_onl%il_varlist(1:nvars) = varlist(1:nvars)

        ! second step: look for "subtraction" of variables groups ("-varname"):
        nsubtract_vars = 0
        DO ivar = 1, nvars
          vname = TRIM(in_varlist(ivar))
          IF (vname(1:1) == "-") THEN
            nsubtract_vars = nsubtract_vars + 1
            varlist(nsubtract_vars) = TRIM(ADJUSTL(vname))
            nsubtract_vars = nsubtract_vars + 1
            varlist(nsubtract_vars) = TRIM(ADJUSTL(vname(2:)))
          END IF
        END DO
        ! remove variables
        CALL difference(in_varlist, nvars, varlist, nsubtract_vars)
      END DO ! i_typ = 1,4
      p_onl => p_onl%next

    END DO ! p_onl

    DEALLOCATE(varlist, grp_vars, new_varlist, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE parse_variable_groups


  !------------------------------------------------------------------------------------------------
  !> Initialize data structures for output.
  !
  !  This routine is called after reading the namelists AND setting up
  !  the domains and variables.
  !
  SUBROUTINE init_name_list_output(sim_step_info, &
    &                              opt_lprintlist, opt_l_is_ocean)

#ifndef NOMPI
#ifdef  __SUNPRO_F95
    INCLUDE "mpif.h"
#else
    USE mpi, ONLY: MPI_ROOT, MPI_PROC_NULL
#endif
! __SUNPRO_F95
#endif
! NOMPI

    !> Data structure containing all necessary data for mapping an
    !  output time stamp onto a corresponding simulation step index.
    TYPE (t_sim_step_info), INTENT(IN) :: sim_step_info

    LOGICAL, OPTIONAL, INTENT(in) :: opt_lprintlist
    LOGICAL, OPTIONAL, INTENT(in) :: opt_l_is_ocean

    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::init_name_list_output"

    LOGICAL                              :: l_print_list ! Flag. Enables  a list of all variables
    INTEGER                              :: i, j, nfiles, i_typ, nvl, vl_list(max_var_lists), &
      &                                     jp, idom, jg, local_i, idom_log,                  &
      &                                     grid_info_mode, ierrstat, jl, idummy, ifile,      &
      &                                     npartitions, ifile_partition
    INTEGER                              :: pe_placement(MAX_NUM_IO_PROCS)
    TYPE (t_output_name_list), POINTER   :: p_onl
    TYPE (t_output_file),      POINTER   :: p_of
    TYPE(t_list_element),      POINTER   :: element
    TYPE(t_fname_metadata)               :: fname_metadata
    TYPE(t_par_output_event),  POINTER   :: ev
    TYPE (t_sim_step_info)               :: dom_sim_step_info
    TYPE(t_cf_var),            POINTER   :: this_cf
    TYPE(timedelta),           POINTER   :: mtime_output_interval, mtime_lower_bound,          &
      &                                     mtime_interval
    TYPE(datetime),            POINTER   :: mtime_datetime, mtime_datetime_start,              &
      &                                     mtime_datetime_end
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: lower_bound_str
    CHARACTER(len=MAX_CHAR_LENGTH)       :: attname                        !< attribute name
    CHARACTER(len=MAX_CHAR_LENGTH)       :: proc_list_str                  !< string (unoccupied I/O ranks)
    LOGICAL                              :: occupied_pes(MAX_NUM_IO_PROCS) !< explicitly placed I/O ranks
    INTEGER                              :: nremaining_io_procs            !< no. of non-placed I/O ranks
    INTEGER                              :: remaining_io_procs(MAX_NUM_IO_PROCS) !< non-placed I/O ranks

    CHARACTER(LEN=MAX_DATETIME_STR_LEN)  :: output_start(MAX_TIME_INTERVALS),     &
      &                                     output_interval(MAX_TIME_INTERVALS)
    INTEGER                              :: idx, istart, iintvl,  nintvls,        &
      &                                     additional_days(MAX_TIME_INTERVALS)
    INTEGER(c_int64_t)                   :: total_ms
    LOGICAL                              :: include_last
#if !defined (__NO_ICON_ATMO__) && !defined (__NO_ICON_OCEAN__)
    CHARACTER(LEN=max_char_length)       :: comp_name
#endif
    l_print_list = .FALSE.
    CALL assign_if_present(l_print_list, opt_lprintlist)

    ! ---------------------------------------------------------------------------

    ! Optional: print list of all variables
    IF (l_print_list) THEN
      DO i = 1, nvar_lists

        IF (my_process_is_stdio()) THEN
          WRITE(message_text,'(3a, i2)') &
            'Var_list name: ',TRIM(var_lists(i)%p%name), &
            ' Patch: ',var_lists(i)%p%patch_id
          CALL message('',message_text)
          element => var_lists(i)%p%first_list_element
          DO
            IF(.NOT. ASSOCIATED(element)) EXIT
            
            IF (element%field%info%post_op%lnew_cf) THEN
              this_cf => element%field%info%post_op%new_cf
            ELSE
              this_cf => element%field%info%cf
            END IF

            WRITE (message_text,'(a,a,l1,a,a)') &
                 &     '    ',element%field%info%name,              &
                 &            element%field%info%loutput, '  ',     &
                 &            TRIM(this_cf%long_name)
            CALL message('',message_text)
            element => element%next_list_element
          ENDDO
        ENDIF
      ENDDO
    ENDIF ! IF (l_print_list)

#ifndef NOMPI
    ! Set broadcast root for intercommunicator broadcasts
    IF(my_process_is_io()) THEN
      ! Root is proc 0 on the compute PEs
      bcast_root = 0
    ELSE
      ! Special root setting for intercommunicators:
      ! The PE really sending must use MPI_ROOT, the others MPI_PROC_NULL
      IF(p_pe_work == 0) THEN
        bcast_root = MPI_ROOT
      ELSE
        bcast_root = MPI_PROC_NULL
      ENDIF
    ENDIF
#else
    ! bcast_root is not used in this case
    bcast_root = 0
#endif
! NOMPI

    ! ---------------------------------------------------------------------------

    ! Set the number of output domains depending on
    ! l_output_phys_patch

    IF(l_output_phys_patch) THEN
      n_dom_out = n_phys_dom
    ELSE
      n_dom_out = n_dom
    ENDIF

    ! Replicate physical domain setup, only the number of domains and
    ! the logical ID is needed
    IF (use_async_name_list_io .AND.  &
      & .NOT. my_process_is_mpi_test()) THEN
      CALL p_bcast(n_phys_dom, bcast_root, p_comm_work_2_io)
      DO jg = 1, n_phys_dom
        CALL p_bcast(p_phys_patch(jg)%logical_id, bcast_root, p_comm_work_2_io)
      ENDDO
    END IF

    ! reset n_dom_out (required on I/O PEs):
    IF(l_output_phys_patch) THEN
      n_dom_out = n_phys_dom
    ELSE
      n_dom_out = n_dom
    ENDIF

    ! allocate patch info data structure for unstructured and regular
    ! grids:
    ALLOCATE(patch_info(n_dom_out), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ALLOCATE(lonlat_info(n_lonlat_grids, n_dom), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! ---------------------------------------------------------------------------

    ! If dom(:) was not specified in namelist input, it is set
    ! completely to -1.  In this case all domains are wanted in the
    ! output, so set it here appropriately - this cannot be done
    ! during reading of the namelists since the number of physical
    ! domains is not known there.
    p_onl => first_output_name_list
    ! Loop over all "output_nml" namelists:
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT
      IF(p_onl%dom(1) <= 0) THEN
        DO i = 1, n_dom_out
          p_onl%dom(i) = i
        ENDDO
      ENDIF
      p_onl => p_onl%next
    ENDDO

    ! ---------------------------------------------------------------------------

    ! Set number of global cells/edges/verts and logical patch ID
    DO jp = 1, n_dom_out
      IF(l_output_phys_patch) THEN
        patch_info(jp)%log_patch_id = p_phys_patch(jp)%logical_id
        IF (.NOT. my_process_is_io()) THEN
          patch_info(jp)%p_pat_c    => p_phys_patch(jp)%comm_pat_gather_c
          patch_info(jp)%nblks_glb_c = (p_phys_patch(jp)%n_patch_cells-1)/nproma + 1
          patch_info(jp)%p_pat_e    => p_phys_patch(jp)%comm_pat_gather_e
          patch_info(jp)%nblks_glb_e = (p_phys_patch(jp)%n_patch_edges-1)/nproma + 1
          patch_info(jp)%p_pat_v    => p_phys_patch(jp)%comm_pat_gather_v
          patch_info(jp)%nblks_glb_v = (p_phys_patch(jp)%n_patch_verts-1)/nproma + 1
          patch_info(jp)%max_cell_connectivity = p_patch(patch_info(jp)%log_patch_id)%cells%max_connectivity
        END IF
      ELSE
        patch_info(jp)%log_patch_id = jp
        IF (.NOT. my_process_is_io()) THEN
          patch_info(jp)%p_pat_c    => p_patch(jp)%comm_pat_gather_c
          patch_info(jp)%nblks_glb_c = (p_patch(jp)%n_patch_cells_g-1)/nproma + 1
          patch_info(jp)%p_pat_e    => p_patch(jp)%comm_pat_gather_e
          patch_info(jp)%nblks_glb_e = (p_patch(jp)%n_patch_edges_g-1)/nproma + 1
          patch_info(jp)%p_pat_v    => p_patch(jp)%comm_pat_gather_v
          patch_info(jp)%nblks_glb_v = (p_patch(jp)%n_patch_verts_g-1)/nproma + 1
          patch_info(jp)%max_cell_connectivity = p_patch(patch_info(jp)%log_patch_id)%cells%max_connectivity
        END IF
      ENDIF
    ENDDO ! jp

    ! ---------------------------------------------------------------------------
    !--- determine the mode how to collect grid information for output:

    patch_info(:)%grid_info_mode    = GRID_INFO_NONE
    lonlat_info(:,:)%grid_info_mode = GRID_INFO_NONE

    p_onl => first_output_name_list
    ! Loop over all "output_nml" namelists:
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT
      
      ! Loop over all domains for which this name list should be used
      DO i = 1, SIZE(p_onl%dom)
        IF(p_onl%dom(i) <= 0) EXIT ! Last one was reached
        idom = p_onl%dom(i)
        ! non-existent domains are simply ignored:
        IF(p_onl%dom(i) > n_dom_out)  CYCLE
        
        IF (p_onl%output_grid) THEN
          grid_info_mode = GRID_INFO_BCAST
          ! For hexagons, we still copy grid info from file; for
          ! triangular grids we have a faster method without file access
          ! IF (max_cell_connectivity == 6)  grid_info_mode = GRID_INFO_FILE
          IF (PRESENT(opt_l_is_ocean)) THEN
            IF (opt_l_is_ocean) grid_info_mode = GRID_INFO_BCAST
          ENDIF
          IF (p_onl%remap==REMAP_REGULAR_LATLON) THEN
            lonlat_info(p_onl%lonlat_id,patch_info(idom)%log_patch_id)%grid_info_mode = grid_info_mode
          ELSE
            patch_info(idom)%grid_info_mode = grid_info_mode
          END IF
        END IF
      ENDDO ! i=1,ndom

      p_onl => p_onl%next
    ENDDO

    ! replicate grid_info_mode on I/O PEs:
    IF (use_async_name_list_io .AND.  &
      & .NOT. my_process_is_mpi_test()) THEN
      ! Go over all output domains
      DO idom = 1, n_dom_out
        CALL p_bcast(patch_info(idom)%grid_info_mode, bcast_root, p_comm_work_2_io)
      END DO
      ! A similar process as above - for the lon-lat grids
      DO jl = 1,n_lonlat_grids
        DO jg = 1,n_dom
          CALL p_bcast(lonlat_info(jl,jg)%grid_info_mode, bcast_root, p_comm_work_2_io)
        END DO
      END DO
    END IF

    ! Prepare the output of grid information: For each
    ! physical/logical patch we must collect the geographical
    ! locations of cells, edges, and vertices

    ! Pure I/O PEs may skip this...
    IF (.NOT. (use_async_name_list_io .AND. my_process_is_io())) THEN
      ! Go over all output domains
      DO idom = 1, n_dom_out
        IF (patch_info(idom)%grid_info_mode == GRID_INFO_BCAST) THEN
          ! logical domain ID
          idom_log = patch_info(idom)%log_patch_id
          CALL collect_all_grid_info(p_patch(idom_log), patch_info(idom))
        END IF
      END DO
    END IF

    ! ---------------------------------------------------------------------------

    ! If async IO is used, replicate data (mainly the variable lists) on IO procs

#ifndef NOMPI
    IF (use_async_name_list_io) THEN
      CALL replicate_data_on_io_procs

      ! Clear patch_info fields clon, clat, etc. (especially on work
      ! PE 0) since they aren't needed there any longer.
      IF ( (.NOT. my_process_is_io()) .AND. &
        &  (.NOT. my_process_is_mpi_test())) THEN
        ! Go over all output domains (deallocation is skipped if data
        ! structures were not allocated)
        DO idom = 1, n_dom_out
          CALL deallocate_all_grid_info(patch_info(idom))
        END DO
      END IF
    END IF
#endif
! NOMPI

    ! Set the number of domains in output and the patch reorder information
    CALL set_patch_info

    ! Date-time computations:
    !
    ! There are two alternative implementations for setting the output
    ! intervals, "output_bounds" and "output_start" / "output_end" /
    ! "output_interval". The former defines the output events relative
    ! to the simulation start (in seconds) and the latter define the
    ! output events by setting ISO8601-conforming date-time strings.
    !
    ! If the user has set values for "output_bounds", we compute the
    ! "output_start" / "output_end" / "output_interval" from this
    ! info:
    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl))  EXIT

      ! there may be multiple "output_bounds" intervals, consider all:        
      DO idx=1,MAX_TIME_INTERVALS
        istart = (idx-1)*3
        IF (p_onl%output_bounds(istart+1) == -1._wp) CYCLE
        
        CALL get_datetime_string(p_onl%output_start(idx), &
          &                      sim_step_info%sim_start, INT(p_onl%output_bounds(istart+1)))
        CALL get_datetime_string(p_onl%output_end(idx),   &
          &                      sim_step_info%sim_start, INT(p_onl%output_bounds(istart+2)))
        CALL get_duration_string(INT(p_onl%output_bounds(istart+3)), p_onl%output_interval(idx), &
          &                      p_onl%additional_days(idx))
        IF (my_process_is_stdio()) THEN
          IF (p_onl%additional_days(idx) == 0) THEN
            WRITE (0,*) "setting output bounds as ", TRIM(p_onl%output_start(idx)), " / ", &
              &                                      TRIM(p_onl%output_end(idx)),   " / ", &
              &                                      TRIM(p_onl%output_interval(idx))
          ELSE
            WRITE (0,*) "setting output bounds as ", TRIM(p_onl%output_start(idx)), " / ", &
              &                                      TRIM(p_onl%output_end(idx)),   " / ", &
              &                                      TRIM(p_onl%output_interval(idx)), " + ", &
              &                                      p_onl%additional_days(idx), " days"
          END IF
        END IF
      END DO

      !--- consistency check: do not allow output intervals < dtime:

      ! there may be multiple "output_bounds" intervals, consider all:        
      INTVL_LOOP : DO idx=1,MAX_TIME_INTERVALS
        IF (p_onl%additional_days(idx) == 0) THEN
          IF (TRIM(p_onl%output_start(idx)) == '') CYCLE INTVL_LOOP

          ! compare start date and end date: if these are equal, then
          ! the interval does not matter and must not be checked.
          mtime_datetime_start => newDatetime(p_onl%output_start(idx))
          mtime_datetime_end   => newDatetime(p_onl%output_end(idx))
          IF (mtime_datetime_end > mtime_datetime_start) THEN
            mtime_output_interval => newTimedelta(TRIM(p_onl%output_interval(idx)))
            
            !Special case for very small time steps
            IF(sim_step_info%dtime .LT. 1._wp)THEN
              CALL get_duration_string_real(REAL(sim_step_info%dtime), &
                &                           lower_bound_str)
              idummy = 0
            ELSE  
              CALL get_duration_string(INT(sim_step_info%dtime), &
                &                      lower_bound_str, idummy)
            END IF

            IF (idummy > 0)  CALL finish(routine, "Internal error: get_duration_string")
            mtime_lower_bound     => newTimedelta(TRIM(lower_bound_str))
            IF (mtime_output_interval < mtime_lower_bound) THEN
              CALL finish(routine, "Output interval "//TRIM(p_onl%output_interval(idx))//" < dtime !")
            END IF
            CALL deallocateTimedelta(mtime_output_interval)
            CALL deallocateTimedelta(mtime_lower_bound)
          END IF
          CALL deallocateDatetime(mtime_datetime_start)
          CALL deallocateDatetime(mtime_datetime_end)
        END IF
      END DO INTVL_LOOP

      p_onl => p_onl%next
    ENDDO

    ! Get the number of output files needed (by counting the domains per name list)

    p_onl => first_output_name_list
    nfiles = 0
    DO

      IF(.NOT.ASSOCIATED(p_onl))  EXIT

      DO i = 1, SIZE(p_onl%dom)
        IF(p_onl%dom(i) <= 0) EXIT ! Last one was reached
        ! non-existent domains are simply ignored:
        IF(p_onl%dom(i) > n_dom_out)  CYCLE

        DO i_typ = 1, 4
          ! Check if name_list has variables of corresponding type,
          ! then increase file counter.
          SELECT CASE(i_typ)
          CASE (level_type_ml)
            IF (p_onl%ml_varlist(1) == ' ') CYCLE
          nfiles = nfiles + p_onl%stream_partitions_ml
          CASE (level_type_pl)
            IF (p_onl%pl_varlist(1) == ' ') CYCLE
          nfiles = nfiles + p_onl%stream_partitions_pl
          CASE (level_type_hl) 
            IF (p_onl%hl_varlist(1) == ' ') CYCLE
            nfiles = nfiles + p_onl%stream_partitions_hl
          CASE (level_type_il)
            IF (p_onl%il_varlist(1) == ' ') CYCLE
            nfiles = nfiles + p_onl%stream_partitions_il
          END SELECT
        ENDDO
      ENDDO

      p_onl => p_onl%next

    ENDDO
    WRITE(message_text,'(a,i4)') 'Number of name list output files: ',nfiles
    CALL message(routine,message_text)

    ! Allocate output_file struct for all output files

    ALLOCATE(output_file(nfiles))

    output_file(:)%cdiFileID  = CDI_UNDEFID ! i.e. not opened
    output_file(:)%cdiVlistId = CDI_UNDEFID ! i.e. not defined


    ! --------------------------------------------------------------------------------------
    ! Loop over all output namelists, set up the output_file struct for all associated files
    ! --------------------------------------------------------------------------------------

    p_onl => first_output_name_list
    ifile = 0
    LOOP_NML : DO

      IF(.NOT.ASSOCIATED(p_onl)) EXIT

      ! Loop over all domains for which this name list should be used

      LOOP_DOM : DO i = 1, SIZE(p_onl%dom)
        IF(p_onl%dom(i) <= 0) EXIT ! Last one was reached
        idom = p_onl%dom(i)
        ! non-existent domains are simply ignored:
        IF(p_onl%dom(i) > n_dom_out)  CYCLE

        ! Loop over model/pressure/height levels

        DO i_typ = 1, 4

          ! Check if name_list has variables of corresponding type
          SELECT CASE(i_typ)
          CASE (level_type_ml)
            IF (p_onl%ml_varlist(1) == ' ') CYCLE
            npartitions     = p_onl%stream_partitions_ml
            pe_placement(:) = p_onl%pe_placement_ml(:)
          CASE (level_type_pl)
            IF (p_onl%pl_varlist(1) == ' ') CYCLE
            npartitions     = p_onl%stream_partitions_pl
            pe_placement(:) = p_onl%pe_placement_pl(:)
          CASE (level_type_hl) 
            IF (p_onl%hl_varlist(1) == ' ') CYCLE
            npartitions     = p_onl%stream_partitions_hl
            pe_placement(:) = p_onl%pe_placement_hl(:)
          CASE (level_type_il)
            IF (p_onl%il_varlist(1) == ' ') CYCLE
            npartitions     = p_onl%stream_partitions_il
            pe_placement(:) = p_onl%pe_placement_il(:)
          END SELECT

          IF (npartitions > 1) THEN
            WRITE(message_text,'(a,i4,a)') "Fork file into: ", npartitions, " concurrent parts."
            CALL message(routine, message_text)
          END IF

          ! Split one namelist into concurrent, alternating files:
          DO ifile_partition = 1,npartitions

              ifile = ifile+1
              p_of => output_file(ifile)
              p_of%ilev_type = i_typ

              ! Fill data members of "t_output_file" data structures
              p_of%filename_pref   = TRIM(p_onl%output_filename)
              p_of%phys_patch_id   = idom
              p_of%log_patch_id    = patch_info(idom)%log_patch_id
              p_of%output_type     = p_onl%filetype
              p_of%name_list       => p_onl
              p_of%remap           = p_onl%remap
              p_of%cdiCellGridID   = CDI_UNDEFID
              p_of%cdiEdgeGridID   = CDI_UNDEFID
              p_of%cdiVertGridID   = CDI_UNDEFID
              p_of%cdiLonLatGridID = CDI_UNDEFID
              p_of%cdiTaxisID      = CDI_UNDEFID
              p_of%cdiZaxisID(:)   = CDI_UNDEFID
              p_of%cdiVlistID      = CDI_UNDEFID

              p_of%npartitions     = npartitions
              p_of%ifile_partition = ifile_partition
              
              ! (optional:) explicitly specified I/O rank
              p_of%io_proc_id      = -1 ! undefined MPI rank
              p_of%pe_placement    = pe_placement(ifile_partition)

              ! Select all var_lists which belong to current logical domain and i_typ
              nvl = 0
              DO j = 1, nvar_lists
                
                IF(.NOT. var_lists(j)%p%loutput) CYCLE
                ! patch_id in var_lists always corresponds to the LOGICAL domain
                IF(var_lists(j)%p%patch_id /= patch_info(idom)%log_patch_id) CYCLE
                
                IF(i_typ /= var_lists(j)%p%vlevel_type) CYCLE
                
                nvl = nvl + 1
                vl_list(nvl) = j
                
              ENDDO

              SELECT CASE(i_typ)
              CASE(level_type_ml)
                CALL add_varlist_to_output_file(p_of,vl_list(1:nvl),p_onl%ml_varlist)
              CASE(level_type_pl)
                CALL add_varlist_to_output_file(p_of,vl_list(1:nvl),p_onl%pl_varlist)
              CASE(level_type_hl)
                CALL add_varlist_to_output_file(p_of,vl_list(1:nvl),p_onl%hl_varlist)
              CASE(level_type_il)
                CALL add_varlist_to_output_file(p_of,vl_list(1:nvl),p_onl%il_varlist)
              END SELECT
        
          END DO ! ifile_partition
            
        ENDDO ! i_typ
          
      ENDDO LOOP_DOM ! i=1,ndom

      p_onl => p_onl%next
      
    ENDDO LOOP_NML

    ! ------------------------------------------------------
    ! Set ID of process doing I/O
    ! ------------------------------------------------------

    ! --- First, set those MPI ranks which were explicitly specified
    !     by the user:
    !
    occupied_pes(:) = .FALSE.
    DO i = 1, nfiles
      IF(use_async_name_list_io) THEN
        ! Asynchronous I/O
        !
        ! MPI ranks "p_io_pe0 ... (p_io_pe0+process_mpi_io_size-1)" are available.
        IF (process_mpi_io_size == 0) CALL finish(routine, "Asynchronous I/O but no IO procs!")
        IF ((output_file(i)%pe_placement /= -1) .AND. &
          & ((output_file(i)%pe_placement < 0) .OR.   &
          &  (output_file(i)%pe_placement > process_mpi_io_size))) THEN
          CALL finish(routine, "Invalid explicit placement of IO rank!")
        END IF

        IF (output_file(i)%pe_placement /= -1) THEN
          output_file(i)%io_proc_id = p_io_pe0 + output_file(i)%pe_placement
          occupied_pes(output_file(i)%pe_placement+1) = .TRUE.
        END IF
      ELSE
        ! Normal I/O done by the standard I/O processor
        !
        ! Only MPI rank "process_mpi_stdio_id" is available.
        IF ((output_file(i)%pe_placement /= -1) .AND. &
          & (output_file(i)%pe_placement /=  0)) &
          &  CALL finish(routine, "Invalid explicit placement of IO rank!")

        IF (p_test_run .AND. .NOT. my_process_is_mpi_test()) THEN
          output_file(i)%io_proc_id = process_mpi_stdio_id + 1
        ELSE
          output_file(i)%io_proc_id = process_mpi_stdio_id
        END IF
      ENDIF
    ENDDO

    ! --- Build a list of MPI ranks that have not yet been occupied:
    !
    nremaining_io_procs = 0
    DO i=1,process_mpi_io_size
      IF (.NOT. occupied_pes(i)) THEN
        nremaining_io_procs = nremaining_io_procs + 1
        remaining_io_procs(nremaining_io_procs) = (i-1)
      END IF
    END DO
    ! status output, if some ranks were explicitly specified
    IF ((process_mpi_io_size /= nremaining_io_procs) .AND. my_process_is_stdio()) THEN
      WRITE (0,*) " "
      WRITE (0,*) "I/O : Explicit placement of I/O ranks:"
      DO i = 1, nfiles
        IF (output_file(i)%pe_placement /= -1) THEN
          WRITE (0,'(a,i0,a,i0)') "    file #", i, " placed on rank #", output_file(i)%io_proc_id
        END IF
      END DO
      IF (nremaining_io_procs > 0) THEN
        CALL sort_and_compress_list(remaining_io_procs(1:nremaining_io_procs), proc_list_str)
        WRITE (0,*) "I/O : Remaining I/O ranks: # ", TRIM(proc_list_str)
      END IF
    END IF

    ! --- Then, set MPI ranks in a Round-Robin fashion for those
    !     namelists which had no explicitly specified "pe_placement":
    !
    ! status print-out only when some PEs were explicitly set.
    IF ((process_mpi_io_size /= nremaining_io_procs) .AND. &
      & my_process_is_stdio()                        .AND. &
      & (nremaining_io_procs > 0)                    .AND. &
      & ANY(output_file(1:nfiles)%pe_placement == -1)) THEN
      WRITE (0,*) " "
      WRITE (0,*) "I/O : Round-Robin placement of I/O ranks:"
    END IF
    j = 0
    DO i = 1, nfiles
      IF (output_file(i)%pe_placement /= -1) CYCLE
      IF(use_async_name_list_io) THEN
        IF (nremaining_io_procs == 0) THEN
          CALL finish(routine, "No I/O proc left after explicit placement!")
        END IF
        ! Asynchronous I/O
        j = j + 1
        output_file(i)%io_proc_id = p_io_pe0 + remaining_io_procs(MOD(j-1,nremaining_io_procs) + 1)
        IF ((process_mpi_io_size /= nremaining_io_procs) .AND. my_process_is_stdio()) THEN
          WRITE (0,'(a,i0,a,i0)') "    file #", i, " placed on rank #", output_file(i)%io_proc_id
        END IF
      ENDIF
    ENDDO
    IF ((process_mpi_io_size /= nremaining_io_procs) .AND. my_process_is_stdio()) THEN
      WRITE (0,*) " "
    END IF

    ! ------------------------------------------------------
    ! Create I/O event data structures:
    ! ------------------------------------------------------
    ! 
    !  Regular output is triggered at
    ! so-called "output event steps". The completion of an output
    ! event step is communicated via non-blocking MPI messages to the
    ! root I/O PE, which keeps track of the overall event status. This
    ! event handling is static, i.e. all event occurrences are
    ! pre-defined during the initialization.
    !
    local_i = 0
    DO i = 1, nfiles
      p_of  => output_file(i)

      IF (use_async_name_list_io  .AND. &
        & my_process_is_io()      .AND. &
        & p_of%io_proc_id /= p_pe) THEN
        NULLIFY(p_of%out_event)
        CYCLE
      END IF

      p_onl => p_of%name_list
      ! pack file-name meta-data into a derived type to pass them
      ! to "new_parallel_output_event":
      fname_metadata%steps_per_file             = p_onl%steps_per_file
      fname_metadata%steps_per_file_inclfirst   = p_onl%steps_per_file_inclfirst
      fname_metadata%file_interval              = p_onl%file_interval
      fname_metadata%phys_patch_id              = p_of%phys_patch_id
      fname_metadata%ilev_type                  = p_of%ilev_type
      fname_metadata%filename_format            = TRIM(p_onl%filename_format)
      fname_metadata%filename_pref              = TRIM(p_of%filename_pref)
      fname_metadata%npartitions                = p_of%npartitions
      fname_metadata%ifile_partition            = p_of%ifile_partition
      ! set user-specified filename extension or use the default
      ! extension:
      IF (TRIM(p_onl%filename_extn) == "default") THEN
        fname_metadata%extn                     = TRIM(get_file_extension(p_onl%filetype))
      ELSE
        fname_metadata%extn                     = TRIM(p_onl%filename_extn)
      END IF

      IF (is_restart_run() .AND. .NOT. time_config%is_relative_time) THEN
        ! Restart case: Get starting index of ouput from restart file
        !               (if there is such an attribute available).
        WRITE(attname,'(a,i2.2)') 'output_jfile_',i
        CALL get_restart_attribute(TRIM(attname), fname_metadata%jfile_offset, opt_default=0)
      ELSE
        fname_metadata%jfile_offset             = 0
      END IF

      ! set model domain start/end time
      dom_sim_step_info = sim_step_info
      CALL get_datetime_string(dom_sim_step_info%dom_start_time, time_config%ini_datetime, NINT(start_time(p_of%log_patch_id)))
      IF (end_time(p_of%log_patch_id) < DEFAULT_ENDTIME) THEN
        CALL get_datetime_string(dom_sim_step_info%dom_end_time,   time_config%ini_datetime, NINT(end_time(p_of%log_patch_id)))
      ELSE
        dom_sim_step_info%dom_end_time = dom_sim_step_info%sim_end
      END IF
      local_i = local_i + 1

      ! special treatment of ocean model: model_date/run_start is the time at
      ! the beginning of the timestep. Output is written at the end of the
      ! timestep
      IF (is_restart_run() .AND. my_process_is_ocean()) THEN
        IF (TRIM(p_onl%output_start(2)) /= '') &
          CALL finish(routine, "Not implemented for ocean model with restart!")
        CALL get_datetime_string(p_onl%output_start(1), &
          &                      time_config%cur_datetime,                &
          &                      opt_td_string=p_onl%output_interval(1))
      ENDIF

      include_last    = p_onl%include_last
      output_interval = p_onl%output_interval
      additional_days = p_onl%additional_days
      output_start    = p_onl%output_start

      ! Handle the case that one namelist has been split into
      ! concurrent, alternating files:
      !
      IF (p_of%npartitions > 1) THEN
        ! count the number of different time intervals for this event (usually 1)
        nintvls = 0
        DO
          IF (TRIM(output_start(nintvls+1)) == '') EXIT
          nintvls = nintvls + 1
          IF (nintvls == MAX_TIME_INTERVALS) EXIT
        END DO
    
        DO iintvl=1,nintvls
          mtime_interval => newTimedelta(output_interval(iintvl))
          mtime_datetime => newDatetime(output_start(iintvl))
          !
          ! - The start_date gets an offset of
          !         "(ifile_partition - 1) * output_interval"
          DO ifile=1,(p_of%ifile_partition-1)
            mtime_datetime = mtime_datetime + mtime_interval
          END DO
          CALL datetimeToString(mtime_datetime, output_start(iintvl))
          ! - The output_interval is replaced by "
          !         "npartitions * output_interval"
          total_ms = getTotalMilliSecondsTimeDelta(mtime_interval, mtime_datetime)
          total_ms = total_ms + additional_days(iintvl)*86400000
          total_ms = total_ms * p_of%npartitions
          CALL get_duration_string(INT(total_ms/1000), output_interval(iintvl), additional_days(iintvl))
          IF (p_of%ifile_partition == 1) THEN
            WRITE(message_text,'(a,a)') "File stream partitioning: total output interval = ", &
              &                         output_interval(iintvl)
            CALL message(routine, message_text)
          END IF
          ! - The "include_last" flag is set to .FALSE.
          include_last                  = .FALSE.
          ! - The "steps_per_file" counter is set to 1
          fname_metadata%steps_per_file = 1        
          ! - The "steps_per_file_inclfirst" flag is set to .FALSE.
          fname_metadata%steps_per_file_inclfirst = .FALSE.
          !
          CALL deallocateTimedelta(mtime_interval)
          CALL deallocateDatetime(mtime_datetime)
        END DO ! iintvl
      END IF

      ! ------------------------------------------------------------------------------------------
      ! --- I/O PEs communicate their event data, the other PEs create
      ! --- the event data only locally for their own event control:
      p_of%out_event => new_parallel_output_event(p_onl%ready_file,                              &
        &                  output_start, p_onl%output_end, output_interval,                      &
        &                  additional_days, include_last,                                        &
        &                  dom_sim_step_info, fname_metadata, compute_matching_sim_steps,        &
        &                  generate_output_filenames, local_i, p_comm_io)
      ! ------------------------------------------------------------------------------------------
      IF (dom_sim_step_info%jstep0 > 0) &
        &  CALL set_event_to_simstep(p_of%out_event, dom_sim_step_info%jstep0 + 1, &
        &                            is_restart_run(), lrecover_open_file=.TRUE.)
    END DO

    ! tell the root I/O process that all output event data structures
    ! have been created:
    CALL complete_event_setup(p_comm_io)

    ! -----------------------------------------------------------
    ! The root I/O MPI rank asks all participating I/O PEs for their
    ! output event info and generates a unified output event,
    ! indicating which PE performs a write process at which step.
    all_events => union_of_all_events(compute_matching_sim_steps, generate_output_filenames, p_comm_io, &
         &                               p_comm_work_io, process_work_io0)

    IF (dom_sim_step_info%jstep0 > 0) &
      &  CALL set_event_to_simstep(all_events, dom_sim_step_info%jstep0 + 1, &
      &                            is_restart_run(), lrecover_open_file=.TRUE.)
    ! print a table with all output events
    IF (.NOT. my_process_is_mpi_test()) THEN
       IF ((      use_async_name_list_io .AND. my_process_is_mpi_ioroot()) .OR.  &
            & (.NOT. use_async_name_list_io .AND. my_process_is_mpi_workroot())) THEN
          CALL print_output_event(all_events)                                       ! screen output
          IF (dom_sim_step_info%jstep0 > 0) THEN
#if !defined (__NO_ICON_ATMO__) && !defined (__NO_ICON_OCEAN__)
             IF ( is_coupled_run() ) THEN
               comp_name = TRIM(get_my_process_name())
               CALL print_output_event(all_events, &
                 ! ASCII file output:
      & opt_filename="output_schedule_"//TRIM(comp_name)//"_steps_"//TRIM(int2string(dom_sim_step_info%jstep0))//"+.txt") 
             ELSE
               CALL print_output_event(all_events, &
      & opt_filename="output_schedule_steps_"//TRIM(int2string(dom_sim_step_info%jstep0))//"+.txt") ! ASCII file output
             ENDIF
#else
             CALL print_output_event(all_events, &
      &        opt_filename="output_schedule_steps_"//TRIM(int2string(dom_sim_step_info%jstep0))//"+.txt") ! ASCII file output
#endif
          ELSE
#if !defined (__NO_ICON_ATMO__) && !defined (__NO_ICON_OCEAN__)
             IF ( is_coupled_run() ) THEN
                comp_name = TRIM(get_my_process_name())
                CALL print_output_event(all_events, opt_filename="output_schedule_"//TRIM(comp_name)//".txt") ! ASCII file output
             ELSE
                CALL print_output_event(all_events, opt_filename="output_schedule.txt") ! ASCII file output
             ENDIF
#else
             CALL print_output_event(all_events, opt_filename="output_schedule.txt") ! ASCII file output
#endif
          END IF
       END IF
    END IF

    ! If async IO is used, initialize the memory window for communication
#ifndef NOMPI
    IF(use_async_name_list_io) CALL init_memory_window
#endif
! NOMPI
    
    ! Initial launch of non-blocking requests to all participating PEs
    ! to acknowledge the completion of the next output event
    IF (.NOT. my_process_is_mpi_test()) THEN
      IF ((      use_async_name_list_io .AND. my_process_is_mpi_ioroot()) .OR.  &
        & (.NOT. use_async_name_list_io .AND. my_process_is_mpi_workroot())) THEN
        ev => all_events
        HANDLE_COMPLETE_STEPS : DO
          IF (.NOT. ASSOCIATED(ev)) EXIT HANDLE_COMPLETE_STEPS
          CALL trigger_output_step_irecv(ev)
          ev => ev%next
        END DO HANDLE_COMPLETE_STEPS
      END IF
    END IF

    CALL message(routine,'Done')

  END SUBROUTINE init_name_list_output


  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE add_varlist_to_output_file(p_of, vl_list, varlist)

    TYPE (t_output_file), INTENT(INOUT) :: p_of
    INTEGER,              INTENT(IN)    :: vl_list(:)
    CHARACTER(LEN=*),     INTENT(IN)    :: varlist(:)
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::add_varlist_to_output_file"
    INTEGER                       :: ivar, nvars, i, iv, tl, grid_of, grid_var
    LOGICAL                       :: found
    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_desc),     TARGET  :: var_desc   !< variable descriptor
    TYPE(t_var_desc),     POINTER :: p_var_desc               !< variable descriptor (pointer)
    TYPE(t_cf_var),       POINTER :: this_cf


    ! Get the number of variables in varlist
    nvars = 0
    DO ivar = 1, SIZE(varlist)
      IF(varlist(ivar) == ' ') EXIT ! Last one reached
      IF (.NOT. is_grid_info_var(varlist(ivar)))  nvars = nvars + 1
    ENDDO

    ! Allocate a list of variable descriptors:
    p_of%max_vars = nvars
    p_of%num_vars = 0
    ALLOCATE(p_of%var_desc(p_of%max_vars))
    DO ivar = 1,nvars
      ! Nullify pointers in p_of%var_desc
      p_of%var_desc(ivar)%r_ptr => NULL()
      p_of%var_desc(ivar)%i_ptr => NULL()
      DO i = 1, max_time_levels
        p_of%var_desc(ivar)%tlev_rptr(i)%p => NULL()
        p_of%var_desc(ivar)%tlev_iptr(i)%p => NULL()
      ENDDO
    END DO ! ivar

    ! Allocate array of variable descriptions
    DO ivar = 1,nvars
      IF (is_grid_info_var(varlist(ivar)))  CYCLE

      found = .FALSE.
      ! Nullify pointers
      var_desc%r_ptr => NULL()
      var_desc%i_ptr => NULL()
      DO i = 1, max_time_levels
        var_desc%tlev_rptr(i)%p => NULL()
        var_desc%tlev_iptr(i)%p => NULL()
      ENDDO

      ! Loop over all var_lists listed in vl_list to find the variable
      ! Please note that there may be several variables with different time levels,
      ! we just add unconditionally all with the name varlist(ivar).
      ! Remark: The different time levels may appear in different lists
      ! or in the same list, the code will accept both

      DO i = 1, SIZE(vl_list)

        iv = vl_list(i)

        element => NULL()
        DO
          IF(.NOT.ASSOCIATED(element)) THEN
            element => var_lists(iv)%p%first_list_element
          ELSE
            element => element%next_list_element
          ENDIF
          IF(.NOT.ASSOCIATED(element)) EXIT

          ! Do not inspect element if output is disabled
          IF(.NOT.element%field%info%loutput) CYCLE

          IF (p_of%name_list%remap==REMAP_REGULAR_LATLON) THEN
            ! If lon-lat variable is requested, skip variable if it
            ! does not correspond to the same lon-lat grid:
            IF (element%field%info%hgrid /= GRID_REGULAR_LONLAT) CYCLE
            grid_of  = p_of%name_list%lonlat_id
            grid_var = element%field%info%hor_interp%lonlat_id
            IF (grid_of /= grid_var) CYCLE
          ELSE
            ! On the other hand: If no lon-lat interpolation is
            ! requested for this output file, skip all variables of
            ! this kind:
            IF (element%field%info%hgrid == GRID_REGULAR_LONLAT) CYCLE
          END IF ! (remap/=REMAP_REGULAR_LATLON)

          ! Do not inspect element if it is a container
          IF(element%field%info%lcontainer) CYCLE

          ! get time level
          tl = get_var_timelevel(element%field)

          ! Check for matching name
          IF(TRIM(varlist(ivar)) /= TRIM(tolower(get_var_name(element%field)))) CYCLE

          ! Found it, add it to the variable list of output file
          p_var_desc => var_desc

          IF(tl == -1) THEN
            ! Not time level dependent
            IF(found) CALL finish(routine,'Duplicate var name: '//TRIM(varlist(ivar)))
            p_var_desc%r_ptr    => element%field%r_ptr
            p_var_desc%i_ptr    => element%field%i_ptr
            p_var_desc%info     =  element%field%info
            p_var_desc%info_ptr => element%field%info
          ELSE
            IF(found) THEN
              ! We have already the info field, make some plausibility checks:
              IF(ANY(p_var_desc%info%used_dimensions(:) /=  &
                element%field%info%used_dimensions(:))) THEN
                CALL message(routine, "Var "//TRIM(element%field%info%name))
                CALL finish(routine,'Dimension mismatch TL variable: '//TRIM(varlist(ivar)))
              END IF
              ! There must not be a TL independent variable with the same name
              IF (ASSOCIATED(p_var_desc%r_ptr) .OR. ASSOCIATED(p_var_desc%i_ptr)) &
                CALL finish(routine,'Duplicate var name: '//TRIM(varlist(ivar)))
              ! Maybe some more members of info should be tested ...
            ELSE
              ! Variable encountered the first time, set info field ...
              p_var_desc%info = element%field%info
              ! ... and set name without .TL# suffix
              p_var_desc%info%name = TRIM(get_var_name(element%field))
            ENDIF

            IF (ASSOCIATED(p_var_desc%tlev_rptr(tl)%p) .OR. ASSOCIATED(p_var_desc%tlev_iptr(tl)%p)) &
              CALL finish(routine, 'Duplicate time level for '//TRIM(element%field%info%name))
            p_var_desc%tlev_rptr(tl)%p => element%field%r_ptr
            p_var_desc%tlev_iptr(tl)%p => element%field%i_ptr
            p_var_desc%info_ptr        => element%field%info
          ENDIF

          found = .TRUE.
        ENDDO

      ENDDO ! i = 1, SIZE(vl_list)

      ! Check that at least one element with this name has been found

      IF (.NOT. found) THEN

        DO i = 1, nvar_lists
          IF (my_process_is_stdio()) THEN
            WRITE(message_text,'(3a, i2)') &
                 'Variable list name: ',TRIM(var_lists(i)%p%name), &
                 ' Patch: ',var_lists(i)%p%patch_id
            CALL message('',message_text)
            element => var_lists(i)%p%first_list_element
            DO
              IF(.NOT. ASSOCIATED(element)) EXIT

              IF (element%field%info%post_op%lnew_cf) THEN
                this_cf => element%field%info%post_op%new_cf
              ELSE
                this_cf => element%field%info%cf
              END IF

              WRITE (message_text,'(a,a,l1,a,a)') &
                   &     '    ',element%field%info%name,              &
                   &            element%field%info%loutput, '  ',     &
                   &            TRIM(this_cf%long_name)
              CALL message('',message_text)
              element => element%next_list_element
            ENDDO
          ENDIF
        ENDDO

        CALL finish(routine,'Output name list variable not found: '//TRIM(varlist(ivar)))
      ENDIF

      ! append variable descriptor to list
      CALL add_var_desc(p_of, var_desc)

    ENDDO ! ivar = 1,nvars

  END SUBROUTINE add_varlist_to_output_file


  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE set_patch_info()
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_patch_info"
    INTEGER :: jp, jl, jg

    DO jp = 1, n_dom_out

      jl = patch_info(jp)%log_patch_id

      IF(.NOT.my_process_is_io()) THEN
        ! Set reorder_info on work and test PE
        CALL set_reorder_info(jp, p_patch(jl)%n_patch_cells_g, p_patch(jl)%n_patch_cells,             &
          &                   p_patch(jl)%cells%decomp_info%owner_mask, p_patch(jl)%cells%phys_id,    &
          &                   p_patch(jl)%cells%decomp_info%glb_index, patch_info(jp)%cells )

        CALL set_reorder_info(jp, p_patch(jl)%n_patch_edges_g, p_patch(jl)%n_patch_edges,             &
          &                   p_patch(jl)%edges%decomp_info%owner_mask, p_patch(jl)%edges%phys_id,    &
          &                   p_patch(jl)%edges%decomp_info%glb_index, patch_info(jp)%edges )
        
        CALL set_reorder_info(jp, p_patch(jl)%n_patch_verts_g, p_patch(jl)%n_patch_verts,             &
          &                   p_patch(jl)%verts%decomp_info%owner_mask, p_patch(jl)%verts%phys_id,    &
          &                   p_patch(jl)%verts%decomp_info%glb_index, patch_info(jp)%verts )
        ! Set grid_filename on work and test PE
        patch_info(jp)%grid_filename = TRIM(p_patch(jl)%grid_filename)
        ! Set UUID on work and test PE
        patch_info(jp)%grid_uuid = p_patch(jl)%grid_uuid
        ! Set information about numberOfGridUsed on work and test PE
        patch_info(jp)%number_of_grid_used = number_of_grid_used(jl)

        patch_info(jp)%max_cell_connectivity = p_patch(jl)%cells%max_connectivity

      ENDIF
#ifndef NOMPI
      IF(use_async_name_list_io .AND. .NOT. my_process_is_mpi_test()) THEN
        ! Transfer reorder_info to IO PEs
        CALL transfer_reorder_info(patch_info(jp)%cells, patch_info(jp)%grid_info_mode)
        CALL transfer_reorder_info(patch_info(jp)%edges, patch_info(jp)%grid_info_mode)
        CALL transfer_reorder_info(patch_info(jp)%verts, patch_info(jp)%grid_info_mode)
        CALL p_bcast(patch_info(jp)%grid_filename, bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(jp)%grid_uuid%data, SIZE(patch_info(jp)%grid_uuid%data),  &
          &          bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(jp)%number_of_grid_used, bcast_root, p_comm_work_2_io)
      ENDIF
#endif
! NOMPI

    ENDDO ! jp

#ifndef __NO_ICON_ATMO__
    ! A similar process as above - for the lon-lat grids
    DO jl = 1,n_lonlat_grids
      DO jg = 1,n_dom
        IF (.NOT. lonlat_grid_list(jl)%l_dom(jg)) CYCLE
        IF(.NOT.my_process_is_io()) THEN
          ! Set reorder_info on work and test PE
          CALL set_reorder_info_lonlat(lonlat_grid_list(jl)%grid,      &
            &                          lonlat_grid_list(jl)%intp(jg),  &
            &                          lonlat_info(jl,jg))
        ENDIF
#ifndef NOMPI
        IF(use_async_name_list_io .AND. .NOT. my_process_is_mpi_test()) THEN
          ! Transfer reorder_info to IO PEs
          CALL transfer_reorder_info(lonlat_info(jl,jg)%ri, lonlat_info(jl,jg)%grid_info_mode)
        ENDIF
#endif
! NOMPI
      END DO ! jg
    ENDDO ! jl
#endif
! #ifndef __NO_ICON_ATMO__
  END SUBROUTINE set_patch_info


  !------------------------------------------------------------------------------------------------
  !> Sets the reorder_info for cells/edges/verts
  !  ATTENTION: This routine must only be called on work and test PE (i.e. not on IO PEs)
  !             The arguments don't make sense on the IO PEs anyways
  !
  SUBROUTINE set_reorder_info(phys_patch_id, n_points_g, n_points, owner_mask, phys_id, &
                              glb_index, p_ri)

    INTEGER, INTENT(IN) :: phys_patch_id   ! Physical patch ID
    INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
    INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
    LOGICAL, INTENT(IN) :: owner_mask(:,:) ! owner_mask for logical patch
    INTEGER, INTENT(IN) :: phys_id(:,:)    ! phys_id for logical patch
    INTEGER, INTENT(IN) :: glb_index(:)    ! glb_index for logical patch
    TYPE(t_reorder_info), INTENT(INOUT) :: p_ri ! Result: reorder info
    ! local variables
    INTEGER :: i, n, il, ib, mpierr
    LOGICAL, ALLOCATABLE :: phys_owner_mask(:) ! owner mask for physical patch
    INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:), reorder_index_log_dom(:)

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_reorder_info"

    ! Just for safety
    IF(my_process_is_io()) CALL finish(routine, 'Must not be called on IO PEs')

    p_ri%n_log = n_points_g ! Total points in logical domain

    ! Set the physical patch owner mask

    ALLOCATE(phys_owner_mask(n_points))
    DO i = 1, n_points
      il = idx_no(i)
      ib = blk_no(i)
      phys_owner_mask(i) = owner_mask(il,ib)
      IF(l_output_phys_patch) &
        phys_owner_mask(i) = phys_owner_mask(i) .AND. (phys_id(il,ib) == phys_patch_id)
    ENDDO

    ! Get number of owned cells/edges/verts (without halos, physical patch only)

    p_ri%n_own = COUNT(phys_owner_mask(:))

    ! Set index arrays to own cells/edges/verts

    ALLOCATE(p_ri%own_idx(p_ri%n_own))
    ALLOCATE(p_ri%own_blk(p_ri%n_own))
    ALLOCATE(glbidx_own(p_ri%n_own)) ! Global index of my own points

    n = 0
    DO i = 1, n_points
      IF(phys_owner_mask(i)) THEN
        n = n+1
        p_ri%own_idx(n) = idx_no(i)
        p_ri%own_blk(n) = blk_no(i)
        glbidx_own(n)   = glb_index(i)
      ENDIF
    ENDDO

    ! Gather the number of own points for every PE into p_ri%pe_own

    ALLOCATE(p_ri%pe_own(0:p_n_work-1))
    ALLOCATE(p_ri%pe_off(0:p_n_work-1))
#ifndef NOMPI
    CALL MPI_Allgather(p_ri%n_own,  1, p_int, &
                       p_ri%pe_own, 1, p_int, &
                       p_comm_work, mpierr)
#else
    p_ri%pe_own(0) = p_ri%n_own
#endif
! NOMPI

    ! Get offset within result array
    p_ri%pe_off(0) = 0
    DO i = 1, p_n_work-1
      p_ri%pe_off(i) = p_ri%pe_off(i-1) + p_ri%pe_own(i-1)
    ENDDO

    ! Get global number of points for current (physical!) patch

    p_ri%n_glb = SUM(p_ri%pe_own(:))

    ! Get the global index numbers of the data when it is gathered on PE 0
    ! exactly in the same order as it is retrieved later during I/O

    ALLOCATE(glbidx_glb(p_ri%n_glb))
#ifndef NOMPI
    CALL MPI_Allgatherv(glbidx_own, p_ri%n_own, p_int, &
                        glbidx_glb, p_ri%pe_own, p_ri%pe_off, p_int, &
                        p_comm_work, mpierr)
#else
    glbidx_glb(:) = glbidx_own(:)
#endif
! NOMPI

    ! Get reorder_index

    ALLOCATE(p_ri%reorder_index(p_ri%n_glb))
    IF (patch_info(phys_patch_id)%grid_info_mode == GRID_INFO_FILE) THEN
      ALLOCATE(p_ri%grid_info%log_dom_index(p_ri%n_glb))
    END IF
    ALLOCATE(reorder_index_log_dom(n_points_g)) ! spans the complete logical domain
    reorder_index_log_dom(:) = 0

    DO i = 1, p_ri%n_glb
      ! reorder_index_log_dom stores where a global point in logical domain comes from.
      ! It is nonzero only at the physical patch locations
      reorder_index_log_dom(glbidx_glb(i)) = i
    ENDDO

    ! Gather the reorder index for the physical domain
    n = 0
    DO i = 1, n_points_g
      IF(reorder_index_log_dom(i)>0) THEN
        n = n+1
        p_ri%reorder_index(reorder_index_log_dom(i)) = n
      ENDIF
    ENDDO

    IF (patch_info(phys_patch_id)%grid_info_mode == GRID_INFO_FILE) THEN
      n = 0
      DO i = 1, n_points_g
        IF(reorder_index_log_dom(i)>0) THEN
          n = n+1
          p_ri%grid_info%log_dom_index(n) = i
        ENDIF
      ENDDO
    END IF

    ! Safety check
    IF(n/=p_ri%n_glb) CALL finish(routine,'Reordering failed')

    ! set trivial destination indices:
    IF(my_process_is_mpi_seq()) THEN
      ALLOCATE(p_ri%own_dst_idx(p_ri%n_own), &
        &      p_ri%own_dst_blk(p_ri%n_own))
      DO i=1,p_ri%n_own
        p_ri%own_dst_idx(i) = idx_no(i)
        p_ri%own_dst_blk(i) = blk_no(i)
      END DO ! i
    END IF

    DEALLOCATE(phys_owner_mask)
    DEALLOCATE(glbidx_own)
    DEALLOCATE(glbidx_glb)
    DEALLOCATE(reorder_index_log_dom)

  END SUBROUTINE set_reorder_info


  !------------------------------------------------------------------------------------------------
  !> Sets the reorder_info for lon-lat-grids
  !
#ifndef __NO_ICON_ATMO__
  SUBROUTINE set_reorder_info_lonlat(grid, intp, patch_info_ll)
    TYPE(t_lon_lat_grid),  INTENT(IN)    :: grid
    TYPE(t_lon_lat_intp),  INTENT(IN)    :: intp
    TYPE(t_patch_info_ll), INTENT(INOUT) :: patch_info_ll      ! Result: reorder info

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_reorder_info_lonlat"
    INTEGER :: ierrstat, i, jc, jb, this_pe, mpierr, &
    &          ioffset, gidx, n_own

    ! Just for safety
    IF(my_process_is_io()) CALL finish(routine, 'Must not be called on IO PEs')
    this_pe = get_my_mpi_work_id()

    n_own                  = intp%nthis_local_pts        ! No. of own points
    patch_info_ll%ri%n_glb = grid%lon_dim * grid%lat_dim ! Total points in lon-lat grid
    patch_info_ll%ri%n_own = n_own
    ! Set index arrays to own cells/edges/verts
    ALLOCATE(patch_info_ll%ri%own_idx(n_own), &
      &      patch_info_ll%ri%own_blk(n_own), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    jc = 0
    jb = 1
    DO i=1,n_own
      jc = jc+1
      IF (jc>nproma) THEN
        jb = jb+1;  jc = 1
      END IF
      patch_info_ll%ri%own_idx(i) = jc
      patch_info_ll%ri%own_blk(i) = jb
    END DO ! i

    ! set destination indices (for sequential/test PEs). This is
    ! important for the case that the local patch is smaller than the
    ! lon-lat grid:
    IF(my_process_is_mpi_seq()) THEN
      ALLOCATE(patch_info_ll%ri%own_dst_idx(n_own), &
        &      patch_info_ll%ri%own_dst_blk(n_own))
      DO i=1,n_own
        gidx = intp%global_idx(i)
        patch_info_ll%ri%own_dst_idx(i) = idx_no(gidx)
        patch_info_ll%ri%own_dst_blk(i) = blk_no(gidx)
      END DO ! i
    END IF

    ! Gather the number of own points for every PE into p_ri%pe_own
    ALLOCATE(patch_info_ll%ri%pe_own(0:p_n_work-1), &
      &      patch_info_ll%ri%pe_off(0:p_n_work-1), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
#ifndef NOMPI
    CALL MPI_Allgather(n_own,  1, p_int,                  &
                       patch_info_ll%ri%pe_own, 1, p_int, &
                       p_comm_work, mpierr)
#else
    patch_info_ll%ri%pe_own(0) = n_own
#endif
! NOMPI

    ! Get offset within result array
    patch_info_ll%ri%pe_off(0) = 0
    DO i = 1, p_n_work-1
      patch_info_ll%ri%pe_off(i) = patch_info_ll%ri%pe_off(i-1) + patch_info_ll%ri%pe_own(i-1)
    ENDDO

    ! Get the global index numbers of the data when it is gathered on PE 0
    ! exactly in the same order as it is retrieved later during I/O
    ALLOCATE(patch_info_ll%ri%reorder_index(patch_info_ll%ri%n_glb), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    IF (patch_info_ll%grid_info_mode == GRID_INFO_FILE) THEN
      ALLOCATE(patch_info_ll%ri%grid_info%log_dom_index(patch_info_ll%ri%n_glb), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    END IF

    ioffset = patch_info_ll%ri%pe_off(this_pe)
    patch_info_ll%ri%reorder_index = -1
    DO i=1,intp%nthis_local_pts
      patch_info_ll%ri%reorder_index(ioffset + i) = intp%global_idx(i)
    END DO
    ! merge all fields across working PEs:
    patch_info_ll%ri%reorder_index = p_max(patch_info_ll%ri%reorder_index, &
      &                                    comm=get_my_mpi_work_communicator())
    IF (patch_info_ll%grid_info_mode == GRID_INFO_FILE) THEN
      ! mapping between logical and physical patch is trivial for
      ! lon-lat grids:
      patch_info_ll%ri%grid_info%log_dom_index(:) = (/ (i, i=1,patch_info_ll%ri%n_glb) /)
    END IF
  END SUBROUTINE set_reorder_info_lonlat
#endif

#ifdef USE_CRAY_POINTER
  !------------------------------------------------------------------------------------------------
  ! Helper routines for setting mem_ptr with the correct size information

  SUBROUTINE set_mem_ptr_sp(arr, len)
    INTEGER          :: len
    REAL(sp), TARGET :: arr(len)
    mem_ptr_sp => arr
  END SUBROUTINE set_mem_ptr_sp
  !------------------------------------------------------------------------------------------------
  SUBROUTINE set_mem_ptr_dp(arr, len)
    INTEGER          :: len
    REAL(dp), TARGET :: arr(len)
    mem_ptr_dp => arr
  END SUBROUTINE set_mem_ptr_dp
#endif
! USE_CRAY_POINTER

  !------------------------------------------------------------------------------------------------
  !> Sets up the vlist for a t_output_file structure
  !
  SUBROUTINE setup_output_vlist(of)
    TYPE(t_output_file), INTENT(INOUT) :: of
    ! local variables
    CHARACTER(LEN=*), PARAMETER     :: routine = modname//"::setup_output_vlist"
    INTEGER                         :: k, i_dom, ll_dim(2), gridtype, idate, itime
    TYPE(t_lon_lat_data), POINTER   :: lonlat
    TYPE(t_datetime)                :: ini_datetime
    CHARACTER(len=1)                :: uuid_string(16)
    REAL(wp)                        :: pi_180
    INTEGER                         :: max_cell_connectivity
    REAL(wp), ALLOCATABLE           :: p_lonlat(:)

    pi_180 = ATAN(1._wp)/45._wp

    IF (of%output_type == FILETYPE_GRB2) THEN
      ! since the current CDI-version does not fully support "GRID_UNSTRUCTURED", the
      ! grid type is changed to "GRID_REFERENCE".
      gridtype = GRID_REFERENCE
    ELSE
      gridtype = GRID_UNSTRUCTURED
    ENDIF

    i_dom = of%phys_patch_id
    max_cell_connectivity = patch_info(i_dom)%max_cell_connectivity

    !
    ! The following sections add the file global properties collected in init_name_list_output
    !
    ! 1. create cdi vlist
    !
    of%cdiVlistID = vlistCreate()
    !
    ! 2. add global attributes for netCDF
    !
    ! define output generating institute
    !
    ! inquire the Institute ID from (center/subcenter)
    !
    of%cdiInstID = institutInq(gribout_config(i_dom)%generatingCenter,          &
      &                        gribout_config(i_dom)%generatingSubcenter, '', '')


    ! define Institute
    CALL vlistDefInstitut(of%cdiVlistID,of%cdiInstID)



    ! 3. add horizontal grid descriptions

    IF(of%name_list%remap == REMAP_REGULAR_LATLON) THEN
#ifndef __NO_ICON_ATMO__

      ! Lon/Lat Interpolation requested

      of%cdiCellGridID = CDI_UNDEFID
      of%cdiEdgeGridID = CDI_UNDEFID
      of%cdiVertGridID = CDI_UNDEFID

      lonlat => lonlat_grid_list(of%name_list%lonlat_id)
      ll_dim(1) = lonlat%grid%lon_dim
      ll_dim(2) = lonlat%grid%lat_dim

      of%cdiLonLatGridID = gridCreate(GRID_LONLAT, ll_dim(1)*ll_dim(2))

      CALL gridDefXsize(of%cdiLonLatGridID, ll_dim(1))
      CALL gridDefXname(of%cdiLonLatGridID, 'lon')
      CALL gridDefXunits(of%cdiLonLatGridID, 'degrees_east')

      CALL gridDefYsize(of%cdiLonLatGridID, ll_dim(2))
      CALL gridDefYname(of%cdiLonLatGridID, 'lat')
      CALL gridDefYunits(of%cdiLonLatGridID, 'degrees_north')

      ALLOCATE(p_lonlat(ll_dim(1)))
      IF (lonlat%grid%reg_lon_def(2) <= threshold_delta_or_intvls) THEN
        DO k=1,ll_dim(1)
          p_lonlat(k) = lonlat%grid%reg_lon_def(1) + REAL(k-1,wp)*lonlat%grid%reg_lon_def(2)
        END DO
      ELSE
        DO k=1,ll_dim(1)
          p_lonlat(k) = (lonlat%grid%start_corner(1) + REAL(k-1,wp)*lonlat%grid%delta(1)) / pi_180
        END DO
      END IF
      CALL gridDefXvals(of%cdiLonLatGridID, p_lonlat)
      DEALLOCATE(p_lonlat)

      ALLOCATE(p_lonlat(ll_dim(2)))
      IF (lonlat%grid%reg_lat_def(2) <= threshold_delta_or_intvls) THEN
        DO k=1,ll_dim(2)
          p_lonlat(k) = lonlat%grid%reg_lat_def(1) + REAL(k-1,wp)*lonlat%grid%reg_lat_def(2)
        END DO
      ELSE
        DO k=1,ll_dim(2)
          p_lonlat(k) = (lonlat%grid%start_corner(2) + REAL(k-1,wp)*lonlat%grid%delta(2)) / pi_180
        END DO
      END IF
      CALL gridDefYvals(of%cdiLonLatGridID, p_lonlat)
      DEALLOCATE(p_lonlat)
#endif
! #ifndef __NO_ICON_ATMO__
    ELSE

      ! Cells

      of%cdiCellGridID = gridCreate(gridtype, patch_info(i_dom)%cells%n_glb)
      CALL gridDefNvertex(of%cdiCellGridID, max_cell_connectivity)
      !
      CALL gridDefXname(of%cdiCellGridID, 'clon')
      CALL gridDefXlongname(of%cdiCellGridID, 'center longitude')
      CALL gridDefXunits(of%cdiCellGridID, 'radian')
      !
      CALL gridDefYname(of%cdiCellGridID, 'clat')
      CALL gridDefYlongname(of%cdiCellGridID, 'center latitude')
      CALL gridDefYunits(of%cdiCellGridID, 'radian')
      !
      CALL uuid2char(patch_info(i_dom)%grid_uuid, uuid_string)
      CALL gridDefUUID(of%cdiCellGridID, uuid_string)
      !
      CALL gridDefNumber(of%cdiCellGridID, patch_info(i_dom)%number_of_grid_used)

      !
      ! not clear whether meta-info GRID_CELL or GRID_UNSTRUCTURED_CELL should be used
      CALL gridDefPosition(of%cdiCellGridID, GRID_CELL)

      ! Single point grid for monitoring
      of%cdiSingleGridID = gridCreate(GRID_LONLAT, 1)
      !
      CALL griddefxsize(of%cdiSingleGridID, 1)                                                                         
      CALL griddefysize(of%cdiSingleGridID, 1)
      CALL griddefxvals(of%cdiSingleGridID, (/0.0_wp/))
      CALL griddefyvals(of%cdiSingleGridID, (/0.0_wp/))


      ! Verts

      of%cdiVertGridID = gridCreate(gridtype, patch_info(i_dom)%verts%n_glb)
      CALL gridDefNvertex(of%cdiVertGridID, 9-max_cell_connectivity)
      !
      CALL gridDefXname(of%cdiVertGridID, 'vlon')
      CALL gridDefXlongname(of%cdiVertGridID, 'vertex longitude')
      CALL gridDefXunits(of%cdiVertGridID, 'radian')
      !
      CALL gridDefYname(of%cdiVertGridID, 'vlat')
      CALL gridDefYlongname(of%cdiVertGridID, 'vertex latitude')
      CALL gridDefYunits(of%cdiVertGridID, 'radian')
      !
      CALL uuid2char(patch_info(i_dom)%grid_uuid, uuid_string)
      CALL gridDefUUID(of%cdiVertGridID, uuid_string)
      !
      CALL gridDefNumber(of%cdiVertGridID, patch_info(i_dom)%number_of_grid_used)

      !
      ! not clear whether meta-info GRID_VERTEX or GRID_UNSTRUCTURED_VERTEX should be used
      CALL gridDefPosition(of%cdiVertGridID, GRID_VERTEX)

      ! Edges

      of%cdiEdgeGridID = gridCreate(gridtype, patch_info(i_dom)%edges%n_glb)
      CALL gridDefNvertex(of%cdiEdgeGridID, 4)
      !
      CALL gridDefXname(of%cdiEdgeGridID, 'elon')
      CALL gridDefXlongname(of%cdiEdgeGridID, 'edge midpoint longitude')
      CALL gridDefXunits(of%cdiEdgeGridID, 'radian')
      !
      CALL gridDefYname(of%cdiEdgeGridID, 'elat')
      CALL gridDefYlongname(of%cdiEdgeGridID, 'edge midpoint latitude')
      CALL gridDefYunits(of%cdiEdgeGridID, 'radian')
      !
      CALL uuid2char(patch_info(i_dom)%grid_uuid, uuid_string)
      CALL gridDefUUID(of%cdiEdgeGridID, uuid_string)
      !
      CALL gridDefNumber(of%cdiEdgeGridID, patch_info(i_dom)%number_of_grid_used)

      !
      ! not clear whether meta-info GRID_EDGE or GRID_UNSTRUCTURED_EDGE should be used
      CALL gridDefPosition(of%cdiEdgeGridID, GRID_EDGE)

      of%cdiLonLatGridID = CDI_UNDEFID

      ! If wanted, set grid info
      IF(of%name_list%output_grid) THEN
        SELECT CASE(of%name_list%filetype)
        CASE (FILETYPE_NC2, FILETYPE_NC4)
          ! encode grid info in NetCDF format:
          SELECT CASE(patch_info(of%phys_patch_id)%grid_info_mode)
          CASE (GRID_INFO_FILE)
            CALL copy_grid_info(of, patch_info)
          CASE (GRID_INFO_BCAST) 
            IF (.NOT. my_process_is_mpi_test()) THEN
              CALL set_grid_info_netcdf(of%cdiCellGridID, patch_info(of%phys_patch_id)%cells%grid_info)
              CALL set_grid_info_netcdf(of%cdiEdgeGridID, patch_info(of%phys_patch_id)%edges%grid_info)
              CALL set_grid_info_netcdf(of%cdiVertGridID, patch_info(of%phys_patch_id)%verts%grid_info)
            END IF
          END SELECT
        CASE (FILETYPE_GRB2)
          ! handled later...
        CASE DEFAULT
          CALL finish(routine, "Unknown grid type")
        END SELECT
      END IF

    ENDIF

    !
    ! 4. add vertical grid descriptions

    IF (iequations/=ihs_ocean) THEN ! atm
      SELECT CASE(of%ilev_type)
      CASE (level_type_ml) 
        CALL setup_ml_axes_atmo(of)
      CASE (level_type_pl) 
        CALL setup_pl_axis_atmo(of)
      CASE (level_type_hl) 
        CALL setup_hl_axis_atmo(of)
      CASE (level_type_il) 
        CALL setup_il_axis_atmo(of)
      CASE DEFAULT
        CALL finish(routine, "Internal error!")
      END SELECT
    ELSE
      CALL setup_zaxes_oce(of)
    END IF

    !
    ! 5. output does contain absolute time
    !
    SELECT CASE (of%name_list%mode)
    CASE (1)  ! forecast mode
     of%cdiTaxisID = taxisCreate(TAXIS_RELATIVE)

     IF (of%name_list%taxis_tunit > 10 .OR. of%name_list%taxis_tunit < 1 ) THEN
       of%name_list%taxis_tunit=TUNIT_MINUTE
       CALL message('','invalid taxis_tunit, reset to TUNIT_MINUTE')
     END IF
     CALL taxisDefTunit (of%cdiTaxisID, of%name_list%taxis_tunit)
     ini_datetime = time_config%ini_datetime
     CALL taxisDefCalendar (of%cdiTaxisID, time_config%calendar)
     idate = cdiEncodeDate(ini_datetime%year, ini_datetime%month, ini_datetime%day)
     itime = cdiEncodeTime(ini_datetime%hour, ini_datetime%minute, &
                           NINT(ini_datetime%second))
     !WRITE(6,'(a,i,a)')'calendar ', time_config%calendar, &
     !                & 'julian_gregorian 0 -  proleptic_gregorian 1 -  cly360 2'
     CALL taxisDefRdate (of%cdiTaxisID, idate )
     CALL taxisDefRtime (of%cdiTaxisID, itime )

     !WRITE(6,'(a,i,a,i)')'idate ',idate,' ',taxisInqRdate(of%cdiTaxisID)
     !WRITE(6,'(a,i,a,i)')'itime ',itime,' ',taxisInqRtime(of%cdiTaxisID)
    CASE (2)  ! climate mode
     of%cdiTaxisID = taxisCreate(TAXIS_ABSOLUTE)
    CASE DEFAULT
     of%cdiTaxisID = taxisCreate(TAXIS_ABSOLUTE)
    END SELECT

    !
    CALL vlistDefTaxis(of%cdiVlistID, of%cdiTaxisID)

    !
    ! add variables
    !
    CALL add_variables_to_vlist(of)

    ! GRB2 format: define geographical longitude, latitude as special
    ! variables "RLON", "RLAT". Note that the grid information may be
    ! contained in the patch_info data structure (for a different
    ! output_file) but still we may do not want to include it into
    ! this output file's (of) data file.
    IF ((of%name_list%output_grid)                                      .AND. &
      & (patch_info(of%phys_patch_id)%grid_info_mode /= GRID_INFO_NONE) .AND. &
      & (of%name_list%filetype == FILETYPE_GRB2)) THEN
      CALL set_grid_info_grb2(of)
    END IF

  END SUBROUTINE setup_output_vlist


  !------------------------------------------------------------------------------------------------
  !> define variables and attributes
  !
  SUBROUTINE add_variables_to_vlist(of)
    TYPE (t_output_file), INTENT(IN), TARGET :: of
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::add_variables_to_vlist"
    TYPE (t_var_metadata), POINTER :: info
    INTEGER                        :: iv, vlistID, varID, gridID, &
      &                               zaxisID
    CHARACTER(LEN=DICT_MAX_STRLEN) :: mapped_name
    TYPE(t_cf_var), POINTER        :: this_cf

    vlistID = of%cdiVlistID

    DO iv = 1, of%num_vars
      !
      info => of%var_desc(iv)%info
      !
      ! set grid ID
      !
      SELECT CASE (info%hgrid)
      CASE(GRID_UNSTRUCTURED_CELL)
        info%cdiGridID = of%cdiCellGridID
      CASE(GRID_LONLAT)
        info%cdiGridID = of%cdiSingleGridID
      CASE(GRID_UNSTRUCTURED_VERT)
        info%cdiGridID = of%cdiVertGridID
      CASE(GRID_UNSTRUCTURED_EDGE)
        info%cdiGridID = of%cdiEdgeGridID
      CASE (GRID_REGULAR_LONLAT)
        info%cdiGridID = of%cdiLonLatGridID
      CASE DEFAULT
        CALL finish(routine, 'GRID definition missing for '//TRIM(info%name))
      END SELECT

      gridID = info%cdiGridID



      !
      ! set z axis ID
      !
      IF (info%cdiZaxisID == CDI_UNDEFID) info%cdiZaxisID = of%cdiZaxisID(info%vgrid)
      zaxisID = info%cdiZaxisID
      IF (zaxisID /= CDI_UNDEFID) THEN

!DR *********** FIXME *************
        ! Re-set
        ! ZA_HYBRID       -> ZA_REFERENCE
        ! ZA_HYBRID_HALF  -> ZA_REFERENCE_HALF
        ! as long as ZA_hybrid/ZA_hybrid_half is used throughout the code.
        ! Should be replaced by ZA_reference/ZA_reference_half for the 
        ! nonhydrostatic model.
        IF (zaxisID == of%cdiZaxisID(ZA_hybrid)) THEN
          zaxisID = of%cdiZaxisID(ZA_reference)
        ELSE IF (zaxisID == of%cdiZaxisID(ZA_hybrid_half)) THEN
          zaxisID = of%cdiZaxisID(ZA_reference_half)
        ELSE IF (zaxisID == of%cdiZaxisID(ZA_hybrid_half_hhl)) THEN
          zaxisID = of%cdiZaxisID(ZA_reference_half_hhl)
        ENDIF
        info%cdiZaxisID = zaxisID
!DR*********WILL BE REMOVED SOON**********

      ELSE
        WRITE (message_text,'(a,i3,a,i3)') &
             &  'Zaxis Nr.: ',info%vgrid,' not defined. zaxisID= ',zaxisID
        CALL finish(routine, message_text)
      ENDIF

      ! Search name mapping for name in NetCDF file
      mapped_name = TRIM(dict_get(out_varnames_dict, info%name, default=info%name))

      ! note that an explicit call of vlistDefVarTsteptype is obsolete, since
      ! isteptype is already defined via vlistDefVar
      varID = vlistDefVar(vlistID, gridID, zaxisID, info%isteptype)
      info%cdiVarID   = varID

      CALL vlistDefVarName(vlistID, varID, TRIM(mapped_name))
      IF (info%post_op%lnew_cf) THEN
        this_cf => info%post_op%new_cf
      ELSE
        this_cf => info%cf
      END IF

      IF (this_cf%long_name /= '')     CALL vlistDefVarLongname(vlistID, varID, this_cf%long_name)
      IF (this_cf%standard_name /= '') CALL vlistDefVarStdname(vlistID, varID, this_cf%standard_name)
      IF (this_cf%units /= '')         CALL vlistDefVarUnits(vlistID, varID, this_cf%units)

      ! Currently only real valued variables are allowed, so we can always use info%missval%rval
      IF (info%lmiss) CALL vlistDefVarMissval(vlistID, varID, info%missval%rval)

      ! Set GRIB2 Triplet
      IF (info%post_op%lnew_grib2) THEN
        CALL vlistDefVarParam(vlistID, varID,                   &
          &  cdiEncodeParam(info%post_op%new_grib2%number,      &
          &                 info%post_op%new_grib2%category,    &
          &                 info%post_op%new_grib2%discipline) )
        IF ( of%output_type == FILETYPE_GRB2 ) THEN
          CALL vlistDefVarDatatype(vlistID, varID, info%post_op%new_grib2%bits)
        END IF
      ELSE
        CALL vlistDefVarParam(vlistID, varID,                   &
          &  cdiEncodeParam(info%grib2%number,                  &
          &                 info%grib2%category,                &
          &                 info%grib2%discipline) )
        IF ( of%output_type == FILETYPE_GRB2 ) THEN
          CALL vlistDefVarDatatype(vlistID, varID, info%grib2%bits)
        END IF
      END IF

      IF ( of%output_type == FILETYPE_GRB2 ) THEN

!!$        ! GRIB2 Quick hack: Set additional GRIB2 keys
!!$        CALL set_additional_GRIB2_keys(vlistID, varID, gribout_config(of%phys_patch_id), &
!!$          &                            get_var_tileidx(TRIM(info%name)) )

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!          ATTENTION                    !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Note that re-setting of the surface types must come AFTER (re)setting
        ! "productDefinitionTemplateNumber" in set_additional_GRIB2_keys. It was observed 
        ! (i.e. for Ensemble output), that the surface-type information is lost again, if 
        ! these settings are performed prior to "productDefinitionTemplateNumber"  

        ! Re-Set "typeOfSecondFixedSurface" for the following set of variables
        !
        ! GRIB_CHECK(grib_set_long(gh, "typeOfSecondFixedSurface", xxx), 0);
        !
        ! HHL     : typeOfSecondFixedSurface = 101
        ! HSURF   : typeOfSecondFixedSurface = 101
        ! HBAS_CON: typeOfSecondFixedSurface = 101
        ! HTOP_CON: typeOfSecondFixedSurface = 101
        ! HTOP_DC : typeOfSecondFixedSurface = 101
        ! HZEROCL : typeOfSecondFixedSurface = 101
        ! CLCL    : typeOfSecondFixedSurface = 1
        ! C_T_LK  : typeOfSecondFixedSurface = 162
        ! H_B1_LK : typeOfSecondFixedSurface = 165
        ! SNOWLMT : typeOfSecondFixedSurface = 101
        IF ( one_of(TRIM(info%name),sfs_name_list) /= -1 ) THEN
          CALL vlistDefVarIntKey(vlistID, varID, "typeOfSecondFixedSurface", &
            &                    second_tos(one_of(TRIM(info%name),sfs_name_list)))
        ENDIF

        ! Re-Set "typeOfFirstFixedSurface" for the following set of variables
        !
        ! GRIB_CHECK(grib_set_long(gh, "typeOfFirstFixedSurface", xxx), 0);
        !
        ! T_MNW_LK: typeOfFirstFixedSurface = 1
        ! DEPTH_LK: typeOfFirstFixedSurface = 1
        ! T_WML_LK: typeOfFirstFixedSurface = 1
        ! H_ML_LK : typeOfFirstFixedSurface = 1
        !
        IF ( one_of(TRIM(info%name),ffs_name_list) /= -1 ) THEN
          CALL vlistDefVarIntKey(vlistID, varID, "typeOfFirstFixedSurface", &
            &                    first_tos(one_of(TRIM(info%name),ffs_name_list)))
        ENDIF


        ! Quick hack: shortName.def should be revised, instead
        IF ( TRIM(info%name)=='qv_s' ) THEN
          CALL vlistDefVarIntKey(vlistID, varID, "scaleFactorOfFirstFixedSurface", 0)
        ENDIF


        ! GRIB2 Quick hack: Set additional GRIB2 keys
        CALL set_additional_GRIB2_keys(vlistID, varID, gribout_config(of%phys_patch_id), &
          &                            get_var_tileidx(TRIM(info%name)) )

      ELSE ! NetCDF
        CALL vlistDefVarDatatype(vlistID, varID, this_cf%datatype)
      ENDIF

    ENDDO
    !
  END SUBROUTINE add_variables_to_vlist


  !------------------------------------------------------------------------------------------------
  !> Transfers reorder_info to IO PEs
  !
  SUBROUTINE transfer_reorder_info(p_ri, grid_info_mode)
    TYPE(t_reorder_info), INTENT(INOUT) :: p_ri ! Result: reorder info
    INTEGER,              INTENT(IN)    :: grid_info_mode

    ! There is nothing to do for the test PE:
    IF(my_process_is_mpi_test()) RETURN

    ! Transfer the global number of points, this is not yet known on IO PEs
    CALL p_bcast(p_ri%n_glb,  bcast_root, p_comm_work_2_io)
    CALL p_bcast(p_ri%n_log,  bcast_root, p_comm_work_2_io)

    IF(my_process_is_io()) THEN

      ! On IO PEs: n_own = 0, own_idx and own_blk are not allocated
      p_ri%n_own = 0

      ! pe_own/pe_off must be allocated for num_work_procs, not for p_n_work
      ALLOCATE(p_ri%pe_own(0:num_work_procs-1))
      ALLOCATE(p_ri%pe_off(0:num_work_procs-1))

      ALLOCATE(p_ri%reorder_index(p_ri%n_glb))
      IF (grid_info_mode == GRID_INFO_FILE) THEN
        ALLOCATE(p_ri%grid_info%log_dom_index(p_ri%n_glb))
      END IF
    ENDIF

    CALL p_bcast(p_ri%pe_own, bcast_root, p_comm_work_2_io)
    CALL p_bcast(p_ri%pe_off, bcast_root, p_comm_work_2_io)

    CALL p_bcast(p_ri%reorder_index, bcast_root, p_comm_work_2_io)
    IF (grid_info_mode == GRID_INFO_FILE) THEN
      CALL p_bcast(p_ri%grid_info%log_dom_index, bcast_root, p_comm_work_2_io)
    END IF

  END SUBROUTINE transfer_reorder_info


  !-------------------------------------------------------------------------------------------------
  !> Replicates data (mainly the variable lists) needed for async I/O on the I/O procs.
  !  ATTENTION: The data is not completely replicated, only as far as needed for I/O.
  !
  !  This routine has to be called by all PEs (work and I/O)
  !
  SUBROUTINE replicate_data_on_io_procs()

    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//"::replicate_data_on_io_procs"
    INTEGER                       :: ivct_len
    INTEGER                       :: info_size, iv, nv, nelems, n, list_info(4)
    INTEGER, ALLOCATABLE          :: info_storage(:,:)
    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_metadata)          :: info
    TYPE(t_var_list)              :: p_var_list
    ! var_list_name should have at least the length of var_list names
    ! (although this doesn't matter as long as it is big enough for every name)
    CHARACTER(LEN=256)            :: var_list_name
    INTEGER                       :: idom

!DR Test
    INTEGER :: nvgrid, ivgrid



    ! There is nothing to do for the test PE:
    IF(my_process_is_mpi_test()) RETURN

    !-----------------------------------------------------------------------------------------------

    ! Replicate vertical coordinate table
#ifndef __NO_ICON_ATMO__
    IF(.NOT.my_process_is_io()) ivct_len = SIZE(vct)
    CALL p_bcast(ivct_len, bcast_root, p_comm_work_2_io)

    IF(my_process_is_io()) ALLOCATE(vct(ivct_len))
    CALL p_bcast(vct, bcast_root, p_comm_work_2_io)
#endif
! #ifndef __NO_ICON_ATMO__
    !-----------------------------------------------------------------------------------------------
    ! Replicate variable lists

    ! Get the size - in default INTEGER words - which is needed to
    ! hold the contents of TYPE(t_var_metadata)

    info_size = SIZE(TRANSFER(info, (/ 0 /)))

    ! Get the number of var_lists
    IF(.NOT.my_process_is_io()) nv = nvar_lists
    CALL p_bcast(nv, bcast_root, p_comm_work_2_io)

    ! For each var list, get its components
    DO iv = 1, nv

      ! Send name
      IF(.NOT.my_process_is_io()) var_list_name = var_lists(iv)%p%name
      CALL p_bcast(var_list_name, bcast_root, p_comm_work_2_io)

      IF(.NOT.my_process_is_io()) THEN

        ! Count the number of variable entries
        element => var_lists(iv)%p%first_list_element
        nelems = 0
        DO
          IF(.NOT.ASSOCIATED(element)) EXIT
          nelems = nelems+1
          element => element%next_list_element
        ENDDO

        ! Gather the components needed for name list I/O and send them.
        ! Please note that not the complete list is replicated, unneeded
        ! entries are left away!

        list_info(1) = nelems
        list_info(2) = var_lists(iv)%p%patch_id
        list_info(3) = var_lists(iv)%p%vlevel_type
        list_info(4) = MERGE(1,0,var_lists(iv)%p%loutput)

      ENDIF

      ! Send basic info:

      CALL p_bcast(list_info, bcast_root, p_comm_work_2_io)

      IF(my_process_is_io()) THEN
        nelems = list_info(1)
        ! Create var list
        CALL new_var_list( p_var_list, var_list_name, patch_id=list_info(2), &
                           vlevel_type=list_info(3), loutput=(list_info(4)/=0) )
      ENDIF

      ! Get the binary representation of all info members of the variables
      ! of the list and send it to the receiver.
      ! Using the Fortran TRANSFER intrinsic may seem like a hack,
      ! but it has the advantage that it is completely independet of the
      ! actual declaration if TYPE(t_var_metadata).
      ! Thus members may added to or removed from TYPE(t_var_metadata)
      ! without affecting the code below and we don't have an additional
      ! cross dependency between TYPE(t_var_metadata) and this module.

      ALLOCATE(info_storage(info_size, nelems))

      IF(.NOT.my_process_is_io()) THEN
        element => var_lists(iv)%p%first_list_element
        nelems = 0
        DO
          IF(.NOT.ASSOCIATED(element)) EXIT
          nelems = nelems+1
          info_storage(:,nelems) = TRANSFER(element%field%info, (/ 0 /))
          element => element%next_list_element
        ENDDO
      ENDIF

      ! Send binary representation of all info members

      CALL p_bcast(info_storage, bcast_root, p_comm_work_2_io)

      IF(my_process_is_io()) THEN

        ! Insert elements into var list

        p_var_list%p%first_list_element => NULL()
        element => NULL() ! Safety only

        DO n = 1, nelems
          IF(.NOT.ASSOCIATED(p_var_list%p%first_list_element)) THEN
            ALLOCATE(p_var_list%p%first_list_element)
            element => p_var_list%p%first_list_element
          ELSE
            ALLOCATE(element%next_list_element)
            element => element%next_list_element
          ENDIF

          element%next_list_element => NULL()

          ! Nullify all pointers in element%field, they don't make sense on the I/O PEs

          element%field%r_ptr => NULL()
          element%field%i_ptr => NULL()
          element%field%l_ptr => NULL()
          element%field%var_base_size = 0 ! Unknown here

          ! Set info structure from binary representation in info_storage
          element%field%info = TRANSFER(info_storage(:, n), info)
        ENDDO

      ENDIF

      DEALLOCATE(info_storage)

    ENDDO

    ! Map the variable groups given in the output namelist onto the
    ! corresponding variable subsets:
    CALL parse_variable_groups()


    ! Go over all output domains
    DO idom = 1, n_dom_out
      CALL p_bcast(gribout_config(idom)%generatingCenter,    bcast_root, p_comm_work_2_io)
      CALL p_bcast(gribout_config(idom)%generatingSubcenter, bcast_root, p_comm_work_2_io)
    ENDDO

    ! allocate vgrid_buffer on asynchronous output PEs, for storing 
    ! the vertical grid UUID
    !
    ! get buffer size and broadcast
    IF (ALLOCATED(vgrid_buffer)) THEN
       nvgrid = SIZE(vgrid_buffer)
    ELSE
       nvgrid = 0
    END IF
    CALL p_bcast(nvgrid, bcast_root, p_comm_work_2_io)
    !
    ! allocate on asynchronous PEs
    IF(my_process_is_io()) THEN
      ALLOCATE(vgrid_buffer(nvgrid))
    ENDIF
    ! broadcast
    DO ivgrid = 1,nvgrid 
      CALL p_bcast(vgrid_buffer(ivgrid)%uuid, bcast_root, p_comm_work_2_io)
    ENDDO

    !-----------------------------------------------------------------------------------------------
    ! Replicate coordinates of cells/edges/vertices:

    ! Go over all output domains
    DO idom = 1, n_dom_out
      CALL p_bcast(patch_info(idom)%nblks_glb_c, bcast_root, p_comm_work_2_io)
      CALL p_bcast(patch_info(idom)%nblks_glb_e, bcast_root, p_comm_work_2_io)
      CALL p_bcast(patch_info(idom)%nblks_glb_v, bcast_root, p_comm_work_2_io)
      CALL p_bcast(patch_info(idom)%max_cell_connectivity, bcast_root, p_comm_work_2_io)

      IF (patch_info(idom)%grid_info_mode == GRID_INFO_BCAST) THEN
        CALL bcast_grid_info(patch_info(idom), bcast_root)
      END IF
    END DO

  END SUBROUTINE replicate_data_on_io_procs


#ifndef NOMPI

  !------------------------------------------------------------------------------------------------
  !> Initializes the memory window for asynchronous IO
  !
  SUBROUTINE init_memory_window

#ifdef __SUNPRO_F95
    INCLUDE "mpif.h"
#else
    USE mpi, ONLY: MPI_ADDRESS_KIND
#endif
! __SUNPRO_F95

    ! local variables
    CHARACTER(LEN=*), PARAMETER     :: routine = modname//"::init_async_name_list_output"
    
    INTEGER                         :: jp, i, iv, nlevs, i_log_dom, &
      &                                n_own, lonlat_id
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_size


    ! There is nothing to do for the test PE:
    IF(my_process_is_mpi_test()) RETURN

    ! Go over all output files
    OUT_FILE_LOOP : DO i = 1, SIZE(output_file)
      
      ! Get size of the data for every output file
      mem_size = 0_i8
      
      ! Go over all name list variables for this output file
      DO iv = 1, output_file(i)%num_vars

        jp = output_file(i)%phys_patch_id

        IF(output_file(i)%var_desc(iv)%info%ndims == 2) THEN
          nlevs = 1
        ELSE
          ! handle the case that a few levels have been selected out of
          ! the total number of levels:
          IF (ASSOCIATED(output_file(i)%level_selection)) THEN
            nlevs = output_file(i)%level_selection%n_selected
          ELSE
            nlevs = output_file(i)%var_desc(iv)%info%used_dimensions(2)
          END IF
        ENDIF

        SELECT CASE (output_file(i)%var_desc(iv)%info%hgrid)
        CASE (GRID_UNSTRUCTURED_CELL)
          mem_size = mem_size + INT(nlevs*patch_info(jp)%cells%n_own,i8)
        CASE (GRID_UNSTRUCTURED_EDGE)
          mem_size = mem_size + INT(nlevs*patch_info(jp)%edges%n_own,i8)
        CASE (GRID_UNSTRUCTURED_VERT)
          mem_size = mem_size + INT(nlevs*patch_info(jp)%verts%n_own,i8)

#ifndef __NO_ICON_ATMO__
        CASE (GRID_REGULAR_LONLAT)
          lonlat_id = output_file(i)%var_desc(iv)%info%hor_interp%lonlat_id
          i_log_dom = output_file(i)%log_patch_id
          n_own     = lonlat_info(lonlat_id, i_log_dom)%ri%n_own
          mem_size  = mem_size + INT(nlevs*n_own,i8)
#endif
! #ifndef __NO_ICON_ATMO__

        CASE DEFAULT
          CALL finish(routine,'unknown grid type')
        END SELECT

      ENDDO ! vars
     
      ! allocate amount of memory needed with MPI_Alloc_mem
#ifdef USE_CRAY_POINTER
      CALL allocate_mem_cray(mem_size, output_file(i))
#else
      CALL allocate_mem_noncray(mem_size, output_file(i))
#endif
      ! USE_CRAY_POINTER

      ! allocate memory window for meta-info communication between
      ! PE#0 and the I/O PEs:
      CALL metainfo_allocate_memory_window(output_file(i)%mem_win, output_file(i)%num_vars)

    ENDDO OUT_FILE_LOOP
    
  END SUBROUTINE init_memory_window
  

#ifdef USE_CRAY_POINTER
  !------------------------------------------------------------------------------------------------
  !> allocate amount of memory needed with MPI_Alloc_mem
  !
  !  @note Implementation for Cray pointers
  !
  SUBROUTINE allocate_mem_cray(mem_size, of)
#ifdef __SUNPRO_F95
    INCLUDE "mpif.h"
#else
    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_INFO_NULL
#endif
! __SUNPRO_F95

    INTEGER (KIND=MPI_ADDRESS_KIND), INTENT(IN)    :: mem_size
    TYPE (t_output_file),            INTENT(INOUT) :: of
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allocate_mem_cray"
    INTEGER (KIND=MPI_ADDRESS_KIND) :: iptr
    REAL(sp)                        :: tmp_sp
    REAL(dp)                        :: tmp_dp
    INTEGER                         :: mpierr
    INTEGER                         :: nbytes_real
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_bytes
    POINTER(tmp_ptr_sp,tmp_sp(*))
    POINTER(tmp_ptr_dp,tmp_dp(*))

    ! Get the amount of bytes per REAL*8 or REAL*4 variable (as used in MPI
    ! communication)
    IF (use_dp_mpi2io) THEN
      CALL MPI_Type_extent(p_real_dp, nbytes_real, mpierr)
    ELSE
      CALL MPI_Type_extent(p_real_sp, nbytes_real, mpierr)
    ENDIF

    ! For the IO PEs the amount of memory needed is 0 - allocate at least 1 word there:
    mem_bytes = MAX(mem_size,1_i8)*INT(nbytes_real,i8)

    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, iptr, mpierr)

    IF (use_dp_mpi2io) THEN
      tmp_ptr_dp = iptr
      CALL set_mem_ptr_dp(tmp_dp, INT(mem_size))
    ELSE
      tmp_ptr_sp = iptr
      CALL set_mem_ptr_sp(tmp_sp, INT(mem_size))
    ENDIF

    ! Create memory window for communication
    IF (use_dp_mpi2io) THEN
      of%mem_win%mem_ptr_dp(:) = 0._dp
      CALL MPI_Win_create( of%mem_win%mem_ptr_dp,mem_bytes,nbytes_real,MPI_INFO_NULL,&
        &                  p_comm_work_io,of%mem_win%mpi_win,mpierr )
    ELSE
      of%mem_win%mem_ptr_sp(:) = 0._sp
      CALL MPI_Win_create( of%mem_win%mem_ptr_sp,mem_bytes,nbytes_real,MPI_INFO_NULL,&
        &                  p_comm_work_io,of%mem_win%mpi_win,mpierr )
    ENDIF
    IF (mpierr /= 0) CALL finish(TRIM(routine), "MPI error!")

  END SUBROUTINE allocate_mem_cray
#endif
! USE_CRAY_POINTER


#ifndef USE_CRAY_POINTER
  !------------------------------------------------------------------------------------------------
  !> allocate amount of memory needed with MPI_Alloc_mem
  !
  !  @note Implementation for non-Cray pointers
  !
  SUBROUTINE allocate_mem_noncray(mem_size, of)
#ifdef __SUNPRO_F95
    INCLUDE "mpif.h"
#else
    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_INFO_NULL
#endif
! __SUNPRO_F95

    INTEGER (KIND=MPI_ADDRESS_KIND), INTENT(IN)    :: mem_size
    TYPE (t_output_file),            INTENT(INOUT) :: of
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allocate_mem_noncray"
    TYPE(c_ptr)                     :: c_mem_ptr
    INTEGER                         :: mpierr
    INTEGER                         :: nbytes_real
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_bytes

    ! Get the amount of bytes per REAL*8 or REAL*4 variable (as used in MPI
    ! communication)
    IF (use_dp_mpi2io) THEN
      CALL MPI_Type_extent(p_real_dp, nbytes_real, mpierr)
    ELSE
      CALL MPI_Type_extent(p_real_sp, nbytes_real, mpierr)
    ENDIF

    ! For the IO PEs the amount of memory needed is 0 - allocate at least 1 word there:
    mem_bytes = MAX(mem_size,1_i8)*INT(nbytes_real,i8)

    ! TYPE(c_ptr) and INTEGER(KIND=MPI_ADDRESS_KIND) do NOT necessarily have the same size!!!
    ! So check if at least c_intptr_t and MPI_ADDRESS_KIND are the same, else we may get
    ! into deep, deep troubles!
    ! There is still a slight probability that TYPE(c_ptr) does not have the size indicated
    ! by c_intptr_t since the standard only requires c_intptr_t is big enough to hold pointers
    ! (so it may be bigger than a pointer), but I hope no vendor screws up its ISO_C_BINDING
    ! in such a way!!!
    ! If c_intptr_t<=0, this type is not defined and we can't do this check, of course.

    IF(c_intptr_t > 0 .AND. c_intptr_t /= MPI_ADDRESS_KIND) &
     & CALL finish(routine,'c_intptr_t /= MPI_ADDRESS_KIND, too dangerous to proceed!')

    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, c_mem_ptr, mpierr)

    ! The NEC requires a standard INTEGER array as 3rd argument for c_f_pointer,
    ! although it would make more sense to have it of size MPI_ADDRESS_KIND.

    NULLIFY(of%mem_win%mem_ptr_sp)
    NULLIFY(of%mem_win%mem_ptr_dp)

    IF (use_dp_mpi2io) THEN

#ifdef __SX__
      CALL C_F_POINTER(c_mem_ptr, of%mem_win%mem_ptr_dp, (/ INT(mem_size) /) )
#else
      CALL C_F_POINTER(c_mem_ptr, of%mem_win%mem_ptr_dp, (/ mem_size /) )
#endif
      ! Create memory window for communication
      of%mem_win%mem_ptr_dp(:) = 0._dp
      CALL MPI_Win_create( of%mem_win%mem_ptr_dp,mem_bytes,nbytes_real,MPI_INFO_NULL,&
        &                  p_comm_work_io,of%mem_win%mpi_win,mpierr )
      IF (mpierr /= 0) CALL finish(TRIM(routine), "MPI error!")

    ELSE

#ifdef __SX__
      CALL C_F_POINTER(c_mem_ptr, of%mem_win%mem_ptr_sp, (/ INT(mem_size) /) )
#else
      CALL C_F_POINTER(c_mem_ptr, of%mem_win%mem_ptr_sp, (/ mem_size /) )
#endif
      ! Create memory window for communication
      of%mem_win%mem_ptr_sp(:) = 0._sp
      CALL MPI_Win_create( of%mem_win%mem_ptr_sp,mem_bytes,nbytes_real,MPI_INFO_NULL,&
        &                  p_comm_work_io,of%mem_win%mpi_win,mpierr )
      IF (mpierr /= 0) CALL finish(TRIM(routine), "MPI error!")

    ENDIF ! use_dp_mpi2io

  END SUBROUTINE allocate_mem_noncray
#endif
! .not. USE_CRAY_POINTER

#endif
! NOMPI

END MODULE mo_name_list_output_init
