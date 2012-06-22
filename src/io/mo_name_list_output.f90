! Define USE_CRAY_POINTER for platforms having problems with ISO_C_BINDING
! BUT understand CRAY pointers

!   #define USE_CRAY_POINTER

!>
!! Module handling synchronous and asynchronous output; supporting
!! multiple I/O PEs and horizontal interpolation.
!!
!! @author R. Johanni
!!
!! @par Revision History
!! Initial implementation  by  R. Johanni  (2011)
!! Removed interpolation on I/O PEs : F. Prill, DWD (2012-03-28)
!!
!! @todo In asynchronous I/O mode, windows are created but not freed
!!
!! @todo Several fields are allocated but not freed at the end of the
!!       simulation. A pseudo-destructor should be implemented!
!!
MODULE mo_name_list_output

  ! Please note: The spelling "name_list" (with underscore) is intented to make
  ! clear that this does not pertain to a FORTRAN namelist but rather
  ! to a list of names of output variables

#ifndef USE_CRAY_POINTER
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_intptr_t, c_f_pointer
#endif

  USE mo_kind,                  ONLY: wp, i8, dp, sp
  USE mo_impl_constants,        ONLY: max_phys_dom, ihs_ocean, zml_soil, MAX_NVARS,   &
    &                                 vname_len, max_dom, SUCCESS, HINTP_TYPE_LONLAT, &
    &                                 min_rlcell_int, min_rledge_int, min_rlvert
  USE mo_grid_config,           ONLY: n_dom, n_phys_dom, global_cell_type, &
    &                                 grid_rescale_factor, start_time, end_time
  USE mo_grid_levels,           ONLY: check_orientation
  USE mo_cdi_constants          ! We need all
  USE mo_io_units,              ONLY: filename_max, nnml, nnml_output, find_next_free_unit
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_var_metadata,          ONLY: t_var_metadata, VARNAME_LEN
  USE mo_linked_list,           ONLY: t_var_list, t_list_element
  USE mo_var_list,              ONLY: nvar_lists, max_var_lists, var_lists,     &
    &                                 new_var_list, get_all_var_names,          &
    &                                 total_number_of_variables, collect_group
  USE mo_var_list_element,      ONLY: level_type_ml, level_type_pl, level_type_hl
  USE mo_util_uuid,             ONLY: t_uuid
  ! MPI Communication routines
  USE mo_mpi,                   ONLY: p_send, p_recv, p_bcast, p_barrier, p_stop, &
    &                                 get_my_mpi_work_id, p_max, get_my_mpi_work_communicator
  ! MPI Communicators
  USE mo_mpi,                   ONLY: p_comm_work, p_comm_work_io, p_comm_work_2_io
  ! MPI Data types
  USE mo_mpi,                   ONLY: p_int, p_int_i8, &
    &                                 p_real_dp, p_real_sp
  ! MPI Process type intrinsics
  USE mo_mpi,                   ONLY: my_process_is_stdio, my_process_is_mpi_test,          &
                                      my_process_is_mpi_workroot, my_process_is_mpi_seq,    &
                                      my_process_is_io
  ! MPI Process IDs
  USE mo_mpi,                   ONLY: process_mpi_all_test_id, process_mpi_all_workroot_id, &
                                      process_mpi_stdio_id
  ! MPI Process group sizes
  USE mo_mpi,                   ONLY: process_mpi_io_size, num_work_procs, p_n_work
  ! Processor numbers
  USE mo_mpi,                   ONLY: p_pe, p_pe_work, p_work_pe0, p_io_pe0
  USE mo_model_domain,          ONLY: t_patch, p_patch, p_phys_patch
  USE mo_parallel_config,       ONLY: nproma, p_test_run, use_sp_output
  USE mo_vertical_coord_table,  ONLY: vct
  USE mo_dynamics_config,       ONLY: iequations, nnow, nnow_rcf
  USE mo_run_config,            ONLY: num_lev, num_levp1, dtime, ldump_states, ldump_dd, &
    &                                 msg_level, output_mode, ltestcase
  USE mo_nh_pzlev_config,       ONLY: nh_pzlev_config
  USE mo_lnd_nwp_config,        ONLY: nlev_snow
  USE mo_datetime,              ONLY: t_datetime
  USE mo_time_config,           ONLY: time_config
  USE mo_lonlat_grid,           ONLY: t_lon_lat_grid
  USE mo_intp_data_strc,        ONLY: t_lon_lat_intp,                         &
    &                                 t_lon_lat_data, get_free_lonlat_grid,   &
    &                                 lonlat_grid_list, n_lonlat_grids
  USE mo_master_nml,            ONLY: model_base_dir
  USE mo_ocean_nml,             ONLY: n_zlev
  USE mo_oce_state,             ONLY: set_zlev

  USE mo_util_string,           ONLY: toupper, t_keyword_list, associate_keyword, &
    &                                 with_keywords, insert_group, MAX_STRING_LEN
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_communication,         ONLY: exchange_data, t_comm_pattern, idx_no, blk_no
  USE mo_math_utilities,        ONLY: t_geographical_coordinates
  USE mo_math_constants,        ONLY: pi
  USE mo_name_list_output_config, ONLY: name_list_output_active, &
  &                                     use_async_name_list_io,  &
  &                                     l_output_phys_patch,     &
  &                                     max_var_ml,              &
  &                                     max_var_pl,              &
  &                                     max_var_hl,              &
  &                                     max_bounds,              &
  &                                     max_levels,              &
  &                                     vname_len,               &  
  &                                     t_output_name_list,      &
  &                                     first_output_name_list,  &
  &                                     max_time_levels,         &
  &                                     t_output_file,           &
  &                                     is_output_nml_active,    &
  &                                     is_output_file_active,   &
  &                                     t_var_desc,              &
  &                                     add_var_desc
  ! meteogram output
  USE mo_meteogram_output,    ONLY: meteogram_init, meteogram_finalize, meteogram_flush_file
  USE mo_meteogram_config,    ONLY: meteogram_output_config


  IMPLICIT NONE

  PRIVATE

  INCLUDE 'netcdf.inc'

  PUBLIC :: read_name_list_output_namelists
  PUBLIC :: init_name_list_output
  PUBLIC :: write_name_list_output
  PUBLIC :: close_name_list_output
  PUBLIC :: istime4name_list_output
  PUBLIC :: name_list_io_main_proc
  PUBLIC :: output_file
  PUBLIC :: parse_variable_groups

  !------------------------------------------------------------------------------------------------

  ! Parameters for naming all used Zaxis ID's in cdiZaxisID in TYPE t_output_file below

  INTEGER, PARAMETER, PUBLIC      :: ZA_surface             =  1
  ! Atmosphere
  INTEGER, PARAMETER, PUBLIC      :: ZA_hybrid              =  2
  INTEGER, PARAMETER, PUBLIC      :: ZA_hybrid_half         =  3
  INTEGER, PARAMETER, PUBLIC      :: ZA_depth_below_land    =  4
  INTEGER, PARAMETER, PUBLIC      :: ZA_depth_below_land_p1 =  5
  INTEGER, PARAMETER, PUBLIC      :: ZA_generic_snow        =  6
  INTEGER, PARAMETER, PUBLIC      :: ZA_generic_snow_p1     =  7
  INTEGER, PARAMETER, PUBLIC      :: ZA_pressure            =  8
  INTEGER, PARAMETER, PUBLIC      :: ZA_height              =  9
  INTEGER, PARAMETER, PUBLIC      :: ZA_altitude            = 10
  INTEGER, PARAMETER, PUBLIC      :: ZA_meansea             = 11
  ! Ocean
  INTEGER, PARAMETER, PUBLIC      :: ZA_depth               = 12
  INTEGER, PARAMETER, PUBLIC      :: ZA_depth_half          = 13
  INTEGER, PARAMETER, PUBLIC      :: ZA_generic_ice         = 14

  ! prefix for group identifier in output namelist
  CHARACTER(len=6) :: GRP_PREFIX = "group:"

  !------------------------------------------------------------------------------------------------

  TYPE(t_output_file), ALLOCATABLE, TARGET :: output_file(:)

  !------------------------------------------------------------------------------------------------
  ! TYPE t_reorder_info describes how local cells/edges/verts
  ! have to be reordered to get the global array.
  ! Below, "points" refers to either cells, edges or verts.
  !
  ! TODO[FP] Note that the "reorder_info" contains fields of *global*
  !          size (reorder_index). On the compute PEs these fields
  !          could be deallocated after the call to
  !          "transfer_reorder_info" in the setup phase!

  TYPE t_reorder_info
    INTEGER :: n_glb  ! Global number of points per physical patch
    INTEGER :: n_log  ! Global number of points in the associated logical patch
    INTEGER :: n_own  ! Number of own points (without halo, only belonging to phyiscal patch)
                      ! Only set on compute PEs, set to 0 on IO PEs
    INTEGER, ALLOCATABLE :: own_idx(:), own_blk(:)
                      ! idx and blk for own points, only set on compute PEs
    INTEGER, ALLOCATABLE :: pe_own(:)
                      ! n_own, gathered for all compute PEs (set on all PEs)
    INTEGER, ALLOCATABLE :: pe_off(:)
                      ! offset of contributions of PEs (set on all PEs)
    INTEGER, ALLOCATABLE :: reorder_index(:)
                      ! Index how to reorder the contributions of all compute PEs
                      ! into the global array (set on all PEs)

    ! only used when copying grid info from file (l_grid_info_from_file = .TRUE.):
    INTEGER, ALLOCATABLE :: log_dom_index(:)
                      ! Index where a point of the physical domains is in the logical domain
  END TYPE t_reorder_info


  TYPE t_grid_info
    REAL(wp), ALLOCATABLE :: lon   (:), lat   (:)
    REAL(wp), ALLOCATABLE :: lonv(:,:), latv(:,:)
  END TYPE t_grid_info

  ! TYPE t_patch_info contains the reordering info for cells, edges and verts
  TYPE t_patch_info
    TYPE(t_reorder_info) :: cells
    TYPE(t_reorder_info) :: edges
    TYPE(t_reorder_info) :: verts
    INTEGER :: log_patch_id

    ! pointer to communication pattern for GATHER operation;
    ! corresponds to physical or logical patch, depending on
    ! "l_output_phys_patch"
    TYPE(t_comm_pattern),  POINTER :: p_pat_c, p_pat_v, p_pat_e 
    
    ! global number of points, corresponds to physical or logical
    ! patch, depending on "l_output_phys_patch"
    INTEGER :: nblks_glb_c, nblks_glb_v, nblks_glb_e

    ! grid information: geographical locations of cells, edges, and
    ! vertices which is first collected on working PE 0 - from where
    ! it will be broadcasted to the pure I/O PEs.
    TYPE (t_grid_info) :: grid_c, grid_e, grid_v

    ! Filename of grid file, needed only if grid information is output
    ! to NetCDF since this information is normally not read and
    ! thus not present in the patch description
    CHARACTER(LEN=filename_max) :: grid_filename

    ! uuid of grid
    TYPE(t_uuid) :: grid_uuid

  END TYPE t_patch_info

  TYPE(t_patch_info),   ALLOCATABLE, TARGET :: patch_info (:)
  TYPE(t_reorder_info), ALLOCATABLE, TARGET :: lonlat_info(:,:)

  ! Number of output domains. This depends on l_output_phys_patch and is either the number
  ! of physical or the number of logical domains.

  INTEGER :: n_dom_out

  ! Total number of variables of all lists that have been tagged with
  ! "loutput=.TRUE."
  INTEGER :: n_allvars

  CHARACTER(LEN=vname_len) :: all_varlist(MAX_NVARS)

  !------------------------------------------------------------------------------------------------
  ! Currently, we use only 1 MPI window for all output files
  INTEGER mpi_win
  REAL(dp), POINTER :: mem_ptr_dp(:) ! Pointer to memory window (REAL*8)
  REAL(sp), POINTER :: mem_ptr_sp(:) ! Pointer to memory window (REAL*4)
  !------------------------------------------------------------------------------------------------
  ! Broadcast root for intercommunicator broadcasts form compute PEs to IO PEs using p_comm_work_2_io
  INTEGER :: bcast_root
  !------------------------------------------------------------------------------------------------

  ! Tags for communication between compute PEs and I/O PEs

  INTEGER, PARAMETER :: msg_io_start    = 12345
  INTEGER, PARAMETER :: msg_io_done     = 54321
  INTEGER, PARAMETER :: msg_io_shutdown = 99999

  ! TYPE t_datetime has no default constructor for setting all members to 0 or a defined value.
  ! Thus we declare a instance here which should have zeros everywhere since it is static.

  TYPE(t_datetime) :: zero_datetime

  ! Flag. If .TRUE. grid info will be copied from grid file, otherwise
  ! geographical locations of cells, edges, and vertices are first
  ! collected on working PE 0 - from where they will be broadcasted to
  ! the pure I/O PEs.
  LOGICAL :: l_grid_info_from_file

  !------------------------------------------------------------------------------------------------
  ! local copy of iadv_rcf, reducing dependencies and core specialities 
  INTEGER :: i_sample
  !------------------------------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output'

CONTAINS

  SUBROUTINE read_name_list_output_namelists( filename )

    ! The name is a bit strange, but if follows the convention to read namelists
    ! with a routine called read_XXX_namelist (plural is used here since several
    ! namelists are read).
    ! Please note the difference between name_list and namelist!

    CHARACTER(LEN=*), INTENT(IN)   :: filename

    INTEGER :: istat, i
    TYPE(t_output_name_list), POINTER :: p_onl
    CHARACTER(LEN=*), PARAMETER :: routine = 'read_name_list_output_namelists'
    INTEGER :: nnamelists
    LOGICAL :: lrewind

    ! Local variables corresponding to members of output_name_list
    INTEGER  :: filetype
    CHARACTER(LEN=8) :: namespace
    CHARACTER(LEN=filename_max) :: map_file
    INTEGER  :: mode
    INTEGER  :: taxis_tunit
    INTEGER  :: dom(max_phys_dom)
    INTEGER  :: output_time_unit
    REAL(wp) :: output_bounds(3,max_bounds)
    INTEGER  :: steps_per_file
    LOGICAL  :: include_last
    LOGICAL  :: output_grid
    CHARACTER(LEN=filename_max) :: output_filename
    LOGICAL  :: lwrite_ready
    CHARACTER(LEN=filename_max) :: ready_directory
    CHARACTER(LEN=vname_len)  :: ml_varlist(max_var_ml)
    CHARACTER(LEN=vname_len)  :: pl_varlist(max_var_pl)
    REAL(wp) :: p_levels(max_levels)
    CHARACTER(LEN=vname_len)  :: hl_varlist(max_var_hl)
    REAL(wp) :: h_levels(max_levels)
    INTEGER  :: remap
    LOGICAL  :: remap_internal
    REAL(wp) :: reg_lon_def(3)
    REAL(wp) :: reg_lat_def(3)
    INTEGER  :: gauss_tgrid_def
    REAL(wp) :: north_pole(2)
    TYPE(t_lon_lat_data), POINTER :: lonlat

    ! The namelist containing all variables above
    NAMELIST /output_nml/ &
      filetype,           &
      namespace,          &
      map_file,           &
      mode,               &
      taxis_tunit,        &
      dom,                &
      output_time_unit,   &
      output_bounds,      &
      steps_per_file,     &
      include_last,       &
      output_grid,        &
      output_filename,    &
      lwrite_ready,       &
      ready_directory,    &
      ml_varlist,         &
      pl_varlist,         &
      p_levels,           &
      hl_varlist,         &
      h_levels,           &
      remap,              &
      remap_internal,     &
      reg_lon_def,        &
      reg_lat_def,        &
      gauss_tgrid_def,    &
      north_pole

    ! Open input file and position to first namelist 'output_nml'

    CALL open_nml(TRIM(filename))

    ! As in COSMO, there may exist several output_nml namelists in the input file
    ! Loop until EOF is reached

    p_onl => NULL()
    first_output_name_list => NULL()
    name_list_output_active = .FALSE.
    nnamelists = 0
    lrewind = .TRUE.

    IF (.NOT. output_mode%l_nml) RETURN ! do not read output namelists if main switch is set to false
    DO
      CALL position_nml ('output_nml', lrewind=lrewind, status=istat)
      IF(istat /= POSITIONED) THEN
        CALL close_nml
        RETURN
      ENDIF
      lrewind = .FALSE.

      ! Set all variables in output_nml to their default values

      filetype           = FILETYPE_NC2 ! NetCDF
      namespace          = ' '
      map_file           = ' '
      mode               = 2
      taxis_tunit        = TUNIT_HOUR
      dom(:)             = -1
      output_time_unit   = 1
      output_bounds(:,:) = 0._wp
      steps_per_file     = 100
      include_last       = .TRUE.
      output_grid        = .FALSE.
      output_filename    = ' '
      lwrite_ready       = .FALSE.
      ready_directory    = ' '
      ml_varlist(:)      = ' '
      pl_varlist(:)      = ' '
      p_levels(:)        = 0._wp
      hl_varlist(:)      = ' '
      h_levels(:)        = 0._wp
      remap              = 0
      remap_internal     = .FALSE.
      reg_lon_def(:)     = 0._wp
      reg_lat_def(:)     = 0._wp
      gauss_tgrid_def    = 0
      north_pole(:)      = (/ 0._wp, 90._wp /)

      !------------------------------------------------------------------
      !  If this is a resumed integration, overwrite the defaults above 
      !  by values used in the previous integration.
      !  RJ: Disabled:
      !  The output_nml namelists need not be written to the restart files.
      !------------------------------------------------------------------
      !IF (is_restart_run()) THEN
      !  funit = open_and_restore_namelist('output_nml')
      !  READ(funit,NML=output_nml)
      !  CALL close_tmpfile(funit)
      !END IF

      ! Read output_nml

      READ (nnml, output_nml, iostat=istat)
      WRITE(message_text,'(a,i0)') 'Read namelist "output_nml", status = ', istat
      CALL message('',message_text)
      IF(istat > 0) THEN
        WRITE(message_text,'(a,i0)') 'Read error in namelist "output_nml", status = ', istat
        CALL finish(routine, message_text)
      ENDIF
      
      nnamelists = nnamelists+1
!LK      WRITE(message_text,'(a,i0)') 'Read namelist "output_nml", number = ', nnamelists
!LK      CALL message('',message_text)

      ! Check input

      ! We need dtime for this check
      IF(dtime<=0._wp) CALL finish(routine, 'dtime must be set before reading output namelists')

      ! Output bounds
      IF(output_bounds(1,1) < 0._wp .OR. &
         output_bounds(2,1) <= output_bounds(1,1) .OR. &
         output_bounds(3,1) <= dtime * grid_rescale_factor ) THEN
        CALL finish(routine,'Illegal output_bounds(:,1)')
      ENDIF

      DO i = 2, max_bounds-1
        IF(output_bounds(3,i) <= 0._wp) EXIT ! The last one
        IF(output_bounds(1,i) <= output_bounds(2,i-1)) &
          CALL finish(routine,'output_bounds not increasing')
        IF(output_bounds(2,i) <= output_bounds(1,i)) &
          CALL finish(routine,'output_bounds end <= start')
        IF(output_bounds(3,i) <  dtime * grid_rescale_factor ) &
          CALL finish(routine,'output_bounds inc < dtime')
      ENDDO

      ! For safety, at least last bounds triple must always be 0
      output_bounds(:,i:) = 0._wp

      ! Allocate next output_name_list

      IF(.NOT.ASSOCIATED(first_output_name_list)) THEN
        ! Allocate first name_list
        ALLOCATE(first_output_name_list)
        p_onl => first_output_name_list
        name_list_output_active = .TRUE.
      ELSE
        ! This is not the first one, p_onl points to the last one which was created
        ALLOCATE(p_onl%next)
        p_onl => p_onl%next
      ENDIF

      ! consistency check
      IF (output_grid .AND. (filetype == FILETYPE_GRB2)) THEN
        WRITE(message_text,*) &
          & 'WARNING: Output file type does not support grid output! => output_grid := .FALSE.'
        CALL message(TRIM(routine), TRIM(message_text))
        output_grid = .FALSE.
      END IF

      ! Set next output_name_list from values read

      p_onl%filetype         = filetype
      p_onl%namespace        = namespace
      p_onl%map_file         = map_file
      p_onl%mode             = mode
      p_onl%taxis_tunit      = taxis_tunit
      p_onl%dom(:)           = dom(:)
      p_onl%output_time_unit = output_time_unit

      ! output_bounds is always in seconds - the question is what to do with months or years
      ! output_time_unit: 1 = second, 2=minute, 3=hour, 4=day, 5=month, 6=year
      SELECT CASE(output_time_unit)
        CASE(1); p_onl%output_bounds(:,:) = output_bounds(:,:)
        CASE(2); p_onl%output_bounds(:,:) = output_bounds(:,:)*60._wp
        CASE(3); p_onl%output_bounds(:,:) = output_bounds(:,:)*3600._wp
        CASE(4); p_onl%output_bounds(:,:) = output_bounds(:,:)*86400._wp
        CASE(5); p_onl%output_bounds(:,:) = output_bounds(:,:)*86400._wp*30._wp  ! Not a real calender month
        CASE(6); p_onl%output_bounds(:,:) = output_bounds(:,:)*86400._wp*365._wp ! Not a real calender year
        CASE DEFAULT
          CALL finish(routine,'Illegal output_time_unit')
      END SELECT

      p_onl%steps_per_file   = steps_per_file
      p_onl%include_last     = include_last
      p_onl%output_grid      = output_grid
      p_onl%output_filename  = output_filename
      p_onl%lwrite_ready     = lwrite_ready
      p_onl%ready_directory  = ready_directory
      p_onl%ml_varlist(:)    = ml_varlist(:)
      p_onl%pl_varlist(:)    = pl_varlist(:)
      p_onl%p_levels         = p_levels
      p_onl%hl_varlist(:)    = hl_varlist(:)
      p_onl%h_levels         = h_levels
      p_onl%remap            = remap
      p_onl%remap_internal   = remap_internal
      p_onl%lonlat_id        = -1

      ! If "remap=1": lon-lat interpolation requested
      IF(remap/=0 .AND. remap/=1) &
        CALL finish(routine,'Unsupported value for remap')
      IF (remap == 1) THEN
        ! Register a lon-lat grid data structure in global list
        p_onl%lonlat_id = get_free_lonlat_grid()
        lonlat => lonlat_grid_list(p_onl%lonlat_id)

        lonlat%grid%reg_lon_def(:) = reg_lon_def(:)
        lonlat%grid%reg_lat_def(:) = reg_lat_def(:)
        lonlat%grid%north_pole(:)  = north_pole(:)

        IF(reg_lon_def(2)<=0._wp) CALL finish(routine,'Illegal LON increment')
        IF(reg_lat_def(2)==0._wp) CALL finish(routine,'Illegal LAT increment')
        IF(reg_lon_def(3)<=reg_lon_def(1)) CALL finish(routine,'end lon <= start lon')
        
        lonlat%grid%lon_dim = INT( (reg_lon_def(3)-reg_lon_def(1))/reg_lon_def(2) ) + 1
        lonlat%grid%lat_dim = INT( (reg_lat_def(3)-reg_lat_def(1))/reg_lat_def(2) ) + 1

        IF(lonlat%grid%lat_dim <= 0) CALL finish(routine,'Illegal LAT grid description')

        ! Flag those domains, which are used for this lon-lat grid:
        !     If dom(:) was not specified in namelist input, it is set
        !     completely to -1.  In this case all domains are wanted in
        !     the output
        IF (p_onl%dom(1) < 0)  THEN
          lonlat%l_dom(:) = .TRUE.
        ELSE
          DOM_LOOP : DO i = 1, max_dom
            IF (p_onl%dom(i) < 0) exit DOM_LOOP
            lonlat%l_dom( p_onl%dom(i) ) = .TRUE.
          ENDDO DOM_LOOP
        END IF
      ENDIF
      p_onl%gauss_tgrid_def      = gauss_tgrid_def

      p_onl%cur_bounds_triple= 1
      p_onl%next_output_time = p_onl%output_bounds(1,1)
      p_onl%n_output_steps   = 0
      p_onl%next => NULL()

      !-----------------------------------------------------
      ! Store the namelist for restart
      ! RJ: Disabled:
      ! The output_nml namelists need not be written to the restart files.
      !-----------------------------------------------------
      !IF(my_process_is_stdio())  THEN
      !  funit = open_tmpfile()
      !  WRITE(funit,NML=output_nml)
      !  CALL store_and_close_namelist(funit, 'output_nml')
      !ENDIF
      !-----------------------------------------------------
      ! write the contents of the namelist to an ASCII file
      !-----------------------------------------------------
      IF(my_process_is_stdio()) WRITE(nnml_output,nml=output_nml)

    ENDDO

    CALL close_nml

  END SUBROUTINE read_name_list_output_namelists


  !------------------------------------------------------------------------------------------------
 
  !> Looks for variable groups ("group:xyz") and replaces them
  !
  SUBROUTINE parse_variable_groups()
    CHARACTER(LEN=*), PARAMETER :: routine = &
      &  TRIM('mo_name_list_output/parse_variable_groups')
    !
    CHARACTER(LEN=VARNAME_LEN), ALLOCATABLE :: varlist(:), grp_vars(:), new_varlist(:)
    CHARACTER(LEN=VARNAME_LEN) :: vname, grp_name
    INTEGER :: nvars, ngrp_vars, i_typ, ierrstat, ivar, ntotal_vars, jvar
    CHARACTER(LEN=vname_len), POINTER :: in_varlist(:)
    TYPE (t_output_name_list), POINTER :: p_onl

    ntotal_vars = total_number_of_variables()
    ! temporary variables needed for variable group parsing
    ALLOCATE(varlist(ntotal_vars), grp_vars(ntotal_vars), &
      &      new_varlist(ntotal_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT
      
      ! process i_typ=ml_varlist, pl_varlist, hl_varlist:
      DO i_typ = 1, 3

        IF (i_typ == 1)  in_varlist => p_onl%ml_varlist
        IF (i_typ == 2)  in_varlist => p_onl%pl_varlist
        IF (i_typ == 3)  in_varlist => p_onl%hl_varlist

        ! Get the number of variables in varlist
        nvars = 1
        DO
          IF (nvars>SIZE(in_varlist))   EXIT
          IF (in_varlist(nvars) == ' ') EXIT
          nvars = nvars + 1
        END DO
        nvars = nvars - 1

        varlist(1:nvars) = in_varlist(1:nvars)
        varlist((nvars+1):ntotal_vars) = " "
        ! look for variable groups ("group:xyz") and replace them:
        DO ivar = 1, nvars
          vname = in_varlist(ivar)
          IF (INDEX(vname, GRP_PREFIX) > 0) THEN
            ! this is a group identifier
            grp_name = vname((LEN(TRIM(GRP_PREFIX))+1) : LEN(vname))
            CALL collect_group(grp_name, grp_vars, ngrp_vars)
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

        IF (i_typ == 1)  p_onl%ml_varlist(1:nvars) = varlist(1:nvars) 
        IF (i_typ == 2)  p_onl%pl_varlist(1:nvars) = varlist(1:nvars)
        IF (i_typ == 3)  p_onl%hl_varlist(1:nvars) = varlist(1:nvars)
      END DO ! i_typ = 1,3
      p_onl => p_onl%next
      
    END DO ! p_onl

    DEALLOCATE(varlist, grp_vars, new_varlist, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE parse_variable_groups


  !------------------------------------------------------------------------------------------------


  SUBROUTINE init_name_list_output(lprintlist, isample)

#ifndef NOMPI
    USE mpi, ONLY: MPI_ROOT, MPI_PROC_NULL
#endif

    LOGICAL, OPTIONAL, INTENT(in) :: lprintlist
    INTEGER, OPTIONAL, INTENT(in) :: isample

    ! For a list of all variables, enable the following:
    LOGICAL :: l_print_list

    INTEGER :: i, j, nfiles, i_typ, nvl, vl_list(max_var_lists), jp
    CHARACTER(LEN=2) :: lev_type
    INTEGER          :: ilev_type, idom, ierrstat, jg, idom_log
    TYPE (t_output_name_list), POINTER :: p_onl
    TYPE (t_output_file), POINTER :: p_of
    TYPE(t_list_element), POINTER :: element

    CHARACTER(LEN=*), PARAMETER :: routine = &
      &  TRIM('mo_name_list_output/init_name_list_output')
    REAL(wp), ALLOCATABLE :: lonv(:,:,:), latv(:,:,:)

    l_print_list = .FALSE.
    IF (PRESENT(lprintlist)) l_print_list = lprintlist

    i_sample = 1
    IF (PRESENT(isample)) i_sample = isample

    ! For hexagons, we still copy grid info from file; for triangular
    ! grids we have a faster method without file access:
    l_grid_info_from_file = (global_cell_type == 6)


    DO i = 1, nvar_lists

      ! print list of all variables
      IF (l_print_list) THEN
        IF (my_process_is_stdio()) THEN
          WRITE(message_text,'(3a, i2)') &
               'Var_list name: ',TRIM(var_lists(i)%p%name), &
               ' Patch: ',var_lists(i)%p%patch_id
          CALL message('',message_text)
          element => var_lists(i)%p%first_list_element
          DO
            IF(.NOT. ASSOCIATED(element)) EXIT
            WRITE (message_text,'(a,a,l1,a,a)') &
                 &     '    ',element%field%info%name,              &
                 &            element%field%info%loutput, '  ',     &
                 &            trim(element%field%info%cf%long_name)
            CALL message('',message_text)
            element => element%next_list_element
          ENDDO
        ENDIF
      ENDIF ! IF (l_print_list)

    ENDDO

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

    ! ---------------------------------------------------------------------------

    ! Set the number of output domains depending on
    ! l_output_phys_patch

    ! Replicate physical domain setup, only the number of domains and
    ! the logical ID is needed
    IF (use_async_name_list_io .AND.  &
      & .NOT. my_process_is_mpi_test()) THEN
      CALL p_bcast(n_phys_dom, bcast_root, p_comm_work_2_io)
      DO jg = 1, n_phys_dom
        CALL p_bcast(p_phys_patch(jg)%logical_id, bcast_root, p_comm_work_2_io)
      ENDDO
    END IF

    IF(l_output_phys_patch) THEN
      n_dom_out = n_phys_dom
    ELSE
      n_dom_out = n_dom
    ENDIF

    ! allocate patch info data structure
    ALLOCATE(patch_info(n_dom_out))
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
        END IF
      ENDIF
    ENDDO ! jp

    ! ---------------------------------------------------------------------------

    ! Prepare the output of grid information: For each
    ! physical/logical patch we must collect the geographical
    ! locations of cells, edges, and vertices

    ! Pure I/O PEs may skip this...
    IF (.NOT. (use_async_name_list_io .AND. my_process_is_io()) .AND.  &
      & .NOT. l_grid_info_from_file ) THEN

      ! Go over all output domains
      DO idom = 1, n_dom_out
        
        ! logical domain ID
        idom_log = patch_info(idom)%log_patch_id
        
        !-- collect domain data on working PE 0
        ! --cells
        ALLOCATE(lonv(nproma, p_patch(idom_log)%nblks_c, global_cell_type), &
          &      latv(nproma, p_patch(idom_log)%nblks_c, global_cell_type), &
          &      STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
        CALL cf_1_1_grid_cells(p_patch(idom_log), lonv, latv)
        CALL collect_grid_info(patch_info(idom)%nblks_glb_c,        &
          &                    p_patch(idom_log)%nblks_c,           &
          &                    p_patch(idom_log)%cells%center,      &
          &                    lonv, latv,                          &
          &                    patch_info(idom)%grid_c,             &
          &                    global_cell_type,                    &
          &                    patch_info(idom)%p_pat_c)
        DEALLOCATE(lonv, latv, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')        

        !-- edges
        ALLOCATE(lonv(nproma, p_patch(idom_log)%nblks_e, 4), &
          &      latv(nproma, p_patch(idom_log)%nblks_e, 4), &
          &      STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
        CALL cf_1_1_grid_edges(p_patch(idom_log), lonv, latv)
        CALL collect_grid_info(patch_info(idom)%nblks_glb_e,        &
          &                    p_patch(idom_log)%nblks_e,           &
          &                    p_patch(idom_log)%edges%center,      &
          &                    lonv, latv,                          &
          &                    patch_info(idom)%grid_e,             &
          &                    4,                                   &
          &                    patch_info(idom)%p_pat_e)
        DEALLOCATE(lonv, latv, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')        

        !-- verts
        ALLOCATE(lonv(nproma, p_patch(idom_log)%nblks_v, 9-global_cell_type), &
          &      latv(nproma, p_patch(idom_log)%nblks_v, 9-global_cell_type), &
          &      STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
        CALL cf_1_1_grid_verts(p_patch(idom_log), lonv, latv)
        CALL collect_grid_info(patch_info(idom)%nblks_glb_v,        &
          &                    p_patch(idom_log)%nblks_v,           &
          &                    p_patch(idom_log)%verts%vertex,      &
          &                    lonv, latv,                          &
          &                    patch_info(idom)%grid_v,             &
          &                    9-global_cell_type,                  &
          &                    patch_info(idom)%p_pat_v)
        DEALLOCATE(lonv, latv, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')        
        
      END DO
    END IF

    ! ---------------------------------------------------------------------------

    ! If async IO is used, replicate data (mainly the variable lists) on IO procs

#ifndef NOMPI
    IF(use_async_name_list_io) THEN
      CALL replicate_data_on_io_procs

      IF (.NOT. l_grid_info_from_file) THEN
        ! Clear patch_info fields clon, clat, etc. (especially on work
        ! PE 0) since they aren't needed there any longer.
        IF ( (.NOT. my_process_is_io()) .AND. &
          &  (.NOT. my_process_is_mpi_test())) THEN
          DO idom = 1, n_dom_out
            DEALLOCATE(patch_info(idom)%grid_c%lon,  patch_info(idom)%grid_c%lat,  &
              &        patch_info(idom)%grid_e%lon,  patch_info(idom)%grid_e%lat,  &
              &        patch_info(idom)%grid_v%lon,  patch_info(idom)%grid_v%lat,  &
              !
              &        patch_info(idom)%grid_c%lonv, patch_info(idom)%grid_c%latv, &
              &        patch_info(idom)%grid_e%lonv, patch_info(idom)%grid_e%latv, &
              &        patch_info(idom)%grid_v%lonv, patch_info(idom)%grid_v%latv, &
              &        STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
          END DO
        END IF
      END IF
    END IF
#endif

    ! Set the number of domains in output and the patch reorder information
    CALL set_patch_info
    
    ! Get the number of output files needed (by counting the domains per name list)

    p_onl => first_output_name_list
    nfiles = 0

    DO

      IF(.NOT.ASSOCIATED(p_onl)) THEN
        EXIT
      ENDIF

      ! If dom(:) was not specified in namelist input, it is set completely to -1.
      ! In this case all domains are wanted in the output, so set it here
      ! appropriately - this cannot be done during reading of the namelists
      ! since the number of physical domains is not known there.

      IF(p_onl%dom(1) <= 0) THEN
        DO i = 1, n_dom_out
          p_onl%dom(i) = i
        ENDDO
      ENDIF

      DO i = 1, SIZE(p_onl%dom)
        IF(p_onl%dom(i) <= 0) EXIT ! Last one was reached
        IF(p_onl%dom(i) > n_dom_out) THEN
          WRITE(message_text,'(a,i6,a)') &
            'Illegal domain number ',p_onl%dom(i),' in name list input'
          CALL finish(routine,message_text)
        ENDIF
        DO i_typ = 1, 3
          ! Check if name_list has variables of corresponding type
          IF(i_typ == 1 .AND. p_onl%ml_varlist(1) == ' ') CYCLE
          IF(i_typ == 2 .AND. p_onl%pl_varlist(1) == ' ') CYCLE
          IF(i_typ == 3 .AND. p_onl%hl_varlist(1) == ' ') CYCLE
          nfiles = nfiles+1
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
    output_file(:)%my_mem_win_off = 0_i8 ! Set if async IO is enabled

    ! Loop over all output namelists, set up the output_file struct for all associated files

    p_onl => first_output_name_list
    nfiles = 0

    DO

      IF(.NOT.ASSOCIATED(p_onl)) EXIT

      ! Loop over all domains for which this name list should be used

      DO i = 1, SIZE(p_onl%dom)
        IF(p_onl%dom(i) <= 0) EXIT ! Last one was reached
        idom = p_onl%dom(i)

        ! Loop over model/pressure/height levels

        DO i_typ = 1, 3

          ! Check if name_list has variables of corresponding type
          IF(i_typ == 1 .AND. p_onl%ml_varlist(1) == ' ') CYCLE
          IF(i_typ == 2 .AND. p_onl%pl_varlist(1) == ' ') CYCLE
          IF(i_typ == 3 .AND. p_onl%hl_varlist(1) == ' ') CYCLE

          nfiles = nfiles+1
          p_of => output_file(nfiles)

          SELECT CASE(i_typ)
            CASE(1); lev_type = 'ML'; ilev_type=level_type_ml
            CASE(2); lev_type = 'PL'; ilev_type=level_type_pl
            CASE(3); lev_type = 'HL'; ilev_type=level_type_hl
          END SELECT

          ! Set prefix of output_file name

          WRITE(p_of%filename_pref,'(a,"_DOM",i2.2,"_",a)') &
            TRIM(p_onl%output_filename),idom,lev_type

          p_of%phys_patch_id = idom
          p_of%log_patch_id  = patch_info(idom)%log_patch_id
          p_of%output_type   = p_onl%filetype
          p_of%name_list     => p_onl
          p_of%remap         = p_onl%remap

          p_of%start_time    = start_time(p_of%log_patch_id)
          p_of%end_time      = end_time(p_of%log_patch_id)
          p_of%initialized   = .FALSE.

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

          ! build a total varlist of variables tagged with
          ! "loutput=.TRUE."
          IF (p_of%remap == 1) THEN
            CALL get_all_var_names(all_varlist, n_allvars, opt_loutput=.TRUE.,    &
              &       opt_vlevel_type=ilev_type, opt_hor_intp_type=HINTP_TYPE_LONLAT,  &
              &       opt_patch_id=patch_info(idom)%log_patch_id)
          ELSE
            CALL get_all_var_names(all_varlist, n_allvars, opt_loutput=.TRUE.,    &
              &       opt_vlevel_type=ilev_type, opt_patch_id=patch_info(idom)%log_patch_id)
          END IF


          SELECT CASE(i_typ)
            CASE(1)
              IF (toupper(TRIM(p_onl%ml_varlist(1))) == "ALL") THEN
                IF (n_allvars > 0) &
                  CALL add_varlist_to_output_file(p_of,vl_list(1:nvl), &
                  &                               all_varlist(1:n_allvars))
              ELSE
                CALL add_varlist_to_output_file(p_of,vl_list(1:nvl),p_onl%ml_varlist)
              END IF
            CASE(2)
              IF (toupper(TRIM(p_onl%pl_varlist(1))) == "ALL") THEN
                IF (n_allvars > 0) &
                  CALL add_varlist_to_output_file(p_of,vl_list(1:nvl), &
                  &                               all_varlist(1:n_allvars))
              ELSE
                CALL add_varlist_to_output_file(p_of,vl_list(1:nvl),p_onl%pl_varlist)
              END IF
            CASE(3)
              IF (toupper(TRIM(p_onl%hl_varlist(1))) == "ALL") THEN
                IF (n_allvars > 0) &
                  CALL add_varlist_to_output_file(p_of,vl_list(1:nvl), &
                  &                               all_varlist(1:n_allvars))
              ELSE
                CALL add_varlist_to_output_file(p_of,vl_list(1:nvl),p_onl%hl_varlist)
              END IF
          END SELECT

        ENDDO

      ENDDO

      p_onl => p_onl%next

    ENDDO

    ! Set ID of process doing I/O
  
    DO i = 1, nfiles
      IF(use_async_name_list_io) THEN
        ! Asynchronous I/O
        IF(process_mpi_io_size==0) CALL finish(routine,'Async IO but no IO procs')
        ! Currently the output files are assigned round robin to the I/O PEs.
        ! Maybe this needs a finer adjustment in the future ...
        output_file(i)%io_proc_id = MOD(i-1,process_mpi_io_size) + p_io_pe0
      ELSE
        ! Normal I/O done by the standard I/O processor
        output_file(i)%io_proc_id = process_mpi_stdio_id
      ENDIF
    ENDDO

    ! If async IO is used, initialize the memory window for communication

#ifndef NOMPI
    IF(use_async_name_list_io) CALL init_memory_window
#endif

    CALL message(routine,'Done')

  END SUBROUTINE init_name_list_output

  !------------------------------------------------------------------------------------------------

  SUBROUTINE add_varlist_to_output_file(p_of, vl_list, varlist)

    TYPE (t_output_file), INTENT(INOUT) :: p_of
    INTEGER, INTENT(IN) :: vl_list(:)
    CHARACTER(LEN=*), INTENT(IN) :: varlist(:)

    INTEGER :: ivar, i, iv, idx, idx_t, idx_x, idx_y, tl, grid_of, grid_var
    LOGICAL :: found, found_1, found_2
    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_desc), TARGET  ::  var_desc_1, var_desc_2   !< variable descriptor
    TYPE(t_var_desc), POINTER ::  p_var_desc               !< variable descriptor (pointer)

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/add_varlist_to_output_file'

    ! Get the number of variables in varlist
    DO ivar = 1, SIZE(varlist)
      IF(varlist(ivar) == ' ') EXIT ! Last one reached
    ENDDO

    ! Allocate a list of variable descriptors: 
    p_of%max_vars = ivar-1
    p_of%num_vars = 0
    ALLOCATE(p_of%var_desc(p_of%max_vars))
    DO ivar = 1,(ivar-1)
      ! Nullify pointers in p_of%var_desc
      p_of%var_desc(ivar)%r_ptr => NULL()
      DO i = 1, max_time_levels
        p_of%var_desc(ivar)%tlev_ptr(i)%p => NULL()
      ENDDO
    END DO ! ivar

    ! Allocate array of variable descriptions
    DO ivar = 1,(ivar-1)

      found_1 = .FALSE.
      found_2 = .FALSE.
      ! Nullify pointers
      var_desc_1%r_ptr => NULL()
      var_desc_2%r_ptr => NULL()
      DO i = 1, max_time_levels
        var_desc_1%tlev_ptr(i)%p => NULL()
        var_desc_2%tlev_ptr(i)%p => NULL()
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

          IF (p_of%name_list%remap==1) THEN
            ! If lon-lat variable is requested, skip variable if it
            ! does not correspond to the same lon-lat grid:
            grid_of  = p_of%name_list%lonlat_id
            grid_var = element%field%info%hor_interp%lonlat_id

            IF (grid_of /= grid_var) CYCLE
          ELSE
            ! On the other hand: If no lon-lat interpolation is
            ! requested for this output file, skip all variables of
            ! this kind:
            IF (element%field%info%hgrid == GRID_REGULAR_LONLAT) CYCLE
          END IF ! (remap/=1)

          ! Do not inspect element if it is a container
          IF(element%field%info%lcontainer) CYCLE

          ! find suffix position for component and time level indices:
          idx_x = INDEX(element%field%info%name,'.X')
          idx_y = INDEX(element%field%info%name,'.Y')
          idx_t = INDEX(element%field%info%name,'.TL')

          idx = vname_len
          IF (idx_t > 0) idx=MIN(idx, idx_t)
          IF (idx_x > 0) idx=MIN(idx, idx_x)
          IF (idx_y > 0) idx=MIN(idx, idx_y)
          IF (idx==vname_len) idx=0
          
          ! Check for matching name
          IF(idx == 0) THEN
            IF(varlist(ivar) /= element%field%info%name) CYCLE
          ELSE
            IF(varlist(ivar) /= element%field%info%name(1:idx-1)) CYCLE
          ENDIF

          ! Found it, add it to the variable list of output file
          
          ! If we are dealing with an edge-based variable: look for
          ! both, the x and y component
          IF ((p_of%name_list%remap/=1) .OR. (idx_y == 0)) THEN
            p_var_desc => var_desc_1
            found      =  found_1
          ELSE
            p_var_desc => var_desc_2
            found      =  found_2
          END IF

          IF(idx_t == 0) THEN
            ! Not time level dependent
            IF(found) CALL finish(routine,'Duplicate var name: '//TRIM(varlist(ivar)))
            p_var_desc%r_ptr => element%field%r_ptr
            p_var_desc%info  = element%field%info
          ELSE
            IF(found) THEN
              ! We have already the info field, make some plausibility checks:
              IF(ANY(p_var_desc%info%used_dimensions(:) /=  &
                element%field%info%used_dimensions(:))) THEN
                CALL message(routine, "Var "//TRIM(element%field%info%name))
                CALL finish(routine,'Dimension mismatch TL variable: '//TRIM(varlist(ivar)))
              END IF
              ! There must not be a TL independent variable with the same name
              IF(ASSOCIATED(p_var_desc%r_ptr)) &
                CALL finish(routine,'Duplicate var name: '//TRIM(varlist(ivar)))
              ! Maybe some more members of info should be tested ...
            ELSE
              ! Variable encountered the first time, set info field ...
              p_var_desc%info = element%field%info
              ! ... and set name without .TL# suffix
              p_var_desc%info%name = element%field%info%name(1:idx_t-1)
            ENDIF
            ! Get time level
            tl = ICHAR(element%field%info%name(idx_t+3:idx_t+3)) - ICHAR('0')
            IF(tl<=0 .OR. tl>max_time_levels) &
              CALL finish(routine, 'Illegal time level in '//TRIM(element%field%info%name))
            IF(ASSOCIATED(p_var_desc%tlev_ptr(tl)%p)) &
              CALL finish(routine, 'Duplicate time level for '//TRIM(element%field%info%name))
            p_var_desc%tlev_ptr(tl)%p => element%field%r_ptr
          ENDIF

          IF ((p_of%name_list%remap/=1) .OR. (idx_y == 0)) THEN
            found_1 = .TRUE.
          ELSE
            found_2 = .TRUE.
          END IF

        ENDDO

      ENDDO ! i = 1, SIZE(vl_list)

      ! Check that at least one element with this name has been found

      IF(.NOT. found_1) &
        CALL finish(routine,'Output name list variable not found: '//TRIM(varlist(ivar)))
      
      ! append variable descriptor to list; append two different
      ! variables (X and Y) if we have a lon-lat interpolated variable
      ! defined on edges:
      IF (found_2) THEN
        var_desc_1%info%name = TRIM(var_desc_1%info%name)//".X"
        CALL add_var_desc(p_of, var_desc_1)
        var_desc_2%info%name = TRIM(var_desc_2%info%name)//".Y"
        CALL add_var_desc(p_of, var_desc_2)
      ELSE
        CALL add_var_desc(p_of, var_desc_1)
      END IF

    ENDDO ! ivar = 1,(ivar-1)

  END SUBROUTINE add_varlist_to_output_file


  !------------------------------------------------------------------------------------------------

  SUBROUTINE set_patch_info

    INTEGER :: jp, jl, ierrstat, jg

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/set_patch_info'

    DO jp = 1, n_dom_out

      jl = patch_info(jp)%log_patch_id

      IF(.NOT.my_process_is_io()) THEN
        ! Set reorder_info on work and test PE
        CALL set_reorder_info(jp, p_patch(jl)%n_patch_cells_g, p_patch(jl)%n_patch_cells, &
                              p_patch(jl)%cells%owner_mask, p_patch(jl)%cells%phys_id,    &
                              p_patch(jl)%cells%glb_index, patch_info(jp)%cells)

        CALL set_reorder_info(jp, p_patch(jl)%n_patch_edges_g, p_patch(jl)%n_patch_edges, &
                              p_patch(jl)%edges%owner_mask, p_patch(jl)%edges%phys_id,    &
                              p_patch(jl)%edges%glb_index, patch_info(jp)%edges)

        CALL set_reorder_info(jp, p_patch(jl)%n_patch_verts_g, p_patch(jl)%n_patch_verts, &
                              p_patch(jl)%verts%owner_mask, p_patch(jl)%verts%phys_id,    &
                              p_patch(jl)%verts%glb_index, patch_info(jp)%verts)
        ! Set grid_filename on work and test PE
        patch_info(jp)%grid_filename = TRIM(p_patch(jl)%grid_filename)
        ! Set UUID on work and test PE
        patch_info(jp)%grid_uuid = p_patch(jl)%grid_uuid
      ENDIF
#ifndef NOMPI
      IF(use_async_name_list_io .AND. .NOT. my_process_is_mpi_test()) THEN
        ! Transfer reorder_info to IO PEs
        CALL transfer_reorder_info(patch_info(jp)%cells)
        CALL transfer_reorder_info(patch_info(jp)%edges)
        CALL transfer_reorder_info(patch_info(jp)%verts)
        CALL p_bcast(patch_info(jp)%grid_filename, bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(jp)%grid_uuid%data, SIZE(patch_info(jp)%grid_uuid%data),  &
          &          bcast_root, p_comm_work_2_io)
      ENDIF
#endif

    ENDDO ! jp

    ! A similar process as above - for the lon-lat grids
    ALLOCATE(lonlat_info(n_lonlat_grids, n_dom), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
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
          CALL transfer_reorder_info(lonlat_info(jl,jg))
        ENDIF
#endif
      END DO ! jg
    ENDDO ! jl

  END SUBROUTINE set_patch_info


  !------------------------------------------------------------------------------------------------
  !> SUBROUTINE collect_grid_info
  !
  !  Prepare the output of grid information: For each physical/logical
  !  patch we must collect the geographical locations of cells, edges,
  !  and vertices first on working PE 0 - from where it will be
  !  broadcasted to the pure I/O PEs.
  !
  SUBROUTINE collect_grid_info(nblks_glb, nblks_loc, in_lonlat, lonv, latv, out_lonlat, &
    &                          dim3, p_pat)
    INTEGER,                          INTENT(IN)    :: nblks_glb,     &  ! global number of blocks
      &                                                nblks_loc         ! local  number of blocks
    TYPE(t_geographical_coordinates), INTENT(IN)    :: in_lonlat(:,:)
    REAL(wp),                         INTENT(IN)    :: lonv(:,:,:), latv(:,:,:)
    TYPE(t_grid_info),                INTENT(INOUT) :: out_lonlat
    INTEGER,                          INTENT(IN)    :: dim3
    TYPE(t_comm_pattern),  POINTER                  :: p_pat
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine  =  TRIM('mo_name_list_output/collect_grid_info')
    INTEGER               :: ierrstat, jb, jc, idim
    REAL(wp), ALLOCATABLE :: r_tmp_lon(:,:), r_tmp_lat(:,:)

    ! skip this on test PE...
    IF (my_process_is_mpi_test()) RETURN

    ! allocate destination (on work root)
    IF ( my_process_is_mpi_workroot() ) THEN
      ALLOCATE(out_lonlat%lon (nproma*nblks_glb),      out_lonlat%lat (nproma*nblks_glb),      &
        &      out_lonlat%lonv(dim3,nproma*nblks_glb), out_lonlat%latv(dim3,nproma*nblks_glb), &
        &      STAT=ierrstat)
    ELSE
      ALLOCATE(out_lonlat%lon(1), out_lonlat%lat (1),            &
        &      out_lonlat%lonv(dim3,1), out_lonlat%latv(dim3,1), &
        &      STAT=ierrstat)
    END IF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! allocate temporary fields:
    IF (my_process_is_mpi_test() .OR. my_process_is_mpi_workroot()) THEN
      ALLOCATE(r_tmp_lon (nproma, nblks_glb), r_tmp_lat (nproma, nblks_glb), STAT=ierrstat)
    ELSE
      ALLOCATE(r_tmp_lon (nproma, nblks_loc), r_tmp_lat (nproma, nblks_loc), STAT=ierrstat)
    ENDIF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    
    !-- part 1: exchange lon/lat coordinates:
    IF (.NOT. my_process_is_mpi_test()) THEN
      ! copy coordinates into sender array:
      DO jb=1,nblks_loc
        DO jc=1,nproma
          r_tmp_lon(jc,jb) = in_lonlat(jc,jb)%lon
          r_tmp_lat(jc,jb) = in_lonlat(jc,jb)%lat
        END DO
      END DO
    ENDIF
    ! gather data on work root
    IF(.NOT. my_process_is_mpi_seq()) THEN
      CALL exchange_data(p_pat, RECV=r_tmp_lon,  SEND=r_tmp_lon )
      CALL exchange_data(p_pat, RECV=r_tmp_lat,  SEND=r_tmp_lat )
    END IF
    ! on work root: reshape into 1D arrays:
    IF ( my_process_is_mpi_workroot() ) THEN
      out_lonlat%lon(:)  = RESHAPE(r_tmp_lon(:,:), (/ nproma*nblks_glb /))
      out_lonlat%lat(:)  = RESHAPE(r_tmp_lat(:,:), (/ nproma*nblks_glb /))
    END IF

    !-- part 2: exchange vertex lon/lat coordinates:
    DO idim=1,dim3
      IF (.NOT. my_process_is_mpi_test()) THEN
        ! copy coordinates into sender array:
        DO jb=1,nblks_loc
          DO jc=1,nproma
            r_tmp_lon(jc,jb) = lonv(jc,jb,idim)
            r_tmp_lat(jc,jb) = latv(jc,jb,idim)
          END DO
        END DO
      ENDIF
      ! gather data on work root
      IF(.NOT. my_process_is_mpi_seq()) THEN
        CALL exchange_data(p_pat, RECV=r_tmp_lon,  SEND=r_tmp_lon )
        CALL exchange_data(p_pat, RECV=r_tmp_lat,  SEND=r_tmp_lat )
      END IF
      ! on work root: reshape into 1D arrays:
      IF ( my_process_is_mpi_workroot() ) THEN
        out_lonlat%lonv(idim,:) = RESHAPE(r_tmp_lon(:,:), (/ nproma*nblks_glb /))
        out_lonlat%latv(idim,:) = RESHAPE(r_tmp_lat(:,:), (/ nproma*nblks_glb /))
      END IF
    END DO ! idim

    ! clean up
    DEALLOCATE(r_tmp_lon, r_tmp_lat, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE collect_grid_info


  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "clonv", "clatv" into a CF1.1 compatible form.
  !
  !  based on SUBROUTINE mo_gridrefinement::write_patch(p)
  !
  SUBROUTINE cf_1_1_grid_cells(p_patch, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER :: jc, jb, j, iidx, iblk,                  &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom

    rl_start   = 1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    lonv(:,:,:) = 0.0_wp
    latv(:,:,:) = 0.0_wp
    DO j = 1, global_cell_type
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          iidx = p_patch%cells%vertex_idx(jc,jb,j)
          iblk = p_patch%cells%vertex_blk(jc,jb,j)
          lonv(jc,jb,j) = p_patch%verts%vertex(iidx,iblk)%lon
          latv(jc,jb,j) = p_patch%verts%vertex(iidx,iblk)%lat
        END DO
      END DO
    END DO
    WHERE (ABS(lonv(:,:,:)) < EPSILON(0.0_wp))
      lonv(:,:,:) = 0.0_wp
    END WHERE
    WHERE (ABS(latv(:,:,:)) < EPSILON(0.0_wp))
      latv(:,:,:) = 0.0_wp
    END WHERE
    DO j = 1, global_cell_type
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          IF (ABS(latv(jc,jb,j)) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
            lonv(jc,jb,j) = p_patch%cells%center(jc,jb)%lon
          ENDIF
        END DO
      END DO
    END DO
    
  END SUBROUTINE cf_1_1_grid_cells


  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "elonv", "elatv" into a CF1.1 compatible form.
  !
  !  based on SUBROUTINE mo_gridrefinement::write_patch(p)
  !
  SUBROUTINE cf_1_1_grid_edges(p_patch, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER  :: jc, jb, j, iidx, iblk,                 &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom
    REAL(wp) :: swap(4)

    rl_start   = 1
    rl_end     = min_rledge_int
    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

    lonv(:,:,:) = 0.0_wp
    latv(:,:,:) = 0.0_wp
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        iidx = p_patch%edges%vertex_idx(jc,jb,1)
        iblk = p_patch%edges%vertex_blk(jc,jb,1)
        lonv(jc,jb,1) = p_patch%verts%vertex(iidx,iblk)%lon
        latv(jc,jb,1) = p_patch%verts%vertex(iidx,iblk)%lat

        iidx = p_patch%edges%vertex_idx(jc,jb,2)
        iblk = p_patch%edges%vertex_blk(jc,jb,2)
        lonv(jc,jb,3) = p_patch%verts%vertex(iidx,iblk)%lon
        latv(jc,jb,3) = p_patch%verts%vertex(iidx,iblk)%lat

        IF (p_patch%edges%cell_idx(jc,jb,1) > 0) THEN
          iidx = p_patch%edges%cell_idx(jc,jb,1)
          iblk = p_patch%edges%cell_blk(jc,jb,1)
          lonv(jc,jb,4) = p_patch%cells%center(iidx,iblk)%lon
          latv(jc,jb,4) = p_patch%cells%center(iidx,iblk)%lat
        ELSE
          lonv(jc,jb,4) = 0._wp
          latv(jc,jb,4) = 0._wp
        ENDIF

        IF (p_patch%edges%cell_idx(jc,jb,2) > 0) THEN
          iidx = p_patch%edges%cell_idx(jc,jb,2)
          iblk = p_patch%edges%cell_blk(jc,jb,2)
          lonv(jc,jb,2) = p_patch%cells%center(iidx,iblk)%lon
          latv(jc,jb,2) = p_patch%cells%center(iidx,iblk)%lat
        ELSE
          lonv(jc,jb,2) = 0._wp
          latv(jc,jb,2) = 0._wp
        END IF
      END DO
    END DO

    WHERE (ABS(lonv(:,:,:)) < EPSILON(0.0_wp))
      lonv(:,:,:) = 0.0_wp
    END WHERE
    WHERE (ABS(latv(:,:,:)) < EPSILON(0.0_wp))
      latv(:,:,:) = 0.0_wp
    END WHERE

    DO j = 1, 4
      DO jb = i_startblk, i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          IF ( ABS(latv(jc,jb,j)) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
            lonv(jc,jb,j) = p_patch%edges%center(jc,jb)%lon
          ENDIF
        ENDDO
      ENDDO
    END DO
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        IF ( check_orientation(p_patch%edges%center(jc,jb)%lon, &
          &                    lonv(jc,jb,:), latv(jc,jb,:), 4) < 0 ) THEN
          swap(1:4) = lonv(jc,jb,4:1:-1)
          lonv(jc,jb,:) = swap(:)
          swap(1:4) = latv(jc,jb,4:1:-1)
          latv(jc,jb,:) = swap(:)
        END IF
      END DO
    ENDDO
    
  END SUBROUTINE cf_1_1_grid_edges


  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "vlonv", "vlatv" into a CF1.1 compatible form.
  !
  !  based on SUBROUTINE mo_gridrefinement::write_patch(p)
  !
  SUBROUTINE cf_1_1_grid_verts(p_patch, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER :: jc, jb, j, iidx, iblk,                  &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom

    rl_start   = 2
    rl_end     = min_rlvert
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%verts%start_blk(rl_start,1)
    i_endblk   = p_patch%verts%end_blk(rl_end,i_nchdom)    

    lonv(:,:,:) = 0.0_wp
    latv(:,:,:) = 0.0_wp
    DO j = 1,(9-global_cell_type)
      DO jb = i_startblk, i_endblk
        CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          IF ((p_patch%verts%cell_idx(jc,jb,j) == 0) .AND. &
            & (p_patch%verts%refin_ctrl(jc,jb) /= 1)) THEN
            iidx = p_patch%verts%cell_idx(jc,jb,5)
            iblk = p_patch%verts%cell_blk(jc,jb,5)
            lonv(jc,jb,7-j) = p_patch%cells%center(iidx,iblk)%lon
            latv(jc,jb,7-j) = p_patch%cells%center(iidx,iblk)%lat
          ELSE IF ((p_patch%verts%cell_idx(jc,jb,j) < 0) .OR. &
            &      (p_patch%verts%refin_ctrl(jc,jb) == 1)) THEN
            lonv(jc,jb,7-j) = 0._wp
            latv(jc,jb,7-j) = 0._wp
          ELSE
            iidx = p_patch%verts%cell_idx(jc,jb,j)
            iblk = p_patch%verts%cell_blk(jc,jb,j)
            lonv(jc,jb,7-j) = p_patch%cells%center(iidx,iblk)%lon
            latv(jc,jb,7-j) = p_patch%cells%center(iidx,iblk)%lat
          ENDIF
        ENDDO
      ENDDO
    END DO
  END SUBROUTINE cf_1_1_grid_verts


  !------------------------------------------------------------------------------------------------
  !>
  !! Sets the reorder_info for cells/edges/verts
  !! ATTENTION: This routine must only be called on work and test PE (i.e. not on IO PEs)
  !!            The arguments don't make sense on the IO PEs anyways

  SUBROUTINE set_reorder_info(phys_patch_id, n_points_g, n_points, owner_mask, phys_id, &
                              glb_index, p_ri)

    INTEGER, INTENT(IN) :: phys_patch_id   ! Physical patch ID
    INTEGER, INTENT(IN) :: n_points_g      ! Global number of cells/edges/verts in logical patch
    INTEGER, INTENT(IN) :: n_points        ! Local number of cells/edges/verts in logical patch
    LOGICAL, INTENT(IN) :: owner_mask(:,:) ! owner_mask for logical patch
    INTEGER, INTENT(IN) :: phys_id(:,:)    ! phys_id for logical patch
    INTEGER, INTENT(IN) :: glb_index(:)    ! glb_index for logical patch

    TYPE(t_reorder_info), INTENT(INOUT) :: p_ri ! Result: reorder info

    INTEGER :: i, n, il, ib, mpierr
    LOGICAL, ALLOCATABLE :: phys_owner_mask(:) ! owner mask for physical patch
    INTEGER, ALLOCATABLE :: glbidx_own(:), glbidx_glb(:), reorder_index_log_dom(:)

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/set_reorder_info'

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

    ! Get reorder_index

    ALLOCATE(p_ri%reorder_index(p_ri%n_glb))
    IF (l_grid_info_from_file) THEN
      ALLOCATE(p_ri%log_dom_index(p_ri%n_glb))
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
        p_ri%reorder_index(n) = reorder_index_log_dom(i)
      ENDIF
    ENDDO

    IF (l_grid_info_from_file) THEN
      n = 0
      DO i = 1, n_points_g
        IF(reorder_index_log_dom(i)>0) THEN
          n = n+1
          p_ri%log_dom_index(n) = i
        ENDIF
      ENDDO
    END IF

    ! Safety check
    IF(n/=p_ri%n_glb) CALL finish(routine,'Reordering failed')

    DEALLOCATE(phys_owner_mask)
    DEALLOCATE(glbidx_own)
    DEALLOCATE(glbidx_glb)
    DEALLOCATE(reorder_index_log_dom)

  END SUBROUTINE set_reorder_info


  !------------------------------------------------------------------------------------------------
  !> Sets the reorder_info for lon-lat-grids
  !
  SUBROUTINE set_reorder_info_lonlat(grid, intp, p_ri)
    TYPE(t_lon_lat_grid), INTENT(IN)    :: grid
    TYPE(t_lon_lat_intp), INTENT(IN)    :: intp
    TYPE(t_reorder_info), INTENT(INOUT) :: p_ri ! Result: reorder info

    CHARACTER(LEN=*), PARAMETER :: routine = &
      &   'mo_name_list_output/set_reorder_info_lonlat'
    INTEGER :: ierrstat, i, jc, jb, this_pe, mpierr, & 
    &          ioffset

    ! Just for safety
    IF(my_process_is_io()) CALL finish(routine, 'Must not be called on IO PEs')
    this_pe = get_my_mpi_work_id()

    p_ri%n_glb = grid%lon_dim * grid%lat_dim ! Total points in lon-lat grid
    p_ri%n_own = intp%nthis_local_pts        ! No. of own points

    ! Set index arrays to own cells/edges/verts
    ALLOCATE(p_ri%own_idx(p_ri%n_own), p_ri%own_blk(p_ri%n_own), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    jc = 0
    jb = 1
    DO i=1,p_ri%n_own
      jc = jc+1
      IF (jc>nproma) THEN
        jb = jb+1;  jc = 1
      END IF

      p_ri%own_idx(i) = jc
      p_ri%own_blk(i) = jb
    END DO ! i

    ! Gather the number of own points for every PE into p_ri%pe_own
    ALLOCATE(p_ri%pe_own(0:p_n_work-1), p_ri%pe_off(0:p_n_work-1), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
#ifndef NOMPI
    CALL MPI_Allgather(p_ri%n_own,  1, p_int, &
                       p_ri%pe_own, 1, p_int, &
                       p_comm_work, mpierr)
#else
    p_ri%pe_own(0) = p_ri%n_own
#endif
    ! Get offset within result array
    p_ri%pe_off(0) = 0
    DO i = 1, p_n_work-1
      p_ri%pe_off(i) = p_ri%pe_off(i-1) + p_ri%pe_own(i-1)
    ENDDO

    ! Get the global index numbers of the data when it is gathered on PE 0
    ! exactly in the same order as it is retrieved later during I/O
    ALLOCATE(p_ri%reorder_index(p_ri%n_glb), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    IF (l_grid_info_from_file) THEN
      ALLOCATE(p_ri%log_dom_index(p_ri%n_glb), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    END IF

    ioffset = p_ri%pe_off(this_pe)
    p_ri%reorder_index = -1
    DO i=1,intp%nthis_local_pts
      p_ri%reorder_index(intp%global_idx(i)) = ioffset + i
    END DO
    ! merge all fields across working PEs:
    p_ri%reorder_index = p_max(p_ri%reorder_index, &
      &                        comm=get_my_mpi_work_communicator())

    IF (l_grid_info_from_file) THEN
      ! mapping between logical and physical patch is trivial for
      ! lon-lat grids:
      p_ri%log_dom_index(:) = (/ (i, i=1,p_ri%n_glb) /)
    END IF

  END SUBROUTINE set_reorder_info_lonlat


  !------------------------------------------------------------------------------------------------
  !> setup_output_vlist:
  !! Sets up the vlist for a t_output_file structure

  SUBROUTINE setup_output_vlist(of)

    TYPE(t_output_file), INTENT(INOUT) :: of

    INTEGER :: k, nlev, nlevp1, nplev, nzlev, nzlevp1, znlev_soil, i_dom
    INTEGER :: ll_dim(2)
    INTEGER :: gridtype
    REAL(wp), ALLOCATABLE :: levels(:), levels_i(:), levels_m(:), p_lonlat(:)

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/setup_output_vlist'
    TYPE(t_lon_lat_data), POINTER :: lonlat
    LOGICAL :: lwrite_pzlev
    INTEGER :: idate, itime
    TYPE(t_datetime) :: ini_datetime

    ! Read map_file - we do it here and not during initialization
    ! since it is enough if only the output PE does the read (and not all)

    IF(of%name_list%map_file /= ' ') CALL read_map_file(of, of%name_list%map_file)


    IF (of%output_type == FILETYPE_GRB2) THEN
      ! since the current CDI-version does not fully support "GRID_UNSTRUCTURED", the 
      ! grid type is changed to "GRID_REFERENCE".
      gridtype = GRID_REFERENCE
    ELSE
      gridtype = GRID_UNSTRUCTURED
    ENDIF

    i_dom = of%phys_patch_id

    ! set initialization flag to true
    of%initialized = .TRUE.

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
    ! get instID
    of%cdiInstID = institutInq(0, 0, TRIM(of%name_list%namespace), "")
    CALL vlistDefInstitut(of%cdiVlistID,of%cdiInstID)


    ! 3. add horizontal grid descriptions

    IF(of%name_list%remap == 1) THEN

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

      ! TODO [FP]

      ! When specifying the north pole of the rotated lon-lat grid,
      ! CDO conversion to GRIB format yields an incorrect longitude
      ! axis, even if the north pole is identical to the lon-lat grid
      ! without rotation.

      ! As a TEMPORARY solution, the position of the rotated poly is
      ! not stored in the NetCDF file, s.t. GRIB conversion works
      ! correctly.

      ! CALL gridDefXpole(of%cdiLonLatGridID, grid%poleN(1)/pi_180)
      ! CALL gridDefYpole(of%cdiLonLatGridID, grid%poleN(2)/pi_180)
      ! Note: CALL gridDefAngle() not yet supported

      ALLOCATE(p_lonlat(ll_dim(1)))
      DO k = 1, ll_dim(1)
        p_lonlat(k) = lonlat%grid%reg_lon_def(1)  &
          &            +  REAL(k-1,wp)*lonlat%grid%reg_lon_def(2)
      ENDDO
      CALL gridDefXvals(of%cdiLonLatGridID, p_lonlat)
      DEALLOCATE(p_lonlat)

      ALLOCATE(p_lonlat(ll_dim(2)))
      DO k = 1, ll_dim(2)
        p_lonlat(k) = lonlat%grid%reg_lat_def(1)  &
          &            +  REAL(k-1,wp)*lonlat%grid%reg_lat_def(2)
      ENDDO
      CALL gridDefYvals(of%cdiLonLatGridID, p_lonlat)
      DEALLOCATE(p_lonlat)

    ELSE

      ! Cells

      of%cdiCellGridID = gridCreate(gridtype, patch_info(i_dom)%cells%n_glb)
      CALL gridDefNvertex(of%cdiCellGridID, global_cell_type)
      !
      CALL gridDefXname(of%cdiCellGridID, 'clon')
      CALL gridDefXlongname(of%cdiCellGridID, 'center longitude')
      CALL gridDefXunits(of%cdiCellGridID, 'radian')
      !
      CALL gridDefYname(of%cdiCellGridID, 'clat')
      CALL gridDefYlongname(of%cdiCellGridID, 'center latitude')
      CALL gridDefYunits(of%cdiCellGridID, 'radian')
      !
      CALL gridDefUUID(of%cdiCellGridID, patch_info(i_dom)%grid_uuid%data)
   
      ! Verts

      of%cdiVertGridID = gridCreate(gridtype, patch_info(i_dom)%verts%n_glb)
      CALL gridDefNvertex(of%cdiVertGridID, 9-global_cell_type)
      !
      CALL gridDefXname(of%cdiVertGridID, 'vlon')
      CALL gridDefXlongname(of%cdiVertGridID, 'vertex longitude')
      CALL gridDefXunits(of%cdiVertGridID, 'radian')
      !
      CALL gridDefYname(of%cdiVertGridID, 'vlat')
      CALL gridDefYlongname(of%cdiVertGridID, 'vertex latitude')
      CALL gridDefYunits(of%cdiVertGridID, 'radian')
      !
      CALL gridDefUUID(of%cdiVertGridID, patch_info(i_dom)%grid_uuid%data)

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
      CALL gridDefUUID(of%cdiEdgeGridID, patch_info(i_dom)%grid_uuid%data)

      of%cdiLonLatGridID = CDI_UNDEFID

      ! If wanted, set grid info
      IF(of%name_list%output_grid) THEN
        IF (l_grid_info_from_file) THEN
          CALL copy_grid_info(of)
        ELSE IF (.NOT. my_process_is_mpi_test()) THEN
          CALL set_grid_info(of)
        END IF
      END IF

    ENDIF

    !
    ! 4. add vertical grid descriptions
    !    RJ: This is copied from mo_io_vlist

    of%cdiZaxisID(:) = CDI_UNDEFID ! not all are set

    ! surface level
    of%cdiZaxisID(ZA_surface) = zaxisCreate(ZAXIS_SURFACE, 1)
    ALLOCATE(levels(1))
    levels(1) = 0.0_wp
    CALL zaxisDefLevels(of%cdiZaxisID(ZA_surface), levels)
    DEALLOCATE(levels)

    ! atm (pressure) height, ocean depth
    IF (iequations/=ihs_ocean) THEN ! atm 

      nlev   = num_lev(of%log_patch_id)
      nlevp1 = num_levp1(of%log_patch_id)
      ! introduce temporary variable znlev_soil, since global variable nlev_soil 
      ! is unknown to the I/O-Processor. Otherwise receive_patch_configuration in 
      ! mo_io_async complains about mismatch of levels. 
      znlev_soil = SIZE(zml_soil)-1

      ! Hybrid

      of%cdiZaxisID(ZA_hybrid)      = zaxisCreate(ZAXIS_HYBRID, nlev)
      of%cdiZaxisID(ZA_hybrid_half) = zaxisCreate(ZAXIS_HYBRID_HALF, nlevp1)

      ALLOCATE(levels(nlev))
      DO k = 1, nlev
        levels(k) = REAL(k,wp)
      END DO
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_hybrid), levels)
      DEALLOCATE(levels)
      CALL zaxisDefVct(of%cdiZaxisID(ZA_hybrid), 2*nlevp1, vct(1:2*nlevp1))
      !
      ALLOCATE(levels(nlevp1))
      DO k = 1, nlevp1
        levels(k) = REAL(k,wp)
      END DO
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_hybrid_half), levels)
      DEALLOCATE(levels)
      CALL zaxisDefVct(of%cdiZaxisID(ZA_hybrid_half), 2*nlevp1, vct(1:2*nlevp1))

      ! Define axis for output on mean sea level
      !
      of%cdiZaxisID(ZA_meansea) = zaxisCreate(ZAXIS_MEANSEA, 1)
      ALLOCATE(levels(1))
      levels(1) = 0.0_wp
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_meansea), levels)
      DEALLOCATE(levels)

      ! Define axes for soil model

      of%cdiZaxisID(ZA_depth_below_land_p1) = &
        & zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+2)
      ALLOCATE(levels(znlev_soil+2))
      levels(1) = 0._wp
      DO k = 1, znlev_soil+1
        levels(k+1) = zml_soil(k)*1000._wp  ! in mm
      END DO
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_depth_below_land_p1), levels)
      CALL zaxisDefUnits(of%cdiZaxisID(ZA_depth_below_land_p1), "mm")
      DEALLOCATE(levels)

      of%cdiZaxisID(ZA_depth_below_land) = &
        & zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+1)
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_depth_below_land), zml_soil*1000._wp) ! in mm
      CALL zaxisDefUnits (of%cdiZaxisID(ZA_depth_below_land), "mm")
      !
      of%cdiZaxisID(ZA_generic_snow_p1) = zaxisCreate(ZAXIS_GENERIC, nlev_snow+1)
      ALLOCATE(levels(nlev_snow+1))
      DO k = 1, nlev_snow+1
        levels(k) = REAL(k,wp)
      END DO
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_generic_snow_p1), levels)
      DEALLOCATE(levels)
      !
      of%cdiZaxisID(ZA_generic_snow) = zaxisCreate(ZAXIS_GENERIC, nlev_snow)
      ALLOCATE(levels(nlev_snow))
      DO k = 1, nlev_snow
        levels(k) = REAL(k,wp)
      END DO
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_generic_snow), levels)
      DEALLOCATE(levels)


      ! Define axes for output on p- and z-levels
      !
      lwrite_pzlev = (of%name_list%pl_varlist(1) /= ' ')  .OR.  &
        &            (of%name_list%hl_varlist(1) /= ' ')
      IF (lwrite_pzlev) THEN
        !
        ! p-axis
        !
        nplev = nh_pzlev_config(of%log_patch_id)%nplev
        of%cdiZaxisID(ZA_pressure) = zaxisCreate(ZAXIS_PRESSURE, nplev)
        ALLOCATE(levels(nplev))
        DO k = 1, nplev
          levels(k) = nh_pzlev_config(of%log_patch_id)%plevels(k)
        END DO
        CALL zaxisDefLevels(of%cdiZaxisID(ZA_pressure), levels)
        CALL zaxisDefVct(of%cdiZaxisID(ZA_pressure), nplev, levels)
        DEALLOCATE(levels)

        !
        ! z-axis (height above ground)
        !
        ! DR: most probably this is wrong, because levels 
        ! nh_pzlev_config(jg)%zlevels(k) are defined as "altitude above mean 
        ! sea level"
        ! 
        nzlev = nh_pzlev_config(of%log_patch_id)%nzlev
        of%cdiZaxisID(ZA_height)  = zaxisCreate(ZAXIS_HEIGHT, nzlev)
        ALLOCATE(levels(nzlev))
        DO k = 1, nzlev
          levels(k) = nh_pzlev_config(of%log_patch_id)%zlevels(k)
        END DO
        CALL zaxisDefLevels(of%cdiZaxisID(ZA_height), levels)
        CALL zaxisDefVct(of%cdiZaxisID(ZA_height), nzlev, levels)
        DEALLOCATE(levels)
        !
        ! Altitude above mean sea level
        !
        nzlev = nh_pzlev_config(of%log_patch_id)%nzlev
        of%cdiZaxisID(ZA_altitude)  = zaxisCreate(ZAXIS_ALTITUDE, nzlev)
        ALLOCATE(levels(nzlev))
        DO k = 1, nzlev
          levels(k) = nh_pzlev_config(of%log_patch_id)%zlevels(k)
        END DO
        CALL zaxisDefLevels(of%cdiZaxisID(ZA_altitude), levels)
        CALL zaxisDefVct(of%cdiZaxisID(ZA_altitude), nzlev, levels)
        DEALLOCATE(levels)
      ENDIF

    ELSE ! oce
      of%cdiZaxisID(ZA_depth)      = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, n_zlev)
      nzlevp1 = n_zlev + 1
      of%cdiZaxisID(ZA_depth_half) = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, nzlevp1)

      ALLOCATE(levels_i(nzlevp1))
      ALLOCATE(levels_m(n_zlev))
      CALL set_zlev(levels_i, levels_m)
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_depth)     , levels_m)
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_depth_half), levels_i)
      DEALLOCATE(levels_i)
      DEALLOCATE(levels_m)
      of%cdiZaxisID(ZA_generic_ice) = zaxisCreate(ZAXIS_GENERIC, 1)
    ENDIF

    !
    ! 5. output does contain absolute time 

!PR    !
    !WRITE(6,'(a,i)') 'mode ',of%name_list%mode
    SELECT CASE (of%name_list%mode)
    CASE (1)  ! forecast mode
     of%cdiTaxisID = taxisCreate(TAXIS_RELATIVE)
     !CALL taxisDefTunit (of%cdiTaxisID, TUNIT_SECOND)
     !CALL taxisDefTunit (of%cdiTaxisID, TUNIT_MINUTE)
     IF (of%name_list%taxis_tunit > 10 .OR. of%name_list%taxis_tunit < 1 ) THEN
       of%name_list%taxis_tunit=TUNIT_HOUR
       CALL message('','invalid taxis_tunit, reset to TUNIT_HOUR')
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

  END SUBROUTINE setup_output_vlist

  !------------------------------------------------------------------------------------------------

  SUBROUTINE read_map_file(p_of, map_file)

    TYPE (t_output_file), INTENT(INOUT) :: p_of
    CHARACTER(LEN=*) :: map_file

    INTEGER :: iunit, ist, nvars
    CHARACTER(LEN=256) :: line
    CHARACTER(LEN=vname_len) :: key, val
    LOGICAL :: valid

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/read_map_file'
    TYPE (t_keyword_list), POINTER :: keywords => NULL()
    CHARACTER(len=MAX_STRING_LEN)  :: filename

    CALL associate_keyword("<path>", TRIM(model_base_dir), keywords)
    filename = TRIM(with_keywords(keywords, map_file))

    iunit = find_next_free_unit(10,99)
    OPEN (unit=iunit,file=filename,access='SEQUENTIAL', &
      &  form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)

    IF(ist/=0)THEN
      CALL finish (routine, 'Open of map_file '//TRIM(map_file)//' failed')
    ENDIF

    nvars = 0

    DO
      line = ' '
      READ(iunit,'(a)',iostat=ist) line
      IF(ist < 0) EXIT ! No more lines
      IF(ist > 0) CALL finish(routine, 'Read error in map_file '//TRIM(map_file))
      CALL parse_line(line,key,val,valid)
      IF(valid) nvars = nvars+1
    ENDDO

    REWIND(iunit)
    ALLOCATE(p_of%name_map(2,nvars))

    nvars = 0

    DO
      line = ' '
      READ(iunit,'(a)',iostat=ist) line
      IF(ist < 0) EXIT ! No more lines
      IF(ist > 0) CALL finish(routine, 'Read error in map_file '//TRIM(map_file))
      CALL parse_line(line,key,val,valid)
      IF(valid) THEN
        nvars = nvars+1
        ! Safety only
        IF(nvars > UBOUND(p_of%name_map,2)) &
          CALL finish(routine,'INTERNAL: Length of name_map calculated wrong')
        p_of%name_map(1,nvars) = key
        p_of%name_map(2,nvars) = val
      ENDIF
    ENDDO

    CLOSE(unit=iunit)

  END SUBROUTINE read_map_file

  !------------------------------------------------------------------------------------------------
  ! parse_line:
  ! Parses a line from map_file and returns key and value part
  ! (first two nonblank words separated by blanks or tabs)

  SUBROUTINE parse_line(line, key, val, valid)

    CHARACTER(LEN=*), INTENT(IN)  :: line
    CHARACTER(LEN=*), INTENT(OUT) :: key, val
    LOGICAL, INTENT(OUT)          :: valid

    INTEGER :: ipos1, ipos2

    valid = .FALSE.
    key = ' '
    val = ' '

    ! Search first nonblank character
    DO ipos1 = 1, LEN(line)
      IF(.NOT.isblank(line(ipos1:ipos1))) EXIT
    ENDDO

    IF(ipos1 > LEN(line)) RETURN ! completely empty line

    IF(line(ipos1:ipos1) == '#') RETURN ! comment line

    ! Search end of key
    DO ipos2 = ipos1, LEN(line)
      IF(isblank(line(ipos2:ipos2))) EXIT
    ENDDO

    key = line(ipos1:ipos2-1)

    ! Search next nonblank character
    DO ipos1 = ipos2, LEN(line)
      IF(.NOT.isblank(line(ipos1:ipos1))) EXIT
    ENDDO

    IF(ipos1 > LEN(line)) THEN ! line contains no value part
      CALL message('mo_name_list_output','Illegal line in name_map:')
      CALL message('mo_name_list_output',TRIM(line))
      RETURN
    ENDIF

    ! Search end of val
    DO ipos2 = ipos1, LEN(line)
      IF(isblank(line(ipos2:ipos2))) EXIT
    ENDDO

    val = line(ipos1:ipos2-1)

    valid = .TRUE.

    CONTAINS

    ! Fortran equivalent to C isblank function
    LOGICAL FUNCTION isblank(c)
      CHARACTER(LEN=1) :: c
      IF(c==' ' .OR. ICHAR(c)==9) THEN
        ! Blank or Tab
        isblank = .TRUE.
      ELSE
        isblank = .FALSE.
      ENDIF
    END FUNCTION isblank

  END SUBROUTINE parse_line

  !------------------------------------------------------------------------------------------------
  !
  ! Sets the grid information in output file
  !
  SUBROUTINE set_grid_info(of)

    TYPE (t_output_file), INTENT(IN) :: of

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/set_grid_info'
    INTEGER :: idom

    idom  = of%phys_patch_id

    ! cell grid
    CALL gridDefXvals  (of%cdiCellGridID, patch_info(idom)%grid_c%lon)
    CALL gridDefYvals  (of%cdiCellGridID, patch_info(idom)%grid_c%lat)
    CALL gridDefXbounds(of%cdiCellGridID, patch_info(idom)%grid_c%lonv)
    CALL gridDefYbounds(of%cdiCellGridID, patch_info(idom)%grid_c%latv)

    ! edge grid
    CALL gridDefXvals  (of%cdiEdgeGridID, patch_info(idom)%grid_e%lon)
    CALL gridDefYvals  (of%cdiEdgeGridID, patch_info(idom)%grid_e%lat)
    CALL gridDefXbounds(of%cdiEdgeGridID, patch_info(idom)%grid_e%lonv)
    CALL gridDefYbounds(of%cdiEdgeGridID, patch_info(idom)%grid_e%latv)

    ! vertex grid
    CALL gridDefXvals  (of%cdiVertGridID, patch_info(idom)%grid_v%lon)
    CALL gridDefYvals  (of%cdiVertGridID, patch_info(idom)%grid_v%lat)
    CALL gridDefXbounds(of%cdiVertGridID, patch_info(idom)%grid_v%lonv)
    CALL gridDefYbounds(of%cdiVertGridID, patch_info(idom)%grid_v%latv)

  END SUBROUTINE set_grid_info


  !------------------------------------------------------------------------------------------------
  !
  ! Copies the grid information from grid file to output file
  !
    SUBROUTINE copy_grid_info(of)

    TYPE (t_output_file), INTENT(IN) :: of

    INTEGER :: ncid, dimid, varid
    INTEGER :: i_nc, i_ne, i_nv
    INTEGER :: i_dom

    REAL(wp), ALLOCATABLE :: clon(:), clat(:), clonv(:,:), clatv(:,:)
    REAL(wp), ALLOCATABLE :: elon(:), elat(:), elonv(:,:), elatv(:,:)
    REAL(wp), ALLOCATABLE :: vlon(:), vlat(:), vlonv(:,:), vlatv(:,:)

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/copy_grid_info'

    i_dom = of%phys_patch_id

    ! Please note: The following is more or less a copy from mo_io_vlist with adaptions
    ! to the data structures used here.
    ! Unfortunately it seems necessary to open the gridfile for reading the information
    ! since it is not read and stored during patch input.

    !---------------------------------------------------------------------------
    ! Open grid file, read dimensions and make a cross check if they match.
    ! This is just for safety and could be skipped, of course.

    CALL nf(nf_open(TRIM(patch_info(i_dom)%grid_filename), NF_NOWRITE, ncid))
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
    CASE (6)
      CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
    END SELECT
    CALL nf(nf_inq_dimlen(ncid, dimid, i_nc))
    !
    CALL nf(nf_inq_dimid(ncid, 'edge', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, i_ne))
    !
    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
    CASE (6)
      CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
    END SELECT
    CALL nf(nf_inq_dimlen(ncid, dimid, i_nv))

    IF(i_nc /= patch_info(i_dom)%cells%n_log) &
      CALL finish(routine,'Number of cells differs in '//TRIM(patch_info(i_dom)%grid_filename))
    IF(i_ne /= patch_info(i_dom)%edges%n_log) &
      CALL finish(routine,'Number of edges differs in '//TRIM(patch_info(i_dom)%grid_filename))
    IF(i_nv /= patch_info(i_dom)%verts%n_log) &
      CALL finish(routine,'Number of verts differs in '//TRIM(patch_info(i_dom)%grid_filename))
    !
    !---------------------------------------------------------------------------
    ! cell grid

    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_varid(ncid, 'clon', varid))
    CASE (6)
      CALL nf(nf_inq_varid(ncid, 'vlon', varid))
    END SELECT
    
    ALLOCATE(clon(i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clon))
    CALL reorder1(patch_info(i_dom)%cells%n_glb, patch_info(i_dom)%cells%log_dom_index, clon)
    CALL gridDefXvals(of%cdiCellGridID, clon)
    DEALLOCATE(clon)


    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_varid(ncid, 'clat', varid))
    CASE (6)
      CALL nf(nf_inq_varid(ncid, 'vlat', varid))
    END SELECT
    
    ALLOCATE(clat(i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clat))
    CALL reorder1(patch_info(i_dom)%cells%n_glb, patch_info(i_dom)%cells%log_dom_index, clat)

    CALL gridDefYvals(of%cdiCellGridID, clat)
    DEALLOCATE(clat)
    

    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_varid(ncid, 'clon_vertices', varid))
    CASE (6)
      CALL nf(nf_inq_varid(ncid, 'vlon_vertices', varid))
    END SELECT
    
    ALLOCATE(clonv(global_cell_type, i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clonv))
    CALL reorder2(patch_info(i_dom)%cells%n_glb, patch_info(i_dom)%cells%log_dom_index, clonv)

    CALL gridDefXbounds(of%cdiCellGridID, clonv)
    DEALLOCATE(clonv)


    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_varid(ncid, 'clat_vertices', varid))
    CASE (6)
      CALL nf(nf_inq_varid(ncid, 'vlat_vertices', varid))
    END SELECT

    ALLOCATE(clatv(global_cell_type, i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clatv))
    CALL reorder2(patch_info(i_dom)%cells%n_glb, patch_info(i_dom)%cells%log_dom_index, clatv)

    CALL gridDefYbounds(of%cdiCellGridID, clatv)
    DEALLOCATE(clatv)

    !-------------------------------------------------------------------------
    ! edge grid

    ALLOCATE(elon(i_ne))
    CALL nf(nf_inq_varid(ncid, 'elon', varid))
    CALL nf(nf_get_var_double(ncid, varid, elon))
    CALL reorder1(patch_info(i_dom)%edges%n_glb, patch_info(i_dom)%edges%log_dom_index, elon)

    CALL gridDefXvals(of%cdiEdgeGridID, elon)
    DEALLOCATE(elon)

    ALLOCATE(elat(i_ne))
    CALL nf(nf_inq_varid(ncid, 'elat', varid))
    CALL nf(nf_get_var_double(ncid, varid, elat))
    CALL reorder1(patch_info(i_dom)%edges%n_glb, patch_info(i_dom)%edges%log_dom_index, elat)

    CALL gridDefYvals(of%cdiEdgeGridID, elat)
    DEALLOCATE(elat)

    ALLOCATE(elonv(4, i_ne))
    CALL nf(nf_inq_varid(ncid, 'elon_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, elonv))
    CALL reorder2(patch_info(i_dom)%edges%n_glb, patch_info(i_dom)%edges%log_dom_index, elonv)

    CALL gridDefXbounds(of%cdiEdgeGridID, elonv)
    DEALLOCATE(elonv)

    ALLOCATE(elatv(4, i_ne))
    CALL nf(nf_inq_varid(ncid, 'elat_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, elatv))
    CALL reorder2(patch_info(i_dom)%edges%n_glb, patch_info(i_dom)%edges%log_dom_index, elatv)

    CALL gridDefYbounds(of%cdiEdgeGridID, elatv)
    DEALLOCATE(elatv)

    !-------------------------------------------------------------------------
    ! vertex grid

    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_varid(ncid, 'vlon', varid))
    CASE (6)
      CALL nf(nf_inq_varid(ncid, 'clon', varid))
    END SELECT

    ALLOCATE(vlon(i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlon))
    CALL reorder1(patch_info(i_dom)%verts%n_glb, patch_info(i_dom)%verts%log_dom_index, vlon)

    CALL gridDefXvals(of%cdiVertGridID, vlon)
    DEALLOCATE(vlon)

    SELECT CASE (global_cell_type)
    CASE (3)
      CALL nf(nf_inq_varid(ncid, 'vlat', varid))
    CASE (6)
      CALL nf(nf_inq_varid(ncid, 'clat', varid))
    END SELECT

    ALLOCATE(vlat(i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlat))
    CALL reorder1(patch_info(i_dom)%verts%n_glb, patch_info(i_dom)%verts%log_dom_index, vlat)

    CALL gridDefYvals(of%cdiVertGridID, vlat)
    DEALLOCATE(vlat)

    IF(global_cell_type==3) THEN
      CALL nf(nf_inq_varid(ncid, 'vlon_vertices', varid))
    ELSEIF(global_cell_type==6) THEN
      CALL nf(nf_inq_varid(ncid, 'clon_vertices', varid))
    ENDIF

    ALLOCATE(vlonv(9-global_cell_type, i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlonv))
    CALL reorder2(patch_info(i_dom)%verts%n_glb, patch_info(i_dom)%verts%log_dom_index, vlonv)

    CALL gridDefXbounds(of%cdiVertGridID, vlonv)
    DEALLOCATE(vlonv)

    IF(global_cell_type==3) THEN
      CALL nf(nf_inq_varid(ncid, 'vlat_vertices', varid))
    ELSEIF(global_cell_type==6) THEN
      CALL nf(nf_inq_varid(ncid, 'clat_vertices', varid))
    ENDIF

    ALLOCATE(vlatv(9-global_cell_type, i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlatv))
    CALL reorder2(patch_info(i_dom)%verts%n_glb, patch_info(i_dom)%verts%log_dom_index, vlatv)

    CALL gridDefYbounds(of%cdiVertGridID, vlatv)
    DEALLOCATE(vlatv)

    !-------------------------------------------------------------------------

    ! Close NetCDF file, it is not needed any more
    CALL nf(nf_close(ncid))

  CONTAINS

    SUBROUTINE nf(status)
      INTEGER, INTENT(in) :: status

      IF (status /= nf_noerr) THEN
        CALL finish(routine, 'NetCDF Error: '//nf_strerror(status))
      ENDIF
    END SUBROUTINE nf

    ! reorder1: get the physical patch points from the logical patch
    ! Note that this works within the array as long as idx is monotonically increasing
    SUBROUTINE reorder1(n, idx, array)
      INTEGER, INTENT(IN)     :: n, idx(:)
      REAL(wp), INTENT(INOUT) :: array(:)
      INTEGER :: i

      DO i = 1, n
        array(i) = array(idx(i))
      ENDDO
    END SUBROUTINE reorder1

    ! reorder2: same as reorder1 for 2D array
    SUBROUTINE reorder2(n, idx, array)
      INTEGER, INTENT(IN)     :: n, idx(:)
      REAL(wp), INTENT(INOUT) :: array(:,:)
      INTEGER :: i

      DO i = 1, n
        array(:,i) = array(:,idx(i))
      ENDDO
    END SUBROUTINE reorder2

  END SUBROUTINE copy_grid_info


  !------------------------------------------------------------------------------------------------
  !
  ! define variables and attributes
  !
  SUBROUTINE add_variables_to_vlist(of)

    TYPE (t_output_file), INTENT(IN), TARGET :: of

    TYPE (t_var_metadata), POINTER :: info
    !
    INTEGER :: iv, vlistID, varID, gridID, zaxisID, nlev, nlevp1, znlev_soil, i
    CHARACTER(LEN=vname_len) :: mapped_name

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/add_variables_to_vlist'

    vlistID = of%cdiVlistID

    nlev   = num_lev(of%log_patch_id)
    nlevp1 = num_levp1(of%log_patch_id)

    ! See above ...
    znlev_soil = SIZE(zml_soil)-1
    !
    DO iv = 1, of%num_vars
      !
      info => of%var_desc(iv)%info
      !
      ! set grid ID
      !
      SELECT CASE (info%hgrid)
      CASE(GRID_UNSTRUCTURED_CELL)
        info%cdiGridID = of%cdiCellGridID
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
      SELECT CASE (info%vgrid)

      CASE (ZAXIS_SURFACE)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_surface)

      CASE (ZAXIS_MEANSEA)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_meansea)

      CASE (ZAXIS_HYBRID)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_hybrid)

      CASE (ZAXIS_HYBRID_HALF)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_hybrid_half)

      CASE (ZAXIS_DEPTH_BELOW_SEA)
        ! RJ: Not sure about this ...
        IF (info%used_dimensions(2) == n_zlev) THEN
          info%cdiZaxisID =  of%cdiZaxisID(ZA_depth)
        ELSE IF (info%used_dimensions(2) == n_zlev+1) THEN
          info%cdiZaxisID =  of%cdiZaxisID(ZA_depth_half)
        ELSE
          WRITE (message_text,'(a,a,a,i0)') &
               &  'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
          CALL message('',message_text)
          CALL finish(routine,'Dimension mismatch for ZAXIS_DEPTH_BELOW_SEA')
        ENDIF

      CASE (ZAXIS_DEPTH_BELOW_LAND)
        IF (info%used_dimensions(2) == znlev_soil+1) THEN
          info%cdiZaxisID =  of%cdiZaxisID(ZA_depth_below_land)
        ELSE IF (info%used_dimensions(2) == znlev_soil+2) THEN
          info%cdiZaxisID =  of%cdiZaxisID(ZA_depth_below_land_p1)
        ELSE
          WRITE (message_text,'(a,a,a,i0)') &
               &  'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
          CALL message('',message_text)
          CALL finish(routine,'Dimension mismatch for ZAXIS_DEPTH_BELOW_LAND')
        ENDIF

      CASE (ZAXIS_HEIGHT)
          ! In all other cases, ZAXIS_HEIGHT seems to be equivalent to ZAXIS_HYBRID/ZAXIS_HYBRID_HALF
          ! TODO: Is there a difference to ZAXIS_HYBRID/ZAXIS_HYBRID_HALF ???
          IF (info%used_dimensions(2) == nlevp1) THEN
            info%cdiZaxisID =  of%cdiZaxisID(ZA_hybrid_half)
          ELSE IF (info%used_dimensions(2) == nlev) THEN
            info%cdiZaxisID =  of%cdiZaxisID(ZA_hybrid)
          ELSE
            WRITE (message_text,'(a,a,a,i0)') &
                 &  'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
            CALL message('',message_text)
            CALL finish(routine,'Dimension mismatch for ZAXIS_HEIGHT')
          ENDIF

      CASE (ZAXIS_GENERIC)
        IF(info%used_dimensions(2) == nlev_snow) THEN
          info%cdiZaxisID =  of%cdiZaxisID(ZA_generic_snow)
        ELSE IF(info%used_dimensions(2) == nlev_snow+1) THEN
          info%cdiZaxisID =  of%cdiZaxisID(ZA_generic_snow_p1)
        ELSE
          WRITE (message_text,'(a,a,a,i0)') &
               &  'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
          CALL message('',message_text)
          CALL finish(routine,'Dimension mismatch for ZAXIS_GENERIC')
        ENDIF

      CASE (ZAXIS_PRESSURE)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_pressure)

      CASE (ZAXIS_ALTITUDE)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_altitude)


      CASE DEFAULT
        WRITE (message_text,'(a,a,a,i0)') &
             &  'Variable: ',TRIM(info%name),' ZAXIS: ',info%vgrid
        CALL message('',message_text)
        CALL finish(routine, 'ZAXIS definition missing for '//TRIM(info%name))

      END SELECT

      zaxisID = info%cdiZaxisID
      ! Search name mapping for name in NetCDF file

      mapped_name = info%name
      IF(ALLOCATED(of%name_map)) THEN
        DO i = 1, UBOUND(of%name_map, 2)
          IF(TRIM(info%name) == TRIM(of%name_map(1,i))) THEN
            mapped_name = of%name_map(2,i)
            EXIT
          ENDIF
        ENDDO
      ENDIF

      varID = vlistDefVar(vlistID, gridID, zaxisID, info%cdiTimeID)
      
      info%cdiVarID   = varID

      CALL vlistDefVarDatatype(vlistID, varID, DATATYPE_FLT32)

      CALL vlistDefVarName(vlistID, varID, TRIM(mapped_name))

      IF (info%cf%long_name /= '') CALL vlistDefVarLongname(vlistID, varID, info%cf%long_name)
      IF (info%cf%units /= '') CALL vlistDefVarUnits(vlistID, varID, info%cf%units)

      ! Currently only real valued variables are allowed, so we can always use info%missval%rval
      IF (info%lmiss) CALL vlistDefVarMissval(vlistID, varID, info%missval%rval)


      IF ( of%output_type == FILETYPE_GRB2 ) THEN
        ! Set GRIB2 Triplet
        CALL vlistDefVarParam(vlistID, varID,                                              &
          &  cdiEncodeParam(info%grib2%number, info%grib2%category, info%grib2%discipline) )
      ENDIF

!DR
!DR Still missing: Set typeOfStatisticalProcessing
!DR This feature is not yet fully supported by CDI
!DR
!DR      CALL  vlistDefVarTsteptype(vlistID, varID, info%istatproc);
        
    ENDDO
    !
  END SUBROUTINE add_variables_to_vlist

  !------------------------------------------------------------------------------------------------
  !> open_output_file:
  !! Opens a output file and sets its vlist

  SUBROUTINE open_output_file(of, jfile)

    TYPE(t_output_file), INTENT(INOUT) :: of
    INTEGER, INTENT(IN) :: jfile ! Number of file set to open

    CHARACTER(LEN=16) :: extn
    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/open_output_file'

    ! Please note that this routine is only executed on one processor (for a specific file)
    ! and thus all calls to message get the all_print=.TRUE. argument so that the messages
    ! really appear in the log

    !
    ! check output file type
    !
    SELECT CASE (of%output_type)       
    CASE (FILETYPE_NC)
      CALL finish(routine,'netCDF classic not supported')
    CASE (FILETYPE_NC2, FILETYPE_NC4)
      ! this is ok, both formats can write more than 2GB files
      extn = '.nc'
    CASE (FILETYPE_GRB)
      CALL finish(routine,'GRIB1 not supported')
    CASE (FILETYPE_GRB2)
      CALL message(routine,'GRIB2 support experimental',all_print=.TRUE.)
      extn = '.grb'
    CASE default
      CALL finish(routine,'unknown output_type')
    END SELECT

    ! Set actual output file name and open file

    IF(my_process_is_mpi_test()) THEN
      WRITE(of%filename,'(a,"_",i4.4,"_TEST",a)') TRIM(of%filename_pref),jfile,TRIM(extn)
    ELSE
      WRITE(of%filename,'(a,"_",i4.4,a)') TRIM(of%filename_pref),jfile,TRIM(extn)
    ENDIF
    of%cdiFileID = streamOpenWrite(TRIM(of%filename), of%output_type)

    IF (of%cdiFileID < 0) THEN
      WRITE(message_text,'(a)') cdiStringError(of%cdiFileID)
      CALL message('',message_text,all_print=.TRUE.)
      CALL finish (routine, 'open failed on '//TRIM(of%filename))
    ELSE
      CALL message (routine, 'opened '//TRIM(of%filename),all_print=.TRUE.)
    END IF

    ! assign the vlist (which must have ben set before)
!PR 
    !WRITE(6,'(a,i)') 'taxisTunit openning ',taxisInqTunit(of%cdiTaxisID)

    CALL streamDefVlist(of%cdiFileID, of%cdiVlistID)

    ! set cdi internal time index to 0 for writing time slices in netCDF

    of%cdiTimeIndex = 0

  END SUBROUTINE open_output_file

  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE close_name_list_output
    !
    ! Close all name_list files
    !
    INTEGER :: i, jg


#ifndef NOMPI
    IF(use_async_name_list_io.AND..NOT.my_process_is_io().AND..NOT.my_process_is_mpi_test()) THEN

      ! write recent samples of meteogram output
      DO jg = 1, n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          CALL meteogram_flush_file(jg)
        END IF
      END DO

      CALL compute_wait_for_async_io()
      CALL compute_shutdown_async_io()

    ELSE
#endif
      DO i = 1, SIZE(output_file)
        CALL close_output_file(output_file(i))
        CALL destroy_output_vlist(output_file(i))
      ENDDO
#ifndef NOMPI
    ENDIF
#endif

    DEALLOCATE(output_file)

  END SUBROUTINE close_name_list_output

  !------------------------------------------------------------------------------------------------
  !
  ! Close output stream and the associated file,
  ! destroy all vlist related data for this file
  !
  SUBROUTINE close_output_file(of)

    TYPE (t_output_file), INTENT(INOUT) :: of

    CHARACTER(LEN=*), PARAMETER :: routine = &
      & TRIM('mo_name_list_output/close_output_file')

    IF(of%cdiFileID /= CDI_UNDEFID) CALL streamClose(of%cdiFileID)

    of%cdiFileID = CDI_UNDEFID

  END SUBROUTINE close_output_file

  !------------------------------------------------------------------------------------------------
  !
  ! Close output stream and the associated file,
  ! destroy all vlist related data for this file
  !
  SUBROUTINE destroy_output_vlist(of)

    TYPE (t_output_file), INTENT(INOUT) :: of

    INTEGER :: j, vlistID


    vlistID = of%cdiVlistID
    IF(vlistID /= CDI_UNDEFID) THEN
      IF(of%cdiCellGridID   /= CDI_UNDEFID) CALL gridDestroy(of%cdiCellGridID)
      IF(of%cdiEdgeGridID   /= CDI_UNDEFID) CALL gridDestroy(of%cdiEdgeGridID)
      IF(of%cdiVertGridID   /= CDI_UNDEFID) CALL gridDestroy(of%cdiVertGridID)
      IF(of%cdiLonLatGridID /= CDI_UNDEFID) CALL gridDestroy(of%cdiLonLatGridID)
      CALL taxisDestroy(of%cdiTaxisID)
      DO j = 1, SIZE(of%cdiZaxisID)
        IF(of%cdiZaxisID(j) /= CDI_UNDEFID) CALL zaxisDestroy(of%cdiZaxisID(j))
      ENDDO
      CALL vlistDestroy(vlistID)
    ENDIF

    of%cdiVlistID      = CDI_UNDEFID
    of%cdiCellGridID   = CDI_UNDEFID
    of%cdiEdgeGridID   = CDI_UNDEFID
    of%cdiVertGridID   = CDI_UNDEFID
    of%cdiLonLatGridID = CDI_UNDEFID
    of%cdiTaxisID      = CDI_UNDEFID
    of%cdiZaxisID(:)   = CDI_UNDEFID

  END SUBROUTINE destroy_output_vlist

  !------------------------------------------------------------------------------------------------
  !
  ! Loop over all output_name_list's, write the ones for which output is due
  ! This routine also cares about opening the output files the first time
  ! and reopening the files after a certain number of steps.
  !
  SUBROUTINE write_name_list_output(datetime, sim_time, last_step)

    TYPE(t_datetime), INTENT(in) :: datetime
    REAL(wp), INTENT(in)         :: sim_time
    LOGICAL, INTENT(IN)          :: last_step

    INTEGER :: i, idate, itime, iret, n, jg
    TYPE(t_output_name_list), POINTER :: p_onl
    CHARACTER(LEN=filename_max+100) :: text
    REAL(wp), PARAMETER :: eps = 1.d-10 ! Tolerance for checking output bounds
    LOGICAL :: lnewly_initialized = .FALSE.

    ! If asynchronous I/O is enabled, the compute PEs have to make sure
    ! that the I/O PEs are ready with the last output step before
    ! writing data into the I/O memory window.
    ! This routine (write_name_list_output) is also called from the I/O PEs,
    ! but in this case the calling routine cares about the flow control.

#ifndef NOMPI
    IF(use_async_name_list_io) THEN
      IF(.NOT.my_process_is_io().AND..NOT.my_process_is_mpi_test()) THEN
        ! write recent samples of meteogram output
        DO jg = 1, n_dom
          IF (meteogram_output_config(jg)%lenabled) THEN
            CALL meteogram_flush_file(jg)
          END IF
        END DO
        CALL compute_wait_for_async_io()
      END IF
    ENDIF
#endif

    ! Check if files have to be (re)opened

    ! Go over all output files
    DO i = 1, SIZE(output_file)

      p_onl => output_file(i)%name_list

      ! Check if output is due for this file
      IF (is_output_file_active(output_file(i), sim_time, dtime, i_sample, last_step)) THEN

        IF (output_file(i)%io_proc_id == p_pe) THEN
          IF (.NOT. output_file(i)%initialized) THEN
            CALL setup_output_vlist(output_file(i))
            lnewly_initialized = .TRUE.
          ELSE
            lnewly_initialized = .FALSE.
          ENDIF
        ENDIF

        IF (lnewly_initialized .OR. MOD(p_onl%n_output_steps,p_onl%steps_per_file) == 0) THEN
          IF (output_file(i)%io_proc_id == p_pe) THEN
            IF(.NOT. lnewly_initialized) THEN
              CALL close_output_file(output_file(i))
            ENDIF
            CALL open_output_file(output_file(i),p_onl%n_output_steps/p_onl%steps_per_file+1)
          ENDIF
        ENDIF

      ENDIF

    ENDDO

    ! Do the output

    idate = cdiEncodeDate(datetime%year, datetime%month, datetime%day)
    itime = cdiEncodeTime(datetime%hour, datetime%minute, NINT(datetime%second))

    ! Go over all output files
    DO i = 1, SIZE(output_file)

      ! Check if output is due for this file
      IF (is_output_file_active(output_file(i), sim_time, dtime, i_sample, last_step)) THEN

        IF (output_file(i)%io_proc_id == p_pe) THEN
          CALL taxisDefVdate(output_file(i)%cdiTaxisID, idate)
          CALL taxisDefVtime(output_file(i)%cdiTaxisID, itime)
!PR
    !WRITE(6,'(a,i,a,i)')'idate forc and ref',idate,' ',taxisInqRdate(output_file(i)%cdiTaxisID)
    !WRITE(6,'(a,i,a,i)')'itime forc and ref',itime,' ',taxisInqRtime(output_file(i)%cdiTaxisID)
    !WRITE(6,'(a,i,i)')'validation date and time ',taxisInqVdate(output_file(i)%cdiTaxisID), &
    !             & taxisInqVtime(output_file(i)%cdiTaxisID)
 
    !WRITE(6,'(a,i)') 'taxisTunit writing ',taxisInqTunit(output_file(i)%cdiTaxisID)

!PR
          iret = streamDefTimestep(output_file(i)%cdiFileId, output_file(i)%cdiTimeIndex)
          output_file(i)%cdiTimeIndex = output_file(i)%cdiTimeIndex + 1

          ! Notify user

#ifndef __SX__
          WRITE(text,'(a,a,a,1pg15.9,a,i6)') &
            'Output to ',TRIM(output_file(i)%filename),' at simulation time ',sim_time, &
             ' by PE ',p_pe
          CALL message('mo_name_list_output',text,all_print=.TRUE.)
#endif

        ENDIF

        IF(my_process_is_io()) THEN
#ifndef NOMPI
          IF(output_file(i)%io_proc_id == p_pe) THEN
            CALL io_proc_write_name_list(output_file(i), &
              &          (MOD(p_onl%n_output_steps,p_onl%steps_per_file) == 0) )
          ENDIF
#endif
        ELSE
          CALL write_name_list(output_file(i), &
            &            (MOD(p_onl%n_output_steps,p_onl%steps_per_file) == 0) )
        ENDIF
      ENDIF

    ENDDO

    ! Increase next_output_time for all output name lists which have been written.
    ! Please note that more than 1 output_file may point to the same name_list,
    ! so this must not be done in the loop above!

    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT

      IF (is_output_nml_active(p_onl, sim_time, dtime, i_sample, last_step)) THEN
        p_onl%n_output_steps = p_onl%n_output_steps + 1
      ENDIF

      ! Switch all name lists to next output time.
      ! This is done in a DO WHILE loop to catch the case
      ! where an output increment less than the time step
      ! or two output_bounds triples which are too close are specified.

      DO WHILE (is_output_nml_active(p_onl, sim_time, dtime, i_sample))
        n = p_onl%cur_bounds_triple
        IF(p_onl%next_output_time + p_onl%output_bounds(3,n) <= p_onl%output_bounds(2,n)+eps) THEN
          ! Next output time will be within current bounds triple
          p_onl%next_output_time = p_onl%next_output_time + p_onl%output_bounds(3,n)
        ELSE
          ! the last bounds triple is always set to 0, so no need to check if n >= max_bounds
          IF(p_onl%output_bounds(3,n+1) > 0._wp) THEN
            ! Start to use next bounds triple
            p_onl%next_output_time = p_onl%output_bounds(1,n+1)
          ELSE
            ! Output is finished for this name list
            p_onl%next_output_time = HUGE(0._wp)
          ENDIF
          p_onl%cur_bounds_triple = p_onl%cur_bounds_triple + 1
        ENDIF
      ENDDO

      p_onl => p_onl%next

    ENDDO

    ! If asynchronous I/O is enabled, the compute PEs can now start the I/O PEs
#ifndef NOMPI
    IF(use_async_name_list_io) THEN
      IF(.NOT.my_process_is_io().AND..NOT.my_process_is_mpi_test()) &
        CALL compute_start_async_io(datetime, sim_time, last_step)
    ENDIF
#endif

    ! Close output file when the related model domain has stopped execution
    DO i = 1, SIZE(output_file)
      IF (sim_time > output_file(i)%end_time) CALL close_output_file(output_file(i))
    ENDDO

  END SUBROUTINE write_name_list_output

  !------------------------------------------------------------------------------------------------
  !
  ! Write a output name list
  !
  SUBROUTINE write_name_list(of, l_first_write)

#ifndef NOMPI
    USE mpi, ONLY: MPI_LOCK_EXCLUSIVE, MPI_MODE_NOCHECK
#endif

    TYPE (t_output_file), INTENT(INOUT), TARGET :: of
    LOGICAL, INTENT(IN) :: l_first_write

    INTEGER :: tl, i_dom, i_log_dom, i, iv, jk, n_points, nlevs, nblks, &
      &        nindex, mpierr, lonlat_id, ierrstat
    INTEGER(i8) :: ioff
    TYPE (t_var_metadata), POINTER :: info
    TYPE(t_reorder_info), POINTER  :: p_ri

    REAL(wp), POINTER :: r_ptr(:,:,:)
    REAL(wp), ALLOCATABLE :: r_tmp(:,:,:), r_out_recv(:,:)
    REAL(sp), ALLOCATABLE :: r_out_sp(:,:)
    REAL(dp), ALLOCATABLE :: r_out_dp(:,:)
    TYPE(t_comm_pattern), POINTER :: p_pat

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/write_name_list'

    ! Offset in memory window for async I/O
    ioff = of%my_mem_win_off

    i_dom = of%phys_patch_id
    i_log_dom = of%log_patch_id

    tl = 0 ! to prevent warning

#ifndef NOMPI
    ! In case of async IO: Lock own window before writing to it
    IF(use_async_name_list_io .AND. .NOT.my_process_is_mpi_test()) &
      CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p_pe_work, MPI_MODE_NOCHECK, mpi_win, mpierr)
#endif

    ! Go over all name list variables for this output file

    DO iv = 1, of%num_vars

      info => of%var_desc(iv)%info

      ! inspect time-constant variables only if we are writing the
      ! first step in this file:
      IF ((info%cdiTimeID == TIME_CONSTANT) .AND. .NOT. l_first_write) CYCLE

      ! Check if first dimension of array is nproma.
      ! Otherwise we got an array which is not suitable for this output scheme.

      IF(info%used_dimensions(1) /= nproma) &
        CALL finish(routine,'1st dim is not nproma: '//TRIM(info%name))

      r_ptr => NULL()

      ! For time level dependent elements: set time level and check if time level is present

      IF(.NOT.ASSOCIATED(of%var_desc(iv)%r_ptr)) THEN
        SELECT CASE (info%tlev_source)
          CASE(0); tl = nnow(i_log_dom)
          CASE(1); tl = nnow_rcf(i_log_dom)
          CASE DEFAULT
            CALL finish(routine,'Unsupported tlev_source')
        END SELECT
        IF(tl<=0 .OR. tl>max_time_levels) &
          CALL finish(routine, 'Illegal time level in nnow()/nnow_rcf()')
        ! Check if present
        IF(.NOT.ASSOCIATED(of%var_desc(iv)%tlev_ptr(tl)%p)) &
          CALL finish(routine,'Actual timelevel not in '//TRIM(info%name))
      ENDIF
      !
      IF (info%lcontained) THEN 
        nindex = info%ncontained
      ELSE
        nindex = 1
      ENDIF

      SELECT CASE (info%ndims)
      CASE (1)
        CALL message(routine, info%name)
        CALL finish(routine,'1d arrays not handled yet.')
      CASE (2)
        ! Make a 3D copy of the array
        ALLOCATE(r_ptr(info%used_dimensions(1),1,info%used_dimensions(2)))
        IF(ASSOCIATED(of%var_desc(iv)%r_ptr)) THEN
          r_ptr(:,1,:) = of%var_desc(iv)%r_ptr(:,:,nindex,1,1)
        ELSE
          r_ptr(:,1,:) = of%var_desc(iv)%tlev_ptr(tl)%p(:,:,nindex,1,1)
        ENDIF
      CASE (3)
        ! Just set a pointer to the array
        IF(ASSOCIATED(of%var_desc(iv)%r_ptr)) THEN
          r_ptr => of%var_desc(iv)%r_ptr(:,:,:,nindex,1)
        ELSE
          r_ptr => of%var_desc(iv)%tlev_ptr(tl)%p(:,:,:,nindex,1)
        ENDIF
      CASE (4)
        CALL message(routine, info%name)
        CALL finish(routine,'4d arrays not handled yet.')
      CASE (5)
        CALL message(routine, info%name)
        CALL finish(routine,'5d arrays not handled yet.')
      CASE DEFAULT 
        CALL message(routine, info%name)
        CALL finish(routine,'dimension not set.')        
      END SELECT

      IF(info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = info%used_dimensions(2)
      ENDIF

      ! Get pointer to appropriate reorder_info

      SELECT CASE (info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        p_ri => patch_info(i_dom)%cells
        IF(l_output_phys_patch) THEN
          p_pat => p_phys_patch(i_dom)%comm_pat_gather_c
        ELSE
          p_pat => p_patch(i_dom)%comm_pat_gather_c
        ENDIF
      CASE (GRID_UNSTRUCTURED_EDGE)
        p_ri => patch_info(i_dom)%edges
        IF(l_output_phys_patch) THEN
          p_pat => p_phys_patch(i_dom)%comm_pat_gather_e
        ELSE
          p_pat => p_patch(i_dom)%comm_pat_gather_e
        ENDIF
      CASE (GRID_UNSTRUCTURED_VERT)
        p_ri => patch_info(i_dom)%verts
        IF(l_output_phys_patch) THEN
          p_pat => p_phys_patch(i_dom)%comm_pat_gather_v
        ELSE
          p_pat => p_patch(i_dom)%comm_pat_gather_v
        ENDIF
      CASE (GRID_REGULAR_LONLAT)
        lonlat_id = info%hor_interp%lonlat_id
        p_ri  => lonlat_info(lonlat_id, i_log_dom)
        p_pat => lonlat_grid_list(lonlat_id)%p_pat(i_log_dom)
      CASE default
        CALL finish(routine,'unknown grid type')
      END SELECT

      IF(.NOT.use_async_name_list_io .OR. my_process_is_mpi_test()) THEN

        ! No asynchronous I/O, gather the array on stdio PE and write it out there

        n_points = p_ri%n_glb
        nblks = (n_points-1)/nproma + 1

        IF(my_process_is_mpi_test() .OR. my_process_is_mpi_workroot()) THEN
          ALLOCATE(r_tmp(nproma,nlevs,nblks))
        ELSE
          ! Dimensions 1 and 2 of r_tmp must always be nproma and nlevs,
          ! otherwise exchange_data doesn't work!
          ALLOCATE(r_tmp(nproma,nlevs,1))
        ENDIF

        ! Gather data on root
        IF(my_process_is_mpi_seq()) THEN
          DO jk = 1, nlevs
            DO i = 1, n_points
              r_tmp(idx_no(i),jk,blk_no(i)) = r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i))
            ENDDO
          ENDDO
        ELSE
          CALL exchange_data(p_pat, RECV=r_tmp, SEND=r_ptr)
        ENDIF

        IF(my_process_is_mpi_test() .OR. my_process_is_mpi_workroot()) THEN

          ! De-block the array
          IF(use_sp_output) THEN
            ALLOCATE(r_out_sp(n_points, nlevs), STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
            DO jk = 1, nlevs
              r_out_sp(:,jk) = REAL(RESHAPE(r_tmp(:,jk,:), (/ n_points /)), sp)
            ENDDO
          ELSE
            ALLOCATE(r_out_dp(n_points, nlevs), STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
            DO jk = 1, nlevs
              r_out_dp(:,jk) = REAL(RESHAPE(r_tmp(:,jk,:), (/ n_points /)), dp)
            ENDDO
          ENDIF

          ! In the case of a test run: Compare results on worker PEs and test PE
          IF(p_test_run .AND. .NOT.use_async_name_list_io .AND. .NOT.use_sp_output) THEN
            ! Currently we don't do the check for REAL*4, we would need
            ! p_send/p_recv for this type
            IF(.NOT. my_process_is_mpi_test()) THEN
              ! Send to test PE
              CALL p_send(r_out_dp, process_mpi_all_test_id, 1)
            ELSE
              ! Receive result from parallel worker PEs and check for correctness
              ALLOCATE(r_out_recv(n_points,nlevs))
              CALL p_recv(r_out_recv, process_mpi_all_workroot_id, 1)
              IF(ANY(r_out_recv(:,:) /= r_out_dp(:,:))) &
                CALL finish(routine,'Sync error test PE/worker PEs for '//TRIM(info%name))
              DEALLOCATE(r_out_recv)
            ENDIF
          ENDIF

        ENDIF
        !
        ! write data
        !
        IF (my_process_is_stdio() .AND. .NOT. my_process_is_mpi_test()) THEN
          IF(use_sp_output) THEN
            CALL streamWriteVarF(of%cdiFileID, info%cdiVarID, r_out_sp, 0)
          ELSE
            CALL streamWriteVar(of%cdiFileID, info%cdiVarID, r_out_dp, 0)
          ENDIF
        ENDIF

        DEALLOCATE(r_tmp)
        IF(my_process_is_mpi_test() .OR. my_process_is_mpi_workroot()) THEN
          IF(use_sp_output) THEN
            DEALLOCATE(r_out_sp, STAT=ierrstat)
          ELSE
            DEALLOCATE(r_out_dp, STAT=ierrstat)
          END IF
          IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
        ENDIF

      ELSE

        ! Asynchronous I/O is used, just copy the OWN data points to the memory window

        DO jk = 1, nlevs
          IF(use_sp_output) THEN
            DO i = 1, p_ri%n_own
              mem_ptr_sp(ioff+INT(i,i8)) = REAL(r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),sp)
            ENDDO
          ELSE
            DO i = 1, p_ri%n_own
              mem_ptr_dp(ioff+INT(i,i8)) = REAL(r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i)),dp)
            ENDDO
          END IF
          ioff = ioff + INT(p_ri%n_own,i8)
        END DO
      END IF

      ! Reset variable if laccu flag is set

      IF(info%laccu) THEN
        SELECT CASE (info%ndims)
        CASE (2)
          IF(ASSOCIATED(of%var_desc(iv)%r_ptr)) THEN
            of%var_desc(iv)%r_ptr(:,:,nindex,1,1)          = info%resetval%rval
          ELSE
            of%var_desc(iv)%tlev_ptr(tl)%p(:,:,nindex,1,1) = info%resetval%rval
          ENDIF
        CASE (3)
          IF(ASSOCIATED(of%var_desc(iv)%r_ptr)) THEN
            of%var_desc(iv)%r_ptr(:,:,:,nindex,1)          = info%resetval%rval
          ELSE
            of%var_desc(iv)%tlev_ptr(tl)%p(:,:,:,nindex,1) = info%resetval%rval
          ENDIF
        END SELECT
      ENDIF

      IF(info%ndims == 2) DEALLOCATE(r_ptr)

    ENDDO

#ifndef NOMPI
    ! In case of async IO: Done writing to memory window, unlock it
    IF(use_async_name_list_io .AND. .NOT.my_process_is_mpi_test()) &
      CALL MPI_Win_unlock(p_pe_work, mpi_win, mpierr)
#endif


  END SUBROUTINE write_name_list

  !------------------------------------------------------------------------------------------------
  !> Returns if it is time for the next output step
  !! Please note:
  !! This function returns .TRUE. whenever the next output time of any name list
  !! is less or equal (sim_time+i_sample*dtime/2._wp).
  !! It DOES NOT check if the step indicated at sim_time is an advection step.
  !! THIS CHECK MUST BE DONE BY THE CALLER !!!!

  FUNCTION istime4name_list_output(sim_time) RESULT(retval)

    REAL(wp), INTENT(IN)  :: sim_time            ! simulation time [s]
    LOGICAL               :: retval

    TYPE(t_output_name_list), POINTER :: p_onl

    retval = .FALSE.

    ! Go through all output name lists and check if next_output_time is reached

    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT
      IF (retval) EXIT
      retval = is_output_nml_active(p_onl, sim_time, dtime, i_sample)
      p_onl => p_onl%next
    ENDDO

  END FUNCTION istime4name_list_output


  !------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------
  ! The following routines are only needed for asynchronous IO
  !------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------

#ifdef NOMPI
  ! Just define the entry point of name_list_io_main_proc, it will never be called

  SUBROUTINE name_list_io_main_proc
  END SUBROUTINE name_list_io_main_proc

#else

  !-------------------------------------------------------------------------------------------------
  !>
  !! Main routine for I/O PEs.
  !! Please note that this routine never returns.

  SUBROUTINE name_list_io_main_proc

    LOGICAL :: done, last_step
    TYPE(t_datetime) :: datetime
    REAL(wp) :: sim_time
    INTEGER :: jg


    print '(a,i0)','============================================ Hello from I/O PE ',p_pe

    ! If ldump_states or ldump_dd is set, the compute PEs will exit after dumping,
    ! there is nothing to do at all for I/O PEs

    IF(ldump_states .OR. ldump_dd) THEN
      CALL p_stop
      STOP
    ENDIF

    ! setup of meteogram output
    DO jg =1,n_dom
      IF (meteogram_output_config(jg)%lenabled) THEN
        CALL meteogram_init(meteogram_output_config(jg), jg)
      END IF
    END DO

    ! Initialize name list output, this is a collective call for all PEs
    CALL init_name_list_output

    ! Tell the compute PEs that we are ready to work
    CALL async_io_send_ready_message

    ! write recent samples of meteogram output
    DO jg = 1, n_dom
      IF (meteogram_output_config(jg)%lenabled) THEN
        CALL meteogram_flush_file(jg)
      END IF
    END DO

    ! Enter I/O loop

    DO
      
      ! Wait for a message from the compute PEs to start
      CALL async_io_wait_for_start(done, datetime, sim_time, last_step)

      IF(done) EXIT ! leave loop, we are done

      ! perform I/O
      CALL write_name_list_output(datetime, sim_time, last_step)

      ! Inform compute PEs that we are done
      CALL async_io_send_ready_message

      ! write recent samples of meteogram output
      DO jg = 1, n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          CALL meteogram_flush_file(jg)
        END IF
      END DO

    ENDDO

    ! Finalization sequence:

    print '(a,i0,a)','============================================ I/O PE ',p_pe,' shutting down'
    CALL close_name_list_output

    ! finalize meteogram output
    DO jg = 1, n_dom
      IF (meteogram_output_config(jg)%lenabled) THEN
        CALL meteogram_finalize(jg)
      END IF
    END DO
    DO jg = 1, max_dom
      DEALLOCATE(meteogram_output_config(jg)%station_list)
    END DO

    ! Shut down MPI
    !
    CALL p_stop

    STOP

  END SUBROUTINE name_list_io_main_proc

  !------------------------------------------------------------------------------------------------
  !>
  !! Transfers reorder_info to IO PEs

  SUBROUTINE transfer_reorder_info(p_ri)

    TYPE(t_reorder_info), INTENT(INOUT) :: p_ri ! Result: reorder info

    ! There is nothing to do for the test PE:
    IF(my_process_is_mpi_test()) RETURN

    ! Transfer the global number of points, this is not yet known on IO PEs
    CALL p_bcast(p_ri%n_glb, bcast_root, p_comm_work_2_io)
    CALL p_bcast(p_ri%n_log, bcast_root, p_comm_work_2_io)

    IF(my_process_is_io()) THEN

      ! On IO PEs: n_own = 0, own_idx and own_blk are not allocated
      p_ri%n_own = 0

      ! pe_own/pe_off must be allocated for num_work_procs, not for p_n_work
      ALLOCATE(p_ri%pe_own(0:num_work_procs-1))
      ALLOCATE(p_ri%pe_off(0:num_work_procs-1))

      ALLOCATE(p_ri%reorder_index(p_ri%n_glb))
      IF (l_grid_info_from_file) THEN
        ALLOCATE(p_ri%log_dom_index(p_ri%n_glb))
      END IF
    ENDIF

    CALL p_bcast(p_ri%pe_own, bcast_root, p_comm_work_2_io)
    CALL p_bcast(p_ri%pe_off, bcast_root, p_comm_work_2_io)

    CALL p_bcast(p_ri%reorder_index, bcast_root, p_comm_work_2_io)
    IF (l_grid_info_from_file) THEN
      CALL p_bcast(p_ri%log_dom_index, bcast_root, p_comm_work_2_io)
    END IF

  END SUBROUTINE transfer_reorder_info

  !-------------------------------------------------------------------------------------------------
  !>
  !! Replicates data (mainly the variable lists) needed for async I/O on the I/O procs.
  !! ATTENTION: The data is not completely replicated, only as far as needed for I/O.
  !!
  !! This routine has to be called by all PEs (work and I/O)

  SUBROUTINE replicate_data_on_io_procs

    INTEGER :: ivct_len
    INTEGER :: info_size, iv, nv, nelems, n, list_info(4)
    INTEGER, ALLOCATABLE :: info_storage(:,:)

    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_metadata) :: info
    TYPE(t_var_list) :: p_var_list
    ! var_list_name should have at least the length of var_list names
    ! (although this doesn't matter as long as it is big enough for every name)
    CHARACTER(LEN=256) :: var_list_name
    INTEGER :: idom, ierrstat, dim_c, dim_e, dim_v

    CHARACTER(len=*), PARAMETER :: routine = 'mo_name_list_output/replicate_data_on_io_procs'

    ! There is nothing to do for the test PE:
    IF(my_process_is_mpi_test()) RETURN

    !-----------------------------------------------------------------------------------------------

    ! Replicate vertical coordinate table

    IF(.NOT.my_process_is_io()) ivct_len = SIZE(vct)
    CALL p_bcast(ivct_len, bcast_root, p_comm_work_2_io)

    IF(my_process_is_io()) ALLOCATE(vct(ivct_len))
    CALL p_bcast(vct, bcast_root, p_comm_work_2_io)

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

    !-----------------------------------------------------------------------------------------------
    ! Replicate coordinates of cells/edges/vertices:
    IF (.NOT. l_grid_info_from_file) THEN

      ! Go over all output domains
      DO idom = 1, n_dom_out
        CALL p_bcast(patch_info(idom)%nblks_glb_c, bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%nblks_glb_e, bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%nblks_glb_v, bcast_root, p_comm_work_2_io)

        dim_c =   global_cell_type
        dim_e =                  4
        dim_v = 9-global_cell_type

        IF(my_process_is_io()) THEN    
          ALLOCATE(patch_info(idom)%grid_c%lon(nproma*patch_info(idom)%nblks_glb_c),  &
            &      patch_info(idom)%grid_c%lat(nproma*patch_info(idom)%nblks_glb_c),  &
            &      patch_info(idom)%grid_c%lonv(nproma*patch_info(idom)%nblks_glb_c, dim_c), &
            &      patch_info(idom)%grid_c%latv(nproma*patch_info(idom)%nblks_glb_c, dim_c), &
                                !
            &      patch_info(idom)%grid_e%lon(nproma*patch_info(idom)%nblks_glb_e), &
            &      patch_info(idom)%grid_e%lat(nproma*patch_info(idom)%nblks_glb_e), &
            &      patch_info(idom)%grid_e%lonv(nproma*patch_info(idom)%nblks_glb_e, dim_e), &
            &      patch_info(idom)%grid_e%latv(nproma*patch_info(idom)%nblks_glb_e, dim_e), &
                                !
            &      patch_info(idom)%grid_v%lon(nproma*patch_info(idom)%nblks_glb_v), &
            &      patch_info(idom)%grid_v%lat(nproma*patch_info(idom)%nblks_glb_v), &
            &      patch_info(idom)%grid_v%lonv(nproma*patch_info(idom)%nblks_glb_v, dim_v), &
            &      patch_info(idom)%grid_v%latv(nproma*patch_info(idom)%nblks_glb_v, dim_v), &
            &      STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
        END IF

        ! cells
        CALL p_bcast(patch_info(idom)%grid_c%lon,  bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%grid_c%lat,  bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%grid_c%lonv, bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%grid_c%latv, bcast_root, p_comm_work_2_io)
        ! edges
        CALL p_bcast(patch_info(idom)%grid_e%lon,  bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%grid_e%lat,  bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%grid_e%lonv, bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%grid_e%latv, bcast_root, p_comm_work_2_io)
        ! vertices
        CALL p_bcast(patch_info(idom)%grid_v%lon,  bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%grid_v%lat,  bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%grid_v%lonv, bcast_root, p_comm_work_2_io)
        CALL p_bcast(patch_info(idom)%grid_v%latv, bcast_root, p_comm_work_2_io)
      END DO
    END IF ! IF (.NOT. l_grid_info_from_file)

  END SUBROUTINE replicate_data_on_io_procs

  !------------------------------------------------------------------------------------------------
  !>
  !! Initializes the memory window for asynchronous IO

  SUBROUTINE init_memory_window

    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_INFO_NULL

    INTEGER :: jp, i, iv, nlevs
    INTEGER :: nbytes_real, mpierr, rma_cache_hint
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_size, mem_bytes
    TYPE (t_var_metadata), POINTER :: info
#ifdef USE_CRAY_POINTER
    INTEGER (KIND=MPI_ADDRESS_KIND) :: iptr
    REAL(sp) :: tmp_sp
    POINTER(tmp_ptr_sp,tmp_sp(*))
    REAL(dp) :: tmp_dp
    POINTER(tmp_ptr_dp,tmp_dp(*))
#else
    TYPE(c_ptr) :: c_mem_ptr
#endif

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/init_async_name_list_output'
    INTEGER :: i_log_dom, n_own, lonlat_id

    ! There is nothing to do for the test PE:
    IF(my_process_is_mpi_test()) RETURN


    ! Get size and offset of the data for every output file

    mem_size = 0_i8

    ! Go over all output files
    DO i = 1, SIZE(output_file)

      output_file(i)%my_mem_win_off = mem_size

      ! Go over all name list variables for this output file

      DO iv = 1, output_file(i)%num_vars

        info => output_file(i)%var_desc(iv)%info
        jp = output_file(i)%phys_patch_id

        IF(info%ndims == 2) THEN
          nlevs = 1
        ELSE
          nlevs = info%used_dimensions(2)
        ENDIF

        SELECT CASE (info%hgrid)
          CASE (GRID_UNSTRUCTURED_CELL)
            mem_size = mem_size + INT(nlevs*patch_info(jp)%cells%n_own,i8)
          CASE (GRID_UNSTRUCTURED_EDGE)
            mem_size = mem_size + INT(nlevs*patch_info(jp)%edges%n_own,i8)
          CASE (GRID_UNSTRUCTURED_VERT)
            mem_size = mem_size + INT(nlevs*patch_info(jp)%verts%n_own,i8)
          CASE (GRID_REGULAR_LONLAT)
            lonlat_id = info%hor_interp%lonlat_id
            i_log_dom = output_file(i)%log_patch_id
            n_own     = lonlat_info(lonlat_id, i_log_dom)%n_own
            mem_size  = mem_size + INT(nlevs*n_own,i8)
          CASE DEFAULT
            CALL finish(routine,'unknown grid type')
        END SELECT

      ENDDO

      ! Get the offset on all PEs
      ALLOCATE(output_file(i)%mem_win_off(0:num_work_procs-1))
      IF(.NOT.my_process_is_io()) THEN
        CALL MPI_Allgather(output_file(i)%my_mem_win_off, 1, p_int_i8, &
                           output_file(i)%mem_win_off, 1, p_int_i8,    &
                           p_comm_work, mpierr)
      ENDIF

      CALL p_bcast(output_file(i)%mem_win_off, bcast_root, p_comm_work_2_io)

    ENDDO

    ! mem_size is calculated as number of variables above, get number of bytes

    ! Get the amount of bytes per REAL*4 or REAL*8 variable (as used in MPI
    ! communication)
    IF(use_sp_output) THEN
      CALL MPI_Type_extent(p_real_sp, nbytes_real, mpierr)
    ELSE
      CALL MPI_Type_extent(p_real_dp, nbytes_real, mpierr)
    ENDIF

    ! For the IO PEs the amount of memory needed is 0 - allocate at least 1 word there:
    mem_bytes = MAX(mem_size,1_i8)*INT(nbytes_real,i8)

    ! allocate amount of memory needed with MPI_Alloc_mem

#ifdef USE_CRAY_POINTER
    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, iptr, mpierr)

    IF(use_sp_output) THEN
      tmp_ptr_sp = iptr
      CALL set_mem_ptr_sp(tmp_sp, INT(mem_size))
    ELSE
      tmp_ptr_dp = iptr
      CALL set_mem_ptr_dp(tmp_dp, INT(mem_size))
    ENDIF
#else
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

    NULLIFY(mem_ptr_sp)
    NULLIFY(mem_ptr_dp)

#ifdef __SX__
    IF(use_sp_output) THEN
      CALL C_F_POINTER(c_mem_ptr, mem_ptr_sp, (/ INT(mem_size) /) )
    ELSE
      CALL C_F_POINTER(c_mem_ptr, mem_ptr_dp, (/ INT(mem_size) /) )
    ENDIF
#else
    IF(use_sp_output) THEN
      CALL C_F_POINTER(c_mem_ptr, mem_ptr_sp, (/ mem_size /) )
    ELSE
      CALL C_F_POINTER(c_mem_ptr, mem_ptr_dp, (/ mem_size /) )
    ENDIF
#endif
#endif

    rma_cache_hint = MPI_INFO_NULL
#ifdef __xlC__
    ! IBM specific RMA hint, that we don't want window caching
    CALL MPI_Info_create(rma_cache_hint, mpierr);
    IF (mpierr /= 0) CALL finish(trim(routine), "MPI error!")
    CALL MPI_Info_set(rma_cache_hint, "IBM_win_cache","0", mpierr)
    IF (mpierr /= 0) CALL finish(trim(routine), "MPI error!")
#endif

    ! Create memory window for communication
    IF(use_sp_output) THEN
      mem_ptr_sp(:) = 0._sp
      CALL MPI_Win_create( mem_ptr_sp,mem_bytes,nbytes_real,MPI_INFO_NULL,&
        &                  p_comm_work_io,mpi_win,mpierr )
    ELSE
      mem_ptr_dp(:) = 0._dp
      CALL MPI_Win_create( mem_ptr_dp,mem_bytes,nbytes_real,MPI_INFO_NULL,&
        &                  p_comm_work_io,mpi_win,mpierr )
    ENDIF
    IF (mpierr /= 0) CALL finish(TRIM(routine), "MPI error!")

#ifdef __xlC__
    CALL MPI_Info_free(rma_cache_hint, mpierr);
    IF (mpierr /= 0) CALL finish(trim(routine), "MPI error!")
#endif

  END SUBROUTINE init_memory_window

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


  !------------------------------------------------------------------------------------------------
  !>
  !! Output routine on the IO PE

  SUBROUTINE io_proc_write_name_list(of, l_first_write)

    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_LOCK_SHARED, MPI_MODE_NOCHECK, MPI_Wtime

    TYPE (t_output_file), INTENT(IN), TARGET :: of
    LOGICAL, INTENT(IN) :: l_first_write

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/io_proc_write_name_list'

    INTEGER nval, nlev_max, iv, jk, i, nlevs, mpierr, nv_off, np, i_dom, &
      &     lonlat_id, i_log_dom, ierrstat
    INTEGER(KIND=MPI_ADDRESS_KIND) :: ioff(0:num_work_procs-1)
    INTEGER :: voff(0:num_work_procs-1)
    REAL(sp), ALLOCATABLE :: var1_sp(:), var2_sp(:), var3_sp(:,:)
    REAL(dp), ALLOCATABLE :: var1_dp(:), var2_dp(:), var3_dp(:,:)
    TYPE (t_var_metadata), POINTER :: info
    TYPE(t_reorder_info), POINTER  :: p_ri
    CHARACTER*10 ctime
    REAL(wp) :: t_get, t_write, t_copy, t_intp, t_0, mb_get, mb_wr

#if defined (__SX__) && !defined (NOMPI)
! It may be necessary that var1 is in global memory on NEC
! (Note: this is only allowed when we compile with MPI.)
!CDIR GM_ARRAY(var1)
#endif

    t_get   = 0._wp
    t_write = 0._wp
    t_copy  = 0._wp
    t_intp  = 0._wp
    mb_get  = 0._wp
    mb_wr   = 0._wp

    CALL date_and_time(TIME=ctime)
    print '(a,i0,a)','#################### I/O PE ',p_pe,' starting I/O at '//ctime

    ! Get maximum number of data points in a slice and allocate tmp variables

    i_dom = of%phys_patch_id
    nval = MAX(patch_info(i_dom)%cells%n_glb, &
               patch_info(i_dom)%edges%n_glb, &
               patch_info(i_dom)%verts%n_glb)
    ! take also the lon-lat grids into account
    DO iv = 1, of%num_vars
      info => of%var_desc(iv)%info
      IF (info%hgrid == GRID_REGULAR_LONLAT) THEN
        lonlat_id = info%hor_interp%lonlat_id
        i_log_dom = of%log_patch_id
        p_ri  => lonlat_info(lonlat_id, i_log_dom)
        nval = MAX(nval, p_ri%n_glb)
      END IF
    END DO

    nlev_max = 1
    DO iv = 1, of%num_vars
      info => of%var_desc(iv)%info
      IF(info%ndims == 3) nlev_max = MAX(nlev_max, info%used_dimensions(2))
    ENDDO

    IF(use_sp_output) THEN
      ALLOCATE(var1_sp(nval*nlev_max), var2_sp(nval), STAT=ierrstat)
    ELSE
      ALLOCATE(var1_dp(nval*nlev_max), var2_dp(nval), STAT=ierrstat)
    ENDIF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ioff(:) = of%mem_win_off(:)


    ! Go over all name list variables for this output file

    DO iv = 1, of%num_vars

      info => of%var_desc(iv)%info

      ! inspect time-constant variables only if we are writing the
      ! first step in this file:
      IF ((info%cdiTimeID == TIME_CONSTANT) .AND. .NOT. l_first_write) CYCLE

      IF(info%ndims == 2) THEN
        nlevs = 1
      ELSE
        nlevs = info%used_dimensions(2)
      ENDIF

      ! Get pointer to appropriate reorder_info

      SELECT CASE (info%hgrid)
        CASE (GRID_UNSTRUCTURED_CELL)
          p_ri => patch_info(of%phys_patch_id)%cells
        CASE (GRID_UNSTRUCTURED_EDGE)
          p_ri => patch_info(of%phys_patch_id)%edges
        CASE (GRID_UNSTRUCTURED_VERT)
          p_ri => patch_info(of%phys_patch_id)%verts
        CASE (GRID_REGULAR_LONLAT)
          lonlat_id = info%hor_interp%lonlat_id
          i_log_dom = of%log_patch_id
          p_ri  => lonlat_info(lonlat_id, i_log_dom)
        CASE DEFAULT
          CALL finish(routine,'unknown grid type')
      END SELECT

      ! Retrieve part of variable from every worker PE using MPI_Get

      nv_off = 0
      DO np = 0, num_work_procs-1

        IF(p_ri%pe_own(np) == 0) CYCLE

        nval = p_ri%pe_own(np)*nlevs ! Number of words to transfer

        t_0 = MPI_Wtime()
        CALL MPI_Win_lock(MPI_LOCK_SHARED, np, MPI_MODE_NOCHECK, mpi_win, mpierr)

        IF(use_sp_output) THEN
          CALL MPI_Get(var1_sp(nv_off+1), nval, p_real_sp, np, ioff(np), &
            &          nval, p_real_sp, mpi_win, mpierr)
        ELSE
          CALL MPI_Get(var1_dp(nv_off+1), nval, p_real_dp, np, ioff(np), &
            &          nval, p_real_dp, mpi_win, mpierr)
        ENDIF

        CALL MPI_Win_unlock(np, mpi_win, mpierr)
        t_get  = t_get  + MPI_Wtime()-t_0
        mb_get = mb_get + nval

        ! Update the offset in var1
        nv_off = nv_off + nval

        ! Update the offset in the memory window on compute PEs
        ioff(np) = ioff(np) + INT(nval,i8)

      ENDDO

      ! compute the total offset for each PE
      nv_off = 0
      DO np = 0, num_work_procs-1
        voff(np) = nv_off
        nval     = p_ri%pe_own(np)*nlevs
        nv_off   = nv_off + nval
      END DO

      ! var1 is stored in the order in which the variable was stored on compute PEs,
      ! get it back into the global storage order

      t_0 = MPI_Wtime()
      IF(use_sp_output) THEN
        ALLOCATE(var3_sp(p_ri%n_glb,nlevs), STAT=ierrstat) ! Must be allocated to exact size
      ELSE
        ALLOCATE(var3_dp(p_ri%n_glb,nlevs), STAT=ierrstat) ! Must be allocated to exact size
      ENDIF
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
 
      DO jk = 1, nlevs
        nv_off = 0
        DO np = 0, num_work_procs-1
          IF(use_sp_output) THEN
            var2_sp(nv_off+1:nv_off+p_ri%pe_own(np)) = var1_sp(voff(np)+1:voff(np)+p_ri%pe_own(np))
          ELSE
            var2_dp(nv_off+1:nv_off+p_ri%pe_own(np)) = var1_dp(voff(np)+1:voff(np)+p_ri%pe_own(np))
          ENDIF
          nv_off = nv_off+p_ri%pe_own(np)
          voff(np) = voff(np)+p_ri%pe_own(np)
        ENDDO
        IF(use_sp_output) THEN
          DO i = 1, p_ri%n_glb
            var3_sp(i,jk) = var2_sp(p_ri%reorder_index(i))
          ENDDO
        ELSE
          DO i = 1, p_ri%n_glb
            var3_dp(i,jk) = var2_dp(p_ri%reorder_index(i))
          ENDDO
        ENDIF
      ENDDO ! Loop over levels
      t_copy = t_copy + MPI_Wtime()-t_0

      t_0 = MPI_Wtime()
      IF(use_sp_output) THEN
        CALL streamWriteVarF(of%cdiFileID, info%cdiVarID, var3_sp, 0)
        mb_wr = mb_wr + SIZE(var3_sp)
      ELSE
        CALL streamWriteVar(of%cdiFileID, info%cdiVarID, var3_dp, 0)
        mb_wr = mb_wr + SIZE(var3_dp)
      ENDIF
      t_write = t_write + MPI_Wtime()-t_0

      IF(use_sp_output) THEN
        DEALLOCATE(var3_sp, STAT=ierrstat)
      ELSE
        DEALLOCATE(var3_dp, STAT=ierrstat)
      ENDIF
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    ENDDO ! Loop over output variables

    IF(use_sp_output) THEN
      DEALLOCATE(var1_sp, var2_sp, STAT=ierrstat)
    ELSE
      DEALLOCATE(var1_dp, var2_dp, STAT=ierrstat)
    ENDIF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    CALL date_and_time(TIME=ctime)
    print '(a,i0,a)','#################### I/O PE ',p_pe,' done at '//ctime
    ! Convert mb_get/mb_wr to MB
    mb_get = mb_get*8*1.d-6
    mb_wr  = mb_wr*4*1.d-6 ! 4 byte since dp output is implicitly converted to sp
    ! writing this message causes a runtime error on the NEC because formatted output to stdio/stderr is limited to 132 chars
#ifndef __SX__
    IF (msg_level >= 12) THEN
      WRITE (message_text,'(10(a,f10.3))') &
           & ' Got ',mb_get,' MB, time get: ',t_get,' s [',mb_get/t_get, &
           & ' MB/s], time write: ',t_write,' s [',mb_wr/t_write,        &
           & ' MB/s], times copy+intp: ',t_copy+t_intp,' s'
      CALL message('',message_text)
    ENDIF
#endif

    ! Convert mb_get/mb_wr to MB
    IF(use_sp_output) THEN
      mb_get = mb_get*4*1.d-6
    ELSE
      mb_get = mb_get*8*1.d-6
    ENDIF
    mb_wr = mb_wr*4*1.d-6 ! always 4 byte since dp output is implicitly converted to sp
    ! PRINT *,' Got ',mb_get,' MB, time get: ',t_get,' s [',mb_get/t_get,&
    !   ' MB/s], time write: ',t_write,' s [',mb_wr/t_write, &
    !   ' MB/s], times copy+intp: ',t_copy+t_intp,' s'

  END SUBROUTINE io_proc_write_name_list

  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  ! Flow control routines between compute and IO procs ...
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------
  ! ... called on IO procs:
  !-------------------------------------------------------------------------------------------------
  !>
  !! async_io_send_ready_message: Send a message to the compute PEs that the I/O is ready
  !! The counterpart on the compute side is compute_wait_for_async_io

  SUBROUTINE async_io_send_ready_message

    REAL(wp) :: msg

    ! make sure all are done
    CALL p_barrier(comm=p_comm_work)

    ! Simply send a message from I/O PE 0 to compute PE 0
    IF(p_pe_work == 0) THEN
      msg = REAL(msg_io_done, wp)
      CALL p_send(msg, p_work_pe0, 0)
    ENDIF

  END SUBROUTINE async_io_send_ready_message

  !-------------------------------------------------------------------------------------------------
  !>
  !! async_io_wait_for_start: Wait for a message from I/O PEs that we should start I/O or finish
  !! The counterpart on the compute side is compute_start_io/compute_shutdown_io

  SUBROUTINE async_io_wait_for_start(done, datetime, sim_time, last_step)

    LOGICAL, INTENT(OUT)          :: done ! flag if we should shut down
    TYPE(t_datetime), INTENT(OUT) :: datetime
    REAL(wp), INTENT(OUT)         :: sim_time
    LOGICAL, INTENT(OUT)          :: last_step

    REAL(wp) :: msg(9)

    ! Set output parameters to default values

    done      = .FALSE.
    datetime  = zero_datetime
    sim_time  = 0._wp
    last_step = .FALSE.

    ! Receive message that we may start I/O (or should finish)

    IF(p_pe_work == 0) CALL p_recv(msg, p_work_pe0, 0)
    CALL p_bcast(msg, 0, comm=p_comm_work)

    SELECT CASE(INT(msg(1)))

    CASE(msg_io_start)

      datetime%year   = INT(msg(2))
      datetime%month  = INT(msg(3))
      datetime%day    = INT(msg(4))
      datetime%hour   = INT(msg(5))
      datetime%minute = INT(msg(6))
      datetime%second = msg(7)

      sim_time = msg(8)

      IF(msg(9) /= 0._wp) last_step = .TRUE.

    CASE(msg_io_shutdown)
      done = .TRUE.

    CASE DEFAULT
      ! Anything else is an error
      CALL finish(modname,'I/O PE: Got illegal I/O tag')

    END SELECT

  END SUBROUTINE async_io_wait_for_start

  !-------------------------------------------------------------------------------------------------
  ! ... called on compute procs:
  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_wait_for_async_io: Wait for a message that the I/O is ready
  !! The counterpart on the I/O side is io_send_ready_message

  SUBROUTINE compute_wait_for_async_io

    REAL(wp) :: msg

    ! First compute PE receives message from I/O leader
    IF(p_pe_work==0) THEN
      CALL p_recv(msg, p_io_pe0, 0)
      ! Just for safety: Check if we got the correct tag
      IF(INT(msg) /= msg_io_done) CALL finish(modname,'Compute PE: Got illegal I/O tag')
    ENDIF

    ! Wait in barrier until message is here
    CALL p_barrier(comm=p_comm_work)

  END SUBROUTINE compute_wait_for_async_io

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_start_async_io: Send a message to I/O PEs that they should start I/O
  !! The counterpart on the I/O side is io_wait_for_start_message

  SUBROUTINE compute_start_async_io(datetime, sim_time, last_step)

    TYPE(t_datetime), INTENT(IN) :: datetime
    REAL(wp), INTENT(IN)         :: sim_time
    LOGICAL, INTENT(IN)          :: last_step

    REAL(wp) :: msg(9)

    CALL p_barrier(comm=p_comm_work) ! make sure all are here

    msg(1) = REAL(msg_io_start,    wp)
    msg(2) = REAL(datetime%year,   wp)
    msg(3) = REAL(datetime%month,  wp)
    msg(4) = REAL(datetime%day,    wp)
    msg(5) = REAL(datetime%hour,   wp)
    msg(6) = REAL(datetime%minute, wp)
    msg(7) = datetime%second
    msg(8) = sim_time
    msg(9) = MERGE(1._wp, 0._wp, last_step)

    IF(p_pe_work==0) CALL p_send(msg, p_io_pe0, 0)

  END SUBROUTINE compute_start_async_io

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_shutdown_async_io: Send a message to I/O PEs that they should shut down
  !! The counterpart on the I/O side is io_wait_for_start_message

  SUBROUTINE compute_shutdown_async_io

    REAL(wp) :: msg(9)

    CALL p_barrier(comm=p_comm_work) ! make sure all are here

    msg(1) = REAL(msg_io_shutdown, wp)
    msg(2:) = 0._wp

    IF(p_pe_work==0) CALL p_send(msg, p_io_pe0, 0)

  END SUBROUTINE compute_shutdown_async_io

  !-------------------------------------------------------------------------------------------------
#endif

END MODULE mo_name_list_output
