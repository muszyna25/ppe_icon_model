
! Define USE_CRAY_POINTER for platforms having problems with ISO_C_BINDING
! BUT understand CRAY pointers

!   #define USE_CRAY_POINTER

MODULE mo_name_list_output

  ! Please note: The spelling "name_list" (with underscore) is intented to make
  ! clear that this does not pertain to a FORTRAN namelist but rather
  ! to a list of names of output variables

#ifndef USE_CRAY_POINTER
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_intptr_t, c_f_pointer
#endif

  USE mo_kind,                  ONLY: wp, i8
  USE mo_impl_constants,        ONLY: max_phys_dom, ihs_ocean, zml_soil
  USE mo_grid_config,           ONLY: n_dom, n_phys_dom, global_cell_type
  USE mo_cdi_constants          ! We need all
  USE mo_io_units,              ONLY: filename_max, nnml, nnml_output
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_var_metadata,          ONLY: t_var_metadata
  USE mo_linked_list            ! we need all
  USE mo_var_list,              ONLY: t_var_list, nvar_lists, max_var_lists, var_lists, &
                                      new_var_list
  ! MPI Communication routines
  USE mo_mpi,                   ONLY: p_send, p_recv, p_bcast, p_barrier, p_stop
  ! MPI Communicators
  USE mo_mpi,                   ONLY: p_comm_work, p_comm_work_io, p_comm_work_2_io
  ! MPI Data types
  USE mo_mpi,                   ONLY: p_real_dp, p_int, p_int_i8
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
  USE mo_communication,         ONLY: blk_no, idx_no, idx_1d
  USE mo_model_domain,          ONLY: t_patch, p_patch, p_phys_patch
  USE mo_parallel_config,       ONLY: nproma, p_test_run
  USE mo_vertical_coord_table,  ONLY: vct
  USE mo_dynamics_config,       ONLY: iequations, nnow, nnow_rcf
  USE mo_io_config,             ONLY: out_expname, lwrite_pzlev
  USE mo_run_config,            ONLY: num_lev, num_levp1, dtime, ldump_states
  USE mo_nh_pzlev_config,       ONLY: nh_pzlev_config
  USE mo_lnd_nwp_config,        ONLY: nlev_snow
  USE mo_datetime,              ONLY: t_datetime
  USE mo_math_utilities,        ONLY: t_lon_lat_grid
  USE mo_intp_data_strc,        ONLY: t_int_state, t_lon_lat_intp, p_int_state
  USE mo_intp_state,            ONLY: allocate_int_state_lonlat_grid, &
                                      deallocate_int_state_lonlat
  USE mo_intp_rbf_coeffs,       ONLY: rbf_setup_interpol_lonlat_grid
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_lonlat_nl, rbf_interpol_lonlat_nl

  USE mo_ocean_nml,             ONLY: n_zlev
  USE mo_oce_state,             ONLY: set_zlev

  USE mo_io_vlist,              ONLY: addGlobAtts, addAtmAtts, addOceAtts
  USE mo_var_list_element,      ONLY: t_var_list_element
  USE mo_communication,         ONLY: exchange_data, t_comm_pattern, idx_no, blk_no
  USE mo_interpol_config,       ONLY: rbf_vec_dim_c, rbf_vec_dim_v, rbf_vec_dim_e, rbf_c2grad_dim
  USE mo_nonhydrostatic_config, ONLY: iadv_rcf



  IMPLICIT NONE

  PRIVATE

  ! Flag whether name_list output is active, i.e. at least one /output_nml/ has been read

  LOGICAL, PUBLIC :: name_list_output_active = .FALSE.

  ! Flag whether async name_list I/O is used, it is set in the main program:

  LOGICAL, PUBLIC :: use_async_name_list_io = .FALSE.

  PUBLIC :: read_name_list_output_namelists
  PUBLIC :: init_name_list_output
  PUBLIC :: write_name_list_output
  PUBLIC :: close_name_list_output
  PUBLIC :: istime4name_list_output
  PUBLIC :: name_list_io_main_proc

  ! The following parameter decides whether physical or logical patches are output
  ! and thus whether the domain number in output name lists pertains to physical
  ! or logical patches.

  LOGICAL, PARAMETER :: l_output_phys_patch = .TRUE.

  INTEGER, PARAMETER :: &
    max_var_ml = 400, & ! maximum number of output model-level variables
    max_var_pl = 400, & ! maximum number of pressure-level variables
    max_var_hl = 400, & ! maximum number of height-level variables
    max_bounds = 100, & ! maximum number of output_bounds
    max_levels = 100, & ! maximum number of pressure/height levels
    vname_len  =  32    ! variable name length in I/O namelists

  TYPE t_output_name_list

    INTEGER  :: filetype            ! One of CDI's FILETYPE_XXX constants
    CHARACTER(LEN=8) :: namespace   ! 'DWD' - DWD short names (or 'MPIM', 'CMIP', 'ECMWF')
    INTEGER  :: mode                ! 1 = forecast mode, 2 = climate mode
    INTEGER  :: dom(max_phys_dom)   ! domains for which this namelist is used, ending with -1
    INTEGER  :: output_time_unit    ! 1 = second, 2=minute, 3=hour, 4=day, 5=month, 6=year
    REAL(wp) :: output_bounds(3,max_bounds) ! post-processing times in units defined by output_time_unit: start, end, increment
    INTEGER  :: steps_per_file      ! Max number of output steps in one output file
    LOGICAL  :: include_last        ! Flag whether to include the last timestep in output
    CHARACTER(LEN=filename_max) :: output_filename   ! output filename prefix
    LOGICAL  :: lwrite_ready        ! Flag. TRUE if a "ready file" (sentinel file) should be written at the end of each output stage
    CHARACTER(LEN=filename_max) :: ready_directory        ! output directory for ready files
    CHARACTER(LEN=vname_len)  :: ml_varlist(max_var_ml)   ! name of model level fields (translation to model by namespace)
    CHARACTER(LEN=vname_len)  :: pl_varlist(max_var_pl)   ! name of pressure level fields (translation to model by namespace)
    REAL(wp) :: p_levels(max_levels)                      ! pressure levels [hPa]
    CHARACTER(LEN=vname_len)  :: hl_varlist(max_var_hl)   ! name of height level fields
    REAL(wp) :: h_levels(max_levels)                      ! height levels
    INTEGER  :: remap               ! interpolate horizontally, 0: none, 1: to regular lat-lon grid, 2: to Gaussian grids, (3:...)
    LOGICAL  :: remap_internal      ! do interpolations online in the model or external (including triggering)
    REAL(wp) :: reg_lon_def(3)      ! if remap=1: start, increment, end longitude in degrees
    REAL(wp) :: reg_lat_def(3)      ! if remap=1: start, increment, end latitude in degrees
    INTEGER  :: gauss_tgrid_def     ! if remap=2: triangular truncation (e.g.63 for T63) for which the Gauss grid should be used
    REAL(wp) :: north_pole(2)       ! definition of north pole for rotated lon-lat grids.

    ! Internal members, not read from input
    INTEGER  :: lon_dim             ! Number of points in lon direction
    INTEGER  :: lat_dim             ! Number of points in lat direction
    INTEGER  :: cur_bounds_triple   ! current output_bounds triple in use
    REAL(wp) :: next_output_time    ! next output time (in seconds simulation time)
    INTEGER  :: n_output_steps
    TYPE(t_output_name_list), POINTER :: next ! Pointer to next output_name_list

  END TYPE t_output_name_list

  ! Pointer to a linked list of output name lists:
  TYPE(t_output_name_list), POINTER :: first_output_name_list => NULL()

  !------------------------------------------------------------------------------------------------
  ! Max number of time levels:
  INTEGER, PARAMETER :: max_time_levels = 5

  ! Unfortunately, Fortran does not allow arrays of pointers, so we have to define an extra type
  TYPE t_rptr_5d
    REAL(wp), POINTER :: p(:,:,:,:,:)
  END TYPE

  TYPE t_var_desc
    REAL(wp), POINTER :: r_ptr(:,:,:,:,:) ! Pointer to time level independent data (or NULL)
    TYPE(t_rptr_5d) :: tlev_ptr(max_time_levels) ! Pointers to time level dependet data
    TYPE(t_var_metadata) :: info          ! Info structure for variable
  END TYPE

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
  ! Ocean
  INTEGER, PARAMETER, PUBLIC      :: ZA_depth               = 10
  INTEGER, PARAMETER, PUBLIC      :: ZA_depth_half          = 11
  INTEGER, PARAMETER, PUBLIC      :: ZA_generic_ice         = 12

  !------------------------------------------------------------------------------------------------
  TYPE t_output_file

    ! The following data must be set before opening the output file:
    CHARACTER(LEN=filename_max) :: filename_pref ! Prefix of output file name
    INTEGER                     :: output_type   ! CDI format
    INTEGER                     :: phys_patch_id ! ID of physical output patch
    INTEGER                     :: log_patch_id  ! ID of logical output patch
    INTEGER                     :: num_vars
    TYPE(t_var_desc), ALLOCATABLE :: var_desc(:)
    TYPE(t_output_name_list), POINTER :: name_list ! Pointer to corresponding output name list

    INTEGER                     :: remap         ! Copy of remap from associated namelist

    !----------------------------
    ! Used for lon/lat interpolation only
    TYPE (t_lon_lat_grid)       :: lonlat_grid
    TYPE(t_lon_lat_intp)        :: int_state_lonlat
    !----------------------------

    INTEGER                     :: io_proc_id    ! ID of process doing I/O on this file

    !----------------------------
    ! Used for async IO only
    INTEGER(i8)                 :: my_mem_win_off
    INTEGER(i8), ALLOCATABLE    :: mem_win_off(:)
    !----------------------------

    ! The following members are set during open
    CHARACTER(LEN=filename_max) :: filename      ! Actual name of output file
    INTEGER                     :: cdiFileId
    INTEGER                     :: cdiVlistId         ! cdi vlist handler
    INTEGER                     :: cdiCellGridID
    INTEGER                     :: cdiVertGridID
    INTEGER                     :: cdiEdgeGridID
    INTEGER                     :: cdiLonLatGridID
    INTEGER                     :: cdiZaxisID(12) ! All types of possible Zaxis ID's
    INTEGER                     :: cdiTaxisID
    INTEGER                     :: cdiTimeIndex

  END TYPE t_output_file
  !------------------------------------------------------------------------------------------------

  TYPE(t_output_file), ALLOCATABLE, TARGET :: output_file(:)

  !------------------------------------------------------------------------------------------------
  ! TYPE t_reorder_info describes how local cells/edges/verts
  ! have to be reordered to get the global array.
  ! Below, "points" refers to either cells, edges or verts.

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
    INTEGER, ALLOCATABLE :: log_dom_index(:)
                      ! Index where a point of the physical domains is in the logical domain
  END TYPE t_reorder_info

  ! TYPE t_patch_info contains the reordering info for cells, edges and verts
  TYPE t_patch_info
    TYPE(t_reorder_info) :: cells
    TYPE(t_reorder_info) :: edges
    TYPE(t_reorder_info) :: verts
    INTEGER :: log_patch_id
  END TYPE t_patch_info

  TYPE(t_patch_info), ALLOCATABLE, TARGET :: patch_info(:)

  ! Number of output domains. This depends on l_output_phys_patch and is either the number
  ! of physical or the number of logical domains.

  INTEGER :: n_dom_out

  !------------------------------------------------------------------------------------------------
  ! Currently, we use only 1 MPI window for all output files
  INTEGER mpi_win
  REAL(wp), POINTER :: mem_ptr(:) ! Pointer to memory window
  !------------------------------------------------------------------------------------------------
  ! Broadcast root for intercommunicator broadcasts form compute PEs to IO PEs using p_comm_work_2_io
  INTEGER :: bcast_root
  !------------------------------------------------------------------------------------------------

  ! Tags for communication between compute PEs and I/O PEs

  INTEGER, PARAMETER :: msg_io_start    = 12345
  INTEGER, PARAMETER :: msg_io_done     = 54321
  INTEGER, PARAMETER :: msg_io_shutdown = 99999

  ! TYPE t_datetime has no default constructor for setting all members to 0 or a defined value.
  ! Thus we declare a instance here which should have zeors everywhere since it is static.

  TYPE(t_datetime) :: zero_datetime

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

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/read_name_list_output_namelists'

    ! Local variables corresponding to members of output_name_list
    INTEGER  :: filetype
    CHARACTER(LEN=8) :: namespace
    INTEGER  :: mode
    INTEGER  :: dom(max_phys_dom)
    INTEGER  :: output_time_unit
    REAL(wp) :: output_bounds(3,max_bounds)
    INTEGER  :: steps_per_file
    LOGICAL  :: include_last
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

    ! The namelist containing all variables above
    NAMELIST /output_nml/ &
      filetype,           &
      namespace,          &
      mode,               &
      dom,                &
      output_time_unit,   &
      output_bounds,      &
      steps_per_file,     &
      include_last,       &
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
    CALL position_nml ('output_nml', status=istat)
    IF(istat /= POSITIONED) THEN
      ! There exist no output_name_lists
      first_output_name_list => NULL()
      CALL close_nml
      RETURN
    ENDIF

    ! As in COSMO, there may exist several output_nml namelists in the input file
    ! Loop until EOF is reached

    p_onl => NULL()
    name_list_output_active = .FALSE.

    DO

      ! Set all variables in output_nml to their default values

      filetype           = FILETYPE_NC2 ! NetCDF
      namespace          = ' '
      mode               = 1
      dom(:)             = -1
      output_time_unit   = 1
      output_bounds(:,:) = 0._wp
      steps_per_file     = 100
      include_last       = .TRUE.
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

      ! Read output_nml

      READ (nnml, output_nml, iostat=istat)
      IF(istat < 0) EXIT ! No more namelists
      IF(istat > 0) CALL finish(routine, 'Read error in namelist "output_nml"')

      ! Check input

      ! We need dtime for this check
      IF(dtime<=0) CALL finish(routine, 'dtime must be set before reading output namelists')

      ! Output bounds
      IF(output_bounds(1,1) < 0 .OR. &
         output_bounds(2,1) <= output_bounds(1,1) .OR. &
         output_bounds(3,1) <= dtime) THEN
        CALL finish(routine,'Illegal output_bounds(:,1)')
      ENDIF

      DO i = 2, max_bounds-1
        IF(output_bounds(3,i) <= 0) EXIT ! The last one
        IF(output_bounds(1,i) <= output_bounds(2,i-1)) &
          CALL finish(routine,'output_bounds not increasing')
        IF(output_bounds(2,i) <= output_bounds(1,i)) &
          CALL finish(routine,'output_bounds end <= start')
        IF(output_bounds(3,i) <  dtime) &
          CALL finish(routine,'output_bounds inc < dtime')
      ENDDO

      ! For safety, at least last bounds triple must always be 0
      output_bounds(:,i:) = 0._wp

      ! Allocate next output_name_list

      IF(.NOT.ASSOCIATED(first_output_name_list)) THEN
        ! Allocate first name_list
        ALLOCATE(first_output_name_list)
        p_onl => first_output_name_list
      ELSE
        ! This is not the first one, p_onl points to the last one which was created
        ALLOCATE(p_onl%next)
        p_onl => p_onl%next
      ENDIF

      ! Set next output_name_list from values read

      p_onl%filetype         = filetype
      p_onl%namespace        = namespace
      p_onl%mode             = mode
      p_onl%dom(:)           = dom(:)
      p_onl%output_time_unit = output_time_unit

      ! output_bounds is always in seconds - the question is what to do with months or years
      ! output_time_unit: 1 = second, 2=minute, 3=hour, 4=day, 5=month, 6=year
      SELECT CASE(output_time_unit)
        CASE(1); p_onl%output_bounds(:,:) = output_bounds(:,:)
        CASE(2); p_onl%output_bounds(:,:) = output_bounds(:,:)*60
        CASE(3); p_onl%output_bounds(:,:) = output_bounds(:,:)*3600
        CASE(4); p_onl%output_bounds(:,:) = output_bounds(:,:)*86400
        CASE(5); p_onl%output_bounds(:,:) = output_bounds(:,:)*86400*30  ! Not a real calender month
        CASE(6); p_onl%output_bounds(:,:) = output_bounds(:,:)*86400*365 ! Not a real calender year
        CASE DEFAULT
          CALL finish(routine,'Illegal output_time_unit')
      END SELECT

      p_onl%steps_per_file   = steps_per_file
      p_onl%include_last     = include_last
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
      p_onl%reg_lon_def(:)   = reg_lon_def(:)
      p_onl%reg_lat_def(:)   = reg_lat_def(:)
      IF(remap/=0 .AND. remap/=1) &
        CALL finish(routine,'Unsupported value for remap')
      IF(remap==1) THEN
        IF(reg_lon_def(2)<=0._wp) CALL finish(routine,'Illegal LON increment')
        IF(reg_lat_def(2)<=0._wp) CALL finish(routine,'Illegal LAT increment')
        IF(reg_lon_def(3)<=reg_lon_def(1)) CALL finish(routine,'end lon <= start lon')
        IF(reg_lat_def(3)<=reg_lat_def(1)) CALL finish(routine,'end lat <= start lat')
        p_onl%lon_dim = INT( (reg_lon_def(3)-reg_lon_def(1))/reg_lon_def(2) ) + 1
        p_onl%lat_dim = INT( (reg_lat_def(3)-reg_lat_def(1))/reg_lat_def(2) ) + 1
      ENDIF
      p_onl%gauss_tgrid_def  = gauss_tgrid_def
      p_onl%north_pole(:)    = north_pole(:)

      p_onl%cur_bounds_triple= 1
      p_onl%next_output_time = p_onl%output_bounds(1,1)
      p_onl%n_output_steps   = 0
      p_onl%next => NULL()

    ENDDO

    CALL close_nml

    IF(ASSOCIATED(first_output_name_list)) name_list_output_active = .TRUE.

  END SUBROUTINE read_name_list_output_namelists

  !------------------------------------------------------------------------------------------------

  SUBROUTINE init_name_list_output

#ifndef NOMPI
    USE mpi, ONLY: MPI_ROOT, MPI_PROC_NULL
#endif

    INTEGER :: i, j, nfiles, i_typ, i_dom, nvl, vl_list(max_var_lists)
    CHARACTER(LEN=2) :: lev_type
    TYPE (t_output_name_list), POINTER :: p_onl
    TYPE (t_output_file), POINTER :: p_of
    TYPE(t_list_element), POINTER :: element

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/init_name_list_output'

    CALL message('init_name_list_output','Start')

    ! Preliminary until this is done correctly during creating the var_lists:
    ! Set loutput to .TRUE. unless we know exactly that we do NOT want to output this list

    DO i = 1, nvar_lists

      var_lists(i)%p%loutput = .TRUE.

      ! Disable output for var_lists which must not be output,
      ! e.g. because of name clashes with other var_lists, e.g.:
      !IF(var_lists(i)%p%name(1:15) == 'ext_data_atm_td' ) var_lists(i)%p%loutput = .FALSE.
      !IF(var_lists(i)%p%name(1:16) == 'nh_state_metrics') var_lists(i)%p%loutput = .FALSE.

! For a list of all variables, enable the following!
IF(.FALSE.) THEN
      IF (my_process_is_stdio()) THEN
        PRINT '(3a, i2)','Var_list name: ',TRIM(var_lists(i)%p%name), &
                         ' Patch: ',var_lists(i)%p%patch_id
        element => var_lists(i)%p%first_list_element
        DO
          IF(.NOT. ASSOCIATED(element)) EXIT
          PRINT *,'    ',element%field%info%name,element%field%info%loutput
          element => element%next_list_element
        ENDDO
      ENDIF
ENDIF

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

    ! If async IO is used, replicate data (mainly the variable lists) on IO procs

#ifndef NOMPI
    IF(use_async_name_list_io) CALL replicate_data_on_io_procs
#endif

    ! Set the number of domains in output and the patch reorder information
    CALL set_patch_info

    ! Get the number of output files needed (by counting the domains per name list)

    p_onl => first_output_name_list
    nfiles = 0

    DO

      IF(.NOT.ASSOCIATED(p_onl)) EXIT

      ! If dom(:) was not specified in namelist input, it is set completely to -1.
      ! In this case all domains are wanted in the output, so set it here
      ! appropriatly - this cannot be done during reading of the namelists
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
    CALL message('init_name_list_output',message_text)

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
        i_dom = p_onl%dom(i)

        ! Loop over model/pressure/height levels

        DO i_typ = 1, 3

          ! Check if name_list has variables of corresponding type
          IF(i_typ == 1 .AND. p_onl%ml_varlist(1) == ' ') CYCLE
          IF(i_typ == 2 .AND. p_onl%pl_varlist(1) == ' ') CYCLE
          IF(i_typ == 3 .AND. p_onl%hl_varlist(1) == ' ') CYCLE

          nfiles = nfiles+1
          p_of => output_file(nfiles)

          SELECT CASE(i_typ)
            CASE(1); lev_type = 'ML'
            CASE(2); lev_type = 'PL'
            CASE(3); lev_type = 'HL'
          END SELECT

          ! Set prefix of output_file name

          WRITE(p_of%filename_pref,'(a,"_DOM",i2.2,"_",a)') &
            TRIM(p_onl%output_filename),i_dom,lev_type

          p_of%phys_patch_id = i_dom
          p_of%log_patch_id  = patch_info(i_dom)%log_patch_id
          p_of%output_type   = p_onl%filetype
          p_of%name_list     => p_onl
          p_of%remap         = p_onl%remap

          ! Select all var_lists which belong to current logical domain and i_typ

          nvl = 0
          DO j = 1, nvar_lists

            IF(.NOT. var_lists(j)%p%loutput) CYCLE
            ! patch_id in var_lists always corresponds to the LOGICAL domain
            IF(var_lists(j)%p%patch_id /= patch_info(i_dom)%log_patch_id) CYCLE

            IF(i_typ /= var_lists(j)%p%level_type) CYCLE

            nvl = nvl + 1
            vl_list(nvl) = j

          ENDDO

          SELECT CASE(i_typ)
            CASE(1)
              CALL add_varlist_to_output_file(p_of,vl_list(1:nvl),p_onl%ml_varlist)
            CASE(2)
              CALL add_varlist_to_output_file(p_of,vl_list(1:nvl),p_onl%pl_varlist)
            CASE(3)
              CALL add_varlist_to_output_file(p_of,vl_list(1:nvl),p_onl%hl_varlist)
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

    ! Set up lonlat interpolation for output files where this is enabled.
    ! This is done after setting io_proc_id so that we know where to
    ! store the lonlat interpolation data.

    DO i = 1, nfiles
      IF(output_file(i)%remap == 1) CALL set_output_file_lonlat_intp(output_file(i))
    ENDDO

    CALL message('init_name_list_output','Done')

  END SUBROUTINE init_name_list_output

  !------------------------------------------------------------------------------------------------

  SUBROUTINE add_varlist_to_output_file(p_of, vl_list, varlist)

    TYPE (t_output_file), INTENT(INOUT) :: p_of
    INTEGER, INTENT(IN) :: vl_list(:)
    CHARACTER(LEN=*), INTENT(IN) :: varlist(:)

    INTEGER :: ivar, i, iv, idx, tl
    LOGICAL :: found
    TYPE(t_list_element), POINTER :: element

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/add_varlist_to_output_file'

    ! Get the number of variables in varlist
    DO ivar = 1, SIZE(varlist)
      IF(varlist(ivar) == ' ') EXIT ! Last one reached
    ENDDO

    p_of%num_vars = ivar-1
    IF(p_of%num_vars == 0) CALL finish(routine,'num_vars == 0')

    ! Allocate array of variable descriptions

    ALLOCATE(p_of%var_desc(p_of%num_vars))

    ! Loop over all variables in varlist

    DO ivar = 1, p_of%num_vars

      ! Nullify pointers in p_of%var_desc
      p_of%var_desc(ivar)%r_ptr => NULL()
      DO i = 1, max_time_levels
        p_of%var_desc(ivar)%tlev_ptr(i)%p => NULL()
      ENDDO

      found = .FALSE.

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

          ! Do not inspect element if it is a container
          IF(element%field%info%lcontainer) CYCLE

          ! Check for matching name
          idx = INDEX(element%field%info%name,'.TL')
          IF(idx == 0) THEN
            IF(varlist(ivar) /= element%field%info%name) CYCLE
          ELSE
            IF(varlist(ivar) /= element%field%info%name(1:idx-1)) CYCLE
          ENDIF

          ! Found it, add it to the variable list of output file

          IF(idx == 0) THEN
            ! Not time level dependent
            IF(found) CALL finish(routine,'Duplicate var name: '//TRIM(varlist(ivar)))
            p_of%var_desc(ivar)%r_ptr => element%field%r_ptr
            p_of%var_desc(ivar)%info  = element%field%info
          ELSE
            IF(found) THEN
              ! We have already the info field, make some plausibilty checks:
              IF(ANY(p_of%var_desc(ivar)%info%used_dimensions(:) /=  &
                     element%field%info%used_dimensions(:))) &
                CALL finish(routine,'Dimension mismatch TL variable: '//TRIM(varlist(ivar)))
              ! There must not be a TL independent variable with the same name
              IF(ASSOCIATED(p_of%var_desc(ivar)%r_ptr)) &
                CALL finish(routine,'Duplicate var name: '//TRIM(varlist(ivar)))
              ! Maybe some more members of info should be tested ...
            ELSE
              ! Variable encountered the first time, set info field ...
              p_of%var_desc(ivar)%info = element%field%info
              ! ... and set name without .TL# suffix
              p_of%var_desc(ivar)%info%name = element%field%info%name(1:idx-1)
            ENDIF
            ! Get time level
            tl = ICHAR(element%field%info%name(idx+3:idx+3)) - ICHAR('0')
            IF(tl<=0 .OR. tl>max_time_levels) &
              CALL finish(routine, 'Illegal time level in '//TRIM(element%field%info%name))
            IF(ASSOCIATED(p_of%var_desc(ivar)%tlev_ptr(tl)%p)) &
              CALL finish(routine, 'Duplicate time level for '//TRIM(element%field%info%name))
            p_of%var_desc(ivar)%tlev_ptr(tl)%p => element%field%r_ptr
          ENDIF

          found = .TRUE.

        ENDDO

      ENDDO

      ! Check that at least one element with this name has been found

      IF(.NOT. found) &
        CALL finish(routine,'Output name list variable not found: '//TRIM(varlist(ivar)))

    ENDDO

  END SUBROUTINE add_varlist_to_output_file

  !------------------------------------------------------------------------------------------------

  SUBROUTINE set_output_file_lonlat_intp(p_of)

    TYPE (t_output_file), INTENT(INOUT) :: p_of

    REAL(wp) :: pi_180
    INTEGER  :: jl

    pi_180 = ATAN(1._wp)/45._wp

    jl = p_of%log_patch_id

    ! Set grid

    p_of%lonlat_grid%delta(1)        = p_of%name_list%reg_lon_def(2) * pi_180
    p_of%lonlat_grid%delta(2)        = p_of%name_list%reg_lat_def(2) * pi_180
    p_of%lonlat_grid%start_corner(1) = p_of%name_list%reg_lon_def(1) * pi_180
    p_of%lonlat_grid%start_corner(2) = p_of%name_list%reg_lat_def(1) * pi_180
    p_of%lonlat_grid%poleN(1:2)      = p_of%name_list%north_pole(1:2) * pi_180
    p_of%lonlat_grid%dimen(1)        = p_of%name_list%lon_dim
    p_of%lonlat_grid%dimen(2)        = p_of%name_list%lat_dim

    p_of%lonlat_grid%total_dim = p_of%lonlat_grid%dimen(1)*p_of%lonlat_grid%dimen(2)
    p_of%lonlat_grid%nblks     = (p_of%lonlat_grid%total_dim - 1)/nproma + 1
    p_of%lonlat_grid%npromz    = p_of%lonlat_grid%total_dim - (p_of%lonlat_grid%nblks-1)*nproma


    ! Allocate interpolation state - this has to be on all PEs
    ! although finally we need it only on one

    CALL allocate_int_state_lonlat_grid(p_of%lonlat_grid, p_of%int_state_lonlat)

    ! Calculate lonlat state

    IF(.NOT.my_process_is_io()) THEN
      CALL rbf_setup_interpol_lonlat_grid(p_of%lonlat_grid, p_patch(jl), &
                                          p_of%int_state_lonlat, p_int_state(jl))
    ENDIF

    ! In the case of async IO: transfer int_state_lonlat to IO procs

    IF(use_async_name_list_io .AND. .NOT.my_process_is_mpi_test()) THEN
      CALL p_bcast(p_of%int_state_lonlat%rbf_vec_coeff,    bcast_root, p_comm_work_2_io)
      CALL p_bcast(p_of%int_state_lonlat%rbf_c2grad_coeff, bcast_root, p_comm_work_2_io)
      CALL p_bcast(p_of%int_state_lonlat%rbf_vec_idx,      bcast_root, p_comm_work_2_io)
      CALL p_bcast(p_of%int_state_lonlat%rbf_vec_blk,      bcast_root, p_comm_work_2_io)
      CALL p_bcast(p_of%int_state_lonlat%rbf_vec_stencil,  bcast_root, p_comm_work_2_io)
      CALL p_bcast(p_of%int_state_lonlat%rbf_c2grad_idx,   bcast_root, p_comm_work_2_io)
      CALL p_bcast(p_of%int_state_lonlat%rbf_c2grad_blk,   bcast_root, p_comm_work_2_io)
      CALL p_bcast(p_of%int_state_lonlat%rdist,            bcast_root, p_comm_work_2_io)
      CALL p_bcast(p_of%int_state_lonlat%tri_idx,          bcast_root, p_comm_work_2_io)
      ! nlocal_pts and owner are not needed and nlocal_pts has a different size on the IO PE,
      ! so don't transfer them!
    ENDIF

    ! Deallocate int_state_lonlat where it is not needed

    IF (.NOT.(p_of%io_proc_id == p_pe .OR. my_process_is_mpi_test()) ) THEN
      CALL deallocate_int_state_lonlat(p_of%int_state_lonlat)
    ENDIF

    ! Just for safety: Set this again - it might have been changed!
    p_of%lonlat_grid%total_dim = p_of%lonlat_grid%dimen(1)*p_of%lonlat_grid%dimen(2)
    p_of%lonlat_grid%nblks     = (p_of%lonlat_grid%total_dim - 1)/nproma + 1
    p_of%lonlat_grid%npromz    = p_of%lonlat_grid%total_dim - (p_of%lonlat_grid%nblks-1)*nproma

  END SUBROUTINE set_output_file_lonlat_intp

  !------------------------------------------------------------------------------------------------

  SUBROUTINE set_patch_info

    INTEGER :: jp, jl

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/set_patch_info'

    ! Set the number of output domains depending on l_output_phys_patch

    IF(l_output_phys_patch) THEN
      n_dom_out = n_phys_dom
    ELSE
      n_dom_out = n_dom
    ENDIF

    ALLOCATE(patch_info(n_dom_out))

    DO jp = 1, n_dom_out

      ! Set number of global cells/edges/verts and logical patch ID

      IF(l_output_phys_patch) THEN
        patch_info(jp)%log_patch_id = p_phys_patch(jp)%logical_id
      ELSE
        patch_info(jp)%log_patch_id = jp
      ENDIF

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
      ENDIF
#ifndef NOMPI
      IF(use_async_name_list_io) THEN
        ! Transfer reorder_info to IO PEs
        CALL transfer_reorder_info(patch_info(jp)%cells)
        CALL transfer_reorder_info(patch_info(jp)%edges)
        CALL transfer_reorder_info(patch_info(jp)%verts)
      ENDIF
#endif

    ENDDO

  END SUBROUTINE set_patch_info


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
    ALLOCATE(p_ri%log_dom_index(p_ri%n_glb))
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
        p_ri%log_dom_index(n) = i
      ENDIF
    ENDDO

    ! Safety check
    IF(n/=p_ri%n_glb) CALL finish(routine,'Reordering failed')

    DEALLOCATE(phys_owner_mask)
    DEALLOCATE(glbidx_own)
    DEALLOCATE(glbidx_glb)
    DEALLOCATE(reorder_index_log_dom)

  END SUBROUTINE set_reorder_info


  !------------------------------------------------------------------------------------------------

  SUBROUTINE open_output_file(of, jfile)

    TYPE(t_output_file), INTENT(INOUT) :: of
    INTEGER, INTENT(IN) :: jfile ! Number of file set to open

    INTEGER :: i, j, k, nlev, nlevp1, nplev, nzlev, nzlevp1, znlev_soil, astatus, iv, jv, i_dom
    INTEGER :: ll_dim(2)
    REAL(wp), ALLOCATABLE :: levels(:), levels_i(:), levels_m(:), p_lonlat(:)

    CHARACTER(LEN=16) :: extn

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/open_output_file'


    i_dom = of%phys_patch_id

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
      CALL message(routine,'GRIB2 support experimental')
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
      CALL message('',message_text)
      CALL finish (routine, 'open failed on '//TRIM(of%filename))
    ELSE
      CALL message (routine, 'opened '//TRIM(of%filename))
    END IF
    !
    ! The following sections add the file global properties collected in init_name_list_output
    !
    ! 1. create cdi vlist 
    !
    of%cdiVlistID = vlistCreate()
    !
    !    set cdi internal time index to 0 for writing time slices in netCDF
    !
    of%cdiTimeIndex = 0
    !
    ! 2. add global attributes for netCDF
    !
    ! 3. add horizontal grid descriptions

    IF(of%name_list%remap == 1) THEN

      ! Lon/Lat Interpolation requested

      of%cdiCellGridID = CDI_UNDEFID
      of%cdiEdgeGridID = CDI_UNDEFID
      of%cdiVertGridID = CDI_UNDEFID

      ll_dim(1) = of%name_list%lon_dim
      ll_dim(2) = of%name_list%lat_dim

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
        p_lonlat(k) = of%name_list%reg_lon_def(1) + k*of%name_list%reg_lon_def(2)
      ENDDO
      CALL gridDefXvals(of%cdiLonLatGridID, p_lonlat)
      DEALLOCATE(p_lonlat)

      ALLOCATE(p_lonlat(ll_dim(2)))
      DO k = 1, ll_dim(2)
        p_lonlat(k) = of%name_list%reg_lat_def(1) + k*of%name_list%reg_lat_def(2)
      ENDDO
      CALL gridDefYvals(of%cdiLonLatGridID, p_lonlat)
      DEALLOCATE(p_lonlat)

    ELSE
      ! Cells

      of%cdiCellGridID = gridCreate(GRID_UNSTRUCTURED, patch_info(i_dom)%cells%n_glb)
      CALL gridDefNvertex(of%cdiCellGridID, global_cell_type)
      !
      CALL gridDefXname(of%cdiCellGridID, 'clon')
      CALL gridDefXlongname(of%cdiCellGridID, 'center longitude')
      CALL gridDefXunits(of%cdiCellGridID, 'radians')
      !
      CALL gridDefYname(of%cdiCellGridID, 'clat')
      CALL gridDefYlongname(of%cdiCellGridID, 'center latitude')
      CALL gridDefYunits(of%cdiCellGridID, 'radians')

      ! Verts

      of%cdiVertGridID = gridCreate(GRID_UNSTRUCTURED, patch_info(i_dom)%verts%n_glb)
      CALL gridDefNvertex(of%cdiVertGridID, 9-global_cell_type)
      !
      CALL gridDefXname(of%cdiVertGridID, 'vlon')
      CALL gridDefXlongname(of%cdiVertGridID, 'vertex longitude')
      CALL gridDefXunits(of%cdiVertGridID, 'radians')
      !
      CALL gridDefYname(of%cdiVertGridID, 'vlat')
      CALL gridDefYlongname(of%cdiVertGridID, 'vertex latitude')
      CALL gridDefYunits(of%cdiVertGridID, 'radians')

      ! Edges

      of%cdiEdgeGridID = gridCreate(GRID_UNSTRUCTURED, patch_info(i_dom)%edges%n_glb)
      CALL gridDefNvertex(of%cdiEdgeGridID, 4)
      !
      CALL gridDefXname(of%cdiEdgeGridID, 'elon')
      CALL gridDefXlongname(of%cdiEdgeGridID, 'edge longitude')
      CALL gridDefXunits(of%cdiEdgeGridID, 'radians')
      !
      CALL gridDefYname(of%cdiEdgeGridID, 'elat')
      CALL gridDefYlongname(of%cdiEdgeGridID, 'edge latitude')
      CALL gridDefYunits(of%cdiEdgeGridID, 'radians')

      of%cdiLonLatGridID = CDI_UNDEFID
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

      ! Define axes for soil model

      of%cdiZaxisID(ZA_depth_below_land_p1) = &
        & zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+2)
      ALLOCATE(levels(znlev_soil+2))
      levels(1) = 0._wp
      DO k = 1, znlev_soil+1
        levels(k+1) = zml_soil(k)*100._wp
      END DO
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_depth_below_land_p1), levels)
      DEALLOCATE(levels)

      of%cdiZaxisID(ZA_depth_below_land) = &
        & zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+1)
      CALL zaxisDefLevels(of%cdiZaxisID(ZA_depth_below_land), zml_soil*100._wp)
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
      IF (lwrite_pzlev) THEN
        nplev = nh_pzlev_config(of%log_patch_id)%nplev
        of%cdiZaxisID(ZA_pressure) = zaxisCreate(ZAXIS_PRESSURE, nplev)
        ALLOCATE(levels(nplev))
        DO k = 1, nplev
          levels(k) = nh_pzlev_config(of%log_patch_id)%plevels(k)
        END DO
        CALL zaxisDefLevels(of%cdiZaxisID(ZA_pressure), levels)
        CALL zaxisDefVct(of%cdiZaxisID(ZA_pressure), nplev, levels)
        DEALLOCATE(levels)

        nzlev = nh_pzlev_config(of%log_patch_id)%nzlev
        of%cdiZaxisID(ZA_height)  = zaxisCreate(ZAXIS_HEIGHT, nzlev)
        ALLOCATE(levels(nzlev))
        DO k = 1, nzlev
          levels(k) = nh_pzlev_config(of%log_patch_id)%zlevels(k)
        END DO
        CALL zaxisDefLevels(of%cdiZaxisID(ZA_height), levels)
        CALL zaxisDefVct(of%cdiZaxisID(ZA_height), nzlev, levels)
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
    !
    of%cdiTaxisID = taxisCreate(TAXIS_ABSOLUTE)
    CALL vlistDefTaxis(of%cdiVlistID, of%cdiTaxisID)

    !
    ! 6. global attributes
    !
    CALL addGlobAtts(of%cdiVlistID,of%log_patch_id,astatus)
    IF (iequations/=ihs_ocean) THEN
      CALL addAtmAtts(of%cdiVlistID,of%log_patch_id,astatus)
    ELSE
      CALL addOceAtts(of%cdiVlistID,astatus)
    END IF

    !
    ! add variables
    !
    CALL add_variables_to_vlist(of)

    CALL streamDefVlist(of%cdiFileID, of%cdiVlistID)

  END SUBROUTINE open_output_file

  !------------------------------------------------------------------------------------------------
  !
  ! define variables and attributes
  !
  SUBROUTINE add_variables_to_vlist(of)

    TYPE (t_output_file), INTENT(IN), TARGET :: of

    TYPE (t_var_metadata), POINTER :: info
    !
    INTEGER :: iv, vlistID, varID, gridID, zaxisID, nlev, nlevp1, znlev_soil, i, nvars

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
      IF(of%name_list%remap == 1) THEN
        info%cdiGridID = of%cdiLonLatGridID
      ELSE
        SELECT CASE (info%hgrid)
        CASE(GRID_UNSTRUCTURED_CELL)
          info%cdiGridID = of%cdiCellGridID
        CASE(GRID_UNSTRUCTURED_VERT)
          info%cdiGridID = of%cdiVertGridID
        CASE(GRID_UNSTRUCTURED_EDGE)
          info%cdiGridID = of%cdiEdgeGridID
        CASE DEFAULT
          CALL finish(routine, 'GRID definition missing for '//TRIM(info%name))
        END SELECT
      ENDIF

      gridID = info%cdiGridID

      !
      ! set z axis ID
      !
      SELECT CASE (info%vgrid)

      CASE (ZAXIS_SURFACE)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_surface)

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
          PRINT *,'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
          CALL finish(routine,'Dimension mismatch for ZAXIS_DEPTH_BELOW_SEA')
        ENDIF

      CASE (ZAXIS_DEPTH_BELOW_LAND)
        IF (info%used_dimensions(2) == znlev_soil+1) THEN
          info%cdiZaxisID =  of%cdiZaxisID(ZA_depth_below_land)
        ELSE IF (info%used_dimensions(2) == znlev_soil+2) THEN
          info%cdiZaxisID =  of%cdiZaxisID(ZA_depth_below_land_p1)
        ELSE
          PRINT *,'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
          CALL finish(routine,'Dimension mismatch for ZAXIS_DEPTH_BELOW_LAND')
        ENDIF

      CASE (ZAXIS_HEIGHT)
        IF(info%name(1:8)=='lnd_prog' .AND. INDEX(info%name,'snow') /= 0) THEN
          ! This is a special case - use ZA_generic_snow[_p1]
          IF(info%used_dimensions(2) == nlev_snow) THEN
            info%cdiZaxisID =  of%cdiZaxisID(ZA_generic_snow)
          ELSE IF(info%used_dimensions(2) == nlev_snow+1) THEN
            info%cdiZaxisID =  of%cdiZaxisID(ZA_generic_snow_p1)
          ELSE
            PRINT *,'SNOW variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
            CALL finish(routine,'Dimension mismatch for ZAXIS_HEIGHT')
          ENDIF
        ELSE
          ! In all other cases, ZAXIS_HEIGHT seems to be equivalent to ZAXIS_HYBRID/ZAXIS_HYBRID_HALF
          ! TODO: Is there a difference to ZAXIS_HYBRID/ZAXIS_HYBRID_HALF ???
          IF (info%used_dimensions(2) == nlevp1) THEN
            info%cdiZaxisID =  of%cdiZaxisID(ZA_hybrid_half)
          ELSE IF (info%used_dimensions(2) == nlev) THEN
            info%cdiZaxisID =  of%cdiZaxisID(ZA_hybrid)
          ELSE
            PRINT *,'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
            CALL finish(routine,'Dimension mismatch for ZAXIS_HEIGHT')
          ENDIF
        ENDIF

      CASE (ZAXIS_PRESSURE)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_pressure)

      CASE (ZAXIS_ALTITUDE)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_height)

      CASE DEFAULT
        PRINT *,'Variable: ',TRIM(info%name),' ZAXIS: ',info%vgrid
        CALL finish(routine, 'ZAXIS definition missing for '//TRIM(info%name))

      END SELECT

      zaxisID = info%cdiZaxisID
      !
      IF(of%name_list%remap == 1 .AND. info%hgrid == GRID_UNSTRUCTURED_EDGE) THEN
        ! lat/lon interpolation and a vector variable => 2 variables are output
        nvars = 2
      ELSE
        nvars = 1
      ENDIF

      DO i = 1, nvars
        varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE)
        IF(i==1) THEN
          info%cdiVarID   = varID
        ELSE
          info%cdiVarID_2 = varID
        ENDIF
        !
        CALL vlistDefVarDatatype(vlistID, varID, DATATYPE_FLT32)
        IF(nvars==2) THEN
          IF(i==1) THEN
            CALL vlistDefVarName(vlistID, varID, TRIM(info%name)//'_LONLAT_X')
          ELSE
            CALL vlistDefVarName(vlistID, varID, TRIM(info%name)//'_LONLAT_Y')
          ENDIF
        ELSE
          IF(of%name_list%remap == 1) THEN
            CALL vlistDefVarName(vlistID, varID, TRIM(info%name)//'_LONLAT')
          ELSE
            CALL vlistDefVarName(vlistID, varID, TRIM(info%name))
          ENDIF
        ENDIF
        !
        IF (info%cf%long_name /= '') CALL vlistDefVarLongname(vlistID, varID, info%cf%long_name)
        IF (info%cf%units /= '') CALL vlistDefVarUnits(vlistID, varID, info%cf%units)

        ! Currently only real valued variables are allowed, so we can always use info%missval%rval
        IF (info%lmiss) CALL vlistDefVarMissval(vlistID, varID, info%missval%rval)
      ENDDO
      !
    ENDDO
    !
  END SUBROUTINE add_variables_to_vlist
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE close_name_list_output
    !
    ! Close all name_list files
    !
    INTEGER :: i

#ifndef NOMPI
    IF(use_async_name_list_io.AND..NOT.my_process_is_io().AND..NOT.my_process_is_mpi_test()) THEN
      CALL compute_wait_for_async_io()
      CALL compute_shutdown_async_io()
    ELSE
#endif
      DO i = 1, SIZE(output_file)
        CALL close_output_file(output_file(i))
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

    INTEGER :: i, j, fileID, vlistID


    fileID = of%cdiFileID
    IF(fileID /= CDI_UNDEFID) CALL streamClose(fileID)

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

    of%cdiFileID       = CDI_UNDEFID
    of%cdiVlistID      = CDI_UNDEFID
    of%cdiCellGridID   = CDI_UNDEFID
    of%cdiEdgeGridID   = CDI_UNDEFID
    of%cdiVertGridID   = CDI_UNDEFID
    of%cdiLonLatGridID = CDI_UNDEFID
    of%cdiTaxisID      = CDI_UNDEFID
    of%cdiZaxisID(:)   = CDI_UNDEFID

  END SUBROUTINE close_output_file

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

    INTEGER :: i, j, idate, itime, iret, n
    TYPE(t_output_name_list), POINTER :: p_onl
    CHARACTER(LEN=filename_max+100) :: text
    REAL(wp), PARAMETER :: eps = 1.d-10 ! Tolerance for checking output bounds

    CALL message('write_name_list_output','Start')

    ! If asynchronous I/O is enabled, the compute PEs have to make sure
    ! that the I/O PEs are ready with the last output step before
    ! writing data into the I/O memory window.
    ! This routine (write_name_list_output) is also called from the I/O PEs,
    ! but in this case the calling routine cares about the flow control.

#ifndef NOMPI
    IF(use_async_name_list_io) THEN
      IF(.NOT.my_process_is_io().AND..NOT.my_process_is_mpi_test()) &
        CALL compute_wait_for_async_io()
    ENDIF
#endif

    ! Check if files have to be (re)opened

    ! Go over all output files
    DO i = 1, SIZE(output_file)

      p_onl => output_file(i)%name_list

      ! Check if output is due for this file
      IF((p_onl%include_last .AND. last_step) .OR. &
          p_onl%next_output_time <= sim_time+iadv_rcf*dtime/2._wp) THEN

        IF(MOD(p_onl%n_output_steps,p_onl%steps_per_file) == 0) THEN
          IF (output_file(i)%io_proc_id == p_pe .OR. my_process_is_mpi_test()) THEN
            IF(p_onl%n_output_steps > 0) CALL close_output_file(output_file(i))
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

      p_onl => output_file(i)%name_list

      ! Check if output is due for this file
      IF((p_onl%include_last .AND. last_step) .OR. &
          p_onl%next_output_time <= sim_time+iadv_rcf*dtime/2._wp) THEN

        IF(.NOT.use_async_name_list_io) THEN
          WRITE(text,'(a,a,a,1pg15.9)') &
            'Output to ',TRIM(output_file(i)%filename),' at simulation time ',sim_time
          CALL message('mo_name_list_output',text)
        ENDIF

        IF (output_file(i)%io_proc_id == p_pe .OR. my_process_is_mpi_test()) THEN
          CALL taxisDefVdate(output_file(i)%cdiTaxisID, idate)
          CALL taxisDefVtime(output_file(i)%cdiTaxisID, itime)

          iret = streamDefTimestep(output_file(i)%cdiFileId, output_file(i)%cdiTimeIndex)
          output_file(i)%cdiTimeIndex = output_file(i)%cdiTimeIndex + 1
        ENDIF

        IF(my_process_is_io()) THEN
#ifndef NOMPI
          IF(output_file(i)%io_proc_id == p_pe) THEN
            WRITE(*,'(a,a,a,1pg15.9,a,i6)') &
              'Output to ',TRIM(output_file(i)%filename),' at simulation time ',sim_time, &
              ' by PE ',p_pe
            CALL io_proc_write_name_list(output_file(i))
          ENDIF
#endif
        ELSE
          CALL write_name_list(output_file(i))
        ENDIF
      ENDIF

    ENDDO

    ! Increase next_output_time for all output name lists which have been written.
    ! Please note that more than 1 output_file may point to the same name_list,
    ! so this must not be done in the loop above!

    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT

      IF((p_onl%include_last .AND. last_step) .OR. &
          p_onl%next_output_time <= sim_time+iadv_rcf*dtime/2._wp) THEN

        p_onl%n_output_steps = p_onl%n_output_steps + 1

      ENDIF

      ! Switch all name lists to next output time.
      ! This is done in a DO WHILE loop to catch the case
      ! where an output increment less than the time step
      ! or two output_bounds triples which are too close are specified.

      DO WHILE(p_onl%next_output_time <= sim_time+iadv_rcf*dtime/2._wp)
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

  END SUBROUTINE write_name_list_output

  !------------------------------------------------------------------------------------------------
  !
  ! Write a output name list
  !
  SUBROUTINE write_name_list(of)

#ifndef NOMPI
    USE mpi, ONLY: MPI_LOCK_EXCLUSIVE, MPI_MODE_NOCHECK
#endif

    TYPE (t_output_file), INTENT(INOUT), TARGET :: of

    INTEGER :: tl, i_dom, i_log_dom, i, iv, jk, n_points, nlevs, nblks, n, nindex
    INTEGER :: mpierr, ll_dim(2)
    INTEGER(i8) :: ioff
    TYPE (t_var_metadata), POINTER :: info
    TYPE(t_reorder_info), POINTER  :: p_ri
    REAL(wp), POINTER :: r_ptr(:,:,:)
    REAL(wp), ALLOCATABLE :: r_tmp(:,:,:), r_out(:,:), r_out_recv(:,:), r_ll(:,:,:)
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
        CALL finish(routine,'4d arrays not handled yet.')
      CASE (5)
        CALL finish(routine,'5d arrays not handled yet.')
      CASE DEFAULT 
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
          ALLOCATE(r_out(n_points, nlevs))
          DO jk = 1, nlevs
            r_out(:,jk) = RESHAPE(r_tmp(:,jk,:), (/ n_points /))
          ENDDO

          ! In the case of a test run: Compare results on worker PEs and test PE
          IF(p_test_run .AND. .NOT.use_async_name_list_io) THEN
            IF(.NOT. my_process_is_mpi_test()) THEN
              ! Send to test PE
              CALL p_send(r_out, process_mpi_all_test_id, 1)
            ELSE
              ! Receive result from parallel worker PEs and check for correctness
              ALLOCATE(r_out_recv(n_points,nlevs))
              CALL p_recv(r_out_recv, process_mpi_all_workroot_id, 1)
              IF(ANY(r_out_recv(:,:) /= r_out(:,:))) &
                CALL finish(routine,'Sync error test PE/worker PEs for '//TRIM(info%name))
              DEALLOCATE(r_out_recv)
            ENDIF
          ENDIF

        ENDIF
        !
        ! write data
        !
        IF (my_process_is_stdio()) THEN
          IF(of%name_list%remap == 1) THEN
            ! Lon/Lat Interpolation requested
            ll_dim(1) = of%name_list%lon_dim
            ll_dim(2) = of%name_list%lat_dim
            ALLOCATE(r_ll(ll_dim(1)*ll_dim(2),nlevs,2))
            CALL interpolate_lonlat(of, info, r_out, r_ll)
            CALL streamWriteVar(of%cdiFileID, info%cdiVarID, r_ll(:,:,1), 0)
            IF(info%hgrid == GRID_UNSTRUCTURED_EDGE) &
              CALL streamWriteVar(of%cdiFileID, info%cdiVarID_2, r_ll(:,:,2), 0)
            DEALLOCATE(r_ll)
          ELSE
            CALL streamWriteVar(of%cdiFileID, info%cdiVarID, r_out, 0)
          ENDIF
        ENDIF

        DEALLOCATE(r_tmp)
        IF(my_process_is_mpi_test() .OR. my_process_is_mpi_workroot()) DEALLOCATE(r_out)

      ELSE

        ! Asynchronous I/O is used, just copy the OWN data points to the memory window

        DO jk = 1, nlevs
          DO i = 1, p_ri%n_own
            mem_ptr(ioff+i) = r_ptr(p_ri%own_idx(i),jk,p_ri%own_blk(i))
          ENDDO
          ioff = ioff + p_ri%n_own
        ENDDO

      ENDIF

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
  !! is less or equal (sim_time+iadv_rcf*dtime/2._wp).
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
      IF(p_onl%next_output_time <= sim_time+iadv_rcf*dtime/2._wp) retval = .TRUE.
      p_onl => p_onl%next
    ENDDO

  END FUNCTION istime4name_list_output

  !------------------------------------------------------------------------------------------------

  SUBROUTINE interpolate_lonlat(of, info, r_grid, r_ll)

    TYPE (t_output_file), INTENT(IN), TARGET :: of
    TYPE (t_var_metadata), INTENT(IN) :: info
    REAL(wp), INTENT(IN) :: r_grid(:,:)    ! Variable on ICON grid
    REAL(wp), INTENT(INOUT) :: r_ll(:,:,:) ! Interploated variable


    INTEGER :: nblks_ll, npromz_ll, n_tot, nlev, nblks, jk, i, iidx, iblk
    REAL(wp), ALLOCATABLE :: r_full(:,:,:), r_out1(:,:,:), r_out2(:,:,:)
    TYPE(t_reorder_info), POINTER  :: p_ri

    nblks_ll  = of%lonlat_grid%nblks
    npromz_ll = of%lonlat_grid%npromz
    n_tot     = of%lonlat_grid%total_dim

    IF(info%ndims == 2) THEN
      nlev = 1
    ELSE
      nlev = info%used_dimensions(2)
    ENDIF

    ! Get pointer to appropriate reorder_info

    SELECT CASE (info%hgrid)
      CASE (GRID_UNSTRUCTURED_CELL)
        p_ri => patch_info(of%phys_patch_id)%cells
      CASE (GRID_UNSTRUCTURED_EDGE)
        p_ri => patch_info(of%phys_patch_id)%edges
      CASE (GRID_UNSTRUCTURED_VERT)
        p_ri => patch_info(of%phys_patch_id)%verts
      CASE DEFAULT
        CALL finish('interpolate_lonlat','unknown grid type')
    END SELECT

    ! The interpolation routines below need the complete logical patch,
    ! r_grid contains only the physical patch - we have to copy back
    ! the physical patch into the logical patch here (and re-block it)!

    nblks = (p_ri%n_log-1)/nproma + 1 ! Total number of blocks in logical domain
    ALLOCATE(r_full(nproma,nlev,nblks))
    r_full(:,:,:) = 0

    DO jk = 1, nlev
      DO i = 1, p_ri%n_glb
        iidx = idx_no(p_ri%log_dom_index(i))
        iblk = blk_no(p_ri%log_dom_index(i))
        r_full(iidx,jk,iblk) = r_grid(i,jk)
      ENDDO
    ENDDO

    ! Allocate output variables
    ALLOCATE(r_out1(nproma,nlev,nblks_ll))
    ALLOCATE(r_out2(nproma,nlev,nblks_ll)) ! Only needed for GRID_UNSTRUCTURED_EDGE

    SELECT CASE (info%hgrid)
    CASE (GRID_UNSTRUCTURED_CELL)

      ! for cell-based variables: interpolate gradients (finite differences) and reconstruct
      CALL rbf_interpol_lonlat_nl(r_full, of%int_state_lonlat, r_out1, nblks_ll, npromz_ll)

    CASE (GRID_UNSTRUCTURED_VERT)

      ! vertex-based variables (vorticity "VOR") are treated in a
      ! special way: They are first interpolated onto the cell centers
      ! and afterwards treated as scalar variables of type "GATHER_C"

      ! Note: The whole process is overly expensive! For the case that
      ! there are other vertex-based variables for output, this
      ! routine must be optimized!

      ! RJ: This doesn't work at the moment !!!
      CALL finish('interpolate_lonlat', 'Vertex based variables currently not supported!')

    CASE (GRID_UNSTRUCTURED_EDGE)

      ! for edge-based variables: simple interpolation
      CALL rbf_vec_interpol_lonlat_nl( r_full, of%int_state_lonlat, r_out1, r_out2, &
                                       nblks_ll, npromz_ll)
    CASE DEFAULT
      CALL finish('interpolate_lonlat','unknown grid type')
    END SELECT

    ! reshape lon-lat variable (nproma,nlev,nblks) -> (ntot,nlev)       
    DO jk = 1, nlev
      r_ll(1:n_tot,jk,1) = RESHAPE(r_out1(:,jk,:), (/n_tot/))
      IF(info%hgrid == GRID_UNSTRUCTURED_EDGE) &
        r_ll(1:n_tot,jk,2) = RESHAPE(r_out2(:,jk,:), (/n_tot/))
    ENDDO

    DEALLOCATE(r_out1, r_out2, r_full)

  END SUBROUTINE interpolate_lonlat

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
!!    INTEGER :: jg


    print '(a,i0)','============================================ Hello from I/O PE ',p_pe

    ! If ldump_states is set, the compute PEs will exit after dumping,
    ! there is nothing to do at all for I/O PEs

    IF(ldump_states) THEN
      CALL p_stop
      STOP
    ENDIF

!!    ! setup of meteogram output
!!    DO jg =1,n_dom
!!      IF (meteogram_output_config(jg)%lenabled) THEN
!!        CALL meteogram_init(meteogram_output_config(jg), jg)
!!      END IF
!!    END DO

    ! Initialize name list output, this is a collective call for all PEs

    CALL init_name_list_output

    ! Tell the compute PEs that we are ready to work

    CALL async_io_send_ready_message

    ! Enter I/O loop

    DO

      ! Wait for a message from the compute PEs to start
      CALL async_io_wait_for_start(done, datetime, sim_time, last_step)

      IF(done) EXIT ! leave loop, we are done

      ! perform I/O
      CALL write_name_list_output(datetime, sim_time, last_step)

!!      ! write recent samples of meteogram output
!!      DO jg = 1, n_dom
!!        IF (meteogram_output_config(jg)%lenabled) THEN
!!          CALL meteogram_flush_file(jg)
!!        END IF
!!      END DO

      ! Inform compute PEs that we are done
      CALL async_io_send_ready_message

    ENDDO

    ! Finalization sequence:

    print '(a,i0,a)','============================================ I/O PE ',p_pe,' shutting down'
    CALL close_name_list_output

!!    ! finalize meteogram output
!!    DO jg = 1, n_dom
!!      IF (meteogram_output_config(jg)%lenabled) THEN
!!        CALL meteogram_finalize(jg)
!!      END IF
!!    END DO
!!    DO jg = 1, max_dom
!!      DEALLOCATE(meteogram_output_config(jg)%station_list)
!!    END DO

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
      ALLOCATE(p_ri%log_dom_index(p_ri%n_glb))

    ENDIF

    CALL p_bcast(p_ri%pe_own, bcast_root, p_comm_work_2_io)
    CALL p_bcast(p_ri%pe_off, bcast_root, p_comm_work_2_io)

    CALL p_bcast(p_ri%reorder_index, bcast_root, p_comm_work_2_io)
    CALL p_bcast(p_ri%log_dom_index, bcast_root, p_comm_work_2_io)

  END SUBROUTINE transfer_reorder_info

  !-------------------------------------------------------------------------------------------------
  !>
  !! Replicates data (mainly the variable lists) needed for async I/O on the I/O procs.
  !! ATTENTION: The data is not completely replicated, only as far as needed for I/O.
  !!
  !! This routine has to be called by all PEs (work and I/O)

  SUBROUTINE replicate_data_on_io_procs

    USE mpi, ONLY: MPI_ROOT, MPI_PROC_NULL

    INTEGER :: ivct_len, jg
    INTEGER :: info_size, iv, nv, nelems, n, list_info(4)
    INTEGER, ALLOCATABLE :: info_storage(:,:)

    TYPE(t_list_element), POINTER :: element
    TYPE(t_var_metadata) :: info
    TYPE(t_var_list) :: p_var_list
    ! var_list_name should have at least the length of var_list names
    ! (although this doesn't matter as long as it is big enough for every name)
    CHARACTER(LEN=256) :: var_list_name

    CHARACTER(len=*), PARAMETER :: routine = 'mo_name_list_output/replicate_data_on_io_procs'

    ! There is nothing to do for the test PE:
    IF(my_process_is_mpi_test()) RETURN

    !-----------------------------------------------------------------------------------------------

    ! Replicate physical domain setup, only the number of domains and the logical ID is needed
    CALL p_bcast(n_phys_dom, bcast_root, p_comm_work_2_io)
    DO jg = 1, n_phys_dom
      CALL p_bcast(p_phys_patch(jg)%logical_id, bcast_root, p_comm_work_2_io)
    ENDDO

    ! Replicate vertical coordinate table

    IF(.NOT.my_process_is_io()) ivct_len = SIZE(vct)
    CALL p_bcast(ivct_len, bcast_root, p_comm_work_2_io)

    IF(my_process_is_io()) ALLOCATE(vct(ivct_len))
    CALL p_bcast(vct, bcast_root, p_comm_work_2_io)

    ! Replicate dimensions needed for lonlat interpolation

    CALL p_bcast(rbf_vec_dim_c,  bcast_root, p_comm_work_2_io)
    CALL p_bcast(rbf_vec_dim_e,  bcast_root, p_comm_work_2_io)
    CALL p_bcast(rbf_vec_dim_v,  bcast_root, p_comm_work_2_io)
    CALL p_bcast(rbf_c2grad_dim, bcast_root, p_comm_work_2_io)

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
        list_info(3) = var_lists(iv)%p%level_type
        list_info(4) = MERGE(1,0,var_lists(iv)%p%loutput)

      ENDIF

      ! Send basic info:

      CALL p_bcast(list_info, bcast_root, p_comm_work_2_io)

      IF(my_process_is_io()) THEN
        nelems = list_info(1)
        ! Create var list
        CALL new_var_list( p_var_list, var_list_name, patch_id=list_info(2), &
                           level_type=list_info(3), loutput=(list_info(4)/=0) )
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

  END SUBROUTINE replicate_data_on_io_procs

  !------------------------------------------------------------------------------------------------
  !>
  !! Initializes the memory window for asynchronous IO

  SUBROUTINE init_memory_window

    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_ROOT, MPI_PROC_NULL, MPI_INFO_NULL

    INTEGER :: jp, jl, n_cells, n_edges, n_verts, i, iv, nlevs
    INTEGER :: nbytes_real, mpierr
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_size, mem_bytes
    TYPE (t_var_metadata), POINTER :: info
#ifdef USE_CRAY_POINTER
    INTEGER (KIND=MPI_ADDRESS_KIND) :: iptr
    REAL(wp) :: tmp
    POINTER(tmp_ptr,tmp(*))
#else
    TYPE(c_ptr) :: c_mem_ptr
#endif

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/init_async_name_list_output'


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
            mem_size = mem_size + INT(nlevs*patch_info(jp)%cells%n_own)
          CASE (GRID_UNSTRUCTURED_EDGE)
            mem_size = mem_size + INT(nlevs*patch_info(jp)%edges%n_own)
          CASE (GRID_UNSTRUCTURED_VERT)
            mem_size = mem_size + INT(nlevs*patch_info(jp)%verts%n_own)
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

    ! Get the amount of bytes per default REAL variable (as used in MPI communication)
    CALL MPI_Type_extent(p_real_dp, nbytes_real, mpierr)

    ! For the IO PEs the amount of memory needed is 0 - allocate at least 1 word there:
    mem_bytes = MAX(mem_size,1_i8)*INT(nbytes_real,i8)

    ! allocate amount of memory needed with MPI_Alloc_mem

#ifdef USE_CRAY_POINTER
    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, iptr, mpierr)

    tmp_ptr = iptr
    CALL set_mem_ptr(tmp, INT(mem_size))
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

#ifdef __SX__
    CALL c_f_pointer(c_mem_ptr, mem_ptr, (/ INT(mem_size) /) )
#else
    CALL c_f_pointer(c_mem_ptr, mem_ptr, (/ mem_size /) )
#endif
#endif

    mem_ptr(:) = 0._wp

    ! Create memory window for communication

    CALL MPI_Win_create(mem_ptr,mem_bytes,nbytes_real,MPI_INFO_NULL,p_comm_work_io,mpi_win,mpierr)

  END SUBROUTINE init_memory_window

#ifdef USE_CRAY_POINTER
  !------------------------------------------------------------------------------------------------
  ! Helper routine for setting mem_ptr with the correct size information

  SUBROUTINE set_mem_ptr(arr, len)

    INTEGER len
    REAL(wp), TARGET :: arr(len)

    mem_ptr => arr

  END SUBROUTINE set_mem_ptr
#endif


  !------------------------------------------------------------------------------------------------
  !>
  !! Output routine on the IO PE

  SUBROUTINE io_proc_write_name_list(of)

    USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_LOCK_SHARED, MPI_MODE_NOCHECK

    TYPE (t_output_file), INTENT(IN), TARGET :: of

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/io_proc_write_name_list'

    INTEGER nval, nlev_max, iv, jk, i, nlevs, mpierr, nv_off, np, i_dom, ll_dim(2)
    INTEGER(KIND=MPI_ADDRESS_KIND) :: ioff(0:num_work_procs-1)
    INTEGER :: voff(0:num_work_procs-1)
    REAL(wp), ALLOCATABLE :: var1(:), var2(:), var3(:,:), r_ll(:,:,:)
    TYPE (t_var_metadata), POINTER :: info
    TYPE(t_reorder_info), POINTER  :: p_ri
    CHARACTER*10 ctime

#if defined (__SX__) && !defined (NOMPI)
! It may be necessary that var1 is in global memory on NEC
! (Note: this is only allowed when we compile with MPI.)
!CDIR GM_ARRAY(var1)
#endif

    CALL date_and_time(TIME=ctime)
    print '(a,i0,a)','#################### I/O PE ',p_pe,' starting I/O at '//ctime

    ! Get maximum number of data points in a slice and allocate tmp variables

    i_dom = of%phys_patch_id
    nval = MAX(patch_info(i_dom)%cells%n_glb, &
               patch_info(i_dom)%edges%n_glb, &
               patch_info(i_dom)%verts%n_glb)
    nlev_max = 1
    DO iv = 1, of%num_vars
      info => of%var_desc(iv)%info
      IF(info%ndims == 3) nlev_max = MAX(nlev_max, info%used_dimensions(2))
    ENDDO

    ALLOCATE(var1(nval*nlev_max), var2(nval))

    ioff(:) = of%mem_win_off(:)


    ! Go over all name list variables for this output file

    DO iv = 1, of%num_vars

      info => of%var_desc(iv)%info

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
        CASE DEFAULT
          CALL finish(routine,'unknown grid type')
      END SELECT

      ! Retrieve part of variable from every worker PE using MPI_Get

      nv_off = 0
      DO np = 0, num_work_procs-1

        voff(np) = nv_off

        IF(p_ri%pe_own(np) == 0) CYCLE

        nval = p_ri%pe_own(np)*nlevs ! Number of words to transfer

        CALL MPI_Win_lock(MPI_LOCK_SHARED, np, MPI_MODE_NOCHECK, mpi_win, mpierr)

        CALL MPI_Get(var1(nv_off+1), nval, p_real_dp, np, ioff(np), &
         &           nval, p_real_dp, mpi_win, mpierr)

        CALL MPI_Win_unlock(np, mpi_win, mpierr)

        ! Update the offset in var1
        nv_off = nv_off + nval

        ! Update the offset in the memory window on compute PEs
        ioff(np) = ioff(np) + INT(nval,i8)

      ENDDO

      ! var1 is stored in the order in which the variable was stored on compute PEs,
      ! get it back into the global storage order

      ALLOCATE(var3(p_ri%n_glb,nlevs)) ! Must be allocated to exact size

      DO jk = 1, nlevs
        nv_off = 0
        DO np = 0, num_work_procs-1
          var2(nv_off+1:nv_off+p_ri%pe_own(np)) = var1(voff(np)+1:voff(np)+p_ri%pe_own(np))
          nv_off = nv_off+p_ri%pe_own(np)
          voff(np) = voff(np)+p_ri%pe_own(np)
        ENDDO
        DO i = 1, p_ri%n_glb
          var3(i,jk) = var2(p_ri%reorder_index(i))
        ENDDO
      ENDDO ! Loop over levels

      IF(of%name_list%remap == 1) THEN
        ! Lon/Lat Interpolation requested
        ll_dim(1) = of%name_list%lon_dim
        ll_dim(2) = of%name_list%lat_dim
        ALLOCATE(r_ll(ll_dim(1)*ll_dim(2),nlevs,2))
        CALL interpolate_lonlat(of, info, var3, r_ll)
        CALL streamWriteVar(of%cdiFileID, info%cdiVarID, r_ll(:,:,1), 0)
        IF(info%hgrid == GRID_UNSTRUCTURED_EDGE) &
          CALL streamWriteVar(of%cdiFileID, info%cdiVarID_2, r_ll(:,:,2), 0)
        DEALLOCATE(r_ll)
      ELSE
        CALL streamWriteVar(of%cdiFileID, info%cdiVarID, var3, 0)
      ENDIF

      DEALLOCATE(var3)

    ENDDO ! Loop over output variables


    DEALLOCATE(var1, var2)

    CALL date_and_time(TIME=ctime)
    print '(a,i0,a)','#################### I/O PE ',p_pe,' done at '//ctime

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
