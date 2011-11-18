MODULE mo_name_list_output

  ! Please note: The spelling "name_list" (with underscore) is intented to make
  ! clear that this does not pertain to a FORTRAN namelist but rather
  ! to a list of names of output variables

  USE mo_kind,                  ONLY: wp
  USE mo_impl_constants,        ONLY: max_phys_dom, ihs_ocean, zml_soil
  USE mo_grid_config,           ONLY: n_dom, n_phys_dom, global_cell_type
  USE mo_cdi_constants          ! We need all
  USE mo_io_units,              ONLY: filename_max, nnml, nnml_output
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_var_metadata,          ONLY: t_var_metadata
  USE mo_linked_list            ! we need all
  USE mo_var_list,              ONLY: t_var_list, nvar_lists, var_lists
  USE mo_mpi,                   ONLY: p_send, p_recv,                                     &
                                      my_process_is_stdio, my_process_is_mpi_test,        &
                                      my_process_is_mpi_workroot, my_process_is_mpi_seq,  &
                                      process_mpi_all_test_id, process_mpi_all_workroot_id
  USE mo_model_domain,          ONLY: t_patch, p_patch, p_phys_patch
  USE mo_parallel_config,       ONLY: nproma, p_test_run
  USE mo_vertical_coord_table,  ONLY: vct
  USE mo_dynamics_config,       ONLY: iequations, nnow
  USE mo_io_config,             ONLY: out_expname, lwrite_pzlev
  USE mo_run_config,            ONLY: num_lev, num_levp1, dtime
  USE mo_nh_pzlev_config,       ONLY: nh_pzlev_config
  USE mo_lnd_nwp_config,        ONLY: nlev_snow
  USE mo_datetime,              ONLY: t_datetime

  USE mo_ocean_nml,             ONLY: n_zlev
  USE mo_oce_state,             ONLY: set_zlev

  USE mo_io_vlist,              ONLY: addGlobAtts, addAtmAtts, addOceAtts
  USE mo_var_list_element,      ONLY: t_var_list_element
  USE mo_communication,         ONLY: exchange_data, t_comm_pattern, idx_no, blk_no


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_name_list_output_namelists
  PUBLIC :: init_name_list_output
  PUBLIC :: write_name_list_output
  PUBLIC :: close_name_list_output
  PUBLIC :: istime4name_list_output

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
    LOGICAL  :: include_start       ! include initial time step
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
    INTEGER  :: cur_bounds_triple   ! current output_bounds triple in use
    REAL(wp) :: next_output_time    ! next output time (in seconds simulation time)
    INTEGER  :: n_output_steps
    TYPE(t_output_name_list), POINTER :: next ! Pointer to next output_name_list

  END TYPE t_output_name_list

  TYPE(t_output_name_list), POINTER :: first_output_name_list => NULL()

  !------------------------------------------------------------------------------------------------
  ! The following structure is like t_list_element (used in var_list's) but contains
  ! the additinal members time_level and info_1st for the handling of different time levels
  TYPE t_olist_element
    TYPE(t_var_list_element)       :: field
    INTEGER                        :: time_level
    TYPE(t_var_metadata), POINTER  :: info_1st
    TYPE(t_olist_element), POINTER :: next_list_element
  END TYPE t_olist_element

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
    INTEGER                     :: patch_id      ! ID of output patch (physical if l_output_phys_patch is set!)
    INTEGER                     :: output_type   ! CDI format
    TYPE(t_output_name_list), POINTER :: name_list ! Pointer to corresponding output name list
    TYPE(t_olist_element), POINTER    :: first_list_element ! Linked list of items to be output

    ! The following members are set during open
    CHARACTER(LEN=filename_max) :: filename      ! Actual name of output file
    INTEGER                     :: cdiFileId
    INTEGER                     :: cdiVlistId         ! cdi vlist handler
    INTEGER                     :: cdiCellGridID
    INTEGER                     :: cdiVertGridID
    INTEGER                     :: cdiEdgeGridID
    INTEGER                     :: cdiZaxisID(12) ! All types of possible Zaxis ID's
    INTEGER                     :: cdiTaxisID
    INTEGER                     :: cdiTimeIndex

  END TYPE t_output_file
  !------------------------------------------------------------------------------------------------

  TYPE(t_output_file), ALLOCATABLE, TARGET :: output_file(:)

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
    LOGICAL  :: include_start
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
      include_start,      &
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

    DO

      ! Set all variables in output_nml to their default values

      filetype           = FILETYPE_NC2 ! NetCDF
      namespace          = ' '
      mode               = 1
      dom(:)             = -1
      output_time_unit   = 1
      output_bounds(:,:) = 0._wp
      steps_per_file     = 100
      include_start      = .FALSE.
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
      north_pole(:)      = 0._wp

      ! Read output_nml

      READ (nnml, output_nml, iostat=istat)
      IF(istat /= 0) EXIT ! No more namelists

      ! Check input

      ! We need dtime for this check
      IF(dtime<=0) CALL finish('read_name_list_output_namelists', &
                               'dtime must be set before reading output namelists')

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
      p_onl%include_start    = include_start
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
      p_onl%gauss_tgrid_def  = gauss_tgrid_def
      p_onl%north_pole(:)    = north_pole(:)

      p_onl%cur_bounds_triple= 1
      p_onl%next_output_time = p_onl%output_bounds(1,1)
      p_onl%n_output_steps   = 0
      p_onl%next => NULL()

    ENDDO

    CALL close_nml

  END SUBROUTINE read_name_list_output_namelists

  !------------------------------------------------------------------------------------------------

  SUBROUTINE init_name_list_output

    INTEGER :: i, j, nfiles, i_typ, i_dom, i_log_dom, nvl, var_list_typ, vl_list(nvar_lists)
    CHARACTER(LEN=2) :: lev_type
    TYPE (t_output_name_list), POINTER :: p_onl
    TYPE (t_output_file), POINTER :: p_of
    TYPE(t_list_element), POINTER :: element

    CALL message('init_name_list_output','Start')

    ! Preliminary until this is done correctly during creating the var_lists:
    ! Set loutput to .TRUE. unless we know exactly that we do NOT want to output this list

    DO i = 1, nvar_lists

      var_lists(i)%p%loutput = .TRUE.

      ! Set disable output for var_lists which must not be output,
      ! e.g. because of name clashes with other var_lists, e.g.:
      !IF(var_lists(i)%p%name(1:15) == 'ext_data_atm_td' ) var_lists(i)%p%loutput = .FALSE.
      !IF(var_lists(i)%p%name(1:16) == 'nh_state_metrics') var_lists(i)%p%loutput = .FALSE.

! For a list of all variables, enable the following!
IF(.TRUE.) THEN
      IF (my_process_is_stdio()) THEN
        PRINT '(3a, i2)','Var_list name: ',TRIM(var_lists(i)%p%name),' Patch: ',var_lists(i)%p%patch_id
        element => var_lists(i)%p%first_list_element
        DO
          IF(.NOT. ASSOCIATED(element)) EXIT
          PRINT *,'    ',element%field%info%name,element%field%info%loutput
          element => element%next_list_element
        ENDDO
      ENDIF
ENDIF

    ENDDO


    ! Get the number of output files needed (by counting the domains per name list)

    p_onl => first_output_name_list
    nfiles = 0

    DO

      IF(.NOT.ASSOCIATED(p_onl)) EXIT

      ! Hack, done here since we need the number of physical domains
      IF(p_onl%dom(1) <= 0) THEN
        ! This means to include all domains
        IF(l_output_phys_patch) THEN
          DO i = 1, n_phys_dom
            p_onl%dom(i) = i
          ENDDO
        ELSE
          DO i = 1, n_dom
            p_onl%dom(i) = i
          ENDDO
        ENDIF
      ENDIF

      DO i = 1, SIZE(p_onl%dom)
        IF(p_onl%dom(i) <= 0) EXIT ! Last one was reached
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

    ! Loop over all output namelists, set up the output_file struct for all associated files

    p_onl => first_output_name_list
    nfiles = 0

    DO

      IF(.NOT.ASSOCIATED(p_onl)) EXIT

      ! Loop over all domains for which this name list should be used

      DO i = 1, SIZE(p_onl%dom)
        IF(p_onl%dom(i) <= 0) EXIT ! Last one was reached
        i_dom = p_onl%dom(i)
        IF(l_output_phys_patch) THEN
          i_log_dom = p_phys_patch(i_dom)%logical_id
        ELSE
          i_log_dom = i_dom
        ENDIF

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

          p_of%patch_id     = i_dom
          p_of%output_type  = p_onl%filetype
          p_of%name_list    => p_onl

          ! Select all var_lists which belong to current logical domain and i_typ

          nvl = 0
          DO j = 1, nvar_lists

            IF(.NOT. var_lists(j)%p%loutput) CYCLE
            IF(var_lists(j)%p%patch_id /= i_log_dom) CYCLE

            ! The following has to be improved, we need a better way to
            ! recognize pressure/height var_lists
            var_list_typ = 1
            IF(var_lists(j)%p%name(1:15) == 'nh_state_diag_p') var_list_typ = 2
            IF(var_lists(j)%p%name(1:15) == 'nh_state_diag_z') var_list_typ = 3

            IF(i_typ /= var_list_typ) CYCLE

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

    CALL message('init_name_list_output','Done')

  END SUBROUTINE init_name_list_output

  !------------------------------------------------------------------------------------------------

  SUBROUTINE add_varlist_to_output_file(p_of, vl_list, varlist)

    TYPE (t_output_file), INTENT(INOUT) :: p_of
    INTEGER, INTENT(IN) :: vl_list(:)
    CHARACTER(LEN=*), INTENT(IN) :: varlist(:)

    INTEGER :: ivar, i, iv, idx
    TYPE(t_var_metadata), POINTER :: p_info_1st
    TYPE(t_list_element), POINTER :: element
    TYPE(t_olist_element), POINTER :: p_elem

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/add_varlist_to_output_file'

    p_of%first_list_element => NULL()
    p_elem => NULL()

    ! Loop over all variables in varlist

    DO ivar = 1, SIZE(varlist)

      IF(varlist(ivar) == ' ') EXIT ! Last one reached

      ! Loop over all var_lists listed in vl_list to find the variable
      ! Please note that there may be several variables with different time levels,
      ! we just add unconditionally all with the name varlist(ivar).
      ! Remark: The different time levels may appear in different lists
      ! or in the same list, the code sill accept both

      p_info_1st => NULL()

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

          IF(.NOT.ASSOCIATED(p_of%first_list_element)) THEN
            ALLOCATE(p_of%first_list_element)
            p_elem => p_of%first_list_element
          ELSE
            ALLOCATE(p_elem%next_list_element)
            p_elem => p_elem%next_list_element
          ENDIF

          ! Please note: The following statement requires that the complete
          ! structure (including pointers etc.) is copied.
          ! If that doesn't work for any compiler, it has to be replaced by
          ! single copies of all members

          p_elem%field = element%field

          ! Additional members in p_elem

          IF(idx == 0) THEN
            p_elem%time_level = -1
          ELSE
            p_elem%field%info%name = element%field%info%name(1:idx-1)
            p_elem%time_level = ICHAR(element%field%info%name(idx+3:idx+3)) - ICHAR('0')
          ENDIF

          ! info_1st: Pointer to the info member of the 1st occurence of this variable name
          ! This is a NULL Pointer for the 1st occurence itself
          p_elem%info_1st => p_info_1st

          p_elem%next_list_element => NULL()

          IF(.NOT.ASSOCIATED(p_info_1st)) p_info_1st => p_elem%field%info

        ENDDO

      ENDDO

      ! Check that at least one element with this name has been found

      IF(.NOT.ASSOCIATED(p_elem)) &
        CALL finish(routine,'Output name list variable not found: '//TRIM(varlist(ivar)))
      IF(p_elem%field%info%name /= varlist(ivar)) &
        CALL finish(routine,'Output name list variable not found: '//TRIM(varlist(ivar)))

    ENDDO

  END SUBROUTINE add_varlist_to_output_file

  !------------------------------------------------------------------------------------------------

  SUBROUTINE open_output_file(of, jfile)

    TYPE(t_output_file), INTENT(INOUT) :: of
    INTEGER, INTENT(IN) :: jfile ! Number of file set to open

    INTEGER :: i, j, k, k_jg, nlev, nlevp1, nplev, nzlev, nzlevp1, znlev_soil, astatus, ip, iv, jv
    INTEGER :: n_cells, n_edges, n_verts
    REAL(wp), ALLOCATABLE :: levels(:), levels_i(:), levels_m(:)

    CHARACTER(LEN=16) :: extn

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/open_output_file'

    ip = of%patch_id
    IF(l_output_phys_patch) THEN
      k_jg = p_phys_patch(ip)%logical_id
      n_cells = p_phys_patch(ip)%n_patch_cells
      n_edges = p_phys_patch(ip)%n_patch_edges
      n_verts = p_phys_patch(ip)%n_patch_verts
    ELSE
      k_jg = ip
      n_cells = p_patch(ip)%n_patch_cells_g
      n_edges = p_patch(ip)%n_patch_edges_g
      n_verts = p_patch(ip)%n_patch_verts_g
    ENDIF

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

    WRITE(of%filename,'(a,"_",i4.4,a)') TRIM(of%filename_pref),jfile,TRIM(extn)
    of%cdiFileID = streamOpenWrite(TRIM(of%filename), of%output_type)

    IF (of%cdiFileID < 0) THEN
      WRITE(message_text,'(a)') cdiStringError(of%cdiFileID)
      CALL message('',message_text)
      CALL finish (routine, 'open failed on '//TRIM(of%filename))
    ELSE
      CALL message (routine, 'opened '//TRIM(of%filename))
    END IF
    !
    ! The following sections add the file global properties collected in init_output
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

    ! Cells

    of%cdiCellGridID = gridCreate(GRID_UNSTRUCTURED, n_cells)
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

    of%cdiVertGridID = gridCreate(GRID_UNSTRUCTURED, n_verts)
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

    of%cdiEdgeGridID = gridCreate(GRID_UNSTRUCTURED, n_edges)
    CALL gridDefNvertex(of%cdiEdgeGridID, 4)
    !
    CALL gridDefXname(of%cdiEdgeGridID, 'elon')
    CALL gridDefXlongname(of%cdiEdgeGridID, 'edge longitude')
    CALL gridDefXunits(of%cdiEdgeGridID, 'radians')
    !
    CALL gridDefYname(of%cdiEdgeGridID, 'elat')
    CALL gridDefYlongname(of%cdiEdgeGridID, 'edge latitude')
    CALL gridDefYunits(of%cdiEdgeGridID, 'radians')

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

      nlev   = num_lev(k_jg)
      nlevp1 = num_levp1(k_jg)
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
        nplev = nh_pzlev_config(k_jg)%nplev
        of%cdiZaxisID(ZA_pressure) = zaxisCreate(ZAXIS_PRESSURE, nplev)
        ALLOCATE(levels(nplev))
        DO k = 1, nplev
          levels(k) = nh_pzlev_config(k_jg)%plevels(k)
        END DO
        CALL zaxisDefLevels(of%cdiZaxisID(ZA_pressure), levels)
        CALL zaxisDefVct(of%cdiZaxisID(ZA_pressure), nplev, levels)
        DEALLOCATE(levels)

        nzlev = nh_pzlev_config(k_jg)%nzlev
        of%cdiZaxisID(ZA_height)  = zaxisCreate(ZAXIS_HEIGHT, nzlev)
        ALLOCATE(levels(nzlev))
        DO k = 1, nzlev
          levels(k) = nh_pzlev_config(k_jg)%zlevels(k)
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
    CALL addGlobAtts(of%cdiVlistID,k_jg,astatus)
    IF (iequations/=ihs_ocean) THEN
      CALL addAtmAtts(of%cdiVlistID,k_jg,astatus)
    ELSE
      CALL addOceAtts(of%cdiVlistID,astatus)
    END IF

    !
    ! add variables
    !
    CALL addVarListToVlist(of)

    CALL streamDefVlist(of%cdiFileID, of%cdiVlistID)

  END SUBROUTINE open_output_file

  !------------------------------------------------------------------------------------------------
  !
  ! define variables and attributes
  !
  SUBROUTINE addVarListToVlist(of)

    TYPE (t_output_file),  INTENT(IN) :: of


    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_olist_element), POINTER :: element
    !
    INTEGER :: vlistID, varID, gridID, zaxisID, k_jg, nlev, nlevp1, znlev_soil
    !
    REAL(wp) :: casted_missval

    vlistID = of%cdiVlistID
    element => NULL()

    IF(l_output_phys_patch) THEN
      k_jg = p_phys_patch(of%patch_id)%logical_id
    ELSE
      k_jg = of%patch_id
    ENDIF
    nlev   = num_lev(k_jg)
    nlevp1 = num_levp1(k_jg)

    ! See above ...
    znlev_soil = SIZE(zml_soil)-1
    !
    for_all_list_elements: DO
      !
      IF(.NOT.ASSOCIATED(element)) THEN
        ! First loop iteration
        element => of%first_list_element
      ELSE
        element => element%next_list_element
      ENDIF
      IF (.NOT.ASSOCIATED(element)) EXIT

      ! Elements for which info_1st is not NULL are not added,
      ! they share the variable description with the corresponding 1st element

      IF(ASSOCIATED(element%info_1st)) CYCLE
      !
      ! retrieve information from actual linked list element
      !
      info => element%field%info
      !
      ! set grid ID
      !
      SELECT CASE (info%hgrid)
      CASE(GRID_UNSTRUCTURED_CELL)
        info%cdiGridID = of%cdiCellGridID
        gridID = info%cdiGridID
      CASE(GRID_UNSTRUCTURED_VERT)
        info%cdiGridID = of%cdiVertGridID
        gridID = info%cdiGridID
      CASE(GRID_UNSTRUCTURED_EDGE)
        info%cdiGridID = of%cdiEdgeGridID
        gridID = info%cdiGridID
      CASE DEFAULT
        CALL finish('addVarListToVlist', 'GRID definition missing for '//TRIM(info%name))
      END SELECT

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
          CALL finish('addVarListToVlist','Dimension mismatch for ZAXIS_DEPTH_BELOW_SEA')
        ENDIF

      CASE (ZAXIS_DEPTH_BELOW_LAND)
        IF (info%used_dimensions(2) == znlev_soil+1) THEN
          info%cdiZaxisID =  of%cdiZaxisID(ZA_depth_below_land)
        ELSE IF (info%used_dimensions(2) == znlev_soil+2) THEN
          info%cdiZaxisID =  of%cdiZaxisID(ZA_depth_below_land_p1)
        ELSE
          PRINT *,'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
          CALL finish('addVarListToVlist','Dimension mismatch for ZAXIS_DEPTH_BELOW_LAND')
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
            CALL finish('addVarListToVlist','Dimension mismatch for ZAXIS_HEIGHT')
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
            CALL finish('addVarListToVlist','Dimension mismatch for ZAXIS_HEIGHT')
          ENDIF
        ENDIF

      CASE (ZAXIS_PRESSURE)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_pressure)

      CASE (ZAXIS_ALTITUDE)
        info%cdiZaxisID =  of%cdiZaxisID(ZA_height)

      CASE DEFAULT
        PRINT *,'Variable: ',TRIM(info%name),' ZAXIS: ',info%vgrid
        CALL finish('addVarListToVlist', 'ZAXIS definition missing for '//TRIM(info%name))

      END SELECT

      zaxisID = info%cdiZaxisID
      !
      info%cdiVarID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE)
      varID = info%cdiVarID 
      !
      CALL vlistDefVarDatatype(vlistID, varID, DATATYPE_FLT32)
      CALL vlistDefVarName(vlistID, varID, info%name)
      !
      IF (info%cf%long_name /= '') CALL vlistDefVarLongname(vlistID, varID, info%cf%long_name)
      IF (info%cf%units /= '') CALL vlistDefVarUnits(vlistID, varID, info%cf%units)

      IF (info%lmiss) THEN
        IF (ASSOCIATED(element%field%r_ptr)) THEN
          casted_missval = info%missval%rval
        ELSE IF (ASSOCIATED(element%field%i_ptr)) THEN
          casted_missval = REAL(info%missval%ival,wp)
        ELSE
          IF (info%missval%lval) THEN
            casted_missval = 1.0_wp
          ELSE
            casted_missval = 0.0_wp
          ENDIF
        ENDIF
        CALL vlistDefVarMissval(vlistID, varID, casted_missval)
      ENDIF
      !
    ENDDO for_all_list_elements
    !
  END SUBROUTINE addVarListToVlist
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE close_name_list_output
    !
    ! Close all name_list files
    !
    INTEGER :: i

    DO i = 1, SIZE(output_file)

      CALL close_output_file(output_file(i))

    ENDDO

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
      CALL gridDestroy(of%cdiCellGridID)
      CALL gridDestroy(of%cdiVertGridID)
      CALL gridDestroy(of%cdiEdgeGridID)
      CALL taxisDestroy(of%cdiTaxisID)
      DO j = 1, SIZE(of%cdiZaxisID)
        IF(of%cdiZaxisID(j) /= CDI_UNDEFID) &
          & CALL zaxisDestroy(of%cdiZaxisID(j))
      ENDDO
      CALL vlistDestroy(vlistID)
    ENDIF

    of%cdiFileID     = CDI_UNDEFID
    of%cdiVlistID    = CDI_UNDEFID
    of%cdiCellGridID = CDI_UNDEFID
    of%cdiVertGridID = CDI_UNDEFID
    of%cdiEdgeGridID = CDI_UNDEFID
    of%cdiTaxisID    = CDI_UNDEFID
    of%cdiZaxisID(:) = CDI_UNDEFID

  END SUBROUTINE close_output_file

  !------------------------------------------------------------------------------------------------
  !
  ! Loop over all output_name_list's, write the ones for which output is due
  ! This routine also cares about opening the output files the first time
  ! and reopening the files after a certain number of steps.
  !
  ! It must be called at the beginning with sim_time=0 in order to write
  ! all name lists for which include_start is set.
  !
  SUBROUTINE write_name_list_output(datetime, sim_time)

    TYPE(t_datetime), INTENT(in) :: datetime
    REAL(wp), INTENT(in)         :: sim_time

    INTEGER :: i, j, idate, itime, iret, n
    TYPE(t_output_name_list), POINTER :: p_onl
    CHARACTER(LEN=filename_max+100) :: text
    REAL(wp), PARAMETER :: eps = 1.d-10 ! Tolerance for checking output bounds

    CALL message('write_name_list_output','Start')

    ! Check if files have to be (re)opened

    ! Go over all output files
    DO i = 1, SIZE(output_file)

      p_onl => output_file(i)%name_list

      ! Check if output is due for this file
      IF((p_onl%include_start .AND. sim_time==0) .OR. &
          p_onl%next_output_time <= sim_time+dtime/2._wp) THEN

        IF(MOD(p_onl%n_output_steps,p_onl%steps_per_file) == 0) THEN
          IF (my_process_is_stdio()) THEN
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
      IF((p_onl%include_start .AND. sim_time==0) .OR. &
          p_onl%next_output_time <= sim_time+dtime/2._wp) THEN

        IF (my_process_is_stdio()) THEN
          CALL taxisDefVdate(output_file(i)%cdiTaxisID, idate)
          CALL taxisDefVtime(output_file(i)%cdiTaxisID, itime)

          iret = streamDefTimestep(output_file(i)%cdiFileId, output_file(i)%cdiTimeIndex)
          output_file(i)%cdiTimeIndex = output_file(i)%cdiTimeIndex + 1
        ENDIF

        WRITE(text,'(a,a,a,1pg15.9)') &
          'Output to ',TRIM(output_file(i)%filename),' at simulation time ',sim_time
        CALL message('mo_name_list_output',text)
        CALL write_name_list(output_file(i))
      ENDIF

    ENDDO

    ! Increase next_output_time for all output name lists which have been written.
    ! Please note that more than 1 output_file may point to the same name_list,
    ! so this must not be done in the loop above!

    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT

      IF((p_onl%include_start .AND. sim_time==0) .OR. &
          p_onl%next_output_time <= sim_time+dtime/2._wp) THEN

        p_onl%n_output_steps = p_onl%n_output_steps + 1

      ENDIF

      ! Switch all name lists to next output time.
      ! This is done in a DO WHILE loop to catch the case
      ! where an output increment less than the time step
      ! or two output_bounds triples which are too close are specified.

      DO WHILE(p_onl%next_output_time <= sim_time+dtime/2._wp)
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

  END SUBROUTINE write_name_list_output

  !------------------------------------------------------------------------------------------------
  !
  ! Write a output name list
  !
  SUBROUTINE write_name_list(of)

    TYPE (t_output_file), INTENT(IN) :: of

    INTEGER :: i_dom, i_log_dom, i, iv, jk, n_points, nlevs, nblks, n, n_points_tot
    INTEGER :: nindex, nmiss, gridtype, fileID, varID
    INTEGER, POINTER :: p_phys(:,:)
    TYPE (t_olist_element), POINTER :: element
    TYPE (t_var_metadata), POINTER :: info
    REAL(wp), POINTER :: r_ptr(:,:,:)
    REAL(wp), ALLOCATABLE :: r_tmp(:,:,:), r_out(:,:), r_out_recv(:,:)
    TYPE(t_comm_pattern), POINTER :: p_pat

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_name_list_output/write_name_list'

    element => NULL()

    ! Just to prevent possible warnings
    n_points_tot = 0
    p_phys => NULL()

    i_dom = of%patch_id
    IF(l_output_phys_patch) THEN
      i_log_dom = p_phys_patch(i_dom)%logical_id
    ELSE
      i_log_dom = i_dom
    ENDIF

    ! Go over all list lelements for this output file

    element => of%first_list_element

    DO

      IF(.NOT.ASSOCIATED(element)) EXIT

      r_ptr => NULL()

      ! For time level dependent elements: check time level

      IF(element%time_level>0 .AND. element%time_level /= nnow(i_log_dom)) THEN
        element => element%next_list_element
        CYCLE
      ENDIF

      ! retrieve information from actual linked list element
      ! if several time levels are present, they share the information of of first one

      IF(ASSOCIATED(element%info_1st)) THEN
        info => element%info_1st
      ELSE
        info => element%field%info
      ENDIF
      !
      IF (info%lcontained) THEN 
        nindex = info%ncontained
      ELSE
        nindex = 1
      ENDIF

      ! Check if first dimension of array is nproma.
      ! Otherwise we got an array which is not suitable for this output scheme.

      IF(info%used_dimensions(1) /= nproma) &
        CALL finish(routine,'1st dim is not nproma: '//TRIM(info%name))

      SELECT CASE (info%ndims)
      CASE (1)
        CALL finish(routine,'1d arrays not handled yet.')
      CASE (2)
        ! Make a 3D copy of the array
        ALLOCATE(r_ptr(info%used_dimensions(1),1,info%used_dimensions(2)))
        r_ptr(:,1,:) = element%field%r_ptr(:,:,nindex,1,1)
      CASE (3)
        ! Just set a pointer to the array
        r_ptr => element%field%r_ptr(:,:,:,nindex,1)
      CASE (4)
        CALL finish(routine,'4d arrays not handled yet.')
      CASE (5)
        CALL finish(routine,'5d arrays not handled yet.')
      CASE DEFAULT 
        CALL finish(routine,'dimension not set.')        
      END SELECT
      !
      gridtype = info%hgrid

      SELECT CASE (gridtype)
      CASE (GRID_UNSTRUCTURED_CELL)
        IF(l_output_phys_patch) THEN
          p_pat => p_phys_patch(i_dom)%comm_pat_gather_c
          n_points = p_phys_patch(i_dom)%n_patch_cells
          n_points_tot = p_patch(i_log_dom)%n_patch_cells_g
          p_phys => p_patch(i_log_dom)%cells%phys_id
        ELSE
          p_pat => p_patch(i_dom)%comm_pat_gather_c
          n_points = p_patch(i_dom)%n_patch_cells_g
        ENDIF
      CASE (GRID_UNSTRUCTURED_EDGE)
        IF(l_output_phys_patch) THEN
          p_pat => p_phys_patch(i_dom)%comm_pat_gather_e
          n_points = p_phys_patch(i_dom)%n_patch_edges
          n_points_tot = p_patch(i_log_dom)%n_patch_edges_g
          p_phys => p_patch(i_log_dom)%edges%phys_id
        ELSE
          p_pat => p_patch(i_dom)%comm_pat_gather_e
          n_points = p_patch(i_dom)%n_patch_edges_g
        ENDIF
      CASE (GRID_UNSTRUCTURED_VERT)
        IF(l_output_phys_patch) THEN
          p_pat => p_phys_patch(i_dom)%comm_pat_gather_v
          n_points = p_phys_patch(i_dom)%n_patch_verts
          n_points_tot = p_patch(i_log_dom)%n_patch_verts_g
          p_phys => p_patch(i_log_dom)%verts%phys_id
        ELSE
          p_pat => p_patch(i_dom)%comm_pat_gather_v
          n_points = p_patch(i_dom)%n_patch_verts_g
        ENDIF
      CASE default
        CALL finish(routine,'unknown grid type')
      END SELECT

      nlevs = UBOUND(r_ptr,2)
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
        IF(l_output_phys_patch) THEN
          n = 0
          DO i = 1, n_points_tot
            IF(p_phys(idx_no(i),blk_no(i)) == i_dom) THEN
              n = n+1
              IF(n>n_points) CALL finish(routine,'n>n_points') ! Safety check only
              r_tmp(idx_no(n),:,blk_no(n)) = r_ptr(idx_no(i),:,blk_no(i))
            ENDIF
          ENDDO
          IF(n/=n_points) CALL finish(routine,'n/=n_points') ! Safety check only
        ELSE
          r_tmp(:,:,1:nblks) = r_ptr(:,:,1:nblks)
        ENDIF
      ELSE
        CALL exchange_data(p_pat, RECV=r_tmp, SEND=r_ptr)
      ENDIF

      ! -----------------------------------------------------
      ! --- Here is the place to do lat/lon interpolation ---
      ! -----------------------------------------------------

      IF(my_process_is_mpi_test() .OR. my_process_is_mpi_workroot()) THEN

        ! De-block the array
        ALLOCATE(r_out(n_points, nlevs))
        DO jk = 1, nlevs
          r_out(:,jk) = RESHAPE(r_tmp(:,jk,:), (/ n_points /))
        ENDDO

        ! In the case of a test run: Compare results on worker PEs and test PE
        IF(p_test_run) THEN
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
        nmiss   = 0
        fileID  = of%cdiFileID
        varID   = info%cdiVarID
        !
        CALL streamWriteVar(fileID, varID, r_out, nmiss)
      END IF
      !
      ! deallocate temporary global arrays
      !
      DEALLOCATE(r_tmp)
      IF(my_process_is_mpi_test() .OR. my_process_is_mpi_workroot()) DEALLOCATE(r_out)
      IF(info%ndims == 2) DEALLOCATE(r_ptr)

      element => element%next_list_element

    ENDDO

  END SUBROUTINE write_name_list

  !------------------------------------------------------------------------------------------------

  FUNCTION istime4name_list_output(sim_time) RESULT(retval)

    REAL(wp), INTENT(IN)  :: sim_time            ! simulation time [s]
    LOGICAL               :: retval

    TYPE(t_output_name_list), POINTER :: p_onl

    retval = .FALSE.

    ! Go through all output name lists and check if next_output_time is reached

    p_onl => first_output_name_list
    DO
      IF(.NOT.ASSOCIATED(p_onl)) EXIT
      IF(p_onl%next_output_time <= sim_time+dtime/2._wp) retval = .TRUE.
      p_onl => p_onl%next
    ENDDO

  END FUNCTION istime4name_list_output

  !------------------------------------------------------------------------------------------------

END MODULE mo_name_list_output
