!>
!! Processing of nudging parameters for the nudging types:
!! - Upper boundary nudging
!! - Global nudging
!! For the lateral boundary nudging, please see: 
!! - src/namelists/mo_limarea_nml
!! - src/namelists/mo_interpol_nml
!!
!! @par Revision History
!! Initial revision by Guenther Zaengl and Sebastian Borchert, DWD (2018-09)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines. 
!!
MODULE mo_nudging_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, inh_atmosphere, inwp
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_nudging_config,      ONLY: indg_type, indg_profile, ithermdyn_type, &
    &                               NDG_VARLIST_STR_LEN, indg_var, cndg_var, &
    &                               nudging_config
  USE mo_util_string,         ONLY: int2string, real2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_nudging_namelist, check_nudging

  ! Module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nudging_nml'

CONTAINS

  !>
  !!
  SUBROUTINE read_nudging_namelist( filename )

    ! In/out variables
    CHARACTER(LEN=*), INTENT(IN) :: filename

    ! Local variables
    INTEGER  :: istat, iunit, funit
    INTEGER  :: jitem, iitem
    LOGICAL  :: lerror, lentry, lapplicable
    REAL(wp), PARAMETER :: eps = dbl_eps * 1000._wp
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: &
      & routine = modname//":read_nudging_namelist"

    ! Defaults, which are reasonable for upper-boundary nudging 
    ! are not necessarily reasonable for global nudging, unfortunately. 
    ! For our approach to handle this, some of the namelist variables 
    ! need to become arrays effectively, although "officially", 
    ! i.e. in the namelist description they go as non-arrays. 
    ! So here, we define some array tools:
    ! The number of entries is the number of nudging types + 2. 
    ! One additional entry is for the actual namelist entry, 
    ! so it gets the index 1 (-> 'inml'). Since we prefer to work 
    ! with the already existing type identifiers 'indg_type%...', 
    ! we introduce the offset 'ishift = 1' (alternatively, we might use 
    ! a start index of zero for the arrays, but that might be more error-prone). 
    ! The other additional entry is for reference, 
    ! so it becomes the last entry (-> 'iref').
    INTEGER, PARAMETER :: initem = indg_type%nitem + 2
    INTEGER, PARAMETER :: ishift = 1
    INTEGER, PARAMETER :: inml   = 1
    INTEGER, PARAMETER :: iref   = initem

    ! Namelist variables:.............................................................................................
    INTEGER  :: nudge_type                           ! Nudging type identifier:
                                                     ! - 0 -> switched off
                                                     ! - 1 -> upper boundary nudging
                                                     ! - 2 -> global nudging
    REAL(wp) :: max_nudge_coeff_vn(initem)           ! Max. nudging coefficient for horizontal wind
    REAL(wp) :: max_nudge_coeff_thermdyn(initem)     ! Max. nudging coefficient for thermodynamic variables
    REAL(wp) :: nudge_start_height(initem)           ! Nuding start height
    INTEGER  :: thermdyn_type(initem)                ! Type of thermodynamic variables (from the perspective
                                                     ! of the nudging-target variables in ICON):
                                                     ! - 1 -> hydrostatic: hydrostatic pressure & temperature are nudged
                                                     !        (diagnostic variables)
                                                     ! - 2 -> nonhydrostatic: density & virtual potential temperature are nudged
                                                     !        (prognostic variables)
    INTEGER  :: nudge_profile(initem)                ! Profile of nudging strength:
                                                     ! - 1 -> squared scaled vertical distance from nudging start height
                                                     ! - 2 -> constant profile (for global nudging only)
                                                     ! - 3 -> hyperbolic tangent profile (for global nudging only)
                                                     ! - 4 -> trapezoidal profile (for global nudging only)
    REAL(wp) :: nudge_end_height(initem)             ! Nudging end height
    REAL(wp) :: nudge_scale_height(initem)           ! Scale height for profiles
    CHARACTER(LEN=NDG_VARLIST_STR_LEN) :: nudge_var(initem)  ! Nudging variable(s) (for global nudging only): 
                                                     ! - 'vn'       -> horizontal wind 
                                                     ! - 'thermdyn' -> thermodynamic variables
                                                     ! - 'qv'       -> water vapour
                                                     ! - string with a comma-separated list with variables, e.g., 'vn,thermdyn'
                                                     ! - 'all'      -> all variables (i.e., 'vn,thermdyn,qv')
    REAL(wp) :: max_nudge_coeff_qv(initem)           ! Max. nudging coefficient for water vapour (for global nudging only)
    INTEGER  :: idiagnose(initem)                    ! Multiple of 'run_nml: dtime' for diagnosing the nudging success/impact 
                                                     ! during runtime (i.e., output quality measures)
    !.................................................................................................................

    
    NAMELIST /nudging_nml/ nudge_type, max_nudge_coeff_vn, max_nudge_coeff_thermdyn,         &
      &                    nudge_start_height, nudge_var, nudge_profile, max_nudge_coeff_qv, &
      &                    nudge_end_height, nudge_scale_height, thermdyn_type, idiagnose    

    !----------------------------------------

    !------------------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------------------

    nudge_type                = indg_type%off    ! 0 -> nudging switched off

    ! For the following variables, 
    ! it might be reasonable to indtroduce 
    ! type-dependent defaults

    ! Initialization 
    max_nudge_coeff_vn(:)       = -999._wp
    max_nudge_coeff_thermdyn(:) = -999._wp
    max_nudge_coeff_qv(:)       = -999._wp
    nudge_start_height(:)       = -999._wp
    nudge_end_height(:)         = -999._wp
    nudge_scale_height(:)       = -999._wp 
    thermdyn_type(:)            = -ithermdyn_type%nitem 
    nudge_profile(:)            = -indg_profile%nitem
    nudge_var(:)                = " "
    idiagnose(:)                = -999 

    ! For upper-boundary nudging
    ! --------------------------
    iitem = indg_type%ubn + ishift
    !
    ! Namelist-controllable
    max_nudge_coeff_vn(iitem)       = 0.04_wp                    ! Max. nudging coefficient for horizontal wind
    max_nudge_coeff_thermdyn(iitem) = 0.075_wp                   ! Max. nudging coefficient for thermodynamic variables
    nudge_start_height(iitem)       = 12000._wp                  ! (m) Nudging start height
    ! Not namelist-controllable
    ! (technically, a namelist entry would be possible, but is disregarded)
    thermdyn_type(iitem)            = ithermdyn_type%hydrostatic ! 1 -> diagnostic variables: hydrostatic pressure 
                                                                 !      and temperature are nudged 
    nudge_profile(iitem)            = indg_profile%sqrddist      ! 1 -> inverse squared scaled vertical distance from start height
    nudge_end_height(iitem)         = -1._wp                     ! (m) Nudging end height,
                                                                 ! negative value means: 'nudge_end_height = top_height' 
                                                                 ! ('sleve_nml: top_height' -> height of model top)
    nudge_scale_height(iitem)       = -1._wp                     ! (m) Scale height for nudging strength profile,
                                                                 ! negative value means: 
                                                                 ! 'nudge_scale_height = nudge_end_height - nudge_start_height'
    max_nudge_coeff_qv(iitem)       = 0._wp                      ! Max. nudging coefficient for water vapour
    nudge_var(iitem)                = "vn,thermdyn"              ! Horizontal wind and thermodynamic variables are nudged
    idiagnose(iitem)                = -1                         ! No analysis of nudging success/impact

    ! For global nudging
    ! ------------------
    iitem = indg_type%globn + ishift
    !
    ! Namelist-controllable
    max_nudge_coeff_vn(iitem)       = 0.016_wp                     ! Max. nudging coefficient for horizontal wind
    max_nudge_coeff_thermdyn(iitem) = 0.03_wp                      ! Max. nudging coefficient for thermodynamic variables
    max_nudge_coeff_qv(iitem)       = 0.008_wp                     ! Max. nudging coefficient for water vapour
    nudge_start_height(iitem)       = 2000._wp                     ! (m) Nudging start height 
                                                                   ! (maybe better omit the planetary boundary layer)
    nudge_end_height(iitem)         = 40000._wp                    ! (m) Nudging end height
                                                                   ! (IFS analysis as driving data: exclude most of 
                                                                   ! that sponge layer of the IFS, which starts above 10 hPa)
    nudge_scale_height(iitem)       = 3000._wp                     ! (m) Scale height for nudging strength profile
    nudge_profile(iitem)            = indg_profile%trapid          ! 4 -> trapezoidal profile
    thermdyn_type(iitem)            = ithermdyn_type%hydrostatic   ! 1 -> diagnostic variables: hydrostatic pressure 
                                                                   !      and temperature are nudged 
    nudge_var(iitem)                = TRIM(cndg_var(indg_var%all)) ! 'all' -> nudge all variables
    idiagnose(iitem)                = -1                           ! Negative number -> no analysis of nudging success/impact

    !------------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !------------------------------------------------------------------------

    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('nudging_nml')
      READ(funit,NML=nudging_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml ('nudging_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, nudging_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, nudging_nml)
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, nudging_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !---------------------------------------------------------
    ! Check and store the namelist entry for the nudging type
    !---------------------------------------------------------
    
    IF (nudge_type /= indg_type%off) THEN
      lerror = .TRUE.
      ! Check 'nudge_type'
      DO jitem = 1, indg_type%nitem
        IF (nudge_type == jitem) THEN
          lerror = .FALSE.
          EXIT
        ENDIF
      ENDDO  !jitem
      IF (lerror) CALL finish(TRIM(routine), 'Invalid nudge_type: '//TRIM(int2string(nudge_type)))
    ENDIF  !IF (nudge_type /= indg_type%off)

    nudging_config%nudge_type = nudge_type

    !---------------------------------------------------------------
    ! Store the other namelist variables in the nudging config type
    !---------------------------------------------------------------

    jitem       = nudge_type + ishift
    lapplicable = .TRUE.
    ! Max. nudging coefficient for horizontal wind
    lentry = lapplicable .AND. (ABS(max_nudge_coeff_vn(inml) - max_nudge_coeff_vn(iref)) > eps)
    iitem  = MERGE( inml,   & ! .TRUE.  -> take namelist entry
      &             jitem,  & ! .FALSE. -> take default value
      &             lentry  ) ! Is there a (valid) namelist entry?
    nudging_config%max_nudge_coeff_vn = max_nudge_coeff_vn(iitem)
    ! Max. nudging coefficient for thermodynamic variables
    lentry = lapplicable .AND. (ABS(max_nudge_coeff_thermdyn(inml) - max_nudge_coeff_thermdyn(iref)) > eps)
    iitem  = MERGE(inml, jitem, lentry)
    nudging_config%max_nudge_coeff_thermdyn = max_nudge_coeff_thermdyn(iitem)
    ! Nudging start height
    lentry = lapplicable .AND. (ABS(nudge_start_height(inml) - nudge_start_height(iref)) > eps)
    iitem  = MERGE(inml, jitem, lentry)
    nudging_config%nudge_start_height = nudge_start_height(iitem)
    ! 
    ! The following namelist variables are for global nudging only
    lapplicable = nudge_type == indg_type%globn
    ! Max. nudging coefficient for water vapour
    lentry = lapplicable .AND. (ABS(max_nudge_coeff_qv(inml) - max_nudge_coeff_qv(iref)) > eps)
    iitem  = MERGE(inml, jitem, lentry)
    nudging_config%max_nudge_coeff_qv = max_nudge_coeff_qv(iitem)
    ! Nudging end height
    lentry = lapplicable .AND. (ABS(nudge_end_height(inml) - nudge_end_height(iref)) > eps)
    iitem  = MERGE(inml, jitem, lentry)    
    nudging_config%nudge_end_height = nudge_end_height(iitem)
    ! Scale height for nudging strength profile
    lentry = lapplicable .AND. (ABS(nudge_scale_height(inml) - nudge_scale_height(iref)) > eps)
    iitem  = MERGE(inml, jitem, lentry)    
    nudging_config%nudge_scale_height = nudge_scale_height(iitem)
    ! Type of thermodynamic variables 
    lentry = lapplicable .AND. (thermdyn_type(inml) /= thermdyn_type(iref))
    iitem  = MERGE(inml, jitem, lentry)    
    nudging_config%thermdyn_type = thermdyn_type(iitem)
    ! Profile of nudging strength
    lentry = lapplicable .AND. (nudge_profile(inml) /= nudge_profile(iref))
    iitem  = MERGE(inml, jitem, lentry)    
    nudging_config%nudge_profile = nudge_profile(iitem)
    ! Nudging variable(s)
    lentry = lapplicable .AND. (TRIM(nudge_var(inml)) /= TRIM(nudge_var(iref)))
    iitem  = MERGE(inml, jitem, lentry)
    nudging_config%nudge_var = nudge_var(iitem)
    ! dtime-sampling interval for analysis of nudging success/impact
    lentry = lapplicable .AND. (idiagnose(inml) /= idiagnose(iref))
    iitem  = MERGE(inml, jitem, lentry)
    nudging_config%idiagnose = idiagnose(iitem) 

    !------------------------------------------------------------------------
    ! Consistency checks
    !------------------------------------------------------------------------

    ! For reasons of clarity, we check nudging-type-wise, 
    ! although this means some code redundancy.

    ! Upper-boundary nudging
    ! ----------------------
    IF (nudge_type == indg_type%ubn) THEN
      !
      IF (nudging_config%max_nudge_coeff_vn < 0._wp) THEN
        CALL finish(TRIM(routine), 'max_nudge_coeff_vn >= 0 required.')
      ELSEIF(nudging_config%max_nudge_coeff_thermdyn < 0._wp) THEN
        CALL finish(TRIM(routine), 'max_nudge_coeff_thermdyn >= 0 required.')
      ELSEIF (nudging_config%nudge_start_height < 0._wp) THEN
        CALL finish(TRIM(routine), 'nudge_start_height >= 0 required.') 
      ENDIF
    !
    ! Global nudging
    ! --------------
    ELSEIF (nudge_type == indg_type%globn) THEN
      !
      ! Check 'max_nudge_coeff_...'
      IF (nudging_config%max_nudge_coeff_vn < 0._wp) THEN
        CALL finish(TRIM(routine), 'max_nudge_coeff_vn >= 0 required.')
      ELSEIF(nudging_config%max_nudge_coeff_thermdyn < 0._wp) THEN
        CALL finish(TRIM(routine), 'max_nudge_coeff_thermdyn >= 0 required.')
      ELSEIF(nudging_config%max_nudge_coeff_qv < 0._wp) THEN
        CALL finish(TRIM(routine), 'max_nudge_coeff_qv >= 0 required.')
      ENDIF
      ! 
      ! For global nudging there is no check of start_height, 
      ! end_height and scale_height, since for reasons of convenience 
      ! negative values mean:
      ! - start_height -> z = 0, nudging starts in the lowermost grid layer
      ! - end_height   -> z = top_height, nudging ends in the uppermost grid layer
      ! - scale_height -> end_height - start_height
      ! But later on, we will check, if start_height < end_height
      !
      ! Check 'thermdyn_type'
      lerror = .TRUE.
      DO jitem = 1, ithermdyn_type%nitem
        IF (nudging_config%thermdyn_type == jitem) THEN
          lerror = .FALSE.
          EXIT
        ENDIF
      ENDDO  !jitem
      IF (lerror) CALL finish(TRIM(routine), 'Invalid thermdyn_type: ' &
        & //TRIM(int2string(nudging_config%thermdyn_type)))
      !
      ! Check 'nudge_profile'
      lerror = .TRUE.
      DO jitem = 1, indg_profile%nitem
        IF (nudging_config%nudge_profile == jitem) THEN
          lerror = .FALSE.
          EXIT
        ENDIF
      ENDDO  !jitem
      IF (lerror) CALL finish(TRIM(routine), 'Invalid nudge_profile: ' &
        & //TRIM(int2string(nudging_config%nudge_profile))) 
      !
      ! Check 'nudge_var'
      lerror = .TRUE.
      VAR_LOOP : DO jitem = 1, indg_var%nitem_2  ! (in contrast to 'nitem', 'nitem_2' includes the entry "all")
        IF (INDEX(TRIM(nudging_config%nudge_var), TRIM(cndg_var(jitem))) > 0) THEN
          lerror = .FALSE.
          EXIT VAR_LOOP
        ENDIF
      ENDDO VAR_LOOP
      IF (lerror) THEN
        CALL finish(TRIM(routine), 'Invalid nudge_var: '//TRIM(nudging_config%nudge_var))            
      ENDIF
    ENDIF  !nudging type

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------

    IF(my_process_is_stdio()) THEN
      funit = open_tmpfile()
      WRITE(funit,NML=nudging_nml)
      CALL store_and_close_namelist(funit, 'nudging_nml')
    ENDIF

    !-----------------------------------------------------
    ! Write the contents of the namelist to an ASCII file
    !-----------------------------------------------------

    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nudging_nml)

  END SUBROUTINE read_nudging_namelist

  !...................................................................................................................

  !>
  !! Subroutine to crosscheck the nudging namelist.
  !!
  SUBROUTINE check_nudging( n_dom, iequations, iforcing, ivctype, top_height,               &
    &                       l_limited_area, num_prefetch_proc, lsparse_latbc, itype_latbc,  &
    &                       nudge_hydro_pres, latbc_varnames_map_file, LATBC_TYPE_CONST,    & 
    &                       LATBC_TYPE_EXT, is_plane_torus, lart,ndyn_substeps, ltransport, &
    &                       nsteps, msg_level                                               )

    ! In/out variables
    ! (In order to avoid circular dependencies, 
    ! we prefer to keep the number of USE-bindings small. 
    ! This is the reason, why there are so many in/out variables.)
    INTEGER,          INTENT(IN) :: n_dom                   ! Number of model domains, 1=primary domain only 
    INTEGER,          INTENT(IN) :: iequations              ! Type of model equations
    INTEGER,          INTENT(IN) :: iforcing                ! Physics package
    INTEGER,          INTENT(IN) :: ivctype                 ! Type of vertical grid
    REAL(wp),         INTENT(IN) :: top_height              ! Height of model top
    LOGICAL,          INTENT(IN) :: l_limited_area          ! .TRUE.: limited-area mode is switched on
    INTEGER,          INTENT(IN) :: num_prefetch_proc       ! Number of PEs for prefetching the driving data
    LOGICAL,          INTENT(IN) :: lsparse_latbc           ! .TRUE.: only nudging data for boundary rows are read in
    INTEGER,          INTENT(IN) :: itype_latbc             ! Type of limited area boundary nudging
    LOGICAL,          INTENT(IN) :: nudge_hydro_pres        ! .TRUE.: use hydrostatic pressure for lateral boundary nudging
    CHARACTER(LEN=*), INTENT(IN) :: latbc_varnames_map_file ! Dictionary which maps internal variable names onto GRIB2
                                                            ! shortnames or NetCDF var names used in lateral boundary nudging
    INTEGER,          INTENT(IN) :: LATBC_TYPE_CONST        ! Identifier for itype_latbc = constant nudging target data
    INTEGER,          INTENT(IN) :: LATBC_TYPE_EXT          ! Identifier for itype_latbc = time-dependent nudging target data
    LOGICAL,          INTENT(IN) :: is_plane_torus          ! .TRUE.: torus mode
    LOGICAL,          INTENT(IN) :: lart                    ! .TRUE.: ART interface cut in
    INTEGER,          INTENT(IN) :: ndyn_substeps           ! Number of dynamics substeps per fast-physics step
    LOGICAL,          INTENT(IN) :: ltransport              ! .TRUE.: tracer transport switched on 
    INTEGER,          INTENT(IN) :: nsteps                  ! Number of time steps
    INTEGER,          INTENT(IN) :: msg_level               ! Verbosity level

    ! Local variables
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: &
      & routine = modname//":check_nudging"

    !----------------------------------------

    !---------------------------------------------
    ! Crosscheck supplementaries
    !---------------------------------------------

    ! Constant-in-time lateral boundary data?
    nudging_config%llatbc_type_const = itype_latbc == LATBC_TYPE_CONST

    ! Main nudging switch
    nudging_config%lnudging = ( ( nudging_config%nudge_type == indg_type%ubn .AND.        &
      &                           l_limited_area                                   ) .OR. &
      &                         ( nudging_config%nudge_type == indg_type%globn )          )

    ! Asynchronous read-in of driving data
    nudging_config%lasync = num_prefetch_proc > 0
    
    !---------------------------------------------
    ! Crosschecks
    !---------------------------------------------

    IF (nudging_config%lnudging) THEN

      ! Abort criteria:

      IF (iequations /= inh_atmosphere) THEN
        ! Nudging only in case of the non-hydrostatic model
        WRITE(message_text,'(a)') "Nudging requires: iequations = "//TRIM(int2string(inh_atmosphere))
        CALL finish(routine, message_text)
      ELSEIF (iforcing /= inwp) THEN
        ! Nudging only in case of the NWP physics package
        WRITE(message_text,'(a)') "Nudging requires: iforcing = "//TRIM(int2string(inwp))
        CALL finish(routine, message_text)
      ELSEIF (ivctype /= 2) THEN
        ! Nudging only in case of a vertical grid of the SLEVE-type
        CALL finish(routine, "Nudging requires: ivctype = 2")
      ELSEIF (is_plane_torus) THEN
        ! Nudging is not foreseen for the torus mode
        CALL finish(routine, "Nudging requires: is_plane_torus = .false.")
      ELSEIF ((nudging_config%nudge_type == indg_type%ubn) .AND. (.NOT. l_limited_area)) THEN
        ! Upper boundary nudging requires the limited-area mode
        ! (Note: this check is actually redundant, but we keep it as a safety net)
        CALL finish(routine, "Upper boundary nudging requires: l_limited_area = .true.")
      ELSEIF ((nudging_config%nudge_type == indg_type%ubn) .AND. lsparse_latbc) THEN
        ! For upper boundary nudging it is essential that "latbc"-data are read in 
        ! on the entire limited-area domain, not only on a small strip along the lateral boundary
        CALL finish(routine, "Upper boundary nudging requires read-in of nudging-data on entire domain.")
      ELSEIF ((nudging_config%nudge_type == indg_type%globn) .AND. lsparse_latbc) THEN
        ! For global nudging it is essential that "latbc"-data are read in on the entire domain
        CALL finish(routine, "Global nudging requires read-in of nudging-data on entire domain.")
      ELSEIF ((nudging_config%nudge_type == indg_type%globn) .AND. &
        &     (itype_latbc /= LATBC_TYPE_EXT)) THEN
        ! Global nudging requires time-dependent driving data
        WRITE(message_text,'(a)') "Global nudging requires: itype_latbc = "//TRIM(int2string(LATBC_TYPE_EXT))
        CALL finish(routine, message_text)
      ELSEIF ((nudging_config%nudge_type == indg_type%globn) .AND. nudging_config%lasync .AND. &
        &     (LEN_TRIM(latbc_varnames_map_file) == 0)) THEN
        ! In case of asynchronous read-in of driving data in a limited-area configuration, 
        ! the presence of a variable dictionary file is not mandatory 
        ! (see 'src/configure_model/mo_limarea_config: configure_latbc'). 
        ! However, when 'src/io/atmo/mo_async_latbc: check_variables' checks for a consistent set 
        ! of variables in the driving data file, the allowed variation of short names is understandably limited 
        ! (e.g. in case of pressure and temperature). So without a dictionary file you may quickly 
        ! run into the situation that this subroutine does not recognize a variable in the data file,  
        ! although it is there (for instance, you may have stored temperature under the short name "T" 
        ! in the driving data files, but 'check_variables' locks for "TEMP"). 
        ! For this reason, we decided that a dictionary file is mandatory in case of global nudging.
        ! Synchronous read-in does not use the dictionary file, but provides different reading subroutines 
        ! in 'src/io/atmo/mo_sync_latbc'.
        WRITE(message_text,'(a)') "Global nudging with num_prefetch_proc = 1 requires to specify " &
          & //"a dictionary file in latbc_varnames_map_file (try icon/run/dict.latbc)."
        CALL finish(routine, message_text)
      ELSEIF ((nudging_config%idiagnose > 0) .AND. (msg_level < 11)) THEN
        ! The analysis of the nudging impact includes the global mean of
        ! the absolute surface pressure tendency. However, this quantity is diagnosed 
        ! only if 'msg_level >= 11'. 
        ! ('msg_level' might be regarded as a "secondary" parameter, 
        ! which we could simply set to the desired value here. 
        ! However, for safety reasons, we prefer not to touch parameters from other namelists.)
        CALL finish(routine, "idiagnose > 0 requires msg_level >= 11.")
      ENDIF

      ! Warning criteria:

      IF (n_dom > 1) THEN
        ! Grid nesting + nudging has not been checked yet
        WRITE(message_text,'(a)') 'WARNING, compatibility of grid nesting and nudging has not been checked!'
        CALL message(TRIM(routine), message_text)
      ENDIF
      IF (lart) THEN
        ! Nudging + enabled ART-interface has not been checked yet
        WRITE(message_text,'(a)') 'WARNING, compatibility of ART and nudging has not been checked!'
        CALL message(TRIM(routine), message_text)        
      ENDIF
      IF (nudging_config%max_nudge_coeff_vn > 1._wp / REAL(ndyn_substeps,wp)) THEN
        ! The nudging coefficient should not be too large
        WRITE(message_text,'(2a)') 'WARNING, max_nudge_coeff_vn is quite large! ', &
          & 'It is recommended, to choose values according to: 0 <= max_nudge_coeff_vn << ' &
          & //TRIM(real2string(1._wp / REAL(ndyn_substeps,wp)))
        CALL message(TRIM(routine), message_text)                
      ENDIF
      IF (nudging_config%max_nudge_coeff_thermdyn > 1._wp / REAL(ndyn_substeps,wp)) THEN
        WRITE(message_text,'(2a)') 'WARNING, max_nudge_coeff_thermdyn is quite large! ', &
          & 'It is recommended, to choose values according to: 0 <= max_nudge_coeff_thermdyn << ' &
          & //TRIM(real2string(1._wp / REAL(ndyn_substeps,wp)))
        CALL message(TRIM(routine), message_text)                
      ENDIF
      IF ((nudging_config%nudge_type == indg_type%globn) .AND.                &
        & (nudging_config%max_nudge_coeff_qv > 1._wp / REAL(ndyn_substeps,wp))) THEN
        ! (Applies to global nudging only)
        WRITE(message_text,'(2a)') 'WARNING, max_nudge_coeff_qv is quite large! ', &
          & 'It is recommended, to choose values according to: 0 <= max_nudge_coeff_qv << ' &
          & //TRIM(real2string(1._wp / REAL(ndyn_substeps,wp)))
        CALL message(TRIM(routine), message_text)                
      ENDIF
      IF ((nudging_config%nudge_type == indg_type%globn) .AND. l_limited_area) THEN
        ! The web of infrastructure related to the limited-area mode has a high degree of complexity  
        ! and we cannot yet exclude negative interferences between the limited-area mode and global nudging. 
        ! Parts of the code that might still require modification are (non-exhaustive list!):
        ! * src/atm_dyn_iconam/mo_nh_diffusion: diffusion
        ! * src/atm_dyn_iconam/mo_solve_nonhydro: solve_nh (damping)
        ! * src/atm_dyn_iconam/mo_nh_deepatmo_solve: solve_nh_deepatmo (damping)
        ! * src/atm_dyn_iconam/mo_vertical_grid: set_nh_metrics ('nudge_c/e_idx/blk')
        ! * src/shr_horizontal/mo_intp_coeffs: init_nudgecoeffs ('nudgecoeff_c/e')
        ! * src/atm_dyn_iconam/mo_nh_stepping: integrate_nh:
        !   We disabled lateral boundary nudging, if global nudging is switched on. 
        !   This is because the measures that would be necessary to avoid double counting 
        !   of nudging tendencies in the lateral boundary nudging zone 
        !   would increase the computational costs and inflate the amount of code.
        WRITE (message_text,'(a)') &
          & 'WARNING, application of global nudging in limited area mode is in experimental stage! ...'
        CALL message(TRIM(routine), message_text)
        WRITE (message_text,'(a)') '... Additional lateral boundary nudging is not applied! ...'
        CALL message(TRIM(routine), message_text)
        WRITE (message_text,'(a)') '... Increased diffusion in lateral boundary nudging zone is applied! ...'
        CALL message(TRIM(routine), message_text)
        WRITE (message_text,'(a)') '... Use at your own risk!'
        CALL message(TRIM(routine), message_text)
      ENDIF

      ! Mixed criteria 
      ! + crosscheck-induced harmonizations:

      ! For reasons of convenience and safety, 
      ! we check nudging-type-wise, 
      ! although this means some code redundancy.

      ! Upper-boundary nudging
      ! ----------------------
      IF (nudging_config%nudge_type == indg_type%ubn) THEN 
        ! The nudging end height is always at the model top
        nudging_config%nudge_end_height = top_height
        !
        ! If the nudging start height is > 0, 
        ! has been checked in 'read_nudging_namelist' above
        !
        ! We have to check, if 'nudge_start_height < nudge_end_height' holds
        IF (nudging_config%nudge_start_height > nudging_config%nudge_end_height) THEN
          WRITE(message_text,'(2a)') 'nudge_start_height > top_height! ', &
            & 'Please, choose a value according to: nudge_start_height < top_height = ' &
            & //TRIM(real2string(top_height))//' m'
          CALL finish(routine, message_text)
        ENDIF  !IF (nudging_config%nudge_start_height > nudging_config%nudge_end_height)
        !
        ! The nudging scale height will be computed in 'src/configure_model/mo_nudging_config: configure_nudging'.
        !
        ! There are two options for the thermodynamic nudging variables:
        ! (1) hydrostatic pressure and temperature 
        ! (2) density and virtual potential temperature
        ! In case of upper boundary nudging, we should be consistent with the choice of 'latbc_config%nudge_hydro_pres'
        ! (see, e.g., the processes controlled by 'latbc_config%nudge_hydro_pres' in 'src/atm_dyn_iconam/mo_nh_stepping')
        nudging_config%thermdyn_type = MERGE( ithermdyn_type%hydrostatic,    & ! .TRUE.  -> option (1)
          &                                   ithermdyn_type%nonhydrostatic, & ! .FALSE. -> option (2)
          &                                   nudge_hydro_pres               ) ! Condition
        ! 
        ! The option to nudge the thermodynamic variables: hydrostatic pressure & temperature 
        ! is not available, if the nudging target data are constant in time
        IF ( nudging_config%llatbc_type_const .AND.                     &
          &  nudging_config%thermdyn_type == ithermdyn_type%hydrostatic ) THEN
          ! If upper boundary nudging is switched on, we have to abort, 
          ! since we should not change 'latbc_config%nudge_hydro_pres' at this place
          CALL finish(routine, "Nudging target data are constant in time. limarea_nml: nudge_hydro_pres = .false. required")
        ENDIF
      !
      ! Global nudging
      ! --------------
      ELSEIF (nudging_config%nudge_type == indg_type%globn) THEN
        IF ( nudging_config%nudge_end_height < 0._wp      .OR. &
          &  nudging_config%nudge_end_height > top_height      ) THEN
          ! We set end_height = top_height
          nudging_config%nudge_end_height = top_height
          WRITE(message_text,'(a)') 'Nudging end height set to nudge_end_height = top_height = ' &
            & //TRIM(real2string(top_height))//' m'
          CALL message(TRIM(routine), message_text)
        ENDIF  !IF (nudging_config%nudge_end_height ...)
        !
        IF ( nudging_config%nudge_start_height < 0._wp      .OR. &
          &  nudging_config%nudge_start_height > top_height      ) THEN
          ! We set start_height = 0 
          nudging_config%nudge_start_height = 0._wp
          WRITE(message_text,'(a)') 'Nudging start height set 0 m' 
          CALL message(TRIM(routine), message_text)
        ENDIF  !IF (nudging_config%nudge_start_height ...)    
        !
        ! Check, if 'nudge_start_height < nudge_end_height' holds
        IF (nudging_config%nudge_start_height > nudging_config%nudge_end_height) THEN
          WRITE(message_text,'(2a)') 'nudge_start_height > nudge_end_height! ',                                &
            & 'Please, choose values according to: 0 < nudge_start_height < nudge_end_height <= top_height = ' &
            & //TRIM(real2string(top_height))//' m'
          CALL finish(routine, message_text)
        ENDIF  !IF (nudging_config%nudge_start_height > nudging_config%nudge_end_height)
        !
        ! The nudging scale height will be computed in 'src/configure_model/mo_nudging_config: configure_nudging'.
        !
        ! Nudging water vapour is not possible, if tracer transport is switched off 
        nudging_config%lqv_nudgable = ltransport
        !
        ! If analysis of nudging success/impact is switched on, 
        ! the analysis interval in units of the time step 'dtime' should be less than 
        ! the total number of time steps
        ! ('nsteps' is not necessarily directly entered via namelist, but might be derived 
        ! from other namelist entries. This potential derivation should have happened before 
        ! the call of this subroutine. To be on the safe side, we check for its default value 
        ! 'nsteps = -999', see 'src/namelists/mo_run_nml'.)
        IF (nudging_config%idiagnose > 0 .AND. nsteps < 0) THEN
          CALL finish(routine, "nsteps seems to be not yet configured") 
        ELSEIF (nudging_config%idiagnose >= nsteps) THEN
          WRITE(message_text,'(a)') "idiagnose < nsteps = "//TRIM(int2string(nsteps))//" required"
          CALL finish(routine, message_text)
        ENDIF
      ENDIF  !nudging type

    ENDIF  !IF (nudging_config%lnudging)

    ! Indicate that crosschecks took place
    nudging_config%lchecked = .TRUE.

    ! If nudging is switch on, we checked above, if the conditions are met 
    ! that will lead to a call of 'src/configure_model/mo_nudging_config: configure_nudging', 
    ! where the switch 'nudging_config%lconfigured' is set to .TRUE. 
    ! If nudging is switched off, it is well possible that 'configure_nudging' is not called at all. 
    ! If a query for 'lconfigured' is executed somewhere in the program, 
    ! this happens typically before a query for 'lnudging'. 
    ! So we have to set 'lconfigured = .TRUE.', otherwise the program might stop inadvertently. 
    IF (.NOT. nudging_config%lnudging) nudging_config%lconfigured = .TRUE.

  END SUBROUTINE check_nudging

END MODULE mo_nudging_nml
