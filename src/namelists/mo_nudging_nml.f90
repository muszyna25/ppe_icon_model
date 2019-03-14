!>
!! Processing of nudging parameters for the nudging types:
!! - Upper boundary nudging
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
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_nudging_config,      ONLY: indg_type, indg_profile, ithermdyn_type, &
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
    INTEGER  :: jitem
    LOGICAL  :: lerror, lexpl_endheight, lexpl_scaleheight
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: &
      & routine = modname//":read_nudging_namelist"

    ! Namelist variables:.............................................................................................
    INTEGER  :: nudge_type                           ! Nudging type identifier:
                                                     ! - 0 -> switched off
                                                     ! - 1 -> upper boundary nudging
    REAL(wp) :: max_nudge_coeff_vn                   ! Max. nudging coefficient for horizontal wind
    REAL(wp) :: max_nudge_coeff_thermdyn             ! Max. nudging coefficient for thermodynamic variables
    REAL(wp) :: nudge_start_height                   ! Nuding start height

    ! Non-namelist parameters 
    ! in case of upper boundary nudging
    INTEGER  :: thermdyn_type                        ! Type of thermodynamic variables (from the perspective
                                                     ! of the nudging-target variables in ICON):
                                                     ! - 1 -> hydrostatic: hydrostatic pressure & temperature are nudged
                                                     !        (diagnostic variables)
                                                     ! - 2 -> nonhydrostatic: density & virtual potential temperature are nudged
                                                     !        (prognostic variables)
    INTEGER  :: nudge_profile                        ! Profile of nudging strength:
                                                     ! - 1 -> squared scaled vertical distance from nudging start height
    REAL(wp) :: nudge_end_height                     ! Nudging end height
    REAL(wp) :: nudge_scale_height                   ! Scale height for profiles
    !.................................................................................................................

    
    NAMELIST /nudging_nml/ nudge_type, max_nudge_coeff_vn, max_nudge_coeff_thermdyn,  &
      &                    nudge_start_height

    !----------------------------------------

    !------------------------------------------------------------------------
    ! Default settings
    !------------------------------------------------------------------------

    nudge_type                = indg_type%off                 ! 0 -> upper boundary nudging switched off
    max_nudge_coeff_vn        = 0.04_wp                       ! Max. nudging coefficient for horizontal wind
    max_nudge_coeff_thermdyn  = 0.075_wp                      ! Max. nudging coefficient for thermodynamic variables
    nudge_start_height        = 12000._wp                     ! (m) Nudging start height

    ! Non-namelist parameters 
    ! in case of upper boundary nudging
    thermdyn_type             = ithermdyn_type%hydrostatic    ! 1 -> diagnostic variables: hydrostatic pressure 
                                                              !      and temperature are nudged
    nudge_profile             = indg_profile%sqrddist         ! 1 -> inverse squared scaled vertical distance from start height
    nudge_end_height          = -1._wp                        ! (m) Nudging end height
                                                              ! negative value means: 'nudge_end_height = top_height' 
                                                              ! ('sleve_nml: top_height' -> height of model top)
    nudge_scale_height        = -1._wp                        ! (m) Scale height for nudging strength profile,
                                                              ! negative value means: 'nudge_scale_height = nudge_end_height - nudge_start_height'

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

    !------------------------------------------------------------------------
    ! Consistency checks
    !------------------------------------------------------------------------

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
      ! Check 'max_nudge_coeff_...'
      IF (max_nudge_coeff_vn < 0._wp) THEN
        CALL finish(TRIM(routine), 'max_nudge_coeff_vn >= 0 required.')
      ELSEIF(max_nudge_coeff_thermdyn < 0._wp) THEN
        CALL finish(TRIM(routine), 'max_nudge_coeff_thermdyn >= 0 required.')
      ENDIF
      ! Check 'nudge_start_height': 
      IF (nudge_start_height < 0._wp) THEN
        CALL finish(TRIM(routine), 'nudge_start_height >= 0 required.') 
      ENDIF
      ! "Check" parameters, which are non-namelist parameters 
      ! in case of upper boundary nudging:
      ! Check 'thermdyn_type'
      lerror = .TRUE.
      DO jitem = 1, ithermdyn_type%nitem
        IF (thermdyn_type == jitem) THEN
          lerror = .FALSE.
          EXIT
        ENDIF
      ENDDO  !jitem
      IF (lerror) CALL finish(TRIM(routine), 'Invalid thermdyn_type: '//TRIM(int2string(thermdyn_type))) 
      ! Check 'nudge_profile'
      lerror = .TRUE.
      DO jitem = 1, indg_profile%nitem
        IF (nudge_profile == jitem) THEN
          lerror = .FALSE.
          EXIT
        ENDIF
      ENDDO  !jitem
      IF (lerror) CALL finish(TRIM(routine), 'Invalid nudge_profile: '//TRIM(int2string(nudge_profile))) 
      ! Check 'nudge_end_height'
      ! If the nudging end height would be a namelist variable, 
      ! there would be an explicit entry for it in the namelist input, 
      ! if 'nudge_end_height > 0' (compare default above)
      IF (nudge_end_height > 0._wp) THEN
        lexpl_endheight = .TRUE. 
      ELSE
        ! This case means: 'nudge_end_height = top_height' (see 'check_nudging' below)
        lexpl_endheight = .FALSE.
      ENDIF
      ! Likewise the nudging scale height, if it would be a namelist variable
      IF (nudge_scale_height > 0._wp) THEN
        lexpl_scaleheight = .TRUE.
      ELSE
        ! This case means: 'nudge_scale_height = nudge_end_height - nudge_start_height' 
        ! (see 'check_nudging' below and 'src/configure_model/mo_nudging_config: configure_nudging')
        lexpl_scaleheight = .FALSE.
      ENDIF
      IF(lexpl_endheight .AND. (nudge_start_height > nudge_end_height)) THEN
        CALL finish(TRIM(routine), 'nudge_start_height < nudge_end_height required.') 
      ENDIF

    ELSE

        lexpl_endheight = .FALSE.
        lexpl_scaleheight = .FALSE.
    ENDIF  !IF (nudge_type /= indg_type%off)

    !-----------------------------------------------------
    ! Store the namelist in the nudging config type
    !-----------------------------------------------------

    nudging_config%nudge_type               = nudge_type
    nudging_config%max_nudge_coeff_vn       = max_nudge_coeff_vn
    nudging_config%max_nudge_coeff_thermdyn = max_nudge_coeff_thermdyn
    nudging_config%nudge_start_height       = nudge_start_height
    ! Non-namelist parameters 
    ! in case of upper boundary nudging
    nudging_config%thermdyn_type            = thermdyn_type
    nudging_config%nudge_profile            = nudge_profile
    nudging_config%nudge_end_height         = nudge_end_height
    nudging_config%nudge_scale_height       = nudge_scale_height
    nudging_config%lexpl_endheight          = lexpl_endheight
    nudging_config%lexpl_scaleheight        = lexpl_scaleheight

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
  SUBROUTINE check_nudging( n_dom, iequations, iforcing, ivctype, top_height, &
    &                       l_limited_area, lsparse_latbc, itype_latbc, nudge_hydro_pres, &
    &                       LATBC_TYPE_CONST, is_plane_torus, lart, ndyn_substeps         )

    ! In/out variables
    INTEGER,  INTENT(IN) :: n_dom             ! Number of model domains, 1=primary domain only 
    INTEGER,  INTENT(IN) :: iequations        ! Type of model equations
    INTEGER,  INTENT(IN) :: iforcing          ! Physics package
    INTEGER,  INTENT(IN) :: ivctype           ! Type of vertical grid
    REAL(wp), INTENT(IN) :: top_height        ! Height of model top
    LOGICAL,  INTENT(IN) :: l_limited_area    ! .true.: limited-area mode is switched on
    LOGICAL,  INTENT(IN) :: lsparse_latbc     ! .true.: only nudging data for boundary rows are read in
    INTEGER,  INTENT(IN) :: itype_latbc       ! Type of limited area boundary nudging
    LOGICAL,  INTENT(IN) :: nudge_hydro_pres  ! .true.: use hydrostatic pressure for lateral boundary nudging
    INTEGER,  INTENT(IN) :: LATBC_TYPE_CONST  ! Identifier for itype_latbc = constant nudging target data
    LOGICAL,  INTENT(IN) :: is_plane_torus    ! .true.: torus mode
    LOGICAL,  INTENT(IN) :: lart              ! .true.: ART interface cut in
    INTEGER,  INTENT(IN) :: ndyn_substeps     ! Number of dynamics substeps per fast-physics step

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
    nudging_config%lnudging = ( nudging_config%nudge_type == indg_type%ubn .AND. &
      &                         l_limited_area                                   )

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

      ! Mixed criteria 
      ! + crosscheck-induced harmonizations:

      IF ((.NOT. nudging_config%lexpl_endheight) .OR. (nudging_config%nudge_end_height > top_height)) THEN
        ! We (re)set the nudging end height to the height of the model top
        nudging_config%nudge_end_height = top_height
        IF (.NOT. (nudging_config%nudge_type == indg_type%ubn)) THEN
          ! This should not be written in case of upper boundary nudging, 
          ! since it is the default behaviour for this nudging type
          WRITE(message_text,'(a)') 'Nudging end height (re)set to nudge_end_height = top_height = ' &
            & //TRIM(real2string(top_height))//' m'
          CALL message(TRIM(routine), message_text)
        ENDIF
        ! Now, we have to check, if 'nudge_start_height < nudge_end_height' still holds
        IF (nudging_config%nudge_start_height > nudging_config%nudge_end_height) THEN
          IF ((nudging_config%nudge_type == indg_type%ubn)) THEN
            ! Message in case of upper boundary nudging
            WRITE(message_text,'(2a)') 'nudge_start_height > top_height! ', &
              & 'Please, choose a value according to: nudge_start_height < top_height = ' &
              & //TRIM(real2string(top_height))//' m'
          ELSE
            ! Message for all other cases
            WRITE(message_text,'(2a)') 'After (re)set of nudge_end_height, nudge_start_height > nudge_end_height! ',&
              & 'Please, choose values according to: nudge_start_height < nudge_end_height <= top_height = ' &
              & //TRIM(real2string(top_height))//' m'
          ENDIF
          CALL finish(routine, message_text)
        ENDIF
        ! The nudging scale height will be computed in 'src/configure_model/mo_nudging_config: configure_nudging'
      ENDIF  !IF ((.NOT. nudging_config%lexpl_endheight) .OR. (nudging_config%nudge_end_height > top_height))

      ! There are two options for the thermodynamic nudging variables:
      ! (1) hydrostatic pressure and temperature 
      ! (2) density and virtual potential temperature
      ! In case of upper boundary nudging, we should be consistent with the choice of 'latbc_config%nudge_hydro_pres'
      ! (see, e.g., the processes controlled by 'latbc_config%nudge_hydro_pres' in 'src/atm_dyn_iconam/mo_nh_stepping')
      IF ((nudging_config%nudge_type == indg_type%ubn) .AND. nudge_hydro_pres) THEN
        ! Option (1)
        nudging_config%thermdyn_type = ithermdyn_type%hydrostatic
      ELSEIF ((nudging_config%nudge_type == indg_type%ubn) .AND. (.NOT. nudge_hydro_pres)) THEN
        ! Option (2)
        nudging_config%thermdyn_type = ithermdyn_type%nonhydrostatic
      ENDIF
      ! The option to nudge the thermodynamic variables: hydrostatic pressure & temperature 
      ! is not available, if the nudging target data are constant in time
      IF ( (nudging_config%nudge_type == indg_type%ubn) .AND. &
        &  nudging_config%llatbc_type_const             .AND. &
        &  nudging_config%thermdyn_type == ithermdyn_type%hydrostatic ) THEN
        ! If upper boundary nudging is switched on, we have to abort, 
        ! since we should not change 'latbc_config%nudge_hydro_pres' at this place
        CALL finish(routine, "Nudging target data are constant in time. limarea_nml: nudge_hydro_pres = .false. required")
      ELSEIF ( nudging_config%llatbc_type_const .AND. &
        &      nudging_config%thermdyn_type == ithermdyn_type%hydrostatic ) THEN
        ! Fallback option: nudge density & virtual potential temperature
        nudging_config%thermdyn_type = ithermdyn_type%nonhydrostatic
        ! In this case a warning should suffice
        WRITE(message_text,'(a)') 'WARNING, the nudging target data are constant in time. &
          &Reset thermdyn_type to '//TRIM(int2string(ithermdyn_type%nonhydrostatic))
        CALL message(TRIM(routine), message_text)
      ENDIF

    ENDIF  !IF (nudging_config%lnudging)

    ! Indicate that crosschecks took place
    nudging_config%lchecked = .TRUE.

    ! If nudging is switch on, we checked above, if the conditions are met 
    ! that will lead to a call of 'src/configure_model/mo_nudging_config: configure_nudging', 
    ! where the switch 'nudging_config%lconfigured' is set to .true. 
    ! If nudging is switched off, it is well possible that 'configure_nudging' is not called at all. 
    ! If a query for 'lconfigured' is executed somewhere in the program, 
    ! this happens typically before a query for 'lnudging'. 
    ! So we have to set 'lconfigured = .true.', otherwise the program might stop inadvertently. 
    IF (.NOT. nudging_config%lnudging) nudging_config%lconfigured = .TRUE.

  END SUBROUTINE check_nudging

END MODULE mo_nudging_nml
