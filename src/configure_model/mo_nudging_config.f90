!>
!! Processing of nudging parameters for the nudging types:
!! - Upper boundary nudging
!! For the lateral boundary nudging, please see: 
!! - src/configure_model/mo_limarea_config
!! - src/configure_model/mo_interpol_config
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
MODULE mo_nudging_config

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH
  USE mo_vertical_coord_table,  ONLY: vct_a
  USE mo_util_string,           ONLY: int2string, real2string
  USE mo_run_config,            ONLY: msg_level

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: indg_type, indg_profile, ithermdyn_type
  PUBLIC :: nudging_config
  PUBLIC :: configure_nudging

  ! Module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nudging_config'

  !---------------------------------------------------------
  !                     Parameter types
  !---------------------------------------------------------

  ! Note: the following parameter types are introduced for convenience: 
  ! * combine thematically related identifiers, to unburden the USE-areas in modules, where they are required
  ! * use the number of list box entries 'ixyz%nitem' for allocation purposes and loop boundaries
  ! * most changes with regard to the identifiers can be confined to this module
  ! * avoid optical overload of 'src/shared/mo_impl_constants'
  ! For the time being, many of them contain only one identifier, 
  ! but more will follow, when the global nudging will be implemented.

  ! Their construction generally follows the rule:
  !
  ! TYPE t_ixyz
  !   INTEGER :: a      ! Identifier 1
  !   INTEGER :: b      ! Identifier 2
  !   INTEGER :: c      ! Identifier 3
  !   INTEGER :: d      ! Identifier 4
  !   ! ...
  !   !
  !   INTEGER :: nitem  ! Number of identifiers in t_ixyz
  ! END TYPE t_ixyz
  ! TYPE(t_ixyz), PARAMETER :: ixyz = t_ixyz( 1, & ! a
  !                                           2, & ! b = a + 1
  !                                           3, & ! c = b + 1
  !                                           4, & ! d = c + 1
  !
  !                                           4  ) ! nitem = d, since d is the last list element
  !
  ! Why this construction rule?
  ! -> We may want to loop over whatever the items of ixyz identify: 
  !
  !    DO jitem = 1, ixyz%nitem
  !      ! jitem successivly takes the values of ixyz%a, ixyz%b, ixyz%c, ... 
  !      ! and can be compared with namelist entries etc.
  !    ENDDO
  ! 
  ! -> The items of ixyz shall be parameters, so a setup during runtime 
  !    by means of some object is ineligible. 
  !    We have to count the items of ixyz manually and store the value in nitem. 
  !
  ! Please, be careful, if you modify the following types, 
  ! since a reasonable check for the adherence to the above construction rule 
  ! during runtime is not possible.

  ! IDENTIFIERS for nudging type:.....................................................................................
  !  
  TYPE t_indg_type
    INTEGER :: off      ! Switched off
    INTEGER :: ubn      ! Upper boundary nudging 
    !   
    INTEGER :: nitem    ! Number of preceding list elements
  END TYPE t_indg_type
  TYPE(t_indg_type), PARAMETER :: indg_type = t_indg_type( 0, & ! off
    &                                                      1, & ! ubn
    !
    &                                                      1  ) ! (not 2) nitem 

  ! IDENTIFIERS for nudging profile:..................................................................................
  !  
  TYPE t_indg_profile
    INTEGER :: sqrddist  ! Squared scaled vertical distance
    ! 
    integer :: nitem     ! Number of list elements
  END TYPE t_indg_profile
  TYPE(t_indg_profile), PARAMETER :: indg_profile = t_indg_profile( 1, & ! sqrddist
    ! 
    &                                                               1  ) ! nitem 

  ! IDENTIFIERS for type of thermodynamic variables:..................................................................
  !  
  TYPE t_ithermdyn_type
    INTEGER :: hydrostatic      ! Hydrostatically balanced
    INTEGER :: nonhydrostatic   ! Non-hydrostatic
    ! 
    INTEGER :: nitem            ! Number of list elements
  END TYPE t_ithermdyn_type
  TYPE(t_ithermdyn_type), PARAMETER ::  ithermdyn_type = t_ithermdyn_type( 1, & ! hydrostatic
    &                                                                      2, & ! nonhydrostatic
    !
    &                                                                      2  ) ! nitem 

  !--------------------------------------------------------
  !         Derived types for nudging parameters
  !--------------------------------------------------------

  ! Thresholds
  TYPE t_thr
    INTEGER :: low
    INTEGER :: med
    integer :: high
  END TYPE t_thr

  ! MAIN CONFIGURATION TYPE:..........................................................................................
  !
  TYPE t_nudging_config

    ! Namelist parameter
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
                                                     ! - 1 -> inverse squared scaled vertical distance from nudging start height
    REAL(wp) :: nudge_end_height                     ! Nudging end height
    REAL(wp) :: nudge_scale_height                   ! Scale height for profiles

    ! Derived parameters
    LOGICAL  :: lnudging                             ! Overall switch for nudging
    LOGICAL  :: ltype(indg_type%nitem)               ! And nudging-type-specific
    LOGICAL  :: lexpl_endheight                      ! .TRUE. -> explicit entry for nudge_end_height in namelist (for global nudging)
    LOGICAL  :: lexpl_scaleheight                    ! .TRUE. -> explicit entry for nudge_scale_height in namelist (for global nudging)
    INTEGER  :: ilev_start                           ! Index of grid layer corresponding to 'nudge_end_height'
    INTEGER  :: ilev_end                             ! Index of grid layer corresponding to 'nudge_start_height'
    
    ! Auxiliary variables
    TYPE(t_thr) :: msg_thr                           ! Message thresholds for print
    TYPE(t_thr) :: tmr_thr                           ! Timer thresholds
    LOGICAL     :: llatbc_type_const                 ! .TRUE. -> constant lateral boundary condition data

    ! Parameters for type status
    LOGICAL  :: lchecked    = .FALSE.                ! .TRUE. -> crosschecks of namelist took place
    LOGICAL  :: lconfigured = .FALSE.                ! .TRUE. -> configuration of this type took place
 
  END TYPE t_nudging_config

  ! Important: 'nudging_config' is uninitialized until the call of 'configure_nudging'. 
  ! So whenever you use 'nudging_config' outside this module at a relatively early point of the program sequence, 
  ! please check its status by: 
  ! 
  !  IF (nudging_config%lconfigured) THEN 
  !    ! Do what you would like to do
  !  ELSE 
  !    CALL finish(routine, "nudging_config still unconfigured.")
  !  ENDIF
  ! 
  ! to be on the safe side.
  TYPE(t_nudging_config), TARGET :: nudging_config

CONTAINS !............................................................................................................

  !>
  !! Subroutine to configure:
  !! - Upper boundary nudging
  !!
  SUBROUTINE configure_nudging(nlev)

    ! In/out variables
    INTEGER, INTENT(IN) :: nlev  ! Number of grid layers of primary domain

    ! Local variables
    REAL(wp) :: height, start_height, end_height
    INTEGER  :: jk, nlevp1
    LOGICAL  :: lfound
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: c4print
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: &
      & routine = modname//":configure_nudging"

    !----------------------------------------

    ! Crosscheck of namelist entries should have taken place beforehand 
    ! in 'src/namelists/mo_nudging_nml: check_nudging'.
    IF (.NOT. nudging_config%lchecked) THEN
      CALL finish(routine, "Crosscheck of nudging_nml still pending. &
        &Please, check the program sequence.")
    ENDIF

    !---------------------------------------------------
    ! Set auxiliary variables in 'nudging_config'
    !---------------------------------------------------
    
    ! Set thresholds for message output
    nudging_config%msg_thr%low  = 4   ! For important output
    nudging_config%msg_thr%med  = 10  ! For less important output
    nudging_config%msg_thr%high = 12  ! For least important output
    
    ! Set thresholds for timers
    nudging_config%tmr_thr%low  = 1
    nudging_config%tmr_thr%med  = 5
    nudging_config%tmr_thr%high = 10

    !---------------------------------------------------
    ! Initialize derived parameters in 'nudging_config'
    !---------------------------------------------------

    nudging_config%ltype(:)   = .FALSE.
    nudging_config%ilev_start = -1
    nudging_config%ilev_end   = -1

    !---------------------------------------------------
    ! Determine derived parameters
    !---------------------------------------------------

    IF (nudging_config%lnudging) THEN

      ! Set switch for selected nudging type
      nudging_config%ltype(nudging_config%nudge_type) = .TRUE.

      ! Evaluate nudging start and end heights:
      ! Check status of 'vct_a'
      IF (.NOT. ALLOCATED(vct_a)) CALL finish(TRIM(routine), 'vct_a still uninitialized.')
      ! Auxiliary variables:
      ! Height of model top
      nlevp1    = nlev + 1
      lfound    = .FALSE.
      ! 1) Start height:
      start_height = nudging_config%nudge_start_height
      DO jk = nlevp1, 1, -1
        ! (No shift of vertical index necessary, because this is done for the primary domain, only)
        ! Nominal height of half level
        ! (Note: we (arbitrarily) defined that those vertical grid layers, 
        ! in which 'start/end_height' lie, should be subject to nudging. 
        ! This is, why we do not use 'height = 0.5_wp * (vct_a(jk) + vct_a(jk+1))' here.)
        height = vct_a(jk)
        IF (height > start_height) THEN
          ! We have found the lowermost grid layer, where nudging is applied
          ! (Because vertical looping starts from model top, 'ilev_end' corresponds 
          ! to the start height, and 'ilev_start' corresponds to the end height)
          nudging_config%ilev_end = jk
          ! Indicate the find
          lfound = .TRUE.
          EXIT 
        ENDIF
      ENDDO  !jk
      IF (.NOT. lfound) CALL finish(TRIM(routine), 'Could not find ilev_end corresponding to nudge_start_height.')
      ! Just to make sure
      nudging_config%ilev_end = MIN( MAX( 1, nudging_config%ilev_end ), nlev )
      ! 2) End height
      end_height = nudging_config%nudge_end_height
      lfound     = .FALSE.
      DO jk = 1, nudging_config%ilev_end + 1
        height = vct_a(jk)
        IF (height < end_height) THEN
          ! We have found the uppermost grid layer, where nudging is applied
          nudging_config%ilev_start = jk - 1
          ! Indicate the find
          lfound = .TRUE.
          EXIT
        ENDIF
      ENDDO
      IF (.NOT. lfound) CALL finish(TRIM(routine), 'Could not find ilev_start corresponding to nudge_end_height.')
      ! Just to make sure
      nudging_config%ilev_start = MIN( MAX( 1, nudging_config%ilev_start ), nlev )
      ! We replace the namelist entry for end height by the effective end height
      nudging_config%nudge_end_height   = 0.5_wp * ( vct_a(nudging_config%ilev_start)     &
        &                                          + vct_a(nudging_config%ilev_start + 1) )
      ! Finally, we can compute the nudging scale height 
      ! (if there is no explicit namelist input for it, or 'nudge_profile = 1' is selected)
      IF ((.NOT. nudging_config%lexpl_scaleheight) .OR. nudging_config%nudge_profile == indg_profile%sqrddist) THEN
        nudging_config%nudge_scale_height = nudging_config%nudge_end_height - nudging_config%nudge_start_height
      ENDIF

      !---------------------------------------------------
      ! Print some info
      !---------------------------------------------------

      IF (msg_level >= nudging_config%msg_thr%med) THEN 

        ! Nudging type
        SELECT CASE(nudging_config%nudge_type)
        CASE(indg_type%ubn)
          c4print = "upper boundary nudging for limited-area mode"
        END SELECT
        WRITE(message_text,'(a)') 'Nudging is switched on.'
        CALL message(TRIM(routine), message_text)
        WRITE(message_text,'(a)') 'Selected nudging type: '//TRIM(c4print)
        CALL message(TRIM(routine), message_text)

        ! Profile of nudging strength
        SELECT CASE(nudging_config%nudge_profile)
        CASE(indg_profile%sqrddist)
          c4print = "inverse squared scaled vertical distance from nudging start height"
        END SELECT
        WRITE(message_text,'(a)') 'Selected profile of nudging strength between start and end height: ' &
          & //TRIM(c4print)
        CALL message(TRIM(routine), message_text)

        ! Nudging start and end heights and the scale height
        WRITE(message_text,'(a)') 'Nudging includes the full levels between:'
        CALL message(TRIM(routine), message_text)
        jk     = nudging_config%ilev_start
        height = 0.5_wp * ( vct_a(jk) + vct_a(jk+1) )  
        WRITE(message_text,'(a)') ' - Level jk = '//TRIM(int2string(jk))//', height z(jk) = ' &
          & //TRIM(real2string(height))//' m'
        CALL message(' ', message_text, adjust_right=.TRUE.)
        jk     = nudging_config%ilev_end
        height = 0.5_wp * ( vct_a(jk) + vct_a(jk+1) )  
        WRITE(message_text,'(a)') ' - Level jk = '//TRIM(int2string(jk))//', height z(jk) = ' &
          & //TRIM(real2string(height))//' m'
        CALL message(' ', message_text, adjust_right=.TRUE.)
        WRITE(message_text,'(a)') 'Scale height for profile of nudging strength: ' &
          & //TRIM(real2string(nudging_config%nudge_scale_height))//' m'
        CALL message(TRIM(routine), message_text)
        
      ENDIF  !IF (msg_level >= nudging_config%msg_thr%med)
      
    ENDIF  !IF (nudging_config%lnudging)
    
    ! Indicate that the setup of 'nudging_config' took place
    nudging_config%lconfigured = .TRUE.

  END SUBROUTINE configure_nudging

END MODULE mo_nudging_config
