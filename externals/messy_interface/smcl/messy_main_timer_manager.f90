MODULE messy_main_timer_manager

  ! MESSy
  USE messy_main_constants_mem, ONLY: dp, OneDay              &
                                    , STRLEN_MEDIUM, STRLEN_ULONG
  USE messy_main_timer

  IMPLICIT NONE
  PRIVATE
  SAVE

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! TIME MANAGER 
!  adopted from mo_time_manager of ECHAM5
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  TYPE time_manager
     PRIVATE     ! no external access possible
     CHARACTER(STRLEN_MEDIUM) :: label = ''   ! short description of the manager
     TYPE (time_days)         :: start_date   ! initial date at pos=0
     LOGICAL  :: init         = .FALSE. ! manager state for access control
     LOGICAL  :: freeze       = .FALSE. ! lock/unlock the manager
     INTEGER  :: pos          = 0       ! present step of the manager
     REAL(dp) :: delta_time   = 0._dp   ! distance between two points in seconds
  END TYPE time_manager
  PUBLIC :: time_manager

  ! *************** INTERFACE SUBROUTINES AND FUNCTIONS
  !
  ! input: label, startdate, increment
  !                different formats
  !        only a  date will reset  the start date
  !        only an increment resets the increment and the start date
  !                           
  PUBLIC :: manager_init    !   set/reset manager state
  !   PARAMETERS(
  !     manager [time_manager] (inp, out) the time axis manager
  !     label   [character]    (inp)      name of the manager
  !     date    [time_days|time_intern|time_native] 
  !                            (inp)      start date
  !     second  [integer]      (inp)      distance between two points
  !     zsecond [real]         (inp)      distance between two points
  !     offset  [integer]      (inp)      point offset moving position
  !     lfreeze [logical]      (inp)      lock/unlock feature
  !   )
  !
  !   separation of functions using the type and number of parameter
  !
  INTERFACE manager_init
    !   manager, label, date, second  -> initialisation
    MODULE PROCEDURE manager_init_days
    !   manager, date                 -> reset start date, adjust position
    MODULE PROCEDURE manager_reinit_days
    !   manager, zsecond              -> reset increment, adjust start date
    MODULE PROCEDURE manager_reinit_incr
    !   manager, offset               -> set new position using offset
    MODULE PROCEDURE manager_step
    !   manager, lfreeze              -> lock/unlock time manager
    MODULE PROCEDURE manager_freeze
  END INTERFACE

  PUBLIC :: timer_rewind_manager
  PUBLIC :: manager_print
  PUBLIC :: manager_state
  !
  ! manager_state [subroutine, interface]
  !   get informations about the manager state
  !   (
  !   manager  [time_manager] (inp)      the time axis manager
  !   date     [time_days]    (inp, out) calendar date
  !   position [integer]      (out)      position on the manager
  !   lequal   [logical]      (out)      date fit to manager position
  !   )
  !   manager, date, position         -> get date for position at manager
  !   manager, date                   -> get present date of manager
  !   manager, date, position, lequal -> get nearest position for date
  !   manager, dtime                  -> get delta time
  !   manager, pos                    -> get time step
  !
  INTERFACE manager_state
    MODULE PROCEDURE manager_any_date
    MODULE PROCEDURE manager_current_date
    MODULE PROCEDURE manager_step_at_date
    MODULE PROCEDURE manager_dtime
    MODULE PROCEDURE manager_pos
  END INTERFACE

CONTAINS

  SUBROUTINE manager_init_days (manager, name, date, incr, message)
    
    IMPLICIT NONE
    INTRINSIC :: LEN, MIN, PRESENT, REAL, TRIM
    !
    ! initialisation - reinitialisation of a time manager
    !
    TYPE (time_manager), INTENT(inout) :: manager
    CHARACTER(len=*),    INTENT(in)    :: name
    TYPE (time_days),    INTENT(in)    :: date
    INTEGER,             INTENT(in)    :: incr
    CHARACTER(LEN=STRLEN_ULONG), OPTIONAL, INTENT(OUT) :: message

    ! LOCAL
    INTEGER    :: is
    CHARACTER(LEN=STRLEN_ULONG) :: message_text

    IF (manager%freeze) THEN
      WRITE(message_text,*) 'Time manager <',&
           TRIM(manager%label),'> is frozen, skip initialization.'

    ELSE IF (manager%init) THEN
      WRITE(message_text,*) 'Time manager <',&
           TRIM(manager%label), '> was initialzed before.'

    ELSE
      DO is=1,STRLEN_MEDIUM
        manager%label(is:is) = ''
      END DO
      is = MIN(LEN(TRIM(name)),STRLEN_MEDIUM)
      manager%label(1:is)  = name(1:is)
      manager%start_date   = date
      manager%pos          = 0
      manager%delta_time   = REAL(incr,dp)
      manager%init         = .TRUE.
      WRITE(message_text,*) &
           'Basic initialization of time manager <',&
           TRIM(manager%label),'> done.'

    END IF
    IF (PRESENT(message)) THEN
       message = message_text
    ELSE
       WRITE (*,*) message_text
    ENDIF

  END SUBROUTINE manager_init_days

  !------------------------------------------------------------
  ! reset the manager start date and correct the present time step
  ! the manager position in real time should not changed for consistency


  SUBROUTINE manager_reinit_days (manager, new_date, ierr, message)

    IMPLICIT NONE
    INTRINSIC :: INT, PRESENT, REAL , TRIM

    ! reset start date
    TYPE (time_manager), INTENT(inout) :: manager
    TYPE (time_days),    INTENT(in)    :: new_date
    INTEGER,             INTENT(OUT)   :: ierr
    CHARACTER(LEN=STRLEN_ULONG), OPTIONAL, INTENT(OUT) :: message

    ! LOCAL
    TYPE (time_days) :: current_date
    REAL(dp)         :: zsecs
    INTEGER          :: iday1, isec1, iday2, isec2, new_step
    CHARACTER(LEN=STRLEN_ULONG) :: message_text

    message_text= ' '

    IF (.NOT. manager%init) THEN
       ierr=3436
       RETURN

    ELSE IF (manager%freeze) THEN
       WRITE(message_text,*) 'Time manager <',&
            TRIM(manager%label),'> is frozen, reinit skiped.'

    ELSE

      ! get current date of manager
      CALL manager_current_date (manager, current_date, ierr)
      IF (ierr /= 0) RETURN

      ! get time since new start date
      CALL date_get (current_date, iday1, isec1, ierr)
      IF (ierr /=0 ) RETURN
      CALL date_get (new_date,     iday2, isec2, ierr)
      IF (ierr /=0 ) RETURN
      zsecs = OneDay*REAL((iday1-iday2),dp) + REAL((isec1-isec2),dp)

      ! get new manager step
      new_step = INT((zsecs+0.0001_dp)/manager%delta_time)

      IF (new_step < 0) &
           WRITE(message_text,*) 'Warning: New time manager position < 0'

      manager %pos        = new_step
      manager %start_date = new_date          ! correct the start date

   END IF
    IF (PRESENT(message)) THEN
       message = message_text
    ELSE
       write (*,*) message_text
    ENDIF

    ierr = 0
    
  END SUBROUTINE manager_reinit_days

  !------------------------------------------------------------

  SUBROUTINE manager_reinit_incr (manager, incr, ierr, message)

    IMPLICIT NONE
    INTRINSIC :: INT, PRESENT, REAL,  TRIM
    ! reset increment between two steps
    TYPE (time_manager), INTENT(inout) :: manager
    REAL(dp),            INTENT(in)    :: incr
    INTEGER,             INTENT(OUT)   :: ierr
    CHARACTER(LEN=STRLEN_ULONG), OPTIONAL, INTENT(OUT) :: message

    ! LOCAL
    REAL(dp)         :: zsecs
    TYPE (time_days) :: current_date
    INTEGER          :: idays, isecs
    CHARACTER(LEN=STRLEN_ULONG) :: message_text

    ierr = 0
    message_text = ' '
    ! find a new start date
    ! the position of the manager will be conserved

    IF (manager%freeze) THEN
      WRITE(message_text,*) 'Time manager <',&
           TRIM(manager%label),'> is frozen, increment not changed.'

    ELSE

      ! distance from present point backward with new incr
      zsecs = REAL(manager%pos ,dp)* incr
      idays = INT((zsecs+0.0001_dp)/OneDay)
      isecs = INT(zsecs - REAL(idays,dp)*OneDay+0.0001_dp)
      idays = -idays
      isecs = -isecs

      ! get the present date
      CALL manager_current_date (manager, current_date, ierr)
      IF (ierr /= 0) RETURN

      ! calculates the new start date
      CALL add_date (idays, isecs, current_date, ierr)
      IF (ierr /=0 ) RETURN

      CALL copy_date(current_date,manager %start_date, ierr)
      IF (ierr /=0 ) RETURN
      
      manager %delta_time = incr     ! set the new increment

    END IF

    IF (PRESENT(message)) THEN
       message = message_text
    ELSE
       write (*,*) message_text
    ENDIF
    
  END SUBROUTINE manager_reinit_incr

  !------------------------------------------------------------
  ! reset position

  SUBROUTINE manager_step (manager, offset, ierr)
    TYPE (time_manager), INTENT(inout) :: manager
    INTEGER,             INTENT(in)    :: offset
    INTEGER,             INTENT(OUT)   :: ierr

    INTEGER :: new_step

    ierr = 0

    new_step = manager%pos + offset
    IF (new_step < 0) THEN
       ierr = 3437
       RETURN
    ENDIF

    manager %pos = new_step

  END SUBROUTINE manager_step

  !------------------------------------------------------------

  SUBROUTINE manager_freeze (manager, lfreeze, ierr, message)

    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM

    TYPE (time_manager), INTENT(inout) :: manager
    LOGICAL,             INTENT(in)    :: lfreeze
    INTEGER,             INTENT(OUT)   :: ierr
    CHARACTER(LEN=STRLEN_ULONG), OPTIONAL, INTENT(OUT) :: message

    ! LOCAL
    CHARACTER(LEN=STRLEN_ULONG) :: message_text = ''

    IF (.NOT. manager%init) THEN
       ierr = 3400
       IF (PRESENT(message)) &
       WRITE(message,*) &
           'Can not freeze state of time manager <',&
           TRIM(manager%label),&
           '>, missing initialization before.'
       RETURN
    ELSE IF (manager%freeze) THEN
      WRITE(message_text,*) 'Warning Time manager <',&
           TRIM(manager%label),'> was frozen before.'
    ELSE IF(lfreeze) THEN
      manager %freeze = .TRUE.
      WRITE(message_text,*) &
           'Initial state of time manager <',&
           TRIM(manager%label),'> is locked now!'
    END IF

    IF (.NOT.lfreeze) THEN
      manager %freeze = .FALSE.
      WRITE(message_text,*) TRIM(message_text),' The time manager <',&
           TRIM(manager%label),'> is unlocked now!'

    END IF

    IF (PRESENT(message)) THEN
       message = message_text
    ELSE
       write (*,*) 'messy_main_timer: freeze:', message_text
    ENDIF

    ierr = 0
    
  END SUBROUTINE manager_freeze

  !------------------------------------------------------------

  SUBROUTINE manager_print (manager, ierr, short, lwrite)

    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM
    ! print out manager settings
    ! call for p_io only
    TYPE (time_manager), INTENT(in) :: manager
    INTEGER,             INTENT(OUT):: ierr
    LOGICAL, OPTIONAL,   INTENT(in) :: short
    LOGICAL, OPTIONAL,   INTENT(in) :: lwrite
    
    TYPE (time_days)   :: current_date
    CHARACTER(len=256) :: my_message
    
    LOGICAL :: my_short

    IF (PRESENT(lwrite)) THEN
       IF (.NOT. lwrite) RETURN
    ENDIF

    my_short = .FALSE. ; IF (PRESENT(short)) my_short = short

    IF (manager%init) THEN

       IF (my_short) THEN         ! short message
          WRITE(*,*) 'State of >>',TRIM(manager%label),'<<'
          WRITE(*,*) 'Step: ',manager%pos,' dtime [s]: ',manager%delta_time
       ELSE                       ! long message
          WRITE(*,*) 'State of manager >>',TRIM(manager%label),'<<'
          WRITE(*,*) 'Step counter : ', manager%pos
          WRITE(*,*) 'Time step [s]: ', manager%delta_time
       END IF

       CALL print_date_components(manager%start_date,ierr,mess=my_message)
       IF (ierr /= 0) RETURN
       WRITE(*,*) 'Initial date: ',TRIM(my_message)
       
       CALL manager_current_date (manager, current_date, ierr)
       IF  (ierr /= 0) RETURN
       CALL print_date_components (current_date,ierr,mess=my_message)
       IF (ierr /=0 ) RETURN
       WRITE(*,*) 'Current date: ',TRIM(my_message)

    ELSE

       WRITE(*,*) 'Warning: Time manager >>',TRIM(manager%label), &
            '<< was not initialized, can not print settings.'
    END IF
    WRITE (*,*) ' '

    ierr = 0
    
  END SUBROUTINE manager_print


  !------------------------------------------------------------
  ! get date from manager

  SUBROUTINE manager_any_date (manager, date, point, ierr)

    IMPLICIT NONE
    INTRINSIC :: INT, REAL
    ! calculates date at a given timestep of the manager
    TYPE (time_manager), INTENT(in)  :: manager
    TYPE (time_days),    INTENT(out) :: date
    INTEGER,             INTENT(in)  :: point
    INTEGER,             INTENT(OUT) :: ierr

    ! LOCAL
    REAL(dp)     :: zsecs
    INTEGER      :: idays, isecs

    ierr = 0

    date  = manager %start_date      ! load start date

    ! calculates new offset
    IF (point < 0) THEN
      zsecs = REAL(point,dp)*manager%delta_time - 0.0001_dp
    ELSE
      zsecs = REAL(point,dp)*manager%delta_time + 0.0001_dp
    END IF
    idays = INT(zsecs/OneDay)

    IF (zsecs < 0.0_dp) THEN
      isecs = INT(zsecs - REAL(idays,dp)*OneDay - 0.0001_dp)
    ELSE
      isecs = INT(zsecs - REAL(idays,dp)*OneDay + 0.0001_dp)
    END IF

    ! final calculation of date at point
    CALL add_date (idays, isecs, date,  ierr)
    IF (ierr /=0 ) RETURN

  END SUBROUTINE manager_any_date

  ! -----------------------------------------------------------

  SUBROUTINE manager_current_date (manager, date, ierr)

    TYPE (time_manager), INTENT(in)  :: manager
    TYPE (time_days),    INTENT(out) :: date
    INTEGER,             INTENT(OUT) :: ierr

    CALL manager_any_date (manager, date, manager%pos, ierr)

  END SUBROUTINE manager_current_date


  !------------------------------------------------------------

  SUBROUTINE manager_step_at_date (manager, current_date, steps, lfit, ierr)

    IMPLICIT NONE
    INTRINSIC :: INT, MOD, REAL

    ! get nearest manager step for special date
    ! get the time step (pointer) for a known date/time
    TYPE (time_manager), INTENT(in)  :: manager
    TYPE (time_days),    INTENT(in)  :: current_date
    INTEGER,             INTENT(out) :: steps
    LOGICAL,             INTENT(out) :: lfit
    INTEGER,             INTENT(OUT) :: ierr

    ! LOCAL
    TYPE (time_days) :: check_date
    INTEGER          :: iday1, isec1, iday2, isec2, iday, isecond
    REAL(dp)         :: zsecs

    ierr = 0
    ! calculate time difference in seconds
    CALL date_get (current_date,       iday1, isec1,ierr)
    IF (ierr /=0 ) RETURN 
    CALL date_get (manager%start_date, iday2, isec2,ierr)
    IF (ierr /=0 ) RETURN 
    ! calculates the difference between the dates in seconds
    zsecs = REAL(iday1-iday2,dp) * OneDay +  REAL(isec1-isec2,dp)
    IF (zsecs < 0.0_dp) THEN
      zsecs = zsecs-0.001_dp
    ELSE
      zsecs = zsecs+0.001_dp
    END IF
    steps = INT(zsecs/manager%delta_time)

    ! recalculate the present day again using the number of steps found
    zsecs = manager%delta_time*REAL(steps,dp)
    IF (zsecs < 0.0_dp) THEN
      zsecs   = zsecs - 0.001_dp
    ELSE
      zsecs   = zsecs + 0.001_dp
    END IF
    iday    = INT(zsecs/OneDay)
    isecond = MOD(INT(zsecs),INT(OneDay))

    CALL copy_date(manager %start_date, check_date, ierr)
    IF (ierr /= 0) RETURN
    CALL add_date (iday, isecond, check_date,ierr)
    IF (ierr /= 0) RETURN

    CALL if_equal(check_date, current_date, lfit, ierr)
    IF (ierr /= 0) RETURN

    IF (.NOT. lfit) THEN
       ierr =1
       write(*,*) 'manager_step_at_date: date does not fit into time grid'
    END IF
  END SUBROUTINE manager_step_at_date

! ------------------------------------------------------------------------
! MANAGER INTERFACE:
! ------------------------------------------------------------------------
  SUBROUTINE manager_dtime (manager,delta_time, ierr)
    ! return the manager increment (time step)

    TYPE (time_manager), INTENT(in)  :: manager
    REAL(dp),            INTENT(out) :: delta_time
    INTEGER,             INTENT(OUT) :: ierr

    IF (manager%init) THEN
      delta_time = manager%delta_time
      ierr = 0
    ELSE
      CALL manager_print(manager,ierr)
      ierr = 3434
      RETURN
    END IF

  END SUBROUTINE manager_dtime

! ------------------------------------------------------------------------

  SUBROUTINE manager_pos (manager,position, ierr)

    ! return the manager position
    TYPE (time_manager), INTENT(in)  :: manager
    INTEGER            , INTENT(out) :: position
    INTEGER            , INTENT(OUT) :: ierr

    IF (manager%init) THEN
      position = manager %pos
      ierr = 0
    ELSE
      CALL manager_print(manager, ierr)
      ierr = 3434
      RETURN
    END IF

  END SUBROUTINE manager_pos

! ------------------------------------------------------------------------

  SUBROUTINE timer_rewind_manager (manager, new_date, fix_start_date, ierr &
       , message)

    ! reset the present date of the time manager

    IMPLICIT NONE
    INTRINSIC :: TRIM

    TYPE (time_manager), INTENT(inout) :: manager    ! the time manager
    TYPE (time_days),    INTENT(in)    :: new_date   ! set position to that date
    LOGICAL,             INTENT(in)    :: fix_start_date  ! T changes time step
    INTEGER, OPTIONAL,   INTENT(OUT)   :: ierr            ! error status
    CHARACTER(LEN=STRLEN_ULONG), OPTIONAL, INTENT(OUT) :: &
         message

    ! LOCAL
    INTEGER          :: isecs, idays, new_step, old_step, offset
    TYPE (time_days) :: wrk_day
    LOGICAL          :: lfit, learlier
    INTEGER          :: ierror

    ierror = 0
    IF (PRESENT(message)) message = ''

    IF (manager%freeze ) THEN
       ierror =3400
       IF (PRESENT(message)) &
       WRITE(message,*) 'Operation not allowed on time manager',&
            TRIM(manager%label),': unfreeze first!'
       RETURN
    ELSE IF (fix_start_date) THEN
       ! *********** reset manager without changing the start date
       CALL if_less(new_date,manager%start_date,learlier,ierror)
       IF (ierror /= 0) THEN
          IF (PRESENT(ierr)) ierr=ierror
          RETURN
       ENDIF
       IF (learlier) THEN
          IF (PRESENT(ierr)) ierr = 3438
          RETURN
       ENDIF
       ! get present position
       CALL manager_state(manager,old_step, ierror)
       IF (ierror /= 0) THEN
          IF (PRESENT(ierr)) ierr=ierror
          RETURN
       ENDIF

       ! get new position
       CALL manager_step_at_date(manager, new_date, new_step, lfit, ierr)
       IF (ierror /= 0) THEN
          IF (PRESENT(ierr)) ierr=ierror
          RETURN
       ENDIF
       IF (.NOT. lfit) THEN
          IF (PRESENT(ierr)) ierr = 3439
          RETURN
       ENDIF

       offset = new_step - old_step
       ! set to new position
       CALL manager_init(manager,offset, ierror)  
       IF (ierror /= 0) THEN
          IF (PRESENT(ierr)) ierr=ierror
          RETURN
       ENDIF
    ELSE
      ! reset the manager with changing the start date

       ! get current_date of manager
      CALL manager_state (manager, wrk_day, ierror)
       IF (ierror /= 0) THEN
          IF (PRESENT(ierr)) ierr=ierror
          RETURN
       ENDIF
      CALL date_get(wrk_day, idays, isecs, ierror)
       IF (ierror /= 0) THEN
          IF (PRESENT(ierr)) ierr=ierror
          RETURN
       ENDIF
      
      CALL copy_date(new_date, wrk_day,ierror)
       IF (ierror /= 0) THEN
          IF (PRESENT(ierr)) ierr=ierror
          RETURN
       ENDIF

      ! wrk_day = new_date - manager %current_date
      CALL add_date (-idays, -isecs, wrk_day, ierror) 
       IF (ierror /= 0) THEN
          IF (PRESENT(ierr)) ierr=ierror
          RETURN
       ENDIF

      IF (ierr /= 0) RETURN
      CALL date_get (wrk_day, idays, isecs, ierror)
       IF (ierror /= 0) THEN
          IF (PRESENT(ierr)) ierr=ierror
          RETURN
       ENDIF

      ! correct start date of manager
      CALL add_date (idays,isecs,manager%start_date, ierror)
       IF (ierror /= 0) THEN
          IF (PRESENT(ierr)) ierr=ierror
          RETURN
       ENDIF

      ! the time step is not changed

    END IF
    
  END SUBROUTINE timer_rewind_manager

END MODULE messy_main_timer_manager
