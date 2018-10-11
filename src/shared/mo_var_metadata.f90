!>
!! Utility funtions for handling field specific meta information
!!
!! Contains utility funtions which are used for defining variable specific
!! meta information. These have nothing to do with var lists itself. That's
!! why they have been moved here from mo_var_list.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2014-01-22)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_var_metadata

  USE mo_kind,               ONLY: wp, sp
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: VINTP_METHOD_LIN, HINTP_TYPE_LONLAT_RBF, &
    &                              MAX_CHAR_LENGTH
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var
  USE mo_var_groups,         ONLY: var_groups_dyn
  USE mo_var_metadata_types, ONLY: t_hor_interp_meta, t_vert_interp_meta, &
    &                              t_union_vals,                          &
    &                              t_post_op_meta,                        &
    &                              VINTP_TYPE_LIST, POST_OP_NONE
  USE mo_action_types,       ONLY: t_var_action_element, t_var_action
  USE mo_util_string,        ONLY: toupper
  USE mo_fortran_tools,      ONLY: assign_if_present
  USE mo_time_config,        ONLY: time_config
  USE mtime,                 ONLY: datetime, newDatetime, deallocateDatetime,    &
    &                              timedelta, newTimedelta, deallocateTimedelta, &
    &                              OPERATOR(+), dateTimeToString, MAX_DATETIME_STR_LEN

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_var_metadata'


  PUBLIC  :: create_hor_interp_metadata
  PUBLIC  :: create_vert_interp_metadata
  PUBLIC  :: post_op
  PUBLIC  :: vintp_types
  PUBLIC  :: vintp_type_id
  PUBLIC  :: new_action
  PUBLIC  :: actions


CONTAINS


  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for horizontal interpolation meta data
  !
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION create_hor_interp_metadata(hor_intp_type, fallback_type, lonlat_id)    &
    RESULT(hor_interp_meta)

    TYPE(t_hor_interp_meta) :: hor_interp_meta
    INTEGER, INTENT(IN), OPTIONAL      :: &
      &  hor_intp_type, fallback_type, lonlat_id

    ! set default values
    hor_interp_meta%hor_intp_type    = HINTP_TYPE_LONLAT_RBF
    hor_interp_meta%fallback_type    = HINTP_TYPE_LONLAT_RBF
    hor_interp_meta%lonlat_id        = 0 ! invalid ID

    ! supersede with user definitions
    CALL assign_if_present(hor_interp_meta%hor_intp_type, hor_intp_type)
    CALL assign_if_present(hor_interp_meta%fallback_type, fallback_type)
    CALL assign_if_present(hor_interp_meta%lonlat_id,     lonlat_id)

  END FUNCTION create_hor_interp_metadata


  !------------------------------------------------------------------------------------------------
  ! HANDLING OF VERTICAL INTERPOLATION MODES
  !------------------------------------------------------------------------------------------------

  !> Implements a (somewhat randomly chosen) one-to-one mapping
  !  between a string and an integer ID number between 1 and
  !  MAX_VINTP_TYPES.
  !
  FUNCTION vintp_type_id(in_str)
    INTEGER                      :: vintp_type_id, ivintp_type
    CHARACTER(LEN=*), INTENT(IN) :: in_str
    CHARACTER(len=*), PARAMETER :: routine = modname//"::vintp_type_id"
    CHARACTER(len=len_trim(in_str)) :: in_str_upper
    INTEGER :: n

    vintp_type_id = 0
    in_str_upper = toupper(in_str)
    n = SIZE(VINTP_TYPE_LIST)
    LOOP_VINTP_TYPES : DO ivintp_type=1,n
      IF (in_str_upper == toupper(VINTP_TYPE_LIST(ivintp_type))) THEN
        vintp_type_id = ivintp_type
        EXIT LOOP_VINTP_TYPES
      END IF
    END DO LOOP_VINTP_TYPES
    ! paranoia:
    IF ((vintp_type_id < 1) .OR. (vintp_type_id > n)) &
      &  CALL finish(routine, "Invalid vertical interpolation type!")
  END FUNCTION vintp_type_id


  !> Utility function with *a lot* of optional string parameters v1,
  !  v2, v3, v4, ...; mapping those onto a
  !  LOGICAL(:) according to the "vintp_type_id"
  !  function.
  !
  FUNCTION vintp_types(v01, v02, v03, v04, v05, v06, v07, v08, v09, v10)
    LOGICAL :: vintp_types(SIZE(VINTP_TYPE_LIST))
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: &
      &   v01, v02, v03, v04, v05, v06, v07, v08, v09, v10

    vintp_types(:) = .FALSE.
    IF (PRESENT(v01)) vintp_types(vintp_type_id(v01)) = .TRUE.
    IF (PRESENT(v02)) vintp_types(vintp_type_id(v02)) = .TRUE.
    IF (PRESENT(v03)) vintp_types(vintp_type_id(v03)) = .TRUE.
    IF (PRESENT(v04)) vintp_types(vintp_type_id(v04)) = .TRUE.
    IF (PRESENT(v05)) vintp_types(vintp_type_id(v05)) = .TRUE.
    IF (PRESENT(v06)) vintp_types(vintp_type_id(v06)) = .TRUE.
    IF (PRESENT(v07)) vintp_types(vintp_type_id(v07)) = .TRUE.
    IF (PRESENT(v08)) vintp_types(vintp_type_id(v08)) = .TRUE.
    IF (PRESENT(v09)) vintp_types(vintp_type_id(v09)) = .TRUE.
    IF (PRESENT(v10)) vintp_types(vintp_type_id(v10)) = .TRUE.
  END FUNCTION vintp_types


  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for vertical interpolation meta data
  !
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION create_vert_interp_metadata(vert_intp_type, vert_intp_method,                     &
    &  l_hires_intp, l_restore_fricred, l_loglin, l_extrapol, l_satlimit, l_restore_pbldev,  &
    &  l_pd_limit, lower_limit)               &
    RESULT(vert_interp_meta)

    TYPE(t_vert_interp_meta) :: vert_interp_meta
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  vert_intp_type(SIZE(VINTP_TYPE_LIST))
    INTEGER, INTENT(IN), OPTIONAL      :: &
      &  vert_intp_method
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_hires_intp, l_restore_fricred, l_loglin, &
      &  l_extrapol, l_satlimit, l_restore_pbldev,  &
      &  l_pd_limit
    REAL(wp), INTENT(IN), OPTIONAL     :: &
      &  lower_limit

    ! set default values
    vert_interp_meta%vert_intp_type(:) = .FALSE.
    vert_interp_meta%vert_intp_method  = VINTP_METHOD_LIN
    vert_interp_meta%l_hires_intp      = .FALSE.
    vert_interp_meta%l_restore_fricred = .FALSE.
    vert_interp_meta%l_loglin          = .FALSE.
    vert_interp_meta%l_extrapol        = .TRUE.
    vert_interp_meta%l_satlimit        = .FALSE.
    vert_interp_meta%l_restore_pbldev  = .FALSE.
    vert_interp_meta%l_pd_limit        = .FALSE.
    vert_interp_meta%lower_limit       = 0._wp
    ! supersede with user definitions
    CALL assign_if_present(vert_interp_meta%vert_intp_type     , vert_intp_type    )
    CALL assign_if_present(vert_interp_meta%vert_intp_method   , vert_intp_method  )
    CALL assign_if_present(vert_interp_meta%l_hires_intp       , l_hires_intp      )
    CALL assign_if_present(vert_interp_meta%l_restore_fricred  , l_restore_fricred )
    CALL assign_if_present(vert_interp_meta%l_loglin           , l_loglin          )
    CALL assign_if_present(vert_interp_meta%l_extrapol         , l_extrapol        )
    CALL assign_if_present(vert_interp_meta%l_satlimit         , l_satlimit        )
    CALL assign_if_present(vert_interp_meta%l_restore_pbldev   , l_restore_pbldev  )
    CALL assign_if_present(vert_interp_meta%l_pd_limit         , l_pd_limit        )
    CALL assign_if_present(vert_interp_meta%lower_limit        , lower_limit       )

  END FUNCTION create_vert_interp_metadata


  !----------------------------------------------------------------------------------------
  !
  FUNCTION post_op(ipost_op_type, new_cf, new_grib2, arg1)
    TYPE(t_post_op_meta) :: post_op
    INTEGER,           INTENT(IN), OPTIONAL :: ipost_op_type    !< type of post-processing operation
    TYPE(t_cf_var),    INTENT(IN), OPTIONAL :: new_cf           !< CF information of modified field
    TYPE(t_grib2_var), INTENT(IN), OPTIONAL :: new_grib2        !< GRIB2 information of modified field
    CLASS(*),          INTENT(IN), OPTIONAL :: arg1             !< post-op argument (e.g. scaling factor)


    post_op%ipost_op_type = POST_OP_NONE
    post_op%lnew_cf       = .FALSE.
    post_op%lnew_grib2    = .FALSE.
    post_op%arg1          = t_union_vals( 0._wp, 0._sp, 0, .FALSE.)

    IF (PRESENT(ipost_op_type)) post_op%ipost_op_type = ipost_op_type

    IF (PRESENT(arg1)) THEN
      SELECT TYPE(arg1)
      TYPE is (INTEGER)
        post_op%arg1 = t_union_vals( 0.0_wp, 0.0_sp, arg1, .FALSE.)
      TYPE is (REAL(wp))
        post_op%arg1 = t_union_vals( arg1  , 0.0_sp,    0, .FALSE.)
      TYPE is (REAL(sp))
        post_op%arg1 = t_union_vals( 0.0_wp,   arg1,    0, .FALSE.)
      END SELECT
    ENDIF

    IF (PRESENT(new_cf)) THEN
      post_op%lnew_cf = .TRUE.
      post_op%new_cf  = new_cf
    END IF
    IF (PRESENT(new_grib2)) THEN
      post_op%lnew_grib2 = .TRUE.
      post_op%new_grib2  = new_grib2
    END IF
  END FUNCTION post_op


  !------------------------------------------------------------------------------------------------
  ! HANDLING OF ACTION EVENTS
  !------------------------------------------------------------------------------------------------
  !>
  !! Initialize single variable specific action
  !!
  !! Initialize single variable specific action. A variable named 'var_action'
  !! of type t_var_action_element is initialized.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-13)
  !! Modification by Daniel Reinert, DWD (2014-12-03)
  !! - add optional start and end time arguments
  !!
  FUNCTION new_action(actionTyp, intvl, opt_start, opt_end, opt_ref) RESULT(var_action)

    INTEGER                   , INTENT(IN) :: actionTyp ! type of action
    CHARACTER(LEN=*)          , INTENT(IN) :: intvl     ! action interval [ISO_8601]
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: opt_start ! action start time [ISO_8601]
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: opt_end   ! action end time [ISO_8601]
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: opt_ref   ! action reference time [ISO_8601]

    ! local variables
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = modname//':new_action'
    TYPE(timedelta), POINTER              :: start_offset, end_offset, ref_offset
    TYPE(datetime), TARGET                :: startdatetime, enddatetime, refdatetime
    TYPE(datetime), POINTER               :: dummy_ptr
    TYPE(t_var_action_element)            :: var_action
    TYPE(datetime), POINTER               :: inidatetime
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: iso8601_ini_datetime      ! ISO_8601
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: iso8601_expstart_datetime ! ISO_8601
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: iso8601_end_datetime      ! ISO_8601
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: start, end, ref      ! start, end, and reference time
                                                                  ! in ISO_8601 format
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)   :: start0, end0, ref0   ! start, end, and reference time
                                                                  ! in ISO_8601 format
    !---------------------------------------------------------------------------------

    ! create model ini_datetime in ISO_8601 format
    CALL dateTimeToString(time_config%tc_startdate, iso8601_ini_datetime)
    ! create experiment start in ISO_8601 format
    CALL dateTimeToString(time_config%tc_exp_startdate, iso8601_expstart_datetime)
    ! create model end_datetime in ISO_8601 format
    CALL dateTimeToString(time_config%tc_stopdate, iso8601_end_datetime)

    ! default start time = model initialization time
    start0 = TRIM(iso8601_ini_datetime)
    ! default end time = model end time
    end0 = TRIM(iso8601_end_datetime)
    ! default reference time = experiment start time
    ref0 = TRIM(iso8601_expstart_datetime)


    ! assign modified start time if offset opt_start is present
    IF (PRESENT(opt_start)) THEN
      !
      ! convert model ini datetime from ISO_8601 format to type datetime
      inidatetime  => newDatetime(TRIM(iso8601_ini_datetime))
      !
      ! convert start offset from ISO_8601 to TYPE timedelta
      start_offset => newTimedelta(TRIM(opt_start))
      !
      ! add start offset to model ini date
      startdatetime = inidatetime + start_offset
      ! transform back from TYPE datetime to ISO_8601
      dummy_ptr => startdatetime
      CALL dateTimeToString(dummy_ptr, start0)
      ! cleanup
      CALL deallocateDatetime(inidatetime)
      CALL deallocateTimeDelta(start_offset)
    ENDIF


    ! assign modified end time if offset opt_end is present
    IF (PRESENT(opt_end)) THEN
      !
      ! convert model ini datetime from ISO_8601 format to type datetime
      inidatetime  => newDatetime(TRIM(iso8601_ini_datetime))
      !
      ! convert end offset from ISO_8601 to TYPE timedelta
      end_offset => newTimedelta(TRIM(opt_end))
      !
      ! add end offset to model ini date
      enddatetime = inidatetime + end_offset
      ! transform back from TYPE datetime to ISO_8601
      dummy_ptr => enddatetime
      CALL dateTimeToString(dummy_ptr, end0)
      ! cleanup
      CALL deallocateDatetime(inidatetime)
      CALL deallocateTimeDelta(end_offset)
    ENDIF

    ! assign modified reference time if offset opt_ref is present
    IF (PRESENT(opt_ref)) THEN
      !
      ! convert model ini datetime from ISO_8601 format to type datetime
      inidatetime  => newDatetime(TRIM(iso8601_ini_datetime))
      !
      ! convert ref offset from ISO_8601 to TYPE timedelta
      ref_offset => newTimedelta(TRIM(opt_ref))
      !
      ! add ref offset to model ini date
      refdatetime = inidatetime + ref_offset
      ! transform back from TYPE datetime to ISO_8601
      dummy_ptr => refdatetime
      CALL dateTimeToString(dummy_ptr, ref0)
      ! cleanup
      CALL deallocateDatetime(inidatetime)
      CALL deallocateTimeDelta(ref_offset)
    ENDIF

    ! default start time = model initialization time
    start = TRIM(iso8601_ini_datetime)
    ! default end time = model end time
    end = TRIM(iso8601_end_datetime)
    ! default reference time = experiment start time
    ref = TRIM(iso8601_expstart_datetime)

    !---------------------------------------------------------------------------------

#ifdef _MTIME_DEBUG

    ! CONSISTENCY CHECK: compares the implementation above with an
    ! mtime-based implementation

    ! assign modified start time if offset opt_start is present
    IF (PRESENT(opt_start)) THEN
      !
      ! convert start offset from ISO_8601 to TYPE timedelta
      start_offset => newTimedelta(TRIM(opt_start))
      !
      ! add start offset to model ini date
      startdatetime = time_config%tc_startdate + start_offset
      ! transform back from TYPE datetime to ISO_8601
      dummy_ptr => startdatetime
      CALL dateTimeToString(dummy_ptr, start)
      ! cleanup
      CALL deallocateTimeDelta(start_offset)
    ENDIF


    ! assign modified end time if offset opt_end is present
    IF (PRESENT(opt_end)) THEN
      !
      ! convert end offset from ISO_8601 to TYPE timedelta
      end_offset => newTimedelta(TRIM(opt_end))
      !
      ! add end offset to model ini date
      enddatetime = time_config%tc_startdate + end_offset
      ! transform back from TYPE datetime to ISO_8601
      dummy_ptr => enddatetime
      CALL dateTimeToString(dummy_ptr, end)
      ! cleanup
      CALL deallocateTimeDelta(end_offset)
    ENDIF

    ! assign modified reference time if offset opt_ref is present
    IF (PRESENT(opt_ref)) THEN
      !
      ! convert ref offset from ISO_8601 to TYPE timedelta
      ref_offset => newTimedelta(TRIM(opt_ref))
      !
      ! add ref offset to model ini date
      refdatetime = time_config%tc_startdate + ref_offset
      ! transform back from TYPE datetime to ISO_8601
      dummy_ptr => refdatetime
      CALL dateTimeToString(dummy_ptr, ref)
      ! cleanup
      CALL deallocateTimeDelta(ref_offset)
    ENDIF

    IF ((TRIM(start0) /= TRIM(start)) .OR.   &
      & (TRIM(end0)   /= TRIM(end))   .OR.   &
      & (TRIM(ref0)   /= TRIM(ref))) THEN
      CALL finish(routine, "Error in mtime consistency check!")
    END IF

#endif

    !---------------------------------------------------------------------------------

    ! define var_action
    var_action%actionTyp  = actionTyp
    var_action%intvl      = TRIM(intvl)               ! interval
    var_action%start      = TRIM(start)               ! start
    var_action%end        = TRIM(end)                 ! end
    var_action%ref        = TRIM(ref)                 ! ref date
    var_action%lastActive = TRIM(start)               ! arbitrary init

    !
    ! convert start datetime from ISO_8601 format to type datetime
    dummy_ptr => newDatetime(TRIM(start))
    IF (.NOT. ASSOCIATED(dummy_ptr)) THEN
      CALL finish(routine, "date/time conversion error: "//TRIM(start))
    END IF
    var_Action%EventLastTriggerDate = dummy_ptr    ! arbitrary init

    ! cleanup
    CALL deallocateDatetime(dummy_ptr)

  END FUNCTION new_action


  !>
  !! Generate list (array) of variable specific actions
  !!
  !! Generate list (array) of variable specific actions.
  !! Creates array 'action_list' of type t_var_action
  !
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-13)
  !!
  FUNCTION actions(a01, a02, a03, a04, a05)  RESULT(action_list)

    TYPE(t_var_action_element), INTENT(IN), OPTIONAL :: a01, a02, a03, a04, a05
    TYPE(t_var_action)             :: action_list

    INTEGER :: n_act             ! action counter

    ! create action list
    !
    n_act = 0
    IF (PRESENT(a01))  THEN
      n_act = n_act + 1
      action_list%action(n_act) = a01
    ENDIF

    IF (PRESENT(a02))  THEN
      n_act = n_act + 1
      action_list%action(n_act) = a02
    ENDIF

    IF (PRESENT(a03))  THEN
      n_act = n_act + 1
      action_list%action(n_act) = a03
    ENDIF

    IF (PRESENT(a04))  THEN
      n_act = n_act + 1
      action_list%action(n_act) = a04
    ENDIF

    IF (PRESENT(a05))  THEN
      n_act = n_act + 1
      action_list%action(n_act) = a05
    ENDIF

    action_list%n_actions = n_act

  END FUNCTION actions

END MODULE mo_var_metadata

