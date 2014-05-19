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
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_var_metadata

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_impl_constants,     ONLY: VINTP_METHOD_LIN, HINTP_TYPE_LONLAT_RBF
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var
  USE mo_var_metadata_types, ONLY: t_hor_interp_meta, t_vert_interp_meta, &
    &                              t_tracer_meta, t_var_metadata,         &
    &                              t_post_op_meta, VAR_GROUPS,            &
    &                              VINTP_TYPE_LIST, VARNAME_LEN, POST_OP_NONE
  USE mo_action_types,       ONLY: t_var_action_element, t_var_action
  USE mo_util_string,        ONLY: toupper
  USE mo_fortran_tools,      ONLY: assign_if_present

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: create_hor_interp_metadata
  PUBLIC  :: create_vert_interp_metadata
  PUBLIC  :: create_tracer_metadata
  PUBLIC  :: groups
  PUBLIC  :: group_id
  PUBLIC  :: post_op
  PUBLIC  :: vintp_types
  PUBLIC  :: vintp_type_id
  PUBLIC  :: new_action
  PUBLIC  :: actions

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS


  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for horizontal interpolation meta data
  ! 
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION create_hor_interp_metadata(hor_intp_type, lonlat_id)    &
    RESULT(hor_interp_meta)

    TYPE(t_hor_interp_meta) :: hor_interp_meta    
    INTEGER, INTENT(IN), OPTIONAL      :: &
      &  hor_intp_type, lonlat_id

    ! set default values
    hor_interp_meta%hor_intp_type    = HINTP_TYPE_LONLAT_RBF
    hor_interp_meta%lonlat_id        = 0 ! invalid ID

    ! supersede with user definitions
    CALL assign_if_present(hor_interp_meta%hor_intp_type, hor_intp_type)
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
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_var_list:vintp_type_id")

    vintp_type_id = 0
    LOOP_VINTP_TYPES : DO ivintp_type=1,SIZE(VINTP_TYPE_LIST)
      IF (toupper(TRIM(in_str)) == toupper(TRIM(VINTP_TYPE_LIST(ivintp_type)))) THEN
        vintp_type_id = ivintp_type
        EXIT LOOP_VINTP_TYPES
      END IF
    END DO LOOP_VINTP_TYPES
    ! paranoia:
    IF ((vintp_type_id < 1) .OR. (vintp_type_id > SIZE(VINTP_TYPE_LIST))) &
      &  CALL finish(routine, "Invalid vertical interpolation type!")
  END FUNCTION vintp_type_id


  !> Utility function with *a lot* of optional string parameters v1,
  !  v2, v3, v4, ...; mapping those onto a
  !  LOGICAL(DIMENSION=MAX_VAR_GROUPS) according to the "group_id"
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
    &  l_pd_limit, l_restore_sfcinv, l_hires_corr, lower_limit, extrapol_dist)               &
    RESULT(vert_interp_meta)

    TYPE(t_vert_interp_meta) :: vert_interp_meta    
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  vert_intp_type(SIZE(VINTP_TYPE_LIST))
    INTEGER, INTENT(IN), OPTIONAL      :: &
      &  vert_intp_method
    LOGICAL, INTENT(IN), OPTIONAL      :: &
      &  l_hires_intp, l_restore_fricred, l_loglin, &
      &  l_extrapol, l_satlimit, l_restore_pbldev,  &
      &  l_pd_limit, l_restore_sfcinv, l_hires_corr
    REAL(wp), INTENT(IN), OPTIONAL     :: &
      &  lower_limit, extrapol_dist

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
    vert_interp_meta%l_restore_sfcinv  = .FALSE.
    vert_interp_meta%l_hires_corr      = .FALSE.
    vert_interp_meta%lower_limit       = 0._wp
    vert_interp_meta%extrapol_dist     = 500._wp
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
    CALL assign_if_present(vert_interp_meta%l_restore_sfcinv   , l_restore_sfcinv  )
    CALL assign_if_present(vert_interp_meta%l_hires_corr       , l_hires_corr      )
    CALL assign_if_present(vert_interp_meta%lower_limit        , lower_limit       )
    CALL assign_if_present(vert_interp_meta%extrapol_dist      , extrapol_dist     )

  END FUNCTION create_vert_interp_metadata


  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for tracer metadata
  ! (public routine. Can be used for two things:
  ! 1.) default settings: If used without any argument, it gives back a variable
  !     of type t_tracer_meta, containing the default settings.
  ! 2.) Setting of metadata: If used with arguments, it gives back a variable 
  !     of type t_tracer_meta, containing the default settings except for those components 
  !     which are given in the argument list.
  !
  ! Comment by DR: Maybe for the future one could define different sets of default values
  ! for different groups of ART species.
  ! 
  FUNCTION create_tracer_metadata(lis_tracer,tracer_class,                              &
    &                            ihadv_tracer, ivadv_tracer, lturb_tracer,              &
    &                            lsed_tracer, ldep_tracer, lconv_tracer,                &
    &                            lwash_tracer,rdiameter_tracer, rrho_tracer,            &
    &                            halflife_tracer,imis_tracer,lifetime_tracer) RESULT(tracer_meta) 

    LOGICAL, INTENT(IN), OPTIONAL :: lis_tracer      ! this is a tracer field (TRUE/FALSE)
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: tracer_class  ! type of tracer (cloud, volcash, radioact,...) 
    INTEGER, INTENT(IN), OPTIONAL :: ihadv_tracer    ! method for horizontal transport
    INTEGER, INTENT(IN), OPTIONAL :: ivadv_tracer    ! method for vertical transport
    LOGICAL, INTENT(IN), OPTIONAL :: lturb_tracer    ! turbulent transport (TRUE/FALSE)
    LOGICAL, INTENT(IN), OPTIONAL :: lsed_tracer     ! sedimentation (TRUE/FALSE)
    LOGICAL, INTENT(IN), OPTIONAL :: ldep_tracer     ! dry deposition (TRUE/FALSE)  
    LOGICAL, INTENT(IN), OPTIONAL :: lconv_tracer    ! convection  (TRUE/FALSE)
    LOGICAL, INTENT(IN), OPTIONAL :: lwash_tracer    ! washout (TRUE/FALSE)
    REAL(wp),INTENT(IN), OPTIONAL :: rdiameter_tracer! particle diameter in m
    REAL(wp),INTENT(IN), OPTIONAL :: rrho_tracer     ! particle density in kg m^-3
    REAL(wp),INTENT(IN), OPTIONAL :: halflife_tracer ! radioactive half-life in s^-1
    INTEGER, INTENT(IN), OPTIONAL :: imis_tracer     ! IMIS number
    REAL(wp),INTENT(IN), OPTIONAL :: lifetime_tracer ! lifetime of a chemical tracer

    TYPE(t_tracer_meta) :: tracer_meta               ! tracer metadata

    ! lis_tracer
    IF ( PRESENT(lis_tracer) ) THEN
      tracer_meta%lis_tracer = lis_tracer
    ELSE
      tracer_meta%lis_tracer = .FALSE.
    ENDIF

    ! tracer_class
    IF ( PRESENT(tracer_class) ) THEN
      tracer_meta%tracer_class = tracer_class
    ELSE
      tracer_meta%tracer_class = "cloud"  ! USE cloud as standard
    ENDIF

    ! ihadv_tracer
    IF ( PRESENT(ihadv_tracer) ) THEN
      tracer_meta%ihadv_tracer = ihadv_tracer
    ELSE
      tracer_meta%ihadv_tracer = 2
    ENDIF

    ! ivadv_tracer
    IF ( PRESENT(ivadv_tracer) ) THEN
      tracer_meta%ivadv_tracer = ivadv_tracer
    ELSE
      tracer_meta%ivadv_tracer = 3
    ENDIF

    ! lturb_tracer  
    IF ( PRESENT(lturb_tracer) ) THEN
      tracer_meta%lturb_tracer = lturb_tracer
    ELSE
      tracer_meta%lturb_tracer = .FALSE.
    ENDIF

    ! lsed_tracer
    IF ( PRESENT(lsed_tracer) ) THEN
      tracer_meta%lsed_tracer = lsed_tracer
    ELSE
      tracer_meta%lsed_tracer = .FALSE.
    ENDIF

    ! ldep_tracer
    IF ( PRESENT(ldep_tracer) ) THEN
      tracer_meta%ldep_tracer = ldep_tracer
    ELSE
      tracer_meta%ldep_tracer = .FALSE.
    ENDIF

    ! lconv_tracer
    IF ( PRESENT(lconv_tracer) ) THEN
      tracer_meta%lconv_tracer = lconv_tracer
    ELSE
      tracer_meta%lconv_tracer = .FALSE.
    ENDIF

    ! lwash_tracer
    IF ( PRESENT(lwash_tracer) ) THEN
      tracer_meta%lwash_tracer = lwash_tracer
    ELSE
      tracer_meta%lwash_tracer = .FALSE.
    ENDIF

    ! rdiameter_tracer
    IF ( PRESENT(rdiameter_tracer) ) THEN
      tracer_meta%rdiameter_tracer = rdiameter_tracer
    ELSE
      tracer_meta%rdiameter_tracer = 1.0e-6_wp ! particle diameter in m
    ENDIF

    ! rrho_tracer
    IF ( PRESENT(rrho_tracer) ) THEN
      tracer_meta%rrho_tracer = rrho_tracer
    ELSE
      tracer_meta%rrho_tracer = 1000.0_wp      ! particle density in kg m^-3
    ENDIF

    ! halflife_tracer
    IF ( PRESENT(halflife_tracer) ) THEN
      tracer_meta%halflife_tracer = halflife_tracer
    ELSE
      tracer_meta%halflife_tracer = 0.0_wp     ! half-life in s^-1
    ENDIF

    ! imis_tracer
    IF ( PRESENT(imis_tracer) ) THEN
      tracer_meta%imis_tracer = imis_tracer
    ELSE
      tracer_meta%imis_tracer = -999           ! IMIS number
    ENDIF
    
    ! lifetime_tracer
    IF ( PRESENT(lifetime_tracer) ) THEN
      tracer_meta%lifetime_tracer = lifetime_tracer
    ELSE
      tracer_meta%lifetime_tracer = -999       ! lifetime of a chemical tracer
    ENDIF

  END FUNCTION create_tracer_metadata


  !------------------------------------------------------------------------------------------------
  ! HANDLING OF VARIABLE GROUPS
  !------------------------------------------------------------------------------------------------
  
  !> Implements a (somewhat randomly chosen) one-to-one mapping
  !  between a string and an integer ID number between 1 and
  !  MAX_VAR_GROUPS.
  !
  FUNCTION group_id(in_str)
    INTEGER                      :: group_id, igrp
    CHARACTER(LEN=*), INTENT(IN) :: in_str
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_var_list:group_id")

    group_id = 0
    LOOP_GROUPS : DO igrp=1,SIZE(VAR_GROUPS)
      IF (toupper(TRIM(in_str)) == toupper(TRIM(VAR_GROUPS(igrp)))) THEN
        group_id = igrp
        EXIT LOOP_GROUPS
      END IF
    END DO LOOP_GROUPS
    ! paranoia:
    IF ((group_id < 1) .OR. (group_id > SIZE(VAR_GROUPS))) &
      &  CALL finish(routine, "Invalid group ID!")

  END FUNCTION group_id




  !----------------------------------------------------------------------------------------
  !
  !> Utility function with *a lot* of optional string parameters g1,
  !  g2, g3, g4, ...; mapping those onto a
  !  LOGICAL(DIMENSION=MAX_VAR_GROUPS) according to the "group_id"
  !  function.
  !
  FUNCTION groups(g01, g02, g03, g04, g05, g06, g07, g08, g09, g10)
    LOGICAL :: groups(SIZE(VAR_GROUPS))
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: &
      &   g01, g02, g03, g04, g05, g06, g07, g08, g09, g10
    
    groups(:) = .FALSE.
    groups(group_id("ALL")) = .TRUE.
    IF (PRESENT(g01)) groups(group_id(g01)) = .TRUE.
    IF (PRESENT(g02)) groups(group_id(g02)) = .TRUE.
    IF (PRESENT(g03)) groups(group_id(g03)) = .TRUE.
    IF (PRESENT(g04)) groups(group_id(g04)) = .TRUE.
    IF (PRESENT(g05)) groups(group_id(g05)) = .TRUE.
    IF (PRESENT(g06)) groups(group_id(g06)) = .TRUE.
    IF (PRESENT(g07)) groups(group_id(g07)) = .TRUE.
    IF (PRESENT(g08)) groups(group_id(g08)) = .TRUE.
    IF (PRESENT(g09)) groups(group_id(g09)) = .TRUE.
    IF (PRESENT(g10)) groups(group_id(g10)) = .TRUE.
  END FUNCTION groups


  !----------------------------------------------------------------------------------------
  !
  !> Utility function with *a lot* of optional string parameters g1,
  !  g2, g3, g4, ...; mapping those onto a
  !  LOGICAL(DIMENSION=MAX_VAR_GROUPS) according to the "group_id"
  !  function.
  !
  FUNCTION post_op(ipost_op_type, new_cf, new_grib2, arg1)
    TYPE(t_post_op_meta) :: post_op
    INTEGER,           INTENT(IN), OPTIONAL :: ipost_op_type    !< type of post-processing operation
    TYPE(t_cf_var),    INTENT(IN), OPTIONAL :: new_cf           !< CF information of modified field
    TYPE(t_grib2_var), INTENT(IN), OPTIONAL :: new_grib2        !< GRIB2 information of modified field
    REAL(wp),          INTENT(IN), OPTIONAL :: arg1             !< post-op argument (e.g. scaling factor)

    post_op%ipost_op_type = POST_OP_NONE
    post_op%lnew_cf       = .FALSE.
    post_op%lnew_grib2    = .FALSE.
    post_op%arg1          = 0._wp

    IF (PRESENT(ipost_op_type)) post_op%ipost_op_type = ipost_op_type
    IF (PRESENT(arg1))          post_op%arg1          = arg1
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
  !!
  FUNCTION new_action(actionID, intvl) RESULT(var_action)

    INTEGER           , INTENT(IN) :: actionID  ! action ID  
    CHARACTER(LEN=*)  , INTENT(IN) :: intvl     ! action interval [h]
    TYPE(t_var_action_element)     :: var_action

    ! define var_action
    var_action%actionID   = actionID
    var_action%intvl      = TRIM(intvl)
    var_action%lastActive ='0000-00-00T00:00:00.000Z'  ! init

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

