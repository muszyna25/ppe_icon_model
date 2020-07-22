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

  USE mo_kind,               ONLY: wp, sp, dp
  USE mo_exception,          ONLY: finish
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var
  USE mo_var_metadata_types, ONLY: t_hor_interp_meta, t_vert_interp_meta, &
    &                              t_post_op_meta, VINTP_TYPE_LIST
  USE mo_util_string,        ONLY: toupper

  IMPLICIT NONE
  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_var_metadata'

  PUBLIC  :: create_hor_interp_metadata
  PUBLIC  :: create_vert_interp_metadata
  PUBLIC  :: post_op
  PUBLIC  :: vintp_types
  PUBLIC  :: vintp_type_id

CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  ! Quasi-constructor for horizontal interpolation meta data
  !
  ! Fills data structure with default values (unless set otherwise).
  FUNCTION create_hor_interp_metadata(hor_intp_type, fallback_type, lonlat_id)    &
    & RESULT(him)
    TYPE(t_hor_interp_meta) :: him
    INTEGER, INTENT(IN), OPTIONAL :: hor_intp_type, fallback_type, lonlat_id

    ! supersede with user definitions
    IF (PRESENT(hor_intp_type)) him%hor_intp_type = hor_intp_type
    IF (PRESENT(fallback_type)) him%fallback_type = fallback_type
    IF (PRESENT(lonlat_id))     him%lonlat_id     = lonlat_id
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
  FUNCTION create_vert_interp_metadata(vert_intp_type, vert_intp_method, l_hires_intp,       &
    & l_restore_fricred, l_loglin, l_extrapol, l_satlimit, l_restore_pbldev, l_pd_limit,     &
    & lower_limit) RESULT(vim)
    TYPE(t_vert_interp_meta) :: vim
    LOGICAL, INTENT(IN), OPTIONAL :: vert_intp_type(SIZE(VINTP_TYPE_LIST))
    INTEGER, INTENT(IN), OPTIONAL :: vert_intp_method 
    LOGICAL, INTENT(IN), OPTIONAL :: l_hires_intp, l_restore_fricred, &
      & l_loglin, l_extrapol, l_satlimit, l_restore_pbldev, l_pd_limit
    REAL(wp), INTENT(IN), OPTIONAL :: lower_limit

    ! supersede with user definitions
    IF (PRESENT(vert_intp_type))    vim%vert_intp_type    = vert_intp_type
    IF (PRESENT(vert_intp_method))  vim%vert_intp_method  = vert_intp_method
    IF (PRESENT(l_hires_intp))      vim%l_hires_intp      = l_hires_intp
    IF (PRESENT(l_restore_fricred)) vim%l_restore_fricred = l_restore_fricred
    IF (PRESENT(l_loglin))          vim%l_loglin          = l_loglin
    IF (PRESENT(l_extrapol))        vim%l_extrapol        = l_extrapol
    IF (PRESENT(l_satlimit))        vim%l_satlimit        = l_satlimit
    IF (PRESENT(l_restore_pbldev))  vim%l_restore_pbldev  = l_restore_pbldev
    IF (PRESENT(l_pd_limit))        vim%l_pd_limit        = l_pd_limit
    IF (PRESENT(lower_limit))       vim%lower_limit       = lower_limit
  END FUNCTION create_vert_interp_metadata

  !----------------------------------------------------------------------------------------
  !
  FUNCTION post_op(ipost_op_type, new_cf, new_grib2, arg1)
    TYPE(t_post_op_meta) :: post_op
    INTEGER,           INTENT(IN), OPTIONAL :: ipost_op_type    !< type of post-processing operation
    TYPE(t_cf_var),    INTENT(IN), OPTIONAL :: new_cf           !< CF information of modified field
    TYPE(t_grib2_var), INTENT(IN), OPTIONAL :: new_grib2        !< GRIB2 information of modified field
    CLASS(*),          INTENT(IN), OPTIONAL :: arg1             !< post-op argument (e.g. scaling factor)

    IF (PRESENT(ipost_op_type)) post_op%ipost_op_type = ipost_op_type
    IF (PRESENT(arg1)) THEN
      SELECT TYPE(arg1)
      TYPE is (INTEGER)
        post_op%arg1%ival = arg1
      TYPE is (REAL(dp))
        post_op%arg1%rval = arg1
      TYPE is (REAL(sp))
        post_op%arg1%sval = arg1
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

END MODULE mo_var_metadata

