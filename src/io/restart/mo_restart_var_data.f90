!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The glue layer between the restart writing code and the variables:
!! The restart writing code only ever sees 3D pointers returned by get_var_3d_ptr.

MODULE mo_restart_var_data

  USE mo_dynamics_config,    ONLY: iequations
  USE mo_exception,          ONLY: finish
  USE mo_fortran_tools,      ONLY: insert_dimension
  USE mo_grid_config,        ONLY: l_limited_area
  USE mo_impl_constants,     ONLY: INH_ATMOSPHERE, TLEV_NNOW, TLEV_NNOW_RCF, SUCCESS
#ifdef DEBUG
  USE mo_io_units,           ONLY: nerr
#endif
  USE mo_kind,               ONLY: dp, sp
  USE mo_util_string,        ONLY: int2string
  USE mo_var,                ONLY: t_var
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_var_metadata,       ONLY: get_var_timelevel
#ifdef _OPENACC
  USE mo_mpi,                       ONLY: i_am_accel_node
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_var_3d_ptr, has_valid_time_level

  INTERFACE get_var_3d_ptr
    MODULE PROCEDURE get_var_3d_ptr_dp
    MODULE PROCEDURE get_var_3d_ptr_sp
    MODULE PROCEDURE get_var_3d_ptr_int
  END INTERFACE get_var_3d_ptr

  CHARACTER(*), PARAMETER :: modname = "mo_restart_var_data"

CONTAINS

  SUBROUTINE get_var_3d_ptr_dp(vd, r_ptr_3d)
    TYPE(t_var), POINTER, INTENT(IN) :: vd
    REAL(dp), POINTER, INTENT(OUT) :: r_ptr_3d(:,:,:)
    REAL(dp), POINTER :: r_ptr_2d(:,:)
    INTEGER :: nindex, nlevs, var_ref_pos
    CHARACTER(*), PARAMETER :: routine = modname//":get_var_3d_ptr_dp"

    nindex = MERGE(vd%info%ncontained, 1, vd%info%lcontained)
    nlevs = MERGE(vd%info%used_dimensions(2), 1, vd%info%ndims /= 2)
    var_ref_pos = MERGE(vd%info%var_ref_pos, vd%info%ndims + 1, vd%info%lcontained)
    SELECT CASE (vd%info%ndims)
    CASE (2)
      SELECT CASE(var_ref_pos)
      CASE (1)
        r_ptr_2d => vd%r_ptr(nindex,:,:,1,1)
        CALL insert_dimension(r_ptr_3d, r_ptr_2d, 2)
      CASE (2)
        r_ptr_3d => vd%r_ptr(:,nindex:nindex,:,1,1)
      CASE (3)
        r_ptr_2d => vd%r_ptr(:,:,nindex,1,1)
        CALL insert_dimension(r_ptr_3d, r_ptr_2d, 2)
      CASE default
        CALL finish(routine, "internal error!")
      END SELECT
    CASE (3)
      SELECT CASE(var_ref_pos)
      CASE (1)
        r_ptr_3d => vd%r_ptr(nindex,:,1:nlevs,:,1)
      CASE (2)
        r_ptr_3d => vd%r_ptr(:,nindex,1:nlevs,:,1)
      CASE (3)
        r_ptr_3d => vd%r_ptr(:,1:nlevs,nindex,:,1)
      CASE (4)
        r_ptr_3d => vd%r_ptr(:,1:nlevs,:,nindex,1)
      CASE default
        CALL finish(routine, "internal error!")
      END SELECT
    CASE DEFAULT
      CALL finish(routine, "'"//TRIM(vd%info%NAME)//"': "//&
           & TRIM(int2string(vd%info%ndims))//"d arrays not handled yet")
    END SELECT
!$ACC UPDATE HOST(r_ptr_3d) IF ( i_am_accel_node )
  END SUBROUTINE get_var_3d_ptr_dp

  SUBROUTINE get_var_3d_ptr_sp(vd, s_ptr_3d)
    TYPE(t_var), POINTER, INTENT(IN) :: vd
    REAL(sp), POINTER, INTENT(OUT) :: s_ptr_3d(:,:,:)
    REAL(sp), POINTER :: s_ptr_2d(:,:)
    INTEGER :: nindex, nlevs, var_ref_pos
    CHARACTER(*), PARAMETER :: routine = modname//":get_var_3d_ptr_sp"

    nindex = MERGE(vd%info%ncontained, 1, vd%info%lcontained)
    nlevs = MERGE(vd%info%used_dimensions(2), 1, vd%info%ndims /= 2)
    var_ref_pos = MERGE(vd%info%var_ref_pos, vd%info%ndims + 1, vd%info%lcontained)
    SELECT CASE (vd%info%ndims)
    CASE (2)
      SELECT CASE(var_ref_pos)
      CASE (1)
        s_ptr_2d => vd%s_ptr(nindex,:,:,1,1)
        CALL insert_dimension(s_ptr_3d, s_ptr_2d, 2)
      CASE (2)
        s_ptr_3d => vd%s_ptr(:,nindex:nindex,:,1,1)
      CASE (3)
        s_ptr_2d => vd%s_ptr(:,:,nindex,1,1)
        CALL insert_dimension(s_ptr_3d, s_ptr_2d, 2)
      CASE default
        CALL finish(routine, "internal error!")
      END SELECT
    CASE (3)
      SELECT CASE(var_ref_pos)
      CASE (1)
        s_ptr_3d => vd%s_ptr(nindex,:,1:nlevs,:,1)
      CASE (2)
        s_ptr_3d => vd%s_ptr(:,nindex,1:nlevs,:,1)
      CASE (3)
        s_ptr_3d => vd%s_ptr(:,1:nlevs,nindex,:,1)
      CASE (4)
        s_ptr_3d => vd%s_ptr(:,1:nlevs,:,nindex,1)
      CASE default
        CALL finish(routine, "internal error!")
      END SELECT
    CASE DEFAULT
      CALL finish(routine, "'"//TRIM(vd%info%NAME)//"': "//&
           & TRIM(int2string(vd%info%ndims))//"d arrays not handled yet")
    END SELECT
!$ACC UPDATE HOST(s_ptr_3d) IF ( i_am_accel_node )
  END SUBROUTINE get_var_3d_ptr_sp

  SUBROUTINE get_var_3d_ptr_int(vd, i_ptr_3d)
    TYPE(t_var), POINTER, INTENT(IN) :: vd
    INTEGER, POINTER, INTENT(OUT) :: i_ptr_3d(:,:,:)
    INTEGER, POINTER :: i_ptr_2d(:,:)
    INTEGER :: nindex, nlevs, var_ref_pos
    CHARACTER(*), PARAMETER :: routine = modname//":get_var_3d_ptr_i"

    nindex = MERGE(vd%info%ncontained, 1, vd%info%lcontained)
    nlevs = MERGE(vd%info%used_dimensions(2), 1, vd%info%ndims /= 2)
    var_ref_pos = MERGE(vd%info%var_ref_pos, vd%info%ndims + 1, vd%info%lcontained)
    SELECT CASE (vd%info%ndims)
    CASE (2)
      SELECT CASE(var_ref_pos)
      CASE (1)
        i_ptr_2d => vd%i_ptr(nindex,:,:,1,1)
        CALL insert_dimension(i_ptr_3d, i_ptr_2d, 2)
      CASE (2)
        i_ptr_3d => vd%i_ptr(:,nindex:nindex,:,1,1)
      CASE (3)
        i_ptr_2d => vd%i_ptr(:,:,nindex,1,1)
        CALL insert_dimension(i_ptr_3d, i_ptr_2d, 2)
      CASE default
        CALL finish(routine, "internal error!")
      END SELECT
    CASE (3)
      SELECT CASE(var_ref_pos)
      CASE (1)
        i_ptr_3d => vd%i_ptr(nindex,:,1:nlevs,:,1)
      CASE (2)
        i_ptr_3d => vd%i_ptr(:,nindex,1:nlevs,:,1)
      CASE (3)
        i_ptr_3d => vd%i_ptr(:,1:nlevs,nindex,:,1)
      CASE (4)
        i_ptr_3d => vd%i_ptr(:,1:nlevs,:,nindex,1)
      CASE default
        CALL finish(routine, "internal error!")
      END SELECT
    CASE DEFAULT
      CALL finish(routine, "'"//TRIM(vd%info%NAME)//"': "//&
           & TRIM(int2string(vd%info%ndims))//"d arrays not handled yet")
    END SELECT
!$ACC UPDATE HOST(i_ptr_3d) IF ( i_am_accel_node )
  END SUBROUTINE get_var_3d_ptr_int

  ! Returns true, if the time level of the given field is valid, else false.
  LOGICAL FUNCTION has_valid_time_level(p_info, patch_id, nnew, nnew_rcf) RESULT(has_vtl)
    TYPE(t_var_metadata), INTENT(IN) :: p_info
    INTEGER, INTENT(in) :: patch_id, nnew, nnew_rcf
    INTEGER :: time_level
    LOGICAL :: lskip_timelev, lskip_extra_timelevs
    CHARACTER(*), PARAMETER :: routine = modname//':has_valid_time_level'

    has_vtl = .FALSE.
    IF (.NOT. p_info%lrestart) RETURN
#ifndef __NO_ICON_ATMO__
    lskip_timelev = .FALSE.
    lskip_extra_timelevs = iequations == INH_ATMOSPHERE .AND. &
      &                    .NOT. (l_limited_area .AND. patch_id == 1)
    ! get time index of the given field
    time_level = get_var_timelevel(p_info%name)
    !TODO: I found the `time_level >= 0` condition IN the async restart code ONLY. Check whether it should be removed OR NOT.
    IF(time_level >= 0) THEN
      ! get information about time level to be skipped for current field
      IF (p_info%tlev_source == TLEV_NNOW) THEN
        IF (time_level == nnew) lskip_timelev = .TRUE.
        ! this is needed to skip the extra time levels allocated for nesting
        IF (lskip_extra_timelevs .AND. time_level > 2) lskip_timelev = .TRUE.
      ELSE IF (p_info%tlev_source == TLEV_NNOW_RCF) THEN
        IF (time_level == nnew_rcf) lskip_timelev = .TRUE.
      ENDIF
    ENDIF
    !
    IF ( lskip_timelev ) RETURN
#endif
    has_vtl = .TRUE.
  END FUNCTION has_valid_time_level

END MODULE mo_restart_var_data
