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
!! The restart writing code only ever sees the bunch of 2D pointers returned by getLevelPointers().

MODULE mo_restart_var_data

  USE mo_dynamics_config,    ONLY: iequations
  USE mo_exception,          ONLY: finish
  USE mo_fortran_tools,      ONLY: t_ptr_2d, t_ptr_2d_sp, t_ptr_2d_int
  USE mo_grid_config,        ONLY: l_limited_area
  USE mo_impl_constants,     ONLY: IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER, INH_ATMOSPHERE, &
    &                              TLEV_NNOW, TLEV_NNOW_RCF, SUCCESS, LEAPFROG_EXPL, LEAPFROG_SI, &
    &                              REAL_T, SINGLE_T, INT_T
#ifdef DEBUG
  USE mo_io_units,           ONLY: nerr
#endif
  USE mo_kind,               ONLY: dp, sp
  USE mo_linked_list,        ONLY: t_list_element, t_var_list
  USE mo_parallel_config,    ONLY: nproma
  USE mo_util_string,        ONLY: int2string
  USE mo_var_list,           ONLY: nvar_lists, var_lists, get_var_timelevel
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_multifile_restart_util, ONLY: dataPtrs_t
#ifndef __NO_ICON_ATMO__
  USE mo_ha_dyn_config, ONLY: ha_dyn_config
#endif
#ifdef _OPENACC
  USE mo_mpi,                       ONLY: i_am_accel_node
#endif

  IMPLICIT NONE

  PUBLIC :: t_RestartVarData
  PUBLIC :: createRestartVarData
  PUBLIC :: getLevelPointers
  PUBLIC :: has_valid_time_level

  INTERFACE getLevelPointers
    MODULE PROCEDURE getLevelPointers_dp
    MODULE PROCEDURE getLevelPointers_sp
    MODULE PROCEDURE getLevelPointers_int
    MODULE PROCEDURE getLevelPointers_hub
  END INTERFACE

  PRIVATE

  ! All the info that's required about a variable for restart purposes.
  TYPE t_RestartVarData
    REAL(dp), POINTER :: r_ptr(:,:,:,:,:)
    REAL(sp), POINTER :: s_ptr(:,:,:,:,:)
    INTEGER,  POINTER :: i_ptr(:,:,:,:,:)
    TYPE(t_var_metadata) :: info
  CONTAINS
    PROCEDURE :: isDoublePrecision => restartVarData_isDoublePrecision
    PROCEDURE :: getDatatype => restartVarData_getDatatype
  END TYPE t_RestartVarData

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_var_data"

CONTAINS

  LOGICAL FUNCTION restartVarData_isDoublePrecision(me) RESULT(resultVar)
    CLASS(t_RestartVarData), INTENT(IN) :: me
    CHARACTER(*), PARAMETER :: routine = modname//":restartVarData_isDoublePrecision"

    SELECT CASE(me%info%data_type)
    CASE(REAL_T) 
      resultVar = .TRUE.
    CASE(SINGLE_T, INT_T)
      resultVar = .FALSE.
    CASE DEFAULT
      CALL finish(routine, "assertion failed: unexpected type of restart variable '"//TRIM(me%info%NAME)//"'")
    END SELECT
  END FUNCTION restartVarData_isDoublePrecision

  INTEGER FUNCTION restartVarData_getDatatype(me) RESULT(resultVar)
    CLASS(t_RestartVarData), INTENT(IN) :: me
    
    resultVar = me%info%data_type
  END FUNCTION restartVarData_getDatatype

  LOGICAL FUNCTION wantVarlist(varlist, patch_id, modelType) RESULT(resultVar)
    TYPE(t_var_list), INTENT(IN) :: varlist
    INTEGER, VALUE :: patch_id
    CHARACTER(LEN = *), INTENT(IN) :: modelType

    resultVar = .FALSE.
    IF(.NOT. varlist%p%lrestart) RETURN
    IF(varlist%p%patch_id /= patch_id) RETURN
    IF(varlist%p%model_type /= modelType) RETURN
    resultVar = .TRUE.
  END FUNCTION wantVarlist

  ! compute the number of  restart variables for the given logical patch.
  INTEGER FUNCTION countRestartVariables(patch_id, modelType) RESULT(resultVar)
    INTEGER, VALUE :: patch_id
    CHARACTER(LEN = *), INTENT(IN) :: modelType
    INTEGER :: i, fld_cnt
    TYPE(t_list_element), POINTER :: element

    resultVar = 0
    DO i = 1, nvar_lists
      IF(.NOT.wantVarlist(var_lists(i), patch_id, modelType)) CYCLE
      ! check, if the list has valid restart fields
      fld_cnt = 0
      element => var_lists(i)%p%first_list_element
      DO
        IF(.NOT. ASSOCIATED(element)) EXIT
        IF(element%field%info%lrestart) fld_cnt = fld_cnt + 1
        element => element%next_list_element
      ENDDO
      resultVar = resultVar + fld_cnt
    ENDDO
  END FUNCTION countRestartVariables

  FUNCTION createRestartVarData(patch_id, modelType, out_restartType) RESULT(resultVar)
    TYPE(t_RestartVarData), POINTER :: resultVar(:)
    INTEGER, VALUE :: patch_id
    CHARACTER(LEN = *), INTENT(IN) :: modelType
    INTEGER, INTENT(OUT) :: out_restartType
    INTEGER :: varCount, varIndex, error, curList
    TYPE(t_list_element), POINTER :: element
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":createRestartVarData"

    ! init. main variables
    resultVar => NULL()
    out_restartType = -1
    ! counts number of restart variables for this file (logical patch ident)
    varCount = countRestartVariables(patch_id, modelType)
#ifdef DEBUG
    WRITE(nerr, *) routine//' numvars = '//TRIM(int2string(varCount))
#endif
    IF(varCount <= 0) RETURN
    ! allocate the array of restart variables
    ALLOCATE(resultVar(varCount), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    ! fill the array of restart variables
    varIndex = 1
    DO curList = 1, nvar_lists
      IF(.NOT.wantVarlist(var_lists(curList), patch_id, modelType)) CYCLE
      IF(out_restartType == -1) THEN
        out_restartType = var_lists(curList)%p%restart_type
      ELSE IF(out_restartType /= var_lists(curList)%p%restart_type) THEN
        CALL finish(routine, "assertion failed: var_lists contains inconsistent restart_type values")
      END IF
      ! check, if the list has valid restart fields
      element => var_lists(curList)%p%first_list_element
      DO
        IF(.NOT. ASSOCIATED(element)) EXIT
        IF(element%field%info%lrestart) THEN
          resultVar(varIndex)%info = element%field%info
          resultVar(varIndex)%r_ptr => element%field%r_ptr
          resultVar(varIndex)%s_ptr => element%field%s_ptr
          resultVar(varIndex)%i_ptr => element%field%i_ptr
          varIndex = varIndex + 1
        END IF
        element => element%next_list_element
      END DO
    END DO
    IF(varIndex /= varCount + 1) CALL finish(routine, "assertion failed: wrong restart variable count")
  END FUNCTION createRestartVarData

  SUBROUTINE getLevelPointers_helper(vMeta, nindex, nlevs, var_ref_pos)
    TYPE(t_var_metadata), INTENT(IN) :: vMeta
    INTEGER, INTENT(OUT) :: nindex, nlevs, var_ref_pos
    INTEGER :: var_ref_pos_max
    CHARACTER(LEN=*), PARAMETER :: routine = modname//":getLevelPointers_helper"

    IF (vMeta%used_dimensions(1) /= nproma) &
      & CALL finish(routine,'1st dim is not nproma: '//TRIM(vMeta%name))
    nindex = MERGE(vMeta%ncontained, 1, vMeta%lcontained)
    var_ref_pos_max = MERGE(3, 4, vMeta%ndims .EQ. 2)
    nlevs = 1
    IF (vMeta%ndims /= 2) nlevs = vMeta%used_dimensions(2)
    var_ref_pos = MERGE(vMeta%var_ref_pos, vMeta%ndims + 1, vMeta%lcontained)
    IF (.NOT.(vMeta%ndims .EQ. 2 .OR. vMeta%ndims .EQ. 3)) &
      & CALL finish(routine, "'"//TRIM(vMeta%NAME)//"': "//&
        &TRIM(int2string(vMeta%ndims))//"d arrays not handled yet")
    IF (var_ref_pos .GT. var_ref_pos_max .OR. var_ref_pos .LT. 1) &
      CALL finish(routine, "internal error!")
  END SUBROUTINE getLevelPointers_helper

  SUBROUTINE getLevelPointers_hub(vMeta, varData, levelPointers, lCnt)
    TYPE(t_var_metadata), INTENT(IN) :: vMeta
    TYPE(t_RestartVarData) :: varData
    TYPE(dataPtrs_t), INTENT(INOUT) :: levelPointers
    INTEGER, OPTIONAL, INTENT(OUT) :: lCnt
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":getLevelPointers"

    SELECT CASE(vMeta%data_type)
    CASE(REAL_T)
      CALL getLevelPointers_dp(vMeta, varData%r_ptr, levelPointers%d)
      IF (ALLOCATED(levelPointers%s)) DEALLOCATE(levelPointers%s)
      IF (ALLOCATED(levelPointers%i)) DEALLOCATE(levelPointers%i)
    CASE(SINGLE_T)
      CALL getLevelPointers_sp(vMeta, varData%s_ptr, levelPointers%s)
      IF (ALLOCATED(levelPointers%d)) DEALLOCATE(levelPointers%d)
      IF (ALLOCATED(levelPointers%i)) DEALLOCATE(levelPointers%i)
     CASE(INT_T)
      CALL getLevelPointers_int(vMeta, varData%i_ptr, levelPointers%i)
      IF (ALLOCATED(levelPointers%s)) DEALLOCATE(levelPointers%s)
      IF (ALLOCATED(levelPointers%d)) DEALLOCATE(levelPointers%d)
    CASE DEFAULT
      CALL finish(routine, "datatype not recognized ")
    END SELECT
    IF(PRESENT(lCnt)) THEN
      lCnt = 1
      IF (vMeta%ndims /= 2) lCnt = vMeta%used_dimensions(2)
    END IF
  END SUBROUTINE getLevelPointers_hub

  SUBROUTINE getLevelPointers_dp(varMetadata, varData, levelPointers)
    TYPE(t_var_metadata), INTENT(IN) :: varMetadata
    REAL(dp), TARGET :: varData(:,:,:,:,:)
    TYPE(t_ptr_2d), ALLOCATABLE, INTENT(INOUT) :: levelPointers(:)
    INTEGER :: nindex, nlevs, var_ref_pos, error, jk
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":getLevelPointers"

    CALL getLevelPointers_helper(varMetadata, nindex, nlevs, var_ref_pos)
    IF(ALLOCATED(levelPointers)) THEN
      IF(SIZE(levelPointers) /= nlevs) DEALLOCATE(levelPointers)
    END IF
    IF(.NOT.ALLOCATED(levelPointers)) THEN
      ALLOCATE(levelPointers(nlevs), STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    END IF
    IF (varMetadata%ndims .EQ. 2) THEN
      SELECT CASE(var_ref_pos)
      CASE (1)
        levelPointers(1)%p => varData(nindex,:,:,1,1)
      CASE (2)
        levelPointers(1)%p => varData(:,nindex,:,1,1)
      CASE (3)
        levelPointers(1)%p => varData(:,:,nindex,1,1)
      END SELECT
!$ACC UPDATE HOST(levelPointers(1)%p) IF ( i_am_accel_node )
    ELSE
      DO jk = 1, nlevs
        SELECT CASE(var_ref_pos)
        CASE (1)
          levelPointers(jk)%p => varData(nindex,:,jk,:,1)
        CASE (2)
          levelPointers(jk)%p => varData(:,nindex,jk,:,1)
        CASE (3)
          levelPointers(jk)%p => varData(:,jk,nindex,:,1)
        CASE (4)
          levelPointers(jk)%p => varData(:,jk,:,nindex,1)
        END SELECT
!$ACC UPDATE HOST(levelPointers(jk)%p) IF ( i_am_accel_node )
      END DO
    END IF
  END SUBROUTINE getLevelPointers_dp

  SUBROUTINE getLevelPointers_sp(varMetadata, varData, levelPointers)
    TYPE(t_var_metadata), INTENT(IN) :: varMetadata
    REAL(sp), TARGET :: varData(:,:,:,:,:)
    TYPE(t_ptr_2d_sp), ALLOCATABLE, INTENT(INOUT) :: levelPointers(:)
    INTEGER :: nindex, nlevs, var_ref_pos, error, jk
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":getLevelPointers"

    CALL getLevelPointers_helper(varMetadata, nindex, nlevs, var_ref_pos)
    IF(ALLOCATED(levelPointers)) THEN
      IF(SIZE(levelPointers) /= nlevs) DEALLOCATE(levelPointers)
    END IF
    IF(.NOT.ALLOCATED(levelPointers)) THEN
      ALLOCATE(levelPointers(nlevs), STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    END IF
    IF(varMetadata%ndims .EQ. 2) THEN
      SELECT CASE(var_ref_pos)
      CASE (1)
        levelPointers(1)%p => varData(nindex,:,:,1,1)
      CASE (2)
        levelPointers(1)%p => varData(:,nindex,:,1,1)
      CASE (3)
        levelPointers(1)%p => varData(:,:,nindex,1,1)
      END SELECT
!$ACC UPDATE HOST(levelPointers(1)%p) IF ( i_am_accel_node )
     ELSE
       DO jk = 1, nlevs
         SELECT CASE(var_ref_pos)
         CASE (1)
           levelPointers(jk)%p => varData(nindex,:,jk,:,1)
         CASE (2)
           levelPointers(jk)%p => varData(:,nindex,jk,:,1)
         CASE (3)
           levelPointers(jk)%p => varData(:,jk,nindex,:,1)
         CASE (4)
           levelPointers(jk)%p => varData(:,jk,:,nindex,1)
         END SELECT
!$ACC UPDATE HOST(levelPointers(jk)%p) IF ( i_am_accel_node )
       END DO
    END IF
  END SUBROUTINE getLevelPointers_sp

  SUBROUTINE getLevelPointers_int(varMetadata, varData, levelPointers)
    TYPE(t_var_metadata), INTENT(IN) :: varMetadata
    INTEGER, TARGET :: varData(:,:,:,:,:)
    TYPE(t_ptr_2d_int), ALLOCATABLE, INTENT(INOUT) :: levelPointers(:)
    INTEGER :: nindex, nlevs, var_ref_pos, error, jk
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":getLevelPointers"

    CALL getLevelPointers_helper(varMetadata, nindex, nlevs, var_ref_pos)
    IF(ALLOCATED(levelPointers)) THEN
      IF(SIZE(levelPointers) /= nlevs) DEALLOCATE(levelPointers)
    END IF
    IF(.NOT.ALLOCATED(levelPointers)) THEN
      ALLOCATE(levelPointers(nlevs), STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    END IF
    IF(varMetadata%ndims .EQ. 2) THEN
      SELECT CASE(var_ref_pos)
      CASE (1)
        levelPointers(1)%p => varData(nindex,:,:,1,1)
      CASE (2)
        levelPointers(1)%p => varData(:,nindex,:,1,1)
      CASE (3)
        levelPointers(1)%p => varData(:,:,nindex,1,1)
      END SELECT
!$ACC UPDATE HOST(levelPointers(1)%p) IF ( i_am_accel_node )
    ELSE
      DO jk = 1, nlevs
        SELECT CASE(var_ref_pos)
        CASE (1)
          levelPointers(jk)%p => varData(nindex,:,jk,:,1)
        CASE (2)
          levelPointers(jk)%p => varData(:,nindex,jk,:,1)
        CASE (3)
          levelPointers(jk)%p => varData(:,jk,nindex,:,1)
        CASE (4)
          levelPointers(jk)%p => varData(:,jk,:,nindex,1)
        END SELECT
!$ACC UPDATE HOST(levelPointers(jk)%p) IF ( i_am_accel_node )
      END DO
    END IF
  END SUBROUTINE getLevelPointers_int

  ! Returns true, if the time level of the given field is valid, else false.
  LOGICAL FUNCTION has_valid_time_level(p_info, domain, nnew, nnew_rcf) RESULT(resultVar)
    TYPE(t_var_metadata), INTENT(IN) :: p_info
    INTEGER, VALUE :: domain, nnew, nnew_rcf
    INTEGER :: time_level
    LOGICAL :: lskip_timelev, lskip_extra_timelevs
    CHARACTER(LEN = *), PARAMETER :: routine = modname//':has_valid_time_level'

    resultVar = .FALSE.
    IF (.NOT. p_info%lrestart) RETURN
#ifndef __NO_ICON_ATMO__
    lskip_timelev = .FALSE.
    lskip_extra_timelevs = iequations == INH_ATMOSPHERE .AND. &
      &                    .NOT. (l_limited_area .AND. domain == 1)
    ! get time index of the given field
    time_level = get_var_timelevel(p_info)
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
    SELECT CASE (iequations)
    CASE(IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER)
      IF ( lskip_timelev                        &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_EXPL &
          & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_SI) &
          & RETURN
    CASE default
      IF ( lskip_timelev ) RETURN
    END SELECT
#endif
    resultVar = .TRUE.
  END FUNCTION has_valid_time_level

END MODULE mo_restart_var_data
