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
    USE mo_dynamics_config, ONLY: iequations
    USE mo_exception, ONLY: finish
    USE mo_fortran_tools, ONLY: t_ptr_2d
    USE mo_grid_config, ONLY: l_limited_area
    USE mo_ha_dyn_config, ONLY: ha_dyn_config
    USE mo_impl_constants, ONLY: IHS_ATM_TEMP, IHS_ATM_THETA, ISHALLOW_WATER, INH_ATMOSPHERE, TLEV_NNOW, TLEV_NNOW_RCF, SUCCESS, &
                               & LEAPFROG_EXPL, LEAPFROG_SI
#ifdef DEBUG
    USE mo_io_units, ONLY: nerr
#endif
    USE mo_kind, ONLY: wp
    USE mo_linked_list, ONLY: t_list_element, t_var_list
    USE mo_parallel_config, ONLY: nproma
    USE mo_util_string, ONLY: int2string
    USE mo_var_list, ONLY: nvar_lists, var_lists, get_var_timelevel
    USE mo_var_metadata_types, ONLY: t_var_metadata

    IMPLICIT NONE

    PUBLIC :: t_RestartVarData
    PUBLIC :: createRestartVarData
    PUBLIC :: getLevelPointers
    PUBLIC :: has_valid_time_level

    PRIVATE

    ! All the info that's required about a variable for restart purposes.
    TYPE t_RestartVarData
        REAL(wp), POINTER :: r_ptr(:,:,:,:,:)
        TYPE(t_var_metadata) :: info
    END TYPE t_RestartVarData

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_var_data"

CONTAINS

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
                    varIndex = varIndex + 1
                END IF
                element => element%next_list_element
            END DO
        END DO
        IF(varIndex /= varCount + 1) CALL finish(routine, "assertion failed: wrong restart variable count")
    END FUNCTION createRestartVarData

    SUBROUTINE getLevelPointers(varMetadata, varData, levelPointers)
        TYPE(t_var_metadata), INTENT(IN) :: varMetadata
        REAL(wp), TARGET :: varData(:,:,:,:,:)
        TYPE(t_ptr_2d), ALLOCATABLE, INTENT(INOUT) :: levelPointers(:)

        INTEGER :: nindex, nlevs, var_ref_pos, error, jk
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":getLevelPointers"

        ! Check if first dimension of array is nproma.
        ! Otherwise we got an array which is not suitable for this output scheme.
        IF (varMetadata%used_dimensions(1) /= nproma) CALL finish(routine,'1st dim is not nproma: '//TRIM(varMetadata%name))

        ! get data index
        nindex = 1
        IF(varMetadata%lcontained) nindex = varMetadata%ncontained

        ! get number of data levels
        nlevs = 1
        IF(varMetadata%ndims /= 2) nlevs = varMetadata%used_dimensions(2)

        var_ref_pos = varMetadata%ndims + 1
        IF(varMetadata%lcontained) var_ref_pos = varMetadata%var_ref_pos

        ! ALLOCATE the POINTER array
        IF(ALLOCATED(levelPointers)) THEN
            ! can we reuse the old allocation?
            IF(SIZE(levelPointers) /= nlevs) DEALLOCATE(levelPointers)
        END IF
        IF(.NOT.ALLOCATED(levelPointers)) THEN
            ! no suitable old allocation that can be reused
            ALLOCATE(levelPointers(nlevs), STAT = error)
            IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
        END IF

        ! get data pointers
        SELECT CASE (varMetadata%ndims)
            CASE (2)
                ! make a 3D copy of the array
                SELECT CASE(var_ref_pos)
                    CASE (1)
                        levelPointers(1)%p => varData(nindex,:,:,1,1)
                    CASE (2)
                        levelPointers(1)%p => varData(:,nindex,:,1,1)
                    CASE (3)
                        levelPointers(1)%p => varData(:,:,nindex,1,1)
                    CASE default
                        CALL finish(routine, "internal error!")
                END SELECT
            CASE (3)
                ! copy the pointer
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
                        CASE default
                            CALL finish(routine, "internal error!")
                    END SELECT
                END DO
            CASE DEFAULT
              CALL finish(routine, "'"//TRIM(varMetadata%NAME)//"': "//&
                                 & TRIM(int2string(varMetadata%ndims))//"d arrays not handled yet")
        END SELECT
    END SUBROUTINE getLevelPointers

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
        lskip_extra_timelevs = iequations == INH_ATMOSPHERE .AND. .NOT. (l_limited_area .AND. domain == 1)

        ! get time index of the given field
        time_level = get_var_timelevel(p_info)

        !TODO[NH]: I found the `time_level >= 0` condition IN the async restart code ONLY. Check whether it should be removed OR NOT.
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
                    & .AND. ha_dyn_config%itime_scheme/=LEAPFROG_SI   ) &
                    & RETURN
            CASE default
                IF ( lskip_timelev ) RETURN
        END SELECT
        resultVar = .TRUE.
#endif
    END FUNCTION has_valid_time_level

END MODULE mo_restart_var_data
