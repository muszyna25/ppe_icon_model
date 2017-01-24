!>
!! This MODULE provides an input instruction list that IS used to determine, which variables may be READ from which input file.
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_input_instructions

    USE ISO_C_BINDING, ONLY: C_CHAR
    USE mo_dictionary, ONLY: dict_get
    USE mo_exception, ONLY: message, finish
    USE mo_impl_constants, ONLY: SUCCESS, MODE_DWDANA, MODE_ICONVREMAP, MODE_IAU, MODE_IAU_OLD, MODE_COMBINED, MODE_COSMODE
    USE mo_initicon_config, ONLY: initicon_config, lread_ana, ltile_coldstart, lp2cintp_incr, lp2cintp_sfcana, &
      &                           l_sst_in
    USE mo_initicon_types, ONLY: ana_varnames_dict
    USE mo_input_request_list, ONLY: t_InputRequestList
    USE mo_lnd_nwp_config, ONLY: lsnowtile
    USE mo_model_domain, ONLY: t_patch
    USE mo_util_string, ONLY: difference, add_to_list, int2string
    USE mo_util_table, ONLY: t_table, initialize_table, add_table_column, set_table_entry, print_table, finalize_table
    USE mo_var_list, ONLY: collect_group
    USE mo_var_metadata_types,  ONLY: VARNAME_LEN

    IMPLICIT NONE

PUBLIC :: t_readInstructionList, readInstructionList_make, t_readInstructionListPtr
PUBLIC :: kInputSourceNone, kInputSourceFg, kInputSourceAna, kInputSourceBoth
PUBLIC :: kStateNoFetch, kStateFailedFetch, kStateRead
    ! The possible RETURN values of readInstructionList_sourceOfVar().
    INTEGER, PARAMETER :: kInputSourceUnset = -1, kInputSourceNone = 0, kInputSourceFg = 1, kInputSourceAna = 2, &
                        & kInputSourceBoth = 3

    ! The possible values for statusFg AND statusAna:
    INTEGER(KIND = C_CHAR), PARAMETER :: kStateNoFetch = 0, kStateFailedFetch = 1, kStateRead = 2

    ! The readInstructionList is used to tell the fetch_dwd*() routines which fields may be read from first guess and/or analysis,
    ! and to signal whether we have found data in the first guess file, so that we can fail correctly.
    ! I. e. for each variable, we expect a pair of calls to
    !
    !     wantVarXXX()
    !     handleErrorXXX()
    !
    ! The wantVarXXX() CALL tells the caller whether they should try to READ the variable,
    ! the handleErrorXXX() CALL tells the ReadInstructionList whether that READ was successfull,
    ! possibly panicking with a `finish()` CALL IN CASE there are no options left to READ the DATA.
    !
    ! If reading of the variable itself is not required, optionalReadResultXXX() has to be used
    ! instead of handleErrorXXX(), this won't panick if no data is available.
    TYPE :: t_readInstructionList
        TYPE(t_readInstruction), POINTER :: list(:)
        INTEGER :: nInstructions

    CONTAINS
        PROCEDURE :: fileRequests => readInstructionList_fileRequests   ! tell an InputRequestList what variables are required

        ! inquire whether an attempt should be made to read a variable
        PROCEDURE :: wantVar => readInstructionList_wantVar
        PROCEDURE :: wantVarFg => readInstructionList_wantVarFg
        PROCEDURE :: wantVarAna => readInstructionList_wantVarAna

        ! inform the ReadInstructionList about the result of a read, possibly triggering a `finish()` call if there is no alternative left to read the variable
        PROCEDURE :: handleError => readInstructionList_handleError
        PROCEDURE :: handleErrorFg => readInstructionList_handleErrorFg
        PROCEDURE :: handleErrorAna => readInstructionList_handleErrorAna

        ! inform the ReadInstructionList about the result of a read (nofail variant)
        PROCEDURE :: optionalReadResult => readInstructionList_optionalReadResult
        PROCEDURE :: optionalReadResultFg => readInstructionList_optionalReadResultFg
        PROCEDURE :: optionalReadResultAna => readInstructionList_optionalReadResultAna

        PROCEDURE :: sourceOfVar => readInstructionList_sourceOfVar ! returns kInputSourceNone, kInputSourceFg, kInputSourceAna, OR kInputSourceBoth
        PROCEDURE :: setSource => readInstructionList_setSource ! overrides the automatic calculation of the input source, argument must be 'kInputSource\(None\|Fg\|Ana\|Both\)'

        PROCEDURE :: fetchStatus => readInstructionList_fetchStatus  ! returns the result of a read attempt for a particular field

        ! print a table with the information obtained via `handleErrorXXX()`, `optionalReadResult()`,
        ! and `setSource()`; THIS DEPENDS ON THE ACTUAL READ ATTEMPTS AND THEIR RESULTS,
        ! NOT ON `lReadFg` or `lReadAna`.
        PROCEDURE :: printSummary => readInstructionList_printSummary

        PROCEDURE :: destruct => readInstructionList_destruct

        PROCEDURE, PRIVATE :: construct => readInstructionList_construct
        PROCEDURE, PRIVATE :: resize => readInstructionList_resize
        PROCEDURE, PRIVATE :: findInstruction => readInstructionList_findInstruction
    END TYPE t_readInstructionList

    TYPE :: t_readInstructionListPtr
        TYPE(t_readInstructionList), POINTER :: ptr
    END TYPE t_readInstructionListPtr

PRIVATE

    TYPE :: t_readInstruction
        CHARACTER(LEN = VARNAME_LEN) :: varName
        LOGICAL :: lReadFg, lReadAna    ! These reflect the result of interpreting the variable groups, they are not used in the table output.
        LOGICAL :: lRequireAna  ! Panick if reading from analysis fails (set according to the `ana_varlist` namelist parameter)
        INTEGER(KIND = C_CHAR) :: statusFg, statusAna   ! These reflect what read attempts really have been made, and what their results were. This is the basis for the table output.
        INTEGER :: sourceOverride   ! IF this IS NOT kInputSourceUnset, it overrides the automatic calculation of the input source.
    CONTAINS
        PROCEDURE :: source => readInstruction_source
    END TYPE t_readInstruction

    CHARACTER(LEN = *), PARAMETER :: modname = 'mo_input_instructions'

CONTAINS

    SUBROUTINE collectGroupFg(outGroup, outGroupSize, init_mode)
        CHARACTER(LEN = VARNAME_LEN), INTENT(INOUT) :: outGroup(:)
        INTEGER, INTENT(OUT) :: outGroupSize
        INTEGER, VALUE :: init_mode

        SELECT CASE(init_mode)
            CASE(MODE_DWDANA, MODE_ICONVREMAP)
                CALL collect_group('mode_dwd_fg_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE(MODE_IAU)
                CALL collect_group('mode_iau_fg_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE(MODE_IAU_OLD)
                CALL collect_group('mode_iau_old_fg_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE(MODE_COMBINED)
                CALL collect_group('mode_combined_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE(MODE_COSMODE)
                CALL collect_group('mode_cosmode_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE DEFAULT
                outGroupSize = 0
        END SELECT
    END SUBROUTINE collectGroupFg

    SUBROUTINE collectGroupAna(outGroup, outGroupSize, init_mode)
        CHARACTER(LEN = VARNAME_LEN), INTENT(INOUT) :: outGroup(:)
        INTEGER, INTENT(OUT) :: outGroupSize
        INTEGER, VALUE :: init_mode

        SELECT CASE(init_mode)
            CASE(MODE_DWDANA, MODE_ICONVREMAP)
                CALL collect_group('mode_dwd_ana_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE(MODE_IAU)
                CALL collect_group('mode_iau_ana_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE(MODE_IAU_OLD)
                CALL collect_group('mode_iau_old_ana_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE DEFAULT
                outGroupSize = 0
        END SELECT
    END SUBROUTINE collectGroupAna

    SUBROUTINE collectGroupAnaAtm(outGroup, outGroupSize, init_mode)
        CHARACTER(LEN = VARNAME_LEN), INTENT(INOUT) :: outGroup(:)
        INTEGER, INTENT(OUT) :: outGroupSize
        INTEGER, VALUE :: init_mode

        SELECT CASE(init_mode)
            CASE(MODE_IAU, MODE_IAU_OLD)
                CALL collect_group('mode_iau_anaatm_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE DEFAULT
                outGroupSize = 0
        END SELECT
    END SUBROUTINE collectGroupAnaAtm

    SUBROUTINE copyGroup(inGroup, inGroupSize, outGroup, outGroupSize)
        CHARACTER(LEN = VARNAME_LEN), INTENT(IN) :: inGroup(:)
        INTEGER, VALUE :: inGroupSize
        CHARACTER(LEN = VARNAME_LEN), INTENT(INOUT) :: outGroup(:)
        INTEGER, INTENT(OUT) :: outGroupSize

        outGroup(1:inGroupSize) = inGroup(1:inGroupSize)
        outGroupSize = inGroupSize
    END SUBROUTINE copyGroup

    SUBROUTINE mergeAnaIntoFg(anaGroup, anaGroupSize, fgGroup, fgGroupSize)
        CHARACTER(LEN = VARNAME_LEN), INTENT(INOUT) :: anaGroup(:), fgGroup(:)
        INTEGER, INTENT(INOUT) :: anaGroupSize, fgGroupSize

        ! fgGroup += anaGroup
        CALL add_to_list(fgGroup, fgGroupSize, anaGroup(1:anaGroupSize), anaGroupSize)

        ! Remove fields 'u', 'v', 'temp', 'pres'
        CALL difference(fgGroup, fgGroupSize, (/'u   ','v   ','temp','pres'/), 4)

        ! anaGroup = --
        anaGroupSize = 0
    END SUBROUTINE mergeAnaIntoFg

    SUBROUTINE collectGroups(p_patch, init_mode, fgGroup, fgGroupSize, anaGroup, anaGroupSize)
        TYPE(t_patch), INTENT(IN) :: p_patch
        INTEGER, VALUE :: init_mode
        CHARACTER(LEN = VARNAME_LEN), INTENT(INOUT) :: anaGroup(:), fgGroup(:)
        INTEGER, INTENT(OUT) :: anaGroupSize, fgGroupSize

        CHARACTER(LEN = *), PARAMETER :: routine = modname//':collectGroups'
        CHARACTER(LEN=VARNAME_LEN), DIMENSION(200) :: anaAtmGroup
        INTEGER :: anaAtmGroupSize, jg
        LOGICAL :: lRemoveSnowfrac
        CHARACTER(LEN = 256) :: message_text

        ! get the raw DATA
        CALL collectGroupFg(fgGroup, fgGroupSize, init_mode)
        CALL collectGroupAna(anaGroup, anaGroupSize, init_mode)
        CALL collectGroupAnaAtm(anaAtmGroup, anaAtmGroupSize, init_mode)
        jg = p_patch%id


        ! integrate the information of these three groups, the init_mode, AND some flags to produce the effective fgGroup AND anaGroup
        SELECT CASE(init_mode)
            CASE(MODE_DWDANA, MODE_ICONVREMAP)
                IF(.NOT.lread_ana) THEN
                    ! lump together fgGroup and anaGroup
                    CALL mergeAnaIntoFg(anaGroup, anaGroupSize, fgGroup, fgGroupSize)
                ENDIF

            CASE(MODE_IAU, MODE_IAU_OLD)
                ! in case of tile coldstart, we can omit snowfrac_lc
                ! Remove field 'snowfrac_lc' from FG list
                lRemoveSnowfrac = ltile_coldstart
                IF(init_mode == MODE_IAU .AND. .NOT. lsnowtile) lRemoveSnowfrac = .TRUE.
                IF(lRemoveSnowfrac) CALL difference(fgGroup, fgGroupSize, (/'snowfrac_lc'/), 1)

                IF (.NOT. lp2cintp_incr(jg) .AND. .NOT. lp2cintp_sfcana(jg) ) THEN
                    ! full ANA read

                ELSE IF (lp2cintp_incr(jg) .AND. .NOT. lp2cintp_sfcana(jg)) THEN
                    ! SFC-ANA read
                    ! atmospheric analysis fieds are interpolated from parent domain, 
                    ! however surface analysis fields are read from file

                    ! Remove fields atmospheric analysis fields from anaGroup
                    CALL difference(anaGroup, anaGroupSize, anaAtmGroup, anaAtmGroupSize)

                ELSE IF (lp2cintp_incr(jg) .AND. lp2cintp_sfcana(jg) ) THEN
                    ! no ANA-read
                    ! lump together fgGroup and anaGroup
                    CALL mergeAnaIntoFg(anaGroup, anaGroupSize, fgGroup, fgGroupSize)

                ELSE
                    WRITE(message_text,'(a,l1,a,l1,a)') 'Combination lp2cintp_incr=',lp2cintp_incr(jg), &
                    &                       ' and lp2cintp_sfcana=',lp2cintp_sfcana(jg),' not allowed'
                    CALL finish(routine, TRIM(message_text))
                ENDIF

            CASE(MODE_COMBINED,MODE_COSMODE)
                ! remove W_SO from default list and replace it by SMI
                CALL difference (fgGroup, fgGroupSize, (/'w_so'/), 1)
                CALL add_to_list(fgGroup, fgGroupSize, (/'smi'/) , 1)

                ! no analysis group
                anaGroupSize = 0

            CASE DEFAULT
                fgGroupSize = 0
                anaGroupSize = 0

        END SELECT

    END SUBROUTINE collectGroups


    !-------------
    !>
    !! SUBROUTINE readInstructionList_make
    !! Generates a mapping from icon variable names to a pair of logicals that define
    !! from which file each variable IS allowed to be READ.
    !!
    !! This IS based on:
    !!  1. the group specifications IN the add_var() calls
    !!  2. user input via initicon_config(:)%ana_varlist,
    !!     ANY of these variables must be READ from the analysis file
    !!  3. some other interesting rules, like the substitution of 'smi' for 'w_so'
    !!     IN the CASE of MODE_COSMODE AND MODE_COMBINED.
    !!
    !! These are the variable groups which are used to determine from which file a variable IS READ:
    !!     MODE_DWDANA    : mode_dwd_fg_in + mode_dwd_ana_in
    !!     MODE_ICONVREMAP: mode_dwd_fg_in + mode_dwd_ana_in
    !!     MODE_IAU       : mode_iau_fg_in + mode_iau_ana_in - mode_iau_anaatm_in
    !!     MODE_IAU_OLD   : mode_iau_old_fg_in + mode_iau_old_ana_in - mode_iau_anaatm_in
    !!     MODE_COMBINED  : mode_combined_in
    !!     MODE_COSMODE   : mode_cosmode_in
    !!     MODE_IFSANA    : <NONE>
    !!
    !! In contrast to the old create_input_groups() SUBROUTINE, this does NOT check the contents of the files itself.
    !! It ONLY provides concise instructions of which fields to READ from which file.
    !!
    !! Returns a new instruction list object.
    FUNCTION readInstructionList_make(p_patch, init_mode) RESULT(resultVar)
        TYPE(t_patch), INTENT(IN) :: p_patch
        INTEGER, VALUE :: init_mode
        ! the resulting list of variable names to be READ together with flags defining which input file may be used
        TYPE(t_readInstructionList), POINTER :: resultVar

        ! local variables
        CHARACTER(LEN = *), PARAMETER :: routine = modname//':readInstructionList_make'
        INTEGER :: ivar, error
        TYPE(t_readInstruction), pointer :: curInstruction

        ! lists of variable names
        CHARACTER(LEN=VARNAME_LEN), DIMENSION(200) :: grp_vars_fg, grp_vars_ana

        ! the corresponding sizes of the lists above
        INTEGER :: ngrp_vars_fg, ngrp_vars_ana

        !-------------------------------------------------------------------------


        ! create a list of instructions according to the current init_mode AND configuration flags
        CALL collectGroups(p_patch, init_mode, grp_vars_fg, ngrp_vars_fg, grp_vars_ana, ngrp_vars_ana)
        ALLOCATE(resultVar, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure");
        CALL resultVar%construct()
        DO ivar = 1, ngrp_vars_fg
            curInstruction => resultVar%findInstruction(grp_vars_fg(ivar))
            curInstruction%lReadFg = .TRUE.
        END DO
        DO ivar = 1, ngrp_vars_ana
            curInstruction => resultVar%findInstruction(grp_vars_ana(ivar))
            curInstruction%lReadFg = .TRUE.
            curInstruction%lReadAna = .TRUE.
        END DO


        ! Allow the user to override the DEFAULT settings via ana_varlist. These variables must be READ from analysis.
        IF( lread_ana ) THEN
            ! translate GRIB2 varname to internal netcdf varname
            ! If requested GRIB2 varname is not found in the dictionary
            ! (i.e. due to typos) -> Model abort
            DO ivar=1,SIZE(initicon_config(p_patch%id)%ana_varlist)
                IF (initicon_config(p_patch%id)%ana_varlist(ivar) == ' ') EXIT

                curInstruction => resultVar%findInstruction(TRIM(dict_get(ana_varnames_dict, &
                &                                                      initicon_config(p_patch%id)%ana_varlist(ivar), &
                &                                                      linverse=.TRUE.)))
                curInstruction%lReadAna = .TRUE.
                curInstruction%lRequireAna = .TRUE.
            ENDDO
        END IF

    END FUNCTION readInstructionList_make

    SUBROUTINE readInstructionList_construct(me)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        INTEGER :: error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":readInstructionList_construct"

        ALLOCATE(me%list(8), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        me%nInstructions = 0
    END SUBROUTINE readInstructionList_construct

    FUNCTION readInstructionList_findInstruction(me, varName) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName
        TYPE(t_readInstruction), POINTER :: resultVar

        INTEGER :: i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":readInstructionList_findInstruction"

        ! try to find it IN the current list
        DO i = 1, me%nInstructions
            resultVar => me%list(i)
            IF(resultVar%varName == varName) RETURN
        END DO

        ! the variable IS NOT found IN the list, expand the list
        ! first increase the allocation IF necessary
        IF(me%nInstructions == SIZE(me%list, 1)) CALL me%resize(2*me%nInstructions)
        ! THEN add a new entry at the END
        me%nInstructions = me%nInstructions + 1
        resultVar => me%list(me%nInstructions)
        resultVar%varName = varName
        resultVar%lReadFg = .FALSE.
        resultVar%lReadAna = .FALSE.
        resultVar%lRequireAna = .FALSE.
        resultVar%statusFg = kStateNoFetch
        resultVar%statusAna = kStateNoFetch
        resultVar%sourceOverride = kInputSourceUnset
    END FUNCTION readInstructionList_findInstruction

    SUBROUTINE readInstructionList_resize(me, newSize)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        INTEGER, VALUE :: newSize

        TYPE(t_readInstruction), POINTER :: tempList(:)
        INTEGER :: i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":readInstructionList_resize"

        ! just IN CASE we are actually shrinking the list, avoid overrunning the new buffer
        IF(newSize < me%nInstructions) me%nInstructions = newSize

        ! get a new allocation
        ALLOCATE(tempList(newSize), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ! copy the list contents over
        DO i = 1, me%nInstructions
            tempList(i)%varName = me%list(i)%varName
            tempList(i)%lReadFg = me%list(i)%lReadFg
            tempList(i)%lReadAna = me%list(i)%lReadAna
            tempList(i)%lRequireAna = me%list(i)%lRequireAna
            tempList(i)%statusFg = me%list(i)%statusFg
            tempList(i)%statusAna = me%list(i)%statusAna
            tempList(i)%sourceOverride = me%list(i)%sourceOverride
        END DO

        ! replace the old list with the new one
        DEALLOCATE(me%list)
        me%list => tempList
        tempList => NULL()
    END SUBROUTINE readInstructionList_resize

    SUBROUTINE readInstructionList_fileRequests(me, requestList, lIsFg)
        CLASS(t_readInstructionList), INTENT(IN) :: me
        TYPE(t_inputRequestList), INTENT(INOUT) :: requestList
        LOGICAL, VALUE :: lIsFg

        INTEGER :: i

        DO i = 1, me%nInstructions
            IF(lIsFg) THEN
                IF(me%list(i)%lReadFg) CALL requestList%request(TRIM(me%list(i)%varName))
            ELSE
                IF(me%list(i)%lReadAna) CALL requestList%request(TRIM(me%list(i)%varName))
            END IF
        END DO
    END SUBROUTINE readInstructionList_fileRequests

    LOGICAL FUNCTION readInstructionList_wantVar(me, varName, lIsFg) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(IN) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName
        LOGICAL, VALUE :: lIsFg

        IF(lIsFg) THEN
            resultVar = me%wantVarFg(varName)
        ELSE
            resultVar = me%wantVarAna(varName)
        END IF
    END FUNCTION readInstructionList_wantVar

    LOGICAL FUNCTION readInstructionList_wantVarFg(me, varName) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(IN) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName

        INTEGER :: i

        resultVar = .FALSE.
        DO i = 1, me%nInstructions
            IF(me%list(i)%varName == varName) THEN
                resultVar = me%list(i)%lReadFg
                RETURN
            END IF
        END DO
    END FUNCTION readInstructionList_wantVarFg

    LOGICAL FUNCTION readInstructionList_wantVarAna(me, varName) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(IN) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName

        INTEGER :: i

        resultVar = .FALSE.
        DO i = 1, me%nInstructions
            IF(me%list(i)%varName == varName) THEN
                resultVar = me%list(i)%lReadAna
                RETURN
            END IF
        END DO
    END FUNCTION readInstructionList_wantVarAna

    SUBROUTINE readInstructionList_handleError(me, lSuccess, varName, caller, lIsFg)
        LOGICAL, VALUE :: lSuccess, lIsFg
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        IF(lIsFg) THEN
            CALL me%handleErrorFg(lSuccess, varName, caller)
        ELSE
            CALL me%handleErrorAna(lSuccess, varName, caller)
        END IF
    END SUBROUTINE readInstructionList_handleError

    SUBROUTINE readInstructionList_handleErrorFg(me, lSuccess, varName, caller)
        LOGICAL, VALUE :: lSuccess
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        CHARACTER(LEN = *), PARAMETER :: routine = modname//"readInstructionList_handleErrorFg"
        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)

        ! sanity check
        IF(.NOT.instruction%lReadFg) CALL finish(routine, "internal error: variable '"//varName//"' was read even though there &
            &are no read instructions for it")

        IF(lSuccess) THEN
            instruction%statusFg = kStateRead
        ELSE
            instruction%statusFg = kStateFailedFetch

            ! check whether we can READ this variable from the analysis file as well
            IF(instruction%lReadAna) THEN
                ! now this variable must be READ from the analysis file
                instruction%lReadFg = .FALSE.
            ELSE
                CALL finish(caller, "failed to read variable '"//varName//"' from the first guess file, &
                &and reading from analysis file is not allowed")
            END IF
        END IF
    END SUBROUTINE readInstructionList_handleErrorFg

    SUBROUTINE readInstructionList_handleErrorAna(me, lSuccess, varName, caller)
        LOGICAL, VALUE :: lSuccess
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        CHARACTER(LEN = *), PARAMETER :: routine = modname//"readInstructionList_handleErrorAna"
        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)

        ! sanity check
        IF(.NOT.instruction%lReadAna) CALL finish(routine, "internal error: variable '"//varName//"' was read even though there &
            &are no read instructions for it")

        IF(lSuccess) THEN
            instruction%statusAna = kStateRead
        ELSE
            instruction%statusAna = kStateFailedFetch

            ! check whether this variable was already READ from the first guess file; IF NOT, we have to error OUT
            IF(.NOT.instruction%lReadFg) THEN
                CALL finish(caller, "failed to read variable '"//varName//"' from the analysis file, &
                    &and reading from first guess is not allowed or failed")
            END IF

            ! check whether we are forced by the user to READ analysis DATA for this variable (ana_varlist)
            IF(instruction%lRequireAna) THEN
                CALL finish(caller, "failed to read variable '"//varName//"' from the analysis file (required by ana_varlist)")
            END IF
        END IF
    END SUBROUTINE readInstructionList_handleErrorAna

    SUBROUTINE readInstructionList_optionalReadResult(me, lSuccess, varName, caller, lIsFg)
        LOGICAL, VALUE :: lSuccess, lIsFg
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        IF(lIsFg) THEN
            CALL me%optionalReadResultFg(lSuccess, varName, caller)
        ELSE
            CALL me%optionalReadResultAna(lSuccess, varName, caller)
        END IF
    END SUBROUTINE readInstructionList_optionalReadResult

    SUBROUTINE readInstructionList_optionalReadResultFg(me, lSuccess, varName, caller)
        LOGICAL, VALUE :: lSuccess
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        CHARACTER(LEN = *), PARAMETER :: routine = modname//"readInstructionList_optionalReadResultFg"
        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)

        ! sanity check
        IF(.NOT.instruction%lReadFg) CALL finish(routine, "internal error: variable '"//varName//"' was read even though there &
            &are no read instructions for it")

        IF(lSuccess) THEN
            instruction%statusFg = kStateRead
        ELSE
            instruction%statusFg = kStateFailedFetch

            ! check whether we can READ this variable from the analysis file as well
            IF(instruction%lReadAna) THEN
                ! now this variable must be READ from the analysis file
                instruction%lReadFg = .FALSE.
            END IF
        END IF
    END SUBROUTINE readInstructionList_optionalReadResultFg

    SUBROUTINE readInstructionList_optionalReadResultAna(me, lSuccess, varName, caller)
        LOGICAL, VALUE :: lSuccess
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        CHARACTER(LEN = *), PARAMETER :: routine = modname//"readInstructionList_optionalReadResultAna"
        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)

        ! sanity check
        IF(.NOT.instruction%lReadAna) CALL finish(routine, "internal error: variable '"//varName//"' was read even though there &
            &are no read instructions for it")

        IF(lSuccess) THEN
            instruction%statusAna = kStateRead
        ELSE
            instruction%statusAna = kStateFailedFetch

            ! check whether we are forced by the user to READ analysis DATA for this variable (ana_varlist)
            IF(instruction%lRequireAna) THEN
                CALL finish(caller, "failed to read variable '"//varName//"' from the analysis file (required by ana_varlist)")
            END IF
        END IF
    END SUBROUTINE readInstructionList_optionalReadResultAna

    INTEGER FUNCTION readInstruction_source(me) RESULT(resultVar)
        CLASS(t_readInstruction), INTENT(IN) :: me

        IF(me%sourceOverride /= kInputSourceUnset) THEN
            resultVar = me%sourceOverride
        ELSE
            IF(me%statusAna == kStateRead) THEN
                resultVar = kInputSourceAna
            ELSE IF(me%statusFg == kStateRead) THEN
                resultVar = kInputSourceFg
            ELSE
                resultVar = kInputSourceNone
            END IF
        END IF
    END FUNCTION readInstruction_source

    INTEGER FUNCTION readInstructionList_sourceOfVar(me, varName) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName

        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)
        resultVar = instruction%source()
    END FUNCTION readInstructionList_sourceOfVar

    SUBROUTINE readInstructionList_setSource(me, varName, source)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName
        INTEGER, VALUE :: source

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":readInstructionList_setSource"
        TYPE(t_readInstruction), POINTER :: instruction

        SELECT CASE(source)
            CASE(kInputSourceNone, kInputSourceFg, kInputSourceAna, kInputSourceBoth)
                instruction => me%findInstruction(varName)
                instruction%sourceOverride = source
            CASE DEFAULT
                CALL finish(routine, "assertion failed: illegal source argument")
        END SELECT
    END SUBROUTINE readInstructionList_setSource


    INTEGER FUNCTION readInstructionList_fetchStatus(me, varName, lIsFg) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName
        LOGICAL, VALUE :: lIsFg

        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)
        IF(lIsFg) THEN
            resultVar = instruction%statusFg
        ELSE
            resultVar = instruction%statusAna
        ENDIF
    END FUNCTION readInstructionList_fetchStatus




    ! The table that is printed by this function deliberately depends on the actual read attempts and their results, not on `lReadFg` or `lReadAna`.
    ! This is due to the fact that there are existing discrepancies between the input groups and the actual read attempts made by the `fetch...()` routines
    ! in `mo_initicon_io`: The table is supposed to show the reality of which data was read from where, and which inputs were used,
    ! not some hypothetical this-is-what-should-have-been-done info.
    SUBROUTINE readInstructionList_printSummary(me, domain)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        INTEGER, VALUE :: domain

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":readInstructionList_printSummary"
        CHARACTER(LEN = 12) :: fgString, anaString
        INTEGER :: i
        TYPE(t_readInstruction), POINTER :: curInstruction
        TYPE(t_table) :: table
        CHARACTER(LEN = *), PARAMETER :: variableCol = "variable", &
                                       & fgAttemptCol = "FG read attempt", &
                                       & fgSuccessCol = "FG data found", &
                                       & anaAttemptCol = "ANA read attempt", &
                                       & anaSuccessCol = "ANA data found", &
                                       & anaRequiredCol = "ANA DATA required (ana_varlist)", &
                                       & useCol = "data used from"

        IF(me%nInstructions == 0) THEN
            ! Don't print empty tables, just give a message that there are no input instructions.
            WRITE(0, *) "no input results available for domain "//TRIM(int2string(domain))
            RETURN
        END IF

        ! Print the title of the list.
        WRITE(0,*) ""
        WRITE(0,*) "input results for domain "//TRIM(int2string(domain))//":"

        CALL initialize_table(table)
        CALL add_table_column(table, variableCol)
        CALL add_table_column(table, fgAttemptCol)
        CALL add_table_column(table, fgSuccessCol)
        CALL add_table_column(table, anaAttemptCol)
        CALL add_table_column(table, anaSuccessCol)
        CALL add_table_column(table, anaRequiredCol)
        CALL add_table_column(table, useCol)

        DO i = 1, me%nInstructions
            curInstruction => me%list(i)

            CALL set_table_entry(table, i, variableCol, TRIM(curInstruction%varName))

            ! Display whether an actual read attempt from FG has been made, and its result.
            ! read attempt == there was a call from a `fetch...()` routine in `mo_initicon_io`
            SELECT CASE(curInstruction%statusFg)
                CASE(kStateNoFetch)
                    CALL set_table_entry(table, i, fgAttemptCol, "no")
                CASE(kStateFailedFetch)
                    CALL set_table_entry(table, i, fgAttemptCol, "yes")
                    CALL set_table_entry(table, i, fgSuccessCol, "no")
                CASE(kStateRead)
                    CALL set_table_entry(table, i, fgAttemptCol, "yes")
                    CALL set_table_entry(table, i, fgSuccessCol, "yes")
                CASE DEFAULT
                    CALL finish(routine, "unexpected value in statusFg")
            END SELECT

            ! As above, but for ANA.
            SELECT CASE(curInstruction%statusAna)
                CASE(kStateNoFetch)
                    CALL set_table_entry(table, i, anaAttemptCol, "no")
                CASE(kStateFailedFetch)
                    CALL set_table_entry(table, i, anaAttemptCol, "yes")
                    CALL set_table_entry(table, i, anaSuccessCol, "no")
                CASE(kStateRead)
                    CALL set_table_entry(table, i, anaAttemptCol, "yes")
                    CALL set_table_entry(table, i, anaSuccessCol, "yes")
                CASE DEFAULT
                    CALL finish(routine, "unexpected value in statusAna")
            END SELECT

            IF(curInstruction%lRequireAna) THEN
                CALL set_table_entry(table, i, anaRequiredCol, "required")
            END IF

            ! Display which data is actually used. If there was a call to `setSource()`, that is what is printed,
            ! otherwise this column is computed from the two `statusXXX` variables printed above.
            SELECT CASE(curInstruction%source())
                CASE(kInputSourceNone)
                    CALL set_table_entry(table, i, useCol, "none")
                CASE(kInputSourceFg)
                    CALL set_table_entry(table, i, useCol, "fg")
                CASE(kInputSourceAna)
                    CALL set_table_entry(table, i, useCol, "ana")
                CASE(kInputSourceBoth)
                    CALL set_table_entry(table, i, useCol, "both")
                CASE DEFAULT
                    CALL finish(routine, "unexpected RESULT from curInstruction%source()")
            END SELECT
        END DO

        CALL print_table(table, opt_delimiter = " | ")

        CALL finalize_table(table)

        WRITE(0,*) ""   ! Separate the table from the following messages.
    END SUBROUTINE readInstructionList_printSummary

    SUBROUTINE readInstructionList_destruct(me)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me

        DEALLOCATE(me%list)
    END SUBROUTINE readInstructionList_destruct

END MODULE mo_input_instructions
