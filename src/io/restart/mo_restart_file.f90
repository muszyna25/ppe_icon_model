!>
!! A class to open, write, and close restart files.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_restart_file
    USE mo_cdi, ONLY: CDI_UNDEFID, FILETYPE_NC2, FILETYPE_NC4, streamWriteVarSlice
    USE mo_exception, ONLY: finish
    USE mo_io_units, ONLY: filename_max
    USE mo_kind, ONLY: dp
    USE mo_linked_list, ONLY: t_var_list
    USE mo_restart_attributes, ONLY: t_RestartAttributeList
    USE mo_restart_patch_description, ONLY: t_restart_patch_description
    USE mo_restart_util, ONLY: getRestartFilename, t_restart_cdi_ids, t_var_data, t_restart_args, has_valid_time_level
    USE mo_var_list, ONLY: nvar_lists, var_lists

    IMPLICIT NONE

    PUBLIC :: t_RestartFile

    PRIVATE

    TYPE t_RestartFile
        CHARACTER(LEN=filename_max) :: filename
        CHARACTER(LEN=32) :: model_type
        TYPE(t_restart_cdi_ids) :: cdiIds
    CONTAINS
        PROCEDURE :: open => restartFile_open
        PROCEDURE :: writeLevel => restartFile_writeLevel
        PROCEDURE :: close => restartFile_close
    END TYPE t_RestartFile

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_file"

CONTAINS

    SUBROUTINE restartFile_open(me, description, varData, restart_args, restartAttributes)
        CLASS(t_RestartFile), INTENT(INOUT) :: me
        TYPE(t_restart_patch_description), INTENT(IN) :: description
        TYPE(t_var_data), INTENT(INOUT) :: varData(:)
        TYPE(t_restart_args), INTENT(IN) :: restart_args
        TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes

        TYPE(t_var_list), POINTER     :: p_re_list
        INTEGER                       :: restart_type, i
        CHARACTER(LEN=*), PARAMETER   :: routine = modname//':restartFile_open'

#ifdef DEBUG
        WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif
        CALL me%cdiIds%init()

        ! get the first restart list from the global lists
        p_re_list => NULL()
        DO i = 1, nvar_lists
            ! skip, if var_list is not required for restart
            IF (.NOT. var_lists(i)%p%lrestart) CYCLE

            ! save pointer
            p_re_list => var_lists(i)
            EXIT
        ENDDO
        ! assume all restart variables have the same model name
        me%model_type = p_re_list%p%model_type

        ! assume all restart variables uses the same file format
        SELECT CASE (p_re_list%p%restart_type)
            CASE (FILETYPE_NC2, FILETYPE_NC4)
                restart_type = p_re_list%p%restart_type
            CASE default
                CALL finish(routine, "file format for restart variables must be NetCDF")
        END SELECT

        me%filename = getRestartFilename(description%base_filename, description%id, restart_args%datetime, me%model_type)

        IF(ALLOCATED(description%opt_pvct)) THEN
            CALL me%cdiIds%openRestartAndCreateIds(TRIM(me%filename), restart_type, restartAttributes, &
                                                  &description%n_patch_cells_g, description%n_patch_verts_g, &
                                                  &description%n_patch_edges_g, description%cell_type, &
                                                  &description%v_grid_defs(1:description%v_grid_count), description%opt_pvct)
        ELSE
            CALL me%cdiIds%openRestartAndCreateIds(TRIM(me%filename), restart_type, restartAttributes, &
                                                  &description%n_patch_cells_g, description%n_patch_verts_g, &
                                                  &description%n_patch_edges_g, description%cell_type, &
                                                  &description%v_grid_defs(1:description%v_grid_count))
        END IF

#ifdef DEBUG
        WRITE (nerr, FORMAT_VALS5)routine,' p_pe=',p_pe,' open netCDF file with ID=',me%cdiIds%file
#endif

        ! go over the all restart variables in the associated array AND define those that have a valid time level
        DO i = 1, SIZE(varData)
            IF(has_valid_time_level(varData(i)%info, description%id, description%nnew, description%nnew_rcf)) THEN
                CALL me%cdiIds%defineVariable(varData(i)%info)
            END IF
        ENDDO

        CALL me%cdiIds%finalizeVlist(restart_args%datetime)
    END SUBROUTINE restartFile_open

    SUBROUTINE restartFile_writeLevel(me, varId, levelId, DATA)
        CLASS(t_RestartFile), INTENT(IN) :: me
        INTEGER, VALUE :: varId, levelId
        REAL(dp), INTENT(IN) :: DATA(:)

        CALL streamWriteVarSlice(me%cdiIds%file, varId, levelId, DATA, 0)
    END SUBROUTINE restartFile_writeLevel

    !------------------------------------------------------------------------------------------------
    !
    ! Closes the given restart file.
    !
    SUBROUTINE restartFile_close(me)
        CLASS(t_RestartFile), INTENT(INOUT) :: me
        CHARACTER(LEN = *), PARAMETER :: routine = modname//':restartFile_close'

#ifdef DEBUG
        WRITE (nerr,FORMAT_VALS3)routine,' p_pe=',p_pe
#endif

        IF (me%cdiIds%file /= CDI_UNDEFID) THEN
#ifdef DEBUG
            WRITE (nerr,'(3a)')routine,' try to close restart file=',TRIM(me%filename)
            WRITE (nerr, FORMAT_VALS5)routine,' p_pe=',p_pe,' close netCDF file with ID=',me%cdiIds%file
#endif
        ENDIF

        CALL me%cdiIds%closeAndDestroyIds()
    END SUBROUTINE restartFile_close

END MODULE mo_restart_file
