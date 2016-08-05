!> Module for reading restart files
!!
!! Note: The asynchronous implementation of the restart output can be
!!       found in the module "mo_async_restart"
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_load_restart
    USE mo_cdi, ONLY: streamOpenRead, streamInqVlist, streamClose, streamOpenRead, streamReadVarSlice, vlistInqTaxis, vlistNvars, &
                    & vlistInqVarName, vlistInqVarGrid, vlistInqVarZaxis, taxisInqVdate, taxisInqVtime, zaxisInqType, &
                    & zaxisInqSize, gridInqSize, ZAXIS_SURFACE
    USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE
    USE mo_communication, ONLY: t_ScatterPattern
    USE mo_exception, ONLY: message, finish
    USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH
    USE mo_kind, ONLY: wp
    USE mo_linked_list, ONLY: t_list_element
    USE mo_model_domain, ONLY: t_patch
    USE mo_mpi, ONLY: p_comm_work, p_comm_rank, p_bcast, my_process_is_mpi_workroot
    USE mo_restart_attributes, ONLY: t_RestartAttributeList, RestartAttributeList_make, setAttributesForRestarting
    USE mo_restart_namelist, ONLY: read_and_bcast_restart_namelists
    USE mo_util_cdi, ONLY: cdiGetStringError
    USE mo_util_hash, ONLY: util_hashword
    USE mo_util_string, ONLY: int2string, separator
    USE mo_var_list, ONLY: nvar_lists, var_lists, find_list_element
    USE mo_var_metadata_types, ONLY: t_var_metadata
    IMPLICIT NONE

    PUBLIC :: read_restart_files, read_restart_header

    CHARACTER(LEN = *), PARAMETER :: modname = "mo_load_restart"

CONTAINS

  !------------------------------------------------------------------------------------------------
  !> Reads attributes and namelists for all available domains from restart file.
  !>
  !> XXX: The code that I found READ the restart attributes from all the files,
  !>      *appending* them to the list of attributes. Since all restart files contain the same attributes,
  !>      this ONLY served to duplicate them. The current code just ignores the attributes IN all but the first file,
  !>      just like the original code ignored the namelists IN all but the first file.
  !>      However, it might be a good idea to add some consistency checking on the attributes/namelists of the other files.
  SUBROUTINE read_restart_header(modeltype_str)
    CHARACTER(LEN=*), INTENT(IN)  :: modeltype_str
    ! local variables
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::read_restart_header"
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: rst_filename
    LOGICAL                        :: lsuccess, lexists
    INTEGER                        :: idom, total_dom, fileID, vlistID, myRank
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText
    INTEGER, PARAMETER :: root_pe = 0
    TYPE(t_RestartAttributeList), SAVE, POINTER :: restartAttributes

    ! rank of broadcast root PE
    myRank = p_comm_rank(p_comm_work)

    idom = 1
    rst_filename = "restart_"//TRIM(modeltype_str)//"_DOM"//TRIM(int2string(idom, "(i2.2)"))//".nc"

    ! test if the domain-dependent restart file exists:
    INQUIRE(file=TRIM(rst_filename), exist=lexists)
    ! otherwise, give a warning and resort to the old naming scheme:
    IF (.NOT. lexists) THEN
        CALL finish(routine, "Restart file not found! Expected name: restart_<model_type>_DOM<2 digit domain number>.nc")
    END IF

    ! Read all namelists used in the previous run
    ! and store them in a buffer. These values will overwrite the
    ! model default, and will later be overwritten if the user has
    ! specified something different for this integration.

    ! Note: We read the namelists AND attributes only once and assume that these
    !       are identical for all domains (which IS guaranteed by the way the restart files are written).

    IF (myRank == root_pe) THEN
        fileID  = streamOpenRead(TRIM(rst_filename))
        ! check if the file could be opened
        IF (fileID < 0) THEN
            CALL cdiGetStringError(fileID, cdiErrorText)
            CALL finish(routine, 'File '//TRIM(rst_filename)//' cannot be opened: '//TRIM(cdiErrorText))
        ENDIF

        vlistID = streamInqVlist(fileID)
    END IF
    CALL read_and_bcast_restart_namelists(vlistID, root_pe, p_comm_work)
    restartAttributes => RestartAttributeList_make(vlistID, root_pe, p_comm_work)
    CALL setAttributesForRestarting(restartAttributes)
    IF (myRank == root_pe) CALL streamClose(fileID)

    CALL message(TRIM(routine), 'read namelists AND attributes from restart file')

    ! since we do not know about the total number of domains yet,
    ! we have to ask the restart file for this information:
    total_dom = restartAttributes%getInteger( 'n_dom' )

    ! check whether we have all the restart files we expect
    IF (myRank == root_pe) THEN
        DO idom = 2, total_dom
            rst_filename = "restart_"//TRIM(modeltype_str)//"_DOM"//TRIM(int2string(idom, "(i2.2)"))//".nc"
            IF (idom > 1) INQUIRE(file=TRIM(rst_filename), exist=lexists)
            IF (lexists) THEN
                CALL message(TRIM(routine), 'read global attributes from restart file')
            ELSE
                CALL message(TRIM(routine), 'Warning: domain not active at restart time')
            ENDIF
        END DO
    END IF

  END SUBROUTINE read_restart_header

  SUBROUTINE read_restart_files(p_patch, opt_ndom)

    TYPE(t_patch), INTENT(in) :: p_patch
    INTEGER,       OPTIONAL, INTENT(in) :: opt_ndom
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::read_restart_files"

    TYPE model_search
      CHARACTER(len=8) :: abbreviation
      INTEGER          :: key
    END type model_search
    TYPE(model_search) :: abbreviations(nvar_lists)

    TYPE (t_list_element), POINTER   :: element
    TYPE (t_var_metadata), POINTER   :: info
    CHARACTER(len=80)                :: restart_filename, name
    INTEGER                          :: fileID, vlistID, gridID, zaxisID, taxisID,    &
      &                                 varID, idate, itime, ic, il, n, nfiles, i,   &
      &                                 iret, istat, key, vgrid, gdims(5), nindex,   &
      &                                 nmiss, nvars, root_pe, var_ref_pos, lev
    CHARACTER(len=8)                 :: model_type
    REAL(wp), POINTER                :: r1d(:), rptr2d(:,:), rptr3d(:,:,:)
    CLASS(t_scatterPattern), POINTER :: scatter_pattern

    ! rank of broadcast root PE
    root_pe = 0

    IF (my_process_is_mpi_workroot()) THEN
      WRITE(0,*) "read_restart_files, nvar_lists=", nvar_lists
    END IF
    abbreviations(1:nvar_lists)%key          = 0
    abbreviations(1:nvar_lists)%abbreviation = ""
    key = 0
    n   = 1
    for_all_model_types: DO i = 1, nvar_lists
      ! skip var_list if it does not match the current patch ID
      IF (var_lists(i)%p%patch_id /= p_patch%id) CYCLE

      key = util_hashword(TRIM(var_lists(i)%p%model_type), LEN_TRIM(var_lists(i)%p%model_type), 0)
      IF (.NOT. ANY(abbreviations(1:n)%key == key)) THEN
        abbreviations(n)%abbreviation = var_lists(i)%p%model_type
        abbreviations(n)%key = key
        n = n+1
      ENDIF
    ENDDO for_all_model_types

    ALLOCATE(r1d(1),STAT=istat)
    IF (istat /= 0) THEN
      CALL finish('','allocation of r1d failed ...')
    ENDIF

    nfiles = n-1

    for_all_files: DO n = 1, nfiles
      model_type =TRIM(abbreviations(n)%abbreviation)
      IF (PRESENT(opt_ndom)) THEN
        IF (opt_ndom > 1) THEN
          restart_filename = 'restart_'//TRIM(model_type)//"_DOM"//TRIM(int2string(p_patch%id, "(i2.2)"))//'.nc'
        ELSE
          restart_filename = 'restart_'//TRIM(model_type)//'_DOM01.nc'
        END IF
      ELSE
        restart_filename = 'restart_'//TRIM(model_type)//'_DOM01.nc'
      END IF
      name = TRIM(restart_filename)//CHAR(0)

      IF (my_process_is_mpi_workroot()) THEN
        WRITE(0,*) "streamOpenRead ", TRIM(restart_filename)

        fileID  = streamOpenRead(TRIM(name))
        IF(fileID < 0) CALL finish(routine, "could not open file '"//TRIM(NAME)//"'")
        vlistID = streamInqVlist(fileID)
        taxisID = vlistInqTaxis(vlistID)

        idate   = taxisInqVdate(taxisID)
        itime   = taxisInqVtime(taxisID)
      END IF
      CALL message('read_restart_files', 'Read restart for: '//TRIM(int2string(idate))//'T'//TRIM(int2string(itime))//'Z '//&
                                       & 'from '//TRIM(restart_filename))

      IF (my_process_is_mpi_workroot()) THEN
        nvars = vlistNvars(vlistID)
      END IF
      CALL p_bcast(nvars, root_pe, comm=p_comm_work)

      for_all_vars: DO varID = 0, (nvars-1)

        IF (my_process_is_mpi_workroot()) THEN
          CALL vlistInqVarName(vlistID, varID, name)
        END IF
        CALL p_bcast(name, root_pe, comm=p_comm_work)

        for_all_lists: DO i = 1, nvar_lists
          ! skip var_list if it does not match the current patch ID
          IF (var_lists(i)%p%patch_id   /= p_patch%id) CYCLE
          IF (var_lists(i)%p%model_type /= model_type) CYCLE

          element => find_list_element(var_lists(i), TRIM(name))
          IF (ASSOCIATED(element)) THEN
            IF (element%field%info%lrestart) THEN

              info => element%field%info

              ! allocate temporary global array on output processor
              ! and gather field from other processors

              NULLIFY(rptr2d)
              NULLIFY(rptr3d)

              IF (my_process_is_mpi_workroot()) THEN
                gridID  = vlistInqVarGrid(vlistID, varID)
                zaxisID = vlistInqVarZaxis(vlistID, varID)
                vgrid   = zaxisInqType(zaxisID)
                ic = gridInqSize(gridID)
                IF (SIZE(r1d, 1) /= ic) THEN
                  DEALLOCATE(r1d)
                  ALLOCATE(r1d(ic),STAT=istat)
                  IF (istat /= 0) THEN
                    CALL finish('','allocation of r1d failed ...')
                  ENDIF
                END IF
                IF (vgrid == ZAXIS_SURFACE) THEN
                  il = 1
                ELSE
                  il = zaxisInqSize(zaxisID)
                ENDIF
              END IF

              CALL p_bcast(il, root_pe, comm=p_comm_work)

              IF (info%lcontained) THEN
                nindex = info%ncontained
              ELSE
                nindex = 1
              ENDIF


              SELECT CASE(info%ndims)
              CASE (2)
                var_ref_pos = 3
                IF (info%lcontained)  var_ref_pos = info%var_ref_pos
                SELECT CASE(var_ref_pos)
                CASE (1)
                  rptr2d => element%field%r_ptr(nindex,:,:,1,1)
                CASE (2)
                  rptr2d => element%field%r_ptr(:,nindex,:,1,1)
                CASE (3)
                  rptr2d => element%field%r_ptr(:,:,nindex,1,1)
                CASE default
                  CALL finish(routine, "internal error!")
                END SELECT
              CASE (3)
                var_ref_pos = 4
                IF (info%lcontained)  var_ref_pos = info%var_ref_pos
                SELECT CASE(var_ref_pos)
                CASE (1)
                  rptr3d => element%field%r_ptr(nindex,:,:,:,1)
                CASE (2)
                  rptr3d => element%field%r_ptr(:,nindex,:,:,1)
                CASE (3)
                  rptr3d => element%field%r_ptr(:,:,nindex,:,1)
                CASE (4)
                  rptr3d => element%field%r_ptr(:,:,:,nindex,1)
                CASE default
                  CALL finish(routine, "internal error!")
                END SELECT
              CASE DEFAULT
                CALL finish(routine, "internal error!")
              END SELECT

              SELECT CASE (info%hgrid)
              CASE (GRID_UNSTRUCTURED_CELL)
                scatter_pattern => p_patch%comm_pat_scatter_c
              CASE (GRID_UNSTRUCTURED_VERT)
                scatter_pattern => p_patch%comm_pat_scatter_v
              CASE (GRID_UNSTRUCTURED_EDGE)
                scatter_pattern => p_patch%comm_pat_scatter_e
              CASE default
                CALL finish('out_stream','unknown grid type')
              END SELECT

              DO lev = 1, il
                IF (my_process_is_mpi_workroot()) THEN
                  CALL streamReadVarSlice(fileID, varID, lev-1, r1d, nmiss)
                ENDIF
                IF (info%ndims == 2) THEN
                  CALL scatter_pattern%distribute(r1d, rptr2d, .FALSE.)
                ELSE
                  CALL scatter_pattern%distribute(r1d, rptr3d(:,lev,:), .FALSE.)
                ENDIF
              END DO


              ! deallocate temporary global arrays


              IF (my_process_is_mpi_workroot()) THEN
                WRITE (0,*) ' ... read ',TRIM(element%field%info%name)
              ENDIF
              CYCLE for_all_vars
            ENDIF
          ENDIF
        ENDDO for_all_lists
        CALL message('reading_restart_file','Variable '//TRIM(name)//' not defined.')
      ENDDO for_all_vars

      IF (my_process_is_mpi_workroot())  CALL streamClose(fileID)

    ENDDO for_all_files

    DEALLOCATE (r1d)

    CALL message('','')
    CALL message('',separator)
    CALL message('','')

  END SUBROUTINE read_restart_files
END MODULE mo_load_restart
