!> Module for reading multifile restart files
!!
!! 2018-08: Major revision / revamp / refactoring : Harald Braun (Atos SE)
!!
!! Note: The single file implementation of the restart input can be
!!       found in the module "mo_load_singlefile_restart"
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! ----------------------------------------------------------------------------------------------------
!!
!! All work processes are used to read the data from the
!! multifile. Each processor handles all domains of its horizontal
!! "chunk", since the domain loop is located outside of this routine,
!! in the calling "src/drivers" routine.
!!
#include "omp_definitions.inc"

MODULE mo_load_multifile_restart

  USE ISO_C_BINDING,             ONLY: C_INT
  USE mo_c_restart_util,         ONLY: checkMultifileDir
  USE mo_impl_constants,         ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_cdi_constants,          ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_VERT
  USE mo_cdi,                    ONLY: streamOpenRead, streamInqVlist, vlistNvars, vlistinqvarname, &
    &                                  streamClose, streamReadVar, streamReadVarF, cdiInqAttInt, &
    &                                  CDI_GLOBAL, streamReadVarSlice, streamReadVarSliceF, &
    &                                  cdi_max_name, CDI_UNDEFID
  USE mo_communication,          ONLY: t_comm_pattern, t_p_comm_pattern, exchange_data
  USE mo_communication_factory,  ONLY: setup_comm_pattern
  USE mo_decomposition_tools,    ONLY: t_glb2loc_index_lookup, init_glb2loc_index_lookup, set_inner_glb_index, &
    &                                  deallocate_glb2loc_index_lookup
  USE mo_dynamics_config,        ONLY: nnew, nnew_rcf
  USE mo_exception,              ONLY: finish, warning
  USE mo_kind,                   ONLY: sp, dp
  USE mo_model_domain,           ONLY: t_patch
  USE mo_mpi,                    ONLY: p_comm_work, p_comm_size, p_comm_rank, my_process_is_work, &
    &                                  p_allreduce, mpi_sum, my_process_is_mpi_workroot, p_alltoall, p_alltoallv
  USE mo_multifile_restart_util, ONLY: multifilePayloadPath, rBuddy, rGroup, vNames_glbIdx
  USE mo_parallel_config,        ONLY: nproma, idx_no, blk_no, restart_load_scale_max
  USE mo_restart_nml_and_att,    ONLY: getAttributesForRestarting, ocean_initFromRestart_OVERRIDE
  USE mo_key_value_store,        ONLY: t_key_value_store
  USE mo_restart_var_data,       ONLY: get_var_3d_ptr, has_valid_time_level
  USE mo_var,                    ONLY: t_var_ptr
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_load_restart_io, timers_level, &
    &                                  timer_load_restart_comm_setup, timer_load_restart_communication, &
    &                                  timer_load_restart_get_var_id

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: multifileCheckRestartFiles, multifileReadPatch

  TYPE t_PayloadFile
    INTEGER :: streamId, vlistId, varCount
    INTEGER :: iCnts(3), iVarIds(3)
  END TYPE t_PayloadFile

  CHARACTER(*), PARAMETER :: modname = "mo_load_multifile_restart"

CONTAINS

  SUBROUTINE multifileCheckRestartFiles(filename)
    CHARACTER(*), INTENT(IN) :: filename
    TYPE(t_key_value_store), POINTER :: restartAttributes
    INTEGER :: n_dom, multifile_file_count
    CHARACTER(*), PARAMETER :: routine = modname//":multifileCheckRestartFiles"

    CALL getAttributesForRestarting(restartAttributes)
    CALL restartAttributes%get('multifile_n_dom_active', n_dom)
    CALL restartAttributes%get('multifile_file_count', multifile_file_count)
    IF(checkMultifileDir(filename, n_dom, multifile_file_count) /= SUCCESS) &
      CALL finish(routine, "'"//filename//"' is not a valid restart multifile")
  END SUBROUTINE multifileCheckRestartFiles

  ! Opens all the payload files that should be opened by this
  ! process, creating t_PayloadFile objects for them, AND
  ! calculating the total numbers of cells/edges/vertices contained
  ! IN these files.
  !
  ! This FUNCTION performs the distribution of the payload files
  ! among the restart processes, subsequent code just has to iterate
  ! over the array returned by this FUNCTION.
  SUBROUTINE openPayloadFiles(mfPath, mfCnt, dom, tCnt, files, at_once)
    CHARACTER(*), INTENT(IN) :: mfPath
    INTEGER, INTENT(IN) :: mfCnt, dom
    INTEGER, INTENT(OUT) :: tCnt(3)
    TYPE(t_PayloadFile), INTENT(INOUT), ALLOCATABLE :: files(:)
    LOGICAL, INTENT(OUT) :: at_once
    INTEGER :: ierr, pCnt, myR, myFCnt, myFF, i
    CHARACTER(*), PARAMETER :: routine = modname//":openPayloadFiles"

    IF(.NOT.my_process_is_work()) CALL finish(routine, "this is not a work process")
    pCnt = p_comm_size(p_comm_work)
    myR = p_comm_rank(p_comm_work)
! number of reader ranks buffering one full 3d field of at max restart_load_scale_max worker ranks
! default : restart_load_scale_max = 1
! if restart_load_scale_max == 0 - - - - - - - -> read data in a LEVEL-based loop
! if restart_load_scale_max <  #ranks / #files -> read data in a LEVEL-based loop
! if restart_load_scale_max >= #ranks / #files -> read data in a VARIABLE-based loop
! if restart_load_scale_max < 0  - - - - - - - -> read data in a VARIABLE-based loop
    at_once = pCnt .LE. mfCnt * &
      & MERGE(restart_load_scale_max, pCnt, restart_load_scale_max .GE. 0)
    myFF = -1
    IF (pCnt .GE. mfCnt) THEN
      myFCnt = MERGE(1, 0, rBuddy(nr_in=mfCnt) .EQ. myR)
      IF (myFCnt .GT. 0) myFF = rGroup(nr_in=mfCnt)
    ELSE
      myFCnt = COUNT( (/ (rGroup(nr_in=pCnt, nw_in=mfCnt, pe_in=i) .EQ. myR, &
        &                    i = myR, mfCnt-1) /) )
      i = myR
      DO WHILE(myFF .LT. 0)
        IF (rGroup(nr_in=pCnt, nw_in=mfCnt, pe_in=i) .EQ. myR) myFF = i
        i = i + 1
      END DO
    END IF
#ifdef DEBUG
    WRITE(0, "(3(a,i4),a)") "restart: rank ", myR, " will start reading from file ", &
      &                      myFF, ", ", myFCnt, " files in total"
#endif
    ALLOCATE(files(myFCnt), STAT = ierr)
    IF(ierr /= SUCCESS) CALL finish(routine, "memory allocation failure")
    tCnt(:) = 0
    DO i = 1, myFCnt
      CALL payloadfile_open(files(i), myFF - 1 + i)
      tCnt(:) = tCnt(:) + files(i)%iCnts(:)
    END DO
  CONTAINS

    SUBROUTINE payloadFile_open(me, partId)
      TYPE(t_PayloadFile), INTENT(INOUT) :: me
      INTEGER, INTENT(in) :: partId
      INTEGER :: varId, trash, icount(1), iG
      CHARACTER(*), PARAMETER :: routine = modname//":payloadFile_init", &
       & attNames(3) = (/"cellCount", "vertCount", "edgeCount"/)
      CHARACTER(:), ALLOCATABLE :: pathname
      CHARACTER(len=cdi_max_name) :: varname
  
      IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
      CALL multifilePayloadPath(mfPath, dom, partId, pathname)
      me%streamId = streamOpenRead(pathname)
      IF (me%streamId < 0) CALL finish( routine, "failed to open file "//pathname)
      me%vlistId = streamInqVlist(me%streamId)
      me%varCount = vlistNvars(me%vlistId)
      me%iVarIds(:) = CDI_UNDEFID
      varId = 0
      DO WHILE(ANY(me%iVarIds(:) .EQ. CDI_UNDEFID) .AND. &
        &      me%varCount .GT. varId)
        varname = ''
        CALL vlistInqVarName(me%vlistId, varId, varname)
        DO iG = 1, 3
          IF (varname == vNames_glbIdx(iG)) THEN
            IF(me%iVarIds(iG) .NE. CDI_UNDEFID) &
              & CALL finish(routine, "corrupted: doubly defined "//vNames_glbIdx(iG))
            trash = cdiInqAttInt(me%vlistId, CDI_GLOBAL, attNames(iG), 1, icount);
            me%iCnts(iG) = icount(1)
            me%iVarIds(iG) = varId
          END IF
        END DO
        varId = varId + 1
      END DO
      IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
      IF (ANY(me%iVarIds(:) .EQ. CDI_UNDEFID)) &
        & CALL finish(routine, "corrupted: index variable missing")
    END SUBROUTINE payloadFile_open
  END SUBROUTINE openPayloadFiles

  ! This asserts that the given pattern redistributes a given array
  ! exactly as described by the global indices arrays.
  SUBROUTINE checkRedistributionPattern(pat, provGlbIdces, reqdGlbIdces)
    CLASS(t_comm_pattern), INTENT(IN), POINTER :: pat
    INTEGER, INTENT(IN) :: provGlbIdces(:), reqdGlbIdces(:)
    INTEGER :: iSize, oSize, i, error 
    CHARACTER(*), PARAMETER :: routine = modname//":checkRedistributionPattern"
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: input, output

    iSize = SIZE(provGlbIdces)
    oSize = SIZE(reqdGlbIdces)
    ALLOCATE(input(nproma, blk_no(iSize)), output(nproma, blk_no(oSize)))
!ICON_OMP PARALLEL PRIVATE(i)
!ICON_OMP DO SCHEDULE(STATIC)
    DO i = 1, blk_no(oSize)
      output(:,i) = REAL(-1, dp)
    END DO
!ICON_OMP DO SCHEDULE(STATIC)
    DO i = 1, iSize
      input(idx_no(i), blk_no(i)) = REAL(provGlbIdces(i), dp)
    END DO
!ICON_OMP END PARALLEL
    CALL exchange_data(pat, output, input)
    error = 0
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO i = 1, oSize
      IF(INT(output(idx_no(i), blk_no(i))) /= reqdGlbIdces(i)) THEN
!ICON_OMP CRITICAL
        error = error + 1
        PRINT "(6(i6,1x))", i, blk_no(i), idx_no(i), oSize, &
          & INT(output(idx_no(i), blk_no(i))), reqdGlbIdces(i)
!ICON_OMP END CRITICAL
      END IF
    END DO
    IF (error .GT. 0) THEN
      PRINT "(2(a,i8))", "check comm pattern: found deviations: ", &
        & error, " out of ", oSize
      CALL finish(routine, "comm pattern does not work as expected")
    END IF
    DEALLOCATE(input, output)
  END SUBROUTINE checkRedistributionPattern

  ! Create a t_comm_pattern that IS able to redistribute the DATA as
  ! it IS READ from the multifile to the domain decomposition of the
  ! work processes.  The tricky part about this IS, that
  ! setup_comm_pattern() requires an owner(:) array as input, which
  ! we have to create first.
  !
  ! The resulting t_comm_pattern will also initialize the halo points.
  FUNCTION makeRedistributionPattern(pIds, rIds) RESULT(pat)
    CLASS(t_comm_pattern), POINTER :: pat
    INTEGER, INTENT(IN) :: pIds(:), rIds(:)
    INTEGER :: glbSize, locSize, cuml, np, rnk, ii
    INTEGER, ALLOCATABLE :: own_src(:), own_dst(:), bd(:)
    TYPE(t_glb2loc_index_lookup) :: lupTbl
    CHARACTER(*), PARAMETER :: routine = modname//":makeRedistributionPattern"

    IF(timers_level >= 7) CALL timer_start(timer_load_restart_comm_setup)
    locSize = SIZE(pIds)
    np = p_comm_size(p_comm_work)
    rnk = p_comm_rank(p_comm_work)
    glbSize = p_allreduce(locSize, mpi_sum, p_comm_work)
    CALL get_owners()
    CALL init_glb2loc_index_lookup(lupTbl, glbSize)
    CALL set_inner_glb_index(lupTbl, pIds, [(ii, ii = 1, locSize)])
    CALL setup_comm_pattern(SIZE(rIds), own_dst, rIds, lupTbl, locSize, &
      & own_src, pIds, pat, inplace=.FALSE., comm=p_comm_work)
    CALL deallocate_glb2loc_index_lookup(lupTbl)
    DEALLOCATE(own_dst, own_src)
#ifdef DEBUG
    CALL checkRedistributionPattern(pat, pGlbIds, rGlbIds)
#endif
    IF(timers_level >= 7) CALL timer_stop(timer_load_restart_comm_setup)
  CONTAINS

    SUBROUTINE get_owners()
      INTEGER, ALLOCATABLE, DIMENSION(:) :: dI_s, dI_c, dI_i, dI_o, o_m, &
        & iPk_s, iPk_c, iPk_i, iPk_o, nI_s, nI_c, nI_i, nI_o, oPk_i, oPk_o, oPk_c
      INTEGER :: i, j

#ifdef DEBUG 
      ! consistency checks
      IF(1 > MINVAL(pIds)) &
        & CALL finish(routine, "global index out of RANGE (<1)")
      IF(glbSize < MAXVAL(pIds)) &
        &  CALL finish(routine, "global index out of range (>globalSize)")
#endif
      ALLOCATE(bd(0:np), nI_i(np), nI_o(np))
      bd(:) = (/(INT(REAL(i,dp)*(REAL(glbSize,dp)/REAL(np,dp))), i=0,np)/)
      bd(np) = glbSize
      CALL get_disps(dI_s, nI_s, iPk_s, pIds)
      CALL get_disps(dI_c, nI_c, iPk_c, rIds)
      CALL p_alltoall(nI_s, nI_i, p_comm_work)
      CALL p_alltoall(nI_c, nI_o, p_comm_work)
      ALLOCATE(dI_o(np+1), dI_i(np+1))
      dI_o(1) = 0
      dI_i(1) = 0
      DO i = 2, np+1
        dI_o(i) = dI_o(i-1) + nI_o(i-1)
        dI_i(i) = dI_i(i-1) + nI_i(i-1)
      END DO
      ALLOCATE(oPk_i(dI_i(np+1)), o_m(bd(rnk)+1:bd(rnk+1)), &
        & oPk_o(dI_o(np+1)), iPk_o(dI_o(np+1)), oPk_c(dI_c(np+1)), &
        & iPk_i(dI_i(np+1)))
      CALL p_alltoallv(iPk_s, nI_s, dI_s, iPk_i, nI_i, dI_i, p_comm_work)
      CALL p_alltoallv(iPk_c, nI_c, dI_c, iPk_o, nI_o, dI_o, p_comm_work)
      DEALLOCATE(dI_s, nI_s, nI_i, iPk_s, iPk_c)
!ICON_OMP PARALLEL PRIVATE(i)
!ICON_OMP DO SCHEDULE(GUIDED)
      DO i = 1, np
        IF (dI_i(i)+1 .LE. dI_i(i+1)) &
          & oPk_i(dI_i(i)+1:dI_i(i+1)) = i - 1
      END DO
!ICON_OMP DO SCHEDULE(GUIDED)
      DO i = 1, dI_i(np+1)
        o_m(iPk_i(i)) = oPk_i(i)
      END DO
!ICON_OMP DO SCHEDULE(GUIDED)
      DO i = 1, dI_o(np+1)
        oPk_o(i) = o_m(iPk_o(i))
      END DO
!ICON_OMP END PARALLEL
      CALL p_alltoallv(oPk_o, nI_o, dI_o, oPk_c, nI_c, dI_c, p_comm_work)
      DEALLOCATE(o_m, oPk_o, iPk_o, iPk_i, dI_i, nI_o, dI_o)
      ALLOCATE(own_dst(SIZE(rIds)), own_src(locSize))
!ICON_OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(cuml,j)
      DO i = 1, np
        IF (dI_c(i+1)-dI_c(i) .GT. 0) THEN
          cuml = dI_c(i)
          DO j = 1, SIZE(rIds)
            IF (rIds(j) .GT. bd(i-1) .AND. rIds(j) .LE. bd(i)) THEN
              cuml = cuml + 1
              own_dst(j) = oPk_c(cuml)
            END IF
          END DO
        END IF
      END DO
      DEALLOCATE(bd, oPk_c, nI_c, dI_c)
      own_src(:) = rnk
    END SUBROUTINE get_owners

    SUBROUTINE get_disps(dId, nId, idPk, id)
      INTEGER, INTENT(INOUT), ALLOCATABLE :: dId(:), nId(:), idPk(:)
      INTEGER, INTENT(IN), CONTIGUOUS :: id(:)
      INTEGER :: i, j

      ALLOCATE(dId(np+1), idPk(SIZE(id)), nId(np))
      cuml = 0
      DO i = 1, np
        dId(i) = cuml
        DO j = 1, SIZE(id)
          IF (id(j) .GT. bd(i-1) .AND. id(j) .LE. bd(i)) THEN
            cuml = cuml + 1
            idPk(cuml) = id(j)
          END IF
        END DO
        nId(i) = cuml - dId(i)
      END DO
      dId(np+1) = cuml
    END SUBROUTINE get_disps
  END FUNCTION makeRedistributionPattern

  SUBROUTINE multifileReadPatch(vDat, p_patch, multifilePath)
    TYPE(t_var_ptr), INTENT(IN) :: vDat(:)
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    CHARACTER(*), INTENT(IN) :: multifilePath
    TYPE(t_p_comm_pattern) :: cpat(3)
    LOGICAL :: load_var_at_once
    INTEGER :: hi
    TYPE(t_PayloadFile), ALLOCATABLE :: files(:)

    CALL construct()
    CALL readData()
    DO hi = 1, 3
      CALL cpat(hi)%p%delete()
    END DO
    IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
    IF(ALLOCATED(files)) THEN
      DO hi = 1, SIZE(files)
        CALL streamClose(files(hi)%streamId)
      END DO
    END IF
    IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
  CONTAINS

  SUBROUTINE construct()
    INTEGER :: mfileCnt, ierr, cFId, tCnt(3), cOff(3), iG, i, n
    INTEGER(C_INT) :: trash
    INTEGER, ALLOCATABLE :: glbidx_read(:)
    INTEGER, POINTER :: glb_index(:)
    REAL(dp), ALLOCATABLE :: buffer(:)
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchReader_construct"

    CALL getAttributesForRestarting(restartAttributes)
    CALL restartAttributes%get('multifile_file_count', mfileCnt)
    IF (my_process_is_mpi_workroot()) &
      & WRITE(0, *) "reading from ", mfileCnt, " files/patch."
    CALL openPayloadFiles(multifilePath, mfileCnt, p_patch%id, tCnt, files, load_var_at_once)
    cOff(:) = 0
    ALLOCATE(glbidx_read(MAXVAL(tCnt(:),1)), STAT = ierr)
    IF (ierr /= SUCCESS) CALL finish(routine, "memory allocation failure")
    DO iG = 1, 3
      DO cFId = 1, SIZE(files)
        n = files(cFId)%iCnts(iG)
        IF (n .LE. 0) CYCLE
        ALLOCATE(buffer(n), STAT = ierr)
        IF (ierr /= SUCCESS) CALL finish(routine, "memory allocation failed")
        IF (timers_level >= 7) CALL timer_start(timer_load_restart_io)
        CALL streamReadVar(files(cFId)%streamId, files(cFId)%iVarIds(iG), &
          &                buffer, trash)
        IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        DO i = 1, n
          glbidx_read(cOff(iG) + i) = INT(buffer(i))
        END DO
        DEALLOCATE(buffer)
        cOff(iG) = cOff(iG) + n
      END DO
      SELECT CASE (ig)
      CASE (GRID_UNSTRUCTURED_CELL)
        glb_index => p_patch%cells%decomp_info%glb_index
      CASE (GRID_UNSTRUCTURED_VERT)
        glb_index => p_patch%verts%decomp_info%glb_index
      CASE (GRID_UNSTRUCTURED_EDGE)
        glb_index => p_patch%edges%decomp_info%glb_index
      END SELECT
      cpat(ig)%p => makeRedistributionPattern(glbidx_read(:cOff(iG)), glb_index)
    END DO
  END SUBROUTINE construct

  SUBROUTINE readData()
    INTEGER :: vId, lId, lcnt, fId, varID, pCt(SIZE(files)), vlID
    INTEGER :: rbuf_size, ofs, nxt, dummy, i, hgrid, mock_nblk, dom
    REAL(KIND=dp), POINTER :: ptr_3d_d(:,:,:), buf_3d_d(:,:,:), buf_d(:)
    REAL(KIND=sp), POINTER :: ptr_3d_s(:,:,:), buf_3d_s(:,:,:), buf_s(:)
    INTEGER,       POINTER :: ptr_3d_i(:,:,:), buf_3d_i(:,:,:), buf_i(:)
    TYPE(t_key_value_store) :: vname_map
    CHARACTER(LEN=cdi_max_name) :: cVname
    CHARACTER(*), PARAMETER :: routine = modname//":multifilePatchReader_readData"
    
    CALL vname_map%init(.FALSE.)
    IF (SIZE(files) .GT. 0) THEN
      vlID = streamInqVlist(files(1)%streamId)
      DO vId = 0, vlistNvars(vlID) - 1
        CALL vlistInqVarName(vlID, vId, cVname)
        CALL vname_map%put(cVname, vId)
      END DO
    END IF
    dom = p_patch%id
    IF (ocean_initFromRestart_OVERRIDE) CALL vname_map%bcast(0, p_comm_work)
    DO vId = 1, SIZE(vDat)
      IF (.NOT.has_valid_time_level(vDat(vId)%p%info, dom, nnew(dom), nnew_rcf(dom))) CYCLE
      IF (timers_level >= 7) CALL timer_start(timer_load_restart_get_var_id)
      CALL vname_map%get(vDat(vId)%p%info%NAME, varId, opt_err=dummy)
      IF (timers_level >= 7) CALL timer_stop(timer_load_restart_get_var_id)
      IF (dummy .NE. 0) THEN
        IF (ocean_initFromRestart_OVERRIDE) THEN
! fatal hack from coding hell to make init_fromRestart=.true. (ocean) work
          CALL warning(routine, "variable not found: '" // TRIM(vDat(vId)%p%info%NAME))
          CALL warning(routine, &
            & "that MAY be intended if initialize_fromRestart=.true.")
          CYCLE
        ELSE IF (SIZE(files) .GT. 0) THEN
          CALL finish(routine, "variable not found: "//TRIM(vDat(vId)%p%info%NAME))
        END IF
      END IF
      hgrid = vDat(vId)%p%info%hgrid
      IF (hgrid < 1 .OR. hgrid > 3) &
        CALL finish(routine, "unexpected varData(varIndex)%info%hgrid")
      pCt = files(:)%iCnts(hgrid)
      rbuf_size = SUM(pCt)
      mock_nblk = blk_no(rbuf_size)
      rbuf_size = mock_nblk * nproma
      IF (load_var_at_once) THEN
        SELECT CASE(vDat(vId)%p%info%data_type)
        CASE(REAL_T)
          CALL get_var_3d_ptr(vDat(vId)%p, ptr_3d_d)
          lcnt = SIZE(ptr_3d_d, 2)
          ALLOCATE(buf_d(MAXVAL(pCt)*lCnt), buf_3d_d(nproma,lcnt,mock_nblk))
          ofs = 0
          IF (timers_level >= 7) CALL timer_start(timer_load_restart_io)
          DO fId = 1, SIZE(files)
            IF (pCt(fId) .LT. 1) CYCLE
            CALL streamReadVar(files(fId)%streamId, varId, buf_d, dummy)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2)
            DO lId = 1, lCnt
              DO i = 1, pCt(fId)
                buf_3d_d(idx_no(i+ofs), lId, blk_no(i+ofs)) = &
                  & buf_d(i+(lId-1)*pCt(fId))
              END DO
            END DO
            ofs = ofs + pCt(fId)
          END DO
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
          IF (timers_level >= 7) CALL timer_start(timer_load_restart_communication)
          CALL exchange_data(cpat(hgrid)%p, ptr_3d_d, buf_3d_d)
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
          DEALLOCATE(buf_3d_d, buf_d)
        CASE(SINGLE_T)
          CALL get_var_3d_ptr(vDat(vId)%p, ptr_3d_s)
          lcnt = SIZE(ptr_3d_s, 2)
          ALLOCATE(buf_s(MAXVAL(pCt)*lCnt), buf_3d_s(nproma,lcnt,mock_nblk))
          ofs = 0
          IF (timers_level >= 7) CALL timer_start(timer_load_restart_io)
          DO fId = 1, SIZE(files)
            IF (pCt(fId) .LT. 1) CYCLE
            CALL streamReadVarF(files(fId)%streamId, varId, buf_s, dummy)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2)
            DO lId = 1, lCnt
              DO i = 1, pCt(fId)
                buf_3d_s(idx_no(i+ofs), lId, blk_no(i+ofs)) = &
                  & buf_s(i+(lId-1)*pCt(fId))
              END DO
            END DO
            ofs = ofs + pCt(fId) 
          END DO
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
          IF (timers_level >= 7) CALL timer_start(timer_load_restart_communication)
          CALL exchange_data(cpat(hgrid)%p, ptr_3d_s, buf_3d_s)
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
          DEALLOCATE(buf_3d_s, buf_s)
        CASE(INT_T)
          CALL get_var_3d_ptr(vDat(vId)%p, ptr_3d_i)
          lcnt = SIZE(ptr_3d_i, 2)
          ALLOCATE(buf_d(MAXVAL(pCt)*lCnt), buf_3d_i(nproma,lcnt,mock_nblk))
          ofs = 0
          IF (timers_level >= 7) CALL timer_start(timer_load_restart_io)
          DO fId = 1, SIZE(files)
            IF (pCt(fId) .LT. 1) CYCLE
            CALL streamReadVar(files(fId)%streamId, varId, buf_d, dummy)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2)
            DO lId = 1, lCnt
              DO i = 1, pCt(fId)
                buf_3d_i(idx_no(i+ofs), lId, blk_no(i+ofs)) = &
                  & INT(buf_d(i+(lId-1)*pCt(fId)))
              END DO
            END DO
            ofs = ofs + pCt(fId)
          END DO
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
          IF (timers_level >= 7) CALL timer_start(timer_load_restart_communication)
          CALL exchange_data(cpat(hgrid)%p, ptr_3d_i, buf_3d_i)
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
          DEALLOCATE(buf_3d_i, buf_d)
        END SELECT
      ELSE
        SELECT CASE(vDat(vId)%p%info%data_type)
        CASE(REAL_T)
          CALL get_var_3d_ptr(vDat(vId)%p, ptr_3d_d)
          lcnt = SIZE(ptr_3d_d, 2)
          ALLOCATE(buf_d(rbuf_size))
          DO lId = 1, lcnt
            ofs = 0
            IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
            DO fId = 1, SIZE(files)
              IF (pCt(fId) .LT. 1) CYCLE
              nxt = ofs + pCt(fId)
              CALL streamReadVarSlice(files(fId)%streamId, varId, lid-1, &
                   buf_d(ofs+1:nxt), dummy)
              ofs = nxt
            END DO
            IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
            IF (timers_level >= 7) CALL timer_start(timer_load_restart_communication)
            buf_3d_d(1:nproma, 1:1, 1:mock_nblk) => buf_d
           CALL exchange_data(cpat(hgrid)%p, ptr_3d_d(:,lid:lid,:), buf_3d_d)
            IF (timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
          END DO
          DEALLOCATE(buf_d)
        CASE(SINGLE_T)
          ALLOCATE(buf_s(rbuf_size))
          CALL get_var_3d_ptr(vDat(vId)%p, ptr_3d_s)
          lcnt = SIZE(ptr_3d_s, 2)
          DO lId = 1, lcnt
            ofs = 0
            IF (timers_level >= 7) CALL timer_start(timer_load_restart_io)
            DO fId = 1, SIZE(files)
              IF (pCt(fId) .LT. 1) CYCLE
              nxt = ofs + pCt(fId)
              CALL streamReadVarSliceF(files(fId)%streamId, varId, lid-1, &
                   buf_s(ofs+1:nxt), dummy)
              ofs = nxt
            END DO
            IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
            IF (timers_level >= 7) CALL timer_start(timer_load_restart_communication)
            buf_3d_s(1:nproma, 1:1, 1:mock_nblk) => buf_s
            CALL exchange_data(cpat(hgrid)%p, ptr_3d_s(:,lid:lid,:), buf_3d_s)
            IF (timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
          END DO
          DEALLOCATE(buf_s)
        CASE(INT_T)
          ALLOCATE(buf_i(rbuf_size), buf_d(MAX(1,MAXVAL(pCt))))
          CALL get_var_3d_ptr(vDat(vId)%p, ptr_3d_i)
          lcnt = SIZE(ptr_3d_i, 2)
          DO lId = 1, lcnt
            ofs = 0
            IF (timers_level >= 7) CALL timer_start(timer_load_restart_io)
            DO fId = 1, SIZE(files)
              IF (pCt(fId) .LT. 1) CYCLE
              CALL streamReadVarSlice(files(fId)%streamId, varId, lid-1, &
                   buf_d(:pCt(fId)), dummy)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
              DO i = 1, pCt(fId)
                buf_i(ofs+i) = INT(buf_d(i))
              END DO
              ofs = ofs + pCt(fId)
            END DO
            IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
            IF (timers_level >= 7) CALL timer_start(timer_load_restart_communication)
            buf_3d_i(1:nproma, 1:1, 1:mock_nblk) => buf_i
            CALL exchange_data(cpat(hgrid)%p, ptr_3d_i(:,lid:lid,:), buf_3d_i)
            IF (timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
          END DO
          DEALLOCATE(buf_d, buf_i)
        END SELECT
      END IF
    END DO
  END SUBROUTINE readData

  END SUBROUTINE multifileReadPatch

END MODULE mo_load_multifile_restart
