!>
!! Contains utility routines for reading NetCDF and GRIB2 files (using
!! the CDI library) and communicating fields in parallel.
!!
!! @author F. Prill, DWD
!!
!!
!! @par Revision History
!! Initial revision: 2013-02-19 : F. Prill, DWD
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_util_cdi

  USE mo_kind,               ONLY: wp, sp, dp, i8
  USE mo_exception,          ONLY: finish
  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_scatter_pattern,    ONLY: t_scatterPattern
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_run_config,         ONLY: msg_level
  USE mo_mpi,                ONLY: p_bcast, p_comm_work, p_comm_work_test,  &
    &                              p_io, my_process_is_stdio, p_mpi_wtime
  USE mo_util_string,        ONLY: tolower
  USE mo_fortran_tools,      ONLY: assign_if_present
  USE mo_dictionary,         ONLY: t_dictionary, dict_get, dict_init, dict_copy, dict_finalize, DICT_MAX_STRLEN
  USE mo_gribout_config,     ONLY: t_gribout_config
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_action,             ONLY: ACTION_RESET
  ! calendar operations
  USE mtime,                 ONLY: timedelta, newTimedelta,                 &
    &                              datetime, newDatetime,                   &
    &                              deallocateTimedelta, deallocateDatetime, &
    &                              MAX_DATETIME_STR_LEN
  USE mo_mtime_extensions,   ONLY: getTimeDeltaFromDateTime

  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PRIVATE

  PUBLIC :: read_cdi_2d, read_cdi_3d
  PUBLIC :: get_cdi_varID
  PUBLIC :: test_cdi_varID
  PUBLIC :: set_additional_GRIB2_keys
  PUBLIC :: set_timedependent_GRIB2_keys
  PUBLIC :: t_inputParameters, makeInputParameters, deleteInputParameters

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_cdi'

  INTERFACE read_cdi_2d
    MODULE PROCEDURE read_cdi_2d_real
    MODULE PROCEDURE read_cdi_2d_int
    MODULE PROCEDURE read_cdi_2d_time
  END INTERFACE

  INTERFACE read_cdi_3d
    MODULE PROCEDURE read_cdi_3d_real
  END INTERFACE

  ! This is a small type that serves two functions:
  ! 1. It encapsulates three parameters to the read functions into one, significantly reducing the hassle to call them.
  ! 2. It allows some information about the file to be distributed to the non-I/O PEs in a single broadcast that would otherwise
  !    require a broadcast in every call.
  !
  ! Create this by calling makeInputParameters(), destroy with deleteInputParameters().
  TYPE :: t_inputParameters
    !PRIVATE! Don't use outside of this file!
    INTEGER :: streamId     !< CDI stream ID
    INTEGER :: glb_arr_len  !< global array length
    TYPE (t_dictionary) :: dict     !< a dictionary that is used to translate variable names
    CLASS(t_scatterPattern), POINTER :: distribution    !< a t_scatterPattern to distribute the data
    LOGICAL :: have_dict    !< whether `dict` is used or not
    CHARACTER(LEN=MAX_CHAR_LENGTH), ALLOCATABLE :: variableNames(:) !< the names of the variables
    INTEGER, ALLOCATABLE :: variableDatatype(:) !< the datatypes that are used by the file to store the variables
    INTEGER, ALLOCATABLE :: variableTileIdx(:)  !< tile indices of the variables
    REAL(dp) :: readDuration    !< statistic on how long we took to read the data
    INTEGER(i8) :: readBytes    !< statistic on how much data we read from the file

  CONTAINS
    PROCEDURE :: findVarId => inputParametersFindVarId  !< determine the ID of an named variable
    PROCEDURE :: findVarDatatype => inputParametersFindVarDatatype  !< determine the type with which this variable is stored on disk
    PROCEDURE :: lookupDatatype => inputParametersLookupDatatype    !< determine datatype based on the variable ID

    GENERIC :: getDatatype => findVarDatatype, lookupDatatype
  END TYPE

CONTAINS


  !---------------------------------------------------------------------------------------------------------------------------------
  !> Creates a t_inputParameters object for the given arguments.
  !---------------------------------------------------------------------------------------------------------------------------------
  FUNCTION makeInputParameters(streamId, glb_arr_len, distribution, opt_dict) RESULT(me)
    IMPLICIT NONE
    TYPE(t_inputParameters) :: me
    INTEGER, INTENT(IN) :: streamID, glb_arr_len
    CLASS(t_scatterPattern), POINTER, INTENT(IN) :: distribution
    TYPE(t_dictionary), INTENT(IN), OPTIONAL :: opt_dict

    CHARACTER(len=*), PARAMETER :: routine = modname//':makeInputParameters'
    INTEGER :: vlistId, variableCount, i, ierrstat

    !first forward the arguments to the object we are building
    me%streamId = streamId
    me%glb_arr_len = glb_arr_len

    me%have_dict = .FALSE.
    IF(PRESENT(opt_dict)) THEN
      CALL dict_init(me%dict, lcase_sensitive=.FALSE.)
      CALL dict_copy(opt_dict, me%dict)
      me%have_dict = .TRUE.
    END IF
    me%distribution => distribution
    CALL me%distribution%resetStatistics()

    !now the interesting part: introspect the file and broadcast the variable info needed to avoid broadcasting it later.
    IF(my_process_is_stdio()) THEN
        vlistId = streamInqVlist(streamId)
        variableCount = vlistNvars(vlistId)
    END IF
    CALL p_bcast(variableCount, p_io, distribution%communicator)
    ALLOCATE(me%variableNames(variableCount), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ALLOCATE(me%variableDatatype(variableCount), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ALLOCATE(me%variableTileIdx(variableCount), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    IF(my_process_is_stdio()) THEN
        do i = 1, variableCount
            CALL vlistInqVarName(vlistId, i-1, me%variableNames(i))
            me%variableDatatype(i) = vlistInqVarDatatype(vlistId, i-1)
            me%variableTileIdx(i) = vlistInqVarIntKey(vlistId, i-1, "localInformationNumber")
        END do
    END IF
    CALL p_bcast(me%variableNames, p_io, distribution%communicator)
    CALL p_bcast(me%variableDatatype, p_io, distribution%communicator)
    CALL p_bcast(me%variableTileIdx, p_io, distribution%communicator)
    me%readDuration = 0.0_dp
    me%readBytes = 0_i8
  END FUNCTION makeInputParameters

  !---------------------------------------------------------------------------------------------------------------------------------
  !> Determine the datatype of the given variable in the input file
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER FUNCTION inputParametersFindVarId(me, name, opt_tileidx) RESULT(result)
    IMPLICIT NONE
    CLASS(t_inputParameters), INTENT(IN) :: me
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, INTENT(IN), OPTIONAL :: opt_tileidx

    CHARACTER(len=*), PARAMETER :: routine = modname//':inputParametersFindVarId'
    CHARACTER(len=DICT_MAX_STRLEN) :: mapped_name
    INTEGER :: i

    mapped_name = trim(name)
    IF(me%have_dict) THEN
      ! Search name mapping for name in NetCDF/GRIB2 file
      mapped_name = TRIM(dict_get(me%dict, TRIM(name), DEFAULT=name))
    END IF
    mapped_name = tolower(trim(mapped_name))

    result = -1
    do i = 1, size(me%variableNames)
        IF(present(opt_tileidx)) THEN
            !can't be the variable we are looking for if it's got the wrong tile number
            IF(me%variableTileIdx(i) /= opt_tileidx) cycle
        END IF
        IF(trim(tolower(trim(me%variableNames(i)))) == trim(mapped_name)) THEN
            result = i-1
            exit
        END IF
    END do

    !sanity check
    IF(result < 0) THEN
        IF(present(opt_tileidx)) THEN
            WRITE (0,*) "tileidx = ", opt_tileidx
        END IF
        CALL finish(routine, "Variable "//trim(name)//" not found!")
    END IF
  END FUNCTION inputParametersFindVarId

  !---------------------------------------------------------------------------------------------------------------------------------
  !> Lookup the datatype of a variable given its cdi varId
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER FUNCTION inputParametersLookupDatatype(me, varId) RESULT(result)
    IMPLICIT NONE
    CLASS(t_inputParameters), INTENT(IN) :: me
    INTEGER, INTENT(IN) :: varId

    result = me%variableDatatype(varId+1)
  END FUNCTION inputParametersLookupDatatype

  !---------------------------------------------------------------------------------------------------------------------------------
  !> Find a variable and determine its datatype
  !---------------------------------------------------------------------------------------------------------------------------------
  INTEGER FUNCTION inputParametersFindVarDatatype(me, name, opt_tileidx) RESULT(result)
    IMPLICIT NONE
    CLASS(t_inputParameters), INTENT(IN) :: me
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, INTENT(IN), OPTIONAL :: opt_tileidx

    result = me%lookupDatatype(me%findVarId(name, opt_tileidx))
  END FUNCTION inputParametersFindVarDatatype

  !---------------------------------------------------------------------------------------------------------------------------------
  !> Destroys a t_inputParameters object
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE deleteInputParameters(me)
    IMPLICIT NONE
    CLASS(t_inputParameters), INTENT(INOUT) :: me

    CHARACTER(len=*), PARAMETER :: routine = modname//':inputParametersFindVarId'

    IF(msg_level >= 5 .and. my_process_is_stdio()) then
        WRITE(0,*) routine, ": Total read statistics for stream ID ", me%streamId
        WRITE(0,'(8X,A,I19,A)')   "amount:    ", me%readBytes, " bytes"
        WRITE(0,'(8X,A,F19.3,A)') "duration:  ", me%readDuration, " seconds"
        WRITE(0,'(8X,A,F19.3,A)') "bandwidth: ", REAL(me%readBytes, dp)/(1048576.0_dp*me%readDuration), " MiB/s"
    END IF
    CALL me%distribution%printStatistics()

    IF(me%have_dict) CALL dict_finalize(me%dict)
    DEALLOCATE(me%variableNames)
    DEALLOCATE(me%variableDatatype)
    DEALLOCATE(me%variableTileIdx)
  END SUBROUTINE deleteInputParameters

  !-------------------------------------------------------------------------
  !> @return CDI variable ID if CDI stream contains a variable of the
  !> given name, -1 otherwise
  !
  !  Uses cdilib for file access.
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  FUNCTION test_cdi_varID(streamID, name, opt_tileidx, opt_dict) RESULT(result_varID)
    INTEGER                                 :: result_varID
    INTEGER,           INTENT(IN)           :: streamID            !< link to file 
    CHARACTER (LEN=*), INTENT(IN)           :: name                !< variable name
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx       !< tile index, encoded as "localInformationNumber"
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict          !< optional: variable name dictionary
    ! local variables
    CHARACTER(len=MAX_CHAR_LENGTH)  :: zname
    INTEGER                         :: nvars, varID, vlistID, tileidx
    CHARACTER(LEN=DICT_MAX_STRLEN)  :: mapped_name

    mapped_name = TRIM(name)
    IF (PRESENT(opt_dict)) THEN
      ! Search name mapping for name in NetCDF/GRIB2 file
      mapped_name = TRIM(dict_get(opt_dict, name, DEFAULT=name))
    END IF
    mapped_name = tolower(TRIM(mapped_name))

    zname   = ""
    vlistID = streamInqVlist(streamID)

    result_varID = -1
    ! total number of available fields:
    nvars = vlistNvars(vlistID)
    ! loop over vlist, find the corresponding varID
    LOOP : DO varID=0,(nvars-1)

      CALL vlistInqVarName(vlistID, varID, zname)
      IF (TRIM(tolower(TRIM(zname))) == TRIM(mapped_name)) THEN

        ! check tile index
        IF (PRESENT(opt_tileidx)) THEN
          tileidx = vlistInqVarIntKey(vlistID, varID, "localInformationNumber")
          IF (tileidx /= opt_tileidx) CYCLE LOOP
        END IF

        result_varID = varID
        EXIT LOOP
      END IF
    END DO LOOP
  END FUNCTION test_cdi_varID


  !-------------------------------------------------------------------------
  !> @return vlist variable ID for a given variable name
  !
  !  Uses cdilib for file access.
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  FUNCTION get_cdi_varID(streamID, name, opt_tileidx, opt_dict) RESULT(result_varID)
    INTEGER                                 :: result_varID
    INTEGER,           INTENT(IN)           :: streamID            !< link to file 
    CHARACTER (LEN=*), INTENT(IN)           :: name                !< variable name
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx       !< tile index, encoded as "localInformationNumber"
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict          !< optional: variable name dictionary
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::get_cdi_varID'

    result_varID = test_cdi_varID(streamID, name, opt_tileidx, opt_dict)
    IF (result_varID== -1) THEN
      IF (PRESENT(opt_tileidx)) THEN
        WRITE (0,*) "tileidx = ", opt_tileidx
      END IF
      CALL finish(routine, "Variable "//TRIM(name)//" not found!")
    END IF
  END FUNCTION get_cdi_varID


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for streamReadVarSlice() that measures the time
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE timeStreamReadVarSlice(parameters, varID, level, buffer, nmiss)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID, level
    REAL(dp), INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(OUT) :: nmiss

    CHARACTER(len=*), PARAMETER :: routine = modname//':timeStreamReadVarSlice'
    REAL(dp) :: startTime, duration
    INTEGER(i8) :: bytes
    REAL(dp), SAVE :: totalTime
    INTEGER(i8), SAVE :: totalBytes

    startTime = p_mpi_wtime()
    CALL streamReadVarSlice(parameters%streamID, varID, level, buffer, nmiss)
    IF(msg_level >= 5 .and. my_process_is_stdio()) THEN
        duration = p_mpi_wtime() - startTime
        bytes = INT(SIZE(buffer, 1), i8) * 8_i8
        parameters%readDuration = parameters%readDuration + duration
        parameters%readBytes = parameters%readBytes + bytes
        IF(msg_level >= 15) then
            WRITE(0,*) routine, ": Read ", bytes, " bytes in ", duration, " seconds (", &
                &      REAL(bytes, dp)/(1048576.0_dp*duration), " MiB/s)"
        END IF
    END IF
  END SUBROUTINE timeStreamReadVarSlice


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for streamReadVarSliceF() that measures the time
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE timeStreamReadVarSliceF(parameters, varID, level, buffer, nmiss)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID, level
    REAL(sp), INTENT(INOUT) :: buffer(:)
    INTEGER, INTENT(OUT) :: nmiss

    CHARACTER(len=*), PARAMETER :: routine = modname//':timeStreamReadVarSliceF'
    REAL(dp) :: startTime, duration
    INTEGER(i8) :: bytes
    REAL(dp), SAVE :: totalTime
    INTEGER(i8), SAVE :: totalBytes

    startTime = p_mpi_wtime()
    CALL streamReadVarSliceF(parameters%streamID, varID, level, buffer, nmiss)
    IF(msg_level >= 5 .and. my_process_is_stdio()) THEN
        duration = p_mpi_wtime() - startTime
        bytes = INT(SIZE(buffer, 1), i8) * 8_i8
        parameters%readDuration = parameters%readDuration + duration
        parameters%readBytes = parameters%readBytes + bytes
        IF(msg_level >= 15) then
            WRITE(0,*) routine, ": Read ", bytes, " bytes in ", duration, " seconds (", &
                &      REAL(bytes, dp)/(1048576.0_dp*duration), " MiB/s)"
        END IF
    END IF
  END SUBROUTINE timeStreamReadVarSliceF


  !---------------------------------------------------------------------------------------------------------------------------------
  !> read and distribute a 3D variable across the processes
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_cdi_3d_wp(parameters, varID, nlevs, levelDimension, var_out, lvalue_add)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID, nlevs, levelDimension
    REAL(wp), INTENT(INOUT) :: var_out(:,:,:) !< output field
    LOGICAL, INTENT(IN) :: lvalue_add         !< If .TRUE., add values to given field

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_3d_wp'
    INTEGER :: jk, ierrstat, nmiss
    REAL(wp), ALLOCATABLE :: tmp_buf(:) ! temporary local array

    ! allocate a buffer for one vertical level
    IF(my_process_is_stdio()) THEN
        ALLOCATE(tmp_buf(parameters%glb_arr_len), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF

    ! initialize output field:
    !XXX: It stands to reason whether we really need to do this. It's only effective iff there are points in var_out that are not
    !     covered by blk_no/idx_no. But the unconditional version was definitely wrong as it entirely negated the effect of
    !     lvalue_add.
    IF(.not. lvalue_add) var_out(:,:,:) = 0._wp

    !FIXME: This code is most likely latency bound, not throughput bound. Needs some asynchronicity to hide the latencies.
    DO jk=1,nlevs
      IF(my_process_is_stdio()) THEN
        ! read record as 1D field
        CALL timeStreamReadVarSlice(parameters, varID, jk-1, tmp_buf(:), nmiss)
      END IF

      SELECT CASE(levelDimension)
          CASE(2)
              CALL parameters%distribution%distribute(tmp_buf(:), var_out(:, jk, :), lvalue_add)
          CASE(3)
              CALL parameters%distribution%distribute(tmp_buf(:), var_out(:, :, jk), lvalue_add)
          CASE DEFAULT
              CALL finish(routine, "Internal error!")
      END SELECT
    END DO ! jk=1,nlevs

    ! clean up
    IF(my_process_is_stdio()) THEN
        DEALLOCATE(tmp_buf, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF

  END SUBROUTINE read_cdi_3d_wp


  !---------------------------------------------------------------------------------------------------------------------------------
  !> read and distribute a 3D variable across the processes
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_cdi_3d_sp(parameters, varID, nlevs, levelDimension, var_out, lvalue_add)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID, nlevs, levelDimension
    REAL(wp), INTENT(INOUT) :: var_out(:,:,:) !< output field
    LOGICAL, INTENT(IN) :: lvalue_add         !< If .TRUE., add values to given field

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_3d_sp'
    INTEGER :: jk, ierrstat, nmiss
    REAL(sp), ALLOCATABLE :: tmp_buf(:) ! temporary local array

    ! allocate a buffer for one vertical level
    IF(my_process_is_stdio()) THEN
        ALLOCATE(tmp_buf(parameters%glb_arr_len), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF

    ! initialize output field:
    !XXX: It stands to reason whether we really need to do this. It's only effective iff there are points in var_out that are not
    !     covered by blk_no/idx_no. But the unconditional version was definitely wrong as it entirely negated the effect of
    !     lvalue_add.
    IF(.not. lvalue_add) var_out(:,:,:) = 0._wp

    !FIXME: This code is most likely latency bound, not throughput bound. Needs some asynchronicity to hide the latencies.
    DO jk=1,nlevs
      IF(my_process_is_stdio()) THEN
        ! read record as 1D field
        CALL timeStreamReadVarSliceF(parameters, varID, jk-1, tmp_buf(:), nmiss)
      END IF

      SELECT CASE(levelDimension)
          CASE(2)
              CALL parameters%distribution%distribute(tmp_buf(:), var_out(:, jk, :), lvalue_add)
          CASE(3)
              CALL parameters%distribution%distribute(tmp_buf(:), var_out(:, :, jk), lvalue_add)
          CASE DEFAULT
              CALL finish(routine, "Internal error!")
      END SELECT
    END DO ! jk=1,nlevs

    ! clean up
    IF(my_process_is_stdio()) THEN
        DEALLOCATE(tmp_buf, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF

  END SUBROUTINE read_cdi_3d_sp

  !-------------------------------------------------------------------------
  !> Read 3D dataset from file.
  ! 
  !  Note: This implementation uses a 2D buffer.
  ! 
  !  @par Revision History
  !  Initial revision by F. Prill, DWD (2013-02-19)
  ! 
  SUBROUTINE read_cdi_3d_real(parameters, varname, nlevs, var_out, opt_tileidx, opt_lvalue_add, opt_lev_dim)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),    INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,             INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    REAL(wp),            INTENT(INOUT) :: var_out(:,:,:) !< output field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx          !< tile index, encoded as "localInformationNumber"
    LOGICAL,             INTENT(IN), OPTIONAL :: opt_lvalue_add       !< If .TRUE., add values to given field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_lev_dim          !< array dimension (of the levels)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_3d_real'
    INTEGER :: vlistId, varId, zaxisId, gridId, levelDimension
    LOGICAL :: lvalue_add

    levelDimension = 2
    CALL assign_if_present(levelDimension, opt_lev_dim)

    lvalue_add = .FALSE.
    CALL assign_if_present(lvalue_add, opt_lvalue_add)

    varId = parameters%findVarId(varname, opt_tileidx)
    IF(my_process_is_stdio()) THEN
        ! sanity check of the variable dimensions
        vlistId = streamInqVlist(parameters%streamId)
        zaxisId = vlistInqVarZaxis(vlistId, varId)
        gridId = vlistInqVarGrid(vlistId, varId)
        IF ((gridInqSize(gridId) /= parameters%glb_arr_len) .OR.  &
            & (zaxisInqSize(zaxisId) /= nlevs)) THEN
            CALL finish(routine, "Incompatible dimensions!")
        END IF
    END IF

    SELECT CASE(parameters%lookupDatatype(varId))
        CASE(DATATYPE_FLT64, DATATYPE_INT32)
            ! int32 is treated as double precision because single precision floats would cut off up to seven bits from the integer
            CALL read_cdi_3d_wp(parameters, varId, nlevs, levelDimension, var_out, lvalue_add)
        CASE DEFAULT
            CALL read_cdi_3d_sp(parameters, varId, nlevs, levelDimension, var_out, lvalue_add)
    END SELECT

  END SUBROUTINE read_cdi_3d_real


  !---------------------------------------------------------------------------------------------------------------------------------
  !> read and distribute a 2D variable across the processes
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_cdi_2d_sp(parameters, varID, var_out)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID
    REAL(wp), INTENT(INOUT) :: var_out(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_sp'
    INTEGER :: nmiss, ierrstat
    REAL(sp), ALLOCATABLE :: tmp_buf(:)

    IF (my_process_is_stdio()) THEN
        ! read record as 1D field
        ALLOCATE(tmp_buf(parameters%glb_arr_len), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        CALL timeStreamReadVarSliceF(parameters, varID, 0, tmp_buf(:), nmiss)
    END IF

    CALL parameters%distribution%distribute(tmp_buf(:), var_out(:, :), .false.)

    IF (my_process_is_stdio()) THEN
        DEALLOCATE(tmp_buf, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE read_cdi_2d_sp

  !---------------------------------------------------------------------------------------------------------------------------------
  !> read and distribute a 2D variable across the processes
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_cdi_2d_wp(parameters, varID, var_out)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: varID
    REAL(wp), INTENT(INOUT) :: var_out(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_wp'
    INTEGER :: nmiss, ierrstat
    REAL(wp), ALLOCATABLE :: tmp_buf(:)

    IF (my_process_is_stdio()) THEN
        ! read record as 1D field
        ALLOCATE(tmp_buf(parameters%glb_arr_len), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        CALL timeStreamReadVarSlice(parameters, varID, 0, tmp_buf(:), nmiss)
    END IF

    CALL parameters%distribution%distribute(tmp_buf(:), var_out(:, :), .false.)

    IF (my_process_is_stdio()) THEN
        DEALLOCATE(tmp_buf, STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    END IF
  END SUBROUTINE read_cdi_2d_wp


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  SUBROUTINE read_cdi_2d_real (parameters, varname, var_out, opt_tileidx)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*), INTENT(IN)    :: varname        !< Var name of field to be read
    REAL(wp),         INTENT(INOUT) :: var_out(:,:)   !< output field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx  !< tile index, encoded as "localInformationNumber"

    ! local variables:
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_real'
    INTEGER       :: varId, vlistId, gridId

    varId = parameters%findVarId(varname, opt_tileidx)
    IF (my_process_is_stdio()) THEN
        !sanity check on the variable dimensions
        vlistId   = streamInqVlist(parameters%streamId)
        gridId    = vlistInqVarGrid(vlistId, varId)
        IF (gridInqSize(gridId) /= parameters%glb_arr_len) CALL finish(routine, "Incompatible dimensions!")
    END IF
    SELECT CASE(parameters%lookupDatatype(varId))
        CASE(DATATYPE_FLT64, DATATYPE_INT32)
            ! int32 is treated as double precision because single precision floats would cut off up to seven bits from the integer
            CALL read_cdi_2d_wp(parameters, varId, var_out)
        CASE DEFAULT
            CALL read_cdi_2d_sp(parameters, varId, var_out)
    END SELECT
  END SUBROUTINE read_cdi_2d_real


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for INTEGER fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !

  SUBROUTINE read_cdi_2d_int(parameters, varname, var_out, opt_tileidx)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*), INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,          INTENT(INOUT) :: var_out(:,:)   !< output field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx  !< tile index, encoded as "localInformationNumber"

    ! local variables:
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_int'
    REAL(wp), ALLOCATABLE :: var_tmp(:,:)
    INTEGER               :: ierrstat

    ! allocate a temporary array:
    ALLOCATE(var_tmp(SIZE(var_out,1), SIZE(var_out,2)), STAT=ierrstat)
    ! Initialization is needed in order to avoid errors in subsequent conversion to INTEGER
    var_tmp(:,:) = 0._wp
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! read the field as a REAL-valued field:
    CALL read_cdi_2d_real (parameters, varname, var_tmp, opt_tileidx)
    ! perform number conversion
    var_out(:,:) = NINT(var_tmp)
    ! clean up
    DEALLOCATE(var_tmp, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE read_cdi_2d_int

  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  SUBROUTINE read_cdi_2d_time (parameters, ntime, varname, var_out, opt_tileidx)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER,          INTENT(IN)    :: ntime          !< time levels of file
    CHARACTER(len=*), INTENT(IN)    :: varname        !< Var name of field to be read
    REAL(wp),         INTENT(INOUT) :: var_out(:,:,:) !< output field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx  !< tile index, encoded as "localInformationNumber"

    ! local variables:
    INTEGER :: jt, nrecs, varId

    ! Get var ID
    varId = parameters%findVarId(varname, opt_tileidx)
    DO jt = 1, ntime
      IF (my_process_is_stdio()) THEN
        nrecs = streamInqTimestep(parameters%streamId, (jt-1))
      END IF
      SELECT CASE(parameters%lookupDatatype(varId))
        CASE(DATATYPE_FLT64, DATATYPE_INT32)
            ! int32 is treated as double precision because single precision floats would cut off up to seven bits from the integer
            CALL read_cdi_2d_wp(parameters, varId, var_out(:,:,jt))
        CASE DEFAULT
            CALL read_cdi_2d_sp(parameters, varId, var_out(:,:,jt))
      END SELECT
    END DO
  END SUBROUTINE read_cdi_2d_time


  !------------------------------------------------------------------------------------------------
  !> Set additional GRIB2 keys
  !
  ! GRIB2 Quick hack
  ! ----------------
  !
  ! Set additional GRIB2 keys. These are added to each single variable, even though
  ! adding it to the vertical or horizontal grid description may be more elegant.
  !
  SUBROUTINE set_additional_GRIB2_keys(vlistID, varID, grib_conf, tileidx)
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf
    INTEGER,                INTENT(IN) :: tileidx

    CHARACTER(len=MAX_CHAR_LENGTH)     :: ydate, ytime
    INTEGER :: cent, year, month, day    ! date
    INTEGER :: hour, minute, second      ! time
    INTEGER :: steptype

    IF (grib_conf%ldate_grib_act) THEN
      ! get date and time
      ! ydate : ccyymmdd, ytime : hhmmss.sss
      CALL date_and_time(ydate,ytime)
      READ(ydate,'(4i2)') cent, year, month, day
      READ(ytime,'(3i2)') hour, minute, second
    ELSE ! set date to "01010101" (for better comparability of GRIB files)
      cent  = 1
      year  = 1
      month = 1
      year  = 1
      hour  = 1
      minute= 1
      second= 1
    ENDIF

    ! inquire steptype (needed below)
    steptype = vlistInqVarTsteptype(vlistID, varID)

    ! Load correct tables and activate section 2
    !
    ! set tablesVersion=11
    CALL vlistDefVarIntKey(vlistID, varID, "tablesVersion", 11)
    ! Initialize section 2
    CALL vlistDefVarIntKey(vlistID, varID, "grib2LocalSectionPresent", 1)
    !
    ! Product definition
    CALL vlistDefVarIntKey(vlistID, varID, "significanceOfReferenceTime",     &
      &                    grib_conf%significanceOfReferenceTime)
    CALL vlistDefVarIntKey(vlistID, varID, "productionStatusOfProcessedData", &
      &                    grib_conf%productionStatusOfProcessedData)
    CALL vlistDefVarIntKey(vlistID, varID, "typeOfProcessedData",             &
      &                    grib_conf%typeOfProcessedData)
    CALL vlistDefVarIntKey(vlistID, varID, "backgroundProcess",               &
      &                    grib_conf%backgroundProcess)
    ! in case of lon-lat output, "1" has to be added to generatingProcessIdentifier
    CALL vlistDefVarIntKey(vlistID, varID, "generatingProcessIdentifier",     &
      &                    grib_conf%generatingProcessIdentifier)


    IF (grib_conf%lspecialdate_invar) THEN
      ! Use special date for invariant and climatological fields
      !
      IF ( steptype == TSTEP_CONSTANT ) THEN
        ! invariant data 
        CALL vlistDefVarTypeOfGeneratingProcess(vlistID, varID, 196)
      ELSE
        CALL vlistDefVarTypeOfGeneratingProcess(vlistID, varID, grib_conf%typeOfGeneratingProcess)
      ENDIF
    ELSE
      ! no special treatment of invariant and climatological fields
      !
      CALL vlistDefVarTypeOfGeneratingProcess(vlistID, varID, grib_conf%typeOfGeneratingProcess)
    ENDIF


    ! Product Generation (local), !! DWD only !!
    ! DWD :: center    = 78
    !        subcenter = 255
    IF ((grib_conf%generatingCenter == 78) .AND. (grib_conf%generatingSubcenter == 255)) THEN
      CALL vlistDefVarIntKey(vlistID, varID, "localDefinitionNumber"  ,         &
        &                    grib_conf%localDefinitionNumber)

      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateYear"  , 100*cent+year)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateMonth" , month)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateDay"   , day)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateHour"  , hour)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateMinute", minute)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateSecond", second)
      ! CALL vlistDefVarIntKey(vlistID, varID, "localValidityDateYear"  , 2013)

      ! preliminary HACK for identifying tile based variables
      CALL vlistDefVarIntKey(vlistID, varID, "localNumberOfExperiment",                 &
        &                    grib_conf%localNumberOfExperiment)

      CALL vlistDefVarIntKey(vlistID, varID, "localInformationNumber" , tileidx)

      IF (grib_conf%localDefinitionNumber == 254) THEN
        !
        ! -------------------------------------------
        ! Local definition for deterministic forecast
        ! -------------------------------------------

        ! store GRIB_API library version
        ! TODO: replace by wrapper call in cdi not available yet (2013-07-29)
        !CALL vlistDefVarIntKey(vlistID, varID, "localVersionNumber" , gribGetAPIVersion())

        !
      ELSE IF (grib_conf%localDefinitionNumber == 253) THEN
        !
        ! --------------------------------------
        ! Local definition for ensemble products
        ! --------------------------------------

        ! Switch between template 
        ! 1: Individual ensemble forecast, control and perturbed, at a horizontal level 
        !    or in a horizontal layer at a point in time.
        ! and 
        ! 11: Individual ensemble forecast, control and perturbed, at a horizontal level 
        !     or in a horizontal layer, in a continuous or non-continuous time interval
        ! depending on variable-specific stepType
        !!! Note that this is only a workaround, since we cannot guarantee that !!!
        !!! no information (keys) set prior to that template-switch gets lost.  !!!
!DR Actually, this should be removed, since CDI is taking care of template switching.
!DR However, for backward compatibility, we keep the following lines for a couple of months.
!DR Otherwise, people would be forced to use the grib-API prerelease 1.12.3 when doing 
!DR ensemble runs. 
        IF ( ANY((/TSTEP_INSTANT,TSTEP_CONSTANT/) == steptype) ) THEN 
          CALL vlistDefVarIntKey(vlistID, varID, "productDefinitionTemplateNumber",1)
        ELSE
          CALL vlistDefVarIntKey(vlistID, varID, "productDefinitionTemplateNumber",11)
        ENDIF


        IF (grib_conf%typeOfEnsembleForecast /= -1)                                       &
          &   CALL vlistDefVarIntKey(vlistID, varID, "typeOfEnsembleForecast" ,           &
          &                          grib_conf%typeOfEnsembleForecast)
        IF (grib_conf%localTypeOfEnsembleForecast /= -1)                                  &
          &   CALL vlistDefVarIntKey(vlistID, varID, "localTypeOfEnsembleForecast" ,      &
          &                          grib_conf%localTypeOfEnsembleForecast)
        IF (grib_conf%numberOfForecastsInEnsemble /= -1)                                  &
          &   CALL vlistDefVarIntKey(vlistID, varID, "numberOfForecastsInEnsemble" ,      &
          &                          grib_conf%numberOfForecastsInEnsemble)
        IF (grib_conf%perturbationNumber /= -1)                                           &
          &   CALL vlistDefVarIntKey(vlistID, varID, "perturbationNumber" ,               &
          &                          grib_conf%perturbationNumber)
      END IF ! localDefinitionNumber
    END IF

    ! SECTION 3
    CALL vlistDefVarIntKey(vlistID, varID, "shapeOfTheEarth", 6)

  END SUBROUTINE set_additional_GRIB2_keys


  !------------------------------------------------------------------------------------------------
  !> Set additional, time-dependent GRIB2 keys
  !!  
  !!  This subroutine sets all GRIB2 keys that may change during the
  !!  simulation. Currently this is the case for the start and end time 
  !!  of the statistical process interval.
  !!
  !!  Description of how the GRIB2 keys 'forecastTime' and 'lengthOfTimeRange' are computed.
  !!
  !!  |===================*======================|================================> time axis
  !!  ^start_time         |                      ^current_time (output)
  !!                      |                      |
  !!                      ^EventLastTriggerTime  |
  !!                      |                      |
  !!                      |                      |
  !!  <-------------------><--------------------->
  !!     forecastTime         lengthOfTimeRange
  !!
  !!  @author D. Reinert, F. Prill (DWD)
  !!
  !!  TODOs
  !!
  !!  - we implicitly assume that we always need to INQUIRE the interval from action 1.
  !!  - we implicitly assume that the time range is always smaller than one month
  !!
  SUBROUTINE set_timedependent_GRIB2_keys(streamID, varID, info, start_date, cur_date)
    INTEGER                             ,INTENT(IN) :: streamID, varID
    TYPE (t_var_metadata)               ,INTENT(IN) :: info
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) ,INTENT(IN) :: start_date, cur_date
    ! local variables
    TYPE(timedelta), POINTER      :: mtime_lengthOfTimeRange       ! length of time range (mtime format)
    INTEGER                       :: ilengthOfTimeRange_secs       ! length of time range in seconds
    INTEGER                       :: ilengthOfTimeRange            ! length of time range in appropriate time units
    INTEGER                       :: taxis_tunit                   ! time axis units
    INTEGER                       :: vlistID, taxisID
    INTEGER                       :: indicatorOfUnitForTimeRange
    TYPE(timedelta), POINTER      :: forecast_delta                ! forecast time (mtime format)
    INTEGER                       :: forecast_secs                 ! forecast time in seconds
    INTEGER                       :: forecast_time                 ! forecast time in appropriate time units
    TYPE(datetime),  POINTER      :: mtime_start                   ! model start (initialization) time
    TYPE(datetime),  POINTER      :: mtime_cur                     ! model current time (rounded)
    TYPE(datetime)                :: statProc_startDateTime        ! statistical process starting DateTime

    !---------------------------------------------------------
    ! Set time-dependent metainfo
    !---------------------------------------------------------
    !
    ! Skip inapplicable fields
    IF (ALL((/TSTEP_MAX, TSTEP_MIN/) /= info%isteptype) ) RETURN

    ! get vlistID. Note that the stream-internal vlistID must be used. 
    ! It is obtained via streamInqVlist(streamID)
    vlistID     = streamInqVlist(streamID)
    taxisID     = vlistInqTaxis(vlistID)
    taxis_tunit = taxisInqTunit(taxisID)

    ! get current DateTime (rounded)
    mtime_cur   => newDatetime(TRIM(cur_date))
    ! get model start date
    mtime_start => newDatetime(TRIM(start_date))

    IF (info%action_list%n_actions > 0) THEN

      ! here we implicitly assume, that action Nr. 1 is the 'nullify'-action.
      ! Currently this is true, but must not be true for all times. A less 
      ! ad-hoc solution would be desirable. A warning is issued, if nr. 1 
      ! is NOT the 'nullify' action.
      IF (info%action_list%action(1)%actionID /= ACTION_RESET) THEN
        write(0,*) 'set_timedependent_GRIB2_keys: actionID of action 1 is not equal to ACTION_RESET.'//&
          &             'lengthOfTimeRange may not be set correctly'
      ENDIF

      ! get latest (intended) triggering time, which is equivalent to 
      ! the statistical process starting time
      statProc_startDateTime = info%action_list%action(1)%EventLastTriggerDate


      ! get time interval, over which statistical process has been performed
      ! It is the time difference between the current time (rounded) mtime_cur and 
      ! the last time the nullify-action took place (rounded) statProc_startDateTime.
      mtime_lengthOfTimeRange  => newTimedelta("P01D")  ! init
      !
      ! mtime_lengthOfTimeRange = mtime_cur - statProc_startDateTime
      CALL getTimeDeltaFromDateTime(mtime_cur, statProc_startDateTime, mtime_lengthOfTimeRange)


      ! time interval over which statistical process has been performed (in secs)    
      ilengthOfTimeRange_secs = 86400 *INT(mtime_lengthOfTimeRange%day)    &
           &                  + 3600  *INT(mtime_lengthOfTimeRange%hour)   &
           &                  + 60    *INT(mtime_lengthOfTimeRange%minute) &
           &                  +        INT(mtime_lengthOfTimeRange%second) 
           

      ! cleanup
      CALL deallocateTimedelta(mtime_lengthOfTimeRange)
    ELSE
      ilengthOfTimeRange_secs = 0
    END IF


    ! get forecast_time: forecast_time = statProc_startDateTime - model_startDateTime
    ! Note that for statistical quantities, the forecast time is the time elapsed between the 
    ! model start time and the start time of the statistical process
    forecast_delta => newTimedelta("P01D")
    CALL getTimeDeltaFromDateTime(statProc_startDateTime, mtime_start, forecast_delta)

    ! forecast time in seconds
    forecast_secs =    forecast_delta%second    +   &
      &             60*(forecast_delta%minute   +   & 
      &                 60*(forecast_delta%hour +   &
      &                       24*forecast_delta%day))


    SELECT CASE (taxis_tunit)
    CASE (TUNIT_SECOND)
       indicatorOfUnitForTimeRange = 13
       ilengthOfTimeRange          = ilengthOfTimeRange_secs
       forecast_time               = forecast_secs
    CASE (TUNIT_MINUTE)
       indicatorOfUnitForTimeRange = 0
       ilengthOfTimeRange          = INT(ilengthOfTimeRange_secs/60)
       forecast_time               = INT(forecast_secs/60)
    CASE (TUNIT_HOUR)
       indicatorOfUnitForTimeRange = 1
       ilengthOfTimeRange          = INT(ilengthOfTimeRange_secs/3600)
       forecast_time               = INT(forecast_secs/3600)
    CASE DEFAULT
    END SELECT

    ! set forecast time: statProc_startDateTime - model_startDateTime
    CALL vlistDefVarIntKey(vlistID, varID, "forecastTime", forecast_time) 

    ! set Indicator of unit for time range over which statistical processing is done
    ! equal to the Indicator of unit of time range. This is not time dependent and may be 
    ! moved to some better place. Note that the CDI-internal numbers differ from the official 
    ! WMO numbers!
    CALL vlistDefVarIntKey(vlistID, varID, "indicatorOfUnitForTimeRange", indicatorOfUnitForTimeRange)
    !
    ! set length of time range: current time - statProc_startDateTime
    CALL vlistDefVarIntKey(vlistID, varID, "lengthOfTimeRange",           ilengthOfTimeRange)


    ! cleanup
    CALL deallocateDatetime(mtime_start)
    CALL deallocateDatetime(mtime_cur)
    CALL deallocateTimedelta(forecast_delta)

  END SUBROUTINE set_timedependent_GRIB2_keys

END MODULE mo_util_cdi
