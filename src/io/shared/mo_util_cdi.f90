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
  USE mo_communication,      ONLY: idx_no, blk_no, t_scatterPattern
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_run_config,         ONLY: msg_level
  USE mo_mpi,                ONLY: p_bcast, p_comm_work, p_comm_work_test,  &
    &                              p_io, my_process_is_stdio, p_mpi_wtime
  USE mo_util_string,        ONLY: tolower, toupper, one_of
  USE mo_fortran_tools,      ONLY: assign_if_present
  USE mo_dictionary,         ONLY: t_dictionary, dict_get, dict_init, dict_copy, dict_finalize, DICT_MAX_STRLEN
  USE mo_gribout_config,     ONLY: t_gribout_config
  USE mo_lnd_nwp_config,     ONLY: getNumberOfTiles, tileid_int2grib, t_tile
  USE mo_var_metadata_types, ONLY: t_var_metadata, CLASS_TILE, CLASS_TILE_LAND
  USE mo_action,             ONLY: ACTION_RESET
  ! calendar operations
  USE mtime,                 ONLY: timedelta, newTimedelta,                 &
    &                              datetime, newDatetime,                   &
    &                              deallocateTimedelta, deallocateDatetime, &
    &                              MAX_DATETIME_STR_LEN
  USE mo_mtime_extensions,   ONLY: getTimeDeltaFromDateTime
  USE mo_cdi_constants,      ONLY: FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4, streamInqVlist, &
    &                              vlistNvars, vlistInqVarDatatype, vlistInqVarIntKey,      &
    &                              vlistInqVarTypeOfGeneratingProcess,                      &
    &                              vlistInqVarZaxis, zaxisInqType, ZAXIS_REFERENCE,         &
    &                              zaxisInqNlevRef, vlistInqVarGrid, gridInqSize,           &
    &                              zaxisInqSize, DATATYPE_FLT64, DATATYPE_INT32,            &
    &                              streamInqTimestep, vlistInqVarTsteptype, TSTEP_CONSTANT, &
    &                              TSTEP_INSTANT, TSTEP_AVG, TSTEP_ACCUM, TSTEP_MAX,        &
    &                              TSTEP_MIN, vlistInqTaxis, taxisInqTunit,                 &
    &                              TUNIT_SECOND, TUNIT_MINUTE, TUNIT_HOUR

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: has_filetype_netcdf
  PUBLIC :: read_cdi_2d, read_cdi_3d
  PUBLIC :: test_cdi_varID
  PUBLIC :: get_cdi_varID
  PUBLIC :: get_cdi_NlevRef
  PUBLIC :: set_GRIB2_additional_keys
  PUBLIC :: set_GRIB2_ensemble_keys
  PUBLIC :: set_GRIB2_local_keys
  PUBLIC :: set_GRIB2_tile_keys
  PUBLIC :: set_GRIB2_timedep_keys
  PUBLIC :: set_GRIB2_timedep_local_keys
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

  !> Provides a common interface to all NetCDF flavors used within
  !> ICON.
  !
  FUNCTION has_filetype_netcdf(cdi_filetype)
    LOGICAL :: has_filetype_netcdf
    INTEGER, INTENT(IN) :: cdi_filetype

    has_filetype_netcdf = ((cdi_filetype == FILETYPE_NC)  .OR. &
      &                    (cdi_filetype == FILETYPE_NC2) .OR. &
      &                    (cdi_filetype == FILETYPE_NC4))
  END FUNCTION has_filetype_netcdf


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



  !-------------------------------------------------------------------------
  !> @return the number of half levels of a generalized Z-axis for a given variable name
  !
  !  Uses cdilib for file access.
  !  Initial revision by D. Reinert, DWD (2014-10-24)
  !
  FUNCTION get_cdi_NlevRef(streamID, name, opt_tileidx, opt_dict) RESULT(result_NlevRef)
    INTEGER                                 :: result_NlevRef
    INTEGER,           INTENT(IN)           :: streamID            !< link to file 
    CHARACTER (LEN=*), INTENT(IN)           :: name                !< variable name
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx       !< tile index, encoded as "localInformationNumber"
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict          !< optional: variable name dictionary
    ! local variables
    INTEGER :: varID                                              ! variable ID
    INTEGER :: vlistID
    INTEGER :: zaxisID 
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::get_cdi_NlevRef'

    varID = test_cdi_varID(streamID, name, opt_tileidx, opt_dict)
    IF (varID== -1) THEN
      IF (PRESENT(opt_tileidx)) THEN
        WRITE (0,*) "tileidx = ", opt_tileidx
      END IF
      CALL finish(routine, "Variable "//TRIM(name)//" not found!")
    END IF
    vlistID = streamInqVlist(streamID)
    zaxisID = vlistInqVarZaxis(vlistID,varID)
    IF (zaxisInqType(zaxisID) /= ZAXIS_REFERENCE) THEN
      CALL finish(routine, "Variable "//TRIM(name)//" has no generalized Z-axis!")
    ENDIF
    ! number of half levels of the generalized Z-axis
    result_NlevRef = zaxisInqNlevRef(zaxisID)
  END FUNCTION get_cdi_NlevRef


  !---------------------------------------------------------------------------------------------------------------------------------
  !> small helper message to output the speed of an individual input operation
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE writeSpeedMessage(routineName, bytes, duration)
    CHARACTER(len=*), INTENT(IN) :: routineName
    INTEGER(i8), VALUE :: bytes
    REAL(dp), VALUE :: duration

    IF(msg_level >= 15) THEN
        IF(duration /= 0.0) THEN
            WRITE(0,*) routineName, ": Read ", bytes, " bytes in ", duration, " seconds (", &
                &      REAL(bytes, dp)/(1048576.0_dp*duration), " MiB/s)"
        ELSE
            WRITE(0,*) routineName, ": Read ", bytes, " bytes in no time. ", &
                &      "Failure to measure time may be due to compilation without MPI."
        END IF
    END IF
  END SUBROUTINE


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
        CALL writeSpeedMessage(routine, bytes, duration)
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
        CALL writeSpeedMessage(routine, bytes, duration)
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
    ELSE
        ALLOCATE(tmp_buf(1), STAT=ierrstat)
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
    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

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
    ELSE
        ALLOCATE(tmp_buf(1), STAT=ierrstat)
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
    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

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
    ELSE
        ALLOCATE(tmp_buf(1), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF

    CALL parameters%distribution%distribute(tmp_buf(:), var_out(:, :), .false.)

    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
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
    ELSE
        ALLOCATE(tmp_buf(1), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    END IF

    CALL parameters%distribution%distribute(tmp_buf(:), var_out(:, :), .false.)

    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
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
  SUBROUTINE set_GRIB2_additional_keys(vlistID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf
    !
    ! Local
    INTEGER :: steptype


    ! inquire steptype (needed below)
    steptype = vlistInqVarTsteptype(vlistID, varID)


    ! SECTION 1: Identification Section
    !
    ! Load correct tables
    !
    ! set tablesVersion=11
    CALL vlistDefVarIntKey(vlistID, varID, "tablesVersion", 11)
    !
    CALL vlistDefVarIntKey(vlistID, varID, "significanceOfReferenceTime",     &
      &                    grib_conf%significanceOfReferenceTime)
    CALL vlistDefVarIntKey(vlistID, varID, "productionStatusOfProcessedData", &
      &                    grib_conf%productionStatusOfProcessedData)
    CALL vlistDefVarIntKey(vlistID, varID, "typeOfProcessedData",             &
      &                    grib_conf%typeOfProcessedData)


    ! SECTION 3: Grid definition Section
    CALL vlistDefVarIntKey(vlistID, varID, "shapeOfTheEarth", 6)



    ! SECTION 4: Product definition
    CALL vlistDefVarIntKey(vlistID, varID, "backgroundProcess",               &
      &                    grib_conf%backgroundProcess)
    ! in case of lon-lat output, "1" has to be added to generatingProcessIdentifier
    CALL vlistDefVarIntKey(vlistID, varID, "generatingProcessIdentifier",     &
      &                    grib_conf%generatingProcessIdentifier)


    IF (ANY((/TSTEP_AVG,TSTEP_ACCUM,TSTEP_MAX,TSTEP_MIN/) == steptype)) THEN
      ! Always set
      !   typeOfTimeIncrement = 2 
      !   "Successive times processed have same start time of forecast, 
      !    forecast time is incremented"
      ! since this is the only type of time processing available in ICON
      CALL vlistDefVarIntKey(vlistID, varID, "typeOfTimeIncrement", 2)
    ENDIF

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



  END SUBROUTINE set_GRIB2_additional_keys



  !------------------------------------------------------------------------------------------------

  !>
  !! Set local GRIB2 keys (SECTION 2)
  !!
  !! Set DWD specific local keys (SECTION 2)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-20)
  !! Modification by Daniel Reinert, DWD (2015-01-20)
  !! - extract local settings from routine set_additional_GRIB2_keys 
  !!
  SUBROUTINE set_GRIB2_local_keys(vlistID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf

    ! Local
    CHARACTER(len=MAX_CHAR_LENGTH)     :: ydate, ytime
    INTEGER :: cent, year, month, day    ! date
    INTEGER :: hour, minute, second      ! time
  !----------------------------------------------------------------

    ! SECTION 2: Initialize local use section
    CALL vlistDefVarIntKey(vlistID, varID, "grib2LocalSectionPresent", 1)


    ! SECTION 2: DWD specific settings (local use)
    !
    ! DWD :: center    = 78
    !        subcenter = 255
    IF ((grib_conf%generatingCenter == 78) .AND. (grib_conf%generatingSubcenter == 255)) THEN

      CALL vlistDefVarIntKey(vlistID, varID, "localDefinitionNumber"  ,         &
        &                    grib_conf%localDefinitionNumber)

      CALL vlistDefVarIntKey(vlistID, varID, "localNumberOfExperiment",         &
        &                    grib_conf%localNumberOfExperiment)



      IF (grib_conf%localDefinitionNumber == 254) THEN
        !
        ! -------------------------------------------
        ! Local definition for deterministic forecast
        ! -------------------------------------------

      ELSE IF (grib_conf%localDefinitionNumber == 253) THEN
        !
        ! --------------------------------------
        ! Local definition for ensemble products
        ! --------------------------------------

        IF (grib_conf%localTypeOfEnsembleForecast /= -1)                                  &
          &   CALL vlistDefVarIntKey(vlistID, varID, "localTypeOfEnsembleForecast" ,      &
          &                          grib_conf%localTypeOfEnsembleForecast)

      END IF ! localDefinitionNumber
    END IF

  END SUBROUTINE set_GRIB2_local_keys




  !------------------------------------------------------------------------------------------------

  !>
  !! Set GRIB2 ensemble keys
  !!
  !! Set GRIB2 ensemble keys (SECTION 4)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-20)
  !! Modification by Daniel Reinert, DWD (2015-01-20)
  !! - extract ensemble settings from routine set_additional_GRIB2_keys
  !!
  !! ATTENTION: To be called AFTER set_GRIB2_additional_keys
  !!            due to its dependency on typeOfGeneratingProcess that may 
  !!            be changed for invariant fields in set_GRIB2_additional_keys
  !!
  SUBROUTINE set_GRIB2_ensemble_keys(vlistID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf

    ! Local
    INTEGER  :: typeOfGeneratingProcess 
  !----------------------------------------------------------------

    ! get typeOfGeneratingProcess
    ! We do not make use of grib_conf%typeOfGeneratingProcess, since 
    ! typeOfGeneratingProcess is modified for invariant fields.
    typeOfGeneratingProcess = vlistInqVarTypeOfGeneratingProcess(vlistID, varID)


    IF (typeOfGeneratingProcess == 4) THEN  ! Ensemble forecast

      ! SECTION 4: Product definition Section

      IF (grib_conf%typeOfEnsembleForecast /= -1)                                       &
        &   CALL vlistDefVarIntKey(vlistID, varID, "typeOfEnsembleForecast" ,           &
        &                          grib_conf%typeOfEnsembleForecast)
      IF (grib_conf%numberOfForecastsInEnsemble /= -1)                                  &
        &   CALL vlistDefVarIntKey(vlistID, varID, "numberOfForecastsInEnsemble" ,      &
        &                          grib_conf%numberOfForecastsInEnsemble)
      IF (grib_conf%perturbationNumber /= -1)                                           &
        &   CALL vlistDefVarIntKey(vlistID, varID, "perturbationNumber" ,               &
        &                          grib_conf%perturbationNumber)
    END IF ! typeOfGeneratingProcess

  END SUBROUTINE set_GRIB2_ensemble_keys


  !>
  !! Set tile-specific keys
  !!
  !! Set GRIB2 keys which are specific to tile-variables.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-20)
  !!
  SUBROUTINE set_GRIB2_tile_keys (vlistID, varID, info, i_lctype)

    INTEGER,                INTENT(IN) :: vlistID, varID
    TYPE (t_var_metadata),  INTENT(IN) :: info
    INTEGER,                INTENT(IN) :: i_lctype  !< Tile classification

    ! local
    INTEGER :: typeOfGeneratingProcess
    INTEGER :: productDefinitionTemplate        ! Tile template number 

    TYPE(t_tile) :: tileID

  !----------------------------------------------------------------

    ! Skip inapplicable fields
    IF ( ALL((/CLASS_TILE,CLASS_TILE_LAND/) /= info%var_class) ) RETURN

    typeOfGeneratingProcess = vlistInqVarTypeOfGeneratingProcess(vlistID, varID)

    ! change product definition template
    !
    IF (typeOfGeneratingProcess == 4) THEN
      ! ensemble
      productDefinitionTemplate = 40056     ! temporary
    ELSE
      ! deterministic
      productDefinitionTemplate = 40055     ! temporary
    ENDIF
    CALL vlistDefVarProductDefinitionTemplate(vlistID, varID, productDefinitionTemplate)

    ! Set tile classification
    CALL vlistDefVarIntKey(vlistID, varID, "tileClassification" , i_lctype)

    ! Set number of used tiles
    CALL vlistDefVarIntKey(vlistID, varID, "numberOfTiles" , &
      &                    getNumberOfTiles(class=info%var_class))

    ! get the following attributes:
    ! - identificationNumberOfTile
    ! - numberOfAttributes
    ! - identificationNumberOfAttribute
    !
    ! Convert internal tile ID into GRIB2 tile ID
    tileID = tileid_int2grib(info%ncontained)

    ! Set tile index
    CALL vlistDefVarIntKey(vlistID, varID, "identificationNumberOfTile" , tileID%GRIB2_tile%itn)

    ! Set total number of attributes for given tile
    CALL vlistDefVarIntKey(vlistID, varID, "numberOfAttributes" , tileID%GRIB2_tile%nat)

    ! Set attribute index
!!$    CALL vlistDefVarIntKey(vlistID, varID, "identificationNumberOfAttribute" , tileID%GRIB2_att%iat)
    ! misuse identificationNumberOfAttribute for storing the total number of records.
    CALL vlistDefVarIntKey(vlistID, varID, "identificationNumberOfAttribute" , info%maxcontained)

    ! Set attribute
    CALL vlistDefVarIntKey(vlistID, varID, "attribute" , tileID%GRIB2_att%attribute)

  END SUBROUTINE set_GRIB2_tile_keys


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
  SUBROUTINE set_GRIB2_timedep_keys(streamID, varID, info, start_date, cur_date)
    INTEGER                             ,INTENT(IN) :: streamID, varID
    TYPE (t_var_metadata)               ,INTENT(IN) :: info
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) ,INTENT(IN) :: start_date, cur_date
    ! local variables
    TYPE(timedelta), POINTER      :: mtime_lengthOfTimeRange       ! length of time range (mtime format)
    INTEGER                       :: ilengthOfTimeRange_secs       ! length of time range in seconds
    INTEGER                       :: ilengthOfTimeRange            ! length of time range in appropriate time units
    INTEGER                       :: taxis_tunit                   ! time axis units
    INTEGER                       :: vlistID, taxisID
    TYPE(timedelta), POINTER      :: forecast_delta                ! forecast time (mtime format)
    INTEGER                       :: forecast_secs                 ! forecast time in seconds
    INTEGER                       :: forecast_time                 ! forecast time in appropriate time units
    TYPE(datetime),  POINTER      :: mtime_start                   ! model start (initialization) time
    TYPE(datetime),  POINTER      :: mtime_cur                     ! model current time (rounded)
    TYPE(datetime)                :: statProc_startDateTime        ! statistical process starting DateTime

    ! special fields for which time-dependent metainfos should be set even though they are not of 
    ! steptype TSTEP_MAX or TSTEP_MIN. These fields are special in the sense that averaging is not 
    ! performed over the entire model run but over only some intervals.
    CHARACTER(LEN=8) :: ana_avg_vars(5) = (/"u_avg   ", "v_avg   ", "pres_avg", "temp_avg", "qv_avg  "/)

    !---------------------------------------------------------
    ! Set time-dependent metainfo
    !---------------------------------------------------------
    !
    ! Skip inapplicable fields
    ! Currently all TSTEP_AVG and TSTEP_ACC fields are skipped, except for special ones 
    ! listed in ana_avg_vars
    IF ((ALL((/TSTEP_MAX, TSTEP_MIN/) /= info%isteptype)) .AND. &
      & (one_of(TRIM(info%name),ana_avg_vars) == -1) ) RETURN


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
       ilengthOfTimeRange          = ilengthOfTimeRange_secs
       forecast_time               = forecast_secs
    CASE (TUNIT_MINUTE)
       ilengthOfTimeRange          = INT(ilengthOfTimeRange_secs/60)
       forecast_time               = INT(forecast_secs/60)
    CASE (TUNIT_HOUR)
       ilengthOfTimeRange          = INT(ilengthOfTimeRange_secs/3600)
       forecast_time               = INT(forecast_secs/3600)
    CASE DEFAULT
    END SELECT

    ! set forecast time: statProc_startDateTime - model_startDateTime
    CALL vlistDefVarIntKey(vlistID, varID, "forecastTime", forecast_time) 

    !
    ! set length of time range: current time - statProc_startDateTime
    CALL vlistDefVarIntKey(vlistID, varID, "lengthOfTimeRange",           ilengthOfTimeRange)
    ! Note that if one of the statistics templates 4.8 or 4.11 is selected, the time unit 
    ! (GRIB2 key "indicatorOfUnitForTimeRange") is set automatically by CDI.
    ! It is always set identical to "indicatorOfUnitOFTimeRange"


    ! cleanup
    CALL deallocateDatetime(mtime_start)
    CALL deallocateDatetime(mtime_cur)
    CALL deallocateTimedelta(forecast_delta)

  END SUBROUTINE set_GRIB2_timedep_keys


  !>
  !! Set time-dependent local GRIB2 keys (SECTION 2)
  !!
  !! Set DWD specific time-dependent local keys (SECTION 2)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-02-06)
  !!
  SUBROUTINE set_GRIB2_timedep_local_keys(streamID, varID, grib_conf)
    INTEGER,                INTENT(IN) :: streamID, varID
    TYPE(t_gribout_config), INTENT(IN) :: grib_conf

    ! Local
    CHARACTER(len=MAX_CHAR_LENGTH)     :: ydate, ytime
    INTEGER :: cent, year, month, day    ! date
    INTEGER :: hour, minute, second      ! time
    INTEGER :: vlistID
  !----------------------------------------------------------------


    ! get vlistID. Note that the stream-internal vlistID must be used. 
    ! It is obtained via streamInqVlist(streamID)
    vlistID = streamInqVlist(streamID)

    ! SECTION 2: DWD specific settings (local use)
    !
    ! DWD :: center    = 78
    !        subcenter = 255
    IF ((grib_conf%generatingCenter == 78) .AND. (grib_conf%generatingSubcenter == 255)) THEN

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
        hour  = 1
        minute= 1
        second= 1
      ENDIF


      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateYear"  , 100*cent+year)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateMonth" , month)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateDay"   , day)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateHour"  , hour)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateMinute", minute)
      CALL vlistDefVarIntKey(vlistID, varID, "localCreationDateSecond", second)
      ! CALL vlistDefVarIntKey(vlistID, varID, "localValidityDateYear"  , 2013)

    END IF

  END SUBROUTINE set_GRIB2_timedep_local_keys

END MODULE mo_util_cdi
