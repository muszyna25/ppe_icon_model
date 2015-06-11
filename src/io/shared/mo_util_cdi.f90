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
  USE mo_communication,      ONLY: t_scatterPattern
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_run_config,         ONLY: msg_level
  USE mo_mpi,                ONLY: p_bcast, p_io, my_process_is_stdio, p_mpi_wtime,  &
    &                              my_process_is_mpi_workroot
  USE mo_util_string,        ONLY: tolower
  USE mo_fortran_tools,      ONLY: assign_if_present
  USE mo_dictionary,         ONLY: t_dictionary, dict_get, dict_init, dict_copy, dict_finalize, DICT_MAX_STRLEN
  USE mo_cdi,                ONLY: FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4, streamInqVlist, vlistNvars, vlistInqVarDatatype, &
                                 & vlistInqVarIntKey, vlistInqVarZaxis, zaxisInqType, ZAXIS_REFERENCE, zaxisInqNlevRef, &
                                 & vlistInqVarGrid, gridInqSize, zaxisInqSize, DATATYPE_FLT64, DATATYPE_INT32, streamInqTimestep, &
                                 & vlistInqVarTsteptype, TSTEP_CONSTANT, TSTEP_INSTANT, TSTEP_MAX, TSTEP_MIN, vlistInqTaxis, &
                                 & taxisInqTunit, TUNIT_SECOND, TUNIT_MINUTE, TUNIT_HOUR, vlistDefVarIntKey, &
                                 & vlistDefVarTypeOfGeneratingProcess, streamReadVarSliceF, streamReadVarSlice, vlistInqVarName, &
                                 & TSTEP_AVG,TSTEP_ACCUM,TSTEP_MAX,TSTEP_MIN, vlistInqVarSubtype, subtypeInqSize

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: has_filetype_netcdf
  PUBLIC :: read_cdi_2d, read_cdi_3d
  PUBLIC :: test_cdi_varID
  PUBLIC :: get_cdi_varID
  PUBLIC :: get_cdi_NlevRef
  PUBLIC :: t_inputParameters, makeInputParameters, deleteInputParameters
  PUBLIC :: t_tileinfo_elt, trivial_tileinfo

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_cdi'

  INTERFACE read_cdi_2d
    MODULE PROCEDURE read_cdi_2d_real
    MODULE PROCEDURE read_cdi_2d_int
    MODULE PROCEDURE read_cdi_2d_time
    MODULE PROCEDURE read_cdi_2d_real_tiles
    MODULE PROCEDURE read_cdi_2d_int_tiles
    MODULE PROCEDURE read_cdi_2d_time_tiles
  END INTERFACE

  INTERFACE read_cdi_3d
    MODULE PROCEDURE read_cdi_3d_real
    MODULE PROCEDURE read_cdi_3d_real_tiles
  END INTERFACE

  TYPE t_tileinfo_elt
    INTEGER :: idx  !< variable specific tile index
    INTEGER :: att  !< variable specific tile attribute
  END TYPE t_tileinfo_elt

  TYPE t_tileinfo
    TYPE(t_tileinfo_elt), ALLOCATABLE :: tile(:)  !< variable specific
    !tile indices > CDI internal, one-dimensional index for the list
    !of (idx/att) > pairs:
    INTEGER,              ALLOCATABLE :: tile_index(:)
  END TYPE t_tileinfo

  ! trivial tile information, denoting a "no-tile" field:
  TYPE(t_tileinfo_elt), PARAMETER :: trivial_tileinfo = t_tileinfo_elt(idx=-99, att=-99)

  ! This is a small type that serves two functions:
  ! 1. It encapsulates three parameters to the read functions into one, significantly reducing the hassle to call them.
  ! 2. It allows some information about the file to be distributed to the non-I/O PEs in a single broadcast that would otherwise
  !    require a broadcast in every call.
  !
  ! Create this by calling makeInputParameters(), destroy with deleteInputParameters().
  TYPE :: t_inputParameters
    !PRIVATE! Don't use outside of this file!
    INTEGER                                     :: streamId     !< CDI stream ID
    INTEGER                                     :: glb_arr_len  !< global array length
    TYPE (t_dictionary)                         :: dict     !< a dictionary that is used to translate variable names
    CLASS(t_scatterPattern), POINTER            :: distribution    !< a t_scatterPattern to distribute the data
    LOGICAL                                     :: have_dict    !< whether `dict` is used or not
    CHARACTER(LEN=MAX_CHAR_LENGTH), ALLOCATABLE :: variableNames(:) !< the names of the variables
    INTEGER, ALLOCATABLE                        :: variableDatatype(:) !< the datatypes that are used by the file to store the variables
    TYPE(t_tileinfo), ALLOCATABLE               :: variableTileinfo(:) !< variable specific tile meta info

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
    INTEGER :: subtypeID
    INTEGER :: ientry
    INTEGER :: cnt

    INTEGER, ALLOCATABLE :: tileIdx_container(:)
    INTEGER, ALLOCATABLE :: tileAtt_container(:)
    INTEGER, ALLOCATABLE :: tileTid_container(:)
    INTEGER, ALLOCATABLE :: subtypeSize(:)

    INTEGER :: idx, att, tile_index


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
    IF(my_process_is_mpi_workroot()) THEN
        vlistId = streamInqVlist(streamId)
        variableCount = vlistNvars(vlistId)
    END IF
    CALL p_bcast(variableCount, p_io, distribution%communicator)

    ALLOCATE(me%variableNames(variableCount), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ALLOCATE(me%variableDatatype(variableCount), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ALLOCATE(me%variableTileinfo(variableCount), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! local arrays
    ALLOCATE(subtypeSize(variableCount), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    IF(my_process_is_mpi_workroot()) THEN
        subtypeSize(1:variableCount) = 0
        do i = 1, variableCount
            CALL vlistInqVarName(vlistId, i-1, me%variableNames(i))
            me%variableDatatype(i) = vlistInqVarDatatype(vlistId, i-1)

            subtypeID = vlistInqVarSubtype(vlistID,i-1)
            subtypeSize(i) = subtypeInqSize(subtypeID)
            !
            ALLOCATE(me%variableTileinfo(i)%tile(subtypeSize(i)), &
              &      me%variableTileinfo(i)%tile_index(subtypeSize(i)), STAT=ierrstat)
            IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

            IF (vlistInqVarIntKey(vlistId, i-1, "totalNumberOfTileAttributePairs") <= 0) THEN
              ! not a tile variable
              me%variableTileinfo(i)%tile(:)       = trivial_tileinfo
              me%variableTileinfo(i)%tile_index(:) = -99
            ELSE
              ! tile
              DO ientry=1, subtypeSize(i)
                CALL subtypeDefActiveIndex(subtypeID,ientry-1)  ! starts with 0

                idx        = vlistInqVarIntKey(vlistId, i-1, "tileIndex")
                att        = vlistInqVarIntKey(vlistId, i-1, "tileAttribute")
                tile_index = ientry-1

                me%variableTileinfo(i)%tile(ientry)       = t_tileinfo_elt( idx, att )
                me%variableTileinfo(i)%tile_index(ientry) = tile_index
                ! reset active index
                CALL subtypeDefActiveIndex(subtypeID,0)
              END DO
            ENDIF  ! totalNumberOfTileAttributePairs <= 0

        END do

    END IF

    CALL p_bcast(subtypeSize, p_io, distribution%communicator)

    ! put tile info into local 1D arrays, for broadcasting
    ! currently direct bradcast of variable variableTileinfo not possible
    ALLOCATE(tileIdx_container(SUM(subtypeSize(1:variableCount))), &
      &      tileAtt_container(SUM(subtypeSize(1:variableCount))), &
      &      tileTId_container(SUM(subtypeSize(1:variableCount))), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    IF(my_process_is_mpi_workroot()) THEN
      cnt = 0
      DO i=1, variableCount
        tileIdx_container(cnt+1:cnt+subtypeSize(i)) = me%variableTileinfo(i)%tile(1:subtypeSize(i))%idx
        tileAtt_container(cnt+1:cnt+subtypeSize(i)) = me%variableTileinfo(i)%tile(1:subtypeSize(i))%att
        tileTid_container(cnt+1:cnt+subtypeSize(i)) = me%variableTileinfo(i)%tile_index(1:subtypeSize(i))
        cnt = cnt + subtypeSize(i)
      END DO
    ENDIF

    CALL p_bcast(me%variableNames   , p_io, distribution%communicator)
    CALL p_bcast(me%variableDatatype, p_io, distribution%communicator)
    CALL p_bcast(tileIdx_container  , p_io, distribution%communicator)
    CALL p_bcast(tileAtt_container  , p_io, distribution%communicator)
    CALL p_bcast(tileTid_container  , p_io, distribution%communicator)

    ! read tile info from broadcasted 1D array and store in array variableTileinfo of TYPE t_tileinfo
    IF (.NOT. my_process_is_mpi_workroot()) THEN
      cnt = 0
      DO i=1, variableCount
        IF (subtypeSize(i) > 0) THEN
          ALLOCATE(me%variableTileinfo(i)%tile(subtypeSize(i)), &
            &      me%variableTileinfo(i)%tile_index(subtypeSize(i)), STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
          me%variableTileinfo(i)%tile(1:subtypeSize(i))%idx        = &
            &      tileIdx_container(cnt+1:cnt+subtypeSize(i))
          me%variableTileinfo(i)%tile(1:subtypeSize(i))%att        = &
            &      tileAtt_container(cnt+1:cnt+subtypeSize(i))
          me%variableTileinfo(i)%tile_index(1:subtypeSize(i))      = &
            &      tileTid_container(cnt+1:cnt+subtypeSize(i))
          cnt = cnt + subtypeSize(i)
        ENDIF
      ENDDO
    ENDIF

    me%readDuration = 0.0_dp
    me%readBytes = 0_i8

    ! cleanup
    DEALLOCATE(tileIdx_container, tileAtt_container, tileTid_container, subtypeSize, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END FUNCTION makeInputParameters


  !---------------------------------------------------------------------------------------------------------------------------------
  !> Determine the datatype of the given variable in the input file
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE inputParametersFindVarId(me, name, tileinfo, varID, tile_index) 
    IMPLICIT NONE
    CLASS(t_inputParameters), INTENT(IN)  :: me
    CHARACTER(len=*),         INTENT(IN)  :: name
    TYPE(t_tileinfo_elt),     INTENT(IN)  :: tileinfo
    INTEGER,                  INTENT(OUT) :: varID, tile_index

    CHARACTER(len=*), PARAMETER :: routine = modname//':inputParametersFindVarId'
    CHARACTER(len=DICT_MAX_STRLEN) :: mapped_name
    INTEGER :: i, j

    mapped_name = trim(name)
    IF(me%have_dict) THEN
      ! Search name mapping for name in NetCDF/GRIB2 file
      mapped_name = TRIM(dict_get(me%dict, TRIM(name), DEFAULT=name))
    END IF
    mapped_name = tolower(trim(mapped_name))

    varID      = -1
    tile_index = -1
    DO i = 1, SIZE(me%variableNames)
      DO j = 1, SIZE(me%variableTileinfo(i)%tile)
        IF ((TRIM(tolower(TRIM(me%variableNames(i)))) == TRIM(mapped_name)) .AND.  & 
          & (me%variableTileinfo(i)%tile(j)%idx == tileinfo%idx)            .AND.  &
          & (me%variableTileinfo(i)%tile(j)%att == tileinfo%att)) THEN
          varID      = i-1
          tile_index = me%variableTileinfo(i)%tile_index(j)
          EXIT
        END IF
      END DO
    END DO


    ! insanity check
    IF(varID < 0) THEN
      CALL finish(routine, "Variable "//TRIM(name)//" not found!")
    END IF
  END SUBROUTINE inputParametersFindVarId

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
  INTEGER FUNCTION inputParametersFindVarDatatype(me, name, tileinfo) RESULT(result)
    IMPLICIT NONE
    CLASS(t_inputParameters), INTENT(IN) :: me
    CHARACTER(len=*),         INTENT(IN) :: name
    TYPE(t_tileinfo_elt),     INTENT(IN)    :: tileinfo
    ! local variables
    INTEGER :: varID, tile_index
    CALL me%findVarId(name, tileinfo, varID, tile_index)
    result = me%lookupDatatype(varID)
  END FUNCTION inputParametersFindVarDatatype


  !---------------------------------------------------------------------------------------------------------------------------------
  !> Destroys a t_inputParameters object
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE deleteInputParameters(me)
    IMPLICIT NONE
    CLASS(t_inputParameters), INTENT(INOUT) :: me

    CHARACTER(len=*), PARAMETER :: routine = modname//':deleteInputParameters'

    IF(msg_level >= 10 .and. my_process_is_stdio()) then
        WRITE(0,*) routine, ": Total read statistics for stream ID ", me%streamId
        WRITE(0,'(8X,A,I19,A)')   "amount:    ", me%readBytes, " bytes"
        WRITE(0,'(8X,A,F19.3,A)') "duration:  ", me%readDuration, " seconds"
        IF (me%readDuration > 0._wp) THEN
          WRITE(0,'(8X,A,F19.3,A)') "bandwidth: ", REAL(me%readBytes, dp)/(1048576.0_dp*me%readDuration), " MiB/s"
        ENDIF
    END IF
    CALL me%distribution%printStatistics()

    IF(me%have_dict) CALL dict_finalize(me%dict)
    DEALLOCATE(me%variableNames)
    DEALLOCATE(me%variableDatatype)
    DEALLOCATE(me%variableTileinfo)
  END SUBROUTINE deleteInputParameters

  !-------------------------------------------------------------------------
  !> @return CDI variable ID if CDI stream contains a variable of the
  !> given name, -1 otherwise
  !
  !  Uses cdilib for file access.
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  FUNCTION test_cdi_varID(streamID, name, opt_dict) RESULT(result_varID)
    INTEGER                                 :: result_varID
    INTEGER,           INTENT(IN)           :: streamID            !< link to file 
    CHARACTER (LEN=*), INTENT(IN)           :: name                !< variable name
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict          !< optional: variable name dictionary
    ! local variables
    CHARACTER(len=MAX_CHAR_LENGTH)  :: zname
    INTEGER                         :: nvars, varID, vlistID
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
  FUNCTION get_cdi_varID(streamID, name, opt_dict) RESULT(result_varID)
    INTEGER                                 :: result_varID
    INTEGER,           INTENT(IN)           :: streamID            !< link to file 
    CHARACTER (LEN=*), INTENT(IN)           :: name                !< variable name
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict          !< optional: variable name dictionary
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::get_cdi_varID'

    result_varID = test_cdi_varID(streamID, name, opt_dict)
    IF (result_varID== -1) THEN
      CALL finish(routine, "Variable "//TRIM(name)//" not found!")
    END IF
  END FUNCTION get_cdi_varID



  !-------------------------------------------------------------------------
  !> @return the number of half levels of a generalized Z-axis for a given variable name
  !
  !  Uses cdilib for file access.
  !  Initial revision by D. Reinert, DWD (2014-10-24)
  !
  FUNCTION get_cdi_NlevRef(streamID, name, opt_dict) RESULT(result_NlevRef)
    INTEGER                                 :: result_NlevRef
    INTEGER,           INTENT(IN)           :: streamID            !< link to file 
    CHARACTER (LEN=*), INTENT(IN)           :: name                !< variable name
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict          !< optional: variable name dictionary
    ! local variables
    INTEGER :: varID                                              ! variable ID
    INTEGER :: vlistID
    INTEGER :: zaxisID 
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::get_cdi_NlevRef'

    varID = test_cdi_varID(streamID, name, opt_dict)
    IF (varID== -1) THEN
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
    IF(my_process_is_mpi_workroot()) THEN
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
      IF(my_process_is_mpi_workroot()) THEN
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
    IF(my_process_is_mpi_workroot()) THEN
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
      IF(my_process_is_mpi_workroot()) THEN
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
  SUBROUTINE read_cdi_3d_real_tiles(parameters, varname, nlevs, var_out, tileinfo, opt_lvalue_add, opt_lev_dim)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,                 INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    REAL(wp),                INTENT(INOUT) :: var_out(:,:,:) !< output field
    TYPE(t_tileinfo_elt),    INTENT(IN)    :: tileinfo
    LOGICAL,             INTENT(IN), OPTIONAL :: opt_lvalue_add       !< If .TRUE., add values to given field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_lev_dim          !< array dimension (of the levels)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_3d_real_tiles'
    INTEGER :: vlistId, varId, zaxisId, gridId, levelDimension, subtypeID, tile_index
    LOGICAL :: lvalue_add

    levelDimension = 2
    CALL assign_if_present(levelDimension, opt_lev_dim)

    lvalue_add = .FALSE.
    CALL assign_if_present(lvalue_add, opt_lvalue_add)

    CALL parameters%findVarId(varname, tileinfo, varID, tile_index)

    IF ((tile_index < 0) .AND. (tileinfo%idx >= 0)) THEN
      CALL finish(routine, "Requested tile not found!")
    END IF

    IF(my_process_is_mpi_workroot()) THEN
      ! set active tile index, if this is a tile-based variable
      vlistId = streamInqVlist(parameters%streamId)
      subtypeID  = vlistInqVarSubtype(vlistID,varID)
      IF (tile_index >= 0)  CALL subtypeDefActiveIndex(subtypeID, tile_index)
      ! sanity check of the variable dimensions
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

    IF(my_process_is_mpi_workroot()) THEN
      ! reset tile index
      IF (tile_index >= 0)  CALL subtypeDefActiveIndex(subtypeID, 0)
    END IF
  END SUBROUTINE read_cdi_3d_real_tiles


  !-------------------------------------------------------------------------
  !> Read 3D dataset from file.
  ! 
  !  Note: This implementation uses a 2D buffer.
  ! 
  !  @par Revision History
  !  Initial revision by F. Prill, DWD (2013-02-19)
  ! 
  SUBROUTINE read_cdi_3d_real(parameters, varname, nlevs, var_out, opt_lvalue_add, opt_lev_dim)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,                 INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    REAL(wp),                INTENT(INOUT) :: var_out(:,:,:) !< output field
    LOGICAL,             INTENT(IN), OPTIONAL :: opt_lvalue_add       !< If .TRUE., add values to given field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_lev_dim          !< array dimension (of the levels)

    CALL read_cdi_3d_real_tiles(parameters, varname, nlevs, var_out, &
      &                         trivial_tileinfo, opt_lvalue_add, opt_lev_dim)
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

    IF (my_process_is_mpi_workroot()) THEN
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

    IF (my_process_is_mpi_workroot()) THEN
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
  SUBROUTINE read_cdi_2d_real_tiles (parameters, varname, var_out, tileinfo)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    REAL(wp),                INTENT(INOUT) :: var_out(:,:)   !< output field
    TYPE(t_tileinfo_elt),    INTENT(IN)    :: tileinfo

    ! local variables:
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_real_tiles'
    INTEGER       :: varId, vlistId, gridId, subtypeID, tile_index

    CALL parameters%findVarId(varname, tileinfo, varID, tile_index)

    IF ((tile_index < 0) .AND. (tileinfo%idx >= 0)) THEN
      CALL finish(routine, "Requested tile not found!")
    END IF

    IF (my_process_is_mpi_workroot()) THEN
      ! set active tile index, if this is a tile-based variable
      vlistId   = streamInqVlist(parameters%streamId)
      subtypeID  = vlistInqVarSubtype(vlistID,varID)
      IF (tile_index >= 0)  CALL subtypeDefActiveIndex(subtypeID, tile_index)
      !sanity check on the variable dimensions
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

    IF(my_process_is_mpi_workroot()) THEN
      ! reset tile index
      IF (tile_index >= 0)  CALL subtypeDefActiveIndex(subtypeID, 0)
    END IF
  END SUBROUTINE read_cdi_2d_real_tiles


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  SUBROUTINE read_cdi_2d_real (parameters, varname, var_out)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    REAL(wp),                INTENT(INOUT) :: var_out(:,:)   !< output field

    CALL read_cdi_2d_real_tiles (parameters, varname, var_out, trivial_tileinfo)
  END SUBROUTINE read_cdi_2d_real


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for INTEGER fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !

  SUBROUTINE read_cdi_2d_int_tiles(parameters, varname, var_out, tileinfo)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,                 INTENT(INOUT) :: var_out(:,:)   !< output field
    TYPE(t_tileinfo_elt),    INTENT(IN)    :: tileinfo

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
    CALL read_cdi_2d_real_tiles (parameters, varname, var_tmp, tileinfo )
    ! perform number conversion
    var_out(:,:) = NINT(var_tmp)
    ! clean up
    DEALLOCATE(var_tmp, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE read_cdi_2d_int_tiles


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for INTEGER fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !

  SUBROUTINE read_cdi_2d_int(parameters, varname, var_out)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,                 INTENT(INOUT) :: var_out(:,:)   !< output field

    CALL read_cdi_2d_int_tiles(parameters, varname, var_out, trivial_tileinfo)
  END SUBROUTINE read_cdi_2d_int


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  SUBROUTINE read_cdi_2d_time_tiles (parameters, ntime, varname, var_out, tileinfo)
    TYPE(t_inputParameters), INTENT(INOUT)  :: parameters
    INTEGER,                 INTENT(IN)     :: ntime          !< time levels of file
    CHARACTER(len=*),        INTENT(IN)     :: varname        !< Var name of field to be read
    REAL(wp),                INTENT(INOUT)  :: var_out(:,:,:) !< output field
    TYPE(t_tileinfo_elt),    INTENT(IN)     :: tileinfo

    ! local variables:
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_cdi_2d_time_tiles'
    INTEGER :: jt, nrecs, varId, subtypeID, tile_index, vlistID

    ! Get var ID
    CALL parameters%findVarId(varname, tileinfo, varID, tile_index)

    IF ((tile_index < 0) .AND. (tileinfo%idx >= 0)) THEN
      CALL finish(routine, "Requested tile not found!")
    END IF

    DO jt = 1, ntime
      IF (my_process_is_mpi_workroot()) THEN
        ! set active tile index, if this is a tile-based variable
        vlistId    = streamInqVlist(parameters%streamId)
        subtypeID  = vlistInqVarSubtype(vlistID, varID)
        IF (tile_index >= 0)  CALL subtypeDefActiveIndex(subtypeID, tile_index)
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
    IF(my_process_is_mpi_workroot()) THEN
      ! reset tile index
      IF (tile_index >= 0)  CALL subtypeDefActiveIndex(subtypeID, 0)
    END IF
  END SUBROUTINE read_cdi_2d_time_tiles


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  SUBROUTINE read_cdi_2d_time (parameters, ntime, varname, var_out)
    TYPE(t_inputParameters), INTENT(INOUT) :: parameters
    INTEGER,                 INTENT(IN)    :: ntime          !< time levels of file
    CHARACTER(len=*),        INTENT(IN)    :: varname        !< Var name of field to be read
    REAL(wp),                INTENT(INOUT) :: var_out(:,:,:) !< output field

    CALL read_cdi_2d_time_tiles (parameters, ntime, varname, var_out, trivial_tileinfo)
  END SUBROUTINE read_cdi_2d_time

END MODULE mo_util_cdi
