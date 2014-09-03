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

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_mpi,                ONLY: p_bcast, p_comm_work, p_comm_work_test,  &
    &                              p_io, p_pe
  USE mo_util_string,        ONLY: tolower
  USE mo_fortran_tools,      ONLY: assign_if_present
  USE mo_dictionary,         ONLY: t_dictionary, dict_get, DICT_MAX_STRLEN
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

  PUBLIC  :: read_cdi_2d, read_cdi_3d
  PUBLIC  :: get_cdi_varID
  PUBLIC  :: test_cdi_varID
  PUBLIC  :: set_additional_GRIB2_keys
  PUBLIC  :: set_timedependent_GRIB2_keys

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_cdi'

  INTERFACE read_cdi_2d
    MODULE PROCEDURE read_cdi_2d_real_id
    MODULE PROCEDURE read_cdi_2d_real
    MODULE PROCEDURE read_cdi_2d_int
    MODULE PROCEDURE read_cdi_2d_time
  END INTERFACE

CONTAINS

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
    LOGICAL                         :: l_found
    INTEGER                         :: nvars, varID, vlistID, tileidx
    CHARACTER(LEN=DICT_MAX_STRLEN)  :: mapped_name

    mapped_name = TRIM(name)
    IF (PRESENT(opt_dict)) THEN
      ! Search name mapping for name in NetCDF/GRIB2 file
      mapped_name = TRIM(dict_get(opt_dict, name, default=name))
    END IF
    mapped_name = tolower(TRIM(mapped_name))

    zname   = ""
    vlistID = streamInqVlist(streamID)

    result_varID = -1
    ! total number of available fields:
    nvars = vlistNvars(vlistID)
    ! loop over vlist, find the corresponding varID
    l_found = .FALSE.
    LOOP : DO varID=0,(nvars-1)

      CALL vlistInqVarName(vlistID, varID, zname)
      IF (TRIM(tolower(TRIM(zname))) == TRIM(mapped_name)) THEN

        ! check tile index
        IF (PRESENT(opt_tileidx)) THEN
          tileidx = vlistInqVarIntKey(vlistID, varID, "localInformationNumber")
          IF (tileidx /= opt_tileidx) CYCLE LOOP
        END IF

        result_varID = varID
        l_found = .TRUE.
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
  !> Read 3D dataset from file.
  ! 
  !  Note: This implementation uses a 2D buffer.
  ! 
  !  @par Revision History
  !  Initial revision by F. Prill, DWD (2013-02-19)
  ! 
  SUBROUTINE read_cdi_3d(streamID, varname, glb_arr_len, loc_arr_len, glb_index, &
    &                     nlevs, var_out, opt_tileidx, opt_lvalue_add, opt_dict, opt_lev_dim)

    INTEGER,             INTENT(IN)    :: streamID       !< ID of CDI file stream
    CHARACTER(len=*),    INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,             INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    INTEGER,             INTENT(IN)    :: glb_arr_len    !< length of 1D field (global)
    INTEGER,             INTENT(IN)    :: loc_arr_len    !< length of 1D field (local)
    INTEGER,             INTENT(IN)    :: glb_index(:)   !< Index mapping local to global
    REAL(wp),            INTENT(INOUT) :: var_out(:,:,:) !< output field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx          !< tile index, encoded as "localInformationNumber"
    LOGICAL,             INTENT(IN), OPTIONAL :: opt_lvalue_add       !< If .TRUE., add values to given field
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict             !< optional: variable name dictionary
    INTEGER,             INTENT(IN), OPTIONAL :: opt_lev_dim          !< array dimension (of the levels)
    ! local constants:
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':read_cdi_3d'
    ! local variables:
    INTEGER                         :: vlistID, varID, zaxisID, gridID,   &
      &                                mpi_comm, j, jl, jb, jk, ierrstat, &
      &                                dimlen(3), nmiss
    REAL(wp), ALLOCATABLE           :: tmp_buf(:) ! temporary local array
    LOGICAL                         :: lvalue_add
    INTEGER                         :: lev_dim

    lev_dim = 2
    CALL assign_if_present(lev_dim, opt_lev_dim)

    lvalue_add = .FALSE.
    CALL assign_if_present(lvalue_add, opt_lvalue_add)

    ! allocate a buffer for one vertical level
    ALLOCATE(tmp_buf(glb_arr_len), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! get var ID
    IF(p_pe == p_io) THEN
      vlistID   = streamInqVlist(streamID)
      varID     = get_cdi_varID(streamID, name=TRIM(varname), opt_tileidx=opt_tileidx, opt_dict=opt_dict)
      zaxisID   = vlistInqVarZaxis(vlistID, varID)
      gridID    = vlistInqVarGrid(vlistID, varID)
      dimlen(1) = gridInqSize(gridID)
      dimlen(2) = zaxisInqSize(zaxisID)

      ! Check variable dimensions:
      IF ((dimlen(1) /= glb_arr_len) .OR.  &
        & (dimlen(2) /= nlevs)) THEN
        CALL finish(routine, "Incompatible dimensions!")
      END IF
    END IF

    ! initialize output field:
    var_out(:,:,:) = 0._wp

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data
   
    DO jk=1,nlevs
      IF(p_pe == p_io) THEN
        ! read record as 1D field
        CALL streamReadVarSlice(streamID, varID, jk-1, tmp_buf(:), nmiss)
      END IF
      
      ! broadcast data: 
      CALL p_bcast(tmp_buf, p_io, mpi_comm)
      ! Set var_out from global data

      SELECT CASE(lev_dim)
      CASE(2) 
        IF (lvalue_add) THEN
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jk,jb) = var_out(jl,jk,jb) + REAL(tmp_buf(glb_index(j)), wp)
          ENDDO
        ELSE
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jk,jb) = REAL(tmp_buf(glb_index(j)), wp)
          ENDDO
        END IF
      CASE(3)
        IF (lvalue_add) THEN
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jb,jk) = var_out(jl,jb,jk) + REAL(tmp_buf(glb_index(j)), wp)
          ENDDO
        ELSE
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jb,jk) = REAL(tmp_buf(glb_index(j)), wp)
          ENDDO
        END IF
      CASE DEFAULT
        CALL finish(routine, "Internal error!")
      END SELECT
    END DO ! jk=1,nlevs
      
    ! clean up
    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE read_cdi_3d


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  SUBROUTINE read_cdi_2d_real_id (streamID, varID, glb_arr_len, loc_arr_len, glb_index, var_out)

    INTEGER,          INTENT(IN)    :: streamID       !< ID of CDI file stream
    INTEGER,          INTENT(IN)    :: varID          !< ID of CDI variable
    INTEGER,          INTENT(IN)    :: glb_arr_len    !< length of 1D field (global)
    INTEGER,          INTENT(IN)    :: loc_arr_len    !< length of 1D field (local)
    INTEGER,          INTENT(IN)    :: glb_index(:)   !< Index mapping local to global
    REAL(wp),         INTENT(INOUT) :: var_out(:,:)   !< output field
    ! local variables:
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':read_cdi_2d_real_id'
    INTEGER       :: mpi_comm, j, jl, jb, nmiss
    REAL(wp)      :: z_dummy_array(glb_arr_len)       !< local dummy array

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF (p_pe == p_io) THEN
      ! read record as 1D field
      CALL streamReadVarSlice(streamID, varID, 0, z_dummy_array(:), nmiss)
    END IF
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:) = 0._wp

    ! Set var_out from global data
    DO j = 1, loc_arr_len
      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch
      var_out(jl,jb) = z_dummy_array(glb_index(j))
    ENDDO
  END SUBROUTINE read_cdi_2d_real_id


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  SUBROUTINE read_cdi_2d_real (streamID, varname, glb_arr_len, loc_arr_len, glb_index, var_out, &
    &                          opt_tileidx, opt_dict)
    INTEGER,          INTENT(IN)    :: streamID       !< ID of CDI file stream
    CHARACTER(len=*), INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,          INTENT(IN)    :: glb_arr_len    !< length of 1D field (global)
    INTEGER,          INTENT(IN)    :: loc_arr_len    !< length of 1D field (local)
    INTEGER,          INTENT(IN)    :: glb_index(:)   !< Index mapping local to global
    REAL(wp),         INTENT(INOUT) :: var_out(:,:)   !< output field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx  !< tile index, encoded as "localInformationNumber"
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict     !< optional: variable name dictionary
    ! local variables:
    CHARACTER(len=max_char_length), PARAMETER :: routine = modname//':read_cdi_2d_real'
    INTEGER       :: varID, vlistID, gridID

    ! Get var ID
    IF (p_pe == p_io) THEN
      vlistID   = streamInqVlist(streamID)
      varID     = get_cdi_varID(streamID, name=TRIM(varname), opt_tileidx=opt_tileidx, opt_dict=opt_dict)
      gridID    = vlistInqVarGrid(vlistID, varID)
      ! Check variable dimensions:
      IF (gridInqSize(gridID) /= glb_arr_len) &
        &  CALL finish(routine, "Incompatible dimensions!")
    END IF
    CALL read_cdi_2d_real_id (streamID, varID, glb_arr_len, loc_arr_len, glb_index, var_out)
  END SUBROUTINE read_cdi_2d_real


  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for INTEGER fields
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  SUBROUTINE read_cdi_2d_int (streamID, varname, glb_arr_len, loc_arr_len, glb_index, var_out, &
    &                         opt_tileidx, opt_dict)

    INTEGER,          INTENT(IN)    :: streamID       !< ID of CDI file stream
    CHARACTER(len=*), INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,          INTENT(IN)    :: glb_arr_len    !< length of 1D field (global)
    INTEGER,          INTENT(IN)    :: loc_arr_len    !< length of 1D field (local)
    INTEGER,          INTENT(IN)    :: glb_index(:)   !< Index mapping local to global
    INTEGER,          INTENT(INOUT) :: var_out(:,:)   !< output field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx  !< tile index, encoded as "localInformationNumber"
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict     !< optional: variable name dictionary
    ! local variables:
    CHARACTER(len=max_char_length), PARAMETER :: routine = modname//':read_cdi_2d_int'
    REAL(wp), ALLOCATABLE :: var_tmp(:,:)
    INTEGER               :: ierrstat

    ! allocate a temporary array:
    ALLOCATE(var_tmp(SIZE(var_out,1), SIZE(var_out,2)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! read the field as a REAL-valued field:
    CALL read_cdi_2d_real (streamID, varname, glb_arr_len, loc_arr_len, glb_index, var_tmp, &
      &                    opt_tileidx, opt_dict=opt_dict)
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
  SUBROUTINE read_cdi_2d_time (streamID, ntime, varname, glb_arr_len, loc_arr_len, glb_index, var_out, &
    &                          opt_tileidx, opt_dict)

    INTEGER,          INTENT(IN)    :: streamID       !< ID of CDI file stream
    INTEGER,          INTENT(IN)    :: ntime          !< time levels of file
    CHARACTER(len=*), INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,          INTENT(IN)    :: glb_arr_len    !< length of 1D field (global)
    INTEGER,          INTENT(IN)    :: loc_arr_len    !< length of 1D field (local)
    INTEGER,          INTENT(IN)    :: glb_index(:)   !< Index mapping local to global
    REAL(wp),         INTENT(INOUT) :: var_out(:,:,:) !< output field
    INTEGER,             INTENT(IN), OPTIONAL :: opt_tileidx  !< tile index, encoded as "localInformationNumber"
    TYPE (t_dictionary), INTENT(IN), OPTIONAL :: opt_dict     !< optional: variable name dictionary
    ! local variables:
    CHARACTER(len=max_char_length), PARAMETER :: routine = modname//':read_cdi_2d_time'
    INTEGER :: jt, nrecs, vlistID, varID

    ! Get var ID
    IF (p_pe == p_io) THEN
      vlistID   = streamInqVlist(streamID)
      varID     = get_cdi_varID(streamID, name=TRIM(varname), opt_tileidx=opt_tileidx, &
        &                       opt_dict=opt_dict)
    END IF

    DO jt = 1, ntime
      IF (p_pe == p_io) THEN
        nrecs = streamInqTimestep(streamID, (jt-1))
      END IF
      CALL read_cdi_2d (streamID, varID, glb_arr_len, loc_arr_len, glb_index, var_out(:,:,jt))
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
