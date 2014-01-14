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
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_util_cdi

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_mpi,                ONLY: my_process_is_stdio, p_bcast,           &
    &                              p_comm_work, p_comm_work_test,          &
    &                              p_io, p_pe
  USE mo_util_string,        ONLY: tolower
  USE mo_fortran_tools,      ONLY: assign_if_present
  USE mo_dictionary,         ONLY: t_dictionary, dict_get, DICT_MAX_STRLEN

  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PRIVATE

  PUBLIC  :: read_cdi_2d, read_cdi_3d
  PUBLIC  :: get_cdi_varID
  PUBLIC  :: test_cdi_varID

  CHARACTER(len=*), PARAMETER :: version = &
    &    '$Id$'
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
    &                     nlevs, var_out, opt_tileidx, opt_lvalue_add, opt_dict)

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
    ! local constants:
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':read_cdi_3d'
    ! local variables:
    INTEGER                         :: vlistID, varID, zaxisID, gridID,   &
      &                                mpi_comm, j, jl, jb, jk, ierrstat, &
      &                                dimlen(3), nmiss
    REAL(wp), ALLOCATABLE           :: tmp_buf(:) ! temporary local array
    LOGICAL                         :: lvalue_add

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

END MODULE mo_util_cdi
