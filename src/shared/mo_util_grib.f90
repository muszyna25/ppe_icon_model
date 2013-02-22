!>
!! Contains routines for reading GRIB2 files (using the CDI library).
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
MODULE mo_util_grib

  USE mo_kind,               ONLY: wp, sp
  USE mo_exception,          ONLY: finish
  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_mpi,                ONLY: my_process_is_stdio, p_bcast,  &
    &                              p_comm_work, p_comm_work_test, &
    &                              p_io, p_pe
  USE mo_util_string,        ONLY: tolower
  USE mo_fortran_tools,      ONLY: assign_if_present

  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PRIVATE

  PUBLIC  :: read_grib_2d, read_grib_3d
  PUBLIC  :: get_varID


  CHARACTER(len=*), PARAMETER :: version = &
    &    '$Id$'
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_grib'

CONTAINS

  !-------------------------------------------------------------------------
  !> @return vlist variable ID for a given variable name
  !
  !  Uses cdilib for file access.
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  FUNCTION get_varID(streamID, name, opt_tileidx) RESULT(result_varID)
    INTEGER                                 :: result_varID
    INTEGER,           INTENT(IN)           :: streamID            !< link to GRIB file 
    CHARACTER (LEN=*), INTENT(IN)           :: name                !< variable name
    INTEGER,           INTENT(IN), OPTIONAL :: opt_tileidx         !< tile index, encoded as "localInformationNumber"
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::get_varID'
    CHARACTER(len=MAX_CHAR_LENGTH) :: zname
    LOGICAL                        :: l_found
    INTEGER                        :: nvars, varID, vlistID, tileidx

    zname   = ""
    vlistID = streamInqVlist(streamID)

    result_varID = -1
    ! total number of available fields:
    nvars = vlistNvars(vlistID)
    ! loop over vlist, find the corresponding varID
    l_found = .FALSE.
    LOOP : DO varID=0,(nvars-1)

      CALL vlistInqVarName(vlistID, varID, zname)


      IF (tolower(TRIM(zname)) == tolower(TRIM(name))) THEN

        ! check tile index
        IF (PRESENT(opt_tileidx)) THEN
          CALL vlistInqVarRawBegin(streamID, varID)
          tileidx = vlistInqVarIntKey(streamID, "localInformationNumber")
          CALL vlistInqVarRawEnd(streamID)
          IF (tileidx /= opt_tileidx) CYCLE LOOP
        END IF

        result_varID = varID
        l_found = .TRUE.
        EXIT LOOP
      END IF
    END DO LOOP
    IF (.NOT. l_found) THEN
      if (present(opt_tileidx)) then
        write (0,*) "tileidx = ", opt_tileidx
      end if
      CALL finish(routine, "Variable "//TRIM(name)//" not found!")
    END IF
  END FUNCTION get_varID


  !-------------------------------------------------------------------------
  !> Read 3D dataset from GRIB2 file.
  ! 
  !  Note: This implementation uses a 2D buffer.
  ! 
  !  @par Revision History
  !  Initial revision by F. Prill, DWD (2013-02-19)
  ! 
  SUBROUTINE read_grib_3d(streamID, varname, glb_arr_len, loc_arr_len, glb_index, &
    &                     nlevs, var_out, opt_tileidx, opt_lvalue_add)

    INTEGER,          INTENT(IN)    :: streamID       !< ID of CDI file stream
    CHARACTER(len=*), INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,          INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    INTEGER,          INTENT(IN)    :: glb_arr_len    !< length of 1D field (global)
    INTEGER,          INTENT(IN)    :: loc_arr_len    !< length of 1D field (local)
    INTEGER,          INTENT(IN)    :: glb_index(:)   !< Index mapping local to global
    INTEGER,          INTENT(IN), OPTIONAL :: opt_tileidx  !< tile index, encoded as "localInformationNumber"
    REAL(wp),         INTENT(INOUT) :: var_out(:,:,:) !< output field
    LOGICAL, INTENT(IN), OPTIONAL   :: opt_lvalue_add !< If .TRUE., add values to given field
    ! local constants:
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':read_grib_3d'
    ! local variables:
    INTEGER               :: vlistID, varID, zaxisID, gridID,   &
      &                      mpi_comm, j, jl, jb, jk, ierrstat, &
      &                      dimlen(3), nmiss
    REAL(wp), ALLOCATABLE :: tmp_buf(:) ! temporary local array
    LOGICAL               :: lvalue_add


    lvalue_add = .FALSE.
    CALL assign_if_present(lvalue_add, opt_lvalue_add)

    ! allocate a buffer for one vertical level
    ALLOCATE(tmp_buf(glb_arr_len), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! get var ID
    IF(p_pe == p_io) THEN
      vlistID   = streamInqVlist(streamID)
      varID     = get_varID(streamID, name=TRIM(varname), opt_tileidx=opt_tileidx)
      zaxisID   = vlistInqVarZaxis(vlistID, varID)
      gridID    = vlistInqVarGrid(vlistID, varID)
      dimlen(1) = gridInqSize(gridID)
      dimlen(2) = zaxisInqSize(zaxisID)

!DR write(0,*) "varname, varID, zaxisID, dimlen(2), nlevs: ", TRIM(varname), &
!DR  & varID, zaxisID, dimlen(2), nlevs
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

  END SUBROUTINE read_grib_3d


  !-------------------------------------------------------------------------
  !> Read 2D dataset from GRIB2 file
  !
  !  @par Revision History
  ! 
  !  Initial revision by F. Prill, DWD (2013-02-19)
  !
  SUBROUTINE read_grib_2d (streamID, varname, glb_arr_len, loc_arr_len, glb_index, var_out, opt_tileidx)

    INTEGER,          INTENT(IN)    :: streamID       !< ID of CDI file stream
    CHARACTER(len=*), INTENT(IN)    :: varname        !< Var name of field to be read
    INTEGER,          INTENT(IN)    :: glb_arr_len    !< length of 1D field (global)
    INTEGER,          INTENT(IN)    :: loc_arr_len    !< length of 1D field (local)
    INTEGER,          INTENT(IN)    :: glb_index(:)   !< Index mapping local to global
    REAL(wp),         INTENT(INOUT) :: var_out(:,:)   !< output field
    INTEGER,          INTENT(IN), OPTIONAL :: opt_tileidx  !< tile index, encoded as "localInformationNumber"
    ! local variables:
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':read_grib_2d'
    INTEGER       :: varID, mpi_comm, j, jl, jb, &
      &              nmiss, vlistID, gridID
    REAL(wp)      :: z_dummy_array(glb_arr_len)       !< local dummy array


    ! Get var ID
    IF (p_pe == p_io) THEN
      vlistID   = streamInqVlist(streamID)
      varID     = get_varID(streamID, name=TRIM(varname), opt_tileidx=opt_tileidx)
      gridID    = vlistInqVarGrid(vlistID, varID)
      ! Check variable dimensions:
      IF (gridInqSize(gridID) /= glb_arr_len) THEN
        CALL finish(routine, "Incompatible dimensions!")
      END IF
    END IF

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
  END SUBROUTINE read_grib_2d

END MODULE mo_util_grib
