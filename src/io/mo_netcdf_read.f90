!>
!! This module provides basic methods for reading 
!! a NetCDF file in parallel or sequential in a transparent way.
!!
!! Contains routines for reading data from netcdf-Files of various shape.
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!!
!! @author Daniel Reinert, DWD
!! @author Leonidas Linardakis, MPIM
!!
!!
!! @par Revision History
!! Moved here from mo_ext_data by Daniel reinert, DWD (2012-02-01)
!! Moved from  mo_util_netcdf, added netcdf_read_oncells_2D, by L. Linardakis, (2013-03-15)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
!!
MODULE mo_netcdf_read

  USE mo_kind
  USE mo_mpi
  USE mo_gather_scatter,     ONLY: scatter_cells
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message_text, message, warning, finish, em_warn
  USE mo_impl_constants,     ONLY: success, max_char_length
  USE mo_parallel_config,    ONLY: nproma
  USE mo_io_units,           ONLY: filename_max

  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work
  USE mo_fortran_tools,      ONLY: assign_if_present
  !-------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'netcdf.inc'

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: read_netcdf_data
  PUBLIC :: read_netcdf_data_single
  PUBLIC :: read_netcdf_lu
  PUBLIC :: nf
  PUBLIC :: netcdf_open_input, netcdf_close
  PUBLIC :: netcdf_read_ONCELLS_2D

  INTERFACE read_netcdf_data
    MODULE PROCEDURE read_netcdf_2d
    MODULE PROCEDURE read_netcdf_2d_int
    MODULE PROCEDURE read_netcdf_3d
    MODULE PROCEDURE read_netcdf_aero
    MODULE PROCEDURE read_netcdf_4d
    MODULE PROCEDURE read_netcdf_time
  END INTERFACE read_netcdf_data

  INTERFACE read_netcdf_data_single
    MODULE PROCEDURE read_netcdf_3d_single
  END INTERFACE read_netcdf_data_single

  INTERFACE netcdf_read_oncells_2D
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_2D_filename
    MODULE PROCEDURE netcdf_read_REAL_ONCELLS_2D_fileid
  END INTERFACE netcdf_read_oncells_2D

  INTEGER, PARAMETER :: MAX_VAR_DIMS = NF_MAX_VAR_DIMS
  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_read_REAL_ONCELLS_2D_filename(filename, variable_name, fill_array, patch)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp), POINTER            :: fill_array(:,:)
    TYPE(t_patch)                :: patch

    INTEGER :: ncid
    INTEGER :: return_status    
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_2D_filename'

    ncid = netcdf_open_input(filename)
    netcdf_read_REAL_ONCELLS_2D_filename = &
      & netcdf_read_REAL_ONCELLS_2D_fileid(ncid, variable_name, fill_array, patch)
    return_status = netcdf_close(ncid)
                              
  END FUNCTION netcdf_read_REAL_ONCELLS_2D_filename
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_read_REAL_ONCELLS_2D_fileid(ncid, variable_name, fill_array, patch)
    INTEGER, INTENT(IN)          :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp), POINTER            :: fill_array(:,:)
    TYPE(t_patch)                :: patch

    INTEGER :: total_number_of_cells
    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    INTEGER :: return_status
    REAL(wp), POINTER :: tmp_array(:)
    
    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_netcdf_read:netcdf_read_REAL_ONCELLS_2D_fileid'


    ! LEONIDAS, CAN YOU PLEASE HAVE A LOOK

    ! trivial return value.
    netcdf_read_REAL_ONCELLS_2D_fileid = 0


    total_number_of_cells = patch%n_patch_cells_g
    
    IF( my_process_is_mpi_workroot()  ) THEN
      CALL nf(netcdf_inq_var(ncid, variable_name, varid, var_type, var_dims, var_size), variable_name)
      
      IF (var_dims /= 1 .OR. var_size(1) /= total_number_of_cells) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=", var_size, &
          & " total_number_of_cells=", total_number_of_cells
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF
    
    ALLOCATE( tmp_array(total_number_of_cells), stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF
    
    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf_get_var_double(ncid, varid, tmp_array(:)), variable_name)
    ENDIF

    IF (.NOT. ASSOCIATED(fill_array)) THEN
      ALLOCATE( fill_array(nproma, patch%nblks_c), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( fill_array )')
      ENDIF
    ENDIF
    
    CALL scatter_cells(tmp_array, fill_array, patch)

    DEALLOCATE(tmp_array)    
                              
  END FUNCTION netcdf_read_REAL_ONCELLS_2D_fileid
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_open_input(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: ncid

    IF( my_process_is_mpi_workroot()  ) THEN
        CALL nf(nf_open(TRIM(filename), nf_nowrite, ncid), TRIM(filename))
    ELSE
        ncid = -1 ! set it to an invalid value
    ENDIF

    netcdf_open_input = ncid
    
  END FUNCTION netcdf_open_input
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_close(ncid)
    INTEGER, INTENT(IN) :: ncid

    netcdf_close = -1
    IF( my_process_is_mpi_workroot()  ) THEN
        netcdf_close = nf_close(ncid)
    ENDIF

  END FUNCTION netcdf_close
  !-------------------------------------------------------------------------
  

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_inq_var(ncid, name, varid, var_type, var_dims, var_size)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: name
    
    INTEGER, INTENT(OUT) :: varid, var_type, var_dims
    INTEGER, INTENT(OUT) :: var_size(MAX_VAR_DIMS)

    INTEGER  :: number_of_attributes  
    INTEGER :: var_dims_reference(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: check_var_name
    INTEGER :: i, return_status

    netcdf_inq_var = -1
    IF ( .NOT. my_process_is_mpi_workroot() ) RETURN


    netcdf_inq_var = nf_inq_varid(ncid, name, varid)
    CALL nf(netcdf_inq_var, name)
    netcdf_inq_var = nf_inq_var (ncid, varid, check_var_name, var_type, var_dims, &
      & var_dims_reference, number_of_attributes)
    CALL nf(netcdf_inq_var, check_var_name)
    DO i=1, var_dims
!       return_status = nf_inq_dimlen(ncid, var_dims_reference(i), var_size(i))
      CALL nf(nf_inq_dimlen(ncid, var_dims_reference(i), var_size(i)), check_var_name)
    ENDDO
!     write(0,*) " Read var_dims, var_size:",  var_dims, var_size
!     write(0,*) " check_var_name:",  check_var_name
!     write(0,*) " name:", name    

  END FUNCTION netcdf_inq_var
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Read dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !!
  SUBROUTINE read_netcdf_2d (ncid, varname, glb_arr_len, loc_arr_len, glb_index, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_2d'

    INTEGER :: varid, mpi_comm, j, jl, jb
    REAL(wp):: z_dummy_array(glb_arr_len)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF( my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:)), routine)
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:) = 0._wp

    ! Set var_out from global data
    DO j = 1, loc_arr_len

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      var_out(jl,jb) = z_dummy_array(glb_index(j))
    ENDDO

  END SUBROUTINE read_netcdf_2d



  !-------------------------------------------------------------------------
  !>
  !! Read dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !!
  SUBROUTINE read_netcdf_2d_int (ncid, varname, glb_arr_len, loc_arr_len, glb_index, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    INTEGER, INTENT(INOUT) :: &          !< output field
      &  var_out(:,:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_2d_int'

    INTEGER :: varid, mpi_comm, j, jl, jb
    INTEGER :: z_dummy_array(glb_arr_len)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF( my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF( my_process_is_stdio()) CALL nf(nf_get_var_int(ncid, varid, z_dummy_array(:)), routine)
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:) = 0

    ! Set var_out from global data
    DO j = 1, loc_arr_len

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      var_out(jl,jb) = z_dummy_array(glb_index(j))
    ENDDO

  END SUBROUTINE read_netcdf_2d_int



  !-------------------------------------------------------------------------
  !>
  !! Read 3D (inlcuding height) dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !! Adapted for 3 D by Guenther Zaengl, DWD (2011-07-11)
  !!
  SUBROUTINE read_netcdf_3d (ncid, varname, glb_arr_len, loc_arr_len, glb_index, &
    &                        nlevs, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: nlevs         !< vertical levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_3d'

    INTEGER :: varid, mpi_comm, j, jl, jb, jk
    REAL(wp):: z_dummy_array(glb_arr_len,nlevs)!< local dummy array

    !-------------------------------------------------------------------------

    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:,:)), routine)
    CALL p_bcast(z_dummy_array, p_io, mpi_comm)

    var_out(:,:,:) = 0._wp

    ! Set var_out from global data
     DO jk = 1, nlevs
       DO j = 1, loc_arr_len

         jb = blk_no(j) ! Block index in distributed patch
         jl = idx_no(j) ! Line  index in distributed patch

         var_out(jl,jk,jb) = z_dummy_array(glb_index(j),jk)

       ENDDO
     ENDDO

  END SUBROUTINE read_netcdf_3d


 !-------------------------------------------------------------------------
  !>
  !! Read 3D dataset from netcdf file in SINGLE PRECISION
  !!
  !! @par Revision History
  !! Initial revision by F. Prill, DWD (2012-02-15)
  !! Optional switch to read 3D field in 2D slices: F. Prill, DWD (2012-12-19)
  !!
  SUBROUTINE read_netcdf_3d_single (ncid, varname, glb_arr_len, loc_arr_len, glb_index, &
    &                               nlevs, var_out, opt_lvalue_add)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: nlevs         !< vertical levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global
    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_lvalue_add !< If .TRUE., add values to given field

    ! local constants:
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_3d_single'
    ! enable this flag to use a 3d buffer (which may be faster)
    LOGICAL, PARAMETER :: luse3dbuffer = .FALSE.
    ! time level (fixed)
    INTEGER, PARAMETER :: itime = 1

    ! local variables:
    INTEGER :: varid, mpi_comm, j, jl, jb, jk, &
      &        istart(3), icount(3), ierrstat, &
      &        dimlen(3), dims(3)
    ! SINGLE PRECISION local array
    REAL(sp), ALLOCATABLE:: tmp_buf(:,:)
    LOGICAL :: lvalue_add

    !-------------------------------------------------------------------------

    lvalue_add = .FALSE.
    CALL assign_if_present(lvalue_add, opt_lvalue_add)

    ! allocate temporary buffer:
    IF (luse3dbuffer) THEN
      ! allocate a buffer for all levels
      ALLOCATE(tmp_buf(glb_arr_len,nlevs), STAT=ierrstat)
    ELSE
      ! allocate a buffer for one vertical level
      ALLOCATE(tmp_buf(glb_arr_len,1), STAT=ierrstat)
    END IF
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! Get var ID
    IF(my_process_is_stdio()) THEN
      CALL nf(nf_inq_varid(ncid, TRIM(varname), varid), routine)

      ! Check variable dimensions:
      CALL nf(NF_INQ_VARDIMID(ncid, varid, dims(:)), routine)
      DO j=1,3
        CALL nf(NF_INQ_DIMLEN  (ncid, dims(j), dimlen(j)), routine)
      END DO
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
    IF (luse3dbuffer) THEN
      !-- 3D buffer implementation

      IF(my_process_is_stdio()) THEN
        CALL nf(nf_get_var_real(ncid, varid, tmp_buf(:,:)), routine)
      END IF
      ! broadcast data:
      CALL p_bcast(tmp_buf, p_io, mpi_comm)
      ! Set var_out from global data
      IF (lvalue_add) THEN
        DO jk = 1, nlevs
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jk,jb) = var_out(jl,jk,jb) + REAL(tmp_buf(glb_index(j),jk), wp)
          ENDDO
        ENDDO
      ELSE
        DO jk = 1, nlevs
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jk,jb) = REAL(tmp_buf(glb_index(j),jk), wp)
          ENDDO
        ENDDO
      END IF

    ELSE
      !-- 2D buffer implementation

      icount = (/ glb_arr_len,1,1 /)
      DO jk=1,nlevs
        istart = (/ 1,jk,itime /)
        IF(my_process_is_stdio()) THEN
          CALL nf(nf_get_vara_real(ncid, varid, &
            &     istart, icount, tmp_buf(:,:)), routine)
        END IF

        ! broadcast data:
        CALL p_bcast(tmp_buf, p_io, mpi_comm)
        ! Set var_out from global data
        IF (lvalue_add) THEN
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jk,jb) = var_out(jl,jk,jb) + REAL(tmp_buf(glb_index(j),1), wp)
          ENDDO
        ELSE
          DO j = 1, loc_arr_len
            jb = blk_no(j) ! Block index in distributed patch
            jl = idx_no(j) ! Line  index in distributed patch
            var_out(jl,jk,jb) = REAL(tmp_buf(glb_index(j),1), wp)
          ENDDO
        END IF
      END DO ! jk=1,nlevs

    END IF

    ! clean up
    DEALLOCATE(tmp_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE read_netcdf_3d_single


  !-------------------------------------------------------------------------

  SUBROUTINE read_netcdf_aero (ncid, ntime, varname, glb_arr_len, &
       &                     loc_arr_len, glb_index, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: ntime         !< time levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_aero'

    INTEGER :: varid, mpi_comm, j, jl, jb, jt
    REAL(wp):: z_dummy_array(glb_arr_len,ntime)!< local dummy array
  !-------------------------------------------------------------------------

    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF


    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:,:)), routine)
    CALL p_bcast(z_dummy_array, p_io , mpi_comm)

    var_out(:,:,:) = 0._wp

    ! Set var_out from global data
    DO jt = 1, ntime
           DO j = 1, loc_arr_len

             jb = blk_no(j) ! Block index in distributed patch
             jl = idx_no(j) ! Line  index in distributed patch

             var_out(jl,jb,jt) = z_dummy_array(glb_index(j),jt)

          ENDDO
    ENDDO

  END SUBROUTINE read_netcdf_aero




  !-------------------------------------------------------------------------
  ! Specific read-routine for LU_CLASS_FRACTION. Probably, read_netcdf_aero
  ! and read_netcdf_lu can be merged into a single routine in the near
  ! future.

  SUBROUTINE read_netcdf_lu (ncid, varname, glb_arr_len,            &
       &                     loc_arr_len, glb_index, nslice, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: nslice        !< slices o netcdf field
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_lu'

    INTEGER :: varid, mpi_comm, j, jl, jb, js
    REAL(wp):: z_dummy_array(glb_arr_len, nslice)!< local dummy array

  !-------------------------------------------------------------------------


    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF



    ! I/O PE reads and broadcasts data

    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:,:)), routine)
    CALL p_bcast(z_dummy_array, p_io , mpi_comm)

    var_out(:,:,:) = 0._wp


    ! Set var_out from global data
    DO js = 1, nslice
      DO j = 1, loc_arr_len

        jb = blk_no(j) ! Block index in distributed patch
        jl = idx_no(j) ! Line  index in distributed patch

        var_out(jl,jb,js) = z_dummy_array(glb_index(j),js)

      ENDDO
    ENDDO

  END SUBROUTINE read_netcdf_lu


  !-------------------------------------------------------------------------
  !>
  !! Read 4D (inlcuding height and time) dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !! Adapted for 4 D by Kristina Froehlich, MPI-M (2011-06-16)
  !!
  SUBROUTINE read_netcdf_4d (ncid, varname, glb_arr_len, &
       &                     loc_arr_len, glb_index, &
       &                     nlevs, ntime,      var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: nlevs         !< vertical levels of netcdf file
    INTEGER, INTENT(IN) :: ntime         !< time levels of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:,:)                !< dimensions: nproma, nlevs, nblks, ntime

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_4d'

    INTEGER :: varid, mpi_comm, j, jl, jb, jk, jt
    REAL(wp):: z_dummy_array(glb_arr_len,nlevs,ntime)!< local dummy array

    !-------------------------------------------------------------------------

    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data

    write(0,*) ' ncep set ',varname,': begin of read - whole time array'
    IF(my_process_is_stdio()) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:,:,:)), routine)
    CALL p_bcast(z_dummy_array, p_io , mpi_comm)

    var_out(:,:,:,:) = 0._wp

    ! Set var_out from global data
    DO jt = 1, ntime
        DO jk = 1, nlevs
           DO j = 1, loc_arr_len

             jb = blk_no(j) ! Block index in distributed patch
             jl = idx_no(j) ! Line  index in distributed patch

             var_out(jl,jk,jb,jt) = z_dummy_array(glb_index(j),jk,jt)

          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE read_netcdf_4d


  !-------------------------------------------------------------------------
  !>
  !! Read 3D (inlcuding a period of time) dataset from netcdf file
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !! Adapted for parallel runs by Rainer Johanni (2010-12-07)
  !! Adapted for time periods by Stephan Lorenz, MPI-M (2012-02-22)
  !!
  SUBROUTINE read_netcdf_time (ncid, varname, glb_arr_len, &
       &                       loc_arr_len, glb_index,     &
       &                       ntime, nstart, ncount, var_out)

    CHARACTER(len=*), INTENT(IN)  ::  &  !< Var name of field to be read
      &  varname

    INTEGER, INTENT(IN) :: ncid          !< id of netcdf file
    INTEGER, INTENT(IN) :: glb_arr_len   !< length of 1D field (global)
    INTEGER, INTENT(IN) :: loc_arr_len   !< length of 1D field (local)
    INTEGER, INTENT(IN) :: glb_index(:)  !< Index mapping local to global
    INTEGER, INTENT(IN) :: ntime         !< number of time steps to read
    INTEGER, INTENT(IN) :: nstart(2)     !< start value for reading in all array dims
    INTEGER, INTENT(IN) :: ncount(2)     !< count value for length of array to read

    REAL(wp), INTENT(INOUT) :: &         !< output field
      &  var_out(:,:,:)                  !< dimensions: nproma, nblks, ntime

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_util_netcdf:read_netcdf_time'

    INTEGER :: varid, mpi_comm, j, jl, jb, jt
    REAL(wp):: z_dummy_array(glb_arr_len,ntime)!< local dummy array

    !-------------------------------------------------------------------------

    ! Get var ID
    IF(my_process_is_stdio()) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid), routine)

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    ! I/O PE reads and broadcasts data
    z_dummy_array(:,:) = 0.0_wp

    !write(0,*) ' Dimensions: glb, ntime: ',glb_arr_len,ntime
    !write(0,*) ' ncep set ',varname,': begin of read - time period'
    !write(0,*) ' nstart, ncount: ',nstart, ncount

    IF(my_process_is_stdio()) CALL nf(nf_get_vara_double(ncid, varid, &
      &                               nstart(:), ncount(:), z_dummy_array(:,:)), routine)
    CALL p_bcast(z_dummy_array, p_io , mpi_comm)

    var_out(:,:,:) = 0.0_wp

    ! Set var_out from global data
    DO jt = 1, ntime
      DO j = 1, loc_arr_len

        jb = blk_no(j) ! Block index in distributed patch
        jl = idx_no(j) ! Line  index in distributed patch

        var_out(jl,jb,jt) = z_dummy_array(glb_index(j),jt)

      ENDDO
    ENDDO

    !write(0,*) ' READ_NETCD_TIME: z_dummy_array stress-x, index 4*64+1,5:'
    !do jt=1,3
    !  write(0,*) 'jt=',jt,' val:',(z_dummy_array(j+4*64,jt),j=1,5)
    !enddo

  END SUBROUTINE read_netcdf_time


  !-------------------------------------------------------------------------


  SUBROUTINE nf(STATUS, routine, warnonly, silent)

    INTEGER, INTENT(in)           :: STATUS
    CHARACTER(len=*), INTENT(in) :: routine
    LOGICAL, INTENT(in), OPTIONAL :: warnonly
    LOGICAL, INTENT(in), OPTIONAL :: silent

    LOGICAL :: lwarnonly, lsilent

    lwarnonly = .FALSE.
    lsilent   = .FALSE.
    IF(PRESENT(warnonly)) lwarnonly = .TRUE.
    IF(PRESENT(silent))   lsilent   = silent

    IF (lsilent) RETURN


    IF (STATUS /= nf_noerr) THEN
      IF (lwarnonly) THEN
        CALL message( TRIM(routine)//' netCDF error', nf_strerror(STATUS), &
          & level=em_warn)
      ELSE
        CALL finish( TRIM(routine)//' netCDF error', nf_strerror(STATUS))
      ENDIF
    ENDIF

  END SUBROUTINE nf


END MODULE mo_netcdf_read
