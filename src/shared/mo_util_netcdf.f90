!>
!! Contains routines for reading netcdf-Files of various shape.
!!
!! Contains routines for reading data from netcdf-Files of various shape.
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Moved here from mo_ext_data by Daniel reinert, DWD (2012-02-01)
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
MODULE mo_util_netcdf

  USE mo_kind,               ONLY: wp, sp
  USE mo_exception,          ONLY: message, finish, em_warn
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work
  USE mo_fortran_tools,      ONLY: assign_if_present

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'


  PRIVATE

  PUBLIC  :: read_netcdf_data
  PUBLIC  :: read_netcdf_data_single
  PUBLIC  :: read_netcdf_lu
  PUBLIC  :: nf


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


  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

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


END MODULE mo_util_netcdf
