!>
!!
!! This module provides methods for reading a NetCDF file in a distributed way.
!! This approach reduces memory consumption
!!
!! @par Revision History
!! Initial version by Moritz Hanke, December 2013
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_read_netcdf_distributed

  USE mo_kind, ONLY: wp, sp
  USE mo_exception, ONLY: finish, message, em_warn
  USE mo_mpi, ONLY: p_n_work, p_pe_work, p_bcast, p_comm_work, p_max
  USE ppm_extents, ONLY: extent
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info, &
    & t_glb2loc_index_lookup, init_glb2loc_index_lookup, &
    & set_inner_glb_index, deallocate_glb2loc_index_lookup, &
    & uniform_partition, partidx_of_elem_uniform_deco
  USE mo_communication, ONLY: t_comm_pattern, idx_no, blk_no, &
    & delete_comm_pattern, exchange_data
  USE mo_parallel_config, ONLY: nproma, io_process_stride, io_process_rotate
  USE mo_communication_factory, ONLY: setup_comm_pattern
  USE mo_fortran_tools, ONLY: t_ptr_2d, t_ptr_2d_int, t_ptr_2d_sp, &
    & t_ptr_3d, t_ptr_3d_int, t_ptr_3d_sp, t_ptr_4d, t_ptr_4d_int, &
    & t_ptr_4d_sp
#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
  USE mpi, ONLY: MPI_INFO_NULL, MPI_UNDEFINED, MPI_Comm_split, MPI_COMM_NULL
#endif

  IMPLICIT NONE
  PRIVATE

  CHARACTER(*), PARAMETER :: modname = 'mo_read_netcdf_distributed'

  PUBLIC :: distrib_read
  PUBLIC :: distrib_nf_open
  PUBLIC :: distrib_nf_inq_varexists
  PUBLIC :: distrib_nf_close
  PUBLIC :: setup_distrib_read
  PUBLIC :: delete_distrib_read
  PUBLIC :: t_distrib_read_data
  PUBLIC :: distrib_inq_var_dims
  PUBLIC :: idx_lvl_blk, idx_blk_time
  PUBLIC :: nf

  INCLUDE 'netcdf.inc'

  INTEGER, PARAMETER :: idx_lvl_blk = 2
  INTEGER, PARAMETER :: idx_blk_time = 3

  !modules interface-------------------------------------------
  !subroutines

  TYPE t_basic_distrib_read_data
    INTEGER :: n_g = -1 ! global number of points (-1 if unused)
    TYPE(extent) :: io_chunk ! io decomposition
    TYPE(t_glb2loc_index_lookup), POINTER :: glb2loc_index => NULL()
    INTEGER :: n_ref = 0 ! number of times this data is referenced
  END TYPE t_basic_distrib_read_data

  TYPE t_distrib_read_data
    INTEGER :: basic_data_index = -1
    CLASS(t_comm_pattern), POINTER :: pat => NULL()
  END TYPE t_distrib_read_data

  TYPE(t_basic_distrib_read_data), TARGET, ALLOCATABLE :: basic_data(:)

  ! This one only depends on n_io_processes and io_process_stride, so it is
  ! save to store it in a module variable.
  INTEGER :: parRootRank = -1
#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
  INTEGER :: io_comm = MPI_COMM_NULL
#endif
  LOGICAL :: this_PE_does_IO = .FALSE.

CONTAINS

  INTEGER FUNCTION distrib_nf_open(path) RESULT(nfid)
    CHARACTER(*), INTENT(in) :: path
#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
    INTEGER              :: ierr, nvars, i
    INTEGER, ALLOCATABLE :: varids(:)
    LOGICAL              :: exists
#endif
    CHARACTER(*), PARAMETER :: routine = modname//'::distrib_nf_open'

    CALL message (routine, path)
    nfid = -1
    IF (this_PE_does_IO) THEN
#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
      ierr = nf_open_par(path, IOR(nf_nowrite, nf_mpiio), io_comm, &
        & MPI_INFO_NULL, nfid)

      ! We do our own error handling here to give the filename to the user if
      ! a file does not exist.
      IF (ierr == nf_noerr) THEN
        ! Switch all vars to collective. Hopefully this is sufficient.
        CALL nf(nf_inq_nvars(nfid, nvars))
        ALLOCATE(varids(nvars))
        CALL nf(nf_inq_varids(nfid, nvars, varids))
        DO i = 1,nvars
          CALL nf(nf_var_par_access(nfid, varids(i), NF_COLLECTIVE))
        ENDDO
      ELSE
        INQUIRE(file=path, exist=exists)
        IF (.NOT. exists) THEN
          CALL finish("mo_read_netcdf_distributed", "File "//TRIM(path)//" does not exist.")
        ELSE
          ! If file exists just do the usual thing.
          CALL nf(nf_open(path, nf_nowrite, nfid))
          CALL message(routine, 'warning: falling back to serial semantics for&
               & opening netcdf file '//path)
        ENDIF
      ENDIF
#else
      CALL nf(nf_open(path, nf_nowrite, nfid))
#endif
    END IF
  END FUNCTION distrib_nf_open

  !-------------------------------------------------------------------------

  LOGICAL FUNCTION distrib_nf_inq_varexists(ncid, vname) result(ret)
    INTEGER, INTENT(in) :: ncid
    CHARACTER(*), INTENT(in) :: vname
    INTEGER :: err, vid

    IF (p_pe_work == parRootRank) err = nf_inq_varid(ncid, vname, vid)
    CALL p_bcast(err, parRootRank, p_comm_work)
    ret = (err == nf_noerr)
  END FUNCTION distrib_nf_inq_varexists

  SUBROUTINE distrib_inq_var_dims(file_id, var_name, var_ndims, var_dimlen)
    INTEGER, INTENT(IN) :: file_id
    CHARACTER(*), INTENT(IN) :: var_name
    INTEGER, INTENT(OUT) :: var_ndims, var_dimlen(:)
    INTEGER :: varid, i
    INTEGER :: temp_var_dimlen(NF_MAX_VAR_DIMS), var_dimids(NF_MAX_VAR_DIMS)

    IF ( p_pe_work .EQ. parRootRank ) THEN
      CALL nf(nf_inq_varid(file_id, var_name, varid))
      CALL nf(nf_inq_varndims(file_id, varid, var_ndims))
      CALL nf(nf_inq_vardimid(file_id, varid, var_dimids))
      DO i=1, var_ndims
        CALL nf(nf_inq_dimlen(file_id, var_dimids(i), temp_var_dimlen(i)))
      ENDDO
    END IF
#ifndef NOMPI
    CALL p_bcast(var_ndims, parRootRank, p_comm_work)
    CALL p_bcast(temp_var_dimlen(1:var_ndims), parRootRank, p_comm_work)
#endif
    IF (SIZE(var_dimlen) .LT. var_ndims) &
      CALL finish("distrib_inq_var_dims", &
        &         "array size of argument var_dimlen is too small")
    var_dimlen(1:var_ndims) = temp_var_dimlen(1:var_ndims)
  END SUBROUTINE distrib_inq_var_dims

  !-------------------------------------------------------------------------

  SUBROUTINE distrib_nf_close(ncid)
    INTEGER, INTENT(in) :: ncid

    IF (this_PE_does_IO) CALL nf(nf_close(ncid))
  END SUBROUTINE distrib_nf_close

  !-------------------------------------------------------------------------

  SUBROUTINE setup_distrib_read(n_g, decomp_info, io_data)
    INTEGER, INTENT(in) :: n_g
    TYPE(t_grid_domain_decomp_info), INTENT(in) :: decomp_info
    TYPE(t_distrib_read_data), INTENT(inout) :: io_data
    INTEGER :: n, n_inner, i
    INTEGER, SAVE :: io_stride = -1, n_io_proc = -1
    INTEGER, ALLOCATABLE :: owner(:)
    TYPE(t_basic_distrib_read_data), POINTER :: bio

    IF (n_io_proc .EQ. -1) CALL init_module_vars()
    io_data%basic_data_index = get_data_idx()
    bio => basic_data(io_data%basic_data_index)
    n = 0
    IF (ALLOCATED(decomp_info%glb_index)) n = SIZE(decomp_info%glb_index)
    ALLOCATE(owner(MAX(1,n)))
    IF (n .GT. 0) &
      owner = (partidx_of_elem_uniform_deco(extent(1, bio%n_g), & 
         n_io_proc,  decomp_info%glb_index(:)) - 1) * io_stride &
         + MODULO(io_process_rotate, io_stride)
    n_inner = 0
    IF (ALLOCATED(bio%glb2loc_index%inner_glb_index)) &
      & n_inner = SIZE(bio%glb2loc_index%inner_glb_index, 1)
    CALL setup_comm_pattern(n, owner(:), decomp_info%glb_index(:), &
      & bio%glb2loc_index, n_inner, &
      & (/(p_pe_work, i = 1, n_inner)/), &
      & bio%glb2loc_index%inner_glb_index, io_data%pat)
  CONTAINS

    SUBROUTINE init_module_vars()
      INTEGER :: rotate, temp_n, temp_stride, ierr, myColor

      rotate = MODULO(io_process_rotate, io_stride)
      IF (io_process_stride .GT. 0) THEN
        io_stride = MAX(1, MODULO(io_process_stride, p_n_work))
        n_io_proc = (p_n_work - rotate + io_stride - 1) / io_stride
      ELSE
        temp_n = NINT(SQRT(REAL(p_n_work)))
        temp_stride = (p_n_work + temp_n - 1) / temp_n
        ! improve io process stride by rounding to the next power of two
        io_stride = 2**CEILING(LOG(REAL(temp_stride))/LOG(2.))
        n_io_proc = (p_n_work + io_stride - 1) / io_stride
      END IF
      this_PE_does_IO = MOD(p_pe_work, io_stride) == rotate
      parRootRank = MERGE(p_pe_work, 0, this_PE_does_IO)
      parRootRank = p_max(parRootRank, p_comm_work)
#if defined (HAVE_PARALLEL_NETCDF) && !defined (NOMPI)
      myColor = MERGE(1, MPI_UNDEFINED, this_PE_does_IO)
      IF (io_comm .EQ. MPI_COMM_NULL) &
        CALL MPI_Comm_split(p_comm_work, myColor, 0, io_comm, ierr)
#endif
    END SUBROUTINE init_module_vars

    INTEGER FUNCTION get_data_idx() RESULT(idx)
      INTEGER :: j, m, j_empty
      TYPE(t_basic_distrib_read_data), ALLOCATABLE :: temp_bio(:)

      m = 0
      j_empty = 0
      idx = 0
      IF (ALLOCATED(basic_data)) m = SIZE(basic_data)
      DO j = 1, m
        idx = MERGE(j, 0, basic_data(j)%n_g .EQ. n_g)
        IF (j_empty .EQ. 0) &
          j_empty = MERGE(j, 0, basic_data(j)%n_g .EQ. -1)
        IF (idx .NE. 0) EXIT
      END DO
      IF (idx .EQ. 0) idx = j_empty
      IF (idx .EQ. 0) THEN
        IF (m .GT. 0) CALL MOVE_ALLOC(basic_data, temp_bio)
        ALLOCATE(basic_data(m + 16))
        IF (m .GT. 0) basic_data(1:m) = temp_bio(1:m)
        idx = m + 1
      END IF
      IF (basic_data(idx)%n_ref .EQ. 0) CALL init_data(idx)
      basic_data(idx)%n_ref = basic_data(idx)%n_ref + 1
    END FUNCTION get_data_idx

    SUBROUTINE init_data(idx)
      INTEGER, INTENT(IN) :: idx
      INTEGER :: j
      INTEGER, ALLOCATABLE :: idxmap(:,:)

      basic_data(idx)%n_g = n_g
      ALLOCATE(basic_data(idx)%glb2loc_index)
      CALL init_glb2loc_index_lookup(basic_data(idx)%glb2loc_index, n_g)
      ! if the process takes part in the reading
      IF (this_PE_does_IO) THEN
        basic_data(idx)%io_chunk = &
          & uniform_partition(extent(1, n_g), n_io_proc, p_pe_work/io_stride + 1)
        ALLOCATE(idxmap(basic_data(idx)%io_chunk%size, 2))
        FORALL(j = 1:basic_data(idx)%io_chunk%size)
          idxmap(j, 1) = j + basic_data(idx)%io_chunk%first - 1
          idxmap(j, 2) = j
        END FORALL
        CALL set_inner_glb_index(basic_data(idx)%glb2loc_index, &
             idxmap(:, 1), idxmap(:, 2))
      ELSE
        basic_data(idx)%io_chunk = extent(1, 0)
      END IF
    END SUBROUTINE init_data
  END SUBROUTINE setup_distrib_read

  !-------------------------------------------------------------------------

  SUBROUTINE delete_distrib_read(io_data)
    TYPE(t_distrib_read_data), INTENT(inout) :: io_data
    INTEGER :: idx
    CHARACTER(*), PARAMETER :: routine = modname//"::delete_distrib_read"

    CALL delete_comm_pattern(io_data%pat)
    idx = io_data%basic_data_index
    IF (.NOT.ALLOCATED(basic_data)) CALL finish(routine, "basic_data not allocated")
    IF (idx .GT. SIZE(basic_data)) CALL finish(routine, "invalid index")
    IF (basic_data(idx)%n_ref .LT. 1) CALL finish(routine, "invalid reference counter")
    basic_data(idx)%n_ref = basic_data(idx)%n_ref - 1
    ! some other distrib_read is still using this basic_distrib_read_data
    IF (basic_data(idx)%n_ref .GT. 0) RETURN
    basic_data(idx)%n_g = -1
    CALL deallocate_glb2loc_index_lookup(basic_data(idx)%glb2loc_index)
    DEALLOCATE(basic_data(idx)%glb2loc_index)
  END SUBROUTINE delete_distrib_read

  !-------------------------------------------------------------------------

  SUBROUTINE distrib_read(ncid, vname, vdata, iod, &
    & edim, dimo, start_ext_dim, end_ext_dim)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: vname
    CLASS(*), INTENT(INOUT) :: vdata(:)
    TYPE(t_distrib_read_data), INTENT(IN) :: iod(:)
    INTEGER, INTENT(IN), OPTIONAL :: edim(:), dimo, start_ext_dim(:), end_ext_dim(:)
    INTEGER, ALLOCATABLE :: bufi_i(:,:,:), bufo_i(:,:,:,:)
    REAL(wp), ALLOCATABLE :: bufi_d(:,:,:), bufo_d(:,:,:,:)
    REAL(sp), ALLOCATABLE :: bufi_s(:,:,:), bufo_s(:,:,:,:)
    INTEGER :: ish(3), osh(4), vid, vtype, strt(3)
    TYPE(t_basic_distrib_read_data), POINTER :: bio
    CHARACTER(*), PARAMETER :: routine = modname//"distrib_read"

    IF (SIZE(iod) == 0) RETURN
    IF (SIZE(vdata) < SIZE(iod)) CALL finish(routine, "var_data too small")
    IF (SIZE(iod) .GT. 1) THEN
      IF (.NOT.ALL(iod(2:)%basic_data_index .EQ. iod(1)%basic_data_index)) &
        & CALL finish(routine, "basic_data_index do not match")
    END IF
    bio => basic_data(iod(1)%basic_data_index)
    ish(:) = [bio%io_chunk%size, 1, 1]
    osh(:) = [nproma, (ish(1)+nproma-1)/nproma, 1, 1]
    strt(:) = [bio%io_chunk%first, 1, 1]
    IF (PRESENT(edim)) THEN
      IF (SIZE(edim) .EQ. 1) THEN
        IF (.NOT.PRESENT(dimo)) CALL finish(routine, "invalid arguments")
        SELECT CASE(dimo)
        CASE(idx_blk_time)
          osh(3) = edim(1)
        CASE(idx_lvl_blk)
          osh(3) = osh(2)
          osh(2) = edim(1)
        CASE DEFAULT
          CALL finish(routine, "invalid dim_order")
        END SELECT
        ish(2) = edim(1)
      ELSE IF (SIZE(edim) .EQ. 2) THEN
        osh(3) = osh(2)
        osh(2) = edim(1)
        osh(4) = edim(2)
        ish(2:3) = edim(1:2)
      END IF
      IF (PRESENT(start_ext_dim)) THEN
        IF (.NOT.PRESENT(end_ext_dim)) CALL finish(routine, "invalid arguments")
        IF (SIZE(edim) .NE. SIZE(start_ext_dim) .OR. SIZE(edim) .NE. SIZE(end_ext_dim)) &
          & CALL finish(routine, "invalid arguments")
        IF (ANY((end_ext_dim - start_ext_dim + 1) /= edim)) &
          & CALL finish(routine, "invalid arguments: edim /= end-start-1")
        strt(2:SIZE(edim)+1) = start_ext_dim(1:SIZE(edim))
      END IF
    END IF
    IF (ish(1) .GT. 0) CALL nf(nf_inq_varid(ncid, vname, vid))
    SELECT TYPE(vdata)
    TYPE IS(t_ptr_2d)
      CALL read_multi_var_2dwp(vdata)
    TYPE IS(t_ptr_2d_sp)
      CALL read_multi_var_2dsp(vdata)
    TYPE IS(t_ptr_2d_int)
      CALL read_multi_var_2dint(vdata)
    TYPE IS(t_ptr_3d)
      CALL read_multi_var_3dwp(vdata, dimo)
    TYPE IS(t_ptr_3d_sp)
      CALL read_multi_var_3dsp(vdata, dimo)
    TYPE IS(t_ptr_3d_int)
      CALL read_multi_var_3dint(vdata, dimo)
    TYPE IS(t_ptr_4d)
      CALL read_multi_var_4dwp(vdata)
    TYPE IS(t_ptr_4d_sp)
      CALL read_multi_var_4dsp(vdata)
    TYPE IS(t_ptr_4d_int)
      CALL read_multi_var_4dint(vdata)
    CLASS DEFAULT
      CALL finish("distrib_read_multi_var", "un-recognized type")
    END SELECT
  CONTAINS

    SUBROUTINE read_multi_var_2dint(vd)
      TYPE(t_ptr_2d_int), INTENT(INOUT) :: vd(:)
      INTEGER :: i

      ALLOCATE(bufi_i(ish(1),ish(2),ish(3)), bufo_i(osh(1),osh(2),osh(3),osh(4)))
      IF (ish(1) > 0) THEN
        CALL nf(nf_inq_vartype(ncid, vid, vtype))
        IF (vtype .NE. NF_INT) CALL finish(routine, "not an NF_INT")
        CALL nf(nf_get_vara_int(ncid, vid, strt(1:1), ish(1:1), bufi_i(:,1,1)))
        DO i = 1, ish(1)
          bufo_i(idx_no(i),blk_no(i),1,1) = bufi_i(i,1,1)
        END DO
      END IF
      DO i = 1, SIZE(iod)
        CALL exchange_data(iod(i)%pat, vd(i)%p, bufo_i(:,:,1,1))
      END DO
    END SUBROUTINE read_multi_var_2dint

    SUBROUTINE read_multi_var_2dwp(vd)
      TYPE(t_ptr_2d), INTENT(INOUT) :: vd(:)
      INTEGER :: i

      ALLOCATE(bufi_d(ish(1),ish(2),ish(3)), bufo_d(osh(1),osh(2),osh(3),osh(4)))
      IF (ish(1) > 0) THEN
        CALL nf(nf_get_vara_double(ncid, vid, strt(1:1), ish(1:1), bufi_d(:,1,1)))
        DO i = 1, ish(1)
          bufo_d(idx_no(i),blk_no(i),1,1) = bufi_d(i,1,1)
        END DO
      END IF
      DO i = 1, SIZE(iod)
        CALL exchange_data(iod(i)%pat, vd(i)%p, bufo_d(:,:,1,1))
      END DO
    END SUBROUTINE read_multi_var_2dwp

    SUBROUTINE read_multi_var_2dsp(vd)
      TYPE(t_ptr_2d_sp), INTENT(INOUT) :: vd(:)
      INTEGER :: i

      ALLOCATE(bufi_s(ish(1),ish(2),ish(3)), bufo_s(osh(1),osh(2),osh(3),osh(4)))
      IF (ish(1) > 0) THEN
        CALL nf(nf_get_vara_real(ncid, vid, strt(1:1), ish(1:1), bufi_s(:,1,1)))
        DO i = 1, ish(1)
          bufo_s(idx_no(i),blk_no(i),1,1) = bufi_s(i,1,1)
        END DO
      END IF
      DO i = 1, SIZE(iod)
        CALL exchange_data(iod(i)%pat, vd(i)%p, bufo_s(:,:,1,1))
      END DO
    END SUBROUTINE read_multi_var_2dsp

    SUBROUTINE read_multi_var_3dint(vd, o)
      TYPE(t_ptr_3d_int), INTENT(INOUT) :: vd(:)
      INTEGER, INTENT(IN) :: o
      INTEGER :: i, j

      ALLOCATE(bufi_i(ish(1),ish(2),ish(3)), bufo_i(osh(1),osh(2),osh(3),osh(4)))
      IF (ish(1) > 0) THEN
        CALL nf(nf_inq_vartype(ncid, vid, vtype))
        IF (vtype .NE. NF_INT) CALL finish(routine, "not an NF_INT")
        CALL nf(nf_get_vara_int(ncid, vid, strt(1:2), ish(1:2), bufi_i(:,:,1)))
        IF (o .EQ. idx_blk_time) THEN
          DO j = 1, ish(2)
            DO i = 1, ish(1)
              bufo_i(idx_no(i),blk_no(i),j,1) = bufi_i(i,j,1)
            END DO
          END DO
        ELSE IF(o .EQ. idx_lvl_blk) THEN
          DO j = 1, ish(2)
            DO i = 1, ish(1)
              bufo_i(idx_no(i),j,blk_no(i),1) = bufi_i(i,j,1)
            END DO
          END DO          
        END IF
      END IF
      DO i = 1, SIZE(iod)
        IF (o .EQ. idx_blk_time) THEN
          DO j = 1, ish(2)
            CALL exchange_data(iod(i)%pat, vd(i)%p(:,:,j), bufo_i(:,:,j,1))
          END DO
        ELSE IF(o .EQ. idx_lvl_blk) THEN
          CALL exchange_data(iod(i)%pat, vd(i)%p, bufo_i(:,:,:,1))
        END IF
      END DO
    END SUBROUTINE read_multi_var_3dint

    SUBROUTINE read_multi_var_3dwp(vd, o)
      TYPE(t_ptr_3d), INTENT(INOUT) :: vd(:)
      INTEGER, INTENT(IN) :: o
      INTEGER :: i, j

      ALLOCATE(bufi_d(ish(1),ish(2),ish(3)), bufo_d(osh(1),osh(2),osh(3),osh(4)))
      IF (ish(1) > 0) THEN
        CALL nf(nf_get_vara_double(ncid, vid, strt(1:2), ish(1:2), bufi_d(:,:,1)))
        IF (o .EQ. idx_blk_time) THEN
          DO j = 1, ish(2)
            DO i = 1, ish(1)
              bufo_d(idx_no(i),blk_no(i),j,1) = bufi_d(i,j,1)
            END DO
          END DO
        ELSE IF(o .EQ. idx_lvl_blk) THEN
          DO j = 1, ish(2)
            DO i = 1, ish(1)
              bufo_d(idx_no(i),j,blk_no(i),1) = bufi_d(i,j,1)
            END DO
          END DO
        END IF
      END IF
      DO i = 1, SIZE(iod)
        IF (o .EQ. idx_blk_time) THEN
          DO j = 1, ish(2)
            CALL exchange_data(iod(i)%pat, vd(i)%p(:,:,j), bufo_d(:,:,j,1))
          END DO
        ELSE IF(o .EQ. idx_lvl_blk) THEN
          CALL exchange_data(iod(i)%pat, vd(i)%p, bufo_d(:,:,:,1))
        END IF
      END DO
    END SUBROUTINE read_multi_var_3dwp

    SUBROUTINE read_multi_var_3dsp(vd, o)
      TYPE(t_ptr_3d_sp), INTENT(INOUT) :: vd(:)
      INTEGER, INTENT(IN) :: o
      INTEGER :: i, j

      ALLOCATE(bufi_s(ish(1),ish(2),ish(3)), bufo_s(osh(1),osh(2),osh(3),osh(4)))
      IF (ish(1) > 0) THEN
        CALL nf(nf_get_vara_real(ncid, vid, strt(1:2), ish(1:2), bufi_s(:,:,1)))
        IF (o .EQ. idx_blk_time) THEN
          DO j = 1, ish(2)
            DO i = 1, ish(1)
              bufo_s(idx_no(i),blk_no(i),j,1) = bufi_s(i,j,1)
            END DO
          END DO
        ELSE IF(o .EQ. idx_lvl_blk) THEN
          DO j = 1, ish(2)
            DO i = 1, ish(1)
              bufo_s(idx_no(i),j,blk_no(i),1) = bufi_s(i,j,1)
            END DO
          END DO
        END IF
      END IF
      DO i = 1, SIZE(iod)
        IF (o .EQ. idx_blk_time) THEN
          DO j = 1, ish(2)
            CALL exchange_data(iod(i)%pat, vd(i)%p(:,:,j), bufo_s(:,:,j,1))
          END DO
        ELSE IF(o .EQ. idx_lvl_blk) THEN
          CALL exchange_data(iod(i)%pat, vd(i)%p, bufo_s(:,:,:,1))
        END IF
      END DO
    END SUBROUTINE read_multi_var_3dsp

    SUBROUTINE read_multi_var_4dint(vd)
      TYPE(t_ptr_4d_int), INTENT(INOUT) :: vd(:)
      INTEGER :: i, j, k

      ALLOCATE(bufi_i(ish(1),ish(2),ish(3)), bufo_i(osh(1),osh(2),osh(3),osh(4)))
      IF (ish(1) > 0) THEN
        CALL nf(nf_inq_vartype(ncid, vid, vtype))
        IF (vtype .NE. NF_INT) CALL finish(routine, "not an NF_INT")
        CALL nf(nf_get_vara_int(ncid, vid, strt, ish, bufi_i(:,:,:)))
        DO k = 1, ish(3)
          DO j = 1, ish(2)
            DO i = 1, ish(1)
              bufo_i(idx_no(i),j,blk_no(i), k) = bufi_i(i, j, k)
            END DO
          END DO
        END DO
      END IF
      DO i = 1, SIZE(iod)
        DO j = 1, ish(3)
          CALL exchange_data(iod(i)%pat, vd(i)%p(:,:,:,j), bufo_i(:,:,:,j))
        END DO
      END DO
    END SUBROUTINE read_multi_var_4dint

    SUBROUTINE read_multi_var_4dwp(vd)
      TYPE(t_ptr_4d), INTENT(INOUT) :: vd(:)
      INTEGER :: i, j, k

      ALLOCATE(bufi_d(ish(1),ish(2),ish(3)), bufo_d(osh(1),osh(2),osh(3),osh(4)))
      IF (ish(1) > 0) THEN
        CALL nf(nf_get_vara_double(ncid, vid, strt, ish, bufi_d(:,:,:)))
        DO k = 1, ish(3)
          DO j = 1, ish(2)
            DO i = 1, ish(1)
              bufo_d(idx_no(i),j,blk_no(i), k) = bufi_d(i, j, k)
            END DO
          END DO
        END DO
      END IF
      DO i = 1, SIZE(iod)
        DO j = 1, ish(3)
          CALL exchange_data(iod(i)%pat, vd(i)%p(:,:,:,j), bufo_d(:,:,:,j))
        END DO
      END DO
    END SUBROUTINE read_multi_var_4dwp

    SUBROUTINE read_multi_var_4dsp(vd)
      TYPE(t_ptr_4d_sp), INTENT(INOUT) :: vd(:)
      INTEGER :: i, j, k

      ALLOCATE(bufi_s(ish(1),ish(2),ish(3)), bufo_s(osh(1),osh(2),osh(3),osh(4)))
      IF (ish(1) > 0) THEN
        CALL nf(nf_get_vara_real(ncid, vid, strt, ish, bufi_s(:,:,:)))
        DO k = 1, ish(3)
          DO j = 1, ish(2)
            DO i = 1, ish(1)
              bufo_s(idx_no(i),j,blk_no(i), k) = bufi_s(i, j, k)
            END DO
          END DO
        END DO
      END IF
      DO i = 1, SIZE(iod)
        DO j = 1, ish(3)
          CALL exchange_data(iod(i)%pat, vd(i)%p(:,:,:,j), bufo_s(:,:,:,j))
        END DO
      END DO
    END SUBROUTINE read_multi_var_4dsp
  END SUBROUTINE distrib_read

  !-------------------------------------------------------------------------

  SUBROUTINE nf(ierrstat, warnonly, silent)
    INTEGER, INTENT(in)           :: ierrstat
    LOGICAL, INTENT(in), OPTIONAL :: warnonly, silent

    IF(PRESENT(silent)) THEN
      IF (silent) RETURN
    END IF
    IF (ierrstat .NE. nf_noerr) THEN
      IF (PRESENT(warnonly)) THEN
        CALL message(modname//' netCDF error', &
          &          nf_strerror(ierrstat), level=em_warn)
      ELSE
        CALL finish(modname//' netCDF error', &
          &         nf_strerror(ierrstat))
      ENDIF
    ENDIF
  END SUBROUTINE nf

END MODULE mo_read_netcdf_distributed
