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

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish, message, em_warn
  USE mo_mpi, ONLY: p_n_work, p_pe_work, p_bcast, p_comm_work
  USE ppm_extents, ONLY: extent
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info, &
    & t_glb2loc_index_lookup, &
    & init_glb2loc_index_lookup, &
    & set_inner_glb_index, &
    & deallocate_glb2loc_index_lookup, &
    & uniform_partition, partidx_of_elem_uniform_deco
  USE mo_communication, ONLY: t_comm_pattern, idx_no, blk_no, &
    & delete_comm_pattern, exchange_data
  USE mo_parallel_config, ONLY: nproma, &
       config_io_process_stride => io_process_stride, &
       config_io_process_rotate => io_process_rotate
  USE mo_communication_factory, ONLY: setup_comm_pattern

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: distrib_read
  PUBLIC :: distrib_nf_open
  PUBLIC :: distrib_nf_close
  PUBLIC :: setup_distrib_read
  PUBLIC :: delete_distrib_read
  PUBLIC :: t_distrib_read_data
  PUBLIC :: var_data_1d_int
  PUBLIC :: var_data_2d_int, var_data_2d_wp
  PUBLIC :: var_data_3d_int, var_data_3d_wp
  PUBLIC :: var_data_4d_int, var_data_4d_wp
  PUBLIC :: distrib_inq_var_dims
  PUBLIC :: idx_lvl_blk, idx_blk_time
  PUBLIC :: nf

  INCLUDE 'netcdf.inc'

  INTERFACE distrib_read
    MODULE PROCEDURE distrib_read_int_2d_multi_var
    MODULE PROCEDURE distrib_read_int_3d_multi_var
    MODULE PROCEDURE distrib_read_int_4d_multi_var
    MODULE PROCEDURE distrib_read_real_2d_multi_var
    MODULE PROCEDURE distrib_read_real_3d_multi_var
    MODULE PROCEDURE distrib_read_real_4d_multi_var
    MODULE PROCEDURE distrib_read_int_2d
    MODULE PROCEDURE distrib_read_int_3d
    MODULE PROCEDURE distrib_read_int_4d
    MODULE PROCEDURE distrib_read_real_2d
    MODULE PROCEDURE distrib_read_real_3d
    MODULE PROCEDURE distrib_read_real_4d
  END INTERFACE distrib_read

  INTEGER, PARAMETER :: idx_lvl_blk = 2
  INTEGER, PARAMETER :: idx_blk_time = 3

  !modules interface-------------------------------------------
  !subroutines

  TYPE t_basic_distrib_read_data
    INTEGER :: n_g ! global number of points (-1 if unused)
    TYPE(extent) :: io_chunk ! io decomposition
    TYPE(t_glb2loc_index_lookup), POINTER :: glb2loc_index
    INTEGER :: io_process_stride, num_io_processes

    INTEGER :: n_ref ! number of times this data is referenced
  END TYPE t_basic_distrib_read_data

  TYPE t_distrib_read_data

    INTEGER :: basic_data_index
    CLASS(t_comm_pattern), POINTER :: redistrib_pattern

  END TYPE t_distrib_read_data

  TYPE var_data_1d_int
    INTEGER, POINTER :: DATA(:)
  END TYPE
  TYPE var_data_2d_int
    INTEGER, POINTER :: DATA(:,:) ! idx, blk
  END TYPE
  TYPE var_data_2d_wp
    REAL(wp), POINTER :: DATA(:,:) ! idx, blk
  END TYPE
  TYPE var_data_3d_int
    INTEGER, POINTER :: DATA(:,:,:) ! idx, lvl, blk / idx, blk, time
  END TYPE
  TYPE var_data_3d_wp
    REAL(wp), POINTER :: DATA(:,:,:) ! idx, lvl, blk / idx, blk, time
  END TYPE
  TYPE var_data_4d_int
    INTEGER, POINTER :: DATA(:,:,:,:) ! idx, lvl, blk, time
  END TYPE
  TYPE var_data_4d_wp
    REAL(wp), POINTER :: DATA(:,:,:,:) ! idx, lvl, blk, time
  END TYPE

  TYPE(t_basic_distrib_read_data), TARGET, ALLOCATABLE :: basic_data(:)

CONTAINS

  FUNCTION get_empty_basic_distrib_read_data()

    INTEGER :: get_empty_basic_distrib_read_data

    TYPE(t_basic_distrib_read_data), ALLOCATABLE :: temp_basic_distrib_read_data(:)
    INTEGER :: i, n

    get_empty_basic_distrib_read_data = 0

    DO i = 1, SIZE(basic_data)
      IF (basic_data(i)%n_g == -1) THEN
        get_empty_basic_distrib_read_data = i
        EXIT
      END IF
    END DO

    IF (get_empty_basic_distrib_read_data == 0) THEN

      IF (ALLOCATED(basic_data)) THEN

        n = SIZE(basic_data)

        ALLOCATE(temp_basic_distrib_read_data(n))

        temp_basic_distrib_read_data(:) = basic_data(:)

        DEALLOCATE(basic_data)

      ELSE
        n = 0
        ALLOCATE(temp_basic_distrib_read_data(0))
      END IF

      ALLOCATE(basic_data(n+64))

      DO i = n+1, n+64
        basic_data(i)%n_g = -1
      END DO

      basic_data(1:n) = temp_basic_distrib_read_data(:)

      DEALLOCATE(temp_basic_distrib_read_data)

      get_empty_basic_distrib_read_data = n + 1
    END IF

  END FUNCTION get_empty_basic_distrib_read_data

  !-------------------------------------------------------------------------

  SUBROUTINE init_basic_distrib_read_data(basic_read_data, n_g)

    TYPE(t_basic_distrib_read_data), INTENT(inout) :: basic_read_data
    INTEGER, INTENT(in) :: n_g

    INTEGER :: n_io_processes, io_process_stride

    ! data required for the communication pattern
    INTEGER :: i, n
    INTEGER, ALLOCATABLE :: idxmap(:,:)

    basic_read_data%n_g = n_g
    CALL distrib_nf_io_rank_distribution(n_io_processes, io_process_stride)

    basic_read_data%io_process_stride = io_process_stride
    basic_read_data%num_io_processes = n_io_processes

    ALLOCATE(basic_read_data%glb2loc_index)

    CALL init_glb2loc_index_lookup(basic_read_data%glb2loc_index, n_g)

    ! if the process takes part in the reading
    IF (distrib_nf_rank_does_io(n_io_processes, io_process_stride)) THEN
      basic_read_data%io_chunk = uniform_partition(extent(1, n_g), &
           n_io_processes, p_pe_work / io_process_stride + 1)
      n = basic_read_data%io_chunk%size
      ALLOCATE(idxmap(n, 2))
      DO i = 1, n
        idxmap(i, 1) = i + basic_read_data%io_chunk%first - 1
        idxmap(i, 2) = i
      END DO
      CALL set_inner_glb_index(basic_read_data%glb2loc_index, &
           idxmap(:, 1), idxmap(:, 2))
    ELSE
      basic_read_data%io_chunk = extent(1, 0)
    END IF

  END SUBROUTINE init_basic_distrib_read_data

  !-------------------------------------------------------------------------

  FUNCTION get_basic_distrib_read_data(n_g)

    INTEGER, INTENT(in) :: n_g
    INTEGER :: get_basic_distrib_read_data

    INTEGER :: i

    get_basic_distrib_read_data = 0

    IF (.NOT. ALLOCATED(basic_data)) THEN

      ALLOCATE(basic_data(64))
      basic_data(:)%n_g = -1
    END IF

    DO i = 1, SIZE(basic_data)
      IF (basic_data(i)%n_g == n_g) THEN
        get_basic_distrib_read_data = i
        EXIT
      END IF
    END DO

    IF (get_basic_distrib_read_data == 0) THEN
      get_basic_distrib_read_data = get_empty_basic_distrib_read_data()
      CALL init_basic_distrib_read_data(basic_data(get_basic_distrib_read_data), n_g)
      basic_data(get_basic_distrib_read_data)%n_ref = 1
    ELSE
      basic_data(get_basic_distrib_read_data)%n_ref = &
        & basic_data(get_basic_distrib_read_data)%n_ref + 1
    END IF

  END FUNCTION get_basic_distrib_read_data

  !-------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE distrib_nf_io_rank_distribution(&
       n_io_processes, io_process_stride)
    INTEGER, INTENT(out) :: n_io_processes, io_process_stride
    INTEGER :: io_process_rotate

    IF (config_io_process_stride > 0) THEN
      io_process_stride = MAX(1, MODULO(config_io_process_stride, p_n_work))
      io_process_rotate = MODULO(config_io_process_rotate, io_process_stride)
      n_io_processes = (p_n_work - io_process_rotate + io_process_stride - 1)&
           / io_process_stride
    ELSE
      n_io_processes = NINT(SQRT(REAL(p_n_work)))
      io_process_stride = (p_n_work + n_io_processes - 1) / n_io_processes
    END IF
  END SUBROUTINE distrib_nf_io_rank_distribution

  FUNCTION distrib_nf_rank_does_io(n_io_processes, io_process_stride) &
       RESULT (do_open)
    INTEGER, INTENT(in) :: n_io_processes, io_process_stride
    LOGICAL :: do_open
    do_open = MOD(p_pe_work, io_process_stride) &
         == MODULO(config_io_process_rotate, io_process_stride)
  END FUNCTION distrib_nf_rank_does_io

  SUBROUTINE distrib_read_compute_owner(n, glb_index, owner, basic_io_data)
    INTEGER, INTENT(in) :: n
    INTEGER, INTENT(in) :: glb_index(n)
    INTEGER, INTENT(out) :: owner(n)
    TYPE(t_basic_distrib_read_data), INTENT(in) :: basic_io_data

    owner = (partidx_of_elem_uniform_deco(extent(1, basic_io_data%n_g), &
         basic_io_data%num_io_processes, glb_index) - 1) &
         * basic_io_data%io_process_stride &
         + MODULO(config_io_process_rotate, basic_io_data%io_process_stride)
  END SUBROUTINE distrib_read_compute_owner

  !-------------------------------------------------------------------------

  INTEGER FUNCTION distrib_nf_open(path)

    CHARACTER(LEN=*), INTENT(in) :: path

    INTEGER :: n_io_processes, io_process_stride

    CALL distrib_nf_io_rank_distribution(n_io_processes, io_process_stride)
    IF (distrib_nf_rank_does_io(n_io_processes, io_process_stride)) THEN
      CALL nf(nf_open(path, nf_nowrite, distrib_nf_open))
    ELSE
      distrib_nf_open = -1
    END IF

  END FUNCTION distrib_nf_open

  !-------------------------------------------------------------------------

  SUBROUTINE distrib_nf_close(ncid)

    INTEGER, INTENT(in) :: ncid

    INTEGER :: n_io_processes, io_process_stride

    CALL distrib_nf_io_rank_distribution(n_io_processes, io_process_stride)

    IF (distrib_nf_rank_does_io(n_io_processes, io_process_stride)) THEN
      CALL nf(nf_close(ncid))
    END IF

  END SUBROUTINE distrib_nf_close

  !-------------------------------------------------------------------------

  SUBROUTINE setup_distrib_read(n_g, decomp_info, io_data)

    INTEGER, INTENT(in) :: n_g
    TYPE(t_grid_domain_decomp_info), INTENT(in) :: decomp_info
    TYPE(t_distrib_read_data), INTENT(inout) :: io_data

    INTEGER :: n, n_inner, i
    INTEGER, ALLOCATABLE :: owner(:)
    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data

    io_data%basic_data_index = get_basic_distrib_read_data(n_g)
    basic_io_data => basic_data(io_data%basic_data_index)

    n = SIZE(decomp_info%glb_index)

    ALLOCATE(owner(n))

    CALL distrib_read_compute_owner(n, decomp_info%glb_index(:), &
      & owner(:), basic_io_data)
    n_inner = SIZE(basic_io_data%glb2loc_index%inner_glb_index, 1)
    CALL setup_comm_pattern(n, owner(:), decomp_info%glb_index(:), &
      & basic_io_data%glb2loc_index, n_inner, &
      & (/(p_pe_work, i = 1, n_inner)/), &
      & basic_io_data%glb2loc_index%inner_glb_index, io_data%redistrib_pattern)

    DEALLOCATE(owner)

  END SUBROUTINE setup_distrib_read

  !-------------------------------------------------------------------------

  SUBROUTINE delete_basic_distrib_read(idx)

    INTEGER, INTENT(in) :: idx

    IF (.NOT. ALLOCATED(basic_data)) &
      & CALL finish("delete_basic_distrib_read", "basic_data is noch allocated")
    IF (idx > SIZE(basic_data)) &
      & CALL finish("delete_basic_distrib_read", "invalid basic_data index")
    IF (basic_data(idx)%n_ref < 1) &
      & CALL finish("delete_basic_distrib_read", "invalid reference counter")

    basic_data(idx)%n_ref = basic_data(idx)%n_ref - 1

    ! some other distrib_read is still using this basic_distrib_read_data
    IF (basic_data(idx)%n_ref > 0) RETURN

    basic_data(idx)%n_g = -1
    CALL deallocate_glb2loc_index_lookup(basic_data(idx)%glb2loc_index)
    DEALLOCATE(basic_data(idx)%glb2loc_index)

  END SUBROUTINE delete_basic_distrib_read

  !-------------------------------------------------------------------------

  SUBROUTINE delete_distrib_read(io_data)

    TYPE(t_distrib_read_data), INTENT(inout) :: io_data

    CALL delete_comm_pattern(io_data%redistrib_pattern)
    CALL delete_basic_distrib_read(io_data%basic_data_index)

  END SUBROUTINE delete_distrib_read

  !-------------------------------------------------------------------------

  ELEMENTAL INTEGER FUNCTION distrib_read_get_buffer_size(io_data)

    TYPE(t_distrib_read_data), INTENT(in) :: io_data

    distrib_read_get_buffer_size &
      = basic_data(io_data%basic_data_index)%io_chunk%size

  END FUNCTION distrib_read_get_buffer_size

  !-------------------------------------------------------------------------

  SUBROUTINE check_basic_data_index(io_data)

    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)

    IF (SIZE(io_data) > 1) THEN

      IF (.NOT. ALL(io_data(2:)%basic_data_index == &
        & io_data(1)%basic_data_index)) &
        & CALL finish("check_basic_data_index", "basic_data_index do not match")
    END IF
  END SUBROUTINE check_basic_data_index

  !-------------------------------------------------------------------------

  SUBROUTINE distrib_read_int_2d(ncid, var_name, var_data, io_data)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, TARGET, INTENT(INOUT) :: var_data(:,:) ! idx, blk
    TYPE(t_distrib_read_data), INTENT(IN) :: io_data

    TYPE(var_data_2d_int) :: var_data_(1)

    var_data_(1)%data => var_data

    CALL distrib_read_int_2d_multi_var(ncid, var_name, var_data_, (/io_data/))

  END SUBROUTINE distrib_read_int_2d

  SUBROUTINE distrib_read_int_2d_multi_var(ncid, var_name, var_data, io_data)

    INTEGER, INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: var_name
    TYPE(var_data_2d_int), INTENT(inout) :: var_data(:)
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)

    INTEGER, ALLOCATABLE :: local_buffer_1d(:)
    INTEGER, ALLOCATABLE :: local_buffer_2d(:, :)
    INTEGER :: buffer_size, buffer_nblks
    INTEGER :: varid, i, var_type

    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data

    IF (SIZE(io_data) == 0) RETURN

    IF (SIZE(var_data) < SIZE(io_data)) &
      & CALL finish("distrib_read_int_2d_multi_var", "var_data too small")

    CALL check_basic_data_index(io_data(:))
    basic_io_data => basic_data(io_data(1)%basic_data_index)

    buffer_size = MAXVAL(distrib_read_get_buffer_size(io_data(:)))
    buffer_nblks = (buffer_size + nproma - 1) / nproma
    ALLOCATE(local_buffer_1d(buffer_size), &
      &      local_buffer_2d(nproma, buffer_nblks))

    ! read data from file into io decomposition
    IF (basic_io_data%io_chunk%size > 0) THEN
      CALL nf(nf_inq_varid(ncid, var_name, varid))
      CALL nf(nf_inq_vartype(ncid, varid, var_type))

      IF (var_type /= NF_INT) &
        CALL finish("distrib_read_int_2d_multi_var", "invalid var_type")

      ! only read io_decomp part
      CALL nf(nf_get_vara_int(ncid, varid, (/basic_io_data%io_chunk%first/), &
        & (/basic_io_data%io_chunk%size/), local_buffer_1d(:)))

      DO i = 1, SIZE(local_buffer_1d(:))
        local_buffer_2d(idx_no(i),blk_no(i)) = local_buffer_1d(i)
      END DO
    END IF

    ! redistribute data into final decomposition
    DO i = 1, SIZE(io_data)
      CALL exchange_data(io_data(i)%redistrib_pattern, var_data(i)%DATA, &
        & local_buffer_2d)
    END DO

    DEALLOCATE(local_buffer_1d, local_buffer_2d)

  END SUBROUTINE distrib_read_int_2d_multi_var

  !-------------------------------------------------------------------------

  SUBROUTINE distrib_read_real_2d(ncid, var_name, var_data, io_data)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    REAL(wp), TARGET, INTENT(INOUT) :: var_data(:,:) ! idx, blk
    TYPE(t_distrib_read_data), INTENT(IN) :: io_data

    TYPE(var_data_2d_wp) :: var_data_(1)

    var_data_(1)%data => var_data

    CALL distrib_read_real_2d_multi_var(ncid, var_name, var_data_, (/io_data/))

  END SUBROUTINE distrib_read_real_2d

  SUBROUTINE distrib_read_real_2d_multi_var(ncid, var_name, var_data, io_data)

    INTEGER, INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: var_name
    TYPE(var_data_2d_wp), INTENT(inout) :: var_data(:)
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)

    REAL(wp), ALLOCATABLE :: local_buffer_1d(:)
    REAL(wp), ALLOCATABLE :: local_buffer_2d(:, :)
    INTEGER :: buffer_size, buffer_nblks
    INTEGER :: varid, i

    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data

    IF (SIZE(io_data) == 0) RETURN

    IF (SIZE(var_data) < SIZE(io_data)) &
      & CALL finish("distrib_read_real_2d_multi_var", "var_data too small")

    CALL check_basic_data_index(io_data(:))
    basic_io_data => basic_data(io_data(1)%basic_data_index)

    buffer_size = MAXVAL(distrib_read_get_buffer_size(io_data(:)))
    buffer_nblks = (buffer_size + nproma - 1) / nproma
    ALLOCATE(local_buffer_1d(buffer_size), &
      &      local_buffer_2d(nproma, buffer_nblks))

    ! read data from file into io decomposition
    IF (basic_io_data%io_chunk%size > 0) THEN
      CALL nf(nf_inq_varid(ncid, var_name, varid))
      ! only read io_decomp part
      CALL nf(nf_get_vara_double(ncid, varid, (/basic_io_data%io_chunk%first/), &
        & (/basic_io_data%io_chunk%size/), local_buffer_1d(:)))
      DO i = 1, SIZE(local_buffer_1d(:))
        local_buffer_2d(idx_no(i),blk_no(i)) = local_buffer_1d(i)
      END DO
    END IF

    ! redistribute data into final decomposition
    DO i = 1, SIZE(io_data)
      CALL exchange_data(io_data(i)%redistrib_pattern, var_data(i)%DATA, &
        & local_buffer_2d)
    END DO

    DEALLOCATE(local_buffer_1d, local_buffer_2d)

  END SUBROUTINE distrib_read_real_2d_multi_var

  !-------------------------------------------------------------------------

  SUBROUTINE distrib_read_int_3d(ncid, var_name, var_data, ext_dim_size, &
    &                            dim_order, io_data, start_ext_dim, &
    &                            end_ext_dim)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, TARGET, INTENT(INOUT) :: var_data(:,:,:)
    INTEGER, INTENT(IN) :: ext_dim_size, dim_order
    TYPE(t_distrib_read_data), INTENT(IN) :: io_data
    INTEGER, OPTIONAL, INTENT(IN) :: start_ext_dim, end_ext_dim
    TYPE(var_data_3d_int) :: var_data_(1)

    var_data_(1)%data => var_data

    CALL distrib_read_int_3d_multi_var(ncid, var_name, var_data_, &
      &                                ext_dim_size, dim_order, (/io_data/), &
      &                                start_ext_dim, end_ext_dim)

  END SUBROUTINE distrib_read_int_3d

  SUBROUTINE distrib_read_int_3d_multi_var(ncid, var_name, var_data, &
    &                                      ext_dim_size, dim_order, io_data, &
    &                                      start_ext_dim, end_ext_dim)

    INTEGER, INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: var_name
    TYPE(var_data_3d_int), INTENT(inout) :: var_data(:)
    INTEGER, INTENT(IN) :: ext_dim_size, dim_order
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)
    INTEGER, OPTIONAL, INTENT(IN) :: start_ext_dim, end_ext_dim

    INTEGER, ALLOCATABLE :: local_buffer_2d(:,:) ! (n io points, ext_dim_size)
    INTEGER, ALLOCATABLE :: local_buffer_3d(:,:,:)
    INTEGER :: buffer_size, buffer_nblks
    INTEGER :: varid, i, j, var_type
    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data

    IF (SIZE(io_data) == 0) RETURN

    IF (SIZE(var_data) < SIZE(io_data)) &
      & CALL finish("distrib_read_int_3d_multi_var", "var_data too small")

    IF (PRESENT(start_ext_dim) .NEQV. PRESENT(end_ext_dim)) &
      CALL finish("distrib_read_real_3d_multi_var", "invalid arguments")

    IF (PRESENT(start_ext_dim)) THEN
      IF ((end_ext_dim - start_ext_dim + 1) /= ext_dim_size) &
        CALL finish("distrib_read_real_3d_multi_var", &
          &         "invalid arguments; (end_ext_dim - start_ext_dim + 1) /= ext_dim_size")
    END IF

    CALL check_basic_data_index(io_data(:))
    basic_io_data => basic_data(io_data(1)%basic_data_index)

    buffer_size = MAXVAL(distrib_read_get_buffer_size(io_data(:)))
    buffer_nblks = (buffer_size + nproma - 1) / nproma
    ALLOCATE(local_buffer_2d(buffer_size, ext_dim_size))
    SELECT CASE(dim_order)
      CASE(idx_blk_time)
        ALLOCATE(local_buffer_3d(nproma, buffer_nblks, ext_dim_size))
      CASE(idx_lvl_blk)
        ALLOCATE(local_buffer_3d(nproma, ext_dim_size, buffer_nblks))
      CASE DEFAULT
        CALL finish("distrib_read_real_3d_multi_var", "invalid dim_order")
    END SELECT

    IF (basic_io_data%io_chunk%size > 0) THEN

      CALL nf(nf_inq_varid(ncid, var_name, varid))
      CALL nf(nf_inq_vartype(ncid, varid, var_type))

      IF (var_type /= NF_INT) &
        CALL finish("distrib_read_int_3d_multi_var", "invalid var_type")

      ! only read io_decomp part
      IF (PRESENT(start_ext_dim)) THEN
        CALL nf(nf_get_vara_int(ncid, varid, &
          &                     (/basic_io_data%io_chunk%first, start_ext_dim/), &
          &                     (/basic_io_data%io_chunk%size, ext_dim_size/), &
          &                     local_buffer_2d(:,:)))
      ELSE
        CALL nf(nf_get_vara_int(ncid, varid, (/basic_io_data%io_chunk%first, 1/), &
          &                     (/basic_io_data%io_chunk%size, ext_dim_size/), &
          &                     local_buffer_2d(:,:)))
      END IF
    END IF

    IF (dim_order == idx_blk_time) THEN
      DO j = 1, ext_dim_size
        DO i = 1, basic_io_data%io_chunk%size
          local_buffer_3d(idx_no(i),blk_no(i), j) = local_buffer_2d(i, j)
        END DO
      END DO

      DO i = 1, SIZE(io_data)
        DO j = 1, ext_dim_size
          CALL exchange_data(io_data(i)%redistrib_pattern, &
            & var_data(i)%DATA(:,:,j), &
            & local_buffer_3d(:,:,j))
!            & var_data(i)%DATA(:,:,LBOUND(var_data(i)%DATA,3)+j-1), &
        END DO
      END DO
    ELSE
      DO j = 1, ext_dim_size
        DO i = 1, basic_io_data%io_chunk%size
          local_buffer_3d(idx_no(i), j, blk_no(i)) = local_buffer_2d(i, j)
        END DO
      END DO

      DO i = 1, SIZE(io_data)
        CALL exchange_data(io_data(i)%redistrib_pattern, &
          & var_data(i)%DATA(:,:,:), &
          & local_buffer_3d(:,:,:))
!           & var_data(i)%DATA(:,LBOUND(var_data(i)%DATA,2): &
!           &                    UBOUND(var_data(i)%DATA,2),:), &
      END DO
    END IF

    DEALLOCATE(local_buffer_2d, local_buffer_3d)

  END SUBROUTINE distrib_read_int_3d_multi_var

  !-------------------------------------------------------------------------

  SUBROUTINE distrib_read_real_3d(ncid, var_name, var_data, ext_dim_size, &
    &                             dim_order, io_data, start_ext_dim, &
    &                             end_ext_dim)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    REAL(wp), TARGET, INTENT(INOUT) :: var_data(:,:,:)
    INTEGER, INTENT(IN) :: ext_dim_size, dim_order
    TYPE(t_distrib_read_data), INTENT(IN) :: io_data
    INTEGER, OPTIONAL, INTENT(IN) :: start_ext_dim, end_ext_dim
    TYPE(var_data_3d_wp) :: var_data_(1)

    var_data_(1)%data => var_data

    CALL distrib_read_real_3d_multi_var(ncid, var_name, var_data_, &
      &                                 ext_dim_size, dim_order, (/io_data/), &
      &                                 start_ext_dim, end_ext_dim)

  END SUBROUTINE distrib_read_real_3d

  SUBROUTINE distrib_read_real_3d_multi_var(ncid, var_name, var_data, &
    &                                       ext_dim_size, dim_order, io_data, &
    &                                       start_ext_dim, end_ext_dim)

    INTEGER, INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: var_name
    TYPE(var_data_3d_wp), INTENT(inout) :: var_data(:)
    INTEGER, INTENT(IN) :: ext_dim_size, dim_order
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)
    INTEGER, OPTIONAL, INTENT(IN) :: start_ext_dim, end_ext_dim

    REAL(wp), ALLOCATABLE :: local_buffer_2d(:,:) ! (n io points, ext_dim_size)
    REAL(wp), ALLOCATABLE :: local_buffer_3d(:,:,:)
    INTEGER :: buffer_size, buffer_nblks
    INTEGER :: varid, i, j
    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data

    IF (SIZE(io_data) == 0) RETURN

    IF (SIZE(var_data) < SIZE(io_data)) &
      CALL finish("distrib_read_real_3d_multi_var", "var_data too small")

    IF (PRESENT(start_ext_dim) .NEQV. PRESENT(end_ext_dim)) &
      CALL finish("distrib_read_real_3d_multi_var", "invalid arguments")

    IF (PRESENT(start_ext_dim)) THEN
      IF ((end_ext_dim - start_ext_dim + 1) /= ext_dim_size) &
        CALL finish("distrib_read_real_3d_multi_var", &
          &         "invalid arguments; (end_ext_dim - start_ext_dim + 1) /= ext_dim_size")
    END IF

    CALL check_basic_data_index(io_data(:))
    basic_io_data => basic_data(io_data(1)%basic_data_index)

    buffer_size = MAXVAL(distrib_read_get_buffer_size(io_data(:)))
    buffer_nblks = (buffer_size + nproma - 1) / nproma
    ALLOCATE(local_buffer_2d(buffer_size, ext_dim_size))
    SELECT CASE(dim_order)
      CASE(idx_blk_time)
        ALLOCATE(local_buffer_3d(nproma, buffer_nblks, ext_dim_size))
      CASE(idx_lvl_blk)
        ALLOCATE(local_buffer_3d(nproma, ext_dim_size, buffer_nblks))
      CASE DEFAULT
        CALL finish("distrib_read_real_3d_multi_var", "invalid dim_order")
    END SELECT

    IF (basic_io_data%io_chunk%size > 0) THEN

      CALL nf(nf_inq_varid(ncid, var_name, varid))
      ! only read io_decomp part
      IF (PRESENT(start_ext_dim)) THEN
        CALL nf(nf_get_vara_double(ncid, varid, &
          &                        (/basic_io_data%io_chunk%first, start_ext_dim/), &
          &                        (/basic_io_data%io_chunk%size, ext_dim_size/), &
          &                        local_buffer_2d(:,:)))
      ELSE
        CALL nf(nf_get_vara_double(ncid, varid, (/basic_io_data%io_chunk%first, 1/), &
          &                        (/basic_io_data%io_chunk%size, ext_dim_size/), &
          &                        local_buffer_2d(:,:)))
      END IF
    END IF

    IF (dim_order == idx_blk_time) THEN
      DO j = 1, ext_dim_size
        DO i = 1, basic_io_data%io_chunk%size
          local_buffer_3d(idx_no(i),blk_no(i), j) = local_buffer_2d(i, j)
        END DO
      END DO

      DO i = 1, SIZE(io_data)
        DO j = 1, ext_dim_size
          CALL exchange_data(io_data(i)%redistrib_pattern, &
            & var_data(i)%DATA(:,:,LBOUND(var_data(i)%DATA,3)+j-1), &
            & local_buffer_3d(:,:,j))
        END DO
      END DO
    ELSE
      DO j = 1, ext_dim_size
        DO i = 1, basic_io_data%io_chunk%size
          local_buffer_3d(idx_no(i), j, blk_no(i)) = local_buffer_2d(i, j)
        END DO
      END DO

      DO i = 1, SIZE(io_data)
        CALL exchange_data(io_data(i)%redistrib_pattern, &
          & var_data(i)%DATA(:,LBOUND(var_data(i)%DATA,2): &
          &                    UBOUND(var_data(i)%DATA,2),:), &
          & local_buffer_3d(:,:,:))
      END DO
    END IF

    DEALLOCATE(local_buffer_2d, local_buffer_3d)

  END SUBROUTINE distrib_read_real_3d_multi_var

  !-------------------------------------------------------------------------

  SUBROUTINE distrib_read_int_4d(ncid, var_name, var_data, ext_dim_size, &
    &                            io_data)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    INTEGER, TARGET, INTENT(INOUT) :: var_data(:,:,:,:)
    INTEGER, INTENT(IN) :: ext_dim_size(2)
    TYPE(t_distrib_read_data), INTENT(IN) :: io_data
    TYPE(var_data_4d_int) :: var_data_(1)

    var_data_(1)%data => var_data

    CALL distrib_read_int_4d_multi_var(ncid, var_name, var_data_, &
      &                                 ext_dim_size, (/io_data/))

  END SUBROUTINE distrib_read_int_4d

  SUBROUTINE distrib_read_int_4d_multi_var(ncid, var_name, var_data, &
    &                                      ext_dim_size, io_data)

    INTEGER, INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: var_name
    TYPE(var_data_4d_int), INTENT(inout) :: var_data(:)
    INTEGER, INTENT(IN) :: ext_dim_size(2)
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)

    INTEGER, ALLOCATABLE :: local_buffer_3d(:,:,:) ! (n io points,
                                                   !  ext_dim_size(1),
                                                   !  ext_dim_size(2))
    INTEGER, ALLOCATABLE :: local_buffer_4d(:,:,:,:)
    INTEGER :: buffer_size, buffer_nblks
    INTEGER :: varid, i, j, k, var_type

    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data

    IF (SIZE(io_data) == 0) RETURN

    IF (SIZE(var_data) < SIZE(io_data)) &
      & CALL finish("distrib_read_int_4d_multi_var", "var_data too small")

    CALL check_basic_data_index(io_data(:))
    basic_io_data => basic_data(io_data(1)%basic_data_index)

    buffer_size = MAXVAL(distrib_read_get_buffer_size(io_data(:)))
    buffer_nblks = (buffer_size + nproma - 1) / nproma
    ALLOCATE(local_buffer_3d(buffer_size, ext_dim_size(1), ext_dim_size(2)))
    ALLOCATE(local_buffer_4d(nproma, ext_dim_size(1), buffer_nblks, &
      &                      ext_dim_size(2)))

    IF (basic_io_data%io_chunk%size > 0) THEN

      CALL nf(nf_inq_varid(ncid, var_name, varid))
      CALL nf(nf_inq_vartype(ncid, varid, var_type))

      IF (var_type /= NF_INT) &
        CALL finish("distrib_read_int_3d_multi_var", "invalid var_type")

      ! only read io_decomp part
      CALL nf(nf_get_vara_int(ncid, varid, (/basic_io_data%io_chunk%first, 1, 1/), &
        & (/basic_io_data%io_chunk%size, ext_dim_size(1), ext_dim_size(2)/), &
        & local_buffer_3d(:,:,:)))
    END IF

    DO k = 1, ext_dim_size(2)
      DO j = 1, ext_dim_size(1)
        DO i = 1, basic_io_data%io_chunk%size
          local_buffer_4d(idx_no(i),j,blk_no(i), k) = local_buffer_3d(i, j, k)
        END DO
      END DO
    END DO

    DO i = 1, SIZE(io_data)
      DO j = 1, ext_dim_size(2)
        CALL exchange_data(io_data(i)%redistrib_pattern, &
          & var_data(i)%DATA(:,:,:,j), local_buffer_4d(:,:,:,j))
      END DO
    END DO

    DEALLOCATE(local_buffer_3d, local_buffer_4d)

  END SUBROUTINE distrib_read_int_4d_multi_var

  !-------------------------------------------------------------------------

  SUBROUTINE distrib_read_real_4d(ncid, var_name, var_data, ext_dim_size, &
    &                             io_data, start_ext_dim, end_ext_dim)

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: var_name
    REAL(wp), TARGET, INTENT(INOUT) :: var_data(:,:,:,:)
    INTEGER, INTENT(IN) :: ext_dim_size(2)
    TYPE(t_distrib_read_data), INTENT(IN) :: io_data
    INTEGER, OPTIONAL, INTENT(IN) :: start_ext_dim(2), end_ext_dim(2)
    TYPE(var_data_4d_wp) :: var_data_(1)

    var_data_(1)%data => var_data

    CALL distrib_read_real_4d_multi_var(ncid, var_name, var_data_, &
      &                                 ext_dim_size, (/io_data/), &
      &                                 start_ext_dim, end_ext_dim)

  END SUBROUTINE distrib_read_real_4d

  SUBROUTINE distrib_read_real_4d_multi_var(ncid, var_name, var_data, &
    &                                       ext_dim_size, io_data, &
    &                                       start_ext_dim, end_ext_dim)

    INTEGER, INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: var_name
    TYPE(var_data_4d_wp), INTENT(inout) :: var_data(:)
    INTEGER, INTENT(IN) :: ext_dim_size(2)
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)
    INTEGER, OPTIONAL, INTENT(IN) :: start_ext_dim(2), end_ext_dim(2)

    REAL(wp), ALLOCATABLE :: local_buffer_3d(:,:,:) ! (n io points,
                                                    !  ext_dim_size(1),
                                                    !  ext_dim_size(2))
    REAL(wp), ALLOCATABLE :: local_buffer_4d(:,:,:,:)
    INTEGER :: buffer_size, buffer_nblks
    INTEGER :: varid, i, j, k

    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data

    IF (SIZE(io_data) == 0) RETURN

    IF (SIZE(var_data) < SIZE(io_data)) &
      & CALL finish("distrib_read_real_4d_multi_var", "var_data too small")

    IF (PRESENT(start_ext_dim) .NEQV. PRESENT(end_ext_dim)) &
      CALL finish("distrib_read_real_4d_multi_var", "invalid arguments")

    IF (PRESENT(start_ext_dim)) THEN
      IF (ANY((end_ext_dim(:) - start_ext_dim(:) + 1) /= ext_dim_size(:))) &
        CALL finish("distrib_read_real_4d_multi_var", &
          "invalid arguments; (end_ext_dim - start_ext_dim + 1) /= ext_dim_size(2)")
    END IF

    CALL check_basic_data_index(io_data(:))
    basic_io_data => basic_data(io_data(1)%basic_data_index)

    buffer_size = MAXVAL(distrib_read_get_buffer_size(io_data(:)))
    buffer_nblks = (buffer_size + nproma - 1) / nproma
    ALLOCATE(local_buffer_3d(buffer_size, ext_dim_size(1), ext_dim_size(2)))
    ALLOCATE(local_buffer_4d(nproma, ext_dim_size(1), buffer_nblks, &
      &                      ext_dim_size(2)))

    IF (basic_io_data%io_chunk%size > 0) THEN

      CALL nf(nf_inq_varid(ncid, var_name, varid))
      ! only read io_decomp part
      IF (PRESENT(start_ext_dim)) THEN
        CALL nf(nf_get_vara_double(ncid, varid, (/basic_io_data%io_chunk%first, &
          &                                       start_ext_dim(1), &
          &                                       start_ext_dim(2)/), &
          & (/basic_io_data%io_chunk%size, ext_dim_size(1), ext_dim_size(2)/), &
          & local_buffer_3d(:,:,:)))
      ELSE
        CALL nf(nf_get_vara_double(ncid, varid, (/basic_io_data%io_chunk%first, 1, 1/), &
          & (/basic_io_data%io_chunk%size, ext_dim_size(1), ext_dim_size(2)/), &
          & local_buffer_3d(:,:,:)))
      END IF
    END IF

    DO k = 1, ext_dim_size(2)
      DO j = 1, ext_dim_size(1)
        DO i = 1, basic_io_data%io_chunk%size
          local_buffer_4d(idx_no(i),j,blk_no(i), k) = local_buffer_3d(i, j, k)
        END DO
      END DO
    END DO

    DO i = 1, SIZE(io_data)
      DO j = 1, ext_dim_size(2)
        CALL exchange_data(io_data(i)%redistrib_pattern, &
          var_data(i)%DATA(:, &
          LBOUND(var_data(i)%DATA,2):UBOUND(var_data(i)%DATA,2),:, &
          LBOUND(var_data(i)%DATA,4)+j-1), local_buffer_4d(:,:,:,j))
      END DO
    END DO

    DEALLOCATE(local_buffer_3d, local_buffer_4d)

  END SUBROUTINE distrib_read_real_4d_multi_var

  !-------------------------------------------------------------------------

  SUBROUTINE distrib_inq_var_dims(file_id, var_name, var_ndims, var_dimlen)
    INTEGER, INTENT(IN) :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: var_name

    INTEGER, INTENT(OUT) :: var_ndims
    INTEGER, INTENT(OUT) :: var_dimlen(:)

    INTEGER :: varid, temp(1), i
    INTEGER :: temp_var_dimlen(NF_MAX_VAR_DIMS), var_dimids(NF_MAX_VAR_DIMS)

    IF ( p_pe_work == 0 ) THEN
      CALL nf(nf_inq_varid(file_id, var_name, varid))
      CALL nf(nf_inq_varndims(file_id, varid, var_ndims))
      CALL nf(nf_inq_vardimid(file_id, varid, var_dimids))
      DO i=1, var_ndims
        CALL nf(nf_inq_dimlen(file_id, var_dimids(i), temp_var_dimlen(i)))
      ENDDO
      temp(1) = var_ndims
    END IF

#ifndef NOMPI
    CALL p_bcast(temp, 0, p_comm_work)
    CALL p_bcast(var_ndims, 0, p_comm_work)
    CALL p_bcast(temp_var_dimlen(1:var_ndims), 0, p_comm_work)
    temp(1) = var_ndims
#endif

    IF (SIZE(var_dimlen) < var_ndims) &
      CALL finish("distrib_inq_var_dims", &
        &         "array size of argument var_dimlen is too small")

    var_dimlen(1:var_ndims) = temp_var_dimlen(1:var_ndims)

  END SUBROUTINE distrib_inq_var_dims

  !-------------------------------------------------------------------------

  SUBROUTINE nf(STATUS, warnonly, silent)

    INTEGER, INTENT(in)           :: STATUS
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
        CALL message('mo_read_netcdf_distributed netCDF error', &
          &          nf_strerror(STATUS), level=em_warn)
      ELSE
        CALL finish('mo_read_netcdf_distributed netCDF error', &
          &         nf_strerror(STATUS))
      ENDIF
    ENDIF

  END SUBROUTINE nf

END MODULE mo_read_netcdf_distributed
