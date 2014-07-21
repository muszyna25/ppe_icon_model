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
  
  USE mo_kind, ONLY: dp, wp
  USE mo_exception, ONLY: finish, message, em_warn
  USE mo_mpi, ONLY: p_n_work, p_pe_work
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info, &
    & t_glb2loc_index_lookup, &
    & init_glb2loc_index_lookup, &
    & set_inner_glb_index
  USE mo_communication, ONLY: t_comm_pattern, idx_no, blk_no, &
    & setup_comm_pattern, delete_comm_pattern, &
    & exchange_data
  USE mo_alloc_patches, ONLY: deallocate_glb2loc_index_lookup
  USE mo_parallel_config, ONLY: nproma
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: distrib_read
  PUBLIC :: distrib_nf_open
  PUBLIC :: distrib_nf_close
  PUBLIC :: setup_distrib_read
  PUBLIC :: delete_distrib_read
  PUBLIC :: distrib_read_get_buffer_size
  PUBLIC :: t_distrib_read_data
  PUBLIC :: var_data_1d_int, var_data_1d_wp
  PUBLIC :: var_data_2d_int, var_data_2d_wp
  PUBLIC :: idx_blk_lvl, idx_lvl_blk
  
  INCLUDE 'netcdf.inc'
  
  INTERFACE distrib_read
    MODULE PROCEDURE distrib_read_int_1d_multi_var
    MODULE PROCEDURE distrib_read_int_2d_multi_var
    MODULE PROCEDURE distrib_read_real_1d_multi_var
    MODULE PROCEDURE distrib_read_real_2d_multi_var
  END INTERFACE distrib_read
  
  INTEGER, PARAMETER :: nf_read = nf_nowrite
  INTEGER, PARAMETER :: idx_lvl_blk = 1
  INTEGER, PARAMETER :: idx_blk_lvl = 2
  
  !modules interface-------------------------------------------
  !subroutines
  
  TYPE t_basic_distrib_read_data
    INTEGER :: n_g ! global number of points (-1 if unused)
    INTEGER :: start, COUNT ! io decomposition
    TYPE(t_glb2loc_index_lookup), POINTER :: glb2loc_index
    INTEGER :: num_points_per_io_process, io_process_stride
    
    INTEGER :: n_ref ! number of times this data is referenced
  END TYPE t_basic_distrib_read_data
  
  TYPE t_distrib_read_data
    
    INTEGER :: basic_data_index
    TYPE(t_comm_pattern) :: redistrib_pattern
    
  END TYPE t_distrib_read_data
  
  TYPE var_data_1d_int
    INTEGER, POINTER :: DATA(:,:) ! idx, blk
  END TYPE
  TYPE var_data_1d_wp
    REAL(wp), POINTER :: DATA(:,:) ! idx, blk
  END TYPE
  TYPE var_data_2d_int
    INTEGER, POINTER :: DATA(:,:,:) ! idx, lvl, blk
  END TYPE
  TYPE var_data_2d_wp
    REAL(wp), POINTER :: DATA(:,:,:) ! idx, lvl, blk
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
    INTEGER :: num_points_per_io_process
    
    ! data required for the communication pattern
    INTEGER :: i
    
    basic_read_data%n_g = n_g
    n_io_processes = NINT(SQRT(REAL(p_n_work)))
    io_process_stride = (p_n_work + n_io_processes - 1) / n_io_processes
    num_points_per_io_process = (n_g + n_io_processes - 1) / n_io_processes
    
    basic_read_data%num_points_per_io_process = num_points_per_io_process
    basic_read_data%io_process_stride = io_process_stride
    
    ALLOCATE(basic_read_data%glb2loc_index)

    CALL init_glb2loc_index_lookup(basic_read_data%glb2loc_index, n_g)
    
    ! if the process takes part in the reading
    IF (MOD(p_pe_work, io_process_stride) == 0) THEN
      basic_read_data%start = (p_pe_work / n_io_processes) * num_points_per_io_process + 1
      basic_read_data%COUNT = MIN(num_points_per_io_process, n_g - basic_read_data%start + 1)
      
      CALL set_inner_glb_index(basic_read_data%glb2loc_index, (/(i + basic_read_data%start - 1, &
        & i = 1, &
        & basic_read_data%COUNT)/), &
        & (/(i, i = 1, basic_read_data%COUNT)/))
    ELSE
      basic_read_data%start = 1
      basic_read_data%COUNT = 0
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
  
  INTEGER FUNCTION distrib_nf_open(path)
    
    CHARACTER(LEN=*), INTENT(in) :: path
    
    INTEGER :: n_io_processes, io_process_stride
    
    n_io_processes = NINT(SQRT(REAL(p_n_work)))
    io_process_stride = (p_n_work + n_io_processes - 1) / n_io_processes
    
    IF (MOD(p_pe_work, io_process_stride) == 0) THEN
      CALL nf(nf_open(path, nf_nowrite, distrib_nf_open))
    ELSE
      distrib_nf_open = -1
    END IF
    
  END FUNCTION distrib_nf_open
  
  !-------------------------------------------------------------------------
  
  SUBROUTINE distrib_nf_close(ncid)
    
    INTEGER, INTENT(in) :: ncid
    
    INTEGER :: n_io_processes, io_process_stride
    
    n_io_processes = NINT(SQRT(REAL(p_n_work)))
    io_process_stride = (p_n_work + n_io_processes - 1) / n_io_processes
    
    IF (MOD(p_pe_work, io_process_stride) == 0) THEN
      CALL nf(nf_close(ncid))
    END IF
    
  END SUBROUTINE distrib_nf_close
  
  !-------------------------------------------------------------------------
  
  SUBROUTINE setup_distrib_read(n_g, decomp_info, io_data)
    
    INTEGER, INTENT(in) :: n_g
    TYPE(t_grid_domain_decomp_info), INTENT(in) :: decomp_info
    TYPE(t_distrib_read_data), INTENT(inout) :: io_data
    
    INTEGER :: n
    INTEGER, ALLOCATABLE :: owner(:)
    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data
    
    io_data%basic_data_index = get_basic_distrib_read_data(n_g)
    basic_io_data => basic_data(io_data%basic_data_index)
    
    n = SIZE(decomp_info%glb_index)
    
    ALLOCATE(owner(n))
    
    CALL distrib_read_compute_owner(n, decomp_info%glb_index(:), &
      & owner(:), basic_io_data)
    CALL setup_comm_pattern(n, owner(:), decomp_info%glb_index(:), &
      & basic_io_data%glb2loc_index, &
      & io_data%redistrib_pattern)
    
    DEALLOCATE(owner)
    
  END SUBROUTINE setup_distrib_read
  
  !-------------------------------------------------------------------------
  
  SUBROUTINE distrib_read_compute_owner(n, glb_index, owner, basic_io_data)
    
    INTEGER, INTENT(in) :: n
    INTEGER, INTENT(in) :: glb_index(n)
    INTEGER, INTENT(out) :: owner(n)
    TYPE(t_basic_distrib_read_data), INTENT(in) :: basic_io_data
    
    owner(:) = ((glb_index(:) - 1) / basic_io_data%num_points_per_io_process) &
      & * basic_io_data%io_process_stride
    
  END SUBROUTINE distrib_read_compute_owner
  
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
  
  INTEGER FUNCTION distrib_read_get_buffer_size(io_data)
    
    TYPE(t_distrib_read_data), INTENT(in) :: io_data
    
    distrib_read_get_buffer_size = basic_data(io_data%basic_data_index)%COUNT
    
  END FUNCTION distrib_read_get_buffer_size
  
  !-------------------------------------------------------------------------
  
  SUBROUTINE check_basic_data_index(io_data)
    
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)
    
    IF (.NOT. ALL(io_data(2:)%basic_data_index == &
      & io_data(1)%basic_data_index)) &
      & CALL finish("check_basic_data_index", "basic_data_index do not match")
  END SUBROUTINE check_basic_data_index
  
  !-------------------------------------------------------------------------
  
  SUBROUTINE distrib_read_int_1d_multi_var(ncid, var_name, local_buffer_1d, &
    & local_buffer_2d, var_data, io_data)
    
    INTEGER, INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: var_name
    INTEGER, INTENT(inout) :: local_buffer_1d(:)
    INTEGER, INTENT(inout) :: local_buffer_2d(:, :)
    TYPE(var_data_1d_int), INTENT(inout) :: var_data(:)
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)
    
    INTEGER :: varid, i
    
    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data
    
    IF (SIZE(io_data) == 0) RETURN
    
    IF (SIZE(var_data) < SIZE(io_data)) &
      & CALL finish("distrib_read_int_1d_multi_var", "var_data too small")
    
    CALL check_basic_data_index(io_data(:))
    
    basic_io_data => basic_data(io_data(1)%basic_data_index)
    
    ! read data from file into io decomposition
    IF (basic_io_data%COUNT > 0) THEN
      CALL nf(nf_inq_varid(ncid, var_name, varid))
      
      ! only read io_decomp part
      CALL nf(nf_get_vara_int(ncid, varid, (/basic_io_data%start/), &
        & (/basic_io_data%COUNT/), local_buffer_1d(:)))
      
      DO i = 1, SIZE(local_buffer_1d(:))
        local_buffer_2d(idx_no(i),blk_no(i)) = local_buffer_1d(i)
      END DO
    END IF
    
    ! redistribute data into final decomposition
    DO i = 1, SIZE(io_data)
      CALL exchange_data(io_data(i)%redistrib_pattern, var_data(i)%DATA, &
        & local_buffer_2d)
    END DO
    
  END SUBROUTINE distrib_read_int_1d_multi_var
  
  !-------------------------------------------------------------------------
  
  SUBROUTINE distrib_read_real_1d_multi_var(ncid, var_name, local_buffer_1d, &
    & local_buffer_2d, var_data, io_data)
    
    INTEGER, INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: var_name
    REAL(wp), INTENT(inout) :: local_buffer_1d(:)
    REAL(wp), INTENT(inout) :: local_buffer_2d(:, :)
    TYPE(var_data_1d_wp), INTENT(inout) :: var_data(:)
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)
    
    INTEGER :: varid, i
    
    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data
    
    IF (SIZE(io_data) == 0) RETURN
    
    IF (SIZE(var_data) < SIZE(io_data)) &
      & CALL finish("distrib_read_real_1d_multi_var", "var_data too small")
    
    CALL check_basic_data_index(io_data(:))
    
    basic_io_data => basic_data(io_data(1)%basic_data_index)
    
    ! read data from file into io decomposition
    IF (basic_io_data%COUNT > 0) THEN
      CALL nf(nf_inq_varid(ncid, var_name, varid))
      
      ! only read io_decomp part
      CALL nf(nf_get_vara_double(ncid, varid, (/basic_io_data%start/), &
        & (/basic_io_data%COUNT/), local_buffer_1d(:)))
      
      DO i = 1, SIZE(local_buffer_1d(:))
        local_buffer_2d(idx_no(i),blk_no(i)) = local_buffer_1d(i)
      END DO
    END IF
    
    ! redistribute data into final decomposition
    DO i = 1, SIZE(io_data)
      CALL exchange_data(io_data(i)%redistrib_pattern, var_data(i)%DATA, &
        & local_buffer_2d)
    END DO
    
  END SUBROUTINE distrib_read_real_1d_multi_var
  
  !-------------------------------------------------------------------------
  
  SUBROUTINE distrib_read_int_2d_multi_var(ncid, var_name, local_buffer_2d, &
    & local_buffer_3d, var_data, dim_order, &
    & io_data)
    
    INTEGER, INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: var_name
    INTEGER, INTENT(inout) :: local_buffer_2d(:,:) ! (n io points, nlev)
    INTEGER, INTENT(inout) :: local_buffer_3d(:,:,:) ! if (dim_order == IDX_BLK_LVL)
    !   dimensions(nproma, n io blk, nlev)
    ! if (dim_order == IDX_LVL_BLK)
    !   dimensions(nproma, nlev, n io blk)
    TYPE(var_data_2d_int), INTENT(inout) :: var_data(:)
    INTEGER, INTENT(in) :: dim_order
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)
    
    INTEGER :: varid, i, j, nlev
    
    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data
    
    IF (SIZE(io_data) == 0) RETURN
    
    IF (SIZE(var_data) < SIZE(io_data)) &
      & CALL finish("distrib_read_int_2d_multi_var", "var_data too small")
    
    IF (dim_order /= idx_blk_lvl .AND. dim_order /= idx_lvl_blk) &
      & CALL finish("distrib_read_int_2d_multi_var", "invalid argument dim_order")
    
    CALL check_basic_data_index(io_data(:))
    
    basic_io_data => basic_data(io_data(1)%basic_data_index)
    
    nlev = SIZE(local_buffer_2d, 2)
    
    IF (basic_io_data%COUNT > 0) THEN
      
      IF ((dim_order == idx_blk_lvl) .AND. &
        & (nlev /= SIZE(local_buffer_3d, 3))) &
        & CALL finish("distrib_read_int_2d_multi_var", &
        & "nlev is not the same for all buffers(a)")
      IF ((dim_order == idx_lvl_blk) .AND. &
        & (nlev /= SIZE(local_buffer_3d, 2))) &
        & CALL finish("distrib_read_int_2d_multi_var", &
        & "nlev is not the same for all buffers(b)")
      
      CALL nf(nf_inq_varid(ncid, var_name, varid))
      
      ! only read io_decomp part
      CALL nf(nf_get_vara_int(ncid, varid, (/basic_io_data%start, 1/), &
        & (/basic_io_data%COUNT, nlev/), &
        & local_buffer_2d(:,:)))
    END IF
    
    IF (dim_order == idx_blk_lvl) THEN
      DO j = 1, nlev
        DO i = 1, basic_io_data%COUNT
          local_buffer_3d(idx_no(i),blk_no(i), j) = local_buffer_2d(i, j)
        END DO
      END DO
      
      DO i = 1, SIZE(io_data)
        DO j = 1, nlev
          CALL exchange_data(io_data(i)%redistrib_pattern, &
            & var_data(i)%DATA(:,:,j), &
            & local_buffer_3d(:,:,j))
        END DO
      END DO
    ELSE
      DO j = 1, nlev
        DO i = 1, basic_io_data%COUNT
          local_buffer_3d(idx_no(i), j, blk_no(i)) = local_buffer_2d(i, j)
        END DO
      END DO
      
      DO i = 1, SIZE(io_data)
        CALL exchange_data(io_data(i)%redistrib_pattern, &
          & var_data(i)%DATA(:,:,:), &
          & local_buffer_3d(:,:,:))
      END DO
    END IF
    
  END SUBROUTINE distrib_read_int_2d_multi_var
  
  !-------------------------------------------------------------------------
  
  SUBROUTINE distrib_read_real_2d_multi_var(ncid, var_name, local_buffer_2d, &
    & local_buffer_3d, var_data, dim_order, &
    & io_data)
    
    INTEGER, INTENT(in) :: ncid
    CHARACTER(LEN=*), INTENT(in) :: var_name
    REAL(wp), INTENT(inout) :: local_buffer_2d(:,:) ! (n io points, nlev)
    REAL(wp), INTENT(inout) :: local_buffer_3d(:,:,:) ! if (dim_order == IDX_BLK_LVL)
    !   dimensions(nproma, n io blk, nlev)
    ! if (dim_order == IDX_LVL_BLK)
    !   dimensions(nproma, nlev, n io blk)
    TYPE(var_data_2d_wp), INTENT(inout) :: var_data(:)
    INTEGER, INTENT(in) :: dim_order
    TYPE(t_distrib_read_data), INTENT(in) :: io_data(:)
    
    INTEGER :: varid, i, j, nlev
    
    TYPE(t_basic_distrib_read_data), POINTER :: basic_io_data
    
    IF (SIZE(io_data) == 0) RETURN
    
    IF (SIZE(var_data) < SIZE(io_data)) &
      & CALL finish("distrib_read_real_2d_multi_var", "var_data too small")
    
    IF (dim_order /= idx_blk_lvl .AND. dim_order /= idx_lvl_blk) &
      & CALL finish("distrib_read_real_2d_multi_var", "invalid argument dim_order")
    
    CALL check_basic_data_index(io_data(:))
    
    basic_io_data => basic_data(io_data(1)%basic_data_index)
    
    nlev = SIZE(local_buffer_2d, 2)
    
    IF (basic_io_data%COUNT > 0) THEN
      
      IF ((dim_order == idx_blk_lvl) .AND. &
        & (nlev /= SIZE(local_buffer_3d, 3))) &
        & CALL finish("distrib_read_real_2d_multi_var", &
        & "nlev is not the same for all buffers(a)")
      IF ((dim_order == idx_lvl_blk) .AND. &
        & (nlev /= SIZE(local_buffer_3d, 2))) &
        & CALL finish("distrib_read_real_2d_multi_var", &
        & "nlev is not the same for all buffers(b)")
      
      CALL nf(nf_inq_varid(ncid, var_name, varid))
      
      ! only read io_decomp part
      CALL nf(nf_get_vara_double(ncid, varid, (/basic_io_data%start, 1/), &
        & (/basic_io_data%COUNT, nlev/), &
        & local_buffer_2d(:,:)))
    END IF
    
    IF (dim_order == idx_blk_lvl) THEN
      DO j = 1, nlev
        DO i = 1, basic_io_data%COUNT
          local_buffer_3d(idx_no(i),blk_no(i), j) = local_buffer_2d(i, j)
        END DO
      END DO
      
      DO i = 1, SIZE(io_data)
        DO j = 1, nlev
          CALL exchange_data(io_data(i)%redistrib_pattern, &
            & var_data(i)%DATA(:,:,j), &
            & local_buffer_3d(:,:,j))
        END DO
      END DO
    ELSE
      DO j = 1, nlev
        DO i = 1, basic_io_data%COUNT
          local_buffer_3d(idx_no(i), j, blk_no(i)) = local_buffer_2d(i, j)
        END DO
      END DO
      
      DO i = 1, SIZE(io_data)
        CALL exchange_data(io_data(i)%redistrib_pattern, &
          & var_data(i)%DATA(:,:,:), &
          & local_buffer_3d(:,:,:))
      END DO
    END IF
    
  END SUBROUTINE distrib_read_real_2d_multi_var
  
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
        CALL message('mo_model_domain_import netCDF error', nf_strerror(STATUS), &
          & level=em_warn)
      ELSE
        CALL finish('mo_model_domain_import netCDF error', nf_strerror(STATUS))
      ENDIF
    ENDIF
    
  END SUBROUTINE nf
  
END MODULE mo_read_netcdf_distributed
