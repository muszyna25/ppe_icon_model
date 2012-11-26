! ----------------------------------------------------------------------
! Main subroutine for reading and remapping.
!
! Initial revision: F. Prill, DWD (2012-11-09)
! ----------------------------------------------------------------------
!
! Restrictions of this version:
!
! - Land-sea mask not taken into account for interpolation weights.
! - Current algorithm is not bit-reproducible (requires unique heap
!   ordering).
! - only 1st order weights for Gaussian->ICON interpolation
! - no internal conversion of Omega to W.
!
! Major TODOs
!
! - Documentation
!
! Minor TODOs:
!
! - don't define local TOL constants.
! - replace "WRITE(0,*)" by "CALL message(...)".
! - replace "tic"/"toc" by timer calls
! - define nproma as a namelist parameter
!
MODULE mo_remap

!$  USE OMP_LIB

  USE mo_kind,              ONLY: wp
  USE mo_parallel_config,   ONLY: nproma
  USE mo_communication,     ONLY: blk_no
  USE mo_timer,             ONLY: tic, toc
  USE mo_mpi,               ONLY: get_my_mpi_work_id, p_n_work
  USE mo_gnat_gridsearch,   ONLY: gnat_init_grid
  USE mo_remap_config,      ONLY: dbg_level, MAX_NAME_LENGTH
  USE mo_remap_sync,        ONLY: t_gather_c, allocate_gather_c,            &
    &                             gather_field2D_c, gather_field3D_c,       &
    &                             finalize_gather_c
  USE mo_remap_weights,     ONLY: prepare_interpolation, consistency_check
  USE mo_remap_grid,        ONLY: load_grid, finalize_grid
  USE mo_remap_shared,      ONLY: t_grid, GRID_TYPE_ICON
  USE mo_remap_intp,        ONLY: t_intp_data, allocate_intp_data,          &
    &                             deallocate_intp_data
  USE mo_remap_input,       ONLY: load_metadata_input, read_input_namelist, &
    &                             n_input_fields, input_field, n_zaxis,     &
    &                             zaxis_metadata, global_metadata,          &
    &                             input_import_data, close_input
  USE mo_remap_subdivision, ONLY: decompose_grid, create_grid_covering,     &
    &                             IMAX, IMIN, get_latitude_range
  USE mo_remap_output,      ONLY: load_metadata_output, open_output,        &
    &                             close_output, store_field, varID,         &
    &                             t_output_grid
  USE mo_remap_io,          ONLY: l_have3dbuffer, t_file_metadata,          &
    &                             open_file, close_file, in_filename,       &
    &                             out_filename, in_grid_filename,           &
    &                             out_grid_filename, in_type, out_type,     &
    &                             t_file_metadata, read_remap_namelist

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: remap_main

  INTEGER,  PARAMETER :: rank0     = 0 !< Input process rank
  INTEGER,  PARAMETER :: rank0_out = 0 !< Output process rank

CONTAINS

  ! Main subroutine for reading data from one file, remapping and writing
  ! to an output file.
  !
  ! @note All input arguments for this subroutine are read from
  !       a namelist file.
  !
  SUBROUTINE remap_main(namelist_filename, opt_input_cfg_filename)
    CHARACTER (LEN=*), INTENT(IN)           :: namelist_filename
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: opt_input_cfg_filename
    ! local variables:
    CHARACTER (len=MAX_NAME_LENGTH) :: input_cfg_filename
    TYPE (t_file_metadata) :: in_data, out_data, in_grid, out_grid
    TYPE (t_grid), TARGET  :: gridA_global, gridB_global,        &
      &                       gridA,        gridB,               &
      &                       gridA_cov,    gridB_cov
    TYPE (t_intp_data)     :: intp_data_A, intp_data_B
    TYPE(t_gather_c)       :: comm_data_c_A_cov, comm_data_c_B
    TYPE (t_output_grid)   :: output_grid

    REAL(wp), ALLOCATABLE  :: fieldB_loc(:,:,:), fieldB_glb(:,:,:)
    INTEGER                :: i, nthreads, ivar, max_nlev,       &
      &                       glb_max_nlev, glb_idx
    REAL(wp)               :: lat_rangeA(2), lat_rangeB(2),      &
      &                       lat_range(2)
    REAL                   :: time_s,                            &
      &                       time_comm,     time_write,         &
      &                       time_comm_tot, time_read_tot,      &
      &                       time_intp_tot, time_write_tot

    WRITE (0,*) "# Conservative Remapping"

    IF (PRESENT(opt_input_cfg_filename)) THEN
      input_cfg_filename = TRIM(opt_input_cfg_filename)
    ELSE
      input_cfg_filename = TRIM(namelist_filename)
    END IF

    nthreads = 1
    !$  nthreads = omp_get_max_threads()
    IF (get_my_mpi_work_id() == rank0) THEN
      WRITE (0,'(a,i4,a,i4,a)') "# running with ", p_n_work, " MPI processes and ", nthreads, " thread(s)."
    END IF

    ! ----------------------------------------------------------------------------
    ! load data
    ! ----------------------------------------------------------------------------

    CALL tic(time_s)  ! performance measurement: start
    ! read namelists
    CALL read_remap_namelist(namelist_filename)
    CALL read_input_namelist(input_cfg_filename, rank0)

    ! load variable meta-data
    CALL open_file(in_filename,       in_type,  in_data,  rank0)
    CALL open_file(in_grid_filename,  in_type,  in_grid,  rank0)
    CALL open_file(out_grid_filename, out_type, out_grid, rank0)
    CALL load_metadata_input (in_data, rank0)
    CALL load_metadata_output(out_grid, output_grid, rank0)

    ! load global grid to internal data structure
    IF (get_my_mpi_work_id() == rank0) THEN
      CALL open_output(out_filename, global_metadata, input_field(:), n_input_fields, &
        &              zaxis_metadata, n_zaxis, output_grid, out_data, out_grid)
      CALL load_grid(gridA_global, in_data%structure, rank0,  opt_file=in_grid)
      CALL load_grid(gridB_global, out_grid%structure, rank0, opt_file=out_grid)
    ELSE
      CALL load_grid(gridA_global, in_data%structure,  rank0)
      CALL load_grid(gridB_global, out_grid%structure, rank0)
    END IF

    ! performance measurement: stop
    IF (dbg_level >= 2)  WRITE (0,*) "# > loading grids: elapsed time: ", toc(time_s), " sec."
    CALL tic(time_s)  ! performance measurement: start

    ! grid subdivision, create local grid coverings:
    CALL decompose_grid(gridA_global, gridA)
    CALL decompose_grid(gridB_global, gridB)
    lat_rangeA = get_latitude_range(gridA)
    lat_rangeB = get_latitude_range(gridB)
    lat_range = (/ MIN(lat_rangeB(IMIN), lat_rangeA(IMIN)), &
      &            MAX(lat_rangeB(IMAX), lat_rangeA(IMAX)) /)
    CALL create_grid_covering(gridA_global, gridA, gridA_cov, lat_range)
    CALL create_grid_covering(gridB_global, gridB, gridB_cov, lat_range)

    ! free global grids, they are no longer needed:
    CALL finalize_grid(gridA_global, gridB_global)

    ! performance measurement: stop
    IF (dbg_level >= 2)  WRITE (0,*) "# > partitioning grids: elapsed time: ", toc(time_s), " sec."

    ! ----------------------------------------------------------------------------
    ! compute weights
    ! ----------------------------------------------------------------------------

    CALL tic(time_s)  ! performance measurement: start

    ! build fast search tree for proximity queries in unstructured triangular grid
    IF (gridA_cov%structure == GRID_TYPE_ICON) THEN
      IF (dbg_level >= 2) WRITE (0,*) "# init GNAT data structure"
      CALL gnat_init_grid(gridA_cov%p_patch, .TRUE., 1, gridA_cov%p_patch%nblks_c)
      IF (dbg_level >= 3) WRITE (0,*) "# done."
    END IF
    IF (gridB_cov%structure == GRID_TYPE_ICON) THEN
      IF (dbg_level >= 3) WRITE (0,*) "# init GNAT data structure"
      CALL gnat_init_grid(gridB_cov%p_patch, .TRUE., 1, gridB_cov%p_patch%nblks_c)
      IF (dbg_level >= 3) WRITE (0,*) "# done."
    END IF

    ! allocate interpolation weights
    CALL allocate_intp_data(intp_data_A, gridA)
    CALL allocate_intp_data(intp_data_B, gridB)

    ! compute interpolation weights
    CALL prepare_interpolation(gridA, gridB, gridA_cov, gridB_cov, &
      &                        intp_data_A, intp_data_B)

    ! performance measurement: stop
    IF (dbg_level >= 2)  WRITE (0,*) "# > weight computation: elapsed time: ", toc(time_s), " sec."

    ! perform some consistency checks:
    IF (dbg_level >= 2) THEN
      CALL consistency_check(gridA, intp_data_A, "list 1")
      CALL consistency_check(gridB, intp_data_B, "list 2")
    END IF

    ! ----------------------------------------------------------------------------
    ! interpolation + output
    ! ----------------------------------------------------------------------------

    CALL tic(time_s)  ! performance measurement: start

    ! for distributed computation: create communication data
    CALL allocate_gather_c(gridA_cov, rank0,     comm_data_c_A_cov)
    CALL allocate_gather_c(gridB,     rank0_out, comm_data_c_B)

    max_nlev     = MAXVAL(input_field(1:n_input_fields)%nlev)
    glb_max_nlev = max_nlev
    IF (.NOT. l_have3dbuffer) glb_max_nlev = 1
    ALLOCATE(fieldB_loc(nproma, max_nlev, gridB%p_patch%nblks_c), &
      &      fieldB_glb(nproma, glb_max_nlev, blk_no(gridB%p_patch%n_patch_cells_g)))

    time_comm_tot  = 0.
    time_read_tot  = 0.
    time_intp_tot  = 0.
    time_write_tot = 0.
    DO ivar=1,n_input_fields
      IF (dbg_level >= 1) &
        &  WRITE (0,*) "# read and interpolate variable '", TRIM(input_field(ivar)%inputname), "'"

      ! load input fields, e.g. IFS GRIB data, interpolate onto ICON grid
      CALL input_import_data(in_data, ivar, comm_data_c_A_cov,          &
        &                    fieldB_loc, gridA_cov, gridB, intp_data_B, &
        &                    time_comm_tot, time_read_tot, time_intp_tot)

      IF (l_have3dbuffer) THEN
        ! for distributed computation: gather interpolation result on rank 0
        CALL tic(time_comm)  ! performance measurement: start
        IF (input_field(ivar)%nlev == 1) THEN
          CALL gather_field2D_c(comm_data_c_B, fieldB_loc(:,1,:), fieldB_glb(:,1,:))
        ELSE
          CALL gather_field3D_c(comm_data_c_B, max_nlev, fieldB_loc(:,:,:), fieldB_glb)
        END IF
        time_comm_tot = time_comm_tot + toc(time_comm)
      END IF

      ! output on PE "rank0"
      DO i=1,input_field(ivar)%nlev
        glb_idx = i
        IF (.NOT. l_have3dbuffer) THEN
          glb_idx = 1
          CALL tic(time_comm)  ! performance measurement: start
          CALL gather_field2D_c(comm_data_c_B, fieldB_loc(:,i,:), fieldB_glb(:,glb_idx,:))
          time_comm_tot = time_comm_tot + toc(time_comm)
        END IF

        IF (get_my_mpi_work_id() == rank0) THEN
          CALL tic(time_write)  ! performance measurement: start
          CALL store_field(out_data, varID(ivar), i, fieldB_glb(:,glb_idx,:))
          time_write_tot = time_write_tot + toc(time_write)
        END IF
      END DO
    END DO ! ivar

    ! performance measurement: stop
    IF (dbg_level >= 2)  WRITE (0,*) "# > interpolation and output: elapsed time: ", toc(time_s), " sec."
    IF (dbg_level >= 1) THEN
      WRITE (0,*) "#     comm:  ", time_comm_tot,  " sec."
      WRITE (0,*) "#     read:  ", time_read_tot,  " sec."
      WRITE (0,*) "#     intp:  ", time_intp_tot,  " sec."
      WRITE (0,*) "#     write: ", time_write_tot, " sec."
    END IF

    ! ----------------------------------------------------------------------------
    ! clean up
    ! ----------------------------------------------------------------------------

    CALL deallocate_intp_data(intp_data_A)
    CALL deallocate_intp_data(intp_data_B)
    CALL finalize_gather_c(comm_data_c_A_cov, comm_data_c_B)
    CALL finalize_grid(gridA, gridB, gridA_cov, gridB_cov)

    DEALLOCATE(fieldB_loc, fieldB_glb)

    IF (get_my_mpi_work_id() == rank0) THEN
      CALL close_input()
      CALL close_output(out_data, output_grid)
      CALL close_file(in_data)
      CALL close_file(in_grid)
      CALL close_file(out_grid)
    END IF
    WRITE (0,*) "# Done."

  END SUBROUTINE remap_main

END MODULE mo_remap
