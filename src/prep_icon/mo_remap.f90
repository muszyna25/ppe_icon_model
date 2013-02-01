! ----------------------------------------------------------------------
! Main subroutine for remapping data between different grids.
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
! - RBF interpolation only for "u","v" (REGULAR) -> "vn"/"vt" (ICON)
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
  USE mo_impl_constants,    ONLY: SUCCESS
  USE mo_communication,     ONLY: blk_no
  USE mo_exception,         ONLY: finish
  USE mo_timer,             ONLY: tic, toc
  USE mo_mpi,               ONLY: get_my_mpi_work_id, p_n_work
  USE mo_grib2,             ONLY: t_grib2_var
  USE mo_cf_convention,     ONLY: t_cf_var
  USE mo_gnat_gridsearch,   ONLY: gnat_init_grid
  USE mo_remap_config,      ONLY: dbg_level, MAX_NSTENCIL_CONS,                &
    &                             MAX_NSTENCIL_RBF, MIN_NFOREIGN,              &
    &                             MAX_NAME_LENGTH                         
  USE mo_remap_sync,        ONLY: t_gather, allocate_gather_c,                 &
    &                             allocate_gather_e,                           &
    &                             gather_field2D, gather_field3D,              &
    &                             finalize_gather                         
  USE mo_remap_weights_cons,ONLY: prepare_interpolation_cons,                  &
    &                             consistency_check
  USE mo_remap_weights_rbf, ONLY: prepare_interpolation_rbf_vec                
  USE mo_remap_intp,        ONLY: t_intp_data, allocate_intp_data,             &
    &                             finalize_intp_data, mask_intp_coeffs,        &
    &                             copy_intp_data_sthreaded,                    &
    &                             interpolate_e_2D_rbf
  USE mo_remap_io,          ONLY: l_have3dbuffer, t_file_metadata,             &
    &                             open_file, in_filename,                      &
    &                             out_filename, in_grid_filename,              &
    &                             out_grid_filename, in_type, out_type,        &
    &                             read_remap_namelist,                         &
    &                             lcompute_vn, lcompute_vt, rbf_vec_scale,     &
    &                             close_file
  USE mo_remap_input,       ONLY: load_metadata_input, read_input_namelist,    &
    &                             n_input_fields, input_field, n_zaxis,        &
    &                             zaxis_metadata, global_metadata,             &
    &                             input_import_data, close_input,              &
    &                             generate_missval_mask, INTP_CONS, INTP_NONE, &
    &                             INTP_RBF, EDGE_GRID, CONST_UNINITIALIZED,    &
    &                             CONST_FALSE, CONST_TRUE, field_id_u,         &
    &                             field_id_v, t_field_adjustment, CELL_GRID,   &
    &                             EDGE_GRID
  USE mo_remap_output,      ONLY: load_metadata_output, open_output,           &
    &                             close_output, store_field, varID,            &
    &                             t_output_grid
  USE mo_remap_grid,        ONLY: load_grid, finalize_grid
  USE mo_remap_shared,      ONLY: t_grid, GRID_TYPE_ICON, GRID_TYPE_REGULAR
  USE mo_remap_hydcorr,     ONLY: remap_horz
  USE mo_remap_subdivision, ONLY: decompose_grid, create_grid_covering,     &
    &                             IMAX, IMIN, get_latitude_range



  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PRIVATE
  PUBLIC :: remap_main

  INTEGER,  PARAMETER :: rank0        = 0        !< Input process rank
  INTEGER,  PARAMETER :: rank0_out    = 0        !< Output process rank

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_remap'

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
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::remap_main'

    TYPE (t_file_metadata)  :: in_data, out_data, in_grid, out_grid
    TYPE (t_grid), TARGET   :: gridA_global, gridB_global,           &
      &                        gridA,        gridB,                  &
      &                        gridA_cov,    gridB_cov
    TYPE (t_intp_data)      :: intp_data_A, intp_data_B,             &
      &                        intp_data_rbf_vn_u_B,                 &
      &                        intp_data_rbf_vn_v_B,                 &
      &                        intp_data_rbf_vt_u_B,                 &
      &                        intp_data_rbf_vt_v_B
    TYPE(t_gather), TARGET  :: comm_data_c_A_cov, comm_data_c_B,     &
      &                        comm_data_e_B
    TYPE(t_gather), POINTER :: comm_data_B
    TYPE (t_output_grid)    :: output_grid
                            
    REAL(wp), ALLOCATABLE   :: fieldB_loc(:,:,:), fieldB_glb(:,:,:), &
      &                        tmp_rfield3D(:,:,:)
    INTEGER                 :: i, nthreads, ivar, max_nlev,          &
      &                        glb_max_nlev, glb_idx, ilev,          &
      &                        fac, max_nforeignA, max_nforeignB,    &
      &                        max_nforeign, igather_size, ierrstat, &
      &                        ihorizontal_size_loc,                 &
      &                        ihorizontal_size_glb
    REAL(wp)                :: lat_rangeA(2), lat_rangeB(2),         &
      &                        lat_range(2)
    REAL                    :: time_s, time_comm, time_write,        &
      &                        time_comm_tot, time_read_tot,         &
      &                        time_intp_tot, time_write_tot,        &
      &                        time_intp

    WRITE (0,*) "# Conservative Remapping"

    IF (PRESENT(opt_input_cfg_filename)) THEN
      input_cfg_filename = TRIM(opt_input_cfg_filename)
    ELSE
      input_cfg_filename = TRIM(namelist_filename)
    END IF

    nthreads = 1
!$  nthreads = omp_get_max_threads()
    IF (get_my_mpi_work_id() == rank0) THEN
      WRITE (0,'(a,i4,a,i4,a)') " # running with ", p_n_work, " MPI processes and ", nthreads, " thread(s)."
    END IF

    ! --------------------------------------------------------------------------
    ! load data
    ! --------------------------------------------------------------------------

    CALL tic(time_s)  ! performance measurement: start
    ! read namelists
    CALL read_remap_namelist(namelist_filename)
    CALL read_input_namelist(input_cfg_filename, rank0, lcompute_vn, lcompute_vt)

    ! load variable meta-data
    CALL open_file(in_filename,       in_type,  in_data,  rank0)
    CALL open_file(in_grid_filename,  in_type,  in_grid,  rank0)
    CALL open_file(out_grid_filename, out_type, out_grid, rank0)
    CALL load_metadata_input (in_data, rank0)
    CALL load_metadata_output(out_grid, output_grid, rank0)

    ! special case of "u","v" -> "vn":
    IF ((field_id_u /= CONST_UNINITIALIZED)  .AND.   &
      & (field_id_v /= CONST_UNINITIALIZED)) THEN
      IF (lcompute_vn) THEN
        CALL insert_v_edge( "vn", &
          &      t_cf_var('normal_velocity', 'm s-1', 'velocity normal to edge', DATATYPE_FLT32), &
          &      t_grib2_var(0, 2, 34, DATATYPE_PACK16, GRID_REFERENCE, EDGE_GRID) )
      END IF
      IF (lcompute_vt) THEN
        CALL insert_v_edge( "vt", &
          &      t_cf_var('tangential_wind', 'm s-1', 'tangential-component of wind', DATATYPE_FLT32), &
          &      t_grib2_var(0, 2, 35, DATATYPE_PACK16, GRID_REFERENCE, EDGE_GRID) )
      END IF
    END IF

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

    ! estimate the number of cells of a covering which are foreign to
    ! the local process (needed to set the dimension of a list which
    ! will be communicated):
    fac = 10 ! safety factor (>= 1)
    max_nforeignA = NINT(fac * REAL(gridB_cov%p_patch%n_patch_cells,wp)) * &
      &   ((lat_range(IMAX) - lat_rangeA(IMAX)) + (lat_rangeA(IMIN) - lat_range(IMIN))) / &
      &   (lat_range(IMAX) - lat_range(IMIN))
    max_nforeignB = NINT(fac * REAL(gridA_cov%p_patch%n_patch_cells,wp)) * &
      &   ((lat_range(IMAX) - lat_rangeB(IMAX)) + (lat_rangeB(IMIN) - lat_range(IMIN))) / &
      &   (lat_range(IMAX) - lat_range(IMIN))
    max_nforeign = MAX(MIN_NFOREIGN, MAX(max_nforeignA, max_nforeignB))

    ! performance measurement: stop
    IF (dbg_level >= 2)  THEN
      WRITE (0,*) "# > partitioning grids: elapsed time: ", toc(time_s), " sec."
      WRITE (0,*) "# allocating ", max_nforeign, " foreign weights."
    END IF

    ! --------------------------------------------------------------------------
    ! compute weights
    ! --------------------------------------------------------------------------

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
    CALL allocate_intp_data(intp_data_A,          gridA%p_patch%nblks_c, MAX_NSTENCIL_CONS)
    CALL allocate_intp_data(intp_data_B,          gridB%p_patch%nblks_c, MAX_NSTENCIL_CONS)

    ! compute interpolation weights: RBF remapping
    IF (ANY(input_field(:)%intp_method == INTP_RBF)) THEN
      IF (lcompute_vn) THEN
        CALL allocate_intp_data(intp_data_rbf_vn_u_B, gridB%p_patch%nblks_e, MAX_NSTENCIL_RBF, opt_smaxsize=1)
        CALL allocate_intp_data(intp_data_rbf_vn_v_B, gridB%p_patch%nblks_e, MAX_NSTENCIL_RBF, opt_smaxsize=1)
        CALL prepare_interpolation_rbf_vec(gridA_cov, gridB, intp_data_rbf_vn_u_B, intp_data_rbf_vn_v_B, &
          &                                .FALSE., rbf_vec_scale)
      END IF
      IF (lcompute_vt) THEN
        CALL allocate_intp_data(intp_data_rbf_vt_u_B, gridB%p_patch%nblks_e, MAX_NSTENCIL_RBF, opt_smaxsize=1)
        CALL allocate_intp_data(intp_data_rbf_vt_v_B, gridB%p_patch%nblks_e, MAX_NSTENCIL_RBF, opt_smaxsize=1)
        CALL prepare_interpolation_rbf_vec(gridA_cov, gridB, intp_data_rbf_vt_u_B, intp_data_rbf_vt_v_B, &
          &                                .TRUE., rbf_vec_scale)
      END IF
    END IF ! if RBF

    ! compute interpolation weights: Conservative remapping
    CALL prepare_interpolation_cons(gridA, gridB, gridA_cov, gridB_cov, &
      &                             intp_data_A, intp_data_B, max_nforeign)

    ! performance measurement: stop
    IF (dbg_level >= 2)  WRITE (0,*) "# > weight computation: elapsed time: ", toc(time_s), " sec."

    ! perform some consistency checks:
    IF (dbg_level >= 2) THEN
      CALL consistency_check(gridA, intp_data_A, "list 1")
      CALL consistency_check(gridB, intp_data_B, "list 2")
    END IF

    ! --------------------------------------------------------------------------
    ! interpolation + output
    ! --------------------------------------------------------------------------

    CALL tic(time_s)  ! performance measurement: start

    ! for distributed computation: create communication data
    CALL allocate_gather_c(gridA_cov, rank0,     comm_data_c_A_cov)
    CALL allocate_gather_c(gridB,     rank0_out, comm_data_c_B)
    CALL allocate_gather_e(gridB,     rank0_out, comm_data_e_B)
    comm_data_B => NULL()

    max_nlev     = MAXVAL(input_field(1:n_input_fields)%nlev)
    glb_max_nlev = max_nlev
    IF (.NOT. l_have3dbuffer) glb_max_nlev = 1

    ! allocate temporary, local field for cell grids:
    ALLOCATE(tmp_rfield3D(nproma, max_nlev, gridA_cov%p_patch%nblks_c), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    time_comm_tot  = 0.
    time_read_tot  = 0.
    time_intp_tot  = 0.
    time_write_tot = 0.
    VAR_LOOP : DO ivar=1,n_input_fields

      IF (dbg_level >= 1)  &
        &  WRITE (0,*) "# read and interpolate variable '", TRIM(input_field(ivar)%inputname), "'"

      ! initialize temporary variables (if necessary):
      ihorizontal_size_loc = 0
      ihorizontal_size_glb = 0
      SELECT CASE(input_field(ivar)%grid_type)
      CASE (CELL_GRID) 
        IF (dbg_level >= 3)  WRITE (0,*) "# cell based field."
        ihorizontal_size_loc = gridB%p_patch%nblks_c
        ihorizontal_size_glb = gridB%p_patch%n_patch_cells_g
        comm_data_B => comm_data_c_B
      CASE (EDGE_GRID)
        IF (dbg_level >= 3)  WRITE (0,*) "# edge based field."
        ihorizontal_size_loc = gridB%p_patch%nblks_e
        ihorizontal_size_glb = gridB%p_patch%n_patch_edges_g
        comm_data_B => comm_data_e_B
      CASE DEFAULT
        CALL finish(routine, "Internal error!")
      END SELECT
      ! deallocate old temporary fields, if it does not match the
      ! required size:
      IF (ALLOCATED(fieldB_loc)) THEN
        IF (ihorizontal_size_loc /= SIZE(fieldB_loc,3)) THEN
          DEALLOCATE(fieldB_loc, STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
        END IF
      END IF
      IF (ALLOCATED(fieldB_glb)) THEN
        IF (ihorizontal_size_glb /= SIZE(fieldB_glb,3)) THEN
          DEALLOCATE(fieldB_glb, STAT=ierrstat)
          IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
        END IF
      END IF
      ! allocate temporary fields:
      IF (.NOT. ALLOCATED(fieldB_loc)) THEN
        ALLOCATE(fieldB_loc(nproma, max_nlev, ihorizontal_size_loc), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      END IF
      IF (.NOT. ALLOCATED(fieldB_glb)) THEN
        ALLOCATE(fieldB_glb(nproma, glb_max_nlev, blk_no(ihorizontal_size_glb)), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      END IF

      ! distinguish between the different interpolation methods:
      SELECT CASE(input_field(ivar)%intp_method)
      CASE(INTP_CONS)

        ! load input field, e.g. IFS GRIB data
        CALL input_import_data(in_data, ivar, comm_data_c_A_cov, &
          &                    tmp_rfield3D, gridA_cov,          &
          &                    time_comm_tot, time_read_tot      )

        ! interpolate onto other (e.g. ICON) grid

        ! perform horizontal interpolation
        CALL tic(time_intp)  ! performance measurement: start
        IF (input_field(ivar)%has_missvals == CONST_FALSE) THEN

          DO ilev=1,input_field(ivar)%nlev
            CALL remap_horz(in_data,                           &
              tmp_rfield3D(:,ilev,:), fieldB_loc  (:,ilev,:),  &
              gridA_cov, gridB, input_field(ivar)%fa,          &
              intp_data_B, comm_data_c_A_cov)
          END DO ! ilev

        ELSE

          ! missing values present:
          CALL handle_missing_values(in_data, fieldB_loc, gridA_cov, gridB,   &
            &                        intp_data_B, input_field(ivar)%fa,       &
            &                        input_field(ivar)%missval, tmp_rfield3D, &
            &                        input_field(ivar)%nlev, comm_data_c_A_cov)

        END IF
        time_intp_tot = time_intp_tot + toc(time_intp)

      CASE (INTP_RBF)

        IF (TRIM(input_field(ivar)%inputname) == "vn") THEN
          CALL interpolate_v_edge(in_data, comm_data_c_A_cov, gridA_cov, gridB,              &
            &                     intp_data_rbf_vn_u_B, intp_data_rbf_vn_v_B,                &
            &                     fieldB_loc, time_comm_tot, time_read_tot, time_intp_tot)
        END IF
        IF (TRIM(input_field(ivar)%inputname) == "vt") THEN
          CALL interpolate_v_edge(in_data, comm_data_c_A_cov, gridA_cov, gridB,              &
            &                     intp_data_rbf_vt_u_B, intp_data_rbf_vt_v_B,                &
            &                     fieldB_loc, time_comm_tot, time_read_tot, time_intp_tot)
        END IF

      CASE (INTP_NONE)
        IF (dbg_level >= 2)  WRITE (0,*) "# skipping"
        CYCLE VAR_LOOP
      CASE DEFAULT
        CALL finish(routine, "Internal error!")

      END SELECT ! intp_method

      !----------------------------
      ! write result to output file
      !----------------------------

      IF (l_have3dbuffer) THEN
        ! for distributed computation: gather interpolation result on rank 0
        CALL tic(time_comm)  ! performance measurement: start
        igather_size = max_nlev
        IF (input_field(ivar)%nlev == 1) igather_size = 1
        CALL gather_field3D(comm_data_B, igather_size, fieldB_loc, fieldB_glb)
        time_comm_tot = time_comm_tot + toc(time_comm)
      END IF

      ! output on PE "rank0"
      DO i=1,input_field(ivar)%nlev
        glb_idx = i
        IF (.NOT. l_have3dbuffer) THEN
          glb_idx = 1
          CALL tic(time_comm)  ! performance measurement: start
          CALL gather_field2D(comm_data_B, fieldB_loc(:,i,:), fieldB_glb(:,glb_idx,:))
          time_comm_tot = time_comm_tot + toc(time_comm)
        END IF

        IF (get_my_mpi_work_id() == rank0) THEN
          CALL tic(time_write)  ! performance measurement: start
          CALL store_field(out_data, varID(ivar), i, fieldB_glb(:,glb_idx,:))
          time_write_tot = time_write_tot + toc(time_write)
        END IF
      END DO
    END DO VAR_LOOP ! ivar
    DEALLOCATE(tmp_rfield3D, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

    ! performance measurement: stop
    IF (dbg_level >= 2) THEN
      WRITE (0,*) "# > interpolation and output: elapsed time: ", toc(time_s), " sec."
      WRITE (0,*) "#     comm:  ", time_comm_tot,  " sec."
      WRITE (0,*) "#     read:  ", time_read_tot,  " sec."
      WRITE (0,*) "#     intp:  ", time_intp_tot,  " sec."
      WRITE (0,*) "#     write: ", time_write_tot, " sec."
    END IF

    ! --------------------------------------------------------------------------
    ! clean up
    ! --------------------------------------------------------------------------

    CALL finalize_intp_data(intp_data_A, intp_data_B)
    IF (lcompute_vn) THEN
      CALL finalize_intp_data(intp_data_rbf_vn_u_B, intp_data_rbf_vn_v_B)
    END IF
    IF (lcompute_vt) THEN
      CALL finalize_intp_data(intp_data_rbf_vt_u_B, intp_data_rbf_vt_v_B)
    END IF
    CALL finalize_gather(comm_data_c_A_cov, comm_data_c_B, comm_data_e_B)
    CALL finalize_grid(gridA, gridB, gridA_cov, gridB_cov)

    DEALLOCATE(fieldB_loc, fieldB_glb, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

    IF (get_my_mpi_work_id() == rank0) THEN
      CALL close_input()
      CALL close_output(out_data, output_grid)
      CALL close_file(in_data, in_grid, out_grid)
    END IF
    WRITE (0,*) "# Done."
  END SUBROUTINE remap_main


  ! -------------------------------------------------------------------------------


  !> Adds a new, edge-based variable "vn" (velocity component normal
  !  to edges) or "vt" (tangential wind component) to the field
  !  list. This "vn"/"vt" is not read from the input file but serves
  !  only as the output field from the interpolation of "u", "v".
  !
  SUBROUTINE insert_v_edge(name, cf_desc, grib2_desc)
    CHARACTER (LEN=*), INTENT(IN) :: name        !< variable name (e.g. "vn" or "vt")
    TYPE(t_cf_var)   , INTENT(IN) :: cf_desc     !< NetCDF variable descriptor
    TYPE(t_grib2_var), INTENT(IN) :: grib2_desc  !< GRIB2 variable descriptor
    ! internal variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::insert_v_edge')
   
    n_input_fields = n_input_fields + 1
    input_field(n_input_fields)%inputname       = TRIM(name) ! unused
    input_field(n_input_fields)%outputname      = TRIM(name)
    input_field(n_input_fields)%code            = 0    ! unused
    input_field(n_input_fields)%type_of_layer   = "-"  ! unused
    input_field(n_input_fields)%grid_type       = EDGE_GRID
    input_field(n_input_fields)%cf              = cf_desc
    input_field(n_input_fields)%grib2           = grib2_desc
    input_field(n_input_fields)%has_missvals    = CONST_FALSE
    input_field(n_input_fields)%intp_method     = INTP_RBF
    input_field(n_input_fields)%varID           = -1   ! unused
    input_field(n_input_fields)%fa%lhydrostatic_correction  = .FALSE.

    ! we copy the some information from the "u" velocity component:
    IF (field_id_u == CONST_UNINITIALIZED) CALL finish(routine, "Internal error!")
    IF ((input_field(field_id_u)%zaxisID == CDI_UNDEFID) .AND.   &
      & (get_my_mpi_work_id() == rank0)) THEN
      CALL finish(routine, "Internal error (zaxis)!")
    END IF
    input_field(n_input_fields)%nlev            = input_field(field_id_u)%nlev
    input_field(n_input_fields)%steptype        = input_field(field_id_u)%steptype
    input_field(n_input_fields)%zaxisID         = input_field(field_id_u)%zaxisID 
    input_field(n_input_fields)%missval         = input_field(field_id_u)%missval
  END SUBROUTINE insert_v_edge


  ! -------------------------------------------------------------------------------


  !> Build a temporary set of interpolation weights masking entries
  !  with missing values.
  !
  SUBROUTINE handle_missing_values(in_data, fieldB_loc, gridA_cov, gridB,  &
    &                              intp_data_B, fa, missval, tmp_rfield3D, &
    &                              nlev, comm_data_c_A_cov)
    TYPE (t_file_metadata),    INTENT(IN)     :: in_data
    REAL(wp),                  INTENT(INOUT)  :: fieldB_loc(:,:,:)
    TYPE (t_grid),             INTENT(IN)     :: gridB, gridA_cov
    TYPE (t_intp_data),        INTENT(IN)     :: intp_data_B
    TYPE (t_field_adjustment), INTENT(IN)     :: fa
    REAL(wp),                  INTENT(IN)     :: missval
    REAL(wp),                  INTENT(INOUT)  :: tmp_rfield3D(:,:,:)
    INTEGER,                   INTENT(IN)     :: nlev
    TYPE(t_gather),            INTENT(IN)     :: comm_data_c_A_cov
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::handle_missing_values'
    LOGICAL, ALLOCATABLE   :: missval_mask2D(:,:), clear_mask2D(:,:)
    INTEGER                :: ierrstat, ilev
    TYPE (t_intp_data)     :: tmp_intp_data

    ! generate a LOGICAL-array mask for missing values based on
    ! the value on level 1:
    ALLOCATE(missval_mask2D(nproma,gridA_cov%p_patch%nblks_c), &
      &      clear_mask2D(nproma,gridB%p_patch%nblks_c), STAT=ierrstat)
    clear_mask2D(:,:) = .FALSE.
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    CALL generate_missval_mask(tmp_rfield3D(:,1,:), missval, missval_mask2D, &
      &                        gridA_cov%p_patch%nblks_c, gridA_cov%p_patch%npromz_c)

    ! copy data structure containing the interpolation
    ! coefficients:
    CALL copy_intp_data_sthreaded(intp_data_B, tmp_intp_data)
    ! modify the interpolation coefficients according to the
    ! missing value mask:
    CALL mask_intp_coeffs(missval_mask2D, gridB, clear_mask2D, tmp_intp_data)

    DO ilev=1,nlev
      CALL remap_horz(in_data,                                  &
        tmp_rfield3D(:,ilev,:), fieldB_loc(:,ilev,:),           &
        gridA_cov, gridB, fa, tmp_intp_data, comm_data_c_A_cov)
      ! clear entries without contributing, non-missing values:
      WHERE (clear_mask2D(:,:)) fieldB_loc(:,ilev,:) = missval
    END DO ! ilev

    ! clean up
    DEALLOCATE(missval_mask2D, clear_mask2D, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    CALL finalize_intp_data(tmp_intp_data)

  END SUBROUTINE handle_missing_values


  ! -------------------------------------------------------------------------------


  !> RBF interpolation for "u","v" (REGULAR) -> "vn" or "vt" (ICON)
  !
  SUBROUTINE interpolate_v_edge(in_data, comm_data_c_A_cov, gridA_cov, gridB,          &
    &                           intp_data_rbf_v_u_B, intp_data_rbf_v_v_B,              &
    &                           fieldB_loc, time_comm_tot, time_read_tot, time_intp_tot)

    TYPE (t_file_metadata), INTENT(IN) :: in_data                !< input file
    TYPE(t_gather)        , INTENT(IN) :: comm_data_c_A_cov      !< communication pattern for source partition
    TYPE (t_grid)         , INTENT(IN) :: gridA_cov, gridB       !< source and destination grid
    TYPE (t_intp_data)    , INTENT(IN) :: intp_data_rbf_v_u_B, & !< interpolation coefficients (u comp.)
      &                                   intp_data_rbf_v_v_B    !< interpolation coefficients (v comp.)
    REAL(wp)              , INTENT(INOUT) :: fieldB_loc(:,:,:)   !< output result (local)
    REAL                  , INTENT(INOUT) :: time_comm_tot, &    !< performance timer: communication
      &                                      time_read_tot, &    !< performance timer: file read
      &                                      time_intp_tot       !< performance timer: interpolation
    ! local variables:
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::interpolate_v_edge'
    INTEGER :: ilev, ierrstat, nlev
    REAL    :: time_intp
    REAL(wp), ALLOCATABLE :: tmp_rfield3D_u(:,:,:), tmp_rfield3D_v(:,:,:)

    IF ((gridA_cov%structure /= GRID_TYPE_REGULAR) .OR.  &
      & (gridB%structure /= GRID_TYPE_ICON)) THEN
      CALL finish(routine, "Not implemented!")
    END IF

    ! --- allocate temporary arrays:
    nlev = input_field(field_id_u)%nlev
    IF (nlev /= input_field(field_id_v)%nlev) &
      & CALL finish(routine, "Inconsistent input data!")
    
    ALLOCATE(tmp_rfield3D_u(nproma, nlev, gridA_cov%p_patch%nblks_c), &
      &      tmp_rfield3D_v(nproma, nlev, gridA_cov%p_patch%nblks_c), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! --- load input field, e.g. IFS GRIB data

    ! u component
    CALL input_import_data(in_data, field_id_u, comm_data_c_A_cov, &
      &                    tmp_rfield3D_u, gridA_cov,              &
      &                    time_comm_tot, time_read_tot      )
    ! v component
    CALL input_import_data(in_data, field_id_v, comm_data_c_A_cov, &
      &                    tmp_rfield3D_v, gridA_cov,              &
      &                    time_comm_tot, time_read_tot      )
    
    ! perform horizontal interpolation onto ICON grid
    CALL tic(time_intp)  ! performance measurement: start
    DO ilev=1,nlev
      CALL interpolate_e_2D_rbf(tmp_rfield3D_u(:,ilev,:), tmp_rfield3D_v(:,ilev,:), &
        &                       fieldB_loc(:,ilev,:),                               &
        &                       gridB, intp_data_rbf_v_u_B, intp_data_rbf_v_v_B)
    END DO ! ilev
    time_intp_tot = time_intp_tot + toc(time_intp)

    ! --- clean up:
    DEALLOCATE(tmp_rfield3D_u, tmp_rfield3D_v, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE interpolate_v_edge

END MODULE mo_remap
