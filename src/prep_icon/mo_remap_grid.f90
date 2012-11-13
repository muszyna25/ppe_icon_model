MODULE mo_remap_grid

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_math_utilities,     ONLY: t_geographical_coordinates
  USE mo_gnat_gridsearch,    ONLY: gnat_tree, UNASSOCIATED
  USE mo_remap_config,       ONLY: dbg_level
  USE mo_remap_shared,       ONLY: t_grid,                                  &
    &                              GRID_TYPE_ICON, GRID_TYPE_REGULAR,       &
    &                              get_containing_cell_generic,             &
    &                              backtransform_lambert_azimuthal,         &
    &                              get_containing_cell_gnat,                &
    &                              LIST_DEFAULT
  USE mo_remap_grid_icon,    ONLY: load_icon_grid
  USE mo_remap_grid_regular, ONLY: load_gaussian_grid,                      &
    &                              get_containing_cell_gauss => get_containing_cell
  USE mo_remap_io,           ONLY: t_file_metadata
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: load_grid
  PUBLIC :: finalize_grid
  PUBLIC :: get_containing_cell

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_grid')  

CONTAINS


  SUBROUTINE load_grid(grid, istructure, rank0, opt_file)
    TYPE (t_grid),    INTENT(INOUT) :: grid
    INTEGER,          INTENT(IN)    :: istructure
    INTEGER,          INTENT(IN)    :: rank0    !< MPI rank where file is actually read
    TYPE (t_file_metadata), INTENT(IN), OPTIONAL :: opt_file
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::load_grid')
    
    grid%structure = istructure
    SELECT CASE(grid%structure) 
    CASE (GRID_TYPE_ICON)
      CALL load_icon_grid(grid, rank0, opt_file)
    CASE (GRID_TYPE_REGULAR)
      CALL load_gaussian_grid(grid, rank0, opt_file)
    CASE DEFAULT
      CALL finish(routine, "Unknown grid type")
    END SELECT
  END SUBROUTINE load_grid


  !> Free one or more grid data structures
  RECURSIVE SUBROUTINE finalize_grid(grid, opt_grid2, opt_grid3, opt_grid4)
    TYPE (t_grid), INTENT(INOUT) :: grid
    TYPE (t_grid), INTENT(INOUT), OPTIONAL :: opt_grid2, opt_grid3, opt_grid4
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::finalize_grid')
    INTEGER :: ierrstat

    ! the following recursion makes this more readable on the caller
    ! side:
    IF (PRESENT(opt_grid4)) THEN
      CALL finalize_grid(opt_grid4)
      CALL finalize_grid(grid, opt_grid2, opt_grid3)
      RETURN
    ELSE IF (PRESENT(opt_grid3)) THEN
      CALL finalize_grid(opt_grid3, opt_grid4)
      CALL finalize_grid(grid, opt_grid2)
      RETURN
    ELSE IF (PRESENT(opt_grid2)) THEN
      CALL finalize_grid(opt_grid2)
      CALL finalize_grid(grid)
      RETURN
    END IF

    IF (dbg_level >= 10) WRITE (0,*) "# finalize grid '", TRIM(grid%name), "'"
    IF ((grid%structure /= GRID_TYPE_ICON) .AND. &
      & (grid%structure /= GRID_TYPE_REGULAR)) THEN
      CALL finish(routine, "Unknown grid type")
    END IF
    
    IF (grid%structure == GRID_TYPE_REGULAR) THEN
      DEALLOCATE(grid%regular_grid%xvals1D, grid%regular_grid%yvals1D,   &
        &        grid%regular_grid%edge_lon, grid%regular_grid%edge_lat, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed (Gaussian)!")
    END IF
    IF (grid%structure == GRID_TYPE_ICON) THEN
      DEALLOCATE(grid%p_patch%verts%cell_idx,   grid%p_patch%verts%cell_blk,       &
        &        grid%p_patch%cells%neighbor_idx, grid%p_patch%cells%neighbor_blk, &
        &        grid%vertex_nb_idx, grid%vertex_nb_blk,                           &
        &        grid%vertex_nb_stencil, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed (ICON)!")
    END IF

    DEALLOCATE(grid%p_patch%verts%vertex,                                    &
      &        grid%p_patch%cells%vertex_idx, grid%p_patch%cells%vertex_blk, &
      &        grid%p_patch%cells%edge_idx,   grid%p_patch%cells%edge_blk,   &
      &        grid%p_patch%cells%center, grid%p_patch%cells%area,           &
      &        grid%p_patch%edges%vertex_idx, grid%p_patch%edges%vertex_blk, &
      &        grid%p_patch%edges%cell_idx,   grid%p_patch%edges%cell_blk,   &
      &        grid%p_patch%cells%glb_index,                                 &
      &        STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed (common)!")

    IF (ALLOCATED(grid%owner_c)) THEN
      ! for "patch coverings": cell index/block owner PE
      DEALLOCATE(grid%owner_c, grid%local_c, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed (cover)!")
    END IF
    IF (ALLOCATED(grid%cov_c)) THEN
      DEALLOCATE(grid%cov_c, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed (cover, inv)!")
    END IF
  END SUBROUTINE finalize_grid


  FUNCTION get_containing_cell(grid, p_in, coord_transform) RESULT(global_idx)
    INTEGER :: global_idx 
    TYPE (t_grid),                     INTENT(INOUT) :: grid
    TYPE (t_geographical_coordinates), INTENT(IN)    :: p_in
    INTEGER,                           INTENT(IN)    :: coord_transform
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::get_containing_cell')
    TYPE (t_geographical_coordinates)                :: p   

    global_idx = -1
    SELECT CASE(grid%structure) 
    CASE (GRID_TYPE_ICON)
      ! general implementation
      IF (gnat_tree == UNASSOCIATED) THEN
        ! sequential search
        global_idx = get_containing_cell_generic(grid, p_in, coord_transform)
      ELSE
        ! using fast search tree
        global_idx = get_containing_cell_gnat(grid, p_in, coord_transform)
      END IF
    CASE (GRID_TYPE_REGULAR)
      IF (coord_transform /= LIST_DEFAULT) THEN
        p = backtransform_lambert_azimuthal(p_in, coord_transform)
      ELSE
        p = p_in
      END IF
      ! exploit regular grid structure
      global_idx = get_containing_cell_gauss(grid%regular_grid, p )
    CASE DEFAULT
      CALL finish(routine, "Unknown grid type")
    END SELECT
  END FUNCTION get_containing_cell

END MODULE mo_remap_grid
