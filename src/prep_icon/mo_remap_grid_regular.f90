!
! @todo cell centers are only necessary for debugging purposes and can
! be removed later.
!
MODULE mo_remap_grid_regular

  USE mo_kind,              ONLY: wp
  USE mo_parallel_config,   ONLY: nproma
  USE mo_communication,     ONLY: blk_no, idx_no, idx_1d
  USE mo_exception,         ONLY: finish
  USE mo_impl_constants,    ONLY: SUCCESS
  USE mo_math_constants,    ONLY: pi_180
  USE mo_math_utilities,    ONLY: t_geographical_coordinates
  USE mo_remap_config,      ONLY: dbg_level, MAX_NAME_LENGTH
  USE mo_util_sort,         ONLY: quicksort
  USE mo_mpi,               ONLY: p_comm_work, p_int, p_real_dp, get_my_mpi_work_id, &
    &                             p_bcast
  USE mo_remap_shared,      ONLY: t_grid, t_regular_grid, normalized_coord
  USE mo_remap_io,          ONLY: t_file_metadata

  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PRIVATE
  PUBLIC :: load_gaussian_grid
  PUBLIC :: allocate_gaussian_grid
  PUBLIC :: get_containing_cell
  PUBLIC :: latlon_compute_area

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_grid_regular')  

  INTEGER, PARAMETER :: IND_UNDEFINED = -1

CONTAINS

  ELEMENTAL INTEGER FUNCTION regular_idx_x(j, xsize)
    INTEGER, INTENT(in) :: j, xsize
    IF(j==0) THEN
      regular_idx_x = 0
    ELSE
      regular_idx_x = SIGN(MOD(ABS(j)-1,xsize)+1, j)
    ENDIF
  END FUNCTION regular_idx_x

  ELEMENTAL INTEGER FUNCTION regular_idx_y(j, xsize)
    INTEGER, INTENT(in) :: j, xsize
    regular_idx_y = MAX((ABS(j)-1)/xsize + 1, 1)
  END FUNCTION regular_idx_y

  
  ! The connectivity information of the regular, global grid is
  ! expressed as follows:
  !
  ! - Cell   indices are in lexicographic order
  ! - Vertex indices are in lexicographic order
  ! - Edge   indices are also in lexicographic order, but we first
  !          number all zonal grid edges (contiguously) and then
  !          follow the indices of all meridional edges.
  ! - In longitudinal direction: periodicity

  !> Returns the 4 global edge indices of a given cell.
  !  counter-clockwise orientation:
  SUBROUTINE edges_of_cell_global(icell, xsize, ysize, xperm, yperm, iedges)
    INTEGER, INTENT(IN)  :: icell, xsize, ysize
    INTEGER, INTENT(OUT) :: iedges(4)
    INTEGER, INTENT(IN)  :: xperm(:), yperm(:)
    ! local variables
    INTEGER :: k1, k2, cell_row, cell_col, offset, cell1D

    cell_row = yperm( (icell-1)/xsize + 1 )
    cell_col = xperm( MOD(icell-1,xsize) + 1 )
    cell1D = (cell_row - 1)*xsize + cell_col

    offset = xsize*(ysize+1)
    k1 = offset + xsize*(cell_row-1) + cell_col
    IF (cell_col < xsize) THEN
      k2 = k1 + 1
    ELSE
      k2 = offset + xsize*(cell_row-1) + 1
    END IF

    iedges = (/ cell1D, k2, cell1D+xsize, k1 /)
  END SUBROUTINE edges_of_cell_global

  !> Returns the 4 global vertex indices of a given cell.
  !  counter-clockwise orientation:    
  SUBROUTINE vertices_of_cell_global(icell, xsize, xperm, yperm, iverts)
    INTEGER, INTENT(IN)  :: icell, xsize
    INTEGER, INTENT(OUT) :: iverts(4)
    INTEGER, INTENT(IN)  :: xperm(:), yperm(:)
    ! local variables
    INTEGER :: k1, cell_row, cell_col, cell1D

    cell_row = yperm( (icell-1)/xsize + 1 )
    cell_col = xperm( MOD(icell-1,xsize) + 1 )
    cell1D = (cell_row - 1)*xsize + cell_col

    IF (cell_col < xsize) THEN
      k1 = cell1D + 1
    ELSE
      k1 = (cell_row-1)*xsize + 1
    END IF
    iverts = (/ cell1D, k1, k1+xsize, cell1D+xsize /)
  END SUBROUTINE vertices_of_cell_global

  !> Returns the 2 (or 1) global cell indices of a given edge.
  SUBROUTINE cells_of_edge_global(iedge, xsize, ysize, xperm, yperm, icells)
    INTEGER, INTENT(IN)  :: iedge, xsize, ysize
    INTEGER, INTENT(OUT) :: icells(2)
    INTEGER, INTENT(IN)  :: xperm(:), yperm(:) !< inverse permutation
    ! local variables
    INTEGER :: m,q, edge_row, edge_col, k1, k2, i, cell_row, cell_col

    m = xsize*(ysize+1)
    IF (iedge <= m) THEN
      ! edge in zonal direction:
      icells(:) = IND_UNDEFINED
      IF (iedge >  xsize)         icells(1) = iedge-xsize
      IF (iedge <= (xsize*ysize)) icells(2) = iedge
    ELSE
      ! edge in meridional direction:
      q = iedge-m
      edge_row = (q-1)/xsize + 1
      edge_col = MOD(q-1,xsize) + 1
      k1 = (edge_row-1)*xsize + edge_col
      IF (edge_col > 1) THEN
        k2 = k1 - 1
      ELSE
        k2 = (edge_row-1)*xsize + xsize
      END IF
      icells = (/ k2, k1 /)
    END IF
    ! translate cell index
    DO i=1,2
      IF (icells(i) /= IND_UNDEFINED) THEN
        cell_row = yperm( (icells(i)-1)/xsize + 1 )
        cell_col = xperm( MOD(icells(i)-1,xsize) + 1 )
        icells(i) = (cell_row - 1)*xsize + cell_col
      END IF
    END DO
  END SUBROUTINE cells_of_edge_global

  !> Returns the 2 global vertex indices of a given edge.
  SUBROUTINE vertices_of_edge_global(iedge, xsize, ysize, iverts)
    INTEGER, INTENT(IN)  :: iedge, xsize, ysize
    INTEGER, INTENT(OUT) :: iverts(2)
    ! local variables
    INTEGER :: m
    
    iverts(:) = IND_UNDEFINED
    m = xsize*(ysize+1)
    IF (iedge <= m) THEN
      ! edge in zonal direction:
      iverts = (/ iedge, 1+MOD(iedge,xsize) + ((iedge-1)/xsize)*xsize /)
    ELSE
      ! edge in meridional direction:
      iverts = (/ iedge-m, iedge-(ysize+1)*xsize+xsize /)
    END IF
  END SUBROUTINE vertices_of_edge_global


  ! --------------------------------------------------------------------
  !> @return global index of containing cell.
  !
  ! @todo enhance performance by
  !       - assuming an equidistant grid, or
  !       - using a bisection method
  !
  ELEMENTAL FUNCTION get_containing_cell(grid, p_in) RESULT(global_idx)
    INTEGER :: global_idx 
    TYPE (t_regular_grid),             INTENT(IN) :: grid
    TYPE (t_geographical_coordinates), INTENT(IN) :: p_in
    ! local variables
    TYPE (t_geographical_coordinates) :: p
    INTEGER :: cell_row, cell_col, i, k1

    p = normalized_coord(p_in)
    cell_row = 1
    cell_col = 1

    k1 = -1
    DO i=1,(grid%nxpoints-1)
      IF ((p%lon >= grid%xvals1D(i)) .AND. &
        & (p%lon <= grid%xvals1D(i+1))) THEN
        k1=i
        IF (ABS(p%lon - grid%xvals1D(i)) < ABS(p%lon - grid%xvals1D(i+1))) THEN
          cell_col=i
        ELSE
          cell_col=i+1
        END IF
      END IF
    END DO
    IF (k1 == -1) THEN
      IF (ABS(p%lon - grid%xvals1D(1)) < ABS(p%lon+360. - grid%xvals1D(grid%nxpoints))) THEN
        cell_col=1
      ELSE
        cell_col=grid%nxpoints
      END IF
    END IF

    IF (grid%yvals1D(1) > p%lat) THEN
      cell_row = 1
    ELSE IF (grid%yvals1D(grid%nypoints) < p%lat) THEN
      cell_row = grid%nypoints
    ELSE
      DO i=1,(grid%nypoints-1)
        IF ((p%lat >= grid%yvals1D(i)) .AND. &
          & (p%lat <= grid%yvals1D(i+1))) THEN
          IF (ABS(p%lat - grid%yvals1D(i)) < ABS(p%lat - grid%yvals1D(i+1))) THEN
            cell_row=i
          ELSE
            cell_row=i+1
          END IF
        END IF
      END DO
    END IF
    global_idx = (grid%yperm1D(cell_row)-1)*grid%nxpoints + grid%xperm1D(cell_col)
  END FUNCTION get_containing_cell


  ! --------------------------------------------------------------------
  !> Load (non-reduced) Gaussian grid to internal data structure,
  !  build topology.
  !
  ! TODO[FP]:
  !
  ! - Is there a possibility to read the list of grids that are
  !   defined in a file without referring to a special variable?
  !   Otherwise: Make sure that all variables share the same gridID.
  !
  ! - When building the grid topology, we assume a global grid with
  !   periodic index definition wrt. to longitude.
  !
  ! - We expect longitude from the interval [  0., +360.] degrees
  !             latitude  from the interval [-90.,  +90.] degrees
  !
  SUBROUTINE load_gaussian_grid(grid, rank0, opt_file)
    TYPE (t_grid),          INTENT(INOUT)        :: grid
    INTEGER,                INTENT(IN)           :: rank0    !< MPI rank where file is actually read
    TYPE (t_file_metadata), INTENT(IN), OPTIONAL :: opt_file
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine    = TRIM(TRIM(modname)//'::load_gaussian_grid')
    ! flag: enable for a simple synthetic test grid
    LOGICAL,          PARAMETER :: synth_test = .FALSE.
    INTEGER :: &
      &  varID, gridID, xsize, ysize, ierrstat, idx_x, idx_y, i,     &
      &  start_idx, end_idx, start_blk, end_blk, global_idx, jb, jc, &
      &  n_patch_cells, n_patch_edges, n_patch_verts, vlistID
    REAL(wp), ALLOCATABLE :: area(:), xvals1D(:), yvals1D(:)
    INTEGER,  ALLOCATABLE :: xperm(:), yperm(:)
    CHARACTER(len=MAX_NAME_LENGTH) :: zname
    INTEGER :: gidx_vertices_of_cell(4), gidx_edges_of_cell(4), &
      &        gidx_vertices_of_edge(2), gidx_cells_of_edge(2)

    IF (synth_test) THEN
      WRITE (0,*) "# load regular grid ('synthetic' grid)."
    ELSE
      WRITE (0,*) "# load regular grid from file."
    END IF

    IF (dbg_level >= 5) WRITE (0,*) "# load 1D coordinate arrays"
    IF (synth_test) THEN
      grid%regular_grid%nxpoints = 50
      grid%regular_grid%nypoints = 50
      ALLOCATE(grid%regular_grid%xvals1D(grid%regular_grid%nxpoints),    &
        &      grid%regular_grid%yvals1D(grid%regular_grid%nypoints),    &
        &      grid%regular_grid%xperm1D(grid%regular_grid%nxpoints),    &
        &      grid%regular_grid%yperm1D(grid%regular_grid%nypoints),    &
        &      grid%regular_grid%edge_lon(grid%regular_grid%nxpoints),   &
        &      grid%regular_grid%edge_lat(grid%regular_grid%nypoints+1), &
        &      STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      xsize = grid%regular_grid%nxpoints
      ysize = grid%regular_grid%nypoints

      grid%regular_grid%xvals1D = (/ ( -180. + 360./(xsize)*(i-1), i=1,xsize ) /)
      grid%regular_grid%yvals1D = (/ (  -90. + 180./(ysize)*(i-1), i=1,ysize ) /)
    ELSE
      ! note: we read only on PE#0 and broadcast the data to the other
      ! working PEs:
      IF (get_my_mpi_work_id() == rank0) THEN
        varID = 0
        IF (.NOT. PRESENT(opt_file))  CALL finish(routine, "Internal error!")
        vlistID = opt_file%vlistID
        CALL vlistInqVarName(vlistID, varID, zname);

        gridID = vlistInqVarGrid(vlistID, varID)
        grid%regular_grid%nxpoints = gridInqXsize(gridID)
        grid%regular_grid%nypoints = gridInqYsize(gridID)
        IF (dbg_level >= 2) THEN
          WRITE (0,*) "# grid size: ", grid%regular_grid%nxpoints, &
            &         ", ", grid%regular_grid%nypoints
        END IF
      END IF
      CALL p_bcast(grid%regular_grid%nxpoints,rank0,p_comm_work)
      CALL p_bcast(grid%regular_grid%nypoints,rank0,p_comm_work)

      ! we expect the grid point number of a linear, non-reduced
      ! Gaussian grid:

      ! wave_m = 255 ! maximum wave number
      ! =>  2*wave_m+1, (2*wave_m+1)/2

      ! read latitude and longitude values:
      ALLOCATE(grid%regular_grid%xvals1D(grid%regular_grid%nxpoints),    &
        &      grid%regular_grid%yvals1D(grid%regular_grid%nypoints),    &
        &      grid%regular_grid%xperm1D(grid%regular_grid%nxpoints),    &
        &      grid%regular_grid%yperm1D(grid%regular_grid%nypoints),    &
        &      grid%regular_grid%edge_lon(grid%regular_grid%nxpoints),   &
        &      grid%regular_grid%edge_lat(grid%regular_grid%nypoints+1), &
        &      STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

      IF (get_my_mpi_work_id() == rank0) THEN
        xsize = gridInqXvals(gridID, grid%regular_grid%xvals1D)
        ysize = gridInqYvals(gridID, grid%regular_grid%yvals1D)
      ELSE
        xsize = grid%regular_grid%nxpoints
        ysize = grid%regular_grid%nypoints
      END IF
      CALL p_bcast(grid%regular_grid%xvals1D, rank0, p_comm_work)
      CALL p_bcast(grid%regular_grid%yvals1D, rank0, p_comm_work)

      IF ((xsize /= grid%regular_grid%nxpoints) .OR. (ysize /= grid%regular_grid%nypoints)) &
        &  CALL finish(routine, "Internal error!")

      !> normalize geographical coordinates (in degrees) to the 
      !  interval [-180., 180]/[-90.,90]
      DO i=1,xsize
        IF (grid%regular_grid%xvals1D(i) > 180._wp) &
          grid%regular_grid%xvals1D(i) = grid%regular_grid%xvals1D(i) - 360._wp
      END DO
    END IF
    ! make a sorted, temporary copy of xvals and yvals
    ALLOCATE(xvals1D(grid%regular_grid%nxpoints),    &
      &      yvals1D(grid%regular_grid%nypoints),    &
      &      xperm(grid%regular_grid%nxpoints),    &
      &      yperm(grid%regular_grid%nypoints),    &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")    
    xvals1D(:) = grid%regular_grid%xvals1D(:)
    yvals1D(:) = grid%regular_grid%yvals1D(:)
    grid%regular_grid%xperm1D = (/ ( i, i=1,grid%regular_grid%nxpoints ) /)
    grid%regular_grid%yperm1D = (/ ( i, i=1,grid%regular_grid%nypoints ) /)
    CALL quicksort(grid%regular_grid%xvals1D, grid%regular_grid%xperm1D)
    CALL quicksort(grid%regular_grid%yvals1D, grid%regular_grid%yperm1D)
    ! compute inverse permutation
    DO i=1,grid%regular_grid%nxpoints
      xperm(grid%regular_grid%xperm1D(i)) = i
    END DO
    DO i=1,grid%regular_grid%nypoints
      yperm(grid%regular_grid%yperm1D(i)) = i
    END DO

    ! initialize grid topology:
    n_patch_cells = xsize*ysize
    n_patch_edges = xsize*(2*ysize + 1) ! assume a *global* grid
    n_patch_verts = xsize*(  ysize + 1)

    CALL allocate_gaussian_grid(grid, "global Gaussian grid", n_patch_cells, n_patch_edges, n_patch_verts)
    grid%p_patch%n_patch_cells_g = n_patch_cells
    grid%p_patch%n_patch_edges_g = n_patch_edges
    grid%p_patch%n_patch_verts_g = n_patch_verts

    ! We will copy the structured grid information to a data structure
    ! which is of the same type as is used for the unstructured
    ! grids. This must not be problematic in terms of performance as
    ! long as we have a systematic numbering of vertices, edges, and
    ! cells, that allows to "jump" directly to each cell's neighbors.
    
    ALLOCATE(area(ysize), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! construct edge positions (located halfway between data point positions):
    IF (dbg_level >= 5) WRITE (0,*) "# construct edge positions"
    DO i=2,xsize
      grid%regular_grid%edge_lon(i) = (grid%regular_grid%xvals1D(i) + grid%regular_grid%xvals1D(i-1))/2._wp
    END DO
    DO i=2,ysize
      grid%regular_grid%edge_lat(i) = (grid%regular_grid%yvals1D(i) + grid%regular_grid%yvals1D(i-1))/2._wp
    END DO
    grid%regular_grid%edge_lon(1) = (grid%regular_grid%xvals1D(1) + 360._wp + &
      &     grid%regular_grid%xvals1D(xsize))/2._wp
    IF (grid%regular_grid%edge_lon(1) > 360._wp) THEN
      grid%regular_grid%edge_lon(1) = grid%regular_grid%edge_lon(1) - 360._wp
    END IF
    ! (works for both scan directions:)
    grid%regular_grid%edge_lat(1)       = SIGN( 90._wp, grid%regular_grid%edge_lat(2)     ) 
    grid%regular_grid%edge_lat(ysize+1) = SIGN( 90._wp, grid%regular_grid%edge_lat(ysize) )

    ! fill array with vertex positions:
    IF (dbg_level >= 5) WRITE (0,*) "# fill array with vertex positions"
    start_blk = 1
    end_blk   = grid%p_patch%nblks_v
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid%p_patch%npromz_v

      DO jc=start_idx,end_idx
        global_idx = idx_1d(jc,jb)
        idx_x = regular_idx_x(global_idx, xsize)
        idx_y = regular_idx_y(global_idx, xsize)
        grid%p_patch%verts%vertex(jc,jb)%lon = grid%regular_grid%edge_lon(idx_x)
        grid%p_patch%verts%vertex(jc,jb)%lat = grid%regular_grid%edge_lat(idx_y)
      END DO ! jc
    END DO !jb

    ! construct cell centers (located at data points)
    IF (dbg_level >= 5) WRITE (0,*) "# construct cell centers"
    start_blk = 1
    end_blk   = grid%p_patch%nblks_c
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid%p_patch%npromz_c

      DO jc=start_idx,end_idx
        global_idx = idx_1d(jc,jb)
        idx_x = regular_idx_x(global_idx, xsize)
        idx_y = regular_idx_y(global_idx, xsize)
        grid%p_patch%cells%center(jc,jb)%lon = xvals1D(idx_x)
        grid%p_patch%cells%center(jc,jb)%lat = yvals1D(idx_y)
      END DO ! jc
    END DO !jb
    ! fill out *cell* data structure with connectivity info:
    IF (dbg_level >= 5) WRITE (0,*) "# fill out *cell* data structure with connectivity info"
    start_blk = 1
    end_blk   = grid%p_patch%nblks_c
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid%p_patch%npromz_c

      DO jc=start_idx,end_idx
        global_idx = idx_1d(jc,jb)
        ! vertex of cell:
        CALL vertices_of_cell_global(global_idx, xsize, xperm, yperm, gidx_vertices_of_cell)
        grid%p_patch%cells%vertex_idx(jc,jb,1:4) = idx_no(gidx_vertices_of_cell(1:4))
        grid%p_patch%cells%vertex_blk(jc,jb,1:4) = blk_no(gidx_vertices_of_cell(1:4))
        ! edges of cell:
        CALL edges_of_cell_global(global_idx, xsize, ysize, xperm, yperm, gidx_edges_of_cell)
        grid%p_patch%cells%edge_idx(jc,jb,1:4) = idx_no(gidx_edges_of_cell(1:4))
        grid%p_patch%cells%edge_blk(jc,jb,1:4) = blk_no(gidx_edges_of_cell(1:4))
      END DO ! jc
    END DO !jb

    ! fill out *edge* data structure with connectivity info:
    IF (dbg_level >= 5) WRITE (0,*) "# fill out *edge* data structure with connectivity info"
    start_blk = 1
    end_blk   = grid%p_patch%nblks_e
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid%p_patch%npromz_e

      DO jc=start_idx,end_idx
        global_idx = idx_1d(jc,jb)
        ! vertex of edge:
        CALL vertices_of_edge_global(global_idx, xsize, ysize, gidx_vertices_of_edge)
        grid%p_patch%edges%vertex_idx(jc,jb,1:2) = idx_no(gidx_vertices_of_edge(1:2))
        grid%p_patch%edges%vertex_blk(jc,jb,1:2) = blk_no(gidx_vertices_of_edge(1:2))
        ! cell of edge:
        CALL cells_of_edge_global(global_idx, xsize, ysize, &
          &                       grid%regular_grid%xperm1D, grid%regular_grid%yperm1D, &
          &                       gidx_cells_of_edge)
        grid%p_patch%edges%cell_idx(jc,jb,1:2) = idx_no(gidx_cells_of_edge(1:2))
        grid%p_patch%edges%cell_blk(jc,jb,1:2) = blk_no(gidx_cells_of_edge(1:2))
      END DO ! jc
    END DO !jb

    ! compute cell areas:
    IF (dbg_level >= 5) WRITE (0,*) "# compute cell areas"
    CALL latlon_compute_area( grid%regular_grid, area )

    start_blk = 1
    end_blk   = grid%p_patch%nblks_c
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid%p_patch%npromz_c

      DO jc=start_idx,end_idx
        global_idx = idx_1d(jc,jb)
        idx_x = MOD(global_idx-1, xsize)+1
        idx_y = MAX((ABS(global_idx)-1)/xsize + 1, 1)
        grid%p_patch%cells%area(jc,jb) = area(yperm(idx_y))
      END DO ! jc
    END DO !jb

    ! set (trivial) local-to-global mapping:
    grid%p_patch%cells%glb_index(:) = (/ ( i, i=1,grid%p_patch%n_patch_cells) /)

    ! set characteristic grid size
    grid%char_length = ABS(grid%regular_grid%edge_lat(yperm(3))-grid%regular_grid%edge_lat(yperm(2)))

    ! clean up
    IF (dbg_level >= 5) WRITE (0,*) "# clean up"
    DEALLOCATE(area, xvals1D, yvals1D, xperm, yperm, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE load_gaussian_grid


  SUBROUTINE allocate_gaussian_grid(grid, name, n_patch_cells, n_patch_edges, n_patch_verts)
    TYPE (t_grid),    INTENT(INOUT) :: grid
    CHARACTER(LEN=*), INTENT(IN)    :: name !< name string (for screen messages)
    INTEGER,          INTENT(IN)    :: n_patch_cells, n_patch_edges, n_patch_verts
    ! local variables
    CHARACTER(LEN=*), PARAMETER   :: routine = TRIM(TRIM(modname)//'::allocate_gaussian_grid')
    INTEGER :: ierrstat

    grid%name = TRIM(name)
    grid%p_patch%cell_type = 4
    grid%p_patch%n_patch_cells = n_patch_cells
    grid%p_patch%n_patch_edges = n_patch_edges
    grid%p_patch%n_patch_verts = n_patch_verts

    ! compute the no. of blocks:
    grid%p_patch%nblks_c  = blk_no(grid%p_patch%n_patch_cells)
    grid%p_patch%npromz_c = idx_no(grid%p_patch%n_patch_cells)
    grid%p_patch%nblks_e  = blk_no(grid%p_patch%n_patch_edges)
    grid%p_patch%npromz_e = idx_no(grid%p_patch%n_patch_edges)
    grid%p_patch%nblks_v  = blk_no(grid%p_patch%n_patch_verts)
    grid%p_patch%npromz_v = idx_no(grid%p_patch%n_patch_verts)

    ! create vertices and topology info:
    IF (dbg_level >=5) WRITE (0,*) "# allocate data structures"
    ALLOCATE(grid%p_patch%verts%vertex    (nproma, grid%p_patch%nblks_v), &
      &      grid%p_patch%cells%vertex_idx(nproma, grid%p_patch%nblks_c, grid%p_patch%cell_type),  &
      &      grid%p_patch%cells%vertex_blk(nproma, grid%p_patch%nblks_c, grid%p_patch%cell_type),  &
      &      grid%p_patch%cells%edge_idx(  nproma, grid%p_patch%nblks_c, grid%p_patch%cell_type),  &
      &      grid%p_patch%cells%edge_blk(  nproma, grid%p_patch%nblks_c, grid%p_patch%cell_type),  &
      &      grid%p_patch%cells%center(   nproma, grid%p_patch%nblks_c),                           &
      &      grid%p_patch%cells%area(   nproma, grid%p_patch%nblks_c),                             &
      &      grid%p_patch%edges%vertex_idx(  nproma, grid%p_patch%nblks_e, 2),                     &
      &      grid%p_patch%edges%vertex_blk(  nproma, grid%p_patch%nblks_e, 2),                     &
      &      grid%p_patch%edges%cell_idx(    nproma, grid%p_patch%nblks_e, 2),                     &
      &      grid%p_patch%edges%cell_blk(    nproma, grid%p_patch%nblks_e, 2),                     &
      &      grid%p_patch%cells%glb_index(grid%p_patch%n_patch_cells),                             &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! set some flags to "not yet initialized":
    grid%index_list%l_initialized = .FALSE.
    grid%lookup_tbl%l_initialized = .FALSE.
  END SUBROUTINE allocate_gaussian_grid


  !> Compute area of a longitudinal strip of a regular grid
  !
  ! @note We assume an equidistant grid in longitudinal direction.
  SUBROUTINE latlon_compute_area( grid, area )
    TYPE(t_regular_grid),  INTENT(IN)    :: grid
    REAL(wp),              INTENT(INOUT) :: area(:)
    ! local variables
    REAL(wp) :: delta_lon, radius, rr_dlon
    REAL(wp) :: delta_lat(grid%nypoints)
    INTEGER  :: k, pole1, pole2, ysize

    ysize = grid%nypoints
    radius = 1.0_wp

    delta_lon   = (grid%xvals1D(2) - grid%xvals1D(1)) * pi_180
    rr_dlon     = radius*radius * delta_lon

    ! for each latitude, compute area of a grid box with a lon-lat
    ! point at its center
    DO k = 1, ysize
      delta_lat(k) = pi_180*(grid%edge_lat(k+1)-grid%edge_lat(k))
      area(k)      = 2._wp*rr_dlon * SIN(delta_lat(k)/2._wp) * COS(pi_180*grid%yvals1D(k))
    END DO
    ! treat the special case of the poles (compute area of "triangle"
    ! with one vertex at the pole and the opposite side with constant
    ! latitude)
    pole1 = ysize
    pole2 = 1
    area(pole1) = rr_dlon*(1._wp-   (SIN(pi_180*grid%edge_lat(pole1))))
    area(pole2) = rr_dlon*(1._wp-ABS(SIN(pi_180*grid%edge_lat(pole2+1))))
  END SUBROUTINE latlon_compute_area

END MODULE mo_remap_grid_regular
