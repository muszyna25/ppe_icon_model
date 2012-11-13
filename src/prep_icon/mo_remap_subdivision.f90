!! subdivision of grids for MPI parallel runs,
!! partly a mock-up for ICON's module.
!!
MODULE mo_remap_subdivision

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_communication,      ONLY: idx_no, blk_no, idx_1d
  USE mo_math_utilities,     ONLY: t_geographical_coordinates
  USE mo_remap_config,       ONLY: dbg_level
  USE mo_util_sort,          ONLY: quicksort
  USE mo_remap_shared,       ONLY: t_grid, GRID_TYPE_REGULAR, GRID_TYPE_ICON
  USE mo_remap_grid_regular, ONLY: allocate_gaussian_grid, latlon_compute_area
  USE mo_remap_grid_icon,    ONLY: allocate_icon_grid
  USE mo_mpi,                ONLY: p_n_work, get_my_mpi_work_id,                &
    &                              p_int, p_comm_work,  p_allreduce_max
  !
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: decompose_grid
  PUBLIC :: create_grid_covering
  PUBLIC :: get_latitude_range
  PUBLIC :: IMIN, IMAX

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_subdivision')  

  ! generic interface for routine translating fields and indices with
  ! global grid numbering to local grid numbering
  INTERFACE translate_g2l
    MODULE PROCEDURE translate_g2l_idx
    MODULE PROCEDURE translate_g2l_2Dreal
    MODULE PROCEDURE translate_g2l_2Dint
    MODULE PROCEDURE translate_g2l_coord
  END INTERFACE

  !> constants for better readability:
  INTEGER, PARAMETER :: IMIN = 1
  INTEGER, PARAMETER :: IMAX = 2
  
CONTAINS

  SUBROUTINE translate_g2l_idx(nlocal, l2g_idx1,l2g_blk1,g2l_idx2,g2l_blk2, &
    &                          in_idx,in_blk, out_idx,out_blk)
    INTEGER, INTENT(IN)    :: nlocal
    INTEGER, INTENT(IN)    :: l2g_idx1(:,:), l2g_blk1(:,:)
    INTEGER, INTENT(IN)    :: g2l_idx2(:,:), g2l_blk2(:,:)
    INTEGER, INTENT(IN)    :: in_idx(:,:,:), in_blk(:,:,:)
    INTEGER, INTENT(INOUT) :: out_idx(:,:,:), out_blk(:,:,:)
    ! local variables
    INTEGER :: jc_local,  jb_local,  jc_global,  jb_global, &
      &        jc_local2, jb_local2, jc_global2, jb_global2, nn, &
      &        ilocal_idx, i

    nn = UBOUND(in_idx,3)
    DO ilocal_idx=1,nlocal
      jc_local  = idx_no(ilocal_idx)
      jb_local  = blk_no(ilocal_idx)
      jc_global = l2g_idx1(jc_local,jb_local)
      jb_global = l2g_blk1(jc_local,jb_local)
      DO i=1,nn
        jc_global2 = in_idx(jc_global,jb_global,i)
        jb_global2 = in_blk(jc_global,jb_global,i)
        ! take care of pentagon cells:
        IF (jc_global2 > 0) THEN
          jc_local2  = g2l_idx2(jc_global2,jb_global2)
          jb_local2  = g2l_blk2(jc_global2,jb_global2)
          out_idx(jc_local,jb_local,i) = jc_local2
          out_blk(jc_local,jb_local,i) = jb_local2
        ELSE
          out_idx(jc_local,jb_local,i) = jc_global2
          out_blk(jc_local,jb_local,i) = jb_global2
        END IF
      END DO
    END DO
  END SUBROUTINE translate_g2l_idx


  SUBROUTINE translate_g2l_2Dreal(nlocal, l2g_idx1,l2g_blk1, &
    &                             in_field,out_field)
    INTEGER, INTENT(IN)    :: nlocal
    INTEGER, INTENT(IN)    :: l2g_idx1(:,:), l2g_blk1(:,:)
    REAL(wp),INTENT(IN)    :: in_field(:,:)
    REAL(wp),INTENT(INOUT) :: out_field(:,:)
    ! local variables
    INTEGER :: jc_local,  jb_local,  jc_global,  jb_global, ilocal_idx

    DO ilocal_idx=1,nlocal
      jc_local  = idx_no(ilocal_idx)
      jb_local  = blk_no(ilocal_idx)
      jc_global = l2g_idx1(jc_local,jb_local)
      jb_global = l2g_blk1(jc_local,jb_local)
      out_field(jc_local,jb_local) = in_field(jc_global,jb_global)
    END DO
  END SUBROUTINE translate_g2l_2Dreal


  SUBROUTINE translate_g2l_2Dint(nlocal, l2g_idx1,l2g_blk1, &
    &                            in_field,out_field)
    INTEGER, INTENT(IN)    :: nlocal
    INTEGER, INTENT(IN)    :: l2g_idx1(:,:), l2g_blk1(:,:)
    INTEGER, INTENT(IN)    :: in_field(:,:)
    INTEGER, INTENT(INOUT) :: out_field(:,:)
    ! local variables
    INTEGER :: jc_local,  jb_local,  jc_global,  jb_global, ilocal_idx

    DO ilocal_idx=1,nlocal
      jc_local = idx_no(ilocal_idx)
      jb_local = blk_no(ilocal_idx)
      jc_global = l2g_idx1(jc_local,jb_local)
      jb_global = l2g_blk1(jc_local,jb_local)
      out_field(jc_local,jb_local) = in_field(jc_global,jb_global)
    END DO
  END SUBROUTINE translate_g2l_2Dint


  SUBROUTINE translate_g2l_coord(nlocal, l2g_idx1,l2g_blk1, &
    &                            in_field,out_field)
    INTEGER, INTENT(IN)    :: nlocal
    INTEGER, INTENT(IN)    :: l2g_idx1(:,:), l2g_blk1(:,:)
    TYPE(t_geographical_coordinates), INTENT(IN)    :: in_field(:,:)
    TYPE(t_geographical_coordinates), INTENT(INOUT) :: out_field(:,:)
    ! local variables
    INTEGER :: jc_local,  jb_local,  jc_global,  jb_global, ilocal_idx

    DO ilocal_idx=1,nlocal
      jc_local  = idx_no(ilocal_idx)
      jb_local  = blk_no(ilocal_idx)
      jc_global = l2g_idx1(jc_local,jb_local)
      jb_global = l2g_blk1(jc_local,jb_local)
      out_field(jc_local,jb_local) = in_field(jc_global,jb_global)
    END DO
  END SUBROUTINE translate_g2l_coord


  !> Partition grid, return a different subset for each MPI task.
  !
  !  @note This routine simply subdivides the global grid by the
  !        latitude of the cell centers.
  !  @todo For the regular grid: Instead of generating a global grid
  !        from which subpartitions are copied, it would be better
  !        to generate only the subpartitions!
  !  
  SUBROUTINE decompose_grid(grid_in, grid_out, opt_lat_range, opt_name)
    TYPE (t_grid),        INTENT(IN)    :: grid_in
    TYPE (t_grid),        INTENT(INOUT) :: grid_out
    REAL(wp), INTENT(IN), OPTIONAL      :: opt_lat_range(2)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_name !< name string (for screen messages)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::decompose_grid')
    REAL(wp),         PARAMETER :: eps     = 1e-12_wp

    INTEGER  :: nparts, n_patch_cells, n_patch_edges,             &
      &         n_patch_verts, ierrstat, jc,jb, ilocal_idx,       &
      &         start_idx, end_idx, start_blk, end_blk, nlocal_c, &
      &         nlocal_e, nlocal_v, jc_global, jb_global,         &
      &         ilocal_idx_e, ilocal_idx_v, iedge, ivert,         &
      &         jc_e, jb_e, jc_v, jb_v, min1D, max1D,             &
      &         nypoints_loc, jc_local, jb_local, i, j
    REAL(wp) :: lat_min, lat_max, lat_delta, lat
    INTEGER,  ALLOCATABLE :: g2l_cells_idx(:,:), g2l_cells_blk(:,:),  &
      &         g2l_edges_idx(:,:), g2l_edges_blk(:,:),               &
      &         g2l_verts_idx(:,:), g2l_verts_blk(:,:),               &
      &         l2g_cells_idx(:,:), l2g_cells_blk(:,:),               &
      &         l2g_edges_idx(:,:), l2g_edges_blk(:,:),               &
      &         l2g_verts_idx(:,:), l2g_verts_blk(:,:)
    LOGICAL,  ALLOCATABLE :: lverts(:,:), ledges(:,:), lvals_glb(:)
    REAL(wp), ALLOCATABLE :: yvals_glb(:)

    ! distribute grid over all working PEs:
    nparts = p_n_work
    IF (dbg_level >= 5) THEN
      IF (.NOT. PRESENT(opt_name)) THEN
        WRITE (0,*) "# partition ", TRIM(grid_in%name), " into ",nparts, " part(s)."
      ELSE
        WRITE (0,*) "# partition ", TRIM(grid_in%name), " into ",nparts, " part(s): ", opt_name
      END IF
    END IF

    IF (.NOT. PRESENT(opt_lat_range)) THEN 
      ! determine range (note: we do not check the last block, because npromz<nproma)
      SELECT CASE(grid_in%structure) 
      CASE (GRID_TYPE_ICON)
        lat_min = MINVAL(grid_in%p_patch%cells%center(:,1:(grid_in%p_patch%nblks_c-1))%lat)
        lat_max = MAXVAL(grid_in%p_patch%cells%center(:,1:(grid_in%p_patch%nblks_c-1))%lat)
        lat_delta = 1.0001_wp * (lat_max - lat_min)/nparts
        lat_min = lat_min + lat_delta*get_my_mpi_work_id()
        lat_max = lat_min + lat_delta
        IF (dbg_level >= 2)  WRITE (0,*) "# partition with lat = [", lat_min, ",", lat_max, "]"
      CASE (GRID_TYPE_REGULAR)
        ! for the regular grid: we cannot use more processors than we
        ! have grid latitudes:
        if (nparts > grid_in%regular_grid%nypoints) &
          &   CALL finish(routine, "Invalid decomposition! More PEs than grid partitions!")

        nypoints_loc = MAX((ABS(grid_in%regular_grid%nypoints)-1)/nparts + 1, 1)
        min1D = get_my_mpi_work_id() * nypoints_loc + 1
        max1D = min1D + nypoints_loc - 1
        IF (dbg_level >= 2) WRITE (0,*) "# min1D, max1D = ", min1D, max1D
        lat_min = grid_in%regular_grid%yvals1D(min1D) - eps
        lat_max = grid_in%regular_grid%yvals1D(max1D) + eps
        IF (dbg_level >= 2) WRITE (0,*) "# partition with lat = [", lat_min, ",", lat_max, "]"
      CASE DEFAULT
        CALL finish(routine, "Unknown grid type")
      END SELECT
    ELSE
      IF (dbg_level >= 5) WRITE (0,*) "# partition range: ", opt_lat_range
      SELECT CASE(grid_in%structure) 
      CASE (GRID_TYPE_ICON)
        lat_min = opt_lat_range(IMIN)
        lat_max = opt_lat_range(IMAX)
      CASE (GRID_TYPE_REGULAR)
        ! determine min, max indices of regular grid that contain the
        ! given range:
        max1D = grid_in%regular_grid%nypoints
        min1D = 1
        DO i=1,grid_in%regular_grid%nypoints
          IF (grid_in%regular_grid%yvals1D(i) < opt_lat_range(IMIN))  min1D = i
          IF ((grid_in%regular_grid%yvals1D(i) > opt_lat_range(IMAX)) .AND. &
            & (max1D == grid_in%regular_grid%nypoints)) THEN
            max1D = i
          END IF
        END DO
        nypoints_loc = max1D - min1D + 1
        lat_min = grid_in%regular_grid%yvals1D(min1D) - eps
        lat_max = grid_in%regular_grid%yvals1D(max1D) + eps
        IF (dbg_level >= 2) WRITE (0,*) "# min1D, max1D = ", min1D, max1D
      CASE DEFAULT
        CALL finish(routine, "Unknown grid type")
      END SELECT
    END IF !  (.NOT. PRESENT(opt_lat_min)) 

    ! create global-to-local mappings for cells, edges, vertices
    ALLOCATE(g2l_cells_idx(nproma,grid_in%p_patch%nblks_c), g2l_cells_blk(nproma,grid_in%p_patch%nblks_c), &
      &      g2l_edges_idx(nproma,grid_in%p_patch%nblks_e), g2l_edges_blk(nproma,grid_in%p_patch%nblks_e), &
      &      g2l_verts_idx(nproma,grid_in%p_patch%nblks_v), g2l_verts_blk(nproma,grid_in%p_patch%nblks_v), &
      &      l2g_cells_idx(nproma,grid_in%p_patch%nblks_c), l2g_cells_blk(nproma,grid_in%p_patch%nblks_c), &
      &      l2g_edges_idx(nproma,grid_in%p_patch%nblks_e), l2g_edges_blk(nproma,grid_in%p_patch%nblks_e), &
      &      l2g_verts_idx(nproma,grid_in%p_patch%nblks_v), l2g_verts_blk(nproma,grid_in%p_patch%nblks_v), &
      &      ledges(nproma,grid_in%p_patch%nblks_e), lverts(nproma,grid_in%p_patch%nblks_v),               &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    g2l_edges_idx(:,:) = -1
    g2l_verts_idx(:,:) = -1

    ilocal_idx = 0
    ! loop over global cells, mark those "owned" by this PE:
    start_blk = 1
    end_blk   = grid_in%p_patch%nblks_c
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid_in%p_patch%npromz_c
      DO jc=start_idx,end_idx
        lat = grid_in%p_patch%cells%center(jc,jb)%lat
        IF ((lat >= lat_min) .AND. (lat < lat_max)) THEN
          ilocal_idx = ilocal_idx + 1
          jc_local = idx_no(ilocal_idx)
          jb_local = blk_no(ilocal_idx)
          g2l_cells_idx(jc,jb) = jc_local
          g2l_cells_blk(jc,jb) = jb_local
          l2g_cells_idx(jc_local,jb_local) = jc
          l2g_cells_blk(jc_local,jb_local) = jb
        END IF
      END DO ! jc
    END DO !jb
    nlocal_c = ilocal_idx

    ! loop over local cells, mark adjacent edges and vertices:
    ledges(:,:) = .FALSE.
    lverts(:,:) = .FALSE.
    DO ilocal_idx=1,nlocal_c
      jc_local = idx_no(ilocal_idx)
      jb_local = blk_no(ilocal_idx)
      jc = l2g_cells_idx(jc_local, jb_local)
      jb = l2g_cells_blk(jc_local, jb_local)
      DO iedge=1,grid_in%p_patch%cell_type
        jc_e = grid_in%p_patch%cells%edge_idx(jc,jb,iedge)
        jb_e = grid_in%p_patch%cells%edge_blk(jc,jb,iedge)
        IF (jc_e > 0) THEN ! take care of pentagons
          ledges(jc_e, jb_e) = .TRUE.
        END IF
      END DO
      DO ivert=1,grid_in%p_patch%cell_type
        jc_v = grid_in%p_patch%cells%vertex_idx(jc,jb,ivert)
        jb_v = grid_in%p_patch%cells%vertex_blk(jc,jb,ivert)
        IF (jc_v > 0) THEN ! take care of pentagons
          lverts(jc_v,jb_v) = .TRUE.
        END IF
      END DO
    END DO

    ! Now loop all marked edges and vertices to the list. This way we
    ! preserve the ordering of edges and vertices.
    ilocal_idx_e = 0
    start_blk    = 1
    end_blk      = grid_in%p_patch%nblks_e
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid_in%p_patch%npromz_e
      DO jc=start_idx,end_idx
        IF (ledges(jc,jb)) THEN
          ilocal_idx_e = ilocal_idx_e + 1
          jc_local = idx_no(ilocal_idx_e)
          jb_local = blk_no(ilocal_idx_e)
          g2l_edges_idx(jc,jb) = jc_local
          g2l_edges_blk(jc,jb) = jb_local
          l2g_edges_idx(jc_local,jb_local) = jc
          l2g_edges_blk(jc_local,jb_local) = jb
        END IF
      END DO ! jc
    END DO !jb
    nlocal_e = ilocal_idx_e

    ilocal_idx_v = 0
    start_blk    = 1
    end_blk      = grid_in%p_patch%nblks_v
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid_in%p_patch%npromz_v
      DO jc=start_idx,end_idx
        IF (lverts(jc,jb)) THEN
          ilocal_idx_v = ilocal_idx_v + 1
          jc_local = idx_no(ilocal_idx_v)
          jb_local = blk_no(ilocal_idx_v)
          g2l_verts_idx(jc,jb) = jc_local
          g2l_verts_blk(jc,jb) = jb_local
          l2g_verts_idx(jc_local,jb_local) = jc
          l2g_verts_blk(jc_local,jb_local) = jb
        END IF
      END DO ! jc
    END DO !jb
    nlocal_v = ilocal_idx_v

    ! now copy the selected cells to the local subdomain:
    n_patch_cells = nlocal_c
    n_patch_edges = nlocal_e
    n_patch_verts = nlocal_v

    ! write some statistics to stdout:
    IF (dbg_level >= 2) THEN
      WRITE (0,*) "# grid subdivision, PE ", get_my_mpi_work_id(), ": "
      WRITE (0,*) "# cells,edges,verts:                         ", &
        &   nlocal_c,nlocal_e,nlocal_v
      WRITE (0,*) "# out of global grid with cells,edges,verts: ", &
        &   grid_in%p_patch%n_patch_cells_g, grid_in%p_patch%n_patch_edges_g, &
        &   grid_in%p_patch%n_patch_verts_g
    END IF

    grid_out%structure = grid_in%structure
    SELECT CASE(grid_in%structure) 
    CASE (GRID_TYPE_ICON)
      IF (PRESENT(opt_name)) THEN
        CALL allocate_icon_grid(grid_out, opt_name, n_patch_cells, &
          &                     n_patch_edges, n_patch_verts)
      ELSE
        CALL allocate_icon_grid(grid_out, "local ICON grid", n_patch_cells, &
          &                     n_patch_edges, n_patch_verts)
      END IF
      CALL translate_g2l(nlocal_v, l2g_verts_idx, l2g_verts_blk, &
        &                g2l_cells_idx, g2l_cells_blk,           &
        &                grid_in%p_patch%verts%cell_idx,         &
        &                grid_in%p_patch%verts%cell_blk,         &
        &                grid_out%p_patch%verts%cell_idx,        &
        &                grid_out%p_patch%verts%cell_blk)
      CALL translate_g2l(nlocal_c, l2g_cells_idx, l2g_cells_blk, &
        &                g2l_cells_idx, g2l_cells_blk,           &
        &                grid_in%p_patch%cells%neighbor_idx,     &
        &                grid_in%p_patch%cells%neighbor_blk,     &
        &                grid_out%p_patch%cells%neighbor_idx,    &
        &                grid_out%p_patch%cells%neighbor_blk)
      grid_out%vertex_nb_idx(:,:,:) = -1
      grid_out%vertex_nb_blk(:,:,:) = -1
      CALL translate_g2l(nlocal_c, l2g_cells_idx, l2g_cells_blk,       &
        &                g2l_cells_idx, g2l_cells_blk,                 &
        &                grid_in%vertex_nb_idx, grid_in%vertex_nb_blk, &
        &                grid_out%vertex_nb_idx,grid_out%vertex_nb_blk)
      CALL translate_g2l(nlocal_c, l2g_cells_idx, l2g_cells_blk,       &
        &                grid_in%vertex_nb_stencil, grid_out%vertex_nb_stencil)

    CASE (GRID_TYPE_REGULAR)
      IF (PRESENT(opt_name)) THEN
        CALL allocate_gaussian_grid(grid_out, opt_name, n_patch_cells, &
          &                         n_patch_edges, n_patch_verts)
      ELSE
        CALL allocate_gaussian_grid(grid_out, "local Gaussian grid", n_patch_cells, &
          &                         n_patch_edges, n_patch_verts)
      END IF
      grid_out%regular_grid%nxpoints = grid_in%regular_grid%nxpoints
      grid_out%regular_grid%nypoints = nypoints_loc
      ALLOCATE(grid_out%regular_grid%xvals1D(grid_out%regular_grid%nxpoints),    &
        &      grid_out%regular_grid%yvals1D(grid_out%regular_grid%nypoints),    &
        &      grid_out%regular_grid%xperm1D(grid_out%regular_grid%nxpoints),    &
        &      grid_out%regular_grid%yperm1D(grid_out%regular_grid%nypoints),    &
        &      grid_out%regular_grid%edge_lon(grid_out%regular_grid%nxpoints),   &
        &      grid_out%regular_grid%edge_lat(grid_out%regular_grid%nypoints+1), &
        &      STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      grid_out%regular_grid%xvals1D(:)  = grid_in%regular_grid%xvals1D(:)
      grid_out%regular_grid%xperm1D(:)  = grid_in%regular_grid%xperm1D(:)
      ! create a new array of latitudes (and a new permutation array:
      ALLOCATE(yvals_glb(grid_in%regular_grid%nypoints), &
        &      lvals_glb(grid_in%regular_grid%nypoints), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
      lvals_glb(:) = .FALSE.
      DO i=min1D,max1D
        yvals_glb(grid_in%regular_grid%yperm1D(i)) = grid_in%regular_grid%yvals1D(i)
        lvals_glb(grid_in%regular_grid%yperm1D(i)) = .TRUE.
      END DO
      j=0
      DO i=1,SIZE(yvals_glb)
        IF (lvals_glb(i)) THEN
          j=j+1
          yvals_glb(j) = yvals_glb(i)
        END IF
      END DO
      IF (j /= grid_out%regular_grid%nypoints) THEN
        CALL finish(routine, "Internal error! Inconsistent sizes!")
      END IF
      grid_out%regular_grid%yvals1D(1:grid_out%regular_grid%nypoints) = yvals_glb(1:grid_out%regular_grid%nypoints)
      grid_out%regular_grid%yperm1D = (/ (i, i=1,grid_out%regular_grid%nypoints) /)
      CALL quicksort(grid_out%regular_grid%yvals1D, grid_out%regular_grid%yperm1D)
      DEALLOCATE(yvals_glb, lvals_glb, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

      grid_out%regular_grid%edge_lon(:) = grid_in%regular_grid%edge_lon(:)
      grid_out%regular_grid%edge_lat(:) = grid_in%regular_grid%edge_lat(min1D:max1D+1)
    CASE DEFAULT
      CALL finish(routine, "Unknown grid type")
    END SELECT

    CALL translate_g2l(nlocal_v, l2g_verts_idx, l2g_verts_blk, &
      &                grid_in%p_patch%verts%vertex, grid_out%p_patch%verts%vertex)
    CALL translate_g2l(nlocal_c, l2g_cells_idx, l2g_cells_blk, &
      &                grid_in%p_patch%cells%center, grid_out%p_patch%cells%center)
    CALL translate_g2l(nlocal_c, l2g_cells_idx, l2g_cells_blk, &
      &                grid_in%p_patch%cells%area, grid_out%p_patch%cells%area)
    CALL translate_g2l(nlocal_c, l2g_cells_idx, l2g_cells_blk, &
      &                g2l_verts_idx, g2l_verts_blk,           &
      &                grid_in%p_patch%cells%vertex_idx,       &
      &                grid_in%p_patch%cells%vertex_blk,       &
      &                grid_out%p_patch%cells%vertex_idx,      &
      &                grid_out%p_patch%cells%vertex_blk)
    CALL translate_g2l(nlocal_c, l2g_cells_idx, l2g_cells_blk, &
      &                g2l_edges_idx, g2l_edges_blk,           &
      &                grid_in%p_patch%cells%edge_idx,         &
      &                grid_in%p_patch%cells%edge_blk,         & 
      &                grid_out%p_patch%cells%edge_idx,        &
      &                grid_out%p_patch%cells%edge_blk)
    CALL translate_g2l(nlocal_e, l2g_edges_idx, l2g_edges_blk, &
      &                g2l_verts_idx, g2l_verts_blk,           &
      &                grid_in%p_patch%edges%vertex_idx,       &
      &                grid_in%p_patch%edges%vertex_blk,       &
      &                grid_out%p_patch%edges%vertex_idx,      &
      &                grid_out%p_patch%edges%vertex_blk)
    CALL translate_g2l(nlocal_e, l2g_edges_idx, l2g_edges_blk, &
      &                g2l_cells_idx, g2l_cells_blk,           &
      &                grid_in%p_patch%edges%cell_idx,         &
      &                grid_in%p_patch%edges%cell_blk,         &
      &                grid_out%p_patch%edges%cell_idx,        &
      &                grid_out%p_patch%edges%cell_blk)

    ! set local-to-global mapping:
    DO ilocal_idx=1,nlocal_c
      jc_local  = idx_no(ilocal_idx)
      jb_local  = blk_no(ilocal_idx)
      jc_global = l2g_cells_idx(jc_local,jb_local)
      jb_global = l2g_cells_blk(jc_local,jb_local)
      grid_out%p_patch%cells%glb_index(ilocal_idx) = idx_1d(jc_global,jb_global)
    END DO

    ! copy local grid sizes:
    grid_out%p_patch%n_patch_cells_g = grid_in%p_patch%n_patch_cells_g
    grid_out%p_patch%n_patch_edges_g = grid_in%p_patch%n_patch_edges_g
    grid_out%p_patch%n_patch_verts_g = grid_in%p_patch%n_patch_verts_g
    ! copy characteristic grid length
    grid_out%char_length = grid_in%char_length
    ! clean up
    DEALLOCATE(g2l_cells_idx,g2l_cells_blk,g2l_edges_idx,g2l_edges_blk, &
      &        g2l_verts_idx,g2l_verts_blk,l2g_cells_idx,l2g_cells_blk, &
      &        l2g_edges_idx,l2g_edges_blk,l2g_verts_idx,l2g_verts_blk, &
      &        lverts, ledges,                                          &
      &        STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE decompose_grid


  !> Decomposes grid and creates a grid covering for each MPI task.
  !
  !  For distributed coefficient computation, we create "patch
  !  coverings", i.e.  local grid partitions used as source grids,
  !  which are larger than the local grid partition.
  !  
  SUBROUTINE create_grid_covering(grid_in, grid_loc, grid_cov, opt_lat_range)
    TYPE (t_grid),        INTENT(IN)    :: grid_in   !< Global input grid
    TYPE (t_grid),        INTENT(INOUT) :: grid_loc  !< local grid partition for this MPI task
    TYPE (t_grid),        INTENT(INOUT) :: grid_cov  !< local grid covering for this MPI task

    !> Optional latitude range for grid covering. If missing, we copy
    !  the global grid as a trivial covering:
    REAL(wp), INTENT(IN), OPTIONAL      :: opt_lat_range(2)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::create_grid_covering')
    INTEGER :: n_patch_cells, n_patch_edges,        &
      &        n_patch_verts, ierrstat, i, iglobal_idx,     &
      &        this_pe, icount
    INTEGER, ALLOCATABLE :: owner_glb(:)

    n_patch_cells = grid_in%p_patch%n_patch_cells
    n_patch_edges = grid_in%p_patch%n_patch_edges
    n_patch_verts = grid_in%p_patch%n_patch_verts
    IF (PRESENT(opt_lat_range)) THEN
      SELECT CASE(grid_in%structure) 
      CASE (GRID_TYPE_ICON)
        CALL decompose_grid(grid_in, grid_cov, opt_lat_range, "local ICON covering")
      CASE (GRID_TYPE_REGULAR)
        CALL decompose_grid(grid_in, grid_cov, opt_lat_range, "local Gaussian covering")
      CASE DEFAULT
        CALL finish(routine, "Unknown grid type")
      END SELECT
    ELSE
      grid_cov%structure = grid_in%structure

      SELECT CASE(grid_in%structure) 
      CASE (GRID_TYPE_ICON)
        CALL allocate_icon_grid(grid_cov, "local ICON covering", n_patch_cells, &
          &                     n_patch_edges, n_patch_verts)
        ! simply copy the data:
        grid_cov%p_patch%verts%cell_idx(:,:,:)     = grid_in%p_patch%verts%cell_idx(:,:,:)    
        grid_cov%p_patch%verts%cell_blk(:,:,:)     = grid_in%p_patch%verts%cell_blk(:,:,:)    
        grid_cov%p_patch%cells%neighbor_idx(:,:,:) = grid_in%p_patch%cells%neighbor_idx(:,:,:)
        grid_cov%p_patch%cells%neighbor_blk(:,:,:) = grid_in%p_patch%cells%neighbor_blk(:,:,:)
        grid_cov%vertex_nb_idx(:,:,:)              = grid_in%vertex_nb_idx(:,:,:)
        grid_cov%vertex_nb_blk(:,:,:)              = grid_in%vertex_nb_blk(:,:,:)
        grid_cov%vertex_nb_stencil(:,:)            = grid_in%vertex_nb_stencil(:,:)
      CASE (GRID_TYPE_REGULAR)
        CALL allocate_gaussian_grid(grid_cov, "local Gaussian covering", n_patch_cells, &
          &                         n_patch_edges, n_patch_verts)
        ! simply copy the data:
        grid_cov%regular_grid%nxpoints = grid_in%regular_grid%nxpoints
        grid_cov%regular_grid%nypoints = grid_in%regular_grid%nypoints
        ALLOCATE(grid_cov%regular_grid%xvals1D(grid_cov%regular_grid%nxpoints), &
          &      grid_cov%regular_grid%yvals1D(grid_cov%regular_grid%nypoints), &
          &      grid_cov%regular_grid%xperm1D(grid_cov%regular_grid%nxpoints), &
          &      grid_cov%regular_grid%yperm1D(grid_cov%regular_grid%nypoints), &
          &      grid_cov%regular_grid%edge_lon(grid_cov%regular_grid%nxpoints),   &
          &      grid_cov%regular_grid%edge_lat(grid_cov%regular_grid%nypoints+1), &
          &      STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
        grid_cov%regular_grid%xvals1D(:) = grid_in%regular_grid%xvals1D(:)
        grid_cov%regular_grid%yvals1D(:) = grid_in%regular_grid%yvals1D(:)
        grid_cov%regular_grid%xperm1D(:) = grid_in%regular_grid%xperm1D(:)
        grid_cov%regular_grid%yperm1D(:) = grid_in%regular_grid%yperm1D(:)
        grid_cov%regular_grid%edge_lon(:) = grid_in%regular_grid%edge_lon(:)
        grid_cov%regular_grid%edge_lat(:) = grid_in%regular_grid%edge_lat(:)
      CASE DEFAULT
        CALL finish(routine, "Unknown grid type")
      END SELECT
      ! simply copy the data:
      grid_cov%p_patch%verts%vertex(:,:)        = grid_in%p_patch%verts%vertex(:,:)        
      grid_cov%p_patch%cells%vertex_idx(:,:,:)  = grid_in%p_patch%cells%vertex_idx(:,:,:)  
      grid_cov%p_patch%cells%vertex_blk(:,:,:)  = grid_in%p_patch%cells%vertex_blk(:,:,:)  
      grid_cov%p_patch%cells%edge_idx(:,:,:)    = grid_in%p_patch%cells%edge_idx(:,:,:)    
      grid_cov%p_patch%cells%edge_blk(:,:,:)    = grid_in%p_patch%cells%edge_blk(:,:,:)    
      grid_cov%p_patch%cells%center(:,:)        = grid_in%p_patch%cells%center(:,:)        
      grid_cov%p_patch%cells%area(:,:)          = grid_in%p_patch%cells%area(:,:)          
      grid_cov%p_patch%edges%vertex_idx(:,:,:)  = grid_in%p_patch%edges%vertex_idx(:,:,:)  
      grid_cov%p_patch%edges%vertex_blk(:,:,:)  = grid_in%p_patch%edges%vertex_blk(:,:,:)  
      grid_cov%p_patch%edges%cell_idx(:,:,:)    = grid_in%p_patch%edges%cell_idx(:,:,:)    
      grid_cov%p_patch%edges%cell_blk(:,:,:)    = grid_in%p_patch%edges%cell_blk(:,:,:)    
      grid_cov%p_patch%cells%glb_index(:)       = grid_in%p_patch%cells%glb_index(:)       
      grid_cov%p_patch%n_patch_cells_g          = grid_in%p_patch%n_patch_cells_g
      grid_cov%p_patch%n_patch_edges_g          = grid_in%p_patch%n_patch_edges_g
      grid_cov%p_patch%n_patch_verts_g          = grid_in%p_patch%n_patch_verts_g
      ! copy characteristic grid length
      grid_cov%char_length = grid_in%char_length
    END IF

    ! create owner information for local covering:
    ! - allocate a global field, set local entries
    this_pe = get_my_mpi_work_id()
    ALLOCATE(owner_glb(grid_in%p_patch%n_patch_cells_g), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    owner_glb(:) = -1
    DO i=1,grid_loc%p_patch%n_patch_cells
      iglobal_idx = grid_loc%p_patch%cells%glb_index(i)
      owner_glb(iglobal_idx) = this_pe
    END DO
    ! - reduce (merge) to a global field with MPI
    !
    !   this is indenpendent from the no. of processors as long as the
    !   grid decomposition is a valid partitioning into disjoint
    !   subsets
    CALL p_allreduce_max(owner_glb, p_comm_work)
    ! - allocate local owner and index fields in covering data structure
    ALLOCATE(grid_cov%owner_c(grid_cov%p_patch%n_patch_cells), &
      &      grid_cov%local_c(grid_cov%p_patch%n_patch_cells), &
      &      STAT=ierrstat)
    ! - allocate index mapping: local-to-covering
    ALLOCATE(grid_loc%cov_c(grid_loc%p_patch%n_patch_cells), & 
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ! - set local owner field
    icount = 0
    DO i=1,grid_cov%p_patch%n_patch_cells
      iglobal_idx = grid_cov%p_patch%cells%glb_index(i)
      grid_cov%owner_c(i) = owner_glb(iglobal_idx)
      IF (grid_cov%owner_c(i) == this_pe) icount = icount + 1
    END DO
    ! consistency check:
    IF (icount /= grid_loc%p_patch%n_patch_cells) CALL finish(routine, "Internal error.")

    ! - set index mapping: covering-to-local, reuse the global field
    !   for this purpose:
    owner_glb(:) = -1
    DO i=1,grid_loc%p_patch%n_patch_cells
      iglobal_idx = grid_loc%p_patch%cells%glb_index(i)
      owner_glb(iglobal_idx) = i
    END DO
    CALL p_allreduce_max(owner_glb, p_comm_work)
    grid_cov%local_c(:) = -1
    DO i=1,grid_cov%p_patch%n_patch_cells
      iglobal_idx = grid_cov%p_patch%cells%glb_index(i)
      grid_cov%local_c(i) = owner_glb(iglobal_idx)
      IF (grid_cov%owner_c(i) == this_pe) THEN
        grid_loc%cov_c(owner_glb(iglobal_idx)) = i
      END IF
    END DO

    ! - clean up
    DEALLOCATE(owner_glb, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE create_grid_covering


  !> @return minimum, maximum latitude of a given grid partition.
  !
  FUNCTION get_latitude_range(grid) RESULT(lat_range)
    TYPE (t_grid),  INTENT(IN)    :: grid           !< input grid
    REAL(wp)                      :: lat_range(2)   !< latitude range
    ! local variables
    INTEGER :: start_idx, end_idx, start_blk, end_blk, jc, jb

    lat_range(IMIN) = grid%p_patch%verts%vertex(1,1)%lat
    lat_range(IMAX) = lat_range(IMIN)

    start_blk = 1
    end_blk   = grid%p_patch%nblks_v
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid%p_patch%npromz_v
      DO jc=start_idx,end_idx
        lat_range(IMIN) = MIN(lat_range(IMIN), grid%p_patch%verts%vertex(jc,jb)%lat)
        lat_range(IMAX) = MAX(lat_range(IMAX), grid%p_patch%verts%vertex(jc,jb)%lat)
      END DO ! jc
    END DO !jb
    ! enlarge the latitude range (for a suitable covering)
    lat_range(IMIN) = lat_range(IMIN) - 2._wp*grid%char_length
    lat_range(IMAX) = lat_range(IMAX) + 2._wp*grid%char_length
  END FUNCTION get_latitude_range

END MODULE mo_remap_subdivision
