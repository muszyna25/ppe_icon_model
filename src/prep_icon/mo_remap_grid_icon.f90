!! Loads an ICON grid into patch variables.
!!
!! @note We use direct NetCDF library calls since the CDI have no
!! sufficient support for INTEGER fields.
!!
!! Possible compiler bug: The following option directive DISABLES VECTORIZATION! Proven to be necessary
!! when compiling with "-Chopt" on the NEC SX9, Rev.451, running with 24x1 MPI tasks. - 2012/12/08 [F. Prill, DWD]:
!option! -Nv
MODULE mo_remap_grid_icon

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_communication,      ONLY: idx_no, blk_no
  USE mo_math_constants,     ONLY: pi, pi_180
  USE mo_physical_constants, ONLY: inverse_earth_radius
  USE mo_util_netcdf,        ONLY: nf
  USE mo_remap_config,       ONLY: dbg_level, N_VNB_STENCIL
  USE mo_mpi,                ONLY: p_comm_work, p_int, p_real_dp,   &
    &                              get_my_mpi_work_id, p_bcast
  USE mo_remap_shared,       ONLY: t_grid, normalized_coord
  USE mo_remap_io,           ONLY: t_file_metadata

  IMPLICIT NONE
  INCLUDE 'netcdf.inc'

  PRIVATE
  PUBLIC :: load_icon_grid
  PUBLIC :: allocate_icon_grid

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_grid_icon')  

CONTAINS

  ! --------------------------------------------------------------------
  !> Load ICON grid to internal data structure
  !
  SUBROUTINE load_icon_grid(grid, rank0, opt_file)
    TYPE (t_grid),          INTENT(INOUT)        :: grid
    INTEGER,                INTENT(IN)           :: rank0    !< MPI rank where file is actually read
    TYPE (t_file_metadata), INTENT(IN), OPTIONAL :: opt_file
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine    = TRIM(TRIM(modname)//'::load_icon_grid')
    INTEGER :: &
      &  i, start_idx, end_idx, start_blk, end_blk, jb, jc,        &
      &  ncid, dimid, j, idx, varID, n_patch_cells, n_patch_edges, &
      &  n_patch_verts
    REAL(wp), ALLOCATABLE :: vlon(:), vlat(:), clon(:), clat(:), area_of_c(:)
    INTEGER,  ALLOCATABLE :: v_of_c(:,:), e_of_c(:,:), c_of_c(:,:), &
      &                      v_of_e(:,:), c_of_e(:,:), c_of_v(:,:)
    REAL(wp)              :: i_rr

    IF (get_my_mpi_work_id() == rank0) THEN
      IF (.NOT. PRESENT(opt_file))  CALL finish(routine, "Internal error!")
      ncid = opt_file%ncfileID
      IF (nf_inq_dimid(ncid, 'cell', dimid) /= nf_noerr) THEN
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid))
      END IF
      CALL nf(nf_inq_dimlen(ncid, dimid, n_patch_cells))
      IF (nf_inq_dimid(ncid, 'edge', dimid) /= nf_noerr) THEN
        CALL nf(nf_inq_dimid(ncid, 'ncells2', dimid))
      END IF
      CALL nf(nf_inq_dimlen(ncid, dimid, n_patch_edges))
      IF (nf_inq_dimid(ncid, 'vertex', dimid) /= nf_noerr) THEN
        CALL nf(nf_inq_dimid(ncid, 'ncells3', dimid))
      END IF
      CALL nf(nf_inq_dimlen(ncid, dimid, n_patch_verts))
      IF (dbg_level >= 5) THEN
        WRITE (0,*) "# n_patch_cells = ", n_patch_cells
        WRITE (0,*) "# n_patch_verts = ", n_patch_verts
        WRITE (0,*) "# n_patch_edges = ", n_patch_edges
      END IF
    END IF
    CALL p_bcast(n_patch_cells,rank0,p_comm_work)
    CALL p_bcast(n_patch_verts,rank0,p_comm_work)
    CALL p_bcast(n_patch_edges,rank0,p_comm_work)

    CALL allocate_icon_grid(grid, "global ICON grid", n_patch_cells, n_patch_edges, n_patch_verts)
    grid%p_patch%n_patch_cells_g = n_patch_cells
    grid%p_patch%n_patch_edges_g = n_patch_edges
    grid%p_patch%n_patch_verts_g = n_patch_verts

    IF (dbg_level >=5) WRITE (0,*) "# load vertex coordinates from file"
    ALLOCATE(vlon(grid%p_patch%n_patch_verts))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'vlon', varID) /= nf_noerr) &
        &  CALL finish(routine, "Field <vlon> missing in grid file.")
      CALL nf(nf_get_var_double(ncid, varid, vlon))
    END IF
    CALL p_bcast(vlon,rank0,p_comm_work)
    DO i=1,grid%p_patch%n_patch_verts
      jc = idx_no(i)
      jb = blk_no(i)
      grid%p_patch%verts%vertex(jc,jb)%lon = vlon(i)/pi_180
    END DO
    DEALLOCATE(vlon)
    ALLOCATE(vlat(grid%p_patch%n_patch_verts))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'vlat', varid) /= nf_noerr) &
        &  CALL finish(routine, "Field <vlat> missing in grid file.")
      CALL nf(nf_get_var_double(ncid, varid, vlat))
    END IF
    CALL p_bcast(vlat,rank0,p_comm_work)
    DO i=1,grid%p_patch%n_patch_verts
      jc = idx_no(i)
      jb = blk_no(i)
      grid%p_patch%verts%vertex(jc,jb)%lat = vlat(i)/pi_180
    END DO    
    DEALLOCATE(vlat)
    ! convert rad->deg, normalize coordinates
    start_blk = 1
    end_blk   = grid%p_patch%nblks_v
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid%p_patch%npromz_v
      DO jc=start_idx,end_idx
        grid%p_patch%verts%vertex(jc,jb) = normalized_coord(grid%p_patch%verts%vertex(jc,jb))
      END DO ! jc
    END DO !jb

    IF (dbg_level >=5) WRITE (0,*) "# load cell-vertex indices from file"
    ALLOCATE(v_of_c(grid%p_patch%n_patch_cells,grid%p_patch%cell_type))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'vertex_of_cell', varid) /= nf_noerr) &
        &  CALL finish(routine, "Field <vertex_of_cell> missing in grid file.")
      CALL nf(nf_get_var_int(ncid, varid, v_of_c))
    END IF
    CALL p_bcast(v_of_c,rank0,p_comm_work)
    DO i=1,grid%p_patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
      DO j=1,grid%p_patch%cell_type
        idx = v_of_c(i,j)
        grid%p_patch%cells%vertex_idx(jc,jb,j) = idx_no(idx)
        grid%p_patch%cells%vertex_blk(jc,jb,j) = blk_no(idx)
      END DO
    END DO
    DEALLOCATE(v_of_c)

    IF (dbg_level >=5) WRITE (0,*) "# load cell-edge indices from file"
    ALLOCATE(e_of_c(grid%p_patch%n_patch_cells, grid%p_patch%cell_type))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'edge_of_cell', varid) /= nf_noerr) &
        &  CALL finish(routine, "Field <edge_of_cell> missing in grid file.")
      CALL nf(nf_get_var_int(ncid, varid, e_of_c))
    END IF
    CALL p_bcast(e_of_c,rank0,p_comm_work)
    DO i=1,grid%p_patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
      DO j=1,grid%p_patch%cell_type
        idx = e_of_c(i,j)
        grid%p_patch%cells%edge_idx(jc,jb,j) = idx_no(idx)
        grid%p_patch%cells%edge_blk(jc,jb,j) = blk_no(idx)
      END DO
    END DO
    DEALLOCATE(e_of_c)

    IF (dbg_level >=5) WRITE (0,*) "# load cell-neighbor indices from file"
    ALLOCATE(c_of_c(grid%p_patch%n_patch_cells,grid%p_patch%cell_type))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'neighbor_cell_index', varid) /= nf_noerr) &
        &  CALL finish(routine, "Field <neighbor_cell_index> missing in grid file.")
      CALL nf(nf_get_var_int(ncid, varid, c_of_c))
    END IF
    CALL p_bcast(c_of_c,rank0,p_comm_work)
    DO i=1,grid%p_patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
      DO j=1,grid%p_patch%cell_type
        idx = c_of_c(i,j)
        grid%p_patch%cells%neighbor_idx(jc,jb,j) = idx_no(idx)
        grid%p_patch%cells%neighbor_blk(jc,jb,j) = blk_no(idx)
      END DO
    END DO
    DEALLOCATE(c_of_c)

    IF (dbg_level >=5) WRITE (0,*) "# load edge-vertex indices from file"
    ALLOCATE(v_of_e(grid%p_patch%n_patch_edges,2))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'edge_vertices', varid) /= nf_noerr) &
        &  CALL finish(routine, "Field <edge_vertices> missing in grid file.")
      CALL nf(nf_get_var_int(ncid, varid, v_of_e))
    END IF
    CALL p_bcast(v_of_e,rank0,p_comm_work )
    DO i=1,grid%p_patch%n_patch_edges
      jc = idx_no(i)
      jb = blk_no(i)
      DO j=1,2
        idx = v_of_e(i,j)
        grid%p_patch%edges%vertex_idx(jc,jb,j) = idx_no(idx)
        grid%p_patch%edges%vertex_blk(jc,jb,j) = blk_no(idx)
      END DO
    END DO
    DEALLOCATE(v_of_e)

    IF (dbg_level >=5) WRITE (0,*) "# load edge-cell indices from file"
    ALLOCATE(c_of_e(grid%p_patch%n_patch_edges,2))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'adjacent_cell_of_edge', varid) /= nf_noerr) &
        &  CALL finish(routine, "Field <adjacent_cell_of_edge> missing in grid file.")
      CALL nf(nf_get_var_int(ncid, varid, c_of_e))
    END IF
    CALL p_bcast(c_of_e,rank0,p_comm_work )
    DO i=1,grid%p_patch%n_patch_edges
      jc = idx_no(i)
      jb = blk_no(i)
      DO j=1,2
        idx = c_of_e(i,j)
        grid%p_patch%edges%cell_idx(jc,jb,j) = idx_no(idx)
        grid%p_patch%edges%cell_blk(jc,jb,j) = blk_no(idx)
      END DO
    END DO
    DEALLOCATE(c_of_e)

    IF (dbg_level >=5) WRITE (0,*) "# load vertex-cell indices from file"
    ALLOCATE(c_of_v(grid%p_patch%n_patch_verts,6))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'cells_of_vertex', varid) /= nf_noerr) &
        &  CALL finish(routine, "Field <cells_of_vertex> missing in grid file.")
      CALL nf(nf_get_var_int(ncid, varid, c_of_v))
    END IF
    CALL p_bcast(c_of_v,rank0,p_comm_work)
    grid%p_patch%verts%cell_idx(:,:,:) = -1
    grid%p_patch%verts%cell_blk(:,:,:) = -1
    DO i=1,grid%p_patch%n_patch_verts
      jc = idx_no(i)
      jb = blk_no(i)
      DO j=1,6
        idx = c_of_v(i,j)
        ! take care of pentagon cells:
        IF (idx > 0) THEN
          grid%p_patch%verts%cell_idx(jc,jb,j) = idx_no(idx)
          grid%p_patch%verts%cell_blk(jc,jb,j) = blk_no(idx)
        END IF
      END DO
    END DO
    DEALLOCATE(c_of_v)

    IF (dbg_level >=5) WRITE (0,*) "# load cell centers from file"
    ALLOCATE(clon(grid%p_patch%n_patch_cells))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'clon', varid) /= nf_noerr) &
        &  CALL finish(routine, "Field <clon> missing in grid file.")
      CALL nf(nf_get_var_double(ncid, varid, clon))
    END IF
    CALL p_bcast(clon,rank0,p_comm_work)
    DO i=1,grid%p_patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
      grid%p_patch%cells%center(jc,jb)%lon = clon(i)/pi_180
    END DO
    DEALLOCATE(clon)
    ALLOCATE(clat(grid%p_patch%n_patch_cells))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'clat', varid) /= nf_noerr) &
        &  CALL finish(routine, "Field <clat> missing in grid file.")
      CALL nf(nf_get_var_double(ncid, varid, clat))
    END IF
    CALL p_bcast(clat,rank0,p_comm_work)
    DO i=1,grid%p_patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
      grid%p_patch%cells%center(jc,jb)%lat = clat(i)/pi_180
    END DO    
    DEALLOCATE(clat)
    ! convert rad->deg, normalize coordinates
    start_blk = 1
    end_blk   = grid%p_patch%nblks_c
    DO jb=start_blk,end_blk
      start_idx = 1
      end_idx   = nproma
      if (jb == end_blk) end_idx = grid%p_patch%npromz_c
      DO jc=start_idx,end_idx
        grid%p_patch%cells%center(jc,jb) = normalized_coord(grid%p_patch%cells%center(jc,jb))
      END DO ! jc
    END DO !jb

    IF (dbg_level >=5) WRITE (0,*) "# load cell areas from file"
    i_rr = inverse_earth_radius*inverse_earth_radius
    ALLOCATE(area_of_c(grid%p_patch%n_patch_cells))
    IF (get_my_mpi_work_id() == rank0) THEN
      IF (nf_inq_varid(ncid, 'cell_area', varid) /= nf_noerr) &
        &  CALL finish(routine, "Field <cell area> missing in grid file.")
      CALL nf(nf_get_var_double(ncid, varid, area_of_c))
    END IF
    CALL p_bcast(area_of_c,rank0,p_comm_work)
    DO i=1,grid%p_patch%n_patch_cells
      jc = idx_no(i)
      jb = blk_no(i)
      grid%p_patch%cells%area(jc,jb) = area_of_c(i)*i_rr
    END DO
    DEALLOCATE(area_of_c)

    ! set (trivial) local-to-global mapping:
    grid%p_patch%cells%glb_index(:) = (/ ( i, i=1,grid%p_patch%n_patch_cells) /)

    ! set characteristic grid size
    grid%char_length = SQRT(4*pi/grid%p_patch%n_patch_cells)/pi_180

    CALL compute_vertex_nb_index(grid)
 
  END SUBROUTINE load_icon_grid


  ! --------------------------------------------------------------------
  !> Load ICON grid to internal data structure
  !
  SUBROUTINE allocate_icon_grid(grid, name, n_patch_cells, n_patch_edges, n_patch_verts)
    TYPE (t_grid),    INTENT(INOUT) :: grid
    CHARACTER(LEN=*), INTENT(IN)    :: name !< name string (for screen messages)
    INTEGER,          INTENT(IN)    :: n_patch_cells, n_patch_edges, n_patch_verts
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::allocate_icon_grid')
    INTEGER :: ierrstat

    grid%name = TRIM(name)
    grid%p_patch%cell_type = 3
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
      &      grid%p_patch%verts%cell_idx(  nproma, grid%p_patch%nblks_v, 6),                       &
      &      grid%p_patch%verts%cell_blk(  nproma, grid%p_patch%nblks_v, 6),                       &
      &      grid%p_patch%cells%neighbor_idx(nproma, grid%p_patch%nblks_c,grid%p_patch%cell_type), &
      &      grid%p_patch%cells%neighbor_blk(nproma, grid%p_patch%nblks_c,grid%p_patch%cell_type), &
      &      grid%p_patch%cells%glb_index(grid%p_patch%n_patch_cells),                             &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ALLOCATE(grid%vertex_nb_idx(nproma, grid%p_patch%nblks_c,N_VNB_STENCIL), &
      &      grid%vertex_nb_blk(nproma, grid%p_patch%nblks_c,N_VNB_STENCIL), &
      &      grid%vertex_nb_stencil(nproma, grid%p_patch%nblks_c),            &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! set some flags to "not yet initialized":
    grid%index_list%l_initialized = .FALSE.
    grid%lookup_tbl%l_initialized = .FALSE.
  END SUBROUTINE allocate_icon_grid


  !> compute stencil indices of vertex-neighbors (cells)
  !
  ! @todo No parallel synchronization implemented!
  !
  SUBROUTINE compute_vertex_nb_index( grid )
    TYPE (t_grid), INTENT(INOUT) :: grid

    INTEGER :: ilc_n(3), ibc_n(3)       ! line and block index for neighbors of direct neighbors
    INTEGER :: ilv(3), ibv(3)           ! vertex line and block indices
    INTEGER :: ilc_v(3,6), ibc_v(3,6)   ! cell line and block indices around each of the three vertices
    INTEGER :: jb, jc, jj, jec, jtri, cnt, nblks_c, end_blk, start_idx, end_idx, &
      &        nb_idx(13), nb_blk(13)

    grid%vertex_nb_idx(:,:,:) = -1
    grid%vertex_nb_blk(:,:,:) = -1

    nblks_c  = grid%p_patch%nblks_c
    end_blk  = nblks_c

!$OMP PARALLEL

    ! The stencil consists of 12 cells surrounding the control
    ! volume that share a vertex with the control volume, plus the
    ! control volume itself.
    !
    ! Note: At pentagon points the size of the stencil reduces to 12.

!$OMP DO PRIVATE(jb,jc,jec,jj,jtri,cnt,ilv,ibv,start_idx,end_idx, &
!$OMP            ilc_v,ibc_v,ilc_n,ibc_n,nb_idx,nb_blk) 
    DO jb = 1,nblks_c
      start_idx = 1
      end_idx   = nproma
      IF (jb == end_blk) end_idx = grid%p_patch%npromz_c

      DO jc = start_idx, end_idx

        nb_idx(:) = -1
        nb_blk(:) = -1

        ! First, add the control volume itself
        cnt = 1
        nb_idx(cnt) = jc
        nb_blk(cnt) = jb

        ! get get line and block indices of cell vertices
        ilv(1:3) = grid%p_patch%cells%vertex_idx(jc,jb,1:3)
        ibv(1:3) = grid%p_patch%cells%vertex_blk(jc,jb,1:3)

        ! for each vertex: get all the cells which share this vertex
        DO jj = 1,3
          ilc_v(jj,:)=grid%p_patch%verts%cell_idx(ilv(jj),ibv(jj),:)
          ibc_v(jj,:)=grid%p_patch%verts%cell_blk(ilv(jj),ibv(jj),:)
        ENDDO

        ! 1. add the 3 direct neighbors to the stencil
        DO jec = 1, 3
          ! get line and block indices of direct neighbors
          ilc_n(jec) = grid%p_patch%cells%neighbor_idx(jc,jb,jec)
          ibc_n(jec) = grid%p_patch%cells%neighbor_blk(jc,jb,jec)

          cnt = cnt + 1
          nb_idx(cnt) = ilc_n(jec)
          nb_blk(cnt) = ibc_n(jec)
        ENDDO

        ! 2. loop over the vertices and add all the cells
        !    that are no direct neighbors and not our CV.
        DO jj = 1,3   ! loop over vertices
          DO jtri=1,6 ! loop over cells around each vertex

            IF (.NOT.( (ilc_v(jj,jtri) == ilc_n(1) .AND. ibc_v(jj,jtri) == ibc_n(1))  &
              &  .OR.  (ilc_v(jj,jtri) == ilc_n(2) .AND. ibc_v(jj,jtri) == ibc_n(2))  &
              &  .OR.  (ilc_v(jj,jtri) == ilc_n(3) .AND. ibc_v(jj,jtri) == ibc_n(3))  &
              &  .OR.  (ilc_v(jj,jtri) == jc       .AND. ibc_v(jj,jtri) == jb)        &
              &  .OR.  (ilc_v(jj,jtri) <= 0) ) ) THEN

              cnt = cnt + 1
              nb_idx(cnt) = ilc_v(jj,jtri)
              nb_blk(cnt) = ibc_v(jj,jtri)
            ENDIF
          ENDDO
        ENDDO
        grid%vertex_nb_idx(jc,jb,1:cnt) = nb_idx(1:cnt)
        grid%vertex_nb_blk(jc,jb,1:cnt) = nb_blk(1:cnt)
        grid%vertex_nb_stencil(jc,jb)   = cnt
      ENDDO ! loop over cells
    ENDDO ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_vertex_nb_index

END MODULE mo_remap_grid_icon
