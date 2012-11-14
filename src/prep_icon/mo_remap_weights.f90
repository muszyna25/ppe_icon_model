!! This algorithm computes the interpolation weights for conservative remapping.
!!
!! See
!! Jones, Philip W.: "First- and Second-Order Conservative Remapping Schemes
!!                    for Grids in Spherical Coordinates"
!! Monthly Weather Review,  1999, 127, pp. 2204-2210
!!
!! @author F. Prill, DWD
!!
!! @note We do not check for a pole correction where the pole is inside an element!
!!

!! Possible compiler bug: The following option directive line has proven to be necessary
!! when compiling with "-Chopt" on the NEC SX9, Rev.451 2012/06/20 [F. Prill, DWD]:
!option! -O nodarg
MODULE mo_remap_weights

!$  USE OMP_LIB

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_communication,      ONLY: idx_1d, blk_no, idx_no
  USE mo_math_constants,     ONLY: pi, pi2, pi_180
  USE mo_math_utilities,     ONLY: t_geographical_coordinates
  USE mo_mpi,                ONLY: get_my_mpi_work_id
  USE mo_util_binheap,       ONLY: t_heap_data, heap_node_init, heap_insert,  &
    &                              node_storage, nnode, get_free_node,        &
    &                              node_storage_init, node_storage_finalize
  USE mo_remap_config,       ONLY: dbg_level
  USE mo_remap_io,           ONLY: s_maxsize
  USE mo_remap_shared,       ONLY: t_line, t_poly, t_grid, inside,            &
    &                              edge2line, cell2poly, intersect,           &
    &                              compute_intersection, dist_deg, dist_cc,   &
    &                              LIST_NPOLE, LIST_SPOLE, LIST_DEFAULT,      &
    &                              transform_lambert_azimuthal,               &
    &                              backtransform_lambert_azimuthal,           &
    &                              allocate_lookup_tables,                    &
    &                              compute_coordinate_transform,              &
    &                              clear_lookup_tables,                       &
    &                              deallocate_lookup_tables,                  &
    &                              copy_vertex_coords,                        &
    &                              compute_length, ccw_orientation,           &
    &                              pole_thresh, compute_index_lists,          &
    &                              finalize_index_lists, LIST_NAME,           &
    &                              finalize_vertex_coords, npole, spole
  USE mo_remap_grid,         ONLY: get_containing_cell
  USE mo_remap_intp,         ONLY: t_intp_data, t_intp_data_mt,               &
    &                              reduce_mthreaded_weights,                  &
    &                              allocate_intp_data,  deallocate_intp_data, &
    &                              merge_heaps, sync_foreign_wgts
 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: prepare_interpolation
  PUBLIC :: consistency_check

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_weights')
 
CONTAINS

  !> Compute intersection of a line with a given edge.
  !  returns indices of next cell and edge to proceed with the remaining segment.
  !
  SUBROUTINE intersect_segment(s, s_beg, cc_idx, cc_blk, ee_idx, ee_blk, grid, coord_transform, &
    &                          thresh, out_s1, out_s2, out_c_k_idx, out_c_k_blk,                &
    &                          out_e_k_idx, out_e_k_blk, istat, lthresh, l_determine_cell)

    TYPE (t_line), INTENT(INOUT) :: s
    TYPE (t_geographical_coordinates), INTENT(IN) :: s_beg  !< segment begin point (in degrees)
    INTEGER,       INTENT(IN)    :: cc_idx, cc_blk, ee_idx, ee_blk
    TYPE(t_grid),  INTENT(IN)    :: grid
    INTEGER,       INTENT(IN)    :: coord_transform
    REAL(wp),      INTENT(IN)    :: thresh

    TYPE (t_line), INTENT(OUT)   :: out_s1, out_s2            !< remaining segments
    INTEGER,       INTENT(OUT)   :: out_c_k_idx, out_c_k_blk  !< new cell indices
    INTEGER,       INTENT(OUT)   :: out_e_k_idx, out_e_k_blk  !< new edge indices
    LOGICAL,       INTENT(OUT)   :: istat                     !< return code (status)
    LOGICAL,       INTENT(INOUT) :: lthresh
    LOGICAL,       INTENT(INOUT) :: l_determine_cell

    ! local variables
    REAL(wp), PARAMETER :: ZERO_THRESHOLD = 1e-12_wp
    INTEGER                           :: k, e_k_idx, e_k_blk, c_k_idx, c_k_blk, &
      &                                  i_ne, nints
    TYPE (t_line)                     :: line_k, full_segment
    TYPE (t_geographical_coordinates) :: p_k, s_end, endp, p_k_2
    TYPE (t_poly)                     :: cell_cc
    REAL(wp)                          :: rdist, mdist, dlon, dlat, roffset
    LOGICAL                           :: l_cut_threshold_n, l_cut_threshold_s

    istat = .FALSE.
    roffset = 0._wp
    ! build cell polygon "cc"
    CALL cell2poly(grid, cc_idx,cc_blk, coord_transform, cell_cc, i_ne)
    out_c_k_idx = cc_idx
    out_c_k_blk = cc_blk
    out_s1 = s

    IF (coord_transform == LIST_DEFAULT) THEN
      endp = s%p(2)
      full_segment%p(1) = s_beg
    ELSE
      endp = backtransform_lambert_azimuthal(s%p(2),coord_transform)
      full_segment%p(1) = transform_lambert_azimuthal(s_beg,coord_transform)
    END IF
    full_segment%p(2) = s%p(2)

    ! Test if second point inside
    s_end = s%p(2)
    IF ((s_end%lon < -90._wp) .AND. (MAXVAL(cell_cc%p(1:i_ne)%lon) >  90._wp)) THEN
      s_end%lon = s_end%lon + 360._wp
    ELSE IF ((s_end%lon >  90._wp) .AND. (MINVAL(cell_cc%p(1:i_ne)%lon) < -90._wp)) THEN
      s_end%lon = s_end%lon - 360._wp
    END IF
    IF (inside(s_end, i_ne, cell_cc)) THEN
      istat = .FALSE.
    ELSE 
      ! Remaining segment s.p1 -- s.p2
      mdist = -1._wp
      nints = 0 ! no. of intersections
      LOOP : DO k=1,grid%p_patch%cell_type
        e_k_idx = grid%p_patch%cells%edge_idx(cc_idx,cc_blk,k)
        e_k_blk = grid%p_patch%cells%edge_blk(cc_idx,cc_blk,k)
        ! looking for intersection with edge e_k
        IF (l_determine_cell .OR.  &
          & ((e_k_idx /= ee_idx) .OR. (e_k_blk /= ee_blk))) THEN
          line_k = edge2line(grid, e_k_idx, e_k_blk, coord_transform)
          IF (intersect(line_k,s)) THEN
            nints = nints + 1
            ! found  intersection with edge e_k
            CALL compute_intersection(full_segment, line_k, p_k);
            ! make sure we don't lie *on* the edge
            IF (coord_transform /= LIST_DEFAULT) THEN
              rdist = dist_cc(s%p(1), p_k)
            ELSE
              rdist = dist_deg(s%p(1), p_k)
            END IF
            IF (rdist>mdist)  THEN
              mdist = rdist
              ! get cell c_k adjacent to edge e_k, neighbouring cc:
              c_k_idx = grid%p_patch%edges%cell_idx(e_k_idx, e_k_blk, 1)
              c_k_blk = grid%p_patch%edges%cell_blk(e_k_idx, e_k_blk, 1)
              IF ((c_k_idx == cc_idx) .AND. (c_k_blk == cc_blk)) THEN
                c_k_idx = grid%p_patch%edges%cell_idx(e_k_idx, e_k_blk, 2)
                c_k_blk = grid%p_patch%edges%cell_blk(e_k_idx, e_k_blk, 2)
              END IF

              ! We may cross the "pole threshold" with this
              ! segment. In this case, divide the segment at the
              ! intersection point with the threshold and return only
              ! the first part.
              ! - if we did this already for the last subsegment, skip this step.

              ! compute endpoint latitude
              IF (coord_transform == LIST_DEFAULT) THEN
                s_end = p_k
              ELSE
                s_end = backtransform_lambert_azimuthal(p_k,coord_transform)
              END IF

              l_cut_threshold_n = .FALSE.
              l_cut_threshold_s = .FALSE.

              IF (.NOT. lthresh) THEN
                ! north pole region, crossing from south to north
                IF ((coord_transform == LIST_DEFAULT)      .AND. &
                  & (s_end%lat > 0._wp) .AND. (ABS( 90._wp - s_end%lat) < thresh)) THEN
                  roffset =  1.0_wp
                  l_cut_threshold_n = .TRUE.
                ! north pole region, crossing from north to south
                ELSE IF ((coord_transform == LIST_NPOLE)   .AND. &
                  &      (ABS( 90._wp - s_end%lat) > thresh)) THEN 
                  roffset = -1.0_wp
                  l_cut_threshold_n = .TRUE.
                ! south pole region, crossing from north to south
                ELSE IF ((coord_transform == LIST_DEFAULT) .AND. &
                  &      (s_end%lat < 0._wp) .AND. (ABS(-90._wp - s_end%lat) < thresh)) THEN
                  roffset = -1.0_wp
                  l_cut_threshold_s = .TRUE.
                ! south pole region, crossing from south to north
                ELSE IF ((coord_transform == LIST_SPOLE)   .AND. &
                  &      (ABS(-90._wp - s_end%lat) > thresh)) THEN 
                  roffset =  1.0_wp
                  l_cut_threshold_s = .TRUE.
                END IF

                lthresh = l_cut_threshold_n .OR. l_cut_threshold_s
                IF (lthresh) THEN
                  ! determine intersection between segment and threshold line:
                  dlat = endp%lat - s_beg%lat
                  dlon = endp%lon - s_beg%lon
                  IF (dlon >  180.0_wp) dlon = dlon - 360._wp
                  IF (dlon < -180.0_wp) dlon = dlon + 360._wp
                  IF (l_cut_threshold_n) THEN
                    p_k_2%lat =  90._wp - thresh + roffset*1.e-12
                  ELSE
                    p_k_2%lat = -90._wp + thresh + roffset*1.e-12
                  END IF
                  IF (ABS(dlat) < ZERO_THRESHOLD) THEN
                    p_k_2%lon = s_beg%lon
                  ELSE
                    p_k_2%lon = s_beg%lon + (p_k_2%lat - s_beg%lat)*dlon/dlat
                  END IF

                  IF (coord_transform /= LIST_DEFAULT) THEN
                    p_k   = transform_lambert_azimuthal(p_k,coord_transform)
                    p_k_2 = transform_lambert_azimuthal(p_k_2,coord_transform)
                  END IF
                END IF
              END IF ! IF (.NOT. lthresh)

              IF (.NOT. l_cut_threshold_n .AND. .NOT. l_cut_threshold_s) THEN
                out_s1%p(1:2) = (/ s%p(1), p_k    /)
                out_s2%p(1:2) = (/ p_k,    s%p(2) /)
                out_e_k_idx = e_k_idx
                out_e_k_blk = e_k_blk
                out_c_k_idx = c_k_idx
                out_c_k_blk = c_k_blk
              ELSE
                out_s1%p(1:2) = (/ s%p(1), p_k_2  /)
                out_s2%p(1:2) = (/ p_k_2,  s%p(2) /)
                out_e_k_idx = ee_idx
                out_e_k_blk = ee_blk
                out_c_k_idx = cc_idx
                out_c_k_blk = cc_blk
              END IF
              istat = (out_c_k_idx >= 0)
              ! note: we cannot simply EXIT LOOP now, since we might
              ! have started our search accidentally on an edge,
              ! s.t. we get 2 intersections.
            END IF ! (rdist>mdist)
          END IF ! intersect
        END IF
      END DO LOOP
    END IF
  END SUBROUTINE intersect_segment


  !> Computes all intersection of a given line @l1 with edges of a grid
  !
  SUBROUTINE intersect_line(edge_v_idx, edge_v_blk, l1_in, grid, coord_transform, &
    &                       out_segments, out_cells, out_fct, icount, thresh)

    INTEGER,       INTENT(INOUT) :: edge_v_idx(2), edge_v_blk(2)
    TYPE (t_line), INTENT(IN)    :: l1_in                !< line
    TYPE(t_grid),  INTENT(INOUT) :: grid
    INTEGER,       INTENT(INOUT) :: coord_transform
    TYPE (t_line), INTENT(OUT)   :: out_segments(:)      !< list of line segments
    INTEGER,       INTENT(OUT)   :: out_cells(:,:)       !< list (idx/blk,:) of cells containing segments
    REAL(wp),      INTENT(OUT)   :: out_fct(:)           !< list of weight factors (+/- 1.0) for each segment
    INTEGER,       INTENT(INOUT) :: icount               !< no. of segments
    REAL(wp),      INTENT(IN)    :: thresh
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::intersect_line')
    TYPE (t_line)               :: l1, s, out_s1, out_s2
    INTEGER                     :: cc_idx, cc_blk, ee_idx, ee_blk,    &
      &                            out_c_k_idx, out_c_k_blk,          &
      &                            out_e_k_idx, out_e_k_blk,          &
      &                            nmax, ne, global_idx, ilist, ithrd
    LOGICAL                     :: istat, lthresh, l_determine_cell, l_lookup, &
      &                            l_firstsubsegment, lskip_segment
    TYPE (t_geographical_coordinates) ::  s_beg, l_beg
    REAL(wp)                    :: distn, dists, fct

    ithrd = 1
!$  ithrd = OMP_GET_THREAD_NUM() + 1

    ! number of cell vertices:
    ne = grid%p_patch%cell_type
    l1 = l1_in

    nmax   = MIN(UBOUND(out_segments,1), UBOUND(out_cells,1))
    IF (UBOUND(out_cells,2) /= 2)  CALL finish(routine, "Invalid dimension!")

    cc_idx = 0
    cc_blk = 0
    ee_idx = -1
    ee_blk = -1
    s      = l1
    fct    = 1.0_wp

    l_determine_cell  = .TRUE.
    lthresh           = .FALSE.
    l_lookup          = .TRUE.
    l_firstsubsegment = .TRUE.
    lskip_segment     = .FALSE.

    ! loop: successively "eat away" the segments to the next
    ! intersection point:
    LOOP : DO 
      ! decide, if we should use the Lambert transformed coords for
      ! this segment:
      s_beg = s%p(1)
      IF (coord_transform /= LIST_DEFAULT) THEN
        s_beg = backtransform_lambert_azimuthal(s_beg, coord_transform)
      END IF

      distn = ABS(npole%lat - s_beg%lat)
      dists = ABS(spole%lat - s_beg%lat)
      IF (distn < thresh) THEN
        ilist = LIST_NPOLE
      ELSE IF (dists < thresh) THEN
        ilist = LIST_SPOLE
      ELSE
        ilist = LIST_DEFAULT
      END IF

      IF ((ilist /= LIST_DEFAULT) .AND. (coord_transform == LIST_DEFAULT)) THEN
        s%p(1) = transform_lambert_azimuthal(s%p(1), ilist)
        s%p(2) = transform_lambert_azimuthal(s%p(2), ilist)
        coord_transform = ilist
      ELSE IF ((ilist == LIST_DEFAULT) .AND. (coord_transform /= LIST_DEFAULT)) THEN
        s%p(1) = s_beg
        s%p(2) = backtransform_lambert_azimuthal(s%p(2), coord_transform)
        coord_transform = ilist
      END IF

      ! determine containing cell
      IF (l_determine_cell) THEN
        ! get cell in grid 1 containing the starting vertex:
        IF (.NOT. l_lookup .OR. &
          & (grid%lookup_tbl%vertex_c_idx(edge_v_idx(1), edge_v_blk(1),coord_transform,ithrd) == -1)) THEN
!! Possible compiler bug: The following compiler directive line has proven to be necessary
!! when compiling with "-Chopt" on the NEC SX9, Rev.451 2012/06/20 [F. Prill, DWD]:
!CDIR NOIEXPAND
          global_idx = get_containing_cell(grid, s%p(1), coord_transform)
          IF (global_idx > 0) THEN
            cc_idx = idx_no(global_idx)
            cc_blk = blk_no(global_idx)
            IF (cc_idx <= 0) &
              CALL finish(routine, "Internal error: source cell could not be identified!")
            ! insert start cell into lookup table (for later use):
            grid%lookup_tbl%vertex_c_idx(edge_v_idx(1), edge_v_blk(1), coord_transform,ithrd) = cc_idx
            grid%lookup_tbl%vertex_c_blk(edge_v_idx(1), edge_v_blk(1), coord_transform,ithrd) = cc_blk
          ELSE
            IF (fct < 0._wp) THEN
              ! skip this segment if we already tried both endpoints:
              lskip_segment = .TRUE.
              EXIT LOOP
            END IF
            ! for local grids: flip edge, if there exists no containing cell:
            edge_v_idx(1:2) = (/ edge_v_idx(2), edge_v_idx(1) /)
            edge_v_blk(1:2) = (/ edge_v_blk(2), edge_v_blk(1) /)
            l1%p(1:2)       = (/ l1%p(2),       l1%p(1)       /)
            s   = l1
            fct = -1.0_wp
            CYCLE LOOP
          END IF
        ELSE
          cc_idx = grid%lookup_tbl%vertex_c_idx(edge_v_idx(1), edge_v_blk(1), coord_transform,ithrd)
          cc_blk = grid%lookup_tbl%vertex_c_blk(edge_v_idx(1), edge_v_blk(1), coord_transform,ithrd)
        END IF
        l_lookup = .FALSE.
      END IF

      IF (l_firstsubsegment) THEN
        l_beg = s_beg
        l_firstsubsegment = .FALSE.
      END IF

      CALL intersect_segment(s, l_beg, cc_idx, cc_blk, ee_idx, ee_blk, grid, coord_transform, &
        &                    thresh,  out_s1, out_s2, out_c_k_idx, out_c_k_blk,               &
        &                    out_e_k_idx, out_e_k_blk, istat, lthresh, l_determine_cell)
      l_determine_cell = .FALSE.

      ! append segment; do not append *very* short segments:
      icount = icount + 1
      IF (icount > nmax) THEN
        WRITE (0,*) "start_cell_idx = ", cc_idx, cc_blk
        CALL finish(routine, "Max size of input field exceeded!")
      END IF

      IF (coord_transform /= LIST_DEFAULT) THEN
        out_s1%p(1) = backtransform_lambert_azimuthal(out_s1%p(1), coord_transform)
        out_s1%p(2) = backtransform_lambert_azimuthal(out_s1%p(2), coord_transform)
      END IF

      out_segments(icount)  = out_s1
      out_cells(icount,1:2) = (/ cc_idx, cc_blk /)
      out_fct(icount)       = fct
      ! update indices:
      s = out_s2

      cc_idx = out_c_k_idx
      cc_blk = out_c_k_blk
      ee_idx = out_e_k_idx
      ee_blk = out_e_k_blk
      IF (.NOT. istat) EXIT LOOP
    END DO LOOP

    ! add the end point to the lookup table
    IF ((.NOT. lskip_segment) .AND. (out_c_k_idx >= 0))  THEN 
      grid%lookup_tbl%vertex_c_idx(edge_v_idx(2), edge_v_blk(2), coord_transform,ithrd) = out_cells(icount,1)
      grid%lookup_tbl%vertex_c_blk(edge_v_idx(2), edge_v_blk(2), coord_transform,ithrd) = out_cells(icount,2)
    END IF
  END SUBROUTINE intersect_line


  !> Store interpolation weight in data structure.
  !  The indices are stored in ascending order wrt. the global index.
  !
  !  @todo Add an INLINE directive here!
  !
  SUBROUTINE store_weight(weight_data, intp_data)
    TYPE(t_heap_data),    INTENT(IN)    :: weight_data
    TYPE(t_intp_data_mt), INTENT(INOUT) :: intp_data
    ! local variables 
    INTEGER  :: ithrd, hn

    ithrd = 1
!$  ithrd = OMP_GET_THREAD_NUM() + 1

    hn = get_free_node(ithrd)
    node_storage(ithrd)%v(hn) = heap_node_init(weight_data)
    CALL heap_insert(intp_data%wgt_heap(ithrd), hn, ithrd)
  END SUBROUTINE store_weight


  !> Computes (an approximation of) the line integrals.
  !
  !  Depending on the accuracy of the mapping procedure we compute
  !  line integrals (mapping weights).
  !
  ! @todo For the time being, only a first order remapping is
  !       implemented.
  !
  SUBROUTINE compute_segment_weight(segt, c1_idx, c1_blk, c2_idx, c2_blk, &
    &                               fct, grid1, grid2, intp_data1, intp_data2)
    TYPE (t_line)       , INTENT(IN)    :: segt                   !< line segment
    INTEGER             , INTENT(IN)    :: c1_idx, c1_blk         !< containing cell (in grid 1)
    INTEGER             , INTENT(IN)    :: c2_idx, c2_blk         !< cell in grid 2
    REAL(wp)            , INTENT(IN)    :: fct                    !< factors (+/- 1.0) for this segment
    TYPE (t_grid),        INTENT(IN)    :: grid1, grid2
    TYPE(t_intp_data_mt), INTENT(INOUT) :: intp_data1, intp_data2 !< interpolation weights
    ! local variables
    REAL(wp) :: sin_lat(2), dlon,  wgt
    INTEGER  :: idx_cov, idx_src, idx_local, jc_src, jb_src, jc_dst, jb_dst, nn

    dlon = pi_180 * (segt%p(2)%lon - segt%p(1)%lon)
    IF ( dlon >  pi ) THEN
      dlon = dlon - pi2
    ELSE IF ( dlon < -1._wp*pi ) THEN
      dlon = dlon + pi2
    END IF
    sin_lat(1:2) = SIN(pi_180 * segt%p(1:2)%lat)

    ! --- 1st order remapping:
    !
    !     Approximate overlap area with trapezoidal rule
    wgt =  fct * (-1._wp) * dlon*( sin_lat(1) + sin_lat(2) )/2._wp

    ! --- Store computed weights in the interpolation data structure:
    !                          [ weight     source       dest ]
    CALL store_weight(t_heap_data(wgt, c1_idx,c1_blk, c2_idx,c2_blk), intp_data2)

    ! index mapping: covering-to-local (for dest cell)
    idx_cov   = idx_1d(c1_idx,c1_blk)
    idx_local = grid1%local_c(idx_cov)
    jc_dst    = idx_no(idx_local)
    jb_dst    = blk_no(idx_local)
    IF (grid1%owner_c(idx_cov) == get_my_mpi_work_id()) THEN
      ! index mapping: local-to-covering (for source cell)
      idx_local = idx_1d(c2_idx,c2_blk)
      idx_src   = grid2%cov_c(idx_local)
      jc_src    = idx_no(idx_src)
      jb_src    = blk_no(idx_src)    
      ! weight for local PE:
      CALL store_weight(t_heap_data(wgt, jc_src,jb_src, jc_dst,jb_dst), intp_data1)
    ELSE
      ! translate source index to global index
      idx_local = idx_1d(c2_idx,c2_blk)
      idx_src   = grid2%p_patch%cells%glb_index(idx_local)
      jc_src    = idx_no(idx_src)
      jb_src    = blk_no(idx_src)    
      ! weight for foreign PE:
      nn = intp_data1%nforeign + 1
      intp_data1%nforeign = nn
      intp_data1%foreign_wgt(nn) = t_heap_data(wgt, jc_src,jb_src, jc_dst,jb_dst)
      intp_data1%foreign_pe (nn) = grid1%owner_c(idx_cov)
    END IF
  END SUBROUTINE compute_segment_weight


  !> Compute intersection segments of edges.
  !
  !  Compute intersection segments of edges in p_grid2 with edges in
  !  p_grid1.
  !
  SUBROUTINE compute_segments(grid1, grid2, intp_data1, intp_data2, thresh)
    TYPE (t_grid),        INTENT(INOUT) :: grid1, grid2
    TYPE(t_intp_data_mt), INTENT(INOUT) :: intp_data1, intp_data2
    REAL(wp),             INTENT(IN)    :: thresh
    ! local variables
    INTEGER,          PARAMETER :: MAX_NSEGMENTS       = 5000
    REAL(wp),         PARAMETER :: TOL_CONST_LONGITUDE = 1.e-12_wp

    INTEGER           jc, jb, ilist, iidx, nidx, c_idx, c_blk,     &
      &               i, iedge, nsegments, &
      &               nsegments0, coord_transform
    TYPE (t_line)  :: l1, out_s(MAX_NSEGMENTS)    !< list of line segments
    INTEGER        :: out_c(MAX_NSEGMENTS,2)      !< list (:,idx/blk) of cells containing segments
    REAL(wp)       :: out_fct(MAX_NSEGMENTS)      !< list of weight factors (+/- 1.0) for each segment
    INTEGER        :: in_e(MAX_NSEGMENTS,2)
    INTEGER        :: edge_v_idx(2), edge_v_blk(2)
    INTEGER        :: tot_nsegments, tot_skip
    
    tot_nsegments = 0
    tot_skip      = 0

    ! loop over list segments: north, south pole region and
    ! complementary ("default") set:
    DO ilist=1,3
      ! no. of entries in index list:
      nidx = grid2%index_list%cnlist(ilist)
      IF (dbg_level >= 2)  &
        &  WRITE (0,*) "# ilist ", ilist, "(", nidx, ", ", LIST_NAME(ilist), ")"

!$OMP PARALLEL DO  &
!$OMP PRIVATE(i, iidx, c_idx, c_blk, nsegments, iedge, jc, jb, coord_transform, l1, &
!$OMP         edge_v_idx, edge_v_blk, nsegments0, out_s, out_c, in_e, out_fct)
      DO iidx=1,nidx
        ! get cell index/block from index list ("ilist"):
        c_idx = grid2%index_list%clist_idx( idx_no(iidx), blk_no(iidx), ilist )
        c_blk = grid2%index_list%clist_blk( idx_no(iidx), blk_no(iidx), ilist )

        nsegments = 0  ! start index in segment list
        coord_transform = LIST_DEFAULT
        ! loop over cell edges:
        DO iedge=1,grid2%p_patch%cell_type

          jc = grid2%p_patch%cells%edge_idx(c_idx,c_blk,iedge)
          jb = grid2%p_patch%cells%edge_blk(c_idx,c_blk,iedge)
          ! get geographical coordinates of current edge
          l1 = edge2line(grid2, jc,jb, coord_transform)

          ! get edge
          edge_v_idx(1:2) = grid2%p_patch%edges%vertex_idx(jc, jb, 1:2)
          edge_v_blk(1:2) = grid2%p_patch%edges%vertex_blk(jc, jb, 1:2)

          IF (ccw_orientation(jc, jb, c_idx, c_blk, grid2) < 0._wp) THEN
            ! swap edge vertices:
            edge_v_idx(1:2) = (/ edge_v_idx(2), edge_v_idx(1) /)
            edge_v_blk(1:2) = (/ edge_v_blk(2), edge_v_blk(1) /)
            l1%p(1:2)       = (/ l1%p(2),       l1%p(1)       /)
          END IF

          ! compute intersection segments
          nsegments0 = nsegments+1
          CALL intersect_line(edge_v_idx, edge_v_blk, l1, grid1, coord_transform, &
            &                 out_s, out_c, out_fct, nsegments, thresh)

          DO i=nsegments0,nsegments
            ! adjust start/end point (this avoids most pole intersection problems)
            out_s(nsegments0)%p(1) = grid2%p_patch%verts%vertex(edge_v_idx(1),edge_v_blk(1))
            out_s(nsegments)%p(2)  = grid2%p_patch%verts%vertex(edge_v_idx(2),edge_v_blk(2))
            in_e(i,1:2) = (/ jc,jb /)
          END DO

          tot_nsegments = tot_nsegments + nsegments - nsegments0 + 1
          IF (nsegments0 > nsegments) tot_skip = tot_skip + 1
        END DO ! cell edges
        DO i=1,nsegments
          ! skip this routine for a constant-longitude segment (no weight
          ! contribution in this case)
          IF (ABS(out_s(i)%p(2)%lon - out_s(i)%p(1)%lon) > TOL_CONST_LONGITUDE) THEN
            edge_v_idx(1:2) = grid2%p_patch%edges%vertex_idx(in_e(i,1), in_e(i,2), 1:2)
            edge_v_blk(1:2) = grid2%p_patch%edges%vertex_blk(in_e(i,1), in_e(i,2), 1:2)
            CALL compute_segment_weight(out_s(i), out_c(i,1), out_c(i,2), c_idx, c_blk, &
              &                         out_fct(i), grid1, grid2, intp_data1, intp_data2)
          END IF
        END DO
      END DO ! iidx
!$OMP END PARALLEL DO
      IF (dbg_level >= 10) THEN
        WRITE (0,*) "# tot_nsegments = ", tot_nsegments
        WRITE (0,*) "# tot_skip      = ", tot_skip
      END IF
    END DO ! ilist

  END SUBROUTINE compute_segments


  !> Computes interpolation weights.
  !
  SUBROUTINE prepare_interpolation(grid1, grid2, grid1_cov, grid2_cov, intp_data1, intp_data2)
    TYPE (t_grid),     INTENT(INOUT) :: grid1, grid2            !< local grid partition
    TYPE (t_grid),     INTENT(INOUT) :: grid1_cov, grid2_cov    !< local grid partition coverings
    TYPE(t_intp_data), INTENT(INOUT) :: intp_data1, intp_data2  !< interpolation coefficients (result)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::prepare_interpolation')
    TYPE (t_intp_data_mt) :: intp_data1_mt, intp_data2_mt
    TYPE (t_heap_data) :: t
    INTEGER  :: nthreads, i, ierrstat, ithrd
    REAL(wp) :: thresh

    nthreads = 1
!$  nthreads = omp_get_max_threads()

    ALLOCATE(node_storage(nthreads), nnode(nthreads), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    DO ithrd=1,nthreads
      CALL node_storage_init(s_maxsize, ithrd)
    END DO

    thresh = pole_thresh*MAX(grid1%char_length, grid2%char_length)
    IF (dbg_level >= 10) &
      &   WRITE (0,*) "# pole threshold = ", thresh
    CALL compute_index_lists(grid1, thresh)
    CALL compute_index_lists(grid2, thresh)

    ! compute coordinate transformation for pole regions
    CALL copy_vertex_coords(grid1)
    CALL copy_vertex_coords(grid2)
    CALL compute_coordinate_transform(grid1, LIST_NPOLE)
    CALL compute_coordinate_transform(grid2, LIST_NPOLE)
    CALL compute_coordinate_transform(grid1, LIST_SPOLE)
    CALL compute_coordinate_transform(grid2, LIST_SPOLE)

    CALL copy_vertex_coords(grid1_cov)
    CALL copy_vertex_coords(grid2_cov)
    CALL compute_coordinate_transform(grid1_cov, LIST_NPOLE)
    CALL compute_coordinate_transform(grid2_cov, LIST_NPOLE)
    CALL compute_coordinate_transform(grid1_cov, LIST_SPOLE)
    CALL compute_coordinate_transform(grid2_cov, LIST_SPOLE)

    ! allocate thread-safe data structure for interpolation weights
    CALL allocate_intp_data(intp_data1_mt, nthreads)
    CALL allocate_intp_data(intp_data2_mt, nthreads)

    ! allocate lookup tables
    CALL allocate_lookup_tables(grid1, grid2, nthreads)
    CALL clear_lookup_tables(grid1, grid2)
    CALL allocate_lookup_tables(grid1_cov, grid2_cov, nthreads)
    CALL clear_lookup_tables(grid1_cov, grid2_cov)

    ! compute line integrals, part I
    IF (dbg_level >= 2) &
      &    WRITE (0,*) "# computing segments of edges for ", TRIM(grid2%name)
    CALL compute_segments(grid1_cov, grid2, intp_data1_mt, intp_data2_mt, thresh)
    ! compute line integrals, part II
    IF (dbg_level >= 2) &
      &    WRITE (0,*) "# computing segments of edges for ", TRIM(grid1%name)
    CALL compute_segments(grid2_cov, grid1, intp_data2_mt, intp_data1_mt, thresh)
    CALL deallocate_lookup_tables(grid1, grid2)
    CALL deallocate_lookup_tables(grid1_cov, grid2_cov)
    ! clear index lists
    CALL finalize_index_lists(grid1)
    CALL finalize_index_lists(grid2)
    ! clear "working copy" with vertex coordinates
    CALL finalize_vertex_coords(grid1)
    CALL finalize_vertex_coords(grid2)
    CALL finalize_vertex_coords(grid1_cov)
    CALL finalize_vertex_coords(grid2_cov)
    ! merge multi-threaded weight lists into a single one
    CALL merge_heaps(intp_data1_mt, intp_data2_mt, nthreads)

    ! MPI computation: communicate lists of weights on PE boundaries:
    CALL sync_foreign_wgts(grid1_cov, intp_data2_mt)
    CALL sync_foreign_wgts(grid2_cov, intp_data1_mt)

!$OMP PARALLEL
!$OMP DO
    DO ithrd=1,2
      SELECT CASE(ithrd)
      CASE(1)
        CALL reduce_mthreaded_weights(intp_data1_mt, intp_data1)
        IF (dbg_level >= 10) &
          &   WRITE (0,*) "# list length 1: ", intp_data1%s_nlist
        ! clean up
        CALL deallocate_intp_data(intp_data1_mt)
      CASE(2)
        CALL reduce_mthreaded_weights(intp_data2_mt, intp_data2)
        IF (dbg_level >= 10) &
          &   WRITE (0,*) "# list length 2: ", intp_data2%s_nlist
        ! clean up
        CALL deallocate_intp_data(intp_data2_mt)
      END SELECT
    END DO
!$OMP END DO
!$OMP END  PARALLEL

    ! divide by area of destination grid cell:
    DO i=1,intp_data1%nstencil
     intp_data1%wgt(i,:,:) = intp_data1%wgt(i,:,:)/intp_data1%area(:,:)
    END DO
    DO i=1,intp_data2%nstencil
      intp_data2%wgt(i,:,:) = intp_data2%wgt(i,:,:)/intp_data2%area(:,:)
    END DO
!$OMP PARALLEL 
!$OMP DO PRIVATE(i,t)
    DO i=1,intp_data1%s_nlist
      t = intp_data1%sl(i)
      intp_data1%sl(i)%wgt  = t%wgt/intp_data1%area(t%didx,t%dblk)
    END DO
!$OMP END DO
!$OMP DO PRIVATE(i,t)
    DO i=1,intp_data2%s_nlist
      t = intp_data2%sl(i)
      intp_data2%sl(i)%wgt  = t%wgt/intp_data2%area(t%didx,t%dblk)
    END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL node_storage_finalize(1)
    DEALLOCATE(node_storage, nnode, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE prepare_interpolation


  !> Perform some consistency checks: Check for negative weights and
  !  large deviations.
  !
  SUBROUTINE consistency_check(grid, intp_data, prefix_str)
    TYPE (t_grid),     INTENT(IN) :: grid
    TYPE(t_intp_data), INTENT(IN) :: intp_data
    CHARACTER(LEN=*),  INTENT(IN) :: prefix_str
    ! local variables
    INTEGER  :: nblks, npromz, i_startidx, i_endidx, jc, jb, nfail
    REAL(wp) :: val

    WRITE (0,*) "# ", TRIM(prefix_str), ": consistency check"
    nblks  = grid%p_patch%nblks_c
    npromz = grid%p_patch%npromz_c
    nfail  = 0
!$OMP PARALLEL PRIVATE(i_startidx, i_endidx,jc,val)
!$OMP DO
    DO jb=1,nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks) i_endidx = npromz

      DO jc=i_startidx, i_endidx
        ! for local grids: cells lying outside have the value -1. (from
        ! initialization).
        IF (intp_data%area(jc,jb) < 0._wp) CYCLE
        val = intp_data%area(jc,jb)/grid%p_patch%cells%area(jc,jb)
        IF (ABS(1._wp - val) > 1e-1) THEN
!$OMP CRITICAL
          nfail = nfail + 1
!$OMP END CRITICAL
          WRITE (0,*) "#A: ", jc,jb, val
        END IF
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL
    WRITE (0,*) "# ", nfail, " cells did not pass consistency check."
  END SUBROUTINE consistency_check

END MODULE mo_remap_weights
