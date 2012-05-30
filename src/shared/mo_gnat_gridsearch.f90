!>
!! Module for efficient promixity queries / nearest neighbor search
!! with "geometric near-neighbor access trees" (GNATs)
!!
!! uses ICON grids, operating on lon-lat coordinate system
!!
!! See
!! Brin, Sergey: "Near Neighbor Search in Large Metric Spaces"
!! VLDB '95 : Proceedings of the 21st International Conference on
!!            Very Large Data Bases,
!! Zurich Switzerland, Sept. 11--15, 1995, pp. 574-584, 
!! Morgan Kaufmann Publishers, 1995. 
!!
!! @author F. Prill, DWD
!!
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2011-08-15)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!! TODO[FP]
!! - Any beneficiary effects of NOT using REAL(wp) for distance measurements?
!! -----------------------------------------------------------------------------------
MODULE mo_gnat_gridsearch

!$  USE OMP_LIB

  USE mo_kind,                ONLY: wp, sp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_math_utilities,      ONLY: t_geographical_coordinates
  USE mo_lonlat_grid,         ONLY: t_lon_lat_grid
  USE mo_model_domain,        ONLY: t_grid_cells, t_grid_vertices, t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_mpi,                 ONLY: p_n_work, get_my_mpi_work_id, &
    &                               p_allreduce_minloc, p_comm_work
  USE mo_kind
  USE mo_grid_config,         ONLY: grid_sphere_radius

  IMPLICIT NONE

  !> REAL kind value for GNAT data structure
  INTEGER, PARAMETER  :: gk           = wp

  !> maximum distance (unit sphere!)
  REAL(gk), PARAMETER :: MAX_RANGE    = 99._gk

  !> marks a free pointer
  INTEGER,  PARAMETER :: UNASSOCIATED = -1

  !> value for invalid result (algorithm failure)
  INTEGER, PARAMETER  :: INVALID_NODE = -1
  INTEGER, PARAMETER  :: UNUSED_NODE  = -2

  !> level of output verbosity
  INTEGER, PARAMETER  :: dbg_level = 0

  !> dimension of coordinate system
  INTEGER, PARAMETER :: icoord_dim = 2

  !> maximum tree depth
  INTEGER :: tree_depth

  ! GNAT degree
  ! value ~ time consumption ~ (query time)**(-1)
  ! On the SX9 architecture, distance computations
  ! are vectorized. Thus it is beneficial to choose
  ! a shallow tree with a large degree.
#ifdef __SX__
  INTEGER, PARAMETER :: gnat_k =  30
#else
  INTEGER, PARAMETER :: gnat_k =  10
#endif

  !> coordinate type
  ! note: we do not use mo_base_geometry::t_geographical_coordinates
  !       here, since icoord_dim appears as a dimension in many lists
  TYPE  t_coord
    REAL(gk)                    :: p (icoord_dim)
    REAL(gk)                    :: sin_p, cos_p ! sin/cos of latitude
    INTEGER                     :: idx(2) ! dimensions: (idx,blk)
  END TYPE t_coord

  !> tree structure based on generalized hyperplanes
  TYPE t_gnat_tree
    INTEGER                     :: isplit_pts
    TYPE (t_coord)              :: p(gnat_k)
    INTEGER                     :: child(gnat_k) ! index pointer
    REAL(gk)                    :: drange(gnat_k, gnat_k, 2)
  END TYPE t_gnat_tree

  !> GNAT-based tree (index of root node)
  INTEGER :: gnat_tree = UNASSOCIATED

  !> the linear array where the tree nodes are actually stored:
  TYPE (t_gnat_tree), TARGET, ALLOCATABLE, SAVE :: node_storage(:)

  !> expected total number of nodes (may be altered during init)
  INTEGER :: expected_num_nodes = 100000
  !> size of array "node_storage"
  INTEGER :: max_num_nodes      = 0
  !> current usage of array "node_storage"
  INTEGER :: num_nodes          = 0

  ! public interface definition
  PRIVATE
  ! functions and subroutines
  PUBLIC :: gnat_init_grid
  PUBLIC :: gnat_query_list_mt
  PUBLIC :: gnat_destroy
  PUBLIC :: gnat_std_radius
  PUBLIC :: gnat_query_containing_triangles
  PUBLIC :: gnat_merge_distributed_queries
  ! data
  PUBLIC :: gk
  PUBLIC :: gnat_k
  PUBLIC :: t_gnat_tree, gnat_tree
  PUBLIC :: icoord_dim
  PUBLIC :: t_coord

CONTAINS

  ! -----------------------------------------------------------------------------------
  !> distance computation, uses precomputed sine, cosine
  PURE FUNCTION dist(p1, p2)
    REAL(gk)                    :: dist
    TYPE (t_coord), INTENT(IN)  :: p1
    REAL(gk),       INTENT(IN)  :: p2(icoord_dim)

    ! spherical distance:
    dist = ACOS( p1%sin_p*SIN(p2(2)) + p1%cos_p*COS(p2(2))*COS(p1%p(1)-p2(1)) )

  END FUNCTION dist


  ! -----------------------------------------------------------------------------------
  !> distance computation
  PURE FUNCTION dist_p(p1, p2)
    REAL(gk)              :: dist_p
    REAL(gk), INTENT(IN)  :: p1(icoord_dim), p2(icoord_dim)

    ! spherical distance:
    dist_p = ACOS( SIN(p1(2))*SIN(p2(2)) + COS(p1(2))*COS(p2(2))*COS(p1(1)-p2(1)) )

  END FUNCTION dist_p


  ! -----------------------------------------------------------------------------------
  !> vector variant of distance computation, uses precomputed sine, cosine
  PURE SUBROUTINE dist_vect(p1, p2, n, pdist)

    INTEGER       , INTENT(IN)     :: n 
    TYPE (t_coord), INTENT(IN)     :: p1(gnat_k)
    REAL(gk)      , INTENT(IN)     :: p2(icoord_dim)
    REAL(gk)      , INTENT(INOUT)  :: pdist(gnat_k)
    INTEGER                        :: i
    REAL(gk)                       :: sin_p2, cos_p2

    sin_p2 = SIN(p2(2))
    cos_p2 = COS(p2(2))

    ! spherical distance:
    FORALL (i=1:n)
      pdist(i) = ACOS((p1(i)%sin_p*sin_p2 +  &
        &              p1(i)%cos_p*cos_p2*COS(p1(i)%p(1)-p2(1))))
    END FORALL

  END SUBROUTINE dist_vect


  ! -----------------------------------------------------------------------------------
  !> Utility function: returns index to a free node, allocates if necessary
  FUNCTION get_free_node()
    INTEGER :: get_free_node

    num_nodes = num_nodes + 1
    get_free_node = num_nodes

  END FUNCTION get_free_node


  ! -----------------------------------------------------------------------------------
  !> Utility function: allocates space for tree nodes
  SUBROUTINE resize_node_storage()

    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_gnat_gridsearch:resize_node_storage")
    INTEGER :: new_max_num_nodes, errstat, i
    TYPE (t_gnat_tree), TARGET, ALLOCATABLE :: tmp(:)

    ! triangle copy
    new_max_num_nodes = max_num_nodes + expected_num_nodes

    IF (max_num_nodes > 0) THEN
      ALLOCATE(tmp(max_num_nodes), STAT=errstat)
      IF (errstat /= 0)  &
        CALL finish (routine, 'Error in ALLOCATE operation!')
      tmp(1:max_num_nodes) = node_storage(1:max_num_nodes)
      DEALLOCATE(node_storage, STAT=errstat)
      IF (errstat /= 0)  &
        CALL finish (routine, 'Error in DEALLOCATE operation!')
    END IF

    IF (dbg_level > 1) THEN
      WRITE(message_text,*) "allocating node storage of size=", new_max_num_nodes
      CALL message(routine, TRIM(message_text))
    END IF
    ALLOCATE(node_storage(new_max_num_nodes), STAT=errstat)
    IF (errstat /= 0)  &
      CALL finish (routine, 'Error in ALLOCATE operation!')
    FORALL (i=(max_num_nodes+1):new_max_num_nodes)
      node_storage(i)%child(:)   = UNASSOCIATED
      node_storage(i)%isplit_pts = 0
    END FORALL

    IF (max_num_nodes > 0) THEN
      node_storage(1:max_num_nodes) = tmp(1:max_num_nodes)
      DEALLOCATE(tmp, STAT=errstat)
      IF (errstat /= 0)  &
        CALL finish (routine, 'Error in DEALLOCATE operation!')
    END IF
    max_num_nodes = new_max_num_nodes

  END SUBROUTINE resize_node_storage


  ! -----------------------------------------------------------------------------------
  !> dynamic insertion of a point into GNAT-based tree
  ! Non-recursive version, suitable for parallel processing
  ! note: sequence of insertion points should be more or less random
  ! in order to provide split points which are fairly far apart.
  SUBROUTINE gnat_insert(tree, p, idx, free_node, lcomplete, count_dist)
    INTEGER,         INTENT(INOUT) :: tree ! index pointer
    REAL(gk),        INTENT(IN)    :: p(icoord_dim)
    INTEGER,         INTENT(IN)    :: idx(2)
    INTEGER,         INTENT(INOUT) :: free_node
    LOGICAL,         INTENT(OUT)   :: lcomplete
    INTEGER,         INTENT(INOUT) :: count_dist ! counter: distance calculations

    INTEGER                        :: i, j, imin(1)
    REAL(gk)                       :: pdist(gnat_k)
    TYPE (t_gnat_tree), POINTER    :: node

    node => node_storage(tree)

    !-- case 1: tree is an empty/incomplete node
    ! compute distances
    CALL dist_vect(node%p, p, node%isplit_pts, pdist)
    count_dist = count_dist + node%isplit_pts

    ! add point p as split point
    IF (node%isplit_pts < gnat_k) THEN
      i = node%isplit_pts + 1
      node%isplit_pts = i
      node%p(i)%p(:)  = p(:)
      node%p(i)%sin_p = SIN(p(2))
      node%p(i)%cos_p = COS(p(2))
      node%p(i)%idx   = idx
      node%drange(i,i,1:2) = 0._gk

      FORALL (j=1:i-1)
        node%drange(j,i,1) = MIN( node%drange(j,i,1), pdist(j) )
        node%drange(j,i,2) = MAX( node%drange(j,i,2), pdist(j) )
        node%drange(i,j,1:2) = pdist(j)
      END FORALL
      lcomplete = .TRUE.
      RETURN
    END IF

    !-- case 2: traverse into subtree for closest split point
    imin  = MINLOC(pdist)
    FORALL (j=1:gnat_k)
      node%drange(j,imin,1) = MIN( node%drange(j,imin,1), pdist(j) )
      node%drange(j,imin,2) = MAX( node%drange(j,imin,2), pdist(j) )
    END FORALL

    ! instead of recursive traversal, like
    ! CALL gnat_insert(tree%child(imin(1))%ptr, p, idx)
    ! we return with status lcomplete == .FALSE.
    tree = node%child(imin(1))
    IF (tree == UNASSOCIATED) THEN
      ! use the free node provided by the calling
      ! function to store the new node child:
      node%child(imin(1)) = free_node
      free_node = UNASSOCIATED
      tree = node%child(imin(1))
      node => node_storage(tree)
      node%isplit_pts    = 0
      node%drange(:,:,1) = MAX_RANGE
      node%drange(:,:,2) = 0._gk
      node%child(:) = UNASSOCIATED
    END IF

    lcomplete = .FALSE.
    
  END SUBROUTINE gnat_insert


  ! -----------------------------------------------------------------------------------
  !>  identify points within given radius
  RECURSIVE SUBROUTINE gnat_recursive_query(tree_idx, v, r, min_dist, &
    &                                       min_node_idx, count_dist)
    INTEGER,           INTENT(IN)    :: tree_idx ! index pointer
    REAL(gk),          INTENT(IN)    :: v(icoord_dim)
    REAL(gk),          INTENT(IN)    :: r
    REAL(gk),          INTENT(INOUT) :: min_dist
    INTEGER,           INTENT(INOUT) :: min_node_idx(2)
    INTEGER,           INTENT(INOUT) :: count_dist ! counter: distance calculations
    ! local variables
    LOGICAL                          :: pflag(gnat_k)
    INTEGER                          :: ip, isplit_pts
    REAL(gk)                         :: dist_vp
    TYPE(t_gnat_tree), POINTER       :: tree
#ifdef __SX__
    REAL(gk)                         :: pdist(gnat_k)
#endif
    
    tree => node_storage(tree_idx)
    isplit_pts = tree%isplit_pts
    IF (isplit_pts == 0) RETURN

    ! include all split points in set P
    ! note that values for indices > tree%isplit_pts are undefined!
    pflag(1:isplit_pts) = .TRUE.

    ! note: On the SX it might be cheaper to compute all distances (vectorized)
    ! outside of the following loop (though many won't be needed then)
#ifdef __SX__
    CALL dist_vect(tree%p, v, isplit_pts, pdist)
    count_dist = count_dist + isplit_pts
#endif
    
    ! remove some elements from P
    P : DO ip=1,isplit_pts
      IF (pflag(ip)) THEN

#ifdef __SX__
        dist_vp = pdist(ip)
#else
        dist_vp = dist(tree%p(ip), v)
        count_dist = count_dist + 1
#endif
        IF ((dist_vp <= r) .AND. (dist_vp < min_dist)) THEN
          min_dist = dist_vp
          min_node_idx(1:2) = tree%p(ip)%idx(1:2)
        END IF

        WHERE (pflag(1:isplit_pts))
          pflag(1:isplit_pts) = (tree%drange(ip,1:isplit_pts,1) <= (dist_vp+r)) .AND. &
            &                   (tree%drange(ip,1:isplit_pts,2) >= (dist_vp-r))
        END WHERE

      END IF
    END DO P

    ! traverse subtrees for remaining ip in P
    DO ip=1,isplit_pts
      IF (pflag(ip) .AND. (tree%child(ip) /= UNASSOCIATED)) THEN
        CALL gnat_recursive_query(tree%child(ip), v, r, min_dist, &
          &                       min_node_idx(:), count_dist)
      END IF
    END DO

  END SUBROUTINE gnat_recursive_query


  ! -----------------------------------------------------------------------------------
  !> queries a sequence of points, exploiting previous search results
  SUBROUTINE gnat_query_list(tree, v, iv_nproma, iv_nblks, iv_npromz, iv, istart, &
    &                        rr, min_dist, min_node_idx, count_dist)

    integer,       INTENT(IN)    :: tree                 ! current tree node (index)
    INTEGER,       INTENT(IN)    :: iv_nproma, iv_nblks, iv_npromz     ! list size
    INTEGER,       INTENT(IN)    :: iv                   ! local chunk size (blocks)
    INTEGER,       INTENT(IN)    :: istart               ! start idx
    REAL(gk),      INTENT(IN)    :: v(iv_nproma, iv_nblks, icoord_dim)    ! search point
    REAL(gk),      INTENT(IN)    :: rr                   ! initial search radius
    REAL(gk),      INTENT(INOUT) :: min_dist(iv_nproma, iv_nblks)         ! minimal distance
    INTEGER,       INTENT(INOUT) :: min_node_idx(iv_nproma, iv_nblks,2)   ! corresponding node
    INTEGER,       INTENT(INOUT) :: count_dist           ! counter: distance calculations
    ! local parameters
    INTEGER                      :: jb, jc, end_idx      ! block, index loop counter
    REAL(gk)                     :: r2, r, min_dist_old
    REAL(gk)                     :: p_old(icoord_dim), p_new(icoord_dim)
    REAL(gk)                     :: vmin_dist
    INTEGER                      :: vmin_node_idx(2)

    ! Description:

    ! Often, proximity search for a sequence of points means that
    ! subsequent points are located next to each other. In this case we can
    ! shorten the tree traversal (and reduce the number of distance
    ! computations) by the following strategies:
    !    the search result from a previous run gives us a true upper
    !    bound for the new search radius that may be smaller than the
    !    previous search radius:
    !    radius_new = min(radius_old, min_dist_old + |v_old, v_new|)

    r = rr
    min_dist(:,istart:istart+iv-1)       = MAX_RANGE
    min_node_idx(:,istart:istart+iv-1,1) = INVALID_NODE

    DO jb=istart, istart+iv-1

      ! set end index in current block:
      end_idx = iv_nproma
      IF (jb == iv_nblks) THEN
        end_idx = iv_npromz
        min_node_idx((end_idx+1):, jb, 1) = UNUSED_NODE
      END IF
      
      DO jc=1,end_idx

        p_new(1:icoord_dim) = v(jc,jb,1:icoord_dim)

        ! first query, block "istart", index "1": no previous result
        ! available, otherwise use the previous search result
        IF ((jb > istart) .OR. (jc>1)) THEN

          ! distance between old and new search point:
          r2 = dist_p(p_old(:), p_new(:))
          count_dist = count_dist + 1

          ! compute a true upper bound for search radius
          ! ("1.05" is just for safety)
          r = MINVAL((/ rr, 1.05_gk*(r2 + min_dist_old) /))
          
        END IF

        ! traverse tree
        vmin_dist        = MAX_RANGE
        vmin_node_idx(1) = INVALID_NODE
        CALL gnat_recursive_query(tree, p_new(:), r,           &
          &                       vmin_dist, vmin_node_idx(:), &
          &                       count_dist)
        min_dist(jc,jb)         = vmin_dist
        min_node_idx(jc,jb,1:2) = vmin_node_idx(1:2)

        p_old(1:icoord_dim)     = p_new(1:icoord_dim)
        min_dist_old = vmin_dist

      END DO
    END DO
    
  END SUBROUTINE gnat_query_list


  ! -----------------------------------------------------------------------------------
  !> multi-threaded insertion of points into GNAT data structure
  !
  ! @note The construction of the GNAT data structure is aimed at
  !       interpolation purposes. Thus we insert only cells with
  !       "ref_ctrl" flags >= 2.

  SUBROUTINE gnat_insert_mt(p_patch, nproc, count_dist)

    TYPE(t_patch), INTENT(IN)  :: p_patch
    INTEGER,       INTENT(IN)  :: nproc
    INTEGER,       INTENT(OUT) :: count_dist ! counter: distance calculations

    ! Local parameters
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_gnat_gridsearch:gnat_insert_mt")
    INTEGER                    :: root ! index pointer 
    INTEGER                    :: idx(2,nproc), work(nproc), depth(nproc)
    INTEGER                    :: node_proc(nproc) ! index pointer
    INTEGER                    :: free_node(nproc) ! short list of available nodes
    REAL(gk)                   :: p(icoord_dim,nproc)
    INTEGER                    :: j, new_idx(2)
    INTEGER                    :: iproc
    LOGICAL                    :: lcomplete
    TYPE (t_geographical_coordinates) :: pcoord
    TYPE (t_gnat_tree), POINTER:: node
    INTEGER                    :: i_startblk, i_endblk, &
      &                           i_startidx, i_endidx, &
      &                           rl_start, rl_end, i_nchdom
    LOGICAL                    :: l_loop_end

    IF (dbg_level > 10) THEN
        WRITE(message_text,*) "Running with ", nproc, "thread(s)."
        CALL message(routine, TRIM(message_text))
    END IF

    ! values for the blocking
    rl_start = 2
    rl_end = min_rlcell_int
    
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)    

    ! initialize cell counter
    CALL get_indices_c(p_patch, i_startblk,  &
      &                i_startblk, i_endblk, &
      &                i_startidx, i_endidx, &
      &                rl_start, rl_end)
    new_idx    =  (/ i_startidx, i_startblk /)

    count_dist = 1

    ! create a root node
    CALL resize_node_storage()
    root = get_free_node()
    node => node_storage(root)
    node%isplit_pts    = 0
    node%drange(:,:,1) = MAX_RANGE
    node%drange(:,:,2) = 0._gk
    node%child(:) = UNASSOCIATED

    idx (:,:)    = -1  ! indices of nodes processed in parallel
    work(:)      =  0  ! counter: no. of pts. inserted by procs

    depth(:)     =  0
    tree_depth   =  0  ! maximum tree depth
    free_node(:) = UNASSOCIATED
    
    ! Short description of parallel processing:

    ! 0 mark all processes as "idle" (idx == -1)
    ! 1 insert one point into the "pipeline" of points processed in parallel.
    ! 2 move all currently processed points one level down the tree
    !   this happens in parallel, since we know that all active
    !   nodes are on a different level of the GNAT structure.
    ! 3 if one point has reached its final position in the
    !   GNAT structure, mark the corresponding process as "idle" again
    ! 4 repeat steps 1-3 until all points have been inserted

    ! disadvantage of this simple pipelining: for a moderate number
    ! of nodes in the GNAT, the tree is rather flat, thus processes
    ! may be often idle since all levels are currently occupied.
    ! In other words: This algorithm does not scale well beyond ~8 threads!

    count_dist = 0
    l_loop_end = .FALSE.

!$OMP PARALLEL num_threads(nproc),                    &
!$OMP          private(lcomplete,iproc),              &
!$OMP          reduction(+:count_dist)

      iproc = 1
!$    iproc = OMP_GET_THREAD_NUM() + 1

    POINTS : DO

!$OMP BARRIER
!$OMP MASTER
      !-- block executed by first thread:
      !-- resize storage for tree nodes if necessary
      IF (num_nodes >= max_num_nodes-2*nproc) THEN
        CALL resize_node_storage()
      END IF

      !-- build a short list of available node
      !   indices for NEW tree nodes
      DO j=1,nproc
        IF (free_node(j) == UNASSOCIATED) THEN
          free_node(j) = get_free_node()
        END IF
      END DO

      !-- assign new work to one of the PEs
      ! find first non-associated PE in node_proc
      j = 1
      FINDFIRST : DO 
        IF (j > nproc)      EXIT FINDFIRST
        IF (idx(1,j) == -1) EXIT FINDFIRST
        j = j + 1
      END DO FINDFIRST

      IF (j<=nproc) THEN
        ! insert next point for this proc
        IF (new_idx(2) <= i_endblk) THEN
          new_idx(1) = new_idx(1) + 1   ! increase index
          IF (new_idx(1) > i_endidx) THEN
            new_idx(2) = new_idx(2) + 1 ! increase block

            CALL get_indices_c(p_patch, new_idx(2),  &
              &                i_startblk, i_endblk, &
              &                i_startidx, i_endidx, &
              &                rl_start, rl_end)

            new_idx(1) = i_startidx
          END IF
          
          IF ((new_idx(2) <= i_endblk) .AND.  &
            & (new_idx(1) <= i_endidx)) THEN
            pcoord       = p_patch%cells%center(new_idx(1), new_idx(2))
            p(1,j)       = REAL( pcoord%lon, gk)
            p(2,j)       = REAL( pcoord%lat, gk)
            idx(1:2,j)   = new_idx(1:2)
            work(j)      = work(j) + 1
            tree_depth   = MAX(tree_depth, depth(j))
            depth(j)     = 0
            node_proc(j) = root
          END IF
        END IF
      END IF
      IF ((new_idx(2) > i_endblk) .AND. (ALL(idx(1,:) == -1))) THEN
        !EXIT POINTS
        l_loop_end = .TRUE.
      END IF
!$OMP END MASTER
!$OMP BARRIER

      IF (l_loop_end) EXIT points

      ! step all active nodes one step down the tree
      IF (idx(1,iproc) /= -1) THEN
        CALL gnat_insert(node_proc(iproc), p(:,iproc), idx(:,iproc), &
          &              free_node(iproc), lcomplete, count_dist)
        ! "lcomplete" : proc is free for new point
        IF (lcomplete) THEN
          idx(1,iproc) = -1
        ELSE
          depth(iproc) = depth(iproc) + 1
        END IF
      END IF
    END DO POINTS
!$OMP END PARALLEL

    gnat_tree = root
    IF (dbg_level > 10) THEN
      WRITE(message_text,*) "distribution of insertions over threads: ", work(1:nproc), &
        &                   ", maximum tree depth:", tree_depth
      CALL message(routine, TRIM(message_text))
    END IF

  END SUBROUTINE gnat_insert_mt


  ! -----------------------------------------------------------------------------------
  !> multi-threaded query for a sequence of points
  ! note: 
  ! - the search radius parameter is optional
  ! - search radius is adapted in case of algorithm failure
  SUBROUTINE gnat_query_list_mt(tree, v, iv_nproma, iv_nblks, iv_npromz,                   &
    &                           min_dist, min_node_idx, count_dist, ladapt_radius, opt_rr)

    ! max. trial steps to adapt search radius, r_new = opt_rr*10**iadapt
    INTEGER, PARAMETER :: nadapt = 15

    INTEGER,            INTENT(IN)    :: tree                 ! tree root node (index)
    INTEGER,            INTENT(IN)    :: iv_nproma, iv_nblks, iv_npromz     ! list size
    REAL(gk),           INTENT(IN)    :: v(iv_nproma, iv_nblks, icoord_dim) ! list of search points
    REAL(gk),           INTENT(INOUT) :: min_dist(iv_nproma, iv_nblks)      ! minimal distance
    INTEGER,            INTENT(INOUT) :: min_node_idx(iv_nproma, iv_nblks,2)! corresponding node
    INTEGER,            INTENT(INOUT) :: count_dist           ! counter: distance calculations
    LOGICAL,            INTENT(IN)    :: ladapt_radius
    REAL(gk), OPTIONAL, INTENT(IN)    :: opt_rr                   ! opt. search radius
    ! local parameters
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_gnat_gridsearch:gnat_query_list_mt")
    INTEGER                           :: istart, ichunksize, iend,     &
    &                                    nproc, iproc, iadapt, jb, jc, &
    &                                    end_idx, nqueries
    REAL(gk)                          :: r                    ! search radius
    REAL(gk)                          :: vmin_dist
    INTEGER                           :: vmin_node_idx(2)

    IF (PRESENT(opt_rr)) THEN
      r = opt_rr
    ELSE
      ! determine a meaningful search radius
      r = gnat_std_radius(tree) ! search radius ~ (query time)**(-1)
    END IF
    count_dist = 0 ! reset performance counter

    nproc = 1
!$  nproc = OMP_GET_MAX_THREADS()

!$OMP PARALLEL private(istart, iend, ichunksize, iproc), &
!$OMP          reduction(+:count_dist)
    
    iproc      = 1
!$  iproc      = OMP_GET_THREAD_NUM() + 1

    ! proc "iproc" will query the entry blocks (istart,...,iend)
    ichunksize = (iv_nblks+nproc-1)/nproc
    istart     = (iproc-1)*ichunksize + 1
    iend       = MIN(iv_nblks, istart+ichunksize-1)

    IF (iend >= istart) THEN
      ichunksize = iend - istart + 1
      CALL gnat_query_list(tree, v, iv_nproma, iv_nblks, iv_npromz, ichunksize, &
        &                  istart, r, min_dist, min_node_idx, count_dist)
    END IF

!$OMP END PARALLEL

    ! for all failed queries (if there are any) repeat query with a
    ! larger search radius:
    IF (.NOT. ladapt_radius) RETURN

    iadapt = 0
    DO
      IF ((iadapt > nadapt) .OR. &
        & (COUNT(min_node_idx(:,:,1) == INVALID_NODE) == 0)) EXIT
      ! adapt search radius
      iadapt = iadapt + 1
      r = r*2._gk
      IF (r > MAX_RANGE) EXIT
      IF (dbg_level > 1) &
        CALL message(routine, "adapting radius")

!$OMP PARALLEL DO private(end_idx,jb,jc,vmin_dist, vmin_node_idx),  &
!$OMP             reduction(+:count_dist)
      DO jb=1,iv_nblks

        ! set end index in current block:
        end_idx = iv_nproma
        IF (jb == iv_nblks) end_idx = iv_npromz

        DO jc=1,end_idx
          IF (min_node_idx(jc,jb,1) == INVALID_NODE) THEN
            vmin_dist        = MAX_RANGE
            vmin_node_idx(1) = INVALID_NODE
            CALL gnat_recursive_query(tree, v(jc,jb,:), r,     &
              &                       vmin_dist, vmin_node_idx(:), count_dist)
            min_dist(jc,jb)         = vmin_dist
            min_node_idx(jc,jb,1:2) = vmin_node_idx(1:2)
          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
    END DO

    IF (iadapt > nadapt) THEN
      WRITE(message_text,*) "failed to adapt search radius! Last radius=",r
      CALL message(routine, TRIM(message_text))
      ! CALL finish (routine, TRIM(message_text))
    END IF

    IF (dbg_level > 10) THEN
      IF (iadapt <= nadapt) THEN
        WRITE(message_text,*) "no. of adaptations: ", iadapt
        CALL message(routine, TRIM(message_text))
      END IF
      nqueries = (iv_nblks-1)*iv_nproma + iv_npromz
      WRITE(message_text,*) "no. of distance computations for tree query: ", count_dist
      CALL message(routine, TRIM(message_text))
      WRITE(message_text,*) "  no. of queries: ", nqueries
      CALL message(routine, TRIM(message_text))
      WRITE(message_text,*) "  i.e. ", count_dist/nqueries, " computations/query."
      CALL message(routine, TRIM(message_text))
    END IF

  END SUBROUTINE gnat_query_list_mt


  ! -----------------------------------------------------------------------------------
  !> determine a meaningful value for the search radius
  FUNCTION gnat_std_radius(tree_idx)

    REAL(gk)                                :: gnat_std_radius
    INTEGER,                    INTENT(IN)  :: tree_idx ! current tree node (index)
    ! local parameters
    TYPE (t_gnat_tree), POINTER             :: node
    INTEGER                                 :: i

    ! move down tree an collect range values
    node => node_storage(tree_idx)
    i = node%isplit_pts
    gnat_std_radius = MINVAL(node%drange(1, 2:i, 1))
    DO
      i = node%isplit_pts
      gnat_std_radius = MIN(gnat_std_radius, MINVAL(node%drange(1, 2:i, 1)))

      IF (node%child(1) == UNASSOCIATED) EXIT
      node => node_storage(node%child(1))
    END DO
    ! increase radius for safety
    gnat_std_radius = 2._gk*gnat_std_radius
  END FUNCTION gnat_std_radius


  ! -----------------------------------------------------------------------------------
  !> reads ICON triangle center points into GNAT data structure
  ! (Possible) improvement for the future:
  ! - implementation of min/max range boundaries if not
  !   all grid points are of interest.
  SUBROUTINE gnat_init_grid(p_patch)

    TYPE(t_patch), INTENT(IN)    :: p_patch
    ! local variables:
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_gnat_gridsearch:gnat_init_grid")
    INTEGER                      :: count_dist, nproc

    nproc = 1
!$  nproc = OMP_GET_MAX_THREADS()

    ! consistency check
    IF (gnat_tree /= UNASSOCIATED) THEN
      CALL message(routine, "Discarding existing data in GNAT structure!")
      CALL gnat_destroy()
    END IF

    ! tell the GNAT how many nodes are expected
    ! ("5" is for safety reasons)
    expected_num_nodes = 5*p_patch%n_patch_cells/gnat_k + 1

    ! build GNAT data structure based on clon, clat:
    CALL gnat_insert_mt(p_patch, nproc, count_dist)

  END SUBROUTINE gnat_init_grid


  ! -----------------------------------------------------------------------------------
  !>  destroy GNAT data structure
  RECURSIVE SUBROUTINE gnat_destroy()
    INTEGER :: errstat

    IF (max_num_nodes > 0) THEN
      DEALLOCATE(node_storage, STAT=errstat)
      IF (errstat /= 0)  &
          CALL finish ('mo_gnat_gridsearch:gnat_destroy', &
          &            'Error in DEALLOCATE operation!')
      num_nodes     = 0
      max_num_nodes = 0
      gnat_tree     = UNASSOCIATED
    END IF

  END SUBROUTINE gnat_destroy


  ! -----------------------------------------------------------------------------------
  !> test for counter-clockwise direction when travelling 
  !  from point p0 to p1 to p2.
  ! @return -1 if clockwise direction or         
  !            if pts are collinear and p0 is between p1 and p2
  ! @return  0 if pts are collinear and p2 is between p0 and p1
  ! @return +1 if counter-clockwise direction or 
  !            if pts are collinear and p1 is between p0 and p2
  SUBROUTINE ccw(p0, p1, p2, ccw_result)
    
    INTEGER, INTENT(OUT)               :: ccw_result
    REAL(gk), DIMENSION(2), INTENT(IN) :: p0, p1, p2
    ! local parameters
    REAL(gk)                           :: d1(2), d2(2), v1, v2

    d1 = p1 - p0
    d2 = p2 - p0

    v1 = d1(1)*d2(2)
    v2 = d1(2)*d2(1)
    IF (v1 > v2) THEN
      ccw_result = +1
      RETURN
    ELSE IF (v1 < v2) THEN
      ccw_result = -1
      RETURN
    END IF

    IF ((d1(1)*d2(1) < 0._wp) .OR. (d1(2)*d2(2) < 0._wp)) THEN
      ccw_result = -1
      RETURN
    ELSE 
      v1 = d1(1)*d1(1) + d1(2)*d1(2)
      v2 = d2(1)*d2(1) + d2(2)*d2(2)
      IF (v1 < v2) THEN
        ccw_result = +1
        RETURN
      END IF
    END IF
    ccw_result = 0
  END SUBROUTINE ccw


  ! -----------------------------------------------------------------------------------
  !> Simple geometric test for "point inside triangle"
  ! See, e.g., the chapter on elementary geometric methods in
  ! Sedgewick, Robert: Algorithms in C. Reading, MA: Addison Wesley Longman.
  SUBROUTINE point_inside_triangle(v, p, flag)
    
    LOGICAL, INTENT(OUT) :: flag
    REAL(gk), INTENT(IN) :: v(2,3), p(2)
    INTEGER              :: ccw1, ccw2, ccw3

    CALL ccw(v(:,1), v(:,2), p, ccw1)
    CALL ccw(v(:,2), v(:,3), p, ccw2)
    CALL ccw(v(:,3), v(:,1), p, ccw3)

    flag = ((ccw1 >= 0)  .AND.  (ccw2 >= 0)  .AND.  (ccw3 >= 0)) .OR. &
      &    ((ccw1 <= 0)  .AND.  (ccw2 <= 0)  .AND.  (ccw3 <= 0))

  END SUBROUTINE point_inside_triangle


  ! -----------------------------------------------------------------------------------
  ! In this routine, the master process merges the results of several
  ! distributed nearest neighbor searches, i.e. we reduce multiple
  ! lists of pairs (distance, index) of the same length to a single
  ! list. All information except the one with minimal distance is
  ! discarded.
  SUBROUTINE gnat_merge_distributed_queries(p_patch, total_dim, iv_nproma, iv_nblks, min_dist, &
    &                                       tri_idx, lonlat_points, global_idx,                &
    &                                       ithis_local_pts)
    TYPE(t_patch),         INTENT(IN)    :: p_patch
    INTEGER,               INTENT(IN)    :: total_dim
    INTEGER,               INTENT(IN)    :: iv_nproma, iv_nblks
    REAL(gk),              INTENT(IN)    :: min_dist(:,:)
    INTEGER,               INTENT(INOUT) :: tri_idx(:,:,:)
    REAL(gk),              INTENT(INOUT) :: lonlat_points(:,:,:)
    INTEGER,               INTENT(INOUT) :: global_idx(:)
    INTEGER,               INTENT(OUT)   :: ithis_local_pts !< no. of points on this PE
    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_gnat_gridsearch:gnat_merge_distributed_queries")
    INTEGER  :: i, j, jc, jb
    INTEGER  :: new_tri_idx(2,total_dim)
    REAL(gk) :: new_lonlat_points(total_dim,2)
    REAL     :: in(2,total_dim)
    INTEGER  :: iowner, idummy_applied
    INTEGER  :: array_shape(2), dummy_idx(2)
    INTEGER  :: i_endblk, i_endidx, &
      &         rl_start, rl_end, i_nchdom, my_id

    my_id = get_my_mpi_work_id()

    ! for the (rather pathological) case of local patches covered by a
    ! larger lon-lat grid, we have to set some "dummy" indices outside
    ! the domain:
    rl_start = 2
    rl_end = min_rlcell_int
    i_nchdom     = MAX(1,p_patch%n_childdom)
    dummy_idx(2) = p_patch%cells%start_blk(rl_start,1)
    i_endblk     = p_patch%cells%end_blk(rl_end,i_nchdom)    
    CALL get_indices_c(p_patch, dummy_idx(2),  &
      &                dummy_idx(2), i_endblk, &
      &                dummy_idx(1), i_endidx, &
      &                rl_start, rl_end)

    ! 1. Perform an MPI_ALLREDUCE operation with the min_dist vector
    !    to find out which process is in charge of which in_point

    in(1,:) = RESHAPE(REAL(min_dist(:,:)), (/ total_dim /) )
    in(2,:) = REAL( my_id )
    IF (p_patch%n_patch_cells == 0) in(2,:) = -1.;
    
    CALL p_allreduce_minloc(in, total_dim, p_comm_work)

    ! 2. If we are a working PE, reduce the list of in_points and the
    !    tri_idx array to points which are actually located on this
    !    portion of the domain.

    ! store list of global indices owned by this PE:
    j=0
    DO i=1,total_dim
      iowner = NINT(in(2,i))
      IF (iowner == my_id) THEN
        j = j + 1 
        global_idx(j) = i
      END IF
    END DO
    ! now, j points are left for proc "get_my_mpi_work_id()"
    ithis_local_pts = j

    idummy_applied = 0
    DO i=1,ithis_local_pts
      ! convert global index into local idx/block pair:
      j = global_idx(i)
      jb = (j-1)/iv_nproma + 1
      jc = j - (jb-1)*iv_nproma

      IF (tri_idx(1,jc,jb) /= INVALID_NODE) THEN
        new_tri_idx(1:2,i) = tri_idx(1:2,jc,jb)          
      ELSE
        new_tri_idx(1:2,i) = dummy_idx
        idummy_applied = idummy_applied + 1
      END IF
      new_lonlat_points(i,1:2) = lonlat_points(jc,jb,1:2)
    END DO
    
    array_shape(:) = (/ iv_nproma, iv_nblks /)
    tri_idx(1,:,:) = RESHAPE( new_tri_idx(1,:), array_shape(:), (/ 0 /) )
    tri_idx(2,:,:) = RESHAPE( new_tri_idx(2,:), array_shape(:), (/ 0 /) )
    lonlat_points(:,:,1) = RESHAPE( new_lonlat_points(:,1), array_shape(:), (/ 0._gk /) )
    lonlat_points(:,:,2) = RESHAPE( new_lonlat_points(:,2), array_shape(:), (/ 0._gk /) )

    IF (dbg_level > 10) THEN
      IF (idummy_applied > 0) THEN
        WRITE(message_text,*) "proc ", my_id, ": dummy applied: ", idummy_applied
        CALL message(routine, TRIM(message_text))
      END IF
    END IF
    
  END SUBROUTINE gnat_merge_distributed_queries


  ! -----------------------------------------------------------------------------------
  ! Performs a nearest-neighbor query for a given list of points and
  ! returns the indices and block indices of the mesh triangles that
  ! contain these points.
  SUBROUTINE gnat_query_containing_triangles(p_patch, tree, v, iv_nproma, iv_nblks,  &
    &                                        iv_npromz, tri_idx, min_dist)

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    INTEGER,  INTENT(IN)    :: tree                                ! tree root node (index)
    INTEGER,  INTENT(IN)    :: iv_nproma, iv_nblks, iv_npromz      ! list size
    REAL(gk), INTENT(IN)    :: v(iv_nproma, iv_nblks, icoord_dim)  ! list of search points
    INTEGER,  INTENT(OUT)   :: tri_idx(2,iv_nproma, iv_nblks)      ! containing triangle (idx,block)
    REAL(gk), INTENT(OUT)   :: min_dist(iv_nproma, iv_nblks)       ! minimal distance
    ! local parameters
    TYPE (t_grid_cells)   , POINTER  :: cells
    TYPE (t_grid_vertices), POINTER  :: verts
    INTEGER                     :: i_nv
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_gnat_gridsearch:gnat_query_containing_triangles")
    INTEGER                :: min_node_idx(iv_nproma, iv_nblks, 2)  ! corresponding GNAT nodes
    INTEGER                :: count_dist                            ! counter: distance calculations

    INTEGER                :: jb, jc, j, k, i_nb,  &
      &                       end_idx, i_end,      &
      &                       tmp_idx(2),          &   ! (idx,blk)
      &                       nb_idx(0:(5*3),2)        ! (idx,blk)
    INTEGER                :: tri_vertex_idx(3,2)   ! (idx,blk)
    REAL(gk)               :: tri_v(2,3), p(2), radius
    LOGICAL                :: l_inside, l_check_point_inside
    INTEGER                :: i_startblk, &
      &                       i_startidx, i_endidx, &
      &                       rl_start, rl_end, i_nchdom

    ! set default value ("failure notice")
    tri_idx(1,:,:) = INVALID_NODE
    min_dist(:,:)  = MAX_RANGE

    IF (p_patch%n_patch_cells == 0) RETURN;

    ! TODO[FP] : for the time being, we find it sufficiently accurate to
    !            perform a simple nearest neighbor query.
    l_check_point_inside = .FALSE.

    cells => p_patch%cells
    verts => p_patch%verts
    i_nv  =  p_patch%cell_type

    IF (i_nv /= 3) THEN
      CALL finish(routine, "Wrong number of cell vertices!")
    END IF

    ! make a sensible guess for search radius (just taking a "randomly
    ! chosen" edge length)
    rl_start = 2
    rl_end   = min_rlcell_int
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    CALL get_indices_e(p_patch, i_startblk,  &
      &                i_startblk, i_startblk, &
      &                i_startidx, i_endidx, &
      &                rl_start, rl_end)
    radius = 3._gk * REAL(p_patch%edges%primal_edge_length(i_startidx,i_startblk) &
      & /grid_sphere_radius, gk)

    ! query list of nearest neighbors
    ! TODO[FP] : For some test cases it might be reasonable to enable
    ! radius adaptation
    count_dist = 0
    CALL gnat_query_list_mt(tree, v, iv_nproma, iv_nblks, iv_npromz, &
      &                     min_dist, min_node_idx, count_dist, .FALSE., radius)

    ! loop over blocks
!$OMP PARALLEL DO private(jb, end_idx, jc, nb_idx, tri_vertex_idx, &
!$OMP                     i_nb, tmp_idx, i_end, tri_v, k, j, p, l_inside)
    DO jb=1,iv_nblks

      ! set end index in current block:
      end_idx = iv_nproma
      if (jb == iv_nblks) end_idx = iv_npromz

      ! loop over indices
      DO jc=1,end_idx

        CHECK_INSIDE : IF (.NOT. l_check_point_inside) THEN

          tri_idx(1:2,jc,jb) = min_node_idx(jc,jb,1:2)

        ELSE

          p = REAL( v(jc, jb, 1:2), gk ) ! current search point

          IF (dbg_level > 10) THEN
            WRITE(message_text,*) "jb,jc = ", jb,",",jc, ": (", v(jc,jb,:), ")"
            CALL message(routine, TRIM(message_text))
          END IF

          ! for each nearest neighbor, retrieve local neighborhood
          ! consisting of the neighbors of all 3 vertices
          ! (without the center cell this makes 15 cells instead of 12
          !  since duplicates are not removed from the set)
          nb_idx(0,1:2) = min_node_idx(jc,jb,1:2)
          IF (nb_idx(0,1) == INVALID_NODE) CYCLE

          ! get all vertices of center cell
          tri_vertex_idx(1:3,1) = cells%vertex_idx(nb_idx(0,1), nb_idx(0,2), 1:3)
          tri_vertex_idx(1:3,2) = cells%vertex_blk(nb_idx(0,1), nb_idx(0,2), 1:3)
          ! get all cells adjacent to these vertices:
          i_nb = 1
          VERTEX : DO j=1,3 ! i_nv
            CELLS_OF_VERTEX : DO k=1,6
              tmp_idx(1) = verts%cell_idx(tri_vertex_idx(j,1), tri_vertex_idx(j,2), k)
              tmp_idx(2) = verts%cell_blk(tri_vertex_idx(j,1), tri_vertex_idx(j,2), k)

              if (tmp_idx(1) == 0) CYCLE CELLS_OF_VERTEX

              ! don't store the central cell of stencil;
              ! also, don't include neighbor with refin_ctrl=1, because
              ! FD stencil is not available.
              IF ((cells%refin_ctrl(tmp_idx(1), tmp_idx(2)) > 1) .AND.  &
                & (tmp_idx(1) /= nb_idx(0,1)) .OR. (tmp_idx(2) /= nb_idx(0,2))) THEN
                nb_idx(i_nb,1:2) = tmp_idx(1:2)
                i_nb = i_nb + 1
              END IF
            END DO CELLS_OF_VERTEX
          END DO VERTEX

          ! loop over neighborhood
          i_end = i_nb-1
          LOOP : DO i_nb=0,i_end
            ! for each neighborhood collect triangle vertices
            tri_vertex_idx(1:3,1) = cells%vertex_idx(nb_idx(i_nb,1), nb_idx(i_nb,2), 1:3)
            tri_vertex_idx(1:3,2) = cells%vertex_blk(nb_idx(i_nb,1), nb_idx(i_nb,2), 1:3)
            FORALL (j=1:3) ! i_nv
              tri_v(1,j) = REAL( verts%vertex(tri_vertex_idx(j,1), tri_vertex_idx(j,2))%lon, gk)
              tri_v(2,j) = REAL( verts%vertex(tri_vertex_idx(j,1), tri_vertex_idx(j,2))%lat, gk)
            END FORALL

            IF (dbg_level > 10) THEN
              WRITE(message_text,*) "i_nb=", i_nb, ", vertices: ", tri_v(:,1), ";", &
                &                   tri_v(:,2), ";", tri_v(:,3)
              CALL message(routine, TRIM(message_text))
            END IF

            ! for each neighborhood find which triangle contains point
            CALL point_inside_triangle(tri_v(:,:), p(:), l_inside)
            IF (l_inside) THEN
              tri_idx(:,jc,jb) = nb_idx(i_nb,:)
              IF (dbg_level > 10) THEN
                WRITE(message_text,*) "Found containing triangle! Index(", i_nb, ") ", &
                  &                   tri_idx(:,jc,jb)
                CALL message(routine, TRIM(message_text))
              END IF
              EXIT LOOP
            END IF
          END DO LOOP

          ! failure: point was not contained in stencil
          IF (tri_idx(1,jc,jb) == -1) THEN
            IF (dbg_level > 1) &        
              CALL message(routine, &
              &      "resorting to neighbor guess (possibly gaps contained in the grid)")
            tri_idx(:,jc,jb) = nb_idx(0,:)
          END IF

        END IF CHECK_INSIDE

      END DO
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE gnat_query_containing_triangles

END MODULE mo_gnat_gridsearch
