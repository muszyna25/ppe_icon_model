!! Binomial heap data structures
!!
!! "A binomial heap is a heap similar to a binary heap but also 
!!  supports quickly merging two heaps." (Wikipedia/en)
!!
!! Literature:
!! - Vuillemin, CACM 21 (1978), 309-315
!! - Cormen, Leiserson, Rivest, Stein: Introduction To Algorithms, pp. 455
!!
!! This implementation is derived from a test program written in C by
!! Bjoern B. Brandenburg <bbb [at] cs.unc.edu>.
!!
!! Initial revision: F. Prill, DWD (2012-09-27)
!!
MODULE mo_util_binheap

  USE mo_kind, ONLY: wp
  IMPLICIT NONE

  INTEGER, PARAMETER   :: NOT_IN_HEAP  = -1, INVALID_NODE = -1

  !> single data item in heap
  TYPE t_heap_data
    SEQUENCE
    !> interpolation weights
    REAL(wp) :: wgt
    !> src indices (index/block)
    INTEGER  :: sidx, sblk
    !> dst indices (index/block)
    INTEGER  :: didx, dblk
  END TYPE t_heap_data

  !> heap node, pointing to a data item
  TYPE t_heap_node 
    SEQUENCE
    INTEGER           :: parent, next, child ! "pointer" indices
    INTEGER           :: degree
    TYPE(t_heap_data) :: rdata
  END TYPE t_heap_node

  !> global heap data structure
  TYPE t_heap
    INTEGER           :: head, min
  END TYPE t_heap

  !> Global heap node storage
  TYPE t_node_storage
    TYPE(t_heap_node), ALLOCATABLE :: v(:)
  END TYPE t_node_storage

  TYPE(t_node_storage), ALLOCATABLE, TARGET :: node_storage(:)
  INTEGER, ALLOCATABLE :: nnode(:)

  PRIVATE
  PUBLIC :: heap_init, heap_node_init
  PUBLIC :: heap_cmp
  PUBLIC :: node_storage_init, node_storage_finalize
  PUBLIC :: resize_node_storage
  PUBLIC :: get_free_node
  PUBLIC :: heap_empty
  PUBLIC :: heap_insert
  PUBLIC :: heap_union
  PUBLIC :: heap_min
  PUBLIC :: heap_take, heap_take_accumulated
  PUBLIC :: heap_add_offset
  PUBLIC :: INVALID_NODE
  ! data types and variables
  PUBLIC :: t_heap_data, t_heap_node, t_heap
  PUBLIC :: t_node_storage, node_storage, nnode

CONTAINS

  !> Generic comparison of two data items.
  !  @return "+1" if a<b, "0", if a==b and "-1" if a>b
  ELEMENTAL FUNCTION heap_cmp(a,b)
    INTEGER :: heap_cmp
    TYPE (t_heap_data), INTENT(IN) :: a,b

    heap_cmp = 0
    IF (a%dblk < b%dblk) THEN
      heap_cmp =  1
      RETURN
    ELSE IF (a%dblk > b%dblk) THEN
      heap_cmp = -1
      RETURN
    END IF
    IF (a%didx < b%didx) THEN
      heap_cmp =  1
      RETURN
    ELSE IF (a%didx > b%didx) THEN
      heap_cmp = -1
      RETURN
    END IF
    IF (a%sblk < b%sblk) THEN
      heap_cmp =  1
      RETURN
    ELSE IF (a%sblk > b%sblk) THEN
      heap_cmp = -1
      RETURN
    END IF
    IF (a%sidx < b%sidx) THEN
      heap_cmp =  1
      RETURN
    ELSE IF (a%sidx > b%sidx) THEN
      heap_cmp = -1
      RETURN
    END IF
  END FUNCTION heap_cmp


  !> Utility function: allocates space for tree nodes
  SUBROUTINE resize_node_storage(expected_size, ithrd)
    INTEGER, INTENT(IN) :: expected_size, ithrd
    ! local variables
    INTEGER :: max_num_nodes, new_max_num_nodes, i
    TYPE (t_heap_node), TARGET, ALLOCATABLE :: tmp(:)

    ! triangle copy
    IF (ALLOCATED(node_storage(ithrd)%v)) THEN
      max_num_nodes = SIZE(node_storage(ithrd)%v)
    ELSE
      max_num_nodes = 0
    ENDIF
    new_max_num_nodes = max_num_nodes + expected_size
    IF (max_num_nodes > 0) THEN
      ALLOCATE(tmp(max_num_nodes))
      tmp(1:max_num_nodes) = node_storage(ithrd)%v(1:max_num_nodes)
      DEALLOCATE(node_storage(ithrd)%v)
    END IF
    ALLOCATE(node_storage(ithrd)%v(new_max_num_nodes))
    DO i=(max_num_nodes+1),new_max_num_nodes
!CDIR IEXPAND
      node_storage(ithrd)%v(i) = heap_node_init(t_heap_data(0._wp,0,0,0,0))
    END DO
    IF (max_num_nodes > 0) THEN
      node_storage(ithrd)%v(1:max_num_nodes) = tmp(1:max_num_nodes)
      DEALLOCATE(tmp)
    END IF
  END SUBROUTINE resize_node_storage


  FUNCTION get_free_node(ithrd)
    INTEGER :: get_free_node
    INTEGER, INTENT(IN) :: ithrd

    nnode(ithrd) = nnode(ithrd) + 1
    IF (nnode(ithrd) > SIZE(node_storage(ithrd)%v)) &
      &  CALL resize_node_storage(SIZE(node_storage(ithrd)%v), ithrd)
    get_free_node = nnode(ithrd)
  END FUNCTION get_free_node


  !> Initializes heap node storage.
  SUBROUTINE node_storage_init(expected_size, ithrd)
    INTEGER, INTENT(IN) :: expected_size, ithrd
    CALL resize_node_storage(expected_size, ithrd)
    nnode(ithrd) = 0
  END SUBROUTINE node_storage_init


  !> Frees heap data structure.
  SUBROUTINE node_storage_finalize(ithrd)
    INTEGER, INTENT(IN) :: ithrd
    DEALLOCATE(node_storage(ithrd)%v)
    nnode(ithrd) = 0
  END SUBROUTINE node_storage_finalize


  !> Initializes empty heap data structure.
  SUBROUTINE heap_init(heap)
    TYPE(t_heap), INTENT(INOUT) :: heap
    heap = t_heap(INVALID_NODE, INVALID_NODE)
  END SUBROUTINE heap_init


  !> Initializes emtpy heap node.
  FUNCTION heap_node_init(rdata) RESULT(h)
    TYPE(t_heap_node) :: h
    TYPE(t_heap_data), INTENT(IN) :: rdata

    h%parent = INVALID_NODE
    h%next   = INVALID_NODE
    h%child  = INVALID_NODE
    h%degree = NOT_IN_HEAP
    h%rdata  = rdata
  END FUNCTION heap_node_init


  !> @return .TRUE. if heap data structure is empty.
  FUNCTION heap_empty(heap)
    LOGICAL :: heap_empty
    TYPE(t_heap), INTENT(INOUT) :: heap

    heap_empty = (heap%head == INVALID_NODE) .AND. &
      &       (heap%min  == INVALID_NODE)
  END FUNCTION heap_empty

  !> Merges two heaps (internal routine)
  FUNCTION heap_merge(a_in,b_in, ithrd) RESULT(head)
    INTEGER :: head
    INTEGER, INTENT(IN) :: a_in,b_in, ithrd
    ! local variables:
    INTEGER :: pos,a,b, deg_a, deg_b
    TYPE(t_heap_node), POINTER :: nodes(:)

    nodes => node_storage(ithrd)%v
    a     = a_in
    deg_a = nodes(a)%degree
    b     = b_in
    deg_b = nodes(b)%degree

    head = INVALID_NODE
    pos = head
    IF ((a /= INVALID_NODE) .AND. (b /= INVALID_NODE)) THEN
      DO
        IF (deg_a < deg_b) THEN
          IF (pos == INVALID_NODE) THEN
            head = a
            pos  = head
          ELSE
            nodes(pos)%next = a
            pos = nodes(pos)%next
          END IF
          a     = nodes(a)%next
          IF (a == INVALID_NODE) EXIT
          deg_a = nodes(a)%degree
        ELSE
          IF (pos == INVALID_NODE) THEN
            head = b
            pos  = head
          ELSE
            nodes(pos)%next = b
            pos = nodes(pos)%next
          END IF
          b     = nodes(b)%next
          IF (b == INVALID_NODE) EXIT
          deg_b = nodes(b)%degree
        END IF
      END DO
    END IF
    IF (a /= INVALID_NODE) THEN
      nodes(pos)%next = a
    ELSE
      nodes(pos)%next = b
    END IF
  END FUNCTION heap_merge


  !> reverse a linked list of nodes. also clears parent pointer (internal routine)
  FUNCTION heap_reverse(h,ithrd) RESULT(res)
    INTEGER :: res
    INTEGER, INTENT(INOUT) :: h
    INTEGER, INTENT(IN)    :: ithrd
    ! local variables
    INTEGER :: tail, next

    tail = INVALID_NODE
    res  = h
    IF (h == INVALID_NODE)  RETURN
    node_storage(ithrd)%v(h)%parent = INVALID_NODE
    DO
      IF (node_storage(ithrd)%v(h)%next == INVALID_NODE) EXIT
      next                            = node_storage(ithrd)%v(h)%next
      node_storage(ithrd)%v(h)%next   = tail
      tail                            = h
      h                               = next
      node_storage(ithrd)%v(h)%parent = INVALID_NODE
    END DO
    node_storage(ithrd)%v(h)%next = tail
    res = h
  END FUNCTION heap_reverse


  SUBROUTINE heap_min(heap, prev, node, ithrd)
    TYPE(t_heap),  INTENT(IN)    :: heap
    INTEGER,       INTENT(INOUT) :: prev, node
    INTEGER,       INTENT(IN)    :: ithrd
    ! local variables
    INTEGER :: lprev, cur, hc
    TYPE (t_heap_node), POINTER :: nodes(:), pnode

    nodes => node_storage(ithrd)%v
    prev = INVALID_NODE
    node  = heap%head
    IF (node == INVALID_NODE)  RETURN
    pnode => nodes(node)
    lprev = node
    cur   = pnode%next
    DO
      IF (cur == INVALID_NODE) EXIT
!CDIR IEXPAND
      hc = heap_cmp(nodes(cur)%rdata, pnode%rdata)
      IF (hc > 0) THEN
        node  =  cur
        pnode => nodes(cur)
        prev  =  lprev
      END IF
      lprev = cur
      cur   = nodes(cur)%next
    END DO
  END SUBROUTINE heap_min


  SUBROUTINE local_heap_union(heap, h2, ithrd)
    TYPE(t_heap), INTENT(INOUT) :: heap
    INTEGER     , INTENT(IN)    :: h2, ithrd
    ! local variables
    INTEGER :: h1, prev, x, next, c
    TYPE(t_heap_node), POINTER :: nn, px, pnext, nodes(:)

    IF (h2 == INVALID_NODE)  RETURN
    h1 = heap%head
    IF (h1 == INVALID_NODE) THEN
      heap%head = h2
      RETURN
    END IF
    h1 = heap_merge(h1,h2,ithrd)
    prev = INVALID_NODE
    nodes => node_storage(ithrd)%v
    x     =  h1
    px    => nodes(x)
    next  =  px%next
    DO
      IF (next == INVALID_NODE) EXIT
      pnext => nodes(next)
      IF (px%degree /= pnext%degree) THEN
        ! nothing to do, advance
        prev = x
        x     =  next
        px    => pnext
        next  =  nodes(x)%next
        CYCLE
      END IF
      nn => nodes(next)
      IF (nn%next /= INVALID_NODE) THEN
        IF (nodes(nn%next)%degree == px%degree) THEN
          ! nothing to do, advance
          prev = x
          x     =  next
          px    => pnext
          next  =  nodes(x)%next
          CYCLE
        END IF
      END IF
!CDIR IEXPAND
      c = heap_cmp(px%rdata, pnext%rdata)
      IF (c > 0) THEN
        ! x becomes the root of next
        px%next = nn%next

        nodes(next)%parent = x
        nodes(next)%next   = nodes(x)%child
        nodes(x)%child     = next
        nodes(x)%degree    = nodes(x)%degree+1
      ELSE 
        ! next becomes the root of x
        IF (prev /= INVALID_NODE) THEN
          nodes(prev)%next = next
        ELSE
          h1 = next
        END IF
        nodes(x)%parent     = next
        nodes(x)%next       = nodes(next)%child
        nodes(next)%child   = x
        nodes(next)%degree  = nodes(next)%degree+1

        x   = next
        px  => nodes(x)
      END IF
      next = nodes(x)%next
    END DO
    heap%head = h1
  END SUBROUTINE local_heap_union


  FUNCTION heap_extract_min(heap, ithrd, opt_min_prev, opt_min_node) RESULT(node)
    INTEGER :: node
    TYPE(t_heap),  INTENT(INOUT)        :: heap
    INTEGER,       INTENT(IN)           :: ithrd
    INTEGER,       INTENT(IN), OPTIONAL :: opt_min_prev, opt_min_node
    ! local variables
    INTEGER :: prev

    IF (PRESENT(opt_min_prev) .AND. PRESENT(opt_min_node)) THEN
      prev = opt_min_prev
      node = opt_min_node
    ELSE
      CALL heap_min(heap, prev, node, ithrd)
    END IF
    IF (node == INVALID_NODE)  RETURN

    IF (prev /= INVALID_NODE) THEN
      node_storage(ithrd)%v(prev)%next = node_storage(ithrd)%v(node)%next
    ELSE
      heap%head = node_storage(ithrd)%v(node)%next
    END IF
    CALL local_heap_union(heap, heap_reverse(node_storage(ithrd)%v(node)%child, ithrd), ithrd)
  END FUNCTION heap_extract_min


  !>  insert (and reinitialize) a node into the heap
  SUBROUTINE heap_insert(heap, inode, ithrd)
    TYPE(t_heap), INTENT(INOUT) :: heap
    INTEGER,      INTENT(IN)    :: inode, ithrd
    ! local variables
    INTEGER :: min, hc
    TYPE(t_heap_node), POINTER :: node

    node => node_storage(ithrd)%v(inode)
    node%child  = INVALID_NODE
    node%parent = INVALID_NODE
    node%next   = INVALID_NODE
    node%degree = 0
    IF (heap%min /= INVALID_NODE) THEN
!CDIR IEXPAND
      hc = heap_cmp(node_storage(ithrd)%v(inode)%rdata, node_storage(ithrd)%v(heap%min)%rdata)
      IF (hc > 0) THEN
        ! swap min cache
        min = heap%min
!CDIR IEXPAND
        node_storage(ithrd)%v(min) = heap_node_init(t_heap_data(0._wp,0,0,0,0))
        node_storage(ithrd)%v(min)%degree = 0
        CALL local_heap_union(heap, min, ithrd)
        heap%min = inode
      ELSE
        CALL local_heap_union(heap, inode, ithrd)
      END IF
    ELSE
      CALL local_heap_union(heap, inode, ithrd)
    END IF
  END SUBROUTINE heap_insert


  SUBROUTINE uncache_min(heap, ithrd)
    TYPE(t_heap), INTENT(INOUT) :: heap
    INTEGER, INTENT(IN) :: ithrd
    ! local variables
    INTEGER :: min

    IF (heap%min /= INVALID_NODE) THEN
      min = heap%min
      heap%min = INVALID_NODE
      CALL heap_insert(heap, min, ithrd)
    END IF
  END SUBROUTINE uncache_min


  !>  merge addition into target
  SUBROUTINE heap_union(dest, addition, ithrd)
    TYPE(t_heap), INTENT(INOUT) :: dest, addition
    INTEGER, INTENT(IN) :: ithrd

    ! first insert any cached minima, if necessary
    CALL uncache_min(dest, ithrd)
    CALL uncache_min(addition, ithrd)
    CALL local_heap_union(dest, addition%head, ithrd)
    ! this is a destructive merge
    addition%head = INVALID_NODE
  END SUBROUTINE heap_union


  !> Remove minimum element from heap. 
  !
  !  The optional arguments @p opt_min_prev, @p opt_min_node may be provided
  !  if the minimum has been computed beforehand.
  FUNCTION heap_take(heap, ithrd, opt_min_prev, opt_min_node) RESULT(node)
    INTEGER :: node
    TYPE(t_heap), INTENT(INOUT)   :: heap
    INTEGER, INTENT(IN)           :: ithrd
    INTEGER, INTENT(IN), OPTIONAL :: opt_min_prev, opt_min_node

    IF (heap%min == INVALID_NODE) &
      &  heap%min = heap_extract_min(heap, ithrd, opt_min_prev, opt_min_node)
    node = heap%min
    heap%min = INVALID_NODE
    IF (node /= INVALID_NODE)  node_storage(ithrd)%v(node)%degree = NOT_IN_HEAP
  END FUNCTION heap_take


  !> Same as heap_take, but sums data of all items with the same key.
  FUNCTION heap_take_accumulated(heap, ithrd) RESULT(pdata)
    TYPE(t_heap_data) :: pdata
    TYPE(t_heap), INTENT(INOUT) :: heap
    INTEGER,      INTENT(IN)    :: ithrd
    ! local variables:
    INTEGER :: hn, prev, node

    hn = heap_take(heap, ithrd)
    pdata = node_storage(ithrd)%v(hn)%rdata
    DO
      CALL heap_min(heap, prev, node, ithrd)
      IF (node == INVALID_NODE) EXIT
      IF (heap_cmp(pdata, node_storage(ithrd)%v(node)%rdata) /= 0)   EXIT
      pdata%wgt = pdata%wgt + node_storage(ithrd)%v(node)%rdata%wgt
      hn = heap_take(heap, ithrd, prev, node)
    END DO
  END FUNCTION heap_take_accumulated


  SUBROUTINE heap_add_offset(nodes, nn, offset)
    TYPE(t_heap_node), INTENT(INOUT) :: nodes(:)
    INTEGER, INTENT(IN)              :: nn, offset
    ! local variables
    INTEGER :: i, parent, next, child

    DO i=1,nn
      parent = nodes(i)%parent
      next   = nodes(i)%next
      child  = nodes(i)%child
      IF (parent /= INVALID_NODE) &
        &  nodes(i)%parent = parent + offset
      IF (next   /= INVALID_NODE) &
        &  nodes(i)%next   = next   + offset
      IF (child  /= INVALID_NODE) &
        &  nodes(i)%child  = child  + offset
    END DO
  END SUBROUTINE heap_add_offset

END MODULE mo_util_binheap
