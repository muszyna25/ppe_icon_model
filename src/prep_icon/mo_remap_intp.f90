!! This module contains subroutines for allocating interpolation coefficient tables,
!! for reduction of multi-threaded coefficient tables and for the interpolation
!! operation itself.
!!
!! @author F. Prill, DWD
!!
MODULE mo_remap_intp

#if !defined(HAVE_NOMPI)
#ifdef __SX__
  USE MPI
#endif
#endif

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_communication,      ONLY: blk_no, idx_no
  USE mo_io_units,           ONLY: nnml
  USE mo_namelist,           ONLY: POSITIONED, position_nml, open_nml, close_nml
  USE mo_util_sort,          ONLY: quicksort
  USE mo_mpi,                ONLY: p_n_work, p_comm_work, p_int, p_real_dp,    &
    &                              p_commit_type_struct, p_alltoall
  USE mo_util_binheap,       ONLY: t_heap_data, t_heap, heap_cmp, heap_init,   &
    &                              heap_empty, heap_take_accumulated,          &
    &                              heap_union, node_storage, nnode,            &
    &                              INVALID_NODE, heap_add_offset,              &
    &                              resize_node_storage, node_storage_finalize, &
    &                              get_free_node, heap_node_init, heap_insert
  USE mo_remap_config,       ONLY: dbg_level, MAX_NSTENCIL
  USE mo_remap_shared,       ONLY: t_grid

  IMPLICIT NONE

  PRIVATE
  ! subroutines
  PUBLIC :: allocate_intp_data
  PUBLIC :: deallocate_intp_data
  PUBLIC :: interpolate_c
  PUBLIC :: merge_heaps
  PUBLIC :: sync_foreign_wgts
  PUBLIC :: reduce_mthreaded_weights
  PUBLIC :: read_interp_namelist
  ! variables and data types
  PUBLIC :: s_maxsize
  PUBLIC :: t_intp_data
  PUBLIC :: t_intp_data_mt

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_intp')

#if !defined(HAVE_NOMPI)
#ifndef __SX__
  INCLUDE "mpif.h"
#endif
#endif
  
  ! Maximum size of sequential list for very large stencils
  INTEGER  :: s_maxsize

  ! Threshold. Weights smaller than this value are neglected.
  REAL(wp), PARAMETER :: W_THRESHOLD      = 1e-10_wp

  ! Distributed computation: max. size of weight storage for other PEs
  ! @todo Compute this value dependent on grid size.
  INTEGER,  PARAMETER :: MAX_NFOREIGN     = 50000

  ! namelist definition: main namelist
  NAMELIST/interp_nml/   s_maxsize

  !> data structure containing interpolation coefficients
  !
  ! @note For vectorization we assume a fixed stencil size. If this
  !       stencil size is chosen too small, we must disregard some of
  !       the interpolation weights (the smallest ones)
  !
  TYPE t_intp_data
    !--- array-structured stencil (processed fast)
    !> fixed stencil size
    INTEGER :: nstencil
    !> interpolation weights
    REAL(wp), ALLOCATABLE :: wgt(:,:,:)                ! (nstencil, nproma, nblks_out)
    !> cell area (for normalization)
    REAL(wp), ALLOCATABLE :: area(:,:)                 ! (nproma, nblks_out)
    !> destination stencil indices (index/block)
    INTEGER, ALLOCATABLE  :: iidx(:,:,:), iblk(:,:,:)  ! (nstencil, nproma, nblks_out)
    !> local stencil size (<= nstencil)
    INTEGER, ALLOCATABLE  :: nidx(:,:)                 ! (nproma, nblks_out)

    !--- sequential list (processed less efficiently)
    !> length of sequential list (for cells with very large stencils):
    INTEGER :: s_nlist
    !> sequential list
    TYPE (t_heap_data), ALLOCATABLE :: sl(:)           ! (1,...,s_nlist)

    !> minimum and maximum index for weights associated with a destination
    !  cell in sequential list (if there are any entries)
    INTEGER, ALLOCATABLE :: smin(:,:), smax(:,:)

    !> maximum stencil of a destination point (standard stencil + seq. list)
    INTEGER :: max_totstencil
   END TYPE t_intp_data

  !> multi-threaded storage of interpolation coefficients
  !
  TYPE t_intp_data_mt
    TYPE (t_heap), ALLOCATABLE :: wgt_heap(:)          ! (thread)
    !> distributed computation: list of weights belonging to other PEs
    TYPE (t_heap_data) :: foreign_wgt(MAX_NFOREIGN)
    INTEGER            :: foreign_pe(MAX_NFOREIGN)
    INTEGER            :: nforeign
  END TYPE t_intp_data_mt

  ! generic interface
  INTERFACE allocate_intp_data
    MODULE PROCEDURE allocate_intp_data_sthreaded
    MODULE PROCEDURE allocate_intp_data_mthreaded
  END INTERFACE

  ! generic interface
  INTERFACE deallocate_intp_data
    MODULE PROCEDURE deallocate_intp_data_sthreaded
    MODULE PROCEDURE deallocate_intp_data_mthreaded
  END INTERFACE

  ! generic interface
  INTERFACE interpolate_c
    MODULE PROCEDURE interpolate_c_2D
  END INTERFACE
  
CONTAINS

  !> Opens namelist file, reads interpolation config.
  !
  SUBROUTINE read_interp_namelist(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename !< main namelist file
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::read_interp_namelist')
    INTEGER :: istat

    IF (dbg_level >= 5) WRITE (0,*) "# read interpolation namelist"
    ! default settings
    s_maxsize        = 500000
    ! read user's (new) specifications
    CALL open_nml(TRIM(filename))
    CALL position_nml ('interp_nml', status=istat)
    IF (istat == POSITIONED) THEN
      READ (nnml, interp_nml)
    ELSE
      CALL finish(routine, "Internal error!")
    END IF
    CALL close_nml
    ! status output
    IF (dbg_level >= 2)  WRITE (0,*) "# stencil size: ", MAX_NSTENCIL, "/", s_maxsize
  END SUBROUTINE read_interp_namelist


  !> Allocate data structure for interpolation coefficients.
  !
  SUBROUTINE allocate_intp_data_sthreaded(intp_data, grid)
    TYPE(t_intp_data),    INTENT(INOUT) :: intp_data   !< data structure with intp. weights
    TYPE (t_grid),        INTENT(IN)    :: grid        !< output grid 
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::allocate_intp_data')
    INTEGER :: ierrstat, nstencil

    nstencil = MAX_NSTENCIL
    intp_data%nstencil       = nstencil
    intp_data%max_totstencil = 0
    ALLOCATE(intp_data%wgt (nstencil, nproma, grid%p_patch%nblks_c), &
      &      intp_data%area(nproma, grid%p_patch%nblks_c),           &
      &      intp_data%iidx(nstencil, nproma, grid%p_patch%nblks_c), &
      &      intp_data%iblk(nstencil, nproma, grid%p_patch%nblks_c), &
      &      intp_data%nidx(nproma, grid%p_patch%nblks_c),           &
      &      intp_data%smin(nproma, grid%p_patch%nblks_c),           &
      &      intp_data%smax(nproma, grid%p_patch%nblks_c),           &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    intp_data%wgt(:,:,:)  =  0._wp
    intp_data%iidx(:,:,:) =  1
    intp_data%iblk(:,:,:) =  1
    intp_data%nidx(:,:)   =  0
    intp_data%area(:,:)   =  -1._wp
    intp_data%smin(:,:)   =  1
    intp_data%smax(:,:)   =  0

    ALLOCATE(intp_data%sl(S_MAXSIZE), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    intp_data%s_nlist = 0

  END SUBROUTINE allocate_intp_data_sthreaded


  !> Allocate MULTI-THREADED data structure for interpolation coefficients.
  !
  SUBROUTINE allocate_intp_data_mthreaded(intp_data, nthreads)
    TYPE(t_intp_data_mt), INTENT(INOUT) :: intp_data   !< data structure with intp. weights
    INTEGER,              INTENT(IN)    :: nthreads    !< no. of threads
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = &
      &  TRIM(TRIM(modname)//'::allocate_intp_data_mthreaded')
    INTEGER :: ierrstat, ithrd

    ALLOCATE(intp_data%wgt_heap(nthreads), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    DO ithrd=1,nthreads
      CALL heap_init(intp_data%wgt_heap(ithrd))
    END DO
    intp_data%nforeign = 0
  END SUBROUTINE allocate_intp_data_mthreaded


  !> Clear data structure for interpolation coefficients.
  !
  SUBROUTINE deallocate_intp_data_sthreaded(intp_data)
    TYPE(t_intp_data), INTENT(INOUT) :: intp_data
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//'::deallocate_intp_data')
    INTEGER :: ierrstat

    intp_data%nstencil = 0
    DEALLOCATE(intp_data%wgt, intp_data%area, intp_data%iidx, intp_data%iblk, &
      &        intp_data%nidx, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    DEALLOCATE(intp_data%sl, intp_data%smin, intp_data%smax, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE deallocate_intp_data_sthreaded


  !> Clear data structure for interpolation coefficients.
  !
  SUBROUTINE deallocate_intp_data_mthreaded(intp_data)
    TYPE(t_intp_data_mt), INTENT(INOUT) :: intp_data
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = &
      &      TRIM(TRIM(modname)//'::deallocate_intp_data_mthreaded')
    INTEGER :: ierrstat

    intp_data%nforeign = 0
    DEALLOCATE(intp_data%wgt_heap, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
  END SUBROUTINE deallocate_intp_data_mthreaded


  !> Merges two binomial heaps into a single heap
  !
  SUBROUTINE merge_heaps2(intp_data_mt1, intp_data_mt2, i1,i2)
    TYPE(t_intp_data_mt), INTENT(INOUT) :: intp_data_mt1, intp_data_mt2 ! modified!
    INTEGER,              INTENT(IN)    :: i1, i2
    ! local variables
    INTEGER :: nnode1, nnode2
    TYPE (t_heap) :: heap

    nnode1 = nnode(i1)
    nnode2 = nnode(i2)
    IF (nnode2 > 0) THEN
      ! copy node from i2 to i1 (correcting the offset)
      IF (dbg_level >= 10) &
        &   WRITE (0,*) "# copy nodes from ", i2, " to ", i1
      CALL heap_add_offset(node_storage(i2)%v, nnode2, nnode1)
      IF ((nnode1+SIZE(node_storage(i2)%v)) > SIZE(node_storage(i1)%v)) &
        &  CALL resize_node_storage(nnode2, i1)
      node_storage(i1)%v(nnode1+1:nnode1+nnode2) = node_storage(i2)%v(1:nnode2)
      nnode(i1)     = nnode1 + nnode2
      CALL node_storage_finalize(i2)
    END IF

    IF (dbg_level >= 10) &
      &   WRITE (0,*) "# merging heaps ", i2, " to ", i1
    heap = intp_data_mt1%wgt_heap(i2)
    IF (heap%head /= INVALID_NODE) heap%head = heap%head + nnode1
    IF (heap%min  /= INVALID_NODE) heap%min  = heap%min  + nnode1
    ! merge the two heaps
    CALL heap_union(intp_data_mt1%wgt_heap(i1), heap, i1)

    heap = intp_data_mt2%wgt_heap(i2)
    IF (heap%head /= INVALID_NODE) heap%head = heap%head + nnode1
    IF (heap%min  /= INVALID_NODE) heap%min  = heap%min  + nnode1
    ! merge the two heaps
    CALL heap_union(intp_data_mt2%wgt_heap(i1), heap, i1)
  END SUBROUTINE merge_heaps2


  !> Merges *two or more* binomial heaps into a single heap
  !
  SUBROUTINE merge_heaps(intp_data_mt1, intp_data_mt2, nthreads)
    TYPE(t_intp_data_mt), INTENT(INOUT) :: intp_data_mt1, intp_data_mt2 ! modified!
    INTEGER,              INTENT(IN)    :: nthreads
    ! local variables:
    INTEGER                             :: nlists, new_nlists, i,                     &
      &                                    listidx(nthreads), new_listidx(nthreads)

    !-- reduce multi-threaded heaps to a single one
    IF (dbg_level >= 2)  WRITE (0,*) "# reducing multi-threaded weights."
    nlists  = nthreads
    listidx = (/ (i, i=1,nthreads) /)
    DO 
      IF (dbg_level >= 10) &
        &   WRITE (0,*) "# remaining sequential lists: ", listidx(1:nlists)
      ! merge pairs of sequential lists
      IF (nlists > 1) THEN
!$OMP PARALLEL DO 
        DO i=1,(nlists-1),2
          CALL merge_heaps2(intp_data_mt1, intp_data_mt2, listidx(i),listidx(i+1))
        END DO
!$OMP END PARALLEL DO
      END IF
      new_nlists  = nlists / 2
      new_listidx = (/ (listidx(i), i=1,nlists,2) /)
      IF (MOD(nlists,2) == 1) THEN
        new_nlists = new_nlists + 1
        new_listidx(new_nlists) = listidx(nlists)
      END IF
      listidx = new_listidx
      nlists  = new_nlists
      IF (nlists == 1) EXIT
    END DO
  END SUBROUTINE merge_heaps


  !> Distributed computation: Communicate lists of weights on PE
  !  boundaries, add them to multi-threaded binomial heap.
  !
  SUBROUTINE sync_foreign_wgts(grid, intp_data_mt)
    TYPE (t_grid),        INTENT(IN)    :: grid !< local grid (covering)
    TYPE(t_intp_data_mt), INTENT(INOUT) :: intp_data_mt
    ! local variables
#if !defined(HAVE_NOMPI)
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(TRIM(modname)//':sync_foreign_wgts')
    INTEGER,          PARAMETER :: thrd0 = 1 ! (we have already reduced data to first thread)
    INTEGER :: ierrstat, nn, i, pe, ierr, nrecv_tot, hn, idx_cov, &
      &        idx_glb, jc_glb, jb_glb, idx_local, &
      &        oldtypes(2),blockcounts(2), mpi_heap_data_type
    INTEGER, ALLOCATABLE :: perm(:), glb2cov(:,:)
    INTEGER :: offset(p_n_work), local_size(p_n_work), recvcounts(p_n_work), &
      &        rdispls(p_n_work)
    TYPE (t_heap_data), ALLOCATABLE :: sendbuf(:), recvbuf(:)
    TYPE (t_heap_data) :: wgt_data

    ! allocate temporary data structures
    nn = intp_data_mt%nforeign
    ALLOCATE(perm(nn), sendbuf(nn), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! commit data type "t_heap_data" to MPI:
    ! cf. setup description of t_heap_data
    oldtypes(:)    = (/ p_real_dp,    p_int /)
    blockcounts(:) = (/         1,        4 /)
    ! define structured type and commit it  
    mpi_heap_data_type = p_commit_type_struct(oldtypes, blockcounts)

    ! sort list of weights for foreign PEs 
    IF (nn > 0) THEN
      perm(:) = (/ ( i, i=1,nn) /)
      CALL quicksort(intp_data_mt%foreign_pe(1:nn), perm, 1,nn)
      DO i=1,nn
        sendbuf(i) = intp_data_mt%foreign_wgt(perm(i))
      END DO
      ! determine the number of entries to send to each PE
      offset(:)     = 0
      local_size(:) = 0
      DO i=1,nn
        pe = intp_data_mt%foreign_pe(i) + 1
        local_size(pe) = local_size(pe) + 1
        IF (offset(pe) == 0) offset(pe)=i
      END DO
      offset(:) = offset(:) - 1
    ELSE
      offset(:)     = 0
      local_size(:) = 0
    END IF

    ! communicate the number of entries to be sent:
    recvcounts(:) = 0
    CALL p_alltoall(local_size, recvcounts, p_comm_work)

    ! compute offsets in receive buffer:
    rdispls(:) = recvcounts(:)
    DO i=2,p_n_work
      rdispls(i) = rdispls(i) + rdispls(i-1)
    END DO
    DO i=p_n_work,2,-1
      rdispls(i) = rdispls(i-1)
    END DO
    rdispls(1) = 0

    ! allocate temporary data structures for sending/receiving
    nrecv_tot = SUM(recvcounts)
    ALLOCATE(recvbuf(nrecv_tot), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! communicate weights
    ! note: we cannot move this call to mo_mpi module since we would need
    !       the definition of "t_heap_data" there.
    CALL MPI_ALLTOALLV(sendbuf,local_size, offset, mpi_heap_data_type,  & ! sendbuf,sendcounts,sdispls,sendtype
      &                recvbuf,recvcounts,rdispls,mpi_heap_data_type,   & ! recvbuf,recvcounts,rdispls,recvtype
      &                p_comm_work, ierr)

    ! communicate mapping global->covering
    ALLOCATE(glb2cov(nproma, blk_no(grid%p_patch%n_patch_cells_g)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    
    ! mapping of source index: global->covering
    glb2cov(:,:) = -1
    DO idx_local=1,grid%p_patch%n_patch_cells
      idx_glb   = grid%p_patch%cells%glb_index(idx_local)
      jc_glb    = idx_no(idx_glb)
      jb_glb    = blk_no(idx_glb)    
      glb2cov(jc_glb,jb_glb) = idx_local ! grid%cov_c(idx_local)
    END DO

    ! insert the received weights in to the binomial heap
    DO i=1,nrecv_tot
      wgt_data = recvbuf(i)
      idx_cov = glb2cov(wgt_data%sidx, wgt_data%sblk)
      wgt_data%sidx = idx_no(idx_cov)
      wgt_data%sblk = blk_no(idx_cov)

      hn = get_free_node(thrd0)
      node_storage(thrd0)%v(hn) = heap_node_init(wgt_data)
      CALL heap_insert(intp_data_mt%wgt_heap(thrd0), hn, thrd0)
    END DO

    ! clean up
    DEALLOCATE(perm, sendbuf, recvbuf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
    intp_data_mt%nforeign = 0
#endif
  END SUBROUTINE sync_foreign_wgts


  !> Reduce a multi-threaded interpolation weight data structure to a
  !> single-threaded one.
  !
  SUBROUTINE reduce_mthreaded_weights(intp_data_mt, intp_data)
    TYPE(t_intp_data_mt), INTENT(INOUT) :: intp_data_mt ! will be modified!
    TYPE(t_intp_data),    INTENT(INOUT) :: intp_data
    ! local variables:
    TYPE (t_heap_data)                  :: wgt_data
    INTEGER                             :: seq_idx, maxstencil, nstencil, didx, dblk, &
      &                                    icount, thrd0

    IF (dbg_level >= 2)  &
      &  WRITE (0,*) "# copy weights to final data structure."
    thrd0 = 1 ! (we have already reduced data to first thread)
    IF (dbg_level >= 5) &
      &   WRITE (0,*) "# fill interpolation stencil with data from heap ", thrd0
    maxstencil = intp_data%nstencil
    icount = 0
    !-- read heap minima successively into standard stencil and sequential list
    DO
      IF (heap_empty(intp_data_mt%wgt_heap(thrd0))) EXIT
      wgt_data = heap_take_accumulated(intp_data_mt%wgt_heap(thrd0), thrd0)
      IF (wgt_data%wgt < W_THRESHOLD) CYCLE

      didx = wgt_data%didx
      dblk = wgt_data%dblk
      nstencil = intp_data%nidx(didx,dblk)
      icount = icount + 1
      IF (nstencil >= maxstencil) THEN
        seq_idx = intp_data%s_nlist + 1
        ! update bounds for minimum and maximum index associated with a destination
        ! cell in sequential list:
        IF (intp_data%smax(didx,dblk) == 0) THEN
          intp_data%smin(didx,dblk) = seq_idx
        END IF
        intp_data%smax(didx,dblk) = seq_idx
        ! update maximum stencil of all destination points
        ! (standard stencil + seq. list)
        intp_data%max_totstencil = MAX(intp_data%max_totstencil, &
          & (maxstencil + intp_data%smax(didx,dblk) - intp_data%smin(didx,dblk) + 1))

        IF (seq_idx > S_MAXSIZE) THEN
          WRITE (0,*) "didx, dblk = ", didx, dblk
          CALL finish("reduce_mthreaded_weights", "List size exceeded!")
        END IF
        intp_data%s_nlist     = seq_idx
        intp_data%sl(seq_idx) = wgt_data
      ELSE
        nstencil = nstencil + 1
        IF (nstencil == 1) intp_data%area(didx,dblk) = 0._wp
        intp_data%nidx(didx,dblk) = nstencil
        intp_data%wgt (nstencil,didx,dblk) = wgt_data%wgt
        intp_data%iidx(nstencil,didx,dblk) = wgt_data%sidx
        intp_data%iblk(nstencil,didx,dblk) = wgt_data%sblk
        intp_data%max_totstencil = MAX(intp_data%max_totstencil, nstencil)
      END IF
      intp_data%area(didx,dblk) = intp_data%area(didx,dblk) + wgt_data%wgt
    END DO
    IF (dbg_level >= 10) WRITE (0,*) "# extracted ", icount, " element(s)."
  END SUBROUTINE reduce_mthreaded_weights


  ! ---------------------------------------------------------------------------


  !> Perform interpolation operation: 2D case.
  !
  SUBROUTINE interpolate_c_2D(field1, field2, out_grid, intp_data)
    REAL(wp),           INTENT(IN)    :: field1(:,:)
    REAL(wp),           INTENT(INOUT) :: field2(:,:)
    TYPE (t_grid),      INTENT(IN)    :: out_grid
    TYPE (t_intp_data), INTENT(IN)    :: intp_data
    ! local variables
    INTEGER  :: i, jc,jb, nblks, npromz, &
      &         i_startidx, i_endidx
    REAL(wp) :: val
    TYPE(t_heap_data) :: t

    IF (dbg_level >= 11)  WRITE (0,*) "# perform horizontal interpolation."

    nblks    = out_grid%p_patch%nblks_c
    npromz   = out_grid%p_patch%npromz_c

!$OMP PARALLEL PRIVATE(i_startidx, i_endidx,jc,val)
!$OMP DO
    DO jb=1,nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks) i_endidx = npromz

      DO jc=i_startidx, i_endidx
        val = 0._wp
!CDIR UNROLL=MAX_NSTENCIL
        DO i=1,MAX_NSTENCIL
          val = val + intp_data%wgt(i,jc,jb) * &
            &         field1(intp_data%iidx(i,jc,jb),intp_data%iblk(i,jc,jb))
        END DO
        field2(jc,jb) = val
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! consider sequential weight list:
!$OMP PARALLEL PRIVATE(t,jc,jb,val)
!$OMP DO
    DO i=1,intp_data%s_nlist
      t   = intp_data%sl(i)
      jc  = t%didx
      jb  = t%dblk
      val = field2(jc,jb) + t%wgt * field1(t%sidx,t%sblk)
!$OMP CRITICAL
      field2(jc,jb) = val
!$OMP END CRITICAL
    END DO
!$OMP END DO
!$OMP END PARALLEL

    IF (dbg_level >= 11) THEN
      WRITE (0,*) "# minmax: ", MINVAL(field2(:,:)), MAXVAL(field2(:,:))
      WRITE (0,*) "# minmax loc: ", MINLOC(field2(:,:)), MAXLOC(field2(:,:))
    END IF
  END SUBROUTINE interpolate_c_2D

END MODULE mo_remap_intp
