!>
!! Routine to invoke the index search for coupling of ICON grids
!!
!! <Description>
!! 
!! The routine takes lists of integers refering to global indices
!! of grid coordinates. The assumption here is that the indices in
!! the local lists are stored in ascending order. It is important to
!! note that this simple search does not perform any geographical
!! search. We assume that a grid point with global index n is on
!! component A is located on the same geographical position as a grid
!! point with global index n on component B. The complete list of
!! global indices in component B may contain only a subset of the
!! global indices which are present in component A and vice versa.
!!
!! Currently the upper and lower global index of each process is
!! distributed to all other participating processes. Common intersections
!! are identified. We assume a bidirectional exchange of data later on
!! therefore the initial search is done on all processes, as they
!! can act later as source or target. By doing it this way each process
!! can determine the total number of messages that have to be exchanged
!! without communicating further information.
!!
!! In a second step each target process sends the subset of points that
!! lie in the common range to the respective source process(es). Each
!! source process now determines all matching locations. This list is
!! stored on the source process (for later sending the data) and transferred
!! back to the target process (for receiving the data). This allows the
!! exchange of compact lists during the exchange of physical field later.  
!!
!! Function location
!!
!! The task to locate the postion global index (coordinate) in a list
!! of global indices of the target on the source side is done with
!! a simple bisection algorithm. (see Numerical Recipies in Fortran, 2nd
!! Edition, Sec. 3.4 How to Search an Ordered Table, page 110 ff.)
!!
!! @author Rene Redler, MPI-M
!!
!! $Id:$
!!
!! @par Revision History
!! first implementation by Rene Redler (2010-02-13)
!!
!! @par Copyright
!! 2010-2011 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and WARRANTY conditions.
!! 
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! &ltol>
!! &ltli> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! &ltli> The code may not be re-distributed without the consent of the authors.
!! &ltli> The copyright notice and statement of authorship must appear in all
!!    copies.
!! &ltli> You accept the warranty conditions (see WARRANTY).
!! &ltli> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!!
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_icon_cpl_search

#ifndef NOMPI

  USE mpi,         ONLY : MPI_SUCCESS, MPI_SOURCE, MPI_INTEGER,   &
       &                  MPI_STATUS_SIZE, MPI_TAG, MPI_ANY_SOURCE

  USE mo_icon_cpl, ONLY : target_locs, t_target_struct,           & 
       &                  source_locs,                            &
       &                  msg_len, initag,                        &
       &                  cplout, debug_coupler, debug_level,           &
       &                  nbr_active_comps, nbr_active_grids,     &
       &                  ICON_root, ICON_comm, ICON_comm_active, &
       &                  ICON_global_rank, ICON_global_size,     &
       &                  ICON_local_rank,                        &
       &                  grids, t_grid

  USE mo_icon_cpl_send_restart, ONLY : ICON_cpl_send_restart

  USE mo_master_control, ONLY: get_my_model_no

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: incr = 256

  INTEGER :: rstatus(MPI_STATUS_SIZE) ! MPI_Irecv status
  INTEGER :: wstatus(MPI_STATUS_SIZE) ! MPI_Wait status

  INTEGER :: i, j, ii, n
  INTEGER :: source_list_len

  INTEGER :: source_rank
  INTEGER :: msgtag

  INTEGER :: comp_id
  INTEGER :: grid_id

  INTEGER :: grid_extent (msg_len)
  INTEGER :: msg_to_src  (msg_len)
  INTEGER :: msg_to_tgt  (msg_len)

  INTEGER :: color
  INTEGER :: key
  INTEGER :: index

  INTEGER :: idx_range(2)

  INTEGER :: n_answers2recv

  INTEGER, ALLOCATABLE :: msg_fm_tgt (:,:)
  INTEGER, ALLOCATABLE :: msg_fm_src (:,:)
  INTEGER, ALLOCATABLE :: lrequests  (:)

  INTEGER, ALLOCATABLE :: all_extents(:,:)
  INTEGER, ALLOCATABLE :: idx(:)
  INTEGER, ALLOCATABLE :: grid_global_index(:)
  INTEGER, ALLOCATABLE :: grid_global_position(:)
  INTEGER, ALLOCATABLE :: global_index_rank(:)

  INTEGER, POINTER     :: srcbuffer     (:)
  INTEGER, POINTER     :: tgtbuffer     (:)
  INTEGER, POINTER     :: new_srcbuffer (:)
  INTEGER, POINTER     :: new_tgtbuffer (:)

  TYPE (t_grid), POINTER           :: gptr
  TYPE (t_target_struct), POINTER  :: tptr

  INTEGER              :: nbr_local_points ! Marker to distinguish internal from halo points

  ! Return code for error handling

  CHARACTER(len=132)               :: err_string
  INTEGER                          :: len
  INTEGER                          :: ierr
  INTEGER                          :: ierror

#else

  IMPLICIT NONE

#endif

  PUBLIC :: ICON_cpl_search

CONTAINS

  SUBROUTINE ICON_cpl_search
#ifndef NOMPI
    ! -------------------------------------------------------------------
    ! Assertion
    ! -------------------------------------------------------------------

    IF ( nbr_active_grids > 1 .OR. nbr_active_comps > 1 ) THEN
       PRINT *, 'More than 1 grid and component per process not yet supported'
       CALL MPI_Abort ( ICON_comm, 1, ierr )
    ENDIF
    !
    ! -------------------------------------------------------------------
    ! Create a communicator which includes only active processes
    ! -------------------------------------------------------------------
    !
    key   = 1
    color = 0

    IF ( nbr_active_comps > 0 .AND. nbr_active_grids > 0 ) color = 1

    CALL MPI_Comm_split ( ICON_comm, color, key, ICON_comm_active, ierr )
    IF ( ierr /= MPI_SUCCESS ) THEN
       CALL MPI_Error_string ( ierr, err_string, len, ierror )
       WRITE  ( * , '(a14,i3,a)' ) 'Error on rank ', ICON_global_rank, err_string
    ENDIF

    ! -------------------------------------------------------------------
    ! Important note: From here onwards for the time being we only
    !                 support one local grid and one local component!
    ! -------------------------------------------------------------------

    grid_id = 1
    comp_id = 1

    gptr => grids(grid_id)

    ! Allocate memory and sort global index list that we got from ICON
    ! First we get rid of the halo points, next we sort the remaining
    ! values. The original position is kept which is later used for the
    ! reordering of the results.

    len = gptr%grid_shape(2) - gptr%grid_shape(1) + 1

    nbr_local_points = 0

    IF ( ASSOCIATED(gptr%glob_index_rank) ) THEN
       DO i = 1, len
          IF ( gptr%glob_index_rank(i) == ICON_local_rank ) nbr_local_points = nbr_local_points + 1
       ENDDO
    ENDIF
   
    ! Subsampling of internal points (w/o halo)

    IF ( ASSOCIATED(gptr%glob_index_rank) .AND. nbr_local_points > 0 ) THEN

       ALLOCATE ( grid_global_index(nbr_local_points), &
            &     grid_global_position(nbr_local_points) )

       ii = 0

       DO i = 1, len
          IF ( gptr%glob_index_rank(i) == ICON_local_rank ) THEN
             ii = ii + 1
             grid_global_index(ii) = gptr%grid_glob_index(i)
             grid_global_position(ii) = i
          ENDIF
       ENDDO

       len = ii

    ELSE

       ALLOCATE( grid_global_index(len), &
            &    grid_global_position(len) )

       grid_global_index = gptr%grid_glob_index

       DO i = 1, len
          grid_global_position(i) = i
       ENDDO

    ENDIF

    ALLOCATE(idx(len))
    DO i = 1, len
       idx(i) = i
    ENDDO

    CALL quicksort_index (grid_global_index, len, idx)
    
    IF ( debug_coupler .AND. debug_level > 1 ) THEN
       DO i = 1, len
          WRITE ( cplout , '(a,4i8)' ) ' Global indices without halo after sorting ', &
               &                         i, grid_global_index(i), idx(i), grid_global_position(i)
       ENDDO
    ENDIF

    !
    ! -------------------------------------------------------------------
    ! Determine local index range
    ! -------------------------------------------------------------------
    !

    IF ( debug_coupler ) &
      WRITE ( cplout , * ) ' nbr_active is ', nbr_active_comps, nbr_active_grids

    DO comp_id = 1, nbr_active_comps
       DO grid_id = 1, nbr_active_grids
          grid_extent(1) = grid_global_index(1)
          grid_extent(2) = grid_global_index(len)
          grid_extent(3) = get_my_model_no() ! Rene: this should not be in here.
          IF ( debug_coupler ) &
               WRITE ( cplout , * ) ' extents ',  grid_extent(1),  grid_extent(2),  grid_extent(3)
       ENDDO
    ENDDO

    ! -------------------------------------------------------------------

    IF ( debug_coupler ) &
         WRITE ( cplout , '(i4,a1,a11,3i8)' ) ICON_global_rank, ':', ' extent is ', grid_extent

    ! -------------------------------------------------------------------
    ! Eventually it is wiser to store only target extents on the source
    !  processes rather that having source extents in all_extents as well.
    ! -------------------------------------------------------------------
    !
    ! Allocate memory to store extents on all processes
    ! -------------------------------------------------------------------
    !
    ALLOCATE ( all_extents(msg_len,ICON_global_size), STAT = ierr )
    IF ( ierr > 0 ) THEN
       PRINT *, ' Error allocating all_extents '
       CALL MPI_Abort ( ICON_comm, 1, ierr )
    ENDIF

    ! Collect all extents

    CALL MPI_Allgather ( grid_extent, msg_len, MPI_Integer, &
         all_extents, msg_len, MPI_Integer, ICON_comm_active, ierr )

    IF ( debug_coupler .AND. debug_level > 0 ) THEN
       IF ( ICON_global_rank == ICON_Root ) THEN
          DO i = 1, ICON_global_size
             WRITE ( cplout, '(i3,a1,a13,3i8)' ) i, ':', ' all_extents ', &
                  all_extents(1,i), &
                  all_extents(2,i), &
                  all_extents(3,i)
          ENDDO
       ENDIF
    ENDIF
    !
    ! -------------------------------------------------------------------
    ! Determine range of common extent
    ! -------------------------------------------------------------------
    !
    n_answers2recv = 0

    DO i = 1, ICON_global_size

       IF ( all_extents(3,i) == get_my_model_no() ) CYCLE

       idx_range(1) = MAX(all_extents(1,i),grid_extent(1))
       idx_range(2) = MIN(all_extents(2,i),grid_extent(2))

       IF ( idx_range(1) <= idx_range(2) ) n_answers2recv =  n_answers2recv + 1

    ENDDO

    IF ( debug_coupler ) &
         WRITE ( cplout , '(a,i3,a,i3)' ) &
         ' Global rank ', ICON_global_rank, ' n_answers2recv ', n_answers2recv

    IF ( n_answers2recv > 0 ) THEN

       ALLOCATE ( source_locs(n_answers2recv), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating source_list '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       ALLOCATE ( target_locs(n_answers2recv), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating target_locs type '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       ! ----------------------------------
       ! Initialise pointer
       ! ----------------------------------

       DO n = 1, n_answers2recv
          NULLIFY(source_locs(n)%source_list)
          NULLIFY(target_locs(n)%source_list)
          NULLIFY(target_locs(n)%target_list)
          target_locs(n)%source_rank     = -999
          target_locs(n)%source_list_len = 0
          target_locs(n)%offset          = 0
       ENDDO

       ALLOCATE ( lrequests(n_answers2recv), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating lrequests '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       ALLOCATE ( msg_fm_tgt(msg_len,n_answers2recv), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating msg_fm_tgt '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       ALLOCATE ( msg_fm_src(msg_len,n_answers2recv), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating msg_fm_src'
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       !
       ! ----------------------------------------------------------------
       ! Start receiving header messgage 1 from target
       ! ----------------------------------------------------------------
       !
       DO n = 1, n_answers2recv
          msgtag = initag + 1
          CALL MPI_Irecv ( msg_fm_tgt(1,n), msg_len, MPI_INTEGER, MPI_ANY_SOURCE, &
               msgtag, ICON_comm_active, lrequests(n), ierr )
       ENDDO

       n = 0

       source_rank = -999

       DO i = 1, ICON_global_size

          IF ( all_extents(3,i) == get_my_model_no() ) CYCLE

          idx_range(1) = MAX(all_extents(1,i),grid_extent(1))
          idx_range(2) = MIN(all_extents(2,i),grid_extent(2))

          IF ( idx_range(1) <= idx_range(2) ) THEN

             n           = n + 1
             source_rank = i - 1

             IF ( debug_coupler ) &
                  WRITE ( cplout , '(a,i3,a,i3)' ) &
                  ' Global rank ', ICON_global_rank, ' found match with global rank ', source_rank

             !
             ! silly search to determine start and end dimensions
             !

             DO ii = 1, len
                IF ( grid_global_index(ii) <  idx_range(1) ) CYCLE
                IF ( grid_global_index(ii) >= idx_range(1) ) THEN
                   idx_range(1) = ii
                   EXIT
                ENDIF
             ENDDO

             DO ii = idx_range(1), len
                IF ( grid_global_index(ii) > idx_range(2) ) EXIT
             ENDDO
             IF ( ii > len ) THEN
                idx_range(2) = len
             ELSE
                idx_range(2) = ii - 1
             ENDIF

             !
             ! ----------------------------------------------------------
             ! Buffered send of list length to the source process 
             ! with rank source_rank
             ! ----------------------------------------------------------
             !
             ! 1st) send header messgage 1 to source
             ! ----------------------------------------------------------
             !
             msg_to_src(1) = idx_range(2) - idx_range(1) + 1
             msg_to_src(2) = idx_range(1)
             msg_to_src(3) = 999

             msgtag = initag + 1
             CALL psmile_bsend ( msg_to_src, msg_len, &
                  MPI_INTEGER, source_rank, msgtag, ICON_comm_active, ierr ) 
             !
             ! ----------------------------------------------------------
             ! 2nd) send data message 1 to source
             ! ----------------------------------------------------------
             !
             IF ( msg_to_src(1) > 0 ) THEN
                msgtag = initag + 2
                CALL psmile_bsend ( grid_global_index(idx_range(1):idx_range(2)), &
                     msg_to_src(1), MPI_INTEGER, source_rank, msgtag, ICON_comm_active, ierr )
             ENDIF

          ENDIF ! ( idx_range(1) <= idx_range(2) )

       ENDDO ! i = 1, ICON_global_size

    ENDIF ! ( n_answers2recv > 0 )

    !
    ! -------------------------------------------------------------------
    ! Receive lists of common index ranges from target processes
    ! -------------------------------------------------------------------
    !

    DO n = 1, n_answers2recv
       !
       ! ----------------------------------------------------------------
       ! 1st) wait for header message 1 from target
       ! ----------------------------------------------------------------
       !
       CALL MPI_Waitany ( n_answers2recv, lrequests, index, wstatus, ierr )

       tptr => target_locs(index)

       tptr%target_list_len = msg_fm_tgt(1,index)
       tptr%offset          = msg_fm_tgt(2,index)
       tptr%target_rank     = wstatus(MPI_SOURCE)

       ALLOCATE ( tptr%target_list(tptr%target_list_len), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating target_list '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       !
       ! ----------------------------------------------------------------
       ! 2nd) receive data message 1 from target
       ! ----------------------------------------------------------------
       !
       source_list_len = 0

       IF ( tptr%target_list_len > 0 ) THEN

          msgtag = wstatus(MPI_TAG) + 1
          CALL MPI_Recv ( tptr%target_list, tptr%target_list_len, MPI_INTEGER, &
               tptr%target_rank, msgtag, ICON_comm_active, rstatus, ierr )

          IF ( debug_coupler .AND. debug_level > 1 ) THEN
             DO i = 1,  tptr%target_list_len
                WRITE ( cplout , '(a,3i8)' ) 'Received target list', i, tptr%target_list(i)
             ENDDO
          ENDIF

          ! =====================
          ! >>> Beginn ASSERTION
          ! =====================

          DO i = 1, len - 1
             IF ( i < idx_range(2) ) THEN
                IF ( grid_global_index(i) > grid_global_index(i+1) ) THEN
                   PRINT *, 'ERROR: Currently ascending order of src list indices is assumed.'
                   PRINT *, i, grid_global_index(i), grid_global_index(i+1)
                   CALL MPI_Abort ( ICON_comm, 1, ierr )
                ENDIF
             ENDIF
          ENDDO

          DO i = 1, tptr%target_list_len - 1
             IF ( tptr%target_list(i) > tptr%target_list(i+1) ) THEN
                PRINT *, 'ERROR: Currently ascending order of tgt list indices is assumed .'
                PRINT *, i, tptr%target_list(i), tptr%target_list(i+1)
                CALL MPI_Abort ( ICON_comm, 1, ierr )
             ENDIF
          ENDDO

          ! ==================
          ! <<< End ASSERTION
          ! ==================
          
          !
          ! ----------------------------------------------------------------
          ! Compare tptr%target_lists with local grid and detect matching indices
          ! ----------------------------------------------------------------
          !
          
          i = 1

          DO WHILE ( grid_global_index(i) < tptr%target_list(1) .AND. i < len )
             i = i + 1
          ENDDO

          ii           = 1
          idx_range(1) = i
          idx_range(2) = len ! gptr%grid_shape(2)

          IF ( debug_coupler ) &
               WRITE ( cplout , '(a,2i8)' ) 'Index range ', idx_range (1), idx_range (2)

          source_list_len = incr

          ALLOCATE ( srcbuffer(incr), tgtbuffer(incr), STAT = ierr )
          IF ( ierr > 0 ) THEN
             PRINT *, ' Error allocating buffers '
             CALL MPI_Abort ( ICON_comm, 1, ierr )
          ENDIF

          j = 0

          DO i = idx_range(1), idx_range(2)

             IF ( grid_global_index(i) > tptr%target_list(tptr%target_list_len) ) EXIT

             ii = location ( tptr%target_list_len, tptr%target_list, grid_global_index(i) )

             IF ( ii > 0 ) THEN

                ! ---------------------------------
                ! Check whether we need more memory
                ! ---------------------------------

                IF ( j + 1 > source_list_len ) THEN

                   ! Allocate a new chunk of memory
                   ! ------------------------------

                   ALLOCATE ( new_srcbuffer(j+incr), new_tgtbuffer(j+incr), STAT = ierr )
                   IF ( ierr > 0 ) THEN
                      PRINT *, ' Error allocating buffers '
                      CALL MPI_Abort ( ICON_comm, 1, ierr )
                   ENDIF

                   ! Transfer data to new buffer
                   ! ------------------------------

                   new_srcbuffer (1:j) = srcbuffer(1:j)
                   new_tgtbuffer (1:j) = tgtbuffer(1:j)

                   ! Deallocate old buffer
                   ! ------------------------------

                   DEALLOCATE ( srcbuffer, tgtbuffer, STAT = ierr )
                   IF ( ierr > 0 ) THEN
                      PRINT *, ' Error deallocating buffers '
                      CALL MPI_Abort ( ICON_comm, 1, ierr )
                   ENDIF

                   ! Assign original pointer
                   ! ------------------------------

                   srcbuffer => new_srcbuffer
                   tgtbuffer => new_tgtbuffer

                ENDIF

                j = j + 1
                srcbuffer(j) = i
                tgtbuffer(j) = ii

                IF ( debug_coupler .AND. debug_level > 1 ) &
                     WRITE ( cplout , '(2i8,a,3i8)' ) i, grid_global_index(i), &
                     ' <-> ', j, tgtbuffer(j),  tptr%target_list_len

             ENDIF

          ENDDO ! i = idx_range(1), idx_range(2)

          source_list_len = j

       ENDIF

       !
       ! ----------------------------------------------------------------
       ! Store matching global indices on the source side
       ! ----------------------------------------------------------------
       !

       ALLOCATE ( source_locs(index)%source_list(source_list_len), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating source_list '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       ALLOCATE ( source_locs(index)%target_list(source_list_len), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating source_list '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       source_locs(index)%target_rank      = tptr%target_rank
       source_locs(index)%source_list_len  = source_list_len

       !
       ! Convert back to unsorted order
       !

       DO i = 1, source_list_len
         source_locs(index)%source_list(i) = grid_global_position(idx(srcbuffer(i)))
         source_locs(index)%target_list(i) = tgtbuffer(i)
       ENDDO

       IF (ASSOCIATED(srcbuffer) .AND. ASSOCIATED(tgtbuffer) ) THEN
          DEALLOCATE ( srcbuffer, tgtbuffer, STAT = ierr )
          IF ( ierr > 0 ) THEN
             PRINT *, ' Error deallocating buffers '
             CALL MPI_Abort ( ICON_comm, 1, ierr )
          ENDIF
       ENDIF

       !
       ! Deallocate target_list as it is not needed anymore
       !
       DEALLOCATE ( tptr%target_list, STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error deallocating target_list '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       ! ----------------------------------------------------------------
       ! Send each list back to target
       ! ----------------------------------------------------------------
       !
       !  1) first send header with number of entries
       ! ----------------------------------------------------------------
       !
       msg_to_tgt(1) = source_list_len
       msg_to_tgt(2) = tptr%offset
       msg_to_tgt(3) = 999

       msgtag = initag + 2 * ICON_global_size + 1

       CALL psmile_bsend ( msg_to_tgt, msg_len, &
            MPI_INTEGER, tptr%target_rank, msgtag, ICON_comm_active, ierr ) 

       IF ( debug_coupler .AND. debug_level > 1 ) &
            WRITE ( cplout , '(a13,i8,a4,i8,a10,i8)' ) 'Sending back ', msg_len, &
            ' to ', tptr%target_rank, ' with tag ', msgtag

       !
       ! ----------------------------------------------------------------
       !  2) send indices for which a match was found
       ! ----------------------------------------------------------------
       !
       IF ( source_list_len > 0 ) THEN
          msgtag = initag  + 2 * ICON_global_size + 2
          CALL psmile_bsend ( source_locs(index)%target_list, source_list_len, &
               MPI_INTEGER, tptr%target_rank, msgtag, ICON_comm_active, ierr ) 
       ENDIF

    ENDDO ! n = 1, n_answers2recv
    !
    ! -------------------------------------------------------------------
    ! Post receive of header messgage from source
    ! -------------------------------------------------------------------
    !
    msgtag = initag + 2 * ICON_global_size + 1

    DO n = 1, n_answers2recv

       CALL MPI_Irecv ( msg_fm_src(1,n), msg_len, MPI_INTEGER, &
            MPI_ANY_SOURCE, msgtag, ICON_comm_active, lrequests(n), ierr )

       IF ( debug_coupler ) &
            WRITE ( cplout , '(a16,i8,a9,i12)' ) 'Posting receive ', msg_len, &
            ' request ', lrequests(n)

    ENDDO

    ii = 0

    DO n = 1, n_answers2recv
       !
       ! ----------------------------------------------------------------
       ! Wait for header messages from source
       ! ----------------------------------------------------------------
       !
       CALL MPI_Waitany ( n_answers2recv, lrequests, index, wstatus, ierr )
       !
       ! ----------------------------------------------------------------
       ! Receive and store lists except empty lists (msg_fm_src(1,index)=0)
       ! ----------------------------------------------------------------
       !
       IF ( msg_fm_src(1,index) > 0 ) THEN

          ii = ii + 1
          tptr => target_locs(ii)

          tptr%source_list_len = msg_fm_src(1,index)
          tptr%offset          = msg_fm_src(2,index)
          tptr%source_rank     = wstatus(MPI_SOURCE)

          IF ( debug_coupler ) &
               WRITE ( cplout , '(a9,i8,a9,i12,a6,i8)' ) 'Received ', msg_len, ' request ', &
               lrequests(index), ' from ', wstatus(MPI_SOURCE)

          ALLOCATE ( tptr%source_list(tptr%source_list_len), &
                     tgtbuffer(tptr%source_list_len), STAT = ierr )
          IF ( ierr > 0 ) THEN
             PRINT *, ' Error allocating source_list '
             CALL MPI_Abort ( ICON_comm, 1, ierr )
          ENDIF

          msgtag = wstatus(MPI_TAG) + 1
          CALL MPI_Recv ( tgtbuffer, tptr%source_list_len, MPI_INTEGER, &
               tptr%source_rank, msgtag, ICON_comm_active, rstatus, ierr )

          ! Add offset on the target side and convert back to original ordering

          DO i = 1, tptr%source_list_len
             tptr%source_list(i) = idx(tgtbuffer(i)+tptr%offset-1)
          ENDDO

          IF ( debug_coupler ) &
               WRITE ( cplout , '(a,i8)' ) 'Offset added: ', tptr%offset - 1 

           DEALLOCATE (tgtbuffer)

       ENDIF

    ENDDO ! n = 1, n_answers2recv

    !
    ! -------------------------------------------------------------------
    ! Check whether we need to provide a restart field
    ! -------------------------------------------------------------------
    !
    CALL ICON_cpl_send_restart
    !
    ! -------------------------------------------------------------------
    ! Deallocate memory
    ! -------------------------------------------------------------------
    !
    IF (  n_answers2recv > 0 ) THEN

       DEALLOCATE ( all_extents, STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error deallocating all_extents '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       DEALLOCATE ( lrequests, STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error deallocating lrequests '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       DEALLOCATE ( msg_fm_tgt, STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error deallocating msg_fm_tgt '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

       DEALLOCATE ( msg_fm_src, STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error deallocating msg_fm_src '
          CALL MPI_Abort ( ICON_comm, 1, ierr )
       ENDIF

    ENDIF
#endif
  END SUBROUTINE ICON_cpl_search

  INTEGER FUNCTION location ( len, list, list_req )

    !
    ! simple bisection algorithm
    ! --------------------------

    INTEGER, INTENT(in)  :: len
    INTEGER, INTENT(in)  :: list(len)
    INTEGER, INTENT(in)  :: list_req

    INTEGER :: j1, j2, j3

    ! Initialize lower and upper limits

    location = -1

    j1 = 0

    j3 = len + 1

    DO WHILE ( j3-j1 > 1 )

       ! Compute midpoint

       j2 = (j3+j1)/2

       ! replace lower or upper limit

       IF ( (list(len) > list(1)) .EQV. (list_req > list(j2)) ) THEN
          j1 = j2
       ELSE
          j3 = j2
       ENDIF

       IF ( list_req == list(j2) ) THEN
          location = j2
          EXIT
       ENDIF

    ENDDO

  END FUNCTION location


  SUBROUTINE quicksort_index(a, n, idx)

    IMPLICIT NONE
    !
    ! !INPUT PARAMETERS:
    !
    INTEGER, INTENT(In)    :: n
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !
    INTEGER, INTENT(InOut) :: a(n)
    INTEGER, INTENT(InOut) :: idx(n)
    !
    ! !LOCAL VARIABLES
    !
    INTEGER, PARAMETER     :: pointer_inc = 64
    INTEGER                :: pointer_size
    INTEGER, POINTER       :: stackl(:), stackr(:)
    INTEGER, POINTER       :: new_stackl(:), new_stackr(:)

    INTEGER                :: i, j, k, l, r, s
    INTEGER                :: ww
    INTEGER                :: w, x
    !
    ! !DESCRIPTION:
    !
    ! Non-recursive stack version of Quicksort from N. Wirth's Pascal Book,
    ! 'Algorithms + Data Structures = Programms'.
    !
    ! taken from:
    !    http://www.nag.com/nagware/examples.asp
    !    http://www.nag.com/nagware/Examples/nur.f90
    !
    ! see also:
    !    http://en.wikipedia.org/wiki/Quicksort
    !
    ! !REVISION HISTORY:
    !
    !   Date      Programmer   Description
    ! ----------  ----------   -----------
    ! 19.07.95    Alan Miller  created
    ! 13.02.08    R. Redler    revised
    !
    !----------------------------------------------------------------------
    !
    !  Initialization
    !
    pointer_size = pointer_inc

    ALLOCATE(new_stackl(pointer_inc), new_stackr(pointer_inc))

    stackl => new_stackl
    stackr => new_stackr

    s = 1
    stackl(1) = 1
    stackr(1) = n
    !
    !  Start sorting
    !        
    ! ... keep taking the top request from the stack until s = 0.

    DO WHILE (s /= 0) ! 10  CONTINUE

       l = stackl(s)
       r = stackr(s)
       s = s - 1

       ! ... keep splitting a(l), ... ,a(r) until l>= r.

       DO WHILE (l <  r) ! 20  CONTINUE

          i = l
          j = r
          k = (l+r) / 2
          x = a(k)

          ! ... repeat until i > j.

          DO
             DO
                IF (a(i) < x) THEN ! Search from lower end
                   i = i + 1
                   CYCLE
                ELSE
                   EXIT
                END IF
             END DO

             DO
                IF (x < a(j)) THEN ! Search from upper end
                   j = j - 1
                   CYCLE
                ELSE
                   EXIT
                END IF
             END DO

             IF (i <= j) THEN ! Swap positions i & j

                w = a(i)
                ww = idx(i)
                a(i) = a(j)
                idx(i) = idx(j)
                a(j) = w
                idx(j) = ww
                i = i + 1
                j = j - 1
                IF (i.GT.j) EXIT

             ELSE

                EXIT

             END IF

          END DO

          IF (j-l >= r-i) THEN

             IF (l < j) THEN

                IF ( s+1 > pointer_size ) THEN
                   pointer_size = pointer_size + pointer_inc
                   ALLOCATE(new_stackl(pointer_size), new_stackr(pointer_size))
                   new_stackl(1:pointer_size-pointer_inc) = &
                        stackl(1:pointer_size-pointer_inc)
                   new_stackr(1:pointer_size-pointer_inc) = &
                        stackr(1:pointer_size-pointer_inc)
                   DEALLOCATE(stackl, stackr)
                   stackl => new_stackl
                   stackr => new_stackr
                ENDIF

                s = s + 1
                stackl(s) = l
                stackr(s) = j

             END IF

             l = i

          ELSE

             IF (i < r) THEN

                IF ( s+1 > pointer_size ) THEN
                   pointer_size = pointer_size + pointer_inc
                   ALLOCATE(new_stackl(pointer_size), new_stackr(pointer_size))
                   new_stackl(1:pointer_size-pointer_inc) = &
                        stackl(1:pointer_size-pointer_inc)
                   new_stackr(1:pointer_size-pointer_inc) = &
                        stackr(1:pointer_size-pointer_inc)
                   DEALLOCATE(stackl, stackr)
                   stackl => new_stackl
                   stackr => new_stackr
                ENDIF

                s = s + 1
                stackl(s) = i
                stackr(s) = r

             END IF

             r = j

          END IF

       ENDDO ! WHILE (l <  r) GO TO 20

    ENDDO ! WHILE (s /= 0) GO TO 10

  END SUBROUTINE quicksort_index

END MODULE mo_icon_cpl_search
