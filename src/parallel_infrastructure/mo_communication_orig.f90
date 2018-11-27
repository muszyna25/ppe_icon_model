!>
!!               This module provides the communication routines.
!!
!!               This module provides the communication routines
!! for parallel runs
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
! (GZ, 2013-08-30): So far, the Cray compiler is the only one for which an OpenMP parallelization
! of copying data into / back from the MPI-buffer seems to give a benefit. Further compilers may
! be added here once the OpenMP implementation is sufficiently efficient
#if ((defined(_CRAYFTN) && !defined(_OPENACC)) || defined(__INTEL_COMPILER))
#define __OMPPAR_COPY__
#endif

!----------------------------
#include "icon_definitions.inc"
#include "crayftn_ptr_fail.inc"
!----------------------------
MODULE mo_communication_orig
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                 ONLY: dp, sp
USE mo_exception,            ONLY: finish, message, message_text
USE mo_mpi,                  ONLY: p_send, p_recv, p_irecv, p_wait, p_isend,                &
     &                             p_comm_work, my_process_is_mpi_seq,                      &
     &                             get_my_mpi_work_communicator, get_my_mpi_work_comm_size, &
     &                             get_my_mpi_work_id, p_gather, p_gatherv,                 &
     &                             p_alltoall, p_comm_rank, p_comm_size,                    &
     &                             p_barrier,                                               &
     &                             p_comm_is_intercomm, p_comm_remote_size
USE mo_parallel_config,      ONLY: iorder_sendrecv, nproma, itype_exch_barrier
USE mo_timer,                ONLY: timer_start, timer_stop, timer_exch_data, &
     &                             timer_barrier, &
     &                             timer_exch_data_wait
USE mo_fortran_tools,        ONLY: t_ptr_3d, t_ptr_3d_sp
USE mo_run_config,           ONLY: msg_level, activate_sync_timers
USE mo_decomposition_tools,  ONLY: t_glb2loc_index_lookup, get_local_index
USE mo_parallel_config,      ONLY: blk_no, idx_no, idx_1d
USE mo_communication_types,  ONLY: t_comm_pattern, t_p_comm_pattern, &
  &                                t_comm_pattern_collection, xfer_list
#ifdef _OPENACC
USE mo_mpi,                  ONLY: i_am_accel_node
#endif


IMPLICIT NONE

PRIVATE

!modules interface-------------------------------------------
PUBLIC :: t_comm_pattern_orig
PUBLIC :: t_comm_pattern_collection_orig
!
!variables

!--------------------------------------------------------------------------------------------------
!
TYPE, EXTENDS(t_comm_pattern) :: t_comm_pattern_orig

  PRIVATE

   ! Number of points we receive in communication,
   ! this is the same as recv_limits

   INTEGER :: n_recv  ! Number of points we receive from other PEs
   INTEGER :: n_pnts  ! Number of points we output into local array;
                      ! this may be bigger than n_recv due to
                      ! duplicate entries
   INTEGER :: n_send  ! Number of points we send to other PEs

   INTEGER :: np_recv ! Number of PEs from which data have to be received
   INTEGER :: np_send ! Number of PEs to which data have to be sent

   !> which communicator to apply this pattern to
   INTEGER :: comm

   ! "recv_limits":
   !
   ! All data that is received from PE np is buffered in the receive
   ! buffer between start index "p_pat%recv_limits(np)+1" and the end
   ! index "p_pat%recv_limits(np+1)".
   INTEGER, ALLOCATABLE :: recv_limits(:)

   ! "recv_src", "recv_dst_blk/idx":
   !
   ! For all points i=1,n_pnts the data received at index recv_src(i)
   ! in the receiver buffer is copied to the destination array at
   ! position recv_dst_idx/blk(i)
   INTEGER, ALLOCATABLE :: recv_src(:)
   INTEGER, ALLOCATABLE :: recv_dst_blk(:)
   INTEGER, ALLOCATABLE :: recv_dst_idx(:)

   ! "send_limits":
   !
   ! All data that is sent to PE np is buffered by the local PE in the
   ! send buffer between start index "p_pat%send_limits(np)+1" and the
   ! end index "p_pat%send_limits(np+1)".
   INTEGER, ALLOCATABLE :: send_limits(:)

   ! "send_src_idx/blk":
   !
   ! For all points i=1,n_send the data in the send buffer at the ith
   ! position is copied from the source array at position
   ! send_src_idx/blk(i)
   INTEGER, ALLOCATABLE :: send_src_blk(:)
   INTEGER, ALLOCATABLE :: send_src_idx(:)

   ! "pelist_send", "pelist_recv":
   !
   ! list of PEs where to send the data to, and from where to receive
   ! the data
   INTEGER, ALLOCATABLE :: pelist_send(:)
   INTEGER, ALLOCATABLE :: pelist_recv(:)

   ! "send_startidx", "send_count":
   !
   ! The local PE sends send_count(i) data items to PE pelist_send(i),
   ! starting at send_startidx(i) in the send buffer.
   INTEGER, ALLOCATABLE :: send_startidx(:)
   INTEGER, ALLOCATABLE :: send_count(:)

   ! "recv_startidx", "recv_count":
   !
   ! The local PE recvs recv_count(i) data items from PE pelist_recv(i),
   ! starting at recv_startidx(i) in the receiver buffer.
   INTEGER, ALLOCATABLE :: recv_startidx(:)
   INTEGER, ALLOCATABLE :: recv_count(:)

  CONTAINS

    PROCEDURE :: setup => setup_comm_pattern
    PROCEDURE :: setup2 => setup_comm_pattern2
    PROCEDURE :: delete => delete_comm_pattern
    PROCEDURE :: exchange_data_r3d => exchange_data_r3d
    PROCEDURE :: exchange_data_s3d => exchange_data_s3d
    PROCEDURE :: exchange_data_i3d => exchange_data_i3d
    PROCEDURE :: exchange_data_r2d => exchange_data_r2d
    PROCEDURE :: exchange_data_s2d => exchange_data_s2d
    PROCEDURE :: exchange_data_i2d => exchange_data_i2d
    PROCEDURE :: exchange_data_mult => exchange_data_mult
    PROCEDURE :: exchange_data_mult_mixprec => exchange_data_mult_mixprec
    PROCEDURE :: exchange_data_4de1 => exchange_data_4de1
    PROCEDURE :: get_np_recv => get_np_recv
    PROCEDURE :: get_np_send => get_np_send
    PROCEDURE :: get_pelist_recv => get_pelist_recv

END TYPE t_comm_pattern_orig

TYPE t_p_comm_pattern_orig
   TYPE(t_comm_pattern_orig), POINTER :: p
END TYPE t_p_comm_pattern_orig

TYPE, EXTENDS(t_comm_pattern_collection) :: t_comm_pattern_collection_orig

   PRIVATE

   TYPE(t_p_comm_pattern_orig), ALLOCATABLE :: patterns(:)

   CONTAINS

   PROCEDURE :: setup => setup_comm_pattern_collection
   PROCEDURE :: delete => delete_comm_pattern_collection
   PROCEDURE :: exchange_data_grf => exchange_data_grf

END TYPE t_comm_pattern_collection_orig

#if defined( _OPENACC )
#if defined(__COMMUNICATION_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
#endif

!
!------------------------------------------------------------------------------------------------
!

CHARACTER(*), PARAMETER :: modname = "mo_communication_orig"

  PUBLIC :: exchange_data_noblk
  INTERFACE exchange_data_noblk
    MODULE PROCEDURE exchange_data_r1d_2d
    MODULE PROCEDURE exchange_data_s1d_2d
    MODULE PROCEDURE exchange_data_i1d_2d
  END INTERFACE exchange_data_noblk

!-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !
  !

  !>
  !! Sets up a communication pattern for exchanging data.
  !!
  !! Note: This setup routine works only for the trivial communication
  !!       patterns in sequential runs.
  !!
  !! dst_n_points     Total number of points in the RECEIVER array,
  !!                  not every point is necessarily set during exchange
  !!                  (see owner!)
  !!
  !! dst_owner        Owner PE number of every point in the RECEIVER array,
  !!                  if owner(.) == -1, this point will not be set during exchange.
  !!                  If owner(.) == p_pe, this point will be exchanged,
  !!                  this is necessary if sender and receiver arrays are
  !!                  different (e.g. feedback, gather, scatter)
  !!
  !! dst_global_index Global index of of every point in the RECEIVER array
  !!                  There may be more than 1 point with the same global index,
  !!                  in this case the point is exchanged only once and
  !!                  locally distributed.
  !!                  - If this argument is not present, we assume global_index=1,2.3,...
  !! inplace          In-place data exchanges are allowed (source and destination
  !!                  arrays can be identically for the exchange)
  !!                  - if inplace == true, the user guarantees that
  !!                    (src_n_points == dst_n_points) and that points, which will
  !!                    have to be sent to other processes, are disjunct from the
  !!                    points that will have to be received
  !!                  - in case the user only provides the receive array to an
  !!                    exchange call and not send array, the exchange will be
  !!                    faster if inplace == true
  !!
  !! send_decomp_info domain decomposition information for the SENDER array
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  SUBROUTINE setup_comm_pattern(p_pat, dst_n_points, dst_owner, &
                                dst_global_index, send_glb2loc_index, &
                                src_n_points, src_owner, src_global_index, &
                                inplace, comm)

    CLASS(t_comm_pattern_orig), INTENT(OUT) :: p_pat

    INTEGER, INTENT(IN) :: dst_n_points        ! Total number of points
    INTEGER, INTENT(IN) :: dst_owner(:)        ! Owner of every point
    INTEGER, INTENT(IN) :: dst_global_index(:) ! Global index of every point
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: send_glb2loc_index
                                               ! global to local index
                                               ! lookup information
                                               ! of the SENDER array
    INTEGER, INTENT(IN) :: src_n_points        ! Total number of points
    INTEGER, INTENT(IN) :: src_owner(:)        ! Owner of every point
    INTEGER, INTENT(IN) :: src_global_index(:) ! Global index of every point

    LOGICAL, OPTIONAL, INTENT(IN) :: inplace
    INTEGER, OPTIONAL, INTENT(in) :: comm

    CHARACTER(len=*), PARAMETER :: routine = modname//"::setup_comm_pattern"
    INTEGER, ALLOCATABLE :: icnt(:), flag(:), global_recv_index(:), send_src(:), num_rcv(:)
    INTEGER :: i, n, np, nr, num_recv, irs, ire, num_send, iss, ise, max_glb, &
      comm_size, comm_rank, recv_idx, abs_dst_idx, n_pnts
    LOGICAL :: any_np_le_0

    !-----------------------------------------------------------------------
    IF (PRESENT(comm)) THEN
      p_pat%comm = comm
    ELSE
      p_pat%comm = p_comm_work
    END IF

    comm_size = p_comm_size(p_pat%comm)
    comm_rank = p_comm_rank(p_pat%comm)
    ALLOCATE(icnt(0:comm_size-1), num_rcv(0:comm_size-1))
    max_glb = MAX(MAXVAL(ABS(dst_global_index(1:dst_n_points)), &
      &                  mask=(dst_owner(1:dst_n_points)>=0)),1)
    ALLOCATE(flag(max_glb))

    ! Count the number of points we want to receive from every PE
    ! and the total number of points to output

    icnt(:) = 0
    flag(:) = 0

    n_pnts = 0

    DO i = 1, dst_n_points
      IF (dst_owner(i) >= 0) THEN
        n_pnts = n_pnts + 1 ! Count total number of points we output
        abs_dst_idx = ABS(dst_global_index(i))
        IF (flag(abs_dst_idx)==0) THEN
          icnt(dst_owner(i)) = icnt(dst_owner(i))+1 ! Number to get from dst_owner(i)
          flag(abs_dst_idx) = 1 ! Flag that this global point is already on the list
        ENDIF
      ENDIF
    ENDDO
    p_pat%n_pnts = n_pnts

    ! Allocate and set up the recv_limits array

    ALLOCATE(p_pat%recv_limits(0:comm_size))

    i = 0
    DO np = 0, comm_size - 1
      p_pat%recv_limits(np) = i
      i = i + icnt(np)
    ENDDO
    p_pat%recv_limits(comm_size) = i

    ! The last entry in recv_limits is the total number of points we receive
    p_pat%n_recv = i

    ! Allocate and set up the recv_src array

    ALLOCATE(p_pat%recv_src(n_pnts))
    ALLOCATE(p_pat%recv_dst_blk(n_pnts))
    ALLOCATE(p_pat%recv_dst_idx(n_pnts))
    ALLOCATE(global_recv_index(p_pat%n_recv))

    DO np = 0, comm_size-1
      icnt(np) = p_pat%recv_limits(np)
    ENDDO

    flag(:) = 0
    n = 0 ! Counts total number of local points

    DO i = 1, dst_n_points
      IF(dst_owner(i)>=0) THEN
        n = n+1
        abs_dst_idx = ABS(dst_global_index(i))
        recv_idx = flag(abs_dst_idx)
        IF (recv_idx == 0) THEN
          recv_idx = icnt(dst_owner(i)) + 1
          icnt(dst_owner(i)) = recv_idx   ! Current index in recv array
          ! Global index of points in receive array
          global_recv_index(recv_idx) = abs_dst_idx
          flag(abs_dst_idx) = recv_idx    ! Store from where to get duplicates
        ENDIF
        p_pat%recv_src(n) = recv_idx      ! From where in the receive array
                                          ! this process receives the local
                                          ! point
        p_pat%recv_dst_blk(n) = blk_no(i) ! Where to put the local point
        p_pat%recv_dst_idx(n) = idx_no(i) ! Where to put the local point
      ENDIF
    ENDDO


    ! Exchange the number of points we want to receive with the respective senders
    DO np = 0, comm_size-1 ! loop over PEs where to send the data
      num_rcv(np) = p_pat%recv_limits(np+1) - p_pat%recv_limits(np)
    ENDDO

    CALL p_alltoall(num_rcv, icnt, p_pat%comm)
    ! Now send the global index of the points we need from PE np
    DO np = 0, comm_size-1 ! loop over PEs where to send the data


      irs = p_pat%recv_limits(np)+1 ! Start index in global_recv_index
      ire = p_pat%recv_limits(np+1) ! End   index in global_recv_index

      IF (np /= comm_rank .AND. num_rcv(np) > 0) &
        CALL p_isend(global_recv_index(irs), np, 1, &
        &            p_count=ire-irs+1, comm=p_pat%comm)

    ENDDO

    DEALLOCATE(num_rcv)
    ! Allocate and set up the send_limits array
    ALLOCATE(p_pat%send_limits(0:comm_size))

    p_pat%send_limits(0) = 0
    i = 0
    DO nr = 0, comm_size-1
      p_pat%send_limits(nr) = i
      i = i + icnt(nr)
    ENDDO
    p_pat%send_limits(comm_size) = i

    ! The last entry in send_limits is the total number of points we receive
    p_pat%n_send = i

    ! Allocate and set up the send_src array

    ALLOCATE(send_src(p_pat%n_send))
    DO nr = 0, comm_size-1
      iss = p_pat%send_limits(nr)   ! Start index in send_src
      ise = p_pat%send_limits(nr+1) ! End   index in send_src
      num_send = ise - iss
      iss = iss + 1
      IF (num_send > 0) THEN
        IF (nr /= comm_rank) THEN
          CALL p_recv(send_src(iss), nr, 1, p_count=ise-iss+1, comm=p_pat%comm)
        ELSE
          irs = p_pat%recv_limits(comm_rank)+1   ! Start index in global_recv_index
          ire = p_pat%recv_limits(comm_rank+1)   ! End   index in global_recv_index
          send_src(iss:ise) = global_recv_index(irs:ire)
        ENDIF
      END IF
    ENDDO

    CALL p_wait

    ALLOCATE(p_pat%send_src_blk(p_pat%n_send))
    ALLOCATE(p_pat%send_src_idx(p_pat%n_send))

    ! The indices in p_pat%send_src are global, convert to local

    any_np_le_0 = .FALSE.
    DO i = 1, p_pat%n_send

      np = get_local_index(send_glb2loc_index, send_src(i))
      IF (np <= 0) THEN
        WRITE (0, '(3(a,i0))') 'problem at i=', i, &
          ', send_src(i)=', send_src(i), &
          ', np=', np
      END IF
      any_np_le_0 = any_np_le_0 .OR. np <= 0
      p_pat%send_src_blk(i) = blk_no(np)
      p_pat%send_src_idx(i) = idx_no(np)
    ENDDO
    IF (any_np_le_0) THEN
      WRITE (0, '(a)') 'send_glb2loc_index%outer_glb_index'
      WRITE (0, '(10(i0," "))') send_glb2loc_index%outer_glb_index
      WRITE (0, '(a)') 'send_glb2loc_index%inner_glb_index'
      WRITE (0, '(10(i0," "))') send_glb2loc_index%inner_glb_index
      FLUSH(0)
      CALL finish(routine, 'Got illegal index')
    END IF
    ! Finally, compute lists of processors for send and receive operations

    num_send = 0
    num_recv = 0

    DO np = 0, comm_size-1 ! loop over PEs

      iss = p_pat%send_limits(np)+1
      ise = p_pat%send_limits(np+1)
      IF(ise >= iss) num_send = num_send + 1

      irs = p_pat%recv_limits(np)+1
      ire = p_pat%recv_limits(np+1)
      IF(ire >= irs) num_recv = num_recv + 1

    ENDDO

    p_pat%np_send = num_send
    p_pat%np_recv = num_recv

    ALLOCATE (p_pat%pelist_send(num_send), p_pat%pelist_recv(num_recv),     &
      p_pat%send_startidx(num_send), p_pat%recv_startidx(num_recv), &
      p_pat%send_count(num_send), p_pat%recv_count(num_recv)        )

    num_send = 0
    num_recv = 0

    DO np = 0, comm_size-1 ! loop over PEs

      iss = p_pat%send_limits(np)+1
      ise = p_pat%send_limits(np+1)
      IF(ise >= iss) THEN
        num_send = num_send + 1
        p_pat%pelist_send(num_send)   = np
        p_pat%send_startidx(num_send) = iss
        p_pat%send_count(num_send)    = ise - iss + 1
      ENDIF

      irs = p_pat%recv_limits(np)+1
      ire = p_pat%recv_limits(np+1)
      IF(ire >= irs) THEN
        num_recv = num_recv + 1
        p_pat%pelist_recv(num_recv)   = np
        p_pat%recv_startidx(num_recv) = irs
        p_pat%recv_count(num_recv)    = ire - irs + 1
      ENDIF

    ENDDO

    DEALLOCATE(icnt, flag, global_recv_index, send_src)

    ! consistency check of communication pattern
#ifndef NOMPI
    IF (msg_level >= 25)  &
      CALL check_comm_pattern(p_pat)
#endif

  END SUBROUTINE setup_comm_pattern


  !-------------------------------------------------------------------------


  SUBROUTINE setup_comm_pattern2(p_pat, comm, recv_msg, send_msg, &
       glb2loc_index_recv, glb2loc_index_send, inplace)
    CLASS(t_comm_pattern_orig), INTENT(out) :: p_pat
    INTEGER, INTENT(in) :: comm
    TYPE(xfer_list), INTENT(in) :: recv_msg(:), send_msg(:)
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: glb2loc_index_recv, &
         glb2loc_index_send
    LOGICAL, OPTIONAL, INTENT(in) :: inplace
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: recv_msg, send_msg
#endif
    INTEGER :: np_recv, np_send, n_send, n_recv, i, comm_rank, &
         n_pnts_recv, n_pnts_send, comm_size
    LOGICAL :: is_inter
    CHARACTER(len=*), PARAMETER :: routine &
         = 'mo_communication::setup_comm_pattern2'

    p_pat%comm = comm
    comm_rank = p_comm_rank(comm)
    is_inter = p_comm_is_intercomm(comm)
    IF (is_inter) THEN
      comm_size = p_comm_remote_size(comm)
    ELSE
      comm_size = p_comm_size(comm)
    END IF
    np_send = SIZE(send_msg)
    p_pat%np_send = np_send
    np_recv = SIZE(recv_msg)
    p_pat%np_recv = np_recv
    ALLOCATE(p_pat%recv_limits(0:comm_size), p_pat%send_limits(0:comm_size), &
         p_pat%pelist_send(np_send), p_pat%pelist_recv(np_recv), &
         p_pat%send_startidx(np_send), p_pat%send_count(np_send), &
         p_pat%recv_startidx(np_recv), p_pat%recv_count(np_recv))
    CALL count_msg_size(np_recv, recv_msg, n_pnts_recv, n_recv, &
         is_inter, comm_rank, p_pat%pelist_recv, p_pat%recv_count, &
         p_pat%recv_startidx)
    CALL count_msg_size(np_send, send_msg, n_pnts_send, n_send, &
         is_inter, comm_rank, p_pat%pelist_send, p_pat%send_count, &
         p_pat%send_startidx)
    p_pat%n_send = n_send
    IF (n_pnts_recv /= n_pnts_send) &
      CALL finish(routine, "inconsistent lists")
    n_pnts_recv = n_pnts_recv + n_recv
    p_pat%n_pnts = n_pnts_recv
    p_pat%n_recv = n_recv
    ALLOCATE(p_pat%recv_src(n_pnts_recv), &
      &      p_pat%recv_dst_blk(n_pnts_recv), p_pat%recv_dst_idx(n_pnts_recv), &
      &      p_pat%send_src_blk(n_send), p_pat%send_src_idx(n_send))
    CALL list2limits(comm_size, p_pat%recv_limits, recv_msg, is_inter)
    CALL list2limits(comm_size, p_pat%send_limits, send_msg, is_inter)
    CALL expand1d2blkidx(recv_msg, p_pat%recv_limits, p_pat%recv_src)
    CALL expandglb2blkidx(recv_msg, glb2loc_index_recv, p_pat%recv_limits, &
         p_pat%recv_dst_blk, p_pat%recv_dst_idx)
    CALL expandglb2blkidx(send_msg, glb2loc_index_send, p_pat%send_limits, &
         p_pat%send_src_blk, p_pat%send_src_idx)
  CONTAINS
    SUBROUTINE count_msg_size(nmsg, msg, nself, nremote, is_inter, comm_rank, &
         ranks, counts, starts)
      INTEGER, INTENT(in) :: nmsg, comm_rank
      TYPE(xfer_list), INTENT(in) :: msg(nmsg)
      INTEGER, INTENT(out) :: nself, nremote, ranks(nmsg), counts(nmsg), &
           starts(nmsg)
      LOGICAL, INTENT(in) :: is_inter
      INTEGER :: nidx_remote, nidx_local, sz, msg_rank, sz_accum
      nidx_remote = 0
      nidx_local = 0
      sz_accum = 1
      DO i = 1, nmsg
        msg_rank = msg(i)%rank
        starts(i) = sz_accum
        sz = SIZE(msg(i)%glob_idx)
        sz_accum = sz_accum + sz
        ranks(i) = msg_rank
        counts(i) = sz
        IF (is_inter .OR. msg_rank /= comm_rank) THEN
          nidx_remote = nidx_remote + sz
        ELSE
          nidx_local = nidx_local + sz
        END IF
      END DO
      nremote = nidx_remote
      nself = nidx_local
    END SUBROUTINE count_msg_size

    SUBROUTINE list2limits(comm_size, limits, msg, is_inter)
      INTEGER, INTENT(in) :: comm_size
      INTEGER, INTENT(out) :: limits(0:comm_size)
      TYPE(xfer_list), INTENT(in) :: msg(:)
      LOGICAL, INTENT(in) :: is_inter
      INTEGER :: i, msg_rank, nmsg, limits_psum
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      CONTIGUOUS :: msg
#endif
      nmsg = SIZE(msg)
      limits(0:comm_size) = 0
      DO i = 1, nmsg
        msg_rank = msg(i)%rank
        limits(msg_rank+1) = limits(msg_rank+1) + SIZE(msg(i)%glob_idx)
      END DO
      limits_psum = 0
      DO i = 1, comm_size
        limits_psum = limits_psum + limits(i)
        limits(i) = limits_psum
      END DO
    END SUBROUTINE list2limits

    SUBROUTINE expand1d2blkidx(msg, limits, recv_src)
      TYPE(xfer_list), INTENT(in) :: msg(:)
      INTEGER, INTENT(in) :: limits(0:)
      INTEGER, INTENT(out) :: recv_src(:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      CONTIGUOUS :: msg, limits, recv_src
#endif
      INTEGER :: nmsg, i, sz, sz_psum, msg_rank, jls, jl, jle
      nmsg = SIZE(msg)
      DO i = 1, nmsg
        msg_rank = msg(i)%rank
        jle = limits(msg_rank+1)
        jls = limits(msg_rank)
        sz = jle - jls
        DO jl = 1, sz
          recv_src(jls+jl) = jls + jl
        END DO
      END DO
    END SUBROUTINE expand1d2blkidx

    SUBROUTINE  expandglb2blkidx(msg, glb2loc_index, limits, &
         a_blk, a_idx)
      TYPE(xfer_list), INTENT(in) :: msg(:)
      TYPE(t_glb2loc_index_lookup), INTENT(IN) :: glb2loc_index
      INTEGER, INTENT(in) :: limits(0:)
      INTEGER, INTENT(out) :: a_blk(:), a_idx(:)
      INTEGER :: i, jl, msg_rank, sz, jls, jle, np, nmsg
      LOGICAL :: any_np_le_0
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      CONTIGUOUS :: limits, msg, a_blk, a_idx
#endif
      nmsg = SIZE(msg)
      any_np_le_0 = .FALSE.
      DO i = 1, nmsg
        msg_rank = msg(i)%rank
        jle = limits(msg_rank+1)
        jls = limits(msg_rank)
        sz = jle - jls
        DO jl = 1, sz
          np = get_local_index(glb2loc_index, msg(i)%glob_idx(jl))
          any_np_le_0 = any_np_le_0 .OR. np <= 0
          a_blk(jls+jl) = blk_no(np)
          a_idx(jls+jl) = idx_no(np)
        END DO
      END DO
      IF (any_np_le_0) CALL finish('setup_comm_pattern2','Got illegal index')
    END SUBROUTINE expandglb2blkidx
  END SUBROUTINE setup_comm_pattern2


  !-------------------------------------------------------------------------


  SUBROUTINE setup_comm_pattern_collection(pattern_collection, patterns)

    CLASS(t_comm_pattern_collection_orig), INTENT(OUT) :: pattern_collection
    TYPE(t_p_comm_pattern), INTENT(IN) :: patterns(:)

    INTEGER :: i

    ALLOCATE(pattern_collection%patterns(SIZE(patterns)))

    DO i = 1, SIZE(patterns)
      SELECT TYPE (pattern_orig => patterns(i)%p)
        TYPE is (t_comm_pattern_orig)
          pattern_collection%patterns(i)%p => pattern_orig
        CLASS DEFAULT
          CALL finish("setup_comm_pattern_collection", &
                      "wrong t_comm_pattern type")
      END SELECT
    END DO

  END SUBROUTINE setup_comm_pattern_collection


  !-------------------------------------------------------------------------
  !
  !>
  !! Deletes a communication pattern, i.e. deallocates all arrays
  !! and sets all other members to 0
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Oct 2011
  !!
  !
  SUBROUTINE delete_comm_pattern(p_pat)

    CLASS(t_comm_pattern_orig), INTENT(INOUT) :: p_pat

    ! deallocate arrays

    IF(ALLOCATED(p_pat%recv_limits))   DEALLOCATE(p_pat%recv_limits)

    IF(ALLOCATED(p_pat%recv_src))      DEALLOCATE(p_pat%recv_src)
    IF(ALLOCATED(p_pat%recv_dst_blk))  DEALLOCATE(p_pat%recv_dst_blk)
    IF(ALLOCATED(p_pat%recv_dst_idx))  DEALLOCATE(p_pat%recv_dst_idx)

    IF(ALLOCATED(p_pat%send_limits))   DEALLOCATE(p_pat%send_limits)

    IF(ALLOCATED(p_pat%send_src_blk))  DEALLOCATE(p_pat%send_src_blk)
    IF(ALLOCATED(p_pat%send_src_idx))  DEALLOCATE(p_pat%send_src_idx)

    IF(ALLOCATED(p_pat%pelist_send))   DEALLOCATE(p_pat%pelist_send)
    IF(ALLOCATED(p_pat%pelist_recv))   DEALLOCATE(p_pat%pelist_recv)

    IF(ALLOCATED(p_pat%send_startidx)) DEALLOCATE(p_pat%send_startidx)
    IF(ALLOCATED(p_pat%recv_startidx)) DEALLOCATE(p_pat%recv_startidx)

    IF(ALLOCATED(p_pat%send_count))    DEALLOCATE(p_pat%send_count)
    IF(ALLOCATED(p_pat%recv_count))    DEALLOCATE(p_pat%recv_count)

    ! Set other members to 0

    p_pat%n_recv  = 0
    p_pat%n_pnts  = 0
    p_pat%n_send  = 0
    p_pat%np_recv = 0
    p_pat%np_send = 0

  END SUBROUTINE delete_comm_pattern


  !-------------------------------------------------------------------------


  SUBROUTINE delete_comm_pattern_collection(pattern_collection)

    CLASS(t_comm_pattern_collection_orig), INTENT(INOUT) :: pattern_collection

    INTEGER :: i

    DO i = 1, SIZE(pattern_collection%patterns)
      CALL pattern_collection%patterns(i)%p%delete()
      DEALLOCATE(pattern_collection%patterns(i)%p)
    END DO
    DEALLOCATE(pattern_collection%patterns)

  END SUBROUTINE delete_comm_pattern_collection


  !-------------------------------------------------------------------------
  !> Consistency check of communication pattern.
  !! Sends pattern info to working PE 0, which checks this data
  !! for consistency wrt. send/receive counts.
  !!
  !! @par Revision History
  !! Initial implementation  by  F. Prill, DWD (2012-01-20)
  !!
  SUBROUTINE check_comm_pattern(p_pat)
    CLASS(t_comm_pattern_orig), INTENT(INOUT) :: p_pat

    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":check_comm_pattern"
    INTEGER, PARAMETER      :: root    = 0
    INTEGER                 :: ierr, npes, i_pe, i_target_pe, ntarget_pes,        &
      &                        target_pe, target_cnt, this_pe, i_source_pe,       &
      &                        nsource_pes, source_pe, source_cnt, p_comm
    INTEGER, ALLOCATABLE    :: recvbuf_send(:), displs(:),                        &
      &                        recvbuf_recv(:), recvbuf_scnt(:), recvbuf_rcnt(:), &
      &                        recvbuf_npsnd(:), recvbuf_nprcv(:)
    LOGICAL                 :: lcheck, lfound_peer

    CALL message(routine, "Consistency check of communication pattern.")
    this_pe = get_my_mpi_work_id()
    p_comm  = get_my_mpi_work_communicator()

    !-- allocate memory -------------------------
    npes = get_my_mpi_work_comm_size()
    ALLOCATE(recvbuf_send(npes*npes), displs(npes), recvbuf_recv(npes*npes), &
      &      recvbuf_scnt(npes*npes), recvbuf_rcnt(npes*npes),               &
      &      recvbuf_npsnd(npes), recvbuf_nprcv(npes),                       &
      &      STAT=ierr)
    IF (ierr /= 0) CALL finish (routine, 'Error in ALLOCATE operation!')

    !-- gather fields at work PE 0 --------------

    ! gather numbers of send/receiver partners
    CALL p_gather(p_pat%np_send, recvbuf_npsnd, root, p_comm)
    CALL p_gather(p_pat%np_recv, recvbuf_nprcv, root, p_comm)

    ! set displacements array
    displs(:)     = (/ ( (i_pe-1)*npes, i_pe=1, npes) /)

    ! field 1: list of target PEs
    CALL p_gatherv(p_pat%pelist_send, p_pat%np_send,    &
      &            recvbuf_send, recvbuf_npsnd, displs, &
      &            root, p_comm)
    ! field 2: list of source PEs
    CALL p_gatherv(p_pat%pelist_recv, p_pat%np_recv,    &
      &            recvbuf_recv, recvbuf_nprcv, displs, &
      &            root, p_comm)
    ! field 3: list of target counts
    CALL p_gatherv(p_pat%send_count, p_pat%np_send,     &
      &            recvbuf_scnt, recvbuf_npsnd, displs, &
      &            root, p_comm)
    ! field 4: list of source counts
    CALL p_gatherv(p_pat%recv_count, p_pat%np_recv,     &
      &            recvbuf_rcnt, recvbuf_nprcv, displs, &
      &            root, p_comm)

    !-- perform consistency checks --------------
    lcheck = .TRUE.
    IF (this_pe == root) THEN
      ! now loop over PEs and check if their send counts match the info
      ! on the receiver side:
      DO i_pe=1, npes
        ntarget_pes = recvbuf_npsnd(i_pe)
        DO i_target_pe=1,ntarget_pes
          target_pe  = recvbuf_send(displs(i_pe) + i_target_pe)
          target_cnt = recvbuf_scnt(displs(i_pe) + i_target_pe)
          ! loop over target PE's receiver array
          lfound_peer = .FALSE.
          nsource_pes = recvbuf_nprcv(target_pe+1)
          DO i_source_pe=1, nsource_pes
            source_pe  = recvbuf_recv(displs(target_pe+1) + i_source_pe)
            IF ((source_pe+1) == i_pe) THEN
              lfound_peer = .TRUE.
              source_cnt = recvbuf_rcnt(displs(target_pe+1) + i_source_pe)
              IF (source_cnt /= target_cnt) THEN
                WRITE (message_text,*) "PE ", i_pe-1, ": sends ", target_cnt, &
                  &         " values to ", target_pe, " but PE ", target_pe,  &
                  &         " receives ", source_cnt, " values from ", source_pe
                CALL message(routine, message_text)
                lcheck = .FALSE.
              END IF
            END IF
          END DO ! i_source_pe
          IF (.NOT. lfound_peer) THEN
            WRITE (message_text,*) "PE ", i_pe-1, ": Missing peer!"
            CALL message(routine, message_text)
            lcheck = .FALSE.
          ELSE
            IF (msg_level >= 25) THEN
              WRITE (message_text,*) "PE ", i_pe-1, ": sends ", target_cnt, &
                &                    " values to ", target_pe
              CALL message(routine, message_text)
            END IF
          END IF
        END DO ! i_target_pe
      END DO ! i_pe

      ! loop over PEs and check if their receive counts match the info
      ! on the sender side:
      DO i_pe=1, npes
        nsource_pes = recvbuf_nprcv(i_pe)
        DO i_source_pe=1,nsource_pes
          source_pe  = recvbuf_recv(displs(i_pe) + i_source_pe)
          source_cnt = recvbuf_rcnt(displs(i_pe) + i_source_pe)
          ! loop over source PE's sender array
          lfound_peer = .FALSE.
          ntarget_pes = recvbuf_npsnd(source_pe+1)
          DO i_target_pe=1, ntarget_pes
            target_pe  = recvbuf_send(displs(source_pe+1) + i_target_pe)
            IF ((target_pe+1) == i_pe) THEN
              lfound_peer = .TRUE.
              target_cnt = recvbuf_scnt(displs(source_pe+1) + i_target_pe)
              IF (source_cnt /= target_cnt) THEN
                WRITE (message_text,*) "PE ", i_pe-1, ": receives ", source_cnt, &
                  &         " values from ", source_pe, " but PE ", source_pe,   &
                  &         " sends ", target_cnt, " values to ", target_pe
                CALL message(routine, message_text)
                lcheck = .FALSE.
              END IF
            END IF
          END DO ! i_target_pe
          IF (.NOT. lfound_peer) THEN
            WRITE (message_text,*) "PE ", i_pe-1, ": Missing peer!"
            CALL message(routine, message_text)
            lcheck = .FALSE.
          ELSE
            IF (msg_level >= 25) THEN
              WRITE (message_text,*) "PE ", i_pe-1, ": receives ", source_cnt, &
                &                    " values from ", source_pe
              CALL message(routine, message_text)
            END IF
          END IF
        END DO ! i_source_pe
      END DO ! i_pe
    END IF ! (this_pe == root)

    ! clean up
    DEALLOCATE(recvbuf_send, displs,                      &
      &        recvbuf_recv, recvbuf_scnt, recvbuf_rcnt,  &
      &        recvbuf_npsnd, recvbuf_nprcv,              &
      &        STAT=ierr)
    IF (ierr /= 0) CALL finish (routine, 'Error in DEALLOCATE operation!')

    ! abort in case of inconsistencies:
    IF (.NOT. lcheck) &
      CALL finish(routine, "Inconsistencies detected in communication pattern!")

    CALL message(routine, "Done.")

  END SUBROUTINE check_comm_pattern

  
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Modified by Guenther Zaengl for vectorization
  !!
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_r3d(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), INTENT(INOUT) :: p_pat
    REAL(dp), INTENT(INOUT), TARGET           :: recv(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET    :: send(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET    :: add (:,:,:)

    REAL(dp)  :: send_buf(SIZE(recv,2),p_pat%n_send), recv_buf(SIZE(recv,2),p_pat%n_recv)
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif
    REAL(dp), POINTER :: send_ptr(:,:,:)

    INTEGER :: i, k, np, irs, iss, pid, icount, ndim2

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_r3d_seq(p_pat, recv, send, add)
      RETURN
    END IF

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      CALL timer_start(timer_barrier)
      CALL p_barrier(p_pat%comm)
      CALL timer_stop(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    IF(SIZE(recv,1) /= nproma) THEN
      CALL finish('exchange_data_r3d','Illegal first dimension of data array')
    ENDIF

    ndim2 = SIZE(recv,2)

#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
#endif

!$ACC DATA CREATE( send_buf, recv_buf ), &
!$ACC      COPYIN( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src ), &
!$ACC      IF ( i_am_accel_node .AND. acc_on )

    IF (iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) THEN
      ! Set up irecv's for receive buffers
      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2
        CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF

    ! Set up send buffer

    IF(PRESENT(send)) THEN
      send_ptr => send
    ELSE
      send_ptr => recv
    ENDIF

    IF (ndim2 == 1) THEN
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
      DO i = 1, n_send
        send_buf(1,i) = send_ptr(tmp_send_src_idx(i),1,tmp_send_src_blk(i))
      ENDDO
!$ACC END PARALLEL
#else
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i))
      ENDDO
#endif
    ELSE

#if defined( __SX__ ) || defined( _OPENACC )
!CDIR UNROLL=6
!$ACC PARALLEL PRESENT( send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO k = 1, ndim2
!$ACC LOOP VECTOR
        DO i = 1, n_send
          send_buf(k,i) = send_ptr(tmp_send_src_idx(i),k,tmp_send_src_blk(i))
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
      DO i = 1, p_pat%n_send
        send_buf(1:ndim2,i) = send_ptr(p_pat%send_src_idx(i),1:ndim2,   &
          &                            p_pat%send_src_blk(i))
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
    ENDIF

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    ! Send our data
    IF (iorder_sendrecv == 1) THEN
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2
        CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2
        CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO

      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2
        CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ELSE IF (iorder_sendrecv == 3) THEN ! use irecv/isend
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2
        CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF

    IF (iorder_sendrecv > 0) THEN
      ! Wait for all outstanding requests to finish
      start_sync_timer(timer_exch_data_wait)
      CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif
      stop_sync_timer(timer_exch_data_wait)

    ELSE IF (iorder_sendrecv == 0) THEN
      ! dummy mode; just fill the receive buffer with "something reasonable"
      IF (p_pat%n_send > 0) THEN
!$ACC KERNELS PRESENT( send_buf, recv_buf ), IF( i_am_accel_node .AND. acc_on )
        DO k = 1, ndim2
          recv_buf(k,:) = send_buf(k,1)
        ENDDO
!$ACC END KERNELS
      ELSE
!$ACC KERNELS PRESENT( recv_buf ), IF( i_am_accel_node .AND. acc_on )
        DO k = 1, ndim2
          recv_buf(k,:) = 1._dp
        ENDDO
!$ACC END KERNELS
      ENDIF
    ENDIF

    IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer

    IF(PRESENT(add)) THEN
      IF (ndim2 == 1) THEN
        k = 1
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( add, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
        DO i = 1, n_pnts
          recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) =   &
            recv_buf(k,tmp_recv_src(i)) +                       &
            add(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i))
        ENDDO
!$ACC END PARALLEL
#else
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
            recv_buf(k,p_pat%recv_src(i)) +                       &
            add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
        ENDDO
#endif
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( recv_buf, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, n_pnts
            recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) =   &
              recv_buf(k,tmp_recv_src(i)) +                       &
              add(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) =   &
            recv_buf(:,p_pat%recv_src(i)) +                       &
            add(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i))
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
      ENDIF
    ELSE
      IF (ndim2 == 1) THEN
        k = 1
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
        DO i = 1, n_pnts
          recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) = &
            recv_buf(k,tmp_recv_src(i))
        ENDDO
!$ACC END PARALLEL
#else
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
            recv_buf(k,p_pat%recv_src(i))
        ENDDO
#endif
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, n_pnts
            recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) = &
              recv_buf(k,tmp_recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) = &
            recv_buf(:,p_pat%recv_src(i))
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
      ENDIF
    ENDIF

!$ACC END DATA
    stop_sync_timer(timer_exch_data)

  END SUBROUTINE exchange_data_r3d


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Modified by Guenther Zaengl for vectorization
  !!
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_s3d(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), INTENT(INOUT) :: p_pat
    REAL(sp), INTENT(INOUT), TARGET           :: recv(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET    :: send(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET    :: add (:,:,:)

    REAL(sp) :: send_buf(SIZE(recv,2),p_pat%n_send), &
      recv_buf(SIZE(recv,2),p_pat%n_recv)
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif
    REAL(sp), POINTER :: send_ptr(:,:,:)

    INTEGER :: i, k, np, irs, iss, pid, icount, ndim2

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_s3d_seq(p_pat, recv, send, add)
      RETURN
    END IF

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    IF(SIZE(recv,1) /= nproma) THEN
      CALL finish('exchange_data_s3d','Illegal first dimension of data array')
    ENDIF

    ndim2 = SIZE(recv,2)

#ifdef _OPENACC
! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
#endif

!$ACC DATA CREATE( send_buf, recv_buf ), &
!$ACC      COPYIN( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src ), &
!$ACC      IF ( i_am_accel_node .AND. acc_on )


    IF (iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) THEN
      ! Set up irecv's for receive buffers
      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2
        CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF

    ! Set up send buffer

    IF(PRESENT(send)) THEN
      send_ptr => send
    ELSE
      send_ptr => recv
    ENDIF

    IF (ndim2 == 1) THEN
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
      DO i = 1, n_send
        send_buf(1,i) = send_ptr(tmp_send_src_idx(i),1,tmp_send_src_blk(i))
      ENDDO
!$ACC END PARALLEL
#else
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i))
      ENDDO
#endif
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!CDIR UNROLL=6
!$ACC PARALLEL PRESENT( send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO k = 1, ndim2
!$ACC LOOP VECTOR
        DO i = 1,n_send
          send_buf(k,i) = send_ptr(tmp_send_src_idx(i),k,tmp_send_src_blk(i))
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
      DO i = 1, p_pat%n_send
        send_buf(1:ndim2,i) = send_ptr(p_pat%send_src_idx(i),1:ndim2,   &
          &                            p_pat%send_src_blk(i))
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
    ENDIF

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    ! Send our data
    IF (iorder_sendrecv == 1) THEN
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2
        CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2
        CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO

      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2
        CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ELSE IF (iorder_sendrecv == 3) THEN ! use irecv/isend
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2
        CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF

    IF (iorder_sendrecv > 0) THEN
      ! Wait for all outstanding requests to finish
      start_sync_timer(timer_exch_data_wait)
      CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif
      stop_sync_timer(timer_exch_data_wait)

    ELSE IF (iorder_sendrecv == 0) THEN
      ! dummy mode; just fill the receive buffer with "something reasonable"
      IF (p_pat%n_send > 0) THEN
!$ACC KERNELS PRESENT( send_buf, recv_buf ), IF( i_am_accel_node .AND. acc_on )
        DO k = 1, ndim2
          recv_buf(k,:) = send_buf(k,1)
        ENDDO
!$ACC END KERNELS
      ELSE
!$ACC KERNELS PRESENT( recv_buf ), IF( i_am_accel_node .AND. acc_on )
        DO k = 1, ndim2
          recv_buf(k,:) = 1._sp
        ENDDO
!$ACC END KERNELS
      ENDIF
    ENDIF

    IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer

    IF(PRESENT(add)) THEN
      IF (ndim2 == 1) THEN
        k = 1
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( add, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
        DO i = 1, n_pnts
          recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) =   &
            recv_buf(k,tmp_recv_src(i)) +                       &
            add(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i))
        ENDDO
!$ACC END PARALLEL
#else
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
            recv_buf(k,p_pat%recv_src(i)) +                       &
            add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
        ENDDO
#endif
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( recv_buf, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, n_pnts
            recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) =   &
              recv_buf(k,tmp_recv_src(i)) +                       &
              add(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) =   &
            recv_buf(:,p_pat%recv_src(i)) +                       &
            add(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i))
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
      ENDIF
    ELSE
      IF (ndim2 == 1) THEN
        k = 1
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
        DO i = 1, n_pnts
          recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) = &
            recv_buf(k,tmp_recv_src(i))
        ENDDO
!$ACC END PARALLEL
#else
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
            recv_buf(k,p_pat%recv_src(i))
        ENDDO
#endif
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, n_pnts
            recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) = &
              recv_buf(k,tmp_recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) = &
            recv_buf(:,p_pat%recv_src(i))
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
      ENDIF
    ENDIF

!$ACC END DATA

    stop_sync_timer(timer_exch_data)

  END SUBROUTINE exchange_data_s3d


  ! SEQUENTIAL version of subroutine "exchange_data_r3d"
  !
  SUBROUTINE exchange_data_r3d_seq(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_r3d_seq"
    INTEGER :: i, ndim2
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif

    ! consistency checks
    ! ------------------

    ! make sure that we are in sequential mode
    IF (.NOT. my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Internal error: sequential routine called in parallel run!")
    END IF
    ! further tests
    IF ( (p_pat%np_recv /= 1) .OR. (p_pat%np_send /= 1) ) THEN
      CALL finish(routine, "Internal error: inconsistent no. send/receive peers!")
    END IF
    IF ( (p_pat%recv_limits(1) - p_pat%recv_limits(0)) /= (p_pat%send_limits(1) - p_pat%send_limits(0)) ) THEN
      CALL finish(routine, "Internal error: inconsistent sender/receiver size!")
    END IF
    IF ( (p_pat%recv_limits(0) /= 0) .OR. (p_pat%send_limits(0) /= 0) ) THEN
      CALL finish(routine, "Internal error: inconsistent sender/receiver start position!")
    END IF
    IF ( (p_pat%recv_limits(1) /= p_pat%n_recv) .OR. (p_pat%n_recv /= p_pat%n_send) ) THEN
      CALL finish(routine, "Internal error: inconsistent counts for sender/receiver!")
    END IF

    ! "communication" (direct copy)
    ! -----------------------------

    ndim2 = SIZE(recv,2)

    ! The next piece of code is a condensed version of the following
    ! (under the assumptions asserted above):
    !
    !     ! fill sender buffer
    !     DO i=1,n_send
    !       send_buf(i) = array_in(send_src_idx(i), send_src_blk(i))
    !     END DO
    !     ! copy sender to receiver buffer
    !     recv_buf(p_pat%recv_limits(0)+1:p_pat%recv_limits(1)) = &
    !       &  send_buf(p_pat%send_limits(0)+1:p_pat%send_limits(1))
    !     ! copy from receiver buffer
    !     DO i=1,n_pnts
    !       array_out( recv_dst_idx(i), recv_dst_blk(i) ) = recv_buf(recv_src(i))
    !     END DO

#ifdef _OPENACC
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
!$ACC UPDATE DEVICE( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src ), IF ( i_am_accel_node .AND. acc_on )
#endif

    IF(PRESENT(add)) THEN
#ifdef _OPENACC
!$ACC KERNELS PRESENT( add, send, recv ), IF( i_am_accel_node .AND. acc_on )
      DO i=1,n_pnts
        recv( tmp_recv_dst_idx(i), 1:ndim2, tmp_recv_dst_blk(i) )  =                    &
          &  add( tmp_recv_dst_idx(i), 1:ndim2, tmp_recv_dst_blk(i) )                +  &
          &  send(tmp_send_src_idx(tmp_recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       tmp_send_src_blk(tmp_recv_src(i)))
      END DO
!$ACC END KERNELS
#else
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
          &  add( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )                +  &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
#endif
    ELSE
#ifdef _OPENACC
!$ACC KERNELS PRESENT( send, recv ), IF( i_am_accel_node .AND. acc_on )
      DO i=1,n_pnts
        recv( tmp_recv_dst_idx(i), 1:ndim2, tmp_recv_dst_blk(i) )  =                    &
          &  send(tmp_send_src_idx(tmp_recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       tmp_send_src_blk(tmp_recv_src(i)))
      END DO
!$ACC END KERNELS
#else
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
#endif
    END IF

  END SUBROUTINE exchange_data_r3d_seq


  ! SEQUENTIAL version of subroutine "exchange_data_s3d"
  !
  SUBROUTINE exchange_data_s3d_seq(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), INTENT(IN) :: p_pat
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_s3d_seq"
    INTEGER :: i, ndim2
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif

    ! consistency checks
    ! ------------------

    ! make sure that we are in sequential mode
    IF (.NOT. my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Internal error: sequential routine called in parallel run!")
    END IF
    ! further tests
    IF ( (p_pat%np_recv /= 1) .OR. (p_pat%np_send /= 1) ) THEN
      CALL finish(routine, "Internal error: inconsistent no. send/receive peers!")
    END IF
    IF ( (p_pat%recv_limits(1) - p_pat%recv_limits(0)) /= (p_pat%send_limits(1) - p_pat%send_limits(0)) ) THEN
      CALL finish(routine, "Internal error: inconsistent sender/receiver size!")
    END IF
    IF ( (p_pat%recv_limits(0) /= 0) .OR. (p_pat%send_limits(0) /= 0) ) THEN
      CALL finish(routine, "Internal error: inconsistent sender/receiver start position!")
    END IF
    IF ( (p_pat%recv_limits(1) /= p_pat%n_recv) .OR. (p_pat%n_recv /= p_pat%n_send) ) THEN
      CALL finish(routine, "Internal error: inconsistent counts for sender/receiver!")
    END IF

    ! "communication" (direct copy)
    ! -----------------------------

    ndim2 = SIZE(recv,2)

#ifdef _OPENACC
! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
#endif

    ! The next piece of code is a condensed version of the following
    ! (under the assumptions asserted above):
    !
    !     ! fill sender buffer
    !     DO i=1,n_send
    !       send_buf(i) = array_in(send_src_idx(i), send_src_blk(i))
    !     END DO
    !     ! copy sender to receiver buffer
    !     recv_buf(p_pat%recv_limits(0)+1:p_pat%recv_limits(1)) = &
    !       &  send_buf(p_pat%send_limits(0)+1:p_pat%send_limits(1))
    !     ! copy from receiver buffer
    !     DO i=1,n_pnts
    !       array_out( recv_dst_idx(i), recv_dst_blk(i) ) = recv_buf(recv_src(i))
    !     END DO

    IF(PRESENT(add)) THEN
#ifdef _OPENACC
!$ACC KERNELS PRESENT( add, send, recv ), IF( i_am_accel_node .AND. acc_on )
      DO i=1,n_pnts
        recv( tmp_recv_dst_idx(i), 1:ndim2, tmp_recv_dst_blk(i) )  =                    &
          &  add( tmp_recv_dst_idx(i), 1:ndim2, tmp_recv_dst_blk(i) )                +  &
          &  send(tmp_send_src_idx(tmp_recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       tmp_send_src_blk(tmp_recv_src(i)))
      END DO
!$ACC END KERNELS
#else
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
          &  add( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )                +  &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
#endif
    ELSE
#ifdef _OPENACC
!$ACC KERNELS PRESENT( send, recv ), IF( i_am_accel_node .AND. acc_on )
      DO i=1,n_pnts
        recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
!$ACC END KERNELS
#else
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
#endif
    END IF

  END SUBROUTINE exchange_data_s3d_seq


  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_i3d(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), INTENT(INOUT) :: p_pat
    INTEGER, INTENT(INOUT), TARGET            :: recv(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET     :: send(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET     :: add (:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_i3d"
    INTEGER :: send_buf(SIZE(recv,2),p_pat%n_send), &
      recv_buf(SIZE(recv,2),p_pat%n_recv)

    INTEGER, POINTER :: send_ptr(:,:,:)

    INTEGER :: i, k, np, irs, iss, pid, icount, ndim2
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif

    IF(my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Not yet implemented!")
    END IF

    !-----------------------------------------------------------------------
    start_sync_timer(timer_exch_data)

    IF(my_process_is_mpi_seq()) &
      CALL finish('exchange_data','must not be called on single PE/test PE')

    IF(SIZE(recv,1) /= nproma) THEN
      CALL finish('exchange_data','Illegal first dimension of data array')
    ENDIF

    ndim2 = SIZE(recv,2)

#ifdef _OPENACC
! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
#endif

!$ACC DATA CREATE( send_buf, recv_buf ), &
!$ACC      COPYIN( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src ), &
!$ACC      IF ( i_am_accel_node .AND. acc_on )

    IF (iorder_sendrecv == 1 .OR. iorder_sendrecv >= 3) THEN
      ! Set up irecv's for receive buffers
      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2
        CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF

    ! Set up send buffer

    IF(PRESENT(send)) THEN
      send_ptr => send
    ELSE
      send_ptr => recv
    ENDIF

    IF (ndim2 == 1) THEN
#if _OPENACC
!$ACC PARALLEL PRESENT( send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
      DO i = 1, n_send
        send_buf(1,i) = send_ptr(tmp_send_src_idx(i),1,tmp_send_src_blk(i))
      ENDDO
!$ACC END PARALLEL
#else
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i))
      ENDDO
#endif
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
      DO k = 1, ndim2
!$ACC LOOP VECTOR
        DO i = 1, n_send
          send_buf(k,i) = send_ptr(tmp_send_src_idx(i),k,tmp_send_src_blk(i))
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
      DO i = 1, p_pat%n_send
        send_buf(1:ndim2,i) = send_ptr(p_pat%send_src_idx(i),1:ndim2,   &
          &                            p_pat%send_src_blk(i))
      ENDDO
#endif
    ENDIF

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    ! Send our data
    IF (iorder_sendrecv == 1) THEN
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2
        CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2
        CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO

      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2
        CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ELSE IF (iorder_sendrecv >= 3) THEN ! use irecv/isend
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2
        CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF

    ! Wait for all outstanding requests to finish

    CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    ! Fill in receive buffer

    IF(PRESENT(add)) THEN
      IF (ndim2 == 1) THEN
        k = 1
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( recv_buf, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
        DO i = 1, n_pnts
          recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) =   &
            recv_buf(k,tmp_recv_src(i)) +                       &
            add(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i))
        ENDDO
!$ACC END PARALLEL
#else
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
            recv_buf(k,p_pat%recv_src(i)) +                       &
            add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
        ENDDO
#endif
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( recv_buf, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, n_pnts
            recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) =   &
              recv_buf(k,tmp_recv_src(i)) +                       &
              add(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) =   &
            recv_buf(:,p_pat%recv_src(i)) +                       &
            add(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i))
        ENDDO
#endif
      ENDIF
    ELSE
      IF (ndim2 == 1) THEN
        k = 1
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
        DO i = 1, n_pnts
          recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) = &
            recv_buf(k,tmp_recv_src(i))
        ENDDO
!$ACC END PARALLEL
#else
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
            recv_buf(k,p_pat%recv_src(i))
        ENDDO
#endif
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, n_pnts
            recv(tmp_recv_dst_idx(i),k,tmp_recv_dst_blk(i)) = &
              recv_buf(k,tmp_recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) = &
            recv_buf(:,p_pat%recv_src(i))
        ENDDO
#endif
      ENDIF
    ENDIF

!$ACC END DATA

    stop_sync_timer(timer_exch_data)

  END SUBROUTINE exchange_data_i3d

  !>
  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Optimized version by Guenther Zaengl to process 4D fields or up to seven 3D fields
  !! in one step
  !!
  SUBROUTINE exchange_data_mult(p_pat, ndim2tot, &
       recv, send, nshift)

    CLASS(t_comm_pattern_orig), INTENT(INOUT) :: p_pat
    INTEGER, INTENT(IN)           :: ndim2tot
    INTEGER, OPTIONAL, INTENT(IN) :: nshift

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_mult"
    TYPE(t_ptr_3d), PTR_INTENT(in) :: recv(:)
    TYPE(t_ptr_3d), OPTIONAL, PTR_INTENT(in) :: send(:)
    INTEGER        :: ndim2(SIZE(recv)), noffset(SIZE(recv))

    REAL(dp) :: send_buf(ndim2tot,p_pat%n_send),recv_buf(ndim2tot,p_pat%n_recv)
#if defined( __SX__ ) || defined( _OPENACC )
    REAL(dp), POINTER :: send_ptr(:,:,:), recv_ptr(:,:,:)  ! Refactoring for OpenACC
#endif
    INTEGER :: nfields
    INTEGER :: i, k, kshift(SIZE(recv)), jb,ik, jl, n, np, irs, iss, pid, icount
    LOGICAL :: lsend
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif

    !-----------------------------------------------------------------------

    nfields = SIZE(recv)
    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    lsend     = PRESENT(send)

    IF (PRESENT(nshift)) THEN
      kshift = nshift
    ELSE
      kshift = 0
    ENDIF

#ifdef _OPENACC
! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
#endif

!$ACC DATA CREATE( send_buf, recv_buf ), &
!$ACC      COPYIN( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src ), &
!$ACC      IF ( i_am_accel_node .AND. acc_on )

    IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) .AND. &
      & .NOT. my_process_is_mpi_seq()) THEN
      ! Set up irecv's for receive buffers
      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2tot
        CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF


    IF(my_process_is_mpi_seq()) THEN
      DO n = 1, nfields
        IF(lsend) THEN
          CALL exchange_data_r3d_seq(p_pat, recv(n)%p(:,:,:), send(n)%p(:,:,:))
        ELSE
          CALL exchange_data_r3d_seq(p_pat, recv(n)%p(:,:,:))
        ENDIF
      ENDDO

    ELSE          ! WS: Removed RETURN in order to properly support OpenACC DATA region

      ! Reset kshift to 0 if 2D fields are passed together with 3D fields
      DO n = 1, nfields
        IF (SIZE(recv(n)%p,2) == 1) kshift(n) = 0
      ENDDO

      noffset(1) = 0
      ndim2(1)   = SIZE(recv(1)%p,2) - kshift(1)
      DO n = 2, nfields
        noffset(n) = noffset(n-1)+ndim2(n-1)
        ndim2(n)   = SIZE(recv(n)%p,2) - kshift(n)
      ENDDO

      ! Set up send buffer
#if defined( __SX__ ) || defined( _OPENACC )
      IF ( lsend ) THEN
        DO n = 1, nfields
          send_ptr => send(n)%p   ! Refactoring for OpenACC
!$ACC PARALLEL PRESENT( send_ptr, send_buf ), COPYIN( kshift, noffset, ndim2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
          DO k = 1, ndim2(n)
!$ACC LOOP VECTOR
            DO i = 1, n_send
              send_buf(k+noffset(n),i) = &
                send_ptr(tmp_send_src_idx(i),k+kshift(n),tmp_send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
      ELSE
        ! Send and receive arrays are identical (for boundary exchange)
        DO n = 1, nfields
          recv_ptr => recv(n)%p
!$ACC PARALLEL PRESENT( recv_ptr, send_buf ), COPYIN( kshift, noffset, ndim2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
          DO k = 1, ndim2(n)
!$ACC LOOP VECTOR
            DO i = 1, n_send
              send_buf(k+noffset(n),i) = &
                recv_ptr(tmp_send_src_idx(i),k+kshift(n),tmp_send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
      ENDIF
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,n,k)
#endif
      DO i = 1, p_pat%n_send
        jb = p_pat%send_src_blk(i)
        jl = p_pat%send_src_idx(i)
        IF ( lsend ) THEN
          DO n = 1, nfields
            DO k = 1, ndim2(n)
              send_buf(k+noffset(n),i) = send(n)%p(jl,k+kshift(n),jb)
            ENDDO
          ENDDO
        ELSE
          DO n = 1, nfields
            DO k = 1, ndim2(n)
              send_buf(k+noffset(n),i) = recv(n)%p(jl,k+kshift(n),jb)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

      ! Send our data
      IF (iorder_sendrecv == 1) THEN
        DO np = 1, p_pat%np_send ! loop over PEs where to send the data

          pid    = p_pat%pelist_send(np) ! ID of sender PE
          iss    = p_pat%send_startidx(np)
          icount = p_pat%send_count(np)*ndim2tot
          CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

        ENDDO
      ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
        DO np = 1, p_pat%np_send ! loop over PEs where to send the data

          pid    = p_pat%pelist_send(np) ! ID of sender PE
          iss    = p_pat%send_startidx(np)
          icount = p_pat%send_count(np)*ndim2tot
          CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

        ENDDO

        DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

          pid    = p_pat%pelist_recv(np) ! ID of receiver PE
          irs    = p_pat%recv_startidx(np)
          icount = p_pat%recv_count(np)*ndim2tot
          CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

        ENDDO
      ELSE IF (iorder_sendrecv == 3) THEN ! use isend/irecv
        DO np = 1, p_pat%np_send ! loop over PEs where to send the data

          pid    = p_pat%pelist_send(np) ! ID of sender PE
          iss    = p_pat%send_startidx(np)
          icount = p_pat%send_count(np)*ndim2tot
          CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

        ENDDO
      ENDIF

      IF (iorder_sendrecv > 0) THEN
        ! Wait for all outstanding requests to finish
        start_sync_timer(timer_exch_data_wait)
        CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif
        stop_sync_timer(timer_exch_data_wait)

      ELSE IF (iorder_sendrecv == 0) THEN
        ! dummy mode; just fill the receive buffer with "something reasonable"
        IF (p_pat%n_send > 0) THEN
!$ACC KERNELS PRESENT( send_buf, recv_buf ), IF( i_am_accel_node .AND. acc_on )
          DO k = 1, ndim2tot
            recv_buf(k,:) = send_buf(k,1)
          ENDDO
!$ACC END KERNELS
        ELSE
!$ACC KERNELS PRESENT( recv_buf ), IF( i_am_accel_node .AND. acc_on )
          DO k = 1, ndim2tot
            recv_buf(k,:) = 1._dp
          ENDDO
!$ACC END KERNELS
        ENDIF
      ENDIF

      IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
        start_sync_timer(timer_barrier)
        CALL p_barrier(p_pat%comm)
        stop_sync_timer(timer_barrier)
      ENDIF

      ! Fill in receive buffer

#if defined( __SX__ ) || defined( _OPENACC )
      DO n = 1, nfields
        recv_ptr => recv(n)%p
!$ACC PARALLEL &
!$ACC PRESENT( recv_buf, recv_ptr ), COPYIN( kshift, noffset, ndim2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2(n)
!$ACC LOOP VECTOR
          DO i = 1, n_pnts
            recv_ptr(tmp_recv_dst_idx(i),k+kshift(n),tmp_recv_dst_blk(i)) =  &
              recv_buf(k+noffset(n),tmp_recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,ik,n,k)
#endif
      DO i = 1, p_pat%n_pnts
        jb = p_pat%recv_dst_blk(i)
        jl = p_pat%recv_dst_idx(i)
        ik  = p_pat%recv_src(i)
        DO n = 1, nfields
          DO k = 1, ndim2(n)
            recv(n)%p(jl,k+kshift(n),jb) = recv_buf(k+noffset(n),ik)
          ENDDO
        ENDDO
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif

    ENDIF  ! .NOT. my_process_is_mpi_seq()

!$ACC END DATA

    stop_sync_timer(timer_exch_data)

  END SUBROUTINE exchange_data_mult

  !>
  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Optimized version by Guenther Zaengl to process 4D fields and 3D fields with either single
  !! precision or double precision
  !!
  SUBROUTINE exchange_data_mult_mixprec(p_pat, nfields_dp, ndim2tot_dp, &
       nfields_sp, ndim2tot_sp, recv_dp, send_dp, recv_sp, send_sp, nshift)

    CLASS(t_comm_pattern_orig), INTENT(INOUT) :: p_pat

    INTEGER, INTENT(IN)           :: nfields_dp, ndim2tot_dp, nfields_sp, ndim2tot_sp
    TYPE(t_ptr_3d), PTR_INTENT(in), OPTIONAL :: recv_dp(:)
    TYPE(t_ptr_3d), PTR_INTENT(in), OPTIONAL :: send_dp(:)
    TYPE(t_ptr_3d_sp), PTR_INTENT(in), OPTIONAL :: recv_sp(:)
    TYPE(t_ptr_3d_sp), PTR_INTENT(in), OPTIONAL :: send_sp(:)

    INTEGER, OPTIONAL, INTENT(IN) :: nshift

    INTEGER             :: ndim2_dp(nfields_dp), noffset_dp(nfields_dp), &
                           ndim2_sp(nfields_sp), noffset_sp(nfields_sp)

    REAL(dp) :: send_buf_dp(ndim2tot_dp,p_pat%n_send),recv_buf_dp(ndim2tot_dp,p_pat%n_recv)
    REAL(sp) :: send_buf_sp(ndim2tot_sp,p_pat%n_send),recv_buf_sp(ndim2tot_sp,p_pat%n_recv)
#if defined( __SX__ ) || defined( _OPENACC )
    ! Refactoring for OpenACC
    REAL(sp), POINTER :: send_fld_sp(:,:,:), recv_fld_sp(:,:,:)
    REAL(dp), POINTER :: send_fld_dp(:,:,:), recv_fld_dp(:,:,:)
#endif
    INTEGER :: i, k, kshift_dp(nfields_dp), kshift_sp(nfields_sp), &
         jb, ik, jl, n, np, irs, iss, pid, icount
    LOGICAL :: lsend
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif

    !-----------------------------------------------------------------------

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    lsend     = PRESENT(send_sp)

    IF (PRESENT(nshift)) THEN
      kshift_dp = nshift
      kshift_sp = nshift
    ELSE
      kshift_dp = 0
      kshift_sp = 0
    ENDIF

#ifdef _OPENACC
! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
#endif

!$ACC DATA CREATE( send_buf_dp, recv_buf_dp, send_buf_sp, recv_buf_sp ), &
!$ACC      COPYIN( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src ), &
!$ACC      IF ( i_am_accel_node .AND. acc_on )

    IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) .AND. &
      & .NOT. my_process_is_mpi_seq()) THEN
      ! Set up irecv's for receive buffers
      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2tot_dp
        IF (icount>0) CALL p_irecv(recv_buf_dp(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

        icount = p_pat%recv_count(np)*ndim2tot_sp
        IF (icount>0) CALL p_irecv(recv_buf_sp(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF

    IF(my_process_is_mpi_seq()) THEN
      DO n = 1, nfields_dp
        IF(lsend) THEN
          CALL exchange_data_r3d_seq(p_pat, recv_dp(n)%p(:,:,:), send_dp(n)%p(:,:,:))
        ELSE
          CALL exchange_data_r3d_seq(p_pat, recv_dp(n)%p(:,:,:))
        ENDIF
      ENDDO
      DO n = 1, nfields_sp
        IF(lsend) THEN
          CALL exchange_data_s3d_seq(p_pat, recv_sp(n)%p(:,:,:), send_sp(n)%p(:,:,:))
        ELSE
          CALL exchange_data_s3d_seq(p_pat, recv_sp(n)%p(:,:,:))
        ENDIF
      ENDDO

    ELSE          ! WS: Removed RETURN in order to properly support OpenACC DATA region

      ! Reset kshift to 0 if 2D fields are passed together with 3D fields
      DO n = 1, nfields_dp
        IF (SIZE(recv_dp(n)%p,2) == 1) kshift_dp(n) = 0
      ENDDO
      DO n = 1, nfields_sp
        IF (SIZE(recv_sp(n)%p,2) == 1) kshift_sp(n) = 0
      ENDDO

      IF (nfields_dp > 0) THEN
        noffset_dp(1) = 0
        ndim2_dp(1)   = SIZE(recv_dp(1)%p,2) - kshift_dp(1)
      ENDIF
      DO n = 2, nfields_dp
        noffset_dp(n) = noffset_dp(n-1)+ndim2_dp(n-1)
        ndim2_dp(n)   = SIZE(recv_dp(n)%p,2) - kshift_dp(n)
      ENDDO
      IF (nfields_sp > 0) THEN
        noffset_sp(1) = 0
        ndim2_sp(1)   = SIZE(recv_sp(1)%p,2) - kshift_sp(1)
      ENDIF
      DO n = 2, nfields_sp
        noffset_sp(n) = noffset_sp(n-1)+ndim2_sp(n-1)
        ndim2_sp(n)   = SIZE(recv_sp(n)%p,2) - kshift_sp(n)
      ENDDO


      ! Set up send buffer
#if defined( __SX__ ) || defined( _OPENACC )
      IF ( lsend ) THEN
        DO n = 1, nfields_dp
          send_fld_dp => send_dp(n)%p   ! Refactoring for OpenACC
!$ACC PARALLEL PRESENT(send_fld_dp,send_buf_dp), COPYIN(kshift_dp,noffset_dp,ndim2_dp), IF(i_am_accel_node .AND. acc_on)
!$ACC LOOP GANG
          DO k = 1, ndim2_dp(n)
!$ACC LOOP VECTOR
            DO i = 1, n_send
              send_buf_dp(k+noffset_dp(n),i) = &
                send_fld_dp(tmp_send_src_idx(i),k+kshift_dp(n),tmp_send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
        DO n = 1, nfields_sp
          send_fld_sp => send_sp(n)%p   ! Refactoring for OpenACC
!$ACC PARALLEL PRESENT(send_fld_sp,send_buf_sp), COPYIN(kshift_sp,noffset_sp,ndim2_sp), IF(i_am_accel_node .AND. acc_on)
!$ACC LOOP GANG
          DO k = 1, ndim2_sp(n)
!$ACC LOOP VECTOR
            DO i = 1, n_send
              send_buf_sp(k+noffset_sp(n),i) = &
                send_fld_sp(tmp_send_src_idx(i),k+kshift_sp(n),tmp_send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
      ELSE
        ! Send and receive arrays are identical (for boundary exchange)
        DO n = 1, nfields_dp
          recv_fld_dp => recv_dp(n)%p
!$ACC PARALLEL PRESENT(recv_fld_dp,send_buf_dp), COPYIN(kshift_dp,noffset_dp,ndim2_dp), IF(i_am_accel_node .AND. acc_on)
!$ACC LOOP GANG
!CDIR UNROLL=6
          DO k = 1, ndim2_dp(n)
!$ACC LOOP VECTOR
            DO i = 1, n_send
              send_buf_dp(k+noffset_dp(n),i) = &
                recv_fld_dp(tmp_send_src_idx(i),k+kshift_dp(n),tmp_send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
        DO n = 1, nfields_sp
          recv_fld_sp => recv_sp(n)%p
!$ACC PARALLEL PRESENT(recv_fld_sp,send_buf_sp), COPYIN(kshift_sp,noffset_sp,ndim2_sp), IF(i_am_accel_node .AND. acc_on)
!$ACC LOOP GANG
!CDIR UNROLL=6
          DO k = 1, ndim2_sp(n)
!$ACC LOOP VECTOR
            DO i = 1, n_send
              send_buf_sp(k+noffset_sp(n),i) = &
                recv_fld_sp(tmp_send_src_idx(i),k+kshift_sp(n),tmp_send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
      ENDIF
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,n,k)
#endif
      DO i = 1, p_pat%n_send
        jb = p_pat%send_src_blk(i)
        jl = p_pat%send_src_idx(i)
        IF ( lsend ) THEN
          DO n = 1, nfields_dp
            DO k = 1, ndim2_dp(n)
              send_buf_dp(k+noffset_dp(n),i) = send_dp(n)%p(jl,k+kshift_dp(n),jb)
            ENDDO
          ENDDO
          DO n = 1, nfields_sp
            DO k = 1, ndim2_sp(n)
              send_buf_sp(k+noffset_sp(n),i) = send_sp(n)%p(jl,k+kshift_sp(n),jb)
            ENDDO
          ENDDO
        ELSE
          DO n = 1, nfields_dp
            DO k = 1, ndim2_dp(n)
              send_buf_dp(k+noffset_dp(n),i) = recv_dp(n)%p(jl,k+kshift_dp(n),jb)
            ENDDO
          ENDDO
          DO n = 1, nfields_sp
            DO k = 1, ndim2_sp(n)
              send_buf_sp(k+noffset_sp(n),i) = recv_sp(n)%p(jl,k+kshift_sp(n),jb)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf_dp, send_buf_sp ), IF ( i_am_accel_node .AND. acc_on )
#endif

      ! Send our data
      IF (iorder_sendrecv == 1) THEN
        DO np = 1, p_pat%np_send ! loop over PEs where to send the data

          pid    = p_pat%pelist_send(np) ! ID of sender PE
          iss    = p_pat%send_startidx(np)
          icount = p_pat%send_count(np)*ndim2tot_dp
          IF (icount>0) CALL p_send(send_buf_dp(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)
          icount = p_pat%send_count(np)*ndim2tot_sp
          IF (icount>0) CALL p_send(send_buf_sp(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

        ENDDO
      ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
        DO np = 1, p_pat%np_send ! loop over PEs where to send the data

          pid    = p_pat%pelist_send(np) ! ID of sender PE
          iss    = p_pat%send_startidx(np)
          icount = p_pat%send_count(np)*ndim2tot_dp
          IF (icount>0) CALL p_isend(send_buf_dp(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)
          icount = p_pat%send_count(np)*ndim2tot_sp
          IF (icount>0) CALL p_isend(send_buf_sp(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

        ENDDO

        DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

          pid    = p_pat%pelist_recv(np) ! ID of receiver PE
          irs    = p_pat%recv_startidx(np)
          icount = p_pat%recv_count(np)*ndim2tot_dp
          IF (icount>0) CALL p_recv(recv_buf_dp(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)
          icount = p_pat%recv_count(np)*ndim2tot_sp
          IF (icount>0) CALL p_recv(recv_buf_sp(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

        ENDDO
      ELSE IF (iorder_sendrecv == 3) THEN ! use isend/irecv
        DO np = 1, p_pat%np_send ! loop over PEs where to send the data

          pid    = p_pat%pelist_send(np) ! ID of sender PE
          iss    = p_pat%send_startidx(np)
          icount = p_pat%send_count(np)*ndim2tot_dp
          IF (icount>0) CALL p_isend(send_buf_dp(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)
          icount = p_pat%send_count(np)*ndim2tot_sp
          IF (icount>0) CALL p_isend(send_buf_sp(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

        ENDDO
      ENDIF

      IF (iorder_sendrecv > 0) THEN
        ! Wait for all outstanding requests to finish
        start_sync_timer(timer_exch_data_wait)
        CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf_dp, recv_buf_sp ), IF ( i_am_accel_node .AND. acc_on )
#endif
        stop_sync_timer(timer_exch_data_wait)

      ELSE IF (iorder_sendrecv == 0) THEN
        ! dummy mode; just fill the receive buffer with "something reasonable"
        IF (p_pat%n_send > 0) THEN
!$ACC KERNELS PRESENT( send_buf_dp, recv_buf_dp, send_buf_sp, recv_buf_sp ), IF( i_am_accel_node .AND. acc_on )
          DO k = 1, ndim2tot_dp
            recv_buf_dp(k,:) = send_buf_dp(k,1)
          ENDDO
          DO k = 1, ndim2tot_sp
            recv_buf_sp(k,:) = send_buf_sp(k,1)
          ENDDO
!$ACC END KERNELS
        ELSE
!$ACC KERNELS PRESENT( recv_buf_dp, recv_buf_sp ), IF( i_am_accel_node .AND. acc_on )
          DO k = 1, ndim2tot_dp
            recv_buf_dp(k,:) = 1._dp
          ENDDO
          DO k = 1, ndim2tot_sp
            recv_buf_sp(k,:) = 1._sp
          ENDDO
!$ACC END KERNELS
        ENDIF
      ENDIF

      IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
        start_sync_timer(timer_barrier)
        CALL p_barrier(p_pat%comm)
        stop_sync_timer(timer_barrier)
      ENDIF

      ! Fill in receive buffer

#if defined( __SX__ ) || defined( _OPENACC )
      DO n = 1, nfields_dp
        recv_fld_dp => recv_dp(n)%p
!$ACC PARALLEL &
!$ACC PRESENT( recv_buf_dp, recv_fld_dp ), COPYIN( kshift_dp, noffset_dp, ndim2_dp ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2_dp(n)
!$ACC LOOP VECTOR
          DO i = 1, n_pnts
            recv_fld_dp(tmp_recv_dst_idx(i),k+kshift_dp(n),tmp_recv_dst_blk(i)) =  &
              recv_buf_dp(k+noffset_dp(n),tmp_recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
      DO n = 1, nfields_sp
        recv_fld_sp => recv_sp(n)%p
!$ACC PARALLEL &
!$ACC PRESENT( recv_buf_sp, recv_fld_sp ), COPYIN( kshift_sp, noffset_sp, ndim2_sp ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2_sp(n)
!$ACC LOOP VECTOR
          DO i = 1, n_pnts
            recv_fld_sp(tmp_recv_dst_idx(i),k+kshift_sp(n),tmp_recv_dst_blk(i)) =  &
              recv_buf_sp(k+noffset_sp(n),tmp_recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,ik,n,k)
#endif
      DO i = 1, p_pat%n_pnts
        jb = p_pat%recv_dst_blk(i)
        jl = p_pat%recv_dst_idx(i)
        ik  = p_pat%recv_src(i)
        DO n = 1, nfields_dp
          DO k = 1, ndim2_dp(n)
            recv_dp(n)%p(jl,k+kshift_dp(n),jb) = recv_buf_dp(k+noffset_dp(n),ik)
          ENDDO
        ENDDO
        DO n = 1, nfields_sp
          DO k = 1, ndim2_sp(n)
            recv_sp(n)%p(jl,k+kshift_sp(n),jb) = recv_buf_sp(k+noffset_sp(n),ik)
          ENDDO
        ENDDO
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif

    ENDIF  ! .NOT. my_process_is_mpi_seq()

!$ACC END DATA

    stop_sync_timer(timer_exch_data)

  END SUBROUTINE exchange_data_mult_mixprec


  !>
  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Optimized version by Guenther Zaengl to process a 4D field whose extra dimension
  !! is on the first index
  !!
  SUBROUTINE exchange_data_4de1(p_pat, nfields, ndim2tot, recv, send)

    CLASS(t_comm_pattern_orig), INTENT(INOUT) :: p_pat

    REAL(dp), INTENT(INOUT)           :: recv(:,:,:,:)
    REAL(dp), INTENT(IN   ), OPTIONAL :: send(:,:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_4de1"
    INTEGER, INTENT(IN)           :: nfields, ndim2tot

    INTEGER :: ndim2, koffset

    REAL(dp) :: send_buf(ndim2tot,p_pat%n_send),recv_buf(ndim2tot,p_pat%n_recv)

    INTEGER :: i, k, ik, jb, jl, n, np, irs, iss, pid, icount
    LOGICAL :: lsend
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif

    !-----------------------------------------------------------------------

    IF(my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Not yet implemented!")
    END IF

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    IF (PRESENT(send)) THEN
      lsend  = .TRUE.
    ELSE
      lsend  = .FALSE.
    ENDIF

    ndim2 = SIZE(recv,3)

#ifdef _OPENACC
! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
#endif

!$ACC DATA CREATE( send_buf, recv_buf ), &
!$ACC      COPYIN( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src ), &
!$ACC      IF ( i_am_accel_node .AND. acc_on )

    IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3)) THEN
      ! Set up irecv's for receive buffers
      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2tot
        CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( send, recv, send_buf ), &
!$ACC IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
    DO i = 1, n_send
      jb = tmp_send_src_blk(i)
      jl = tmp_send_src_idx(i)
      IF ( lsend ) THEN
!$ACC LOOP VECTOR
        DO k = 1, ndim2
          koffset = (k-1)*nfields
          DO n = 1, nfields
            send_buf(n+koffset,i) = send(n,jl,k,jb)
          ENDDO
        ENDDO
      ELSE
!$ACC LOOP VECTOR
        DO k = 1, ndim2
          koffset = (k-1)*nfields
          DO n = 1, nfields
            send_buf(n+koffset,i) = recv(n,jl,k,jb)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__ 
!$OMP PARALLEL DO PRIVATE(jb,jl,koffset,k,n)
#endif
    DO i = 1, p_pat%n_send
      jb = p_pat%send_src_blk(i)
      jl = p_pat%send_src_idx(i)
      IF ( lsend ) THEN
        DO k = 1, ndim2
          koffset = (k-1)*nfields
          DO n = 1, nfields
            send_buf(n+koffset,i) = send(n,jl,k,jb)
          ENDDO
        ENDDO
      ELSE
        DO k = 1, ndim2
          koffset = (k-1)*nfields
          DO n = 1, nfields
            send_buf(n+koffset,i) = recv(n,jl,k,jb)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif

#endif

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    ! Send our data
    IF (iorder_sendrecv == 1) THEN
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2tot
        CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2tot
        CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO

      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2tot
        CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ELSE IF (iorder_sendrecv == 3) THEN ! use isend/irecv
      DO np = 1, p_pat%np_send ! loop over PEs where to send the data

        pid    = p_pat%pelist_send(np) ! ID of sender PE
        iss    = p_pat%send_startidx(np)
        icount = p_pat%send_count(np)*ndim2tot
        CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF

    IF (iorder_sendrecv > 0) THEN
      ! Wait for all outstanding requests to finish
      start_sync_timer(timer_exch_data_wait)
      CALL p_wait
      stop_sync_timer(timer_exch_data_wait)

#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    ELSE IF (iorder_sendrecv == 0) THEN
      ! dummy mode; just fill the receive buffer with "something reasonable"
      IF (p_pat%n_send > 0) THEN
!$ACC KERNELS PRESENT( recv_buf, send_buf ), IF ( i_am_accel_node .AND. acc_on )
        DO k = 1, ndim2tot
          recv_buf(k,:) = send_buf(k,1)
        ENDDO
!$ACC END KERNELS
      ELSE
!$ACC KERNELS PRESENT( recv_buf ), IF ( i_am_accel_node .AND. acc_on )
        DO k = 1, ndim2tot
          recv_buf(k,:) = 1._dp
        ENDDO
!$ACC END KERNELS
      ENDIF
    ENDIF

    IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer
#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( recv_buf, recv ), &
!$ACC IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
    DO i = 1, n_pnts
      jb = tmp_recv_dst_blk(i)
      jl = tmp_recv_dst_idx(i)
      ik  = tmp_recv_src(i)
!$ACC LOOP VECTOR
      DO k = 1, ndim2
        koffset = (k-1)*nfields
        DO n = 1, nfields
          recv(n,jl,k,jb) = recv_buf(n+koffset,ik)
        ENDDO
      ENDDO
    ENDDO
!$ACC END PARALLEL
#else
#if defined( __OMPPAR_COPY__ )
!$OMP PARALLEL DO PRIVATE(jb,jl,ik,koffset,k,n)
#endif
    DO i = 1, p_pat%n_pnts
      jb = p_pat%recv_dst_blk(i)
      jl = p_pat%recv_dst_idx(i)
      ik  = p_pat%recv_src(i)
      DO k = 1, ndim2
        koffset = (k-1)*nfields
        DO n = 1, nfields
          recv(n,jl,k,jb) = recv_buf(n+koffset,ik)
        ENDDO
      ENDDO
    ENDDO
#if defined( __OMPPAR_COPY__ )
!$OMP END PARALLEL DO
#endif

#endif

!$ACC END DATA

    stop_sync_timer(timer_exch_data)

  END SUBROUTINE exchange_data_4de1


  !>
  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Optimized version by Guenther Zaengl to process up to two 4D fields or up to six 3D fields
  !! for an array-sized communication pattern (as needed for boundary interpolation) in one step
  !!
  SUBROUTINE exchange_data_grf(p_pat_coll, nfields, ndim2tot, recv1, send1, &
    recv2, send2, recv3, send3, recv4, send4, &
    recv5, send5, recv6, send6, recv4d1, send4d1, &
    recv4d2, send4d2)

    CLASS(t_comm_pattern_collection_orig), INTENT(INOUT), TARGET :: p_pat_coll

    REAL(dp), INTENT(INOUT), TARGET, OPTIONAL ::  &
      recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4d1(:,:,:,:), &
      recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), recv4d2(:,:,:,:)
    ! Note: the last index of the send fields corresponds to the dimension of p_pat
    ! On the other hand, they are not blocked and have the vertical index first
    REAL(dp), INTENT(IN   ), TARGET, OPTIONAL ::  &
      send1(:,:,:), send2(:,:,:), send3(:,:,:), send4d1(:,:,:,:), &
      send4(:,:,:), send5(:,:,:), send6(:,:,:), send4d2(:,:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_grf"
    INTEGER, INTENT(IN)           :: nfields  ! total number of input fields
    INTEGER, INTENT(IN)           :: ndim2tot ! sum of vertical levels of input fields

    INTEGER           :: nsendtot ! total number of send points
    INTEGER           :: nrecvtot ! total number of receive points

    TYPE t_fieldptr_recv
      REAL(dp), POINTER :: fld(:,:,:)
    END TYPE t_fieldptr_recv

    TYPE t_fieldptr_send
      REAL(dp), POINTER :: fld(:,:)
    END TYPE t_fieldptr_send

    TYPE(t_fieldptr_recv) :: recv(nfields)
    TYPE(t_fieldptr_send) :: send(nfields*SIZE(p_pat_coll%patterns))

    INTEGER        :: ndim2(nfields), noffset(nfields), &
      ioffset_s(SIZE(p_pat_coll%patterns)), &
      ioffset_r(SIZE(p_pat_coll%patterns))

    REAL(dp), ALLOCATABLE :: send_buf(:,:),recv_buf(:,:), &
      auxs_buf(:,:),auxr_buf(:,:)

    INTEGER :: i, j, k, ik, jb, jl, n, np, irs, ire, iss, ise, &
      npats, isum, ioffset, isum1, n4d, pid, num_send, num_recv, &
      comm_size
    INTEGER, ALLOCATABLE :: pelist_send(:), pelist_recv(:)

#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER, ALLOCATABLE   :: tmp_send_src_idx(:), tmp_send_src_blk(:)
    INTEGER, ALLOCATABLE   :: tmp_recv_dst_idx(:), tmp_recv_dst_blk(:)
    INTEGER, ALLOCATABLE   :: tmp_recv_src(:)
    INTEGER   :: n_pnts, n_send
#endif

    TYPE(t_p_comm_pattern_orig), POINTER :: p_pat(:)

    !-----------------------------------------------------------------------

    p_pat => p_pat_coll%patterns

    !-----------------------------------------------------------------------

    nsendtot = 0
    nrecvtot = 0
    DO i = 1, SIZE(p_pat)
      nsendtot = nsendtot + p_pat(i)%p%n_send
      nrecvtot = nrecvtot + p_pat(i)%p%n_recv
    END DO
    comm_size = p_comm_size(p_pat_coll%patterns(1)%p%comm)

    ALLOCATE(send_buf(ndim2tot,nsendtot),recv_buf(ndim2tot,nrecvtot), &
      auxs_buf(ndim2tot,nsendtot),auxr_buf(ndim2tot,nrecvtot))

    !-----------------------------------------------------------------------

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat_coll%patterns(1)%p%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    npats = SIZE(p_pat)  ! Number of communication patterns provided on input

    !-----------------------------------------------------------------------

    ! some adjustmens to the standard communication patterns in order to make
    ! them work in this routine

    num_send = 0
    num_recv = 0

    DO np = 0, comm_size-1 ! loop over PEs

      DO n = 1, npats  ! loop over communication patterns
        iss = p_pat(n)%p%send_limits(np)+1
        ise = p_pat(n)%p%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, npats  ! loop over communication patterns
        irs = p_pat(n)%p%recv_limits(np)+1
        ire = p_pat(n)%p%recv_limits(np+1)
        IF(ire >= irs) THEN
          num_recv = num_recv + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

    ENDDO

    ALLOCATE(pelist_send(num_send), pelist_recv(num_recv))

    num_send = 0
    num_recv = 0

    ! Now compute "envelope PE lists" for all communication patterns
    DO np = 0, comm_size-1 ! loop over PEs

      DO n = 1, npats  ! loop over communication patterns
        iss = p_pat(n)%p%send_limits(np)+1
        ise = p_pat(n)%p%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          DO j = 1, npats
            pelist_send(num_send) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, npats  ! loop over communication patterns
        irs = p_pat(n)%p%recv_limits(np)+1
        ire = p_pat(n)%p%recv_limits(np+1)
        IF(ire >= irs) THEN
          num_recv = num_recv + 1
          DO j = 1, npats
            pelist_recv(num_recv) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

    ENDDO

    !-----------------------------------------------------------------------

!$ACC DATA CREATE( send_buf, recv_buf, auxs_buf, auxr_buf), &
!$ACC      IF ( i_am_accel_node .AND. acc_on )

    ! Set up irecv's for receive buffers
    ! Note: the dummy mode (iorder_sendrecv=0) does not work for nest boundary communication
    ! because there may be PEs receiving but not sending data. Therefore, iorder_sendrecv=0
    ! is treated as iorder_sendrecv=1 in this routine.
    IF ((iorder_sendrecv <= 1 .OR. iorder_sendrecv >= 3) .AND. .NOT. my_process_is_mpi_seq()) THEN

      ioffset = 0
      DO np = 1, num_recv ! loop over PEs from where to receive the data

        pid = pelist_recv(np) ! ID of receiver PE

        ! Sum up receive points over all communication patterns to be processed
        isum = ioffset
        DO n = 1, npats
          isum = isum + p_pat(n)%p%recv_limits(pid+1) - &
            p_pat(n)%p%recv_limits(pid)
        ENDDO

        IF(isum > ioffset) &
          CALL p_irecv(auxr_buf(1,ioffset+1), pid, 1, &
          &            p_count=(isum-ioffset)*ndim2tot, &
          &            comm=p_pat_coll%patterns(1)%p%comm)
        ioffset = isum

      ENDDO

    ENDIF

    ! Set pointers to input fields
    IF (PRESENT(recv4d1) .AND. .NOT. PRESENT(recv4d2)) THEN
      DO n = 1, nfields
        recv(n)%fld => recv4d1(:,:,:,n)
        DO np = 1, npats
          send(np+(n-1)*npats)%fld => send4d1(:,:,np,n)
        ENDDO
      ENDDO
    ELSE IF (PRESENT(recv4d1) .AND. PRESENT(recv4d2)) THEN
      n4d = nfields/2
      DO n = 1, n4d
        recv(n)%fld => recv4d1(:,:,:,n)
        DO np = 1, npats
          send(np+(n-1)*npats)%fld => send4d1(:,:,np,n)
        ENDDO
      ENDDO
      DO n = 1, n4d
        recv(n4d+n)%fld => recv4d2(:,:,:,n)
        DO np = 1, npats
          send(np+(n4d+n-1)*npats)%fld => send4d2(:,:,np,n)
        ENDDO
      ENDDO
    ELSE
      IF (PRESENT(recv1)) THEN
        recv(1)%fld => recv1
        DO np = 1, npats
          send(np)%fld => send1(:,:,np)
        ENDDO
      ENDIF
      IF (PRESENT(recv2)) THEN
        recv(2)%fld => recv2
        DO np = 1, npats
          send(np+npats)%fld => send2(:,:,np)
        ENDDO
      ENDIF
      IF (PRESENT(recv3)) THEN
        recv(3)%fld => recv3
        DO np = 1, npats
          send(np+2*npats)%fld => send3(:,:,np)
        ENDDO
      ENDIF
      IF (PRESENT(recv4)) THEN
        recv(4)%fld => recv4
        DO np = 1, npats
          send(np+3*npats)%fld => send4(:,:,np)
        ENDDO
      ENDIF
      IF (PRESENT(recv5)) THEN
        recv(5)%fld => recv5
        DO np = 1, npats
          send(np+4*npats)%fld => send5(:,:,np)
        ENDDO
      ENDIF
      IF (PRESENT(recv6)) THEN
        recv(6)%fld => recv6
        DO np = 1, npats
          send(np+5*npats)%fld => send6(:,:,np)
        ENDDO
      ENDIF
    ENDIF

    noffset(1) = 0
    ndim2(1)   = SIZE(recv(1)%fld,2)
    DO n = 2, nfields
      noffset(n) = noffset(n-1)+ndim2(n-1)
      ndim2(n)   = SIZE(recv(n)%fld,2)
    ENDDO

    ioffset_r(1) = 0
    ioffset_s(1) = 0
    DO np = 2, npats
      ioffset_r(np) = ioffset_r(np-1) + p_pat(np-1)%p%n_recv
      ioffset_s(np) = ioffset_s(np-1) + p_pat(np-1)%p%n_send
    ENDDO


    IF (my_process_is_mpi_seq()) THEN

#ifdef _OPENACC
      DO np = 1, npats
! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
        n_pnts           = p_pat(np)%p%n_pnts
        n_send           = p_pat(np)%p%n_send
        ALLOCATE( tmp_send_src_idx(n_send), tmp_send_src_blk(n_send), &
                  tmp_recv_dst_idx(n_pnts), tmp_recv_dst_blk(n_pnts), tmp_recv_src(n_pnts) )
        tmp_send_src_idx = p_pat(np)%p%send_src_idx
        tmp_send_src_blk = p_pat(np)%p%send_src_blk
        tmp_recv_dst_idx = p_pat(np)%p%recv_dst_idx
        tmp_recv_dst_blk = p_pat(np)%p%recv_dst_blk
        tmp_recv_src     = p_pat(np)%p%recv_src
!$ACC DATA COPYIN( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src )

!$ACC PARALLEL &
!$ACC PRESENT( send, recv ), COPYIN( ndim2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
        DO i = 1, n_pnts
          DO n = 1, nfields
            DO k = 1, ndim2(n)
              recv(n)%fld( tmp_recv_dst_idx(i), k, &
                tmp_recv_dst_blk(i) ) =            &
                send(np+(n-1)*npats)%fld( k,                             &
                  idx_1d(tmp_send_src_idx(         &
                           tmp_recv_src(i)),       &
                         tmp_send_src_blk(           &
                           tmp_recv_src(i))))
            ENDDO
          ENDDO
        ENDDO
!$ACC END PARALLEL

!$ACC END DATA
        DEALLOCATE( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src )
      ENDDO
#else
      DO np = 1, npats
        DO i = 1, p_pat(np)%p%n_pnts
          DO n = 1, nfields
            DO k = 1, ndim2(n)
              recv(n)%fld( p_pat(np)%p%recv_dst_idx(i), k, &
                p_pat(np)%p%recv_dst_blk(i) ) =            &
                send(np+(n-1)*npats)%fld( k,                             &
                  idx_1d(p_pat(np)%p%send_src_idx(         &
                           p_pat(np)%p%recv_src(i)),       &
                         p_pat(np)%p%send_src_blk(           &
                           p_pat(np)%p%recv_src(i))))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
#endif

    ELSE


      ! Set up send buffer
#if defined( __SX__ ) || defined( _OPENACC )
      DO np = 1, npats

! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
        n_send           = p_pat(np)%p%n_send
        tmp_send_src_idx = p_pat(np)%p%send_src_idx
        tmp_send_src_blk = p_pat(np)%p%send_src_blk
        ALLOCATE( tmp_send_src_idx(n_send), tmp_send_src_blk(n_send) )
!$ACC DATA COPYIN( tmp_send_src_idx, tmp_send_src_blk )

!$ACC PARALLEL &
!$ACC PRESENT( p_pat, send, send_buf ), COPYIN( ndim2, noffset ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP SEQ
        DO n = 1, nfields
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=6
          DO k = 1, ndim2(n)
            DO i = 1, n_send
              send_buf(k+noffset(n),i+ioffset_s(np)) =                &
                & send(np+(n-1)*npats)%fld(k, &
                &   idx_1d(tmp_send_src_idx(i), &
                &          tmp_send_src_blk(i)))
            ENDDO
          ENDDO
        ENDDO
!$ACC END PARALLEL
!$ACC END DATA 
        DEALLOCATE( tmp_send_src_idx, tmp_send_src_blk )
      ENDDO
#else

#ifdef __OMPPAR_COPY__
!$OMP PARALLEL
#endif
      DO np = 1, npats
#ifdef __OMPPAR_COPY__
!$OMP DO PRIVATE(jl,n,k)
#endif
        DO i = 1, p_pat(np)%p%n_send
          jl = idx_1d(p_pat(np)%p%send_src_idx(i), &
            p_pat(np)%p%send_src_blk(i))
          DO n = 1, nfields
            DO k = 1, ndim2(n)
              send_buf(k+noffset(n),i+ioffset_s(np)) = &
                send(np+(n-1)*npats)%fld(k,jl)
            ENDDO
          ENDDO
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END DO
#endif
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL
#endif

#endif

      IF (iorder_sendrecv <= 1) THEN

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

        ! Send our data
        ioffset = 0
        DO np = 1, num_send ! loop over PEs where to send the data

          pid = pelist_send(np) ! ID of sender PE

          ! Copy send points for all communication patterns into one common send buffer
          isum = ioffset
          DO n = 1, npats
            iss = p_pat(n)%p%send_limits(pid)+1 + ioffset_s(n)
            ise = p_pat(n)%p%send_limits(pid+1) + ioffset_s(n)
            isum1 = ise - iss + 1
            IF (isum1 > 0) THEN
!CDIR COLLAPSE
!$ACC KERNELS PRESENT( auxs_buf, send_buf ), IF ( i_am_accel_node .AND. acc_on )
              auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
!$ACC END KERNELS
              isum = isum+isum1
            ENDIF
          ENDDO

#ifndef __USE_G2G
!$ACC UPDATE HOST( auxs_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

          IF(isum > ioffset) CALL p_send(auxs_buf(1,ioffset+1), pid, 1,             &
            p_count=(isum-ioffset)*ndim2tot, comm=p_pat_coll%patterns(1)%p%comm)

          ioffset = isum

        ENDDO
      ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
        ioffset = 0
        DO np = 1, num_send ! loop over PEs where to send the data

          pid = pelist_send(np) ! ID of sender PE

          ! Copy send points for all communication patterns into one common send buffer
          isum = ioffset
          DO n = 1, npats
            iss = p_pat(n)%p%send_limits(pid)+1 + ioffset_s(n)
            ise = p_pat(n)%p%send_limits(pid+1) + ioffset_s(n)
            isum1 = ise - iss + 1
            IF (isum1 > 0) THEN
!CDIR COLLAPSE
!$ACC KERNELS PRESENT( auxs_buf, send_buf ), IF ( i_am_accel_node .AND. acc_on )
              auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
!$ACC END KERNELS
              isum = isum+isum1
            ENDIF
          ENDDO

#ifndef __USE_G2G
!$ACC UPDATE HOST( auxs_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif
          IF(isum > ioffset) CALL p_isend(auxs_buf(1,ioffset+1), pid, 1,            &
            p_count=(isum-ioffset)*ndim2tot, comm=p_pat_coll%patterns(1)%p%comm)

          ioffset = isum

        ENDDO

        ioffset = 0
        DO np = 1, num_recv ! loop over PEs from where to receive the data

          pid = pelist_recv(np) ! ID of receiver PE

          ! Sum up receive points over all communication patterns to be processed
          isum = ioffset
          DO n = 1, npats
            isum = isum + p_pat(n)%p%recv_limits(pid+1) - &
              p_pat(n)%p%recv_limits(pid)
          ENDDO

          IF(isum > ioffset) CALL p_recv(auxr_buf(1,ioffset+1), pid, 1,             &
            p_count=(isum-ioffset)*ndim2tot, comm=p_pat_coll%patterns(1)%p%comm)
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( auxs_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif
          ioffset = isum

        ENDDO
      ELSE IF (iorder_sendrecv >= 3) THEN ! use isend/recv
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL PRIVATE(ioffset,pid,isum,iss,ise,isum1,np,n,i)
#endif
        ioffset = 0
        DO np = 1, num_send ! loop over PEs where to send the data

          pid = pelist_send(np) ! ID of sender PE

          ! Copy send points for all communication patterns into one common send buffer
          isum = ioffset
          DO n = 1, npats
            iss = p_pat(n)%p%send_limits(pid)+1 + ioffset_s(n)
            ise = p_pat(n)%p%send_limits(pid+1) + ioffset_s(n)
            isum1 = ise - iss + 1
            IF (isum1 > 0) THEN
#ifdef __OMPPAR_COPY__
!$OMP DO
#endif
!$ACC KERNELS PRESENT( auxs_buf, send_buf ), IF ( i_am_accel_node .AND. acc_on )
              DO i = 1, isum1
                auxs_buf(:,isum+i) = send_buf(:,iss-1+i)
              ENDDO
!$ACC END KERNELS
#ifdef __OMPPAR_COPY__
!$OMP END DO
#endif
              isum = isum+isum1
            ENDIF
          ENDDO

#ifndef __USE_G2G
!$ACC UPDATE HOST( auxs_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif
!$OMP MASTER
          IF(isum > ioffset) CALL p_isend(auxs_buf(1,ioffset+1), pid, 1,            &
            p_count=(isum-ioffset)*ndim2tot, comm=p_pat_coll%patterns(1)%p%comm)
!$OMP END MASTER

          ioffset = isum

        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL
#endif
      ENDIF

      ! Wait for all outstanding requests to finish
      start_sync_timer(timer_exch_data_wait)
      CALL p_wait
      stop_sync_timer(timer_exch_data_wait)


      IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
        start_sync_timer(timer_barrier)
        CALL p_barrier(p_pat_coll%patterns(1)%p%comm)
        stop_sync_timer(timer_barrier)
      ENDIF

      ! Copy exchanged data back to receive buffer

#ifdef __OMPPAR_COPY__
!$OMP PARALLEL PRIVATE(ioffset,pid,isum,irs,ire,isum1,n,i)
#endif

      ioffset = 0
      DO np = 1, num_recv ! loop over PEs from where to receive the data

        pid = pelist_recv(np) ! ID of receiver PE

        isum = ioffset
        DO n = 1, npats
          irs = p_pat(n)%p%recv_limits(pid)+1 + ioffset_r(n)
          ire = p_pat(n)%p%recv_limits(pid+1) + ioffset_r(n)
          isum1 = ire - irs + 1
          IF (isum1 > 0) THEN
#ifdef __OMPPAR_COPY__
!$OMP DO
#endif
!$ACC KERNELS PRESENT( recv_buf, auxr_buf ), IF ( i_am_accel_node .AND. acc_on )
            DO i = 1, isum1
              recv_buf(:,irs-1+i) = auxr_buf(:,isum+i)
            ENDDO
!$ACC END KERNELS
#ifdef __OMPPAR_COPY__
!$OMP END DO
#endif
            isum = isum + isum1
          ENDIF
        ENDDO

        ioffset = isum

      ENDDO

      ! Fill in receive buffer

#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

#if defined( __SX__ ) || defined( _OPENACC )
      DO np = 1, npats

! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
        n_pnts           = p_pat(np)%p%n_pnts
        ALLOCATE( tmp_recv_dst_idx(n_pnts), tmp_recv_dst_blk(n_pnts), tmp_recv_src(n_pnts) )
        tmp_recv_dst_idx = p_pat(np)%p%recv_dst_idx
        tmp_recv_dst_blk = p_pat(np)%p%recv_dst_blk
        tmp_recv_src     = p_pat(np)%p%recv_src
!$ACC DATA COPYIN( tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src )

!$ACC PARALLEL PRESENT( recv_buf, recv ), COPYIN( ndim2, noffset, ioffset_r ), IF( i_am_accel_node .AND. acc_on )
        DO n = 1, nfields
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=6
          DO k = 1, ndim2(n)
            DO i = 1, n_pnts
              recv(n)%fld(tmp_recv_dst_idx(i),k, &
                          tmp_recv_dst_blk(i)) =   &
                recv_buf(k+noffset(n),tmp_recv_src(i)+ioffset_r(np))
            ENDDO
          ENDDO
        ENDDO
!$ACC END PARALLEL
!$ACC END DATA
        DEALLOCATE( tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src )
      ENDDO
#else
      DO np = 1, npats
#ifdef __OMPPAR_COPY__
!$OMP DO PRIVATE(jb,jl,ik,n,k)
#endif
        DO i = 1, p_pat(np)%p%n_pnts
          jb = p_pat(np)%p%recv_dst_blk(i)
          jl = p_pat(np)%p%recv_dst_idx(i)
          ik  = p_pat(np)%p%recv_src(i)+ioffset_r(np)
          DO n = 1, nfields
            DO k = 1, ndim2(n)
              recv(n)%fld(jl,k,jb) = recv_buf(k+noffset(n),ik)
            ENDDO
          ENDDO
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END DO
#endif
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL
#endif
#endif

      !---------------------------------------------------------
    ENDIF  ! .NOT. my_process_is_mpi_seq()

!$ACC END DATA
    stop_sync_timer(timer_exch_data)

  END SUBROUTINE exchange_data_grf


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Interface for 2D arrays for exchange_data.
  !!
  !! Just reshapes the arrays and calls exchange_data.
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_r2d(p_pat, recv, send, add, l_recv_exists)
    !
    CLASS(t_comm_pattern_orig), INTENT(INOUT), TARGET :: p_pat
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_r2d"
    REAL(dp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_r2d_seq(p_pat, recv, send, add, l_recv_exists)
      RETURN
    END IF

!$ACC DATA CREATE( tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
    IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
!$ACC KERNELS PRESENT( tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
      tmp_recv(:,1,:) = 0._dp
!$ACC END KERNELS
    ELSE
!$ACC KERNELS PRESENT( recv, tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
      tmp_recv(:,1,:) = recv(:,:)
!$ACC END KERNELS
    ENDIF

    IF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
        CALL exchange_data_r3d(p_pat, tmp_recv, &
          send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
          add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
        CALL exchange_data_r3d(p_pat, tmp_recv, &
          send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
      ENDIF
    ELSE
      IF (PRESENT(add)) THEN
        CALL exchange_data_r3d(p_pat, tmp_recv, &
          add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
        CALL exchange_data_r3d(p_pat, tmp_recv)
      ENDIF
    ENDIF

!$ACC KERNELS PRESENT( recv, tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
    recv(:,:) = tmp_recv(:,1,:)
!$ACC END KERNELS

!$ACC END DATA

  END SUBROUTINE exchange_data_r2d


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Interface for 2D arrays for exchange_data.
  !!
  !! Just reshapes the arrays and calls exchange_data.
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_s2d(p_pat, recv, send, add, l_recv_exists)
    !
    CLASS(t_comm_pattern_orig), INTENT(INOUT), TARGET :: p_pat
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_s2d"
    REAL(sp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_s2d_seq(p_pat, recv, send, add, l_recv_exists)
      RETURN
    END IF

!$ACC DATA CREATE( tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
    IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
!$ACC KERNELS PRESENT( tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
      tmp_recv(:,1,:) = 0._sp
!$ACC END KERNELS
    ELSE
!$ACC KERNELS PRESENT( recv, tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
      tmp_recv(:,1,:) = recv(:,:)
!$ACC END KERNELS
    ENDIF

    IF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
        CALL exchange_data_s3d(p_pat, tmp_recv, &
          send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
          add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
        CALL exchange_data_s3d(p_pat, tmp_recv, &
          send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
      ENDIF
    ELSE
      IF (PRESENT(add)) THEN
        CALL exchange_data_s3d(p_pat, tmp_recv, &
          add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
        CALL exchange_data_s3d(p_pat, tmp_recv)
      ENDIF
    ENDIF

!$ACC KERNELS PRESENT( recv, tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
    recv(:,:) = tmp_recv(:,1,:)
!$ACC END KERNELS

!$ACC END DATA

  END SUBROUTINE exchange_data_s2d


  ! SEQUENTIAL version of subroutine "exchange_data_r3d"
  !
  SUBROUTINE exchange_data_r2d_seq(p_pat, recv, send, add, l_recv_exists)

    CLASS(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_r2d_seq"
    INTEGER :: i
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif

    ! consistency checks
    ! ------------------

    ! make sure that we are in sequential mode
    IF (.NOT. my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Internal error: sequential routine called in parallel run!")
    END IF
    ! further tests
    IF ( (p_pat%np_recv /= 1) .OR. (p_pat%np_send /= 1) ) THEN
      CALL finish(routine, "Internal error: inconsistent no. send/receive peers!")
    END IF
    IF ( (p_pat%recv_limits(1) - p_pat%recv_limits(0)) /= (p_pat%send_limits(1) - p_pat%send_limits(0)) ) THEN
      CALL finish(routine, "Internal error: inconsistent sender/receiver size!")
    END IF
    IF ( (p_pat%recv_limits(0) /= 0) .OR. (p_pat%send_limits(0) /= 0) ) THEN
      CALL finish(routine, "Internal error: inconsistent sender/receiver start position!")
    END IF
    IF ( (p_pat%recv_limits(1) /= p_pat%n_recv) .OR. (p_pat%n_recv /= p_pat%n_send) ) THEN
      CALL finish(routine, "Internal error: inconsistent counts for sender/receiver!")
    END IF

    ! "communication" (direct copy)
    ! -----------------------------

#ifdef _OPENACC
! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
#endif

!$ACC DATA COPYIN( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src ), &
!$ACC      IF ( i_am_accel_node .AND. acc_on )

    IF(PRESENT(add)) THEN
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( p_pat, send, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
      DO i=1,n_pnts
        recv( tmp_recv_dst_idx(i), tmp_recv_dst_blk(i) )  =                    &
          &  add( tmp_recv_dst_idx(i), tmp_recv_dst_blk(i) )                +  &
          &  send(tmp_send_src_idx(tmp_recv_src(i)),                                    &
          &       tmp_send_src_blk(tmp_recv_src(i)))
      END DO
!$ACC END PARALLEL
#else
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  add( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )                +  &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
#endif
    ELSE
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( send, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
      DO i=1,n_pnts
        recv( tmp_recv_dst_idx(i), tmp_recv_dst_blk(i) )  =                    &
          &  send(tmp_send_src_idx(tmp_recv_src(i)),                                    &
          &       tmp_send_src_blk(tmp_recv_src(i)))
      END DO
!$ACC END PARALLEL
#else
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
#endif
    END IF

!$ACC END DATA

  END SUBROUTINE exchange_data_r2d_seq


  ! SEQUENTIAL version of subroutine "exchange_data_r3d"
  !
  SUBROUTINE exchange_data_s2d_seq(p_pat, recv, send, add, l_recv_exists)

    CLASS(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_s2d_seq"
    INTEGER :: i
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif

    ! consistency checks
    ! ------------------

    ! make sure that we are in sequential mode
    IF (.NOT. my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Internal error: sequential routine called in parallel run!")
    END IF
    ! further tests
    IF ( (p_pat%np_recv /= 1) .OR. (p_pat%np_send /= 1) ) THEN
      CALL finish(routine, "Internal error: inconsistent no. send/receive peers!")
    END IF
    IF ( (p_pat%recv_limits(1) - p_pat%recv_limits(0)) /= (p_pat%send_limits(1) - p_pat%send_limits(0)) ) THEN
      CALL finish(routine, "Internal error: inconsistent sender/receiver size!")
    END IF
    IF ( (p_pat%recv_limits(0) /= 0) .OR. (p_pat%send_limits(0) /= 0) ) THEN
      CALL finish(routine, "Internal error: inconsistent sender/receiver start position!")
    END IF
    IF ( (p_pat%recv_limits(1) /= p_pat%n_recv) .OR. (p_pat%n_recv /= p_pat%n_send) ) THEN
      CALL finish(routine, "Internal error: inconsistent counts for sender/receiver!")
    END IF

    ! "communication" (direct copy)
    ! -----------------------------

#ifdef _OPENACC
! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
#endif

!$ACC DATA COPYIN( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src ), &
!$ACC      IF ( i_am_accel_node .AND. acc_on )

    IF(PRESENT(add)) THEN
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( send, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i=1,n_pnts
        recv( tmp_recv_dst_idx(i), tmp_recv_dst_blk(i) )  =                    &
          &  add( tmp_recv_dst_idx(i), tmp_recv_dst_blk(i) )                +  &
          &  send(tmp_send_src_idx(tmp_recv_src(i)),                                    &
          &       tmp_send_src_blk(tmp_recv_src(i)))
      END DO
!$ACC END PARALLEL
#else
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  add( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )                +  &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
#endif
    ELSE
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( send, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i=1,n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
!$ACC END PARALLEL
#else
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
#endif
    END IF

!$ACC END DATA

  END SUBROUTINE exchange_data_s2d_seq


  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_i2d(p_pat, recv, send, add, l_recv_exists)
    !
    CLASS(t_comm_pattern_orig), INTENT(INOUT), TARGET :: p_pat
    INTEGER, INTENT(INOUT), TARGET        :: recv(:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_i2d"
    INTEGER :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_i2d_seq(p_pat, recv, send, add)
      RETURN
    END IF

!$ACC DATA CREATE( tmp_recv ), IF ( i_am_accel_node .AND. acc_on )

    IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
!$ACC KERNELS PRESENT( tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
      tmp_recv(:,1,:) = 0
!$ACC END KERNELS
    ELSE
!$ACC KERNELS PRESENT( recv, tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
      tmp_recv(:,1,:) = recv(:,:)
!$ACC END KERNELS
    ENDIF

    IF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
        CALL exchange_data_i3d(p_pat, tmp_recv, &
          send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
          add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
        CALL exchange_data_i3d(p_pat, tmp_recv, &
          send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
      ENDIF
    ELSE
      IF (PRESENT(add)) THEN
        CALL exchange_data_i3d(p_pat, tmp_recv, &
          add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
        CALL exchange_data_i3d(p_pat, tmp_recv)
      ENDIF
    ENDIF

!$ACC KERNELS PRESENT( recv, tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
    recv(:,:) = tmp_recv(:,1,:)
!$ACC END KERNELS

!$ACC END DATA

  END SUBROUTINE exchange_data_i2d


  ! SEQUENTIAL version of subroutine "exchange_data_r3d"
  !
  SUBROUTINE exchange_data_i2d_seq(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    INTEGER, INTENT(INOUT), TARGET        :: recv(:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_i2d_seq"
    INTEGER :: i
#ifdef _OPENACC
! This is an embarrassing temporary hack to circumvent the OPENACC problems from the "EXTEND" of comm_pattern
    INTEGER   :: tmp_send_src_idx(p_pat%n_send), tmp_send_src_blk(p_pat%n_send)
    INTEGER   :: tmp_recv_dst_idx(p_pat%n_pnts), tmp_recv_dst_blk(p_pat%n_pnts)
    INTEGER   :: tmp_recv_src(p_pat%n_pnts)
    INTEGER   :: n_pnts, n_send
#endif

    ! consistency checks
    ! ------------------

    ! make sure that we are in sequential mode
    IF (.NOT. my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Internal error: sequential routine called in parallel run!")
    END IF
    ! further tests
    IF ( (p_pat%np_recv /= 1) .OR. (p_pat%np_send /= 1) ) THEN
      CALL finish(routine, "Internal error: inconsistent no. send/receive peers!")
    END IF
    IF ( (p_pat%recv_limits(1) - p_pat%recv_limits(0)) /= (p_pat%send_limits(1) - p_pat%send_limits(0)) ) THEN
      CALL finish(routine, "Internal error: inconsistent sender/receiver size!")
    END IF
    IF ( (p_pat%recv_limits(0) /= 0) .OR. (p_pat%send_limits(0) /= 0) ) THEN
      CALL finish(routine, "Internal error: inconsistent sender/receiver start position!")
    END IF
    IF ( (p_pat%recv_limits(1) /= p_pat%n_recv) .OR. (p_pat%n_recv /= p_pat%n_send) ) THEN
      CALL finish(routine, "Internal error: inconsistent counts for sender/receiver!")
    END IF

    ! "communication" (direct copy)
    ! -----------------------------

#ifdef _OPENACC
! 2018.11.13  Part of nasty hack to avoid the t_comm_pattern "EXTEND" problem for OpenACC
   n_pnts           = p_pat%n_pnts
   n_send           = p_pat%n_send
   tmp_send_src_idx = p_pat%send_src_idx
   tmp_send_src_blk = p_pat%send_src_blk
   tmp_recv_dst_idx = p_pat%recv_dst_idx
   tmp_recv_dst_blk = p_pat%recv_dst_blk
   tmp_recv_src     = p_pat%recv_src
#endif

!$ACC DATA COPYIN( tmp_send_src_idx, tmp_send_src_blk, tmp_recv_dst_idx, tmp_recv_dst_blk, tmp_recv_src ), &
!$ACC      IF ( i_am_accel_node .AND. acc_on )

    IF(PRESENT(add)) THEN
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( send, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
      DO i=1,n_pnts
        recv( tmp_recv_dst_idx(i), tmp_recv_dst_blk(i) )  =                    &
          &  add( tmp_recv_dst_idx(i), tmp_recv_dst_blk(i) )                +  &
          &  send(tmp_send_src_idx(tmp_recv_src(i)),                                    &
          &       tmp_send_src_blk(tmp_recv_src(i)))
      END DO
!$ACC END PARALLEL
#else
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  add( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )                +  &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
#endif
    ELSE
#ifdef _OPENACC
!$ACC PARALLEL PRESENT( p_pat, send, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
      DO i=1,n_pnts
        recv( tmp_recv_dst_idx(i), tmp_recv_dst_blk(i) )  =                    &
          &  send(tmp_send_src_idx(tmp_recv_src(i)),                                    &
          &       tmp_send_src_blk(tmp_recv_src(i)))
      END DO
!$ACC END PARALLEL
#else
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i)))
      END DO
#endif
    END IF
!$ACC END DATA

  END SUBROUTINE exchange_data_i2d_seq


  !-------------------------------------------------------------------------


  FUNCTION get_np_recv(comm_pat)
    CLASS (t_comm_pattern_orig), INTENT(IN) :: comm_pat
    INTEGER :: get_np_recv

    get_np_recv = comm_pat%np_recv
  END FUNCTION get_np_recv


  !-------------------------------------------------------------------------


  FUNCTION get_np_send(comm_pat)
    CLASS (t_comm_pattern_orig), INTENT(IN) :: comm_pat
    INTEGER :: get_np_send

    get_np_send = comm_pat%np_send
  END FUNCTION get_np_send


  !-------------------------------------------------------------------------


  SUBROUTINE get_pelist_recv(comm_pat, pelist_recv)
    CLASS (t_comm_pattern_orig), INTENT(IN) :: comm_pat
    INTEGER, INTENT(OUT) :: pelist_recv(:)

    pelist_recv = comm_pat%pelist_recv
  END SUBROUTINE get_pelist_recv

  !================================================================================================
  !
  ! Variant of exchange routine that expects input as a 1D unblocked vector and provides
  ! output as a blocked (nproma,nblks) 2D array.
  !
  SUBROUTINE exchange_data_r1d_2d(p_pat, recv, send)

    TYPE(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    REAL(dp), INTENT(INOUT), TARGET          :: recv(:,:)
    REAL(dp), INTENT(IN)                     :: send(:)

    REAL(dp) :: send_buf(p_pat%n_send), recv_buf(p_pat%n_recv)

    INTEGER :: i, il, k, np, irs, iss, pid, icount

    !-----------------------------------------------------------------------

    ! Set up irecv's for receive buffers
    DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

      pid    = p_pat%pelist_recv(np) ! ID of receiver PE
      irs    = p_pat%recv_startidx(np)
      icount = p_pat%recv_count(np)
      CALL p_irecv(recv_buf(irs), pid, 1, p_count=icount, comm=p_comm_work)

    ENDDO

    ! Set up send buffer
!$OMP PARALLEL DO PRIVATE(il)
    DO i = 1, p_pat%n_send
      il = (p_pat%send_src_blk(i)-1)*nproma + p_pat%send_src_idx(i)
      send_buf(i) = send(il)
    ENDDO
!$OMP END PARALLEL DO

    ! Send our data
    DO np = 1, p_pat%np_send ! loop over PEs where to send the data

      pid    = p_pat%pelist_send(np) ! ID of sender PE
      iss    = p_pat%send_startidx(np)
      icount = p_pat%send_count(np)
      CALL p_isend(send_buf(iss), pid, 1, p_count=icount, comm=p_comm_work)

    ENDDO

    CALL p_wait

    ! Fill in receive buffer
!$OMP PARALLEL DO
    DO i = 1, p_pat%n_pnts
      recv(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i)) = recv_buf(p_pat%recv_src(i))
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE exchange_data_r1d_2d

  !================================================================================================
  !
  ! Variant of exchange routine that expects input as a 1D unblocked vector and provides
  ! output as a blocked (nproma,nblks) 2D array.
  !
  SUBROUTINE exchange_data_s1d_2d(p_pat, recv, send)

    TYPE(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    REAL(sp), INTENT(INOUT), TARGET          :: recv(:,:)
    REAL(sp), INTENT(IN)                     :: send(:)

    REAL(sp) :: send_buf(p_pat%n_send), recv_buf(p_pat%n_recv)

    INTEGER :: i, il, k, np, irs, iss, pid, icount

    !-----------------------------------------------------------------------

    ! Set up irecv's for receive buffers
    DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

      pid    = p_pat%pelist_recv(np) ! ID of receiver PE
      irs    = p_pat%recv_startidx(np)
      icount = p_pat%recv_count(np)
      CALL p_irecv(recv_buf(irs), pid, 1, p_count=icount, comm=p_comm_work)

    ENDDO


    ! Set up send buffer
!$OMP PARALLEL DO PRIVATE(il)
    DO i = 1, p_pat%n_send
      il = (p_pat%send_src_blk(i)-1)*nproma + p_pat%send_src_idx(i)
      send_buf(i) = send(il)
    ENDDO
!$OMP END PARALLEL DO


    ! Send our data
    DO np = 1, p_pat%np_send ! loop over PEs where to send the data

      pid    = p_pat%pelist_send(np) ! ID of sender PE
      iss    = p_pat%send_startidx(np)
      icount = p_pat%send_count(np)
      CALL p_isend(send_buf(iss), pid, 1, p_count=icount, comm=p_comm_work)

    ENDDO

    CALL p_wait

    ! Fill in receive buffer
!$OMP PARALLEL DO
    DO i = 1, p_pat%n_pnts
      recv(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i)) = recv_buf(p_pat%recv_src(i))
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE exchange_data_s1d_2d

  !================================================================================================
  !
  ! Variant of exchange routine that expects input as a 1D unblocked vector and provides
  ! output as a blocked (nproma,nblks) 2D array.
  !
  SUBROUTINE exchange_data_i1d_2d(p_pat, recv, send)

    TYPE(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    INTEGER, INTENT(INOUT), TARGET          :: recv(:,:)
    INTEGER, INTENT(IN)                     :: send(:)

    INTEGER :: send_buf(p_pat%n_send), recv_buf(p_pat%n_recv)

    INTEGER :: i, il, k, np, irs, iss, pid, icount

    !-----------------------------------------------------------------------

    ! Set up irecv's for receive buffers
    DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

      pid    = p_pat%pelist_recv(np) ! ID of receiver PE
      irs    = p_pat%recv_startidx(np)
      icount = p_pat%recv_count(np)
      CALL p_irecv(recv_buf(irs), pid, 1, p_count=icount, comm=p_comm_work)

    ENDDO


    ! Set up send buffer
!$OMP PARALLEL DO PRIVATE(il)
    DO i = 1, p_pat%n_send
      il = (p_pat%send_src_blk(i)-1)*nproma + p_pat%send_src_idx(i)
      send_buf(i) = send(il)
    ENDDO
!$OMP END PARALLEL DO


    ! Send our data
    DO np = 1, p_pat%np_send ! loop over PEs where to send the data

      pid    = p_pat%pelist_send(np) ! ID of sender PE
      iss    = p_pat%send_startidx(np)
      icount = p_pat%send_count(np)
      CALL p_isend(send_buf(iss), pid, 1, p_count=icount, comm=p_comm_work)

    ENDDO

    CALL p_wait

    ! Fill in receive buffer
!$OMP PARALLEL DO
    DO i = 1, p_pat%n_pnts
      recv(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i)) = recv_buf(p_pat%recv_src(i))
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE exchange_data_i1d_2d

END MODULE mo_communication_orig
!
! Local Variables:
! f90-continuation-indent: 2
! End:
!
