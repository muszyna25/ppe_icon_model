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
USE mo_fortran_tools,        ONLY: t_ptr_3d, t_ptr_3d_sp, t_ptr_2d, t_ptr_1d_int, &
     &                             insert_dimension
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
    PROCEDURE :: exchange_data_l3d => exchange_data_l3d
    PROCEDURE :: exchange_data_r2d => exchange_data_r2d
    PROCEDURE :: exchange_data_s2d => exchange_data_s2d
    PROCEDURE :: exchange_data_i2d => exchange_data_i2d
    PROCEDURE :: exchange_data_l2d => exchange_data_l2d
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
  PUBLIC :: t_p_comm_pattern_orig

TYPE, EXTENDS(t_comm_pattern_collection) :: t_comm_pattern_collection_orig

   PRIVATE

   TYPE(t_p_comm_pattern_orig), ALLOCATABLE :: patterns(:)

   CONTAINS

   PROCEDURE :: setup => setup_comm_pattern_collection
   PROCEDURE :: delete => delete_comm_pattern_collection
   PROCEDURE :: exchange_data_grf => exchange_data_grf

END TYPE t_comm_pattern_collection_orig

#if defined( _OPENACC )
#define ACC_DEBUG NO_ACC
! #define __COMMUNICATION_NOACC
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

    CLASS(t_comm_pattern_orig), TARGET, INTENT(OUT) :: p_pat

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

!$ACC ENTER DATA COPYIN(p_pat, p_pat%recv_src, &
!$ACC                   p_pat%send_src_idx, p_pat%send_src_blk, &
!$ACC                   p_pat%recv_dst_idx, p_pat%recv_dst_blk) IF (acc_on)

  END SUBROUTINE setup_comm_pattern


  !-------------------------------------------------------------------------


  SUBROUTINE setup_comm_pattern2(p_pat, comm, recv_msg, send_msg, &
       glb2loc_index_recv, glb2loc_index_send, inplace)
    CLASS(t_comm_pattern_orig), TARGET, INTENT(out) :: p_pat
    INTEGER, INTENT(in) :: comm
    TYPE(xfer_list), INTENT(in) :: recv_msg(:), send_msg(:)
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: glb2loc_index_recv, &
         glb2loc_index_send
    LOGICAL, OPTIONAL, INTENT(in) :: inplace
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: recv_msg, send_msg
#endif
    INTEGER :: np_recv, np_send, n_send, n_recv, comm_rank, &
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

!$ACC ENTER DATA COPYIN(p_pat, p_pat%recv_src, &
!$ACC                   p_pat%send_src_idx, p_pat%send_src_blk, &
!$ACC                   p_pat%recv_dst_idx, p_pat%recv_dst_blk) IF (acc_on)

  CONTAINS
    SUBROUTINE count_msg_size(nmsg, msg, nself, nremote, is_inter, comm_rank, &
         ranks, counts, starts)
      INTEGER, INTENT(in) :: nmsg, comm_rank
      TYPE(xfer_list), INTENT(in) :: msg(nmsg)
      INTEGER, INTENT(out) :: nself, nremote, ranks(nmsg), counts(nmsg), &
           starts(nmsg)
      LOGICAL, INTENT(in) :: is_inter
      INTEGER :: nidx_remote, nidx_local, sz, msg_rank, sz_accum, i
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
      INTEGER :: nmsg, i, sz, msg_rank, jls, jl, jle
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

    CLASS(t_comm_pattern_orig), TARGET, INTENT(INOUT) :: p_pat

!$ACC EXIT DATA DELETE(p_pat%send_src_idx, p_pat%send_src_blk, &
!$ACC                  p_pat%recv_dst_idx, p_pat%recv_dst_blk, &
!$ACC                  p_pat%recv_src, p_pat) IF (acc_on)

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
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_r3d(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), TARGET, INTENT(INOUT) :: p_pat
    REAL(dp), INTENT(INOUT), TARGET           :: recv(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET    :: send(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET    :: add (:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_r3d"
    REAL(dp) :: send_buf(SIZE(recv,2),p_pat%n_send), &
      recv_buf(SIZE(recv,2),p_pat%n_recv)

    REAL(dp), POINTER :: send_ptr(:,:,:)
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)

    INTEGER :: i, k, np, irs, iss, pid, icount, ndim2
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)

    !-----------------------------------------------------------------------
    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_r3d_seq(p_pat, recv, send, add)
      RETURN
    END IF

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    IF(SIZE(recv,1) /= nproma) THEN
      CALL finish(routine,'Illegal first dimension of data array')
    ENDIF

    ndim2 = SIZE(recv,2)

!$ACC DATA CREATE( send_buf, recv_buf )                                                      &
!$ACC      PRESENT( recv, recv_src, recv_dst_blk, recv_dst_idx, send_src_blk, send_src_idx ) &
!$ACC      IF (use_gpu)

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
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(send_src_idx(i),1,send_src_blk(i))
      ENDDO
!$ACC END PARALLEL
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$NEC outerloop_unroll(4)
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k = 1, ndim2
        DO i = 1, p_pat%n_send
          send_buf(k,i) = send_ptr(send_src_idx(i),k,send_src_blk(i))
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
      DO i = 1, p_pat%n_send
        send_buf(1:ndim2,i) = send_ptr(send_src_idx(i),1:ndim2, send_src_blk(i))
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
    ENDIF

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ) IF (use_gpu)
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

    ! Wait for all outstanding requests to finish
    start_sync_timer(timer_exch_data_wait)
    CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ) IF (use_gpu)
#endif
    stop_sync_timer(timer_exch_data_wait)

    IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer

    IF(PRESENT(add)) THEN
      IF (ndim2 == 1) THEN
        k = 1
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),k,recv_dst_blk(i)) = &
            recv_buf(k,recv_src(i)) + add(recv_dst_idx(i),k,recv_dst_blk(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
        DO k = 1, ndim2
          DO i = 1, p_pat%n_pnts
            recv(recv_dst_idx(i),k,recv_dst_blk(i)) = &
              recv_buf(k,recv_src(i)) + add(recv_dst_idx(i),k,recv_dst_blk(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),:,recv_dst_blk(i)) = &
            recv_buf(:,recv_src(i)) + add(recv_dst_idx(i),1:ndim2,recv_dst_blk(i))
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
      ENDIF
    ELSE
      IF (ndim2 == 1) THEN
        k = 1
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),k,recv_dst_blk(i)) = recv_buf(k,recv_src(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
        DO k = 1, ndim2
          DO i = 1, p_pat%n_pnts
            recv(recv_dst_idx(i),k,recv_dst_blk(i)) = recv_buf(k,recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),:,recv_dst_blk(i)) = recv_buf(:,recv_src(i))
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
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_s3d(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), TARGET, INTENT(INOUT) :: p_pat
    REAL(sp), INTENT(INOUT), TARGET           :: recv(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET    :: send(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET    :: add (:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_s3d"
    REAL(sp) :: send_buf(SIZE(recv,2),p_pat%n_send), &
      recv_buf(SIZE(recv,2),p_pat%n_recv)

    REAL(sp), POINTER :: send_ptr(:,:,:)
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)

    INTEGER :: i, k, np, irs, iss, pid, icount, ndim2
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)

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
      CALL finish(routine,'Illegal first dimension of data array')
    ENDIF

    ndim2 = SIZE(recv,2)

!$ACC DATA CREATE( send_buf, recv_buf )                                                      &
!$ACC      PRESENT( recv, recv_src, recv_dst_blk, recv_dst_idx, send_src_blk, send_src_idx ) &
!$ACC      IF (use_gpu)

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
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(send_src_idx(i),1,send_src_blk(i))
      ENDDO
!$ACC END PARALLEL
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$NEC outerloop_unroll(4)
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k = 1, ndim2
        DO i = 1, p_pat%n_send
          send_buf(k,i) = send_ptr(send_src_idx(i),k,send_src_blk(i))
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
      DO i = 1, p_pat%n_send
        send_buf(1:ndim2,i) = send_ptr(send_src_idx(i),1:ndim2, send_src_blk(i))
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
    ENDIF

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ) IF (use_gpu)
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

    ! Wait for all outstanding requests to finish
    start_sync_timer(timer_exch_data_wait)
    CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ) IF (use_gpu)
#endif
    stop_sync_timer(timer_exch_data_wait)

    IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer

    IF(PRESENT(add)) THEN
      IF (ndim2 == 1) THEN
        k = 1
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),k,recv_dst_blk(i)) = &
            recv_buf(k,recv_src(i)) + add(recv_dst_idx(i),k,recv_dst_blk(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
        DO k = 1, ndim2
          DO i = 1, p_pat%n_pnts
            recv(recv_dst_idx(i),k,recv_dst_blk(i)) = &
              recv_buf(k,recv_src(i)) + add(recv_dst_idx(i),k,recv_dst_blk(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),:,recv_dst_blk(i)) = &
            recv_buf(:,recv_src(i)) + add(recv_dst_idx(i),1:ndim2,recv_dst_blk(i))
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
      ENDIF
    ELSE
      IF (ndim2 == 1) THEN
        k = 1
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),k,recv_dst_blk(i)) = recv_buf(k,recv_src(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
        DO k = 1, ndim2
          DO i = 1, p_pat%n_pnts
            recv(recv_dst_idx(i),k,recv_dst_blk(i)) = recv_buf(k,recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),:,recv_dst_blk(i)) = recv_buf(:,recv_src(i))
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
    INTEGER :: i, k, ndim2
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)

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

    IF(PRESENT(add)) THEN
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=1,ndim2
        DO i=1,p_pat%n_pnts
          recv( recv_dst_idx(i), k, recv_dst_blk(i) )  =                    &
          &  add( recv_dst_idx(i), k, recv_dst_blk(i) )                +  &
          &  send(send_src_idx(recv_src(i)),                                    &
          &       k,                                                                  &
          &       send_src_blk(recv_src(i)))
        ENDDO
      ENDDO
!$ACC END PARALLEL
    ELSE
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=1,ndim2
        DO i=1,p_pat%n_pnts
          recv( recv_dst_idx(i), k, recv_dst_blk(i) )  =                    &
          &  send(send_src_idx(recv_src(i)),                                    &
          &       k,                                                                  &
          &       send_src_blk(recv_src(i)))
        ENDDO
      ENDDO
!$ACC END PARALLEL
    END IF

  END SUBROUTINE exchange_data_r3d_seq


  ! SEQUENTIAL version of subroutine "exchange_data_s3d"
  !
  SUBROUTINE exchange_data_s3d_seq(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_s3d_seq"
    INTEGER :: i, k, ndim2
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)

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

    IF(PRESENT(add)) THEN
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=1,ndim2
        DO i=1,p_pat%n_pnts
          recv( recv_dst_idx(i), k, recv_dst_blk(i) )  =                    &
          &  add( recv_dst_idx(i), k, recv_dst_blk(i) )                +  &
          &  send(send_src_idx(recv_src(i)),                                    &
          &       k,                                                                  &
          &       send_src_blk(recv_src(i)))
        ENDDO
      ENDDO
!$ACC END PARALLEL
    ELSE
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=1,ndim2
        DO i=1,p_pat%n_pnts
          recv( recv_dst_idx(i), k, recv_dst_blk(i) )  =                    &
          &  send(send_src_idx(recv_src(i)),                                    &
          &       k,                                                                  &
          &       send_src_blk(recv_src(i)))
        ENDDO
      ENDDO
!$ACC END PARALLEL
    END IF

  END SUBROUTINE exchange_data_s3d_seq


  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_i3d(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), TARGET, INTENT(INOUT) :: p_pat
    INTEGER, INTENT(INOUT), TARGET           :: recv(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET    :: send(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET    :: add (:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_i3d"
    INTEGER :: send_buf(SIZE(recv,2),p_pat%n_send), &
      recv_buf(SIZE(recv,2),p_pat%n_recv)

    INTEGER, POINTER :: send_ptr(:,:,:)
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)

    INTEGER :: i, k, np, irs, iss, pid, icount, ndim2
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)

    !-----------------------------------------------------------------------
    IF(my_process_is_mpi_seq()) THEN
      CALL finish(routine, 'must not be called on single PE/test PE')
    END IF

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    IF(SIZE(recv,1) /= nproma) THEN
      CALL finish(routine,'Illegal first dimension of data array')
    ENDIF

    ndim2 = SIZE(recv,2)

!$ACC DATA CREATE( send_buf, recv_buf )                                                      &
!$ACC      PRESENT( recv, recv_src, recv_dst_blk, recv_dst_idx, send_src_blk, send_src_idx ) &
!$ACC      IF (use_gpu)

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
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(send_src_idx(i),1,send_src_blk(i))
      ENDDO
!$ACC END PARALLEL
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$NEC outerloop_unroll(4)
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k = 1, ndim2
        DO i = 1, p_pat%n_send
          send_buf(k,i) = send_ptr(send_src_idx(i),k,send_src_blk(i))
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
      DO i = 1, p_pat%n_send
        send_buf(1:ndim2,i) = send_ptr(send_src_idx(i),1:ndim2, send_src_blk(i))
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
    ENDIF

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ) IF (use_gpu)
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

    ! Wait for all outstanding requests to finish
    start_sync_timer(timer_exch_data_wait)
    CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ) IF (use_gpu)
#endif
    stop_sync_timer(timer_exch_data_wait)

    IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer

    IF(PRESENT(add)) THEN
      IF (ndim2 == 1) THEN
        k = 1
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),k,recv_dst_blk(i)) = &
            recv_buf(k,recv_src(i)) + add(recv_dst_idx(i),k,recv_dst_blk(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
        DO k = 1, ndim2
          DO i = 1, p_pat%n_pnts
            recv(recv_dst_idx(i),k,recv_dst_blk(i)) = &
              recv_buf(k,recv_src(i)) + add(recv_dst_idx(i),k,recv_dst_blk(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),:,recv_dst_blk(i)) = &
            recv_buf(:,recv_src(i)) + add(recv_dst_idx(i),1:ndim2,recv_dst_blk(i))
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
      ENDIF
    ELSE
      IF (ndim2 == 1) THEN
        k = 1
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),k,recv_dst_blk(i)) = recv_buf(k,recv_src(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
        DO k = 1, ndim2
          DO i = 1, p_pat%n_pnts
            recv(recv_dst_idx(i),k,recv_dst_blk(i)) = recv_buf(k,recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),:,recv_dst_blk(i)) = recv_buf(:,recv_src(i))
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
      ENDIF
    ENDIF

!$ACC END DATA

    stop_sync_timer(timer_exch_data)

  END SUBROUTINE exchange_data_i3d


  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_l3d(p_pat, recv, send)

    CLASS(t_comm_pattern_orig), TARGET, INTENT(INOUT) :: p_pat
    LOGICAL, INTENT(INOUT), TARGET           :: recv(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL, TARGET    :: send(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_l3d"
    LOGICAL :: send_buf(SIZE(recv,2),p_pat%n_send), &
      recv_buf(SIZE(recv,2),p_pat%n_recv)

    LOGICAL, POINTER :: send_ptr(:,:,:)
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)

    INTEGER :: i, k, np, irs, iss, pid, icount, ndim2
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)

    !-----------------------------------------------------------------------
    IF(my_process_is_mpi_seq()) THEN
      CALL finish(routine, 'must not be called on single PE/test PE')
    END IF

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    IF(SIZE(recv,1) /= nproma) THEN
      CALL finish(routine,'Illegal first dimension of data array')
    ENDIF

    ndim2 = SIZE(recv,2)

!$ACC DATA CREATE( send_buf, recv_buf )                                                      &
!$ACC      PRESENT( recv, recv_src, recv_dst_blk, recv_dst_idx, send_src_blk, send_src_idx ) &
!$ACC      IF (use_gpu)

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
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(send_src_idx(i),1,send_src_blk(i))
      ENDDO
!$ACC END PARALLEL
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$NEC outerloop_unroll(4)
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k = 1, ndim2
        DO i = 1, p_pat%n_send
          send_buf(k,i) = send_ptr(send_src_idx(i),k,send_src_blk(i))
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
      DO i = 1, p_pat%n_send
        send_buf(1:ndim2,i) = send_ptr(send_src_idx(i),1:ndim2, send_src_blk(i))
      ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
    ENDIF

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ) IF (use_gpu)
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

    ! Wait for all outstanding requests to finish
    start_sync_timer(timer_exch_data_wait)
    CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ) IF (use_gpu)
#endif
    stop_sync_timer(timer_exch_data_wait)

    IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer

    IF (ndim2 == 1) THEN
      k = 1
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i = 1, p_pat%n_pnts
        recv(recv_dst_idx(i),k,recv_dst_blk(i)) = recv_buf(k,recv_src(i))
      ENDDO
!$ACC END PARALLEL
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
      DO k = 1, ndim2
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),k,recv_dst_blk(i)) = recv_buf(k,recv_src(i))
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
        DO i = 1, p_pat%n_pnts
          recv(recv_dst_idx(i),:,recv_dst_blk(i)) = recv_buf(:,recv_src(i))
        ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
    ENDIF

!$ACC END DATA

    stop_sync_timer(timer_exch_data)

  END SUBROUTINE exchange_data_l3d

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

    CLASS(t_comm_pattern_orig), TARGET, INTENT(INOUT) :: p_pat
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
    INTEGER :: nfields, accum
    INTEGER :: i, k, kshift(SIZE(recv)), jb,ik, jl, n, np, irs, iss, pid, icount
    LOGICAL :: lsend
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)
    INTEGER :: n_send, n_pnts
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)
    n_send = p_pat%n_send
    n_pnts = p_pat%n_pnts

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

!$ACC DATA CREATE(send_buf, recv_buf) COPYIN(kshift)                                         &
!$ACC      PRESENT( recv_src, recv_dst_blk, recv_dst_idx, send_src_blk, send_src_idx ) &
!$ACC      IF (use_gpu)

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

      accum = 0
      DO n = 1, nfields
        noffset(n) = accum
        ndim2(n) = SIZE(recv(n)%p,2) - kshift(n)
        accum = accum + ndim2(n)
      ENDDO

!$ACC DATA COPYIN(noffset, ndim2) IF (use_gpu)

      ! Set up send buffer
#if defined( __SX__ ) || defined( _OPENACC )
      IF ( lsend ) THEN
        DO n = 1, nfields
          send_ptr => send(n)%p
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
          DO k = 1, ndim2(n)
            DO i = 1, n_send
              send_buf(k+noffset(n),i) = &
                send_ptr(send_src_idx(i),k+kshift(n),send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
      ELSE
        ! Send and receive arrays are identical (for boundary exchange)
        DO n = 1, nfields
          recv_ptr => recv(n)%p
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
          DO k = 1, ndim2(n)
            DO i = 1, n_send
              send_buf(k+noffset(n),i) = &
                recv_ptr(send_src_idx(i),k+kshift(n),send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
      ENDIF
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,n,k)
#endif
      DO i = 1, n_send
        jb = send_src_blk(i)
        jl = send_src_idx(i)
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
!$ACC UPDATE HOST( send_buf ) IF (use_gpu)
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

      ! Wait for all outstanding requests to finish
      start_sync_timer(timer_exch_data_wait)
      CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ) IF (use_gpu)
#endif
      stop_sync_timer(timer_exch_data_wait)

      IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
        start_sync_timer(timer_barrier)
        CALL p_barrier(p_pat%comm)
        stop_sync_timer(timer_barrier)
      ENDIF

      ! Fill in receive buffer

#if defined( __SX__ ) || defined( _OPENACC )
      DO n = 1, nfields
        recv_ptr => recv(n)%p
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
        DO k = 1, ndim2(n)
          DO i = 1, n_pnts
            recv_ptr(recv_dst_idx(i),k+kshift(n),recv_dst_blk(i)) =  &
              recv_buf(k+noffset(n),recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,ik,n,k)
#endif
      DO i = 1, n_pnts
        jb = recv_dst_blk(i)
        jl = recv_dst_idx(i)
        ik  = recv_src(i)
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

!$ACC END DATA

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

    CLASS(t_comm_pattern_orig), TARGET, INTENT(INOUT) :: p_pat

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
         jb, ik, jl, n, np, irs, iss, pid, icount, accum
    LOGICAL :: lsend
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)
    INTEGER :: n_send, n_pnts
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)
    n_send = p_pat%n_send
    n_pnts = p_pat%n_pnts

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

!$ACC DATA CREATE( send_buf_dp, recv_buf_dp, send_buf_sp, recv_buf_sp )                &
!$ACC      PRESENT( recv_src, recv_dst_blk, recv_dst_idx, send_src_blk, send_src_idx ) &
!$ACC      IF (use_gpu)

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

      accum = 0
      DO n = 1, nfields_dp
        noffset_dp(n) = accum
        ndim2_dp(n)   = SIZE(recv_dp(n)%p,2) - kshift_dp(n)
        accum = accum + ndim2_dp(n)
      ENDDO
      accum = 0
      DO n = 1, nfields_sp
        noffset_sp(n) = accum
        ndim2_sp(n)   = SIZE(recv_sp(n)%p,2) - kshift_sp(n)
        accum = accum + ndim2_sp(n)
      ENDDO

!$ACC DATA COPYIN(kshift_dp,noffset_dp,ndim2_dp,kshift_sp,noffset_sp,ndim2_sp)

      ! Set up send buffer
#if defined( __SX__ ) || defined( _OPENACC )
      IF ( lsend ) THEN
        DO n = 1, nfields_dp
          send_fld_dp => send_dp(n)%p   ! Refactoring for OpenACC
!$ACC PARALLEL DEFAULT(PRESENT) IF(use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
          DO k = 1, ndim2_dp(n)
            DO i = 1, n_send
              send_buf_dp(k+noffset_dp(n),i) = &
                send_fld_dp(send_src_idx(i),k+kshift_dp(n),send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
        DO n = 1, nfields_sp
          send_fld_sp => send_sp(n)%p   ! Refactoring for OpenACC
!$ACC PARALLEL DEFAULT(PRESENT) IF(use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
          DO k = 1, ndim2_sp(n)
            DO i = 1, n_send
              send_buf_sp(k+noffset_sp(n),i) = &
                send_fld_sp(send_src_idx(i),k+kshift_sp(n),send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
      ELSE
        ! Send and receive arrays are identical (for boundary exchange)
        DO n = 1, nfields_dp
          recv_fld_dp => recv_dp(n)%p
!$ACC PARALLEL DEFAULT(PRESENT) IF(use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
          DO k = 1, ndim2_dp(n)
            DO i = 1, n_send
              send_buf_dp(k+noffset_dp(n),i) = &
                recv_fld_dp(send_src_idx(i),k+kshift_dp(n),send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
        DO n = 1, nfields_sp
          recv_fld_sp => recv_sp(n)%p
!$ACC PARALLEL DEFAULT(PRESENT) IF(use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
          DO k = 1, ndim2_sp(n)
            DO i = 1, n_send
              send_buf_sp(k+noffset_sp(n),i) = &
                recv_fld_sp(send_src_idx(i),k+kshift_sp(n),send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
      ENDIF
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,n,k)
#endif
      DO i = 1, n_send
        jb = send_src_blk(i)
        jl = send_src_idx(i)
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
!$ACC UPDATE HOST( send_buf_dp, send_buf_sp ) IF (use_gpu)
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

      ! Wait for all outstanding requests to finish
      start_sync_timer(timer_exch_data_wait)
      CALL p_wait
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf_dp, recv_buf_sp ) IF (use_gpu)
#endif
      stop_sync_timer(timer_exch_data_wait)

      IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
        start_sync_timer(timer_barrier)
        CALL p_barrier(p_pat%comm)
        stop_sync_timer(timer_barrier)
      ENDIF

      ! Fill in receive buffer

#if defined( __SX__ ) || defined( _OPENACC )
      DO n = 1, nfields_dp
        recv_fld_dp => recv_dp(n)%p
!$ACC PARALLEL DEFAULT(PRESENT) IF(use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
!$NEC outerloop_unroll(4)
        DO k = 1, ndim2_dp(n)
          DO i = 1, n_pnts
            recv_fld_dp(recv_dst_idx(i),k+kshift_dp(n),recv_dst_blk(i)) =  &
              recv_buf_dp(k+noffset_dp(n),recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
      DO n = 1, nfields_sp
        recv_fld_sp => recv_sp(n)%p
!$ACC PARALLEL DEFAULT(PRESENT) IF(use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2) 
!$NEC outerloop_unroll(4)
        DO k = 1, ndim2_sp(n)
          DO i = 1, n_pnts
            recv_fld_sp(recv_dst_idx(i),k+kshift_sp(n),recv_dst_blk(i)) =  &
              recv_buf_sp(k+noffset_sp(n),recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,ik,n,k)
#endif
      DO i = 1, n_pnts
        jb = recv_dst_blk(i)
        jl = recv_dst_idx(i)
        ik  = recv_src(i)
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

!$ACC END DATA

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

    CLASS(t_comm_pattern_orig), TARGET, INTENT(INOUT) :: p_pat

    REAL(dp), INTENT(INOUT)           :: recv(:,:,:,:)
    REAL(dp), INTENT(IN   ), OPTIONAL :: send(:,:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_4de1"
    INTEGER, INTENT(IN)           :: nfields, ndim2tot

    INTEGER :: ndim2, koffset

    REAL(dp) :: send_buf(ndim2tot,p_pat%n_send),recv_buf(ndim2tot,p_pat%n_recv)

    INTEGER :: i, k, ik, jb, jl, n, np, irs, iss, pid, icount
    LOGICAL :: lsend
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)
    INTEGER :: n_send, n_pnts
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)
    n_send = p_pat%n_send
    n_pnts = p_pat%n_pnts

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

!$ACC DATA CREATE( send_buf, recv_buf )                                                      &
!$ACC      PRESENT( recv, recv_src, recv_dst_blk, recv_dst_idx, send_src_blk, send_src_idx ) &
!$ACC      IF (use_gpu)

    IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3)) THEN
      ! Set up irecv's for receive buffers
      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2tot
        CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_pat%comm)

      ENDDO
    ENDIF

#ifdef __SX__
    IF ( lsend ) THEN
      DO k = 1, ndim2
        koffset = (k-1)*nfields
!$NEC novector
        DO n = 1, nfields
          DO i = 1, n_send
            send_buf(n+koffset,i) = send(n,send_src_idx(i),k,send_src_blk(i))
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO k = 1, ndim2
        koffset = (k-1)*nfields
!$NEC novector
        DO n = 1, nfields
          DO i = 1, n_send
            send_buf(n+koffset,i) = recv(n,send_src_idx(i),k,send_src_blk(i))
          ENDDO
        ENDDO
      ENDDO
    ENDIF
#else
#if defined( __OMPPAR_COPY__ ) && !defined( _OPENACC )
!$OMP PARALLEL DO PRIVATE(jb,jl,koffset,k,n)
#endif
    DO i = 1, n_send
      jb = send_src_blk(i)
      jl = send_src_idx(i)
      IF ( lsend ) THEN
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k = 1, ndim2
          DO n = 1, nfields
            koffset = (k-1)*nfields
            send_buf(n+koffset,i) = send(n,jl,k,jb)
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ELSE
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k = 1, ndim2
          DO n = 1, nfields
            koffset = (k-1)*nfields
            send_buf(n+koffset,i) = recv(n,jl,k,jb)
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDIF
    ENDDO
#if defined( __OMPPAR_COPY__ ) && !defined( _OPENACC )
!$OMP END PARALLEL DO
#endif
#endif

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ) IF (use_gpu)
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

    ! Wait for all outstanding requests to finish
    start_sync_timer(timer_exch_data_wait)
    CALL p_wait
    stop_sync_timer(timer_exch_data_wait)

#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ) IF (use_gpu)
#endif

    IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL p_barrier(p_pat%comm)
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer
#ifdef __SX__
    DO k = 1, ndim2
      koffset = (k-1)*nfields
!$NEC novector
      DO n = 1, nfields
        DO i = 1, n_pnts
          recv(n,recv_dst_idx(i),k,recv_dst_blk(i)) = recv_buf(n+koffset,recv_src(i))
        ENDDO
      ENDDO
    ENDDO
#else
#if defined( __OMPPAR_COPY__ ) && !defined( _OPENACC )
!$OMP PARALLEL DO PRIVATE(jb,jl,ik,koffset,k,n)
#else
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG
#endif
    DO i = 1, n_pnts
      jb = recv_dst_blk(i)
      jl = recv_dst_idx(i)
      ik  = recv_src(i)
!$ACC LOOP VECTOR COLLAPSE(2)
      DO k = 1, ndim2
        DO n = 1, nfields
          koffset = (k-1)*nfields
          recv(n,jl,k,jb) = recv_buf(n+koffset,ik)
        ENDDO
      ENDDO
    ENDDO
#if defined( __OMPPAR_COPY__ ) && !defined( _OPENACC )
!$OMP END PARALLEL DO
#else
!$ACC END PARALLEL
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

    TYPE(t_ptr_3d) :: recv(nfields), send(nfields)

    TYPE(t_ptr_1d_int) :: p_send_src_idx(SIZE(p_pat_coll%patterns))
    TYPE(t_ptr_1d_int) :: p_send_src_blk(SIZE(p_pat_coll%patterns))
    TYPE(t_ptr_1d_int) :: p_recv_dst_idx(SIZE(p_pat_coll%patterns))
    TYPE(t_ptr_1d_int) :: p_recv_dst_blk(SIZE(p_pat_coll%patterns))
    TYPE(t_ptr_1d_int) :: p_recv_src(SIZE(p_pat_coll%patterns))
    INTEGER :: n_pnts(SIZE(p_pat_coll%patterns))
    INTEGER :: n_send(SIZE(p_pat_coll%patterns))

    INTEGER        :: ndim2(nfields), noffset(nfields), &
      ioffset_s(SIZE(p_pat_coll%patterns)), &
      ioffset_r(SIZE(p_pat_coll%patterns))

    REAL(dp), ALLOCATABLE :: send_buf(:,:),recv_buf(:,:), &
      auxs_buf(:,:),auxr_buf(:,:)

    INTEGER :: i, j, k, ik, jb, jl, n, np, irs, ire, iss, ise, &
      npats, isum, ioffset, isum1, n4d, pid, num_send, num_recv, &
      comm_size, idx_1d_i, accum, accum2
    INTEGER, ALLOCATABLE :: pelist_send(:), pelist_recv(:)

    TYPE(t_p_comm_pattern_orig), POINTER :: p_pat(:)
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

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

    ! Set pointers to input fields
    IF (PRESENT(recv4d1) .AND. .NOT. PRESENT(recv4d2)) THEN
      DO n = 1, nfields
        recv(n)%p => recv4d1(:,:,:,n)
        send(n)%p => send4d1(:,:,:,n)
      ENDDO
    ELSE IF (PRESENT(recv4d1) .AND. PRESENT(recv4d2)) THEN
      n4d = nfields/2
      DO n = 1, n4d
        recv(n)%p => recv4d1(:,:,:,n)
        send(n)%p => send4d1(:,:,:,n)
      ENDDO
      DO n = 1, n4d
        recv(n4d+n)%p => recv4d2(:,:,:,n)
        send(n4d+n)%p => send4d2(:,:,:,n)
      ENDDO
    ELSE
      IF (PRESENT(recv1)) THEN
        recv(1)%p => recv1
        send(1)%p => send1
      ENDIF
      IF (PRESENT(recv2)) THEN
        recv(2)%p => recv2
        send(2)%p => send2
      ENDIF
      IF (PRESENT(recv3)) THEN
        recv(3)%p => recv3
        send(3)%p => send3
      ENDIF
      IF (PRESENT(recv4)) THEN
        recv(4)%p => recv4
        send(4)%p => send4
      ENDIF
      IF (PRESENT(recv5)) THEN
        recv(5)%p => recv5
        send(5)%p => send5
      ENDIF
      IF (PRESENT(recv6)) THEN
        recv(6)%p => recv6
        send(6)%p => send6
      ENDIF
    ENDIF

    accum = 0
    DO n = 1, nfields
      noffset(n) = accum
      ndim2(n) = SIZE(recv(n)%p,2)
      accum = accum + ndim2(n)
    ENDDO

    accum = 0
    accum2 = 0
    DO np = 1, npats
      ioffset_r(np) = accum
      accum = accum + p_pat(np)%p%n_recv
      ioffset_s(np) = accum2
      accum2 = accum2 + p_pat(np)%p%n_send
    ENDDO

    DO np = 1, npats
      p_send_src_idx(np)%p => p_pat(np)%p%send_src_idx
      p_send_src_blk(np)%p => p_pat(np)%p%send_src_blk
      p_recv_dst_idx(np)%p => p_pat(np)%p%recv_dst_idx
      p_recv_dst_blk(np)%p => p_pat(np)%p%recv_dst_blk
      p_recv_src(np)%p => p_pat(np)%p%recv_src
      n_pnts(np) = p_pat(np)%p%n_pnts
      n_send(np) = p_pat(np)%p%n_send
    END DO

!$ACC DATA CREATE( send_buf, recv_buf, auxs_buf, auxr_buf ) &
!$ACC      COPYIN( ndim2, noffset, ioffset_s, ioffset_r, &
!$ACC              send, recv, n_pnts, n_send, p_recv_src, &
!$ACC              p_send_src_idx, p_send_src_blk, &
!$ACC              p_recv_dst_idx, p_recv_dst_blk ) IF (use_gpu)

#ifdef _OPENACC
     DO n = 1, nfields
       DO np = 1, npats
!$ACC ENTER DATA ATTACH( send(np+(n-1)*npats)%p )
       ENDDO
!$ACC ENTER DATA ATTACH( recv(n)%p )
     ENDDO

     DO np = 1, npats
!$ACC ENTER DATA ATTACH( p_send_src_idx(np)%p, p_send_src_blk(np)%p, &
!$ACC                    p_recv_dst_idx(np)%p, p_recv_dst_blk(np)%p, &
!$ACC                    p_recv_src(np)%p ) IF (use_gpu)
     ENDDO
#endif

    IF (my_process_is_mpi_seq()) THEN

!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG
      DO np = 1, npats
!$ACC LOOP VECTOR
        DO i = 1, n_pnts(np)
          idx_1d_i = idx_1d(p_send_src_idx(np)%p(p_recv_src(np)%p(i)),       &
                            p_send_src_blk(np)%p(p_recv_src(np)%p(i)))
          DO n = 1, nfields
            DO k = 1, ndim2(n)
              recv(n)%p( p_recv_dst_idx(np)%p(i), k, &
                p_recv_dst_blk(np)%p(i) ) =            &
                send(n)%p(k, idx_1d_i, np)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL

    ELSE    ! WS: removed RETURN statement to allow for OpenACC DATA region

      !-----------------------------------------------------------------------
      ! Set up irecv's for receive buffers
      IF (iorder_sendrecv <= 1 .OR. iorder_sendrecv >= 3) THEN

        ioffset = 0
        DO np = 1, num_recv ! loop over PEs from where to receive the data

          pid = pelist_recv(np) ! ID of receiver PE

          ! Sum up receive points over all communication patterns to be processed
          isum = ioffset
          DO n = 1, npats
            isum = isum + p_pat(n)%p%recv_limits(pid+1) - &
              p_pat(n)%p%recv_limits(pid)
          ENDDO

!TODO: this will probably not work with __USE_G2G
          IF(isum > ioffset) &
            CALL p_irecv(auxr_buf(1,ioffset+1), pid, 1, &
            &            p_count=(isum-ioffset)*ndim2tot, &
            &            comm=p_pat_coll%patterns(1)%p%comm)
          ioffset = isum

        ENDDO

      ENDIF

      ! Set up send buffer
#if defined( __SX__ ) || defined( _OPENACC )

!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG
      DO n = 1, nfields
!$ACC LOOP VECTOR
        DO np = 1, npats
!$NEC novector
          DO k = 1, ndim2(n)
            DO i = 1, n_send(np)
              send_buf(k+noffset(n),i+ioffset_s(np)) =                &
                & send(n)%p(k, idx_1d(p_send_src_idx(np)%p(i),    &
                &                     p_send_src_blk(np)%p(i)), np)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL PRIVATE(np)
#endif
      DO np = 1, npats
#ifdef __OMPPAR_COPY__
!$OMP DO PRIVATE(idx_1d_i,n,k)
#endif
        DO i = 1, n_send(np)
          idx_1d_i = idx_1d(p_send_src_idx(np)%p(i), p_send_src_blk(np)%p(i))
          DO n = 1, nfields
            DO k = 1, ndim2(n)
              send_buf(k+noffset(n),i+ioffset_s(np)) = send(n)%p(k, idx_1d_i, np)
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
!$ACC KERNELS DEFAULT(PRESENT) IF (use_gpu)
!
!  TODO:  Makes sure this is set up correctly
              auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
!$ACC END KERNELS
              isum = isum+isum1
            ENDIF
          ENDDO

#ifndef __USE_G2G
!$ACC UPDATE HOST( auxs_buf(:,ioffset+1:ioffset+isum) ) IF (use_gpu)
#endif

!TODO: this will probably not work with __USE_G2G
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
!$ACC KERNELS DEFAULT(PRESENT) IF (use_gpu)
              auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
!$ACC END KERNELS
              isum = isum+isum1
            ENDIF
          ENDDO

#ifndef __USE_G2G
!$ACC UPDATE HOST( auxs_buf(:,ioffset+1:ioffset+isum) ) IF (use_gpu)
#endif

!TODO: this will probably not work with __USE_G2G
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

!TODO: this will probably not work with __USE_G2G
          IF(isum > ioffset) CALL p_recv(auxr_buf(1,ioffset+1), pid, 1,             &
            p_count=(isum-ioffset)*ndim2tot, comm=p_pat_coll%patterns(1)%p%comm)
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
!$ACC KERNELS DEFAULT(PRESENT) IF (use_gpu)
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
!$ACC UPDATE HOST( auxs_buf(:,ioffset+1:ioffset+isum) ) IF (use_gpu)
#endif
!$OMP MASTER
!TODO: this will probably not work with __USE_G2G
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

#ifndef __USE_G2G
!$ACC UPDATE DEVICE( auxr_buf ) IF (use_gpu)
#endif

      ! Copy exchanged data back to receive buffer

#ifdef __OMPPAR_COPY__
!$OMP PARALLEL PRIVATE(ioffset,pid,isum,irs,ire,isum1,n,i)
#endif

      ioffset = 0
      DO np = 1, num_recv ! loop over PEs from where to receive the data

        pid = pelist_recv(np) ! ID of receiver PE

        isum = ioffset
        DO n = 1, npats
          irs = p_pat(n)%p%recv_limits(pid)   + ioffset_r(n)
          ire = p_pat(n)%p%recv_limits(pid+1) + ioffset_r(n)
          isum1 = ire - irs
          IF (isum1 > 0) THEN
#ifdef __OMPPAR_COPY__
!$OMP DO
#endif
!
!  TODO:  check compiler output:   KERNELS may decide to perform this on the CPU...
!         replace KERNELS with PARALLEL
!
!$ACC KERNELS DEFAULT(PRESENT) IF (use_gpu)
            DO i = 1, isum1
              recv_buf(:,irs+i) = auxr_buf(:,isum+i)
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

#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG
      DO n = 1, nfields
        DO np = 1, npats
!$ACC LOOP VECTOR
!$NEC novector
          DO k = 1, ndim2(n)
            DO i = 1, n_pnts(np)
              recv(n)%p(p_recv_dst_idx(np)%p(i),k, &
                        p_recv_dst_blk(np)%p(i)) =   &
                recv_buf(k+noffset(n),p_recv_src(np)%p(i)+ioffset_r(np))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
      DO np = 1, npats
#ifdef __OMPPAR_COPY__
!$OMP DO PRIVATE(jb,jl,ik,n,k)
#endif
        DO i = 1, n_pnts(np)
          jb = p_recv_dst_blk(np)%p(i)
          jl = p_recv_dst_idx(np)%p(i)
          ik  = p_recv_src(np)%p(i)+ioffset_r(np)
          DO n = 1, nfields
            DO k = 1, ndim2(n)
              recv(n)%p(jl,k,jb) = recv_buf(k+noffset(n),ik)
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

#ifdef _OPENACC
     DO n = 1, nfields
!$ACC EXIT DATA DETACH( recv(n)%p ) IF (use_gpu)
       DO np = 1, npats
!$ACC EXIT DATA DETACH( send(np+(n-1)*npats)%p ) IF (use_gpu)
       ENDDO
     ENDDO

     DO np = 1, npats
!$ACC EXIT DATA DETACH( p_send_src_idx(np)%p, p_send_src_blk(np)%p, &
!$ACC                   p_recv_dst_idx(np)%p, p_recv_dst_blk(np)%p, &
!$ACC                   p_recv_src(np)%p ) IF (use_gpu)
     ENDDO
#endif

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
    REAL(dp), POINTER :: recv3d(:,:,:)
    REAL(dp), POINTER :: send3d(:,:,:)
    REAL(dp), POINTER :: add3d(:,:,:)

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_r2d_seq(p_pat, recv, send, add)
      RETURN
    END IF

    IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) recv = 0._dp

    CALL insert_dimension(recv3d, recv, 2)
    IF (PRESENT(send)) CALL insert_dimension(send3d, send, 2)
    IF (PRESENT(add))  CALL insert_dimension(add3d, add, 2)

    IF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
        CALL exchange_data_r3d(p_pat, recv3d, send=send3d, add=add3d)
      ELSE
        CALL exchange_data_r3d(p_pat, recv3d, send=send3d)
      ENDIF
    ELSE
      IF (PRESENT(add)) THEN
        CALL exchange_data_r3d(p_pat, recv3d, add=add3d)
      ELSE
        CALL exchange_data_r3d(p_pat, recv3d)
      ENDIF
    ENDIF

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
    REAL(sp), POINTER :: recv3d(:,:,:)
    REAL(sp), POINTER :: send3d(:,:,:)
    REAL(sp), POINTER :: add3d(:,:,:)

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_s2d_seq(p_pat, recv, send, add)
      RETURN
    END IF

    IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) recv = 0._sp

    CALL insert_dimension(recv3d, recv, 2)
    IF (PRESENT(send)) CALL insert_dimension(send3d, send, 2)
    IF (PRESENT(add))  CALL insert_dimension(add3d, add, 2)

    IF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
        CALL exchange_data_s3d(p_pat, recv3d, send=send3d, add=add3d)
      ELSE
        CALL exchange_data_s3d(p_pat, recv3d, send=send3d)
      ENDIF
    ELSE
      IF (PRESENT(add)) THEN
        CALL exchange_data_s3d(p_pat, recv3d, add=add3d)
      ELSE
        CALL exchange_data_s3d(p_pat, recv3d)
      ENDIF
    ENDIF

  END SUBROUTINE exchange_data_s2d


  ! SEQUENTIAL version of subroutine "exchange_data_r3d"
  !
  SUBROUTINE exchange_data_r2d_seq(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_r2d_seq"
    INTEGER :: i
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)
    INTEGER :: n_pnts
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)
    n_pnts = p_pat%n_pnts

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

    IF(PRESENT(add)) THEN
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i=1,n_pnts
        recv( recv_dst_idx(i), recv_dst_blk(i) )  =    &
          &  add( recv_dst_idx(i), recv_dst_blk(i) ) + &
          &  send(send_src_idx(recv_src(i)),           &
          &       send_src_blk(recv_src(i)))
      END DO
!$ACC END PARALLEL
    ELSE
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i=1,n_pnts
        recv( recv_dst_idx(i), recv_dst_blk(i) )  = &
          &  send(send_src_idx(recv_src(i)),        &
          &       send_src_blk(recv_src(i)))
      END DO
!$ACC END PARALLEL
    END IF

  END SUBROUTINE exchange_data_r2d_seq


  ! SEQUENTIAL version of subroutine "exchange_data_r3d"
  !
  SUBROUTINE exchange_data_s2d_seq(p_pat, recv, send, add)

    CLASS(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_s2d_seq"
    INTEGER :: i
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)
    INTEGER :: n_pnts
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)
    n_pnts = p_pat%n_pnts

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

    IF(PRESENT(add)) THEN
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i=1,n_pnts
        recv( recv_dst_idx(i), recv_dst_blk(i) )  =    &
          &  add( recv_dst_idx(i), recv_dst_blk(i) ) + &
          &  send(send_src_idx(recv_src(i)),           &
          &       send_src_blk(recv_src(i)))
      END DO
!$ACC END PARALLEL
    ELSE
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i=1,n_pnts
        recv( recv_dst_idx(i), recv_dst_blk(i) )  = &
          &  send(send_src_idx(recv_src(i)),        &
          &       send_src_blk(recv_src(i)))
      END DO
!$ACC END PARALLEL
    END IF

  END SUBROUTINE exchange_data_s2d_seq


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
    INTEGER, POINTER :: recv3d(:,:,:)
    INTEGER, POINTER :: send3d(:,:,:)
    INTEGER, POINTER :: add3d(:,:,:)

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_i2d_seq(p_pat, recv, send, add)
      RETURN
    END IF

    IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) recv = 0

    CALL insert_dimension(recv3d, recv, 2)
    IF (PRESENT(send)) CALL insert_dimension(send3d, send, 2)
    IF (PRESENT(add))  CALL insert_dimension(add3d, add, 2)

    IF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
        CALL exchange_data_i3d(p_pat, recv3d, send=send3d, add=add3d)
      ELSE
        CALL exchange_data_i3d(p_pat, recv3d, send=send3d)
      ENDIF
    ELSE
      IF (PRESENT(add)) THEN
        CALL exchange_data_i3d(p_pat, recv3d, add=add3d)
      ELSE
        CALL exchange_data_i3d(p_pat, recv3d)
      ENDIF
    ENDIF

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
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)
    INTEGER :: n_pnts
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)
    n_pnts = p_pat%n_pnts

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

    IF(PRESENT(add)) THEN
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i=1,n_pnts
        recv( recv_dst_idx(i), recv_dst_blk(i) )  =    &
          &  add( recv_dst_idx(i), recv_dst_blk(i) ) + &
          &  send(send_src_idx(recv_src(i)),           &
          &       send_src_blk(recv_src(i)))
      END DO
!$ACC END PARALLEL
    ELSE
!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
      DO i=1,n_pnts
        recv( recv_dst_idx(i), recv_dst_blk(i) )  = &
          &  send(send_src_idx(recv_src(i)),        &
          &       send_src_blk(recv_src(i)))
      END DO
!$ACC END PARALLEL
    END IF

  END SUBROUTINE exchange_data_i2d_seq


  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_l2d(p_pat, recv, send, l_recv_exists)
    !
    CLASS(t_comm_pattern_orig), INTENT(INOUT), TARGET :: p_pat
    LOGICAL, INTENT(INOUT), TARGET        :: recv(:,:)
    LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    LOGICAL, OPTIONAL :: l_recv_exists

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_l2d"
    LOGICAL, POINTER :: recv3d(:,:,:)
    LOGICAL, POINTER :: send3d(:,:,:)

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_l2d_seq(p_pat, recv, send)
      RETURN
    END IF

    IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) recv = .FALSE.

    CALL insert_dimension(recv3d, recv, 2)
    IF (PRESENT(send)) CALL insert_dimension(send3d, send, 2)

    IF (PRESENT(send)) THEN
      CALL exchange_data_l3d(p_pat, recv3d, send=send3d)
    ELSE
      CALL exchange_data_l3d(p_pat, recv3d)
    ENDIF

  END SUBROUTINE exchange_data_l2d


  ! SEQUENTIAL version of subroutine "exchange_data_l3d"
  !
  SUBROUTINE exchange_data_l2d_seq(p_pat, recv, send)

    CLASS(t_comm_pattern_orig), INTENT(IN), TARGET :: p_pat
    LOGICAL, INTENT(INOUT), TARGET        :: recv(:,:)
    LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_l2d_seq"
    INTEGER :: i
    INTEGER, POINTER :: recv_src(:)
    INTEGER, POINTER :: recv_dst_blk(:)
    INTEGER, POINTER :: recv_dst_idx(:)
    INTEGER, POINTER :: send_src_blk(:)
    INTEGER, POINTER :: send_src_idx(:)
    INTEGER :: n_pnts
#ifdef _OPENACC
    LOGICAL :: use_gpu

    use_gpu = i_am_accel_node .AND. acc_on
#endif

    recv_src => p_pat%recv_src(:)
    recv_dst_blk => p_pat%recv_dst_blk(:)
    recv_dst_idx => p_pat%recv_dst_idx(:)
    send_src_blk => p_pat%send_src_blk(:)
    send_src_idx => p_pat%send_src_idx(:)
    n_pnts = p_pat%n_pnts

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

!$ACC PARALLEL DEFAULT(PRESENT) IF (use_gpu)
!$ACC LOOP GANG VECTOR
    DO i=1,n_pnts
      recv( recv_dst_idx(i), recv_dst_blk(i) )  = &
        &  send(send_src_idx(recv_src(i)),        &
        &       send_src_blk(recv_src(i)))
    END DO
!$ACC END PARALLEL

  END SUBROUTINE exchange_data_l2d_seq


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

END MODULE mo_communication_orig
!
! Local Variables:
! f90-continuation-indent: 2
! End:
!
