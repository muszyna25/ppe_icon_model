!>
!!               This module provides the communication routines.
!!
!!               This module provides the communication routines
!! for parallel runs
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_communication
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

USE mo_kind,               ONLY: wp
USE mo_exception,          ONLY: finish
USE mo_mpi,                ONLY: p_send, p_recv, p_irecv, p_wait, p_isend, &
     & p_real_dp, p_int, p_bool, p_comm_work, &
     & my_process_is_mpi_seq, &
     & p_pe_work, p_n_work
USE mo_parallel_config, ONLY: iorder_sendrecv, nproma
USE mo_timer,           ONLY: timer_start, timer_stop, activate_sync_timers, &
  & timer_exch_data, timer_exch_data_rv, timer_exch_data_async


IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

!modules interface-------------------------------------------
!subroutines
PUBLIC :: blk_no, idx_no, idx_1d
PUBLIC :: setup_comm_pattern, delete_comm_pattern, exchange_data, exchange_data_reverse,   &
          exchange_data_mult, exchange_data_grf, exchange_data_gm,    &
          start_delayed_exchange, do_delayed_exchange,                &
          start_async_comm, complete_async_comm, exchange_data_4de3
PUBLIC :: t_comm_pattern
!
!variables

!--------------------------------------------------------------------------------------------------
!
TYPE t_comm_pattern

   ! Number of points we receive in communication,
   ! this is the same as recv_limits

   INTEGER :: n_recv ! Number of points we receive from other PEs
   INTEGER :: n_pnts ! Number of points we output into local array;
                     ! this may be bigger than n_recv due to
                     ! duplicate entries
   INTEGER :: n_send ! Number of points we send to other PEs

   INTEGER :: np_recv ! Number of PEs from which data have to be received
   INTEGER :: np_send ! Number of PEs to which data have to be sent

   INTEGER, ALLOCATABLE :: recv_limits(:)

   INTEGER, ALLOCATABLE :: recv_src(:)
   INTEGER, ALLOCATABLE :: recv_dst_blk(:)
   INTEGER, ALLOCATABLE :: recv_dst_idx(:)

   INTEGER, ALLOCATABLE :: send_limits(:)

   INTEGER, ALLOCATABLE :: send_src_blk(:)
   INTEGER, ALLOCATABLE :: send_src_idx(:)

   INTEGER, ALLOCATABLE :: pelist_send(:)
   INTEGER, ALLOCATABLE :: pelist_recv(:)

   INTEGER, ALLOCATABLE :: send_startidx(:)
   INTEGER, ALLOCATABLE :: recv_startidx(:)

   INTEGER, ALLOCATABLE :: send_count(:)
   INTEGER, ALLOCATABLE :: recv_count(:)

END TYPE t_comm_pattern

!--------------------------------------------------------------------------------------------------
!
TYPE t_buffer_r
   INTEGER :: nelems
   REAL(wp), ALLOCATABLE :: buf(:)
END TYPE t_buffer_r

TYPE(t_buffer_r), ALLOCATABLE :: send_bufs_r(:)

TYPE t_buffer_i
   INTEGER :: nelems
   INTEGER, ALLOCATABLE :: buf(:)
END TYPE t_buffer_i

TYPE(t_buffer_i), ALLOCATABLE :: send_bufs_i(:)

TYPE t_buffer_l
   INTEGER :: nelems
   LOGICAL, ALLOCATABLE :: buf(:)
END TYPE t_buffer_l

TYPE(t_buffer_l), ALLOCATABLE :: send_bufs_l(:)

!--------------------------------------------------------------------------------------------------
!
TYPE t_request_r
   TYPE(t_comm_pattern), POINTER :: p_pat
   INTEGER :: ndim2
   LOGICAL :: reverse
   REAL(wp), POINTER :: recv2(:,:)
   REAL(wp), POINTER :: recv3(:,:,:)
   REAL(wp), POINTER :: add2(:,:)
   REAL(wp), POINTER :: add3(:,:,:)
END TYPE t_request_r

TYPE t_request_i
   TYPE(t_comm_pattern), POINTER :: p_pat
   INTEGER :: ndim2
   LOGICAL :: reverse
   INTEGER, POINTER :: recv2(:,:)
   INTEGER, POINTER :: recv3(:,:,:)
   INTEGER, POINTER :: add2(:,:)
   INTEGER, POINTER :: add3(:,:,:)
END TYPE t_request_i

TYPE t_request_l
   TYPE(t_comm_pattern), POINTER :: p_pat
   INTEGER :: ndim2
   LOGICAL :: reverse
   LOGICAL, POINTER :: recv2(:,:)
   LOGICAL, POINTER :: recv3(:,:,:)
END TYPE t_request_l

INTEGER, PARAMETER :: max_delayed_requests = 10000

TYPE(t_request_r) :: delayed_request_r(max_delayed_requests)
TYPE(t_request_i) :: delayed_request_i(max_delayed_requests)
TYPE(t_request_l) :: delayed_request_l(max_delayed_requests)

INTEGER :: n_delayed_requests_r = 0
INTEGER :: n_delayed_requests_i = 0
INTEGER :: n_delayed_requests_l = 0

LOGICAL :: use_exchange_delayed = .FALSE.

!--------------------------------------------------------------------------------------------------
!

INTERFACE exchange_data
   MODULE PROCEDURE exchange_data_r3d
   MODULE PROCEDURE exchange_data_i3d
   MODULE PROCEDURE exchange_data_l3d
   MODULE PROCEDURE exchange_data_r2d
   MODULE PROCEDURE exchange_data_i2d
   MODULE PROCEDURE exchange_data_l2d
END INTERFACE

INTERFACE exchange_data_reverse
   MODULE PROCEDURE exchange_data_reverse_3
   MODULE PROCEDURE exchange_data_reverse_2
END INTERFACE

INTERFACE buffer_data
  MODULE PROCEDURE buffer_data_r
  MODULE PROCEDURE buffer_data_i
  MODULE PROCEDURE buffer_data_l
END INTERFACE buffer_data

!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
! The following functions are for conversion of 1D to 2D indices and vice versa
!
! Treatment of 0 (important for empty patches) and negative numbers:
!
! Converting 1D => 2D:
!
! 0 always is mapped to blk_no = 1, idx_no = 0
! negative numbers: Convert usings ABS(j) and negate idx_no
!
! Thus: blk_no >= 1 always!
!       idx_no > 0  for j > 0
!       idx_no = 0  for j = 0
!       idx_no < 0  for j < 0
!
! This mimics mostly the behaviour of reshape_idx in mo_model_domimp_patches
! with a difference for nproma=1 and j=0 (where reshape_idx returns blk_no=0, idx_no=1)
!
! The consisten treatment of 0 in the above way is very important for empty patches
! where start_index=1, end_index=0
!
! Converting 2D => 1D:
! Trying to invert the above and catching cases with blk_no < 1
!
!-------------------------------------------------------------------------
ELEMENTAL INTEGER FUNCTION blk_no(j)
  INTEGER, INTENT(in) :: j
  blk_no = MAX((ABS(j)-1)/nproma + 1, 1) ! i.e. also 1 for j=0, nproma=1
END FUNCTION blk_no
!-------------------------------------------------------------------------
ELEMENTAL INTEGER FUNCTION idx_no(j)
  INTEGER, INTENT(in) :: j
  IF(j==0) THEN
    idx_no = 0
  ELSE
    idx_no = SIGN(MOD(ABS(j)-1,nproma)+1, j)
  ENDIF
END FUNCTION idx_no
!-------------------------------------------------------------------------
ELEMENTAL INTEGER FUNCTION idx_1d(jl,jb)
  INTEGER, INTENT(in) :: jl, jb
  IF(jb<=0) THEN
    idx_1d = 0 ! This covers the special case nproma==1,jb=0,jl=1
               ! All other cases are invalid and get also a 0 returned
  ELSE
    idx_1d = SIGN((jb-1)*nproma + ABS(jl), jl)
  ENDIF
END FUNCTION idx_1d
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!
!

!>
!! Sets up a communication pattern for exchanging data.
!!
!! n_points       Total number of points in the RECEIVER array,
!!                not every point is necessarily set during exchange
!!                (see owner!)
!!
!! owner          Owner PE number of every point in the RECEIVER array,
!!                if owner(.) == -1, this point will not be set during exchange.
!!                If owner(.) == p_pe, this point will be exchanged,
!!                this is necessary if sender and receiver arrays are
!!                different (e.g. feedback, gather, scatter)
!!
!! global_index   Global index of of every point in the RECEIVER array
!!                There may be more tha 1 point with the same global index,
!!                in this case the point is exchanged only once and
!!                locally distributed.
!!
!! local_index    Local index in the SENDER array.
!!                This array must have the local index at the global
!!                position for points which are local and a negative
!!                value at every position which is not owned by the local PE.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
SUBROUTINE setup_comm_pattern(n_points, owner, global_index, local_index, p_pat)

!

   INTEGER, INTENT(IN) :: n_points         ! Total number of points
   INTEGER, INTENT(IN) :: owner(:)         ! Owner of every point
   INTEGER, INTENT(IN) :: global_index(:)  ! Global index of every point
   INTEGER, INTENT(IN) :: local_index(:)   ! Array mapping global indices to local ones
                                           ! valid on the remote side!

   TYPE(t_comm_pattern), INTENT(INOUT) :: p_pat



   INTEGER, ALLOCATABLE :: icnt(:), flag(:), global_recv_index(:), send_src(:)

   INTEGER :: i, n, np, nr, num_recv, irs, ire, num_send, iss, ise, max_glb

!-----------------------------------------------------------------------

   if(my_process_is_mpi_seq()) &
      CALL finish('setup_comm_pattern','must not be called on single PE/test PE')

   ALLOCATE(icnt(0:p_n_work-1))
   max_glb = MAX(MAXVAL(ABS(global_index(1:n_points)),mask=(owner(1:n_points)>=0)),1)
   ALLOCATE(flag(max_glb))

   ! Count the number of points we want to receive from every PE
   ! and the total number of points to output

   icnt(:) = 0
   flag(:) = 0

   p_pat%n_pnts = 0

   DO i = 1, n_points
      IF(owner(i)>=0) THEN
         p_pat%n_pnts = p_pat%n_pnts + 1 ! Count total number of points we output
         IF(flag(ABS(global_index(i)))==0) THEN
            icnt(owner(i)) = icnt(owner(i))+1 ! Number to get from owner(i)
            flag(ABS(global_index(i))) = 1 ! Flag that this global point is already on the list
         ENDIF
      ENDIF
   ENDDO

   ! Allocate and set up the recv_limits array

   ALLOCATE(p_pat%recv_limits(0:p_n_work))

   p_pat%recv_limits(0) = 0
   DO np = 0, p_n_work-1
      p_pat%recv_limits(np+1) = p_pat%recv_limits(np) + icnt(np)
   ENDDO

   ! The last entry in recv_limits is the total number of points we receive

   p_pat%n_recv = p_pat%recv_limits(p_n_work)

   ! Allocate and set up the recv_src array

   ALLOCATE(p_pat%recv_src(p_pat%n_pnts))
   ALLOCATE(p_pat%recv_dst_blk(p_pat%n_pnts))
   ALLOCATE(p_pat%recv_dst_idx(p_pat%n_pnts))
   ALLOCATE(global_recv_index(p_pat%n_recv))

   DO np = 0, p_n_work-1
      icnt(np) = p_pat%recv_limits(np)
   ENDDO

   flag(:) = 0
   n = 0 ! Counts total number of local points

   DO i = 1, n_points
      IF(owner(i)>=0) THEN
         n = n+1
         IF(flag(ABS(global_index(i)))==0) THEN
            icnt(owner(i)) = icnt(owner(i)) + 1    ! Current index in recv array
            global_recv_index(icnt(owner(i))) = ABS(global_index(i))
                                                   ! Global index of points in receive array
            p_pat%recv_src(n) = icnt(owner(i))     ! From where in the receive array we get
                                                   ! the local point
            p_pat%recv_dst_blk(n) = blk_no(i)      ! Where to put the local point
            p_pat%recv_dst_idx(n) = idx_no(i)      ! Where to put the local point
            flag(ABS(global_index(i))) = icnt(owner(i)) ! Store from where to get duplicates
         ELSE
            p_pat%recv_src(n) = flag(ABS(global_index(i)))
            p_pat%recv_dst_blk(n) = blk_no(i)
            p_pat%recv_dst_idx(n) = idx_no(i)
         ENDIF
      ENDIF
   ENDDO


   ! Exchange the number of points we want to receive with the respective senders

   DO np = 0, p_n_work-1 ! loop over PEs where to send the data

      num_recv = p_pat%recv_limits(np+1) - p_pat%recv_limits(np)
      irs = p_pat%recv_limits(np)+1 ! Start index in global_recv_index
      ire = p_pat%recv_limits(np+1) ! End   index in global_recv_index

      IF(np/=p_pe_work) THEN
         ! Just send the global index of the points we need from PE np
         CALL p_send(num_recv, np, 1, comm=p_comm_work)
         if(num_recv>0) &
           CALL p_send(global_recv_index(irs:ire), np, 2, comm=p_comm_work)
      ELSE

         ! In this turn, we receive the info about which points are needed from us

         ! First get the number of points only
         DO nr = 0, p_n_work-1
            IF(nr /= p_pe_work) THEN
               CALL p_recv(icnt(nr), nr, 1,  comm=p_comm_work)
            ELSE
               icnt(nr) = num_recv
            ENDIF
         ENDDO

         ! Allocate and set up the send_limits array

         ALLOCATE(p_pat%send_limits(0:p_n_work))

         p_pat%send_limits(0) = 0
         DO nr = 0, p_n_work-1
            p_pat%send_limits(nr+1) = p_pat%send_limits(nr) + icnt(nr)
         ENDDO

         ! The last entry in send_limits is the total number of points we receive

         p_pat%n_send = p_pat%send_limits(p_n_work)

         ! Allocate and set up the send_src array

         ALLOCATE(send_src(p_pat%n_send))

         DO nr = 0, p_n_work-1
            num_send = p_pat%send_limits(nr+1) - p_pat%send_limits(nr)
            iss = p_pat%send_limits(nr)+1 ! Start index in send_src
            ise = p_pat%send_limits(nr+1) ! End   index in send_src
            IF(nr /= p_pe_work) THEN
               IF(num_send>0) &
                  CALL p_recv(send_src(iss:ise), nr, 2, comm=p_comm_work)
            ELSE
               IF(num_send>0) &
                  send_src(iss:ise) = global_recv_index(irs:ire)
            ENDIF
         ENDDO

      ENDIF

   ENDDO

   ALLOCATE(p_pat%send_src_blk(p_pat%n_send))
   ALLOCATE(p_pat%send_src_idx(p_pat%n_send))

   ! The indices in p_pat%send_src are global, convert ot local

   DO i = 1, p_pat%n_send
      np = send_src(i)
      if(np<1 .or. np>SIZE(local_index)) CALL finish('setup_comm_pattern','Got illegal index')
      np = local_index(np)
      if(np<0) CALL finish('setup_comm_pattern','Got illegal index')
      p_pat%send_src_blk(i) = blk_no(np)
      p_pat%send_src_idx(i) = idx_no(np)
   ENDDO

   ! Finally, compute lists of processors for send and receive operations

   num_send = 0
   num_recv = 0

   DO np = 0, p_n_work-1 ! loop over PEs

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

   DO np = 0, p_n_work-1 ! loop over PEs

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

END SUBROUTINE setup_comm_pattern
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

   TYPE(t_comm_pattern), INTENT(INOUT) :: p_pat

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

   p_pat%n_recv = 0
   p_pat%n_pnts = 0
   p_pat%n_send = 0
   p_pat%np_recv = 0
   p_pat%np_send = 0

END SUBROUTINE delete_comm_pattern

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
SUBROUTINE exchange_data_r3d(p_pat, recv, send, add, send_lbound3)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   REAL(wp), INTENT(INOUT), TARGET        :: recv(:,:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
   INTEGER, OPTIONAL :: send_lbound3

   REAL(wp) :: send_buf(SIZE(recv,2),p_pat%n_send), &
               recv_buf(SIZE(recv,2),p_pat%n_recv)

   REAL(wp), POINTER :: send_ptr(:,:,:)

   INTEGER :: i, k, np, irs, iss, pid, icount, ndim2, lbound3

   !-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_exch_data)

   IF(my_process_is_mpi_seq()) &
     CALL finish('exchange_data','must not be called on single PE/test PE')

   IF(SIZE(recv,1) /= nproma) THEN
     CALL finish('exchange_data','Illegal first dimension of data array')
   ENDIF

   ndim2 = SIZE(recv,2)

   IF((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) .AND. .NOT. use_exchange_delayed) THEN
     ! Set up irecv's for receive buffers
     DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_recv(np) ! ID of receiver PE
       irs    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2
       CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF

   IF(PRESENT(send) .AND. PRESENT(send_lbound3)) THEN
     lbound3 = send_lbound3
   ELSE
     lbound3 = 1
   ENDIF

   ! Set up send buffer

   IF(PRESENT(send)) THEN
     send_ptr => send
   ELSE
     send_ptr => recv
   ENDIF

   IF (ndim2 == 1) THEN
     DO i = 1, p_pat%n_send
       send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i)-lbound3+1)
     ENDDO
   ELSE
#ifdef __SX__
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_send
         send_buf(k,i) = send_ptr(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)-lbound3+1)
       ENDDO
     ENDDO
#else
     DO i = 1, p_pat%n_send
       send_buf(1:ndim2,i) = send_ptr(p_pat%send_src_idx(i),1:ndim2,   &
         &                            p_pat%send_src_blk(i)-lbound3+1)
     ENDDO
#endif
   ENDIF

   IF(use_exchange_delayed) THEN
     CALL add_delayed_request(p_pat, ndim2, .FALSE., p_real_dp)
     delayed_request_r(n_delayed_requests_r)%recv3 => recv
     IF(PRESENT(add)) delayed_request_r(n_delayed_requests_r)%add3 => add

     CALL buffer_data(p_pat%send_limits, send_buf)
     RETURN
   ENDIF


   ! Send our data
   IF (iorder_sendrecv == 1) THEN
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2
       CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO

     DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_recv(np) ! ID of receiver PE
       irs    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2
       CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 3) THEN ! use irecv/isend
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

   IF(PRESENT(add)) THEN
     IF (ndim2 == 1) THEN
       k = 1
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
           recv_buf(k,p_pat%recv_src(i)) +                       &
           add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
       ENDDO
     ELSE
#ifdef __SX__
!CDIR UNROLL=6
       DO k = 1, ndim2
         DO i = 1, p_pat%n_pnts
           recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
             recv_buf(k,p_pat%recv_src(i)) +                       &
             add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
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
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
           recv_buf(k,p_pat%recv_src(i))
       ENDDO
     ELSE
#ifdef __SX__
!CDIR UNROLL=6
       DO k = 1, ndim2
         DO i = 1, p_pat%n_pnts
           recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
             recv_buf(k,p_pat%recv_src(i))
         ENDDO
       ENDDO
#else
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) = &
           recv_buf(:,p_pat%recv_src(i))
       ENDDO
#endif
     ENDIF
   ENDIF

#ifdef __BOUNDCHECK
   DO i = 1, p_pat%n_pnts
     IF (p_pat%recv_dst_idx(i)>nproma .OR. p_pat%recv_dst_blk(i)>SIZE(recv,3)) &
       WRITE(0,*)'exch_r3d',i,p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i),SIZE(recv,3)
   ENDDO
#endif

   IF (activate_sync_timers) CALL timer_stop(timer_exch_data)
   
END SUBROUTINE exchange_data_r3d
!================================================================================================
! INTEGER SECTION -------------------------------------------------------------------------------
! 
SUBROUTINE exchange_data_i3d(p_pat, recv, send, add, send_lbound3)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   INTEGER, INTENT(INOUT), TARGET        :: recv(:,:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
   INTEGER, OPTIONAL :: send_lbound3

   INTEGER :: send_buf(SIZE(recv,2),p_pat%n_send), &
               recv_buf(SIZE(recv,2),p_pat%n_recv)

   INTEGER, POINTER :: send_ptr(:,:,:)

   INTEGER :: i, k, np, irs, iss, pid, icount, ndim2, lbound3

   !-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_exch_data)

   IF(my_process_is_mpi_seq()) &
     CALL finish('exchange_data','must not be called on single PE/test PE')

   IF(SIZE(recv,1) /= nproma) THEN
     CALL finish('exchange_data','Illegal first dimension of data array')
   ENDIF

   ndim2 = SIZE(recv,2)

   IF((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) .AND. .NOT. use_exchange_delayed) THEN
     ! Set up irecv's for receive buffers
     DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_recv(np) ! ID of receiver PE
       irs    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2
       CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF

   IF(PRESENT(send) .AND. PRESENT(send_lbound3)) THEN
     lbound3 = send_lbound3
   ELSE
     lbound3 = 1
   ENDIF

   ! Set up send buffer

   IF(PRESENT(send)) THEN
     send_ptr => send
   ELSE
     send_ptr => recv
   ENDIF

   IF (ndim2 == 1) THEN
     DO i = 1, p_pat%n_send
       send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i)-lbound3+1)
     ENDDO
   ELSE
#ifdef __SX__
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_send
         send_buf(k,i) = send_ptr(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)-lbound3+1)
       ENDDO
     ENDDO
#else
     DO i = 1, p_pat%n_send
       send_buf(1:ndim2,i) = send_ptr(p_pat%send_src_idx(i),1:ndim2,   &
         &                            p_pat%send_src_blk(i)-lbound3+1)
     ENDDO
#endif
   ENDIF

   IF(use_exchange_delayed) THEN
     CALL add_delayed_request(p_pat, ndim2, .FALSE., p_int)
     delayed_request_i(n_delayed_requests_i)%recv3 => recv
     IF(PRESENT(add)) delayed_request_i(n_delayed_requests_i)%add3 => add

     CALL buffer_data(p_pat%send_limits, send_buf)
     RETURN
   ENDIF


   ! Send our data
   IF (iorder_sendrecv == 1) THEN
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2
       CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO

     DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_recv(np) ! ID of receiver PE
       irs    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2
       CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 3) THEN ! use irecv/isend
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

   IF(PRESENT(add)) THEN
     IF (ndim2 == 1) THEN
       k = 1
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
           recv_buf(k,p_pat%recv_src(i)) +                       &
           add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
       ENDDO
     ELSE
#ifdef __SX__
!CDIR UNROLL=6
       DO k = 1, ndim2
         DO i = 1, p_pat%n_pnts
           recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
             recv_buf(k,p_pat%recv_src(i)) +                       &
             add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
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
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
           recv_buf(k,p_pat%recv_src(i))
       ENDDO
     ELSE
#ifdef __SX__
!CDIR UNROLL=6
       DO k = 1, ndim2
         DO i = 1, p_pat%n_pnts
           recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
             recv_buf(k,p_pat%recv_src(i))
         ENDDO
       ENDDO
#else
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) = &
           recv_buf(:,p_pat%recv_src(i))
       ENDDO
#endif
     ENDIF
   ENDIF

#ifdef __BOUNDCHECK
   DO i = 1, p_pat%n_pnts
     IF (p_pat%recv_dst_idx(i)>nproma .OR. p_pat%recv_dst_blk(i)>SIZE(recv,3)) &
       WRITE(0,*)'exch_i3d',i,p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i),SIZE(recv,3)
   ENDDO
#endif
   IF (activate_sync_timers) CALL timer_stop(timer_exch_data)

END SUBROUTINE exchange_data_i3d
!================================================================================================
! LOGICAL SECTION -------------------------------------------------------------------------------
! 
SUBROUTINE exchange_data_l3d(p_pat, recv, send, send_lbound3)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   LOGICAL, INTENT(INOUT), TARGET        :: recv(:,:,:)
   LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
   INTEGER, OPTIONAL :: send_lbound3

   LOGICAL :: send_buf(SIZE(recv,2),p_pat%n_send), &
               recv_buf(SIZE(recv,2),p_pat%n_recv)

   LOGICAL, POINTER :: send_ptr(:,:,:)

   INTEGER :: i, k, np, irs, ire, iss, ise, ndim2, lbound3

   !-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_exch_data)

   IF(my_process_is_mpi_seq()) &
     CALL finish('exchange_data','must not be called on single PE/test PE')

   IF(SIZE(recv,1) /= nproma) THEN
     CALL finish('exchange_data','Illegal first dimension of data array')
   ENDIF

   ndim2 = SIZE(recv,2)

   IF((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) .AND. .NOT. use_exchange_delayed) THEN
     ! Set up irecv's for receive buffers
     DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

       irs = p_pat%recv_limits(np)+1
       ire = p_pat%recv_limits(np+1)
       IF(ire >= irs) &
         CALL p_irecv(recv_buf(1,irs), np, 1, p_count=(ire-irs+1)*ndim2, comm=p_comm_work)

     ENDDO
   ENDIF

   IF(PRESENT(send) .AND. PRESENT(send_lbound3)) THEN
     lbound3 = send_lbound3
   ELSE
     lbound3 = 1
   ENDIF

   ! Set up send buffer

   IF(PRESENT(send)) THEN
     send_ptr => send
   ELSE
     send_ptr => recv
   ENDIF

   IF (ndim2 == 1) THEN
     DO i = 1, p_pat%n_send
       send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i)-lbound3+1)
     ENDDO
   ELSE
#ifdef __SX__
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_send
         send_buf(k,i) = send_ptr(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)-lbound3+1)
       ENDDO
     ENDDO
#else
     DO i = 1, p_pat%n_send
       send_buf(1:ndim2,i) = send_ptr(p_pat%send_src_idx(i),1:ndim2,   &
         &                            p_pat%send_src_blk(i)-lbound3+1)
     ENDDO
#endif
   ENDIF

   IF(use_exchange_delayed) THEN
     CALL add_delayed_request(p_pat, ndim2, .FALSE., p_bool)
     delayed_request_l(n_delayed_requests_l)%recv3 => recv

     CALL buffer_data(p_pat%send_limits, send_buf)
     RETURN
   ENDIF


   ! Send our data
   IF (iorder_sendrecv == 1) THEN
     DO np = 0, p_n_work-1 ! loop over PEs where to send the data

       iss = p_pat%send_limits(np)+1
       ise = p_pat%send_limits(np+1)

       IF(ise >= iss) &
         CALL p_send(send_buf(1,iss), np, 1, p_count=(ise-iss+1)*ndim2, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
     DO np = 0, p_n_work-1 ! loop over PEs where to send the data

       iss = p_pat%send_limits(np)+1
       ise = p_pat%send_limits(np+1)

       IF(ise >= iss) &
         CALL p_isend(send_buf(1,iss), np, 1, p_count=(ise-iss+1)*ndim2, comm=p_comm_work)

     ENDDO

     DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

       irs = p_pat%recv_limits(np)+1
       ire = p_pat%recv_limits(np+1)
       IF(ire >= irs) &
         CALL p_recv(recv_buf(1,irs), np, 1, p_count=(ire-irs+1)*ndim2, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 3) THEN ! use irecv/isend
     DO np = 0, p_n_work-1 ! loop over PEs where to send the data

       iss = p_pat%send_limits(np)+1
       ise = p_pat%send_limits(np+1)

       IF(ise >= iss) &
         CALL p_isend(send_buf(1,iss), np, 1, p_count=(ise-iss+1)*ndim2, comm=p_comm_work)

     ENDDO
   ENDIF

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

   IF (ndim2 == 1) THEN
     k = 1
     DO i = 1, p_pat%n_pnts
       recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
            recv_buf(k,p_pat%recv_src(i))
     ENDDO
   ELSE
#ifdef __SX__
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
              recv_buf(k,p_pat%recv_src(i))
       ENDDO
     ENDDO
#else
     DO i = 1, p_pat%n_pnts
       recv(p_pat%recv_dst_idx(i),:,p_pat%recv_dst_blk(i)) = &
            recv_buf(:,p_pat%recv_src(i))
     ENDDO
#endif
   ENDIF

#ifdef __BOUNDCHECK
   DO i = 1, p_pat%n_pnts
     IF (p_pat%recv_dst_idx(i)>nproma .OR. p_pat%recv_dst_blk(i)>SIZE(recv,3)) &
       WRITE(0,*)'exch_l3d',i,p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i),SIZE(recv,3)
   ENDDO
#endif
   IF (activate_sync_timers) CALL timer_stop(timer_exch_data)

END SUBROUTINE exchange_data_l3d
!
!================================================================================================
!
!--------------------------------------------------------------------------------------------------
!>
!! Does a reverse data exchange according to a communication pattern (in p_pat),
!! i.e. if exchange_data transfers data from source point x to destination point y
!! this routine transfers from data from y to x.
!!
!! Please note:
!! The communication pattern used should have exactly one destination point
!! for every source point, otherwise it is undefined which of the different
!! destination data points is transferred to the source.
!!
!! !REVISION HISTORY:
!! Initial version by Rainer Johanni, Dec 2009
!!
SUBROUTINE exchange_data_reverse_3(p_pat, recv, send)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   REAL(wp), INTENT(INOUT), TARGET        :: recv(:,:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)

! !LOCAL VARIABLES:

   ! Note reversed second dimensions
   REAL(wp) :: send_buf(SIZE(recv,2),p_pat%n_recv), &
               recv_buf(SIZE(recv,2),p_pat%n_send)

   REAL(wp), POINTER :: send_ptr(:,:,:)

   INTEGER :: i, k, jb, jl, np, irs, iss, pid, icount, ndim2

!EOP
!-----------------------------------------------------------------------
!BOC
   IF (activate_sync_timers) CALL timer_start(timer_exch_data_rv)

   IF(my_process_is_mpi_seq()) &
      CALL finish('exchange_data_reverse','must not be called on single PE/test PE')

   IF(SIZE(recv,1) /= nproma) THEN
      CALL finish('exchange_data_reverse','Illegal first dimension of data array')
   ENDIF

   ndim2 = SIZE(recv,2)

   IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) .AND. .NOT. use_exchange_delayed) THEN
     ! Set up irecv's for receive buffers
     DO np = 1, p_pat%np_send ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_send(np) ! ID of receiver PE
       irs    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2
       CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF

   ! Set up send buffer

   IF(PRESENT(send)) THEN
     send_ptr => send
   ELSE
     send_ptr => recv
   ENDIF

   IF (ndim2 == 1) THEN
     DO i = 1, p_pat%n_pnts
       send_buf(1,p_pat%recv_src(i)) = send_ptr(p_pat%recv_dst_idx(i),1,p_pat%recv_dst_blk(i))
     ENDDO
   ELSE
#ifndef __SX__
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_pnts
         send_buf(k,p_pat%recv_src(i)) = send_ptr(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
       ENDDO
     ENDDO
#else
     DO i = 1, p_pat%n_pnts
       send_buf(1:ndim2,p_pat%recv_src(i)) = &
         & send_ptr(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i))
     ENDDO
#endif
   ENDIF

   IF(use_exchange_delayed) THEN
     CALL add_delayed_request(p_pat, ndim2, .TRUE., p_real_dp)
     delayed_request_r(n_delayed_requests_r)%recv3 => recv

     CALL buffer_data(p_pat%recv_limits, send_buf)
     RETURN
   ENDIF


   ! Send our data
   IF (iorder_sendrecv == 1) THEN
     DO np = 1, p_pat%np_recv ! loop over PEs where to send the data

       pid    = p_pat%pelist_recv(np) ! ID of sender PE
       iss    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2
       CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
     DO np = 1, p_pat%np_recv ! loop over PEs where to send the data

       pid    = p_pat%pelist_recv(np) ! ID of sender PE
       iss    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO

     DO np = 1, p_pat%np_send ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_send(np) ! ID of receiver PE
       irs    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2
       CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 3) THEN ! use isend/irecv
     DO np = 1, p_pat%np_recv ! loop over PEs where to send the data

       pid    = p_pat%pelist_recv(np) ! ID of sender PE
       iss    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

#ifndef __SX__
   DO i = 1, p_pat%n_send
      jb = p_pat%send_src_blk(i)
      jl = p_pat%send_src_idx(i)
      recv(jl,:,jb) = recv_buf(:,i)
   ENDDO
#else
   IF (ndim2 == 1) THEN
     k = 1
     DO i = 1, p_pat%n_send
       recv(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)) = recv_buf(k,i)
     ENDDO
   ELSE
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_send
         recv(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)) = recv_buf(k,i)
       ENDDO
     ENDDO
   ENDIF
#endif

#ifdef __BOUNDCHECK
   DO i = 1, p_pat%n_send
     IF (p_pat%send_src_idx(i)>nproma .OR.p_pat%send_src_blk(i)>SIZE(recv,3)) &
       WRITE(0,*)'exch_3r',i,p_pat%send_src_idx(i),p_pat%send_src_blk(i),SIZE(recv,3)
   ENDDO
#endif
   IF (activate_sync_timers) CALL timer_stop(timer_exch_data_rv)

END SUBROUTINE exchange_data_reverse_3

!-------------------------------------------------------------------------
!
!

!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process 4D fields or up to seven 3D fields
!! in one step
!!
SUBROUTINE exchange_data_mult(p_pat, nfields, ndim2tot, recv1, send1, add1, recv2, send2,   &
                              add2, recv3, send3, add3, recv4, send4, add4, recv5, send5,   &
                              add5, recv6, send6, add6, recv7, send7, add7, recv4d, send4d, &
                              add4d, lpar, nshift)

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat

   REAL(wp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), &
     recv7(:,:,:), recv4d(:,:,:,:)
   REAL(wp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     send1(:,:,:), send2(:,:,:), send3(:,:,:), send4(:,:,:), send5(:,:,:), send6(:,:,:), &
     send7(:,:,:), send4d(:,:,:,:)
   REAL(wp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     add1(:,:,:), add2(:,:,:), add3(:,:,:), add4(:,:,:), add5(:,:,:), add6(:,:,:),       &
     add7(:,:,:), add4d(:,:,:,:)

   INTEGER, INTENT(IN)           :: nfields, ndim2tot
   LOGICAL, OPTIONAL, INTENT(IN) :: lpar
   INTEGER, OPTIONAL, INTENT(IN) :: nshift

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   TYPE(t_fieldptr) :: recv(nfields), send(nfields), add(nfields)
   INTEGER        :: ndim2(nfields), noffset(nfields)

   REAL(wp) :: send_buf(ndim2tot,p_pat%n_send),recv_buf(ndim2tot,p_pat%n_recv)

   INTEGER :: i, k, kshift, jb,ik, jl, n, np, irs, iss, pid, icount
   LOGICAL :: lsend, ladd, l_par

!-----------------------------------------------------------------------


   lsend     = .FALSE.
   ladd      = .FALSE.

   IF (PRESENT(nshift)) THEN
     kshift = nshift
   ELSE
     kshift = 0
   ENDIF

   IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) .AND. .NOT. use_exchange_delayed) THEN
     ! Set up irecv's for receive buffers
     DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_recv(np) ! ID of receiver PE
       irs    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2tot
       CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF

   ! Set pointers to input fields
   IF (PRESENT(recv4d)) THEN
     DO i = 1, nfields
       recv(i)%fld => recv4d(:,:,:,i)
     ENDDO
     IF (PRESENT(send4d)) THEN ! all 4D fields must have the same dimensions
       DO i = 1, nfields
         send(i)%fld => send4d(:,:,:,i)
       ENDDO
       lsend = .TRUE.
     ENDIF
     IF (PRESENT(add4d)) THEN
       DO i = 1, nfields
         add(i)%fld => add4d(:,:,:,i)
       ENDDO
       ladd = .TRUE.
     ENDIF
   ELSE
     IF (PRESENT(recv1)) THEN
       recv(1)%fld => recv1
       IF (PRESENT(send1)) THEN
         send(1)%fld => send1
         lsend = .TRUE.
       ENDIF
       IF (PRESENT(add1)) THEN
         add(1)%fld => add1
         ladd = .TRUE.
       ENDIF
       IF (PRESENT(recv2)) THEN
         recv(2)%fld => recv2
         IF (lsend) send(2)%fld => send2
         IF (ladd)  add(2)%fld  => add2
         IF (PRESENT(recv3)) THEN
           recv(3)%fld => recv3
           IF (lsend) send(3)%fld => send3
           IF (ladd)  add(3)%fld  => add3
           IF (PRESENT(recv4)) THEN
             recv(4)%fld => recv4
             IF (lsend) send(4)%fld => send4
             IF (ladd)  add(4)%fld  => add4
             IF (PRESENT(recv5)) THEN
               recv(5)%fld => recv5
               IF (lsend) send(5)%fld => send5
               IF (ladd)  add(5)%fld  => add5
               IF (PRESENT(recv6)) THEN
                 recv(6)%fld => recv6
                 IF (lsend) send(6)%fld => send6
                 IF (ladd)  add(6)%fld  => add6
                 IF (PRESENT(recv7)) THEN
                   recv(7)%fld => recv7
                   IF (lsend) send(7)%fld => send7
                   IF (ladd)  add(7)%fld  => add7
                 ENDIF
               ENDIF
             ENDIF
           ENDIF
         ENDIF
       ENDIF
     ENDIF
   ENDIF

   IF(use_exchange_delayed) THEN
     DO n = 1, nfields
       IF(lsend) THEN
         IF(ladd) THEN
           CALL exchange_data(p_pat, recv(n)%fld(:,:,:), send(n)%fld(:,:,:), add(n)%fld(:,:,:))
         ELSE
           CALL exchange_data(p_pat, recv(n)%fld(:,:,:), send(n)%fld(:,:,:))
         ENDIF
       ELSE
         IF(ladd) THEN
           CALL exchange_data(p_pat, recv(n)%fld(:,:,:), add=add(n)%fld(:,:,:))
         ELSE
           CALL exchange_data(p_pat, recv(n)%fld(:,:,:))
         ENDIF
       ENDIF
     ENDDO
     RETURN
   ENDIF

   IF (PRESENT(lpar)) THEN
     l_par = lpar
   ELSE
     l_par = .FALSE.
   ENDIF

   noffset(1) = 0
   ndim2(1)   = SIZE(recv(1)%fld,2) - kshift
   DO n = 2, nfields
     noffset(n) = noffset(n-1)+ndim2(n-1)
     ndim2(n)   = SIZE(recv(n)%fld,2) - kshift
   ENDDO

   ! Set up send buffer
#ifdef __SX__
   IF ( lsend .AND. l_par ) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(n,k,i)
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_send
           send_buf(k+noffset(n),i) = &
             send(n)%fld(p_pat%send_src_idx(i),k+kshift,p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
!$OMP END DO
!$OMP END PARALLEL
   ELSE IF ( lsend ) THEN
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_send
           send_buf(k+noffset(n),i) = &
             send(n)%fld(p_pat%send_src_idx(i),k+kshift,p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ELSE
       ! Send and receive arrays are identical (for boundary exchange)
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_send
           send_buf(k+noffset(n),i) = &
             recv(n)%fld(p_pat%send_src_idx(i),k+kshift,p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ENDIF
#else
   DO i = 1, p_pat%n_send
     jb = p_pat%send_src_blk(i)
     jl = p_pat%send_src_idx(i)
     IF ( lsend ) THEN
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),i) = send(n)%fld(jl,k+kshift,jb)
         ENDDO
       ENDDO
     ELSE
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),i) = recv(n)%fld(jl,k+kshift,jb)
         ENDDO
       ENDDO
     ENDIF
   ENDDO
#endif


   ! Send our data
   IF (iorder_sendrecv == 1) THEN
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2tot
       CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2tot
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO

     DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_recv(np) ! ID of receiver PE
       irs    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2tot
       CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 3) THEN ! use isend/irecv
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2tot
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

#ifdef __SX__
   IF (ladd .AND. l_par) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(n,k,i)
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           recv(n)%fld(p_pat%recv_dst_idx(i),k+kshift,p_pat%recv_dst_blk(i)) =  &
             recv_buf(k+noffset(n),p_pat%recv_src(i)) +                         &
             add(n)%fld(p_pat%recv_dst_idx(i),k+kshift,p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
     ENDDO
!$OMP END DO
!$OMP END PARALLEL
   ELSE IF (ladd) THEN
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           recv(n)%fld(p_pat%recv_dst_idx(i),k+kshift,p_pat%recv_dst_blk(i)) =  &
             recv_buf(k+noffset(n),p_pat%recv_src(i)) +                         &
             add(n)%fld(p_pat%recv_dst_idx(i),k+kshift,p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ELSE
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           recv(n)%fld(p_pat%recv_dst_idx(i),k+kshift,p_pat%recv_dst_blk(i)) =  &
             recv_buf(k+noffset(n),p_pat%recv_src(i))
         ENDDO
       ENDDO
     ENDDO
   ENDIF
#else
   DO i = 1, p_pat%n_pnts
     jb = p_pat%recv_dst_blk(i)
     jl = p_pat%recv_dst_idx(i)
     ik  = p_pat%recv_src(i)
     IF (ladd) THEN
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           recv(n)%fld(jl,k+kshift,jb) = recv_buf(k+noffset(n),ik)+add(n)%fld(jl,k+kshift,jb)
         ENDDO
       ENDDO
     ELSE
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           recv(n)%fld(jl,k+kshift,jb) = recv_buf(k+noffset(n),ik)
         ENDDO
       ENDDO
     ENDIF
   ENDDO
#endif

#ifdef __BOUNDCHECK
   DO n = 1, nfields
     DO i = 1, p_pat%n_pnts
       IF (p_pat%recv_dst_idx(i)>nproma .OR.p_pat%recv_dst_blk(i)>SIZE(recv(n)%fld,3)) &
         WRITE(0,*)'exch_mult',n,i,p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i),SIZE(recv(n)%fld,3)
     ENDDO
   ENDDO
#endif
     
END SUBROUTINE exchange_data_mult

!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process a 4D field whose extra dimension
!! is on the third index
!!
SUBROUTINE exchange_data_4de3(p_pat, nfields, ndim2tot, recv, send)

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat

   REAL(wp), INTENT(INOUT)           :: recv(:,:,:,:)
   REAL(wp), INTENT(IN   ), OPTIONAL :: send(:,:,:,:)

   INTEGER, INTENT(IN)           :: nfields, ndim2tot

   INTEGER :: ndim2, noffset

   REAL(wp) :: send_buf(ndim2tot,p_pat%n_send),recv_buf(ndim2tot,p_pat%n_recv)

   INTEGER :: i, k, ik, jb, jl, n, np, irs, iss, pid, icount
   LOGICAL :: lsend

!-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_exch_data)

   IF (PRESENT(send)) THEN
     lsend  = .TRUE.
   ELSE
     lsend  = .FALSE.
   ENDIF

   ndim2 = SIZE(recv,2)

   IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3)) THEN
     ! Set up irecv's for receive buffers
     DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_recv(np) ! ID of receiver PE
       irs    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2tot
       CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF


   ! Set up send buffer
#ifdef __SX__
   IF ( lsend ) THEN
     DO n = 1, nfields
       noffset = (n-1)*ndim2
!CDIR UNROLL=6
       DO k = 1, ndim2
         DO i = 1, p_pat%n_send
           send_buf(k+noffset,i) = &
             send(p_pat%send_src_idx(i),k,n,p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ELSE
       ! Send and receive arrays are identical (for boundary exchange)
     DO n = 1, nfields
       noffset = (n-1)*ndim2
!CDIR UNROLL=6
       DO k = 1, ndim2
         DO i = 1, p_pat%n_send
           send_buf(k+noffset,i) = &
             recv(p_pat%send_src_idx(i),k,n,p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ENDIF
#else
   DO i = 1, p_pat%n_send
     jb = p_pat%send_src_blk(i)
     jl = p_pat%send_src_idx(i)
     IF ( lsend ) THEN
       DO n = 1, nfields
         noffset = (n-1)*ndim2
         DO k = 1, ndim2
           send_buf(k+noffset,i) = send(jl,k,n,jb)
         ENDDO
       ENDDO
     ELSE
       DO n = 1, nfields
         noffset = (n-1)*ndim2
         DO k = 1, ndim2
           send_buf(k+noffset,i) = recv(jl,k,n,jb)
         ENDDO
       ENDDO
     ENDIF
   ENDDO
#endif


   ! Send our data
   IF (iorder_sendrecv == 1) THEN
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2tot
       CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2tot
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO

     DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_recv(np) ! ID of receiver PE
       irs    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2tot
       CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 3) THEN ! use isend/irecv
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2tot
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF

   ! Wait for all outstanding requests to finish

   CALL p_wait

   ! Fill in receive buffer

#ifdef __SX__
   DO n = 1, nfields
     noffset = (n-1)*ndim2
!CDIR UNROLL=6
     DO k = 1, ndim2
       DO i = 1, p_pat%n_pnts
         recv(p_pat%recv_dst_idx(i),k,n,p_pat%recv_dst_blk(i)) =  &
           recv_buf(k+noffset,p_pat%recv_src(i))
       ENDDO
     ENDDO
   ENDDO
#else
   DO i = 1, p_pat%n_pnts
     jb = p_pat%recv_dst_blk(i)
     jl = p_pat%recv_dst_idx(i)
     ik  = p_pat%recv_src(i)
     DO n = 1, nfields
       noffset = (n-1)*ndim2
       DO k = 1, ndim2
         recv(jl,k,n,jb) = recv_buf(k+noffset,ik)
       ENDDO
     ENDDO
   ENDDO
#endif
   IF (activate_sync_timers) CALL timer_stop(timer_exch_data)

END SUBROUTINE exchange_data_4de3


!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process 4D fields or up to seven 3D fields
!! in one step
!! Note: this routine has to be called from within an OpenMP-parallel section
!!
SUBROUTINE exchange_data_gm(p_pat, nfields, ndim2tot, send_buf, recv_buf, recv1, send1, add1, &
                            recv2, send2, add2, recv3, send3, add3, recv4, send4, add4,       &
                            recv5, send5, add5, recv6, send6, add6, recv7, send7, add7,       &
                            recv4d, send4d, add4d)

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat

   REAL(wp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), &
     recv7(:,:,:), recv4d(:,:,:,:)
   REAL(wp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     send1(:,:,:), send2(:,:,:), send3(:,:,:), send4(:,:,:), send5(:,:,:), send6(:,:,:), &
     send7(:,:,:), send4d(:,:,:,:)
   REAL(wp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     add1(:,:,:), add2(:,:,:), add3(:,:,:), add4(:,:,:), add5(:,:,:), add6(:,:,:),       &
     add7(:,:,:), add4d(:,:,:,:)

   INTEGER, INTENT(IN)           :: nfields, ndim2tot

   REAL(wp) , INTENT(INOUT) :: send_buf(:,:),recv_buf(:,:)

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   TYPE(t_fieldptr) :: recv(nfields), send(nfields), add(nfields)
   INTEGER        :: ndim2(nfields), noffset(nfields)

   INTEGER :: i, k, ik, jb, jl, n, np, irs, iss, pid, icount
   LOGICAL :: lsend, ladd

!-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_exch_data)

   lsend     = .FALSE.
   ladd      = .FALSE.

   IF (iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) THEN
!$OMP MASTER
     ! Set up irecv's for receive buffers
     DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_recv(np) ! ID of receiver PE
       irs    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2tot
       CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
!$OMP END MASTER
   ENDIF

   ! Note: there is no OpenMP barrier set here because the subsequent filling of
   ! the send buffer is supposed to be executed in parallel to the irecv calls


   ! Set pointers to input fields
   IF (PRESENT(recv4d)) THEN
     DO i = 1, nfields
       recv(i)%fld => recv4d(:,:,:,i)
     ENDDO
     IF (PRESENT(send4d)) THEN ! all 4D fields must have the same dimensions
       DO i = 1, nfields
         send(i)%fld => send4d(:,:,:,i)
       ENDDO
       lsend = .TRUE.
     ENDIF
     IF (PRESENT(add4d)) THEN
       DO i = 1, nfields
         add(i)%fld => add4d(:,:,:,i)
       ENDDO
       ladd = .TRUE.
     ENDIF
   ELSE
     IF (PRESENT(recv1)) THEN
       recv(1)%fld => recv1
       IF (PRESENT(send1)) THEN
         send(1)%fld => send1
         lsend = .TRUE.
       ENDIF
       IF (PRESENT(add1)) THEN
         add(1)%fld => add1
         ladd = .TRUE.
       ENDIF
       IF (PRESENT(recv2)) THEN
         recv(2)%fld => recv2
         IF (lsend) send(2)%fld => send2
         IF (ladd)  add(2)%fld  => add2
         IF (PRESENT(recv3)) THEN
           recv(3)%fld => recv3
           IF (lsend) send(3)%fld => send3
           IF (ladd)  add(3)%fld  => add3
           IF (PRESENT(recv4)) THEN
             recv(4)%fld => recv4
             IF (lsend) send(4)%fld => send4
             IF (ladd)  add(4)%fld  => add4
             IF (PRESENT(recv5)) THEN
               recv(5)%fld => recv5
               IF (lsend) send(5)%fld => send5
               IF (ladd)  add(5)%fld  => add5
               IF (PRESENT(recv6)) THEN
                 recv(6)%fld => recv6
                 IF (lsend) send(6)%fld => send6
                 IF (ladd)  add(6)%fld  => add6
                 IF (PRESENT(recv7)) THEN
                   recv(7)%fld => recv7
                   IF (lsend) send(7)%fld => send7
                   IF (ladd)  add(7)%fld  => add7
                 ENDIF
               ENDIF
             ENDIF
           ENDIF
         ENDIF
       ENDIF
     ENDIF
   ENDIF


   noffset(1) = 0
   ndim2(1)   = SIZE(recv(1)%fld,2)
   DO n = 2, nfields
     noffset(n) = noffset(n-1)+ndim2(n-1)
     ndim2(n)   = SIZE(recv(n)%fld,2)
   ENDDO

   ! Set up send buffer
#ifdef __SX__
   IF ( lsend ) THEN
!$OMP DO SCHEDULE(guided)
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_send
           send_buf(k+noffset(n),i) = &
             send(n)%fld(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
!$OMP END DO NOWAIT
   ELSE
       ! Send and receive arrays are identical (for boundary exchange)
!$OMP DO SCHEDULE(guided)
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_send
           send_buf(k+noffset(n),i) = &
             recv(n)%fld(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
!$OMP END DO NOWAIT
   ENDIF
#else
!$OMP DO SCHEDULE(guided)
   DO i = 1, p_pat%n_send
     jb = p_pat%send_src_blk(i)
     jl = p_pat%send_src_idx(i)
     IF ( lsend ) THEN
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),i) = send(n)%fld(jl,k,jb)
         ENDDO
       ENDDO
     ELSE
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),i) = recv(n)%fld(jl,k,jb)
         ENDDO
       ENDDO
     ENDIF
   ENDDO
!$OMP END DO NOWAIT
#endif

!$OMP BARRIER

!$OMP MASTER
   ! Send our data
   IF (iorder_sendrecv == 1) THEN
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2tot
       CALL p_send(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2tot
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO

     DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

       pid    = p_pat%pelist_recv(np) ! ID of receiver PE
       irs    = p_pat%recv_startidx(np)
       icount = p_pat%recv_count(np)*ndim2tot
       CALL p_recv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ELSE IF (iorder_sendrecv == 3) THEN ! use isend/irecv
     DO np = 1, p_pat%np_send ! loop over PEs where to send the data

       pid    = p_pat%pelist_send(np) ! ID of sender PE
       iss    = p_pat%send_startidx(np)
       icount = p_pat%send_count(np)*ndim2tot
       CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

     ENDDO
   ENDIF

   ! Wait for all outstanding requests to finish

   CALL p_wait

!$OMP END MASTER
!$OMP BARRIER

   ! Fill in receive buffer

#ifdef __SX__
   IF (ladd) THEN
!$OMP DO
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           recv(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
             recv_buf(k+noffset(n),p_pat%recv_src(i)) + &
             add(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
     ENDDO
!$OMP END DO
   ELSE
!$OMP DO
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           recv(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
             recv_buf(k+noffset(n),p_pat%recv_src(i))
         ENDDO
       ENDDO
     ENDDO
!$OMP END DO
   ENDIF
#else
!$OMP DO
   DO i = 1, p_pat%n_pnts
     jb = p_pat%recv_dst_blk(i)
     jl = p_pat%recv_dst_idx(i)
     ik  = p_pat%recv_src(i)
     IF (ladd) THEN
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           recv(n)%fld(jl,k,jb) = recv_buf(k+noffset(n),ik)+add(n)%fld(jl,k,jb)
         ENDDO
       ENDDO
     ELSE
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           recv(n)%fld(jl,k,jb) = recv_buf(k+noffset(n),ik)
         ENDDO
       ENDDO
     ENDIF
   ENDDO
!$OMP END DO
#endif

#ifdef __BOUNDCHECK
   DO n = 1, nfields
     DO i = 1, p_pat%n_pnts
       IF (p_pat%recv_dst_idx(i)>nproma .OR.p_pat%recv_dst_blk(i)>SIZE(recv(n)%fld,3)) &
         WRITE(0,*)'exch_gm',n,i,p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i),SIZE(recv(n)%fld,3)
     ENDDO
   ENDDO
#endif
   IF (activate_sync_timers) CALL timer_stop(timer_exch_data)

END SUBROUTINE exchange_data_gm


!>
!! Starts asynchronous halo communication by filling the send buffer and issueing
!! non-blocking send and receive calls
!!
!! @par Revision History
!! Initial version by Guenther Zaengl, Feb 2011, based on work by Rainer Johanni
!!
SUBROUTINE start_async_comm(p_pat, nfields, ndim2tot, send_buf, recv_buf, recv1, recv2, recv3, &
                            recv4, recv5, recv4d)

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat

   REAL(wp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv4d(:,:,:,:)

   INTEGER, INTENT(IN)           :: nfields, ndim2tot

   REAL(wp) , INTENT(INOUT) :: send_buf(:,:),recv_buf(:,:)

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   TYPE(t_fieldptr) :: recv(nfields)
   INTEGER        :: ndim2(nfields), noffset(nfields)

   INTEGER :: i, k, jb, jl, n, np, irs, iss, pid, icount

!-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_exch_data_async)

   IF(my_process_is_mpi_seq()) RETURN

   ! Start with non-blocking receive calls
   DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

     pid    = p_pat%pelist_recv(np) ! ID of receiver PE
     irs    = p_pat%recv_startidx(np)
     icount = p_pat%recv_count(np)*ndim2tot
     CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

   ENDDO

   ! Set pointers to input fields
   IF (PRESENT(recv4d)) THEN
     DO i = 1, nfields
       recv(i)%fld => recv4d(:,:,:,i)
     ENDDO
   ELSE
     IF (PRESENT(recv1)) THEN
       recv(1)%fld => recv1
       IF (PRESENT(recv2)) THEN
         recv(2)%fld => recv2
         IF (PRESENT(recv3)) THEN
           recv(3)%fld => recv3
           IF (PRESENT(recv4)) THEN
             recv(4)%fld => recv4
             IF (PRESENT(recv5)) THEN
               recv(5)%fld => recv5
             ENDIF
           ENDIF
         ENDIF
       ENDIF
     ENDIF
   ENDIF

   noffset(1) = 0
   ndim2(1)   = SIZE(recv(1)%fld,2)
   DO n = 2, nfields
     noffset(n) = noffset(n-1)+ndim2(n-1)
     ndim2(n)   = SIZE(recv(n)%fld,2)
   ENDDO

   ! Set up send buffer
#ifdef __SX__
   DO n = 1, nfields
!CDIR UNROLL=6
     DO k = 1, ndim2(n)
       DO i = 1, p_pat%n_send
         send_buf(k+noffset(n),i) = &
           recv(n)%fld(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i))
       ENDDO
     ENDDO
   ENDDO
#else
   DO i = 1, p_pat%n_send
     jb = p_pat%send_src_blk(i)
     jl = p_pat%send_src_idx(i)
     DO n = 1, nfields
       DO k = 1, ndim2(n)
         send_buf(k+noffset(n),i) = recv(n)%fld(jl,k,jb)
       ENDDO
     ENDDO
   ENDDO
#endif

   ! Finally issue the non-blocking send calls
   DO np = 1, p_pat%np_send ! loop over PEs where to send the data

     pid    = p_pat%pelist_send(np) ! ID of sender PE
     iss    = p_pat%send_startidx(np)
     icount = p_pat%send_count(np)*ndim2tot
     CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

   ENDDO
   
   IF (activate_sync_timers) CALL timer_stop(timer_exch_data_async)

END SUBROUTINE start_async_comm



!>
!! Completes asynchronous halo communication by filling the receive buffer after completion
!! of the exchange process
!!
!! @par Revision History
!! Initial version by Guenther Zaengl, Feb 2011, based on work by Rainer Johanni
!!
SUBROUTINE complete_async_comm(p_pat, nfields, ndim2tot, recv_buf, recv1, recv2, recv3, &
                               recv4, recv5, recv4d)

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat

   REAL(wp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv4d(:,:,:,:)

   INTEGER, INTENT(IN)      :: nfields, ndim2tot

   REAL(wp) , INTENT(INOUT) :: recv_buf(:,:)

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   TYPE(t_fieldptr) :: recv(nfields)
   INTEGER        :: ndim2(nfields), noffset(nfields)

   INTEGER :: i, k, ik, jb, jl, n

!-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_exch_data_async)

   IF(my_process_is_mpi_seq()) RETURN

   ! Set pointers to output fields
   IF (PRESENT(recv4d)) THEN
     DO i = 1, nfields
       recv(i)%fld => recv4d(:,:,:,i)
     ENDDO
   ELSE
     IF (PRESENT(recv1)) THEN
       recv(1)%fld => recv1
       IF (PRESENT(recv2)) THEN
         recv(2)%fld => recv2
         IF (PRESENT(recv3)) THEN
           recv(3)%fld => recv3
           IF (PRESENT(recv4)) THEN
             recv(4)%fld => recv4
             IF (PRESENT(recv5)) THEN
               recv(5)%fld => recv5
             ENDIF
           ENDIF
         ENDIF
       ENDIF
     ENDIF
   ENDIF

   noffset(1) = 0
   ndim2(1)   = SIZE(recv(1)%fld,2)
   DO n = 2, nfields
     noffset(n) = noffset(n-1)+ndim2(n-1)
     ndim2(n)   = SIZE(recv(n)%fld,2)
   ENDDO

 
   ! Wait for all outstanding requests (issued in start_async_comm) to finish

   CALL p_wait

   ! Fill in receive buffer

#ifdef __SX__
   DO n = 1, nfields
!CDIR UNROLL=6
     DO k = 1, ndim2(n)
       DO i = 1, p_pat%n_pnts
         recv(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
           recv_buf(k+noffset(n),p_pat%recv_src(i))
       ENDDO
     ENDDO
   ENDDO
#else
   DO i = 1, p_pat%n_pnts
     jb = p_pat%recv_dst_blk(i)
     jl = p_pat%recv_dst_idx(i)
     ik  = p_pat%recv_src(i)
     DO n = 1, nfields
       DO k = 1, ndim2(n)
         recv(n)%fld(jl,k,jb) = recv_buf(k+noffset(n),ik)
       ENDDO
     ENDDO
   ENDDO
#endif

   IF (activate_sync_timers) CALL timer_stop(timer_exch_data_async)

END SUBROUTINE complete_async_comm


!>
!! Does data exchange according to a communication pattern (in p_pat).
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl to process 4D fields or up to three 3D fields
!! for an array-sized communication pattern (as needed for boundary interpolation) in one step
!!
SUBROUTINE exchange_data_grf(p_pat, nfields, ndim2tot, nsendtot, nrecvtot, recv1, send1, &
                             recv2, send2, recv3, send3, recv4d, send4d, send_lbound3)

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat(:)

   REAL(wp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4d(:,:,:,:)
   ! Note: the send fields have one additional dimension in this case because
   ! the fourth index corresponds to the dimension of p_pat
   REAL(wp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     send1(:,:,:,:), send2(:,:,:,:), send3(:,:,:,:), send4d(:,:,:,:,:)

   INTEGER, INTENT(IN)           :: nfields  ! total number of input fields
   INTEGER, INTENT(IN)           :: ndim2tot ! sum of vertical levels of input fields
   INTEGER, INTENT(IN)           :: nsendtot ! total number of send points
   INTEGER, INTENT(IN)           :: nrecvtot ! total number of receive points
   INTEGER, OPTIONAL, INTENT(IN) :: send_lbound3

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   TYPE(t_fieldptr) :: recv(nfields), send(nfields*SIZE(p_pat))
   INTEGER        :: ndim2(nfields), noffset(nfields),            &
                     ioffset_s(SIZE(p_pat)), ioffset_r(SIZE(p_pat))

   REAL(wp) :: send_buf(ndim2tot,nsendtot),recv_buf(ndim2tot,nrecvtot), &
               auxs_buf(ndim2tot,nsendtot),auxr_buf(ndim2tot,nrecvtot)

   INTEGER :: i, k, ik, jb, jl, n, np, irs, ire, iss, ise, &
              npats, isum, ioffset, isum1

!-----------------------------------------------------------------------

    npats = SIZE(p_pat)  ! Number of communication patterns provided on input


   ! Set up irecv's for receive buffers
   IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) .AND. .NOT. use_exchange_delayed) THEN

     ioffset = 0
     DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

       ! Sum up receive points over all communication patterns to be processed
       isum = ioffset
       DO n = 1, npats
         isum = isum + p_pat(n)%recv_limits(np+1) - p_pat(n)%recv_limits(np)
       ENDDO

       IF(isum > ioffset) &
         CALL p_irecv(auxr_buf(1,ioffset+1), np, 1, p_count=(isum-ioffset)*ndim2tot, &
                      comm=p_comm_work)
       ioffset = isum

     ENDDO

   ELSE IF (use_exchange_delayed) THEN! use_exchange_delayed
     IF(.NOT.PRESENT(send_lbound3)) CALL finish('exchange_data_grf','Missing send_lbound3')
     IF (PRESENT(recv4d)) THEN
       DO n = 1, UBOUND(recv4d, 4)
         DO np = 1, npats
           CALL exchange_data(p_pat(np), recv4d(:,:,:,n), send4d(:,:,:,np,n), &
             &                  send_lbound3=send_lbound3)
         ENDDO
       ENDDO
     ENDIF
     IF (PRESENT(recv1)) THEN
       DO np = 1, npats
         CALL exchange_data(p_pat(np), recv1(:,:,:), send1(:,:,:,np),send_lbound3=send_lbound3)
       ENDDO
     ENDIF
     IF (PRESENT(recv2)) THEN
       DO np = 1, npats
         CALL exchange_data(p_pat(np), recv2(:,:,:), send2(:,:,:,np),send_lbound3=send_lbound3)
       ENDDO
     ENDIF
     IF (PRESENT(recv3)) THEN
       DO np = 1, npats
         CALL exchange_data(p_pat(np), recv3(:,:,:), send3(:,:,:,np),send_lbound3=send_lbound3)
       ENDDO
     ENDIF

     RETURN
   ENDIF

   ! Set pointers to input fields
   IF (PRESENT(recv4d)) THEN
     DO n = 1, nfields
       recv(n)%fld => recv4d(:,:,:,n)
       DO np = 1, npats
         send(np+(n-1)*npats)%fld => send4d(:,:,:,np,n)
       ENDDO
     ENDDO
   ELSE
     IF (PRESENT(recv1)) THEN
       recv(1)%fld => recv1
       DO np = 1, npats
         send(np)%fld => send1(:,:,:,np)
       ENDDO
     ENDIF
     IF (PRESENT(recv2)) THEN
       recv(2)%fld => recv2
       DO np = 1, npats
         send(np+npats)%fld => send2(:,:,:,np)
       ENDDO
     ENDIF
     IF (PRESENT(recv3)) THEN
       recv(3)%fld => recv3
       DO np = 1, npats
         send(np+2*npats)%fld => send3(:,:,:,np)
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
     ioffset_r(np) = ioffset_r(np-1) + p_pat(np-1)%n_recv
     ioffset_s(np) = ioffset_s(np-1) + p_pat(np-1)%n_send
   ENDDO

   ! Set up send buffer
#ifdef __SX__
   DO n = 1, nfields
     DO np = 1, npats
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
!CDIR NODEP
         DO i = 1, p_pat(np)%n_send
           send_buf(k+noffset(n),i+ioffset_s(np)) =                                            &
             & send(np+(n-1)*npats)%fld(p_pat(np)%send_src_idx(i),k,p_pat(np)%send_src_blk(i)- &
             & send_lbound3+1)
         ENDDO
       ENDDO
     ENDDO
   ENDDO
#else
   DO np = 1, npats
     DO i = 1, p_pat(np)%n_send
       jb = p_pat(np)%send_src_blk(i)
       jl = p_pat(np)%send_src_idx(i)
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),i+ioffset_s(np)) = &
             send(np+(n-1)*npats)%fld(jl,k,jb-send_lbound3+1)
         ENDDO
       ENDDO
     ENDDO
   ENDDO
#endif

   IF (iorder_sendrecv == 1) THEN
     ! Send our data
     ioffset = 0
     DO np = 0, p_n_work-1 ! loop over PEs where to send the data

       ! Copy send points for all communication patterns into one common send buffer
       isum = ioffset
       DO n = 1, npats
         iss = p_pat(n)%send_limits(np)+1 + ioffset_s(n)
         ise = p_pat(n)%send_limits(np+1) + ioffset_s(n)
         isum1 = ise - iss + 1
         IF (isum1 > 0) THEN
!CDIR COLLAPSE
           auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
           isum = isum+isum1
         ENDIF
       ENDDO

       IF(isum > ioffset) CALL p_send(auxs_buf(1,ioffset+1), np, 1,             &
                               p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)

       ioffset = isum

     ENDDO
   ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
     ioffset = 0
     DO np = 0, p_n_work-1 ! loop over PEs where to send the data

       ! Copy send points for all communication patterns into one common send buffer
       isum = ioffset
       DO n = 1, npats
         iss = p_pat(n)%send_limits(np)+1 + ioffset_s(n)
         ise = p_pat(n)%send_limits(np+1) + ioffset_s(n)
         isum1 = ise - iss + 1
         IF (isum1 > 0) THEN
!CDIR COLLAPSE
           auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
           isum = isum+isum1
         ENDIF
       ENDDO

       IF(isum > ioffset) CALL p_isend(auxs_buf(1,ioffset+1), np, 1,            &
                               p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)

       ioffset = isum

     ENDDO

     ioffset = 0
     DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

       ! Sum up receive points over all communication patterns to be processed
       isum = ioffset
       DO n = 1, npats
         isum = isum + p_pat(n)%recv_limits(np+1) - p_pat(n)%recv_limits(np)
       ENDDO

       IF(isum > ioffset) CALL p_recv(auxr_buf(1,ioffset+1), np, 1,             &
                               p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)
       ioffset = isum

     ENDDO
   ELSE IF (iorder_sendrecv == 3) THEN ! use isend/recv
     ioffset = 0
     DO np = 0, p_n_work-1 ! loop over PEs where to send the data

       ! Copy send points for all communication patterns into one common send buffer
       isum = ioffset
       DO n = 1, npats
         iss = p_pat(n)%send_limits(np)+1 + ioffset_s(n)
         ise = p_pat(n)%send_limits(np+1) + ioffset_s(n)
         isum1 = ise - iss + 1
         IF (isum1 > 0) THEN
!CDIR COLLAPSE
           auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
           isum = isum+isum1
         ENDIF
       ENDDO

       IF(isum > ioffset) CALL p_isend(auxs_buf(1,ioffset+1), np, 1,            &
                               p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)

       ioffset = isum

     ENDDO
   ENDIF

   ! Wait for all outstanding requests to finish
   CALL p_wait

   ! Copy exchanged data back to receive buffer
   ioffset = 0
   DO np = 0, p_n_work-1

     isum = ioffset
     DO n = 1, npats
       irs = p_pat(n)%recv_limits(np)+1 + ioffset_r(n)
       ire = p_pat(n)%recv_limits(np+1) + ioffset_r(n)
       isum1 = ire - irs + 1
       IF (isum1 > 0) THEN
!CDIR COLLAPSE
         recv_buf(:,irs:ire) = auxr_buf(:,isum+1:isum+isum1)
         isum = isum + isum1
       ENDIF
     ENDDO

     ioffset = isum

   ENDDO

   ! Fill in receive buffer

#ifdef __SX__
   DO n = 1, nfields
     DO np = 1, npats
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat(np)%n_pnts
           recv(n)%fld(p_pat(np)%recv_dst_idx(i),k,p_pat(np)%recv_dst_blk(i)) =   &
             recv_buf(k+noffset(n),p_pat(np)%recv_src(i)+ioffset_r(np))
         ENDDO
       ENDDO
     ENDDO
   ENDDO
#else
   DO np = 1, npats
     DO i = 1, p_pat(np)%n_pnts
       jb = p_pat(np)%recv_dst_blk(i)
       jl = p_pat(np)%recv_dst_idx(i)
       ik  = p_pat(np)%recv_src(i)+ioffset_r(np)
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           recv(n)%fld(jl,k,jb) = recv_buf(k+noffset(n),ik)
         ENDDO
       ENDDO
     ENDDO
   ENDDO
#endif

#ifdef __BOUNDCHECK
   DO n = 1, nfields
     DO np = 1, npats
       DO i = 1, p_pat(np)%n_pnts
         IF (p_pat(np)%recv_dst_idx(i)>nproma .OR.p_pat(np)%recv_dst_blk(i)>SIZE(recv(n)%fld,3)) &
           WRITE(0,*)'exch_grf',n,np,i,p_pat(np)%recv_dst_idx(i),p_pat(np)%recv_dst_blk(i),     &
           SIZE(recv(n)%fld,3)
       ENDDO
     ENDDO
   ENDDO
#endif

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
SUBROUTINE exchange_data_r2d(p_pat, recv, send, add, send_lbound2, l_recv_exists)
   !
   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   REAL(wp), INTENT(INOUT), TARGET        :: recv(:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
   INTEGER, OPTIONAL :: send_lbound2
   LOGICAL, OPTIONAL :: l_recv_exists

   REAL(wp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))
   REAL(wp) :: send_buf(1,p_pat%n_send)
   INTEGER :: i, lbound2

   !-----------------------------------------------------------------------

   IF (use_exchange_delayed) THEN

     CALL add_delayed_request(p_pat, 1, .FALSE., p_real_dp)

     delayed_request_r(n_delayed_requests_r)%recv2 => recv
     IF (PRESENT(add)) delayed_request_r(n_delayed_requests_r)%add2 => add

     IF (PRESENT(send) .AND. PRESENT(send_lbound2)) THEN
       lbound2 = send_lbound2
     ELSE
       lbound2 = 1
     ENDIF

     IF ( PRESENT(send)) THEN
       DO i = 1, p_pat%n_send
         send_buf(1,i) = send(p_pat%send_src_idx(i),p_pat%send_src_blk(i)-lbound2+1)
       ENDDO
     ELSE
       DO i = 1, p_pat%n_send
         send_buf(1,i) = recv(p_pat%send_src_idx(i),p_pat%send_src_blk(i))
       ENDDO
     ENDIF

     CALL buffer_data(p_pat%send_limits, send_buf)
     RETURN

   ENDIF

   IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
     tmp_recv(:,1,:) = 0._wp
   ELSE
     tmp_recv(:,1,:) = recv(:,:)
   ENDIF

   IF (PRESENT(send) .AND. PRESENT(send_lbound2)) THEN
      IF (PRESENT(add)) THEN
         CALL exchange_data(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
                              add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)), &
                              send_lbound3=send_lbound2 )
      ELSE
         CALL exchange_data(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
                              send_lbound3=send_lbound2)
      ENDIF
   ELSEIF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
         CALL exchange_data(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
                              add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
      ENDIF
   ELSE
      IF (PRESENT(add)) THEN
         CALL exchange_data(p_pat, tmp_recv, &
                              add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data(p_pat, tmp_recv)
      ENDIF
   ENDIF

   recv(:,:) = tmp_recv(:,1,:)

END SUBROUTINE exchange_data_r2d
!================================================================================================
! INTEGER SECTION -------------------------------------------------------------------------------
! 
SUBROUTINE exchange_data_i2d(p_pat, recv, send, add, send_lbound2, l_recv_exists)
   !
   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   INTEGER, INTENT(INOUT), TARGET        :: recv(:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:)
   INTEGER, OPTIONAL :: send_lbound2
   LOGICAL, OPTIONAL :: l_recv_exists

   INTEGER :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))
   INTEGER :: send_buf(1,p_pat%n_send)
   INTEGER :: i, lbound2

   !-----------------------------------------------------------------------

   IF (use_exchange_delayed) THEN

     CALL add_delayed_request(p_pat, 1, .FALSE., p_int)

     delayed_request_i(n_delayed_requests_i)%recv2 => recv
     IF (PRESENT(add)) delayed_request_i(n_delayed_requests_i)%add2 => add

     IF (PRESENT(send) .AND. PRESENT(send_lbound2)) THEN
       lbound2 = send_lbound2
     ELSE
       lbound2 = 1
     ENDIF

     IF ( PRESENT(send)) THEN
       DO i = 1, p_pat%n_send
         send_buf(1,i) = send(p_pat%send_src_idx(i),p_pat%send_src_blk(i)-lbound2+1)
       ENDDO
     ELSE
       DO i = 1, p_pat%n_send
         send_buf(1,i) = recv(p_pat%send_src_idx(i),p_pat%send_src_blk(i))
       ENDDO
     ENDIF

     CALL buffer_data(p_pat%send_limits, send_buf)
     RETURN

   ENDIF

   IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
     tmp_recv(:,1,:) = 0
   ELSE
     tmp_recv(:,1,:) = recv(:,:)
   ENDIF

   IF (PRESENT(send) .AND. PRESENT(send_lbound2)) THEN
      IF (PRESENT(add)) THEN
         CALL exchange_data(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
                              add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)), &
                              send_lbound3=send_lbound2 )
      ELSE
         CALL exchange_data(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
                              send_lbound3=send_lbound2)
      ENDIF
   ELSEIF (PRESENT(send)) THEN
      IF (PRESENT(add)) THEN
         CALL exchange_data(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
                              add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data(p_pat, tmp_recv, &
                              send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
      ENDIF
   ELSE
      IF (PRESENT(add)) THEN
         CALL exchange_data(p_pat, tmp_recv, &
                              add =RESHAPE(add, (/SIZE(add ,1),1,SIZE(add ,2)/)))
      ELSE
         CALL exchange_data(p_pat, tmp_recv)
      ENDIF
   ENDIF

   recv(:,:) = tmp_recv(:,1,:)

END SUBROUTINE exchange_data_i2d
!================================================================================================
! LOGICAL SECTION -------------------------------------------------------------------------------
! 
SUBROUTINE exchange_data_l2d(p_pat, recv, send, send_lbound2, l_recv_exists)
   !
   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   LOGICAL, INTENT(INOUT), TARGET        :: recv(:,:)
   LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
   INTEGER, OPTIONAL :: send_lbound2
   LOGICAL, OPTIONAL :: l_recv_exists

   LOGICAL :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))
   LOGICAL :: send_buf(1,p_pat%n_send)
   INTEGER :: i, lbound2

   !-----------------------------------------------------------------------

   IF (use_exchange_delayed) THEN

     CALL add_delayed_request(p_pat, 1, .FALSE., p_bool)

     delayed_request_l(n_delayed_requests_l)%recv2 => recv

     IF (PRESENT(send) .AND. PRESENT(send_lbound2)) THEN
       lbound2 = send_lbound2
     ELSE
       lbound2 = 1
     ENDIF

     IF ( PRESENT(send)) THEN
       DO i = 1, p_pat%n_send
         send_buf(1,i) = send(p_pat%send_src_idx(i),p_pat%send_src_blk(i)-lbound2+1)
       ENDDO
     ELSE
       DO i = 1, p_pat%n_send
         send_buf(1,i) = recv(p_pat%send_src_idx(i),p_pat%send_src_blk(i))
       ENDDO
     ENDIF

     CALL buffer_data(p_pat%send_limits, send_buf)
     RETURN

   ENDIF

   IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
     tmp_recv(:,1,:) = .FALSE.
   ELSE
     tmp_recv(:,1,:) = recv(:,:)
   ENDIF

   IF (PRESENT(send) .AND. PRESENT(send_lbound2)) THEN
     CALL exchange_data(p_pat, tmp_recv, &
          &             send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)), &
          &             send_lbound3=send_lbound2)
   ELSEIF (PRESENT(send)) THEN
     CALL exchange_data(p_pat, tmp_recv, &
          &             send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
   ELSE
     CALL exchange_data(p_pat, tmp_recv)
   ENDIF

   recv(:,:) = tmp_recv(:,1,:)

END SUBROUTINE exchange_data_l2d
!
!================================================================================================

!-------------------------------------------------------------------------

!>
!! Interface for 2D arrays for exchange\_data\_reverse.
!!
!! Just reshapes the arrays and calls exchange\_data\_reverse\_3.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Dec 2009
!!
SUBROUTINE exchange_data_reverse_2(p_pat, recv, send)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   REAL(wp), INTENT(INOUT), TARGET        :: recv(:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)


! !LOCAL VARIABLES:

   REAL(wp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))
   REAL(wp) :: send_buf(1,p_pat%n_recv)
   INTEGER :: i

!EOP
!-----------------------------------------------------------------------
!BOC

   IF(use_exchange_delayed) THEN

     CALL add_delayed_request(p_pat, 1, .TRUE., p_real_dp)
     delayed_request_r(n_delayed_requests_r)%recv2 => recv

     IF(PRESENT(send)) THEN
       DO i = 1, p_pat%n_pnts
         send_buf(1,p_pat%recv_src(i)) = send(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i))
       ENDDO
     ELSE
       DO i = 1, p_pat%n_pnts
         send_buf(1,p_pat%recv_src(i)) = recv(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i))
       ENDDO
     ENDIF

     CALL buffer_data(p_pat%recv_limits, send_buf)
     RETURN

   ENDIF

   tmp_recv(:,1,:) = recv(:,:)

   IF(PRESENT(send)) THEN
      CALL exchange_data_reverse_3(p_pat, tmp_recv, &
                           send=RESHAPE(send,(/SIZE(send,1),1,SIZE(send,2)/)))
   ELSE
      CALL exchange_data_reverse_3(p_pat, tmp_recv)
   ENDIF

   recv(:,:) = tmp_recv(:,1,:)

END SUBROUTINE exchange_data_reverse_2

!--------------------------------------------------------------------------------------------------

!> Starts a delayed exchange, i.e. all following exchange_data_xxx calls are
!! buffered until do_delayed_exchange() is called which actually does the
!! communication calls and distributes the received data.
!!
!! ATTENTION:
!! Since a pointer to the recv and add arrays is stored internally,
!! these arrays must be a real reference (either the target array or a pointer)
!! and not a copy to the real data and  must neither go out of scope nor
!! be deallocated between the call to exchange_data_xxx and the corresponding
!! call to do_delayed_exchange().
!! Same holds for p_pat (but the communication patterns shouldn't be
!! changed after initialization anyways).
!!
!! Initial version by Rainer Johanni, Oct 2010

SUBROUTINE start_delayed_exchange

  INTEGER :: np

  ! Just ignore this call for 1 processor runs
  IF (my_process_is_mpi_seq()) RETURN

  ! If the first time here: initialize send_bufs_r

  IF(.NOT.ALLOCATED(send_bufs_r)) THEN

    ALLOCATE(send_bufs_r(0:p_n_work-1))

    DO np = 0, p_n_work-1
      ALLOCATE(send_bufs_r(np)%buf(1)) ! will be enlarged later
      send_bufs_r(np)%nelems = 0 ! nothing filled in yet
    ENDDO

  ENDIF

  use_exchange_delayed = .TRUE.

END SUBROUTINE start_delayed_exchange

!--------------------------------------------------------------------------------------------------
!
!> Adds a delayed request, i.e. increments n_delayed_requests_r, checks if the limit is not exceeded
!! sets the p_pat, ndim2, reverse members and nullifies the remaining pointers
!!
SUBROUTINE add_delayed_request(p_pat, ndim2, reverse, typ)

  TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
  INTEGER, INTENT(IN) :: ndim2
  LOGICAL, INTENT(IN) :: reverse
  INTEGER, INTENT(in) :: typ

  IF (typ == p_real_dp) THEN

    ! Check for limit
    IF (n_delayed_requests_r >= max_delayed_requests) &
         CALL finish('mo_communication','max_delayed_requests exceeded')

    n_delayed_requests_r = n_delayed_requests_r+1

    delayed_request_r(n_delayed_requests_r)%p_pat => p_pat
    delayed_request_r(n_delayed_requests_r)%ndim2 =  ndim2
    delayed_request_r(n_delayed_requests_r)%reverse = reverse
    delayed_request_r(n_delayed_requests_r)%recv2 => NULL()
    delayed_request_r(n_delayed_requests_r)%recv3 => NULL()
    delayed_request_r(n_delayed_requests_r)%add2  => NULL()
    delayed_request_r(n_delayed_requests_r)%add3  => NULL()

  ELSEIF (typ == p_int) THEN

    ! Check for limit
    IF (n_delayed_requests_i >= max_delayed_requests) &
         CALL finish('mo_communication','max_delayed_requests exceeded')
    
    n_delayed_requests_i = n_delayed_requests_i+1

    delayed_request_i(n_delayed_requests_i)%p_pat => p_pat
    delayed_request_i(n_delayed_requests_i)%ndim2 =  ndim2
    delayed_request_i(n_delayed_requests_i)%reverse = reverse
    delayed_request_i(n_delayed_requests_i)%recv2 => NULL()
    delayed_request_i(n_delayed_requests_i)%recv3 => NULL()
    delayed_request_i(n_delayed_requests_i)%add2  => NULL()
    delayed_request_i(n_delayed_requests_i)%add3  => NULL()

  ELSEIF (typ == p_bool) THEN

    ! Check for limit

    IF (n_delayed_requests_l >= max_delayed_requests) &
         CALL finish('mo_communication','max_delayed_requests exceeded')

    n_delayed_requests_l = n_delayed_requests_l+1

    delayed_request_l(n_delayed_requests_l)%p_pat => p_pat
    delayed_request_l(n_delayed_requests_l)%ndim2 =  ndim2
    delayed_request_l(n_delayed_requests_l)%reverse = reverse
    delayed_request_l(n_delayed_requests_l)%recv2 => NULL()
    delayed_request_l(n_delayed_requests_l)%recv3 => NULL()
    
  ENDIF

END SUBROUTINE add_delayed_request
!
!================================================================================================
!
!--------------------------------------------------------------------------------------------------
!
!> Backend for buffering data for delayed exchanges
!!
!! Initial version by Rainer Johanni, Oct 2010
!!
!================================================================================================
! REAL SECTION ----------------------------------------------------------------------------------
! 
SUBROUTINE buffer_data_r(limits, send)

  INTEGER, INTENT(IN) :: limits(0:)
  REAL(wp), INTENT(IN) :: send(:,:)

  REAL(wp), ALLOCATABLE :: buf(:)
  INTEGER :: ndim2, nitems, n, np, iss, ise

  ndim2 = UBOUND(send,1)

  ! Store data in send_bufs_r

  DO np = 0, p_n_work-1 ! loop over PEs where to send the data

    iss = limits(np)+1
    ise = limits(np+1)

    nitems = ndim2*(ise-iss+1)

    IF(ise<iss) CYCLE

    n = send_bufs_r(np)%nelems ! Number of elements currently in buffer

    ! Make sure that send_bufs_r(np)%buf has enough capacity.
    ! There is some deallocating/allocating of the buffers at the beginning
    ! but after the first timestep the final size should be reached.

    IF(n + nitems > SIZE(send_bufs_r(np)%buf)) THEN
      IF(n>0) THEN
        ALLOCATE(buf(n))
        buf(1:n) = send_bufs_r(np)%buf(1:n)
      ENDIF
      DEALLOCATE(send_bufs_r(np)%buf)
      ALLOCATE(send_bufs_r(np)%buf(n + nitems))
      IF(n>0) THEN
        send_bufs_r(np)%buf(1:n) = buf(1:n)
        DEALLOCATE(buf)
      ENDIF
    ENDIF

    send_bufs_r(np)%buf(n+1:n+nitems) = RESHAPE(send(:,iss:ise), (/ nitems /))

    send_bufs_r(np)%nelems = n+nitems ! new number of elements

  ENDDO

END SUBROUTINE buffer_data_r
!
!================================================================================================
! INTEGER SECTION -------------------------------------------------------------------------------
!
SUBROUTINE buffer_data_i(limits, send)

  INTEGER, INTENT(IN) :: limits(0:)
  INTEGER, INTENT(IN) :: send(:,:)

  INTEGER, ALLOCATABLE :: buf(:)
  INTEGER :: ndim2, nitems, n, np, iss, ise

  ndim2 = UBOUND(send,1)

  ! Store data in send_bufs_i

  DO np = 0, p_n_work-1 ! loop over PEs where to send the data

    iss = limits(np)+1
    ise = limits(np+1)

    nitems = ndim2*(ise-iss+1)

    IF(ise<iss) CYCLE

    n = send_bufs_i(np)%nelems ! Number of elements currently in buffer

    ! Make sure that send_bufs_i(np)%buf has enough capacity.
    ! There is some deallocating/allocating of the buffers at the beginning
    ! but after the first timestep the final size should be reached.

    IF(n + nitems > SIZE(send_bufs_i(np)%buf)) THEN
      IF(n>0) THEN
        ALLOCATE(buf(n))
        buf(1:n) = send_bufs_i(np)%buf(1:n)
      ENDIF
      DEALLOCATE(send_bufs_i(np)%buf)
      ALLOCATE(send_bufs_i(np)%buf(n + nitems))
      IF(n>0) THEN
        send_bufs_i(np)%buf(1:n) = buf(1:n)
        DEALLOCATE(buf)
      ENDIF
    ENDIF

    send_bufs_i(np)%buf(n+1:n+nitems) = RESHAPE(send(:,iss:ise), (/ nitems /))

    send_bufs_i(np)%nelems = n+nitems ! new number of elements

  ENDDO

END SUBROUTINE buffer_data_i
!
!================================================================================================
! LOGICAL SECTION -------------------------------------------------------------------------------
!
SUBROUTINE buffer_data_l(limits, send)

  INTEGER, INTENT(IN) :: limits(0:)
  LOGICAL, INTENT(IN) :: send(:,:)

  LOGICAL, ALLOCATABLE :: buf(:)
  INTEGER :: ndim2, nitems, n, np, iss, ise

  ndim2 = UBOUND(send,1)

  ! Store data in send_bufs_l

  DO np = 0, p_n_work-1 ! loop over PEs where to send the data

    iss = limits(np)+1
    ise = limits(np+1)

    nitems = ndim2*(ise-iss+1)

    IF(ise<iss) CYCLE

    n = send_bufs_l(np)%nelems ! Number of elements currently in buffer

    ! Make sure that send_bufs_l(np)%buf has enough capacity.
    ! There is some deallocating/allocating of the buffers at the beginning
    ! but after the first timestep the final size should be reached.

    IF(n + nitems > SIZE(send_bufs_l(np)%buf)) THEN
      IF(n>0) THEN
        ALLOCATE(buf(n))
        buf(1:n) = send_bufs_l(np)%buf(1:n)
      ENDIF
      DEALLOCATE(send_bufs_l(np)%buf)
      ALLOCATE(send_bufs_l(np)%buf(n + nitems))
      IF(n>0) THEN
        send_bufs_l(np)%buf(1:n) = buf(1:n)
        DEALLOCATE(buf)
      ENDIF
    ENDIF

    send_bufs_l(np)%buf(n+1:n+nitems) = RESHAPE(send(:,iss:ise), (/ nitems /))

    send_bufs_l(np)%nelems = n+nitems ! new number of elements

  ENDDO

END SUBROUTINE buffer_data_l
!
!================================================================================================
!
!--------------------------------------------------------------------------------------------------

SUBROUTINE do_delayed_exchange()

!
! !DESCRIPTION:
! Does data exchange according to a communication pattern (in p_pat).
!
! !REVISION HISTORY:
! Initial version by Rainer Johanni, Oct 2010
!

! !INPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

! !LOCAL VARIABLES:

  TYPE(t_buffer_r), ALLOCATABLE :: recv_bufs_r(:)
  TYPE(t_buffer_i), ALLOCATABLE :: recv_bufs_i(:)
  TYPE(t_buffer_l), ALLOCATABLE :: recv_bufs_l(:)
!
  REAL(wp), ALLOCATABLE :: recv_buf(:,:)
  TYPE(t_comm_pattern), POINTER :: p_pat

  INTEGER :: i, n, np, nr, irs, ire, ndim2
#ifdef __SX__
  INTEGER :: k
#endif

!EOP
!-----------------------------------------------------------------------

  ! Just ignore this call for 1 processor runs
  IF (my_process_is_mpi_seq()) RETURN

  ALLOCATE(recv_bufs_r(0:p_n_work-1))
  ALLOCATE(recv_bufs_i(0:p_n_work-1))
  ALLOCATE(recv_bufs_l(0:p_n_work-1))

  ! Set up irecv's for receive buffers

  DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

    n = 0
    DO nr = 1, n_delayed_requests_r
      IF(delayed_request_r(nr)%reverse) THEN
        irs = delayed_request_r(nr)%p_pat%send_limits(np)+1
        ire = delayed_request_r(nr)%p_pat%send_limits(np+1)
      ELSE
        irs = delayed_request_r(nr)%p_pat%recv_limits(np)+1
        ire = delayed_request_r(nr)%p_pat%recv_limits(np+1)
      ENDIF
      n = n + (ire-irs+1)*delayed_request_r(nr)%ndim2
    ENDDO

    IF(n>0) THEN
      ALLOCATE(recv_bufs_r(np)%buf(n))
      CALL p_irecv(recv_bufs_r(np)%buf, np, 1, p_count=n, comm=p_comm_work)
    ENDIF

  ENDDO

  ! Send our data

  DO np = 0, p_n_work-1 ! loop over PEs where to send the data

     IF(send_bufs_r(np)%nelems > 0) &
        CALL p_send(send_bufs_r(np)%buf, np, 1, p_count=send_bufs_r(np)%nelems, comm=p_comm_work)

  ENDDO

  ! Wait for all outstanding requests to finish

  CALL p_wait

  recv_bufs_r(:)%nelems = 0 ! counts consumed elements below
  recv_bufs_i(:)%nelems = 0 ! counts consumed elements below
  recv_bufs_l(:)%nelems = 0 ! counts consumed elements below

  ! Fill in receive buffer

  DO nr = 1, n_delayed_requests_r

    p_pat => delayed_request_r(nr)%p_pat
    ndim2 =  delayed_request_r(nr)%ndim2

    IF(delayed_request_r(nr)%reverse) THEN
      ALLOCATE(recv_buf(ndim2,p_pat%n_send))
    ELSE
      ALLOCATE(recv_buf(ndim2,p_pat%n_recv))
    ENDIF

    DO np = 0, p_n_work-1 ! loop over PEs from where to receive the data

      IF(delayed_request_r(nr)%reverse) THEN
        irs = p_pat%send_limits(np)+1
        ire = p_pat%send_limits(np+1)
      ELSE
        irs = p_pat%recv_limits(np)+1
        ire = p_pat%recv_limits(np+1)
      ENDIF
      IF(ire<irs) CYCLE

      n = recv_bufs_r(np)%nelems
      recv_buf(1:ndim2,irs:ire) = RESHAPE(recv_bufs_r(np)%buf(n+1:n+ndim2*(ire-irs+1)), &
        &                                 (/ndim2,ire-irs+1/))
      recv_bufs_r(np)%nelems = recv_bufs_r(np)%nelems + ndim2*(ire-irs+1)
    ENDDO

    IF(delayed_request_r(nr)%reverse) THEN
      IF(ASSOCIATED(delayed_request_r(nr)%recv3)) THEN
#ifdef __SX__
!CDIR UNROLL=6
        DO k = 1, ndim2
          DO i = 1, p_pat%n_send
            delayed_request_r(nr)%recv3(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)) = &
              recv_buf(k,i)
          ENDDO
        ENDDO
#else
        DO i = 1, p_pat%n_send
          delayed_request_r(nr)%recv3(p_pat%send_src_idx(i),1:ndim2,p_pat%send_src_blk(i)) = &
            recv_buf(1:ndim2,i)
        ENDDO
#endif
      ELSE
        DO i = 1, p_pat%n_send
          delayed_request_r(nr)%recv2(p_pat%send_src_idx(i),p_pat%send_src_blk(i)) = &
            recv_buf(1,i)
        ENDDO
      ENDIF
    ELSE
      IF(ASSOCIATED(delayed_request_r(nr)%recv3)) THEN
        IF(ASSOCIATED(delayed_request_r(nr)%add3)) THEN
#ifdef __SX__
!CDIR UNROLL=6
          DO k = 1, ndim2
            DO i = 1, p_pat%n_pnts
              delayed_request_r(nr)%recv3(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
                recv_buf(k,p_pat%recv_src(i)) +                                          &
                delayed_request_r(nr)%add3(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
            ENDDO
          ENDDO
#else
          DO i = 1, p_pat%n_pnts
            delayed_request_r(nr)%recv3(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i)) = &
              recv_buf(1:ndim2,p_pat%recv_src(i)) +                                          &
              delayed_request_r(nr)%add3(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i))
          ENDDO
#endif
        ELSE
#ifdef __SX__
!CDIR UNROLL=6
          DO k = 1, ndim2
            DO i = 1, p_pat%n_pnts
              delayed_request_r(nr)%recv3(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
                recv_buf(k,p_pat%recv_src(i))
            ENDDO
          ENDDO
#else
          DO i = 1, p_pat%n_pnts
            delayed_request_r(nr)%recv3(p_pat%recv_dst_idx(i),1:ndim2,p_pat%recv_dst_blk(i)) = &
              recv_buf(1:ndim2,p_pat%recv_src(i))
          ENDDO
#endif
        ENDIF
      ELSE
        IF(ASSOCIATED(delayed_request_r(nr)%add2)) THEN
          DO i = 1, p_pat%n_pnts
            delayed_request_r(nr)%recv2(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i)) = &
              recv_buf(1,p_pat%recv_src(i)) +                                        &
              delayed_request_r(nr)%add2(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i))
          ENDDO
        ELSE
          DO i = 1, p_pat%n_pnts
            delayed_request_r(nr)%recv2(p_pat%recv_dst_idx(i),p_pat%recv_dst_blk(i)) = &
              recv_buf(1,p_pat%recv_src(i))
          ENDDO
        ENDIF
      ENDIF
    ENDIF

    DEALLOCATE(recv_buf)

  ENDDO

  DO np = 0, p_n_work-1
     IF(ALLOCATED(recv_bufs_r(np)%buf)) DEALLOCATE(recv_bufs_r(np)%buf)
     IF(ALLOCATED(recv_bufs_i(np)%buf)) DEALLOCATE(recv_bufs_i(np)%buf)
     IF(ALLOCATED(recv_bufs_l(np)%buf)) DEALLOCATE(recv_bufs_l(np)%buf)
  ENDDO

  DEALLOCATE(recv_bufs_r)
  DEALLOCATE(recv_bufs_i)
  DEALLOCATE(recv_bufs_l)

  n_delayed_requests_r = 0

  DO np = 0, p_n_work-1
     send_bufs_r(np)%nelems = 0
     send_bufs_i(np)%nelems = 0
     send_bufs_l(np)%nelems = 0
  ENDDO

  use_exchange_delayed = .FALSE.

END SUBROUTINE do_delayed_exchange

END MODULE mo_communication
