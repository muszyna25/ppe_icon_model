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
#if (defined(_CRAYFTN) )
#define __OMPPAR_COPY__
#endif
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
USE mo_exception,          ONLY: finish, message, message_text
USE mo_mpi,                ONLY: p_send, p_recv, p_irecv, p_wait, p_isend, &
     & p_comm_work, my_process_is_mpi_seq, p_pe_work, p_n_work, &
     & get_my_mpi_work_communicator, get_my_mpi_work_comm_size, &
     & get_my_mpi_work_id, p_gather, p_gatherv, work_mpi_barrier, &
     & p_alltoallv, p_alltoall, process_mpi_root_id, p_bcast
USE mo_parallel_config, ONLY: iorder_sendrecv, nproma, itype_exch_barrier
USE mo_timer,           ONLY: timer_start, timer_stop, activate_sync_timers, &
  & timer_exch_data, timer_exch_data_async, timer_barrier, timer_exch_data_wait
USE mo_run_config,      ONLY: msg_level
USE mo_decomposition_tools, ONLY: t_glb2loc_index_lookup, get_local_index
USE mo_util_sort,          ONLY: quicksort


IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

!modules interface-------------------------------------------
!subroutines
PUBLIC :: blk_no, idx_no, idx_1d
PUBLIC :: setup_comm_pattern, delete_comm_pattern, exchange_data,  &
          exchange_data_mult, exchange_data_grf,                   &
          start_async_comm, complete_async_comm,                   &
          exchange_data_4de3, delete_comm_gather_pattern
PUBLIC :: t_comm_pattern
PUBLIC :: reorder_comm_pattern
PUBLIC :: reorder_comm_pattern_snd
PUBLIC :: reorder_comm_pattern_rcv
PUBLIC :: reorder_comm_gather_pattern

PUBLIC :: t_comm_gather_pattern
PUBLIC :: setup_comm_gather_pattern

PUBLIC :: ASSIGNMENT(=)
!
!variables

!--------------------------------------------------------------------------------------------------
!
TYPE t_comm_pattern

   ! Number of points we receive in communication,
   ! this is the same as recv_limits

   INTEGER :: n_recv  ! Number of points we receive from other PEs
   INTEGER :: n_pnts  ! Number of points we output into local array;
                      ! this may be bigger than n_recv due to
                      ! duplicate entries
   INTEGER :: n_send  ! Number of points we send to other PEs

   INTEGER :: np_recv ! Number of PEs from which data have to be received
   INTEGER :: np_send ! Number of PEs to which data have to be sent

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
   ! The local PE recvs recv_count(i) data items to PE pelist_recv(i),
   ! starting at recv_startidx(i) in the receiver buffer.
   INTEGER, ALLOCATABLE :: recv_startidx(:)
   INTEGER, ALLOCATABLE :: recv_count(:)

END TYPE t_comm_pattern

!
!------------------------------------------------------------------------------------------------
!

TYPE t_comm_gather_pattern
  INTEGER, ALLOCATABLE :: collector_pes(:) ! ranks of collector processes
  INTEGER, ALLOCATABLE :: collector_size(:) ! total number of points per
                                            ! collector
  INTEGER, ALLOCATABLE :: collector_send_size(:) ! local number of points per
                                                 ! collector
  INTEGER, ALLOCATABLE :: loc_index(:) ! local indices of all points that
                                       ! need to be sent to a collector
  INTEGER, ALLOCATABLE :: recv_buffer_reorder(:) ! once the data is received
                                                 ! on the collectors, it has
                                                 ! to be reordered according
                                                 ! to this array

  INTEGER, ALLOCATABLE :: recv_pes(:) ! ranks from which data is to be received
  INTEGER, ALLOCATABLE :: recv_size(:) ! number of remote points received
END TYPE t_comm_gather_pattern

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE copy_t_comm_gather_pattern
END INTERFACE

!--------------------------------------------------------------------------------------------------
!

INTERFACE exchange_data
   MODULE PROCEDURE exchange_data_r3d
   MODULE PROCEDURE exchange_data_i3d
   MODULE PROCEDURE exchange_data_l3d
   MODULE PROCEDURE exchange_data_r2d
   MODULE PROCEDURE exchange_data_i2d
   MODULE PROCEDURE exchange_data_l2d
   MODULE PROCEDURE gather_r_2d_deblock
   MODULE PROCEDURE gather_r_1d_deblock
   MODULE PROCEDURE gather_i_2d_deblock
   MODULE PROCEDURE gather_i_1d_deblock
END INTERFACE

INTERFACE exchange_data_seq
   MODULE PROCEDURE exchange_data_r3d_seq
   MODULE PROCEDURE exchange_data_r2d_seq
END INTERFACE


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
!! Note: This setup routine works only for the trivial communication
!!       patterns in sequential runs.
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
!!                There may be more than 1 point with the same global index,
!!                in this case the point is exchanged only once and
!!                locally distributed.
!!                - If this argument is not present, we assume global_index=1,2.3,...
!!
!! send_decomp_info domain decomposition information for the SENDER array
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
SUBROUTINE setup_comm_pattern(n_points, owner, opt_global_index, &
                              send_glb2loc_index, p_pat)

!

   INTEGER, INTENT(IN)           :: n_points             ! Total number of points
   INTEGER, INTENT(IN)           :: owner(:)             ! Owner of every point
   INTEGER, INTENT(IN), OPTIONAL :: opt_global_index(:)  ! Global index of every point
   TYPE(t_glb2loc_index_lookup), INTENT(IN) :: send_glb2loc_index
                                                         ! global to local index
                                                         ! lookup information
                                                         ! of the SENDER array

   TYPE(t_comm_pattern), INTENT(INOUT) :: p_pat


   INTEGER, ALLOCATABLE :: icnt(:), flag(:), global_recv_index(:), send_src(:), num_rcv(:)
   INTEGER              :: global_index(n_points)
   INTEGER :: i, n, np, nr, num_recv, irs, ire, num_send, iss, ise, max_glb

!-----------------------------------------------------------------------

   IF (PRESENT(opt_global_index)) THEN
     global_index(:) = opt_global_index(1:n_points)
   ELSE
     global_index(:) = (/ (i, i=1,n_points) /)
   END IF

   ALLOCATE(icnt(0:p_n_work-1), num_rcv(0:p_n_work-1))
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
     num_rcv(np) = p_pat%recv_limits(np+1) - p_pat%recv_limits(np)
     ! First send the number of points to be received
     IF (np /= p_pe_work) CALL p_isend(num_rcv(np), np, 1, comm=p_comm_work)
   ENDDO

   ! Now, we receive the number of points are needed from us
   DO nr = 0, p_n_work-1
     
     IF(nr /= p_pe_work) THEN
       CALL p_recv(icnt(nr), nr, 1,  comm=p_comm_work)
     ELSE
       icnt(nr) = num_rcv(nr)
     ENDIF
     
   ENDDO
   CALL p_wait

   ! Now send the global index of the points we need from PE np
   DO np = 0, p_n_work-1 ! loop over PEs where to send the data

     IF (np == p_pe_work) CYCLE

     irs = p_pat%recv_limits(np)+1 ! Start index in global_recv_index
     ire = p_pat%recv_limits(np+1) ! End   index in global_recv_index

     IF(num_rcv(np)>0) CALL p_isend(global_recv_index(irs), np, 1, &
        p_count=ire-irs+1, comm=p_comm_work)

   ENDDO

   irs = p_pat%recv_limits(p_pe_work)+1   ! Start index in global_recv_index
   ire = p_pat%recv_limits(p_pe_work+1)   ! End   index in global_recv_index

   DEALLOCATE(num_rcv)
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
       IF(num_send>0) CALL p_recv(send_src(iss), nr, 1, &
         p_count=ise-iss+1, comm=p_comm_work)
     ELSE
       IF(num_send>0) send_src(iss:ise) = global_recv_index(irs:ire)
     ENDIF
   ENDDO

   CALL p_wait

   ALLOCATE(p_pat%send_src_blk(p_pat%n_send))
   ALLOCATE(p_pat%send_src_idx(p_pat%n_send))

   ! The indices in p_pat%send_src are global, convert to local

   DO i = 1, p_pat%n_send

      np = get_local_index(send_glb2loc_index, send_src(i))
      IF(np <= 0) CALL finish('setup_comm_pattern','Got illegal index')
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

   ! consistency check of communication pattern
#ifndef NOMPI
   IF (msg_level >= 25)  &
     CALL check_comm_pattern(p_pat)
#endif

END SUBROUTINE setup_comm_pattern

!-------------------------------------------------------------------------

SUBROUTINE setup_comm_gather_pattern(global_size, owner_local, glb_index, &
  &                                  gather_pattern)
  INTEGER, INTENT(IN) :: global_size, owner_local(:), glb_index(:)
  TYPE(t_comm_gather_pattern), INTENT(INOUT) :: gather_pattern

  INTEGER :: num_collectors
  LOGICAL, ALLOCATABLE :: pack_mask(:)
  INTEGER :: num_local_points, num_points_per_coll
  INTEGER, ALLOCATABLE :: packed_glb_index(:)
  INTEGER :: coll_stride
  INTEGER :: num_send_per_process(p_n_work), num_recv_per_process(p_n_work)
  INTEGER :: send_displ(p_n_work+1), recv_displ(p_n_work)
  INTEGER, ALLOCATABLE :: send_buffer(:), recv_buffer(:)
  INTEGER :: num_recv, num_recv_points
  INTEGER :: i, j, n

  ! determine collector ranks and the data associated to each collector
  num_collectors = NINT(SQRT(REAL(p_n_work)))
  IF (ALLOCATED(gather_pattern%collector_pes)) &
    DEALLOCATE(gather_pattern%collector_pes)
  IF (ALLOCATED(gather_pattern%collector_size)) &
    DEALLOCATE(gather_pattern%collector_size)
  IF (ALLOCATED(gather_pattern%collector_send_size)) &
    DEALLOCATE(gather_pattern%collector_send_size)
  ALLOCATE(gather_pattern%collector_pes(num_collectors), &
    &      gather_pattern%collector_size(num_collectors), &
    &      gather_pattern%collector_send_size(num_collectors))
  coll_stride = (p_n_work + num_collectors - 1) / &
    &           num_collectors
  num_points_per_coll = (global_size + num_collectors - 1) / num_collectors
  DO i = 1, num_collectors
    ! set collector ranks
    gather_pattern%collector_pes(i) = (i-1) * coll_stride
  END DO

  ! mask for all locally owned points
  ALLOCATE(pack_mask(SIZE(owner_local)))
  pack_mask(:) = owner_local(:) == p_pe_work
  num_local_points = COUNT(pack_mask(:))

  ! determine local indices of all points that need to be sent to a collector
  IF (ALLOCATED(gather_pattern%loc_index)) &
    DEALLOCATE(gather_pattern%loc_index)
  ALLOCATE(gather_pattern%loc_index(num_local_points), &
    &      packed_glb_index(num_local_points))
  packed_glb_index(:) = PACK(glb_index(:), pack_mask(:))
  gather_pattern%loc_index(:) = PACK((/(i, i = 1, &
    &                                SIZE(owner_local))/), pack_mask(:))

  DEALLOCATE(pack_mask)

  ! sort loc_index according to the respective global indices
  CALL quicksort(packed_glb_index(:), gather_pattern%loc_index(:))

  ! determine number of points that need to be sent to each collector
  gather_pattern%collector_send_size(:) = 0
  DO i = 1, num_local_points
    n = 1 + (packed_glb_index(i) - 1) / num_points_per_coll
    gather_pattern%collector_send_size(n) = &
      gather_pattern%collector_send_size(n) + 1
  END DO

  ! generate send and receive counts for all processes
  num_send_per_process(:) = 0
  DO i = 1, num_collectors
    num_send_per_process(gather_pattern%collector_pes(i)+1) = &
      gather_pattern%collector_send_size(i)
  END DO
  CALL p_alltoall(num_send_per_process(:), num_recv_per_process(:), &
    &             p_comm_work)
  num_recv_points = SUM(num_recv_per_process(:))

  ! exchange number of points per collector
  IF (p_pe_work == gather_pattern%collector_pes(1)) THEN
    gather_pattern%collector_size(1) = num_recv_points
    DO i = 2, num_collectors
      CALL p_recv(gather_pattern%collector_size(i), &
        &         gather_pattern%collector_pes(i), 0, 1, p_comm_work)
    END DO
  ELSE IF (ANY(gather_pattern%collector_pes(:) == p_pe_work)) THEN
    CALL p_send(num_recv_points, gather_pattern%collector_pes(1), 0, 1, &
      &         p_comm_work)
  END IF
  CALL p_bcast(gather_pattern%collector_size(:), &
    &          gather_pattern%collector_pes(1), p_comm_work)

  ! number of messages to be received (is 0 on non-collector processes)
  num_recv = COUNT(num_recv_per_process(:) /= 0)
  IF (ALLOCATED(gather_pattern%recv_pes)) &
    DEALLOCATE(gather_pattern%recv_pes)
  IF (ALLOCATED(gather_pattern%recv_size)) &
    DEALLOCATE(gather_pattern%recv_size)
  ALLOCATE(gather_pattern%recv_pes(num_recv), &
    &      gather_pattern%recv_size(num_recv))
  num_recv = 0
  DO i = 1, p_n_work
    IF (num_recv_per_process(i) /= 0) THEN
      num_recv = num_recv + 1
      gather_pattern%recv_pes(num_recv) = i - 1
      gather_pattern%recv_size(num_recv) = num_recv_per_process(i)
    END IF
  END DO

  ! generate send and receive displacement and fill send buffer
  ! remark: at first the content of send_displ is shifted by one element
  !         to the back, this is done in order to ease the copying of data
  !         into the send buffer
  send_displ(1:2) = 0
  recv_displ(1) = 0
  DO i = 2, p_n_work
    send_displ(i+1) = send_displ(i) + num_send_per_process(i-1)
    recv_displ(i)   = recv_displ(i-1) + num_recv_per_process(i-1)
  END DO
  ALLOCATE(send_buffer(SUM(num_send_per_process(:))), &
    &      recv_buffer(num_recv_points))

  DO i = 1, num_local_points
    j = 2 + gather_pattern%collector_pes(1 + (packed_glb_index(i)-1) / &
      &                                  num_points_per_coll)
    send_displ(j) = send_displ(j) + 1
    send_buffer(send_displ(j)) = packed_glb_index(i)
  END DO

  DEALLOCATE(packed_glb_index)

  ! collect the global indices from all processes
  CALL p_alltoallv(send_buffer, num_send_per_process, send_displ, &
    &              recv_buffer, num_recv_per_process, recv_displ, &
    &              p_comm_work)

  ! compute the final position of received data on the collectors
  IF (ALLOCATED(gather_pattern%recv_buffer_reorder)) &
    DEALLOCATE(gather_pattern%recv_buffer_reorder)
  ALLOCATE(gather_pattern%recv_buffer_reorder(num_recv_points))
  gather_pattern%recv_buffer_reorder(:) = (/(i, i = 1, num_recv_points)/)
  CALL quicksort(recv_buffer(:), gather_pattern%recv_buffer_reorder(:))

  DEALLOCATE(send_buffer, recv_buffer)

END SUBROUTINE setup_comm_gather_pattern

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

SUBROUTINE delete_comm_gather_pattern(gather_pattern)
  TYPE(t_comm_gather_pattern), INTENT(INOUT) :: gather_pattern

  IF (ALLOCATED(gather_pattern%collector_pes)) &
    DEALLOCATE(gather_pattern%collector_pes)
  IF (ALLOCATED(gather_pattern%collector_size)) &
    DEALLOCATE(gather_pattern%collector_size)
  IF (ALLOCATED(gather_pattern%collector_send_size)) &
    DEALLOCATE(gather_pattern%collector_send_size)
  IF (ALLOCATED(gather_pattern%loc_index)) &
    DEALLOCATE(gather_pattern%loc_index)
  IF (ALLOCATED(gather_pattern%recv_buffer_reorder)) &
    DEALLOCATE(gather_pattern%recv_buffer_reorder)
  IF (ALLOCATED(gather_pattern%recv_pes)) &
    DEALLOCATE(gather_pattern%recv_pes)
  IF (ALLOCATED(gather_pattern%recv_size)) &
    DEALLOCATE(gather_pattern%recv_size)
END SUBROUTINE delete_comm_gather_pattern

ELEMENTAL SUBROUTINE copy_t_comm_gather_pattern(out_arg, in_arg)

  TYPE(t_comm_gather_pattern), INTENT(OUT) :: out_arg
  TYPE(t_comm_gather_pattern), INTENT(IN) :: in_arg

  IF (ALLOCATED(in_arg%collector_pes)) THEN
    ALLOCATE(out_arg%collector_pes(SIZE(in_arg%collector_pes)))
    out_arg%collector_pes(:) = in_arg%collector_pes(:)
  END IF
  IF (ALLOCATED(in_arg%collector_size)) THEN
    ALLOCATE(out_arg%collector_size(SIZE(in_arg%collector_size)))
    out_arg%collector_size(:) = in_arg%collector_size(:)
  END IF
  IF (ALLOCATED(in_arg%collector_send_size)) THEN
    ALLOCATE(out_arg%collector_send_size(SIZE(in_arg%collector_send_size)))
    out_arg%collector_send_size(:) = in_arg%collector_send_size(:)
  END IF
  IF (ALLOCATED(in_arg%loc_index)) THEN
    ALLOCATE(out_arg%loc_index(SIZE(in_arg%loc_index)))
    out_arg%loc_index(:) = in_arg%loc_index(:)
  END IF
  IF (ALLOCATED(in_arg%recv_buffer_reorder)) THEN
    ALLOCATE(out_arg%recv_buffer_reorder(SIZE(in_arg%recv_buffer_reorder)))
    out_arg%recv_buffer_reorder(:) = in_arg%recv_buffer_reorder(:)
  END IF
  IF (ALLOCATED(in_arg%recv_pes)) THEN
    ALLOCATE(out_arg%recv_pes(SIZE(in_arg%recv_pes)))
    out_arg%recv_pes(:) = in_arg%recv_pes(:)
  END IF
  IF (ALLOCATED(in_arg%recv_size)) THEN
    ALLOCATE(out_arg%recv_size(SIZE(in_arg%recv_size)))
    out_arg%recv_size(:) = in_arg%recv_size(:)
  END IF

END SUBROUTINE copy_t_comm_gather_pattern

!-------------------------------------------------------------------------
!> Consistency check of communication pattern.
!! Sends pattern info to working PE 0, which checks this data
!! for consistency wrt. send/receive counts.
!!
!! @par Revision History
!! Initial implementation  by  F. Prill, DWD (2012-01-20)
!!
SUBROUTINE check_comm_pattern(p_pat)
  TYPE(t_comm_pattern), INTENT(INOUT) :: p_pat

  ! local variables
  CHARACTER(*), PARAMETER :: routine = TRIM("mo_communication:check_comm_pattern")
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

   ! special treatment for trivial communication patterns of
   ! sequential runs
   IF(my_process_is_mpi_seq()) THEN
     CALL exchange_data_r3d_seq(p_pat, recv, send, add, send_lbound3)
     RETURN
   END IF

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL work_mpi_barrier()
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   IF (activate_sync_timers) CALL timer_start(timer_exch_data)

   IF(SIZE(recv,1) /= nproma) THEN
     CALL finish('exchange_data','Illegal first dimension of data array')
   ENDIF

   ndim2 = SIZE(recv,2)

   IF (iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) THEN
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
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO
#endif
     DO i = 1, p_pat%n_send
       send_buf(1:ndim2,i) = send_ptr(p_pat%send_src_idx(i),1:ndim2,   &
         &                            p_pat%send_src_blk(i)-lbound3+1)
     ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
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

   IF (iorder_sendrecv > 0) THEN
     ! Wait for all outstanding requests to finish
     IF (activate_sync_timers) CALL timer_start(timer_exch_data_wait)
     CALL p_wait
     IF (activate_sync_timers) CALL timer_stop(timer_exch_data_wait)

   ELSE IF (iorder_sendrecv == 0) THEN
     ! dummy mode; just fill the receive buffer with "something reasonable"
     IF (p_pat%n_send > 0) THEN
       DO k = 1, ndim2
         recv_buf(k,:) = send_buf(k,1)
       ENDDO
     ELSE
       DO k = 1, ndim2
         recv_buf(k,:) = 1._wp
       ENDDO
     ENDIF
   ENDIF

   IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL work_mpi_barrier()
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

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

   IF (activate_sync_timers) CALL timer_stop(timer_exch_data)
   
END SUBROUTINE exchange_data_r3d



! SEQUENTIAL version of subroutine "exchange_data_r3d"
!
SUBROUTINE exchange_data_r3d_seq(p_pat, recv, send, add, send_lbound3)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   REAL(wp), INTENT(INOUT), TARGET        :: recv(:,:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
   INTEGER, OPTIONAL :: send_lbound3
   ! local variables
    CHARACTER(*), PARAMETER :: routine = "mo_communication:exchange_data_r3d_seq"
   INTEGER :: i, lbound3, ndim2

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
   IF(PRESENT(send) .AND. PRESENT(send_lbound3)) THEN
     lbound3 = send_lbound3 - 1
   ELSE
     lbound3 = 0
   ENDIF

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
     DO i=1,p_pat%n_pnts
       recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
         &  add( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )                +  &
         &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
         &       1:ndim2,                                                                  &
         &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound3)
     END DO
   ELSE
     DO i=1,p_pat%n_pnts
       recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
         &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
         &       1:ndim2,                                                                  &
         &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound3)
     END DO
   END IF

 END SUBROUTINE exchange_data_r3d_seq



!================================================================================================
! INTEGER SECTION -------------------------------------------------------------------------------
! 
SUBROUTINE exchange_data_i3d(p_pat, recv, send, add, send_lbound3)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   INTEGER, INTENT(INOUT), TARGET        :: recv(:,:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
   INTEGER, OPTIONAL :: send_lbound3

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_i3d"
   INTEGER :: send_buf(SIZE(recv,2),p_pat%n_send), &
               recv_buf(SIZE(recv,2),p_pat%n_recv)

   INTEGER, POINTER :: send_ptr(:,:,:)

   INTEGER :: i, k, np, irs, iss, pid, icount, ndim2, lbound3

   IF(my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Not yet implemented!")
   END IF

   !-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_exch_data)

   IF(my_process_is_mpi_seq()) &
     CALL finish('exchange_data','must not be called on single PE/test PE')

   IF(SIZE(recv,1) /= nproma) THEN
     CALL finish('exchange_data','Illegal first dimension of data array')
   ENDIF

   ndim2 = SIZE(recv,2)

   IF (iorder_sendrecv == 1 .OR. iorder_sendrecv >= 3) THEN
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
   ELSE IF (iorder_sendrecv >= 3) THEN ! use irecv/isend
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

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_l3d"
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

   IF (iorder_sendrecv == 1 .OR. iorder_sendrecv >= 3) THEN
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
   ELSE IF (iorder_sendrecv >= 3) THEN ! use irecv/isend
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

   IF (activate_sync_timers) CALL timer_stop(timer_exch_data)

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
SUBROUTINE exchange_data_mult(p_pat, nfields, ndim2tot, recv1, send1, add1, recv2, send2,   &
                              add2, recv3, send3, add3, recv4, send4, add4, recv5, send5,   &
                              add5, recv6, send6, add6, recv7, send7, add7, recv4d, send4d, &
                              add4d, nshift)

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
   INTEGER, OPTIONAL, INTENT(IN) :: nshift

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_mult"
   TYPE(t_fieldptr) :: recv(nfields), send(nfields), add(nfields)
   INTEGER        :: ndim2(nfields), noffset(nfields)

   REAL(wp) :: send_buf(ndim2tot,p_pat%n_send),recv_buf(ndim2tot,p_pat%n_recv)

   INTEGER :: i, k, kshift(nfields), jb,ik, jl, n, np, irs, iss, pid, icount, nf4d
   LOGICAL :: lsend, ladd

!-----------------------------------------------------------------------

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL work_mpi_barrier()
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   IF (activate_sync_timers) CALL timer_start(timer_exch_data)

   lsend     = .FALSE.
   ladd      = .FALSE.

   IF (PRESENT(nshift)) THEN
     kshift = nshift
   ELSE
     kshift = 0
   ENDIF

   IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) .AND. &
     & .NOT. my_process_is_mpi_seq()) THEN
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
     nf4d = SIZE(recv4d,4)
     DO i = 1, nf4d
       recv(i)%fld => recv4d(:,:,:,i)
     ENDDO
     IF (PRESENT(send4d)) THEN ! all 4D fields must have the same dimensions
       DO i = 1, nf4d
         send(i)%fld => send4d(:,:,:,i)
       ENDDO
       lsend = .TRUE.
     ENDIF
     IF (PRESENT(add4d)) THEN
       DO i = 1, nf4d
         add(i)%fld => add4d(:,:,:,i)
       ENDDO
       ladd = .TRUE.
     ENDIF
   ELSE
     nf4d = 0
   ENDIF
   IF (PRESENT(recv1)) THEN
     recv(nf4d+1)%fld => recv1
     IF (PRESENT(send1)) THEN
       send(nf4d+1)%fld => send1
       lsend = .TRUE.
     ENDIF
     IF (PRESENT(add1)) THEN
       add(nf4d+1)%fld => add1
       ladd = .TRUE.
     ENDIF
     IF (PRESENT(recv2)) THEN
       recv(nf4d+2)%fld => recv2
       IF (lsend) send(nf4d+2)%fld => send2
       IF (ladd)  add(nf4d+2)%fld  => add2
       IF (PRESENT(recv3)) THEN
         recv(nf4d+3)%fld => recv3
         IF (lsend) send(nf4d+3)%fld => send3
         IF (ladd)  add(nf4d+3)%fld  => add3
         IF (PRESENT(recv4)) THEN
           recv(nf4d+4)%fld => recv4
           IF (lsend) send(nf4d+4)%fld => send4
           IF (ladd)  add(nf4d+4)%fld  => add4
           IF (PRESENT(recv5)) THEN
             recv(nf4d+5)%fld => recv5
             IF (lsend) send(nf4d+5)%fld => send5
             IF (ladd)  add(nf4d+5)%fld  => add5
             IF (PRESENT(recv6)) THEN
               recv(nf4d+6)%fld => recv6
               IF (lsend) send(nf4d+6)%fld => send6
               IF (ladd)  add(nf4d+6)%fld  => add6
               IF (PRESENT(recv7)) THEN
                 recv(nf4d+7)%fld => recv7
                 IF (lsend) send(nf4d+7)%fld => send7
                 IF (ladd)  add(nf4d+7)%fld  => add7
               ENDIF
             ENDIF
           ENDIF
         ENDIF
       ENDIF
     ENDIF
   ENDIF

   IF(my_process_is_mpi_seq()) THEN
     DO n = 1, nfields
       IF(lsend) THEN
         IF(ladd) THEN
           CALL exchange_data_seq(p_pat, recv(n)%fld(:,:,:), send(n)%fld(:,:,:), add(n)%fld(:,:,:))
         ELSE
           CALL exchange_data_seq(p_pat, recv(n)%fld(:,:,:), send(n)%fld(:,:,:))
         ENDIF
       ELSE
         IF(ladd) THEN
           CALL exchange_data_seq(p_pat, recv(n)%fld(:,:,:), add=add(n)%fld(:,:,:))
         ELSE
           CALL exchange_data_seq(p_pat, recv(n)%fld(:,:,:))
         ENDIF
       ENDIF
     ENDDO

     IF (activate_sync_timers) CALL timer_stop(timer_exch_data)
     RETURN
     !---------------------------------------------------------------
   ENDIF

   ! Reset kshift to 0 if 2D fields are passed together with 3D fields
   DO n = 1, nfields
     IF (SIZE(recv(n)%fld,2) == 1) kshift(n) = 0
   ENDDO

   noffset(1) = 0
   ndim2(1)   = SIZE(recv(1)%fld,2) - kshift(1)
   DO n = 2, nfields
     noffset(n) = noffset(n-1)+ndim2(n-1)
     ndim2(n)   = SIZE(recv(n)%fld,2) - kshift(n)
   ENDDO

   ! Set up send buffer
#ifdef __SX__
   IF ( lsend ) THEN
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_send
           send_buf(k+noffset(n),i) = &
             send(n)%fld(p_pat%send_src_idx(i),k+kshift(n),p_pat%send_src_blk(i))
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
             recv(n)%fld(p_pat%send_src_idx(i),k+kshift(n),p_pat%send_src_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ENDIF
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl)
#endif
   DO i = 1, p_pat%n_send
     jb = p_pat%send_src_blk(i)
     jl = p_pat%send_src_idx(i)
     IF ( lsend ) THEN
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),i) = send(n)%fld(jl,k+kshift(n),jb)
         ENDDO
       ENDDO
     ELSE
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           send_buf(k+noffset(n),i) = recv(n)%fld(jl,k+kshift(n),jb)
         ENDDO
       ENDDO
     ENDIF
   ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
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

   IF (iorder_sendrecv > 0) THEN
     ! Wait for all outstanding requests to finish
     IF (activate_sync_timers) CALL timer_start(timer_exch_data_wait)
     CALL p_wait
     IF (activate_sync_timers) CALL timer_stop(timer_exch_data_wait)

   ELSE IF (iorder_sendrecv == 0) THEN
     ! dummy mode; just fill the receive buffer with "something reasonable"
     IF (p_pat%n_send > 0) THEN
       DO k = 1, ndim2tot
         recv_buf(k,:) = send_buf(k,1)
       ENDDO
     ELSE
       DO k = 1, ndim2tot
         recv_buf(k,:) = 1._wp
       ENDDO
     ENDIF
   ENDIF

   IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL work_mpi_barrier()
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   ! Fill in receive buffer

#ifdef __SX__
   IF (ladd) THEN
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           recv(n)%fld(p_pat%recv_dst_idx(i),k+kshift(n),p_pat%recv_dst_blk(i)) =  &
             recv_buf(k+noffset(n),p_pat%recv_src(i)) +                         &
             add(n)%fld(p_pat%recv_dst_idx(i),k+kshift(n),p_pat%recv_dst_blk(i))
         ENDDO
       ENDDO
     ENDDO
   ELSE
     DO n = 1, nfields
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat%n_pnts
           recv(n)%fld(p_pat%recv_dst_idx(i),k+kshift(n),p_pat%recv_dst_blk(i)) =  &
             recv_buf(k+noffset(n),p_pat%recv_src(i))
         ENDDO
       ENDDO
     ENDDO
   ENDIF
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,ik)
#endif
   DO i = 1, p_pat%n_pnts
     jb = p_pat%recv_dst_blk(i)
     jl = p_pat%recv_dst_idx(i)
     ik  = p_pat%recv_src(i)
     IF (ladd) THEN
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           recv(n)%fld(jl,k+kshift(n),jb)= recv_buf(k+noffset(n),ik)+add(n)%fld(jl,k+kshift(n),jb)
         ENDDO
       ENDDO
     ELSE
       DO n = 1, nfields
         DO k = 1, ndim2(n)
           recv(n)%fld(jl,k+kshift(n),jb) = recv_buf(k+noffset(n),ik)
         ENDDO
       ENDDO
     ENDIF
   ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
     
   IF (activate_sync_timers) CALL timer_stop(timer_exch_data)

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

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_4de3"
   INTEGER, INTENT(IN)           :: nfields, ndim2tot

   INTEGER :: ndim2, noffset

   REAL(wp) :: send_buf(ndim2tot,p_pat%n_send),recv_buf(ndim2tot,p_pat%n_recv)

   INTEGER :: i, k, ik, jb, jl, n, np, irs, iss, pid, icount
   LOGICAL :: lsend

!-----------------------------------------------------------------------

   IF(my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Not yet implemented!")
   END IF

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL work_mpi_barrier()
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

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
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,noffset)
#endif
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
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
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

   IF (iorder_sendrecv > 0) THEN
     ! Wait for all outstanding requests to finish
     IF (activate_sync_timers) CALL timer_start(timer_exch_data_wait)
     CALL p_wait
     IF (activate_sync_timers) CALL timer_stop(timer_exch_data_wait)

   ELSE IF (iorder_sendrecv == 0) THEN
     ! dummy mode; just fill the receive buffer with "something reasonable"
     IF (p_pat%n_send > 0) THEN
       DO k = 1, ndim2tot
         recv_buf(k,:) = send_buf(k,1)
       ENDDO
     ELSE
       DO k = 1, ndim2tot
         recv_buf(k,:) = 1._wp
       ENDDO
     ENDIF
   ENDIF

   IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL work_mpi_barrier()
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

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
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,ik,noffset)
#endif
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
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL DO
#endif
#endif
   IF (activate_sync_timers) CALL timer_stop(timer_exch_data)

END SUBROUTINE exchange_data_4de3


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
SUBROUTINE complete_async_comm(p_pat, nfields, recv_buf, recv1, recv2, recv3, &
                               recv4, recv5, recv4d)

   TYPE(t_comm_pattern), INTENT(IN) :: p_pat

   REAL(wp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv4d(:,:,:,:)

   INTEGER, INTENT(IN)      :: nfields

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
!! Optimized version by Guenther Zaengl to process up to two 4D fields or up to six 3D fields
!! for an array-sized communication pattern (as needed for boundary interpolation) in one step
!!
SUBROUTINE exchange_data_grf(p_pat, nfields, ndim2tot, nsendtot, nrecvtot, recv1, send1,   &
                             recv2, send2, recv3, send3, recv4, send4, recv5, send5,       &
                             recv6, send6, recv4d1, send4d1, recv4d2, send4d2)

  TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat(:)

   REAL(wp), INTENT(INOUT), TARGET, OPTIONAL ::  &
     recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4d1(:,:,:,:), &
     recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), recv4d2(:,:,:,:)
   ! Note: the last index of the send fields corresponds to the dimension of p_pat
   ! On the other hand, they are not blocked and have the vertical index first
   REAL(wp), INTENT(IN   ), TARGET, OPTIONAL ::  &
     send1(:,:,:), send2(:,:,:), send3(:,:,:), send4d1(:,:,:,:), &
     send4(:,:,:), send5(:,:,:), send6(:,:,:), send4d2(:,:,:,:)

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_grf"
   INTEGER, INTENT(IN)           :: nfields  ! total number of input fields
   INTEGER, INTENT(IN)           :: ndim2tot ! sum of vertical levels of input fields
   INTEGER, INTENT(IN)           :: nsendtot ! total number of send points
   INTEGER, INTENT(IN)           :: nrecvtot ! total number of receive points

   TYPE t_fieldptr_recv
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr_recv

   TYPE t_fieldptr_send
     REAL(wp), POINTER :: fld(:,:)
   END TYPE t_fieldptr_send

   TYPE(t_fieldptr_recv) :: recv(nfields)
   TYPE(t_fieldptr_send) :: send(nfields*SIZE(p_pat))

   INTEGER        :: ndim2(nfields), noffset(nfields),            &
                     ioffset_s(SIZE(p_pat)), ioffset_r(SIZE(p_pat))

   REAL(wp) :: send_buf(ndim2tot,nsendtot),recv_buf(ndim2tot,nrecvtot), &
               auxs_buf(ndim2tot,nsendtot),auxr_buf(ndim2tot,nrecvtot)

   INTEGER :: i, k, ik, jb, jl, n, np, irs, ire, iss, ise, &
              npats, isum, ioffset, isum1, n4d, pid

!-----------------------------------------------------------------------

   IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL work_mpi_barrier()
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   IF (activate_sync_timers) CALL timer_start(timer_exch_data)

   npats = SIZE(p_pat)  ! Number of communication patterns provided on input

   ! Set up irecv's for receive buffers
   ! Note: the dummy mode (iorder_sendrecv=0) does not work for nest boundary communication
   ! because there may be PEs receiving but not sending data. Therefore, iorder_sendrecv=0
   ! is treated as iorder_sendrecv=1 in this routine.
   IF ((iorder_sendrecv <= 1 .OR. iorder_sendrecv >= 3) .AND. .NOT. my_process_is_mpi_seq()) THEN

     ioffset = 0
     DO np = 1, p_pat(1)%np_recv ! loop over PEs from where to receive the data

       pid = p_pat(1)%pelist_recv(np) ! ID of receiver PE

       ! Sum up receive points over all communication patterns to be processed
       isum = ioffset
       DO n = 1, npats
         isum = isum + p_pat(n)%recv_limits(pid+1) - p_pat(n)%recv_limits(pid)
       ENDDO

       IF(isum > ioffset) &
         CALL p_irecv(auxr_buf(1,ioffset+1), pid, 1, p_count=(isum-ioffset)*ndim2tot, &
                      comm=p_comm_work)
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
     ioffset_r(np) = ioffset_r(np-1) + p_pat(np-1)%n_recv
     ioffset_s(np) = ioffset_s(np-1) + p_pat(np-1)%n_send
   ENDDO


   IF (my_process_is_mpi_seq()) THEN
     DO np = 1, npats
       DO i = 1, p_pat(np)%n_pnts
         DO n = 1, nfields
           DO k = 1, ndim2(n) 
             recv(n)%fld( p_pat(np)%recv_dst_idx(i), k, p_pat(np)%recv_dst_blk(i) ) =  &
               send(np+(n-1)*npats)%fld( k, p_pat(np)%send_src_idx(p_pat(np)%recv_src(i)) )
           ENDDO
         ENDDO
       ENDDO
     ENDDO
     IF (activate_sync_timers) CALL timer_stop(timer_exch_data)

     RETURN
     !---------------------------------------------------------
   ENDIF


   ! Set up send buffer
#ifdef __SX__
   DO n = 1, nfields
     DO np = 1, npats
!CDIR UNROLL=6
       DO k = 1, ndim2(n)
         DO i = 1, p_pat(np)%n_send
           send_buf(k+noffset(n),i+ioffset_s(np)) =                &
             & send(np+(n-1)*npats)%fld(k,p_pat(np)%send_src_idx(i))
         ENDDO
       ENDDO
     ENDDO
   ENDDO
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL
#endif
   DO np = 1, npats
#ifdef __OMPPAR_COPY__
!$OMP DO PRIVATE(jl)
#endif
     DO i = 1, p_pat(np)%n_send
       jl = p_pat(np)%send_src_idx(i)
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
     ! Send our data
     ioffset = 0
     DO np = 1, p_pat(1)%np_send ! loop over PEs where to send the data

       pid = p_pat(1)%pelist_send(np) ! ID of sender PE

       ! Copy send points for all communication patterns into one common send buffer
       isum = ioffset
       DO n = 1, npats
         iss = p_pat(n)%send_limits(pid)+1 + ioffset_s(n)
         ise = p_pat(n)%send_limits(pid+1) + ioffset_s(n)
         isum1 = ise - iss + 1
         IF (isum1 > 0) THEN
!CDIR COLLAPSE
           auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
           isum = isum+isum1
         ENDIF
       ENDDO

       IF(isum > ioffset) CALL p_send(auxs_buf(1,ioffset+1), pid, 1,             &
                               p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)

       ioffset = isum

     ENDDO
   ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
     ioffset = 0
     DO np = 1, p_pat(1)%np_send ! loop over PEs where to send the data

       pid = p_pat(1)%pelist_send(np) ! ID of sender PE

       ! Copy send points for all communication patterns into one common send buffer
       isum = ioffset
       DO n = 1, npats
         iss = p_pat(n)%send_limits(pid)+1 + ioffset_s(n)
         ise = p_pat(n)%send_limits(pid+1) + ioffset_s(n)
         isum1 = ise - iss + 1
         IF (isum1 > 0) THEN
!CDIR COLLAPSE
           auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
           isum = isum+isum1
         ENDIF
       ENDDO

       IF(isum > ioffset) CALL p_isend(auxs_buf(1,ioffset+1), pid, 1,            &
                               p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)

       ioffset = isum

     ENDDO

     ioffset = 0
     DO np = 1, p_pat(1)%np_recv ! loop over PEs from where to receive the data

       pid = p_pat(1)%pelist_recv(np) ! ID of receiver PE

       ! Sum up receive points over all communication patterns to be processed
       isum = ioffset
       DO n = 1, npats
         isum = isum + p_pat(n)%recv_limits(pid+1) - p_pat(n)%recv_limits(pid)
       ENDDO

       IF(isum > ioffset) CALL p_recv(auxr_buf(1,ioffset+1), pid, 1,             &
                               p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)
       ioffset = isum

     ENDDO
   ELSE IF (iorder_sendrecv >= 3) THEN ! use isend/recv
     ioffset = 0
     DO np = 1, p_pat(1)%np_send ! loop over PEs where to send the data

       pid = p_pat(1)%pelist_send(np) ! ID of sender PE

       ! Copy send points for all communication patterns into one common send buffer
       isum = ioffset
       DO n = 1, npats
         iss = p_pat(n)%send_limits(pid)+1 + ioffset_s(n)
         ise = p_pat(n)%send_limits(pid+1) + ioffset_s(n)
         isum1 = ise - iss + 1
         IF (isum1 > 0) THEN
!CDIR COLLAPSE
           auxs_buf(:,isum+1:isum+isum1) = send_buf(:,iss:ise)
           isum = isum+isum1
         ENDIF
       ENDDO

       IF(isum > ioffset) CALL p_isend(auxs_buf(1,ioffset+1), pid, 1,            &
                               p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)

       ioffset = isum

     ENDDO
   ENDIF

   ! Wait for all outstanding requests to finish
   IF (activate_sync_timers) CALL timer_start(timer_exch_data_wait)
   CALL p_wait
   IF (activate_sync_timers) CALL timer_stop(timer_exch_data_wait)


   IF (itype_exch_barrier == 2 .OR. itype_exch_barrier == 3) THEN
     IF (activate_sync_timers) CALL timer_start(timer_barrier)
     CALL work_mpi_barrier()
     IF (activate_sync_timers) CALL timer_stop(timer_barrier)
   ENDIF

   ! Copy exchanged data back to receive buffer
   ioffset = 0
   DO np = 1, p_pat(1)%np_recv ! loop over PEs from where to receive the data

     pid = p_pat(1)%pelist_recv(np) ! ID of receiver PE

     isum = ioffset
     DO n = 1, npats
       irs = p_pat(n)%recv_limits(pid)+1 + ioffset_r(n)
       ire = p_pat(n)%recv_limits(pid+1) + ioffset_r(n)
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
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL
#endif
   DO np = 1, npats
#ifdef __OMPPAR_COPY__
!$OMP DO PRIVATE(jb,jl,ik)
#endif
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
#ifdef __OMPPAR_COPY__
!$OMP END DO
#endif
   ENDDO
#ifdef __OMPPAR_COPY__
!$OMP END PARALLEL
#endif
#endif

   IF (activate_sync_timers) CALL timer_stop(timer_exch_data)

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

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_r2d"
   REAL(wp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

   !-----------------------------------------------------------------------

   ! special treatment for trivial communication patterns of
   ! sequential runs
   IF(my_process_is_mpi_seq()) THEN
     CALL exchange_data_r2d_seq(p_pat, recv, send, add, send_lbound2,l_recv_exists)
     RETURN
   END IF

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


! SEQUENTIAL version of subroutine "exchange_data_r3d"
!
SUBROUTINE exchange_data_r2d_seq(p_pat, recv, send, add, send_lbound2, l_recv_exists)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   REAL(wp), INTENT(INOUT), TARGET        :: recv(:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
   REAL(wp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
   INTEGER, OPTIONAL :: send_lbound2
   LOGICAL, OPTIONAL :: l_recv_exists
   ! local variables
    CHARACTER(*), PARAMETER :: routine = "mo_communication:exchange_data_r2d_seq"
   INTEGER :: i, lbound2

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

   IF(PRESENT(send) .AND. PRESENT(send_lbound2)) THEN
     lbound2 = send_lbound2 - 1
   ELSE
     lbound2 = 0
   ENDIF

   IF(PRESENT(add)) THEN
     DO i=1,p_pat%n_pnts
       recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
         &  add( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )                +  &
         &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
         &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound2)
     END DO
   ELSE
     DO i=1,p_pat%n_pnts
       recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
         &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
         &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound2)
     END DO
   END IF

 END SUBROUTINE exchange_data_r2d_seq


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

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_i2d"
   INTEGER :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

   !-----------------------------------------------------------------------

   ! special treatment for trivial communication patterns of
   ! sequential runs
   IF(my_process_is_mpi_seq()) THEN
     CALL exchange_data_i2d_seq(p_pat, recv, send, add, send_lbound2)
     RETURN
   END IF

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

! SEQUENTIAL version of subroutine "exchange_data_r3d"
!
SUBROUTINE exchange_data_i2d_seq(p_pat, recv, send, add, send_lbound2)

   TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
   INTEGER, INTENT(INOUT), TARGET        :: recv(:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
   INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:)
   INTEGER, OPTIONAL :: send_lbound2
   ! local variables
    CHARACTER(*), PARAMETER :: routine = "mo_communication:exchange_data_i2d_seq"
   INTEGER :: i, lbound2

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

   IF(PRESENT(send) .AND. PRESENT(send_lbound2)) THEN
     lbound2 = send_lbound2 - 1
   ELSE
     lbound2 = 0
   ENDIF

   IF(PRESENT(add)) THEN
     DO i=1,p_pat%n_pnts
       recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
         &  add( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )                +  &
         &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
         &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound2)
     END DO
   ELSE
     DO i=1,p_pat%n_pnts
       recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
         &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
         &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound2)
     END DO
   END IF

 END SUBROUTINE exchange_data_i2d_seq

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

   CHARACTER(len=*), PARAMETER :: routine = "mo_communication::exchange_data_l2d"
   LOGICAL :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

   !-----------------------------------------------------------------------

   IF(my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Not yet implemented!")
   END IF


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


!> In-situ reordering of communication gather pattern.
!
!  @author M. Hanke, DKRZ (2013-11-13)
!
SUBROUTINE reorder_comm_gather_pattern(gather_pattern, idx_old2new)
  TYPE (t_comm_gather_pattern), INTENT(INOUT) :: gather_pattern
  INTEGER,                      INTENT(IN)    :: idx_old2new(:) ! permutation array

  gather_pattern%loc_index(:) = idx_old2new(gather_pattern%loc_index(:))

END SUBROUTINE reorder_comm_gather_pattern



!> In-situ reordering of communication pattern.
!
!  @author F. Prill, DWD (2013-07-30)
!
SUBROUTINE reorder_comm_pattern(comm_pat, idx_old2new)
  TYPE (t_comm_pattern), INTENT(INOUT) :: comm_pat
  INTEGER,               INTENT(IN)    :: idx_old2new(:) ! permutation array

  CALL reorder_comm_pattern_rcv(comm_pat, idx_old2new)
  CALL reorder_comm_pattern_snd(comm_pat, idx_old2new)

END SUBROUTINE reorder_comm_pattern


!> In-situ reordering of communication pattern.
!
!  @author F. Prill, DWD (2013-07-30)
!
SUBROUTINE reorder_comm_pattern_snd(comm_pat, idx_old2new)
  TYPE (t_comm_pattern), INTENT(INOUT) :: comm_pat
  INTEGER,               INTENT(IN)    :: idx_old2new(:) ! permutation array
  ! local variables
  INTEGER :: i, iidx

  ! return, if communication pattern has not yet been initialized
  IF (.NOT. ALLOCATED(comm_pat%recv_limits))  RETURN

  ! "send_src_idx/blk(i)" : For all points i=1,n_send the data in the
  ! send buffer is copied from the source array at position
  ! send_src_idx/blk(i)

  ! translate indices, overwriting old values
  DO i=1,comm_pat%n_send
    iidx            = idx_1d(comm_pat%send_src_idx(i), comm_pat%send_src_blk(i))
    iidx            = idx_old2new(iidx)
    comm_pat%send_src_idx(i) = idx_no(iidx)
    comm_pat%send_src_blk(i) = blk_no(iidx)
  END DO
END SUBROUTINE reorder_comm_pattern_snd

!> In-situ reordering of communication pattern.
!
!  @author F. Prill, DWD (2013-07-30)
!
SUBROUTINE reorder_comm_pattern_rcv(comm_pat, idx_old2new)
  TYPE (t_comm_pattern), INTENT(INOUT) :: comm_pat
  INTEGER,               INTENT(IN)    :: idx_old2new(:) ! permutation array
  ! local variables
  INTEGER :: i, iidx

  ! return, if communication pattern has not yet been initialized
  IF (.NOT. ALLOCATED(comm_pat%recv_limits))  RETURN

  ! "recv_dst_idx/blk(i)": For all points i=1,n_pnts the data received
  ! is copied to the destination array at position recv_dst_idx/blk(i)

  ! translate indices, overwriting old values
  DO i=1,comm_pat%n_pnts
    iidx            = idx_1d(comm_pat%recv_dst_idx(i), comm_pat%recv_dst_blk(i))
    iidx            = idx_old2new(iidx)
    comm_pat%recv_dst_idx(i) = idx_no(iidx)
    comm_pat%recv_dst_blk(i) = blk_no(iidx)
  END DO

END SUBROUTINE reorder_comm_pattern_rcv

SUBROUTINE gather_r_1d_deblock(in_array, out_array, gather_pattern)
  ! dimension (nproma, nblk)
  REAL(wp), INTENT(IN) :: in_array(:,:)
  ! dimension (global length); only required on root
  REAL(wp), INTENT(INOUT) :: out_array(:)
  TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

  REAL(wp), ALLOCATABLE :: send_buffer(:,:), recv_buffer(:,:)
  INTEGER :: i, num_send_points, idx, blk

  num_send_points = SUM(gather_pattern%collector_send_size(:))
  ALLOCATE(send_buffer(1, num_send_points))

  IF (p_pe_work == process_mpi_root_id) THEN
    ALLOCATE(recv_buffer(1, SUM(gather_pattern%collector_size(:))))
  ELSE
    ALLOCATE(recv_buffer(0,0))
  END IF

  DO i = 1, SIZE(gather_pattern%loc_index(:))
    idx = idx_no(gather_pattern%loc_index(i))
    blk = blk_no(gather_pattern%loc_index(i))
    send_buffer(1,i) = in_array(idx, blk)
  END DO

  CALL two_phase_gather(send_buffer_r=send_buffer, recv_buffer_r=recv_buffer, &
    &                   gather_pattern=gather_pattern)

  IF (p_pe_work == process_mpi_root_id) &
    out_array(1:SIZE(recv_buffer, 2)) = recv_buffer(1,:)

  DEALLOCATE(send_buffer,recv_buffer)

END SUBROUTINE gather_r_1d_deblock

SUBROUTINE gather_i_1d_deblock(in_array, out_array, gather_pattern)
  ! dimension (nproma, nblk)
  INTEGER, INTENT(IN) :: in_array(:,:)
  ! dimension (global length); only required on root
  INTEGER, INTENT(INOUT) :: out_array(:)
  TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

  INTEGER, ALLOCATABLE :: send_buffer(:,:), recv_buffer(:,:)
  INTEGER :: i, num_send_points, idx, blk

  num_send_points = SUM(gather_pattern%collector_send_size(:))
  ALLOCATE(send_buffer(1, num_send_points))

  IF (p_pe_work == process_mpi_root_id) THEN
    ALLOCATE(recv_buffer(1, SUM(gather_pattern%collector_size(:))))
  ELSE
    ALLOCATE(recv_buffer(0,0))
  END IF

  DO i = 1, SIZE(gather_pattern%loc_index(:))
    idx = idx_no(gather_pattern%loc_index(i))
    blk = blk_no(gather_pattern%loc_index(i))
    send_buffer(1,i) = in_array(idx, blk)
  END DO

  CALL two_phase_gather(send_buffer_i=send_buffer, recv_buffer_i=recv_buffer, &
    &                   gather_pattern=gather_pattern)

  IF (p_pe_work == process_mpi_root_id) &
    out_array(1:SIZE(recv_buffer, 2)) = recv_buffer(1,:)

  DEALLOCATE(send_buffer,recv_buffer)

END SUBROUTINE gather_i_1d_deblock

SUBROUTINE gather_r_2d_deblock(in_array, out_array, gather_pattern)
  ! dimension (nproma, nlev, nblk)
  REAL(wp), INTENT(IN) :: in_array(:,:,:)
  ! dimension (global length, nlev); only required on root
  REAL(wp), INTENT(INOUT) :: out_array(:,:)
  TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

  REAL(wp), ALLOCATABLE :: send_buffer(:,:), recv_buffer(:,:)
  INTEGER :: i, num_send_points, nlev, idx, blk

  nlev = SIZE(in_array, 2)

  IF (SIZE(in_array, 1) /= nproma) &
    CALL finish("gather_r_2d_deblock", &
      &         "size of first dimension of in_array is not nproma")

  IF (nlev /= SIZE(out_array, 2) .AND. p_pe_work == process_mpi_root_id) &
    CALL finish("gather_r_2d_deblock", &
      &         "second size of in_array and out_array are not the same")

  num_send_points = SUM(gather_pattern%collector_send_size(:))
  IF (SIZE(in_array, 1) * SIZE(in_array, 3) < num_send_points) &
    CALL finish("gather_r_2d_deblock", "in_array is too small")

  ALLOCATE(send_buffer(nlev, num_send_points))
  IF (p_pe_work == process_mpi_root_id) THEN
    ALLOCATE(recv_buffer(nlev, SUM(gather_pattern%collector_size(:))))
  ELSE
    ALLOCATE(recv_buffer(0,0))
  END IF

  DO i = 1, SIZE(gather_pattern%loc_index(:))
    idx = idx_no(gather_pattern%loc_index(i))
    blk = blk_no(gather_pattern%loc_index(i))
    send_buffer(:,i) = in_array(idx, :, blk)
  END DO

  CALL two_phase_gather(send_buffer_r=send_buffer, recv_buffer_r=recv_buffer, &
    &                   gather_pattern=gather_pattern)

  IF (p_pe_work == process_mpi_root_id) &
    out_array(1:SIZE(recv_buffer, 2),1:SIZE(recv_buffer, 1)) = &
      TRANSPOSE(recv_buffer(:,:))
END SUBROUTINE gather_r_2d_deblock

SUBROUTINE gather_i_2d_deblock(in_array, out_array, gather_pattern)
  ! dimension (nproma, nlev, nblk)
  INTEGER, INTENT(IN) :: in_array(:,:,:)
  ! dimension (global length, nlev); only required on root
  INTEGER, INTENT(INOUT) :: out_array(:,:)
  TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

  INTEGER, ALLOCATABLE :: send_buffer(:,:), recv_buffer(:,:)
  INTEGER :: i, num_send_points, nlev, idx, blk

  nlev = SIZE(in_array, 2)

  IF (SIZE(in_array, 1) /= nproma) &
    CALL finish("gather_i_2d_deblock", &
      &         "size of first dimension of in_array is not nproma")

  IF (nlev /= SIZE(out_array, 2) .AND. p_pe_work == process_mpi_root_id) &
    CALL finish("gather_i_2d_deblock", &
      &         "second size of in_array and out_array are not the same")

  num_send_points = SUM(gather_pattern%collector_send_size(:))
  IF (SIZE(in_array, 1) * SIZE(in_array, 3) < num_send_points) &
    CALL finish("gather_r_2d_deblock", "in_array is too small")

  ALLOCATE(send_buffer(nlev, num_send_points))
  IF (p_pe_work == process_mpi_root_id) THEN
    ALLOCATE(recv_buffer(nlev, SUM(gather_pattern%collector_size(:))))
  ELSE
    ALLOCATE(recv_buffer(0,0))
  END IF

  DO i = 1, SIZE(gather_pattern%loc_index(:))
    idx = idx_no(gather_pattern%loc_index(i))
    blk = blk_no(gather_pattern%loc_index(i))
    send_buffer(:,i) = in_array(idx, :, blk)
  END DO

  CALL two_phase_gather(send_buffer_i=send_buffer, recv_buffer_i=recv_buffer, &
    &                   gather_pattern=gather_pattern)

  IF (p_pe_work == process_mpi_root_id) &
    out_array(1:SIZE(recv_buffer, 2),1:SIZE(recv_buffer, 1)) = &
      TRANSPOSE(recv_buffer(:,:))
END SUBROUTINE gather_i_2d_deblock

SUBROUTINE two_phase_gather(send_buffer_r, send_buffer_i, &
                            recv_buffer_r, recv_buffer_i, gather_pattern)
  ! dimension (:, length), prepared according to gather pattern
  REAL(wp), OPTIONAL, INTENT(IN) :: send_buffer_r(:,:)
  INTEGER, OPTIONAL, INTENT(IN) :: send_buffer_i(:,:)
  ! dimension (:, global length); only required on root
  REAL(wp), OPTIONAL, INTENT(INOUT) :: recv_buffer_r(:,:)
  INTEGER, OPTIONAL, INTENT(INOUT) :: recv_buffer_i(:,:)
  TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

  REAL(wp), ALLOCATABLE :: collector_buffer_r(:,:)
  INTEGER, ALLOCATABLE :: collector_buffer_i(:,:)
  INTEGER :: num_send_per_process(p_n_work), send_displ(p_n_work)
  INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
  INTEGER :: i

  IF ((PRESENT(send_buffer_r) .NEQV. PRESENT(recv_buffer_r)) .OR. &
      (PRESENT(send_buffer_i) .NEQV. PRESENT(recv_buffer_i))) &
    CALL finish("two_phase_gather", "invalid arguments")

  IF (PRESENT(send_buffer_r)) &
    ALLOCATE(collector_buffer_r(SIZE(send_buffer_r, 1), &
      &                         SUM(gather_pattern%recv_size(:))))
  IF (PRESENT(send_buffer_i)) &
    ALLOCATE(collector_buffer_i(SIZE(send_buffer_i, 1), &
      &                         SUM(gather_pattern%recv_size(:))))

  num_send_per_process(:) = 0
  num_send_per_process(gather_pattern%collector_pes(:)+1) = &
      gather_pattern%collector_send_size(:)
  num_recv_per_process(:) = 0
  num_recv_per_process(gather_pattern%recv_pes(:)+1) = &
    gather_pattern%recv_size(:)
  send_displ(1) = 0
  recv_displ(1) = 0
  DO i = 2, p_n_work
    send_displ(i) = send_displ(i-1) + num_send_per_process(i-1)
    recv_displ(i) = recv_displ(i-1) + num_recv_per_process(i-1)
  END DO

  IF (PRESENT(send_buffer_r)) &
    CALL p_alltoallv(send_buffer_r, num_send_per_process, send_displ, &
      &              collector_buffer_r, num_recv_per_process, recv_displ, &
      &              p_comm_work)
  IF (PRESENT(send_buffer_i)) &
    CALL p_alltoallv(send_buffer_i, num_send_per_process, send_displ, &
      &              collector_buffer_i, num_recv_per_process, recv_displ, &
      &              p_comm_work)

  ! reorder collector_buffer
  IF (PRESENT(send_buffer_r)) &
    collector_buffer_r(:,:) = &
      collector_buffer_r(:, gather_pattern%recv_buffer_reorder(:))
  IF (PRESENT(send_buffer_i)) &
    collector_buffer_i(:,:) = &
      collector_buffer_i(:, gather_pattern%recv_buffer_reorder(:))

  num_recv_per_process(:) = 0
  IF (p_pe_work == process_mpi_root_id) &
    num_recv_per_process(gather_pattern%collector_pes(:)+1) = &
      gather_pattern%collector_size(:)
  recv_displ(1) = 0
  DO i = 2, p_n_work
    recv_displ(i) = recv_displ(i-1) + num_recv_per_process(i-1)
  END DO

  IF (PRESENT(send_buffer_r)) &
    CALL p_gatherv(collector_buffer_r, SIZE(collector_buffer_r, 2), &
      &            recv_buffer_r, num_recv_per_process, recv_displ, &
      &            process_mpi_root_id, p_comm_work)
  IF (PRESENT(send_buffer_i)) &
    CALL p_gatherv(collector_buffer_i, SIZE(collector_buffer_i, 2), &
      &            recv_buffer_i, num_recv_per_process, recv_displ, &
      &            process_mpi_root_id, p_comm_work)

  IF (ALLOCATED(collector_buffer_r)) DEALLOCATE(collector_buffer_r)
  IF (ALLOCATED(collector_buffer_i)) DEALLOCATE(collector_buffer_i)
END SUBROUTINE two_phase_gather

END MODULE mo_communication


