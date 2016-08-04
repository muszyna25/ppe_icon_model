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
#if (defined(_CRAYFTN) && !defined(_OPENACC) )
#define __OMPPAR_COPY__
#endif

!----------------------------
#include "icon_definitions.inc"
!----------------------------
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
USE mo_impl_constants,       ONLY: SUCCESS
USE mo_scatter_pattern_base, ONLY: t_ScatterPattern, t_ScatterPatternPtr, deleteScatterPattern
USE mo_kind,                 ONLY: dp, sp
USE mo_exception,            ONLY: finish, message, message_text
USE mo_mpi,                  ONLY: p_send, p_recv, p_irecv, p_wait, p_isend,                               &
     &                             p_comm_work, my_process_is_mpi_seq, p_pe_work, p_n_work,                &
     &                             get_my_mpi_work_communicator, get_my_mpi_work_comm_size,                &
     &                             get_my_mpi_work_id, p_gather, p_gatherv, work_mpi_barrier,              &
     &                             p_alltoallv, p_alltoall, process_mpi_root_id, p_bcast,                  &
     &                             p_comm_is_intercomm, p_comm_remote_size, p_allgather, p_allgatherv,     &
     &                             MPI_COMM_NULL
USE mo_parallel_config,      ONLY: iorder_sendrecv, nproma, itype_exch_barrier
USE mo_timer,                ONLY: timer_start, timer_stop,                                                &
  &                                timer_exch_data, timer_exch_data_async, timer_barrier, timer_exch_data_wait
USE mo_run_config,           ONLY: msg_level
USE mo_decomposition_tools,  ONLY: t_glb2loc_index_lookup, get_local_index
USE mo_util_sort,            ONLY: quicksort
USE mo_util_string,          ONLY: int2string
USE mo_parallel_config,      ONLY: blk_no, idx_no, idx_1d
USE mo_fortran_tools,        ONLY: t_ptr_3d
#ifdef _OPENACC
USE mo_mpi,                  ONLY: i_am_accel_node
#endif


IMPLICIT NONE

PRIVATE

!modules interface-------------------------------------------
!subroutines
PUBLIC :: blk_no, idx_no, idx_1d
PUBLIC :: setup_comm_pattern, delete_comm_pattern, exchange_data,  &
          exchange_data_mult, exchange_data_grf,                   &
          start_async_comm, complete_async_comm,                   &
          exchange_data_4de1, exchange_data_mult_mixprec,          &
          get_np_recv, get_np_send, get_pelist_recv, exchange_data_noblk
PUBLIC :: t_comm_pattern

PUBLIC :: t_comm_gather_pattern
PUBLIC :: setup_comm_gather_pattern
PUBLIC :: delete_comm_gather_pattern

PUBLIC :: t_comm_allgather_pattern
PUBLIC :: setup_comm_allgather_pattern
PUBLIC :: delete_comm_allgather_pattern

PUBLIC :: t_ScatterPattern, t_ScatterPatternPtr, makeScatterPattern, deleteScatterPattern

PUBLIC :: ASSIGNMENT(=)
!
!variables

!--------------------------------------------------------------------------------------------------
!
TYPE t_comm_pattern

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

  PRIVATE

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
  INTEGER, ALLOCATABLE :: recv_buffer_reorder_fill(:) ! once the data is received
                                                      ! on the collectors, it has
                                                      ! to be reordered according
                                                      ! to this array (it leaves
                                                      ! holes for missing values)

  INTEGER, ALLOCATABLE :: recv_pes(:) ! ranks from which data is to be received
  INTEGER, ALLOCATABLE :: recv_size(:) ! number of remote points received
  INTEGER :: global_size ! global size of the array that is to be gathered
END TYPE t_comm_gather_pattern

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE copy_t_comm_gather_pattern
END INTERFACE

!
!------------------------------------------------------------------------------------------------
!

TYPE t_comm_allgather_pattern
  TYPE(t_comm_gather_pattern), POINTER :: gather_pattern
  INTEGER :: intercomm
END TYPE t_comm_allgather_pattern

!--------------------------------------------------------------------------------------------------
!

INTERFACE exchange_data
   MODULE PROCEDURE exchange_data_r3d
   MODULE PROCEDURE exchange_data_s3d
   MODULE PROCEDURE exchange_data_i3d
   MODULE PROCEDURE exchange_data_l3d
   MODULE PROCEDURE exchange_data_r2d
   MODULE PROCEDURE exchange_data_s2d
   MODULE PROCEDURE exchange_data_i2d
   MODULE PROCEDURE exchange_data_l2d
   MODULE PROCEDURE gather_r_2d_deblock
   MODULE PROCEDURE gather_r_1d_deblock
   MODULE PROCEDURE gather_s_1d_deblock
   MODULE PROCEDURE gather_i_2d_deblock
   MODULE PROCEDURE gather_i_1d_deblock
   MODULE PROCEDURE allgather_r_1d_deblock
   MODULE PROCEDURE allgather_i_1d_deblock
END INTERFACE

INTERFACE exchange_data_noblk
   MODULE PROCEDURE exchange_data_r1d_2d
   MODULE PROCEDURE exchange_data_s1d_2d
   MODULE PROCEDURE exchange_data_i1d_2d
END INTERFACE

INTERFACE exchange_data_seq
   MODULE PROCEDURE exchange_data_r3d_seq
   MODULE PROCEDURE exchange_data_s3d_seq
   MODULE PROCEDURE exchange_data_r2d_seq
   MODULE PROCEDURE exchange_data_s2d_seq
END INTERFACE

INTERFACE two_phase_gather_first
  MODULE PROCEDURE two_phase_gather_first_r
  MODULE PROCEDURE two_phase_gather_first_s
  MODULE PROCEDURE two_phase_gather_first_i
END INTERFACE

INTERFACE two_phase_gather_second
  MODULE PROCEDURE two_phase_gather_second_r
  MODULE PROCEDURE two_phase_gather_second_s
  MODULE PROCEDURE two_phase_gather_second_i
END INTERFACE


#if defined( _OPENACC )
#define ACC_DEBUG NO_ACC
! OpenACC Not currently enabled
#define __COMMUNICATION_NOACC
#if defined(__COMMUNICATION_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
#endif

CHARACTER(*), PARAMETER :: modname = "mo_communication"

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
    &                           send_glb2loc_index, p_pat)

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
    ENDDO

    CALL p_alltoall(num_rcv, icnt, p_pat%comm)
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
    &                                  gather_pattern, disable_consistency_check)
    INTEGER, INTENT(IN) :: global_size, owner_local(:), glb_index(:)
    TYPE(t_comm_gather_pattern), INTENT(INOUT) :: gather_pattern
    LOGICAL, INTENT(IN), OPTIONAL :: disable_consistency_check

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
    LOGICAL :: consistency_check
    CHARACTER(*), PARAMETER :: routine = modname//":setup_comm_gather_pattern"

    gather_pattern%global_size = global_size

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
    IF (ALLOCATED(gather_pattern%recv_buffer_reorder_fill)) &
      DEALLOCATE(gather_pattern%recv_buffer_reorder_fill)
    ALLOCATE(gather_pattern%recv_buffer_reorder_fill(num_recv_points))
    gather_pattern%recv_buffer_reorder_fill(:) = &
      MOD(recv_buffer-1, num_points_per_coll)+1

    IF (PRESENT(disable_consistency_check)) THEN
      consistency_check = .NOT. disable_consistency_check
    ELSE
      consistency_check = .TRUE.
    END IF

    ! consistency check
    IF (consistency_check) THEN

      ! check whether a point is owned by multiple processes
      DO i = 2, num_recv_points
        IF (recv_buffer(i) == recv_buffer(i-1)) &
          CALL finish(routine, &
          &         "One or more points are owned by multiple processes")
      END DO

      ! check if all points have an owner
      IF (SUM(gather_pattern%collector_size(:)) /= global_size) &
        CALL finish(routine, int2string(SUM(gather_pattern%collector_size(:)))//" out of "//int2string(global_size)//&
        &" points have no owner!")
    END IF

    DEALLOCATE(send_buffer, recv_buffer)

  END SUBROUTINE setup_comm_gather_pattern


  !-------------------------------------------------------------------------


  SUBROUTINE setup_comm_allgather_pattern(gather_pattern, intercomm, &
    &                                     allgather_pattern)
    TYPE(t_comm_gather_pattern), TARGET, INTENT(INOUT) :: gather_pattern
    INTEGER, OPTIONAL, INTENT(IN) :: intercomm
    TYPE(t_comm_allgather_pattern), INTENT(INOUT) :: allgather_pattern

    IF (PRESENT(intercomm)) THEN
      ! check whether intercomm is really a intercommunicator
      IF (.NOT. p_comm_is_intercomm(intercomm)) &
        CALL finish("setup_comm_allgather_pattern", "invalid intercomm")
      allgather_pattern%intercomm = intercomm
    ELSE
      allgather_pattern%intercomm = MPI_COMM_NULL
    END IF
    allgather_pattern%gather_pattern => gather_pattern
  END SUBROUTINE setup_comm_allgather_pattern


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

    p_pat%n_recv  = 0
    p_pat%n_pnts  = 0
    p_pat%n_send  = 0
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
    IF (ALLOCATED(gather_pattern%recv_buffer_reorder_fill)) &
      DEALLOCATE(gather_pattern%recv_buffer_reorder_fill)
    IF (ALLOCATED(gather_pattern%recv_pes)) &
      DEALLOCATE(gather_pattern%recv_pes)
    IF (ALLOCATED(gather_pattern%recv_size)) &
      DEALLOCATE(gather_pattern%recv_size)
  END SUBROUTINE delete_comm_gather_pattern


  !-------------------------------------------------------------------------


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


  SUBROUTINE delete_comm_allgather_pattern(allgather_pattern)
    TYPE(t_comm_allgather_pattern), INTENT(INOUT) :: allgather_pattern

    RETURN

  END SUBROUTINE delete_comm_allgather_pattern


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
  SUBROUTINE exchange_data_r3d(p_pat, recv, send, add, send_lbound3)

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
    INTEGER, OPTIONAL :: send_lbound3

    REAL(dp) :: send_buf(SIZE(recv,2),p_pat%n_send), &
      recv_buf(SIZE(recv,2),p_pat%n_recv)

    REAL(dp), POINTER :: send_ptr(:,:,:)

    INTEGER :: i, k, np, irs, iss, pid, icount, ndim2, lbound3

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_r3d_seq(p_pat, recv, send, add, send_lbound3)
      RETURN
    END IF

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL work_mpi_barrier()
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    IF(SIZE(recv,1) /= nproma) THEN
      CALL finish('exchange_data_r3d','Illegal first dimension of data array')
    ENDIF

    ndim2 = SIZE(recv,2)

!$ACC DATA CREATE( send_buf, recv_buf ) IF ( i_am_accel_node .AND. acc_on )

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
!$ACC PARALLEL PRESENT( p_pat, send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i)-lbound3+1)
      ENDDO
!$ACC END PARALLEL
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!CDIR UNROLL=6
!$ACC PARALLEL PRESENT( p_pat, send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO k = 1, ndim2
!$ACC LOOP VECTOR
        DO i = 1, p_pat%n_send
          send_buf(k,i) = send_ptr(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)-lbound3+1)
        ENDDO
      ENDDO
!$ACC END PARALLEL
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

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

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
      CALL work_mpi_barrier()
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer

    IF(PRESENT(add)) THEN
      IF (ndim2 == 1) THEN
        k = 1
!$ACC PARALLEL PRESENT( p_pat, add, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
            recv_buf(k,p_pat%recv_src(i)) +                       &
            add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( p_pat, recv_buf, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, p_pat%n_pnts
            recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
              recv_buf(k,p_pat%recv_src(i)) +                       &
              add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
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
!$ACC PARALLEL PRESENT( p_pat, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
            recv_buf(k,p_pat%recv_src(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( p_pat, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, p_pat%n_pnts
            recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
              recv_buf(k,p_pat%recv_src(i))
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
  SUBROUTINE exchange_data_s3d(p_pat, recv, send, add, send_lbound3)

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
    INTEGER, OPTIONAL :: send_lbound3

    REAL(sp) :: send_buf(SIZE(recv,2),p_pat%n_send), &
      recv_buf(SIZE(recv,2),p_pat%n_recv)

    REAL(sp), POINTER :: send_ptr(:,:,:)

    INTEGER :: i, k, np, irs, iss, pid, icount, ndim2, lbound3

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_s3d_seq(p_pat, recv, send, add, send_lbound3)
      RETURN
    END IF

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL work_mpi_barrier()
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    IF(SIZE(recv,1) /= nproma) THEN
      CALL finish('exchange_data_s3d','Illegal first dimension of data array')
    ENDIF

    ndim2 = SIZE(recv,2)

!$ACC DATA CREATE( send_buf, recv_buf ) IF ( i_am_accel_node .AND. acc_on )

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
!$ACC PARALLEL PRESENT( p_pat, send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i)-lbound3+1)
      ENDDO
!$ACC END PARALLEL
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!CDIR UNROLL=6
!$ACC PARALLEL PRESENT( p_pat, send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO k = 1, ndim2
!$ACC LOOP VECTOR
        DO i = 1, p_pat%n_send
          send_buf(k,i) = send_ptr(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)-lbound3+1)
        ENDDO
      ENDDO
!$ACC END PARALLEL
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

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

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
      CALL work_mpi_barrier()
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer

    IF(PRESENT(add)) THEN
      IF (ndim2 == 1) THEN
        k = 1
!$ACC PARALLEL PRESENT( p_pat, add, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
            recv_buf(k,p_pat%recv_src(i)) +                       &
            add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( p_pat, recv_buf, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, p_pat%n_pnts
            recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
              recv_buf(k,p_pat%recv_src(i)) +                       &
              add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
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
!$ACC PARALLEL PRESENT( p_pat, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
            recv_buf(k,p_pat%recv_src(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( p_pat, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, p_pat%n_pnts
            recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
              recv_buf(k,p_pat%recv_src(i))
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
  SUBROUTINE exchange_data_r3d_seq(p_pat, recv, send, add, send_lbound3)

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
    INTEGER, OPTIONAL :: send_lbound3
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_r3d_seq"
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
!$ACC KERNELS PRESENT( p_pat, add, send, recv ), IF( i_am_accel_node .AND. acc_on )
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
          &  add( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )                +  &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound3)
      END DO
!$ACC END KERNELS
    ELSE
!$ACC KERNELS PRESENT( p_pat, send, recv ), IF( i_am_accel_node .AND. acc_on )
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound3)
      END DO
!$ACC END KERNELS
    END IF

  END SUBROUTINE exchange_data_r3d_seq


  ! SEQUENTIAL version of subroutine "exchange_data_s3d"
  !
  SUBROUTINE exchange_data_s3d_seq(p_pat, recv, send, add, send_lbound3)

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
    INTEGER, OPTIONAL :: send_lbound3
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_s3d_seq"
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
!$ACC KERNELS PRESENT( p_pat, add, send, recv ), IF( i_am_accel_node .AND. acc_on )
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
          &  add( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )                +  &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound3)
      END DO
!$ACC END KERNELS
    ELSE
!$ACC KERNELS PRESENT( p_pat, send, recv ), IF( i_am_accel_node .AND. acc_on )
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), 1:ndim2, p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       1:ndim2,                                                                  &
          &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound3)
      END DO
!$ACC END KERNELS
    END IF

  END SUBROUTINE exchange_data_s3d_seq


  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_i3d(p_pat, recv, send, add, send_lbound3)

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
    INTEGER, INTENT(INOUT), TARGET        :: recv(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
    INTEGER, OPTIONAL :: send_lbound3

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_i3d"
    INTEGER :: send_buf(SIZE(recv,2),p_pat%n_send), &
      recv_buf(SIZE(recv,2),p_pat%n_recv)

    INTEGER, POINTER :: send_ptr(:,:,:)

    INTEGER :: i, k, np, irs, iss, pid, icount, ndim2, lbound3

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

!$ACC DATA CREATE( send_buf, recv_buf ) IF ( i_am_accel_node .AND. acc_on )

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
!$ACC PARALLEL PRESENT( p_pat, send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i)-lbound3+1)
      ENDDO
!$ACC END PARALLEL
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( p_pat, send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
      DO k = 1, ndim2
!$ACC LOOP VECTOR
        DO i = 1, p_pat%n_send
          send_buf(k,i) = send_ptr(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)-lbound3+1)
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
      DO i = 1, p_pat%n_send
        send_buf(1:ndim2,i) = send_ptr(p_pat%send_src_idx(i),1:ndim2,   &
          &                            p_pat%send_src_blk(i)-lbound3+1)
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
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    ! Fill in receive buffer

    IF(PRESENT(add)) THEN
      IF (ndim2 == 1) THEN
        k = 1
!$ACC PARALLEL PRESENT( p_pat, recv_buf, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
            recv_buf(k,p_pat%recv_src(i)) +                       &
            add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( p_pat, recv_buf, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, p_pat%n_pnts
            recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
              recv_buf(k,p_pat%recv_src(i)) +                       &
              add(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i))
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
!$ACC PARALLEL PRESENT( p_pat, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
            recv_buf(k,p_pat%recv_src(i))
        ENDDO
!$ACC END PARALLEL
      ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( p_pat, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2
!$ACC LOOP VECTOR
          DO i = 1, p_pat%n_pnts
            recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
              recv_buf(k,p_pat%recv_src(i))
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


  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_l3d(p_pat, recv, send, send_lbound3)

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
    LOGICAL, INTENT(INOUT), TARGET        :: recv(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    INTEGER, OPTIONAL :: send_lbound3

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_l3d"
    LOGICAL :: send_buf(SIZE(recv,2),p_pat%n_send), &
      recv_buf(SIZE(recv,2),p_pat%n_recv)

    LOGICAL, POINTER :: send_ptr(:,:,:)

    INTEGER :: i, k, np, irs, ire, iss, ise, ndim2, lbound3

    !-----------------------------------------------------------------------
    start_sync_timer(timer_exch_data)

    IF(my_process_is_mpi_seq()) &
      CALL finish('exchange_data','must not be called on single PE/test PE')

    IF(SIZE(recv,1) /= nproma) THEN
      CALL finish('exchange_data','Illegal first dimension of data array')
    ENDIF

    ndim2 = SIZE(recv,2)

!$ACC DATA CREATE( send_buf, recv_buf ) IF ( i_am_accel_node .AND. acc_on )

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
!$ACC PARALLEL PRESENT( p_pat, send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i = 1, p_pat%n_send
        send_buf(1,i) = send_ptr(p_pat%send_src_idx(i),1,p_pat%send_src_blk(i)-lbound3+1)
      ENDDO
!$ACC END PARALLEL
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( p_pat, send_ptr, send_buf ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
      DO k = 1, ndim2
!$ACC LOOP VECTOR
        DO i = 1, p_pat%n_send
          send_buf(k,i) = send_ptr(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i)-lbound3+1)
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
      DO i = 1, p_pat%n_send
        send_buf(1:ndim2,i) = send_ptr(p_pat%send_src_idx(i),1:ndim2,   &
          &                            p_pat%send_src_blk(i)-lbound3+1)
      ENDDO
#endif
    ENDIF

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

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
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    ! Fill in receive buffer

    IF (ndim2 == 1) THEN
      k = 1
!$ACC PARALLEL PRESENT( p_pat, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i = 1, p_pat%n_pnts
        recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
          recv_buf(k,p_pat%recv_src(i))
      ENDDO
!$ACC END PARALLEL
    ELSE
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL PRESENT( p_pat, recv_buf, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
      DO k = 1, ndim2
!$ACC LOOP VECTOR
        DO i = 1, p_pat%n_pnts
          recv(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) = &
            recv_buf(k,p_pat%recv_src(i))
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
  SUBROUTINE exchange_data_mult(p_pat, nfields, ndim2tot, recv1, send1, recv2, send2,   &
    recv3, send3, recv4, send4,  recv5, send5, recv6, send6, recv7, send7,              &
    recv4d, send4d, nshift, recv3d_arr, send3d_arr)

    TYPE(t_comm_pattern), INTENT(IN) :: p_pat

    REAL(dp), INTENT(INOUT), TARGET, OPTIONAL ::  &
      recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), &
      recv7(:,:,:), recv4d(:,:,:,:)
    REAL(dp), INTENT(IN   ), TARGET, OPTIONAL ::  &
      send1(:,:,:), send2(:,:,:), send3(:,:,:), send4(:,:,:), send5(:,:,:), send6(:,:,:), &
      send7(:,:,:), send4d(:,:,:,:)

    INTEGER, INTENT(IN)           :: nfields, ndim2tot
    TYPE(t_ptr_3d), INTENT(   IN), TARGET, OPTIONAL :: recv3d_arr(:)
    TYPE(t_ptr_3d), INTENT(INOUT), TARGET, OPTIONAL :: send3d_arr(:)
    INTEGER, OPTIONAL, INTENT(IN) :: nshift

    TYPE t_fieldptr
      REAL(dp), POINTER :: fld(:,:,:)
    END TYPE t_fieldptr

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_mult"
    TYPE(t_fieldptr) :: recv(nfields), send(nfields)
    INTEGER        :: ndim2(nfields), noffset(nfields)

    REAL(dp) :: send_buf(ndim2tot,p_pat%n_send),recv_buf(ndim2tot,p_pat%n_recv)
    REAL(dp), POINTER :: send_fld(:,:,:), recv_fld(:,:,:)  ! Refactoring for OpenACC

    INTEGER :: i, k, kshift(nfields), jb,ik, jl, n, np, irs, iss, pid, icount, nf4d
    LOGICAL :: lsend

    !-----------------------------------------------------------------------

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL work_mpi_barrier()
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    lsend     = .FALSE.

    IF (PRESENT(nshift)) THEN
      kshift = nshift
    ELSE
      kshift = 0
    ENDIF

!$ACC DATA CREATE( send_buf, recv_buf ) IF ( i_am_accel_node .AND. acc_on )

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
    ELSE
      nf4d = 0
    ENDIF


    ! Set pointers to input fields
    IF (PRESENT(recv3d_arr)) THEN
      DO i = 1, SIZE(recv3d_arr)
        recv(i+nf4d)%fld => recv3d_arr(i)%p
      ENDDO
      IF (PRESENT(send3d_arr)) THEN ! all 4D fields must have the same dimensions
        DO i = 1, SIZE(recv3d_arr)
          send(i+nf4d)%fld => send3d_arr(i)%p
        ENDDO
        lsend = .TRUE.
      ENDIF
      nf4d = nf4d + SIZE(recv3d_arr)
    ENDIF


    IF (PRESENT(recv1)) THEN
      recv(nf4d+1)%fld => recv1
      IF (PRESENT(send1)) THEN
        send(nf4d+1)%fld => send1
        lsend = .TRUE.
      ENDIF
      IF (PRESENT(recv2)) THEN
        recv(nf4d+2)%fld => recv2
        IF (lsend) send(nf4d+2)%fld => send2
        IF (PRESENT(recv3)) THEN
          recv(nf4d+3)%fld => recv3
          IF (lsend) send(nf4d+3)%fld => send3
          IF (PRESENT(recv4)) THEN
            recv(nf4d+4)%fld => recv4
            IF (lsend) send(nf4d+4)%fld => send4
            IF (PRESENT(recv5)) THEN
              recv(nf4d+5)%fld => recv5
              IF (lsend) send(nf4d+5)%fld => send5
              IF (PRESENT(recv6)) THEN
                recv(nf4d+6)%fld => recv6
                IF (lsend) send(nf4d+6)%fld => send6
                IF (PRESENT(recv7)) THEN
                  recv(nf4d+7)%fld => recv7
                  IF (lsend) send(nf4d+7)%fld => send7
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
          CALL exchange_data_seq(p_pat, recv(n)%fld(:,:,:), send(n)%fld(:,:,:))
        ELSE
          CALL exchange_data_seq(p_pat, recv(n)%fld(:,:,:))
        ENDIF
      ENDDO

    ELSE          ! WS: Removed RETURN in order to properly support OpenACC DATA region

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
#if defined( __SX__ ) || defined( _OPENACC )
      IF ( lsend ) THEN
        DO n = 1, nfields
          send_fld => send(n)%fld   ! Refactoring for OpenACC
!$ACC PARALLEL PRESENT( p_pat, send_fld, send_buf ), COPYIN( kshift, noffset, ndim2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
          DO k = 1, ndim2(n)
!$ACC LOOP VECTOR
            DO i = 1, p_pat%n_send
              send_buf(k+noffset(n),i) = &
                send(n)%fld(p_pat%send_src_idx(i),k+kshift(n),p_pat%send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
      ELSE
        ! Send and receive arrays are identical (for boundary exchange)
        DO n = 1, nfields
          recv_fld => recv(n)%fld
!$ACC PARALLEL PRESENT( p_pat, recv_fld, send_buf ), COPYIN( kshift, noffset, ndim2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
          DO k = 1, ndim2(n)
!$ACC LOOP VECTOR
            DO i = 1, p_pat%n_send
              send_buf(k+noffset(n),i) = &
                recv(n)%fld(p_pat%send_src_idx(i),k+kshift(n),p_pat%send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
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

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
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
        CALL work_mpi_barrier()
        stop_sync_timer(timer_barrier)
      ENDIF

      ! Fill in receive buffer

#if defined( __SX__ ) || defined( _OPENACC )
      DO n = 1, nfields
        recv_fld => recv(n)%fld
!$ACC PARALLEL &
!$ACC PRESENT( p_pat, recv_buf, recv_fld ), COPYIN( kshift, noffset, ndim2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2(n)
!$ACC LOOP VECTOR
          DO i = 1, p_pat%n_pnts
            recv(n)%fld(p_pat%recv_dst_idx(i),k+kshift(n),p_pat%recv_dst_blk(i)) =  &
              recv_buf(k+noffset(n),p_pat%recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,ik)
#endif
      DO i = 1, p_pat%n_pnts
        jb = p_pat%recv_dst_blk(i)
        jl = p_pat%recv_dst_idx(i)
        ik  = p_pat%recv_src(i)
        DO n = 1, nfields
          DO k = 1, ndim2(n)
            recv(n)%fld(jl,k+kshift(n),jb) = recv_buf(k+noffset(n),ik)
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
  SUBROUTINE exchange_data_mult_mixprec(p_pat, nfields_dp, ndim2tot_dp, nfields_sp, ndim2tot_sp,        &
    recv1_dp, send1_dp, recv2_dp, send2_dp, recv3_dp, send3_dp, recv4_dp, send4_dp, recv5_dp, send5_dp, &
    recv1_sp, send1_sp, recv2_sp, send2_sp, recv3_sp, send3_sp, recv4_sp, send4_sp, recv5_sp, send5_sp, &
    recv4d_dp, send4d_dp, recv4d_sp, send4d_sp, nshift)

    TYPE(t_comm_pattern), INTENT(IN) :: p_pat

    REAL(dp), INTENT(INOUT), TARGET, OPTIONAL ::  &
      recv1_dp(:,:,:), recv2_dp(:,:,:), recv3_dp(:,:,:), recv4_dp(:,:,:), recv5_dp(:,:,:), recv4d_dp(:,:,:,:)
    REAL(dp), INTENT(IN   ), TARGET, OPTIONAL ::  &
      send1_dp(:,:,:), send2_dp(:,:,:), send3_dp(:,:,:), send4_dp(:,:,:), send5_dp(:,:,:), send4d_dp(:,:,:,:)

    REAL(sp), INTENT(INOUT), TARGET, OPTIONAL ::  &
      recv1_sp(:,:,:), recv2_sp(:,:,:), recv3_sp(:,:,:), recv4_sp(:,:,:), recv5_sp(:,:,:), recv4d_sp(:,:,:,:)
    REAL(sp), INTENT(IN   ), TARGET, OPTIONAL ::  &
      send1_sp(:,:,:), send2_sp(:,:,:), send3_sp(:,:,:), send4_sp(:,:,:), send5_sp(:,:,:), send4d_sp(:,:,:,:)

    INTEGER, INTENT(IN)           :: nfields_dp, ndim2tot_dp, nfields_sp, ndim2tot_sp
    INTEGER, OPTIONAL, INTENT(IN) :: nshift

    TYPE t_fieldptr_dp
      REAL(dp), POINTER :: fld(:,:,:)
    END TYPE t_fieldptr_dp
    TYPE t_fieldptr_sp
      REAL(sp), POINTER :: fld(:,:,:)
    END TYPE t_fieldptr_sp

    TYPE(t_fieldptr_dp) :: recv_dp(nfields_dp), send_dp(nfields_dp)
    TYPE(t_fieldptr_sp) :: recv_sp(nfields_sp), send_sp(nfields_sp)
    INTEGER             :: ndim2_dp(nfields_dp), noffset_dp(nfields_dp), &
                           ndim2_sp(nfields_sp), noffset_sp(nfields_sp)

    REAL(dp) :: send_buf_dp(ndim2tot_dp,p_pat%n_send),recv_buf_dp(ndim2tot_dp,p_pat%n_recv)
    REAL(dp), POINTER :: send_fld_dp(:,:,:), recv_fld_dp(:,:,:)  ! Refactoring for OpenACC
    REAL(sp) :: send_buf_sp(ndim2tot_sp,p_pat%n_send),recv_buf_sp(ndim2tot_sp,p_pat%n_recv)
    REAL(sp), POINTER :: send_fld_sp(:,:,:), recv_fld_sp(:,:,:)  ! Refactoring for OpenACC

    INTEGER :: i, k, kshift_dp(nfields_dp), kshift_sp(nfields_sp), jb, ik, jl, n, np, irs, iss, pid, &
      icount, nf4d_dp, nf4d_sp
    LOGICAL :: lsend

    !-----------------------------------------------------------------------

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL work_mpi_barrier()
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    lsend     = .FALSE.

    IF (PRESENT(nshift)) THEN
      kshift_dp = nshift
      kshift_sp = nshift
    ELSE
      kshift_dp = 0
      kshift_sp = 0
    ENDIF

!$ACC DATA CREATE( send_buf_dp, recv_buf_dp, send_buf_sp, recv_buf_sp ) IF ( i_am_accel_node .AND. acc_on )

    IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3) .AND. &
      & .NOT. my_process_is_mpi_seq()) THEN
      ! Set up irecv's for receive buffers
      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2tot_dp
        IF (icount>0) CALL p_irecv(recv_buf_dp(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

        icount = p_pat%recv_count(np)*ndim2tot_sp
        IF (icount>0) CALL p_irecv(recv_buf_sp(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

      ENDDO
    ENDIF

    ! Set pointers to input fields
    IF (PRESENT(recv4d_dp)) THEN
      nf4d_dp = SIZE(recv4d_dp,4)
      DO i = 1, nf4d_dp
        recv_dp(i)%fld => recv4d_dp(:,:,:,i)
      ENDDO
      IF (PRESENT(send4d_dp)) THEN ! all 4D fields must have the same dimensions
        DO i = 1, nf4d_dp
          send_dp(i)%fld => send4d_dp(:,:,:,i)
        ENDDO
        lsend = .TRUE.
      ENDIF
    ELSE
      nf4d_dp = 0
    ENDIF
    IF (PRESENT(recv4d_sp)) THEN
      nf4d_sp = SIZE(recv4d_sp,4)
      DO i = 1, nf4d_sp
        recv_sp(i)%fld => recv4d_sp(:,:,:,i)
      ENDDO
      IF (PRESENT(send4d_sp)) THEN ! all 4D fields must have the same dimensions
        DO i = 1, nf4d_sp
          send_sp(i)%fld => send4d_sp(:,:,:,i)
        ENDDO
        lsend = .TRUE.
      ENDIF
    ELSE
      nf4d_sp = 0
    ENDIF

    IF (PRESENT(recv1_dp)) THEN
      recv_dp(nf4d_dp+1)%fld => recv1_dp
      IF (PRESENT(send1_dp)) THEN
        send_dp(nf4d_dp+1)%fld => send1_dp
        lsend = .TRUE.
      ENDIF
      IF (PRESENT(recv2_dp)) THEN
        recv_dp(nf4d_dp+2)%fld => recv2_dp
        IF (lsend) send_dp(nf4d_dp+2)%fld => send2_dp
        IF (PRESENT(recv3_dp)) THEN
          recv_dp(nf4d_dp+3)%fld => recv3_dp
          IF (lsend) send_dp(nf4d_dp+3)%fld => send3_dp
          IF (PRESENT(recv4_dp)) THEN
            recv_dp(nf4d_dp+4)%fld => recv4_dp
            IF (lsend) send_dp(nf4d_dp+4)%fld => send4_dp
            IF (PRESENT(recv5_dp)) THEN
              recv_dp(nf4d_dp+5)%fld => recv5_dp
              IF (lsend) send_dp(nf4d_dp+5)%fld => send5_dp
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    IF (PRESENT(recv1_sp)) THEN
      recv_sp(nf4d_sp+1)%fld => recv1_sp
      IF (PRESENT(send1_sp)) THEN
        send_sp(nf4d_sp+1)%fld => send1_sp
        lsend = .TRUE.
      ENDIF
      IF (PRESENT(recv2_sp)) THEN
        recv_sp(nf4d_sp+2)%fld => recv2_sp
        IF (lsend) send_sp(nf4d_sp+2)%fld => send2_sp
        IF (PRESENT(recv3_sp)) THEN
          recv_sp(nf4d_sp+3)%fld => recv3_sp
          IF (lsend) send_sp(nf4d_sp+3)%fld => send3_sp
          IF (PRESENT(recv4_sp)) THEN
            recv_sp(nf4d_sp+4)%fld => recv4_sp
            IF (lsend) send_sp(nf4d_sp+4)%fld => send4_sp
            IF (PRESENT(recv5_sp)) THEN
              recv_sp(nf4d_sp+5)%fld => recv5_sp
              IF (lsend) send_sp(nf4d_sp+5)%fld => send5_sp
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF


    IF(my_process_is_mpi_seq()) THEN
      DO n = 1, nfields_dp
        IF(lsend) THEN
          CALL exchange_data_seq(p_pat, recv_dp(n)%fld(:,:,:), send_dp(n)%fld(:,:,:))
        ELSE
          CALL exchange_data_seq(p_pat, recv_dp(n)%fld(:,:,:))
        ENDIF
      ENDDO
      DO n = 1, nfields_sp
        IF(lsend) THEN
          CALL exchange_data_seq(p_pat, recv_sp(n)%fld(:,:,:), send_sp(n)%fld(:,:,:))
        ELSE
          CALL exchange_data_seq(p_pat, recv_sp(n)%fld(:,:,:))
        ENDIF
      ENDDO

    ELSE          ! WS: Removed RETURN in order to properly support OpenACC DATA region

      ! Reset kshift to 0 if 2D fields are passed together with 3D fields
      DO n = 1, nfields_dp
        IF (SIZE(recv_dp(n)%fld,2) == 1) kshift_dp(n) = 0
      ENDDO
      DO n = 1, nfields_sp
        IF (SIZE(recv_sp(n)%fld,2) == 1) kshift_sp(n) = 0
      ENDDO

      IF (nfields_dp > 0) THEN
        noffset_dp(1) = 0
        ndim2_dp(1)   = SIZE(recv_dp(1)%fld,2) - kshift_dp(1)
      ENDIF
      DO n = 2, nfields_dp
        noffset_dp(n) = noffset_dp(n-1)+ndim2_dp(n-1)
        ndim2_dp(n)   = SIZE(recv_dp(n)%fld,2) - kshift_dp(n)
      ENDDO
      IF (nfields_sp > 0) THEN
        noffset_sp(1) = 0
        ndim2_sp(1)   = SIZE(recv_sp(1)%fld,2) - kshift_sp(1)
      ENDIF
      DO n = 2, nfields_sp
        noffset_sp(n) = noffset_sp(n-1)+ndim2_sp(n-1)
        ndim2_sp(n)   = SIZE(recv_sp(n)%fld,2) - kshift_sp(n)
      ENDDO


      ! Set up send buffer
#if defined( __SX__ ) || defined( _OPENACC )
      IF ( lsend ) THEN
        DO n = 1, nfields_dp
          send_fld_dp => send_dp(n)%fld   ! Refactoring for OpenACC
!$ACC PARALLEL PRESENT(p_pat,send_fld_dp,send_buf_dp), COPYIN(kshift_dp,noffset_dp,ndim2_dp), IF(i_am_accel_node .AND. acc_on)
!$ACC LOOP GANG
          DO k = 1, ndim2_dp(n)
!$ACC LOOP VECTOR
            DO i = 1, p_pat%n_send
              send_buf_dp(k+noffset_dp(n),i) = &
                send_dp(n)%fld(p_pat%send_src_idx(i),k+kshift_dp(n),p_pat%send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
        DO n = 1, nfields_sp
          send_fld_sp => send_sp(n)%fld   ! Refactoring for OpenACC
!$ACC PARALLEL PRESENT(p_pat,send_fld_sp,send_buf_sp), COPYIN(kshift_sp,noffset_sp,ndim2_sp), IF(i_am_accel_node .AND. acc_on)
!$ACC LOOP GANG
          DO k = 1, ndim2_sp(n)
!$ACC LOOP VECTOR
            DO i = 1, p_pat%n_send
              send_buf_sp(k+noffset_sp(n),i) = &
                send_sp(n)%fld(p_pat%send_src_idx(i),k+kshift_sp(n),p_pat%send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
      ELSE
        ! Send and receive arrays are identical (for boundary exchange)
        DO n = 1, nfields_dp
          recv_fld_dp => recv_dp(n)%fld
!$ACC PARALLEL PRESENT(p_pat,recv_fld_dp,send_buf_dp), COPYIN(kshift_dp,noffset_dp,ndim2_dp), IF(i_am_accel_node .AND. acc_on)
!$ACC LOOP GANG
!CDIR UNROLL=6
          DO k = 1, ndim2_dp(n)
!$ACC LOOP VECTOR
            DO i = 1, p_pat%n_send
              send_buf_dp(k+noffset_dp(n),i) = &
                recv_dp(n)%fld(p_pat%send_src_idx(i),k+kshift_dp(n),p_pat%send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
        DO n = 1, nfields_sp
          recv_fld_sp => recv_sp(n)%fld
!$ACC PARALLEL PRESENT(p_pat,recv_fld_sp,send_buf_sp), COPYIN(kshift_sp,noffset_sp,ndim2_sp), IF(i_am_accel_node .AND. acc_on)
!$ACC LOOP GANG
!CDIR UNROLL=6
          DO k = 1, ndim2_sp(n)
!$ACC LOOP VECTOR
            DO i = 1, p_pat%n_send
              send_buf_sp(k+noffset_sp(n),i) = &
                recv_sp(n)%fld(p_pat%send_src_idx(i),k+kshift_sp(n),p_pat%send_src_blk(i))
            ENDDO
          ENDDO
!$ACC END PARALLEL
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
          DO n = 1, nfields_dp
            DO k = 1, ndim2_dp(n)
              send_buf_dp(k+noffset_dp(n),i) = send_dp(n)%fld(jl,k+kshift_dp(n),jb)
            ENDDO
          ENDDO
          DO n = 1, nfields_sp
            DO k = 1, ndim2_sp(n)
              send_buf_sp(k+noffset_sp(n),i) = send_sp(n)%fld(jl,k+kshift_sp(n),jb)
            ENDDO
          ENDDO
        ELSE
          DO n = 1, nfields_dp
            DO k = 1, ndim2_dp(n)
              send_buf_dp(k+noffset_dp(n),i) = recv_dp(n)%fld(jl,k+kshift_dp(n),jb)
            ENDDO
          ENDDO
          DO n = 1, nfields_sp
            DO k = 1, ndim2_sp(n)
              send_buf_sp(k+noffset_sp(n),i) = recv_sp(n)%fld(jl,k+kshift_sp(n),jb)
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
          IF (icount>0) CALL p_send(send_buf_dp(1,iss), pid, 1, p_count=icount, comm=p_comm_work)
          icount = p_pat%send_count(np)*ndim2tot_sp
          IF (icount>0) CALL p_send(send_buf_sp(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

        ENDDO
      ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
        DO np = 1, p_pat%np_send ! loop over PEs where to send the data

          pid    = p_pat%pelist_send(np) ! ID of sender PE
          iss    = p_pat%send_startidx(np)
          icount = p_pat%send_count(np)*ndim2tot_dp
          IF (icount>0) CALL p_isend(send_buf_dp(1,iss), pid, 1, p_count=icount, comm=p_comm_work)
          icount = p_pat%send_count(np)*ndim2tot_sp
          IF (icount>0) CALL p_isend(send_buf_sp(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

        ENDDO

        DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

          pid    = p_pat%pelist_recv(np) ! ID of receiver PE
          irs    = p_pat%recv_startidx(np)
          icount = p_pat%recv_count(np)*ndim2tot_dp
          IF (icount>0) CALL p_recv(recv_buf_dp(1,irs), pid, 1, p_count=icount, comm=p_comm_work)
          icount = p_pat%recv_count(np)*ndim2tot_sp
          IF (icount>0) CALL p_recv(recv_buf_sp(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

        ENDDO
      ELSE IF (iorder_sendrecv == 3) THEN ! use isend/irecv
        DO np = 1, p_pat%np_send ! loop over PEs where to send the data

          pid    = p_pat%pelist_send(np) ! ID of sender PE
          iss    = p_pat%send_startidx(np)
          icount = p_pat%send_count(np)*ndim2tot_dp
          IF (icount>0) CALL p_isend(send_buf_dp(1,iss), pid, 1, p_count=icount, comm=p_comm_work)
          icount = p_pat%send_count(np)*ndim2tot_sp
          IF (icount>0) CALL p_isend(send_buf_sp(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

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
        CALL work_mpi_barrier()
        stop_sync_timer(timer_barrier)
      ENDIF

      ! Fill in receive buffer

#if defined( __SX__ ) || defined( _OPENACC )
      DO n = 1, nfields_dp
        recv_fld_dp => recv_dp(n)%fld
!$ACC PARALLEL &
!$ACC PRESENT( p_pat, recv_buf_dp, recv_fld_dp ), COPYIN( kshift_dp, noffset_dp, ndim2_dp ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2_dp(n)
!$ACC LOOP VECTOR
          DO i = 1, p_pat%n_pnts
            recv_dp(n)%fld(p_pat%recv_dst_idx(i),k+kshift_dp(n),p_pat%recv_dst_blk(i)) =  &
              recv_buf_dp(k+noffset_dp(n),p_pat%recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
      DO n = 1, nfields_sp
        recv_fld_sp => recv_sp(n)%fld
!$ACC PARALLEL &
!$ACC PRESENT( p_pat, recv_buf_sp, recv_fld_sp ), COPYIN( kshift_sp, noffset_sp, ndim2_sp ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
!CDIR UNROLL=6
        DO k = 1, ndim2_sp(n)
!$ACC LOOP VECTOR
          DO i = 1, p_pat%n_pnts
            recv_sp(n)%fld(p_pat%recv_dst_idx(i),k+kshift_sp(n),p_pat%recv_dst_blk(i)) =  &
              recv_buf_sp(k+noffset_sp(n),p_pat%recv_src(i))
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL DO PRIVATE(jb,jl,ik)
#endif
      DO i = 1, p_pat%n_pnts
        jb = p_pat%recv_dst_blk(i)
        jl = p_pat%recv_dst_idx(i)
        ik  = p_pat%recv_src(i)
        DO n = 1, nfields_dp
          DO k = 1, ndim2_dp(n)
            recv_dp(n)%fld(jl,k+kshift_dp(n),jb) = recv_buf_dp(k+noffset_dp(n),ik)
          ENDDO
        ENDDO
        DO n = 1, nfields_sp
          DO k = 1, ndim2_sp(n)
            recv_sp(n)%fld(jl,k+kshift_sp(n),jb) = recv_buf_sp(k+noffset_sp(n),ik)
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

    TYPE(t_comm_pattern), INTENT(IN) :: p_pat

    REAL(dp), INTENT(INOUT)           :: recv(:,:,:,:)
    REAL(dp), INTENT(IN   ), OPTIONAL :: send(:,:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_4de1"
    INTEGER, INTENT(IN)           :: nfields, ndim2tot

    INTEGER :: ndim2, koffset

    REAL(dp) :: send_buf(ndim2tot,p_pat%n_send),recv_buf(ndim2tot,p_pat%n_recv)

    INTEGER :: i, k, ik, jb, jl, n, np, irs, iss, pid, icount
    LOGICAL :: lsend

    !-----------------------------------------------------------------------

    IF(my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Not yet implemented!")
    END IF

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL work_mpi_barrier()
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    IF (PRESENT(send)) THEN
      lsend  = .TRUE.
    ELSE
      lsend  = .FALSE.
    ENDIF

    ndim2 = SIZE(recv,3)

!$ACC DATA CREATE( send_buf, recv_buf ), IF ( i_am_accel_node .AND. acc_on )

    IF ((iorder_sendrecv == 1 .OR. iorder_sendrecv == 3)) THEN
      ! Set up irecv's for receive buffers
      DO np = 1, p_pat%np_recv ! loop over PEs from where to receive the data

        pid    = p_pat%pelist_recv(np) ! ID of receiver PE
        irs    = p_pat%recv_startidx(np)
        icount = p_pat%recv_count(np)*ndim2tot
        CALL p_irecv(recv_buf(1,irs), pid, 1, p_count=icount, comm=p_comm_work)

      ENDDO
    ENDIF

#if defined( __OMPPAR_COPY__ ) && !defined( _OPENACC )
!$OMP PARALLEL DO PRIVATE(jb,jl,koffset)
#else
!$ACC PARALLEL &
!$ACC PRESENT( p_pat, send, recv, send_buf ), &
!$ACC IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
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
#if defined( __OMPPAR_COPY__ ) && !defined( _OPENACC )
!$OMP END PARALLEL DO
#else
!$ACC END PARALLEL
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
      CALL work_mpi_barrier()
      stop_sync_timer(timer_barrier)
    ENDIF

    ! Fill in receive buffer
#if defined( __OMPPAR_COPY__ ) && !defined( _OPENACC )
!$OMP PARALLEL DO PRIVATE(jb,jl,ik,koffset)
#else
!$ACC PARALLEL &
!$ACC PRESENT( p_pat, recv_buf, recv ), &
!$ACC IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
#endif
    DO i = 1, p_pat%n_pnts
      jb = p_pat%recv_dst_blk(i)
      jl = p_pat%recv_dst_idx(i)
      ik  = p_pat%recv_src(i)
!$ACC LOOP VECTOR
      DO k = 1, ndim2
        koffset = (k-1)*nfields
        DO n = 1, nfields
          recv(n,jl,k,jb) = recv_buf(n+koffset,ik)
        ENDDO
      ENDDO
    ENDDO
#if defined( __OMPPAR_COPY__ ) && !defined( _OPENACC )
!$OMP END PARALLEL DO
#else
!$ACC END PARALLEL
#endif

!$ACC END DATA

    stop_sync_timer(timer_exch_data)

  END SUBROUTINE exchange_data_4de1


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

    REAL(dp), INTENT(INOUT), TARGET, OPTIONAL ::  &
      recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv4d(:,:,:,:)

    INTEGER, INTENT(IN)           :: nfields, ndim2tot

    REAL(dp) , INTENT(INOUT) :: send_buf(:,:),recv_buf(:,:)

    TYPE t_fieldptr
      REAL(dp), POINTER :: fld(:,:,:)
    END TYPE t_fieldptr

    TYPE(t_fieldptr) :: recv(nfields)
    INTEGER        :: ndim2(nfields), noffset(nfields)

    INTEGER :: i, k, jb, jl, n, np, irs, iss, pid, icount

    !-----------------------------------------------------------------------
    start_sync_timer(timer_exch_data_async)

    IF(my_process_is_mpi_seq()) RETURN

!$ACC DATA CREATE( send_buf, recv_buf, recv ), IF ( i_am_accel_node .AND. acc_on )

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
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL &
!$ACC PRESENT( p_pat, recv, send_buf ), COPYIN( noffset, ndim2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
    DO n = 1, nfields
!CDIR UNROLL=6
      DO k = 1, ndim2(n)
        DO i = 1, p_pat%n_send
          send_buf(k+noffset(n),i) = &
            recv(n)%fld(p_pat%send_src_idx(i),k,p_pat%send_src_blk(i))
        ENDDO
      ENDDO
    ENDDO
!$ACC END PARALLEL
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

#ifndef __USE_G2G
!$ACC UPDATE HOST( send_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    ! Finally issue the non-blocking send calls
    DO np = 1, p_pat%np_send ! loop over PEs where to send the data

      pid    = p_pat%pelist_send(np) ! ID of sender PE
      iss    = p_pat%send_startidx(np)
      icount = p_pat%send_count(np)*ndim2tot
      CALL p_isend(send_buf(1,iss), pid, 1, p_count=icount, comm=p_comm_work)

    ENDDO

!$ACC END DATA

    stop_sync_timer(timer_exch_data_async)

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

    REAL(dp), INTENT(INOUT), TARGET, OPTIONAL ::  &
      recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv4d(:,:,:,:)

    INTEGER, INTENT(IN)      :: nfields

    REAL(dp) , INTENT(INOUT) :: recv_buf(:,:)

    TYPE t_fieldptr
      REAL(dp), POINTER :: fld(:,:,:)
    END TYPE t_fieldptr

    TYPE(t_fieldptr) :: recv(nfields)
    INTEGER        :: ndim2(nfields), noffset(nfields)

    INTEGER :: i, k, ik, jb, jl, n

    !-----------------------------------------------------------------------
    start_sync_timer(timer_exch_data_async)

    IF(my_process_is_mpi_seq()) RETURN

!$ACC DATA CREATE( recv ), IF ( i_am_accel_node .AND. acc_on )

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
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( recv_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif

    ! Fill in receive buffer

#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL &
!$ACC PRESENT( p_pat, recv_buf, recv ), COPYIN( noffset, ndim2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
    DO n = 1, nfields
!$ACC LOOP VECTOR
!CDIR UNROLL=6
      DO k = 1, ndim2(n)
        DO i = 1, p_pat%n_pnts
          recv(n)%fld(p_pat%recv_dst_idx(i),k,p_pat%recv_dst_blk(i)) =   &
            recv_buf(k+noffset(n),p_pat%recv_src(i))
        ENDDO
      ENDDO
    ENDDO
!$ACC END PARALLEL
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

!$ACC END DATA

    stop_sync_timer(timer_exch_data_async)

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
  SUBROUTINE exchange_data_grf(p_pat, nfields, ndim2tot, recv1, send1, &
    recv2, send2, recv3, send3, recv4, send4, &
    recv5, send5, recv6, send6, recv4d1, send4d1, &
    recv4d2, send4d2)

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat(:)

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
    TYPE(t_fieldptr_send) :: send(nfields*SIZE(p_pat))

    INTEGER        :: ndim2(nfields), noffset(nfields),            &
      ioffset_s(SIZE(p_pat)), ioffset_r(SIZE(p_pat))

    REAL(dp), ALLOCATABLE :: send_buf(:,:),recv_buf(:,:), &
      auxs_buf(:,:),auxr_buf(:,:)

    INTEGER :: i, k, ik, jb, jl, n, np, irs, ire, iss, ise, &
      npats, isum, ioffset, isum1, n4d, pid, num_send, num_recv, j
    INTEGER, ALLOCATABLE :: pelist_send(:), pelist_recv(:)

    !-----------------------------------------------------------------------

    nsendtot = SUM(p_pat(:)%n_send)
    nrecvtot = SUM(p_pat(:)%n_recv)

    ALLOCATE(send_buf(ndim2tot,nsendtot),recv_buf(ndim2tot,nrecvtot), &
      auxs_buf(ndim2tot,nsendtot),auxr_buf(ndim2tot,nrecvtot))

    !-----------------------------------------------------------------------

    IF (itype_exch_barrier == 1 .OR. itype_exch_barrier == 3) THEN
      start_sync_timer(timer_barrier)
      CALL work_mpi_barrier()
      stop_sync_timer(timer_barrier)
    ENDIF

    start_sync_timer(timer_exch_data)

    npats = SIZE(p_pat)  ! Number of communication patterns provided on input

    !-----------------------------------------------------------------------

    ! some adjustmens to the standard communication patterns in order to make
    ! them work in this routine

    num_send = 0
    num_recv = 0

    DO np = 0, p_n_work-1 ! loop over PEs

      DO n = 1, npats  ! loop over communication patterns
        iss = p_pat(n)%send_limits(np)+1
        ise = p_pat(n)%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, npats  ! loop over communication patterns
        irs = p_pat(n)%recv_limits(np)+1
        ire = p_pat(n)%recv_limits(np+1)
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
    DO np = 0, p_n_work-1 ! loop over PEs

      DO n = 1, npats  ! loop over communication patterns
        iss = p_pat(n)%send_limits(np)+1
        ise = p_pat(n)%send_limits(np+1)
        IF(ise >= iss) THEN
          num_send = num_send + 1
          DO j = 1, npats
            pelist_send(num_send) = np
          ENDDO
          EXIT ! count each processor only once
        ENDIF
      ENDDO

      DO n = 1, npats  ! loop over communication patterns
        irs = p_pat(n)%recv_limits(np)+1
        ire = p_pat(n)%recv_limits(np+1)
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

!$ACC DATA CREATE( send_buf, recv_buf, auxs_buf, auxr_buf), IF ( i_am_accel_node .AND. acc_on )

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

!$ACC PARALLEL &
!$ACC PRESENT( p_pat, send, recv ), COPYIN( ndim2 ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO np = 1, npats
!$ACC LOOP VECTOR
        DO i = 1, p_pat(np)%n_pnts
          DO n = 1, nfields
            DO k = 1, ndim2(n)
              recv(n)%fld( p_pat(np)%recv_dst_idx(i), k, p_pat(np)%recv_dst_blk(i) ) =  &
                send(np+(n-1)*npats)%fld( k, &
                idx_1d(p_pat(np)%send_src_idx(p_pat(np)%recv_src(i)), &
                p_pat(np)%send_src_blk(p_pat(np)%recv_src(i))))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL

      stop_sync_timer(timer_exch_data)

    ELSE    ! WS: removed RETURN statement to allow for OpenACC DATA region


      ! Set up send buffer
#if defined( __SX__ ) || defined( _OPENACC )
!$ACC PARALLEL &
!$ACC PRESENT( p_pat, send, send_buf ), COPYIN( ndim2, noffset ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO n = 1, nfields
!$ACC LOOP VECTOR
        DO np = 1, npats
!CDIR UNROLL=6
          DO k = 1, ndim2(n)
            DO i = 1, p_pat(np)%n_send
              send_buf(k+noffset(n),i+ioffset_s(np)) =                &
                & send(np+(n-1)*npats)%fld(k,idx_1d(p_pat(np)%send_src_idx(i), &
                &                                   p_pat(np)%send_src_blk(i)))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL
#endif
      DO np = 1, npats
#ifdef __OMPPAR_COPY__
!$OMP DO PRIVATE(jl)
#endif
        DO i = 1, p_pat(np)%n_send
          jl = idx_1d(p_pat(np)%send_src_idx(i), p_pat(np)%send_src_blk(i))
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
            iss = p_pat(n)%send_limits(pid)+1 + ioffset_s(n)
            ise = p_pat(n)%send_limits(pid+1) + ioffset_s(n)
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
            p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)

          ioffset = isum

        ENDDO
      ELSE IF (iorder_sendrecv == 2) THEN ! use isend/recv
        ioffset = 0
        DO np = 1, num_send ! loop over PEs where to send the data

          pid = pelist_send(np) ! ID of sender PE

          ! Copy send points for all communication patterns into one common send buffer
          isum = ioffset
          DO n = 1, npats
            iss = p_pat(n)%send_limits(pid)+1 + ioffset_s(n)
            ise = p_pat(n)%send_limits(pid+1) + ioffset_s(n)
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
            p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)

          ioffset = isum

        ENDDO

        ioffset = 0
        DO np = 1, num_recv ! loop over PEs from where to receive the data

          pid = pelist_recv(np) ! ID of receiver PE

          ! Sum up receive points over all communication patterns to be processed
          isum = ioffset
          DO n = 1, npats
            isum = isum + p_pat(n)%recv_limits(pid+1) - p_pat(n)%recv_limits(pid)
          ENDDO

          IF(isum > ioffset) CALL p_recv(auxr_buf(1,ioffset+1), pid, 1,             &
            p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)
#ifndef __USE_G2G
!$ACC UPDATE DEVICE( auxs_buf ), IF ( i_am_accel_node .AND. acc_on )
#endif
          ioffset = isum

        ENDDO
      ELSE IF (iorder_sendrecv >= 3) THEN ! use isend/recv
#ifdef __OMPPAR_COPY__
!$OMP PARALLEL PRIVATE(ioffset,pid,isum,iss,ise,isum1)
#endif
        ioffset = 0
        DO np = 1, num_send ! loop over PEs where to send the data

          pid = pelist_send(np) ! ID of sender PE

          ! Copy send points for all communication patterns into one common send buffer
          isum = ioffset
          DO n = 1, npats
            iss = p_pat(n)%send_limits(pid)+1 + ioffset_s(n)
            ise = p_pat(n)%send_limits(pid+1) + ioffset_s(n)
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
            p_count=(isum-ioffset)*ndim2tot, comm=p_comm_work)
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
        CALL work_mpi_barrier()
        stop_sync_timer(timer_barrier)
      ENDIF

      ! Copy exchanged data back to receive buffer

#ifdef __OMPPAR_COPY__
!$OMP PARALLEL PRIVATE(ioffset,pid,isum,irs,ire,isum1)
#endif

      ioffset = 0
      DO np = 1, num_recv ! loop over PEs from where to receive the data

        pid = pelist_recv(np) ! ID of receiver PE

        isum = ioffset
        DO n = 1, npats
          irs = p_pat(n)%recv_limits(pid)+1 + ioffset_r(n)
          ire = p_pat(n)%recv_limits(pid+1) + ioffset_r(n)
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
!$ACC PARALLEL PRESENT( p_pat, recv_buf, recv ), COPYIN( ndim2, noffset, ioffset_r ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO n = 1, nfields
        DO np = 1, npats
!$ACC LOOP VECTOR
!CDIR UNROLL=6
          DO k = 1, ndim2(n)
            DO i = 1, p_pat(np)%n_pnts
              recv(n)%fld(p_pat(np)%recv_dst_idx(i),k,p_pat(np)%recv_dst_blk(i)) =   &
                recv_buf(k+noffset(n),p_pat(np)%recv_src(i)+ioffset_r(np))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL
#else
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
  SUBROUTINE exchange_data_r2d(p_pat, recv, send, add, send_lbound2, l_recv_exists)
    !
    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    INTEGER, OPTIONAL :: send_lbound2
    LOGICAL, OPTIONAL :: l_recv_exists

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_r2d"
    REAL(dp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_r2d_seq(p_pat, recv, send, add, send_lbound2,l_recv_exists)
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
  SUBROUTINE exchange_data_s2d(p_pat, recv, send, add, send_lbound2, l_recv_exists)
    !
    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    INTEGER, OPTIONAL :: send_lbound2
    LOGICAL, OPTIONAL :: l_recv_exists

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_s2d"
    REAL(sp) :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_s2d_seq(p_pat, recv, send, add, send_lbound2,l_recv_exists)
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

!$ACC KERNELS PRESENT( recv, tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
    recv(:,:) = tmp_recv(:,1,:)
!$ACC END KERNELS

!$ACC END DATA

  END SUBROUTINE exchange_data_s2d


  ! SEQUENTIAL version of subroutine "exchange_data_r3d"
  !
  SUBROUTINE exchange_data_r2d_seq(p_pat, recv, send, add, send_lbound2, l_recv_exists)

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    INTEGER, OPTIONAL :: send_lbound2
    LOGICAL, OPTIONAL :: l_recv_exists
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_r2d_seq"
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
!$ACC PARALLEL PRESENT( p_pat, send, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  add( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )                +  &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound2)
      END DO
!$ACC END PARALLEL
    ELSE
!$ACC PARALLEL PRESENT( p_pat, send, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound2)
      END DO
!$ACC END PARALLEL
    END IF

  END SUBROUTINE exchange_data_r2d_seq


  ! SEQUENTIAL version of subroutine "exchange_data_r3d"
  !
  SUBROUTINE exchange_data_s2d_seq(p_pat, recv, send, add, send_lbound2, l_recv_exists)

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    INTEGER, OPTIONAL :: send_lbound2
    LOGICAL, OPTIONAL :: l_recv_exists
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_s2d_seq"
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
!$ACC PARALLEL PRESENT( p_pat, send, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  add( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )                +  &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound2)
      END DO
!$ACC END PARALLEL
    ELSE
!$ACC PARALLEL PRESENT( p_pat, send, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound2)
      END DO
!$ACC END PARALLEL
    END IF

  END SUBROUTINE exchange_data_s2d_seq


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

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_i2d"
    INTEGER :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

    !-----------------------------------------------------------------------

    ! special treatment for trivial communication patterns of
    ! sequential runs
    IF(my_process_is_mpi_seq()) THEN
      CALL exchange_data_i2d_seq(p_pat, recv, send, add, send_lbound2)
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

!$ACC KERNELS PRESENT( recv, tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
    recv(:,:) = tmp_recv(:,1,:)
!$ACC END KERNELS

!$ACC END DATA

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
    CHARACTER(*), PARAMETER :: routine = modname//":exchange_data_i2d_seq"
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
!$ACC PARALLEL PRESENT( p_pat, send, add, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  add( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )                +  &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound2)
      END DO
!$ACC END PARALLEL
    ELSE
!$ACC PARALLEL PRESENT( p_pat, send, recv ), IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG
      DO i=1,p_pat%n_pnts
        recv( p_pat%recv_dst_idx(i), p_pat%recv_dst_blk(i) )  =                    &
          &  send(p_pat%send_src_idx(p_pat%recv_src(i)),                                    &
          &       p_pat%send_src_blk(p_pat%recv_src(i))-lbound2)
      END DO
!$ACC END PARALLEL
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

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_l2d"
    LOGICAL :: tmp_recv(SIZE(recv,1),1,SIZE(recv,2))

    !-----------------------------------------------------------------------

    IF(my_process_is_mpi_seq()) THEN
      CALL finish(routine, "Not yet implemented!")
    END IF

!$ACC DATA CREATE( tmp_recv ), IF ( i_am_accel_node .AND. acc_on )

    IF (PRESENT(send) .AND. .NOT. PRESENT(l_recv_exists)) THEN
!$ACC KERNELS PRESENT( tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
      tmp_recv(:,1,:) = .FALSE.
!$ACC END KERNELS
    ELSE
!$ACC KERNELS PRESENT( recv, tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
      tmp_recv(:,1,:) = recv(:,:)
!$ACC END KERNELS
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

!$ACC KERNELS PRESENT( recv, tmp_recv ), IF ( i_am_accel_node .AND. acc_on )
    recv(:,:) = tmp_recv(:,1,:)
!$ACC END KERNELS

!$ACC END DATA

  END SUBROUTINE exchange_data_l2d


  !================================================================================================
  !
  ! Variant of exchange routine that expects input as a 1D unblocked vector and provides
  ! output as a blocked (nproma,nblks) 2D array.
  !
  SUBROUTINE exchange_data_r1d_2d(p_pat, recv, send)

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
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

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
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

    TYPE(t_comm_pattern), INTENT(IN), TARGET :: p_pat
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


  !-------------------------------------------------------------------------


  FUNCTION get_np_recv(comm_pat)
    TYPE (t_comm_pattern), INTENT(IN) :: comm_pat
    INTEGER :: get_np_recv

    get_np_recv = comm_pat%np_recv
  END FUNCTION get_np_recv


  !-------------------------------------------------------------------------


  FUNCTION get_np_send(comm_pat)
    TYPE (t_comm_pattern), INTENT(IN) :: comm_pat
    INTEGER :: get_np_send

    get_np_send = comm_pat%np_send
  END FUNCTION get_np_send


  !-------------------------------------------------------------------------


  SUBROUTINE get_pelist_recv(comm_pat, pelist_recv)
    TYPE (t_comm_pattern), INTENT(IN) :: comm_pat
    INTEGER, INTENT(OUT) :: pelist_recv(:)

    pelist_recv = comm_pat%pelist_recv
  END SUBROUTINE get_pelist_recv


  !-------------------------------------------------------------------------


  SUBROUTINE gather_r_1d_deblock(in_array, out_array, fill_value, gather_pattern)
    ! dimension (nproma, nblk)
    REAL(dp), INTENT(IN) :: in_array(:,:)
    ! dimension (global length); only required on root
    REAL(dp), INTENT(INOUT) :: out_array(:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

    REAL(dp), ALLOCATABLE :: send_buffer(:,:)
    REAL(dp), POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, idx, blk

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    num_send_points = SUM(gather_pattern%collector_send_size(:))
    ALLOCATE(send_buffer(1, num_send_points))

    DO i = 1, SIZE(gather_pattern%loc_index(:))
      idx = idx_no(gather_pattern%loc_index(i))
      blk = blk_no(gather_pattern%loc_index(i))
      send_buffer(1,i) = in_array(idx, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_r=send_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_r=collector_buffer)
    CALL out_array_to_2d(out_array, SIZE(out_array))

    DEALLOCATE(send_buffer)

  CONTAINS

    SUBROUTINE out_array_to_2d(out_array_2d, n)
      INTEGER, INTENT(IN) :: n
      REAL(dp), INTENT(INOUT) :: out_array_2d(1,n)

      CALL two_phase_gather_second(recv_buffer_r=out_array_2d, &
        fill_value=fill_value, &
        gather_pattern=gather_pattern, &
        collector_buffer_r=collector_buffer)
    END SUBROUTINE out_array_to_2d

  END SUBROUTINE gather_r_1d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE gather_s_1d_deblock(in_array, out_array, fill_value, gather_pattern)
    ! dimension (nproma, nblk)
    REAL(sp), INTENT(IN) :: in_array(:,:)
    ! dimension (global length); only required on root
    REAL(sp), INTENT(INOUT) :: out_array(:)
    REAL(sp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

    REAL(sp), ALLOCATABLE :: send_buffer(:,:)
    REAL(sp), POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, idx, blk

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    num_send_points = SUM(gather_pattern%collector_send_size(:))
    ALLOCATE(send_buffer(1, num_send_points))

    DO i = 1, SIZE(gather_pattern%loc_index(:))
      idx = idx_no(gather_pattern%loc_index(i))
      blk = blk_no(gather_pattern%loc_index(i))
      send_buffer(1,i) = in_array(idx, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_r=send_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_r=collector_buffer)
    CALL out_array_to_2d(out_array, SIZE(out_array))

    DEALLOCATE(send_buffer)

  CONTAINS

    SUBROUTINE out_array_to_2d(out_array_2d, n)
      INTEGER, INTENT(IN) :: n
      REAL(sp), INTENT(INOUT) :: out_array_2d(1,n)

      CALL two_phase_gather_second(recv_buffer_r=out_array_2d, &
        fill_value=fill_value, &
        gather_pattern=gather_pattern, &
        collector_buffer_r=collector_buffer)
    END SUBROUTINE out_array_to_2d

  END SUBROUTINE gather_s_1d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE gather_i_1d_deblock(in_array, out_array, fill_value, gather_pattern)
    ! dimension (nproma, nblk)
    INTEGER, INTENT(IN) :: in_array(:,:)
    ! dimension (global length); only required on root
    INTEGER, INTENT(INOUT) :: out_array(:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

    INTEGER, ALLOCATABLE :: send_buffer(:,:)
    INTEGER, POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, idx, blk

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    num_send_points = SUM(gather_pattern%collector_send_size(:))
    ALLOCATE(send_buffer(1, num_send_points))

    DO i = 1, SIZE(gather_pattern%loc_index(:))
      idx = idx_no(gather_pattern%loc_index(i))
      blk = blk_no(gather_pattern%loc_index(i))
      send_buffer(1,i) = in_array(idx, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_i=send_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_i=collector_buffer)
    CALL out_array_to_2d(out_array, SIZE(out_array))

    DEALLOCATE(send_buffer)

  CONTAINS

    SUBROUTINE out_array_to_2d(out_array_2d, n)
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(INOUT) :: out_array_2d(1,n)

      CALL two_phase_gather_second(recv_buffer_i=out_array_2d, &
        fill_value=fill_value, &
        gather_pattern=gather_pattern, &
        collector_buffer_i=collector_buffer)
    END SUBROUTINE out_array_to_2d

  END SUBROUTINE gather_i_1d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE gather_r_2d_deblock(in_array, out_array, fill_value, gather_pattern)
    ! dimension (nproma, nlev, nblk)
    REAL(dp), INTENT(IN) :: in_array(:,:,:)
    ! dimension (global length, nlev); only required on root
    REAL(dp), INTENT(INOUT) :: out_array(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

    REAL(dp), ALLOCATABLE :: send_buffer(:,:), recv_buffer(:,:)
    REAL(dp), POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, nlev, idx, blk

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

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
      ALLOCATE(recv_buffer(nlev, MERGE(gather_pattern%global_size, &
        &                              SUM(gather_pattern%collector_size(:)), &
        &                              PRESENT(fill_value))))
    ELSE
      ALLOCATE(recv_buffer(0,0))
    END IF

    DO i = 1, SIZE(gather_pattern%loc_index(:))
      idx = idx_no(gather_pattern%loc_index(i))
      blk = blk_no(gather_pattern%loc_index(i))
      send_buffer(:,i) = in_array(idx, :, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_r=send_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_r=collector_buffer)
    CALL two_phase_gather_second(recv_buffer_r=recv_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_r=collector_buffer)

    IF (p_pe_work == process_mpi_root_id) &
      out_array(1:SIZE(recv_buffer, 2),1:SIZE(recv_buffer, 1)) = &
      TRANSPOSE(recv_buffer(:,:))

    DEALLOCATE(send_buffer,recv_buffer)

  END SUBROUTINE gather_r_2d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE gather_i_2d_deblock(in_array, out_array, fill_value, gather_pattern)
    ! dimension (nproma, nlev, nblk)
    INTEGER, INTENT(IN) :: in_array(:,:,:)
    ! dimension (global length, nlev); only required on root
    INTEGER, INTENT(INOUT) :: out_array(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

    INTEGER, ALLOCATABLE :: send_buffer(:,:), recv_buffer(:,:)
    INTEGER, POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, nlev, idx, blk

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    nlev = SIZE(in_array, 2)

    IF (SIZE(in_array, 1) /= nproma) &
      CALL finish("gather_i_2d_deblock", &
      &         "size of first dimension of in_array is not nproma")

    IF (nlev /= SIZE(out_array, 2) .AND. p_pe_work == process_mpi_root_id) &
      CALL finish("gather_i_2d_deblock", &
      &         "second size of in_array and out_array are not the same")

    num_send_points = SUM(gather_pattern%collector_send_size(:))
    IF (SIZE(in_array, 1) * SIZE(in_array, 3) < num_send_points) &
      CALL finish("gather_i_2d_deblock", "in_array is too small")

    ALLOCATE(send_buffer(nlev, num_send_points))
    IF (p_pe_work == process_mpi_root_id) THEN
      ALLOCATE(recv_buffer(nlev, MERGE(gather_pattern%global_size, &
        &                              SUM(gather_pattern%collector_size(:)), &
        &                              PRESENT(fill_value))))
    ELSE
      ALLOCATE(recv_buffer(0,0))
    END IF

    DO i = 1, SIZE(gather_pattern%loc_index(:))
      idx = idx_no(gather_pattern%loc_index(i))
      blk = blk_no(gather_pattern%loc_index(i))
      send_buffer(:,i) = in_array(idx, :, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_i=send_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_i=collector_buffer)
    CALL two_phase_gather_second(recv_buffer_i=recv_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_i=collector_buffer)

    IF (p_pe_work == process_mpi_root_id) &
      out_array(1:SIZE(recv_buffer, 2),1:SIZE(recv_buffer, 1)) = &
      TRANSPOSE(recv_buffer(:,:))

    DEALLOCATE(send_buffer,recv_buffer)

  END SUBROUTINE gather_i_2d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE allgather_r_1d_deblock(in_array, out_array, fill_value, &
    &                               allgather_pattern)
    ! dimension (nproma, nblk)
    REAL(dp), INTENT(IN) :: in_array(:,:)
    ! dimension (global length); only required on root
    REAL(dp), INTENT(INOUT) :: out_array(:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_allgather_pattern), INTENT(IN) :: allgather_pattern

    REAL(dp), ALLOCATABLE :: send_buffer(:,:)
    REAL(dp), POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, idx, blk, n_procs, comm
    INTEGER, ALLOCATABLE :: collector_buffer_sizes(:)

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    num_send_points = SUM(allgather_pattern%gather_pattern%collector_send_size(:))
    ALLOCATE(send_buffer(1, num_send_points))

    DO i = 1, SIZE(allgather_pattern%gather_pattern%loc_index(:))
      idx = idx_no(allgather_pattern%gather_pattern%loc_index(i))
      blk = blk_no(allgather_pattern%gather_pattern%loc_index(i))
      send_buffer(1,i) = in_array(idx, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_r=send_buffer, fill_value=fill_value,&
      gather_pattern=allgather_pattern%gather_pattern, &
      collector_buffer_r=collector_buffer)
    DEALLOCATE(send_buffer)
    IF (allgather_pattern%intercomm /= MPI_COMM_NULL) THEN
      n_procs = p_comm_remote_size(allgather_pattern%intercomm)
      comm = allgather_pattern%intercomm
    ELSE
      n_procs = p_n_work
      comm = p_comm_work
    END IF
    ALLOCATE(collector_buffer_sizes(n_procs))
    CALL p_allgather(SIZE(collector_buffer, 2), collector_buffer_sizes, comm)
    IF (SIZE(out_array, 1) < SUM(collector_buffer_sizes)) &
      CALL finish("allgather_r_1d_deblock", "invalid out_array size")
    CALL p_allgatherv(collector_buffer(1,:), out_array, collector_buffer_sizes,&
      comm)

    DEALLOCATE(collector_buffer_sizes)
  END SUBROUTINE allgather_r_1d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE allgather_i_1d_deblock(in_array, out_array, fill_value, &
    &                               allgather_pattern)
    ! dimension (nproma, nblk)
    INTEGER, INTENT(IN) :: in_array(:,:)
    ! dimension (global length); only required on root
    INTEGER, INTENT(INOUT) :: out_array(:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_allgather_pattern), INTENT(IN) :: allgather_pattern

    INTEGER, ALLOCATABLE :: send_buffer(:,:)
    INTEGER, POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, idx, blk, n_procs, comm
    INTEGER, ALLOCATABLE :: collector_buffer_sizes(:)

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    num_send_points = SUM(allgather_pattern%gather_pattern%collector_send_size(:))
    ALLOCATE(send_buffer(1, num_send_points))

    DO i = 1, SIZE(allgather_pattern%gather_pattern%loc_index(:))
      idx = idx_no(allgather_pattern%gather_pattern%loc_index(i))
      blk = blk_no(allgather_pattern%gather_pattern%loc_index(i))
      send_buffer(1,i) = in_array(idx, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_i=send_buffer, fill_value=fill_value,&
      gather_pattern=allgather_pattern%gather_pattern, &
      collector_buffer_i=collector_buffer)
    DEALLOCATE(send_buffer)
    IF (allgather_pattern%intercomm /= MPI_COMM_NULL) THEN
      n_procs = p_comm_remote_size(allgather_pattern%intercomm)
      comm = allgather_pattern%intercomm
    ELSE
      n_procs = p_n_work
      comm = p_comm_work
    END IF
    ALLOCATE(collector_buffer_sizes(n_procs))
    CALL p_allgather(SIZE(collector_buffer, 2), collector_buffer_sizes, comm)
    IF (SIZE(out_array, 1) < SUM(collector_buffer_sizes)) &
      CALL finish("allgather_i_1d_deblock", "invalid out_array size")
    CALL p_allgatherv(collector_buffer(1,:), out_array, collector_buffer_sizes,&
      comm)

    DEALLOCATE(collector_buffer_sizes)
  END SUBROUTINE allgather_i_1d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE two_phase_gather_first_r(send_buffer_r, fill_value,                &
    &                                 gather_pattern, collector_buffer_r)
    ! dimension (:, length), prepared according to gather pattern
    REAL(dp), INTENT(IN) :: send_buffer_r(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    REAL(dp), POINTER, INTENT(OUT) :: collector_buffer_r(:,:)

    REAL(dp), POINTER :: collector_buffer_nofill_r(:,:)
    REAL(dp), POINTER :: collector_buffer_fill_r(:,:)
    INTEGER :: num_send_per_process(p_n_work), send_displ(p_n_work)
    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
    INTEGER :: i, collector_idx, collector_buffer_nofill_size, &
      &        collector_buffer_fill_size, num_collectors, num_points_per_coll
    LOGICAL :: use_fill_value

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    IF (ANY(gather_pattern%collector_pes == p_pe_work)) THEN
      collector_idx = MINLOC(ABS(gather_pattern%collector_pes - p_pe_work), 1)
    ELSE
      collector_idx = -1
    END IF

    IF (collector_idx /= -1) THEN
      IF (PRESENT(fill_value)) THEN
        num_collectors = SIZE(gather_pattern%collector_pes)
        num_points_per_coll = (gather_pattern%global_size + num_collectors - 1) &
          &                   / num_collectors
        collector_buffer_fill_size = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (collector_idx- 1)))
        use_fill_value = (collector_buffer_fill_size /= &
          &               gather_pattern%collector_size(collector_idx)) .AND. &
          &              (collector_buffer_fill_size > 0)
      ELSE
        collector_buffer_fill_size = gather_pattern%collector_size(collector_idx)
        use_fill_value = .FALSE.
      END IF
      collector_buffer_nofill_size = &
        gather_pattern%collector_size(collector_idx)
    ELSE
      use_fill_value = .FALSE.
      collector_buffer_nofill_size = 0
      collector_buffer_fill_size = 0
    END IF

    ALLOCATE(collector_buffer_nofill_r(SIZE(send_buffer_r, 1), &
      &      collector_buffer_nofill_size))

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

    CALL p_alltoallv(send_buffer_r, num_send_per_process, send_displ, &
      &              collector_buffer_nofill_r, num_recv_per_process, &
      &              recv_displ, p_comm_work)

    ! reorder collector_buffer
    collector_buffer_nofill_r(:,:) = &
      collector_buffer_nofill_r(:, gather_pattern%recv_buffer_reorder)

    IF (use_fill_value) THEN
      ALLOCATE(collector_buffer_fill_r(SIZE(send_buffer_r, 1), &
        &      collector_buffer_fill_size))
      collector_buffer_fill_r = fill_value
      collector_buffer_fill_r(:,gather_pattern%recv_buffer_reorder_fill) = &
        collector_buffer_nofill_r
      DEALLOCATE(collector_buffer_nofill_r)
      collector_buffer_r => collector_buffer_fill_r
    ELSE
      collector_buffer_r => collector_buffer_nofill_r
    END IF
  END SUBROUTINE two_phase_gather_first_r


  !-------------------------------------------------------------------------


  SUBROUTINE two_phase_gather_first_s(send_buffer_r, fill_value,                &
    &                                 gather_pattern, collector_buffer_r)
    ! dimension (:, length), prepared according to gather pattern
    REAL(sp), INTENT(IN) :: send_buffer_r(:,:)
    REAL(sp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    REAL(sp), POINTER, INTENT(OUT) :: collector_buffer_r(:,:)

    REAL(sp), POINTER :: collector_buffer_nofill_r(:,:)
    REAL(sp), POINTER :: collector_buffer_fill_r(:,:)
    INTEGER :: num_send_per_process(p_n_work), send_displ(p_n_work)
    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
    INTEGER :: i, collector_idx, collector_buffer_nofill_size, &
      &        collector_buffer_fill_size, num_collectors, num_points_per_coll
    LOGICAL :: use_fill_value

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    IF (ANY(gather_pattern%collector_pes == p_pe_work)) THEN
      collector_idx = MINLOC(ABS(gather_pattern%collector_pes - p_pe_work), 1)
    ELSE
      collector_idx = -1
    END IF

    IF (collector_idx /= -1) THEN
      IF (PRESENT(fill_value)) THEN
        num_collectors = SIZE(gather_pattern%collector_pes)
        num_points_per_coll = (gather_pattern%global_size + num_collectors - 1) &
          &                   / num_collectors
        collector_buffer_fill_size = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (collector_idx- 1)))
        use_fill_value = (collector_buffer_fill_size /= &
          &               gather_pattern%collector_size(collector_idx)) .AND. &
          &              (collector_buffer_fill_size > 0)
      ELSE
        collector_buffer_fill_size = gather_pattern%collector_size(collector_idx)
        use_fill_value = .FALSE.
      END IF
      collector_buffer_nofill_size = &
        gather_pattern%collector_size(collector_idx)
    ELSE
      use_fill_value = .FALSE.
      collector_buffer_nofill_size = 0
      collector_buffer_fill_size = 0
    END IF

    ALLOCATE(collector_buffer_nofill_r(SIZE(send_buffer_r, 1), &
      &      collector_buffer_nofill_size))

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

    CALL p_alltoallv(send_buffer_r, num_send_per_process, send_displ, &
      &              collector_buffer_nofill_r, num_recv_per_process, &
      &              recv_displ, p_comm_work)

    ! reorder collector_buffer
    collector_buffer_nofill_r(:,:) = &
      collector_buffer_nofill_r(:, gather_pattern%recv_buffer_reorder)
    IF (use_fill_value) THEN
      ALLOCATE(collector_buffer_fill_r(SIZE(send_buffer_r, 1), &
        &      collector_buffer_fill_size))
      collector_buffer_fill_r = fill_value
      collector_buffer_fill_r(:,gather_pattern%recv_buffer_reorder_fill) = &
        collector_buffer_nofill_r
      DEALLOCATE(collector_buffer_nofill_r)
      collector_buffer_r => collector_buffer_fill_r
    ELSE
      collector_buffer_r => collector_buffer_nofill_r
    END IF
  END SUBROUTINE two_phase_gather_first_s


  !-------------------------------------------------------------------------


  SUBROUTINE two_phase_gather_first_i(send_buffer_i, fill_value,                &
    &                                 gather_pattern, collector_buffer_i)
    ! dimension (:, length), prepared according to gather pattern
    INTEGER, INTENT(IN) :: send_buffer_i(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    INTEGER, POINTER, INTENT(OUT) :: collector_buffer_i(:,:)

    INTEGER, POINTER :: collector_buffer_nofill_i(:,:)
    INTEGER, POINTER :: collector_buffer_fill_i(:,:)
    INTEGER :: num_send_per_process(p_n_work), send_displ(p_n_work)
    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
    INTEGER :: i, collector_idx, collector_buffer_nofill_size, &
      &        collector_buffer_fill_size, num_collectors, num_points_per_coll
    LOGICAL :: use_fill_value

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    IF (ANY(gather_pattern%collector_pes == p_pe_work)) THEN
      collector_idx = MINLOC(ABS(gather_pattern%collector_pes - p_pe_work), 1)
    ELSE
      collector_idx = -1
    END IF

    IF (collector_idx /= -1) THEN
      IF (PRESENT(fill_value)) THEN
        num_collectors = SIZE(gather_pattern%collector_pes)
        num_points_per_coll = (gather_pattern%global_size + num_collectors - 1) &
          &                   / num_collectors
        collector_buffer_fill_size = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (collector_idx- 1)))
        use_fill_value = (collector_buffer_fill_size /= &
          &               gather_pattern%collector_size(collector_idx)) .AND. &
          &              (collector_buffer_fill_size > 0)
      ELSE
        collector_buffer_fill_size = gather_pattern%collector_size(collector_idx)
        use_fill_value = .FALSE.
      END IF
      collector_buffer_nofill_size = &
        gather_pattern%collector_size(collector_idx)
    ELSE
      use_fill_value = .FALSE.
      collector_buffer_nofill_size = 0
      collector_buffer_fill_size = 0
    END IF

    ALLOCATE(collector_buffer_nofill_i(SIZE(send_buffer_i, 1), &
      &      collector_buffer_nofill_size))

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

    CALL p_alltoallv(send_buffer_i, num_send_per_process, send_displ, &
      &              collector_buffer_nofill_i, num_recv_per_process, &
      &              recv_displ, p_comm_work)

    ! reorder collector_buffer
    collector_buffer_nofill_i(:,:) = &
      collector_buffer_nofill_i(:, gather_pattern%recv_buffer_reorder)
    IF (use_fill_value) THEN
      ALLOCATE(collector_buffer_fill_i(SIZE(send_buffer_i, 1), &
        &      collector_buffer_fill_size))
      collector_buffer_fill_i = fill_value
      collector_buffer_fill_i(:,gather_pattern%recv_buffer_reorder_fill) = &
        collector_buffer_nofill_i
      DEALLOCATE(collector_buffer_nofill_i)
      collector_buffer_i => collector_buffer_fill_i
    ELSE
      collector_buffer_i => collector_buffer_nofill_i
    END IF
  END SUBROUTINE two_phase_gather_first_i


  !-------------------------------------------------------------------------


  SUBROUTINE two_phase_gather_second_r(recv_buffer_r, fill_value,                &
    &                                  gather_pattern, collector_buffer_r)
    ! dimension (:, global length); only required on root
    REAL(dp), INTENT(INOUT) :: recv_buffer_r(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    REAL(dp), POINTER, INTENT(INOUT) :: collector_buffer_r(:,:)

    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
    INTEGER :: i, collector_idx, collector_buffer_fill_size, num_collectors, &
      &        num_points_per_coll
    LOGICAL :: use_fill_value

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    IF (ANY(gather_pattern%collector_pes == p_pe_work)) THEN
      collector_idx = MINLOC(ABS(gather_pattern%collector_pes - p_pe_work), 1)
    ELSE
      collector_idx = -1
    END IF

    IF (collector_idx /= -1) THEN
      IF (PRESENT(fill_value)) THEN
        num_collectors = SIZE(gather_pattern%collector_pes)
        num_points_per_coll = (gather_pattern%global_size + num_collectors - 1) &
          &                   / num_collectors
        collector_buffer_fill_size = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (collector_idx- 1)))
        use_fill_value = (collector_buffer_fill_size /= &
          &               gather_pattern%collector_size(collector_idx)) .AND. &
          &              (collector_buffer_fill_size > 0)
      ELSE
        use_fill_value = .FALSE.
      END IF
    ELSE
      use_fill_value = .FALSE.
    END IF

    num_recv_per_process(:) = 0
    IF (p_pe_work == process_mpi_root_id) THEN
      IF (use_fill_value) THEN
        num_recv_per_process(gather_pattern%collector_pes(:)+1) = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (/(i, i = 0, num_collectors - 1)/)))
      ELSE
        num_recv_per_process(gather_pattern%collector_pes(:)+1) = &
          gather_pattern%collector_size(:)
      END IF
    END IF
    recv_displ(1) = 0
    DO i = 2, p_n_work
      recv_displ(i) = recv_displ(i-1) + num_recv_per_process(i-1)
    END DO

    CALL p_gatherv(collector_buffer_r, SIZE(collector_buffer_r, 2), &
      &            recv_buffer_r, num_recv_per_process, recv_displ, &
      &            process_mpi_root_id, p_comm_work)
    DEALLOCATE(collector_buffer_r)
  END SUBROUTINE two_phase_gather_second_r


  !-------------------------------------------------------------------------


  SUBROUTINE two_phase_gather_second_s(recv_buffer_r, fill_value,                &
    &                                  gather_pattern, collector_buffer_r)
    ! dimension (:, global length); only required on root
    REAL(sp), INTENT(INOUT) :: recv_buffer_r(:,:)
    REAL(sp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    REAL(sp), POINTER, INTENT(INOUT) :: collector_buffer_r(:,:)

    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
    INTEGER :: i, collector_idx, collector_buffer_fill_size, num_collectors, &
      &        num_points_per_coll
    LOGICAL :: use_fill_value

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    IF (ANY(gather_pattern%collector_pes == p_pe_work)) THEN
      collector_idx = MINLOC(ABS(gather_pattern%collector_pes - p_pe_work), 1)
    ELSE
      collector_idx = -1
    END IF

    IF (collector_idx /= -1) THEN
      IF (PRESENT(fill_value)) THEN
        num_collectors = SIZE(gather_pattern%collector_pes)
        num_points_per_coll = (gather_pattern%global_size + num_collectors - 1) &
          &                   / num_collectors
        collector_buffer_fill_size = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (collector_idx- 1)))
        use_fill_value = (collector_buffer_fill_size /= &
          &               gather_pattern%collector_size(collector_idx)) .AND. &
          &              (collector_buffer_fill_size > 0)
      ELSE
        use_fill_value = .FALSE.
      END IF
    ELSE
      use_fill_value = .FALSE.
    END IF

    num_recv_per_process(:) = 0
    IF (p_pe_work == process_mpi_root_id) THEN
      IF (use_fill_value) THEN
        num_recv_per_process(gather_pattern%collector_pes(:)+1) = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (/(i, i = 0, num_collectors - 1)/)))
      ELSE
        num_recv_per_process(gather_pattern%collector_pes(:)+1) = &
          gather_pattern%collector_size(:)
      END IF
    END IF
    recv_displ(1) = 0
    DO i = 2, p_n_work
      recv_displ(i) = recv_displ(i-1) + num_recv_per_process(i-1)
    END DO

    CALL p_gatherv(collector_buffer_r, SIZE(collector_buffer_r, 2), &
      &            recv_buffer_r, num_recv_per_process, recv_displ, &
      &            process_mpi_root_id, p_comm_work)
    DEALLOCATE(collector_buffer_r)
  END SUBROUTINE two_phase_gather_second_s


  !-------------------------------------------------------------------------


  SUBROUTINE two_phase_gather_second_i(recv_buffer_i, fill_value,       &
    &                                  gather_pattern,                  &
    &                                  collector_buffer_i)
    ! dimension (:, global length); only required on root
    INTEGER, INTENT(INOUT) :: recv_buffer_i(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    INTEGER, POINTER, INTENT(INOUT) :: collector_buffer_i(:,:)

    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
    INTEGER :: i, collector_idx, collector_buffer_fill_size, num_collectors, &
      &        num_points_per_coll
    LOGICAL :: use_fill_value

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    IF (ANY(gather_pattern%collector_pes == p_pe_work)) THEN
      collector_idx = MINLOC(ABS(gather_pattern%collector_pes - p_pe_work), 1)
    ELSE
      collector_idx = -1
    END IF

    IF (collector_idx /= -1) THEN
      IF (PRESENT(fill_value)) THEN
        num_collectors = SIZE(gather_pattern%collector_pes)
        num_points_per_coll = (gather_pattern%global_size + num_collectors - 1) &
          &                   / num_collectors
        collector_buffer_fill_size = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (collector_idx- 1)))
        use_fill_value = (collector_buffer_fill_size /= &
          &               gather_pattern%collector_size(collector_idx)) .AND. &
          &              (collector_buffer_fill_size > 0)
      ELSE
        use_fill_value = .FALSE.
      END IF
    ELSE
      use_fill_value = .FALSE.
    END IF

    num_recv_per_process(:) = 0
    IF (p_pe_work == process_mpi_root_id) THEN
      IF (use_fill_value) THEN
        num_recv_per_process(gather_pattern%collector_pes(:)+1) = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (/(i, i = 0, num_collectors - 1)/)))
      ELSE
        num_recv_per_process(gather_pattern%collector_pes(:)+1) = &
          gather_pattern%collector_size(:)
      END IF
    END IF
    recv_displ(1) = 0
    DO i = 2, p_n_work
      recv_displ(i) = recv_displ(i-1) + num_recv_per_process(i-1)
    END DO

    CALL p_gatherv(collector_buffer_i, SIZE(collector_buffer_i, 2), &
      &            recv_buffer_i, num_recv_per_process, recv_displ, &
      &            process_mpi_root_id, p_comm_work)
    DEALLOCATE(collector_buffer_i)
  END SUBROUTINE two_phase_gather_second_i


  !-----------------------------------------------------------------------------
  !> Factory method for t_ScatterPattern. Destroy with deleteScatterPattern().
  !-----------------------------------------------------------------------------
  FUNCTION makeScatterPattern(jg, loc_arr_len, glb_index, communicator)
    USE mo_scatter_pattern_scatter
    IMPLICIT NONE
    CLASS(t_ScatterPattern), POINTER :: makeScatterPattern
    INTEGER, VALUE :: jg, loc_arr_len, communicator
    INTEGER, INTENT(IN) :: glb_index(:)

    CHARACTER(*), PARAMETER :: routine = modname//":makeScatterPattern"
    INTEGER :: ierr

    ALLOCATE(t_ScatterPatternScatter::makeScatterPattern, stat = ierr)
    IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
    CALL makeScatterPattern%construct(jg, loc_arr_len, glb_index, communicator)
  END FUNCTION makeScatterPattern

END MODULE mo_communication


