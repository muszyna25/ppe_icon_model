!>
!! This module provides basic infrastracture for
!! dynamic sync of two openmp treads.
!!
!! @author
!! Leonidas Linardakis, MPI-M, 2011-08
!!
!! @par Revision History
!! Initial release by Leonidas Linardakis, MPI-M, 2011-08
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ompthreads

#ifdef _OPENMP
  USE omp_lib,       ONLY: omp_get_thread_num, omp_get_max_threads, &
                          omp_in_parallel, omp_get_num_threads
!                          omp_get_dynamic, omp_set_dynamic
#endif

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: operator(.syncto.) ! Note: this is not a commutative operator
  PUBLIC :: t_ompthread
  PUBLIC :: init_ompthreads
  PUBLIC :: new_ompthread
  PUBLIC :: end_ompthread
  PUBLIC :: request_end_ompthread
  PUBLIC :: delete_ompthread
  PUBLIC :: ompthread_has_endrequest
  PUBLIC :: ompthread_is_alive
  PUBLIC :: ompthread_is_dead

  INTEGER, PARAMETER :: ompthreads_both_end = -2
  
  INTEGER, PARAMETER :: ompthread_notused = -1
  INTEGER, PARAMETER :: ompthread_null = 0
  INTEGER, PARAMETER :: ompthread_busy = 1
  INTEGER, PARAMETER :: ompthread_waits = 2
  INTEGER, PARAMETER :: ompthread_in_barrier = 3
  INTEGER, PARAMETER :: ompthread_ready = 4
  INTEGER, PARAMETER :: ompthread_ended = 5
  INTEGER, PARAMETER :: ompthread_endrequest = 6
  
  INTEGER, PARAMETER :: ompthread_syncrequest = 7
  INTEGER, PARAMETER :: ompthread_syncaknw = 8
  
  TYPE t_ompthread_strc
    INTEGER :: status
    INTEGER :: received_request
    INTEGER :: send_request
!     INTEGER :: comm_request
!     INTEGER :: from_thread
    INTEGER :: sync_status
  END TYPE

  TYPE t_ompthread
    INTEGER :: id
  END TYPE

  ! Note: this is not a commutative operator
  INTERFACE OPERATOR (.syncto.)
    MODULE PROCEDURE sync_ompthread_to 
  END INTERFACE


  INTEGER, PARAMETER ::  max_ompthreads = 8
  TYPE(t_ompthread_strc), TARGET :: ompthread_of(max_ompthreads)


CONTAINS

  !-----------------------------------------
  !>
  INTEGER FUNCTION init_ompthreads()

    DO init_ompthreads=1, max_ompthreads
      ompthread_of(init_ompthreads)%status = ompthread_notused
    ENDDO
    init_ompthreads = max_ompthreads
    
  END FUNCTION init_ompthreads
  !-----------------------------------------

  !-----------------------------------------
  !>
  ! Returns a new ompthread 
  ! Important note: it must be called outside an openmp region
  TYPE(t_ompthread) FUNCTION new_ompthread()

    INTEGER :: thread_id
#ifdef _OPENMP
!$ IF (omp_in_parallel()) THEN
!$   CALL finish("new_ompthread","cannot be called from an omp parallel region")
!$  ENDIF
#endif
    DO thread_id=1, max_ompthreads
      IF (ompthread_of(thread_id)%status == ompthread_notused) THEN
        ompthread_of(thread_id)%status           = ompthread_ready
        ompthread_of(thread_id)%sync_status      = ompthread_null
        ompthread_of(thread_id)%received_request = ompthread_null
        ompthread_of(thread_id)%send_request     = ompthread_null
        new_ompthread%id = thread_id
        RETURN
      ENDIF
    ENDDO
    
  END FUNCTION new_ompthread
  !-----------------------------------------

  !-----------------------------------------
  !>
  LOGICAL FUNCTION ompthread_is_alive(thread)

    TYPE (t_ompthread), INTENT(in) :: thread

    ompthread_is_alive = &
      & (ompthread_of(thread%id)%status /= ompthread_ended)

  END FUNCTION ompthread_is_alive
  !-----------------------------------------
  
  !-----------------------------------------
  !>
  LOGICAL FUNCTION ompthread_is_dead(thread)

    TYPE (t_ompthread), INTENT(in) :: thread

    ompthread_is_dead = &
      & (ompthread_of(thread%id)%status == ompthread_ended)

  END FUNCTION ompthread_is_dead
  !-----------------------------------------
  
  !-----------------------------------------
  !>
  LOGICAL FUNCTION ompthread_has_endrequest(thread)

    TYPE (t_ompthread), INTENT(in) :: thread

    ompthread_has_endrequest = &
      & (ompthread_of(thread%id)%received_request == ompthread_endrequest)

  END FUNCTION ompthread_has_endrequest
  !-----------------------------------------


  !-----------------------------------------
  !>
  ! Note: this is not a commutative operator
  INTEGER FUNCTION sync_ompthread_to(thread_from, thread_to)

    TYPE (t_ompthread), INTENT(in) :: thread_from,  thread_to
    
    TYPE(t_ompthread_strc), POINTER :: my_thread, to_thread
    INTEGER  :: my_thread_id, to_thread_id
    INTEGER wait_cnt
    CHARACTER(len=*), PARAMETER :: method_name = "sync_ompthread_to"
    
    my_thread_id = thread_from%id
    to_thread_id = thread_to%id
    my_thread => ompthread_of(my_thread_id)
    to_thread => ompthread_of(to_thread_id)
    
    IF (ompthread_of(my_thread_id)%status == ompthread_ended) THEN
      write(0,*) my_thread_id, ' is dead, but requests sync from thread ', to_thread_id
      sync_ompthread_to = ompthread_ended
      my_thread%received_request = sync_ompthread_to
      RETURN
    ENDIF
    
    my_thread%sync_status = ompthread_syncrequest
!$OMP FLUSH(ompthread_of)

    !-----------------------------
    ! if the other thread is dead then cancel sync
    IF (ompthread_of(to_thread_id)%status == ompthread_ended) THEN
      write(0,*) my_thread_id, ' requests sync from dead thread ', to_thread_id
      sync_ompthread_to = ompthread_ended
      my_thread%received_request = sync_ompthread_to
      RETURN
    ENDIF

    !-----------------------------
    ! sync the threads
    write(0,*) my_thread_id, ' waits for sync from ', to_thread_id
    wait_cnt = wait_for_syncrequest_ompthread(to_thread_id)
    write(0,*) my_thread_id, ' received sync from ', to_thread_id
    !-----------------------------
    ! make sure we communicate with the right thread.
    ! ..................

    !-----------------------------
    ! get request
    sync_ompthread_to = to_thread%send_request    
    IF (sync_ompthread_to == ompthread_endrequest) THEN
      write(0,*) my_thread_id, ' received ompthread_endrequest from ', to_thread_id
    ENDIF
    
    ! -----------------------------------------
    ! at this point both threads exchanged requests
    ! aknowledge the sync by sending to the opposite site thread_ready
    ! (otherwise we cannot be sure if the other thread got the right info)
    to_thread%sync_status = ompthread_syncaknw
!$OMP FLUSH(ompthread_of)
    ! wait until we get the ready signal
    wait_cnt = wait_for_syncaknw_ompthread(my_thread_id)
    write(0,*) my_thread_id, ' ackowledged sync from ', to_thread_id
   
    !----------------------------------------------
    ! update my request
    my_thread%received_request = sync_ompthread_to
    
    RETURN
    
  END FUNCTION sync_ompthread_to
  !-----------------------------------------

  !-----------------------------------------
  !>
  INTEGER FUNCTION wait_for_syncrequest_ompthread(thread_id) result(wait_cnt)

    INTEGER, INTENT(in) :: thread_id
    
    wait_cnt=0
    DO WHILE (.true.)
      IF ( ompthread_of(thread_id)%sync_status == ompthread_syncrequest) EXIT
!       write(0,*) "thread_is_busy(thread_status):", thread_status
      wait_cnt = wait_cnt + 1
    ENDDO    
    RETURN
    
  END FUNCTION wait_for_syncrequest_ompthread
  !-----------------------------------------

  !-----------------------------------------
  !>
  INTEGER FUNCTION wait_for_syncaknw_ompthread(thread_id) result(wait_cnt)

    INTEGER, INTENT(in) :: thread_id

    wait_cnt=0
    DO WHILE (.true.)
      IF ( ompthread_of(thread_id)%sync_status == ompthread_syncaknw) EXIT
!       write(0,*) "thread_is_busy(thread_status):", thread_status
      wait_cnt = wait_cnt + 1
    ENDDO
    RETURN

  END FUNCTION wait_for_syncaknw_ompthread
  !-----------------------------------------


  !-----------------------------------------
  !>
  ! Note: if an ompthread ends, it will never be alive again
  ! ( a new ompthread has to be created )
  ! An ended ompthread will not communicate, its dead
  SUBROUTINE end_ompthread(thread)

    TYPE (t_ompthread), INTENT(in) :: thread

    ompthread_of(thread%id)%status = ompthread_ended
    write(0,*) 'Reached end_ompthread', thread%id

  END SUBROUTINE end_ompthread
  !-----------------------------------------
  

  !-----------------------------------------
  !>
  SUBROUTINE request_end_ompthread(thread)

    TYPE (t_ompthread), INTENT(in) :: thread

    ompthread_of(thread%id)%send_request = ompthread_endrequest
    write(0,*) 'set ompthread_endrequest for', thread%id

  END SUBROUTINE request_end_ompthread
  !-----------------------------------------
  

  !-----------------------------------------
  !>
  SUBROUTINE delete_ompthread(thread)

    TYPE (t_ompthread), INTENT(in) :: thread

    ompthread_of(thread%id)%status = ompthread_notused
    write(0,*) 'Reached end_ompthread', thread%id

  END SUBROUTINE delete_ompthread
  !-----------------------------------------
  

END MODULE mo_ompthreads

