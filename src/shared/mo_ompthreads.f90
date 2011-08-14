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
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_ompthreads

!  USE mo_atm_phy_nwp_nml,      ONLY: inwp_radiation, dt_rad, dt_radheat
  USE mo_exception,            ONLY: message, warning, finish !message_tex
  USE mo_run_config,           ONLY: msg_level
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
  PUBLIC :: ompthread_has_endrequest
  PUBLIC :: ompthread_is_alive

  INTEGER, PARAMETER :: ompthreads_both_end = -2
  
  INTEGER, PARAMETER :: ompthread_notused = -1
  INTEGER, PARAMETER :: ompthread_null = 0
  INTEGER, PARAMETER :: ompthread_busy = 1
  INTEGER, PARAMETER :: ompthread_waits = 2
  INTEGER, PARAMETER :: ompthread_in_barrier = 3
  INTEGER, PARAMETER :: ompthread_ready = 4
  INTEGER, PARAMETER :: ompthread_ends = 5
  INTEGER, PARAMETER :: ompthread_syncrequest = 6
  
  TYPE t_ompthread_strc
    INTEGER :: status
    INTEGER :: request
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
!$ IF (omp_in_parallel()) THEN
!$   CALL finish("new_ompthread","cannot be called form an omp parallel region")
!$  ENDIF
    DO thread_id=1, max_ompthreads
      IF (ompthread_of(thread_id)%status == ompthread_notused) THEN
        ompthread_of(thread_id)%status      = ompthread_ready
        ompthread_of(thread_id)%sync_status = ompthread_ready
        ompthread_of(thread_id)%request     = ompthread_null
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
      & (ompthread_of(thread%id)%status /= ompthread_ends)

  END FUNCTION ompthread_is_alive
  !-----------------------------------------
  
  !-----------------------------------------
  !>
  LOGICAL FUNCTION ompthread_has_endrequest(thread)

    TYPE (t_ompthread), INTENT(in) :: thread

    ompthread_has_endrequest = &
      & (ompthread_of(thread%id)%request == ompthread_ends)

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

    sync_ompthread_to = ompthread_ready ! default value
    
    my_thread_id = thread_from%id
    to_thread_id = thread_to%id
    my_thread => ompthread_of(my_thread_id)
    to_thread => ompthread_of(to_thread_id)
    
!     ompthread_of(to_thread_id)%comm_request = ompthread_sync
!     ompthread_of(to_thread_id)%from_thread = thread_id
    my_thread%sync_status = ompthread_syncrequest
!$OMP FLUSH(ompthread_of)

    write(0,*) my_thread_id, ' waits for sync from ', to_thread_id
    wait_cnt = wait_for_syncrequest_ompthread(to_thread_id)
    write(0,*) my_thread_id, ' recieved sync from ', to_thread_id
    !-----------------------------
    ! make sure we communicate with the right thread.

    !----------------------------- 
    ! Note: if an ompthread ends, it will never be alive again
    ! ( a new ompthread has to be created )
    IF (to_thread%status == ompthread_ends) THEN
      ! the other thread requests to end
      IF (my_thread%status == ompthread_ends) THEN
        CALL warning(method_name, 'Both threads request end')
        sync_ompthread_to = ompthread_ends
        
      ELSE
        ! if the other thread requests to end
        ! then return the request and return
        CALL message(method_name, ' Received ompthread_ends')
        sync_ompthread_to = ompthread_ends        
        
      ENDIF  
    ENDIF
    
    IF (sync_ompthread_to == ompthread_ends) THEN
      write(0,*) to_thread_id, ' requested an end to ', my_thread_id
    ENDIF
    ! -----------------------------------------
    ! at this point both threads are not busy
    ! aknowledge the sync by sendind to the opposite site thread_ready
    to_thread%sync_status = ompthread_ready
!$OMP FLUSH(ompthread_of)

    ! wait until we get the ready signal
    wait_cnt = wait_for_ready_ompthread(my_thread_id)
   
    my_thread%request = sync_ompthread_to
    IF (sync_ompthread_to == ompthread_ends) THEN
      write(0,*) to_thread_id, ' ackowledged ompthread_ends'
    ENDIF
    
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
  INTEGER FUNCTION wait_for_ready_ompthread(thread_id) result(wait_cnt)

    INTEGER, INTENT(in) :: thread_id

    wait_cnt=0
    DO WHILE (.true.)
      IF ( ompthread_of(thread_id)%sync_status == ompthread_ready) EXIT
!       write(0,*) "thread_is_busy(thread_status):", thread_status
      wait_cnt = wait_cnt + 1
    ENDDO
    RETURN

  END FUNCTION wait_for_ready_ompthread
  !-----------------------------------------


  !-----------------------------------------
  !>
  SUBROUTINE end_ompthread(thread)

    TYPE (t_ompthread), INTENT(in) :: thread

    ompthread_of(thread%id)%status = ompthread_ends
!$OMP FLUSH(ompthread_of)
    write(0,*) 'Reached end_ompthread', thread%id

  END SUBROUTINE end_ompthread
  !-----------------------------------------
  

END MODULE mo_ompthreads

