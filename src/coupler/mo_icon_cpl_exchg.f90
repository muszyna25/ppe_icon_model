!>
!! Sending and receiving ICON fields for coupling
!!
!! <Description>
!!
!! Data are send to remote processes using the non-blocking psmile_bsend.
!! Thus, the ICON_cpl_put routine returns. The user is allowed to reuse
!! the send buffer from the ICON_cpl_put while the psmile_bsend is keeping
!! trac of the memory. The ICON_cpl_get is blocking, it only returns when
!! data heve been received.
!! 
!! @author Rene Redler, Max-Planck Institute for Meteorology, Germany
!!
!! $Id:$
!!
!! @par Revision History
!! first implementation by Rene Redler (2010-02-13)
!!
!! @par Copyright
!! 2010 by MPI-M
!! This software is provided for non-commerncial use only.
!! See the LICENSE and WARRANTY conditions.
!! 
!! @par License
!!
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
MODULE mo_icon_cpl_exchg

  USE mo_kind, ONLY          : wp
  USE mo_event_manager, ONLY : event_check, events

#ifndef NOMPI
  USE mpi, ONLY : MPI_INTEGER, MPI_TAG, MPI_SOURCE, MPI_STATUS_SIZE

  USE mo_icon_cpl, ONLY : t_cpl_field, cpl_fields,  &
   &                      nbr_ICON_fields,  &
   &                      t_source_struct, source_locs, &
   &                      t_target_struct, target_locs, &
   &                      l_debug, debug_level, cplout, &
   &                      ICON_global_rank,             &
   &                      ICON_comm, ICON_comm_active,  &
   &                      initag,           &
   &                      msg_len,          &
   &                      datatype,         &
   &                      datatype,         &
   &                      cpl_field_none,   &
   &                      cpl_field_avg
  
  USE mo_icon_cpl_write_restart, ONLY : icon_cpl_write_restart

  IMPLICIT NONE

  PRIVATE

  TYPE(t_cpl_field), POINTER         :: fptr
  TYPE(t_source_struct), POINTER :: sptr
  TYPE(t_target_struct), POINTER :: tptr

  INTEGER                        :: rstatus(MPI_STATUS_SIZE) ! MPI_Irecv status
  INTEGER                        :: wstatus(MPI_STATUS_SIZE) ! MPI_Wait status

  INTEGER, ALLOCATABLE           :: lrequests(:)
  INTEGER                        :: msgtag
  INTEGER                        :: index

  INTEGER                        :: msg_to_tgt(msg_len)
  INTEGER, ALLOCATABLE           :: msg_fm_src(:,:)
  INTEGER                        :: n_send
  INTEGER                        :: n_recv
  INTEGER                        :: source_rank

  INTEGER                        :: i, n, m
  INTEGER                        :: nbr_bundles

  REAL (wp), ALLOCATABLE         :: recv_buffer(:,:)
  REAL (wp), ALLOCATABLE         :: send_buffer(:,:)

  ! Return code for error handling

  INTEGER                        :: len
  INTEGER                        :: ierr

  ! Logical to activate MPI_Send/MPI_Recv behind the put/get

  LOGICAL                        :: l_action

#else

  USE mo_icon_cpl, ONLY : nbr_ICON_fields, t_cpl_field, cpl_fields, target_locs

  IMPLICIT NONE

  PRIVATE

#endif

  PUBLIC :: ICON_cpl_get, ICON_cpl_put

CONTAINS

  SUBROUTINE ICON_cpl_get ( field_id, field_shape, recv_field, info, ierror )

    INTEGER, INTENT(in)    :: field_id         !<  field id
    INTEGER, INTENT(in)    :: field_shape(3)   !<  shape of recv field
    REAL(wp), INTENT(out)  :: recv_field (field_shape(1):field_shape(2),field_shape(3))
    INTEGER, INTENT(out)   :: info             !<  action performed
    INTEGER, INTENT(out)   :: ierror           !<  returned error code

    ierror = 0

    info   = 0

    ! -------------------------------------------------------------------
    ! Check field id and return if field was not declared.
    ! -------------------------------------------------------------------

    IF ( field_id < 1 .OR.  field_id > nbr_ICON_fields ) RETURN

    IF ( .NOT. cpl_fields(field_id)%l_field_status ) RETURN

#ifndef NOMPI

    ! -------------------------------------------------------------------
    ! Check event
    ! -------------------------------------------------------------------

    fptr => cpl_fields(field_id)

    l_action = event_check ( fptr%event_id )

    IF ( l_debug .AND. debug_level > 0 ) THEN
       WRITE ( cplout , * ) ICON_global_rank, ' : get action for event ', &
                              fptr%event_id,                              &
                              events(fptr%event_id)%time_step,            &
                              events(fptr%event_id)%elapsed_time,         &
                              events(fptr%event_id)%delta_time, l_action
    ENDIF

    IF ( .NOT. l_action ) RETURN

    ! ----------------------------------------------------------------------
    ! First check whether this process has to receive data from someone else
    ! ----------------------------------------------------------------------

    n_recv = 0

    IF ( ASSOCIATED (target_locs) ) THEN
       DO  i = 1, SIZE(target_locs)
          IF ( target_locs(i)%source_list_len > 0 ) n_recv = n_recv + 1
       ENDDO
    ELSE
       RETURN
    ENDIF

    IF ( n_recv == 0 ) RETURN

    ALLOCATE ( lrequests(n_recv), STAT = ierr )
    IF ( ierr > 0 ) THEN
       WRITE ( * , * ) ' Error allocating lrequests for field id ', field_id
       CALL MPI_Abort ( ICON_comm, 1, ierr )
    ENDIF

    ALLOCATE ( msg_fm_src(msg_len, n_recv), STAT = ierr )
    IF ( ierr > 0 ) THEN
       WRITE ( * , * ) ' Error allocating msg_fm_src for field id ', field_id
       CALL MPI_Abort ( ICON_comm, 1, ierr )
    ENDIF

    ! -------------------------------------------------------------------
    ! Post Receives of header messages
    ! -------------------------------------------------------------------

    DO n = 1, size(target_locs) ! n_recv

       tptr => target_locs(n)

       IF ( tptr%source_list_len > 0 ) THEN

          msgtag = initag + 1000 * fptr%global_field_id

          IF ( l_debug .AND. debug_level > 0 ) &
               WRITE ( cplout , * ) ICON_global_rank, ' irecv : tag ', msgtag, &
               ' length ', msg_len, ' from ', tptr%source_rank

          CALL MPI_Irecv ( msg_fm_src(1,n), msg_len, MPI_INTEGER, &
               tptr%source_rank, msgtag, ICON_comm_active, lrequests(n), ierr )

       ENDIF

    ENDDO

    ! -------------------------------------------------------------------
    ! Loop over the number of send operation (determined during the search)
    ! -------------------------------------------------------------------

    DO n = 1, n_recv

       CALL MPI_Waitany ( n_recv, lrequests, index, wstatus, ierr )

       len         = msg_fm_src(1, index)
       nbr_bundles = msg_fm_src(2, index)

       ! ----------------------------------------------------------------
       ! Receive and store lists except empty lists
       ! ----------------------------------------------------------------

       IF ( len > 0 ) THEN

          tptr => target_locs(index)

          ALLOCATE ( recv_buffer(len,nbr_bundles), STAT = ierr )
          IF ( ierr > 0 ) THEN
             WRITE ( * , '(a,i4)') ' Error allocating receive buffer for field id ', field_id
             CALL MPI_Abort ( ICON_comm, 1, ierr )
          ENDIF

          source_rank = wstatus(MPI_SOURCE)
          msgtag      = wstatus(MPI_TAG) + 1

          IF ( tptr%source_rank /= msg_fm_src(3, index) ) THEN
             WRITE ( * , '(a,2i6)') ' Error: Messages got mixed up ', &
                   tptr%source_rank, msg_fm_src(3, index)
                   CALL MPI_Abort ( ICON_comm, 1, ierr )
          ENDIF

          CALL MPI_Recv ( recv_buffer, len*nbr_bundles, datatype, &
               source_rank, msgtag, ICON_comm_active, rstatus, ierr )

          ! -------------------------------------------------------------
          ! Scatter data into user buffer
          ! -------------------------------------------------------------

          IF ( nbr_bundles /= field_shape(3) ) THEN
             WRITE ( * , '(a,i4,a4,i4)' ) 'Number of bundles does not match!', &
                nbr_bundles, ' /= ', field_shape(3)
             CALL MPI_Abort ( ICON_comm, 1, ierr )
          ENDIF

          DO m = 1, nbr_bundles
             DO i = 1, len
                recv_field(tptr%source_list(i),m) = recv_buffer(i,m)
                IF ( l_debug .AND. debug_level > 1 ) &
                     WRITE ( cplout, '(i4,a,i4,a,i4,f13.6)' ) ICON_global_rank, ' extract from ', &
                     source_rank, ' : ', tptr%source_list(i), recv_buffer(i,m) 
             ENDDO
          ENDDO

          DEALLOCATE (recv_buffer)

       ENDIF

    ENDDO

    info = 1

    DEALLOCATE ( lrequests )
    DEALLOCATE ( msg_fm_src )
#else
    recv_field(:,:) = 0.0_wp
#endif

  END SUBROUTINE ICON_cpl_get


  SUBROUTINE ICON_cpl_put ( field_id, field_shape, send_field, ierror )

    INTEGER, INTENT(in)    :: field_id         !<  field id
    INTEGER, INTENT(in)    :: field_shape(3)   !<  shape of send field

    REAL (wp), INTENT(in)  :: send_field (field_shape(1):field_shape(2),field_shape(3))

    INTEGER, INTENT(out)   :: ierror           !<  returned error code

#ifndef NOMPI

    ! Local variables

    LOGICAL                :: l_coupling
    LOGICAL                :: l_end_of_run
    REAL (wp)              :: weight

    ierror = 0

    l_coupling   = .TRUE.
    l_end_of_run = .FALSE.

    ! -------------------------------------------------------------------
    ! Check field id and return if field was not declared.
    ! -------------------------------------------------------------------

    IF ( field_id < 1 .OR.  field_id > nbr_ICON_fields ) RETURN

    IF ( .NOT. cpl_fields(field_id)%l_field_status ) RETURN

    ! -------------------------------------------------------------------
    ! First check whether this process has to send data to someone else
    ! -------------------------------------------------------------------
    IF ( ASSOCIATED (target_locs) ) THEN
       n_send = SIZE(target_locs)
    ELSE
       RETURN
    ENDIF

    IF ( n_send == 0 ) RETURN

    fptr => cpl_fields(field_id)

    ! -------------------------------------------------------------------
    ! Store data for averaging and/or accumulation
    ! -------------------------------------------------------------------

    IF ( fptr%coupling%time_operation /= cpl_field_none ) THEN

       IF ( .NOT. ASSOCIATED(fptr%send_field_acc) ) THEN
          ALLOCATE ( fptr%send_field_acc(field_shape(1):field_shape(2), &
                                           field_shape(3)), STAT = ierr )
          IF ( ierr > 0 ) THEN
             WRITE ( * , * ) ' Error allocating storage for accumulation of field id ', field_id
             CALL MPI_Abort ( ICON_comm, 1, ierr )
          ENDIF
          fptr%send_field_acc     = 0.0_wp
          fptr%accumulation_count = 0
       ENDIF
       
       fptr%accumulation_count   = fptr%accumulation_count + 1
       fptr%send_field_acc(:,:)  = fptr%send_field_acc(:,:) + send_field(:,:)

    ENDIF

    ! -------------------------------------------------------------------
    ! Check event
    ! -------------------------------------------------------------------

    l_action = event_check ( fptr%event_id )

    IF ( l_action .AND. l_end_of_run .AND. fptr%coupling%lag > 0 ) l_coupling = .FALSE.

    IF ( l_debug .AND. debug_level > 0 ) THEN
       WRITE ( cplout , * ) ICON_global_rank, ' : put action for event ', &
                                  fptr%event_id,                      &
                                  events(fptr%event_id)%time_step,    &
                                  events(fptr%event_id)%elapsed_time, &
                                  events(fptr%event_id)%delta_time, l_action
    ENDIF

    IF ( .NOT. l_action ) RETURN

    ! -------------------------------------------------------------------
    ! Loop over the number of send operation (determined during the search)
    ! -------------------------------------------------------------------

    DO n = 1, n_send

       sptr => source_locs(n)

       len         = sptr%source_list_len
       nbr_bundles = field_shape(3)

       ! ----------------------------------------------------------------
       ! Send data except of empty lists
       ! ----------------------------------------------------------------

       IF ( len > 0 ) THEN

          IF ( l_coupling ) THEN

             ALLOCATE ( send_buffer(len,nbr_bundles), STAT = ierr )
             IF ( ierr > 0 ) THEN
                WRITE ( * , * ) ' Error allocating send buffer for field id ', field_id
                CALL MPI_Abort ( ICON_comm, 1, ierr )
             ENDIF

             ! -------------------------------------------------------------
             ! Gather data into a compact send buffer
             ! -------------------------------------------------------------

             IF (  fptr%coupling%time_operation == cpl_field_none ) THEN

                DO m = 1, nbr_bundles
                   DO i = 1, len
                      send_buffer(i,m) = send_field(sptr%source_list(i),m)
                      IF ( l_debug .AND. debug_level > 1) &
                           WRITE ( cplout, '(i4,a,i4,a,f13.6)' ) ICON_global_rank, &
                               ' extract for ', sptr%target_rank, ' : ', send_buffer(i,m) 
                   ENDDO
                ENDDO

             ELSE

                IF (  fptr%coupling%time_operation == cpl_field_avg ) THEN
                   weight = 1.0_wp / REAL(fptr%accumulation_count,wp)
                ELSE
                   weight = 1.0_wp
                ENDIF

                DO m = 1, nbr_bundles
                   DO i = 1, len
                      send_buffer(i,m) = fptr%send_field_acc(sptr%source_list(i),m) * weight
                      IF ( l_debug .AND. debug_level > 1 ) &
                           WRITE ( cplout, '(i4,a,i4,a,f13.6)' ) ICON_global_rank, &
                               ' extract for ', sptr%target_rank, ' : ', send_buffer(i,m) 
                   ENDDO
                ENDDO

             ENDIF

             ! -------------------------------------------------------------
             ! Prepare and send header messages
             ! -------------------------------------------------------------

             msgtag = initag + 1000 * fptr%global_field_id

             msg_to_tgt(1) = len
             msg_to_tgt(2) = nbr_bundles
             msg_to_tgt(3) = ICON_global_rank

             CALL psmile_bsend ( msg_to_tgt, msg_len, &
                  MPI_INTEGER, sptr%target_rank, msgtag, ICON_comm_active, ierr ) 

             ! -------------------------------------------------------------
             ! Send data 
             ! -------------------------------------------------------------

             msgtag = msgtag + 1

             CALL psmile_bsend ( send_buffer, len*nbr_bundles, &
                  datatype, sptr%target_rank, msgtag, ICON_comm_active, ierr ) 

             ! -------------------------------------------------------------
             ! Deallocate send buffer as the memory management is done inside
             ! psmile_bsend 
             ! -------------------------------------------------------------

             DEALLOCATE (send_buffer)

          ELSE  ! l_coupling

             IF ( fptr%coupling%time_operation == cpl_field_none ) THEN

                CALL ICON_cpl_write_restart ( field_id, field_shape, &
                      send_field, 1, ierror )

             ELSE

                CALL ICON_cpl_write_restart ( field_id, field_shape, &
                      fptr%send_field_acc, fptr%accumulation_count, ierror )

             ENDIF

          ENDIF ! l_coupling

       ENDIF ! len > 0

    ENDDO

    ! -------------------------------------------------------------
    ! Reset accumulation buffer and counter
    ! -------------------------------------------------------------

    IF ( fptr%coupling%time_operation /= cpl_field_none ) THEN
       fptr%send_field_acc     = 0.0_wp
       fptr%accumulation_count = 0
    ENDIF

#else

    ierror = 0

#endif

  END SUBROUTINE ICON_cpl_put

END MODULE mo_icon_cpl_exchg
