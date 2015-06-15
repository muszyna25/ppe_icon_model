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
!!
!! @par Revision History
!! first implementation by Rene Redler (2010-02-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_icon_cpl_exchg

  USE mo_kind, ONLY          : wp
  USE mo_event_manager, ONLY : event_check, events
  USE mo_run_config, ONLY    : dtime, nsteps

#ifndef NOMPI
  USE mpi, ONLY : MPI_INTEGER, MPI_TAG, MPI_SOURCE, MPI_STATUS_SIZE, &
                  MPI_MIN, MPI_MAX, MPI_SUM

  USE mo_icon_cpl, ONLY : t_cpl_field, cpl_fields,  &
   &                      nbr_ICON_fields,  &
   &                      t_source_struct, source_locs, &
   &                      t_target_struct, target_locs, &
   &                      debug_coupler_level, cplout,  &
   &                      ICON_global_rank,             &
   &                      ICON_comm, ICON_comm_active,  &
   &                      ICON_comp_comm,   &
   &                      initag,           &
   &                      msg_len,          &
   &                      datatype,         &
   &                      cpl_field_none,   &
   &                      cpl_field_avg,    &
   &                      NOTHING,          &
   &                      RESTART,          &
   &                      INITIAL,          &
   &                      XCHANGE

  USE mo_master_control,  ONLY  : is_restart_run
  USE mo_io_config, ONLY        : dt_checkpoint
  USE mo_time_config, ONLY      : time_config
  USE mo_datetime, ONLY         : iso8601

  IMPLICIT NONE

  PRIVATE

  TYPE(t_cpl_field), POINTER     :: fptr
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

  REAL (wp), PARAMETER           :: dummy = -151169.0_wp

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

  PUBLIC :: ICON_cpl_get, ICON_cpl_get_field, ICON_cpl_put

CONTAINS
  !>
  !!  This routine is used on the receiver side to get a coupling field
  !!  from a remote component. It receives messages from those remote
  !!  processes which share common grid points with the caller. Received
  !!  data are gather on the local array and provided to the caller of
  !!  ICON_cpl_get. By calling ICON_cpl_get the internal event handler
  !!  gets updated. 
  !!
  SUBROUTINE ICON_cpl_get ( field_id,          &! in
                            field_shape,       &! in
                            recv_field,        &! out
                            info,              &! out
                            ierror )            ! out

    INTEGER, INTENT(in)    :: field_id         !<  field id
    INTEGER, INTENT(in)    :: field_shape(3)   !<  shape of recv field
    REAL(wp), INTENT(out)  :: recv_field (field_shape(1):field_shape(2),field_shape(3))
    INTEGER, INTENT(out)   :: info             !<  action performed
    INTEGER, INTENT(out)   :: ierror           !<  returned error code

#ifndef NOMPI

    ! Local variables and fields

    REAL(wp)               :: recv_buf(field_shape(3))
    REAL(wp)               :: recv_min(field_shape(3))
    REAL(wp)               :: recv_max(field_shape(3))
    REAL(wp)               :: recv_avg(field_shape(3))

    INTEGER                :: nsum
    LOGICAL                :: l_end_of_run = .FALSE.

    ! -------------------------------------------------------------------
    ! Initialisation
    ! -------------------------------------------------------------------

    ierror = 0
    info   = NOTHING

    ! -------------------------------------------------------------------
    ! Check field id and return if field was not declared.
    ! -------------------------------------------------------------------

    IF ( field_id < 1 .OR.  field_id > nbr_ICON_fields ) RETURN

    IF ( .NOT. cpl_fields(field_id)%l_field_status ) RETURN

    ! -------------------------------------------------------------------
    ! Check event
    ! -------------------------------------------------------------------

    fptr => cpl_fields(field_id)

    IF ( .NOT. fptr%coupling%l_activated ) THEN
       IF ( debug_coupler_level > 0 ) THEN
          WRITE ( cplout , * ) ICON_global_rank, ' : field ', &
               TRIM(fptr%field_name), ' is not activatd for coupling!'
       ENDIF
       RETURN
    ENDIF

    l_action = event_check ( fptr%event_id )

    IF ( debug_coupler_level > 1 ) THEN
       WRITE ( cplout , * ) ICON_global_rank, ' : get action for event ', &
                            fptr%event_id,                                &
                            l_action,                           &
                            events(fptr%event_id)%time_step,    &
                            events(fptr%event_id)%delta_time,   &
                            events(fptr%event_id)%elapsed_time, &
                            events(fptr%event_id)%event_time,   &
                            events(fptr%event_id)%restart_time, &
                            events(fptr%event_id)%lag
    ENDIF

    l_end_of_run = events(fptr%event_id)%elapsed_time >  &
                   events(fptr%event_id)%restart_time

    IF ( .NOT. l_action ) RETURN

    IF ( debug_coupler_level > 1 ) THEN
       WRITE ( cplout , '(A10,A8,A10,A16,L1)' ) 'Receiving ', &
            TRIM(cpl_fields(field_id)%field_name), ' for date ', &
            iso8601(time_config%cur_datetime), l_end_of_run
    ENDIF

    ! No get action at the end of the run

    IF ( l_end_of_run ) info = RESTART

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

          IF ( debug_coupler_level > 2 ) &
               WRITE ( cplout , * ) ICON_global_rank, ' irecv : tag ', msgtag, &
               ' length ', msg_len, ' from ', tptr%source_rank

          CALL MPI_Irecv ( msg_fm_src(1,n), msg_len, MPI_INTEGER, &
               tptr%source_rank, msgtag, ICON_comm_active, lrequests(n), ierr )

       ENDIF

    ENDDO

    ! ------------------------------------------------------------------------
    ! Loop over the number of receive operation (determined during the search)
    ! ------------------------------------------------------------------------

    recv_field = dummy

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

          IF ( debug_coupler_level > 1 ) THEN
             WRITE ( cplout , '(A10,A8,A10,A16)' ) 'Receiving ', &
                  TRIM(cpl_fields(field_id)%field_name), ' for date ', &
                  iso8601(time_config%cur_datetime)
          ENDIF
          
          CALL MPI_Recv ( recv_buffer, len*nbr_bundles, datatype, &
               source_rank, msgtag, ICON_comm_active, rstatus, ierr )

          ! -------------------------------------------------------------
          ! Scatter data into user buffer. Note that only internal values
          ! get updated. The coupling routines do not work on the halos!
          ! -------------------------------------------------------------

          IF ( nbr_bundles /= field_shape(3) ) THEN
             WRITE ( * , '(a,i4,a4,i4)' ) 'Number of bundles does not match!', &
                nbr_bundles, ' /= ', field_shape(3)
             CALL MPI_Abort ( ICON_comm, 1, ierr )
          ENDIF

          DO m = 1, nbr_bundles
             DO i = 1, len
                recv_field(tptr%source_list(i),m) = recv_buffer(i,m)
             ENDDO
          ENDDO

          IF ( debug_coupler_level > 2 ) THEN
             DO m = 1, nbr_bundles
                DO i = 1, len
                   WRITE ( cplout, '(i4,a,i4,a,i4,f13.6)' ) ICON_global_rank, ' extract from ', &
                           source_rank, ' : ', tptr%source_list(i), recv_buffer(i,m) 
                ENDDO
             ENDDO
             PRINT *, 'Control ', TRIM(fptr%field_name), '(1,1)', recv_buffer(1,1)
          ENDIF
 
          IF ( fptr%coupling%l_diagnostic ) THEN

             DO m = 1, nbr_bundles
                recv_min(m) = MINVAL(recv_buffer(:,m))
                recv_max(m) = MAXVAL(recv_buffer(:,m))
                recv_avg(m) = 0.0_wp
                DO i = 1, len
                   recv_avg(m) = recv_avg(m) + recv_buffer(i,m)
                ENDDO
             ENDDO
          ENDIF

          DEALLOCATE (recv_buffer)

       ENDIF

    ENDDO

    IF ( fptr%coupling%l_diagnostic ) THEN

       CALL MPI_Allreduce ( recv_min, recv_buf, nbr_bundles, datatype, &
            MPI_MIN, ICON_comp_comm, ierror )
       recv_min(:) = recv_buf(:)

       CALL MPI_Allreduce ( recv_max, recv_buf, nbr_bundles, datatype, &
            MPI_MAX, ICON_comp_comm, ierror )
       recv_max(:) = recv_buf(:)

       CALL MPI_Allreduce ( recv_avg, recv_buf, nbr_bundles, datatype, &
            MPI_SUM, ICON_comp_comm, ierror )
       recv_avg(:) = recv_buf(:)

       CALL MPI_Allreduce ( len, nsum, 1, MPI_INTEGER, &
            MPI_SUM, ICON_comp_comm, ierror )

       recv_avg(:) = recv_avg(:) / REAL(nsum,wp)

       DO i = 1, nbr_bundles
          WRITE ( cplout, '(a32,a3,3(a6,f13.6))' ) fptr%field_name, ' : ', &
            ' Min: ', recv_min(i), &      
            ' Avg: ', recv_avg(i), &
            ' Max: ', recv_max(i)
       ENDDO

    ENDIF

    info = MAX(INITIAL,info)

    DEALLOCATE ( lrequests )
    DEALLOCATE ( msg_fm_src )
#else

    PRINT *, ' Restart requires MPI! '
    PRINT *, ' field ID    ', field_id
    PRINT *, ' field shape ', field_shape

    recv_field(:,:) = 0.0_wp

    ierror = -1
    info   = 0

#endif

  END SUBROUTINE ICON_cpl_get

  ! --------------------------------------------------------------------
  !>
  !!  This routine is used on the sender side to get accumulated fields
  !!  out of the coupling layer in order to stored them in a file, possibly
  !!  a restart file. This routine does not invoke any MPI or other communication
  !!  and is just a processor local operation.
  !!
  SUBROUTINE ICON_cpl_get_field ( field_id,     &! in
                                  field_shape,  &! in
                                  data,         &! out
                                  count,        &! out
                                  info,         &! out
                                  ierror )       ! out

    INTEGER, INTENT(in)           :: field_id         !<  field id
    INTEGER, INTENT(in)           :: field_shape(3)   !<  shape of outgoing data
    REAL(wp), INTENT(out)         :: data(field_shape(1):field_shape(2),field_shape(3))
    INTEGER, INTENT(out)          :: count            !<  number of accumulations in data
    INTEGER, INTENT(out)          :: info             !<  returned info code
    INTEGER, INTENT(out)          :: ierror           !<  returned error code

    ! -------------------------------------------------------------------
    ! Initialise variables
    ! -------------------------------------------------------------------

#ifndef NOMPI

    ierror = 0
    info   = NOTHING

    ! -------------------------------------------------------------------
    ! Check field id and return if field was not declared.
    ! -------------------------------------------------------------------

    IF ( field_id < 1 .OR. field_id > nbr_ICON_fields ) RETURN

    IF ( .NOT. cpl_fields(field_id)%l_field_status ) RETURN

    fptr => cpl_fields(field_id)

    ! -------------------------------------------------------------------
    ! Check whether we have something accumulated for this field id
    ! -------------------------------------------------------------------

    IF ( fptr%coupling%time_operation /= cpl_field_none ) RETURN

    IF ( .NOT. ASSOCIATED(fptr%send_field_acc) ) RETURN

    IF ( fptr%accumulation_count == 0 ) RETURN

    info = INITIAL

    ! -------------------------------------------------------------------
    ! Get averaged and/or accumulated data
    ! -------------------------------------------------------------------

    count = fptr%accumulation_count

    data  = fptr%send_field_acc

#else

    PRINT *, ' Restart requires MPI! '
    PRINT *, ' field ID    ', field_id
    PRINT *, ' field shape ', field_shape

    data   = 0.0_wp
    count  = 0

    ierror = -1
    info   = 0

#endif
  END SUBROUTINE ICON_cpl_get_field

  ! --------------------------------------------------------------------
  !>
  !!  This routine is used on the sender side to either accumulate/average
  !!  coupling field or send instant/accumulated/averaged fields to a
  !!  remote component. By calling ICON_cpl_put the internal event handler
  !!  gets updated.
  !!
  SUBROUTINE ICON_cpl_put ( field_id,    &! in
                            field_shape, &! in
                            send_field,  &! in
                            info,        &! out
                            ierror )      ! out

    INTEGER, INTENT(in)    :: field_id         !<  field id
    INTEGER, INTENT(in)    :: field_shape(3)   !<  shape of send field

    REAL (wp), INTENT(in)  :: send_field (field_shape(1):field_shape(2),field_shape(3))

    INTEGER, INTENT(out)   :: info             !<  returned error code
    INTEGER, INTENT(out)   :: ierror           !<  returned error code

#ifndef NOMPI

    ! Local variables

    LOGICAL                :: l_restart    = .FALSE.
    LOGICAL                :: l_checkpoint = .FALSE.
    REAL (wp)              :: weight
    !
    ! for coupling diagnostic
    !
    REAL(wp)               :: send_buf(field_shape(3))
    REAL(wp)               :: send_min(field_shape(3))
    REAL(wp)               :: send_max(field_shape(3))
    REAL(wp)               :: send_avg(field_shape(3))
    INTEGER                :: j, nsum

    ierror = 0
    info   = NOTHING

    fptr => cpl_fields(field_id)

    ! -------------------------------------------------------------------
    ! Check field id and return if field was not declared.
    ! -------------------------------------------------------------------

    IF ( field_id < 1 .OR.  field_id > nbr_ICON_fields ) RETURN

    IF ( .NOT. cpl_fields(field_id)%l_field_status ) RETURN

    IF ( .NOT. cpl_fields(field_id)%coupling%l_activated ) THEN
       IF ( debug_coupler_level > 0 ) THEN
          WRITE ( cplout , * ) ICON_global_rank, ' : field ', &
               TRIM(cpl_fields(field_id)%field_name), ' is not activatd for coupling!'
       ENDIF
       RETURN
    ENDIF

    IF ( debug_coupler_level > 1 ) THEN
       WRITE ( cplout , '(A10,A8,A10,A16)' ) 'Sending   ', &
            TRIM(cpl_fields(field_id)%field_name), ' for date ', &
            iso8601(time_config%cur_datetime)
    ENDIF

    ! -------------------------------------------------------------------
    ! First check whether this process has to send data to someone else
    ! -------------------------------------------------------------------

    IF ( ASSOCIATED (target_locs) ) THEN
       n_send = SIZE(target_locs)
    ELSE
       RETURN
    ENDIF

    IF ( n_send == 0 ) RETURN

    ! -------------------------------------------------------------------
    ! Check event
    !
    ! The differentiation between end of run and restart is necessary as
    ! the non-hydrostatic atmosphere performs one more time step than requested.
    !
    ! -------------------------------------------------------------------

    l_action = event_check ( fptr%event_id )

    IF ( debug_coupler_level > 0 ) THEN
       WRITE ( cplout , * ) 'restart check: ', &
                            nsteps*dtime,                         &
                            events(fptr%event_id)%elapsed_time,   &
                            events(fptr%event_id)%restart_time
    ENDIF

    l_restart    = events(fptr%event_id)%elapsed_time == &
                   events(fptr%event_id)%restart_time

    IF ( .NOT. l_restart ) &
       l_checkpoint = events(fptr%event_id)%elapsed_time == &
                         events(fptr%event_id)%check_time

    IF ( debug_coupler_level > 1 ) THEN
       WRITE ( cplout , * ) ICON_global_rank, ' : put action for event ', &
                              fptr%event_id, l_action,                    &
                              events(fptr%event_id)%time_step,            &
                              events(fptr%event_id)%delta_time,           &
                              events(fptr%event_id)%elapsed_time,         &
                              events(fptr%event_id)%event_time,           &
                              events(fptr%event_id)%restart_time
    ENDIF

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

       IF ( debug_coupler_level > 1 ) &
         WRITE ( cplout , * ) 'Accumulation for ', TRIM(cpl_fields(field_id)%field_name)

       fptr%accumulation_count  = fptr%accumulation_count + 1
       fptr%send_field_acc(:,:) = fptr%send_field_acc(:,:) + send_field(:,:)

    ENDIF

    ! Currently, we assume that data are coupled on the last time step.
    ! Nevertheless, we have to make them available for the next run as well.

    IF ( l_checkpoint ) THEN
       ! This is bad coding style, but ...
       events(fptr%event_id)%check_time = events(fptr%event_id)%check_time &
                                        + NINT(dt_checkpoint)
    ENDIF

    ! set restart flag

    IF ( l_restart .OR. l_checkpoint ) info = RESTART

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
                   IF ( debug_coupler_level > 2 ) &
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
                   IF ( debug_coupler_level > 2 ) &
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
          msg_to_tgt(4) = XCHANGE

          CALL psmile_bsend ( msg_to_tgt, msg_len, &
               MPI_INTEGER, sptr%target_rank, msgtag, ICON_comm_active, ierr ) 

          ! -------------------------------------------------------------
          ! Send data 
          ! -------------------------------------------------------------

          msgtag = msgtag + 1

          IF ( debug_coupler_level > 1 ) THEN
             WRITE ( cplout , '(A10,A8,A10,A16)' ) 'Sending   ', &
                  TRIM(cpl_fields(field_id)%field_name), ' for date ', &
                  iso8601(time_config%cur_datetime)
          ENDIF

          CALL psmile_bsend ( send_buffer, len*nbr_bundles, &
               datatype, sptr%target_rank, msgtag, ICON_comm_active, ierr ) 

          info = XCHANGE

          IF ( fptr%coupling%l_diagnostic ) THEN

             send_avg(:) = 0.0_wp

             DO i = 1, nbr_bundles
                send_min(i) = MINVAL(send_buffer(:,i))
                send_max(i) = MAXVAL(send_buffer(:,i))
                DO j = 1, len
                   send_avg(i) = send_avg(i) + send_buffer(j,i)
                ENDDO
             ENDDO

             CALL MPI_Allreduce ( send_min, send_buf, nbr_bundles, datatype, &
                  MPI_MIN, ICON_comp_comm, ierror )
             send_min(:) = send_buf(:)

             CALL MPI_Allreduce ( send_max, send_buf, nbr_bundles, datatype, &
                  MPI_MAX, ICON_comp_comm, ierror )
             send_max(:) = send_buf(:)

             CALL MPI_Allreduce ( send_avg, send_buf, nbr_bundles, datatype, &
                  MPI_SUM, ICON_comp_comm, ierror )
             send_avg(:) = send_buf(:)

             CALL MPI_Allreduce ( len, nsum, 1, MPI_INTEGER, &
                  MPI_SUM, ICON_comp_comm, ierror )

             send_avg(:) = send_avg(:) / REAL(nsum,wp)

             DO i = 1, nbr_bundles
                WRITE ( cplout, '(a32,a3,3(a6,f13.6))' ) fptr%field_name, ' : ', &
                     ' Min: ', send_min(i), &      
                     ' Avg: ', send_avg(i), &
                     ' Max: ', send_max(i)
             ENDDO

          ENDIF

          ! -------------------------------------------------------------
          ! Deallocate send buffer as the memory management is done inside
          ! psmile_bsend 
          ! -------------------------------------------------------------

          DEALLOCATE (send_buffer)

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

    PRINT *, ' Restart requires MPI! '
    PRINT *, ' field ID    ', field_id
    PRINT *, ' field shape ', field_shape
    PRINT *, ' send field  ', send_field(1,1)

    ierror = -1 
    info   = 0

#endif

  END SUBROUTINE ICON_cpl_put

END MODULE mo_icon_cpl_exchg
