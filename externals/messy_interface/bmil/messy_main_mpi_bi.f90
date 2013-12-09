#define _DIAGONALEX_
! **************************************************************************
MODULE messy_main_mpi_bi
! **************************************************************************

#if defined(ECHAM5)
  ! ECHAM5
  USE mo_mpi,           ONLY: p_parallel_io               &  ! LOGICAL
                            , p_io, p_pe, p_nprocs        &  ! INTEGER
                            , p_bcast                     &  ! SUBROUTINES
                            , p_send, p_recv, p_sendrecv  &  ! SUBROUTINES
                            , p_all_comm                  &  ! INTEGER
                            , p_abort, p_parallel         &
                            , p_sum                       &  ! SUBROUTINE
                            , p_set_communicator          &  ! mz_ab_20100307+
                            , p_barrier                      ! mz_ab_20100518
#ifdef MPIOM_2000
  USE mo_mpi,           ONLY: init_mpi_datatypes             !mz_ap_20101126
#endif
! p_set_communicator required for exchg_boundaries (CMAT submodel)
! p_barrier required for ... ???
!                                                            ! mz_ab_20100307-

  USE mo_exception,     ONLY: finish, message                ! SUBROUTINES

! um_ak_20080702+
  ! DECOMPOSITION
  USE mo_decomposition,    ONLY: dcg => global_decomposition       &
                               , dcl => local_decomposition

  ! TRANSFORMATION
  USE mo_transpose,        ONLY: gather_gp, scatter_gp, indx, reorder &
                               , gather_sa, scatter_sa  &
                               , gather_sp, scatter_sp
  ! op_pj_20110205+
  USE mo_tr_gather,        ONLY: gather_field
  ! op_pj_20110205-

! um_ak_20080702-

  IMPLICIT NONE

  PUBLIC :: exchg_boundaries_cg        ! mz_ab_20100503 for submodel CMAT
  PUBLIC :: exchg_boundaries_tracer_cg ! mz_ab_20100509 for submodel CMAT

  CONTAINS                    ! mz_ab_20100503: exchg_boundaries_cg
#endif

! mz_ab_20100223+
#ifdef MBM_CMAT
  USE messy_main_constants_mem, ONLY: dp
  ! ECHAM5
  USE mo_mpi,           ONLY: p_parallel_io               &  ! LOGICAL
                            , p_io, p_pe, p_nprocs        &  ! INTEGER
                            , p_bcast                     &  ! SUBROUTINES
                            , p_send, p_recv, p_sendrecv  &  ! SUBROUTINES
                            , p_all_comm                  &  ! INTEGER
                            , p_abort, p_parallel         &
                            , p_sum                       &  ! SUBROUTINE
                            , p_set_communicator          &  ! SUBROUTINE
                            , p_start, p_stop             &  ! SUBROUTINES
                            , p_barrier                      ! mz_ab_20100518

  IMPLICIT NONE

  PUBLIC :: exchg_boundaries_cg ! exchange boundaries of DC_CG 
  PUBLIC :: exchg_boundaries_tracer_cg  ! exchange tracer boundaries of DC_CG 
  PUBLIC :: reorder           ! reorder arrays in gridpoint space 

  TYPE decomp
     LOGICAL :: lreg = .TRUE.
  END TYPE decomp
  TYPE(decomp), POINTER :: dcg
  TYPE(decomp), SAVE    :: dcl

  INTERFACE reorder
    MODULE PROCEDURE reorder2
    MODULE PROCEDURE reorder3
    MODULE PROCEDURE reorder4
  END INTERFACE

CONTAINS
#endif
! mz_ab_20100223-

! mz_ab_20100503+
#if defined(ECHAM5) || defined(MBM_CMAT)
SUBROUTINE exchg_boundaries_cg                                               &
               (sendbuf, isendbuflen,                                        &
                 idim, jdim, kdim, jstartpar, jendpar, nlines,               &
                 neighbors, neighbors_pole, do_pole, ntag, ntype, ierror,    &
                 var01, var02, var03, var04, var05, var06, var07, var08,     &
                 var09, var10, var11, var12, var13, var14, var15, var16,     &
                 var17, var18, var19, var20, var21, var22, var23, var24 )
!------------------------------------------------------------------------------
!
! Description:
!   This subroutine performs the boundary exchange of up to 20 variables. Only
!   one variable has to occur, the others are optional. 
!
! Method:
!   At the moment there are 3 different MPI-communications implemented:
!     1) immediate send, blocking receive and wait
!     2) immediate receive, blocking send and wait
!     3) MPI_SendRecv
!   Also there is the choice of an explicit buffering (putbuf, getbuf) or 
!   implicit buffering (MPI-Datatypes) of the data to be send.
!
!------------------------------------------------------------------------------

  USE messy_main_constants_mem, ONLY: dp

#ifndef NOMPI
  INCLUDE 'mpif.h'
#endif

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    isendbuflen,        & ! length of sendbuffer
 !   imp_type,           & ! determines the REAL type used
 !   icomm,              & ! communicator for virtual cartesian topology
    idim, jdim,         & ! horizontal dimensions of the fields
    kdim,               & ! array for the vertical dimensions of var1..var20
    jstartpar,          & ! start index in j-direction
    jendpar,            & ! end index in j-direction
    nlines,             & ! number of lines that have to be exchanged
                          ! (<= nboundlines)
 !   nboundlines,        & ! number of overlapping boundary lines
#ifdef _DIAGONALEX_
    neighbors(8),       & ! process-id's of the neighbors in the grid
    neighbors_pole(8),  & ! process-id's of the neighbors in the grid
#else
    neighbors(4),       & ! process-id's of the neighbors in the grid
    neighbors_pole(4),  & ! process-id's of the neighbors in the grid
#endif
    ntag,               & ! tag of the message
    ntype                 ! indicates how the communication should be done


  LOGICAL, INTENT(IN) :: do_pole
!!$#ifdef _DIAGONALEX_
!!$  LOGICAL, INTENT(IN) :: neighbors_pole(8)
!!$#else
!!$  LOGICAL, INTENT(IN) :: neighbors_pole(4)
!!$#endif

  INTEGER :: icomm 

  INTEGER, INTENT (OUT)        ::    &
    ierror                ! error status variable

!!$  CHARACTER (LEN=*),        INTENT(OUT)  ::       &
!!$    yerrmsg               ! for MPI error message
  CHARACTER(len=*), PARAMETER :: substr='exchg_boundaries_cg' 

  REAL (dp),       INTENT (INOUT)      ::    &
    sendbuf (:,:) !,& ! send buffer
!    var01 (kdim, jdim, idim)   ! first field that has to occur
 REAL (dp), DIMENSION(:,:,:), POINTER :: var01
  

  REAL (dp), DIMENSION(:,:,:), POINTER,  OPTIONAL  ::    &
    var02,& ! additional optional fields
    var03,& ! additional optional fields
    var04,& ! additional optional fields
    var05,& ! additional optional fields
    var06,& ! additional optional fields
    var07,& ! additional optional fields
    var08,& ! additional optional fields
    var09,& ! additional optional fields
    var10,& ! additional optional fields
    var11,& ! additional optional fields
    var12,& ! additional optional fields
    var13,& ! additional optional fields
    var14,& ! additional optional fields
    var15,& ! additional optional fields
    var16,& ! additional optional fields
    var17,& ! additional optional fields
    var18,& ! additional optional fields
    var19,& ! additional optional fields
    var20,& ! additional optional fields
    var21,& ! additional optional fields
    var22,& ! additional optional fields
    var23,& ! additional optional fields
    var24   ! additional optional fields

#ifndef NOMPI

! Local variables

  INTEGER   ::       &
    ! the following numbers are for filling/emptying the buffers for each
    ! neighbor
    izlo_lr, izup_lr, jzlo_lr, jzup_lr,     & ! left , receive
    izlo_rr, izup_rr, jzlo_rr, jzup_rr,     & ! right, receive
    izlo_ur, izup_ur, jzlo_ur, jzup_ur,     & ! upper, receive
    izlo_dr, izup_dr, jzlo_dr, jzup_dr,     & ! down , receive
    izlo_ls, izup_ls, jzlo_ls, jzup_ls,     & ! left , send
    izlo_rs, izup_rs, jzlo_rs, jzup_rs,     & ! right, send
    izlo_us, izup_us, jzlo_us, jzup_us,     & ! upper, send
    izlo_ds, izup_ds, jzlo_ds, jzup_ds,     & ! down , send
#ifdef _DIAGONALEX_
    izlo_lur, izup_lur, jzlo_lur, jzup_lur,     & ! left + upper, receive
    izlo_rdr, izup_rdr, jzlo_rdr, jzup_rdr,     & ! right + down, receive
    izlo_urr, izup_urr, jzlo_urr, jzup_urr,     & ! upper + right, receive
    izlo_dlr, izup_dlr, jzlo_dlr, jzup_dlr,     & ! down + left, receive
    izlo_lus, izup_lus, jzlo_lus, jzup_lus,     & ! left + upper, send
    izlo_rds, izup_rds, jzlo_rds, jzup_rds,     & ! right + down, send
    izlo_urs, izup_urs, jzlo_urs, jzup_urs,     & ! upper + right, send
    izlo_dls, izup_dls, jzlo_dls, jzup_dls,     & ! down + left, send
#endif
    nzcount_ls, nzcount_rs,     & ! counting the values
    nzcount_us, nzcount_ds,     & ! counting the values
    nzcount_lr, nzcount_rr,     & ! counting the values
    nzcount_ur, nzcount_dr,     & ! counting the values
    nzcount_usp, nzcount_dsp,   & ! counting the values for poles
    nzcount_urp, nzcount_drp,   & ! counting the values for poles
#ifdef _DIAGONALEX_
    nzcount_lus, nzcount_rds,     & ! counting the values
    nzcount_urs, nzcount_dls,     & ! counting the values
    nzcount_lur, nzcount_rdr,     & ! counting the values
    nzcount_urr, nzcount_dlr,     & ! counting the values
#endif
    nzrequest(MPI_STATUS_SIZE), & ! for MPI-receive
    nzstatus (MPI_STATUS_SIZE), & ! for MPI-WAIT
    ncount, type_handle,        & ! return values from setup_data_type
#ifdef _DIAGONALEX_
    MPI_neighbors(8), i,        & ! same as neighbors, if neighbor exists
    MPI_neighbors_pole(8),      & ! same as neighbors, if neighbor exists
#else
    MPI_neighbors(4), i,        & ! same as neighbors, if neighbor exists
    MPI_neighbors_pole(4),      & ! same as neighbors, if neighbor exists
#endif
    ilocalreq(4)                  ! the local requests for the ISEND and IRECV

  INTEGER   ::       &
    izmplcode                   ! for MPI error code

  INTEGER, PARAMETER :: imp_type = MPI_DOUBLE_PRECISION
!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------
  icomm = p_all_comm

  ierror     = 0
  izmplcode  = 0

  ! Determine the start- and end-indices for routines putbuf, getbuf
  ! e.g.   sendbuf, isendbuflen, idim, jdim, kdim,                &
  !        izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )

  ! izlo_ls
  !   ||
  !   lo/up index (lo:up)
  !      |
  !      left/right/up/down
  !       |
  !       send/receive


! possibly swap ur/dr and us/ds
! mz_ab_20100514 update: count latitudes from south to north, 
! so order here is correct!

  izlo_ls = 1
  izup_ls = nlines
  jzlo_ls = jstartpar
  jzup_ls = jendpar

  izlo_lr = 1 - nlines
  izup_lr = 0
  jzlo_lr = jstartpar
  jzup_lr = jendpar

  izlo_us = 1
  izup_us = idim
  jzlo_us = jdim - nlines + 1
  jzup_us = jdim 

  izlo_ur = 1
  izup_ur = idim
  jzlo_ur = jdim + 1
  jzup_ur = jdim + nlines

  izlo_rs = idim - nlines + 1
  izup_rs = idim
  jzlo_rs = jstartpar
  jzup_rs = jendpar

  izlo_rr = idim + 1
  izup_rr = idim + nlines
  jzlo_rr = jstartpar
  jzup_rr = jendpar

  izlo_ds = 1
  izup_ds = idim
  jzlo_ds = 1
  jzup_ds = nlines

  izlo_dr = 1
  izup_dr = idim
  jzlo_dr = 1-nlines
  jzup_dr = 0

  nzcount_lr = 0
  nzcount_rr = 0
  nzcount_ur = 0
  nzcount_dr = 0
  nzcount_ls = 0
  nzcount_rs = 0
  nzcount_us = 0
  nzcount_ds = 0

  nzcount_usp = 0
  nzcount_dsp = 0
  nzcount_urp = 0
  nzcount_drp = 0

#ifdef _DIAGONALEX_
  ! Diagonal exchange
  ! lu = left + up
  ! ld = left + down 
  ! ...
  ! (k,j,i); jstartpar=1; jendpar=lat_dim
  izlo_lus = 1
  izup_lus = nlines
  jzlo_lus = jdim - nlines + 1
  jzup_lus = jdim

  izlo_lur = 1 - nlines
  izup_lur = 0
  jzlo_lur = jdim + 1
  jzup_lur = jdim + nlines

  izlo_urs = idim - nlines + 1
  izup_urs = idim
  jzlo_urs = jdim - nlines + 1
  jzup_urs = jdim

  izlo_urr = idim + 1
  izup_urr = idim + nlines
  jzlo_urr = jdim + 1
  jzup_urr = jdim + nlines

  izlo_rds = idim - nlines + 1
  izup_rds = idim
  jzlo_rds = 1
  jzup_rds = nlines

  izlo_rdr = idim + 1
  izup_rdr = idim + nlines
  jzlo_rdr = 1-nlines
  jzup_rdr = 0

  izlo_dls = 1
  izup_dls = nlines
  jzlo_dls = 1
  jzup_dls = nlines

  izlo_dlr = 1-nlines
  izup_dlr = 0
  jzlo_dlr = 1-nlines
  jzup_dlr = 0

  nzcount_lur = 0
  nzcount_rdr = 0
  nzcount_urr = 0
  nzcount_dlr = 0
  nzcount_lus = 0
  nzcount_rds = 0
  nzcount_urs = 0
  nzcount_dls = 0
  ! (end of diagonal exchange)
#endif

  ! Fix list of neighbors (use MPI_PROC_NULL rather than -1 to indicate
  ! missing neighbor).
#ifdef _DIAGONALEX_
  DO i= 1, 8
#else
  DO i= 1, 4
#endif
    IF ( neighbors(i) /= -1 ) THEN
      MPI_neighbors(i) = neighbors(i)
    ELSE
      MPI_neighbors(i) = MPI_PROC_NULL
    ENDIF
    IF ( neighbors_pole(i) /= -1 ) THEN
      MPI_neighbors_pole(i) = neighbors_pole(i)
    ELSE
      MPI_neighbors_pole(i) = MPI_PROC_NULL
    ENDIF
  ENDDO


!------------------------------------------------------------------------------
!- Section 3: Exchange with immediate Send and blocking Recv
!------------------------------------------------------------------------------

  IF   (ntype == 1) THEN

     ! mz_ab_20100630+
     ! IN ORDER TO IMPLEMENT POLE EXCHANGE, ADD ADDITIONAL BLOCK
     ! AS IN ntype=3 BECAUSE SEND+RECEIVE USES SAME DESTINATION CPU     
     ! mz_ab_20100630-

      !------------------------------------------------------------------------
      !- Section 3.3: exchange with left and right neigh. using explict buff.
      !------------------------------------------------------------------------

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
  
        ! determine start- and end-indices for routine putbuf
        nzcount_ls = 0
        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),   &
                         ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_ISEND')
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! right neighbor is present

        ! determine start- and end-indices for routine putbuf
        nzcount_rs = 0
        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),   &
                         ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_ISEND')
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! receive the data
        CALL MPI_RECV ( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_RECV')
        ENDIF

        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
        ! receive the data

        CALL MPI_RECV ( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_RECV')
        ENDIF

        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! wait for the completion of the last send to neighbors(1)
        ! to safely reuse sendbuf(1,1)
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! wait for the completion of the last send to neighbors(3)
        ! to safely reuse sendbuf(1,3)
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            CALL p_abort !(substr,'MPI_WAIT')
          ENDIF
      ENDIF

      !------------------------------------------------------------------------
      !- Section 3.4: exchange with lower and upper neigh. using explict buff.
      !------------------------------------------------------------------------

!!$      IF ((neighbors(2) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(2))) THEN
      IF (neighbors(2) /= -1) THEN
        ! upper neighbor is present

        nzcount_us = 0
        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),   &
                         ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_ISEND')
         ENDIF
      ENDIF

!!$      IF ((neighbors(4) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(4))) THEN
      IF (neighbors(4) /= -1) THEN
        ! lower neighbor is present

        nzcount_ds = 0
        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),   &
                         ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_ISEND')
        ENDIF
      ENDIF

!!$      IF ((neighbors(4) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(4))) THEN
      IF (neighbors(4) /= -1) THEN
        ! lower neighbor is present
        ! receive the data
 
        CALL MPI_RECV ( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_RECV')
        ENDIF
  
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF
  
!!$      IF ((neighbors(2) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(2))) THEN
      IF (neighbors(2) /= -1) THEN
        ! upper neighbor is present
        ! receive the data
  
        CALL MPI_RECV ( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_RECV')
        ENDIF
  
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF
  
!!$      IF ((neighbors(2) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(2))) THEN
      IF (neighbors(2) /= -1) THEN
        ! wait for the completion of the last send to neighbors(2)
        ! to safely reuse sendbuf(1,2)
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF
      ENDIF
  
!!$      IF ((neighbors(4) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(4))) THEN
      IF (neighbors(4) /= -1) THEN
        ! wait for the completion of the last send to neighbors(4)
        ! to safely reuse sendbuf(1,4)
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF
      ENDIF

 
!------------------------------------------------------------------------------
!- Section 4: Exchange with immediate Recv and blocking Send
!------------------------------------------------------------------------------

  ELSEIF (ntype == 2) THEN


      !------------------------------------------------------------------------
      !- Section 4.3: exchange with left and right neighbors
      !------------------------------------------------------------------------

      IF (neighbors(3) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),  &
                        ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_IRECV')
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),  &
                        ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_IRECV')
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        nzcount_ls = 0
        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )

        CALL MPI_SEND ( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_SEND')
         ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        nzcount_rs = 0
        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )

        CALL MPI_SEND ( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_SEND')
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF

        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF

        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      !------------------------------------------------------------------------
      !- Section 4.4: exchange with upper and lower neighbors
      !------------------------------------------------------------------------

!!$      IF ((neighbors(2) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(2))) THEN
      IF (neighbors(2) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),  &
                        ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_IRECV')
        ENDIF
      ENDIF

!!$      IF ((neighbors(4) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(4))) THEN
      IF (neighbors(4) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),  &
                        ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_IRECV')
        ENDIF
      ENDIF

!!$      IF ((neighbors(4) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(4))) THEN
      IF (neighbors(4) /= -1) THEN

        nzcount_ds = 0
        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )

        CALL MPI_SEND ( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_SEND')
        ENDIF
      ENDIF

!!$      IF ((neighbors(2) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(2))) THEN
      IF (neighbors(2) /= -1) THEN

        nzcount_us = 0
        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )
 
        CALL MPI_SEND ( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_SEND')
        ENDIF
      ENDIF

!!$      IF ((neighbors(2) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(2))) THEN
      IF (neighbors(2) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF

        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF

 
!!$      IF ((neighbors(4) /= -1) .AND. (do_pole .OR. .NOT.neighbors_pole(4))) THEN
      IF (neighbors(4) /= -1) THEN
        ! Now wait until the data have arrived
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF

        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF


!------------------------------------------------------------------------------
!- Section 5: Exchange with SendRecv
!------------------------------------------------------------------------------

  ELSEIF (ntype == 3) THEN


      !--------------------------------------------------------------------------
      !- Section 5.5: Send data to the left and receive from the right neighbor
      !--------------------------------------------------------------------------

      nzcount_ls = 0
      IF (MPI_neighbors(1) /= MPI_PROC_NULL) THEN
        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,1), nzcount_ls,  imp_type, MPI_neighbors(1), ntag,    &
             sendbuf(1,7), isendbuflen, imp_type, MPI_neighbors(3), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(3) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.6: Send data to the right and receive from the left neighbor
      !--------------------------------------------------------------------------

      nzcount_rs = 0
      IF (MPI_neighbors(3) /= MPI_PROC_NULL) THEN
        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,3), nzcount_rs,  imp_type, MPI_neighbors(3), ntag,    &
             sendbuf(1,5), isendbuflen, imp_type, MPI_neighbors(1), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(1) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.7: Send data to the upper and receive from the lower neighbor
      !--------------------------------------------------------------------------
 
      nzcount_us = 0
      IF (MPI_neighbors(2) /= MPI_PROC_NULL) THEN

        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,2), nzcount_us,  imp_type, MPI_neighbors(2), ntag,    &
             sendbuf(1,8), isendbuflen, imp_type, MPI_neighbors(4), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(4) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.8: Send data to the lower and receive from the upper neighbor
      !--------------------------------------------------------------------------
 
      nzcount_ds = 0
      IF (MPI_neighbors(4) /= MPI_PROC_NULL) THEN

        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,4), nzcount_ds,  imp_type, MPI_neighbors(4), ntag,    &
             sendbuf(1,6), isendbuflen, imp_type, MPI_neighbors(2), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(2) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.9: Send data to neighbors across the poles
      !--------------------------------------------------------------------------
  
      ! Testing if neighbours correct: grep "11:" xmessy_cmat.48.941.log | grep -A8 Neighbor
      ! util/bnd_exchg_eval.jnl
      IF_DO_POLE : IF (do_pole) THEN

         nzcount_usp = 0
         IF (MPI_neighbors_pole(2) /= MPI_PROC_NULL) THEN

            CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                 var09, var10, var11, var12, var13, var14, var15 ,var16,&
                 var17, var18, var19, var20, var21, var22, var23 ,var24,&
                 sendbuf, isendbuflen, idim, jdim, kdim,                &
                 izlo_us, izup_us, jzlo_us, jzup_us, nzcount_usp, 9 )
         ENDIF
 
         CALL MPI_SENDRECV                                                      &
              ( sendbuf(1,9), nzcount_usp,  imp_type, MPI_neighbors_pole(2), ntag,    &
              sendbuf(1,10), isendbuflen, imp_type, MPI_neighbors_pole(2), ntag,    &
              icomm, nzstatus, izmplcode)
         IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            CALL p_abort !(substr,'MPI_SENDRECV')
         ENDIF

         nzcount_urp = 0
         IF (MPI_neighbors_pole(2) /= MPI_PROC_NULL) THEN
            CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                 var09, var10, var11, var12, var13, var14, var15 ,var16,&
                 var17, var18, var19, var20, var21, var22, var23 ,var24,&
                 sendbuf, isendbuflen, idim, jdim, kdim,                &
                 izlo_ur, izup_ur, jzup_ur, jzlo_ur, nzcount_urp, 10, .TRUE. )
         ENDIF

         nzcount_dsp = 0
         IF (MPI_neighbors_pole(4) /= MPI_PROC_NULL) THEN

            CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                 var09, var10, var11, var12, var13, var14, var15 ,var16,&
                 var17, var18, var19, var20, var21, var22, var23 ,var24,&
                 sendbuf, isendbuflen, idim, jdim, kdim,                &
                 izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_dsp, 11 )
         ENDIF

         CALL MPI_SENDRECV                                                      &
              ( sendbuf(1,11), nzcount_dsp,  imp_type, MPI_neighbors_pole(4), ntag,    &
              sendbuf(1,12), isendbuflen, imp_type, MPI_neighbors_pole(4), ntag,    &
              icomm, nzstatus, izmplcode)
         IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            CALL p_abort !(substr,'MPI_SENDRECV')
         ENDIF

         nzcount_drp = 0
         IF (MPI_neighbors_pole(4) /= MPI_PROC_NULL) THEN
            CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                 var09, var10, var11, var12, var13, var14, var15 ,var16,&
                 var17, var18, var19, var20, var21, var22, var23 ,var24,&
                 sendbuf, isendbuflen, idim, jdim, kdim,                &
                 izlo_dr, izup_dr, jzup_dr, jzlo_dr, nzcount_drp, 12, .TRUE. )
         ENDIF

      END IF IF_DO_POLE

#ifdef _DIAGONALEX_
      !-------------------------------------------------------------------------------------
      !- Section 5.10: Send data to the lower-left and receive from the upper right neighbor
      !-------------------------------------------------------------------------------------
 
      nzcount_dls = 0
      IF (MPI_neighbors(6) /= MPI_PROC_NULL) THEN

        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dls, izup_dls, jzlo_dls, jzup_dls, nzcount_dls, 13 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,13), nzcount_dls,  imp_type, MPI_neighbors(6), ntag,    &
             sendbuf(1,14),  isendbuflen, imp_type, MPI_neighbors(7), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(7) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_urr, izup_urr, jzlo_urr, jzup_urr, nzcount_urr, 14 )
      ENDIF

      !-------------------------------------------------------------------------------------
      !- Section 5.11: Send data to the lower-right and receive from the upper left neighbor
      !-------------------------------------------------------------------------------------
 
      nzcount_rds = 0
      IF (MPI_neighbors(8) /= MPI_PROC_NULL) THEN

        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rds, izup_rds, jzlo_rds, jzup_rds, nzcount_rds, 15 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,15), nzcount_rds,  imp_type, MPI_neighbors(8), ntag,    &
             sendbuf(1,16),  isendbuflen, imp_type, MPI_neighbors(5), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(5) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lur, izup_lur, jzlo_lur, jzup_lur, nzcount_lur, 16 )
      ENDIF

      !-------------------------------------------------------------------------------------
      !- Section 5.12: Send data to the upper-left and receive from the lower right neighbor
      !-------------------------------------------------------------------------------------
 
      nzcount_lus = 0
      IF (MPI_neighbors(5) /= MPI_PROC_NULL) THEN

        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lus, izup_lus, jzlo_lus, jzup_lus, nzcount_lus, 17 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,17), nzcount_lus,  imp_type, MPI_neighbors(5), ntag,    &
             sendbuf(1,18),  isendbuflen, imp_type, MPI_neighbors(8), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(8) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rdr, izup_rdr, jzlo_rdr, jzup_rdr, nzcount_rdr, 18 )
      ENDIF

      !-------------------------------------------------------------------------------------
      !- Section 5.13: Send data to the upper-right and receive from the lower left neighbor
      !-------------------------------------------------------------------------------------
 
      nzcount_urs = 0
      IF (MPI_neighbors(7) /= MPI_PROC_NULL) THEN

        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_urs, izup_urs, jzlo_urs, jzup_urs, nzcount_urs, 19 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,19), nzcount_urs,  imp_type, MPI_neighbors(7), ntag,    &
             sendbuf(1,20),  isendbuflen, imp_type, MPI_neighbors(6), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(6) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dlr, izup_dlr, jzlo_dlr, jzup_dlr, nzcount_dlr, 20 )
      ENDIF

      !-------------------------------------------------------------------------------------
      !- Section 5.14: Diagonal pole exchange
      !-------------------------------------------------------------------------------------
 
      IF_DO_POLE2 : IF (do_pole) THEN

      ! 7 - 5
      nzcount_urs = 0
      IF (MPI_neighbors_pole(7) /= MPI_PROC_NULL) THEN

        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_urs, izup_urs, jzlo_urs, jzup_urs, nzcount_urs, 21 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,21), nzcount_urs,  imp_type, MPI_neighbors_pole(7), ntag,    &
             sendbuf(1,22),  isendbuflen, imp_type, MPI_neighbors_pole(5), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      nzcount_lur = 0
      IF (MPI_neighbors_pole(5) /= MPI_PROC_NULL) THEN

        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lur, izup_lur, jzup_lur, jzlo_lur, nzcount_lur, 22, .TRUE. )
      ENDIF

      !-------------------------------------------------------------------------------------
      ! 5 - 7
      !-------------------------------------------------------------------------------------
      nzcount_lus = 0
      IF (MPI_neighbors_pole(5) /= MPI_PROC_NULL) THEN

        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lus, izup_lus, jzlo_lus, jzup_lus, nzcount_lus, 23 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,23), nzcount_lus,  imp_type, MPI_neighbors_pole(5), ntag,    &
             sendbuf(1,24),  isendbuflen, imp_type, MPI_neighbors_pole(7), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      nzcount_urr = 0
      IF (MPI_neighbors_pole(7) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_urr, izup_urr, jzup_urr, jzlo_urr, nzcount_urr, 24, .TRUE. )
      ENDIF

      !-------------------------------------------------------------------------------------
      ! 8 - 6
      !-------------------------------------------------------------------------------------
      nzcount_rds = 0
      IF (MPI_neighbors_pole(8) /= MPI_PROC_NULL) THEN

        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rds, izup_rds, jzlo_rds, jzup_rds, nzcount_rds, 25 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,25), nzcount_rds,  imp_type, MPI_neighbors_pole(8), ntag,    &
             sendbuf(1,26),  isendbuflen, imp_type, MPI_neighbors_pole(6), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      nzcount_dlr = 0
      IF (MPI_neighbors_pole(6) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dlr, izup_dlr, jzup_dlr, jzlo_dlr, nzcount_dlr, 26, .TRUE. )
      ENDIF

      !-------------------------------------------------------------------------------------
      ! 6 - 8
      !-------------------------------------------------------------------------------------
      nzcount_dls = 0
      IF (MPI_neighbors_pole(6) /= MPI_PROC_NULL) THEN

        CALL putbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dls, izup_dls, jzlo_dls, jzup_dls, nzcount_dls, 27 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,27), nzcount_dls,  imp_type, MPI_neighbors_pole(6), ntag,    &
             sendbuf(1,28),  isendbuflen, imp_type, MPI_neighbors_pole(8), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      nzcount_rdr = 0
      IF (MPI_neighbors_pole(8) /= MPI_PROC_NULL) THEN
        CALL getbuf ( var01, var02, var03, var04, var05, var06, var07, var08,&
                      var09, var10, var11, var12, var13, var14, var15 ,var16,&
                      var17, var18, var19, var20, var21, var22, var23 ,var24,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rdr, izup_rdr, jzup_rdr, jzlo_rdr, nzcount_rdr, 28, .TRUE. )
      ENDIF

      END IF IF_DO_POLE2
#endif

  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
#endif
END SUBROUTINE exchg_boundaries_cg

!==============================================================================
!==============================================================================
!+ This subroutine puts all necessary values into sendbuf
!------------------------------------------------------------------------------

SUBROUTINE putbuf  (var01, var02, var03, var04, var05, var06, var07, var08, &
                    var09, var10, var11, var12, var13, var14, var15, var16, &
                    var17, var18, var19, var20, var21, var22, var23, var24, &
                    sendbuf, isendbuflen, idim, jdim, kdim,                 &
                    ilo, iup, jlo, jup, ncount, nentry )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine puts the necessary values from the present variables
!   (determined by ilo, iup, jlo, jup) into sendbuf.
!
! Method:
!   Check which variables are present.
!
!------------------------------------------------------------------------------

  USE messy_main_constants_mem, ONLY: dp

! Subroutine arguments
  INTEGER, INTENT (IN)      ::    &
    isendbuflen,                  & ! length of sendbuffer
    idim, jdim, kdim,             & ! dimensions of the fields
    ilo, iup, jlo, jup,           & ! start- and end-indices
    nentry                          ! specifies the row of the sendbuf

  INTEGER, INTENT (INOUT)   ::    &
    ncount                          ! counts the variables

  REAL (dp), INTENT (INOUT) ::    &
    sendbuf (:,:)!,     & ! send buffer
 !   var01  (kdim, jdim, idim)   ! first field that has to occur
 REAL (dp), DIMENSION(:,:,:), POINTER :: var01

  REAL (dp), DIMENSION(:,:,:), POINTER, OPTIONAL  ::    &
    var02,& ! additional optional fields
    var03,& ! additional optional fields
    var04,& ! additional optional fields
    var05,& ! additional optional fields
    var06,& ! additional optional fields
    var07,& ! additional optional fields
    var08,& ! additional optional fields
    var09,& ! additional optional fields
    var10,& ! additional optional fields
    var11,& ! additional optional fields
    var12,& ! additional optional fields
    var13,& ! additional optional fields
    var14,& ! additional optional fields
    var15,& ! additional optional fields
    var16,& ! additional optional fields
    var17,& ! additional optional fields
    var18,& ! additional optional fields
    var19,& ! additional optional fields
    var20,& ! additional optional fields
    var21,& ! additional optional fields
    var22,& ! additional optional fields
    var23,& ! additional optional fields
    var24   ! additional optional fields

! Local variables

  INTEGER   ::       &
    i, j, k, nzc

  LOGICAL                    ::       &
    lpres02, lpres03, lpres04, lpres05, lpres06, lpres07, lpres08, lpres09, &
    lpres10, lpres11, lpres12, lpres13, lpres14, lpres15, lpres16, lpres17, &
    lpres18, lpres19, lpres20, lpres21, lpres22, lpres23, lpres24

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ! check which variables are present
  lpres02 = PRESENT (var02)
  lpres03 = PRESENT (var03)
  lpres04 = PRESENT (var04)
  lpres05 = PRESENT (var05)
  lpres06 = PRESENT (var06)
  lpres07 = PRESENT (var07)
  lpres08 = PRESENT (var08)
  lpres09 = PRESENT (var09)
  lpres10 = PRESENT (var10)
  lpres11 = PRESENT (var11)
  lpres12 = PRESENT (var12)
  lpres13 = PRESENT (var13)
  lpres14 = PRESENT (var14)
  lpres15 = PRESENT (var15)
  lpres16 = PRESENT (var16)
  lpres17 = PRESENT (var17)
  lpres18 = PRESENT (var18)
  lpres19 = PRESENT (var19)
  lpres20 = PRESENT (var20)
  lpres21 = PRESENT (var21)
  lpres22 = PRESENT (var22)
  lpres23 = PRESENT (var23)
  lpres24 = PRESENT (var24)

!------------------------------------------------------------------------------
!- Section 2: Put data into the buffer
!------------------------------------------------------------------------------

  ! use nzc as a local counter  (based on a work from Mike O'Neill to 
  ! improve vectorization of putbuf and getbuf)
  nzc = ncount

  ! first variable that has to be present
  DO k = 1, kdim
    DO j = jlo, jup
      DO i = ilo, iup
        nzc = nzc + 1
        sendbuf (nzc,nentry) = var01(k,j,i)
      ENDDO
    ENDDO
  ENDDO

  ! optional variables that are present
  IF (lpres02 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var02(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres03 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var03(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres04 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var04(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres05 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var05(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres06 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var06(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres07 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var07(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres08 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var08(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres09 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var09(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres10 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var10(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres11 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var11(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres12 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var12(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres13 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var13(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres14 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var14(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres15 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var15(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres16 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var16(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres17 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var17(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres18 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var18(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres19 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var19(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres20 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var20(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres21 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var21(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres22 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var22(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres23 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var23(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres24 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          sendbuf (nzc,nentry) = var24(k,j,i)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! put nzc to global counter
  ncount = nzc

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  
END SUBROUTINE putbuf 

!==============================================================================
!==============================================================================
!+ This subroutine gets all necessary values from sendbuf
!------------------------------------------------------------------------------

SUBROUTINE getbuf  (var01, var02, var03, var04, var05, var06, var07, var08, &
                    var09, var10, var11, var12, var13, var14, var15, var16, &
                    var17, var18, var19, var20, var21, var22, var23, var24, &
                    sendbuf, isendbuflen, idim, jdim, kdim,                 &
                    ilo, iup, jlo, jup, ncount, nentry, reverse_jorder )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine gets the necessary values for the present variables
!   (determined by ilo, iup, jlo, jup) from sendbuf.
!
! Method:
!   Check which variables are present.
!
!------------------------------------------------------------------------------

  USE messy_main_constants_mem, ONLY: dp

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    isendbuflen,                  & ! length of sendbuffer
    idim, jdim, kdim,             & ! dimensions of the fields
    ilo, iup, jlo, jup,           & ! start- and end-indices
    nentry                          ! specifies the row of sendbuf to be used

  INTEGER, INTENT (INOUT)      ::    &
    ncount                          ! counts the variables

  REAL(dp), INTENT (INOUT)            ::    &
    sendbuf (:, :)!,     & ! send buffer
    !var01  (kdim, jdim, idim)       ! first field that has to occur
  REAL (dp), DIMENSION(:,:,:), POINTER :: var01

  REAL(dp), DIMENSION(:,:,:), POINTER, OPTIONAL  ::    &
    var02,& ! additional optional fields
    var03,& ! additional optional fields
    var04,& ! additional optional fields
    var05,& ! additional optional fields
    var06,& ! additional optional fields
    var07,& ! additional optional fields
    var08,& ! additional optional fields
    var09,& ! additional optional fields
    var10,& ! additional optional fields
    var11,& ! additional optional fields
    var12,& ! additional optional fields
    var13,& ! additional optional fields
    var14,& ! additional optional fields
    var15,& ! additional optional fields
    var16,& ! additional optional fields
    var17,& ! additional optional fields
    var18,& ! additional optional fields
    var19,& ! additional optional fields
    var20,& ! additional optional fields
    var21,& ! additional optional fields
    var22,& ! additional optional fields
    var23,& ! additional optional fields
    var24   ! additional optional fields

  ! mz_ab_20100915+
  LOGICAL, OPTIONAL :: reverse_jorder
  ! mz_ab_20100915-

! Local variables

  INTEGER    ::    i, j, k, nzc
  INTEGER    ::    jorder         ! mz_ab_20100915
  LOGICAL                    ::       &
    lpres02, lpres03, lpres04, lpres05, lpres06, lpres07, lpres08, lpres09, &
    lpres10, lpres11, lpres12, lpres13, lpres14, lpres15, lpres16, lpres17, &
    lpres18, lpres19, lpres20, lpres21, lpres22, lpres23, lpres24

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ! check which variables are present
  lpres02 = PRESENT (var02)
  lpres03 = PRESENT (var03)
  lpres04 = PRESENT (var04)
  lpres05 = PRESENT (var05)
  lpres06 = PRESENT (var06)
  lpres07 = PRESENT (var07)
  lpres08 = PRESENT (var08)
  lpres09 = PRESENT (var09)
  lpres10 = PRESENT (var10)
  lpres11 = PRESENT (var11)
  lpres12 = PRESENT (var12)
  lpres13 = PRESENT (var13)
  lpres14 = PRESENT (var14)
  lpres15 = PRESENT (var15)
  lpres16 = PRESENT (var16)
  lpres17 = PRESENT (var17)
  lpres18 = PRESENT (var18)
  lpres19 = PRESENT (var19)
  lpres20 = PRESENT (var20)
  lpres21 = PRESENT (var21)
  lpres22 = PRESENT (var22)
  lpres23 = PRESENT (var23)
  lpres24 = PRESENT (var24)

  ! mz_ab_20100915+
  jorder = 1
  IF (PRESENT(reverse_jorder)) THEN
     IF (reverse_jorder) jorder = -1
  ENDIF
  ! mz_ab_20100915-

!------------------------------------------------------------------------------
!- Section 2: Get data from the buffer
!------------------------------------------------------------------------------

  ! use nzc as a local counter  (based on a work from Mike O'Neill to 
  ! improve vectorization of putbuf and getbuf)
  nzc = ncount

  ! first variable that has to be present
  DO k = 1, kdim
    DO j = jlo, jup, jorder ! mz_ab_20100915
      DO i = ilo, iup
        nzc = nzc + 1
        var01(k,j,i) = sendbuf (nzc,nentry)
      ENDDO
    ENDDO
  ENDDO

  ! optional variables that are present
  IF (lpres02.EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var02(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres03.EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var03(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres04.EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var04(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres05.EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var05(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres06.EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var06(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres07.EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var07(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres08.EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var08(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres09 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var09(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres10 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var10(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres11 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var11(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres12 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var12(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres13 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var13(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres14 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var14(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres15 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var15(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres16 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var16(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres17 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var17(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres18 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var18(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres19 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var19(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres20 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup
        DO i = ilo, iup
          nzc = nzc + 1
          var20(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres21 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var21(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres22 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var22(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres23 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var23(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! optional variables that are present
  IF (lpres24 .EQV. .TRUE.) THEN
    DO k = 1, kdim
      DO j = jlo, jup, jorder ! mz_ab_20100915
        DO i = ilo, iup
          nzc = nzc + 1
          var24(k,j,i) = sendbuf (nzc,nentry)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! put nzc to global counter
  ncount = nzc

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  
END SUBROUTINE getbuf 
  ! mz_ab_20100210-

SUBROUTINE exchg_boundaries_tracer_cg                                        &
               (sendbuf, isendbuflen,                                        &
                 idim, jdim, kdim, odim, jstartpar, jendpar, nlines,         &
                 neighbors, ntag, ntype, ierror,                             &
                 var01)
! Description: see above exchg_boundaries_cg

  USE messy_main_constants_mem, ONLY: dp

#ifndef NOMPI
  INCLUDE 'mpif.h'
#endif

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    isendbuflen,        & ! length of sendbuffer
    idim, jdim,         & ! horizontal dimensions of the fields
    kdim, odim,         & ! array for the vertical dimensions of var1..var20
    jstartpar,          & ! start index in j-direction
    jendpar,            & ! end index in j-direction
    nlines,             & ! number of lines that have to be exchanged
                          ! (<= nboundlines)
    neighbors(4),       & ! process-id's of the neighbors in the grid
    ntag,               & ! tag of the message
    ntype                 ! indicates how the communication should be done

  INTEGER :: icomm 

  INTEGER, INTENT (OUT)        ::    &
    ierror                ! error status variable

  CHARACTER(len=*), PARAMETER :: substr='exchg_boundaries_tracer_cg' 

  REAL (dp),       INTENT (INOUT)      ::    &
    sendbuf (:,:) !,& ! send buffer
  REAL (dp), DIMENSION(:,:,:,:), POINTER :: var01
  
#ifndef NOMPI

! Local variables

  INTEGER   ::       &
    ! the following numbers are for filling/emptying the buffers for each
    ! neighbor
    izlo_lr, izup_lr, jzlo_lr, jzup_lr,     & ! left , receive
    izlo_rr, izup_rr, jzlo_rr, jzup_rr,     & ! right, receive
    izlo_ur, izup_ur, jzlo_ur, jzup_ur,     & ! upper, receive
    izlo_dr, izup_dr, jzlo_dr, jzup_dr,     & ! down , receive
    izlo_ls, izup_ls, jzlo_ls, jzup_ls,     & ! left , send
    izlo_rs, izup_rs, jzlo_rs, jzup_rs,     & ! right, send
    izlo_us, izup_us, jzlo_us, jzup_us,     & ! upper, send
    izlo_ds, izup_ds, jzlo_ds, jzup_ds,     & ! down , send
    nzcount_ls, nzcount_rs,     & ! counting the values
    nzcount_us, nzcount_ds,     & ! counting the values
    nzcount_lr, nzcount_rr,     & ! counting the values
    nzcount_ur, nzcount_dr,     & ! counting the values
    nzrequest(MPI_STATUS_SIZE), & ! for MPI-receive
    nzstatus (MPI_STATUS_SIZE), & ! for MPI-WAIT
    ncount, type_handle,        & ! return values from setup_data_type
    MPI_neighbors(4), i,        & ! same as neighbors, if neighbor exists
    ilocalreq(4)                  ! the local requests for the ISEND and IRECV

  INTEGER   ::       &
    izmplcode                   ! for MPI error code

  INTEGER, PARAMETER :: imp_type = MPI_DOUBLE_PRECISION
!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------
  icomm = p_all_comm

  ierror     = 0
  izmplcode  = 0

  ! Determine the start- and end-indices for routines putbuf, getbuf
  ! e.g.   sendbuf, isendbuflen, idim, jdim, kdim,                &
  !        izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )

  ! izlo_ls
  !   ||
  !   lo/up index (lo:up)
  !      |
  !      left/right/up/down
  !       |
  !       send/receive


  izlo_ls = 1
  izup_ls = nlines
  jzlo_ls = jstartpar
  jzup_ls = jendpar

  izlo_lr = 1 - nlines
  izup_lr = 0
  jzlo_lr = jstartpar
  jzup_lr = jendpar

  izlo_us = 1
  izup_us = idim
  jzlo_us = jdim - nlines + 1
  jzup_us = jdim

  izlo_ur = 1
  izup_ur = idim
  jzlo_ur = jdim + 1
  jzup_ur = jdim + nlines

  izlo_rs = idim - nlines + 1
  izup_rs = idim
  jzlo_rs = jstartpar
  jzup_rs = jendpar

  izlo_rr = idim + 1
  izup_rr = idim + nlines
  jzlo_rr = jstartpar
  jzup_rr = jendpar

  izlo_ds = 1
  izup_ds = idim
  jzlo_ds = 1
  jzup_ds = nlines

  izlo_dr = 1
  izup_dr = idim
  jzlo_dr = 1-nlines
  jzup_dr = 0

  nzcount_lr = 0
  nzcount_rr = 0
  nzcount_ur = 0
  nzcount_dr = 0
  nzcount_ls = 0
  nzcount_rs = 0
  nzcount_us = 0
  nzcount_ds = 0

  ! Fix list of neighbors (use MPI_PROC_NULL rather than -1 to indicate
  ! missing neighbor).
  DO i= 1, 4
    IF ( neighbors(i) /= -1 ) THEN
      MPI_neighbors(i) = neighbors(i)
    ELSE
      MPI_neighbors(i) = MPI_PROC_NULL
    ENDIF
  ENDDO


!------------------------------------------------------------------------------
!- Section 3: Exchange with immediate Send and blocking Recv
!------------------------------------------------------------------------------

  IF   (ntype == 1) THEN

      !------------------------------------------------------------------------
      !- Section 3.3: exchange with left and right neigh. using explict buff.
      !------------------------------------------------------------------------

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
  
        ! determine start- and end-indices for routine putbuf
        nzcount_ls = 0
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),   &
                         ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_ISEND')
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! right neighbor is present

        ! determine start- and end-indices for routine putbuf
        nzcount_rs = 0
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),   &
                         ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_ISEND')
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! receive the data
        CALL MPI_RECV ( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_RECV')
        ENDIF

        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
        ! receive the data

        CALL MPI_RECV ( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_RECV')
        ENDIF

        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! wait for the completion of the last send to neighbors(1)
        ! to safely reuse sendbuf(1,1)
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! wait for the completion of the last send to neighbors(3)
        ! to safely reuse sendbuf(1,3)
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            CALL p_abort !(substr,'MPI_WAIT')
          ENDIF
      ENDIF

      !------------------------------------------------------------------------
      !- Section 3.4: exchange with lower and upper neigh. using explict buff.
      !------------------------------------------------------------------------

      IF (neighbors(2) /= -1) THEN
        nzcount_us = 0
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),   &
                         ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_ISEND')
         ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! lower neighbor is present

        nzcount_ds = 0
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),   &
                         ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_ISEND')
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! lower neighbor is present
        ! receive the data
  
        CALL MPI_RECV ( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_RECV')
        ENDIF
  
        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF
  
      IF (neighbors(2) /= -1) THEN
        ! upper neighbor is present
        ! receive the data
  
        CALL MPI_RECV ( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_RECV')
                  ENDIF
  
        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF
  
      IF (neighbors(2) /= -1) THEN
        ! wait for the completion of the last send to neighbors(2)
        ! to safely reuse sendbuf(1,2)
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF
      ENDIF
  
      IF (neighbors(4) /= -1) THEN
        ! wait for the completion of the last send to neighbors(4)
        ! to safely reuse sendbuf(1,4)
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF
      ENDIF

 
!------------------------------------------------------------------------------
!- Section 4: Exchange with immediate Recv and blocking Send
!------------------------------------------------------------------------------

  ELSEIF (ntype == 2) THEN


      !------------------------------------------------------------------------
      !- Section 4.3: exchange with left and right neighbors
      !------------------------------------------------------------------------

      IF (neighbors(3) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),  &
                        ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_IRECV')
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),  &
                        ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_IRECV')
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        nzcount_ls = 0
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )

        CALL MPI_SEND ( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_SEND')
         ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        nzcount_rs = 0
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )

        CALL MPI_SEND ( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_SEND')
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF

        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,  &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF

        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      !------------------------------------------------------------------------
      !- Section 4.4: exchange with upper and lower neighbors
      !------------------------------------------------------------------------

      IF (neighbors(2) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),  &
                        ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_IRECV')
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),  &
                        ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_IRECV')
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        nzcount_ds = 0
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )

        CALL MPI_SEND ( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_SEND')
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        nzcount_us = 0
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )
 
        CALL MPI_SEND ( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_SEND')
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF

        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF

 
      IF (neighbors(4) /= -1) THEN
        ! Now wait until the data have arrived
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          CALL p_abort !(substr,'MPI_WAIT')
        ENDIF

        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF


!------------------------------------------------------------------------------
!- Section 5: Exchange with SendRecv
!------------------------------------------------------------------------------

  ELSEIF (ntype == 3) THEN


      !--------------------------------------------------------------------------
      !- Section 5.5: Send data to the left and receive from the right neighbor
      !--------------------------------------------------------------------------

      nzcount_ls = 0
      IF (MPI_neighbors(1) /= MPI_PROC_NULL) THEN
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,1), nzcount_ls,  imp_type, MPI_neighbors(1), ntag,    &
             sendbuf(1,7), isendbuflen, imp_type, MPI_neighbors(3), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(3) /= MPI_PROC_NULL) THEN
        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.6: Send data to the right and receive from the left neighbor
      !--------------------------------------------------------------------------

      nzcount_rs = 0
      IF (MPI_neighbors(3) /= MPI_PROC_NULL) THEN
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,3), nzcount_rs,  imp_type, MPI_neighbors(3), ntag,    &
             sendbuf(1,5), isendbuflen, imp_type, MPI_neighbors(1), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(1) /= MPI_PROC_NULL) THEN
        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.7: Send data to the upper and receive from the lower neighbor
      !--------------------------------------------------------------------------
 
      nzcount_us = 0
      IF (MPI_neighbors(2) /= MPI_PROC_NULL) THEN
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,2), nzcount_us,  imp_type, MPI_neighbors(2), ntag,    &
             sendbuf(1,8), isendbuflen, imp_type, MPI_neighbors(4), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(4) /= MPI_PROC_NULL) THEN
        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.8: Send data to the lower and receive from the upper neighbor
      !--------------------------------------------------------------------------
 
      nzcount_ds = 0
      IF (MPI_neighbors(4) /= MPI_PROC_NULL) THEN
        CALL putbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,4), nzcount_ds,  imp_type, MPI_neighbors(4), ntag,    &
             sendbuf(1,6), isendbuflen, imp_type, MPI_neighbors(2), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        CALL p_abort !(substr,'MPI_SENDRECV')
      ENDIF

      IF (MPI_neighbors(2) /= MPI_PROC_NULL) THEN
        CALL getbuf4 ( var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,           &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF


  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

CONTAINS
  SUBROUTINE putbuf4(var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,       &
                    ilo, iup, jlo, jup, ncount, nentry )

    USE messy_main_constants_mem, ONLY: dp
 
    ! Subroutine arguments
    INTEGER, INTENT (IN)      ::    &
         isendbuflen,                  & ! length of sendbuffer
         idim, jdim, kdim, odim,       & ! dimensions of the fields
         ilo, iup, jlo, jup,           & ! start- and end-indices
         nentry                          ! specifies the row of the sendbuf
    
    INTEGER, INTENT (INOUT)   ::    &
         ncount                          ! counts the variables
    
    REAL (dp), INTENT (INOUT) ::    &
         sendbuf (:,:)!,     & ! send buffer
    !   var01  (kdim, jdim, idim)   ! first field that has to occur
    REAL (dp), DIMENSION(:,:,:,:), POINTER :: var01
    
    ! Local variables
    INTEGER   :: i, j, k, o, nzc
    
    !------------------------------------------------------------------------------
    !- Section 2: Put data into the buffer
    !------------------------------------------------------------------------------
    
    ! use nzc as a local counter  (based on a work from Mike O'Neill to 
    ! improve vectorization of putbuf and getbuf)
    nzc = ncount
    
    ! first variable that has to be present
    DO k = 1, kdim
       DO j = jlo, jup
          DO o = 1, odim
             DO i = ilo, iup
                nzc = nzc + 1
                sendbuf (nzc,nentry) = var01(k,j,o,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 
    ! put nzc to global counter
    ncount = nzc
    
  END SUBROUTINE putbuf4

  SUBROUTINE getbuf4(var01, sendbuf, isendbuflen, idim, jdim, kdim, odim,   &
                     ilo, iup, jlo, jup, ncount, nentry )
    ! See description of getbuf (above)

    USE messy_main_constants_mem, ONLY: dp
    
    ! Subroutine arguments
    INTEGER, INTENT (IN)         ::    &
         isendbuflen,                  & ! length of sendbuffer
         idim, jdim, kdim, odim,       & ! dimensions of the fields
         ilo, iup, jlo, jup,           & ! start- and end-indices
         nentry                          ! specifies the row of sendbuf to be used
    
    INTEGER, INTENT (INOUT)      ::    &
         ncount                          ! counts the variables
    
    REAL(dp), INTENT (INOUT)            ::    &
         sendbuf (:,:)!,     & ! send buffer
    !var01  (kdim, jdim, idim)       ! first field that has to occur
    REAL (dp), DIMENSION(:,:,:,:), POINTER :: var01
    
    ! Local variables
    INTEGER    ::  i, j, k, o, nzc
    
    !------------------------------------------------------------------------------
    !- Section 2: Get data from the buffer
    !------------------------------------------------------------------------------
    
    ! use nzc as a local counter  (based on a work from Mike O'Neill to 
    ! improve vectorization of putbuf and getbuf)
    nzc = ncount
    
    ! first variable that has to be present
    DO k = 1, kdim
       DO j = jlo, jup
          DO o = 1, odim
             DO i = ilo, iup
                nzc = nzc + 1
              !  write(*,*) k,j,o,i,nzc, lbound(var01,2),ubound(sendbuf,1),ubound(sendbuf,2)
                var01(k,j,o,i) = sendbuf (nzc,nentry)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    ! put nzc to global counter
    ncount = nzc
  
  END SUBROUTINE getbuf4

#endif
END SUBROUTINE exchg_boundaries_tracer_cg

#endif

#if defined(MBM_CMAT)
!==============================================================================
  SUBROUTINE reorder2 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:)

    y = RESHAPE (x,(/SIZE(y,1),SIZE(y,2)/),(/0._dp/))

  END SUBROUTINE reorder2
!------------------------------------------------------------------------------
  SUBROUTINE reorder3 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:,:)
    INTEGER :: k

    DO k=1,SIZE(x,2)
      y(:,k,:) = RESHAPE (x(:,k,:),(/SIZE(y,1),SIZE(y,3)/),(/0._dp/))
    END DO

  END SUBROUTINE reorder3
!------------------------------------------------------------------------------
  SUBROUTINE reorder4 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:,:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:,:,:)
    INTEGER :: k, l

    DO l=1,SIZE(x,3)
      DO k=1,SIZE(x,2)
        y(:,k,l,:) = RESHAPE (x(:,k,l,:),(/SIZE(y,1),SIZE(y,4)/),(/0._dp/))
      END DO
    END DO

  END SUBROUTINE reorder4
!==============================================================================

#endif
! mz_ab_20100503-

#if (defined COSMO) || defined(BLANK)

#ifdef COSMO
  ! COSMO
  USE data_parallel,    ONLY: imp_reals, imp_grib,  imp_integers     &
                            , imp_byte, imp_character, imp_logical   &
                            , icomm_world, p_all_comm => icomm_world &
                            , p_nprocs => nproc, p_pe => my_world_id &
                            , nprocx, nprocy, isubpos, icomm_compute &
                            , num_compute, nboundlines, my_cart_id   &
                            , icomm_cart
  USE parallel_utilities, ONLY: distribute_field, gather_field         &
                              , gather_values, global_values, ij_local !&
                             ! , MPI_SUM, MPI_MAX, MPI_MIN
#endif

  ! MESSY
  USE messy_main_constants_mem, ONLY: dp, sp, i4, i8

  IMPLICIT NONE
  PUBLIC
  SAVE
#ifdef COSMO
  PUBLIC :: exchg_trac_boundaries
  PUBLIC :: exchg_trac_boundaries2
  EXTERNAL :: MPI_BCAST

#ifndef NOMPI
  INCLUDE 'mpif.h'
#endif

#endif
  
! NOTE: The following is used as a dummy in order to mimic the
!       vectorisation procedures of ECHAM5 for COSMO,
!       e.g., bi_vector, bi_decompose etc.

  TYPE decomp
     LOGICAL :: lreg = .TRUE.
  END TYPE decomp
  !
  TYPE(decomp), POINTER :: dcg
  TYPE(decomp)          :: dcl 

  INTEGER, PARAMETER :: p_io          = 0       ! CPU writing LOG-File

#ifdef BLANK
  ! DUMMY FOR (NON-)PARALLEL ENVIRONMENT
  INTEGER, PARAMETER :: p_pe         = 0       ! CPU number
  INTEGER, PARAMETER :: p_nprocs     = 1       ! number of parallel CPUs
  INTEGER, PARAMETER :: p_all_comm = 1         ! communicator dummy
  INTEGER(I4), PARAMETER :: icomm_world   = 1  ! communicator dummy
  INTEGER, PARAMETER :: imp_integers = 1       ! MPI dummy
  INTEGER, PARAMETER :: imp_reals    = 2       ! MPI dummy
  INTEGER, PARAMETER :: imp_logical  = 3       ! MPI dummy
  INTEGER, PARAMETER :: imp_character = 4      ! MPI dummy
  LOGICAL            :: p_parallel_io = .TRUE. ! TRUE if p_pe = p_io
  LOGICAL            :: p_parallel    = .FALSE.  ! parallel environment
#else
  LOGICAL            :: p_parallel_io = .FALSE. ! TRUE if p_pe = p_io
  LOGICAL            :: p_parallel    = .TRUE.  ! parallel environment
#endif
  
  INTEGER            :: itag = 88

! Global variables
#ifndef BLANK
  INTEGER, PARAMETER :: iexchg_MPI_type_len    = 200
       ! length of global vector iexchg_MPI_type used in environment.f90


  INTEGER :: iexchg_MPI_types(iexchg_MPI_type_len) = MPI_DATATYPE_NULL
                  ! List of MPI data types used in exchg_datatypes
                  ! routine. If set to MPI_DATATYPE_NULL, the
                  ! corresponding entry has not been properly set up.
  INTEGER ::  iexchg_counts(iexchg_MPI_type_len)
                  ! List of counts of these data types
                  ! used in exchg_datatypes routine.
                  ! Has meaningful value only if the corresponding
                  ! vector element of iexchg_MPI_types is not
                  ! MPI_DATATYPE_NULL.
#endif

  ! INTERFACE
  INTERFACE gather_gp
    MODULE PROCEDURE gather_gp432 ! gather gridp. field (nlon,nlev,ntrac,nlat)
                                  !                  or (nlon,nlev,nlat,1)
                                  !                  or (nlon,nlat)
    MODULE PROCEDURE gather_gp32  ! gather gridp. field (nlon,nlev,nlat)
                                  !                 or  (nlon,nlat,1)
    MODULE PROCEDURE gather_gp2   ! gather only m=0 wave number (nlon,nlat)
  END INTERFACE

  INTERFACE scatter_gp
    MODULE PROCEDURE scatter_gp432! scatter gridp. field (nlon,nlev,ntrac,nlat)
                                  !                   or (nlon,nlev,nlat,1)
    MODULE PROCEDURE scatter_gp32 ! scatter gridp. field (nlon,nlev,nlat)
                                  !                   or (nlon,nlat,1)
    MODULE PROCEDURE scatter_gp2  ! scatter gridp. field (nlon,nlat)
  END INTERFACE


  INTERFACE p_bcast
     MODULE PROCEDURE              &
          distribute_kind8,        &
          distribute_kind4,        &
          distribute_oneinteger,   &
          distribute_dp,           &
          distribute_sp,           &
          distribute_onedouble,    &
          distribute_onesingle,    &
          distribute_logical,      &
          distribute_onelogical,   &
          distribute_character,    &
          distribute_onecharacter, &
          bcast_4d
  END INTERFACE

#ifndef BLANK
  INTERFACE global_fields
     MODULE PROCEDURE      &
     global_1d_intfield,   &
     global_1d_realfield,  &
     global_2d_intfield,   &
     global_2d_realfield,  &
     global_3d_intfield,   &
     global_3d_realfield
  END INTERFACE
#endif

! SUBROUTINES
! PUBLIC :: messy_mpi_initialize
! PUBLIC :: p_abort
! PUBLIC :: bcast_4d

CONTAINS

SUBROUTINE messy_mpi_initialize
!
! Author: Astrid Kerkweg, Uni-Mainz, Mar 2008
!
! This subroutines initializes some variables needed in MESSy
!
  IF (p_nprocs > 1 ) THEN
     p_parallel = .TRUE.
  ELSE
     p_parallel = .FALSE.
  ENDIF

  p_parallel_io = (p_pe == p_io)   

END SUBROUTINE messy_mpi_initialize

!------------------------------------------------------------------------------
SUBROUTINE p_abort(pstr,qstr)
!
! Author: Astrid Kerkweg, Uni-Mainz, Mar 2008
!
! This subroutines calls the COSMO model_abort routine
!
#ifdef COSMO
  USE environment,      ONLY: model_abort
#endif

  IMPLICIT NONE
  INTRINSIC :: PRESENT, TRIM

  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: pstr
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: qstr

  CHARACTER(LEN=256) :: lpstr
  CHARACTER(LEN=256) :: lqstr

  IF (PRESENT(pstr) ) THEN
     lpstr=pstr
  ELSE
     lpstr='messy_main_mpi_bi'
  ENDIF
  IF (PRESENT(qstr) ) THEN
     lqstr=qstr
  ELSE
     lqstr='forced exit'
  ENDIF

#ifdef COSMO
  CALL model_abort(p_pe,42,TRIM(lpstr),TRIM(lqstr))
#endif
END SUBROUTINE p_abort

!-----------------------------------------------------------

SUBROUTINE bcast_4d (buffer, sender)

  IMPLICIT NONE
  INTRINSIC :: INT, SIZE

  REAL(dp), INTENT(INOUT) :: buffer(:,:,:,:)
  INTEGER,  INTENT(IN)    :: sender

  ! LOCAL
  INTEGER(i4) :: isender
  INTEGER     :: i,j,k

  isender=INT(sender,i4)

  DO i=1,SIZE(buffer,4)
     DO j=1,SIZE(buffer,3)
        DO k=1,SIZE(buffer,2)
           CALL p_bcast(buffer(:,k,j,i), isender)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE bcast_4d

!------------------------------------------------------------------------------
!
! SUBROUTINE p_bcast(buffer, isender, icomm, ibufferlen, idatatype, 
!                                ierrorcode)
!
!------------------------------------------------------------------------------
!
! Description:
!  p_bcast is a generic name for several subroutines that distribute
!  values from one processor to all others. Depending on the type of the 
!  first argument, the appropriate procedure is chosen.
!
! Method:
!  With the MPI_BCAST command the buffer is distributed to all other 
!  processors.
!  
!  This routines are -in principle- copies of the COSMO routines 
!  distribute_values (parallel_utilities.f90) but they arguments that are
!  not given by the p_bcast calls in the MESSy submodel interface layer are
!  made OPTIONAL parameter here.
!
!------------------------------------------------------------------------------

! Following are the different subroutines

!------------------------------------------------------------------------------

!==============================================================================

!+ Subroutine for array of kind=8 integers

SUBROUTINE distribute_kind8 (buffer, isender, ibufferlen, icomm,  &
                             idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
INTEGER (KIND=i8), INTENT(INOUT) :: buffer(:) ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender   ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):

INTEGER (KIND=i4), INTENT(OUT), OPTIONAL ::                         &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4)  :: p_ierrorcode ! local error code
INTEGER (KIND=i4) :: p_ibufferlen  ! local buffer length           
INTEGER (KIND=i4) :: p_idatatype   ! local data type         
INTEGER (KIND=i4) :: p_icomm       ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_integers 
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,              &
                  p_icomm, p_ierrorcode)  
#else
p_ierrorcode = 0
#endif
IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_kind8

!==============================================================================

!==============================================================================

!+ Subroutine for array of kind=4 integers

SUBROUTINE distribute_kind4 (buffer, isender, ibufferlen, icomm,  &
                             idatatype, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE

INTEGER (KIND=i4), INTENT(INOUT) :: buffer(:) ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender   ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL ::                         &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4)  :: p_ierrorcode       ! local error code
INTEGER (KIND=i4) :: p_ibufferlen       ! local buffer length           
INTEGER (KIND=i4) :: p_idatatype        ! local data type         
INTEGER (KIND=i4) :: p_icomm    ! local communicator
!- End of header
!------------------------------------------------------------------------------
IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_integers 
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,                 &
                  p_icomm, p_ierrorcode)  
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_kind4

!==============================================================================

!==============================================================================

!+ Subroutine for one model integers

SUBROUTINE distribute_oneinteger(buffer, isender, ibufferlen, icomm,  &
                             idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
INTEGER (KIND=i4),  INTENT(INOUT) :: buffer  ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)     :: isender ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

INTEGER (KIND=i4), INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode  ! local error code
INTEGER (KIND=i4) :: p_idatatype   ! local data type         
INTEGER (KIND=i4) :: p_icomm       ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_integers 
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST ( buffer, 1, p_idatatype, isender, p_icomm &
                 , p_ierrorcode)  
#else
p_ierrorcode = 0
#endif
IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_oneinteger

!==============================================================================

!==============================================================================

!+ Subroutine for array of doubles

SUBROUTINE distribute_dp (buffer, isender, ibufferlen, icomm,  &
                             idatatype, ierrorcode) 

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
REAL (KIND=dp),  INTENT(INOUT) :: buffer(:)  ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)  :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen, & ! length of the buffer
  idatatype,  & ! type of buffer
  icomm         ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode    ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode ! local error code
INTEGER (KIND=i4) :: p_ibufferlen ! local buffer length           
INTEGER (KIND=i4) :: p_idatatype  ! local data type         
INTEGER (KIND=i4) :: p_icomm      ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_reals 
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm, p_ierrorcode)
#else
p_ierrorcode = 0
#endif
IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_dp

!==============================================================================

!==============================================================================

!+ Subroutine for one double

SUBROUTINE distribute_onedouble (buffer, isender, ibufferlen, icomm,  &
                             idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
REAL (KIND=dp),  INTENT(INOUT) :: buffer     ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)  :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode       ! local error code
INTEGER (KIND=i4) :: p_ibufferlen       ! local buffer length           
INTEGER (KIND=i4) :: p_idatatype        ! local data type         
INTEGER (KIND=i4) :: p_icomm    ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = 1
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_reals 
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
 CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm, p_ierrorcode)
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_onedouble 

!==============================================================================

!==============================================================================

!+ Subroutine for array of singles

SUBROUTINE distribute_sp (buffer, isender, ibufferlen  &
                        , icomm, idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
REAL (KIND=sp),    INTENT(INOUT) ::  buffer(:) ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode     ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode  ! local error code
INTEGER (KIND=i4) :: p_ibufferlen  ! local buffer length           
INTEGER (KIND=i4) :: p_idatatype   ! local data type         
INTEGER (KIND=i4) :: p_icomm       ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_reals 
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm, p_ierrorcode)
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_sp

!==============================================================================

!==============================================================================

!+ Subroutine for one single

SUBROUTINE distribute_onesingle (buffer, isender, ibufferlen  &
                              , icomm, idatatype, ierrorcode)
  IMPLICIT NONE

  INTRINSIC :: PRESENT

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
REAL (KIND=sp),    INTENT(INOUT) :: buffer     ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode   ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode   ! local error code
INTEGER (KIND=i4) :: p_ibufferlen   ! local buffer length           
INTEGER (KIND=i4) :: p_idatatype    ! local data type         
INTEGER (KIND=i4) :: p_icomm        ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = 1
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_reals 
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm, p_ierrorcode)
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_onesingle

!==============================================================================

!==============================================================================

!+ Subroutine for array of default logicals

SUBROUTINE distribute_logical  (buffer, isender, ibufferlen  &
                              , icomm, idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
LOGICAL,           INTENT(INOUT) :: buffer(:)  ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode  ! local error code
INTEGER (KIND=i4) :: p_ibufferlen  ! local buffer length           
INTEGER (KIND=i4) :: p_idatatype   ! local data type         
INTEGER (KIND=i4) :: p_icomm       ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_logical
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm, p_ierrorcode)  
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
  ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_logical

!==============================================================================

!==============================================================================

!+ Subroutine for one default logical

SUBROUTINE distribute_onelogical (buffer, isender, ibufferlen  &
                                , icomm, idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
LOGICAL,           INTENT(INOUT) :: buffer     ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL ::  ierrorcode   ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode ! local error code
INTEGER (KIND=i4) :: p_idatatype  ! local data type         
INTEGER (KIND=i4) :: p_icomm      ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_logical
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF
#ifdef COSMO
 CALL MPI_BCAST (buffer, 1, p_idatatype, isender, p_icomm &
               , p_ierrorcode)  
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_onelogical

!==============================================================================

!==============================================================================

!+ Subroutine for array of characters

SUBROUTINE distribute_character (buffer, isender, ibufferlen  &
                               , icomm, idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: CHAR, ICHAR, PRESENT, SIZE
#ifdef COSMO
  EXTERNAL  :: MPI_COMM_RANK
#endif

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments


! Array arguments with intent(inout):
! op_pj_20110313+
!!$CHARACTER (LEN=100), INTENT(INOUT) :: buffer(:)  ! buffer to be broadcasted
CHARACTER (LEN=*), INTENT(INOUT) :: buffer(:)  ! buffer to be broadcasted
! op_pj_20110313-

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)      :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode  ! error code

!------------------------------------------------------------------------------

! Local Scalars and Arrays
INTEGER (KIND=i4) :: p_ierrorcode       ! local error code
INTEGER (KIND=i4) :: p_ibufferlen       ! local buffer length           
INTEGER (KIND=i4) :: p_icomm            ! local communicator
INTEGER (KIND=i4) ::  my_comm_id, i, j

INTEGER    :: intbuf(100)      ! Standard integer
!- End of header
!------------------------------------------------------------------------------


 IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF
#ifdef COSMO
  CALL MPI_COMM_RANK(p_icomm, my_comm_id, p_ierrorcode)
#else 
  my_comm_id=0
#endif
  DO i=1,p_ibufferlen
    IF (my_comm_id == isender) THEN
       ! um_ak_20110721+
       !DO j=1,100
       DO j=1,LEN(buffer(i))
       ! um_ak_20110721-
        intbuf(j) = ICHAR ( buffer(i)(j:j) )
      ENDDO
    ENDIF

#ifdef COSMO
    CALL MPI_BCAST (intbuf, 100, imp_integers, isender  &
                  , p_icomm, p_ierrorcode)
#else
    p_ierrorcode = 0
#endif

    IF (my_comm_id /= isender ) THEN
      !um_ak_20110721+
      !DO j=1,100
      DO j=1,LEN(buffer(i))
      !um_ak_20110721-
        buffer(i)(j:j) = CHAR (intbuf(j) )
      ENDDO
    ENDIF
  ENDDO

! and this would be the normal way
! CALL MPI_BCAST (buffer, ibufferlen, MPI_CHARACTER, isender,   &
!                 icomm, implcode)  

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = p_ierrorcode
  ENDIF

END SUBROUTINE distribute_character

!==============================================================================

!==============================================================================

!+ Subroutine for one word of characters

SUBROUTINE distribute_onecharacter (buffer, isender, ibufferlen  &
                              , icomm, idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: LEN, PRESENT

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments

! Array arguments with intent(inout):
CHARACTER (LEN=*), INTENT(INOUT) ::  buffer     ! character to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors


! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars and Arrays
INTEGER (KIND=i4) :: p_ierrorcode ! local error code
INTEGER (KIND=i4) :: p_idatatype  ! local data type         
INTEGER (KIND=i4) :: p_icomm      ! local communicator

!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_character
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
 CALL MPI_BCAST (buffer, len(buffer), p_idatatype, isender,              &
                 p_icomm, p_ierrorcode)  
#else
 p_ierrorcode = 0
#endif

 IF (PRESENT(ierrorcode)) THEN
    ierrorcode = p_ierrorcode
 ENDIF

END SUBROUTINE distribute_onecharacter

!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================


!==============================================================================
  SUBROUTINE gather_gp432 (gl, lc, gl_dc, source, lg32)

  !
  ! gather global grid point field from pe's (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   ) 
  !                                       or (nlon,nlat,1    ,1   )
  !

    IMPLICIT NONE
    INTRINSIC :: SIZE

    ! INPUT
    REAL(dp),             POINTER     :: gl   (:,:,:,:) ! global field
    REAL(dp), TARGET,     INTENT(in)  :: lc   (:,:,:,:) ! local  field
    TYPE(decomp),OPTIONAL,INTENT(in)  :: gl_dc          ! global decomposition
    INTEGER, OPTIONAL    ,INTENT(in)  :: source         ! source to gather from
    !                                                 ! -1=all;0=p_io;1=not p_io
    LOGICAL, OPTIONAL,    INTENT(in)  :: lg32           ! um_ak_20110603
    ! LOCAL
    INTEGER(i4) :: size4 ! size of 4th dimension
    !
    INTEGER :: i
    REAL(dp), POINTER :: gl3(:,:,:)
    !
    ! call 2D gather routine if 4th dimension size is 1
    ! else loop over 3th index
    !
#ifdef COSMO
    IF (p_parallel_io) size4 = SIZE(gl,4) 
    CALL p_bcast (size4, p_io)
    DO i=1,size4
       NULLIFY (gl3)
       gl3 => gl(:,:,:,i)
       CALL gather_gp32 (gl3, lc(:,:,:,i), gl_dc, &
            source=source, lg32=lg32)! um_ak_20110603 lg32 added
    END DO
#endif
    !
  END SUBROUTINE gather_gp432
!------------------------------------------------------------------------------
  SUBROUTINE gather_gp32 (gl, lc, gl_dc,source, lg32)
  !
  ! gather global grid point field from pe's (nlon,nlev,nlat) or (nlon,nlat,1)
  !
    IMPLICIT NONE

    INTRINSIC :: SIZE

    REAL(dp)             ,POINTER     :: gl   (:,:,:) ! global field
    REAL(dp)    , TARGET ,INTENT(in)  :: lc   (:,:,:) ! local  field
    TYPE(decomp),OPTIONAL,INTENT(in)  :: gl_dc        ! global decomposition
    INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                                 ! -1=all;0=p_io;1=not p_io
    LOGICAL, OPTIONAL    , INTENT(in) :: lg32         ! um_ak_20110603

    ! LOCAL
    INTEGER(i4) :: isize ! size of 3rd dimension
    INTEGER     :: i
    REAL(dp),POINTER :: gl2(:,:)
    LOGICAL     :: l3
    !
    ! call 2D gather routine if 3rd dimension size is 1
    ! else call 3D gather routine
    !
#ifdef COSMO
    l3=.TRUE.
    IF (PRESENT(lg32)) THEN
       IF (lg32) l3=.FALSE. 
    ENDIF

    IF (l3) THEN
       IF (p_parallel_io) isize = SIZE(gl,3)
       CALL p_bcast(isize, p_io)
       NULLIFY (gl2)
       DO i=1,isize 
          !IF (p_pe == p_io) gl2 => gl(:,:,i)
          gl2 => gl(:,:,i)
          CALL gather_gp2 (gl2, lc(:,:,i))
       ENDDO
    ELSE
       IF (p_parallel_io) isize = SIZE(gl,2)
       CALL p_bcast(isize, p_io)
       NULLIFY (gl2)
       DO i=1,isize
          !IF (p_pe == p_io) gl2 => gl(:,:,i)
          gl2 => gl(:,i,:)
          CALL gather_gp2 (gl2, lc(:,i,:))
       ENDDO
       
    ENDIF

#endif
  END SUBROUTINE gather_gp32
!------------------------------------------------------------------------------
  SUBROUTINE gather_gp2 (gl, lc)
 
    IMPLICIT NONE
    INTRINSIC :: SIZE

    REAL(dp),POINTER              :: gl   (:,:) ! global field
    REAL(dp), TARGET ,INTENT(in)  :: lc   (:,:) ! local  field
    INTEGER(i4) :: l_size(2), g_size(2)
    INTEGER     :: ierror = 0
    !
#ifdef COSMO
    IF (p_pe == p_io) THEN
       g_size = (/ SIZE(gl,1), SIZE(gl,2) /)
    END IF
    CALL p_bcast (g_size, p_io)
    l_size = (/ SIZE(lc,1), SIZE(lc,2) /)
  
    CALL gather_field(lc,l_size(1),l_size(2),gl,g_size(1),g_size(2) &
         ,p_io, ierror)
    
    IF (ierror /= 0) THEN
       CALL p_abort('gather_gp2','messy_main_mpi_bi')
       RETURN
    ENDIF
#endif    
  END SUBROUTINE gather_gp2
!==============================================================================


!==============================================================================
  SUBROUTINE scatter_gp432 (gl, lc, gl_dc, lg32)
  !
  ! scatter global grid point field from pe's (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   ) 
  !                                       or (nlon,nlat,1    ,1   )
  !
    IMPLICIT NONE

    INTRINSIC :: SIZE

    REAL(dp), POINTER                 :: gl   (:,:,:,:) ! global field
    REAL(dp), TARGET ,INTENT(out)     :: lc   (:,:,:,:) ! local  field
    TYPE(DECOMP),OPTIONAL,INTENT(in)  :: gl_dc          ! global decomposition
    LOGICAL, OPTIONAL,    INTENT(in)  :: lg32           ! um_ak_20110603
    ! LOCAL
    INTEGER(i4)      :: size4 ! size of 4th dimension
    INTEGER          :: i
    REAL(dp),POINTER :: gl3(:,:,:)
    !
    ! call 3D scatter routine if 4th dimension size is 1
    ! else loop over 3th index
    !
#ifdef COSMO
    IF (p_parallel_io) size4 = (SIZE(gl,4))
    CALL p_bcast (size4, p_io)
    NULLIFY(gl3)
    DO i=1,size4
!       IF (p_pe == p_io) gl3 => gl(:,:,:,i)
       gl3 => gl(:,:,:,i)
       CALL scatter_gp32 (gl3, lc(:,:,:,i), gl_dc, lg32=lg32)! um_ak_20110603 lg32 added)
    END DO
#else
    lc(:,:,:,:) = gl(1:SIZE(lc,1),1:SIZE(lc,2),1:SIZE(lc,3),1:SIZE(lc,4))
#endif
    ! 
  END SUBROUTINE scatter_gp432
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp32 (gl, lc, gl_dc, lg32)

    IMPLICIT NONE

    INTRINSIC :: SIZE
  !
  ! send global grid point field to pe's (nlon,nlev,nlat) or (nlon,nlat,1)
  !

  REAL(dp), POINTER                 :: gl   (:,:,:) ! global field
  REAL(dp), TARGET ,INTENT(out)     :: lc   (:,:,:) ! local  field
  TYPE(decomp),OPTIONAL,INTENT(in)  :: gl_dc        ! global decomposition
  LOGICAL, OPTIONAL    , INTENT(in) :: lg32         ! um_ak_20110603

  ! LOCAL
  REAL(dp),POINTER :: gl2(:,:)
  INTEGER(i4)      :: isize    ! size of 3rd dimension
  INTEGER          :: i
  LOGICAL          :: l3
  !
  ! call 2D scatter routine if 3rd dimension size is 1
  ! else call 3D scatter routine
  !
#ifdef COSMO
    l3=.TRUE.
    IF (PRESENT(lg32)) THEN
       IF (lg32) l3=.FALSE. 
    ENDIF

    IF (l3) THEN
       IF (p_pe == p_io) isize = (SIZE(gl,3))
       CALL p_bcast (isize, p_io)
       DO i=1,isize
          NULLIFY(gl2)
          !       IF (p_pe == p_io) gl2 => gl(:,:,i)
          gl2 => gl(:,:,i)
          CALL scatter_gp2 (gl2, lc(:,:,i), gl_dc)
       ENDDO
    ELSE
       IF (p_pe == p_io) isize = (SIZE(gl,2))
       CALL p_bcast (isize, p_io)
       DO i=1,isize
          NULLIFY(gl2)
          !       IF (p_pe == p_io) gl2 => gl(:,:,i)
          gl2 => gl(:,i,:)
          CALL scatter_gp2 (gl2, lc(:,i,:), gl_dc)
       ENDDO
    ENDIF
#else
    lc(:,:,:) = gl(1:SIZE(lc,1),1:SIZE(lc,2),1:SIZE(lc,3))
#endif
    !
  END SUBROUTINE scatter_gp32
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp2 (gl, lc, gl_dc)
  !
  ! send global 2D grid point field to local pe's (nlon,nlat)
  !
    IMPLICIT NONE
    INTRINSIC :: SIZE

  REAL(dp)                 ,POINTER     :: gl   (:,:) ! global field
  REAL(dp)         ,TARGET ,INTENT(out) :: lc   (:,:) ! local  field
  TYPE(decomp),   OPTIONAL ,INTENT(in)  :: gl_dc !global decomposition
    !
  INTEGER(i4) :: l_size(2), g_size(2)
  INTEGER(i4) :: ierror = 0
  !
#ifdef COSMO
  IF (p_pe == p_io) THEN
     g_size = (/ SIZE(gl,1), SIZE(gl,2) /)
  END IF
  CALL p_bcast (g_size, p_io)
  l_size = (/ SIZE(lc,1), SIZE(lc,2) /)

  CALL distribute_field(gl,g_size(1),g_size(2),lc,l_size(1),l_size(2) &
       , ierror )

  IF (ierror /= 0) THEN
     CALL p_abort('scatter_gp2','messy_main_mpi_bi')
     RETURN
  ENDIF
#else
    lc(:,:) = gl(1:SIZE(lc,1),1:SIZE(lc,2))
#endif
  END SUBROUTINE scatter_gp2
!------------------------------------------------------------------------------
!==============================================================================
  SUBROUTINE reorder (y,x)

    IMPLICIT NONE
    INTRINSIC :: RESHAPE, SIZE

    REAL(dp) ,INTENT(out) :: y (:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:)

#if (defined CRAY) || (defined sun) || (defined NAG) || (defined __SX__) || defined(__PGI)
    CALL util_reshape(y, x, SIZE(y,1)*SIZE(y,2), SIZE(x,1)*SIZE(x,2))
#else
    y = RESHAPE (x,(/SIZE(y,1),SIZE(y,2)/),(/0._dp/))
#endif

  END SUBROUTINE reorder
!------------------------------------------------------------------------------

#ifdef I2CINC
! ----------------------------------------------------------------------
  SUBROUTINE switch_par_utilities(flag)

    ! SWITCH BETWEEN COSMO AND INT2COSMO DIMENSIONS OF FIELDS
    ! FOR parallel_utilities

    ! INT2COSMO
    USE data_grid_lm,         ONLY: istartpar_i2c   => istartpar    &
                                  , iendpar_i2c     => iendpar      &
                                  , jstartpar_i2c   => jstartpar    &
                                  , jendpar_i2c     => jendpar      &
                                  , ie2lm_tot, je2lm_tot, kelm_tot  &
                                  , kelm, ie2lm_max, je2lm_max      &
                                  , ie2lm, je2lm
    
    USE data_int2lm_parallel, ONLY: nboundlines_i2c => nboundlines &
                                  , isubpos_i2c     => isubpos

    ! COSMO
    USE data_modelconfig, ONLY: istartpar_c4   => istartpar &
                              , iendpar_c4     => iendpar   &  
                              , jstartpar_c4   => jstartpar &
                              , jendpar_c4     => jendpar   &
                              , ke, ke_tot, ie, je, ie_tot  &
                              , je_tot, ie_max, je_max             
    USE data_parallel,    ONLY: nboundlines_c4 => nboundlines &
                              , isubpos_c4     => isubpos     &
    ! COSMO/INT2COSMO
                              , icomm_cart, imp_integers       &
                              , nprocx, nprocy, nproc, nprocio &
                              , my_cart_id

    USE parallel_utilities, ONLY: init_par_utilities

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='switch_par_utilities'

    SELECT CASE (flag)
    CASE(1)
       ! SWITCH FROM COSMO TO INT2COSMO
       CALL init_par_utilities (ie2lm, je2lm, kelm     &
            , ie2lm_tot, je2lm_tot, kelm_tot           &
            , ie2lm_max, je2lm_max                     &
            , istartpar_i2c, iendpar_i2c               &
            , jstartpar_i2c, jendpar_i2c               &
            , nproc, nprocx, nprocy, nprocio           &
            , isubpos_i2c, nboundlines_i2c, icomm_cart &
            , my_cart_id, imp_reals, imp_integers      )
    CASE(2)
       ! SWITCH FROM INT2COSMO TO COSMO
       CALL init_par_utilities (ie, je, ke, ie_tot, je_tot, ke_tot     &
            , ie_max, je_max, istartpar_c4, iendpar_c4                 &
            , jstartpar_c4, jendpar_c4, nproc, nprocx, nprocy, nprocio &
            , isubpos_c4, nboundlines_c4, icomm_cart, my_cart_id       &
            , imp_reals, imp_integers)

    CASE DEFAULT
       ! NO FURTHER SWITCHING OPTIONS
       write (0,*) ' THIS SWITCHING OPTION IS NOT DEFINED! '
    END SELECT

  END SUBROUTINE switch_par_utilities

! -------------------------------------------------------------------------
#endif

#ifndef BLANK
! -------------------------------------------------------------------------
  SUBROUTINE global_1d_intfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)  :: dims
    INTEGER, INTENT(IN)  :: infield(:)
    INTEGER, INTENT(OUT) :: outfield(:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF
    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_integers, noper, icomm_cart, ierror)

  END SUBROUTINE global_1d_intfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE global_2d_intfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)  :: dims
    INTEGER, INTENT(IN)  :: infield(:,:)
    INTEGER, INTENT(OUT) :: outfield(:,:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF

    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_integers, noper, icomm_cart, ierror)

  END SUBROUTINE global_2d_intfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE global_3d_intfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)  :: dims
    INTEGER, INTENT(IN)  :: infield(:,:,:)
    INTEGER, INTENT(OUT) :: outfield(:,:,:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF

    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_integers, noper, icomm_cart, ierror)

  END SUBROUTINE global_3d_intfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE global_1d_realfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)   :: dims
    REAL(dp), INTENT(IN)  :: infield(:)
    REAL(dp), INTENT(OUT) :: outfield(:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF

    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_reals, noper, icomm_cart, ierror)

  END SUBROUTINE global_1d_realfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE global_2d_realfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)   :: dims
    REAL(dp), INTENT(IN)  :: infield(:,:)
    REAL(dp), INTENT(OUT) :: outfield(:,:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF

    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_reals, noper, icomm_cart, ierror)

  END SUBROUTINE global_2d_realfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE global_3d_realfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)   :: dims
    REAL(dp), INTENT(IN)  :: infield(:,:,:)
    REAL(dp), INTENT(OUT) :: outfield(:,:,:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF

    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_reals, noper, icomm_cart, ierror)


  END SUBROUTINE global_3d_realfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE exchange_boundaries(kdims, field, ierror, yerrmsg )

    USE environment,      ONLY: exchg_boundaries
    USE data_parallel,    ONLY: ldatatypes, ncomm_type, my_cart_neigh &
                              , sendbuf, isendbuflen, icomm_cart, nboundlines
    USE data_modelconfig, ONLY: jstartpar, jendpar, ie, je
    USE data_runcontrol,  ONLY: nnow

    IMPLICIT NONE

    INTEGER, INTENT(IN)                       :: kdims(24)
    REAL(dp), INTENT(INOUT), DIMENSION(:,:,:) :: field
    INTEGER, INTENT(OUT)                      :: ierror
    CHARACTER(LEN=*), INTENT(OUT)             :: yerrmsg
    
    itag = itag+1

    CALL exchg_boundaries                                                   &
             (nnow+39, sendbuf, isendbuflen, imp_reals, icomm_cart, ie, je  &
             , kdims, jstartpar, jendpar, nboundlines, nboundlines          &
             , my_cart_neigh, itag, ldatatypes, ncomm_type                  &
             , ierror, yerrmsg                                             &
             , field(:,:,:))


  END SUBROUTINE exchange_boundaries
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------

  SUBROUTINE exchg_trac_boundaries2 ( idim, jdim, kdim, trac_dim,  &
                 jstartpar, jendpar, nlines, nboundlines,         &
                 neighbors, ntag, ierror, yerrmsg, trac_field)

    !----------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine performs the boundary exchange for the full tracer field.
    !
    !----------------------------------------------------------------------------

    IMPLICIT NONE

#ifndef NOMPI
  INCLUDE 'mpif.h'
#endif

    ! Subroutine arguments
    INTEGER, INTENT(IN)         ::    &
         idim, jdim,         & ! horizontal dimensions of the fields
         kdim,               & ! vertical dimensions  of the tracer field
         trac_dim,           & ! number dimension of the tracer field
         jstartpar,          & ! start index in j-direction
         jendpar,            & ! end index in j-direction
         nlines,             & ! number of lines that have to be exchanged
                               ! (<= nboundlines)
         nboundlines,        & ! number of overlapping boundary lines
         neighbors(4),       & ! process-id's of the neighbors in the grid
         ntag                  ! tag of the message

    INTEGER, INTENT (OUT) :: ierror    ! error status variable
    
    CHARACTER (LEN=*), INTENT(OUT)  :: yerrmsg ! for MPI error message
    
    REAL (dp), DIMENSION(:,:,:,:), POINTER :: trac_field 
    
    ! LOCAL
    INTEGER :: izmplcode            ! for MPI error code
    INTEGER :: i,j,k, jt, ind, datasize
    INTEGER :: ireq(4)
    INTEGER :: status (MPI_STATUS_SIZE)
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: sbuffer
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: rbuffer
  
!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ierror     = 0
  izmplcode  = 0
  yerrmsg    = '    '

  ! check whether nlines <= nboundlines
  IF (nlines > nboundlines) THEN
    ierror  = 9011
    yerrmsg     = ' *** nlines > nboundlines *** '
    RETURN
  ENDIF

  ! Fix list of neighbors (use MPI_PROC_NULL rather than -1 to indicate
  ! missing neighbor).
  ! Horizontal exchange

  datasize = nlines * (jendpar - jstartpar + 1) * kdim * trac_dim

  ALLOCATE (sbuffer(datasize))
  ALLOCATE (rbuffer(datasize))

  ! left neighbor is present
  IF (neighbors(1) /= -1) THEN
     ind = 1
     DO jt = 1, trac_dim
        DO i = nboundlines+1, nboundlines + nlines
           DO j = jstartpar, jendpar
              DO k= 1, kdim
                 sbuffer(ind) = trac_field(i,j,jt,k)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO

     CALL MPI_ISEND (sbuffer, datasize,  imp_reals, &
          neighbors(1), ntag, icomm_cart, ireq(1), ierror)
     IF (ierror /= 0) RETURN

  END IF

  ! right neighbor exists
  IF (neighbors(3) /= -1) THEN
     ! 2.) send buffer
     ind = 1
     DO jt = 1, trac_dim
        DO i = idim - nboundlines - nlines + 1, idim-nboundlines
           DO j = jstartpar, jendpar
              DO k= 1, kdim
                 sbuffer(ind) = trac_field(i,j,jt,k)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
     CALL MPI_ISEND (sbuffer, datasize,  imp_reals, &
          neighbors(3), ntag+1, icomm_cart, ireq(3), ierror)
     IF (ierror /= 0) RETURN 

     ! 1.) receive buffer of right neighbor 
     rbuffer(:) = 0._dp

      CALL MPI_RECV(rbuffer, datasize,  imp_reals &
          , neighbors(3), ntag, icomm_cart, MPI_STATUS_IGNORE, ierror)
     IF (ierror /= 0) RETURN

     ind = 1
     DO jt = 1, trac_dim
        DO i = idim - nboundlines + 1,idim - nboundlines + nlines
           DO j = jstartpar, jendpar
              DO k= 1, kdim
                 trac_field(i,j,jt,k) = rbuffer(ind)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO

  END IF
 
  ! left neighbors exists recv buffer
  IF (neighbors(1) /= -1) THEN

     rbuffer(:) = 0._dp

     CALL MPI_RECV(rbuffer, datasize,  imp_reals &
          , neighbors(1), ntag+1, icomm_cart, MPI_STATUS_IGNORE, ierror)
     IF (ierror /= 0) RETURN

     ind = 1
     DO jt = 1, trac_dim
        DO i = nboundlines + 1 - nlines, nboundlines 
           DO j = jstartpar, jendpar
              DO k= 1, kdim
                 trac_field(i,j,jt,k) = rbuffer(ind)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
  END IF
  
  IF (neighbors(1) /= -1) THEN
     ! wait for the completion of the last send to neighbors(1)
     CALL MPI_WAIT (ireq(1), status, izmplcode)
     IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_WAIT 1'
        RETURN
     ENDIF
  ENDIF
  IF (neighbors(3) /= -1) THEN
     ! wait for the completion of the last send to neighbors(1)
     CALL MPI_WAIT (ireq(3), status, izmplcode)
     IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_WAIT 3'
        RETURN
     ENDIF
  ENDIF
  


  DEALLOCATE (rbuffer, sbuffer)

  ! *************************************************
  ! EXCHANGE with northern and southern neighbors
  ! *************************************************

  ! DIMENSION buffer
  datasize =  nlines * idim * kdim * trac_dim
  ALLOCATE (sbuffer(datasize))
  ALLOCATE (rbuffer(datasize))
  
  ! upper neighbor is present
  IF (neighbors(2) /= -1) THEN
      ! send buffer
     ind = 1
     DO jt = 1, trac_dim
        DO i = 1 , idim
           DO j = jdim - nboundlines - nlines + 1, jdim - nboundlines
              DO k= 1, kdim
                 sbuffer(ind) = trac_field(i,j,jt,k)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
     CALL MPI_ISEND (sbuffer, datasize,  imp_reals, &
          neighbors(2), ntag+2, icomm_cart, ireq(2), ierror)
     IF (ierror /= 0) RETURN
 
  ENDIF

  ! lower neighbor is present
  IF (neighbors(4) /= -1) THEN
      ! send buffer to lower neighbor
     ind = 1
     DO jt = 1, trac_dim
        DO i = 1 , idim
           DO j = nboundlines + 1, nboundlines + nlines
              DO k= 1, kdim
                 sbuffer(ind) = trac_field(i,j,jt,k)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
     CALL MPI_ISEND (sbuffer, datasize,  imp_reals, &
          neighbors(4), ntag+3, icomm_cart, ireq(4), ierror)
     IF (ierror /= 0) RETURN
    
     rbuffer(:) = 0._dp

     ! receive buffer of lower neighbor
     CALL MPI_RECV(rbuffer, datasize,  imp_reals &
          , neighbors(4), ntag+2, icomm_cart, MPI_STATUS_IGNORE, ierror)
     IF (ierror /= 0) RETURN

     ind = 1
     DO jt = 1, trac_dim
        DO i = 1, idim
           DO j = nboundlines-nlines+1, nboundlines
              DO k= 1, kdim
                 trac_field(i,j,jt,k) = rbuffer(ind)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
  ENDIF

  ! upper neighbor is present
  IF (neighbors(2) /= -1) THEN
     ! receive buffer of upper neighbor
     rbuffer(:) = 0._dp

     CALL MPI_RECV(rbuffer, datasize,  imp_reals &
          , neighbors(2), ntag+3, icomm_cart, MPI_STATUS_IGNORE, ierror)
     IF (ierror /= 0) RETURN

     ind = 1
     DO jt = 1, trac_dim
        DO i = 1, idim
!           DO j = jdim -nboundlines - nlines +1 , jdim - nboundlines
           DO j = jdim -nboundlines + 1 , jdim - nboundlines + nlines
              DO k= 1, kdim
                 trac_field(i,j,jt,k) = rbuffer(ind)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
  ENDIF

  IF (neighbors(2) /= -1) THEN
     ! wait for the completion of the last send to neighbors(1)
     write (0,*) 'WAIT2 01'
     CALL MPI_WAIT (ireq(2), status, izmplcode)
     write (0,*) 'WAIT2 01'
     IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_WAIT 2'
        RETURN
     ENDIF
  ENDIF
  IF (neighbors(4) /= -1) THEN
     ! wait for the completion of the last send to neighbors(1)
     write (0,*) 'WAIT4 01'
     CALL MPI_WAIT (ireq(4), status, izmplcode)
     write (0,*) 'WAIT4 02'
     IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_WAIT 4'
        RETURN
     ENDIF
  ENDIF
   DEALLOCATE (rbuffer, sbuffer)


END SUBROUTINE exchg_trac_boundaries2

!==============================================================================
!==============================================================================
!+ This subroutine performs the data exchange between boundaries
!------------------------------------------------------------------------------

SUBROUTINE exchg_trac_boundaries                                             &
               ( icase, imp_type, icomm, idim, jdim,   &
                 kdim, trac_dim, jstartpar, jendpar, nlines, nboundlines,    &
                 neighbors, ntag, lmpi_types, ntype, ierror, yerrmsg,        &
                 trac_field)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine performs the boundary exchange of up to 20 variables. Only
!   one variable has to occur, the others are optional. 
!
! Method:
!   At the moment there are 3 different MPI-communications implemented:
!     1) immediate send, blocking receive and wait
!     2) immediate receive, blocking send and wait
!     3) MPI_SendRecv
!   Also there is the choice of an explicit buffering (putbuf, getbuf) or 
!   implicit buffering (MPI-Datatypes) of the data to be send.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    icase,              & ! tag for exchange scenario
!    isendbuflen,        & ! length of sendbuffer
    imp_type,           & ! determines the REAL type used
    icomm,              & ! communicator for virtual cartesian topology
    idim, jdim,         & ! horizontal dimensions of the fields
    kdim,               & ! array for the vertical dimensions of var1..var20
    trac_dim,           & ! number of tracers
    jstartpar,          & ! start index in j-direction
    jendpar,            & ! end index in j-direction
    nlines,             & ! number of lines that have to be exchanged
                          ! (<= nboundlines)
    nboundlines,        & ! number of overlapping boundary lines
    neighbors(4),       & ! process-id's of the neighbors in the grid
    ntag,               & ! tag of the message
    ntype                 ! indicates how the communication should be done

  LOGICAL, INTENT(IN) :: lmpi_types ! whether implicit (with MPI-Datatypes) 
                                    ! or explicit
                                    ! (putbuf, getbuf) buffering of data is used

  INTEGER,         INTENT (OUT) :: ierror      ! error status variable

  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg     ! for MPI error message

  REAL(dp),       INTENT (INOUT)      ::    &
!    sendbuf (isendbuflen, 8),& ! send buffer
    trac_field(idim,jdim,trac_dim,kdim)


  ! LOCAL
  INTEGER ::       &
    ! the following numbers are for filling/emptying the buffers for each
    ! neighbor
    izlo_lr, izup_lr, jzlo_lr, jzup_lr,     & ! left , receive
    izlo_rr, izup_rr, jzlo_rr, jzup_rr,     & ! right, receive
    izlo_ur, izup_ur, jzlo_ur, jzup_ur,     & ! upper, receive
    izlo_dr, izup_dr, jzlo_dr, jzup_dr,     & ! down , receive
    izlo_ls, izup_ls, jzlo_ls, jzup_ls,     & ! left , send
    izlo_rs, izup_rs, jzlo_rs, jzup_rs,     & ! right, send
    izlo_us, izup_us, jzlo_us, jzup_us,     & ! upper, send
    izlo_ds, izup_ds, jzlo_ds, jzup_ds,     & ! down , send
    nzcount_ls, nzcount_rs,     & ! counting the values
    nzcount_us, nzcount_ds,     & ! counting the values
    nzcount_lr, nzcount_rr,     & ! counting the values
    nzcount_ur, nzcount_dr,     & ! counting the values
    nzrequest(MPI_STATUS_SIZE), & ! for MPI-receive
    nzstatus (MPI_STATUS_SIZE), & ! for MPI-WAIT
    ncount, type_handle,        & ! return values from setup_data_type
    MPI_neighbors(4), i,        & ! same as neighbors, if neighbor exists
    ilocalreq(4)                  ! the local requests for the ISEND and IRECV

  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: sendbuf
  INTEGER  :: isendbuflen
  INTEGER  :: izmplcode                   ! for MPI error code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ierror     = 0
  izmplcode  = 0
  yerrmsg    = '    '

  isendbuflen = MAX(idim,jdim)* kdim*trac_dim*nlines
  ALLOCATE (sendbuf(isendbuflen, 8))
 ! check whether nlines <= nboundlines
  IF (nlines > nboundlines) THEN
    ierror  = 9011
    yerrmsg     = ' *** nlines > nboundlines *** '
    RETURN
  ENDIF

  ! Determine the start- and end-indices for routines putbuf, getbuf
  izlo_ls = nboundlines + 1
  izup_ls = nboundlines + nlines
  jzlo_ls = jstartpar
  jzup_ls = jendpar

  izlo_lr = nboundlines + 1 - nlines
  izup_lr = nboundlines
  jzlo_lr = jstartpar
  jzup_lr = jendpar

  izlo_us = 1
  izup_us = idim
  jzlo_us = jdim - nboundlines - nlines + 1
  jzup_us = jdim - nboundlines

  izlo_ur = 1
  izup_ur = idim
  jzlo_ur = jdim - nboundlines + 1
  jzup_ur = jdim - nboundlines + nlines

  izlo_rs = idim - nboundlines - nlines + 1
  izup_rs = idim - nboundlines
  jzlo_rs = jstartpar
  jzup_rs = jendpar

  izlo_rr = idim - nboundlines + 1
  izup_rr = idim - nboundlines + nlines
  jzlo_rr = jstartpar
  jzup_rr = jendpar

  izlo_ds = 1
  izup_ds = idim
  jzlo_ds = nboundlines + 1
  jzup_ds = nboundlines + nlines

  izlo_dr = 1
  izup_dr = idim
  jzlo_dr = nboundlines + 1 - nlines
  jzup_dr = nboundlines

  nzcount_lr = 0
  nzcount_rr = 0
  nzcount_ur = 0
  nzcount_dr = 0
  nzcount_ls = 0
  nzcount_rs = 0
  nzcount_us = 0
  nzcount_ds = 0

  ! Fix list of neighbors (use MPI_PROC_NULL rather than -1 to indicate
  ! missing neighbor).
  DO i= 1, 4
    IF ( neighbors(i) /= -1 ) THEN
      MPI_neighbors(i) = neighbors(i)
    ELSE
      MPI_neighbors(i) = MPI_PROC_NULL
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- Section 2: Determine the necessary datatypes
!------------------------------------------------------------------------------

  IF (lmpi_types) THEN
    ! Exchange with left and right neighbor
    IF ( iexchg_MPI_types(2*icase-1) == MPI_DATATYPE_NULL ) THEN
      IF ( MPI_neighbors(1) /= MPI_PROC_NULL ) THEN
        CALL setup_data_type( trac_field, trac_dim,                        &
           idim, jdim, kdim, izlo_ls, izup_ls, jzlo_ls, jzup_ls,           &
           ierror, yerrmsg, imp_type, ncount, type_handle )
        iexchg_MPI_types(2*icase-1) = type_handle
        iexchg_counts   (2*icase-1) = ncount
      ELSE
        CALL setup_data_type(trac_field, trac_dim,                         &
           idim, jdim, kdim, izlo_rr, izup_rr, jzlo_rr, jzup_rr,           &
           ierror, yerrmsg, imp_type, ncount, type_handle )
        iexchg_MPI_types(2*icase-1) = type_handle
        iexchg_counts   (2*icase-1) = ncount
      ENDIF
    ENDIF
  
    ! Exchange with upper and lower neighbor
    IF ( iexchg_MPI_types(2*icase) == MPI_DATATYPE_NULL ) THEN
      IF ( MPI_neighbors(2) /= MPI_PROC_NULL ) THEN
        CALL setup_data_type(trac_field, trac_dim,                         &
           idim, jdim, kdim, izlo_us, izup_us, jzlo_us, jzup_us,           &
           ierror, yerrmsg, imp_type, ncount, type_handle )
        iexchg_MPI_types(2*icase) = type_handle
        iexchg_counts   (2*icase) = ncount
      ELSE
        CALL setup_data_type(trac_field, trac_dim,                         &
           idim, jdim, kdim, izlo_dr, izup_dr, jzlo_dr, jzup_dr,           &
           ierror, yerrmsg, imp_type, ncount, type_handle )
        iexchg_MPI_types(2*icase) = type_handle
        iexchg_counts   (2*icase) = ncount
      ENDIF
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- Section 3: Exchange with immediate Send and blocking Recv
!------------------------------------------------------------------------------

  IF   (ntype == 1) THEN

    IF (lmpi_types) THEN

      !------------------------------------------------------------------------
      !- Section 3.1: exchange with left and right neighbors using datatypes
      !------------------------------------------------------------------------

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
        CALL MPI_ISEND ( trac_field(izlo_ls,jzlo_ls,1,1) &
                       , iexchg_counts(2*icase-1), &
                         iexchg_MPI_types(2*icase-1), MPI_neighbors(1),      &
                         ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! right neighbor is present
        ! send the data
        CALL MPI_ISEND ( trac_field(izlo_rs,jzlo_rs,1,1) &
                       , iexchg_counts(2*icase-1), &
                         iexchg_MPI_types(2*icase-1), MPI_neighbors(3),      &
                         ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! right neighbor is present
        ! receive the data
        CALL MPI_RECV (trac_field(izlo_rr,jzlo_rr,1,1) &
                     , iexchg_counts(2*icase-1),   &
                       iexchg_MPI_types(2*icase-1), MPI_neighbors(3),        &
                       ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
        ! receive the data
        CALL MPI_RECV ( trac_field(izlo_lr,jzlo_lr,1,1) &
                      , iexchg_counts(2*icase-1),  &
                        iexchg_MPI_types(2*icase-1), MPI_neighbors(1),       &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! wait for the completion of the last send to neighbors(1)
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! wait for the completion of the last send to neighbors(3)
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      !------------------------------------------------------------------------
      !- Section 3.2: exchange with upper and lower neighbors using datatypes
      !------------------------------------------------------------------------

      IF (neighbors(2) /= -1) THEN
        ! upper neighbor is present
        ! send the data
        CALL MPI_ISEND (trac_field(izlo_us,jzlo_us,1,1) &
                      , iexchg_counts(2*icase),    &
                        iexchg_MPI_types(2*icase), MPI_neighbors(2),         &
                        ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! send the data
        CALL MPI_ISEND (trac_field(izlo_ds,jzlo_ds,1,1), iexchg_counts(2*icase),    &
                        iexchg_MPI_types(2*icase), MPI_neighbors(4),         &
                        ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! lower neighbor is present
        CALL MPI_RECV (trac_field(izlo_dr,jzlo_dr,1,1), iexchg_counts(2*icase),    &
                       iexchg_MPI_types(2*icase), MPI_neighbors(4),         &
                       ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        ! upper neighbor is present
        CALL MPI_RECV (trac_field(izlo_ur,jzlo_ur,1,1), iexchg_counts(2*icase),     &
                        iexchg_MPI_types(2*icase), MPI_neighbors(2),         &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF
      ENDIF
  
      IF (neighbors(2) /= -1) THEN
        ! wait for the completion of the last send to neighbors(2)
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF
  
      IF (neighbors(4) /= -1) THEN
        ! wait for the completion of the last send to neighbors(4)
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

    ELSE

      !------------------------------------------------------------------------
      !- Section 3.3: exchange with left and right neigh. using explict buff.
      !------------------------------------------------------------------------

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
  
        ! determine start- and end-indices for routine putbuf
        nzcount_ls = 0
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),   &
                         ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! right neighbor is present

        ! determine start- and end-indices for routine putbuf
        nzcount_rs = 0
        CALL putbuf ( trac_field, trac_dim,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),   &
                         ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! receive the data
        CALL MPI_RECV ( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF

        CALL getbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
        ! receive the data

        CALL MPI_RECV ( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF

        CALL getbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! wait for the completion of the last send to neighbors(1)
        ! to safely reuse sendbuf(1,1)
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! wait for the completion of the last send to neighbors(3)
        ! to safely reuse sendbuf(1,3)
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
      ENDIF

      !------------------------------------------------------------------------
      !- Section 3.4: exchange with lower and upper neigh. using explict buff.
      !------------------------------------------------------------------------

      IF (neighbors(2) /= -1) THEN
        nzcount_us = 0
        CALL putbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),   &
                         ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! lower neighbor is present

        nzcount_ds = 0
        CALL putbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),   &
                         ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! lower neighbor is present
        ! receive the data
  
        CALL MPI_RECV ( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF
  
        CALL getbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF
  
      IF (neighbors(2) /= -1) THEN
        ! upper neighbor is present
        ! receive the data
  
        CALL MPI_RECV ( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF
  
        CALL getbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF
  
      IF (neighbors(2) /= -1) THEN
        ! wait for the completion of the last send to neighbors(2)
        ! to safely reuse sendbuf(1,2)
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF
  
      IF (neighbors(4) /= -1) THEN
        ! wait for the completion of the last send to neighbors(4)
        ! to safely reuse sendbuf(1,4)
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

    ENDIF

!------------------------------------------------------------------------------
!- Section 4: Exchange with immediate Recv and blocking Send
!------------------------------------------------------------------------------

  ELSEIF (ntype == 2) THEN

    IF (lmpi_types) THEN

      !------------------------------------------------------------------------
      !- Section 4.1: exchange with left and right neighbors
      !------------------------------------------------------------------------

      IF (neighbors(3) /= -1) THEN
        CALL MPI_IRECV (trac_field(izlo_rr,jzlo_rr,1,1) &
                     , iexchg_counts(2*icase-1), &
                        iexchg_MPI_types(2*icase-1), MPI_neighbors(3),      &
                        ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_IRECV (trac_field(izlo_lr,jzlo_lr,1,1), iexchg_counts(2*icase-1),&
                        iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                        ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! now send the data to the left neighbor
        CALL MPI_SEND ( trac_field(izlo_ls,jzlo_ls,1,1), iexchg_counts(2*icase-1),&
                        iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! now send the data to the right neighbor
        CALL MPI_SEND ( trac_field(izlo_rs,jzlo_rs,1,1), iexchg_counts(2*icase-1),&
                        iexchg_MPI_types(2*icase-1), MPI_neighbors(3),     &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      !------------------------------------------------------------------------
      !- Section 4.2: exchange with upper and lower neighbors
      !------------------------------------------------------------------------

      IF (neighbors(2) /= -1) THEN
        CALL MPI_IRECV (trac_field(izlo_ur,jzlo_ur,1,1), iexchg_counts(2*icase),  &
                        iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                        ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        CALL MPI_IRECV (trac_field(izlo_dr,jzlo_dr,1,1), iexchg_counts(2*icase),  &
                        iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                        ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        CALL MPI_SEND (trac_field(izlo_ds,jzlo_ds,1,1), iexchg_counts(2*icase),  &
                       iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                       ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        CALL MPI_SEND (trac_field(izlo_us,jzlo_us,1,1), iexchg_counts(2*icase),  &
                       iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                       ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

    ELSE

      !------------------------------------------------------------------------
      !- Section 4.3: exchange with left and right neighbors
      !------------------------------------------------------------------------

      IF (neighbors(3) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),  &
                        ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),  &
                        ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        nzcount_ls = 0
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )

        CALL MPI_SEND ( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        nzcount_rs = 0
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )

        CALL MPI_SEND ( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF

        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF

        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      !------------------------------------------------------------------------
      !- Section 4.4: exchange with upper and lower neighbors
      !------------------------------------------------------------------------

      IF (neighbors(2) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),  &
                        ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),  &
                        ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        nzcount_ds = 0
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )

        CALL MPI_SEND ( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        nzcount_us = 0
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )
 
        CALL MPI_SEND ( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF

        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF

 
      IF (neighbors(4) /= -1) THEN
        ! Now wait until the data have arrived
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF

        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF

    ENDIF

!------------------------------------------------------------------------------
!- Section 5: Exchange with SendRecv
!------------------------------------------------------------------------------

  ELSEIF (ntype == 3) THEN

    IF (lmpi_types) THEN

      !--------------------------------------------------------------------------
      !- Section 5.1: Send data to the left and receive from the right neighbor
      !--------------------------------------------------------------------------
 
      CALL MPI_SENDRECV ( trac_field(izlo_ls,jzlo_ls,1,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                          ntag,                                              &
                          trac_field(izlo_rr,jzlo_rr,1,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(3),     &
                          ntag,                                              &
                          icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV-1'
        RETURN
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.2: Receive data from the left and send to the right neighbor
      !--------------------------------------------------------------------------

      CALL MPI_SENDRECV ( trac_field(izlo_rs,jzlo_rs,1,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(3),     &
                          ntag+1,                                            &
                          trac_field(izlo_lr,jzlo_lr,1,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                          ntag+1,                                            &
                          icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV-2'
        RETURN
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.3: Send data to the upper and receive from the lower neighbor
      !--------------------------------------------------------------------------

      CALL MPI_SENDRECV ( trac_field(izlo_us,jzlo_us,1,1), iexchg_counts(2*icase),  &
                          iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                          ntag+2,                                            &
                          trac_field(izlo_dr,jzlo_dr,1,1), iexchg_counts(2*icase),  &
                          iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                          ntag+2,                                            &
                          icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV-3'
        RETURN
      ENDIF
 

      !--------------------------------------------------------------------------
      !- Section 5.4: Receive data from the upper and send to the lower neighbor
      !--------------------------------------------------------------------------
 
      CALL MPI_SENDRECV ( trac_field(izlo_ds,jzlo_ds,1,1), iexchg_counts(2*icase),  &
                          iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                          ntag+3,                                            &
                          trac_field(izlo_ur,jzlo_ur,1,1), iexchg_counts(2*icase),  &
                          iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                          ntag+3,                                            &
                          icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV-4'
        RETURN
      ENDIF

    ELSE

      !--------------------------------------------------------------------------
      !- Section 5.5: Send data to the left and receive from the right neighbor
      !--------------------------------------------------------------------------

      nzcount_ls = 0
      IF (MPI_neighbors(1) /= MPI_PROC_NULL) THEN
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,1), nzcount_ls,  imp_type, MPI_neighbors(1), ntag,    &
             sendbuf(1,7), isendbuflen, imp_type, MPI_neighbors(3), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV'
        RETURN
      ENDIF

      IF (MPI_neighbors(3) /= MPI_PROC_NULL) THEN
        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.6: Send data to the right and receive from the left neighbor
      !--------------------------------------------------------------------------

      nzcount_rs = 0
      IF (MPI_neighbors(3) /= MPI_PROC_NULL) THEN
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,3), nzcount_rs,  imp_type, MPI_neighbors(3), ntag,    &
             sendbuf(1,5), isendbuflen, imp_type, MPI_neighbors(1), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV'
        RETURN
      ENDIF

      IF (MPI_neighbors(1) /= MPI_PROC_NULL) THEN
        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.7: Send data to the upper and receive from the lower neighbor
      !--------------------------------------------------------------------------
 
      nzcount_us = 0
      IF (MPI_neighbors(2) /= MPI_PROC_NULL) THEN
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,2), nzcount_us,  imp_type, MPI_neighbors(2), ntag,    &
             sendbuf(1,8), isendbuflen, imp_type, MPI_neighbors(4), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV'
        RETURN
      ENDIF

      IF (MPI_neighbors(4) /= MPI_PROC_NULL) THEN
        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.8: Send data to the lower and receive from the upper neighbor
      !--------------------------------------------------------------------------
 
      nzcount_ds = 0
      IF (MPI_neighbors(4) /= MPI_PROC_NULL) THEN
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,4), nzcount_ds,  imp_type, MPI_neighbors(4), ntag,    &
             sendbuf(1,6), isendbuflen, imp_type, MPI_neighbors(2), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV'
        RETURN
      ENDIF

      IF (MPI_neighbors(2) /= MPI_PROC_NULL) THEN
        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF

    ENDIF

  ENDIF

  DEALLOCATE (sendbuf)
!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE exchg_trac_boundaries

!==============================================================================
!==============================================================================
!+ This subroutine puts all necessary values into sendbuf
!------------------------------------------------------------------------------

SUBROUTINE putbuf  (trac_field, trac_dim,                    &
                    sendbuf, isendbuflen, idim, jdim, kdim,  &
                    ilo, iup, jlo, jup, ncount, nentry )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine puts the necessary values from the present variables
!   (determined by ilo, iup, jlo, jup) into sendbuf.
!
! Method:
!   Check which variables are present.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    isendbuflen,                  & ! length of sendbuffer
    idim, jdim, kdim,         & ! dimensions of the fields
    trac_dim, &
    ilo, iup, jlo, jup,           & ! start- and end-indices
    nentry                          ! specifies the row of the sendbuf

  INTEGER, INTENT (INOUT)      ::    &
    ncount                          ! counts the variables

  REAL (dp), INTENT (INOUT)            ::    &
    sendbuf (isendbuflen, 8),     & ! send buffer
    trac_field (idim, jdim, trac_dim, kdim)   ! first field that has to occur
! Local variables

  INTEGER :: i, j, jt, k, nzc

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 2: Put data into the buffer
!------------------------------------------------------------------------------

  ! use nzc as a local counter  (based on a work from Mike O'Neill to 
  ! improve vectorization of putbuf and getbuf)
  nzc = ncount

  ! first variable that has to be present
  DO k = 1, kdim
     DO jt =1, trac_dim
        DO j = jlo, jup
           DO i = ilo, iup
              nzc = nzc + 1
              sendbuf (nzc,nentry) = trac_field(i,j,jt,k)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! put nzc to global counter
  ncount = nzc

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  
END SUBROUTINE putbuf 

!==============================================================================
!==============================================================================
!+ This subroutine gets all necessary values from sendbuf
!------------------------------------------------------------------------------

SUBROUTINE getbuf  (trac_field, trac_dim, &
                    sendbuf, isendbuflen, idim, jdim, kdim,                 &
                    ilo, iup, jlo, jup, ncount, nentry )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine gets the necessary values for the present variables
!   (determined by ilo, iup, jlo, jup) from sendbuf.
!
! Method:
!   Check which variables are present.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    isendbuflen,                  & ! length of sendbuffer
    idim, jdim, kdim,         & ! dimensions of the fields
    trac_dim,     &
    ilo, iup, jlo, jup,           & ! start- and end-indices
    nentry                          ! specifies the row of sendbuf to be used

  INTEGER, INTENT (INOUT)  :: ncount      ! counts the variables

  REAL (dp), INTENT (INOUT)            ::    &
    sendbuf (isendbuflen, 8),     & ! send buffer
    trac_field(idim, jdim, trac_dim,kdim)   ! first field that has to occur

! Local variables

  INTEGER ::   i, j, jt, k, nzc

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!- Section 2: Get data from the buffer
!------------------------------------------------------------------------------

  ! use nzc as a local counter  (based on a work from Mike O'Neill to 
  ! improve vectorization of putbuf and getbuf)
  nzc = ncount

  ! first variable that has to be present
  DO k = 1, kdim
     DO jt = 1, trac_dim
        DO j = jlo, jup
           DO i = ilo, iup
              nzc = nzc + 1
              trac_Field(i,j,jt,k) = sendbuf (nzc,nentry)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
    
  ! put nzc to global counter
  ncount = nzc

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  
END SUBROUTINE getbuf 

!==============================================================================
!==============================================================================
!+ Calls MPI_BARRIER
!------------------------------------------------------------------------------

SUBROUTINE comm_barrier (icomm, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!   MPI-routine MPI_BARRIER
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)  :: icomm            ! communicator to be used

  INTEGER, INTENT (OUT) :: ierror           ! error-status variable

  CHARACTER (LEN=*), INTENT (OUT) :: yerrmsg  ! for MPI error message

! Local variables

  INTEGER ::  izmplcode                     ! for MPI error code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

ierror    = 0
yerrmsg   = '   '
izmplcode = 0

CALL MPI_BARRIER (icomm, izmplcode)

IF (izmplcode /= 0) THEN
  ierror = izmplcode
  yerrmsg = 'MPI_BARRIER'
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE comm_barrier

!==============================================================================
!==============================================================================
!+ This subroutine defines and allocates MPI data types for exchg_boundaries
!------------------------------------------------------------------------------

SUBROUTINE setup_data_type(trac_field, trac_dim,                             &
       idim, jdim, kdim, ilo, iup, jlo, jup,                                 &
       ierror, yerrmsg, imp_type, ncount, type_handle)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine allocates and commits a data type consisting of a
!   certain subgrid of up to 14 variables. The subgrid is specified
!   via  idim, jdim, kdim, ilo, iup, jlo, jup, klo, kup . Only
!   one variable has to occur, the other are optional.
!
!   The values of "var1(ilo, jlo, klo), ncount, datatype" 
!   should constitute the right entries for the dummy variables 
!   "BUF, COUNT, DATATYPE" in a call to
!   MPI_SEND/MPI_RECV. This should describe the subgrid for all
!   up to 14 variables.
!
!   As a consequence, the same data type can only be used, if the same
!   variables (same position in memory !) are used. 
!
!  Author: C. Pospiech, IBM
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    imp_type,                     & ! determines the REAL type used
    idim, jdim,                   & ! horizontal dimensions of the fields
    kdim,                     & ! vertical dimensions of var01..var20
    trac_dim, &
    ilo, iup, jlo, jup              ! start- and end-indices

  INTEGER, INTENT (OUT)      ::    &
    ncount,                       & ! how many copies of type_handle
    type_handle                     ! handle for MPI data type

  INTEGER,           INTENT (OUT) ::  ierror       ! error status variable

  CHARACTER (LEN=*), INTENT(OUT)  :: yerrmsg       ! for MPI error message

  REAL(dp), INTENT (INOUT) ::  trac_field(idim,jdim,trac_dim, kdim) 

! Local variables

  INTEGER ::               &
    nzc,                   &
    sect1d, sect2d,        &   ! Variables to hold intermediate 
    sect3d, sect4d,        &
    meta_vect,             &   ! MPI data types
    meta_addr,             &   ! Vector of addresses of the varxx
    meta_disp,             &   ! Displacements of varxx in memory
    meta_blklen,           &   ! some intermediate variable
    num_meta_entries,      &   ! how many varxx are present
    disp(2), blocklen(2),  &   ! Variables needed to define
    vartype(2),            &   ! MPI data type of a certain extent
    sizeofreal                 ! size of data type in byte

  INTEGER :: izmplcode                   ! for MPI error code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  sect3d = MPI_DATATYPE_NULL
  sect4d = MPI_DATATYPE_NULL

!------------------------------------------------------------------------------
!- Section 2: Set up of MPI data types *** subarrays
!------------------------------------------------------------------------------

  ! set up 1-dimensional section
  nzc = iup - ilo + 1
  CALL MPI_TYPE_CONTIGUOUS  (nzc, imp_type, sect1d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_CONTIGUOUS TRAC'
    RETURN
  ENDIF

  ! set up 2-dimensional section
  nzc = jup - jlo +1
  CALL MPI_TYPE_EXTENT(imp_type, sizeofreal, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_EXTENT TRAC'
    RETURN
  ENDIF
  CALL MPI_TYPE_HVECTOR    (nzc, 1, idim*sizeofreal,         &
                            sect1d, sect2d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-2 TRAC'
    RETURN
  ENDIF 

!US: this must be done for every optional entry
  ! set up 3-dimensional section
  CALL MPI_TYPE_HVECTOR    (trac_dim, 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-3 TRAC'
    RETURN
  ENDIF 
  CALL MPI_TYPE_HVECTOR    (kdim, 1, idim*jdim*sizeofreal,    &
                              sect3d, type_handle, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-4 TRAC'
    RETURN
  ENDIF 

!------------------------------------------------------------------------------
!- Section 3: Set up of MPI data types *** meta structure from all varxx
!------------------------------------------------------------------------------

  meta_addr = 0 ! initialize 

  CALL MPI_ADDRESS (trac_field, meta_addr, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_ADDRESS-01 TRAC'
    RETURN
  ENDIF

!!$  meta_disp(:)   = meta_addr(:) - meta_addr(1)
!!$  meta_blklen(:) = 1
!!$  CALL MPI_TYPE_STRUCT  (num_meta_entries, meta_blklen, meta_disp, sect3d,  &
!!$                         meta_vect, izmplcode)
!!$  IF (izmplcode /= 0) THEN
!!$    ierror  = izmplcode
!!$    yerrmsg = 'MPI_TYPE_STRUCT'
!!$    RETURN
!!$  ENDIF
 
!------------------------------------------------------------------------------
!- Section 4: Reset extent of this new data type by defining upper bound
!------------------------------------------------------------------------------

!US Bin mir nicht sicher, wozu das gut sein soll
!!$  blocklen(:)   = 1
!!$  disp(1)     = 0
!!$  disp(2)     = sizeofreal*(iup - ilo + 2)
!!$  vartype(1)  = meta_vect
!!$  vartype(2)  = MPI_UB
!!$  CALL MPI_TYPE_STRUCT     (2, blocklen, disp, vartype, type_handle, izmplcode)
!!$  IF (izmplcode /= 0) THEN
!!$    ierror  = izmplcode
!!$    yerrmsg = 'MPI_TYPE_STRUCT'
!!$    RETURN
!!$  ENDIF

!------------------------------------------------------------------------------
!- Section 5: Commit the data type
!------------------------------------------------------------------------------

  CALL MPI_TYPE_COMMIT       (type_handle,izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_COMMIT'
    RETURN
  ENDIF
  ncount = 1

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE setup_data_type

!==============================================================================
#endif
! ifndef BLANK

#endif
! if (defined COSMO) || defined(BLANK)

! op_bk_20130820+
#if defined (__ICON__)
  USE mo_mpi,           ONLY: my_process_is_stdio              &  ! LOGICAL FUNC
       &                    , p_io=>process_mpi_stdio_id       &  ! INTEGER
       &                    , get_my_mpi_all_id                &  ! p_pe
       &                    , get_my_mpi_all_comm_size         &  ! p_nprocs
       &                    , p_bcast                          &  ! SUBROUTINES
       &                    , p_send, p_recv, p_sendrecv       &  ! SUBROUTINES
       &                    , p_all_comm=>process_mpi_all_comm &  ! INTEGER
       &                    , p_abort                          &
       &                    , my_process_is_mpi_all_parallel   &  ! p_parallel
       &                    , p_sum                            &  ! SUBROUTINE
       &                    , get_my_mpi_all_id                &  ! p_pe
       &                    , my_process_is_mpi_all_parallel   &  ! p_parallel
       &                    , my_process_is_stdio                 ! p_parallel_io
  ! op_bk_20130828-



  USE mo_exception,     ONLY: finish, message                ! SUBROUTINES

  USE messy_main_constants_mem, ONLY: dp, sp, i4, i8

  IMPLICIT NONE

#ifndef NOMPI
  INCLUDE 'mpif.h'
#endif

  INTEGER, PUBLIC, SAVE :: p_pe, p_nprocs
  LOGICAL, PUBLIC, SAVE :: p_parallel, p_parallel_io

! NOTE: The following is used as a dummy in order to mimic the
!       vectorisation procedures of ECHAM5 for COSMO,
!       e.g., bi_vector, bi_decompose etc.

  TYPE decomp
     LOGICAL :: lreg = .TRUE.
  END TYPE decomp
  !
  TYPE(decomp), PUBLIC, POINTER :: dcg
  TYPE(decomp), PUBLIC, SAVE    :: dcl 

  INTEGER, PARAMETER :: iexchg_MPI_type_len    = 200
       ! length of global vector iexchg_MPI_type used in environment.f90


  INTEGER :: iexchg_MPI_types(iexchg_MPI_type_len) = MPI_DATATYPE_NULL
                  ! List of MPI data types used in exchg_datatypes
                  ! routine. If set to MPI_DATATYPE_NULL, the
                  ! corresponding entry has not been properly set up.
  INTEGER ::  iexchg_counts(iexchg_MPI_type_len)
                  ! List of counts of these data types
                  ! used in exchg_datatypes routine.
                  ! Has meaningful value only if the corresponding
                  ! vector element of iexchg_MPI_types is not
                  ! MPI_DATATYPE_NULL.

  ! INTERFACE
  INTERFACE reorder
    MODULE PROCEDURE reorder2
    MODULE PROCEDURE reorder3
    MODULE PROCEDURE reorder4
  END INTERFACE

  PUBLIC :: p_io
  PUBLIC :: p_abort
  PUBLIC :: p_bcast
  PUBLIC :: p_all_comm

! SUBROUTINES
  PUBLIC reorder
  PUBLIC :: messy_mpi_initialize
! PUBLIC :: p_abort
! PUBLIC :: bcast_4d

CONTAINS

SUBROUTINE messy_mpi_initialize
!
! This subroutines initializes some variables needed in MESSy
!
  p_parallel = my_process_is_mpi_all_parallel()
  p_parallel_io = my_process_is_stdio()
  p_pe = get_my_mpi_all_id()
  p_nprocs = get_my_mpi_all_comm_size()

END SUBROUTINE messy_mpi_initialize

!==============================================================================
!==============================================================================
!+ This subroutine puts all necessary values into sendbuf
!------------------------------------------------------------------------------

SUBROUTINE putbuf  (trac_field, trac_dim,                    &
                    sendbuf, isendbuflen, idim, jdim, kdim,  &
                    ilo, iup, jlo, jup, ncount, nentry )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine puts the necessary values from the present variables
!   (determined by ilo, iup, jlo, jup) into sendbuf.
!
! Method:
!   Check which variables are present.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    isendbuflen,                  & ! length of sendbuffer
    idim, jdim, kdim,         & ! dimensions of the fields
    trac_dim, &
    ilo, iup, jlo, jup,           & ! start- and end-indices
    nentry                          ! specifies the row of the sendbuf

  INTEGER, INTENT (INOUT)      ::    &
    ncount                          ! counts the variables

  REAL (dp), INTENT (INOUT)            ::    &
    sendbuf (isendbuflen, 8),     & ! send buffer
    trac_field (idim, jdim, trac_dim, kdim)   ! first field that has to occur
! Local variables

  INTEGER :: i, j, jt, k, nzc

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 2: Put data into the buffer
!------------------------------------------------------------------------------

  ! use nzc as a local counter  (based on a work from Mike O'Neill to 
  ! improve vectorization of putbuf and getbuf)
  nzc = ncount

  ! first variable that has to be present
  DO k = 1, kdim
     DO jt =1, trac_dim
        DO j = jlo, jup
           DO i = ilo, iup
              nzc = nzc + 1
              sendbuf (nzc,nentry) = trac_field(i,j,jt,k)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! put nzc to global counter
  ncount = nzc

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  
END SUBROUTINE putbuf 

!==============================================================================
!==============================================================================
!+ This subroutine gets all necessary values from sendbuf
!------------------------------------------------------------------------------

SUBROUTINE getbuf  (trac_field, trac_dim, &
                    sendbuf, isendbuflen, idim, jdim, kdim,                 &
                    ilo, iup, jlo, jup, ncount, nentry )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine gets the necessary values for the present variables
!   (determined by ilo, iup, jlo, jup) from sendbuf.
!
! Method:
!   Check which variables are present.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    isendbuflen,                  & ! length of sendbuffer
    idim, jdim, kdim,         & ! dimensions of the fields
    trac_dim,     &
    ilo, iup, jlo, jup,           & ! start- and end-indices
    nentry                          ! specifies the row of sendbuf to be used

  INTEGER, INTENT (INOUT)  :: ncount      ! counts the variables

  REAL (dp), INTENT (INOUT)            ::    &
    sendbuf (isendbuflen, 8),     & ! send buffer
    trac_field(idim, jdim, trac_dim,kdim)   ! first field that has to occur

! Local variables

  INTEGER ::   i, j, jt, k, nzc

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!- Section 2: Get data from the buffer
!------------------------------------------------------------------------------

  ! use nzc as a local counter  (based on a work from Mike O'Neill to 
  ! improve vectorization of putbuf and getbuf)
  nzc = ncount

  ! first variable that has to be present
  DO k = 1, kdim
     DO jt = 1, trac_dim
        DO j = jlo, jup
           DO i = ilo, iup
              nzc = nzc + 1
              trac_Field(i,j,jt,k) = sendbuf (nzc,nentry)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
    
  ! put nzc to global counter
  ncount = nzc

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  
END SUBROUTINE getbuf 

!==============================================================================
!==============================================================================
!+ Calls MPI_BARRIER
!------------------------------------------------------------------------------

SUBROUTINE comm_barrier (icomm, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!   MPI-routine MPI_BARRIER
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)  :: icomm            ! communicator to be used

  INTEGER, INTENT (OUT) :: ierror           ! error-status variable

  CHARACTER (LEN=*), INTENT (OUT) :: yerrmsg  ! for MPI error message

! Local variables

  INTEGER ::  izmplcode                     ! for MPI error code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

ierror    = 0
yerrmsg   = '   '
izmplcode = 0

CALL MPI_BARRIER (icomm, izmplcode)

IF (izmplcode /= 0) THEN
  ierror = izmplcode
  yerrmsg = 'MPI_BARRIER'
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE comm_barrier

!==============================================================================
!==============================================================================
!+ This subroutine defines and allocates MPI data types for exchg_boundaries
!------------------------------------------------------------------------------

SUBROUTINE setup_data_type(trac_field, trac_dim,                             &
       idim, jdim, kdim, ilo, iup, jlo, jup,                                 &
       ierror, yerrmsg, imp_type, ncount, type_handle)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine allocates and commits a data type consisting of a
!   certain subgrid of up to 14 variables. The subgrid is specified
!   via  idim, jdim, kdim, ilo, iup, jlo, jup, klo, kup . Only
!   one variable has to occur, the other are optional.
!
!   The values of "var1(ilo, jlo, klo), ncount, datatype" 
!   should constitute the right entries for the dummy variables 
!   "BUF, COUNT, DATATYPE" in a call to
!   MPI_SEND/MPI_RECV. This should describe the subgrid for all
!   up to 14 variables.
!
!   As a consequence, the same data type can only be used, if the same
!   variables (same position in memory !) are used. 
!
!  Author: C. Pospiech, IBM
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    imp_type,                     & ! determines the REAL type used
    idim, jdim,                   & ! horizontal dimensions of the fields
    kdim,                     & ! vertical dimensions of var01..var20
    trac_dim, &
    ilo, iup, jlo, jup              ! start- and end-indices

  INTEGER, INTENT (OUT)      ::    &
    ncount,                       & ! how many copies of type_handle
    type_handle                     ! handle for MPI data type

  INTEGER,           INTENT (OUT) ::  ierror       ! error status variable

  CHARACTER (LEN=*), INTENT(OUT)  :: yerrmsg       ! for MPI error message

  REAL(dp), INTENT (INOUT) ::  trac_field(idim,jdim,trac_dim, kdim) 

! Local variables

  INTEGER ::               &
    nzc,                   &
    sect1d, sect2d,        &   ! Variables to hold intermediate 
    sect3d, sect4d,        &
    meta_vect,             &   ! MPI data types
    meta_addr,             &   ! Vector of addresses of the varxx
    meta_disp,             &   ! Displacements of varxx in memory
    meta_blklen,           &   ! some intermediate variable
    num_meta_entries,      &   ! how many varxx are present
    disp(2), blocklen(2),  &   ! Variables needed to define
    vartype(2),            &   ! MPI data type of a certain extent
    sizeofreal                 ! size of data type in byte

  INTEGER :: izmplcode                   ! for MPI error code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  sect3d = MPI_DATATYPE_NULL
  sect4d = MPI_DATATYPE_NULL

!------------------------------------------------------------------------------
!- Section 2: Set up of MPI data types *** subarrays
!------------------------------------------------------------------------------

  ! set up 1-dimensional section
  nzc = iup - ilo + 1
  CALL MPI_TYPE_CONTIGUOUS  (nzc, imp_type, sect1d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_CONTIGUOUS TRAC'
    RETURN
  ENDIF

  ! set up 2-dimensional section
  nzc = jup - jlo +1
  CALL MPI_TYPE_EXTENT(imp_type, sizeofreal, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_EXTENT TRAC'
    RETURN
  ENDIF
  CALL MPI_TYPE_HVECTOR    (nzc, 1, idim*sizeofreal,         &
                            sect1d, sect2d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-2 TRAC'
    RETURN
  ENDIF 

!US: this must be done for every optional entry
  ! set up 3-dimensional section
  CALL MPI_TYPE_HVECTOR    (trac_dim, 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-3 TRAC'
    RETURN
  ENDIF 
  CALL MPI_TYPE_HVECTOR    (kdim, 1, idim*jdim*sizeofreal,    &
                              sect3d, type_handle, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-4 TRAC'
    RETURN
  ENDIF 

!------------------------------------------------------------------------------
!- Section 3: Set up of MPI data types *** meta structure from all varxx
!------------------------------------------------------------------------------

  meta_addr = 0 ! initialize 

  CALL MPI_ADDRESS (trac_field, meta_addr, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_ADDRESS-01 TRAC'
    RETURN
  ENDIF

!!$  meta_disp(:)   = meta_addr(:) - meta_addr(1)
!!$  meta_blklen(:) = 1
!!$  CALL MPI_TYPE_STRUCT  (num_meta_entries, meta_blklen, meta_disp, sect3d,  &
!!$                         meta_vect, izmplcode)
!!$  IF (izmplcode /= 0) THEN
!!$    ierror  = izmplcode
!!$    yerrmsg = 'MPI_TYPE_STRUCT'
!!$    RETURN
!!$  ENDIF
 
!------------------------------------------------------------------------------
!- Section 4: Reset extent of this new data type by defining upper bound
!------------------------------------------------------------------------------

!US Bin mir nicht sicher, wozu das gut sein soll
!!$  blocklen(:)   = 1
!!$  disp(1)     = 0
!!$  disp(2)     = sizeofreal*(iup - ilo + 2)
!!$  vartype(1)  = meta_vect
!!$  vartype(2)  = MPI_UB
!!$  CALL MPI_TYPE_STRUCT     (2, blocklen, disp, vartype, type_handle, izmplcode)
!!$  IF (izmplcode /= 0) THEN
!!$    ierror  = izmplcode
!!$    yerrmsg = 'MPI_TYPE_STRUCT'
!!$    RETURN
!!$  ENDIF

!------------------------------------------------------------------------------
!- Section 5: Commit the data type
!------------------------------------------------------------------------------

  CALL MPI_TYPE_COMMIT       (type_handle,izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_COMMIT'
    RETURN
  ENDIF
  ncount = 1

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE setup_data_type

!==============================================================================

!==============================================================================
  SUBROUTINE reorder2 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:)

    y = RESHAPE (x,(/SIZE(y,1),SIZE(y,2)/),(/0._dp/))

  END SUBROUTINE reorder2
!------------------------------------------------------------------------------
  SUBROUTINE reorder3 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:,:)
    INTEGER :: k

    DO k=1,SIZE(x,2)
      y(:,k,:) = RESHAPE (x(:,k,:),(/SIZE(y,1),SIZE(y,3)/),(/0._dp/))
    END DO

  END SUBROUTINE reorder3
!------------------------------------------------------------------------------
  SUBROUTINE reorder4 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:,:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:,:,:)
    INTEGER :: k, l

    DO l=1,SIZE(x,3)
      DO k=1,SIZE(x,2)
        y(:,k,l,:) = RESHAPE (x(:,k,l,:),(/SIZE(y,1),SIZE(y,4)/),(/0._dp/))
      END DO
    END DO

  END SUBROUTINE reorder4
!==============================================================================


#endif
! ifdef __ICON__

! op_bk_20130820-

! **************************************************************************
END MODULE messy_main_mpi_bi
! **************************************************************************
