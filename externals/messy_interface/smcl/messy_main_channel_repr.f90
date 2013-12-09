! **********************************************************************
MODULE messy_main_channel_repr
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Patrick Joeckel, MPICH, May 2005

!#define GENERIC

  USE messy_main_constants_mem,     ONLY: DP, STRLEN_MEDIUM

  USE messy_main_channel_dimensions, ONLY: t_dimension_ptr

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  ! INTERNAL RANK OF ALL REPRESENTATIONS FOR MEMORY MANAGEMENT
  INTEGER, PARAMETER, PUBLIC :: IRANK = 4

  ! UNDEFINED REPRESENTATION
!  INTEGER, PARAMETER, PUBLIC :: UNDEFINED = -2
  INTEGER, PARAMETER, PUBLIC :: REPR_UNDEF = -2

  ! AUTOMATIC CHOICE OF DIMENSION LENGTH
  INTEGER, PARAMETER, PUBLIC :: AUTO = -1   ! LOCAL = GLOBAL

  ! FOR PARALLEL DECOMPOSITION
  INTEGER, PARAMETER, PUBLIC :: PIOTYPE_SGL = 0  ! single PE
  INTEGER, PARAMETER, PUBLIC :: PIOTYPE_IND = 1  ! independent
  INTEGER, PARAMETER, PUBLIC :: PIOTYPE_COL = 2  ! collective

  TYPE t_repr_pdecomp
     ! HOLDS INFORMATION ABOUT PRALLEL DECOMPOSITION
     LOGICAL                  :: lpdecomp = .FALSE.  ! info available ?
     !LOGICAL                 :: l_this_pe = .TRUE.  ! ACTIVE ON THIS PE ?
     !
     ! (UN-DEFORMED (e.g. DE-VECTORISED)) SHAPE IN MEMORY
     INTEGER,DIMENSION(IRANK) :: shape_mem = 0 ! shape in memory (local)
     INTEGER,DIMENSION(IRANK) :: shape_out = 0 ! shape for output (local)
     !
     ! MAPPING BETWEEN GLOBAL AND LOCAL
     INTEGER                          :: nseg = 0        ! number of segments
     INTEGER, DIMENSION(:,:), POINTER :: ml    => NULL() ! lower bounds in mem
     INTEGER, DIMENSION(:,:), POINTER :: mu    => NULL() ! upper bounds in mem
     INTEGER, DIMENSION(:,:), POINTER :: start => NULL() ! start in output
     INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL() ! count in output
     !INTEGER, DIMENSION(:,:), POINTER :: stride => NULL()
     !INTEGER, DIMENSION(:,:), POINTER :: map => NULL()
     ! 
     ! SPECIAL FOR PERFORMANCE TUNING ...
     INTEGER :: piotype = PIOTYPE_IND
     !
  END TYPE t_repr_pdecomp

  ! mz_ab_20090921+
  TYPE t_repr_boundary
     ! HOLDS INFORMATION ABOUT BOUNDARIES
     LOGICAL                   :: lbounds = .FALSE.  ! info available ?
     !
     ! number of boundary indices
     INTEGER, DIMENSION(IRANK) :: nbounds = 0
     !
  END TYPE t_repr_boundary
  ! mz_ab_20090921-

  TYPE t_representation
     ! IDENTIFICATION
     CHARACTER(LEN=STRLEN_MEDIUM) :: name    = ''    ! NAME
     INTEGER                      :: id      = 0     ! ID
     !
     ! op_bk_20130906+
     INTEGER                      :: patch_id = -1   ! patch
     ! op_bk_20130906-
     INTEGER                      :: rank    = 0     ! RANK
     CHARACTER(LEN=IRANK)         :: link    = ''    ! LINK-STRING
     CHARACTER(LEN=IRANK)         :: axis    = ''    ! AXIS-STRING
     INTEGER, DIMENSION(IRANK)    :: gdimlen = 0     ! GLOBAL DIMENSION LENGTH
     !
     ! DECOMPOSITION INFORMATION
     INTEGER, DIMENSION(IRANK)    :: ldimlen = 0     ! LOCAL DIMENSION LENGHT
     !                                               ! (ON THIS PE)
     INTEGER                      :: dctype  = 0     ! DECOMPOSITION TYPE
     !
     ! FULL DIMENSION INFORMATION (POINTER TO DIMENSION; GLOBAL !!!)
     TYPE(t_dimension_ptr), DIMENSION(IRANK) :: dim
     !
     ! INPUT/OUTPUT CONVERSION
     ! - PERMUTATION OF DIMENSIONS BEFORE OUTPUT
     LOGICAL                  :: lperm          = .FALSE.
     !INTEGER,DIMENSION(IRANK) :: shape_mem     = 0 ! shape in memory = gdimlen
     INTEGER,DIMENSION(IRANK) :: order_mem2out  = 0 ! output
     INTEGER,DIMENSION(IRANK) :: shape_out      = 0 ! shape in output
     INTEGER,DIMENSION(IRANK) :: order_out2mem  = 0 ! input
     !
#ifndef GENERIC
     INTEGER :: i_m2o
     INTEGER :: i_o2m
#endif
     !
     ! PARALLEL DECOMPOSITION INFORMATION
     TYPE(t_repr_pdecomp)     :: pdecomp
     !
     ! mz_ab_20090921+
     ! BOUNDARY MEMORY MANAGEMENT
     TYPE(t_repr_boundary)    :: bounds
     ! mz_ab_20090921-
     !
  END TYPE t_representation

  TYPE t_representation_list
     TYPE(t_representation)               :: this
     TYPE(t_representation_list), POINTER :: next => NULL()
  END TYPE t_representation_list

  ! ====================================================================
  ! CONCT. LIST OF REPRESENTATIONS
  TYPE(t_representation_list), POINTER, SAVE         :: GREPRLIST => NULL()
  ! NUMBER OF REPRESENTATIONS
  INTEGER,                              SAVE, PUBLIC :: NREP = 0
  ! ====================================================================

  ! TYPES
  PUBLIC :: t_representation         ! TYPE: representation
  PUBLIC :: t_repr_pdecomp           ! TYPE: parallel decomposition table
  PUBLIC :: t_repr_boundary          ! TYPE: boundaries ! mz_ab_20090921

  ! INTERFACES
  PUBLIC :: new_representation       ! add representation to concat. list
  PUBLIC :: set_representation_decomp ! add decomposition information
  INTERFACE write_representation
     MODULE PROCEDURE write_representation_one
     MODULE PROCEDURE write_representation_all
  END INTERFACE
  PUBLIC :: write_representation     ! write representation (list) to 
  !                                  ! standard output
  PUBLIC :: write_representation_dc  ! write representation decomposition
  INTERFACE get_representation
     MODULE PROCEDURE get_representation_by_name
     MODULE PROCEDURE get_representation_by_id
  END INTERFACE
  PUBLIC :: get_representation       ! get pointer to representation
  PUBLIC :: get_representation_id    ! get representation id from name
  PUBLIC :: clean_representations    ! cleanup memory
  ! um_ak_20090729+
  INTERFACE get_representation_info
     MODULE PROCEDURE get_representation_info
  END INTERFACE
  PUBLIC :: get_representation_info
  ! um_ak_20090729-
  INTERFACE loc_representation
     MODULE PROCEDURE loc_representation_by_name
     MODULE PROCEDURE loc_representation_by_id
  END INTERFACE
  !PRIVATE :: loc_representation
  !
  PUBLIC :: repr_reorder
  PUBLIC :: repr_getptr

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_representation(status, id, name, rank, link, dctype &
       , dimension_ids, ldimlen, output_order, axis, nbounds, patch_id ) 
    ! um_ak_20090805 axis added
    ! mz_ab_20090921 nbounds added
    USE messy_main_channel_dimensions, ONLY: get_dimension
   
    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, PRESENT, SUM, TRIM, ALL

    ! I/O
    INTEGER,                            INTENT(OUT) :: status
    INTEGER,                            INTENT(OUT) :: id
    CHARACTER(LEN=*),                   INTENT(IN)  :: name
    INTEGER,                            INTENT(IN)  :: rank
    CHARACTER(LEN=IRANK),               INTENT(IN)  :: link
    INTEGER,                            INTENT(IN)  :: dctype
    INTEGER,          DIMENSION(rank),  INTENT(IN)  :: dimension_ids
    INTEGER,          DIMENSION(rank),  INTENT(IN)  :: ldimlen
    INTEGER,          DIMENSION(rank),  INTENT(IN), OPTIONAL :: output_order
    ! um_ak_20090805+
    CHARACTER(LEN=IRANK),               INTENT(IN), OPTIONAL :: axis
    ! um_ak_20090805-
    ! mz_ab_20090921+
    INTEGER,          DIMENSION(rank),  INTENT(IN), OPTIONAL :: nbounds
    ! mz_ab_20090921-
    ! op_bk_20130906+
    INTEGER,                            INTENT(IN), OPTIONAL :: patch_id
    ! op_bk_20130906-

    ! LOCAL
    TYPE(t_representation_list), POINTER :: ai   => NULL()
    TYPE(t_representation_list), POINTER :: ae   => NULL()
    TYPE(t_representation),      POINTER :: repr => NULL()
    INTEGER                              :: zstat
    INTEGER                              :: i, j
    CHARACTER(LEN=IRANK)                 :: zlink
    INTEGER                              :: nx
    INTEGER :: ja      ! active rank counter
    INTEGER :: ji      ! inactive rank counter
    INTEGER, PARAMETER :: checksum = 10 ! 1+2+3+4
    INTEGER :: bnd     ! mz_ab_20090921

    ! INIT
    id = 0
    nx = 0
    zlink = ''

    CALL loc_representation_by_name(zstat, GREPRLIST, name, repr)
    IF (zstat /= 2003) THEN  ! REPRESENTATION DOES NOT EXIST (IS OK HERE !)
       IF (zstat == 0) THEN  ! REPRESENTATION EXISTS ALREADY
          status = 2002      ! REPRESENTATION EXISTS ALREADY
       ELSE
          status = zstat    ! ERROR
       END IF
       RETURN
    END IF

    ! CHECK RANK
    IF (rank < 0) THEN
       status = 2016   ! MINIMUM RANK LIMIT EXCEEDED
       RETURN
    END IF
    IF (rank > IRANK) THEN
       status = 2010   ! MAXIMUM RANK LIMIT EXCEEDED
       RETURN
    END IF
    ! 
    ! CHECK 'link'
    DO i=1, IRANK
       SELECT CASE(link(i:i))
          CASE('x','X')
             zlink(i:i) = 'x'
             nx = nx + 1
          CASE('-')
             zlink(i:i) = '-'
          CASE DEFAULT
             status = 2011  ! SYNTAX ERROR IN LINK-STRING
             RETURN
       END SELECT
    END DO
    IF (nx /= rank) THEN
       status = 2012  ! LINK STRING NOT CONFORM WITH RANK
       RETURN
    ENDIF

    ! GOTO END OF LIST
    ai => GREPRLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       ae => ai
       ai => ai%next
    END DO

    ! ADD NEW
   ALLOCATE(ai)
    NULLIFY(ai%next)
    IF (.NOT. ASSOCIATED(GREPRLIST)) THEN
       GREPRLIST => ai                 ! SET POINTER TO FIRST ELEMENT
    ELSE
       ae%next => ai                   ! SET NEXT POINTER OF LAST ELEMENT
       !                               ! TO NEW ELEMENT
    END IF
       
    ! SET VALUES
    ai%this%name   = TRIM(ADJUSTL(name))
    ai%this%rank   = rank
    ai%this%link   = zlink
    ai%this%dctype = dctype
    !
    !
    ! FOR INPUT/OUTPUT CONVERSION
    ! INITIALISE ...
    !ai%this%shape_mem(:)     = 0
    ai%this%order_mem2out(:) = 0
    ai%this%shape_out(:)     = 0
    ai%this%order_out2mem(:) = 0

    ! SET DIMENSION INFORMATION
    ja  = 0
    ji = rank

    ! mz_ab_20090921+
    IF (PRESENT(nbounds)) THEN
       ai%this%bounds%lbounds = .TRUE.
    ELSE
       ai%this%bounds%lbounds = .FALSE.
    END IF
    ! mz_ab_20090921-

    dimension_loop1: DO i=1, IRANK

       IF (zlink(i:i) == '-') THEN
          ! THIS DIMENSION IS INACTIVE
          !
          ji = ji + 1
          !
          ai%this%gdimlen(i) = 1
          ai%this%ldimlen(i) = 1
          ai%this%dim(i)%ptr => NULL()
          !
          ai%this%order_mem2out(ji) = i
          !
       ELSE
          ! THIS DIMENSION IS ACTIVE
          !
          ja = ja + 1
          !
          IF (dimension_ids(ja) <= 0) THEN
             status = 2009 ! DIMENSION ID IS <= ZERO
             RETURN
          END IF
          CALL get_dimension(status, dimension_ids(ja), ai%this%dim(i)%ptr)
          IF (status /= 0) RETURN
          !
          ! SET GLOBAL DIMENSION LENGTH
          ai%this%gdimlen(i) = ai%this%dim(i)%ptr%len
          ! SET LOCAL DIMENSION LENGTH (=GLOBAL OR DECOMPOSED DIMENSION)
          ! mz_ab_20090921+
          IF (ai%this%bounds%lbounds) THEN
             ai%this%bounds%nbounds(ja) = nbounds(ja)
             bnd = 2 * nbounds(ja)
          ELSE
             ai%this%bounds%nbounds(ja) = 0
             bnd = 0
          END IF
          ! mz_ab_20090921-
          IF (ldimlen(ja) == AUTO) THEN
!!$             ai%this%ldimlen(i) = ai%this%gdimlen(i)       ! mz_ab_20090921+
             ai%this%ldimlen(i) = ai%this%gdimlen(i) + bnd    ! mz_ab_20090921+
          ELSE
!!$             ai%this%ldimlen(i) = ldimlen(ja)              ! mz_ab_20090921+
             ai%this%ldimlen(i) = ldimlen(ja) + bnd           ! mz_ab_20090921+
          END IF
          !
          ! SET ORDER OF DIMENSIONS FOR OUTPUT
          IF (PRESENT(output_order)) THEN
             ai%this%order_mem2out(ja) = output_order(ja)
          ELSE
             ai%this%order_mem2out(ja) = i
          END IF
          !
       END IF

    END DO dimension_loop1

    ! SET DIMENSION INFORMATION
    ja  = 0
    ji = rank

    dimension_loop2: DO i=1, IRANK
       ai%this%shape_out(i) = ai%this%gdimlen(ai%this%order_mem2out(i))
       !
       j = ai%this%order_mem2out(i)
       IF (zlink(j:j) == '-') THEN
          ! THIS DIMENSION IS INACTIVE
          ai%this%order_out2mem(j) = i
       ELSE
          ja = ja + 1
          ai%this%order_out2mem(j) = ja
       END IF
    END DO dimension_loop2

    ! CHECKS
    IF ( (SUM(ai%this%order_mem2out) /= checksum) .OR. &
         (SUM(ai%this%order_out2mem) /= checksum) )THEN
       status = 2014 ! INVALID ORDER OF OUTPUT DIMENSIONS
       RETURN
    END IF

    ! PERMUTATION(S) REQUIRED FOR I/O ?
    ai%this%lperm =  .NOT. ( ALL( ai%this%order_mem2out == (/1,2,3,4/) ) &
                       .AND. ALL( ai%this%order_out2mem == (/1,2,3,4/) ) )
    
    ! COUNT AND SET ID
    NREP       = NREP + 1
    ai%this%id = NREP
    id         = NREP

#ifndef GENERIC
    ai%this%i_m2o = 0
    ai%this%i_o2m = 0
    DO i=1, IRANK
       ai%this%i_m2o = ai%this%i_m2o + ai%this%order_mem2out(i)*10**(IRANK-i)
       ai%this%i_o2m = ai%this%i_o2m + ai%this%order_out2mem(i)*10**(IRANK-i)
    END DO
#endif

    ! um_ak_20090729+
    IF (PRESENT(axis)) THEN
        ai%this%axis = axis
    ELSE
       ai%this%axis = '----'
    ENDIF
    ! um_ak_20090729-

    ! op_bk_20130906+
    IF (PRESENT(patch_id)) THEN
       ai%this%patch_id = patch_id
    ENDIF
    ! op_bk_20130906-

    status = 0

  END SUBROUTINE new_representation
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_representation_decomp(status, id, start, cnt, mu, ml &
       , lchk, piotype)

    IMPLICIT NONE

    INTRINSIC :: PRODUCT, MINVAL, MAXVAL, TRIM

    ! I/O
    INTEGER,                   INTENT(OUT) :: status
    INTEGER,                   INTENT(IN)  :: id        ! representation ID
    INTEGER, DIMENSION(:,:),   INTENT(IN)  :: start     ! start vector
    INTEGER, DIMENSION(:,:),   INTENT(IN)  :: cnt       ! count vector
    INTEGER, DIMENSION(:,:),   INTENT(IN)  :: ml        ! memory lower bound
    INTEGER, DIMENSION(:,:),   INTENT(IN)  :: mu        ! memory upper bound
    LOGICAL, OPTIONAL,         INTENT(IN)  :: lchk      ! check memory size
    INTEGER, OPTIONAL,         INTENT(IN)  :: piotype   ! type of parallel IO

    ! LOCAL
    TYPE(t_representation), POINTER  :: repr
    INTEGER                          :: nseg
    INTEGER                          :: nrank
    INTEGER                          :: jr, i
    INTEGER                          :: cs1, cs2, cs3, cs4
    LOGICAL                          :: zlchk

    CALL loc_representation_by_id(status, GREPRLIST, id, repr)
    IF (status /= 0) RETURN

    IF (repr%pdecomp%lpdecomp) THEN
       ! REPRESENTATION DECOMPOSITION TABLE EXISTS ALREADY
       status = 2023
       RETURN
    END IF
    
    ! INIT
    nrank = repr%rank
    IF (PRESENT(lchk)) THEN
       zlchk = lchk
    ELSE
       zlchk = .TRUE. ! DEFAULT
    END IF

    ! NUMBER OF SEGMENTS
    nseg =  SIZE(start,1)
    repr%pdecomp%nseg = nseg

    ! ALLOCATE SPACE
    ALLOCATE(repr%pdecomp%start(nseg,IRANK))
    ALLOCATE(repr%pdecomp%cnt(nseg,IRANK))
    ALLOCATE(repr%pdecomp%ml(nseg,IRANK))
    ALLOCATE(repr%pdecomp%mu(nseg,IRANK))

    ! INITIALIZE
    repr%pdecomp%start(:,:) = 0
    repr%pdecomp%cnt(:,:) = 0
    repr%pdecomp%ml(:,:) = 0
    repr%pdecomp%mu(:,:) = 0

    rank_loop: DO jr=1, nrank
       cs1 = 0
       cs2 = 0
       segment_loop: DO i=1, nseg
          repr%pdecomp%start(i,jr) = start(i,repr%order_mem2out(jr))
          repr%pdecomp%cnt(i,jr)   = cnt(i,repr%order_mem2out(jr))
          cs1 = cs1 + repr%pdecomp%cnt(i,jr)
          repr%pdecomp%ml(i,jr)  = ml(i,repr%order_mem2out(jr))
          repr%pdecomp%mu(i,jr)  = mu(i,repr%order_mem2out(jr))
          IF (repr%pdecomp%cnt(i,jr) /= 0) &
               cs2 = cs2 + ( repr%pdecomp%mu(i,jr) - repr%pdecomp%ml(i,jr) + 1)
       END DO segment_loop
       IF (cs1 /= cs2) THEN
          WRITE(*,*) 'ERROR: (REPRESENTATION ',TRIM(repr%name) &
               ,', RANK ',jr,'): ',cs1,' =/= ',cs2
          status = 2024 ! REPRESENTATION SEGMENTATION MISMATCH
          RETURN
       END IF
    END DO rank_loop

    ! CHECK MEMORY SIZE
    IF (zlchk) THEN
       cs4 = PRODUCT(repr%ldimlen)
       cs3 = 0
       !
       DO i=1, nseg
          cs3 = cs3 + PRODUCT(repr%pdecomp%cnt(i,1:nrank))
       END DO
       !
       IF (cs3 /= cs4) THEN
          WRITE(*,*) 'ERROR: (REPRESENTATION ',TRIM(repr%name) &
               ,'): ',cs3,' =/= ',cs4
          status = 2026 ! REPRESENTATION DECOMPOSITION MISMATCH
          RETURN
       END IF
    END IF

    ! DETERMINE SHAPE IN MEMORY FOR OUTPUT (LOCAL)
    repr%pdecomp%shape_mem(:) = 1
    DO jr=1, IRANK
       repr%pdecomp%shape_mem(jr) = MAXVAL(mu(:,jr))-MINVAL(ml(:,jr))+1
    END DO

    ! DETERMINE SHAPE FOR OUTPUT (LOCAL)
    DO jr=1, IRANK
       repr%pdecomp%shape_out(jr) = &
            repr%pdecomp%shape_mem( repr%order_mem2out(jr) )
    END DO

    ! PARALLEL IO TYPE
    IF (PRESENT(piotype)) THEN
       repr%pdecomp%piotype = piotype
    END IF

    ! UPDATE
    repr%pdecomp%lpdecomp = .TRUE.
    status = 0

  END SUBROUTINE set_representation_decomp
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_representation_by_name(status, name, repr)

    IMPLICIT NONE

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    CHARACTER(LEN=*),       INTENT(IN)  :: name
    TYPE(t_representation), POINTER     :: repr

    CALL loc_representation_by_name(status, GREPRLIST, name, repr)

  END SUBROUTINE get_representation_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_representation_by_id(status, id, repr)

    IMPLICIT NONE

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    INTEGER,                INTENT(IN)  :: id
    TYPE(t_representation), POINTER     :: repr

    CALL loc_representation_by_id(status, GREPRLIST, id, repr)

  END SUBROUTINE get_representation_by_id
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_representation_id(status, name, reprid)

    IMPLICIT NONE

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    CHARACTER(LEN=*),       INTENT(IN)  :: name
    INTEGER,                INTENT(OUT) :: reprid

    ! LOCAL
    TYPE(t_representation), POINTER     :: repr

    reprid = REPR_UNDEF

    CALL loc_representation_by_name(status, GREPRLIST, name, repr)
    IF (status /= 0) RETURN

    reprid = repr%id
    status = 0

  END SUBROUTINE get_representation_id
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE clean_representations(status)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER, INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_representation_list), POINTER :: ai => NULL()
    TYPE(t_representation_list), POINTER :: ae => NULL()

    ai => GREPRLIST
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       !
       ae => ai
       ai => ai%next
       !
       IF (ae%this%pdecomp%lpdecomp) THEN
          DEALLOCATE(ae%this%pdecomp%start) ; NULLIFY(ae%this%pdecomp%start)
          DEALLOCATE(ae%this%pdecomp%cnt)   ; NULLIFY(ae%this%pdecomp%cnt)
          DEALLOCATE(ae%this%pdecomp%ml)    ; NULLIFY(ae%this%pdecomp%ml)
          DEALLOCATE(ae%this%pdecomp%mu)    ; NULLIFY(ae%this%pdecomp%mu)
       END IF
       !
       DEALLOCATE(ae)
       NULLIFY(ae)
       !
       ! COUNT
       NREP = NREP - 1
    END DO    

    NULLIFY(GREPRLIST)

    IF (NREP /= 0) THEN
       status = 1003
       RETURN
    END IF
    
    status = 0
    
  END SUBROUTINE clean_representations
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_representation_all(status)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER, INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_representation_list), POINTER :: ai => NULL()

    WRITE(*,*) '### REPRESENTATIONS: ###############################'

    ! EMPTY LIST ?
    IF (.NOT. ASSOCIATED(GREPRLIST)) THEN

       WRITE(*,*) '*** REPRESENTATION LIST IS EMPTY ***'

    ELSE

       ai => GREPRLIST
       DO
          IF (.NOT. ASSOCIATED(ai)) EXIT
          !
          CALL write_representation_one(status, ai%this)
          IF (status /= 0) RETURN
          !
          ai => ai%next
       END DO
       
    END IF

    WRITE(*,*) '####################################################'
    
    status = 0

  END SUBROUTINE write_representation_all
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_representation_one(status, repr)

    USE messy_main_channel_dimensions, ONLY: write_dimension

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    TYPE(t_representation), INTENT(IN)  :: repr

    ! LOCAL
    INTEGER :: i

    WRITE(*,*) '# NAME                   : ',repr%name
    WRITE(*,*) '# ID                     : ',repr%id
    WRITE(*,*) '# RANK                   : ',repr%rank
    WRITE(*,*) '# LINK                   : ',repr%link
    WRITE(*,*) '# AXIS                   : ',repr%axis
    WRITE(*,*) '# DECOMPOSITION TYPE     : ',repr%dctype
    WRITE(*,*) '# GLOBAL DIMENSION LENGTH: ',repr%gdimlen
    WRITE(*,*) '#     SHAPE IN OUTPUT    : ',repr%shape_out
    WRITE(*,*) '#     ORDER mem2out      : ',repr%order_mem2out
    WRITE(*,*) '#     ORDER out2mem      : ',repr%order_out2mem 
    WRITE(*,*) '#     PERMUTATION(S)     : ',repr%lperm
    WRITE(*,*) '# LOCAL  DIMENSION LENGTH: ',repr%ldimlen
    !
    IF (repr%pdecomp%lpdecomp) THEN
       WRITE(*,*) '# PARALLEL DECOMP.-INFO  : YES'
       SELECT CASE(repr%pdecomp%piotype)
       CASE(PIOTYPE_SGL)
          WRITE(*,*) '# PARALLEL I/O TYPE      : SINGLE'
       CASE(PIOTYPE_IND)
          WRITE(*,*) '# PARALLEL I/O TYPE      : INDEPENDENT'
       CASE(PIOTYPE_COL)
          WRITE(*,*) '# PARALLEL I/O TYPE      : COLLECTIVE'
       CASE DEFAULT
       END SELECT
    ELSE
       WRITE(*,*) '# PARALLEL DECOMP.-INFO  : NO'
    ENDIF
    !
    ! mz_ab_20090921+
    WRITE(*,*) '# BOUNDARIES             : ',repr%bounds%nbounds(1:repr%rank)
    ! mz_ab_20090921-
    !
    WRITE(*,*) '# DIMENSIONS:'
    DO i=1, IRANK
       IF (ASSOCIATED(repr%dim(i)%ptr)) THEN
          CALL write_dimension(status, repr%dim(i)%ptr)
          IF (status /= 0) RETURN
       ELSE
          WRITE(*,*) '   |-> *** UNDEFINED DIMENSION ***'
       END IF
    END DO

    status = 0

  END SUBROUTINE write_representation_one
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE write_representation_dc(status, p_pe)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: p_pe

    ! LOCAL
    TYPE(t_representation_list), POINTER :: ai => NULL()
    INTEGER :: i

    WRITE(*,*) '### REPRESENTATIONS DECOMPOSITION TABLE: ###########'

    ! EMPTY LIST ?
    IF (.NOT. ASSOCIATED(GREPRLIST)) THEN

       WRITE(*,*) '*** REPRESENTATION LIST IS EMPTY ***'

    ELSE

       ai => GREPRLIST
       DO
          IF (.NOT. ASSOCIATED(ai)) EXIT
          !
          IF (ai%this%pdecomp%lpdecomp) &
               WRITE(*,*) 'DECOMPOSITION OF REPRESENTATION ' &
               ,TRIM(ai%this%name), &
               ' ON PE ',p_pe,': ', &
               ' [shape_mem] ',ai%this%pdecomp%shape_mem, &
               ' [shape_out] ',ai%this%pdecomp%shape_out, &
               ( &
               ' [START] ',ai%this%pdecomp%start(i,:), &
               ' [CNT  ] ',ai%this%pdecomp%cnt(i,:),   &
               ' [ML   ] ',ai%this%pdecomp%ml(i,:),  &
               ' [MU   ] ',ai%this%pdecomp%mu(i,:),  &
               i=1, ai%this%pdecomp%nseg )
          !
          ai => ai%next
       END DO
       
    END IF

    WRITE(*,*) '####################################################'
    
    status = 0

  END SUBROUTINE write_representation_dc
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE repr_reorder(status, flag, lparallel, repr, mem, out)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, RESHAPE

    ! I/O
    INTEGER,                      INTENT(OUT)   :: status
    INTEGER,                      INTENT(IN)    :: flag
    LOGICAL,                      INTENT(IN)    :: lparallel
    TYPE(t_representation),       POINTER       :: repr
    REAL(DP), DIMENSION(:,:,:,:), POINTER       :: mem   ! MEMORY
    REAL(DP), DIMENSION(:,:,:,:), POINTER       :: out   ! OUTPUT

    ! LOCAL
    INTEGER :: zstat
    INTEGER,DIMENSION(IRANK) :: shape_mem    ! shape in memory = gdimlen
    INTEGER,DIMENSION(IRANK) :: order_mem2out   ! output
    INTEGER,DIMENSION(IRANK) :: shape_out       ! shape in output
    INTEGER,DIMENSION(IRANK) :: order_out2mem   ! input

    ! INIT
    IF (.NOT. lparallel) THEN
       shape_mem(:)     = repr%gdimlen(:)
       shape_out(:)     = repr%shape_out(:)
    ELSE
       shape_mem(:)     = repr%pdecomp%shape_mem(:)
       shape_out(:)     = repr%pdecomp%shape_out(:)
    END IF
    !
    order_mem2out(:) = repr%order_mem2out(:)
    order_out2mem(:) = repr%order_out2mem(:)

    SELECT CASE(flag)
    CASE(1)
       ! ################# MEMORY -> OUTPUT ##############################
       !
       ! INIT
       IF (ASSOCIATED(out)) THEN
          DEALLOCATE(out)
          NULLIFY(out)
       END IF
       !
       ALLOCATE(out(shape_out(1),shape_out(2),shape_out(3),shape_out(4)) &
            , STAT=zstat)
       IF (zstat /= 0) THEN
          status = 1000 ! MEMORY ALLOCATION FAILED
          RETURN
       END IF
       !
       IF (repr%lperm) THEN
#ifdef GENERIC
          out = RESHAPE(mem, shape=shape_out, order=order_out2mem)
#else
          CALL FAST_REORDER(status, out, mem, repr%i_m2o)
          IF (status /= 0) RETURN
#endif
       ELSE
!!$          out(:,:,:,:) = mem(:,:,:,:)
          out = mem
       END IF
       !
       status = 0
       !
       ! #################################################################
    CASE(-1)
       ! ################# MEMORY <- OUTPUT ##############################
       !
       ! INIT
       IF (ASSOCIATED(mem)) THEN
          DEALLOCATE(mem)
          NULLIFY(mem)
       END IF
       !
       ALLOCATE(mem(shape_mem(1), shape_mem(2), shape_mem(3), shape_mem(4)) &
            , STAT=zstat)
       !
       IF (zstat /= 0) THEN
          status = 1000 ! MEMORY ALLOCATION FAILED
          RETURN
       END IF
       !
       IF (repr%lperm) THEN
#ifdef GENERIC
          mem = RESHAPE(out, shape=shape_mem, order=order_mem2out)
#else
          CALL FAST_REORDER(status, mem, out, repr%i_o2m)
          IF (status /= 0) RETURN
#endif
       ELSE
!!$          mem(:,:,:,:) = out(:,:,:,:)
          mem = out
       END IF
       !
       status = 0
       !
       ! #################################################################
    CASE DEFAULT
       !
       status = 2015 ! INVALID REPR_REORDER FLAG
       !
    END SELECT

#ifndef GENERIC
    CONTAINS

      SUBROUTINE FAST_REORDER(status, out, inp, order)

        IMPLICIT NONE

        ! I/O
        INTEGER,                      INTENT(OUT) :: status
        REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: out
        REAL(DP), DIMENSION(:,:,:,:), INTENT(IN)  :: inp
        INTEGER,                      INTENT(IN)  :: order

        ! LOCAL
        INTEGER :: i1, i2, i3, i4
        INTEGER :: n1, n2, n3, n4

        n1 = SIZE(inp, 1)
        n2 = SIZE(inp, 2)
        n3 = SIZE(inp, 3)
        n4 = SIZE(inp, 4)

        ! no checks of size(out) (performance !)

        SELECT CASE(order)

        CASE(1324)
           DO i1=1, n1
              DO i2=1, n2
                 DO i3=1, n3
                    DO i4=1, n4
                       out(i1, i3, i2, i4) = inp(i1, i2, i3, i4)
                    END DO
                 END DO
              END DO
           END DO

        CASE(2314)
           DO i1=1, n1
              DO i2=1, n2
                 DO i3=1, n3
                    DO i4=1, n4
                       out(i2, i3, i1, i4) = inp(i1, i2, i3, i4)
                    END DO
                 END DO
              END DO
           END DO

        CASE(3124)
           DO i1=1, n1
              DO i2=1, n2
                 DO i3=1, n3
                    DO i4=1, n4
                       out(i3, i1, i2, i4) = inp(i1, i2, i3, i4)
                    END DO
                 END DO
              END DO
           END DO

        CASE(3142)
           DO i1=1, n1
              DO i2=1, n2
                 DO i3=1, n3
                    DO i4=1, n4
                       out(i3, i1, i4, i2) = inp(i1, i2, i3, i4)
                    END DO
                 END DO
              END DO
           END DO

        CASE(3214)
           DO i1=1, n1
              DO i2=1, n2
                 DO i3=1, n3
                    DO i4=1, n4
                       out(i3, i2, i1, i4) = inp(i1, i2, i3, i4)
                    END DO
                 END DO
              END DO
           END DO           

        CASE(2134)
           DO i1=1, n1
              DO i2=1, n2
                 DO i3=1, n3
                    DO i4=1, n4
                       out(i2, i1, i3, i4) = inp(i1, i2, i3, i4)
                    END DO
                 END DO
              END DO
           END DO     

           ! op_bk_20131118+
        CASE(3421)
           DO i1=1, n1
              DO i2=1, n2
                 DO i3=1, n3
                    DO i4=1, n4
                       out(i3, i4, i2, i1) = inp(i1, i2, i3, i4)
                    END DO
                 END DO
              END DO
           END DO     
           ! op_bk_20131118-

        CASE(2413)
           DO i1=1, n1
              DO i2=1, n2
                 DO i3=1, n3
                    DO i4=1, n4
                       out(i2, i4, i1, i3) = inp(i1, i2, i3, i4)
                    END DO
                 END DO
              END DO
           END DO     

        CASE DEFAULT
           status = 2100 ! NON-GENERIC FAST REORDER NOT IMPLEMENTED
                         ! FOR SPECIFIED PERMUTATION
           RETURN
        END SELECT

        status = 0

      END SUBROUTINE FAST_REORDER
#endif

  END SUBROUTINE repr_reorder
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE repr_getptr(status, repr, in, p0, p1, p2, p3, p4, linner)

    USE messy_main_tools, ONLY: remap_bounds  ! mz_ab_20100610

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRESENT

    ! I/O
    INTEGER,                      INTENT(OUT)          :: status
    REAL(DP), DIMENSION(:,:,:,:), POINTER              :: in
    TYPE(t_representation),       POINTER              :: repr
    REAL(DP),                     POINTER,    OPTIONAL :: p0
    REAL(DP), DIMENSION(:),       POINTER,    OPTIONAL :: p1
    REAL(DP), DIMENSION(:,:),     POINTER,    OPTIONAL :: p2
    REAL(DP), DIMENSION(:,:,:),   POINTER,    OPTIONAL :: p3
    REAL(DP), DIMENSION(:,:,:,:), POINTER,    OPTIONAL :: p4
    ! mz_ab_20100211+
    LOGICAL,                      INTENT(IN), OPTIONAL :: linner
    
    LOGICAL :: zlinner 
    INTEGER :: lb1, lb2, lb3, lb4

    IF (PRESENT(linner)) THEN
       zlinner = linner
    ELSE
       zlinner = .FALSE.   ! default
    ENDIF   
    ! mz_ab_20100211-

    ! POINTER ASSOCIATION
 
    IF ((.NOT. repr%bounds%lbounds) .OR. zlinner) THEN ! mz_ab_20100222
       ! .NOT. (A .AND. (.NOT. B) ) => (.NOT. A) .OR. B
       SELECT CASE(repr%link)
       CASE('----')
          IF (PRESENT(p0)) p0 => in(1,1,1,1)
          IF (PRESENT(p1)) p1 => in(:,1,1,1)
          IF (PRESENT(p2)) p2 => in(:,:,1,1)
          IF (PRESENT(p3)) p3 => in(:,:,:,1)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)

       CASE('x---')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => in(:,1,1,1)
          IF (PRESENT(p2)) p2 => in(:,:,1,1)
          IF (PRESENT(p3)) p3 => in(:,:,:,1)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('-x--')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => in(1,:,1,1)
          IF (PRESENT(p2)) p2 => in(:,:,1,1)
          IF (PRESENT(p3)) p3 => in(:,:,:,1)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('--x-')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => in(1,1,:,1)
          IF (PRESENT(p2)) p2 => in(:,1,:,1)
          IF (PRESENT(p3)) p3 => in(:,:,:,1)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('---x')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => in(1,1,1,:)
          IF (PRESENT(p2)) p2 => in(:,1,1,:)
          IF (PRESENT(p3)) p3 => in(:,:,1,:)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)

       CASE('xx--')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => in(:,:,1,1)
          IF (PRESENT(p3)) p3 => in(:,:,:,1)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('x-x-')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => in(:,1,:,1)
          IF (PRESENT(p3)) p3 => in(:,:,:,1)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('x--x')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => in(:,1,1,:)
          IF (PRESENT(p3)) p3 => in(:,:,1,:)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('-xx-')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => in(1,:,:,1)
          IF (PRESENT(p3)) p3 => in(:,:,:,1)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('-x-x')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => in(1,:,1,:)
          IF (PRESENT(p3)) p3 => in(:,:,1,:)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('--xx')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => in(1,1,:,:)
          IF (PRESENT(p3)) p3 => in(:,1,:,:)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)

       CASE('-xxx')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => NULL()
          IF (PRESENT(p3)) p3 => in(1,:,:,:)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('x-xx')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => NULL()
          IF (PRESENT(p3)) p3 => in(:,1,:,:)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('xx-x')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => NULL()
          IF (PRESENT(p3)) p3 => in(:,:,1,:)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)
       CASE('xxx-')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => NULL()
          IF (PRESENT(p3)) p3 => in(:,:,:,1)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)

       CASE('xxxx')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => NULL()
          IF (PRESENT(p3)) p3 => NULL()
          IF (PRESENT(p4)) p4 => in(:,:,:,:)

       CASE DEFAULT
       END SELECT

       ! mz_ab_20100211+
    ELSE ! RETURN FIELD INCL BOUNDS AS E.G. (-1:...)

       lb1 = 1-repr%bounds%nbounds(1) 
       lb2 = 1-repr%bounds%nbounds(2) 
       lb3 = 1-repr%bounds%nbounds(3) 
       lb4 = 1-repr%bounds%nbounds(4) 

       SELECT CASE(repr%link)
       CASE('----')
          IF (PRESENT(p0)) p0 => in(1,1,1,1)
          IF (PRESENT(p1)) p1 => in(:,1,1,1)
          IF (PRESENT(p2)) p2 => in(:,:,1,1)
          IF (PRESENT(p3)) p3 => in(:,:,:,1)
          IF (PRESENT(p4)) p4 => in(:,:,:,:)

       CASE('x---')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => remap_bounds(lb1,                in(:,1,1,1))
          IF (PRESENT(p2)) p2 => remap_bounds(lb1, lb2,           in(:,:,1,1))
          IF (PRESENT(p3)) p3 => remap_bounds(lb1, lb2, lb3,      in(:,:,:,1))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1, lb2, lb3, lb4, in(:,:,:,:))
       CASE('-x--')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => remap_bounds(lb1,                in(1,:,1,1))
          IF (PRESENT(p2)) p2 => remap_bounds(lb1, lb2,           in(:,:,1,1))
          IF (PRESENT(p3)) p3 => remap_bounds(lb1, lb2, lb3,      in(:,:,:,1))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1, lb2, lb3, lb4, in(:,:,:,:))
       CASE('--x-')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => remap_bounds(lb1,                in(1,1,:,1))
          IF (PRESENT(p2)) p2 => remap_bounds(lb1, lb2,           in(:,1,:,1))
          IF (PRESENT(p3)) p3 => remap_bounds(lb1, lb2, lb3,      in(:,:,:,1))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1, lb2, lb3, lb4, in(:,:,:,:))
       CASE('---x')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => remap_bounds(lb1,                in(1,1,1,:))
          IF (PRESENT(p2)) p2 => remap_bounds(lb1, lb2,           in(:,1,1,:))
          IF (PRESENT(p3)) p3 => remap_bounds(lb1, lb2, lb3,      in(:,:,1,:))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1, lb2, lb3, lb4, in(:,:,:,:))

       CASE('xx--')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => remap_bounds(lb1,lb2,        in(:,:,1,1))
          IF (PRESENT(p3)) p3 => remap_bounds(lb1,lb2,lb3,    in(:,:,:,1))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:))
       CASE('x-x-')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => remap_bounds(lb1,lb2,        in(:,1,:,1))
          IF (PRESENT(p3)) p3 => remap_bounds(lb1,lb2,lb3,    in(:,:,:,1))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:))
       CASE('x--x')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => remap_bounds(lb1,lb2,        in(:,1,1,:))
          IF (PRESENT(p3)) p3 => remap_bounds(lb1,lb2,lb3,    in(:,:,1,:))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:))
       CASE('-xx-')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => remap_bounds(lb1,lb2,        in(1,:,:,1))
          IF (PRESENT(p3)) p3 => remap_bounds(lb1,lb2,lb3,    in(:,:,:,1))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:))
       CASE('-x-x')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => remap_bounds(lb1,lb2,        in(1,:,1,:))
          IF (PRESENT(p3)) p3 => remap_bounds(lb1,lb2,lb3,    in(:,:,1,:))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:))
       CASE('--xx')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => remap_bounds(lb1,lb2,        in(1,1,:,:))
          IF (PRESENT(p3)) p3 => remap_bounds(lb1,lb2,lb3,    in(:,1,:,:))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:)) 
       CASE('-xxx')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => NULL()
          IF (PRESENT(p3)) p3 => remap_bounds(lb1,lb2,lb3,    in(1,:,:,:))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:))
       CASE('x-xx')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => NULL()
          IF (PRESENT(p3)) p3 => remap_bounds(lb1,lb2,lb3,    in(:,1,:,:))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:))
       CASE('xx-x')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => NULL()
          IF (PRESENT(p3)) p3 => remap_bounds(lb1,lb2,lb3,    in(:,:,1,:))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:))
       CASE('xxx-')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => NULL()
          IF (PRESENT(p3)) p3 => remap_bounds(lb1,lb2,lb3,    in(:,:,:,1))
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:))

       CASE('xxxx')
          IF (PRESENT(p0)) p0 => NULL()
          IF (PRESENT(p1)) p1 => NULL()
          IF (PRESENT(p2)) p2 => NULL()
          IF (PRESENT(p3)) p3 => NULL()
          IF (PRESENT(p4)) p4 => remap_bounds(lb1,lb2,lb3,lb4,in(:,:,:,:))

       CASE DEFAULT
       END SELECT

    ENDIF
    ! mz_ab_20100211-

    IF (PRESENT(p4)) THEN
       IF (.NOT.ASSOCIATED(p4)) THEN
          status = 2021 ! POINTER p4 COULD NOT BE ASSOCIATED
          RETURN
       END IF
    END IF

    IF (PRESENT(p3)) THEN
       IF (.NOT.ASSOCIATED(p3)) THEN
          status = 2020 ! POINTER p3 COULD NOT BE ASSOCIATED
          RETURN
       END IF
    END IF

    IF (PRESENT(p2)) THEN
       IF (.NOT.ASSOCIATED(p2)) THEN
          status = 2019 ! POINTER p2 COULD NOT BE ASSOCIATED
          RETURN
       END IF
    END IF

    IF (PRESENT(p1)) THEN
       IF (.NOT.ASSOCIATED(p1)) THEN
          status = 2018 ! POINTER p1 COULD NOT BE ASSOCIATED
          RETURN
       END IF
    END IF

    IF (PRESENT(p0)) THEN
       IF (.NOT.ASSOCIATED(p0)) THEN
          status = 2017 ! POINTER p0 COULD NOT BE ASSOCIATED
          RETURN
       END IF
    END IF

    status = 0

  END SUBROUTINE repr_getptr
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_representation_by_name(status, list, name, repr)

    IMPLICIT NONE

    INTRINSIC :: ADJUSTL, ASSOCIATED, LEN_TRIM, TRIM

    ! I/O
    INTEGER,                     INTENT(OUT)          :: status
    TYPE(t_representation_list), POINTER              :: list
    CHARACTER(LEN=*),            INTENT(IN)           :: name
    TYPE(t_representation),      POINTER              :: repr

    ! LOCAL
    TYPE(t_representation_list), POINTER :: ai  => NULL()
    TYPE(t_representation_list), POINTER :: ae  => NULL()
    LOGICAL                              :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(repr)

    ! CHECKS
    IF (LEN_TRIM(ADJUSTL(name)) > STRLEN_MEDIUM) THEN
       WRITE(*,*) "MESSY !! ", name
       status = 2001  ! REPRESENTATION NAME TOO LONG
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(list)) THEN
       status = 2003  ! REPRESENTATION DOES NOT EXIST
       RETURN
    END IF

    ! CHECK, IF IT EXISTS
    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF (TRIM(ADJUSTL(name)) == TRIM(ai%this%name)) THEN
          lex = .TRUE.
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO    

    IF (lex) THEN
       repr => ai%this
    ELSE
       status = 2003  ! REPRESENTATION DOES NOT EXIST
       RETURN
    END IF

    status = 0


  END SUBROUTINE loc_representation_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE loc_representation_by_id(status, list, id, repr)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,                     INTENT(OUT)          :: status
    TYPE(t_representation_list), POINTER              :: list
    INTEGER,                     INTENT(IN)           :: id
    TYPE(t_representation),      POINTER              :: repr

    ! LOCAL
    TYPE(t_representation_list), POINTER :: ai  => NULL()
    TYPE(t_representation_list), POINTER :: ae  => NULL()
    LOGICAL                              :: lex

    ! INIT
    lex = .FALSE.
    NULLIFY(repr)

    ! CHECKS
    IF (id <= 0) THEN
       status = 2013  ! INVALID REPRESENTATION ID
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(list)) THEN
       status = 2003  ! REPRESENTATION DOES NOT EXIST
       RETURN
    END IF

    ! CHECK, IF IT EXISTS
    ai => list
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT
       IF (id == ai%this%id) THEN
          lex = .TRUE.
          EXIT
       END IF
       ae => ai
       ai => ai%next
    END DO    

    IF (lex) THEN
       repr => ai%this
    ELSE
       status = 2003  ! REPRESENTATION DOES NOT EXIST
       RETURN
    END IF

    status = 0

  END SUBROUTINE loc_representation_by_id
  ! -------------------------------------------------------------------

  ! um_ak_20090729+
  ! -------------------------------------------------------------------
  SUBROUTINE get_representation_info(status, inpname, Id, rank, link, axis &
                                    , gdimlen, ldimlen, dctype, name       &
                                    , nbounds) ! op_pj_20100121

    IMPLICIT NONE

    INTEGER,                      INTENT(OUT)             :: status
    CHARACTER(LEN=*),             INTENT(IN)              :: inpname
    INTEGER,                      INTENT(INOUT), OPTIONAL :: ID
    INTEGER,                      INTENT(OUT),   OPTIONAL :: rank
    CHARACTER(LEN=IRANK),         INTENT(OUT),   OPTIONAL :: link
    CHARACTER(LEN=IRANK),         INTENT(OUT),   OPTIONAL :: axis 
    INTEGER, DIMENSION(IRANK),    INTENT(OUT),   OPTIONAL :: gdimlen
    INTEGER, DIMENSION(IRANK),    INTENT(OUT),   OPTIONAL :: ldimlen
    INTEGER,                      INTENT(OUT),   OPTIONAL :: dctype 
    CHARACTER(LEN=STRLEN_MEDIUM), INTENT(OUT),   OPTIONAL :: name
    ! op_pj_20100121+
    INTEGER, DIMENSION(IRANK),    INTENT(OUT),   OPTIONAL :: nbounds
    ! op_pj_20100121-

    ! LOCAL
    TYPE(t_representation), POINTER  :: repr

    IF (TRIM(ADJUSTL(inpname)) /= '') THEN
       CALL loc_representation_by_name(status, GREPRLIST, inpname, repr)
       IF (status/=0) RETURN
       IF (PRESENT(id)) id = repr%id    
    ELSE IF (PRESENT(id)) THEN
       CALL loc_representation_by_id(status, GREPRLIST, ID, repr)
       IF (status/=0) RETURN
    ELSE
       status = 2027
       RETURN
    ENDIF

    IF (PRESENT(rank))    rank    = repr%rank
    IF (PRESENT(link))    link    = repr%link
    IF (PRESENT(axis))    axis    = repr%axis
    IF (PRESENT(gdimlen)) gdimlen = repr%gdimlen
    IF (PRESENT(ldimlen)) ldimlen = repr%ldimlen
    IF (PRESENT(dctype))  dctype  = repr%dctype
    IF (PRESENT(name))    name    = repr%name
    ! op_pj_20100121+
    IF (PRESENT(nbounds)) nbounds = repr%bounds%nbounds
    ! op_pj_20100121-

    status = 0

  END SUBROUTINE get_representation_info
  ! -------------------------------------------------------------------
  ! um_ak_20090729-

! **********************************************************************
END MODULE messy_main_channel_repr
! **********************************************************************
