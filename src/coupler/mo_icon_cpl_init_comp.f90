!>
!! Initialisation of the ICON coupler
!!
!! <Description>
!! 
!! This routine returns a grid_id and sets up a MPI communicator
!! which can be used fo local communication inside the model component.
!! This communicator can be accessed via 
!!
!! @author Rene Redler, MPI-M
!!
!! $Id:$
!!
!! @par Revision History
!! first implementation by Rene Redler (2010-02-13)
!!
!! @par Copyright
!! 2010-2011 by MPI-M
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
#if ! defined (__INTEL_COMPILER) && ! defined (__SX__) && ! defined (__PGI)
#define __NAG
#endif
#if ! defined (__xlC)            && ! defined (__sun)  && ! defined (__SUNPRO_F95)
#define __NAG
#endif
#if ! defined (__GFORTRAN__)
#define __NAG
#endif

MODULE mo_icon_cpl_init_comp

#ifndef NOMPI
  USE mpi, ONLY : MPI_SUCCESS, MPI_COMM_NULL
#endif

  USE mo_icon_cpl, ONLY : set_cpl_local_comm, maxchar,        &
   &                      l_MPI_was_initialized,              &
   &                      l_debug, cplout,                    &
   &                      grids, comps, cpl_fields,           &
   &                      nbr_active_comps,                   &
   &                      nbr_active_grids,                   &
   &                      nbr_active_fields,                  &
   &                      nbr_ICON_fields, nbr_ICON_comps,    &
   &                      nbr_ICON_grids, nbr_ICON_inc,       &
   &                      ICON_comm,                          &
   &                      ICON_global_rank, ICON_global_size, &
   &                      ICON_local_rank, ICON_local_size

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER    :: version = '$Id$'

  LOGICAL                        :: l_MPI_is_initialized ! to check whether MPI_init was called.

  ! Return code for error handling

  CHARACTER(len=132)             :: err_string ! Error string for MPI
  INTEGER                        :: len        ! length 
  INTEGER                        :: ierr       ! returned error from MPI functions

  PRIVATE

  PUBLIC :: icon_cpl_init_comp, get_my_local_comp_id, icon_cpl_redirect_stdout

CONTAINS
  

  SUBROUTINE icon_cpl_init_comp ( comp_name, global_comp_no, global_comp_type, &
    comp_id, ierror )

    CHARACTER(len=*), INTENT(in) :: comp_name
    INTEGER, INTENT(in)          :: global_comp_no   ! the component unique number
    INTEGER, INTENT(in)          :: global_comp_type ! the component type (ocean, atmo, etc)
    
    INTEGER, INTENT(out)         :: comp_id
    INTEGER, INTENT(out)         :: ierror

    INTEGER                      :: comp_comm
    INTEGER                      :: key, color

    INTEGER                      :: i               !< loop count


    ! -------------------------------------------------------------------
    ! Initialise variables
    ! -------------------------------------------------------------------

    ierror  = 0

#ifndef NOMPI

    ! -------------------------------------------------------------------
    ! Check MPI Initialization
    ! -------------------------------------------------------------------

    l_MPI_was_initialized = .FALSE.

    CALL MPI_Initialized ( l_MPI_is_initialized, ierr )

    IF ( .NOT. l_MPI_is_initialized ) THEN

       CALL MPI_Init ( ierr )

       IF ( ierr /= MPI_SUCCESS ) THEN
          CALL MPI_Error_string ( ierr, err_string, len, ierror )
          WRITE  ( * , '(A,A)' ) 'Error in MPI_Init ', err_string
       ENDIF

       l_MPI_was_initialized = .TRUE.

    ENDIF

    ! -------------------------------------------------------------------
    ! Get a local component ID
    ! -------------------------------------------------------------------

    DO comp_id = 1, nbr_ICON_comps
       IF ( .NOT. comps(comp_id)%l_comp_status ) EXIT
    ENDDO

    IF ( comp_id > nbr_ICON_comps ) THEN
       PRINT *, 'number of requested components exceeds maximum of nbr_ICON_comps'
       CALL MPI_Abort ( ICON_comm, 1, ierr )
    ENDIF

    ! -------------------------------------------------------------------
    ! Derive component communicators
    ! -------------------------------------------------------------------

    color = global_comp_no
    key   = 0

    comp_comm = MPI_COMM_NULL

    CALL MPI_Comm_split ( ICON_comm, color, key, comp_comm, ierr )
    IF ( ierr /= MPI_SUCCESS ) THEN
       CALL MPI_Error_string ( ierr, err_string, len, ierror )
       WRITE  ( * , '(A14,I3,A)' ) 'Error on rank ', ICON_global_rank, err_string
    ENDIF

    IF ( l_debug ) &
         WRITE ( cplout , '(I3,A1,A)' ) ICON_global_rank, ':', ' returned from split '

    ! -------------------------------------------------------------------
    ! Size of component and rank in component
    ! -------------------------------------------------------------------

    CALL MPI_Comm_rank ( comp_comm, ICON_local_rank, ierr )
    IF ( ierr /= MPI_SUCCESS ) THEN
       CALL MPI_Error_string ( ierr, err_string, len, ierror )
       WRITE  ( * , '(A,A)' ) 'Error in getting local rank ', err_string
    ENDIF

    CALL MPI_Comm_size ( comp_comm, ICON_local_size, ierr )
    IF ( ierr /= MPI_SUCCESS ) THEN
       CALL MPI_Error_string ( ierr, err_string, len, ierror )
       WRITE  ( * , '(A,A)' ) 'Error in getting local size ', err_string
    ENDIF

    ! -------------------------------------------------------------------
    ! Initialize components 
    ! -------------------------------------------------------------------

    ! Are there compilers that associate pointer at initialisation?

    IF ( .NOT. ASSOCIATED(comps) ) THEN

       nbr_active_comps  = 0

       ALLOCATE ( comps(nbr_ICON_comps), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating comps '
          CALL MPI_Abort ( ICON_COMM, 1, ierr )
       ENDIF

       comps(:)%l_comp_status = .FALSE.

    ENDIF

    comps(comp_id)%comp_name     = TRIM(comp_name)
    comps(comp_id)%l_comp_status = .TRUE.

    nbr_active_comps = nbr_active_comps + 1

    CALL icon_cpl_redirect_stdout ( comp_id )

    ! -------------------------------------------------------------------
    ! Initialize grids
    ! -------------------------------------------------------------------

    IF ( .NOT. ASSOCIATED(grids) ) THEN

       nbr_active_grids  = 0

       nbr_ICON_grids = nbr_ICON_inc

       ALLOCATE ( grids(nbr_ICON_grids), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating grids '
          CALL MPI_Abort ( ICON_COMM, 1, ierr )
       ENDIF

       grids(:)%comp_id       = -1
       grids(:)%l_grid_status = .FALSE.

    ENDIF

    ! -------------------------------------------------------------------
    ! Initialize fields
    ! -------------------------------------------------------------------

    IF ( .NOT. ASSOCIATED(cpl_fields) ) THEN

       nbr_active_fields = 0

       nbr_ICON_fields = nbr_ICON_inc

       ALLOCATE ( cpl_fields(nbr_ICON_fields), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating cpl_fields '
          CALL MPI_Abort ( ICON_COMM, 1, ierr )
       ENDIF

       cpl_fields(:)%comp_id        = -1
       cpl_fields(:)%grid_id        = -1
       cpl_fields(:)%event_id       = -1
       cpl_fields(:)%l_field_status = .FALSE.

       cpl_fields(:)%coupling%lag            = 0
       cpl_fields(:)%coupling%time_operation = 0
       cpl_fields(:)%coupling%frequency      = 0
       cpl_fields(:)%coupling%time_step      = 0

       DO i = 1, nbr_ICON_fields
          NULLIFY ( cpl_fields(i)%send_field_acc )
       ENDDO

    ENDIF

#else

    ! -------------------------------------------------------------------
    ! Initialise variables
    ! -------------------------------------------------------------------

    i                        = 0
    color                    = global_comp_no
    key                      = 0
    ierror                   = 0
    comp_id                  = 1
    comp_comm                = 0
    comps(comp_id)%comp_name = TRIM(comp_name)

#endif

    ierror = set_cpl_local_comm ( comp_comm )

  END SUBROUTINE icon_cpl_init_comp

  ! -------------------------------------------------------------------
  ! Provide the local component id
  ! -------------------------------------------------------------------

  INTEGER FUNCTION get_my_local_comp_id ( comp_process )

     INTEGER, INTENT (IN) :: comp_process

     INTEGER :: i

     DO i = 1, nbr_active_comps
       IF ( comps(i)%comp_process == comp_process ) EXIT
     ENDDO

     IF ( i > nbr_active_comps ) THEN
        ! TODO:
        ! We have an error and need to stop
     ENDIF

     get_my_local_comp_id = i
     
  END FUNCTION get_my_local_comp_id

  ! -------------------------------------------------------------------
  ! Redirect standard output and standard error output if requested
  ! -------------------------------------------------------------------

  SUBROUTINE icon_cpl_redirect_stdout ( comp_id )

    INTEGER, INTENT (IN) :: comp_id

    CHARACTER(len=maxchar)       :: comp_name
    INTEGER                      :: lenstr
#if ! defined (__NAG)
    INTEGER                      :: istat
    INTEGER                      :: length_of_integer
    INTEGER                      :: ndibuf
    INTEGER                      :: ipos
    INTEGER, ALLOCATABLE         :: ibuf(:)
#endif

    ! TODO: replace with user namelist input

    LOGICAL                      :: l_redirect_stdout 

    INTEGER, PARAMETER           :: parallel_io = 1 !< switch for writing cplout
                                                    !<  - 1 write in parallel
                                                    !<  - 0 all in one

    l_redirect_stdout = .FALSE.

    IF ( l_redirect_stdout ) THEN

       comp_name = comps(comp_id)%comp_name

       lenstr = LEN_TRIM (comp_name)

#if defined (__NAG)
       CALL psmile_redirstdout ( comp_name(1:lenstr), lenstr, &
            parallel_io, ICON_global_rank, ICON_global_size, ierr )
#else
       length_of_integer = BIT_SIZE (istat) / 8
       ndibuf = lenstr / length_of_integer + 1

       ALLOCATE (ibuf(1:ndibuf), STAT = ierr)
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating ibuf '
          CALL MPI_Abort ( ICON_COMM, 1, ierr )
       ENDIF

       ipos = 0

       CALL psmile_char2buf (ibuf, ndibuf, ipos, comp_name(1:lenstr))

       CALL psmile_redirstdout ( ibuf, lenstr, &
            parallel_io, ICON_global_rank, ICON_global_size, ierr )
#endif
    ENDIF

  END SUBROUTINE icon_cpl_redirect_stdout

END MODULE mo_icon_cpl_init_comp
