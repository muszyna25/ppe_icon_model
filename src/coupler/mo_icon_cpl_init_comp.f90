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

  USE mo_cpl_nml,  ONLY  : cpl_nml_setup, l_redirect_stdout

  USE mo_icon_cpl, ONLY : set_cpl_local_comm

#ifndef NOMPI

  USE mo_icon_cpl, ONLY : MPI_SUCCESS, MPI_COMM_NULL,         &
   &                      l_MPI_was_initialized,              &
   &                      l_debug, cplout,                    &
   &                      fieldname,                          &
   &                      grids, comps,                       &
   &                      fields, complist,                   &
   &                      nbr_active_comps,                   &
   &                      nbr_active_grids,                   &
   &                      nbr_active_fields,                  &
   &                      nbr_ICON_fields, nbr_ICON_comps,    &
   &                      nbr_ICON_grids, nbr_ICON_inc,       &
   &                      cpl_field_none,                     &
   &                      cpl_field_avg,                      &
   &                      cpl_field_acc,                      &
   &                      ICON_comm,                          &
   &                      ICON_global_rank, ICON_global_size, &
   &                      ICON_local_rank, ICON_local_size

!  USE mo_master_control,      ONLY: get_my_process_component

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER    :: version = '$Id$'

  LOGICAL                        :: l_MPI_is_initialized ! to check whether MPI_init was called.

  ! Return code for error handling

  CHARACTER(len=132)             :: err_string ! Error string for MPI
  INTEGER                        :: len        ! length 
  INTEGER                        :: ierr       ! returned error from MPI functions

#else

  IMPLICIT NONE

#endif

  PUBLIC :: icon_cpl_init_comp
!  PRIVATE

CONTAINS

  SUBROUTINE icon_cpl_init_comp ( comp_name, comp_id, ierror )

    CHARACTER(len=*), INTENT(in) :: comp_name
    INTEGER, INTENT(out)         :: comp_id
    INTEGER, INTENT(out)         :: ierror

    INTEGER                      :: comp_comm

    ! -------------------------------------------------------------------
    ! for redirecting stdout
    ! -------------------------------------------------------------------

    INTEGER                      :: lenstr
#if ! defined (__NAG)
    INTEGER                      :: istat
    INTEGER                      :: length_of_integer
    INTEGER                      :: ndibuf
    INTEGER                      :: ipos
    INTEGER, ALLOCATABLE         :: ibuf(:)
#endif

    ! TODO: replace with mo_namelist

    INTEGER, PARAMETER           :: nnml = 2
    INTEGER, PARAMETER           :: parallel_io = 1 !< switch for writing cplout
                                                    !<  - 1 write in parallel
                                                    !<  - 0 all in one

    INTEGER                      :: i               !< loop count
    INTEGER                      :: key, color


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
    ! Set up list of allowed coupling fields
    ! -------------------------------------------------------------------

    fieldname(:) = ''
    fieldname(1) = 'TEST1'
    fieldname(2) = 'TEST2'
    fieldname(3) = 'TEST3'
    fieldname(4) = 'TEST4'
    fieldname(5) = 'TEST5'
    fieldname(6) = 'TEST6'
    fieldname(7) = 'TEST7'
    fieldname(8) = 'TEST8'

    ! -------------------------------------------------------------------
    ! Derive component communicators
    ! -------------------------------------------------------------------

    color = 1!get_my_process_component() ! Rene: this should not be in here.
    key   = 0

    CALL cpl_nml_setup(comp_id)
    
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
    ! Redirect standard output and standard error output if requested
    ! -------------------------------------------------------------------

    IF ( l_redirect_stdout ) THEN

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

    DO comp_id = 1, nbr_ICON_comps
       IF ( .NOT. comps(comp_id)%l_comp_status ) EXIT
    ENDDO

    IF ( comp_id > nbr_ICON_comps ) THEN
       PRINT *, 'number of requested components exceeds maximum of nbr_ICON_comps'
       CALL MPI_Abort ( ICON_comm, 1, ierr )
    ENDIF

    comps(comp_id)%comp_name     = TRIM(comp_name)
    comps(comp_id)%l_comp_status = .TRUE.

    nbr_active_comps = nbr_active_comps + 1

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

    IF ( .NOT. ASSOCIATED(fields) ) THEN

       nbr_active_fields = 0

       nbr_ICON_fields = nbr_ICON_inc

       ALLOCATE ( fields(nbr_ICON_fields), STAT = ierr )
       IF ( ierr > 0 ) THEN
          PRINT *, ' Error allocating fields '
          CALL MPI_Abort ( ICON_COMM, 1, ierr )
       ENDIF

       fields(:)%comp_id        = -1
       fields(:)%grid_id        = -1
       fields(:)%event_id       = -1
       fields(:)%lag            =  0
       fields(:)%l_field_status = .FALSE.

       DO i = 1, nbr_ICON_fields
          NULLIFY ( fields(i)%send_field_acc )
       ENDDO

    ENDIF

    ! -------------------------------------------------------------------
    ! Initialize coupling, substitute for the OASIS4 XML reading and storage
    ! -------------------------------------------------------------------

    DO i = 1, nbr_ICON_fields
       IF ( complist(comp_id)%l_time_accumulation ) THEN
          fields(i)%coupling%time_operation = cpl_field_acc
       ELSE IF ( complist(comp_id)%l_time_average ) THEN
          fields(i)%coupling%time_operation = cpl_field_avg
       ELSE
          fields(i)%coupling%time_operation = cpl_field_none
       ENDIF

       fields(i)%coupling%coupling_freq = complist(comp_id)%coupling_freq
       fields(i)%coupling%time_step     = complist(comp_id)%time_step
    ENDDO
#else

    ! -------------------------------------------------------------------
    ! Initialise variables
    ! -------------------------------------------------------------------

    i         = 0
    ierror    = 0
    comp_id   = 0
    comp_comm = 0

#endif

    ierror = set_cpl_local_comm ( comp_comm )

  END SUBROUTINE icon_cpl_init_comp

END MODULE mo_icon_cpl_init_comp
