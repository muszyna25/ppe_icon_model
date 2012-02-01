!>
!! Initialisation of the ICON coupler
!!
!! <Description>
!! 
!! This routine returns a grid_id and an MPI communicator
!! The MPI communicator can be used fo local communication
!! inside the model component.
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
!! This software is provided for non-commercial use only.
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

MODULE mo_icon_cpl_init

  USE mo_kind, ONLY     : wp, dp

#ifndef NOMPI

  USE mpi

  USE mo_icon_cpl, ONLY : comps, nbr_ICON_comps, datatype,     &
   &                      PRISM_CHARACTER, PRISM_INTEGER,      &
   &                      PRISM_LOGICAL, PRISM_REAL,           &
   &                      PRISM_DOUBLE_PRECISION,              &
   &                      PRISM_COMPLEX, PRISM_DOUBLE_COMPLEX, &
   &                      l_MPI_was_initialized,               &
   &                      debug_coupler_level,                 &
   &                      cplout, maxchar,                     &
   &                      ICON_comm,                           &
   &                      ICON_global_rank, ICON_global_size,  &
   &                      nbr_ICON_comps,                      &
   &                      initial_date

  USE mo_time_base,     ONLY : set_julian_day
  USE mo_event_manager, ONLY : event_init

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER    :: version = '$Id$'

  LOGICAL                        :: l_MPI_is_initialized ! to check whether MPI_init was called.

  CHARACTER(len=maxchar)         :: filename
  LOGICAL                        :: unit_is_occupied
  INTEGER                        :: i          ! loop count
  ! Return code for error handling

  CHARACTER(len=132)             :: err_string ! Error string for MPI
  INTEGER                        :: len        ! length 
  INTEGER                        :: ierr       ! returned error from MPI functions
  INTEGER                        :: ierror     ! returned error from MPI functions

#else
  USE mo_icon_cpl, ONLY : comps, nbr_ICON_comps

  IMPLICIT NONE

  INTEGER                        :: ierror     ! returned error from MPI functions

#endif

  PUBLIC :: ICON_cpl_init

CONTAINS

  SUBROUTINE ICON_cpl_init(debug_level)

     INTEGER, OPTIONAL, INTENT(in) :: debug_level

    ! -------------------------------------------------------------------
    ! Initialise variables
    ! -------------------------------------------------------------------

    ALLOCATE ( comps(nbr_ICON_comps) )

    comps(:)%l_comp_status = .FALSE.

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
    ! Determine datatype for physical fields
    ! -------------------------------------------------------------------

    IF ( wp == dp ) THEN
       datatype = MPI_DOUBLE_PRECISION
    ELSE
       datatype = MPI_REAL
    ENDIF

    ! -------------------------------------------------------------------
    ! Create own communicator for coupling
    ! -------------------------------------------------------------------

    CALL MPI_Comm_dup (MPI_COMM_WORLD, ICON_comm, ierr)

    IF ( ierr /= MPI_SUCCESS ) THEN
       CALL MPI_Error_string ( ierr, err_string, len, ierror )
       WRITE  ( * , '(A,A)' ) 'Error in MPI_Comm_dub ', err_string
    ENDIF

    ! -------------------------------------------------------------------
    ! Size of coupled application and rank in the coupled system
    ! -------------------------------------------------------------------

    CALL MPI_Comm_rank ( ICON_comm, ICON_global_rank, ierr)
    IF ( ierr /= MPI_SUCCESS ) THEN
       CALL MPI_Error_string ( ierr, err_string, len, ierror )
       WRITE  ( * , '(A,A)' ) 'Error in getting global rank ', err_string
    ENDIF

    CALL MPI_Comm_size ( ICON_comm, ICON_global_size, ierr)
    IF ( ierr /= MPI_SUCCESS ) THEN
       CALL MPI_Error_string ( ierr, err_string, len, ierror )
       WRITE  ( * , '(A,A)' ) 'Error in getting global rank ', err_string
    ENDIF

    ! -------------------------------------------------------------------
    ! Preparation for couler debug output
    ! -------------------------------------------------------------------

    debug_coupler_level=0

    IF (PRESENT(debug_level)) debug_coupler_level = debug_level
    
    IF ( debug_coupler_level > 0 ) THEN

       ! Find a free unit

       DO i = 10, 100
          INQUIRE (UNIT=i, OPENED=unit_is_occupied)
          IF ( .NOT. unit_is_occupied ) EXIT
       ENDDO

       IF ( i > 100 ) THEN
          WRITE ( * , '(A)' ) 'WARNING: No free unit for coupler debug output!'
          WRITE ( * , '(A)' ) '         Debug output is turned off.'
          debug_coupler_level = 0
       ELSE
          cplout = i
       ENDIF

    ENDIF

    IF ( debug_coupler_level > 0 ) THEN
       WRITE ( filename, '(a5,I4.4,A4)' )  'ICON_', ICON_global_rank, '.log'
       OPEN  ( unit = cplout, file = filename, status = 'unknown', form = 'formatted' )
    ENDIF

    ! -------------------------------------------------------------------
    ! Initialise datatypes for psmile_bsend_init
    ! -------------------------------------------------------------------

    CALL init_bsend_datatypes ( ierr )

    ! -------------------------------------------------------------------
    ! Initialise list of events (move to some other place later)
    ! -------------------------------------------------------------------

    ! From somewhere we need to get the time_config%ini_datetime
    ! t_julian_date format.

    call set_julian_day ( 0, 1, 1, 0, initial_date )

    CALL event_init

#else

    ! -------------------------------------------------------------------
    ! Initialise variables
    ! -------------------------------------------------------------------

    ierror = 0

    IF (PRESENT(debug_level)) RETURN

#endif


  END SUBROUTINE ICON_cpl_init
  
  !! @par License
  !!
  !!-----------------------------------------------------------------------
  !! Copyright 2010, NEC Europe Ltd., Germany.
  !! All rights reserved. Use is subject to OASIS4 license terms.
  !!-----------------------------------------------------------------------
  
  SUBROUTINE init_bsend_datatypes ( ierror )
    !
    IMPLICIT NONE
    !
    ! !OUTPUT PARAMETERS:
    !
    INTEGER, INTENT (Out)               :: ierror
    !
    !     Returns the error code of init_bsend_datatypes;
    !             ierror = 0 : No error
    !             ierror > 0 : Severe error
    !
    ! !LOCAL VARIABLES
    !
    INTEGER            :: i
    !
    !     ... Lengths of Fortran MPI datatypes
    !
    INTEGER, PARAMETER :: noDatatypes = 9
    INTEGER            :: datatypes2mpi    (noDatatypes)
    INTEGER            :: lengths_of_types (noDatatypes)
    INTEGER            :: length_of_integer
    !
    ! !DESCRIPTION:
    !
    ! Copy of psmile_init_datatypes from the OASIS4 distribution.
    ! Subroutine init_bsend_datatypes initializes the handling of
    ! predefined (MPI) datatypes.
    !
    ! !REVISION HISTORY:
    !
    !   Date      Programmer   Description
    ! ----------  ----------   -----------
    ! 03.06.20    H. Ritzdorf  created
    ! 10.02.13    R. Redler    modified for use within ICON
    !
    !EOP
    !----------------------------------------------------------------------
    !
    ! Initialisation
    !
    ierror = 0
    !
    ! Transfer-table from PRISM datatypes into MPI datatypes
    !

#ifndef NOMPI

    datatypes2mpi (1:noDatatypes)          = MPI_DATATYPE_NULL
    !
    datatypes2mpi (PRISM_Character)        = MPI_CHARACTER
    datatypes2mpi (PRISM_Integer)          = MPI_INTEGER
    datatypes2mpi (PRISM_LogicaL)          = MPI_LOGICAL
    datatypes2mpi (PRISM_Real)             = MPI_REAL
    datatypes2mpi (PRISM_Double_precision) = MPI_DOUBLE_PRECISION
    datatypes2mpi (PRISM_Complex)          = MPI_COMPLEX
    datatypes2mpi (PRISM_Double_complex)   = MPI_DOUBLE_COMPLEX
    !
    ! Get size of integer`s in bytes
    !
    length_of_integer = BIT_SIZE (datatypes2mpi (1)) / 8
    !
    ! Initialize datatypes for PSMILe_bsend
    !
    DO i = 1, noDatatypes
       IF (datatypes2mpi (i) /= MPI_DATATYPE_NULL) THEN
          CALL MPI_Type_size (datatypes2mpi (i), lengths_of_types (i), ierror)
          !
          !           if (ierror /= MPI_SUCCESS) then
          !           end if
       ELSE
          lengths_of_types (i) = 0
       ENDIF
    ENDDO
    !
    CALL psmile_bsend_init ( datatypes2mpi, lengths_of_types, &
         noDatatypes, ierror )

#else

    length_of_integer = 0
    datatypes2mpi     = 0
    lengths_of_types  = 0
    i      = 0
    ierror = 0

#endif

  END SUBROUTINE init_bsend_datatypes

END MODULE mo_icon_cpl_init
