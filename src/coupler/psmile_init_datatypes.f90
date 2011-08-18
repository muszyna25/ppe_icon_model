!-----------------------------------------------------------------------
! Copyright 2010, NEC Europe Ltd., Germany.
! All rights reserved. Use is subject to OASIS4 license terms.
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: PSMILe_Init_datatypes
!
! !INTERFACE:

SUBROUTINE psmile_init_datatypes ( ierror )
  !
  ! !USES:
  !
#ifndef NOMPI
  USE mpi
#endif
  USE mo_icon_cpl, ONLY : PRISM_CHARACTER, PRISM_INTEGER,     &
   &                      PRISM_LOGICAL, PRISM_REAL,          &
   &                      PRISM_DOUBLE_PRECISION,             &
   &                      PRISM_COMPLEX, PRISM_DOUBLE_COMPLEX

  IMPLICIT NONE
  !
  ! !OUTPUT PARAMETERS:
  !
  INTEGER, INTENT (Out)               :: ierror
  !
  !     Returns the error code of PSMILe_Init_datatypes;
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
  ! Taken from the OASIS4 distribution.
  ! Subroutine PSMILe_Init_datatypes initializes the handling of
  ! predefined (MPI) datatypes.
  !
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
  ! $Id:$
  ! $Author:$
  !
  !----------------------------------------------------------------------
  !
  !===> Initialisation
  !
  ierror = 0
  !
  !===> Transfer-table from PRISM datatypes into MPI datatypes
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
  !===> Get size of integer`s in bytes
  !
  length_of_integer = BIT_SIZE (datatypes2mpi (1)) / 8
  !
  !===> Initialize datatypes for PSMILe_bsend
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

END SUBROUTINE PSMILe_Init_datatypes
