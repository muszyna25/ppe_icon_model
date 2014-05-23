!>
!! Finalisation of the ICON coupler
!!
!! <Description>
!!
!! This routine finalises the ICON coupling, call MPI_Finalize
!! if required and delaaocte memory.
!!
!! 
!! @author Rene Redler, MPI-M
!!
!! $Id:$
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
MODULE mo_icon_cpl_finalize

#ifndef NOMPI
  USE mpi,         ONLY : MPI_SUCCESS
  USE mo_icon_cpl, ONLY : ICON_comm_active,       &
   &                      cplout,                 &
   &                      debug_coupler_level,    &
   &                      l_MPI_was_initialized,  &
   &                      target_locs,            &
   &                      source_locs

#else

  USE mo_icon_cpl, ONLY : ICON_comm_active,       &
   &                      cplout,                 &
   &                      debug_coupler_level,    &
   &                      l_MPI_was_initialized,  &
   &                      target_locs,            &
   &                      source_locs

#endif

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER    :: version = '$Id$'

  INTEGER                        :: i
  LOGICAL                        :: l_MPI_is_finalized ! to check whether MPI_finalize was called.

  ! Return code for error handling

  CHARACTER(len=132)             :: err_string ! Error string for MPI
  INTEGER                        :: len        ! length 
  INTEGER                        :: ierr       ! returned error from MPI functions 
  INTEGER                        :: ierror     ! returned error from MPI functions 

  PUBLIC :: ICON_cpl_finalize

CONTAINS

  SUBROUTINE ICON_cpl_finalize

    IF ( debug_coupler_level > 0 ) &
       WRITE ( cplout , '(a)') 'Calling ICON_Finalize '

    ! -------------------------------------------------------------------
    ! Deallocate memory
    ! -------------------------------------------------------------------
    !
    IF ( ASSOCIATED(target_locs) ) THEN
       DO i = 1, SIZE(target_locs)
          IF ( ASSOCIATED ( target_locs(i)%target_list ) ) &
               DEALLOCATE ( target_locs(i)%target_list )
          IF ( ASSOCIATED ( target_locs(i)%source_list ) ) &
               DEALLOCATE ( target_locs(i)%source_list )
       ENDDO
       DEALLOCATE ( target_locs )
    ENDIF
    !
    IF ( ASSOCIATED(source_locs) ) THEN
       DO i = 1, SIZE(source_locs)
          IF ( ASSOCIATED ( source_locs(i)%source_list ) ) &
               DEALLOCATE ( source_locs(i)%source_list )
       ENDDO
       DEALLOCATE ( source_locs )
    ENDIF
#ifndef NOMPI

    CALL MPI_BARRIER ( ICON_comm_active, ierr )

    ! -------------------------------------------------------------------
    ! Check MPI Initialization and finalize MPI
    ! -------------------------------------------------------------------

    IF ( l_MPI_was_initialized ) THEN

       CALL MPI_Finalized ( l_MPI_is_finalized, ierr )

       IF ( .NOT. l_MPI_is_finalized ) THEN

          IF ( debug_coupler_level > 0 ) &
             WRITE ( cplout , '(a)') 'Calling MPI_Finalize '

          CALL MPI_Finalize ( ierr )

          IF ( ierr /= MPI_SUCCESS ) THEN
             CALL MPI_Error_string ( ierr, err_string, len, ierror )
             WRITE  ( * , '(a,a)' ) 'Error in MPI_Init ', err_string
          ENDIF

       ENDIF

    ENDIF

#endif

    IF ( debug_coupler_level > 0 ) CLOSE ( unit = cplout )

  END SUBROUTINE ICON_cpl_finalize

END MODULE mo_icon_cpl_finalize
