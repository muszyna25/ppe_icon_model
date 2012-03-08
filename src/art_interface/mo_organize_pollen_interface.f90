!>
!! Provides interface to ART-routine organize_pollen
!!
!! This module provides an interface to the ART-routine organize_pollen.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2012-01-27)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_organize_pollen_interface

USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH
#ifdef __ICON_ART
!DR  USE mo_art,               ONLY: organize_pollen
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: organize_pollen_ifc


CONTAINS


  !>
  !! Interface for ART-routine organize_pollen
  !!
  !! This interface calls the ART-routine organize_pollen, if ICON has been 
  !! built including the ART-package. Otherwise, this is simply a dummy 
  !! routine.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-01-27)
  !!
  SUBROUTINE organize_pollen_ifc()

    ! CHARACTER (LEN=MAX_CHAR_LENGTH) :: yaction
    ! CHARACTER (LEN=10)              :: ydate_ini
    ! INTEGER                         :: ierror
    ! CHARACTER (LEN=MAX_CHAR_LENGTH) :: yerrmsg

    !-----------------------------------------------------------------------

#ifdef __ICON_ART
!DR    CALL organize_pollen (yaction, ydate_ini, ierror, yerrmsg)
#endif


  END SUBROUTINE organize_pollen_ifc


END MODULE mo_organize_pollen_interface

