!>
!! Provides interface to the ART-routine for initializing ART type structures. 
!!
!! This module provides an interface to the ART-routine art_init_all_dom
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Daniel Rieger, KIT
!!
!! @par Revision History
!! Initial revision by Daniel Rieger (2013-09-13)
!!
!! @par Copyright
!! 2002-2010 by DWD, MPI-M, and KIT.
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
!!    an according license agreement with DWD, MPI-M, and KIT.
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
MODULE mo_art_init_interface

#ifdef __ICON_ART
  USE mo_art_init_all_dom,        ONLY: art_init_all_dom
  USE mo_art_clean_up,            ONLY: art_clean_up
#endif
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_init_interface

CONTAINS

  SUBROUTINE art_init_interface(n_dom,defcase)

    INTEGER,intent(in) :: n_dom     !< number of model domains
    CHARACTER(LEN=*),intent(in) :: defcase  !< construction or destruction?
    
#ifdef __ICON_ART
    
    
    if (TRIM(defcase) == 'construct') then
      CALL art_init_all_dom(n_dom)
    end if
    
    if (TRIM(defcase) == 'destruct') then
      CALL art_clean_up(n_dom)
    end if
    
#endif

  END SUBROUTINE art_init_interface

END MODULE mo_art_init_interface
