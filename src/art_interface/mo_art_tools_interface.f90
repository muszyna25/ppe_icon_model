!>
!! Provides interface to the ART-routine for using tools
!!
!! This module provides an interface to the ART-routine art_ini_tracer.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Daniel Rieger, KIT
!!
!! @par Revision History
!! Initial revision by Daniel Rieger (2013-12-16)
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

MODULE mo_art_tools_interface
  USE mo_nonhydro_types,                ONLY: t_nh_state
#ifdef __ICON_ART
  USE mo_art_unit_conversion,           ONLY: art_massmix2density
#endif

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC  :: art_tools_interface 

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_tools_interface(defcase,p_nh_state,jg)
  !>
  !! Interface for ART tools
  !!
  !! @par Revision History
  !! Initial revision by Daniel Rieger, KIT (2013-12-16)
  
  CHARACTER(len=*),INTENT(in)        :: & 
    &  defcase                            !< definition of case 
  TYPE(t_nh_state),TARGET,INTENT(in) :: &
    &  p_nh_state                         !< prognostic state
  INTEGER,INTENT(in)                 :: &
    &  jg                                 !< domain index
  
#ifdef __ICON_ART

  IF (TRIM(defcase) .EQ. 'unit_conversion') THEN
    CALL art_massmix2density(p_nh_state,jg)
  END IF

#endif
    
END SUBROUTINE art_tools_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_tools_interface
