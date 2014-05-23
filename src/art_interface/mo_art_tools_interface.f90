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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
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
