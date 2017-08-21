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
  USE mo_kind,                          ONLY: wp
  USE mo_run_config,                    ONLY: lart
  USE mo_linked_list,                   ONLY: t_var_list
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_toolInt
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
SUBROUTINE art_tools_interface(defcase, prog_list, tracer_now, tracer_new, rho)
  !>
  !! Interface for ART tools
  !!
  !! @par Revision History
  !! Initial revision by Daniel Rieger, KIT (2013-12-16)
  
  CHARACTER(len=*),INTENT(in)        :: & 
    &  defcase                            !< definition of case 
  TYPE(t_var_list),TARGET,INTENT(in) :: &
    &  prog_list                          !< prognostic state list
  REAL(wp),INTENT(in)                :: &
    &  tracer_now(:,:,:,:),             & !< tracer concentrations at timelevel nnow
    &  rho(:,:,:)                         !< Density
  REAL(wp),INTENT(inout)             :: &
    &  tracer_new(:,:,:,:)                !< tracer concentrations at timelevel nnew
  
#ifdef __ICON_ART
  IF (lart) THEN
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_toolInt)

    IF (TRIM(defcase) .EQ. 'unit_conversion') THEN
      CALL art_massmix2density(prog_list, tracer_now, tracer_new, rho)
    ENDIF

    IF (timers_level > 3) CALL timer_stop(timer_art_toolInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  ENDIF
#endif
    
END SUBROUTINE art_tools_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_tools_interface
