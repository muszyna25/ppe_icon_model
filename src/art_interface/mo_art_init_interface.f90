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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_init_interface

  USE mo_run_config,                    ONLY: lart
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
    timer_art_initInt
#ifdef __ICON_ART
  USE mo_art_init_all_dom,              ONLY: art_init_all_dom
  USE mo_art_clean_up,                  ONLY: art_clean_up
  USE mo_art_read_xml,                  ONLY: art_get_childnumber_xml,   &
                                          &   art_open_xml_file,         &
                                          &   art_close_xml_file,        &
                                          &   t_xml_file
#endif
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_init_interface, art_calc_number_of_art_tracers_xml

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_init_interface(n_dom,defcase)

  INTEGER,intent(in)          :: &
    &  n_dom                        !< number of model domains
  CHARACTER(LEN=*),intent(in) :: &
    &  defcase                      !< construction or destruction?
    
#ifdef __ICON_ART
  if (lart) then
    IF (timers_level > 3) CALL timer_start(timer_art_initInt)

    if (TRIM(defcase) == 'construct') then
      CALL art_init_all_dom(n_dom)
    end if
      
    if (TRIM(defcase) == 'destruct') then
      CALL art_clean_up(n_dom)
    end if

    IF (timers_level > 3) CALL timer_stop(timer_art_initInt)

  end if
#endif

END SUBROUTINE art_init_interface
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_number_of_art_tracers_xml(xml_filename,auto_ntracer)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in) :: &
     &   xml_filename                    !< name of the xml file
  INTEGER, INTENT(out) ::         &
     &   auto_ntracer                    !< automatically computed number of
                                         !   tracer within one xml file
#ifdef __ICON_ART
  TYPE(t_xml_file) :: tixi_file
#endif

  auto_ntracer = 0

#ifdef __ICON_ART
  CALL art_open_xml_file(TRIM(xml_filename),tixi_file)

  CALL art_get_childnumber_xml(tixi_file,"/tracers",auto_ntracer)

  CALL art_close_xml_file(tixi_file)
#endif

END SUBROUTINE art_calc_number_of_art_tracers_xml
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_init_interface
