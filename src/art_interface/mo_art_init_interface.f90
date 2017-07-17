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
  USE mo_art_config,                    ONLY: t_art_config
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

  PUBLIC :: art_init_interface, art_calc_number_of_art_tracers

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
SUBROUTINE art_calc_number_of_art_tracers(art_config_jg)
  IMPLICIT NONE
  TYPE(t_art_config), INTENT(inout) :: art_config_jg
#ifdef __ICON_ART
  INTEGER :: number_elem
  TYPE(t_xml_file) :: tixi_file
#endif

  art_config_jg%iart_ntracer = 0

#ifdef __ICON_ART

  IF (art_config_jg%lart_chem) THEN
    CALL art_open_xml_file(TRIM(art_config_jg%cart_chemistry_xml),tixi_file)

    CALL art_get_childnumber_xml(tixi_file,"/tracers",number_elem)

    CALL art_close_xml_file(tixi_file)

    art_config_jg%iart_ntracer = art_config_jg%iart_ntracer + number_elem
  END IF

  IF (art_config_jg%lart_aerosol) THEN
    CALL art_open_xml_file(TRIM(art_config_jg%cart_aerosol_xml),tixi_file)

    CALL art_get_childnumber_xml(tixi_file,"/tracers",number_elem)

    CALL art_close_xml_file(tixi_file)

    art_config_jg%iart_ntracer = art_config_jg%iart_ntracer + number_elem
  END IF

  IF (art_config_jg%lart_passive) THEN
    CALL art_open_xml_file(TRIM(art_config_jg%cart_passive_xml),tixi_file)

    CALL art_get_childnumber_xml(tixi_file,"/tracers",number_elem)

    CALL art_close_xml_file(tixi_file)

    art_config_jg%iart_ntracer = art_config_jg%iart_ntracer + number_elem
  END IF
#endif

END SUBROUTINE art_calc_number_of_art_tracers
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_init_interface
