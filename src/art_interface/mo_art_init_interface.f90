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
                                          &   timer_art_initInt
  USE mo_storage,                       ONLY: t_storage
#ifdef __ICON_ART
  USE mo_art_init_all_dom,              ONLY: art_init_all_dom
  USE mo_art_clean_up,                  ONLY: art_clean_up
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN
  USE mo_art_tagging,                   ONLY: get_number_tagged_tracer
  USE mo_art_read_xml,                  ONLY: art_get_childnumber_xml,   &
                                          &   art_read_elements_xml,     &
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
      CALL art_write_vcs_info
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
  INTEGER ::          &
     &   idx_tracer,  &                  !< index of the tracer in XML
     &   ntags                           !< number of tags for the current tracer
  CHARACTER(LEN = 5) :: &
     &   idx_tracer_str                  !< string of the index
  TYPE(t_storage) :: &
     &   storage                         !< temporally created storage for the tracer

#ifdef __ICON_ART
  TYPE(t_xml_file) :: tixi_file          !< tracer XML file
#endif

  auto_ntracer = 0

#ifdef __ICON_ART
  CALL art_open_xml_file(TRIM(xml_filename),tixi_file)

  CALL art_get_childnumber_xml(tixi_file,"/tracers",auto_ntracer)

  IF (auto_ntracer > 0) THEN
    DO idx_tracer = 1,auto_ntracer
      ! Create a storage container
      CALL storage%init(lcase_sensitivity=.FALSE.)

      WRITE(idx_tracer_str,'(I5)') idx_tracer
      CALL art_read_elements_xml(tixi_file,'/tracers/*['     &
               &               //TRIM(ADJUSTL(idx_tracer_str))//']/',storage)

      ntags = get_number_tagged_tracer(storage)

      IF (ntags > 1) THEN
        ! the first tag is already included in the calculation of the number of
        ! tracer so add only the tagged tracer with number greater than 1
        ! (tag002 to tag999)
        auto_ntracer = auto_ntracer + ntags - 1
      END IF

      ! Set metadata storage free again
      CALL storage%destruct

    END DO
  END IF
  

  CALL art_close_xml_file(tixi_file)
#endif

END SUBROUTINE art_calc_number_of_art_tracers_xml
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_write_vcs_info
#ifdef __ICON_ART

!<
! SUBROUTINE art_write_vcs_info                   
! This subroutine writes information about repository
! branch etc. used in the .out fie
! Part of Module: mo_art_init_interface.f90
! Author: Jennifer Schroeter, KIT
! Initial Release: 2018-01-18                
! Modifications:
!>

  USE mo_art_util_vcs,     ONLY: art_util_repository_url, art_util_branch_name, &
                           & art_util_revision_key
  USE mo_exception,    ONLY: message_text, message, finish
  USE mo_mpi,          ONLY: my_process_is_global_root
                       

  CHARACTER(len=256) :: art_repository  = ''
  CHARACTER(len=256) :: art_branch      = ''
  CHARACTER(len=256) :: art_revision    = ''
  

  INTEGER :: nlen
  nlen = 256
  call art_util_repository_url(art_repository, nlen)
  nlen = 256
  call art_util_branch_name(art_branch, nlen)
  nlen = 256
  call art_util_revision_key(art_revision, nlen)
  IF (my_process_is_global_root()) THEN


  WRITE(message_text,'(a,a)') 'ART Repository: ', TRIM(art_repository)
    CALL message('',message_text)
  WRITE(message_text,'(a,a)') 'ART Branch    : ', TRIM(art_branch)
    CALL message('',message_text)
  WRITE(message_text,'(a,a)') 'ART Revision  : ', TRIM(art_revision)
    CALL message('',message_text)
  ENDIF
#endif
  
END SUBROUTINE art_write_vcs_info


END MODULE mo_art_init_interface
