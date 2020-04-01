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

  USE mo_kind,                          ONLY: wp
  USE mo_run_config,                    ONLY: lart, ntracer
  USE mo_exception,                     ONLY: finish
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art_initInt
  USE mo_storage,                       ONLY: t_storage
  USE mo_linked_list,                   ONLY: t_var_list
  USE mo_nonhydro_types,                ONLY: t_nh_prog, t_nh_state
  USE mo_ext_data_types,                ONLY: t_external_data
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag

  USE mo_art_config,                    ONLY: ctracer_art
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH

  USE mtime,                            ONLY: datetime
#ifdef __ICON_ART
  USE mo_art_collect_atmo_state,        ONLY: art_collect_atmo_state_nwp,   &
                                          &   art_update_atmo_state_nwp,    &
                                          &   art_collect_atmo_state_echam, &
                                          &   art_update_atmo_state_echam,  &
                                          &   art_init_tracer_values_nwp,   &
                                          &   art_init_tracer_values_echam

  USE mo_art_init_all_dom,              ONLY: art_init_all_dom
  USE mo_art_clean_up,                  ONLY: art_clean_up
  USE mo_art_tagging,                   ONLY: get_number_tagged_tracer
  USE mo_art_read_xml,                  ONLY: art_get_childnumber_xml,   &
                                          &   art_read_elements_xml,     &
                                          &   art_open_xml_file,         &
                                          &   art_close_xml_file,        &
                                          &   t_xml_file
#endif
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_init_interface, art_calc_ntracer_and_names
  PUBLIC :: art_init_atmo_tracers_nwp, art_init_atmo_tracers_echam
  PUBLIC :: art_update_atmo_phy

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
  IF (lart) THEN
    IF (timers_level > 3) CALL timer_start(timer_art_initInt)

    IF (TRIM(defcase) == 'construct') THEN
      CALL art_write_vcs_info
      CALL art_init_all_dom(n_dom)
    END IF
      
    IF (TRIM(defcase) == 'destruct') THEN
      CALL art_clean_up(n_dom)
    END IF

    IF (timers_level > 3) CALL timer_stop(timer_art_initInt)

   END IF 
#endif

END SUBROUTINE art_init_interface
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_number_of_art_tracers_xml(xml_filename,auto_ntracer, tracer_names)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in) :: &
     &   xml_filename                    !< name of the xml file
  INTEGER, INTENT(out) ::         &
     &   auto_ntracer                    !< automatically computed number of
                                         !   tracer within one xml file
  CHARACTER(LEN = *), POINTER, INTENT(inout) :: &
     &   tracer_names(:)

  ! Local variables
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
    ALLOCATE(tracer_names(auto_ntracer))

    DO idx_tracer = 1,auto_ntracer
      ! Create a storage container
      CALL storage%init(lcase_sensitivity=.FALSE.)

      WRITE(idx_tracer_str,'(I5)') idx_tracer
      CALL art_read_elements_xml(tixi_file,'/tracers/*['     &
               &               //TRIM(ADJUSTL(idx_tracer_str))//']/',storage)

      CALL storage%get('name',tracer_names(idx_tracer))

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
! branch etc. used in the .out file
! Part of Module: mo_art_init_interface.f90
! Author: Jennifer Schroeter, KIT
! Initial Release: 2018-01-18                
! Modifications:
!>

  USE mo_art_util_vcs,   ONLY: art_util_repository_url, art_util_branch_name, &
                           &   art_util_revision_key
  USE mo_exception,      ONLY: message_text, message, finish
  USE mo_mpi,            ONLY: my_process_is_global_root
                       

  CHARACTER(len=256) :: art_repository  = ''
  CHARACTER(len=256) :: art_branch      = ''
  CHARACTER(len=256) :: art_revision    = ''
  

  INTEGER :: nlen
  nlen = 256
  CALL art_util_repository_url(art_repository, nlen)
  nlen = 256
  CALL art_util_branch_name(art_branch, nlen)
  nlen = 256
  CALL art_util_revision_key(art_revision, nlen)

  IF (my_process_is_global_root()) THEN

    WRITE(message_text,'(a,a)') 'ART Repository: ', TRIM(art_repository)
      CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'ART Branch    : ', TRIM(art_branch)
      CALL message('',message_text)
    WRITE(message_text,'(a,a)') 'ART Revision  : ', TRIM(art_revision)
      CALL message('',message_text)

  END IF
#endif
  
END SUBROUTINE art_write_vcs_info
!!
!!-------------------------------------------------------------------------
!!

SUBROUTINE art_calc_ntracer_and_names(auto_ntracer,                           &
               &      cart_chemistry_xml, cart_aerosol_xml, cart_passive_xml, &
               &      lart_chem, lart_aerosol, lart_passive)
  INTEGER, INTENT(inout) :: auto_ntracer      !< automatically computed number of tracers
  CHARACTER(LEN = *), INTENT(in) :: &
       &      cart_chemistry_xml,   &
       &      cart_aerosol_xml,     &
       &      cart_passive_xml
  LOGICAL, INTENT(in) ::            &
       &      lart_chem,            &
       &      lart_aerosol,         &
       &      lart_passive

  ! Local variables
  INTEGER  ::                       &
       &   auto_ntracer_chemistry,  &
       &   auto_ntracer_aerosol,    &
       &   auto_ntracer_passive          !< art ntracer from one xml file
  CHARACTER(LEN = MAX_CHAR_LENGTH), POINTER ::   &
       &   tracer_names(:),           &
       &   tracer_names_chemistry(:), &
       &   tracer_names_passive(:),   &
       &   tracer_names_aerosol(:)
  LOGICAL ::                 &
       &   l_exist

  auto_ntracer = 0
  ! chemistry xml file
  IF (lart_chem) THEN
    IF (TRIM(cart_chemistry_xml) == '') THEN
      CALL finish('mo_art_init_interface:art_calc_ntracer_and_names',  &
             &    'namelist parameter cart_chemistry_xml'              &
             &  //' has to be given for lart_chem == .TRUE.')
    ELSE
      INQUIRE(file = TRIM(cart_chemistry_xml), EXIST = l_exist)
  
      IF (l_exist) THEN
        CALL art_calc_number_of_art_tracers_xml(TRIM(cart_chemistry_xml),  &
                               &                auto_ntracer_chemistry, tracer_names_chemistry)
        auto_ntracer = auto_ntracer + auto_ntracer_chemistry
      ELSE
        CALL finish('mo_art_init_interface:art_calc_ntracer_and_names',  &
                    TRIM(cart_chemistry_xml)//  &
                    & ' could not be found. Check cart_chemistry_xml.')
      END IF
    END IF
  ELSE
    auto_ntracer_chemistry = 0
    NULLIFY(tracer_names_chemistry)
  END IF
  
  ! aerosol xml file
  IF (lart_aerosol) THEN
    IF (TRIM(cart_aerosol_xml) == '') THEN
      CALL finish('mo_art_init_interface:art_calc_ntracer_and_names',  &
              &   'namelist parameter cart_aerosol_xml'                &
              & //' has to be given for lart_aerosol == .TRUE.')
    ELSE
      INQUIRE(file = TRIM(cart_aerosol_xml), EXIST = l_exist)
  
      IF (l_exist) THEN
        CALL art_calc_number_of_art_tracers_xml(TRIM(cart_aerosol_xml),  &
                               &                auto_ntracer_aerosol, tracer_names_aerosol)
        auto_ntracer = auto_ntracer + auto_ntracer_aerosol
      ELSE
        CALL finish('mo_art_init_interface:art_calc_ntracer_and_names',  &
                    TRIM(cart_aerosol_xml)//  &
                    & ' could not be found. Check cart_aerosol_xml.')
      END IF
    END IF
  ELSE
    auto_ntracer_aerosol = 0
    NULLIFY(tracer_names_aerosol)
  END IF
  
  ! passive xml file
  IF (lart_passive) THEN
    IF (TRIM(cart_passive_xml) == '') THEN
      CALL finish('mo_art_init_interface:art_calc_ntracer_and_names',  &
             &    'namelist parameter cart_passive_xml'                &
             &  //' has to be given for lart_passive == .TRUE.')
    ELSE
      INQUIRE(file = TRIM(cart_passive_xml), EXIST = l_exist)
  
      IF (l_exist) THEN
        CALL art_calc_number_of_art_tracers_xml(TRIM(cart_passive_xml),  &
                               &                auto_ntracer_passive, tracer_names_passive)
        auto_ntracer = auto_ntracer + auto_ntracer_passive
      ELSE
        CALL finish('mo_art_init_interface:art_calc_ntracer_and_names',  &
                    TRIM(cart_passive_xml)//  &
                    & ' could not be found. Check cart_passive_xml.')
      END IF
    END IF
  ELSE
    auto_ntracer_passive = 0
    NULLIFY(tracer_names_passive)
  END IF


  IF (auto_ntracer > 0) THEN
    ALLOCATE(tracer_names(auto_ntracer))

    ! This setting of the tracer names actually requires the correct order of
    ! including the ART tracers in mo_art_tracer: first aerosols, then chemical,
    ! then passive tracers
    IF (auto_ntracer_aerosol > 0) THEN
      tracer_names(1:auto_ntracer_aerosol) = tracer_names_aerosol
    END IF
  
    IF (auto_ntracer_chemistry > 0) THEN
      tracer_names(auto_ntracer_aerosol+1:auto_ntracer_aerosol + auto_ntracer_chemistry) &
         &    = tracer_names_chemistry
    END IF
  
    IF (auto_ntracer_passive > 0) THEN
      tracer_names(auto_ntracer_aerosol+auto_ntracer_chemistry+1:) = tracer_names_passive
    END IF
  
    ALLOCATE(ctracer_art(ntracer+auto_ntracer))
    ctracer_art(ntracer+1:) = tracer_names
  
    DEALLOCATE(tracer_names)
    IF (ASSOCIATED(tracer_names_passive)) DEALLOCATE(tracer_names_passive)
    IF (ASSOCIATED(tracer_names_chemistry)) DEALLOCATE(tracer_names_chemistry)
    IF (ASSOCIATED(tracer_names_aerosol)) DEALLOCATE(tracer_names_aerosol)
  ELSE
    ALLOCATE(ctracer_art(ntracer))
  END IF

END SUBROUTINE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_init_atmo_tracers_nwp(jg, mtime_current, p_nh_state, ext_data, &
                 &                   prm_diag, p_prog, tracer, p_prog_list)
  IMPLICIT NONE
  INTEGER, INTENT(in) ::  &
    &  jg                   !< patch id
  TYPE(datetime), POINTER :: &
    &  mtime_current        !< current model date
  TYPE(t_nh_state), INTENT(in) :: &
    &  p_nh_state           !< state variables of ICON
  TYPE(t_external_data), INTENT(in) :: &
    &  ext_data             !< external fields for NWP physics (boundary fields etc.)
  TYPE(t_nwp_phy_diag), INTENT(in) :: &
    &  prm_diag             !< physics fields for NWP physics
  TYPE(t_nh_prog), INTENT(in) :: &
    &  p_prog               !< prognostic variables of ICON
  REAL(wp), POINTER :: &
    &  tracer(:,:,:,:)      !< tracer values (ICON kg/kg, Aerosols ??, chemistry  mol/mol, passive none)
  TYPE(t_var_list), INTENT(in) :: &
    &  p_prog_list          !< list of prognostic variables

#ifdef __ICON_ART
  IF (lart) THEN
    CALL art_collect_atmo_state_nwp(jg, mtime_current, p_nh_state,  &
                   &                ext_data, prm_diag, p_prog)

    CALL art_init_tracer_values_nwp(jg, tracer, mtime_current, p_prog_list)
  END IF
#endif

END SUBROUTINE art_init_atmo_tracers_nwp
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_init_atmo_tracers_echam(jg, mtime_current, p_nh_state, &
                 &                     p_prog, tracer, p_prog_list)
  IMPLICIT NONE
  INTEGER, INTENT(in) ::  &
    &  jg                   !< patch id
  TYPE(datetime), POINTER :: &
    &  mtime_current        !< current model date
  TYPE(t_nh_state), INTENT(in) :: &
    &  p_nh_state           !< state variables of ICON
  TYPE(t_nh_prog), INTENT(in) :: &
    &  p_prog               !< prognostic variables of ICON
  REAL(wp), POINTER  :: &
    &  tracer(:,:,:,:)      !< tracer values (ICON kg/kg, Aerosols ??, chemistry  mol/mol, passive none)
  TYPE(t_var_list), INTENT(in) :: &
    &  p_prog_list          !< list of prognostic variables

#ifdef __ICON_ART
  IF (lart) THEN
    CALL art_collect_atmo_state_echam(jg, mtime_current, p_nh_state, p_prog)

    CALL art_init_tracer_values_echam(jg, tracer, mtime_current, p_prog_list)
  END IF
#endif

END SUBROUTINE art_init_atmo_tracers_echam
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_update_atmo_phy(jg, mtime_current, p_prog, prm_diag)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: &
    &  jg                   !< patch id
  TYPE(datetime), POINTER :: &
    &  mtime_current        !< current model date
  TYPE(t_nh_prog), INTENT(in) :: &
    &  p_prog               !< prognostic variables of ICON
  TYPE(t_nwp_phy_diag), INTENT(in), OPTIONAL :: &
    &  prm_diag             !< phyics fields for NWP physics


#ifdef __ICON_ART
  IF (lart) THEN
    IF (PRESENT(prm_diag)) THEN
      CALL art_update_atmo_state_nwp(jg,mtime_current, p_prog, prm_diag)
    ELSE
      CALL art_update_atmo_state_echam(jg,mtime_current, p_prog)
    END IF
  END IF
#endif
END SUBROUTINE art_update_atmo_phy

END MODULE mo_art_init_interface
