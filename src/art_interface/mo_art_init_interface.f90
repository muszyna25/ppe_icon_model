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

  USE mo_kind,                          ONLY: wp, i8
  USE mo_run_config,                    ONLY: lart, ntracer
  USE mo_exception,                     ONLY: finish
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art_initInt
  USE mo_time_config,                   ONLY: time_config
  USE mo_grid_config,                   ONLY: start_time
  USE mo_master_config,                 ONLY: isRestart
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_linked_list,                   ONLY: t_var_list
  USE mo_nonhydro_types,                ONLY: t_nh_prog, t_nh_state
  USE mo_ext_data_types,                ONLY: t_external_data
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag

  USE mo_art_config,                    ONLY: art_config, ctracer_art
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH

  USE mtime,                            ONLY: datetime, timedelta,           &
                                          &   MAX_TIMEDELTA_STR_LEN,         &
                                          &   newTimedelta,                  &
                                          &   getTotalMilliSecondsTimeDelta, &
                                          &   getPTStringFromMS,             &
                                          &   deallocateTimedelta
#ifdef __ICON_ART
  USE mo_art_collect_atmo_state,        ONLY: art_collect_atmo_state_nwp,   &
                                          &   art_update_atmo_state_nwp,    &
                                          &   art_collect_atmo_state_echam, &
                                          &   art_update_atmo_state_echam,  &
                                          &   art_init_tracer_values_nwp,   &
                                          &   art_init_tracer_values_echam

  USE mo_art_init_all_dom,              ONLY: art_init_all_dom
  USE mo_art_init,                      ONLY: art_init
  USE mo_art_clean_up,                  ONLY: art_clean_up
  USE mo_art_tagging,                   ONLY: get_number_tagged_tracer
  USE mo_art_read_xml,                  ONLY: art_check_tracer_children_xml,  &
                                          &   art_read_elements_xml,          &
                                          &   art_open_xml_file,              &
                                          &   art_close_xml_file,             &
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
SUBROUTINE art_calc_number_of_art_tracers_xml(xml_filename,auto_ntracer,   &
                          &                   tracer_element, tracer_names)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in) :: &
     &   xml_filename                    !< name of the xml file
  INTEGER, INTENT(out) ::         &
     &   auto_ntracer                    !< automatically computed number of
                                         !   tracer within one xml file
  CHARACTER(LEN = *), INTENT(in) :: &
    &    tracer_element
  CHARACTER(LEN = *), POINTER, INTENT(inout) :: &
     &   tracer_names(:)                 !< names of the tracers in XML

  ! Local variables
  INTEGER ::          &
     &   idx_tracer,  &                  !< index of the tracer in XML
     &   ntags                           !< number of tags for the current tracer
  CHARACTER(LEN = 5) :: &
     &   idx_tracer_str                  !< string of the index
  TYPE(t_key_value_store) :: &
     &   storage                         !< temporally created storage for the tracer
  CHARACTER(:), ALLOCATABLE :: &
     &   tracer_name

#ifdef __ICON_ART
  TYPE(t_xml_file) :: tixi_file          !< tracer XML file
  
#endif

  auto_ntracer = 0

#ifdef __ICON_ART
  CALL art_open_xml_file(TRIM(xml_filename),tixi_file)

  CALL art_check_tracer_children_xml(tixi_file,"/tracers",TRIM(tracer_element), auto_ntracer)

  IF (auto_ntracer > 0) THEN
    ALLOCATE(tracer_names(auto_ntracer))

    DO idx_tracer = 1,auto_ntracer
      ! Create a storage container
      CALL storage%init(.FALSE.)

      WRITE(idx_tracer_str,'(I5)') idx_tracer
      CALL art_read_elements_xml(tixi_file,'/tracers/*['     &
               &               //TRIM(ADJUSTL(idx_tracer_str))//']/',storage)

      CALL storage%get('name',tracer_name)
      tracer_names(idx_tracer) = TRIM(tracer_name)

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
  USE mo_exception,      ONLY: message_text, message
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

SUBROUTINE art_calc_ntracer_and_names()

  ! Local variables
  CHARACTER(LEN = MAX_CHAR_LENGTH) ::             &
       &      cart_chemtracer_xml,  &
       &      cart_mecca_xml,       &
       &      cart_aerosol_xml
  LOGICAL ::                        &
       &      lart_chem,            &
       &      lart_chemtracer,      &
       &      lart_mecca,           &
       &      lart_aerosol

  INTEGER  ::                       &
       &   auto_ntracer,            &
       &   auto_ntracer_chemtracer, &
       &   auto_ntracer_mecca,      &
       &   auto_ntracer_aerosol        !< art ntracer from one xml file
  CHARACTER(LEN = MAX_CHAR_LENGTH), POINTER ::   &
       &   tracer_names(:),            &
       &   tracer_names_chemtracer(:), &
       &   tracer_names_mecca(:),      &
       &   tracer_names_aerosol(:)
  LOGICAL ::                 &
       &   l_exist

  
  cart_chemtracer_xml = art_config(1)%cart_chemtracer_xml
  cart_mecca_xml      = art_config(1)%cart_mecca_xml
  cart_aerosol_xml    = art_config(1)%cart_aerosol_xml

  lart_chem          = art_config(1)%lart_chem
  lart_aerosol       = art_config(1)%lart_aerosol
  lart_chemtracer    = art_config(1)%lart_chemtracer
  lart_mecca         = art_config(1)%lart_mecca

  auto_ntracer = 0

  ! chemtracer xml file
  IF ((lart_chem) .AND. (lart_chemtracer)) THEN
    IF (TRIM(cart_chemtracer_xml) == '') THEN
      CALL finish('mo_art_init_interface:art_calc_ntracer_and_names',   &
             &    'namelist parameter cart_chemtracer_xml'              &
             &  //' has to be given for lart_chem == .TRUE. and lart_chemtracer = .TRUE.')
    ELSE
      INQUIRE(file = TRIM(cart_chemtracer_xml), EXIST = l_exist)
    
      IF (l_exist) THEN
        CALL art_calc_number_of_art_tracers_xml(TRIM(cart_chemtracer_xml),  &
                               &                auto_ntracer_chemtracer,    &
                               &                'chemtracer',               &
                               &                tracer_names_chemtracer)
        auto_ntracer = auto_ntracer + auto_ntracer_chemtracer
      ELSE
        CALL finish('mo_art_init_interface:art_calc_ntracer_and_names',  &
                &   TRIM(cart_chemtracer_xml)                            &
                & //' could not be found. Check cart_chemtracer_xml.')
      END IF
    END IF
  ELSE
    auto_ntracer_chemtracer = 0
    NULLIFY(tracer_names_chemtracer)
  END IF


  IF ((lart_chem) .AND. (lart_mecca)) THEN
    IF (TRIM(cart_mecca_xml) == '') THEN
      CALL finish('mo_art_init_interface:art_calc_ntracer_and_names',   &
             &    'namelist parameter cart_mecca_xml'                   &
             &  //' has to be given for lart_chem == .TRUE. and lart_mecca = .TRUE.')
    ELSE
      INQUIRE(file = TRIM(cart_mecca_xml), EXIST = l_exist)
    
      IF (l_exist) THEN
        CALL art_calc_number_of_art_tracers_xml(TRIM(cart_mecca_xml),  &
                               &                auto_ntracer_mecca,    &
                               &                'meccatracer',         &
                               &                tracer_names_mecca)
        auto_ntracer = auto_ntracer + auto_ntracer_mecca
      ELSE
        CALL finish('mo_art_init_interface:art_calc_ntracer_and_names',  &
               &    TRIM(cart_mecca_xml)                                 &
               &  //' could not be found. Check cart_mecca_xml.')
      END IF
    END IF
  ELSE
    auto_ntracer_mecca = 0
    NULLIFY(tracer_names_mecca)
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
                               &                auto_ntracer_aerosol,    &
                               &                'aerosol',               &
                               &                tracer_names_aerosol)
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
  
  IF (auto_ntracer > 0) THEN
    ALLOCATE(tracer_names(auto_ntracer))

    ! This setting of the tracer names actually requires the correct order of
    ! including the ART tracers in mo_art_tracer: first aerosols, then chemical
    ! tracers
    IF (auto_ntracer_aerosol > 0) THEN
      tracer_names(1:auto_ntracer_aerosol) = tracer_names_aerosol
    END IF
  
    IF (auto_ntracer_chemtracer > 0) THEN
      tracer_names(auto_ntracer_aerosol+1:auto_ntracer_aerosol + auto_ntracer_chemtracer) &
         &    = tracer_names_chemtracer
    END IF
  
    IF (auto_ntracer_mecca > 0) THEN
      tracer_names(auto_ntracer_chemtracer+1:auto_ntracer_chemtracer + auto_ntracer_mecca) &
         &    = tracer_names_mecca
    END IF
  
    ALLOCATE(ctracer_art(ntracer+auto_ntracer))
    ctracer_art(ntracer+1:) = tracer_names
  
    DEALLOCATE(tracer_names)
    IF (ASSOCIATED(tracer_names_chemtracer)) DEALLOCATE(tracer_names_chemtracer)
    IF (ASSOCIATED(tracer_names_mecca)) DEALLOCATE(tracer_names_mecca)
    IF (ASSOCIATED(tracer_names_aerosol)) DEALLOCATE(tracer_names_aerosol)
  ELSE
    ALLOCATE(ctracer_art(ntracer))
  END IF

  ! Setting iart_ntracer to determined number of tracer
  art_config(1)%iart_ntracer = auto_ntracer

END SUBROUTINE art_calc_ntracer_and_names
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_init_one_dom(jg, p_prog_list, tracer, nest_level)
  IMPLICIT NONE
  INTEGER, INTENT(in) ::  &
    &  jg                   !< patch id
  TYPE(t_var_list), INTENT(in) :: &
    &  p_prog_list          !< list of prognostic variables
  REAL(wp), POINTER :: &
    &  tracer(:,:,:,:)      !< tracer values (ICON kg/kg, Aerosols ??, chemistry  mol/mol)
  INTEGER, INTENT(in) :: &
    &  nest_level           !< beginning with zero in global domain
  ! local variables
#ifdef __ICON_ART
  INTEGER(i8) ::   &
    &   dtime_ms
  REAL(wp) ::  &
    &   dtime_real
  CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) ::  &
    &   dtime_string
  TYPE(timedelta), POINTER ::  &
    &   dt_model
  
  ! set timedelta according to nest level
  dtime_ms = getTotalMilliSecondsTimeDelta(time_config%tc_dt_model,  &
    &                                      time_config%tc_exp_refdate)
  dtime_real = REAL(dtime_ms, wp) / 1000._wp         ! in seconds
  dtime_real = dtime_real / 2._wp**nest_level
  dtime_ms   = NINT(dtime_real*1000, i8)
  CALL getPTStringFromMS(dtime_ms, dtime_string)
  dt_model => newTimedelta(TRIM(dtime_string))
  CALL art_init(jg, dt_model, time_config%tc_exp_refdate,  &
           &    p_prog_list, tracer)

  CALL deallocateTimedelta(dt_model)
#endif
END SUBROUTINE art_init_one_dom

SUBROUTINE art_init_atmo_tracers_nwp(jg, mtime_current, p_nh_state, ext_data, &
                 &                   prm_diag, p_prog, tracer, p_prog_list,   &
                 &                   nest_level)
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
    &  tracer(:,:,:,:)      !< tracer values (ICON kg/kg, Aerosols ??, chemistry  mol/mol)
  TYPE(t_var_list), INTENT(in) :: &
    &  p_prog_list          !< list of prognostic variables
  INTEGER, INTENT(in) :: &
    &  nest_level           !< beginning with zero in global domain

#ifdef __ICON_ART
  IF (lart) THEN
    CALL art_collect_atmo_state_nwp(jg, mtime_current, p_nh_state,  &
                   &                ext_data, prm_diag, p_prog)
    CALL art_init_one_dom(jg, p_prog_list, tracer, nest_level)

    IF ((start_time(jg) <= 0.0_wp) .AND. (.NOT. isRestart())) THEN
      CALL art_init_tracer_values_nwp(jg, tracer, mtime_current, p_prog_list)
    END IF
  END IF
#endif

END SUBROUTINE art_init_atmo_tracers_nwp
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_init_atmo_tracers_echam(jg, mtime_current, p_nh_state, &
                 &                     p_prog, tracer, p_prog_list,   &
                 &                     nest_level)
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
    &  tracer(:,:,:,:)      !< tracer values (ICON kg/kg, Aerosols ??, chemistry  mol/mol)
  TYPE(t_var_list), INTENT(in) :: &
    &  p_prog_list          !< list of prognostic variables
  INTEGER, INTENT(in) :: &
    &  nest_level           !< beginning with zero in global domain

#ifdef __ICON_ART
  IF (lart) THEN
    CALL art_collect_atmo_state_echam(jg, mtime_current, p_nh_state, p_prog)

    CALL art_init_one_dom(jg, p_prog_list, tracer, nest_level)

    IF ((start_time(jg) <= 0.0_wp) .AND. (.NOT. isRestart())) THEN
      CALL art_init_tracer_values_echam(jg, tracer, mtime_current, p_prog_list)
    END IF
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
