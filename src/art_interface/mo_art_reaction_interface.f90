!>
!! Provides interface to ART-routines dealing with chemical reactions and radioactive decay
!!
!! This module provides an interface to the ART-routines.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Max Bangert, KIT
!!
!! @par Revision History
!! Initial revision by Max Bangert, KIT (2013-02-15)
!! Modifications by Daniel Rieger, KIT (2014-05-22)
!! - Adaption to changes in ART data structure
!! Modifications by Roland Ruhnke, Jennifer Schröter, KIT (2015-08-06)
!! - Splitting between simple chemtracer and full gas phase chemistry
!!   called by iart_chem_mechanism == 0 or 1
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_reaction_interface

  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_linked_list,                   ONLY: t_var_list
  USE mo_run_config,                    ONLY: lart
  USE mtime,                            ONLY: datetime
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_reacInt,            &
                                          &   timer_art_losschem, timer_art_photo
  USE mo_radiation_config,              ONLY: irad_o3
#ifdef __ICON_ART
  USE mo_art_data,                      ONLY: p_art_data, t_art_atmo, t_art_chem
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_decay_radioact,            ONLY: art_decay_radioact
  USE mo_art_chemtracer,                ONLY: art_loss_chemtracer
  USE mo_art_gasphase,                  ONLY: art_loss_gasphase
  USE mo_art_photolysis,                ONLY: art_photolysis
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_radio
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_psc_state,                 ONLY: art_psc_main
  USE mo_art_feedback_icon,             ONLY: art_feedback_o3,     &
                                          &   art_dry_freezing_H2O
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_reaction_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_reaction_interface(jg,current_date,p_dtime,p_prog_list,tracer)

  !>
  !! Interface for ART-routines treating reactions of any kind (chemistry, radioactive decay)
  !!
  !! This interface calls the ART-routines, if ICON has been 
  !! built including the ART-package. Otherwise, this is simply a dummy 
  !! routine.
  !!
  !! @par Revision History
  !! Initial revision by Max Bangert, KIT (2013-02-25)
  !!
  ! atmosphere external data                                
  INTEGER, INTENT(IN)                 :: & 
    &  jg                                !< patch id
  TYPE(datetime), POINTER, INTENT(IN) :: &
    &  current_date                      !< current time and date
  REAL(wp), INTENT(IN)                :: &
    &  p_dtime                           !< time step
  TYPE(t_var_list), INTENT(INOUT)     :: &
    &  p_prog_list                       !< current prognostic state list
  REAL(wp), INTENT(INOUT)           :: &
    &  tracer(:,:,:,:)                   !< tracer mixing ratios (mol/mol)
! Local variables
  INTEGER                           :: &
    &  jb,                             & !< loop index
    &  istart, iend                      !< Start and end of nproma loop
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode
  TYPE(t_art_chem), POINTER :: &
    &  art_chem
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo
  
  !-----------------------------------------------------------------------
 
  IF(lart) THEN

    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_reacInt)

    art_atmo => p_art_data(jg)%atmo
    art_chem => p_art_data(jg)%chem

    IF (art_config(jg)%lart_aerosol) THEN
      ! ----------------------------------
      ! --- Radioactive particles
      ! ----------------------------------
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
    
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          TYPE IS (t_fields_radio)
            DO jb = art_atmo%i_startblk, art_atmo%i_endblk
              CALL art_get_indices_c(jg, jb, istart, iend)

              CALL art_decay_radioact(p_dtime,istart,iend,art_atmo%nlev,     &
                &                     fields%halflife,tracer(:,:,jb,fields%itr))
            ENDDO
        END SELECT
        this_mode => this_mode%next_mode
      END DO
    
      NULLIFY(this_mode)
      
    ENDIF !lart_aerosol
    
    ! ----------------------------------
    ! --- chemical tracer reactions
    ! ----------------------------------

    IF (art_config(jg)%lart_chem) THEN

      ! ----------------------------------
      ! ---  Treat PSCs
      ! ---------------------------------
    
      IF (art_config(jg)%lart_psc) THEN
       CALL art_psc_main(p_art_data(jg)%PSC_meta,    &
                 &       p_dtime,                    &
                 &       art_atmo%temp,              &
                 &       art_atmo%pres,              &
                 &       art_atmo%rho,               &
                 &       art_atmo%dz,                &
                 &       art_atmo%cell_area,         &
                 &       jg, tracer)
      END IF

      ! Dry freezing of H2O
      IF (art_chem%iTRH2O /= 0) THEN
        CALL art_dry_freezing_H2O(jg, tracer(:,:,:,art_chem%iTRH2O))
      END IF


      SELECT CASE(art_config(jg)%iart_chem_mechanism)
        CASE(0,1)

          IF (p_art_data(jg)%OH_chem_meta%is_init) THEN
            IF (timers_level > 3) CALL timer_start(timer_art_photo)

            CALL art_photolysis(jg,                     &
                   &            current_date,           &
                   &            p_dtime,                &
                   &            tracer)
            IF (timers_level > 3) CALL timer_stop(timer_art_photo)
          END IF


          IF (timers_level > 3) CALL timer_start(timer_art_losschem)

          CALL art_loss_chemtracer(jg,             &
                                 & current_date,   &
                                 & p_dtime,        &
                                 & p_prog_list,    &
                                 & tracer)
          IF (timers_level > 3) CALL timer_stop(timer_art_losschem)
        CASE(2)
          IF (timers_level > 3) CALL timer_start(timer_art_photo)

          CALL art_photolysis(jg,                 &
                            & current_date,       &
                            & p_dtime,            &
                            & tracer)

          IF (timers_level > 3) CALL timer_stop(timer_art_photo)
          IF (timers_level > 3) CALL timer_start(timer_art_losschem)

          CALL art_loss_gasphase(current_date,    &
                               & jg,              &
                               & p_dtime,         &
                               & tracer)

          IF (timers_level > 3) CALL timer_stop(timer_art_losschem)

        CASE DEFAULT
          CALL finish('mo_art_reaction_interface:art_reaction_interface', &
               &      'ART: Unknown iart_chem_mechanism')
      END SELECT


      SELECT CASE ( irad_o3 )
        CASE (10)
          IF (art_config(jg)%O3_feedback == 1) THEN 
      
            IF (art_chem%iTRO3 /= 0) THEN
              CALL art_feedback_o3(jg,art_atmo%o3_field_icon,tracer(:,:,:,art_chem%iTRO3))
            ELSE
               CALL finish('mo_art_reaction_interface:art_reaction_interface', &
                    &      'You have chosen ART-ozone feedback, irad_o3 = 10, '&
                    &       //'but O3 is not present')
         
            ENDIF
          ENDIF
      
      END SELECT

      
    ENDIF !lart_chem

    IF (timers_level > 3) CALL timer_stop(timer_art_reacInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  ENDIF ! lart
#endif


END SUBROUTINE art_reaction_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_reaction_interface
