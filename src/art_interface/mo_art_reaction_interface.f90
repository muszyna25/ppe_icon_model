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
  USE mo_linked_list,                   ONLY: t_var_list_ptr
  USE mo_run_config,                    ONLY: lart
  USE mtime,                            ONLY: datetime
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_reacInt,            &
                                          &   timer_art_losschem, timer_art_photo
  USE mo_radiation_config,              ONLY: irad_o3
#ifdef __ICON_ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_chem_data,                 ONLY: t_art_chem_param,       &
                                          &   t_art_chem_indices
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_decay_radioact,            ONLY: art_decay_radioact
  USE mo_art_chemtracer,                ONLY: art_loss_chemtracer
  USE mo_art_full_chemistry,            ONLY: art_loss_full_chemistry
  USE mo_art_photolysis,                ONLY: art_photolysis
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_radio
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_psc_state,                 ONLY: art_psc_main
  USE mo_art_feedback_icon,             ONLY: art_feedback_o3,     &
                                          &   art_dry_freezing_H2O
  USE mo_art_cell_loop,                 ONLY: art_loop_cell_array
  USE mo_art_chem_utils,                ONLY: art_calc_vmr2Nconc,            &
                                          &   art_convert_tracers_mmr_Nconc, &
                                          &   art_convert_tracers_Nconc_mmr
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
  TYPE(t_var_list_ptr), INTENT(INOUT)     :: &
    &  p_prog_list                       !< current prognostic state list
  REAL(wp), INTENT(INOUT)           :: &
    &  tracer(:,:,:,:)                   !< tracer mixing ratios (kg/kg)
! Local variables
  INTEGER                           :: &
    &  jb,                             & !< loop index
    &  istart, iend                      !< Start and end of nproma loop
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode
  TYPE(t_art_chem_indices), POINTER :: &
    &  art_indices
  TYPE(t_art_chem_param), POINTER :: &
    &  art_param
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo
  
  !-----------------------------------------------------------------------
 
  IF (lart) THEN

    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_reacInt)

    art_atmo  => p_art_data(jg)%atmo
    art_param => p_art_data(jg)%chem%param
    art_indices => p_art_data(jg)%chem%indices

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

      ! Convert mass mixing ratio into number concentration in #/cm3
      CALL art_loop_cell_array(jg,p_art_data(jg)%chem%vmr2Nconc,art_calc_vmr2Nconc)

      CALL art_convert_tracers_mmr_Nconc(jg, tracer, p_prog_list, p_art_data(jg)%chem%vmr2Nconc)

      ! ----------------------------------
      ! ---  Treat PSCs
      ! ---------------------------------
    
      IF (art_config(jg)%lart_psc) THEN
       CALL art_psc_main(p_art_data(jg)%chem%PSC_meta,  &
                 &       p_dtime,                       &
                 &       art_atmo%temp,                 &
                 &       art_atmo%pres,                 &
                 &       art_atmo%dz,                   &
                 &       art_atmo%cell_area,            &
                 &       jg, tracer)
      END IF

      ! Dry freezing of H2O
      IF (art_indices%iTRH2O /= 0) THEN
        CALL art_dry_freezing_H2O(jg, tracer(:,:,:,art_indices%iTRH2O),  &
                     &            p_art_data(jg)%chem%vmr2Nconc)
      END IF


      IF ((p_art_data(jg)%chem%param%OH_chem_meta%is_init)  &
        &  .OR. art_config(jg)%lart_mecca) THEN
        IF (timers_level > 3) CALL timer_start(timer_art_photo)
  
        CALL art_photolysis(jg,                           &
               &            tracer,                       &
               &            p_art_data(jg)%chem%vmr2Nconc)
        IF (timers_level > 3) CALL timer_stop(timer_art_photo)
      END  IF
        

      IF (art_config(jg)%lart_chemtracer) THEN
        ! Only calculate the parametrised tracers if there are any of them
        IF (art_param%number_param_tracers > 0) THEN
          IF (timers_level > 3) CALL timer_start(timer_art_losschem)
  
          CALL art_loss_chemtracer(jg,             &
                                 & current_date,   &
                                 & p_dtime,        &
                                 & p_prog_list,    &
                                 & tracer)
          IF (timers_level > 3) CALL timer_stop(timer_art_losschem)
        END IF
      END IF

      IF (art_config(jg)%lart_mecca) THEN
        IF (timers_level > 3) CALL timer_start(timer_art_losschem)

        CALL art_loss_full_chemistry(current_date,  &
                             & jg,                  &
                             & p_dtime,             &
                             & tracer)

        IF (timers_level > 3) CALL timer_stop(timer_art_losschem)
      END IF

      ! Convert number concentration in #/cm3 back to mass mixing ratio
      CALL art_convert_tracers_Nconc_mmr(jg, tracer, p_prog_list, p_art_data(jg)%chem%vmr2Nconc)


      SELECT CASE ( irad_o3 )
        CASE (10)
          IF (art_config(jg)%O3_feedback == 1) THEN 
      
            IF (art_indices%iTRO3 /= 0) THEN
              CALL art_feedback_o3(jg,art_atmo%o3_field_icon,tracer(:,:,:,art_indices%iTRO3))
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
