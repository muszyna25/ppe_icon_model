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
!!
!!  chemical tracer section
!!  by Roland Ruhnke (RR), Jennifer Schroeter (JS), KIT
!!
!! @par Revision History
!! Initial revision by Max Bangert, KIT (2013-02-15)
!! Modifications by Daniel Rieger, KIT (2014-05-22)
!! - Adaption to changes in ART data structure
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
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_linked_list,                   ONLY: t_var_list
  USE mo_nonhydro_types,                ONLY: t_nh_diag
  USE mo_run_config,                    ONLY: lart
  USE mo_ext_data_types,         ONLY: t_external_data
  USE mo_datetime,               ONLY: t_datetime
  USE mo_nonhydro_types,         ONLY: t_nh_metrics, t_nh_prog, t_nh_diag
#ifdef __ICON_ART
  USE mo_art_radioactive,               ONLY: art_decay_radioact
  USE mo_art_chemtracer,                ONLY: art_loss_chemtracer !<JS
  USE mo_art_gasphase,                  ONLY: art_loss_gasphase !< JS
  USE mo_art_photolysis,                ONLY: art_photolysis !< JS
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_radio
  USE mo_art_config,                    ONLY: art_config
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_reaction_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_reaction_interface(ext_data, p_patch,datetime,p_dtime,p_prog_list,p_prog,p_metrics,p_diag,p_tracer_now)

  !>
  !! Interface for ART-routines treating reactions of any kind (chemistry, radioactive decay)
  !!
  !! This interface calls the ART-routines, if ICON has been 
  !! built including the ART-package. Otherwise, this is simply a dummy 
  !! routine.
  !!
  !! @par Revision History
  !! Initial revision by Max Bangert, KIT (2013-02-25)
    ! atmosphere external data                                
    TYPE(t_external_data),               INTENT(IN) :: ext_data
  TYPE(t_patch), TARGET, INTENT(IN) ::  & 
    &  p_patch                             !< patch on which computation is performed
!<JS
    TYPE(t_datetime), INTENT(IN) :: datetime
  REAL(wp), INTENT(IN)              ::  &
    &  p_dtime                             !< time step
  TYPE(t_nh_diag), INTENT(IN)       ::  &
    &  p_diag                              !< list of diagnostic fields
  TYPE(t_var_list), INTENT(IN)      ::  &
    &  p_prog_list                         !< current prognostic state list

    TYPE(t_nh_prog), INTENT(IN) ::p_prog


    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
  REAL(wp), INTENT(INOUT)           ::  &
    &  p_tracer_now(:,:,:,:)               !< tracer mixing ratios (specific concentrations)

  REAL(wp), POINTER ::  &  !< density of air 
    &  p_rho(:,:,:)              !< [kg/m3]
    
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode
#endif
  INTEGER  :: jg                           !< domain index
  !-----------------------------------------------------------------------
 
#ifdef __ICON_ART
  jg  = p_patch%id
  p_rho => p_prog%rho
  IF(lart) THEN
    IF (art_config(jg)%lart_aerosol) THEN
      ! ----------------------------------
      ! --- Radioactive particles
      ! ----------------------------------
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
    
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        select type (fields=>this_mode%fields)
          type is (t_fields_radio)
            CALL  art_decay_radioact(p_patch,p_dtime,p_tracer_now(:,:,:,fields%itracer),fields%halflife) 
        end select                  
        this_mode => this_mode%next_mode
      END DO
    
      NULLIFY(this_mode)
      
    ENDIF !lart_aerosol
    
    ! ----------------------------------
    ! --- chemical tracer reactions
    ! ----------------------------------
  !< JS added splitting between simple chemtracer and full gas phase chemistry
  !< called by iart_chem_mechanism == 0 or 1

    IF (art_config(jg)%lart_chem) THEN
    
      IF (art_config(jg)%iart_chem_mechanism == 0) THEN
        CALL art_loss_chemtracer(p_patch,           &
               & p_dtime,                           &
               & p_prog_list,                       &
               & p_diag,                            &
               & p_tracer_now)
      ELSEIF (art_config(jg)%iart_chem_mechanism == 1) THEN
        CALL art_photolysis(ext_data,               &
               & p_patch,                           &
               & datetime,                          &
               & p_dtime,                           &
               & p_prog_list,                       &
               & p_prog,                            &
               & p_diag,                            &
               & p_rho,                             &
               & p_metrics,                         &
               & p_tracer_now)

        CALL art_loss_gasphase(p_patch,p_dtime,p_prog_list,p_diag,p_tracer_now)
      ENDIF 
    ENDIF

  ENDIF ! lart
#endif

END SUBROUTINE art_reaction_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_reaction_interface
