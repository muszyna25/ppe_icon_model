!>
!! Provides interface to ART-routines dealing with washout
!!
!! This module provides an interface to the ART-routines.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Daniel Rieger, KIT
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Kristina Lundgren, KIT (2012-06-15)
!! Rewritten by Daniel Rieger, KIT (2013-30-09)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_washout_interface

  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_exception,                     ONLY: finish
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag
  USE mo_run_config,                    ONLY: lart,iqr,iqs
  USE mo_nonhydro_state,                ONLY: p_nh_state

#ifdef __ICON_ART
! Infrastructure Routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_volc
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_aerosol_utilities,         ONLY: art_air_properties
  USE mo_art_clipping,                  ONLY: art_clip_lt
! Washout Routines
  USE mo_art_washout_volc,              ONLY: art_washout_volc
  USE mo_art_radioactive,               ONLY: art_washout_radioact
  USE mo_art_washout_aerosol,           ONLY: art_aerosol_washout
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_washout_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_washout_interface(            & !>in
    &          dt_phy_jg,                    & !>in
    &          p_patch,                      & !>in
    &          prm_diag,                     & !>in
    &          p_rho,                        & !>in
    &          p_tracer_new)                   !>inout
!>
!! Interface for ART-routines dealing with washout
!!
!! This interface calls the ART-routines, if ICON has been 
!! built including the ART-package. Otherwise, this is simply a dummy 
!! routine.
!!
  TYPE(t_patch), TARGET, INTENT(IN) :: &
    &  p_patch                             !< patch on which computation is performed
  REAL(wp), INTENT(IN)              :: &
    &  dt_phy_jg                           !< time interval, fast physics
  TYPE(t_nwp_phy_diag), INTENT(IN)  :: &
    &  prm_diag                            !< diagnostic variables
  REAL(wp), INTENT(IN)              :: &          
    &  p_rho(:,:,:)                        !< density of air  [kg/m3]
  REAL(wp), INTENT(INOUT)           :: &
    &  p_tracer_new(:,:,:,:)               !< tracer mixing ratios after transport  [kg/kg]
  INTEGER  :: jg                           !< patch id
    
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode
#endif
    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART

  jg  = p_patch%id

  IF(lart) THEN 
    IF (art_config(jg)%lart_aerosol) THEN
      CALL art_air_properties(p_patch,p_art_data(jg))
      
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
     
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        select type (fields=>this_mode%fields)
        
          class is (t_fields_2mom)
            ! Before washout, the modal parameters have to be calculated
            call fields%modal_param(p_art_data(jg),p_patch,p_tracer_new)
            call art_aerosol_washout(fields,p_patch,dt_phy_jg,prm_diag,p_nh_state(jg), &
                              &      p_tracer_new,p_rho,p_art_data(jg))
                              
          class is (t_fields_volc)
            call art_washout_volc(p_patch,dt_phy_jg,prm_diag,   &
                           &      p_tracer_new,p_rho,fields%itracer)
                           
          class is (t_fields_radio)
            call art_washout_radioact(fields,p_patch,dt_phy_jg,prm_diag,  &
                               &      p_tracer_new,p_rho,p_art_data(jg))
                               
          class default
            call finish('mo_art_washout_interface:art_washout_interface', &
                 &      'ART: Unknown mode field type')
        end select
                                  
        this_mode => this_mode%next_mode
      END DO
    
      ! ----------------------------------
      ! --- Clip the tracers
      ! ----------------------------------
    
      CALL art_clip_lt(p_tracer_new,0.0_wp)
    ENDIF
  ENDIF
#endif

END SUBROUTINE art_washout_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_washout_interface

