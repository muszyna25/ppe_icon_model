!>
!! Provides interface to ART-routines dealing with emissions
!!
!! This module provides an interface to the ART emission routines.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Daniel Rieger, KIT
!! @author Daniel Reinert, DWD
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2012-01-27)
!! Modification by Kristina Lundgren, KIT (2012-01-30)
!! - Modification for dealing with the ART-routine emission_volc.
!! Rewritten by Daniel Rieger, KIT (2013-09-30)
!! - Usage of the generalized ART infrastructure
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_emission_interface

  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_exception,                     ONLY: finish
  USE mo_nonhydro_types,                ONLY: t_nh_diag
  USE mo_ext_data_types,                ONLY: t_external_data
  USE mo_run_config,                    ONLY: lart,                         &
                                          &   iCS137,iI131,iTE132,          &
                                          &   iZR95,iXE133,iI131g,          &
                                          &   iI131o,iBA140,iRU103
    
#ifdef __ICON_ART
! Infrastructure Routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_volc
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_aerosol_utilities,         ONLY: art_air_properties
! Emission Routines
  USE mo_art_emission_volc,             ONLY: art_organize_emission_volc
  USE mo_art_radioactive,               ONLY: art_emiss_radioact
  USE mo_art_emission_seas,             ONLY: art_emission_seas
  USE mo_art_emission_dust,             ONLY: art_emission_dust
  USE mo_art_chemtracer,                ONLY: art_emiss_chemtracer
  USE mo_art_config,                    ONLY: art_config
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_emission_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_emission_interface(ext_data,p_patch,p_dtime,p_rho,p_diag,p_tracer_now)
  !! Interface for ART: Emissions
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-01-27)
  !! Modification by Kristina Lundgren, KIT (2012-01-30)
  !! Rewritten by Daniel Rieger, KIT (2013-09-30)
  TYPE(t_external_data), INTENT(IN) ::  &
    &  ext_data                            !< atmosphere external data
  TYPE(t_patch), TARGET, INTENT(IN) ::  & 
    &  p_patch                             !< patch on which computation is performed
  REAL(wp), INTENT(IN)              ::  &
    &  p_dtime                             !< time step
  REAL(wp), INTENT(INOUT)           ::  &
    &  p_rho(:,:,:)                        !< density of air [kg/m3]
  TYPE(t_nh_diag), INTENT(IN)       ::  &
    &  p_diag                              !< list of diagnostic fields
  REAL(wp), INTENT(INOUT)           ::  &
    &  p_tracer_now(:,:,:,:)               !< tracer mixing ratios [kg/kg]
  ! Local variables
  INTEGER  ::  jg                          !< patch id
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode     !< pointer to current aerosol mode

  jg  = p_patch%id
     
  IF (lart) THEN

    CALL art_air_properties(p_patch,p_art_data(jg))
       
    IF (art_config(jg)%lart_aerosol) THEN
       
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
     
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        select type (fields=>this_mode%fields)
          type is (t_fields_2mom)
            ! Before emissions, the modal parameters have to be calculated
            call fields%modal_param(p_art_data(jg),p_patch,p_tracer_now)
            ! Now the according emission routine has to be found
            select case(TRIM(fields%info%name))
              case ('seasa')
                call art_emission_seas(fields,p_patch,p_dtime,p_rho,p_tracer_now)
              case ('seasb')
                call art_emission_seas(fields,p_patch,p_dtime,p_rho,p_tracer_now)
              case ('seasc')
                call art_emission_seas(fields,p_patch,p_dtime,p_rho,p_tracer_now)
              case ('dusta')
                call art_emission_dust(fields,p_patch,p_dtime,p_rho,p_tracer_now)
              case ('dustb')
                call art_emission_dust(fields,p_patch,p_dtime,p_rho,p_tracer_now)
              case ('dustc')
                call art_emission_dust(fields,p_patch,p_dtime,p_rho,p_tracer_now)
              case default
                call finish('mo_art_emission_interface:art_emission_interface', &
                     &      'No according emission routine to mode'//TRIM(fields%info%name))
            end select
          class is (t_fields_radio)
            select case(TRIM(fields%info%name))
              case ('CS137')
                call art_emiss_radioact(p_patch,p_dtime,p_rho,p_tracer_now(:,:,:,iCS137),373)
              case ('I131')
                call art_emiss_radioact(p_patch,p_dtime,p_rho,p_tracer_now(:,:,:,iI131),340)
              case ('TE132')
                call art_emiss_radioact(p_patch,p_dtime,p_rho,p_tracer_now(:,:,:,iTE132),325)
              case ('ZR95')
                call art_emiss_radioact(p_patch,p_dtime,p_rho,p_tracer_now(:,:,:,iZR95),184)
              case ('XE133')
                call art_emiss_radioact(p_patch,p_dtime,p_rho,p_tracer_now(:,:,:,iXE133),355)
              case ('I131g')
                call art_emiss_radioact(p_patch,p_dtime,p_rho,p_tracer_now(:,:,:,iI131g),870)
              case ('I131o')
                call art_emiss_radioact(p_patch,p_dtime,p_rho,p_tracer_now(:,:,:,iI131o),880)
              case ('BA140')
                call art_emiss_radioact(p_patch,p_dtime,p_rho,p_tracer_now(:,:,:,iBA140),384)
              case ('RU103')
                call art_emiss_radioact(p_patch,p_dtime,p_rho,p_tracer_now(:,:,:,iRU103),220)
              ! And Default...
              case default
                call finish('mo_art_emission_interface:art_emission_interface', &
                     &      'No according emission routine to mode')
            end select
          class is (t_fields_volc)
            ! nothing to do here, see below
          class default
            call finish('mo_art_emission_interface:art_emission_interface', &
                 &      'ART: Unknown mode field type')
        end select
        this_mode => this_mode%next_mode
      ENDDO
  
      ! ----------------------------------
      ! --- volcano emissions
      ! ----------------------------------
    
      IF (art_config(jg)%iart_volcano == 1) THEN
        CALL art_organize_emission_volc(p_patch,p_dtime,p_rho,p_tracer_now) 
      ENDIF
    ENDIF !lart_aerosol
    
    ! ----------------------------------
    ! --- emissions of chemical tracer
    ! ----------------------------------
  
    IF (art_config(jg)%lart_chem) THEN
    
      IF (art_config(jg)%iart_chem_mechanism == 0) THEN
        CALL art_emiss_chemtracer(ext_data,p_patch,p_dtime,p_diag,p_tracer_now)
      ENDIF
      
    ENDIF
    
  ENDIF !lart
       
#endif

END SUBROUTINE art_emission_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emission_interface

