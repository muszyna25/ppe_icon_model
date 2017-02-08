!>
!! Provides interface to ART-routines dealing with aerosol-cloud-interactions
!!
!! This module provides an interface to a version of the Seifert and Beheng
!! two-moment cloud microphysics scheme which incorporates prognostic aerosol
!! as calculated by the ART routines.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Daniel Rieger, KIT
!!
!! @par Revision History
!! Initial revision by Daniel Rieger, KIT (2014-11-10)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_clouds_interface

  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_exception,                     ONLY: finish
  USE mo_run_config,                    ONLY: lart
#ifdef __ICON_ART
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_2mom_driver,               ONLY: art_2mom_mcrph,               &
                                          &   art_2mom_mcrph_init
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_volc
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_prepare_aerosol,           ONLY: art_prepare_dust_KL06
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_clouds_interface_2mom
  PUBLIC  :: art_clouds_interface_2mom_init

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clouds_interface_2mom(isize, ke, jg, jb, is, ie, ks, dt, &
                                   & dz, rho, pres, tke, p_trac, tk,    &
                                   & w, prec_r, prec_i, prec_s,         &
                                   & prec_g, prec_h, tkvh, msg_level, l_cv)
  !! Interface for ART: Aerosol-Cloud-Interactions
  !! @par Revision History
  !! Initial revision by Daniel Rieger, KIT (2014-11-10)
  ! Setup variables (Grid, timestep, looping)
  INTEGER,            INTENT (in) :: &
    &  isize, ke,                    & !< grid sizes
    &  jg, jb,                       & !< domain index (p_patch%id)
    &  is, ie, ks                      !< start/end indices
  REAL(wp), INTENT(in)            :: &
    &  dt                              !< time step
  ! Dynamical core variables
  REAL(wp), INTENT(in), TARGET    :: &
    &  dz(:,:),                      & !< Vertical layer thickness
    &  rho(:,:),                     & !< Density
    &  pres(:,:),                    & !< Pressure
    &  tke(:,:),                     & !< Turbulent kinetic energy
    &  w(:,:),                       & !< Vertical velocity
    &  tkvh(:,:)                       !< Turbulent diffusion coefficient for heat
  REAL(wp), INTENT(inout), TARGET :: &
    &  tk(:,:)                         !< Temperature
  ! Tracer fields
  REAL(wp), INTENT(inout), TARGET :: &
    &  p_trac(:,:,:)                   !< Tracer fields
  ! Precip rates, vertical profiles
  REAL(wp), INTENT (inout)        :: &
    &  prec_r(:),                    & !< Precipitation rate for rain
    &  prec_i(:),                    & !< Precipitation rate for ice
    &  prec_s(:),                    & !< Precipitation rate for snow
    &  prec_g(:),                    & !< Precipitation rate for graupel
    &  prec_h(:)                       !< Precipitation rate for hail
  ! Switches
  INTEGER, INTENT (in)            :: &
    &  msg_level                       !< Message level
  LOGICAL, INTENT (in)            :: &
    &  l_cv                            !< Use c_v (true) or c_p (false)
    
#ifdef __ICON_ART
  
  IF (lart) THEN
    
    ! ----------------------------------
    ! --- Call of the coupled ART-twomoment microphysics
    ! ----------------------------------
    IF (art_config(jg)%iart_aci_cold == 6 .OR. art_config(jg)%iart_aci_cold == 7) THEN
      CALL art_prepare_dust_KL06(jg,jb,is,ie,ks,ke,rho,p_trac)
    ENDIF
    CALL art_2mom_mcrph(isize, ke, jg, jb, is, ie, ks, dt,           &
                        & dz, rho, pres, tke, p_trac(:,:,:), tk,    &
                        & w, prec_r, prec_i, prec_s,         &
                        & prec_g, prec_h, tkvh, msg_level, l_cv)
  ELSE
    call finish('mo_art_clouds_interface:art_clouds_interface_2mom', &
         &      'Two moment micophysics with ART aerosol chosen (inwp_gscp=6), but lart=.FALSE.')
  ENDIF !lart
#endif
END SUBROUTINE art_clouds_interface_2mom
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_clouds_interface_2mom_init(msg_level)
  !! Interface for ART: Aerosol-Cloud-Interactions Initialization
  !! @par Revision History
  !! Initial revision by Daniel Rieger, KIT (2014-11-10)
  INTEGER, INTENT(IN) :: &
    &  msg_level           !< message level

  
#ifdef __ICON_ART
  IF (lart) THEN
    CALL art_2mom_mcrph_init(msg_level)
  ELSE
    call finish('mo_art_clouds_interface:art_clouds_interface_2mom_init', &
         &      'Two moment micophysics with ART aerosol chosen (inwp_gscp=6), but lart=.FALSE.')
  ENDIF !lart
#endif

END SUBROUTINE art_clouds_interface_2mom_init
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_clouds_interface
