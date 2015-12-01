!>
!! mo_art_diagnostics_interface
!! This module provides and interface to routines calculating diagnostic
!! properties of aerosols and trace gases
!!
!! Author: Daniel Rieger, KIT
!!
!! Revision History
!! Initial Daniel Rieger, KIT (2014-08-04)
!!
!! Copyright
!! Institute for Meteorology and Climate Research, KIT, Karlsruhe
!! This software is provided for non-commercial use only.
!!
!! See the LICENSE conditions:
!! http://icon-art.imk-tro.kit.edu/
!!
!!-----------------------------------------------------------------------------
MODULE mo_art_diagnostics_interface

! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_impl_constants,                ONLY: min_rlcell_int
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c
  USE mo_loopindices,                   ONLY: get_indices_c
  USE mo_run_config,                    ONLY: lart
! ART
#ifdef __ICON_ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_aero_optical_props,        ONLY: art_calc_aod550_dust, art_calc_aod550_seas
  USE mo_art_diagnostics,               ONLY: art_volc_diagnostics
  USE mo_art_clipping,                  ONLY: art_clip_lt
#endif

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_diagnostics_interface
    
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_diagnostics_interface(p_patch, rho, pres, p_trac, dz, hml, jg)
  
  TYPE(t_patch), TARGET, INTENT(in) ::  & 
    &  p_patch                !< Patch on which computation is performed
  REAL(wp), INTENT(in)   :: &
    &  rho(:,:,:),          & !< Air density
    &  pres(:,:,:),         & !< Air pressure
    &  dz(:,:,:),           & !< Layer thickness
    &  hml(:,:,:)             !< Height of main layer
  REAL(wp),INTENT(inout) :: &
    &  p_trac(:,:,:,:)        !< Tracer mixing ratios [kg/kg]
  INTEGER, INTENT(in)    :: &
    & jg                      !< Patch id

#ifdef __ICON_ART
  INTEGER                :: &
    &  jb,                   & !< Counter for block loop
    &  i_startblk, i_endblk, & !< Start and end of block loop
    &  istart, iend,         & !< Start and end index of nproma loop
    &  i_rlstart, i_rlend,   & !< Relaxation start and end
    &  i_nchdom,             & !< Number of child domains
    &  nlev                    !< Number of levels (equals index of lowest full level)

  IF (lart) THEN
    IF (art_config(jg)%lart_diag_out) THEN
      ! --- Get the loop indizes
      i_nchdom   = MAX(1,p_patch%n_childdom)
      nlev       = p_patch%nlev
      i_rlstart  = grf_bdywidth_c+1
      i_rlend    = min_rlcell_int
      i_startblk = p_patch%cells%start_blk(i_rlstart,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

      ! ----------------------------------
      ! --- Clipping as convection might produce negative values
      ! ----------------------------------

      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)

        CALL art_clip_lt(p_trac(:,:,jb,:),0.0_wp)

        ! --------------------------------------
        ! --- Calculate aerosol optical depths
        ! --------------------------------------
        IF (art_config(jg)%iart_dust > 0) THEN
          CALL art_calc_aod550_dust(rho(:,:,jb), p_trac(:,:,jb,:), dz(:,:,jb), istart, iend, nlev, jb, p_art_data(jg))
        ENDIF
        IF (art_config(jg)%iart_seasalt > 0) THEN
          CALL art_calc_aod550_seas(rho(:,:,jb), p_trac(:,:,jb,:), dz(:,:,jb), istart, iend, nlev, jb, p_art_data(jg))
        ENDIF

        ! -------------------------------------
        ! --- Calculate volcanic ash products
        ! -------------------------------------
        IF (art_config(jg)%iart_volcano > 0) THEN
          CALL art_volc_diagnostics( rho(:,:,jb), pres(:,:,jb), p_trac(:,:,jb,:), dz(:,:,jb), hml(:,:,jb),  &
            &                        istart, iend, nlev, jb, p_art_data(jg), art_config(jg)%iart_volcano )
        ENDIF

      ENDDO
    ENDIF
  ENDIF
#endif

END SUBROUTINE art_diagnostics_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_diagnostics_interface
