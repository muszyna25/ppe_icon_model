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
  USE mo_impl_constants,                ONLY: min_rlcell_int
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c
  USE mo_loopindices,                   ONLY: get_indices_c
  USE mo_parallel_config,               ONLY: nproma
  USE mo_exception,                     ONLY: finish
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag
  USE mo_run_config,                    ONLY: lart,iqr,iqnr,iqs
  USE mo_nonhydro_types,                ONLY: t_nh_prog, t_nh_diag

#ifdef __ICON_ART
! Infrastructure Routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_volc
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_aerosol_utilities,         ONLY: art_air_properties
  USE mo_art_clipping,                  ONLY: art_clip_lt
  USE mo_art_integration,               ONLY: art_integrate_explicit
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
SUBROUTINE art_washout_interface(pt_prog,pt_diag, dtime, p_patch, &
              &                  prm_diag, rho, p_trac)
!>
!! Interface for ART-routines dealing with washout
!!
!! This interface calls the ART-routines, if ICON has been 
!! built including the ART-package. Otherwise, this is simply a dummy 
!! routine.
!!
  TYPE(t_nh_prog), TARGET, INTENT(inout) :: & 
    &  pt_prog                           !<the prognostic variables
  TYPE(t_nh_diag), TARGET, INTENT(inout) :: &
    &  pt_diag                           !<the diagnostic variables
  REAL(wp), INTENT(IN)              :: &
    &  dtime                             !< time interval, fast physics
  TYPE(t_patch), TARGET, INTENT(IN) :: &
    &  p_patch                           !< patch on which computation is performed
  TYPE(t_nwp_phy_diag), INTENT(IN)  :: &
    &  prm_diag                          !< diagnostic variables
  REAL(wp), INTENT(IN)              :: &          
    &  rho(:,:,:)                        !< density of air  [kg/m3]
  REAL(wp), INTENT(INOUT)           :: &
    &  p_trac(:,:,:,:)                   !< tracer mixing ratios after transport  [kg/kg]
  ! Local Variables
  INTEGER                 :: & 
    &  jg, jb, ijsp,         & !< patch id, counter for block loop, conuter for jsp loop
    &  i_startblk, i_endblk, & !< Start and end of block loop
    &  istart, iend,         & !< Start and end of nproma loop
    &  i_rlstart, i_rlend,   & !< 
    &  i_nchdom,             & !< 
    &  nlev                    !< Number of levels (equals index of lowest full level)
  REAL(wp),ALLOCATABLE    :: &
    &  washout_rate(:,:)       !< Washout rates [UNIT kg-1 s-1] or [UNIT m-3 s-1], UNIT might be mug, kg or just a number
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode
    !-----------------------------------------------------------------------
    
  ! --- Get the loop indizes
  i_nchdom   = MAX(1,p_patch%n_childdom)
  jg         = p_patch%id
  nlev       = p_patch%nlev
  i_rlstart  = grf_bdywidth_c+1
  i_rlend    = min_rlcell_int
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)
  
  
  IF(lart) THEN 
    IF (art_config(jg)%lart_aerosol) THEN
      ALLOCATE(washout_rate(nproma,nlev))
      CALL art_air_properties(p_patch,p_art_data(jg))
      
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
     
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        select type (fields=>this_mode%fields)
        
          class is (t_fields_2mom)
            ! Before washout, the modal parameters have to be calculated
            call fields%modal_param(p_art_data(jg),p_patch,p_trac)
            
            ! This DO loop will be outside the DO WHILE ASSOCIATED loop as soon as modal_param is rewritten for single blocks
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                &                istart, iend, i_rlstart, i_rlend)
              
              !Washout rate
              IF (.FALSE.) THEN ! Check if qnr is present
                CALL art_aerosol_washout(pt_diag%temp(:,:,jb), pt_diag%pres(:,:,jb),                         &
                   &                p_trac(:,:,jb,fields%info%i_number_conc), fields%density(:,:,jb),        &
                   &                fields%diameter(:,:,jb),fields%info%sg_ini, p_trac(:,:,jb,iqr),          &
                   &                rho(:,:,jb), p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb),               &
                   &                p_art_data(jg)%air_prop%art_free_path(:,:,jb), istart, iend, 15, nlev,   &
                   &                .TRUE., washout_rate(:,:), rrconv_3d=prm_diag%rain_con_rate_3d(:,:,jb),  &
                   &                qnr=p_trac(:,:,jb,iqnr))
              ELSE
                CALL art_aerosol_washout(pt_diag%temp(:,:,jb), pt_diag%pres(:,:,jb),                         &
                   &                p_trac(:,:,jb,fields%info%i_number_conc), fields%density(:,:,jb),        &
                   &                fields%diameter(:,:,jb),fields%info%sg_ini, p_trac(:,:,jb,iqr),          &
                   &                rho(:,:,jb), p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb),               &
                   &                p_art_data(jg)%air_prop%art_free_path(:,:,jb), istart, iend, 15, nlev,   &
                   &                .TRUE., washout_rate(:,:), rrconv_3d=prm_diag%rain_con_rate_3d(:,:,jb))
              ENDIF
!              ! Update mass mixing ratios
!              DO ijsp = 1, fields%info%njsp
!                CALL art_integrate_explicit(p_trac(:,:,jb,fields%info%jsp(ijsp)),  washout_rate(:,:), dtime,          &
!                  &                         istart,iend, nlev, opt_rho = rho(:,:,jb),                                 &
!                  &                         opt_fac = (1._wp/(fields%info%mode_fac * fields%info%factnum)))
!              ENDDO
!              ! Update mass-specific number
!              CALL art_integrate_explicit(p_trac(:,:,jb,fields%info%i_number_conc), washout_rate(:,:), dtime,         &
!                &                         istart,iend, nlev, opt_rho = rho(:,:,jb))
            ENDDO
              
          class is (t_fields_volc)
            call art_washout_volc(p_patch,dtime,prm_diag,   &
                           &      p_trac,rho,fields%itracer)
                           
          class is (t_fields_radio)
            call art_washout_radioact(fields,p_patch,dtime,prm_diag,  &
                               &      p_trac,rho,p_art_data(jg))
                               
          class default
            call finish('mo_art_washout_interface:art_washout_interface', &
                 &      'ART: Unknown mode field type')
        end select
                                  
        this_mode => this_mode%next_mode
      END DO
    
      ! ----------------------------------
      ! --- Clip the tracers
      ! ----------------------------------
    
      CALL art_clip_lt(p_trac,0.0_wp)
      DEALLOCATE(washout_rate)
    ENDIF
  ENDIF
#endif

END SUBROUTINE art_washout_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_washout_interface

