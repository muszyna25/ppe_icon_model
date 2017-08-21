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
  USE mo_nonhydro_types,                ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_washoutInt

#ifdef __ICON_ART
! Infrastructure Routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_pollen, t_fields_volc
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_aerosol_utilities,         ONLY: art_air_properties
  USE mo_art_clipping,                  ONLY: art_clip_lt
  USE mo_art_integration,               ONLY: art_integrate_explicit
! Washout Routines
  USE mo_art_washout_volc,              ONLY: art_washout_volc
  USE mo_art_washout_radioact,          ONLY: art_washout_radioact
  USE mo_art_washout_aerosol,           ONLY: art_aerosol_washout
  USE omp_lib
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_washout_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_washout_interface(pt_prog,pt_diag, dtime, p_patch, &
              &                  prm_diag, p_metrics, tracer)
!>
!! Interface for ART-routines dealing with washout
!!
!! This interface calls the ART-routines, if ICON has been 
!! built including the ART-package. Otherwise, this is simply a dummy 
!! routine.
!!
  TYPE(t_nh_prog), TARGET, INTENT(inout) :: & 
    &  pt_prog                           !< Prognostic variables
  TYPE(t_nh_diag), TARGET, INTENT(inout) :: &
    &  pt_diag                           !< Diagnostic variables
  REAL(wp), INTENT(IN)              :: &
    &  dtime                             !< Time step (fast physics)
  TYPE(t_patch), TARGET, INTENT(IN) :: &
    &  p_patch                           !< Patch on which computation is performed
  TYPE(t_nwp_phy_diag), INTENT(IN)  :: &
    &  prm_diag                          !< Diagnostic variables (Physics)
  TYPE(t_nh_metrics)                :: &
    &  p_metrics                         !< Metrics (dz, ...)
  REAL(wp), INTENT(INOUT)           :: &
    &  tracer(:,:,:,:)                   !< Tracer mixing ratios [kg/kg]
  ! Local Variables
  INTEGER                 :: & 
    &  jg, jb, ijsp,         & !< patch id, counter for block loop, conuter for jsp loop
    &  i_startblk, i_endblk, & !< Start and end of block loop
    &  istart, iend,         & !< Start and end of nproma loop
    &  i_rlstart, i_rlend,   & !< Relaxation start and end
    &  i_nchdom,             & !< Number of child domains
    &  nlev,                 & !< Number of levels (equals index of lowest full level)
    &  num_radioact            !< counter for number of radionuclides
  REAL(wp),POINTER        :: &
    &  rho(:,:,:)              !< Pointer to air density [kg m-3]
  REAL(wp),ALLOCATABLE    :: &
    &  wash_rate_m0(:,:),    & !< Washout rates [# m-3 s-1]
    &  wash_rate_m3(:,:)       !< Washout rates [UNIT m-3 s-1], UNIT might be mug, kg
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode
  !-----------------------------------------------------------------------
    
  rho => pt_prog%rho
  
  ! --- Get the loop indizes
  i_nchdom   = MAX(1,p_patch%n_childdom)
  jg         = p_patch%id
  nlev       = p_patch%nlev
  i_rlstart  = grf_bdywidth_c+1
  i_rlend    = min_rlcell_int
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)
  
  IF(lart) THEN 
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_washoutInt)

    IF (art_config(jg)%lart_aerosol) THEN
      ALLOCATE(wash_rate_m0(nproma,nlev))
      ALLOCATE(wash_rate_m3(nproma,nlev))
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)
        CALL art_air_properties(pt_diag%pres(:,:,jb),pt_diag%temp(:,:,jb), &
          &                     istart,iend,1,nlev,jb,p_art_data(jg))
      ENDDO
      
      num_radioact = 0

      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
     
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          CLASS IS (t_fields_2mom)

!$omp parallel do default(shared) private(jb, istart, iend, wash_rate_m0, wash_rate_m3)            
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                &                istart, iend, i_rlstart, i_rlend)
              ! Before washout, the modal parameters have to be calculated
              CALL fields%modal_param(p_art_data(jg)%air_prop%art_free_path(:,:,jb),                &
                &                     istart, iend, nlev, jb, tracer(:,:,jb,:))
              !Washout rate
              IF (.FALSE.) THEN ! Check if qnr is present
                CALL art_aerosol_washout(pt_diag%temp(:,:,jb),                                               &
                   &                tracer(:,:,jb,fields%itr0), fields%density(:,:,jb),                      &
                   &                fields%diameter(:,:,jb),fields%info%sg_ini, tracer(:,:,jb,iqr),          &
                   &                rho(:,:,jb), p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb),               &
                   &                p_art_data(jg)%air_prop%art_free_path(:,:,jb), istart, iend, 15, nlev,   &
                   &                .TRUE., wash_rate_m0(:,:), wash_rate_m3(:,:),                            &
                   &                rrconv_3d=prm_diag%rain_con_rate_3d(:,:,jb), qnr=tracer(:,:,jb,iqnr))
              ELSE
                CALL art_aerosol_washout(pt_diag%temp(:,:,jb),                                               &
                   &                tracer(:,:,jb,fields%itr0), fields%density(:,:,jb),                      &
                   &                fields%diameter(:,:,jb),fields%info%sg_ini, tracer(:,:,jb,iqr),          &
                   &                rho(:,:,jb), p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb),               &
                   &                p_art_data(jg)%air_prop%art_free_path(:,:,jb), istart, iend, 15, nlev,   &
                   &                .TRUE., wash_rate_m0(:,:), wash_rate_m3(:,:),                            &
                   &                rrconv_3d=prm_diag%rain_con_rate_3d(:,:,jb))
              ENDIF
              ! Update mass mixing ratios
              DO ijsp = 1, fields%ntr-1
                CALL art_integrate_explicit(tracer(:,:,jb,fields%itr3(ijsp)),  wash_rate_m3(:,:), dtime,     &
                  &                         istart,iend, nlev, opt_rho = rho(:,:,jb), opt_fac=(1._wp/fields%info%mode_fac))
              ENDDO
              ! Update mass-specific number
              CALL art_integrate_explicit(tracer(:,:,jb,fields%itr0), wash_rate_m0(:,:), dtime,              &
                &                         istart,iend, nlev, opt_rho = rho(:,:,jb))
            ENDDO
!$omp end parallel do
              
          CLASS IS (t_fields_pollen)
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                &                istart, iend, i_rlstart, i_rlend)
              ! CAREFUL: For the time being we are using this washout routine designed for 1-moment volcanic ash
              !          in order to calculate pollen washout. We need to replace this routine by the proper
              !          pollen washout routine from COSMO-ART (see issue #18 in ICON-ART redmine)
              CALL art_washout_volc(dtime,istart, iend, nlev, tracer(:,:,jb,iqr), prm_diag%rain_gsp_rate(:,jb), &
                &                   prm_diag%rain_con_rate(:,jb), prm_diag%rain_con_rate_3d(:,:,jb),            &
                &                   tracer(:,:,jb,fields%itr))
            ENDDO
          CLASS IS (t_fields_volc)
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                &                istart, iend, i_rlstart, i_rlend)
              CALL art_washout_volc(dtime,istart, iend, nlev, tracer(:,:,jb,iqr), prm_diag%rain_gsp_rate(:,jb), &
                &                   prm_diag%rain_con_rate(:,jb), prm_diag%rain_con_rate_3d(:,:,jb),            &
                &                   tracer(:,:,jb,fields%itr))
            ENDDO
          CLASS IS (t_fields_radio)
            num_radioact = num_radioact+1
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                &                istart, iend, i_rlstart, i_rlend)
              CALL art_washout_radioact(rho(:,:,jb), p_metrics%ddqz_z_full(:,:,jb),                          &
                &                       tracer(:,:,jb,iqr),tracer(:,:,jb,iqs),                               &
                &                       prm_diag%rain_gsp_rate(:,jb), prm_diag%rain_con_rate(:,jb),          &
                &                       prm_diag%rain_con_rate_3d(:,:,jb),                                   &
                &                       prm_diag%snow_gsp_rate(:,jb), prm_diag%snow_con_rate(:,jb),          &
                &                       prm_diag%snow_con_rate_3d(:,:,jb), num_radioact,                     &
                &                       fields%fac_wetdep, fields%exp_wetdep, istart, iend, nlev, jb, dtime, &
                &                       tracer(:,:,jb,fields%itr), p_art_data(jg))
            ENDDO
          CLASS DEFAULT
            CALL finish('mo_art_washout_interface:art_washout_interface', &
                 &      'ART: Unknown mode field type')
        END SELECT
                                  
        this_mode => this_mode%next_mode
      ENDDO !associated(this_mode)
    
      ! ----------------------------------
      ! --- Clip the tracers
      ! ----------------------------------
!$omp parallel do default(shared) private(jb, istart, iend)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)
        CALL art_clip_lt(tracer(istart:iend,1:nlev,jb,:),0.0_wp)
      ENDDO
!$omp end parallel do
      
      DEALLOCATE(wash_rate_m0)
      DEALLOCATE(wash_rate_m3)
    ENDIF !lart_aerosol

    IF (timers_level > 3) CALL timer_stop(timer_art_washoutInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  ENDIF !lart
#endif

END SUBROUTINE art_washout_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_washout_interface

