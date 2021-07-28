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
  USE mo_nwp_phy_state,                 ONLY: phy_params
  USE mo_nonhydro_state,                ONLY: p_nh_state_lists
  USE mo_nonhydrostatic_config,         ONLY: kstart_tracer, kstart_moist

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
#ifdef _OPENMP
  USE omp_lib
#endif
  USE mo_art_diagnostics,               ONLY: art_save_aerosol_wet_deposition
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
    &  nlev, nlevp1,         & !< Number of levels (equals index of lowest full level)
    &  num_radioact            !< counter for number of radionuclides
  REAL(wp),POINTER        :: &
    &  rho(:,:,:)              !< Pointer to air density [kg m-3]
  REAL(wp),ALLOCATABLE    :: &
    &  wash_rate_m0(:,:),    & !< Washout rates [# m-3 s-1]
    &  wash_rate_m3(:,:),    & !< Washout rates [UNIT m-3 s-1], UNIT might be mug, kg
    &  rain_con_rate_wo(:,:)   !< convective rain rate for 2mom aerosol washout (eg dust)
  INTEGER                 :: &
    &  kstart_wo               !< start level for washout routine
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode
  !-----------------------------------------------------------------------

  rho => pt_prog%rho

  ! --- Get the loop indizes
  i_nchdom   = MAX(1,p_patch%n_childdom)
  jg         = p_patch%id
  nlev       = p_patch%nlev
  nlevp1     = p_patch%nlevp1
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
      ALLOCATE(rain_con_rate_wo(nproma,nlevp1))

!$omp parallel do default(shared) private(jb, istart, iend)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)
        CALL art_air_properties(pt_diag%pres(:,:,jb),pt_diag%temp(:,:,jb), &
          &                     istart,iend,1,nlev,jb,p_art_data(jg))
      ENDDO
!$omp end parallel do

      num_radioact = 0

      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          CLASS IS (t_fields_2mom)

          kstart_wo = MAX(kstart_tracer(jg,fields%itr0),kstart_moist(jg))

!$omp parallel do default(shared) private(jb, istart, iend, ijsp, wash_rate_m0, wash_rate_m3, rain_con_rate_wo)
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
                   &                p_art_data(jg)%air_prop%art_free_path(:,:,jb), istart, iend,             &
                   &                kstart_wo, nlev, .TRUE., wash_rate_m0(:,:), wash_rate_m3(:,:),           &
                   &                rrconv_3d=prm_diag%rain_con_rate_3d(:,:,jb), qnr=tracer(:,:,jb,iqnr))

                ! Update mass mixing ratios
                DO ijsp = 1, fields%ntr-1
                  CALL art_integrate_explicit(tracer(:,:,jb,fields%itr3(ijsp)),  wash_rate_m3(:,:), dtime,     &
                    &                         istart,iend, nlev, opt_rho = rho(:,:,jb), opt_fac=(1._wp/fields%info%mode_fac))
                ENDDO
                ! Update mass-specific number
                CALL art_integrate_explicit(tracer(:,:,jb,fields%itr0), wash_rate_m0(:,:), dtime,              &
                  &                         istart,iend, nlev, opt_rho = rho(:,:,jb))

              ELSE
                !iart_aero_washout=0: washout with gscp+con precipitation, only 1 call to washout routine
                !iart_aero_washout=1: washout with gscp, con precipitation separately, 2 calls to washout routine
                !iart_aero_washout=2: washout with gscp, con precipitation separately with con/rcucov, 2 calls to washout routine
                !-------------------
                !1st call:
                CALL art_aerosol_washout(pt_diag%temp(:,:,jb),                                               &
                   &                tracer(:,:,jb,fields%itr0), fields%density(:,:,jb),                      &
                   &                fields%diameter(:,:,jb),fields%info%sg_ini, tracer(:,:,jb,iqr),          &
                   &                rho(:,:,jb), p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb),               &
                   &                p_art_data(jg)%air_prop%art_free_path(:,:,jb), istart, iend,             &
                   &                kstart_wo, nlev, .TRUE., wash_rate_m0(:,:), wash_rate_m3(:,:),           &
                   &                rrconv_3d=prm_diag%rain_con_rate_3d(:,:,jb),                             &
                   &                iart_aero_washout=art_config(jg)%iart_aero_washout)


                ! Update mass mixing ratios
                DO ijsp = 1, fields%ntr-1
                   CALL art_integrate_explicit(tracer(:,:,jb,fields%itr3(ijsp)),  wash_rate_m3(:,:), dtime,  &
                        & istart,iend, nlev, opt_rho = rho(:,:,jb), opt_fac=(1._wp/fields%info%mode_fac))
                  ! DIAGNOSTIC: acc_wetdepo_gscp/acc_wetdepo_con/acc_wetdepo_rrsfc of art-tracer
                   CALL art_save_aerosol_wet_deposition(p_art_data(jg),                                      &
                        & wash_rate_m3(:,:), art_config(jg)%iart_aero_washout,                               &
                        & prm_diag%rain_gsp_rate(:,jb)+prm_diag%rain_con_rate(:,jb),                         &
                        & p_metrics%ddqz_z_full(:,:,:), dtime, fields%itr3(ijsp), jb,                        &
                        & istart, iend, kstart_wo, nlev, opt_fac=(1._wp/fields%info%mode_fac))
                ENDDO !ijsp

                ! Update mass-specific number
                CALL art_integrate_explicit(tracer(:,:,jb,fields%itr0), wash_rate_m0(:,:), dtime,            &
                     &                         istart,iend, nlev, opt_rho = rho(:,:,jb))
                ! DIAGNOSTIC: acc_wetdepo_gscp/acc_wetdepo_con/acc_wetdepo_rrsfc of art-tracer
                CALL art_save_aerosol_wet_deposition(p_art_data(jg),                                         &
                     & wash_rate_m0(:,:), art_config(jg)%iart_aero_washout,                                  &
                     & prm_diag%rain_gsp_rate(:,jb)+prm_diag%rain_con_rate(:,jb),                            &
                     & p_metrics%ddqz_z_full(:,:,:), dtime, fields%itr0, jb,                                 &
                     & istart, iend, kstart_wo, nlev )


                IF (art_config(jg)%iart_aero_washout .GT. 0) THEN
                !-------------------
                !2nd call: (washout due to convective precipitation)

                  rain_con_rate_wo(:,:) = prm_diag%rain_con_rate_3d(:,:,jb)
                  IF (art_config(jg)%iart_aero_washout == 2) THEN
                    !scale convective precipitation rate with rcucov
                    rain_con_rate_wo(:,:) = rain_con_rate_wo(:,:)/phy_params(jg)%rcucov
                  ENDIF

                  !add 100 to iart_aero_washout to indicate second call to washout routine
                  CALL art_aerosol_washout(pt_diag%temp(:,:,jb),                                              &
                     &                tracer(:,:,jb,fields%itr0), fields%density(:,:,jb),                     &
                     &                fields%diameter(:,:,jb),fields%info%sg_ini, tracer(:,:,jb,iqr),         &
                     &                rho(:,:,jb), p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb),              &
                     &                p_art_data(jg)%air_prop%art_free_path(:,:,jb), istart, iend,            &
                     &                kstart_wo, nlev, .TRUE., wash_rate_m0(:,:), wash_rate_m3(:,:),          &
                     &                rrconv_3d=rain_con_rate_wo,                                             &
                     &                iart_aero_washout=art_config(jg)%iart_aero_washout+100)

                  IF (art_config(jg)%iart_aero_washout == 2) THEN
                    !rescale washout rate due to convective precipitation with rcucov
                    wash_rate_m3(:,:) = wash_rate_m3(:,:)*phy_params(jg)%rcucov
                    wash_rate_m0(:,:) = wash_rate_m0(:,:)*phy_params(jg)%rcucov
                  ENDIF

                  ! Update mass mixing ratios
                  DO ijsp = 1, fields%ntr-1
                    CALL art_integrate_explicit(tracer(:,:,jb,fields%itr3(ijsp)),  wash_rate_m3(:,:), dtime, &
                      &  istart,iend, nlev, opt_rho = rho(:,:,jb), opt_fac=(1._wp/fields%info%mode_fac))
                   ! DIAGNOSTIC: acc_wetdepo_gscp/acc_wetdepo_con/acc_wetdepo_rrsfc of art-tracer
                    CALL art_save_aerosol_wet_deposition(p_art_data(jg),                                     &
                         & wash_rate_m3(:,:), art_config(jg)%iart_aero_washout+100,                          &
                         & prm_diag%rain_gsp_rate(:,jb)+prm_diag%rain_con_rate(:,jb),                        &
                         & p_metrics%ddqz_z_full(:,:,:), dtime, fields%itr3(ijsp), jb,                       &
                         & istart, iend, kstart_wo, nlev, opt_fac=(1._wp/fields%info%mode_fac))
                  ENDDO
                  ! Update mass-specific number
                  CALL art_integrate_explicit(tracer(:,:,jb,fields%itr0), wash_rate_m0(:,:), dtime,     &
                       &                         istart,iend, nlev, opt_rho = rho(:,:,jb))
                  ! DIAGNOSTIC: acc_wetdepo_gscp/acc_wetdepo_con/acc_wetdepo_rrsfc of art-tracer
                  CALL art_save_aerosol_wet_deposition(p_art_data(jg),                                  &
                       & wash_rate_m0(:,:), art_config(jg)%iart_aero_washout+100,                       &
                       & prm_diag%rain_gsp_rate(:,jb)+prm_diag%rain_con_rate(:,jb),                     &
                       & p_metrics%ddqz_z_full(:,:,:), dtime, fields%itr0, jb,                          &
                       & istart, iend, kstart_wo, nlev )
               ENDIF !2nd call
              ENDIF
            ENDDO
!$omp end parallel do

          CLASS IS (t_fields_pollen)
            kstart_wo = MAX(kstart_tracer(jg,fields%itr),kstart_moist(jg))
!$omp parallel do default(shared) private(jb, istart, iend)
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                &                istart, iend, i_rlstart, i_rlend)
              ! CAREFUL: For the time being we are using this washout routine designed for 1-moment volcanic ash
              !          in order to calculate pollen washout. We need to replace this routine by the proper
              !          pollen washout routine from COSMO-ART (see issue #18 in ICON-ART redmine)
              CALL art_washout_volc(dtime,istart, iend, kstart_wo, nlev, tracer(:,:,jb,iqr),            &
                &                   prm_diag%rain_gsp_rate(:,jb),                                       &
                &                   prm_diag%rain_con_rate(:,jb), prm_diag%rain_con_rate_3d(:,:,jb),    &
                &                   tracer(:,:,jb,fields%itr))
            ENDDO
!$omp end parallel do
          CLASS IS (t_fields_volc)
            kstart_wo = MAX(kstart_tracer(jg,fields%itr),kstart_moist(jg))
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                &                istart, iend, i_rlstart, i_rlend)
              CALL art_washout_volc(dtime,istart, iend, kstart_wo, nlev, tracer(:,:,jb,iqr),             &
                &                   prm_diag%rain_gsp_rate(:,jb),                                        &
                &                   prm_diag%rain_con_rate(:,jb), prm_diag%rain_con_rate_3d(:,:,jb),     &
                &                   tracer(:,:,jb,fields%itr))
            ENDDO
          CLASS IS (t_fields_radio)
            num_radioact = num_radioact+1
            kstart_wo = MAX(kstart_tracer(jg,fields%itr),kstart_moist(jg))
!$omp parallel do default(shared) private(jb, istart, iend)
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                &                istart, iend, i_rlstart, i_rlend)
              CALL art_washout_radioact(rho(:,:,jb), p_metrics%ddqz_z_full(:,:,jb),                          &
                &                       tracer(:,:,jb,iqr),tracer(:,:,jb,iqs),                               &
                &                       prm_diag%rain_gsp_rate(:,jb), prm_diag%rain_con_rate(:,jb),          &
                &                       prm_diag%rain_con_rate_3d(:,:,jb),                                   &
                &                       prm_diag%snow_gsp_rate(:,jb), prm_diag%snow_con_rate(:,jb),          &
                &                       prm_diag%snow_con_rate_3d(:,:,jb), num_radioact,                     &
                &                       fields%fac_wetdep, fields%exp_wetdep, istart, iend, kstart_wo, nlev, &
                &                       jb, dtime, tracer(:,:,jb,fields%itr), p_art_data(jg))
            ENDDO
!$omp end parallel do
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
      DEALLOCATE(rain_con_rate_wo)

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

