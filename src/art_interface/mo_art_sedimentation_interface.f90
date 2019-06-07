!>
!! Provides interface to ART-routines dealing with sedimentation
!!
!! This module provides an interface to the ART-routine sedi_volc.
!! The interface is written in such a way, that ICON will compile and run
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Kristina Lundgren, KIT (2012-06-01)
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
MODULE mo_art_sedi_interface

  USE mo_kind,                          ONLY: wp
  USE mo_parallel_config,               ONLY: nproma
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_impl_constants,                ONLY: min_rlcell_int
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c
  USE mo_nonhydro_types,                ONLY: t_nh_prog, t_nh_metrics,t_nh_diag
  USE mo_nonhydro_state,                ONLY: p_nh_state_lists
  USE mo_nonhydrostatic_config,         ONLY: kstart_tracer
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag
  USE mo_run_config,                    ONLY: lart,iqr,iqc,iqv
  USE mo_exception,                     ONLY: finish
  USE mo_advection_vflux,               ONLY: upwind_vflux_ppm, implicit_sedim_tracer
  USE mo_loopindices,                   ONLY: get_indices_c
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_sedInt
#ifdef __ICON_ART
! infrastructure routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_pollen,t_fields_volc
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_clipping,                  ONLY: art_clip_lt
  USE mo_art_config,                    ONLY: art_config
! sedimentation and deposition routines
  USE mo_art_sedi_1mom,                 ONLY: art_sedi_1mom
  USE mo_art_sedi_2mom,                 ONLY: art_calc_v_sed, art_calc_sed_flx
  USE mo_art_depo_2mom,                 ONLY: art_calc_v_dep, art_store_v_dep
  USE mo_art_drydepo_radioact,          ONLY: art_drydepo_radioact
  USE mo_art_diagnostics,               ONLY: art_get_diag_tracer_index
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_sedi_interface


CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_sedi_interface(p_patch, p_dtime, p_prog, p_metrics, rho, p_diag, &
              &                   prm_diag, tracer, lprint_cfl)
!! Interface for ART routines for calculation of
!! sedimentation and deposition velocities
!!
!! @par Revision History
!! Initial revision by Kristina Lundgren, KIT (2012-06-01)
!! Modifications by Daniel Rieger, KIT (2014-05-22)
!! - Adaption to changes in ART data structure
!! - Calculation of deposition velocities
  TYPE(t_patch), TARGET, INTENT(IN) :: &
    &  p_patch                           !< Patch on which computation is performed
  TYPE(t_nh_prog), INTENT(IN)       :: &
    &  p_prog                            !< Current prognostic state
  TYPE(t_nh_metrics), INTENT(IN)    :: &
    &  p_metrics                         !< Metrical fields
  TYPE(t_nh_diag), INTENT(IN)       :: &
    &  p_diag                            !< diagnostic variables
  TYPE(t_nwp_phy_diag), INTENT(in)  :: &
    &  prm_diag                          !< list of diagnostic fields (physics)
  REAL(wp), INTENT(IN)              :: &
    &  p_dtime,                        & !< Time step (dynamics)
    &  rho(:,:,:)                        !< density of air at full levels
  REAL(wp), INTENT(INOUT)           :: &
    &  tracer(:,:,:,:)                   !< tracer mixing ratio [kg/kg]
  LOGICAL, INTENT(IN)               :: &               
    &  lprint_cfl                        !< determines if vertical CFL number shall be printed in upwind_vflux_ppm
! Local Variables
  REAL(wp), ALLOCATABLE        :: &
    &  p_upflux_sed(:,:,:),       & !< upwind flux at half levels due to sedimentation
    &  rhodz_new(:,:,:),          & !< density * height of full layer
    !MB:
    &  rho_inv(:,:,:),            & !< 1 / density [m^3/kg]
    &  vsed0(:,:),                & !< Sedimentation velocities 0th moment [m s-1]
    &  vsed3(:,:),                & !< Sedimentation velocities 3th moment [m s-1]
    &  vdep0(:),                  & !< Deposition velocities 0th moment [m s-1]
    &  vdep3(:)                     !< Deposition velocities 3th moment [m s-1]
  REAL(wp),POINTER             :: &
    &  flx_contra_vsed(:,:,:)       !< Flux due to sedimentation (can be mass or number)
  !MB:
  REAL(wp),POINTER             :: &
    &  vsed(:,:,:)                  !< sedimentation velocity (can be for mass or number)
  REAL(wp) ::       &
    &  sedim_update,              & !< tracer tendency due to sedimentation
    &  dt_sub                       !< integration time step of one substep
  INTEGER,ALLOCATABLE          :: &
    &  jsp_ar(:)                    !< Contains index of mass mixing ratios in tracer container
  INTEGER                      :: &
    &  jc, jk, jkp1, jb,          & !< loop index for: cell, level full, level half, block
    &  nlev, nlevp1,              & !< number of full / half levels
    &  nblks,                     & !< number of blocks
    &  i_nchdom,                  & !< Number of child domains
    &  i_rlstart, i_rlend,        & !< Relaxation start and end
    &  i_startblk, i_endblk,      & !< Start and end of block loop
    &  istart,iend,               & !< Start and end of nproma loop
    &  jg, jsp, i, n,             & !< jg: patch id, jsp/i/n: counters
    &  iubc=0,                    & !< upper boundary condition 0 = none
    &  itype_vlimit=2               !< Monotone flux limiter
  LOGICAL                      :: &
    &  lcompute_gt,               & !< compute geometrical terms
    &  lcleanup_gt                  !< clean up geometrical terms.
  INTEGER                      :: &
    &  idx_diag,                  & !< Index of tracer in diagnostics container
    &  nart_substeps                !< Number of substeps for sedimentation

#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode

  lcompute_gt=.TRUE. ! compute geometrical terms
  lcleanup_gt=.TRUE. ! clean up geometrical terms. obs. this i currently done for all components. improvement:
                     ! compute values for first component, cleanup after last component.
  NULLIFY(flx_contra_vsed)
  NULLIFY(vsed)

  jg     = p_patch%id
  nlevp1 = p_patch%nlevp1
  nlev   = p_patch%nlev
  nblks  = p_patch%nblks_c

  nart_substeps = art_config(jg)%nart_substeps_sedi
  dt_sub = p_dtime/REAL(nart_substeps)

  !Get all cell entities, except halos
  i_nchdom  = MAX(1,p_patch%n_childdom)
  i_rlstart = grf_bdywidth_c+1
  i_rlend   = min_rlcell_int
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

  IF (lart) THEN
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_sedInt)

    IF (art_config(jg)%lart_aerosol) THEN
      ALLOCATE(vsed0(nproma,nlev),vsed3(nproma,nlev))
      ALLOCATE(vdep0(nproma),vdep3(nproma))
      ALLOCATE(p_upflux_sed(nproma,nlevp1,nblks))
      ALLOCATE(rhodz_new(nproma,nlev,nblks))

      !MB:
      IF ( art_config(jg)%cart_type_sedim == "impl" ) THEN
        ALLOCATE( rho_inv(nproma,nlev,nblks) )
        !rho_inv(:,:,:) = 0.0_wp  ! might be necessary only for debugging purposes 
      END IF

!$omp parallel do default(shared) private(jb, jk, jc, istart, iend)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
          &                istart, iend, i_rlstart, i_rlend)
        DO jk = 1, nlev
          DO jc = istart, iend
            rhodz_new(jc,jk,jb) = rho(jc,jk,jb) * p_metrics%ddqz_z_full(jc,jk,jb)
          ENDDO
        ENDDO

        !MB:
        IF ( art_config(jg)%cart_type_sedim == "impl" ) THEN
          DO jk = 1, nlev
            DO jc = istart, iend
              rho_inv(jc,jk,jb) = 1.0_wp / rho(jc,jk,jb)
            END DO
          END DO
        END IF

      ENDDO
!$omp end parallel do


      !drieg: Necessary depending on parameterizations that will be used
      !    CALL art_air_properties(p_patch,p_art_data(jg))
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)

          CLASS IS (t_fields_2mom)
!$omp parallel do default(shared) private(jb, i, istart, iend, vsed0, vsed3, vdep0, vdep3, jsp_ar)
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
              &                istart, iend, i_rlstart, i_rlend)
            ! Before sedimentation/deposition velocity calculation, the modal parameters have to be calculated
            CALL fields%modal_param(p_art_data(jg)%air_prop%art_free_path(:,:,jb),                &
              &                     istart, iend, nlev, jb, tracer(:,:,jb,:))
            ! Calculate sedimentation velocities for 0th and 3rd moment
            CALL art_calc_v_sed(p_metrics%ddqz_z_full(:,:,jb),                                    &
              &     p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb), fields%density(:,:,jb),         &
              &     fields%diameter(:,:,jb), fields%info%exp_aero, fields%knudsen_nr(:,:,jb),     &
              &     istart, iend, nlev, vsed0(:,:), vsed3(:,:))

            !MB>>
            IF ( art_config(jg)%cart_type_sedim == "impl" ) THEN
              ! This is a bit ugly: flx_contra_vsed is misused as storage for vsed.
              ! The following code is slightly modified from subr. 'art_calc_sed_flx'.
 
              !fields%flx_contra_vsed0(:,:,jb) = 0.0_wp    ! only necessary for debugging
              !fields%flx_contra_vsed3(:,:,jb) = 0.0_wp    ! only necessary for debugging
  
              ! --- Interpolate on half levels and
  
              !MB: remember the index convention: w(k) = w_{k-1/2}
              DO jk = 2, nlev ! all except top and bottom level
                DO jc = istart, iend
                  fields%flx_contra_vsed0(jc,jk,jb) =  p_metrics%wgtfac_c(jc,jk,jb)  * vsed0(jc,jk)   &
                    &                         + (1._wp-p_metrics%wgtfac_c(jc,jk,jb)) * vsed0(jc,jk-1)
                  fields%flx_contra_vsed3(jc,jk,jb) =  p_metrics%wgtfac_c(jc,jk,jb)  * vsed3(jc,jk)   &
                    &                         + (1._wp-p_metrics%wgtfac_c(jc,jk,jb)) * vsed3(jc,jk-1)
                ENDDO  ! jc
              ENDDO ! jk
  
              ! top and bottom level
              DO jc =  istart, iend
                fields%flx_contra_vsed0(jc,1     ,jb) = vsed0(jc,1)
                fields%flx_contra_vsed0(jc,nlevp1,jb) = vsed0(jc,nlev)
                fields%flx_contra_vsed3(jc,1     ,jb) = vsed3(jc,1)
                fields%flx_contra_vsed3(jc,nlevp1,jb) = vsed3(jc,nlev)
              ENDDO ! jc
              !MB<<

            ELSE IF ( art_config(jg)%cart_type_sedim == "expl" ) THEN

              ! Calculate massflux due to sedimentation for 0th and 3rd moment
              CALL art_calc_sed_flx(vsed0(:,:),vsed3(:,:),p_metrics%wgtfac_c(:,:,jb),             &
                &     p_metrics%wgtfacq_c(:,:,jb),rho(:,:,jb),                                    &
                &     p_diag%rho_ic(:,:,jb), istart, iend, nlev,                                  &
                &     fields%flx_contra_vsed0(:,:,jb),fields%flx_contra_vsed3(:,:,jb))

            ELSE
              CALL finish('mo_art_sedimentation_interface:art_sedimentation_interface', &
                 &         'ART: Unknown value for cart_type_sedim')
            END IF
           
            ! Calculate deposition velocities for 0th and 3rd moment
            CALL art_calc_v_dep(p_diag%temp(:,nlev,jb), p_diag%temp_ifc(:,nlev+1,jb),             &
              &     p_diag%u(:,nlev,jb), p_diag%v(:,nlev,jb), rho(:,nlev,jb),prm_diag%tcm(:,jb),  &
              &     prm_diag%tch(:,jb),tracer(:,nlev,jb,iqv),tracer(:,nlev,jb,iqc),               &
              &     tracer(:,nlev,jb,iqr), p_prog%theta_v(:,nlev,jb), prm_diag%gz0(:,jb),         &
              &     p_metrics%ddqz_z_full(:,nlev,jb),                                             &
              &     p_art_data(jg)%air_prop%art_dyn_visc(:,nlev,jb), fields%diameter(:,nlev,jb),  &
              &     fields%info%exp_aero,fields%knudsen_nr(:,nlev,jb), vsed0(:,nlev),             &
              &     vsed3(:,nlev), istart, iend, nlev, vdep0(:), vdep3(:))
            ! Store deposition velocities for the use in turbulence scheme
            ALLOCATE(jsp_ar(fields%ntr-1))
            DO i =1, fields%ntr-1
              jsp_ar(i) = fields%itr3(i)
            ENDDO
            CALL art_store_v_dep(vdep0(:), vdep3(:), (fields%ntr-1), jsp_ar,                     &
              &     fields%itr0, art_config(jg)%nturb_tracer,                                    &
              &     p_prog%turb_tracer(jb,:),istart,iend,p_art_data(jg)%turb_fields%vdep(:,jb,:))
            DEALLOCATE(jsp_ar)

          ENDDO !jb
!$omp end parallel do

          CLASS IS (t_fields_pollen)
!$omp parallel do default(shared) private(jb, istart, iend)
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
              &                istart, iend, i_rlstart, i_rlend)
            CALL art_sedi_1mom(p_diag%temp(:,:,jb), p_diag%pres(:,:,jb), rho(:,:,jb),    &
              &                p_diag%rho_ic(:,:,jb), p_metrics%wgtfac_c(:,:,jb),        &
              &                p_metrics%wgtfacq_c(:,:,jb), fields%diam, fields%rho,     &
              &                p_prog%turb_tracer(jb,:), istart, iend, nlev, jb,         &
              &                fields%itr, art_config(jg), p_art_data(jg),               &
              &                fields%flx_contra_vsed(:,:,jb))
          ENDDO
!$omp end parallel do

          CLASS IS (t_fields_volc)
!$omp parallel do default(shared) private(jb, istart, iend)
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
              &                istart, iend, i_rlstart, i_rlend)
            CALL art_sedi_1mom(p_diag%temp(:,:,jb), p_diag%pres(:,:,jb), rho(:,:,jb),    &
              &                p_diag%rho_ic(:,:,jb), p_metrics%wgtfac_c(:,:,jb),        &
              &                p_metrics%wgtfacq_c(:,:,jb), fields%diam, fields%rho,     &
              &                p_prog%turb_tracer(jb,:), istart, iend, nlev, jb,         &
              &                fields%itr, art_config(jg), p_art_data(jg),               &
              &                fields%flx_contra_vsed(:,:,jb))
          ENDDO
!$omp end parallel do

          CLASS IS (t_fields_radio)
          ! A constant deposition velocity is used / required
!$omp parallel do default(shared) private(jb, istart, iend)
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
              &                istart, iend, i_rlstart, i_rlend)
            CALL art_drydepo_radioact(p_prog%turb_tracer(jb,:), art_config(jg),        &
              &                       fields%itr, fields%vdep_const, istart, iend, jb, &
              &                       p_art_data(jg))
          ENDDO
!$omp end parallel do

          CLASS DEFAULT
          CALL finish('mo_art_sedimentation_interface:art_sedimentation_interface', &
            &         'ART: Unknown mode field type')

        END SELECT

        !MB:
        IF ( art_config(jg)%cart_type_sedim == "expl" ) THEN

          DO n = 1, nart_substeps

            !drieg: Necessary depending on parameterizations that will be used
            !    CALL art_air_properties(p_patch,p_art_data(jg))
            ! Select type of mode
            SELECT TYPE (fields=>this_mode%fields)

              CLASS IS (t_fields_2mom)

              DO i=1, fields%ntr            !< loop through the tracer contained in the mode
                IF (i .NE. fields%ntr) THEN
                  jsp = fields%itr3(i)
                  flx_contra_vsed => fields%flx_contra_vsed3
                ELSE
                  jsp = fields%itr0
                  flx_contra_vsed => fields%flx_contra_vsed0
                ENDIF

                ! upwind_vflux_ppm is internally OpenMP parallelized
                CALL upwind_vflux_ppm(p_patch, tracer(:,:,:,jsp),             &
                  &                   iubc, flx_contra_vsed, dt_sub,          &
                  &                   lcompute_gt, lcleanup_gt, itype_vlimit, &
                  &                   p_metrics%ddqz_z_full,                  &
                  &                   rhodz_new, lprint_cfl,                  &
                  &                   p_upflux_sed(:,:,:),                    &
                  &                   opt_rlstart=i_rlstart,                  &
                  &                   opt_rlend=i_rlend,                      &
                  &                   opt_elev=nlevp1,                        &
                  &                   opt_slev=kstart_tracer(jg,jsp),         &
                  &                   opt_ti_slev=kstart_tracer(jg,jsp))

                ! ----------------------------------
                ! --- update mixing ratio after sedimentation
                ! ----------------------------------

                ! get tracer index in diagnostics container (for diagnostic ACC_SEDIM)
                CALL art_get_diag_tracer_index(p_nh_state_lists(jg)%diag_list, 'acc_sedim', jsp, idx_diag)

                IF ( idx_diag > 0 .AND. n == 1 ) THEN
!$omp parallel do default(shared) private(jb, jk, jc, istart, iend)
                  DO jb = i_startblk, i_endblk
                    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                      &                istart, iend, i_rlstart, i_rlend)
                    DO jk = kstart_tracer(jg,jsp), nlev
                      DO jc = istart, iend
                        ! save sedimentation in lowest model level (ACC_SEDIM)
                        ! integrate sedim_update vertically
                        p_art_data(jg)%diag%acc_sedim(jc,jb,idx_diag) =         &
                          & p_art_data(jg)%diag%acc_sedim(jc,jb,idx_diag)       &
                          & + tracer(jc,jk,jb,jsp) * rhodz_new(jc,jk,jb)
                      ENDDO!jc
                    ENDDO !jk
                  ENDDO !jb
!$omp end parallel do
                END IF

!$omp parallel do default(shared) private(jb, jk, jc, jkp1, sedim_update, istart, iend)
                DO jb = i_startblk, i_endblk
                  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                    &                istart, iend, i_rlstart, i_rlend)
                  DO jk = kstart_tracer(jg,jsp), nlev
                    ! index of bottom half level
                    jkp1 = jk + 1
                    DO jc = istart, iend
                      sedim_update = (dt_sub * ( p_upflux_sed(jc,jk,  jb)   &
                        &                      - p_upflux_sed(jc,jkp1,jb) ) &
                        &         / rhodz_new(jc,jk,jb))
                      tracer(jc,jk,jb,jsp) =   tracer(jc,jk,jb,jsp) - sedim_update
                    ENDDO!jc
                  ENDDO !jk
                ENDDO !jb
!$omp end parallel do

                IF ( idx_diag > 0 .AND. n == nart_substeps ) THEN
!$omp parallel do default(shared) private(jb, jk, jc, istart, iend)
                  DO jb = i_startblk, i_endblk
                    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                      &                istart, iend, i_rlstart, i_rlend)
                    DO jk = kstart_tracer(jg,jsp), nlev
                      DO jc = istart, iend
                        ! save sedimentation in lowest model level (ACC_SEDIM)
                        ! integrate sedim_update vertically
                        p_art_data(jg)%diag%acc_sedim(jc,jb,idx_diag) =         &
                          & p_art_data(jg)%diag%acc_sedim(jc,jb,idx_diag)       &
                          & - tracer(jc,jk,jb,jsp) * rhodz_new(jc,jk,jb)
                      ENDDO!jc
                    ENDDO !jk
                  ENDDO !jb
!$omp end parallel do
                END IF

              ENDDO !i

              CLASS IS (t_fields_pollen)

              flx_contra_vsed => fields%flx_contra_vsed
              jsp = fields%itr

              ! upwind_vflux_ppm is internally OpenMP parallelized
              CALL upwind_vflux_ppm(p_patch, tracer(:,:,:,jsp),             &
                &                   iubc, flx_contra_vsed, dt_sub,          &
                &                   lcompute_gt, lcleanup_gt, itype_vlimit, &
                &                   p_metrics%ddqz_z_full,                  &
                &                   rhodz_new, lprint_cfl,                  &
                &                   p_upflux_sed(:,:,:),                    &
                &                   opt_rlstart=i_rlstart,                  &
                &                   opt_rlend=i_rlend,                      &
                &                   opt_elev=nlevp1,                        &
                &                   opt_slev=kstart_tracer(jg,jsp),         &
                &                   opt_ti_slev=kstart_tracer(jg,jsp))

              ! ----------------------------------
              ! --- update mixing ratio after sedimentation
              ! ----------------------------------

!$omp parallel do default(shared) private(jb, jk, jc, jkp1, sedim_update, istart, iend)
              DO jb = i_startblk, i_endblk
                CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                  &                istart, iend, i_rlstart, i_rlend)
                DO jk = kstart_tracer(jg,jsp), nlev
                  ! index of bottom half level
                  jkp1 = jk + 1
                  DO jc = istart, iend
                    sedim_update = (dt_sub * ( p_upflux_sed(jc,jk,  jb)   &
                      &                      - p_upflux_sed(jc,jkp1,jb) )  &
                      &         / rhodz_new(jc,jk,jb))
                    tracer(jc,jk,jb,jsp) =   tracer(jc,jk,jb,jsp) - sedim_update
                  ENDDO!jc
                ENDDO !jk
              ENDDO !jb
!$omp end parallel do

              CLASS IS (t_fields_volc)

              flx_contra_vsed => fields%flx_contra_vsed
              jsp = fields%itr

              ! upwind_vflux_ppm is internally OpenMP parallelized
              CALL upwind_vflux_ppm(p_patch, tracer(:,:,:,jsp),             &
                &                   iubc, flx_contra_vsed, dt_sub,          &
                &                   lcompute_gt, lcleanup_gt, itype_vlimit, &
                &                   p_metrics%ddqz_z_full,                  &
                &                   rhodz_new, lprint_cfl,                  &
                &                   p_upflux_sed(:,:,:),                    &
                &                   opt_rlstart=i_rlstart,                  &
                &                   opt_rlend=i_rlend,                      &
                &                   opt_elev=nlevp1,                        &
                &                   opt_slev=kstart_tracer(jg,jsp),         &
                &                   opt_ti_slev=kstart_tracer(jg,jsp))

              ! ----------------------------------
              ! --- update mixing ratio after sedimentation
              ! ----------------------------------

!$omp parallel do default(shared) private(jb, jk, jc, jkp1, sedim_update, istart, iend)
              DO jb = i_startblk, i_endblk
                CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                  &                istart, iend, i_rlstart, i_rlend)
                DO jk = kstart_tracer(jg,jsp), nlev
                  ! index of bottom half level
                  jkp1 = jk + 1
                  DO jc = istart, iend
                    sedim_update = (dt_sub * ( p_upflux_sed(jc,jk,  jb)   &
                      &                      - p_upflux_sed(jc,jkp1,jb) )  &
                      &         / rhodz_new(jc,jk,jb))
                    tracer(jc,jk,jb,jsp) =   tracer(jc,jk,jb,jsp) - sedim_update
                  ENDDO!jc
                ENDDO !jk
              ENDDO !jb
!$omp end parallel do

            END SELECT

          ENDDO !nart_substeps


        ELSE IF ( art_config(jg)%cart_type_sedim == "impl" ) THEN
          !MB>>

          ! Select type of mode
          SELECT TYPE (fields=>this_mode%fields)

            CLASS IS (t_fields_2mom)

            DO i=1, fields%ntr            !< loop through the tracer contained in the mode
              IF (i .NE. fields%ntr) THEN
                jsp = fields%itr3(i)
                vsed => fields%flx_contra_vsed3
              ELSE
                jsp = fields%itr0
                vsed => fields%flx_contra_vsed0
              ENDIF

              ! get tracer index in diagnostics container (for diagnostic ACC_SEDIM)
              CALL art_get_diag_tracer_index(p_nh_state_lists(jg)%diag_list, 'acc_sedim', jsp, idx_diag)

              IF ( idx_diag > 0 ) THEN
!$omp parallel do default(shared) private(jb, jk, jc, istart, iend)
                DO jb = i_startblk, i_endblk
                  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                    &                istart, iend, i_rlstart, i_rlend)
                  DO jk = 1, nlev
                    DO jc = istart, iend
                      ! save sedimentation in lowest model level (ACC_SEDIM)
                      ! integrate sedim_update vertically
                      p_art_data(jg)%diag%acc_sedim(jc,jb,idx_diag) =         &
                        & p_art_data(jg)%diag%acc_sedim(jc,jb,idx_diag)       &
                        & + tracer(jc,jk,jb,jsp) * rhodz_new(jc,jk,jb)
                    ENDDO!jc
                  ENDDO !jk
                ENDDO !jb
!$omp end parallel do
              END IF

              ! implicit_sedim_tracer is internally OpenMP parallelized
              CALL implicit_sedim_tracer( tracer(:,:,:,jsp),  &
                &                  rho, rho_inv,              &
              !!&                  vsed, vsed_old,            &
                &                  vsed, vsed,           &
                &                  p_dtime,              &
                &                  p_patch, p_metrics,   &
                &                  i_rlstart, i_rlend  )

              IF ( idx_diag > 0 ) THEN
!$omp parallel do default(shared) private(jb, jk, jc, istart, iend)
                DO jb = i_startblk, i_endblk
                  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                    &                istart, iend, i_rlstart, i_rlend)
                  DO jk = 1, nlev
                    DO jc = istart, iend
                      ! save sedimentation in lowest model level (ACC_SEDIM)
                      ! integrate sedim_update vertically
                      p_art_data(jg)%diag%acc_sedim(jc,jb,idx_diag) =         &
                        & p_art_data(jg)%diag%acc_sedim(jc,jb,idx_diag)       &
                        & - tracer(jc,jk,jb,jsp) * rhodz_new(jc,jk,jb)
                    ENDDO!jc
                  ENDDO !jk
                ENDDO !jb
!$omp end parallel do
              END IF

            ENDDO !i

            CLASS IS (t_fields_pollen)

            jsp = fields%itr
            vsed => fields%flx_contra_vsed

            ! implicit_sedim_tracer is internally OpenMP parallelized
            CALL implicit_sedim_tracer( tracer(:,:,:,jsp),  &
              &                  rho, rho_inv,              &
            !!&                  vsed, vsed_old,            &
              &                  vsed, vsed,           &
              &                  p_dtime,              &
              &                  p_patch, p_metrics,   &
              &                  i_rlstart, i_rlend  )

            CLASS IS (t_fields_volc)

            jsp = fields%itr
            vsed => fields%flx_contra_vsed

            ! implicit_sedim_tracer is internally OpenMP parallelized
            CALL implicit_sedim_tracer( tracer(:,:,:,jsp),  &
              &                  rho, rho_inv,              &
            !!&                  vsed, vsed_old,            &
              &                  vsed, vsed,           &
              &                  p_dtime,              &
              &                  p_patch, p_metrics,   &
              &                  i_rlstart, i_rlend  )

          END SELECT
          !MB<<

        ELSE
          CALL finish('mo_art_sedimentation_interface:art_sedimentation_interface', &
            &         'ART: Unknown value for cart_type_sedim')
        END IF

        this_mode => this_mode%next_mode
      ENDDO ! end do while associated

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

      DEALLOCATE(p_upflux_sed)
      DEALLOCATE(vsed0,vsed3,vdep0,vdep3)
      DEALLOCATE(rhodz_new)

      IF ( art_config(jg)%cart_type_sedim == "impl" ) THEN
        DEALLOCATE( rho_inv )
      END IF

    ENDIF !lart_aerosol

    IF (timers_level > 3) CALL timer_stop(timer_art_sedInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  ENDIF !lart

#endif

END SUBROUTINE art_sedi_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_sedi_interface

