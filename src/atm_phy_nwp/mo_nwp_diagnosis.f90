!>
!! @brief diagnosis of physics after physic's call 
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <Pilar Ripodas, DWD>
!!
!!
!! @par Revision History
!! first implementation by Pilar Ripodas, DWD (2011-03)
!! generalized overlap by Martin Koehler, DWD (2014-04)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
#include "consistent_fma.inc"
#include "icon_contiguous_defines.inc"
!----------------------------

MODULE mo_nwp_diagnosis


  USE mo_kind,               ONLY: wp

  USE mo_impl_constants,     ONLY: itccov, itconv, itradheat, itturb, itsfc, itfastphy, &
    &                              min_rlcell_int
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_exception,          ONLY: message, message_text
  USE mo_model_domain,       ONLY: t_patch
  USE mo_run_config,         ONLY: iqv, iqc, iqi, iqr, iqs,  &
                                   iqni, iqg, iqh, iqnc, iqm_max, iqgl, iqhl
  USE mo_grid_config,        ONLY: n_dom, n_dom_start
  USE mo_timer,              ONLY: ltimer, timer_start, timer_stop, timer_nh_diagnostics
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_nh_state
  USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_math_divrot,        ONLY: rot_vertex
  USE mo_intp,               ONLY: verts2cells_scalar
  USE mo_parallel_config,    ONLY: nproma, proc0_offloading
  USE mo_lnd_nwp_config,     ONLY: nlev_soil, ntiles_total
  USE mo_nwp_lnd_types,      ONLY: t_lnd_diag, t_wtr_prog, t_lnd_prog
  USE mo_physical_constants, ONLY: tmelt, grav, cpd, vtmpc1
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_advection_config,   ONLY: advection_config
  USE mo_io_config,          ONLY: lflux_avg, t_var_in_output, uh_max_zmin, uh_max_zmax, &
    &                              luh_max_out, uh_max_nlayer
  USE mo_sync,               ONLY: global_max, global_min
  USE mo_vertical_coord_table,  ONLY: vct_a
  USE mo_satad,              ONLY: sat_pres_water, spec_humi
  USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp
  USE mo_opt_nwp_diagnostics,ONLY: calsnowlmt, cal_cape_cin, maximize_field_lpi, compute_field_tcond_max, &
                                   compute_field_uh_max, compute_field_vorw_ctmax, compute_field_w_ctmax, &
                                   compute_field_dbz3d_lin, maximize_field_dbzctmax, &
                                   compute_field_echotop, compute_field_echotopinm, compute_field_dursun
  USE mo_nwp_ww,             ONLY: ww_diagnostics, ww_datetime
  USE mtime,                 ONLY: datetime, timeDelta, getTimeDeltaFromDateTime,  &
    &                              deallocateTimedelta, newTimeDelta, newDatetime, &
    &                              deallocateDatetime, event
  USE mo_util_mtime,         ONLY: is_event_active
  USE mo_exception,          ONLY: finish
  USE mo_math_constants,     ONLY: pi
  USE mo_statistics,         ONLY: time_avg
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_nwp_parameters,     ONLY: t_phy_params
  USE mo_time_config,        ONLY: time_config
  USE mo_nwp_tuning_config,  ONLY: lcalib_clcov, max_calibfac_clcl
  USE mo_upatmo_impl_const,  ONLY: idamtr
  USE mo_mpi,                ONLY: p_io, p_comm_work, p_bcast

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: nwp_statistics
  PUBLIC  :: nwp_diag_for_output
  PUBLIC  :: nwp_opt_diagnostics
  PUBLIC  :: nwp_opt_diagnostics_2
  PUBLIC  :: nwp_diag_output_1
  PUBLIC  :: nwp_diag_output_2
  PUBLIC  :: nwp_diag_output_minmax_micro


 !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_diagnosis'



CONTAINS

  !>
  !! Computation of time averages, accumulated variables and vertical integrals
  !!
  !! Computation of time averages, accumulated variables and vertical integrals 
  !! for output. The statistics are valid from the beginning of the forecast 
  !! to the output time.
  !!
  !! @par Revision History
  !! Add calculation of high-, mid-, and low-level cloud cover, height
  !! of base and top of convection  by Helmut Frank, DWD (2013-01-17)
  !! Add height of 0 deg C level    by Helmut Frank, DWD (2013-03-11)
  !!
  !!
  SUBROUTINE nwp_statistics(lcall_phy_jg,                 & !in
                            & dt_phy_jg, p_sim_time,      & !in
                            & ext_data, kstart_moist,     & !in
                            & ih_clch, ih_clcm,           & !in
                            & pt_patch, p_metrics,        & !in
                            & pt_prog, pt_prog_rcf,       & !in
                            & pt_diag,                    & !inout
                            & prm_diag, lnd_diag,         & !inout 
                            & linit                       ) !in  
                            

    LOGICAL,            INTENT(IN)   :: lcall_phy_jg(:) !< physics package time control (switches)
                                                        !< for domain jg
    LOGICAL, OPTIONAL,  INTENT(IN)   :: linit           !< initialization flag
    REAL(wp),           INTENT(IN)   :: dt_phy_jg(:)    !< time interval for all physics
                                                        !< packages on domain jg
    REAL(wp),           INTENT(IN)   :: p_sim_time

    TYPE(t_patch),      INTENT(IN)   :: pt_patch    !<grid/patch info.
    TYPE(t_nh_diag),    INTENT(INOUT):: pt_diag     !<the diagnostic variables
    TYPE(t_nh_prog),    INTENT(IN)   :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog),    INTENT(IN)   :: pt_prog_rcf !<the prognostic variables (with
                                                    !< red. calling frequency for tracers!
    TYPE(t_nh_metrics), INTENT(in)   :: p_metrics
    TYPE(t_external_data),INTENT(IN) :: ext_data    !< external data

    TYPE(t_nwp_phy_diag), INTENT(inout):: prm_diag
    TYPE(t_lnd_diag),     INTENT(inout):: lnd_diag      !< diag vars for sfc
 
    INTEGER,           INTENT(IN)  :: kstart_moist
    INTEGER,           INTENT(IN)  :: ih_clch, ih_clcm

    ! Local
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices

    REAL(wp):: r_sim_time
    REAL(wp):: t_wgt                   !< weight for running time average

    INTEGER :: jc,jk,jb,jg      ! block index
    INTEGER :: jt               ! tracer loop index
    LOGICAL :: lacc             ! OpenACC flag


  !-----------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_nh_diagnostics)

    IF(PRESENT(linit)) THEN
      lacc = .NOT. linit
    ELSE
      lacc = .FALSE.
    ENDIF

    jg        = pt_patch%id

    ! Inverse of simulation time
    r_sim_time = 1._wp/MAX(1.e-6_wp, p_sim_time)

    ! time average weight
    t_wgt = dt_phy_jg(itfastphy)/MAX(1.e-6_wp, p_sim_time)


    ! exclude nest boundary interpolation zone
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)
    

    ! Calculate vertical integrals of moisture quantities and cloud cover if 
    ! averages of these variables are requested for output
    IF (atm_phy_nwp_config(jg)%lcalc_moist_integral_avg) THEN
      CALL calc_moist_integrals(pt_patch, p_metrics,        & !in
                              & pt_prog, pt_prog_rcf,       & !in
                              & ext_data, kstart_moist,     & !in
                              & ih_clch, ih_clcm,           & !in
                              & pt_diag, prm_diag,          & !inout
                              & lacc                        ) !in
    ENDIF

    ! Calculation of average/accumulated values since model start
    !
    ! Compute
    !
    ! wind
    !-----------
    ! - maximum gust (including convective contribution)
    ! - max/min 2m temperature
    !
    ! cloud/rain
    !-----------
    ! - total precipitation amount
    ! - time averaged precipitation rates (total, grid-scale, convective)
    ! - time averaged total cloud cover
    ! - time averaged TQV, TQC, TQI, TQR, TQS
    ! - time averaged TQV_DIA, TQC_DIA, TQI_DIA
    !
    ! soil
    !-----
    ! - surface water runoff; sum over forecast
    ! - soil water runoff; sum over forecast
    !
    ! turbulent fluxes
    !-----------------
    ! - surface latent heat flux
    ! - surface latent heat flux from bare soil 
    ! - surface sensible heat flux
    ! - surface moisture flux
    ! - surface u/v-momentum flux (turbulent, sso, resolved)
    !
    ! radiative fluxes
    !------------------
    ! - top net solar radiation
    ! - top down solar radiation
    ! - top net thermal radiation
    ! - surface net solar radiation
    ! - surface net thermal radiation
    ! - surface shortwave diffuse downward radiation
    ! - surface shortwave diffuse upward radiation
    ! - surface shortwave direct downward radiation
    ! - surface downward photosynthetically active flux

!$OMP PARALLEL
    IF ( p_sim_time <= 1.e-6_wp) THEN

      ! ensure that extreme value fields are equal to instantaneous fields 
      ! at VV=0.
      ! In addition, set extreme value fields to instantaneous fields prior 
      ! to first regular time step (i.e. for IAU) 

!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

!DIR$ IVDEP
        !$acc parallel default(present) if(lacc)
        !$acc loop gang vector
        DO jc = i_startidx, i_endidx

          ! set to instantaneous values
          prm_diag%gust10(jc,jb)  = prm_diag%dyn_gust(jc,jb) + prm_diag%con_gust(jc,jb)
          prm_diag%tmax_2m(jc,jb) = prm_diag%t_2m(jc,jb)
          prm_diag%tmin_2m(jc,jb) = prm_diag%t_2m(jc,jb)
        ENDDO
        !$acc end parallel

      ENDDO  ! jb
!$OMP END DO

    ELSE  ! regular time steps
  
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

!DIR$ IVDEP
        !$acc parallel default(present) if(lacc)
        !$acc loop gang vector
        DO jc = i_startidx, i_endidx

          ! maximum 10m gust, including convective contribution
          ! (reset is done on a regular basis in reset_action)
          prm_diag%gust10(jc,jb) = MAX(prm_diag%gust10(jc,jb),                       &
            &                    prm_diag%dyn_gust(jc,jb) + prm_diag%con_gust(jc,jb) )
          
          ! total precipitation
          prm_diag%tot_prec(jc,jb) = prm_diag%prec_gsp(jc,jb) + prm_diag%prec_con(jc,jb)

          ! time averaged total precipitation rate
          prm_diag%tot_prec_rate_avg(jc,jb) = prm_diag%tot_prec(jc,jb) &
            &                               * r_sim_time

          ! time averaged grid scale precipitation rate
          prm_diag%prec_gsp_rate_avg(jc,jb) = prm_diag%prec_gsp(jc,jb) &
            &                               * r_sim_time

          ! time averaged convective precipitation rate
          prm_diag%prec_con_rate_avg(jc,jb) = prm_diag%prec_con(jc,jb) &
            &                               * r_sim_time

        ENDDO  ! jc
        !$acc end parallel

        IF (atm_phy_nwp_config(jg)%lcalc_moist_integral_avg) THEN
!DIR$ IVDEP
          !$acc parallel default(present) if(lacc)
          !$acc loop gang vector
          DO jc = i_startidx, i_endidx
            ! time averaged total cloud cover
            prm_diag%clct_avg(jc,jb) = time_avg(prm_diag%clct_avg(jc,jb), &
              &                                 prm_diag%clct    (jc,jb), &
              &                                 t_wgt)
          ENDDO  ! jc
          !$acc end parallel

          ! time averaged tracer vertical integrals (mass concentrations only)  
          !$acc parallel default(present) if(lacc)
          !$acc loop gang vector collapse(2)
          DO jt = 1, iqm_max
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              pt_diag%tracer_vi_avg(jc,jb,jt) = (1._wp - t_wgt)*pt_diag%tracer_vi_avg(jc,jb,jt) &
                &                              + t_wgt * pt_diag%tracer_vi(jc,jb,jt)
            ENDDO  ! jc
          ENDDO  ! jt
          !$acc end parallel

         ! time averaged TQV_DIA, TQC_DIA, TQI_DIA
          !$acc parallel default(present) if(lacc)
          !$acc loop gang vector collapse(2)
          DO jt = 1, 3
!DIR$ IVDEP
           DO jc = i_startidx, i_endidx
             prm_diag%tot_cld_vi_avg(jc,jb,jt) = (1._wp - t_wgt)*prm_diag%tot_cld_vi_avg(jc,jb,jt) &
               &                                + t_wgt * prm_diag%tot_cld_vi(jc,jb,jt)
            ENDDO
          ENDDO  ! jt
          !$acc end parallel
        ENDIF

        IF (lcall_phy_jg(itsfc)) THEN
          !$acc parallel default(present) if(lacc)
          !$acc loop gang vector collapse(2)
          DO jt=1,ntiles_total
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              lnd_diag%runoff_s_t(jc,jb,jt) = lnd_diag%runoff_s_t(jc,jb,jt) + lnd_diag%runoff_s_inst_t(jc,jb,jt)
              lnd_diag%runoff_g_t(jc,jb,jt) = lnd_diag%runoff_g_t(jc,jb,jt) + lnd_diag%runoff_g_inst_t(jc,jb,jt)
            END DO
          END DO
          !$acc end parallel
        END IF


        ! max/min 2m temperature
        !
        ! note that we do not use the instantaneous aggregated 2m temperature prm_diag%t_2m, 
        ! but the instantaneous max/min over all tiles. In case of no tiles both are equivalent.
        IF (lcall_phy_jg(itturb)) THEN
!DIR$ IVDEP
          !$acc parallel default(present) if(lacc)
          !$acc loop gang vector
          DO jc = i_startidx, i_endidx
            prm_diag%tmax_2m(jc,jb) = MAX(prm_diag%t_tilemax_inst_2m(jc,jb), prm_diag%tmax_2m(jc,jb) )
            prm_diag%tmin_2m(jc,jb) = MIN(prm_diag%t_tilemin_inst_2m(jc,jb), prm_diag%tmin_2m(jc,jb) )
          END DO
          !$acc end parallel
        ENDIF


        IF (lflux_avg) THEN

          IF (lcall_phy_jg(itturb)) THEN
!DIR$ IVDEP
            !$acc parallel default(present) if(lacc)
            !$acc loop gang vector
            DO jc = i_startidx, i_endidx
              ! ATTENTION:
              ! the sign, in the output all fluxes must be positive downwards

              ! time averaged surface latent heat flux
              prm_diag%alhfl_s(jc,jb) = time_avg(prm_diag%alhfl_s(jc,jb), &
                &                                prm_diag%lhfl_s (jc,jb), & 
                &                                t_wgt) 

              ! time averaged surface latent heat flux from bare soil
              prm_diag%alhfl_bs(jc,jb)= time_avg(prm_diag%alhfl_bs(jc,jb),& 
                &                                prm_diag%lhfl_bs (jc,jb),& 
                &                                t_wgt)

              ! time averaged surface sensible heat flux
              prm_diag%ashfl_s(jc,jb) = time_avg(prm_diag%ashfl_s(jc,jb), & 
                &                                prm_diag%shfl_s (jc,jb), & 
                &                                t_wgt)

              ! time averaged surface moisture flux
              prm_diag%aqhfl_s(jc,jb) = time_avg(prm_diag%aqhfl_s(jc,jb), & 
                &                                prm_diag%qhfl_s (jc,jb), & 
                &                                t_wgt )

              ! time averaged surface u-momentum flux turbulence
              prm_diag%aumfl_s(jc,jb) = time_avg(prm_diag%aumfl_s(jc,jb), &
                &                                prm_diag%umfl_s (jc,jb), &
                &                                t_wgt )

              ! time averaged surface v-momentum flux turbulence
              prm_diag%avmfl_s(jc,jb) = time_avg(prm_diag%avmfl_s(jc,jb), &
                &                                prm_diag%vmfl_s (jc,jb), &
                &                                t_wgt )

            ENDDO  ! jc
            !$acc end parallel

            !$acc parallel default(present) if(lacc)
            !$acc loop gang vector collapse(2)
            DO jk = 1, nlev_soil
!DIR$ IVDEP
              DO jc = i_startidx, i_endidx
              prm_diag%alhfl_pl(jc,jk,jb) = time_avg(prm_diag%alhfl_pl(jc,jk,jb), &
                &                                    prm_diag%lhfl_pl (jc,jk,jb), &
                &                                    t_wgt)
              ENDDO  ! jc
            ENDDO  ! jk
            !$acc end parallel

            IF (atm_phy_nwp_config(jg)%lcalc_extra_avg) THEN
!DIR$ IVDEP
              !$acc parallel default(present) if(lacc)
              !$acc loop gang vector
              DO jc = i_startidx, i_endidx
                ! time averaged surface u-momentum flux SSO
                prm_diag%astr_u_sso(jc,jb) = time_avg(prm_diag%astr_u_sso(jc,jb), &
                  &                                   prm_diag%str_u_sso (jc,jb), &
                  &                                   t_wgt )

                ! time averaged surface v-momentum flux SSO
                prm_diag%astr_v_sso(jc,jb) = time_avg(prm_diag%astr_v_sso(jc,jb), &
                  &                                   prm_diag%str_v_sso (jc,jb), &
                  &                                   t_wgt )

                ! time averaged surface u-momentum flux resolved
                prm_diag%adrag_u_grid(jc,jb) = time_avg(prm_diag%adrag_u_grid(jc,jb), &
                  &                                prm_diag%drag_u_grid (jc,jb), &
                  &                                t_wgt )

                ! time averaged surface v-momentum flux resolved
                prm_diag%adrag_v_grid(jc,jb) = time_avg(prm_diag%adrag_v_grid(jc,jb), &
                  &                                prm_diag%drag_v_grid (jc,jb), &
                  &                                t_wgt )
              ENDDO  ! jc
              !$acc end parallel

            ENDIF  ! lcalc_extra_avg

          ENDIF  ! inwp_turb > 0


          IF ( lcall_phy_jg(itradheat) ) THEN
            !sum up for averaged fluxes
            !T.R.: this is not correct for output after 1st timestep,
            !e.g. dt_phy_jg(itradheat) may then be greater than p_sim_time
            !leading to wrong averaging.
!DIR$ IVDEP
            !$acc parallel default(present) if(lacc)
            !$acc loop gang vector
            DO jc = i_startidx, i_endidx

              ! time averaged shortwave net flux at surface
              prm_diag%swflxsfc_a(jc,jb) = time_avg(prm_diag%swflxsfc_a(jc,jb), &
                &                                   prm_diag%swflxsfc  (jc,jb), &
                &                                   t_wgt)

              ! time averaged clear-sky shortwave net flux at surface
              prm_diag%swflxclrsfc_a(jc,jb) = time_avg(prm_diag%swflxclrsfc_a(jc,jb), &
                &                                      prm_diag%swflxclr_sfc (jc,jb), &
                &                                      t_wgt)

              ! time averaged shortwave diffuse downward flux at surface
              prm_diag%asodifd_s (jc,jb) = time_avg(prm_diag%asodifd_s        (jc,jb), &
                &                                   prm_diag%swflx_dn_sfc_diff(jc,jb), &
                &                                   t_wgt)

              ! time averaged shortwave diffuse upward flux at surface
              prm_diag%asodifu_s (jc,jb) = time_avg(prm_diag%asodifu_s   (jc,jb), &
                &                                   prm_diag%swflx_up_sfc(jc,jb), &
                &                                   t_wgt)

              ! time averaged longwave net flux at surface
              prm_diag%lwflxsfc_a(jc,jb) = time_avg(prm_diag%lwflxsfc_a(jc,jb), &
                &                                   prm_diag%lwflxsfc  (jc,jb), &
                &                                   t_wgt)

              ! time averaged clear-sky longwave net flux at surface
              prm_diag%lwflxclrsfc_a(jc,jb) = time_avg(prm_diag%lwflxclrsfc_a(jc,jb), &
                &                                      prm_diag%lwflxclr_sfc (jc,jb), &
                &                                      t_wgt)

              ! time averaged longwave upward flux at surface
              prm_diag%athu_s    (jc,jb) = time_avg(prm_diag%athu_s      (jc,jb), &
                &                                   prm_diag%lwflx_up_sfc(jc,jb), &
                &                                   t_wgt)

              ! time averaged longwave downward flux at surface
              prm_diag%athd_s(jc,jb) = prm_diag%lwflxsfc_a(jc,jb) + prm_diag%athu_s(jc,jb)

              ! time averaged shortwave net flux at TOA
              prm_diag%swflxtoa_a(jc,jb) = time_avg(prm_diag%swflxtoa_a(jc,jb), &
                &                                   prm_diag%swflxtoa  (jc,jb), &
                &                                   t_wgt)

              ! time averaged longwave net flux at TOA
              prm_diag%lwflxtoa_a(jc,jb) = time_avg(prm_diag%lwflxtoa_a(jc,jb), &
                &                                   prm_diag%lwflxall(jc,1,jb), &
                &                                   t_wgt)

              ! time averaged top down solar radiation
              prm_diag%asod_t    (jc,jb) = time_avg(prm_diag%asod_t    (jc,jb), &
                &                                   prm_diag%flxdwswtoa(jc,jb), &
                &                                   t_wgt)

              ! time averaged solar upward flux at TOA
              prm_diag%asou_t(jc,jb) = prm_diag%asod_t(jc,jb) - prm_diag%swflxtoa_a(jc,jb)

              ! time averaged shortwave direct downward flux at surface
              prm_diag%asodird_s (jc,jb) = MAX(0._wp, prm_diag%swflxsfc_a(jc,jb) &
                &                        -            prm_diag%asodifd_s (jc,jb) &
                &                        +            prm_diag%asodifu_s (jc,jb) )

              ! downward solar radiation = sum of direct + diffuse
              prm_diag%asod_s(jc,jb) = prm_diag%asodifd_s(jc,jb) + prm_diag%asodird_s(jc,jb)

              ! time averaged downward photosynthetically active flux at surface
              prm_diag%aswflx_par_sfc(jc,jb) = time_avg(prm_diag%aswflx_par_sfc(jc,jb), &
                &                                       prm_diag%swflx_par_sfc(jc,jb),  &
                &                                       t_wgt)
            ENDDO
            !$acc end parallel

          ENDIF  ! lcall_phy_jg(itradheat)

        ELSEIF (.NOT. lflux_avg) THEN

#ifdef _OPENACC
          CALL finish('mo_nwp_diagnosis:','accumulating surface not available on GPU. ')
#endif

          IF (lcall_phy_jg(itturb)) THEN

!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              ! ATTENTION:
              ! the sign, in the output all fluxes must be positive downwards

              ! accumulated surface latent heat flux
              prm_diag%alhfl_s(jc,jb) =  prm_diag%alhfl_s(jc,jb)       &
                                 &  + prm_diag%lhfl_s(jc,jb)           & 
                                 &  * dt_phy_jg(itfastphy) 

              ! accumulated surface latent heat flux from bare soil
              prm_diag%alhfl_bs(jc,jb) =  prm_diag%alhfl_bs(jc,jb)     &
                                 &  + prm_diag%lhfl_bs(jc,jb)          & 
                                 &  * dt_phy_jg(itfastphy) 

              ! accumulated surface sensible heat flux
              prm_diag%ashfl_s(jc,jb) =  prm_diag%ashfl_s(jc,jb)       &
                                 &  + prm_diag%shfl_s(jc,jb)           & 
                                 &  * dt_phy_jg(itfastphy) 

              ! accumulated surface moisture flux
              prm_diag%aqhfl_s(jc,jb) =  prm_diag%aqhfl_s(jc,jb)       &
                                 &  + prm_diag%qhfl_s(jc,jb)           & 
                                 &  * dt_phy_jg(itfastphy)

              ! accumulated surface u-momentum flux turbulence
              prm_diag%aumfl_s(jc,jb) = prm_diag%aumfl_s(jc,jb)        &
                                &   + prm_diag%umfl_s(jc,jb)           &
                                &   * dt_phy_jg(itfastphy)

              ! accumulated surface v-momentum flux turbulence
              prm_diag%avmfl_s(jc,jb) = prm_diag%avmfl_s(jc,jb)        &
                                &   + prm_diag%vmfl_s(jc,jb)           &
                                &   * dt_phy_jg(itfastphy)
            ENDDO

            DO jk = 1, nlev_soil
!DIR$ IVDEP
              DO jc = i_startidx, i_endidx
                prm_diag%alhfl_pl(jc,jk,jb) =  prm_diag%alhfl_pl(jc,jk,jb)&
                                 &  + prm_diag%lhfl_pl(jc,jk,jb)          & 
                                 &  * dt_phy_jg(itfastphy) 
              ENDDO  ! jc
            ENDDO  ! jk


            IF (atm_phy_nwp_config(jg)%lcalc_extra_avg) THEN
!DIR$ IVDEP
              DO jc = i_startidx, i_endidx
                ! accumulated surface u-momentum flux SSO
                prm_diag%astr_u_sso(jc,jb) = prm_diag%astr_u_sso(jc,jb)     &
                                       &   + prm_diag%str_u_sso(jc,jb)      &
                                       &   * dt_phy_jg(itfastphy)

                ! accumulated surface v-momentum flux SSO
                prm_diag%astr_v_sso(jc,jb) = prm_diag%astr_v_sso(jc,jb)     &
                                       &   + prm_diag%str_v_sso(jc,jb)      &
                                       &   * dt_phy_jg(itfastphy)

                ! accumulated surface u-momentum flux resolved
                prm_diag%adrag_u_grid(jc,jb) = prm_diag%adrag_u_grid(jc,jb) &
                                         &   + prm_diag%drag_u_grid(jc,jb)  &
                                         &   * dt_phy_jg(itfastphy)

                ! accumulated surface v-momentum flux resolved
                prm_diag%adrag_v_grid(jc,jb) = prm_diag%adrag_v_grid(jc,jb) &
                                         &   + prm_diag%drag_v_grid(jc,jb)  &
                                         &   * dt_phy_jg(itfastphy)
              ENDDO  ! jc

            ENDIF  ! lcalc_extra_avg

          ENDIF  ! inwp_turb > 0


          IF ( lcall_phy_jg(itradheat) ) THEN
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              ! accumulated shortwave net flux at surface
              prm_diag%swflxsfc_a(jc,jb) = prm_diag%swflxsfc_a(jc,jb) &
             &                         + prm_diag%swflxsfc(jc,jb)     &
             &                         * dt_phy_jg(itfastphy)

              ! accumulated clear-sky shortwave net flux at surface
              prm_diag%swflxclrsfc_a(jc,jb) = prm_diag%swflxclrsfc_a(jc,jb) &
             &                           + prm_diag%swflxclr_sfc(jc,jb)     &
             &                           * dt_phy_jg(itfastphy)

              ! accumulated shortwave diffuse downward flux at surface
              prm_diag%asodifd_s (jc,jb) = prm_diag%asodifd_s        (jc,jb)  &
             &                           + prm_diag%swflx_dn_sfc_diff(jc,jb)  &
             &                           * dt_phy_jg(itfastphy)

              ! accumulated shortwave diffuse upward flux at surface
              prm_diag%asodifu_s (jc,jb) = prm_diag%asodifu_s   (jc,jb)  &
             &                           + prm_diag%swflx_up_sfc(jc,jb)  &
             &                           * dt_phy_jg(itfastphy)

              ! accumulated longwave net flux at surface
              prm_diag%lwflxsfc_a(jc,jb) = prm_diag%lwflxsfc_a(jc,jb) &
                                   &   + prm_diag%lwflxsfc(jc,jb)     &
                                   &   * dt_phy_jg(itfastphy)

              ! accumulated clear-sky longwave net flux at surface
              prm_diag%lwflxclrsfc_a(jc,jb) = prm_diag%lwflxclrsfc_a(jc,jb) &
                                     &   + prm_diag%lwflxclr_sfc(jc,jb)     &
                                     &   * dt_phy_jg(itfastphy)

              ! accumulated shortwave net flux at TOA
              prm_diag%swflxtoa_a(jc,jb) = prm_diag%swflxtoa_a(jc,jb) &
                                   &   + prm_diag%swflxtoa(jc,jb)     &
                                   &   * dt_phy_jg(itfastphy)

              ! accumulated longwave net flux at TOA
              prm_diag%lwflxtoa_a(jc,jb) = prm_diag%lwflxtoa_a(jc,jb) &
                                   &   + prm_diag%lwflxall(jc,1,jb)   &
                                   &   * dt_phy_jg(itfastphy)

              ! accumulated longwave upward flux at surface
              prm_diag%athu_s    (jc,jb) = prm_diag%athu_s(jc,jb) &
                &                  + prm_diag%lwflx_up_sfc(jc,jb) &
                                   &   * dt_phy_jg(itfastphy)

              ! accumulated longwave downward flux at surface
              prm_diag%athd_s(jc,jb) = prm_diag%lwflxsfc_a(jc,jb) + prm_diag%athu_s(jc,jb)

              ! accumulated top down solar radiation
              prm_diag%asod_t    (jc,jb) = prm_diag%asod_t(jc,jb)     &
                                   &   + prm_diag%flxdwswtoa(jc,jb)   &
                                   &   * dt_phy_jg(itfastphy)

              ! accumulated solar upward flux at TOA
              prm_diag%asou_t(jc,jb) = prm_diag%asod_t(jc,jb) - prm_diag%swflxtoa_a(jc,jb)

              ! accumulated shortwave direct downward flux at surface
              prm_diag%asodird_s (jc,jb) = MAX(0._wp, prm_diag%swflxsfc_a(jc,jb) &
                &                        -            prm_diag%asodifd_s (jc,jb) &
                &                        +            prm_diag%asodifu_s (jc,jb) )

              ! downward solar radiation = sum of direct + diffuse
              prm_diag%asod_s(jc,jb) = prm_diag%asodifd_s(jc,jb) + prm_diag%asodird_s(jc,jb)

              ! accumulated downward photosynthetically active flux at surface
              prm_diag%aswflx_par_sfc(jc,jb) = prm_diag%aswflx_par_sfc(jc,jb)  &
                &                            + prm_diag%swflx_par_sfc(jc,jb)   &
                &                            * dt_phy_jg(itfastphy)

            END DO
          ENDIF  ! lcall_phy_jg(itradheat)

        ENDIF  ! lflux_avg

      ENDDO ! nblks
!$OMP END DO NOWAIT

    END IF  ! p_sim_time

!$OMP END PARALLEL  

    IF (ltimer) CALL timer_stop(timer_nh_diagnostics)

  END SUBROUTINE nwp_statistics


  !>
  !! Computation of vertical integrals of moisture and cloud cover
  !!
  !!
  !! @par Revision History
  !! Separated from module nwp_statistics by Guenther Zaengl, DWD (2014-07-11)
  !! Includes calculation of high-, mid-, and low-level cloud cover, height
  !! of base and top of convection  by Helmut Frank, DWD (2013-01-17)
  !! Add height of 0 deg C level    by Helmut Frank, DWD (2013-03-11)
  !!
  !!
  SUBROUTINE calc_moist_integrals(pt_patch, p_metrics,        & !in
                                & pt_prog, pt_prog_rcf,       & !in
                                & ext_data, kstart_moist,     & !in
                                & ih_clch, ih_clcm,           & !in
                                & pt_diag, prm_diag,          & !inout
                                & lacc                        ) !in


    TYPE(t_patch),      INTENT(IN)   :: pt_patch    !<grid/patch info.
    TYPE(t_nh_diag),    INTENT(INOUT):: pt_diag     !<the diagnostic variables
    TYPE(t_nh_prog),    INTENT(IN)   :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog),    INTENT(IN)   :: pt_prog_rcf !<the prognostic variables (with red. calling frequency for tracers!)
    TYPE(t_nh_metrics), INTENT(in)   :: p_metrics
    TYPE(t_external_data),INTENT(IN) :: ext_data    !< external data
    TYPE(t_nwp_phy_diag), INTENT(inout):: prm_diag

    INTEGER,           INTENT(IN)  :: kstart_moist
    INTEGER,           INTENT(IN)  :: ih_clch, ih_clcm
    LOGICAL, OPTIONAL, INTENT(IN)  :: lacc

    ! Local
    INTEGER :: nlev                    !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices

    REAL(wp):: rhodz(nproma,pt_patch%nlev)   ! rho times delta z
    REAL(wp):: z_help

    INTEGER :: jc,jk,jb,jg,jk1      ! block index
    INTEGER :: jt               ! tracer loop index

    REAL(wp):: clearsky(nproma)
    REAL(wp):: ccmax, ccran, alpha(nproma,pt_patch%nlev), clcl_mod, clcm_mod, clct_fac
    LOGICAL :: lland
    LOGICAL :: lzacc

    REAL(wp), PARAMETER :: eps_clc = 1.e-7_wp

    INTEGER,  PARAMETER :: i_overlap = 2       ! 1: maximum-random overlap
                                               ! 2: generalized overlap (Hogan, Illingworth, 2000) 
    REAL(wp) :: zdecorr(pt_patch%nlev)         ! decorrelation length scale del(z0) 

  !-----------------------------------------------------------------

    IF(PRESENT(lacc)) THEN
      lzacc = lacc
    ELSE
      lzacc = .FALSE.
    ENDIF

    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev

    ! exclude nest boundary interpolation zone
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)

    ! set height-dependent decorrelation length scale
    !$acc data create (zdecorr, rhodz) if(lzacc)
    !$acc kernels if(lzacc)
    zdecorr(:) = 2000._wp
    !$acc end kernels
    !$acc parallel copyin (vct_a) default(present) if(lzacc)
    !$acc loop private(jk1,z_help)
    DO jk = nlev, 1, -1
      jk1 = jk + pt_patch%nshift_total
      z_help = 0.5_wp*(vct_a(jk1)+vct_a(jk1+1))
      IF (z_help < 3000._wp) THEN
        zdecorr(jk) = 800._wp + 0.4_wp*z_help
      ELSE
        EXIT
      ENDIF
    ENDDO
    !$acc end parallel

!$OMP PARALLEL
    IF ( atm_phy_nwp_config(jg)%lenabled(itccov) ) THEN

!$OMP DO PRIVATE(jc,jk,jb,z_help,i_startidx,i_endidx,clearsky,ccmax,ccran,alpha,clcl_mod,clcm_mod,clct_fac,lland)
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)


        ! if cloud cover is called, vertical integration of cloud content
        ! (for iqv, iqc, iqi)

        !$acc kernels if(lzacc)
        prm_diag%tot_cld_vi(i_startidx:i_endidx,jb,1:3) = 0.0_wp
        !$acc end kernels

        !$acc parallel default(present) if(lzacc)
        !$acc loop seq
        DO jk = kstart_moist, nlev
!DIR$ IVDEP
          !$acc loop gang vector private(z_help)
          DO jc = i_startidx, i_endidx

           ! (deep-atmosphere modification applied: height-dependence of grid cell volume)
           z_help = p_metrics%ddqz_z_full(jc,jk,jb) * pt_prog%rho(jc,jk,jb) & 
             &    * p_metrics%deepatmo_t1mc(jk,idamtr%t1mc%vol)

           ! TQV, TQC, TQI
           prm_diag%tot_cld_vi(jc, jb,iqv) = prm_diag%tot_cld_vi(jc, jb,iqv)    + &
                                             z_help * prm_diag%tot_cld(jc,jk,jb,iqv)
           prm_diag%tot_cld_vi(jc, jb,iqc) = prm_diag%tot_cld_vi(jc, jb,iqc)    + &
                                             z_help * prm_diag%tot_cld(jc,jk,jb,iqc)
           prm_diag%tot_cld_vi(jc, jb,iqi) = prm_diag%tot_cld_vi(jc, jb,iqi)    + &
                                             z_help * prm_diag%tot_cld(jc,jk,jb,iqi)
          ENDDO
        ENDDO
        !$acc end parallel


        ! cloud cover calculation
        ! note: the conversion into % is done within the internal output postprocessing

        SELECT CASE ( i_overlap )
 
        CASE ( 1 )      ! maximum-random overlap

          !$acc data create(clearsky)
          !$acc parallel default(present) if(lzacc)
          !$acc loop gang vector
          DO jc = i_startidx, i_endidx
            clearsky(jc) = 1._wp - prm_diag%clc(jc,kstart_moist,jb)
          ENDDO
          !$acc end parallel
          
          !$acc parallel default(present) if(lzacc)
          !$acc loop seq
          DO jk = kstart_moist+1, ih_clch
            !$acc loop gang vector
            DO jc = i_startidx, i_endidx
              clearsky(jc) = clearsky(jc)*    &
              &  ( 1._wp - MAX( prm_diag%clc(jc,jk  ,jb), prm_diag%clc(jc,jk-1,jb))) &
              & /( 1._wp - MIN( prm_diag%clc(jc,jk-1,jb), 1._wp - eps_clc) )
            ENDDO
          ENDDO
          !$acc end parallel
          
          ! store high-level clouds
          !$acc parallel default(present) if(lzacc)
          !$acc loop gang vector
          DO jc = i_startidx, i_endidx
            prm_diag%clch(jc,jb) = MAX( 0._wp, 1._wp - clearsky(jc) - eps_clc)
          ENDDO
          !$acc end parallel
          
          ! continue downward for total cloud cover
          !$acc parallel default(present) if(lzacc)
          !$acc loop seq
          DO jk = ih_clch+1, nlev
            !$acc loop gang vector
            DO jc = i_startidx, i_endidx
              clearsky(jc) = clearsky(jc)*    &
              &  ( 1._wp - MAX( prm_diag%clc(jc,jk,jb), prm_diag%clc(jc,jk-1,jb))) &
              & /( 1._wp - MIN( prm_diag%clc(jc,jk-1,jb), 1._wp - eps_clc) )
            ENDDO
          ENDDO
          !$acc end parallel
          
          ! store total cloud cover, start for mid-level clouds
!DIR$ IVDEP
          !$acc parallel default(present) if(lzacc)
          !$acc loop gang vector
          DO jc = i_startidx, i_endidx
            prm_diag%clct(jc,jb) = MAX( 0._wp, 1._wp - clearsky(jc) - eps_clc)
            clearsky(jc) = 1._wp - prm_diag%clc(jc,ih_clch+1,jb)
          ENDDO
          !$acc end parallel
          
          ! mid-level clouds
          !$acc parallel default(present) if(lzacc)
          !$acc loop seq
          DO jk = ih_clch+2, ih_clcm
            !$acc loop gang vector
            DO jc = i_startidx, i_endidx
              clearsky(jc) = clearsky(jc)*    &
              &  ( 1._wp - MAX( prm_diag%clc(jc,jk,jb), prm_diag%clc(jc,jk-1,jb))) &
              & /( 1._wp - MIN( prm_diag%clc(jc,jk-1,jb), 1._wp - eps_clc) )
            ENDDO
          ENDDO
          !$acc end parallel
          
          ! store mid-level cloud cover, start for low-level clouds
!DIR$ IVDEP
          !$acc parallel default(present) if(lzacc)
          !$acc loop gang vector
          DO jc = i_startidx, i_endidx
            prm_diag%clcm(jc,jb) = MAX( 0._wp, 1._wp - clearsky(jc) - eps_clc)
          
            clearsky(jc) = 1._wp - prm_diag%clc(jc,ih_clcm+1,jb)
          ENDDO
          !$acc end parallel
          
          ! continue downward for mid-level clouds
          !$acc parallel default(present) if(lzacc)
          !$acc loop seq
          DO jk = ih_clcm+2, nlev
            !$acc loop gang vector
            DO jc = i_startidx, i_endidx
              clearsky(jc) = clearsky(jc)*    &
              &  ( 1._wp - MAX( prm_diag%clc(jc,jk,jb), prm_diag%clc(jc,jk-1,jb))) &
              & /( 1._wp - MIN( prm_diag%clc(jc,jk-1,jb), 1._wp - eps_clc) )
            ENDDO
          ENDDO
          !$acc end parallel
          
          ! store low-level clouds
          !$acc parallel default(present) if(lzacc)
          !$acc loop gang vector
          DO jc = i_startidx, i_endidx
            prm_diag%clcl(jc,jb) = MAX( 0._wp, 1._wp - clearsky(jc) - eps_clc)
          ENDDO
          !$acc end parallel
          !$acc end data

        CASE ( 2 )      ! generalized overlap (Hogan, Illingworth, 2000)

          !$acc data create(alpha)
          !$acc parallel default(present) if(lzacc)
          !$acc loop gang vector
          DO jc = i_startidx, i_endidx
            prm_diag%clct(jc,jb) = prm_diag%clc(jc,kstart_moist,jb)
            prm_diag%clch(jc,jb) = prm_diag%clc(jc,kstart_moist,jb)
            prm_diag%clcm(jc,jb) = 0.0_wp 
            prm_diag%clcl(jc,jb) = 0.0_wp 
          ENDDO

!PREVENT_INCONSISTENT_IFORT_FMA
          !$acc loop seq
          DO jk = kstart_moist+1, nlev
!DIR$ IVDEP
            !$acc loop gang vector private(ccmax,ccran)
            DO jc = i_startidx, i_endidx   ! total cloud cover
              ccmax = MAX( prm_diag%clc(jc,jk,jb),  prm_diag%clct(jc,jb) )
              ccran =      prm_diag%clc(jc,jk,jb) + prm_diag%clct(jc,jb) - &
                       & ( prm_diag%clc(jc,jk,jb) * prm_diag%clct(jc,jb) )
              alpha(jc,jk) = MIN( EXP( - (p_metrics%z_mc(jc,jk-1,jb)-p_metrics%z_mc(jc,jk,jb)) / zdecorr(jk) ), &
                             prm_diag%clc(jc,jk-1,jb)/MAX(eps_clc,prm_diag%clc(jc,jk,jb)) )
              prm_diag%clct(jc,jb) = alpha(jc,jk) * ccmax + (1._wp-alpha(jc,jk)) * ccran
            ENDDO

            !$acc loop gang vector private(ccmax,ccran)
            DO jc = i_startidx, i_endidx
              IF (jk <= prm_diag%k400(jc,jb)-1) THEN    ! high cloud cover
                ccmax = MAX( prm_diag%clc(jc,jk,jb),  prm_diag%clch(jc,jb) )
                ccran =      prm_diag%clc(jc,jk,jb) + prm_diag%clch(jc,jb) - &
                         & ( prm_diag%clc(jc,jk,jb) * prm_diag%clch(jc,jb) )
                prm_diag%clch(jc,jb) = alpha(jc,jk) * ccmax + (1._wp-alpha(jc,jk)) * ccran
              
              ELSE IF (jk <= prm_diag%k800(jc,jb)-1) THEN  ! midlevel cloud cover
                ccmax = MAX( prm_diag%clc(jc,jk,jb),  prm_diag%clcm(jc,jb) )
                ccran =      prm_diag%clc(jc,jk,jb) + prm_diag%clcm(jc,jb) - &
                         & ( prm_diag%clc(jc,jk,jb) * prm_diag%clcm(jc,jb) )
                prm_diag%clcm(jc,jb) = alpha(jc,jk) * ccmax + (1._wp-alpha(jc,jk)) * ccran
              
              ELSE  ! low cloud cover
                ccmax = MAX( prm_diag%clc(jc,jk,jb),  prm_diag%clcl(jc,jb) )
                ccran =      prm_diag%clc(jc,jk,jb) + prm_diag%clcl(jc,jb) - &
                         & ( prm_diag%clc(jc,jk,jb) * prm_diag%clcl(jc,jb) )
                prm_diag%clcl(jc,jb) = alpha(jc,jk) * ccmax + (1._wp-alpha(jc,jk)) * ccran
              ENDIF
            ENDDO

          ENDDO
          !$acc end parallel

          ! calibration of layer-wise cloud cover fields
          IF (lcalib_clcov) THEN
            !$acc parallel default(present) if(lzacc)
            !$acc loop gang vector private(lland,clcl_mod,clcm_mod,clct_fac)
            DO jc = i_startidx, i_endidx
              lland = ext_data%atm%fr_land(jc,jb)+ext_data%atm%fr_lake(jc,jb) > 0._wp
              clcl_mod = MIN(MERGE(max_calibfac_clcl,1.5_wp,lland)*prm_diag%clcl(jc,jb), &
                EXP((1._wp+prm_diag%clcl(jc,jb))/2._wp*LOG(MAX(eps_clc,prm_diag%clcl(jc,jb)))))
              clcm_mod = MIN(MERGE(3._wp,1.5_wp,lland)*prm_diag%clcm(jc,jb), &
                EXP((1._wp+prm_diag%clcm(jc,jb))/2._wp*LOG(MAX(eps_clc,prm_diag%clcm(jc,jb)))))
              clct_fac = (clcl_mod+clcm_mod+prm_diag%clch(jc,jb)) /                        &
                MAX(eps_clc,prm_diag%clcl(jc,jb)+prm_diag%clcm(jc,jb)+prm_diag%clch(jc,jb))
              clct_fac = MIN(clct_fac, SQRT(1._wp/MAX(0.05_wp,prm_diag%clct(jc,jb))) )
              prm_diag%clct(jc,jb) = clct_fac*prm_diag%clct(jc,jb)
              prm_diag%clcm(jc,jb) = clcm_mod
              prm_diag%clcl(jc,jb) = clcl_mod
            ENDDO
            !$acc end parallel
          ENDIF
          !$acc end data

        END SELECT

      ENDDO ! nblks
!$OMP END DO NOWAIT

    END IF !cloud cover



    ! Calculate vertically integrated values of the grid-scale tracers 
    ! Vertical integrals are computed for all mass concentrations.
    ! Number concentrations are skipped.
    !
!$OMP DO PRIVATE(jt,jc,jk,jb,i_startidx,i_endidx,rhodz) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      ! pre-computation of rho * \Delta z
      ! (deep-atmosphere modification applied: height-dependence of grid cell volume)
      !$acc parallel default(present) if(lzacc)
      !$acc loop collapse(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx 
          rhodz(jc,jk) = p_metrics%ddqz_z_full(jc,jk,jb) * pt_prog%rho(jc,jk,jb) & 
            &          * p_metrics%deepatmo_t1mc(jk,idamtr%t1mc%vol)  
        ENDDO
      ENDDO
      !$acc end parallel

      DO jt = 1, iqm_max
        !$acc kernels default(present) if(lzacc)
        pt_diag%tracer_vi(i_startidx:i_endidx,jb,jt) = 0.0_wp
        !$acc end kernels

        !$acc parallel default(present) if(lzacc)
        !$acc loop seq
        DO jk = advection_config(jg)%iadv_slev(jt), nlev

!DIR$ IVDEP
          !$acc loop gang vector
          DO jc = i_startidx, i_endidx 

            pt_diag%tracer_vi(jc,jb,jt) = pt_diag%tracer_vi(jc,jb,jt)   &
              &                + rhodz(jc,jk) * pt_prog_rcf%tracer(jc,jk,jb,jt) 

          ENDDO  ! jc
        ENDDO  ! jk
        !$acc end parallel
      ENDDO  ! jt

    ENDDO ! nblks   
!$OMP END DO
!$OMP END PARALLEL  
    
    !$acc end data

  END SUBROUTINE calc_moist_integrals




  !>
  !! Diagnostics which are only required for output
  !!
  !! Diagnostics which are only required for output. Gathers  
  !! computations which are purely diagnostic and only required for 
  !! (meteogram) output.
  !!
  !! Available diagnostics:
  !! - height of convection base and top: hbas_con, htop_con
  !! - height of the top of dry convection: htop_dc
  !! - height of 0 deg C level: hzerocl
  !! - height of snow fall limit above MSL
  !! - CLDEPTH: modified cloud depth for media
  !! - CLCT_MOD: modified total cloud cover (between 0 and 1) 
  !! - t_ice is filled with t_so(0) for non-ice points (h_ice=0)
  !! - instantaneous 10m wind speed (resolved scales)
  !!
  !! @par Revision History
  !! Add calculation of high-, mid-, and low-level cloud cover, height
  !! of base and top of convection  by Helmut Frank, DWD (2013-01-17)
  !! Add height of 0 deg C level    by Helmut Frank, DWD (2013-03-11)
  !! Modification by Daniel Reinert (2014-02-27)
  !! - separated all those diagnostics which are only required at output 
  !!   times and moved them into a separate routine.
  !!
  !!
  SUBROUTINE nwp_diag_for_output(mtime_current,           & !in
                            & kstart_moist,               & !in
                            & ih_clch, ih_clcm,           & !in
                            & phy_params,                 & !in
                            & pt_patch, p_metrics,        & !in
                            & pt_prog, pt_prog_rcf,       & !in
                            & pt_diag,                    & !in
                            & lnd_diag,                   & !in
                            & p_prog_lnd_now,             & !in
                            & p_prog_wtr_now,             & !in
                            & ext_data,                   & !in 
                            & prm_diag                    ) !inout
              
    TYPE(datetime),   POINTER     :: mtime_current     ! current datetime (mtime)
    INTEGER,         INTENT(IN)   :: kstart_moist
    INTEGER,         INTENT(IN)   :: ih_clch, ih_clcm

    TYPE(t_phy_params),INTENT(IN) :: phy_params
    TYPE(t_patch),   INTENT(IN)   :: pt_patch    !<grid/patch info.
    TYPE(t_nh_prog), INTENT(IN)   :: pt_prog     !<the prognostic variables 
    TYPE(t_nh_prog), INTENT(IN)   :: pt_prog_rcf !<the prognostic variables (with
                                                 !< red. calling frequency for tracers!
    TYPE(t_nh_metrics)  ,INTENT(IN) :: p_metrics

    TYPE(t_nh_diag),     INTENT(INOUT):: pt_diag     ! the diagnostic variables
    TYPE(t_lnd_diag),    INTENT(IN)   :: lnd_diag    ! land diag state
    TYPE(t_lnd_prog),    INTENT(IN)   :: p_prog_lnd_now ! land prognostic state (now)
    TYPE(t_wtr_prog),    INTENT(INOUT):: p_prog_wtr_now ! water prognostic state (now)
    TYPE(t_external_data),INTENT(IN)  :: ext_data       !< external data
    TYPE(t_nwp_phy_diag),INTENT(INOUT):: prm_diag

    ! Local
    INTEGER :: jc,jk,jb,jg             !< loop index
    INTEGER :: nlev, nlevp1            !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    REAL(wp):: zbuoy, zqsat, zcond
    INTEGER :: mtop_min
    REAL(wp):: ztp(nproma), zqp(nproma)
    INTEGER :: mlab(nproma)

    REAL(wp), PARAMETER :: grav_o_cpd = grav/cpd

    REAL(wp), PARAMETER :: zundef = -999._wp   ! undefined value for 0 deg C level

    TYPE(timeDelta), POINTER :: time_diff

  !-----------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_nh_diagnostics)

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1    


    ! exclude nest boundary interpolation zone
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


    ! minimum top index for dry convection
    mtop_min = (ih_clch+ih_clcm)/2    

    CALL calc_moist_integrals(pt_patch, p_metrics,        & !in
                            & pt_prog, pt_prog_rcf,       & !in
                            & ext_data, kstart_moist,     & !in
                            & ih_clch, ih_clcm,           & !in
                            & pt_diag, prm_diag           ) !inout

    ! time difference since last call of ww_diagnostics
    time_diff => newTimedelta("PT0S")
    time_diff =  getTimeDeltaFromDateTime(mtime_current, ww_datetime(jg))

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,mlab,ztp,zqp,zbuoy,zqsat,zcond) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      !
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      !
      ! Calculation of grid scale (gsp) and total (gsp+con) instantaneous precipitation rates:
      !
      SELECT CASE (atm_phy_nwp_config(jg)%inwp_gscp)
      CASE(4,5,6,7)
        DO jc =  i_startidx, i_endidx
          prm_diag%prec_gsp_rate(jc,jb) = prm_diag%rain_gsp_rate(jc,jb)  &
               &                        + prm_diag%ice_gsp_rate(jc,jb)   &
               &                        + prm_diag%snow_gsp_rate(jc,jb)  &
               &                        + prm_diag%hail_gsp_rate(jc,jb)  &
               &                        + prm_diag%graupel_gsp_rate(jc,jb)
          prm_diag%tot_prec_rate(jc,jb) = prm_diag%prec_gsp_rate(jc,jb)
        ENDDO
      CASE(2)
        DO jc =  i_startidx, i_endidx
          prm_diag%prec_gsp_rate(jc,jb) = prm_diag%rain_gsp_rate(jc,jb)  &
               ! not sure what to do with ice. To be consistent to prm_diag%prec_gsp, where ice is neglected
               ! because it predominantly is made of blowing snow, we neglect it also here:
!               &                        + prm_diag%ice_gsp_rate(jc,jb)   &
               &                        + prm_diag%snow_gsp_rate(jc,jb)  &
               &                        + prm_diag%graupel_gsp_rate(jc,jb)
          prm_diag%tot_prec_rate(jc,jb) = prm_diag%prec_gsp_rate(jc,jb)
        ENDDO
      CASE (1)
        DO jc =  i_startidx, i_endidx
          prm_diag%prec_gsp_rate(jc,jb) = prm_diag%rain_gsp_rate(jc,jb)  &
               ! not sure what to do with ice. To be consistent to prm_diag%prec_gsp, where ice is neglected
               ! because it predominantly is made of blowing snow, we neglect it also here:
!               &                        + prm_diag%ice_gsp_rate(jc,jb)   &
               &                        + prm_diag%snow_gsp_rate(jc,jb)
          prm_diag%tot_prec_rate(jc,jb) = prm_diag%prec_gsp_rate(jc,jb)
        ENDDO
      CASE (9)
        DO jc =  i_startidx, i_endidx
          prm_diag%prec_gsp_rate(jc,jb) = prm_diag%rain_gsp_rate(jc,jb)
          prm_diag%tot_prec_rate(jc,jb) = prm_diag%prec_gsp_rate(jc,jb)
        ENDDO
      CASE default
        DO jc =  i_startidx, i_endidx
          prm_diag%prec_gsp_rate(jc,jb) = 0.0_wp
          prm_diag%tot_prec_rate(jc,jb) = 0.0_wp
        ENDDO
      END SELECT
      !
      ! Add convective contributions to the total precipitation rate:
      !
      IF (atm_phy_nwp_config(jg)%inwp_convection > 0) THEN
        DO jc = i_startidx, i_endidx
          prm_diag%tot_prec_rate(jc,jb) = prm_diag%tot_prec_rate(jc,jb) + prm_diag%rain_con_rate(jc,jb) + &
               &                          prm_diag%snow_con_rate(jc,jb)
        ENDDO  ! jc
      END IF

   
      IF (atm_phy_nwp_config(jg)%lenabled(itconv))THEN !convection parameterization switched on
        !
        ! height of convection base and top, hbas_con, htop_con
        ! 
        DO jc = i_startidx, i_endidx
          IF ( prm_diag%locum(jc,jb) ) THEN
            prm_diag%hbas_con(jc,jb) = p_metrics%z_ifc( jc, prm_diag%mbas_con(jc,jb), jb)
            prm_diag%htop_con(jc,jb) = p_metrics%z_ifc( jc, prm_diag%mtop_con(jc,jb), jb)
!           Do not allow diagnostic depth of convection to be thinner than 100m or one model layer
            IF ( prm_diag%htop_con(jc,jb) - prm_diag%hbas_con(jc,jb) < 100._wp ) THEN
              prm_diag%hbas_con(jc,jb) = -500._wp
              prm_diag%htop_con(jc,jb) = -500._wp
            END IF
          ELSE
            prm_diag%hbas_con(jc,jb) = -500._wp
            prm_diag%htop_con(jc,jb) = -500._wp
          END IF
        ENDDO  ! jc


        !
        ! height of the top of dry convection
        !
        DO jc = i_startidx, i_endidx 
          prm_diag%htop_dc(jc,jb) = zundef
          mlab(jc) = 1
          ztp (jc) = pt_diag%temp(jc,nlev,jb) + 0.25_wp
          zqp (jc) = pt_prog_rcf%tracer(jc,nlev,jb,iqv)
        ENDDO

        DO jk = nlev-1, mtop_min, -1
          DO jc = i_startidx, i_endidx 
            IF ( mlab(jc) == 1) THEN
              ztp(jc) = ztp(jc)  - grav_o_cpd*( p_metrics%z_mc(jc,jk,jb)    &
             &                                 -p_metrics%z_mc(jc,jk+1,jb) )
              zbuoy = ztp(jc)*( 1._wp + vtmpc1*zqp(jc) ) - pt_diag%tempv(jc,jk,jb)
              zqsat = spec_humi( sat_pres_water(ztp(jc)), pt_diag%pres(jc,jk,jb) )
              zcond = zqp(jc) - zqsat

              IF ( zcond < 0._wp .AND. zbuoy > 0._wp) THEN
                prm_diag%htop_dc(jc,jb) = p_metrics%z_ifc(jc,jk,jb)
              ELSE
                mlab(jc) = 0
              END IF
            END IF
          ENDDO
        ENDDO

        DO jc = i_startidx, i_endidx 
          IF ( prm_diag%htop_dc(jc,jb) > zundef) THEN
            prm_diag%htop_dc(jc,jb) = MIN( prm_diag%htop_dc(jc,jb),        &
           &                p_metrics%z_ifc(jc,nlevp1,jb) + 3000._wp )
            IF ( prm_diag%hbas_con(jc,jb) /= -500._wp) THEN
              prm_diag%htop_dc(jc,jb) = MIN( prm_diag%htop_dc(jc,jb),      &
             &                               prm_diag%hbas_con(jc,jb) )
            END IF
          ELSE
            prm_diag%htop_dc(jc,jb) = MIN( 0._wp, p_metrics%z_ifc(jc,nlevp1,jb) )
          END IF
        ENDDO

      END IF !convection parameterization on


      !
      ! height of 0 deg C level "hzerocl". Take uppermost freezing level in case of multiple 
      ! occurrences, use orography height if temperature is below freezing in all levels
      !
      ! Initialization with orography height
      prm_diag%hzerocl(i_startidx:i_endidx,jb) = p_metrics%z_ifc(i_startidx:i_endidx,nlevp1,jb)

      DO jk = kstart_moist+1, nlev
        DO jc = i_startidx, i_endidx 
          IF ( prm_diag%hzerocl(jc,jb) > p_metrics%z_ifc(jc,nlevp1,jb)) THEN ! freezing level found
            CYCLE
          ELSE IF (pt_diag%temp(jc,jk-1,jb) < tmelt .AND. pt_diag%temp(jc,jk,jb) >= tmelt) THEN
            prm_diag%hzerocl(jc,jb) = p_metrics%z_mc(jc,jk-1,jb) -            &
           &      ( p_metrics%z_mc(jc,jk-1,jb) - p_metrics%z_mc(jc,jk,jb) )*  &
           &      (    pt_diag%temp(jc,jk-1,jb) - tmelt ) /                   &
           &      (    pt_diag%temp(jc,jk-1,jb) - pt_diag%temp(jc,jk,jb) )
          END IF
        ENDDO
      ENDDO


      !
      !  Height of snow fall limit above MSL (snow line)
      !
      CALL calsnowlmt ( snowlmt = prm_diag%snowlmt(:,jb)        , & !inout
        &               temp    = pt_diag%temp(:,:,jb)          , & !in
        &               pres    = pt_diag%pres(:,:,jb)          , & !in
        &               qv      = pt_prog_rcf%tracer(:,:,jb,iqv), & !in
        &               hhl     = p_metrics%z_ifc(:,:,jb)       , & !in
        &               hhlr    = vct_a(pt_patch%nshift_total+1:),& !in
        &               istart  = i_startidx                    , & !in
        &               iend    = i_endidx                      , & !in
        &               wbl     = 1.3_wp )


      ! Fill t_ice with t_so(1) for ice-free points (h_ice<=0)
      ! This was demanded by FE14 (surface analysis)
      !
      ! Note, that t_ice contains ice temperature information from 
      ! the sea ice model as well as the lake model.
      !
      ! Furthermore, note that filling t_ice with t_so(1) only makes 
      ! sense when running without tiles. When using tiles, t_ice contains 
      ! the temperatures of sea-ice tiles and frozen lake tiles. Mixing this field 
      ! with aggeregated t_so values makes no sense from my point of view.
      IF ( (ntiles_total == 1) .AND. (atm_phy_nwp_config(jg)%inwp_surface > 0)) THEN
        DO jc = i_startidx, i_endidx 
          p_prog_wtr_now%t_ice(jc,jb) = MERGE(                               &
            &                           lnd_diag%t_so(jc,1,jb),              &
            &                           p_prog_wtr_now%t_ice(jc,jb),         &
            &                           p_prog_wtr_now%h_ice(jc,jb) <= 0._wp )
        ENDDO  !jc
      ENDIF



      ! Compute resolved surface drag: ps * del(orog)
 
      DO jc = i_startidx, i_endidx
         prm_diag%drag_u_grid(jc,jb) = pt_diag%pres_ifc(jc,nlevp1,jb) * ext_data%atm%grad_topo(1,jc,jb)
         prm_diag%drag_v_grid(jc,jb) = pt_diag%pres_ifc(jc,nlevp1,jb) * ext_data%atm%grad_topo(2,jc,jb)
      ENDDO


      IF (atm_phy_nwp_config(jg)%inwp_gscp > 0 ) THEN

        CALL ww_diagnostics( nproma, nlev, nlevp1, i_startidx, i_endidx, jg,             &
            &                pt_diag%temp(:,:,jb), pt_prog_rcf%tracer(:,:,jb,iqv),       &
            &                pt_prog_rcf%tracer(:,:,jb,iqc),                             &
            &                pt_diag%u   (:,:,jb), pt_diag%v         (:,:,jb),           &
            &                prm_diag%clc(:,:,jb),                                       &
            &                pt_diag%pres(:,:,jb), pt_diag%pres_ifc  (:,:,jb),           &
            &                prm_diag%t_2m     (:,jb), prm_diag%td_2m   (:,jb),          &
            &                p_prog_lnd_now%t_g(:,jb),                                   &
            &                prm_diag%clct     (:,jb), prm_diag%clcm    (:,jb),          &
            &                prm_diag%u_10m    (:,jb), prm_diag%v_10m   (:,jb),          &
            &                prm_diag%rain_gsp0(:,jb), prm_diag%rain_gsp(:,jb),          &
            &                prm_diag%rain_con0(:,jb), prm_diag%rain_con(:,jb),          &
            &                prm_diag%snow_gsp0(:,jb), prm_diag%snow_gsp(:,jb),          &
            &                prm_diag%snow_con0(:,jb), prm_diag%snow_con(:,jb),          &
            &                prm_diag%mbas_con (:,jb), prm_diag%mtop_con(:,jb),          &
            &                time_diff, prm_diag%iww(:,jb) )
!       Save precipitation and time until next call of ww_diagnostics
        DO jc = i_startidx, i_endidx
          prm_diag%rain_gsp0(jc,jb) = prm_diag%rain_gsp(jc,jb)
          prm_diag%rain_con0(jc,jb) = prm_diag%rain_con(jc,jb)
          prm_diag%snow_gsp0(jc,jb) = prm_diag%snow_gsp(jc,jb)
          prm_diag%snow_con0(jc,jb) = prm_diag%snow_con(jc,jb)
        ENDDO
      ENDIF

      !
      !  CAPE and CIN of mean surface layer parcel
      !
      !  start level (kmoist) is limited to pressure heights above p=60hPa, 
      !  in order to avoid unphysically low test parcel temperature.
      !  Otherwise computation crashes in sat_pres_water  
      CALL cal_cape_cin( i_startidx, i_endidx,                     &
        &                kmoist  = MAX(kstart_moist,phy_params%k060), & !in
        &                te      = pt_diag%temp(:,:,jb)          , & !in
        &                qve     = pt_prog_rcf%tracer(:,:,jb,iqv), & !in
        &                prs     = pt_diag%pres(:,:,jb)          , & !in
        &                hhl     = p_metrics%z_ifc(:,:,jb)       , & !in
        &                cape_ml = prm_diag%cape_ml(:,jb)        , & !in
        &                cin_ml  = prm_diag%cin_ml(:,jb) )

    ENDDO  ! jb
!$OMP END DO

!$OMP END PARALLEL  
    ww_datetime(jg) = time_config%tc_current_date

    ! compute modified cloud parameters for TV presentation
    CALL calcmod( pt_patch, pt_diag, prm_diag )

    IF (ltimer) CALL timer_stop(timer_nh_diagnostics)
    CALL deallocateTimedelta(time_diff)

  END SUBROUTINE nwp_diag_for_output


  !-------------------------------------------------------------------------
  !>
  !! Calculates modified cloud parameters for TV presentation.
  !!
  !! Calculates modified cloud parameters for TV presentation, namely
  !! CLDEPTH : modified cloud depth (scaled between 0 and 1)
  !! CLCT_MOD: modified total cloud cover (between 0 and 1) 
  !! 
  !! Both quantities are derived from the cloud cover "clc" on each
  !! model layer by neglecting cirrus clouds if they are the only
  !! clouds at this grid point. The reason for this treatment is that
  !! the general public does not regard transparent cirrus clouds as
  !! "real" clouds.
  !!
  !! @par Revision History
  !! Developed by D. Majewski, DWD (2004).
  !! Modification by Daniel Reinert, DWD (2014-09-10)
  !! - Adapted to and implemented into ICON
  !!
  !!
  SUBROUTINE calcmod( pt_patch, pt_diag, prm_diag )    
              
    TYPE(t_patch)       ,INTENT(IN)   :: pt_patch  !<grid/patch info.
    TYPE(t_nh_diag)     ,INTENT(IN)   :: pt_diag
    TYPE(t_nwp_phy_diag),INTENT(INOUT):: prm_diag

    ! Local
    INTEGER :: jc,jk,jb                !< loop index
    INTEGER :: nlev                    !< number of full levels
    INTEGER :: nlevp1                  !< number of half levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    INTEGER :: jk_top, jk_bot         !< start and end level for linear interpolation
    INTEGER :: iclbas(nproma)         !< base level of significant cloudiness
    REAL(wp):: p_clbas(nproma)        !< pressure (Pa) at base of significant cloudiness
    REAL(wp):: zred                   !< cloud cover reduction factor

    REAL(wp), PARAMETER :: p_clbas_min = 200.0E2_wp ! lower bound for reduction factor
    REAL(wp), PARAMETER :: p_clbas_max = 600.0E2_wp ! upper bound for reduction factor
    REAL(wp), PARAMETER :: clct_min    = 0.5_wp     ! threshold for significant cloudiness

  !--------------------------------------------------------------------

    i_nchdom  = MAX(1,pt_patch%n_childdom)

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    ! exclude nest boundary interpolation zone
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,jk_top,jk_bot,iclbas,p_clbas,zred) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      !
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      ! 
      ! modified cloud depth for media
      !
      ! calculation of the normalized cloud depth 'cldepth' as a modified cloud parameter 
      ! for TV presentation. The vertical integral of cloud cover in pressure units is 
      ! normalized by 700hPa. Thus, cldepth=1 for a cloud extending vertically over a 
      ! range of 700 hPa. Only used for visualization purpose (i.e. gray-scale pictures)
      !
      prm_diag%cldepth(i_startidx:i_endidx,jb) = 0._wp
      !
      DO jk=1, nlev
        DO jc = i_startidx, i_endidx 
           prm_diag%cldepth(jc,jb) = prm_diag%cldepth(jc,jb)   & 
             &                     + prm_diag%clc(jc,jk,jb) * pt_diag%dpres_mc(jc,jk,jb)
        ENDDO  ! jc
      ENDDO  ! jk
      !
      ! Normalize:
      DO jc = i_startidx, i_endidx 
        prm_diag%cldepth(jc,jb) = MIN(1._wp,prm_diag%cldepth(jc,jb)/700.E2_wp)
      ENDDO


      !
      ! modified total cloud cover for media
      !
      ! do not take high clouds into account, if they are the only clouds present 
      ! at this grid point. The computation of the cloud cover uses maximum overlapping. 
      !
      ! initialize
      prm_diag%clct_mod(i_startidx:i_endidx,jb) = 0._wp  ! modified cloud cover
      p_clbas(i_startidx:i_endidx)              = 0._wp  ! pressure at base of significant cloudiness
      iclbas(i_startidx:i_endidx)               = 1      ! level at base of significant cloudiness

      ! Determine base level of significant cloudiness
      ! Cloudiness is assumed to be significant, if clc>clct_min (=0.5)
      ! If there is no significant cloudiness within a column: iclbas = 1
      DO jk=1, nlev
        DO jc = i_startidx, i_endidx 
          IF ( prm_diag%clc(jc,jk,jb) >= clct_min ) THEN
            ! half-level index at base of significant cloudiness
            iclbas(jc) = jk+1
          ENDIF
        ENDDO  ! jc
      ENDDO  ! jk


      ! compute pressure at base of significant cloudiness, i.e. pressure at 
      ! height where clct_min is reached (linear interpolation is performed 
      ! between pressure at upper and lower main level)
      !
      ! setup for linear interpolation
      ! |
      ! |--------------------------------------
      ! |
      ! |               X  pres(jk_top); clc >= 0.5
      ! |
      ! |------------ iclbas -------------------
      ! |
      ! |               X  pres(jk_bot); clc < 0.5
      ! |
      ! |--------------------------------------
      !
      DO jc = i_startidx, i_endidx
        IF (iclbas(jc) == 1) THEN     ! no cloud at this grid point
          p_clbas(jc) = 0._wp
        ELSE IF (iclbas(jc) == nlevp1) THEN  ! set to surface pressure
          p_clbas(jc) = pt_diag%pres_ifc(jc,nlevp1,jb)
        ELSE                          ! Interpolate base pressure
          jk_top = iclbas(jc) - 1
          jk_bot = iclbas(jc)
          p_clbas(jc) = (prm_diag%clc(jc,jk_top,jb) - clct_min) &
            &         / MAX(0.001_wp,prm_diag%clc(jc,jk_top,jb) - prm_diag%clc(jc,jk_bot,jb)) &
            &         * (pt_diag%pres(jc,jk_bot,jb) - pt_diag%pres(jc,jk_top,jb))  &
            &         + pt_diag%pres(jc,jk_top,jb)
        ENDIF
      ENDDO

      ! compute cloud cover using maximum overlapping
      DO jk = 1,nlev
        DO jc = i_startidx, i_endidx
          prm_diag%clct_mod(jc,jb) = MAX (prm_diag%clct_mod(jc,jb), prm_diag%clc(jc,jk,jb))
        ENDDO
      ENDDO  ! jk

      ! Compute the modified total cloud cover; do not take the high clouds
      ! into account if they are the only clouds present at this grid
      ! point
      !
      !
      ! Profile of the reduction factor
      ! |
      ! |      zred = 0
      ! |
      ! |----------- 200 hPa -------------
      ! |
      ! |      zred = COS (0.5*pi*(...)/(...))
      ! |
      ! |----------- 600 hPa -------------
      ! |
      ! |      zred = 1
      ! |__________________________________>

      DO jc = i_startidx, i_endidx
        zred = 1._wp
        IF (p_clbas(jc) < p_clbas_min) THEN
          zred  = 0._wp
        ELSE IF (p_clbas(jc) < p_clbas_max) THEN
          zred  = MAX (0._wp, COS(0.5_wp*pi/(p_clbas_min - p_clbas_max)  &
            &     * (p_clbas(jc) - p_clbas_max)) )
        ENDIF
        prm_diag%clct_mod(jc,jb) = zred * prm_diag%clct_mod(jc,jb)
      ENDDO  ! jc


    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE calcmod



  !>
  !! Subroutine collecting the calls for computing optional diagnostic output variables
  !! needed primarily for convection-permitting model configurations
  !!
  !! Moved from nh_stepping for better code structure
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2020-02-14)
  !!
  !!
  SUBROUTINE nwp_opt_diagnostics(p_patch, p_patch_lp, p_int_lp, p_nh, p_int, prm_diag, &
     l_output, nnow, nnow_rcf, var_in_output, &
     lpi_max_Event, celltracks_Event, dbz_Event, mtime_current,  plus_slack)

    TYPE(t_patch)       ,INTENT(IN)   :: p_patch(:), p_patch_lp(:)  ! patches and their local parents
    TYPE(t_int_state)   ,INTENT(IN)   :: p_int_lp(:)                ! interpolation state for local parents
    TYPE(t_nh_state)    ,INTENT(INOUT):: p_nh(:)                    ! nonhydro state
    TYPE(t_int_state)   ,INTENT(IN)   :: p_int(:)                   ! interpolation state
    TYPE(t_nwp_phy_diag),INTENT(INOUT):: prm_diag(:)                ! physics diagnostics

    TYPE(event),     POINTER, INTENT(INOUT) :: lpi_max_Event, celltracks_Event, dbz_Event
    TYPE(datetime),  POINTER, INTENT(IN   ) :: mtime_current  !< current_datetime
    TYPE(timedelta), POINTER, INTENT(IN   ) :: plus_slack

    TYPE(t_var_in_output),    INTENT(IN   ) :: var_in_output(:)

    LOGICAL, INTENT(IN) :: l_output(:)
    INTEGER, INTENT(IN) :: nnow(:), nnow_rcf(:)

    LOGICAL :: l_active(3), l_lpimax_event_active, l_celltracks_event_active, l_dbz_event_active, &
               l_need_dbz3d, l_need_temp, l_need_pres
    INTEGER :: jg, k

    l_active(1) = is_event_active(lpi_max_Event,    mtime_current, proc0_offloading, plus_slack, opt_lasync=.TRUE.)
    l_active(2) = is_event_active(celltracks_Event, mtime_current, proc0_offloading, plus_slack, opt_lasync=.TRUE.)
    l_active(3) = is_event_active(dbz_Event,        mtime_current, proc0_offloading, plus_slack, opt_lasync=.TRUE.)

    ! In NEC hybrid mode, mtime is called on p_io only, so result needs to be broadcasted
    IF (proc0_offloading)  CALL p_bcast(l_active, p_io, p_comm_work)

    l_lpimax_event_active     = l_active(1)
    l_celltracks_event_active = l_active(2)
    l_dbz_event_active        = l_active(3)

    IF (ltimer) CALL timer_start(timer_nh_diagnostics)

    DO jg = 1, n_dom
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE 

      ! First check if pressure and/or temperature need to be diagnosed before calling the routines below

      l_need_dbz3d = l_output(jg) .AND. (             &
           &         var_in_output(jg)%dbz       .OR. &
           &         var_in_output(jg)%dbz850    .OR. &
           &         var_in_output(jg)%dbzcmax   .OR. &
           &         var_in_output(jg)%dbzctmax  .OR. &
           &         var_in_output(jg)%echotop   .OR. &
           &         var_in_output(jg)%echotopinm     )
      
      l_need_dbz3d = l_need_dbz3d .OR. (             &
           &         l_dbz_event_active .AND. (      &
           &         var_in_output(jg)%dbzctmax .OR. &
           &         var_in_output(jg)%echotop  .OR. &
           &         var_in_output(jg)%echotopinm  ) )

      ! Remark: at output dates, temperature and pressure have already been diagnosed
      l_need_temp = var_in_output(jg)%lpi_max     .AND. l_lpimax_event_active     .OR. &
                    var_in_output(jg)%tcond10_max .AND. l_celltracks_event_active .OR. &
                    l_need_dbz3d .AND. .NOT. l_output(jg)
      l_need_pres = l_need_dbz3d .AND. .NOT. l_output(jg)

      IF (l_need_temp .OR. l_need_pres) THEN
        CALL diagnose_pres_temp(p_nh(jg)%metrics, p_nh(jg)%prog(nnow(jg)), p_nh(jg)%prog(nnow_rcf(jg)),         &
                                p_nh(jg)%diag, p_patch(jg), opt_calc_temp=l_need_temp, opt_calc_pres=l_need_pres)
      ENDIF


      ! maximization of LPI_MAX (LPI max. during the time interval "celltracks_interval") if required
      IF ( var_in_output(jg)%lpi_max .AND. (l_output(jg) .OR. l_lpimax_event_active ) ) THEN
        IF ( jg >= n_dom_start+1 ) THEN
          ! p_patch_local_parent(jg) seems to exist
          CALL maximize_field_lpi( p_patch(jg), jg, p_patch_lp(jg), p_int_lp(jg), p_nh(jg)%metrics,      &
               &                   p_nh(jg)%prog(nnow(jg)), p_nh(jg)%prog(nnow_rcf(jg)), p_nh(jg)%diag,  &
               &                   prm_diag(jg)%lpi_max  )
        ELSE
          CALL message( "perform_nh_timeloop", "WARNING: LPI_MAX cannot be computed since no reduced grid is available" )
        END IF
      END IF

      ! update of TCOND_MAX (total column-integrated condensate, max. during the time interval "celltracks_interval") if required
      IF ( ( var_in_output(jg)%tcond_max .OR. var_in_output(jg)%tcond10_max ) .AND.  &
           ( l_output(jg) .OR. l_celltracks_event_active) ) THEN
        CALL compute_field_tcond_max( p_patch(jg), jg, p_nh(jg)%metrics,   &
             &                        p_nh(jg)%prog(nnow(jg)),  p_nh(jg)%prog(nnow_rcf(jg)), p_nh(jg)%diag, &
             &                        var_in_output(jg)%tcond_max, var_in_output(jg)%tcond10_max,     &
             &                        prm_diag(jg)%tcond_max, prm_diag(jg)%tcond10_max  )
      END IF

      ! update vorticity for the calculation of uh_max / vorw_ctmax;
      ! otherwise, diag%vor is only diagnosed every output time step
      IF ( (ANY(luh_max_out(jg,:)) .OR. var_in_output(jg)%vorw_ctmax) .AND. &
            l_celltracks_event_active .AND. .NOT. l_output(jg) ) THEN
        CALL rot_vertex (p_nh(jg)%prog(nnow(jg))%vn, p_patch(jg), p_int(jg), p_nh(jg)%diag%omega_z)
        ! Diagnose relative vorticity on cells
        CALL verts2cells_scalar(p_nh(jg)%diag%omega_z, p_patch(jg), &
           p_int(jg)%verts_aw_cells, p_nh(jg)%diag%vor)
      END IF

      DO k = 1,uh_max_nlayer

        ! update of UH_MAX (updraft helicity, max.  during the time interval "celltracks_interval") if required
        IF ( luh_max_out(jg,k) .AND. (l_output(jg) .OR. l_celltracks_event_active ) ) THEN
          CALL compute_field_uh_max( p_patch(jg), p_nh(jg)%metrics, p_nh(jg)%prog(nnow(jg)), p_nh(jg)%diag,  &
               &                     uh_max_zmin(k), uh_max_zmax(k), prm_diag(jg)%uh_max_3d(:,:,k) )
        END IF

      END DO

      ! update of VORW_CTMAX (Maximum rotation amplitude during the time interval "celltracks_interval") if required
      IF ( var_in_output(jg)%vorw_ctmax .AND. (l_output(jg) .OR. l_celltracks_event_active ) ) THEN
        CALL compute_field_vorw_ctmax( p_patch(jg), p_nh(jg)%metrics, p_nh(jg)%diag,  &
             &                         prm_diag(jg)%vorw_ctmax  )
      END IF

      ! update of W_CTMAX (Maximum updraft track during the ime interval "celltracks_interval") if required
      IF ( var_in_output(jg)%w_ctmax .AND. (l_output(jg) .OR. l_celltracks_event_active ) ) THEN
        CALL compute_field_w_ctmax( p_patch(jg), p_nh(jg)%metrics, p_nh(jg)%prog(nnow(jg)),  &
             &                      prm_diag(jg)%w_ctmax  )
      END IF


      ! Compute diagnostic 3D radar reflectivity (in linear units) for a specific domain
      ! if some derived output variables are present in any namelist
      ! for this domain, or if it is needed for statistical variables between output time steps on this domain.
      ! Has to be computed before pp_scheduler_process(simulation_status) and before statistical processing between timesteps below!

      IF (l_need_dbz3d) THEN
        CALL compute_field_dbz3D_lin( jg, p_patch(jg), p_nh(jg)%prog(nnow(jg)), p_nh(jg)%prog(nnow_rcf(jg)), &
             &                        p_nh(jg)%diag, prm_diag(jg), prm_diag(jg)%dbz3d_lin )
      END IF

      ! output of dbz_ctmax (column maximum reflectivity during a time interval (namelist param. celltracks_interval) is required
      IF ( var_in_output(jg)%dbzctmax .AND. (l_output(jg) .OR. l_dbz_event_active ) ) THEN
        CALL maximize_field_dbzctmax( p_patch(jg), jg, prm_diag(jg)%dbz3d_lin, prm_diag(jg)%dbz_ctmax )
      END IF

      ! output of echotop (minimum pressure where reflectivity exceeds threshold(s) 
      ! during a time interval (namelist param. echotop_meta(jg)%time_interval) is required:
      IF ( var_in_output(jg)%echotop .AND. (l_output(jg) .OR. l_dbz_event_active ) ) THEN
        CALL compute_field_echotop ( p_patch(jg), jg, p_nh(jg)%diag, prm_diag(jg)%dbz3d_lin, prm_diag(jg)%echotop )
      END IF

      ! output of echotopinm (maximum height where reflectivity exceeds threshold(s) 
      ! during a time interval (namelist param. echotop_meta(jg)%time_interval) is required:
      IF ( var_in_output(jg)%echotopinm .AND. (l_output(jg) .OR. l_dbz_event_active ) ) THEN
        CALL compute_field_echotopinm ( p_patch(jg), jg, p_nh(jg)%metrics, &
                                        prm_diag(jg)%dbz3d_lin, prm_diag(jg)%echotopinm )
      END IF
    END DO

    IF (ltimer) CALL timer_stop(timer_nh_diagnostics)

  END SUBROUTINE nwp_opt_diagnostics

  !>
  !! Subroutine collecting the calls for computing optional diagnostic output variables
  !! that should be called directly from the physics interface in order to consider
  !! the intermediate values of nesting time steps
  !!
  !! @par Revision History
  !! Developed by Daniel Rieger, DWD (2021-07-06)
  !!
  !!
  SUBROUTINE nwp_opt_diagnostics_2(p_patch, prm_diag, cosmu0, var_in_output, dt_phy)

    TYPE(t_patch),         INTENT(IN)    :: p_patch         !< current patch
    TYPE(t_nwp_phy_diag),  INTENT(INOUT) :: prm_diag        !< physics diagnostics
    REAL(wp),              INTENT(IN)    :: cosmu0(:,:)     !< Cosine of solar zenith angle
    TYPE(t_var_in_output), INTENT(IN)    :: var_in_output   !< variable contained in output nml
    REAL(wp),              INTENT(IN)    :: dt_phy          !< time interval for fast physics

    IF (ltimer) CALL timer_start(timer_nh_diagnostics)

    IF (var_in_output%dursun) THEN
      ! sunshine duration is an accumulative value like precipitation or runoff
      CALL compute_field_dursun( p_patch, dt_phy,                                           &
           &                     prm_diag%dursun, prm_diag%swflxsfc, prm_diag%swflx_up_sfc, &
           &                     prm_diag%swflx_dn_sfc_diff, cosmu0 )
    ENDIF

    IF (ltimer) CALL timer_stop(timer_nh_diagnostics)

  END SUBROUTINE nwp_opt_diagnostics_2

  !-------------------------------------------------------------------------
  !>
  !! Extended diagnostics for NWP physics interface - part 1
  !! Was included in mo_nh_interface_nwp before
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2013-01-07)
  !!
  SUBROUTINE nwp_diag_output_1(p_patch, p_diag, p_prog_rcf)

    TYPE(t_patch),   INTENT(in) :: p_patch     !< grid/patch info.
    TYPE(t_nh_diag), INTENT(in) :: p_diag      !< NH diagnostic state
    TYPE(t_nh_prog), INTENT(in) :: p_prog_rcf  !< state for tracer variables


    ! Local variables
    REAL(wp), DIMENSION(p_patch%nblks_c,p_patch%nlev) ::             &
      maxabs_u, maxabs_v, maxtemp, mintemp, maxqv, minqv, maxqc, minqc
    REAL(wp), DIMENSION(p_patch%nlev) ::               &
      umax, vmax, tmax, tmin, qvmax, qvmin, qcmax, qcmin

    ! loop indices
    INTEGER :: jc,jk,jb,jg

    INTEGER :: nlev                    !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    CALL message('mo_nwp_diagnosis:','Initial diagnostic output')

    maxabs_u(:,:) = 0._wp
    maxabs_v(:,:) = 0._wp
    maxtemp(:,:)  = 0._wp
    mintemp(:,:)  = 1.e20_wp
    maxqv(:,:)    = 0._wp
    minqv(:,:)    = 1.e20_wp
    maxqc(:,:)    = 0._wp
    minqc(:,:)    = 1.e20_wp

    nlev = p_patch%nlev
    jg   = p_patch%id

    ! Exclude the nest boundary zone
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_nchdom  = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          maxabs_u(jb,jk) = MAX(maxabs_u(jb,jk),ABS(p_diag%u(jc,jk,jb)))
          maxabs_v(jb,jk) = MAX(maxabs_v(jb,jk),ABS(p_diag%v(jc,jk,jb)))
          maxtemp(jb,jk)  = MAX(maxtemp(jb,jk),p_diag%temp(jc,jk,jb))
          mintemp(jb,jk)  = MIN(mintemp(jb,jk),p_diag%temp(jc,jk,jb))
          maxqv(jb,jk)    = MAX(maxqv(jb,jk),p_prog_rcf%tracer(jc,jk,jb,iqv))
          minqv(jb,jk)    = MIN(minqv(jb,jk),p_prog_rcf%tracer(jc,jk,jb,iqv))
          maxqc(jb,jk)    = MAX(maxqc(jb,jk),p_prog_rcf%tracer(jc,jk,jb,iqc))
          minqc(jb,jk)    = MIN(minqc(jb,jk),p_prog_rcf%tracer(jc,jk,jb,iqc))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DO jk = 1, nlev
      umax(jk)  = MAXVAL(maxabs_u(:,jk))
      vmax(jk)  = MAXVAL(maxabs_v(:,jk))
      tmax(jk)  = MAXVAL(maxtemp(:,jk))
      tmin(jk)  = MINVAL(mintemp(:,jk))
      qvmax(jk) = MAXVAL(maxqv(:,jk))
      qvmin(jk) = MINVAL(minqv(:,jk))
      qcmax(jk) = MAXVAL(maxqc(:,jk))
      qcmin(jk) = MINVAL(minqc(:,jk))
    ENDDO

    ! Finally take maximum/minimum over all PEs
    umax  = global_max(umax)
    vmax  = global_max(vmax)
    tmax  = global_max(tmax)
    tmin  = global_min(tmin)
    qvmax = global_max(qvmax)
    qvmin = global_min(qvmin)
    qcmax = global_max(qcmax)
    qcmin = global_min(qcmin)

    WRITE(message_text,'(a,i2)') 'max |U|, max |V|, min/max T, min/max QV,&
      & max QC per level in domain ',jg
    CALL message('', TRIM(message_text))
    DO jk = 1, nlev
      WRITE(message_text,'(a,i3,7(a,e12.5))') 'level ',jk,': u =',umax(jk),', v =',vmax(jk), &
        ', t =', tmin(jk),' ', tmax(jk),', qv =', qvmin(jk),' ', qvmax(jk), &
        ', qc =', qcmax(jk)   !,' ',qcmin(jk)
      CALL message('', TRIM(message_text))
    ENDDO

  END SUBROUTINE nwp_diag_output_1


  !-------------------------------------------------------------------------
  !>
  !! Extended diagnostics for NWP physics interface - part 2
  !! Was included in mo_nh_interface_nwp before
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2013-01-07)
  !!
  SUBROUTINE nwp_diag_output_2(p_patch, p_prog_rcf, prm_nwp_tend)

    TYPE(t_patch), TARGET,INTENT(in) :: p_patch      !< grid/patch info.
    TYPE(t_nh_prog),      INTENT(in) :: p_prog_rcf   !< state for TKE
    TYPE(t_nwp_phy_tend), INTENT(in) :: prm_nwp_tend !< physics tendencies

    ! Local variables

    ! variables for turbulence diagnostics
    REAL(wp) :: maxtke(p_patch%nblks_c,p_patch%nlevp1),tkemax(p_patch%nlevp1)
    REAL(wp), DIMENSION(p_patch%nblks_c,p_patch%nlev) :: maxtturb, maxuturb, maxvturb
    REAL(wp), DIMENSION(p_patch%nlev) :: tturbmax, uturbmax, vturbmax

    ! loop indices
    INTEGER :: jc,jk,jb,jg

    INTEGER :: nlev, nlevp1            !< number of model levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1
    jg     = p_patch%id

    ! Exclude the nest boundary zone
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_nchdom  = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    ! initialization
    maxtke(:,:)   = 0._wp
    maxtturb(:,:) = 0._wp
    maxuturb(:,:) = 0._wp
    maxvturb(:,:) = 0._wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlevp1
        DO jc = i_startidx, i_endidx
          maxtke(jb,jk) = MAX(maxtke(jb,jk),p_prog_rcf%tke(jc,jk,jb))
        ENDDO
      ENDDO

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          maxtturb(jb,jk) = MAX(maxtturb(jb,jk),ABS(prm_nwp_tend%ddt_temp_turb(jc,jk,jb)))
          maxuturb(jb,jk) = MAX(maxuturb(jb,jk),ABS(prm_nwp_tend%ddt_u_turb(jc,jk,jb)))
          maxvturb(jb,jk) = MAX(maxvturb(jb,jk),ABS(prm_nwp_tend%ddt_v_turb(jc,jk,jb)))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DO jk = 1, nlevp1
      tkemax(jk) = MAXVAL(maxtke(:,jk))
    ENDDO
    DO jk = 1, nlev
      tturbmax(jk) = MAXVAL(maxtturb(:,jk))
      uturbmax(jk) = MAXVAL(maxuturb(:,jk))
      vturbmax(jk) = MAXVAL(maxvturb(:,jk))
    ENDDO

    ! Take maximum over all PEs
    tkemax   = global_max(tkemax)
    tturbmax = global_max(tturbmax)
    uturbmax = global_max(uturbmax)
    vturbmax = global_max(vturbmax)

    WRITE(message_text,'(a,i2)') 'Extended turbulence diagnostic for domain ',jg
    CALL message('nwp_diag_output_2: ', TRIM(message_text))
    WRITE(message_text,'(a)') 'maximum TKE [m**2/s**2] and U,V,T-tendencies/s per level'
    CALL message('', TRIM(message_text))

    DO jk = 1, nlev
      WRITE(message_text,'(a,i3,4(a,e13.5))') 'level ',jk,': TKE =',tkemax(jk), &
        ', utend =',uturbmax(jk),', vtend =',vturbmax(jk),', ttend =',tturbmax(jk)
      CALL message('', TRIM(message_text))
    ENDDO
    jk = nlevp1
    WRITE(message_text,'(a,i3,a,e13.5)') 'level ',jk,': TKE =',tkemax(jk)
    CALL message('', TRIM(message_text))

  END SUBROUTINE nwp_diag_output_2

  !-------------------------------------------------------------------------
  !>
  !! Extended diagnostics for NWP physics interface
  !! for run-time min/max output of microphysics variables
  !!

  SUBROUTINE nwp_diag_output_minmax_micro(p_patch, p_prog, p_diag, ptr_tracer)

    TYPE(t_nh_prog), INTENT(in) :: p_prog      !< the dyn prog vars
    TYPE(t_patch),   INTENT(in) :: p_patch     !< grid/patch info.
    TYPE(t_nh_diag), INTENT(in) :: p_diag      !< NH diagnostic state
    !> tracer variables
    REAL(wp), CONTIGUOUS_ARGUMENT(in) :: ptr_tracer(:,:,:,:)


    ! Local variables
    REAL(wp), DIMENSION(p_patch%nblks_c) ::                                              &
         & qvmax, qcmax, qrmax, qimax, qsmax, qhmax, qgmax, tmax, wmax, qncmax, qnimax, qglmax, qhlmax, &
         & qvmin, qcmin, qrmin, qimin, qsmin, qhmin, qgmin, tmin, wmin, qncmin, qnimin, qglmin, qhlmin
    REAL(wp) ::                                                                          &
         & qvmaxi, qcmaxi, qrmaxi, qimaxi, qsmaxi, qhmaxi, qgmaxi, tmaxi, wmaxi, qncmaxi, qnimaxi, qglmaxi, qhlmaxi,  &
         & qvmini, qcmini, qrmini, qimini, qsmini, qhmini, qgmini, tmini, wmini, qncmini, qnimini, qglmini, qhlmini
    
    ! loop indices
    INTEGER :: jc,jk,jb,jg

    INTEGER :: nlev                    !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    CALL message('mo_nwp_diagnosis:','output min/max values of microphysics')

    nlev = p_patch%nlev
    jg   = p_patch%id

    ! Exclude the nest boundary zone 

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_nchdom  = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    ! Find local min/max
    wmax  = 0.0_wp
    wmin  = 0.0_wp
    tmax  = 0.0_wp
    tmin  = 0.0_wp
    qvmax = 0.0_wp
    qvmin = 0.0_wp
    qcmax = 0.0_wp
    qcmin = 0.0_wp
    qrmax = 0.0_wp
    qrmin = 0.0_wp
    qimax = 0.0_wp
    qimin = 0.0_wp
    qsmax = 0.0_wp
    qsmin = 0.0_wp
    qgmax = 0.0_wp
    qgmin = 0.0_wp
    qhmax = 0.0_wp
    qhmin = 0.0_wp
    
    qglmax = 0.0_wp
    qglmin = 0.0_wp
    qhlmax = 0.0_wp
    qhlmin = 0.0_wp

    qncmax = 0.0_wp
    qncmin = 0.0_wp
    qnimax = 0.0_wp
    qnimin = 0.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
         DO jc = i_startidx, i_endidx
            wmax(jb)  = MAX(wmax(jb), p_prog%w(jc,jk,jb))
            wmin(jb)  = MIN(wmin(jb), p_prog%w(jc,jk,jb))
            tmax(jb)  = MAX(tmax(jb), p_diag%temp(jc,jk,jb))
            tmin(jb)  = MIN(tmin(jb), p_diag%temp(jc,jk,jb))
            qvmax(jb) = MAX(qvmax(jb),ptr_tracer(jc,jk,jb,iqv))
            qvmin(jb) = MIN(qvmin(jb),ptr_tracer(jc,jk,jb,iqv))
            qcmax(jb) = MAX(qcmax(jb),ptr_tracer(jc,jk,jb,iqc))
            qcmin(jb) = MIN(qcmin(jb),ptr_tracer(jc,jk,jb,iqc))
            qrmax(jb) = MAX(qrmax(jb),ptr_tracer(jc,jk,jb,iqr))
            qrmin(jb) = MIN(qrmin(jb),ptr_tracer(jc,jk,jb,iqr))
            qimax(jb) = MAX(qimax(jb),ptr_tracer(jc,jk,jb,iqi))
            qimin(jb) = MIN(qimin(jb),ptr_tracer(jc,jk,jb,iqi))
            qsmax(jb) = MAX(qsmax(jb),ptr_tracer(jc,jk,jb,iqs))
            qsmin(jb) = MIN(qsmin(jb),ptr_tracer(jc,jk,jb,iqs))
            
            IF(atm_phy_nwp_config(jg)%inwp_gscp==4 &
                 & .OR.atm_phy_nwp_config(jg)%inwp_gscp==5 .OR. atm_phy_nwp_config(jg)%inwp_gscp==7)THEN
               qgmax(jb) = MAX(qgmax(jb),ptr_tracer(jc,jk,jb,iqg))
               qgmin(jb) = MIN(qgmin(jb),ptr_tracer(jc,jk,jb,iqg))
               qhmax(jb) = MAX(qhmax(jb),ptr_tracer(jc,jk,jb,iqh))
               qhmin(jb) = MIN(qhmin(jb),ptr_tracer(jc,jk,jb,iqh))
            END IF
            IF(atm_phy_nwp_config(jg)%inwp_gscp==7)THEN
               qglmax(jb) = MAX(qglmax(jb),ptr_tracer(jc,jk,jb,iqgl))
               qglmin(jb) = MIN(qglmin(jb),ptr_tracer(jc,jk,jb,iqgl))
               qhlmax(jb) = MAX(qhlmax(jb),ptr_tracer(jc,jk,jb,iqhl))
               qhlmin(jb) = MIN(qhlmin(jb),ptr_tracer(jc,jk,jb,iqhl))
            END IF
            IF(atm_phy_nwp_config(jg)%inwp_gscp==5)THEN
               qncmax(jb) = MAX(qncmax(jb),ptr_tracer(jc,jk,jb,iqnc))
               qncmin(jb) = MIN(qncmin(jb),ptr_tracer(jc,jk,jb,iqnc))
               qnimax(jb) = MAX(qnimax(jb),ptr_tracer(jc,jk,jb,iqni))
               qnimin(jb) = MIN(qnimin(jb),ptr_tracer(jc,jk,jb,iqni))
            END IF
         ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Take maximum/minimum over blocks
    wmaxi = MAXVAL(wmax(i_startblk:i_endblk))
    wmini = MINVAL(wmin(i_startblk:i_endblk))
    tmaxi  = MAXVAL(tmax(i_startblk:i_endblk))
    tmini  = MINVAL(tmin(i_startblk:i_endblk))
    qvmaxi = MAXVAL(qvmax(i_startblk:i_endblk))
    qvmini = MINVAL(qvmin(i_startblk:i_endblk))
    qcmaxi = MAXVAL(qcmax(i_startblk:i_endblk))
    qcmini = MINVAL(qcmin(i_startblk:i_endblk))
    qrmaxi = MAXVAL(qrmax(i_startblk:i_endblk))
    qrmini = MINVAL(qrmin(i_startblk:i_endblk))
    qimaxi = MAXVAL(qimax(i_startblk:i_endblk))
    qimini = MINVAL(qimin(i_startblk:i_endblk))
    qsmaxi = MAXVAL(qsmax(i_startblk:i_endblk))
    qsmini = MINVAL(qsmin(i_startblk:i_endblk))
    IF(atm_phy_nwp_config(jg)%inwp_gscp==4 &
         & .OR.atm_phy_nwp_config(jg)%inwp_gscp==5 .OR.atm_phy_nwp_config(jg)%inwp_gscp==7)THEN
       qgmaxi = MAXVAL(qgmax(i_startblk:i_endblk))
       qgmini = MINVAL(qgmin(i_startblk:i_endblk))
       qhmaxi = MAXVAL(qhmax(i_startblk:i_endblk))
       qhmini = MINVAL(qhmin(i_startblk:i_endblk))
    END IF
    IF(atm_phy_nwp_config(jg)%inwp_gscp==7)THEN
       qglmaxi = MAXVAL(qglmax(i_startblk:i_endblk))
       qglmini = MINVAL(qglmin(i_startblk:i_endblk))
       qhlmaxi = MAXVAL(qhlmax(i_startblk:i_endblk))
       qhlmini = MINVAL(qhlmin(i_startblk:i_endblk))
    END IF
    IF(atm_phy_nwp_config(jg)%inwp_gscp==5)THEN
       qncmaxi = MAXVAL(qncmax(i_startblk:i_endblk))
       qncmini = MINVAL(qncmin(i_startblk:i_endblk))
       qnimaxi = MAXVAL(qnimax(i_startblk:i_endblk))
       qnimini = MINVAL(qnimin(i_startblk:i_endblk))
    END IF

    ! Take maximum/minimum over all PEs
    wmaxi  = global_max(wmaxi)
    wmini  = global_min(wmini)
    tmaxi  = global_max(tmaxi)
    tmini  = global_min(tmini)
    qvmaxi = global_max(qvmaxi)
    qvmini = global_min(qvmini)
    qcmaxi = global_max(qcmaxi)
    qcmini = global_min(qcmini)
    qrmaxi = global_max(qrmaxi)
    qrmini = global_min(qrmini)
    qimaxi = global_max(qimaxi)
    qimini = global_min(qimini)
    qsmaxi = global_max(qsmaxi)
    qsmini = global_min(qsmini)
    IF(atm_phy_nwp_config(jg)%inwp_gscp==4 &
         & .OR.atm_phy_nwp_config(jg)%inwp_gscp==5 .OR. atm_phy_nwp_config(jg)%inwp_gscp==7)THEN
       qgmaxi = global_max(qgmaxi)
       qgmini = global_min(qgmini)
       qhmaxi = global_max(qhmaxi)
       qhmini = global_min(qhmini)
    END IF
    IF(atm_phy_nwp_config(jg)%inwp_gscp==7)THEN
       qglmaxi = global_max(qglmaxi)
       qglmini = global_min(qglmini)
       qhlmaxi = global_max(qhlmaxi)
       qhlmini = global_min(qhlmini)
    END IF
    IF(atm_phy_nwp_config(jg)%inwp_gscp==5)THEN
       qncmaxi = global_max(qncmaxi)
       qncmini = global_min(qncmini)
       qnimaxi = global_max(qnimaxi)
       qnimini = global_min(qnimini)
    END IF

    ! Standard output
    SELECT CASE (atm_phy_nwp_config(jg)%inwp_gscp)
    CASE(1)
       WRITE(message_text,'(A10,8A11)')   '  var: ', 'w','qv','qc','qr','qi','qs'
       CALL message("",TRIM(message_text))
       WRITE(message_text,'(A10,8E11.3)') '  max: ', wmaxi,qvmaxi,qcmaxi,qrmaxi,qimaxi,qsmaxi
       CALL message("",TRIM(message_text))
       WRITE(message_text,'(A10,8E11.3)') '  min: ', wmini,qvmini,qcmini,qrmini,qimini,qsmini
       CALL message("",TRIM(message_text))
    CASE(4)
       WRITE(message_text,'(A10,9A11)')   '  var: ', 'w','qv','qc','qr','qi','qs','qg','qh','temp'
       CALL message("",TRIM(message_text))
       WRITE(message_text,'(A10,9E11.3)') '  max: ', wmaxi,qvmaxi,qcmaxi,qrmaxi,qimaxi,qsmaxi,qgmaxi,qhmaxi,tmaxi
       CALL message("",TRIM(message_text))
       WRITE(message_text,'(A10,9E11.3)') '  min: ', wmini,qvmini,qcmini,qrmini,qimini,qsmini,qgmini,qhmini,tmini
       CALL message("",TRIM(message_text))
    CASE(5)
       WRITE(message_text,'(A10,10A11)')   '  var: ', 'w','qv','qc','qr','qi','qs','qg','qh','qnc','qni'
       CALL message("",TRIM(message_text))
       WRITE(message_text,'(A10,10E11.3)') '  max: ', wmaxi,qvmaxi,qcmaxi,qrmaxi,qimaxi,qsmaxi,qgmaxi,qhmaxi,qncmaxi,qnimaxi
       CALL message("",TRIM(message_text))
       WRITE(message_text,'(A10,10E11.3)') '  min: ', wmini,qvmini,qcmini,qrmini,qimini,qsmini,qgmini,qhmini,qncmini,qnimini
       CALL message("",TRIM(message_text))       
    CASE(7)
       WRITE(message_text,'(A10,10A11)')   '  var: ', 'w','qv','qc','qr','qi','qs','qg','qh','qgl','qhl'
       CALL message("",TRIM(message_text))
       WRITE(message_text,'(A10,10E11.3)') '  max: ', wmaxi,qvmaxi,qcmaxi,qrmaxi,qimaxi,qsmaxi,qgmaxi,qhmaxi,qglmaxi,qhlmaxi
       CALL message("",TRIM(message_text))
       WRITE(message_text,'(A10,10E11.3)') '  min: ', wmini,qvmini,qcmini,qrmini,qimini,qsmini,qgmini,qhmini,qhlmini,qhlmini
       CALL message("",TRIM(message_text))       
    CASE DEFAULT       
          CALL finish('nwp_diag_output_minmax_micro', 'Cloud microphysics scheme not yet known in diagnostics.')
    END SELECT


  END SUBROUTINE nwp_diag_output_minmax_micro


END MODULE mo_nwp_diagnosis

