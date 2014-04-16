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
!!
!! $Id: n/a$
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_diagnosis


  USE mo_kind,               ONLY: wp

  USE mo_impl_constants,     ONLY: itccov, itconv, itradheat, itturb, itfastphy, &
    &                              min_rlcell_int
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_exception,          ONLY: message, message_text
  USE mo_model_domain,       ONLY: t_patch
  USE mo_run_config,         ONLY: msg_level, iqv, iqc, iqi
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_parallel_config,    ONLY: nproma
  USE mo_lnd_nwp_config,     ONLY: nlev_soil, ntiles_total
  USE mo_nwp_lnd_types,      ONLY: t_lnd_diag, t_wtr_prog
  USE mo_physical_constants, ONLY: tmelt, grav, cpd, vtmpc1
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_io_config,          ONLY: lflux_avg
  USE mo_sync,               ONLY: global_max, global_min
  USE mo_satad,              ONLY: sat_pres_water, spec_humi

  IMPLICIT NONE

  PRIVATE

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  PUBLIC  :: nwp_statistics
  PUBLIC  :: nwp_diag_for_output
  PUBLIC  :: nwp_diag_output_1
  PUBLIC  :: nwp_diag_output_2

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
                            & kstart_moist,               & !in
                            & ih_clch, ih_clcm,           & !in
                            & pt_patch, p_metrics,        & !in
                            & pt_prog, pt_prog_rcf,       & !in
                            & pt_diag,                    & !inout
                            & prm_diag                    ) !inout   
                            

    LOGICAL,            INTENT(IN)   :: lcall_phy_jg(:) !< physics package time control (switches)
                                                        !< for domain jg
    REAL(wp),           INTENT(IN)   :: dt_phy_jg(:)    !< time interval for all physics
                                                        !< packages on domain jg
    REAL(wp),           INTENT(IN)   :: p_sim_time

    TYPE(t_patch),      INTENT(IN)   :: pt_patch    !<grid/patch info.
    TYPE(t_nh_diag),    INTENT(INOUT):: pt_diag     !<the diagnostic variables
    TYPE(t_nh_prog),    INTENT(IN)   :: pt_prog     !<the prognostic variables
    TYPE(t_nh_prog),    INTENT(IN)   :: pt_prog_rcf !<the prognostic variables (with
                                                         !< red. calling frequency for tracers!
    TYPE(t_nh_metrics), INTENT(in)   :: p_metrics

    TYPE(t_nwp_phy_diag), INTENT(inout):: prm_diag

    ! Local
    INTEGER :: nlev, nlevp1            !< number of full levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    REAL(wp):: z_help, r_sim_time
    REAL(wp):: t_wgt                   !< weight for running time average

    INTEGER :: jc,jk,jb,jg      !block index
    INTEGER :: jk_max
    INTEGER :: kstart_moist
    INTEGER :: ih_clch, ih_clcm


    REAL(wp):: clearsky(nproma)
    REAL(wp):: d_theta_dz, d_theta_dz_max


    REAL(wp), PARAMETER:: eps_clc = 1.e-7_wp

  !-----------------------------------------------------------------


    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1    

    ! Inverse of simulation time
    r_sim_time = 1._wp/MAX(1.e-6_wp, p_sim_time)

    ! time average weight
    t_wgt = dt_phy_jg(itfastphy)/MAX(1.e-6_wp, p_sim_time)


    ! exclude nest boundary interpolation zone
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
    

!$OMP PARALLEL
    IF ( lcall_phy_jg(itccov) ) THEN

!$OMP DO PRIVATE(jc,jk,jb,z_help,i_startidx,i_endidx,clearsky) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)


        ! if cloud cover is called, vertical integration of cloud content
        ! (for iqv, iqc, iqi)
        !
        prm_diag%tot_cld_vi(i_startidx:i_endidx,jb,1:3) = 0.0_wp

        DO jk = kstart_moist, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

           z_help = p_metrics%ddqz_z_full(jc,jk,jb) * pt_prog%rho(jc,jk,jb)  

           ! TQV, TQC, TQI
           prm_diag%tot_cld_vi(jc, jb,iqv) = prm_diag%tot_cld_vi(jc, jb,iqv)    + &
                                             z_help * prm_diag%tot_cld(jc,jk,jb,iqv)
           prm_diag%tot_cld_vi(jc, jb,iqc) = prm_diag%tot_cld_vi(jc, jb,iqc)    + &
                                             z_help * prm_diag%tot_cld(jc,jk,jb,iqc)
           prm_diag%tot_cld_vi(jc, jb,iqi) = prm_diag%tot_cld_vi(jc, jb,iqi)    + &
                                             z_help * prm_diag%tot_cld(jc,jk,jb,iqi)
          ENDDO  ! jc
        ENDDO  ! jk



        ! cloud cover calculation
        ! note: the conversion into % is done within the internal output postprocessing
        !
        DO jc = i_startidx, i_endidx
          clearsky(jc) = 1._wp - prm_diag%clc(jc,kstart_moist,jb)
        ENDDO

        DO jk = kstart_moist+1, ih_clch
          DO jc = i_startidx, i_endidx
            clearsky(jc) = clearsky(jc)*    &
            &  ( 1._wp - MAX( prm_diag%clc(jc,jk  ,jb), prm_diag%clc(jc,jk-1,jb))) &
            & /( 1._wp - MIN( prm_diag%clc(jc,jk-1,jb), 1._wp - eps_clc) )
          ENDDO
        ENDDO

        ! store high-level clouds
        DO jc = i_startidx, i_endidx
          prm_diag%clch(jc,jb) = MAX( 0._wp, 1._wp - clearsky(jc) - eps_clc)
        ENDDO

        ! continue downward for total cloud cover
        DO jk = ih_clch+1, nlev
          DO jc = i_startidx, i_endidx
            clearsky(jc) = clearsky(jc)*    &
            &  ( 1._wp - MAX( prm_diag%clc(jc,jk,jb), prm_diag%clc(jc,jk-1,jb))) &
            & /( 1._wp - MIN( prm_diag%clc(jc,jk-1,jb), 1._wp - eps_clc) )
          ENDDO
        ENDDO

        ! store total cloud cover, start for mid-level clouds
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          prm_diag%clct(jc,jb) = MAX( 0._wp, 1._wp - clearsky(jc) - eps_clc)
          clearsky(jc) = 1._wp - prm_diag%clc(jc,ih_clch+1,jb)
        ENDDO

        ! mid-level clouds
        DO jk = ih_clch+2, ih_clcm
          DO jc = i_startidx, i_endidx
            clearsky(jc) = clearsky(jc)*    &
            &  ( 1._wp - MAX( prm_diag%clc(jc,jk,jb), prm_diag%clc(jc,jk-1,jb))) &
            & /( 1._wp - MIN( prm_diag%clc(jc,jk-1,jb), 1._wp - eps_clc) )
          ENDDO
        ENDDO

        ! store mid-level cloud cover, start for low-level clouds
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          prm_diag%clcm(jc,jb) = MAX( 0._wp, 1._wp - clearsky(jc) - eps_clc)

          clearsky(jc) = 1._wp - prm_diag%clc(jc,ih_clcm+1,jb)
        ENDDO

        ! continue downward for mid-level clouds
        DO jk = ih_clcm+2, nlev
          DO jc = i_startidx, i_endidx
            clearsky(jc) = clearsky(jc)*    &
            &  ( 1._wp - MAX( prm_diag%clc(jc,jk,jb), prm_diag%clc(jc,jk-1,jb))) &
            & /( 1._wp - MIN( prm_diag%clc(jc,jk-1,jb), 1._wp - eps_clc) )
          ENDDO
        ENDDO

        ! store low-level clouds
        DO jc = i_startidx, i_endidx
          prm_diag%clcl(jc,jb) = MAX( 0._wp, 1._wp - clearsky(jc) - eps_clc)
        ENDDO

      ENDDO ! nblks
!$OMP END DO

    END IF !cloud cover



 !! Calculate vertically integrated values of the grid-scale tracers q1, q2, q3, q4 and q5
 
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,z_help) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      pt_diag%tracer_vi(i_startidx:i_endidx,jb,1:5) = 0.0_wp
      DO jk = kstart_moist, nlev
        DO jc = i_startidx, i_endidx 

          z_help = p_metrics%ddqz_z_full(jc,jk,jb) * pt_prog%rho(jc,jk,jb)  
          pt_diag%tracer_vi(jc, jb,1:5) = pt_diag%tracer_vi(jc, jb,1:5)    + &
                                 z_help * pt_prog_rcf%tracer(jc,jk,jb,1:5) 

        ENDDO
      ENDDO

    ENDDO ! nblks   
!$OMP END DO



    ! Calculation of average/accumulated values since model start
    !
    ! Compute
    !
    ! wind
    !-----------
    ! - maximum gust (including convective contribution)
    !
    ! cloud/rain
    !-----------
    ! - time averaged precipitation rates (total, grid-scale, convective)
    ! - time averaged total cloud cover
    ! - time averaged TQV, TQC, TQI, TQR, TQS
    ! - time averaged TQV_DIA, TQC_DIA, TQI_DIA
    !
    ! turbulent fluxes
    !-----------------
    ! - surface latent heat flux
    ! - surface latent heat flux from bare soil 
    ! - surface sensible heat flux
    ! - surface moisture flux
    ! - surface u/v-momentum flux
    !
    ! radiative fluxes
    !------------------
    ! - top net solar radiation
    ! - top down solar radiation
    ! - top net thermal radiation
    ! - surface net solar radiation
    ! - surface net thermal radiation

    IF ( p_sim_time > 1.e-6_wp) THEN

!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

!DIR$ IVDEP
        DO jc = i_startidx, i_endidx

          ! maximum 10m gust, including convective contribution
          ! (reset is done on a regular basis in reset_action)
          prm_diag%gust10(jc,jb) = MAX(prm_diag%gust10(jc,jb),                       &
            &                    prm_diag%dyn_gust(jc,jb) + prm_diag%con_gust(jc,jb) )

          ! time averaged total precipitation rate
          prm_diag%tot_prec_rate_avg(jc,jb) = prm_diag%tot_prec(jc,jb) &
            &                               * r_sim_time

          ! time averaged grid scale precipitation rate
          prm_diag%gsp_prec_rate_avg(jc,jb) = (prm_diag%rain_gsp(jc,jb) &
            &                               + prm_diag%snow_gsp(jc,jb)) &
            &                               * r_sim_time

          ! time averaged convective precipitation rate
          prm_diag%con_prec_rate_avg(jc,jb) = (prm_diag%rain_con(jc,jb) & 
            &                               + prm_diag%snow_con(jc,jb)) &
            &                               * r_sim_time

          ! time averaged total cloud cover
          prm_diag%clct_avg(jc,jb) = time_avg(prm_diag%clct_avg(jc,jb), &
            &                                 prm_diag%clct    (jc,jb), &
            &                                 t_wgt)

          ! time averaged TQV, TQC, TQI, TQR, TQS  
          pt_diag%tracer_vi_avg(jc,jb,1:5) = (1._wp - t_wgt)*pt_diag%tracer_vi_avg(jc,jb,1:5) &
            &                              + t_wgt * pt_diag%tracer_vi(jc,jb,1:5)

          ! time averaged TQV_DIA, TQC_DIA, TQI_DIA
          prm_diag%tot_cld_vi_avg(jc,jb,1:3) = (1._wp - t_wgt)*prm_diag%tot_cld_vi_avg(jc,jb,1:3) &
            &                                + t_wgt * prm_diag%tot_cld_vi(jc,jb,1:3)
        ENDDO

        !Add more precip vars in case of two moment microphysics
        IF(atm_phy_nwp_config(jg)%inwp_gscp==4)THEN
         
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            prm_diag%gsp_prec_rate_avg(jc,jb) = prm_diag%gsp_prec_rate_avg(jc,jb) + &
              &                               ( prm_diag%ice_gsp(jc,jb)  +    &
              &                                 prm_diag%hail_gsp(jc,jb) +    &
              &                                 prm_diag%graupel_gsp(jc,jb) ) &
              &                               * r_sim_time
          END DO

        END IF!IF inwp_gscp==4
        

        IF (lflux_avg) THEN

          IF (lcall_phy_jg(itturb)) THEN
!DIR$ IVDEP
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

              ! time averaged surface u-momentum flux
              prm_diag%aumfl_s(jc,jb) = time_avg(prm_diag%aumfl_s(jc,jb), &
                &                                prm_diag%umfl_s (jc,jb), &
                &                                t_wgt )

              ! time averaged surface v-momentum flux
              prm_diag%avmfl_s(jc,jb) = time_avg(prm_diag%avmfl_s(jc,jb), &
                &                                prm_diag%vmfl_s (jc,jb), &
                &                                t_wgt )

            ENDDO  ! jc
            DO jk = 1, nlev_soil
!DIR$ IVDEP
              DO jc = i_startidx, i_endidx
              prm_diag%alhfl_pl(jc,jk,jb) = time_avg(prm_diag%alhfl_pl(jc,jk,jb), &
                &                                    prm_diag%lhfl_pl (jc,jk,jb), &
                &                                    t_wgt)
              ENDDO  ! jc
            ENDDO  ! jk

          ENDIF  ! inwp_turb > 0


          IF ( lcall_phy_jg(itradheat) ) THEN
            !sum up for averaged fluxes
            !T.R.: this is not correct for output after 1st timestep,
            !e.g. dt_phy_jg(itradheat) may then be greater than p_sim_time
            !leading to wrong averaging.
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              ! time averaged shortwave net flux at surface
              prm_diag%swflxsfc_a(jc,jb) = time_avg(prm_diag%swflxsfc_a(jc,jb), &
                &                                   prm_diag%swflxsfc  (jc,jb), &
                &                                   t_wgt)

              ! time averaged longwave net flux at surface
              prm_diag%lwflxsfc_a(jc,jb) = time_avg(prm_diag%lwflxsfc_a(jc,jb), &
                &                                   prm_diag%lwflxsfc  (jc,jb), &
                &                                   t_wgt)

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

            ENDDO

          ENDIF  ! lcall_phy_jg(itradheat)

        ELSEIF (.NOT. lflux_avg) THEN


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

              ! accumulated surface u-momentum flux
              prm_diag%aumfl_s(jc,jb) = prm_diag%aumfl_s(jc,jb)        &
                                &   + prm_diag%umfl_s(jc,jb)           &
                                &   * dt_phy_jg(itfastphy)

              ! accumulated surface v-momentum flux
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
          ENDIF  ! inwp_turb > 0


          IF ( lcall_phy_jg(itradheat) ) THEN
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              ! accumulated shortwave net flux at surface
              prm_diag%swflxsfc_a(jc,jb) = prm_diag%swflxsfc_a(jc,jb) &
                                   &   + prm_diag%swflxsfc(jc,jb)     &
                                   &   * dt_phy_jg(itfastphy)

              ! accumulated longwave net flux at surface
              prm_diag%lwflxsfc_a(jc,jb) = prm_diag%lwflxsfc_a(jc,jb) &
                                   &   + prm_diag%lwflxsfc(jc,jb)     &
                                   &   * dt_phy_jg(itfastphy)

              ! accumulated shortwave net flux at TOA
              prm_diag%swflxtoa_a(jc,jb) = prm_diag%swflxtoa_a(jc,jb) &
                                   &   + prm_diag%swflxtoa(jc,jb)     &
                                   &   * dt_phy_jg(itfastphy)

              ! accumulated longwave net flux at TOA
              prm_diag%lwflxtoa_a(jc,jb) = prm_diag%lwflxtoa_a(jc,jb) &
                                   &   + prm_diag%lwflxall(jc,1,jb)   &
                                   &   * dt_phy_jg(itfastphy)

              ! accumulated top down solar radiation
              prm_diag%asod_t    (jc,jb) = prm_diag%asod_t(jc,jb)     &
                                   &   + prm_diag%flxdwswtoa(jc,jb)   &
                                   &   * dt_phy_jg(itfastphy)
            END DO
          ENDIF  ! lcall_phy_jg(itradheat)

        ENDIF  ! lflux_avg

      ENDDO ! nblks     
!$OMP END DO NOWAIT

    END IF  ! p_sim_time





!  Calculation of average/accumulated u and v has been moved here from nwp_interface. 
!  Thus there is no need anymore for a LES-specific code. (Daniel Reinert, Jan 2014)
!  -included calculation of boundary layer height (Anurag Dipankar, MPI Octo 2013).
!  -soon all these LES diagnostics have to be moved to different module.

    IF( atm_phy_nwp_config(jg)%is_les_phy )THEN

!     !2D Boundary layer height
!     Switch to a more general formulation applicable for stable case as well (AD,2013-11-16)
       

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,d_theta_dz,jk_max,d_theta_dz_max) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
              & i_startidx, i_endidx, rl_start, rl_end)
            DO jc = i_startidx, i_endidx           
              d_theta_dz_max = -1.0_wp
              jk_max = -1
              DO jk = 5 , nlev-1
                d_theta_dz  = (pt_diag%temp(jc,jk-1,jb)/pt_prog%exner(jc,jk-1,jb) - &
                      pt_diag%temp(jc,jk,jb)/pt_prog%exner(jc,jk,jb))*p_metrics%inv_ddqz_z_half(jc,jk,jb) 

                IF(d_theta_dz > d_theta_dz_max)THEN
                  jk_max = jk
                  d_theta_dz_max = d_theta_dz
                END IF

              END DO !jk

              prm_diag%z_pbl(jc,jb) = p_metrics%z_mc(jc,jk_max,jb)

            ENDDO
          ENDDO ! jb
!$OMP END DO

    END IF !is_les_phy

!$OMP END PARALLEL  


  END SUBROUTINE nwp_statistics





  !>
  !! Diagnostics which are only required for output
  !!
  !! Diagnostics which are only required for output. Should contain all 
  !! those computations which are purely diagnostic and only required for 
  !! (meteogram) output.
  !!
  !! Available diagnostics:
  !! - height of convection base and top: hbas_con, htop_con
  !! - height of the top of dry convection: htop_dc
  !! - height of 0 deg C level: hzerocl
  !! - t_ice is filled with t_so(0) for non-ice points (h_ice=0)
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
  SUBROUTINE nwp_diag_for_output(kstart_moist,            & !in
                            & ih_clch, ih_clcm,           & !in
                            & pt_patch, p_metrics,        & !in
                            & pt_prog_rcf,                & !in
                            & pt_diag,                    & !in
                            & lnd_diag,                   & !in
                            & p_prog_wtr_now,             & !in
                            & prm_diag                    ) !inout    
              
    INTEGER,         INTENT(IN)   :: kstart_moist
    INTEGER,         INTENT(IN)   :: ih_clch, ih_clcm

    TYPE(t_patch),   INTENT(IN)   :: pt_patch    !<grid/patch info.
    TYPE(t_nh_prog), INTENT(IN)   :: pt_prog_rcf !<the prognostic variables (with
                                                 !< red. calling frequency for tracers!
    TYPE(t_nh_metrics)  ,INTENT(IN) :: p_metrics

    TYPE(t_nh_diag),     INTENT(IN)   :: pt_diag     ! the diagnostic variables
    TYPE(t_lnd_diag),    INTENT(IN)   :: lnd_diag    ! land diag state
    TYPE(t_wtr_prog),    INTENT(INOUT):: p_prog_wtr_now ! water prognostic state (now)
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

  !-----------------------------------------------------------------


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

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,mlab,ztp,zqp,zbuoy,zqsat,zcond) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      !
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)


      IF (atm_phy_nwp_config(jg)%lproc_on(itconv)) THEN  ! convection parameterization switched on

        !
        ! height of convection base and top, hbas_con, htop_con
        ! 
        DO jc = i_startidx, i_endidx
          IF ( prm_diag%locum(jc,jb)) THEN
            prm_diag%hbas_con(jc,jb) = p_metrics%z_ifc( jc, prm_diag%mbas_con(jc,jb), jb)
            prm_diag%htop_con(jc,jb) = p_metrics%z_mc ( jc, prm_diag%mtop_con(jc,jb), jb)
          ELSE
            prm_diag%hbas_con(jc,jb) = 0._wp
            prm_diag%htop_con(jc,jb) = 0._wp
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
            IF ( prm_diag%locum(jc,jb)) THEN
              prm_diag%htop_dc(jc,jb) = MIN( prm_diag%htop_dc(jc,jb),      &
             &                               prm_diag%hbas_con(jc,jb) )
            END IF
          ELSE
            prm_diag%htop_dc(jc,jb) = MIN( 0._wp, p_metrics%z_ifc(jc,nlevp1,jb) )
          END IF
        ENDDO

      END IF !convection parameterization switched on


      !
      ! height of 0 deg C level "hzerocl". Not higher than htop_moist_proc
      !
      ! Surface temperature below 0 deg C
      WHERE( pt_diag%temp(i_startidx:i_endidx,nlev,jb) < tmelt)
        prm_diag%hzerocl(i_startidx:i_endidx,jb) = zundef
      ELSEWHERE
        prm_diag%hzerocl(i_startidx:i_endidx,jb) = 0._wp
      END WHERE

      !AD(MPIM): ending the loop at kstart_moist+1 to avoid runtime error in dry case
      DO jk = nlev, kstart_moist+1, -1
        DO jc = i_startidx, i_endidx 
          IF ( prm_diag%hzerocl(jc,jb) /= 0._wp) THEN
            CYCLE
          ELSE IF ( pt_diag%temp(jc,jk  ,jb) >= tmelt .AND. &
           &        pt_diag%temp(jc,jk-1,jb) <  tmelt ) THEN
            prm_diag%hzerocl(jc,jb) = p_metrics%z_ifc(jc,jk-1,jb) -  &
           &      ( p_metrics%z_ifc(jc,jk-1,jb) - p_metrics%z_ifc(jc,jk,jb) )*  &
           &      (    pt_diag%temp(jc,jk-1,jb) - tmelt ) /                     &
           &      (    pt_diag%temp(jc,jk-1,jb) - pt_diag%temp(jc,jk,jb) )
          END IF
        ENDDO
      ENDDO


      ! Fill t_ice with t_so(1) for ice-free points (h_ice<=0)
      ! This was demanded by FE14 (surface analysis)
      !
      ! Note, that t_ice contains ice temperature information from 
      ! the sea ice model as well as the lake model.
      !
      ! Furthermore, note that filling t_ice with t_so(1) only makes 
      ! sense when running without tiles. When using tiles, contains 
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

    ENDDO  ! jb
!$OMP END DO

!$OMP END PARALLEL  


  END SUBROUTINE nwp_diag_for_output


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
  SUBROUTINE nwp_diag_output_2(p_patch, p_prog_rcf, prm_nwp_tend, lcall_turb)

    TYPE(t_patch), TARGET,INTENT(in) :: p_patch      !< grid/patch info.
    TYPE(t_nh_prog),      INTENT(in) :: p_prog_rcf   !< state for TKE
    TYPE(t_nwp_phy_tend), INTENT(in) :: prm_nwp_tend !< physics tendencies

    LOGICAL,  INTENT(in) :: lcall_turb ! switch if turbulence has been called

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


    ! In case that turbulence diagnostics are computed
    IF (msg_level >= 18) THEN
      maxtke(:,:)   = 0._wp
      maxtturb(:,:) = 0._wp
      maxuturb(:,:) = 0._wp
      maxvturb(:,:) = 0._wp
    ENDIF

!$OMP PARALLEL

    ! Extended turbulence diagnostics if msg_level >= 18
    IF (lcall_turb .AND. msg_level >= 18) THEN

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

    ENDIF
!$OMP END PARALLEL

    IF (msg_level >= 18 .AND. lcall_turb) THEN ! extended turbulence diagnostic
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
    ENDIF

  END SUBROUTINE nwp_diag_output_2
!
!===========================================================================
!
  !>
  !! Computes updated time average
  !!
  !! Computes updated time average for a particular field
  !!
  !! @Literature
  !! Based on proposal found in 
  !! Jochen Froehlich, 2006:Large Eddy Simulation turbulenter Str�mungen, Teubner
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-17)
  !!
  FUNCTION time_avg (psi_avg_old, psi_inst, wgt)  RESULT (psi_avg_new)

    REAL(wp), INTENT(IN) :: psi_avg_old       !< time average at t(n-1)
    REAL(wp), INTENT(IN) :: psi_inst          !< instantaneous value
    REAL(wp), INTENT(IN) :: wgt               !< weight (=dt/sim_time)

    ! Result
    REAL(wp) :: psi_avg_new                   !< updated time average

    !--------------------------------------------------------------------

    ! compute updated time average
    psi_avg_new = (1._wp - wgt)*psi_avg_old + wgt*psi_inst

  END FUNCTION time_avg

END MODULE mo_nwp_diagnosis

