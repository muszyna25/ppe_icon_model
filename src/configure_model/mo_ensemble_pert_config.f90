!>
!! @brief Ensemble perturbations of nwp physics
!!
!! configuration setup for ensemble physics perturbations
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! Initial revision by Guenther Zaengl, DWD (2015-04-23)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ensemble_pert_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom
  USE mo_nwp_tuning_config,  ONLY: tune_gkwake, tune_gkdrag, tune_gfluxlaun, tune_zvz0i,  &
    &                        tune_entrorg, tune_capdcfac_et, tune_box_liq
  USE mo_turbdiff_config,    ONLY: turbdiff_config
  USE mo_gribout_config,     ONLY: gribout_config
  USE mo_exception,          ONLY: message_text, message

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: use_ensemble_pert, configure_ensemble_pert
  PUBLIC :: range_gkwake, range_gkdrag, range_gfluxlaun, range_zvz0i, range_entrorg, range_capdcfac_et, &
            range_box_liq, range_tkhmin, range_tkmmin, range_rlam_heat

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for ensemble perturbations
  !!--------------------------------------------------------------------------


    ! namelist variables
  REAL(wp) :: &                    !< low level wake drag constant
    &  range_gkwake

  REAL(wp) :: &                    !< gravity wave drag constant
    &  range_gkdrag

  REAL(wp) :: &                    !< gravity wave flux emission
    &  range_gfluxlaun

  REAL(wp) :: &                    !< Terminal fall velocity of ice 
    &  range_zvz0i

  REAL(wp) :: &                    !< Entrainment parameter for deep convection valid at dx=20 km 
    &  range_entrorg

  REAL(wp) :: &                    !< Fraction of CAPE diurnal cycle correction applied in the extratropics
    &  range_capdcfac_et            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< Box width for liquid clouds assumed in the cloud cover scheme
    &  range_box_liq                ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Minimum vertical diffusion for heat/moisture 
    &  range_tkhmin

  REAL(wp) :: &                    !< Minimum vertical diffusion for momentum 
    &  range_tkmmin

  REAL(wp) :: &                    !< Laminar transport resistance parameter 
    &  range_rlam_heat

  LOGICAL :: use_ensemble_pert     !< main switch

  CONTAINS


  !>
  !! Application of the ensemble perturbation to the config/namelist variables 
  !!
  !! This is done based on randum numbers determined by the ensemble member ID
  !!
  !! @par Revision History
  !! Initial revision by Guenther Zaengl, DWD (2015-04-23)
  !!
  SUBROUTINE configure_ensemble_pert

    INTEGER, ALLOCATABLE :: rnd_seed(:)
    INTEGER  :: rnd_size, i
    REAL(wp) :: rnd_num, rnd_fac


    IF (use_ensemble_pert) THEN

      CALL RANDOM_SEED(rnd_size)
      ALLOCATE(rnd_seed(rnd_size))

      ! Initialize randum number generator with an integer sequence depending on the ensemble member ID
      DO i = 1, rnd_size
        rnd_seed(i) = (gribout_config(1)%perturbationNumber - 1) * 100 + i
      ENDDO
      CALL RANDOM_SEED(PUT=rnd_seed)


      ! Apply perturbations to physics tuning parameters

      CALL RANDOM_NUMBER(rnd_num)
      ! perturbations for gkwake and gkdrag must be anticorrelated
      tune_gkwake = tune_gkwake + 2._wp*(rnd_num-0.5_wp)*range_gkwake
      tune_gkdrag = tune_gkdrag - 2._wp*(rnd_num-0.5_wp)*range_gkdrag

      CALL RANDOM_NUMBER(rnd_num)
      ! perturbations for zvz0i and entrorg must be anticorrelated
      tune_zvz0i   = tune_zvz0i   + 2._wp*(rnd_num-0.5_wp)*range_zvz0i
      tune_entrorg = tune_entrorg - 2._wp*(rnd_num-0.5_wp)*range_entrorg

      CALL RANDOM_NUMBER(rnd_num)
      tune_gfluxlaun = tune_gfluxlaun + 2._wp*(rnd_num-0.5_wp)*range_gfluxlaun

      CALL RANDOM_NUMBER(rnd_num)
      ! Perturbation for CAPE diurnal cycle must be one-sided and is chosen to be zero for half of the
      ! ensemble members because applying this perturbation degrades scores in some cases
      tune_capdcfac_et = tune_capdcfac_et + MAX(0._wp,2._wp*(rnd_num-0.5_wp))*range_capdcfac_et

      CALL RANDOM_NUMBER(rnd_num)
      tune_box_liq = tune_box_liq + 2._wp*(rnd_num-0.5_wp)*range_box_liq

      CALL RANDOM_NUMBER(rnd_num)
      turbdiff_config(1:max_dom)%tkhmin = turbdiff_config(1:max_dom)%tkhmin + 2._wp*(rnd_num-0.5_wp)*range_tkhmin

      CALL RANDOM_NUMBER(rnd_num)
      turbdiff_config(1:max_dom)%tkmmin = turbdiff_config(1:max_dom)%tkmmin + 2._wp*(rnd_num-0.5_wp)*range_tkmmin

      CALL RANDOM_NUMBER(rnd_num)
      rnd_fac = range_rlam_heat**(2._wp*(rnd_num-0.5_wp))
      ! The product rlam_heat*rat_sea must stay unchanged in order to keep evaporation over the oceans constant
      turbdiff_config(1:max_dom)%rlam_heat = turbdiff_config(1:max_dom)%rlam_heat * rnd_fac
      turbdiff_config(1:max_dom)%rat_sea   = turbdiff_config(1:max_dom)%rat_sea   / rnd_fac


      ! control output
      WRITE(message_text,'(2f8.4,e11.4)') tune_gkwake, tune_gkdrag, tune_gfluxlaun
      CALL message('Perturbed values, gkwake, gkdrag, gfluxlaun', TRIM(message_text))

      WRITE(message_text,'(3f8.4,e11.4)') tune_box_liq, tune_zvz0i, tune_capdcfac_et, tune_entrorg
      CALL message('Perturbed values, box_liq, zvz0i, capdcfac_et, entrorg', TRIM(message_text))

      WRITE(message_text,'(3f8.4,f8.3)') turbdiff_config(1)%tkhmin, turbdiff_config(1)%tkmmin, &
        turbdiff_config(1)%rlam_heat , turbdiff_config(1)%rat_sea
      CALL message('Perturbed values, tkhmin, tkmmin, rlam_heat, rat_sea', TRIM(message_text))

      DEALLOCATE(rnd_seed)
    ENDIF

  END SUBROUTINE configure_ensemble_pert

END MODULE mo_ensemble_pert_config
