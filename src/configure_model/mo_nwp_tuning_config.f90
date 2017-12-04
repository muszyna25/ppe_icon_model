!>
!! @brief Tuning and/or perturbing nwp physics
!!
!! configuration setup for NWP physics tuning
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2014-09-25)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nwp_tuning_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tune_gkwake
  PUBLIC :: tune_gkdrag
  PUBLIC :: tune_gfrcrit
  PUBLIC :: tune_grcrit
  PUBLIC :: tune_gfluxlaun
  PUBLIC :: tune_zceff_min
  PUBLIC :: tune_v0snow
  PUBLIC :: tune_zvz0i
  PUBLIC :: tune_entrorg
  PUBLIC :: tune_capdcfac_et
  PUBLIC :: tune_rhebc_land
  PUBLIC :: tune_rhebc_ocean
  PUBLIC :: tune_rcucov
  PUBLIC :: tune_rhebc_land_trop
  PUBLIC :: tune_rhebc_ocean_trop
  PUBLIC :: tune_rcucov_trop
  PUBLIC :: tune_texc
  PUBLIC :: tune_qexc
  PUBLIC :: tune_minsnowfrac
  PUBLIC :: tune_box_liq
  PUBLIC :: tune_dust_abs
  PUBLIC :: itune_albedo
  PUBLIC :: lcalib_clcov
  PUBLIC :: max_freshsnow_inc


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for physics tuning
  !!--------------------------------------------------------------------------

!  TYPE :: t_nwp_tuning_config

    ! namelist variables
  REAL(wp) :: &                    !< low level wake drag constant
    &  tune_gkwake(max_dom)

  REAL(wp) :: &                    !< gravity wave drag constant
    &  tune_gkdrag(max_dom)

  REAL(wp) :: &                    !< critical Froude number in SSO scheme
    &  tune_gfrcrit(max_dom)

  REAL(wp) :: &                    !< critical Richardson number in SSO scheme
    &  tune_grcrit(max_dom)

  REAL(wp) :: &                    !< total launch momentum flux in each azimuth (rho_o x F_o)
    &  tune_gfluxlaun

  REAL(wp) :: &                    !< Minimum value for sticking efficiency
    &  tune_zceff_min

  REAL(wp) :: &                    !< factor in the terminal velocity for snow
    &  tune_v0snow

  REAL(wp) :: &                    !< Terminal fall velocity of ice 
    &  tune_zvz0i

  REAL(wp) :: &                    !< Entrainment parameter for deep convection valid at dx=20 km 
    &  tune_entrorg

  REAL(wp) :: &                    !< Fraction of CAPE diurnal cycle correction applied in the extratropics
    &  tune_capdcfac_et            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< RH threshold for onset of evaporation below cloud base over land
    &  tune_rhebc_land

  REAL(wp) :: &                    !< RH threshold for onset of evaporation below cloud base over sea
    &  tune_rhebc_ocean

  REAL(wp) :: &                    !< Convective area fraction
    &  tune_rcucov

  REAL(wp) :: &                    !< RH threshold for onset of evaporation below cloud base over tropical land
    &  tune_rhebc_land_trop        !  (relevant only if smaller than rhebc_land)

  REAL(wp) :: &                    !< RH threshold for onset of evaporation below cloud base over tropical sea
    &  tune_rhebc_ocean_trop       !  (relevant only if smaller than rhebc_ocean)

  REAL(wp) :: &                    !< Convective area fraction in the tropics
    &  tune_rcucov_trop            !  (relevant only if smaller than rcucov)

  REAL(wp) :: &                    !< Excess value for temperature used in test parcel ascent
    &  tune_texc

  REAL(wp) :: &                    !< Excess fraction of grid-scale QV used in test parcel ascent
    &  tune_qexc

  REAL(wp) :: &                    !< Minimum value to which the snow cover fraction is artificially reduced
    &  tune_minsnowfrac            !  in case of melting show (in case of idiag_snowfrac = 20/30/40)

  REAL(wp) :: &                    !< Box width for liquid clouds assumed in the cloud cover scheme
    &  tune_box_liq                ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Tuning factor for enhanced LW absorption of mineral dust in the Saharan region
    &  tune_dust_abs               !

  INTEGER :: &                     !< (MODIS) albedo tuning
    &  itune_albedo                ! 1: dimmed Sahara
                                   ! 2: dimmed Sahara and brighter Antarctica

  LOGICAL :: &                     ! cloud cover calibration over land points
    &  lcalib_clcov

  REAL(wp) :: &                    !< maximum allowed positive freshsnow increment
    &  max_freshsnow_inc

!  END TYPE t_nwp_tuning_config


!CONTAINS


END MODULE mo_nwp_tuning_config
