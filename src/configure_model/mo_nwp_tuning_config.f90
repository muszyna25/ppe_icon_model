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


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tune_gkwake
  PUBLIC :: tune_gkdrag
  PUBLIC :: tune_gfluxlaun
  PUBLIC :: tune_zceff_min
  PUBLIC :: tune_v0snow
  PUBLIC :: tune_zvz0i
  PUBLIC :: tune_entrorg
  PUBLIC :: tune_capdcfac_et
  PUBLIC :: tune_box_liq
  PUBLIC :: itune_albedo
  PUBLIC :: max_freshsnow_inc


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for physics tuning
  !!--------------------------------------------------------------------------

!  TYPE :: t_nwp_tuning_config

    ! namelist variables
  REAL(wp) :: &                    !< low level wake drag constant
    &  tune_gkwake

  REAL(wp) :: &                    !< gravity wave drag constant
    &  tune_gkdrag

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

  REAL(wp) :: &                    !< Box width for liquid clouds assumed in the cloud cover scheme
    &  tune_box_liq                ! (in case of inwp_cldcover = 1)

  INTEGER :: &                     !< (MODIS) albedo tuning
    &  itune_albedo                ! 1: dimmed Sahara
                                   ! 2: dimmed Sahara and brighter Antarctica


  REAL(wp) :: &                    !< maximum allowed positive freshsnow increment
    &  max_freshsnow_inc

!  END TYPE t_nwp_tuning_config


!CONTAINS


END MODULE mo_nwp_tuning_config
