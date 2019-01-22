!>
!! @brief Namelist for tuning and/or perturbing nwp physics
!!
!! These subroutines are called by read_atmo_namelists and do some 
!! nwp physics tuning 
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
MODULE mo_nwp_tuning_nml

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_nwp_tuning_config,   ONLY: config_tune_gkwake    => tune_gkwake,    &
    &                               config_tune_gkdrag    => tune_gkdrag,    &
    &                               config_tune_gfrcrit   => tune_gfrcrit,   &
    &                               config_tune_grcrit    => tune_grcrit,    &
    &                               config_tune_gfluxlaun => tune_gfluxlaun, &
    &                               config_tune_zceff_min => tune_zceff_min, &
    &                               config_tune_v0snow    => tune_v0snow,    &
    &                               config_tune_zvz0i     => tune_zvz0i,     &
    &                               config_tune_icesedi_exp => tune_icesedi_exp, &
    &                               config_tune_entrorg     => tune_entrorg,     &  
    &                               config_tune_capdcfac_et => tune_capdcfac_et, &
    &                               config_tune_capdcfac_tr => tune_capdcfac_tr, &
    &                               config_tune_rhebc_land  => tune_rhebc_land,  &  
    &                               config_tune_rhebc_ocean => tune_rhebc_ocean, &  
    &                               config_tune_rcucov      => tune_rcucov,      &  
    &                               config_tune_rhebc_land_trop  => tune_rhebc_land_trop,  &  
    &                               config_tune_rhebc_ocean_trop => tune_rhebc_ocean_trop, &  
    &                               config_tune_rcucov_trop      => tune_rcucov_trop,      &  
    &                               config_tune_texc        => tune_texc,        &  
    &                               config_tune_qexc        => tune_qexc,        &  
    &                               config_tune_minsnowfrac => tune_minsnowfrac, &  
    &                               config_tune_box_liq   => tune_box_liq,       &  
    &                               config_tune_box_liq_asy => tune_box_liq_asy, &
    &                               config_tune_dust_abs  => tune_dust_abs,      &  
    &                               config_itune_albedo   => itune_albedo,       &
    &                               config_lcalib_clcov   => lcalib_clcov,       &
    &                               config_max_freshsnow_inc => max_freshsnow_inc 
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_nwp_tuning_namelist


  !-----------------------------------!
  ! nwp_tuning_nml namelist variables !
  !-----------------------------------!

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

  REAL(wp) :: &                    !< Exponent for density correction of cloud ice sedimentation
    &  tune_icesedi_exp

  REAL(wp) :: &                    !< Entrainment parameter for deep convection valid at dx=20 km 
    &  tune_entrorg

  REAL(wp) :: &                    !< Fraction of CAPE diurnal cycle correction applied in the extratropics
    &  tune_capdcfac_et            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< Fraction of CAPE diurnal cycle correction applied in the tropics
    &  tune_capdcfac_tr            ! (relevant only if icapdcycl = 3)

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

  REAL(wp) :: &                    !< Asymmetry factor liquid cloud parameterization
    &  tune_box_liq_asy            ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Tuning factor for enhanced LW absorption of mineral dust in the Saharan region
    &  tune_dust_abs               !

  INTEGER :: &                     !< (MODIS) albedo tuning
    &  itune_albedo                ! 0: no tuning
                                   ! 1: dimmed Sahara
                                   ! 2: dimmed Sahara and brighter Antarctica

  LOGICAL :: &                     ! cloud cover calibration over land points
    &  lcalib_clcov

  REAL(wp) :: &                    !< maximum allowed positive freshsnow increment
    &  max_freshsnow_inc

  NAMELIST/nwp_tuning_nml/ tune_gkwake, tune_gkdrag, tune_gfluxlaun,        &
    &                      tune_zceff_min, tune_v0snow, tune_zvz0i,         &
    &                      tune_entrorg, itune_albedo, max_freshsnow_inc,   &
    &                      tune_capdcfac_et, tune_box_liq, tune_rhebc_land, &
    &                      tune_rhebc_ocean, tune_rcucov, tune_texc,        &
    &                      tune_qexc, tune_minsnowfrac,tune_rhebc_land_trop,&
    &                      tune_rhebc_ocean_trop, tune_rcucov_trop,         &
    &                      tune_dust_abs, tune_gfrcrit, tune_grcrit,        &
    &                      lcalib_clcov, tune_box_liq_asy, tune_capdcfac_tr,&
    &                      tune_icesedi_exp

CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for NWP physics tuning. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP physics tuning
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !! Initial Revision by Daniel Reinert, DWD (2014-09-25)
  !!
  SUBROUTINE read_nwp_tuning_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit, jg

    REAL(wp) :: gkwake_def, gkdrag_def, gfrcrit_def, grcrit_def

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nwp_tuning_nml: read_tuning_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    ! Comment: In case we want to draw from a normal distribution, the namelist 
    ! parameters could be extended to arrays of size 2. The first value is the mean, 
    ! while the second one is the standard deviation. 

    ! SSO tuning
    gkwake_def  = 1.5_wp       ! original COSMO value 0.5
    gkdrag_def  = 0.075_wp     ! original COSMO value 0.075
    gfrcrit_def = 0.4_wp       ! original COSMO value 0.5
    grcrit_def  = 0.25_wp      ! original COSMO value 0.25

    tune_gkwake(:)  = gkwake_def
    tune_gkdrag(:)  = gkdrag_def
    tune_gfrcrit(:) = gfrcrit_def
    tune_grcrit(:)  = grcrit_def
    !
    ! GWD tuning
    tune_gfluxlaun  = 2.50e-3_wp   ! original IFS value 3.75e-3
    !
    ! grid scale microphysics
    tune_zceff_min   = 0.01_wp
    tune_v0snow      = 25.0_wp      ! previous ICON value was 20
    tune_zvz0i       = 1.25_wp      ! original value of Heymsfield+Donner 1990: 3.29
    tune_icesedi_exp = 0.33_wp      ! exponent for density correction of cloud ice sedimentation

    !
    ! convection
    !

    !> entrainment parameter for deep convection:
    tune_entrorg     = 1.95e-3_wp   

    !> fraction of CAPE diurnal cycle correction applied in the extratropics
    tune_capdcfac_et = 0.5_wp

    !> fraction of CAPE diurnal cycle correction applied in the tropics
    tune_capdcfac_tr = 0.5_wp

    !> RH threshold for onset of evaporation below cloud base over land (original IFS value 0.7)
    tune_rhebc_land  = 0.75_wp

    !> RH threshold for onset of evaporation below cloud base over sea (original IFS value 0.9)
    tune_rhebc_ocean = 0.85_wp

    !> Convective area fraction used for computing evaporation below cloud base (original IFS value 0.05)
    tune_rcucov      = 0.05_wp

    !> Excess value for temperature used in test parcel ascent (K) (original IFS value 0.2 K)
    tune_texc        = 0.125_wp

    !> Excess fraction of grid-scale QV used in test parcel ascent
    !  (original IFS value 0.1 g/kg independent of grid-scale QV))
    tune_qexc        = 1.25e-2_wp

    ! The following switches allow separate tuning for evaporation below cloud base in the tropics
    tune_rhebc_land_trop  = 0.70_wp
    tune_rhebc_ocean_trop = 0.76_wp
    tune_rcucov_trop      = 0.03_wp

    !
    ! snow cover diagnosis
    tune_minsnowfrac = 0.20_wp     ! Minimum value to which the snow cover fraction is artificially reduced
                                   ! in case of melting show (in case of idiag_snowfrac = 20/30/40)
    !
    ! cloud cover
    tune_box_liq     = 0.05_wp     ! box width scale of liquid clouds
    tune_box_liq_asy = 2.5_wp      ! asymmetry factor for liquid cloud parameterization
    lcalib_clcov     = .TRUE.      ! use calibration of layer-wise cloud cover diagnostics over land

    tune_dust_abs   = 0._wp        ! no tuning of LW absorption of mineral dust
    itune_albedo    = 0            ! original (measured) albedo
    !
    ! IAU increment tuning
    max_freshsnow_inc = 0.025_wp   ! maximum allowed positive freshsnow increment

    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, nwp_tuning_nml)   ! write defaults to temporary text file
    END IF

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('nwp_tuning_nml')
      READ(funit,NML=nwp_tuning_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('nwp_tuning_nml', STATUS=istat)

    SELECT CASE (istat)
    CASE (POSITIONED)

      ! Set array parameters to dummy values to determine which ones are actively set in the namelist
      tune_gkwake(:)  = -1._wp
      tune_gkdrag(:)  = -1._wp
      tune_gfrcrit(:) = -1._wp
      tune_grcrit(:)  = -1._wp

      READ (nnml, nwp_tuning_nml)    ! overwrite default settings

      ! Reset first element to default values if not specified in the namelist
      IF (tune_gkwake(1)  < 0._wp) tune_gkwake(1)  = gkwake_def
      IF (tune_gkdrag(1)  < 0._wp) tune_gkdrag(1)  = gkdrag_def
      IF (tune_gfrcrit(1) < 0._wp) tune_gfrcrit(1) = gfrcrit_def
      IF (tune_grcrit(1)  < 0._wp) tune_grcrit(1)  = grcrit_def

      ! Fill remaining array elements with entry of parent domain if not specified in the namelist
      DO jg = 2, max_dom
        IF (tune_gkwake(jg)  < 0._wp) tune_gkwake(jg)  = tune_gkwake(jg-1)
        IF (tune_gkdrag(jg)  < 0._wp) tune_gkdrag(jg)  = tune_gkdrag(jg-1)
        IF (tune_gfrcrit(jg) < 0._wp) tune_gfrcrit(jg) = tune_gfrcrit(jg-1)
        IF (tune_grcrit(jg)  < 0._wp) tune_grcrit(jg)  = tune_grcrit(jg-1)
      ENDDO

      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, nwp_tuning_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    config_tune_gkwake           = tune_gkwake
    config_tune_gkdrag           = tune_gkdrag
    config_tune_gfrcrit          = tune_gfrcrit
    config_tune_grcrit           = tune_grcrit
    config_tune_gfluxlaun        = tune_gfluxlaun
    config_tune_zceff_min        = tune_zceff_min 
    config_tune_v0snow           = tune_v0snow
    config_tune_zvz0i            = tune_zvz0i
    config_tune_icesedi_exp      = tune_icesedi_exp
    config_tune_entrorg          = tune_entrorg
    config_tune_capdcfac_et      = tune_capdcfac_et
    config_tune_capdcfac_tr      = tune_capdcfac_tr
    config_tune_rhebc_land       = tune_rhebc_land
    config_tune_rhebc_ocean      = tune_rhebc_ocean
    config_tune_rcucov           = tune_rcucov
    config_tune_rhebc_land_trop  = tune_rhebc_land_trop
    config_tune_rhebc_ocean_trop = tune_rhebc_ocean_trop
    config_tune_rcucov_trop      = tune_rcucov_trop
    config_tune_texc             = tune_texc
    config_tune_qexc             = tune_qexc
    config_tune_minsnowfrac      = tune_minsnowfrac
    config_tune_box_liq          = tune_box_liq
    config_tune_box_liq_asy      = tune_box_liq_asy
    config_tune_dust_abs         = tune_dust_abs
    config_itune_albedo          = itune_albedo
    config_lcalib_clcov          = lcalib_clcov
    config_max_freshsnow_inc     = max_freshsnow_inc


    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=nwp_tuning_nml)                    
      CALL store_and_close_namelist(funit, 'nwp_tuning_nml')             
    ENDIF

    !--------------------------------------------------------
    ! 7. write the contents of the namelist to an ASCII file
    !--------------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nwp_tuning_nml)


  END SUBROUTINE read_nwp_tuning_namelist

END MODULE mo_nwp_tuning_nml
