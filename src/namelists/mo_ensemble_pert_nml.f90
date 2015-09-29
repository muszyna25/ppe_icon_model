!>
!! @brief Namelist for perturbing nwp physics
!!
!! These subroutines are called by read_atmo_namelists and set the ranges
!! for ensemble physics perturbations
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
MODULE mo_ensemble_pert_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_ensemble_pert_config,ONLY: config_range_gkwake    => range_gkwake,    &
    &                               config_range_gkdrag    => range_gkdrag,    &
    &                               config_range_gfluxlaun => range_gfluxlaun, &
    &                               config_range_zvz0i     => range_zvz0i,     &  
    &                               config_range_entrorg   => range_entrorg,   &  
    &                               config_range_capdcfac_et => range_capdcfac_et, &  
    &                               config_range_minsnowfrac => range_minsnowfrac, &  
    &                               config_range_rhebc     => range_rhebc,     &  
    &                               config_range_texc      => range_texc,      &  
    &                               config_range_box_liq   => range_box_liq,   &  
    &                               config_range_tkhmin    => range_tkhmin,    &  
    &                               config_range_tkmmin    => range_tkmmin,    &  
    &                               config_range_rlam_heat => range_rlam_heat, &
    &                               config_use_ensemble_pert => use_ensemble_pert

  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_ensemble_pert_namelist


  !--------------------------------------!
  ! ensemble_pert_nml namelist variables !
  !--------------------------------------!

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

  REAL(wp) :: &                    !< RH thresholds for evaporation below cloud base
    &  range_rhebc

  REAL(wp) :: &                    !< Excess value for temperature used in test parcel ascent
    &  range_texc

  REAL(wp) :: &                    !< Minimum value to which the snow cover fraction is artificially reduced
    &  range_minsnowfrac           !  in case of melting show (in case of idiag_snowfrac = 20/30/40)

  REAL(wp) :: &                    !< Box width for liquid clouds assumed in the cloud cover scheme
    &  range_box_liq                ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Minimum vertical diffusion for heat/moisture 
    &  range_tkhmin

  REAL(wp) :: &                    !< Minimum vertical diffusion for momentum 
    &  range_tkmmin

  REAL(wp) :: &                    !< Laminar transport resistance parameter 
    &  range_rlam_heat

  LOGICAL :: use_ensemble_pert     !< main switch

  NAMELIST/ensemble_pert_nml/ use_ensemble_pert, range_gkwake, range_gkdrag, range_gfluxlaun, range_zvz0i, &
    &                         range_entrorg, range_capdcfac_et, range_box_liq, range_tkhmin, range_tkmmin, &
    &                         range_rlam_heat, range_rhebc, range_texc, range_minsnowfrac

CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for NWP ensemble perturbations. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP ensemble perturbations
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !! Initial Revision by Guenther Zaengl, DWD (2015-04-23)
  !!
  SUBROUTINE read_ensemble_pert_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_ensemble_pert_nml: read_ensemble_pert_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    ! Ranges for ensemble perturbations:

    ! SSO tuning
    range_gkwake     = 1._wp/6._wp   ! low-level blocking parameter
    range_gkdrag     = 0.02_wp       ! parameter for wave drag deposition in the atmosphere
    !
    ! GWD tuning
    range_gfluxlaun  = 0.50e-3_wp   ! scaling parameter for GWD flux production
    !
    ! grid scale microphysics
    range_zvz0i      = 0.125_wp     ! scaling for cloud ice sedimentation speed
    !
    ! convection
    range_entrorg    = 0.125e-3_wp  ! entrainment parameter for deep convection
    range_capdcfac_et = 0.125_wp    ! fraction of CAPE diurnal cycle correction applied in the extratropics
    range_rhebc      = 0.05_wp      ! RH thresholds for evaporation below cloud base
    range_texc       = 0.025_wp     ! Excess value for temperature used in test parcel ascent
    !
    ! cloud cover
    range_box_liq    = 0.01_wp      ! box width scale of liquid clouds
    !
    ! turbulence scheme
    range_tkhmin     = 0.15_wp      ! minimum vertical diffusion for heat/moisture
    range_tkmmin     = 0.15_wp      ! minimum vertical diffusion for momentum
    range_rlam_heat  = 1.5_wp       ! multiplicative change of laminar transport resistance parameter
                                    ! (compensated by an inverse change of rat_sea)
    !
    ! snow cover diagnosis
    range_minsnowfrac = 0.05_wp     ! Minimum value to which the snow cover fraction is artificially reduced
                                    ! in case of melting show (in case of idiag_snowfrac = 20/30/40)

    use_ensemble_pert = .FALSE.     ! Usage of ensemble perturbations must be turned on explicitly

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('ensemble_pert_nml')
      READ(funit,NML=ensemble_pert_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('ensemble_pert_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ensemble_pert_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, ensemble_pert_nml)             ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ensemble_pert_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    config_range_gkwake       = range_gkwake
    config_range_gkdrag       = range_gkdrag
    config_range_gfluxlaun    = range_gfluxlaun
    config_range_zvz0i        = range_zvz0i
    config_range_entrorg      = range_entrorg
    config_range_capdcfac_et  = range_capdcfac_et
    config_range_rhebc        = range_rhebc
    config_range_texc         = range_texc
    config_range_minsnowfrac  = range_minsnowfrac
    config_range_box_liq      = range_box_liq
    config_range_tkhmin       = range_tkhmin
    config_range_tkmmin       = range_tkmmin
    config_range_rlam_heat    = range_rlam_heat
    config_use_ensemble_pert  = use_ensemble_pert


    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=ensemble_pert_nml)                    
      CALL store_and_close_namelist(funit, 'ensemble_pert_nml')             
    ENDIF

    !--------------------------------------------------------
    ! 7. write the contents of the namelist to an ASCII file
    !--------------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=ensemble_pert_nml)


  END SUBROUTINE read_ensemble_pert_namelist

END MODULE mo_ensemble_pert_nml
