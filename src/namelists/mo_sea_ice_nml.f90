!>
!! @brief namelist setup for the sea-ice model
!!
!! Namelist setup for the sea-ice model
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Einar Olason, MPI-M (2013-01-29)
!!
!!
!! @par Revision History
!! New file based on mo_lnd_jsbach_nml.f90 by Einar Olason, MPI-M (2013-01-29)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_sea_ice_nml

  USE mo_kind,                ONLY: wp
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist, &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_exception,           ONLY: finish, message
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_run_config,          ONLY: dtime

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_sea_ice_namelist

  INTEGER,PUBLIC :: kice                !< Number of ice classes
  INTEGER,PUBLIC :: i_ice_therm         !< Thermodynamic model switch:
                                        ! 0: No thermodynamics for test, full ice model initialization and cleanup
                                        ! 1: Zero layers
                                        ! 2: Winton's two layer model
                                        ! 3: Zero layers with analytical fluxes
                                        ! 4: Zero layers - only calculates surface temperature
  INTEGER,PUBLIC :: i_ice_albedo        !< Albedo model (no albedo model implemented yet)
  INTEGER,PUBLIC :: i_ice_dyn           !< Dynamical model switch:
                                        ! 0: No dynamics
                                        ! 1: FEM dynamics (from AWI)
  INTEGER,PUBLIC :: i_Qio_type          !< Methods to calculate ice-ocean heatflux
                                        ! 1: Proportional to temperature difference below ice-covered part only
                                        ! 2: Proportional to speed difference btwn. ice and ocean
                                        ! 3: Proportional to temperature difference of whole grid-cell (like MPI-OM)

  INTEGER ,PUBLIC :: i_therm_slo = 0    ! cleanup switch for thermodynamic model - 0: corrected model, old style; 1: cleaned model

  REAL(wp),PUBLIC :: hnull              !< Hibler's h_0 for new ice formation
  REAL(wp),PUBLIC :: hmin               !< Minimum ice thickness allowed in the model
  REAL(wp),PUBLIC :: ramp_wind          !< Time (in days) that the wind stress is increased over.
                                        !  This is only necessary for runs which start off from a
                                        !  still ocean, i.e. not re-start
  REAL(wp),PUBLIC :: hci_layer          !< Thickness of stabilizing constant heat capacity layer
  REAL(wp),PUBLIC :: leadclose_1        !< Hibler's leadclose parameter for lateral melting
  REAL(wp),PUBLIC :: leadclose_2n       !< MPIOM's leadclose parameters par_3/par_2 to push new ice together

  ! some analytic initialization parameters
  REAL(wp),PUBLIC :: init_analytic_temp_under_ice= -1.6_wp
  REAL(wp),PUBLIC :: init_analytic_conc_param    = 0.9_wp
  REAL(wp),PUBLIC :: init_analytic_hi_param      = 2.0_wp
  REAL(wp),PUBLIC :: init_analytic_hs_param      = 0.2_wp
  REAL(wp),PUBLIC :: t_heat_base                 = -5.0_wp  !  arbitrary temperature used as basis for
                                                            !  heat content in energy calculations
  LOGICAL, PUBLIC :: use_constant_tfreez         = .TRUE.   !  constant freezing temperature for ocean water (=Tf)
  LOGICAL, PUBLIC :: use_IceInitialization_fromTemperature = .true.
  LOGICAL, PUBLIC :: stress_ice_zero             = .TRUE.   !  set stress below sea ice to zero
  LOGICAL, PUBLIC :: use_calculated_ocean_stress = .FALSE.  !  calculate ocean stress instead of reading from OMIP
  LOGICAL, PUBLIC :: use_no_flux_gradients       = .TRUE.   !  simplified ice_fast without flux gradients

  ! Rheology
  LOGICAL,PUBLIC  :: ice_VP_rheology=.false.  ! VP=.true., EVP=.false.
  LOGICAL,PUBLIC  :: ice_EVP1=.true.      ! standard EVP .true., modified EVP .false.
  REAL(wp),PUBLIC :: delta_min            ! (1/s) Limit for minimum divergence; valid for both VP and EVP
  ! VP, not implemented
!  INTEGER,PUBLIC  :: ice_VP_iter=10       ! The numberof Pickard iterations
!  REAL(wp),PUBLIC :: ice_VP_soltol=1.0e-10
  ! EVP
  INTEGER,PUBLIC  :: evp_rheol_steps   ! the number of sybcycling steps in EVP
  REAL(wp),PUBLIC :: Tevp_inv          ! Time scale of the standard EVP, inverse of the relaxation time
                                       ! Tevp_inv=3/dtime is the default value
!  REAL(wp),PUBLIC :: alpha_evp        ! Parameters  of modified EVP formulation in Bouillon (2013)
!  REAL(wp),PUBLIC :: beta_evp
  ! Motion
  INTEGER ,PUBLIC :: ice_advection         ! type of ice advection
  REAL(wp),PUBLIC :: theta_io          ! ice/ocean rotation angle. Implemented in EVP, can beadded to VP

  INTEGER         :: iunit

  NAMELIST /sea_ice_nml/ &
    &  kice, &
    &  i_ice_therm, &
    &  i_ice_albedo, &
    &  i_ice_dyn, &
    &  hnull, & 
    &  hmin, &
    &  ramp_wind, &
    &  i_Qio_type, &
    &  i_therm_slo, &
    &  hci_layer, &
    &  leadclose_1, &
    &  leadclose_2n, &
    &  t_heat_base, &
    &  use_IceInitialization_fromTemperature, &
    &  use_constant_tfreez, &
    &  stress_ice_zero,    &
    &  use_calculated_ocean_stress, &
    &  use_no_flux_gradients, &
    &  init_analytic_temp_under_ice, &
    &  init_analytic_conc_param , &
    &  init_analytic_hi_param, &
    &  init_analytic_hs_param, &

    &  ice_VP_rheology, &
    &  ice_EVP1, &
    &  delta_min, &
!    &  ice_VP_iter, &
!    &  ice_VP_soltol, &
    &  evp_rheol_steps, &
    &  Tevp_inv, &
!    &  alpha_evp, &
!    &  beta_evp, &
    &  ice_advection, &
    &  theta_io

CONTAINS
  !>
  !!
  SUBROUTINE read_sea_ice_namelist(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_sea_ice_nml:read_sea_ice_namelist'

    !------------------------------------------------------------------
    ! Set default values
    !------------------------------------------------------------------
    kice        = 1
    i_ice_therm = 1
    i_ice_albedo= 1
    i_ice_dyn   = 0
    i_Qio_type  = 3

    hnull       = 0.5_wp
    hmin        = 0.05_wp
    hci_layer   = 0.10_wp
    leadclose_1 = 0.5_wp
    leadclose_2n = 0.0_wp

    ramp_wind    = 1.0_wp

    delta_min    = 2.0e-11_wp ! Hibler, Hunke normally use 2.0e-9, which does much stronger limiting
    evp_rheol_steps = 120
    Tevp_inv     = 3.0_wp/dtime
!    alpha_evp=500            ! Parameters  of modified EVP formulation in Bouillon (2013)
!    beta_evp=1000
    ice_advection=1           ! 1 switches on FCT advection, and
                                  ! 0 switches on Backward Euler advection
    theta_io=0._wp    !0.436

    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('sea_ice_nml')
      READ(funit,NML=sea_ice_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------
    ! Read user's (new) specifications (done by all MPI processes)
    !------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('sea_ice_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, sea_ice_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (positioned)
      READ (nnml, sea_ice_nml)                                        ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, sea_ice_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------------
    ! Sanity Check
    !------------------------------------------------------------------
    IF (kice /= 1 ) THEN
      CALL finish(TRIM(routine), 'Currently, kice must be 1.')
    END IF

    IF (i_ice_therm < 0 .OR. i_ice_therm > 4) THEN
      CALL finish(TRIM(routine), 'i_ice_therm must be between 0 and 4.')
    END IF

    IF (i_ice_therm == 2) THEN
      CALL finish(TRIM(routine), 'i_ice_therm = 2 is not allowed - Winton thermodynamics not active anymore')
    END IF

    IF (i_ice_albedo < 1 .OR. i_ice_albedo > 2 ) THEN
      CALL finish(TRIM(routine), 'i_ice_albedo must be either 1 or 2.')
    END IF

    IF (i_ice_dyn < 0 .OR. i_ice_dyn > 1) THEN
      CALL finish(TRIM(routine), 'i_ice_dyn must be either 0 or 1.')
    END IF

    IF (i_ice_dyn == 1 ) THEN
      CALL message(TRIM(routine), 'WARNING: i_ice_dyn is 1 - BUT SEA ICE DYNAMICS INCLUDE ERRORS')
      ! TODO: This can be changed when we start advecting T1 and T2
   !  CALL message(TRIM(routine), 'WARNING: i_ice_therm set to 1 because i_ice_dyn is 1')
   !  i_ice_therm = 1   !  no Winton thermodynamics allowed, switched off by default

      ! When using routine ice_ocean_stress, ocean stress below sea ice is considered accordingly
      CALL message(TRIM(routine), 'WARNING: stress_ice_zero=FALSE because i_ice_dyn is 1')
      stress_ice_zero = .TRUE.
    ENDIF


    IF (i_Qio_type < 1 .OR. i_Qio_type > 3) THEN
      CALL finish(TRIM(routine), 'i_Qio_type must be either 1,2 or 3.')
    END IF

    IF (i_ice_dyn == 0 .AND. i_Qio_type == 2) THEN
      CALL message(TRIM(routine), 'i_Qio_type set to 3 because i_ice_dyn is 0')
      i_Qio_type = 3
    ENDIF

    IF (hmin > hnull) THEN
      CALL message(TRIM(routine), 'hmin cannot be larger than hnull')
      CALL message(TRIM(routine), 'setting hmin to hnull')
      hmin = hnull
    ENDIF

    IF (ramp_wind <= 0) THEN
      CALL message(TRIM(routine), 'ramp_wind cannot be smaller than or equal to zero')
      CALL message(TRIM(routine), 'setting ramp_wind to TINY(1._wp)')
      ramp_wind = TINY(1._wp)
    ENDIF

    IF (hci_layer < 0) THEN
      CALL message(TRIM(routine), 'hci_layer < 0, setting it equal to zero')
    ENDIF
      
  ! IF (i_ice_dyn == 1 ) THEN
  !   i_ice_dyn = 0
  !   CALL warning(method_name,"Disable sea-ice dynamics. It does not work.")
  ! ENDIF

    !------------------------------------------------------------------
    ! Store the namelist for restart
    !------------------------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=sea_ice_nml)
      CALL store_and_close_namelist(funit, 'sea_ice_nml')
    ENDIF

    !------------------------------------------------------------------
    ! Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) WRITE(nnml_output,nml=sea_ice_nml)

    !------------------------------------------------------------------
    ! Fill the configuration state
    !------------------------------------------------------------------

  END SUBROUTINE read_sea_ice_namelist

END MODULE mo_sea_ice_nml
