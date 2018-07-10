!>
!! Data type containing basic control variables for a model integration. 
!!
!! Note that in a coupled simulation, each model component (e.g.,
!! atmosphere, ocean) will have its own run-configuration.
!! 
!! @author Kristina Froehlich, MPI-M (2011-07-12)
!! @author Hui Wan, MPI-M (2011-07-12)
!!
!! @par Revision History
!! Initial version by Hui Wan (MPI-M, 2011-07-05)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_run_config

  USE mo_kind,           ONLY: wp
  USE mo_memory_log,     ONLY: memory_log_initialize
  USE mo_impl_constants, ONLY: MAX_DOM, IHELDSUAREZ, INWP, IECHAM, ILDF_ECHAM, &
                               IMPIOM, INOFORCING, ILDF_DRY, MAX_CHAR_LENGTH,  &
                               TIMER_MODE_AGGREGATED, TIMER_MODE_DETAILED
  USE mtime,             ONLY: MAX_TIMEDELTA_STR_LEN
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ltestcase, ldynamics, iforcing, lforcing, logmaxrss, logmaxrss_all
  PUBLIC :: ltransport, ntracer, nlev, nlevm1, nlevp1
  PUBLIC :: lart
  PUBLIC :: lvert_nest, num_lev, nshift, nsteps, dtime
  PUBLIC :: ltimer, timers_level, activate_sync_timers, msg_level
  PUBLIC :: iqv, iqc, iqi, iqs, iqr, iqtvar, nqtendphy, iqt, ico2, ich4, in2o, io3
  PUBLIC :: iqni, iqni_nuc, iqg, iqm_max
  PUBLIC :: iqh, iqnh, iqnr, iqns, iqng, iqnc, inccn, ininpot, ininact
  PUBLIC :: iqtke
  PUBLIC :: grid_generatingCenter     ! non-namelist variables
  PUBLIC :: grid_generatingSubcenter  ! non-namelist variables
  PUBLIC :: number_of_grid_used       ! non-namelist variables
  PUBLIC :: test_mode
  PUBLIC :: configure_run
  PUBLIC :: output, t_output_mode, output_mode, max_output_modes
  PUBLIC :: debug_check_level
  PUBLIC :: restart_filename
  PUBLIC :: profiling_output, TIMER_MODE_AGGREGATED, TIMER_MODE_DETAILED
  PUBLIC :: check_uuid_gracefully
  PUBLIC :: modelTimeStep

    ! Namelist variables
    !
    LOGICAL :: ltestcase       !< Run idealized test case
    LOGICAL :: ldynamics       !< Switch on model dynamics
    INTEGER :: iforcing        !< Choice of diabatic forcing

    LOGICAL :: ltransport      !< Switch on tracer transport
    INTEGER :: ntracer         !< Total number of advected tracers

    LOGICAL :: lart            !< switch for ICON-ART (Treatment of Aerosols and Trace Gases)

    LOGICAL :: lvert_nest         !< switch for vertical nesting
    INTEGER :: num_lev  (MAX_DOM) !< number of full levels for each domain
    INTEGER :: nshift   (MAX_DOM) !< half level of parent domain which coincides 
                                  !< with the upper boundary of the current domain jg

    INTEGER :: nsteps          !< number of time steps to integrate
    REAL(wp):: dtime           !< [s] length of a time step
    
    LOGICAL :: ltimer          !< if .TRUE.,  the timer is switched on
    INTEGER :: timers_level    !< what level of timers to run
    LOGICAL :: activate_sync_timers
    INTEGER :: profiling_output = TIMER_MODE_AGGREGATED  !< switch defining the kind of timer output

    LOGICAL :: check_uuid_gracefully !< Flag. If .TRUE. then we give only warnings for non-matching UUIDs
  
    INTEGER :: test_mode = 0   !< 0= run the model, /=0 run in test mode
    INTEGER :: debug_check_level = 10  ! Define debug checks level. This is not related to the debug output in
                                      ! mo_dbg_nml, it only controls the activation of internal checks

    INTEGER :: msg_level       !< how much printout is generated during runtime
    LOGICAL :: logmaxrss
    LOGICAL :: logmaxrss_all

 
    !> output mode (string)
    !  one or multiple of "none", "nml", "totint"
    INTEGER, PARAMETER :: max_output_modes = 5
    CHARACTER(len=32) :: output(max_output_modes)

    ! Non-Namelist variables
    ! These are read from the grid file in mo_model_domimp_patches/read_basic_patch
    ! 
    INTEGER :: grid_generatingCenter   (0:MAX_DOM)   !< patch generating center
    INTEGER :: grid_generatingSubcenter(0:MAX_DOM)   !< patch generating subcenter
    INTEGER :: number_of_grid_used(0:MAX_DOM)  !< Number of grid used (GRIB2 key)
    

    ! Derived variables
    !
    ! Tracer indices of water species
    INTEGER :: iqv        !< water vapor
    INTEGER :: iqc        !< cloud water
    INTEGER :: iqi        !< cloud ice
    !CK>
    INTEGER :: iqni       !< cloud ice number
    INTEGER :: iqni_nuc   !< activated ice nuclei
    INTEGER :: iqg        !< graupel
    !CK<    
    INTEGER :: iqr        !< rain water
    INTEGER :: iqs        !< snow
    INTEGER :: iqtvar     !< total water variance
    INTEGER :: nqtendphy  !< number of water species for which physical tendencies are stored
    INTEGER :: iqm_max    !< highest tracer index carrying a mass-related moisture variable
  
    !For 2 moment microphysics
    INTEGER :: iqh        !<hail
    INTEGER :: iqnh       !<hail number
    INTEGER :: iqnr       !<rain number
    INTEGER :: iqns       !<snow number
    INTEGER :: iqng       !<graupel number
    INTEGER :: iqnc       !<cloud number
    INTEGER :: inccn      !<ccn number
    INTEGER :: ininpot    !<number of aerosol particles which are potential IN
    INTEGER :: ininact    !<number of activated IN

    ! For TKE advection
    INTEGER :: iqtke      !< turbulent kinetic energy

    ! Tracer indices of other species
    INTEGER :: iqt        !< start index of other tracers than hydrometeors
    INTEGER :: ico2       !< CO2
    INTEGER :: ich4       !< CH4
    INTEGER :: in2o       !< N2O
    INTEGER :: io3        !< O3


    INTEGER :: nlev               !< number of full levels for each domain
    INTEGER :: nlevm1             !< number of half levels for each domain without boundaries
    INTEGER :: nlevp1             !< number of half levels for each domain with    boundaries

    LOGICAL :: lforcing           !< diabatic forcing TRUE/FALSE

    !> output mode (logicals)
    !  one or multiple of "none", "nml", "totint", "maxwinds"
    TYPE t_output_mode
      LOGICAL :: l_none, l_nml, l_totint, l_maxwinds
    END TYPE t_output_mode

    TYPE (t_output_mode) output_mode

    !> file name for restart/checkpoint files (containg keyword
    !> substition patterns)
    CHARACTER(len=MAX_CHAR_LENGTH) :: restart_filename

    !> namelist parameter (as raw character string):
    CHARACTER(len=max_timedelta_str_len) :: modelTimeStep

CONTAINS
  !>
  !!
  !! Assign value to components of the run configuration state that have no
  !! corresponding namelist variable.
  !!
  !! Exceptions: grid_generatingCenter, grid_generatingSubcenter and number_of_grid_used 
  !!             are set in mo_model_domimp_patches/read_basic_patch 
  !!
  SUBROUTINE configure_run( )

    CHARACTER(LEN=*),PARAMETER :: routine = 'mo_run_config:configure_run'
    

    !----------------------------
    ! Number of vertical levels

    IF (.NOT.lvert_nest) THEN
      num_lev(:) = nlev
      nshift (:) = 0
    END IF

    nlevm1       = nlev - 1
    nlevp1       = nlev + 1

    !-------------------------------------
    ! Logical switch for diabatic forcing

    SELECT CASE (iforcing)
    CASE(IHELDSUAREZ,INWP,IECHAM,ILDF_ECHAM,IMPIOM)
      lforcing = .TRUE.

    CASE(INOFORCING,ILDF_DRY)
      lforcing = .FALSE.
    END SELECT

    ! activate memory logging if needed
    IF (logmaxrss .OR. logmaxrss_all) CALL memory_log_initialize(logmaxrss_all)

  END SUBROUTINE configure_run
  !-------------------------------------------------------------

 
!  !---------------------------------------
!  !>
!  LOGICAL FUNCTION get_ltestcase()
!    get_ltestcase = ltestcase 
!  END FUNCTION get_ltestcase
!  !---------------------------------------
!  !>
!  LOGICAL FUNCTION get_ldynamics()
!    get_ldynamics = ldynamics 
!  END FUNCTION get_ldynamics
!  !---------------------------------------
!  !>
!  LOGICAL FUNCTION get_ltransport()
!    get_ltransport = ltransport 
!  END FUNCTION get_ltransport
!  !---------------------------------------
!  !>
!  INTEGER FUNCTION get_ntracer()
!    get_ntracer = ntracer 
!  END FUNCTION get_ntracer
!  !---------------------------------------
!  !>
!  INTEGER FUNCTION get_iforcing()
!    get_iforcing = iforcing 
!  END FUNCTION get_iforcing
!  !---------------------------------------
!  !>
!  REAL(wp) FUNCTION get_dtime()
!    get_dtime = dtime 
!  END FUNCTION get_dtime
!  !---------------------------------------
!  !>
!  INTEGER FUNCTION get_nsteps()
!    get_nsteps = nsteps 
!  END FUNCTION get_nsteps
!  !---------------------------------------
!

END MODULE mo_run_config


