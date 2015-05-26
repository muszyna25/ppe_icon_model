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
  USE mo_impl_constants, ONLY: MAX_DOM, IHELDSUAREZ, INWP, IECHAM, ILDF_ECHAM, &
                               IMPIOM, INOFORCING, ILDF_DRY, MAX_CHAR_LENGTH,  &
                               TIMER_MODE_AGGREGATED, TIMER_MODE_DETAILED

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ltestcase, ldynamics, iforcing, lforcing
  PUBLIC :: ltransport, ntracer, nlev, nlevm1, nlevp1
  PUBLIC :: lart
  PUBLIC :: lvert_nest, num_lev, nshift, nsteps, dtime
  PUBLIC :: ltimer, timers_level, activate_sync_timers, msg_level
  PUBLIC :: iqv, iqc, iqi, iqs, iqr, iqtvar, nqtendphy, iqt, ico2
  PUBLIC :: iqni, iqni_nuc, iqg, iqm_max
  PUBLIC :: iqh, iqnh, iqnr, iqns, iqng, iqnc, inccn, ininpot, ininact
  PUBLIC :: iqtke
  PUBLIC :: iash1,iash2,iash3,iash4,iash5,iash6                          !Running index for Volcanic Ash in ICON-ART 
  PUBLIC :: iCS137,iI131,iTE132,iZR95,iXE133,iI131g,iI131o,iBA140,iRU103 !Running index for radioactive nuclides  in ICON-ART
  PUBLIC :: iseasa,iseasb,iseasc,iseasa0,iseasb0,iseasc0                 !Running index for sea salt in ICON-ART
  PUBLIC :: idusta,idustb,idustc,idusta0,idustb0,idustc0,idust_act       !Running index for mineral dust in ICON-ART
  PUBLIC :: iTRCHBR3,iTRCH2BR2,iTRBRy                                    !Running index for chemical tracer in ICON-ART - VSLS-BRy
  PUBLIC :: iTRCH4,iTRCO2,iTRCO,iTRH2O,iTRO3                             !Running index for chemical tracer in ICON-ART - CH4-CO-CO2-H2O-O3
  PUBLIC :: iTRSF6l,iTRSF6r,iTRSF6d                                      !Running index for chemical tracer in ICON-ART - SF6
  PUBLIC :: iTRN2O                                                       !Running index for chemical tracer in ICON-ART - N2O
  PUBLIC :: iTR1                                                         !Running index for chemical tracer in ICON-ART - artificial tracer
                                                                         ! RR JS  !Running index for chemical tracer in ICON-ART - VSLS
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
  PUBLIC :: irad_type

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
    INTEGER :: ico2       !< CO2
    INTEGER :: iqt        !< start index of other tracers than hydrometeors
    ! Tracer indices of ICON-ART species
    INTEGER :: iash1        !< Volcanic ash, first class
    INTEGER :: iash2        !< Volcanic ash, second class
    INTEGER :: iash3        !< Volcanic ash, third class
    INTEGER :: iash4        !< Volcanic ash, fourth class
    INTEGER :: iash5        !< Volcanic ash, fifth class
    INTEGER :: iash6        !< Volcanic ash, sixth class
    INTEGER :: iCS137       !< radioactive nuclides
    INTEGER :: iI131        !<
    INTEGER :: iTE132       !< 
    INTEGER :: iZR95        !<
    INTEGER :: iXE133       !<
    INTEGER :: iI131g       !<
    INTEGER :: iI131o       !<
    INTEGER :: iBA140       !< 
    INTEGER :: iRU103       !<
    INTEGER :: iseasa       !< Sea Salt Aerosol Mode A Mass Density
    INTEGER :: iseasb       !< Sea Salt Aerosol Mode B Mass Density
    INTEGER :: iseasc       !< Sea Salt Aerosol Mode C Mass Density
    INTEGER :: iseasa0      !< Sea Salt Aerosol Mode A Number Density
    INTEGER :: iseasb0      !< Sea Salt Aerosol Mode B Number Density
    INTEGER :: iseasc0      !< Sea Salt Aerosol Mode C Number Density
    INTEGER :: idusta       !< Mineral Dust Aerosol Mode A Mass Density
    INTEGER :: idustb       !< Mineral Dust Aerosol Mode B Mass Density
    INTEGER :: idustc       !< Mineral Dust Aerosol Mode C Mass Density
    INTEGER :: idusta0      !< Mineral Dust Aerosol Mode A Number Density
    INTEGER :: idustb0      !< Mineral Dust Aerosol Mode B Number Density
    INTEGER :: idustc0      !< Mineral Dust Aerosol Mode C Number Density
    INTEGER :: idust_act    !< Activated Mineral Dust Number Density
    INTEGER :: iTRCHBR3     !< chemical tracer in ICON-ART
    INTEGER :: iTRCH2BR2    !< chemical tracer in ICON-ART
    INTEGER :: iTRBRy       !< chemical tracer in ICON-ART
    INTEGER :: iTRCH4       !< chemical tracer in ICON-ART
    INTEGER :: iTRCO2       !< chemical tracer in ICON-ART
    INTEGER :: iTRCO        !< chemical tracer in ICON-ART
    INTEGER :: iTRH2O       !< chemical tracer in ICON-ART
    INTEGER :: iTRO3        !< chemical tracer in ICON-ART
    INTEGER :: iTRSF6l      !< chemical tracer in ICON-ART
    INTEGER :: iTRSF6r      !< chemical tracer in ICON-ART
    INTEGER :: iTRSF6d      !< chemical tracer in ICON-ART
    INTEGER :: iTRN2O       !< chemical tracer in ICON-ART
    INTEGER :: iTR1         !< chemical tracer in ICON-ART
                            !< RR JS

    INTEGER :: nlev               !< number of full levels for each domain
    INTEGER :: nlevm1             !< number of half levels for each domain without boundaries
    INTEGER :: nlevp1             !< number of half levels for each domain with    boundaries

    LOGICAL :: lforcing           !< diabatic forcing TRUE/FALSE

    !> output mode (logicals)
    !  one or multiple of "none", "nml", "totint"
    TYPE t_output_mode
      LOGICAL :: l_none, l_nml, l_totint
    END TYPE t_output_mode

    TYPE (t_output_mode) output_mode

    !> file name for restart/checkpoint files (containg keyword
    !> substition patterns)
    CHARACTER(len=MAX_CHAR_LENGTH) :: restart_filename

    !> variable irad_type determines choice of radiation flux scheme
    !> irad_type=1: rrtm, irad_type=2: psrad
    INTEGER :: irad_type

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


