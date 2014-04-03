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
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_run_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: MAX_DOM, IHELDSUAREZ, INWP, IECHAM, ILDF_ECHAM, &
                               IMPIOM, INOFORCING, ILDF_DRY, MAX_CHAR_LENGTH,  &
                               TIMER_MODE_AGGREGATED, TIMER_MODE_DETAILED
  USE mo_io_units,       ONLY: filename_max
  USE mo_grid_config,    ONLY: get_grid_rescale_factor

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ltestcase, ldynamics, iforcing, lforcing
  PUBLIC :: ltransport, ntracer, nlev, nlevm1, nlevp1, nvclev
  PUBLIC :: lvert_nest, num_lev, num_levp1, nshift, nsteps, dtime, dtime_adv
  PUBLIC :: ltimer, timers_level, activate_sync_timers, msg_level
  PUBLIC :: iqv, iqc, iqi, iqs, iqr, iqtvar, nqtendphy, iqt, ico2
  PUBLIC :: iqni, iqni_nuc, iqg, iqm_max
  PUBLIC :: iqh, iqnh, iqnr, iqns, iqng
  PUBLIC :: iash1,iash2,iash3,iash4,iash5,iash6 !Running index for Volcanic Ash in ICON-ART 
  PUBLIC :: iCS137,iI131,iTE132,iZR95,iXE133,iI131g,iI131o,iBA140,iRU103 !Running index for radioactive nuclides  in ICON-ART
  PUBLIC :: iseasa,iseasb,iseasc,iseasa0,iseasb0,iseasc0 !Running index for sea salt in ICON-ART
  PUBLIC :: idusta,idustb,idustc,idusta0,idustb0,idustc0 !Running index for mineral dust in ICON-ART
  PUBLIC :: iTRCHBR3,iTRCH2BR2,iTRBRy             !Running index for chemical tracer in ICON-ART - VSLS
  PUBLIC :: iTRCH4,iTRCO2,iTRCO,iTRH2O,iTRO3      !Running index for chemical tracer in ICON-ART - CH4-CO-CO2-H20-O3
  PUBLIC :: iTRSF6l,iTRSF6r,iTRSF6d               !Running index for chemical tracer in ICON-ART - SF6
  PUBLIC :: iTRN2O,iTR1                           !Running index for chemical tracer in ICON-ART - other
  PUBLIC :: grid_generatingCenter     ! non-namelist variables
  PUBLIC :: grid_generatingSubcenter  ! non-namelist variables
  PUBLIC :: number_of_grid_used       ! non-namelist variables
  PUBLIC :: check_epsilon, test_mode
  PUBLIC :: configure_run
  PUBLIC :: output, t_output_mode, output_mode, max_output_modes
  PUBLIC :: debug_check_level
  PUBLIC :: restart_filename
  PUBLIC :: profiling_output, TIMER_MODE_AGGREGATED, TIMER_MODE_DETAILED


  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

    ! Namelist variables
    !
    LOGICAL :: ltestcase       !< Run idealized test case
    LOGICAL :: ldynamics       !< Switch on model dynamics
    INTEGER :: iforcing        !< Choice of diabatic forcing

    LOGICAL :: ltransport      !< Switch on tracer transport
    INTEGER :: ntracer         !< Total number of advected tracers

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
  
    REAL(wp):: check_epsilon   !< small value for checks
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
    INTEGER :: idusta       !< Sea Salt Aerosol Mode A Mass Density
    INTEGER :: idustb       !< Sea Salt Aerosol Mode B Mass Density
    INTEGER :: idustc       !< Sea Salt Aerosol Mode C Mass Density
    INTEGER :: idusta0      !< Sea Salt Aerosol Mode A Number Density
    INTEGER :: idustb0      !< Sea Salt Aerosol Mode B Number Density
    INTEGER :: idustc0      !< Sea Salt Aerosol Mode C Number Density
    INTEGER :: iTRCHBR3        !< chemical tracer in ICON-ART
    INTEGER :: iTRCH2BR2       !< chemical tracer in ICON-ART
    INTEGER :: iTRBRy          !< chemical tracer in ICON-ART
    INTEGER :: iTRCH4          !< chemical tracer in ICON-ART
    INTEGER :: iTRCO2          !< chemical tracer in ICON-ART
    INTEGER :: iTRCO           !< chemical tracer in ICON-ART
    INTEGER :: iTRH2O          !< chemical tracer in ICON-ART
    INTEGER :: iTRO3           !< chemical tracer in ICON-ART
    INTEGER :: iTRSF6l         !< chemical tracer in ICON-ART
    INTEGER :: iTRSF6r         !< chemical tracer in ICON-ART
    INTEGER :: iTRSF6d         !< chemical tracer in ICON-ART
    INTEGER :: iTRN2O          !< chemical tracer in ICON-ART
    INTEGER :: iTR1            !< chemical tracer in ICON-ART

    REAL(wp) :: dtime_adv = 0.0_wp!< advective timestep on global patch (iadv_rcf*dtime) [s]

    INTEGER :: num_levp1(MAX_DOM) !< number of half levels for each domain
    INTEGER :: nlev               !< number of full levels for each domain
    INTEGER :: nlevm1             !< number of half levels for each domain without boundaries
    INTEGER :: nlevp1             !< number of half levels for each domain with    boundaries
    INTEGER :: nvclev             !< number of levels at which the coeffs A, B are given

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

CONTAINS
  !>
  !!
  !! Assign value to components of the run configuration state that have no
  !! corresponding namelist variable.
  !!
  !! Exceptions: grid_generatingCenter, grid_generatingSubcenter and number_of_grid_used 
  !!             are set in mo_model_domimp_patches/read_basic_patch 
  !!
  SUBROUTINE configure_run( opt_iadv_rcf )

    INTEGER, OPTIONAL, INTENT(IN) :: opt_iadv_rcf  !< reduced calling freq. for advection

    CHARACTER(LEN=*),PARAMETER :: routine = 'mo_run_config:configure_run'
    
    !----------------------------
    ! advective timestep on global patch
    IF ( PRESENT(opt_iadv_rcf) ) THEN
      dtime_adv = REAL(opt_iadv_rcf,wp) * dtime  ! NH-atm
    ELSE
      dtime_adv = dtime  ! oce, H-atm
      WRITE(0,*) routine, ': dtime_adv initialized with', dtime
    ENDIF


    !----------------------------
    ! Number of vertical levels

    IF (.NOT.lvert_nest) THEN
      num_lev(:) = nlev
      nshift (:) = 0
    END IF

    nlevm1       = nlev - 1
    nlevp1       = nlev + 1
    nvclev       = nlev + 1
    num_levp1(:) = num_lev(:) + 1

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

