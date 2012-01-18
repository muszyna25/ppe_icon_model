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
                               IMPIOM, INOFORCING, ILDF_DRY

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ldump_states, lrestore_states, ltestcase, ldynamics, iforcing, lforcing
  PUBLIC :: ltransport, ntracer, ntracer_static, nlev, nlevp1, nvclev
  PUBLIC :: lvert_nest, num_lev, num_levp1, nshift, nsteps, dtime, dtime_adv
  PUBLIC :: ltimer, timers_level, activate_sync_timers, msg_level
  PUBLIC :: iqv, iqc, iqi, iqs, iqr, iqhydro, nqtendphy, iqt, io3, ico2
  PUBLIC :: check_epsilon, testbed_mode
  PUBLIC :: configure_run, l_one_file_per_patch, ldump_dd, lread_dd, nproc_dd

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

    ! Namelist variables
    !
    LOGICAL :: ldump_states    !< Compute interpolation coefficients and stop.
    LOGICAL :: lrestore_states !< Read interpolation coefficients from external file.
    LOGICAL :: l_one_file_per_patch !< Use one or several files per patch for dump/restore
    LOGICAL :: ldump_dd        !< Compute domain decomposition and stop
    LOGICAL :: lread_dd        !< Read domain decomposition from external file
    INTEGER :: nproc_dd        !< Number of procs for domain decomposition (when done on 1 CPU)

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
  
    REAL(wp):: check_epsilon   !< small value for checks
    INTEGER :: testbed_mode

    INTEGER :: msg_level       !< how much printout is generated during runtime


    ! Derived variables
    !
    ! Tracer indices of water species
    INTEGER :: iqv        !< water vapor
    INTEGER :: iqc        !< cloud water
    INTEGER :: iqi        !< cloud ice
    INTEGER :: iqr        !< rain water
    INTEGER :: iqs        !< snow
    INTEGER :: iqhydro    !< index of last hydrometeor to ease summation over all of them
    INTEGER :: nqtendphy  !< number of water species for which physical tendencies are stored
  
    ! Tracer indices of other species
    INTEGER :: io3        !< O3
    INTEGER :: ico2       !< CO2
    INTEGER :: iqt        !< start index of other tracers than hydrometeors


    REAL(wp) :: dtime_adv !< advective timestep on global patch (iadv_rcf*dtime) [s]

    INTEGER :: num_levp1(MAX_DOM) !< number of half levels for each domain
    INTEGER :: nlev               !< number of full levels for each domain
    INTEGER :: nlevp1             !< number of half levels for each domain
    INTEGER :: nvclev             !< number of levels at which the coeffs A, B are given

    INTEGER :: ntracer_static     !< Total number of non-advected tracers

    LOGICAL :: lforcing           !< diabatic forcing TRUE/FALSE

CONTAINS
  !>
  !!
  !! Assign value to components of the run configuration state that have no
  !! corresponding namelist variable. 
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


!  !---------------------------------------
!  !>
!  LOGICAL FUNCTION get_ldump_states()
!    get_ldump_states = ldump_states 
!  END FUNCTION get_ldump_states
!  !---------------------------------------
!  !>
!  LOGICAL FUNCTION get_lrestore_states()
!    get_lrestore_states = lrestore_states 
!  END FUNCTION get_lrestore_states
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
!  INTEGER FUNCTION get_ntracer_static()
!    get_ntracer_static = ntracer_static
!  END FUNCTION get_ntracer_static
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

