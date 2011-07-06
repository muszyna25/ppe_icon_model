!>
!! Data type containing basic control variables for a model integration. 
!!
!! Note that in a coupled simulation, each model component (e.g.,
!! atmosphere, ocean) will have its own run-configuration.
!! 
!! @author <name, affiliation>
!! @author <name, affiliation>
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
  USE mo_impl_constants, ONLY: max_dom


  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !>
  !! Derived type containing variables for time control. 
  !!
  TYPE :: t_run_config

    LOGICAL :: ldump_states    !< Compute interpolation coefficients and stop.
    LOGICAL :: lrestore_states !< Read interpolation coefficients from external file.

    LOGICAL :: ltestcase       !< Run idealized test case
    LOGICAL :: ldynamics       !< Switch on model dynamics
    LOGICAL :: ltheta_dyn      !< if .true., use potential temperature
                               !< times delta p as prognostic variable
    LOGICAL :: ltransport      !< Switch on tracer transport
    INTEGER :: ntracer         !< Total number of advected tracers
    INTEGER :: ntracer_static  !< Total number of non-advected tracers

    INTEGER :: iforcing        !< Choice of diabatic forcing
    INTEGER :: iequations      ! equation system:
    LOGICAL :: lcorio          ! if .TRUE.,  the Coriolis force is switched on,
                               ! if .FALSE., the Coriolis force is switched off
                               ! timer
    LOGICAL :: ltimer          ! if .TRUE.,  the timer is switched on
    INTEGER :: timers_level    ! what level of timers to run


    INTEGER   :: num_lev  ! number of full levels for each domain
    INTEGER   :: num_levp1! number of half levels for each domain
    INTEGER   :: nshift   ! half level of parent domain which coincides 
                                ! with the upper boundary of the current domain jg
    LOGICAL   :: lvert_nest !< switches on vertical nesting (.TRUE.)
    INTEGER   :: nvclev              ! no. of levels at which the coeffs A, B are given
    
    INTEGER   :: run_day              ! run length
    INTEGER   :: run_hour, run_minute ! - in day,hr,min,sec
    REAL(wp)  :: run_second

    INTEGER   :: nsteps              ! number of time steps
    REAL(wp)  :: dtime               ! [s] length of a time step
    REAL(wp)  :: dtrk(3)             ! [s] Runge Kutta 3 time steps [s]

    INTEGER   :: itopo     ! flag for topography handling
                         ! 0: corresponds to analytical topography,
                         ! 1: corresponds to netcdf files provided by
                         !    Herrmann Asensio

                           ! messages
    INTEGER  :: msg_level  ! Determines how much printout is generated during runtime

    INTEGER :: inextra_2d        !> number of extra output fields for debugging
    INTEGER :: inextra_3d        !> number of extra output fields for debugging
 
  END TYPE t_run_config 

  !>
  !! The actual variable
  !!
  TYPE(t_run_config) :: run_config(max_dom)

! CONTAINS


!  !---------------------------------------
!  !>
!  LOGICAL FUNCTION get_ldump_states()
!    get_ldump_states = run_config%ldump_states 
!  END FUNCTION get_ldump_states
!  !---------------------------------------
!  !>
!  LOGICAL FUNCTION get_lrestore_states()
!    get_lrestore_states = run_config%lrestore_states 
!  END FUNCTION get_lrestore_states
!  !---------------------------------------
!  !>
!  LOGICAL FUNCTION get_ltestcase()
!    get_ltestcase = run_config%ltestcase 
!  END FUNCTION get_ltestcase
!  !---------------------------------------
!  !>
!  LOGICAL FUNCTION get_ldynamics()
!    get_ldynamics = run_config%ldynamics 
!  END FUNCTION get_ldynamics
!  !---------------------------------------
!  !>
!  LOGICAL FUNCTION get_ltransport()
!    get_ltransport = run_config%ltransport 
!  END FUNCTION get_ltransport
!  !---------------------------------------
!  !>
!  INTEGER FUNCTION get_ntracer()
!    get_ntracer = run_config%ntracer 
!  END FUNCTION get_ntracer
!  !---------------------------------------
!  !>
!  INTEGER FUNCTION get_ntracer_static()
!    get_ntracer_static = run_config%ntracer_static
!  END FUNCTION get_ntracer_static
!  !---------------------------------------
!  !>
!  INTEGER FUNCTION get_iforcing()
!    get_iforcing = run_config%iforcing 
!  END FUNCTION get_iforcing
!  !---------------------------------------
!  !>
!  REAL(wp) FUNCTION get_dtime()
!    get_dtime = run_config%dtime 
!  END FUNCTION get_dtime
!  !---------------------------------------
!  !>
!  INTEGER FUNCTION get_nsteps()
!    get_nsteps = run_config%nsteps 
!  END FUNCTION get_nsteps
!  !---------------------------------------
!

END MODULE mo_run_config

