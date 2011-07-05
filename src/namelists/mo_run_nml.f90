!>
!!        
!! @par Revision History
!!   Revision History in mo_global_variables.f90 (r3513)
!!   Modification by Constantin Junk (2010-02-22)
!!     - separated namelist run_ctl from mo_global_variables
!!     - therefore, added mo_run_nml
!!     - restructured variable declaration section
!!     - introduced initialization variables 
!!       (iinit,ianalytic,irestart)
!!     - added cell geomtry variables itri and ihex
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
!!
MODULE mo_run_nml

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: max_char_length, max_dom, max_ntracer,     &
    &                              ianalytic, ihs_atm_temp, itri, ihex,       &
    &                              inoforcing, ildf_echam, ishallow_water,    &
    &                              iecham, ihs_atm_theta, inh_atmosphere,     &
    &                              ihs_ocean, iheldsuarez, inwp, impiom,      &
    &                              ildf_dry
  USE mo_physical_constants, ONLY: grav
  USE mo_datetime,           ONLY: t_datetime, proleptic_gregorian,          &
                                 & date_to_time, add_time, print_datetime_all
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_master_nml,         ONLY: lrestart
  USE mo_io_restart_attributes, ONLY: get_restart_attribute
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC
  PRIVATE :: run_ctl, calendar,                                              &
    &        ini_year, ini_month, ini_day, ini_hour, ini_minute, ini_second, &
    &        end_year, end_month, end_day, end_hour, end_minute, end_second, &
    &                             run_day, run_hour, run_minute, run_second


  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters
  ! ------------------------------------------------------------------------

  ! initialization
  ! --------------
  INTEGER          :: iinit               ! model initialization:


  ! computing setup
  ! ---------------

  ! number of levels
  ! ----------------
  INTEGER          :: nlev                ! number of full levels = number of layers
  INTEGER          :: nlevp1              ! number of half levels = nlev+1
  INTEGER          :: num_lev(max_dom)  ! number of full levels for each domain
  INTEGER          :: num_levp1(max_dom)! number of half levels for each domain
  INTEGER          :: nshift(max_dom)   ! half level of parent domain which coincides 
                                ! with the upper boundary of the current domain jg
  LOGICAL          :: lvert_nest !< switches on vertical nesting (.TRUE.)
  INTEGER          :: nvclev              ! no. of levels at which the coeffs A, B are given
  INTEGER          :: ntracer             ! number of advected tracers
  INTEGER          :: ntracer_static      ! number of non-advected tracers

  ! time information
  ! ----------------
  !
  ! calendar type
  INTEGER          :: calendar
  !
  ! initial date and time
  ! - namelist variables
  INTEGER          :: ini_year, ini_month, ini_day
  INTEGER          :: ini_hour, ini_minute
  REAL(wp)         :: ini_second
  ! - data and time structure
  TYPE(t_datetime) :: ini_datetime
  !
  ! current model time, not a namelist variable
  TYPE(t_datetime) :: current_datetime
  !
  ! end date and time
  ! - namelist variables
  INTEGER          :: end_year, end_month, end_day
  INTEGER          :: end_hour, end_minute
  REAL(wp)         :: end_second
  ! - data and time structure
  TYPE(t_datetime) :: end_datetime
  !
  ! run length
  ! - in day,hr,min,sec
  INTEGER          :: run_day
  INTEGER          :: run_hour, run_minute
  REAL(wp)         :: run_second
  ! - in time steps
  INTEGER          :: nsteps              ! number of time steps
  REAL(wp)         :: dtime               ! [s] length of a time step
  REAL(wp)         :: dtrk(3)    ! [s] Runge Kutta 3 time steps [s]

  ! restart interval
  ! ----------------
  REAL(wp)         :: dt_restart          ! [s] length of a restart cycle 

  ! equations to be solved
  ! ----------------------
  INTEGER          :: iequations          ! equation system:


  ! cell geometry
  ! -------------
  INTEGER          :: i_cell_type         ! cell type:


  ! dynamics
  ! --------
  LOGICAL  :: ldynamics  ! if .TRUE.=default then
                         ! - read &dynamics_ctl namelist and
                         ! - compute adiabatic dynamical tendencies
                         ! else skip dynamics.

  ! parameterized forcing (right hand side) of dynamics
  INTEGER          :: iforcing            ! forcing package

  
  ! auxiliary variables to access single fields of the 4D tracer array
  ! - H2O
  INTEGER          :: iqv        !> water vapor
  INTEGER          :: iqc        !! cloud water
  INTEGER          :: iqi        !! cloud ice
  INTEGER          :: iqr        !! rain water
  INTEGER          :: iqs        !! snow
  INTEGER          :: iqcond     !! index of last hydrometeor to ease summation over all of them
  ! - other species
  INTEGER          :: io3        !< O3
  INTEGER          :: ico2       !< CO2

  INTEGER          :: iqt        !< start index of other tracers than hydrometeors

  ! transport
  ! ---------
  LOGICAL  :: ltransport ! if .TRUE. then
                         ! - read &transport_ctl namelist and
                         ! - compute tracer tendencies by transport
                         ! else skip

  ! testcases
  !----------
  LOGICAL  :: ltestcase  ! if .TRUE. then
                         ! - read &testcase_ctl namelist,
                         ! - compute analytical initial state,
                         !   depending on the specified test case,
                         ! - compute analytical boundary conditions,
                         ! - and compute analytical forcing
                         ! else
                         ! - initialize model from external data
                         ! - compute boundary conditions from external data
                         ! - compute processes as specified elsewhere

  LOGICAL  :: lcorio     ! if .TRUE.,  the Coriolis force is switched on,
                         ! if .FALSE., the Coriolis force is switched off

  INTEGER  :: itopo      ! flag for topography handling
                         ! 0: corresponds to analytical topography,
                         ! 1: corresponds to netcdf files provided by
                         !    Herrmann Asensio

  ! messages
  INTEGER  :: msg_level  ! Determines how much printout is generated during runtime

  INTEGER :: inextra_2d        !> number of extra output fields for debugging
  INTEGER :: inextra_3d        !> number of extra output fields for debugging


  ! timer
  ! -----
  LOGICAL :: ltimer     ! if .TRUE.,  the timer is switched on
  INTEGER :: timers_level = 1  ! what level of timers to run

  ! dump/restore
  ! ------------
  LOGICAL :: ldump_states    ! Dump patch/interpolation/grid refinement state of every
                             ! patch (after subdivision in case of a parallel run)
                             ! end exit program
  LOGICAL :: lrestore_states ! Restore patch/interpolation/grid refinement states
                             ! from dump files instead of calculating them
  !
  LOGICAL :: ltheta_dyn ! if .true., use potential temperature times delta p as prognostic variable

  NAMELIST /run_ctl/ iinit,                                &
    &                nlev, num_lev, nshift,                &
    &                lvert_nest, ntracer, calendar,        &
    &                ini_year, ini_month,  ini_day,        &
    &                ini_hour, ini_minute, ini_second,     &
    &                end_year, end_month,  end_day,        &
    &                end_hour, end_minute, end_second,     &
    &                run_day,                              &
    &                run_hour, run_minute, run_second,     &
    &                nsteps, dt_restart, dtime,            &
    &                iequations, i_cell_type,              &
    &                ldynamics, ltransport, iforcing,      &
    &                ltestcase, lcorio, itopo, msg_level,  &
    &                ldump_states, lrestore_states, ltimer,&
    &                timers_level, inextra_2d, inextra_3d

!--------------------------------------------------------------------
! for external parameters
!--------------------------------------------------------------------

  ! Namelist variables
  !
  REAL(wp):: fac_smooth_topo
  INTEGER :: n_iter_smooth_topo

  NAMELIST/ext_par_ctl/ fac_smooth_topo,n_iter_smooth_topo


  !
  ! -----------------------------------------------------------------------
  ! 2.0 Declaration of dependent control variables 
  ! -----------------------------------------------------------------------
  !
  LOGICAL  :: lshallow_water! if .TRUE., the model runs in shallow water mode
  LOGICAL  :: latmosphere   ! if .TRUE., the model runs in atmosphere mode
  LOGICAL  :: locean        ! if .TRUE., the model runs in ocean mode
  LOGICAL  :: lhydrostatic  ! if .TRUE., the model runs in hydrostatic mode
  LOGICAL  :: lforcing      ! if .TRUE., the model runs with parameterized forcing
  !


CONTAINS
  !-------------------------------------------------------------------------
  !>
  !!  Initialization of variables that contain general information.
  !!
  !!  Initialization of variables that contain general information
  !!  about the model run. The configuration is read from
  !!  namelist 'run_ctl'.
  !!
  !! @par Revision History
  !!  Reading of the namelist and checking of the validity were
  !!  in some other modules in the earlier version of the shallow water model.
  !!  Moved to this module by Hui Wan, MPI-M (2007-02-23)
  !!  The character parameter <i>routine</i> was introduced and used
  !!  for error information by Hui Wan, MPI-M (2007-02-23).
  !!  Modified by Almut Gassmann, MPI-M (2008-09-23)
  !!  - introduced i_cell_type, lidealized and lshallow_water
  !!  Modified by Marco Giorgetta, MPI-M (2009-02-23)
  !!  - lidealized replaced by ltestcase
  !!  Modification by Constantin Junk, MPI-M (2010-02-22)
  !!  - changes to consistency checks
  !!
  SUBROUTINE run_nml_setup
                                               
   INTEGER  :: istat, funit, calendar_old
   INTEGER  :: ini_year_old, ini_month_old, ini_day_old, ini_hour_old, ini_minute_old
   INTEGER  :: restart_year, restart_month, restart_day, restart_hour, restart_minute
   REAL(wp) :: ini_second_old
   REAL(wp) :: restart_second
   REAL(wp) :: cur_datetime_calsec, end_datetime_calsec, length_sec

   CHARACTER(len=max_char_length), PARAMETER ::   &
            &  routine = 'mo_run_nml/run_nml_setup'

   !------------------------------------------------------------
   ! 1. set up the default values for run_ctl
   !------------------------------------------------------------
   ! initialization
   iinit          = ianalytic

   ! dimensions for new, initialized experiments
   nlev           = 31  ! number of full levels
   num_lev(:)     = nlev ! number of full levels for each domain
   nshift(:)      = 0    ! please do not change the default.
                         ! otherwise the initialization of 
                         ! p_patch(jg)%nshift in "import patches" 
                         ! will not work properly.
   lvert_nest = .FALSE.  ! no vertical nesting
   ntracer        = 0   ! number of advected tracers
   ntracer_static = 0   ! number of non-advected tracers

   ! initial date and time
   calendar       = proleptic_gregorian
   ini_year       = 2008
   ini_month      = 9
   ini_day        = 1
   ini_hour       = 0
   ini_minute     = 0
   ini_second     = 0.0_wp
   !
   ! end date and time
   end_year       = 2008
   end_month      = 9
   end_day        = 1
   end_hour       = 1
   end_minute     = 40
   end_second     = 0.0_wp
   !
   ! length of integration = (number of timesteps)*(length of timestep)
   ! - If nsteps is set to a non-zero positive value, then the end date is computed
   !   from the initial date and time, the time step dtime, and nsteps.
   ! - Else if run_day, run_hour, run_minute or run_second is set to a non-zero,
   !   positive value, then the initial date and time and the run_... variables are
   !   used to compute the end date and time and, using dtime, nsteps.
   !   Else nsteps is computed from the initial and end date and time and dtime.
   !
   ! initialize run_... variables with zero
   run_day        = 0
   run_hour       = 0
   run_minute     = 0
   run_second     = 0.0_wp
   !
   ! initialize nsteps with zero
   nsteps         = 0
   !
   ! length of restart cycle
   dt_restart     = 86400._wp*30._wp   ! = 30 days
   !
   ! time step
   dtime          = 600._wp   ! [s] for R2B04 + semi-implicit time steppping

   ! select model and numerics
   iequations     = ihs_atm_temp
   i_cell_type    = itri

   ! switches for tendency computation
   ldynamics      = .TRUE.
   ltransport     = .FALSE.
   iforcing       = inoforcing

   ! switch for running a predefined testcase
   ! details of the testcase are controled by 'testcase_ctl'
   ltestcase      = .TRUE.

   lcorio         = .TRUE.
   itopo          = 0
   msg_level      = 10
   ltimer         = .TRUE.

   ! dump/restore
   ldump_states    = .FALSE.
   lrestore_states = .FALSE.

   ! The following values are deduced from the namelist variables:
   ! auxiliary switches set as function of iequations
   ! (not included in run_ctl)
   lshallow_water = .FALSE.
   latmosphere    = .FALSE.
   locean         = .FALSE.
   lhydrostatic   = .FALSE.
   lforcing       = .FALSE.

   ! For debugging purposes define number of output variables
   inextra_2d      = 0           !> no extra output 2D fields
   inextra_3d      = 0           !> no extra output 3D fields

    !------------------------------------------------------------------------                  
    ! 2. If this is a resumed integration...
    !------------------------------------------------------------------------                  
    IF (lrestart) THEN                                                                 

      ! 2.1 Overwrite the defaults above by values in the restart file

      funit = open_and_restore_namelist('run_ctl')
      READ(funit,NML=run_ctl)
      CALL close_tmpfile(funit)

      ! 2.2 Save the calendar and initial date/time of the old run

      calendar_old    = calendar
      ini_year_old    = ini_year
      ini_month_old   = ini_month
      ini_day_old     = ini_day
      ini_hour_old    = ini_hour
      ini_minute_old  = ini_minute
      ini_second_old  = ini_second

      ! 2.2 Inquire the date/time at which the previous run stopped

      CALL get_restart_attribute( 'current_year'  , restart_year   )
      CALL get_restart_attribute( 'current_month' , restart_month  )
      CALL get_restart_attribute( 'current_day'   , restart_day    )
      CALL get_restart_attribute( 'current_hour'  , restart_hour   )
      CALL get_restart_attribute( 'current_minute', restart_minute )
      CALL get_restart_attribute( 'current_second', restart_second )

    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL position_nml('run_ctl', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, run_ctl)
    END SELECT

    !---------------------------------------------------------------
    ! 4. Check whether the namelist varibles have reasonable values
    !---------------------------------------------------------------

    IF (nlev < 1)  CALL finish(TRIM(routine),'"nlev" must be positive')
    nlevp1 = nlev+1
    nvclev = nlevp1

    ! vertical nesting
    IF (.NOT. lvert_nest) THEN
      ! overwrite num_lev with nlev
      num_lev(1:max_dom) = nlev
      ! set nshift to 0
      nshift(1:max_dom) = 0 
    ENDIF
    num_levp1(1:max_dom) = num_lev(1:max_dom) + 1

    IF (ANY(num_lev < 0)) CALL finish(TRIM(routine),'"num_lev" must be positive')
    IF (ANY(nshift < 0)) CALL finish(TRIM(routine),'"nshift" must be positive')

    SELECT CASE (i_cell_type)
    CASE (itri,ihex)
      ! ok
    CASE default
      CALL finish( TRIM(routine),'wrong cell type specifier, "i_cell_type" must be 3 or 6')
    END SELECT

    SELECT CASE (iequations)
    CASE (ishallow_water)
      lshallow_water = .TRUE.
      IF ( nlev/=1 ) THEN
        CALL finish(TRIM(routine),'Shallow water model needs nlev=1')
      ENDIF
      ltheta_dyn     = .FALSE.
    CASE (ihs_atm_temp)
      latmosphere    = .TRUE.
      lhydrostatic   = .TRUE.
      ltheta_dyn     = .FALSE.
    CASE (ihs_atm_theta)
      latmosphere    = .TRUE.
      lhydrostatic   = .TRUE.
      ltheta_dyn     = .TRUE.
    CASE (inh_atmosphere)
      latmosphere    = .TRUE.
      lhydrostatic   = .FALSE.
    CASE (ihs_ocean)
      locean         = .TRUE.
      lhydrostatic   = .TRUE.
    CASE default
      CALL finish( TRIM(routine),'wrong equation specifier iequations')
    END SELECT

    SELECT CASE (iforcing)
    CASE (iheldsuarez, inwp, iecham, ildf_echam, impiom)
     lforcing = .TRUE.
    CASE (inoforcing, ildf_dry)
     lforcing = .FALSE.
    CASE DEFAULT
      CALL finish( TRIM(routine),'wrong forcing specifier iforcing')
    END SELECT

    IF(ntracer<0 .OR. ntracer>max_ntracer) THEN
      CALL finish( TRIM(routine),'wrong number of tracers. Valid range: 0<= ntracer <=20')
    ENDIF

    SELECT CASE(iforcing)
    CASE (iecham,ildf_echam)
      IF (ntracer<3) CALL finish( TRIM(routine),'ECHAM forcing needs at least 3 tracers')
      iqv    = 1     !> water vapour
      iqc    = 2     !! cloud water
      iqi    = 3     !! ice
      iqcond = iqi   !! index of last hydrometeor to ease summation over all of them
      iqt    = 4     !! starting index of non-water species 
      io3    = 5     !! O3
      ico2   = 6     !! CO2
    CASE (inwp)
      iqv    = 1     !> water vapour
      iqc    = 2     !! cloud water
      iqi    = 3     !! ice
      iqr    = 4     !! rain water
      iqs    = 5     !! snow
      iqcond = iqs   !! index of last hydrometeor to ease summation over all of them
      io3    = 6     !! O3
      ico2   = 7     !! CO2
      iqt    = 6     !! start index of other tracers than hydrometeors
    END SELECT

    ! time step
    IF (dtime  <= 0._wp) CALL finish(routine,'"dtime" must be positive')

    !---------------------------------------------------------------
    ! 5. Set up model time
    !---------------------------------------------------------------
    ! 5.1 Initial date and time

    ini_datetime%calendar = calendar
    ini_datetime%year     = ini_year
    ini_datetime%month    = ini_month
    ini_datetime%day      = ini_day
    ini_datetime%hour     = ini_hour
    ini_datetime%minute   = ini_minute
    ini_datetime%second   = ini_second

    CALL date_to_time(ini_datetime) ! fill date time structure
    CALL message(' ',' ')
    CALL message(routine,'Initial date and time')
    CALL message(routine,'---------------------')
    CALL print_datetime_all(ini_datetime)  ! print all date and time components

    ! 5.2 Current date and time:

    IF (lrestart) THEN
      ! In a resumed integration, if the calendar or initial date/time 
      ! is different from those in the restart file,
      ! we regard this integration as a new one with its own calendar. 
      ! Model time at which the previous run stopped is thus not relevant. 
      ! Simulation will start from the user-specified initial date/time,
      ! which is also the current model date/time.

      IF (calendar  /=calendar_old   .OR.                                 &
          ini_year  /=ini_year_old   .OR. ini_month  /=ini_month_old .OR. &
          ini_day   /=ini_day_old    .OR. ini_hour   /=ini_hour_old  .OR. &
          ini_minute/=ini_minute_old .oR. ini_second /=ini_second_old     ) THEN

        current_datetime = ini_datetime

      ELSE
      ! Otherwise we start from the point when the previous integration stopped.

        current_datetime%calendar = calendar
        current_datetime%year     = restart_year
        current_datetime%month    = restart_month
        current_datetime%day      = restart_day
        current_datetime%hour     = restart_hour
        current_datetime%minute   = restart_minute
        current_datetime%second   = restart_second

        CALL date_to_time(current_datetime) ! fill date time structure
      END IF

    ELSE
      ! In an initial run, current date/time is, naturally, the initial date/time
      current_datetime = ini_datetime

    END IF !lrestart

    CALL message(' ',' ')
    CALL message(' ',' ')
    CALL message(routine,'Current date and time')
    CALL message(routine,'---------------------')
    CALL print_datetime_all(current_datetime)  ! print all date and time components

    ! 5.3 End date and time, and length of integration
    !     Here we define "nsteps" as the number of time steps THIS integration
    !     will last, regardless of the restart status.

    IF (nsteps/=0) THEN

      IF (nsteps < 0) CALL finish(routine,'"nsteps" must not be negative')
      length_sec   = REAL(nsteps,wp)*dtime
      end_datetime = current_datetime
      CALL add_time(length_sec,0,0,0,end_datetime)
      !
    ELSE IF (run_day/=0 .OR. run_hour/=0 .OR. run_minute/=0 .OR. run_second/=0.0_wp) THEN
      IF (run_day    < 0    ) CALL finish(routine,'"run_day" must not be negative')
      IF (run_hour   < 0    ) CALL finish(routine,'"run_hour" must not be negative')
      IF (run_minute < 0    ) CALL finish(routine,'"run_minute" must not be negative')
      IF (run_second < 0._wp) CALL finish(routine,'"run_second" must not be negative')
      !
      end_datetime = current_datetime
      CALL add_time(run_second,run_minute,run_hour,run_day,end_datetime)
      !
      cur_datetime_calsec = (REAL(current_datetime%calday,wp)+current_datetime%caltime) &
        &                   *REAL(current_datetime%daylen,wp)
      end_datetime_calsec = (REAL(end_datetime%calday,wp)+end_datetime%caltime) &
        &                   *REAL(end_datetime%daylen,wp)
      nsteps=INT((end_datetime_calsec-cur_datetime_calsec)/dtime)
      !
    ELSE
      ! compute nsteps from current_datetime, end_datetime and dtime
      end_datetime%calendar = calendar
      end_datetime%year     = end_year
      end_datetime%month    = end_month
      end_datetime%day      = end_day
      end_datetime%hour     = end_hour
      end_datetime%minute   = end_minute
      end_datetime%second   = end_second
      CALL date_to_time      (end_datetime) ! fill date time structure
      !
      cur_datetime_calsec = (REAL(current_datetime%calday,wp)+current_datetime%caltime) &
        &                   *REAL(current_datetime%daylen,wp)
      end_datetime_calsec = (REAL(end_datetime%calday,wp)+end_datetime%caltime) &
        &                   *REAL(end_datetime%daylen,wp)
      IF (end_datetime_calsec < cur_datetime_calsec) &
        & CALL finish(routine,'The end date and time must not be '// &
        &            'before the current date and time')
      !
      nsteps=INT((end_datetime_calsec-cur_datetime_calsec)/dtime)
      !
    END IF

    nsteps = MIN(nsteps,INT(dt_restart/dtime))

    IF (lhydrostatic) THEN
     ! If running the HYDROSTATIC version,
     ! let the model integrate one more step after the desired end of
     ! simulation in order to get the proper output. This additional step is
     ! necessary because the HYDROSTATIC model writes out values of step N
     ! after the integration from N to N+1 is finished. Also note that
     ! this additional step is done only for the regular output, and is 
     ! ignored for restart.

     nsteps = nsteps + 1

     ! The additional step is not needed in the NON-hydrostatic version because
     ! in this case the model writes out values of step N
     ! after the integration from N-1 to N is finished.
    ENDIF
    !..................................................................
    !
    CALL message(' ',' ')
    CALL message(routine,'End date and time')
    CALL message(routine,'-----------------')
    CALL print_datetime_all(end_datetime)  ! print all date and time components

    CALL message(' ',' ')
    CALL message(routine,'Length of restart cycle')
    CALL message(routine,'-----------------------')
    WRITE(message_text,'(a,f10.2,a,f16.10,a)') &
         &'dt_restart :',dt_restart,' seconds =', dt_restart/86400._wp, ' days'
    CALL message(routine,message_text)

    CALL message(' ',' ')
    CALL message(routine,'Length of this run')
    CALL message(routine,'------------------')
    WRITE (message_text,'(a,f7.2)') 'dtime [s] :',dtime
    CALL message(routine,message_text)
    WRITE (message_text,'(a,i7)')   'nsteps    :',nsteps
    CALL message(routine,message_text)
    CALL message(' ',' ')


    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=run_ctl)
    CALL store_and_close_namelist(funit, 'run_ctl')

    ! write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=run_ctl)

    !-----------------------------------------------------------
    ! Topography
    !-----------------------------------------------------------
    SELECT CASE (itopo)
    CASE (0)
      ! ok
    CASE (1)

      ! set up the default values
      fac_smooth_topo    = 0.015625_wp
      n_iter_smooth_topo = 35

      ! If this is a resumed integration, overwrite the defaults above 
      ! by values in the previous integration.
      IF (lrestart) THEN
        funit = open_and_restore_namelist('ext_par_ctl')
        READ(funit,NML=ext_par_ctl)
        CALL close_tmpfile(funit)
      END IF

      ! read namelist for external parameters
      CALL position_nml ('ext_par_ctl', status=istat)
      SELECT CASE (istat)
      CASE (POSITIONED)
        READ (nnml, ext_par_ctl)
      END SELECT

      ! Store the namelist for restart
      funit = open_tmpfile()
      WRITE(funit,NML=ext_par_ctl)
      CALL store_and_close_namelist(funit, 'ext_par_ctl')

      !write the contents of the namelist to an ASCII file
      IF(p_pe == p_io) WRITE(nnml_output,nml=ext_par_ctl)

    CASE default
       CALL finish( TRIM(routine),'wrong topography specifier, itopo must be in {0,1}]')
    END SELECT

 END SUBROUTINE run_nml_setup

END MODULE mo_run_nml
