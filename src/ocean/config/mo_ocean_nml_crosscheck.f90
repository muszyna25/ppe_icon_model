!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ocean_nml_crosscheck

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, message_text, finish, print_value, warning
  USE mo_grid_config,       ONLY: init_grid_configuration
  USE mo_parallel_config,   ONLY: check_parallel_configuration, p_test_run, l_fast_sum
  USE mo_run_config,        ONLY: nsteps, dtime, nlev
  USE mo_time_config,       ONLY: time_config, restart_experiment
  USE mo_datetime,          ONLY: add_time, print_datetime_all
  USE mo_io_config,         ONLY: dt_checkpoint, write_initial_state
  USE mo_grid_config,       ONLY: grid_rescale_factor, use_duplicated_connectivity
  USE mo_ocean_nml
  USE mo_master_config,     ONLY: isRestart

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ocean_crosscheck

CONTAINS

  !>
  !! Check and, if necessary, adapt simulation lengths
  !! This partially duplicates the code resize_atmo_simulation_length
  !!  not a nice way around it though, needs to be cleaned
  !!
  !! @par Revision History
  !! Initial revision by Hui Wan, MPI (2011-07)
  !!
!<Optimize:inUse>
  SUBROUTINE resize_ocean_simulation_length()

    REAL(wp):: cur_datetime_calsec, end_datetime_calsec, length_sec
    CHARACTER(len=*), PARAMETER :: method_name =  'mo_ocean_nml_crosscheck:resize_ocean_simulation_length'

    !----------------------------
    ! rescale timestep
    dtime     = dtime     * grid_rescale_factor
    !---------------------------------
    ! Check length of this integration
    !---------------------------------
    !
    IF (nsteps>=0) THEN   ! User specified a value

      length_sec = REAL(nsteps,wp)*dtime
      time_config%end_datetime = time_config%cur_datetime
      CALL add_time(length_sec,0,0,0,time_config%end_datetime)

      !HW (2011-07-17): run_day/hour/... not implemented in the restructured version ------
      !ELSE IF (run_day/=0 .OR. run_hour/=0 .OR. run_minute/=0 .OR. run_second/=0.0_wp) THEN
      !  IF (run_day    < 0    ) CALL finish(method_name,'"run_day" must not be negative')
      !  IF (run_hour   < 0    ) CALL finish(method_name,'"run_hour" must not be negative')
      !  IF (run_minute < 0    ) CALL finish(method_name,'"run_minute" must not be negative')
      !  IF (run_second < 0._wp) CALL finish(method_name,'"run_second" must not be negative')
      !  !
      !  end_datetime = cur_datetime
      !  CALL add_time(run_second,run_minute,run_hour,run_day,end_datetime)
      !  !
      !  cur_datetime_calsec = (REAL(cur_datetime%calday,wp)+cur_datetime%caltime) &
      !    &                   *REAL(cur_datetime%daylen,wp)
      !  end_datetime_calsec = (REAL(end_datetime%calday,wp)+end_datetime%caltime) &
      !    &                   *REAL(end_datetime%daylen,wp)
      !  nsteps=INT((end_datetime_calsec-cur_datetime_calsec)/dtime)
      !-------------------------

    ELSE
      ! Compute nsteps from cur_datetime, end_datetime and dtime
      !
      cur_datetime_calsec = (REAL(time_config%cur_datetime%calday,wp)  &
        +time_config%cur_datetime%caltime   ) &
        * REAL(time_config%cur_datetime%daylen,wp)
      end_datetime_calsec = (REAL(time_config%end_datetime%calday,wp)  &
        +time_config%end_datetime%caltime   ) &
        * REAL(time_config%end_datetime%daylen,wp)

      IF (end_datetime_calsec < cur_datetime_calsec) &
        & CALL finish(TRIM(method_name),'The end date and time must not be '// &
        &            'before the current date and time')

      nsteps=INT((end_datetime_calsec-cur_datetime_calsec)/dtime)
    END IF


    ! Length of this integration is limited by length of the restart cycle.
    !
    IF (nsteps > INT(time_config%dt_restart/dtime)) THEN
      nsteps = INT(time_config%dt_restart/dtime)
      restart_experiment = .TRUE.
    ELSE
      restart_experiment = .FALSE.
    ENDIF
!     nsteps = MIN(nsteps,INT(time_config%dt_restart/dtime))


    CALL message(' ',' ')
    CALL message(method_name,'Initial date and time')
    CALL message(method_name,'---------------------')
    CALL print_datetime_all(time_config%ini_datetime)  ! print all date and time components

    CALL message(' ',' ')
    CALL message(method_name,'End date and time')
    CALL message(method_name,'-----------------')
    CALL print_datetime_all(time_config%end_datetime)  ! print all date and time components

    CALL message(' ',' ')
    CALL message(method_name,'Length of restart cycle')
    CALL message(method_name,'-----------------------')
    WRITE(message_text,'(a,f10.2,a,f16.10,a)') &
         &'dt_restart :',time_config%dt_restart,' seconds =', &
         & time_config%dt_restart/86400._wp, ' days'
    CALL message(method_name,message_text)


    ! Reset the value of dt_checkpoint if it is longer than dt_restart
    ! so that at least one restart file is generated at the end of the cycle.
    !
    dt_checkpoint = MIN(dt_checkpoint,time_config%dt_restart)

    WRITE(message_text,'(a,f10.2,a,f16.10,a)')          &
         &'dt_checkpoint :',dt_checkpoint,' seconds =', &
         & dt_checkpoint/86400._wp, ' days'
    CALL message(method_name,message_text)

  END SUBROUTINE resize_ocean_simulation_length

  

!<Optimize:inUse>
  SUBROUTINE check_thicknesses

    ! ensure, that all used thicknesses are non-zero in case of non-shallow-water.
    !For shallow-water this test makes no sense since dzlev==0 
    IF(iswm_oce==0)THEN
      IF (MINVAL(dzlev_m(1:n_zlev)) <= 0.0_wp) THEN
        CALL finish("check_thicknesses","Found zero or negative thicknesses")
      END IF
    ENDIF  
  END SUBROUTINE check_thicknesses

!<Optimize:inUse>
  SUBROUTINE ocean_crosscheck()

    CHARACTER(len=*), PARAMETER :: method_name =  'mo_ocean_nml_crosscheck:ocean_crosscheck'

    CALL check_parallel_configuration()
    CALL resize_ocean_simulation_length()
    CALL init_grid_configuration

    ! set the patch-related nlev variable to the ocean setup n_ zlev
    nlev = n_zlev

    IF (p_test_run .AND. l_fast_sum ) THEN                      
       CALL warning(method_name, "p_test_run sets l_fast_sum=.f alse.")
       l_fast_sum = .false.                                     
    ENDIF                                                       
    
    SELECT CASE (select_solver)
      CASE (select_gmres)

      CASE (select_restart_gmres, select_restart_mixedPrecision_gmres)

        IF (p_test_run .OR. .NOT. l_fast_sum ) THEN
           CALL warning(method_name, "p_test_run .OR. .NOT. l_fast_sum cannot be used by the restart gmres solver")
           CALL message(method_name, "Using the standard gmres solver")
           select_solver = select_gmres
!        ELSE
!           use_absolute_solver_tolerance = .true.
        ENDIF

      CASE default
        CALL finish(method_name, "Unknown solver")

    END SELECT
    
    
    IF (no_tracer < 1) THEN
      CALL warning("ocean_crosscheck", "no_tracer < 1, use_constant_mixing")
      physics_parameters_type = physics_parameters_Constant_type
    ENDIF

    CALL check_thicknesses

    IF  (RichardsonDiffusion_threshold < convection_InstabilityThreshold) &
      CALL finish (method_name, "RichardsonDiffusion_threshold < convection_InstabilityThreshold")

     
    IF (l_rigid_lid .AND. iswm_oce /= 1) THEN
      CALL finish(method_name, "l_rigid_lid .AND. iswm_oce /= 1")
    ENDIF

    IF (use_duplicated_connectivity) THEN
      use_duplicated_connectivity = .FALSE.
      CALL message(method_name, "Set use_duplicated_connectivity to FALSE")
    ENDIF
    
    IF (isRestart() .AND. write_initial_state) THEN
      CALL warning(method_name, "write_initial_state is disbaled for restarts")
      write_initial_state = .false.
    ENDIF


  END SUBROUTINE ocean_crosscheck


END MODULE mo_ocean_nml_crosscheck
