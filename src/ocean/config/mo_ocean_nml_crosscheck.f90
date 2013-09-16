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
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_ocean_nml_crosscheck

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, message_text, finish, print_value, warning
  USE mo_grid_config,       ONLY: init_grid_configuration
  USE mo_parallel_config,   ONLY: check_parallel_configuration, p_test_run, l_fast_sum
  USE mo_run_config,        ONLY: nsteps, dtime
  USE mo_time_config,       ONLY: time_config, restart_experiment
  USE mo_datetime,          ONLY: add_time, print_datetime_all
  USE mo_io_config,         ONLY: dt_checkpoint
  USE mo_grid_config,       ONLY: grid_rescale_factor
  USE mo_ocean_nml,         ONLY: select_solver, select_restart_gmres, select_gmres, &
    & use_absolute_solver_tolerance

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: oce_crosscheck

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

  !>
  !! Check and, if necessary, adapt simulation lengths
  !! This partially duplicates the code resize_atmo_simulation_length
  !!  not a nice way around it though, needs to be cleaned
  !!
  !! @par Revision History
  !! Initial revision by Hui Wan, MPI (2011-07)
  !!
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
    IF (nsteps/=0) THEN   ! User specified a value

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


  SUBROUTINE oce_crosscheck()

    CHARACTER(len=*), PARAMETER :: method_name =  'mo_ocean_nml_crosscheck:oce_crosscheck'

    CALL check_parallel_configuration()
    CALL resize_ocean_simulation_length()
    CALL init_grid_configuration

    SELECT CASE (select_solver)
      CASE (select_gmres)

      CASE (select_restart_gmres)

        IF (p_test_run .OR. .NOT. l_fast_sum ) THEN
           CALL message(method_name, "p_test_run .OR. .NOT. l_fast_sum cannot be used by the restart gmres solver")
           CALL message(method_name, "Using the standard gmres solver")
           select_solver = select_gmres
        ELSE
           use_absolute_solver_tolerance = .true.
        ENDIF

      CASE default
        CALL finish(method_name, "Unknown solver")

    END SELECT

  END SUBROUTINE oce_crosscheck


END MODULE mo_ocean_nml_crosscheck
