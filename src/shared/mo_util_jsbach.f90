!>
!! Contains modules that provide interfaces to ICON infrastructure for JSBACH4
!!
!! @par Revision History
!! Moved from adapters directory of JSBACH4 to ICON     by Reiner Schnur (2019-01-30)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#ifndef __NO_JSBACH__
!>
!! @brief Contains interfaces to ICON parallel infrastructure for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version             by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_parallel_iface

  USE mo_mpi,             ONLY: p_comm_work_test, p_comm_work, my_process_is_stdio, my_process_is_mpi_parallel, &
                                get_my_global_mpi_id
  USE mo_parallel_config, ONLY: p_test_run
  USE mo_mpi, ONLY: p_io,    & !< Processor performing I/O
                    p_bcast    !< Broadcasting routine

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_parallel_iface, my_process_is_mpi_parallel, my_process_is_stdio, get_my_global_mpi_id
  PUBLIC :: p_io, mpi_comm, p_bcast

  INTEGER :: mpi_comm

CONTAINS

  SUBROUTINE init_parallel_iface

    IF (p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    END IF

  END SUBROUTINE init_parallel_iface

END MODULE mo_jsb_parallel_iface

!! ==============================================================================================================================
!>
!! @brief Contains interfaces to ICON spatial domains for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                              by Reiner Schnur (2013-04-04)
!! New subroutines and some adaptations for jsb4 standalone   by Julia Nabel (2016-02-20)
!!
MODULE mo_jsb_domain_iface

  USE mo_model_domain,       ONLY: t_patch
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_impl_constants,     ONLY: min_rlcell_int
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_physical_constants, ONLY: earth_radius
  USE mo_math_constants,     ONLY: pi

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_patch, get_host_patch_id
  PUBLIC :: get_grid_filename, get_ntotal, get_ntotal_g, get_dims_g, get_nlat_g, get_nproma, get_nblks
  PUBLIC :: get_blk_start, get_blk_end, get_col_start, get_col_end

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_domain_iface'

CONTAINS

  INTEGER FUNCTION get_host_patch_id(patch)

    TYPE(t_patch), INTENT(in)   :: patch

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_host_patch_id'

    get_host_patch_id = patch%id

  END FUNCTION get_host_patch_id

  FUNCTION get_grid_filename(patch) RESULT(filename)

    USE mo_io_units, ONLY: filename_max

    TYPE(t_patch), INTENT(in) :: patch
    CHARACTER(len=filename_max)               :: filename

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_grid_filename'

    filename = patch%grid_filename

  END FUNCTION get_grid_filename

  INTEGER FUNCTION get_ntotal(patch)

    TYPE(t_patch), INTENT(in) :: patch

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_ntotal'

    get_ntotal = patch%n_patch_cells

  END FUNCTION get_ntotal

  INTEGER FUNCTION get_ntotal_g(patch)

    TYPE(t_patch), INTENT(in) :: patch

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_ntotal_g'

    get_ntotal_g = patch%n_patch_cells_g

  END FUNCTION get_ntotal_g

  FUNCTION get_dims_g(patch)

    TYPE(t_patch), INTENT(in) :: patch

    INTEGER                                   :: get_dims_g(2)

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_dims_g'

    get_dims_g(1) = patch%n_patch_cells_g
    get_dims_g(2) = 1

  END FUNCTION get_dims_g

  FUNCTION get_nlat_g(patch)

    TYPE(t_patch), INTENT(in) :: patch

    INTEGER                                   :: get_nlat_g

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nlat_g'

    ! Derive an effective number of latitudes from characteristic mesh size
    get_nlat_g = pi * earth_radius / patch%geometry_info%mean_characteristic_length

  END FUNCTION get_nlat_g

  INTEGER FUNCTION get_nblks(patch)

    TYPE(t_patch), INTENT(in) :: patch

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nblks'

    get_nblks = patch%nblks_c

  END FUNCTION get_nblks

  INTEGER FUNCTION get_nproma(patch)

    USE mo_parallel_config,    ONLY: nproma

    TYPE(t_patch), INTENT(in) :: patch

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nproma'

    IF (.FALSE.) PRINT*, patch%nblks_c ! only to avoid compiler warnings about unused argument

    get_nproma = nproma
!!$      IF (patch%nblks_c > 1) THEN
!!$        get_nproma = INT( (patch%n_patch_cells - patch%npromz_c)/(patch%nblks_c-1))
!!$      ELSE
!!$        get_nproma = patch%n_patch_cells
!!$      END IF

  END FUNCTION get_nproma

  INTEGER FUNCTION get_blk_start(patch)

    TYPE(t_patch), INTENT(in) :: patch

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_blk_start'

    get_blk_start = patch%cells%start_blk(grf_bdywidth_c + 1, 1)

  END FUNCTION get_blk_start

  INTEGER FUNCTION get_blk_end(patch)

    TYPE(t_patch), INTENT(in) :: patch

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_blk_end'

    get_blk_end = patch%cells%end_blk(min_rlcell_int, MAX(1, patch%n_childdom))

  END FUNCTION get_blk_end

  INTEGER FUNCTION get_col_start(jb, patch)

    INTEGER,       INTENT(in) :: jb
    TYPE(t_patch), INTENT(in) :: patch

    INTEGER :: i_startblk, i_endblk, jcs, jce, rl_start, rl_end

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_col_start'

    i_startblk = get_blk_start(patch)
    i_endblk   = get_blk_end(patch)
    rl_start   = grf_bdywidth_c + 1
    rl_end     = min_rlcell_int

    CALL get_indices_c(patch, jb, i_startblk, i_endblk, jcs, jce, rl_start, rl_end)

    get_col_start = jcs

  END FUNCTION get_col_start

  INTEGER FUNCTION get_col_end(jb, patch)

    INTEGER,       INTENT(in) :: jb
    TYPE(t_patch), INTENT(in) :: patch

    INTEGER :: i_startblk, i_endblk, jcs, jce, rl_start, rl_end

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_col_end'

    i_startblk = get_blk_start(patch)
    i_endblk   = get_blk_end(patch)
    rl_start   = grf_bdywidth_c + 1
    rl_end     = min_rlcell_int

    CALL get_indices_c(patch, jb, i_startblk, i_endblk, jcs, jce, rl_start, rl_end)

    get_col_end = jce

  END FUNCTION get_col_end

END MODULE mo_jsb_domain_iface

!! ==============================================================================================================================
!>
!! @brief Contains interface to ICON namelist handling for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version             by Reiner Schnur (2013-02-06)
!!
MODULE mo_jsb_namelist_iface

  USE mo_exception, ONLY: finish

  USE mo_namelist, ONLY: open_nml_icon => open_nml, close_nml_icon => close_nml, &
                         position_nml_icon => position_nml
  USE mo_io_units, ONLY: nnml
  USE mo_namelist, ONLY: POSITIONED

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: open_nml, POSITIONED, position_nml, close_nml

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_namelist_iface'

CONTAINS

  !>
  !! Wrapper function to make ICON subroutine open_nml a function, as in ECHAM.
  !! Should be replaced when open_nml in ICON gets updated.
  INTEGER FUNCTION open_nml(namelist_filename)

    CHARACTER(len=*), INTENT(in) :: namelist_filename

    LOGICAL :: unit_is_occupied

    CHARACTER(len=*), PARAMETER :: routine = modname//':open_nml'

    INQUIRE(UNIT=nnml, OPENED=unit_is_occupied)
    IF (unit_is_occupied) CALL finish(TRIM(routine), 'namelist unit is already open')

    CALL open_nml_icon(TRIM(namelist_filename))
    open_nml = 1

  END FUNCTION open_nml

  !>
  !! Wrapper function to make ICON subroutine close_nml a function take unit as argument,
  !! as in ECHAM. Should be replaced when close_nml in ICON gets updated.
  SUBROUTINE close_nml(nml_handler)

    INTEGER, INTENT(in) :: nml_handler

    CHARACTER(len=*), PARAMETER :: routine = modname//':close_nml'

    IF (nml_handler /= 1) CALL finish(TRIM(routine), 'Cannot close namelist unit')

    CALL close_nml_icon()

  END SUBROUTINE close_nml

  !>
  !! Wrapper function to make ICON subroutine position_nml a function, as in ECHAM.
  !! Should be replaced when open_nml in ICON gets updated.
  INTEGER FUNCTION position_nml(name, nml_handler, lrewind, status)

    CHARACTER(len=*), INTENT(in)            :: name         ! namelist group name
    INTEGER,          INTENT(in)            :: nml_handler  ! file handler
    LOGICAL,          INTENT(in)  ,OPTIONAL :: lrewind      ! default: true
    INTEGER,          INTENT(out) ,OPTIONAL :: status       ! error return value

    CHARACTER(len=*), PARAMETER :: routine = modname//':position_nml'

    CALL position_nml_icon(name, lrewind=lrewind, status=status)

    position_nml = nnml
    ! avoid compiler warnings about dummy arguments not being used
    IF (nml_handler > 0) CONTINUE

  END FUNCTION position_nml

END MODULE mo_jsb_namelist_iface

!! ==============================================================================================================================
!>
!! @brief Contains interfaces to ICON time control for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                              by Reiner Schnur (2013-02-06)
!! New subroutines and some adaptations for jsb4 standalone   by Julia Nabel (2016-02-20)
!! Added ltimer although this doesn't really belong here      by Reiner Schnur (2018-06-04)
!!
MODULE mo_jsb_time_iface

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message, message_text

  USE mtime,                     ONLY: t_datetime => datetime, newDatetime, deallocateDatetime,      &
    &                                  OPERATOR(+), OPERATOR(*), OPERATOR(==),                       &
    &                                  OPERATOR(<=), OPERATOR(>), OPERATOR(-),                       &
    &                                  divisionquotienttimespan, getDayOfYearFromDateTime,           &
    &                                  getNoOfDaysInMonthDateTime, getNoOfDaysInYearDateTime,        &
    &                                  getTotalMilliSecondsTimeDelta, no_of_sec_in_a_day,            &
    &                                  getNoOfSecondsElapsedInDayDateTime, getTotalSecondsTimeDelta, &
    &                                  divideDatetimeDifferenceInSeconds !, isCurrentEventActive
  USE mo_time_config,            ONLY: time_config !configure_time
  USE mo_time_nml,               ONLY: read_time_namelist
  USE mo_dynamics_config,        ONLY: iequations
  USE mo_impl_constants,         ONLY: ihs_atm_temp, ihs_atm_theta, ishallow_water, &
    &                                  leapfrog_expl, leapfrog_si
  USE mo_ha_dyn_config,          ONLY: ha_dyn_config
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights,              &
    &                                  calculate_time_interpolation_weights
  USE mo_echam_phy_config,       ONLY: echam_phy_tc, dt_zero
  USE mo_run_config,             ONLY: l_timer_host => ltimer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_datetime
  PUBLIC :: get_year_length, get_month_length, get_day_length, get_year_day
  PUBLIC :: get_time_dt, get_time_nsteps, &
            get_time_start, get_time_stop, &
            get_time_current, get_time_previous, is_time_experiment_start, is_time_restart, &
            is_time_ltrig_rad_m1
  PUBLIC :: read_time_namelist !, configure_time
  PUBLIC :: start_timestep, finish_timestep
  PUBLIC :: configure_time_and_events
  PUBLIC :: get_date_components
  PUBLIC :: get_time_interpolation_weights
  PUBLIC :: get_asselin_coef
  PUBLIC :: l_timer_host

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_time_iface'

CONTAINS

  SUBROUTINE configure_time_and_events(model_shortname)

    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: model_shortname

    ! CALL configure_time
    ! avoid compiler warnings about dummy arguments not being used
    IF (PRESENT(model_shortname)) CONTINUE

    RETURN

  END SUBROUTINE configure_time_and_events

  !>
  !! Get time step length
  !!
  !! Function returns the length of the integration time step
  !!
  REAL(wp) FUNCTION get_time_dt()

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_time_dt'
    REAL(wp) :: ztime

    TYPE(t_datetime), POINTER :: reference_dt
    reference_dt => newDatetime("1980-06-01T00:00:00.000")
    ztime = 0.001_wp*getTotalMilliSecondsTimeDelta(time_config%tc_dt_model, reference_dt)
    CALL deallocateDatetime(reference_dt)

    IF (ztime <= 0.0_wp) &
      CALL finish(routine, 'time step not configured yet.')

    get_time_dt = ztime

  END FUNCTION get_time_dt

  !>
  !! Get number of time steps
  !!
  !! Function returns the number of time steps for current model integration
  !!
  INTEGER FUNCTION get_time_nsteps()

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_time_nsteps'

    INTEGER :: isteps

    TYPE(divisionquotienttimespan) :: tq

    CALL divideDatetimeDifferenceInSeconds(time_config%tc_stopdate,    &
         &                                 time_config%tc_startdate,   &
         &                                 time_config%tc_dt_model, tq)
    isteps = INT(tq%quotient)

!@todo: error message is misleading! should be: error in istep calculation or comparable!?
    IF (isteps <= 0) &
      CALL finish(routine, 'number of time steps not configured yet.')

    get_time_nsteps = isteps

  END FUNCTION get_time_nsteps

  TYPE(t_datetime) FUNCTION get_time_start()

    get_time_start = time_config%tc_startdate

  END FUNCTION get_time_start

  TYPE(t_datetime) FUNCTION get_time_stop()

    get_time_stop = time_config%tc_stopdate

  END FUNCTION get_time_stop

  !>
  !! @brief Get current time
  !!
  TYPE(t_datetime) FUNCTION get_time_current()

    get_time_current = time_config%tc_current_date

  END FUNCTION get_time_current

  !>
  !! @brief Get time at previous time step
  !!
  FUNCTION get_time_previous() RESULT(previous)

    TYPE(t_datetime) :: previous

    previous = time_config%tc_current_date + (-1) * time_config%tc_dt_model

  END FUNCTION get_time_previous

  !>
  !! @brief Get time at next time step
  !!
  FUNCTION get_time_next() RESULT(next)

    TYPE(t_datetime) :: next

    next = time_config%tc_current_date + (1) * time_config%tc_dt_model

  END FUNCTION get_time_next

  !
  ! Check whether we are one time step before radiation is calculated in atmosphere
  !
  FUNCTION is_time_ltrig_rad_m1(model_id) RESULT(ltrig_rad_m1)

    INTEGER, INTENT(in) :: model_id
    LOGICAL :: ltrig_rad_m1
    !
    ! In the ICON MPI atmosphere, radiation is triggered based on the radiation time delta and the old  time step.
    ! We therefore use the current time step to check whether we're one time step before the radiation is triggered
    ! in the atmosphere
    TYPE(t_datetime) :: current
    INTEGER          :: dt_rad

    current = get_time_current()
    IF (echam_phy_tc(model_id)%dt_rad > dt_zero) THEN
      dt_rad = getTotalSecondsTimeDelta(echam_phy_tc(model_id)%dt_rad, current)
      ltrig_rad_m1 = MOD(getNoOfSecondsElapsedInDayDateTime(current), dt_rad) == 0
    ELSE
      ltrig_rad_m1 = .FALSE.
    END IF

  END FUNCTION is_time_ltrig_rad_m1

  LOGICAL FUNCTION is_time_experiment_start()

    ! In ICON, the first computed time step is one time step after the start date
    is_time_experiment_start = get_time_previous() == time_config%tc_exp_startdate

  END FUNCTION is_time_experiment_start

  LOGICAL FUNCTION is_time_restart()

    ! In ICON, the first computed time step is one time step after the restart date
    is_time_restart = .NOT. is_time_experiment_start() .AND. (get_time_previous() == time_config%tc_startdate)

  END FUNCTION is_time_restart

  SUBROUTINE start_timestep(istep)

    INTEGER, INTENT(in) :: istep

    CHARACTER(len=*), PARAMETER :: routine = modname//':start_timestep'

    !
    ! Luis:
    !
    ! The following commented region is not allowed as time_config is
    ! a protected varaible. Reiner told me that this is meant for
    ! later use in JSBACH standalone.
    ! time_config%tc_current_date = time_config%tc_current_date &
    !      &                       +time_config%tc_dt_model
    !
    WRITE(message_text,'(a,i10)') 'JSBACH: Start TIME STEP n: ', istep
    CALL message(TRIM(routine),message_text)

  END SUBROUTINE start_timestep

  SUBROUTINE finish_timestep(istep)

    INTEGER, INTENT(in) :: istep

    CHARACTER(len=*), PARAMETER :: routine = modname//':finish_timestep'

    WRITE(message_text,'(a,i10)') 'JSBACH: Finish TIME STEP n: ', istep
    CALL message(TRIM(routine),message_text)

  END SUBROUTINE finish_timestep

  SUBROUTINE get_date_components(this_datetime, year, month, day, hour, minute, second)

    TYPE(t_datetime),  INTENT(in)  :: this_datetime
    INTEGER, OPTIONAL, INTENT(out) :: year, month, day, hour, minute, second

    IF (PRESENT(year))   year   = this_datetime%date%year
    IF (PRESENT(month))  month  = this_datetime%date%month
    IF (PRESENT(day))    day    = this_datetime%date%day
    IF (PRESENT(hour))   hour   = this_datetime%time%hour
    IF (PRESENT(minute)) minute = this_datetime%time%minute
    IF (PRESENT(second)) second = this_datetime%time%second

  END SUBROUTINE get_date_components

  !>
  !! Get length of year
  !!
  !! Function returns the length of a year in days
  !!
  INTEGER FUNCTION get_year_length(yr)
    INTEGER, INTENT(in) :: yr

    TYPE(t_datetime), POINTER :: datetime_tmp

    datetime_tmp => newDatetime(yr, 1, 1, 0, 0, 0, 0)

    get_year_length = getNoOfDaysInYearDateTime(datetime_tmp)

    CALL deallocateDatetime(datetime_tmp)

  END FUNCTION get_year_length

  !>
  !! Get length of month
  !!
  !! Function returns the length of a month in days
  !!
  INTEGER FUNCTION get_month_length(yr, mo)
    INTEGER, INTENT(in) :: yr, mo

    TYPE(t_datetime), POINTER :: datetime_tmp

    datetime_tmp => newDatetime(yr, mo, 1, 0, 0, 0, 0)

    get_month_length = getNoOfDaysInMonthDateTime(datetime_tmp)

    CALL deallocateDatetime(datetime_tmp)

  END FUNCTION get_month_length

  !>
  !! Get length of day
  !!
  !! Function returns the length of a day in seconds
  !!
  REAL(wp) FUNCTION get_day_length()

    get_day_length = REAL(no_of_sec_in_a_day,wp)

  END FUNCTION get_day_length

  !>
  !! Get day in year
  !!
  !! Function returns the day of the year for a given date.
  !! The seconds in the day are given as fraction of the day.
  !!
  REAL(wp) FUNCTION get_year_day(date)

    ! t_datetime points to mtimes datetime
    TYPE(t_datetime),  INTENT(in)  :: date

    get_year_day = REAL(getDayOfYearFromDateTime(date), wp) &
         &        +REAL(getNoOfSecondsElapsedInDayDateTime(date),wp) &
         &        /REAL(no_of_sec_in_a_day,wp)

  END FUNCTION get_year_day

  SUBROUTINE  get_time_interpolation_weights(w1, w2, n1, n2)

    REAL(wp), INTENT(out) :: w1, w2
    INTEGER, INTENT(out)  :: n1, n2

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_time_interpolation_weights'

    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights

!    CALL message('', '')
!    CALL message(routine, 'generate time interpolation weights.')
!    CALL message('', '')
    current_time_interpolation_weights = calculate_time_interpolation_weights(time_config%tc_current_date)
    w1 = current_time_interpolation_weights%weight1
    w2 = current_time_interpolation_weights%weight2
    n1 = current_time_interpolation_weights%month1_index
    n2 = current_time_interpolation_weights%month2_index

  END SUBROUTINE get_time_interpolation_weights

  REAL(wp) FUNCTION get_asselin_coef()

    SELECT CASE(iequations)
    CASE(ishallow_water,ihs_atm_temp,ihs_atm_theta)
      IF (ha_dyn_config%itime_scheme == leapfrog_expl .OR. ha_dyn_config%itime_scheme == leapfrog_si) THEN
        get_asselin_coef = ha_dyn_config%asselin_coeff
      END IF
    CASE DEFAULT
      get_asselin_coef = 0._wp
    END SELECT

  END FUNCTION get_asselin_coef

END MODULE mo_jsb_time_iface

!! ==============================================================================================================================
!>
!! @brief Contains interfaces to ICON io for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                 by Reiner Schnur (2013-04-04)
!! IO namelist subroutine for jsb4 standalone    by Julia Nabel (2016-02-20)
!!
MODULE mo_jsb_io_iface

  USE mo_kind,               ONLY: wp
  USE mo_cf_convention, ONLY: t_cf_var
  USE mo_grib2,         ONLY: t_grib2_var, grib2_var
  USE mo_io_config,     ONLY: lnetcdf_flt64_output
  USE mo_cdi,           ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24, &
                            & ZAXIS_SURFACE, ZAXIS_GENERIC, ZAXIS_DEPTH_BELOW_LAND, &
                            & FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2, &
                            & GRID_UNSTRUCTURED, &
                            & TSTEP_CONSTANT, TSTEP_INSTANT, &
                            & cdiDefMissval
  USE mo_cdi_constants, ONLY: GRID_CELL, GRID_UNSTRUCTURED_CELL
  USE mo_zaxis_type,    ONLY: zaxisTypeList, ZA_SURFACE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_cf_var, t_grib2_var, grib2_var
  PUBLIC :: Get_netcdf_precision
  PUBLIC :: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, DATATYPE_PACK24,  &
            ZAXIS_SURFACE, ZAXIS_GENERIC, ZAXIS_DEPTH_BELOW_LAND,               &
            FILETYPE_NC2, FILETYPE_NC4, FILETYPE_GRB, FILETYPE_GRB2, &
            GRID_CELL, GRID_UNSTRUCTURED, GRID_UNSTRUCTURED_CELL, &
            TSTEP_CONSTANT, TSTEP_INSTANT, &
            Create_zaxis, cdiDefMissval, read_io_namelist
            ! Create_zaxis, Destroy_zaxis, cdiDefMissval


  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_iface'

CONTAINS

  FUNCTION Get_netcdf_precision() RESULT(return_value)

    INTEGER :: return_value

    IF (lnetcdf_flt64_output) THEN
      return_value = DATATYPE_FLT64
    ELSE
      return_value = DATATYPE_FLT32
    END IF

  END FUNCTION Get_netcdf_precision

  SUBROUTINE Create_zaxis(ZaxisID, echamZaxisIdx, cdi_axis_type, length, name, &
    & levels, units, longname, lbounds, ubounds)

    INTEGER,           INTENT(out)   :: ZaxisID
    INTEGER,           INTENT(out)   :: echamZaxisIdx
    INTEGER,           INTENT(in)    :: cdi_axis_type
    INTEGER,           INTENT(in)    :: length
    CHARACTER(len=*),  INTENT(in)    :: name
    REAL(wp),          INTENT(in)    :: levels(:)
    CHARACTER(len=*),  INTENT(in)    :: units
    CHARACTER(len=*),  INTENT(in)    :: longname
    REAL(wp),          INTENT(in)    :: lbounds(:)
    REAL(wp),          INTENT(in)    :: ubounds(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Create_zaxis'

    ! Note: This is only here to suppress compiler warnings about unused dummy arguments
    IF (.FALSE.) THEN
      IF (SIZE(levels)  > 1) CONTINUE
      IF (SIZE(lbounds) > 1) CONTINUE
      IF (SIZE(ubounds) > 1) CONTINUE
      IF (LEN(units)    > 1) CONTINUE
      IF (LEN(longname) > 1) CONTINUE
    END IF

    echamZaxisIdx = -1

    IF (TRIM(name) == 'surface') THEN
      ZaxisID = ZA_SURFACE
      RETURN
    ELSE
      ZaxisID = zaxisTypeList%register(cdi_zaxis_type=cdi_axis_type, is_2D=(length==1))
    END IF

  END SUBROUTINE Create_zaxis

  ! Dummy subroutine, currently not needed for ICON
  SUBROUTINE read_io_namelist(filename)

    CHARACTER(LEN=*), INTENT(in) :: filename
    CHARACTER(len=:), ALLOCATABLE :: filename_loc

    ! Note: This is only here to suppress compiler warnings about unused dummy argument
    !       This function is not implemented for ICON, yet.
    IF (.FALSE.) filename_loc = filename

  END SUBROUTINE read_io_namelist

END MODULE mo_jsb_io_iface

!! ==============================================================================================================================
!>
!! @brief Contains interfaces to ICON NetCDF io for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version             by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_io_netcdf_iface

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message
  USE mo_io_units,           ONLY: filename_max
  USE mo_jsb_domain_iface,   ONLY: t_patch
  USE mo_jsb_parallel_iface, ONLY: my_process_is_mpi_parallel, my_process_is_stdio, p_bcast, p_io, mpi_comm
  USE mo_read_interface, ONLY: read_1D, read_2D, read_2D_time, read_2D_extdim, openInputFile, closeFile, on_cells, &
                               t_stream_id, read_netcdf_broadcast_method

  IMPLICIT NONE

  PUBLIC

  TYPE t_input_file
    INTEGER           :: TYPE                ! 1: ICON infrastructure, 2: ECHAM infrastructure
    CHARACTER(len=filename_max) :: filename  ! Name of input file
    LOGICAL           :: is_open = .FALSE.
    TYPE(t_stream_id) :: stream_id  ! ICON stream_id structure (for > 1D)
    INTEGER :: file_id              ! netcdf file id (for 1D)
  CONTAINS
    PROCEDURE :: Close          => netcdf_close_file
    PROCEDURE :: Read_1d        => netcdf_read_real_1d
    PROCEDURE :: Read_2d        => netcdf_read_real_2d
    PROCEDURE :: Read_2d_time   => netcdf_read_real_2d_time
    PROCEDURE :: Read_2d_extdim => netcdf_read_real_2d_extdim
    PROCEDURE :: Has_dim        => netcdf_file_has_dim
    PROCEDURE :: Has_var        => netcdf_file_has_var
  END TYPE t_input_file

  INTERFACE netcdf_open_input
    MODULE PROCEDURE netcdf_open_input_generic
    MODULE PROCEDURE netcdf_open_input_patch
  END INTERFACE netcdf_open_input

  ! INTERFACE netcdf_read_2d
  !   MODULE PROCEDURE netcdf_read_real_2d
  ! END INTERFACE netcdf_read_2d

  ! INTERFACE netcdf_read_2d_time
  !   MODULE PROCEDURE netcdf_read_real_2d_time
  ! END INTERFACE netcdf_read_2d_time

  ! INTERFACE netcdf_read_2d_extdim
  !   MODULE PROCEDURE netcdf_read_real_2d_extdim
  ! END INTERFACE netcdf_read_2d_extdim

  INCLUDE 'netcdf.inc'

  INTEGER, PARAMETER :: MAX_VAR_DIMS = 16 ! NF_MAX_VAR_DIMS

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_io_netcdf_iface'

CONTAINS

  TYPE(t_input_file) FUNCTION netcdf_open_input_generic(filename)

    CHARACTER(LEN=*),    INTENT(in) :: filename

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_open_input_generic'

    netcdf_open_input_generic%filename = TRIM(filename)

    netcdf_open_input_generic%type      = 1
    ! netcdf_open_input%stream_id = openInputFile(TRIM(filename), patch, read_netcdf_broadcast_method)
    netcdf_open_input_generic%file_id = openInputFile(TRIM(filename))
    netcdf_open_input_generic%is_open = netcdf_open_input_generic%file_id > 0

  END FUNCTION netcdf_open_input_generic

  TYPE(t_input_file) FUNCTION netcdf_open_input_patch(filename, patch)

    CHARACTER(LEN=*),    INTENT(in) :: filename
    TYPE(t_patch),       INTENT(in) :: patch

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_open_input_patch'

    netcdf_open_input_patch%filename = TRIM(filename)

    IF (patch%id > 0) THEN
      netcdf_open_input_patch%type      = 1
      ! netcdf_open_input_patch%stream_id = openInputFile(TRIM(filename), patch, read_netcdf_broadcast_method)
      netcdf_open_input_patch%stream_id = openInputFile(TRIM(filename), patch)
      netcdf_open_input_patch%file_id   = netcdf_open_input_patch%stream_id%file_id
      netcdf_open_input_patch%is_open   = netcdf_open_input_patch%file_id > 0
    ELSE
      CALL finish(TRIM(routine), 'No patch')
    END IF

  END FUNCTION netcdf_open_input_patch

  SUBROUTINE netcdf_close_file(input_file)

    CLASS(t_input_file), INTENT(inout) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_close_file'

    IF (ALLOCATED(input_file%stream_id%read_info)) THEN
      CALL closeFile(input_file%stream_id)
    ELSE
      CALL closeFile(input_file%file_id)
    END IF

      input_file%is_open = .FALSE.

  END SUBROUTINE netcdf_close_file

  FUNCTION netcdf_read_real_1d(input_file, variable_name, fill_array)

    CLASS(t_input_file), INTENT(inout) :: input_file
    CHARACTER(LEN=*),   INTENT(in)     :: variable_name
    REAL(wp), TARGET, OPTIONAL         :: fill_array(:)
    REAL(wp), POINTER                  :: netcdf_read_real_1d(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_read_real_1d'

    NULLIFY(netcdf_read_real_1d)

    IF (input_file%type == 1) THEN
      CALL read_1D(input_file%file_id, TRIM(variable_name), fill_array, netcdf_read_real_1d)
    ELSE IF (input_file%type == 2) THEN
      CALL finish(TRIM(routine), 'Incompatible input file type')
    ELSE
      CALL finish(TRIM(routine), 'Input file type not recognized.')
    END IF

  END FUNCTION netcdf_read_real_1d

  FUNCTION netcdf_read_real_2d(input_file, variable_name, fill_array)

    CLASS(t_input_file), INTENT(inout) :: input_file
    CHARACTER(LEN=*),   INTENT(in)    :: variable_name
    REAL(wp), TARGET, OPTIONAL        :: fill_array(:,:)
    REAL(wp), POINTER                 :: netcdf_read_real_2d(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_read_real_2d'

    NULLIFY(netcdf_read_real_2d)

    IF (input_file%type == 1) THEN
      CALL read_2D(input_file%stream_id, on_cells, TRIM(variable_name), fill_array, netcdf_read_real_2d)
    ELSE IF (input_file%type == 2) THEN
      CALL finish(TRIM(routine), 'Incompatible input file type')
    ELSE
      CALL finish(TRIM(routine), 'Input file type not recognized.')
    END IF

  END FUNCTION netcdf_read_real_2d

  FUNCTION netcdf_read_real_2d_time(input_file, variable_name, fill_array, &
    start_time_step, end_time_step)

    CLASS(t_input_file), INTENT(inout) :: input_file
    CHARACTER(LEN=*),   INTENT(in)    :: variable_name
    REAL(wp), TARGET, OPTIONAL        :: fill_array(:,:,:)
    INTEGER, OPTIONAL,  INTENT(in)    :: start_time_step, end_time_step
    REAL(wp), POINTER                 :: netcdf_read_real_2d_time(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_read_real_2d_time'

    netcdf_read_real_2d_time => netcdf_read_real_2d_extdim( &
      input_file, variable_name, fill_array,    &
      start_extdim=start_time_step, end_extdim=end_time_step, extdim_name="time")

  END FUNCTION netcdf_read_real_2d_time

  FUNCTION netcdf_read_real_2d_extdim(input_file, variable_name, fill_array, &
    start_extdim, end_extdim, extdim_name)

    CLASS(t_input_file),         INTENT(inout) :: input_file
    CHARACTER(LEN=*),           INTENT(in)    :: variable_name
    REAL(wp), TARGET, OPTIONAL                :: fill_array(:,:,:)
    INTEGER, OPTIONAL,          INTENT(in)    :: start_extdim, end_extdim
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: extdim_name
    REAL(wp), POINTER                         :: netcdf_read_real_2d_extdim(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_read_real_2d_extdim'

    NULLIFY(netcdf_read_real_2d_extdim)

    IF (input_file%type == 1) THEN
      CALL read_2D_extdim(input_file%stream_id, on_cells, TRIM(variable_name), fill_array, netcdf_read_real_2d_extdim, &
                          start_extdim=start_extdim, end_extdim=end_extdim, extdim_name=extdim_name)
    ELSE IF (input_file%type == 2) THEN
      CALL finish(TRIM(routine), 'Incompatible input file type')
    ELSE
      CALL finish(TRIM(routine), 'Input file type not recognized.')
    END IF

  END FUNCTION netcdf_read_real_2d_extdim

  LOGICAL FUNCTION netcdf_file_has_dim(input_file, dimname)

    CLASS(t_input_file), INTENT(in) :: input_file
    CHARACTER(LEN=*),    INTENT(in) :: dimname

    INTEGER :: dimid, status

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_file_has_dim'

    IF (my_process_is_stdio()) THEN
      IF (.NOT. input_file%is_open) CALL finish(TRIM(routine), 'NetCDF file not open')

      status = nf_inq_dimid(input_file%file_id, TRIM(dimname), dimid)
      netcdf_file_has_dim = (status == NF_NOERR)
    END IF

    CALL p_bcast(netcdf_file_has_dim, p_io, mpi_comm)

  END FUNCTION netcdf_file_has_dim

  LOGICAL FUNCTION netcdf_file_has_var(input_file, varname)

    CLASS(t_input_file), INTENT(in) :: input_file
    CHARACTER(LEN=*),    INTENT(in) :: varname

    INTEGER :: varid, status

    CHARACTER(len=*), PARAMETER :: routine = modname//':netcdf_file_has_var'

    IF (my_process_is_stdio()) THEN
      IF(.NOT. input_file%is_open) CALL finish(TRIM(routine), 'NetCDF file not open')

      status = nf_inq_varid(input_file%file_id, TRIM(varname), varid)
      netcdf_file_has_var = (status == NF_NOERR)
    END IF

    CALL p_bcast(netcdf_file_has_var, p_io, mpi_comm)

  END FUNCTION netcdf_file_has_var

  SUBROUTINE nf(status, routine, warnonly, silent)

    USE mo_exception, ONLY: em_warn

    INTEGER, INTENT(in)           :: status
    CHARACTER(len=*), INTENT(in)  :: routine
    LOGICAL, INTENT(in), OPTIONAL :: warnonly
    LOGICAL, INTENT(in), OPTIONAL :: silent

    LOGICAL :: lwarnonly, lsilent

    lwarnonly = .FALSE.
    lsilent   = .FALSE.
    IF(PRESENT(warnonly)) lwarnonly = .TRUE.
    IF(PRESENT(silent))   lsilent   = silent

    IF (lsilent) RETURN

    IF (status /= nf_noerr) THEN
      IF (lwarnonly) THEN
        CALL message( TRIM(routine)//' netCDF error', NF_STRERROR(status), &
          & level=em_warn)
      ELSE
        CALL finish( TRIM(routine)//' netCDF error', NF_STRERROR(status))
      ENDIF
    ENDIF

  END SUBROUTINE nf

END MODULE mo_jsb_io_netcdf_iface

!! ==============================================================================================================================
!>
!! @brief Contains grid utilities for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                                   by Reiner Schnur (2013-04-04)
!! New subroutines and some adaptations for jsb4 standalone        by Julia Nabel   (2016-02-20)
!! Separated from mo_jsb_domain_iface to avoid cyclic dependencies by Reiner Schnur (2019-02-04)
!!
MODULE mo_jsb_grid_iface

  USE mo_kind,               ONLY: wp
  USE mo_jsb_domain_iface,   ONLY: t_patch, get_nproma, get_nblks

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_lon, get_lat, get_area

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_grid_iface'

CONTAINS

  FUNCTION get_lon(patch) RESULT(lon)

    TYPE(t_patch), INTENT(in)  :: patch
    REAL(wp), POINTER          :: lon(:,:)

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_lon'

    nproma = get_nproma(patch)
    nblks  = get_nblks (patch)
    ALLOCATE(lon(nproma,nblks))

    lon(:,:) = patch%cells%center(:,:)%lon

  END FUNCTION get_lon

  FUNCTION get_lat(patch) RESULT(lat)

    TYPE(t_patch), INTENT(in)  :: patch
    REAL(wp),      POINTER     :: lat(:,:)

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_lat'

    nproma = get_nproma(patch)
    nblks  = get_nblks (patch)
    ALLOCATE(lat(nproma,nblks))

    lat(:,:) = patch%cells%center(:,:)%lat

  END FUNCTION get_lat

  FUNCTION get_area(patch) RESULT(area)

    TYPE(t_patch), INTENT(in)  :: patch
    REAL(wp),      POINTER     :: area(:,:)

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_area'

    nproma = get_nproma(patch)
    nblks  = get_nblks (patch)
    ALLOCATE(area(nproma,nblks))

    area(:,:) = patch%cells%area(:,:)

  END FUNCTION get_area

END MODULE mo_jsb_grid_iface

!! ==============================================================================================================================
!>
!! @brief Contains interfaces to different ICON utilities for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version               by Reiner Schnur (2016-02-29)
!!
MODULE mo_jsb_utils_iface

  USE mo_fortran_tools, ONLY: assign_if_present

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: assign_if_present

END MODULE mo_jsb_utils_iface

!! ==============================================================================================================================
!>
!! @brief Contains interfaces to ICON varlists for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version               by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_varlist_iface

  USE mo_kind,               ONLY: wp, dp
  USE mo_exception,          ONLY: finish
  USE mo_var_list,           ONLY: new_var_list_icon => new_var_list, &
                                   get_var_list,                      &
                                   add_var_icon => add_var,           &
                                   find_list_element
  USE mo_var_groups,         ONLY: groups
  USE mo_var_metadata_types, ONLY: t_var_metadata, VARNAME_LEN
  USE mo_linked_list,        ONLY: t_var_list, t_list_element

  USE mo_jsb_io_iface, ONLY: t_cf_var, t_grib2_var

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: VARNAME_LEN
  PUBLIC :: t_var_list, t_var_metadata, t_list_element, get_var_list
  PUBLIC :: new_var_list
  PUBLIC :: add_var_list_element_r2d, add_var_list_element_r3d

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_varlist_iface'

CONTAINS

  SUBROUTINE new_var_list (this_list, name, patch_id, output_type, restart_type, &
       &                       post_suf, rest_suf, init_suf, loutput, lrestart,  &
       &                       linitial, table)
    !
    TYPE(t_var_list), POINTER, INTENT(inout) :: this_list    ! anchor
    CHARACTER(len=*), INTENT(in)             :: name         ! name of output var_list
    INTEGER,          INTENT(in)             :: patch_id     ! patch ID
    INTEGER,          INTENT(in), OPTIONAL   :: output_type  ! 'GRIB1', 'GRIB2' or 'NetCDF[12]'
    INTEGER,          INTENT(in), OPTIONAL   :: restart_type ! 'NetCDF[12]'
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: post_suf     ! suffix of output file
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: rest_suf     ! suffix of restart file
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: init_suf     ! suffix of initial file
    LOGICAL,          INTENT(in), OPTIONAL   :: loutput      ! write to  output file
    LOGICAL,          INTENT(in), OPTIONAL   :: lrestart     ! write to restart file
    LOGICAL,          INTENT(in), OPTIONAL   :: linitial     ! read from initial file
    INTEGER,          INTENT(in), OPTIONAL   :: table        ! used only for ECHAM

    IF (PRESENT(table)) CONTINUE ! Only here to avoid compiler warning about "table" not being used

    CALL new_var_list_icon(this_list, name,                                    &
                      output_type=output_type, restart_type=restart_type,      &
                      post_suf=post_suf, rest_suf=rest_suf, init_suf=init_suf, &
                      loutput=loutput, lrestart=lrestart, linitial=linitial,   &
                      patch_id=patch_id                                        &
                     )

  END SUBROUTINE new_var_list

  SUBROUTINE add_var_list_element_r2d(this_list, name, ptr,                             &
    hgrid, vgrid, cf, grib2, code, table, ldims, gdims, levelindx, loutput, lcontainer, &
    lrestart, lrestart_cont, initval_r, isteptype,                                      &
    resetval_r, lmiss, missval_r, tlev_source, info, p5,                                &
    in_groups, verbose, new_element                                                     &
    )

    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:)            ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type
                                                                      ! ICON: zaxis type (ZA_*)
                                                                      ! ECHAM: cdi zaxis ID
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in)           :: code                ! GRIB1 code number
    INTEGER,              INTENT(in)           :: table               ! GRIB1 table number
    INTEGER,              INTENT(in), OPTIONAL :: ldims(2)            ! local dimensions
    INTEGER,              INTENT(in), OPTIONAL :: gdims(2)            ! global dimensions
    INTEGER,              INTENT(in), OPTIONAL :: levelindx           ! info%levelindx for ECHAM
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      ,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    CHARACTER(len=VARNAME_LEN), INTENT(in), OPTIONAL :: in_groups(:)  ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element           ! pointer to new var list element

    TYPE (t_list_element), POINTER :: element

    CHARACTER(len=*), PARAMETER :: routine = modname//':add_var_list_element_r2d'

    ! These variables are not used for ICON, but avoid compiler warnings about dummy arguments not being used
    IF (code > 0) CONTINUE
    IF (table > 0) CONTINUE
    IF (PRESENT(gdims)) CONTINUE
    IF (PRESENT(levelindx)) CONTINUE

    IF (PRESENT(p5)) THEN
      IF (SIZE(p5, DIM=5) > 1) &
        CALL finish(TRIM(routine), 'p5: only four dimensions allowed currently because of ECHAM compatibility')
    END IF

    IF (PRESENT(in_groups)) THEN
      CALL add_var_icon(this_list, TRIM(name), ptr, hgrid, vgrid, cf, grib2, &
        ldims=ldims, loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont,     &
        initval=initval_r, isteptype=isteptype, resetval=resetval_r, lmiss=lmiss, missval=missval_r,             &
        tlev_source=tlev_source, info=info, p5=p5, in_group=groups(in_groups), verbose=verbose, new_element=new_element)
    ELSE
      CALL add_var_icon(this_list, TRIM(name), ptr, hgrid, vgrid, cf, grib2, &
        ldims=ldims, loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont,     &
        initval=initval_r, isteptype=isteptype, resetval=resetval_r, lmiss=lmiss, missval=missval_r,             &
        tlev_source=tlev_source, info=info, p5=p5, verbose=verbose, new_element=new_element)
    END IF
    element => find_list_element(this_list, TRIM(name))
    element%field%info%ndims = 2
    element%field%info%used_dimensions(1:2) = ldims(1:2)

  END SUBROUTINE add_var_list_element_r2d

  SUBROUTINE add_var_list_element_r3d(this_list, name, ptr,                             &
    hgrid, vgrid, cf, grib2, code, table, ldims, gdims, levelindx, loutput, lcontainer, &
    lrestart, lrestart_cont, initval_r, isteptype,                                      &
    resetval_r, lmiss, missval_r, tlev_source, info, p5,                                &
    in_groups, verbose, new_element                                                     &
    )

    TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
    CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
    REAL(wp),             POINTER              :: ptr(:,:,:)          ! reference to field
    INTEGER,              INTENT(in)           :: hgrid               ! horizontal grid type
    INTEGER,              INTENT(in)           :: vgrid               ! vertical grid type
                                                                      ! ICON: zaxis type (ZA_*)
                                                                      ! ECHAM: cdi zaxis ID
    TYPE(t_cf_var),       INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var),    INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER,              INTENT(in)           :: code                ! GRIB1 code number
    INTEGER,              INTENT(in)           :: table               ! GRIB1 table number
    INTEGER,              INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions
    INTEGER,              INTENT(in), OPTIONAL :: gdims(3)            ! global dimensions
    INTEGER,              INTENT(in), OPTIONAL :: levelindx           ! info%levelindx for ECHAM
    LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
    LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
    REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
    INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
    REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
    LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
    REAL(dp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
    INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
    REAL(wp),             TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     ,    OPTIONAL :: p5(:,:,:,:,:)       ! provided pointer
    CHARACTER(len=VARNAME_LEN), INTENT(in), OPTIONAL :: in_groups(:)  ! groups to which a variable belongs
    LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
    TYPE(t_list_element), POINTER, OPTIONAL  :: new_element           ! pointer to new var list element

    TYPE (t_list_element), POINTER :: element

    CHARACTER(len=*), PARAMETER :: routine = modname//':add_var_list_element_r3d'

    ! These variables are not used for ICON, but avoid compiler warnings about dummy arguments not being used
    IF (code > 0) CONTINUE
    IF (table > 0) CONTINUE
    IF (PRESENT(gdims)) CONTINUE
    IF (PRESENT(levelindx)) CONTINUE

    IF (PRESENT(p5)) THEN
      IF (SIZE(p5, DIM=5) > 1) &
        CALL finish(TRIM(routine), 'p5: only four dimensions allowed currently because of ECHAM compatibility')
    END IF

    IF (PRESENT(in_groups)) THEN
      CALL add_var_icon(this_list, TRIM(name), ptr, hgrid, vgrid, cf, grib2, &
        ldims=ldims, loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont,     &
        initval=initval_r, isteptype=isteptype, resetval=resetval_r, lmiss=lmiss, missval=missval_r,             &
        tlev_source=tlev_source, info=info, p5=p5, verbose=verbose, in_group=groups(in_groups), new_element=new_element)
    ELSE
      CALL add_var_icon(this_list, TRIM(name), ptr, hgrid, vgrid, cf, grib2, &
        ldims=ldims, loutput=loutput, lcontainer=lcontainer, lrestart=lrestart, lrestart_cont=lrestart_cont,     &
        initval=initval_r, isteptype=isteptype, resetval=resetval_r, lmiss=lmiss, missval=missval_r,             &
        tlev_source=tlev_source, info=info, p5=p5, verbose=verbose, new_element=new_element)
    END IF
    element => find_list_element(this_list, TRIM(name))
    element%field%info%ndims = 3
    element%field%info%used_dimensions(1:3) = ldims(1:3)

  END SUBROUTINE add_var_list_element_r3d

END MODULE mo_jsb_varlist_iface

!! ==============================================================================================================================
!>
!! @brief Contains interface to ICON mo_jsb_vertical_axes for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version               by Reiner Schnur (2017-03-14)
!!
MODULE mo_jsb_vertical_axes_iface

  USE mo_kind, ONLY: wp !, dp

  USE mo_name_list_output_zaxes_types, ONLY: t_verticalAxis, t_verticalAxisList
  USE mo_cdi_ids,                      ONLY: set_vertical_grid, t_Vgrid
  USE mo_zaxis_type,                   ONLY: zaxisTypeList

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_verticalAxisList, Setup_jsb_vertical_axis, t_Vgrid, Set_jsb_restart_vgrid

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_vertical_axes_iface'

CONTAINS

  SUBROUTINE Setup_jsb_vertical_axis(verticalAxisList, zaxis_id, name, length, longname, units, levels, lbounds, ubounds)

    TYPE(t_verticalAxisList), INTENT(inout) :: verticalAxisList
    INTEGER,                  INTENT(in)    :: zaxis_id
    CHARACTER(len=*),         INTENT(in)    :: name
    INTEGER,                  INTENT(in)    :: length
    CHARACTER(len=*),         INTENT(in)    :: longname
    CHARACTER(len=*),         INTENT(in)    :: units
    REAL(wp),                 INTENT(in)    :: levels(:), lbounds(:), ubounds(:)

    CALL verticalAxisList%append( &
      & t_verticalAxis(zaxisTypeList%getEntry(zaxis_id), &
      &                length,                           &
      &                zaxisLevels   = levels,           &
      &                zaxisLbounds  = lbounds,          &
      &                zaxisUbounds  = ubounds,          &
      &                zaxisUnits    = units,            &
      &                zaxisName     = name,             &
      &                zaxisLongname = longname          &
      &               ) &
      & )

  END SUBROUTINE Setup_jsb_vertical_axis

  SUBROUTINE Set_jsb_restart_vgrid(vgrid_defs, count, zaxis_id, levels)

    TYPE(t_Vgrid), INTENT(inout) :: vgrid_defs(:)
    INTEGER,       INTENT(inout) :: count
    INTEGER, VALUE               :: zaxis_id
    REAL(wp),      INTENT(in)    :: levels(:)

    CALL set_vertical_grid(vgrid_defs, count, zaxis_id, levels)

  END SUBROUTINE Set_jsb_restart_vgrid

END MODULE mo_jsb_vertical_axes_iface

!! ==============================================================================================================================
!>
!! @brief Contains interface to ICON physcial constants for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version               by Reiner Schnur (2016-03-18)
!!
MODULE mo_physical_constants_iface

  USE mo_kind, ONLY: wp

  USE mo_physical_constants, ONLY: &
    grav,            & !< Avg. gravitational acceleration        [m/s2]
    stbo,            & !< Stephan-Boltzmann constant             [W/(m2 K4)]
    argas,           & !< Molar/universal/ideal gas constant     [J/(K mol)]
    zemiss_def,      & !< Long-wave surface emissivity factor    []
    rd,              & !< Gas constant for dry air               [J/(K kg)]
    rv,              & !< Gas constant for water vapor           [J/(K*kg)]
    cpd,             & !< Specific heat at constant pressure     [J/(K kg)]
    cvd,             & !< Specific heat at constant volume       [J/(K kg)]
    rhoh2o,          & !< Density of liquid water                [kg/m3]
    alv,             & !< Latent heat for vaporisation           [J/kg]
    als,             & !< Latent heat for sublimation            [J/kg]
    alf,             & !< Latent heat for fusion                 [J/kg]
    tmelt,           & !< Melting temperature of ice/snow        [K]
    amd,             & !< Molar weight of dry air                [g/mol]
    amco2,           & !< Molar weight of CO2                    [g/mol]
    clw,             & !< Specific heat capacity of liquid water [J/(K kg)]
    ci,              & !< Specific heat capacity of ice          [J/(K kg)]
    !cs,              & !< Specific heat capacity of snow         [J/(K kg)]
    rhoi,            & !< Density of (sea) ice                   [kg/m3]
    rhos,            & !< Density of snow                        [kg/m3]
    ki,              & !< Heat conductivity of ice               [J/(m s K)]
    ks,              & !< Heat conductivity of snow              [J/(m s K)]
    ! Auxiliary constants
    rvd1  => vtmpc1, & !< = rv/rd-1                              []
    cpvd1 => vtmpc2    !< = cpv/cpd-1                            []

  USE mo_echam_vdiff_params, ONLY: cvdifts

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER :: &
    cs=2090._wp        !< Specific heat capacity of snow         [J/(K kg)]

  REAL(wp), PARAMETER :: tpfac1 = cvdifts
  REAL(wp), PARAMETER :: tpfac2 = 1._wp / tpfac1
  REAL(wp), PARAMETER :: tpfac3 = 1._wp - tpfac2

END MODULE mo_physical_constants_iface

!! ==============================================================================================================================
#else

! nag does not like empty files
MODULE util_jsbach
END MODULE util_jsbach

#endif
