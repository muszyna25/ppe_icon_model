!>
!! @author G. Zaengl
!!
!! @par Revision History
!! Moved configure state from namelists/mo_initicon_nml:
!! F. Prill, DWD (2012-01-31)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_initicon_config

  USE mo_kind,               ONLY: wp
  USE mo_util_string,        ONLY: t_keyword_list, associate_keyword, with_keywords, &
    &                              int2string
  USE mo_io_units,           ONLY: filename_max
  USE mo_impl_constants,     ONLY: max_dom, vname_len, max_var_ml, MAX_CHAR_LENGTH,  &
    &                              MODE_IFSANA, MODE_COMBINED, MODE_COSMODE,         &
    &                              MODE_IAU, MODE_IAU_OLD, MODE_ICONVREMAP
  USE mo_grid_config,        ONLY: l_limited_area
  USE mo_time_config,        ONLY: time_config
  USE mo_datetime,           ONLY: t_datetime
  USE mtime,                 ONLY: timedelta, newTimedelta, deallocateTimedelta,     &
    &                              max_timedelta_str_len, datetime, newDatetime,     &
    &                              deallocateDatetime, OPERATOR(+),                  &
    &                              MAX_DATETIME_STR_LEN, OPERATOR(<=), OPERATOR(>=), &
    &                              getPTStringFromSeconds
  USE mo_mtime_extensions,   ONLY: get_datetime_string
  USE mo_exception,          ONLY: finish, message_text, message

  IMPLICIT NONE

  PRIVATE 


  ! Types
  PUBLIC :: t_initicon_config

  ! Variables
  PUBLIC :: init_mode, nlevatm_in, nlevsoil_in, zpbl1, zpbl2
  PUBLIC :: dt_iau
  PUBLIC :: type_iau_wgt
  PUBLIC :: iterate_iau
  PUBLIC :: l_sst_in
  PUBLIC :: lread_ana
  PUBLIC :: lread_vn
  PUBLIC :: lconsistency_checks
  PUBLIC :: l_coarse2fine_mode
  PUBLIC :: lp2cintp_incr, lp2cintp_sfcana
  PUBLIC :: ltile_coldstart
  PUBLIC :: ltile_init
  PUBLIC :: lvert_remap_fg
  PUBLIC :: lcalc_avg_fg
  PUBLIC :: start_time_avg_fg
  PUBLIC :: end_time_avg_fg
  PUBLIC :: interval_avg_fg
  PUBLIC :: iso8601_start_timedelta_avg_fg
  PUBLIC :: iso8601_end_timedelta_avg_fg
  PUBLIC :: iso8601_interval_avg_fg
  PUBLIC :: ifs2icon_filename
  PUBLIC :: dwdfg_filename
  PUBLIC :: dwdana_filename
  PUBLIC :: filetype
  PUBLIC :: ana_varnames_map_file
  PUBLIC :: init_mode_soil
  PUBLIC :: is_iau_active
  PUBLIC :: iau_wgt_dyn, iau_wgt_adv
  PUBLIC :: rho_incr_filter_wgt
  PUBLIC :: niter_divdamp, niter_diffu
  PUBLIC :: t_timeshift
  PUBLIC :: timeshift
  PUBLIC :: initicon_config
  PUBLIC :: aerosol_fg_present

  ! Subroutines
  PUBLIC :: configure_initicon

  ! Functions
  PUBLIC :: generate_filename
  PUBLIC :: is_avgFG_time


  TYPE t_timeshift
    REAL(wp)                 :: dt_shift
    TYPE(timedelta)          :: mtime_shift
  END TYPE t_timeshift

  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the init_icon preprocessing program
  ! ----------------------------------------------------------------------------
  !
  TYPE :: t_initicon_config
    CHARACTER(LEN=vname_len) :: ana_varlist(max_var_ml) ! list of mandatory analysis fields. 
                                                        ! This list can include a subset or the 
                                                        ! entire set of default analysis fields.
  END TYPE t_initicon_config
  !
  ! probably those which are domain-dependent could be included into aboves type lateron
  !
  INTEGER  :: init_mode     ! initialization mode
  INTEGER  :: nlevsoil_in   ! number of soil levels of input data

  REAL(wp) :: zpbl1, zpbl2  ! AGL heights used for vertical gradient computation
  LOGICAL  :: l_sst_in      ! logical switch, if sea surface temperature is provided as input
  LOGICAL  :: lread_ana     ! If .TRUE., read analysis fields are read from analysis file
                            ! dwdana_filename. If .FALSE., ICON is soleyly started 
                            ! from first guess fields.   

  LOGICAL  :: lconsistency_checks    ! check validity of input fields (FG and ANA)

  LOGICAL  :: l_coarse2fine_mode(max_dom)  ! If true, apply special corrections for interpolation from coarse
                                           ! to fine resolutions over mountainous terrain
  LOGICAL  :: lp2cintp_incr(max_dom) ! If true, perform parent-to-child interpolation of atmospheric data
                                     ! assimilation increments
  LOGICAL  :: lp2cintp_sfcana(max_dom) ! If true, perform parent-to-child interpolation of
                                       ! surface analysis data
  LOGICAL  :: ltile_coldstart  ! If true, initialize tile-based surface fields from first guess with tile-averaged fields

  LOGICAL  :: ltile_init       ! If true, initialize tile-based surface fields from first guess without tiles

  LOGICAL  :: lvert_remap_fg   ! If true, vertical remappting of first guess input is performed

  ! Variables controlling computation of temporally averaged first guess fields for DA
  ! The calculation is switched on by setting end_time > start_time
  REAL(wp) :: start_time_avg_fg   ! start time [s]
  REAL(wp) :: end_time_avg_fg     ! end time [s]
  REAL(wp) :: interval_avg_fg     ! averaging interval [s]
  TYPE(datetime), TARGET :: startdatetime_avgFG ! as start_time_avg_fg but mtime compatible full datetime
  TYPE(datetime), TARGET :: enddatetime_avgFG   ! as end_time_avg_fg   but mtime compatible full datetime
  CHARACTER(len=max_timedelta_str_len):: iso8601_start_timedelta_avg_fg  ! start time in ISO8601 Format (relative)
  CHARACTER(len=max_timedelta_str_len):: iso8601_end_timedelta_avg_fg    ! end time in ISO8601 Format (relative)
  CHARACTER(len=max_timedelta_str_len):: iso8601_interval_avg_fg      ! averaging interval in ISO8601 Format

  INTEGER  :: filetype      ! One of CDI's FILETYPE\_XXX constants. Possible values: 2 (=FILETYPE\_GRB2), 4 (=FILETYPE\_NC2)

  REAL(wp) :: dt_iau        ! Time interval during which incremental analysis update (IAU) is performed [s]. 
                            ! Only required for init_mode=MODE_IAU, MODE_IAU_OLD

  TYPE(t_timeshift) :: &    ! Allows IAU runs to start earlier than the nominal simulation start date 
    &  timeshift            ! without showing up in the output metadata

  INTEGER  :: type_iau_wgt  ! Type of weighting function for IAU.
                            ! 1: Top-hat
                            ! 2: SIN2
                            ! Only required for init_mode=MODE_IAU, MODE_IAU_OLD
  LOGICAL  :: iterate_iau   ! if .TRUE., iterate IAU phase with halved dt_iau in first iteration
  REAL(wp) :: rho_incr_filter_wgt  ! Vertical filtering weight for density increments 
                                   ! Only applicable for init_mode=MODE_IAU, MODE_IAU_OLD

  INTEGER  :: niter_divdamp ! number of divergence damping iterations on wind increment from DA
  INTEGER  :: niter_diffu   ! number of diffusion iterations on wind increment from DA

  ! IFS2ICON input filename, may contain keywords, by default
  ! ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: ifs2icon_filename

  ! DWD-FG input filename, may contain keywords, by default
  ! dwdfg_filename = "<path>dwdFG_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: dwdfg_filename

  ! DWD-ANA input filename, may contain keywords, by default
  ! dwdana_filename = "<path>dwdana_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: dwdana_filename

  ! analysis file: dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names.
  CHARACTER(LEN=filename_max) :: ana_varnames_map_file      

   
  ! ----------------------------------------------------------------------------
  ! Derived variables / variables based on input file contents
  ! ----------------------------------------------------------------------------

  INTEGER :: nlevatm_in(max_dom) = 0  !< number of atmospheric model levels of input data
  LOGICAL :: lread_vn  = .FALSE. !< control variable that specifies if u/v or vn are read as wind field input

  INTEGER :: init_mode_soil     !< initialization mode of soil model (coldstart, warmstart, warmstart+IAU)

  LOGICAL :: is_iau_active = .FALSE.  !< determines whether IAU is active at current time

  LOGICAL :: lcalc_avg_fg           !< determines whether temporally averaged first guess fields are computed

  REAL(wp):: iau_wgt_dyn = 0._wp    !< IAU weight for dynamics fields 
  REAL(wp):: iau_wgt_adv = 0._wp    !< IAU weight for tracer fields

  LOGICAL :: aerosol_fg_present(max_dom) = .FALSE. !< registers if aerosol fields have been read from the first-guess data

  TYPE(t_initicon_config), TARGET :: initicon_config(0:max_dom)

CONTAINS

  !>
  !! setup additional initicon control variables
  !!
  !! Setup of additional initicon control variables depending on the 
  !! initicon-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2013-07-11)
  !!
  SUBROUTINE configure_initicon(dtime)
    !
    REAL(wp), INTENT(IN)        :: dtime       ! advection/fast physics time step
    !
    CHARACTER(len=*), PARAMETER :: routine = 'mo_initicon_config:configure_initicon'
    !
    CHARACTER(len=max_timedelta_str_len) :: PTshift
    TYPE(timedelta), POINTER             :: mtime_shift_local, td_start_time_avg_fg, td_end_time_avg_fg
    CHARACTER(len=max_timedelta_str_len) :: str_start_time_avg_fg, str_end_time_avg_fg
    !
    TYPE(datetime), POINTER              :: inidatetime          ! in mtime format
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)  :: iso8601_ini_datetime ! ISO_8601

    REAL(wp)                             :: zdt_shift            ! rounded dt_shift
    !
    !-----------------------------------------------------------------------
    !


    IF ( ANY((/MODE_IFSANA,MODE_COMBINED,MODE_COSMODE/) == init_mode) ) THEN
       init_mode_soil = 1   ! full coldstart is executed
       ! i.e. w_so_ice and h_snow are re-diagnosed
    ELSE IF (l_limited_area .AND. init_mode == MODE_ICONVREMAP .AND. .NOT. lread_ana) THEN
       init_mode_soil = 1   ! same initialization for limited-area cold start
    ELSE IF ( ANY((/MODE_IAU, MODE_IAU_OLD/) == init_mode) ) THEN
       init_mode_soil = 3  ! warmstart (within assimilation cycle) with analysis increments for h_snow
    ELSE
       init_mode_soil = 2  ! warmstart with full fields for h_snow from snow analysis
    ENDIF

    !
    ! timeshift-operations
    !

    ! Round dt_shift to the nearest integer multiple of the advection time step
    !
    IF (timeshift%dt_shift < 0._wp) THEN
      zdt_shift = REAL(NINT(timeshift%dt_shift/dtime),wp)*dtime
      IF (ABS((timeshift%dt_shift-zdt_shift)/zdt_shift) > 1.e-10_wp) THEN
        WRITE(message_text,'(a,f10.3,a)') '*** WARNING: dt_shift adjusted to ', zdt_shift, &
          &                               ' s in order to be a multiple of the advection time step ***'
        CALL message('',message_text)
      ENDIF
      timeshift%dt_shift = zdt_shift
    ELSE
      iterate_iau = .FALSE. ! IAU iteration is meaningless if the model starts without backward time shift
    END IF
    !
    ! transform timeshift to mtime-format
    !
    CALL getPTStringFromSeconds(timeshift%dt_shift, PTshift)

    !*******************************************************
    ! can be removed, once the new libmtime is available (timedeltaToString)
    IF (TRIM(PTshift)=="-P00.000S") THEN
      PTshift = "-PT00.000S"
    ELSE IF (TRIM(PTshift)=="+P00.000S") THEN
      PTshift = "+PT00.000S"
    ENDIF 
    !********************************************************

    mtime_shift_local => newTimedelta(TRIM(PTshift))
    timeshift%mtime_shift = mtime_shift_local

    ! cleanup
    CALL deallocateTimedelta(mtime_shift_local)

    ! Preparations for first guess averaging
    !
    IF (end_time_avg_fg > start_time_avg_fg) THEN
      lcalc_avg_fg  = .TRUE.   ! determines allocation
    ELSE
      lcalc_avg_fg  = .FALSE.
    ENDIF
    !
    ! transform end_time_avg_fg, start_time_avg_fg to mtime format
    !
    ! create model ini_datetime in ISO_8601 format
    CALL get_datetime_string(iso8601_ini_datetime, time_config%ini_datetime)
    !
    ! convert model ini datetime from ISO_8601 format to type datetime
    inidatetime => newDatetime(TRIM(iso8601_ini_datetime))
    !
    ! get start and end datetime in mtime-format
    CALL getPTStringFromSeconds(start_time_avg_fg, str_start_time_avg_fg)
    td_start_time_avg_fg => newTimedelta(str_start_time_avg_fg)
    CALL getPTStringFromSeconds(end_time_avg_fg, str_end_time_avg_fg)
    td_end_time_avg_fg   => newTimedelta(str_end_time_avg_fg)
    !
    startdatetime_avgFG = inidatetime + td_start_time_avg_fg
    enddatetime_avgFG   = inidatetime + td_end_time_avg_fg
    !
    ! get start and end datetime in ISO_8601 format relative to inidatetime
    ! start time
    CALL getPTStringFromSeconds(start_time_avg_fg, iso8601_start_timedelta_avg_fg)
    ! end time
    CALL getPTStringFromSeconds(end_time_avg_fg, iso8601_end_timedelta_avg_fg)

    !******************************************************* 
    ! can be removed, once the new libmtime is available (timedeltaToString)
    IF (TRIM(iso8601_start_timedelta_avg_fg)=="-P00.000S") THEN
      iso8601_start_timedelta_avg_fg = "-PT00.000S"
    ELSE IF (TRIM(iso8601_start_timedelta_avg_fg)=="+P00.000S") THEN
      iso8601_start_timedelta_avg_fg = "+PT00.000S"
    ENDIF 
    IF (TRIM(iso8601_end_timedelta_avg_fg)=="-P00.000S") THEN
      iso8601_end_timedelta_avg_fg = "-PT00.000S"
    ELSE IF (TRIM(iso8601_end_timedelta_avg_fg)=="+P00.000S") THEN
      iso8601_end_timedelta_avg_fg = "+PT00.000S"
    ENDIF 
    !********************************************************


    !
    ! transform averaging interval to ISO_8601 format
    !
    CALL getPTStringFromSeconds(interval_avg_fg, iso8601_interval_avg_fg)

    ! cleanup
    CALL deallocateDatetime(inidatetime)

  END SUBROUTINE configure_initicon



  FUNCTION generate_filename(input_filename, model_base_dir, &
    &                        nroot, jlev, idom)  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)   :: input_filename, &
      &                               model_base_dir
    INTEGER,          INTENT(IN)   :: nroot, jlev, idom
    CHARACTER(len=MAX_CHAR_LENGTH) :: result_str
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    CALL associate_keyword("<path>",   TRIM(model_base_dir),             keywords)
    CALL associate_keyword("<nroot>",  TRIM(int2string(nroot,"(i1)")),   keywords)
    CALL associate_keyword("<nroot0>", TRIM(int2string(nroot,"(i2.2)")), keywords)
    CALL associate_keyword("<jlev>",   TRIM(int2string(jlev, "(i2.2)")), keywords)
    CALL associate_keyword("<idom>",   TRIM(int2string(idom, "(i2.2)")), keywords)
    ! replace keywords in "input_filename", which is by default
    ! ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
    result_str = TRIM(with_keywords(keywords, TRIM(input_filename)))

  END FUNCTION generate_filename


  !>
  !! Determines, whether it is time for first guess averaging
  !!
  !! Determines, whether it is time for first guess averaging.
  !! The fields u, v, temp, pres and qv can be averaged over a user-specified 
  !! time interval in order to provide them to the data assimilation 
  !! as 'filtered/averaged' first guess fields. 
  !! This function indicates, whether the current model time is within the 
  !! user-specified averaging interval.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-12-17)
  !!
  FUNCTION is_avgFG_time(cur_datetime)
    TYPE(t_datetime), INTENT(INOUT)  :: cur_datetime

    ! local variables
    TYPE(datetime), POINTER :: curdatetime     ! current datetime in mtime format
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: iso8601_cur_datetime ! ISO_8601

    LOGICAL :: is_avgFG_time
    !---------------------------------------------------------------------

    ! create model cur_datetime in ISO_8601 format
    CALL get_datetime_string(iso8601_cur_datetime, cur_datetime)

    !
    ! convert model cur_datetime from ISO_8601 format to type datetime
    curdatetime  => newDatetime(TRIM(iso8601_cur_datetime))

    ! check whether startdatetime_avgFG <= curdatetime <= enddatetime_avgFG
    !
    is_avgFG_time = (curdatetime >= startdatetime_avgFG) .AND.              &
                    (curdatetime <= enddatetime_avgFG    .AND. lcalc_avg_fg )

    ! cleanup
    CALL deallocateDatetime(curdatetime)

  END FUNCTION is_avgFG_time


END MODULE mo_initicon_config
