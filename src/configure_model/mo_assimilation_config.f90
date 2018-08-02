!>
!! @brief configuration setup for data assimilation
!!
!! @author Klaus Stephan, DWD
!!
!!
!! @par Revision History
!! Initial revision by Klaus Stephan , DWD (2014-12-18)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_assimilation_config
  USE mo_kind,                 ONLY: wp, i8
  USE mo_impl_constants,       ONLY: max_dom, MODE_IAU, MODE_IAU_OLD
  USE mo_exception,            ONLY: message, message_text
  USE mo_run_config,           ONLY: dtime,nsteps
  USE mo_mpi,                  ONLY: my_process_is_stdio
  USE mo_initicon_config,      ONLY: timeshift, init_mode
  USE mo_phy_events,           ONLY: t_phyProcFast, t_phyProcGroup
  USE mo_grid_config,          ONLY: start_time, end_time, DEFAULT_ENDTIME
  USE mo_time_config,          ONLY: time_config
  USE mtime,                   ONLY: datetime, timedelta, newTimedelta, &
    &                                getPTStringFromMS, MAX_TIMEDELTA_STR_LEN, &
    &                                deallocateTimedelta, OPERATOR(+), OPERATOR(>)
  USE mo_model_domain,         ONLY: p_patch 

  IMPLICIT NONE


  PUBLIC 


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for ICON-LHN
  !!--------------------------------------------------------------------------
!  INTEGER, PARAMETER  :: n_noobs  = 36 !Maximum number of times with missing observations

  TYPE :: t_assimilation_config 

    ! Namelist variables for LHN

    LOGICAL                          ::           &
      llhn             ,& ! on/off switch for latent heat nudging (lhn)
      llhnverif        ,& ! on/off switch for verification against radar
      luse_rad         ,& ! read radar data
      lhn_artif        ,& ! apply artificial Gaussian like profile
      lhn_filt         ,& ! vertical filtering of lhn t-increments
      lhn_relax        ,& ! horizontal filtering of lhn t-increments
      lhn_limit        ,& ! impose absolute limit on lhn t-increm. (abs_lhn_lim)
      lhn_no_ttend     ,& ! do not apply the t increments (except for artif points), but only the q increments
      lhn_hum_adj      ,& ! apply a humidity adjustment along with t-increments
      lhn_satad        ,& ! apply saturation adjustment after LHN
      lhn_incloud      ,& ! apply the LHN-scaling in cloudy layers only
      lhn_spqual       ,& ! switch for use of a spatial quality function
      lhn_black        ,& ! use blacklist for radar data
      lhn_qrs          ,& ! calculate the integrated precipitation flux
      lhn_logscale     ,& ! apply logarithmic scaling factors
      lhn_wweight      ,& ! apply a weighting with respect to the mean horizontal wind
      lhn_height       ,& ! use height infos for radar data
      lhn_bright       ,& ! apply bright band detection
      lhn_diag         ,& ! produce more detailed diagnostic output during lhn
      lhn_artif_only      ! apply only artificial temperature profile instead of applying modelled tt_lheat profile
  
    INTEGER ::  &
      nlhn_start       ,& ! start of latent heat nudging period in timesteps
      nlhn_end         ,& ! end of latent heat nudging period in timesteps
      nlhnverif_start  ,& ! start of latent heat nudging period in timesteps
      nlhnverif_end    ,& ! end of latent heat nudging period in timesteps
      nlhn_relax       ,& ! number of interations of horizontal filtering
      nradar           ,& ! max. number of radar stations within input data
      nobs_times          ! number of observation times (i.e. records in radar data file)

    REAL (KIND=wp)                   ::           &
      lhn_coef          ,& ! factor for reduction of lhn t-increments
      lhn_dt_obs        ,& ! time step of input data in minutes
      abs_lhn_lim       ,& ! absolute limit for lhn t-increments (used if lhn_limit)
      fac_lhn_artif     ,& ! factor when artificial profile will applied
      fac_lhn_up        ,& ! limiting factor for upscaling of model heating profile
      fac_lhn_down      ,& ! limiting factor for downscaling model heating profile
      thres_lhn         ,& ! threshold of rain rates to be consinderd within lhn approach
      rqrsgmax          ,& ! ratio of maximum of qrsgflux, needed for reference precipitation
      tt_artif_max      ,& ! maximum of latent heat used as artificial profile
      zlev_artif_max    ,& ! altidude of maximum of artificial profile
      std_artif_max     ,& ! parameter to define vertical width of artifical temperature profile
      start_fadeout        ! time relative to lhn_end, when the lhn coefficient in decreased toward 0

  
   CHARACTER (LEN=100)              ::           &
      radar_in             ,& ! directory for reading radar-files
      radardata_file       ,& ! filename of blacklist for radar data
      blacklist_file       ,& ! filename of blacklist for radar data
      height_file             ! dxheight_file_name
  
   TYPE(t_phyProcFast)  :: dass_lhn   !> event for LHN
   TYPE(t_phyProcFast)  :: dass_lhn_verif   !> event for LHN
   TYPE(t_phyProcGroup) :: dass_g        

  END TYPE t_assimilation_config

  TYPE(t_assimilation_config), TARGET :: assimilation_config(0:max_dom)


  CONTAINS

  SUBROUTINE configure_lhn(jg)
   INTEGER, INTENT(IN) :: jg          !< patch
   CHARACTER (LEN=100)              ::           &
     filepath
   LOGICAL  :: lf_exist,lb_exist,lh_exist
   INTEGER  :: nobs,nt_end,nt_start, idt_shift


    ! local
    TYPE(timedelta), POINTER        :: eventInterval    => NULL()
    TYPE(datetime)                  :: eventEndDate_proc     ! process-specific end date
                                                             ! set to startDate, if process is disabled
    TYPE(timedelta), POINTER        :: plusSlack    => NULL()
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: dt_ass_str       ! physics timestep (PT-format)

    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_start_str
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_end_str
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_dt_str
    TYPE(datetime)           :: domStartDate, domEndDate
    TYPE(datetime)           :: eventStartDate, eventEndDate
    TYPE(timedelta), POINTER :: td_start, td_end, td_dt   => NULL()

    REAL(wp)            :: dt_ass            !< assimilation time intervals



    IF (ANY((/MODE_IAU,MODE_IAU_OLD/)==init_mode)) THEN
       idt_shift=INT(timeshift%dt_shift/dtime)
    ELSE
       idt_shift=0
    ENDIF

    IF (jg > 1) THEN
       dt_ass = (dtime/2._wp**(p_patch(jg)%level -  p_patch(1)%level))            !seconds
    ELSE
       dt_ass = dtime
    ENDIF

    IF (assimilation_config(jg)%llhn .OR. assimilation_config(jg)%llhnverif) THEN
       filepath=TRIM(assimilation_config(jg)%radar_in)//TRIM(assimilation_config(jg)%radardata_file)
       INQUIRE (FILE=filepath, EXIST=lf_exist)
       IF (.NOT. lf_exist) THEN
          WRITE (message_text,'(3a)') "radar data file",TRIM(filepath)," not found, llhn/llhnverif set to false!"
          CALL message('configure_lhn: ',message_text)
          assimilation_config(jg)%llhn = lf_exist
          assimilation_config(jg)%llhnverif = lf_exist
       ENDIF

       IF (assimilation_config(jg)%lhn_black) THEN
          filepath=assimilation_config(jg)%radar_in(1:LEN_TRIM(assimilation_config(jg)%radar_in))&
                  &//assimilation_config(jg)%blacklist_file(1:LEN_TRIM(assimilation_config(jg)%blacklist_file))
          INQUIRE(file=filepath,EXIST=lb_exist)
          IF (.not. lb_exist) THEN
             WRITE (message_text,'(3a)') "radar blacklist file",TRIM(filepath)," not found, lhn_black set to false!"
             CALL message('configure_lhn: ',message_text)
             assimilation_config(jg)%lhn_black = lb_exist
          ENDIF
       ENDIF
       IF (assimilation_config(jg)%lhn_height) THEN
          filepath=assimilation_config(jg)%radar_in(1:LEN_TRIM(assimilation_config(jg)%radar_in))&
                  &//assimilation_config(jg)%height_file(1:LEN_TRIM(assimilation_config(jg)%height_file))
          INQUIRE(file=filepath,EXIST=lh_exist)
          IF (.not. lh_exist) THEN
             WRITE (message_text,'(3a)') "radar height file",TRIM(filepath)," not found, lhn_height set to false!"
             CALL message('configure_lhn: ',message_text)
             assimilation_config(jg)%lhn_height = lh_exist
          ENDIF
       ENDIF
   
       assimilation_config(jg)%lhn_bright = (assimilation_config(jg)%lhn_height .AND. assimilation_config(jg)%lhn_bright)
    ENDIF

    IF (assimilation_config(jg)%lhn_logscale) THEN
        assimilation_config(jg)%fac_lhn_down   = 1.0_wp + log(assimilation_config(jg)%fac_lhn_down)
        assimilation_config(jg)%fac_lhn_up     = 1.0_wp + log(assimilation_config(jg)%fac_lhn_up)
        assimilation_config(jg)%fac_lhn_artif  = 1.0_wp + log(assimilation_config(jg)%fac_lhn_artif)
    ENDIF


    nt_end = MAX(assimilation_config(jg)%nlhn_end,assimilation_config(jg)%nlhnverif_end)
    nt_end = MIN(nt_end,nsteps-idt_shift)
    nt_start = MIN(assimilation_config(jg)%nlhn_start,assimilation_config(jg)%nlhnverif_start)
    nt_start = MIN(nt_end - 1,nt_start)
    nobs = NINT((REAL(nt_end-nt_start)+3600._wp)/60._wp/assimilation_config(jg)%lhn_dt_obs) ! consider one hour more to be safe
    IF (nobs > 0) THEN
      assimilation_config(jg)%nobs_times = nobs
    ELSE
      assimilation_config(jg)%llhn=.false.
      assimilation_config(jg)%llhnverif=.false.
    ENDIF

    assimilation_config(jg)%luse_rad=(assimilation_config(jg)%llhn.OR.assimilation_config(jg)%llhnverif)

    if (assimilation_config(jg)%llhnverif) assimilation_config(jg)%lhn_diag=.true.

    ! LHN event

    ! Events are triggered between [actual_trigger_time, actual_trigger_time + plus_slack]
    plusSlack =>newTimedelta("PT0S")

    ! compute event start date
    IF (jg > 1) THEN
      CALL getPTStringFromMS(INT(assimilation_config(jg)%nlhn_start*1000._wp,i8), td_start_str)
      td_start => newTimedelta(td_start_str)
      domStartDate = time_config%tc_exp_startdate + td_start
      CALL deallocateTimedelta(td_start)
    ELSE
      domStartDate = time_config%tc_exp_startdate
      ! take care of possibe IAU-Timeshift
      IF (timeshift%dt_shift < 0._wp) THEN
        domStartDate = domStartDate + timeshift%mtime_shift
      ENDIF
    ENDIF
    ! Note that the model time is updated at the beginning of a timestep.
    !
    ! by adding td_dt we make sure, that all events are triggered 
    ! during the first integration step of the given patch.
    CALL getPTStringFromMS(INT(dt_ass*1000._wp,i8), td_dt_str)
    td_dt => newTimedelta(td_dt_str)
   !
    eventStartDate = domStartDate + td_dt
!    eventStartDate = domStartDate 
    CALL deallocateTimedelta(td_dt)

    ! compute event end date
    IF (jg > 1) THEN
      IF (end_time(jg) /= DEFAULT_ENDTIME) THEN
!        CALL getPTStringFromMS(INT(end_time(jg)*1000._wp,i8), td_end_str)
        CALL getPTStringFromMS(INT(assimilation_config(jg)%nlhn_end*1000._wp,i8), td_end_str)
        td_end => newTimedelta(td_end_str)
        domEndDate = time_config%tc_exp_startdate + td_end
        CALL deallocateTimedelta(td_end)
        ! make sure that eventEndDate<=tc_exp_stopdate
        IF (domEndDate > time_config%tc_exp_stopdate) &
          &  domEndDate = time_config%tc_exp_stopdate
      ELSE
        domEndDate = time_config%tc_exp_stopdate
      ENDIF
    ELSE
      domEndDate = time_config%tc_exp_stopdate
    ENDIF
    eventEndDate = domEndDate

    eventEndDate_proc = MERGE(eventEndDate, eventStartDate, assimilation_config(jg)%llhn)
    !
    CALL getPTStringFromMS(INT(dt_ass*1000._wp,i8), dt_ass_str)
    eventInterval=>newTimedelta(dt_ass_str)

    CALL assimilation_config(jg)%dass_g%construct(grpName=TRIM("DASS_G"), pid=jg, grpSize=2)
    !
    CALL assimilation_config(jg)%dass_lhn%initialize(                               &
      &                     name        = 'dass_lhn',                           & !in
      &                     id          = 1,                                    & !in
      &                     is_enabled  = assimilation_config(jg)%llhn,             & !in
      &                     startDate   = eventStartDate,                       & !in
      &                     endDate     = eventEndDate_proc,                    & !in
      &                     dt          = eventInterval,                        & !in
      &                     plusSlack   = plusSlack,                            & !in
      &                     optInclStart= .TRUE.                                ) !in
    !
    CALL assimilation_config(jg)%dass_g%addToGroup(                          &
      &                     phyProc = assimilation_config(jg)%dass_lhn)

    ! compute event start date
    IF (jg > 1) THEN
      CALL getPTStringFromMS(INT(assimilation_config(jg)%nlhnverif_start*1000._wp,i8), td_start_str)
      td_start => newTimedelta(td_start_str)
      domStartDate = time_config%tc_exp_startdate + td_start
      CALL deallocateTimedelta(td_start)
    ELSE
      domStartDate = time_config%tc_exp_startdate
      ! take care of possibe IAU-Timeshift
      IF (timeshift%dt_shift < 0._wp) THEN
        domStartDate = domStartDate + timeshift%mtime_shift
      ENDIF
    ENDIF
    ! Note that the model time is updated at the beginning of a timestep.
    !
    ! by adding td_dt we make sure, that all events are triggered 
    ! during the first integration step of the given patch.
    CALL getPTStringFromMS(INT(dt_ass*1000._wp,i8), td_dt_str)
    td_dt => newTimedelta(td_dt_str)
   !
    eventStartDate = domStartDate + td_dt
!    eventStartDate = domStartDate 
    CALL deallocateTimedelta(td_dt)

    ! compute event end date
    IF (jg > 1) THEN
      IF (end_time(jg) /= DEFAULT_ENDTIME) THEN
!        CALL getPTStringFromMS(INT(end_time(jg)*1000._wp,i8), td_end_str)
        CALL getPTStringFromMS(INT(assimilation_config(jg)%nlhnverif_end*1000._wp,i8), td_end_str)
        td_end => newTimedelta(td_end_str)
        domEndDate = time_config%tc_exp_startdate + td_end
        CALL deallocateTimedelta(td_end)
        ! make sure that eventEndDate<=tc_exp_stopdate
        IF (domEndDate > time_config%tc_exp_stopdate) &
          &  domEndDate = time_config%tc_exp_stopdate
      ELSE
        domEndDate = time_config%tc_exp_stopdate
      ENDIF
    ELSE
      domEndDate = time_config%tc_exp_stopdate
    ENDIF
    eventEndDate = domEndDate

    eventEndDate_proc = MERGE(eventEndDate, eventStartDate, assimilation_config(jg)%llhnverif)

    CALL assimilation_config(jg)%dass_lhn_verif%initialize(                         &
      &                     name        = 'dass_lhn_verif',                     & !in
      &                     id          = 2,                                    & !in
      &                     is_enabled  = assimilation_config(jg)%llhnverif,        & !in
      &                     startDate   = eventStartDate,                       & !in
      &                     endDate     = eventEndDate_proc,                    & !in
      &                     dt          = eventInterval,                        & !in
      &                     plusSlack   = plusSlack,                            & !in
      &                     optInclStart= .TRUE.                                ) !in
!    ! add to physics group
    CALL assimilation_config(jg)%dass_g%addToGroup(                          &
      &                     phyProc = assimilation_config(jg)%dass_lhn_verif)
    CALL deallocateTimedelta(eventInterval)

    CALL assimilation_config(jg)%dass_g%printSetup()

  END SUBROUTINE configure_lhn

END MODULE mo_assimilation_config
