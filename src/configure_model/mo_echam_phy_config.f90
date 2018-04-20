!>
!! Configuration of the ECHAM physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!     First version by Marco Giorgetta, MPI-M (2017-04)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_phy_config

  USE mo_exception     ,ONLY: message, print_value, finish

  USE mtime            ,ONLY: OPERATOR(<), OPERATOR(>), OPERATOR(==),                                                 &
       &                      datetime , newDatetime , datetimeToString , max_datetime_str_len ,                      &
       &                      timedelta, newTimedelta, timedeltaToString, max_timedelta_str_len, deallocateTimeDelta, &
       &                      event    , newEvent    , eventGroup         , addEventToEventGroup
  USE mo_event_manager ,ONLY: addEventGroup, getEventGroup, printEventGroup

  USE mo_impl_constants,ONLY: max_dom
  USE mo_grid_config   ,ONLY: n_dom
  USE mo_master_config ,ONLY: experimentStartDate, experimentStopDate

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::                    name    !< name for this unit

  PUBLIC ::         echam_phy_config   !< user specified configuration parameters
  PUBLIC ::         echam_phy_tc       !< derived time control (tc) parameters
  PUBLIC ::    init_echam_phy_config   !< allocate and initialize echam_phy_config
  PUBLIC ::    eval_echam_phy_config   !< evaluate echam_phy_config
  PUBLIC ::    eval_echam_phy_tc       !< evaluate echam_phy_tc
  PUBLIC ::   print_echam_phy_config   !< print out

  PUBLIC ::                  dt_zero   !< a zero (0 sec) time interval 'PT0S'

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_phy'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the ECHAM physics package
  !!
  TYPE t_echam_phy_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! dynamics physics coupling
     LOGICAL  :: ldcphycpl  !< determines the coupling between the dynamics and the phyiscs
     !                      !  .FALSE.: dynamics and physics update sequentially
     !                      !  .TRUE. : dynamics uses physics forcing for updating
     !
     LOGICAL  :: lparamcpl  !< determines the coupling between the parameterizations
     !                      !  .FALSE.: parameterizations do not update the physics state
     !                      !           so that all params get the same input state
     !                      !  .TRUE. : parameterizations update the physics state
     !                      !           so that each param. provides a new prov. state
     !
     LOGICAL  :: ldrymoist  !  .true. : use dry   air mass as conserved reference air mass
     !                      !  .false.: use moist air mass as conserved reference air mass
     !
     !
     ! time control of processes
     ! - sd = start date
     ! - ed = end date
     ! - dt = time interval
     !
     ! forcing control of processes
     ! (to replace later ldcphycpl, which controls all params. together, by a param. specific control)
     ! - fc
     !   fc = 0: diagnostic, tendencies of the param. are not used in integration
     !   fc = 1: prognostic, tendencies of the param. are used to update the model state
     !   fc = 2: prognostic, tendencies of the param. are used as forcing in the dynamical core 
     !
     ! atmospheric physics
     CHARACTER(len=max_timedelta_str_len) :: dt_rad  !< time  step of LW radiation
     CHARACTER(len=max_datetime_str_len ) :: sd_rad  !< start time of LW radiation
     CHARACTER(len=max_datetime_str_len ) :: ed_rad  !< end   time of LW radiation
     INTEGER                              :: fc_rht
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_vdf  !< time  step of vertical diffusion
     CHARACTER(len=max_datetime_str_len ) :: sd_vdf  !< start time of vertical diffusion
     CHARACTER(len=max_datetime_str_len ) :: ed_vdf  !< end   time of vertical diffusion
     INTEGER                              :: fc_vdf
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_cnv  !< time  step of cumulus convection
     CHARACTER(len=max_datetime_str_len ) :: sd_cnv  !< start time of cumulus convection
     CHARACTER(len=max_datetime_str_len ) :: ed_cnv  !< end   time of cumulus convection
     INTEGER                              :: fc_cnv
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_cld  !< time  step of cloud microphysics
     CHARACTER(len=max_datetime_str_len ) :: sd_cld  !< start time of cloud microphysics
     CHARACTER(len=max_datetime_str_len ) :: ed_cld  !< end   time of cloud microphysics
     INTEGER                              :: fc_cld
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_gwd  !< time  step of atmospheric gravity wave drag
     CHARACTER(len=max_datetime_str_len ) :: sd_gwd  !< start time of atmospheric gravity wave drag
     CHARACTER(len=max_datetime_str_len ) :: ed_gwd  !< end   time of atmospheric gravity wave drag
     INTEGER                              :: fc_gwd
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_sso  !< time  step of sub grid scale orogr. effects
     CHARACTER(len=max_datetime_str_len ) :: sd_sso  !< start time of sub grid scale orogr. effects
     CHARACTER(len=max_datetime_str_len ) :: ed_sso  !< end   time of sub grid scale orogr. effects
     INTEGER                              :: fc_sso
     !
     ! atmospheric chemistry
     CHARACTER(len=max_timedelta_str_len) :: dt_mox  !< time  step of methan oxidation and water vapor photolysis
     CHARACTER(len=max_datetime_str_len ) :: sd_mox  !< start time of methan oxidation and water vapor photolysis
     CHARACTER(len=max_datetime_str_len ) :: ed_mox  !< end   time of methan oxidation and water vapor photolysis
     INTEGER                              :: fc_mox
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_car  !< time  step of lin. Cariolle ozone chemistry
     CHARACTER(len=max_datetime_str_len ) :: sd_car  !< start time of lin. Cariolle ozone chemistry
     CHARACTER(len=max_datetime_str_len ) :: ed_car  !< end   time of lin. Cariolle ozone chemistry
     INTEGER                              :: fc_car
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_art  !< time  step of ART chemistry
     CHARACTER(len=max_datetime_str_len ) :: sd_art  !< start time of ART chemistry
     CHARACTER(len=max_datetime_str_len ) :: ed_art  !< end   time of ART chemistry
     INTEGER                              :: fc_art
     !
     ! surface
     LOGICAL                              :: lmlo    !< .true. for mixed layer ocean
     LOGICAL                              :: lice    !< .true. for sea-ice temperature calculation
     LOGICAL                              :: ljsb    !< .true. for calculating the JSBACH land surface
     LOGICAL                              :: llake   !< .true. for using lakes in JSBACH
     LOGICAL                              :: lamip   !< .true. for AMIP simulations
     LOGICAL                              :: lcpl_co2_atmoce !< .true. for coupling of co2 atmo/ocean
     !
  END TYPE t_echam_phy_config


  TYPE t_echam_phy_tc
     !
     ! mtime datetime and time delta, and events
     ! -----------------------------------------
     !
     ! - sd = start date
     ! - ed = end date
     ! - dt = time interval
     ! - ev = event
     !
     ! atmospheric physics
     TYPE(timedelta), POINTER :: dt_rad
     TYPE(datetime ), POINTER :: sd_rad   
     TYPE(datetime ), POINTER :: ed_rad   
     TYPE(event    ), POINTER :: ev_rad
     !
     TYPE(timedelta), POINTER :: dt_vdf
     TYPE(datetime ), POINTER :: sd_vdf
     TYPE(datetime ), POINTER :: ed_vdf
     TYPE(event    ), POINTER :: ev_vdf
     !
     TYPE(timedelta), POINTER :: dt_cnv
     TYPE(datetime ), POINTER :: sd_cnv
     TYPE(datetime ), POINTER :: ed_cnv
     TYPE(event    ), POINTER :: ev_cnv
     !
     TYPE(timedelta), POINTER :: dt_cld
     TYPE(datetime ), POINTER :: sd_cld
     TYPE(datetime ), POINTER :: ed_cld
     TYPE(event    ), POINTER :: ev_cld
     !
     TYPE(timedelta), POINTER :: dt_gwd
     TYPE(datetime ), POINTER :: sd_gwd
     TYPE(datetime ), POINTER :: ed_gwd
     TYPE(event    ), POINTER :: ev_gwd
     !
     TYPE(timedelta), POINTER :: dt_sso
     TYPE(datetime ), POINTER :: sd_sso
     TYPE(datetime ), POINTER :: ed_sso
     TYPE(event    ), POINTER :: ev_sso
     !
     ! atmospheric chemistry
     TYPE(timedelta), POINTER :: dt_mox
     TYPE(datetime ), POINTER :: sd_mox
     TYPE(datetime ), POINTER :: ed_mox
     TYPE(event    ), POINTER :: ev_mox
     !
     TYPE(timedelta), POINTER :: dt_car
     TYPE(datetime ), POINTER :: sd_car
     TYPE(datetime ), POINTER :: ed_car
     TYPE(event    ), POINTER :: ev_car
     !
     TYPE(timedelta), POINTER :: dt_art
     TYPE(datetime ), POINTER :: sd_art
     TYPE(datetime ), POINTER :: ed_art
     TYPE(event    ), POINTER :: ev_art
     !
  END TYPE t_echam_phy_tc

  !>
  !! Configuration/logicals/timecontrol state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_phy_config), TARGET :: echam_phy_config  (max_dom)
  TYPE(t_echam_phy_tc)    , TARGET :: echam_phy_tc      (max_dom)
  
  !>
  !! Events and event group
  !!
  INTEGER                   :: echam_phy_events
  TYPE(eventGroup), POINTER :: echam_phy_event_group

  !>
  !! For convenience
  !!
  TYPE(timedelta) , POINTER :: dt_zero

CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_phy_config
    !
    dt_zero =>  newTimedelta ('PT0S')
    !
    ! ECHAM physics configuration
    ! ---------------------------
    !
    ! dynamics physics coupling
    echam_phy_config(:)%ldcphycpl = .FALSE.
    echam_phy_config(:)%lparamcpl = .TRUE.
    echam_phy_config(:)%ldrymoist = .FALSE.
    !
    ! time control parameters
    echam_phy_config(:)% dt_rad = ''
    echam_phy_config(:)% sd_rad = ''
    echam_phy_config(:)% ed_rad = ''
    echam_phy_config(:)% fc_rht = 1
    !
    echam_phy_config(:)% dt_vdf = ''
    echam_phy_config(:)% sd_vdf = ''
    echam_phy_config(:)% ed_vdf = ''
    echam_phy_config(:)% fc_vdf = 1
    !
    echam_phy_config(:)% dt_cnv = ''
    echam_phy_config(:)% sd_cnv = ''
    echam_phy_config(:)% ed_cnv = ''
    echam_phy_config(:)% fc_cnv = 1
    !
    echam_phy_config(:)% dt_cld = ''
    echam_phy_config(:)% sd_cld = ''
    echam_phy_config(:)% ed_cld = ''
    echam_phy_config(:)% fc_cld = 1
    !
    echam_phy_config(:)% dt_gwd = ''
    echam_phy_config(:)% sd_gwd = ''
    echam_phy_config(:)% ed_gwd = ''
    echam_phy_config(:)% fc_gwd = 1
    !
    echam_phy_config(:)% dt_sso = ''
    echam_phy_config(:)% sd_sso = ''
    echam_phy_config(:)% ed_sso = ''
    echam_phy_config(:)% fc_sso = 1
    !
    echam_phy_config(:)% dt_mox = ''
    echam_phy_config(:)% sd_mox = ''
    echam_phy_config(:)% ed_mox = ''
    echam_phy_config(:)% fc_mox = 1
    !
    echam_phy_config(:)% dt_car = ''
    echam_phy_config(:)% sd_car = ''
    echam_phy_config(:)% ed_car = ''
    echam_phy_config(:)% fc_car = 1
    !
    echam_phy_config(:)% dt_art = ''
    echam_phy_config(:)% sd_art = ''
    echam_phy_config(:)% ed_art = ''
    echam_phy_config(:)% fc_art = 1
    !
    ! logical switches
    echam_phy_config(:)% ljsb  = .FALSE.
    echam_phy_config(:)% llake = .FALSE.
    echam_phy_config(:)% lamip = .FALSE.
    echam_phy_config(:)% lmlo  = .FALSE.
    echam_phy_config(:)% lice  = .FALSE.
    echam_phy_config(:)% lcpl_co2_atmoce  = .FALSE.
    !
  END SUBROUTINE init_echam_phy_config

  !----

  !>
  !! Check the echam_phy_config state
  !!
  SUBROUTINE eval_echam_phy_config
    !
    INTEGER                 :: jg
    CHARACTER(LEN=2)        :: cg
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL eval_echam_phy_config_details(TRIM(cg),                'rad' ,&
            &                             echam_phy_config (jg)% dt_rad  ,&
            &                             echam_phy_config (jg)% sd_rad  ,&
            &                             echam_phy_config (jg)% ed_rad  ,&
            &                             echam_phy_config (jg)% fc_rht  )
       !
       CALL eval_echam_phy_config_details(TRIM(cg),                'vdf' ,&
            &                             echam_phy_config (jg)% dt_vdf  ,&
            &                             echam_phy_config (jg)% sd_vdf  ,&
            &                             echam_phy_config (jg)% ed_vdf  ,&
            &                             echam_phy_config (jg)% fc_vdf  )
       !
       CALL eval_echam_phy_config_details(TRIM(cg),                'cnv' ,&
            &                             echam_phy_config (jg)% dt_cnv  ,&
            &                             echam_phy_config (jg)% sd_cnv  ,&
            &                             echam_phy_config (jg)% ed_cnv  ,&
            &                             echam_phy_config (jg)% fc_cnv  )
       !
       CALL eval_echam_phy_config_details(TRIM(cg),                'cld' ,&
            &                             echam_phy_config (jg)% dt_cld  ,&
            &                             echam_phy_config (jg)% sd_cld  ,&
            &                             echam_phy_config (jg)% ed_cld  ,&
            &                             echam_phy_config (jg)% fc_cld  )
       !
       CALL eval_echam_phy_config_details(TRIM(cg),                'gwd' ,&
            &                             echam_phy_config (jg)% dt_gwd  ,&
            &                             echam_phy_config (jg)% sd_gwd  ,&
            &                             echam_phy_config (jg)% ed_gwd  ,&
            &                             echam_phy_config (jg)% fc_gwd  )
       !
       CALL eval_echam_phy_config_details(TRIM(cg),                'sso' ,&
            &                             echam_phy_config (jg)% dt_sso  ,&
            &                             echam_phy_config (jg)% sd_sso  ,&
            &                             echam_phy_config (jg)% ed_sso  ,&
            &                             echam_phy_config (jg)% fc_sso  )
       !
       CALL eval_echam_phy_config_details(TRIM(cg),                'mox' ,&
            &                             echam_phy_config (jg)% dt_mox  ,&
            &                             echam_phy_config (jg)% sd_mox  ,&
            &                             echam_phy_config (jg)% ed_mox  ,&
            &                             echam_phy_config (jg)% fc_mox  )
       !
       CALL eval_echam_phy_config_details(TRIM(cg),                'car' ,&
            &                             echam_phy_config (jg)% dt_car  ,&
            &                             echam_phy_config (jg)% sd_car  ,&
            &                             echam_phy_config (jg)% ed_car  ,&
            &                             echam_phy_config (jg)% fc_car  )
       !
       CALL eval_echam_phy_config_details(TRIM(cg),                'art' ,&
            &                             echam_phy_config (jg)% dt_art  ,&
            &                             echam_phy_config (jg)% sd_art  ,&
            &                             echam_phy_config (jg)% ed_art  ,&
            &                             echam_phy_config (jg)% fc_art  )
       !
    END DO
    !
  CONTAINS
    !
    SUBROUTINE eval_echam_phy_config_details(cg,process ,&
         &                                   config_dt  ,&
         &                                   config_sd  ,&
         &                                   config_ed  ,&
         &                                   config_fc  )
      !
      CHARACTER(LEN=*),PARAMETER  :: method_name ='eval_echam_phy_config_details'
      !
      ! grid and name of evaluated configuration
      CHARACTER(len=*)                    , INTENT(in)    :: cg
      CHARACTER(len=*)                    , INTENT(in)    :: process
      !
      ! sd, ed and tc arguments are empty strings or 'P...'  strings
      CHARACTER(len=max_timedelta_str_len), INTENT(inout) :: config_dt
      CHARACTER(len=max_datetime_str_len ), INTENT(inout) :: config_sd
      CHARACTER(len=max_datetime_str_len ), INTENT(inout) :: config_ed
      !
      ! forcing control
      INTEGER                             , INTENT(in)    :: config_fc
      !
      ! mtime time control (TC) variables
      TYPE(timedelta), POINTER :: tc_dt
      !
      ! 1. if dt='' or dt contains only blanks, then use dt='PT0S',
      !    because MTIME cannot digest empty strings
      !
      IF (TRIM(config_dt)=='') THEN
         config_dt='PT0S'
      END IF
      !
      !
      ! 2. if dt<0 then stop
      !
      tc_dt => newTimedelta (config_dt)
      IF (tc_dt < dt_zero) THEN
         CALL finish(method_name,'negative echam_phy_config('//TRIM(cg)//')% dt_'//TRIM(process)//' is not allowed')
      END IF
      !
      !
      ! 3. if dt is zero in any format, then set 'PT0S'
      !
      IF (tc_dt == dt_zero) THEN
         config_dt = 'PT0S'
      END IF
      !
      !
      ! 4. if dt>0 check start and end dates
      !
      IF (tc_dt > dt_zero) THEN
         !
         ! if start and end dates are empty strings or contain only blanks
         ! then use the start and stop dates of the experiment
         !
         IF (TRIM(config_sd) == '') config_sd = experimentStartDate
         IF (TRIM(config_ed) == '') config_ed = experimentStopDate
         !
      END IF
      !
      !
      ! 5. fc must be 0, 1 or 2
      !
      SELECT CASE(config_fc)
      CASE(0,1,2)
         ! OK
      CASE DEFAULT
         ! not allowed
         CALL finish(method_name,'echam_phy_config('//TRIM(cg)//')% fc_'//TRIM(process)//' must be 0, 1 or 2')
      END SELECT
      !
      CALL deallocateTimeDelta(tc_dt)
      !
    END SUBROUTINE eval_echam_phy_config_details
    !
  END SUBROUTINE eval_echam_phy_config

  !----

  !>
  !! Evaluate the configuration state
  !!
  SUBROUTINE eval_echam_phy_tc
    !
    INTEGER                 :: jg
    CHARACTER(LEN=2)        :: cg
    !
    ! ECHAM physics timecontrol
    ! -------------------------
    !
    ! mtime events
    !
    echam_phy_events        =  addEventGroup("echam_phy_events_group")
    echam_phy_event_group   => getEventGroup( echam_phy_events )
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL eval_echam_phy_tc_details(cg,                     'rad' ,&
            &                         echam_phy_config(jg)% dt_rad  ,&
            &                         echam_phy_config(jg)% sd_rad  ,&
            &                         echam_phy_config(jg)% ed_rad  ,&
            &                         echam_phy_tc    (jg)% dt_rad  ,&
            &                         echam_phy_tc    (jg)% sd_rad  ,&
            &                         echam_phy_tc    (jg)% ed_rad  ,&
            &                         echam_phy_tc    (jg)% ev_rad  )
       !
       CALL eval_echam_phy_tc_details(cg,                     'vdf' ,&
            &                         echam_phy_config(jg)% dt_vdf  ,&
            &                         echam_phy_config(jg)% sd_vdf  ,&
            &                         echam_phy_config(jg)% ed_vdf  ,&
            &                         echam_phy_tc    (jg)% dt_vdf  ,&
            &                         echam_phy_tc    (jg)% sd_vdf  ,&
            &                         echam_phy_tc    (jg)% ed_vdf  ,&
            &                         echam_phy_tc    (jg)% ev_vdf  )
       !
       CALL eval_echam_phy_tc_details(cg,                     'cnv' ,&
            &                         echam_phy_config(jg)% dt_cnv  ,&
            &                         echam_phy_config(jg)% sd_cnv  ,&
            &                         echam_phy_config(jg)% ed_cnv  ,&
            &                         echam_phy_tc    (jg)% dt_cnv  ,&
            &                         echam_phy_tc    (jg)% sd_cnv  ,&
            &                         echam_phy_tc    (jg)% ed_cnv  ,&
            &                         echam_phy_tc    (jg)% ev_cnv  )
       !
       CALL eval_echam_phy_tc_details(cg,                     'cld' ,&
            &                         echam_phy_config(jg)% dt_cld  ,&
            &                         echam_phy_config(jg)% sd_cld  ,&
            &                         echam_phy_config(jg)% ed_cld  ,&
            &                         echam_phy_tc    (jg)% dt_cld  ,&
            &                         echam_phy_tc    (jg)% sd_cld  ,&
            &                         echam_phy_tc    (jg)% ed_cld  ,&
            &                         echam_phy_tc    (jg)% ev_cld  )
       !
       CALL eval_echam_phy_tc_details(cg,                     'gwd' ,&
            &                         echam_phy_config(jg)% dt_gwd  ,&
            &                         echam_phy_config(jg)% sd_gwd  ,&
            &                         echam_phy_config(jg)% ed_gwd  ,&
            &                         echam_phy_tc    (jg)% dt_gwd  ,&
            &                         echam_phy_tc    (jg)% sd_gwd  ,&
            &                         echam_phy_tc    (jg)% ed_gwd  ,&
            &                         echam_phy_tc    (jg)% ev_gwd  )
       !
       CALL eval_echam_phy_tc_details(cg,                     'sso' ,&
            &                         echam_phy_config(jg)% dt_sso  ,&
            &                         echam_phy_config(jg)% sd_sso  ,&
            &                         echam_phy_config(jg)% ed_sso  ,&
            &                         echam_phy_tc    (jg)% dt_sso  ,&
            &                         echam_phy_tc    (jg)% sd_sso  ,&
            &                         echam_phy_tc    (jg)% ed_sso  ,&
            &                         echam_phy_tc    (jg)% ev_sso  )
       !
       CALL eval_echam_phy_tc_details(cg,                     'mox' ,&
            &                         echam_phy_config(jg)% dt_mox  ,&
            &                         echam_phy_config(jg)% sd_mox  ,&
            &                         echam_phy_config(jg)% ed_mox  ,&
            &                         echam_phy_tc    (jg)% dt_mox  ,&
            &                         echam_phy_tc    (jg)% sd_mox  ,&
            &                         echam_phy_tc    (jg)% ed_mox  ,&
            &                         echam_phy_tc    (jg)% ev_mox  )
       !
       CALL eval_echam_phy_tc_details(cg,                     'car' ,&
            &                         echam_phy_config(jg)% dt_car  ,&
            &                         echam_phy_config(jg)% sd_car  ,&
            &                         echam_phy_config(jg)% ed_car  ,&
            &                         echam_phy_tc    (jg)% dt_car  ,&
            &                         echam_phy_tc    (jg)% sd_car  ,&
            &                         echam_phy_tc    (jg)% ed_car  ,&
            &                         echam_phy_tc    (jg)% ev_car  )
       !
       CALL eval_echam_phy_tc_details(cg,                     'art' ,&
            &                         echam_phy_config(jg)% dt_art  ,&
            &                         echam_phy_config(jg)% sd_art  ,&
            &                         echam_phy_config(jg)% ed_art  ,&
            &                         echam_phy_tc    (jg)% dt_art  ,&
            &                         echam_phy_tc    (jg)% sd_art  ,&
            &                         echam_phy_tc    (jg)% ed_art  ,&
            &                         echam_phy_tc    (jg)% ev_art  )
       !
    END DO
    !
  CONTAINS
    !
    SUBROUTINE eval_echam_phy_tc_details(cg, process ,&
         &                               config_dt   ,&
         &                               config_sd   ,&
         &                               config_ed   ,&
         &                               tc_dt       ,&
         &                               tc_sd       ,&
         &                               tc_ed       ,&
         &                               tc_ev       )
      !
      ! grid and name of evaluated configuration
      CHARACTER(len=*)                    , INTENT(in)    :: cg
      CHARACTER(len=*)                    , INTENT(in)    :: process
      !
      ! configuration strings
      CHARACTER(len=max_timedelta_str_len), INTENT(in) :: config_dt
      CHARACTER(len=max_datetime_str_len ), INTENT(in) :: config_sd
      CHARACTER(len=max_datetime_str_len ), INTENT(in) :: config_ed
      !
      ! mtime time control (TC) variables
      TYPE(timedelta), POINTER :: tc_dt
      TYPE(datetime ), POINTER :: tc_sd
      TYPE(datetime ), POINTER :: tc_ed
      TYPE(event    ), POINTER :: tc_ev
      !
      LOGICAL :: lret
      !
      tc_dt => newTimedelta (config_dt)
      IF (tc_dt > dt_zero) THEN
         tc_sd => newDatetime  (config_sd)
         tc_ed => newDatetime  (config_ed)
         tc_ev => newEvent(process//'_d'//cg, &
              &            tc_sd,             & ! <- start date as reference date!
              &            tc_sd,             &
              &            tc_ed,             &
              &            tc_dt)
         lret = addEventToEventGroup(tc_ev, echam_phy_event_group)
      END IF
      !
    END SUBROUTINE eval_echam_phy_tc_details
    !
  END SUBROUTINE eval_echam_phy_tc

  !----

  !>
  !! Print out the user controlled configuration state and the derived logicals and time controls
  !!
  SUBROUTINE print_echam_phy_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','ECHAM physics configuration')
    CALL message    ('','===========================')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL message    ('','User controlled parameters')
       CALL message    ('','..........................')
       CALL message    ('','')
       CALL message    ('','dynamics physics coupling')
       CALL print_value('    echam_phy_config('//TRIM(cg)//')% ldcphycpl  ',echam_phy_config(jg)% ldcphycpl  )
       CALL message    ('','')
       CALL message    ('','parameterization coupling: .FALSE.: no updates of physics state, .TRUE.: updates phyiscs state')
       CALL print_value('    echam_phy_config('//TRIM(cg)//')% lparamcpl  ',echam_phy_config(jg)% lparamcpl  )
       CALL message    ('','')
       CALL message    ('','tracer reference air mass: .FALSE.: moist air mass, .TRUE.: dry air mass')
       CALL print_value('    echam_phy_config('//TRIM(cg)//')% ldrymoist  ',echam_phy_config(jg)% ldrymoist  )
       CALL message    ('','')
       CALL message    ('','time control parameters')
       CALL message    ('','')
       CALL print_echam_phy_config_details(cg,                     'rad' ,&
            &                              echam_phy_config(jg)% dt_rad  ,&
            &                              echam_phy_config(jg)% sd_rad  ,&
            &                              echam_phy_config(jg)% ed_rad  ,&
            &                              echam_phy_config(jg)% fc_rht  )
       !
       CALL print_echam_phy_config_details(cg,                     'vdf' ,&
            &                              echam_phy_config(jg)% dt_vdf  ,&
            &                              echam_phy_config(jg)% sd_vdf  ,&
            &                              echam_phy_config(jg)% ed_vdf  ,&
            &                              echam_phy_config(jg)% fc_vdf  )
       !
       CALL print_echam_phy_config_details(cg,                     'cnv' ,&
            &                              echam_phy_config(jg)% dt_cnv  ,&
            &                              echam_phy_config(jg)% sd_cnv  ,&
            &                              echam_phy_config(jg)% ed_cnv  ,&
            &                              echam_phy_config(jg)% fc_cnv  )
       !
       CALL print_echam_phy_config_details(cg,                     'cld' ,&
            &                              echam_phy_config(jg)% dt_cld  ,&
            &                              echam_phy_config(jg)% sd_cld  ,&
            &                              echam_phy_config(jg)% ed_cld  ,&
            &                              echam_phy_config(jg)% fc_cld  )
       !
       CALL print_echam_phy_config_details(cg,                     'gwd' ,&
            &                              echam_phy_config(jg)% dt_gwd  ,&
            &                              echam_phy_config(jg)% sd_gwd  ,&
            &                              echam_phy_config(jg)% ed_gwd  ,&
            &                              echam_phy_config(jg)% fc_gwd  )
       !
       CALL print_echam_phy_config_details(cg,                     'sso' ,&
            &                              echam_phy_config(jg)% dt_sso  ,&
            &                              echam_phy_config(jg)% sd_sso  ,&
            &                              echam_phy_config(jg)% ed_sso  ,&
            &                              echam_phy_config(jg)% fc_sso  )
       !
       CALL print_echam_phy_config_details(cg,                     'mox' ,&
            &                              echam_phy_config(jg)% dt_mox  ,&
            &                              echam_phy_config(jg)% sd_mox  ,&
            &                              echam_phy_config(jg)% ed_mox  ,&
            &                              echam_phy_config(jg)% fc_mox  )
       !
       CALL print_echam_phy_config_details(cg,                     'car' ,&
            &                              echam_phy_config(jg)% dt_car  ,&
            &                              echam_phy_config(jg)% sd_car  ,&
            &                              echam_phy_config(jg)% ed_car  ,&
            &                              echam_phy_config(jg)% fc_car  )
       !
       CALL print_echam_phy_config_details(cg,                     'art' ,&
            &                              echam_phy_config(jg)% dt_art  ,&
            &                              echam_phy_config(jg)% sd_art  ,&
            &                              echam_phy_config(jg)% ed_art  ,&
            &                              echam_phy_config(jg)% fc_art  )
       !
       CALL message    ('','logical switches')
       CALL print_value('    echam_phy_config('//TRIM(cg)//')% lmlo ',echam_phy_config(jg)% lmlo  )
       CALL print_value('    echam_phy_config('//TRIM(cg)//')% lice ',echam_phy_config(jg)% lice  )
       CALL print_value('    echam_phy_config('//TRIM(cg)//')% ljsb ',echam_phy_config(jg)% ljsb  )
       CALL print_value('    echam_phy_config('//TRIM(cg)//')% llake',echam_phy_config(jg)% llake )
       CALL print_value('    echam_phy_config('//TRIM(cg)//')% lamip',echam_phy_config(jg)% lamip )
       CALL print_value('    echam_phy_config('//TRIM(cg)//')% lcpl_co2_atmoce',echam_phy_config(jg)% lcpl_co2_atmoce)
       !
       CALL message    ('','')
       CALL message    ('','Derived time control')
       CALL message    ('','....................')
       CALL message    ('','')
       !
       CALL print_echam_phy_tc_details(cg,                 'rad' ,&
            &                          echam_phy_tc(jg)% dt_rad  ,&
            &                          echam_phy_tc(jg)% sd_rad  ,&
            &                          echam_phy_tc(jg)% ed_rad  )
       !
       CALL print_echam_phy_tc_details(cg,                 'vdf' ,&
            &                          echam_phy_tc(jg)% dt_vdf  ,&
            &                          echam_phy_tc(jg)% sd_vdf  ,&
            &                          echam_phy_tc(jg)% ed_vdf  )
       !
       CALL print_echam_phy_tc_details(cg,                 'cnv' ,&
            &                          echam_phy_tc(jg)% dt_cnv  ,&
            &                          echam_phy_tc(jg)% sd_cnv  ,&
            &                          echam_phy_tc(jg)% ed_cnv  )
       !
       CALL print_echam_phy_tc_details(cg,                 'cld' ,&
            &                          echam_phy_tc(jg)% dt_cld  ,&
            &                          echam_phy_tc(jg)% sd_cld  ,&
            &                          echam_phy_tc(jg)% ed_cld  )
       !
       CALL print_echam_phy_tc_details(cg,                 'gwd' ,&
            &                          echam_phy_tc(jg)% dt_gwd  ,&
            &                          echam_phy_tc(jg)% sd_gwd  ,&
            &                          echam_phy_tc(jg)% ed_gwd  )
       !
       CALL print_echam_phy_tc_details(cg,                 'sso' ,&
            &                          echam_phy_tc(jg)% dt_sso  ,&
            &                          echam_phy_tc(jg)% sd_sso  ,&
            &                          echam_phy_tc(jg)% ed_sso  )
       !
       CALL print_echam_phy_tc_details(cg,                 'mox' ,&
            &                          echam_phy_tc(jg)% dt_mox  ,&
            &                          echam_phy_tc(jg)% sd_mox  ,&
            &                          echam_phy_tc(jg)% ed_mox  )
       !
       CALL print_echam_phy_tc_details(cg,                 'car' ,&
            &                          echam_phy_tc(jg)% dt_car  ,&
            &                          echam_phy_tc(jg)% sd_car  ,&
            &                          echam_phy_tc(jg)% ed_car  )
       !
       CALL print_echam_phy_tc_details(cg,                 'art' ,&
            &                          echam_phy_tc(jg)% dt_art  ,&
            &                          echam_phy_tc(jg)% sd_art  ,&
            &                          echam_phy_tc(jg)% ed_art  )
       !
       CALL message    ('','')
       CALL message    ('','------------------------------------------------------------------------')
       CALL message    ('','')
       !
    END DO
    !
    CALL message    ('','Events on all domains')
    CALL message    ('','.....................')
    CALL message    ('','')
    CALL printEventGroup(echam_phy_events)
    CALL message    ('','')
    CALL message    ('','========================================================================')
    !
  CONTAINS
    !
    SUBROUTINE print_echam_phy_config_details(cg, process ,&
         &                                    config_dt   ,&
         &                                    config_sd   ,&
         &                                    config_ed   ,&
         &                                    config_fc   )
      !
      ! grid and name of evaluated configuration
      CHARACTER(len=*)                    , INTENT(in) :: cg
      CHARACTER(len=*)                    , INTENT(in) :: process
      !
      ! configuration strings
      CHARACTER(len=max_timedelta_str_len), INTENT(in) :: config_dt
      CHARACTER(len=max_datetime_str_len ), INTENT(in) :: config_sd
      CHARACTER(len=max_datetime_str_len ), INTENT(in) :: config_ed
      !
      ! forcing control
      INTEGER                             , INTENT(in) :: config_fc
      !
      CALL message       ('    echam_phy_config('//cg//')% dt_'//process,config_dt )
      IF (config_dt /= 'PT0S') THEN
         CALL message    ('    echam_phy_config('//cg//')% sd_'//process,config_sd )
         CALL message    ('    echam_phy_config('//cg//')% ed_'//process,config_ed )
         CALL print_value('    echam_phy_config('//cg//')% fc_'//process,config_fc )
         SELECT CASE(config_fc)
         CASE(0)
            CALL message ('',process//' tendencies are diagnostic')
         CASE(1)
            CALL message ('',process//' tendencies are used to update the model state')
         CASE(2)
            CALL message ('',process//' tendencies are used as forcing in the dynamics')
         END SELECT
      END IF
      CALL message   ('','')
      !
    END SUBROUTINE print_echam_phy_config_details
    !
    !
    SUBROUTINE print_echam_phy_tc_details(cg, process ,&
         &                                tc_dt       ,&
         &                                tc_sd       ,&
         &                                tc_ed       )
      !
      ! grid and name of evaluated configuration
      CHARACTER(len=*)                    , INTENT(in)    :: cg
      CHARACTER(len=*)                    , INTENT(in)    :: process
      !
      ! mtime time control (TC) variables
      TYPE(timedelta), POINTER :: tc_dt
      TYPE(datetime ), POINTER :: tc_sd
      TYPE(datetime ), POINTER :: tc_ed
      !
      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_string
      CHARACTER(LEN=MAX_DATETIME_STR_LEN ) :: dt_string
      !
      CALL timedeltaToString  (tc_dt, td_string)
      CALL message   ('    echam_phy_tc('//cg//')% dt_'//process,td_string)
      IF (tc_dt > dt_zero) THEN
         CALL datetimeToString(tc_sd, dt_string)
         CALL message('    echam_phy_tc('//cg//')% sd_'//process,dt_string)
         CALL datetimeToString(tc_ed, dt_string)
         CALL message('    echam_phy_tc('//cg//')% ed_'//process,dt_string)
      END IF
      CALL message   ('','')
      !
    END SUBROUTINE print_echam_phy_tc_details
    !
  END SUBROUTINE print_echam_phy_config

  !----

END MODULE mo_echam_phy_config
