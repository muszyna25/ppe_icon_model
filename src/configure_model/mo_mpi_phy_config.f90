!>
!! Configuration of the MPI physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!! First version by Marco Giorgetta, MPI (2017-04)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_mpi_phy_config

  USE mo_exception     ,ONLY: message, print_value, finish

  USE mtime            ,ONLY: OPERATOR(<), OPERATOR(>), OPERATOR(==),                                                 &
       &                      datetime , newDatetime , datetimeToString , max_datetime_str_len ,                      &
       &                      timedelta, newTimedelta, timedeltaToString, max_timedelta_str_len, deallocateTimeDelta, &
       &                      event    , newEvent    , eventGroup         , addEventToEventGroup
  USE mo_event_manager ,ONLY: addEventGroup, getEventGroup, printEventGroup

  USE mo_grid_config   ,ONLY: n_dom
  USE mo_master_config ,ONLY: experimentStartDate, experimentStopDate

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::         mpi_phy_config   !< user specified configuration parameters
  PUBLIC ::         mpi_phy_tc       !< derived time control (tc) parameters
  PUBLIC ::    init_mpi_phy_config   !< allocate and initialize mpi_phy_config
  PUBLIC ::    eval_mpi_phy_config   !< evaluate mpi_phy_config
  PUBLIC ::    eval_mpi_phy_tc       !< evaluate mpi_phy_tc
  PUBLIC ::   print_mpi_phy_config   !< print out
  PUBLIC :: dealloc_mpi_phy_config   !< deallocate

  PUBLIC ::               dt_zero    !< a zero (0 sec) time interval 'PT0S'

  !>
  !! Configuration type containing parameters and switches for the configuration of the MPI physics package
  !!
  TYPE t_mpi_phy_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! dynamics physics coupling
     INTEGER  :: idcphycpl  !< determines the coupling between the dynamical core and the
     !                      !  phyiscs package
     !                      !  1: dynamics and physics update sequentially
     !                      !  2: dynamics uses physics forcing for updating
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
     ! atmospheric physics
     CHARACTER(len=max_timedelta_str_len) :: dt_rad  !< time  step of LW radiation
     CHARACTER(len=max_datetime_str_len ) :: sd_rad  !< start time of LW radiation
     CHARACTER(len=max_datetime_str_len ) :: ed_rad  !< end   time of LW radiation
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_vdf  !< time  step of vertical diffusion
     CHARACTER(len=max_datetime_str_len ) :: sd_vdf  !< start time of vertical diffusion
     CHARACTER(len=max_datetime_str_len ) :: ed_vdf  !< end   time of vertical diffusion
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_cnv  !< time  step of cumulus convection
     CHARACTER(len=max_datetime_str_len ) :: sd_cnv  !< start time of cumulus convection
     CHARACTER(len=max_datetime_str_len ) :: ed_cnv  !< end   time of cumulus convection
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_cld  !< time  step of cloud microphysics
     CHARACTER(len=max_datetime_str_len ) :: sd_cld  !< start time of cloud microphysics
     CHARACTER(len=max_datetime_str_len ) :: ed_cld  !< end   time of cloud microphysics
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_gwd  !< time  step of atmospheric gravity wave drag
     CHARACTER(len=max_datetime_str_len ) :: sd_gwd  !< start time of atmospheric gravity wave drag
     CHARACTER(len=max_datetime_str_len ) :: ed_gwd  !< end   time of atmospheric gravity wave drag
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_sso  !< time  step of sub grid scale orogr. effects
     CHARACTER(len=max_datetime_str_len ) :: sd_sso  !< start time of sub grid scale orogr. effects
     CHARACTER(len=max_datetime_str_len ) :: ed_sso  !< end   time of sub grid scale orogr. effects
     !
     ! atmospheric chemistry
     CHARACTER(len=max_timedelta_str_len) :: dt_mox  !< time  step of methan oxidation and water vapor photolysis
     CHARACTER(len=max_datetime_str_len ) :: sd_mox  !< start time of methan oxidation and water vapor photolysis
     CHARACTER(len=max_datetime_str_len ) :: ed_mox  !< end   time of methan oxidation and water vapor photolysis
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_car  !< time  step of lin. Cariolle ozone chemistry
     CHARACTER(len=max_datetime_str_len ) :: sd_car  !< start time of lin. Cariolle ozone chemistry
     CHARACTER(len=max_datetime_str_len ) :: ed_car  !< end   time of lin. Cariolle ozone chemistry
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_art  !< time  step of ART chemistry
     CHARACTER(len=max_datetime_str_len ) :: sd_art  !< start time of ART chemistry
     CHARACTER(len=max_datetime_str_len ) :: ed_art  !< end   time of ART chemistry
     !
     ! surface
     LOGICAL                              :: lmlo    !< .true. for mixed layer ocean
     LOGICAL                              :: lice    !< .true. for sea-ice temperature calculation
     LOGICAL                              :: ljsb    !< .true. for calculating the JSBACH land surface
     LOGICAL                              :: llake   !< .true. for using lakes in JSBACH
     LOGICAL                              :: lamip   !< .true. for AMIP simulations
     !
  END TYPE t_mpi_phy_config


  TYPE t_mpi_phy_tc
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
  END TYPE t_mpi_phy_tc

  !>
  !! Configuration/logicals/timecontrol state vectors, for multiple domains/grids.
  !!
  TYPE(t_mpi_phy_config) , ALLOCATABLE :: mpi_phy_config  (:)
  TYPE(t_mpi_phy_tc)     , ALLOCATABLE :: mpi_phy_tc      (:)
  
  !>
  !! Events and event group
  !!
  INTEGER                   :: mpi_phy_events
  TYPE(eventGroup), POINTER :: mpi_phy_event_group

  !>
  !! For convenience
  !!
  TYPE(timedelta) , POINTER :: dt_zero

CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_mpi_phy_config

    dt_zero =>  newTimedelta ('PT0S')

    ! MPI physics configuration
    ! -------------------------

    ALLOCATE(mpi_phy_config(n_dom))

    ! dynamics physics coupling
    mpi_phy_config(:)%idcphycpl = 1
    mpi_phy_config(:)%ldrymoist = .TRUE.
    !
    ! time control parameters
    mpi_phy_config(:)% dt_rad = ''
    mpi_phy_config(:)% sd_rad = ''
    mpi_phy_config(:)% ed_rad = ''
    !
    mpi_phy_config(:)% dt_vdf = ''
    mpi_phy_config(:)% sd_vdf = ''
    mpi_phy_config(:)% ed_vdf = ''
    !
    mpi_phy_config(:)% dt_cnv = ''
    mpi_phy_config(:)% sd_cnv = ''
    mpi_phy_config(:)% ed_cnv = ''
    !
    mpi_phy_config(:)% dt_cld = ''
    mpi_phy_config(:)% sd_cld = ''
    mpi_phy_config(:)% ed_cld = ''
    !
    mpi_phy_config(:)% dt_gwd = ''
    mpi_phy_config(:)% sd_gwd = ''
    mpi_phy_config(:)% ed_gwd = ''
    !
    mpi_phy_config(:)% dt_sso = ''
    mpi_phy_config(:)% sd_sso = ''
    mpi_phy_config(:)% ed_sso = ''
    !
    mpi_phy_config(:)% dt_mox = ''
    mpi_phy_config(:)% sd_mox = ''
    mpi_phy_config(:)% ed_mox = ''
    !
    mpi_phy_config(:)% dt_car = ''
    mpi_phy_config(:)% sd_car = ''
    mpi_phy_config(:)% ed_car = ''
    !
    mpi_phy_config(:)% dt_art = ''
    mpi_phy_config(:)% sd_art = ''
    mpi_phy_config(:)% ed_art = ''
    !
    ! logical switches
    mpi_phy_config(:)% ljsb  = .FALSE.
    mpi_phy_config(:)% llake = .FALSE.
    mpi_phy_config(:)% lamip = .FALSE.
    mpi_phy_config(:)% lmlo  = .FALSE.
    mpi_phy_config(:)% lice  = .FALSE.

  END SUBROUTINE init_mpi_phy_config

  !----

  !>
  !! Check the mpi_phy_config state
  !!
  SUBROUTINE eval_mpi_phy_config

    INTEGER                 :: jg
    CHARACTER(LEN=2)        :: cg

    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL eval_mpi_phy_config_dt_sd_ed(TRIM(cg),              'rad', &
            &                            mpi_phy_config (jg)% dt_rad,  &
            &                            mpi_phy_config (jg)% sd_rad,  &
            &                            mpi_phy_config (jg)% ed_rad)
       !
       CALL eval_mpi_phy_config_dt_sd_ed(TRIM(cg),              'vdf', &
            &                            mpi_phy_config (jg)% dt_vdf,  &
            &                            mpi_phy_config (jg)% sd_vdf,  &
            &                            mpi_phy_config (jg)% ed_vdf)
       !
       CALL eval_mpi_phy_config_dt_sd_ed(TRIM(cg),              'cnv', &
            &                            mpi_phy_config (jg)% dt_cnv,  &
            &                            mpi_phy_config (jg)% sd_cnv,  &
            &                            mpi_phy_config (jg)% ed_cnv)
       !
       CALL eval_mpi_phy_config_dt_sd_ed(TRIM(cg),              'cld', &
            &                            mpi_phy_config (jg)% dt_cld,  &
            &                            mpi_phy_config (jg)% sd_cld,  &
            &                            mpi_phy_config (jg)% ed_cld)
       !
       CALL eval_mpi_phy_config_dt_sd_ed(TRIM(cg),              'gwd', &
            &                            mpi_phy_config (jg)% dt_gwd,  &
            &                            mpi_phy_config (jg)% sd_gwd,  &
            &                            mpi_phy_config (jg)% ed_gwd)
       !
       CALL eval_mpi_phy_config_dt_sd_ed(TRIM(cg),              'sso', &
            &                            mpi_phy_config (jg)% dt_sso,  &
            &                            mpi_phy_config (jg)% sd_sso,  &
            &                            mpi_phy_config (jg)% ed_sso)
       !
       CALL eval_mpi_phy_config_dt_sd_ed(TRIM(cg),              'mox', &
            &                            mpi_phy_config (jg)% dt_mox,  &
            &                            mpi_phy_config (jg)% sd_mox,  &
            &                            mpi_phy_config (jg)% ed_mox)
       !
       CALL eval_mpi_phy_config_dt_sd_ed(TRIM(cg),              'car', &
            &                            mpi_phy_config (jg)% dt_car,  &
            &                            mpi_phy_config (jg)% sd_car,  &
            &                            mpi_phy_config (jg)% ed_car)
       !
       CALL eval_mpi_phy_config_dt_sd_ed(TRIM(cg),              'art', &
            &                            mpi_phy_config (jg)% dt_art,  &
            &                            mpi_phy_config (jg)% sd_art,  &
            &                            mpi_phy_config (jg)% ed_art)
       !
    END DO

  CONTAINS

    SUBROUTINE eval_mpi_phy_config_dt_sd_ed(cg,process,config_dt,config_sd,config_ed)

      CHARACTER(LEN=*),PARAMETER  :: method_name ='eval_mpi_phy_config_dt_sd_ed'

      ! grid and name of evaluated configuration
      CHARACTER(len=*)                    , INTENT(in)    :: cg
      CHARACTER(len=*)                    , INTENT(in)    :: process

      ! sd, ed and tc arguments are empty strings or 'P...'  strings
      CHARACTER(len=max_timedelta_str_len), INTENT(inout) :: config_dt
      CHARACTER(len=max_datetime_str_len ), INTENT(inout) :: config_sd
      CHARACTER(len=max_datetime_str_len ), INTENT(inout) :: config_ed

      ! mtime time control (TC) variables
      TYPE(timedelta), POINTER :: tc_dt

      ! 1. if dt='' or dt contains only blanks, then use dt='PT0S',
      !    because MTIME cannot digest empty strings
      !
      IF (TRIM(config_dt)=='') THEN
         config_dt='PT0S'
      END IF


      ! 2. if dt<0 then stop
      !
      tc_dt => newTimedelta (config_dt)
      IF (tc_dt < dt_zero) THEN
         CALL finish(method_name,'negative mpi_phy_config('//TRIM(cg)//')% dt_'//TRIM(process)//' is not allowed')
      END IF


      ! 3. if dt is zero in any format, then set 'PT0S'
      !
      IF (tc_dt == dt_zero) THEN
         config_dt = 'PT0S'
      END IF


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

      CALL deallocateTimeDelta(tc_dt)

    END SUBROUTINE eval_mpi_phy_config_dt_sd_ed

  END SUBROUTINE eval_mpi_phy_config

  !----

  !>
  !! Evaluate the configuration state
  !!
  SUBROUTINE eval_mpi_phy_tc

    INTEGER                 :: jg
    CHARACTER(LEN=2)        :: cg
    LOGICAL                 :: lret

    ! MPI physics timecontrol
    ! -----------------------
    !
    ! mtime events
    !
    mpi_phy_events        =  addEventGroup("mpi_phy_events_group")
    mpi_phy_event_group   => getEventGroup( mpi_phy_events )

    ALLOCATE(mpi_phy_tc(n_dom))
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL eval_mpi_phy_tc_dt_sd_ed_ev(cg,                   'rad', &
            &                           mpi_phy_config(jg)% dt_rad,  &
            &                           mpi_phy_config(jg)% sd_rad,  &
            &                           mpi_phy_config(jg)% ed_rad,  &
            &                           mpi_phy_tc    (jg)% dt_rad,  &
            &                           mpi_phy_tc    (jg)% sd_rad,  &
            &                           mpi_phy_tc    (jg)% ed_rad,  &
            &                           mpi_phy_tc    (jg)% ev_rad   )
       !
       CALL eval_mpi_phy_tc_dt_sd_ed_ev(cg,                   'vdf', &
            &                           mpi_phy_config(jg)% dt_vdf,  &
            &                           mpi_phy_config(jg)% sd_vdf,  &
            &                           mpi_phy_config(jg)% ed_vdf,  &
            &                           mpi_phy_tc    (jg)% dt_vdf,  &
            &                           mpi_phy_tc    (jg)% sd_vdf,  &
            &                           mpi_phy_tc    (jg)% ed_vdf,  &
            &                           mpi_phy_tc    (jg)% ev_vdf   )
       !
       CALL eval_mpi_phy_tc_dt_sd_ed_ev(cg,                   'cnv', &
            &                           mpi_phy_config(jg)% dt_cnv,  &
            &                           mpi_phy_config(jg)% sd_cnv,  &
            &                           mpi_phy_config(jg)% ed_cnv,  &
            &                           mpi_phy_tc    (jg)% dt_cnv,  &
            &                           mpi_phy_tc    (jg)% sd_cnv,  &
            &                           mpi_phy_tc    (jg)% ed_cnv,  &
            &                           mpi_phy_tc    (jg)% ev_cnv   )
       !
       CALL eval_mpi_phy_tc_dt_sd_ed_ev(cg,                   'cld', &
            &                           mpi_phy_config(jg)% dt_cld,  &
            &                           mpi_phy_config(jg)% sd_cld,  &
            &                           mpi_phy_config(jg)% ed_cld,  &
            &                           mpi_phy_tc    (jg)% dt_cld,  &
            &                           mpi_phy_tc    (jg)% sd_cld,  &
            &                           mpi_phy_tc    (jg)% ed_cld,  &
            &                           mpi_phy_tc    (jg)% ev_cld   )
       !
       CALL eval_mpi_phy_tc_dt_sd_ed_ev(cg,                   'gwd', &
            &                           mpi_phy_config(jg)% dt_gwd,  &
            &                           mpi_phy_config(jg)% sd_gwd,  &
            &                           mpi_phy_config(jg)% ed_gwd,  &
            &                           mpi_phy_tc    (jg)% dt_gwd,  &
            &                           mpi_phy_tc    (jg)% sd_gwd,  &
            &                           mpi_phy_tc    (jg)% ed_gwd,  &
            &                           mpi_phy_tc    (jg)% ev_gwd   )
       !
       CALL eval_mpi_phy_tc_dt_sd_ed_ev(cg,                   'sso', &
            &                           mpi_phy_config(jg)% dt_sso,  &
            &                           mpi_phy_config(jg)% sd_sso,  &
            &                           mpi_phy_config(jg)% ed_sso,  &
            &                           mpi_phy_tc    (jg)% dt_sso,  &
            &                           mpi_phy_tc    (jg)% sd_sso,  &
            &                           mpi_phy_tc    (jg)% ed_sso,  &
            &                           mpi_phy_tc    (jg)% ev_sso   )
       !
       CALL eval_mpi_phy_tc_dt_sd_ed_ev(cg,                   'mox', &
            &                           mpi_phy_config(jg)% dt_mox,  &
            &                           mpi_phy_config(jg)% sd_mox,  &
            &                           mpi_phy_config(jg)% ed_mox,  &
            &                           mpi_phy_tc    (jg)% dt_mox,  &
            &                           mpi_phy_tc    (jg)% sd_mox,  &
            &                           mpi_phy_tc    (jg)% ed_mox,  &
            &                           mpi_phy_tc    (jg)% ev_mox   )
       !
       CALL eval_mpi_phy_tc_dt_sd_ed_ev(cg,                   'car', &
            &                           mpi_phy_config(jg)% dt_car,  &
            &                           mpi_phy_config(jg)% sd_car,  &
            &                           mpi_phy_config(jg)% ed_car,  &
            &                           mpi_phy_tc    (jg)% dt_car,  &
            &                           mpi_phy_tc    (jg)% sd_car,  &
            &                           mpi_phy_tc    (jg)% ed_car,  &
            &                           mpi_phy_tc    (jg)% ev_car   )
       !
       CALL eval_mpi_phy_tc_dt_sd_ed_ev(cg,                   'art', &
            &                           mpi_phy_config(jg)% dt_art,  &
            &                           mpi_phy_config(jg)% sd_art,  &
            &                           mpi_phy_config(jg)% ed_art,  &
            &                           mpi_phy_tc    (jg)% dt_art,  &
            &                           mpi_phy_tc    (jg)% sd_art,  &
            &                           mpi_phy_tc    (jg)% ed_art,  &
            &                           mpi_phy_tc    (jg)% ev_art   )
       !
    END DO

  CONTAINS

    SUBROUTINE eval_mpi_phy_tc_dt_sd_ed_ev(cg, process, &
         &                                 config_dt,   &
         &                                 config_sd,   &
         &                                 config_ed,   &
         &                                 tc_dt,       &
         &                                 tc_sd,       &
         &                                 tc_ed,       &
         &                                 tc_ev        )

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

      LOGICAL :: lret

      tc_dt => newTimedelta (config_dt)
      IF (tc_dt > dt_zero) THEN
         tc_sd => newDatetime  (config_sd)
         tc_ed => newDatetime  (config_ed)
         tc_ev => newEvent(process//'_d'//cg, &
              &            tc_sd,             & ! <- start date as reference date!
              &            tc_sd,             &
              &            tc_ed,             &
              &            tc_dt)
         lret = addEventToEventGroup(tc_ev, mpi_phy_event_group)
      END IF

    END SUBROUTINE eval_mpi_phy_tc_dt_sd_ed_ev

  END SUBROUTINE eval_mpi_phy_tc

  !----

  !>
  !! Print out the user controlled configuration state and the derived logicals and time controls
  !!
  SUBROUTINE print_mpi_phy_config

    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg

    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','MPI physics configuration')
    CALL message    ('','=========================')
    CALL message    ('','')

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
!!$       CALL message    ('','dynamics physics coupling')
!!$       CALL print_value('    mpi_phy_config('//TRIM(cg)//')% idcphycpl  ',mpi_phy_config(jg)% idcphycpl  )
!!$       CALL print_value('    mpi_phy_config('//TRIM(cg)//')% ldrymoist  ',mpi_phy_config(jg)% ldrymoist  )
!!$       CALL message    ('','')
       CALL message    ('','time control parameters')
       CALL message    ('','')
       CALL print_mpi_phy_config_dt_sd_ed(cg,                   'rad', &
            &                             mpi_phy_config(jg)% dt_rad,  &
            &                             mpi_phy_config(jg)% sd_rad,  &
            &                             mpi_phy_config(jg)% ed_rad   )
       !
       CALL print_mpi_phy_config_dt_sd_ed(cg,                   'vdf', &
            &                             mpi_phy_config(jg)% dt_vdf,  &
            &                             mpi_phy_config(jg)% sd_vdf,  &
            &                             mpi_phy_config(jg)% ed_vdf   )
       !
       CALL print_mpi_phy_config_dt_sd_ed(cg,                   'cnv', &
            &                             mpi_phy_config(jg)% dt_cnv,  &
            &                             mpi_phy_config(jg)% sd_cnv,  &
            &                             mpi_phy_config(jg)% ed_cnv   )
       !
       CALL print_mpi_phy_config_dt_sd_ed(cg,                   'cld', &
            &                             mpi_phy_config(jg)% dt_cld,  &
            &                             mpi_phy_config(jg)% sd_cld,  &
            &                             mpi_phy_config(jg)% ed_cld   )
       !
       CALL print_mpi_phy_config_dt_sd_ed(cg,                   'gwd', &
            &                             mpi_phy_config(jg)% dt_gwd,  &
            &                             mpi_phy_config(jg)% sd_gwd,  &
            &                             mpi_phy_config(jg)% ed_gwd   )
       !
       CALL print_mpi_phy_config_dt_sd_ed(cg,                   'sso', &
            &                             mpi_phy_config(jg)% dt_sso,  &
            &                             mpi_phy_config(jg)% sd_sso,  &
            &                             mpi_phy_config(jg)% ed_sso   )
       !
       CALL print_mpi_phy_config_dt_sd_ed(cg,                   'mox', &
            &                             mpi_phy_config(jg)% dt_mox,  &
            &                             mpi_phy_config(jg)% sd_mox,  &
            &                             mpi_phy_config(jg)% ed_mox   )
       !
       CALL print_mpi_phy_config_dt_sd_ed(cg,                   'car', &
            &                             mpi_phy_config(jg)% dt_car,  &
            &                             mpi_phy_config(jg)% sd_car,  &
            &                             mpi_phy_config(jg)% ed_car   )
       !
       CALL print_mpi_phy_config_dt_sd_ed(cg,                   'art', &
            &                             mpi_phy_config(jg)% dt_art,  &
            &                             mpi_phy_config(jg)% sd_art,  &
            &                             mpi_phy_config(jg)% ed_art   )
       !
       CALL message    ('','logical switches')
       CALL print_value('    mpi_phy_config('//TRIM(cg)//')% lmlo ',mpi_phy_config(jg)% lmlo         )
       CALL print_value('    mpi_phy_config('//TRIM(cg)//')% lice ',mpi_phy_config(jg)% lice         )
       CALL print_value('    mpi_phy_config('//TRIM(cg)//')% ljsb ',mpi_phy_config(jg)% ljsb         )
       CALL print_value('    mpi_phy_config('//TRIM(cg)//')% llake',mpi_phy_config(jg)% llake        )
       CALL print_value('    mpi_phy_config('//TRIM(cg)//')% lamip',mpi_phy_config(jg)% lamip        )
       CALL message    ('','')
       CALL message    ('','')
       !
       CALL message    ('','')
       CALL message    ('','Derived time control')
       CALL message    ('','....................')
       CALL message    ('','')
       !
       CALL print_mpi_phy_tc_dt_sd_ed(cg,               'rad', &
            &                         mpi_phy_tc(jg)% dt_rad,  &
            &                         mpi_phy_tc(jg)% sd_rad,  &
            &                         mpi_phy_tc(jg)% ed_rad   )
       !
       CALL print_mpi_phy_tc_dt_sd_ed(cg,               'vdf', &
            &                         mpi_phy_tc(jg)% dt_vdf,  &
            &                         mpi_phy_tc(jg)% sd_vdf,  &
            &                         mpi_phy_tc(jg)% ed_vdf   )
       !
       CALL print_mpi_phy_tc_dt_sd_ed(cg,               'cnv', &
            &                         mpi_phy_tc(jg)% dt_cnv,  &
            &                         mpi_phy_tc(jg)% sd_cnv,  &
            &                         mpi_phy_tc(jg)% ed_cnv   )
       !
       CALL print_mpi_phy_tc_dt_sd_ed(cg,               'cld', &
            &                         mpi_phy_tc(jg)% dt_cld,  &
            &                         mpi_phy_tc(jg)% sd_cld,  &
            &                         mpi_phy_tc(jg)% ed_cld   )
       !
       CALL print_mpi_phy_tc_dt_sd_ed(cg,               'gwd', &
            &                         mpi_phy_tc(jg)% dt_gwd,  &
            &                         mpi_phy_tc(jg)% sd_gwd,  &
            &                         mpi_phy_tc(jg)% ed_gwd   )
       !
       CALL print_mpi_phy_tc_dt_sd_ed(cg,               'sso', &
            &                         mpi_phy_tc(jg)% dt_sso,  &
            &                         mpi_phy_tc(jg)% sd_sso,  &
            &                         mpi_phy_tc(jg)% ed_sso   )
       !
       CALL print_mpi_phy_tc_dt_sd_ed(cg,               'mox', &
            &                         mpi_phy_tc(jg)% dt_mox,  &
            &                         mpi_phy_tc(jg)% sd_mox,  &
            &                         mpi_phy_tc(jg)% ed_mox   )
       !
       CALL print_mpi_phy_tc_dt_sd_ed(cg,               'car', &
            &                         mpi_phy_tc(jg)% dt_car,  &
            &                         mpi_phy_tc(jg)% sd_car,  &
            &                         mpi_phy_tc(jg)% ed_car   )
       !
       CALL print_mpi_phy_tc_dt_sd_ed(cg,               'art', &
            &                         mpi_phy_tc(jg)% dt_art,  &
            &                         mpi_phy_tc(jg)% sd_art,  &
            &                         mpi_phy_tc(jg)% ed_art   )
       !
       CALL message    ('','')
       CALL message    ('','------------------------------------------------------------------------')
       CALL message    ('','')

    END DO

    CALL message    ('','Events on all domains')
    CALL message    ('','.....................')
    CALL message    ('','')
    CALL printEventGroup(mpi_phy_events)
    CALL message    ('','')
    CALL message    ('','========================================================================')

  CONTAINS

    SUBROUTINE print_mpi_phy_config_dt_sd_ed(cg, process, &
         &                                   config_dt,   &
         &                                   config_sd,   &
         &                                   config_ed    )

      ! grid and name of evaluated configuration
      CHARACTER(len=*), INTENT(in) :: cg
      CHARACTER(len=*), INTENT(in) :: process
      !
      ! configuration strings
      CHARACTER(len=max_timedelta_str_len), INTENT(in) :: config_dt
      CHARACTER(len=max_datetime_str_len ), INTENT(in) :: config_sd
      CHARACTER(len=max_datetime_str_len ), INTENT(in) :: config_ed

      CALL message   ('    mpi_phy_config('//cg//')% dt_'//process,config_dt )
      IF (config_dt /= 'PT0S') THEN
         CALL message('    mpi_phy_config('//cg//')% sd_'//process,config_sd )
         CALL message('    mpi_phy_config('//cg//')% ed_'//process,config_ed )
      END IF
      CALL message   ('','')

    END SUBROUTINE print_mpi_phy_config_dt_sd_ed


    SUBROUTINE print_mpi_phy_tc_dt_sd_ed(cg, process, &
         &                               tc_dt,       &
         &                               tc_sd,       &
         &                               tc_ed        )

      ! grid and name of evaluated configuration
      CHARACTER(len=*)                    , INTENT(in)    :: cg
      CHARACTER(len=*)                    , INTENT(in)    :: process
      !
      ! mtime time control (TC) variables
      TYPE(timedelta), POINTER :: tc_dt
      TYPE(datetime ), POINTER :: tc_sd
      TYPE(datetime ), POINTER :: tc_ed

      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_string
      CHARACTER(LEN=MAX_DATETIME_STR_LEN ) :: dt_string

      CALL timedeltaToString  (tc_dt, td_string)
      CALL message   ('    mpi_phy_tc('//cg//')% dt_'//process,td_string)
      IF (tc_dt > dt_zero) THEN
         CALL datetimeToString(tc_sd, dt_string)
         CALL message('    mpi_phy_tc('//cg//')% sd_'//process,dt_string)
         CALL datetimeToString(tc_ed, dt_string)
         CALL message('    mpi_phy_tc('//cg//')% ed_'//process,dt_string)
      END IF
      CALL message   ('','')

    END SUBROUTINE print_mpi_phy_tc_dt_sd_ed

  END SUBROUTINE print_mpi_phy_config

  !----

  !>
  !! Deallocate the mpi_phy_xyz states
  !!
  SUBROUTINE dealloc_mpi_phy_config

    CHARACTER(LEN=*),PARAMETER  :: method_name ='mo_mpi_phy_config:dealloc_mpi_phy_config'

    DEALLOCATE(mpi_phy_config)
    DEALLOCATE(mpi_phy_tc)

    CALL message    ('','')
    CALL message    ('','------------------------------------------------------------------------')
    CALL message    ('','')
    CALL message    (method_name,'MPI physics configuration states deallocated')
    CALL message    ('','')
    CALL message    ('','------------------------------------------------------------------------')
    CALL message    ('','')

  END SUBROUTINE dealloc_mpi_phy_config

  !----

END MODULE mo_mpi_phy_config
