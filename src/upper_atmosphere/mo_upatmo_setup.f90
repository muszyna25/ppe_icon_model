!>
!! Construction and destruction of the upper atmosphere.
!!
!! @author Sebastian Borchert, DWD
!!
!! @par Revision History
!! Initial revision by Sebastian Borchert, DWD, 2016-09-06
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_upatmo_setup

  USE mo_exception,               ONLY: message, finish
  USE mo_impl_constants,          ONLY: MAX_CHAR_LENGTH
  USE mo_upatmo_impl_const,       ONLY: imsg_thr, itmr_thr, iUpatmoPrcStat
  USE mo_model_domain,            ONLY: t_patch
  USE mo_master_config,           ONLY: getModelBaseDir, isRestart
  USE mo_grid_config,             ONLY: n_dom_start, n_dom, &
    &                                   start_time, end_time
  USE mo_dynamics_config,         ONLY: ldeepatmo
  USE mo_initicon_config,         ONLY: init_mode
  USE mo_run_config,              ONLY: dtime, iforcing, msg_level
  USE mo_nonhydrostatic_config,   ONLY: ndyn_substeps
  USE mo_sleve_config,            ONLY: flat_height
  USE mo_time_config,             ONLY: time_config
  USE mo_echam_rad_config,        ONLY: echam_rad_config
  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config
  USE mo_name_list_output_config, ONLY: first_output_name_list
  USE mo_parallel_config,         ONLY: nproma
  USE mo_vertical_coord_table,    ONLY: vct_a
  USE mo_timer,                   ONLY: timers_level, timer_start, timer_stop,   & 
    &                                   timer_upatmo_constr, timer_upatmo_destr, &
    &                                   timer_upatmo
  USE mo_upatmo_config,           ONLY: upatmo_config, upatmo_phy_config, &
    &                                   configure_upatmo, destruct_upatmo
  USE mo_upatmo_state,            ONLY: construct_upatmo_state, &
    &                                   destruct_upatmo_state
  USE mo_upatmo_phy_setup,        ONLY: finalize_upatmo_phy_nwp

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: upatmo_initialize
  PUBLIC :: upatmo_finalize

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_upatmo_setup'
  
CONTAINS

  !>
  !! Initialize upper atmosphere.
  !!
  SUBROUTINE upatmo_initialize( p_patch )

    ! In/out variables
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(n_dom_start:)

    ! Local variables
    LOGICAL  :: lmessage, ltimer, lrestart
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: model_base_dir
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':upatmo_initialize'
    
    !---------------------------------------------------------

    lmessage = msg_level >= imsg_thr%high
    ltimer   = timers_level > itmr_thr%med

    IF (ltimer) THEN 
      CALL timer_start(timer_upatmo)
      CALL timer_start(timer_upatmo_constr)
    ENDIF

    ! 'upatmo_config' should have been allocated in 'src/namelists/mo_upatmo_nml: check_upatmo', 
    ! which should have been called before this subroutine
    IF (.NOT. ALLOCATED(upatmo_config)) THEN
      CALL finish(TRIM(routine), "Check calling sequence: upatmo_config is not allocated.")
    ENDIF

    IF (lmessage) CALL message(TRIM(routine), 'Initialization of upper atmosphere started')

    !---------------------------------------------------------------------
    !          Construct the upper-atmosphere configuration type
    !---------------------------------------------------------------------

    model_base_dir = getModelBaseDir()
    lrestart       = isRestart()

    CALL configure_upatmo( n_dom_start            = n_dom_start,                       & !in
      &                    n_dom                  = n_dom,                             & !in 
      &                    p_patch                = p_patch(n_dom_start:),             & !in
      &                    lrestart               = lrestart,                          & !in
      &                    ldeepatmo              = ldeepatmo,                         & !in
      &                    lupatmo_phy            = atm_phy_nwp_config(:)%lupatmo_phy, & !in
      &                    init_mode              = init_mode,                         & !in
      &                    iforcing               = iforcing,                          & !in
      &                    tc_exp_startdate       = time_config%tc_exp_startdate,      & !in
      &                    tc_exp_stopdate        = time_config%tc_exp_stopdate,       & !in
      &                    start_time             = start_time(:),                     & !in 
      &                    end_time               = end_time(:),                       & !in
      &                    dtime                  = dtime,                             & !in
      &                    dt_rad_nwp             = atm_phy_nwp_config(:)%dt_rad,      & !in
      &                    ndyn_substeps          = ndyn_substeps,                     & !in
      &                    flat_height            = flat_height,                       & !in
      &                    l_orbvsop87            = echam_rad_config(:)%l_orbvsop87,   & !in
      &                    cecc                   = echam_rad_config(:)%cecc,          & !in
      &                    cobld                  = echam_rad_config(:)%cobld,         & !in
      &                    clonp                  = echam_rad_config(:)%clonp,         & !in
      &                    lyr_perp               = echam_rad_config(:)%lyr_perp,      & !in
      &                    yr_perp                = echam_rad_config(:)%yr_perp,       & !in
      &                    model_base_dir         = model_base_dir,                    & !in
      &                    msg_level              = msg_level,                         & !in
      &                    timers_level           = timers_level,                      & !in
      &                    first_output_name_list = first_output_name_list,            & !in
      &                    vct_a                  = vct_a                              ) !(opt)in

    !---------------------------------------------------------------------
    !              Construct the upper-atmosphere data types
    !---------------------------------------------------------------------
    
    CALL construct_upatmo_state( n_dom             = n_dom,                 & !in
      &                          nproma            = nproma,                & !in
      &                          p_patch           = p_patch(1:),           & !in
      &                          upatmo_config     = upatmo_config(1:),     & !in
      &                          upatmo_phy_config = upatmo_phy_config(1:), & !in
      &                          vct_a             = vct_a                  ) !(opt)in

    !---------------------------------------------------------------------
    !       Construct the upper-atmosphere physics parameterizations
    !---------------------------------------------------------------------

    ! * NWP forcing:
    ! Please see 'src/atm_phy_nwp/mo_nwp_phy_init: init_nwp_phy' 
    ! for call of 'src/upper_atmosphere/mo_upatmo_phy_setup: init_upatmo_phy_nwp'
    
    IF (lmessage) CALL message(TRIM(routine), 'Initialization of upper atmosphere finished')

    IF (ltimer) THEN 
      CALL timer_stop(timer_upatmo_constr)
      CALL timer_stop(timer_upatmo)
    ENDIF

  END SUBROUTINE upatmo_initialize

  !====================================================================================

  !>
  !! Finalize upper atmosphere.
  !!
  SUBROUTINE upatmo_finalize( p_patch ) 

    ! In/out variables
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(n_dom_start:)

    ! Local variables
    INTEGER  :: jg
    LOGICAL  :: lmessage, ltimer
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':upatmo_finalize'
    
    !---------------------------------------------------------

    lmessage = msg_level >= imsg_thr%high
    ltimer   = timers_level > itmr_thr%med

    IF (ltimer) THEN 
      CALL timer_start(timer_upatmo)
      CALL timer_start(timer_upatmo_destr)
    ENDIF

    ! 'upatmo_config' should have been allocated
    IF (.NOT. ALLOCATED(upatmo_config)) THEN
      CALL finish(TRIM(routine), "Error: upatmo_config is not allocated.")
    ENDIF

    IF (lmessage) CALL message(TRIM(routine), 'Finalization of upper atmosphere started')

    !---------------------------------------------------------------------
    !      Finalize the upper-atmosphere physics parameterizations
    !---------------------------------------------------------------------
    
    ! This is required for NWP forcing only. 
    ! For ECHAM forcing, the following will likely be done in 
    ! 'src/atm_phy_echam/mo_echam_phy_cleanup: cleanup_echam_phy'
    DO jg = 1, n_dom
      IF (upatmo_config( jg )%nwp_phy%l_phy_stat( iUpatmoPrcStat%enabled )) THEN
        CALL finalize_upatmo_phy_nwp( p_patch( jg ) ) !in
      ENDIF
    ENDDO  !jg

    !---------------------------------------------------------------------
    !             Destruct the upper-atmosphere data types
    !---------------------------------------------------------------------

    CALL destruct_upatmo_state( n_dom         = n_dom,            & !in
      &                         upatmo_config = upatmo_config(1:) ) !in

    !---------------------------------------------------------------------
    !          Destruct the upper-atmosphere configuration type
    !---------------------------------------------------------------------

    ! After the following call, 'upatmo_config' cannot be used anymore!
    CALL destruct_upatmo( n_dom_start = n_dom_start, & !in
      &                   n_dom       = n_dom        ) !in

    IF (lmessage) CALL message(TRIM(routine), 'Finalization of upper atmosphere finished')

    IF (ltimer) THEN 
      CALL timer_stop(timer_upatmo_destr)
      CALL timer_stop(timer_upatmo)
    ENDIF

  END SUBROUTINE upatmo_finalize

END MODULE mo_upatmo_setup
