!>
!! This module configures:
!! - Upper-atmosphere physics
!! - Deep-atmosphere dynamics
!! - Upper-atmosphere extrapolation
!!
!! @author Guidi Zhou, MPI-M, 2016-03-03
!!         Sebastian Borchert, DWD, 2016-03-03
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_upatmo_config

  USE mo_kind,                     ONLY: wp
  USE mo_exception,                ONLY: message, message_text, finish
  USE mo_impl_constants,           ONLY: max_dom, MAX_CHAR_LENGTH,   &
    &                                    MODE_IFSANA, MODE_COMBINED, &
    &                                    inwp, iecham, inoforcing,   &
    &                                    SUCCESS
  USE mo_model_domain,             ONLY: t_patch
  USE mo_util_string,              ONLY: int2string, logical2string, &
    &                                    real2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: imsg_thr
  PUBLIC :: itmr_thr
  PUBLIC :: idamtr
  PUBLIC :: istatus
  PUBLIC :: t_upatmo_dyn_config
  PUBLIC :: upatmo_dyn_config
  PUBLIC :: upatmo_exp_config
  PUBLIC :: t_upatmo_config
  PUBLIC :: upatmo_config
  PUBLIC :: configure_upatmo
  PUBLIC :: destruct_upatmo

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_upatmo_config'

  ! Note: upper-atmosphere physics, deep-atmosphere dynamics 
  ! and upper-atmosphere extrapolation are regarded as three elements 
  ! of the upper-atmosphere extension of ICON. 
  ! This is why we treat them in one module (as well as code economy).

  !------------------------------------------------------------
  !                      Parameter types
  !------------------------------------------------------------

  ! (The like you find in 'src/shared/mo_impl_constants'.)

  ! Thresholds...
  !
  TYPE t_ithr
    INTEGER :: low    ! Low threshold
    INTEGER :: med    ! Medium threshold
    INTEGER :: high   ! High threshold
  END TYPE t_ithr
  !
  ! ... for message output:
  TYPE(t_ithr), PARAMETER :: imsg_thr = t_ithr(  8, &  ! Low threshold for message output
    &                                           10, &  ! Med(ium) threshold for message output
    &                                           15  )  ! High threshold for message output  
  !
  ! ... for timers:
  TYPE(t_ithr), PARAMETER :: itmr_thr = t_ithr(  2, &  ! Low threshold for timer-call
    &                                            5, &  ! Med(ium) threshold for timer-call
    &                                            8  )  ! High threshold for timer-call  

  !------------------------------------------------------------

  ! Identifiers for deep-atmosphere metrical modification factors.
  !
  ! 1) Full levels, index order (jk, jtype)
  !
  TYPE t_idamtr_idxlist_type_1_mc
    INTEGER :: gradh     ! Horizontal derivatives 
    INTEGER :: divh      ! Horizontal part of divergence
    INTEGER :: vol       ! Cell volume
    INTEGER :: invr      ! = 1 / ( a + z )
    INTEGER :: centri    ! Centrifugal acceleration
    ! 
    INTEGER :: nitem     ! Number of identifiers
  END TYPE t_idamtr_idxlist_type_1_mc
  !
  ! 2) Half levels, index order (jk, jtype)
  !
  TYPE t_idamtr_idxlist_type_1_ifc
    INTEGER :: gradh     ! Horizontal derivatives 
    INTEGER :: invr      ! = 1 / ( a + z )   
    INTEGER :: centri    ! Centrifugal acceleration
    ! 
    INTEGER :: nitem     ! Number of identifiers
  END TYPE t_idamtr_idxlist_type_1_ifc
  !
  ! 3) Full levels, index order (jtype, jk)
  !
  TYPE t_idamtr_idxlist_type_2_mc
    INTEGER :: divzU     ! Vertical part of divergence (Upper interface of cell)
    INTEGER :: divzL     ! Vertical part of divergence (Lower interface of cell)
    !
    INTEGER :: nitem     ! Number of identifiers
  END TYPE t_idamtr_idxlist_type_2_mc
  !
  ! Collector
  !
  TYPE t_idamtr
    TYPE(t_idamtr_idxlist_type_1_mc)  :: t1mc
    TYPE(t_idamtr_idxlist_type_1_ifc) :: t1ifc
    TYPE(t_idamtr_idxlist_type_2_mc)  :: t2mc
  END TYPE t_idamtr
  !
  ! Assign values 
  ! (Please, update the 'nitem', if you modify the identifier lists. Thank you!) 
  !
  TYPE(t_idamtr), PARAMETER :: idamtr = t_idamtr(  &
    &                                   t_idamtr_idxlist_type_1_mc(  1,     &  ! idamtr%t1mc%gradh
    &                                                                2,     &  ! idamtr%t1mc%divh 
    &                                                                3,     &  ! idamtr%t1mc%vol 
    &                                                                4,     &  ! idamtr%t1mc%invr 
    &                                                                5,     &  ! idamtr%t1mc%centri 
    !
    &                                                                5  ),  &  ! idamtr%t1mc%nitem
    !-------------------------------------------------------------------------------
    &                                   t_idamtr_idxlist_type_1_ifc( 1,     &  ! idamtr%t1ifc%gradh
    &                                                                2,     &  ! idamtr%t1ifc%invr
    &                                                                3,     &  ! idamtr%t1ifc%centri
    !
    &                                                                3  ),  &  ! idamtr%t1ifc%nitem 
    !-------------------------------------------------------------------------------
    &                                   t_idamtr_idxlist_type_2_mc(  1,     &  ! idamtr%t2mc%divzU 
    &                                                                2,     &  ! idamtr%t2mc%divzL 
    !
    &                                                                2  )   &  ! idamtr%t2mc%nitem
    &                                              )  

  !------------------------------------------------------------

  ! Identifiers for configuration status.
  !
  TYPE t_istatus
    INTEGER :: checked      ! Upper-atmosphere namelist settings crosschecked?
    INTEGER :: configured   ! Upper-atmosphere configured?
    INTEGER :: required     ! Upper-atmosphere settings required at all?
    !
    INTEGER :: nitem        ! Number of identifiers
  END TYPE t_istatus
  TYPE(t_istatus), PARAMETER :: istatus = t_istatus( 1, &  ! checked
    &                                                2, &  ! configured
    &                                                3, &  ! required
    !
    &                                                3  )  ! nitem

  !------------------------------------------------------------
  !                    Configuration types
  !------------------------------------------------------------

  ! Namelist parameters that control the deep-atmosphere dynamics.
  !
  TYPE t_upatmo_dyn_config
    LOGICAL :: lnontrad        ! Switch for non-traditional deep-atmosphere terms 
                               ! in components of momentum equation
    LOGICAL :: lconstgrav      ! .TRUE. -> gravitational acceleration is const. 
                               ! like in case of the shallow atmosphere
    LOGICAL :: lcentrifugal    ! .TRUE. -> explicit centrifugal acceleration is switched on
    LOGICAL :: ldeepatmo2phys  ! .TRUE. -> the input fields to the ECHAM physics parameterizations 
                               ! are modified for the deep atmosphere, if required 
                               ! .FALSE. -> the input fields are computed in accordance 
                               ! with the shallow-atmosphere approximation (standard) in any case
  END TYPE t_upatmo_dyn_config
  
  !------------------------------------------------------------

  ! Namelist parameters that control the vertical extrapolation 
  ! of the initial atmospheric state into the upper atmosphere.
  !
  TYPE t_upatmo_exp_config
    REAL(wp) :: expol_start_height    ! [m] Height above which extrapolation (blending) of initial data starts
    REAL(wp) :: expol_blending_scale  ! [m] Blending scale height
    REAL(wp) :: expol_vn_decay_scale  ! [m] Scale height for (exponential) decay of extrapolated 
                                      ! horizontal wind component (for stability reasons)
    REAL(wp) :: expol_temp_infty      ! [K] Climatological temperature of exosphere (for z -> infinity)
    LOGICAL  :: lexpol_sanitycheck    ! .TRUE. -> Apply sanity check to extrapolated fields
  END TYPE t_upatmo_exp_config

  !------------------------------------------------------------

  ! For reasons of easier generalizability 
  ! (e.g., for the upcoming upper-atmosphere physics extension) 
  ! the variables of the above types become domain-dependent. 
  ! Since they are used, to store the namelist input, 
  ! i.e. at a point in the program sequence before the actual number of domains is known, 
  ! they have to be allocated for the maximum permissible domain range '0:max_dom'.
  ! A significant number of parameters are derived from the namelist settings 
  ! (especially in the context of the upcoming upper-atmosphere physics extension),  
  ! so they are not gathered in the above, but in the following type(s), 
  ! which can be allocated for the actual number of domains 'n_dom_start:n_dom', 
  ! in order to save some memory in standard simulations without upper-atmosphere.
  !
  TYPE t_phy
    LOGICAL :: l_constgrav    ! Const. gravitational acceleration for ECHAM physics
    LOGICAL :: l_shallowatmo  ! Shallow-atmosphere metrics for ECHAM physics
    ! Status variables
    LOGICAL :: l_status(istatus%nitem)
  END TYPE t_phy
  !
  TYPE t_dyn
    LOGICAL :: l_constgrav    ! .TRUE. -> gravitational acceleration is assumed to be constant
    LOGICAL :: l_centrifugal  ! .TRUE. -> centrifugal acceleration is switched on
    LOGICAL :: l_initonzgpot  ! .TRUE. -> initial data are living on geopotential heights
    ! Status variables
    LOGICAL :: l_status(istatus%nitem)
  END TYPE t_dyn
  !
  TYPE t_exp
    LOGICAL :: l_expol       ! .TRUE. -> upper-atmosphere-specific vertical extrapolation 
                             ! of initial data towards climatological values is switched on    
    INTEGER :: nexpollev     ! Index of grid layer for and above which the extrapolation 
                             ! of initial data will take place
    LOGICAL :: l_initicon_config = .FALSE. ! .TRUE. -> first step of configuration of extrapolation 
                             ! after namelist read-in has taken place
    ! Status variables
    LOGICAL :: l_status(istatus%nitem)
  END TYPE t_exp
  !
  ! Collector
  !
  TYPE t_upatmo_config
    TYPE(t_phy) :: phy
    TYPE(t_dyn) :: dyn
    TYPE(t_exp) :: exp
    ! Miscellaneous
    REAL(wp) :: dt_fastphy             ! On several occasions we need the domain-specific dtime 
                                       ! and 'atm_phy_nwp_config(jg)%dt_fastphy' may not be available
    REAL(wp) :: dt_dyn_nom             ! Domain-specific nominal dynamics time step
    ! Status variables
    ! (Although desirable, the creation 
    ! of a status object is currently too complicated 
    ! and not worth the effort.)
    LOGICAL :: l_status(istatus%nitem)
  END TYPE t_upatmo_config

  !------------------------------------------------------------

  TYPE(t_upatmo_dyn_config), TARGET              :: upatmo_dyn_config(0:max_dom)
  TYPE(t_upatmo_exp_config), TARGET              :: upatmo_exp_config(0:max_dom)
  TYPE(t_upatmo_config),     TARGET, ALLOCATABLE :: upatmo_config(:) 

  !------------------------------------------------------------

CONTAINS !..................................................................................

  !------------------------------------------------------------
  !                       Subroutines
  !------------------------------------------------------------

  !>
  !! This subroutine configures:
  !! - Upper-atmosphere physics
  !! - Deep-atmosphere dynamics
  !! - Upper-atmosphere extrapolation
  !!
  !! (Called in 'src/drivers/mo_atmo_nonhydrostatic: construct_atmo_nonhydrostatic')
  !!
  SUBROUTINE configure_upatmo( n_dom_start,   & !in
    &                          n_dom,         & !in
    &                          p_patch,       & !in
    &                          ldeepatmo,     & !in
    &                          init_mode,     & !in
    &                          iforcing,      & !in
    &                          dtime,         & !in
    &                          ndyn_substeps, & !in
    &                          flat_height,   & !in
    &                          msg_level,     & !in
    &                          vct_a          ) !(opt)in

    ! In/out variables
    INTEGER,            INTENT(IN) :: n_dom_start            ! Start index of domains
    INTEGER,            INTENT(IN) :: n_dom                  ! End index of domains
    TYPE(t_patch),      INTENT(IN) :: p_patch(n_dom_start:)  ! Domain properties    
    LOGICAL,            INTENT(IN) :: ldeepatmo              ! Switch for deep-atmosphere dynamics
    INTEGER,            INTENT(IN) :: init_mode              ! Initialization mode
    INTEGER,            INTENT(IN) :: iforcing               ! Switch for physics package (nwp, echam etc.) 
    REAL(wp),           INTENT(IN) :: dtime                  ! Fast-physics/advective time step for primary domain
    INTEGER,            INTENT(IN) :: ndyn_substeps          ! Number of dynamics' substeps per fast-physics time step
    REAL(wp),           INTENT(IN) :: flat_height            ! Below 'flat_height' grid layer interfaces follow topography
    INTEGER,            INTENT(IN) :: msg_level              ! Message level
    REAL(wp), OPTIONAL, INTENT(IN) :: vct_a(:)               ! Nominal heights of grid layer interfaces

    ! Local variables
    INTEGER :: jg, jg_ordered, jg_ref
    INTEGER :: nlev, nlevp1, nshift_total, n_dom_shift
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':configure_upatmo'

    !---------------------------------------------------------

    !-----------------------------------------------------
    !             Domain-independent checks
    !-----------------------------------------------------

    ! (We assume that the setup of 'p_patch' took already place.)

    IF (.NOT. ALLOCATED(upatmo_config)) THEN
      ! 'upatmo' should have been allocated in 'src/namelists/mo_upatmo_nml: check_upatmo'
      CALL finish(TRIM(routine), "Check calling sequence: upatmo_config is not allocated.")
    ELSEIF (.NOT. ANY((/0, 1/) == n_dom_start)) THEN
      ! For domain-dependent fields the max. range of allocation 
      ! with which we can reckon is '0:max_dom'
      !                              -
      CALL finish(TRIM(routine), 'Something has changed regarding n_dom_start.' )
    ELSEIF (.NOT. PRESENT(vct_a)) THEN
      ! For the setup of the upper-atmosphere extrapolation, we need 'vct_a', 
      ! and we take its absence as an indicator that it has not been allocated yet
      CALL finish(TRIM(routine), 'vct_a still uninitialized.')
    ENDIF

    ! Some of the settings below (especially for the upper-atmosphere extrapolation) 
    ! are more conveniently done, if the domain loop follows the order: 
    ! jg = 1, 2, 3, ..., n_dom, n_dom_start, 
    ! provided that 'n_dom_start = 0'. 
    n_dom_shift = 1 - n_dom_start

    ! Domain loop
    DO jg_ordered = n_dom_start, n_dom

      ! First, 'jg = 1' to 'n_dom' ... 
      jg = jg_ordered + n_dom_shift
      ! ... and then 'jg = 0'
      IF (jg > n_dom) jg = n_dom + 1 - jg

      ! Number of grid layers
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1
      
      ! Shift of grid layer index, 
      ! to account for vertical nesting
      IF (jg >= 1) THEN
        nshift_total = p_patch(jg)%nshift_total
      ELSE
        ! 'nshift_total' is not initialized for 'jg = n_dom_start', 
        ! if 'n_dom_start = 0', so we take the value from the primary domain 
        nshift_total = p_patch(1)%nshift_total
      ENDIF

      !-----------------------------------------------------
      !               Domain-dependent checks
      !-----------------------------------------------------
      
      IF ((jg < n_dom_start) .OR. (jg > n_dom)) THEN
        ! Some rudimentary checks of the reordering 
        ! of the domain sequence may be in order
        CALL finish(TRIM(routine), "Domain index jg has unexpected value.")
      ELSEIF ((jg_ordered == n_dom_start) .AND. (jg /= 1)) THEN
        CALL finish(TRIM(routine), "First domain to be configured has to be jg=1.")
      ELSEIF (.NOT. upatmo_config(jg)%l_status(istatus%checked)) THEN
        ! (Just to make sure. Actually it should have been set to .true. 
        ! right after the allocation of 'upatmo_config')
        CALL finish(TRIM(routine), "Check calling sequence: check_upatmo -> configure_upatmo.")  
      ELSEIF (.NOT. upatmo_config(jg)%exp%l_initicon_config) THEN
        ! For the final configuration of the upper-atmosphere extrapolation it is required 
        ! that the preliminary configuration in 'src/configure_model/mo_initicon_config' has taken place
        CALL finish(TRIM(routine), "Check calling sequence: configure_initicon -> configure_upatmo.")
      ENDIF

      !-----------------------------------------------------
      !                   Configuration
      !-----------------------------------------------------
      
      !---------------
      !   Dynamics
      !---------------

      CALL configure_upatmo_dynamics( ldeepatmo         = ldeepatmo,             & !in
        &                             init_mode         = init_mode,             & !in
        &                             upatmo_dyn_config = upatmo_dyn_config(jg), & !in
        &                             upatmo_config     = upatmo_config(jg)      ) !inout

      !---------------
      !    Physics
      !---------------

      CALL configure_upatmo_physics( ldeepatmo         = ldeepatmo,             & !in
        &                            iforcing          = iforcing,              & !in
        &                            upatmo_dyn_config = upatmo_dyn_config(jg), & !in
        &                            upatmo_config     = upatmo_config(jg)      ) !inout

      !---------------
      ! Extrapolation
      !---------------
      
      ! The primary domain is our reference
      ! (-> One reason for the reordering of the domain loop)
      jg_ref = 1

      ! Because of the consistency check within 'configure_upatmo_extrapolation', 
      ! the initialization of 'l_status' has to be done here
      CALL init_logical_1d( variable=upatmo_config(jg)%exp%l_status, value=.FALSE., &
        &                   opt_ilist=(/istatus%checked/), opt_mask="list"          )

      IF (jg == jg_ref) THEN
        CALL configure_upatmo_extrapolation( jg                = jg,                    & !in
          &                                  jg_ref            = jg_ref,                & !in
          &                                  nlev              = nlev,                  & !in
          &                                  nlevp1            = nlevp1,                & !in
          &                                  nshift_total      = nshift_total,          & !in
          &                                  flat_height       = flat_height,           & !in
          &                                  vct_a             = vct_a,                 & !(opt)in
          &                                  upatmo_exp_config = upatmo_exp_config(jg), & !in
          &                                  upatmo_config     = upatmo_config(jg)      ) !inout
      ELSE
        CALL configure_upatmo_extrapolation( jg                    = jg,                    & !in
          &                                  jg_ref                = jg_ref,                & !in
          &                                  nlev                  = nlev,                  & !in
          &                                  nlevp1                = nlevp1,                & !in
          &                                  nshift_total          = nshift_total,          & !in
          &                                  flat_height           = flat_height,           & !in
          &                                  vct_a                 = vct_a,                 & !(opt)in
          &                                  upatmo_exp_config     = upatmo_exp_config(jg), & !in
          &                                  upatmo_config         = upatmo_config(jg),     & !inout
          &                                  opt_upatmo_config_ref = upatmo_config(jg_ref)  ) !optin
      ENDIF

      !---------------
      ! Miscellaneous
      !---------------
      
      ! On several occasions we need the domain-specific value of dtime, 
      ! but 'iforcing /= inwp', so that 'atm_phy_nwp_config(jg)%dt_fastphy' is not available. 
      ! So we compute it here following its computation in 
      ! 'src/configure_model/mo_atm_phy_nwp_config: configure_atm_phy_nwp'. 
      ! A potential rescaling of 'dtime' with 'src/configure_model/mo_grid_config: grid_rescale_factor' 
      ! took already place in 'src/shared/mo_time_management: compute_timestep_settings', 
      ! which is called in 'src/configure_model/mo_nml_crosscheck: atm_crosscheck'.
      upatmo_config(jg)%dt_fastphy = dtime / 2._wp**(p_patch(jg)%level-p_patch(1)%level)
      
      ! We take the opportunity, to compute the nominal dynamics time step
      ! ("nominal", because 'ndyn_substeps' may change its value during runtime)
      upatmo_config(jg)%dt_dyn_nom = upatmo_config(jg)%dt_fastphy / REAL(ndyn_substeps, wp)

      ! 'upatmo_config' is allocated in any case, but not necessarily required
      upatmo_config(jg)%l_status(istatus%required) = upatmo_config(jg)%dyn%l_status(istatus%required) .OR. &
        &                                            upatmo_config(jg)%phy%l_status(istatus%required) .OR. &
        &                                            upatmo_config(jg)%exp%l_status(istatus%required)
      
      ! Indicate that configuration has taken place
      upatmo_config(jg)%l_status(istatus%configured) = upatmo_config(jg)%dyn%l_status(istatus%configured) .OR. &
        &                                              upatmo_config(jg)%phy%l_status(istatus%configured) .OR. &
        &                                              upatmo_config(jg)%exp%l_status(istatus%configured)

      !---------------
      !   Messages
      !---------------
      
      IF (upatmo_config(jg)%l_status(istatus%required)) THEN
        CALL print_message( jg,                & !in
          &                 iforcing,          & !in
          &                 nshift_total,      & !in
          &                 msg_level,         & !in
          &                 upatmo_config(jg), & !in
          &                 vct_a              ) !(opt)in
      ENDIF
      
    ENDDO  !jg_ordered
    
  END SUBROUTINE configure_upatmo

  !====================================================================================

  !>
  !! Configure deep-atmosphere dynamics.
  !!
  SUBROUTINE configure_upatmo_dynamics( ldeepatmo,         & !in
    &                                   init_mode,         & !in
    &                                   upatmo_dyn_config, & !in
    &                                   upatmo_config      ) !inout

    ! In/out variables
    LOGICAL,                   INTENT(IN)    :: ldeepatmo         ! Main deep-atmosphere switch
    INTEGER,                   INTENT(IN)    :: init_mode         ! Initialization mode
    TYPE(t_upatmo_dyn_config), INTENT(IN)    :: upatmo_dyn_config ! Namelist parameters
    TYPE(t_upatmo_config),     INTENT(INOUT) :: upatmo_config     ! Upper-atmosphere configuration

    ! Local variables
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':configure_upatmo_dynamics'

    !---------------------------------------------------------

    !-----------------------------------------------------
    !         Initialization with default values
    !-----------------------------------------------------

    upatmo_config%dyn%l_constgrav   = .TRUE.   ! Constant gravitational acceleration
    upatmo_config%dyn%l_centrifugal = .FALSE.  ! No explicit centrifugal acceleration
    upatmo_config%dyn%l_initonzgpot = .FALSE.  ! Initial data are living on geometric heights 
    CALL init_logical_1d( variable=upatmo_config%dyn%l_status, value=.FALSE., &
      &                   opt_ilist=(/istatus%checked/), opt_mask="list"      )

    !-----------------------------------------------------
    !                   Configuration
    !-----------------------------------------------------

    ! Check, if gravitational acceleration shall be constant in case of deep atmosphere, too
    upatmo_config%dyn%l_constgrav = MERGE( upatmo_dyn_config%lconstgrav, & ! True-case
      &                                    .TRUE.,                       & ! False-case
      &                                    ldeepatmo                     ) ! Condition
    
    ! Explicit centrifugal acceleration switched on?
    ! (If we assume that the parameter 'src/shared/mo_physical_constants: grav' contains 
    ! a contribution from the centrifugal acceleration, it is not subtracted out!)
    upatmo_config%dyn%l_centrifugal = MERGE( upatmo_dyn_config%lcentrifugal, & 
      &                                      .FALSE.,                        & 
      &                                      ldeepatmo                       ) 
    
    ! Should initial data be interpreted to live on geopotential instead of geometric heights?
    ! (In case of the shallow atmosphere, geometric and geopotential height would coincide,
    !  but this is unfortunately not the case for the deep atmosphere)
    ! Currently this variable is set to .true., if
    ! - the gravitational acceleration varies with height, 
    ! - and initial data are taken from IFS (MODE_IFSANA or MODE_COMBINED)
    upatmo_config%dyn%l_initonzgpot = .NOT. upatmo_config%dyn%l_constgrav .AND. &
      &                               ( init_mode == MODE_IFSANA .OR.           &
      &                                 init_mode == MODE_COMBINED              )
    
    ! 'upatmo_config' is allocated in any case, but not necessarily required
    upatmo_config%dyn%l_status(istatus%required) = ldeepatmo
    
    ! Indicate that configuration has taken place
    upatmo_config%dyn%l_status(istatus%configured) = .TRUE.

  END SUBROUTINE configure_upatmo_dynamics

  !====================================================================================

  !>
  !! Configure (upper-atmosphere) physics.
  !!
  SUBROUTINE configure_upatmo_physics( ldeepatmo,         & !in
    &                                  iforcing,          & !in
    &                                  upatmo_dyn_config, & !in
    &                                  upatmo_config      ) !inout

    ! In/out variables
    LOGICAL,                   INTENT(IN)    :: ldeepatmo         ! Main deep-atmosphere switch
    INTEGER,                   INTENT(IN)    :: iforcing          ! Switch for physics package (nwp, echam etc.) 
    TYPE(t_upatmo_dyn_config), INTENT(IN)    :: upatmo_dyn_config ! Namelist parameters
    TYPE(t_upatmo_config),     INTENT(INOUT) :: upatmo_config     ! Upper-atmosphere configuration

    ! Local variables
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':configure_upatmo_physics'

    !---------------------------------------------------------

    !-----------------------------------------------------
    !         Initialization with default values
    !-----------------------------------------------------

    upatmo_config%phy%l_constgrav   = .TRUE.  ! Constant gravitational acceleration in physics interface 
                                              ! (applies to ECHAM-physics only!)
    upatmo_config%phy%l_shallowatmo = .TRUE.  ! Shallow-atmosphere metrics in physics interface 
                                              ! (applies to ECHAM-physics only!)
    ! Initialize the status switches, 
    ! but skip 'istatus%checked', so as not to overwrite 
    ! its assignment in 'src/namelists/mo_upatmo_nml: check_upatmo'
    CALL init_logical_1d( variable=upatmo_config%phy%l_status, value=.FALSE., &
      &                   opt_ilist=(/istatus%checked/), opt_mask="list"      )

    !-----------------------------------------------------
    !                   Configuration
    !-----------------------------------------------------

    ! Should the input fields to the physics parameterizations 
    ! be modified for the deep atmosphere?
    ! Note: this applies to ECHAM-physics only
    IF (iforcing == iecham) THEN
      ! Gravitational acceleration:
      upatmo_config%phy%l_constgrav = MERGE( upatmo_config%dyn%l_constgrav,   & 
        &                                    .TRUE.,                          & 
        &                                    upatmo_dyn_config%ldeepatmo2phys ) 
      ! Metrics (concerns especially the cell volume):
      upatmo_config%phy%l_shallowatmo = MERGE( .NOT. ldeepatmo,                 &
        &                                      .TRUE.,                          &
        &                                      upatmo_dyn_config%ldeepatmo2phys )
    ELSE  ! (E.g., iforcing == inwp)
      ! Adopt the settings for the dynamics in any case
      upatmo_config%phy%l_constgrav   = upatmo_config%dyn%l_constgrav
      upatmo_config%phy%l_shallowatmo = .NOT. ldeepatmo
    ENDIF  !IF (iforcing == iecham)
    
    ! 'upatmo_config' is allocated in any case, but not necessarily required
    upatmo_config%phy%l_status(istatus%required) = ldeepatmo
    
    ! Indicate that configuration has taken place
    upatmo_config%phy%l_status(istatus%configured) = .TRUE.
    
  END SUBROUTINE configure_upatmo_physics

  !====================================================================================

  !>
  !! Configure upper-atmosphere extrapolation.
  !!
  SUBROUTINE configure_upatmo_extrapolation( jg,                   & !in
    &                                        jg_ref,               & !in
    &                                        nlev,                 & !in
    &                                        nlevp1,               & !in
    &                                        nshift_total,         & !in
    &                                        flat_height,          & !in
    &                                        vct_a,                & !in
    &                                        upatmo_exp_config,    & !in
    &                                        upatmo_config,        & !inout
    &                                        opt_upatmo_config_ref ) !optin

    ! In/out variables
    INTEGER,                             INTENT(IN)    :: jg, jg_ref            ! Domain indices
    INTEGER,                             INTENT(IN)    :: nlev, nlevp1          ! Number of grid layers and interfaces
    INTEGER,                             INTENT(IN)    :: nshift_total          ! Index shift for vertical nesting
    REAL(wp),                            INTENT(IN)    :: flat_height           ! Below 'flat_height' grid layer interfaces 
                                                                                ! follow topography
    REAL(wp),                  OPTIONAL, INTENT(IN)    :: vct_a(:)              ! Nominal heights of grid layer interfaces
    TYPE(t_upatmo_exp_config),           INTENT(IN)    :: upatmo_exp_config     ! Namelist parameters
    TYPE(t_upatmo_config),               INTENT(INOUT) :: upatmo_config         ! Upper-atmosphere configuration
    TYPE(t_upatmo_config),     OPTIONAL, INTENT(IN)    :: opt_upatmo_config_ref ! Configuration of reference domain 

    ! Local variables
    INTEGER :: jk, jks
    LOGICAL :: l_found
    LOGICAL :: l_present_ref, l_configured_ref, l_expol_ref
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':configure_upatmo_extrapolation'

    !---------------------------------------------------------

    ! The initialization of the extrapolation switch 'l_expol' 
    ! takes place in 'src/configure_model/mo_initicon_config', 
    ! depending on 'initicon_nml: itype_vert_expol'.
    ! Here, the final setting takes place.

    !-----------------------------------------------------
    !                 Consistency check
    !-----------------------------------------------------

    ! Actually, vct_a is only optional, because it is already optional one level higher
    IF (.NOT. PRESENT(vct_a)) CALL finish(TRIM(routine), 'vct_a has to be present.')

    IF (PRESENT(opt_upatmo_config_ref)) THEN
      l_present_ref    = .TRUE.
      l_configured_ref = opt_upatmo_config_ref%exp%l_status(istatus%configured)
      l_expol_ref      = opt_upatmo_config_ref%exp%l_expol
    ELSE
      l_present_ref    = .FALSE.
      l_configured_ref = .FALSE.
      l_expol_ref      = .FALSE.
    ENDIF

    IF (l_present_ref .EQV. (jg == jg_ref)) THEN
      ! The reference configuration has to be present, 
      ! if the current domain is not the reference domain, but only then
      CALL finish(TRIM(routine), "Presence/absence of optional reference has to coincide with jg/=jg_ref/jg=jg_ref.")
    ELSEIF (l_present_ref .AND. (.NOT. l_configured_ref)) THEN
      ! The reference domain should be the first one, 
      ! which enters this subroutine
      CALL finish(TRIM(routine), "Optional reference is present, but not yet configured.")
    ENDIF

    !-----------------------------------------------------
    !         Initialization with default values
    !-----------------------------------------------------

    ! Please, do not set 'upatmo_config%exp%l_expol' here, 
    ! otherwise, you overwrite the assignment 
    ! in 'src/configure_model/mo_initicon_config: configure_initicon'
    upatmo_config%exp%nexpollev = -1         

    !-----------------------------------------------------
    !                   Configuration
    !-----------------------------------------------------

    IF ((jg /= jg_ref) .AND. (.NOT. l_expol_ref)) THEN
      ! If extrapolation has been switched off for domain 1 (reference),    
      ! it will be switched off for all other domains
      upatmo_config%exp%l_expol = .FALSE. 
    ELSEIF (upatmo_exp_config%expol_start_height < flat_height .OR.     &
      &     upatmo_exp_config%expol_start_height > vct_a(1+nshift_total)) THEN
      ! If extrapolation start heigth 'expol_start_height' is below 'flat_height' or
      ! above the domain top, no extrapolation of the upper-atmosphere type should be necessary
      upatmo_config%exp%l_expol = .FALSE. 
      ! In all other cases the initial setting of 'l_expol' remains
    ENDIF
    
    ! Determine vertical index of grid layer for and above which 
    ! the extrapolation of the upper-atmosphere type of initial data will take place 
    IF (upatmo_config%exp%l_expol .AND. (jg==jg_ref)) THEN
      l_found = .FALSE.
      DO jk=1, nlevp1
        jks = jk + nshift_total  ! (Not really necessary, but for completeness)
        IF (vct_a(jks) < upatmo_exp_config%expol_start_height) THEN
          ! The grid cell layer within which 'expol_start_height' is located 
          ! shall be the lowermost layer to which the upper-atmosphere extrapolation is applied, 
          ! and the first interface above 'expol_start_height' is the lowermost to which 
          ! the extrapolation is applied, i.e.'nexpollev' applies to full and half(!) levels
          upatmo_config%exp%nexpollev = jk - 1
          ! Indicate the find
          l_found = .TRUE.
          EXIT
        ENDIF
      ENDDO  !jk
      IF (.NOT. l_found) CALL finish(TRIM(routine), 'Could not find nexpollev.')
      ! Just to make sure
      upatmo_config%exp%nexpollev = MIN( MAX( 1, upatmo_config%exp%nexpollev ), nlev )
    ELSEIF (upatmo_config%exp%l_expol .AND. (jg/=jg_ref)) THEN
      ! In the current configuration of the upper-atmosphere namelist, 
      ! no domain-wise specification of the extrapolation parameters is possible. 
      ! So the settings for jg /= jg_ref should be derivable 
      ! from the setting for the primary domain
      upatmo_config%exp%nexpollev = opt_upatmo_config_ref%exp%nexpollev - nshift_total
    ENDIF  !IF (upatmo_config%exp%l_expol .AND. (jg==jg_ref))
    
    ! 'upatmo_config' is allocated in any case, but not necessarily required
    upatmo_config%exp%l_status(istatus%required) = upatmo_config%exp%l_expol
    
    ! Indicate that configuration has taken place
    upatmo_config%exp%l_status(istatus%configured) = .TRUE.
    
  END SUBROUTINE configure_upatmo_extrapolation

  !====================================================================================

  !>
  !! Print message on upper-atmosphere configuration.
  !!
  SUBROUTINE print_message( jg,            & !in
    &                       iforcing,      & !in
    &                       nshift_total,  & !in
    &                       msg_level,     & !in
    &                       upatmo_config, & !in
    &                       vct_a          ) !(opt)in

    ! In/out variables
    INTEGER,                         INTENT(IN) :: jg            ! Domain index
    INTEGER,                         INTENT(IN) :: iforcing      ! Switch for physics package (nwp, echam etc.)  
    INTEGER,                         INTENT(IN) :: nshift_total  ! Shift of vertical grid index for vertical nesting
    INTEGER,                         INTENT(IN) :: msg_level     ! Message level
    TYPE(t_upatmo_config),           INTENT(IN) :: upatmo_config ! Upper-atmosphere configuration
    REAL(wp),              OPTIONAL, INTENT(IN) :: vct_a(:)      ! Nominal heights of grid layer interfaces

    ! Local variables
    LOGICAL :: l_onlyPrimDom
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: msg_prefix
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':print_message'

    !---------------------------------------------------------

    ! Some print is only required for the primary domain
    l_onlyPrimDom = (jg==1)

    !-----------------------------------------------------
    !                "Any-case" output
    !-----------------------------------------------------

    ! Deep atmosphere:
    IF (l_onlyPrimDom .AND. upatmo_config%dyn%l_status(istatus%required)) THEN
      CALL message(TRIM(routine), "Deep-atmosphere modification of non-hydrostatic atmosphere switched on.")
      CALL message(TRIM(routine), "Please note: for efficiency reasons and code economy"// &
        & " the deep-atmosphere modification of the dynamical core disregards:")
      CALL message("", " - horizontal variation of grid layer heights due to terrain")
      CALL message("", " - any kind of diffusion, damping and the like (including LES physics)")
      CALL message("", " - special numerical 'tricks' beyond the main dynamics line,"//&
        & " such as sub-stepping for tracer advection")
      CALL message("", " - the feedback procedures for state relaxation between domains")
      IF (iforcing /= inoforcing) CALL message(TRIM(routine), "Please note: no physics parameterization"// &
        & " is modified for the deep atmosphere!")
    ENDIF

    ! Miscellaneous:
    IF (l_onlyPrimDom .AND. upatmo_config%l_status(istatus%required)) THEN
      CALL message(TRIM(routine), "(Info: most upper-atmosphere-related message output requires msg_level >= " & 
        & //TRIM(int2string(imsg_thr%high))//")")
      CALL message(TRIM(routine), "(Info: most upper-atmosphere-related timers require timers_level >= " &
        & //TRIM(int2string(itmr_thr%med))//")")
    ENDIF

    IF (msg_level >= imsg_thr%low) THEN

      !-----------------------------------------------------
      !              "High-priority" output
      !-----------------------------------------------------

      IF (msg_level >= imsg_thr%med) THEN

        !-----------------------------------------------------
        !             "Medium-priority" output
        !-----------------------------------------------------
        
        IF (msg_level >= imsg_thr%high) THEN 

          !-----------------------------------------------------
          !               "Low-priority" output
          !-----------------------------------------------------

          ! Deep atmosphere:
          IF (upatmo_config%dyn%l_status(istatus%required)) THEN
            msg_prefix = 'upatmo_config('//TRIM(int2string(jg))//')%dyn%'
            ! 
            WRITE(message_text,'(a)') TRIM(msg_prefix)//'l_constgrav: '// &
              & TRIM(logical2string(upatmo_config%dyn%l_constgrav))
            CALL message(TRIM(routine), TRIM(message_text))
            !
            WRITE(message_text,'(a)') TRIM(msg_prefix)//'l_centrifugal: '// &
              & TRIM(logical2string(upatmo_config%dyn%l_centrifugal))
            CALL message(TRIM(routine), TRIM(message_text))
            !
            WRITE(message_text,'(a)') TRIM(msg_prefix)//'l_initonzgpot: '// &
              & TRIM(logical2string(upatmo_config%dyn%l_initonzgpot))
            CALL message(TRIM(routine), TRIM(message_text))
          ENDIF

          ! Physics:
          IF (upatmo_config%phy%l_status(istatus%required)) THEN
            msg_prefix = 'upatmo_config('//TRIM(int2string(jg))//')%phy%'
            !
            WRITE(message_text,'(a)') TRIM(msg_prefix)//'l_constgrav: '// &
              & TRIM(logical2string(upatmo_config%phy%l_constgrav))
            CALL message(TRIM(routine), TRIM(message_text))
            !
            WRITE(message_text,'(a)') TRIM(msg_prefix)//'l_shallowatmo: '// &
              & TRIM(logical2string(upatmo_config%phy%l_shallowatmo))
            CALL message(TRIM(routine), TRIM(message_text))
          ENDIF

          ! Upper-atmosphere extrapolation:
          IF (upatmo_config%exp%l_status(istatus%required) .AND. PRESENT(vct_a)) THEN
            msg_prefix = 'upatmo-expol('//TRIM(int2string(jg))//'): '
            WRITE(message_text,'(a)') TRIM(msg_prefix)//'nexpollev: '// &
              & TRIM(int2string(upatmo_config%exp%nexpollev))
            CALL message(TRIM(routine), TRIM(message_text))
            WRITE(message_text,'(a)') TRIM(msg_prefix)//'interface height above which extrapolation '// &
              & 'potentially takes place: '//                                                           &
              & TRIM(real2string(vct_a(upatmo_config%exp%nexpollev + nshift_total)))
            CALL message(TRIM(routine), TRIM(message_text))
          ENDIF

        ENDIF  !imsg_thr%high
      ENDIF  !imsg_thr%med
    ENDIF  !imsg_thr%low

  END SUBROUTINE print_message

  !====================================================================================

  !>
  !! Initialize logical 1d-array.
  !! (Introduced, because 'src/shared/mo_fortran_tools: init_contiguous_l' 
  !! does not suit our purposes.)
  !!
  SUBROUTINE init_logical_1d( variable,   & !inout
    &                         value,      & !in
    &                         opt_ilist,  & !optin
    &                         opt_istart, & !optin
    &                         opt_mask    ) !optin
    ! In/out variables
    LOGICAL,                    INTENT(INOUT) :: variable(:)  ! Logical array to be assigned with 'value'
    LOGICAL,                    INTENT(IN)    :: value         
    INTEGER,          OPTIONAL, INTENT(IN)    :: opt_ilist(:) ! Optional list with indices of 'variable' 
                                                              ! that shall or shall not be assigned with 'value'. 
                                                              ! The indices are assumed to be 
                                                              ! in '[opt_istart, opt_istart+size(variable)-1]', 
                                                              ! if present, or in '[1, SIZE(variable)]' otherwise.
    INTEGER,          OPTIONAL, INTENT(IN)    :: opt_istart   ! Optional input, if "true" index range 
                                                              ! of 'variable' does not start with 1
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: opt_mask     ! "list" -> those indices of 'variable' stored 
                                                              ! in 'opt_ilist' are not assigned with 'value'
                                                              ! "complement" -> those indices of 'variable' 
                                                              ! not stored in 'opt_ilist' are not assigned with 'value'. 
                                                              ! The case that 'opt_ilist' is present, 
                                                              ! while 'opt_mask' is absent, is interpreted 
                                                              ! as 'opt_mask = "list"'.

    ! Local variables
    LOGICAL, ALLOCATABLE :: mask(:)
    LOGICAL :: lmask
    INTEGER :: varsize, istart, iend, ishift, jloop, istat
    INTEGER, PARAMETER :: MASKLEN = 20
    CHARACTER(LEN=MASKLEN), PARAMETER :: mask_list       = "list"
    CHARACTER(LEN=MASKLEN), PARAMETER :: mask_complement = "complement"
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':init_logical_1d'

    !---------------------------------------------------------

    varsize = SIZE(variable)

    IF (PRESENT(opt_istart)) THEN
      istart = opt_istart
    ELSE
      istart = 1
    ENDIF

    iend   = istart + varsize - 1
    ishift = 1 - istart

    IF (PRESENT(opt_ilist)) THEN
      lmask = .TRUE.
      IF ( MINVAL(opt_ilist) < istart .OR. &
        &  MAXVAL(opt_ilist) > iend        ) THEN
        CALL finish(TRIM(routine), "Index in opt_ilist outside index range of variable.")
      ENDIF
      ALLOCATE(mask(varsize), STAT=istat)
      IF (istat /= SUCCESS) CALL finish(TRIM(routine), "Allocation of mask failed.")

      IF (PRESENT(opt_mask)) THEN
        SELECT CASE(TRIM(opt_mask))
        CASE(TRIM(mask_list))
          mask(:) = .FALSE.
        CASE(TRIM(mask_complement))
          mask(:) = .TRUE.
        CASE default
          CALL finish(TRIM(routine), "Invalid opt_mask.")
        END SELECT
      ELSE
        mask(:) = .FALSE.
      ENDIF

      DO jloop = 1, SIZE(opt_ilist)
        mask(opt_ilist(jloop) - ishift) = .NOT. mask(opt_ilist(jloop) - ishift)
      ENDDO
    ELSE
      lmask = .FALSE.
    ENDIF
    
    IF (lmask) THEN
      DO jloop = 1, varsize
        IF (.NOT. mask(jloop)) variable(jloop) = value
      ENDDO
    ELSE
      DO jloop = 1, varsize
        variable(jloop) = value
      ENDDO
    ENDIF

    IF (lmask) THEN
      DEALLOCATE(mask, STAT=istat)
      IF (istat /= SUCCESS) CALL finish(TRIM(routine), "Deallocation of mask failed.")
    ENDIF

  END SUBROUTINE init_logical_1d

  !====================================================================================

  !>
  !! Destruct upper-atmosphere configuration.
  !!
  SUBROUTINE destruct_upatmo() 

    ! Local variables
    INTEGER :: istat
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':destruct_upatmo'

    !---------------------------------------------------------

    ! Deallocate upatmo_config
    ! (Allocated in 'src/namelists/mo_upatmo_nml: check_upatmo')
    IF (ALLOCATED(upatmo_config)) THEN
      DEALLOCATE(upatmo_config, STAT=istat)
      IF (istat /= SUCCESS) CALL finish(TRIM(routine), "Deallocation of upatmo_config failed.")
    ENDIF

  END SUBROUTINE destruct_upatmo

END MODULE mo_upatmo_config
