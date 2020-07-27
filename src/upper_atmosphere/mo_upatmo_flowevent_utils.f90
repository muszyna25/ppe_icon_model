!>
!! Tools for the upper-atmosphere 
!! that could not be placed in 'mo_upatmo_utils', 
!! due to circular dependencies, unfortunately.
!!
!! @author Sebastian Borchert (DWD)
!!
!! @par Revision History
!! Initial release by Sebastian Borchert, DWD (2019-10-28)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_upatmo_flowevent_utils

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message
  USE mo_impl_constants,        ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_upatmo_impl_const,     ONLY: iUpatmoStat
  USE mo_upatmo_types,          ONLY: t_upatmo
  USE mo_restart_attributes,    ONLY: t_RestartAttributeList, &
    &                                 getAttributesForRestarting
  USE mo_upatmo_config,         ONLY: upatmo_config
  USE mo_upatmo_utils,          ONLY: t_varstate_set
  USE mtime,                    ONLY: datetime
  USE mo_packed_message,        ONLY: t_PackedMessage
  USE mo_util_string,           ONLY: int2string
  USE mo_fortran_tools,         ONLY: assign_if_present_allocatable, &
    &                                 DO_DEALLOCATE

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_upatmoRestartAttributes
  PUBLIC :: upatmoRestartAttributesPrepare
  PUBLIC :: upatmoRestartAttributesPack
  PUBLIC :: upatmoRestartAttributesAssign
  PUBLIC :: upatmoRestartAttributesSet
  PUBLIC :: upatmoRestartAttributesGet
  PUBLIC :: upatmoRestartAttributesDeallocate

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_upatmo_flowevent_utils'

  !====================================================================================

  ! For a restart:
  ! The procedure to transfer metadata to the restart file is extremely complex. 
  ! We were not yet able to deduce how to use it correctly 
  ! (let alone how it works). With the following type and subroutines,
  ! we try to emulate examples of usage in 'src/io/restart' 
  ! in the hope that we obtain a sufficent result.

  ! All attemps to make direct use of src/upper_atmosphere/mo_upatmo_utils: t_varstate_set 
  ! for tendStateSet in t_upatmoRestartAttributes below were unsuccessful 
  ! (e.g., t_PackedMessage raised trouble). So we had no choice but to explicitly repeat 
  ! the content of t_varstate_set in t_upatmoRestartAttributes 
  ! and in all subroutines below.

  ! Important: if more than one nest is used, 
  ! some seemingly arbitrary entries of upatmoRestartAttributes%tendStateSet_...
  ! are no longer written to the restart files. 
  ! We were not yet able to identify the reason for that. 
  ! Since the restart infrastructure is far too complicated
  ! to figure out a workaround, we cannot apply a restart 
  ! for upatmoRestartAttributes%tendStateSet_... of domains with index jg > 2, 
  ! but have to fall back on the "cold start" initialization.
  ! I.e. simulations with n_dom > 2 are not restart-safe!

  ! Please note: if the content of an instance of t_upatmoRestartAttributes 
  ! is allocated, no explicit deallocation takes place in most cases.

  TYPE :: t_upatmoRestartAttributes
    REAL(wp), ALLOCATABLE :: elapsedTimePhy(:)
    REAL(wp), ALLOCATABLE :: elapsedTimeExtdat(:)
    INTEGER,  ALLOCATABLE :: tendStateSet_i_old(:)
    INTEGER,  ALLOCATABLE :: tendStateSet_i_new(:)
    INTEGER,  ALLOCATABLE :: tendStateSet_n_swap(:)
    INTEGER,  ALLOCATABLE :: tendStateSet_n_state(:)
    INTEGER,  ALLOCATABLE :: tendStateSet_n_statep1(:)
    LOGICAL,  ALLOCATABLE :: tendStateSet_l_swapped(:)
    LOGICAL,  ALLOCATABLE :: tendStateSet_l_updated(:)
    LOGICAL,  ALLOCATABLE :: tendStateSet_l_locked(:)
    LOGICAL,  ALLOCATABLE :: tendStateSet_l_locking(:)
    LOGICAL,  ALLOCATABLE :: tendStateSet_l_unlockable(:)
    LOGICAL,  ALLOCATABLE :: tendStateSet_l_final(:)
    LOGICAL,  ALLOCATABLE :: tendStateSet_l_finish_on_error(:)
    LOGICAL,  ALLOCATABLE :: tendStateSet_l_initialized(:)
  END TYPE t_upatmoRestartAttributes

  CHARACTER(LEN=*), PARAMETER :: keyStrElapsedTimePhy                 = "upatmo_elapsedTimePhy_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrElapsedTimeExtdat              = "upatmo_elapsedTimeExtdat_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_i_old             = "upatmo_tendStateSet_i_old_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_i_new             = "upatmo_tendStateSet_i_new_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_n_swap            = "upatmo_tendStateSet_n_swap_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_n_state           = "upatmo_tendStateSet_n_state_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_n_statep1         = "upatmo_tendStateSet_n_statep1_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_l_swapped         = "upatmo_tendStateSet_l_swapped_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_l_updated         = "upatmo_tendStateSet_l_updated_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_l_locked          = "upatmo_tendStateSet_l_locked_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_l_locking         = "upatmo_tendStateSet_l_locking_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_l_unlockable      = "upatmo_tendStateSet_l_unlockable_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_l_final           = "upatmo_tendStateSet_l_final_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_l_finish_on_error = "upatmo_tendStateSet_l_finish_on_error_DOM"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSet_l_initialized     = "upatmo_tendStateSet_l_initialized_DOM"
  !
  CHARACTER(LEN=*), PARAMETER :: keyStrElapsedTimeSuffix  = "_PHY"
  CHARACTER(LEN=*), PARAMETER :: keyStrTendStateSetSuffix = "_TEND"

  INTEGER, PARAMETER :: domRestartLimit = 2

  INTERFACE setRestartAttributes
    MODULE PROCEDURE setRestartAttributes_R1D
    MODULE PROCEDURE setRestartAttributes_I1D
    MODULE PROCEDURE setRestartAttributes_L1D
  END INTERFACE

CONTAINS !..................................................................................

  !====================================================================================

  !>
  !! Prepare metadata for restart file.
  !!
  !! Called in: src/atm_dyn_iconam/mo_nh_stepping: perform_nh_timeloop
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert (DWD) (2019-10-28)
  !!
  SUBROUTINE upatmoRestartAttributesPrepare( jg,                      & !in
    &                                        upatmoRestartAttributes, & !inout
    &                                        prm_upatmo,              & !inout
    &                                        mtime_current            ) !in
    
    ! In/out variables
    INTEGER,                         INTENT(IN)    :: jg
    TYPE(t_upatmoRestartAttributes), INTENT(INOUT) :: upatmoRestartAttributes
    TYPE(t_upatmo),                  INTENT(INOUT) :: prm_upatmo
    TYPE(datetime), POINTER,         INTENT(IN)    :: mtime_current

    ! Local variables
    TYPE(t_varstate_set) :: tendStateSet
    INTEGER :: i, istat, nsize
    LOGICAL :: lmessage
    CHARACTER(LEN=2) :: domStr
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':upatmoRestartAttributesPrepare'

    !----------------------------------------------

    ! Message output desired?
    lmessage = upatmo_config(jg)%l_status(iUpatmoStat%message)

    domStr = TRIM(int2string(jg))

    IF (lmessage) CALL message(TRIM(routine), &
      & 'Start preparation of metadata for restart file on domain '//domStr)

    ! If the current simulation is a multi-domain application, 
    ! the argument upatmoRestartAttributes may have been used several times 
    ! on the invocation site. So we should deallocate its content, if required.
    CALL upatmoRestartAttributesDeallocate(upatmoRestartAttributes)
    
    ! Event management object: get elapsed time since last trigger date 
    CALL upatmo_config(jg)%nwp_phy%event_mgmt_grp%serialize(mtime_current, &
      & upatmoRestartAttributes%elapsedTimePhy                             )
    CALL upatmo_config(jg)%nwp_phy%event_mgmt_extdat%serialize(mtime_current, &
      & upatmoRestartAttributes%elapsedTimeExtdat                             )
    
    ! The restart mechanism can only be applied to upatmoRestartAttributes%tendStateSet_... for jg < 3. 
    ! For jg >= 3, seemingly arbitrary entries of upatmoRestartAttributes%tendStateSet_... are missing 
    ! from the attribute lists of the restart files for not yet identified reasons.
    IF (jg <= domRestartLimit) THEN

      ! Get info about state of accumulative tendencies
      nsize = SIZE(prm_upatmo%tend%ddt%state)
      IF (nsize < 1) CALL finish(TRIM(routine), 'SIZE(prm_upatmo%tend%ddt%state) < 1')
      ALLOCATE( upatmoRestartAttributes%tendStateSet_i_old(nsize),             &
        &       upatmoRestartAttributes%tendStateSet_i_new(nsize),             &
        &       upatmoRestartAttributes%tendStateSet_n_swap(nsize),            &
        &       upatmoRestartAttributes%tendStateSet_n_state(nsize),           &
        &       upatmoRestartAttributes%tendStateSet_n_statep1(nsize),         &
        &       upatmoRestartAttributes%tendStateSet_l_swapped(nsize),         &
        &       upatmoRestartAttributes%tendStateSet_l_updated(nsize),         &
        &       upatmoRestartAttributes%tendStateSet_l_locked(nsize),          &
        &       upatmoRestartAttributes%tendStateSet_l_locking(nsize),         &
        &       upatmoRestartAttributes%tendStateSet_l_unlockable(nsize),      &
        &       upatmoRestartAttributes%tendStateSet_l_final(nsize),           &
        &       upatmoRestartAttributes%tendStateSet_l_finish_on_error(nsize), &
        &       upatmoRestartAttributes%tendStateSet_l_initialized(nsize),     &
        &       STAT=istat                                                     )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
        & 'Allocation of upatmoRestartAttributes%tendStateSet... failed')
      
      DO i = 1, nsize
        tendStateSet = prm_upatmo%tend%ddt%state(i)%getSet(optWhichSet="set4use")
        upatmoRestartAttributes%tendStateSet_i_old(i)             = tendStateSet%i_old
        upatmoRestartAttributes%tendStateSet_i_new(i)             = tendStateSet%i_new
        upatmoRestartAttributes%tendStateSet_n_swap(i)            = tendStateSet%n_swap
        upatmoRestartAttributes%tendStateSet_n_state(i)           = tendStateSet%n_state
        upatmoRestartAttributes%tendStateSet_n_statep1(i)         = tendStateSet%n_statep1
        upatmoRestartAttributes%tendStateSet_l_swapped(i)         = tendStateSet%l_swapped
        upatmoRestartAttributes%tendStateSet_l_updated(i)         = tendStateSet%l_updated
        upatmoRestartAttributes%tendStateSet_l_locked(i)          = tendStateSet%l_locked
        upatmoRestartAttributes%tendStateSet_l_locking(i)         = tendStateSet%l_locking
        upatmoRestartAttributes%tendStateSet_l_unlockable(i)      = tendStateSet%l_unlockable
        upatmoRestartAttributes%tendStateSet_l_final(i)           = tendStateSet%l_final
        upatmoRestartAttributes%tendStateSet_l_finish_on_error(i) = tendStateSet%l_finish_on_error
        upatmoRestartAttributes%tendStateSet_l_initialized(i)     = tendStateSet%l_initialized
      ENDDO

    ENDIF  !IF (jg <= domRestartLimit)

    IF (lmessage) CALL message(TRIM(routine), &
      & 'Finish preparation of metadata for restart file on domain '//domStr)

  END SUBROUTINE upatmoRestartAttributesPrepare

  !***************************************************************

  !>
  !! This subroutine is a structural copy of  
  !! src/io/restart/mo_restart_patch_description: restartPatchDescription_packer, 
  !! from where it is called.
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert (DWD) (2019-10-28)
  !!
  SUBROUTINE upatmoRestartAttributesPack( jg,                      & !in
    &                                     upatmoRestartAttributes, & !inout
    &                                     packedMessage,           & !inout
    &                                     operation                ) !value

    ! In/out variables
    INTEGER,                         INTENT(IN)    :: jg
    TYPE(t_upatmoRestartAttributes), INTENT(INOUT) :: upatmoRestartAttributes
    TYPE(t_PackedMessage),           INTENT(INOUT) :: packedMessage
    INTEGER,                         VALUE         :: operation

    !----------------------------------------------

    CALL packedMessage%packer(operation, upatmoRestartAttributes%elapsedTimePhy)
    CALL packedMessage%packer(operation, upatmoRestartAttributes%elapsedTimeExtdat)
    IF (jg <= domRestartLimit) THEN
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_i_old)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_i_new)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_n_swap)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_n_state)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_n_statep1)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_l_swapped)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_l_updated)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_l_locked)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_l_locking)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_l_unlockable)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_l_final)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_l_finish_on_error)
      CALL packedMessage%packer(operation, upatmoRestartAttributes%tendStateSet_l_initialized)
    ENDIF

  END SUBROUTINE upatmoRestartAttributesPack

  !***************************************************************

  !>
  !! Called in: src/io/restart/mo_restart_patch_description: restartPatchDescription_update
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert (DWD) (2019-10-28)
  !!
  SUBROUTINE upatmoRestartAttributesAssign( jg,                        & !in
    &                                       upatmoRestartAttributes,   & !inout
    &                                       optUpatmoRestartAttributes ) !optin

    ! In/out variables
    INTEGER,                                   INTENT(IN)    :: jg
    TYPE(t_upatmoRestartAttributes),           INTENT(INOUT) :: upatmoRestartAttributes
    TYPE(t_upatmoRestartAttributes), OPTIONAL, INTENT(IN)    :: optUpatmoRestartAttributes

    !----------------------------------------------

    IF (PRESENT(optUpatmoRestartAttributes)) THEN

      CALL assign_if_present_allocatable( upatmoRestartAttributes%elapsedTimePhy,   &
        &                                 optUpatmoRestartAttributes%elapsedTimePhy )
      CALL assign_if_present_allocatable( upatmoRestartAttributes%elapsedTimeExtdat,   &
        &                                 optUpatmoRestartAttributes%elapsedTimeExtdat )
      IF (jg <= domRestartLimit) THEN
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_i_old,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_i_old )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_i_new,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_i_new )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_n_swap,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_n_swap )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_n_state,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_n_state )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_n_statep1,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_n_statep1 )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_l_swapped,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_l_swapped )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_l_updated,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_l_updated )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_l_locked,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_l_locked )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_l_locking,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_l_locking )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_l_unlockable,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_l_unlockable )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_l_final,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_l_final )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_l_finish_on_error,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_l_finish_on_error )
        CALL assign_if_present_allocatable( upatmoRestartAttributes%tendStateSet_l_initialized,   &
          &                                 optUpatmoRestartAttributes%tendStateSet_l_initialized )
      ENDIF

    ENDIF

  END SUBROUTINE upatmoRestartAttributesAssign

  !***************************************************************

  !>
  !! Called in: src/io/restart/mo_restart_patch_description: restartPatchDescription_setRestartAttributes
  !! and src/io/restart/mo_sync_restart: syncRestartDescriptor_defineRestartAttributes
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert (DWD) (2019-10-28)
  !!
  SUBROUTINE upatmoRestartAttributesSet( jg,                      & !in
    &                                    upatmoRestartAttributes, & !in
    &                                    restartAttributes        ) !inout

    ! In/out variables
    INTEGER,                         INTENT(IN)    :: jg
    TYPE(t_upatmoRestartAttributes), INTENT(IN)    :: upatmoRestartAttributes
    TYPE(t_RestartAttributeList),    INTENT(INOUT) :: restartAttributes

    ! Local variables
    CHARACTER(LEN=2) :: domStr
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: keyStr

    !----------------------------------------------

    domStr = TRIM(int2string(jg, "(i2.2)"))
    
    ! elapsedTimePhy:
    keyStr = TRIM(keyStrElapsedTimePhy)//TRIM(domStr)//TRIM(keyStrElapsedTimeSuffix)
    CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%elapsedTimePhy, keyStr)
    ! elapsedTimeExtdat:
    keyStr = TRIM(keyStrElapsedTimeExtdat)//TRIM(domStr)//TRIM(keyStrElapsedTimeSuffix)
    CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%elapsedTimeExtdat, keyStr)
    IF (jg <= domRestartLimit) THEN
      ! i_old:
      keyStr = TRIM(keyStrTendStateSet_i_old)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_i_old, keyStr)
      ! i_new:
      keyStr = TRIM(keyStrTendStateSet_i_new)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_i_new, keyStr)
      ! n_swap:
      keyStr = TRIM(keyStrTendStateSet_n_swap)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_n_swap, keyStr)
      ! n_state:
      keyStr = TRIM(keyStrTendStateSet_n_state)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_n_state, keyStr)
      ! n_statep1:
      keyStr = TRIM(keyStrTendStateSet_n_statep1)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_n_statep1, keyStr)
      ! l_swapped:
      keyStr = TRIM(keyStrTendStateSet_l_swapped)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_l_swapped, keyStr)
      ! l_updated:
      keyStr = TRIM(keyStrTendStateSet_l_updated)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_l_updated, keyStr)
      ! l_locked:
      keyStr = TRIM(keyStrTendStateSet_l_locked)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_l_locked, keyStr)
      ! l_locking:
      keyStr = TRIM(keyStrTendStateSet_l_locking)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_l_locking, keyStr)
      ! l_unlockable:
      keyStr = TRIM(keyStrTendStateSet_l_unlockable)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_l_unlockable, keyStr)
      ! l_final:
      keyStr = TRIM(keyStrTendStateSet_l_final)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_l_final, keyStr)
      ! l_finish_on_error:
      keyStr = TRIM(keyStrTendStateSet_l_finish_on_error)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_l_finish_on_error, keyStr)
      ! l_initialized:
      keyStr = TRIM(keyStrTendStateSet_l_initialized)//TRIM(domStr)//TRIM(keyStrTendStateSetSuffix)
      CALL setRestartAttributes(restartAttributes, upatmoRestartAttributes%tendStateSet_l_initialized, keyStr)
    ENDIF

  END SUBROUTINE upatmoRestartAttributesSet

  !***************************************************************

  !>
  !! Called in: src/atm_dyn_iconam/mo_nh_stepping: allocate_nh_stepping
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert (DWD) (2019-10-28)
  !!
  SUBROUTINE upatmoRestartAttributesGet( jg,           & !in
    &                                    prm_upatmo,   & !inout
    &                                    mtime_current ) !in
    
    ! In/out variables
    INTEGER,                 INTENT(IN)    :: jg
    TYPE(t_upatmo),          INTENT(INOUT) :: prm_upatmo
    TYPE(datetime), POINTER, INTENT(IN)    :: mtime_current

    ! Local variables
    TYPE(t_RestartAttributeList), POINTER :: restartAttributes
    TYPE(t_varstate_set), TARGET :: tendStateSet
    INTEGER :: i
    LOGICAL :: lmessage
    CHARACTER(LEN=2) :: domStr, iStr
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: keyStr
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':upatmoRestartAttributesGet'

    !----------------------------------------------

    IF (.NOT. ALLOCATED(prm_upatmo%tend%ddt%state)) THEN
      CALL finish(TRIM(routine), 'prm_upatmo%tend%ddt%state is not allocated')
    ENDIF

    ! Message output desired?
    lmessage = upatmo_config(jg)%l_status( iUpatmoStat%message )

    domStr = TRIM(int2string(jg, '(i2.2)'))

    IF (lmessage) CALL message(TRIM(routine), &
      & 'Start to get metadata from restart file on domain '//domStr)

    ! Event management object
    CALL upatmo_config(jg)%nwp_phy%event_mgmt_grp%deserialize(mtime_current, optAttnamePrefix=keyStrElapsedTimePhy)
    CALL upatmo_config(jg)%nwp_phy%event_mgmt_extdat%deserialize(mtime_current, optAttnamePrefix=keyStrElapsedTimeExtdat)
        
    ! State of accumulative tendencies
    restartAttributes => getAttributesForRestarting()

    ! If src/configure_model/mo_master_config: isRestart() => .FALSE., 
    ! restartAttributes should point to NULL()
    IF (ASSOCIATED(restartAttributes) .AND. jg <= domRestartLimit) THEN
      
      DO i = 1, SIZE(prm_upatmo%tend%ddt%state)
        iStr = TRIM(int2string(i, '(i2.2)'))
        ! i_old:
        keyStr = TRIM(keyStrTendStateSet_i_old)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%i_old = restartAttributes%getInteger(TRIM(keyStr))
        ! i_new:
        keyStr = TRIM(keyStrTendStateSet_i_new)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%i_new = restartAttributes%getInteger(TRIM(keyStr))
        ! n_swap:
        keyStr = TRIM(keyStrTendStateSet_n_swap)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%n_swap = restartAttributes%getInteger(TRIM(keyStr))
        ! n_state:
        keyStr = TRIM(keyStrTendStateSet_n_state)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%n_state = restartAttributes%getInteger(TRIM(keyStr))
        ! n_statep1:
        keyStr = TRIM(keyStrTendStateSet_n_statep1)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%n_statep1 = restartAttributes%getInteger(TRIM(keyStr))
        ! l_swapped:
        keyStr = TRIM(keyStrTendStateSet_l_swapped)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%l_swapped = restartAttributes%getLogical(TRIM(keyStr))
        ! l_updated:
        keyStr = TRIM(keyStrTendStateSet_l_updated)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%l_updated = restartAttributes%getLogical(TRIM(keyStr))
        ! l_locked:
        keyStr = TRIM(keyStrTendStateSet_l_locked)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%l_locked = restartAttributes%getLogical(TRIM(keyStr))
        ! l_locking:
        keyStr = TRIM(keyStrTendStateSet_l_locking)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%l_locking = restartAttributes%getLogical(TRIM(keyStr))
        ! l_unlockable:
        keyStr = TRIM(keyStrTendStateSet_l_unlockable)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%l_unlockable = restartAttributes%getLogical(TRIM(keyStr))
        ! l_final:
        keyStr = TRIM(keyStrTendStateSet_l_final)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%l_final = restartAttributes%getLogical(TRIM(keyStr))
        ! l_finish_on_error:
        keyStr = TRIM(keyStrTendStateSet_l_finish_on_error)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%l_finish_on_error = restartAttributes%getLogical(TRIM(keyStr))
        ! l_initialized:
        keyStr = TRIM(keyStrTendStateSet_l_initialized)//TRIM(domStr) &
          & //TRIM(keyStrTendStateSetSuffix)//TRIM(iStr)
        tendStateSet%l_initialized = restartAttributes%getLogical(TRIM(keyStr))
        ! 
        CALL prm_upatmo%tend%ddt%state(i)%reset(optSet4Reset=tendStateSet)
      ENDDO
      
    ENDIF

    restartAttributes => NULL()

    ! In case of jg > domRestartLimit the initial values of prm_upatmo%tend%ddt%state(i), 
    ! set in src/upper_atmosphere/mo_upatmo_state: new_upatmo_tend_list, remain.
    
    IF (lmessage) CALL message(TRIM(routine), &
      & 'Finish to get metadata from restart file on domain '//domStr)
    
  END SUBROUTINE upatmoRestartAttributesGet

  !***************************************************************

  !>
  !! Called in: src/atm_dyn_iconam/mo_nh_stepping: perform_nh_timeloop
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert (DWD) (2019-10-28)
  !!
  SUBROUTINE upatmoRestartAttributesDeallocate( upatmoRestartAttributes )

    ! In/out variables
    TYPE(t_upatmoRestartAttributes), INTENT(INOUT) :: upatmoRestartAttributes

    !----------------------------------------------

    ! elapsedTimePhy:
    CALL DO_DEALLOCATE(upatmoRestartAttributes%elapsedTimePhy)
    ! elapsedTimeExtdat:
    CALL DO_DEALLOCATE(upatmoRestartAttributes%elapsedTimeExtdat)
    ! i_old:
    CALL DO_DEALLOCATE(upatmoRestartAttributes%tendStateSet_i_old)
    ! i_new:
    CALL DO_DEALLOCATE(upatmoRestartAttributes%tendStateSet_i_new)
    ! n_swap:
    CALL DO_DEALLOCATE(upatmoRestartAttributes%tendStateSet_n_swap)
    ! n_state:
    CALL DO_DEALLOCATE(upatmoRestartAttributes%tendStateSet_n_state)
    ! n_statep1:
    CALL DO_DEALLOCATE(upatmoRestartAttributes%tendStateSet_n_statep1)
    ! l_swapped:
    CALL DO_DEALLOCATE_l1D(upatmoRestartAttributes%tendStateSet_l_swapped)
    ! l_updated:
    CALL DO_DEALLOCATE_l1D(upatmoRestartAttributes%tendStateSet_l_updated)
    ! l_locked:
    CALL DO_DEALLOCATE_l1D(upatmoRestartAttributes%tendStateSet_l_locked)
    ! l_locking:
    CALL DO_DEALLOCATE_l1D(upatmoRestartAttributes%tendStateSet_l_locking)
    ! l_unlockable:
    CALL DO_DEALLOCATE_l1D(upatmoRestartAttributes%tendStateSet_l_unlockable)
    ! l_final:
    CALL DO_DEALLOCATE_l1D(upatmoRestartAttributes%tendStateSet_l_final)
    ! l_finish_on_error:
    CALL DO_DEALLOCATE_l1D(upatmoRestartAttributes%tendStateSet_l_finish_on_error)
    ! l_initialized:
    CALL DO_DEALLOCATE_l1D(upatmoRestartAttributes%tendStateSet_l_initialized)    

  END SUBROUTINE upatmoRestartAttributesDeallocate

  !***************************************************************

  !>
  !! Auxiliary subroutine.
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert (DWD) (2019-10-28)
  !!
  SUBROUTINE setRestartAttributes_R1D( restartAttributes, attribute, key )

    ! In/out variables
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    REAL(wp), ALLOCATABLE,        INTENT(IN)    :: attribute(:)
    CHARACTER(LEN=*),             INTENT(IN)    :: key

    ! Local variables
    INTEGER :: i
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':setRestartAttributes_R1D'

    !----------------------------------------------

    IF (LEN_TRIM(key) == 0) CALL finish(TRIM(routine), 'Invalid key')
    IF (ALLOCATED(attribute)) THEN
      DO i = 1, SIZE(attribute)
        CALL restartAttributes%setReal(TRIM(key)//TRIM(int2string(i, '(i2.2)')), attribute(i))
      ENDDO
    ENDIF

  END SUBROUTINE setRestartAttributes_R1D

  !***************************************************************

  !>
  !! Auxiliary subroutine.
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert (DWD) (2019-10-28)
  !!
  SUBROUTINE setRestartAttributes_I1D( restartAttributes, attribute, key )

    ! In/out variables
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    INTEGER, ALLOCATABLE,         INTENT(IN)    :: attribute(:)
    CHARACTER(LEN=*),             INTENT(IN)    :: key

    ! Local variables
    INTEGER :: i
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':setRestartAttributes_I1D'

    !----------------------------------------------

    IF (LEN_TRIM(key) == 0) CALL finish(TRIM(routine), 'Invalid key')
    IF (ALLOCATED(attribute)) THEN
      DO i = 1, SIZE(attribute)
        CALL restartAttributes%setInteger(TRIM(key)//TRIM(int2string(i, '(i2.2)')), attribute(i))
      ENDDO
    ENDIF

  END SUBROUTINE setRestartAttributes_I1D

  !***************************************************************

  !>
  !! Auxiliary subroutine.
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert (DWD) (2019-10-28)
  !!
  SUBROUTINE setRestartAttributes_L1D( restartAttributes, attribute, key )

    ! In/out variables
    TYPE(t_RestartAttributeList), INTENT(INOUT) :: restartAttributes
    LOGICAL, ALLOCATABLE,         INTENT(IN)    :: attribute(:)
    CHARACTER(LEN=*),             INTENT(IN)    :: key

    ! Local variables
    INTEGER :: i
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':setRestartAttributes_L1D'

    !----------------------------------------------

    IF (LEN_TRIM(key) == 0) CALL finish(TRIM(routine), 'Invalid key')
    IF (ALLOCATED(attribute)) THEN
      DO i = 1, SIZE(attribute)
        CALL restartAttributes%setLogical(TRIM(key)//TRIM(int2string(i, '(i2.2)')), attribute(i))
      ENDDO
    ENDIF

  END SUBROUTINE setRestartAttributes_L1D

  !***************************************************************

  !>
  !! Auxiliary subroutine.
  !! Structural copy of src/shared/mo_fortran_tools: DO_DEALLOCATE
  !!
  !! For safety reasons we refrain from integrating this subroutine 
  !! into the interface DO_DEALLOCATE. 
  !! The latter is too critical infrastructure.
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert (DWD) (2019-10-28)
  !!
  SUBROUTINE DO_DEALLOCATE_l1D( object )

    ! In/out variables
    LOGICAL, ALLOCATABLE, INTENT(INOUT) :: object(:)

    ! Local variables
    INTEGER :: ierrstat
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':DO_DEALLOCATE_l1D'

    !----------------------------------------------

    IF (ALLOCATED(object)) THEN
      DEALLOCATE(object, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation failed!')
    END IF

  END SUBROUTINE DO_DEALLOCATE_l1D

END MODULE mo_upatmo_flowevent_utils
