!>
!! Type definition and variable declaration for the dynamical core of ICOHAM.
!!
!! This module contains the variables ("state vectors") used for
!! storing the prognostic and diagnostic variables of the hydrostatic
!! atmospheric dynamical core, as well as the tendency of the prognostic
!! variables.
!! Constructors and destructors for these data structures are also
!! defined here.
!!
!! This module is the substitute and extended version of the old
!! "mo_hydro_state".
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
!!     violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!     copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!      an according license agreement with DWD and MPI-M.
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
!!
MODULE mo_icoham_dyn_memory

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS
  USE mo_exception,           ONLY: message,finish
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm, t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_model_domain,        ONLY: t_patch
  USE mo_run_nml,             ONLY: ltheta_dyn
  USE mo_run_nml,             ONLY: nlev, nlevp1, ntracer, nproma

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: p_hydro_state
  PUBLIC :: construct_icoham_dyn_state, destruct_icoham_dyn_state

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_icoham_dyn_memory'

  !!----------------------------------------------------------------------------
  !! Memory buffer
  !!----------------------------------------------------------------------------

  TYPE(t_hydro_atm),TARGET,ALLOCATABLE :: p_hydro_state(:) !< state vector on
                                                           !< different grid levels
                                                           !< shape: (n_dom)

CONTAINS

  !!----------------------------------------------------------------------------
  !! Subroutines for allocating/deallocating memory
  !!----------------------------------------------------------------------------
  !>
  !! Subroutine that allocates memory for the state vector on ALL grid levels
  !!
  SUBROUTINE construct_icoham_dyn_state( ntimelevel, p_patch )

    INTEGER,INTENT(IN)     :: ntimelevel
    TYPE(t_patch),INTENT(IN) :: p_patch (:)

    INTEGER :: ndomain
    INTEGER :: jg    !< grid level/domain index
    INTEGER :: jt    !< time level index
    INTEGER :: ist   !< system status code
    !---

    CALL message(TRIM(thismodule),'Construction of 3D dynamics state vector started.')

    ndomain = SIZE(p_patch)

    DO jg = 1,ndomain

      !----------------------------
      ! 1.  For time integration:
      !----------------------------
      ! 1.1 Prognostic variables

      ALLOCATE(p_hydro_state(jg)%prog(1:ntimelevel), STAT=ist)
      IF (ist/=SUCCESS) &
      CALL finish(TRIM(thismodule),'allocation of prognostic state array failed')

      DO jt = 1,ntimelevel
        CALL construct_hydro_state_prog( p_patch(jg), p_hydro_state(jg)%prog(jt) )
      END DO

      ! 1.2 Diagnostic variables
      CALL construct_hydro_state_diag( p_patch(jg), p_hydro_state(jg)%diag )

      ! 1.3 Tendencies
      CALL construct_hydro_state_prog( p_patch(jg), p_hydro_state(jg)%tend_dyn )
      CALL construct_hydro_state_prog( p_patch(jg), p_hydro_state(jg)%tend_phy )

      !----------------------------
      ! 2.  For organizing output
      !----------------------------
      CALL construct_hydro_state_prog( p_patch(jg), p_hydro_state(jg)%prog_out )
      CALL construct_hydro_state_diag( p_patch(jg), p_hydro_state(jg)%diag_out )

    ENDDO
    CALL message(TRIM(thismodule),'Construction of 3D dynamics state vector finished.')

  END SUBROUTINE construct_icoham_dyn_state

  !>
  !! Subroutine that released memory used by the state vector
  !!
  SUBROUTINE destruct_icoham_dyn_state

    INTEGER :: ntimelevel, ndomain
    INTEGER :: jt    !< time level index
    INTEGER :: jg    !< grid level/domain index
    INTEGER :: ist   !< system status code
    !---

    CALL message(TRIM(thismodule),'Destruction of 3D dynamics state vector started.')

    ndomain    = SIZE(p_hydro_state)
    ntimelevel = SIZE(p_hydro_state(1)%prog)

    DO jg = 1,ndomain

      ! Prognostic variables
      DO jt = 1,ntimelevel
        CALL destruct_hydro_state_prog( p_hydro_state(jg)%prog(jt) )
      END DO

      DEALLOCATE( p_hydro_state(jg)%prog, STAT=ist )
      IF (ist/=SUCCESS) &
      CALL finish(TRIM(thismodule),'deallocation of prognostic state array failed')

      ! Diagnostic variables
      CALL destruct_hydro_state_diag( p_hydro_state(jg)%diag )

      ! Tendencies
      CALL destruct_hydro_state_prog( p_hydro_state(jg)%tend_dyn )
      CALL destruct_hydro_state_prog( p_hydro_state(jg)%tend_phy )

      ! Memory used for organizing output
      CALL destruct_hydro_state_prog( p_hydro_state(jg)%prog_out )
      CALL destruct_hydro_state_diag( p_hydro_state(jg)%diag_out )

    ENDDO
    CALL message(TRIM(thismodule),'Destruction of 3D dynamics state vector finished.')

  END SUBROUTINE destruct_icoham_dyn_state

  !>
  !!----------------------------------------------------------------
  !! For variables of type "t_icoham_dyn_prog"
  !!----------------------------------------------------------------
  !!
  SUBROUTINE construct_hydro_state_prog(p_patch, p_prog)

    TYPE(t_patch),           INTENT(IN) :: p_patch
    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: p_prog

    INTEGER :: nblks_c, nblks_e
    INTEGER :: ist

    ! Size of the horizontal dimension

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! Start allocation

    ALLOCATE( p_prog% pres_sfc(nproma,     nblks_c),          &
              p_prog% vn      (nproma,nlev,nblks_e),          &
              p_prog% temp    (nproma,nlev,nblks_c), STAT=ist )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation failed.')

    IF (ltheta_dyn) THEN
      ALLOCATE( p_prog% theta(nproma,nlev,nblks_c), STAT=ist )
      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation of theta failed.')
    ENDIF

    IF (ntracer > 0) THEN
      ALLOCATE( p_prog% tracer(nproma,nlev,nblks_c,ntracer),STAT=ist )
      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation of tracer field failed.')
    ENDIF

    ! Initialize all fields with zero

    p_prog% pres_sfc = 0.0_wp
    p_prog% vn       = 0.0_wp
    p_prog% temp     = 0.0_wp
    p_prog% vn       = 0.0_wp

    IF (ltheta_dyn)  p_prog% theta  = 0.0_wp
    IF (ntracer > 0) p_prog% tracer = 0.0_wp

  END SUBROUTINE construct_hydro_state_prog

  !>
  !!
  SUBROUTINE destruct_hydro_state_prog( p_prog )

    TYPE(t_hydro_atm_prog),TARGET,INTENT(INOUT) :: p_prog

    INTEGER :: ist

    DEALLOCATE( p_prog% pres_sfc, &
                p_prog% vn      , &
                p_prog% temp    , STAT=ist )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'Deallocation failed.')

    IF (ltheta_dyn) THEN
      DEALLOCATE( p_prog% theta, STAT=ist )
      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'deallocation of theta failed.')
    ENDIF

    IF (ntracer > 0) THEN
      DEALLOCATE( p_prog% tracer, STAT=ist )
      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'deallocation of tracer failed.')
    ENDIF

  END SUBROUTINE destruct_hydro_state_prog

  !>
  !!----------------------------------------------------------------
  !! For variables of type "t_icoham_dyn_diag"
  !!----------------------------------------------------------------
  !!
  SUBROUTINE construct_hydro_state_diag( p_patch, p_diag )

    TYPE(t_patch),           INTENT(IN)    :: p_patch
    TYPE(t_hydro_atm_diag),INTENT(INOUT) :: p_diag

    INTEGER :: nblks_c, nblks_e, nblks_v
    INTEGER :: ist

    ! Size of the horizontal dimension

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

   ! Allocate memory

    ALLOCATE( p_diag% qx          (nproma,nlev  ,nblks_c), &
            & p_diag% u           (nproma,nlev  ,nblks_c), &
            & p_diag% v           (nproma,nlev  ,nblks_c), &
            & p_diag% vt          (nproma,nlev  ,nblks_e), &
            & p_diag% rel_vort    (nproma,nlev  ,nblks_v), &
            & p_diag% rel_vort_e  (nproma,nlev  ,nblks_e), &
            & p_diag% rel_vort_c  (nproma,nlev  ,nblks_c), &
            & p_diag% div         (nproma,nlev  ,nblks_c), &
            & p_diag% e_kin       (nproma,nlev  ,nblks_c), &
            & p_diag% geo_ic      (nproma,nlevp1,nblks_c), &
            & p_diag% geo_mc      (nproma,nlev  ,nblks_c), &
            & p_diag% wpres_mc    (nproma,nlev  ,nblks_c), &
            & p_diag% wpres_ic    (nproma,nlevp1,nblks_c), &
            & p_diag% weta        (nproma,nlevp1,nblks_c), &
            & p_diag% pres_ic     (nproma,nlevp1,nblks_c), &
            & p_diag% pres_ic_new (nproma,nlevp1,nblks_c), &
            & p_diag% pres_mc     (nproma,nlev  ,nblks_c), &
            & p_diag% delp_c      (nproma,nlev  ,nblks_c), &
            & p_diag% delp_c_new  (nproma,nlev  ,nblks_c), &
            & p_diag% rdelp_c     (nproma,nlev  ,nblks_c), &
            & p_diag% rdelp_c_new (nproma,nlev  ,nblks_c), &
            & p_diag% delp_e      (nproma,nlev  ,nblks_e), &
            & p_diag% delp_v      (nproma,nlev  ,nblks_v), &
            & p_diag% virt_incr   (nproma,nlev  ,nblks_c), &
            & p_diag% tempv       (nproma,nlev  ,nblks_c), &
            & p_diag% rdlnpr_c    (nproma,nlev  ,nblks_c), &
            & p_diag% rdalpha_c   (nproma,nlev  ,nblks_c), &
            & p_diag% lnp_ic      (nproma,nlevp1,nblks_c), &
            & p_diag% mass_flux_e (nproma,nlev  ,nblks_e), &
            & STAT=ist )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation of diag state failed.')

    IF (ltheta_dyn) THEN
      ALLOCATE( p_diag% exner (nproma,nlev,nblks_c), STAT=ist )
      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation of exner failed.')
    ENDIF

    IF (ntracer > 0) THEN

      ALLOCATE( p_diag% hfl_tracer (nproma,nlev  ,nblks_e,ntracer), &
              & p_diag% vfl_tracer (nproma,nlevp1,nblks_c,ntracer), &
              & STAT=ist )

      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'allocation of tracer flux failed.')
    ENDIF

    ! Initialize all components with zero

    p_diag% qx          = 0.0_wp
    p_diag% u           = 0.0_wp
    p_diag% v           = 0.0_wp
    p_diag% vt          = 0.0_wp
    p_diag% rel_vort    = 0.0_wp
    p_diag% rel_vort_e  = 0.0_wp
    p_diag% rel_vort_c  = 0.0_wp
    p_diag% div         = 0.0_wp
    p_diag% e_kin       = 0.0_wp
    p_diag% geo_ic      = 0.0_wp
    p_diag% geo_mc      = 0.0_wp
    p_diag% wpres_mc    = 0.0_wp
    p_diag% wpres_ic    = 0.0_wp
    p_diag% weta        = 0.0_wp
    p_diag% pres_ic     = 0.0_wp
    p_diag% pres_ic_new = 0.0_wp
    p_diag% pres_mc     = 0.0_wp
    p_diag% delp_c      = 0.0_wp
    p_diag% delp_c_new  = 0.0_wp
    p_diag% rdelp_c     = 0.0_wp
    p_diag% rdelp_c_new = 0.0_wp
    p_diag% delp_e      = 0.0_wp
    p_diag% delp_v      = 0.0_wp
    p_diag% virt_incr   = 0.0_wp
    p_diag% tempv       = 0.0_wp
    p_diag% rdlnpr_c    = 0.0_wp
    p_diag% rdalpha_c   = 0.0_wp
    p_diag% lnp_ic      = 0.0_wp
    p_diag% mass_flux_e = 0.0_wp

    IF (ltheta_dyn) p_diag% exner =0.0_wp

    IF (ntracer > 0) THEN
      p_diag% hfl_tracer = 0.0_wp
      p_diag% vfl_tracer = 0.0_wp
    END IF

  END SUBROUTINE construct_hydro_state_diag
  !>
  !!
  !!
  SUBROUTINE destruct_hydro_state_diag( p_diag )

    TYPE(t_hydro_atm_diag),TARGET,INTENT(INOUT) :: p_diag
    INTEGER :: ist

   ! Deallocate memory

    DEALLOCATE( p_diag% qx          , &
              & p_diag% u           , &
              & p_diag% v           , &
              & p_diag% vt          , &
              & p_diag% rel_vort    , &
              & p_diag% rel_vort_e  , &
              & p_diag% rel_vort_c  , &
              & p_diag% div         , &
              & p_diag% e_kin       , &
              & p_diag% geo_ic      , &
              & p_diag% geo_mc      , &
              & p_diag% wpres_mc    , &
              & p_diag% wpres_ic    , &
              & p_diag% weta        , &
              & p_diag% pres_ic     , &
              & p_diag% pres_ic_new , &
              & p_diag% pres_mc     , &
              & p_diag% delp_c      , &
              & p_diag% delp_c_new  , &
              & p_diag% rdelp_c     , &
              & p_diag% rdelp_c_new , &
              & p_diag% delp_e      , &
              & p_diag% delp_v      , &
              & p_diag% virt_incr   , &
              & p_diag% tempv       , &
              & p_diag% rdlnpr_c    , &
              & p_diag% rdalpha_c   , &
              & p_diag% lnp_ic      , &
              & p_diag% mass_flux_e , &
              & STAT=ist )

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'deallocation of diag state failed.')

    IF (ltheta_dyn) THEN
      DEALLOCATE( p_diag% exner, STAT=ist )
      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'deallocation of exner failed.')
    ENDIF
    
    IF (ntracer > 0) THEN
      DEALLOCATE( p_diag% hfl_tracer, p_diag% vfl_tracer, STAT=ist )
      IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),'deallocation of tracer flux failed.')
    ENDIF

  END SUBROUTINE destruct_hydro_state_diag

END module mo_icoham_dyn_memory
