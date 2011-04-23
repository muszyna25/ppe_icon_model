MODULE mo_nwp_lnd_state
!>
!!  !MODULE:  mo_nwp_phy_state\\
!!
!! Description:  Contains the data structures
!!  to store the physical model state and other auxiliary variables
!!  in order to run the NWP land physics.
!!  Constructors and destructors for these data structures as well as
!!  initialization of fields are also defined here.
!!  This module should be an analogon to 'mo_hydro_state.f90'

!!  TODO/To think about:
!     - should physics be called before or after dynamics?
!     - allocate fluxes at edges instead at the centers?
!     - horizontal/vertical tracer flux (reconstruct q'v_n' into q'u' and q'v') ?
!     - provide the "virt_inc" with meaning
!     - where to provide the lat/lon info for radiation?
!     - how to implement the echam-modules - rewriting them or "capsulate"?
!     - revision of fields if there are needed or tp be replaced
!     - fill the physics tendency construction/destruction subroutine
!     - later implement already calculated icon gradients for echam physics
!     - think about variables for flexible time steps
!!
!! @author Kristina Froehlich, DWD
!! @author Marco Giorgetta, MPI-M
!!
!! $Id: n/a$
!!
!! @par Revision History
!! Initial  by Kristina Froehlich (2010-11-09)
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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

! !USES:

USE mo_kind,                ONLY: wp
USE mo_impl_constants,      ONLY: SUCCESS
!!$USE mo_global_variables,    ONLY: l_nest_rcf
USE mo_dynamics_nml,        ONLY: nsav1, nsav2
USE mo_run_nml,             ONLY: nproma
USE mo_exception,           ONLY: message, finish
USE mo_model_domain,        ONLY: t_patch
USE mo_model_domain_import, ONLY: n_dom, l_limited_area
USE mo_atm_phy_nwp_nml,     ONLY: inwp_surface

IMPLICIT NONE
PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

!public interface
!
! subroutines
PUBLIC :: construct_nwp_lnd_state
PUBLIC :: destruct_nwp_lnd_state
!
!variables
PUBLIC :: t_lnd_state  !> state vector  for land scheme
PUBLIC :: t_lnd_prog   !!       for prognostic variables
PUBLIC :: t_lnd_diag   !!       for diagnostic variables
!




! prognostic variables state vector
  TYPE t_lnd_prog

  REAL(wp), ALLOCATABLE :: &
         t_snow       (:,:,:)     , & ! temperature of the snow-surface               (  K  )
         t_snow_mult  (:,:,:,:,:) , & ! temperature of the snow-surface               (  K  )
         t_s          (:,:,:,:)   , & ! temperature of the ground surface             (  K  )
         t_g          (:,:)       , & ! weighted surface temperature                  (  K  )
         w_snow       (:,:,:,:)   , & ! water content of snow                         (m H2O)
         rho_snow     (:,:,:)     , & ! snow density                                  (kg/m**3)
         rho_snow_mult(:,:,:,:,:)   , & ! snow density                                  (kg/m**3)
         w_i          (:,:,:,:)     , & ! water content of interception water           (m H2O)
         t_so         (:,:,:,:,:) , & ! soil temperature (main level)                 (  K  )
         w_so         (:,:,:,:,:) , & ! total water conent (ice + liquid water)       (m H20)
         w_so_ice     (:,:,:,:,:)     ! ice content                                   (m H20)

  END TYPE t_lnd_prog

  TYPE t_lnd_diag

  REAL(wp), ALLOCATABLE :: &
!
      qv_s (:,:)       , & ! specific humidity at the surface              (kg/kg)
      h_snow(:,:,:,:)    , & ! snow height                                   (  m  
      freshsnow(:,:,:) , & ! indicator for age of snow in top of snow layer(  -  )
      wliq_snow(:,:,:,:) , & ! liquid water content in the snow              (m H2O)
      wtot_snow(:,:,:,:) , & ! total (liquid + solid) water content of snow  (m H2O)
      dzh_snow(:,:,:,:)  , & ! layer thickness between half levels in snow   (  m  )
      runoff_s(:,:,:)  , & ! surface water runoff; sum over forecast      (kg/m2)
      runoff_g(:,:,:)  , & ! soil water runoff; sum over forecast         (kg/m2)
      rstom(:,:,:)     , & ! stomata resistance                           ( s/m )
      lhfl_bs(:,:,:)   , & ! average latent heat flux from bare soil evap.( W/m2)
      lhfl_pl(:,:,:)   , & ! average latent heat flux from plants         ( W/m2)
      subsfrac(:,:,:)   , & ! 
      fr_seaice(:,:)        !< fraction of sea ice                [ ]   
                            !< as partition of total area of the
                            !< grid element, but set to 0 or 1
                            !< index1=1,nproma, index2=1,nblks_c
  END TYPE t_lnd_diag


! complete state vector
  TYPE t_lnd_state

    TYPE(t_lnd_prog), ALLOCATABLE :: prog_lnd(:)

    TYPE(t_lnd_diag) :: diag_lnd

  END TYPE t_lnd_state



  CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!
  !>
  !! Constructor for prognostic and diagnostic states.
  !!
  !! It calls constructors to
  !! single time level prognostic states, and diagnostic states.
  !! Initialization of all components with zero.
  !!
  !! @par Revision History
  !! Initial release by Kristina Froehlich (2010-11-09)
  !!
  SUBROUTINE construct_nwp_lnd_state(p_patch, p_lnd_state, n_timelevels)
!
  TYPE(t_patch), TARGET, INTENT(IN)     :: p_patch(n_dom) ! patch
  INTEGER, OPTIONAL, INTENT(IN)         :: n_timelevels ! number of timelevels

  TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(n_dom)
                                           ! nh state at different grid levels
  INTEGER     :: ntl, &! local number of timelevels
                 ist, &! status
                 jg,  &! grid level counter
                 jt    ! time level counter

!-----------------------------------------------------------------------

  DO jg = 1, n_dom

    IF(PRESENT(n_timelevels))THEN
      ntl = n_timelevels
    ELSE
      ntl = 1
    ENDIF

!!$    ! If grid nesting is not called at every dynamics time step, an extra time
!!$    ! level is needed for full-field interpolation and boundary-tendency calculation
!!$    IF (l_nest_rcf .AND. n_dom > 1) THEN
!!$      ntl = ntl + 1
!!$      nsav1(jg) = ntl
!!$    ENDIF

    ! In the presence of grid nesting, another extra time level is needed to save
    ! the feedback increments
    ! This extra time level is also used to store the driving-model data in the
    ! limited-area mode
    IF (l_limited_area .OR. jg > 1) ntl = ntl + 1
    nsav2(jg) = ntl

    !create state arrays
    ALLOCATE(p_lnd_state(jg)%prog_lnd(1:ntl), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish ('mo_nonhydro_state:construct_lnd_state', &
           'allocation of land prognostic state array failed')
    ENDIF

    DO jt = 1, ntl
            !construct prognostic state
      CALL construct_lnd_state_prog(p_patch(jg),p_lnd_state(jg)%prog_lnd(jt))

    ENDDO

    !construct diagnostic state
    CALL construct_lnd_state_diag(p_patch(jg), p_lnd_state(jg)%diag_lnd)

  ENDDO !ndom

  CALL message ('mo_lnd_state:construct_lnd_state', &
                'land state construction completed')

  END SUBROUTINE construct_nwp_lnd_state

!-------------------------------------------------------------------------
!
!
  !>
  !! Destructor for prognostic and diagnostic states.
  !!
  !! It calls destructors to
  !! single time level prognostic states, and diagnostic states.
  !!
  !! @par Revision History
  !! Initial release by Kristina Froehlich (2010-11-09)
  !!
  SUBROUTINE destruct_nwp_lnd_state(p_lnd_state)
!

  TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(n_dom)
                                           ! nh state at different grid levels
  INTEGER     :: ntl, &! local number of timelevels
                 ist, &! status
                 jg,  &! grid level counter
                 jt    ! time level counter

!-----------------------------------------------------------------------

  DO jg = 1, n_dom

    ntl = SIZE(p_lnd_state(jg)%prog_lnd(:))
    IF(ntl==0)THEN
      CALL finish('mo_lnad_state:destruct_lnd_state', &
                  'prognostic array has no timelevels')
    ENDIF

    !destruct diagnostic state
!    CALL destruct_lnd_state_diag(lnd_state(jg)%diag)

    DO jt = 1, ntl

      !destruct prognostic state
      CALL destruct_lnd_state_prog(p_lnd_state(jg)%prog_lnd(jt))

    ENDDO

    DEALLOCATE(p_lnd_state(jg)%prog_lnd, STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_lnd_state:destruct_lnd_state', &
             'deallocation of land prognostic state array failed')
      ENDIF

      CALL destruct_lnd_state_diag(p_lnd_state(jg)%diag_lnd )

  ENDDO

  END SUBROUTINE destruct_nwp_lnd_state

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
  !>
  !! Allocation of components of prognostic state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Kristina Froehlich (2010-11-09)
  !!
  SUBROUTINE construct_lnd_state_prog (p_patch, p_prog_lnd)
!
  TYPE(t_patch), TARGET, INTENT(IN)     :: p_patch    ! current patch

  TYPE(t_lnd_prog), TARGET, INTENT(INOUT):: p_prog_lnd     ! current prognostic state

  INTEGER     :: nblks_c, & ! number of cell blocks to allocate
                 ist        ! status
  INTEGER     :: nztlev ! timelevel scheme
  INTEGER     :: nlevs, nlev_snow, nsfc_type    !< number of soil levels and tiles

!-----------------------------------------------------------------------

  !determine size of arrays
  nblks_c = p_patch%nblks_c


  ! T_g
  ALLOCATE(p_prog_lnd%t_g(nproma,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for weighted surface temperature T_G failed')
  ENDIF
  p_prog_lnd%t_g(:,:) =0.0_wp

  IF(inwp_surface > 0) THEN

  ! T_snow
  ALLOCATE(p_prog_lnd%t_snow(nproma,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for surface temperature T_SNOW failed')
  ENDIF
  p_prog_lnd%t_snow(:,:,:) =0.0_wp


  ! T_snow_mult
  ALLOCATE(p_prog_lnd%t_snow_mult(nproma,nlev_snow,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for surface temperature T_SNOW_MULT failed')
  ENDIF
  p_prog_lnd%t_snow_mult(:,:,:,:,:) =0.0_wp


  ! T_s
  ALLOCATE(p_prog_lnd%t_s(nproma,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for surface temperature T_S failed')
  ENDIF
  p_prog_lnd%t_s(:,:,:,:) =0.0_wp


 ! W_snow
  ALLOCATE(p_prog_lnd%w_snow(nproma,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for   water content of snow W_SNOW failed')
  ENDIF
  p_prog_lnd%w_snow(:,:,:,:) =0.0_wp

 ! Rho_snow
  ALLOCATE(p_prog_lnd%rho_snow(nproma,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for  snow density RHO_SNOW failed')
  ENDIF
  p_prog_lnd%rho_snow(:,:,:) =0.0_wp

 ! Rho_snow_mult
  ALLOCATE(p_prog_lnd%rho_snow_mult(nproma,nlev_snow,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for  snow density RHO_SNOW_MULT failed')
  ENDIF
  p_prog_lnd%rho_snow_mult(:,:,:,:,:) =0.0_wp

  ! W_I
  ALLOCATE(p_prog_lnd%w_i(nproma,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for interception store W_I failed')
  ENDIF
  p_prog_lnd%w_i(:,:,:,:) =0.0_wp

  ! T_SO
  ALLOCATE(p_prog_lnd%t_so(nproma,nlevs,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for soil temperature T_SO failed')
  ENDIF
  p_prog_lnd%t_so(:,:,:,:,:) =0.0_wp


  ! W_SO
  ALLOCATE(p_prog_lnd%w_so(nproma,nlevs,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for soil moisture W_SO failed')
  ENDIF
  p_prog_lnd%w_so(:,:,:,:,:) =0.0_wp

  ! W_SO_ICE
  ALLOCATE(p_prog_lnd%w_so_ice(nproma,nlevs,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_prog', &
                'allocation for soil ice W_SO_ICE failed')
  ENDIF
  p_prog_lnd%w_so_ice(:,:,:,:,:) =0.0_wp

ENDIF !inwp_surface

END  SUBROUTINE construct_lnd_state_prog

!-------------------------------------------------------------------------
  !>
  !! Allocation of components of diagnostic state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Kristina Froehlich (2010-11-09)
  !!
  SUBROUTINE construct_lnd_state_diag (p_patch, p_diag_lnd)
!
  TYPE(t_patch), TARGET, INTENT(IN)     :: p_patch    !  patch

  TYPE(t_lnd_diag), TARGET, INTENT(INOUT):: p_diag_lnd !diagnostic state

  INTEGER     :: nblks_c, & ! number of cell blocks to allocate
                 ist        ! status
  INTEGER     :: nztlev ! timelevel scheme
  INTEGER     :: nlevs, nsfc_type    !< number of soil levels and tiles
!-----------------------------------------------------------------------

  !determine size of arrays
  nblks_c = p_patch%nblks_c

 ! qv_s
  ALLOCATE(p_diag_lnd%qv_s(nproma,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for surface moisture QV_S failed')
  ENDIF
  p_diag_lnd%qv_s(:,:) =0.0_wp


  IF(inwp_surface > 0) THEN
  ! H_SNOW
  ALLOCATE(p_diag_lnd%h_snow(nproma,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for snow height H_SNOW failed')
  ENDIF
  p_diag_lnd%h_snow(:,:,:,:) =0.0_wp

 
  ! FRESHSNOW
  ALLOCATE(p_diag_lnd%freshsnow(nproma,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for snow age FRESHSNOW failed')
  ENDIF
  p_diag_lnd%freshsnow(:,:,:) =0.0_wp

  ! WLIQ_SNOW
  ALLOCATE(p_diag_lnd%wliq_snow(nproma,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for snow  liquid water content WLIQ_SNOW failed')
  ENDIF
  p_diag_lnd%wliq_snow(:,:,:,:) =0.0_wp


  ! WTOT_SNOW
  ALLOCATE(p_diag_lnd%wtot_snow(nproma,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for snow total  water content WTOT_SNOW failed')
  ENDIF
  p_diag_lnd%wtot_snow(:,:,:,:) =0.0_wp

  ! DZH_SNOW
  ALLOCATE(p_diag_lnd%dzh_snow(nproma,nztlev,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for snow layer thickness DZH_SNOW failed')
  ENDIF
  p_diag_lnd%dzh_snow(:,:,:,:) =0.0_wp

  ! RUNOFF_S
  ALLOCATE(p_diag_lnd%runoff_s(nproma,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for surface runoff RUNOFF_S failed')
  ENDIF
  p_diag_lnd%runoff_s(:,:,:) =0.0_wp

  ! RUNOFF_G
  ALLOCATE(p_diag_lnd%runoff_g(nproma,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for ground runoff RUNOFF_G failed')
  ENDIF
  p_diag_lnd%runoff_g(:,:,:) =0.0_wp


  ! RSTOM
  ALLOCATE(p_diag_lnd%rstom(nproma,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for  stomata resistance RSTOM  failed')
  ENDIF
  p_diag_lnd%rstom(:,:,:) =0.0_wp

  ! LHFL_BS
  ALLOCATE(p_diag_lnd%lhfl_bs(nproma,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for latent heat flux from bare soil LHFL_BS  failed')
  ENDIF
  p_diag_lnd%lhfl_bs(:,:,:) =0.0_wp


  ! LHFL_PL
  ALLOCATE(p_diag_lnd%lhfl_pl(nproma,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for average latent heat flux from plants  LHFL_PL  failed')
  ENDIF
  p_diag_lnd%lhfl_pl(:,:,:) =0.0_wp

  ! SUBSFRAC
  ALLOCATE(p_diag_lnd%subsfrac(nproma,nsfc_type,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for TILE fraction SUBSFRAC  failed')
  ENDIF
  p_diag_lnd%subsfrac(:,:,:) =0.0_wp

  ALLOCATE(p_diag_lnd%fr_seaice(nproma,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:construct_lnd_state_diag', &
                'allocation for land variables at surface FR_SEAICE failed')
  ENDIF
  p_diag_lnd%fr_seaice(:,:) = 0._wp

  ENDIF

END  SUBROUTINE construct_lnd_state_diag


!-------------------------------------------------------------------------
!
!
  !>
  !! Deallocation of components of prognostic state.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE destruct_lnd_state_prog (p_prog_lnd)
!
  TYPE(t_lnd_prog), TARGET, INTENT(INOUT) :: p_prog_lnd     ! current prognostic state

  INTEGER     :: ist        ! status

!-----------------------------------------------------------------------

  ! T_g
  DEALLOCATE(p_prog_lnd%t_g, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for weighted surface temperature T_G failed')
  ENDIF

  IF(inwp_surface > 0) THEN

  ! T_snow
  DEALLOCATE(p_prog_lnd%t_snow, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for surface temperature T_SNOW failed')
  ENDIF



  ! T_snow_mult
  DEALLOCATE(p_prog_lnd%t_snow_mult, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for surface temperature T_SNOW_MULT failed')
  ENDIF



  ! T_s
  DEALLOCATE(p_prog_lnd%t_s, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for surface temperature T_S failed')
  ENDIF



 ! W_snow
  DEALLOCATE(p_prog_lnd%w_snow, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for   water content of snow W_SNOW failed')
  ENDIF


 ! Rho_snow
  DEALLOCATE(p_prog_lnd%rho_snow, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for  snow density RHO_SNOW failed')
  ENDIF


 ! Rho_snow_mult
  DEALLOCATE(p_prog_lnd%rho_snow_mult, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for  snow density RHO_SNOW_MULT failed')
  ENDIF


  ! W_I
  DEALLOCATE(p_prog_lnd%w_i, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for interception store W_I failed')
  ENDIF


  ! T_SO
  DEALLOCATE(p_prog_lnd%t_so, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for soil temperature T_SO failed')
  ENDIF



  ! W_SO
  DEALLOCATE(p_prog_lnd%w_so, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for soil moisture W_SO failed')
  ENDIF


  ! W_SO_ICE
  DEALLOCATE(p_prog_lnd%w_so_ice, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_prog', &
                'deallocation for soil ice W_SO_ICE failed')
  ENDIF

ENDIF! surface

END  SUBROUTINE destruct_lnd_state_prog
!-------------------------------------------------------------------------
!


  !>
  !! Deallocation of components of prognostic state.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE destruct_lnd_state_diag (p_diag_lnd)
!
  TYPE(t_lnd_diag), TARGET, INTENT(INOUT) :: p_diag_lnd     ! current prognostic state

  INTEGER     :: ist        ! status

!-----------------------------------------------------------------------




 ! qv_s
  DEALLOCATE(p_diag_lnd%qv_s, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for surface moisture QV_S failed')
  ENDIF

  IF(inwp_surface > 0) THEN
  
  ! H_SNOW
  DEALLOCATE(p_diag_lnd%h_snow, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for snow height H_SNOW failed')
  ENDIF


 
  ! FRESHSNOW
  DEALLOCATE(p_diag_lnd%freshsnow, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for snow age FRESHSNOW failed')
  ENDIF


  ! WLIQ_SNOW
  DEALLOCATE(p_diag_lnd%wliq_snow, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for snow  liquid water content WLIQ_SNOW failed')
  ENDIF



  ! WTOT_SNOW
  DEALLOCATE(p_diag_lnd%wtot_snow, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for snow total  water content WTOT_SNOW failed')
  ENDIF


  ! DZH_SNOW
  DEALLOCATE(p_diag_lnd%dzh_snow, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for snow layer thickness DZH_SNOW failed')
  ENDIF


  ! RUNOFF_S
  DEALLOCATE(p_diag_lnd%runoff_s, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for surface runoff RUNOFF_S failed')
  ENDIF


  ! RUNOFF_G
  DEALLOCATE(p_diag_lnd%runoff_g, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for ground runoff RUNOFF_G failed')
  ENDIF



  ! RSTOM
  DEALLOCATE(p_diag_lnd%rstom, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for  stomata resistance RSTOM  failed')
  ENDIF


  ! LHFL_BS
  DEALLOCATE(p_diag_lnd%lhfl_bs, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for latent heat flux from bare soil LHFL_BS  failed')
  ENDIF



  ! LHFL_PL
  DEALLOCATE(p_diag_lnd%lhfl_pl, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for average latent heat flux from plants  LHFL_PL  failed')
  ENDIF
 
  ! SUBSFRAC
  DEALLOCATE(p_diag_lnd%subsfrac, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for SUBSFRAC failed')
  ENDIF

  DEALLOCATE(p_diag_lnd%fr_seaice, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_lnd_state:destruct_lnd_state_diag', &
                'deallocation for land variables at surface  failed')
  ENDIF

ENDIF
END  SUBROUTINE destruct_lnd_state_diag
!-------------------------------------------------------------------------

!
!-------------------------------------------------------------------------

END MODULE mo_nwp_lnd_state
!<
