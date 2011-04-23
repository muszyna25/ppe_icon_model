!>
!! Defines, constructs and destructs the state vector of the nonhydrostatic.
!!
!! Defines, constructs and destructs the state vector of the nonhydrostatic
!! model variables. They are subdivided in several classes: prognostics
!! and diagnostics.
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2009-03-06)
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
!!
MODULE mo_nonhydro_state
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS
  USE mo_exception,           ONLY: message,message_text, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_model_domain_import, ONLY: n_dom, l_limited_area
  USE mo_nonhydrostatic_nml,  ONLY: l_nest_rcf
  USE mo_dynamics_nml,        ONLY: nsav1, nsav2, itime_scheme
  USE mo_run_nml,             ONLY: nproma, i_cell_type, iforcing,               &
    &                               inwp, ltransport, ntracer, ntracer_static,   &
    &                               inextra_2d, inextra_3d
  USE mo_io_nml,              ONLY: lwrite_extra
!  USE mo_atm_phy_nwp_nml,     ONLY: inwp_turb

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: t_nh_prog               ! state vector of prognostic variables
  PUBLIC :: t_nh_diag               ! state vector of diagnostic variables
  PUBLIC :: t_nh_metrics            ! state vector of metrics variables
  PUBLIC :: t_nh_state              ! state vector of nonhydrostatic variables

  PUBLIC :: construct_nh_state    ! Constructor for the nonhydrostatic state
  PUBLIC :: destruct_nh_state     ! Destructor for the nonhydrostatic state

  PUBLIC :: t_buffer_memory, bufr

  ! prognostic variables state vector
  TYPE t_nh_prog

    REAL(wp), ALLOCATABLE :: &
      w(:,:,:),          & !> orthogonal vertical wind (nproma,nlevp1,nblks_c)     [m/s]
      vn(:,:,:),         & !! orthogonal normal wind (nproma,nlev,nblks_e)         [m/s]
      rho(:,:,:),        & !! density (nproma,nlev,nblks_c)                     [kg/m^3]
      exner(:,:,:),      & !! Exner pressure (nproma,nlev,nblks_c)                   [-]
      theta_v(:,:,:),    & !! virtual potential temperature (nproma,nlev,nblks_c)    [K]
      rhotheta_v(:,:,:), & !! rho*theta_v (nproma,nlev,nblks_c)               [K*kg/m^3]
      tracer(:,:,:,:),   & !! tracer concentration (nproma,nlev,nblks_c,ntracer) [kg/kg]
      tke   (:,:,:)        !! SQRT(2 * turbulent kinetik energy)                 [ m/s ]
                           !! (defined on half levels) with 2 time levels

  END TYPE t_nh_prog

  ! diagnostic variables state vector
  TYPE t_nh_diag

    REAL(wp), ALLOCATABLE :: &
    ! a) variables needed for triangles and hexagons
    &  u(:,:,:),            & ! zonal wind (nproma,nlev,nblks_c)               [m/s]
    &  v(:,:,:),            & ! meridional wind (nproma,nlev,nblks_c)          [m/s]
    &  vt(:,:,:),           & ! tangential wind (nproma,nlev,nblks_e)          [m/s]
    &  e_kin(:,:,:),        & ! spec. kinetic energy (nproma,nlev,nblks_c) [m^2/s^2]
    &  omega_z(:,:,:),      & ! vertical vorticity at dual grid
                              ! (nproma,nlev,nblks_v or nblks_e)               [1/s]
    &  ddt_vn(:,:,:),       & ! normal wind tendency from forcing
    &  ddt_vn_phy(:,:,:),       & ! normal wind tendency from forcing
                              ! (nproma,nlev,nblks_e)   [m/s^2]
    &  ddt_w(:,:,:),        & ! vert. wind tendency from forcing
                              ! (nproma,nlevp1,nblks_c)  [m/s^2]
    &  ddt_exner(:,:,:),    & ! exner pressure tendency from forcing (nproma,nlev,nblks_c)  [1/s]
    &  ddt_exner_phy(:,:,:),& ! exner pressure tendency from physical forcing (nproma,nlev,nblks_c)  [1/s]
    &  ddt_tracer_phy(:,:,:,:), &! physics tendency of tracers
                              ! (nproma,nlev,nblks_c,ntracer)             [kg/kg/s] !T.R.
    &  ddt_tracer_adv(:,:,:,:), &! advective tendency of tracers          [kg/kg/s]
    &  exner_old(:,:,:),    & ! exner pres from previous step (nproma,nlev,nblks_c)
                            ! *** needs to be saved for restart ***
    &  w_con(:,:,:),        & ! contravariant vert wind (nproma,nlevp1,nblks_c)[1/s]
    &  temp(:,:,:),         & ! temperature (nproma,nlev,nblks_c)                 [K]
    &  tempv(:,:,:),        & ! virtual temperature (nproma,nlev,nblks_c)         [K]
    &  temp_ifc(:,:,:),     & ! temperature at half levels (nproma,nlevp1,nblks_c)[K]
    &  pres(:,:,:),         & ! pressure (nproma,nlev,nblks_c)                  [Pa]
    &  pres_ifc(:,:,:),     & ! pressure at interfaces (nproma,nlevp1,nblks_c)  [Pa]
    &  pres_sfc(:,:),       & ! diagnosed surface pressure (nproma,nblks_c)     [Pa]
    &  dpres_mc(:,:,:),     & ! pressure thickness at masspoints(nproma,nlevp,nblks_c)  [Pa]
    &  hfl_tracer(:,:,:,:), & ! horizontal tracer flux at edges             [kg/m/s]
                              ! (nproma,nlev,nblks_e,ntracer)
    &  vfl_tracer(:,:,:,:), & ! vertical tracer flux at cells               [kg/m/s]
                              ! (nproma,nlevp1,nblks_c,ntracer)
    &  div(:,:,:),          & ! divergence(nproma,nlev,nblks_c)     [1/s]
    &  mass_fl_e(:,:,:),    & ! horizontal mass flux at edges (nproma,nlev,nblks_e) [kg/m/s]
    &  rho_ic(:,:,:),       & ! density at half levels (nproma,nlevp1,nblks_c)     [kg/m^3]
    &  theta_v_ic(:,:,:),   & ! theta_v at half levels (nproma,nlevp1,nblks_c)         [K]
    &  w_concorr_c(:,:,:),  & ! contravariant vert correction (nproma,nlevp1,nblks_c)[m/s]
                            ! *** needs to be saved for restart for triangular grid ***
    &  e_kinh(:,:,:),       & ! horizontal spec. kinetic energy
                              ! (nproma,nlev,nblks_c)                            [m^2/s^2]
                              !  
    ! b) variables needed for the triangular grid only
    &  vn_ie(:,:,:),        & ! normal wind at half levels (nproma,nlevp1,nblks_e)   [m/s]
    &  ddt_vn_adv(:,:,:,:), & ! normal wind tendency from advection
                              ! (nproma,nlev,nblks_e,1:3)                    [m/s^2]
    &  ddt_w_adv(:,:,:,:),  & ! vert. wind tendency from advection
                              ! (nproma,nlevp1,nblks_c,1:3)                  [m/s^2]
    &  grf_tend_vn(:,:,:),  & ! vn tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_e)                        [m/s^2]
    &  grf_tend_w(:,:,:),   & ! w tendency field for use in grid refinement
                              ! (nproma,nlevp1,nblks_c)                      [m/s^2]
    &  grf_tend_rho(:,:,:), & ! rho tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_c)                     [kg/m^3/s]
    &  grf_tend_thv(:,:,:), & ! theta_v tendency field for use in grid refinement
                              ! (nproma,nlev,nblks_c)                          [K/s]
    &  grf_tend_tracer(:,:,:,:), & ! tracer tendency field for use in grid refinement
                                   ! (nproma,nlev,nblks_c,ntracer)          [kg/kg/s]
    &  dvn_ie_int(:,:),    & ! Storage field for vertical nesting: vn at parent interface level
    &  dvn_ie_ubc(:,:),    & ! Storage field for vertical nesting: vn at child upper boundary
    &  drho_ic_int(:,:),   & ! Storage field for vertical nesting: rho at parent interface level
    &  drho_ic_ubc(:,:),   & ! Storage field for vertical nesting: rho at child upper boundary
    &  dtheta_v_ic_int(:,:),& ! Storage field for vertical nesting: theta at parent interface level
    &  dtheta_v_ic_ubc(:,:),& ! Storage field for vertical nesting: theta at child upper boundary
    &  dw_int(:,:),        & ! Storage field for vertical nesting: w at parent interface level
    &  dw_ubc(:,:),        & ! Storage field for vertical nesting: w at child upper boundary
    &  q_int(:,:,:),       & ! Storage field for vertical nesting: q at parent interface level
    &  q_ubc(:,:,:),       & ! Storage field for vertical nesting: q at child upper boundary
    &  thermal_exp_fastphy(:,:), & ! Contribution of fast physics processes to thermal expansion
                                   ! (needed for open upper boundary condition)
    !
    ! c) variables needed for the hexagonal grid only
    &  theta_v_impl(:,:,:), & ! (nnow+nnew)/2 from impl. vert. adv. of theta_v   [K]
    &  horpgrad(:,:,:),     & ! covariant horizontal pressur gradient term   [m/s^2]
    &  vn_cov(:,:,:),       & ! covariant normal wind (nproma,nlev,nblks_e)    [m/s]
    &  w_cov(:,:,:),        & ! covariant vert wind (nproma,nlevp1,nblks_c)  [m^2/s]
    &  omega_z_con(:,:,:),  & ! vertical contrav. vorticity (nproma,nlev,nblks_v)
                              ! at vertices                                    [1/s]
    &  omega_t_con(:,:,:),  & ! tangential horiz. contravariant vorticity
                              ! (nproma,nlev,nblks_e)[1/s]
    &  omega_x(:,:,:),      & ! zonal vorticity (nproma,nlev,nblks_c)          [1/s]
    &  omega_y(:,:,:),      & ! meridional vorticity (nproma,nlev,nblks_c)     [1/s]
    &  ddt_vn_vort(:,:,:),  & ! normal wind tendency from vorticity flux term
                              ! (nproma,nlev,nblks_e,1:3)                    [m/s^2]
    &  ddt_w_vort(:,:,:)      ! vert. wind tendency from vorticity flux term
                              ! (nproma,nlevp1,nblks_c,1:3)                  [m/s^2]

    REAL(wp), ALLOCATABLE :: & !
     &  extra_2d(:,:,:)  ,   & !> extra debug output in 2d and
     &  extra_3d(:,:,:,:)      !!                       3d

    REAL(wp), POINTER :: &
      vn_con(:,:,:),     &! contravariant normal wind (nproma,nlev,nblks_e)[m/s]
      omega_t(:,:,:)      ! tangent. horiz. vorticity (nproma,nlev,nblks_e)[1/s]

  END TYPE t_nh_diag

  TYPE t_nh_metrics

   ! a) Variables needed for triangles and hexagons
   !
   ! geometric height at the vertical interface of cells (nproma,nlevp1,nblks_c)
   REAL(wp), ALLOCATABLE :: z_ifc(:,:,:)
   ! geometric height at full levels (nproma,nlev,nblks_c)
   REAL(wp), ALLOCATABLE :: z_mc(:,:,:)
   ! geometric height at full level edges (nproma,nlev,nblks_e)
   REAL(wp), ALLOCATABLE :: z_mc_e(:,:,:)

   ! slope of the terrain in normal direction (nproma,nlevp1,nblks_e)
   REAL(wp), ALLOCATABLE :: ddxn_z_half(:,:,:)
   ! slope of the terrain in normal direction (nproma,nlev,nblks_e)
   REAL(wp), ALLOCATABLE :: ddxn_z_full(:,:,:)
   ! slope of the terrain in tangential direction (nproma,nlevp1,nblks_e)
   REAL(wp), ALLOCATABLE :: ddxt_z_half(:,:,:)

   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_c)
   REAL(wp), ALLOCATABLE :: ddqz_z_full(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_e)
   REAL(wp), ALLOCATABLE :: ddqz_z_full_e(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlevp1,nblks_c)
   REAL(wp), ALLOCATABLE :: ddqz_z_half(:,:,:)

   ! 1/dz(k-1)-1/dz(k) (nproma,nlevp1,nblks_c)
   REAL(wp), ALLOCATABLE :: diff_1_o_dz(:,:,:)
   ! 1/(dz(k-1)*dz(k)) (nproma,nlevp1,nblks_c)
   REAL(wp), ALLOCATABLE :: mult_1_o_dz(:,:,:)

   ! geopotential at cell center (nproma,nlev,nblks_c)
   REAL(wp), ALLOCATABLE :: geopot(:,:,:)
   ! geopotential at interfaces at cell center (nproma,nlevp1,nblks_c)
   REAL(wp), ALLOCATABLE :: geopot_ifc(:,:,:)
   ! geopotential above groundlevel at cell center (nproma,nlev,nblks_c)
   REAL(wp), ALLOCATABLE :: geopot_agl(:,:,:)
   ! geopotential above groundlevel at interfaces and cell center (nproma,nlevp1,nblks_c)
   REAL(wp), ALLOCATABLE :: geopot_agl_ifc(:,:,:)
   ! geopotential at cell center (nproma,nlev,nblks_c)
   REAL(wp), ALLOCATABLE :: dgeopot_mc(:,:,:)


   ! Rayleigh damping on the vertical velocity
   REAL(wp), ALLOCATABLE :: rayleigh_w(:,:,:)

   ! b) Variables needed for the triangular grid only
   !
   ! weighting factor for interpolation from full to half levels (nproma,nlevp1,nblks_c)
   REAL(wp), ALLOCATABLE :: wgtfac_c(:,:,:)
   ! weighting factor for interpolation from full to half levels (nproma,nlevp1,nblks_e)
   REAL(wp), ALLOCATABLE :: wgtfac_e(:,:,:)
   ! weighting factor for quadratic interpolation to surface (nproma,3,nblks_c)
   REAL(wp), ALLOCATABLE :: wgtfacq_c(:,:,:)
   ! weighting factor for quadratic interpolation to surface (nproma,3,nblks_e)
   REAL(wp), ALLOCATABLE :: wgtfacq_e(:,:,:)
   ! weighting factor for quadratic interpolation to model top (nproma,3,nblks_c)
   REAL(wp), ALLOCATABLE :: wgtfacq1_c(:,:,:)
   ! weighting factor for quadratic interpolation to model top (nproma,3,nblks_e)
   REAL(wp), ALLOCATABLE :: wgtfacq1_e(:,:,:)
   ! Inverse layer thickness of full levels (nproma,nlev,nblks_c)
   REAL(wp), ALLOCATABLE :: inv_ddqz_z_full(:,:,:)
   ! Inverse distance between full levels jk+1 and jk-1 (nproma,nlev,nblks_c)
   REAL(wp), ALLOCATABLE :: inv_ddqz_z_half2(:,:,:)
   ! Vertical index of neighbor points needed for Taylor-expansion-based pressure gradient
   INTEGER,  ALLOCATABLE :: vertidx_gradp(:,:,:,:) ! (2,nproma,nlev,nblks_e)
   ! Height differences between local edge point and neighbor cell points used for
   ! pressure gradient computation (2,nproma,nlev,nblks_e)
   REAL(wp), ALLOCATABLE :: zdiff_gradp(:,:,:,:)
   ! extrapolation factor for Exner pressure (slope-dependent for stability optimization) (nproma,nlev,nblks_c)
   REAL(wp), ALLOCATABLE :: exner_exfac(:,:,:)
   ! explicit weight in vertical wind solver (nproma,nblks_c)
   REAL(wp), ALLOCATABLE :: vwind_expl_wgt(:,:)
   ! implicit weight in vertical wind solver (nproma,nblks_c)
   REAL(wp), ALLOCATABLE :: vwind_impl_wgt(:,:)
   ! Fields for reference atmosphere
   REAL(wp), ALLOCATABLE :: theta_ref_mc(:,:,:)
   REAL(wp), ALLOCATABLE :: theta_ref_ic(:,:,:)
   REAL(wp), ALLOCATABLE :: exner_ref_mc(:,:,:)
   REAL(wp), ALLOCATABLE :: d_exner_dz_ref_ic(:,:,:)
   REAL(wp), ALLOCATABLE :: d_exner_dz_ref_mc(:,:,:)
   REAL(wp), ALLOCATABLE :: d2_exner_dz2_ref_mc(:,:,:)
   REAL(wp), ALLOCATABLE :: rho_refcorr_ic(:,:,:)
   ! Fields for truly horizontal temperature diffusion
   INTEGER               :: zd_listdim
   INTEGER,  ALLOCATABLE :: zd_indlist(:,:)
   INTEGER,  ALLOCATABLE :: zd_blklist(:,:)
   INTEGER,  ALLOCATABLE :: zd_edgeidx(:,:)
   INTEGER,  ALLOCATABLE :: zd_edgeblk(:,:)
   INTEGER,  ALLOCATABLE :: zd_vertidx(:,:)
   REAL(wp), ALLOCATABLE :: zd_intcoef(:,:)
   REAL(wp), ALLOCATABLE :: zd_geofac(:,:)
   REAL(wp), ALLOCATABLE :: zd_e2cell(:,:)
   REAL(wp), ALLOCATABLE :: zd_diffcoef(:)
   ! Fields for igradp_method = 3
   INTEGER               :: pg_listdim
   INTEGER,  ALLOCATABLE :: pg_edgeidx(:)
   INTEGER,  ALLOCATABLE :: pg_edgeblk(:)
   INTEGER,  ALLOCATABLE :: pg_vertidx(:)
   REAL(wp), ALLOCATABLE :: pg_exdist (:)
   ! Index lists for grid points on which lateral boundary nudging is applied
   INTEGER               :: nudge_c_dim, nudge_e_dim
   INTEGER,  ALLOCATABLE :: nudge_c_idx(:)
   INTEGER,  ALLOCATABLE :: nudge_e_idx(:)
   INTEGER,  ALLOCATABLE :: nudge_c_blk(:)
   INTEGER,  ALLOCATABLE :: nudge_e_blk(:)
   ! Index lists and mask fields needed to minimize the number of halo communications
   ! a) index lists for halo points belonging to the nest boundary region
   INTEGER               :: bdy_halo_c_dim  
   INTEGER,  ALLOCATABLE :: bdy_halo_c_idx(:)
   INTEGER,  ALLOCATABLE :: bdy_halo_c_blk(:)
   ! b) a mask field that excludes boundary halo points
   LOGICAL,  ALLOCATABLE :: mask_prog_halo_c(:,:)
   ! c) index lists for halo points belonging to a nest overlap region
   !    the additional dimension is n_childdom
   INTEGER,  ALLOCATABLE :: ovlp_halo_c_dim(:)
   INTEGER,  ALLOCATABLE :: ovlp_halo_c_idx(:,:)
   INTEGER,  ALLOCATABLE :: ovlp_halo_c_blk(:,:)


   ! c) Variables needed for the hexagonal grid only
   !
   ! slope of the coordinate lines towards the North (nproma,nlev,nblks_e)
   REAL(wp), ALLOCATABLE :: ddnorth_z(:,:,:)
   ! slope of the coordinate lines towards the East (nproma,nlev,nblks_e)
   REAL(wp), ALLOCATABLE :: ddeast_z(:,:,:)

   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_v)
   REAL(wp), ALLOCATABLE :: ddqz_z_full_v(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlevp1,nblks_e)
   REAL(wp), ALLOCATABLE :: ddqz_z_half_e(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlevp1,nblks_e)
   REAL(wp), ALLOCATABLE :: ddqz_z_half_r(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlevp1,nblks_v)
   REAL(wp), ALLOCATABLE :: ddqz_z_half_v(:,:,:)

  END TYPE t_nh_metrics

  ! complete state vector
  TYPE t_nh_state

    !array of prognostic states at different timelevels
    TYPE(t_nh_prog),     ALLOCATABLE :: prog(:)

    TYPE(t_nh_diag)              :: diag

    TYPE(t_nh_metrics)           :: metrics

  END TYPE t_nh_state

  TYPE :: t_buffer_memory

    REAL(wp), ALLOCATABLE :: send_c1(:,:)
    REAL(wp), ALLOCATABLE :: recv_c1(:,:)
    REAL(wp), ALLOCATABLE :: send_c3(:,:)
    REAL(wp), ALLOCATABLE :: recv_c3(:,:)
    REAL(wp), ALLOCATABLE :: send_e1(:,:)
    REAL(wp), ALLOCATABLE :: recv_e1(:,:)
    REAL(wp), ALLOCATABLE :: send_e2(:,:)
    REAL(wp), ALLOCATABLE :: recv_e2(:,:)
    REAL(wp), ALLOCATABLE :: send_e3(:,:)
    REAL(wp), ALLOCATABLE :: recv_e3(:,:)
    REAL(wp), ALLOCATABLE :: send_v2(:,:)
    REAL(wp), ALLOCATABLE :: recv_v2(:,:)

  END TYPE t_buffer_memory

  TYPE (t_buffer_memory), ALLOCATABLE :: bufr(:)

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
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE construct_nh_state(p_patch, p_nh_state, n_timelevels)
!
  TYPE(t_patch), TARGET, INTENT(IN)       :: p_patch(n_dom) ! patch
  INTEGER, OPTIONAL, INTENT(IN)         :: n_timelevels ! number of timelevels

  TYPE(t_nh_state), TARGET, INTENT(INOUT) :: p_nh_state(n_dom)
                                           ! nh state at different grid levels
  INTEGER     :: ntl, &! local number of timelevels
                 ist, &! status
                 jg,  &! grid level counter
                 jt    ! time level counter

  LOGICAL     :: l_alloc_tracer

!-----------------------------------------------------------------------

  DO jg = 1, n_dom

    IF(PRESENT(n_timelevels))THEN
      ntl = n_timelevels
    ELSE
      ntl = 1
    ENDIF

    ! If grid nesting is not called at every dynamics time step, an extra time
    ! level is needed for full-field interpolation and boundary-tendency calculation
    IF (l_nest_rcf .AND. n_dom > 1) THEN
      ntl = ntl + 1
      nsav1(jg) = ntl
    ENDIF

    ! In the presence of grid nesting, another extra time level is needed to save
    ! the feedback increments
    ! This extra time level is also used to store the driving-model data in the
    ! limited-area mode
    IF (l_limited_area .OR. jg > 1) ntl = ntl + 1
    nsav2(jg) = ntl

    !create state arrays
    ALLOCATE(p_nh_state(jg)%prog(1:ntl), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish ('mo_nonhydro_state:construct_nh_state', &
                   'allocation of prognostic state array failed')
    ENDIF

    DO jt = 1, ntl

      ! Tracer fields do not need extra time levels because feedback is not incremental
      ! and the nest-call frequency is always synchronized with the advection time step
      IF (jt <= n_timelevels) THEN
        l_alloc_tracer = .TRUE.
      ELSE
        l_alloc_tracer = .FALSE.
      ENDIF

      !construct prognostic state
      CALL construct_nh_state_prog(p_patch(jg), p_nh_state(jg)%prog(jt), l_alloc_tracer)

    ENDDO

    !construct diagnostic state
    CALL construct_nh_state_diag(p_patch(jg), p_nh_state(jg)%diag)

    !construct metrics state
    CALL construct_nh_metrics(p_patch(jg), p_nh_state(jg)%metrics)

  ENDDO

  CALL message ('mo_nonhydro_state:construct_nh_state', &
                'NH state construction completed')

  END SUBROUTINE construct_nh_state

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
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE destruct_nh_state(p_nh_state)
!

  TYPE(t_nh_state), TARGET, INTENT(INOUT) :: p_nh_state(n_dom)
                                           ! nh state at different grid levels
  INTEGER     :: ntl, &! local number of timelevels
                 ist, &! status
                 jg,  &! grid level counter
                 jt    ! time level counter

!-----------------------------------------------------------------------

  DO jg = 1, n_dom

    ntl = SIZE(p_nh_state(jg)%prog(:))
    IF(ntl==0)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state', &
                  'prognostic array has no timelevels')
    ENDIF

    !destruct diagnostic state
    CALL destruct_nh_state_diag(p_nh_state(jg)%diag)

    !destruct metrics state
    CALL destruct_nh_metrics (p_nh_state(jg)%metrics)

    DO jt = 1, ntl

      !destruct prognostic state
      CALL destruct_nh_state_prog(p_nh_state(jg)%prog(jt))

    ENDDO

    !destruct state array
    DEALLOCATE(p_nh_state(jg)%prog, STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish ('mo_nonhydro_state:destruct_nh_state', &
                   'deallocation of prognostic state array failed')
    ENDIF

  ENDDO

  END SUBROUTINE destruct_nh_state

!-------------------------------------------------------------------------
!
!
  !>
  !! Allocation of components of prognostic state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE construct_nh_state_prog (p_patch, p_prog, l_alloc_tracer)
!
  TYPE(t_patch), TARGET, INTENT(IN)     :: p_patch    ! current patch

  TYPE(t_nh_prog), TARGET, INTENT(INOUT):: p_prog     ! current prognostic state

  LOGICAL, INTENT(IN) :: l_alloc_tracer  ! allocate tracer fields if true

  INTEGER :: nblks_c, &    ! number of cell blocks to allocate
             nblks_e, &    ! number of edge blocks to allocate
             ist           ! status
  INTEGER :: nlev, nlevp1  ! number of full and half levels

!-----------------------------------------------------------------------

  !determine size of arrays
  nblks_c = p_patch%nblks_c
  nblks_e = p_patch%nblks_e

  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  ! vn
  ALLOCATE(p_prog%vn(nproma,nlev,nblks_e), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_prog', &
                'allocation for normal wind failed')
  ENDIF

  ! w
  ALLOCATE(p_prog%w(nproma,nlevp1,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_prog', &
                'allocation for vertical wind failed')
  ENDIF

  ! rho
  ALLOCATE(p_prog%rho(nproma,nlev,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_prog', &
                'allocation for density failed')
  ENDIF

  ! exner
  ALLOCATE(p_prog%exner(nproma,nlev,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_prog', &
                'allocation for Exner pressure failed')
  ENDIF

  ! theta_v
  ALLOCATE(p_prog%theta_v(nproma,nlev,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_prog', &
                'allocation for theta_v failed')
  ENDIF

  ! rhotheta_v
  ALLOCATE(p_prog%rhotheta_v(nproma,nlev,nblks_c), STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_prog', &
                'allocation for rho*theta_v failed')
  ENDIF

   !tracer
  IF ( (ltransport .OR. iforcing == inwp) .AND. l_alloc_tracer  ) THEN
     ALLOCATE(p_prog%tracer(nproma,nlev,nblks_c,ntracer+ntracer_static),STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_prog', &
       'allocation for tracer failed')
    ENDIF
  ENDIF





  !tke
  IF ( iforcing == inwp ) THEN
    ALLOCATE(p_prog%tke(nproma,nlevp1,nblks_c),STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_prog', &
       'allocation for tke failed')
    ENDIF
  ENDIF

  p_prog%vn(:,:,:)         = 0.0_wp
  p_prog%w(:,:,:)          = 0.0_wp
  p_prog%rho(:,:,:)        = 0.0_wp
  p_prog%exner(:,:,:)      = 0.0_wp
  p_prog%rhotheta_v(:,:,:) = 0.0_wp
  p_prog%theta_v(:,:,:)    = 0.0_wp
  IF ( (ltransport .OR. iforcing == inwp ) .AND. l_alloc_tracer ) THEN
    p_prog%tracer(:,:,:,:) = 0.0_wp
  ENDIF
 IF ( iforcing == inwp ) THEN
   p_prog%tke(:,:,:) = 1.e-4_wp
 ENDIF

  END SUBROUTINE construct_nh_state_prog

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
  SUBROUTINE destruct_nh_state_prog (p_prog)
!
  TYPE(t_nh_prog), TARGET, INTENT(INOUT) :: p_prog     ! current prognostic state

  INTEGER     :: ist        ! status

!-----------------------------------------------------------------------

  ! vn
  DEALLOCATE(p_prog%vn, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_prog', &
                'deallocation for normal wind failed')
  ENDIF

  ! w
  DEALLOCATE(p_prog%w, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_prog', &
                'deallocation for vertical wind failed')
  ENDIF

  ! rho
  DEALLOCATE(p_prog%rho, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_prog', &
                'deallocation for density failed')
  ENDIF

  ! exner
  DEALLOCATE(p_prog%exner, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_prog', &
                'deallocation for Exner pressure failed')
  ENDIF

  ! theta_v
  DEALLOCATE(p_prog%theta_v, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_prog', &
                'deallocation for theta_v failed')
  ENDIF

  ! rhotheta_v
  DEALLOCATE(p_prog%rhotheta_v, STAT = ist)
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_prog', &
                'deallocation for rho*theta_v failed')
  ENDIF

  IF ( (ltransport.OR. iforcing == inwp ) .AND. ALLOCATED(p_prog%tracer) ) THEN
    DEALLOCATE( p_prog%tracer, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_prog', &
        'deallocation for tracer failed')
    ENDIF
  ENDIF

  IF ( iforcing == 3 .AND. ALLOCATED(p_prog%tke) ) THEN
    DEALLOCATE( p_prog%tke, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_prog', &
        'deallocation for tke failed')
    ENDIF
  ENDIF

  END SUBROUTINE destruct_nh_state_prog

!
  !>
  !! Allocation of components of diagnostic state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-03-06)
  !!
  !! Modification by Kristina Fröhlich, DWD, (2010-10-22)
  !! - added pressure on interfaces

  SUBROUTINE construct_nh_state_diag (p_patch, p_diag)
!
  TYPE(t_patch), TARGET, INTENT(IN)     :: p_patch    ! current patch

  TYPE(t_nh_diag), TARGET, INTENT(INOUT):: p_diag     ! current prognostic state

  INTEGER :: nblks_c, &    ! number of cell blocks to allocate
             nblks_e, &    ! number of edge blocks to allocate
             nblks_v, &    ! number of vertex blocks to allocate
             ist           ! status
  INTEGER :: nlev, nlevp1  ! number of full and half levels

  INTEGER :: n_timlevs     ! number of time levels for advection tendency fields

!-----------------------------------------------------------------------

  !determine size of arrays
  nblks_c = p_patch%nblks_c
  nblks_e = p_patch%nblks_e
  nblks_v = p_patch%nblks_v

  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  IF (itime_scheme == 6) THEN
   n_timlevs = 2
  ELSE
   n_timlevs = 1
  ENDIF

  ! u
  ALLOCATE(p_diag%u(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for zonal wind failed')
  ENDIF

  ! v
  ALLOCATE(p_diag%v(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for meridional wind failed')
  ENDIF

  ! vt
  ALLOCATE(p_diag%vt(nproma,nlev,nblks_e), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for vt failed')
  ENDIF

  ! e_kin
  ALLOCATE(p_diag%e_kin(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for specific kinetic energy failed')
  ENDIF

  ! omega_z
  ALLOCATE(p_diag%omega_z(nproma,nlev,nblks_v), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for vertical voritcity failed')
  ENDIF

  ! ddt_vn
  ALLOCATE(p_diag%ddt_vn(nproma,nlev,nblks_e), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for normal wind tendency failed')
  ENDIF

  ! ddt_vn_phy
  ALLOCATE(p_diag%ddt_vn_phy(nproma,nlev,nblks_e), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for normal wind physical tendency failed')
  ENDIF

  ! ddt_w
  ALLOCATE(p_diag%ddt_w(nproma,nlevp1,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for vertical wind tendency failed')
  ENDIF

  ! ddt_exner
  ALLOCATE(p_diag%ddt_exner(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for exner pressure tendency failed')
  ENDIF

  ! ddt_exner_phy
  ALLOCATE(p_diag%ddt_exner_phy(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for exner pressure physical tendency failed')
  ENDIF


  ! exner_old
  ALLOCATE(p_diag%exner_old(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for old exner pressure failed')
  ENDIF

  ! w_con
  ALLOCATE(p_diag%w_con(nproma,nlevp1,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for w_con failed')
  ENDIF

  ! pres_sfc
  ALLOCATE(p_diag%pres_sfc(nproma,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for surface pressure failed')
  ENDIF

  ! temp
  ALLOCATE(p_diag%temp(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for temperature failed')
  ENDIF

  ! tempv
  ALLOCATE(p_diag%tempv(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for virtual temperature failed')
  ENDIF

  ! temp_ifc
  ALLOCATE(p_diag%temp_ifc(nproma,nlevp1,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for temperature failed')
  ENDIF

  ! pres
  ALLOCATE(p_diag%pres(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for pressure failed')
  ENDIF

  ! pres_ifc
  ALLOCATE(p_diag%pres_ifc(nproma,nlevp1,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for pressure on interfaces failed')
  ENDIF

  ! dpres_mc
  ALLOCATE(p_diag%dpres_mc(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for pressure thickness failed')
  ENDIF

  ! div
  ALLOCATE(p_diag%div(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for divergence')
  ENDIF

  ! mass_fl_e
  ALLOCATE(p_diag%mass_fl_e(nproma,nlev,nblks_e), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for mass_fl_e failed')
  ENDIF

  ! rho_ic
  ALLOCATE(p_diag%rho_ic(nproma,nlevp1,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for rho_ic failed')
  ENDIF

  ! w_concorr_c
  ALLOCATE(p_diag%w_concorr_c(nproma,nlevp1,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for w_concorr_c failed')
  ENDIF

  ! e_kinh
  ALLOCATE(p_diag%e_kinh(nproma,nlev,nblks_c), STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                'allocation for horizontal specific kinetic energy failed')
  ENDIF

  IF (i_cell_type == 3) THEN

    ! vn_ie
    ALLOCATE(p_diag%vn_ie(nproma,nlevp1,nblks_e), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for vn_ie failed')
    ENDIF

    ! theta_v_ic
    ALLOCATE(p_diag%theta_v_ic(nproma,nlevp1,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for theta_v_ic failed')
    ENDIF

    ! ddt_vn_adv
    ALLOCATE(p_diag%ddt_vn_adv(nproma,nlev,nblks_e,n_timlevs), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for normal wind tendency failed')
    ENDIF
    ! ddt_w_adv
    ALLOCATE(p_diag%ddt_w_adv(nproma,nlevp1,nblks_c,n_timlevs), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for vertical wind tendency failed')
    ENDIF

    ! grf_tend_vn
    ALLOCATE(p_diag%grf_tend_vn(nproma,nlev,nblks_e), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for grf_tend_vn failed')
    ENDIF

    ! grf_tend_w
    ALLOCATE(p_diag%grf_tend_w(nproma,nlevp1,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for grf_tend_w failed')
    ENDIF

    ! grf_tend_rho
    ALLOCATE(p_diag%grf_tend_rho(nproma,nlev,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for grf_tend_rho failed')
    ENDIF

    ! grf_tend_thv
    ALLOCATE(p_diag%grf_tend_thv(nproma,nlev,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for grf_tend_thv failed')
    ENDIF

    ! Storage fields for vertical nesting; the middle index (2) addresses the field and its temporal tendency

    ! dvn_ie_int
    ALLOCATE(p_diag%dvn_ie_int(nproma,nblks_e), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for dvn_ie_int failed')
    ENDIF

    ! dvn_ie_ubc
    ALLOCATE(p_diag%dvn_ie_ubc(nproma,nblks_e), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for dvn_ie_ubc failed')
    ENDIF

    ! drho_ic_int
    ALLOCATE(p_diag%drho_ic_int(nproma,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for drho_ic_int failed')
    ENDIF

    ! drho_ic_ubc
    ALLOCATE(p_diag%drho_ic_ubc(nproma,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for drho_ic_ubc failed')
    ENDIF

    ! dtheta_v_ic_int
    ALLOCATE(p_diag%dtheta_v_ic_int(nproma,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for dtheta_v_ic_int failed')
    ENDIF

    ! dtheta_v_ic_ubc
    ALLOCATE(p_diag%dtheta_v_ic_ubc(nproma,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for dtheta_v_ic_ubc failed')
    ENDIF

    ! dw_int
    ALLOCATE(p_diag%dw_int(nproma,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for dw_int failed')
    ENDIF

    ! dw_ubc
    ALLOCATE(p_diag%dw_ubc(nproma,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for dw_ubc failed')
    ENDIF

    ! q_int
    ALLOCATE(p_diag%q_int(nproma,nblks_c,ntracer), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for q_int failed')
    ENDIF

    ! q_ubc
    ALLOCATE(p_diag%q_ubc(nproma,nblks_c,ntracer), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for q_ubc failed')
    ENDIF

    ! thermal_exp_fastphy
    ALLOCATE(p_diag%thermal_exp_fastphy(nproma,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for thermal_exp_fastphy failed')
    ENDIF

  ELSE IF (i_cell_type == 6) THEN

    ! theta_v_ic
    ALLOCATE(p_diag%theta_v_ic(nproma,nlevp1,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for theta_v_ic failed')
    ENDIF

    ! theta_v_impl
    ALLOCATE(p_diag%theta_v_impl(nproma,nlev,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for theta_v_impl failed')
    ENDIF

    ! horpgrad
    ALLOCATE(p_diag%horpgrad(nproma,nlev,nblks_e), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for horpgrad failed')
    ENDIF

    ! vn_cov
    ALLOCATE(p_diag%vn_cov(nproma,nlev,nblks_e), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for vn_cov failed')
    ENDIF

    ! w_cov
    ALLOCATE(p_diag%w_cov(nproma,nlevp1,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for w_cov failed')
    ENDIF

    ! omega_t_con
    ALLOCATE(p_diag%omega_t_con(nproma,nlevp1,nblks_e), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for tangential contrav. voritcity failed')
    ENDIF

    !omega_x
    ALLOCATE(p_diag%omega_x(nproma,nlev,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for zonal vorticity failed')
    ENDIF

    !omega_y
    ALLOCATE(p_diag%omega_y(nproma,nlev,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for meridional vorticity failed')
    ENDIF

    ! omega_z_con
    ALLOCATE(p_diag%omega_z_con(nproma,nlev,nblks_v), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for vertical contrav. voritcity failed')
    ENDIF

    ! ddt_vn_vort
    ALLOCATE(p_diag%ddt_vn_vort(nproma,nlev,nblks_e), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for normal wind tendency failed')
    ENDIF

    ! ddt_w_vort
    ALLOCATE(p_diag%ddt_w_vort(nproma,nlevp1,nblks_c), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for vertical wind tendency failed')
    ENDIF

  ENDIF

  !tracers
  IF ( ltransport ) THEN
    ! grf_tend_tracer
    ALLOCATE(p_diag%grf_tend_tracer(nproma,nlev,nblks_c,ntracer), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for grf_tend_tracer failed')
    ENDIF

    ! hfl_tracer
    ALLOCATE(p_diag%hfl_tracer(nproma,nlev,nblks_e,ntracer), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for horizontal tracer flux at edges failed')
    ENDIF

    ! vfl_tracer
    ALLOCATE(p_diag%vfl_tracer(nproma,nlevp1,nblks_c,ntracer), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for vertical tracer flux at half levels failed')
    ENDIF

    ! ddt_tracer_adv
    ALLOCATE(p_diag%ddt_tracer_adv(nproma,nlev,nblks_c,ntracer), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
                  'allocation for ddt_tracer_adv failed')
    ENDIF

  ENDIF


  IF( iforcing== inwp) THEN  !T.R
    ALLOCATE(p_diag%ddt_tracer_phy(nproma,nlev,nblks_c,ntracer), STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
           'allocation of p_diag%ddt_tracer_phy failed')
    ENDIF
  ENDIF

  IF( lwrite_extra) THEN
    WRITE(0,*)'inextra_2d=',inextra_2d
  
    IF(inextra_2d > 0) &
    ALLOCATE(p_diag%extra_2d(nproma,     nblks_c,inextra_2d), STAT=ist )
    IF(inextra_3d > 0) &
    ALLOCATE(p_diag%extra_3d(nproma,nlev,nblks_c,inextra_3d), STAT=ist )

    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
           'allocation of extra fields for debugging failed')
    ENDIF
  ENDIF


  p_diag%u(:,:,:)             = 0.0_wp
  p_diag%v(:,:,:)             = 0.0_wp

  p_diag%vt(:,:,:)            = 0.0_wp
  p_diag%e_kin(:,:,:)         = 0.0_wp
  p_diag%omega_z(:,:,:)       = 0.0_wp
  p_diag%ddt_vn(:,:,:)        = 0.0_wp
  p_diag%ddt_w(:,:,:)         = 0.0_wp
  p_diag%ddt_exner(:,:,:)     = 0.0_wp
  p_diag%exner_old(:,:,:)     = 0.0_wp
  p_diag%w_con(:,:,:)         = 0.0_wp
  p_diag%temp(:,:,:)          = 0.0_wp
  p_diag%tempv(:,:,:)         = 0.0_wp
  p_diag%temp_ifc(:,:,:)      = 0.0_wp
  p_diag%pres(:,:,:)          = 0.0_wp
  p_diag%pres_ifc(:,:,:)      = 0.0_wp
  p_diag%dpres_mc(:,:,:)      = 0.0_wp
  p_diag%pres_sfc(:,:)        = 0.0_wp
  p_diag%div(:,:,:)           = 0.0_wp
  p_diag%mass_fl_e(:,:,:)     = 0.0_wp
  p_diag%rho_ic(:,:,:)        = 0.0_wp
  p_diag%w_concorr_c(:,:,:)   = 0.0_wp
  p_diag%e_kinh(:,:,:)        = 0.0_wp
  IF (i_cell_type == 3) THEN
    p_diag%vn_ie(:,:,:)         = 0.0_wp
    p_diag%theta_v_ic(:,:,:)    = 0.0_wp
    p_diag%ddt_vn_adv(:,:,:,:)  = 0.0_wp
    p_diag%ddt_w_adv(:,:,:,:)   = 0.0_wp
    p_diag%grf_tend_vn(:,:,:)   = 0.0_wp
    p_diag%grf_tend_w(:,:,:)    = 0.0_wp
    p_diag%grf_tend_rho(:,:,:)  = 0.0_wp
    p_diag%grf_tend_thv(:,:,:)  = 0.0_wp
    p_diag%dvn_ie_int(:,:)      = 0.0_wp
    p_diag%dvn_ie_ubc(:,:)      = 0.0_wp
    p_diag%drho_ic_int(:,:)     = 0.0_wp
    p_diag%drho_ic_ubc(:,:)     = 0.0_wp
    p_diag%dtheta_v_ic_int(:,:) = 0.0_wp
    p_diag%dtheta_v_ic_ubc(:,:) = 0.0_wp
    p_diag%dw_int(:,:)          = 0.0_wp
    p_diag%dw_ubc(:,:)          = 0.0_wp
    p_diag%q_int(:,:,:)         = 0.0_wp
    p_diag%q_ubc(:,:,:)         = 0.0_wp
    p_diag%thermal_exp_fastphy(:,:) = 0.0_wp
  ELSE IF (i_cell_type == 6) THEN
    p_diag%theta_v_impl(:,:,:)  = 0.0_wp
    p_diag%theta_v_ic(:,:,:)    = 0.0_wp
    p_diag%horpgrad(:,:,:)      = 0.0_wp
    p_diag%vn_cov(:,:,:)        = 0.0_wp
    p_diag%w_cov(:,:,:)         = 0.0_wp
    p_diag%omega_x(:,:,:)       = 0.0_wp
    p_diag%omega_y(:,:,:)       = 0.0_wp
    p_diag%omega_t_con(:,:,:)   = 0.0_wp
    p_diag%omega_z_con(:,:,:)   = 0.0_wp
    p_diag%ddt_vn_vort(:,:,:)   = 0.0_wp
    p_diag%ddt_w_vort(:,:,:)    = 0.0_wp
  ENDIF
  IF ( ltransport ) THEN
    p_diag%ddt_tracer_adv(:,:,:,:) = 0.0_wp
    p_diag%grf_tend_tracer(:,:,:,:)= 0.0_wp
    p_diag%hfl_tracer(:,:,:,:)     = 0.0_wp
    p_diag%vfl_tracer(:,:,:,:)     = 0.0_wp
  ENDIF
  IF ( iforcing == inwp) THEN
    p_diag%ddt_tracer_phy(:,:,:,:)= 0.0_wp
  ENDIF
  p_diag%ddt_vn_phy    (:,:,:)  = 0.0_wp
  p_diag%ddt_exner_phy (:,:,:)  = 0.0_wp

  IF ( lwrite_extra) THEN
    IF(inextra_2d > 0) &
    p_diag%extra_2d   (:,:,:)   = 0.0_wp
    IF(inextra_3d > 0) &
    p_diag%extra_3d   (:,:,:,:) = 0.0_wp
  ENDIF
  END SUBROUTINE construct_nh_state_diag
!-------------------------------------------------------------------------
!
!
  !>
  !! Deallocation of components of diagnostic state.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE destruct_nh_state_diag (p_diag)
!
  TYPE(t_nh_diag), TARGET, INTENT(INOUT) :: p_diag     ! current diagnostic state

  INTEGER     :: ist        ! status

!-----------------------------------------------------------------------

  ! u
  DEALLOCATE(p_diag%u, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for zonal wind failed')
  ENDIF

  ! v
  DEALLOCATE(p_diag%v, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for meridional wind failed')
  ENDIF

  ! vt
  DEALLOCATE(p_diag%vt, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for vt failed')
  ENDIF

  ! e_kin
  DEALLOCATE(p_diag%e_kin, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for specific kinetic energy failed')
  ENDIF

  ! omega_z
  DEALLOCATE(p_diag%omega_z, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for vertical voritcity failed')
  ENDIF

  ! ddt_vn
  DEALLOCATE(p_diag%ddt_vn, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for normal wind tendency failed')
  ENDIF

  ! ddt_w
  DEALLOCATE(p_diag%ddt_w, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for vertical wind tendency failed')
  ENDIF

  ! ddt_exner
  DEALLOCATE(p_diag%ddt_exner, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for exner pressure tendency failed')
  ENDIF

  ! exner_old
  DEALLOCATE(p_diag%exner_old, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for old exner pressure failed')
  ENDIF

  ! w_con
  DEALLOCATE(p_diag%w_con, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for w_con failed')
  ENDIF

  ! pres_sfc
  DEALLOCATE(p_diag%pres_sfc, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for surface pressure failed')
  ENDIF

  ! temp
  DEALLOCATE(p_diag%temp, STAT=ist )
  ! tempv
  DEALLOCATE(p_diag%tempv, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for temperature failed')
  ENDIF

  ! temp_ifc
  DEALLOCATE(p_diag%temp_ifc, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for temperature on interfaces failed')
  ENDIF

  ! pres
  DEALLOCATE(p_diag%pres, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for pressure failed')
  ENDIF

  ! pres_ifc
  DEALLOCATE(p_diag%pres_ifc, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for pressure failed')
  ENDIF

  ! dpres_mc
  DEALLOCATE(p_diag%dpres_mc, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for pressure thickness failed')
  ENDIF


  ! div
  DEALLOCATE(p_diag%div, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for divergence')
  ENDIF

  ! mass_fl_e
  DEALLOCATE(p_diag%mass_fl_e, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for mass_fl_e failed')
  ENDIF

  ! rho_ic
  DEALLOCATE(p_diag%rho_ic, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for rho_ic failed')
  ENDIF

  ! w_concorr_c
  DEALLOCATE(p_diag%w_concorr_c, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for w_concorr_c failed')
  ENDIF

  ! e_kinh
  DEALLOCATE(p_diag%e_kinh, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                'deallocation for horizontal specific kinetic energy failed')
  ENDIF

  IF (i_cell_type == 3) THEN

    ! vn_ie
    DEALLOCATE(p_diag%vn_ie, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for vn_ie failed')
    ENDIF

    ! theta_v_ic
    DEALLOCATE(p_diag%theta_v_ic, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for theta_v_ic failed')
    ENDIF

    ! ddt_vn_adv
    DEALLOCATE(p_diag%ddt_vn_adv, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for normal wind tendency failed')
    ENDIF

    ! ddt_w_adv
    DEALLOCATE(p_diag%ddt_w_adv, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for vertical wind tendency failed')
    ENDIF

    ! grf_tend_vn
    DEALLOCATE(p_diag%grf_tend_vn, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for grf_tend_vn failed')
    ENDIF

    ! grf_tend_w
    DEALLOCATE(p_diag%grf_tend_w, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for grf_tend_w failed')
    ENDIF

    ! grf_tend_rho
    DEALLOCATE(p_diag%grf_tend_rho, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for grf_tend_rho failed')
    ENDIF

    ! grf_tend_thv
    DEALLOCATE(p_diag%grf_tend_thv, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for grf_tend_thv failed')
    ENDIF

    ! dvn_ie_int
    DEALLOCATE(p_diag%dvn_ie_int, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for dvn_ie_int failed')
    ENDIF

    ! dvn_ie_ubc
    DEALLOCATE(p_diag%dvn_ie_ubc, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for dvn_ie_ubc failed')
    ENDIF

    ! drho_ic_int
    DEALLOCATE(p_diag%drho_ic_int, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for drho_ic_int failed')
    ENDIF

    ! drho_ic_ubc
    DEALLOCATE(p_diag%drho_ic_ubc, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for drho_ic_ubc failed')
    ENDIF

    ! dtheta_v_ic_int
    DEALLOCATE(p_diag%dtheta_v_ic_int, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for dtheta_v_ic_int failed')
    ENDIF

    ! dtheta_v_ic_ubc
    DEALLOCATE(p_diag%dtheta_v_ic_ubc, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for dtheta_v_ic_ubc failed')
    ENDIF

    ! dw_int
    DEALLOCATE(p_diag%dw_int, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for dw_int failed')
    ENDIF

    ! dw_ubc
    DEALLOCATE(p_diag%dw_ubc, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for dw_ubc failed')
    ENDIF

    ! q_int
    DEALLOCATE(p_diag%q_int, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for q_int failed')
    ENDIF

    ! q_ubc
    DEALLOCATE(p_diag%q_ubc, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for q_ubc failed')
    ENDIF

    ! thermal_exp_fastphy
    DEALLOCATE(p_diag%thermal_exp_fastphy, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for thermal_exp_fastphy failed')
    ENDIF

  ELSE IF (i_cell_type == 6) THEN

    ! theta_v_ic
    DEALLOCATE(p_diag%theta_v_ic, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for theta_v_ic failed')
    ENDIF

    ! theta_v_impl
    DEALLOCATE(p_diag%theta_v_impl, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for theta_v_impl failed')
    ENDIF

    ! horpgrad
    DEALLOCATE(p_diag%horpgrad, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for horpgrad failed')
    ENDIF

    ! vn_cov
    DEALLOCATE(p_diag%vn_cov, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for vn_cov failed')
    ENDIF

    ! w_cov
    DEALLOCATE(p_diag%w_cov, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for w_cov failed')
    ENDIF

    ! omega_t_con
    DEALLOCATE(p_diag%omega_t_con, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for tangential contrav. voritcity failed')
    ENDIF

    !omega_x
    DEALLOCATE(p_diag%omega_x, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for zonal vorticity failed')
    ENDIF

    !omega_y
    DEALLOCATE(p_diag%omega_y, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for meridional vorticity failed')
    ENDIF

    ! omega_z_con
    DEALLOCATE(p_diag%omega_z_con, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for vertical contrav. voritcity failed')
    ENDIF

    ! ddt_vn_vort
    DEALLOCATE(p_diag%ddt_vn_vort, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for normal wind tendency failed')
    ENDIF

    ! ddt_w_vort
    DEALLOCATE(p_diag%ddt_w_vort, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for vertical wind tendency failed')
    ENDIF

  ENDIF

  !tracers
  IF ( ltransport ) THEN
    ! grf_tend_tracer
    DEALLOCATE(p_diag%grf_tend_tracer, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for grf_tend_tracer failed')
    ENDIF

    ! hfl_tracer
    DEALLOCATE(p_diag%hfl_tracer, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for horizontal tracer flux at edges failed')
    ENDIF

    ! vfl_tracer
    DEALLOCATE(p_diag%vfl_tracer, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for vertical tracer flux at half levels failed')
    ENDIF

    ! ddt_tracer_adv
    DEALLOCATE(p_diag%ddt_tracer_adv, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:destruct_nh_state_diag', &
                  'deallocation for ddt_tracer_adv failed')
    ENDIF

  ENDIF


  IF( iforcing== inwp) THEN  !T.R
    DEALLOCATE(p_diag%ddt_tracer_phy, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
           'deallocation of p_diag%ddt_tracer_phy failed')
    ENDIF
  ENDIF 
  DEALLOCATE(p_diag%ddt_vn_phy, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
         'deallocation of p_diag%ddt_vn_phy failed')
  ENDIF
  DEALLOCATE(p_diag%ddt_exner_phy, STAT=ist )
  IF (ist/=SUCCESS)THEN
    CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
         'deallocation of p_diag%ddt_exner_phy failed')
  ENDIF

  IF( lwrite_extra) THEN  
    IF(inextra_2d > 0) &
    DEALLOCATE(p_diag%extra_2d, STAT=ist )
    IF(inextra_3d > 0) &
    DEALLOCATE(p_diag%extra_3d, STAT=ist )
    IF (ist/=SUCCESS)THEN
      CALL finish('mo_nonhydro_state:construct_nh_state_diag', &
           'deallocation of debugging extra fields failed')
    ENDIF
  ENDIF 


  END SUBROUTINE destruct_nh_state_diag

  !---------------------------------------------------------------------------
  !>
  !! Allocates all metric coefficients defined in type metrics_3d of the patch.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-04-14)
  !! Modification by Daniel Reinert, DWD, (2010-04-22)
  !! - added geometric height at full levels


  SUBROUTINE construct_nh_metrics(p_patch, p_metrics)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch !< patch

    TYPE(t_nh_metrics), TARGET, INTENT(INOUT):: p_metrics ! metrics state

    INTEGER :: nblks_c, nblks_e, nblks_v, ist
    INTEGER :: nlev, nlevp1        ! number of full levels
    !------------------------------------------------------------------

      nblks_c = p_patch%nblks_c
      nblks_e = p_patch%nblks_e
      nblks_v = p_patch%nblks_v

      ! number of vertical levels
      nlev   = p_patch%nlev
      nlevp1 = p_patch%nlevp1

      ! geometric height at the vertical interface of cells
      ALLOCATE(p_metrics%z_ifc(nproma,nlevp1,nblks_c), STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for z_ifc failed')
      ENDIF

      ! geometric height at full levels
      ALLOCATE(p_metrics%z_mc(nproma,nlev,nblks_c), STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for z_mc failed')
      ENDIF

      ! geometric height at full level edges
      ALLOCATE(p_metrics%z_mc_e(nproma,nlev,nblks_e), STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for z_mc_e failed')
      ENDIF

      ! slope of the terrain in normal direction
      ALLOCATE(p_metrics%ddxn_z_half(nproma,nlevp1,nblks_e),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for ddxn_z_half failed')
      ENDIF
      ! slope of the terrain
      ALLOCATE(p_metrics%ddxn_z_full(nproma,nlev,nblks_e),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for ddxn_z_full failed')
      ENDIF
      ! slope of the terrain in tangential direction
      ALLOCATE(p_metrics%ddxt_z_half(nproma,nlevp1,nblks_e),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for ddxt_z_half failed')
      ENDIF

      ! functional determinant of the metrics [sqrt(gamma)]
      ALLOCATE(p_metrics%ddqz_z_full(nproma,nlev,nblks_c), STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for ddqz_z_full failed')
      ENDIF
      ! functional determinant of the metrics [sqrt(gamma)]
      ALLOCATE(p_metrics%ddqz_z_full_e(nproma,nlev,nblks_e),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for ddqz_z_full_e failed')
      ENDIF
      ! functional determinant of the metrics [sqrt(gamma)]
      ALLOCATE(p_metrics%ddqz_z_half(nproma,nlevp1,nblks_c),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for ddqz_z_half failed')
      ENDIF

      ! 1/dz(k-1)-1/dz(k) 
      ALLOCATE(p_metrics%diff_1_o_dz(nproma,nlevp1,nblks_c),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for diff_1_o_dz failed')
      ENDIF
      ! 1/(dz(k-1)*dz(k))
      ALLOCATE(p_metrics%mult_1_o_dz(nproma,nlevp1,nblks_c),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for mult_1_o_dz failed')
      ENDIF

      ! geopotential at cell center
      ALLOCATE(p_metrics%geopot(nproma,nlev,nblks_c), STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for geopot failed')
      ENDIF

      ! geopotential on interfaces at cell center
      ALLOCATE(p_metrics%geopot_ifc(nproma,nlevp1,nblks_c), STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for geopot_ifc failed')
      ENDIF

      ! geopotential above groundlevel at cell center
      ALLOCATE(p_metrics%geopot_agl    (nproma,nlev  ,nblks_c), STAT = ist)
      ALLOCATE(p_metrics%geopot_agl_ifc(nproma,nlevp1,nblks_c), STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for geopot_agl failed')
      ENDIF

      ! geopotential at cell center
      ALLOCATE(p_metrics%dgeopot_mc(nproma,nlev,nblks_c), STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for dgeopot_mc failed')
      ENDIF


      ! Rayleigh damping
      ALLOCATE(p_metrics%rayleigh_w(nproma,nlevp1,nblks_c),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for rayleigh_w failed')
      ENDIF

      ! Explicit weight in vertical wind solver
      ALLOCATE(p_metrics%vwind_expl_wgt(nproma,nblks_c),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for vwind_expl_wgt failed')
      ENDIF

      ! Implicit weight in vertical wind solver
      ALLOCATE(p_metrics%vwind_impl_wgt(nproma,nblks_c),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for vwind_impl_wgt failed')
      ENDIF

! These fields are needed for triangles only once the initialization in
! mo_nh_testcases is properly rewritten for hexagons
!      IF (i_cell_type== 3) THEN
        ! weighting factor for interpolation from full to half levels
        ALLOCATE(p_metrics%wgtfac_c(nproma,nlevp1,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for wgtfac_c failed')
        ENDIF

       ! weighting factor for interpolation from full to half levels
        ALLOCATE(p_metrics%wgtfac_e(nproma,nlevp1,nblks_e),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for wgtfac_e failed')
        ENDIF

       ! weighting factor for quadratic interpolation to surface
        ALLOCATE(p_metrics%wgtfacq_c(nproma,3,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for wgtfacq_c failed')
        ENDIF

       ! weighting factor for quadratic interpolation to surface
        ALLOCATE(p_metrics%wgtfacq_e(nproma,3,nblks_e),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for wgtfacq_e failed')
        ENDIF

       ! weighting factor for quadratic interpolation to model top
        ALLOCATE(p_metrics%wgtfacq1_c(nproma,3,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for wgtfacq1_c failed')
        ENDIF

       ! weighting factor for quadratic interpolation to model top
        ALLOCATE(p_metrics%wgtfacq1_e(nproma,3,nblks_e),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for wgtfacq1_e failed')
        ENDIF
       ! Inverse layer thickness of full levels
        ALLOCATE(p_metrics%inv_ddqz_z_full(nproma,nlev,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for inv_ddqz_z_full failed')
        ENDIF

       ! Inverse distance between full levels jk+1 and jk-1
        ALLOCATE(p_metrics%inv_ddqz_z_half2(nproma,nlev,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for inv_ddqz_z_half2 failed')
        ENDIF

      IF (i_cell_type== 3) THEN
       ! Vertical index of neighbor points needed for Taylor-expansion-based pressure gradient
        ALLOCATE(p_metrics%vertidx_gradp(2,nproma,nlev,nblks_e),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for vertidx_gradp failed')
        ENDIF

       ! Height differences between local edge point and neighbor cell points used for
       ! pressure gradient computation
        ALLOCATE(p_metrics%zdiff_gradp(2,nproma,nlev,nblks_e),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for zdiff_gradp failed')
        ENDIF

       ! Extrapolation factor for Exner pressure
        ALLOCATE(p_metrics%exner_exfac(nproma,nlev,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for exner_exfac failed')
        ENDIF

        ! Reference atmosphere fields
        ALLOCATE(p_metrics%theta_ref_mc(nproma,nlev,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for theta_ref_mc failed')
        ENDIF

        ALLOCATE(p_metrics%theta_ref_ic(nproma,nlevp1,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for theta_ref_ic failed')
        ENDIF

        ALLOCATE(p_metrics%exner_ref_mc(nproma,nlev,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for exner_ref_mc failed')
        ENDIF

        ALLOCATE(p_metrics%d_exner_dz_ref_ic(nproma,nlevp1,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for d_exner_dz_ref_ic failed')
        ENDIF

        ALLOCATE(p_metrics%d_exner_dz_ref_mc(nproma,nlev,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for d_exner_dz_ref_mc failed')
        ENDIF

        ALLOCATE(p_metrics%d2_exner_dz2_ref_mc(nproma,nlev,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for d2_exner_dz2_ref_mc failed')
        ENDIF

        ALLOCATE(p_metrics%rho_refcorr_ic(nproma,nlevp1,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for rho_refcorr_ic failed')
        ENDIF

        ALLOCATE(p_metrics%mask_prog_halo_c(nproma,nblks_c),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for mask_prog_halo_c failed')
        ENDIF

      ELSE IF (i_cell_type== 6) THEN
!      IF (i_cell_type== 6) THEN
        ! slope of the coordinate lines in Northern direction
        ALLOCATE(p_metrics%ddnorth_z(nproma,nlev,nblks_v),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for ddnorth_z failed')
        ENDIF
        ! slope of the coordinate lines in Eastern direction
        ALLOCATE(p_metrics%ddeast_z(nproma,nlev,nblks_v),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for ddeast_z failed')
        ENDIF
        ! functional determinant of the metrics [sqrt(gamma)]
        ALLOCATE(p_metrics%ddqz_z_full_v(nproma,nlev,nblks_v),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for ddqz_z_full_v failed')
        ENDIF

        ! functional determinant of the metrics [sqrt(gamma)]
        ALLOCATE(p_metrics%ddqz_z_half_e(nproma,nlevp1,nblks_e),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for ddqz_z_half_e failed')
        ENDIF
        ! functional determinant of the metrics [sqrt(gamma)]
        ALLOCATE(p_metrics%ddqz_z_half_r(nproma,nlevp1,nblks_e),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for ddqz_z_half_r failed')
        ENDIF
        ! functional determinant of the metrics [sqrt(gamma)]
        ALLOCATE(p_metrics%ddqz_z_half_v(nproma,nlevp1,nblks_v),STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                      'allocation for ddqz_z_half_v failed')
        ENDIF

      ENDIF

  END SUBROUTINE construct_nh_metrics

  !>
  !! Deallocation of components of diagnostic state.
  !!
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE destruct_nh_metrics (p_metrics)
!
    TYPE(t_nh_metrics), TARGET, INTENT(INOUT) :: p_metrics     ! current metrics state

    INTEGER     :: ist        ! status

      ! geometric height at the vertical interface of cells
      DEALLOCATE(p_metrics%z_ifc, STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for z_ifc failed')
      ENDIF

      ! geometric height at full levels
      DEALLOCATE(p_metrics%z_mc, STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for z_mc failed')
      ENDIF

      ! geometric height at full level edges
      DEALLOCATE(p_metrics%z_mc_e, STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for z_mc_e failed')
      ENDIF

      ! slope of the terrain in normal direction
      DEALLOCATE(p_metrics%ddxn_z_half,STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for ddxn_z_half failed')
      ENDIF
      ! slope of the terrain
      DEALLOCATE(p_metrics%ddxn_z_full,STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for ddxn_z_full failed')
      ENDIF
      ! slope of the terrain in tangential direction
      DEALLOCATE(p_metrics%ddxt_z_half,STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for ddxt_z_half failed')
      ENDIF

      ! functional determinant of the metrics [sqrt(gamma)]
      DEALLOCATE(p_metrics%ddqz_z_full, STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for ddqz_z_full failed')
      ENDIF
      ! functional determinant of the metrics [sqrt(gamma)]
      DEALLOCATE(p_metrics%ddqz_z_full_e,STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for ddqz_z_full_e failed')
      ENDIF
      ! functional determinant of the metrics [sqrt(gamma)]
      DEALLOCATE(p_metrics%ddqz_z_half,STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for ddqz_z_half failed')
      ENDIF

      ! 1/dz(k-1)-1/dz(k)
      DEALLOCATE(p_metrics%diff_1_o_dz,STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for diff_1_o_dz failed')
      ENDIF
      ! 1/(dz(k-1)*dz(k))
      DEALLOCATE(p_metrics%mult_1_o_dz,STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for mult_1_o_dz failed')
      ENDIF

      ! geopotential at cell center
      DEALLOCATE(p_metrics%geopot, STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for geopot failed')
      ENDIF

      ! geopotential at cell center
      DEALLOCATE(p_metrics%geopot_ifc, STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for geopot_ifc failed')
      ENDIF

      ! geopotential at cell center above groundlevel
      DEALLOCATE(p_metrics%geopot_agl,     STAT = ist)
      DEALLOCATE(p_metrics%geopot_agl_ifc, STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for geopot_agl failed')
      ENDIF

      ! Rayleigh damping
      DEALLOCATE(p_metrics%rayleigh_w,STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for rayleigh_w failed')
      ENDIF

      ! Explicit weight in vertical wind solver
      DEALLOCATE(p_metrics%vwind_expl_wgt,STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for vwind_expl_wgt failed')
      ENDIF

      ! Implicit weight in vertical wind solver
      DEALLOCATE(p_metrics%vwind_impl_wgt,STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                    'deallocation for vwind_impl_wgt failed')
      ENDIF

! These fields are needed for triangles only once the initialization in
! mo_nh_testcases is properly rewritten for hexagons
!      IF (i_cell_type== 3) THEN
        ! weighting factor for interpolation from full to half levels
        DEALLOCATE(p_metrics%wgtfac_c,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for wgtfac_c failed')
        ENDIF

       ! weighting factor for interpolation from full to half levels
        DEALLOCATE(p_metrics%wgtfac_e,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for wgtfac_e failed')
        ENDIF

       ! weighting factor for quadratic interpolation to surface
        DEALLOCATE(p_metrics%wgtfacq_c,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for wgtfacq_c failed')
        ENDIF

       ! weighting factor for quadratic interpolation to surface
        DEALLOCATE(p_metrics%wgtfacq_e,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for wgtfacq_e failed')
        ENDIF

       ! weighting factor for quadratic interpolation to model top
        DEALLOCATE(p_metrics%wgtfacq1_c,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for wgtfacq1_c failed')
        ENDIF

       ! weighting factor for quadratic interpolation to model top
        DEALLOCATE(p_metrics%wgtfacq1_e,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for wgtfacq1_e failed')
        ENDIF

       ! Inverse layer thickness of full levels
        DEALLOCATE(p_metrics%inv_ddqz_z_full,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for inv_ddqz_z_full failed')
        ENDIF

       ! Inverse distance between full levels jk+1 and jk-1
        DEALLOCATE(p_metrics%inv_ddqz_z_half2,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for inv_ddqz_z_half2 failed')
        ENDIF

      IF (i_cell_type== 3) THEN
       ! Vertical index of neighbor points needed for Taylor-expansion-based pressure gradient
        DEALLOCATE(p_metrics%vertidx_gradp,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for vertidx_gradp failed')
        ENDIF

       ! Height differences between local edge point and neighbor cell points used for
       ! pressure gradient computation
        DEALLOCATE(p_metrics%zdiff_gradp,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for zdiff_gradp failed')
        ENDIF

       ! Extrapolation factor for Exner pressure
        DEALLOCATE(p_metrics%exner_exfac,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for exner_exfac failed')
        ENDIF

        ! Reference atmosphere fields
        DEALLOCATE(p_metrics%theta_ref_mc,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for theta_ref_mc failed')
        ENDIF

        DEALLOCATE(p_metrics%theta_ref_ic,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for theta_ref_ic failed')
        ENDIF

        DEALLOCATE(p_metrics%exner_ref_mc,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for exner_ref_mc failed')
        ENDIF

        DEALLOCATE(p_metrics%d_exner_dz_ref_ic,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for d_exner_dz_ref_ic failed')
        ENDIF

        DEALLOCATE(p_metrics%d_exner_dz_ref_mc,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for d_exner_dz_ref_mc failed')
        ENDIF

        DEALLOCATE(p_metrics%d2_exner_dz2_ref_mc,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for d2_exner_dz2_ref_mc failed')
        ENDIF

        DEALLOCATE(p_metrics%rho_refcorr_ic,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for rho_refcorr_ic failed')
        ENDIF

        DEALLOCATE(p_metrics%mask_prog_halo_c,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for mask_prog_halo_c failed')
        ENDIF

      ELSE IF (i_cell_type== 6) THEN
!      IF (i_cell_type== 6) THEN
        ! slope of the coordinate lines in Northern direction
        DEALLOCATE(p_metrics%ddnorth_z,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for ddnorth_z failed')
        ENDIF
        ! slope of the coordinate lines in Eastern direction
        DEALLOCATE(p_metrics%ddeast_z,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for ddeast_z failed')
        ENDIF
        ! functional determinant of the metrics [sqrt(gamma)]
        DEALLOCATE(p_metrics%ddqz_z_full_v,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for ddqz_z_full_v failed')
        ENDIF

        ! functional determinant of the metrics [sqrt(gamma)]
        DEALLOCATE(p_metrics%ddqz_z_half_e,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for ddqz_z_half_e failed')
        ENDIF
        ! functional determinant of the metrics [sqrt(gamma)]
        DEALLOCATE(p_metrics%ddqz_z_half_r,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for ddqz_z_half_r failed')
        ENDIF
        ! functional determinant of the metrics [sqrt(gamma)]
        DEALLOCATE(p_metrics%ddqz_z_half_v,STAT=ist)
        IF (ist/=SUCCESS)THEN
          CALL finish('mo_nonhydro_state:destruct_nh_metrics', &
                      'deallocation for ddqz_z_half_v failed')
        ENDIF

      ENDIF


  END SUBROUTINE destruct_nh_metrics


END MODULE mo_nonhydro_state


