!>
!! Defines, constructs and destructs the state vector of the nonhydrostatic.
!!
!! Defines, constructs and destructs the state vector of the nonhydrostatic
!! model variables. They are subdivided in several classes: prognostics
!! and diagnostics.
!!
!! @author Almut Gassmann (MPI-M)
!! @author Daniel Reiner©t (DWD-M)
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2009-03-06)
!! Modification by Daniel Reinert, DWD (2011-05-02)
!! - Memory allocation method changed from explicit allocation to Luis'
!!   infrastructure
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
!
  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_model_domain_import, ONLY: n_dom, l_limited_area
  USE mo_nonhydrostatic_nml,  ONLY: l_nest_rcf
  USE mo_dynamics_nml,        ONLY: nsav1, nsav2, itime_scheme
!  USE mo_advection_nml,       ONLY: ctracer_list
  USE mo_run_nml,             ONLY: nproma, i_cell_type, iforcing,             &
    &                               inwp, ltransport, ntracer, ntracer_static, &
    &                               inextra_2d, inextra_3d,&
    &                               iqv, iqc, iqi, iqr, iqs, io3
  USE mo_radiation_nml,       ONLY: irad_o3
  USE mo_io_nml,              ONLY: lwrite_extra
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: default_var_list_settings, &
    &                               add_var,  add_ref,          &
    &                               new_var_list,              &
    &                               delete_var_list
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: t_nh_prog             ! state vector of prognostic variables
  PUBLIC :: t_nh_diag             ! state vector of diagnostic variables
  PUBLIC :: t_nh_metrics          ! state vector of metrics variables
  PUBLIC :: t_nh_state            ! state vector of nonhydrostatic variables

  PUBLIC :: construct_nh_state    ! Constructor for the nonhydrostatic state
  PUBLIC :: destruct_nh_state     ! Destructor for the nonhydrostatic state

  PUBLIC :: t_buffer_memory, bufr
  PUBLIC :: t_ptr_nh

  !>
  !! Derived data type for building pointer arrays
  !!
  TYPE t_ptr_nh
    REAL(wp),POINTER :: p_3d(:,:,:)  ! pointer to 3D (spatial) array
    REAL(wp),POINTER :: p_2d(:,:)    ! pointer to 2D (spatial) array
  END TYPE t_ptr_nh


  ! prognostic variables state vector
  TYPE t_nh_prog

    REAL(wp), POINTER :: &
      w(:,:,:),          & !> orthogonal vertical wind (nproma,nlevp1,nblks_c)     [m/s]
      vn(:,:,:),         & !! orthogonal normal wind (nproma,nlev,nblks_e)         [m/s]
      rho(:,:,:),        & !! density (nproma,nlev,nblks_c)                     [kg/m^3]
      exner(:,:,:),      & !! Exner pressure (nproma,nlev,nblks_c)                   [-]
      theta_v(:,:,:),    & !! virtual potential temperature (nproma,nlev,nblks_c)    [K]
      rhotheta_v(:,:,:), & !! rho*theta_v (nproma,nlev,nblks_c)               [K*kg/m^3]
      tracer(:,:,:,:),   & !! tracer concentration (nproma,nlev,nblks_c,ntracer) [kg/kg]
      tke   (:,:,:)        !! SQRT(2 * turbulent kinetik energy)                 [ m/s ]
                           !! (defined on half levels) with 2 time levels
    TYPE(t_ptr_nh),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer

  END TYPE t_nh_prog


  ! diagnostic variables state vector
  TYPE t_nh_diag

    REAL(wp), POINTER ::    &
    ! a) variables needed for triangles and hexagons
    &  u(:,:,:),            & ! zonal wind (nproma,nlev,nblks_c)               [m/s]
    &  v(:,:,:),            & ! meridional wind (nproma,nlev,nblks_c)          [m/s]
    &  vt(:,:,:),           & ! tangential wind (nproma,nlev,nblks_e)          [m/s]
    &  e_kin(:,:,:),        & ! spec. kinetic energy (nproma,nlev,nblks_c) [m^2/s^2]
    &  omega_z(:,:,:),      & ! vertical vorticity at dual grid
                              ! (nproma,nlev,nblks_v or nblks_e)               [1/s]
    &  ddt_vn(:,:,:),       & ! normal wind tendency from forcing              [m/s^2]
    &  ddt_vn_phy(:,:,:),   & ! normal wind tendency from forcing
                              ! (nproma,nlev,nblks_e)                          [m/s^2]
    &  ddt_w(:,:,:),        & ! vert. wind tendency from forcing
                              ! (nproma,nlevp1,nblks_c)                        [m/s^2]
    &  ddt_exner(:,:,:),    & ! exner pressure tendency from forcing (nproma,nlev,nblks_c)  [1/s]
    &  ddt_exner_phy(:,:,:),& ! exner pressure tendency from physical forcing 
                              ! (nproma,nlev,nblks_c)                     [1/s]
    &  ddt_tracer_phy(:,:,:,:), &! physics tendency of tracers
                              ! (nproma,nlev,nblks_c,ntracer)             [kg/kg/s]
    &  ddt_tracer_adv(:,:,:,:), &! advective tendency of tracers          [kg/kg/s]
    &  tracer_vi(:,:,:),    & ! vertically integrated tracers(for q1,q2,q3) [kg/m**2]
    &  tracer_vi_avg(:,:,:),& ! average since last output of tracer_vi [kg/m**2]
    &  exner_old(:,:,:),    & ! exner pres from previous step (nproma,nlev,nblks_c)
                            ! *** needs to be saved for restart ***
    &  w_con(:,:,:),        & ! contravariant vert wind (nproma,nlevp1,nblks_c)[1/s]
    &  temp(:,:,:),         & ! temperature (nproma,nlev,nblks_c)                 [K]
    &  tempv(:,:,:),        & ! virtual temperature (nproma,nlev,nblks_c)         [K]
    &  temp_ifc(:,:,:),     & ! temperature at half levels (nproma,nlevp1,nblks_c)[K]
    &  pres(:,:,:),         & ! pressure (nproma,nlev,nblks_c)                  [Pa]
    &  pres_ifc(:,:,:),     & ! pressure at interfaces (nproma,nlevp1,nblks_c)  [Pa]
    &  pres_sfc(:,:),       & ! diagnosed surface pressure (nproma,nblks_c)     [Pa]
    &  pres_sfc_s6avg(:,:), & ! 6 hourly sample  surface pressure average       [Pa]
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
    &  theta_v_ave(:,:,:),  & ! time average from horiz. adv. of theta_v         [K]
    &  horpgrad(:,:,:),     & ! covariant horizontal pressure gradient term   [m/s^2]
    &  vn_cov(:,:,:),       & ! covariant normal wind (nproma,nlev,nblks_e)    [m/s]
    &  w_cov(:,:,:),        & ! covariant vert wind (nproma,nlevp1,nblks_c)  [m^2/s]
    &  rho_e(:,:,:),        & ! density at edges nnow                       [kg/m^3]
    &  rho_star_e(:,:,:),   & ! density at edges estim. step                [kg/m^3]
    &  omega_z_con(:,:,:),  & ! vertical contrav. vorticity (nproma,nlev,nblks_v)
                              ! at vertices                                    [1/s]
    &  omega_t_con(:,:,:),  & ! tangential horiz. contravariant vorticity
                              ! (nproma,nlev,nblks_e)[1/s]
    &  omega_x(:,:,:),      & ! zonal vorticity (nproma,nlev,nblks_c)          [1/s]
    &  omega_y(:,:,:),      & ! meridional vorticity (nproma,nlev,nblks_c)     [1/s]
    &  ddt_vn_vort(:,:,:),  & ! normal wind tendency from vorticity flux term
                              ! (nproma,nlev,nblks_e,1:3)                    [m/s^2]
    &  ddt_w_vort(:,:,:),   & ! vert. wind tendency from vorticity flux term
                              ! (nproma,nlevp1,nblks_c,1:3)                  [m/s^2]
    &  ddt_w_phy(:,:,:)       ! vert. wind tendency from phyiscs
                              ! (nproma,nlevp1,nblks_c,1:3)                  [m/s^2]

    REAL(wp), POINTER ::    & !
     &  extra_2d(:,:,:)  ,  & !> extra debug output in 2d and
     &  extra_3d(:,:,:,:)     !!                       3d

    REAL(wp), POINTER :: &
      vn_con(:,:,:),     &! contravariant normal wind (nproma,nlev,nblks_e)[m/s]
      omega_t(:,:,:)      ! tangent. horiz. vorticity (nproma,nlev,nblks_e)[1/s]

    TYPE(t_ptr_nh),ALLOCATABLE :: ddt_grf_trc_ptr(:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_nh),ALLOCATABLE :: hfl_trc_ptr    (:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_nh),ALLOCATABLE :: vfl_trc_ptr    (:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_nh),ALLOCATABLE :: ddt_trc_adv_ptr(:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_nh),ALLOCATABLE :: ddt_trc_phy_ptr(:)  !< pointer array: one pointer for each tracer

    TYPE(t_ptr_nh),ALLOCATABLE :: ddt_vn_adv_ptr(:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_nh),ALLOCATABLE :: ddt_w_adv_ptr (:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_nh),ALLOCATABLE :: q_int_ptr     (:)
    TYPE(t_ptr_nh),ALLOCATABLE :: q_ubc_ptr     (:)
  END TYPE t_nh_diag


  TYPE t_nh_metrics

   ! a) Variables needed for triangles and hexagons
   !
   ! geometric height at the vertical interface of cells (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: z_ifc(:,:,:)
   ! geometric height at full levels (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: z_mc(:,:,:)
   ! geometric height at full level edges (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: z_mc_e(:,:,:)

   ! slope of the terrain in normal direction (nproma,nlevp1,nblks_e)
   REAL(wp), POINTER :: ddxn_z_half(:,:,:)
   ! slope of the terrain in normal direction (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddxn_z_full(:,:,:)
   ! slope of the terrain in tangential direction (nproma,nlevp1,nblks_e)
   REAL(wp), POINTER :: ddxt_z_half(:,:,:)

   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: ddqz_z_full(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddqz_z_full_e(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddqz_z_full_r(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: ddqz_z_half(:,:,:)

   ! 1/dz(k-1)-1/dz(k) (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: diff_1_o_dz(:,:,:)
   ! 1/(dz(k-1)*dz(k)) (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: mult_1_o_dz(:,:,:)

   ! geopotential at cell center (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: geopot(:,:,:)
   ! geopotential at interfaces at cell center (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: geopot_ifc(:,:,:)
   ! geopotential above groundlevel at cell center (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: geopot_agl(:,:,:)
   ! geopotential above groundlevel at interfaces and cell center (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: geopot_agl_ifc(:,:,:)
   ! geopotential at cell center (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: dgeopot_mc(:,:,:)


   ! Rayleigh damping on the vertical velocity
   REAL(wp), POINTER :: rayleigh_w(:,:,:)

   ! b) Variables needed for the triangular grid only
   !
   ! weighting factor for interpolation from full to half levels (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: wgtfac_c(:,:,:)
   ! weighting factor for interpolation from full to half levels (nproma,nlevp1,nblks_e)
   REAL(wp), POINTER :: wgtfac_e(:,:,:)
   ! weighting factor for quadratic interpolation to surface (nproma,3,nblks_c)
   REAL(wp), POINTER :: wgtfacq_c(:,:,:)
   ! weighting factor for quadratic interpolation to surface (nproma,3,nblks_e)
   REAL(wp), POINTER :: wgtfacq_e(:,:,:)
   ! weighting factor for quadratic interpolation to model top (nproma,3,nblks_c)
   REAL(wp), POINTER :: wgtfacq1_c(:,:,:)
   ! weighting factor for quadratic interpolation to model top (nproma,3,nblks_e)
   REAL(wp), POINTER :: wgtfacq1_e(:,:,:)
   ! Inverse layer thickness of full levels (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: inv_ddqz_z_full(:,:,:)
   ! Inverse distance between full levels jk+1 and jk-1 (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: inv_ddqz_z_half2(:,:,:)
   ! Vertical index of neighbor points needed for Taylor-expansion-based pressure gradient
   INTEGER,  POINTER :: vertidx_gradp(:,:,:,:) ! (2,nproma,nlev,nblks_e)
   ! Height differences between local edge point and neighbor cell points used for
   ! pressure gradient computation (2,nproma,nlev,nblks_e)
   REAL(wp), POINTER :: zdiff_gradp(:,:,:,:)
   ! extrapolation factor for Exner pressure (slope-dependent for stability optimization) 
   ! (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: exner_exfac(:,:,:)
   ! explicit weight in vertical wind solver (nproma,nblks_c)
   REAL(wp), POINTER :: vwind_expl_wgt(:,:)
   ! implicit weight in vertical wind solver (nproma,nblks_c)
   REAL(wp), POINTER :: vwind_impl_wgt(:,:)
   ! Fields for reference atmosphere
   REAL(wp), POINTER :: theta_ref_mc(:,:,:)
   REAL(wp), POINTER :: theta_ref_ic(:,:,:)
   REAL(wp), POINTER :: exner_ref_mc(:,:,:)
   REAL(wp), POINTER :: d_exner_dz_ref_ic(:,:,:)
   REAL(wp), POINTER :: d_exner_dz_ref_mc(:,:,:)
   REAL(wp), POINTER :: d2_exner_dz2_ref_mc(:,:,:)
   REAL(wp), POINTER :: rho_refcorr_ic(:,:,:)
   ! Fields for truly horizontal temperature diffusion
   INTEGER           :: zd_listdim
   INTEGER,  POINTER :: zd_indlist(:,:)
   INTEGER,  POINTER :: zd_blklist(:,:)
   INTEGER,  POINTER :: zd_edgeidx(:,:)
   INTEGER,  POINTER :: zd_edgeblk(:,:)
   INTEGER,  POINTER :: zd_vertidx(:,:)
   REAL(wp), POINTER :: zd_intcoef(:,:)
   REAL(wp), POINTER :: zd_geofac(:,:)
   REAL(wp), POINTER :: zd_e2cell(:,:)
   REAL(wp), POINTER :: zd_diffcoef(:)
   ! Fields for igradp_method = 3
   INTEGER           :: pg_listdim
   INTEGER,  POINTER :: pg_edgeidx(:)
   INTEGER,  POINTER :: pg_edgeblk(:)
   INTEGER,  POINTER :: pg_vertidx(:)
   REAL(wp), POINTER :: pg_exdist (:)
   ! Index lists for grid points on which lateral boundary nudging is applied
   INTEGER           :: nudge_c_dim, nudge_e_dim
   INTEGER,  POINTER :: nudge_c_idx(:)
   INTEGER,  POINTER :: nudge_e_idx(:)
   INTEGER,  POINTER :: nudge_c_blk(:)
   INTEGER,  POINTER :: nudge_e_blk(:)
   ! Index lists and mask fields needed to minimize the number of halo communications
   ! a) index lists for halo points belonging to the nest boundary region
   INTEGER           :: bdy_halo_c_dim  
   INTEGER,  POINTER :: bdy_halo_c_idx(:)
   INTEGER,  POINTER :: bdy_halo_c_blk(:)
   ! b) a mask field that excludes boundary halo points
   LOGICAL,  POINTER :: mask_prog_halo_c(:,:)
   ! c) index lists for halo points belonging to a nest overlap region
   !    the additional dimension is n_childdom
   INTEGER,  POINTER :: ovlp_halo_c_dim(:)
   INTEGER,  POINTER :: ovlp_halo_c_idx(:,:)
   INTEGER,  POINTER :: ovlp_halo_c_blk(:,:)


   ! c) Variables needed for the hexagonal grid only
   !
   ! slope of the coordinate lines towards the North (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddnorth_z(:,:,:)
   ! slope of the coordinate lines towards the East (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddeast_z(:,:,:)

   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_v)
   REAL(wp), POINTER :: ddqz_z_full_v(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlevp1,nblks_e)
   REAL(wp), POINTER :: ddqz_z_half_e(:,:,:)

  END TYPE t_nh_metrics

  TYPE :: t_buffer_memory

    REAL(wp), POINTER :: send_c1(:,:)
    REAL(wp), POINTER :: recv_c1(:,:)
    REAL(wp), POINTER :: send_c3(:,:)
    REAL(wp), POINTER :: recv_c3(:,:)
    REAL(wp), POINTER :: send_e1(:,:)
    REAL(wp), POINTER :: recv_e1(:,:)
    REAL(wp), POINTER :: send_e2(:,:)
    REAL(wp), POINTER :: recv_e2(:,:)
    REAL(wp), POINTER :: send_e3(:,:)
    REAL(wp), POINTER :: recv_e3(:,:)
    REAL(wp), POINTER :: send_v2(:,:)
    REAL(wp), POINTER :: recv_v2(:,:)

  END TYPE t_buffer_memory

  TYPE (t_buffer_memory), POINTER :: bufr(:)


!-------------------------------------------------------------------------
!                      STATE VECTORS AND LISTS
!-------------------------------------------------------------------------
  TYPE t_nh_state

    !array of prognostic states at different timelevels
    TYPE(t_nh_prog),  ALLOCATABLE :: prog(:)       !< shape: (timelevels)
    TYPE(t_var_list), ALLOCATABLE :: prog_list(:)  !< shape: (timelevels)

    TYPE(t_nh_diag)    :: diag
    TYPE(t_var_list)   :: diag_list

    TYPE(t_nh_metrics) :: metrics
    TYPE(t_var_list)   :: metrics_list

  END TYPE t_nh_state



  CONTAINS

!-------------------------------------------------------------------------
!!            SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
!-------------------------------------------------------------------------
!
!
  !>
  !! Constructor for prognostic and diagnostic states.
  !!
  !! Top-level procedure for building the prognostic and diagnostic states.
  !! It calls constructors to single time level prognostic states, and 
  !! diagnostic states.
  !! Initialization of all components with zero.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE construct_nh_state(p_patch, p_nh_state, n_timelevels)
!
    TYPE(t_patch),     INTENT(IN)   ::  & ! patch
      &  p_patch(n_dom)

    TYPE(t_nh_state),  INTENT(INOUT)::  & ! nh state at different grid levels
      &  p_nh_state(n_dom)

    INTEGER, OPTIONAL, INTENT(IN)   ::  & ! number of timelevels
      &  n_timelevels    

    INTEGER  :: ntl, &    ! local number of timelevels
                ist, &    ! status
                jg,  &    ! grid level counter
                jt        ! time level counter

    LOGICAL  :: l_alloc_tracer

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname, varname_prefix

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nonhydro_state:construct_nh_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'Construction of NH state started')

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

      !
      ! Allocate pointer array p_nh_state(jg)%prog, as well as the 
      ! corresponding list array for each grid level.
      !
      ! create state arrays
      ALLOCATE(p_nh_state(jg)%prog(1:ntl), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine),                               &
          &          'allocation of prognostic state array failed')
      ENDIF

      ! create state list
      ALLOCATE(p_nh_state(jg)%prog_list(1:ntl), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine),                                    &
          &          'allocation of prognostic state list array failed')
      ENDIF

      !
      ! Build prog list for every timelevel
      ! includes memory allocation
      ! 
      DO jt = 1, ntl

        ! Tracer fields do not need extra time levels because feedback is not incremental
        ! and the nest-call frequency is always synchronized with the advection time step
        IF (jt <= n_timelevels) THEN
          l_alloc_tracer = .TRUE.
        ELSE
          l_alloc_tracer = .FALSE.
        ENDIF

        WRITE(listname,'(a,i2.2,a,i2.2)') 'nh_state_prog_of_domain_',jg, &
          &                               '_and_timelev_',jt
        WRITE(varname_prefix,'(a,i2.2,a)') 'nh_prog_TL',jt,'_'
        CALL construct_nh_state_prog_list(p_patch(jg), p_nh_state(jg)%prog(jt), &
          &  p_nh_state(jg)%prog_list(jt), listname, TRIM(varname_prefix)     , &
          &  l_alloc_tracer)

      ENDDO

      !
      ! Build diag state list
      ! includes memory allocation
      !
      WRITE(listname,'(a,i2.2)') 'nh_state_diag_of_domain_',jg
      CALL construct_nh_state_diag_list(p_patch(jg), p_nh_state(jg)%diag, &
        &  p_nh_state(jg)%diag_list, listname)

      !
      ! Build metrics state list
      ! includes memory allocation
      !
      WRITE(listname,'(a,i2.2)') 'nh_state_metrics_of_domain_',jg
      CALL construct_nh_metrics_list(p_patch(jg), p_nh_state(jg)%metrics, &
        &  p_nh_state(jg)%metrics_list, listname )

    ENDDO

    CALL message (TRIM(routine), 'NH state construction completed')

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
    TYPE(t_nh_state), INTENT(INOUT) :: & ! nh state at different grid levels
      &  p_nh_state(n_dom)
                                             
    INTEGER  :: ntl, &    ! local number of timelevels
                ist, &    ! status
                jg,  &    ! grid level counter
                jt        ! time level counter

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nonhydro_state:destruct_nh_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of NH state started')


    DO jg = 1, n_dom

      ntl = SIZE(p_nh_state(jg)%prog(:))
      IF(ntl==0)THEN
        CALL finish(TRIM(routine), 'prognostic array has no timelevels')
      ENDIF

      ! delete diagnostic state list elements
      CALL delete_var_list( p_nh_state(jg)%diag_list )

      ! delete metrics state list elements
      CALL delete_var_list( p_nh_state(jg)%metrics_list )

      ! delete prognostic state list elements
      DO jt = 1, ntl
        CALL delete_var_list( p_nh_state(jg)%prog_list(jt) )
      ENDDO

      ! destruct state lists and arrays
      DEALLOCATE(p_nh_state(jg)%prog_list, STAT=ist )
      IF(ist/=SUCCESS) CALL finish (TRIM(routine),&
        & 'deallocation of prognostic state list array failed')

      DEALLOCATE(p_nh_state(jg)%prog )
      IF(ist/=SUCCESS) CALL finish (TRIM(routine),&
        & 'deallocation of prognostic state array failed')

    ENDDO

    CALL message (TRIM(routine), 'NH state destruction completed')

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
  SUBROUTINE construct_nh_state_prog_list ( p_patch, p_prog, p_prog_list,  &
    &                                       listname,vname_prefix, l_alloc_tracer)
!
    TYPE(t_patch), TARGET, INTENT(IN) :: & !< current patch
      &  p_patch

    TYPE(t_nh_prog),  INTENT(INOUT)   :: & !< current prognostic state
      &  p_prog 

    TYPE(t_var_list), INTENT(INOUT)   :: p_prog_list !< current prognostic state list

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  listname, vname_prefix

    LOGICAL, INTENT(IN) :: l_alloc_tracer  !< allocate tracer fields if true


    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e       !< number of edge blocks to allocate

    INTEGER :: nlev, nlevp1, ktracer

    INTEGER :: shape3d_c(3), shape3d_e(3), shape3d_chalf(3), &
      &        shape4d_c(4)

    INTEGER :: ientr         !< "entropy" of horizontal slice
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape3d_e     = (/nproma, nlev,   nblks_e  /)
    shape3d_c     = (/nproma, nlev,   nblks_c  /)
    shape3d_chalf = (/nproma, nlevp1, nblks_c  /)
    shape4d_c     = (/nproma, nlev,   nblks_c, ntracer+ntracer_static/)

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_prog_list, TRIM(listname) )
    CALL default_var_list_settings( p_prog_list,               &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )


    !------------------------------
    ! Meteorological quantities
    !------------------------------

    ! vn           p_prog%vn(nproma,nlev,nblks_e)
    cf_desc    = t_cf_var('normal_velocity', 'm s-1', 'velocity normal to edge')
    grib2_desc = t_grib2_var(0, 2, 197, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'vn', p_prog%vn,             &
      &           GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
      &           ldims=shape3d_e )


    ! rho          p_prog%rho(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('density', 'kg m-3', 'density')
    grib2_desc = t_grib2_var(0, 3, 10, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'rho', p_prog%rho,             &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
      &           ldims=shape3d_c )


    ! exner        p_prog%exner(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('exner_pressure', '-', 'exner pressure')
    grib2_desc = t_grib2_var(0, 3, 195, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'exner', p_prog%exner,        &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,     &
      &           ldims=shape3d_c )


    ! theta_v      p_prog%theta_v(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('virtual_potential_temperature', 'K', &
      &                    'virtual potential temperature')
    grib2_desc = t_grib2_var(0, 0, 1, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'theta_v', p_prog%theta_v,   &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
      &           ldims=shape3d_c )


    ! rhotheta_v   p_prog%rhotheta_v(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('rho_virt_pot_temp', 'K', 'rho virt pot temp')
    grib2_desc = t_grib2_var(0, 19, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'rhotheta_v', p_prog%rhotheta_v, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,     &
      &           ldims=shape3d_c )

    ! Tracer array for (model) internal use

    ! tracer         p_prog%tracer(nproma,nlev,nblks_c,ntracer+ntracer_static)
    IF ( (ltransport .OR. iforcing == inwp) .AND. l_alloc_tracer  ) THEN
      cf_desc    = t_cf_var('tracer', 'kg kg-1', 'tracer')
      grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_prog_list, 'tracer', p_prog%tracer,                       &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
        &           ldims=shape4d_c ,                                           &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)


      ! Reference to individual tracer, for I/O

    ktracer=ntracer+ntracer_static
    ALLOCATE( p_prog%tracer_ptr(ktracer) )

           !QV
        CALL add_ref( p_prog_list, 'tracer',                                   &
                    & TRIM(vname_prefix)//'qv', p_prog%tracer_ptr(iqv)%p_3d,   &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                    &
                    & t_cf_var(TRIM(vname_prefix)//'qv',                       &
                    &  'kg kg-1','specific_humidity'),                         &
                    & t_grib2_var(0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d_c)
           !QC
        CALL add_ref( p_prog_list, 'tracer',&
                    & TRIM(vname_prefix)//'qc', p_prog%tracer_ptr(iqc)%p_3d,       &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                        &
                    & t_cf_var(TRIM(vname_prefix)//'qc',                           &
                    &  'kg kg-1', 'specific_cloud_water_content'),                 &
                    & t_grib2_var(192, 201, 31, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)
           !QI
        CALL add_ref( p_prog_list, 'tracer',                                       &
                    & TRIM(vname_prefix)//'qi', p_prog%tracer_ptr(iqi)%p_3d,       &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                        &
                    & t_cf_var(TRIM(vname_prefix)//'qi',                           & 
                    &  'kg kg-1','specific_cloud_ice_content'),                    &
                    & t_grib2_var(192, 201, 33, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)
           !QR
        CALL add_ref( p_prog_list, 'tracer',                                    &
                    & TRIM(vname_prefix)//'qr', p_prog%tracer_ptr(iqr)%p_3d,    &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                     &
                    & t_cf_var(TRIM(vname_prefix)//'qr',                        &       
                    &  'kg kg-1','rain_mixing_ratio'),                          &
                    & t_grib2_var(0, 1, 24, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d_c)
           !QS
        CALL add_ref( p_prog_list, 'tracer',                                   &
                    & TRIM(vname_prefix)//'qs', p_prog%tracer_ptr(iqs)%p_3d,   &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                    &
                    & t_cf_var(TRIM(vname_prefix)//'qs',                       &
                    &  'kg kg-1','snow_mixing_ratio'),                         &
                    & t_grib2_var(0, 1, 25, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)

        IF(irad_o3 == 3) THEN
           !O3
          CALL add_ref( p_prog_list, 'tracer',                         &
            & TRIM(vname_prefix)//'O3', p_prog%tracer_ptr(io3)%p_3d,   &
            & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                    &
            & t_cf_var(TRIM(vname_prefix)//'O3',                       &
            &  'kg kg-1','ozone_mass_mixing_ratio'),                   &
            & t_grib2_var(0, 14, 1, ientr, GRID_REFERENCE, GRID_CELL), &
            & ldims=shape3d_c)
        ENDIF

    ENDIF ! ltransport or iforcing
    !
    ! variables defined at half levels 
    !

    ! w            p_prog%w(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('upward air velocity', 'm s-1', 'upward air velocity')
    grib2_desc = t_grib2_var(0, 2, 9, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'w', p_prog%w,               &
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,     &
      &          ldims=shape3d_chalf )


    ! tke            p_prog%tke(nproma,nlevp1,nblks_c)
    IF ( iforcing == inwp ) THEN
      cf_desc    = t_cf_var('turbulent_kinetic_energy', 'm2 s-2', 'turbulent kinetic energy')
      grib2_desc = t_grib2_var(0, 19, 11, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_prog_list, TRIM(vname_prefix)//'tke', p_prog%tke, &
           &        GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, &
        &           cf_desc, grib2_desc, ldims=shape3d_chalf )
    ENDIF

  END SUBROUTINE construct_nh_state_prog_list




  !-------------------------------------------------------------------------
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
  !!
  SUBROUTINE construct_nh_state_diag_list ( p_patch, p_diag, p_diag_list,  &
    &                                       listname )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_diag),  INTENT(INOUT)   :: &  !< diagnostic state
      &  p_diag 

    TYPE(t_var_list), INTENT(INOUT)   :: p_diag_list  !< diagnostic state list

    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e, &    !< number of edge blocks to allocate
               nblks_v       !< number of vertex blocks to allocate

    INTEGER :: nlev, nlevp1

    INTEGER :: n_timlevs     !< number of time levels for advection 
                             !< tendency fields

    INTEGER :: shape2d_c(2), shape2d_e(2), shape3d_c(3),           &
      &        shape3d_e(3), shape3d_v(3), shape3d_chalf(3),       &
      &        shape3d_ehalf(3), shape4d_chalf(4), shape4d_e(4),   &
      &        shape4d_entl(4), shape4d_chalfntl(4), shape4d_c(4), &
      &        shape3d_ctra(3), shape2d_extra(3), shape3d_extra(4)
 
    INTEGER :: ientr         !< "entropy" of horizontal slice
    INTEGER :: jt

    CHARACTER(LEN=2) :: ctrc
    !--------------------------------------------------------------

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

    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c     = (/nproma,          nblks_c    /)
    shape2d_e     = (/nproma,          nblks_e    /)
    shape2d_extra = (/nproma, nblks_c, inextra_2d /)
    shape3d_c     = (/nproma, nlev   , nblks_c    /)
    shape3d_e     = (/nproma, nlev   , nblks_e    /)
    shape3d_v     = (/nproma, nlev   , nblks_v    /)
    shape3d_chalf = (/nproma, nlevp1 , nblks_c    /)
    shape3d_ehalf = (/nproma, nlevp1 , nblks_e    /)
    shape3d_ctra  = (/nproma, nblks_c, ntracer    /)
    shape3d_extra = (/nproma, nlev   , nblks_c, inextra_3d  /)
    shape4d_c     = (/nproma, nlev   , nblks_c, ntracer     /)
    shape4d_chalf = (/nproma, nlevp1 , nblks_c, ntracer     /)
    shape4d_e     = (/nproma, nlev   , nblks_e, ntracer     /)
    shape4d_entl  = (/nproma, nlev   , nblks_e, n_timlevs   /)
    shape4d_chalfntl = (/nproma, nlevp1, nblks_c, n_timlevs /)


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_diag_list, TRIM(listname) )
    CALL default_var_list_settings( p_diag_list,               &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    ! u           p_diag%u(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'u-component of wind')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'u', p_diag%u,                                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! v           p_diag%v(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'v-component of wind')
    grib2_desc = t_grib2_var(0, 2, 3, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'v', p_diag%v,                                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )

    ! vt           p_diag%vt(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('tangential_wind', 'm s-1', 'tangential-component of wind')
    grib2_desc = t_grib2_var(0, 2, 198, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'vt', p_diag%vt,                                 &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e )


    ! e_kin        p_diag%e_kin(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('specific_kinetic_energy', 'm2 s-2', 'specific kinetic energy')
    grib2_desc = t_grib2_var(0, 2, 196, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'e_kin', p_diag%e_kin,                           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! omega_z      p_diag%omega_z(nproma,nlev,nblks_v)
    !
    cf_desc    = t_cf_var('vertical_vorticity', 'm s-1', 'vertical voritcity')
    grib2_desc = t_grib2_var(0, 2, 197, ientr, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( p_diag_list, 'omega_z', p_diag%omega_z,                       &
                & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_v )


    ! ddt_vn       p_diag%ddt_vn(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2', 'normal wind tendency')
    grib2_desc = t_grib2_var(0, 2, 198, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'ddt_vn', p_diag%ddt_vn,                         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e )


    ! ddt_vn_phy   p_diag%ddt_vn_phy(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('normal_wind_physical_tendency', 'm s-2',             &
      &                   'normal wind physical tendency')
    grib2_desc = t_grib2_var(0, 2, 199, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'ddt_vn_phy', p_diag%ddt_vn_phy,                 &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e )


    ! ddt_w        p_diag%ddt_w(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('vertical_wind_tendency', 'm s-2', 'vertical wind tendency')
    grib2_desc = t_grib2_var(0, 2, 200, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_w', p_diag%ddt_w,                           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! ddt_exner    p_diag%ddt_exner(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('exner_pressure_tendency', 's-1', 'exner pressure tendency')
    grib2_desc = t_grib2_var(0, 3, 196, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_exner', p_diag%ddt_exner,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! ddt_exner_phy  p_diag%ddt_exner_phy(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('exner_pressure_physical_tendency', 's-1',            &
      &                   'exner pressure physical tendency')
    grib2_desc = t_grib2_var(0, 3, 197, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_exner_phy', p_diag%ddt_exner_phy,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! exner_old    p_diag%exner_old(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('old_exner_pressure', '-', 'old exner pressure')
    grib2_desc = t_grib2_var(0, 3, 196, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'exner_old', p_diag%exner_old,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! w_con        p_diag%w_con(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('contravariant_vertical_wind', 'm s-1',               &
      &                   'contravariant_vertical_wind')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'w_con', p_diag%w_con,                           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! pres_sfc     p_diag%pres_sfc(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_sfc', p_diag%pres_sfc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d_c )

    ! pres_sfc_s6avg     p_diag%pres_sfc_s6avg(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_sfc_s6avg', p_diag%pres_sfc_s6avg,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d_c )

    ! temp         p_diag%temp(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('temperature', 'K', 'temperature')
    grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'temp', p_diag%temp,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! tempv        p_diag%tempv(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('virtual_temperature', 'K', 'virtual temperature')
    grib2_desc = t_grib2_var(0, 0, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'tempv', p_diag%tempv,                           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! temp_ifc     p_diag%temp_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('temperature', 'K', 'temperature at half level')
    grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'temp_ifc', p_diag%temp_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! pres         p_diag%pres(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres', p_diag%pres,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! pres_ifc     p_diag%pres_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at half level')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_ifc', p_diag%pres_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! dpres_mc     p_diag%dpres_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('pressure_thickness', 'Pa', 'pressure thickness')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'dpres_mc', p_diag%dpres_mc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c )


    ! div          p_diag%div(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('divergence', 's-1', 'divergence')
    grib2_desc = t_grib2_var( 0, 2, 13, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'div', p_diag%div,                               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c )


    ! mass_fl_e    p_diag%mass_fl_e(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('horizontal_mass_flux_at_edges', 'kg m-1 s-1',        &
       &         'horizontal mass flux at edges')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'mass_fl_e', p_diag%mass_fl_e,                   &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_e )


    ! rho_ic       p_diag%rho_ic(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('density', 'kg m-3', 'density at half level')
    grib2_desc = t_grib2_var( 0, 3, 10, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'rho_ic', p_diag%rho_ic,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! w_concorr_c  p_diag%w_concorr_c(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('contravariant_vertical_correction', 'm s-1',         &
      &                   'contravariant vertical correction')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'w_concorr_c', p_diag%w_concorr_c,               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! e_kinh       p_diag%e_kinh(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('horizontal specific kinetic energy', 'm2 s-2',       &
      &                   'horizontal specific kinetic energy')
    grib2_desc = t_grib2_var( 0, 2, 196, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'e_kinh', p_diag%e_kinh,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c )


    IF (i_cell_type == 3) THEN
      ! vn_ie        p_diag%vn_ie(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_at_half_level', 'm s-1',               &
        &                   'normal wind at half level')
      grib2_desc = t_grib2_var( 0, 2, 197, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'vn_ie', p_diag%vn_ie,                         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf )


      ! theta_v_ic   p_diag%theta_v_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('virtual_potential_temperature_at_half_levels', 'K',&
        &                   'virtual_potential temperature at half levels')
      grib2_desc = t_grib2_var( 0, 0, 1, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'theta_v_ic', p_diag%theta_v_ic,               &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )


      ! ddt_vn_adv   p_diag%ddt_vn_adv(nproma,nlev,nblks_e,n_timlevs)
      !
      cf_desc    = t_cf_var('advective_normal_wind_tendency', 'm s-2',          &
        &                   'advective normal wind tendency')
      grib2_desc = t_grib2_var( 0, 2, 201, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_adv', p_diag%ddt_vn_adv,               &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape4d_entl ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

      ALLOCATE(p_diag%ddt_vn_adv_ptr(n_timlevs))
      DO jt =1,n_timlevs
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'ddt_vn_adv',                                   &
                    & 'ddt_adv_vn'//ctrc, p_diag%ddt_vn_adv_ptr(jt)%p_3d,          &
                    & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT,                        &
                    & t_cf_var('ddt_adv_vn'//ctrc, 'm s-2',''),                    &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape3d_e)
      ENDDO


      ! ddt_w_adv    p_diag%ddt_w_adv(nproma,nlevp1,nblks_c,n_timlevs)
      !
      cf_desc    = t_cf_var('advective_vertical_wind_tendency', 'm s-2',        &
        &                   'advective vertical wind tendency')
      grib2_desc = t_grib2_var( 0, 2, 202, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_w_adv', p_diag%ddt_w_adv,                 &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_chalfntl ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

      ALLOCATE(p_diag%ddt_w_adv_ptr(n_timlevs))
      DO jt =1,n_timlevs
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'ddt_w_adv',                                     &
                    & 'ddt_w_adv'//ctrc, p_diag%ddt_w_adv_ptr(jt)%p_3d,             &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                         &
                    & t_cf_var('ddt_adv_w'//ctrc, 'm s-2',''),                      &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_chalf)
      ENDDO

      ! grf_tend_vn  p_diag%grf_tend_vn(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2',                    &
        &                   'normal wind tendency (grid refinement)')
      grib2_desc = t_grib2_var( 0, 2, 203, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'grf_tend_vn', p_diag%grf_tend_vn,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_e )


      ! grf_tend_w  p_diag%grf_tend_w(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vertical_wind_tendency', 'm s-2',                  &
        &                   'vertical wind tendency (grid refinement)')
      grib2_desc = t_grib2_var( 0, 2, 204, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_w', p_diag%grf_tend_w,               &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )

      ! grf_tend_rho   p_diag%grf_tend_rho(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('density_tendency', 'kg m-3 s-1',                   &
        &                   'density tendency (grid refinement)')
      grib2_desc = t_grib2_var( 0, 3, 198, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_rho', p_diag%grf_tend_rho,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )

      ! grf_tend_thv   p_diag%grf_tend_thv(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('virtual_potential_temperature_tendency', 'K s-1',  &
        &                   'virtual potential temperature tendency (grid refinement)')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_thv', p_diag%grf_tend_thv,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )


      ! Storage fields for vertical nesting; the middle index (2) addresses 
      ! the field and its temporal tendency

      ! dvn_ie_int   p_diag%dvn_ie_int(nproma,nblks_e)
      !
      cf_desc    = t_cf_var('normal_velocity_parent_interface_level', 'm s-1',  &
        &                   'normal velocity at parent interface level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'dvn_ie_int', p_diag%dvn_ie_int,               &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_e )


      ! dvn_ie_ubc   p_diag%dvn_ie_ubc(nproma,nblks_e)
      !
      cf_desc    = t_cf_var('normal_velocity_child_upper_boundary', 'm s-1',    &
        &                   'normal velocity at child upper boundary')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'dvn_ie_ubc', p_diag%dvn_ie_ubc,               &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_e )


      ! drho_ic_int  p_diag%drho_ic_int(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('rho_at_parent_interface_level', 'kg m-3',          &
        &                   'rho at parent interface level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'drho_ic_int', p_diag%drho_ic_int,             &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c )


      ! drho_ic_ubc  p_diag%drho_ic_ubc(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('density_at_child_upper_boundary', 'kg m-3',        &
        &                   'density at child upper boundary')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'drho_ic_ubc', p_diag%drho_ic_ubc,             &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c )


      ! dtheta_v_ic_int    p_diag%dtheta_v_ic_int(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('theta_at_parent_interface_level', 'K',             &
        &                   'potential temperature at parent interface level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dtheta_v_ic_int', p_diag%dtheta_v_ic_int,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c )


      ! dtheta_v_ic_ubc    p_diag%dtheta_v_ic_ubc(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('theta_at_child_upper_boundary', 'K',               &
        &                   'potential temperature at child upper boundary')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dtheta_v_ic_ubc', p_diag%dtheta_v_ic_ubc,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c )


      ! dw_int       p_diag%dw_int(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('w_at_parent_interface_level', 'm s-1',             &
        &                   'vertical velocity at parent interface level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dw_int', p_diag%dw_int,                       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape2d_c )


      ! dw_ubc       p_diag%dw_ubc(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('w at child upper boundary', 'm s-1',               &
        &                   'vertical velocity at child upper boundary')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dw_ubc', p_diag%dw_ubc,                       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape2d_c )


      ! q_int        p_diag%q_int(nproma,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('q_at_parent_interface_level', 'kg kg-1',           &
        &                   'q at parent interface level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'q_int', p_diag%q_int,                         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape3d_ctra ,                                        &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

      ALLOCATE(p_diag%q_int_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'q_int',                                         &
                    & 'q_int'//ctrc, p_diag%q_int_ptr(jt)%p_2d,                     &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
                    & t_cf_var('q_int'//ctrc, 'kg kg-1',''),                        &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape2d_c)
      ENDDO


      ! q_ubc        p_diag%q_ubc(nproma,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('q_at_child_upper_boundary', 'kg kg-1',             &
        &                   'q at child upper boundary')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'q_ubc', p_diag%q_ubc,                         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape3d_ctra,                                         &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

      ALLOCATE(p_diag%q_ubc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'q_ubc',                                         &
                    & 'q_ubc'//ctrc, p_diag%q_ubc_ptr(jt)%p_2d,                     &
                    & GRID_UNSTRUCTURED_CELL,ZAXIS_SURFACE,                         &
                    & t_cf_var('q_ubc'//ctrc, 'kg kg-1',''),                        &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape2d_c)
      ENDDO


      ! thermal_exp_fastphy     p_diag%thermal_exp_fastphy(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('thermal_expansion_due_to_fast_physics', '',        &
        &                   'thermal expansion due to fast physics')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'thermal_exp_fastphy', p_diag%thermal_exp_fastphy, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,     &
                  & ldims=shape2d_c )


    ELSE IF (i_cell_type == 6) THEN

      ! theta_v_ic   p_diag%theta_v_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('virtual potential_temperature_at_half_level', 'K', &
        &                   'virtual potential temperature at half level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'theta_v_ic', p_diag%theta_v_ic,               &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )


      ! theta_v_impl   p_diag%theta_v_impl(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('(nnow+nnew)/2_from_impl._vert._adv._of_theta_v', 'K', &
        &                   '(nnow+nnew)/2 from impl. vert. adv. of theta_v')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'theta_v_impl', p_diag%theta_v_impl,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )

      ! theta_v_ave   p_diag%theta_v_ave(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('time average from horiz. adv. of theta_v', 'K', &
        &                   'time average from horiz. adv. of theta_v')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'theta_v_ave', p_diag%theta_v_ave,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,&
                  & ldims=shape3d_c )


      ! horpgrad     p_diag%horpgrad(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('covariant_horizontal_pressure_gradient', 'Pa m-1', &
        &                   'covariant horizontal pressure gradient')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'horpgrad', p_diag%horpgrad,                   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )


      ! vn_cov       p_diag%vn_cov(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('covariant_normal_wind', 'm s-1',                   &
        &                   'covariant normal wind')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'vn_cov', p_diag%vn_cov,                       &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )


      ! w_cov        p_diag%w_cov(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('covariant_vertical_wind', 'm s-1',                 &
        &                   'covariant vertical wind at half level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'w_cov', p_diag%w_cov,                         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  & 
                  & ldims=shape3d_chalf )

      ! rho_e        p_diag%rho_e(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('density_at_edges_nnow', 'm s-1',                   &
        &                   'density_at_edges_nnow')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'rho_e', p_diag%rho_e,                         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )

      ! rho_star_e   p_diag%rho_star_e(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('density_at_edges_estim_step', 'm s-1',              &
        &                   'density_at_edges_estim_step')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'rho_star_e', p_diag%rho_star_e,                &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,   &
                  & ldims=shape3d_e )

      ! omega_t_con  p_diag%omega_t_con(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('tangential_horiz._contravariant_vorticity', 's-1', &
        &                   'tangential horiz. contravariant vorticity')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'omega_t_con', p_diag%omega_t_con,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf )


      ! omega_x      p_diag%omega_x(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('zonal_vorticity', 's-1', 'zonal vorticity')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'omega_x', p_diag%omega_x,                     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )


      ! omega_y      p_diag%omega_y(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('meridional_vorticity', 's-1', 'meridional vorticity')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'omega_y', p_diag%omega_y,                     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )


      ! omega_z_con  p_diag%omega_z_con(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('vertical_vorticity', 's-1',                        &
        &                   'contravariant vertical vorticity')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_diag_list, 'omega_z_con', p_diag%omega_z_con,             &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_v )


      ! ddt_vn_vort  p_diag%ddt_vn_vort(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2',                    &
        &                   'normal wind tendency from vorticity flux term')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_vort', p_diag%ddt_vn_vort,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )


      ! ddt_w_vort   p_diag%ddt_w_vort(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vert._wind_tendency', 'm s-2',                     &
        &                   'vert. wind tendency from vorticity flux term')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_w_vort', p_diag%ddt_w_vort,               &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )

      ! ddt_w_phy   p_diag%ddt_w_phy(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vert_wind_tendency_phy', 'm s-2',                   &
        &                   'vertical wind tendency from physics')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_w_phy', p_diag%ddt_w_phy,                  &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,   &
                  & ldims=shape3d_chalf )

    ENDIF


    !
    ! tracers
    !
    IF ( ltransport ) THEN
      ! grf_tend_tracer   p_diag%grf_tend_tracer(nproma,nlev,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('tracer_tendency', 'kg kg-1 s-1',                   &
        &                   'tracer_tendency for grid refinement')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_tracer', p_diag%grf_tend_tracer,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_c ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

      ALLOCATE(p_diag%ddt_grf_trc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'grf_tend_tracer',                              &
                    & 'ddt_grf_q'//ctrc, p_diag%ddt_grf_trc_ptr(jt)%p_3d,             &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                        &
                    & t_cf_var('ddt_grf_q'//ctrc, 'kg kg-1 s**-1',''),             &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape3d_c)
      ENDDO


      ! hfl_tracer   p_diag%hfl_tracer(nproma,nlev,nblks_e,ntracer)
      !
      cf_desc    = t_cf_var('horizontal tracer flux', 'kg m-1 s-1',               &
        &                   'horizontal tracer flux')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'hfl_tracer', p_diag%hfl_tracer,                 &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape4d_e ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

      ALLOCATE(p_diag%hfl_trc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'hfl_tracer',                                    &
                    & 'hfl_q'//ctrc, p_diag%hfl_trc_ptr(jt)%p_3d,                   &
                    & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT,                         &
                    & t_cf_var('hfl_q'//ctrc, 'kg m-1 s-1',''),                     &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_e)
      ENDDO


      ! vfl_tracer   p_diag%vfl_tracer(nproma,nlevp1,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('vertical_tracer_flux', 'kg m-1 s-1',               &
        &                   'vertical tracer flux')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'vfl_tracer', p_diag%vfl_tracer,               &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_chalf ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

      ALLOCATE(p_diag%vfl_trc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'vfl_tracer',                                   &
                    & 'vfl_q'//ctrc, p_diag%vfl_trc_ptr(jt)%p_3d,                     &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                        &
                    & t_cf_var('vfl_q'//ctrc, 'kg m-1 s-1',''),                    &
                    & t_grib2_var(255,255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_chalf)
      ENDDO


      ! ddt_tracer_adv   p_diag%ddt_tracer_adv(nproma,nlev,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('advective tracer tendency', 'kg kg-1 s-1',         &
        &                   'advective tracer tendency')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_tracer_adv', p_diag%ddt_tracer_adv,       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_c ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

      ALLOCATE(p_diag%ddt_trc_adv_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'ddt_tracer_adv',                               &
                    & 'ddt_adv_q'//ctrc, p_diag%ddt_trc_adv_ptr(jt)%p_3d,             &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                        &
                    & t_cf_var('ddt_adv_q'//ctrc, 'kg kg-1 s-1',''),               &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape3d_c)
      ENDDO


      ! tracer_vi(nproma,nblks_c,3), only Q1, Q2, Q3
      cf_desc    = t_cf_var('tracer_vi', '', 'tracer_vi')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'tracer_vi', p_diag%tracer_vi,                 &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=(/nproma, nblks_c,3/), lrestart=.FALSE. )

      ! tracer_vi_avg(nproma,nblks_c,3), only Q1, Q2, Q3
      cf_desc    = t_cf_var('tracer_vi_avg', '', 'tracer_vi_avg')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'tracer_vi_avg', p_diag%tracer_vi_avg,          &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,   &
                  & ldims=(/nproma, nblks_c,3/) , lrestart=.FALSE.)
    ENDIF



    IF( iforcing== inwp) THEN  !T.R
      ! ddt_tracer_phy   p_diag%ddt_tracer_phy(nproma,nlev,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('physical tracer tendency', 'kg kg-1 s-1',          &
        &                   'physical tracer tendency')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_tracer_phy', p_diag%ddt_tracer_phy,       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_c ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.)

      ALLOCATE(p_diag%ddt_trc_phy_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)') jt
        CALL add_ref( p_diag_list, 'ddt_tracer_phy',                               &
                    & 'ddt_phy_q'//ctrc, p_diag%ddt_trc_phy_ptr(jt)%p_3d,             &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                        &
                    & t_cf_var('ddt_phy_q'//ctrc, 'kg kg-1 s-1',''),               &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape3d_c)
      ENDDO

    ENDIF


    IF( lwrite_extra) THEN
      WRITE(0,*)'inextra_2d=',inextra_2d

      IF(inextra_2d > 0) THEN

        ! extra_2d   p_diag%extra_2d(nproma,nblks_c,inextra_2d)
        !
        cf_desc    = t_cf_var('extra_field_2D', '-', 'extra field 2D')
        grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_diag_list, 'extra_2d', p_diag%extra_2d,                   &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                    & ldims=shape2d_extra )
      ENDIF

      IF(inextra_3d > 0) THEN

        ! extra_3d   p_diag%extra_3d(nproma,nlev,nblks_c,inextra_3d)
        !
        cf_desc    = t_cf_var('extra_fields_3D', '-', 'extra fields 3D')
        grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_diag_list, 'extra_3d', p_diag%extra_3d,                   &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                    & ldims=shape3d_extra )

      ENDIF
    ENDIF

  END SUBROUTINE construct_nh_state_diag_list




  !---------------------------------------------------------------------------
  !>
  !! Allocates all metric coefficients defined in type metrics_3d of the patch.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-04-14)
  !! Modification by Daniel Reinert, DWD, (2010-04-22)
  !! - added geometric height at full levels
  SUBROUTINE construct_nh_metrics_list ( p_patch, p_metrics, p_metrics_list,  &
    &                                    listname )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_metrics),  INTENT(INOUT):: &  !< diagnostic state
      &  p_metrics 

    TYPE(t_var_list), INTENT(INOUT) :: p_metrics_list   !< diagnostic state list

    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e, &    !< number of edge blocks to allocate
               nblks_v       !< number of vertex blocks to allocate

    INTEGER :: nlev, nlevp1

    INTEGER :: shape2d_c(2), shape3d_c(3), shape3d_e(3),               &
      &        shape3d_v(3), shape3d_chalf(3), shape3d_ehalf(3),       &
      &        shape2d_ccubed(3), shape2d_ecubed(3), shape3d_vhalf(3), & 
      &        shape3d_esquared(4) 
    INTEGER :: ientr         !< "entropy" of horizontal slice
    !--------------------------------------------------------------

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c        = (/nproma,          nblks_c    /)     
    shape2d_ccubed   = (/nproma, 3      , nblks_c    /)     
    shape2d_ecubed   = (/nproma, 3      , nblks_e    /)     
    shape3d_c        = (/nproma, nlev   , nblks_c    /)     
    shape3d_chalf    = (/nproma, nlevp1 , nblks_c    /)      
    shape3d_e        = (/nproma, nlev   , nblks_e    /)     
    shape3d_ehalf    = (/nproma, nlevp1 , nblks_e    /)     
    shape3d_esquared = (/2     , nproma , nlev   , nblks_e /)
    shape3d_v        = (/nproma, nlev   , nblks_v    /)     
    shape3d_vhalf    = (/nproma, nlevp1 , nblks_v    /)


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_metrics_list, TRIM(listname) )
    CALL default_var_list_settings( p_metrics_list,            &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    ! geometric height at the vertical interface of cells
    ! z_ifc        p_metrics%z_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('geometric_height_at_half_level_center', 'm',         &
      &                   'geometric height at half level center')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'z_ifc', p_metrics%z_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! geometric height at full levels
    ! z_mc         p_metrics%z_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geometric_height_at_full_level_center', 'm',         &
      &                   'geometric height at full level center')
    grib2_desc = t_grib2_var( 0, 3, 6, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'z_mc', p_metrics%z_mc,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c )


    ! geometric height at full level edges
    ! z_mc_e       p_metrics%z_mc_e(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('geometric_height_at_full_level_edge', 'm',           &
      &                   'geometric height at full level edge')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'z_mc_e', p_metrics%z_mc_e,                   &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_e )


    ! slope of the terrain in normal direction (half level)
    ! ddxn_z_half  p_metrics%ddxn_z_half(nproma,nlevp1,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',             &
      &                   'terrain slope in normal direction')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxn_z_half', p_metrics%ddxn_z_half,         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_ehalf )


    ! slope of the terrain in normal direction (full level)
    ! ddxn_z_full  p_metrics%ddxn_z_full(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',             &
      &                   'terrain slope in normal direction')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxn_z_full', p_metrics%ddxn_z_full,         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_e )


    ! slope of the terrain in tangential direction (half level)
    ! ddxt_z_half  p_metrics%ddxt_z_half(nproma,nlevp1,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_tangential_direction', '-',         &
      &                   'terrain slope in tangential direction')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxt_z_half', p_metrics%ddxt_z_half,         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_ehalf )


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_full  p_metrics%ddqz_z_full(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'ddqz_z_full', p_metrics%ddqz_z_full,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c )


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_full_e  p_metrics%ddqz_z_full_e(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant (edge)')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddqz_z_full_e', p_metrics%ddqz_z_full_e,     &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_e )


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_half  p_metrics%ddqz_z_half(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'ddqz_z_half', p_metrics%ddqz_z_half,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! 1/dz(k-1)-1/dz(k)
    ! diff_1_o_dz  p_metrics%diff_1_o_dz(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('difference_1_over_dz', 'm-1', 'difference 1 over dz')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'diff_1_o_dz', p_metrics%diff_1_o_dz,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! 1/(dz(k-1)*dz(k))
    ! mult_1_o_dz  p_metrics%mult_1_o_dz(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('mult_1_over_dz', 'm-2', 'mult 1 over dz')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'mult_1_o_dz', p_metrics%mult_1_o_dz,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! geopotential at full level cell center
    ! geopot       p_metrics%geopot(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
       &                  'geopotential at full level cell centre')
    grib2_desc = t_grib2_var( 0, 3, 4, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot', p_metrics%geopot,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c )


    ! geopotential at half level cell center
    ! geopot_ifc   p_metrics%geopot_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential at half level cell centre')
    grib2_desc = t_grib2_var( 0, 3, 4, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot_ifc', p_metrics%geopot_ifc,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! geopotential above groundlevel at cell center
    ! geopot_agl   p_metrics%geopot_agl(nproma,nlev  ,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential above groundlevel at cell center')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot_agl', p_metrics%geopot_agl,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c )


    ! geopotential above groundlevel at cell center
    ! geopot_agl_ifc  p_metrics%geopot_agl_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential above groundlevel at cell center')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot_agl_ifc', p_metrics%geopot_agl_ifc,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! geopotential at cell center
    ! dgeopot_mc   p_metrics%dgeopot_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential at cell center')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'dgeopot_mc', p_metrics%dgeopot_mc,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c )


    ! Rayleigh damping
    ! rayleigh_w   p_metrics%rayleigh_w(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('Rayleigh_damping', '-', 'Rayleigh damping')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'rayleigh_w', p_metrics%rayleigh_w,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! Explicit weight in vertical wind solver
    ! vwind_expl_wgt   p_metrics%vwind_expl_wgt(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('Explicit_weight_in_vertical_wind_solver', '-',       &
      &                   'Explicit weight in vertical wind solver')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'vwind_expl_wgt', p_metrics%vwind_expl_wgt,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,        &
                & ldims=shape2d_c )


    ! Implicit weight in vertical wind solver
    ! vwind_impl_wgt  p_metrics%vwind_impl_wgt(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('Implicit_weight_in_vertical_wind_solver', '-',       &
      &                   'Implicit weight in vertical wind solver')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'vwind_impl_wgt', p_metrics%vwind_impl_wgt,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,        &
                & ldims=shape2d_c )



! These fields are needed for triangles only once the initialization in
! mo_nh_testcases is properly rewritten for hexagons
!    IF (i_cell_type== 3) THEN
      ! weighting factor for interpolation from full to half levels
      ! wgtfac_c     p_metrics%wgtfac_c(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for interpolation from full to half levels')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfac_c', p_metrics%wgtfac_c,             &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )


      ! weighting factor for interpolation from full to half levels
      ! wgtfac_e     p_metrics%wgtfac_e(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for interpolation from full to half levels')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfac_e', p_metrics%wgtfac_e,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf )


      ! weighting factor for quadratic interpolation to surface
      ! wgtfacq_c    p_metrics%wgtfacq_c(nproma,3,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for quadratic interpolation to surface')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfacq_c', p_metrics%wgtfacq_c,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_ccubed )


      ! weighting factor for quadratic interpolation to surface
      ! wgtfacq_e    p_metrics%wgtfacq_e(nproma,3,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for quadratic interpolation to surface')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfacq_e', p_metrics%wgtfacq_e,           &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_ecubed )


      ! weighting factor for quadratic interpolation to model top
      ! wgtfacq1_c    p_metrics%wgtfacq1_c(nproma,3,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for quadratic interpolation to model top')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfacq1_c', p_metrics%wgtfacq1_c,         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_ccubed )


      ! weighting factor for quadratic interpolation to model top
      ! wgtfacq1_e   p_metrics%wgtfacq1_e(nproma,3,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for quadratic interpolation to model top')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfacq1_e', p_metrics%wgtfacq1_e,         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_ecubed )


      ! Inverse layer thickness of full levels
      ! inv_ddqz_z_full   p_metrics%inv_ddqz_z_full(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Inverse_layer_thickness', 'm-1',                   &
      &                     'Inverse layer thickness of full levels')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_full', p_metrics%inv_ddqz_z_full, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )


      ! Inverse distance between full levels jk+1 and jk-1
      ! inv_ddqz_z_half2  p_metrics%inv_ddqz_z_half2(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Inverse_distance_between_full_levels', 'm-1',      &
      &                     'Inverse distance between full levels jk+1 and jk-1')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_half2', p_metrics%inv_ddqz_z_half2, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )



    IF (i_cell_type== 3) THEN

      ! Vertical index of neighbor points needed for Taylor-expansion-based pressure gradient
      ! vertidx_gradp  p_metrics%vertidx_gradp(2,nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('Vertical_index', '-',                              &
      &                     'Vertical index')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'vertidx_gradp', p_metrics%vertidx_gradp,   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_esquared, lrestart=.FALSE. )


      ! Height differences between local edge point and neighbor cell points used for
      ! pressure gradient computation
      ! zdiff_gradp  p_metrics%zdiff_gradp(2,nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('Height_differences', 'm',                          &
      &                     'Height differences')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'zdiff_gradp', p_metrics%zdiff_gradp,       &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_esquared, lrestart=.FALSE.  )


      ! Extrapolation factor for Exner pressure
      ! exner_exfac  p_metrics%exner_exfac(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Extrapolation_factor_for_Exner_pressure', '-',     &
      &                     'Extrapolation factor for Exner pressure')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'exner_exfac', p_metrics%exner_exfac,       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )


      ! Reference atmosphere field theta
      ! theta_ref_mc  p_metrics%theta_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_theta', 'K',            &
      &                     'Reference atmosphere field theta')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'theta_ref_mc', p_metrics%theta_ref_mc,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )


      ! Reference atmosphere field theta
      ! theta_ref_ic  p_metrics%theta_ref_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_theta', 'K',            &
      &                     'Reference atmosphere field theta')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'theta_ref_ic', p_metrics%theta_ref_ic,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )


      ! Reference atmosphere field exner
      ! exner_ref_mc  p_metrics%exner_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
      &                     'Reference atmosphere field exner')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'exner_ref_mc', p_metrics%exner_ref_mc,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )


      ! Reference atmosphere field exner
      ! d_exner_dz_ref_ic  p_metrics%d_exner_dz_ref_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
      &                     'Reference atmosphere field exner')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'd_exner_dz_ref_ic', p_metrics%d_exner_dz_ref_ic, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )


      ! Reference atmosphere field exner
      ! d_exner_dz_ref_mc  p_metrics%d_exner_dz_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
      &                     'Reference atmosphere field exner')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'd_exner_dz_ref_mc', p_metrics%d_exner_dz_ref_mc, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )


      ! Reference atmosphere field exner
      ! d2_exner_dz2_ref_mc  p_metrics%d2_exner_dz2_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
      &                     'Reference atmosphere field exner')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'd2_exner_dz2_ref_mc', p_metrics%d2_exner_dz2_ref_mc, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )


      ! Reference atmosphere field rho
      ! rho_refcorr_ic  p_metrics%rho_refcorr_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_rho', '-',              &
      &                     'Reference atmosphere field rho')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'rho_refcorr_ic', p_metrics%rho_refcorr_ic, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )


      ! mask field that excludes boundary halo points
      ! mask_prog_halo_c  p_metrics%mask_prog_halo_c(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('mask_field', '-',                                  &
      &                     'mask field that excludes boundary halo points')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'mask_prog_halo_c', p_metrics%mask_prog_halo_c, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c, lrestart=.FALSE. )


    ELSE IF (i_cell_type== 6) THEN


      ! slope of the coordinate lines in Northern direction
      ! ddnorth_z    p_metrics%ddnorth_z(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('slope_of_coordinate_lines_in_Northern_direction', '-', &
      &                     'slope of coordinate lines in Northern direction')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddnorth_z', p_metrics%ddnorth_z,               &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,           &
                  & ldims=shape3d_v )


      ! slope of the coordinate lines in Eastern direction
      ! ddeast_z     p_metrics%ddeast_z(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('slope_of_coordinate_lines_in_Eastern_direction', '-', &
      &                     'slope of coordinate lines in Eastern direction')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddeast_z', p_metrics%ddeast_z,                &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,          &
                  & ldims=shape3d_v )


      ! functional determinant of the metrics [sqrt(gamma)]
      ! ddqz_z_full_v   p_metrics%ddqz_z_full_v(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',              &
      &                     'metrics functional determinant')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddqz_z_full_v', p_metrics%ddqz_z_full_v,   &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_v )

      ! functional determinant of the metrics [sqrt(gamma)]
      ! ddqz_z_full_r   p_metrics%ddqz_z_full_r(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',              &
      &                     'metrics functional determinant')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'ddqz_z_full_r', p_metrics%ddqz_z_full_r,   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )

      ! functional determinant of the metrics [sqrt(gamma)]
      ! ddqz_z_half_e   p_metrics%ddqz_z_half_e(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',              &
      &                     'metrics functional determinant')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'ddqz_z_half_e', p_metrics%ddqz_z_half_e,   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf )

    ENDIF

  END SUBROUTINE construct_nh_metrics_list


END MODULE mo_nonhydro_state


