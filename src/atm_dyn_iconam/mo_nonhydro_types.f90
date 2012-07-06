#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>
!! Type definition for the dynamical sore of ICONAM.
!!
!! Type definition for the dynamical sore of ICONAM.
!!
!! @author Almut Gassmann (MPI-M)
!! @author Daniel Reinert (DWD-M)
!!
!! @par Revision History
!! Initial release by Daniel Reinert, DWD (2012-02-07)
!! - Moved here from mo_nonhydro_state to avoid circular dependencies
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
MODULE mo_nonhydro_types

  USE mo_kind,                 ONLY: wp
  USE mo_linked_list,          ONLY: t_var_list


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = &
    & '$Id$'

  PUBLIC :: t_nh_prog             ! state vector of prognostic variables (type)
  PUBLIC :: t_nh_diag             ! state vector of diagnostic variables (type)
                                  ! on p- and/or z-levels
  PUBLIC :: t_nh_ref              ! state vector of reference state (type)
  PUBLIC :: t_nh_metrics          ! state vector of metrics variables (type)
  PUBLIC :: t_nh_state            ! state vector of nonhydrostatic variables (type)

  PUBLIC :: t_buffer_memory

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
      tke   (:,:,:)        !! turbulent kinetic energy                         [m^2/s^2]
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
    &  omega_z(:,:,:),      & ! vertical vorticity at dual grid
                              ! (nproma,nlev,nblks_v or nblks_e)               [1/s]
    &  ddt_vn_phy(:,:,:),   & ! normal wind tendency from forcing
                              ! (nproma,nlev,nblks_e)                          [m/s^2]
    &  ddt_exner(:,:,:),    & ! exner pressure tendency from forcing (nproma,nlev,nblks_c)  [1/s]
    &  ddt_exner_phy(:,:,:),& ! exner pressure tendency from physical forcing 
                              ! (nproma,nlev,nblks_c)                     [1/s]
    &  ddt_temp_dyn(:,:,:), & ! rediagnosed temperature tendency from dynamics [K/s]
    &  ddt_tracer_adv(:,:,:,:), &! advective tendency of tracers          [kg/kg/s]
    &  tracer_vi(:,:,:),    & ! vertically integrated tracers(for q1,q2,q3) [kg/m**2]
    &  tracer_vi_avg(:,:,:),& ! average since last output of tracer_vi [kg/m**2]
    &  exner_old(:,:,:),    & ! exner pres from previous step (nproma,nlev,nblks_c)
    &  exner_dyn_incr(:,:,:), & ! exner pres dynamics increment (nproma,nlev,nblks_c)
    &  temp(:,:,:),         & ! temperature (nproma,nlev,nblks_c)                 [K]
    &  tempv(:,:,:),        & ! virtual temperature (nproma,nlev,nblks_c)         [K]
    &  temp_ifc(:,:,:),     & ! temperature at half levels (nproma,nlevp1,nblks_c)[K]
    &  pres(:,:,:),         & ! pressure (nproma,nlev,nblks_c)                  [Pa]
    &  pres_ifc(:,:,:),     & ! pressure at interfaces (nproma,nlevp1,nblks_c)  [Pa]
    &  pres_sfc(:,:),       & ! diagnosed surface pressure (nproma,nblks_c)     [Pa]
    &  pres_sfc_old(:,:),   & ! diagnosed surface pressure at previous timestep (nproma,nblks_c) [Pa]
    &  pres_sfc_s6avg(:,:), & ! 6 hourly sample  surface pressure average       [Pa]
    &  pres_msl(:,:),       & ! diagnosed mean sea level pressure (nproma,nblks_c)  [Pa]
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

    !
    ! c) variables needed for the hexagonal grid only
    &  e_kin(:,:,:),        & ! spec. kinetic energy (nproma,nlev,nblks_c) [m^2/s^2]
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
    &  ddt_vn(:,:,:),       & ! normal wind tendency from forcing              [m/s^2]
    &  ddt_vn_vort(:,:,:),  & ! normal wind tendency from vorticity flux term
                              ! (nproma,nlev,nblks_e,1:3)                    [m/s^2]
    &  ddt_w(:,:,:),        & ! vert. wind tendency from forcing
                              ! (nproma,nlevp1,nblks_c)                        [m/s^2]
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

    TYPE(t_ptr_nh),ALLOCATABLE :: ddt_vn_adv_ptr(:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_nh),ALLOCATABLE :: ddt_w_adv_ptr (:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_nh),ALLOCATABLE :: q_int_ptr     (:)
    TYPE(t_ptr_nh),ALLOCATABLE :: q_ubc_ptr     (:)

    TYPE(t_ptr_nh),ALLOCATABLE :: tracer_vi_ptr(:)      !< pointer array: one pointer for each tracer
    TYPE(t_ptr_nh),ALLOCATABLE :: tracer_vi_avg_ptr(:)  !< pointer array: one pointer for each tracer
    TYPE(t_ptr_nh),ALLOCATABLE :: extra_2d_ptr(:)
    TYPE(t_ptr_nh),ALLOCATABLE :: extra_3d_ptr(:)

  END TYPE t_nh_diag


  TYPE t_nh_ref
    REAL(wp), POINTER ::    &
      vn_ref(:,:,:),        & !! orthogonal normal wind (nproma,nlev,nblks_e)         [m/s]
      w_ref(:,:,:)             !> orthogonal vertical wind (nproma,nlevp1,nblks_c)    [m/s]
  END TYPE t_nh_ref


  TYPE t_nh_metrics

   ! a) Variables needed for triangles and hexagons
   !
   ! geometric height at the vertical interface of cells (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: z_ifc(:,:,:)
   ! geometric height at full levels (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: z_mc(:,:,:)

   ! slope of the terrain in normal direction (nproma,nlevp1,nblks_e)
   REAL(wp), POINTER :: ddxn_z_half(:,:,:)
   ! slope of the terrain in normal direction (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddxn_z_full(:,:,:)
   ! slope of the terrain in tangential direction (nproma,nlevp1,nblks_e)
   REAL(wp), POINTER :: ddxt_z_half(:,:,:)
   ! slope of the terrain in tangential direction (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddxt_z_full(:,:,:)

   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: ddqz_z_full(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddqz_z_full_e(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddqz_z_full_r(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: ddqz_z_half(:,:,:)

   ! geopotential at cell center (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: geopot(:,:,:)
   ! geopotential above ground level at cell center (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: geopot_agl(:,:,:)
   ! geopotential above ground level at interfaces and cell center (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: geopot_agl_ifc(:,:,:)
   ! geopotential at cell center (nproma,nlev,nblks_c)
   REAL(wp), POINTER :: dgeopot_mc(:,:,:)


   ! Rayleigh damping on the vertical velocity
   REAL(wp), POINTER :: rayleigh_w(:)
   ! Rayleigh damping on the normal velocity
   REAL(wp), POINTER :: rayleigh_vn(:)
   ! Enhancement factor for nabla4 background diffusion
   REAL(wp), POINTER :: enhfac_diffu(:)

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
   ! Interpolation coefficients for cubic interpolation of Exner pressure (8,nproma,nlev,nblks_e)
   REAL(wp), POINTER :: coeff_gradp(:,:,:,:)
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
   REAL(wp), POINTER :: tsfc_ref(:,:)
   REAL(wp), POINTER :: exner_ref_mc(:,:,:)
   REAL(wp), POINTER :: rho_ref_mc  (:,:,:)
   REAL(wp), POINTER :: d_exner_dz_ref_ic(:,:,:)
   REAL(wp), POINTER :: d2dexdz2_fac1_mc(:,:,:)
   REAL(wp), POINTER :: d2dexdz2_fac2_mc(:,:,:)
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
   ! Correction term needed to use perturbation density for lateral boundary nudging
   ! (note: this field is defined on the local parent grid in case of MPI parallelization)
   REAL(wp), POINTER :: rho_ref_corr(:,:,:)
   ! Area of subdomain for which feedback is performed; dim: (nlev)
   REAL(wp), POINTER :: fbk_dom_volume(:)

   ! c) Variables needed for the hexagonal grid only
   !
   ! geometric height at full level edges (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: z_mc_e(:,:,:)

   ! slope of the coordinate lines towards the North (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddnorth_z(:,:,:)
   ! slope of the coordinate lines towards the East (nproma,nlev,nblks_e)
   REAL(wp), POINTER :: ddeast_z(:,:,:)

   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlev,nblks_v)
   REAL(wp), POINTER :: ddqz_z_full_v(:,:,:)
   ! functional determinant of the metrics [sqrt(gamma)] (nproma,nlevp1,nblks_e)
   REAL(wp), POINTER :: ddqz_z_half_e(:,:,:)

   ! 1/dz(k-1)-1/dz(k) (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: diff_1_o_dz(:,:,:)
   ! 1/(dz(k-1)*dz(k)) (nproma,nlevp1,nblks_c)
   REAL(wp), POINTER :: mult_1_o_dz(:,:,:)

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


!-------------------------------------------------------------------------
!                      STATE VECTORS AND LISTS
!-------------------------------------------------------------------------
  TYPE t_nh_state

    !array of prognostic states at different timelevels
    TYPE(t_nh_prog),  ALLOCATABLE :: prog(:)       !< shape: (timelevels)
    TYPE(t_var_list), ALLOCATABLE :: prog_list(:)  !< shape: (timelevels)

    TYPE(t_nh_diag)    :: diag
    TYPE(t_var_list)   :: diag_list

    TYPE(t_nh_ref)     :: ref
    TYPE(t_var_list)   :: ref_list

    TYPE(t_nh_metrics) :: metrics
    TYPE(t_var_list)   :: metrics_list

    TYPE(t_var_list), ALLOCATABLE :: tracer_list(:) !< shape: (timelevels)

  END TYPE t_nh_state


END MODULE mo_nonhydro_types




