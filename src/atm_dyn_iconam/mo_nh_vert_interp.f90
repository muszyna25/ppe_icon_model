!>
!! This module contains routines for the vertical interpolation of
!! atmospheric data provided by external analyses to the ICON grid
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-06-29)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_vert_interp

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag,  &
    &                               t_nh_metrics
  USE mo_opt_diagnostics,     ONLY: t_vcoeff, vcoeff_allocate
  USE mo_intp_data_strc,      ONLY: t_int_state, p_int_state
  USE mo_intp,                ONLY: edges2cells_scalar, cell_avg, &
    &                               cells2edges_scalar, cells2verts_scalar
  USE mo_parallel_config,     ONLY: nproma
  USE mo_physical_constants,  ONLY: grav, rd, rdv, o_m_rdv, dtdz_standardatm
  USE mo_grid_config,         ONLY: n_dom
  USE mo_run_config,          ONLY: iforcing, num_lev
  USE mo_io_config,           ONLY: itype_pres_msl
  USE mo_impl_constants,      ONLY: inwp, iecham, PRES_MSL_METHOD_GME, PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_DWD, &
    &                               PRES_MSL_METHOD_IFS_CORR, MODE_IFSANA, MODE_COMBINED, MODE_ICONVREMAP
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_initicon_config,     ONLY: zpbl1, zpbl2, l_coarse2fine_mode, init_mode, lread_vn, lvert_remap_fg
  USE mo_initicon_types,      ONLY: t_init_state, t_initicon_state
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_nh_init_utils,       ONLY: interp_uv_2_vn, adjust_w, convert_thdvars, &
                                    convert_omega2w
  USE mo_util_phys,           ONLY: virtual_temp, vap_pres
  USE mo_math_gradients,      ONLY: grad_fd_norm, grad_fd_tang
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf, interpol2_vec_grf
  USE mo_sync,                ONLY: sync_patch_array, SYNC_C, SYNC_E, sync_patch_array_mult
  USE mo_satad,               ONLY: sat_pres_water
  USE mo_nwp_sfc_interp,      ONLY: process_sfcfields

  IMPLICIT NONE
  PRIVATE


  ! Threshold for switching between analytical formulas for constant temperature and
  ! constant vertical gradient of temperature, respectively
  REAL(wp), PARAMETER :: dtvdz_thresh = 1.e-4_wp ! 0.1 K/km

  ! Artificial limits on temperature profile used for extrapolation below the ground
  REAL(wp), PARAMETER :: t_low  = 255.0_wp
  REAL(wp), PARAMETER :: t_high = 290.5_wp

  ! Height above ground level from where the temperature is used for downward extrapolation
  REAL(wp), PARAMETER :: zagl_extrap_1 = 10._wp
  REAL(wp), PARAMETER :: topo_extrap_1 = 1500._wp  ! 10 m AGL valid for grid points below 1500 m
  REAL(wp), PARAMETER :: zagl_extrap_2 = 150._wp
  REAL(wp), PARAMETER :: trans_depth   = 500._wp   ! depth of transition zone 500 m


  PUBLIC :: vert_interp
  PUBLIC :: vert_interp_atm
  PUBLIC :: vert_interp_sfc
  PUBLIC :: prepare_lin_intp
  PUBLIC :: prepare_vert_interp_z
  PUBLIC :: prepare_vert_interp_p
  PUBLIC :: prepare_vert_interp_i
  PUBLIC :: prepare_extrap, prepare_extrap_ifspp
  PUBLIC :: lin_intp, uv_intp, qv_intp, temperature_intp, z_at_plevels
  PUBLIC :: diagnose_pmsl, diagnose_pmsl_gme, diagnose_pmsl_ifs


CONTAINS



  !-------------
  !>
  !! SUBROUTINE vert_interp_atm
  !! Outer driver routine for vertical interpolation of analysis
  !! data (atmosphere only) interpolated horizontally by IFS2ICON to
  !! the ICON grid
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE vert_interp_atm(p_patch, p_nh_state, p_int, p_grf, initicon)

    TYPE(t_patch),          INTENT(IN)       :: p_patch(:)
    TYPE(t_nh_state),       INTENT(IN)       :: p_nh_state(:)
    TYPE(t_int_state),      INTENT(IN)       :: p_int(:)
    TYPE(t_gridref_state),  INTENT(IN)       :: p_grf(:)
    TYPE(t_initicon_state), INTENT(INOUT)    :: initicon(:)

    ! LOCAL VARIABLES
    INTEGER :: jg, jn, jgc

!-------------------------------------------------------------------------

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE ! skip model domains not active at initial time

      WRITE(message_text,'(a,i4)') 'Vertical interpolation of analysis data, domain ',jg
      CALL message('vert_interp_atm', TRIM(message_text))

      IF (p_patch(jg)%n_patch_cells==0) CYCLE ! skip empty patches

      CALL vert_interp(p_patch(jg), p_int(jg), p_nh_state(jg)%metrics, initicon(jg))

      ! Apply boundary interpolation for u and v because the outer nest boundary
      ! points would remain undefined otherwise
      DO jn = 1, p_patch(jg)%n_childdom

        jgc = p_patch(jg)%child_id(jn)

        CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_grf(jg)%p_dom(jn), &
                                1, initicon(jg)%atm%w, initicon(jgc)%atm%w )

        CALL interpol2_vec_grf (p_patch(jg), p_patch(jgc), p_grf(jg)%p_dom(jn), &
                                1, initicon(jg)%atm%vn, initicon(jgc)%atm%vn )

      ENDDO
    ENDDO

  END SUBROUTINE vert_interp_atm



  !-------------
  !>
  !! SUBROUTINE vert_interp_sfc
  !! Outer driver routine for vertical interpolation of analysis
  !! data (surface only) interpolated horizontally by ICONREMAP to
  !! the ICON grid
  !!
  !! @par Revision History
  !! Initial version by Daniel Reinert, DWD(2012-12-19)
  !!
  !!
  SUBROUTINE vert_interp_sfc(p_patch, initicon)

    TYPE(t_patch),       INTENT(IN)       :: p_patch(:)
    CLASS(t_init_state), INTENT(INOUT)    :: initicon(:)

    ! LOCAL VARIABLES
    INTEGER :: jg

    !-------------------------------------------------------------------------

    DO jg = 1, n_dom

      IF (p_patch(jg)%n_patch_cells==0) CYCLE ! skip empty patches
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE ! skip model domains not active at initial time

      !
      ! process surface fields
      !
      CALL process_sfcfields(p_patch(jg), initicon(jg))
    ENDDO

  END SUBROUTINE vert_interp_sfc


  !-------------
  !>
  !! SUBROUTINE vert_interp
  !! Domain-wise driver routine for vertical interpolation of analysis
  !! data (atmosphere only) interpolated horizontally by ICONREMAP to
  !! the ICON grid
  !!
  !! SOURCE data taken from initicon%atm_in: 
  !! z3d (psfc or log(psfc), phi_sfc), temp, vn (u,v), qv, qc, qi, qr, qs, pres, w (omega), tke
  !! Alternatives are given in brackets.
  !! - if z3d is not available, it is derived from (psfc or log(psfc), phi_sfc)
  !! - if omega rather than w is given, it is converted to w prior to the interpolation step.
  !!
  !! Fields are vertically interpolated onto the ICON vertical grid 
  !! and stored in the initicon%atm state. They are finally converted into the following 
  !! set of prognostic variables:
  !! vn, w, qv, qc, qi, qr, qs, rho, exner, theta_v
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE vert_interp(p_patch, p_int, p_metrics, initicon, opt_use_vn)

    TYPE(t_patch),          INTENT(IN)       :: p_patch
    TYPE(t_int_state),      INTENT(IN)       :: p_int
    TYPE(t_nh_metrics),     INTENT(IN)       :: p_metrics

    ! (Note: This is a CLASS parameter, st. we can feed
    ! t_initicon_state into this subroutine)
    CLASS(t_init_state),    INTENT(INOUT)    :: initicon

    LOGICAL, OPTIONAL,      INTENT(IN)       :: opt_use_vn


    ! LOCAL VARIABLES
    CHARACTER(LEN=*), PARAMETER :: routine = 'vert_interp'

    INTEGER :: jb
    INTEGER :: nlen, nlev, nlevp1
    INTEGER :: nlev_in         ! number of vertical levels in source vgrid 
    LOGICAL :: lc2f, l_use_vn

    ! Auxiliary fields for input data
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev,p_patch%nblks_c) :: temp_v_in

    ! Auxiliary field for output data
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: &
      z_tempv

    ! CELLS !
    !
    ! Auxiliary fields for coefficients
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: &
      wfac_lin, coef1, coef2, coef3

    REAL(wp), DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_c) :: &
      wfac_lin_w

    INTEGER , DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: &
      idx0_lin, idx0_cub

    INTEGER , DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_c) :: &
      idx0_lin_w

    ! extrapolation weights and indices for cells
    REAL(wp), DIMENSION(nproma,p_patch%nblks_c) :: &
      wfacpbl1, wfacpbl2, slope

    INTEGER , DIMENSION(nproma,p_patch%nblks_c) :: &
      bot_idx_lin, bot_idx_lin_w, bot_idx_cub, kpbl1, kpbl2



    ! EDGES !
    !
    ! Auxiliary fields for coefficients at edges
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: &
      wfac_lin_e, coef1_e, coef2_e, coef3_e

    INTEGER , DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: &
      idx0_lin_e, idx0_cub_e

    ! extrapolation weights
    REAL(wp), DIMENSION(nproma,p_patch%nblks_e) :: &
      wfacpbl1_e, wfacpbl2_e

    INTEGER , DIMENSION(nproma,p_patch%nblks_e) :: &
      bot_idx_lin_e, bot_idx_cub_e, kpbl1_e, kpbl2_e

    ! full level heights of initicon input fields at edge midpoints
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev,p_patch%nblks_e) :: atm_in_z_me

    ! full level heights of ICON vertical grid at edge midpoints
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_me

!-------------------------------------------------------------------------

    nlev_in = initicon%atm_in%nlev

    IF (nlev_in == 0) THEN
       CALL finish(routine, "Number of input levels <nlev_in> not yet initialized.")
    END IF

    IF (PRESENT(opt_use_vn)) THEN
      l_use_vn = opt_use_vn
    ELSE
      l_use_vn = lread_vn ! use vn field if available as input from file
    ENDIF

    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Switch to determine if corrections for coarse-to-fine-mesh interpolation are to be used
    lc2f   = l_coarse2fine_mode(p_patch%id)

    ! Compute virtual temperature of input data
    CALL virtual_temp(p_patch, initicon%atm_in%temp, initicon%atm_in%qv, initicon%atm_in%qc, &
                      initicon%atm_in%qi, initicon%atm_in%qr, initicon%atm_in%qs,            &
                      temp_v=temp_v_in)

    ! Prepare interpolation coefficients for cells
    IF (.NOT. ASSOCIATED(initicon%const%z_mc)) THEN
       CALL finish(routine, "Internal error!")
    END IF
    IF (.NOT. ASSOCIATED(initicon%const%z_mc_in)) THEN
       CALL finish(routine, "Internal error: z_mc_in not associated!")
    END IF
    CALL prepare_lin_intp(initicon%const%z_mc_in, initicon%const%z_mc,      &
                          p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                          wfac_lin, idx0_lin, bot_idx_lin)

    CALL prepare_extrap(initicon%const%z_mc_in,                     &
                        p_patch%nblks_c, p_patch%npromz_c, nlev_in, &
                        kpbl1, wfacpbl1, kpbl2, wfacpbl2 )



    CALL prepare_cubic_intp(initicon%const%z_mc_in, initicon%const%z_mc,      &
                            p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                            coef1, coef2, coef3, idx0_cub, bot_idx_cub)


    ! Perform vertical interpolation

    ! Temperature
    CALL compute_slope(p_patch, p_int, initicon%const%topography_c, slope)

    CALL temperature_intp(initicon%atm_in%temp, initicon%atm%temp,                &
                          initicon%const%z_mc_in, initicon%const%z_mc,            &
                          p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,       &
                          coef1, coef2, coef3, wfac_lin,                          &
                          idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,           &
                          wfacpbl1, kpbl1, wfacpbl2, kpbl2,                       &
                          l_restore_sfcinv=.TRUE., l_hires_corr=lc2f,             &
                          extrapol_dist=-1500._wp, l_pz_mode=.FALSE., slope=slope )


    ! horizontal wind components
    IF (l_use_vn) THEN

      !
      ! prepare interpolation coefficients for edges
      !

      ! full level heights at edges for initicon input
      CALL cells2edges_scalar(initicon%const%z_mc_in, p_patch, p_int%c_lin_e, atm_in_z_me, opt_fill_latbc=.TRUE.)


      ! full level heights at edges for ICON vertical grid
      CALL cells2edges_scalar(p_metrics%z_mc, p_patch, p_int%c_lin_e, z_me, opt_fill_latbc=.TRUE.)


      ! compute extrapolation coefficients for edges
      CALL prepare_extrap(atm_in_z_me, p_patch%nblks_e, p_patch%npromz_e, &
        &                 nlev_in, kpbl1_e, wfacpbl1_e, kpbl2_e, wfacpbl2_e )


      ! compute linear interpolation coefficients for edges
      CALL prepare_lin_intp(atm_in_z_me, z_me,                                &
                            p_patch%nblks_e, p_patch%npromz_e, nlev_in, nlev, &
                            wfac_lin_e, idx0_lin_e, bot_idx_lin_e)

      ! compute cubic interpolation coefficients for edges
      CALL prepare_cubic_intp(atm_in_z_me, z_me,                                &
                            p_patch%nblks_e, p_patch%npromz_e, nlev_in, nlev,   &
                            coef1_e, coef2_e, coef3_e, idx0_cub_e, bot_idx_cub_e)

      ! vertically interpolate vn to ICON grid
      CALL uv_intp(initicon%atm_in%vn, initicon%atm%vn,                  &
                   atm_in_z_me, z_me,                                    &
                   p_patch%nblks_e, p_patch%npromz_e, nlev_in, nlev,     &
                   coef1_e, coef2_e, coef3_e, wfac_lin_e,                &
                   idx0_cub_e, idx0_lin_e, bot_idx_cub_e, bot_idx_lin_e, &
                   wfacpbl1_e, kpbl1_e, wfacpbl2_e, kpbl2_e,             &
                   l_hires_intp=lc2f                                     )

    ELSE
      CALL uv_intp(initicon%atm_in%u, initicon%atm%u,                &
                   initicon%const%z_mc_in, initicon%const%z_mc,      &
                   p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                   coef1, coef2, coef3, wfac_lin,                    &
                   idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                   wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                   l_hires_intp=lc2f                                 )
      CALL uv_intp(initicon%atm_in%v, initicon%atm%v,                &
                   initicon%const%z_mc_in, initicon%const%z_mc,      &
                   p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                   coef1, coef2, coef3, wfac_lin,                    &
                   idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                   wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                   l_hires_intp=lc2f                                 )

      ! Convert u and v on cell points to vn at edge points
      CALL interp_uv_2_vn(p_patch, p_int, initicon%atm%u, initicon%atm%v, initicon%atm%vn)
    ENDIF

    CALL sync_patch_array(SYNC_E,p_patch,initicon%atm%vn)


    ! Preliminary interpolation of QV: this is needed to compute virtual temperature
    ! below, which in turn is required to integrate the hydrostatic equation
    ! A lower limit of 2.5 ppm is imposed
    CALL lin_intp(initicon%atm_in%qv, initicon%atm%qv,                 &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,    &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,    &
                  wfacpbl2, kpbl2, l_loglin=.TRUE., l_extrapol=.TRUE., &
                  l_pd_limit=.TRUE., lower_limit=2.5e-6_wp             )

    ! Cloud and precipitation variables - linear interpolation only because cubic may
    ! cause negative values, and no-gradient condition for downward extrapolation
    ! Positive definite limiter is used in all cases
    CALL lin_intp(initicon%atm_in%qc, initicon%atm%qc,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE.)

    CALL lin_intp(initicon%atm_in%qi, initicon%atm%qi,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE.)

    CALL lin_intp(initicon%atm_in%qr, initicon%atm%qr,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE.)

    CALL lin_intp(initicon%atm_in%qs, initicon%atm%qs,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE.)

    ! Compute virtual temperature with preliminary QV
    CALL virtual_temp(p_patch, initicon%atm%temp, initicon%atm%qv, initicon%atm%qc,    &
                      initicon%atm%qi, initicon%atm%qr, initicon%atm%qs, temp_v=z_tempv)

    ! Interpolate pressure on ICON grid
    CALL pressure_intp_initmode(initicon%atm_in%pres, temp_v_in, initicon%const%z_mc_in,    &
      &                        initicon%atm%pres, z_tempv, initicon%const%z_mc,&
      &                        p_patch%nblks_c, p_patch%npromz_c, nlev,        &
      &                        wfac_lin, idx0_lin, bot_idx_lin)


    CALL qv_intp(initicon%atm_in%qv, initicon%atm%qv, initicon%const%z_mc_in,  &
                 initicon%const%z_mc, initicon%atm_in%temp, initicon%atm_in%pres, &
                 initicon%atm%temp, initicon%atm%pres,                      &
                 p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,          &
                 coef1, coef2, coef3, wfac_lin,                             &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,              &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_satlimit=.TRUE.,       &
                 lower_limit=2.5e-6_wp, l_restore_pbldev=.TRUE.,            &
                 opt_qc=initicon%atm%qc ) ! ... for consistency checking

    ! Compute virtual temperature with final QV
    CALL virtual_temp(p_patch, initicon%atm%temp, initicon%atm%qv, initicon%atm%qc,    &
                      initicon%atm%qi, initicon%atm%qr, initicon%atm%qs, temp_v=z_tempv)

    ! Final interpolation of pressure on ICON grid
    CALL pressure_intp_initmode(initicon%atm_in%pres, temp_v_in, initicon%const%z_mc_in, &
      &                        initicon%atm%pres, z_tempv, initicon%const%z_mc,          &
      &                        p_patch%nblks_c, p_patch%npromz_c, nlev,                  &
      &                        wfac_lin, idx0_lin, bot_idx_lin)


    ! Convert thermodynamic variables into set of NH prognostic variables
    CALL convert_thdvars(p_patch, initicon%atm%pres, z_tempv,    &
                         initicon%atm%rho, initicon%atm%exner,   &
                         initicon%atm%theta_v                    )


      ! Compute coefficients for w interpolation
      CALL prepare_lin_intp(initicon%const%z_mc_in, initicon%const%z_ifc,       &
                            p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlevp1, &
                            wfac_lin_w, idx0_lin_w, bot_idx_lin_w)

      ! Perform linear interpolation of w
      ! Note: the coefficients for gradient computation (*pbl*) do not have to be recomputed
      ! because of l_extrapol=.FALSE.,
      CALL lin_intp(initicon%atm_in%w, initicon%atm%w,                     &
                    p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlevp1,    &
                    wfac_lin_w, idx0_lin_w, bot_idx_lin_w, wfacpbl1, kpbl1,&
                    wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                    l_pd_limit=.FALSE.)

    ! Impose appropriate lower boundary condition on vertical wind field
    CALL adjust_w(p_patch, p_int, initicon%atm%vn, initicon%const%z_ifc, initicon%atm%w)

    CALL sync_patch_array(SYNC_C,p_patch,initicon%atm%w)


    IF ((init_mode == MODE_ICONVREMAP .OR. lvert_remap_fg) .AND. ASSOCIATED(initicon%atm_in%tke)) THEN
      CALL lin_intp(initicon%atm_in%tke, initicon%atm%tke,                 &
                    p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlevp1,    &
                    wfac_lin_w, idx0_lin_w, bot_idx_lin_w, wfacpbl1, kpbl1,&
                    wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                    l_pd_limit=.TRUE., lower_limit=1.e-5_wp)
    ENDIF

  END SUBROUTINE vert_interp




  !-------------
  !>
  !! Compute coefficients for vertical interpolation of model-level
  !! fields to constant-height levels.
  !!
  !! @note Computation of coefficients already performs interpolation
  !! of temperature, QV (prognostic+diagnostic) and pressure onto z-levels.
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-18)
  !! Modification by Daniel Reinert, DWD (2011-09-05)
  !! - routine duplicated for use in ICON
  !! Modification by F. Prill, DWD (2012-02-29)
  !! - Separated coefficient computation from interpolation
  !!
  SUBROUTINE prepare_vert_interp_z(p_patch, p_diag, p_metrics, intp_hrz, nzlev,  &
    &                              temp_z_out, pres_z_out, p_z3d_out, vcoeff_z)

    TYPE(t_patch),        TARGET,      INTENT(IN)    :: p_patch
    TYPE(t_nh_diag),      POINTER                    :: p_diag
    TYPE(t_nh_metrics),                INTENT(IN)    :: p_metrics
    TYPE(t_int_state),                 INTENT(IN)    :: intp_hrz          !< pointer to data structure for interpolation
    INTEGER,                           INTENT(IN)    :: nzlev             !< number of output levels (height)
    REAL(wp),                          INTENT(INOUT) :: temp_z_out(:,:,:),          &
      &                                                 pres_z_out(:,:,:)
    REAL(wp) ,            POINTER                    :: p_z3d_out(:,:,:)  !< height field
    TYPE(t_vcoeff),                    INTENT(INOUT) :: vcoeff_z

    ! LOCAL VARIABLES
!!$    CHARACTER(*), PARAMETER :: routine = TRIM("mo_nh_vert_interp:prepare_vert_interp_z")

    INTEGER :: nlev, nlevp1, jg, i_endblk, nblks_c,nblks_e, npromz_c,npromz_e !< blocking parameters
    ! Auxiliary field for output data
    REAL(wp), DIMENSION(nproma,nzlev,p_patch%nblks_c)        :: z_auxz
    REAL(wp), DIMENSION(nproma,nzlev,p_patch%nblks_e)        :: p_z3d_edge
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_me
    ! Pointer to virtual temperature / temperature, depending on whether the run is moist or dry
    REAL(wp), POINTER, DIMENSION(:,:,:) :: ptr_tempv

    !-------------------------------------------------------------------------

    vcoeff_z%l_initialized = .TRUE.
    IF (p_patch%n_patch_cells==0) RETURN

    nlev     = p_patch%nlev
    nlevp1   = p_patch%nlevp1
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c
    nblks_e  = p_patch%nblks_e
    npromz_e = p_patch%npromz_e
    jg       = p_patch%id

    !--- Coefficients: Interpolation to z-level fields

    ! allocate coefficient table:
    CALL vcoeff_allocate(nblks_c, nblks_e, nzlev, vcoeff_z)

    ! Prepare interpolation coefficients
    CALL prepare_lin_intp(p_metrics%z_mc, p_z3d_out,                                         & !in
      &                   nblks_c, npromz_c, nlev, nzlev,                                    & !in
      &                   vcoeff_z%lin_cell%wfac_lin,                                        & !out
      &                   vcoeff_z%lin_cell%idx0_lin,                                        & !out
      &                   vcoeff_z%lin_cell%bot_idx_lin  )                                     !out
    IF ( ANY((/PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR, PRES_MSL_METHOD_DWD/)== itype_pres_msl) ) THEN
      CALL prepare_extrap_ifspp(p_metrics%z_ifc, p_metrics%z_mc, nblks_c, npromz_c, nlev,    & !in
        &                       vcoeff_z%lin_cell%kpbl1, vcoeff_z%lin_cell%zextrap,          & !out
        &                       vcoeff_z%lin_cell%wfacpbl1)                                    !out
      vcoeff_z%lin_cell%kpbl2(:,:) = 0
      vcoeff_z%lin_cell%wfacpbl2(:,:) = 0._wp
    ELSE
      CALL prepare_extrap(p_metrics%z_mc, nblks_c, npromz_c, nlev,                           & !in
        &                 vcoeff_z%lin_cell%kpbl1, vcoeff_z%lin_cell%wfacpbl1,               & !out
        &                 vcoeff_z%lin_cell%kpbl2, vcoeff_z%lin_cell%wfacpbl2 )                !out
    ENDIF
    CALL prepare_cubic_intp(p_metrics%z_mc, p_z3d_out, nblks_c, npromz_c, nlev, nzlev,       & !in
      &                     vcoeff_z%cub_cell%coef1, vcoeff_z%cub_cell%coef2,                & !out
      &                     vcoeff_z%cub_cell%coef3,                                         & !out
      &                     vcoeff_z%cub_cell%idx0_cub, vcoeff_z%cub_cell%bot_idx_cub )        !out

    ! Perform vertical interpolation
    CALL temperature_intp(p_diag%temp, z_auxz, p_metrics%z_mc, p_z3d_out,                   & !in,out
      &                   nblks_c, npromz_c, nlev, nzlev,                                   & !in
      &                   vcoeff_z%cub_cell%coef1, vcoeff_z%cub_cell%coef2,                 & !in
      &                   vcoeff_z%cub_cell%coef3, vcoeff_z%lin_cell%wfac_lin,              & !in
      &                   vcoeff_z%cub_cell%idx0_cub, vcoeff_z%lin_cell%idx0_lin,           & !in
      &                   vcoeff_z%cub_cell%bot_idx_cub, vcoeff_z%lin_cell%bot_idx_lin,     & !in
      &                   vcoeff_z%lin_cell%wfacpbl1, vcoeff_z%lin_cell%kpbl1,              & !in
      &                   vcoeff_z%lin_cell%wfacpbl2, vcoeff_z%lin_cell%kpbl2,              & !in
      &                   l_restore_sfcinv=.FALSE., l_hires_corr=.FALSE.,                   & !in
      &                   extrapol_dist=0._wp, l_pz_mode=.TRUE., zextrap=vcoeff_z%lin_cell%zextrap) !in


    IF (jg > 1) THEN ! copy outermost nest boundary row in order to avoid missing values
      i_endblk = p_patch%cells%end_blk(1,1)
      temp_z_out(:,:,1:i_endblk) = z_auxz(:,:,1:i_endblk)
    ENDIF
    CALL cell_avg(z_auxz, p_patch, p_int_state(jg)%c_bln_avg, temp_z_out)

    IF (  iforcing == inwp .OR. iforcing == iecham  ) THEN
      ptr_tempv => p_diag%tempv(:,:,:)
    ELSE
      ptr_tempv => p_diag%temp(:,:,:)
    END IF

    ! Interpolate pressure on z-levels
    CALL pressure_intp(p_diag%pres, ptr_tempv, p_metrics%z_mc,                              & !in
      &                z_auxz, p_z3d_out,                                                   & !out,in
      &                nblks_c, npromz_c, nlev, nzlev,                                      & !in
      &                vcoeff_z%lin_cell%wfac_lin,    vcoeff_z%lin_cell%idx0_lin,           & !in
      &                vcoeff_z%lin_cell%bot_idx_lin, vcoeff_z%lin_cell%wfacpbl1,           & !in
      &                vcoeff_z%lin_cell%kpbl1, vcoeff_z%lin_cell%wfacpbl2,                 & !in
      &                vcoeff_z%lin_cell%kpbl2, vcoeff_z%lin_cell%zextrap)                    !in

    IF (jg > 1) THEN ! copy outermost nest boundary row in order to avoid missing values
      i_endblk = p_patch%cells%end_blk(1,1)
      pres_z_out(:,:,1:i_endblk) = z_auxz(:,:,1:i_endblk)
    ENDIF
    CALL cell_avg(z_auxz, p_patch, p_int_state(jg)%c_bln_avg, pres_z_out)

    IF ( ANY((/PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR, PRES_MSL_METHOD_DWD/)== itype_pres_msl) ) THEN
      CALL prepare_extrap(p_metrics%z_mc, nblks_c, npromz_c, nlev,                           & !in
        &                 vcoeff_z%lin_cell%kpbl1, vcoeff_z%lin_cell%wfacpbl1,               & !out
        &                 vcoeff_z%lin_cell%kpbl2, vcoeff_z%lin_cell%wfacpbl2 )                !out
    ENDIF

    !--- Interpolation data for the vertical interface of cells, "nlevp1"

    CALL prepare_lin_intp(p_metrics%z_ifc, p_z3d_out, nblks_c, npromz_c, nlevp1, nzlev,     & !in
      &                   vcoeff_z%lin_cell_nlevp1%wfac_lin,                                & !out
      &                   vcoeff_z%lin_cell_nlevp1%idx0_lin,                                & !out
      &                   vcoeff_z%lin_cell_nlevp1%bot_idx_lin  )                             !out
    CALL prepare_extrap(p_metrics%z_ifc, nblks_c, npromz_c, nlevp1,                         & !in
      &                 vcoeff_z%lin_cell_nlevp1%kpbl1, vcoeff_z%lin_cell_nlevp1%wfacpbl1,  & !out
      &                 vcoeff_z%lin_cell_nlevp1%kpbl2, vcoeff_z%lin_cell_nlevp1%wfacpbl2 )   !out

    !--- Compute data for interpolation of edge-based fields

    ! Compute geometric height at edge points
    CALL cells2edges_scalar(p_metrics%z_mc, p_patch, intp_hrz%c_lin_e, z_me, opt_fill_latbc=.TRUE.)
    CALL cells2edges_scalar(p_z3d_out, p_patch, intp_hrz%c_lin_e, p_z3d_edge, opt_fill_latbc=.TRUE.)
    CALL sync_patch_array_mult(SYNC_E,p_patch,2,z_me,p_z3d_edge)

    CALL prepare_lin_intp(z_me, p_z3d_edge, nblks_e, npromz_e, nlev, nzlev,                 & !in
      &                   vcoeff_z%lin_edge%wfac_lin, vcoeff_z%lin_edge%idx0_lin,           & !out
      &                   vcoeff_z%lin_edge%bot_idx_lin  )                                    !out
    CALL prepare_extrap(z_me, nblks_e, npromz_e, nlev,                                      & !in
      &                 vcoeff_z%lin_edge%kpbl1, vcoeff_z%lin_edge%wfacpbl1,                & !out
      &                 vcoeff_z%lin_edge%kpbl2, vcoeff_z%lin_edge%wfacpbl2 )                 !out
    CALL prepare_cubic_intp(z_me, p_z3d_edge, nblks_e, npromz_e, nlev, nzlev,               & !in
      &                     vcoeff_z%cub_edge%coef1, vcoeff_z%cub_edge%coef2,               & !out
      &                     vcoeff_z%cub_edge%coef3,                                        & !out
      &                     vcoeff_z%cub_edge%idx0_cub, vcoeff_z%cub_edge%bot_idx_cub )       !out

  END SUBROUTINE prepare_vert_interp_z


  !-------------
  !>
  !! Compute coefficients for vertical interpolation of model-level
  !! fields to constant-pressure levels.
  !!
  !! @par Revision History
  !! Initial version: see SR prepare_vert_interp_z
  !!
  SUBROUTINE prepare_vert_interp_p(p_patch, p_diag, p_metrics, intp_hrz, nplev,     &
    &                              gh_p_out, temp_p_out, p_p3d_out, vcoeff_p)

    TYPE(t_patch),        TARGET,      INTENT(IN)    :: p_patch
    TYPE(t_nh_diag),      POINTER                    :: p_diag
    TYPE(t_nh_metrics),                INTENT(IN)    :: p_metrics
    TYPE(t_int_state),                 INTENT(IN)    :: intp_hrz          !< pointer to data structure for interpolation
    INTEGER,                           INTENT(IN)    :: nplev             !< number of output levels (pres)
    REAL(wp),                          INTENT(INOUT) :: gh_p_out(:,:,:), &
      &                                                 temp_p_out(:,:,:)
    REAL(wp),             POINTER                    :: p_p3d_out(:,:,:)  !< pressure field
    TYPE(t_vcoeff),                    INTENT(INOUT) :: vcoeff_p          !< out

    ! LOCAL VARIABLES
!!$    CHARACTER(*), PARAMETER :: routine = TRIM("mo_nh_vert_interp:prepare_vert_interp_p")

    INTEGER :: nlev, nlevp1, jg, i_endblk, nblks_c,nblks_e, npromz_c,npromz_e !< blocking parameters
    ! Auxiliary field for smoothing temperature and geopotential on pressure levels
    REAL(wp), DIMENSION(nproma,nplev,p_patch%nblks_c)        :: z_auxp
    REAL(wp), DIMENSION(nproma,nplev,p_patch%nblks_e)        :: gh_p_edge
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_me
    ! Pointer to virtual temperature / temperature, depending on whether the run is moist or dry
    REAL(wp), POINTER, DIMENSION(:,:,:) :: ptr_tempv

    vcoeff_p%l_initialized = .TRUE.
    IF (p_patch%n_patch_cells==0) RETURN

    nlev     = p_patch%nlev
    nlevp1   = p_patch%nlevp1
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c
    nblks_e  = p_patch%nblks_e
    npromz_e = p_patch%npromz_e
    jg       = p_patch%id

    ! allocate coefficient table:
    CALL vcoeff_allocate(nblks_c, nblks_e, nplev, vcoeff_p)

    !--- Coefficients: Interpolation to pressure-level fields
    IF ( ANY((/PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR, PRES_MSL_METHOD_DWD/)== itype_pres_msl) ) THEN
      CALL prepare_extrap_ifspp(p_metrics%z_ifc, p_metrics%z_mc, nblks_c, npromz_c, nlev, & !in
        &                       vcoeff_p%lin_cell%kpbl1, vcoeff_p%lin_cell%zextrap,       & !out
        &                       vcoeff_p%lin_cell%wfacpbl1)                                 !out
      vcoeff_p%lin_cell%kpbl2(:,:) = 0
      vcoeff_p%lin_cell%wfacpbl2(:,:) = 0._wp
    ELSE
      CALL prepare_extrap(p_metrics%z_mc, nblks_c, npromz_c, nlev,                          & !in
        &                 vcoeff_p%lin_cell%kpbl1, vcoeff_p%lin_cell%wfacpbl1,              & !out
        &                 vcoeff_p%lin_cell%kpbl2, vcoeff_p%lin_cell%wfacpbl2 )               !out
    ENDIF

    ! Compute height at pressure levels (i.e. geopot/g); this height
    ! field is afterwards also used as target coordinate for vertical
    ! interpolation
    IF (  iforcing == inwp .OR. iforcing == iecham  ) THEN
      ptr_tempv => p_diag%tempv
    ELSE
      ptr_tempv => p_diag%temp
    ENDIF
    CALL z_at_plevels(p_diag%pres, ptr_tempv, p_metrics%z_mc,                               & !in
      &               p_p3d_out, z_auxp, nblks_c, npromz_c, nlev, nplev,                    & !in,out,in
      &               vcoeff_p%lin_cell%kpbl1, vcoeff_p%lin_cell%wfacpbl1,                  & !in
      &               vcoeff_p%lin_cell%kpbl2, vcoeff_p%lin_cell%wfacpbl2,                  & !in
      &               vcoeff_p%lin_cell%zextrap)                                              !in

    IF (jg > 1) THEN ! copy outermost nest boundary row in order to avoid missing values
      i_endblk = p_patch%cells%end_blk(1,1)
      gh_p_out(:,:,1:i_endblk) = z_auxp(:,:,1:i_endblk)
    ENDIF
    CALL cell_avg(z_auxp, p_patch, p_int_state(jg)%c_bln_avg, gh_p_out)

    ! Prepare again interpolation coefficients (now for pressure levels)
    ! Note: the coefficients are partly identical to the z-level
    !       setup, but they are computed independently to avoid some
    !       complicated code.
    CALL prepare_lin_intp(p_metrics%z_mc, gh_p_out, nblks_c, npromz_c, nlev, nplev,         & !in
      &                   vcoeff_p%lin_cell%wfac_lin, vcoeff_p%lin_cell%idx0_lin,           & !out
      &                   vcoeff_p%lin_cell%bot_idx_lin )                                     !out
    CALL prepare_cubic_intp(p_metrics%z_mc, gh_p_out, nblks_c, npromz_c, nlev, nplev,       & !in
      &                     vcoeff_p%cub_cell%coef1, vcoeff_p%cub_cell%coef2,               & !out
      &                     vcoeff_p%cub_cell%coef3,                                        & !out
      &                     vcoeff_p%cub_cell%idx0_cub, vcoeff_p%cub_cell%bot_idx_cub )       !out

    ! Perform vertical interpolation
    CALL temperature_intp(p_diag%temp, z_auxp, p_metrics%z_mc, gh_p_out,                    & !in,out
      &                   nblks_c, npromz_c, nlev, nplev,                                   & !in
      &                   vcoeff_p%cub_cell%coef1, vcoeff_p%cub_cell%coef2,                 & !in
      &                   vcoeff_p%cub_cell%coef3,                                          & !in
      &                   vcoeff_p%lin_cell%wfac_lin, vcoeff_p%cub_cell%idx0_cub,           & !in
      &                   vcoeff_p%lin_cell%idx0_lin, vcoeff_p%cub_cell%bot_idx_cub,        & !in
      &                   vcoeff_p%lin_cell%bot_idx_lin, vcoeff_p%lin_cell%wfacpbl1,        & !in
      &                   vcoeff_p%lin_cell%kpbl1, vcoeff_p%lin_cell%wfacpbl2,              & !in
      &                   vcoeff_p%lin_cell%kpbl2, l_restore_sfcinv=.FALSE.,                & !in
      &                   l_hires_corr=.FALSE., extrapol_dist=0._wp, l_pz_mode=.TRUE.,      & !in
      &                   zextrap=vcoeff_p%lin_cell%zextrap)

    IF (jg > 1) THEN ! copy outermost nest boundary row in order to avoid missing values
      i_endblk = p_patch%cells%end_blk(1,1)
      temp_p_out(:,:,1:i_endblk) = z_auxp(:,:,1:i_endblk)
    ENDIF
    CALL cell_avg(z_auxp, p_patch, p_int_state(jg)%c_bln_avg, temp_p_out)

    IF ( ANY((/PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR, PRES_MSL_METHOD_DWD/)== itype_pres_msl) ) THEN
      CALL prepare_extrap(p_metrics%z_mc, nblks_c, npromz_c, nlev,                          & !in
        &                 vcoeff_p%lin_cell%kpbl1, vcoeff_p%lin_cell%wfacpbl1,              & !out
        &                 vcoeff_p%lin_cell%kpbl2, vcoeff_p%lin_cell%wfacpbl2 )               !out
    ENDIF

    !--- Interpolation data for the vertical interface of cells, "nlevp1"

    CALL prepare_lin_intp(p_metrics%z_ifc, gh_p_out, nblks_c, npromz_c, nlevp1, nplev,      & !in
      &                   vcoeff_p%lin_cell_nlevp1%wfac_lin,                                & !out
      &                   vcoeff_p%lin_cell_nlevp1%idx0_lin,                                & !out
      &                   vcoeff_p%lin_cell_nlevp1%bot_idx_lin  )                             !out
    CALL prepare_extrap(p_metrics%z_ifc, nblks_c, npromz_c, nlevp1,                         & !in
      &                 vcoeff_p%lin_cell_nlevp1%kpbl1, vcoeff_p%lin_cell_nlevp1%wfacpbl1,  & !out
      &                 vcoeff_p%lin_cell_nlevp1%kpbl2, vcoeff_p%lin_cell_nlevp1%wfacpbl2 )   !out

    !--- Compute data for interpolation of edge-based fields

    CALL cells2edges_scalar(p_metrics%z_mc, p_patch, intp_hrz%c_lin_e, z_me, opt_fill_latbc=.TRUE.)
    CALL cells2edges_scalar(gh_p_out, p_patch, intp_hrz%c_lin_e, gh_p_edge, opt_fill_latbc=.TRUE.)
    CALL sync_patch_array_mult(SYNC_E,p_patch,2,z_me,gh_p_edge)

    CALL prepare_lin_intp(z_me, gh_p_edge, nblks_e, npromz_e, nlev, nplev,                   & !in
      &                   vcoeff_p%lin_edge%wfac_lin, vcoeff_p%lin_edge%idx0_lin,            & !out
      &                   vcoeff_p%lin_edge%bot_idx_lin )                                      !out
    CALL prepare_extrap(z_me, nblks_e, npromz_e, nlev,                                       & !in
      &                 vcoeff_p%lin_edge%kpbl1, vcoeff_p%lin_edge%wfacpbl1,                 & !out
      &                 vcoeff_p%lin_edge%kpbl2, vcoeff_p%lin_edge%wfacpbl2 )                  !out
    CALL prepare_cubic_intp(z_me, gh_p_edge, nblks_e, npromz_e, nlev, nplev,                 & !in
      &                     vcoeff_p%cub_edge%coef1, vcoeff_p%cub_edge%coef2,                & !out
      &                     vcoeff_p%cub_edge%coef3,                                         & !out
      &                     vcoeff_p%cub_edge%idx0_cub, vcoeff_p%cub_edge%bot_idx_cub )        !out

  END SUBROUTINE prepare_vert_interp_p


  !-------------
  !>
  !! Compute coefficients for vertical interpolation of model-level
  !! fields to isentropes.
  !!
  !! @par Revision History
  !! Initial version: see SR prepare_vert_interp_i
  !!
  SUBROUTINE prepare_vert_interp_i(p_patch, p_prog, p_diag, p_metrics, intp_hrz, nilev,     &
    &                              gh_i_out, temp_i_out, p_i3d_out, vcoeff_i)

    TYPE(t_patch),        TARGET,      INTENT(IN)    :: p_patch
    TYPE(t_nh_prog),      POINTER                    :: p_prog
    TYPE(t_nh_diag),      POINTER                    :: p_diag
    TYPE(t_nh_metrics),                INTENT(IN)    :: p_metrics
    TYPE(t_int_state),                 INTENT(IN)    :: intp_hrz          !< pointer to data structure for interpolation
    INTEGER,                           INTENT(IN)    :: nilev             !< number of output levels (isentropes)
    REAL(wp),                          INTENT(INOUT) :: gh_i_out(:,:,:), &
      &                                                 temp_i_out(:,:,:)
    REAL(wp),             POINTER                    :: p_i3d_out(:,:,:)  !< theta field
    TYPE(t_vcoeff),                    INTENT(INOUT) :: vcoeff_i          !< out

    ! LOCAL VARIABLES
!!$    CHARACTER(*), PARAMETER :: routine = TRIM("mo_nh_vert_interp:prepare_vert_interp_i")
    INTEGER :: nlev, nlevp1, nblks_c, nblks_e, npromz_c,npromz_e !< blocking parameters
    REAL(wp), DIMENSION(nproma,nilev,p_patch%nblks_e)        :: gh_i_edge
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_me

    vcoeff_i%l_initialized = .TRUE.

    nlev     = p_patch%nlev
    nlevp1   = p_patch%nlevp1
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c
    nblks_e  = p_patch%nblks_e
    npromz_e = p_patch%npromz_e

    ! allocate coefficient table:
    CALL vcoeff_allocate(nblks_c, nblks_e, nilev, vcoeff_i)

    !--- Coefficients: Interpolation to isentropes

    ! Compute height at isentropic levels (i.e. geopot/g); this height
    ! field is afterwards also used as target coordinate for vertical
    ! interpolation

    ! Note: the following coefficients are identical to the z-level
    !       setup, but they are computed independently to avoid some
    !       complicated code:
    CALL prepare_extrap(p_metrics%z_mc, nblks_c, npromz_c, nlev,                            & !in
      &                 vcoeff_i%lin_cell%kpbl1, vcoeff_i%lin_cell%wfacpbl1,                & !out
      &                 vcoeff_i%lin_cell%kpbl2, vcoeff_i%lin_cell%wfacpbl2 )                 !out

    CALL z_at_theta_levels(p_prog%theta_v, p_metrics%z_mc,                                  & !in
      &                    p_i3d_out, gh_i_out,                                             & !out
      &                    vcoeff_i%lin_cell%wfacpbl1, vcoeff_i%lin_cell%kpbl1,             & !in
      &                    vcoeff_i%lin_cell%wfacpbl2, vcoeff_i%lin_cell%kpbl2, nblks_c,    & !in
      &                    npromz_c, nlev, nilev)

    ! Prepare again interpolation coefficients (now for isentropic levels)
    CALL prepare_lin_intp(p_metrics%z_mc, gh_i_out, nblks_c, npromz_c, nlev, nilev,         & !in
      &                   vcoeff_i%lin_cell%wfac_lin, vcoeff_i%lin_cell%idx0_lin,           & !out
      &                   vcoeff_i%lin_cell%bot_idx_lin )                                     !out

    CALL prepare_cubic_intp(p_metrics%z_mc, gh_i_out, nblks_c, npromz_c, nlev, nilev,       & !in
      &                     vcoeff_i%cub_cell%coef1, vcoeff_i%cub_cell%coef2,               & !out
      &                     vcoeff_i%cub_cell%coef3,                                        & !out
      &                     vcoeff_i%cub_cell%idx0_cub, vcoeff_i%cub_cell%bot_idx_cub )       !out

    ! Perform vertical interpolation
    CALL temperature_intp(p_diag%temp, temp_i_out, p_metrics%z_mc, gh_i_out,                & !in,out
      &                   nblks_c, npromz_c, nlev, nilev,                                   & !in
      &                   vcoeff_i%cub_cell%coef1, vcoeff_i%cub_cell%coef2,                 & !in
      &                   vcoeff_i%cub_cell%coef3,                                          & !in
      &                   vcoeff_i%lin_cell%wfac_lin, vcoeff_i%cub_cell%idx0_cub,           & !in
      &                   vcoeff_i%lin_cell%idx0_lin, vcoeff_i%cub_cell%bot_idx_cub,        & !in
      &                   vcoeff_i%lin_cell%bot_idx_lin, vcoeff_i%lin_cell%wfacpbl1,        & !in
      &                   vcoeff_i%lin_cell%kpbl1, vcoeff_i%lin_cell%wfacpbl2,              & !in
      &                   vcoeff_i%lin_cell%kpbl2, l_restore_sfcinv=.TRUE.,                 & !in
      &                   l_hires_corr=.FALSE., extrapol_dist=0._wp, l_pz_mode=.TRUE.)        !in

    !--- Interpolation data for the vertical interface of cells, "nlevp1"

    CALL prepare_lin_intp(p_metrics%z_ifc, gh_i_out,                                        & !in
      &                   nblks_c, npromz_c, nlevp1, nilev,                                 & !in
      &                   vcoeff_i%lin_cell_nlevp1%wfac_lin,                                & !out
      &                   vcoeff_i%lin_cell_nlevp1%idx0_lin,                                & !out
      &                   vcoeff_i%lin_cell_nlevp1%bot_idx_lin  )                             !out
    CALL prepare_extrap(p_metrics%z_ifc, nblks_c, npromz_c, nlevp1,                         & !in
      &                 vcoeff_i%lin_cell_nlevp1%kpbl1, vcoeff_i%lin_cell_nlevp1%wfacpbl1,  & !out
      &                 vcoeff_i%lin_cell_nlevp1%kpbl2, vcoeff_i%lin_cell_nlevp1%wfacpbl2 )   !out

    !--- Compute data for interpolation of edge-based fields

    CALL cells2edges_scalar(p_metrics%z_mc, p_patch, intp_hrz%c_lin_e, z_me, opt_fill_latbc=.TRUE.)
    CALL cells2edges_scalar(gh_i_out, p_patch, intp_hrz%c_lin_e, gh_i_edge, opt_fill_latbc=.TRUE.)
    CALL sync_patch_array_mult(SYNC_E,p_patch,2,z_me,gh_i_edge)

    CALL prepare_lin_intp(z_me, gh_i_edge, nblks_e, npromz_e, nlev, nilev,                   & !in
      &                   vcoeff_i%lin_edge%wfac_lin, vcoeff_i%lin_edge%idx0_lin,            & !out
      &                   vcoeff_i%lin_edge%bot_idx_lin )                                      !out
    CALL prepare_extrap(z_me, nblks_e, npromz_e, nlev,                                       & !in
      &                 vcoeff_i%lin_edge%kpbl1, vcoeff_i%lin_edge%wfacpbl1,                 & !out
      &                 vcoeff_i%lin_edge%kpbl2, vcoeff_i%lin_edge%wfacpbl2 )                  !out
    CALL prepare_cubic_intp(z_me, gh_i_edge, nblks_e, npromz_e, nlev, nilev,                 & !in
      &                     vcoeff_i%cub_edge%coef1, vcoeff_i%cub_edge%coef2,                & !out
      &                     vcoeff_i%cub_edge%coef3,                                         & !out
      &                     vcoeff_i%cub_edge%idx0_cub, vcoeff_i%cub_edge%bot_idx_cub )        !out

  END SUBROUTINE prepare_vert_interp_i



  !-------------
  !>
  !! SUBROUTINE prepare_lin_intp
  !! Computes coefficient fields for linear vertical interpolation
  !!
  !! Required input fields: 3D coordinate fields of input and output data
  !! Output: weighting coefficient, index of upper level,
  !! index of lowest level for which interpolation is possible
  !! (as opposed to extrapolation)
  !!
  !! It is assumed that the highest level of the input data is at least
  !! as high as the highest level of the output data
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE prepare_lin_intp(z3d_in, z3d_out,                    &
                              nblks, npromz, nlevs_in, nlevs_out, &
                              wfac, idx0, bot_idx, lextrap        )

    ! Input fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:) ! height coordinate field of input data (m)
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)! height coordinate field of input data (m)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Switch for extrapolation beyond top of input data
    LOGICAL, INTENT(IN), OPTIONAL :: lextrap ! if true, apply linear extrapolation

    ! Output fields
    REAL(wp), INTENT(OUT) :: wfac(:,:,:)       ! weighting factor of upper level
    INTEGER , INTENT(OUT) :: idx0(:,:,:)       ! index of upper level
    INTEGER , INTENT(OUT) :: bot_idx(:,:)      ! index of lowest level for which interpolation is possible

    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jk1, jk_start
    INTEGER :: nlen, ierror(nblks), nerror
    LOGICAL :: l_found(nproma), lfound_all, l_extrap

!-------------------------------------------------------------------------

    IF (PRESENT(lextrap)) THEN
      l_extrap = lextrap
    ELSE
      l_extrap = .TRUE.
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk1,jk_start,l_found,lfound_all) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        bot_idx(nlen+1:nproma,jb) = nlevs_out
        wfac(nlen+1:nproma,:,jb)  = 0.5_wp
        idx0(nlen+1:nproma,:,jb)  = nlevs_in-1
      ENDIF
      ierror(jb) = 0

      jk_start = 1
      DO jk = 1, nlevs_out
        l_found(:) = .FALSE.
        lfound_all = .FALSE.
        DO jk1 = jk_start,nlevs_in-1
          DO jc = 1, nlen
            IF (z3d_out(jc,jk,jb) <= z3d_in(jc,jk1,jb) .AND. &
                z3d_out(jc,jk,jb) >  z3d_in(jc,jk1+1,jb)) THEN
              idx0(jc,jk,jb) = jk1
              wfac(jc,jk,jb) = (z3d_out(jc,jk,jb)-z3d_in(jc,jk1+1,jb))/&
                               (z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+1,jb))
              bot_idx(jc,jb) = jk
              l_found(jc) = .TRUE.
            ELSE IF (z3d_out(jc,jk,jb) <= z3d_in(jc,nlevs_in,jb)) THEN
              l_found(jc) = .TRUE.
              idx0(jc,jk,jb) = nlevs_in
              IF (jk == 1)  bot_idx(jc,jb) = 0
            ELSE IF (z3d_out(jc,jk,jb) > z3d_in(jc,1,jb)) THEN ! linear extrapolation
              idx0(jc,jk,jb) = 1
              IF (l_extrap) THEN
                wfac(jc,jk,jb) = (z3d_out(jc,jk,jb)-z3d_in(jc,2,jb))/&
                                 (z3d_in(jc,1,jb)-z3d_in(jc,2,jb))
              ELSE
                wfac(jc,jk,jb) = 1._wp ! use constant values above top of input data
              ENDIF
              bot_idx(jc,jb) = jk
              l_found(jc) = .TRUE.
            ENDIF
          ENDDO
          IF (ALL(l_found(1:nlen))) THEN
            lfound_all = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (lfound_all) THEN
          jk_start = MIN(MINVAL(idx0(1:nlen,jk,jb)),nlevs_in-1)
        ELSE
          ierror(jb) = ierror(jb) + 1
          ! Write extra debug output before finishing
          WRITE(0,*) 'prepare_lin_intp',jb,jk,nlevs_out,jk_start,jk1,nlevs_in
          DO jc = 1, nlen
            IF(.NOT.l_found(jc)) THEN
              WRITE(0,*)'prepare_lin_intp',z3d_in(jc,jk_start:jk1,jb),&
                z3d_in(jc,nlevs_in,jb),z3d_out(jc,jk,jb)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DO jk = MINVAL(bot_idx(1:nlen,jb))+1, nlevs_out
        DO jc = 1, nlen
          IF (jk >= bot_idx(jc,jb)+1) THEN
            ! Store extrapolation distance on wfac if target point is below the
            ! surface of the input data (note: this is a negative quantity)
            idx0(jc,jk,jb) = nlevs_in
            wfac(jc,jk,jb) = z3d_out(jc,jk,jb)-z3d_in(jc,nlevs_in,jb)
          ENDIF
        ENDDO
      ENDDO


    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    nerror = SUM(ierror)

    IF (nerror > 0) CALL finish("prepare_lin_intp:",&
      "Error in computing linear interpolation coefficients")

  END SUBROUTINE prepare_lin_intp

  !-------------
  !>
  !! SUBROUTINE prepare_extrap
  !! Computes coefficient fields for vertical extrapolation below the surface level
  !! of the input model
  !!
  !! Required input fields: 3D coordinate fields of input data
  !! Output: index and coefficient fields to compute field values at
  !! 500 m and 1000 m above ground
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE prepare_extrap(z3d_in, nblks, npromz, nlevs_in, &
                            kpbl1, wfacpbl1, kpbl2, wfacpbl2 )

    ! Input fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:) ! height coordinate field of input data (m)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! Output fields
    INTEGER , INTENT(OUT) :: kpbl1(:,:), & ! Indices of model levels lying immediately above
                             kpbl2(:,:)    ! (by default) 500 m / 1000 m AGL
    REAL(wp), INTENT(OUT) :: wfacpbl1(:,:), & ! Corresponding interpolation coefficients
                             wfacpbl2(:,:)

    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jk_start
    INTEGER :: nlen


!-------------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk_start) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        kpbl1(nlen+1:nproma,jb) = nlevs_in
        kpbl2(nlen+1:nproma,jb) = nlevs_in

        wfacpbl1(nlen+1:nproma,jb) = 0.5_wp
        wfacpbl2(nlen+1:nproma,jb) = 0.5_wp
      ENDIF

      DO jk = 1, nlevs_in
        IF (MINVAL(z3d_in(1:nlen,jk,jb)-z3d_in(1:nlen,nlevs_in,jb)) <= zpbl2) THEN
          jk_start = jk - 1
          EXIT
        ENDIF
      ENDDO

      DO jk = jk_start, nlevs_in-1
        DO jc = 1, nlen

          IF (z3d_in(jc,jk,jb)  >= z3d_in(jc,nlevs_in,jb)+zpbl2 .AND. &
              z3d_in(jc,jk+1,jb) < z3d_in(jc,nlevs_in,jb)+zpbl2) THEN
            kpbl2(jc,jb) = jk
            wfacpbl2(jc,jb) = (z3d_in(jc,nlevs_in,jb)+zpbl2 - z3d_in(jc,jk+1,jb)) / &
                              (z3d_in(jc,jk,jb)            - z3d_in(jc,jk+1,jb))
          ENDIF

          IF (z3d_in(jc,jk,jb)  >= z3d_in(jc,nlevs_in,jb)+zpbl1 .AND. &
              z3d_in(jc,jk+1,jb) < z3d_in(jc,nlevs_in,jb)+zpbl1) THEN
            kpbl1(jc,jb) = jk
            wfacpbl1(jc,jb) = (z3d_in(jc,nlevs_in,jb)+zpbl1 - z3d_in(jc,jk+1,jb)) / &
                              (z3d_in(jc,jk,jb)            - z3d_in(jc,jk+1,jb))
          ENDIF

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE prepare_extrap


  !-------------
  !>
  !! SUBROUTINE prepare_extrap_ifspp
  !! Computes coefficient fields for vertical extrapolation below the surface level
  !! of the input model for IFS postprocessing method
  !!
  !! Required input fields: 3D coordinate fields of input data
  !! Output: index and coefficient fields to compute field values at
  !! 150 m above ground
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2013-04-11)
  !!
  !!
  !!
  SUBROUTINE prepare_extrap_ifspp(z3d_h_in, z3d_in, nblks, npromz, nlevs_in, kextrap, zextrap, wfac_extrap)

    ! Input fields
    REAL(wp), INTENT(IN) :: z3d_h_in(:,:,:) ! half-level height coordinate field of input data (m)
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:)   ! full-level height coordinate field of input data (m)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! Output fields
    INTEGER , INTENT(OUT) :: kextrap(:,:) ! Indices of model levels lying immediately above zextrap
    REAL(wp), INTENT(OUT) :: zextrap(:,:) ! AGL height from which downward extrapolation starts (between 10 m and 150 m)
    REAL(wp), INTENT(OUT) :: wfac_extrap(:,:) ! Corresponding interpolation coefficients

    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jk_start
    INTEGER :: nlen
    REAL(wp) :: zagl_extrap


!-------------------------------------------------------------------------

    ! Use extrapolation from 150 m AGL if one of the IFS methods is selected, and 
    ! orography-height dependent blending between 10 m and 150 m AGL if the new DWD
    ! method is selected (which constitutes a mixture between the IFS method and the old GME method)
    SELECT CASE (itype_pres_msl)
    CASE (PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR)
      zagl_extrap = zagl_extrap_2
    CASE (PRES_MSL_METHOD_DWD)
      zagl_extrap = zagl_extrap_1
    END SELECT

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk_start) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      kextrap(:,jb) = -1
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        kextrap(nlen+1:nproma,jb) = nlevs_in
        wfac_extrap(nlen+1:nproma,jb) = 0.5_wp
      ENDIF

      ! Compute start height above ground for downward extrapolation, depending on topography height
      DO jc = 1, nlen
        IF (z3d_h_in(jc,nlevs_in+1,jb) <= topo_extrap_1) THEN
          zextrap(jc,jb) = zagl_extrap
        ELSE IF (z3d_h_in(jc,nlevs_in+1,jb) >= topo_extrap_1 + trans_depth) THEN
          zextrap(jc,jb) = zagl_extrap_2
        ELSE
          zextrap(jc,jb) = zagl_extrap + (zagl_extrap_2 - zagl_extrap) *            &
                           (z3d_h_in(jc,nlevs_in+1,jb) - topo_extrap_1) / trans_depth
        ENDIF
      ENDDO

      jk_start = nlevs_in-1
      DO jk = 1, nlevs_in
        IF (MINVAL(z3d_in(1:nlen,jk,jb)-z3d_h_in(1:nlen,nlevs_in+1,jb)) <= zagl_extrap_2) THEN
          jk_start = jk - 1
          EXIT
        ENDIF
      ENDDO

      DO jk = jk_start, nlevs_in-1
        DO jc = 1, nlen

          IF (z3d_in(jc,jk,jb)  >= z3d_h_in(jc,nlevs_in+1,jb)+zextrap(jc,jb) .AND. &
              z3d_in(jc,jk+1,jb) < z3d_h_in(jc,nlevs_in+1,jb)+zextrap(jc,jb) .OR.  & 
              kextrap(jc,jb) == -1 .AND. jk == nlevs_in-1) THEN
            kextrap(jc,jb) = jk
            wfac_extrap(jc,jb) = (z3d_h_in(jc,nlevs_in+1,jb)+zextrap(jc,jb) - z3d_in(jc,jk+1,jb)) / &
                                 (z3d_in(jc,jk,jb)                          - z3d_in(jc,jk+1,jb))
          ENDIF

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE prepare_extrap_ifspp

  !>
  !! SUBROUTINE prepare_cubic_intp
  !! Computes coefficient fields for cubic vertical interpolation
  !!
  !! Required input fields: 3D coordinate fields of input and output data
  !! Output: weighting coefficients, index of upper level,
  !! index of lowest level for which cubic interpolation is possible
  !!
  !! It is assumed that the highest level of the input data is at least
  !! as high as the highest level of the output data
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE prepare_cubic_intp(z3d_in, z3d_out,                    &
                                nblks, npromz, nlevs_in, nlevs_out, &
                                coef1, coef2, coef3, idx0, bot_idx)

    ! Input fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:) ! height coordinate field of input data (m)
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)! height coordinate field of input data (m)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Output fields
    REAL(wp), INTENT(OUT) :: coef1(:,:,:)      ! coefficient for linear term
    REAL(wp), INTENT(OUT) :: coef2(:,:,:)      ! coefficient for quadratic term
    REAL(wp), INTENT(OUT) :: coef3(:,:,:)      ! coefficient for cubic term
    INTEGER , INTENT(OUT) :: idx0(:,:,:)       ! index of upper level
    INTEGER , INTENT(OUT) :: bot_idx(:,:)      ! index of lowest level for which interpolation is possible

    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jk1, jk_start
    INTEGER :: nlen, ierror(nblks), nerror
    LOGICAL :: l_found(nproma),lfound_all

!-------------------------------------------------------------------------

    IF (nlevs_in <= 3) CALL finish("prepare_cubic_intp:",&
      "Error, number of levels too small for cubic interpolation")

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk1,jk_start,l_found,lfound_all) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        bot_idx(nlen+1:nproma,jb) = nlevs_out
        coef1(nlen+1:nproma,:,jb) = 0.5_wp
        coef2(nlen+1:nproma,:,jb) = 0.5_wp
        coef3(nlen+1:nproma,:,jb) = 0.5_wp
        idx0(nlen+1:nproma,:,jb)  = nlevs_in-1
      ENDIF
      ierror(jb) = 0

      jk_start = 2
      DO jk = 1, nlevs_out
        l_found(:) = .FALSE.
        lfound_all = .FALSE.
        DO jk1 = jk_start,nlevs_in-2 ! cubic interpolation requires jk1-1, jk1, jk1+1 and jk1+2
          DO jc = 1, nlen
            IF (z3d_out(jc,jk,jb) <= z3d_in(jc,jk1,jb) .AND. &
                z3d_out(jc,jk,jb) >  z3d_in(jc,jk1+1,jb)) THEN
              idx0(jc,jk,jb)  = jk1-1
              coef1(jc,jk,jb) = z3d_out(jc,jk,jb)-z3d_in(jc,jk1-1,jb)
              coef2(jc,jk,jb) = (z3d_out(jc,jk,jb)-z3d_in(jc,jk1-1,jb))*&
                                (z3d_out(jc,jk,jb)-z3d_in(jc,jk1  ,jb))
              coef3(jc,jk,jb) = (z3d_out(jc,jk,jb)-z3d_in(jc,jk1-1,jb))*&
                                (z3d_out(jc,jk,jb)-z3d_in(jc,jk1  ,jb))*&
                                (z3d_out(jc,jk,jb)-z3d_in(jc,jk1+1,jb))
              bot_idx(jc,jb) = jk
              l_found(jc) = .TRUE.
            ELSE IF (z3d_out(jc,jk,jb) <= z3d_in(jc,nlevs_in-1,jb)) THEN
              l_found(jc) = .TRUE.
              idx0(jc,jk,jb) = nlevs_in-1
            ELSE IF (z3d_out(jc,jk,jb) >  z3d_in(jc,2,jb)) THEN
              ! linear interpolation between two upper input levels or extrapolation beyond top level
              idx0(jc,jk,jb) = 1
              coef1(jc,jk,jb) = z3d_out(jc,jk,jb)-z3d_in(jc,1,jb)
              coef2(jc,jk,jb) = 0._wp
              coef3(jc,jk,jb) = 0._wp
              bot_idx(jc,jb) = jk
              l_found(jc) = .TRUE.
            ENDIF
          ENDDO
          IF (ALL(l_found(1:nlen))) THEN
            lfound_all = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (lfound_all) THEN
          jk_start = MIN(MINVAL(idx0(1:nlen,jk,jb))+1,nlevs_in-2)
        ELSE
          ierror(jb) = ierror(jb) + 1
        ENDIF
      ENDDO

      DO jk = MINVAL(bot_idx(1:nlen,jb))+1, nlevs_out
        DO jc = 1, nlen
          IF (jk >= bot_idx(jc,jb)+1) THEN
            coef1(jc,jk,jb) = 0._wp
            coef2(jc,jk,jb) = 0._wp
            coef3(jc,jk,jb) = 0._wp
            idx0 (jc,jk,jb) = nlevs_in
           ENDIF
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    nerror = SUM(ierror)

    IF (nerror > 0) CALL finish("prepare_cubic_intp:",&
      "Error in computing cubic interpolation coefficients")

  END SUBROUTINE prepare_cubic_intp


  !-------------
  !>
  !! SUBROUTINE lin_intp
  !! Performs linear vertical interpolation of a 3D field
  !!
  !! Required input fields: 3D input field to be interpolated,
  !! coefficient fields from prepare_lin_intp and prepare_extrap
  !! Output: interpolated 3D field
  !!
  !! Note: field values below the surface level of the input field are
  !! computed by linear extrapolation, using the gradient of the lowest (by default)
  !! 500 m above ground. Do not use for temperature!
  !!
  !! Setting l_loglin=.TRUE. activates logarithmic interpolation
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE lin_intp(f3d_in, f3d_out,                      &
                      nblks, npromz, nlevs_in, nlevs_out,   &
                      wfac, idx0, bot_idx, wfacpbl1, kpbl1, &
                      wfacpbl2, kpbl2, l_loglin, l_pd_limit,&
                      l_extrapol, lower_limit               )


    ! Atmospheric fields
    REAL(wp), INTENT(IN)  :: f3d_in (:,:,:) ! input field
    REAL(wp), INTENT(OUT) :: f3d_out(:,:,:) ! output field (on ICON vertical grid)

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: wfac(:,:,:)    ! weighting factor of upper level
    INTEGER , INTENT(IN) :: idx0(:,:,:)    ! index of upper level
    INTEGER , INTENT(IN) :: bot_idx(:,:)   ! index of lowest level for which interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)     ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)  ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)     ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)  ! corresponding interpolation coefficient

    ! Control switches
    LOGICAL,  INTENT(IN) :: l_loglin    ! switch for logarithmic interpolation
    LOGICAL,  INTENT(IN) :: l_pd_limit  ! switch for use of positive definite limiter
    LOGICAL,  INTENT(IN) :: l_extrapol  ! switch for use of downward extrapolation (no-gradient condition otherwise)

    REAL(wp), INTENT(IN), OPTIONAL :: lower_limit ! lower limit of variable

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jc
    INTEGER  :: nlen
    REAL(wp) :: zf_in_tr(nproma,nlevs_in), zf_in_lim(nproma,nlevs_in), z_limit, f3d_z1, f3d_z2, vgrad_f3d

!-------------------------------------------------------------------------

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0))) RETURN

    IF (PRESENT(lower_limit)) THEN
      z_limit = lower_limit
    ELSE
      z_limit = 0._wp
    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,zf_in_tr,zf_in_lim,f3d_z1,f3d_z2,vgrad_f3d) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        f3d_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

      IF (l_pd_limit) THEN
        zf_in_lim(1:nlen,1:nlevs_in) = MAX(z_limit,f3d_in(1:nlen,1:nlevs_in,jb))
      ELSE
        zf_in_lim(1:nlen,1:nlevs_in) = f3d_in(1:nlen,1:nlevs_in,jb)
      ENDIF

      IF (l_loglin) THEN
        zf_in_tr(1:nlen,1:nlevs_in) = LOG(MAX(1.e-20_wp,zf_in_lim(1:nlen,1:nlevs_in)))
      ELSE
        zf_in_tr(1:nlen,1:nlevs_in) = zf_in_lim(1:nlen,1:nlevs_in)
      ENDIF

      DO jk = 1, nlevs_out
        DO jc = 1, nlen
          IF (jk <= bot_idx(jc,jb)) THEN

            ! linear interpolation
            f3d_out(jc,jk,jb) = wfac(jc,jk,jb)*zf_in_tr(jc,idx0(jc,jk,jb)) + &
              (1._wp-wfac(jc,jk,jb))*zf_in_tr(jc,idx0(jc,jk,jb)+1)

            IF (l_loglin) f3d_out(jc,jk,jb) = EXP(f3d_out(jc,jk,jb))

          ELSE IF (l_extrapol) THEN

            ! linear extrapolation, using the gradient between (by default) 1000 m and 500 m AGL
            ! Logarithmic computation is not used here because it would be numerically unstable for extrapolation

            ! Field value at height zpbl1
            f3d_z1 = wfacpbl1(jc,jb) *zf_in_lim(jc,kpbl1(jc,jb)  ) +  &
              (1._wp-wfacpbl1(jc,jb))*zf_in_lim(jc,kpbl1(jc,jb)+1)

            ! Field value at height zpbl2
            f3d_z2 = wfacpbl2(jc,jb) *zf_in_lim(jc,kpbl2(jc,jb)  ) +  &
              (1._wp-wfacpbl2(jc,jb))*zf_in_lim(jc,kpbl2(jc,jb)+1)

            ! vertical gradient
            vgrad_f3d = (f3d_z2-f3d_z1)/(zpbl2-zpbl1)

            ! wfac carries the (negative) extrapolation distance (in m) in this case
            f3d_out(jc,jk,jb) = zf_in_lim(jc,nlevs_in) + vgrad_f3d*wfac(jc,jk,jb)

          ELSE ! use no-gradient condition for extrapolation

            f3d_out(jc,jk,jb) = zf_in_lim(jc,nlevs_in)

          ENDIF
        ENDDO
      ENDDO

      IF (l_pd_limit) THEN
        f3d_out(1:nlen,1:nlevs_out,jb) = MAX(z_limit,f3d_out(1:nlen,1:nlevs_out,jb))
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE lin_intp


  !-------------
  !>
  !! SUBROUTINE pressure_intp
  !! Performs vertical interpolation of pressure
  !!
  !! Required input fields: pressure, virtual temperature and 3D height coordinate
  !! field of input data, virtual temperature and 3D height coordinate of output data
  !! Output: pressure field of output data
  !!
  !! Method: piecewise analytical integration of the hydrostatic equation
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE pressure_intp(pres_in, tempv_in, z3d_in, pres_out, z3d_out,                 &
                           nblks, npromz, nlevs_in, nlevs_out,                           &
                           wfac, idx0, bot_idx, wfacpbl1, kpbl1, wfacpbl2, kpbl2, zextrap)


    ! Input fields
    REAL(wp), INTENT(IN)  :: pres_in  (:,:,:) ! pressure field of input data
    REAL(wp), INTENT(IN)  :: tempv_in (:,:,:) ! virtual temperature of input data
    REAL(wp), INTENT(IN)  :: z3d_in   (:,:,:) ! 3D height coordinate field of input data
    REAL(wp), INTENT(IN)  :: z3d_out  (:,:,:) ! 3D height coordinate field of output data

    ! Output
    REAL(wp), INTENT(OUT) :: pres_out (:,:,:) ! pressure field of output data

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: wfac(:,:,:)    ! weighting factor of upper level
    INTEGER , INTENT(IN) :: idx0(:,:,:)    ! index of upper level
    INTEGER , INTENT(IN) :: bot_idx(:,:)   ! index of lowest level for which interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)     ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)  ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)     ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)  ! corresponding interpolation coefficient

    REAL(wp), OPTIONAL, INTENT(IN) :: zextrap(:,:)   ! AGL height from which downward extrapolation starts (in postprocesing mode)

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jc, jk1, nlen
    REAL(wp), DIMENSION(nproma,nlevs_out) :: dtvdz_down
    REAL(wp), DIMENSION(nproma)           :: tmsl, tsfc_mod, tempv1, tempv2, vtgrad_up, sfc_inv
    REAL(wp) :: p_up, p_down, t_extr

!-------------------------------------------------------------------------

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0))) RETURN

    IF (.NOT. PRESENT(zextrap) .AND. itype_pres_msl >= 3 ) THEN
      CALL finish("pressure_intp:", "zextrap missing in argument list")
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jk1,jc,dtvdz_down,p_up,p_down,tmsl,tsfc_mod,tempv1,tempv2,&
!$OMP            vtgrad_up,sfc_inv,t_extr) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        pres_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

      IF (itype_pres_msl == PRES_MSL_METHOD_GME) THEN

        ! Preparations for extrapolation below the ground: calculate temperature
        ! profile with artificial limits like in IFS. This temperature profile
        ! is NOT consistent with the temperature on pressure levels computed
        ! afterwards, implying that temperature and geopotential on pressure
        ! levels are not in hydrostatic balance. We are not happy with this
        ! solution, but it is used for the time being in order to be consistent
        ! with IFS and GME.

        DO jc = 1, nlen
          tsfc_mod(jc) = tempv_in(jc,nlevs_in,jb)
          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z3d_in(jc,nlevs_in,jb)
          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF
        ENDDO

      ELSE IF ( itype_pres_msl >= 3 ) THEN

        ! Similar method to option 1, but we use the temperature at 150 m AGL
        ! for extrapolation (kpbl1 contains the required index in this case)

        DO jc = 1, nlen

          tsfc_mod(jc) = wfacpbl1(jc,jb) *tempv_in(jc,kpbl1(jc,jb),jb  ) + &
                  (1._wp-wfacpbl1(jc,jb))*tempv_in(jc,kpbl1(jc,jb)+1,jb) - &
                  dtdz_standardatm*(zextrap(jc,jb)-0.5_wp*vct_a(num_lev(1)))

          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z3d_in(jc,nlevs_in,jb)
          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF
        ENDDO

      ELSE

        ! Use a temperature profile that is (roughly) consistent with
        ! that used for temperature extrapolation

        DO jc = 1, nlen
          ! Virtual temperature at height zpbl1
          tempv1(jc) = wfacpbl1(jc,jb) *tempv_in(jc,kpbl1(jc,jb),jb  ) + &
                (1._wp-wfacpbl1(jc,jb))*tempv_in(jc,kpbl1(jc,jb)+1,jb)

          ! Virtual temperature at height zpbl2
          tempv2(jc) = wfacpbl2(jc,jb) *tempv_in(jc,kpbl2(jc,jb),jb  ) + &
                (1._wp-wfacpbl2(jc,jb))*tempv_in(jc,kpbl2(jc,jb)+1,jb)

          ! Vertical gradient between zpbl1 and zpbl2
          vtgrad_up(jc) = (tempv2(jc) - tempv1(jc))/(zpbl2 - zpbl1)

          ! Modified (extrapolated) surface temperature without inversion
          tsfc_mod(jc) = tempv1(jc) - zpbl1*vtgrad_up(jc)

          ! "surface inversion", defined by the difference between the extrapolated
          ! extrapolated temperature from above and the original input temperature
          sfc_inv(jc) = tsfc_mod(jc) - tempv_in(jc,nlevs_in,jb)

          ! Reduction of the surface inversion depending on the extrapolation
          ! distance. The surface inversion is fully restored for extrapolation distances
          ! up to zpbl1 and disregarded for distances larger than 3*zpbl1
          IF (z3d_in(jc,nlevs_in,jb) > 3._wp*zpbl1) THEN
            sfc_inv(jc) = 0._wp
          ELSE IF (z3d_in(jc,nlevs_in,jb) > zpbl1) THEN
            sfc_inv(jc) = sfc_inv(jc)*(1._wp - (z3d_in(jc,nlevs_in,jb)-zpbl1)/(2._wp*zpbl1))
          ENDIF

          tsfc_mod(jc) = tsfc_mod(jc) - sfc_inv(jc)

          ! Limitation of very cold temperatures according to GME method
          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))

          ! Estimated temperature at mean sea level
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z3d_in(jc,nlevs_in,jb)

          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF

        ENDDO

      ENDIF

      ! Compute vertical gradients of virtual potential temperature
      DO jk = 1, nlevs_out

        DO jc = 1, nlen           ! The vertical gradient is needed for downward extrapolation only

          IF (z3d_in(jc,nlevs_in,jb) > 1._wp) THEN
            dtvdz_down(jc,jk) = (tsfc_mod(jc)-tmsl(jc))/z3d_in(jc,nlevs_in,jb)

          ELSE ! avoid pathological results at grid points below sea level
            dtvdz_down(jc,jk) = dtdz_standardatm

          ENDIF

        ENDDO

      ENDDO

      ! Now compute pressure on target grid
      DO jk = 1, nlevs_out

        DO jc = 1, nlen
          IF (jk <= bot_idx(jc,jb)) THEN

            ! interpolation based on piecewise analytical integration of the hydrostatic equation:
            ! we use only the formula for constant (layer-averaged) temperature in order not to
            ! require an already interpolated (virtual) temperature on the output grid;
            ! the linear weighting between the upward and downward integrals nevertheless
            ! ensures negligibly small errors
            jk1 = idx0(jc,jk,jb)

            p_up = pres_in(jc,jk1,jb)*EXP(-grav*(z3d_out(jc,jk,jb)-z3d_in(jc,jk1,jb)) / &
              (rd*0.5_wp*(tempv_in(jc,jk1,jb)+tempv_in(jc,jk1+1,jb))) )

            p_down = pres_in(jc,jk1+1,jb)*EXP(-grav*(z3d_out(jc,jk,jb)-z3d_in(jc,jk1+1,jb)) / &
              (rd*0.5_wp*(tempv_in(jc,jk1,jb)+tempv_in(jc,jk1+1,jb))) )

            ! apply inverse-distance weighting between top-down and bottom-up integrated value
            ! except in the case of extrapolation above the top of the input data
            IF (jk1 == 1 .AND. wfac(jc,jk,jb) > 1._wp) THEN
              pres_out(jc,jk,jb) = p_up
            ELSE
              pres_out(jc,jk,jb) = wfac(jc,jk,jb)*p_up + (1._wp-wfac(jc,jk,jb))*p_down
            ENDIF

          ELSE ! downward extrapolation based on the vertical gradient computed above
            t_extr = tsfc_mod(jc) + (z3d_out(jc,jk,jb)-z3d_in(jc,nlevs_in,jb))*dtvdz_down(jc,jk)
            IF (ABS(dtvdz_down(jc,jk)) > dtvdz_thresh) THEN
              p_down = pres_in(jc,nlevs_in,jb)*EXP(-grav/(rd*dtvdz_down(jc,jk))* &
                LOG(t_extr/tsfc_mod(jc)) )
            ELSE
              p_down = pres_in(jc,nlevs_in,jb)*EXP(-grav*(z3d_out(jc,jk,jb)-z3d_in(jc,nlevs_in,jb)) / &
                (rd*0.5_wp*(t_extr+tsfc_mod(jc))) )
            ENDIF
            pres_out(jc,jk,jb) = p_down
          ENDIF
        ENDDO

      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE pressure_intp


  !-------------
  !>
  !! SUBROUTINE pressure_intp_initmode
  !! Performs vertical interpolation of pressure
  !!
  !! Required input fields: pressure, virtual temperature and 3D height coordinate
  !! field of input data, virtual temperature and 3D height coordinate of output data
  !! Output: pressure field of output data
  !!
  !! Method: piecewise analytical integration of the hydrostatic equation
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE pressure_intp_initmode(pres_in, tempv_in, z3d_in, pres_out, tempv_out, z3d_out, &
    &                               nblks, npromz, nlevs_out, wfac, idx0, bot_idx)


    ! Input fields
    REAL(wp), INTENT(IN)  :: pres_in  (:,:,:) ! pressure field of input data
    REAL(wp), INTENT(IN)  :: tempv_in (:,:,:) ! virtual temperature of input data
    REAL(wp), INTENT(IN)  :: z3d_in   (:,:,:) ! 3D height coordinate field of input data
    REAL(wp), INTENT(IN)  :: tempv_out(:,:,:) ! virtual temperature of output data
    REAL(wp), INTENT(IN)  :: z3d_out  (:,:,:) ! 3D height coordinate field of output data

    ! Output
    REAL(wp), INTENT(OUT) :: pres_out (:,:,:) ! pressure field of output data

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: wfac(:,:,:)    ! weighting factor of upper level
    INTEGER , INTENT(IN) :: idx0(:,:,:)    ! index of upper level
    INTEGER , INTENT(IN) :: bot_idx(:,:)   ! index of lowest level for which interpolation is possible

    ! LOCAL CONSTANTS

    REAL(wp), PARAMETER :: TOL = 1e-12_wp

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jc, jk1, nlen
    REAL(wp), DIMENSION(nproma,nlevs_out) :: dtvdz_up, dtvdz_down
    REAL(wp) :: p_up, p_down

!-------------------------------------------------------------------------

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0))) RETURN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jk1,jc,dtvdz_up,dtvdz_down,p_up, &
!$OMP            p_down) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        pres_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

      ! Compute vertical gradients of virtual potential temperature
      DO jk = 1, nlevs_out

        DO jc = 1, nlen
          IF (jk <= bot_idx(jc,jb)) THEN
            jk1 = idx0(jc,jk,jb)

            IF (ABS(z3d_in  (jc,jk1,jb)-z3d_out  (jc,jk,jb)) > TOL) THEN
              dtvdz_up(jc,jk) = (tempv_in(jc,jk1,jb)-tempv_out(jc,jk,jb)) / &
                &               (z3d_in  (jc,jk1,jb)-z3d_out  (jc,jk,jb))
            ELSE
              dtvdz_up(jc,jk) = 0._wp
            ENDIF

            ! Paranoia
            IF (ABS(z3d_in  (jc,jk1+1,jb)-z3d_out  (jc,jk,jb)) > TOL) THEN
              dtvdz_down(jc,jk) = (tempv_in(jc,jk1+1,jb)-tempv_out(jc,jk,jb)) / &
                &                 (z3d_in  (jc,jk1+1,jb)-z3d_out  (jc,jk,jb))
            ELSE
              dtvdz_down(jc,jk) = 0._wp
            END IF

          ELSE ! downward extrapolation; only dtvdz_down is needed

            dtvdz_down(jc,jk) = (tempv_out(jc,jk-1,jb)-tempv_out(jc,jk,jb)) / &
              (z3d_out(jc,jk-1,jb)-z3d_out(jc,jk,jb))

          ENDIF
        ENDDO

      ENDDO

      ! Now compute pressure on target grid
      DO jk = 1, nlevs_out

        DO jc = 1, nlen
          IF (jk <= bot_idx(jc,jb)) THEN

            ! interpolation based on piecewise analytical integration of the hydrostatic equation
            jk1 = idx0(jc,jk,jb)

            IF (ABS(dtvdz_up(jc,jk)) > dtvdz_thresh) THEN
              p_up = pres_in(jc,jk1,jb)*EXP(-grav/(rd*dtvdz_up(jc,jk))* &
                LOG(tempv_out(jc,jk,jb)/tempv_in(jc,jk1,jb)) )
            ELSE
              p_up = pres_in(jc,jk1,jb)*EXP(-grav*(z3d_out(jc,jk,jb)-z3d_in(jc,jk1,jb)) / &
                (rd*0.5_wp*(tempv_out(jc,jk,jb)+tempv_in(jc,jk1,jb))) )
            ENDIF

            IF (ABS(dtvdz_down(jc,jk)) > dtvdz_thresh) THEN
              p_down = pres_in(jc,jk1+1,jb)*EXP(-grav/(rd*dtvdz_down(jc,jk))* &
                LOG(tempv_out(jc,jk,jb)/tempv_in(jc,jk1+1,jb)) )
            ELSE
              p_down = pres_in(jc,jk1+1,jb)*EXP(-grav*(z3d_out(jc,jk,jb)-z3d_in(jc,jk1+1,jb)) / &
                (rd*0.5_wp*(tempv_out(jc,jk,jb)+tempv_in(jc,jk1+1,jb))) )
            ENDIF

            ! Finally, apply inverse-distance weighting between top-down and bottom-up integrated value
            ! except in the case of extrapolation above the top of the input data
            IF (jk1 == 1 .AND. wfac(jc,jk,jb) > 1._wp) THEN
              pres_out(jc,jk,jb) = p_up
            ELSE
              pres_out(jc,jk,jb) = wfac(jc,jk,jb)*p_up + (1._wp-wfac(jc,jk,jb))*p_down
            ENDIF

          ELSE ! downward extrapolation

            IF (ABS(dtvdz_down(jc,jk)) > dtvdz_thresh) THEN
              p_down = pres_out(jc,jk-1,jb)*EXP(-grav/(rd*dtvdz_down(jc,jk))* &
                LOG(tempv_out(jc,jk,jb)/tempv_out(jc,jk-1,jb)) )
            ELSE
              p_down = pres_out(jc,jk-1,jb)*EXP(-grav*(z3d_out(jc,jk,jb)-z3d_out(jc,jk-1,jb)) / &
                (rd*0.5_wp*(tempv_out(jc,jk,jb)+tempv_out(jc,jk-1,jb))) )
            ENDIF

            pres_out(jc,jk,jb) = p_down
          ENDIF
        ENDDO

      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE pressure_intp_initmode


  !-------------
  !>
  !! SUBROUTINE z_at_plevels
  !! Computes height for a given set of pressure levels based on the same extrapolation
  !! assumptions as used in ECMWF's IFS.
  !! The purpose of this is to prepare diagnostic output on pressure levels.
  !!
  !! Required input fields: pressure, virtual temperature and 3D height coordinate
  !! on model levels, pressure of output levels
  !! Output: 3D height field of pressure levels
  !!
  !! Method: piecewise analytical integration of the hydrostatic equation
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2012-11-13)
  !!
  !!
  !!
  SUBROUTINE z_at_plevels(pres_ml, tempv_ml, z3d_ml, pres_pl, z3d_pl, &
                          nblks, npromz, nlevs_ml, nlevs_pl,          &
                          kpbl1, wfacpbl1, kpbl2, wfacpbl2, zextrap   )


    ! Input fields
    REAL(wp), INTENT(IN)  :: pres_ml  (:,:,:) ! pressure field on model levels
    REAL(wp), INTENT(IN)  :: tempv_ml (:,:,:) ! virtual temperature on model levels
    REAL(wp), INTENT(IN)  :: z3d_ml   (:,:,:) ! 3D height coordinate field on model levels
    REAL(wp), INTENT(IN)  :: pres_pl  (:,:,:) ! pressure of output levels

    ! Comment: for z3d_zl and pres_pl, 1D fields would actually be sufficient, but having
    ! everything as 3D fields simplifies programming

    ! Output
    REAL(wp), INTENT(OUT) :: z3d_pl  (:,:,:) ! 3D height coordinate field of pressure levels

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_ml   ! Number of model levels
    INTEGER , INTENT(IN) :: nlevs_pl   ! Number of pressure levels

    ! Coefficients needed for removing the surface inversion before extrapolating downward
    INTEGER , INTENT(IN) :: kpbl1(:,:)     ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)  ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)     ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)  ! corresponding interpolation coefficient

    REAL(wp), OPTIONAL, INTENT(IN) :: zextrap(:,:)   ! AGL height from which downward extrapolation starts (in postprocesing mode)

    ! LOCAL VARIABLES

    ! Threshold for switching between analytical formulas for constant
    ! temperature and constant vertical gradient of temperature,
    ! respectively
    REAL(wp), PARAMETER :: dtvdz_thresh = 1.e-4_wp ! 0.1 K/km

    INTEGER  :: jb, jkm, jkp, jc, jkm_start
    INTEGER  :: nlen, ierror(nblks)
    INTEGER , DIMENSION(nproma)          :: bot_idx_ml
    INTEGER , DIMENSION(nproma,nlevs_pl) :: idx0_ml

    REAL(wp) :: z_up, z_down
    REAL(wp), DIMENSION(nproma,nlevs_pl) :: wfac_ml, dtvdz
    ! temporary to store pre-computed inverse of differences
    REAL(wp), DIMENSION(nproma,nlevs_ml) :: z3d_ml_di
    REAL(wp), DIMENSION(nproma)          :: tmsl, tsfc_mod, tempv1, tempv2, vtgrad_up, sfc_inv

    LOGICAL :: l_found(nproma),lfound_all,lzextrap

!-------------------------------------------------------------------------

    IF (PRESENT(zextrap)) THEN
      lzextrap = .TRUE.
    ELSE
      lzextrap = .FALSE.
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jkm,jkp,jc,nlen,jkm_start,bot_idx_ml,idx0_ml,z_up,z_down,&
!$OMP wfac_ml,tmsl,tsfc_mod,dtvdz,tempv1,tempv2,vtgrad_up,sfc_inv,l_found,&
!$OMP lfound_all,z3d_ml_di) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        z3d_pl(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF
      ierror(jb) = 0

      ! First compute coefficients for those target levels/points for which interpolation
      ! between model-level data is possible
      jkm_start = 1
      DO jkp = 1, nlevs_pl
        l_found(:) = .FALSE.
        lfound_all = .FALSE.
        DO jkm = jkm_start,nlevs_ml-1
          DO jc = 1, nlen
            IF (pres_pl(jc,jkp,jb) >= pres_ml(jc,jkm,jb) .AND. &
                pres_pl(jc,jkp,jb) <  pres_ml(jc,jkm+1,jb)) THEN
              idx0_ml(jc,jkp) = jkm
              wfac_ml(jc,jkp) = (pres_pl(jc,jkp,jb)-pres_ml(jc,jkm+1,jb))/&
                                (pres_ml(jc,jkm,jb)-pres_ml(jc,jkm+1,jb))
              bot_idx_ml(jc) = jkp
              l_found(jc) = .TRUE.
            ELSE IF (pres_pl(jc,jkp,jb) >= pres_ml(jc,nlevs_ml,jb)) THEN
              l_found(jc) = .TRUE.
              idx0_ml(jc,jkp) = nlevs_ml
              IF (jkp == 1) bot_idx_ml(jc) = 0
            ELSE IF (pres_pl(jc,jkp,jb) < pres_ml(jc,1,jb)) THEN
              idx0_ml(jc,jkp) = 1
              wfac_ml(jc,jkp) = (pres_pl(jc,jkp,jb)-pres_ml(jc,2,jb))/&
                                (pres_ml(jc,1  ,jb)-pres_ml(jc,2,jb))
              bot_idx_ml(jc) = jkp
              l_found(jc) = .TRUE.
            ENDIF
          ENDDO
          IF (ALL(l_found(1:nlen))) THEN
            lfound_all = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (lfound_all) THEN
          jkm_start = MIN(MINVAL(idx0_ml(1:nlen,jkp)),nlevs_ml-1)
        ELSE
          ierror(jb) = ierror(jb) + 1
        ENDIF
      ENDDO

      IF (ierror(jb) > 0) CALL finish("z_at_plevels:",&
        "Error in computing interpolation coefficients")

      IF (itype_pres_msl == PRES_MSL_METHOD_GME) THEN

        ! Preparations for extrapolation below the ground: calculate temperature
        ! profile with artificial limits like in IFS. This temperature profile
        ! is NOT consistent with the temperature on pressure levels computed
        ! afterwards, implying that temperature and geopotential on pressure
        ! levels are not in hydrostatic balance. We are not happy with this
        ! solution, but it is used for the time being in order to be consistent
        ! with IFS and GME.

        DO jc = 1, nlen
          tsfc_mod(jc) = tempv_ml(jc,nlevs_ml,jb)
          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z3d_ml(jc,nlevs_ml,jb)
          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF
        ENDDO

      ELSE IF ( lzextrap .AND. itype_pres_msl >= 3 ) THEN

        ! Similar method to option 1, but we use the temperature at 150 m AGL
        ! for extrapolation (kpbl1 contains the required index in this case)

        DO jc = 1, nlen

          tsfc_mod(jc) = wfacpbl1(jc,jb) *tempv_ml(jc,kpbl1(jc,jb),jb  ) + &
                  (1._wp-wfacpbl1(jc,jb))*tempv_ml(jc,kpbl1(jc,jb)+1,jb) - &
                  dtdz_standardatm*(zextrap(jc,jb)-0.5_wp*vct_a(num_lev(1)))

          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z3d_ml(jc,nlevs_ml,jb)
          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF
        ENDDO

      ELSE
        ! Use a temperature profile that is (roughly) consistent with
        ! that used for temperature extrapolation

        DO jc = 1, nlen
          ! Virtual temperature at height zpbl1
          tempv1(jc) = wfacpbl1(jc,jb) *tempv_ml(jc,kpbl1(jc,jb),jb  ) + &
                (1._wp-wfacpbl1(jc,jb))*tempv_ml(jc,kpbl1(jc,jb)+1,jb)

          ! Virtual temperature at height zpbl2
          tempv2(jc) = wfacpbl2(jc,jb) *tempv_ml(jc,kpbl2(jc,jb),jb  ) + &
                (1._wp-wfacpbl2(jc,jb))*tempv_ml(jc,kpbl2(jc,jb)+1,jb)

          ! Vertical gradient between zpbl1 and zpbl2
          vtgrad_up(jc) = (tempv2(jc) - tempv1(jc))/(zpbl2 - zpbl1)

          ! Modified (extrapolated) surface temperature without inversion
          tsfc_mod(jc) = tempv1(jc) - zpbl1*vtgrad_up(jc)

          ! "surface inversion", defined by the difference between the extrapolated
          ! extrapolated temperature from above and the original input temperature
          sfc_inv(jc) = tsfc_mod(jc) - tempv_ml(jc,nlevs_ml,jb)

          ! Reduction of the surface inversion depending on the extrapolation
          ! distance. The surface inversion is fully restored for extrapolation distances
          ! up to zpbl1 and disregarded for distances larger than 3*zpbl1
          IF (z3d_ml(jc,nlevs_ml,jb) > 3._wp*zpbl1) THEN
            sfc_inv(jc) = 0._wp
          ELSE IF (z3d_ml(jc,nlevs_ml,jb) > zpbl1) THEN
            sfc_inv(jc) = sfc_inv(jc)*(1._wp - (z3d_ml(jc,nlevs_ml,jb)-zpbl1)/(2._wp*zpbl1))
          ENDIF

          tsfc_mod(jc) = tsfc_mod(jc) - sfc_inv(jc)

          ! Limitation of very cold temperatures according to GME method
          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))

          ! Estimated temperature at mean sea level
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z3d_ml(jc,nlevs_ml,jb)

          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF

        ENDDO

      ENDIF

      DO jkp = MINVAL(bot_idx_ml(1:nlen))+1, nlevs_pl
        DO jc = 1, nlen
          IF (jkp >= bot_idx_ml(jc)+1) THEN
            ! Store extrapolation distance on wfac_ml if target point is below the
            ! surface of the input data (note: this is a positive quantity)
            idx0_ml(jc,jkp) = nlevs_ml
            wfac_ml(jc,jkp) = pres_pl(jc,jkp,jb)-pres_ml(jc,nlevs_ml,jb)
          ENDIF
        ENDDO
      ENDDO

      DO jkm = 1, nlevs_ml - 1
        DO jc = 1, nlen
          z3d_ml_di(jc, jkm) = 1.0_wp &
               / (z3d_ml(jc,jkm,jb) - z3d_ml(jc,jkm+1,jb))
        END DO
      END DO

      ! Now, compute vertical gradients of virtual potential temperature
      DO jkp = 1, nlevs_pl
        DO jc = 1, nlen
          IF (jkp <= bot_idx_ml(jc)) THEN
            jkm = idx0_ml(jc,jkp)
            dtvdz(jc,jkp) = (tempv_ml(jc,jkm,jb)-tempv_ml(jc,jkm+1,jb)) &
                 * z3d_ml_di(jc, jkm)
          ELSE IF (z3d_ml(jc,nlevs_ml,jb) > 1._wp) THEN ! extrapolation below lowest model level required
            dtvdz(jc,jkp) = (tsfc_mod(jc)-tmsl(jc))/z3d_ml(jc,nlevs_ml,jb)

          ELSE ! avoid pathological results at grid points below sea level
            dtvdz(jc,jkp) = dtdz_standardatm
          ENDIF
        ENDDO
      ENDDO

      ! Finally, compute height of pressure levels
      DO jkp = 1, nlevs_pl
        DO jc = 1, nlen
          IF (jkp <= bot_idx_ml(jc)) THEN

            ! interpolation based on piecewise analytical integration of the hydrostatic equation
            jkm = idx0_ml(jc,jkp)
            IF (ABS(dtvdz(jc,jkp)) > dtvdz_thresh) THEN
              z_up   = z3d_ml(jc,jkm,jb) + tempv_ml(jc,jkm,jb)*((pres_pl(jc,jkp,jb) /   &
                       pres_ml(jc,jkm,jb))**(-rd*dtvdz(jc,jkp)/grav)-1._wp)/dtvdz(jc,jkp)
              z_down = z3d_ml(jc,jkm+1,jb) + tempv_ml(jc,jkm+1,jb)*((pres_pl(jc,jkp,jb) / &
                       pres_ml(jc,jkm+1,jb))**(-rd*dtvdz(jc,jkp)/grav)-1._wp)/dtvdz(jc,jkp)

            ELSE
              z_up   = z3d_ml(jc,jkm,jb) + (rd*0.5_wp*(tempv_ml(jc,jkm,jb) +                 &
                       tempv_ml(jc,jkm+1,jb))/grav)*LOG(pres_ml(jc,jkm,jb)/pres_pl(jc,jkp,jb))
              z_down = z3d_ml(jc,jkm+1,jb) + (rd*0.5_wp*(tempv_ml(jc,jkm,jb) +                 &
                       tempv_ml(jc,jkm+1,jb))/grav)*LOG(pres_ml(jc,jkm+1,jb)/pres_pl(jc,jkp,jb))
            ENDIF

            ! Finally, apply inverse-distance weighting between top-down and bottom-up integrated value
            ! except in the case of extrapolation above the top of the input data
            IF (jkm == 1 .AND. wfac_ml(jc,jkp) > 1._wp) THEN
              z3d_pl(jc,jkp,jb) = z_up
            ELSE
              z3d_pl(jc,jkp,jb) = wfac_ml(jc,jkp)*z_up + (1._wp-wfac_ml(jc,jkp))*z_down
            ENDIF

          ELSE !  extrapolation below lowest height level

            jkm = nlevs_ml
            IF (ABS(dtvdz(jc,jkp)) > dtvdz_thresh) THEN
              z_up   = z3d_ml(jc,jkm,jb) + tsfc_mod(jc)*((pres_pl(jc,jkp,jb) /       &
                       pres_ml(jc,jkm,jb))**(-rd*dtvdz(jc,jkp)/grav)-1._wp)/dtvdz(jc,jkp)
            ELSE
              z_up   = z3d_ml(jc,jkm,jb) + (rd*tsfc_mod(jc)/grav) * &
                       LOG(pres_ml(jc,jkm,jb)/pres_pl(jc,jkp,jb))
            ENDIF

            z3d_pl(jc,jkp,jb) = z_up

          ENDIF

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE z_at_plevels

  !-------------
  !>
  !! SUBROUTINE z_at_theta_levels
  !! Computes height for a given set of potential temperature levels
  !! The purpose of this is to prepare diagnostic output on isentropic levels.
  !!
  !! Required input fields: potential temperature and 3D height coordinate
  !! on model levels, coefficients for computing the boundary layer temperature gradients
  !! Output: 3D height field of potential temperature levels
  !!
  !! Method: linear vertical interpolation / extrapolation
  !!
  !! Note: this routine should not be used for large extrapolation distances below
  !! the model orography. The target values of the interpolation (theta_thl)
  !! must be in descending order. It the input theta field does not increase monotonically
  !! with height, the first occurrence of the target value - scanning from top to bottom - is selected
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2012-07-19)
  !!
  !!
  SUBROUTINE z_at_theta_levels(theta_ml, z3d_ml, theta_thl, z3d_thl, wfacpbl1, kpbl1, &
                               wfacpbl2, kpbl2, nblks, npromz, nlevs_ml, nlevs_thl)


    ! Input fields
    REAL(wp), INTENT(IN)  :: theta_ml (:,:,:) ! potential temperature on model levels
    REAL(wp), INTENT(IN)  :: z3d_ml   (:,:,:) ! 3D height coordinate field on model levels
    REAL(wp), INTENT(IN)  :: theta_thl(:,:,:) ! Potential temperature on theta levels
                                              ! (target field of interpolation)

    ! Comment: for theta_thl, a 1D field would actually be sufficient, but having
    ! everything as 3D fields simplifies programming

    ! Input coefficients
    INTEGER , INTENT(IN) :: kpbl1(:,:)      ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)   ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)      ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)   ! corresponding interpolation coefficient

    ! Output
    REAL(wp), INTENT(OUT)  :: z3d_thl  (:,:,:) ! 3D height coordinate field on theta levels

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_ml   ! Number of model levels
    INTEGER , INTENT(IN) :: nlevs_thl  ! Number of theta levels

    ! LOCAL VARIABLES


    INTEGER  :: jb, jkm, jkt, jc, jkm_start
    INTEGER  :: nlen, ierror(nblks)
    INTEGER , DIMENSION(nproma)          :: bot_idx_ml
    INTEGER , DIMENSION(nproma,nlevs_thl) :: idx0_ml

    REAL(wp), DIMENSION(nproma,nlevs_thl) :: wfac_ml

    REAL(wp) :: theta1, theta2, vthgrad

    LOGICAL :: l_found(nproma),lfound_all

!-------------------------------------------------------------------------


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jkm,jkt,jc,nlen,jkm_start,bot_idx_ml,idx0_ml,wfac_ml,theta1,theta2,vthgrad, &
!$OMP            l_found,lfound_all) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        z3d_thl(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF
      ierror(jb) = 0

      ! Compute coefficients for those target levels/points for which interpolation
      ! between model-level data is possible
      jkm_start = 1
      bot_idx_ml(:) = 0
      DO jkt = 1, nlevs_thl
        l_found(:) = .FALSE.
        lfound_all = .FALSE.
        DO jkm = jkm_start,nlevs_ml-1
          DO jc = 1, nlen
            IF (theta_thl(jc,jkt,jb) <= theta_ml(jc,jkm,jb) .AND. &
                theta_thl(jc,jkt,jb) >  theta_ml(jc,jkm+1,jb)) THEN
              idx0_ml(jc,jkt) = jkm
              wfac_ml(jc,jkt) = (theta_thl(jc,jkt,jb)-theta_ml(jc,jkm+1,jb))/&
                                (theta_ml(jc,jkm,jb) -theta_ml(jc,jkm+1,jb))
              bot_idx_ml(jc) = jkt
              l_found(jc) = .TRUE.
            ! The second criterion here is needed to prevent running into the extrapolation branch
            ! if a pair of levels fulfiling the interpolation condition has already been found
            ! (relevant if theta does not increase monotonically with height)
            ELSE IF (theta_thl(jc,jkt,jb) <= theta_ml(jc,nlevs_ml,jb) .AND. bot_idx_ml(jc) < jkt ) THEN
              l_found(jc) = .TRUE.
              idx0_ml(jc,jkt) = nlevs_ml
              ! Store "extrapolation distance" on wfac_ml if target point is below the
              ! surface of the input data
              wfac_ml(jc,jkt) = theta_thl(jc,jkt,jb)-theta_ml(jc,nlevs_ml,jb)
            ELSE IF (theta_thl(jc,jkt,jb) > theta_ml(jc,1,jb)) THEN
              idx0_ml(jc,jkt) = 1
              wfac_ml(jc,jkt) = (theta_thl(jc,jkt,jb)-theta_ml(jc,2,jb))/&
                                (theta_ml (jc,1  ,jb)-theta_ml(jc,2,jb))
              bot_idx_ml(jc) = jkt
              l_found(jc) = .TRUE.
            ENDIF
          ENDDO
          IF (ALL(l_found(1:nlen))) THEN
            lfound_all = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF (lfound_all) THEN
          jkm_start = MIN(MINVAL(idx0_ml(1:nlen,jkt)),nlevs_ml-1)
        ELSE
          ierror(jb) = ierror(jb) + 1
        ENDIF
      ENDDO

      IF (ierror(jb) > 0) CALL finish("z_at_theta_levels:",&
        "Calculation of interpolation coefficients failed")

      DO jkt = 1, nlevs_thl
        DO jc = 1, nlen
          IF (jkt <= bot_idx_ml(jc)) THEN

            ! linear interpolation
            z3d_thl(jc,jkt,jb) = wfac_ml(jc,jkt)*z3d_ml(jc,idx0_ml(jc,jkt),jb) + &
              (1._wp-wfac_ml(jc,jkt))*z3d_ml(jc,idx0_ml(jc,jkt)+1,jb)

          ELSE

            ! linear extrapolation, using the gradient between (by default) 1000 m and 500 m AGL

            ! Potential temperature at height zpbl1
            theta1 = wfacpbl1(jc,jb) *theta_ml(jc,kpbl1(jc,jb),jb  ) + &
              (1._wp-wfacpbl1(jc,jb))*theta_ml(jc,kpbl1(jc,jb)+1,jb)

            ! Potential temperature at height zpbl2
            theta2 = wfacpbl2(jc,jb) *theta_ml(jc,kpbl2(jc,jb),jb  ) + &
              (1._wp-wfacpbl2(jc,jb))*theta_ml(jc,kpbl2(jc,jb)+1,jb)

            ! Vertical gradient between zpbl1 and zpbl2
            vthgrad = (theta2 - theta1)/(zpbl2 - zpbl1)

            ! Set reasonable limits
            vthgrad = MAX(vthgrad,3.0e-3_wp)
            vthgrad = MIN(vthgrad,8.5e-3_wp)

            ! wfac carries the theta difference in this case
            z3d_thl(jc,jkt,jb) = z3d_ml(jc,nlevs_ml,jb) + wfac_ml(jc,jkt)/vthgrad

          ENDIF

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE z_at_theta_levels


  !-------------
  !>
  !! SUBROUTINE compute_slope
  !! Computes the slope of the mass points, which is an optional input
  !! to the temperature interpolation routine
  !!
  !! Required input fields: topography fields, patch and interpolation state
  !! Output: slope
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE compute_slope(p_patch, p_int, topo_c, slope_c)

    TYPE(t_patch),          INTENT(IN)       :: p_patch
    TYPE(t_int_state),      INTENT(IN)       :: p_int

    ! Topography data
    REAL(wp), INTENT(IN) :: topo_c(:,:)  ! topography height of mass points

    ! Output
    REAL(wp), INTENT(OUT) :: slope_c(:,:) ! surface slope of mass points

    ! LOCAL VARIABLES

    INTEGER  :: jb, jc, je
    INTEGER  :: i_startblk, i_startidx, i_endidx, nblks_c, nblks_e
    REAL(wp), DIMENSION(nproma,1,p_patch%nblks_c) :: z_topo_c, slope_abs_c
    REAL(wp), DIMENSION(nproma,1,p_patch%nblks_v) :: z_topo_v
    REAL(wp), DIMENSION(nproma,1,p_patch%nblks_e) :: slope_norm, slope_tang, slope_abs_e

!-------------------------------------------------------------------------

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! Initialization - no slope can be computed on lateral boundary points
    slope_c(:,:)    = 0._wp

    slope_norm (:,1,:) = 0._wp
    slope_tang (:,1,:) = 0._wp
    slope_abs_e(:,1,:) = 0._wp
    z_topo_v   (:,1,:) = 0._wp

    z_topo_c(:,1,:) = topo_c(:,:)

    ! Compute auxiliary topography at verices
    CALL cells2verts_scalar(z_topo_c, p_patch, p_int%cells_aw_verts, z_topo_v, 1, 1)

    ! Compute slopes
    CALL grad_fd_norm ( z_topo_c, p_patch, slope_norm, 1, 1)
    CALL grad_fd_tang ( z_topo_v, p_patch, slope_tang, 1, 1)

    i_startblk = p_patch%edges%start_blk(2,1)

    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)

      DO je = i_startidx, i_endidx
        slope_abs_e(je,1,jb) = SQRT(slope_norm(je,1,jb)**2+slope_tang(je,1,jb)**2)
      ENDDO
    ENDDO

    ! Interpolate absolute slope to mass points
    CALL edges2cells_scalar(slope_abs_e, p_patch, p_int%e_bln_c_s, slope_abs_c, 1, 1)

    i_startblk = p_patch%cells%start_blk(2,1)

    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, i_startidx, i_endidx, 2)

      DO jc = i_startidx, i_endidx
        slope_c(jc,jb) = slope_abs_c(jc,1,jb)
      ENDDO
    ENDDO

    CALL sync_patch_array(SYNC_C,p_patch,slope_c)

  END SUBROUTINE compute_slope


  !-------------
  !>
  !! SUBROUTINE temperature_intp
  !! Performs vertical interpolation and extrapolation of temperature with
  !! special treatment of surface inversions
  !!
  !! Required input fields: 3D input field to be interpolated,
  !! coefficient fields from prepare_lin/cubic_intp and prepare_extrap
  !! Output: interpolated 3D field
  !!
  !! Performs cubic interpolation where possible, turning to linear interpolation
  !! close to the surface
  !! The most important ingredient of the refined temperature extrapolation
  !! is to remove the surface inversion, if present, before the extrapolation
  !! and to add it again afterwards with a variable weighting coefficient,
  !! accounting for the slope of the target grid point and for the height
  !! difference between source and target grid
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE temperature_intp(temp_in, temp_out, z3d_in, z3d_out,            &
                             nblks, npromz, nlevs_in, nlevs_out,             &
                             coef1, coef2, coef3, wfac_lin,                  &
                             idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,   &
                             wfacpbl1, kpbl1, wfacpbl2, kpbl2,               &
                             l_hires_corr, l_restore_sfcinv, extrapol_dist,  &
                             l_pz_mode, zextrap, slope)


    ! Atmospheric fields
    REAL(wp), INTENT(IN)  :: temp_in (:,:,:) ! input temperature field
    REAL(wp), INTENT(OUT) :: temp_out(:,:,:) ! output temperature field

    ! Coordinate fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:)   ! 3D height coordinate field of input data
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)  ! 3D height coordinate field of output data

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: coef1(:,:,:)    ! coefficient for linear term
    REAL(wp), INTENT(IN) :: coef2(:,:,:)    ! coefficient for quadratic term
    REAL(wp), INTENT(IN) :: coef3(:,:,:)    ! coefficient for cubic term
    REAL(wp), INTENT(IN) :: wfac_lin(:,:,:) ! weighting factor for linear interpolation

    INTEGER , INTENT(IN) :: idx0_cub(:,:,:) ! index of upper level, cubic interpolation
    INTEGER , INTENT(IN) :: idx0_lin(:,:,:) ! index of upper level, linear interpolation

    INTEGER , INTENT(IN) :: bot_idx_cub(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: bot_idx_lin(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)      ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)   ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)      ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)   ! corresponding interpolation coefficient

    REAL(wp), OPTIONAL, INTENT(IN) :: zextrap(:,:)    ! AGL height from which downward extrapolation starts (in postprocesing mode)

    ! Logical switch if slope-based reduction of surface inversion is to be performed
    ! (recommended when target model has a much finer resolution than source model)
    LOGICAL, INTENT(IN) :: l_hires_corr

    ! Logical switch if surface inversion is to be restored after downward extrapolation
    ! (may be set to .FALSE. when interpolating on constant height/pressure levels)
    LOGICAL, INTENT(IN) :: l_restore_sfcinv

    ! Logical switch if interpolation mode to height/pressure levels is to be used
    LOGICAL, INTENT(IN) :: l_pz_mode

    REAL(wp), INTENT(IN) :: extrapol_dist ! Maximum extrapolation distance using the local vertical gradient

    REAL(wp), INTENT(IN), OPTIONAL :: slope(:,:)  ! slope of mass points

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jk1, jc, nlen, jk_start, jk_start_in, jk_start_out, ik1(nproma)
    REAL(wp) :: wfac, sfcinv, tmsl_max

    REAL(wp), DIMENSION(nproma) :: temp1, temp2, vtgrad_up, zdiff_inout, &
                                   redinv1, redinv2, tmsl, tmsl_mod
    LOGICAL , DIMENSION(nproma) :: l_found

    REAL(wp), DIMENSION(nproma,nlevs_in)  :: zalml_in, sfc_inv, temp_mod, g1, g2, g3
    REAL(wp), DIMENSION(nproma,nlevs_in-1) :: zalml_in_d
    REAL(wp), DIMENSION(nproma,nlevs_out) :: zalml_out

!-------------------------------------------------------------------------

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0))) RETURN

    IF (l_hires_corr .AND. .NOT. PRESENT(slope)) CALL finish("temperature_intp:",&
      "slope correction requires slope data as input")

    IF (.NOT. PRESENT(zextrap) .AND. l_pz_mode .AND. itype_pres_msl >= 3) THEN
      CALL finish("pressure_intp:", "zextrap missing in argument list")
    ENDIF

    ! Artificial upper limit on sea-level temperature for so-called plateau correction
    tmsl_max = 298._wp

!$OMP PARALLEL private(jc, jk, zalml_in, zalml_out)


    do jc = 1, nproma
      do jk = 1, nlevs_in
        zalml_in(jc, jk) = (jc - 1) * nlevs_in + jk
      end do
    end do
    do jc = 1, nproma
      do jk = 1, nlevs_out
        zalml_out(jc, jk) = (jc - 1) * nlevs_out + jk
      end do
    end do

!-------------------------------------------------------------------------

!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,jk_start,jk_start_in,jk_start_out,ik1,wfac,sfcinv,    &
!$OMP            temp1,temp2,tmsl,tmsl_mod,vtgrad_up,zdiff_inout,redinv1,redinv2,l_found,&
!$OMP            zalml_in_d,temp_mod,sfc_inv,g1,g2,g3) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        temp_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

      IF (l_pz_mode .AND. itype_pres_msl < 3) THEN

        DO jk1 = 1, nlevs_in
          DO jc = 1, nlen

            ! just copy temp_in to temp_mod
            temp_mod(jc,jk1) = temp_in(jc,jk1,jb)

          ENDDO
        ENDDO

      ELSE IF (l_pz_mode) THEN ! new IFS method: extrapolate downward from 150 m AGL
                               ! with the standard atmosphere gradient

        DO jc = 1, nlen
          ! Temperature at height zpbl1 (150 m AGL in this case)
          temp1(jc) = wfacpbl1(jc,jb) *temp_in(jc,kpbl1(jc,jb),jb  ) + &
               (1._wp-wfacpbl1(jc,jb))*temp_in(jc,kpbl1(jc,jb)+1,jb)
        ENDDO

        DO jk1 = 1, nlevs_in
          DO jc = 1, nlen
            IF (jk1 <= kpbl1(jc,jb)) THEN ! just use the input temperature
              temp_mod(jc,jk1) = temp_in(jc,jk1,jb)
            ELSE ! extrapolate downward from 150 m AGL with 6.5 K/km
              temp_mod(jc,jk1) = temp1(jc) - dtdz_standardatm*(zextrap(jc,jb) -       &
                0.5_wp*vct_a(num_lev(1)) - (z3d_in(jc,jk1,jb) - z3d_in(jc,nlevs_in,jb)) )
            ENDIF
          ENDDO
        ENDDO

      ELSE
        DO jc = 1, nlen
          ! Temperature at height zpbl1
          temp1(jc) = wfacpbl1(jc,jb) *temp_in(jc,kpbl1(jc,jb),jb  ) + &
               (1._wp-wfacpbl1(jc,jb))*temp_in(jc,kpbl1(jc,jb)+1,jb)

          ! Temperature at height zpbl2
          temp2(jc) = wfacpbl2(jc,jb) *temp_in(jc,kpbl2(jc,jb),jb  ) + &
               (1._wp-wfacpbl2(jc,jb))*temp_in(jc,kpbl2(jc,jb)+1,jb)

          ! Vertical gradient between zpbl1 and zpbl2
          vtgrad_up(jc) = (temp2(jc) - temp1(jc))/(zpbl2 - zpbl1)

          ! Set reasonable limits
          vtgrad_up(jc) = MAX(vtgrad_up(jc),-8.5e-3_wp)
          vtgrad_up(jc) = MIN(vtgrad_up(jc),-1.5e-3_wp)

          ! height distance between lowest input and output grid point
          ! (negative if extrapolation takes place)
          zdiff_inout(jc) = z3d_out(jc,nlevs_out,jb) - z3d_in(jc,nlevs_in,jb)
        ENDDO

        DO jk1 = 1, nlevs_in
          DO jc = 1, nlen

            ! Height above lowest model level
            zalml_in(jc,jk1) = z3d_in(jc,jk1,jb) - z3d_in(jc,nlevs_in,jb)

            ! Modified temperature with surface inversion removed
            IF (zalml_in(jc,jk1) < zpbl1) THEN
              temp_mod(jc,jk1) = temp1(jc)+vtgrad_up(jc)*(zalml_in(jc,jk1)-zpbl1)
            ELSE
              temp_mod(jc,jk1) = temp_in(jc,jk1,jb)
            ENDIF

            ! "surface inversion", defined by the difference between the extrapolated
            ! extrapolated temperature from above and the original input temperature
            !
            ! Note: the extrapolated temperature may also be slightly colder than the original
            ! one in the presence of (super-)adiabatic temperature gradients. In such
            ! cases, using vtgrad_up for downward extrapolation is still appropriate because
            ! the atmosphere in narrow mountain valleys usually does not become adiabatic
            ! because of slope heating
            !
            sfc_inv(jc,jk1) = temp_mod(jc,jk1) - temp_in(jc,jk1,jb)

          ENDDO
        ENDDO

        DO jk1 = 1, nlevs_in
          IF (MINVAL(zalml_in(1:nlen,jk1)) < zpbl1) THEN
            jk_start_in = jk1
            EXIT
          ENDIF
        ENDDO
      ENDIF

      ! Reduce the surface inversion strength depending on the slope of the target point and
      ! on the height difference between the source and target topography.
      !
      ! The following empirical modifications are designed for the case that the target
      ! model has a significantly finer spatial resolution than the source model,
      ! implying that nocturnal surface inversions should be removed on grid points lying
      ! higher than on the source grid (i.e. mountain ranges not resolved in the source model)
      ! and on slope points because cold air drainage in reality lets the cold air
      ! accumulate at the valley bottom.

      IF (l_hires_corr) THEN
        DO jc = 1, nlen

          IF (slope(jc,jb) <= 5.e-2_wp) THEN ! 50 m/km
            redinv1(jc) = 1._wp
          ELSE IF (slope(jc,jb) <= 3.5e-1_wp) THEN ! 350 m/km
            redinv1(jc) = 1._wp - LOG(slope(jc,jb)/5.e-2_wp)/LOG(10._wp)
          ELSE
            redinv1(jc) = 1._wp - LOG(7._wp)/LOG(10._wp)
          ENDIF

          ! zdiff_inout > 0 means that the target grid point is higher than the source point
          IF (zdiff_inout(jc) <= 100._wp) THEN
            redinv2(jc) = 1._wp + MAX(0._wp,MIN(1.5_wp,-1.e-3_wp*zdiff_inout(jc)))
          ELSE IF (zdiff_inout(jc) <= 300._wp) THEN
            redinv2(jc) = 1._wp - (zdiff_inout(jc)-100._wp)/200._wp
          ELSE
            redinv2(jc) = 0._wp
          ENDIF
        ENDDO

        ! Reduce surface inversion if sfc_inv > 0 (i.e. there is really enhanced static stability)
        DO jk1 = jk_start_in, nlevs_in
          DO jc = 1, nlen
            IF (sfc_inv(jc,jk1) > 0._wp) THEN
              sfc_inv(jc,jk1) = sfc_inv(jc,jk1)*MIN(1._wp, redinv1(jc)*redinv2(jc))
            ENDIF
          ENDDO
        ENDDO

      ENDIF


      ! Compute vertical gradients of input data
      DO jk1 = 1, nlevs_in-1
        DO jc = 1, nlen
          g1(jc,jk1) = (temp_mod(jc,jk1 )-temp_mod(jc,jk1+1 ))/ &
                       (z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+1,jb))
        ENDDO
      ENDDO

      ! Compute vertical gradients of gradients
      DO jk1 = 1, nlevs_in-2
        DO jc = 1, nlen
          g2(jc,jk1) = (g1(jc,jk1)-g1(jc,jk1+1))/(z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+2,jb))
        ENDDO
      ENDDO

      ! Compute third-order vertical gradients
      DO jk1 = 1, nlevs_in-3
        DO jc = 1, nlen
          g3(jc,jk1) = (g2(jc,jk1)-g2(jc,jk1+1))/(z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+3,jb))
        ENDDO
      ENDDO

      ! Now perform vertical interpolation, based on the modified temperature field
      IF (.NOT. l_pz_mode) THEN ! mode for interpolating initial data

        DO jk = 1, nlevs_out
          DO jc = 1, nlen
            IF (jk <= bot_idx_cub(jc,jb)) THEN

              ! cubic interpolation
              jk1 = idx0_cub(jc,jk,jb)
              temp_out(jc,jk,jb) = temp_mod(jc,jk1)           + coef1(jc,jk,jb)*g1(jc,jk1) + &
                                   coef2(jc,jk,jb)*g2(jc,jk1) + coef3(jc,jk,jb)*g3(jc,jk1)

            ELSE IF (jk <= bot_idx_lin(jc,jb)) THEN

              ! linear interpolation
              jk1 = idx0_lin(jc,jk,jb)
              temp_out(jc,jk,jb) = wfac_lin(jc,jk,jb)*temp_mod(jc,jk1) + &
                (1._wp-wfac_lin(jc,jk,jb))*temp_mod(jc,jk1+1)

            ELSE

              ! linear extrapolation using the upper temperature gradient; wfac_lin
              ! carries the extrapolation distance
              IF (wfac_lin(jc,jk,jb) > extrapol_dist) THEN
                temp_out(jc,jk,jb) = temp_mod(jc,nlevs_in) + wfac_lin(jc,jk,jb)*vtgrad_up(jc)
              ELSE
                temp_out(jc,jk,jb) = temp_mod(jc,nlevs_in) + extrapol_dist*vtgrad_up(jc) +  &
                                     (wfac_lin(jc,jk,jb)-extrapol_dist)*dtdz_standardatm
              ENDIF

            ENDIF

            ! Height above lowest model level - needed for restoring the surface inversion
            zalml_out(jc,jk) = z3d_out(jc,jk,jb) - z3d_out(jc,nlevs_out,jb)

          ENDDO
        ENDDO

      ELSE ! mode for interpolation to pressure or height levels

        ! Some precomputations that do not depend on the model level
        DO jc = 1, nlen

          ! extrapolated mean sea level temperature
          tmsl(jc) = temp_mod(jc,nlevs_in) - dtdz_standardatm * z3d_in(jc,nlevs_in,jb)

          ! "Plateau correction"
          IF (z3d_in(jc,nlevs_in,jb) > 2000._wp .AND. tmsl(jc) > tmsl_max) THEN
            IF (z3d_in(jc,nlevs_in,jb) < 2500._wp) THEN
              tmsl_mod(jc) = 0.002_wp*((2500._wp-z3d_in(jc,nlevs_in,jb))*tmsl(jc)+  &
                (z3d_in(jc,nlevs_in,jb)-2000._wp)*tmsl_max)
            ELSE
              tmsl_mod(jc) = tmsl_max
            ENDIF
          ELSE
            tmsl_mod(jc) = tmsl(jc)
          ENDIF
        ENDDO

        DO jk = 1, nlevs_out
          DO jc = 1, nlen
            IF (jk <= bot_idx_cub(jc,jb)) THEN

              ! cubic interpolation
              jk1 = idx0_cub(jc,jk,jb)
              temp_out(jc,jk,jb) = temp_mod(jc,jk1)           + coef1(jc,jk,jb)*g1(jc,jk1) + &
                                   coef2(jc,jk,jb)*g2(jc,jk1) + coef3(jc,jk,jb)*g3(jc,jk1)

            ELSE IF (jk <= bot_idx_lin(jc,jb)) THEN

              ! linear interpolation
              jk1 = idx0_lin(jc,jk,jb)
              temp_out(jc,jk,jb) = wfac_lin(jc,jk,jb)*temp_mod(jc,jk1) + &
                (1._wp-wfac_lin(jc,jk,jb))*temp_mod(jc,jk1+1)

            ELSE ! apply empiricial correction computed above for extrapolation below the ground

              temp_out(jc,jk,jb) = temp_mod(jc,nlevs_in) + wfac_lin(jc,jk,jb)* &
                (temp_mod(jc,nlevs_in)-tmsl_mod(jc))/z3d_in(jc,nlevs_in,jb)

            ENDIF

          ENDDO
        ENDDO

      ENDIF

      ! Finally, subtract surface inversion from preliminary temperature field
      IF (l_restore_sfcinv .AND. .NOT. l_pz_mode) THEN

        DO jk = 1, nlevs_out
          IF (MINVAL(zalml_out(1:nlen,jk)) < zpbl1) THEN
            jk_start_out = jk
            EXIT
          ENDIF
        ENDDO

        zalml_in_d = 1.0 / (zalml_in(1:nproma,1:nlevs_in-1) - &
          &                 zalml_in(1:nproma,2:nlevs_in))

        jk_start = jk_start_in - 1
        DO jk = jk_start_out, nlevs_out
          l_found(:) = .FALSE.
          DO jk1 = jk_start, nlevs_in-1
            DO jc = 1, nlen
              IF (zalml_out(jc,jk) >= zpbl1) THEN
                l_found(jc) = .TRUE.
                ik1(jc)     = jk_start
              ELSE IF (zalml_out(jc,jk) <  zalml_in(jc,jk1) .AND. &
                       zalml_out(jc,jk) >= zalml_in(jc,jk1+1)) THEN

                wfac = (zalml_out(jc,jk)-zalml_in(jc,jk1+1))*&
                       zalml_in_d(jc,jk1)
                sfcinv = wfac*sfc_inv(jc,jk1) + (1._wp-wfac)*sfc_inv(jc,jk1+1)

                l_found(jc) = .TRUE.
                ik1(jc)     = jk1

                temp_out(jc,jk,jb) =  temp_out(jc,jk,jb) - sfcinv

              ELSE IF (zalml_out(jc,jk) > zalml_in(jc,jk_start)) THEN
                l_found(jc) = .TRUE.
                ik1(jc)     = jk_start
              ENDIF
            ENDDO
            IF (ALL(l_found(1:nlen))) EXIT
          ENDDO
          jk_start = MINVAL(ik1(1:nlen))
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE temperature_intp


  !-------------
  !>
  !! SUBROUTINE uv_intp
  !! Performs vertical interpolation and extrapolation of horizontal wind components
  !! with special treatment of boundary layer effects
  !!
  !! Required input fields: 3D input field to be interpolated,
  !! coefficient fields from prepare_lin/cubic_intp and prepare_extrap
  !! Output: interpolated 3D field
  !!
  !! Performs cubic interpolation where possible, turning to linear interpolation
  !! close to the surface
  !! Boundary-layer treatment follows a similar reasoning as for temperature,
  !! but there are different consistency checks and limitations.
  !! In particular, wind components are requested not to change sign when being
  !! extrapolated downward
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-18)
  !!
  !!
  !!
  SUBROUTINE uv_intp(uv_in, uv_out, z3d_in, z3d_out,               &
                     nblks, npromz, nlevs_in, nlevs_out,           &
                     coef1, coef2, coef3, wfac_lin,                &
                     idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin, &
                     wfacpbl1, kpbl1, wfacpbl2, kpbl2,             &
                     l_hires_intp, l_restore_fricred, extrap_limit )


    ! Atmospheric fields
    REAL(wp), INTENT(IN)  :: uv_in (:,:,:) ! input wind (component) field
    REAL(wp), INTENT(OUT) :: uv_out(:,:,:) ! output wind (component) field

    ! Coordinate fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:)   ! 3D height coordinate field of input data
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)  ! 3D height coordinate field of output data

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: coef1(:,:,:)    ! coefficient for linear term
    REAL(wp), INTENT(IN) :: coef2(:,:,:)    ! coefficient for quadratic term
    REAL(wp), INTENT(IN) :: coef3(:,:,:)    ! coefficient for cubic term
    REAL(wp), INTENT(IN) :: wfac_lin(:,:,:) ! weighting factor for linear interpolation

    INTEGER , INTENT(IN) :: idx0_cub(:,:,:) ! index of upper level, cubic interpolation
    INTEGER , INTENT(IN) :: idx0_lin(:,:,:) ! index of upper level, linear interpolation

    INTEGER , INTENT(IN) :: bot_idx_cub(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: bot_idx_lin(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)      ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)   ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)      ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)   ! corresponding interpolation coefficient

    LOGICAL,  INTENT(IN) :: l_hires_intp ! mode for interpolation to (much) finer grid
    LOGICAL,  INTENT(IN), OPTIONAL :: l_restore_fricred ! subtract/restore frictional reduction of wind speed
    REAL(wp), INTENT(IN), OPTIONAL :: extrap_limit  ! multiplicative limit in case of downward extrapolation

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jk1, jc, nlen, jk_start, jk_start_in, jk_start_out, ik1(nproma)
    REAL(wp) :: wfac, fricred, mult_limit


    REAL(wp), DIMENSION(nproma) :: uv1, uv2, dudz_up
    LOGICAL , DIMENSION(nproma) :: l_found

    LOGICAL :: lrestore_fricred

    REAL(wp), DIMENSION(nproma,nlevs_in)  :: zalml_in, fric_red, uv_mod, g1, g2, g3
    REAL(wp), DIMENSION(nproma,nlevs_in-1) :: zalml_in_d
    REAL(wp), DIMENSION(nproma,nlevs_out) :: zalml_out, zdiff_inout, red_speed

!-------------------------------------------------------------------------

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0)))  RETURN

    ! consistency check:
    IF (UBOUND(uv_out,2) < nlevs_out) CALL finish("Wrong size of output field!")

    IF (PRESENT(l_restore_fricred)) THEN
      lrestore_fricred = l_restore_fricred
    ELSE
      lrestore_fricred = .TRUE.
    ENDIF

    ! Multiplicative limit for downward extrapolation; not used in high-resolution mode
    ! (wind speed in narrow valleys not resolved in the source model may be reduced
    ! more strongly)
    IF (PRESENT(extrap_limit)) THEN
      mult_limit = extrap_limit
    ELSE IF (lrestore_fricred) THEN
      mult_limit = 0.5_wp
    ELSE
      mult_limit = 1.0_wp
    ENDIF

!$OMP PARALLEL private(jc, jk, zalml_in,zalml_out)

    do jc = 1, nproma
      do jk = 1, nlevs_in
        zalml_in(jc, jk) = (jc - 1) * nlevs_in + jk
      end do
    end do
    do jc = 1, nproma
      do jk = 1, nlevs_out
        zalml_out(jc, jk) = (jc - 1) * nlevs_out + jk
      end do
    end do

!-------------------------------------------------------------------------

!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,jk_start,jk_start_in,jk_start_out,ik1,wfac,fricred,&
!$OMP            uv1,uv2,dudz_up,zdiff_inout,red_speed,l_found,    &
!$OMP            uv_mod,zalml_in_d, fric_red,g1,g2,g3) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        uv_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

      DO jc = 1, nlen
        ! Wind component at height zpbl1
        uv1(jc) = wfacpbl1(jc,jb) *uv_in(jc,kpbl1(jc,jb),jb  ) + &
           (1._wp-wfacpbl1(jc,jb))*uv_in(jc,kpbl1(jc,jb)+1,jb)

        ! Wind component at height zpbl2
        uv2(jc) = wfacpbl2(jc,jb) *uv_in(jc,kpbl2(jc,jb),jb  ) + &
           (1._wp-wfacpbl2(jc,jb))*uv_in(jc,kpbl2(jc,jb)+1,jb)

        ! Vertical gradient between zpbl1 and zpbl2
        dudz_up(jc) = (uv2(jc) - uv1(jc))/(zpbl2 - zpbl1)
      ENDDO

      DO jk1 = 1, nlevs_in
        DO jc = 1, nlen

          ! Height above lowest model level
          zalml_in(jc,jk1) = z3d_in(jc,jk1,jb) - z3d_in(jc,nlevs_in,jb)

          ! Modified wind speed with frictional reduction removed
          IF (lrestore_fricred .AND. zalml_in(jc,jk1) < zpbl1) THEN
            uv_mod(jc,jk1) = uv1(jc)+dudz_up(jc)*(zalml_in(jc,jk1)-zpbl1)
          ELSE
            uv_mod(jc,jk1) = uv_in(jc,jk1,jb)
          ENDIF

          ! frictional reduction of wind speed
          fric_red(jc,jk1) = uv_mod(jc,jk1) - uv_in(jc,jk1,jb)

        ENDDO
      ENDDO

      DO jk1 = 1, nlevs_in
        IF (MINVAL(zalml_in(1:nlen,jk1)) < zpbl1) THEN
          jk_start_in = jk1
          EXIT
        ENDIF
      ENDDO

      ! Compute factor for artificial wind speed reduction when interpolating
      ! into deep valleys, whose presence is inferred from zdiff_inout
      ! (of course, this inference assumes that the target grid has a much finer
      ! resolution than the source grid; thus, the reduction is used for l_hires_intp only

      IF (l_hires_intp) THEN

        DO jk = 1, nlevs_out
          DO jc = 1, nlen

            ! height distance between lowest input and output grid points
            ! (negative if extrapolation takes place)
            zdiff_inout(jc,jk) = z3d_out(jc,jk,jb) - z3d_in(jc,nlevs_in,jb)

            IF (zdiff_inout(jc,jk) <= -500._wp) THEN
              red_speed(jc,jk) = 0.25_wp
            ELSE IF (zdiff_inout(jc,jk) <= -100._wp) THEN
              red_speed(jc,jk) = 1._wp + (zdiff_inout(jc,jk) + 100._wp)*1.875e-3_wp
            ELSE
              red_speed(jc,jk) = 1._wp
            ENDIF

          ENDDO
        ENDDO

      ELSE !  no artificial reduction

        red_speed(:,:) = 1._wp

      ENDIF


      ! Compute vertical gradients of input data
      DO jk1 = 1, nlevs_in-1
        DO jc = 1, nlen
          g1(jc,jk1) = (uv_mod(jc,jk1 )-uv_mod(jc,jk1+1 ))/ &
                       (z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+1,jb))
        ENDDO
      ENDDO

      ! Compute vertical gradients of gradients
      DO jk1 = 1, nlevs_in-2
        DO jc = 1, nlen
          g2(jc,jk1) = (g1(jc,jk1)-g1(jc,jk1+1))/(z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+2,jb))
        ENDDO
      ENDDO

      ! Compute third-order vertical gradients
      DO jk1 = 1, nlevs_in-3
        DO jc = 1, nlen
          g3(jc,jk1) = (g2(jc,jk1)-g2(jc,jk1+1))/(z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+3,jb))
        ENDDO
      ENDDO

      ! Now perform vertical interpolation, based on the modified temperature field
      DO jk = 1, nlevs_out
        DO jc = 1, nlen
          IF (jk <= bot_idx_cub(jc,jb)) THEN

            ! cubic interpolation
            jk1 = idx0_cub(jc,jk,jb)
            uv_out(jc,jk,jb) = uv_mod(jc,jk1)             + coef1(jc,jk,jb)*g1(jc,jk1) + &
                               coef2(jc,jk,jb)*g2(jc,jk1) + coef3(jc,jk,jb)*g3(jc,jk1)

          ELSE IF (jk <= bot_idx_lin(jc,jb)) THEN

            ! linear interpolation
            jk1 = idx0_lin(jc,jk,jb)
            uv_out(jc,jk,jb) = wfac_lin(jc,jk,jb)*uv_mod(jc,jk1) + &
              (1._wp-wfac_lin(jc,jk,jb))*uv_mod(jc,jk1+1)

          ELSE IF (lrestore_fricred) THEN

            ! downward extrapolation
            uv_out(jc,jk,jb) = uv_mod(jc,nlevs_in) + wfac_lin(jc,jk,jb)*dudz_up(jc)

            ! ensure that the extrapolated wind does not change sign and
            ! stays within the specified limit of speed reduction
            IF (uv_out(jc,jk,jb)*uv_in(jc,nlevs_in,jb) >= 0._wp) THEN ! still correct sign

              uv_out(jc,jk,jb) = SIGN(MAX(ABS(uv_out(jc,jk,jb)),  &
                ABS(mult_limit*uv_mod(jc,nlevs_in))), uv_mod(jc,nlevs_in))

              ! For extrapolation to z-levels or p-levels, this limitation is needed in addition
              uv_out(jc,jk,jb) = SIGN(MIN(ABS(uv_out(jc,jk,jb)),  &
                ABS(1.25_wp*uv_mod(jc,nlevs_in))), uv_mod(jc,nlevs_in))

            ELSE  ! i.e. spurious sign change due to extrapolation

              uv_out(jc,jk,jb) = mult_limit*uv_mod(jc,nlevs_in)

            ENDIF

            ! Artificial reduction for large downward extrapolation distances
            ! if high-resolution mode is selected
            uv_out(jc,jk,jb) = uv_out(jc,jk,jb)*red_speed(jc,jk)

          ELSE ! in post-processing mode, just use constant profile below ground

            uv_out(jc,jk,jb) = uv_mod(jc,nlevs_in)

          ENDIF

          ! Height above lowest model level - needed for restoring the frictional layer
          zalml_out(jc,jk) = z3d_out(jc,jk,jb) - z3d_out(jc,nlevs_out,jb)

        ENDDO
      ENDDO

      DO jk = 1, nlevs_out
        IF (MINVAL(zalml_out(1:nlen,jk)) < zpbl1) THEN
          jk_start_out = jk
          EXIT
        ENDIF
      ENDDO

      ! subtract frictional reduction from preliminary wind field
      IF (lrestore_fricred) THEN

        zalml_in_d = 1.0 / (zalml_in(1:nproma,1:nlevs_in-1) - &
          &                 zalml_in(1:nproma,2:nlevs_in))

        jk_start = jk_start_in - 1
        DO jk = jk_start_out, nlevs_out
          l_found(:) = .FALSE.
          DO jk1 = jk_start, nlevs_in-1
            DO jc = 1, nlen
              IF (zalml_out(jc,jk) >= zpbl1) THEN
                l_found(jc) = .TRUE.
                ik1(jc)     = jk_start
              ELSE IF (zalml_out(jc,jk) <  zalml_in(jc,jk1) .AND. &
                       zalml_out(jc,jk) >= zalml_in(jc,jk1+1)) THEN

                wfac = (zalml_out(jc,jk)-zalml_in(jc,jk1+1))*&
                       zalml_in_d(jc,jk1)
                fricred = wfac*fric_red(jc,jk1) + (1._wp-wfac)*fric_red(jc,jk1+1)

                l_found(jc) = .TRUE.
                ik1(jc)     = jk1

                uv_out(jc,jk,jb) = (uv_out(jc,jk,jb)/red_speed(jc,jk)-fricred)*red_speed(jc,jk)

              ELSE IF (zalml_out(jc,jk) > zalml_in(jc,jk_start)) THEN
                l_found(jc) = .TRUE.
                ik1(jc)     = jk_start
              ENDIF
            ENDDO
            IF (ALL(l_found(1:nlen))) EXIT
          ENDDO
          jk_start = MINVAL(ik1(1:nlen))
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE uv_intp


  !-------------
  !>
  !! SUBROUTINE qv_intp
  !! Performs vertical interpolation of specific humidity,
  !! including checks for supersaturation and unrealistically small values
  !!
  !! Required input fields: 3D input field to be interpolated,
  !! coefficient fields from prepare_lin/cubic_intp and prepare_extrap
  !! Output: interpolated 3D field
  !!
  !! Performs cubic interpolation where possible, turning to linear interpolation
  !! close to the surface
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-18)
  !!
  !!
  !!
  SUBROUTINE qv_intp(qv_in, qv_out, z3d_in, z3d_out,                  &
                     temp_in, pres_in, temp_out, pres_out,            &
                     nblks, npromz, nlevs_in, nlevs_out,              &
                     coef1, coef2, coef3, wfac_lin,                   &
                     idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,    &
                     wfacpbl1, kpbl1, wfacpbl2, kpbl2,                &
                     lower_limit, l_satlimit, l_restore_pbldev, opt_qc)


    ! Specific humidity fields
    REAL(wp), INTENT(IN)  :: qv_in (:,:,:) ! input field
    REAL(wp), INTENT(OUT) :: qv_out(:,:,:) ! output field

    ! Optional cloud water field: if provided, consistency checks between QV and QC will
    ! be performed at the end of this routine
    REAL(wp), INTENT(INOUT), OPTIONAL :: opt_qc(:,:,:) ! specific cloud water

    ! Additional atmospheric fields on input and output grids (needed for supersaturation check)
    ! Note: these are all input data
    REAL(wp), INTENT(IN) :: temp_in(:,:,:)  ! temperature on input grid
    REAL(wp), INTENT(IN) :: pres_in(:,:,:)  ! pressure on input grid
    REAL(wp), INTENT(IN) :: temp_out(:,:,:) ! temperature on output grid
    REAL(wp), INTENT(IN) :: pres_out(:,:,:) ! pressure on output grid

    ! Coordinate fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:)   ! 3D height coordinate field of input data
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)  ! 3D height coordinate field of output data

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: coef1(:,:,:)    ! coefficient for linear term
    REAL(wp), INTENT(IN) :: coef2(:,:,:)    ! coefficient for quadratic term
    REAL(wp), INTENT(IN) :: coef3(:,:,:)    ! coefficient for cubic term
    REAL(wp), INTENT(IN) :: wfac_lin(:,:,:) ! weighting factor for linear interpolation

    INTEGER , INTENT(IN) :: idx0_cub(:,:,:) ! index of upper level, cubic interpolation
    INTEGER , INTENT(IN) :: idx0_lin(:,:,:) ! index of upper level, linear interpolation

    INTEGER , INTENT(IN) :: bot_idx_cub(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: bot_idx_lin(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)      ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)   ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)      ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)   ! corresponding interpolation coefficient

    REAL(wp), INTENT(IN) :: lower_limit     ! lower limit of QV
    LOGICAL , INTENT(IN) :: l_satlimit       ! limit input field to water saturation
    LOGICAL , INTENT(IN) :: l_restore_pbldev ! restore PBL deviation of QV from extrapolated profile

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jk1, jc, nlen, jk_start, jk_start_in, jk_start_out, ik1(nproma)
    REAL(wp) :: wfac, pbldev, rhum, qtot
    REAL(wp) :: e_vapor                     ! vapor pressure


    REAL(wp), DIMENSION(nproma) :: qv1, qv2, dqvdz_up
    LOGICAL , DIMENSION(nproma) :: l_found
    LOGICAL                     :: l_check_qv_qc

    REAL(wp), DIMENSION(nproma,nlevs_in)  :: zalml_in, pbl_dev, qv_mod, g1, g2, g3, qsat_in, qv_in_lim
    REAL(wp), DIMENSION(nproma,nlevs_in-1) :: zalml_in_d
    REAL(wp), DIMENSION(nproma,nlevs_out) :: zalml_out, qsat_out

!-------------------------------------------------------------------------

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0))) RETURN

    IF (PRESENT(opt_qc)) THEN
      l_check_qv_qc = .TRUE.
    ELSE
      l_check_qv_qc = .FALSE.
    ENDIF

!$OMP PARALLEL private(jc, jk, zalml_in, zalml_out)

    do jc = 1, nproma
      do jk = 1, nlevs_in
        zalml_in(jc, jk) = (jc - 1) * nlevs_in + jk
      end do
    end do
    do jc = 1, nproma
      do jk = 1, nlevs_out
        zalml_out(jc, jk) = (jc - 1) * nlevs_out + jk
      end do
    end do

!-------------------------------------------------------------------------

!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,jk_start,jk_start_in,jk_start_out,ik1,wfac,pbldev,&
!$OMP            rhum,qtot,qv1,qv2,dqvdz_up,l_found,zalml_in_d,qv_in_lim,    &
!$OMP            qv_mod,pbl_dev,g1,g2,g3,e_vapor,qsat_in,qsat_out) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        qv_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF


      DO jk1 = 1, nlevs_in
        DO jc = 1, nlen

          ! Limit input data to the 'lower_limit' specified in the argument list
          qv_in_lim(jc,jk1) = MAX(qv_in(jc,jk1,jb),lower_limit)

          ! compute vapour pressure e_vapor=f(qv_in,pres_in)
          e_vapor = vap_pres(qv_in_lim(jc,jk1),pres_in(jc,jk1,jb))

          ! saturation specific humidity of input data
          qsat_in(jc,jk1) = rdv*sat_pres_water(temp_in(jc,jk1,jb)) /    &
                            (pres_in(jc,jk1,jb)-o_m_rdv*e_vapor)

          IF (l_satlimit) THEN
            ! limit input data to water saturation when processing interpolated IFS data:
            ! This is needed to remove supersaturations generated (primarily) by the interpolation
            ! from the spherical harmonics to the Gaussain grid; without this limitation, the
            ! QV-QC-adjustment at the end of this routine would generate nonsensically large
            ! cloud water peaks
            qv_in_lim(jc,jk1) = MIN(qv_in_lim(jc,jk1),qsat_in(jc,jk1))
          ENDIF
        ENDDO
      ENDDO



      DO jc = 1, nlen
        ! QV at height zpbl1
        qv1(jc) = wfacpbl1(jc,jb) *qv_in_lim(jc,kpbl1(jc,jb)  ) + &
           (1._wp-wfacpbl1(jc,jb))*qv_in_lim(jc,kpbl1(jc,jb)+1)

        ! QV at height zpbl2
        qv2(jc) = wfacpbl2(jc,jb) *qv_in_lim(jc,kpbl2(jc,jb)  ) + &
           (1._wp-wfacpbl2(jc,jb))*qv_in_lim(jc,kpbl2(jc,jb)+1)

        ! Vertical gradient between zpbl1 and zpbl2
        dqvdz_up(jc) = (qv2(jc) - qv1(jc))/(zpbl2 - zpbl1)
      ENDDO

      DO jk1 = 1, nlevs_in
        DO jc = 1, nlen

          ! Height above lowest model level
          zalml_in(jc,jk1) = z3d_in(jc,jk1,jb) - z3d_in(jc,nlevs_in,jb)

          ! Modified QV with boundary-layer deviation from the extrapolated profile removed
          IF (zalml_in(jc,jk1) < zpbl1) THEN
            qv_mod(jc,jk1) = qv1(jc)+dqvdz_up(jc)*(zalml_in(jc,jk1)-zpbl1)
          ELSE
            qv_mod(jc,jk1) = qv_in_lim(jc,jk1)
          ENDIF

          ! boundary-layer deviation of QV, converted to RH
          pbl_dev(jc,jk1) = (qv_mod(jc,jk1) - qv_in_lim(jc,jk1)) / qsat_in(jc,jk1)

        ENDDO
      ENDDO

      DO jk1 = 1, nlevs_in
        IF (MINVAL(zalml_in(1:nlen,jk1)) < zpbl1) THEN
          jk_start_in = jk1
          EXIT
        ENDIF
      ENDDO


      ! Compute vertical gradients of input data
      DO jk1 = 1, nlevs_in-1
        DO jc = 1, nlen
          g1(jc,jk1) = (qv_mod(jc,jk1 )-qv_mod(jc,jk1+1 ))/ &
                       (z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+1,jb))
        ENDDO
      ENDDO

      ! Compute vertical gradients of gradients
      DO jk1 = 1, nlevs_in-2
        DO jc = 1, nlen
          g2(jc,jk1) = (g1(jc,jk1)-g1(jc,jk1+1))/(z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+2,jb))
        ENDDO
      ENDDO

      ! Compute third-order vertical gradients
      DO jk1 = 1, nlevs_in-3
        DO jc = 1, nlen
          g3(jc,jk1) = (g2(jc,jk1)-g2(jc,jk1+1))/(z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+3,jb))
        ENDDO
      ENDDO

      ! Now perform vertical interpolation, based on the modified temperature field
      DO jk = 1, nlevs_out
        DO jc = 1, nlen

          IF (jk <= bot_idx_cub(jc,jb)) THEN

            ! cubic interpolation
            jk1 = idx0_cub(jc,jk,jb)
            qv_out(jc,jk,jb) = qv_mod(jc,jk1)             + coef1(jc,jk,jb)*g1(jc,jk1) + &
                               coef2(jc,jk,jb)*g2(jc,jk1) + coef3(jc,jk,jb)*g3(jc,jk1)

          ELSE IF (jk <= bot_idx_lin(jc,jb)) THEN

            ! linear interpolation
            jk1 = idx0_lin(jc,jk,jb)
            qv_out(jc,jk,jb) = wfac_lin(jc,jk,jb)*qv_mod(jc,jk1) + &
              (1._wp-wfac_lin(jc,jk,jb))*qv_mod(jc,jk1+1)

          ELSE

            ! if extrapolation is needed, maintain relative humidity (with a slight
            ! approximation of the qsat expression to avoid iterations)
            qv_out(jc,jk,jb) = qv_mod(jc,nlevs_in)/qsat_in(jc,nlevs_in)* &
              rdv*sat_pres_water(temp_out(jc,jk,jb))/pres_out(jc,jk,jb)

          ENDIF

          ! Height above lowest model level - needed for restoring the boundary-layer deviation
          zalml_out(jc,jk) = z3d_out(jc,jk,jb) - z3d_out(jc,nlevs_out,jb)

          ! compute vapour pressure e_vapor=f(qv_out,pres_out)
          e_vapor = vap_pres(qv_out(jc,jk,jb),pres_out(jc,jk,jb))

          ! saturation specific humidity of output data
          qsat_out(jc,jk) = rdv*sat_pres_water(temp_out(jc,jk,jb)) /    &
                            (pres_out(jc,jk,jb)-o_m_rdv*e_vapor)

        ENDDO
      ENDDO

      DO jk = 1, nlevs_out
        IF (MINVAL(zalml_out(1:nlen,jk)) < zpbl1) THEN
          jk_start_out = jk
          EXIT
        ENDIF
      ENDDO

      zalml_in_d = 1.0 / (zalml_in(1:nproma,1:nlevs_in-1) - &
        &                 zalml_in(1:nproma,2:nlevs_in))

      ! subtract boundary-layer deviation from preliminary QV
      IF (l_restore_pbldev) THEN
        jk_start = jk_start_in - 1
        DO jk = jk_start_out, nlevs_out
          l_found(:) = .FALSE.
          DO jk1 = jk_start, nlevs_in-1
            DO jc = 1, nlen
              IF (zalml_out(jc,jk) >= zpbl1) THEN
                l_found(jc) = .TRUE.
                ik1(jc)     = jk_start
              ELSE IF (zalml_out(jc,jk) <  zalml_in(jc,jk1) .AND. &
                       zalml_out(jc,jk) >= zalml_in(jc,jk1+1)) THEN

                wfac = (zalml_out(jc,jk)-zalml_in(jc,jk1+1))*&
                       zalml_in_d(jc,jk1)

                pbldev = wfac*pbl_dev(jc,jk1) + (1._wp-wfac)*pbl_dev(jc,jk1+1)

                l_found(jc) = .TRUE.
                ik1(jc)     = jk1

                qv_out(jc,jk,jb) =  qv_out(jc,jk,jb) - pbldev*qsat_out(jc,jk)

              ELSE IF (zalml_out(jc,jk) > zalml_in(jc,jk_start)) THEN
                l_found(jc) = .TRUE.
                ik1(jc)     = jk_start
              ENDIF
            ENDDO
            IF (ALL(l_found(1:nlen))) EXIT
          ENDDO
          jk_start = MINVAL(ik1(1:nlen))
        ENDDO
      ENDIF

      IF (l_check_qv_qc) THEN ! apply consistency checks between QV and QC
        DO jk = 1, nlevs_out
          DO jc = 1, nlen

            ! Impose limiter to avoid negative or unreasonably small values,
            ! and ensure that RH > 1% below 200 hPa
            qv_out(jc,jk,jb) = MAX(lower_limit,qv_out(jc,jk,jb))

            IF (pres_out(jc,jk,jb) >= 20000._wp) THEN
              qv_out(jc,jk,jb) = MAX(qv_out(jc,jk,jb),0.01_wp*qsat_out(jc,jk))
            ENDIF

            rhum = qv_out(jc,jk,jb)/qsat_out(jc,jk)

            IF (rhum < 0.75_wp) THEN ! remove all the cloud water if present
              opt_qc(jc,jk,jb) = 0._wp
            ELSE ! conserve QV plus QC
              qtot = qv_out(jc,jk,jb) + opt_qc(jc,jk,jb)
              opt_qc(jc,jk,jb) = MAX(0._wp,qtot-qsat_out(jc,jk))
              qv_out(jc,jk,jb) = qtot - opt_qc(jc,jk,jb)
            ENDIF

          ENDDO
        ENDDO

      ELSE ! impose only some obvious constraints on QV
        DO jk = 1, nlevs_out
          DO jc = 1, nlen

            ! Impose limiter to avoid negative or unreasonably small values
            qv_out(jc,jk,jb) = MAX(lower_limit,qv_out(jc,jk,jb))

            ! in addition, limit output data to water saturation,
            ! ensure that RH > 1% below 200 hPa
            qv_out(jc,jk,jb) = MIN(qv_out(jc,jk,jb),qsat_out(jc,jk))

            IF (pres_out(jc,jk,jb) >= 20000._wp) THEN
              qv_out(jc,jk,jb) = MAX(qv_out(jc,jk,jb),0.01_wp*qsat_out(jc,jk))
            ENDIF

          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE qv_intp


  !------------------------------------------------------------------
  !> SUBROUTINE diagnose_pmsl_gme
  !  Compute mean sea level pressure from the surface pressure.
  !
  ! @par Revision History
  ! Initial version by F. Prill, DWD (2012-07-16)
  !
  ! Based on GME's "pp_extrap.f90" which is a modified GM routine by
  ! J.Haseler (ECMWF).
  !
  ! Cf. also
  ! Trenberth, K.E., Berry, J.C., Buja, L.E. 1993:
  ! "Vertical Interpolation and Truncation of Model-coordinate Data."
  !
  ! We assume an atmosphere in hydrostatic equilibrium and having
  ! a constant lapse rate.
  !
  ! Method
  !
  ! For temperature, geopotential and mean sea-level pressure, the
  ! surface temperature, "ztstar", and the msl temperature, "ztmsl"
  ! are calculated, assuming a lapse rate of 6.5 degrees *C/km.  If
  ! both "ztstar" and "ztmsl" are greater than 290.5, "ztmsl" is
  ! corrected according to the formula:
  !
  !    ztmsl = 290.5-0.005*(ztstar-290.5)**2
  !
  ! Comment on this method (from Trenberth et al.):
  ! "The above method seems to work reasonably well as long as the
  ! surface height is less than about 2000m elevation. Problems remain
  ! especially in the Himalayan-Tibetan Plateau region. At very high
  ! elevations, the critical dependence on the local surface
  ! temperature that is used to set the temperature of the whole
  ! artificial below-ground column becomes amplified and leads to
  ! considerable noise in the mean sea level pressure."
  !
  SUBROUTINE diagnose_pmsl_gme(pres3d_in, pres_sfc_in, temp3d_in, z3d_in, pres_out, &
    &                          nblks, npromz, nlevs_in)

    ! Input fields
    REAL(wp), INTENT(IN)  :: pres3d_in   (:,:,:) ! pressure field of input data
    REAL(wp), INTENT(IN)  :: pres_sfc_in (:,:)   ! surface pressure field (input data)
    REAL(wp), INTENT(IN)  :: temp3d_in   (:,:,:) ! temperature of input data
    REAL(wp), INTENT(IN)  :: z3d_in      (:,:,:) ! 3D height coordinate field of input data

    ! Output
    REAL(wp), INTENT(OUT) :: pres_out   (:,:)    ! mean sea level pressure field

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! local variables
!!$    CHARACTER(*), PARAMETER :: routine = TRIM("mo_nh_vert_interp:diagnose_pmsl_gme")
    REAL(wp), PARAMETER :: zthird = 1._wp/3._wp, &
      &                    zlapse = 6.5E-3_wp
    INTEGER             :: nlev, nlevp1, nlen, jb, jc
    REAL(wp)            :: zrg, ztstar, ztmsl, zalph, zprt, zprtal, &
      &                    geop_sfc, temp_in, pres_in, pres_sfc

    !  INITIALISATION
    !  --------------

    nlev   = nlevs_in
    nlevp1 = nlev + 1
    zrg    = 1.0_wp/grav

    ! TEMPERATURE-RELATED FIELDS
    ! --------------------------

    ! initialize output field
    pres_out(:,nblks) = 0._wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jc,geop_sfc,temp_in,pres_in,pres_sfc,ztstar,ztmsl,zalph,zprt,zprtal) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
      ENDIF

      DO jc = 1, nlen

        geop_sfc = z3d_in(jc,nlevp1,jb)*grav    ! surface geopotential
        temp_in  = temp3d_in(jc,nlev,jb)        ! temperature at lowest (full) model level
        pres_in  = pres3d_in(jc,nlev,jb)        ! pressure at lowest (full) level
        pres_sfc = pres_sfc_in(jc,jb)           ! surface pressure

        ! Surface pressure is given on half levels, temperature is given
        ! on full levels. Therefore extrapolate temperature with
        ! zlapse=6.5E-3 to model surface.
        !
        ! This formula follows from a Taylor series expansion of
        !   temp/ztstar  =  (pres_sfc/pres_in)**(R*lapse/grav)
        ! for small increments of p:
        ztstar = ( 1._wp + zlapse*rd*zrg*(pres_sfc/pres_in-1._wp) )*temp_in

        ! slight modifications for exceptionally warm or cold
        ! conditions:
        IF ( ztstar < t_low) THEN
          ztstar = 0.5_wp*(t_low+ztstar)
        END IF

        ztmsl = ztstar + zlapse*zrg*geop_sfc

        IF ( ztmsl > t_high) THEN
          IF ( ztstar > t_high) THEN
            ztstar = 0.5_wp*(t_high+ztstar)
            ztmsl  = ztstar
          ELSE
            ztmsl  = t_high
          END IF
        END IF

        IF ( ABS(ztmsl-ztstar) < 1.E-6_wp) THEN
          zalph = 0._wp
        ELSE IF ( ABS(zrg*geop_sfc) > 1.E-4_wp) THEN
          zalph= rd*(ztmsl-ztstar)/geop_sfc
        ELSE
          zalph= rd*zlapse*zrg
        END IF

        ! MEAN SEA LEVEL PRESSURE
        ! -----------------------

        ! The formula
        !   pres_msl = pres_sfc ( 1 + (zalph*geop_sfc)/(rd*ztstar) )**(1/zalph)
        ! is approximated (using the log-power series) as follows:
        zprt   = geop_sfc/(rd*ztstar)
        zprtal = zprt*zalph
        pres_out(jc,jb) = pres_sfc*EXP( zprt * ( 1.0_wp-zprtal*(0.5_wp-zthird*zprtal) ) )

      END DO ! jc
    END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE diagnose_pmsl_gme



  !------------------------------------------------------------------
  !> SUBROUTINE diagnose_pmsl_ifs
  !  Compute mean sea level pressure from the surface pressure.
  !
  ! @par Revision History
  ! Initial version by Guenther Zaengl, DWD (2013-04-11)
  !
  ! Based on IFS's "ppmer.f90"
  !
  ! The method is basically the same as for the GME routine above, with the
  ! exception that the surface temperature used for extrapolating the pressure
  ! is based on the temperature at 150 m AGL rather than that at the lowest
  ! model level. Moreover, the artificial corrections applied to extremely cold
  ! temperatures are slightly different (remark: they are not self-consistent
  ! for geopotential and sea-level pressure, and are again different to the
  ! temperatures used for temperature extrapolation)
  !
  !
  SUBROUTINE diagnose_pmsl_ifs(pres_sfc_in, temp3d_in, z3d_in, pmsl_out,      &
    &                          nblks, npromz, nlevs_in, wfac_extrap, kextrap, &
    &                          zextrap, method )

    ! Input fields
    REAL(wp), INTENT(IN)  :: pres_sfc_in (:,:)   ! surface pressure field (input data)
    REAL(wp), INTENT(IN)  :: temp3d_in   (:,:,:) ! temperature of input data
    REAL(wp), INTENT(IN)  :: z3d_in      (:,:,:) ! 3D height coordinate field of input data

    ! Output
    REAL(wp), INTENT(OUT) :: pmsl_out   (:,:)    ! mean sea level pressure field

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! Coefficients
    INTEGER , INTENT(IN) :: kextrap(:,:)     ! index of model level immediately above (by default) 150 m AGL
    REAL(wp), INTENT(IN) :: wfac_extrap(:,:) ! corresponding interpolation coefficient
    REAL(wp), INTENT(IN) :: zextrap(:,:)     ! AGL height from which downward extrapolation starts (in postprocesing mode)

    ! Method (PRES_MSL_METHOD_IFS or PRES_MSL_METHOD_IFS_CORR)
    INTEGER , INTENT(IN) :: method           ! valid input: PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR

    ! local variables
    CHARACTER(*), PARAMETER :: routine = TRIM("mo_nh_vert_interp:diagnose_pmsl_ifs")

    INTEGER             :: nlev, nlevp1, nlen, jb, jc

    REAL(wp)            :: temp_extrap, vtgrad, dtdz_thresh
    REAL(wp), DIMENSION(nproma) :: tsfc, tmsl, tsfc_mod, tmsl_mod

!-------------------------------------------------------------------------

    nlev   = nlevs_in
    nlevp1 = nlev + 1

    ! Threshold for switching between analytical formulas for constant temperature and
    ! constant vertical gradient of temperature, respectively
    dtdz_thresh = 1.e-4_wp ! 0.1 K/km


    ! sanity check
    !
    IF (method < 3) THEN
      WRITE(message_text,'(a)') 'invalid value for itype_pres_msl! must be >= 3'
      CALL finish(routine, message_text)
    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jc,temp_extrap,vtgrad,tsfc,tmsl,tsfc_mod,tmsl_mod) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        pmsl_out(nlen+1:nproma,jb)  = 0.0_wp
      ENDIF

      DO jc = 1, nlen

        ! Auxiliary temperature at 150 m AGL from which the sfc temp is extrapolated
        temp_extrap = wfac_extrap(jc,jb) *temp3d_in(jc,kextrap(jc,jb),jb  ) + &
               (1._wp-wfac_extrap(jc,jb))*temp3d_in(jc,kextrap(jc,jb)+1,jb)

        ! extrapolated surface temperature (tstar in IFS nomenclature)
        tsfc(jc) = temp_extrap - dtdz_standardatm*zextrap(jc,jb)

        ! extrapolated sea-level temperature (t0 in IFS nomenclature)
        tmsl(jc) = tsfc(jc) - dtdz_standardatm*z3d_in(jc,nlevp1,jb)

        ! artificial modifications in the presence of extreme temperatures
        IF (tsfc(jc) < t_low) THEN
          tsfc_mod(jc) = 0.5_wp*(tsfc(jc)+t_low)
        ELSE IF (tsfc(jc) < t_high) THEN
          tsfc_mod(jc) = tsfc(jc)
        ELSE
          tsfc_mod(jc) = 0.5_wp*(tsfc(jc)+t_high)
        ENDIF

        tmsl_mod(jc) = tsfc_mod(jc) - dtdz_standardatm*z3d_in(jc,nlevp1,jb)

        IF (tsfc_mod(jc) < t_high .AND. tmsl_mod(jc) > t_high) THEN
          tmsl_mod(jc) = t_high
        ELSE IF (tsfc_mod(jc) >= t_high .AND. tmsl_mod(jc) > tsfc_mod(jc)) THEN
          tmsl_mod(jc) = tsfc_mod(jc)
        ELSE
          IF (method == PRES_MSL_METHOD_IFS) THEN
            tmsl_mod(jc) = tmsl(jc) ! this is missing in the geopotential computation!!!
          ENDIF
          ! tmsl_mod unmodified for method == PRES_MSL_METHOD_IFS_CORR
        ENDIF

        vtgrad = (tsfc_mod(jc) - tmsl_mod(jc))/SIGN(MAX(1.e-4_wp,ABS(z3d_in(jc,nlevp1,jb))),z3d_in(jc,nlevp1,jb))

        IF (ABS(vtgrad) > dtdz_thresh) THEN
          pmsl_out(jc,jb) = pres_sfc_in(jc,jb)*EXP(-grav/(rd*vtgrad)* &
                            LOG(tmsl_mod(jc)/tsfc_mod(jc)) )
        ELSE
          pmsl_out(jc,jb) = pres_sfc_in(jc,jb)*EXP(grav*z3d_in(jc,nlevp1,jb)/ &
                           (rd*0.5_wp*(tsfc_mod(jc)+tmsl_mod(jc))) )
        ENDIF


      END DO ! jc
    END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE diagnose_pmsl_ifs


  !-------------
  !>
  !! SUBROUTINE diagnose_pmsl
  !! Diagnoses sea-level pressure
  !!
  !! Required input fields: pressure, virtual temperature and 3D height coordinate field
  !! Output: sea-level pressure
  !!
  !! Method: stepwise analytical integration of the hydrostatic equation
  !!  step 1: from surface level of input data 500 m (zpbl1) upward
  !!  step 2: from there to zpbl1 meters above sea level
  !!  step 3: from there down to sea level
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2012-07-17)
  !!
  !!
  !!
  SUBROUTINE diagnose_pmsl(pres_in, tempv_in, z3d_in, pmsl_out, nblks, npromz, nlevs_in, &
                           wfacpbl1, kpbl1, wfacpbl2, kpbl2, z_target, extrapol_dist)


    ! Input fields
    REAL(wp), INTENT(IN)  :: pres_in  (:,:,:) ! pressure field on main levels
    REAL(wp), INTENT(IN)  :: tempv_in (:,:,:) ! virtual temperature on main levels
    REAL(wp), INTENT(IN)  :: z3d_in   (:,:,:) ! 3D height coordinate field of main levels

    ! Output
    REAL(wp), INTENT(OUT) :: pmsl_out (:,:) ! sea-level pressure

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! Coefficients
    INTEGER , INTENT(IN) :: kpbl1(:,:)     ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)  ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)     ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)  ! corresponding interpolation coefficient

    REAL(wp), INTENT(IN) :: z_target       ! target height of pressure diagnosis (usually 0 m = sea level)
    REAL(wp), INTENT(IN) :: extrapol_dist  ! maximum extrapolation distance up to which the local vertical
                                           ! temperature gradient is used

    ! LOCAL VARIABLES

    INTEGER  :: jb, jc
    INTEGER  :: nlen
    REAL(wp), DIMENSION(nproma) :: tempv1, tempv2, vtgrad_pbl_in, vtgrad_up, zdiff_inout, &
                                   sfc_inv, tempv_out_pbl1, vtgrad_up_out, tempv_out,     &
                                   vtgrad_pbl_out, p_pbl1, p_pbl1_out
    REAL(wp) :: dtvdz_thresh

!-------------------------------------------------------------------------

    ! Threshold for switching between analytical formulas for constant temperature and
    ! constant vertical gradient of temperature, respectively
    dtvdz_thresh = 1.e-4_wp ! 0.1 K/km

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jc,nlen,tempv1,tempv2,vtgrad_pbl_in,vtgrad_up,zdiff_inout,sfc_inv, &
!$OMP            tempv_out_pbl1,vtgrad_up_out,tempv_out,vtgrad_pbl_out,p_pbl1,         &
!$OMP            p_pbl1_out) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        pmsl_out(nlen+1:nproma,jb)  = 0.0_wp
      ENDIF

      DO jc = 1, nlen
        ! Virtual temperature at height zpbl1
        tempv1(jc) = wfacpbl1(jc,jb) *tempv_in(jc,kpbl1(jc,jb),jb  ) + &
              (1._wp-wfacpbl1(jc,jb))*tempv_in(jc,kpbl1(jc,jb)+1,jb)

        ! Virtual temperature at height zpbl2
        tempv2(jc) = wfacpbl2(jc,jb) *tempv_in(jc,kpbl2(jc,jb),jb  ) + &
              (1._wp-wfacpbl2(jc,jb))*tempv_in(jc,kpbl2(jc,jb)+1,jb)

        ! Vertical gradient between surface level zpbl1
        vtgrad_pbl_in(jc) = (tempv1(jc) - tempv_in(jc,nlevs_in,jb))/zpbl1

        ! Vertical gradient between zpbl1 and zpbl2
        vtgrad_up(jc) = (tempv2(jc) - tempv1(jc))/(zpbl2 - zpbl1)

        ! Set reasonably restrictive limits
        vtgrad_up(jc) = MAX(vtgrad_up(jc),-8.0e-3_wp)
        vtgrad_up(jc) = MIN(vtgrad_up(jc),-5.0e-3_wp)

        ! height distance between lowest input level
        ! (negative if model grid point lies above z_target)
        zdiff_inout(jc) = z_target - z3d_in(jc,nlevs_in,jb)

        ! "surface inversion", defined by the difference between the extrapolated
        ! extrapolated temperature from above and the original input temperature
        sfc_inv(jc) = tempv1(jc)-vtgrad_up(jc)*zpbl1 - tempv_in(jc,nlevs_in,jb)

        ! Reduction of the surface inversion depending on the extrapolation
        ! distance. The surface inversion is fully restored for extrapolation distances
        ! up to zpbl1 and disregarded for distances larger than 3*zpbl1
        IF (ABS(zdiff_inout(jc)) > 3._wp*zpbl1) THEN
          sfc_inv(jc) = 0._wp
        ELSE IF (ABS(zdiff_inout(jc)) > zpbl1) THEN
          sfc_inv(jc) = sfc_inv(jc)*(1._wp - (ABS(zdiff_inout(jc))-zpbl1)/(2._wp*zpbl1))
        ENDIF

        ! Linear extrapolation to zpbl1 above target height using the upper temperature gradient
        ! extrapol_dist is the maximum extrapolation distance up to which the local
        ! vertical gradient is considered; for larger distance, the standard atmosphere
        ! vertical gradient is taken
        IF (zdiff_inout(jc) > extrapol_dist) THEN
          tempv_out_pbl1(jc) = tempv1(jc) + zdiff_inout(jc)*vtgrad_up(jc)
          vtgrad_up_out(jc)  = vtgrad_up(jc)
        ELSE
          tempv_out_pbl1(jc) = tempv1(jc) + extrapol_dist*vtgrad_up(jc) +   &
                            (zdiff_inout(jc)-extrapol_dist)*dtdz_standardatm

          ! Artificial limitation analogous to GME method
          IF (tempv_out_pbl1(jc) < 255.65_wp) THEN
            tempv_out_pbl1(jc) = 0.5_wp*(tempv_out_pbl1(jc)+255.65_wp)
          ELSE IF (tempv_out_pbl1(jc) > 290.65_wp) THEN
            tempv_out_pbl1(jc) = 0.5_wp*(tempv_out_pbl1(jc)+290.65_wp)
            tempv_out_pbl1(jc) = MIN(298.15_wp,tempv_out_pbl1(jc))
            tempv_out_pbl1(jc) = MAX(248.15_wp,tempv_out_pbl1(jc))
          ENDIF

          ! Averaged vertical temperature gradient
          vtgrad_up_out(jc) = (tempv_out_pbl1(jc) - tempv1(jc))/zdiff_inout(jc)
        ENDIF

        ! Final extrapolation of temperature to target height, including restored surface inversion
        tempv_out(jc) = tempv_out_pbl1(jc) - vtgrad_up_out(jc)*zpbl1 - sfc_inv(jc)

        ! Boundary-layer vertical gradient of extrapolation profile
        vtgrad_pbl_out(jc) = (tempv_out_pbl1(jc)-tempv_out(jc))/zpbl1


        ! Three-step vertical integration of pressure
        ! step 1: from surface level of input data to zpbl1 (to get rid of surface inversion)

        IF (ABS(vtgrad_pbl_in(jc)) > dtvdz_thresh) THEN
          p_pbl1(jc) = pres_in(jc,nlevs_in,jb)*EXP(-grav/(rd*vtgrad_pbl_in(jc))* &
                       LOG(tempv1(jc)/tempv_in(jc,nlevs_in,jb)) )
        ELSE
          p_pbl1(jc) = pres_in(jc,nlevs_in,jb)*EXP(-grav*zpbl1 / &
                      (rd*0.5_wp*(tempv1(jc)+tempv_in(jc,nlevs_in,jb))) )
        ENDIF

        ! step 2: from zpbl1 of input data to zpbl1 above target level
        IF (ABS(vtgrad_up_out(jc)) > dtvdz_thresh) THEN
          p_pbl1_out(jc) = p_pbl1(jc)*EXP(-grav/(rd*vtgrad_up_out(jc))* &
                           LOG(tempv_out_pbl1(jc)/tempv1(jc)) )
        ELSE
          p_pbl1_out(jc) = p_pbl1(jc)*EXP(-grav*zdiff_inout(jc)/ &
                          (rd*0.5_wp*(tempv1(jc)+tempv_out_pbl1(jc))) )
        ENDIF

        ! step 3: from zpbl1 above target level to target level
        IF (ABS(vtgrad_pbl_out(jc)) > dtvdz_thresh) THEN
          pmsl_out(jc,jb) = p_pbl1_out(jc)*EXP(-grav/(rd*vtgrad_pbl_out(jc))* &
                            LOG(tempv_out(jc)/tempv_out_pbl1(jc)) )
        ELSE
          pmsl_out(jc,jb) = p_pbl1_out(jc)*EXP(grav*zpbl1/ &
                           (rd*0.5_wp*(tempv_out(jc)+tempv_out_pbl1(jc))) )
        ENDIF

      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE diagnose_pmsl


END MODULE mo_nh_vert_interp
