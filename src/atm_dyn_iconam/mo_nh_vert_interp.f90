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
MODULE mo_nh_vert_interp

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_interpolation,       ONLY: t_int_state, edges2cells_scalar
  USE mo_parallel_config,     ONLY: nproma 
  USE mo_physical_constants,  ONLY: grav, rd, rdv, o_m_rdv
  USE mo_grid_config,         ONLY: n_dom
  USE mo_exception,           ONLY: finish
  USE mo_prepicon_nml,        ONLY: nlev_in, zpbl1, zpbl2, &
                                    i_oper_mode, l_w_in, l_sfc_in
  USE mo_prepicon_utils,      ONLY: t_prepicon_state, nzplev
  USE mo_ifs_coord,           ONLY: half_level_pressure, full_level_pressure, &
                                    auxhyb, geopot
  USE mo_nh_init_utils,       ONLY: interp_uv_2_vn, init_w, adjust_w, convert_thdvars, & 
                                    virtual_temp, convert_omega2w
  USE mo_math_operators,      ONLY: grad_fd_norm, grad_fd_tang
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c
  USE mo_grf_interpolation,   ONLY: t_gridref_state
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf, interpol2_vec_grf
  USE mo_sync,                ONLY: sync_patch_array, SYNC_C, SYNC_E
  USE mo_satad,               ONLY: sat_pres_water
  USE mo_nwp_sfc_interp,      ONLY: process_sfcfields

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: vertical_interpolation, interpolate_to_p_and_z_levels

CONTAINS


  !-------------
  !>
  !! SUBROUTINE interpolate_to_p_and_z_levels
  !! Outer driver routine for vertical interpolation of analysis 
  !! data interpolated horizontally by IFS2ICON to the ICON grid
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE interpolate_to_p_and_z_levels(p_patch, p_int, prepicon)

    TYPE(t_patch),          INTENT(IN)       :: p_patch(:)
    TYPE(t_int_state),      INTENT(IN)       :: p_int(:)
    TYPE(t_prepicon_state), INTENT(INOUT)    :: prepicon(:)

    ! LOCAL VARIABLES
    INTEGER :: jg

!-------------------------------------------------------------------------

    DO jg = 1, n_dom

      CALL intp2pzlevs(p_patch(jg), p_int(jg), prepicon(jg))

    ENDDO

  END SUBROUTINE interpolate_to_p_and_z_levels

  !-------------
  !>
  !! SUBROUTINE vertical_interpolation
  !! Outer driver routine for vertical interpolation of analysis 
  !! data interpolated horizontally by IFS2ICON to the ICON grid
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE vertical_interpolation(p_patch, p_int, p_grf, prepicon)

    TYPE(t_patch),          INTENT(IN)       :: p_patch(:)
    TYPE(t_int_state),      INTENT(IN)       :: p_int(:)
    TYPE(t_gridref_state),  INTENT(IN)       :: p_grf(:)
    TYPE(t_prepicon_state), INTENT(INOUT)    :: prepicon(:)

    ! LOCAL VARIABLES
    INTEGER :: jg, jn, jgc

!-------------------------------------------------------------------------

    DO jg = 1, n_dom

      IF (i_oper_mode == 2) nlev_in = p_patch(jg)%nlev

      CALL vert_interp(p_patch(jg), p_int(jg), prepicon(jg))

      ! Apply boundary interpolation for u and v because the outer nest boundary
      ! points would remain undefined otherwise
      DO jn = 1, p_patch(jg)%n_childdom

        jgc = p_patch(jg)%child_id(jn)

        CALL interpol_scal_grf (p_patch(jg), p_patch(jgc), p_int(jg), p_grf(jg)%p_dom(jn), &
                                jn, 1, prepicon(jg)%atm%w, prepicon(jgc)%atm%w )

        CALL interpol2_vec_grf (p_patch(jg), p_patch(jgc), p_int(jg), p_grf(jg)%p_dom(jn), &
                                jn, prepicon(jg)%atm%vn, prepicon(jgc)%atm%vn )

      ENDDO
    ENDDO

  END SUBROUTINE vertical_interpolation

  !-------------
  !>
  !! SUBROUTINE vert_interp
  !! Domain-wise driver routine for vertical interpolation of analysis 
  !! data interpolated horizontally by IFS2ICON to the ICON grid
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE vert_interp(p_patch, p_int, prepicon)


    TYPE(t_patch),          INTENT(IN)       :: p_patch
    TYPE(t_int_state),      INTENT(IN)       :: p_int
    TYPE(t_prepicon_state), INTENT(INOUT)    :: prepicon


    ! LOCAL VARIABLES

    INTEGER :: jb, jc, jk
    INTEGER :: nlen, nlev, nlevp1

    ! Auxiliary fields for input data
    REAL(wp), DIMENSION(nproma,nlev_in+1) :: pres_ic, lnp_ic, geop_ic
    REAL(wp), DIMENSION(nproma,nlev_in  ) :: delp, rdelp, rdlnpr, rdalpha, geop_mc
    REAL(wp), DIMENSION(nproma,nlev_in,p_patch%nblks_c) :: temp_v_in

    ! Auxiliary field for output data
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: &
      z_tempv

    ! Auxiliary fields for coefficients
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: &
      wfac_lin, coef1, coef2, coef3

    REAL(wp), DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_c) :: &
      wfac_lin_w

    INTEGER , DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: &
      idx0_lin, idx0_cub

    INTEGER , DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_c) :: &
      idx0_lin_w

    REAL(wp), DIMENSION(nproma,p_patch%nblks_c) :: &
      wfacpbl1, wfacpbl2, slope

    INTEGER , DIMENSION(nproma,p_patch%nblks_c) :: &
      bot_idx_lin, bot_idx_lin_w, bot_idx_cub, kpbl1, kpbl2


!-------------------------------------------------------------------------

    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Compute virtual temperature of input data
    CALL virtual_temp(p_patch, prepicon%atm_in%temp, prepicon%atm_in%qv, prepicon%atm_in%qc, &
                      prepicon%atm_in%qi, prepicon%atm_in%qr, prepicon%atm_in%qs, temp_v_in  )

    ! 1. Compute pressure and height of input data, using the IFS routines

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, nlen, pres_ic, lnp_ic, geop_ic, delp, rdelp, rdlnpr, &
!$OMP            rdalpha, geop_mc)

    DO jb = 1,p_patch%nblks_c

      IF (jb /= p_patch%nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = p_patch%npromz_c
      ENDIF

      ! Check if psfc is really psfc or LOG(psfc)
      IF (MAXVAL(prepicon%atm_in%psfc(1:nlen,jb)) <= 100._wp) THEN
        prepicon%atm_in%psfc(1:nlen,jb) = EXP(prepicon%atm_in%psfc(1:nlen,jb))
      ENDIF

      CALL half_level_pressure(prepicon%atm_in%psfc(:,jb), nproma, nlen, pres_ic)

      CALL full_level_pressure(pres_ic,nproma, nlen, prepicon%atm_in%pres(:,:,jb))

      CALL auxhyb(pres_ic, nproma, nlen,              & ! in
                  delp, rdelp, lnp_ic, rdlnpr, rdalpha) ! out

      CALL geopot(temp_v_in(:,:,jb), rdlnpr, rdalpha, prepicon%atm_in%phi_sfc(:,jb), & ! in
                  nproma, 1, nlen, geop_mc, geop_ic                                  ) ! inout

      ! Compute 3D height coordinate field
      prepicon%atm_in%z3d(1:nlen,1:nlev_in,jb) = geop_mc(1:nlen,1:nlev_in)/grav

    ENDDO
!$OMP END DO
!$OMP END PARALLEL


    ! Prepare interpolation coefficients
    CALL prepare_lin_intp(prepicon%atm_in%z3d, prepicon%z_mc,               &
                          p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                          wfac_lin, idx0_lin, bot_idx_lin)

    CALL prepare_extrap(prepicon%atm_in%z3d,                        &
                        p_patch%nblks_c, p_patch%npromz_c, nlev_in, &
                        kpbl1, wfacpbl1, kpbl2, wfacpbl2 )


    CALL prepare_cubic_intp(prepicon%atm_in%z3d, prepicon%z_mc,               &
                            p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                            coef1, coef2, coef3, idx0_cub, bot_idx_cub)

    ! Perform vertical interpolation

    ! Temperature 
    CALL compute_slope(p_patch, p_int, prepicon%topography_c, prepicon%topography_v, slope)

    CALL temperature_intp(prepicon%atm_in%temp, prepicon%atm%temp,          &
                          prepicon%atm_in%z3d, prepicon%z_mc,               &
                          p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                          coef1, coef2, coef3, wfac_lin,                    &
                          idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                          wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                          l_restore_sfcinv=.TRUE., l_hires_corr=.FALSE.,    &
                          extrapol_dist=-1500._wp, slope=slope              )

    ! horizontal wind components
    CALL uv_intp(prepicon%atm_in%u, prepicon%atm%u,                &
                 prepicon%atm_in%z3d, prepicon%z_mc,               &
                 p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                 coef1, coef2, coef3, wfac_lin,                    &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                 l_hires_intp=.FALSE.                              )
    CALL uv_intp(prepicon%atm_in%v, prepicon%atm%v,                &
                 prepicon%atm_in%z3d, prepicon%z_mc,               &
                 p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                 coef1, coef2, coef3, wfac_lin,                    &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                 l_hires_intp=.FALSE.                              )

    ! Preliminary interpolation of QV: this is needed to compute virtual temperature
    ! below, which in turn is required to integrate the hydrostatic equation
    ! A lower limit of 2.5 ppm is imposed
    CALL lin_intp(prepicon%atm_in%qv, prepicon%atm%qv,                 &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,    &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,    &
                  wfacpbl2, kpbl2, l_loglin=.TRUE., l_extrapol=.TRUE., &
                  l_pd_limit=.TRUE., lower_limit=2.5e-6_wp             )

    ! Cloud and precipitation variables - linear interpolation only because cubic may 
    ! cause negative values, and no-gradient condition for downward extrapolation
    ! Positive definite limiter is used in all cases
    CALL lin_intp(prepicon%atm_in%qc, prepicon%atm%qc,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE.)

    CALL lin_intp(prepicon%atm_in%qi, prepicon%atm%qi,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE.)

    CALL lin_intp(prepicon%atm_in%qr, prepicon%atm%qr,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE.)

    CALL lin_intp(prepicon%atm_in%qs, prepicon%atm%qs,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE.)

    ! Compute virtual temperature with preliminary QV
    CALL virtual_temp(p_patch, prepicon%atm%temp, prepicon%atm%qv, prepicon%atm%qc, &
                      prepicon%atm%qi, prepicon%atm%qr, prepicon%atm%qs, z_tempv)

    ! Interpolate pressure on ICON grid
    CALL pressure_intp(prepicon%atm_in%pres, temp_v_in, prepicon%atm_in%z3d, &
                       prepicon%atm%pres, z_tempv, prepicon%z_mc,            &
                       p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,     &
                       wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1      )
 
    CALL qv_intp(prepicon%atm_in%qv, prepicon%atm%qv, prepicon%atm_in%z3d,  &
                 prepicon%z_mc, prepicon%atm_in%temp, prepicon%atm_in%pres, &
                 prepicon%atm%temp, prepicon%atm%pres,                      &
                 p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,          &
                 coef1, coef2, coef3, wfac_lin,                             &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,              &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2,                          &
                 lower_limit=2.5e-6_wp, l_restore_pbldev=.TRUE.,            &
                 opt_qc=prepicon%atm%qc ) ! ... for consistency checking

    ! Compute virtual temperature with final QV
    CALL virtual_temp(p_patch, prepicon%atm%temp, prepicon%atm%qv, prepicon%atm%qc, &
                      prepicon%atm%qi, prepicon%atm%qr, prepicon%atm%qs, z_tempv)

    ! Final interpolation of pressure on ICON grid
    CALL pressure_intp(prepicon%atm_in%pres, temp_v_in, prepicon%atm_in%z3d, &
                       prepicon%atm%pres, z_tempv, prepicon%z_mc,            &
                       p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,     &
                       wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1      )


    ! Convert thermodynamic variables into set of NH prognostic variables
    CALL convert_thdvars(p_patch, prepicon%atm%pres, z_tempv,    &
                         prepicon%atm%rho, prepicon%atm%exner,   &
                         prepicon%atm%theta_v                    )

    ! Convert u and v on cell points to vn at edge points
    CALL interp_uv_2_vn(p_patch, p_int, prepicon%atm%u, prepicon%atm%v, prepicon%atm%vn)

    CALL sync_patch_array(SYNC_E,p_patch,prepicon%atm%vn)

    IF (l_w_in) THEN
      ! convert omega to w
      CALL convert_omega2w(prepicon%atm_in%omega, prepicon%atm_in%w,   &
                           prepicon%atm_in%pres, prepicon%atm_in%temp, &
                           p_patch%nblks_c, p_patch%npromz_c, nlev_in  )

      ! Compute coefficients for w interpolation
      CALL prepare_lin_intp(prepicon%atm_in%z3d, prepicon%z_ifc,                &
                            p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlevp1, &
                            wfac_lin_w, idx0_lin_w, bot_idx_lin_w)

      ! Perform linear interpolation of w
      ! Note: the coefficients for gradient computation (*pbl*) do not have to be recomputed
      ! because of l_extrapol=.FALSE., 
      CALL lin_intp(prepicon%atm_in%w, prepicon%atm%w,                     &
                    p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlevp1,    &
                    wfac_lin_w, idx0_lin_w, bot_idx_lin_w, wfacpbl1, kpbl1,&
                    wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                    l_pd_limit=.FALSE.)

      ! Impose appropriate lower boundary condition on vertical wind field
      CALL adjust_w(p_patch, p_int, prepicon%atm%vn, prepicon%z_ifc, prepicon%atm%w)
    ELSE
      ! Initialize vertical wind field
      CALL init_w(p_patch, p_int, prepicon%atm%vn, prepicon%z_ifc, prepicon%atm%w)
    ENDIF

    CALL sync_patch_array(SYNC_C,p_patch,prepicon%atm%w)

    IF (l_sfc_in) THEN 
      ! process surface fields
      CALL process_sfcfields(p_patch, prepicon)
    ENDIF

  END SUBROUTINE vert_interp

  !-------------
  !>
  !! SUBROUTINE intp2pzlevs
  !! Domain-wise driver routine for vertical interpolation of model-level 
  !! fields to constant-height and constant-pressure levels
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-18)
  !!
  !!
  SUBROUTINE intp2pzlevs(p_patch, p_int, prepicon)


    TYPE(t_patch),          INTENT(IN)       :: p_patch
    TYPE(t_int_state),      INTENT(IN)       :: p_int
    TYPE(t_prepicon_state), INTENT(INOUT)    :: prepicon


    ! LOCAL VARIABLES

    INTEGER :: jb, jk
    INTEGER :: nlen, nlev

    ! Auxiliary field for input data
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: z_tempv_in

    ! Auxiliary field for output data
    REAL(wp), DIMENSION(nproma,nzplev,p_patch%nblks_c) :: z_tempv

    ! Auxiliary fields for coefficients
    REAL(wp), DIMENSION(nproma,nzplev,p_patch%nblks_c) :: &
      wfac_lin, coef1, coef2, coef3

    INTEGER , DIMENSION(nproma,nzplev,p_patch%nblks_c) :: &
      idx0_lin, idx0_cub

    REAL(wp), DIMENSION(nproma,p_patch%nblks_c) :: &
      wfacpbl1, wfacpbl2

    INTEGER , DIMENSION(nproma,p_patch%nblks_c) :: &
      bot_idx_lin, bot_idx_cub, kpbl1, kpbl2


!-------------------------------------------------------------------------

    nlev = p_patch%nlev

    ! Fill z3d field of pressure-level data and pressure field of height-level data

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1,p_patch%nblks_c

      IF (jb /= p_patch%nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = p_patch%npromz_c
      ENDIF

      DO jk = 1, nzplev
        prepicon%plev%pres(1:nlen,jk,jb) = prepicon%plev%levels(jk)
        prepicon%zlev%z3d(1:nlen,jk,jb)  = prepicon%zlev%levels(jk)
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Part 1: Interpolation to z-level fields

    ! Prepare interpolation coefficients
    CALL prepare_lin_intp(prepicon%z_mc, prepicon%zlev%z3d,                &
                          p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev, &
                          wfac_lin, idx0_lin, bot_idx_lin)

    CALL prepare_extrap(prepicon%z_mc,                           &
                        p_patch%nblks_c, p_patch%npromz_c, nlev, &
                        kpbl1, wfacpbl1, kpbl2, wfacpbl2 )


    CALL prepare_cubic_intp(prepicon%z_mc, prepicon%zlev%z3d,                &
                            p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev, &
                            coef1, coef2, coef3, idx0_cub, bot_idx_cub)

    ! Perform vertical interpolation


    CALL temperature_intp(prepicon%atm%temp, prepicon%zlev%temp,           &
                          prepicon%z_mc, prepicon%zlev%z3d,                &
                          p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev, &
                          coef1, coef2, coef3, wfac_lin,                   &
                          idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,    &
                          wfacpbl1, kpbl1, wfacpbl2, kpbl2,                &
                          l_restore_sfcinv=.FALSE., l_hires_corr=.FALSE.,  &
                          extrapol_dist=-500._wp                           )


    ! horizontal wind components
    CALL uv_intp(prepicon%atm%u, prepicon%zlev%u,                  &
                 prepicon%z_mc, prepicon%zlev%z3d,                 &
                 p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev,  &
                 coef1, coef2, coef3, wfac_lin,                    &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                 l_hires_intp=.FALSE., l_restore_fricred=.FALSE.   )
    CALL uv_intp(prepicon%atm%v, prepicon%zlev%v,                  &
                 prepicon%z_mc, prepicon%zlev%z3d,                 &
                 p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev,  &
                 coef1, coef2, coef3, wfac_lin,                    &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                 l_hires_intp=.FALSE., l_restore_fricred=.FALSE.   )

    ! Preliminary interpolation of QV; a lower limit of 2.5 ppm is imposed
    CALL lin_intp(prepicon%atm%qv, prepicon%zlev%qv,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev,     &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,    &
                  wfacpbl2, kpbl2, l_loglin=.TRUE., l_extrapol=.TRUE., &
                  l_pd_limit=.TRUE., lower_limit=2.5e-6_wp             )

    ! Compute virtual temperature for model-level and z-level data
    CALL virtual_temp(p_patch, prepicon%atm%temp, prepicon%atm%qv, temp_v=z_tempv_in)
    CALL virtual_temp(p_patch, prepicon%zlev%temp, prepicon%zlev%qv, temp_v=z_tempv)

    ! Interpolate pressure on z-levels
    CALL pressure_intp(prepicon%atm%pres, z_tempv_in, prepicon%z_mc,     &
                       prepicon%zlev%pres, z_tempv, prepicon%zlev%z3d,   &
                       p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev,  &
                       wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1  )

    ! Final interpolation of QV, including supersaturation limiting
    CALL qv_intp(prepicon%atm%qv, prepicon%zlev%qv, prepicon%z_mc,        &
                 prepicon%zlev%z3d, prepicon%atm%temp, prepicon%atm%pres, &
                 prepicon%zlev%temp, prepicon%zlev%pres,                  &
                 p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev,         &
                 coef1, coef2, coef3, wfac_lin,                           &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,            &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2,                        &
                 lower_limit=2.5e-6_wp, l_restore_pbldev=.FALSE.          )

    ! Part 2: Interpolation to pressure-level fields


    ! Compute height at pressure levels; this height field is afterwards also
    ! used as target coordinate for vertical interpolation
    CALL z_at_plevels(prepicon%atm%pres, z_tempv_in, prepicon%z_mc,           &
                      prepicon%zlev%pres, z_tempv, prepicon%zlev%z3d,         &
                      prepicon%plev%pres, prepicon%plev%z3d, p_patch%nblks_c, &
                      p_patch%npromz_c, nlev, nzplev, nzplev                  )


    ! Prepare again interpolation coefficients (now for pressure levels)
    CALL prepare_lin_intp(prepicon%z_mc, prepicon%plev%z3d,                &
                          p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev, &
                          wfac_lin, idx0_lin, bot_idx_lin)

    CALL prepare_cubic_intp(prepicon%z_mc, prepicon%plev%z3d,                &
                            p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev, &
                            coef1, coef2, coef3, idx0_cub, bot_idx_cub)

    ! Perform vertical interpolation


    CALL temperature_intp(prepicon%atm%temp, prepicon%plev%temp,           &
                          prepicon%z_mc, prepicon%plev%z3d,                &
                          p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev, &
                          coef1, coef2, coef3, wfac_lin,                   &
                          idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,    &
                          wfacpbl1, kpbl1, wfacpbl2, kpbl2,                &
                          l_restore_sfcinv=.FALSE., l_hires_corr=.FALSE.,  &
                          extrapol_dist=-500._wp                           )

    ! horizontal wind components
    CALL uv_intp(prepicon%atm%u, prepicon%plev%u,                  &
                 prepicon%z_mc, prepicon%plev%z3d,                 &
                 p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev,  &
                 coef1, coef2, coef3, wfac_lin,                    &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                 l_hires_intp=.FALSE., l_restore_fricred=.FALSE.   )
    CALL uv_intp(prepicon%atm%v, prepicon%plev%v,                  &
                 prepicon%z_mc, prepicon%plev%z3d,                 &
                 p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev,  &
                 coef1, coef2, coef3, wfac_lin,                    &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                 l_hires_intp=.FALSE., l_restore_fricred=.FALSE.   )

    ! Interpolation of QV, including supersaturation limiting
    CALL qv_intp(prepicon%atm%qv, prepicon%plev%qv, prepicon%z_mc,        &
                 prepicon%plev%z3d, prepicon%atm%temp, prepicon%atm%pres, &
                 prepicon%plev%temp, prepicon%plev%pres,                  &
                 p_patch%nblks_c, p_patch%npromz_c, nlev, nzplev,         &
                 coef1, coef2, coef3, wfac_lin,                           &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,            &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2,                        &
                 lower_limit=2.5e-6_wp, l_restore_pbldev=.FALSE.          )

  END SUBROUTINE intp2pzlevs


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
                              wfac, idx0, bot_idx                 )

    ! Input fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:) ! height coordinate field of input data (m)
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)! height coordinate field of input data (m)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Output fields
    REAL(wp), INTENT(OUT) :: wfac(:,:,:)       ! weighting factor of upper level
    INTEGER , INTENT(OUT) :: idx0(:,:,:)       ! index of upper level
    INTEGER , INTENT(OUT) :: bot_idx(:,:)      ! index of lowest level for which interpolation is possible

    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jk1, jk_start
    INTEGER :: nlen
    LOGICAL :: l_found(nproma)

!-------------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk1,jk_start,l_found)

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        bot_idx(nlen+1:nproma,jb) = nlevs_out
        wfac(nlen+1:nproma,:,jb)  = 0.5_wp
        idx0(nlen+1:nproma,:,jb)  = nlevs_in-1
      ENDIF

      jk_start = 1
      DO jk = 1, nlevs_out
        l_found(:) = .FALSE.
        DO jk1 = jk_start,nlevs_in-1
          DO jc = 1, nlen
            IF (z3d_out(jc,jk,jb) <= z3d_in(jc,jk1,jb) .AND. &
                z3d_out(jc,jk,jb) >  z3d_in(jc,jk1+1,jb)) THEN
              idx0(jc,jk,jb) = jk1
              wfac(jc,jk,jb) = (z3d_out(jc,jk,jb)-z3d_in(jc,jk1+1,jb))/&
                               (z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+1,jb))
              bot_idx(jc,jb) = jk
              l_found(jc) = .TRUE.
            ELSE IF (z3d_out(jc,jk,jb) < z3d_in(jc,nlevs_in,jb)) THEN
              l_found(jc) = .TRUE.
              idx0(jc,jk,jb) = nlevs_in
            ENDIF
          ENDDO
          IF (ALL(l_found(1:nlen))) EXIT
        ENDDO
        jk_start = MINVAL(idx0(1:nlen,jk,jb))
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
!$OMP END DO
!$OMP END PARALLEL

    IF (MINVAL(bot_idx) == 0) CALL finish("prepare_lin_intp:",&
      "ICON top higher than top of input data")

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
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk_start)

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

      DO jk = jk_start, nlevs_in
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
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE prepare_extrap

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
    INTEGER :: nlen
    LOGICAL :: l_found(nproma)

!-------------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk1,jk_start,l_found)

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

      jk_start = 2
      DO jk = 1, nlevs_out
        l_found(:) = .FALSE.
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
            ELSE IF (z3d_out(jc,jk,jb) < z3d_in(jc,nlevs_in-1,jb)) THEN
              l_found(jc) = .TRUE.
              idx0(jc,jk,jb) = nlevs_in-1
            ENDIF
          ENDDO
          IF (ALL(l_found(1:nlen))) EXIT
        ENDDO
        jk_start = MINVAL(idx0(1:nlen,jk,jb))+1
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
!$OMP END DO
!$OMP END PARALLEL

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
    REAL(wp), INTENT(INOUT) :: f3d_in (:,:,:) ! input field (INOUT because of limiter)
    REAL(wp), INTENT(OUT)   :: f3d_out(:,:,:) ! output field (on ICON vertical grid)

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
    REAL(wp) :: zf_in(nproma,nlevs_in), z_limit, f3d_z1, f3d_z2, vgrad_f3d

!-------------------------------------------------------------------------

    IF (PRESENT(lower_limit)) THEN
      z_limit = lower_limit
    ELSE
      z_limit = 0._wp
    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,zf_in,f3d_z1,f3d_z2,vgrad_f3d)

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        f3d_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

      IF (l_pd_limit) THEN
        f3d_in(1:nlen,1:nlevs_in,jb) = MAX(z_limit,f3d_in(1:nlen,1:nlevs_in,jb))
      ENDIF

      IF (l_loglin) THEN
        zf_in(1:nlen,1:nlevs_in) = LOG(f3d_in(1:nlen,1:nlevs_in,jb))
      ELSE
        zf_in(1:nlen,1:nlevs_in) =     f3d_in(1:nlen,1:nlevs_in,jb)
      ENDIF

      DO jk = 1, nlevs_out
        DO jc = 1, nlen
          IF (jk <= bot_idx(jc,jb)) THEN

            ! linear interpolation
            f3d_out(jc,jk,jb) = wfac(jc,jk,jb)*zf_in(jc,idx0(jc,jk,jb)) + &
              (1._wp-wfac(jc,jk,jb))*zf_in(jc,idx0(jc,jk,jb)+1)

            IF (l_loglin) f3d_out(jc,jk,jb) = EXP(f3d_out(jc,jk,jb))

          ELSE IF (l_extrapol) THEN

            ! linear extrapolation, using the gradient between (by default) 1000 m and 500 m AGL
            ! Logarithmic computation is not used here because it would be numerically unstable for extrapolation

            ! Field value at height zpbl1
            f3d_z1 = wfacpbl1(jc,jb) *f3d_in(jc,kpbl1(jc,jb)  ,jb) +  &
              (1._wp-wfacpbl1(jc,jb))*f3d_in(jc,kpbl1(jc,jb)+1,jb)

            ! Field value at height zpbl2
            f3d_z2 = wfacpbl2(jc,jb) *f3d_in(jc,kpbl2(jc,jb)  ,jb) +  &
              (1._wp-wfacpbl2(jc,jb))*f3d_in(jc,kpbl2(jc,jb)+1,jb)

            ! vertical gradient
            vgrad_f3d = (f3d_z2-f3d_z1)/(zpbl2-zpbl1)

            ! wfac carries the (negative) extrapolation distance (in m) in this case
            f3d_out(jc,jk,jb) = f3d_in(jc,nlevs_in,jb) + vgrad_f3d*wfac(jc,jk,jb)

          ELSE ! use no-gradient condition for extrapolation

            f3d_out(jc,jk,jb) = f3d_in(jc,nlevs_in,jb)

          ENDIF
        ENDDO
      ENDDO

      IF (l_pd_limit) THEN
        f3d_out(1:nlen,1:nlevs_out,jb) = MAX(z_limit,f3d_out(1:nlen,1:nlevs_out,jb))
      ENDIF

    ENDDO
!$OMP END DO
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
  SUBROUTINE pressure_intp(pres_in, tempv_in, z3d_in, pres_out, tempv_out, z3d_out, &
                           nblks, npromz, nlevs_in, nlevs_out,                      &
                           wfac, idx0, bot_idx, wfacpbl1, kpbl1)


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
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: wfac(:,:,:)    ! weighting factor of upper level
    INTEGER , INTENT(IN) :: idx0(:,:,:)    ! index of upper level
    INTEGER , INTENT(IN) :: bot_idx(:,:)   ! index of lowest level for which interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)     ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)  ! corresponding interpolation coefficient

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jc, jk1
    INTEGER  :: nlen
    REAL(wp), DIMENSION(nproma,nlevs_out) :: dtvdz_up, dtvdz_down
    REAL(wp) :: dtvdz_thresh, p_up, p_down

!-------------------------------------------------------------------------


    ! Threshold for switching between analytical formulas for constant temperature and
    ! constant vertical gradient of temperature, respectively
    dtvdz_thresh = 1.e-4_wp ! 0.1 K/km

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jk1,jc,dtvdz_up,dtvdz_down,p_up,p_down)

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        pres_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

      ! First, compute vertical gradients of virtual potential temperature
      DO jk = 1, nlevs_out
        DO jc = 1, nlen
          IF (jk <= bot_idx(jc,jb)) THEN
            jk1 = idx0(jc,jk,jb)
            dtvdz_up(jc,jk) = (tempv_in(jc,jk1,jb)-tempv_out(jc,jk,jb)) / &
                              (z3d_in  (jc,jk1,jb)-z3d_out  (jc,jk,jb))

            dtvdz_down(jc,jk) = (tempv_in(jc,jk1+1,jb)-tempv_out(jc,jk,jb)) / &
                                (z3d_in  (jc,jk1+1,jb)-z3d_out  (jc,jk,jb))

          ELSE ! for downward extrapolation, only dtvdz_down is needed

            dtvdz_down(jc,jk) = (tempv_out(jc,jk-1,jb)-tempv_out(jc,jk,jb)) / &
                                (z3d_out  (jc,jk-1,jb)-z3d_out  (jc,jk,jb))

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
            pres_out(jc,jk,jb) = wfac(jc,jk,jb)*p_up + (1._wp-wfac(jc,jk,jb))*p_down

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
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE pressure_intp


  !-------------
  !>
  !! SUBROUTINE z_at_plevels
  !! Computes height for a given set of pressure levels (basically the inversion 
  !! of pressure_intp). The purpose of this is to prepare diagnostic output on pressure levels.
  !!
  !! Required input fields: pressure, virtual temperature and 3D height coordinate
  !! on model levels and height levels (which have to be interpolated before 
  !! calling this routine), pressure of output levels
  !! Output: 3D height field of pressure levels
  !!
  !! Method: piecewise analytical integration of the hydrostatic equation
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE z_at_plevels(pres_ml, tempv_ml, z3d_ml, pres_zl, tempv_zl, z3d_zl,   &
                          pres_pl, z3d_pl, nblks, npromz, nlevs_ml, nlevs_zl,     &
                          nlevs_pl                                                )


    ! Input fields
    REAL(wp), INTENT(IN)  :: pres_ml  (:,:,:) ! pressure field on model levels
    REAL(wp), INTENT(IN)  :: tempv_ml (:,:,:) ! virtual temperature on model levels
    REAL(wp), INTENT(IN)  :: z3d_ml   (:,:,:) ! 3D height coordinate field on model levels
    REAL(wp), INTENT(IN)  :: pres_zl  (:,:,:) ! pressure field on height levels
    REAL(wp), INTENT(IN)  :: tempv_zl (:,:,:) ! virtual temperature on height levels
    REAL(wp), INTENT(IN)  :: z3d_zl   (:,:,:) ! 3D height coordinate field on height levels
    REAL(wp), INTENT(IN)  :: pres_pl  (:,:,:) ! pressure of output levels

    ! Comment: for z3d_zl and pres_pl, 1D fields would actually be sufficient, but having
    ! everything as 3D fields simplifies programming

    ! Output
    REAL(wp), INTENT(OUT) :: z3d_pl  (:,:,:) ! 3D height coordinate field of pressure levels

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_ml   ! Number of model levels
    INTEGER , INTENT(IN) :: nlevs_zl   ! Number of height levels
    INTEGER , INTENT(IN) :: nlevs_pl   ! Number of pressure levels

    ! LOCAL VARIABLES

    REAL(wp) :: dtvdz_thresh

    INTEGER  :: jb, jkm, jkp, jkz, jc, jkm_start, jkz_start 
    INTEGER  :: nlen
    INTEGER , DIMENSION(nproma)          :: bot_idx_ml, bot_idx_zl
    INTEGER , DIMENSION(nproma,nlevs_pl) :: idx0_ml, idx0_zl

    REAL(wp) :: z_up, z_down
    REAL(wp), DIMENSION(nproma,nlevs_pl) :: wfac_ml, wfac_zl, dtvdz

    LOGICAL  :: l_found(nproma)

!-------------------------------------------------------------------------


    ! Threshold for switching between analytical formulas for constant temperature and
    ! constant vertical gradient of temperature, respectively
    dtvdz_thresh = 1.e-4_wp ! 0.1 K/km

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jkm,jkp,jkz,jc,nlen,jkm_start,jkz_start,bot_idx_ml,bot_idx_zl, &
!$OMP            idx0_ml,idx0_zl,z_up,z_down,wfac_ml,wfac_zl,dtvdz,l_found)

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        z3d_pl(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

      ! First compute coefficients for those target levels/points for which interpolation
      ! between model-level data is possible
      jkm_start = 1
      DO jkp = 1, nlevs_pl
        l_found(:) = .FALSE.
        DO jkm = jkm_start,nlevs_ml-1
          DO jc = 1, nlen
            IF (pres_pl(jc,jkp,jb) >= pres_ml(jc,jkm,jb) .AND. &
                pres_pl(jc,jkp,jb) <  pres_ml(jc,jkm+1,jb)) THEN
              idx0_ml(jc,jkp) = jkm
              wfac_ml(jc,jkp) = (pres_pl(jc,jkp  ,jb)-pres_ml(jc,jkm,jb))/&
                                (pres_ml(jc,jkm+1,jb)-pres_ml(jc,jkm,jb))
              bot_idx_ml(jc) = jkp
              l_found(jc) = .TRUE.
            ELSE IF (pres_pl(jc,jkp,jb) > pres_ml(jc,nlevs_ml,jb)) THEN
              l_found(jc) = .TRUE.
              idx0_ml(jc,jkp) = nlevs_ml
            ENDIF
          ENDDO
          IF (ALL(l_found(1:nlen))) EXIT
        ENDDO
        jkm_start = MINVAL(idx0_ml(1:nlen,jkp))
      ENDDO

      ! Now compute coefficients for interpolation between height levels
      ! The temperatures that have already been extrapolated to height levels
      ! will afterwards be used for integrating the hydrostatic equation
      jkz_start = 1
      DO jkp = 1, nlevs_pl
        l_found(:) = .FALSE.
        DO jkz = jkz_start,nlevs_zl-1
          DO jc = 1, nlen
            IF (pres_pl(jc,jkp,jb) >= pres_zl(jc,jkz,jb) .AND. &
                pres_pl(jc,jkp,jb) <  pres_zl(jc,jkz+1,jb)) THEN
              idx0_zl(jc,jkp) = jkz
              wfac_zl(jc,jkp) = (pres_pl(jc,jkp  ,jb)-pres_zl(jc,jkz,jb))/&
                                (pres_zl(jc,jkz+1,jb)-pres_zl(jc,jkz,jb))
              bot_idx_zl(jc) = jkp
              l_found(jc) = .TRUE.
            ELSE IF (pres_pl(jc,jkp,jb) > pres_zl(jc,nlevs_zl,jb)) THEN
              l_found(jc) = .TRUE.
              idx0_zl(jc,jkp) = nlevs_zl
            ELSE IF (pres_pl(jc,jkp,jb) < pres_zl(jc,1,jb)) THEN
              l_found(jc) = .TRUE.
              idx0_zl(jc,jkp) = 1
            ENDIF
          ENDDO
          IF (ALL(l_found(1:nlen))) EXIT
        ENDDO
        jkz_start = MINVAL(idx0_zl(1:nlen,jkp))
      ENDDO

      ! Special treatment for target points lying lower than the lowest available height level
      ! (e.g. below sea level)
      DO jkp = MINVAL(bot_idx_zl(1:nlen))+1, nlevs_pl
        DO jc = 1, nlen
          IF (jkp >= bot_idx_zl(jc)+1) THEN
            ! Store extrapolation distance on wfac_zl if target point is below the 
            ! surface of the input data (note: this is a positive quantity)
            idx0_zl(jc,jkp) = nlevs_zl
            wfac_zl(jc,jkp) = pres_pl(jc,jkp,jb)-pres_zl(jc,nlevs_zl,jb)
          ENDIF
        ENDDO
      ENDDO


      ! Now, compute vertical gradients of virtual potential temperature
      DO jkp = 1, nlevs_pl
        DO jc = 1, nlen
          IF (jkp <= bot_idx_ml(jc)) THEN
            jkm = idx0_ml(jc,jkp)
            dtvdz(jc,jkp) = (tempv_ml(jc,jkm,jb)-tempv_ml(jc,jkm+1,jb)) / &
                            (z3d_ml  (jc,jkm,jb)-z3d_ml  (jc,jkm+1,jb))

          ELSE IF (jkp <= bot_idx_zl(jc)) THEN
            jkz = idx0_zl(jc,jkp)
            dtvdz(jc,jkp) = (tempv_zl(jc,jkz,jb)-tempv_zl(jc,jkz+1,jb)) / &
                            (z3d_zl  (jc,jkz,jb)-z3d_zl  (jc,jkz+1,jb))

          ELSE ! extrapolation below lowest height level required
            dtvdz(jc,jkp) = (tempv_zl(jc,nlevs_zl-1,jb)-tempv_zl(jc,nlevs_zl,jb)) / &
                            (z3d_zl  (jc,nlevs_zl-1,jb)-z3d_zl  (jc,nlevs_zl,jb))

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
            z3d_pl(jc,jkp,jb) = wfac_ml(jc,jkp)*z_up + (1._wp-wfac_ml(jc,jkp))*z_down

          ELSE IF (jkp <= bot_idx_zl(jc)) THEN

            ! interpolation based on piecewise analytical integration of the hydrostatic equation
            jkz = idx0_zl(jc,jkp)
            IF (ABS(dtvdz(jc,jkp)) > dtvdz_thresh) THEN
              z_up   = z3d_zl(jc,jkz,jb) + tempv_zl(jc,jkz,jb)*((pres_pl(jc,jkp,jb) /   &
                       pres_zl(jc,jkz,jb))**(-rd*dtvdz(jc,jkp)/grav)-1._wp)/dtvdz(jc,jkp)
              z_down = z3d_zl(jc,jkz+1,jb) + tempv_zl(jc,jkz+1,jb)*((pres_pl(jc,jkp,jb) / &
                       pres_zl(jc,jkz+1,jb))**(-rd*dtvdz(jc,jkp)/grav)-1._wp)/dtvdz(jc,jkp)

            ELSE
              z_up   = z3d_zl(jc,jkz,jb) + (rd*0.5_wp*(tempv_zl(jc,jkz,jb) +                 &
                       tempv_zl(jc,jkz+1,jb))/grav)*LOG(pres_zl(jc,jkz,jb)/pres_pl(jc,jkp,jb))
              z_down = z3d_zl(jc,jkz+1,jb) + (rd*0.5_wp*(tempv_zl(jc,jkz,jb) +                 &
                       tempv_zl(jc,jkz+1,jb))/grav)*LOG(pres_zl(jc,jkz+1,jb)/pres_pl(jc,jkp,jb))
            ENDIF

            ! Finally, apply inverse-distance weighting between top-down and bottom-up integrated value
            z3d_pl(jc,jkp,jb) = wfac_zl(jc,jkp)*z_up + (1._wp-wfac_zl(jc,jkp))*z_down

          ELSE !  extrapolation below lowest height level

            jkz = nlevs_zl
            IF (ABS(dtvdz(jc,jkp)) > dtvdz_thresh) THEN
              z_up   = z3d_zl(jc,jkz,jb) + tempv_zl(jc,jkz,jb)*((pres_pl(jc,jkp,jb) /   &
                       pres_zl(jc,jkz,jb))**(-rd*dtvdz(jc,jkp)/grav)-1._wp)/dtvdz(jc,jkp)
            ELSE
              z_up   = z3d_zl(jc,jkz,jb) + (rd*tempv_zl(jc,jkz,jb)/grav) * &
                       LOG(pres_zl(jc,jkz,jb)/pres_pl(jc,jkp,jb))
            ENDIF

            z3d_pl(jc,jkp,jb) = z_up

          ENDIF

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE z_at_plevels


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
  SUBROUTINE compute_slope(p_patch, p_int, topo_c, topo_v, slope_c)

    TYPE(t_patch),          INTENT(IN)       :: p_patch
    TYPE(t_int_state),      INTENT(IN)       :: p_int

    ! Topography data
    REAL(wp), INTENT(IN) :: topo_c(:,:)  ! topography height of mass points
    REAL(wp), INTENT(IN) :: topo_v(:,:)  ! topography height of vertices

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

    z_topo_c(:,1,:) = topo_c(:,:)
    z_topo_v(:,1,:) = topo_v(:,:)

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
                             wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_hires_corr, &
                             l_restore_sfcinv, extrapol_dist, slope          )


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

    ! Logical switch if slope-based reduction of surface inversion is to be performed
    ! (recommended when target model has a much finer resolution than source model)
    LOGICAL, INTENT(IN) :: l_hires_corr

    ! Logical switch if surface inversion is to be restored after downward extrapolation
    ! (may be set to .FALSE. when interpolating on constant height/pressure levels)
    LOGICAL, INTENT(IN) :: l_restore_sfcinv

    REAL(wp), INTENT(IN) :: extrapol_dist ! Maximum extrapolation distance using the local vertical gradient

    REAL(wp), INTENT(IN), OPTIONAL :: slope(:,:)  ! slope of mass points

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jk1, jc, nlen, jk_start, jk_start_in, jk_start_out, ik1(nproma)
    REAL(wp) :: wfac, sfcinv, vtgrad_standardatm

    REAL(wp), DIMENSION(nproma) :: temp1, temp2, vtgrad_up, zdiff_inout, &
                                   redinv1, redinv2
    LOGICAL , DIMENSION(nproma) :: l_found

    REAL(wp), DIMENSION(nproma,nlevs_in)  :: zalml_in, sfc_inv, temp_mod, g1, g2, g3
    REAL(wp), DIMENSION(nproma,nlevs_out) :: zalml_out



!-------------------------------------------------------------------------

    IF (l_hires_corr .AND. .NOT. PRESENT(slope)) CALL finish("temperature_intp:",&
      "slope correction requires slope data as input")

    ! Standard atmosphere vertical temperature gradient used for large extrapolation distances
    vtgrad_standardatm = -6.5e-3_wp

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,jk_start,jk_start_in,jk_start_out,ik1,wfac,sfcinv,&
!$OMP            temp1,temp2,vtgrad_up,zdiff_inout,redinv1,redinv2,l_found,          &
!$OMP            zalml_in,zalml_out,temp_mod,sfc_inv,g1,g2,g3)

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        temp_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

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
                                   (wfac_lin(jc,jk,jb)-extrapol_dist)*vtgrad_standardatm
            ENDIF

          ENDIF

          ! Height above lowest model level - needed for restoring the surface inversion
          zalml_out(jc,jk) = z3d_out(jc,jk,jb) - z3d_out(jc,nlevs_out,jb)

        ENDDO
      ENDDO

      DO jk = 1, nlevs_out
        IF (MINVAL(zalml_out(1:nlen,jk)) < zpbl1) THEN
          jk_start_out = jk
          EXIT
        ENDIF
      ENDDO

      ! Finally, subtract surface inversion from preliminary temperature field
      IF (l_restore_sfcinv) THEN
        jk_start = jk_start_in - 1
        DO jk = jk_start_out, nlevs_out
          l_found(:) = .FALSE.
          DO jk1 = jk_start, nlevs_in-1
            DO jc = 1, nlen
              IF (zalml_in(jc,jk1) >= zpbl1 .OR. zalml_out(jc,jk) >= zpbl1) THEN
                l_found(jc) = .TRUE.
                ik1(jc)     = jk_start
              ELSE IF (zalml_out(jc,jk) <  zalml_in(jc,jk1) .AND. &
                       zalml_out(jc,jk) >= zalml_in(jc,jk1+1)) THEN

                wfac = (zalml_out(jc,jk)-zalml_in(jc,jk1+1))/&
                       (zalml_in(jc,jk1)-zalml_in(jc,jk1+1))
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
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE temperature_intp


  !-------------
  !>
  !! SUBROUTINE uv_intp
  !! Performs vertical interpolation and extrapolation of horizontal wind components with 
  !! components with special treatment of boundary layer effects
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
    LOGICAL,  INTENT(IN), OPTIONAL :: l_restore_fricred ! restore frictional reduction of wind speed
    REAL(wp), INTENT(IN), OPTIONAL :: extrap_limit  ! multiplicative limit in case of downward extrapolation

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jk1, jc, nlen, jk_start, jk_start_in, jk_start_out, ik1(nproma)
    REAL(wp) :: wfac, fricred, mult_limit


    REAL(wp), DIMENSION(nproma) :: uv1, uv2, dudz_up
    LOGICAL , DIMENSION(nproma) :: l_found

    LOGICAL :: lrestore_fricred

    REAL(wp), DIMENSION(nproma,nlevs_in)  :: zalml_in, fric_red, uv_mod, g1, g2, g3
    REAL(wp), DIMENSION(nproma,nlevs_out) :: zalml_out, zdiff_inout, red_speed

!-------------------------------------------------------------------------

    ! Multiplicative limit for downward extrapolation; not used in high-resolution mode
    ! (wind speed in narrow valleys not resolved in the source model may be reduced
    ! more strongly)
    IF (PRESENT(extrap_limit)) THEN
      mult_limit = extrap_limit
    ELSE
      mult_limit = 0.5_wp
    ENDIF

    IF (PRESENT(l_restore_fricred)) THEN
      lrestore_fricred = l_restore_fricred
    ELSE
      lrestore_fricred = .TRUE.
    ENDIF

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,jk_start,jk_start_in,jk_start_out,ik1,wfac,fricred,&
!$OMP            uv1,uv2,dudz_up,zdiff_inout,red_speed,l_found,zalml_in,zalml_out,    &
!$OMP            uv_mod,fric_red,g1,g2,g3)

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
          IF (zalml_in(jc,jk1) < zpbl1) THEN
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

          ELSE

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
        jk_start = jk_start_in - 1
        DO jk = jk_start_out, nlevs_out
          l_found(:) = .FALSE.
          DO jk1 = jk_start, nlevs_in-1
            DO jc = 1, nlen
              IF (zalml_in(jc,jk1) >= zpbl1 .OR. zalml_out(jc,jk) >= zpbl1) THEN
                l_found(jc) = .TRUE.
                ik1(jc)     = jk_start
              ELSE IF (zalml_out(jc,jk) <  zalml_in(jc,jk1) .AND. &
                       zalml_out(jc,jk) >= zalml_in(jc,jk1+1)) THEN

                wfac = (zalml_out(jc,jk)-zalml_in(jc,jk1+1))/&
                       (zalml_in(jc,jk1)-zalml_in(jc,jk1+1))
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
!$OMP END DO
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
  SUBROUTINE qv_intp(qv_in, qv_out, z3d_in, z3d_out,                &
                     temp_in, pres_in, temp_out, pres_out,          &    
                     nblks, npromz, nlevs_in, nlevs_out,            &
                     coef1, coef2, coef3, wfac_lin,                 &
                     idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,  &
                     wfacpbl1, kpbl1, wfacpbl2, kpbl2,              &
                     lower_limit, l_restore_pbldev, opt_qc          )


    ! Specific humidity fields
    REAL(wp), INTENT(INOUT) :: qv_in (:,:,:) ! input field (INOUT because of limiter)
    REAL(wp), INTENT(OUT)   :: qv_out(:,:,:) ! output field

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
    LOGICAL , INTENT(IN) :: l_restore_pbldev ! restore PBL deviation of QV from extrapolated profile

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jk1, jc, nlen, jk_start, jk_start_in, jk_start_out, ik1(nproma)
    REAL(wp) :: wfac, pbldev, rhum, qtot


    REAL(wp), DIMENSION(nproma) :: qv1, qv2, dqvdz_up
    LOGICAL , DIMENSION(nproma) :: l_found
    LOGICAL                     :: l_check_qv_qc

    REAL(wp), DIMENSION(nproma,nlevs_in)  :: zalml_in, pbl_dev, qv_mod, g1, g2, g3, qsat_in
    REAL(wp), DIMENSION(nproma,nlevs_out) :: zalml_out, qsat_out

!-------------------------------------------------------------------------

    IF (PRESENT(opt_qc)) THEN
      l_check_qv_qc = .TRUE.
    ELSE
      l_check_qv_qc = .FALSE.
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,jk_start,jk_start_in,jk_start_out,ik1,wfac,pbldev,&
!$OMP            rhum,qtot,qv1,qv2,dqvdz_up,l_found,zalml_in,zalml_out,              &
!$OMP            qv_mod,pbl_dev,g1,g2,g3,qsat_in,qsat_out)

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        qv_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

      DO jk1 = 1, nlevs_in
        DO jc = 1, nlen

          ! saturation specific humidity of input data
          qsat_in(jc,jk1) = rdv*sat_pres_water(temp_in(jc,jk1,jb)) /    &
                            (pres_in(jc,jk1,jb)-o_m_rdv*qv_in(jc,jk1,jb))

          ! limit input data to water saturation:
          ! This is needed to remove supersaturations generated (primarily) by the interpolation
          ! from the spherical harmonics to the Gaussain grid; without this limitation, the 
          ! QV-QC-adjustment at the end of this routine would generate nonsensically large
          ! cloud water peaks
          qv_in(jc,jk1,jb) = MIN(qv_in(jc,jk1,jb),qsat_in(jc,jk1))

        ENDDO
      ENDDO


      DO jc = 1, nlen
        ! QV at height zpbl1
        qv1(jc) = wfacpbl1(jc,jb) *qv_in(jc,kpbl1(jc,jb),jb  ) + &
           (1._wp-wfacpbl1(jc,jb))*qv_in(jc,kpbl1(jc,jb)+1,jb)

        ! QV at height zpbl2
        qv2(jc) = wfacpbl2(jc,jb) *qv_in(jc,kpbl2(jc,jb),jb  ) + &
           (1._wp-wfacpbl2(jc,jb))*qv_in(jc,kpbl2(jc,jb)+1,jb)

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
            qv_mod(jc,jk1) = qv_in(jc,jk1,jb)
          ENDIF

          ! boundary-layer deviation of QV, converted to RH
          pbl_dev(jc,jk1) = (qv_mod(jc,jk1) - qv_in(jc,jk1,jb)) / qsat_in(jc,jk1)

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

          ! saturation specific humidity of output data
          qsat_out(jc,jk) = rdv*sat_pres_water(temp_out(jc,jk,jb)) /    &
                            (pres_out(jc,jk,jb)-o_m_rdv*qv_out(jc,jk,jb))

        ENDDO
      ENDDO

      DO jk = 1, nlevs_out
        IF (MINVAL(zalml_out(1:nlen,jk)) < zpbl1) THEN
          jk_start_out = jk
          EXIT
        ENDIF
      ENDDO

      ! subtract boundary-layer deviation from preliminary QV
      IF (l_restore_pbldev) THEN
        jk_start = jk_start_in - 1
        DO jk = jk_start_out, nlevs_out
          l_found(:) = .FALSE.
          DO jk1 = jk_start, nlevs_in-1
            DO jc = 1, nlen
              IF (zalml_in(jc,jk1) >= zpbl1 .OR. zalml_out(jc,jk) >= zpbl1) THEN
                l_found(jc) = .TRUE.
                ik1(jc)     = jk_start
              ELSE IF (zalml_out(jc,jk) <  zalml_in(jc,jk1) .AND. &
                       zalml_out(jc,jk) >= zalml_in(jc,jk1+1)) THEN

                wfac = (zalml_out(jc,jk)-zalml_in(jc,jk1+1))/&
                       (zalml_in(jc,jk1)-zalml_in(jc,jk1+1))
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
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE qv_intp

END MODULE mo_nh_vert_interp
