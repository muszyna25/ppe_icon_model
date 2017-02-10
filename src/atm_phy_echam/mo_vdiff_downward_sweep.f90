!>
!! @brief First half of the driver routine for turbulent mixing.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2011-04)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
#endif

MODULE mo_vdiff_downward_sweep

  USE mo_kind,               ONLY: wp
  USE mo_echam_vdiff_params, ONLY: tpfac1, itop, eps_corio
  USE mo_turbulence_diag,    ONLY: atm_exchange_coeff, sfc_exchange_coeff
  USE mo_vdiff_solver,       ONLY: nvar_vdiff, nmatrix, ih, iqv, imh, imqv, &
                                 & matrix_setup_elim, rhs_setup, rhs_elim

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: vdiff_down

CONTAINS
  !>
  !!
  !!
  SUBROUTINE vdiff_down( kproma, kbdim, klev, klevm1, klevp1, ktrac,    &! in
                       & ksfc_type, idx_wtr, idx_ice, idx_lnd,          &! in
                       & pdtime,  pcoriol,                              &! in
                       & pzf, pzh,                                      &! in
                       & pfrc,                                          &! in
                       & ptsfc_tile, pocu,      pocv,       ppsfc,      &! in
                       & pum1,       pvm1,      ptm1,       pqm1,       &! in
                       & pxlm1,      pxim1,     pxm1,       pxtm1,      &! in
                       & pmair,      pmdry,                             &! in
                       & paphm1,     papm1,                             &! in
                       & ptvm1,      paclc,     pxt_emis,   pthvvar,    &! in
                       & pxvar,      pz0m_tile,                         &! in
                       & ptkem1,                                        &! in
                       & pustar,     pwstar,    pwstar_tile,            &! inout
                       & pqsat_tile, pghpbl,                            &! out
                       & pri,        pri_tile,  pmixlen,                &! out
                       & pcfm,       pcfm_tile, pcfh,       pcfh_tile,  &! out
                       & pcfv,       pcftke,    pcfthv,                 &! out
                       & aa,         aa_btm,    bb,         bb_btm,     &! out
                       & pfactor_sfc, pcpt_tile,                        &! out
                       & pcptgz,                                        &! out
                       & pzthvvar,   pthvsig,   pztkevn,                &! out
                       & pch_tile,                                      &! out, for "nsurf_diag"
                       & pbn_tile,   pbhn_tile,                         &! out, for "nsurf_diag"
                       & pbm_tile,   pbh_tile,                          &! out, for "nsurf_diag"
                       & pcsat,                                         &! in
                       & pcair,                                         &! in
                       & paz0lh)


    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN) :: pdtime

    REAL(wp),INTENT(IN) ::          &
      & pcoriol   (kbdim)          ,&!< Coriolis parameter: 2*omega*sin(lat)
      & pzf       (kbdim,klev)     ,&!< geopotential height above sea level, full level   
      & pzh       (kbdim,klevp1)   ,&!< geopotential height above sea level, half level   
      & pfrc      (kbdim,ksfc_type),&!< area fraction of each surface type
      & ptsfc_tile(kbdim,ksfc_type),&!< surface temperature
      & pocu      (kbdim)          ,&!< eastward  velocity of ocean sfc current
      & pocv      (kbdim)          ,&!< northward velocity of ocean sfc current
      & ppsfc     (kbdim)            !< surface pressure

    REAL(wp),INTENT(IN) ::        &
      & pum1    (kbdim,klev)     ,&!< u-wind at step t-dt
      & pvm1    (kbdim,klev)     ,&!< q-wind at step t-dt
      & ptm1    (kbdim,klev)     ,&!< temperature at step t-dt
      & pqm1    (kbdim,klev)     ,&!< specific humidity at step t-dt
      & pxlm1   (kbdim,klev)     ,&!< cloud water concentration at step t-dt
      & pxim1   (kbdim,klev)     ,&!< cloud ice   concentration at step t-dt
      & pxm1    (kbdim,klev)     ,&!< cloud water + cloud ice at step t-dt
      & pxtm1   (kbdim,klev,ktrac) !< specific density of other tracers at t-dt

    REAL(wp),INTENT(IN) ::        &
      & pmair   (kbdim,klev)     ,&!<     air mass [kg/m2]
      & pmdry   (kbdim,klev)     ,&!< dra air mass [kg/m2]
      & paphm1  (kbdim,klevp1)   ,&!< half level pressure [Pa]
      & papm1   (kbdim,klev)     ,&!< full level pressure [Pa]
      & ptvm1   (kbdim,klev)     ,&!< virtual temperature
      & paclc   (kbdim,klev)     ,&!< cloud fraction
      & pxt_emis(kbdim,ktrac)      !< tracer tendency due to surface emission
                                   !< and dry deposition

    REAL(wp),INTENT(IN) ::         &
      & pthvvar  (kbdim,klev)     ,&!< variance of virtual pot. temp. at step t-dt
      & pxvar    (kbdim,klev)     ,&!< step t-dt
      & pz0m_tile(kbdim,ksfc_type)  !< roughness length at step t-dt

    REAL(wp),INTENT(IN)  :: ptkem1(kbdim,klev)    !< TKE at step t-dt

    ! Grid-box mean friction velocity.
    ! In: value at step t-2dt computed in the previous time step,
    ! used in the computation of PBL height (then mixing length);
    ! Out: computed in sfc_exchange_coeff at step t-dt.

    REAL(wp),INTENT(INOUT) :: pustar (kbdim),      &
                            & pwstar (kbdim),      &
                            & pwstar_tile(kbdim,ksfc_type)

    ! Variables with intent(out)

    REAL(wp),INTENT(INOUT) :: pqsat_tile(kbdim,ksfc_type) !< saturation specific     out
                                                          !< humidity at sfc.
                                                          !< (step t-dt)

    REAL(wp),INTENT(INOUT) :: pghpbl(kbdim)  !< geopotential height of PBL top

    REAL(wp),INTENT(INOUT) ::      &   ! out
      & pri      (kbdim,klev)     ,&!< Richardson number
      & pri_tile (kbdim,ksfc_type),&!< Richardson number
      & pmixlen  (kbdim,klev)     ,&!< mixing length
      & pcfm     (kbdim,klev)     ,&!< exchange coeff. for u, v
      & pcfm_tile(kbdim,ksfc_type),&!< exchange coeff. for u, v
      & pcfh     (kbdim,klev)     ,&!< exchange coeff. for heat and tracers
      & pcfh_tile(kbdim,ksfc_type),&!< exchange coeff. for heat and tracers
      & pcfv     (kbdim,klev)     ,&!< exchange coeff. for variance of qx
      & pcftke   (kbdim,klev)     ,&!< exchange coeff. for TKE
      & pcfthv   (kbdim,klev)       !< exchange coeff. for variance of theta_v

    ! Coefficient matrices and right-hand-side vectors.
    ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)

    REAL(wp),INTENT(INOUT) ::           &  ! out
      & aa     (kbdim,klev,3,nmatrix)  ,&!< coeff. matrices, all variables
      & aa_btm (kbdim,3,ksfc_type,imh:imqv),&!< last row of coeff. matrix of heat and moisture
      & bb     (kbdim,klev,nvar_vdiff) ,&!< r.h.s., all variables
      & bb_btm (kbdim,ksfc_type,ih:iqv)  !< last row of r.h.s. of heat and moisture

    ! Other variables to be passed on to the second part of turbulence solver

    REAL(wp),INTENT(INOUT) ::       &  ! out
      & pfactor_sfc(kbdim)         ,&!< prefactor for the exchange coeff.
      & pcpt_tile (kbdim,ksfc_type),&!< dry static energy at surface
      & pcptgz    (kbdim,klev)     ,&!< dry static energy
      & pzthvvar  (kbdim,klev)     ,&!<
      & pthvsig   (kbdim)          ,&
      & pztkevn   (kbdim,klev)       !< intermediate value of TKE

    REAL(wp), INTENT(OUT) :: pch_tile(kbdim,ksfc_type)    ! out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbn_tile(kbdim,ksfc_type)    ! out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbhn_tile(kbdim,ksfc_type)   ! out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbm_tile(kbdim,ksfc_type)    ! out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbh_tile(kbdim,ksfc_type)    ! out, for "nsurf_diag"

    REAL(wp), OPTIONAL, INTENT(IN) ::          &
      & pcsat     (kbdim)          ,&!< area fraction with wet land surface
      & pcair     (kbdim)          ,&!< area fraction with wet land surface
      & paz0lh    (kbdim)            !< surface roughness length over land for heat

    ! Local variables

    REAL(wp) :: zfcor  (kbdim)        !< absolute value of Coriolis param., > epsilon
    REAL(wp) :: zghf   (kbdim,klev)   !< geopotential height above ground, full level
    REAL(wp) :: zghh   (kbdim,klevp1) !< geopotential height above ground, full level

    REAL(wp) :: zfactor(kbdim,klev)   !< prefactor for the exchange coefficients
    REAL(wp) :: zrmairm(kbdim,klev)
    REAL(wp) :: zrmairh(kbdim,klevm1)
    REAL(wp) :: zrmdrym(kbdim,klev)

    ! _b denotes value at the bottom level (the klev-th full level)

    REAL(wp) :: ztheta_b (kbdim)  !< potential temperature
    REAL(wp) :: zthetav_b(kbdim)  !< virtual potential temperature
    REAL(wp) :: zthetal_b(kbdim)  !< liquid (and ice?) pot. temp.
    REAL(wp) :: zqsat_b  (kbdim)  !< specific humidity at saturation
    REAL(wp) :: zlh_b    (kbdim)  !< latent heat

    REAL(wp) :: zconst

    !----------------------------------------------------------------------
    ! 0. Compute useful local fields
    !----------------------------------------------------------------------

    ! absolut coriolis parameter > epsilon
    zfcor(:)               = MAX(ABS(pcoriol(:)),eps_corio)

    ! geopotential height above ground
    zghf (:,1:klev)        = pzf(:,1:klev)  -SPREAD(pzh(:,klevp1),2,klev  )
    zghh (:,1:klevp1)      = pzh(:,1:klevp1)-SPREAD(pzh(:,klevp1),2,klevp1)

    ! reciprocal layer mass
    zrmairm(1:kproma,:) = 1._wp/ pmair(1:kproma,:)
    zrmairh(1:kproma,:) = 2._wp/(pmair(1:kproma,1:klevm1)+pmair(1:kproma,2:klev))
    zrmdrym(1:kproma,:) = 1._wp/ pmdry(1:kproma,:)

    !----------------------------------------------------------------------
    ! 1. Compute various thermodynamic variables; Diagnose PBL extension;
    !    Compute exchange coefficients for half levels [1+1/2, klev-1/2];
    !    Get TKE and variance of theta_v at intermediate time step.
    !----------------------------------------------------------------------

    CALL atm_exchange_coeff( kproma, kbdim, klev, klevm1, klevp1,     &! in
                           & pdtime, zfcor,                           &! in
                           & zghf, zghh,                              &! in
                           & pum1, pvm1, ptm1, ptvm1,                 &! in
                           & pqm1, pxm1,                              &! in
                           & papm1, paphm1, paclc, pustar,            &! in
                           & pthvvar, ptkem1,                         &! in
                           & pcptgz, pghpbl,                          &! out
                           & pzthvvar(:,1:klevm1),pztkevn(:,1:klevm1),&! out
                           & pcfm    (:,1:klevm1), pcfh  (:,1:klevm1),&! out
                           & pcfv    (:,1:klevm1), pcftke(:,1:klevm1),&! out
                           & pcfthv  (:,1:klevm1),zfactor(:,1:klevm1),&! out
                           & ztheta_b, zthetav_b, zthetal_b,          &! out, for "sfc_exchange_coeff"
                           & zqsat_b,  zlh_b,                         &! out, for "sfc_exchange_coeff"
                           & pri(:,1:klevm1), pmixlen(:,1:klevm1)     )! out, for output

    pmixlen(:,klev) = -999._wp ! dummy value, as not defined for jk=klev in atm_exchange_coeff

    !-----------------------------------------------------------------------
    ! 2. Compute exchange coefficients at the air-sea/ice/land interface.
    !    Get boundary condition for TKE and variance of theta_v.
    !-----------------------------------------------------------------------

    CALL sfc_exchange_coeff( kproma, kbdim, ksfc_type,              &! in
                           & idx_wtr, idx_ice, idx_lnd,             &! in
                           & pz0m_tile(:,:),  ptsfc_tile(:,:),      &! in
                           & pfrc(:,:),       pghpbl(:),            &! in
                           & pocu(:),         pocv(:),   ppsfc(:),  &! in
                           & zfcor,           zghf(:,klev),         &! in
                           & pum1(:,klev),    pvm1  (:,klev),       &! in
                           & ptm1(:,klev),                          &! in
                           & pqm1(:,klev),    pxm1  (:,klev),       &! in
                           & zqsat_b  (:),    zlh_b    (:),         &! in
                           & ztheta_b (:),    zthetav_b(:),         &! in
                           & zthetal_b(:),    paclc (:,klev),       &! in
                           & ptkem1(:,klev),  pzthvvar(:,klevm1),   &! in
                           & pthvsig(:),                            &! inout
                           & pwstar(:),       pwstar_tile(:,:),     &! inout
                           & pqsat_tile(:,:), pcpt_tile(:,:),       &! out
                           & pri    (:,klev), pri_tile(:,:),        &! out
                           & pcfm   (:,klev), pcfm_tile(:,:),       &! out
                           & pcfh   (:,klev), pcfh_tile(:,:),       &! out
                           & pcfv   (:,klev),                       &! out
                           & pcftke (:,klev), pcfthv  (:,klev),     &! out
                           & zfactor(:,klev),                       &! out
                           & pztkevn(:,klev), pzthvvar(:,klev),     &! out
                           & pustar(:),                             &! out, for "atm_exchange_coeff" at next time step
                           & pch_tile(:,:),                         &! out, for "nsurf_diag"
                           & pbn_tile(:,:),   pbhn_tile(:,:),       &! out, for "nsurf_diag"
                           & pbm_tile(:,:),   pbh_tile(:,:),        &! out, for "nsurf_diag"
                           & paz0lh(:),                             &! in, optional
                           & pcsat(:),                              &! in, optional
                           & pcair(:))                               ! in, optional


    !-----------------------------------------------------------------------
    ! 3. Set up coefficient matrix of the tri-diagonal system, then perform
    !    Gauss elimination for it. The matrix is built from
    !    - the exchange coefficients;
    !    - the prefactor "zfactor" and some additional constants ("zconst")
    !      which are determined by the spatial and temporal discretization
    !      employed for vertical diffusion;
    !    - the assumption about upper and lower boundaries, especially
    !      whether there is turbulent flux at the lower boundary for each
    !      quantity subject to turbulent mixing.
    !-----------------------------------------------------------------------

    zconst = tpfac1*pdtime
    zfactor(1:kproma,1:klevm1) = zfactor(1:kproma,1:klevm1)*zconst
    zfactor(1:kproma,  klev)   = zfactor(1:kproma,  klev)  *zconst

    CALL matrix_setup_elim( kproma, kbdim, klev, klevm1, ksfc_type, itop, &! in
                          & pcfm     (:,:),   pcfh  (:,1:klevm1),         &! in
                          & pcfh_tile(:,:),   pcfv  (:,:),                &! in
                          & pcftke   (:,:),   pcfthv(:,:),                &! in
                          & zfactor  (:,:),                               &! in
                          & zrmairm, zrmairh, zrmdrym,                    &! in
                          & aa, aa_btm                                    )! out

    ! Save for output, to be used in "update_surface"
    pfactor_sfc(1:kproma) = zfactor(1:kproma,klev)

    !-----------------------------------------------------------------------
    ! 4. Set up right-hand side of the tri-diagonal system and perform
    !    Gauss elimination. Factors that determine the r.h.s. include
    !    - time stepping scheme used for vertical diffusion
    !    - whether there is any other process (e.g., tracer emission)
    !      solved together with vertical diffusion.
    !-----------------------------------------------------------------------

    CALL rhs_setup( kproma, kbdim, itop, klev, klevm1,    &! in
                  & ksfc_type, ktrac, pdtime,             &! in
                  & pum1, pvm1, pcptgz, pqm1,             &! in
                  & pxlm1, pxim1, pxvar, pxtm1, pxt_emis, &! in
                  & zrmdrym, pztkevn, pzthvvar, aa,       &! in
                  & bb, bb_btm                            )! out

    CALL rhs_elim ( kproma, kbdim, itop, klev, klevm1, &! in
                  & aa, bb                             )! in, inout

  END SUBROUTINE vdiff_down
  !-------------

END MODULE mo_vdiff_downward_sweep
