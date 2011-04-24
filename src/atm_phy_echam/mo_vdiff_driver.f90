!>
!! @brief Driver routine for turbulent mixing.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2010-09)
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
#ifdef __xlC__
@PROCESS HOT
#endif

MODULE mo_vdiff_driver

  USE mo_kind,               ONLY: wp
  USE mo_turbulence_diag,    ONLY: atm_exchange_coeff, sfc_exchange_coeff
  USE mo_vdiff_solver,       ONLY: nvar_vdiff, nmatrix, ih, iqv,          &
                                 & matrix_setup_elim, rhs_setup, rhs_elim,&
                                 & sfc_solve, rhs_bksub, vdiff_tendencies
#ifdef __ICON__
  USE mo_physical_constants, ONLY: grav, rd
  USE mo_echam_vdiff_params, ONLY: tpfac1, tpfac2, itop,        &
                                 & lsfc_mom_flux, lsfc_heat_flux
#else
  USE mo_constants, ONLY: grav=>g, rd
  USE mo_physc2,    ONLY: tpfac1, tpfac2, itop
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: vdiff

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
  !!
  SUBROUTINE vdiff ( kproma, kbdim, klev, klevm1, klevp1, ktrac,       &! in
                     ksfc_type, idx_wtr, idx_ice, idx_lnd, idx_gbm,    &! in
                     pdtime, pstep_len,                                &! in
                     pcoriol,    pfrc,        ptsfc,                   &! in
                     pocu,       pocv,        ppsfc,                   &! in
                     pum1,       pvm1,        ptm1,        pqm1,       &! in
                     pxlm1,      pxim1,       pxm1,        pxtm1,      &! in
                     paphm1,     papm1,       pdelpm1,     pgeom1,     &! in
#ifdef __ICON__
                     ptvm1,      paclc,       ptkem1,      pxt_emis,   &! in
#else
                     ptvm1,      paclc,    ptkem1, ptkem0, pxt_emis,   &! in/inout
#endif
                     pxvar,      pthvvar,                              &! inout
                     pustar,     pz0m,        pkedisp,                 &! inout
                     pute,       pvte,        ptte,        pqte,       &! inout
                     pxlte,      pxite,       pxtte,                   &! inout
                     pute_vdf,   pvte_vdf,    ptte_vdf,                &! out
                     pqte_vdf,   pxlte_vdf,   pxite_vdf,   pxtte_vdf,  &! out
                     pqsat_tile,                                       &! out
                     pxvarprod,  pvmixtau,    pqv_mflux_sfc,           &! out
                     pthvsig,    ptke,        ihpbl,       pghpbl,     &! out
                     pri,        pmixlen,                              &! out
                     pcfm,       pcfm_tile,                            &! out
                     pcfh,       pcfh_tile,                            &! out
                     pcfv,       pcftke,      pcfthv                   )! out

  INTEGER, INTENT(IN)  :: kproma, kbdim, klev, klevm1, klevp1, ktrac
  INTEGER, INTENT(IN)  :: ksfc_type, idx_wtr, idx_ice, idx_lnd, idx_gbm
  REAL(wp),INTENT(IN)  :: pdtime, pstep_len

  REAL(wp),INTENT(IN)  :: pcoriol   (kbdim)           !< Coriolis parameter: 2*omega*sin(lat)
  REAL(wp),INTENT(IN)  :: pfrc      (kbdim,ksfc_type) !< area fraction of each surface type
  REAL(wp),INTENT(IN)  :: ptsfc     (kbdim,ksfc_type) !< surface temperature
  REAL(wp),INTENT(IN)  :: pocu      (kbdim)           !< eastward  velocity of ocean surface current
  REAL(wp),INTENT(IN)  :: pocv      (kbdim)           !< northward velocity of ocean surface current
  REAL(wp),INTENT(IN)  :: ppsfc     (kbdim)           !<  surface pressure

  REAL(wp),INTENT(IN)  :: pum1    (kbdim,klev)  !< u-wind at step t-dt
  REAL(wp),INTENT(IN)  :: pvm1    (kbdim,klev)  !< q-wind at step t-dt
  REAL(wp),INTENT(IN)  :: ptm1    (kbdim,klev)  !< temperature at step t-dt
  REAL(wp),INTENT(IN)  :: pqm1    (kbdim,klev)  !< specific humidity at step t-dt
  REAL(wp),INTENT(IN)  :: pxlm1   (kbdim,klev)  !< cloud water concentration at step t-dt
  REAL(wp),INTENT(IN)  :: pxim1   (kbdim,klev)  !< cloud ice   concentration at step t-dt
  REAL(wp),INTENT(IN)  :: pxm1    (kbdim,klev)  !< cloud water + cloud ice at step t-dt
  REAL(wp),INTENT(IN)  :: pxtm1   (kbdim,klev,ktrac) !< specific density of other tracers at step t-dt

  REAL(wp),INTENT(IN)  :: paphm1 (kbdim,klevp1) !< half level pressure [Pa]
  REAL(wp),INTENT(IN)  :: papm1  (kbdim,klev)   !< full level pressure [Pa]
  REAL(wp),INTENT(IN)  :: pdelpm1(kbdim,klev)   !< layer thickness [Pa]
  REAL(wp),INTENT(IN)  :: pgeom1 (kbdim,klev)   !< geopotential above ground
  REAL(wp),INTENT(IN)  :: ptvm1  (kbdim,klev)   !< virtual temperature
  REAL(wp),INTENT(IN)  :: paclc  (kbdim,klev)   !< cloud fraction
#ifdef __ICON__
  REAL(wp),INTENT(IN)  :: ptkem1(kbdim,klev)    !< TKE at time step t-dt 
#else
  REAL(wp),INTENT(INOUT) :: ptkem1(kbdim,klev)
  REAL(wp),INTENT(INOUT) :: ptkem0(kbdim,klev)
#endif

  REAL(wp),INTENT(IN)  :: pxt_emis(kbdim,ktrac) !< tracer tendency due to surface emission
                                                !< and dry deposition

  REAL(wp),INTENT(INOUT) :: pxvar  (kbdim,klev) !<distribution width (b-a) 
                                                !< in: step t-dt, out: modified due to vertical diffusion
  REAL(wp),INTENT(INOUT) :: pthvvar(kbdim,klev) !< variance of virtual potential temperature
                                                !< in: step t-dt, out: modified due to vertical diffusion,
                                                !< shear production and dissipation
  REAL(wp),INTENT(INOUT) :: pustar(kbdim)  !< grid-box mean friction velocity
                                           !< in: value at step t-2dt computed in the previous time step,
                                           !< used in the computation of PBL height (then
                                           !< mixing length; out: computed in sfc_exchange_coeff
                                           !< at step t-dt.

  REAL(wp),INTENT(INOUT) :: pz0m(kbdim,idx_gbm:ksfc_type) !< roughness length
                                                  !< in: values over each surface type computed
                                                  !< during the previous time step;
                                                  !< out: new values computed from the temporally
                                                  !< weighted "hat" value of u and v and
                                                  !< the step t-dt value of some other variables.

  REAL(wp),INTENT(INOUT) :: pkedisp(kbdim) !< temporally and vertically integrated
                                           !< dissipation of kinetic energy

  REAL(wp),INTENT(INOUT) :: pute (kbdim,klev)
  REAL(wp),INTENT(INOUT) :: pvte (kbdim,klev)
  REAL(wp),INTENT(INOUT) :: ptte (kbdim,klev)
  REAL(wp),INTENT(INOUT) :: pqte (kbdim,klev)
  REAL(wp),INTENT(INOUT) :: pxlte(kbdim,klev)
  REAL(wp),INTENT(INOUT) :: pxite(kbdim,klev)
  REAL(wp),INTENT(INOUT) :: pxtte(kbdim,klev,ktrac)

  REAL(wp),INTENT(OUT) :: pute_vdf (kbdim,klev)
  REAL(wp),INTENT(OUT) :: pvte_vdf (kbdim,klev)
  REAL(wp),INTENT(OUT) :: ptte_vdf (kbdim,klev)
  REAL(wp),INTENT(OUT) :: pqte_vdf (kbdim,klev)
  REAL(wp),INTENT(OUT) :: pxlte_vdf(kbdim,klev)
  REAL(wp),INTENT(OUT) :: pxite_vdf(kbdim,klev)
  REAL(wp),INTENT(OUT) :: pxtte_vdf(kbdim,klev,ktrac)

  REAL(wp),INTENT(OUT) :: pqsat_tile(kbdim,ksfc_type)!< surface specific humidity at saturation

  REAL(wp),INTENT(OUT) :: pxvarprod    (kbdim,klev) !< shear production of the variance of total water
  REAL(wp),INTENT(OUT) :: pvmixtau     (kbdim,klev) !< vertical mixing time scale
  REAL(wp),INTENT(OUT) :: pqv_mflux_sfc(kbdim)      !< surface mass flux of water vapour
  REAL(wp),INTENT(OUT) :: pthvsig      (kbdim)      !< sqrt( variance of theta_v )

  REAL(wp),INTENT(OUT) :: ptke  (kbdim,klev)
  REAL(wp),INTENT(OUT) :: pghpbl(kbdim)
  INTEGER, INTENT(OUT) :: ihpbl (kbdim)

  REAL(wp),INTENT(OUT) :: pri      (kbdim,klev) !< Richardson number
  REAL(wp),INTENT(OUT) :: pmixlen  (kbdim,klev) !< mixing length
  REAL(wp),INTENT(OUT) :: pcfm     (kbdim,klev) !< exchange coeff. for u, v
  REAL(wp),INTENT(OUT) :: pcfm_tile(kbdim,ksfc_type) !< exchange coeff.
  REAL(wp),INTENT(OUT) :: pcfh     (kbdim,klev) !< exchange coeff. for heat and tracers
  REAL(wp),INTENT(OUT) :: pcfh_tile(kbdim,ksfc_type) !< exchange coeff. for heat and tracers
  REAL(wp),INTENT(OUT) :: pcfv    (kbdim,klev) !< exchange coeff. for variance of qx
  REAL(wp),INTENT(OUT) :: pcftke  (kbdim,klev) !< exchange coeff. for TKE
  REAL(wp),INTENT(OUT) :: pcfthv  (kbdim,klev) !< exchange coeff. for variance of theta_v

  ! Local variables

  REAL(wp) :: zcpt_tile(kbdim,ksfc_type)
  REAL(wp) :: zthvvar (kbdim,klev)
  REAL(wp) :: ztkevn  (kbdim,klev)
  REAL(wp) :: zqshear (kbdim,klev)
  REAL(wp) :: zprfac  (kbdim,klev) !< prefactor for the exchange coefficients
  REAL(wp) :: zrhoh   (kbdim,klev)  !< air density at half levels
  REAL(wp) :: zcptgz  (kbdim,klev)
  REAL(wp) :: zrdpm   (kbdim,klev)
  REAL(wp) :: zrdph   (kbdim,klevm1)

  ! _b denotes value at the bottom level (the klev-th full level)

  REAL(wp) :: ztheta_b (kbdim)  !< potential temperature
  REAL(wp) :: zthetav_b(kbdim)  !< virtual potential temperature
  REAL(wp) :: zthetal_b(kbdim)  !< liquid (and ice?) pot. temp.
  REAL(wp) :: zqsat_b  (kbdim)  !< specific humidity at saturation
  REAL(wp) :: zlh_b    (kbdim)  !< latent heat

  REAL(wp) :: zconst

  ! Coefficient matrices and right-hand-side vectors.
  ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)
 
  REAL(wp) :: aa    (kbdim,klev,3,nmatrix)   !< (nproma,nlev,3,nmatrix)
  REAL(wp) :: aa_btm(kbdim,3,ksfc_type)      !< (nproma,3,ksfc_type), for heat and moisture
  REAL(wp) :: bb    (kbdim,klev,nvar_vdiff)  !< (nproma,nlev,nvar_vdiff)
  REAL(wp) :: bb_btm(kbdim,ksfc_type,ih:iqv) !< (nproma,ksfc_type,ih:iqv), for heat and moisture
  
  !----------------------------------------------------------------------
  ! 0. Reciprocal of layer thickness. It will be used repeatedly.
  !----------------------------------------------------------------------

  zrdpm(1:kproma,:) = 1._wp/pdelpm1(1:kproma,:)
  zrdph(1:kproma,:) = 1._wp/(papm1(1:kproma,2:klev)-papm1(1:kproma,1:klev-1))

  !----------------------------------------------------------------------
  ! 1. Compute various thermodynamic variables; Diagnose PBL extension;
  !    Compute exchange coefficients for half levels [1+1/2, klev-1/2];
  !    Get TKE and variance of theta_v at intermediate time step.
  !----------------------------------------------------------------------

  CALL atm_exchange_coeff( kproma, kbdim, klev, klevm1, klevp1,     &! in
                         & pstep_len, pcoriol,                      &! in
                         & pum1, pvm1, ptm1, ptvm1, pgeom1,         &! in
                         & pqm1, pxm1,                              &! in
                         & papm1, paphm1, paclc, pustar,            &! in
#ifdef __ICON__
                         & pthvvar, ptkem1,                         &! in
#else
                         & pthvvar, ptkem1, ptkem0,                 &! in, inout, inout
#endif
                         & zcptgz (:,1:klev),   ihpbl(:), pghpbl(:),&! out
                         & zqshear(:,1:klevm1),                     &! out
                         & zthvvar(:,1:klevm1), ztkevn(:,1:klevm1), &! out
                         & pcfm   (:,1:klevm1), pcfh  (:,1:klevm1), &! out
                         & pcfv   (:,1:klevm1), pcftke(:,1:klevm1), &! out
                         & pcfthv (:,1:klevm1), zprfac(:,1:klevm1), &! out
                         & zrhoh  (:,1:klevm1),                     &! out, for "vdiff_tendencies"
                         & ztheta_b(:), zthetav_b(:), zthetal_b(:), &! out, for "sfc_exchange_coeff"
                         & zqsat_b(:),  zlh_b(:),                   &! out, for "sfc_exchange_coeff"
                         & pri(:,1:klevm1), pmixlen(:,1:klevm1)     )! out, for output

  !-----------------------------------------------------------------------
  ! 2. Compute exchange coefficients at the air-sea/ice/land interface.
  !    Get boundary condition for TKE and variance of theta_v.
  !-----------------------------------------------------------------------

  CALL sfc_exchange_coeff( kproma, kbdim, ksfc_type,              &! in
                         & idx_wtr, idx_ice, idx_lnd,             &! in
                         & lsfc_mom_flux, lsfc_heat_flux,         &! in
                         & pz0m(:,1:), ptsfc(:,:),                &! in
                         & pfrc(:,:), pghpbl(:),                  &! in
                         & pocu(:),         pocv(:),   ppsfc(:),  &! in
                         & pum1(:,klev),    pvm1  (:,klev),       &! in
                         & ptm1(:,klev),    pgeom1(:,klev),       &! in
                         & pqm1(:,klev),    pxm1  (:,klev),       &! in
                         & zqsat_b  (:),    zlh_b    (:),         &! in
                         & ztheta_b (:),    zthetav_b(:),         &! in
                         & zthetal_b(:),    paclc (:,klev),       &! in
                         & zthvvar(:,klevm1),                     &! in
#ifdef __ICON__
#else
                         & ptkem1(:,klev),  ptkem0(:,klev),       &! inout
#endif
                         & pqsat_tile(:,:), zcpt_tile(:,:),       &! out
                         & pri    (:,klev),                       &! out
                         & pcfm   (:,klev), pcfm_tile(:,:),       &! out
                         & pcfh   (:,klev), pcfh_tile(:,:),       &! out
                         & pcfv   (:,klev),                       &! out
                         & pcftke (:,klev), pcfthv (:,klev),      &! out
                         & zprfac (:,klev), zrhoh  (:,klev),      &! out
                         & ztkevn (:,klev), zthvvar(:,klev),      &! out
                         & zqshear(:,klev),                       &! out, for "vdiff_tendencies"
                         &  pustar(:)                             )! out, for "atm_exchange_coeff"
                                                                   ! at next time step

  !-----------------------------------------------------------------------
  ! 3. Set up coefficient matrix of the tri-diagonal system, then perform
  !    Gaussian elimination for the matrix. The matrix is built from
  !    - the exchange coefficients;
  !    - the prefactor "zprfac" and some additional constants ("zconst")
  !      which are determined by the spatial and temporal discretization
  !      employed for vertical diffusion;
  !    - the assumption about upper and lower boundaries, especially
  !      whether there is turbulent flux at the lower boundary for each
  !      quantity subject to turbulent mixing.
  !-----------------------------------------------------------------------

  zconst = tpfac1*pstep_len*grav*grav
  zprfac(1:kproma,1:klevm1) = zprfac(1:kproma,1:klevm1)*zconst

  zconst = tpfac1*pstep_len*grav/rd
  zprfac(1:kproma,  klev)   = zprfac(1:kproma,  klev)  *zconst

  CALL matrix_setup_elim( kproma, kbdim, klev, klevm1, ksfc_type, itop, &! in
                        & pcfm     (:,:),   pcfh  (:,1:klevm1),         &! in
                        & pcfh_tile(:,:),   pcfv  (:,:),                &! in
                        & pcftke   (:,:),   pcfthv(:,:),                &! in
                        & zprfac   (:,:),   zrdpm, zrdph,               &! in
                        & aa, aa_btm                                    )! out

  !-----------------------------------------------------------------------
  ! 4. Set up right-hand side of the tri-diagonal system and perform 
  !    Gaussian elimination. Factors that determine the r.h.s. include
  !    - time stepping scheme used for vertical diffusion
  !    - whether there is any other process (e.g., tracer emission)
  !      solved together with vertical diffusion.
  !-----------------------------------------------------------------------

  CALL rhs_setup( kproma, kbdim, itop, klev, klevm1,    &! in
                & ksfc_type, ktrac, tpfac2, pstep_len,  &! in
                & pum1, pvm1, zcptgz, pqm1,             &! in
                & pxlm1, pxim1, pxvar, pxtm1, pxt_emis, &! in
                & zrdpm, ztkevn, zthvvar, aa,           &! in
                & bb, bb_btm                            )! out

  CALL rhs_elim ( kproma, kbdim, itop, klev, klevm1, &! in
                & aa, bb                             )! in, inout

  !-----------------------------------------------------------------------
  ! 5. Handle fractional surfaces (tiles) and obtain solution of 
  !    static energy and water vapor concentration in the lowest model 
  !    layer (the klev-th full level). 
  !-----------------------------------------------------------------------

  CALL sfc_solve( kproma, kbdim, klev, klevm1,    &! in
                & ksfc_type, idx_wtr, idx_ice,    &! in
                & lsfc_heat_flux, tpfac2, pfrc,   &! in
                & pocu, pocv, zcpt_tile,          &! in
                & pqsat_tile,                     &! in
                & pcfh_tile, zprfac(:,klev),      &! in
                & aa, aa_btm,                     &! in
                & bb, bb_btm                      )! inout

  !-----------------------------------------------------------------------
  ! 6. Obtain solution of the tri-diagonal system by back-substitution. 
  !    Then compute tendencies and diagnose moisture flux etc.
  !-----------------------------------------------------------------------
  CALL rhs_bksub( kproma, kbdim, itop, klev, aa, bb ) ! in,...,in, inout

  CALL vdiff_tendencies( kproma, kbdim, itop, klev, klevm1, klevp1,   &! in
                       & ktrac, ksfc_type, idx_lnd, idx_wtr,          &! in
                       & pdtime, pstep_len,                           &! in
                       & pum1, pvm1, ptm1, pqm1, pxlm1, pxim1,        &! in
                       & pxtm1, pgeom1, pdelpm1, zcptgz,              &! in
#ifdef __ICON__
                       & ptkem1, ztkevn, zthvvar, zrhoh,              &! in
#else
                       & ptkem1, ptkem0, ztkevn, zthvvar, zrhoh,      &! in
#endif
                       & zqshear, ihpbl, pcfh_tile, pqsat_tile,       &! in
                       & pcfm_tile, pfrc, bb,                         &! in
                       & pkedisp(:),                                  &! inout ("pvdis" in echam)
                       & pxvar(:,:), pz0m(:,1:ksfc_type),             &! inout
                       & pute, pvte, ptte, pqte, pxlte, pxite, pxtte, &! inout
                       & pute_vdf, pvte_vdf, ptte_vdf, pqte_vdf,      &! out
                       & pxlte_vdf, pxite_vdf, pxtte_vdf,             &! out
                       & pxvarprod,                                   &! out ("pvdiffp" in echam)
                       & pz0m(:,idx_gbm),                             &! out
                       & ptke, pthvvar, pthvsig, pvmixtau,            &! out
                       & pqv_mflux_sfc                                )! out ("pqhfla" in echam)

  ! Note: computation of additional diagnostics, e.g., surface sensible heat flux,
  !       wind stress, 10m wind, 2m temperature etc., has not been implemented yet.

  END SUBROUTINE vdiff
  !-------------

END MODULE mo_vdiff_driver

