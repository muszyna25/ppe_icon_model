!>
!! @brief First half of the driver routine for turbulent mixing. 
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2011-04)
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

MODULE mo_vdiff_downward_sweep

  USE mo_kind,               ONLY: wp
  USE mo_turbulence_diag,    ONLY: atm_exchange_coeff, sfc_exchange_coeff
  USE mo_vdiff_solver,       ONLY: nvar_vdiff, nmatrix, ih, iqv,         &
                                 & matrix_setup_elim, rhs_setup, rhs_elim
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
  PUBLIC :: vdiff_down

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
  !!
  SUBROUTINE vdiff_down( kproma, kbdim, klev, klevm1, klevp1, ktrac,    &! in
                       & ksfc_type, idx_wtr, idx_ice, idx_lnd,          &! in
                       & pstep_len,  pcoriol,   pfrc,                   &! in
                       & ptsfc_tile, pocu,      pocv,       ppsfc,      &! in
                       & pum1,       pvm1,      ptm1,       pqm1,       &! in
                       & pxlm1,      pxim1,     pxm1,       pxtm1,      &! in
                       & paphm1,     papm1,     pdelpm1,    pgeom1,     &! in
                       & ptvm1,      paclc,     pxt_emis,   pthvvar,    &! in
                       & pxvar,      pz0m_tile,                         &! in
#ifdef __ICON__
                       & ptkem1,                                        &! in
#else
                       & ptkem1,     ptkem0,                            &! inout
#endif
                       & pustar,                                        &! inout
                       & pqsat_tile, ihpbl,     pghpbl,                 &! out
                       & pri,        pmixlen,                           &! out
                       & pcfm,       pcfm_tile, pcfh,       pcfh_tile,  &! out
                       & pcfv,       pcftke,    pcfthv,                 &! out
                       & aa,         aa_btm,    bb,         bb_btm,     &! out
                       & pprfac_sfc, pcpt_tile,                         &! out
                       & pcptgz,     prhoh,     pqshear,                &! out
                       & pzthvvar,   pztkevn                            )! out

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN) :: pstep_len

    REAL(wp),INTENT(IN) ::          &
      & pcoriol   (kbdim)          ,&!< Coriolis parameter: 2*omega*sin(lat)
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
      & paphm1  (kbdim,klevp1)   ,&!< half level pressure [Pa]
      & papm1   (kbdim,klev)     ,&!< full level pressure [Pa]
      & pdelpm1 (kbdim,klev)     ,&!< layer thickness [Pa]
      & pgeom1  (kbdim,klev)     ,&!< geopotential above ground
      & ptvm1   (kbdim,klev)     ,&!< virtual temperature
      & paclc   (kbdim,klev)     ,&!< cloud fraction
      & pxt_emis(kbdim,ktrac)      !< tracer tendency due to surface emission
                                   !< and dry deposition

    REAL(wp),INTENT(IN) ::         &
      & pthvvar  (kbdim,klev)     ,&!< variance of virtual pot. temp. at step t-dt
      & pxvar    (kbdim,klev)     ,&!< step t-dt
      & pz0m_tile(kbdim,ksfc_type)  !< roughness length at step t-dt

#ifdef __ICON__
    REAL(wp),INTENT(IN)  :: ptkem1(kbdim,klev)    !< TKE at step t-dt 
#else
    REAL(wp),INTENT(INOUT) :: ptkem1(kbdim,klev)  !< TKE at step t-dt
    REAL(wp),INTENT(INOUT) :: ptkem0(kbdim,klev)  !< TKE at step t
#endif

    ! Grid-box mean friction velocity.
    ! In: value at step t-2dt computed in the previous time step,
    ! used in the computation of PBL height (then mixing length);
    ! Out: computed in sfc_exchange_coeff at step t-dt.

    REAL(wp),INTENT(INOUT) :: pustar (kbdim)

    ! Variables with intent(out)

    REAL(wp),INTENT(OUT) :: pqsat_tile(kbdim,ksfc_type) !< saturation specific 
                                                        !< humidity at sfc.
                                                        !< (step t-dt)

    INTEGER, INTENT(OUT) :: ihpbl (kbdim)  !< PBL height given as level index
    REAL(wp),INTENT(OUT) :: pghpbl(kbdim)  !< geopotential height of PBL top

    REAL(wp),INTENT(OUT) ::        &
      & pri      (kbdim,klev)     ,&!< Richardson number
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
 
    REAL(wp),INTENT(OUT) ::             &
      & aa     (kbdim,klev,3,nmatrix)  ,&!< coeff. matrices, all variables
      & aa_btm (kbdim,3,ksfc_type)     ,&!< last row of coeff. matrix of heat and moisture
      & bb     (kbdim,klev,nvar_vdiff) ,&!< r.h.s., all variables
      & bb_btm (kbdim,ksfc_type,ih:iqv)  !< last row of r.h.s. of heat and moisture

    ! Other variables to be passed on to the second part of turbulence solver

    REAL(wp),INTENT(OUT) ::         &
      & pprfac_sfc(kbdim)          ,&!< prefactor for the exchange coeff.
      & pcpt_tile (kbdim,ksfc_type),&!< dry static energy at surface
      & pcptgz    (kbdim,klev)     ,&!< dry static energy
      & prhoh     (kbdim,klev)     ,&!< air density at half levels
      & pqshear   (kbdim,klev)     ,&!<
      & pzthvvar  (kbdim,klev)     ,&!<
      & pztkevn   (kbdim,klev)       !< intermediate value of TKE

    ! Local variables

    REAL(wp) :: zprfac (kbdim,klev)   !< prefactor for the exchange coefficients
    REAL(wp) :: zrdpm  (kbdim,klev)
    REAL(wp) :: zrdph  (kbdim,klevm1)

    ! _b denotes value at the bottom level (the klev-th full level)

    REAL(wp) :: ztheta_b (kbdim)  !< potential temperature
    REAL(wp) :: zthetav_b(kbdim)  !< virtual potential temperature
    REAL(wp) :: zthetal_b(kbdim)  !< liquid (and ice?) pot. temp.
    REAL(wp) :: zqsat_b  (kbdim)  !< specific humidity at saturation
    REAL(wp) :: zlh_b    (kbdim)  !< latent heat

    REAL(wp) :: zconst

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
                           & pcptgz (:,1:klev),   ihpbl(:), pghpbl(:),&! out
                           & pqshear(:,1:klevm1),                     &! out
                           & pzthvvar(:,1:klevm1),pztkevn(:,1:klevm1),&! out
                           & pcfm   (:,1:klevm1), pcfh  (:,1:klevm1), &! out
                           & pcfv   (:,1:klevm1), pcftke(:,1:klevm1), &! out
                           & pcfthv (:,1:klevm1), zprfac(:,1:klevm1), &! out
                           & prhoh  (:,1:klevm1),                     &! out, for "vdiff_tendencies"
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
                           & pz0m_tile(:,:),  ptsfc_tile(:,:),      &! in
                           & pfrc(:,:),       pghpbl(:),            &! in
                           & pocu(:),         pocv(:),   ppsfc(:),  &! in
                           & pum1(:,klev),    pvm1  (:,klev),       &! in
                           & ptm1(:,klev),    pgeom1(:,klev),       &! in
                           & pqm1(:,klev),    pxm1  (:,klev),       &! in
                           & zqsat_b  (:),    zlh_b    (:),         &! in
                           & ztheta_b (:),    zthetav_b(:),         &! in
                           & zthetal_b(:),    paclc (:,klev),       &! in
                           & pzthvvar(:,klevm1),                    &! in
#ifdef __ICON__
#else
                           & ptkem1(:,klev),  ptkem0(:,klev),       &! inout
#endif
                           & pqsat_tile(:,:), pcpt_tile(:,:),       &! out
                           & pri    (:,klev),                       &! out
                           & pcfm   (:,klev), pcfm_tile(:,:),       &! out
                           & pcfh   (:,klev), pcfh_tile(:,:),       &! out
                           & pcfv   (:,klev),                       &! out
                           & pcftke (:,klev), pcfthv  (:,klev),     &! out
                           & zprfac (:,klev), prhoh   (:,klev),     &! out
                           & pztkevn(:,klev), pzthvvar(:,klev),     &! out
                           & pqshear(:,klev),                       &! out, for "vdiff_tendencies"
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
 
    ! Save for output, to be used in "surface_solve"
    pprfac_sfc(1:kproma) = zprfac(1:kproma,klev)

    !-----------------------------------------------------------------------
    ! 4. Set up right-hand side of the tri-diagonal system and perform 
    !    Gaussian elimination. Factors that determine the r.h.s. include
    !    - time stepping scheme used for vertical diffusion
    !    - whether there is any other process (e.g., tracer emission)
    !      solved together with vertical diffusion.
    !-----------------------------------------------------------------------
  
    CALL rhs_setup( kproma, kbdim, itop, klev, klevm1,    &! in
                  & ksfc_type, ktrac, tpfac2, pstep_len,  &! in
                  & pum1, pvm1, pcptgz, pqm1,             &! in
                  & pxlm1, pxim1, pxvar, pxtm1, pxt_emis, &! in
                  & zrdpm, pztkevn, pzthvvar, aa,         &! in
                  & bb, bb_btm                            )! out
  
    CALL rhs_elim ( kproma, kbdim, itop, klev, klevm1, &! in
                  & aa, bb                             )! in, inout

  END SUBROUTINE vdiff_down
  !-------------

END MODULE mo_vdiff_downward_sweep

