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
  USE mo_echam_vdiff_params, ONLY: tpfac1, itop
  USE mo_turbulence_diag,    ONLY: atm_exchange_coeff, sfc_exchange_coeff
  USE mo_vdiff_solver,       ONLY: nvar_vdiff, nmatrix, ih, imh, imqv, &
                                 & matrix_setup_elim, rhs_setup, rhs_elim
  USE mo_exception,          ONLY: message
  USE mo_sma_turbulence_diag,ONLY: atm_exchange_coeff3d, diffuse_hori_velocity, &
                                 & diffuse_vert_velocity,diffuse_scalar
  USE mo_model_domain,       ONLY: t_patch
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_impl_constants,     ONLY: min_rlcell_int, min_rlcell
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_run_config,         ONLY: iqv, iqc, iqi
  USE mo_nh_testcases_nml,   ONLY: isrfc_type

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: vdiff_down

CONTAINS
  !>
  !!
  !!
  SUBROUTINE vdiff_down( jg, kbdim, nblks_c, nblks_v, nblks_e,          &! in
                       & klev, klevm1, klevp1, ktrac,                   &! in
                       & ksfc_type, idx_wtr, idx_ice, idx_lnd,          &! in
                       & pdtime,  pcoriol,                              &! in
                       & turb,                                          &! in
                       & patch,                                         &! in
                       & pzf, pzh,                                      &! in
                       & pfrc,                                          &! in
                       & ptsfc_tile, pocu,      pocv,       ppsfc,      &! in
                       & pum1,       pvm1,      pwp1,                   &! in
                       & ptm1,       pqm1,                              &! in
                       & pxlm1,      pxim1,     pxm1,       pxtm1,      &! in
                       & pmair,      pmref,     rho,                    &! in
                       & paphm1,     papm1,                             &! in
                       & ptvm1,      paclc,     pxt_emis,   pthvvar,    &! in
                       & pxvar,      pz0m_tile,                         &! in
                       & ptottem1,                                      &! in
                       & pustar,     pwstar,    pwstar_tile,            &! inout, out, inout
                       & pqsat_tile, phdtcbl,                           &! out
                       & pri,        pri_tile,  pmixlen,                &! out
                       & pcfm,       pcfm_tile, pcfh,       pcfh_tile,  &! out
                       & pcfv,       pcftotte,  pcfthv,                 &! out
                       & aa,         aa_btm,    bb,         bb_btm,     &! out
                       & ddt_u, ddt_v,                                  &! out
                       & ta_hori_tend,                                  &! out
                       & qv_hori_tend,                                  &! out
                       & ql_hori_tend,                                  &! out
                       & qi_hori_tend,                                  &! out
                       & pfactor_sfc, pcpt_tile,                        &! out
                       & pcptgz,                                        &! out
                       & pzthvvar,   pthvsig,   pztottevn,              &! out
                       & pch_tile,                                      &! out, for "nsurf_diag"
                       & pbn_tile,   pbhn_tile,                         &! out, for "nsurf_diag"
                       & pbm_tile,   pbh_tile,                          &! out, for "nsurf_diag"
                       & pcsat,                                         &! in
                       & pcair,                                         &! in
                       & paz0lh)                                         ! in


    INTEGER, INTENT(IN) :: jg
    INTEGER, INTENT(IN) :: kbdim, klev, klevm1, klevp1, ktrac, nblks_c, nblks_v, nblks_e
    INTEGER, INTENT(IN) :: turb      !< 1: TTE scheme, 2: 3D Smagorisnky
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    REAL(wp),INTENT(IN) :: pdtime

    TYPE(t_patch)   ,TARGET ,INTENT(inout) :: patch

    REAL(wp),INTENT(IN) ::          &
      & pcoriol   (:,:)   ,&!< (kbdim) Coriolis parameter: 2*omega*sin(lat)
      & pzf       (:,:,:) ,&!< (kbdim,klev) geopotential height above sea level, full level
      & pzh       (:,:,:) ,&!< (kbdim,klevp1) geopotential height above sea level, half level
      & pfrc      (:,:,:) ,&!< (kbdim,ksfc_type) area fraction of each surface type
      & ptsfc_tile(:,:,:) ,&!< (kbdim,ksfc_type) surface temperature
      & pocu      (:,:)   ,&!< (kbdim) eastward  velocity of ocean sfc current
      & pocv      (:,:)   ,&!< (kbdim) northward velocity of ocean sfc current
      & ppsfc     (:,:)     !< (kbdim) surface pressure

    REAL(wp),INTENT(IN) ::        &
      & ptm1    (:,:,:)   ,&!< (kbdim,klev) temperature at step t-dt
      & pqm1    (:,:,:)   ,&!< (kbdim,klev) specific humidity at step t-dt
      & pxlm1   (:,:,:)   ,&!< (kbdim,klev) cloud water concentration at step t-dt
      & pxim1   (:,:,:)   ,&!< (kbdim,klev) cloud ice   concentration at step t-dt
      & pxm1    (:,:,:)   ,&!< (kbdim,klev) cloud water + cloud ice at step t-dt
      & pxtm1   (:,:,:,:)   !< (kbdim,klev,ktrac) specific density of other tracers at t-dt

    REAL(wp),INTENT(IN) ::        &
      & pmair   (:,:,:)   ,&!< (kbdim,klev)     air mass [kg/m2]
      & pmref   (:,:,:)   ,&!< (kbdim,klev) dra air mass [kg/m2]
      & paphm1  (:,:,:)   ,&!< (kbdim,klevp1) half level pressure [Pa]
      & papm1   (:,:,:)   ,&!< (kbdim,klev) full level pressure [Pa]
      & ptvm1   (:,:,:)   ,&!< (kbdim,klev) virtual temperature
      & paclc   (:,:,:)   ,&!< (kbdim,klev) cloud fraction
      & pxt_emis(:,:,:)     !< (kbdim,ktrac) tracer tendency due to surface emission
                          !< and dry deposition

    REAL(wp),INTENT(IN) ::         &
      & pthvvar  (:,:,:)  ,&!< (kbdim,klev) variance of virtual pot. temp. at step t-dt
      & pxvar    (:,:,:)  ,&!< (kbdim,klev) step t-dt
      & pz0m_tile(:,:,:)    !< (kbdim,ksfc_type) roughness length at step t-dt

    REAL(wp),INTENT(INOUT) ::        &
      & rho     (:,:,:)   ,&!< (kbdim,klev) air density [kg/m3]
      & pum1    (:,:,:)   ,&!< (kbdim,klev) u-wind at step t-dt
      & pvm1    (:,:,:)     !< (kbdim,klev) q-wind at step t-dt

    REAL(wp),INTENT(IN)  :: ptottem1(:,:,:)    !< (kbdim,klev) TTE at step t-dt

    ! Grid-box mean friction velocity.
    ! In: value at step t-2dt computed in the previous time step,
    ! used in the computation of PBL height (then mixing length);
    ! Out: computed in sfc_exchange_coeff at step t-dt.

    REAL(wp),INTENT(INOUT) :: pustar (:,:)         !< (kbdim)
    REAL(wp),INTENT(OUT)   :: pwstar (:,:)         !< (kbdim)
    REAL(wp),INTENT(INOUT) :: pwstar_tile(:,:,:)   !< (kbdim,ksfc_type)
    REAL(wp),INTENT(INOUT) :: pwp1    (:,:,:)      !< (kbdim,klevp1) vertical wind in m/s
    REAL(wp),INTENT(OUT)   :: ddt_u (:,:,:),      &
                            & ddt_v (:,:,:)
    REAL(wp),INTENT(OUT)   :: ta_hori_tend (:,:,:)
    REAL(wp),INTENT(OUT)   :: qv_hori_tend (:,:,:)
    REAL(wp),INTENT(OUT)   :: ql_hori_tend (:,:,:)
    REAL(wp),INTENT(OUT)   :: qi_hori_tend (:,:,:)

    ! Variables with intent(out)

    REAL(wp),INTENT(OUT) :: pqsat_tile(:,:,:)   !< (kbdim,ksfc_type) saturation specific
                                              !< humidity at sfc.
                                              !< (step t-dt)

    REAL(wp),INTENT(OUT) :: phdtcbl(:,:)    !< (kbdim) height of the top of the atmospheric dry
                                          !< convective boundary layer

    REAL(wp),INTENT(OUT) ::      &   ! out
      & pri      (:,:,:)  ,&!< (kbdim,klev) Richardson number
      & pri_tile (:,:,:)  ,&!< (kbdim,ksfc_type) Richardson number
      & pmixlen  (:,:,:)  ,&!< (kbdim,klev) mixing length
      & pcfm     (:,:,:)  ,&!< (kbdim,klev) exchange coeff. for u, v
      & pcfm_tile(:,:,:)  ,&!< (kbdim,ksfc_type) exchange coeff. for u, v
      & pcfh     (:,:,:)  ,&!< (kbdim,klev) exchange coeff. for heat and tracers
      & pcfh_tile(:,:,:)  ,&!< (kbdim,ksfc_type) exchange coeff. for heat and tracers
      & pcfv     (:,:,:)  ,&!< (kbdim,klev) exchange coeff. for variance of qx
      & pcftotte (:,:,:)  ,&!< (kbdim,klev) exchange coeff. for TTE
      & pcfthv   (:,:,:)    !< (kbdim,klev) exchange coeff. for variance of theta_v

    ! Coefficient matrices and right-hand-side vectors.
    ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)

    REAL(wp),INTENT(OUT) ::           &  ! out
      & aa     (:,:,:,:,:)    ,&!< (kbdim,klev,3,nmatrix) coeff. matrices, all variables
      & aa_btm (:,:,:,imh:,:) ,&!< (kbdim,3,ksfc_type,imh:imqv) last row of coeff. matrix of heat and moisture
      & bb     (:,:,:,:)      ,&!< (kbdim,klev,nvar_vdiff) r.h.s., all variables
      & bb_btm (:,:,ih:,:)      !< (kbdim,ksfc_type,ih:iqv) last row of r.h.s. of heat and moisture

    ! Other variables to be passed on to the second part of turbulence solver

    REAL(wp),INTENT(OUT) ::       &  ! out
      & pfactor_sfc(:,:)  ,&!< (kbdim) prefactor for the exchange coeff.
      & pcpt_tile (:,:,:) ,&!< (kbdim,ksfc_type) dry static energy at surface
      & pcptgz    (:,:,:) ,&!< (kbdim,klev) dry static energy
      & pzthvvar  (:,:,:) ,&!< (kbdim,klev)
      & pthvsig   (:,:)   ,&!< (kbdim)
      & pztottevn (:,:,:)   !< (kbdim,klev) intermediate value of TTE
    REAL(wp) :: jztottevn(kbdim,nblks_c) 

    REAL(wp), INTENT(OUT) :: pch_tile(:,:,:)    !< (kbdim,ksfc_type) out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbn_tile(:,:,:)    !< (kbdim,ksfc_type) out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbhn_tile(:,:,:)   !< (kbdim,ksfc_type) out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbm_tile(:,:,:)    !< (kbdim,ksfc_type) out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbh_tile(:,:,:)    !< (kbdim,ksfc_type) out, for "nsurf_diag"

    REAL(wp), OPTIONAL, INTENT(IN) ::          &
      & pcsat     (:,:)          ,&!< (kbdim) area fraction with wet land surface
      & pcair     (:,:)          ,&!< (kbdim) area fraction with wet land surface
      & paz0lh    (:,:)            !< (kbdim) surface roughness length over land for heat

    ! Local variables

    REAL(wp) :: zghf   (kbdim,klev,nblks_c)   !< geopotential height above ground, full level
    REAL(wp) :: zghh   (kbdim,klevp1,nblks_c) !< geopotential height above ground, full level

    REAL(wp) :: zfactor(kbdim,klev,nblks_c)   !< prefactor for the exchange coefficients
    REAL(wp) :: zrmairm(kbdim,klev,nblks_c)
    REAL(wp) :: zrmairh(kbdim,klevm1,nblks_c)
    REAL(wp) :: zrmrefm(kbdim,klev,nblks_c)

    REAL(wp), DIMENSION(kbdim,klev,nblks_c)   :: km_c
    REAL(wp), DIMENSION(kbdim,klevp1,nblks_v) :: km_iv
    REAL(wp), DIMENSION(kbdim,klevp1,nblks_e) :: km_ie, kh_ie
    REAL(wp), DIMENSION(kbdim,klevp1,nblks_c) :: kh_ic, km_ic
    REAL(wp), DIMENSION(kbdim,klev,nblks_e)   :: vn

    REAL(wp) :: u_vert(kbdim,klev,nblks_v), v_vert(kbdim,klev,nblks_v)
    REAL(wp) :: div_c(kbdim,klev,nblks_c)
    REAL(wp) :: rho_ic(kbdim,klevp1,nblks_c)
    REAL(wp) :: w_vert(kbdim,klevp1,nblks_v)
    REAL(wp) :: w_ie(kbdim,klevp1,nblks_e)

    ! _b denotes value at the bottom level (the klev-th full level)

    REAL(wp) :: ztheta_b (kbdim,nblks_c)  !< potential temperature
    REAL(wp) :: zthetav_b(kbdim,nblks_c)  !< virtual potential temperature
    REAL(wp) :: zthetal_b(kbdim,nblks_c)  !< liquid (and ice?) pot. temp.
    REAL(wp) :: zqsat_b  (kbdim,nblks_c)  !< specific humidity at saturation
    REAL(wp) :: zlh_b    (kbdim,nblks_c)  !< latent heat

    REAL(wp) :: zconst

    INTEGER  :: jl, jk
    INTEGER,parameter :: itte=1, isma=2
    INTEGER, PARAMETER :: tracer_dry_static = 1
    INTEGER, PARAMETER :: tracer_water = 2

    ! OMP variables
    INTEGER                             :: jb,jbs,jbe,jcs,jce,ncd,rls,rle

    !$ACC DATA PRESENT( pxtm1, pxt_emis ) IF( ktrac > 0 )
    !$ACC DATA &
    !$ACC PRESENT(pcoriol,pzf,pzh,pfrc,ptsfc_tile,pocu,pocv,ppsfc) &
    !$ACC PRESENT(pum1,pvm1,pwp1,ptm1,pqm1,pxlm1,pxim1,pxm1) &
    !$ACC PRESENT(pmair,pmref,paphm1,papm1,ptvm1,paclc)   &
    !$ACC PRESENT(pthvvar,pxvar,pz0m_tile,ptottem1)   &
    !---- Argument arrays - intent(inout)
    !$ACC PRESENT(pustar,pwstar,pwstar_tile,pqsat_tile,phdtcbl)    &
    !$ACC PRESENT(pri,pri_tile,pmixlen,pcfm,pcfm_tile,pcfh,pcfh_tile,pcfv,pcftotte)    &
    !$ACC PRESENT(pcfthv,aa,aa_btm,bb,bb_btm,pfactor_sfc,pcpt_tile)    &
    !$ACC PRESENT(pcptgz,pzthvvar,pthvsig,pztottevn)    &
    !---- Argument arrays - intent(out)
    !$ACC PRESENT(pch_tile,pbn_tile,pbhn_tile,pbm_tile,pbh_tile, ddt_u, ddt_v) &
    !---- Optional Argument arrays - intent(in)
    !$ACC PRESENT(pcsat,pcair,paz0lh) &
    !---- Local variables
    !$ACC CREATE(zghf,zghh,zfactor,zrmairm,zrmairh,zrmrefm,jztottevn) &
    !$ACC CREATE(ztheta_b,zthetav_b,zthetal_b,zqsat_b,zlh_b)


    rls = grf_bdywidth_c+1
    rle = min_rlcell_int

    ncd = MAX(1,patch%n_childdom)
    jbs = patch%cells%start_blk(rls,  1)
    jbe = patch%cells%end_blk  (rle,ncd)

!##############################################################################
!## jb loop1
!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = jbs, jbe
      CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
      IF (jcs>jce) CYCLE
!##############################################################################
    !----------------------------------------------------------------------
    ! 0. Compute useful local fields
    !----------------------------------------------------------------------

    ! geopotential height above ground

    !$ACC PARALLEL LOOP DEFAULT(NONE) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jk = 1,klev
         DO jl = jcs,jce
            zghf(jl,jk,jb) = pzf(jl,jk,jb) - pzh(jl,klevp1,jb)
         END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(NONE) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jk = 1,klevp1
         DO jl = jcs,jce
            zghh (jl,jk,jb) = pzh(jl,jk,jb) - pzh(jl,klevp1,jb)
         END DO
    END DO
    !$ACC END PARALLEL LOOP
    
    ! reciprocal layer mass
    !$ACC PARALLEL LOOP DEFAULT(NONE) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jk = 1,klev
         DO jl = jcs,jce
            zrmairm(jl,jk,jb) = 1._wp / pmair(jl,jk,jb)
            zrmrefm(jl,jk,jb) = 1._wp / pmref(jl,jk,jb)
         END DO
    END DO
    !$ACC END PARALLEL LOOP
    
    !$ACC PARALLEL LOOP DEFAULT(NONE) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jk = 1,klevm1
      DO jl = jcs,jce
        zrmairh(jl,jk,jb) = 2._wp / (pmair(jl,jk,jb) + pmair(jl,jk+1,jb))
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(NONE) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jk = 1,klev
       DO jl = jcs,jce
         ddt_u(jl,jk,jb) = 0._wp
         ddt_v(jl,jk,jb) = 0._wp
       END DO
    END DO
    !$ACC END PARALLEL LOOP
    
!##############################################################################
!## jb end loop1
   END DO
!$OMP END PARALLEL DO
!##############################################################################

SELECT CASE ( turb ) ! select turbulent scheme
CASE ( itte ) ! TTE scheme
    !----------------------------------------------------------------------
    ! 1. Compute various thermodynamic variables; Diagnose PBL extension;
    !    Compute exchange coefficients for half levels [1+1/2, klev-1/2];
    !    Get TTE and variance of theta_v at intermediate time step.
    !----------------------------------------------------------------------

    rls = grf_bdywidth_c+1
    rle = min_rlcell_int

    ncd = MAX(1,patch%n_childdom)
    jbs = patch%cells%start_blk(rls,  1)
    jbe = patch%cells%end_blk  (rle,ncd)
!##############################################################################
!## jb loop2
!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = jbs, jbe
      CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
      IF (jcs>jce) CYCLE
!##############################################################################

    ! DA: this routine is async aware, so it's safe not not wait here
    CALL atm_exchange_coeff( jg,                                                         &! in
                           & jb,                                                         &! in, for debugging only
                           & jcs, jce, kbdim, klev, klevm1, klevp1,                      &! in
                           & pdtime, pcoriol(:,jb),                                      &! in
                           & zghf(:,:,jb), zghh(:,:,jb),                                 &! in
                           & pum1(:,:,jb), pvm1(:,:,jb), ptm1(:,:,jb), ptvm1(:,:,jb),    &! in
                           & pqm1(:,:,jb), pxm1(:,:,jb),                                 &! in
                           & papm1(:,:,jb), paphm1(:,:,jb), paclc(:,:,jb), pustar(:,jb), &! in
                           & pthvvar(:,:,jb), ptottem1(:,:,jb),                          &! in
                           & pcptgz(:,:,jb), phdtcbl(:,jb),                              &! out
                           & pzthvvar(:,1:klevm1,jb),                                    &! out
                           & pztottevn(:,1:klevm1,jb),                                   &! out
                           & pcfm    (:,1:klevm1,jb), pcfh  (:,1:klevm1,jb),             &! out
                           & pcfv    (:,1:klevm1,jb),                                    &! out
                           & pcftotte(:,1:klevm1,jb),                                    &! out
                           & pcfthv  (:,1:klevm1,jb),zfactor(:,1:klevm1,jb),             &! out
                           & ztheta_b(:,jb), zthetav_b(:,jb), zthetal_b(:,jb),           &! out, for "sfc_exchange_coeff"
                           & zqsat_b(:,jb),  zlh_b(:,jb),                                &! out, for "sfc_exchange_coeff"
                           & pri(:,1:klevm1,jb), pmixlen(:,1:klevm1,jb)                  )! out, for output

    !$ACC PARALLEL LOOP DEFAULT(NONE) GANG VECTOR ASYNC(1)
    DO jl = jcs,jce
      pmixlen(jl,klev,jb) = -999._wp
    END DO
    !$ACC END PARALLEL LOOP

    !-----------------------------------------------------------------------
    ! 2. Compute exchange coefficients at the air-sea/ice/land interface.
    !    Get boundary condition for TTE and variance of theta_v.
    !-----------------------------------------------------------------------

    ! DA: this routine is async, no need to wait
    CALL sfc_exchange_coeff( jg,                                            &! in
                           & jcs, jce, kbdim, ksfc_type,                    &! in
                           & idx_wtr, idx_ice, idx_lnd,                     &! in
                           & pz0m_tile(:,jb,:),  ptsfc_tile(:,jb,:),        &! in
                           & pfrc(:,jb,:),       phdtcbl(:,jb),             &! in
                           & pocu(:,jb),         pocv(:,jb),   ppsfc(:,jb), &! in
                           & zghf(:,klev,jb),                               &! in
                           & pum1(:,klev,jb),    pvm1  (:,klev,jb),         &! in
                           & ptm1(:,klev,jb),                               &! in
                           & pqm1(:,klev,jb),    pxm1  (:,klev,jb),         &! in
                           & zqsat_b  (:,jb),    zlh_b    (:,jb),           &! in
                           & ztheta_b (:,jb),    zthetav_b(:,jb),           &! in
                           & zthetal_b(:,jb),    paclc (:,klev,jb),         &! in
                           & ptottem1(:,klev,jb),pzthvvar(:,klevm1,jb),     &! in
                           & pthvsig(:,jb),                                 &! out
                           & pwstar(:,jb),       pwstar_tile(:,jb,:),       &! out, inout
                           & pqsat_tile(:,jb,:), pcpt_tile(:,jb,:),         &! out
                           & pri    (:,klev,jb), pri_tile(:,jb,:),          &! out
                           & pcfm   (:,klev,jb), pcfm_tile(:,jb,:),         &! out
                           & pcfh   (:,klev,jb), pcfh_tile(:,jb,:),         &! out
                           & pcfv   (:,klev,jb),                            &! out
                           & pcftotte(:,klev,jb),pcfthv  (:,klev,jb),       &! out
                           & zfactor(:,klev,jb),                            &! out
                           & pztottevn(:,klev,jb),                          &! out
                           & pzthvvar(:,klev,jb),                           &! out
                           & jztottevn(:,jb),                               &! out
                           & pustar(:,jb),                                  &! out, for "atm_exchange_coeff" at next time step
                           & pch_tile(:,jb,:),                              &! out, for "nsurf_diag"
                           & pbn_tile(:,jb,:),   pbhn_tile(:,jb,:),         &! out, for "nsurf_diag"
                           & pbm_tile(:,jb,:),   pbh_tile(:,jb,:),          &! out, for "nsurf_diag"
                           & paz0lh(:,jb),                                  &! in, optional
                           & pcsat(:,jb),                                   &! in, optional
                           & pcair(:,jb))                                    ! in, optional

    IF ( isrfc_type == 1 ) THEN
      !$ACC PARALLEL LOOP DEFAULT(NONE) GANG VECTOR ASYNC(1)
      DO jl = jcs,jce
        pztottevn(jl,klev,jb) = jztottevn(jl,jb)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

!##############################################################################
!## jb end loop2
   END DO
!$OMP END PARALLEL DO
!##############################################################################

CASE ( isma ) ! 3D Smagorinksy scheme

    CALL atm_exchange_coeff3d ( jg, kbdim, nblks_c, nblks_v, nblks_e,                &! in
                           & klev, klevm1, klevp1,                                   &! in
                           & ksfc_type, idx_lnd,                                     &! in
                           & patch,                                                  &! in
                           & pz0m_tile(:,:,:), ptsfc_tile(:,:,:), pfrc(:,:,:),       &! in
                           & ppsfc(:,:),                                             &! in
                           & zghf(:,:,:), zghh(:,:,:),                               &! in
                           & pum1(:,:,:), pvm1(:,:,:), pwp1(:,:,:),                  &! in
                           & ptm1(:,:,:), ptvm1(:,:,:),                              &! in
                           & pqm1(:,:,:), pxm1(:,:,:),                               &! in
                           & rho(:,:,:),                                             &! in
                           & papm1(:,:,:), paphm1(:,:,:),                            &! in
                           & pri_tile(:,:,:),                                        &! out
                           & pthvsig(:,:),                                           &! out
                           & pcfm_tile(:,:,:),                                       &! out
                           & pcfh_tile(:,:,:),                                       &! out
                           & pqsat_tile(:,:,:), pcpt_tile(:,:,:),                    &! out
                           & pcptgz(:,:,:),                                          &! out
                           & pzthvvar(:,:,:),                                        &! out
                           & pztottevn(:,:,:), pmixlen(:,:,:),                       &! out
                           & pcfm    (:,:,:), pcfh  (:,:,:),                         &! out
                           & pcfv    (:,:,:),                                        &! out
                           & pcftotte(:,:,:),                                        &! out
                           & pcfthv  (:,:,:),                                        &! out
                           & km_c(:,:,:), km_iv(:,:,:), km_ie(:,:,:),                &! out
                           & km_ic(:,:,:), kh_ic(:,:,:),                             &! out,
                           & zfactor(:,:,:),                                         &! out
                           & u_vert(:,:,:), v_vert(:,:,:), div_c(:,:,:),             &! out, for "sfc_exchange_coeff"
                           & rho_ic(:,:,:), w_vert(:,:,:), w_ie(:,:,:),              &! out
                           & vn(:,:,:),                                              &! out,
                           & pch_tile(:,:,:),                                        &! out, for "nsurf_diag"
                           & pbn_tile(:,:,:), pbhn_tile(:,:,:),                      &! out
                           & pbm_tile(:,:,:), pbh_tile(:,:,:),                       &! out
                           & paz0lh=paz0lh(:,:),                                     &! in, optional
                           & pcsat=pcsat(:,:), pcair=pcair(:,:)                      &! in, optional
                           & )


    CALL diffuse_hori_velocity( kbdim,                                            &
                              & patch,                                            &
                              & km_c(:,:,:), km_iv(:,:,:), km_ie(:,:,:),          &
                              & u_vert(:,:,:), v_vert(:,:,:), div_c(:,:,:),       &
                              & rho(:,:,:), pum1(:,:,:), pvm1(:,:,:), vn(:,:,:),  &
                              & ddt_u(:,:,:), ddt_v(:,:,:), pdtime)

    CALL diffuse_vert_velocity( kbdim,                                            &
                              & patch,                                            &
                              & rho_ic(:,:,:), w_vert(:,:,:), w_ie(:,:,:),        &
                              & km_c(:,:,:), km_iv(:,:,:), km_ic(:,:,:),          &
                              & u_vert(:,:,:), v_vert(:,:,:), div_c(:,:,:),       &
                              & pum1(:,:,:), pvm1(:,:,:), pwp1(:,:,:), vn(:,:,:), &
                              & pdtime)


! dry static energy
    call diffuse_scalar( kbdim, ptm1(:,:,:),               &
                       & patch,                            &
                       & kh_ic(:,:,:), km_ie(:,:,:),       &
                       & ta_hori_tend(:,:,:),              &
                       & rho,                              &
                       & tracer_dry_static)

    call diffuse_scalar( kbdim, pqm1(:,:,:),               &
                       & patch,                            &
                       & kh_ic(:,:,:), km_ie(:,:,:),       &
                       & qv_hori_tend(:,:,:),              &
                       & rho,                              &
                       & tracer_water)

    call diffuse_scalar( kbdim, pxlm1(:,:,:),              &
                       & patch,                            &
                       & kh_ic(:,:,:), km_ie(:,:,:),       &
                       & ql_hori_tend(:,:,:),              &
                       & rho,                              &
                       & tracer_water)

    call diffuse_scalar( kbdim, pxim1(:,:,:),              &
                       & patch,                            &
                       & kh_ic(:,:,:), km_ie(:,:,:),       &
                       & qi_hori_tend(:,:,:),              &
                       & rho,                              &
                       & tracer_water)

END SELECT    !select turbulent scheme

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

    rls = grf_bdywidth_c+1
    rle = min_rlcell_int

    ncd = MAX(1,patch%n_childdom)
    jbs = patch%cells%start_blk(rls,  1)
    jbe = patch%cells%end_blk  (rle,ncd)
!##############################################################################
!## jb loop3
!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = jbs, jbe
      CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
      IF (jcs>jce) CYCLE
!##############################################################################

    zconst = tpfac1*pdtime

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jl = jcs,jce
        zfactor(jl,jk,jb) = zfactor(jl,jk,jb)*zconst
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs,jce
      zfactor(jl,klev,jb) = zfactor(jl, klev,jb)*zconst
    END DO
    !$ACC END PARALLEL

    CALL matrix_setup_elim( jcs, jce, kbdim, klev, klevm1, ksfc_type, itop, &! in
                          & pcfm     (:,:,jb),   pcfh  (:,1:klevm1,jb),         &! in
                          & pcfh_tile(:,jb,:),   pcfv  (:,:,jb),                &! in
                          & pcftotte (:,:,jb),   pcfthv(:,:,jb),                &! in
                          & zfactor  (:,:,jb),                                  &! in
                          & zrmairm(:,:,jb), zrmairh(:,:,jb), zrmrefm(:,:,jb),  &! in
                          & aa(:,:,:,:,jb), aa_btm(:,:,:,:,jb)                  )! out

    ! Save for output, to be used in "update_surface"
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs,jce
      pfactor_sfc(jl,jb) = zfactor(jl,klev,jb)
    END DO
    !$ACC END PARALLEL

    !-----------------------------------------------------------------------
    ! 4. Set up right-hand side of the tri-diagonal system and perform
    !    Gauss elimination. Factors that determine the r.h.s. include
    !    - time stepping scheme used for vertical diffusion
    !    - whether there is any other process (e.g., tracer emission)
    !      solved together with vertical diffusion.
    !-----------------------------------------------------------------------

    CALL rhs_setup( jcs, jce, kbdim, itop, klev, klevm1,                                            &! in
                  & ksfc_type, ktrac, pdtime,                                                       &! in
                  & pum1(:,:,jb), pvm1(:,:,jb), pcptgz(:,:,jb), pqm1(:,:,jb),                       &! in
                  & pxlm1(:,:,jb), pxim1(:,:,jb), pxvar(:,:,jb), pxtm1(:,:,jb,:), pxt_emis(:,:,jb), &! in
                  & zrmrefm(:,:,jb), pztottevn(:,:,jb), pzthvvar(:,:,jb), aa(:,:,:,:,jb),           &! in
                  & bb(:,:,:,jb), bb_btm(:,:,:,jb)                                                  )! out

    CALL rhs_elim ( jcs, jce, kbdim, itop, klev, klevm1, &! in
                  & aa(:,:,:,:,jb), bb(:,:,:,jb)         )! in, inout

!##############################################################################
!## jb end loop3
   END DO
!$OMP END PARALLEL DO
!##############################################################################

    !$ACC WAIT
    !$ACC END DATA
    !$ACC END DATA

  END SUBROUTINE vdiff_down
  !-------------

END MODULE mo_vdiff_downward_sweep
