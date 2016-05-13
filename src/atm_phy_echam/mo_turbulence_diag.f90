#include "fsel.inc"
!>
!! @brief Subroutines for computing turbulent exchange coefficients.
!!
!! @author Authors of ECHAM code
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! Code originates from ECHAM5/6
!! Re-organized by Hui Wan (2010-09).
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_turbulence_diag

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
  USE mo_convect_tables,    ONLY: prepare_ua_index_spline, lookup_ua_spline, &
    &                             compute_qsat
  USE mo_echam_vdiff_params,ONLY: clam, cgam, ckap, cb,cc, chneu, shn, smn, &
    &                             da1, custf, cwstf, cfreec, tpfac1,        &
    &                             eps_shear, eps_corio, tke_min,            &
    &                             cons2, cons25, cons5
  USE mo_physical_constants,ONLY: grav, rd, cpd, cpv, rd_o_cpd, rv,         &
    &                             vtmpc1, tmelt, alv, als, p0ref
#ifdef __ICON__
  USE mo_echam_phy_config,  ONLY: phy_config => echam_phy_config
#else
  USE mo_time_control,      ONLY: lstart
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: atm_exchange_coeff, sfc_exchange_coeff

CONTAINS
  !>
  !! Compute various thermodynamic variables for all (full) vertical levels;
  !! Diagnose PBL extension;
  !! Diagnose wind shear, buoyancy, Ri-number, mixing length, then compute
  !! the turbulent exchange coefficients of momentum, dry static energy,
  !! tracers, TKE, variance of virtual optential temperation at half levels
  !! [1+1/2, klev-1/2].
  !! Note that
  !! - for all coeffcient arrays, vertical index k in this subroutine
  !!   correspond to interface (half level) k+1/2;
  !! - the exchange coefficient at model top (level 1/2) is zero, thus does
  !!   not need computing;
  !! - the exchange coefficients at the Earth's surface are computed in
  !!   subroutine "sfc_exchanged_coeff".
  !!
  !! @par Revision History
  !! Separated from vdiff.f90 of ECHAM6 and re-organized by Hui Wan (2010-09).
  !!  updated to echam-6.3.01 by Monika Esch (2014-11)
  !!
  SUBROUTINE atm_exchange_coeff( kproma, kbdim, klev, klevm1, klevp1,     &! in
#ifndef __ICON__
                               & lstart,                                  &! in
#endif
                               & pstep_len, pcoriol,                      &! in
                               & pum1, pvm1, ptm1, ptvm1, pgeom1,         &! in
                               & pqm1, pxm1,                              &! in
                               & papm1, paphm1, paclc,                    &! in
                               & pustarm, pthvvar,                        &! in
#ifdef __ICON__
                               & ptkem1,                                  &! in
#else
                               & ptkem1, ptkem0,                          &! inout
#endif
                               & pcptgz, ihpbl, pghabl,                   &! out
                               & pqshear, pzthvvar, ptkevn,               &! out
                               & pcfm, pcfh, pcfv, pcftke, pcfthv, pprfac,&! out
                               & prhoh, ptheta_b, pthetav_b, pthetal_b,   &! out
                               & pqsat_b, plh_b,                          &! out
                               & pri, pmixlen                             )
    ! Arguments

    INTEGER, INTENT(IN) :: kproma, kbdim
    INTEGER, INTENT(IN) :: klev, klevm1, klevp1
#ifndef __ICON__
    LOGICAL, INTENT(IN) :: lstart
#endif
    REAL(wp),INTENT(IN) :: pstep_len
    REAL(wp),INTENT(IN) :: pcoriol(kbdim)   !< Coriolis parameter: 2*omega*sin(lat)
    REAL(wp),INTENT(IN) :: pxm1(kbdim,klev)
    REAL(wp),INTENT(IN) :: ptvm1(kbdim,klev)
    REAL(wp),INTENT(IN) :: pqm1(kbdim,klev),   pum1(kbdim,klev),  pvm1(kbdim,klev)
    REAL(wp),INTENT(IN) :: ptm1(kbdim,klev)
    REAL(wp),INTENT(IN) :: paclc(kbdim,klev)
    REAL(wp),INTENT(IN) :: papm1(kbdim,klev),  paphm1(kbdim,klevp1)
    REAL(wp),INTENT(IN) :: pthvvar(kbdim,klev)
    REAL(wp),INTENT(IN) :: pgeom1(kbdim,klev)
    REAL(wp),INTENT(IN) :: pustarm(kbdim)

#ifdef __ICON__
    REAL(wp),INTENT(IN) :: ptkem1 (kbdim,klev)
#else
    REAL(wp),INTENT(INOUT) :: ptkem1 (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: ptkem0 (kbdim,klev)
#endif

    INTEGER, INTENT(OUT) :: ihpbl   (kbdim)        !< grid level index of PBL top
    REAL(wp),INTENT(OUT) :: pghabl  (kbdim)        !< geopotential of PBL top
    REAL(wp),INTENT(OUT) :: ptkevn  (kbdim,klevm1) !< TKE at intermediate time step
    REAL(wp),INTENT(OUT) :: pcftke  (kbdim,klevm1) !< exchange coeff. for TKE
    REAL(wp),INTENT(OUT) :: pcfthv  (kbdim,klevm1) !< exchange coeff. for variance of theta_v
    REAL(wp),INTENT(OUT) :: pqshear (kbdim,klevm1) !< vertical gradient of qv
    REAL(wp),INTENT(OUT) :: pcfm    (kbdim,klevm1) !< exchange coeff. for u, v
    REAL(wp),INTENT(OUT) :: pcfh    (kbdim,klevm1) !< exchange coeff. for cptgz and tracers
    REAL(wp),INTENT(OUT) :: pcfv    (kbdim,klevm1) !< exchange coeff. for variance of qx
    REAL(wp),INTENT(OUT) :: pzthvvar(kbdim,klevm1) !< variance of theta_v at intermediate step
    REAL(wp),INTENT(OUT) :: pcptgz  (kbdim,klev)   !< dry static energy
    REAL(wp),INTENT(OUT) :: pprfac  (kbdim,klevm1) !< prefactor for the exchange coefficients
    REAL(wp),INTENT(OUT) :: prhoh   (kbdim,klevm1) !< air density at half levels

    ! _b denotes the values at the bottom level (the klev-th full level)
    REAL(wp),INTENT(OUT) :: ptheta_b (kbdim)  !< potential temperature
    REAL(wp),INTENT(OUT) :: pthetav_b(kbdim)  !< virtual potential temperature
    REAL(wp),INTENT(OUT) :: pthetal_b(kbdim)  !< liquid (and ice?) potential temperature
    REAL(wp),INTENT(OUT) :: pqsat_b  (kbdim)  !< specific humidity at saturation
    REAL(wp),INTENT(OUT) :: plh_b    (kbdim)  !< latent heat

    ! Just for output
    REAL(wp),INTENT(OUT) :: pri     (kbdim,klevm1) !< moist Richardson number at half levels
    REAL(wp),INTENT(OUT) :: pmixlen (kbdim,klevm1) !< mixing length

    ! Local variables
    ! - Variables defined at full levels

    REAL(wp) :: zlh    (kbdim,klev)  !< latent heat at full levels
    REAL(wp) :: ztheta (kbdim,klev)  !< potential temperature
    REAL(wp) :: zthetav(kbdim,klev)  !< virtual potential temperature
    REAL(wp) :: zthetal(kbdim,klev)  !< liquid (and ice?) potential temperature
    REAL(wp) :: zqsat  (kbdim,klev)  !< specific humidity at saturation

    ! - Variables defined at half levels

    REAL(wp) :: zlhh  (kbdim,klevm1) !< latent heat at half levels
    REAL(wp) :: zdgh  (kbdim,klevm1) !< geopotential difference between two full levels

    REAL(wp) :: zccover (kbdim,klevm1), zqxmit  (kbdim,klevm1)       &
               ,zqmit   (kbdim,klevm1), zqsatm  (kbdim,klevm1)       &
               ,zthetah (kbdim,klevm1), zthetavh(kbdim,klevm1)       &
               ,ztmitte (kbdim,klevm1)

    ! - 1D variables and scalars

    INTEGER  :: ihpblc(kbdim),          ihpbld(kbdim),          idx(kbdim)
    INTEGER  :: jk, jl, it
    REAL(wp) :: za(kbdim),              zhdyn(kbdim)                          &
               ,zpapm1i(kbdim),         zua(kbdim)
    REAL(wp) :: z2geomf, zalf, zalh2, zbet, zbuoy
    REAL(wp) :: zcons18, zcons23
    REAL(wp) :: zcor, zdisl, zdusq, zdvsq
    REAL(wp) :: zdqtot, zds, zdus1, zdus2, zdz
    REAL(wp) :: zes, zfox, zfux, zhexp
    REAL(wp) :: zktest, zlam
    REAL(wp) :: zmix, zmult1, zmult2, zmult3, zmult4
    REAL(wp) :: zmult5, zqddif, zqtmit, zrdp
    REAL(wp) :: zrdrv, zri, zrvrd, zsdep1
    REAL(wp) :: zsdep2, zsh, zshear, zsm
    REAL(wp) :: zteldif, ztkesq
    REAL(wp) :: zthvprod, zthvdiss, zthvirdif
    REAL(wp) :: zucf, zusus1
    REAL(wp) :: zz2geo, zzb, ztvm

    REAL(wp) :: zonethird

    !-------------------------------------
    ! 1. Some constants
    !-------------------------------------
    zrvrd     = rv/rd
    zrdrv     = rd/rv
    zonethird = 1._wp/3._wp
    zcons18   = tpfac1*pstep_len*grav**2

    !-------------------------------------
    ! 2. NEW THERMODYNAMIC VARIABLES
    !-------------------------------------

    DO 212 jk=1,klev
      CALL prepare_ua_index_spline('vdiff (1)',kproma,ptm1(1,jk),idx(1),za(1))
      CALL lookup_ua_spline(kproma,idx(1),za(1),zua(1))

      zpapm1i(1:kproma) = 1._wp/papm1(1:kproma,jk)
      ztheta(1:kproma,jk) = (p0ref*zpapm1i(1:kproma))**rd_o_cpd

      DO 211 jl=1,kproma

        ! Virtual dry static energy, potential temperature, virtual potential
        ! temperature

        pcptgz (jl,jk) = pgeom1(jl,jk)+ptm1(jl,jk)*(cpd+(cpv-cpd)*pqm1(jl,jk))
        ztheta (jl,jk) = ptm1(jl,jk)*ztheta(jl,jk)
        zthetav(jl,jk) = ztheta(jl,jk)*(1._wp+vtmpc1*pqm1(jl,jk)-pxm1(jl,jk))

        ! Latent heat, liquid (and ice) potential temperature

        zlh(jl,jk)     = FSEL(ptm1(jl,jk)-tmelt,alv,als) ! latent heat
        zusus1         = zlh(jl,jk)/cpd*ztheta(jl,jk)/ptm1(jl,jk)*pxm1(jl,jk)
        zthetal(jl,jk) = ztheta(jl,jk)-zusus1

        zes=zua(jl)*zpapm1i(jl)              ! (sat. vapour pressure)*Rd/Rv/ps
        zes=MIN(zes,0.5_wp)
        zqsat(jl,jk)=zes/(1._wp-vtmpc1*zes)  ! specific humidity at saturation
211   END DO
212 END DO

    ! Copy bottom-level values to dummy arguments

    ptheta_b (1:kproma) = ztheta (1:kproma,klev)
    pthetav_b(1:kproma) = zthetav(1:kproma,klev)
    pthetal_b(1:kproma) = zthetal(1:kproma,klev)
    pqsat_b  (1:kproma) = zqsat  (1:kproma,klev)
    plh_b    (1:kproma) = zlh    (1:kproma,klev)

    !---------------------------------------------------------
    ! 3. Preparation for computation of exchange coefficients
    !---------------------------------------------------------
    ! Vertical interpolation from full levels to half levels
    ! using linear interpolation in pressure coordinate

    DO 214 jk=1,klevm1
      DO 213 jl=1,kproma
        zdgh(jl,jk) = pgeom1(jl,jk) - pgeom1(jl,jk+1)

        ! interpolation coefficients
        zrdp   = 1._wp/(paphm1(jl,jk) - paphm1(jl,jk+2))
        zsdep1 = (paphm1(jl,jk)  - paphm1(jl,jk+1))*zrdp
        zsdep2 = (paphm1(jl,jk+1)- paphm1(jl,jk+2))*zrdp

        ! vertical interpolation for various variables
        zqsatm(jl,jk)  =zsdep1*zqsat(jl,jk)+zsdep2*zqsat(jl,jk+1)
        ztmitte(jl,jk)=zsdep1*ptm1(jl,jk)+zsdep2*ptm1(jl,jk+1)
        zthetavh(jl,jk)=zsdep1*zthetav(jl,jk)+zsdep2*zthetav(jl,jk+1)
        zlhh(jl,jk)=zsdep1*zlh(jl,jk)+zsdep2*zlh(jl,jk+1)
        zqxmit(jl,jk)=zsdep1*pxm1(jl,jk)+zsdep2*pxm1(jl,jk+1)
        zqmit(jl,jk)=zsdep1*pqm1(jl,jk)+zsdep2*pqm1(jl,jk+1)
        zthetah(jl,jk)=zsdep1*ztheta(jl,jk)+zsdep2*ztheta(jl,jk+1)
        zccover(jl,jk)=paclc(jl,jk)*zsdep1+paclc(jl,jk+1)*zsdep2
213   END DO
214 END DO

    !-----------------------------------------------
    ! 4. Compute planetary boundary layer extension
    !-----------------------------------------------
    ! (ECHAM6) JSBACH note: in standard ECHAM5, ustarm is computed before zhdyn for
    ! the current time step. Here, ustarm comes from the call to
    ! mo_surface at the previous timestep.
    ! But this should have only a minor effect on ihpbl and ghabl.
    ! ICON: Initial value of ustar has been set in subroutine "init_phy_memory".
    ! ECHAM: Initial value of ustar has been set in subroutine "physc".

    DO jl = 1,kproma
      zcor=MAX(ABS(pcoriol(jl)),eps_corio)
      zhdyn(jl)=MIN(pgeom1(jl,1)/grav,chneu*pustarm(jl)/zcor)
      ihpblc(jl)=klev
      ihpbld(jl)=klev
    END DO

    DO jk=klevm1,1,-1
      DO jl=1,kproma
        zds=pcptgz(jl,jk)-pcptgz(jl,klev)
        zdz=pgeom1(jl,jk)/grav-zhdyn(jl)
        ihpblc(jl)=MERGE(jk,ihpblc(jl),ihpblc(jl).EQ.klev.AND.zds.GT.0._wp)
        ihpbld(jl)=MERGE(jk,ihpbld(jl),ihpbld(jl).EQ.klev.AND.zdz.GE.0._wp)
      END DO
    END DO

    DO jl=1,kproma
      ihpbl (jl) = MIN(ihpblc(jl),ihpbld(jl))
      pghabl(jl) = MIN(50000._wp,pgeom1(jl,ihpbl(jl)))
    END DO

    !--------------------------------------------
    ! 5. Compute exchange coefficients
    !--------------------------------------------
    DO 372 jk=1,klevm1
      DO 361 jl=1,kproma

        ! gradient of specific humidity, wind shear, buoyancy, Ri-number
        ! according to Brinkop and Roeckner (1995, Tellus A)

        zqtmit=zqxmit(jl,jk)+zqmit(jl,jk)                               ! qt
        zfux=zlhh(jl,jk)/(cpd*ztmitte(jl,jk))                           ! L/(cpd*T)
        zfox=zlhh(jl,jk)/(rd*ztmitte(jl,jk))                            ! L/(Rd*T)
        zmult1=1._wp+vtmpc1*zqtmit                                      ! (1+0.61*qt) = A in clear sky
        zmult2=zfux*zmult1-zrvrd
        zmult3=zrdrv*zfox*zqsatm(jl,jk)/(1._wp+zrdrv*zfux*zfox*zqsatm(jl,jk))
        zmult5=zmult1-zmult2*zmult3                                     ! A in cloud
        zmult4=zfux*zmult5-1._wp                                        ! D in cloud
        zdus1=zccover(jl,jk)*zmult5+(1._wp-zccover(jl,jk))*zmult1       ! A avg
        zdus2=zccover(jl,jk)*zmult4+(1._wp-zccover(jl,jk))*vtmpc1       ! D avg
        zteldif  =(zthetal(jl,jk)-zthetal(jl,jk+1))/zdgh(jl,jk)*grav    ! d theta_l
        zthvirdif=(zthetav(jl,jk)-zthetav(jl,jk+1))/zdgh(jl,jk)*grav    ! d theta_v
        zdqtot=(pqm1(jl,jk)+pxm1(jl,jk))-(pqm1(jl,jk+1)+pxm1(jl,jk+1))  ! d qt
        zqddif=zdqtot/zdgh(jl,jk)*grav                                  ! (d qt)/(d z)
        zbuoy  = (zteldif*zdus1+zthetah(jl,jk)*zdus2*zqddif)            &
                *grav/zthetavh(jl,jk)
        zdusq  = (pum1(jl,jk)-pum1(jl,jk+1))**2
        zdvsq  = (pvm1(jl,jk)-pvm1(jl,jk+1))**2
        zshear = (zdusq+zdvsq)*(grav/zdgh(jl,jk))**2
        zri    = zbuoy/MAX(zshear,eps_shear)

        pqshear(jl,jk) = zqddif   ! store for variance production
        pri(jl,jk) = zri          ! save for output

       ! Asymptotic mixing length for momentum and
       ! heat (zlam) above the PBL as a function of height
       ! according to Holtslag and Boville (1992), J. CLIMATE.

        zhexp=EXP(1._wp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))
        zlam=1._wp+(clam-1._wp)*zhexp
        IF(jk.GE.ihpbl(jl)) THEN
           zcons23=cons25        ! (1/lambda)*(0.5/grav)
        ELSE
           zcons23=cons2/zlam    ! (1/lambda)*(0.5/grav)
        END IF

        ! Mixing length (Blackadar)

        z2geomf = pgeom1(jl,jk)+pgeom1(jl,jk+1)  ! half-level value of gz
        zz2geo  = cons2*z2geomf                  ! z*(Karman constant)
        zmix    = zz2geo/(1._wp+zcons23*z2geomf) ! mixing length

        pmixlen(jl,jk) = zmix   ! save for output

        ! Stability functions (Louis, 1979)

        IF(zri.LT.0._wp) THEN  ! unstable condition
           zalh2=zmix*zmix
           zucf=1._wp/                                              &
                (1._wp+cons5*zalh2*SQRT(ABS(zri)*(((pgeom1(jl,jk)   &
                    /pgeom1(jl,jk+1))**zonethird-1._wp)/(pgeom1(jl,jk) &
                    -pgeom1(jl,jk+1)))**3/(pgeom1(jl,jk+1))))
           zsh=shn*(1._wp-3._wp*cb*zri*zucf)*zmix
           zsm=smn*(1._wp-2._wp*cb*zri*zucf)*zmix

        ELSE  ! stable condition
           zsh=shn/(1._wp+2._wp*cb*zri*SQRT(1._wp+zri))*zmix
           zsm=smn/(1._wp+2._wp*cb*zri/SQRT(1._wp+zri))*zmix
        END IF

        ! TKE at intermediate time step, obtained by solving a prognostic equation
        ! of TKE, the right-hand side of which consists of shear production,
        ! buoyancy production, and dissipation of TKE.
        ! See Appendix A in Brinkop and Roeckner (1995, Tellus) for the numerics.

        zzb=zshear*zsm-zbuoy*zsh
        zdisl=da1*zmix/pstep_len
        zktest=1._wp+(zzb*pstep_len+SQRT(ptkem1(jl,jk))*2._wp)/zdisl
        IF (zktest.LE.1._wp) THEN
           ptkevn(jl,jk)=tke_min
        ELSE
           ptkevn(jl,jk)=MAX(tke_min,(zdisl*(SQRT(zktest)-1._wp))**2)
        END IF

#ifdef __ICON__
        ! Should we do the same as in echam, or should we do the initialization
        ! in, e.g., init_vdiff_xxx?
#else
        ! For ECHAM: initialize tkem1 and tkem0 at the first time step

        IF (lstart) THEN
           ptkem1(jl,jk) = ptkevn(jl,jk)
           ptkem0(jl,jk) = ptkevn(jl,jk)
        END IF
#endif

        ! Square root of TKE at the old time step

        ztkesq = SQRT(MAX(tke_min,ptkem1(jl,jk)))

        ! Virtual potential temperautre variance at intermediate step,
        ! obtained by solving a prognostic equation the variance, the right-hand side
        ! of which consists of a production term and dissipation.
        ! An explicit time stepping method is used here.

        zthvprod = 2._wp*zsh*ztkesq*zthvirdif**2       ! production rate
        zthvdiss = pthvvar(jl,jk)*ztkesq/(da1*zmix)    ! dissipation rate
        pzthvvar(jl,jk) = pthvvar(jl,jk)+(zthvprod-zthvdiss)*pstep_len
        pzthvvar(jl,jk) = MAX(tke_min,pzthvvar(jl,jk))

        ! Exchange coefficients for
        ! - momentum  (variable pcfm),
        ! - heat and tracers (variable pcfh),
        ! - variance of hydrometeors (variable pcfv).
        ! These are proportional to the square root of TKE at the old step.

        pcfm(jl,jk) = zsm*ztkesq
        pcfh(jl,jk) = zsh*ztkesq
        pcfv(jl,jk) = zsh*ztkesq*0.5_wp

        ! Exchange coefficients for
        ! - TKE (variable pfctke),
        ! - and variance of virtual potential temperature (variable pcfthv),
        ! which are proportional to the square root of TKE at the
        ! intermediate time step

        pcftke(jl,jk) = zsm*SQRT(ptkevn(jl,jk))
        pcfthv(jl,jk) = zsh*SQRT(ptkevn(jl,jk))

        ! Air density at half levels, and the prefactor that will be multiplied
        ! later to the exchange coeffcients to build a linear algebraic equation set.

        ztvm = (ptvm1(jl,jk)+ptvm1(jl,jk+1))*0.5_wp   ! Tv at half level k+1/2
        prhoh (jl,jk) = paphm1(jl,jk+1)/(ztvm*rd)     ! air density
        pprfac(jl,jk) = prhoh (jl,jk)/zdgh(jl,jk)     ! air density/dz/g

361   END DO
372 END DO

  END SUBROUTINE atm_exchange_coeff
  !-------------
  !>
  !!
  SUBROUTINE sfc_exchange_coeff( kproma, kbdim, ksfc_type,          &! in
                               & idx_wtr, idx_ice, idx_lnd,              &! in
                               & lsfc_mom_flux, lsfc_heat_flux,          &! in
                               & pz0m, ptsfc,                            &! in
                               & pfrc, pghabl,                           &! in
                               & pocu, pocv, ppsfc,                      &! in
                               & pum1_b, pvm1_b,                         &! in
                               & ptm1_b, pgeom1_b,                       &! in
                               & pqm1_b, pqxm1_b,                        &! in
                               & pqsat_b, plh_b,                         &! in
                               & ptheta_b, pthetav_b,                    &! in
                               & pthetal_b, paclc_b,                     &! in
                               & pthvvar_b,                              &! in
                               & pqsat_sfc, pcpt_sfc,                    &! out
                               & pri_gbm, pri_sfc,                       &! out
                               & pcfm_gbm, pcfm_sfc,                     &! out
                               & pcfh_gbm, pcfh_sfc,                     &! out
                               & pcfv_sfc,                               &! out
                               & pcftke_sfc, pcfthv_sfc,                 &! out
                               & pprfac_sfc, prho_sfc,                   &! out
                               & ptkevn_sfc, pthvvar_sfc,                &! out
                               & pqshear_sfc, pustarm,                   &! out
                               & pch_sfc,                                &! out
!                               & pch_sfc, pchn_sfc, pcdn_sfc, pcfnc_sfc, &! out
                               & pbn_sfc, pbhn_sfc, pbm_sfc, pbh_sfc,    &! out
                               & paz0lh,                                 &! in, optional
                               & pcsat, pcair                            &! in, optional
                               & )

    INTEGER, INTENT(IN) :: kproma, kbdim
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr, idx_ice, idx_lnd

    LOGICAL, INTENT(IN) :: lsfc_mom_flux   !< switch on/off surface momentum flux
    LOGICAL, INTENT(IN) :: lsfc_heat_flux  !< switch on/off surface fluxes of
                                           !< sensible and latent heat

    REAL(wp),INTENT(IN) :: pz0m     (kbdim,ksfc_type) !< aerodynamic roughness length
    REAL(wp),INTENT(IN) :: ptsfc    (kbdim,ksfc_type) !< temp. at surface
    REAL(wp),INTENT(IN) :: pfrc     (kbdim,ksfc_type) !< fraction of the grid box occupied by
                                                      !< each surface type
    REAL(wp),INTENT(IN) :: pghabl   (kbdim)  !< geopotential of PBL top
    REAL(wp),INTENT(IN) :: pocu     (kbdim)  !< ocean surface velocity
    REAL(wp),INTENT(IN) :: pocv     (kbdim)  !< ocean surface velocity
    REAL(wp),INTENT(IN) :: ppsfc    (kbdim)  !< surface pressure

    ! "_b" denotes value at the bottom level (the klev-th full level)

    REAL(wp),INTENT(IN) :: pum1_b   (kbdim)  !< u-wind
    REAL(wp),INTENT(IN) :: pvm1_b   (kbdim)  !< v-wind
    REAL(wp),INTENT(IN) :: ptm1_b   (kbdim)  !< temperature
    REAL(wp),INTENT(IN) :: pgeom1_b (kbdim)  !< geopotential
    REAL(wp),INTENT(IN) :: pqm1_b   (kbdim)  !< specific humidity
    REAL(wp),INTENT(IN) :: pqxm1_b  (kbdim)  !< total concentration of hydrometeors
    REAL(wp),INTENT(IN) :: pqsat_b  (kbdim)  !< saturation specific humidity
    REAL(wp),INTENT(IN) :: plh_b    (kbdim)  !< latent heat
    REAL(wp),INTENT(IN) :: ptheta_b (kbdim)  !< potential temp.
    REAL(wp),INTENT(IN) :: pthetav_b(kbdim)  !< virtual potential temp.
    REAL(wp),INTENT(IN) :: pthetal_b(kbdim)  !< liquid water (?) pot. temp.
    REAL(wp),INTENT(IN) :: paclc_b  (kbdim)  !< cloud cover at lowest model level

    ! For the variance of theta_v, "_b" denotes the lowest computational level
    ! above surface, i.e., the interface between full levels klev-1 and klev.

    REAL(wp),INTENT(IN) :: pthvvar_b (kbdim)  !< variance of theta_v

    ! "_sfc" denotes value at surface

    REAL(wp),INTENT(OUT) :: pqsat_sfc (kbdim,ksfc_type) !< saturation specific humidity at surface
    REAL(wp),INTENT(OUT) :: pcpt_sfc  (kbdim,ksfc_type) !< dry static energy
    REAL(wp),INTENT(OUT) :: pri_gbm   (kbdim)           !< moist Richardson number
    REAL(wp),INTENT(OUT) :: pri_sfc   (kbdim,ksfc_type) !< moist Richardson number

    REAL(wp),INTENT(OUT) :: pcfm_gbm  (kbdim)           !< exchange coeff. of momentum
    REAL(wp),INTENT(OUT) :: pcfm_sfc  (kbdim,ksfc_type) !< exchange coeff. of momentum,
                                                        !< for each type of surface
    REAL(wp),INTENT(OUT) :: pcfh_gbm  (kbdim)           !< exchange coeff. of heat and vapor
    REAL(wp),INTENT(OUT) :: pcfh_sfc  (kbdim,ksfc_type) !< exchange coeff. of heat and vapour
                                                        !< for each type of surface
    REAL(wp),INTENT(OUT) :: pcfv_sfc   (kbdim)  !< exchange coeff. of total water variance
    REAL(wp),INTENT(OUT) :: pcftke_sfc (kbdim)  !< exchange coeff. of TKE
    REAL(wp),INTENT(OUT) :: pcfthv_sfc (kbdim)  !< exchange coeff. of the variance of theta_v
    REAL(wp),INTENT(OUT) :: pprfac_sfc (kbdim)  !< prefactor for exchange coefficients
    REAL(wp),INTENT(OUT) :: prho_sfc   (kbdim)  !< air density
    REAL(wp),INTENT(OUT) :: ptkevn_sfc (kbdim)  !< boundary condition (sfc value) of TKE
    REAL(wp),INTENT(OUT) :: pthvvar_sfc(kbdim)  !< boundary condition (sfc value)
                                                !< of the variance of theta_v
    REAL(wp),INTENT(OUT) :: pqshear_sfc(kbdim)  !< vertical shear of total water concentration
    REAL(wp),INTENT(OUT) :: pustarm    (kbdim)  !< friction velocity, grid-box mean
    REAL(wp),INTENT(OUT) :: pbn_sfc  (kbdim,ksfc_type) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pbhn_sfc (kbdim,ksfc_type) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pbm_sfc  (kbdim,ksfc_type) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pbh_sfc  (kbdim,ksfc_type) !< for diagnostics
!    REAL(wp),INTENT(OUT) :: pchn_sfc (kbdim,ksfc_type) !<
!    REAL(wp),INTENT(OUT) :: pcdn_sfc (kbdim,ksfc_type) !<
!    REAL(wp),INTENT(OUT) :: pcfnc_sfc(kbdim,ksfc_type) !<
    REAL(wp),INTENT(OUT) :: pch_sfc  (kbdim,ksfc_type) !< for TKE boundary condition
!
! optional arguments for use with jsbach
    REAL(wp),OPTIONAL,INTENT(IN) :: paz0lh (kbdim)  !< roughness length for heat over land
    REAL(wp),OPTIONAL,INTENT(IN) :: pcsat  (kbdim)  !< area fraction with wet land surface
    REAL(wp),OPTIONAL,INTENT(IN) :: pcair  (kbdim)  !< area fraction with wet land surface (air)

!    REAL(wp),INTENT(INOUT) :: ptkem1_sfc (kbdim)  !< boundary condition (surface value) for TKE
!    REAL(wp),INTENT(INOUT) :: ptkem0_sfc (kbdim)  !< boundary condition (surface value) for TKE

!    REAL(wp) :: pbn_sfc  (kbdim,ksfc_type) !< for diagnostics
!    REAL(wp) :: pbhn_sfc (kbdim,ksfc_type) !< for diagnostics
!    REAL(wp) :: pbm_sfc  (kbdim,ksfc_type) !< for diagnostics
!    REAL(wp) :: pbh_sfc  (kbdim,ksfc_type) !< for diagnostics
    REAL(wp) :: pchn_sfc (kbdim,ksfc_type) !<
    REAL(wp) :: pcdn_sfc (kbdim,ksfc_type) !<
    REAL(wp) :: pcfnc_sfc(kbdim,ksfc_type) !<

    ! Local variables

    REAL(wp) :: zdu2   (kbdim,ksfc_type) !<
    REAL(wp) :: zcfnch (kbdim,ksfc_type) !<
    REAL(wp) :: zcr    (kbdim)           !< for open water only
    REAL(wp) :: zwst   (kbdim,ksfc_type) !<
    REAL(wp) :: zcsat  (kbdim,ksfc_type) !<
    REAL(wp) :: zustar (kbdim,ksfc_type) !< friction velocity
    REAL(wp) :: ztvsfc (kbdim)           !< virtual temperature at surface
    REAL(wp) :: zqts   (kbdim,ksfc_type)
    INTEGER  :: loidx  (kbdim,ksfc_type) !< counter for masks
    INTEGER  :: is     (ksfc_type)       !< counter for masks

    REAL(wp) :: zrdrv, zrvrd, zrgam, zonethird, ztwothirds
    REAL(wp) :: z2b, z3b, z3bc, zcons17, zepsr, zepdu2, zepsec
    REAL(wp) :: zqtl, zqtmit, zdqt, zqsmit
    REAL(wp) :: ztmit, ztheta, zthetav, zthetamit, zthetavmit, zfux, zfox
    REAL(wp) :: zmult1, zmult2, zmult3, zmult4, zmult5
    REAL(wp) :: zdus1, zdus2, zbuoy, zalo, zaloh, ztkev, zstabf
    REAL(wp) :: zdthetal, zdthv, z0h, zscf,  zconvs, zucf, zucfh
    REAL(wp) :: z2m, zcdn2m, zcdnr, zcfm2m, zust, zrrd
    LOGICAL  :: lhighz0
    INTEGER  :: jsfc, jl, jls, js

    !-------------------
    ! Some constants
    !-------------------
    zepsec     = 1.E-2_wp
    zrrd       = 1._wp/rd
    zrvrd      = rv/rd
    zrdrv      = rd/rv
    zrgam      = 1._wp/cgam
    z2b        = 2._wp*cb
    z3b        = 3._wp*cb
    z3bc       = 3._wp*cb*cc
    zcons17    = 1._wp / ckap**2
    zepsr      = 1.E-10_wp
    zepdu2     = 1._wp
    zonethird  = 1._wp/3._wp
    ztwothirds = 2._wp/3._wp

    !------------------------------------------------------------------------------
    ! The prefactor (= air density * Rd) that will be multiplied to the exchange
    ! coefficients when building the linear algebraic equations. The density here is
    ! computed using air temperature of the lowest model level at time step n-1.
    !------------------------------------------------------------------------------
    pprfac_sfc(1:kproma) =  ppsfc(1:kproma)                                     &
                         & /( ptm1_b(1:kproma)                                  &
                         &   *(1._wp+vtmpc1*pqm1_b(1:kproma)-pqxm1_b(1:kproma)) )

    !-------------------------------------------------------------
    ! COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
    ! RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
    ! AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
    ! COMMON PART OF THE DRAG COEFFICIENTS.
    !-------------------------------------------------------------
!TODO:    preset values to zero
     pcpt_sfc(1:kproma,1:ksfc_type) = 0._wp
     pqsat_sfc(1:kproma,1:ksfc_type) = 0._wp
     pri_sfc(1:kproma,1:ksfc_type) = 0._wp

     DO jsfc = 1,ksfc_type

! check for masks
!
      is(jsfc) = 0
      DO jl = 1,kproma
        IF(pfrc(jl,jsfc).GT.0.0_wp) THEN
          is(jsfc) = is(jsfc) + 1
          loidx(is(jsfc),jsfc) = jl
        ENDIF
      ENDDO

      CALL compute_qsat( kproma, is(jsfc), loidx(1,jsfc), ppsfc, &
                              ptsfc(1,jsfc), pqsat_sfc(1,jsfc) )

! loop over mask only
!
      DO jls = 1,is(jsfc)
        ! set index
        js=loidx(jls,jsfc)
        ! dry static energy pcpt_sfc
        !
        IF(jsfc == idx_lnd) THEN
          zqts(js,jsfc) = pcsat(js) * pqsat_sfc(js,jsfc) + (1._wp - pcair(js))*pqm1_b(js) ! q_total at land surface
        ELSE
          zqts(js,jsfc) = pqsat_sfc(js,jsfc)                                              ! q_total at surface
        END IF
        pcpt_sfc(js,jsfc) = ptsfc(js,jsfc) * (cpd + (cpv - cpd) * zqts(js,jsfc))

        ztheta      = ptsfc(js,jsfc)*(p0ref/ppsfc(js))**rd_o_cpd
        zthetav     = ztheta*(1._wp+vtmpc1*zqts(js,jsfc))

        zqtl       = pqm1_b(js) + pqxm1_b(js)                          ! q_total at lowest model level

        zqtmit     = 0.5_wp*( zqtl + zqts(js,jsfc) )                   ! q_total, vertical average

        zqsmit     = 0.5_wp*( pqsat_b  (js) + pqsat_sfc(js,jsfc) )  ! qs
        ztmit      = 0.5_wp*( ptm1_b   (js) + ptsfc    (js,jsfc) )  ! temp.
        zthetamit  = 0.5_wp*( ptheta_b (js) + ztheta  )  ! potential temp.
        zthetavmit = 0.5_wp*( pthetav_b(js) + zthetav )  ! virtual potential temp.

        zfux = plh_b(js)/(cpd*ztmit)
        zfox = plh_b(js)/(rd*ztmit)

        zmult1 = 1._wp+vtmpc1*zqtmit   ! A in clear sky
        zmult2 = zfux*zmult1-zrvrd
        zmult3 = zrdrv*zfox*zqsmit/(1._wp+zrdrv*zfox*zfux*zqsmit)
        zmult5 = zmult1-zmult2*zmult3  ! A in cloud
        zmult4 = zfux*zmult5-1._wp     ! D in cloud

        zdus1 = paclc_b(js)*zmult5 + (1._wp-paclc_b(js))*zmult1   ! A avg
        zdus2 = paclc_b(js)*zmult4 + (1._wp-paclc_b(js))*vtmpc1   ! D avg

        zdqt     = zqtl - zqts(js,jsfc)                           ! d qt
        zdthetal = pthetal_b(js) - ztheta                         ! d theta_l
        IF (jsfc == idx_lnd) THEN                                 ! over land
          zdu2(js,jsfc) = MAX(zepdu2,(pum1_b(js)**2+pvm1_b(js)**2))
        ELSE                                                      ! over water or ice
          zdu2(js,jsfc) = MAX(zepdu2,(pum1_b(js)-pocu(js))**2 &   ! (d u)^2
                                    +(pvm1_b(js)-pocv(js))**2)    ! (d v)^2
        END IF

        zbuoy        = zdus1*zdthetal + zdus2*zthetamit*zdqt
        pri_sfc(js,jsfc) = pgeom1_b(js)*zbuoy/(zthetavmit*zdu2(js,jsfc))

        zalo = LOG( 1._wp + pgeom1_b(js)/(grav*pz0m(js,jsfc)) )  ! ln[ 1 + zL/z0m ]
        pcdn_sfc(js,jsfc) = (ckap/zalo)**2

        pcfnc_sfc(js,jsfc)= SQRT(zdu2(js,jsfc))*pcdn_sfc(js,jsfc)
 !  REMARK:
 !  compared to the original (precalc_land,ocean,ice) the factor zcons is missing
 !  this factor is included as "prefactor for the exchange coefficients" in
 !  subroutine matrix_setup_elim

        zdthv         = MAX(0._wp,(zthetav-pthetav_b(js)))
        zwst(js,jsfc) = zdthv*SQRT(zdu2(js,jsfc))/zthetavmit

        IF ( jsfc == idx_wtr ) THEN         ! over water
          z0h        =pz0m(js,jsfc)*EXP(2._wp-86.276_wp*pz0m(js,jsfc)**0.375_wp)
          zaloh      =LOG(1._wp+pgeom1_b(js)/(grav*z0h))
          pchn_sfc  (js,jsfc)=ckap**2/(zalo*zaloh)
          zcfnch(js,jsfc)=SQRT(zdu2(js,jsfc))*pchn_sfc(js,jsfc)
          zcr   (js)=(cfreec/(pchn_sfc(js,jsfc)*SQRT(zdu2(js,jsfc))))*ABS(zbuoy)**zonethird
        ELSE IF ( jsfc == idx_ice ) THEN    ! over ice
          pchn_sfc  (js,jsfc) = pcdn_sfc (js,jsfc)
          zcfnch(js,jsfc) = pcfnc_sfc(js,jsfc)  ! coeff. for scalar is the same as for momentum
        ELSE IF ( jsfc == idx_lnd ) THEN    ! over  land
          pchn_sfc  (js,jsfc) = (ckap / LOG(1._wp+pgeom1_b(js)/(grav*paz0lh(js))))**2
          zcfnch(js,jsfc) = SQRT(zdu2(js,jsfc))*pchn_sfc(js,jsfc)
        ENDIF
      ENDDO ! 1:is
    ENDDO   ! 1:ksfc_type

    !-------------------------------------------------------------------------
    ! Compute the exchange coefficients for momentum, heat and water vapour,
    ! for each type of surface
    !-------------------------------------------------------------------------

    IF (lsfc_mom_flux.OR.lsfc_heat_flux) THEN  ! Surface flux is considered

!TODO:    preset values to zero
     pcfh_sfc(1:kproma,1:ksfc_type) = 0._wp
     pcfm_sfc(1:kproma,1:ksfc_type) = 0._wp

      ! stable case              pri_sfc(js,jsfc) > 0.

      DO jsfc = 1,ksfc_type
        DO jls = 1,is(jsfc)
          ! set index
          js=loidx(jls,jsfc)
          IF ( pri_sfc(js,jsfc) > 0._wp ) THEN
            zscf = SQRT(1._wp+pri_sfc(js,jsfc))
            pcfm_sfc(js,jsfc) = pcfnc_sfc (js,jsfc)/(1._wp+z2b*pri_sfc(js,jsfc)/zscf)   ! 5.2? 5.5
            pcfh_sfc(js,jsfc) = zcfnch(js,jsfc)/(1._wp+z2b*pri_sfc(js,jsfc)*zscf)   ! 5.2? 5.6
            pch_sfc (js,jsfc) = pchn_sfc  (js,jsfc)/(1._wp+z2b*pri_sfc(js,jsfc)*zscf)   ! for (5.22)
          ENDIF
        ENDDO
      ENDDO

      ! unstable case           pri_sfc(js,jsfc) <= 0.

      IF (idx_wtr<=ksfc_type) THEN
        jsfc = idx_wtr  ! water
        DO jls = 1,is(jsfc)
          ! set index
          js=loidx(jls,jsfc)
          IF ( pri_sfc(js,jsfc) <= 0._wp ) THEN
            zucf =  SQRT( -pri_sfc(js,jsfc)*(1._wp+ pgeom1_b(js)/(grav*pz0m(js,jsfc))) ) ! sqrt in (5.4)
            zucf =  1._wp+z3bc*pcdn_sfc(js,jsfc)*zucf                          ! denominator in (5.4)
            zucf =  1._wp/zucf
            pcfm_sfc(js,jsfc) = pcfnc_sfc (js,jsfc)*(1._wp-z2b*pri_sfc(js,jsfc)*zucf)  ! (5.2), (5.4)
            pcfh_sfc(js,jsfc) = zcfnch(js,jsfc)*(1._wp+zcr(js)**cgam)**zrgam       ! (5.9)
            pch_sfc (js,jsfc) = pchn_sfc  (js,jsfc)*(1._wp+zcr(js)**cgam)**zrgam
          ENDIF
        ENDDO
      ENDIF

      IF (idx_ice<=ksfc_type) THEN
        jsfc = idx_ice  ! ice
        DO jls = 1,is(jsfc)
          ! set index
          js=loidx(jls,jsfc)
          IF ( pri_sfc(js,jsfc) <= 0._wp ) THEN
            zucf =  SQRT( -pri_sfc(js,jsfc)*(1._wp+ pgeom1_b(js)/(grav*pz0m(js,jsfc))) ) ! sqrt in (5.4)
            zucf =  1._wp+z3bc*pcdn_sfc(js,jsfc)*zucf                   ! denominator in (5.4)
            zucf =  1._wp/zucf
            pcfm_sfc(js,jsfc) = pcfnc_sfc (js,jsfc)*(1._wp-z2b*pri_sfc(js,jsfc)*zucf)  ! (5.2), (5.4)
            pcfh_sfc(js,jsfc) = zcfnch(js,jsfc)*(1._wp-z3b*pri_sfc(js,jsfc)*zucf)  ! (5.2), (5.4)
            pch_sfc (js,jsfc) = pchn_sfc  (js,jsfc)*(1._wp-z3b*pri_sfc(js,jsfc)*zucf)
          ENDIF
        ENDDO
      ENDIF

      IF (idx_lnd<=ksfc_type) THEN
        jsfc = idx_lnd  ! land
        DO jls = 1,is(jsfc)
          ! set index
          js=loidx(jls,jsfc)
          IF ( pri_sfc(js,jsfc) <= 0._wp ) THEN
            zucf =  SQRT( -pri_sfc(js,jsfc)*(1._wp+ pgeom1_b(js)/(grav*pz0m(js,jsfc))) ) ! sqrt in (5.4)
            zucf =  1._wp+z3bc*pcdn_sfc(js,jsfc)*zucf                   ! denominator in (5.4)
            zucf =  1._wp/zucf
            zucfh=  SQRT( -pri_sfc(js,jsfc)*(1._wp+ pgeom1_b(js)/(grav*paz0lh(js))) ) ! sqrt in (5.4)
            zucfh=  1._wp+z3bc*pchn_sfc(js,jsfc)*zucfh                   ! denominator in (5.4)
            zucfh=  1._wp/zucfh
            pcfm_sfc(js,jsfc) = pcfnc_sfc (js,jsfc)*(1._wp-z2b*pri_sfc(js,jsfc)*zucf)  ! (5.2), (5.4)
            pcfh_sfc(js,jsfc) = zcfnch(js,jsfc)*(1._wp-z3b*pri_sfc(js,jsfc)*zucfh) ! (5.2), (5.4)
            pch_sfc (js,jsfc) = pchn_sfc  (js,jsfc)*(1._wp-z3b*pri_sfc(js,jsfc)*zucfh)
          ENDIF
        ENDDO
      ENDIF

    END IF  ! lsfc_mom_flux.OR.lsfc_heat_flux

    IF (.NOT.lsfc_mom_flux) THEN  ! Surface momentum flux is switched off
      pcfm_sfc(1:kproma,1:ksfc_type) = 0._wp
    END IF

    IF (.NOT.lsfc_heat_flux) THEN  ! Surface heat flux is switched off
      pcfh_sfc(1:kproma,1:ksfc_type) = 0._wp
      zwst    (1:kproma,1:ksfc_type) = 0._wp   ! affects TKE at surface
    END IF

    !------------------------------------------------------------------------------
    ! Store values
    ! to be used in subroutine "nsurf_diag" to compute
    ! new t2m, 2m dew point, 10m wind components
    !------------------------------------------------------------------------------

    pbn_sfc (1:kproma,1:ksfc_type) = 0._wp  !
    pbhn_sfc(1:kproma,1:ksfc_type) = 0._wp  !
    pbm_sfc (1:kproma,1:ksfc_type) = 0._wp  !
    pbh_sfc (1:kproma,1:ksfc_type) = 0._wp  !
    DO jsfc = 1,ksfc_type
      DO jls = 1,is(jsfc)
! set index
      js=loidx(jls,jsfc)
        pbn_sfc(js,jsfc) = ckap / SQRT(pcdn_sfc(js,jsfc))
        pbhn_sfc(js,jsfc) = ckap / SQRT(pchn_sfc(js,jsfc))
        pbm_sfc(js,jsfc) = MAX(zepsec, SQRT(pcfm_sfc(js,jsfc)*pcdn_sfc(js,jsfc) *  &
                               zcons17 / pcfnc_sfc(js,jsfc)))
        pbh_sfc(js,jsfc) = MAX(zepsec, pch_sfc(js,jsfc)/pbm_sfc(js,jsfc)*zcons17)
        pbm_sfc(js,jsfc) = 1._wp / pbm_sfc(js,jsfc)
        pbh_sfc(js,jsfc) = 1._wp / pbh_sfc(js,jsfc)
      END DO
    END DO

    !-------------------------------------------------------------------------
    ! Get the aggregated exchange coefficient for momentum
    !-------------------------------------------------------------------------
    ! For sensible heat and water vapour, it is not the exchange coefficient
    ! but the solution at the lowest model level that is aggregated.
    ! The exchange coeffcients (cfh_sfc) computed in this subroutine
    ! are returned to the calling subroutine for each surface type.
    ! They are used later for solving the discretized vertical diffusion
    ! equation at the lowest grid level (klev) for each surface type
    ! separately. Then the solutions are aggregated using
    ! (fraction of type)*(cfh_sfc of type) as the weighting factor, which
    ! ensures conservation of the area-weighted total flux.
    !   Here we compute the aggregated exchange coefficient for output
    ! and for solving the vertical diffusion equation of the variance of
    ! virtual potential temperature (theta_v).
    !-------------------------------------------------------------------------
    ! Add aggregated Richardson number for the surface
    !-------------------------------------------------------------------------

    pcfm_gbm(1:kproma) = 0._wp
    pcfh_gbm(1:kproma) = 0._wp
    pri_gbm (1:kproma) = 0._wp
    DO jsfc = 1,ksfc_type
      DO jls = 1,is(jsfc)
! set index
      js=loidx(jls,jsfc)
        pcfm_gbm(js) = pcfm_gbm(js) + pfrc(js,jsfc)*pcfm_sfc(js,jsfc)
        pcfh_gbm(js) = pcfh_gbm(js) + pfrc(js,jsfc)*pcfh_sfc(js,jsfc)
        pri_gbm (js) = pri_gbm (js) + pfrc(js,jsfc)*pri_sfc (js,jsfc)
      ENDDO
    ENDDO

    !-------------------------------------------------------------------------
    ! Hydrometeors and the other tracers share the same exchange coefficient
    ! with heat and moisture, but have no turbulence-induced surface flux.
    ! These are taken care of in subroutine matrix_setup.
    !   The total water variance has a different exchange coefficient
    ! (variable cfv), and no surface flux. Set the surface exchange coefficient
    ! to zero.
    !-------------------------------------------------------------------------
    pcfv_sfc(1:kproma) = 0._wp

    !-------------------------------------------------------------------------
    ! Compute vertical shear of total water
    !-------------------------------------------------------------------------
    zcsat(1:kproma,1:ksfc_type) = 1._wp
    IF (idx_lnd<=ksfc_type) zcsat(1:kproma,idx_lnd)     = pcsat(1:kproma)

    pqshear_sfc(1:kproma) = 0._wp ! initialization for weighted q_total at surface

    DO jsfc = 1,ksfc_type
      DO jls = 1,is(jsfc)
! set index
      js=loidx(jls,jsfc)
        pqshear_sfc(js) =  pqshear_sfc(js)      &
                              & +  pqsat_sfc(js,jsfc) &
                              &       *zcsat(js,jsfc) &
                              &        *pfrc(js,jsfc)
      ENDDO
    ENDDO

    pqshear_sfc(1:kproma) = ( pqm1_b(1:kproma)+pqxm1_b(1:kproma)          &
                          &  -pqshear_sfc(1:kproma) )*grav/pgeom1_b(1:kproma)

    !-------------------------------------------------------------------------
    ! Diagnose friction velocity. The values of each individual surface type
    ! are used immediately below for computing the surface value of TKE;
    ! The grid-box mean is used in the next time step in subroutine
    ! "atm_exchange_coeff" for finding the PBL height.
    !-------------------------------------------------------------------------
    IF (lsfc_mom_flux) THEN  ! Surface momentum flux is switched on

      z2m = 2._wp            ! 2-m height
      DO jsfc = 1,ksfc_type
        DO jls = 1,is(jsfc)
! set index
        js=loidx(jls,jsfc)
          lhighz0 = pz0m(js,jsfc).GT.z2m
          zcdn2m = MERGE((ckap/LOG(1._wp+pgeom1_b(js)/(grav*z2m)))**2, &
                 &        pcdn_sfc(js,jsfc),lhighz0 )
          zcdnr  = zcdn2m/pcdn_sfc(js,jsfc)

          zucf   = SQRT( ABS(pri_sfc(js,jsfc))*(1._wp+pgeom1_b(js)/(grav*z2m)) )  ! sqrt in (5.4)
          zucf  = 1._wp + z3bc*zcdn2m*zucf                              ! denomenator in (5.4)
          zucf = 1._wp - z2b*pri_sfc(js,jsfc)/zucf                      ! (5.4)
          zcfm2m = MERGE( pcfnc_sfc(js,jsfc)*zcdnr*zucf, pcfm_sfc(js,jsfc)*zcdnr, &
                 &        lhighz0.AND.pri_sfc(js,jsfc).LT.0._wp           )
          zust   = zcfm2m*SQRT(zdu2(js,jsfc))
          zustar(js,jsfc) = SQRT(zust)
        END DO
      END DO
 !  REMARK:
 !  compared to the original (precalc_land,ocean,ice) the factor 1/zcons is missing
 !  this factor is included as "prefactor for the exchange coefficients" in
 !  subroutine matrix_setup_elim

      pustarm(1:kproma) = 0._wp
      DO jsfc = 1,ksfc_type
        DO jls = 1,is(jsfc)
! set index
        js=loidx(jls,jsfc)
          pustarm(js) = pustarm(js) + pfrc(js,jsfc)*zustar(js,jsfc)
        END DO
      END DO

    ELSE ! Surface momentum flux is off. Friction velocity is by definition zero.
      zustar (1:kproma,1:ksfc_type) = 0._wp
      pustarm(1:kproma)             = 0._wp
    END IF

    !---------------------------------------------------
    ! Surface value of TKE and its exchange coefficient
    !---------------------------------------------------
    ! The exchange coefficient of TKE is set to the same value as for momentum
    pcftke_sfc(1:kproma) = pcfm_gbm(1:kproma)

    ! TKE at the surface (formulation of Mailhot and Benoit (1982))

    ptkevn_sfc(1:kproma) = 0._wp  ! initialize the weighted average

    DO jsfc = 1,ksfc_type
      DO jls = 1,is(jsfc)
! set index
      js=loidx(jls,jsfc)
        ztkev = custf*(zustar(js,jsfc)**2)
        IF(zwst(js,jsfc).GT.zepsr) THEN
           zconvs = (zwst(js,jsfc)*pch_sfc(js,jsfc)*pghabl(js))**zonethird
           zstabf = (pgeom1_b(js)*ckap*zwst(js,jsfc)*pch_sfc(js,jsfc))**ztwothirds
           zstabf = MIN(custf*3._wp*zustar(js,jsfc)**2,zstabf)
           ztkev = ztkev + zstabf + cwstf*(zconvs**2)
        END IF
        ptkevn_sfc(js) = ptkevn_sfc(js) + ztkev*pfrc(js,jsfc)
      END DO
    END DO
    ptkevn_sfc(1:kproma) = MAX( tke_min,ptkevn_sfc(1:kproma) )

    !----------------------------------------------------------------
    ! Surface value and exchange coefficient of theta_v variance
    !----------------------------------------------------------------
    ! The exchange coefficient is set to the aggregated coefficient
    ! of heat and moisture.
    pcfthv_sfc(1:kproma) = pcfh_gbm(1:kproma)

    ! thvvar at the surface
    pthvvar_sfc(1:kproma) = pthvvar_b(1:kproma)

    !------------------------------------------------------------------------------
    ! Compute the surface air density using surface temperature and humidity,
    ! to be used in subroutine "vdiff_tendencies" for diagnosing the
    ! turbulence-induced production of total water variance.
    !------------------------------------------------------------------------------

    prho_sfc(1:kproma) = 0._wp  ! Initialize the area weighted average
    DO jsfc = 1,ksfc_type
      DO jls = 1,is(jsfc)
        js=loidx(jls,jsfc) ! set index

        ztvsfc(js) = ptsfc(js,jsfc)*(1._wp + vtmpc1*zqts(js,jsfc))
        prho_sfc(js) = prho_sfc(js) + zrrd*ppsfc(js)/ztvsfc(js)
      END DO
    END DO
! extra variable for land points:
!    IF (idx_lnd<=ksfc_type) THEN
!     jsfc = idx_lnd  ! land
!      DO jls = 1,is(jsfc)
!! set index
!      js=loidx(jls,jsfc)
!        pbhn_sfc(js,jsfc) = ckap / SQRT(pchn_sfc(js,jsfc))
!      END DO
!    END IF

  END SUBROUTINE sfc_exchange_coeff
  !-------------

END MODULE mo_turbulence_diag
