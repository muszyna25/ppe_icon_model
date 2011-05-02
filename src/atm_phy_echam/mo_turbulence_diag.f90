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
MODULE mo_turbulence_diag

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
  USE mo_convect_tables,    ONLY: tlucua, jptlucu1, jptlucu2,  &
                                & lookuperror, lookupoverflow, &
                                & compute_qsat

#ifdef __ICON__
!  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_physical_constants,ONLY: g=>grav, rd, cpd, rd_o_cpd, rv,           &
                                & vtmpc1,vtmpc2,tmelt,alv,als
  USE mo_echam_vdiff_params,ONLY: clam, cgam, ckap, cb,cc, chneu, shn, smn, &
                                & da1, custf, cwstf, cfreec,                &
                                & epshr=>eps_shear, epcor=>eps_corio,       &
                                & tkemin=>tke_min, cons2, cons25, cons5
#else
  USE mo_constants,ONLY: g, rd, cpd, rd_o_cpd, rv,                 &
                       & vtmpc1,vtmpc2,tmelt,alv,als
  USE mo_physc2,   ONLY: clam, cgam, ckap, cb,cc, chneu, shn, smn, &
                       & da1, custf, cwstf, cfreec,                &
                       & epshr, epcor, tkemin, cons2, cons25, cons5
  USE mo_time_control, ONLY: lstart
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: atm_exchange_coeff, sfc_exchange_coeff

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

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
  !!
  SUBROUTINE atm_exchange_coeff( kproma, kbdim, klev, klevm1, klevp1,     &! in
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

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, klevp1
    REAL(wp),INTENT(IN) :: pstep_len
    REAL(wp),INTENT(IN) :: pcoriol(kbdim)   !< Coriolis parameter: 2*omega*sin(lat)
    REAL(wp),INTENT(IN) :: pum1(kbdim,klev),   pvm1(kbdim,klev)
    REAL(wp),INTENT(IN) :: ptm1(kbdim,klev),   ptvm1(kbdim,klev)
    REAL(wp),INTENT(IN) :: pgeom1(kbdim,klev), pqm1(kbdim,klev)
    REAL(wp),INTENT(IN) :: pxm1(kbdim,klev)
    REAL(wp),INTENT(IN) :: papm1(kbdim,klev),  paphm1(kbdim,klevp1)
    REAL(wp),INTENT(IN) :: paclc(kbdim,klev)

    REAL(wp),INTENT(IN) :: pthvvar(kbdim,klev)

#ifdef __ICON__
    REAL(wp),INTENT(IN) :: pustarm(kbdim)
    REAL(wp),INTENT(IN) :: ptkem1 (kbdim,klev)
#else
    REAL(wp),INTENT(INOUT) :: pustarm(kbdim)
    REAL(wp),INTENT(INOUT) :: ptkem1 (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: ptkem0 (kbdim,klev)
#endif

    REAL(wp),INTENT(OUT) :: pcptgz(kbdim,klev)   !< dry static energy
    INTEGER, INTENT(OUT) :: ihpbl(kbdim)         !< grid level index of PBL top
    REAL(wp),INTENT(OUT) :: pghabl(kbdim)        !< geopotential of PBL top

    REAL(wp),INTENT(OUT) :: pqshear (kbdim,klevm1) !< vertical gradient of qv
    REAL(wp),INTENT(OUT) :: pzthvvar(kbdim,klevm1) !< variance of theta_v at intermediate step
    REAL(wp),INTENT(OUT) :: ptkevn  (kbdim,klevm1) !< TKE at intermediate time step

    REAL(wp),INTENT(OUT) :: pcfm    (kbdim,klevm1) !< exchange coeff. for u, v
    REAL(wp),INTENT(OUT) :: pcfh    (kbdim,klevm1) !< exchange coeff. for cptgz and tracers
    REAL(wp),INTENT(OUT) :: pcfv    (kbdim,klevm1) !< exchange coeff. for variance of qx
    REAL(wp),INTENT(OUT) :: pcftke  (kbdim,klevm1) !< exchange coeff. for TKE
    REAL(wp),INTENT(OUT) :: pcfthv  (kbdim,klevm1) !< exchange coeff. for variance of theta_v
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

    REAL(wp) :: ztheta (kbdim,klev)  !< potential temperature
    REAL(wp) :: zthetav(kbdim,klev)  !< virtual potential temperature
    REAL(wp) :: zthetal(kbdim,klev)  !< liquid (and ice?) potential temperature
    REAL(wp) :: zqsat  (kbdim,klev)  !< specific humidity at saturation
    REAL(wp) :: zlh    (kbdim,klev)  !< latent heat at full levels

    ! - Variables defined at half levels

    REAL(wp) :: zlhh (kbdim,klevm1)   !< latent heat at half levels
    REAL(wp) :: zdgh (kbdim,klevm1)   !< geopotential difference between two full levels

    REAL(wp) :: zqsatm  (kbdim,klevm1)
    REAL(wp) :: ztmitte (kbdim,klevm1)
    REAL(wp) :: zqmit   (kbdim,klevm1)
    REAL(wp) :: zthetavh(kbdim,klevm1)
    REAL(wp) :: zthetah (kbdim,klevm1)
    REAL(wp) :: zccover (kbdim,klevm1)
    REAL(wp) :: zqxmit  (kbdim,klevm1)

    ! - 1D variables and scalars

    REAL(wp) :: zhdyn (kbdim)
    INTEGER  :: ihpblc(kbdim), ihpbld(kbdim)
    REAL(wp) :: zrvrd, zrdrv, zonethird
    REAL(wp) :: zusus1, zes, zsdep1, zsdep2, zcor, zcons23
    REAL(wp) :: zrdp, zds, zdz, zri, zqtmit, zfux, zfox
    REAL(wp) :: zmult1, zmult2, zmult3, zmult4, zmult5
    REAL(wp) :: zdus1, zdus2, zteldif, zthvirdif, zdqtot, zqddif
    REAL(wp) :: zbuoy, zdusq, zdvsq, zshear, zhexp, zlam, ztvm
    REAL(wp) :: z2geomf, zz2geo, zmix, zalh2, zucf, zsh, zsm, zzb, zdisl
    REAL(wp) :: zktest, ztkesq, zthvprod, zthvdiss
    INTEGER  :: jk, jl, it

    !-------------------------------------
    ! 1. Some constants
    !-------------------------------------
    zrvrd     = rv/rd
    zrdrv     = rd/rv
    zonethird = 1._wp/3._wp

    !-------------------------------------
    ! 2. NEW THERMODYNAMIC VARIABLES
    !-------------------------------------
    lookupoverflow = .FALSE.

    DO 212 jk=1,klev
      DO 211 jl=1,kproma

        ! Virtual dry static energy, potential temperature, virtual potential
        ! temperature

        pcptgz (jl,jk) = pgeom1(jl,jk)+ptm1(jl,jk)*cpd*(1._wp+vtmpc2*pqm1(jl,jk))
        ztheta (jl,jk) = ptm1(jl,jk)*(100000._wp/papm1(jl,jk))**rd_o_cpd
        zthetav(jl,jk) = ztheta(jl,jk)*(1._wp+vtmpc1*pqm1(jl,jk)-pxm1(jl,jk))

        ! Latent heat, liquid (and ice) potential temperature

        zlh(jl,jk)     = MERGE(alv,als,ptm1(jl,jk).GE.tmelt)  ! latent heat
        zusus1         = zlh(jl,jk)/cpd*ztheta(jl,jk)/ptm1(jl,jk)*pxm1(jl,jk)
        zthetal(jl,jk) = ztheta(jl,jk)-zusus1

        it = NINT(ptm1(jl,jk)*1000._wp)
        IF (it<jptlucu1 .OR. it>jptlucu2) lookupoverflow = .TRUE.
        it = MAX(MIN(it,jptlucu2),jptlucu1)
        zes=tlucua(it)/papm1(jl,jk)          ! (sat. vapour pressure)*Rd/Rv/ps
        zes=MIN(zes,0.5_wp)
        zqsat(jl,jk)=zes/(1._wp-vtmpc1*zes)  ! specific humidity at saturation
211   END DO
212 END DO
    IF (lookupoverflow) CALL lookuperror ('vdiff (1)   ')

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

#ifdef __ICON__
    ! Initial value of ustar has been set in subroutine "init_phy_memory".
#else
    IF (lstart) THEN
      pustarm(1:kproma) = 1._wp
    END IF
#endif
    DO jl = 1,kproma
      zcor=MAX(ABS(pcoriol(jl)),epcor)
      zhdyn(jl)=MIN(pgeom1(jl,1)/g,chneu*pustarm(jl)/zcor)
      ihpblc(jl)=klev
      ihpbld(jl)=klev
    END DO

    DO jk=klevm1,1,-1
      DO jl=1,kproma
        zds=pcptgz(jl,jk)-pcptgz(jl,klev)
        zdz=pgeom1(jl,jk)/g-zhdyn(jl)
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

        zqtmit=zqxmit(jl,jk)+zqmit(jl,jk)       ! qt
        zfux=zlhh(jl,jk)/(cpd*ztmitte(jl,jk))   ! L/(cpd*T)
        zfox=zlhh(jl,jk)/(rd*ztmitte(jl,jk))    ! L/(Rd*T)
        zmult1=1._wp+vtmpc1*zqtmit              ! (1+0.61*qt) = A in clear sky
        zmult2=zfux*zmult1-zrvrd
        zmult3=zrdrv*zfox*zqsatm(jl,jk)             &
               /(1._wp+zrdrv*zfux*zfox*zqsatm(jl,jk))
        zmult5=zmult1-zmult2*zmult3             ! A in cloud
        zmult4=zfux*zmult5-1._wp                ! D in cloud
        zdus1=zccover(jl,jk)*zmult5+(1._wp-zccover(jl,jk))*zmult1       ! A avg
        zdus2=zccover(jl,jk)*zmult4+(1._wp-zccover(jl,jk))*vtmpc1       ! D avg
        zteldif  =(zthetal(jl,jk)-zthetal(jl,jk+1))/zdgh(jl,jk)*g       ! d theta_l
        zthvirdif=(zthetav(jl,jk)-zthetav(jl,jk+1))/zdgh(jl,jk)*g       ! d theta_v
        zdqtot=(pqm1(jl,jk)+pxm1(jl,jk))-(pqm1(jl,jk+1)+pxm1(jl,jk+1))  ! d qt
        zqddif=zdqtot/zdgh(jl,jk)*g                                     ! (d qt)/(d z)
        zbuoy  = (zteldif*zdus1+zthetah(jl,jk)*zdus2*zqddif)  &
                *g/zthetavh(jl,jk)
        zdusq  = (pum1(jl,jk)-pum1(jl,jk+1))**2
        zdvsq  = (pvm1(jl,jk)-pvm1(jl,jk+1))**2
        zshear = (zdusq+zdvsq)*(g/zdgh(jl,jk))**2
        zri    = zbuoy/MAX(zshear,epshr)

        pqshear(jl,jk) = zqddif   ! store for variance production
        pri(jl,jk) = zri          ! save for output

        ! ASYMPTOTIC MIXING LENGTH FOR MOMENTUM AND
        ! HEAT (ZLAM) ABOVE THE PBL AS A FUNCTION OF HEIGHT
        ! ACCORDING TO HOLTSLAG AND BOVILLE (1992), J. CLIMATE.

        zhexp=EXP(1._wp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))
        zlam=1._wp+(clam-1._wp)*zhexp
        IF(jk.GE.ihpbl(jl)) THEN
           zcons23=cons25        ! (1/lambda)*(0.5/grav)
        ELSE
           zcons23=cons2/zlam    ! (1/lambda)*(0.5/grav)
        END IF

        ! MIXING LENGTH (BLACKADAR)

        z2geomf = pgeom1(jl,jk)+pgeom1(jl,jk+1)  ! half-level value of gz
        zz2geo  = cons2*z2geomf                  ! z*(Karman constant)
        zmix    = zz2geo/(1._wp+zcons23*z2geomf) ! mixing length

        pmixlen(jl,jk) = zmix   ! save for output

        ! STABILITY FUNCTIONS (LOUIS, 1979)

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
           ptkevn(jl,jk)=tkemin
        ELSE
           ptkevn(jl,jk)=MAX(tkemin,(zdisl*(SQRT(zktest)-1._wp))**2)
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

        ztkesq = SQRT(MAX(tkemin,ptkem1(jl,jk)))

        ! Virtual potential temperautre variance at intermediate step,
        ! obtained by solving a prognostic equation the variance, the right-hand side
        ! of which consists of a production term and dissipation.
        ! An explicit time stepping method is used here.

        zthvprod = 2._wp*zsh*ztkesq*zthvirdif**2       ! production rate
        zthvdiss = pthvvar(jl,jk)*ztkesq/(da1*zmix)    ! dissipation rate
        pzthvvar(jl,jk) = pthvvar(jl,jk)+(zthvprod-zthvdiss)*pstep_len
        pzthvvar(jl,jk) = MAX(tkemin,pzthvvar(jl,jk))

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
  !!
  SUBROUTINE sfc_exchange_coeff( kproma, kbdim, ksfc_type,               &! in
                               & idx_wtr, idx_ice, idx_lnd,              &! in
                               & lsfc_mom_flux, lsfc_heat_flux,          &! in
                               & pz0m, ptsfc,                            &! in
                               & pfrc, pghabl, pocu, pocv,               &! in
                               & ppsfc, pum1_b, pvm1_b,                  &! in
                               & ptm1_b, pgeom1_b, pqm1_b, pqxm1_b,      &! in
                               & pqsat_b, plh_b, ptheta_b,               &! in
                               & pthetav_b, pthetal_b, paclc_b,          &! in
                               & pthvvar_b,                              &! in
#ifdef __ICON__
#else
                               & ptkem1_sfc, ptkem0_sfc,                 &! inout
#endif
                               & pqsat_sfc, pcpt_sfc, pri_sfc,           &! out
                               & pcfm_gbm, pcfm_sfc, pcfh_gbm, pcfh_sfc, &! out
                               & pcfv_sfc, pcftke_sfc, pcfthv_sfc,       &! out
                               & pprfac_sfc, prho_sfc,                   &! out
                               & ptkevn_sfc, pthvvar_sfc,                &! out
                               & pqshear_sfc, pustarm                    )! out

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

#ifdef __ICON__
#else
    REAL(wp),INTENT(INOUT) :: ptkem1_sfc (kbdim)  !< boundary condition (surface value) for TKE
    REAL(wp),INTENT(INOUT) :: ptkem0_sfc (kbdim)  !< boundary condition (surface value) for TKE
#endif

    ! Local variables

    REAL(wp) :: zdu2   (kbdim,ksfc_type) !<
    REAL(wp) :: zchn   (kbdim,ksfc_type) !<
    REAL(wp) :: zcfnch (kbdim,ksfc_type) !<
    REAL(wp) :: zcfnc  (kbdim,ksfc_type) !<
    REAL(wp) :: zcdn   (kbdim,ksfc_type) !<
    REAL(wp) :: zcr    (kbdim)           !< for open water only
    REAL(wp) :: zch    (kbdim,ksfc_type) !< for TKE boundary condition
    REAL(wp) :: zwst   (kbdim,ksfc_type) !<
    REAL(wp) :: zcsat  (kbdim,ksfc_type) !<
    REAL(wp) :: zustar (kbdim,ksfc_type) !< friction velocity
    REAL(wp) :: ztvsfc (kbdim)           !< virtual temperature at surface

    REAL(wp) :: zrdrv, zrvrd, zrgam, zonethird, ztwothirds
    REAL(wp) :: z2b, z3b, z3bc, zepsr, zepdu2
    REAL(wp) :: zqtl, zqts, zqtmit, zdqt, zqsmit
    REAL(wp) :: ztmit, ztheta, zthetav, zthetamit, zthetavmit, zfux, zfox
    REAL(wp) :: zmult1, zmult2, zmult3, zmult4, zmult5
    REAL(wp) :: zdus1, zdus2, zbuoy, zalo, zaloh, ztkev, zstabf
    REAL(wp) :: zdthetal, zdthv, z0h, zscf,  zconvs, zucf
    REAL(wp) :: z2m, zcdn2m, zcdnr, zcfm2m, zust, zrrd
    LOGICAL  :: lhighz0
    INTEGER  :: jsfc, jl

    !-------------------
    ! Some constants
    !-------------------
    zrrd       = 1._wp/rd
    zrvrd      = rv/rd
    zrdrv      = rd/rv
    zrgam      = 1._wp/cgam
    z2b        = 2._wp*cb
    z3b        = 3._wp*cb
    z3bc       = 3._wp*cb*cc
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

 
     DO jsfc = 1,ksfc_type
      IF ( jsfc == idx_lnd ) CYCLE ! computation below is valid only over water and ice

      CALL compute_qsat( kproma, kbdim, ppsfc, ptsfc(:,jsfc), pqsat_sfc(:,jsfc) )

      DO jl = 1,kproma

        pcpt_sfc(jl,jsfc) = ptsfc(jl,jsfc)*cpd*(1._wp+vtmpc2*pqsat_sfc(jl,jsfc))

        ztheta      = ptsfc(jl,jsfc)*(1.e5_wp/ppsfc(jl))**rd_o_cpd
        zthetav     = ztheta*(1._wp+vtmpc1*pqsat_sfc(jl,jsfc))

        zqtl       = pqm1_b(jl) + pqxm1_b(jl)  ! q_total at lowest model level
        zqts       = pqsat_sfc(jl,jsfc)        ! q_total at surface
        zqtmit     = 0.5_wp*( zqtl + zqts )    ! q_total, vertical average

        zqsmit     = 0.5_wp*( pqsat_b  (jl) + pqsat_sfc(jl,jsfc) )  ! qs
        ztmit      = 0.5_wp*( ptm1_b   (jl) + ptsfc    (jl,jsfc) )  ! temp.
        zthetamit  = 0.5_wp*( ptheta_b (jl) + ztheta  )  ! potential temp.
        zthetavmit = 0.5_wp*( pthetav_b(jl) + zthetav )  ! virtual potential temp.

        zfux = plh_b(jl)/(cpd*ztmit)
        zfox = plh_b(jl)/(rd*ztmit)

        zmult1 = 1._wp+vtmpc1*zqtmit   ! A in clear sky
        zmult2 = zfux*zmult1-zrvrd
        zmult3 = zrdrv*zfox*zqsmit/(1._wp+zrdrv*zfox*zfux*zqsmit)
        zmult5 = zmult1-zmult2*zmult3  ! A in cloud
        zmult4 = zfux*zmult5-1._wp     ! D in cloud

        zdus1 = paclc_b(jl)*zmult5 + (1._wp-paclc_b(jl))*zmult1   ! A avg
        zdus2 = paclc_b(jl)*zmult4 + (1._wp-paclc_b(jl))*vtmpc1   ! D avg

        zdqt     = zqtl - zqts                                    ! d qt
        zdthetal = pthetal_b(jl) - ztheta                         ! d theta_l
        zdu2(jl,jsfc) = MAX(zepdu2,(pum1_b(jl)-pocu(jl))**2 &     ! (d u)^2
                                  +(pvm1_b(jl)-pocv(jl))**2)      ! (d v)^2

        zbuoy        = zdus1*zdthetal + zdus2*zthetamit*zdqt
        pri_sfc(jl,jsfc) = pgeom1_b(jl)*zbuoy/(zthetavmit*zdu2(jl,jsfc))

        zalo = LOG( 1._wp + pgeom1_b(jl)/(g*pz0m(jl,jsfc)) )  ! ln[ 1 + zL/z0m ]
        zcdn(jl,jsfc) = (ckap/zalo)**2

        zcfnc(jl,jsfc)= SQRT(zdu2(jl,jsfc))*zcdn(jl,jsfc)

        zdthv         = MAX(0._wp,(zthetav-pthetav_b(jl)))
        zwst(jl,jsfc) = zdthv*SQRT(zdu2(jl,jsfc))/zthetavmit

        IF ( jsfc == idx_wtr ) THEN
          z0h        =pz0m(jl,jsfc)*EXP(2._wp-86.276_wp*pz0m(jl,jsfc)**0.375_wp)
          zaloh      =LOG(1._wp+pgeom1_b(jl)/(g*z0h))
          zchn  (jl,jsfc)=ckap**2/(zalo*zaloh)
          zcfnch(jl,jsfc)=SQRT(zdu2(jl,jsfc))*zchn(jl,jsfc)
          zcr   (jl)=(cfreec/(zchn(jl,jsfc)*SQRT(zdu2(jl,jsfc))))*ABS(zbuoy)**zonethird
        ELSE ! over ice
          zchn  (jl,jsfc) = zcdn (jl,jsfc)
          zcfnch(jl,jsfc) = zcfnc(jl,jsfc)  ! coeff. for scalar is the same as for momentum
        ENDIF
      ENDDO ! 1:kproma
    ENDDO   ! 1:ksfc_type

    IF ( idx_lnd.LE.ksfc_type ) THEN ! There is land surface in this simulation
      CALL finish('mo_turbulence_diag','computation over land surface not implemented')
    ENDIF

    !-------------------------------------------------------------------------
    ! Compute vertical shear of total water
    !-------------------------------------------------------------------------
    zcsat(1:kproma,1:ksfc_type) = 1._wp
   !IF ( idx_lnd.LE.ksfc_type ) zcsat(1:kproma,idx_lnd) = ...

    pqshear_sfc(1:kproma) = 0._wp ! initialization for weighted q_total at surface

    DO jsfc = 1,ksfc_type
      pqshear_sfc(1:kproma) =  pqshear_sfc(1:kproma)      &
                            & +  pqsat_sfc(1:kproma,jsfc) &
                            &       *zcsat(1:kproma,jsfc) &
                            &        *pfrc(1:kproma,jsfc)
    ENDDO

    pqshear_sfc(1:kproma) = ( pqm1_b(1:kproma)+pqxm1_b(1:kproma)          &
                          &  -pqshear_sfc(1:kproma) )*g/pgeom1_b(1:kproma)

    !-------------------------------------------------------------------------
    ! Compute the exchange coefficients for momentum, heat and water vapour,
    ! for each type of surface
    !-------------------------------------------------------------------------
    IF (lsfc_mom_flux.OR.lsfc_heat_flux) THEN  ! Surface flux is considered 

      ! stable case
  
      DO jsfc = 1,ksfc_type
        DO jl = 1,kproma
          IF ( pri_sfc(jl,jsfc) > 0._wp ) THEN
            zscf = SQRT(1._wp+pri_sfc(jl,jsfc))
            pcfm_sfc(jl,jsfc) = zcfnc (jl,jsfc)/(1._wp+z2b*pri_sfc(jl,jsfc)/zscf)   ! 5.2? 5.5
            pcfh_sfc(jl,jsfc) = zcfnch(jl,jsfc)/(1._wp+z2b*pri_sfc(jl,jsfc)*zscf)   ! 5.2? 5.6
            zch     (jl,jsfc) = zchn  (jl,jsfc)/(1._wp+z2b*pri_sfc(jl,jsfc)*zscf)   ! for (5.22)
          ENDIF
        ENDDO
      ENDDO
  
      ! unstable case
  
      IF (idx_wtr<=ksfc_type) THEN
        jsfc = idx_wtr  ! water
        DO jl = 1,kproma
          IF ( pri_sfc(jl,jsfc) <= 0._wp ) THEN
            zucf =  SQRT( -pri_sfc(jl,jsfc)*(1._wp+ pgeom1_b(jl)/(g*pz0m(jl,jsfc))) ) ! sqrt in (5.4)
            zucf =  1._wp+z3bc*zcdn(jl,jsfc)*zucf                          ! denominator in (5.4)
            zucf =  1._wp/zucf
            pcfm_sfc(jl,jsfc) = zcfnc (jl,jsfc)*(1._wp-z2b*pri_sfc(jl,jsfc)*zucf)  ! (5.2), (5.4)
            pcfh_sfc(jl,jsfc) = zcfnch(jl,jsfc)*(1._wp+zcr(jl)**cgam)**zrgam       ! (5.9)
            zch     (jl,jsfc) = zchn  (jl,jsfc)*(1._wp+zcr(jl)**cgam)**zrgam
          ENDIF
        ENDDO
      ENDIF
  
      IF (idx_ice<=ksfc_type) THEN
        jsfc = idx_ice  ! ice
        DO jl = 1,kproma
          IF ( pri_sfc(jl,jsfc) <= 0._wp ) THEN
            zucf =  SQRT( -pri_sfc(jl,jsfc)*(1._wp+ pgeom1_b(jl)/(g*pz0m(jl,jsfc))) ) ! sqrt in (5.4)
            zucf =  1._wp+z3bc*zcdn(jl,jsfc)*zucf                   ! denominator in (5.4)
            zucf =  1._wp/zucf
            pcfm_sfc(jl,jsfc) = zcfnc (jl,jsfc)*(1._wp-z2b*pri_sfc(jl,jsfc)*zucf)  ! (5.2), (5.4)
            pcfh_sfc(jl,jsfc) = zcfnch(jl,jsfc)*(1._wp-z3b*pri_sfc(jl,jsfc)*zucf)  ! (5.2), (5.4)
            zch     (jl,jsfc) = zchn  (jl,jsfc)*(1._wp-z3b*pri_sfc(jl,jsfc)*zucf)
          ENDIF
        ENDDO
      ENDIF
     
      jsfc = idx_lnd ! land, not implemented
      ! z0 for heat over land is different, thus zucf is also different

    END IF  ! lsfc_mom_flux.OR.lsfc_heat_flux

    IF (.NOT.lsfc_mom_flux) THEN  ! Surface momentum flux is switched off 
      pcfm_sfc(1:kproma,1:ksfc_type) = 0._wp
    END IF

    IF (.NOT.lsfc_heat_flux) THEN  ! Surface heat flux is switched off 
      pcfh_sfc(1:kproma,1:ksfc_type) = 0._wp
      zwst    (1:kproma,1:ksfc_type) = 0._wp   ! affects TKE at surface
    END IF

    !-------------------------------------------------------------------------
    ! Get the aggregated exchange coefficient for momentum
    !-------------------------------------------------------------------------
    pcfm_gbm(1:kproma) = 0._wp
    DO jsfc = 1,ksfc_type
      pcfm_gbm(1:kproma) = pcfm_gbm(1:kproma) + pfrc(1:kproma,jsfc)*pcfm_sfc(1:kproma,jsfc)
    ENDDO

    !-------------------------------------------------------------------------
    ! For sensible heat and water vapour, it is not the exchange coefficient
    ! but the solution at the lowest model level that is aggregated.
    ! The exchange coeffcients (cfh_sfc) computed in this subroutine
    ! are returned to the calling subroutine for each surface type.
    ! They are used later for solving the discretized veritical diffusion
    ! equation at the lowest grid level (klev) for each surface type
    ! separately. Then the solutions are aggregated using
    ! (fraction of type)*(cfh_sfc of type) as the weighting factor, which 
    ! ensures conservation of the area-weighted total flux.
    !   Here we compute the aggregated the exchange coefficient for output
    ! and for solving the vertical diffusion equation of the variance of 
    ! virtual potential temperature (theta_v).
    !-------------------------------------------------------------------------
    pcfh_gbm(1:kproma) = 0._wp                                                         
    DO jsfc = 1,ksfc_type                                                                 
      pcfh_gbm(1:kproma) = pcfh_gbm(1:kproma) + pfrc(1:kproma,jsfc)*pcfh_sfc(1:kproma,jsfc)
    ENDDO     

    !-------------------------------------------------------------------------
    ! Hydrometeors and the other tracers share the same exchange coefficient 
    ! and heat and moisture, but have no turbulence-induced surface flux.
    ! These are taken care of in subroutine matrix_setup.
    !   The total water variance has a different exchange coefficient 
    ! (variable cfv), and no surface flux. Set the surfce exchange coefficient
    ! to zero. 
    !-------------------------------------------------------------------------
    pcfv_sfc(1:kproma) = 0._wp

    !-------------------------------------------------------------------------
    ! Diagnose friction velocity. The values of each individual surface type
    ! are used immediately below for computing the surface value of TKE;
    ! The grid-box mean is used in the next time step in subroutine
    ! "atm_exchange_coeff" for finding the PBL height.
    !-------------------------------------------------------------------------
    IF (lsfc_mom_flux) THEN  ! Surface momentum flux is switched on

      z2m = 2._wp            ! 2-m height
      DO jsfc = 1,ksfc_type
        DO jl = 1,kproma
          lhighz0 = pz0m(jl,jsfc).GT.z2m
          zcdn2m = MERGE((ckap/LOG(1._wp+pgeom1_b(jl)/(g*z2m)))**2, &
                 &        zcdn(jl,jsfc),lhighz0 )
          zcdnr  = zcdn2m/zcdn(jl,jsfc)
  
          zucf   = SQRT( ABS(pri_sfc(jl,jsfc))*(1._wp+pgeom1_b(jl)/(g*z2m)) )  ! sqrt in (5.4)
          zucf  = 1._wp + z3bc*zcdn2m*zucf                              ! denomenator in (5.4)
          zucf = 1._wp - z2b*pri_sfc(jl,jsfc)/zucf                      ! (5.4)
          zcfm2m = MERGE( zcfnc(jl,jsfc)*zcdnr*zucf, pcfm_sfc(jl,jsfc)*zcdnr, &
                 &        lhighz0.AND.pri_sfc(jl,jsfc).LT.0._wp           )
          zust   = zcfm2m*SQRT(zdu2(jl,jsfc))
          zustar(jl,jsfc) = SQRT(zust)
        END DO
      END DO

      pustarm(1:kproma) = 0._wp
      DO jsfc = 1,ksfc_type
        pustarm(1:kproma) = pustarm(1:kproma) + pfrc(1:kproma,jsfc)*zustar(1:kproma,jsfc)
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

    DO jl = 1,kproma
      DO jsfc = 1,ksfc_type
        ztkev = custf*(zustar(jl,jsfc)**2)
        IF(zwst(jl,jsfc).GT.zepsr) THEN
           zconvs = (zwst(jl,jsfc)*zch(jl,jsfc)*pghabl(jl))**zonethird
           zstabf = (pgeom1_b(jl)*ckap*zwst(jl,jsfc)*zch(jl,jsfc))**ztwothirds
           zstabf = MIN(custf*3._wp*zustar(jl,jsfc)**2,zstabf)
           ztkev = ztkev + zstabf + cwstf*(zconvs**2)
        END IF
        ptkevn_sfc(jl) = ptkevn_sfc(jl) + ztkev*pfrc(jl,jsfc)
      END DO
      ptkevn_sfc(jl) = MAX( tkemin,ptkevn_sfc(jl) )
    END DO

#ifdef __ICON__
    ! should the same be done as in echam?
#else
     IF (lstart) THEN
        DO 345 jl=1,kproma
           ptkem1_sfc(jl) = ptkevn_sfc(jl)
           ptkem0_sfc(jl) = ptkevn_sfc(jl)
345     END DO
     END IF
#endif

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
      IF ( jsfc == idx_lnd ) CALL finish('sfc_exchange_coeff','land surface not implemented')
      ztvsfc(1:kproma) = ptsfc(1:kproma,jsfc)*(1._wp + vtmpc1*pqsat_sfc(1:kproma,jsfc))
      prho_sfc(1:kproma) = prho_sfc(1:kproma) + zrrd*ppsfc(1:kproma)/ztvsfc(1:kproma)
    END DO

  END SUBROUTINE sfc_exchange_coeff
  !-------------

END MODULE mo_turbulence_diag
