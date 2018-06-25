!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
#endif
MODULE mo_ssodrag

  USE mo_kind,                  ONLY: wp
  USE mo_math_constants,        ONLY: pi
  USE mo_physical_constants,    ONLY: grav, rd, cpd
  USE mo_echam_sso_config,      ONLY: echam_sso_config,               &
       &                              gfrcrit, grcrit, grahilo,       &
       &                              gsigcr , gssec , gtsec , gvsec

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ssodrag

CONTAINS

SUBROUTINE ssodrag ( jg            ,& ! in,  grid index
  &                  kproma        ,& ! in,  loop length in block of cells/columns
  &                  kbdim         ,& ! in,  dimension of block of cells/columns
  &                  klev          ,& ! in,  number of levels
  !
  &                  pdtime        ,& ! in,  length of timestep (s)
  &                  pcoriol       ,& ! in,  Coriolis parameter (1/s)
  &                  pzf           ,& ! in,  full level height (m)
  &                  pzs           ,& ! in,  surface height    (m)
  !
  &                  paphm1        ,& ! in,  p at half levels
  &                  papm1         ,& ! in,  p at full levels
  &                  pmair         ,& ! in,  air mass (kg/m2)
  &                  ptm1          ,& ! in,  T
  &                  pum1          ,& ! in,  u
  &                  pvm1          ,& ! in,  v
  !
  &                  pmea          ,& ! in,  Mean Orography (m)
  &                  pstd          ,& ! in,  SSO standard deviation (m)
  &                  psig          ,& ! in,  SSO slope
  &                  pgam          ,& ! in,  SSO Anisotropy
  &                  pthe          ,& ! in,  SSO Angle
  &                  ppic          ,& ! in,  SSO Peaks elevation (m)
  &                  pval          ,& ! in,  SSO Valleys elevation (m)
  !
  &                  psftlf        ,& ! in,  area fraction of land incl. lakes
  !                                          where the SSO params are valid
  !
  &                  pustrgw       ,& ! out, u-gravity wave stress
  &                  pvstrgw       ,& ! out, v-gravity wave stress
  &                  pvdisgw       ,& ! out, dissipation by gravity wave drag
  !
  &                  pdis_sso      ,& ! out, sso energy dissipation per mass
  &                  pdu_sso       ,& ! out, sso tendency of zonal wind
  &                  pdv_sso        ) ! out, sso tendency of meridional wind

  !
  ! Description:
  !
  ! Does the SSO drag following LOTT &MILLER (1997)
  !
  !   *ssodrag* is called from *physc*
  !
  ! Method:
  ! Authors:
  !
  ! m.miller + b.ritter   e.c.m.w.f.     15/06/1986.
  ! f.lott + m. miller    e.c.m.w.f.     22/11/1994
  ! e.manzini             mpi            05.09.2000
  ! m.giorgetta           mpi            31.01.2017
  !
  ! for more details see file AUTHORS
  !

  ! scalar arguments with intent(IN):
  INTEGER,  INTENT(in)    :: jg
  INTEGER,  INTENT(in)    :: kproma, kbdim, klev
  REAL(wp), INTENT(in)    :: pdtime               ! length oftimestep (s)

  ! array arguments with intent(IN):
  ! 1D
  REAL(wp), INTENT(in)    :: pcoriol(kbdim)       ! Coriolis parameter (1/s)
  !
  REAL(wp), INTENT(in)    :: pmea(kbdim)          ! Mean Orography (m)
  REAL(wp), INTENT(in)    :: pstd(kbdim)          ! SSO standard deviation (m)
  REAL(wp), INTENT(in)    :: psig(kbdim)          ! SSO slope
  REAL(wp), INTENT(in)    :: pgam(kbdim)          ! SSO Anisotropy
  REAL(wp), INTENT(in)    :: pthe(kbdim)          ! SSO Angle
  REAL(wp), INTENT(in)    :: ppic(kbdim)          ! SSO Peacks elevation (m)
  REAL(wp), INTENT(in)    :: pval(kbdim)          ! SSO Valleys elevation (m)
  !
  REAL(wp), INTENT(in)    :: psftlf(kbdim)        ! area fraction of land incl. lakes
  !
  ! 2D
  REAL(wp), INTENT(in)    :: paphm1(kbdim,klev+1) ! half level pressure (t-dt)
  REAL(wp), INTENT(in)    :: papm1(kbdim,klev)    ! full level pressure (t-dt)
  REAL(wp), INTENT(in)    :: pmair(kbdim,klev)    ! air mass (kg/m2)
  REAL(wp), INTENT(in)    :: pzf  (kbdim,klev)    ! full level height (m)
  REAL(wp), INTENT(in)    :: pzs  (kbdim)         ! surface height    (m)
  REAL(wp), INTENT(in)    :: ptm1(kbdim,klev)     ! temperature (t-dt)
  REAL(wp), INTENT(in)    :: pum1(kbdim,klev)     ! zonal wind (t-dt)
  REAL(wp), INTENT(in)    :: pvm1(kbdim,klev)     ! meridional wind (t-dt)

  ! array arguments with intent(OUT):
  ! 1D
  REAL(wp), INTENT(out)   :: pustrgw(kbdim)       ! u-gravity wave stress
  REAL(wp), INTENT(out)   :: pvstrgw(kbdim)       ! v-gravity wave stress
  REAL(wp), INTENT(out)   :: pvdisgw(kbdim)       ! dissipation by gravity wave drag

  ! array arguments with intent(OUT):
  ! 2D
  REAL(wp), INTENT(out)   :: pdis_sso(kbdim,klev) ! sso energy dissipation rate
  REAL(wp), INTENT(out)   :: pdu_sso(kbdim,klev)  ! sso tendency of zonal wind
  REAL(wp), INTENT(out)   :: pdv_sso(kbdim,klev)  ! sso tendency of wind

  ! Local scalars:
  INTEGER  :: igwd, jk, jl, ji

  ! Local arrays:
  INTEGER  :: idx(kbdim), itest(kbdim)
  REAL(wp) :: zhgeo(kbdim,klev)     ! geopot. height above ground (m)
  REAL(wp) :: zdu_oro(kbdim,klev)   ! tendency due to ORO GW DRAG  (m/s)
  REAL(wp) :: zdv_oro(kbdim,klev)   ! tendency due to ORO GW DRAG  (m/s)
  REAL(wp) :: zdis_oro(kbdim,klev)  ! energy dissipation due to ORO GW DRAG  (?)
  REAL(wp) :: zdu_lif(kbdim,klev)   ! tendency due to MOUNTAIN LIFT(m/s)
  REAL(wp) :: zdv_lif(kbdim,klev)   ! tendency due to MOUNTAIN LIFT(m/s)
  REAL(wp) :: zdis_lif(kbdim,klev)  ! energy dissipation due to MOUNTAIN LIFT(?)

  ! Shortcuts to components of echam_sso_config
  !
  REAL(wp), POINTER :: gpicmea, gstd
  REAL(wp), POINTER :: gkdrag, gkwake, gklift
  !
  gpicmea => echam_sso_config(jg)% gpicmea
  gstd    => echam_sso_config(jg)% gstd
  gkwake  => echam_sso_config(jg)% gkwake
  gkdrag  => echam_sso_config(jg)% gkdrag
  gklift  => echam_sso_config(jg)% gklift

  !
  !*         1.    initialization
  !                --------------

  pustrgw (:)   = 0.0_wp
  pvstrgw (:)   = 0.0_wp
  pvdisgw (:)   = 0.0_wp

  pdis_sso(:,:) = 0.0_wp
  pdu_sso (:,:) = 0.0_wp
  pdv_sso (:,:) = 0.0_wp

  zhgeo   (:,:) = pzf(:,:)-SPREAD(pzs(:),2,klev)
  
  zdu_oro (:,:) = 0.0_wp
  zdv_oro (:,:) = 0.0_wp
  zdis_oro(:,:) = 0.0_wp

  zdu_lif (:,:) = 0.0_wp
  zdv_lif (:,:) = 0.0_wp
  zdis_lif(:,:) = 0.0_wp


  !  SELECTION  POINTS WHERE THE SCHEME IS ACTIVE

  igwd=0
  idx(:) = 0
  DO jl=1,kproma
     itest(jl)=0
     IF (((ppic(jl)-pmea(jl)) > gpicmea).AND.(pstd(jl) > gstd)) THEN
        itest(jl)=1
        igwd=igwd+1
        idx(igwd)=jl
     ENDIF
  ENDDO


  IF (.NOT.((gkwake == 0.0_wp).AND.(gkdrag == 0.0_wp))) THEN
     !
     !*         2.    orographic gravity wave drag
     !                -----------------------------
     !
     CALL orodrag( jg, kproma,  kbdim,   klev,                       &
          &        pdtime,                                           &
          &        igwd,    idx,                                     &
          &        zhgeo,   paphm1,  papm1,                          &
          &        pmair,                                            &
          &        ptm1,    pum1,    pvm1,                           &
          &        pmea,    pstd,    psig,  pgam, pthe, ppic, pval,  &
          &        zdu_oro, zdv_oro, zdis_oro)
     !
  END IF


  IF (.NOT.(gklift == 0.0_wp)) THEN
     !
     !*         3.    mountain lift
     !                --------------
     CALL orolift( jg, kproma,  kbdim,   klev,                       &
          &        pcoriol,                                          &
          &        pdtime,                                           &
          &        itest,                                            &
          &        zhgeo,   paphm1,                                  &
          &        pmair,                                            &
          &        ptm1,    pum1,    pvm1,                           &
          &        pmea,    pstd,    ppic,                           &
          &        zdu_lif, zdv_lif, zdis_lif )
     !
  END IF

  ! STRESS FROM TENDENCIES

  ! Scale with the area fraction of land incl. lakes psftlf, because:
  ! - the SSO parameters are derived for this part of the surface only
  ! - the effects happen only above this part of the surface, assuming
  !   perfect vertical propagation of gravity waves as used here.

  DO jk = 1, klev
!CDIR NODEP
     DO jl = 1, igwd
        ji=idx(jl)
        pustrgw(ji) = pustrgw(ji)+( zdu_oro(ji,jk)+ zdu_lif(ji,jk))*pmair(ji,jk)*psftlf(ji)
        pvstrgw(ji) = pvstrgw(ji)+( zdv_oro(ji,jk)+ zdv_lif(ji,jk))*pmair(ji,jk)*psftlf(ji)
        pvdisgw(ji) = pvdisgw(ji)+(zdis_oro(ji,jk)+zdis_lif(ji,jk))*pmair(ji,jk)*psftlf(ji)
     ENDDO
  ENDDO
  !
  !*         4.    total quantities
  !                ----------------

  do jk=1,klev
!CDIR NODEP
    do jl=1,igwd
      ji=idx(jl)
      pdis_sso(ji,jk)= (zdis_oro(ji,jk) +zdis_lif(ji,jk))*psftlf(ji)
      pdu_sso(ji,jk) = ( zdu_oro(ji,jk) + zdu_lif(ji,jk))*psftlf(ji)
      pdv_sso(ji,jk) = ( zdv_oro(ji,jk) + zdv_lif(ji,jk))*psftlf(ji)
    enddo
  enddo

END SUBROUTINE ssodrag

SUBROUTINE orodrag( jg, kproma, kbdim,  klev,                        &
                    pdtime,                                          &
                    kgwd,   kdx,                                     &
                    phgeo,  paphm1, papm1,                           &
                    pmair,                                           &
                    ptm1,   pum1,   pvm1,                            &
                    pmea,   pstd,   psig,   pgam, pthe, ppic, pval,  &
                    pvom,   pvol,   pdis)
  !
  !
  !**** *orodrag* - does the SSO drag  parametrization.
  !
  !     purpose.
  !     --------
  !
  !     this routine computes the physical tendencies of the
  !     prognostic variables u,v  and t due to  vertical transports by
  !     subgridscale orographically excited gravity waves, and to
  !     low level blocked flow drag.
  !
  !          called from *ssodrag*.
  !
  !     author.
  !     -------
  !     m.miller + b.ritter   e.c.m.w.f.     15/06/86.
  !
  !     f.lott + m. miller    e.c.m.w.f.     22/11/94

  IMPLICIT NONE

  ! scalar arguments with intent(IN):
  INTEGER,  INTENT(in)  :: jg, kproma, kbdim, klev
  INTEGER,  INTENT(in)  :: kgwd      ! Total points where oro scheme is active
  REAL(wp), INTENT(in)  :: pdtime


  ! array arguments with intent(IN):
  ! Input 1D
  INTEGER,  INTENT(in)  :: kdx(kbdim)    ! physical location of an active point

  REAL(wp), INTENT(in)  :: pmea(kbdim)   ! Mean Orography (m)
  REAL(wp), INTENT(in)  :: pstd(kbdim)   ! SSO standard deviation (m)
  REAL(wp), INTENT(in)  :: psig(kbdim)   ! SSO slope
  REAL(wp), INTENT(in)  :: pgam(kbdim)   ! SSO Anisotropy
  REAL(wp), INTENT(in)  :: pthe(kbdim)   ! SSO Angle
  REAL(wp), INTENT(in)  :: ppic(kbdim)   ! SSO Peacks elevation (m)
  REAL(wp), INTENT(in)  :: pval(kbdim)   ! SSO Valleys elevation (m)
  ! Input 2D
  REAL(wp), INTENT(in)  :: paphm1(kbdim,klev+1) ! half level pressure (t-dt)
  REAL(wp), INTENT(in)  :: papm1(kbdim,klev)    ! full level pressure (t-dt)
  REAL(wp), INTENT(in)  :: phgeo(kbdim,klev)    ! geopotential height above surface
  REAL(wp), INTENT(in)  :: pmair(kbdim,klev)    ! air mass (kg/m2)
  REAL(wp), INTENT(in)  :: ptm1(kbdim,klev)     ! temperature (t-dt)
  REAL(wp), INTENT(in)  :: pum1(kbdim,klev)     ! zonal wind (t-dt)
  REAL(wp), INTENT(in)  :: pvm1(kbdim,klev)     ! meridional wind (t-dt)

  ! Output 2D
  REAL(wp), INTENT(out) :: pdis(kbdim,klev)     ! energy dissipation rate per mass
  REAL(wp), INTENT(out) :: pvol(kbdim,klev)     ! tendency of meridional wind
  REAL(wp), INTENT(out) :: pvom(kbdim,klev)     ! tendency of zonal wind

  !  Local scalars:
  INTEGER :: jl, jk, ji
  REAL(wp)    :: ztemp, zforc, ztend, rover, &
       zb, zc, zconb, zabsv, zzd1, ratio, zbet,        &
       zust, zvst, zdis, zred

  !  Local arrays:
  INTEGER  ::    icrit(kbdim),     ikcrith(kbdim), &
       ikenvh(kbdim),   iknu(kbdim),      iknu2(kbdim),          &
       ikcrit(kbdim)

  REAL(wp) ::  zdudt(kbdim),    zdvdt(kbdim),                       &
       &       znu(kbdim),      zd1(kbdim),       zd2(kbdim),       &
       &       zdmod(kbdim),    ztau0(kbdim)

  REAL(wp) ::  ztau(kbdim,klev+1),  zstab(kbdim,klev+1),  zvph(kbdim,klev+1),  &
       &       zrho(kbdim,klev+1),  zri(kbdim,klev+1),    zpsi(kbdim,klev+1),  &
       &       zzdep(kbdim,klev),   pulow(kbdim),         pvlow(kbdim)

  ! Shortcuts to components of echam_sso_config
  !
  REAL(wp), POINTER :: gkdrag, gkwake
  !
  gkdrag => echam_sso_config(jg)% gkdrag
  gkwake => echam_sso_config(jg)% gkwake


  !  Executable statements

  !*         1.1   computational constants
  !*                ------------ ----------
  !
  !     ------------------------------------------------------------------
  !
  !*         2.     precompute basic state variables.
  !*                ---------- ----- ----- ----------
  !*                define low level wind, project winds in plane of
  !*                low level wind, determine sector in which to take
  !*                the variance and set indicator for critical levels.
  !
  CALL orosetup( jg, kproma, kbdim,  klev,   kgwd,   kdx,               &
       &         ikcrit, ikcrith,icrit,  ikenvh, iknu,   iknu2,         &
       &         paphm1, papm1,  pmair,  pum1,   pvm1,   ptm1,   phgeo, &
       &         zrho  , zri,    zstab,  ztau,   zvph,   zpsi,   zzdep, &
       &         pulow,  pvlow,                                         &
       &         pthe,   pgam,   pmea,   ppic,   pval,                  &
       &         znu, zd1, zd2,  zdmod )
  !
  !***********************************************************


  IF (.NOT.(gkdrag == 0.0_wp)) THEN
     !
     !*         3.      compute low level stresses using subcritical and
     !*                 supercritical forms.computes anisotropy coefficient
     !*                 as measure of orographic twodimensionality.
     !
     CALL gwstress( jg, kbdim,  klev,                                      &
          &         kgwd,   kdx,                                           &
          &         ikenvh,                                                &
          &         pstd,   psig,   ppic,   pval,                          &
          &         zrho,   zstab,  zvph,   zdmod,                         &
          &         ztau0 )
     !
     !*         4.      compute stress profile including
     !                  trapped waves, wave breaking,
     !                  linear decay in stratosphere.
     !
     CALL gwprofil( jg, kbdim,  klev,                                      &
          &         kgwd,   kdx,                                           &
          &         ikcrith,icrit,                                         &
          &         pstd,   psig,                                          &
          &         paphm1, zrho,   zri,    zstab,  zvph,   zdmod,         &
          &         ztau0,  ztau )
     !
  END IF
  
  !
  !*         5.      Compute tendencies from waves stress profile.
  !                  Compute low level blocked flow drag.
  !*                 --------------------------------------------
  !
  !  explicit solution at all levels for the gravity wave
  !  implicit solution for the blocked levels
  !
  zdudt(:)   = 0.0_wp
  zdvdt(:)   = 0.0_wp
  !
  pvom (:,:) = 0.0_wp
  pvol (:,:) = 0.0_wp
  pdis (:,:) = 0.0_wp  
  !
  DO 524 jk=1,klev
!CDIR NODEP
!DIR$ CONCURRENT
     DO 523 jl=1,kgwd
        ji=kdx(jl)
        !
        !  WAVE STRESS
        !-------------
        !
        IF (.NOT.(gkdrag == 0._wp)) THEN
           !
           ztemp     = -(ztau(ji,jk+1)-ztau(ji,jk))/(zvph(ji,klev+1)*pmair(ji,jk))
           zdudt(ji) = (pulow(ji)*zd1(ji)-pvlow(ji)*zd2(ji))*ztemp/zdmod(ji)
           zdvdt(ji) = (pvlow(ji)*zd1(ji)+pulow(ji)*zd2(ji))*ztemp/zdmod(ji)
           !
           ! Control Overshoots
           !
           zforc=SQRT(zdudt(ji)**2+zdvdt(ji)**2)
           ztend=SQRT(pum1(ji,jk)**2+pvm1(ji,jk)**2)/pdtime
           rover=0.25_wp
           IF(zforc >= rover*ztend)THEN
              zdudt(ji)=rover*ztend/zforc*zdudt(ji)
              zdvdt(ji)=rover*ztend/zforc*zdvdt(ji)
           ENDIF
           !
        ELSE
           !
           zdudt(ji) = 0.0_wp
           zdvdt(ji) = 0.0_wp
           !
        END IF
        !
        !
        ! BLOCKED FLOW DRAG:
        ! -----------------
        !
        IF(jk > ikenvh(ji)) THEN
           !
           IF (.NOT.(gkwake == 0._wp)) THEN
              !
              zb    = 1.0_wp-0.18_wp*pgam(ji)-0.04_wp*pgam(ji)**2
              zc    =        0.48_wp*pgam(ji)+0.30_wp*pgam(ji)**2
              zconb = 2._wp*pdtime*gkwake*psig(ji)/(4._wp*pstd(ji))
              zabsv = SQRT(pum1(ji,jk)**2+pvm1(ji,jk)**2)/2._wp
              zzd1  = zb*COS(zpsi(ji,jk))**2+zc*SIN(zpsi(ji,jk))**2
              ratio = (COS(zpsi(ji,jk))**2+pgam(ji)*SIN(zpsi(ji,jk))**2)/      &
                   &  (pgam(ji)*COS(zpsi(ji,jk))**2+SIN(zpsi(ji,jk))**2)
              zbet=MAX(0.0_wp,2.0_wp-1.0_wp/ratio)*zconb*zzdep(ji,jk)*zzd1*zabsv
              !
              ! OPPOSED TO THE WIND
              !
              zdudt(ji)=-pum1(ji,jk)/pdtime
              zdvdt(ji)=-pvm1(ji,jk)/pdtime
              !
              ! PERPENDICULAR TO THE SSO MAIN AXIS:
              !
              !mod     zdudt(ji)=-(pum1(ji,jk)*cos(pthe(ji)*pi/180.)
              !mod *              +pvm1(ji,jk)*sin(pthe(ji)*pi/180.))
              !mod *              *cos(pthe(ji)*pi/180.)/pdtime
              !mod     zdvdt(ji)=-(pum1(ji,jk)*cos(pthe(ji)*pi/180.)
              !mod *              +pvm1(ji,jk)*sin(pthe(ji)*pi/180.))
              !mod *              *sin(pthe(ji)*pi/180.)/pdtime
              !
              zdudt(ji) = zdudt(ji)*(zbet/(1._wp+zbet))
              zdvdt(ji) = zdvdt(ji)*(zbet/(1._wp+zbet))
              !
           ELSE
              !
              zdudt(ji) = 0.0_wp
              zdvdt(ji) = 0.0_wp
              !
           END IF
           !
        END IF
        !
        pvom(ji,jk) = zdudt(ji)
        pvol(ji,jk) = zdvdt(ji)
        zust = pum1(ji,jk)+pdtime*zdudt(ji)
        zvst = pvm1(ji,jk)+pdtime*zdvdt(ji)
        zdis = 0.5_wp*(pum1(ji,jk)**2+pvm1(ji,jk)**2-zust**2-zvst**2)
        IF(zdis<0.0_wp) THEN
           zred=SQRT((pum1(ji,jk)**2+pvm1(ji,jk)**2)/(zust**2+zvst**2))
           zust=zust*zred
           zvst=zvst*zred
           pvom(ji,jk)=(zust-pum1(ji,jk))/pdtime
           pvol(ji,jk)=(zvst-pvm1(ji,jk))/pdtime
           zdis=0.5_wp*(pum1(ji,jk)**2+pvm1(ji,jk)**2-zust**2-zvst**2)
        END IF
        pdis(ji,jk)=zdis/pdtime
523  END DO
524 END DO

END SUBROUTINE orodrag

SUBROUTINE orosetup                                           &
     ( jg, kproma, kbdim, klev, kgwd, kdx                     &
     , kkcrit, kkcrith,kcrit                                  &
     , kkenvh, kknu  , kknu2                                  &
     , paphm1, papm1 , pmair , pum1  , pvm1  , ptm1  , phgeo  &
     , prho  , pri   , pstab , ptau  , pvph  , ppsi  , pzdep  &
     , pulow , pvlow                                          &
     , ptheta, pgam, pmea, ppic, pval                         &
     , pnu   , pd1  ,  pd2  ,pdmod  )
  !
  !**** *gwsetup*
  !
  !     purpose.
  !     --------
  !     SET-UP THE ESSENTIAL PARAMETERS OF THE SSO DRAG SCHEME:
  !     DEPTH OF LOW WBLOCKED LAYER, LOW-LEVEL FLOW, BACKGROUND
  !     STRATIFICATION.....
  !
  !    called from *orodrag*
  !
  !     ==== inputs ===
  !
  !
  ! ptsphy--input-R-Time-step (s)
  ! paphm1--input-R: pressure at model 1/2 layer
  ! papm1---input-R: pressure at model layer
  ! pmair---input-R: air mass (kg/m2)
  ! phgeo---input-R: Altitude of layer above ground
  ! ptm1, pum1, pvm1--R-: t, u and v
  ! pmea----input-R-Mean Orography (m)
  ! psig----input-R-SSO slope
  ! pgam----input-R-SSO Anisotropy
  ! pthe----input-R-SSO Angle
  ! ppic----input-R-SSO Peacks elevation (m)
  ! pval----input-R-SSO Valleys elevation (m)
  !
  !     ==== outputs ===
  ! pulow, pvlow -output-R: Low-level wind
  ! kkcrit----I-: Security value for top of low level flow
  ! kcrit-----I-: Critical level
  ! kkenvh----I-: Top of blocked flow layer
  ! kknu------I-: Layer that sees mountain peacks
  ! kknu2-----I-: Layer that sees mountain peacks above mountain mean
  ! kknub-----I-: Layer that sees mountain mean above valleys
  ! prho------R-: Density at 1/2 layers
  ! pri-------R-: Background Richardson Number, Wind shear measured along
  !               GW stress
  ! pstab-----R-: Brunt-Vaisala freq. at 1/2 layers
  ! pvph------R-: Wind in  plan of GW stress, Half levels.
  ! ppsi------R-: Angle between low level wind and SS0 main axis.
  ! pd1-------R-| Compared the ratio of the stress
  ! pd2-------R-| that is along the wind to that Normal to it.
  !               pdi define the plane of low level stress
  !               compared to the low level wind.
  !
  ! see p. 108 Lott & Miller (1997).
  !
  ! pdmod-----R-: Norme of pdi
  !
  !     === local arrays ===
  !
  ! zvpf------R-: Wind projected in the plan of the low-level stress.
  !
  !        implicit arguments :   none
  !        --------------------
  !
  !     reference.
  !     ----------
  !
  !        see ecmwf research department documentation of the "i.f.s."
  !
  !     author.
  !     -------
  !
  !     modifications.
  !     --------------
  !     f.lott  for the new-gwdrag scheme november 1993
  !
  !-----------------------------------------------------------------------
  !

  IMPLICIT NONE

  ! scalar arguments with intent(IN):
  INTEGER, INTENT(in) :: jg, kproma, kbdim, klev, kgwd

  INTEGER :: kkcrit(kbdim), kkcrith(kbdim), kcrit(kbdim),    &
       kdx(kbdim), kkenvh(kbdim)

  REAL(wp) :: pulow(kbdim), pvlow(kbdim), ptheta(kbdim), pgam(kbdim), pnu(kbdim),  &
       pd1(kbdim),   pd2(kbdim),   pdmod(kbdim),                          &
       pmea(kbdim),  ppic(kbdim),  pval(kbdim)

  ! array arguments with intent(IN):
  ! Input 2D
  REAL(wp), INTENT(in)  :: paphm1(kbdim,klev+1) ! half level pressure (t-dt)
  REAL(wp), INTENT(in)  :: papm1(kbdim,klev)    ! full level pressure (t-dt)
  REAL(wp), INTENT(in)  :: pmair(kbdim,klev)    ! air mass (kg/m2)
  REAL(wp), INTENT(in)  :: phgeo(kbdim,klev)    ! geopotential height above surface
  REAL(wp), INTENT(in)  :: ptm1(kbdim,klev)     ! temperature (t-dt)
  REAL(wp), INTENT(in)  :: pum1(kbdim,klev)     ! zonal wind (t-dt)
  REAL(wp), INTENT(in)  :: pvm1(kbdim,klev)     ! meridional wind (t-dt)

  REAL(wp), INTENT(out) :: prho(kbdim,klev+1), pri(kbdim,klev+1),  pstab(kbdim,klev+1),    &
       &                   ptau(kbdim,klev+1), pvph(kbdim,klev+1), ppsi(kbdim,klev+1),     &
       &                   pzdep(kbdim,klev)
  !
  !  Local scalars:
  INTEGER  :: ilevh , jl, jk, ji
  REAL(wp) :: zcons1, zcons2, zu,   zphi,    &
       &      zvt1,   zvt2,   zdwind, zwind, &
       &      zstabm, zstabp, zrhom,  zrhop
  LOGICAL  :: lo
  !
  ! local arrays
  INTEGER  :: kknu(kbdim),  kknu2(kbdim), kknub(kbdim), kknul(kbdim)
  REAL(wp) :: znorm(kbdim), zb(kbdim),   zc(kbdim), znup(kbdim), znum(kbdim)
  REAL(wp) :: zhcrit(kbdim,klev), zvpf(kbdim,klev),   zdp(kbdim,klev)
  REAL(wp) :: zmair(kbdim), zmair_inv
  LOGICAL  :: ll1(kbdim,klev+1)

  ! Shortcuts to components of echam_sso_config
  !
  INTEGER , POINTER :: nktopg
  !
  nktopg => echam_sso_config(jg)% nktopg


  !     ------------------------------------------------------------------
  !*         1.1   computational constants
  !                -----------------------
  !
  ilevh =klev/3
  !
  zcons1=1._wp/rd
  zcons2=grav**2/cpd
  !
  !     ------------------------------------------------------------------
  !*         2.1     define low level wind, project winds in plane of
  !*                 low level wind, determine sector in which to take
  !*                 the variance and set indicator for critical levels.
  !
  DO jl=1,kproma
     kknu(jl)    =klev
     kknu2(jl)   =klev
     kknub(jl)   =klev
     kknul(jl)   =klev
     pgam(jl) =MAX(pgam(jl),gtsec)
     ll1(jl,klev+1)=.FALSE.
  END DO
  !
  ! Ajouter une initialisation (L. Li, le 23fev99):
  !
  DO jk=klev,ilevh,-1
     DO jl=1,kproma
        ll1(jl,jk)= .TRUE.
     END DO
  END DO
  !
  !*      define top of low level flow
  !       ----------------------------
  DO 2002 jk=klev,ilevh,-1
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
    DO 2003 ji=1,kgwd
       jl = kdx(ji)
       lo=(paphm1(jl,jk)/paphm1(jl,klev+1)) >= gsigcr
       IF(lo) THEN
          kkcrit(jl)=jk
       ENDIF
       zhcrit(jl,jk)=ppic(jl)-pval(jl)
       ll1(jl,jk)=(phgeo(jl,jk) > zhcrit(jl,jk))
       IF(ll1(jl,jk).NEQV.ll1(jl,jk+1)) THEN
          kknu(jl)=jk
       ENDIF
       IF(.NOT.ll1(jl,ilevh))kknu(jl)=ilevh
2003 END DO
2002 END DO

  DO 2004 jk=klev,ilevh,-1
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
     DO 2005 ji=1,kgwd
        jl = kdx(ji)
        zhcrit(jl,jk)=ppic(jl)-pmea(jl)
        ll1(jl,jk)=(phgeo(jl,jk) > zhcrit(jl,jk))
        IF(ll1(jl,jk).NEQV.ll1(jl,jk+1)) THEN
           kknu2(jl)=jk
        ENDIF
        IF(.NOT.ll1(jl,ilevh))kknu2(jl)=ilevh
2005 END DO
2004 END DO

  DO 2006 jk=klev,ilevh,-1
!IBM* ASSERT(NODEPS)
     DO 2007 ji=1,kgwd
        jl = kdx(ji)
        zhcrit(jl,jk)=MIN(ppic(jl)-pmea(jl),pmea(jl)-pval(jl))
        ll1(jl,jk)=(phgeo(jl,jk) > zhcrit(jl,jk))
        IF(ll1(jl,jk).NEQV.ll1(jl,jk+1)) THEN
           kknub(jl)=jk
        ENDIF
        IF(.NOT.ll1(jl,ilevh))kknub(jl)=ilevh
2007 END DO
2006 END DO
  !
!IBM* ASSERT(NODEPS)
  DO 2010 ji=1,kgwd
     jl = kdx(ji)
     kknu (jl)=MIN(kknu (jl),nktopg)
     kknu2(jl)=MIN(kknu2(jl),nktopg)
     kknub(jl)=MIN(kknub(jl),nktopg)
     kknul(jl)=klev
2010 END DO
  !
  !     initialize various arrays
  !
  zmair  (:)        = 0.0_wp
  prho   (:,klev+1) = 0.0_wp
  pstab  (:,1)      = 0.0_wp
  pstab  (:,klev+1) = 0.0_wp
  pri    (:,1)      = 0.0_wp
  pri    (:,klev+1) = 9999.0_wp
  pvph   (:,1)      = 0.0_wp
  pvph   (:,klev+1) = 0.0_wp
  ppsi   (:,klev+1) = 0.0_wp
  pulow  (:)        = 0.0_wp
  pvlow  (:)        = 0.0_wp
  kkcrith(:)        = klev
  kkenvh (:)        = klev
  kcrit  (:)        = 1
  ll1    (:,klev+1) = .FALSE.
  !
  !*     define flow density and stratification (rho and N2)
  !      at semi layers.
  !      -------------------------------------------------------
  !
  DO 223 jk=klev,2,-1
!IBM* NOVECTOR
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
     DO 222 ji=1,kgwd
        jl = kdx(ji)
        zdp(jl,jk)=papm1(jl,jk)-papm1(jl,jk-1)
        prho(jl,jk)=2._wp*paphm1(jl,jk)*zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
        pstab(jl,jk)=2._wp*zcons2/(ptm1(jl,jk)+ptm1(jl,jk-1))*             &
             (1._wp-cpd*prho(jl,jk)*(ptm1(jl,jk)-ptm1(jl,jk-1))/zdp(jl,jk))
        pstab(jl,jk)=MAX(pstab(jl,jk),gssec)
222  END DO
223 END DO
  !
  !********************************************************************
  !
  !*     define Low level flow (between ground and peacks-valleys)
  !      ---------------------------------------------------------
  DO 2115 jk=klev,ilevh,-1
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
     DO 2116 ji=1,kgwd
        jl = kdx(ji)
        IF(jk >= kknu2(jl).AND.jk <= kknul(jl)) THEN
           pulow(jl)        = pulow(jl)        + pum1 (jl,jk)*pmair(jl,jk)
           pvlow(jl)        = pvlow(jl)        + pvm1 (jl,jk)*pmair(jl,jk)
           pstab(jl,klev+1) = pstab(jl,klev+1) + pstab(jl,jk)*pmair(jl,jk)
           prho(jl,klev+1)  = prho (jl,klev+1) + prho (jl,jk)*pmair(jl,jk)
           zmair(jl)        = zmair(jl)        +              pmair(jl,jk)
        END IF
2116 END DO
2115 END DO

!IBM* ASSERT(NODEPS)
  DO 2110 ji=1,kgwd
     jl = kdx(ji)
     zmair_inv        = 1._wp/zmair(jl)
     pulow(jl)        = pulow(jl)        *zmair_inv
     pvlow(jl)        = pvlow(jl)        *zmair_inv
     pstab(jl,klev+1) = pstab(jl,klev+1) *zmair_inv
     prho (jl,klev+1) = prho (jl,klev+1) *zmair_inv
     !
     znorm(jl)        = MAX(SQRT(pulow(jl)**2+pvlow(jl)**2),gvsec)
     pvph (jl,klev+1) = znorm(jl)
2110 END DO
  !
  !*******  setup orography orientation relative to the low level
  !       wind and define parameters of the Anisotropic wave stress.
  !
!IBM* ASSERT(NODEPS)
  DO 2112 ji=1,kgwd
     jl = kdx(ji)
     lo=(pulow(jl) < gvsec).AND.(pulow(jl) >= -gvsec)
     IF(lo) THEN
        zu=pulow(jl)+2._wp*gvsec
     ELSE
        zu=pulow(jl)
     ENDIF
     zphi=ATAN(pvlow(jl)/zu)
     ppsi(jl,klev+1)=ptheta(jl)*pi/180._wp-zphi
     zb(jl)=1._wp-0.18_wp*pgam(jl)-0.04_wp*pgam(jl)**2
     zc(jl)=0.48_wp*pgam(jl)+0.3_wp*pgam(jl)**2
     pd1(jl)=zb(jl)-(zb(jl)-zc(jl))*(SIN(ppsi(jl,klev+1))**2)
     pd2(jl)=(zb(jl)-zc(jl))*SIN(ppsi(jl,klev+1))                    &
          *COS(ppsi(jl,klev+1))
     pdmod(jl)=SQRT(pd1(jl)**2+pd2(jl)**2)
2112 END DO
  !
  !  ************ projet flow in plane of lowlevel stress *************
  !  ************ Find critical levels...                 *************
  !
  DO 213 jk=1,klev
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
     DO 212 ji=1,kgwd
        jl = kdx(ji)
        zvt1       =pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk)
        zvt2       =-pvlow(jl)*pum1(jl,jk)+pulow(jl)*pvm1(jl,jk)
        zvpf(jl,jk)=(zvt1*pd1(jl)+zvt2*pd2(jl))/(znorm(jl)*pdmod(jl))
212  END DO
213 END DO
  ptau (:,:) = 0.0_wp
  pzdep(:,:) = 0.0_wp
  ppsi (:,:) = 0.0_wp
  ll1  (:,:) = .FALSE.

!!  DO 215 jk=2,klev-1  ! BUG FIX FOR NaN (undefined values)
  DO 215 jk=2,klev
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
     DO 214 ji=1,kgwd
        jl = kdx(ji)
        zdp (jl,jk)=  papm1 (jl,jk)-papm1 (jl,jk-1)
        pvph(jl,jk)=((paphm1(jl,jk)-papm1 (jl,jk-1))*zvpf(jl,jk)+    &
             &       (papm1 (jl,jk)-paphm1(jl,jk)  )*zvpf(jl,jk-1))  &
             &      /zdp(jl,jk)
        IF(pvph(jl,jk) < gvsec) THEN
           pvph(jl,jk)=gvsec
           IF(jk.LT.klev) kcrit(jl)=jk
        ENDIF
214  END DO
215 END DO
  !
  !*         2.3     mean flow richardson number.
  !
  DO 232 jk=2,klev
 !DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
    DO 231 ji=1,kgwd
       jl = kdx(ji)
       zdwind     = MAX(ABS(zvpf(jl,jk)-zvpf(jl,jk-1)),gvsec)
       pri(jl,jk) = pstab(jl,jk)*(zdp(jl,jk)/(grav*prho(jl,jk)*zdwind))**2
       pri(jl,jk) = MAX(pri(jl,jk),grcrit)
231  END DO
232 END DO
  !
  !*      define top of 'envelope' layer
  !       ----------------------------
  pnu(:)  = 0.0_wp
  znum(:) = 0.0_wp
  !
  DO jk=2,klev-1
!DIR$ CONCURRENT
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
!IBM* ASSERT(NODEPS)
     DO ji=1,kgwd
        jl = kdx(ji)
        IF(jk >= kknu(jl)) THEN
           znum(jl)= pnu(jl)
           zwind   = (pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/        &
                &    MAX(SQRT(pulow(jl)**2+pvlow(jl)**2),gvsec)
           zwind   = MAX(SQRT(zwind**2),gvsec)
           zstabm  = SQRT(MAX(pstab(jl,jk  ),gssec))
           zstabp  = SQRT(MAX(pstab(jl,jk+1),gssec))
           zrhom   = prho(jl,jk  )
           zrhop   = prho(jl,jk+1)
           pnu(jl) = pnu(jl) + pmair(jl,jk)*((zstabp/zrhop+zstabm/zrhom)/2._wp)/zwind
           IF ((znum(jl) <= gfrcrit).AND.(pnu(jl) > gfrcrit).AND.(kkenvh(jl) == klev))  kkenvh(jl)=jk
        ENDIF
     END DO
  END DO
  !
  !  calculation of a dynamical mixing height for when the waves
  !  BREAK AT LOW LEVEL: The drag will be repartited over
  !  a depths that depends on waves vertical wavelength,
  !  not just between two adjacent model layers.
  !  of gravity waves:
  !
  znup(:) = 0.0_wp
  znum(:) = 0.0_wp
  !
  DO jk=klev-1,2,-1
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
    DO ji=1,kgwd
       jl = kdx(ji)
       znum(jl)= znup(jl)
       zwind   = (pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/        &
            &    MAX(SQRT(pulow(jl)**2+pvlow(jl)**2),gvsec)
       zwind   = MAX(SQRT(zwind**2),gvsec)
       zstabm  = SQRT(MAX(pstab(jl,jk  ),gssec))
       zstabp  = SQRT(MAX(pstab(jl,jk+1),gssec))
       zrhom   = prho(jl,jk  )
       zrhop   = prho(jl,jk+1)
       znup(jl)= znup(jl) + pmair(jl,jk)*((zstabp/zrhop+zstabm/zrhom)/2._wp)/zwind
       IF ((znum(jl) <= pi/2._wp).AND.(znup(jl) > pi/2._wp).AND.(kkcrith(jl) == klev))  kkcrith(jl)=jk
    END DO
  END DO

!IBM* ASSERT(NODEPS)
  DO  ji=1,kgwd
     jl = kdx(ji)
     kkcrith(jl)=min0(kkcrith(jl),kknu(jl))
     kkcrith(jl)=max0(kkcrith(jl),ilevh*2)
     IF(kcrit(jl) >= kkcrith(jl))kcrit(jl)=1
  END DO
  !
  !     directional info for flow blocking *************************
  !
  DO 251 jk=1,klev
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
     DO 252 ji=1,kgwd
        jl = kdx(ji)
        lo=(pum1(jl,jk) < gvsec).AND.(pum1(jl,jk) >= -gvsec)
        IF(lo) THEN
           zu=pum1(jl,jk)+2._wp*gvsec
        ELSE
           zu=pum1(jl,jk)
        ENDIF
        zphi=ATAN(pvm1(jl,jk)/zu)
        ppsi(jl,jk)=ptheta(jl)*pi/180._wp-zphi
252  END DO
251 END DO

  !      forms the vertical 'leakiness' **************************

  DO 254  jk=ilevh,klev
!IBM* ASSERT(NODEPS)
     DO 253  ji=1,kgwd
        jl = kdx(ji)
        pzdep(jl,jk)=0._wp
        IF(jk >= kkenvh(jl).AND.kkenvh(jl) /= klev) THEN
           pzdep(jl,jk)=(phgeo(jl,kkenvh(jl)-1)-phgeo(jl,jk))/      &
                &       (phgeo(jl,kkenvh(jl)  )-phgeo(jl,klev))
        END IF
253  END DO
254 END DO

END SUBROUTINE orosetup

SUBROUTINE gwstress( jg, kbdim,  klev,              &
     &               kgwd,   kdx,    kkenvh,        &
     &               pstd,   psig,   ppic,   pval,  &
     &               prho,   pstab,  pvph,   pdmod, &
     &               ptau0 )
  !
  !**** *gwstress*
  !
  !     purpose.
  !     --------
  !  Compute the surface stress due to Gravity Waves, according
  !  to the Phillips (1979) theory of 3-D flow above
  !  anisotropic elliptic ridges.

  !  The stress is reduced two account for cut-off flow over
  !  hill.  The flow only see that part of the ridge located
  !  above the blocked layer (see zeff).
  !

  !     called from orodrag*
  !
  !
  !     reference.
  !     ----------
  !
  !   LOTT and MILLER (1997)  &  LOTT (1999)
  !
  !     author.
  !     -------
  !     f. lott put the new gwd on ifs      22/11/93

  IMPLICIT NONE

  INTEGER,  INTENT(in)  :: jg, kbdim, klev
  INTEGER,  INTENT(in)  :: kgwd, kdx(kbdim)
  INTEGER,  INTENT(in)  :: kkenvh(kbdim)
  !
  REAL(wp), INTENT(in)  :: pstd(kbdim), psig(kbdim), ppic(kbdim), pval(kbdim), pdmod(kbdim)
  REAL(wp), INTENT(in)  :: prho(kbdim,klev+1), pstab(kbdim,klev+1), pvph(kbdim,klev+1)
  !
  REAL(wp), INTENT(out) :: ptau0(kbdim)

  !  Local scalars:
  INTEGER  :: jl, ji
  REAL(wp) :: zeff  ! effective height seen by the flow when there is blocking

  ! Shortcuts to components of echam_sso_config
  !
  REAL(wp), POINTER :: gkdrag
  !
  gkdrag => echam_sso_config(jg)% gkdrag

  !
  !*         1.1     gravity wave stress.
  !

  ptau0(:) = 0.0_wp

!IBM* ASSERT(NODEPS)
  DO 301 ji=1,kgwd
     jl = kdx(ji)

     !  effective mountain height above the blocked flow

     zeff=ppic(jl)-pval(jl)
     IF(kkenvh(jl) < klev)THEN
        zeff=MIN(gfrcrit*pvph(jl,klev+1)/SQRT(pstab(jl,klev+1)),zeff)
     ENDIF

     ptau0(jl) = gkdrag                                    &
          &     *prho(jl,klev+1)                           &
          &     *psig(jl)*pdmod(jl)/4._wp/pstd(jl)         &
          &     *pvph(jl,klev+1)*SQRT(pstab(jl,klev+1))    &
          &     *zeff**2


        !  too small value of stress or  low level flow include critical level
        !  or low level flow:  gravity wave stress nul.

        !       lo=(ptau0(jl).lt.gtsec).or.(kcrit(jl).ge.kknu(jl))
        !    *      .or.(pvph(jl,klev+1).lt.gvcrit)
        !       if(lo) ptau0(jl)=0.0_wp

301 END DO

END SUBROUTINE gwstress


SUBROUTINE gwprofil( jg, kbdim,  klev,                                  &
     &               kgwd,   kdx,                                       &
     &               kkcrith,kcrit,                                     &
     &               pstd,   psig,                                      &
     &               paphm1, prho,   pri,    pstab,  pvph,   pdmod,     &
     &               ptau0,  ptau )

  ! *gwprofil called from orodrag

  !     method:
  !     -------
  !     the stress profile for gravity waves is computed as follows:
  !     it decreases linearly with heights from the ground
  !     to the low-level indicated by kkcrith,
  !     to simulates lee waves or
  !     low-level gravity wave breaking.
  !     above it is constant, except when the waves encounter a critical
  !     level (kcrit) or when they break.
  !     The stress is also uniformly distributed above the level
  !     ntop.
  !

  IMPLICIT NONE

  INTEGER,  INTENT(in)  :: jg, kbdim, klev
  INTEGER,  INTENT(in)  :: kgwd        ! Number  of points where oro scheme is active
  INTEGER,  INTENT(in)  :: kdx(kbdim)  ! Indices of points where oro scheme is active
  INTEGER,  INTENT(in)  :: kkcrith(kbdim), kcrit(kbdim)
  !
  REAL(wp), INTENT(in)  :: pdmod (kbdim), psig(kbdim), pstd(kbdim)
  !
  REAL(wp), INTENT(in)  :: paphm1(kbdim,klev+1), pstab(kbdim,klev+1),    &
       &                   prho  (kbdim,klev+1), pvph (kbdim,klev+1),    &
       &                   pri   (kbdim,klev+1)
  !
  REAL(wp), INTENT(in)  :: ptau0(kbdim)
  REAL(wp), INTENT(out) :: ptau(kbdim,klev+1)

  !  Local scalars:
  INTEGER  :: ji, jl, jk
  REAL(wp) :: zsqr, zalfa, zriw, zdel, zb, zalpha, zdz2n, zdelp, zdelpt

  !  Local arrays:
  REAL(wp) :: zdz2 (kbdim,klev) , znorm(kbdim) , zoro(kbdim), ztau (kbdim,klev+1)

  ! Shortcuts to components of echam_sso_config
  !
  INTEGER , POINTER :: ntop
  !
  ntop   => echam_sso_config(jg)% ntop

  !
  !  Executable statements

  ptau(:,:) = 0.0_wp
  
!CDIR NODEP
!DIR$ CONCURRENT
  DO 400 ji=1,kgwd
     jl=kdx(ji)
     zoro(jl)=psig(jl)*pdmod(jl)/4._wp/pstd(jl)
     ztau(jl,klev+1)=ptau0(jl)
     ztau(jl,kkcrith(jl))=grahilo*ptau0(jl)
400 END DO

  !
  DO 430 jk=klev+1,2,-1
     !
     !
     !*         4.1    constant shear stress until top of the
     !                 low-level breaking/trapped layer

     !
!DIR$ CONCURRENT
!DIR$ PREFERVECTOR
!DIR$ PREFERSTREAM
     DO 411 ji=1,kgwd
        jl=kdx(ji)
        IF(jk > kkcrith(jl)) THEN
           zdelp=paphm1(jl,jk)-paphm1(jl,klev+1)
           zdelpt=paphm1(jl,kkcrith(jl))-paphm1(jl,klev+1)
           ptau(jl,jk)=ztau(jl,klev+1)+zdelp/zdelpt*                    &
                (ztau(jl,kkcrith(jl))-ztau(jl,klev+1))
        ELSE
           ptau(jl,jk)=ztau(jl,kkcrith(jl))
        ENDIF
411  END DO

430 END DO

  !
  !*         4.15   constant shear stress until the top of the
  !                 low level flow layer.

  !
  !
  !*         4.2    wave displacement at next level.
  !
  !

  !
  !*         4.4    wave richardson number, new wave displacement
  !*                and stress:  breaking evaluation and critical
  !                 level
  !

  DO 440 jk=klev,2,-1
!DIR$ CONCURRENT
     DO 441 ji=1,kgwd
        jl=kdx(ji)
        znorm(jl)=prho(jl,jk)*SQRT(pstab(jl,jk))*pvph(jl,jk)
        zdz2(jl,jk)=ptau(jl,jk)/MAX(znorm(jl),gssec)/zoro(jl)
441  END DO

!CDIR NODEP
!DIR$ CONCURRENT
     DO 442 ji=1,kgwd
        jl=kdx(ji)
        IF(jk < kkcrith(jl)) THEN
           IF((ptau(jl,jk+1) < gtsec).OR.(jk <= kcrit(jl))) THEN
              ptau(jl,jk)=0.0_wp
           ELSE
              zsqr=SQRT(pri(jl,jk))
              zalfa=SQRT(pstab(jl,jk)*zdz2(jl,jk))/pvph(jl,jk)
              zriw=pri(jl,jk)*(1._wp-zalfa)/(1+zalfa*zsqr)**2
              IF(zriw < grcrit) THEN
                 zdel=4._wp/zsqr/grcrit+1._wp/grcrit**2+4._wp/grcrit
                 zb=1._wp/grcrit+2._wp/zsqr
                 zalpha=0.5_wp*(-zb+SQRT(zdel))
                 zdz2n=(pvph(jl,jk)*zalpha)**2/pstab(jl,jk)
                 ptau(jl,jk)=znorm(jl)*zdz2n*zoro(jl)
              ENDIF

              ptau(jl,jk)=MIN(ptau(jl,jk),ptau(jl,jk+1))

           ENDIF
        ENDIF
442  END DO
440 END DO



  !  REORGANISATION OF THE STRESS PROFILE AT LOW LEVEL

!DIR$ CONCURRENT
  DO 530 ji=1,kgwd
     jl=kdx(ji)
     ztau(jl,kkcrith(jl))=ptau(jl,kkcrith(jl))
     ztau(jl,ntop)=ptau(jl,ntop)
530 END DO

  DO 531 jk=1,klev
!DIR$ CONCURRENT
     DO 532 ji=1,kgwd
        jl=kdx(ji)

        IF(jk > kkcrith(jl))THEN

           zdelp=paphm1(jl,jk)-paphm1(jl,klev+1    )
           zdelpt=paphm1(jl,kkcrith(jl))-paphm1(jl,klev+1    )
           ptau(jl,jk)=ztau(jl,klev+1    ) +                             &
                (ztau(jl,kkcrith(jl))-ztau(jl,klev+1    ) )*      &
                zdelp/zdelpt

        ENDIF

532  END DO

     !  REORGANISATION AT THE MODEL TOP....

!DIR$ CONCURRENT
     DO 533 ji=1,kgwd
        jl=kdx(ji)

        IF (jk < ntop) THEN

           !         zdelp =paphm1(jl,ntop)
           !         zdelpt=paphm1(jl,jk)
           !         ptau(jl,jk)=ztau(jl,ntop)*zdelpt/zdelp
           ptau(jl,jk)=ztau(jl,ntop)

        ENDIF

533  END DO


531 END DO


END SUBROUTINE gwprofil

SUBROUTINE orolift( jg, kproma, kbdim, klev,  &
  &                 pcoriol,                  &
  &                 pdtime,                   &
  &                 ktest,                    &
  &                 phgeo,  paphm1,           &
  &                 pmair,                    &
  &                 ptm1,   pum1,  pvm1,      &
  &                 pmea,   pstd,  ppic,      &
  &                 pvom,   pvol,  pdis )
  !
  !**** *orolift: simulate the geostrophic lift.
  !
  !     purpose.
  !     --------
  ! this routine computes the physical tendencies of the
  ! prognostic variables u,v  when enhanced vortex stretching
  ! is needed.
  !
  !          called from *ssodrag*.
  !     author.
  !     -------
  !     f.lott  lmd 22/11/95
  !
  !

  IMPLICIT NONE

  ! scalar arguments with intent(IN):
  INTEGER,  INTENT(in)  :: jg, kproma, kbdim, klev
  REAL(wp), INTENT(in)  :: pdtime

  ! array arguments with intent(IN):
  ! Input 1D
  INTEGER,  INTENT(in)  :: ktest(kbdim)  ! Flags to indicate active points

  REAL(wp), INTENT(in)  :: pmea(kbdim)   ! Mean Orography (m)
  REAL(wp), INTENT(in)  :: pstd(kbdim)   ! SSO standard deviation (m)
  REAL(wp), INTENT(in)  :: ppic(kbdim)   ! SSO Peacks elevation (m)
  
  REAL(wp), INTENT(in)  :: pcoriol(kbdim)       ! Coriolis parameter (1/s)
  ! Input 2D
  REAL(wp), INTENT(in)  :: paphm1(kbdim,klev+1) ! half level pressure (t-dt)
  REAL(wp), INTENT(in)  :: pmair(kbdim,klev)    ! air mass (kg/m2)
  REAL(wp), INTENT(in)  :: phgeo(kbdim,klev)    ! geopotential height above surface
  REAL(wp), INTENT(in)  :: ptm1(kbdim,klev)     ! temperature (t-dt)
  REAL(wp), INTENT(in)  :: pum1(kbdim,klev)     ! zonal wind (t-dt)
  REAL(wp), INTENT(in)  :: pvm1(kbdim,klev)     ! meridional wind (t-dt)


  ! array arguments with intent(OUT):
  ! Input 2D
  REAL(wp), INTENT(out) :: pdis(kbdim,klev)    ! energy dissipation rate per mass
  REAL(wp), INTENT(out) :: pvol(kbdim,klev)    ! tendency of meridional wind
  REAL(wp), INTENT(out) :: pvom(kbdim,klev)    ! tendency of zonal wind

  ! Local scalars:
  INTEGER  :: jl, jk
  REAL(wp) :: zslow, zsqua, zscav, zbet
  REAL(wp) :: zcons1
  LOGICAL  :: lifthigh

  ! Local arrays:
  INTEGER  :: iknub(kbdim), iknul(kbdim)
  REAL(wp) :: zdudt(kbdim), zdvdt(kbdim)

  REAL(wp) :: pulow(kbdim), pvlow(kbdim)
  REAL(wp) :: ztau(kbdim,klev+1), ztav(kbdim,klev+1)
  REAL(wp) :: zrho(kbdim,klev+1), zmair(kbdim), zmair_inv
  REAL(wp) :: zhcrit(kbdim,klev)
  LOGICAL  :: ll1(kbdim,klev+1)

  ! Shortcuts to components of echam_sso_config
  !
  REAL(wp), POINTER :: gklift
  INTEGER , POINTER :: nktopg
  !
  nktopg => echam_sso_config(jg)% nktopg
  gklift => echam_sso_config(jg)% gklift

  !-----------------------------------------------------------------------
  !
  !*         1.1  initialisations
  !               ---------------

  lifthigh=.FALSE.

  zcons1=1._wp/rd

  iknub(:)        = klev
  iknul(:)        = klev
  ll1  (:,klev+1) = .FALSE.

  zrho (:,klev+1) = 0.0_wp
  zmair(:)        = 0.0_wp
  pulow(:)        = 0.0_wp
  pvlow(:)        = 0.0_wp

  pvom (:,:)      = 0.0_wp
  pvol (:,:)      = 0.0_wp
  pdis (:,:)      = 0.0_wp

  !
  !*         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
  !*                 LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
  !*                 THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.
  !

  DO 2006 jk=klev,1,-1
     DO 2007 jl=1,kproma
        IF(ktest(jl) == 1) THEN
           zhcrit(jl,jk)=MAX(ppic(jl)-pmea(jl),100.0_wp)
           ll1(jl,jk)=(phgeo(jl,jk) > zhcrit(jl,jk))
           IF(ll1(jl,jk).NEQV.ll1(jl,jk+1)) THEN
              iknub(jl)=jk
           ENDIF
        ENDIF
2007 END DO
2006 END DO

  DO 2010 jl=1,kproma
     IF(ktest(jl) == 1) THEN
        iknub(jl)=MAX(iknub(jl),klev/2)
        iknul(jl)=MAX(iknul(jl),2*klev/3)
        IF(iknub(jl) >  nktopg) iknub(jl)=nktopg
        IF(iknub(jl) == nktopg) iknul(jl)=klev
        IF(iknub(jl) == iknul(jl)) iknub(jl)=iknul(jl)-1
     ENDIF
2010 END DO

  DO 223 jk=klev,2,-1
     DO 222 jl=1,kproma
        zrho(jl,jk)=2._wp*paphm1(jl,jk)*zcons1/(ptm1(jl,jk)+ptm1(jl,jk-1))
222  END DO
223 END DO

  !     print *,'  dans orolift: 223'
  !********************************************************************
  !
  !*     define low level flow
  !      -------------------
  DO 2115 jk=klev,1,-1
     DO 2116 jl=1,kproma
        IF(ktest(jl) == 1) THEN
           IF(jk >= iknub(jl).AND.jk <= iknul(jl)) THEN
              pulow(jl)        = pulow(jl)        + pum1(jl,jk)*pmair(jl,jk)
              pvlow(jl)        = pvlow(jl)        + pvm1(jl,jk)*pmair(jl,jk)
              zrho (jl,klev+1) = zrho (jl,klev+1) + zrho(jl,jk)*pmair(jl,jk)
              zmair(jl)        = zmair(jl)        +             pmair(jl,jk)
           ENDIF
        ENDIF
2116 END DO
2115 END DO

  DO 2110 jl=1,kproma
     IF(ktest(jl) == 1) THEN
        zmair_inv        = 1._wp/zmair(jl)
        pulow(jl)        = pulow(jl)        *zmair_inv
        pvlow(jl)        = pvlow(jl)        *zmair_inv
        zrho (jl,klev+1) = zrho (jl,klev+1) *zmair_inv
     ENDIF
2110 END DO



  !***********************************************************
  !
  !*         3.      COMPUTE MOUNTAIN LIFT
  !
  DO 301 jl=1,kproma
     IF(ktest(jl) == 1) THEN
        ztau(jl,klev+1)= - gklift                       &
             &            *zrho(jl,klev+1)*pcoriol(jl)  &
            !&            *(2._wp*pstd(jl)+pmea(jl))    &
             &            * 2._wp*pstd(jl)              &
             &            * pvlow(jl)
        ztav(jl,klev+1)=   gklift                       &
             &            *zrho(jl,klev+1)*pcoriol(jl)  &
            !&            *(2._wp*pstd(jl)+pmea(jl))    &
             &            * 2._wp*pstd(jl)              &
             &            * pulow(jl)
     ELSE
        ztau(jl,klev+1)=0.0_wp
        ztav(jl,klev+1)=0.0_wp
     ENDIF
301 END DO
  !
  !*         4.      COMPUTE LIFT PROFILE
  !*                 --------------------
  !

  DO 401 jk=1,klev
     DO 402 jl=1,kproma
        IF(ktest(jl) == 1) THEN
           ztau(jl,jk)=ztau(jl,klev+1)*paphm1(jl,jk)/paphm1(jl,klev+1)
           ztav(jl,jk)=ztav(jl,klev+1)*paphm1(jl,jk)/paphm1(jl,klev+1)
        ELSE
           ztau(jl,jk)=0.0_wp
           ztav(jl,jk)=0.0_wp
        ENDIF
402  END DO
401 END DO
  !
  !
  !*         5.      COMPUTE TENDENCIES.
  !*                 -------------------
  IF(lifthigh)THEN


     !
     !  EXPLICIT SOLUTION AT ALL LEVELS
     !
     DO 524 jk=1,klev
        DO 523 jl=1,kproma
           IF(ktest(jl) == 1) THEN
              zmair_inv = 1._wp/pmair(jl,jk)
              zdudt(jl) = -(ztau(jl,jk+1)-ztau(jl,jk))*zmair_inv
              zdvdt(jl) = -(ztav(jl,jk+1)-ztav(jl,jk))*zmair_inv
           ENDIF
523     END DO
524  END DO
     !
     !  PROJECT PERPENDICULARLY TO U NOT TO DESTROY ENERGY
     !
     DO 530 jk=1,klev
        DO 531 jl=1,kproma
           IF(ktest(jl) == 1) THEN
              zslow=SQRT(pulow(jl)**2+pvlow(jl)**2)
              zsqua=MAX(SQRT(pum1(jl,jk)**2+pvm1(jl,jk)**2),gvsec)
              zscav=-zdudt(jl)*pvm1(jl,jk)+zdvdt(jl)*pum1(jl,jk)
              IF(zsqua > gvsec)THEN
                 pvom(jl,jk)=-zscav*pvm1(jl,jk)/zsqua**2
                 pvol(jl,jk)= zscav*pum1(jl,jk)/zsqua**2
              ELSE
                 pvom(jl,jk)=0.0_wp
                 pvol(jl,jk)=0.0_wp
              ENDIF
              zsqua=SQRT(pum1(jl,jk)**2+pum1(jl,jk)**2)
              IF(zsqua < zslow)THEN
                 pvom(jl,jk)=zsqua/zslow*pvom(jl,jk)
                 pvol(jl,jk)=zsqua/zslow*pvol(jl,jk)
              ENDIF
           ENDIF
531     END DO
530  END DO
     !
     !  6.  LOW LEVEL LIFT, SEMI IMPLICIT:
     !  ----------------------------------

  ELSE

     DO 601 jl=1,kproma
        IF(ktest(jl) == 1) THEN
           DO jk=klev,iknub(jl),-1
              zbet =  gklift*pcoriol(jl)*pdtime               &
                   & *(phgeo(jl,iknub(jl)-1)-phgeo(jl,  jk))  &
                   & /(phgeo(jl,iknub(jl)-1)-phgeo(jl,klev))
              zdudt(jl)=-pum1(jl,jk)/pdtime/(1+zbet**2)
              zdvdt(jl)=-pvm1(jl,jk)/pdtime/(1+zbet**2)
              pvom(jl,jk)= zbet**2*zdudt(jl) - zbet   *zdvdt(jl)
              pvol(jl,jk)= zbet   *zdudt(jl) + zbet**2*zdvdt(jl)
           ENDDO
        ENDIF
601  END DO

  ENDIF
  !           PRINT *,' out orolift'
END SUBROUTINE orolift

END MODULE mo_ssodrag
