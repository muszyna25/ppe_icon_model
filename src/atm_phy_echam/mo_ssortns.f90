!>
!! @brief  Does the SSO drag following LOTT &MILLER (1997)
!!
!!  !!!!!!!ATTENTION! This scheme is already part of ICON in its version which has been run
!!  in GME and the COSMO model at DWD. A comparison between the two version is needed
!!  and a decision has to be made to stick to one version.!!!!!!!
!!
!! @author  m.miller + b.ritter   e.c.m.w.f.     15/06/1986.
!! @author f.lott + m. miller    e.c.m.w.f.     22/11/1994
!! 
!!
!! @par Revision History
!! @par  e.manzini    mpi  05.09.2000, MPI-M
!! @par K. Froehlich, MPI  2011/10/04 first implementation into ICON to be compared 
!!                                     with DWD version of SSO
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
MODULE mo_ssortns

  USE mo_kind,              ONLY: wp
  USE mo_ssodrag,           ONLY: sugwd, nktopg, ntop, gstd, gkdrag,gkwake,gfrcrit,grcrit ,&
    &                             gklift, gsigcr, gssec, gtsec, gvsec, grahilo, gpicmea
  USE mo_math_constants,    ONLY: pi    
  USE mo_physical_constants,ONLY: grav, cpd,rd, omega
  USE mo_ssodrag

!  USE mo_time_control,    ONLY: time_step_len


  IMPLICIT NONE

CONTAINS

SUBROUTINE ssodrag ( kproma,  kbdim,   klev,  krow,                   &
                     paphm1,  papm1,   pgeom1,                        &
                     ptm1,    pum1,    pvm1,                          &
                     pmea,    pstd,    psig,  pgam, pthe, ppic, pval, &
                     pustrgw, pvstrgw, pvdisgw,                       &
                     ptte,     pvol,    pvom, time_step_len)
  !
  ! Description: see above
  !
  !   *ssodrag* is called from *mo_echam_phy_main*
  !
  ! for more details see file AUTHORS
  !
 

  ! scalar arguments with intent(IN):
  INTEGER, INTENT(in) :: kproma, kbdim, klev, krow
  REAL(wp), INTENT(in) :: time_step_len

  ! array arguments with intent(IN):
  ! Input 1D
  REAL(wp), INTENT(in) :: pmea(kbdim)   ! Mean Orography (m)
  REAL(wp), INTENT(in) :: pstd(kbdim)   ! SSO standard deviation (m)
  REAL(wp), INTENT(in) :: psig(kbdim)   ! SSO slope
  REAL(wp), INTENT(in) :: pgam(kbdim)   ! SSO Anisotropy
  REAL(wp), INTENT(in) :: pthe(kbdim)   ! SSO Angle
  REAL(wp), INTENT(in) :: ppic(kbdim)   ! SSO Peacks elevation (m)
  REAL(wp), INTENT(in) :: pval(kbdim)   ! SSO Valleys elevation (m)
  ! Input 2D
  REAL(wp), INTENT(in) :: paphm1(kbdim,klev+1) ! half level pressure (t-dt)
  REAL(wp), INTENT(in) :: papm1(kbdim,klev)    ! full level pressure (t-dt)
  REAL(wp), INTENT(in) :: pgeom1(kbdim,klev)   ! geopotential above surface (t-dt)
  REAL(wp), INTENT(in) :: ptm1(kbdim,klev)     ! temperature (t-dt)
  REAL(wp), INTENT(in) :: pum1(kbdim,klev)     ! zonal wind (t-dt)
  REAL(wp), INTENT(in) :: pvm1(kbdim,klev)     ! meridional wind (t-dt)

  ! array arguments with intent(INOUT):
  ! Input 1D
  REAL(wp), INTENT(inout) :: pustrgw(kbdim)    ! u-gravity wave stress
  REAL(wp), INTENT(inout) :: pvstrgw(kbdim)    ! v-gravity wave stress
  REAL(wp), INTENT(inout) :: pvdisgw(kbdim)    ! dissipation by gravity wave drag
  ! Input 2D
  REAL(wp), INTENT(inout) :: ptte(kbdim,klev) ! tendency of temperature
  REAL(wp), INTENT(inout) :: pvol(kbdim,klev) ! tendency of meridional wind
  REAL(wp), INTENT(inout) :: pvom(kbdim,klev) ! tendency of zonal wind

  ! Local scalars:
  INTEGER :: igwd, jk, jl, ji

  ! Local arrays:
  INTEGER :: idx(kbdim), itest(kbdim)
  REAL(wp)    :: zdu_oro(kbdim,klev) ! tendency due to ORO GW DRAG  (m/s)
  REAL(wp)    :: zdv_oro(kbdim,klev) ! tendency due to ORO GW DRAG  (m/s)
  REAL(wp)    :: zdt_oro(kbdim,klev) ! tendency due to ORO GW DRAG  (K)
  REAL(wp)    :: zdu_lif(kbdim,klev) ! tendency due to MOUNTAIN LIFT(m/s)
  REAL(wp)    :: zdv_lif(kbdim,klev) ! tendency due to MOUNTAIN LIFT(m/s)
  REAL(wp)    :: zdt_lif(kbdim,klev) ! tendency due to MOUNTAIN LIFT(K)

  !
  !*         1.    initialization
  !                --------------
  !  INITIALIZE CONSTANT FOR THE GWD SCHEME
  !  ON-LINE SHOULD ONLY BE CALLED AT THE MODEL START UP.

  CALL sugwd(kproma,klev)

  pustrgw(:) = 0.0_wp
  pvstrgw(:) = 0.0_wp
  pvdisgw(:) = 0.0_wp

  zdu_oro(:,:) = 0.0_wp
  zdv_oro(:,:) = 0.0_wp
  zdt_oro(:,:) = 0.0_wp

  zdu_lif(:,:) = 0.0_wp
  zdv_lif(:,:) = 0.0_wp
  zdt_lif(:,:) = 0.0_wp

  idx(:) = 0

  !
  !*         2.    orographic gravity wave drag
  !                -----------------------------

  !  SELECTION  POINTS WHERE THE SCHEME IS ACTIVE

  igwd=0
  DO jl=1,kproma
     itest(jl)=0
     IF (((ppic(jl)-pmea(jl)) > gpicmea).AND.(pstd(jl) > gstd)) THEN
        itest(jl)=1
        igwd=igwd+1
        idx(igwd)=jl
     ENDIF
  ENDDO

  !
  CALL orodrag( kproma,  kbdim,   klev,                           &
                igwd,             idx,   itest,                   &
                paphm1,  papm1,   pgeom1,                         &
                ptm1,    pum1,    pvm1,                           &
                pmea,    pstd,    psig,  pgam, pthe, ppic, pval,  &
                zdu_oro, zdv_oro, zdt_oro, time_step_len)
  !
  !*         3.    mountain lift
  !                --------------
  CALL orolift( kproma, kbdim, klev,krow,                          &
                itest,                                             &
                paphm1,  pgeom1,                                   &
                ptm1,    pum1,    pvm1,                            &
                pmea,    pstd,    ppic,                            &
                zdu_lif, zdv_lif, zdt_lif )

  ! STRESS FROM TENDENCIES

  DO jk = 1, klev
!CDIR NODEP
     DO jl = 1, igwd
        ji=idx(jl)
        pustrgw(ji) = pustrgw(ji)                                       &
             +(zdu_oro(ji,jk)+zdu_lif(ji,jk))                           &
             *(paphm1(ji,jk+1)-paphm1(ji,jk))/grav
        pvstrgw(ji) = pvstrgw(ji)                                       &
             +(zdv_oro(ji,jk)+zdv_lif(ji,jk))                           &
             *(paphm1(ji,jk+1)-paphm1(ji,jk))/grav
        pvdisgw(ji) = pvdisgw(ji)                                       &
             +(zdt_oro(ji,jk)+zdt_lif(ji,jk))                           &
             *(paphm1(ji,jk+1)-paphm1(ji,jk))/grav*cpd
     ENDDO
  ENDDO
  !
  !*         4.    total quantities
  !                ----------------

  do jk=1,klev
!CDIR NODEP
    do jl=1,igwd
      ji=idx(jl)
      ptte(ji,jk) =  ptte(ji,jk)+zdt_oro(ji,jk)+zdt_lif(ji,jk)
      pvol(ji,jk) =  pvol(ji,jk)+zdv_oro(ji,jk)+zdv_lif(ji,jk)
      pvom(ji,jk) =  pvom(ji,jk)+zdu_oro(ji,jk)+zdu_lif(ji,jk)
    enddo
  enddo

  RETURN
END SUBROUTINE ssodrag

SUBROUTINE orodrag( kproma,  kbdim,  klev,            &
     kgwd,           kdx,    ktest,                   &
     paphm1, papm1,  pgeom1,                          &
     ptm1,   pum1,   pvm1,                            &
     pmea,   pstd,   psig,   pgam, pthe, ppic, pval,  &
     pvom,   pvol,   pte, time_step_len )
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
  INTEGER, INTENT(in) :: kproma, kbdim, klev
  INTEGER, INTENT(in) :: kgwd      ! Total points where oro scheme is active
  REAL(wp),INTENT(in) :: time_step_len

  ! array arguments with intent(IN):
  ! Input 1D
  INTEGER, INTENT(in) :: kdx(kbdim)   ! physical location of an active point
  INTEGER, INTENT(in) :: ktest(kbdim) ! Flags to indicate active points

  REAL(wp), INTENT(in) :: pmea(kbdim)   ! Mean Orography (m)
  REAL(wp), INTENT(in) :: pstd(kbdim)   ! SSO standard deviation (m)
  REAL(wp), INTENT(in) :: psig(kbdim)   ! SSO slope
  REAL(wp), INTENT(in) :: pgam(kbdim)   ! SSO Anisotropy
  REAL(wp), INTENT(in) :: pthe(kbdim)   ! SSO Angle
  REAL(wp), INTENT(in) :: ppic(kbdim)   ! SSO Peacks elevation (m)
  REAL(wp), INTENT(in) :: pval(kbdim)   ! SSO Valleys elevation (m)
  ! Input 2D
  REAL(wp), INTENT(in) :: paphm1(kbdim,klev+1) ! half level pressure (t-dt)
  REAL(wp), INTENT(in) :: papm1(kbdim,klev)    ! full level pressure (t-dt)
  REAL(wp), INTENT(in) :: pgeom1(kbdim,klev)   ! geopotential above surface (t-dt)
  REAL(wp), INTENT(in) :: ptm1(kbdim,klev)     ! temperature (t-dt)
  REAL(wp), INTENT(in) :: pum1(kbdim,klev)     ! zonal wind (t-dt)
  REAL(wp), INTENT(in) :: pvm1(kbdim,klev)     ! meridional wind (t-dt)

  ! array arguments with intent(INOUT):
  ! Input 2D
  REAL(wp), INTENT(inout) :: pte(kbdim,klev) ! tendency of temperature
  REAL(wp), INTENT(inout) :: pvol(kbdim,klev) ! tendency of meridional wind
  REAL(wp), INTENT(inout) :: pvom(kbdim,klev) ! tendency of zonal wind

  !  Local scalars:
  INTEGER :: klevm1, jl, jk, ji
  REAL(wp)    :: ztmst, zdelp, ztemp, zforc, ztend, rover, &
       zb, zc, zconb, zabsv, zzd1, ratio, zbet,        &
       zust, zvst, zdis, zred

  !  Local arrays:
  INTEGER  ::    icrit(kbdim),     ikcrith(kbdim), &
       ikenvh(kbdim),   iknu(kbdim),      iknu2(kbdim),          &
       ikcrit(kbdim)

  REAL(wp) ::  zdudt(kbdim),    zdvdt(kbdim),     ztfr(kbdim),  &
       znu(kbdim),      zd1(kbdim),       zd2(kbdim),       &
       zdmod(kbdim)

  REAL(wp) ::  ztau(kbdim,klev+1),  zstab(kbdim,klev+1),  zvph(kbdim,klev+1),  &
       zrho(kbdim,klev+1),   zri(kbdim,klev+1),    zpsi(kbdim,klev+1),     &
       zzdep(kbdim,klev),   pulow(kbdim),         pvlow(kbdim)

  !  Executable statements

  !*         1.1   computational constants
  !*                ------------ ----------

  ztmst = time_step_len

  klevm1=klev-1

  !     ------------------------------------------------------------------
  !
  !*         2.     precompute basic state variables.
  !*                ---------- ----- ----- ----------
  !*                define low level wind, project winds in plane of
  !*                low level wind, determine sector in which to take
  !*                the variance and set indicator for critical levels.
  !
  CALL orosetup( kproma,  kbdim,  klev,  kgwd, kdx, ktest,             &
       ikcrit, ikcrith, icrit,   ikenvh, iknu,   iknu2,                &
       paphm1, papm1,   pum1,    pvm1,  ptm1,   pgeom1,                &
       zrho  , zri,     zstab,   ztau,  zvph,   zpsi,   zzdep,         &
       pulow,  pvlow,                                                  &
       pthe,   pgam,  pmea, ppic, pval, znu, zd1, zd2,  zdmod )
  !
  !***********************************************************
  !
  !*         3.      compute low level stresses using subcritical and
  !*                 supercritical forms.computes anisotropy coefficient
  !*                 as measure of orographic twodimensionality.
  !
  CALL gwstress( kproma, kbdim, klev, kgwd,                            &
       kdx, ktest,  icrit,  ikenvh, iknu,                              &
       zrho,   zstab,  zvph,   pstd,  psig,    ppic, pval,             &
       ztau,   zdmod)
  !
  !*         4.      compute stress profile including
  !                  trapped waves, wave breaking,
  !                  linear decay in stratosphere.
  !
  CALL gwprofil(  kproma, kbdim, klev,                         &
       kgwd,    kdx  ,                                         &
       ikcrith, icrit,                                         &
       paphm1,  zrho,  zstab,  ztfr, zvph,                     &
       zri,     ztau,                                          &
       zdmod,   psig,  pstd)
  !
  !*         5.      Compute tendencies from waves stress profile.
  !                  Compute low level blocked flow drag.
  !*                 --------------------------------------------
  !
  !  explicit solution at all levels for the gravity wave
  !  implicit solution for the blocked levels
  !
     zdudt(:)  = 0.0_wp
     zdvdt(:)  = 0.0_wp
  !
  DO 524 jk=1,klev
     !
     !  WAVE STRESS
     !-------------
     !
!CDIR NODEP
!DIR$ CONCURRENT
     DO 523 jl=1,kgwd
        ji=kdx(jl)
        zdelp=paphm1(ji,jk+1)-paphm1(ji,jk)
        ztemp=-grav*(ztau(ji,jk+1)-ztau(ji,jk))/(zvph(ji,klev+1)*zdelp)
        zdudt(ji)=(pulow(ji)*zd1(ji)-pvlow(ji)*zd2(ji))*ztemp/zdmod(ji)
        zdvdt(ji)=(pvlow(ji)*zd1(ji)+pulow(ji)*zd2(ji))*ztemp/zdmod(ji)
        !
        ! Control Overshoots
        !
        zforc=SQRT(zdudt(ji)**2+zdvdt(ji)**2)
        ztend=SQRT(pum1(ji,jk)**2+pvm1(ji,jk)**2)/ztmst
        rover=0.25_wp
        IF(zforc >= rover*ztend)THEN
           zdudt(ji)=rover*ztend/zforc*zdudt(ji)
           zdvdt(ji)=rover*ztend/zforc*zdvdt(ji)
        ENDIF
        !
        ! BLOCKED FLOW DRAG:
        ! -----------------
        !
        IF(jk > ikenvh(ji)) THEN
           zb=1.0_wp-0.18_wp*pgam(ji)-0.04_wp*pgam(ji)**2
           zc=0.48_wp*pgam(ji)+0.3_wp*pgam(ji)**2
           zconb=2._wp*ztmst*gkwake*psig(ji)/(4._wp*pstd(ji))
           zabsv=SQRT(pum1(ji,jk)**2+pvm1(ji,jk)**2)/2._wp
           zzd1=zb*COS(zpsi(ji,jk))**2+zc*SIN(zpsi(ji,jk))**2
           ratio=(COS(zpsi(ji,jk))**2+pgam(ji)*SIN(zpsi(ji,jk))**2)/      &
                (pgam(ji)*COS(zpsi(ji,jk))**2+SIN(zpsi(ji,jk))**2)
           zbet=MAX(0.0_wp,2.0_wp-1.0_wp/ratio)*zconb*zzdep(ji,jk)*zzd1*zabsv
           !
           ! OPPOSED TO THE WIND
           !
           zdudt(ji)=-pum1(ji,jk)/ztmst
           zdvdt(ji)=-pvm1(ji,jk)/ztmst
           !
           ! PERPENDICULAR TO THE SSO MAIN AXIS:
           !
           !mod     zdudt(ji)=-(pum1(ji,jk)*cos(pthe(ji)*pi/180.)
           !mod *              +pvm1(ji,jk)*sin(pthe(ji)*pi/180.))
           !mod *              *cos(pthe(ji)*pi/180.)/ztmst
           !mod     zdvdt(ji)=-(pum1(ji,jk)*cos(pthe(ji)*pi/180.)
           !mod *              +pvm1(ji,jk)*sin(pthe(ji)*pi/180.))
           !mod *              *sin(pthe(ji)*pi/180.)/ztmst
           !
           zdudt(ji)=zdudt(ji)*(zbet/(1._wp+zbet))
           zdvdt(ji)=zdvdt(ji)*(zbet/(1._wp+zbet))
        END IF
        pvom(ji,jk)=zdudt(ji)
        pvol(ji,jk)=zdvdt(ji)
        zust=pum1(ji,jk)+ztmst*zdudt(ji)
        zvst=pvm1(ji,jk)+ztmst*zdvdt(ji)
        zdis=0.5_wp*(pum1(ji,jk)**2+pvm1(ji,jk)**2-zust**2-zvst**2)
        IF(zdis<0.0_wp) THEN
           zred=SQRT((pum1(ji,jk)**2+pvm1(ji,jk)**2)/(zust**2+zvst**2))
           zust=zust*zred
           zvst=zvst*zred
           pvom(ji,jk)=(zust-pum1(ji,jk))/ztmst
           pvol(ji,jk)=(zvst-pvm1(ji,jk))/ztmst
           zdis=0.5_wp*(pum1(ji,jk)**2+pvm1(ji,jk)**2-zust**2-zvst**2)
        END IF
        pte(ji,jk)=zdis/(ztmst*cpd)
523  END DO
524 END DO

  RETURN
END SUBROUTINE orodrag

SUBROUTINE orosetup                                           &
     ( kproma, kbdim, klev, kgwd, kdx, ktest                  &
     , kkcrit, kkcrith, kcrit                                 &
     , kkenvh, kknu  , kknu2                                  &
     , paphm1, papm1 , pum1   , pvm1 , ptm1  , pgeom1         &
     , prho  , pri   , pstab  , ptau , pvph  ,ppsi, pzdep     &
     , pulow , pvlow                                          &
     , ptheta, pgam, pmea, ppic, pval                         &
     , pnu  ,  pd1  ,  pd2  ,pdmod  )
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
  ! ktest--input-I: Flags to indicate active points
  !
  ! ptsphy--input-R-Time-step (s)
  ! paphm1--input-R: pressure at model 1/2 layer
  ! papm1---input-R: pressure at model layer
  ! pgeom1--input-R: Altitude of layer above ground
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
  INTEGER, INTENT(in) :: kproma, kbdim, klev, kgwd

  INTEGER :: kkcrit(kbdim), kkcrith(kbdim), kcrit(kbdim),    &
       kdx(kbdim), ktest(kbdim),   kkenvh(kbdim)

  REAL(wp) :: pulow(kbdim), pvlow(kbdim), ptheta(kbdim), pgam(kbdim), pnu(kbdim),  &
       pd1(kbdim),   pd2(kbdim),   pdmod(kbdim),                          &
       pmea(kbdim),  ppic(kbdim),  pval(kbdim)

  ! array arguments with intent(IN):
  ! Input 2D
  REAL(wp), INTENT(in) :: paphm1(kbdim,klev+1) ! half level pressure (t-dt)
  REAL(wp), INTENT(in) :: papm1(kbdim,klev)    ! full level pressure (t-dt)
  REAL(wp), INTENT(in) :: pgeom1(kbdim,klev)   ! geopotential above surface (t-dt)
  REAL(wp), INTENT(in) :: ptm1(kbdim,klev)     ! temperature (t-dt)
  REAL(wp), INTENT(in) :: pum1(kbdim,klev)     ! zonal wind (t-dt)
  REAL(wp), INTENT(in) :: pvm1(kbdim,klev)     ! meridional wind (t-dt)

  REAL(wp) ::  prho(kbdim,klev+1), pri(kbdim,klev+1),  pstab(kbdim,klev+1),    &
       ptau(kbdim,klev+1), pvph(kbdim,klev+1), ppsi(kbdim,klev+1),     &
       pzdep(kbdim,klev)
  !
  !  Local scalars:
  INTEGER :: ilevh , jl, jk, ji
  REAL(wp)    :: zcons1, zcons2, zhgeo,  zu,   zphi,   &
       zvt1,   zvt2,   zdwind, zwind, zdelp, &
       zstabm, zstabp, zrhom,  zrhop
  LOGICAL :: lo
  !
  ! local arrays
  INTEGER :: kknu(kbdim),  kknu2(kbdim), kknub(kbdim), kknul(kbdim),    &
       kentp(kbdim), ncount(kbdim)
  REAL(wp) :: znorm(kbdim), zb(kbdim),   zc(kbdim),                        &
       zulow(kbdim), zvlow(kbdim),znup(kbdim), znum(kbdim)
  REAL(wp)    :: zhcrit(kbdim,klev), zvpf(kbdim,klev),   zdp(kbdim,klev)
  LOGICAL :: ll1(kbdim,klev+1)
  !
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
       zhgeo=pgeom1(jl,jk)/grav
       ll1(jl,jk)=(zhgeo > zhcrit(jl,jk))
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
        zhgeo=pgeom1(jl,jk)/grav
        ll1(jl,jk)=(zhgeo > zhcrit(jl,jk))
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
        zhgeo=pgeom1(jl,jk)/grav
        ll1(jl,jk)=(zhgeo > zhcrit(jl,jk))
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
     kknu(jl)=MIN(kknu(jl),nktopg)
     kknu2(jl)=MIN(kknu2(jl),nktopg)
     kknub(jl)=MIN(kknub(jl),nktopg)
     kknul(jl)=klev
2010 END DO
  !
  !     initialize various arrays
  !
  DO 2107 jl=1,kproma
     prho(jl,klev+1)  =0.0_wp
     pstab(jl,1)      =0.0_wp
     pstab(jl,klev+1) =0.0_wp
     pri(jl,1)        =0.0_wp
     pri(jl,klev+1)   =9999.0_wp
     pvph(jl,1)       =0.0_wp
     pvph(jl,klev+1)  =0.0_wp
     ppsi(jl,klev+1)  =0.0_wp
     pulow(jl)        =0.0_wp
     pvlow(jl)        =0.0_wp
     zulow(jl)        =0.0_wp
     zvlow(jl)        =0.0_wp
     kkcrith(jl)      =klev
     kkenvh(jl)       =klev
     kentp(jl)        =klev
     kcrit(jl)        =1
     ncount(jl)       =0
     ll1(jl,klev+1)   =.FALSE.
2107 END DO
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
           pulow(jl)=pulow(jl)+pum1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
           pvlow(jl)=pvlow(jl)+pvm1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
           pstab(jl,klev+1)=pstab(jl,klev+1)                               &
                +pstab(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
           prho(jl,klev+1)=prho(jl,klev+1)                                 &
                +prho(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
        END IF
2116 END DO
2115 END DO

!IBM* ASSERT(NODEPS)
  DO 2110 ji=1,kgwd
     jl = kdx(ji)
     pulow(jl)=pulow(jl)/(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknu2(jl)))
     pvlow(jl)=pvlow(jl)/(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknu2(jl)))
     znorm(jl)=MAX(SQRT(pulow(jl)**2+pvlow(jl)**2),gvsec)
     pvph(jl,klev+1)=znorm(jl)
     pstab(jl,klev+1)=pstab(jl,klev+1)                                 &
          /(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknu2(jl)))
     prho(jl,klev+1)=prho(jl,klev+1)                                   &
          /(paphm1(jl,kknul(jl)+1)-paphm1(jl,kknu2(jl)))
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
     ptau(1:kproma,jk)  =0.0_wp
     pzdep(1:kproma,jk) =0.0_wp
     ppsi(1:kproma,jk)  =0.0_wp
     ll1(1:kproma,jk)   =.FALSE.
213 END DO

!!  DO 215 jk=2,klev-1  ! BUG FIX FOR NaN (undefined values)
  DO 215 jk=2,klev
!DIR$ CONCURRENT
!IBM* ASSERT(NODEPS)
     DO 214 ji=1,kgwd
        jl = kdx(ji)
        zdp(jl,jk)=papm1(jl,jk)-papm1(jl,jk-1)
        pvph(jl,jk)=((paphm1(jl,jk)-papm1(jl,jk-1))*zvpf(jl,jk)+        &
             (papm1(jl,jk)-paphm1(jl,jk))*zvpf(jl,jk-1))           &
             /zdp(jl,jk)
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
       zdwind=MAX(ABS(zvpf(jl,jk)-zvpf(jl,jk-1)),gvsec)
       pri(jl,jk)=pstab(jl,jk)*(zdp(jl,jk)                             &
            /(grav*prho(jl,jk)*zdwind))**2
       pri(jl,jk)=MAX(pri(jl,jk),grcrit)
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
           znum(jl)=pnu(jl)
           zwind=(pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/        &
                MAX(SQRT(pulow(jl)**2+pvlow(jl)**2),gvsec)
           zwind=MAX(SQRT(zwind**2),gvsec)
           zdelp=paphm1(jl,jk+1)-paphm1(jl,jk)
           zstabm=SQRT(MAX(pstab(jl,jk  ),gssec))
           zstabp=SQRT(MAX(pstab(jl,jk+1),gssec))
           zrhom=prho(jl,jk  )
           zrhop=prho(jl,jk+1)
           pnu(jl) = pnu(jl) + (zdelp/grav)*                             &
                ((zstabp/zrhop+zstabm/zrhom)/2._wp)/zwind
           IF((znum(jl) <= gfrcrit).AND.(pnu(jl) > gfrcrit)            &
                .AND.(kkenvh(jl) == klev))              &
                kkenvh(jl)=jk
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
       znum(jl)=znup(jl)
       zwind=(pulow(jl)*pum1(jl,jk)+pvlow(jl)*pvm1(jl,jk))/        &
            MAX(SQRT(pulow(jl)**2+pvlow(jl)**2),gvsec)
       zwind=MAX(SQRT(zwind**2),gvsec)
       zdelp=paphm1(jl,jk+1)-paphm1(jl,jk)
       zstabm=SQRT(MAX(pstab(jl,jk  ),gssec))
       zstabp=SQRT(MAX(pstab(jl,jk+1),gssec))
       zrhom=prho(jl,jk  )
       zrhop=prho(jl,jk+1)
       znup(jl) = znup(jl) + (zdelp/grav)*                           &
            ((zstabp/zrhop+zstabm/zrhom)/2._wp)/zwind
       IF((znum(jl) <= pi/2._wp).AND.(znup(jl) > pi/2._wp)      &
            .AND.(kkcrith(jl) == klev))             &
            kkcrith(jl)=jk
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
           pzdep(jl,jk)=(pgeom1(jl,kkenvh(jl)-1)-pgeom1(jl,jk))/         &
                (pgeom1(jl,kkenvh(jl)  )-pgeom1(jl,klev))
        END IF
253  END DO
254 END DO

  RETURN
END SUBROUTINE orosetup

SUBROUTINE gwstress( kproma, kbdim, klev, kgwd,                     &
     kdx, ktest,  kcrit, kkenvh,  kknu,                             &
     prho,   pstab,  pvph,   pstd,  psig,    ppic,   pval,          &
     ptau,   pdmod )
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

!  USE mo_ssodrag

  IMPLICIT NONE

  ! scalar arguments with intent(IN):
  INTEGER, INTENT(in) :: kproma, kbdim, klev, kgwd

  INTEGER :: kcrit(kbdim), kdx(kbdim), ktest(kbdim),   kkenvh(kbdim), kknu(kbdim)
  !
  REAL(wp) :: psig(kbdim), ppic(kbdim), pval(kbdim), pdmod(kbdim)
  !
  REAL(wp) :: prho(kbdim,klev+1), pstab(kbdim,klev+1), ptau(kbdim,klev+1),  &
       pvph(kbdim,klev+1), ptfr(kbdim), pstd(kbdim)

  !  Local scalars: 
  INTEGER :: jl, ji
  REAL(wp)    :: zeff  ! effective height seen by the flow when there is blocking

  !
  !*         1.1     gravity wave stress.
  !

  ptau(1:kproma,klev+1)=0.0_wp

!IBM* ASSERT(NODEPS)
  DO 301 ji=1,kgwd
     jl = kdx(ji)

     !  effective mountain height above the blocked flow

     zeff=ppic(jl)-pval(jl)
     IF(kkenvh(jl) < klev)THEN
        zeff=MIN(gfrcrit*pvph(jl,klev+1)/SQRT(pstab(jl,klev+1)),zeff)
     ENDIF

     ptau(jl,klev+1)=gkdrag*prho(jl,klev+1)                          &
          *psig(jl)*pdmod(jl)/4._wp/pstd(jl)                         &
          *pvph(jl,klev+1)*SQRT(pstab(jl,klev+1))                    &
          *zeff**2


        !  too small value of stress or  low level flow include critical level
        !  or low level flow:  gravity wave stress nul.

        !       lo=(ptau(jl,klev+1).lt.gtsec).or.(kcrit(jl).ge.kknu(jl))
        !    *      .or.(pvph(jl,klev+1).lt.gvcrit)
        !       if(lo) ptau(jl,klev+1)=0.0_wp

        !      print *,jl,ptau(jl,klev+1)

301 END DO

  !      write(21)(ptau(jl,klev+1),jl=1,kbdim)

  RETURN
END SUBROUTINE gwstress


SUBROUTINE gwprofil( kproma, kbdim, klev,                    &
     kgwd ,kdx  ,                                            &
     kkcrith, kcrit ,                                        &
     paphm1, prho   , pstab , ptfr , pvph , pri , ptau,      &
     pdmod,  psig,    pstd)

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

  ! scalar arguments with intent(IN):
  INTEGER, INTENT(in) :: kproma, kbdim, klev
  INTEGER, INTENT(in) :: kgwd     ! Total points where oro scheme is active    
  !
  INTEGER :: kkcrith(kbdim), kcrit(kbdim), kdx(kbdim)

  REAL(wp) :: pdmod (kbdim), psig(kbdim), pstd(kbdim)
  !
  REAL(wp) :: paphm1(kbdim,klev+1), pstab(kbdim,klev+1),               &
       prho  (kbdim,klev+1), pvph (kbdim,klev+1),               &
       pri   (kbdim,klev+1), ptfr (kbdim), ptau(kbdim,klev+1)

  !  Local scalars: 
  INTEGER :: ji, jl, jk
  REAL(wp) :: zsqr, zalfa, zriw, zdel, zb, zalpha, zdz2n, zdelp, zdelpt

  !  Local arrays: 
  REAL(wp) :: zdz2 (kbdim,klev) , znorm(kbdim) , zoro(kbdim), ztau (kbdim,klev+1)
  !
  !  Executable statements  

!CDIR NODEP
!DIR$ CONCURRENT
  DO 400 ji=1,kgwd
     jl=kdx(ji)
     zoro(jl)=psig(jl)*pdmod(jl)/4._wp/pstd(jl)
     ztau(jl,klev+1)=ptau(jl,klev+1)
     !     print *,jl,ptau(jl,klev+1)
     ztau(jl,kkcrith(jl))=grahilo*ptau(jl,klev+1)
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
              zriw=pri(jl,jk)*(1._wp-zalfa)/(1._wp+zalfa*zsqr)**2
              IF(zriw < grcrit) THEN
                 !                 print *,' breaking!!!',ptau(jl,jk)
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

        IF(jk < ntop)THEN

           zdelp =paphm1(jl,ntop)
           zdelpt=paphm1(jl,jk)
           !         ptau(jl,jk)=ztau(jl,ntop)*zdelpt/zdelp
           ptau(jl,jk)=ztau(jl,ntop)

        ENDIF

533  END DO


531 END DO


  !      do ji=1,kgwd
  !      jl=kdx(ji)
  !      print 123, jl ,(ptau(jl,jk),jk=klev+1,1,-1)
  !      enddo
  !
  ! 123   FORMAT(i4,1x,20(f6.3,1x))


  RETURN
END SUBROUTINE gwprofil

SUBROUTINE orolift( kproma, kbdim, klev,krow,          &
     ktest,                                            &
     paphm1, pgeom1,                                   &
     ptm1,   pum1,  pvm1,                              &
     pmea,   pstd,  ppic,                              &
     pvom,   pvol,  pte )
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
!  USE mo_time_control,      ONLY: delta_time
!  USE mo_geoloc,            ONLY: twomu_2d
!  USE mo_ssodrag
!  USE mo_constants,         ONLY: rd, g, omega

  IMPLICIT NONE
!! KF ATTENETION
  REAL(wp):: delta_time

  ! scalar arguments with intent(IN):
  INTEGER, INTENT(in) :: kproma, kbdim, klev, krow

  ! array arguments with intent(IN):
  ! Input 1D
  INTEGER, INTENT(in) :: ktest(kbdim) ! Flags to indicate active points

  REAL(wp), INTENT(in) :: pmea(kbdim)   ! Mean Orography (m)
  REAL(wp), INTENT(in) :: pstd(kbdim)   ! SSO standard deviation (m)
  REAL(wp), INTENT(in) :: ppic(kbdim)   ! SSO Peacks elevation (m)
  ! Input 2D
  REAL(wp), INTENT(in) :: paphm1(kbdim,klev+1) ! half level pressure (t-dt)
  REAL(wp), INTENT(in) :: pgeom1(kbdim,klev)   ! geopotential above surface (t-dt)
  REAL(wp), INTENT(in) :: ptm1(kbdim,klev)     ! temperature (t-dt)
  REAL(wp), INTENT(in) :: pum1(kbdim,klev)     ! zonal wind (t-dt)
  REAL(wp), INTENT(in) :: pvm1(kbdim,klev)     ! meridional wind (t-dt)


  ! array arguments with intent(INOUT):
  ! Input 2D
  REAL(wp), INTENT(inout) :: pte(kbdim,klev) ! tendency of temperature
  REAL(wp), INTENT(inout) :: pvol(kbdim,klev) ! tendency of meridional wind
  REAL(wp), INTENT(inout) :: pvom(kbdim,klev) ! tendency of zonal wind

  ! Local scalars:
  INTEGER :: jl, ilevh, jk
  REAL(wp)    :: zhgeo, zdelp, zslow, zsqua, zscav, zbet
  REAL(wp)    :: zcons1,ztmst
  LOGICAL :: lifthigh

  ! Local arrays:   
  INTEGER :: iknub(kbdim),  iknul(kbdim)
  REAL(wp)    ::  zdudt(kbdim), zdvdt(kbdim)
 
  REAL(wp)    :: plat(kbdim)                ! Latitude (radiants)
  REAL(wp)    :: pulow(kbdim), pvlow(kbdim)
  REAL(wp)    :: ztau(kbdim,klev+1), ztav(kbdim,klev+1), zrho(kbdim,klev+1)
  REAL(wp)    :: zhcrit(kbdim,klev)
  LOGICAL :: ll1(kbdim,klev+1)

  !-----------------------------------------------------------------------
  !
  !*         1.1  initialisations
  !               ---------------

  lifthigh=.FALSE.

  zcons1=1._wp/rd
  ztmst=2._wp*delta_time

!! KF ATTENTION
  DO jl=1,kproma
    plat(jl)= 0._wp !ASIN(0.5_wp*twomu_2d(jl,krow))
  END DO
  !
  DO 1001 jl=1,kproma
     zrho(jl,klev+1)  =0.0_wp
     pulow(jl)        =0.0_wp
     pvlow(jl)        =0.0_wp
     iknub(jl)   =klev
     iknul(jl)   =klev
     ilevh=klev/3
     ll1(jl,klev+1)=.FALSE.
     DO 1000 jk=1,klev
        pvom(jl,jk)=0.0_wp
        pvol(jl,jk)=0.0_wp
        pte (jl,jk)=0.0_wp
1000 END DO
1001 END DO

  !
  !*         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
  !*                 LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
  !*                 THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.
  !

  DO 2006 jk=klev,1,-1
     DO 2007 jl=1,kproma
        IF(ktest(jl) == 1) THEN
           zhcrit(jl,jk)=MAX(ppic(jl)-pmea(jl),100.0_wp)
           zhgeo=pgeom1(jl,jk)/grav
           ll1(jl,jk)=(zhgeo > zhcrit(jl,jk))
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
        IF(iknub(jl) > nktopg) iknub(jl)=nktopg
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
              pulow(jl)=pulow(jl)+pum1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
              pvlow(jl)=pvlow(jl)+pvm1(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
              zrho(jl,klev+1)=zrho(jl,klev+1)                          &
                   +zrho(jl,jk)*(paphm1(jl,jk+1)-paphm1(jl,jk))
           ENDIF
        ENDIF
2116 END DO
2115 END DO

  DO 2110 jl=1,kproma
     IF(ktest(jl) == 1) THEN
        pulow(jl)=pulow(jl)/(paphm1(jl,iknul(jl)+1)-paphm1(jl,iknub(jl)))
        pvlow(jl)=pvlow(jl)/(paphm1(jl,iknul(jl)+1)-paphm1(jl,iknub(jl)))
        zrho(jl,klev+1)=zrho(jl,klev+1)                                &
             /(paphm1(jl,iknul(jl)+1)-paphm1(jl,iknub(jl)))
     ENDIF
2110 END DO



  !***********************************************************
  !
  !*         3.      COMPUTE MOUNTAIN LIFT
  !
  DO 301 jl=1,kproma
     IF(ktest(jl) == 1) THEN
        ztau(jl,klev+1)= - gklift*zrho(jl,klev+1)*2._wp*omega*        &
             !                (2*pstd(jl)+pmea(jl))*   &
2       *pstd(jl)*                                      &
             SIN(plat(jl))*pvlow(jl)
        ztav(jl,klev+1)=   gklift*zrho(jl,klev+1)*2._wp*omega*        &
             !                (2*pstd(jl)+pmea(jl))*   &
2       *pstd(jl)*                                       &
             SIN(plat(jl))*pulow(jl)
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
              zdelp=paphm1(jl,jk+1)-paphm1(jl,jk)
              zdudt(jl)=-grav*(ztau(jl,jk+1)-ztau(jl,jk))/zdelp
              zdvdt(jl)=-grav*(ztav(jl,jk+1)-ztav(jl,jk))/zdelp
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
              zbet=gklift*2._wp*omega*SIN(plat(jl))*ztmst*          &
                   (pgeom1(jl,iknub(jl)-1)-pgeom1(jl,  jk))/      &
                   (pgeom1(jl,iknub(jl)-1)-pgeom1(jl,klev))
              zdudt(jl)=-pum1(jl,jk)/ztmst/(1._wp+zbet**2)
              zdvdt(jl)=-pvm1(jl,jk)/ztmst/(1._wp+zbet**2)
              pvom(jl,jk)= zbet**2*zdudt(jl) - zbet   *zdvdt(jl)
              pvol(jl,jk)= zbet   *zdudt(jl) + zbet**2*zdvdt(jl)
           ENDDO
        ENDIF
601  END DO

  ENDIF
  !           PRINT *,' out orolift'
  RETURN 
END SUBROUTINE orolift

END MODULE mo_ssortns
