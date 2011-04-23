!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
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
@PROCESS XLF90(NOSIGNEDZERO)
#else
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#define SWDIV_NOCHK(a,b) ((a)/(b))
#endif

MODULE mo_cover

  USE mo_kind,               ONLY : wp
  USE mo_physical_constants, ONLY : rd, cpd, vtmpc1
#ifdef __ibmspline__
  USE mo_convect_tables,     ONLY : prepare_ua_index_spline, lookup_ua_eor_uaw_spline
#else
  USE mo_convect_tables,     ONLY : prepare_ua_index, lookup_ua_eor_uaw
#endif
  USE mo_echam_cloud_params, ONLY : jbmin1, ncctop, cqtmin, cbeta_pq      &
                                  , jbmin, jbmax                          &
                                  , csatsc, ccwmin, cbeta_pq_max, nbetaq  &
                                  , cbetaqs, rbetak, nbetax, tbetai       &
                                  , cvarmin, cmmrmax, crt, crs            &
                                  , nex
#ifdef _PROFILE
  USE mo_profile,         ONLY : trace_start, trace_stop
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cover

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS

SUBROUTINE cover (         kproma,   kbdim, ktdia, klev, klevp1, lcover  &
!
! - INPUT  1D .
                         , ktype,    pfrw                                &
! - INPUT  2D .
                         , paphm1,   papm1                               &
                         , pvervel,  ptm1                                &
                         , pqm1,     pxlm1,    pxim1                     &
! - INPUT AND OUTPUT, 2D
                         , paclc,    pxvar,    pxskew                    &
! - OUTPUT 1D .
                         , knvb,     printop                             &
! - OUTPUT 2D .
                         , pbetaa,   pbetab,   pbetass                   &
                 )
!-----------------------------------------------------------------------
!     *Cover*
!          Diagnoses cloud cover for current timestep
!
!     Subject.
!     --------
!     In ECHAM4 all cloud cover calculations were performed
!     in routine CLOUD.  The routine diagnosed the cloud cover
!     from the previous timestep using m1 variables, but then
!     used a "first guess" calculation of T and Q to give a
!     cloud cover estimate for the next timestep.  This meant
!     that the radiation scheme used cloud cover values that were
!     from a different timestep to the temperature and water vapour
!     values, and also that the T and Q values were anyway
!     preliminary.  Finally, the cover calculation was performed
!     twice each timestep, when one calculation suffices.
!
!     This scheme calculates cover diagnostically and is called
!     once at the beginning of each timestep.  It uses the
!     standard relative humidity calculation from Lohmann and
!     Roeckner (96), or the method from the new prognostic
!     scheme of Tompkins.  The choice of which scheme to use is
!     controlled by the parameter switch LCOVER, which is set in
!     namelist PHYSCTL along with lsurf etc... Note that even if
!     lcover=.false. (RH scheme) you can't restart this model version
!     from restart files saved from a different model version, since
!     the two extra prognostic equations, pxvar and pxskew are still
!     stored even though they are not actively used.  However, this means
!     that once you have restart files from this version, you are able
!     change lcover at will.
!
!     In the new scheme the variable xskew is provided
!     as outlined in the reference, this variable represents
!     directly the Beta distribution shape parameter "q"
!     The shape parameter "p" (zbetap) a tunable parameter and it is
!     recommended that this be set to a low constant 1.5<p<2.0
!     (This may be changed later to prognostic to allow negative skewness
!     from downdraft detrainment, see ref. for details).
!
!     from xi,xl,q,xskew and zbetap, the Beta distribution is defined
!     and cloud cover is diagnosable.  For the iteration, Ridders' method
!     is used (see Numerical Recipes).
!
!     Attention:
!     In the current version the advective tendencies of skewness
!     and variance are set to zero.
!
!     INTERFACE.
!     ----------
!
!     *Call cover*
!
!     Input arguments.
!     ----- ----------
!  - 1D
!  ktype    : type of convection
!  pfrw     :
!
!  - 2D
!  paphm1   : pressure at half levels                              (n-1)
!  papm1    : pressure at full levels                              (n-1)
!  pvervel  : vertical velocity in pressure coordiante             (n)
!  ptm1     : temperature                                          (n-1)
!  pqm1     : specific humidity                                    (n-1)
!  pxlm1    : cloud liquid water                                   (n-1)
!  pxim1    : cloud ice                                            (n-1)
!
!     Input and output arguments.
!     ---------------------------
!  - 2D
!  paclc    : cloud cover
!  pxvar    : the beta distribution width "b-a"                    (n-1)
!  pxskew   : the beta distribution shape parameter "q"            (n-1)
!
!     Output arguments.
!     ------ ----------
!  - 1D
!  knvb     :
!  pfrw     :
!
!  - 2D
!  pbetaa   : the beta distribution minimum a                      (n-1)
!  pbetab   : the beta distribution maximum b                      (n-1)
!  pbetass  :
!
!     Method.
!     -------
!     see References
!
!     References.
!     ----------
!     Diagnostic CC scheme: Lohmann and Roeckner 96, Clim. Dyn.
!     Prognostic CC scheme: Tompkins 2002, J. Atmos. Sci.
!
!     Authors.
!     -------
!     A. Tompkins    MPI-Hamburg  2000
!     K. Ketelesen   NEC, April 2002
!     L. Kornblueh   MPI, April 2002
!
!     Modifications.
!     --------------
!     view cover
!       v2: first working version
!       v5: lookup table added
!       v8: zriddr and functions replaced for vectorization
!       v9: optimizations, longer vector loop and less indirect addressing
!          - introduction of additional arrays
!          - change structure if "beta function scheme" IF block
!          - scattered loops "ictit" and "ictdg" are collected over "kproma" and "klev"
!          - intoducing additional arrays to hold data in "ictit" and "ictdg"
!            addressing scheme
!
  INTEGER, INTENT(IN)    :: kproma, kbdim, ktdia, klev, klevp1
  LOGICAL, INTENT(IN)    :: lcover
  INTEGER, INTENT(IN)    ::  ktype(kbdim)
  REAL(wp),INTENT(IN)    ::   pfrw(kbdim)
  REAL(wp),INTENT(IN)    :: paphm1(kbdim,klevp1),   papm1(kbdim,klev)
  REAL(wp),INTENT(IN)    ::   ptm1(kbdim,klev),   pvervel(kbdim,klev)
  REAL(wp),INTENT(IN)    ::   pqm1(kbdim,klev)
  REAL(wp),INTENT(IN)    ::  pxlm1(kbdim,klev),     pxim1(kbdim,klev)

  REAL(wp),INTENT(INOUT) ::  paclc(kbdim,klev)
  REAL(wp),INTENT(INOUT) ::  pxvar(kbdim,klev),    pxskew(kbdim,klev)

  INTEGER, INTENT(OUT)   ::  knvb(kbdim)
  REAL(wp),INTENT(OUT)   ::  printop(kbdim)
  REAL(wp),INTENT(OUT)   ::  pbetaa(kbdim,klev),   pbetab(kbdim,klev)
  REAL(wp),INTENT(OUT)   ::  pbetass(kbdim,klev)
!
! Local variables -----------------------------------------------
!
  INTEGER  :: jl, jk, kbeta, jb, jt
  INTEGER  :: locnt, nl
  REAL(wp) :: rdcpd, zdthdp, zcor, zbetai0, zbetai1
  REAL(wp) :: zrhc, zsat, zqr

  INTEGER  ::  itv1(kproma*klev), itv2(kproma*klev)
  REAL(wp) ::  zdthmin(kbdim),    ztheta(kbdim,klev)

#ifdef __ibmspline__
  REAL(wp) ::  za(kbdim)
#endif

!
!   Pointers and counters for iteration and diagnostic loop:
!
  INTEGER :: iptit(2,kproma*klev), iptdg(2,kproma*klev)                &
           , iqidx_l1(kproma*klev)
  INTEGER :: ictdg, ictit, ix, ix1, ix2, iq, nphase
  INTEGER :: iqidx(kproma,klev)
!
!   variables required for zriddr iteration scheme:
!
  REAL(wp) :: fh(kproma*klev), fl(kproma*klev), fm(kproma*klev)        &
            , fnew(kproma*klev), xh(kproma*klev), xl(kproma*klev)      &
            , xm(kproma*klev), xnew(kproma*klev)

  REAL(wp) :: unused = -1.11e30_wp
  REAL(wp) :: ztt, zx1, zx2, zvar, zvartarget, zss, zqt, zpp, zqq      &
            , zaa, zbb, zjk, zlo2, zpq, zpqi, zppi             &
            , zqqi, zaa2, zbb2, paclc1
  REAL(wp) :: zqsm1(kproma*klev), zskew1(kproma*klev)                  &
            , zskew2(kproma*klev), zskew(kbdim)
  REAL(wp) :: zbetaqt(kproma,klev), zbetacl(kproma,klev), ua(kproma)

  LOGICAL :: lo1, lo2, lao, lao1
  INTEGER :: liter(kproma*klev)

  INTEGER :: icond(kproma)
  REAL(wp) :: pqm1_l1(kproma*klev), pqm1_l2(kproma*klev), zpapm1i(kbdim,klev) &
            , ztmp(kproma*klev), zalpha1(kproma*klev), zalpha2(kproma*klev)
  REAL(wp) :: zbetaqt_l1(kproma*klev), zbetacl_l1(kproma*klev)         &
            , zbetass_l1(kproma*klev)
  REAL(wp) :: zbetaqt_l2(kproma*klev), zbetacl_l2(kproma*klev)         &
            , zbetass_l2(kproma*klev)

  REAL(wp) :: zknvb(kbdim),zphase(kbdim)

  INTEGER :: lo(kbdim,klev),loidx(kproma*klev)

!   number of iteration loops.  In a small test, most gridpoints
!   converged in 2 to 4 iterations, with occassionally as many as
!   8 needed for an accuracy of 0.01*(ql+qi). Don't set too high since
!   the loops are all performed, albeit often as "empty loops" (i.e.
!   LITER is false because of convergence), due to vectorisation needs.
!
  INTEGER, PARAMETER :: niter=60            ! max loop number
  REAL(wp), PARAMETER :: ziter_acc=0.01_wp  ! iteration accuracy*(ql+qi)
!
#ifdef _PROFILE
  CALL trace_start ('cover', 9)
#endif
!
!   Computational constants
!
  rdcpd  = rd/cpd
!
!   Initialize variables
!
  DO jl=1,kproma
    zdthmin(jl)=0.0_wp
    zknvb(jl)  =1.0_wp
    printop(jl)=1.0_wp
  END DO
!
  DO jk = ktdia,klev
     DO jl = 1,kproma
        zpapm1i(jl,jk) = SWDIV_NOCHK(1._wp,papm1(jl,jk))
        ztheta(jl,jk) = ptm1(jl,jk)*(1.0e5_wp*zpapm1i(jl,jk))**rdcpd
     END DO
  END DO
!
  IF (lcover) THEN
     kbeta = MAX(ncctop,ktdia)
  ELSE
     kbeta = klev+1
  END IF ! lcover
!
!       1.3   Checking occurrence of low-level inversion
!             (below 1000 m, sea points only, no convection)
!

  locnt = 0
  DO jl = 1,kproma
     IF (pfrw(jl).GT.0._wp.AND.ktype(jl).EQ.0) THEN
        locnt = locnt + 1
        loidx(locnt) = jl
     END IF
  END DO

  IF (locnt.GT.0) THEN
     DO jk = klev,jbmin1,-1

!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* novector
        DO nl = 1,locnt
           jl = loidx(nl)
           ztmp(nl) = (ztheta(jl,jk)-ztheta(jl,jk-1)) / (papm1(jl,jk)-papm1(jl,jk-1))
        END DO

        zjk = REAL(jk,wp)
!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO nl = 1,locnt
           jl = loidx(nl)
           zdthdp      = ztmp(nl)
           zknvb(jl)   = FSEL(zdthdp-zdthmin(jl),zknvb(jl),zjk)
           !        zknvb(jl)    = MERGE(zjk,zknvb(jl),zdthdp .LT. zdthmin(jl))
           zdthmin(jl) = MIN(zdthdp,zdthmin(jl))
        END DO
     END DO
  END IF
  knvb(1:kproma) = INT(zknvb(1:kproma))


!
!   Tunable parameters now in mo_cloud.f90
!
  IF (ncctop > 1) THEN
    DO jk=1,ncctop-1
      DO jl=1,kproma
        pbetaa(jl,jk)=0.0_wp
        pbetab(jl,jk)=0.0_wp
        pxvar(jl,jk)=cqtmin
        pxskew(jl,jk)=cbeta_pq
      END DO
    END DO
  END IF
!
  ictit=0
  ictdg=0
  DO jk=ktdia,klev
    IF (lcover .AND. jk >= ncctop) THEN  ! beta function scheme
!
!       1.   Calculate the saturation mixing ratio
!
#ifdef __ibmspline__
        CALL prepare_ua_index_spline('cover (1)',kproma,ptm1(1,jk),itv1(1),  &
                                    & za(1),pxim1(1,jk),nphase,zphase,itv2(1))
        CALL lookup_ua_eor_uaw_spline(kproma,itv1(1),za(1),nphase,itv2(1),ua(1))
#else
        CALL prepare_ua_index('cover (1)',kproma,ptm1(1,jk),itv1(1),pxim1(1,jk), &
                             & nphase,zphase,itv2(1))
        CALL lookup_ua_eor_uaw(kproma,itv1(1),nphase,itv2(1),ua(1))
#endif

!IBM* NOVECTOR
        DO jl=1,kproma
           zqsm1(jl) = MIN(ua(jl)*zpapm1i(jl,jk),0.5_wp)
           zcor      = 1._wp/(1._wp-vtmpc1*zqsm1(jl))
           zqsm1(jl) = zqsm1(jl)*zcor
           icond(jl) = INT(FSEL(-pvervel(jl,jk),0._wp,1._wp))
        END DO

        IF (kproma > 0) THEN
           DO jl=1,kproma
              jb=knvb(jl)
              lo2=(jb.GE.jbmin .AND. jb.LE.jbmax .AND. icond(jl).GT.0)
              lo1=(jk.EQ.jb .OR. jk.EQ.jb+1)
              IF (lo2.AND.lo1) THEN
                 zqsm1(jl)=zqsm1(jl)*csatsc
                 printop(jl)=REAL(jb,wp)
              END IF
           ENDDO
        END IF
        !
        !       2.    calculate cloud cover
        !
        !   Don't need to iterate at every gridpoint, thus make pointers
        !
        DO jl=1,kproma
           zskew(jl) = MAX(MIN(pxskew(jl,jk),cbeta_pq_max),cbeta_pq)
           ztmp(jl)  = (zskew(jl)-cbeta_pq)/rbetak+1._wp
        END DO

        ztmp(1:kproma) = LOG(ztmp(1:kproma))
        iqidx(1:kproma,jk) = INT((REAL(nbetaq,wp)/cbetaqs)*ztmp(1:kproma)+0.5_wp)

        IF (kproma > 0) THEN
           DO jl=1,kproma
              zbetacl(jl,jk)=MAX(0._wp,pxlm1(jl,jk))+MAX(0._wp,pxim1(jl,jk))

              zbetaqt(jl,jk)=MAX(cqtmin,pqm1(jl,jk))+zbetacl(jl,jk)
              pbetass(jl,jk)=MAX(pqm1(jl,jk),zqsm1(jl)) !safety

              ! lo=(pxim1(jl,jk)>ccwmin .OR. pxlm1(jl,jk)>ccwmin) .AND. pqm1(jl,jk)<pbetass(jl,jk)

              zlo2 = FSEL(ccwmin-pxim1(jl,jk),0._wp,1._wp)
              zlo2 = FSEL(ccwmin-pxlm1(jl,jk),zlo2,1._wp)
              zlo2 = FSEL(pqm1(jl,jk)-pbetass(jl,jk),0._wp,zlo2)
              lo(jl,1) = INT(zlo2)
           END DO
        END IF

        IF (ktdia > 0) THEN
           DO jl=1,kproma

              IF (lo(jl,1).EQ.0) THEN
                 ! mpuetz: this is the more probable path
                 ictdg=ictdg+1
                 iptdg(1,ictdg)=jl
                 iptdg(2,ictdg)=jk
                 pqm1_l2(ictdg)    = pqm1(jl,jk)
                 zbetaqt_l2(ictdg) = zbetaqt(jl,jk)
                 zbetacl_l2(ictdg) = zbetacl(jl,jk)
                 zbetass_l2(ictdg) = pbetass(jl,jk)
                 zskew2(ictdg)     = zskew(jl)
              ELSE
                 ictit=ictit+1
                 iptit(1,ictit)=jl
                 iptit(2,ictit)=jk
                 pqm1_l1(ictit)    = pqm1(jl,jk)
                 zbetaqt_l1(ictit) = zbetaqt(jl,jk)
                 zbetacl_l1(ictit) = zbetacl(jl,jk)
                 zbetass_l1(ictit) = pbetass(jl,jk)
                 zskew1(ictit)     = zskew(jl)
                 iqidx_l1(ictit)   = iqidx(jl,jk)
              ENDIF
              !
!
!   setup index for skewness, q, in lookup table (doesn't change)
!
           ENDDO
        END IF

    ENDIF !lcover
  ENDDO
!
  IF (lcover) THEN
!
!   Partially cloudy: Iterative gridpoints
!   uses func using ridders' method, return the root of a function
!   Method: See Numerical Recipes
!
!
!IBM* novector
  DO jl=1,ictit
     !
     !   Lower bound < 0 to catch occassional overflow
     !
     zx1 =-0.1_wp
     xl(jl)=zx1

     ztt=(zbetass_l1(jl)-zx1)*cbeta_pq /                     &
          ((zbetaqt_l1(jl)-zx1)*(cbeta_pq+zskew1(jl)))
     ztt=REAL(nbetax,wp)*MAX(MIN(ztt,1.0_wp),0.0_wp)
     zalpha1(jl) = ztt - AINT(ztt,wp)
     itv1(jl)    = INT(ztt)
  END DO

!IBM* novector
  DO jl=1,ictit

     zx2 =(cbeta_pq+zskew1(jl))*zbetaqt_l1(jl) - cbeta_pq*zbetass_l1(jl)
     zx2 = zx2 / zskew1(jl)
     zx2 = MIN(zx2,zbetaqt_l1(jl),zbetass_l1(jl))
     xh(jl)=zx2
  END DO

!IBM* novector
  DO jl=1,ictit
     zx2 = xh(jl)
     ztt = (zbetass_l1(jl)-zx2)*cbeta_pq
     ztt = ztt / ((zbetaqt_l1(jl)-zx2)*(cbeta_pq+zskew1(jl)))
     ztt = REAL(nbetax,wp)*MAX(MIN(ztt,1.0_wp),0.0_wp)
     zalpha2(jl) = ztt - DINT(ztt)
     itv2(jl)    = INT(ztt)

     xnew(jl)=unused

     !    if (xl(jl)==xh(jl)) liter(jl)=.false.
     liter(jl) = INT(FSEL(-ABS(xl(jl)-xh(jl)),0._wp,1._wp))

  END DO

  IF (ictit > 0) THEN ! mpuetz : don't fuse with previous loop
!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO jl=1,ictit
        ix1 = itv1(jl)
        ix2 = itv2(jl)
        zx1 = xl(jl)
        zx2 = xh(jl)

        iq = iqidx_l1(jl)

        zbetai0=zalpha1(jl)*tbetai(0,iq,ix1+1) + (1._wp - zalpha1(jl))*tbetai(0,iq,ix1)
        zbetai1=zalpha1(jl)*tbetai(1,iq,ix1+1) + (1._wp - zalpha1(jl))*tbetai(1,iq,ix1)
        fl(jl)=(zbetaqt_l1(jl)-zx1)*zbetai1 -                              &
             (zbetass_l1(jl)-zx1)*zbetai0 +                                &
             zbetass_l1(jl)-MAX(cqtmin,pqm1_l1(jl))
        !
        !   3 conditions for the maximum iteration bound: a<qs,a<qt,b>qs
        !
        zbetai0=zalpha2(jl)*tbetai(0,iq,ix2+1) + (1._wp - zalpha2(jl))*tbetai(0,iq,ix2)
        zbetai1=zalpha2(jl)*tbetai(1,iq,ix2+1) + (1._wp - zalpha2(jl))*tbetai(1,iq,ix2)
        fh(jl)=(zbetaqt_l1(jl)-zx2)*zbetai1 -                              &
             (zbetass_l1(jl)-zx2)*zbetai0 +                                &
             zbetass_l1(jl)-MAX(cqtmin,pqm1_l1(jl))
        !    PRINT '(A6,3 I7,2 E26.18,I2)','liter',jl,itv1(jl),itv2(jl),fl(jl),fh(jl)
     END DO
  END IF

!  IF (ictit > 0) STOP 'debug'

!IBM* NOVECTOR
  DO jt=1,niter   ! short iteration loop

     locnt = 1
     DO jl=1,ictit
        loidx(locnt) = jl
        locnt = locnt + liter(jl)
     END DO
     locnt = locnt - 1

     IF (locnt.EQ.0) exit   ! all sites converged

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO nl=1,locnt   ! longer loop over indices
        jl = loidx(nl)
        xm(jl)=0.5_wp*(xl(jl)+xh(jl))
        ztt=(zbetass_l1(jl)-xm(jl))*cbeta_pq
        ztt=SWDIV_NOCHK(ztt,((zbetaqt_l1(jl)-xm(jl))*(cbeta_pq+zskew1(jl))))
        ztt=REAL(nbetax,wp)*MAX(MIN(ztt,1.0_wp),0.0_wp)
        zalpha1(nl) = ztt - DINT(ztt)
        itv1(nl)    = INT(ztt)
     END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO nl=1,locnt   ! longer loop over indices
        jl = loidx(nl)
        iq = iqidx_l1(jl)
        ix = itv1(nl)

        zbetai0=zalpha1(nl)*tbetai(0,iq,ix+1) + (1._wp - zalpha1(nl))*tbetai(0,iq,ix)
        zbetai1=zalpha1(nl)*tbetai(1,iq,ix+1) + (1._wp - zalpha1(nl))*tbetai(1,iq,ix)

        fm(jl)=(zbetaqt_l1(jl)-xm(jl))*zbetai1-                        &
             (zbetass_l1(jl)-xm(jl))*zbetai0+                          &
             zbetass_l1(jl)-MAX(cqtmin,pqm1_l1(jl))
        ztmp(nl)=MAX(fm(jl)**2-fl(jl)*fh(jl),0._wp)
     END DO

     ztmp(1:locnt) = SQRT(ztmp(1:locnt))

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO nl=1,locnt   ! longer loop over indices
        jl = loidx(nl)

        ztt = ztmp(nl)
        ztt = xm(jl)+(xm(jl)-xl(jl))                                   &
            * (sign(1._wp,fl(jl)-fh(jl))*SWDIV_NOCHK(fm(jl),ztt)) !update formula
        ztt = FSEL(-ztmp(nl),xnew(jl),ztt)

        liter(jl) = INT(FSEL(ziter_acc*zbetacl_l1(jl)-ABS(xnew(jl)-ztt),0._wp,1._wp))
!        IF (abs(xnew(jl)-ztt).le.ziter_acc*zbetacl_l1(jl)) liter(jl)=.false.
        xnew(jl)=ztt
        ztt=(zbetass_l1(jl)-xnew(jl))*cbeta_pq/                        &
             ((zbetaqt_l1(jl)-xnew(jl))*                               &
             (cbeta_pq+zskew1(jl)))
        ztt=REAL(nbetax,wp)*MAX(MIN(ztt,1.0_wp),0.0_wp)
        zalpha1(nl) = ztt - DINT(ztt)
        itv1(nl) = INT(ztt)
     END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO nl=1,locnt   ! longer loop over indices
        jl = loidx(nl)
        iq = iqidx_l1(jl)
        ix = itv1(nl)

        zbetai0=zalpha1(nl)*tbetai(0,iq,ix+1) + (1._wp - zalpha1(nl))*tbetai(0,iq,ix)
        zbetai1=zalpha1(nl)*tbetai(1,iq,ix+1) + (1._wp - zalpha1(nl))*tbetai(1,iq,ix)

        fnew(jl)=(zbetaqt_l1(jl)-xnew(jl))*zbetai1-                    &
             (zbetass_l1(jl)-xnew(jl))*zbetai0+                        &
             zbetass_l1(jl)-MAX(cqtmin,pqm1_l1(jl))
     END DO

!CDIR NODEP
!IBM* ASSERT(NODEPS)
     DO nl=1,locnt   ! longer loop over indices
        jl = loidx(nl)
!
!   bookkeeping to keep the root bracketed on next iteration.
!
        IF (sign(fm(jl),fnew(jl)).ne.fm(jl)) THEN
          xl(jl)=xm(jl)
          fl(jl)=fm(jl)
          xh(jl)=xnew(jl)
          fh(jl)=fnew(jl)
        ELSE IF(sign(fl(jl),fnew(jl)).ne.fl(jl)) THEN
          xh(jl)=xnew(jl)
          fh(jl)=fnew(jl)
        ELSE IF (sign(fh(jl),fnew(jl)).ne.fh(jl)) THEN
          xl(jl)=xnew(jl)
          fl(jl)=fnew(jl)
        ENDIF
    END DO

  ENDDO !niter
!
!   Set a and diagnose b
!
!DIR$ CONCURRENT
!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* novector
  DO nl=1,ictit
     jl = iptit(1,nl)
     jk = iptit(2,nl)
     zvartarget=MAX(cqtmin,cvarmin*pqm1_l1(nl))
     zss=zbetass_l1(nl)
     zqt=zbetaqt_l1(nl)
     zpp=cbeta_pq
     zqq=zskew1(nl)
     zaa=MAX(xnew(nl),0.0_wp)
     zaa=MIN(zaa,zqt-zvartarget*(zpp/(zpp+zqq)))
     zbb=(zqt-zaa)*((zpp+zqq)/zpp)+zaa
     zbb=MAX(zbb,zaa+zvartarget)
     pbetaa(jl,jk)=zaa
     pbetab(jl,jk)=zbb
 ENDDO
!
!   Overcast or clear sky: Diagnostic gridpoints
!
!DIR$ CONCURRENT
!CDIR NODEP
!IBM* ASSERT(NODEPS)
!IBM* novector
  DO nl=1,ictdg
     jl = iptdg(1,nl)
     jk = iptdg(2,nl)
     zvartarget=MAX(cqtmin,cvarmin*pqm1_l2(nl))
     zvar=MAX(pxvar(jl,jk),zvartarget) ! b-a width
!
!   Defined to make equations easier to understand:
!
    zss=zbetass_l2(nl)  ! qsat
    zqt=zbetaqt_l2(nl)  ! qtotal
    zpp=cbeta_pq        ! p shape factor
    zqq=zskew2(nl)      ! q shape factor
!
!   Limit a>0 and ensure (a,b)<qs or (a,b)>qs
!
    zpq  = (zpp+zqq)
    zppi = 1._wp/zpp
    zqqi = 1._wp/zqq
    zpqi = 1._wp/zpq

    zaa = MAX(0._wp,zqt-zvar*zpp*zpqi)    !a from eqn(12) limited
    zbb = (zqt-zaa)*zpq*zppi+zaa          !b from eqn(12) again

    zaa2 = MAX(0.0_wp,(zqt*zpq-zss*zpp)*zqqi)
    zbb2 = (zqt-zss)*zpq*zppi+zss          !b from eqn(12) again

    zaa2 = FSEL(zqt - zss,zss,zaa2)
    zbb2 = FSEL(zqt - zss,zbb2,zss)

    zlo2 = FSEL(zss - zbb,0._wp,1._wp)
    zlo2 = FSEL(zaa - zss,0._wp,zlo2)

    zaa = FSEL(-zlo2,zaa,zaa2)
    zbb = FSEL(-zlo2,zbb,zbb2)

    zbb=MAX(zbb,zaa+zvartarget)
    pbetaa(jl,jk)=zaa
    pbetab(jl,jk)=zbb
  ENDDO
!
  IF (kbeta <= klev) THEN

     ! beta function scheme

     nl = 0
     DO jk=kbeta,klev
        DO jl=1,kproma
!
!   Define variance
!
        pxvar(jl,jk) = pbetab(jl,jk)-pbetaa(jl,jk)
!
!       set cloud fraction
!
           paclc(jl,jk) = FSEL(pqm1(jl,jk) - pbetass(jl,jk),1.0_wp,2.0_wp)
           paclc1       = FSEL(ccwmin-pxim1(jl,jk),0.0_wp,paclc(jl,jk))
           paclc(jl,jk) = FSEL(ccwmin-pxlm1(jl,jk),paclc1,paclc(jl,jk))

           ! create an integer flag for all cases where paclc = 2

           lo(jl,jk)       = INT(0.5_wp*paclc(jl,jk))
        END DO
     END DO

     ! compute an index for all (jl,jk) where paclc = 2

     locnt = 1
     IF (kproma > 0) THEN
     DO jk=kbeta,klev
        DO jl=1,kproma
           iptit(1,locnt) = jl
           iptit(2,locnt) = jk
           locnt = locnt + lo(jl,jk)
        END DO
     END DO
     END IF
     locnt = locnt - 1

     IF (locnt.GT.0) THEN

        ! compute paclc using beta function scheme for all (jl,jk) where  0 < pcalc < 1

!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO nl=1,locnt
           jl = iptit(1,nl)
           jk = iptit(2,nl)

           ztt=SWDIV_NOCHK((pbetass(jl,jk)-pbetaa(jl,jk)),(pbetab(jl,jk)-pbetaa(jl,jk)))
           ztt=REAL(nbetax,wp)*MAX(MIN(ztt,1.0_wp),0.0_wp)
           zalpha1(nl) = ztt - DINT(ztt)

!#define __power6opt__
#ifdef __power6opt__
           loidx(nl)   = INT(ztt)
        END DO

        ! power6 optimization : split loop to hide float-to-int coversion load-hit-store latency

!CDIR NODEP
!IBM* ASSERT(NODEPS)
        DO nl=1,locnt
           jl = iptit(1,nl)
           jk = iptit(2,nl)
           ix = loidx(nl)
#else
           ix = INT(ztt)
#endif
           iq = iqidx(jl,jk)

           zbetai0 = zalpha1(nl)*tbetai(0,iq,ix+1) +(1._wp - zalpha1(nl))*tbetai(0,iq,ix)
           paclc(jl,jk) = MAX(zbetacl(jl,jk)/cmmrmax,1.0_wp - zbetai0)
           !    Fractional cloud cover > 0.01 required for middle atmosphere
           paclc(jl,jk) = MAX(MIN(paclc(jl,jk),1.0_wp),0.01_wp)
        END DO

     END IF
  END IF
 END IF   ! lcover
!
!       1.   Calculate the saturation mixing ratio
!
  IF (ktdia < kbeta) THEN

     DO jk=ktdia,kbeta-1

#ifdef __ibmspline__
        CALL prepare_ua_index_spline('cover (2)',kproma,ptm1(1,jk),itv1(1),za(1), &
          &                          pxim1(1,jk),nphase,zphase,itv2)
        CALL lookup_ua_eor_uaw_spline(kproma,itv1(1),za(1),nphase,itv2(1),ua(1))
#else
        CALL prepare_ua_index('cover (2)',kproma,ptm1(1,jk),itv1(1), &
          &                   pxim1(1,jk),nphase,zphase,itv2)
        CALL lookup_ua_eor_uaw(kproma,itv1(1),nphase,itv2(1),ua(1))
#endif

!IBM* novector
        DO jl=1,kproma
           zqsm1(jl) = MIN(ua(jl)*zpapm1i(jl,jk),0.5_wp)
           zcor      = 1._wp/(1._wp-vtmpc1*zqsm1(jl))
           zqsm1(jl) = zqsm1(jl)*zcor
        ENDDO
!
!       Threshold relative humidity, qsat and cloud cover
!       This is from cloud, and is the original calculation for
!       cloud cover, based on relative humidity
!       (Lohmann and Roeckner, Clim. Dyn.  96)
!
      DO jl=1,kproma
!
        zrhc=crt+(crs-crt)*EXP(1._wp-(paphm1(jl,klevp1)                &
             /papm1(jl,jk))**nex)
        zsat=1._wp
        jb=knvb(jl)
        lao=(jb.GE.jbmin .AND. jb.LE.jbmax .AND. pvervel(jl,jk).GT.0._wp)
        lao1=(jk.EQ.jb .OR. jk.EQ.jb+1)
        IF (lao .AND. lao1) THEN
          printop(jl)=REAL(jb,wp)
          zsat=csatsc
        END IF
        zqr=pqm1(jl,jk)/(zqsm1(jl)*zsat)
        paclc(jl,jk)=(zqr-zrhc)/(1.0_wp-zrhc)
        paclc(jl,jk)=MAX(MIN(paclc(jl,jk),1.0_wp),0.0_wp)
        paclc(jl,jk)=1._wp-SQRT(1._wp-paclc(jl,jk))
        IF (pxim1(jl,jk)<=ccwmin .AND. pxlm1(jl,jk)<=ccwmin) &
          paclc(jl,jk)=0.0_wp
      END DO !jl
     END DO  !jk
  END IF ! "pseudo lcover=false"

#ifdef _PROFILE
  CALL trace_stop ('cover', 9)
#endif
!
!
END SUBROUTINE cover

END MODULE mo_cover
