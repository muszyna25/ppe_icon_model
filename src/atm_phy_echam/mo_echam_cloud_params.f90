!>
!! Contains all the relevant constants for Tompkins' cloud cover scheme.
!! and the tunable parameters for the cloud microphysics.
!! Subroutines for intializing these values are also included.
!!
!! @author Adrian Tompkins
!!
!! @par Revision History
!!  Originally the module "mo_cloud" from ECHAM, by A. Tompkins (2000-07)
!!  Transfered to and modified for ICON by Hui Wan, MPI (2010-07)
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
MODULE mo_echam_cloud_params

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: tmelt, g=>grav
  USE mo_datetime,           ONLY: rdaylen
  USE mo_exception,          ONLY: finish, print_value

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_cloud_params'

  PUBLIC :: sucloud, init_cloud_tables

  PUBLIC :: cbeta_cs
  PUBLIC :: cbeta_pq,cbeta_pq_max,cvarmin
  PUBLIC :: ctaus,ctaul,ctauk
  PUBLIC :: cmmrmax,ceffmax
  PUBLIC :: nbetax,nbetaq,cbetaqs,rbetak
  PUBLIC :: tbetai,tbetai0,tbetai1
  PUBLIC :: cthomi,cn0s,crhoi,crhosno
  PUBLIC :: ccsaut
  PUBLIC :: clmin,clmax
  PUBLIC :: crs,crt,nex
  PUBLIC :: cqtmin
  PUBLIC :: cptop, cpbot, ncctop, nccbot
  PUBLIC :: jbmin, jbmin1, jbmax, fjbmin, fjbmin1, fjbmax
  PUBLIC :: lonacc

  PUBLIC :: csatsc
  PUBLIC :: ceffmin
  PUBLIC :: csecfrl
  PUBLIC :: ccraut,ccsacl,ccracl,cauloc
  PUBLIC :: cvtfall
  PUBLIC :: ccwmin

  !----------------------------------------
  ! default values for cloud microphysics
  !----------------------------------------
  REAL(wp)            :: cthomi  = tmelt-35.0_wp
  REAL(wp), PARAMETER :: cn0s    = 3.e6_wp
  REAL(wp), PARAMETER :: crhoi   = 500.0_wp
  REAL(wp), PARAMETER :: crhosno = 100.0_wp
  REAL(wp)            :: ccsaut  = 95.0_wp
  REAL(wp), PARAMETER :: clmax   = 0.5_wp
  REAL(wp), PARAMETER :: clmin   = 0.0_wp
  REAL(wp), PARAMETER :: crs     = 0.9_wp
  REAL(wp), PARAMETER :: crt     = 0.7_wp
  INTEGER,  PARAMETER :: nex     = 4
  REAL(wp), PARAMETER :: ceffmax = 150.0_wp   ! max eff.radius for ice cloud

  !---------------------------------------
  ! default values for cloud cover scheme
  !---------------------------------------

  REAL(wp), PARAMETER :: cbeta_cs     = 10.0_wp                  ! K1: conv source of skew
  REAL(wp), PARAMETER :: ctaus        = 1.0_wp/( 0.5_wp*rdaylen) ! htau shortest timescale
  REAL(wp), PARAMETER :: ctaul        = 1.0_wp/(20.0_wp*rdaylen) ! htau longest timescale
  REAL(wp), PARAMETER :: ctauk        = 0.091625_wp              ! htau K = sqrt(3)*Cs(=0.23)^2.
  REAL(wp), PARAMETER :: cbeta_pq     = 2.0_wp                   ! q_0: target value for q
  REAL(wp), PARAMETER :: cbeta_pq_max = 50.0_wp                  ! max values for q
  REAL(wp), PARAMETER :: cvarmin      = 0.1_wp                   ! b-a_0: min dist width *qv
  REAL(wp), PARAMETER :: cmmrmax      = 0.005_wp                 ! max mmr of cld in cldy region
  REAL(wp), PARAMETER :: cqtmin       = 1.e-12_wp                ! total water minimum

  REAL(wp), PARAMETER :: cptop        = 1000.0_wp   ! min. pressure level for cond.
  REAL(wp), PARAMETER :: cpbot        = 50000.0_wp  ! max. pressure level for tropopause calc.

  LOGICAL,  PARAMETER :: lonacc = .TRUE.

  !-------------------------------------------------------------
  ! parameters initialized in subroutine sucloud of this module
  !-------------------------------------------------------------
  REAL(wp) :: csatsc
  REAL(wp) :: ceffmin
  REAL(wp) :: csecfrl
  REAL(wp) :: ccraut,ccsacl,ccracl,cauloc
  REAL(wp) :: cvtfall
  REAL(wp) :: ccwmin

  INTEGER  :: ncctop           ! max. level for condensation
  INTEGER  :: nccbot           ! lowest level for tropopause calculation
  INTEGER  :: jbmin
  INTEGER  :: jbmin1
  INTEGER  :: jbmax
  REAL(wp) :: fjbmin
  REAL(wp) :: fjbmin1
  REAL(wp) :: fjbmax

  !-----------------------------------------------------------------------
  ! lookup table (set in subroutine init_lookup_table_cld of this module)
  !-----------------------------------------------------------------------
  INTEGER, PARAMETER :: nbetax = 400         ! lookup table size for ibeta
  INTEGER, PARAMETER :: nbetaq = 50          ! lookup table size for ibeta
  REAL(wp),PARAMETER :: cbetaqs = 6.0_wp     ! stretch factor for q

  REAL(wp) :: tbetai(0:1,0:nbetaq,0:nbetax+1) ! betai table for q=cbeta_pq
  REAL(wp) :: tbetai0(0:nbetaq,0:nbetax+1)    ! betai table for q=cbeta_pq
  REAL(wp) :: tbetai1(0:nbetaq,0:nbetax+1)    ! betai table for q=cbeta_pq
  REAL(wp) :: rbetak

CONTAINS
  !>
  !!
  !! @author E. Roeckner, MPI, October 2001
  !!
  SUBROUTINE sucloud ( nlev, vct  &
!!$    &                , lmidatm    &
    &                , lcouple    &
    &                , lipcc      &
!!$    &                , lham       &
    &                )

    INTEGER,INTENT(IN)  :: nlev  !< total # of vertical layers
    REAL(wp),INTENT(IN) :: vct(2*(nlev+1))
    LOGICAL, INTENT(IN) :: lcouple, lipcc
!!$    LOGICAL, INTENT(IN) :: lmidatm, lham

    REAL(wp) :: za, zb, zph(nlev+1), zp(nlev), zh(nlev)
    INTEGER  :: jk

    ! Resolution dependent parameters.
    ! A single set of values are hardwired here. Tuning needed for ICON

    ccsacl  = 0.1_wp
    ccracl  = 6.0_wp
    csecfrl = 5.e-7_wp
    ccraut  = 15.0_wp
    cvtfall = 3.29_wp
    ceffmin = 10.0_wp    ! min eff.radius for ice cloud
    ccwmin  = 1.e-7_wp   ! cloud water limit for cover>0
    csatsc  = 1.0_wp
    cauloc  = 5.0_wp

    IF (lcouple .OR. lipcc) THEN
       csatsc  = 1.0_wp
       cauloc  = 0.0_wp
    ENDIF

!!echam!
!!echam! Special 11-Level values
!!echam!
!!echam
!!echam  IF (nlev == 11) THEN
!!echam    ccsacl  = 0.5_wp
!!echam    ccracl  = 6.0_wp
!!echam    csecfrl = 1.e-6_wp
!!echam    ccraut  = 30.0_wp
!!echam    cvtfall = 7.0_wp
!!echam    ceffmin = 30.0_wp    ! min eff.radius for ice cloud
!!echam    ccwmin  = 5.e-7_wp  ! cloud water limit for cover>0
!!echam    csatsc  = 1.0_wp
!!echam  ELSE
!!echam    ccsacl  = 0.10_wp
!!echam    ccracl  = 6.0_wp
!!echam    csecfrl = 5.e-7_wp
!!echam    IF(lham) csecfrl = 1.e-8_wp
!!echam    ccraut  = 15.0_wp
!!echam    cvtfall = 3.29_wp
!!echam    ceffmin = 10.0_wp    ! min eff.radius for ice cloud
!!echam    ccwmin  = 1.e-7_wp  ! cloud water limit for cover>0
!!echam    csatsc  = 1.0_wp
!!echam  ENDIF
!!echam!
!!echam!                19 Level, no middle atmosphere
!!echam!
!!echam  IF (nlev == 11  .AND. .NOT. lmidatm) THEN
!!echam    IF (nn == 21) THEN
!!echam      cauloc  = 5.0_wp
!!echam    ELSE IF (nn == 31) THEN
!!echam      cauloc  = 5.0_wp
!!echam    ELSE
!!echam      CALL finish (TRIM(thismodule), 'Truncation not supported.')
!!echam    ENDIF
!!echam
!!echam  ELSE IF (nlev == 19 .AND. .NOT. lmidatm) THEN
!!echam    IF (nn == 31) THEN
!!echam      cauloc  = 1.0_wp
!!echam    ELSE IF (nn == 42) THEN
!!echam      cauloc  = 2.0_wp
!!echam      IF(lham) cauloc  = 1.0_wp
!!echam    ELSE
!!echam      CALL finish (TRIM(thismodule), 'Truncation not supported.')
!!echam    ENDIF
!!echam!
!!echam!                31 Level, no middle atmosphere
!!echam!
!!echam  ELSE IF(nlev == 31  .AND. .NOT. lmidatm) THEN
!!echam    IF (nn == 63) THEN
!!echam      cauloc  = 5.0_wp
!!echam    ELSE
!!echam      CALL finish (TRIM(thismodule), 'Truncation not supported.')
!!echam    ENDIF
!!echam!
!!echam!                39 Level, middle atmosphere
!!echam!
!!echam  ELSE IF (nlev == 39 .AND. lmidatm) THEN
!!echam    IF (nn == 31) THEN
!!echam      cauloc  = 3.0_wp
!!echam    ELSE
!!echam      CALL finish (TRIM(thismodule), 'Truncation not supported.')
!!echam    ENDIF
!!echam!
!!echam!                47 Level, middle atmosphere
!!echam!
!!echam  ELSE IF(nlev == 47  .AND. lmidatm) THEN
!!echam    IF (nn == 63) THEN
!!echam      cauloc  = 5.0_wp
!!echam    ELSE
!!echam      CALL finish (TRIM(thismodule), 'Truncation not supported.')
!!echam    ENDIF
!!echam!
!!echam!                95 Level, middle atmosphere
!!echam!
!!echam  ELSE IF(nlev == 95  .AND. lmidatm) THEN
!!echam    IF (nn == 127) THEN
!!echam      cauloc  = 5.0_wp
!!echam    ELSE
!!echam      CALL finish (TRIM(thismodule), 'Truncation not supported.')
!!echam    ENDIF
!!echam!
!!echam!                199 Level, middle atmosphere
!!echam!
!!echam  ELSE IF(nlev == 199  .AND. lmidatm) THEN
!!echam    IF (nn == 255) THEN
!!echam      cauloc  = 5.0_wp
!!echam    ELSE
!!echam      CALL finish (TRIM(thismodule), 'Truncation not supported.')
!!echam    ENDIF
!!echam  ELSE
!!echam    CALL finish (TRIM(thismodule), 'Truncation not supported.')
!!echam  ENDIF
!!echam!
!!echam!-- overwrite values for coupled runs
!!echam!
!!echam!!echam  IF(lcouple .OR. lipcc) THEN
!!echam    csatsc  = 1.0_wp
!!echam    cauloc  = 0.0_wp
!!echam  ENDIF
!
!
!-- half level pressure values, assuming 101320. Pa surface pressure

  DO jk=1,nlev+1
    za=vct(jk)
    zb=vct(jk+nlev+1)
    zph(jk)=za+zb*101320.0_wp
  END DO
!
! -- full level pressure
!
  DO jk = 1, nlev
    zp(jk)=(zph(jk)+zph(jk+1))*0.5_wp
  END DO
!
  DO jk = 1, nlev
    zh(jk)=(zph(nlev+1)-zp(jk))/(g*1.25_wp)
  END DO
!
! -- search for highest inversion level (first full level below 1000 m)
!
  DO jk = 1, nlev
    jbmin=jk
    IF(zh(jk).LT.1000.0_wp) EXIT
  END DO
  fjbmin = REAL(jbmin,wp)
!
! -- search for lowest inversion level (first full level below 500 m)
!
  DO jk = 1, nlev
    jbmax=jk
    IF(zh(jk).LT.500.0_wp) EXIT
  END DO
  fjbmax = REAL(jbmax,wp)
!
  jbmin1  = jbmin-1
  fjbmin1 = REAL(jbmin1,wp)
!
! -- search for pressure level cptop (Pa)
!
  DO jk = 1, nlev
    ncctop=jk
    IF(zp(jk).GE.cptop) EXIT
  END DO
!
  CALL print_value ('highest level for condensation: ncctop = ', ncctop)
!
! -- search for pressure level cpbot (Pa)
!
  DO jk = 1, nlev
    nccbot=jk
    IF(zp(jk).GE.cpbot) EXIT
  END DO
!
  CALL print_value('lowest level for tropopause calc.: nccbot = ', nccbot)

  END SUBROUTINE sucloud
  !-------------
  !>
  !! @par Revision History
  !! Computation taken from subroutine "setphys" of ECHAM6 by Hui Wan (2010-07-20)
  !!
  SUBROUTINE init_cloud_tables

    REAL(wp) :: zx, zq
    INTEGER  :: it, iq

    rbetak = (cbeta_pq_max-cbeta_pq)/(EXP(cbetaqs)-1.0_wp)
    DO it = 0, nbetax
      zx = REAL(it,wp)/REAL(nbetax,wp)
      DO iq = 0, nbetaq
        zq = rbetak*(EXP(cbetaqs*REAL(iq,wp)/REAL(nbetaq,wp))-1.0_wp)+cbeta_pq
        tbetai0(iq,it) = betai(cbeta_pq       ,zq,zx)
        tbetai1(iq,it) = betai(cbeta_pq+1.0_wp,zq,zx)
        tbetai(0,iq,it) = betai(cbeta_pq       ,zq,zx)
        tbetai(1,iq,it) = betai(cbeta_pq+1.0_wp,zq,zx)
      ENDDO
    ENDDO

    ! mpuetz: add nbetax+1 entry to avoid checking for nebtax in spline code
    tbetai(:,:,nbetax+1) = tbetai(:,:,nbetax)
    tbetai0(:,:) = tbetai(0,:,:)
    tbetai1(:,:) = tbetai(1,:,:)

  END SUBROUTINE init_cloud_tables
  !-------------
  !>
  !! @par Uses betacf, gammln; returns the incomplete beta function I x (a; b).
  !!
  !! @par Method
  !!   See Numerical Recipes (Fortran)
  !!
  !! @author  A. Tompkins, MPI, 2000
  !!
  FUNCTION betai(p,q,x)

    REAL(wp)             :: betai
    REAL(wp), INTENT(in) :: p,q, & ! beta shape parameters
                            x      ! integration limit
    !  local scalars:
    REAL(wp) :: bt !,betacf,gammln

    IF (x > 0.0_wp .AND. x < 1.0_wp ) THEN  ! factors in front of the continued fraction.
      bt = EXP(gammln(p+q)-gammln(p)-gammln(q)+p*LOG(x)+q*LOG(1.0_wp-x))
    ELSE
      bt = 0.0_wp
    ENDIF
    IF (x < (p+1.0_wp)/(p+q+2.0_wp)) THEN ! use continued fraction directly.
      betai = bt*betacf(p,q,x)/p
    ELSE     ! use continued fraction after making the symmetry transformation.
      betai = 1.0_wp-bt*betacf(q,p,1.0_wp-x)/q !
    ENDIF

  END FUNCTION betai
  !-----------
  !>
  !! @par Description
  !!  used by betai: evaluates continued fraction for incomplete
  !!  beta function by modi ed lentz's method ( x 5.2).
  !!  first step of lentz's method.
  !!
  !! @par Method
  !!   See Numerical Recipes (Fortran)
  !!
  !! @author  A. Tompkins, MPI, 2000
  !!
  FUNCTION betacf(p,q,x)

    REAL(wp)             :: betacf
    REAL(wp), INTENT(in) :: p,q, & ! beta shape parameters
                            x      ! integration limit

    INTEGER :: maxit = 100, m, m2
    REAL(wp) :: zeps = 3.e-7_wp, fpmin = 1.e-30_wp, &
                aa, c, d, del, h, qab, qam, qap

    qab = p+q

    ! these q's will be used in factors that occur in the coe cients (6.4.6).

    qap = p+1.0_wp
    qam = p-1.0_wp
    c = 1.0_wp
    d = 1.0_wp-qab*x/qap
    IF (ABS(d) < fpmin) d = fpmin
    d = 1.0_wp/d
    h = d
    m = 1
    del = 2.0_wp
    DO WHILE (ABS(del-1.0_wp) > zeps)
      m2 = 2*m
      aa = REAL(m,wp)*(q-REAL(m,wp))*x/((qam+REAL(m2,wp))*(p+REAL(m2,wp)))
      d = 1.0_wp+aa*d  ! one step (the even one) of the recurrence.
      IF (ABS(d) < fpmin) d = fpmin
      c = 1.0_wp+aa/c
      IF (ABS(c) < fpmin) c = fpmin
      d = 1.0_wp/d
      h = h*d*c
      aa = -(p+REAL(m,wp))*(qab+REAL(m,wp))*x/((p+REAL(m2,wp))*(qap+REAL(m2,wp)))
      d = 1.0_wp+aa*d ! next step of the recurrence (the odd one).
      IF (ABS(d) < fpmin) d = fpmin
      c = 1.0_wp+aa/c
      IF (ABS(c) < fpmin) c = fpmin
      d = 1.0_wp/d
      del = d*c
      h = h*del
      m = m+1            ! AMT
      IF (m > maxit) THEN
        CALL finish ('betacf','a or b too big, or maxit too small in betacf')
        del = 1.0_wp
      ENDIF
    ENDDO
    betacf = h

  END FUNCTION betacf
  !-----------
  !>
  !! @par Description
  !!  Gamma function calculation. Returns the value ln[g(xx)] for xx > 0.
  !!
  !! @par Method
  !!   See Numerical Recipes (Fortran)
  !!
  !! @author  A. Tompkins, MPI, 2000
  FUNCTION gammln(xx)

    REAL(wp)             :: gammln
    REAL(wp), INTENT(in) :: xx

    INTEGER :: j

    REAL(wp) :: ser, tmp, x, y
    REAL(wp), PARAMETER :: cof(6) = (/ &
         76.18009172947146_wp, -86.50532032941677_wp, &
         24.01409824083091_wp, -1.231739572450155_wp, &
         0.1208650973866179e-2_wp, -0.5395239384953e-5_wp /)
    REAL(wp), PARAMETER :: stp = 2.5066282746310005_wp

    x = xx
    y = x
    tmp = x+5.5_wp
    tmp = (x+0.5_wp)*LOG(tmp)-tmp
    ser = 1.000000000190015_wp
    DO j =1, 6
      y = y+1.0_wp
      ser = ser+cof(j)/y
    ENDDO
    gammln = tmp+LOG(stp*ser/x)

  END FUNCTION gammln
  !-----------

END MODULE mo_echam_cloud_params
