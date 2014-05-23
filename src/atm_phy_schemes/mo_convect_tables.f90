!>
!! @brief Lookup tables for convective adjustment code
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!!  D. Salmond, CRAY (UK), August 1991, original code
!!  L. Kornblueh, MPI, April 2003, cleanup and move of table setup code
!!                                 from setphys.f90 in module
!!  M. Puetz, IBM, April 2009, added spline optio
!!  H. Wan, MPI, 2010-07-16, transfer to ICON
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#ifdef __xlC__
@PROCESS HOT
#else
#define FSEL(a,b,c) MERGE(b,c,(a).GE.0._wp)
#endif

MODULE mo_convect_tables

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: alv, als, cpd, rd, rv, tmelt, vtmpc1
  USE mo_echam_cloud_params, ONLY: csecfrl, cthomi

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: &
    &  version = '$Id$'

  ! variables public

  PUBLIC :: jptlucu1          ! lookup table lower bound
  PUBLIC :: jptlucu2          ! lookup table upper bound
  PUBLIC :: fjptlucu1         ! lookup table lower bound (real)
  PUBLIC :: fjptlucu2         ! lookup table upper bound (real)
  PUBLIC :: tlucua            ! table -- e_s*Rd/Rv
  PUBLIC :: tlucub            ! table -- for derivative calculation: d es/ d t
  PUBLIC :: tlucuc            ! table -- lv/cp
  PUBLIC :: tlucuaw           ! table -- e_s*Rd/Rv - force saturation wrt water
  Public :: tlucubw           ! table -- for derivative calculation: d es/ d t - forcing saturation wrt water
  PUBLIC :: tlucucw           ! table -- lv/cp - forcing saturation wrt water
  PUBLIC :: tlucu             ! fused table
  PUBLIC :: lookupoverflow
  PUBLIC :: flookupoverflow

  PUBLIC :: c1es, c2es, c3les, c3ies, c4les, c4ies, c5les, c5ies, &
    &       c5alvcp, c5alscp, alvdcp, alsdcp

  ! subroutines public

  PUBLIC :: compute_qsat
  PUBLIC :: compute_qsat_amip
  PUBLIC :: init_convect_tables ! initialize LUTs
  PUBLIC :: lookuperror         ! error handling routine
  PUBLIC :: prepare_ua_index_spline
  PUBLIC :: lookup_ua_spline
  PUBLIC :: lookup_uaw_spline
  PUBLIC :: lookup_ua_eor_uaw_spline

  PUBLIC :: lookup_ua_list_spline
  PUBLIC :: lookup_ua_list_spline_2
  PUBLIC :: lookup_uaw_list_spline

  PUBLIC :: prepare_ua_index
  PUBLIC :: lookup_ua
  PUBLIC :: lookup_uaw
  PUBLIC :: lookup_ua_eor_uaw

  PUBLIC :: lookup_ua_list
  PUBLIC :: lookup_uaw_list

  PUBLIC :: lookup_ubc
  PUBLIC :: lookup_ubc_list

  ! Constants used for the computation of lookup tables of the saturation
  ! mixing ratio over liquid water (*c_les*) or ice(*c_ies*)
  !
  REAL (wp), PARAMETER :: c1es  = 610.78_wp              !
  REAL (wp), PARAMETER :: c2es  = c1es*rd/rv             !
  REAL (wp), PARAMETER :: c3les = 17.269_wp              !
  REAL (wp), PARAMETER :: c3ies = 21.875_wp              !
  REAL (wp), PARAMETER :: c4les = 35.86_wp               !
  REAL (wp), PARAMETER :: c4ies = 7.66_wp                !
  REAL (wp), PARAMETER :: c5les = c3les*(tmelt-c4les)    !
  REAL (wp), PARAMETER :: c5ies = c3ies*(tmelt-c4ies)    !
  REAL (wp), PARAMETER :: c5alvcp = c5les*alv/cpd        !
  REAL (wp), PARAMETER :: c5alscp = c5ies*als/cpd        !
  REAL (wp), PARAMETER :: alvdcp  = alv/cpd              !
  REAL (wp), PARAMETER :: alsdcp  = als/cpd              !

  INTEGER, PARAMETER :: jptlucu1 =  50000  ! lookup table lower bound
  INTEGER, PARAMETER :: jptlucu2 = 400000  ! lookup table upper bound

  ! basic table size is 7000, 1/50 of the original size
  ! for the spline tables the melting point tmelt = 273.15
  ! is resembled twice at index = 273150/50 = 5463 and index = 5462
  ! to be able to interpolate values close to the melting point correctly

  INTEGER, PARAMETER :: lucupmin = 1000  ! lookup table lower bound
  INTEGER, PARAMETER :: lucupmax = 8000  ! lookup table upper bound

  ! add float point constants for accelerated boundary guard code
  ! original big table
  REAL(wp), PARAMETER :: fjptlucu1 =  50000._wp
  REAL(wp), PARAMETER :: fjptlucu2 = 400000._wp
  ! new basic table bounds
  REAL(wp), PARAMETER :: flucupmin = 1000._wp
  REAL(wp), PARAMETER :: flucupmax = 8000._wp

  ! logical is fine and logical, but fp faster
  LOGICAL  :: lookupoverflow = .FALSE.          ! preset with false
  REAL(wp) :: flookupoverflow = 0._wp           ! preset with false

  REAL(wp) :: tlucua(jptlucu1-1:jptlucu2+1)     ! table - e_s*Rd/Rv
  REAL(wp) :: tlucub(jptlucu1-1:jptlucu2+1)     ! table - for derivative calculation
  REAL(wp) :: tlucuc(jptlucu1-1:jptlucu2+1)     ! table - l/cp
  REAL(wp) :: tlucuaw(jptlucu1-1:jptlucu2+1)    ! table
  REAL(wp) :: tlucubw(jptlucu1-1:jptlucu2+1)    ! table
  REAL(wp) :: tlucucw(jptlucu1-1:jptlucu2+1)    ! table
  REAL(wp) :: tlucuad(jptlucu1-1:jptlucu2+1)    ! table - e_s*Rd/Rv
  REAL(wp) :: tlucuawd(jptlucu1-1:jptlucu2+1)   ! table

  ! fused tables for splines
  REAL(wp) :: tlucu(1:2,lucupmin-2:lucupmax+1)     ! fused table
  REAL(wp) :: tlucuw(1:2,lucupmin-2:lucupmax+1)    ! fused table

#ifdef __SPLINE_TEST__
  REAL(wp) :: ztemp(jptlucu1+1:jptlucu2-1)
  REAL(wp) :: ua(jptlucu1+1:jptlucu2-1)
  REAL(wp) :: dua(jptlucu1+1:jptlucu2-1)
  REAL(wp) :: uaw(jptlucu1+1:jptlucu2-1)
  REAL(wp) :: duaw(jptlucu1+1:jptlucu2-1)
  REAL(wp) :: sa(jptlucu1+1:jptlucu2-1)
  REAL(wp) :: dsa(jptlucu1+1:jptlucu2-1)
  REAL(wp) :: saw(jptlucu1+1:jptlucu2-1)
  REAL(wp) :: dsaw(jptlucu1+1:jptlucu2-1)
  INTEGER  :: tmp_idx(jptlucu1+1:jptlucu2-1)
  REAL(wp) :: za(jptlucu1+1:jptlucu2-1)
#endif

  !------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE init_convect_tables

    ! Reference: Sonntag D., 1990: Important new values of the physical
    ! constants of 1986, vapour pressure formulations based on ITS-90,
    ! and psychrometer formulae. Z. Meteor. 70, pp 340-344.

    REAL(wp), PARAMETER :: zavl1 = -6096.9385_wp
    REAL(wp), PARAMETER :: zavl2 =    21.2409642_wp
    REAL(wp), PARAMETER :: zavl3 =    -2.711193_wp
    REAL(wp), PARAMETER :: zavl4 =     1.673952_wp
    REAL(wp), PARAMETER :: zavl5 =     2.433502_wp

    REAL(wp), PARAMETER :: zavi1 = -6024.5282_wp
    REAL(wp), PARAMETER :: zavi2 =    29.32707_wp
    REAL(wp), PARAMETER :: zavi3 =     1.0613868_wp
    REAL(wp), PARAMETER :: zavi4 =    -1.3198825_wp
    REAL(wp), PARAMETER :: zavi5 =    -0.49382577_wp

    REAL(wp) :: z5alvcp, z5alscp, zalvdcp, zalsdcp
    REAL(wp) :: ztt, zldcp
    REAL(wp) :: zcvm3, zcvm4, zcvm5
    REAL(wp) :: zavm1, zavm2, zavm3, zavm4, zavm5
    REAL(wp) :: zminner,zdminner,zlinner,zdlinner
#ifdef __SPLINE_TEST__
    REAL(wp) :: max_dua,max_ua,ua_eps,dua_eps
    REAL(wp) :: max_duaw,max_uaw,uaw_eps,duaw_eps
#endif
    INTEGER :: it

    z5alvcp = c5les*alv/cpd
    z5alscp = c5ies*als/cpd

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    ! extend the table by one to safe guard code for accessing it+1, when it is already checked

!IBM* NOVECTOR
    DO it = jptlucu1-1, jptlucu2+1

      ztt = 0.001_wp*REAL(it,wp)
      IF ((ztt-tmelt) > 0.0_wp) THEN
        zcvm3 = c3les
        zcvm4 = c4les
        zcvm5 = z5alvcp
        zldcp = zalvdcp
        zavm1 = zavl1
        zavm2 = zavl2
        zavm3 = zavl3
        zavm4 = zavl4
        zavm5 = zavl5
      ELSE
        zcvm3 = c3ies
        zcvm4 = c4ies
        zcvm5 = z5alscp
        zldcp = zalsdcp
        zavm1 = zavi1
        zavm2 = zavi2
        zavm3 = zavi3
        zavm4 = zavi4
        zavm5 = zavi5
      END IF

      zminner  = (zavm1/ztt+zavm2+zavm3*0.01_wp*ztt+zavm4*ztt*ztt*1.e-5_wp+zavm5*LOG(ztt))
      zdminner = (-zavm1/(ztt*ztt) + zavm3*0.01_wp + zavm4*ztt*2.e-5_wp + zavm5/ztt)

      zlinner  = (zavl1/ztt+zavl2+zavl3*0.01_wp*ztt+zavl4*ztt*ztt*1.e-5_wp+zavl5*LOG(ztt))
      zdlinner = (-zavl1/(ztt*ztt) + zavl3*0.01_wp + zavl4*ztt*2.e-5_wp + zavl5/ztt)

      tlucua(it)  = EXP(zminner)*rd/rv
      tlucuaw(it) = EXP(zlinner)*rd/rv

      tlucuad(it)  = zdminner*EXP(zminner)*rd/rv
      tlucuawd(it) = zdlinner*EXP(zlinner)*rd/rv

      tlucub(it)  = zcvm5*(1.0_wp/(ztt-zcvm4))**2
      tlucubw(it)  = z5alvcp*(1.0_wp/(ztt-c4les))**2

      tlucuc(it)  = zldcp
      tlucucw(it) = zalvdcp

    END DO

!IBM* NOVECTOR
    DO it = lucupmin-2, 27315/5-1
      ztt = 0.05_wp*REAL(it+1,wp)
      zcvm3 = c3ies
      zcvm4 = c4ies
      zcvm5 = z5alscp
      zldcp = zalsdcp
      zavm1 = zavi1
      zavm2 = zavi2
      zavm3 = zavi3
      zavm4 = zavi4
      zavm5 = zavi5

      zminner  = (zavm1/ztt+zavm2+zavm3*0.01_wp*ztt+zavm4*ztt*ztt*1.e-5_wp+zavm5*LOG(ztt))
      zdminner = (-zavm1/(ztt*ztt) + zavm3*0.01_wp + zavm4*ztt*2.e-5_wp + zavm5/ztt)

      zlinner  = (zavl1/ztt+zavl2+zavl3*0.01_wp*ztt+zavl4*ztt*ztt*1.e-5_wp+zavl5*LOG(ztt))
      zdlinner = (-zavl1/(ztt*ztt) + zavl3*0.01_wp + zavl4*ztt*2.e-5_wp + zavl5/ztt)

      tlucu(1,it)  = EXP(zminner)*rd/rv
      tlucu(2,it)  = 0.05_wp*zdminner*EXP(zminner)*rd/rv

      tlucuw(1,it) = EXP(zlinner)*rd/rv
      tlucuw(2,it) = 0.05_wp*zdlinner*EXP(zlinner)*rd/rv

    END DO

!IBM* NOVECTOR
    DO it = 27315/5,lucupmax
      ztt = 0.05_wp*REAL(it,wp)
      zcvm3 = c3les
      zcvm4 = c4les
      zcvm5 = z5alvcp
      zldcp = zalvdcp
      zavm1 = zavl1
      zavm2 = zavl2
      zavm3 = zavl3
      zavm4 = zavl4
      zavm5 = zavl5

      zminner  = (zavm1/ztt+zavm2+zavm3*0.01_wp*ztt+zavm4*ztt*ztt*1.e-5_wp+zavm5*LOG(ztt))
      zdminner = (-zavm1/(ztt*ztt) + zavm3*0.01_wp + zavm4*ztt*2.e-5_wp + zavm5/ztt)

      zlinner  = (zavl1/ztt+zavl2+zavl3*0.01_wp*ztt+zavl4*ztt*ztt*1.e-5_wp+zavl5*LOG(ztt))
      zdlinner = (-zavl1/(ztt*ztt) + zavl3*0.01_wp + zavl4*ztt*2.e-5_wp + zavl5/ztt)

      tlucu(1,it)  = EXP(zminner)*rd/rv
      tlucu(2,it)  = 0.05_wp*zdminner*EXP(zminner)*rd/rv

      tlucuw(1,it) = EXP(zlinner)*rd/rv
      tlucuw(2,it) = 0.05_wp*zdlinner*EXP(zlinner)*rd/rv
    END DO

#ifdef __SPLINE_TEST__

    DO it = jptlucu1-1, jptlucu2-1
      ztemp(it) = 0.001_wp*it
    END DO

    CALL prepare_ua_index('setup1',jptlucu2-jptlucu1-2,ztemp(jptlucu1+1),tmp_idx(jptlucu1+1))
    CALL lookup_ua(jptlucu2-jptlucu1-2,tmp_idx(jptlucu1+1),ua(jptlucu1+1),dua(jptlucu1+1))
    CALL lookup_uaw(jptlucu2-jptlucu1-2,tmp_idx(jptlucu1+1),uaw(jptlucu1+1),duaw(jptlucu1+1))

    max_ua   = 1.e-20_wp
    max_dua  = 1.e-20_wp
    max_uaw  = 1.e-20_wp
    max_duaw = 1.e-20_wp

    DO it = jptlucu1+1, jptlucu2-1

      ua_eps = ABS(tlucua(it) - ua(it))/tlucua(it)
      max_ua = MAX(max_ua,ua_eps)

      uaw_eps = ABS(tlucuaw(it) - uaw(it))/tlucuaw(it)
      max_uaw = MAX(max_uaw,uaw_eps)

      dua_eps = ABS(tlucuad(it) - dua(it))/tlucuad(it)
      max_dua = MAX(max_dua,dua_eps)

      duaw_eps = ABS(tlucuawd(it) - duaw(it))/tlucuawd(it)
      max_duaw = MAX(max_duaw,duaw_eps)

      WRITE (1,*) 0.001_wp*it, tlucua(it),   ua(it),   ua_eps
      WRITE (2,*) 0.001_wp*it, tlucuaw(it),  uaw(it),  uaw_eps
      WRITE (3,*) 0.001_wp*it, tlucuad(it),  dua(it),  dua_eps
      WRITE (4,*) 0.001_wp*it, tlucuawd(it), duaw(it), duaw_eps

    END DO

    WRITE (0,*) 'max_ua   epsilon',max_ua
    WRITE (0,*) 'max_dua  epsilon',max_dua
    WRITE (0,*) 'max_uaw  epsilon',max_uaw
    WRITE (0,*) 'max_duaw epsilon',max_duaw

    !-------------------------------------------------------------------------

    CALL prepare_ua_index_spline('setup2',jptlucu2-jptlucu1-2,ztemp(jptlucu1+1),tmp_idx(jptlucu1+1),za(jptlucu1+1))
    CALL lookup_ua_spline(jptlucu2-jptlucu1-2,tmp_idx(jptlucu1+1),za(jptlucu1+1),sa(jptlucu1+1), dsa(jptlucu1+1))
    CALL lookup_uaw_spline(jptlucu2-jptlucu1-2,tmp_idx(jptlucu1+1),za(jptlucu1+1),saw(jptlucu1+1),dsaw(jptlucu1+1))

    max_ua   = 1.e-20_wp
    max_dua  = 1.e-20_wp
    max_uaw  = 1.e-20_wp
    max_duaw = 1.e-20_wp

    DO it = jptlucu1+1, jptlucu2-1

      ua_eps = ABS(tlucua(it) - sa(it))/tlucua(it)
      max_ua = MAX(max_ua,ua_eps)

      uaw_eps = ABS(tlucuaw(it) - saw(it))/tlucuaw(it)
      max_uaw = MAX(max_uaw,uaw_eps)

      dua_eps = ABS(tlucuad(it) - dsa(it))/tlucuad(it)
      max_dua = MAX(max_dua,dua_eps)

      duaw_eps = ABS(tlucuawd(it) - dsaw(it))/tlucuawd(it)
      max_duaw = MAX(max_duaw,duaw_eps)

      WRITE (11,*) 0.001_wp*it, tlucua(it),   sa(it),   ua_eps
      WRITE (12,*) 0.001_wp*it, tlucuaw(it),  saw(it),  uaw_eps
      WRITE (13,*) 0.001_wp*it, tlucuad(it),  dsa(it),  dua_eps
      WRITE (14,*) 0.001_wp*it, tlucuawd(it), dsaw(it), duaw_eps

    END DO

    WRITE (0,*) 'max_ua   epsilon',max_ua
    WRITE (0,*) 'max_dua  epsilon',max_dua
    WRITE (0,*) 'max_uaw  epsilon',max_uaw
    WRITE (0,*) 'max_duaw epsilon',max_duaw

    DO it = lucupmin, lucupmax
      IF (it < 27315/5) THEN
        WRITE (21,*) 0.05*(it+1), tlucu(1,it)
        WRITE (23,*) 0.05*(it+1), 20.0_wp*tlucu(2,it)
      ELSE
        WRITE (21,*) 0.05*it, tlucu(1,it)
        WRITE (23,*) 0.05*it, 20.0_wp*tlucu(2,it)
      ENDIF
    ENDDO

    STOP 'lookup tables printed'
#endif

  END SUBROUTINE init_convect_tables


  SUBROUTINE lookup_ubc(size,temp,ub,uc)

    INTEGER,            INTENT(in) :: size
    REAL(wp),           INTENT(in) :: temp(size)

    REAL(wp),           INTENT(out) :: ub(size)
    REAL(wp), OPTIONAL, INTENT(out) :: uc(size)

    REAL(wp) :: z5alvcp, z5alscp, zalvdcp, zalsdcp, zcvm4, zcvm5

    INTEGER :: jl

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    z5alvcp = c5les*zalvdcp
    z5alscp = c5ies*zalsdcp

    IF (PRESENT(uc)) THEN
!IBM* NOVECTOR
      DO jl = 1,size
        zcvm4  = FSEL(tmelt-temp(jl),c4ies,c4les)
        zcvm5  = FSEL(tmelt-temp(jl),z5alscp,z5alvcp)
        uc(jl) = FSEL(tmelt-temp(jl),zalsdcp,zalvdcp)
        ub(jl) = zcvm5/(temp(jl)-zcvm4)**2
      END DO
    ELSE
!IBM* NOVECTOR
      DO jl = 1,size
        zcvm4  = FSEL(tmelt-temp(jl),c4ies,c4les)
        zcvm5  = FSEL(tmelt-temp(jl),z5alscp,z5alvcp)
        ub(jl) = zcvm5/(temp(jl)-zcvm4)**2
      END DO
    END IF

  END SUBROUTINE lookup_ubc

  SUBROUTINE lookup_ubc_list(size,kidx,list,temp,ub,uc)

    INTEGER,           INTENT(in) :: size, kidx

    INTEGER,           INTENT(in) :: list(kidx)
    REAL(wp),          INTENT(in) :: temp(size)

    REAL(wp),          INTENT(out) :: ub(size)
    REAL(wp),OPTIONAL, INTENT(out) :: uc(size)

    REAL(wp) :: z5alvcp, z5alscp, zalvdcp, zalsdcp, zcvm4, zcvm5

    INTEGER :: nl, jl

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    z5alvcp = c5les*zalvdcp
    z5alscp = c5ies*zalsdcp

    IF (PRESENT(uc)) THEN
!IBM* NOVECTOR
!CDIR NODEP,VOVERTAKE,VOB
!IBM* ASSERT(NODEPS)
      DO nl = 1,kidx
        jl = list(nl)
        zcvm4  = FSEL(tmelt-temp(jl),c4ies,c4les)
        zcvm5  = FSEL(tmelt-temp(jl),z5alscp,z5alvcp)
        uc(nl) = FSEL(tmelt-temp(jl),zalsdcp,zalvdcp)
        ub(nl) = zcvm5/(temp(jl)-zcvm4)**2
      END DO
    ELSE
!IBM* NOVECTOR
!CDIR NODEP,VOVERTAKE,VOB
!IBM* ASSERT(NODEPS)
      DO nl = 1,kidx
        jl = list(nl)
        zcvm4  = FSEL(tmelt-temp(jl),c4ies,c4les)
        zcvm5  = FSEL(tmelt-temp(jl),z5alscp,z5alvcp)
        ub(nl) = zcvm5/(temp(jl)-zcvm4)**2
      END DO
    END IF
  END SUBROUTINE lookup_ubc_list

  SUBROUTINE fetch_ua_spline(size,idx,zalpha,table,ua,dua)

    INTEGER,            INTENT(in) :: size
    INTEGER,            INTENT(in) :: idx(size)
    REAL(wp),           INTENT(in) :: zalpha(size)
    REAL(wp),           INTENT(in) :: table(1:2,lucupmin-2:lucupmax+1)

    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    REAL(wp) :: a,b,c,d,dx,ddx,x,bxa
    INTEGER  ::  jl

    IF (PRESENT(ua) .AND. .NOT. PRESENT(dua)) THEN

      DO jl = 1,size

        x = zalpha(jl)

        ! derivative and second derivative approximations (2 flops)

        dx   = table(1,idx(jl)+1) - table(1,idx(jl))
        ddx  = table(2,idx(jl)+1) + table(2,idx(jl))

        ! determine coefficients (2 fma + 1 flop)

        a = ddx - 2._wp*dx
        b = 3._wp*dx - ddx - table(2,idx(jl))
        c = table(2,idx(jl))
        d = table(1,idx(jl))

        ! Horner's scheme to compute the spline functions (3 fmas)

        bxa = b + x*a
        ua(jl) = d + x*(c + x*bxa)
      END DO

    ELSE IF (PRESENT(ua) .AND. PRESENT(dua)) THEN

      DO jl = 1,size

        x = zalpha(jl)

        ! derivate and second derivate approximations (2 flops)

        dx   = table(1,idx(jl)+1) - table(1,idx(jl))
        ddx  = table(2,idx(jl)+1) + table(2,idx(jl))

        ! determine coefficients (2 fma + 1 flop)

        a = ddx - 2._wp*dx
        b = 3._wp*dx - ddx - table(2,idx(jl))
        c = table(2,idx(jl))
        d = table(1,idx(jl))

        ! Horner's scheme to compute the spline functions (5 fmas + 1 flop)

        bxa = b + x*a
        ua(jl)  = d + x*(c + x*bxa)
        dua(jl) = 20._wp*(c + x*(3._wp*bxa - b))
      END DO
    END IF

  END SUBROUTINE fetch_ua_spline


  SUBROUTINE fetch_ua(size,idx,table,u)

    INTEGER,  INTENT(in) :: size
    INTEGER,  INTENT(in) :: idx(size)
    REAL(wp), INTENT(in) :: table(jptlucu1-1:jptlucu2+1)

    REAL(wp), INTENT(out) :: u(size)

    INTEGER  ::  jl

    DO jl = 1,size
      u(jl) = table(idx(jl))
    END DO

  END SUBROUTINE fetch_ua


  SUBROUTINE fetch_ua_list_spline(size,store_idx,lookup_idx,zalpha,table,ua,dua)

    INTEGER,            INTENT(in) :: size
    INTEGER,            INTENT(in) :: store_idx(size), lookup_idx(size)
    REAL(wp),           INTENT(in) :: zalpha(size)
    REAL(wp),           INTENT(in) :: table(1:2,lucupmin-1:lucupmax+1)

    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    REAL(wp) :: a,b,c,d,dx,ddx,x,bxa
    INTEGER  ::  jl,nl

    IF (PRESENT(ua) .AND. .NOT. PRESENT(dua)) THEN
!CDIR NODEP,VOVERTAKE,VOB
!IBM* ASSERT(NODEPS)
      DO nl = 1,size

        jl = store_idx(nl)

        x = zalpha(jl)

        ! derivate and second derivate approximations (2 flops)

        dx   = table(1,lookup_idx(jl)+1) - table(1,lookup_idx(jl))
        ddx  = table(2,lookup_idx(jl)+1) + table(2,lookup_idx(jl))

        ! determine coefficients (2 fma + 1 flop)

        a = ddx - 2._wp*dx
        b = 3._wp*dx - ddx - table(2,lookup_idx(jl))
        c = table(2,lookup_idx(jl))
        d = table(1,lookup_idx(jl))

        ! Horner's scheme to compute the spline functions (3 fmas)

        bxa = b + x*a
        ua(jl) = d + x*(c + x*bxa)
      END DO
    ELSE IF (PRESENT(ua) .AND. PRESENT(dua)) THEN
!CDIR NODEP,VOVERTAKE,VOB
!IBM* ASSERT(NODEPS)
      DO nl = 1,size

        jl = store_idx(nl)

        x = zalpha(jl)

        ! derivate and second derivate approximations (2 flops)

        dx   = table(1,lookup_idx(jl)+1) - table(1,lookup_idx(jl))
        ddx  = table(2,lookup_idx(jl)+1) + table(2,lookup_idx(jl))

        ! determine coefficients (2 fma + 1 flop)

        a = ddx - 2._wp*dx
        b = 3._wp*dx - ddx - table(2,lookup_idx(jl))
        c = table(2,lookup_idx(jl))
        d = table(1,lookup_idx(jl))

        ! Horner's scheme to compute the spline functions (5 fmas + 1 flop)

        bxa = b + x*a
        ua(jl)  = d + x*(c + x*bxa)
        dua(jl) = 20._wp*(c + x*(3._wp*bxa - b))
      END DO
    END IF

  END SUBROUTINE fetch_ua_list_spline


  SUBROUTINE fetch_ua_list(size,kidx,store_idx,lookup_idx,table,u)

    INTEGER,  INTENT(in)    :: size, kidx
    INTEGER,  INTENT(in)    :: store_idx(kidx), lookup_idx(size)
    REAL(wp), INTENT(in)    :: table(jptlucu1-1:jptlucu2+1)

    REAL(wp), INTENT(inout) :: u(size)

    INTEGER  ::  jl,nl

!CDIR NODEP,VOVERTAKE,VOB
!IBM* ASSERT(NODEPS)
    DO nl = 1, kidx

      jl = store_idx(nl)
      u(jl) = table(lookup_idx(jl))

    END DO

  END SUBROUTINE fetch_ua_list


  SUBROUTINE lookup_ua_eor_uaw_spline(size,idx,zalpha,nphase,iphase,ua,dua)

    INTEGER,             INTENT(in) :: size, nphase
    INTEGER,             INTENT(in) :: idx(size), iphase(size)
    REAL(wp),            INTENT(in) :: zalpha(size)

    REAL(wp), OPTIONAL , INTENT(out) :: ua(size), dua(size)

    INTEGER :: tmpidx(size)
    INTEGER ::  jl,iw,inw

    IF (nphase == 0) THEN

      !       CALL fetch_ua_spline(size,idx,zalpha,tlucuw,ua,dua) (causes consistency problems)
      DO jl = 1,size
        tmpidx(jl) = jl
      END DO
      CALL fetch_ua_list_spline(size,tmpidx,idx,zalpha,tlucuw,ua,dua)

    ELSE IF (nphase == size) THEN

      !       CALL fetch_ua_spline(size,idx,zalpha,tlucu,ua,dua) (causes consistency problems)
      DO jl = 1,size
        tmpidx(jl) = jl
      END DO
      CALL fetch_ua_list_spline(size,tmpidx,idx,zalpha,tlucu,ua,dua)

    ELSE
      ! mixed case, must build store index

      iw = 1
      inw = size
      DO jl = 1,size
        tmpidx(iw)  = jl  ! lower part of tmpidx() filled with cond = .true.
        tmpidx(inw) = jl  ! upper part of tmpidx() filled with conf = .false.
        iw = iw + iphase(jl)
        inw = inw - (1 - iphase(jl))
      END DO
      iw = iw - 1

      CALL fetch_ua_list_spline(iw     ,tmpidx(1)   ,idx,zalpha,tlucu ,ua,dua)
      CALL fetch_ua_list_spline(size-iw,tmpidx(iw+1),idx,zalpha,tlucuw,ua,dua)
    END IF

  END SUBROUTINE lookup_ua_eor_uaw_spline

  SUBROUTINE lookup_ua_eor_uaw(size,idx,nphase,iphase,ua,dua)

    INTEGER,            INTENT(in) :: size, nphase
    INTEGER,            INTENT(in) :: idx(size), iphase(size)

    REAL(wp),           INTENT(inout) :: ua(size)
    REAL(wp), OPTIONAL, INTENT(inout) :: dua(size)

    INTEGER :: tmpidx(size)
    INTEGER ::  jl, iw, inw

    IF (nphase == 0) THEN
      ! only water

      CALL fetch_ua(size,idx,tlucuaw,ua)
      IF (PRESENT(dua)) CALL fetch_ua(size,idx,tlucuawd,dua)

    ELSE IF (nphase == size) THEN
      ! only ice

      CALL fetch_ua(size,idx,tlucua,ua)
      IF (PRESENT(dua)) CALL fetch_ua(size,idx,tlucuad,dua)

    ELSE
      ! mixed case, must build store index

      iw = 1
      inw = size
      DO jl = 1,size
        tmpidx(iw)  = jl  ! lower part of tmpidx() filled with cond = .true.
        tmpidx(inw) = jl  ! upper part of tmpidx() filled with cond = .false.
        iw = iw + iphase(jl)
        inw = inw - (1 - iphase(jl))
      END DO
      iw = iw - 1

      CALL fetch_ua_list(size, iw     ,tmpidx(1)   ,idx,tlucua, ua)
      CALL fetch_ua_list(size, size-iw,tmpidx(iw+1),idx,tlucuaw, ua)
      IF (PRESENT(dua)) THEN
        CALL fetch_ua_list(size, iw     ,tmpidx(1)   ,idx,tlucuad, dua)
        CALL fetch_ua_list(size, size-iw,tmpidx(iw+1),idx,tlucuawd, dua)
      END IF
    END IF

  END SUBROUTINE lookup_ua_eor_uaw


  SUBROUTINE lookup_ua_spline(size,idx,zalpha,ua,dua)

    INTEGER,            INTENT(in) :: size
    INTEGER,            INTENT(in) :: idx(size)
    REAL(wp),           INTENT(in) :: zalpha(size)

    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    CALL fetch_ua_spline(size,idx,zalpha,tlucu,ua,dua)

  END SUBROUTINE lookup_ua_spline


  SUBROUTINE lookup_ua(size,idx,ua,dua)

    INTEGER,            INTENT(in) :: size
    INTEGER,            INTENT(in) :: idx(size)

    REAL(wp),           INTENT(inout) :: ua(size)
    REAL(wp), OPTIONAL, INTENT(inout) :: dua(size)

    CALL fetch_ua(size,idx,tlucua,ua)
    IF (PRESENT(dua)) CALL fetch_ua(size,idx,tlucuad,dua)

  END SUBROUTINE lookup_ua


  SUBROUTINE lookup_uaw_spline(size,idx,zalpha,ua,dua)

    INTEGER,            INTENT(in) :: size
    INTEGER,            INTENT(in) :: idx(size)
    REAL(wp),           INTENT(in) :: zalpha(size)

    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    CALL fetch_ua_spline(size,idx,zalpha,tlucuw,ua,dua)

  END SUBROUTINE lookup_uaw_spline


  SUBROUTINE lookup_uaw(size,idx,ua,dua)

    INTEGER,            INTENT(in) :: size
    INTEGER,            INTENT(in) :: idx(size)

    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    CALL fetch_ua(size,idx,tlucuaw,ua)
    IF (PRESENT(dua)) CALL fetch_ua(size,idx,tlucuawd,dua)

  END SUBROUTINE lookup_uaw


  SUBROUTINE prepare_ua_index_spline(name,size,temp,idx,zalpha,xi,nphase,zphase,iphase)

    CHARACTER(len=*),   INTENT(in) :: name

    INTEGER,            INTENT(in) :: size
    REAL(wp),           INTENT(in) :: temp(size)
    REAL(wp), OPTIONAL, INTENT(in) :: xi(size)

    INTEGER,            INTENT(out) :: idx(size)
    REAL(wp),           INTENT(out) :: zalpha(size)

    INTEGER,  OPTIONAL, INTENT(out) :: nphase
    REAL(wp), OPTIONAL, INTENT(out) :: zphase(size)
    INTEGER,  OPTIONAL, INTENT(out) :: iphase(size)

    REAL(wp) :: ztt, ztshft, zinbounds, ztmin,ztmax,znphase,ztest
    INTEGER  ::  jl

    zinbounds = 1._wp
    ztmin = flucupmin
    ztmax = flucupmax

    IF (PRESENT(xi)) THEN

      znphase = 0.0_wp

      DO jl = 1,size

        ztshft = FSEL(temp(jl)-tmelt,0._wp,1._wp)
        ztt = 20._wp*temp(jl)
        zalpha(jl) = ztt - DINT(ztt)
        idx(jl) = INT(ztt-ztshft)

        zinbounds = FSEL(ztmin-ztt,0._wp,zinbounds)
        zinbounds = FSEL(ztt-ztmax,0._wp,zinbounds)

        ! check dual phase conditions

        !          lo2         = (ptm1(jl,jk) .LT. cthomi) .OR.                   &
        !                      (ptm1(jl,jk) .LT. tmelt .AND. zxised .GT. csecfrl)

        ztest = FSEL(temp(jl)-tmelt ,0._wp,1.0_wp)
        ztest = FSEL(csecfrl-xi(jl) ,0._wp,ztest)
        ztest = FSEL(temp(jl)-cthomi,ztest,1.0_wp)

        ! normalize ztest to 0 and 1

        iphase(jl) = INT(ztest)
        zphase(jl) = ztest-0.5_wp
        znphase = znphase + ztest
      END DO

      nphase = INT(znphase)

    ELSE
      DO jl = 1,size

        ztshft = FSEL(temp(jl)-tmelt,0._wp,1._wp)
        ztt = 20._wp*temp(jl)
        zalpha(jl) = ztt - DINT(ztt)
        idx(jl) = INT(ztt-ztshft)

        zinbounds = FSEL(ztmin-ztt,0._wp,zinbounds)
        zinbounds = FSEL(ztt-ztmax,0._wp,zinbounds)
      END DO
    END IF

    ! if one index was out of bounds -> print error and exit

    IF (zinbounds == 0._wp) CALL lookuperror(name)

  END SUBROUTINE prepare_ua_index_spline


  SUBROUTINE prepare_ua_index(name,size,temp,idx,xi,nphase,zphase,iphase)


    CHARACTER(len=*),   INTENT(in) :: name

    INTEGER,            INTENT(in) :: size
    REAL(wp),           INTENT(in) :: temp(size)
    REAL(wp), OPTIONAL, INTENT(in) :: xi(size)

    INTEGER,            INTENT(out)  :: idx(size)

    INTEGER,  OPTIONAL, INTENT(out) :: nphase
    REAL(wp), OPTIONAL, INTENT(out) :: zphase(size)
    INTEGER,  OPTIONAL, INTENT(out) :: iphase(size)

    REAL(wp) :: ztt, zinbounds(size), ztmin,ztmax,znphase,ztest(size)
    INTEGER  ::  jl

    ! first compute all lookup indices and check if they are all within allowed bounds

    zinbounds(:) = 1._wp
    ztmin = fjptlucu1
    ztmax = fjptlucu2

    IF (PRESENT(xi)) THEN

      znphase = 0.0_wp

      DO jl = 1,size

        ztt = DNINT(1000._wp*temp(jl))
        idx(jl) = INT(ztt)

        zinbounds(jl) = FSEL(ztmin-ztt,0._wp,zinbounds(jl))
        zinbounds(jl) = FSEL(ztt-ztmax,0._wp,zinbounds(jl))

        ! check dual phase conditions

        !          lo2         = (ptm1(jl,jk) .LT. cthomi) .OR.                   &
        !                      (ptm1(jl,jk) .LT. tmelt .AND. zxised .GT. csecfrl)

        ztest(jl) = FSEL(temp(jl)-tmelt ,0._wp,1.0_wp)
        ztest(jl) = FSEL(csecfrl-xi(jl) ,0._wp,ztest(jl))
        ztest(jl) = FSEL(temp(jl)-cthomi,ztest(jl),1.0_wp)

        ! normalize ztest to 0 and 1

        iphase(jl) = INT(ztest(jl))
        zphase(jl) = ztest(jl)-0.5_wp
        znphase = znphase + ztest(jl)
      END DO

      nphase = INT(znphase)

    ELSE

      DO jl = 1,size

        ztt = DNINT(1000.0_wp*temp(jl))
        idx(jl) = INT(ztt)

        zinbounds(jl) = FSEL(ztmin-ztt,0._wp,zinbounds(jl))
        zinbounds(jl) = FSEL(ztt-ztmax,0._wp,zinbounds(jl))

      END DO
    END IF

    ! if one index was out of bounds -> print error and exit

    IF (ANY(zinbounds(:) == 0._wp)) CALL lookuperror(name)

  END SUBROUTINE prepare_ua_index


  SUBROUTINE lookup_ua_list_spline(name,size,list,temp,ua,dua)

    CHARACTER(len=*),   INTENT(in) :: name
    INTEGER,            INTENT(in) :: size

    INTEGER,            INTENT(in) :: list(size)
    REAL(wp),           INTENT(in) :: temp(size)

    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    INTEGER  :: idx(size)
    REAL(wp) :: zalpha(size)

    REAL(wp) :: ztt, ztshft, zinbounds, ztmax, ztmin
    INTEGER  :: nl, jl

    zinbounds = 1._wp
    ztmin = flucupmin
    ztmax = flucupmax

    ! first compute all lookup indices and check if they are all within allowed bounds

!CDIR NODEP,VOVERTAKE,VOB
!IBM* ASSERT(NODEPS)
    DO nl = 1,size
      jl = list(nl)

      ztshft = FSEL(temp(jl)-tmelt,0._wp,1._wp)
      ztt = 20._wp*temp(jl)
      zalpha(nl) = ztt - DINT(ztt)
      idx(nl) = INT(ztt-ztshft)

      zinbounds = FSEL(ztmin-ztt,0._wp,zinbounds)
      zinbounds = FSEL(ztt-ztmax,0._wp,zinbounds)
    END DO

    ! if one index was out of bounds -> print error and exit

    IF (zinbounds == 0._wp) CALL lookuperror(name)

    CALL fetch_ua_spline(size,idx,zalpha,tlucu,ua,dua)

  END SUBROUTINE lookup_ua_list_spline

  SUBROUTINE lookup_ua_list_spline_2(name, size, kidx, list, temp, ua, dua)
    CHARACTER(len=*),   INTENT(in)  :: name
    INTEGER,            INTENT(in)  :: size, kidx
    INTEGER,            INTENT(in)  :: list(kidx)
    REAL(wp),           INTENT(in)  :: temp(size)
    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)
 
    INTEGER :: idx(size)
    REAL(wp) :: zalpha(size)

    REAL(wp) :: ztt, ztshft, zinbounds, ztmax, ztmin
    INTEGER :: nl, jl

    zinbounds = 1.0_wp
    ztmin = flucupmin
    ztmax = flucupmax

    ! first compute all lookup indices and check if they are all within allowed bounds

!IBM* ASSERT(NODEPS)
    DO nl = 1, kidx
      jl = list(nl)
      ztshft = FSEL(tmelt-temp(jl),1.0_wp,0.0_wp)
! TODO: original code ztt = rsdeltat*temp(jl) and rsdeltat=40._wp
      ztt = 20._wp*temp(jl)
      zalpha(nl) = ztt - AINT(ztt)
      idx(nl) = INT(ztt-ztshft)
      zinbounds = FSEL(ztmin-ztt,0.0_wp,zinbounds)
      zinbounds = FSEL(ztt-ztmax,0.0_wp,zinbounds)
    END DO
    ! if one index was out of bounds -> print error and exit
    IF (zinbounds == 0.0_wp) CALL lookuperror(name)
    CALL fetch_ua_spline(kidx, idx, zalpha, tlucu, ua, dua)

  END SUBROUTINE lookup_ua_list_spline_2

  SUBROUTINE lookup_ua_list(name,size,kidx,list,temp,ua,dua)

    CHARACTER(len=*),   INTENT(in) :: name
    INTEGER,            INTENT(in) :: size, kidx

    INTEGER,            INTENT(in) :: list(kidx)
    REAL(wp),           INTENT(in) :: temp(size)

    REAL(wp), OPTIONAL, INTENT(out) :: ua(size), dua(size)

    INTEGER  :: idx(size)

    REAL(wp) :: ztt, zinbounds(size), ztmax, ztmin
    INTEGER  :: nl, jl

    zinbounds(:) = 1._wp
    ztmin = fjptlucu1
    ztmax = fjptlucu2

    ! first compute all lookup indices and check if they are all within allowed bounds

!CDIR NODEP,VOVERTAKE,VOB
!IBM* ASSERT(NODEPS)
    DO nl = 1,kidx
      jl = list(nl)

      ztt = DNINT(1000._wp*temp(jl))
      idx(nl) = INT(ztt)

      zinbounds(jl) = FSEL(ztmin-ztt,0._wp,zinbounds(jl))
      zinbounds(jl) = FSEL(ztt-ztmax,0._wp,zinbounds(jl))
    END DO

    ! if one index was out of bounds -> print error and exit

    IF (ANY(zinbounds(:) == 0._wp)) CALL lookuperror(name)

    CALL fetch_ua(kidx,idx,tlucua,ua)
    IF (PRESENT(dua)) CALL fetch_ua(kidx,idx,tlucuad,dua)

  END SUBROUTINE lookup_ua_list


  SUBROUTINE lookup_uaw_list_spline(name,size,list,temp,uaw,duaw)

    CHARACTER(len=*),   INTENT(in) :: name
    INTEGER,            INTENT(in) :: size

    INTEGER,            INTENT(in) :: list(size)
    REAL(wp),           INTENT(in) :: temp(size)

    REAL(wp), OPTIONAL, INTENT(out) :: uaw(size), duaw(size)

    INTEGER  :: idx(size)
    REAL(wp) :: zalpha(size)

    REAL(wp) :: ztt, ztshft, zinbounds, ztmax, ztmin
    INTEGER  :: nl, jl

    zinbounds = 1._wp
    ztmin = flucupmin
    ztmax = flucupmax

    ! first compute all lookup indices and check if they are all within allowed bounds

!CDIR NODEP,VOVERTAKE,VOB
!IBM* ASSERT(NODEPS)
    DO nl = 1,size

      jl = list(nl)

      ztshft = FSEL(temp(jl)-tmelt,0._wp,1._wp)
      ztt = 20._wp*temp(jl)
      zalpha(nl) = ztt - DINT(ztt)
      idx(nl) = INT(ztt-ztshft)

      zinbounds = FSEL(ztmin-ztt,0._wp,zinbounds)
      zinbounds = FSEL(ztt-ztmax,0._wp,zinbounds)
    END DO

    ! if one index was out of bounds -> print error and exit

    IF (zinbounds == 0._wp) CALL lookuperror(name)

    CALL fetch_ua_spline(size,idx,zalpha,tlucuw,uaw,duaw)

  END SUBROUTINE lookup_uaw_list_spline


  SUBROUTINE lookup_uaw_list(name,size,list,temp,uaw,duaw)

    CHARACTER(len=*),   INTENT(in) :: name
    INTEGER,            INTENT(in) :: size

    INTEGER,            INTENT(in) :: list(size)
    REAL(wp),           INTENT(in) :: temp(size)

    REAL(wp), OPTIONAL, INTENT(out) :: uaw(size), duaw(size)

    INTEGER  :: idx(size)

    REAL(wp) :: ztt, zinbounds, ztmax, ztmin
    INTEGER  :: nl, jl

    zinbounds = 1._wp
    ztmin = fjptlucu1
    ztmax = fjptlucu2

    ! first compute all lookup indices and check if they are all within allowed bounds

!CDIR NODEP,VOVERTAKE,VOB
!IBM* ASSERT(NODEPS)
    DO nl = 1,size

      jl = list(nl)

      ztt = DNINT(1000._wp*temp(jl))
      idx(nl) = INT(ztt)

      zinbounds = FSEL(ztmin-ztt,0._wp,zinbounds)
      zinbounds = FSEL(ztt-ztmax,0._wp,zinbounds)
    END DO

    ! if one index was out of bounds -> print error and exit

    IF (zinbounds == 0._wp) CALL lookuperror(name)

    CALL fetch_ua(size,idx,tlucuaw,uaw)
    IF (PRESENT(duaw)) CALL fetch_ua(size,idx,tlucuaw,duaw)

  END SUBROUTINE lookup_uaw_list


  SUBROUTINE lookuperror (name)

    USE mo_exception,  ONLY: finish

    CHARACTER(len=*), INTENT(in) :: name

    CALL finish (name, ' lookup table overflow')

  END SUBROUTINE lookuperror

  !>
  !! Compute saturation specific humidity 
  !! from the given temperature and pressure.
  !!
  SUBROUTINE compute_qsat( kproma, kbdim, ppsfc, ptsfc, pqs )

    INTEGER, INTENT(IN)  :: kbdim, kproma
    REAL(wp),INTENT(IN)  :: ppsfc (kbdim)   !< surface pressure
    REAL(wp),INTENT(IN)  :: ptsfc (kbdim)   !< SST
    REAL(wp),INTENT(OUT) :: pqs   (kbdim)   !< saturation specific humidity

    INTEGER  :: itemp  !< temperature*1000
    INTEGER  :: jc     !< column index
    REAL(wp) :: zes    !< (saturation vapour pressure)*Rd/Rv/ps

    !-----
    lookupoverflow = .FALSE.

    DO jc = 1,kproma
      itemp   = NINT(ptsfc(jc)*1000._wp)
      IF (itemp<jptlucu1 .OR. itemp>jptlucu2) lookupoverflow = .TRUE.
      itemp   = MAX(MIN(itemp,jptlucu2),jptlucu1)
      zes     = tlucua(itemp)/ppsfc(jc)
      pqs(jc) = zes/(1._wp-vtmpc1*zes)
    ENDDO

    IF (lookupoverflow) CALL lookuperror ('compute_qsat')

  END SUBROUTINE compute_qsat
  !-------------
  !>
  !! Compute saturation specific humidity
  !! from the given temperature and pressure.
  !!
  SUBROUTINE compute_qsat_amip( kbdim, is, loidx, ppsfc, ptsfc, pqs )

    INTEGER, INTENT(IN)  :: kbdim, is
    INTEGER ,INTENT(IN)  :: loidx(kbdim)!<
    REAL(wp),INTENT(IN)  :: ppsfc (kbdim)   !< surface pressure
    REAL(wp),INTENT(IN)  :: ptsfc (kbdim)   !< SST
    REAL(wp),INTENT(INOUT) :: pqs   (kbdim)   !< saturation specific humidity

    INTEGER  :: jc, jl      !< column index
    REAL(wp) :: zes         !< (saturation vapour pressure)*Rd/Rv/ps
    REAL(wp) :: zpap, zcor
    REAL(wp) :: ua(kbdim)

    !-----
    lookupoverflow = .FALSE.

      CALL lookup_ua_list_spline_2('compute_qsat',kbdim,is,loidx(1), ptsfc(1), ua(1))
!
    DO jc = 1,is
      jl = loidx(jc)
      zpap    = 1._wp/ppsfc(jl)
      zes     = ua(jc)*zpap
      zcor    = 1._wp/(1._wp-vtmpc1*zes)
      pqs(jl) = zes*zcor
    ENDDO
!
    IF (lookupoverflow) CALL lookuperror ('compute_qsat_amip')

  END SUBROUTINE compute_qsat_amip
  !-------------


END MODULE mo_convect_tables
