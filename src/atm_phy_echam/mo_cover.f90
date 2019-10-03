#ifdef __xlC__
@PROCESS HOT
@PROCESS XLF90(NOSIGNEDZERO)
#endif
#include "fsel.inc"

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module diagnoses cloud cover for current timestep
!!
!! @remarks
!!     In ECHAM4 all cloud cover calculations were performed
!!     in routine CLOUD.  The routine diagnosed the cloud cover
!!     from the previous timestep using m1 variables, but then
!!     used a "first guess" calculation of T and Q to give a
!!     cloud cover estimate for the next timestep.  This meant
!!     that the radiation scheme used cloud cover values that were
!!     from a different timestep to the temperature and water vapour
!!     values, and also that the T and Q values were anyway
!!     preliminary.  Finally, the cover calculation was performed
!!     twice each timestep, when one calculation suffices.
!!
!!     This scheme calculates cover diagnostically and is called
!!     once at the beginning of each timestep.  It uses the
!!     standard relative humidity calculation from Lohmann and
!!     Roeckner (96).
!!
!! @references.
!!     Diagnostic CC scheme: Lohmann and Roeckner 96, Clim. Dyn.
!!
!! @author A. Tompkins    MPI-Hamburg        2000
!!         K. Ketelesen   NEC,         April 2002
!!         L. Kornblueh   MPI-Hamburg, April 2002
!!
!! @par Revision History
!!    v2: first working version
!!    v5: lookup table added
!!    v8: zriddr and functions replaced for vectorization
!!    v9: optimizations, longer vector loop and less indirect addressing
!!       - introduction of additional arrays
!!       - change structure if "beta function scheme" IF block
!!       - scattered loops "ictit" and "ictdg" are collected over "kproma" and "klev"
!!       - intoducing additional arrays to hold data in "ictit" and "ictdg"
!!         addressing scheme
!! - Taken from ECHAM6.2, wrapped in module and modified for ICON
!!   by Monika Esch, MPI-M (2013-11)
!! - updated to ECHAM.6.3 and unified for use in ICON; removed Tompkins scheme
!!   by Monika Esch, MPI-M (2015-05)
!!
!
MODULE mo_cover

  USE mo_kind,                 ONLY : wp
  USE mo_physical_constants,   ONLY : vtmpc1, cpd, grav
  USE mo_echam_convect_tables, ONLY : prepare_ua_index_spline,       lookup_ua_eor_uaw_spline,     &
                                      prepare_ua_index_spline_batch, lookup_ua_eor_uaw_spline_batch
USE mo_echam_cov_config,     ONLY : echam_cov_config
USE mo_index_list,           ONLY : generate_index_list

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cover

CONTAINS
  !>
  !!
  SUBROUTINE cover ( jg                                                    & !in
       &           , jb                                                    & !in
       &           , jcs,      kproma,   kbdim, klev, klevp1               & !in
       &           , ktype,    pfrw,     pfri                              & !in
       &           , zf                                                    & !in
       &           , paphm1,   papm1                                       & !in
       &           , ptm1,     pqm1,     pxlm1, pxim1                      & !in
       &           , paclc                                                 & !out
       &           , printop                                               & !out
       &           )
    !---------------------------------------------------------------------------------
    !
    INTEGER, INTENT(in)    :: jg
    INTEGER, INTENT(in)    :: jb                    !< number of block
    INTEGER, INTENT(in)    :: kbdim, klevp1, klev, jcs, kproma
    INTEGER, INTENT(in)    :: ktype(kbdim)          !< type of convection
    REAL(wp),INTENT(in)    :: pfrw(kbdim)         ,&!< water mask
         &                    pfri(kbdim)           !< ice mask
    REAL(wp),INTENT(in)    :: zf(kbdim,klev)      ,&!< geometric height thickness [m]
         &                    paphm1(kbdim,klevp1),&!< pressure at half levels
         &                    papm1(kbdim,klev)   ,&!< pressure at full levels
         &                    pqm1(kbdim,klev)    ,&!< specific humidity
         &                    ptm1(kbdim,klev)    ,&!< temperature
         &                    pxlm1(kbdim,klev)   ,&!< cloud water
         &                    pxim1(kbdim,klev)     !< cloud ice
    REAL(wp),INTENT(out)   :: paclc(kbdim,klev)     !< cloud cover
    REAL(wp),INTENT(out)   :: printop(kbdim)

    INTEGER :: jl, jk, jbm
    INTEGER :: locnt, nl
    REAL(wp):: zdtdz, zcor, zrhc, zsat, zqr, zqx
    INTEGER :: itv1(kproma), itv2(kproma)

    !
    !   Temporary arrays
    !
    REAL(wp)   ::  zdtmin(kbdim), za(kbdim)
    !
    !   Pointers and counters for iteration and diagnostic loop:
    !
    INTEGER :: nphase
    !
    !   variables required for zriddr iteration scheme:
    !
    REAL(wp) :: zjk, zgam
    REAL(wp) :: zqsm1(kbdim,klev)
    REAL(wp) :: ua(kproma)

    LOGICAL :: lao, lao1

    REAL(wp) :: zpapm1i,     ztmp

    REAL(wp) :: zknvb(kbdim),            zphase(kbdim)

    INTEGER :: knvb(kbdim), loidx(kbdim), lomask(kbdim)

#ifdef _OPENACC
    ! DA: for GPU it's much better to collapse the jk and jl loops,
    !     therefore we need to use 2D temporary arrays
    INTEGER  :: idx_batch(kproma,klev), iphase_batch(kproma,klev), nphase_batch(klev)
    REAL(wp) :: ua_batch (kproma,klev)

    REAL(wp) :: za_batch(kbdim,klev), zphase_batch(kbdim,klev)
    REAL(wp) :: zqsm1_s
#endif

    ! Shortcuts to components of echam_cov_config
    !
    INTEGER :: nex, icov, jkscov, jksinv, jkeinv
    REAL(wp):: csatsc, csat, crt, crs, cinv, cqx, clcon
    !
    icov   = echam_cov_config(jg)% icov
    jkscov = echam_cov_config(jg)% jkscov
    jksinv = echam_cov_config(jg)% jksinv
    jkeinv = echam_cov_config(jg)% jkeinv
    csatsc = echam_cov_config(jg)% csatsc
    csat   = echam_cov_config(jg)% csat
    crs    = echam_cov_config(jg)% crs
    crt    = echam_cov_config(jg)% crt
    nex    = echam_cov_config(jg)% nex
    cinv   = echam_cov_config(jg)% cinv
    cqx    = echam_cov_config(jg)% cqx
    clcon  = echam_cov_config(jg)% clcon

    !$ACC DATA PRESENT( ktype, pfrw, pfri, zf, paphm1, papm1, pqm1, ptm1, pxlm1, pxim1, paclc, printop )  &
    !$ACC       CREATE( itv1, itv2, zdtmin, za, zqsm1, ua, zknvb, zphase,                                 &
    !$ACC               knvb, loidx, lomask, idx_batch, iphase_batch, za_batch, zphase_batch,             &
    !$ACC               ua_batch, nphase_batch )

    ! Initialize output arrays
    !
    !   Initialize variables
    !
    !$ACC KERNELS DEFAULT(NONE) ASYNC(1)
    paclc(jcs:kproma,1:klev) = 0.0_wp
    !$ACC END KERNELS
    !
    !$ACC KERNELS DEFAULT(NONE) ASYNC(1)
    printop(jcs:kproma) = 0.0_wp
    !$ACC END KERNELS


    ! Preparations if needed
    !
    SELECT CASE (icov)
       !
    CASE(1,2) ! Calculate the saturation mixing ratio
       !        for relative humidity based schemes
       !
#ifndef _OPENACC
       DO jk = jkscov,klev
          !
          CALL prepare_ua_index_spline(jg,'cover (2)',jcs,kproma,ptm1(:,jk),itv1(:),   &
               &                       za(:),pxim1(:,jk),nphase,zphase,itv2,       &
               &                       klev=jk,kblock=jb,kblock_size=kbdim)
          ! output: itv1=idx, za=zalpha, nphase, zphase, itv2
          CALL lookup_ua_eor_uaw_spline(jcs,kproma,itv1(:),za(:),nphase,itv2(:),ua(:))
          ! output: ua
          !
!IBM* novector

          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG VECTOR PRIVATE( zcor )
          DO jl = jcs,kproma
             zpapm1i      = SWDIV_NOCHK(1._wp,papm1(jl,jk))
             zqsm1(jl,jk) = MIN(ua(jl)*zpapm1i,0.5_wp)
             zcor         = 1._wp/(1._wp-vtmpc1*zqsm1(jl,jk))
             zqsm1(jl,jk) = zqsm1(jl,jk)*zcor       ! qsat
          END DO
          !$ACC END PARALLEL
          !
       END DO   !jk
#else
      CALL prepare_ua_index_spline_batch( jg,'cover (2)',jcs,kproma,klev-jkscov+1,                                    &
                                          ptm1(:,jkscov:),idx_batch(:,jkscov:),za_batch(:,jkscov:),pxim1(:,jkscov:),  &
                                          nphase_batch(jkscov:), zphase_batch(:,jkscov:), iphase_batch(:,jkscov:),    &
                                          kblock=jb, kblock_size=kbdim, opt_need_host_nphase=.FALSE. )

      CALL lookup_ua_eor_uaw_spline_batch( jcs, kproma, klev-jkscov+1,                    &
                                           idx_batch(:,jkscov:), iphase_batch(:,jkscov:), & 
                                           za_batch(:,jkscov:),  ua_batch(:,jkscov:) )


      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = jkscov,klev
        DO jl = jcs,kproma
          zpapm1i = SWDIV_NOCHK(1._wp,papm1(jl,jk))
          zqsm1(jl,jk) = MIN(ua_batch(jl,jk)*zpapm1i,0.5_wp)
          zcor      = 1._wp/(1._wp-vtmpc1*zqsm1(jl,jk))
          zqsm1(jl,jk) = zqsm1(jl,jk)*zcor       ! qsat
        END DO
      END DO
      !$ACC END PARALLEL
#endif
       !
    END SELECT

    
    ! Calculate the cloud cover
    !
    SELECT CASE (icov)
       !
    CASE(0) ! constant cloud cover scheme
       !      Cloud cover is set to the constant value clcon.
       !
       !$ACC KERNELS DEFAULT(NONE) ASYNC(1)
       paclc(jcs:kproma,jkscov:klev) = clcon
       !$ACC END KERNELS
       !
    CASE(1) ! Fractional cloud cover scheme
       !      The relative humidity dependence follows Sundqvist et al. (1989), Eqs.3.11-3.13,
       !      see also Lohmann and Roeckner (1996), Eq.4, and Roeckner et al. (1996), Eq.55 and 56.
       !      The modification of the fractional cloud cover below inversions, to capture marine
       !      stratocumulus clouds follows Mauritsen et al. (2019), Eq.3.
       !      The vertical dependence of the critical relative humidity is fitted to results
       !      obtained by Xu and Krueger (1991). The vertical dependence on pressure is given
       !      by Roeckner et al. (1996), Eq.57.
       !      - Lohmann and Roeckner, Clim. Dyn., 12, 557–572, 1996.
       !      - Mauritsen et al, JAMES, 11, 998-1038, 2019
       !      - Roeckner et al., MPI-Rport 218, 90pp., 1996.
       !      - Sundqvist et al., Mon. Wea. Rev., 117, 1641-1657, 1989.
       !      - Xu and Krueger., Mon. Wea. Rev., 119, 342-367, 1991.
       !
       !   Initialize variables
       !
       !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs,kproma
          zdtmin(jl) = -cinv * grav/cpd   ! fraction of dry adiabatic lapse rate
          zknvb(jl)= 1.0_wp
       END DO
       !$ACC END PARALLEL
       !
       !   Checking occurrence of low-level inversion
       !   (below 2000 m, sea points only, no convection)
       !
       ! Build index list for columns, which have >50% sea surface, practically no sea ice, and
       ! which are not convective.
       ! Indixes of these columns are stored in the index list loidx, with list indices jcs:locnt.
       !
       !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs,kproma
          lomask(jl) = MERGE( 1, 0, (pfrw(jl).GT.0.5_wp.AND.pfri(jl).LT.1.e-12_wp.AND.ktype(jl).EQ.0) )
       END DO
       !$ACC END PARALLEL

       CALL generate_index_list(lomask, loidx, jcs, kproma, locnt, 1)

       !
       ! For these columns, search from the lowermost layer jk=klev up to jk=jksinv at ~2000m height.
       ! Find the level index with the least negative lapse rate towards
       ! the adjacent upper layer, supposed to be the layer below the inversion.
       IF (locnt.GT.jcs-1) THEN
          zjk = REAL(jk,wp)
          !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
          !$ACC LOOP SEQ
          DO jk = klev,jksinv,-1
             !
!IBM* ASSERT(NODEPS)
            !$ACC LOOP GANG VECTOR
             DO nl = jcs,locnt
                jl = loidx(nl)
                ! Lapse rate dT/dz (K/m) between layers k-1 and k.
                ztmp        = (ptm1(jl,jk-1)-ptm1(jl,jk))/(zf(jl,jk-1)-zf(jl,jk))
                ! Truncate lapse rate dT/dz (K/m) to value <= 0.
                zdtdz       = MIN(0.0_wp, ztmp)
                ! Update inversion level index if the lapse rate is weaker than zdtmin,
                ! otherwise keep old value (initial value = 1)
                zknvb(jl)   = FSEL(zdtmin(jl)-zdtdz,zknvb(jl),zjk)
                ! Update minimum lapse rate, if the lapse rate is less negative,
                ! otherwise keep the old value (initial value = -cinv*grav/cpd)
                zdtmin(jl)  = MAX(zdtdz,zdtmin(jl))
             END DO
          END DO
          !$ACC END PARALLEL
       END IF
       !
       !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs,kproma
          knvb(jl) = INT(zknvb(jl))
       END DO
       !$ACC END PARALLEL
       !
       !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
       !$ACC LOOP GANG VECTOR COLLAPSE(2)
       DO jk = jkscov,klev
          DO jl = jcs,kproma
             !
             ! Scaling factor zsat for the saturation mass mixing ratio, with range [csatsc,1].
             ! Scaling factors < 1 are computed for cells below an inversion in conditions
             ! allowing for marine stratocumulus clouds.
             !
             jbm=knvb(jl)                           ! if inversion was found: level below inversion, otherwise = 1 
             lao=(jksinv <= jbm .AND. jbm <= jkeinv)! true if jbm of this column is in valid range of height
             lao1=(jk == jbm)                       ! true if jk is the level just below the inversion
             IF (lao .AND. lao1) THEN               ! this is a layer below an inversion, where zsat needs to be modified
                printop(jl)=REAL(jk,wp)             ! store level index
                zdtdz = (ptm1(jl,jbm-1)-ptm1(jl,jbm))/(zf(jl,jk-1)-zf(jl,jk)) ! lapse rate dT/dz (K/m)
                zgam  = MAX(0.0_wp,-zdtdz*cpd/grav) ! ratio (dT/dz)/(dry adiab. dT/dz), truncated >= 0
                zsat  = MIN(1.0_wp,csatsc+zgam)     ! scaling factor for saturation mixing ratio, with range [csatsc,1]
                !                                     Mauritsen et al. (2019), Eq.3
             ELSE
                !
                zsat  = 1.0_wp
                !
             END IF
             !
             ! Relative humidity
             zqr=pqm1(jl,jk)/(zqsm1(jl,jk)*zsat)    ! r in Lohmann-scheme (grid-mean rel hum)
             !                                        but scaled by zsat in layer below inversion
             !                                        following Mauritsen et al. (2019)
             !
             ! Critical relative humidity for cloud formation
             zrhc=crt+(crs-crt)*EXP(1._wp-(paphm1(jl,klevp1)/papm1(jl,jk))**nex)
             !
             ! Compute fractional cloud cover
             paclc(jl,jk)=(zqr-zrhc)/(csat-zrhc)    ! = b_o in Eq.4 of Lohmann and Roeckner (1996), linear in zqr
             !                                        = 0 for zqr=zrhc
             !                                        = 1 for zqr=csat
             !
             paclc(jl,jk)=MAX(MIN(paclc(jl,jk),1.0_wp),0.0_wp) ! limit to range [0,1]
             !                                                   = 0                     , zqr<=zrhc
             !                                                   = (zqr-zrhc)/(csat-zrhc), zrhc<zqr<csat
             !                                                   = 1                     , csat<=zqr
             paclc(jl,jk)=1._wp-SQRT(1._wp-paclc(jl,jk))       ! = b in Eq.4 in LR (1996)
             !                                                   = b in Eq.3.13 of Sundqvist et al. (1989)
             !
          END DO  !jl
       END DO   !jk
       !$ACC END PARALLEL
       !
    CASE(2) ! 0/1 cloud cover scheme based on relative humidity.
       !      Cloud cover is 1 if the relative humidity is >= csat, and 0 otherwise
       !
       !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
       !$ACC LOOP GANG VECTOR COLLAPSE(2)
       DO jk = jkscov,klev
          DO jl = jcs,kproma
             !
             zqr = pqm1(jl,jk)/zqsm1(jl,jk)
             !
             IF (zqr >= csat) THEN
                paclc(jl,jk) = 1.0_wp
             END IF
             !
          END DO  !jl
       END DO   !jk
       !$ACC END PARALLEL
       !
    CASE(3) ! 0/1 cloud cover scheme based on cloud condensate.
       !      Cloud cover is 1 if the cloud condensate mixing ration is >= cqx, and 0 otherwise
       !
       !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
       !$ACC LOOP GANG VECTOR COLLAPSE(2)
       DO jk = jkscov,klev
          DO jl = jcs,kproma
             !
             zqx = pxlm1(jl,jk)+pxim1(jl,jk)
             !
             IF (zqx >= cqx) THEN
                paclc(jl,jk) = 1.0_wp
             END IF
             !
          END DO  !jl
       END DO   !jk
       !$ACC END PARALLEL
       !
    END SELECT
    !$ACC END DATA
    !
    !
  END SUBROUTINE cover

END MODULE mo_cover
