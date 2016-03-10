!>
!! @file chemcon.f90
!! @brief Inorganic carbon cycle.
!!
!! Computes chemical constants in the surface layer (aksurf, former chemc)
!! and in the water column (ak13, ak23, akb3, aksp)
!!
#include "hamocc_omp_definitions.inc"

SUBROUTINE CHEMCON ( start_idx, end_idx, klevs, psao, ptho,  &
     &               pddpo, ptiestu, kldtday)

  USE mo_kind, ONLY        : wp
  USE mo_carbch, ONLY      : akb3, akw3, ak13, ak23, aksp, &
       &                     satn2o, satoxy, satn2, solco2,&
       &                     aksurf
  USE mo_biomod, ONLY      : rrrcl
  USE mo_sedmnt, ONLY      : calcon
  USE mo_control_bgc, ONLY : bgc_nproma, bgc_zlevs


  USE mo_bgc_constants

  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in)  :: start_idx            !< start index for j loop (ICON cells, MPIOM lat dir)  
  INTEGER, INTENT(in)  :: end_idx              !< end index  for j loop  (ICON cells, MPIOM lat dir) 
  INTEGER, INTENT(in)  :: klevs(bgc_nproma)    !<  vertical levels

  REAL(wp) :: psao(bgc_nproma,bgc_zlevs)  !< salinity [psu]
  REAL(wp) :: ptho(bgc_nproma,bgc_zlevs)  !< potential temperature [deg C]
  REAL(wp) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp) :: ptiestu(bgc_nproma,bgc_zlevs)       !< depth of scalar grid cell [m]

  !! Local variables

  INTEGER ::  jc, k, kldtday, kpke
  REAL(wp) :: t, q, s, log_t, log_q, sqrt_s, ti, qi, q2, tc2
  REAL(wp) :: cek0, ckb, ck1, ck2, ckw, oxy, ani
  REAL(wp) :: ak1, ak2, akb, akw, ak0, aksp0, log10ksp
  REAL(wp) :: p, cp, tc
  REAL(wp) :: rs

  !     -----------------------------------------------------------------
  !*            SET MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG)
  !             (SEE BROECKER A. PENG, 1982, P. 26)
  !             ([CA++](MOLES/KG)=1.026E-2*(S/35.) AFTER
  !             CULKIN(1965), CF. BROECKER ET AL. 1982)
  !             ------------- --- -------- -- --- -----
  !
  calcon = 1.03e-2_wp




  rrrcl = salchl * 1.025_wp * bor1 * bor2

  !
  !     -----------------------------------------------------------------
  !*        21. CHEMICAL CONSTANTS - SEA SURFACE
  !             --------------------------------
  !
!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(jc,T,log_t,ti,q,log_q,qi,q2,s,sqrt_s,&
!HAMOCC_OMP            cek0,ak0,rs,oxy,ani,ck1,ck2,ckb,ckw,ak1,ak2,&
!HAMOCC_OMP            akw,akb) HAMOCC_OMP_DEFAULT_SCHEDULE

 DO jc = start_idx, end_idx

       IF (pddpo(jc, 1) > 0.5_wp) THEN           ! wet cell
              !
              !*        21.1 SET ABSOLUTE TEMPERATURE
              !              ------------------------
              T = ptho(jc,1) + tmelt                  ! degC to K
              log_t = LOG(t)
              ti = 1.0_wp / t
              Q = T * PERC                             ! perc=0.01
              log_q = LOG(q)
              qi =  1.0_wp / q
              q2 = q**2
              s = MAX(25._wp, psao(jc, 1))           ! minimum salinity 25
              sqrt_s = SQRT(s)
              !
              !
              !*        21.2 CHLORINITY (WOOSTER ET AL., 1969)
              !              ---------------------------------

              !
              !*        21.3 SOLUBILITY OF GASES
              !              -------------------------------------------------
              !
              !      CO2 CARBON DIOXIDE, SEA SURFACE
              !      --------------------------------------
              !         LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
              !

              cek0 = c00*ti + c01 + c02 * LOG_q + s * (c03 + c04*t + c05*q2)
              AK0 = EXP(CEK0)*SMICR
              solco2(jc) = AK0


              !      N2O LAUGHING GAS, SEA SURFACE  (WEISS, 1974)
              !      --------------------------------------
              rs = a1 + a2 * (100._wp *ti) + a3 * LOG(t*0.01_wp) &
                 + s * ( b1 +b2*(t*0.01_wp) + b3*(t*0.01_wp)**2)

              satn2o(jc) = atn2o*EXP(rs)
              !
              !      O2 OXYGEN, SEA SURFACE(EQ. 4, WEISS, 1970)
              !      -----------------------------------------------

              oxy = ox0 + ox1 * qi + ox2 * log_q + ox3 * q &
                  + s * (ox4 + ox5 * q + ox6 * q2)
              ! satoxy is a 3d-field
              satoxy(jc,1) = EXP(OXY)*OXYCO


              !      N2 DINITROGEN, SEA SURFACE
              !      -----------------------------------------------
              !        WEISS, R. F. (1970), DEEP-SEA RESEARCH, VOL. 17, 721-735.

              ani = an0 + an1 * qi + an2 * log_q + an3 * q &
                  + s * (an4 + an5 * q + an6 * q2)
              satn2(jc)  = EXP(ANI)*OXYCO

              !
              !
              !*        21.4 DISSOZIATION CONSTANTS
              !              -----------------------------------------
              !
              !*      PK1 and PK2  OF CARBONIC ACID
              !       -------------------------------

              ck1 = c10  + c11*ti  + c12*log_t -(c13 + c14*ti)*sqrt_s &
                   + c15*s + c16*(s*sqrt_s) + LOG(1.-c17*s)
              ck2 = c20  + c21*ti  + c22*log_t -(c23 + c24*ti)*sqrt_s &
                   + c25*s + c26*(s*sqrt_s) + LOG(1.-c27*s)

              !*      PKB OF BORIC ACID
              !       -------------------------------

              CKB= (cb0 + cb1*sqrt_s + cb2*s + cb3*(s*sqrt_s) + cb4*(s**2)) * ti &
                   + cb5 + cb6*sqrt_s + cb7*s                                &
                   -(cb8+cb9*sqrt_s + cb10*s) * log_t                       &
                   + cb11*sqrt_s*t
              !
              !*      PKW OF WATER (H2O)
              !       ------------------------------------

              CKW = cw0  + cw1*ti + cw2*log_t                                &
                    +SQRT_S*(cw3*ti + cw4 +cw5*log_t) + cw6*s

              !*        21.5 K1, K2 OF CARB. ACID, KB OF BORIC ACID, KW of WATER
              !              ----------------------------------------------------------------
              AK1 = EXP(CK1)
              AK2 = EXP(CK2)
              AKB = EXP(CKB)
              AKW = EXP(CKW)
              !
              !*       21.11 STORE CHEMICAL CONSTANTS AT SEA SURFACE
              !              pressure independent
               aksurf(jc,1) = ak1
               aksurf(jc,2) = ak2
               aksurf(jc,3) = akb
               aksurf(jc,4) = AKW

           ENDIF

     END DO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
  !
  !     -----------------------------------------------------------------
  !*        22. CHEMICAL CONSTANTS - DEEP OCEAN
  !     ----------------------------------------------------------------

  IF ( kldtday == 1 ) THEN
!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(jc,kpke,p,t,ti,log_t,q,log_q,qi,s,sqrt_s,&
!HAMOCC_OMP            oxy,ck1,ck2,ckb,ak1,ak2,akb,akw,cp,tc,tc2,&
!HAMOCC_OMP            log10ksp,aksp0,ckw) HAMOCC_OMP_DEFAULT_SCHEDULE

     !
     !*     22.1 APPROX. SEAWATER PRESSURE AT U-POINT DEPTH (BAR)
     !  ----------------------------------------------------------------


  DO jc = start_idx, end_idx
     kpke=klevs(jc)
      DO k = 1, kpke
        p = 1.025e-1_wp * ptiestu(jc,k)   ! pressure
           IF(kpke.ne.0)THEN
           IF (pddpo(jc, k) > 0.5_wp) THEN           ! wet cell
           !
           !*    22.1 SET ABSOLUTE TEMPERATURE
           ! ----------------------------------------------------------------

           t = ptho(jc,k) + tmelt
           ti = 1.0_wp / t
           log_t = LOG(t)
           q = t * perc
           log_q = LOG(q)
           qi = 1.0_wp / q
           !
           !*    22.2 SET LIMITS FOR SALINITY
           ! ----------------------------------------------------------------
           !   (THIS IS DONE TO AVOID COMPUTATIONAL CRASH AT DRY
           !   POINTS DURING CALCULATION OF CHEMICAL CONSTANTS)

           s = MAX(25._wp, psao(jc, k))
           sqrt_s = SQRT(s)
           !
           !*    22.3 CHLORINITY (WOOSTER ET AL., 1969)
           ! ----------------------------------------------------------------


           !
           !*    22.4 SOLUBILITY OF O2 (EQ. 4, WEISS, 1970)
           ! ----------------------------------------------------------------
           oxy = ox0 + ox1 * qi + ox2 * log_q + ox3 * q &
                + s * (ox4 + ox5 * q + ox6 * q**2)

           satoxy(jc,k) = EXP(oxy) * oxyco

           !
           !*    22.5 PK1, PK2 OF CARBONIC ACID, PKB OF BORIC ACID, PKW OF WATER
           ! ----------------------------------------------------------------
           ck1 = c10  + c11*ti  + c12*log_t -(c13 + c14*ti)*sqrt_s &
                + c15*s + c16*(s*sqrt_s) + LOG(1.-c17*s)
           ck2 = c20  + c21*ti  + c22*log_t -(c23 + c24*ti)*sqrt_s &
                + c25*s + c26*(s*sqrt_s) + LOG(1.-c27*s)

           CKB= (cb0 + cb1*sqrt_s + cb2*s + cb3*(s*sqrt_s) + cb4*(s**2)) * ti &
                   + cb5 + cb6*sqrt_s + cb7*s                                &
                   -(cb8+cb9*sqrt_s + cb10*s) * log_t                       &
                   + cb11*sqrt_s*t

           CKW = cw0  + cw1*ti + cw2*log_t                                &
                    +sqrt_s*(cw3*ti + cw4 +cw5*log_t) + cw6*s

           !
           !*    22.6 K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O)
           ! ----------------------------------------------------------------
           AK1 = EXP(CK1)
           AK2 = EXP(CK2)
           AKB = EXP(CKB)
           AKW = EXP(CKW)

           !
           !       FORMULA FOR CP
           !       ----------------------------------------------

           CP = P / (RGAS*T)
           TC = ptho(jc,k)
           tc2 = tc**2

           !
           !       PRESSURE DEPENDENT KB, K1,K2 AND KW
           !       ----------------------------------------------

           AKB3(JC,K) = AKB * EXP(CP*(DEVKB-DEVKBT*TC-DEVKBT2*tc2))
           AK13(JC,K) = AK1 * EXP(CP*(DEVK1-DEVK1T*TC))
           AK23(JC,K) = AK2 * EXP(CP*(DEVK2-DEVK2T*TC))
           AKW3(JC,K) = AKW * EXP(CP*(DEVKW-DEVKWT*TC-DEVKWT2*tc2))

           !
           !*    22.7 APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
           ! ----------------------------------------------------------------
           !     Mucci (1983)

           log10ksp = akcc1 + akcc2*t + akcc3*ti + akcc4*log_t*invlog10       &
                 + (akcc5 +akcc6*t + akcc7*ti) * sqrt_s                 &
                 + akcc8*s + akcc9*(s*sqrt_s)
           aksp0 = 10**(log10ksp)

           !
           !     PRESSURE DEPENDENT APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE (OR ARAGONITE)
           !     ----------------------------------------------------------------

           AKSP(JC,K) = ARACAL * AKSP0 * EXP(CP*(DEVKS-DEVKST*TC))

          END IF
          END IF
        END DO
  END DO
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

  ENDIF  ! kldtay == 1

END SUBROUTINE CHEMCON
