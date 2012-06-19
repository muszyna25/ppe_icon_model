SUBROUTINE update_surfacetemp_icon(klon, time_step_len,              &
  &            cvdifts, pcp,                                         &
  &            pescoe, pfscoe, peqcoe, pfqcoe,                       &
  &            psold, pqsold, pdqsold,                               &
  &            pnetrad, pgrdfl,                                      &
  &            pcfh, pcair, pcsat, pfracsu, pgrdcap,                 &
  &            psnew, pqsnew)


!!$ TR  USE mo_time_control,     ONLY: time_step_len
!!$ TR#ifndef STANDALONE
!!$ TR    USE mo_physc2, ONLY: cvdifts
!!$ TR#endif
!!$ TR  USE mo_radiation_parameters,        ONLY: cemiss
!!$ TR  USE mo_constants,        ONLY: stbo, cpd, vtmpc2, alv, als
  USE mo_kind,             ONLY: wp

  IMPLICIT NONE

  INTEGER,  INTENT(in)    :: klon
  REAL(wp), INTENT(in)    :: time_step_len
  REAL(wp), INTENT(in)    :: cvdifts
  REAL(wp),     INTENT(in)    :: pcp(klon), pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
  REAL(wp),     INTENT(in)    :: psold(klon), pqsold(klon), pdqsold(klon)
  REAL(wp),     INTENT(in)    :: pnetrad(klon), pgrdfl(klon)
  REAL(wp),     INTENT(in)    :: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
  REAL(wp),     INTENT(in)    :: pgrdcap(klon)
  REAL(wp),     INTENT(out)   :: psnew(klon), pqsnew(klon)
!
  INTEGER :: jl
  REAL(wp) :: zcolin(klon), zcohfl(klon), zcoind(klon), zicp(klon), zca(klon), zcs(klon)
  REAL(wp) :: pdt, pemi, pboltz
  REAL(wp) :: pc16, platev, platsu
  REAL(wp) :: ztpfac1, ztmst
  REAL(wp) :: cpd, cpv

!!$ TR#ifdef STANDALONE
!!$ TR    REAL(wp), PARAMETER :: cvdifts = 1.0_wp
!!$ TR#endif


!-------------------------------------------------------------------------------------
! Constants
  cpd = 1005.46_wp                         ! specific heat of dry air at constant pressure in J/K/kg
  cpv = 1869.46_wp                         ! specific heat of water vapour at constant pressure in J/K/kg
  ztpfac1 = cvdifts
  ztmst   = time_step_len
  pdt     = ztpfac1*ztmst                  ! zcons29 in 'old' vdiff
  pemi    = 0.996_wp                       ! emissivity: compare wirh cemiss in echam/mo_radiation_parameters.f90
  pboltz  = 5.67E-8_wp                     ! Stefan Boltzman constant: compare with stbo in echam/mo_constants.f90
  pc16    = cpd * (cpv / cpd - 1.0_wp)     ! cpd * vtmpc2 : compare with echam/mo_constants.f90
  platev  = 2.5008e6_wp                    ! latent heat for vaporisation in J/kg: compare with alv in echam/mo_constants.f90
  platsu  = 2.8345e6_wp                    ! latent heat for sublimation in J/kg: compare with als in echam/mo_constants.f90

!************************************************************************************
!
     zicp(:) = 1._wp / pcp(:)
!
     zca(:)    = platsu * pfracsu(:) +  platev * (pcair(:) - pfracsu(:))
     zcs(:)    = platsu * pfracsu(:) +  platev * (pcsat(:) - pfracsu(:))
!
     zcolin(:) = pgrdcap(:)*zicp(:) +                                             &
          &      pdt * (zicp(:) * 4._wp * pemi * pboltz * ((zicp(:) * psold(:))**3) -     &
          &      pcfh(:) * (zca(:) * peqcoe(:) - zcs(:) -                      &
          &      zicp(:) * pc16 * psold(:) *                        &
          &      (pcair(:) * peqcoe(:) - pcsat(:)))*        &
          &      zicp(:) * pdqsold(:))
!
     zcohfl(:) = -pdt * pcfh(:) * (pescoe(:) - 1._wp)
!
     zcoind(:) = pdt * (pnetrad(:) + pcfh(:) * pfscoe(:) +  pcfh(:)*          &
          &      ((zca(:) * peqcoe(:) - zcs(:)) * pqsold(:) + zca(:) * pfqcoe(:) -   &
          &      zicp(:) * pc16 * psold(:) *                                &
          &      ((pcair(:) * peqcoe(:) - pcsat(:)) * pqsold(:) +    &
          &      pcair(:) * pfqcoe(:))) + pgrdfl(:))
!
    psnew(:)  = (zcolin(:) * psold(:) + zcoind(:)) / (zcolin(:) + zcohfl(:))
    pqsnew(:) = pqsold(:) + zicp(:) * pdqsold(:) * (psnew(:) - psold(:))

END SUBROUTINE update_surfacetemp_icon
