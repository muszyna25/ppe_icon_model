!===============================================================================!
!
! This was originally part of mo_2mom_mcrph_driver, but ..
!
!  Work-around for Intel 14.0.3 optimizer bug.
!  Host-associated variables are incorrectly propagated at -O2.
!  (SVN Comment by Thomas Jahns in icon-hdcp2-20150604, rev 22867)
!
!===============================================================================!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!===============================================================================!

MODULE mo_2mom_prepare

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message, message_text
  USE mo_2mom_mcrph_main, ONLY: particle, particle_lwf, atmosphere
  IMPLICIT NONE
  PUBLIC :: prepare_twomoment, post_twomoment

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_prepare'

CONTAINS

  SUBROUTINE prepare_twomoment(atmo, cloud, rain, ice, snow, graupel, hail, &
       rho, rhocorr, rhocld, pres, w, tk, hhl, &
       nccn, ninpot, ninact, &
       qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, qgl, qhl, &
       lprogccn, lprogin, lprogmelt, its, ite, kts, kte)

    TYPE(atmosphere), INTENT(inout)   :: atmo
    CLASS(particle),  INTENT(inout)   :: cloud, rain, ice, snow
    CLASS(particle),  INTENT(inout)   :: graupel, hail
    REAL(wp), TARGET, DIMENSION(:, :), INTENT(in) :: &
         rho, rhocorr, rhocld, pres, w, tk, hhl
    REAL(wp), DIMENSION(:,:), INTENT(inout) , TARGET :: &
         &               qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh
    LOGICAL, INTENT(in) :: lprogccn, lprogin, lprogmelt
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               nccn, ninpot, ninact
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               qgl,qhl
    INTEGER, INTENT(in) :: its, ite, kts, kte
    INTEGER :: ii, kk

    ! ... Transformation of microphysics variables to densities
    DO kk = kts, kte
      DO ii = its, ite

        ! ... concentrations --> number densities
        qnc(ii,kk) = rho(ii,kk) * qnc(ii,kk)
        qnr(ii,kk) = rho(ii,kk) * qnr(ii,kk)
        qni(ii,kk) = rho(ii,kk) * qni(ii,kk)
        qns(ii,kk) = rho(ii,kk) * qns(ii,kk)
        qng(ii,kk) = rho(ii,kk) * qng(ii,kk)
        qnh(ii,kk) = rho(ii,kk) * qnh(ii,kk)

        ! ... mixing ratios -> mass densities
        qv(ii,kk) = rho(ii,kk) * qv(ii,kk)
        qc(ii,kk) = rho(ii,kk) * qc(ii,kk)
        qr(ii,kk) = rho(ii,kk) * qr(ii,kk)
        qi(ii,kk) = rho(ii,kk) * qi(ii,kk)
        qs(ii,kk) = rho(ii,kk) * qs(ii,kk)
        qg(ii,kk) = rho(ii,kk) * qg(ii,kk)
        qh(ii,kk) = rho(ii,kk) * qh(ii,kk)

        ninact(ii,kk)  = rho(ii,kk) * ninact(ii,kk)

        IF (lprogccn) THEN
          nccn(ii,kk) = rho(ii,kk) * nccn(ii,kk)
        end if
        if (lprogin) then
          ninpot(ii,kk)  = rho(ii,kk) * ninpot(ii,kk)
        end if
        IF (lprogmelt) THEN
          qgl(ii,kk)  = rho(ii,kk) * qgl(ii,kk)
          qhl(ii,kk)  = rho(ii,kk) * qhl(ii,kk)
        END IF

      END DO
    END DO

    IF (lprogmelt.AND.(.not.PRESENT(qgl).or..not.PRESENT(qhl))) THEN
      CALL finish(TRIM(routine),'Error in prepare_twomoment, something wrong with qgl or qhl')
    END IF
        
    ! set pointers
    atmo%w   => w
    atmo%T   => tk
    atmo%p   => pres
    atmo%qv  => qv
    atmo%rho => rho
    atmo%zh  => hhl

    cloud%rho_v   => rhocld
    rain%rho_v    => rhocorr
    ice%rho_v     => rhocorr
    graupel%rho_v => rhocorr
    snow%rho_v    => rhocorr
    hail%rho_v    => rhocorr

    cloud%q   => qc
    cloud%n   => qnc
    rain%q    => qr
    rain%n    => qnr
    ice%q     => qi
    ice%n     => qni
    snow%q    => qs
    snow%n    => qns
    graupel%q => qg
    graupel%n => qng
    hail%q    => qh
    hail%n    => qnh

    SELECT TYPE (graupel)
    CLASS IS (particle_lwf) 
       graupel%l => qgl
    END SELECT

    SELECT TYPE (hail)
    CLASS IS (particle_lwf) 
       hail%l    => qhl
    END SELECT

    ! enforce upper and lower bounds for number concentrations
    ! (may not be necessary or only at initial time)
    DO kk=kts,kte
      DO ii=its,ite
        rain%n(ii,kk) = MIN(rain%n(ii,kk), rain%q(ii,kk)/rain%x_min)
        rain%n(ii,kk) = MAX(rain%n(ii,kk), rain%q(ii,kk)/rain%x_max)
        ice%n(ii,kk) = MIN(ice%n(ii,kk), ice%q(ii,kk)/ice%x_min)
        ice%n(ii,kk) = MAX(ice%n(ii,kk), ice%q(ii,kk)/ice%x_max)
        snow%n(ii,kk) = MIN(snow%n(ii,kk), snow%q(ii,kk)/snow%x_min)
        snow%n(ii,kk) = MAX(snow%n(ii,kk), snow%q(ii,kk)/snow%x_max)
        graupel%n(ii,kk) = MIN(graupel%n(ii,kk), graupel%q(ii,kk)/graupel%x_min)
        graupel%n(ii,kk) = MAX(graupel%n(ii,kk), graupel%q(ii,kk)/graupel%x_max)
        hail%n(ii,kk) = MIN(hail%n(ii,kk), hail%q(ii,kk)/hail%x_min)
        hail%n(ii,kk) = MAX(hail%n(ii,kk), hail%q(ii,kk)/hail%x_max)
      END DO
    END DO

  END SUBROUTINE prepare_twomoment

  SUBROUTINE post_twomoment(atmo, cloud, rain, ice, snow, graupel, hail, &
       rho_r, qnc, nccn, ninpot, ninact, &
       qv, qc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, qgl, qhl,  &
       lprogccn, lprogin, lprogmelt, its, ite, kts, kte)

    TYPE(atmosphere), INTENT(inout)   :: atmo
    CLASS(particle), INTENT(inout)    :: cloud, rain, ice, snow
    CLASS(particle), INTENT(inout)    :: graupel, hail
    REAL(wp), INTENT(in) :: rho_r(:, :)
    REAL(wp), DIMENSION(:,:), INTENT(inout) :: &
         &           qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &           nccn, ninpot, ninact
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               qgl,qhl
    LOGICAL, INTENT(in) :: lprogccn, lprogin, lprogmelt
    INTEGER, INTENT(in) :: its, ite, kts, kte
    INTEGER :: ii, kk
    REAL(wp) :: hlp

    IF (lprogmelt.AND.(.not.PRESENT(qgl).or..not.PRESENT(qhl))) THEN
      CALL finish(TRIM(routine),'Error in post_twomoment, something wrong with qgl or qhl')
    END IF

    ! nullify pointers
    atmo%w   => null()
    atmo%T   => null()
    atmo%p   => null()
    atmo%qv  => null()
    atmo%rho => null()
    atmo%zh  => null()

    cloud%rho_v   => null()
    rain%rho_v    => null()
    ice%rho_v     => null()
    graupel%rho_v => null()
    snow%rho_v    => null()
    hail%rho_v    => null()

    cloud%q   => null()
    cloud%n   => null()
    rain%q    => null()
    rain%n    => null()
    ice%q     => null()
    ice%n     => null()
    snow%q    => null()
    snow%n    => null()
    graupel%q => null()
    graupel%n => null()
    hail%q    => null()
    hail%n    => null()

    SELECT TYPE (graupel)
    CLASS IS (particle_lwf) 
       graupel%l => null()
    END SELECT

    SELECT TYPE (hail)
    CLASS IS (particle_lwf) 
       hail%l    => null()
    END SELECT

    ! ... Transformation of variables back to ICON standard variables
    DO kk = kts, kte
      DO ii = its, ite

        hlp = rho_r(ii,kk)

        ! ... from mass densities back to mixing ratios
        qv(ii,kk) = hlp * qv(ii,kk)
        qc(ii,kk) = hlp * qc(ii,kk)
        qr(ii,kk) = hlp * qr(ii,kk)
        qi(ii,kk) = hlp * qi(ii,kk)
        qs(ii,kk) = hlp * qs(ii,kk)
        qg(ii,kk) = hlp * qg(ii,kk)
        qh(ii,kk) = hlp * qh(ii,kk)

        ! ... number concentrations
        qnc(ii,kk) = hlp * qnc(ii,kk)
        qnr(ii,kk) = hlp * qnr(ii,kk)
        qni(ii,kk) = hlp * qni(ii,kk)
        qns(ii,kk) = hlp * qns(ii,kk)
        qng(ii,kk) = hlp * qng(ii,kk)
        qnh(ii,kk) = hlp * qnh(ii,kk)

        ninact(ii,kk)  = hlp * ninact(ii,kk)

        if (lprogccn) THEN
          nccn(ii,kk) = hlp * nccn(ii,kk)
        end if
        if (lprogin) THEN
          ninpot(ii,kk)  = hlp * ninpot(ii,kk)
        end if
        IF (lprogmelt) THEN
          qgl(ii,kk)  = hlp * qgl(ii,kk)
          qhl(ii,kk)  = hlp * qhl(ii,kk)
        END IF

      ENDDO
    ENDDO

  END SUBROUTINE post_twomoment

END MODULE mo_2mom_prepare
