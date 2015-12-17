MODULE mo_2mom_prepare
  USE mo_kind,                 ONLY: wp
  USE mo_2mom_mcrph_types,     ONLY: particle, atmosphere
  IMPLICIT NONE
  PUBLIC :: prepare_twomoment, post_twomoment
CONTAINS
  SUBROUTINE prepare_twomoment(atmo, cloud, rain, ice, snow, graupel, hail, &
       rho, rhocorr, rhocld, pres, w, tk, &
       nccn, ninpot, ninact, &
       qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, &
       lprogccn, lprogin, its, ite, kts, kte)
    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud, rain, ice, snow, graupel, hail
    REAL(wp), TARGET, DIMENSION(:, :), INTENT(in) :: &
         rho, rhocorr, rhocld, pres, w, tk
    REAL(wp), DIMENSION(:,:), INTENT(inout) , TARGET :: &
         &               qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh
    LOGICAL, INTENT(in) :: lprogccn, lprogin
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &               nccn, ninpot, ninact
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
      END DO
    END DO

    ! set pointers
    atmo%w   => w
    atmo%T   => tk
    atmo%p   => pres
    atmo%qv  => qv
    atmo%rho => rho

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

  END SUBROUTINE prepare_twomoment

  SUBROUTINE post_twomoment(atmo, cloud, rain, ice, snow, graupel, hail, &
       rho_r, qnc, nccn, ninpot, ninact, &
       qv, qc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh, &
       lprogccn, lprogin, its, ite, kts, kte)

    TYPE(atmosphere), INTENT(inout) :: atmo
    TYPE(particle), INTENT(inout) :: cloud, rain, ice, snow, graupel, hail
    REAL(wp), INTENT(in) :: rho_r(:, :)
    REAL(wp), DIMENSION(:,:), INTENT(inout) :: &
         &           qv, qc, qnc, qr, qnr, qi, qni, qs, qns, qg, qng, qh, qnh
    REAL(wp), DIMENSION(:,:), INTENT(INOUT), TARGET, OPTIONAL :: &
         &           nccn, ninpot, ninact
    LOGICAL, INTENT(in) :: lprogccn, lprogin
    INTEGER, INTENT(in) :: its, ite, kts, kte
    INTEGER :: ii, kk
    REAL(wp) :: hlp

    ! nullify pointers
    atmo%w   => null()
    atmo%T   => null()
    atmo%p   => null()
    atmo%qv  => null()
    atmo%rho => null()

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
      ENDDO
    ENDDO

  END SUBROUTINE post_twomoment

END MODULE mo_2mom_prepare
