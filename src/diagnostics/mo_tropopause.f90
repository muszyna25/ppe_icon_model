!>
!! This module contains parameters and routines needed for the
!! WMO defined tropopause height.
!!
!! Hauke: Is this routine HAMMONIA save?
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
#endif
MODULE mo_tropopause

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: rd, cpd, g=>grav
  USE mo_echam_wmo_config,   ONLY: echam_wmo_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: WMO_tropopause

CONTAINS
  !>
  !! WMO_tropopause - calculation of the tropopause height.
  !!
  !! @par Revision History
  !!  programmed by Dieter Nodorp                 V1.0  March 95
  !!  modified by Thomas Reichler                 V1.2  May 95
  !!  Rewritten for use in ECHAM  Christine Land  V2.0  Jan 96
  !!  Rewritten for use in ECHAM5 Monika Esch           Sep 2003
  !!  Debugged                    Luis Kornblueh        Nov 2003
  !!
  !! @par Purpose
  !!  WMO_tropopause computes the tropopause height following the
  !!  definition of the height of the tropopause as postulated
  !!  by the WMO (1957).
  !!
  !! @par Method
  !!   - p**k interpolation for the temperature Tm and
  !!     the temperature gradient Gamma of each model layer
  !!   - p**kappa interpolation for p-WMO
  !!
  !! @par References
  !!
  !!     WMO (1992): International meteorological vocabulary, Genf, 784pp.
  !!
  !!   " 1. The first tropopause is defined as the lowest level at which
  !!        the lapse rate decreases to 2 deg C per kilometer or less,
  !!        provided also the average lapse rate between this level and
  !!        all higher levels within 2 kilometers does not exceed 2 deg C
  !!        per kilometer."
  !!
  !!     Reichler (1995): Eine globale Klimatologie der Tropopausenhoehe
  !!                 auf der Basis von ECMWF-Analysen, Diplomarbeit
  !!                 Universitaet Augsburg, 145pp.
  !!
  !!     Dameris, M., D. Nodorp and R.Sausen, 1995: Correlation between
  !!                 Tropopause Height Pressure and TOMS-Data for the
  !!                 EASOE-Winter 1991/1992, Beitr. Phys. Atmosph., 68,
  !!                 227-232.
  !!
  SUBROUTINE WMO_tropopause( jg,                         &
                             jcs, kproma, kbdim, klev,   &
                             ptm1, papm1,                &
    ! Routines of ART use this subroutine also in combination with NWP physics.
    ! Therefore, the parameters iplimb and iplimt must not be limited to
    ! ECHAM physics in this case.
                             ptropo, iplimb, iplimt      )

    ! scalar arguments
    INTEGER, INTENT(in) :: jg
    INTEGER, INTENT(in) :: jcs, kproma, kbdim, klev

    INTEGER, INTENT(inout), OPTIONAL :: iplimb     ! bottom level to search for the tropopause
    INTEGER, INTENT(inout), OPTIONAL :: iplimt     !   top  level to search for the tropopause

    ! array arguments
    REAL(wp), INTENT(in)    :: ptm1(kbdim,klev), papm1(kbdim,klev)
    REAL(wp), INTENT(inout) :: ptropo(kbdim)

    ! local arrays
    REAL(wp) :: ztropo(kproma)

    REAL(wp) :: zpmk(kproma,klev), zpm(kproma,klev),   za(kproma,klev)
    REAL(wp) :: zb(kproma,klev),   ztm(kproma,klev),   zdtdz(kproma,klev)
    REAL(wp) :: zplimb(kproma), zplimt(kproma)

    REAL(wp) :: zptph, zp2km, zkappa, zzkap, zfaktor, zgwmo, zdeltaz
    REAL(wp) :: zag, zbg

    INTEGER :: jl, jk, jj

    INTEGER :: kcount

    REAL(wp) :: zasum, zamean
    REAL(wp) :: zpapm1(kproma,klev)

    zkappa   = rd/cpd                   ! 0.286
    zzkap    = 1.0_wp/zkappa
    zfaktor  = (-1.0_wp)*g/rd           ! -9.81/287.0
    zgwmo    = -0.002_wp
    zdeltaz  = 2000.0_wp

    IF (.NOT. PRESENT(iplimb)) THEN
      iplimb   = echam_wmo_config(jg)%jkewmo
    END IF

    IF (.NOT. PRESENT(iplimt)) THEN
      iplimt   = echam_wmo_config(jg)%jkswmo+2
    END IF

    !$ACC DATA PRESENT( ptm1, papm1, ptropo ) &
    !$ACC       CREATE( ztropo, zpmk, zpm, za, zb, ztm, zdtdz, zplimb, zplimt, zpapm1 )

    ! Calculate the height of the tropopause

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs, kproma
      ztropo(jl) = -999.0_wp
      zplimb(jl) = papm1(jl,iplimb)
      zplimt(jl) = papm1(jl,iplimt)
    ENDDO
    !$ACC END PARALLEL

    ! compute dt/dz

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jk = iplimt-2, iplimb+1
      !$ACC LOOP GANG VECTOR
      DO jl = jcs, kproma
        zpapm1(jl,jk)=papm1(jl,jk)**zkappa
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jk = iplimt-1, iplimb+1
      !$ACC LOOP GANG VECTOR
      DO jl = jcs, kproma

        ! ztm   lineare Interpolation in p**kappa
        ! gamma         dt/dp = a * kappa + papm1(jl,jk)**(kappa-1.)

        zpmk(jl,jk) = 0.5_wp*(zpapm1(jl,jk-1)+zpapm1(jl,jk))
        zpm(jl,jk)  = zpmk(jl,jk)**zzkap                   ! p centre
        za(jl,jk)   = (ptm1(jl,jk-1)-ptm1(jl,jk)) &
             /(zpapm1(jl,jk-1)-zpapm1(jl,jk))
        zb(jl,jk)   = ptm1(jl,jk)-(za(jl,jk)*zpapm1(jl,jk))
        ztm(jl,jk)  = za(jl,jk)*zpmk(jl,jk)+zb(jl,jk)      ! T centre

        zdtdz(jl,jk)=zfaktor*zkappa*za(jl,jk)*zpmk(jl,jk)/ztm(jl,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL


    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE( zag, zbg, zptph, zp2km, zasum, kcount, zamean )
    nproma_loop: DO jl = jcs, kproma
      !$ACC LOOP SEQ
      vertical_loop: DO jk = iplimb+1, iplimt-1, -1
        ! First test: valid dt/dz ?
        IF (zdtdz(jl,jk) >  zgwmo .AND.  &     ! dt/dz > -2K/km
              zpm(jl,jk) <= zplimb(jl)) THEN   ! zpm not too low

          ! dtdz is valid - calculatye p_WMO by interpolation between
          ! actual and layer above linear in p^kappa (Dieters new method)

          zag = (zdtdz(jl,jk)-zdtdz(jl,jk+1))/(zpmk(jl,jk)-zpmk(jl,jk+1))
          zbg = zdtdz(jl,jk+1)-zag*zpmk(jl,jk+1)
          IF ((zgwmo-zbg)/zag >= 0.0_wp) THEN
            zptph = ABS((zgwmo-zbg)/zag)**zzkap
          ELSE
            zptph = 0.0_wp
          END IF
          IF (zdtdz(jl,jk+1) >= zgwmo) zptph = zpm(jl,jk)
          zp2km   = zptph+zdeltaz*zfaktor*zpm(jl,jk)/ztm(jl,jk)   ! p at ptph + 2km
          zasum   = 0.0_wp                                        ! zdtdz above
          kcount  = 0                                             ! number of levels above
          ! 2nd test: dt/dz above 2 km must not be lower than -2 K/km
          IF (zptph >= zplimt(jl)) THEN
            !$ACC LOOP SEQ
            vertical_sub_loop: DO jj = jk, iplimt-1, -1

              IF (zpm(jl,jj) <= zptph .AND. zpm(jl,jj) >= zp2km )THEN
                zasum = zasum+zdtdz(jl,jj)
                kcount= kcount+1
                zamean = zasum/REAL(kcount,wp)            ! dt/dz mean
                IF (zamean <= zgwmo) CYCLE vertical_loop  ! dt/dz above < 2K/km
              ELSE IF (zpm(jl,jj) < zp2km) THEN           ! ptropo valid
                ztropo(jl) = zptph
                EXIT vertical_loop
              ENDIF

            END DO vertical_sub_loop

          ENDIF
        ENDIF
      END DO vertical_loop
    END DO nproma_loop
    !$ACC END PARALLEL

    ! if tropopause not found use previous value in time

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs,kproma
      IF (ztropo(jl) > 0.0_wp) THEN
        ptropo(jl) = ztropo(jl)
      END IF
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END SUBROUTINE WMO_tropopause

END MODULE mo_tropopause
