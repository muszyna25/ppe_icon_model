!>
!!  Description:
!!  This module provides service utilities for meteorological calculations.
!!
!! Routines (module procedure)
!!
!!     - satad_V_3D
!!       Corrects the temperature, the specific humidity and the cloud water
!!      content for condensation/evaporation.
!!
!!    - qsat_rho
!!      Specific humidity at water saturation (with respect to flat surface)
!!      depending on the temperature "temp" and the total density "rhotot")
!!
!!    - dqsatdT_rho
!!       Partial derivative of the specific humidity at water saturation with
!!       respect to the temperature at constant total density.
!!
!!      @author Ulrich Blahak
!!
!!    the following functions should  later be replaced by lookup tables
!!     from mo_convect_tables!!!
!!
!!     - pres_sat_water
!!       Saturation water vapour pressure
!!
!!     - pres_sat_ice
!!       Saturation water vapour pressure
!!
!!     - spec_humi
!!       Specific humidity at saturation pressure
!!

!!
!! @par Revision History
!! first implementation and modification for ICON by Kristina Froehlich, DWD (2010-09-15)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_satad


USE mo_kind,               ONLY: ireals=>wp     , &
                                 iintegers=>i4
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                               rvd_m_o => vtmpc1 , & !! rv/rd-1._wp
                                 o_m_rdv        , & !! 1 - r_d/r_v
                                 rdv            , & !! r_d / r_v
                                 cvd            , & !!
                                 cl    => clw   , & !! specific heat of water
                                 lwd   => alv   , & !! latent heat of vapourization
                                 b3    => tmelt , & !!
                                 tmelt

USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b2i   => c3ies , & !!               -- " --
                                 b4w   => c4les , & !!               -- " --
                                 b4i   => c4ies , & !!               -- " --
                                 b234w => c5les , & !!               -- " --
                                 b234i => c5ies     !!               -- " --

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: satad_v_3D
  PUBLIC  :: satad_v_3D_gpu
  PUBLIC  :: spec_humi
  PUBLIC  :: sat_pres_water
  PUBLIC  :: sat_pres_ice
  PUBLIC  :: dqsatdT
  PUBLIC  :: dqsatdT_ice
  PUBLIC  :: qsat_rho



CONTAINS

SUBROUTINE satad_v_3D (maxiter, tol, te, qve, qce,    & ! IN, INOUT
                       rhotot,                        & ! IN
                       qtvar,                         & ! IN optional
     idim, kdim, ilo, iup, klo, kup,                  & ! IN
                                           errstat)     ! optional

  !-------------------------------------------------------------------------------
  !
  ! Description:
  !   This routine corrects the temperature (te), the specific humidity (qve),
  !   the cloud water content (qce) and the pressure (pe) for condensation/evaporation.
  !   Pressure adapts itself in ICON but has to be rediagnosed in COSMO
  !
  ! Method:
  !   Saturation adjustment at constant total density (adjustment of T and p accordingly)
  !   assuming chemical equilibrium of water and vapor. For the heat capacity of
  !   of the total system (dry air, vapor, and hydrometeors) the value of dry air
  !   is taken, which is a common approximation and introduces only a small error.
  !
  ! Inout fields: te, qve, qce (both COSMO and ICON); ppe (COSMO only)
  !
  ! Input only fields: rhotot (ICON only); qle, qie, p0e (COSMO only)
  !
  ! Outputs:  - count  :  number of iterations needed, maximum of all iterated gridpoints
  !           - errstat:  error status, /= 0 denotes an error; in this case, the calling routine
  !                       has to abort the model run immediately!!!
  !
  ! Description of input/inout - fields:
  !
  ! te:           abs. temperature (adjusted during satad)           [K]
  ! qve:          specific vapor content (adjusted during satad)     [-]
  ! qce:          specific cloud content (adjusted during satad)     [-]
  ! rhotot:       total density, assumed constant during satad       [kg/m^3]   (ICON only)
  ! qle:          non-adjusted liquid water (i.e., rain drops)       [-]        (COSMO only)
  ! qie:          non-adjusted total ice (i.e., sum of qi, qs, ...)  [-]        (COSMO only)
  ! p0e:          base state pressure                                [Pa]       (COSMO only)
  ! ppe:          pressure deviation                                 [Pa]       (COSMO only)
  !
  !-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! Subroutine arguments:
  ! --------------------
  INTEGER (KIND=iintegers), INTENT (IN)    ::  &
       maxiter,                 & !  Max. number of iterations in the numerical scheme
       idim, kdim,              & !  Dimension of I/O-fields
       ilo, iup, klo, kup         !  start- and end-indices for the computations

  REAL    (KIND=ireals),    INTENT (IN) ::  &
       tol                 ! Desired abs. accuracy in K of the adjusted temperature

  REAL    (KIND=ireals),    INTENT (INOUT), DIMENSION(:,:) ::  &  !  dim (idim,kdim)
       te      , & ! Temperature on input/ouput
       qve     , & ! Specific humidity on input/output
       qce         ! Specific cloud water content on input/output

  REAL    (KIND=ireals),    INTENT (IN),  DIMENSION(:,:) ::  &  !  dim (idim,kdim)
       rhotot    ! density containing dry air and water constituents

  REAL    (KIND=ireals),    INTENT (IN), OPTIONAL, DIMENSION(:,:)       ::  &  !  dim (idim,kdim)
       qtvar     ! total water variance - needed only for EDMF

!KF error status temporarly set to optional
  INTEGER (KIND=iintegers), INTENT (OUT),  OPTIONAL  ::  &
       errstat                ! Error status of the saturation adjustment
                                !(errstat /= 0 means an error)
  INTEGER (KIND=iintegers)   ::  &
       count                    ! number of iterations actually needed;
                                ! maximum over all iterated gridpoints
  ! Local variables:
  ! -------------
  INTEGER (KIND=iintegers) ::  &
       i, k                      !  Loop indices
  INTEGER (KIND=iintegers) ::  &
       nsat,                   & !  Number of saturated gridpoints
       iwrk(idim*kdim),        & !  i-index of saturated gridpoints
       kwrk(idim*kdim),        & !  k-index of saturated gridpoints
       indx                        !, & !  loop index

  REAL    (KIND=ireals   ) ::  &
       Ttest(idim,kdim), qtest(idim,kdim), qw(idim,kdim), qwd, qwa, dqwd, fT, dfT, & !, cvvmcl, qd,
       zqwmin              ! Minimum cloud water content for adjustment

  REAL    (KIND=ireals),  DIMENSION(idim,kdim) ::  &
       lwdocvd  ! (Temperature-dependent) latent heat of vaporization over cv

  REAL (KIND=ireals), DIMENSION(idim*kdim) ::  &
       twork, tworkold


  REAL (KIND=ireals), PARAMETER :: cp_v = 1850._ireals ! specific heat of water vapor J
                                                       !at constant pressure
                                                       ! (Landolt-Bornstein)
  LOGICAL :: ll_satad, lqtvar


  !------------ End of header ----------------------------------------------------

  !-------------------------------------------------------------------------------
  ! Begin Subroutine satad
  !-------------------------------------------------------------------------------

  ! Initialization

  IF (PRESENT(errstat)) errstat = 0_iintegers

  IF (PRESENT(qtvar)) THEN
    lqtvar = .TRUE.
  ELSE
    lqtvar = .FALSE.
  ENDIF

  zqwmin = 1.0E-20_ireals

  ! Counter for the number of gridpoints which need the Newton iteration:
  nsat = 0

  ! KF old uncorrect version
  ! Some constants:
  ! lwdocvd = lwd / cvd

!!!=============================================================================================


    DO k = klo, kup
      DO i = ilo , iup

        ! total content of the species which are changed by the adjustment:
        qw(i,k) = qve(i,k) + qce(i,k)

        ! check, which points will still be subsaturated even
        ! if all the cloud water would have been evaporated.
        ! At such points, the Newton iteration is not necessary and the
        ! adjusted values of T, p, qv and qc can be obtained directly.

        lwdocvd(i,k) = ( lwd + (cp_v - cl)*(te(i,k)-tmelt) - r_v*te(i,k) )/ cvd

        Ttest(i,k) = te(i,k) - lwdocvd(i,k)*qce(i,k)

        qtest(i,k) = qsat_rho(Ttest(i,k), rhotot(i,k))

      END DO
    END DO

    DO k = klo, kup
      DO i = ilo , iup

       ll_satad = .TRUE.
       IF (lqtvar) THEN
         IF ( qtvar(i,k) > 0.0001_ireals * (qve(i,k)+qce(i,k))**2 ) THEN ! satad not in EDMF boundary layer
           ll_satad = .false.                                            ! sqrt(qtvar) > 0.001
         ENDIF
       ENDIF

       IF ( ll_satad ) THEN

        IF (qw(i,k) <= qtest(i,k) ) THEN
          ! In this case, all the cloud water evaporates and there is still (sub)saturation.
          ! The resulting state depends only on the available cloud water and is
          ! not saturated, which enables direct computation of the adjusted variables:
          qve(i,k)  = qw(i,k)
          qce(i,k)  = 0.0_ireals
          te(i,k)   = Ttest(i,k)
        ELSE
          ! In this case, the Newton interation is needed
          nsat       = nsat+1
          iwrk(nsat) = i
          kwrk(nsat) = k

          ! Field for the iterated temperature, here set the starting value for the below iteration
          ! to the "old" temperature (as an alternative, the arithmetic mean
          ! between "old" temperature and dew point has been tested, but did
          ! not significantly increase the convergence speed of the iteration):
          twork(nsat)  = te(i,k)

          ! And this is the storage variable for the "old" values in the below iteration:
          ! Add some nonesense increment to the starting, which is sufficient to trigger the
          ! iteration below:
          tworkold(nsat) = twork(nsat) + 10.0_ireals

        END IF
        
       END IF

      END DO
    END DO

    ! Do the Newton iteration at gridpoints which need it:
    ! ---------------------------------------------------

!!!=============================================================================
!!! ..This is the version for the NEC. On a scalar system, it might be better
!!!   to do the iteration within the "indx"-loop, not outside!
!!!=============================================================================


    IF (nsat.gt.0) THEN

      count = 0
      DO WHILE (ANY(ABS(twork(1:nsat)-tworkold(1:nsat)) > tol) .AND. count < maxiter)
        DO indx = 1, nsat
          ! The following if-clause is necessary to achieve reproducible results.
          ! If it is not applied, all grid points are iterated until the "worst"
          ! gridpoint has reached convergence, so the result on one gridpoint
          ! depends on some other gridpoint on the same processor grid.
          IF (ABS(twork(indx)-tworkold(indx)) > tol) THEN
            ! Here we still have to iterate ...
            i = iwrk(indx)
            k = kwrk(indx)
            tworkold(indx) = twork(indx)
            qwd  = qsat_rho(twork(indx), rhotot(i,k))
            dqwd = dqsatdT_rho(qwd, twork(indx))
            ! Newton:
            fT = twork(indx) - te(i,k) + lwdocvd(i,k)*(qwd - qve(i,k))
            dfT = 1.0_ireals + lwdocvd(i,k)*dqwd
            twork(indx) = twork(indx) - fT / dfT;
          END IF
        END DO
        count = count + 1
      END DO

!      IF (ANY(ABS(twork(1:nsat)-tworkold(1:nsat)) > tol)) THEN
!        WRITE(*,*) 'SATAD_V_3D: convergence error of Newton iteration'
!        errstat = 2
!        RETURN
!      END IF

      ! Distribute the results back to gridpoints:
      ! ------------------------------------------

!$NEC ivdep
!DIR$ IVDEP
      DO indx = 1, nsat
        i = iwrk(indx)
        k = kwrk(indx)
        te (i,k) = twork(indx)

        ! We disregard here the extrapolation of qsat from the second-last iteration
        ! step, which is done in the original routine to exactly preserve the internal energy.
        ! This introduces a small error (see the COSMO-Documentation, Part III).
        qwa = qsat_rho(te(i,k), rhotot(i,k))
        qce(i,k) = MAX( qce(i,k) + qve(i,k) - qwa,zqwmin)
        qve(i,k) = qwa
      ENDDO

    END IF

    ! ending the pcreate

!!!=============================================================================================

!!!=============================================================================================

!  END SELECT

!KF
! since ICON re-diagnoses the Exner pressure after satad we do not need
! the following:

END SUBROUTINE satad_v_3D

SUBROUTINE satad_v_3D_gpu (maxiter, tol, te, qve, qce,    & ! IN, INOUT
                       rhotot,                        & ! IN
                       qtvar,                         & ! IN optional
     idim, kdim, ilo, iup, klo, kup,                  & ! IN
                                           errstat)     ! optional

  !-------------------------------------------------------------------------------
  !
  ! Description:
  !   This routine corrects the temperature (te), the specific humidity (qve),
  !   the cloud water content (qce) and the pressure (pe) for condensation/evaporation.
  !   Pressure adapts itself in ICON but has to be rediagnosed in COSMO
  !
  ! Method:
  !   Saturation adjustment at constant total density (adjustment of T and p accordingly)
  !   assuming chemical equilibrium of water and vapor. For the heat capacity of
  !   of the total system (dry air, vapor, and hydrometeors) the value of dry air
  !   is taken, which is a common approximation and introduces only a small error.
  !
  ! Inout fields: te, qve, qce (both COSMO and ICON); ppe (COSMO only)
  !
  ! Input only fields: rhotot (ICON only); qle, qie, p0e (COSMO only)
  !
  ! Outputs:  - count  :  number of iterations needed, maximum of all iterated gridpoints
  !           - errstat:  error status, /= 0 denotes an error; in this case, the calling routine
  !                       has to abort the model run immediately!!!
  !
  ! Description of input/inout - fields:
  !
  ! te:           abs. temperature (adjusted during satad)           [K]
  ! qve:          specific vapor content (adjusted during satad)     [-]
  ! qce:          specific cloud content (adjusted during satad)     [-]
  ! rhotot:       total density, assumed constant during satad       [kg/m^3]   (ICON only)
  ! qle:          non-adjusted liquid water (i.e., rain drops)       [-]        (COSMO only)
  ! qie:          non-adjusted total ice (i.e., sum of qi, qs, ...)  [-]        (COSMO only)
  ! p0e:          base state pressure                                [Pa]       (COSMO only)
  ! ppe:          pressure deviation                                 [Pa]       (COSMO only)
  !
  !-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! Subroutine arguments:
  ! --------------------
  INTEGER (KIND=iintegers), INTENT (IN)    ::  &
       maxiter,                 & !  Max. number of iterations in the numerical scheme
       idim, kdim,              & !  Dimension of I/O-fields
       ilo, iup, klo, kup         !  start- and end-indices for the computations

  REAL    (KIND=ireals),    INTENT (IN) ::  &
       tol                 ! Desired abs. accuracy in K of the adjusted temperature


  REAL    (KIND=ireals),    INTENT (INOUT), DIMENSION(:,:) ::  &  !  dim (idim,kdim)
       te      , & ! Temperature on input/ouput
       qve     , & ! Specific humidity on input/output
       qce         ! Specific cloud water content on input/output

  REAL    (KIND=ireals),    INTENT (IN),  DIMENSION(:,:) ::  &  !  dim (idim,kdim)
       rhotot    ! density containing dry air and water constituents

  REAL    (KIND=ireals),    INTENT (IN), OPTIONAL, DIMENSION(:,:)       ::  &  !  dim (idim,kdim)
       qtvar     ! total water variance - needed only for EDMF

!KF error status temporarly set to optional
  INTEGER (KIND=iintegers), INTENT (OUT),  OPTIONAL  ::  &
       errstat                ! Error status of the saturation adjustment
                                !(errstat /= 0 means an error)
  INTEGER (KIND=iintegers)   ::  &
       count                    ! number of iterations actually needed;
                                ! maximum over all iterated gridpoints
  ! Local variables:
  ! -------------
  INTEGER (KIND=iintegers) ::  &
       i, k                      !  Loop indices

  REAL    (KIND=ireals   ) ::  &
       Ttest, qtest, qw, qwd, qwa, dqwd, fT, dfT, & !, cvvmcl, qd,
       zqwmin              ! Minimum cloud water content for adjustment

  REAL    (KIND=ireals) ::  &
       lwdocvd  ! (Temperature-dependent) latent heat of vaporization over cv

  REAL (KIND=ireals) ::  &
       twork, tworkold

  LOGICAL :: iter_mask

  REAL (KIND=ireals), PARAMETER :: cp_v = 1850._ireals ! specific heat of water vapor J
                                                       !at constant pressure
                                                       ! (Landolt-Bornstein)
  LOGICAL :: ll_satad, lqtvar


  !------------ End of header ----------------------------------------------------

  !-------------------------------------------------------------------------------
  ! Begin Subroutine satad
  !-------------------------------------------------------------------------------

  ! Initialization

  IF (PRESENT(errstat)) errstat = 0_iintegers

  IF (PRESENT(qtvar)) THEN
    lqtvar = .TRUE.
  ELSE
    lqtvar = .FALSE.
  ENDIF

  zqwmin = 1.0E-20_ireals

  ! KF old uncorrect version
  ! Some constants:
  ! lwdocvd = lwd / cvd

!!!=============================================================================================

  !$ACC PARALLEL DEFAULT(PRESENT)
  !$ACC LOOP GANG VECTOR TILE(128,*)
  DO k = klo, kup
    DO i = ilo , iup
      ! total content of the species which are changed by the adjustment:
      qw = qve(i,k) + qce(i,k)

      ! check, which points will still be subsaturated even
      ! if all the cloud water would have been evaporated.
      ! At such points, the Newton iteration is not necessary and the
      ! adjusted values of T, p, qv and qc can be obtained directly.

      lwdocvd = ( lwd + (cp_v - cl)*(te(i,k)-tmelt) - r_v*te(i,k) )/ cvd
      Ttest = te(i,k) - lwdocvd*qce(i,k)
      qtest = qsat_rho(Ttest, rhotot(i,k))

      ll_satad = .TRUE.
      IF (lqtvar) THEN
        IF ( qtvar(i,k) > 0.0001_ireals * (qve(i,k)+qce(i,k))**2 ) THEN ! satad not in EDMF boundary layer
          ll_satad = .false.                                            ! sqrt(qtvar) > 0.001
        ENDIF
      ENDIF

      iter_mask = .FALSE.
      twork = 0.0_ireals
      IF ( ll_satad ) THEN
        IF (qw <= qtest ) THEN
          ! In this case, all the cloud water evaporates and there is still (sub)saturation.
          ! The resulting state depends only on the available cloud water and is
          ! not saturated, which enables direct computation of the adjusted variables:
          qve(i,k)  = qw
          qce(i,k)  = 0.0_ireals
          te(i,k)   = Ttest
        ELSE
          ! Field for the iterated temperature, here set the starting value for the below iteration
          ! to the "old" temperature (as an alternative, the arithmetic mean
          ! between "old" temperature and dew point has been tested, but did
          ! not significantly increase the convergence speed of the iteration):
          twork  = te(i,k)
          ! And this is the storage variable for the "old" values in the below iteration:
          ! Add some nonesense increment to the starting, which is sufficient to trigger the
          ! iteration below:
          tworkold = twork + 10.0_ireals
          iter_mask = .TRUE.
        END IF
      END IF

      ! Do the Newton iteration at gridpoints which need it:
      ! ---------------------------------------------------
      !$ACC LOOP SEQ
      DO count = 1, maxiter
        IF (ABS(twork-tworkold) > tol .AND. iter_mask) THEN
          ! Here we still have to iterate ...
          tworkold = twork
          qwd  = qsat_rho(twork, rhotot(i,k))
          dqwd = dqsatdT_rho(qwd, twork)
          ! Newton:
          fT = twork - te(i,k) + lwdocvd*(qwd - qve(i,k))
          dfT = 1.0_ireals + lwdocvd*dqwd
          twork = twork - fT / dfT;
        END IF
      END DO !while

      ! Distribute the results back to gridpoints:
      ! ------------------------------------------
      IF(iter_mask) THEN
        te (i,k) = twork
        ! We disregard here the extrapolation of qsat from the second-last iteration
        ! step, which is done in the original routine to exactly preserve the internal energy.
        ! This introduces a small error (see the COSMO-Documentation, Part III).
        qwa = qsat_rho(te(i,k), rhotot(i,k))
        qce(i,k) = MAX( qce(i,k) + qve(i,k) - qwa,zqwmin)
        qve(i,k) = qwa
      ENDIF

    ENDDO !i
  ENDDO !k
  !$ACC END PARALLEL

!!!=============================================================================================
!!!=============================================================================================

END SUBROUTINE satad_v_3D_gpu

ELEMENTAL FUNCTION sat_pres_water(temp)
IMPLICIT NONE

REAL (KIND=ireals)              :: sat_pres_water
REAL (KIND=ireals), INTENT(IN)  :: temp

!$ACC ROUTINE SEQ

sat_pres_water = b1*EXP( b2w*(temp-b3)/(temp-b4w) )

END FUNCTION sat_pres_water

!==============================================================================

ELEMENTAL FUNCTION sat_pres_ice(temp)
IMPLICIT NONE

REAL (KIND=ireals)              :: sat_pres_ice
REAL (KIND=ireals), INTENT(IN)  :: temp

!$ACC ROUTINE SEQ

sat_pres_ice = b1*EXP( b2i*(temp-b3)/(temp-b4i) )

END FUNCTION sat_pres_ice

!==============================================================================

ELEMENTAL FUNCTION spec_humi(pvap,pres)
IMPLICIT NONE

REAL (KIND=ireals)              :: spec_humi
REAL (KIND=ireals), INTENT(IN)  :: pres,pvap

!$ACC ROUTINE SEQ

spec_humi = rdv*pvap/( pres - o_m_rdv*pvap )

END FUNCTION spec_humi

! Routine with ACC ROUTINE SEQ cannot be elemental
#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION qsat_rho(temp, rhotot)
  !$ACC ROUTINE SEQ

    !-------------------------------------------------------------------------------
    !
    ! Description:
    !   Specific humidity at water saturation (with respect to flat surface)
    !   depending on the temperature "temp" and the total density "rhotot")
    !
    !-------------------------------------------------------------------------------

  ! input variables: temperature [K], total density [kg/m^2]:
  REAL (KIND=ireals)            :: qsat_rho
  REAL(kind=ireals), INTENT(IN) :: temp, rhotot

  qsat_rho   = sat_pres_water(temp) / (rhotot * r_v * temp)

END FUNCTION qsat_rho

! Routine with ACC ROUTINE SEQ cannot be elemental
#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION dqsatdT_rho(zqsat, temp)
  !$ACC ROUTINE SEQ

    !-------------------------------------------------------------------------------
    !
    ! Description:
    !   Partial derivative of the specific humidity at water saturation with
    !   respect to the temperature at constant total density.
    !   Depends on the temperature "temp" and the
    !   saturation specific humidity "zqsat".
    !
    !-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! input variables: specific saturation humidity [-], temperature [K]:
  REAL (KIND=ireals)            :: dqsatdT_rho
  REAL(kind=ireals), INTENT(in) :: zqsat, temp

  ! local variables:
  REAL(kind=ireals) :: beta

  beta        = b234w/(temp-b4w)**2_iintegers - 1.0_ireals / temp
  dqsatdT_rho = beta * zqsat

END FUNCTION dqsatdT_rho

! UB_20100525<<

ELEMENTAL  FUNCTION dqsatdT (zqsat, temp)

    !-------------------------------------------------------------------------------
    !>
    !! Description:
    !!   Partial derivative of the specific humidity at water saturation with
    !!   respect to the temperature
    !-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! input variables: specific saturation humidity [-], temperature [K]:
  REAL (KIND=ireals)            :: dqsatdT
  REAL(kind=ireals), INTENT(in) :: zqsat, temp

  dqsatdT = b234w * ( 1.0_ireals + rvd_m_o*zqsat ) * zqsat &
                             / (temp-b4w)**2

END FUNCTION dqsatdT

ELEMENTAL  FUNCTION dqsatdT_ice (zqsat, temp)

    !-------------------------------------------------------------------------------
    !>
    !! Description:
    !!   Partial derivative of the specific humidity at ice saturation with
    !!   respect to the temperature
    !-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! input variables: specific saturation humidity [-], temperature [K]:
  REAL (KIND=ireals)            :: dqsatdT_ice
  REAL(kind=ireals), INTENT(in) :: zqsat, temp
!$ACC ROUTINE SEQ

  dqsatdT_ice = b234i * ( 1.0_ireals + rvd_m_o*zqsat ) * zqsat &
                             / (temp-b4i)**2

END FUNCTION dqsatdT_ice


END MODULE mo_satad

