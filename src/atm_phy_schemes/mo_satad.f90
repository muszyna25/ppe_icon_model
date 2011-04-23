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
                                 b234w => c5les     !!               -- " --

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: satad_v_3D
  PUBLIC  :: spec_humi
  PUBLIC  :: sat_pres_water
  PUBLIC  :: sat_pres_ice
  PUBLIC  :: dqsatdT
  PUBLIC  :: qsat_rho

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

SUBROUTINE satad_v_3D (maxiter, tol, te, qve, qce,    & ! IN, INOUT
#ifdef __ICON__
                       rhotot,                        & ! IN
#endif
#ifdef __COSMO__
                       qle, qie, p0e, ppe,            & ! IN
#endif
     idim, jdim, kdim, ilo, iup, jlo, jup, klo, kup,  & ! IN
!      count,
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
       idim, jdim, kdim,        & !  Dimension of I/O-fields
       ilo, iup, jlo, jup, klo, kup !  start- and end-indices for the computations

  REAL    (KIND=ireals),    INTENT (IN) ::  &
       tol                 ! Desired abs. accuracy in K of the adjusted temperature

#ifdef __COSMO__
  REAL    (KIND=ireals),    INTENT (IN),  DIMENSION(idim,jdim,kdim) ::  &
       p0e      , & ! Reference pressure
       qle      , & ! Liquid hydrometeors which are not affected by saturation adj.
                    ! (i.e., rain drops)
       qie          ! Solid hydrometeors (sum of qi, qs, qg and, if present, qh)
  REAL    (KIND=ireals),    INTENT (INOUT), DIMENSION(idim,jdim,kdim) ::  &
       ppe          ! Pressure deviation from reference pressure needed in COSMO
#endif

  REAL    (KIND=ireals),    INTENT (INOUT), DIMENSION(idim,jdim,kdim) ::  &
       te      , & ! Temperature on input/ouput
       qve     , & ! Specific humidity on input/output
       qce         ! Specific cloud water content on input/output

#ifdef __ICON__
  REAL    (KIND=ireals),    INTENT (IN),  DIMENSION(idim,jdim,kdim) ::  &
       rhotot    ! density containing dry air and water constituents
#endif

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
       i, j, k,                  & !  Loop indices
       nsat,                     & !  Number of saturated gridpoints
       iwrk(idim*jdim*kdim),     & !  i-index of saturated gridpoints
       jwrk(idim*jdim*kdim),     & !  j-index of saturated gridpoints
       kwrk(idim*jdim*kdim),     & !  k-index of saturated gridpoints
       indx                        !, & !  loop index

  REAL    (KIND=ireals   ) ::  &
       Ttest, qtest, qw, qwd, dqwd, fT, dfT, & !, cvvmcl, qd,
       zqwmin              ! Minimum cloud water content for adjustment

#ifdef __COSMO__
  REAL    (KIND=ireals),  DIMENSION(idim,jdim,kdim) ::  &
       rhotot              ! Total density
#endif

  REAL    (KIND=ireals),  DIMENSION(idim,jdim,kdim) ::  &
       lwdocvd  ! (Temperature-dependent) latent heat of vaporization over cv

  REAL (KIND=ireals), DIMENSION(idim*jdim*kdim) ::  &
       twork, tworkold

  REAL (KIND=ireals), PARAMETER :: cp_v = 1850._ireals ! specific heat of water vapor J
                                                       !at constant pressure
                                                       ! (Landolt-Bornstein)


  !------------ End of header ----------------------------------------------------

  !-------------------------------------------------------------------------------
  ! Begin Subroutine satad
  !-------------------------------------------------------------------------------

  ! Initialization

  IF (PRESENT(errstat)) errstat = 0_iintegers

  zqwmin = 1.0E-20_ireals

  ! Counter for the number of gridpoints which need the Newton iteration:
  nsat = 0

  ! KF old uncorrect version
  ! Some constants:
  ! lwdocvd = lwd / cvd

#ifdef __ICON__
  ! j is always 1 for ICON
  j = 1
#endif

!!!=============================================================================================

#ifdef __COSMO__
   ! diagnostically (re)compute density of moist air for time-level nnow
   ! ... using specific moisture quantities
   CALL calrho( te, ppe,qve,qce,qle+qie, p0e, rhotot, idim, jdim, kdim, r_d, rvd_m_o )
#endif

      DO k = klo, kup
#ifndef __ICON__
      DO j = jlo , jup
#endif
        DO i = ilo , iup

          ! total content of the species which are changed by the adjustment:
          qw = qve(i,j,k) + qce(i,j,k)

          ! check, which points will still be subsaturated even
          ! if all the cloud water would have been evaporated.
          ! At such points, the Newton iteration is not necessary and the
          ! adjusted values of T, p, qv and qc can be obtained directly.

          lwdocvd(i,j,k) = ( lwd + (cp_v - cl)*(te(i,j,k)-tmelt) - r_v*te(i,j,k) )/ cvd

          Ttest = te(i,j,k) - lwdocvd(i,j,k)*qce(i,j,k)

          qtest = qsat_rho(Ttest, rhotot(i,j,k))

          IF (qw <= qtest ) THEN
            ! In this case, all the cloud water evaporates and there is still (sub)saturation.
            ! The resulting state depends only on the available cloud water and is
            ! not saturated, which enables direct computation of the adjusted variables:
!            qve(i,j,k)  = qve(i,j,k) + qce(i,j,k)
            qve(i,j,k)  = qw
            qce(i,j,k)  = 0.0_ireals
            te(i,j,k)   = Ttest
          ELSE
            ! In this case, the Newton interation is needed
            nsat       = nsat+1
            iwrk(nsat) = i
            jwrk(nsat) = j
            kwrk(nsat) = k
            ! Field for the iterated temperature, here set the starting value for the below iteration
            ! to the "old" temperature (as an alternative, the arithmetic mean
            ! between "old" temperature and dew point has been tested, but did
            ! not significantly increase the convergence speed of the iteration):
            twork(nsat)  = te(i,j,k)
            ! And this is the storage variable for the "old" values in the below iteration:
            ! Add some nonesense increment to the starting, which is sufficient to trigger the
            ! iteration below:
            tworkold(nsat) = twork(nsat) + 10.0_ireals
          END IF

        END DO
#ifndef __ICON__
      END DO
#endif
    END DO

    ! Do the Newton iteration at gridpoints which need it:
    ! ---------------------------------------------------

!!!=============================================================================
!!! ..This is the version for the NEC. On a scalar system, it might be better
!!!   to do the iteration within the "indx"-loop, not outside!
!!!=============================================================================


    IF (nsat.gt.0) THEN
      count = 0
      DO WHILE (MAXVAL(ABS(twork(1:nsat)-tworkold(1:nsat))) > tol .AND. count < maxiter)
        DO indx = 1, nsat
          ! The following if-clause is necessary to achieve reproducible results.
          ! If it is not applied, all grid points are iterated until the "worst"
          ! gridpoint has reached convergence, so the result on one gridpoint
          ! depends on some other gridpoint on the same processor grid.
          IF (ABS(twork(indx)-tworkold(indx)) > tol) THEN
            ! Here we still have to iterate ...
            i = iwrk(indx)
#ifndef __ICON__
            j = jwrk(indx)
#endif
            k = kwrk(indx)
            tworkold(indx) = twork(indx)
            qwd  = qsat_rho(twork(indx), rhotot(i,j,k))
            dqwd = dqsatdT_rho(qwd, twork(indx))
            ! Newton:
            fT = twork(indx) - te(i,j,k) + lwdocvd(i,j,k)*(qwd - qve(i,j,k))
            dfT = 1.0_ireals + lwdocvd(i,j,k)*dqwd
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

!CDIR NODEP,VOVERTAKE,VOB
      DO indx = 1, nsat
        i = iwrk(indx)
#ifndef __ICON__
        j = jwrk(indx)
#endif
        k = kwrk(indx)
        te (i,j,k) = twork(indx)
        ! We disregard here the extrapolation of qsat from the second-last iteration
        ! step, which is done in the original routine to exactly preserve the internal energy.
        ! This introduces a small error (see the COSMO-Documentation, Part III).
        qw = qsat_rho(te(i,j,k), rhotot(i,j,k))
        qce(i,j,k) = MAX( qce(i,j,k) + qve(i,j,k) - qw,zqwmin)
        qve(i,j,k) = qw
      ENDDO

    END IF

!!!=============================================================================================

!!!=============================================================================================

!  END SELECT

!KF
! since ICON re-diagnoses the Exner pressure after satad we do not need
! the following:
#ifdef __COSMO__
  ! Re-diagnose pressure at all gridpoints:
  ppe(ilo:iup,jlo:jup,:) = rhotot(ilo:iup,jlo:jup,:) * r_d * te(ilo:iup,jlo:jup,:) * &
       ( 1.0_ireals + rvd_m_o*qve(ilo:iup,jlo:jup,:)   &
       - qce(ilo:iup,jlo:jup,:) &
       - qle(ilo:iup,jlo:jup,:) &
       - qie(ilo:iup,jlo:jup,:) )     -    p0e(ilo:iup,jlo:jup,:)
#endif

END SUBROUTINE satad_v_3D

ELEMENTAL FUNCTION sat_pres_water(temp)
IMPLICIT NONE

REAL (KIND=ireals)              :: sat_pres_water
REAL (KIND=ireals), INTENT(IN)  :: temp

sat_pres_water = b1*EXP( b2w*(temp-b3)/(temp-b4w) )

END FUNCTION sat_pres_water

!==============================================================================

ELEMENTAL FUNCTION sat_pres_ice(temp)
IMPLICIT NONE

REAL (KIND=ireals)              :: sat_pres_ice
REAL (KIND=ireals), INTENT(IN)  :: temp

sat_pres_ice = b1*EXP( b2i*(temp-b3)/(temp-b4i) )

END FUNCTION sat_pres_ice

!==============================================================================

ELEMENTAL FUNCTION spec_humi(pvap,pres)
IMPLICIT NONE

REAL (KIND=ireals)              :: spec_humi
REAL (KIND=ireals), INTENT(IN)  :: pres,pvap

spec_humi = rdv*pvap/( pres - o_m_rdv*pvap )

END FUNCTION spec_humi

ELEMENTAL FUNCTION qsat_rho(temp, rhotot)

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

ELEMENTAL  FUNCTION dqsatdT_rho(zqsat, temp)

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

END MODULE mo_satad

