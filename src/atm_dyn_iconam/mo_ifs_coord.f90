!>
!! Contains subroutines for calculating the pressure values and
!! some other auxiliary variables related to the pressure-sigma
!! hybrid vertical coordinate (the eta coordinate).
!!
!! @par Revision History
!!   Original version from ECHAM5.3.01
!!   Adapted for ICOHDC by Hui Wan, 2006-02
!!   Modifications include:
!!   - Calculation of half-level geopotential added to *geopot*
!!   - Calculation of logorithm of half-level pressure added to *auxhyb*
!!   - Subroutine <i>full_level_pressure</i> added
!!   Modifications by Hui Wan (MPI-M, 2010-01-29)
!!   - Renamed subroutines and variables.
!!   - Changed the sequence of items in the argument lists.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ifs_coord


  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: success, max_char_length
  USE mo_physical_constants, ONLY: grav, rcpd, rd

  IMPLICIT NONE

  PRIVATE


  INTEGER :: nlevm1         ! (number of levels)-1.
  INTEGER :: nplev          ! *number of pressure levels.
  INTEGER :: nplvp1         ! *nplev+1.*
  INTEGER :: nplvp2         ! *nplev+2.*
  INTEGER :: nlmsgl         ! *nlev* - (number of sigma levels).
  INTEGER :: nlmslp         ! *nlmsgl+1.*
  INTEGER :: nlmsla         ! *nlmslp,* or 2 if *nlmslp=1.*

  REAL(wp) :: apzero        ! *reference pressure for computation of the
  !                         !  hybrid vertical levels.
  REAL(wp) :: t0icao        ! *surface temperatur of reference atmosphere
  REAL(wp) :: tsticao       ! *stratospheric temperature of reference atmosphere
  REAL(wp) :: rdlnp0i       ! *rd*ln(surface pressure) of reference atmosphere
  REAL(wp) :: alrrdic       ! *lapse-rate parameter of reference atmosphere
  REAL(wp) :: rdlnpti       ! *rd*ln(ptricao)

  REAL(wp), ALLOCATABLE :: ralpha(:) ! rd*alpha at pressure and sigma levels.
  REAL(wp), ALLOCATABLE :: rlnpr(:)  ! rd*ln(p(k+.5)/p(k-.5))
  REAL(wp), ALLOCATABLE :: dela(:)   ! a(k+.5)-a(k-.5).
  REAL(wp), ALLOCATABLE :: delb(:)   ! b(k+.5)-b(k-.5).
  REAL(wp), ALLOCATABLE :: rddelb(:) ! rd*delb.
  REAL(wp), ALLOCATABLE :: cpg(:)    ! a(k+.5)*b(k-.5)-b(k+.5)*a(k-.5).
  REAL(wp), ALLOCATABLE :: delpr(:)  ! p(k+.5)-p(k-.5) of
                                     ! the refrence surface pressure.
  REAL(wp), ALLOCATABLE :: rdelpr(:) ! *reciprocal of *delpr.*
  REAL(wp), ALLOCATABLE :: alpham(:) ! *constant array for use by dyn.
  REAL(wp), ALLOCATABLE :: ardprc(:) ! *constant array for use by dyn.
  REAL(wp), ALLOCATABLE :: ceta(:)   ! *full hybrid vertical levels.
  REAL(wp), ALLOCATABLE :: cetah(:)  ! *half hybrid vertical levels.

  REAL(wp), ALLOCATABLE :: vct_a(:) ! param. A of the vertical coordinte
  REAL(wp), ALLOCATABLE :: vct_b(:) ! param. B of the vertical coordinate
  REAL(wp), ALLOCATABLE :: vct  (:) ! param. A and B of the vertical coordinate

  PUBLIC :: half_level_pressure, full_level_pressure
  PUBLIC :: auxhyb, geopot, alloc_vct, init_vct, vct, vct_a, vct_b

CONTAINS

  !>
  !!
  SUBROUTINE alloc_vct(nlev)

    INTEGER, INTENT(IN) :: nlev
    INTEGER :: ist
    INTEGER :: nlevp1

    CHARACTER(len=max_char_length),PARAMETER :: routine = &
         & 'mo_vertical_coord_table:alloc_vct'

    !-----------------------------------------------------------------------
    !BOC

    nlevp1 = nlev+1

    ALLOCATE (ralpha(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of ralpha failed')
    ENDIF

    ALLOCATE (rlnpr(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rlnpr failed')
    ENDIF

    ALLOCATE (dela(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of dela failed')
    ENDIF

    ALLOCATE (delb(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of delb failed')
    ENDIF

    ALLOCATE (rddelb(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rddelb failed')
    ENDIF

    ALLOCATE (cpg(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of cpg failed')
    ENDIF

    ALLOCATE (delpr(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rdelpr failed')
    ENDIF

    ALLOCATE (rdelpr(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of rdelpr failed')
    ENDIF

    ALLOCATE (alpham(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of alpham failed')
    ENDIF

    ALLOCATE (ardprc(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of ardprc failed')
    ENDIF

    ALLOCATE (ceta(nlev), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of ceta failed')
    ENDIF

    ALLOCATE (cetah(nlevp1), STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), ' allocation of cetah failed')
    ENDIF

  END SUBROUTINE alloc_vct

  !EOC
  !-------------------------------------------------------------------------
  !
  !

  !>
  !!  Initializes constants for vertical coordinate calculations.
  !!
  !!  Method:
  !!    Compute loop indices and surface-pressure independent
  !!    variables associated with the vertical finite-difference scheme.
  !!    Output is in module *mo_hyb*
  !!
  !! @par Revision History
  !!  A. J. Simmons, ECMWF, November 1981, original source
  !!  L. Kornblueh, MPI, May 1998, f90 rewrite
  !!  U. Schulzweida, MPI, May 1998, f90 rewrite
  !!  A. Rhodin, MPI, Jan 1999, subroutine inihyb -> module mo_hyb
  !!  H. Wan, MPI-M, Feb 2006, rewrite: "goto" removed
  !!  H. Wan, MPI-M, Aug 2007, new name: init_hyb_params
  !!  A. Gassmann, MPI-M, (2008-04-23), change apzero to 10^5Pa
  !! @par
  !!  for more details see file AUTHORS
  !!
  SUBROUTINE init_vct(nlev)

    INTEGER, INTENT(IN) :: nlev

    !  Local scalars:
    REAL(wp) :: za, zb, zetam, zetap, zp, zp0icao, zpp, zrd, zs, zsm
    INTEGER  :: ilev, ilevp1, iplev, iplvp1, is, ism, ist, &
      &         jk, jlev, nvclev

    !  Intrinsic functions
    INTRINSIC EXP, LOG

    !-----------------------------------------------------------------------
    !BOC

    !  Executable statements

    !-- 1. Initialize variables
    nvclev = nlev+1
!ag    apzero    = 101325._wp ! changed for NCAR summer colloquium!
    apzero    = 100000._wp
    zrd       = rd
    ralpha(1) = zrd*LOG(2._wp)
    rlnpr(1)  = 2._wp*ralpha(1)
    ilev      = nlev
    ilevp1    = ilev + 1
    nlevm1    = ilev - 1
    iplev     = 0
    iplvp1    = 1
    is        = nvclev + ilevp1
    ism       = is - 1
    zpp       = vct(1)
    zsm       = vct(is)

    t0icao  = 288._wp
    tsticao = 216.5_wp
    zp0icao = 101320._wp
    rdlnp0i = rd*LOG(zp0icao)
    alrrdic = 0.0065_wp/grav
    rdlnpti = rdlnp0i + (LOG(tsticao/t0icao))/alrrdic

    zb      = vct(nvclev+iplvp1+1)

    !-- 2. Calculate pressure-level values

    DO WHILE ( zb == 0._wp )

      iplev  = iplvp1
      iplvp1 = iplev + 1
      IF (iplvp1==ilevp1) EXIT    ! if all levels are pressure levels

      zp            = zpp
      zpp           = vct(iplvp1)
      delpr(iplev)  = zpp - zp
      rdelpr(iplev) = 1._wp/delpr(iplev)

      IF ( iplev>1 ) THEN
        rlnpr(iplev)  = zrd*LOG(zpp/zp)
        ralpha(iplev) = zrd - zp*rlnpr(iplev)/delpr(iplev)
      END IF

      alpham(iplev) = ralpha(iplev)*rcpd
      ardprc(iplev) = rlnpr(iplev)*rdelpr(iplev)*rcpd
      zb            = vct(nvclev+iplvp1+1)

    ENDDO


    IF (iplvp1/=ilevp1) THEN   ! All levels are not pressure-levels

      nplev  = iplev
      nplvp1 = iplvp1
      nplvp2 = iplvp1 + 1


      !-- 3. Calculate sigma-level values

      za = vct(ism-nvclev)

      DO WHILE ( za == 0._wp )

        is  = ism
        ism = is - 1
        ist = is - nvclev
        zs  = zsm
        zsm = vct(is)
        IF (ist==1) THEN
          nlmsgl = 0
          nlmslp = 1
          nlmsla = 2
          EXIT
        ELSE
          rlnpr(ist)  = zrd*LOG(zs/zsm)
          ralpha(ist) = zrd - zsm*rlnpr(ist)/(zs-zsm)
        END IF
        za = vct(ism-nvclev)

      END DO

      IF (za>0._wp) THEN
        nlmsgl = ism - nvclev
        nlmslp = nlmsgl + 1
        nlmsla = nlmslp
      END IF

      !-- 4. Calculate dela, delb, rddelb, cpg, and complete alphdb

      DO jk = 1, nlev
        dela(jk)   = vct(jk+1) - vct(jk)
        delb(jk)   = vct(nvclev+jk+1) - vct(nvclev+jk)
        rddelb(jk) = rd*delb(jk)
        cpg(jk)    = vct(nvclev+jk)*vct(jk+1) - vct(nvclev+jk+1)*vct(jk)
      END DO

      DO jk = nlmslp, nlev
        alpham(jk) = ralpha(jk)*delb(jk)
      END DO

    ENDIF  ! If all levels are not pressure-levels

    !-- 5. Compute full level values of the hybrid coordinate

    zetam    = vct(1)/apzero + vct(nvclev+1)
    cetah(1) = zetam

    DO jlev = 1, nlev
      zetap         = vct(jlev+1)/apzero + vct(nvclev+1+jlev)
      ceta(jlev)    = (zetam+zetap)*.5_wp
      cetah(jlev+1) = zetap
      zetam = zetap
    END DO

  END SUBROUTINE init_vct

!-------------------------------------------------------------------------
  !>
  !! Calculate half-level pressures at all model levels
  !! for a given surface pressure.
  !!
  !! @par Method
  !!  Calculations are performed separately for pressure,
  !!  hybrid and sigma levels.
  !!
  !! @par Arguments
  !!   *ps*        surface pressure.
  !!   *kdimp*     first dimension of 2-d array *ph.*
  !!   *klen*      number of points for which calculation is
  !!               performed.
  !!   *ph*        computed half-level pressures.
  !!
  !! @par Parameters
  !!  Required constants are obtained from module *mo_hyb*.
  !!  The latter must have been initialized by a call of
  !!  subroutine init_vertical_coord.
  !!
  !! @par Results
  !!  Results are computed for *klen* consecutive points at
  !!  each model half level.
  !!
  !! @see
  !!  External documentation of the model equations and the
  !!  organization of the vertical calculation.
  !!
  !! @par Revision History
  !!    A. J. Simmons, ECMWF, November 1981, original source
  !!    L. Kornblueh, MPI, May 1998, f90 rewrite
  !!    U. Schulzweida, MPI, May 1998, f90 rewrite
  !!    H. Wan, MPI, Feb 2006, adapted for ICOHDC
  !!    H. Wan, MPI, Jan 2010, renamed the interface
  !!
  SUBROUTINE half_level_pressure( ps,kdimp,klen,nlev_in, ph)

    INTEGER ,INTENT(in)  :: kdimp
    REAL(wp),INTENT(in)  :: ps(kdimp)   !< surface pressure
    INTEGER ,INTENT(in)  :: klen, nlev_in

    REAL(wp),INTENT(inout) :: ph(kdimp,nlev_in+1) !< half-level pressure

    REAL(wp):: zb, zp
    INTEGER :: jk, jl, nvclev, nlevp1

    nvclev = nlev_in+1
    nlevp1 = nvclev

    ! Transfer pressure level values

    DO jk = 1, nplvp1
      zp = vct(jk)
      DO jl = 1, klen
        ph(jl,jk) = zp
      END DO
    END DO

    ! Compute hybrid level values

    DO jk = nplvp2, nlmsgl
      zp = vct(jk)
      zb = vct(jk+nvclev)
      DO jl = 1, klen
        ph(jl,jk) = zp + zb*ps(jl)
      END DO
    END DO

    ! Compute sigma-level values

    DO jk = nlmslp, nlevp1
      zb = vct(jk+nvclev)
      DO jl = 1, klen
        ph(jl,jk) = zb*ps(jl)
      END DO
    END DO

  END SUBROUTINE half_level_pressure

  !>
  !! Calculate full-level pressures for all vertical layers
  !! Method: Simmons and Burridge (Mon.Wea.Rev.,1981,p761,Eqn.(3.18))
  !!
  !! @par Parameters
  !!    *pres_i*    half-level pressure values.
  !!    *pres_m*    computed full-level pressure values.
  !!    *kdimp*     first dimension of 2-d array *pres_i*
  !!    *klen*      number of points for which calculation is
  !!                performed.
  !!
  !!  Required constants are obtained from module *mo_hyb*.
  !!  The latter must have been initialiazed
  !!  by a call of subroutine *inihyb*.
  !!
  !! @par Revision History
  !!    H. Wan, MPI, 2006-08-17
  !!
  SUBROUTINE full_level_pressure( pres_i, kdimp, klen, nlev_in, pres_m)

    INTEGER ,INTENT(in) :: kdimp, klen, nlev_in    !< dimension parameters
    REAL(wp),INTENT(in) :: pres_i(kdimp,nlev_in+1) !< half-level pressure

    REAL(wp),INTENT(inout) :: pres_m(kdimp,nlev_in) !< full(/mid)-level pressure

    REAL(wp):: ztmp, zpres_i_top_min
    INTEGER :: jk, jl, ikp, ik_top, nlev

    !-----

    nlev = nlev_in
    zpres_i_top_min = vct(1)
    IF ( zpres_i_top_min > 0._wp ) THEN
      ik_top = 1
    ELSE
      ik_top = 2
      pres_m(1:klen,1) = pres_i(1:klen,2)*0.5_wp
    END IF

    DO jk = ik_top, nlev
       ikp = jk+1
       DO jl = 1, klen
         ztmp = ( pres_i(jl,ikp)*LOG(pres_i(jl,ikp))   &
         &       -pres_i(jl,jk )*LOG(pres_i(jl,jk )) ) &
         &     /( pres_i(jl,ikp)-pres_i(jl,jk) )
         pres_m(jl,jk) = EXP(ztmp-1._wp)
       END DO
    END DO

  END SUBROUTINE full_level_pressure

  !-------------------------------------------------------------------------
  !>
  !! Calculates auxiliary variables connected with
  !! the vertical finite-difference scheme.
  !!
  !! @par Arguments
  !!
  !!   *ph*          *specified half-level pressures.
  !!   *kdim*        *first dimension of 2-d arrays *pdelp,*
  !!                  *plnpr,* *palpha,* and *ph.*
  !!   *klen*        *number of points for which calculation is performed.
  !!   *pdelp*       *computed pressure difference across layers.
  !!   *prdelp*      *reciprocal of *pdelp.*
  !!   *plnpr*       *computed logarithm of ratio of pressures.
  !!   *palpha*      *computed alphas for integration of the
  !!                  hydrostatic equation and related terms.
  !!
  !!  Required constants are obtained from modules
  !!  *mo_constants* and *mo_hyb*. The latter should have been
  !!  initialized by a call of subroutine *inihyb*.
  !!
  !!  Results are computed for *klen* consecutive points for
  !!  all required levels.
  !!
  !!  Calculations are performed separately for pressure,
  !!  hybrid and sigma levels.
  !!
  !!  External documentation of the model equations and the
  !!  organization of the vertical calculation can be found
  !!  in MPI-M technical report No. 349
  !!
  !! @par Revision History
  !!    A. J. Simmons, ECMWF, November 1981, original source
  !!    L. Kornblueh, MPI, May 1998, f90 rewrite
  !!    U. Schulzweida, MPI, May 1998, f90 rewrite
  !!
  SUBROUTINE auxhyb( ph,kdim,klen,nlev_in,                   &
                     pdelp,prdelp,plnph,plnpr,palpha )

    INTEGER ,INTENT(in)  :: kdim, klen, nlev_in
    REAL(wp),INTENT(in)  :: ph(kdim, nlev_in+1)

    REAL(wp), INTENT(inout) :: pdelp (kdim, nlev_in  ), prdelp(kdim, nlev_in)
    REAL(wp), INTENT(inout) :: plnph (kdim, nlev_in+1), plnpr (kdim, nlev_in)
    REAL(wp), INTENT(inout) :: palpha(kdim, nlev_in  )

    REAL(wp) :: za, zd, zl, zr
    INTEGER  :: jk, jl, nlev, nlevp1

    nlev = nlev_in
    nlevp1 = nlev+1

    !-----
    ! Set pressure-level values or other top-level values

    DO jk = 1, nplev
      zd = delpr(jk)
      zr = rdelpr(jk)
      zl = rlnpr(jk)
      za = ralpha(jk)
      DO jl = 1, klen
        pdelp(jl,jk) = zd
        prdelp(jl,jk) = zr
        plnpr(jl,jk) = zl
        palpha(jl,jk) = za
      END DO
    END DO

    ! Calculate hybrid-level values

    DO jk = nplvp1, nlmsgl

      IF (jk == 1 .AND. vct(1) == 0.0_wp ) THEN

        DO jl = 1, klen
          pdelp(jl,jk)  = ph(jl,jk+1) - ph(jl,jk)
          prdelp(jl,jk) = 1._wp/pdelp(jl,jk)
          palpha(jl,jk) = rd*LOG(2._wp)
          plnpr(jl,jk)  = 2._wp*palpha(jl,jk)
        END DO

      ELSE
        DO jl = 1, klen
          pdelp(jl,jk)  = ph(jl,jk+1) - ph(jl,jk)
          prdelp(jl,jk) = 1._wp/pdelp(jl,jk)
          plnpr(jl,jk)  = rd*LOG(ph(jl,jk+1)/ph(jl,jk))
          palpha(jl,jk) = rd - ph(jl,jk)*plnpr(jl,jk)*prdelp(jl,jk)
        END DO

      ENDIF
    END DO

    ! Set sigma-level values

    DO jk = nlmsla, nlev
      zl = rlnpr(jk)
      za = ralpha(jk)
      DO jl = 1, klen
        pdelp(jl,jk)  = ph(jl,jk+1) - ph(jl,jk)
        prdelp(jl,jk) = 1._wp/pdelp(jl,jk)
        plnpr(jl,jk)  = zl
        palpha(jl,jk) = za
      END DO
    END DO

    DO jk = 2,nlevp1
      DO jl = 1, klen
        plnph(jl,jk) = LOG(ph(jl,jk))
      END DO
    END DO

  END SUBROUTINE auxhyb

  !-------------------------------------------------------------------------
  !>
  !! Calculates full- and half-level geopotential
  !!
  !! Method: Integrate the hydrostatic equation in the vertical
  !! to obtain full- or half-level values of geopotential.
  !!
  !! *geopot* is called during the calculation of adiabatic
  !! tendencies, prior to the calculation of the physical paramet-
  !! erizations, and in the post-processing.
  !!
  !! Parameters are
  !!  *pgeop_sfc*   *surface geopotential.
  !!  *pgeop_m*     *computed geopotential on full-levels
  !!  *pgeop_i*     *computed geopotential on half-levels
  !!  *ptv*         *virtual temperature.
  !!  *plnpr*       *logarithm of ratio of pressures, computed by *auxhyb*.
  !!  *palpha*      *for full-level values use *alpha* as computed by *auxhyb*.
  !!  *kdim*        *first dimension of 2-d arrays *phi,*
  !!                 *ptv,* *plnpr,* and *palpha.*
  !!  *kstart,kend* *start and end points for which calculation is performed.
  !!
  !! Required constants are obtained from module *mo_hyb*.
  !! The latter should have been initialized by a call of subroutine *inihyb*.
  !!
  !! Results are computed for *klen* consecutive points for *nlev* levels.
  !!
  !! The choice of full- or half-level values is determined
  !! by the specification of the input array *alpha.*
  !!
  !! External documentation of the model equations and the
  !! organization of the vertical calculation.
  !!
  !! @par Revision History
  !!    A. J. Simmons, ECMWF, January 1982, original source
  !!    L. Kornblueh, MPI, May 1998, f90 rewrite
  !!    U. Schulzweida, MPI, May 1998, f90 rewrite
  !!    H. Wan, MPI-M, 2006-08-10, included in m_vertical
  !!    H. Wan, MPI-M, 2006-08-15, half-level geopotential added to the output
  !!    H. Wan, MPI-M, 2007-07-19, calculation of the full-level geopotential
  !!                               simplified.
  !!    G. Zaengl, DWD, 2009-07-06, replace klen with kstart/kend for correct
  !!                                execution on nested domains
  !!    A. Seifert, DWD, 2010-06-21, add missing lower boundary by restructuring
  !!                                 the loop for half and full levels

  SUBROUTINE geopot( ptv,plnpr,palpha,pgeop_sfc,kdim,kstart,kend,nlev_in, &
                     pgeop_m, pgeop_i )

    INTEGER ,INTENT(in) :: kdim, kstart, kend, nlev_in

    REAL(wp) ,INTENT(in)    :: ptv   (kdim,nlev_in),     plnpr(kdim,nlev_in)
    REAL(wp) ,INTENT(in)    :: palpha(kdim,nlev_in), pgeop_sfc(kdim)
    REAL(wp) ,INTENT(inout) :: pgeop_m(kdim,nlev_in),   pgeop_i(kdim,nlev_in+1)

    INTEGER :: jk, jl, jkp, nlev, nlevp1

    nlev = nlev_in
    nlevp1 = nlev+1

    ! Integrate hydrostatic equation

    DO jl = kstart, kend
      pgeop_i(jl,nlevp1) = pgeop_sfc(jl)
      pgeop_m(jl,nlev)   = palpha(jl,nlev)*ptv(jl,nlev) + pgeop_sfc(jl)
    END DO

    ! half levels
    DO jk = nlev, 1, -1
      jkp = jk + 1
      DO jl = kstart, kend
        pgeop_i(jl,jk) = plnpr(jl,jk)*ptv(jl,jk) + pgeop_i(jl,jkp)
      END DO
    END DO

    ! full levels
    DO jk = nlevm1, 1, -1
      jkp = jk + 1
      DO jl = kstart, kend
        pgeop_m(jl,jk)  = pgeop_i(jl,jkp) + palpha(jl,jk)*ptv(jl,jk)
      END DO
    END DO

  END SUBROUTINE geopot

END MODULE mo_ifs_coord

