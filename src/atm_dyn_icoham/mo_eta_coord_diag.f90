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
MODULE mo_eta_coord_diag
!RJ: For some strange reasons the Intel compiler produces wrong code when
!optimizing convert_theta2t_lin, therefore we switch off optimization here:
!DEC$ NOOPTIMIZE

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: rd
  USE mo_run_config,         ONLY: nlevp1, nvclev, nlev
  USE mo_vertical_coord_table, ONLY: nlmsgl, nlmslp, nplvp1, nplvp2, vct, &
    &                                delpr,  nlmsla, nplev,  ralpha,      &
    &                                rdelpr, rlnpr,  nlevm1

  USE mo_fast_math_lib,      ONLY: vec_log, vec_exp
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: half_level_pressure, full_level_pressure
  PUBLIC :: auxhyb, geopot

CONTAINS
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
  SUBROUTINE half_level_pressure( ps,kdimp,klen,ph)

    INTEGER ,INTENT(in)  :: kdimp
    REAL(wp),INTENT(in)  :: ps(kdimp)   !< surface pressure
    INTEGER ,INTENT(in)  :: klen

    REAL(wp),INTENT(inout) :: ph(kdimp,nlevp1) !< half-level pressure

    REAL(wp):: zb, zp
    INTEGER :: jk, jl

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
  SUBROUTINE full_level_pressure( pres_i, kdimp, klen, pres_m)

    INTEGER ,INTENT(in) :: kdimp, klen          !< dimension parameters
    REAL(wp),INTENT(in) :: pres_i(kdimp,nlevp1) !< half-level pressure

    REAL(wp),INTENT(inout) :: pres_m(kdimp,nlev) !< full(/mid)-level pressure

    REAL(wp):: zpres_i_top_min
    REAL(wp):: ztmp(kdimp,nlev)
    REAL(wp):: log_press(kdimp,nlevp1)
    INTEGER :: jk, jl, ikp, ik_top

    !-----

    zpres_i_top_min = vct(1)
    IF ( zpres_i_top_min > 0._wp ) THEN
      ik_top = 1
    ELSE
      ik_top = 2
      pres_m(1:klen,1) = pres_i(1:klen,2)*0.5_wp
    END IF

    DO jk = ik_top, nlev+1
      CALL vec_log(pres_i(1:klen,jk), log_press(1:klen,jk), klen)
    ENDDO
    
    DO jk = ik_top, nlev
       ikp = jk+1
       DO jl = 1, klen
         ztmp(jl, jk) = (( pres_i(jl,ikp)*log_press(jl,ikp)   &
         &       -pres_i(jl,jk )*log_press(jl,jk ) ) &
         &     /( pres_i(jl,ikp)-pres_i(jl,jk) )) - 1._wp
       END DO
    END DO
    
    DO jk = ik_top, nlev
      CALL vec_exp(ztmp(1:klen, jk), pres_m(1:klen,jk), klen)
    ENDDO

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
  SUBROUTINE auxhyb( ph,kdim,klen,                   &
                     pdelp,prdelp,plnph,plnpr,palpha )

    INTEGER ,INTENT(in)  :: kdim, klen
    REAL(wp),INTENT(in)  :: ph(kdim, nlevp1)

    REAL(wp), INTENT(inout) :: pdelp (kdim, nlev  ), prdelp(kdim, nlev)
    REAL(wp), INTENT(inout) :: plnph (kdim, nlevp1), plnpr (kdim, nlev)
    REAL(wp), INTENT(inout) :: palpha(kdim, nlev  )

    REAL(wp) :: za, zd, zl, zr
    INTEGER  :: jk, jl

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

  SUBROUTINE geopot( ptv,plnpr,palpha,pgeop_sfc,kdim,kstart,kend, &
                     pgeop_m, pgeop_i )

    INTEGER ,INTENT(in) :: kdim, kstart, kend

    REAL(wp) ,INTENT(in)    :: ptv   (kdim,nlev),     plnpr(kdim,nlev)
    REAL(wp) ,INTENT(in)    :: palpha(kdim,nlev), pgeop_sfc(kdim)
    REAL(wp) ,INTENT(inout) :: pgeop_m(kdim,nlev),   pgeop_i(kdim,nlevp1)

    INTEGER :: jk, jl, jkp

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

END MODULE mo_eta_coord_diag

