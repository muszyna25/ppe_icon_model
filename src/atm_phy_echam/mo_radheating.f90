!>
!! @brief Module to provide interface to radiation routines.
!!
!! @remarks
!!   This module contains routines that provide the interface between ECHAM
!!   and the radiation code.  Mostly it organizes and calculates the
!!   information necessary to call the radiative transfer solvers.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19):
!!
!!         Hauke Schmidt, MPI-M, Hamburg (2009-12-18): Few modifications to
!!              allow specific solar irradiance for AMIP-type and preindustrial
!!              simulations.
!!         Luis Kornblueh, MPI-M, Hamburg (2010-04-06): Never ever use write
!!              directly
!!         Martin Schultz, FZJ, Juelich (2010-04-13):
!!              Extracted public parameters into new module mo_radiation_parameters
!!              to avoid circular dependencies in submodels
!!                                      (2010-06-03):
!!              Added submodel calls, decl_sun_cur
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Major segments of this code combines and rewrites (for the ICON standard)
!!   code previously contained in the ECHAM5 routines rad_int.f90,
!!   radiation.f90 and prerad.f90.  Modifications were also made to provide
!!   a cleaner interface to the aerosol and cloud properties. Contributors to
!!   the code from which the present routines were derived include:  M. Jarraud,
!!   ECMWF (1983-06); M.A. Giorgetta, MPI-M (2002-05); U. Schulzweida,  MPI-M
!!   (2002-05); P. Stier MPI-M \& Caltech (2004-04, 2006-07), M. Thomas MPI-M
!!   (2007-06); U. Schlese, MPI-M (2007-06); M. Esch, MPI-M (2007-06); S.J.
!!   Lorenz, MPI-M (2007-11); T. Raddatz, MPI-M (2006-05); I. Kirchner.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!
MODULE mo_radheating

  USE mo_kind                       , ONLY: wp
  USE mo_physical_constants         , ONLY: stbo
  USE mo_psrad_radiation_parameters , ONLY: diff, psctm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: radheating


CONTAINS

  !-----------------------------------------------------------------------------
  !>
  !! Compute shortwave and longwave heating rates
  !!
  !! The radheat subroutine computes the radiative heating rates resulting from
  !! the divergence of the vertical profiles of longwave and shortwave net fluxes.
  !!
  !! - Shortwave net flux profiles are computed from:
  !!   - the vertical profiles of net transmissivity
  !!   - the solar incoming flux at TOA
  !! - Longwave net flux profiles are given as input
  !! - Specific heat depends on the moisture in the air
  !!
  !! @author Marco Giorgetta, Max Planck Institute for Meteorology
  !!
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!

  SUBROUTINE radheating ( &
       !
       ! input
       ! -----
       !
       & jcs            ,&
       & jce            ,&
       & kbdim          ,&
       & klev           ,&
       & klevp1         ,&
       !
       & cosmu0_rad     ,&
       & cosmu0         ,&
       !
       & prsdt          ,&
       & pemiss         ,&
       & ptsfc          ,&
       & ptsfctrad      ,&
       & lwflx_up_sfc_rs,&
       !
       & rsd            ,&
       & rsu            ,&
       & rld            ,&
       & rlu            ,&
       !
       & rsdcs          ,&
       & rsucs          ,&
       & rldcs          ,&
       & rlucs          ,&
       !
       & ptrmsw         ,&
       & rln            ,&
       & ptrmswclr      ,&
       & rlncs          ,&
       !
       ! output
       ! ------
       !
       & pq_rsw         ,&
       & pq_rlw         ,&
       !
       & rsns           ,&
       & rlns           ,&
       & rsnt           ,&
       & rlnt           ,&
       & lwflx_up_sfc   )

    INTEGER,  INTENT(in)  ::    &
      &     jcs, jce, kbdim,    &
      &     klev,   klevp1

    REAL(wp), INTENT(in)  ::          &
      &     cosmu0_rad (kbdim)       ,&! solar zenith angle at radiation time
      &     cosmu0     (kbdim)       ,&! solar zenith angle at current   time
      !
      &     prsdt      (kbdim)       ,&! local solar incoming flux at TOA         [W/m2]
      &     pemiss     (kbdim)       ,&! lw sfc emissivity
      &     ptsfc      (kbdim)       ,&! surface temperature at t                 [K]
      &     ptsfctrad  (kbdim)       ,&! surface temperature at trad              [K]
      &     lwflx_up_sfc_rs(kbdim),   &! longwave upward flux at surface calculated at radiation time steps
      !
      &     rsd   (kbdim,klevp1)     ,&! all-sky   shortwave downward flux at last radiation step [W/m2]
      &     rsu   (kbdim,klevp1)     ,&! all-sky   shortwave upward   flux at last radiation step [W/m2]
      &     rld   (kbdim,klevp1)     ,&! all-sky   longwave  downward flux at last radiation step [W/m2]
      &     rlu   (kbdim,klevp1)     ,&! all-sky   longwave  upward   flux at last radiation step [W/m2]
      !
      &     rsdcs (kbdim,klevp1)     ,&! clear-sky shortwave downward flux at last radiation step [W/m2]
      &     rsucs (kbdim,klevp1)     ,&! clear-sky shortwave upward   flux at last radiation step [W/m2]
      &     rldcs (kbdim,klevp1)     ,&! clear-sky longwave  downward flux at last radiation step [W/m2]
      &     rlucs (kbdim,klevp1)       ! clear-sky longwave  upward   flux at last radiation step [W/m2]
      !
    REAL(wp), INTENT(inout)  ::       &
      &     ptrmsw     (kbdim,klevp1),&! shortwave transmissivity at trad         []
      &     rln        (kbdim,klevp1),&! longwave net flux at trad                [W/m2]
      &     ptrmswclr  (kbdim,klevp1),&! shortwave net transmissivity at last rad. step clear sky []
      &     rlncs      (kbdim,klevp1)  ! longwave net flux at last rad. step clear sky [W/m2]
   
    REAL(wp), INTENT(inout) ::        &
      &     pq_rsw (kbdim,klev)      ,&! shortwave temperature heating  [W/m2]
      &     pq_rlw (kbdim,klev)        ! longwave  temperature heating  [W/m2]

    REAL(wp), INTENT(inout) ::        &
      &     rsns (kbdim)             ,&! shortwave surface net flux [W/m2]
      &     rlns (kbdim)             ,&! longwave  surface net flux [W/m2]
      &     rsnt (kbdim)             ,&! shortwave toa net flux [W/m2]
      &     rlnt (kbdim)             ,&! longwave  toa net flux [W/m2]
      &     lwflx_up_sfc(kbdim)        ! longwave upward flux at surface [W/m2]

    ! Local arrays
    REAL(wp) ::                     &
      &     zflxsw (kbdim,klevp1)  ,&
      &     zflxlw (kbdim,klevp1)  ,&
      &     zflxswclr(kbdim,klevp1),&
      &     zflxlwclr(kbdim,klevp1),&
      &     dlwem_o_dtg(kbdim)     ,&
      &     flux_factor(kbdim)

    ! local scalars
    INTEGER :: jk

    ! lev == 1        => TOA
    ! lev in [2,klev] => Atmosphere
    ! lev == klevp1   => Surface

!!$    flux_factor(jcs:jce) = 1._wp / (psctm*cosmu0_rad(jcs:jce))
    
    ! Shortwave fluxes
    DO jk = 1, klevp1
!!$      ptrmsw   (jcs:jce,jk)  = (rsd  (jcs:jce,jk) - rsu  (jcs:jce,jk)) * flux_factor(jcs:jce)
!!$      ptrmswclr(jcs:jce,jk)  = (rsdcs(jcs:jce,jk) - rsucs(jcs:jce,jk)) * flux_factor(jcs:jce)
      !
      zflxsw   (jcs:jce,jk)  = ptrmsw   (jcs:jce,jk) * prsdt(jcs:jce)
      zflxswclr(jcs:jce,jk)  = ptrmswclr(jcs:jce,jk) * prsdt(jcs:jce)
    END DO

    ! Longwave fluxes

      ! Longwave fluxes: For now keep fluxes fixed at TOA and in atmosphere,
      ! but adjust flux from surface to the current surface temperature.
      ! - TOA
      zflxlw(jcs:jce,1)      = rln(jcs:jce,1)
      ! - Atmosphere
      zflxlw(jcs:jce,2:klev) = rln(jcs:jce,2:klev)

      ! - Surface
      !   Adjust net sfc longwave radiation for changed surface temperature (ptsfc) with respect to the
      !   surface temperature used for the longwave flux computation (ptsfctrad).
      !   --> modifies heating in lowermost layer only (is this smart?)
      !   This assumes that downward sfc longwave radiation is constant between radiation time steps and
      !   upward and net sfc longwave radiation are updated between radiation time steps
      dlwem_o_dtg(jcs:jce) = pemiss(jcs:jce)*4._wp*stbo*ptsfc(jcs:jce)**3    ! Derivative of upward sfc rad wrt to sfc temperature
      lwflx_up_sfc(jcs:jce) = lwflx_up_sfc_rs(jcs:jce)                 & ! Upward longwave sfc rad at radiation time step
        &   + dlwem_o_dtg(jcs:jce) * (ptsfc(jcs:jce) - ptsfctrad(jcs:jce))   ! Correction for new sfc temp between radiation time steps
      zflxlw(jcs:jce,klevp1) = rln(jcs:jce,klevp1)                      & ! Net longwave sfc rad at radiation time step
        & - dlwem_o_dtg(jcs:jce) * (ptsfc(jcs:jce) - ptsfctrad(jcs:jce))     ! Correction for new sfc temp between radiation time steps
!!$      zflxlw(jcs:jce,klevp1) = rln(jcs:jce,klevp1)                      &
!!$        &                   + pemiss(jcs:jce)*stbo * ptsfctrad(jcs:jce)**4 &
!!$        &                   - pemiss(jcs:jce)*stbo * ptsfc    (jcs:jce)**4
        ! Longwave fluxes clear sky: For now keep fluxes fixed at TOA and in atmosphere,
        ! but adjust flux from surface to the current surface temperature.
        ! - TOA
        zflxlwclr(jcs:jce,1)      = rlncs(jcs:jce,1)
        ! - Atmosphere
        zflxlwclr(jcs:jce,2:klev) = rlncs(jcs:jce,2:klev)

        ! - Surface
        !   Adjust net sfc longwave radiation for changed surface temperature (ptsfc) with respect to the
        !   surface temperature used for the longwave flux computation (ptsfctrad).
        !   --> modifies heating in lowermost layer only (is this smart?)
        !   This assumes that downward sfc longwave radiation is constant between radiation time steps and
        !   upward and net sfc longwave radiation are updated between radiation time steps
        dlwem_o_dtg(jcs:jce) = pemiss(jcs:jce)*4._wp*stbo*ptsfc(jcs:jce)**3    ! Derivative of upward sfc rad wrt to sfc temperature
        zflxlwclr(jcs:jce,klevp1) = rlncs(jcs:jce,klevp1)                & ! Net longwave sfc rad at radiation time step
        & - dlwem_o_dtg(jcs:jce) * (ptsfc(jcs:jce) - ptsfctrad(jcs:jce))       ! Correction for new sfc temp between radiation time steps

    !
    !
    !     4.2  Fluxes and heating rates except for lowest layer
    !
    pq_rsw(jcs:jce,1:klev) = zflxsw(jcs:jce,1:klev)-zflxsw(jcs:jce,2:klev+1)
    pq_rlw(jcs:jce,1:klev) = zflxlw(jcs:jce,1:klev)-zflxlw(jcs:jce,2:klev+1)

    !
    !     4.3 net fluxes at surface
    !
    rsns(jcs:jce) = zflxsw(jcs:jce,klevp1)
    rlns(jcs:jce) = zflxlw(jcs:jce,klevp1)

    !
    !     4.4 net sw flux at toa
    !
    rsnt(jcs:jce) = zflxsw(jcs:jce,1)
    rlnt(jcs:jce) = zflxlw(jcs:jce,1)

  END SUBROUTINE radheating

END MODULE mo_radheating
