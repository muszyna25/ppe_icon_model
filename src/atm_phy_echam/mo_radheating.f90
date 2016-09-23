!>
!! @brief Module to compute the radiative heating in W/m2.
!!
!! @remarks
!!   This code is derived from the overloaded "radheat" subroutine
!!   of mo_radiation.
!!
!! @author Marco Giorgetta, MPI-M, Hamburg (2009-09-192016-09-19):
!!
!! @par Origin
!!   This code is derived from the overloaded "radheat" subroutine
!!   of atm_phy_schemes/mo_radiation.f90.
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

  USE mo_kind,                 ONLY: wp
  USE mo_physical_constants,   ONLY: stbo

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: radheating


  ! --- radiative transfer parameters
  !
  REAL(wp), PARAMETER :: diff   = 1.66_wp   !< LW Diffusivity Factor

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

  SUBROUTINE radheating (jcs, jce, kbdim, &
    &                    klev  , klevp1,  &
    &                    pi0           ,  &
    &                    pemiss        ,  &
    &                    ptsfc         ,  &
    &                    ptsfctrad     ,  &
    &                    ptrmsw        ,  &
    &                    pflxlw        ,  &
    &                    lwflx_up_sfc_rs, & ! optional: longwave upward flux at surface
    &                    ptrmswclr     ,  & ! optional: shortwave net transmissivity at last rad. step clear sky []
    &                    pflxlwclr     ,  & ! optional: longwave net flux at last rad. step clear sky [W/m2]
    &                    pq_rsw        ,  &
    &                    pq_rlw        ,  &
    &                    pflxsfcsw     ,  &
    &                    pflxsfclw     ,  &
    &                    lwflx_up_sfc  ,  &
    &                    pflxtoasw     ,  &
    &                    pflxtoalw      )

    INTEGER,  INTENT(in)  ::    &
      &     jcs, jce, kbdim,    &
      &     klev,   klevp1

    REAL(wp), INTENT(in)  ::           &
      &     pi0        (kbdim),        & ! local solar incoming flux at TOA         [W/m2]
      &     pemiss     (kbdim),        & ! lw sfc emissivity
      &     ptsfc      (kbdim),        & ! surface temperature at t                 [K]
      &     ptsfctrad  (kbdim),        & ! surface temperature at trad              [K]
      &     ptrmsw     (kbdim,klevp1), & ! shortwave transmissivity at trad         []
      &     pflxlw     (kbdim,klevp1)    ! longwave net flux at trad                [W/m2]

    REAL(wp), INTENT(in)  ::           &
      &     lwflx_up_sfc_rs(kbdim)       ! longwave upward flux at surface calculated at radiation time steps

    REAL(wp), INTENT(in)  ::           &
      &     ptrmswclr  (kbdim,klevp1), & ! shortwave net transmissivity at last rad. step clear sky []
      &     pflxlwclr  (kbdim,klevp1)    ! longwave net flux at last rad. step clear sky [W/m2]
   
    REAL(wp), INTENT(inout) ::   &
      &     pq_rsw (kbdim,klev), & ! heating by shortwave radiation  [W/m2]
      &     pq_rlw (kbdim,klev)    ! heating by longwave  radiation  [W/m2]         [K/s]

    REAL(wp), INTENT(inout) :: &
      &     pflxsfcsw (kbdim), &       ! shortwave surface net flux [W/m2]
      &     pflxsfclw (kbdim), &       ! longwave  surface net flux [W/m2]
      &     pflxtoasw (kbdim), &       ! shortwave toa net flux [W/m2]
      &     pflxtoalw (kbdim), &       ! longwave  toa net flux [W/m2]
      &     lwflx_up_sfc(kbdim)        ! longwave upward flux at surface [W/m2]

    ! Local arrays
    REAL(wp) ::                    &
      &     zflxsw (kbdim,klevp1), &
      &     zflxlw (kbdim,klevp1), &
      &     zflxswclr(kbdim,klevp1),&
      &     zflxlwclr(kbdim,klevp1),&
      &     dlwem_o_dtg(kbdim)

    INTEGER :: jk

    ! lev == 1        => TOA
    ! lev in [2,klev] => Atmosphere
    ! lev == klevp1   => Surface
    DO jk = 1, klevp1
      zflxsw   (jcs:jce,jk)  = ptrmsw   (jcs:jce,jk)*pi0(jcs:jce)
      zflxswclr(jcs:jce,jk)  = ptrmswclr(jcs:jce,jk)*pi0(jcs:jce)
    END DO

    ! Longwave fluxes: For now keep fluxes fixed at TOA and in atmosphere,
    ! but adjust flux from surface to the current surface temperature.
    !
    ! - TOA
    zflxlw(jcs:jce,1)      = pflxlw(jcs:jce,1)
    !
    ! - Atmosphere
    zflxlw(jcs:jce,2:klev) = pflxlw(jcs:jce,2:klev)
    !
    ! - Surface
    !   Adjust net sfc longwave radiation for changed surface temperature (ptsfc) with respect to the
    !   surface temperature used for the longwave flux computation (ptsfctrad).
    !   --> modifies heating in lowermost layer only (is this smart?)
    !   This assumes that downward sfc longwave radiation is constant between radiation time steps and
    !   only upward and net sfc longwave radiation are updated between radiation time steps
    dlwem_o_dtg(jcs:jce)   = pemiss(jcs:jce)*4._wp*stbo*ptsfc(jcs:jce)**3    ! Derivative of upward sfc rad wrt to sfc temperature
    lwflx_up_sfc(jcs:jce)  = lwflx_up_sfc_rs(jcs:jce)                 & ! Upward longwave sfc rad at radiation time step
      &                    + dlwem_o_dtg(jcs:jce) * (ptsfc(jcs:jce) - ptsfctrad(jcs:jce))   ! Correction for new sfc temp between radiation time steps
    zflxlw(jcs:jce,klevp1) = pflxlw(jcs:jce,klevp1)                      & ! Net longwave sfc rad at radiation time step
      &                    - dlwem_o_dtg(jcs:jce) * (ptsfc(jcs:jce) - ptsfctrad(jcs:jce))     ! Correction for new sfc temp between radiation time steps
    !
    !     4.2  Fluxes and heating rates except for lowest layer
    !
    pq_rsw(jcs:jce,1:klev) = (zflxsw(jcs:jce,1:klev)-zflxsw(jcs:jce,2:klev+1))
    pq_rlw(jcs:jce,1:klev) = (zflxlw(jcs:jce,1:klev)-zflxlw(jcs:jce,2:klev+1))
    !
    !     4.3 net fluxes at surface
    !
    pflxsfcsw(jcs:jce) = zflxsw(jcs:jce,klevp1)
    pflxsfclw(jcs:jce) = zflxlw(jcs:jce,klevp1)

    !
    !     4.4 net sw flux at toa
    !
    pflxtoasw(jcs:jce) = zflxsw(jcs:jce,1)
    pflxtoalw(jcs:jce) = zflxlw(jcs:jce,1)

  END SUBROUTINE radheating

END MODULE mo_radheating
