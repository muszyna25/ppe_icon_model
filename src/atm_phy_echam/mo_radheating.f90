!>
!! @brief Module to provide SW and LW fluxes and heating rates
!!
!! @remarks
!!   This module contains the "radheating" routine that diagnoses
!!   SW and LW fluxes at TOA and at the surface and the heating
!!   in the atmosphere.
!!
!! @author Marco Giorgetta, MPI-M, Hamburg (2016-11-02)
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

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: radheating


CONTAINS

  !-----------------------------------------------------------------------------
  !>
  !! Compute shortwave and longwave heating rates
  !!
  !! The radheating subroutine computes the radiative heating rates in the
  !! shortwave (SW) and longwave (LW) part of the spectrum.
  !!
  !! SW and LW fluxes are computed only every radiation time step. The resulting
  !! upward and downward fluxes at all half-levels are used here to approximate
  !! fluxes at the current time, as follows:
  !!
  !! SW fluxes are rescaled by the SW indicent radiation at the top of the
  !! atmosphere. In order to make this work the radiative transfer uses a zenith
  !! angle limited to 84deg so that also at nighttime non-zero fluxes exist,
  !! as needed for the scaling.
  !!
  !! LW fluxes are kept constant except for the upward flux from the surface,
  !! which is corrected for the new radiative surface temperature.
  !!
  !! The corrected fluxes are used to diagnose various SW and LW flux components
  !! at the top of the atmosphere and at the surface, and to compute the SW and LW
  !! heating from the vertical flux convergences.
  !!
  !! The vertical ordering is assumed to be from the top of the atmosphere
  !! downward to the surface.
  !!
  !! @author Marco Giorgetta, Max Planck Institute for Meteorology
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!

  SUBROUTINE radheating ( &
       !
       ! input
       ! -----
       !
       & jcs        ,&
       & jce        ,&
       & kbdim      ,&
       & klev       ,&
       & klevp1     ,&
       !
       & rsdt0      ,&
       & cosmu0     ,&
       !
       & emiss      ,&
       & tsr        ,&
       & tsr_rt     ,&
       !
       & rsd        ,&
       & rsu        ,&
       & rld        ,&
       & rlu        ,&
       !
       & rsdcs      ,&
       & rsucs      ,&
       & rldcs      ,&
       & rlucs      ,&
       !
       ! output
       ! ------
       !
       & rsdt       ,&
       & rsut       ,&
       & rsnt       ,&
       & rsds       ,&
       & rsus       ,&
       & rsns       ,&
       !
       & rsutcs     ,&
       & rsntcs     ,&
       & rsdscs     ,&
       & rsuscs     ,&
       & rsnscs     ,&
       !
       & rlut       ,&
       & rlnt       ,&
       & rlds       ,&
       & rlus       ,&
       & rlns       ,&
       !
       & rlutcs     ,&
       & rlntcs     ,&
       & rldscs     ,&
       & rlnscs     ,&
       !
       & q_rsw      ,&
       & q_rlw      )

    
    INTEGER,  INTENT(in)  :: &
      &     jcs, jce, kbdim, &
      &     klev, klevp1

    REAL(wp), INTENT(in)  ::        &
      &     rsdt0                  ,&! indicent SW flux for sun in zenith
      &     cosmu0(kbdim)          ,&! cosine of solar zenith angle at current time
      &     emiss (kbdim)          ,&! lw sfc emissivity
      &     tsr   (kbdim)          ,&! radiative surface temperature at current   time [K]
      &     tsr_rt(kbdim)          ,&! radiative surface temperature at radiation time [K]
      !
      &     rsd   (kbdim,klevp1)   ,&! all-sky   shortwave downward flux at radiation time [W/m2]
      &     rsu   (kbdim,klevp1)   ,&! all-sky   shortwave upward   flux at radiation time [W/m2]
      !
      &     rsdcs (kbdim,klevp1)   ,&! clear-sky shortwave downward flux at radiation time [W/m2]
      &     rsucs (kbdim,klevp1)   ,&! clear-sky shortwave upward   flux at radiation time [W/m2]
      !
      &     rld   (kbdim,klevp1)   ,&! all-sky   longwave  downward flux at radiation time [W/m2]
      &     rlu   (kbdim,klevp1)   ,&! all-sky   longwave  upward   flux at radiation time [W/m2]
      !
      &     rldcs (kbdim,klevp1)   ,&! clear-sky longwave  downward flux at radiation time [W/m2]
      &     rlucs (kbdim,klevp1)     ! clear-sky longwave  upward   flux at radiation time [W/m2]
      !
    REAL(wp), INTENT(out) ::        &
      &     rsdt  (kbdim)          ,&! all-sky   shortwave downward flux at current   time [W/m2]
      &     rsut  (kbdim)          ,&! all-sky   shortwave upward   flux at current   time [W/m2]
      &     rsnt  (kbdim)          ,&! all-sky   shortwave net      flux at current   time [W/m2]
      &     rsds  (kbdim)          ,&! all-sky   shortwave downward flux at current   time [W/m2]
      &     rsus  (kbdim)          ,&! all-sky   shortwave upward   flux at current   time [W/m2]
      &     rsns  (kbdim)          ,&! all-sky   shortwave net      flux at current   time [W/m2]
      !
      &     rsutcs(kbdim)          ,&! clear-sky shortwave upward   flux at current   time [W/m2]
      &     rsntcs(kbdim)          ,&! clear-sky shortwave net      flux at current   time [W/m2]
      &     rsdscs(kbdim)          ,&! clear-sky shortwave downward flux at current   time [W/m2]
      &     rsuscs(kbdim)          ,&! clear-sky shortwave upward   flux at current   time [W/m2]
      &     rsnscs(kbdim)          ,&! clear-sky shortwave net      flux at current   time [W/m2]
      !
      &     rlut  (kbdim)          ,&! all-sky   longwave  upward   flux at current   time [W/m2]
      &     rlnt  (kbdim)          ,&! all-sky   longwave  net      flux at current   time [W/m2]
      &     rlds  (kbdim)          ,&! all-sky   longwave  downward flux at current   time [W/m2]
      &     rlus  (kbdim)          ,&! all-sky   longwave  upward   flux at current   time [W/m2]
      &     rlns  (kbdim)          ,&! all-sky   longwave  net      flux at current   time [W/m2]
      !
      &     rlutcs(kbdim)          ,&! clear-sky longwave  upward   flux at current   time [W/m2]
      &     rlntcs(kbdim)          ,&! clear-sky longwave  net      flux at current   time [W/m2]
      &     rldscs(kbdim)          ,&! clear-sky longwave  downward flux at current   time [W/m2]
      &     rlnscs(kbdim)          ,&! clear-sky longwave  net      flux at current   time [W/m2]
      !
      &     q_rsw (kbdim,klev)     ,&! radiative shortwave heating  [W/m2]
      &     q_rlw (kbdim,klev)       ! radiative longwave  heating  [W/m2]

    ! Local arrays
    REAL(wp) ::                     &
      &     xsdt  (kbdim)          ,&
      &     rsn   (kbdim,klevp1)   ,&
      &     rln   (kbdim,klevp1)   ,&
      &     drlus_dtsr(kbdim)      ,&
      &     dtsr  (kbdim)

    ! Shortwave fluxes
    ! ----------------
    !
    ! The original downward and upward fluxes form the radiative transfer (rt) calculation
    ! are scaled by the ratio of the incident solar fluxes of the current time and the
    ! rt-time. This assumes that the incident solar flux at rt-time is non-zero in all
    ! columns.
    ! - incident solar radiation at rt-time     : rsdt_rt = rsd(jk=1)
    ! - incident solar radiation at current time: rsdt    = rsdt0*MAX(0,cosmu0)
    ! - scaling ratio for fluxes at current time: xsdt    = rsdt / rsdt_rt
    !
    ! top of atmophere
    rsdt  (jcs:jce)   = rsdt0*MAX(0._wp,cosmu0(jcs:jce))
    xsdt  (jcs:jce)   = rsdt  (jcs:jce) / rsd  (jcs:jce,1)
    !
    rsut  (jcs:jce)   = rsu   (jcs:jce,1) * xsdt(jcs:jce)
    rsutcs(jcs:jce)   = rsucs (jcs:jce,1) * xsdt(jcs:jce)
    !
    rsnt  (jcs:jce)   = rsdt  (jcs:jce) - rsut  (jcs:jce)
    rsntcs(jcs:jce)   = rsdt  (jcs:jce) - rsutcs(jcs:jce)
    !
    ! all half levels
    rsn   (jcs:jce,:) = (rsd  (jcs:jce,:) - rsu  (jcs:jce,:)) * SPREAD(xsdt(jcs:jce),2,klevp1)
    !
    ! surface
    rsds  (jcs:jce)   = rsd   (jcs:jce,klevp1) * xsdt(jcs:jce)
    rsdscs(jcs:jce)   = rsdcs (jcs:jce,klevp1) * xsdt(jcs:jce)
    !
    rsus  (jcs:jce)   = rsu   (jcs:jce,klevp1) * xsdt(jcs:jce)
    rsuscs(jcs:jce)   = rsucs (jcs:jce,klevp1) * xsdt(jcs:jce)
    !
    rsns  (jcs:jce)   = rsds  (jcs:jce) - rsus  (jcs:jce)
    rsnscs(jcs:jce)   = rsdscs(jcs:jce) - rsuscs(jcs:jce)
    

    ! Longwave fluxes
    ! ---------------
    !
    ! The original downward and upward fluxes form the radiative transfer (rt) calculation
    ! are kept constant, except for the upward flux from the surface, which is corrected
    ! for the change in radiative surface temperature using a 1st order Taylor expansion.
    ! - surface upward flux at rt-time      : rlus_rt = rlu(jk=klevp1)
    ! - rad. surface temp.  at rt-time      : tsr_rt
    ! - rad. surface temp.  at current time : tsr
    !
    ! top of atmophere
    rlut  (jcs:jce)   = rlu  (jcs:jce,1)
    rlutcs(jcs:jce)   = rlucs(jcs:jce,1)
    !
    rlnt  (jcs:jce)   = -rlut  (jcs:jce)
    rlntcs(jcs:jce)   = -rlutcs(jcs:jce)
    
    ! all half levels
    rln   (jcs:jce,:) = (rld  (jcs:jce,:) - rlu  (jcs:jce,:))
    !
    ! surface
    rlds  (jcs:jce)  = rld  (jcs:jce,klevp1)
    rldscs(jcs:jce)  = rldcs(jcs:jce,klevp1)
    !
    ! - correct upward flux for changed radiative surface temperature
    drlus_dtsr(jcs:jce) = emiss(jcs:jce)*4._wp*stbo*tsr(jcs:jce)**3 ! derivative
    dtsr      (jcs:jce) = tsr(jcs:jce) - tsr_rt(jcs:jce)            ! change in tsr
    rlus(jcs:jce)       = rlu(jcs:jce,klevp1)                     & ! rlus = rlus_rt
      &                 + drlus_dtsr(jcs:jce) * dtsr(jcs:jce)       !       + correction
    !
    rlns  (jcs:jce) = rlds  (jcs:jce) - rlus  (jcs:jce)
    rlnscs(jcs:jce) = rldscs(jcs:jce) - rlus  (jcs:jce)

    
    ! Heating rates in atmosphere
    !----------------------------
    q_rsw(jcs:jce,1:klev) = rsn(jcs:jce,1:klev)-rsn(jcs:jce,2:klev+1)
    q_rlw(jcs:jce,1:klev) = rln(jcs:jce,1:klev)-rln(jcs:jce,2:klev+1)

  END SUBROUTINE radheating

END MODULE mo_radheating
