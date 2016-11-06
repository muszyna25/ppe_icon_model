!>
!! @brief Module to provide SW and LW fluxes and heating rates
!!
!! @remarks
!!   This module contains the "radheating" routine to 
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
  USE mo_psrad_radiation_parameters , ONLY: diff, psct

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
    REAL(wp), INTENT(inout) ::      &
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
      &     rsncs (kbdim,klevp1)   ,&
      &     rln   (kbdim,klevp1)   ,&
      &     rlncs (kbdim,klevp1)   ,&
      &     drlus_dtsr(kbdim)      ,&
      &     dtsr  (kbdim)

    ! local scalars
    INTEGER :: jk

    ! Shortwave fluxes
    ! ----------------
    !
    ! top of atmophere (toa)
    ! - toa incident SW radiation at the radiation time step: rsd(jk=1)
    ! - toa incident SW radiation at the current   time step: rsdt
    ! --> use rsdt/rsd(jk=1) for scaling of SW fluxes
    rsdt(jcs:jce) = MAX(0._wp,cosmu0(jcs:jce)) * psct
    xsdt(jcs:jce) = rsdt(jcs:jce) / rsd(jcs:jce,1)
    !
    rsut  (jcs:jce) = xsdt(jcs:jce) * rsu  (jcs:jce,1)
    rsutcs(jcs:jce) = xsdt(jcs:jce) * rsucs(jcs:jce,1)
    !
    rsnt  (jcs:jce) = rsdt  (jcs:jce) - rsut  (jcs:jce)
    rsntcs(jcs:jce) = rsdt  (jcs:jce) - rsutcs(jcs:jce)
    !
    ! all half levels
    DO jk = 1, klevp1
      rsn  (jcs:jce,jk) = xsdt(jcs:jce) * (rsd  (jcs:jce,jk) - rsu  (jcs:jce,jk))
      rsncs(jcs:jce,jk) = xsdt(jcs:jce) * (rsdcs(jcs:jce,jk) - rsucs(jcs:jce,jk))
    END DO
    !
    ! surface
    rsds  (jcs:jce) = xsdt(jcs:jce) * rsd  (jcs:jce,klevp1)
    rsdscs(jcs:jce) = xsdt(jcs:jce) * rsdcs(jcs:jce,klevp1)
    !
    rsus  (jcs:jce) = xsdt(jcs:jce) * rsu  (jcs:jce,klevp1)
    rsuscs(jcs:jce) = xsdt(jcs:jce) * rsucs(jcs:jce,klevp1)
    !
    rsns  (jcs:jce) = rsds  (jcs:jce) - rsus  (jcs:jce)
    rsnscs(jcs:jce) = rsdscs(jcs:jce) - rsuscs(jcs:jce)
    

    ! Longwave fluxes
    ! ---------------
    
    ! top of atmophere (toa)
    rlut  (jcs:jce) = rlu  (jcs:jce,1)
    rlutcs(jcs:jce) = rlucs(jcs:jce,1)
    !
    rlnt  (jcs:jce) = -rlut  (jcs:jce)
    rlntcs(jcs:jce) = -rlutcs(jcs:jce)
    
    ! all half levels
    DO jk = 1, klevp1
      rln  (jcs:jce,jk) = (rld  (jcs:jce,jk) - rlu  (jcs:jce,jk))
      rlncs(jcs:jce,jk) = (rldcs(jcs:jce,jk) - rlucs(jcs:jce,jk))
    END DO
    !
    ! surface
    rlds  (jcs:jce) = rld  (jcs:jce,klevp1)
    rldscs(jcs:jce) = rldcs(jcs:jce,klevp1)
    !
    ! - adjust upward sfc longwave radiation for changed surface temperature
    drlus_dtsr(jcs:jce) = emiss(jcs:jce)*4._wp*stbo*tsr(jcs:jce)**3 ! derivative of rlus wrt. to tsr
    dtsr      (jcs:jce) = tsr(jcs:jce) - tsr_rt(jcs:jce)            ! change in tsr
    rlus(jcs:jce)       = rlu(jcs:jce,klevp1)                     & ! current rlus
      &                 + drlus_dtsr(jcs:jce) * dtsr(jcs:jce)
    !
    rlns  (jcs:jce) = rlds  (jcs:jce) - rlus  (jcs:jce)
    rlnscs(jcs:jce) = rldscs(jcs:jce) - rlus  (jcs:jce)

    
    ! Heating rates in atmosphere
    !----------------------------
    q_rsw(jcs:jce,1:klev) = rsn(jcs:jce,1:klev)-rsn(jcs:jce,2:klev+1)
    q_rlw(jcs:jce,1:klev) = rln(jcs:jce,1:klev)-rln(jcs:jce,2:klev+1)

  END SUBROUTINE radheating

END MODULE mo_radheating
