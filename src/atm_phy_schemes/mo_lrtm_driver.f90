!>
!! @brief Module to provide interface to rrtmg longwave radiation
!!
!! @remarks
!!   This module contains routines that provide the interface between ECHAM
!!   and the AER RRTMG radiation code.  Mostly it organizes and calculates the
!!   information necessary to call the radiative transfer solvers.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-11-29)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Major segments of this code combines and rewrites (for the ICON standard)
!!   code previously provided by AER and copyrighed by them.  The authors of the
!!   original AER code are: Eli J. Mlawer, Jennifer S. Delamere, Michael J.
!!   Iacono and Shepard A. Clough with acknowledgments to Steven J. Taubman,
!!   Karen Cady-Pereira, Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert
!!   Bergstrom. The rewrites were designed to better interface with the structure
!!   of the ICON family of models and elements of the ICON programming standard.
!!
!! @par Copyright
!!   2002-2009 by the Deutsche Wetterdienst (DWD) and the Max-Planck-Institut
!!   for Meteorology (MPI-M).  This software is provided for non-commerical
!!   use only.  See the LICENSE and the WARRANTY conditions.  The AER copyright
!!   on the original code is as follows: Copyright 2002-2009, Atmospheric and
!!   Environmental Research, Inc. (AER). This software may be used, copied, or
!!   redistributed as long as it is not sold and this copyright notice is
!!   reproduced on each copy made.  This model is provided as is without any
!!   express or implied warranties. (http://www.rtweb.aer.com/)
!!
!! @par License
!!   The use of ICON is hereby granted free of charge for an unlimited time,
!!   provided:
!!   <ol>
!!    <li> Its use is limited to own non-commercial and non-violent purposes;
!!    <li> The code is not re-distributed without the consent of DWD and MPI-M;
!!    <li> This header appears in all copies of the code;
!!    <li> You accept the warranty conditions (see WARRANTY).
!!   </ol>
!!   Commericial use of the code is allowed subject to a separate licensing
!!   agreement with the DWD and MPI-M
!!
!! @par Warranty
!!   This code is distributed in the home that it will be useful, but WITHOUT
!!   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!!   FITNESS FOR A PARTICULAR PURPOSE.
!
MODULE mo_lrtm

  USE mo_kind,         ONLY: wp

#ifdef __ICON__
  USE mo_physical_constants, ONLY: amd, amw, grav
#else
  USE mo_constants,          ONLY: amd, amw, grav=>g
#endif

  USE mo_lrtm_par,     ONLY: nbndlw, nmol, ngptlw, ngb
  USE mo_lrtm_rtrnmr,  ONLY: lrtm_rtrnmr
  USE mo_lrtm_coeffs,  ONLY: lrtm_coeffs
  USE mo_lrtm_taumol,  ONLY: lrtm_taumol

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lrtm

CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! @brief Prepares information for radiation call
  !!
  !! @remarks: This program is the driver subroutine for RRTMG_LW, the AER LW
  !! radiation model for application to GCMs, that has been adapted from RRTM_LW
  !! for improved efficiency.
  !! This routine:
  !!    1) calls INATM to read in the atmospheric profile from GCM;
  !!       all layering in RRTMG is ordered from surface to toa.
  !!    2) calls COEFFS to calculate various quantities needed for
  !!       the radiative transfer algorithm
  !!    3) calls TAUMOL to calculate gaseous optical depths for each
  !!       of the 16 spectral bands
  !!    4) calls RTRNMR (for both clear and cloudy profiles) to perform the
  !!       radiative transfer calculation with a maximum-random cloud
  !!       overlap method, or calls RTRN to use random cloud overlap.
  !!    5) passes the necessary fluxes and cooling rates back to GCM
  !!
  !! The optional calculation of the change in upward flux as a function of
  !! surface temperature is available (controlled by input flag idrv).  This
  !! can be utilized  to approximate adjustments to the upward flux profile
  !! caused only by a change in surface temperature between full radiation
  !! calls.  This feature uses the pre-calculated derivative of the Planck
  !! function with respect to surface temperature. To calculate these
  !! derivatives idrv must be set to 1.
  !
  SUBROUTINE lrtm     (                            &
       & kproma     ,klev     ,play     ,psfc     ,&
       & tlay       ,tlev     ,tsfc     ,wkl_2d   ,&
       & wx_2d      ,coldry_2d,emis     ,cldfr    ,&
       & taucld     ,tauaer   ,uflx     ,dflx     ,&
       & uflxc      ,dflxc    ,duflx_dt ,duflxc_dt)

    INTEGER, INTENT(in) :: &
         & kproma,        & !< Number of horizontal columns
         & klev             !< Number of model layers

    REAL(wp), INTENT(in) :: &
         & play(:,:),     & !< Layer pressures [hPa, mb] (kbdim,klev)
         & psfc(:),       & !< Surface pressure [hPa, mb] (kbdim)
         & tlay(:,:),     & !< Layer temperatures [K] (kbdim,klev)
         & tlev(:,:),     & !< Interface temperatures [K] (kbdim,klev+1)
         & tsfc(:),       & !< Surface temperature [K] (kbdim)
         & wkl_2d(:,:,:), & !< Gas volume mixing ratios
         & wx_2d(:,:,:),  & !< CFC type gas volume mixing ratios
         & coldry_2d(:,:),& !< Column dry amount
         & emis(:,:),     & !< Surface emissivity  (kbdim,nbndlw)
         & cldfr(:,:),    & !< Cloud fraction  (kbdim,klev)
         & taucld(:,:,:), & !< Coud optical depth (kbdim,klev,nbndlw)
         & tauaer(:,:,:)    !< Qerosol optical depth (kbdim,klev,nbndlw)

    REAL(wp), INTENT(out) :: &
         uflx(:,:) , & !< Tot sky longwave upward flux [W/m2], (kbdim,klev+1)
         dflx(:,:) , & !< Tot sky longwave downward flux [W/m2], (kbdim,klev+1)
         uflxc(:,:), & !< Clr sky longwave upward flux [W/m2], (kbdim,klev+1)
         dflxc(:,:)    !< Clr sky longwave downward flux [W/m2], (kbdim,klev+1)

    REAL(wp), INTENT(out), OPTIONAL :: &
         & duflx_dt(:,:), & !< change in uflx wrt Temp [w/m2/k] (kbdim/klev)
         & duflxc_dt(:,:)   !< change in uflxc wrt Temp [w/m2/k] (kbdim/klev)

    INTEGER :: ncbands(kproma)   !< number of cloud spectral bands

    REAL(wp) ::               &
         & tz(kproma,0:klev),           & !< level (interface) temperatures [K]
         & wbrodl(kproma,klev),         & !< broadening gas column density (mol/cm2)
         & pwvcm(kproma),               & !< precipitable water vapor [cm]
         & fracs(kproma,klev,ngptlw),   & !< layer cloud fraction
         & taug(kproma,klev,ngptlw),    & !< gas optical depth
         & taut(kproma,klev,ngptlw),    & !< gaseous + aerosol optical depths
         & taucloud(kproma,klev,nbndlw),& !< layer in-cloud optical depth
         & totuflux(kproma,0:klev),     & !< upward longwave flux (w/m2)
         & totdflux(kproma,0:klev),     & !< downward longwave flux (w/m2)
         & fnet(kproma,0:klev),         & !< net longwave flux (w/m2)
         & totuclfl(kproma,0:klev),     & !< clear sky upward longwave flux (w/m2)
         & totdclfl(kproma,0:klev),     & !< clear sky downward longwave flux (w/m2)
         & fnetc(kproma,0:klev),        & !< clear sky net longwave flux (w/m2)
         & dtotuflux_dt(kproma,0:klev), & !< change in uflxx (w/m2/k) wrt Temp
         & dtotuclfl_dt(kproma,0:klev), & !< change in cuflx  (w/m2/k) wrt Temp
         & tauctot(kproma,klev)             !< band integrated cloud optical depth

    INTEGER ::            &
         & laytrop(kproma),             & !< tropopause layer index
         & jp(kproma,klev),             & !< lookup table index
         & jt(kproma,klev),             & !< lookup table index
         & jt1(kproma,klev),            & !< lookup table index
         & indself(kproma,klev),        &
         & indfor(kproma,klev),         &
         & indminor(kproma,klev)

    REAL(wp) ::                  &
         & planklay(kproma,klev,nbndlw),   & !
         & planklev(kproma,0:klev,nbndlw), & !
         & plankbnd(kproma,nbndlw),        & !
         & dplankbnd_dt(kproma,nbndlw),    & !
         & colh2o(kproma,klev),            & !< column amount (h2o)
         & colco2(kproma,klev),            & !< column amount (co2)
         & colo3(kproma,klev),             & !< column amount (o3)
         & coln2o(kproma,klev),            & !< column amount (n2o)
         & colco(kproma,klev),             & !< column amount (co)
         & colch4(kproma,klev),            & !< column amount (ch4)
         & colo2(kproma,klev),             & !< column amount (o2)
         & colbrd(kproma,klev),            & !< column amount (broadening gases)
         & selffac(kproma,klev),           &
         & selffrac(kproma,klev),          &
         & forfac(kproma,klev),            &
         & forfrac(kproma,klev),           &
         & minorfrac(kproma,klev),         &
         & scaleminor(kproma,klev),        &
         & scaleminorn2(kproma,klev)

    REAL(wp) ::      &
         fac00(kproma,klev),        &
         fac01(kproma,klev),        &
         fac10(kproma,klev),        &
         fac11(kproma,klev),        &
         rat_h2oco2(kproma,klev),   &
         rat_h2oco2_1(kproma,klev), &
         rat_h2oo3(kproma,klev),    &
         rat_h2oo3_1(kproma,klev),  &
         rat_h2on2o(kproma,klev),   &
         rat_h2on2o_1(kproma,klev), &
         rat_h2och4(kproma,klev),   &
         rat_h2och4_1(kproma,klev), &
         rat_n2oco2(kproma,klev),   &
         rat_n2oco2_1(kproma,klev), &
         rat_o3co2(kproma,klev),    &
         rat_o3co2_1(kproma,klev)

    INTEGER :: jl, jk, ib, ig ! loop indicies

    INTEGER, PARAMETER :: idrv   =  0 !< Flag for calculating derivs
    INTEGER, PARAMETER :: iout   =  0 !< option output flag
    INTEGER, PARAMETER :: iend   = 16 !< last band
    INTEGER, PARAMETER :: istart =  1 !< first band
    REAL(wp),    PARAMETER :: cldmin = 1.e-20_wp ! minimum val for clouds

    !
    ! 1.0 Convert 2D data into 1D column data an prepare some auxilliary info
    ! --------------------------------
    !
    CALL inatm(                                          &
         & kproma       ,klev           ,psfc           ,&
         & tlev         ,wkl_2d         ,&
         & coldry_2d    ,&
         & tz           ,&
         & wbrodl       ,pwvcm)

    DO jl = 1, kproma  ! loop over columns
      ncbands(jl)      = 1
      tauctot(jl,:)    = 0.0_wp
      taucloud(jl,:,:) = 0.0_wp
    ENDDO

    DO ib = 1,nbndlw
      DO jk = 1, klev
        DO jl = 1, kproma  ! loop over columns
          tauctot(jl,jk)  = tauctot(jl,jk) + taucld(jl,jk,ib)
        END DO
      END DO
    ENDDO

    DO jk = 1, klev
      DO jl = 1, kproma  ! loop over columns
        IF (cldfr(jl,jk) .GE. cldmin .AND. tauctot(jl,jk) .GE. cldmin) THEN
          ncbands(jl) = 16
!CDIR EXPAND=nbndlw
          DO ib = 1,nbndlw
            taucloud(jl,jk,ib) = taucld(jl,jk,ib)
          END DO
        END IF
      END DO
    ENDDO

    !
    ! 2.0  Calculate information needed by the radiative transfer routine
    ! that is specific to this atmosphere, especially some of the
    ! coefficients and indices needed to compute the optical depths
    ! by interpolating data from stored reference atmospheres.
    ! --------------------------------
    !
    CALL lrtm_coeffs(                                                &
         & kproma, klev ,istart       ,play         ,tlay           ,&
         & tz           ,tsfc         ,emis         ,coldry_2d      ,&
         & wkl_2d       ,wbrodl       ,laytrop      ,jp             ,&
         & jt           ,jt1          ,planklay     ,planklev       ,&
         & plankbnd     ,idrv         ,dplankbnd_dt ,colh2o         ,&
         & colco2       ,colo3        ,coln2o       ,colco          ,&
         & colch4       ,colo2        ,colbrd       ,fac00          ,&
         & fac01        ,fac10        ,fac11        ,rat_h2oco2     ,&
         & rat_h2oco2_1 ,rat_h2oo3    ,rat_h2oo3_1  ,rat_h2on2o     ,&
         & rat_h2on2o_1 ,rat_h2och4   ,rat_h2och4_1 ,rat_n2oco2     ,&
         & rat_n2oco2_1 ,rat_o3co2    ,rat_o3co2_1  ,selffac        ,&
         & selffrac     ,indself      ,forfac       ,forfrac        ,&
         & indfor       ,minorfrac    ,scaleminor   ,scaleminorn2   ,&
         & indminor     )

    !
    !  3.0 Calculate the gaseous optical depths and Planck fractions for
    !  each longwave spectral band.
    ! --------------------------------
    !
    CALL lrtm_taumol(                                                &
         & kproma, klev ,play         ,wx_2d        ,coldry_2d      ,&
         & laytrop      ,jp           ,jt           ,jt1            ,&
         & colh2o       ,colco2       ,colo3        ,coln2o         ,&
         & colco        ,colch4       ,colo2        ,colbrd         ,&
         & fac00        ,fac01        ,fac10        ,fac11          ,&
         & rat_h2oco2   ,rat_h2oco2_1 ,rat_h2oo3    ,rat_h2oo3_1    ,&
         & rat_h2on2o   ,rat_h2on2o_1 ,rat_h2och4   ,rat_h2och4_1   ,&
         & rat_n2oco2   ,rat_n2oco2_1 ,rat_o3co2    ,rat_o3co2_1    ,&
         & selffac      ,selffrac     ,indself      ,forfac         ,&
         & forfrac      ,indfor       ,minorfrac    ,scaleminor     ,&
         & scaleminorn2 ,indminor     ,fracs        ,taug)
    !
    ! --- Combine gaseous and aerosol optical depths, if aerosol active
    !
    DO ig = 1, ngptlw
      DO jk = 1, klev
        DO jl = 1, kproma  ! loop over columns
          taut(jl,jk,ig) = taug(jl,jk,ig) + tauaer(jl,jk,ngb(ig))
        ENDDO
      ENDDO
    ENDDO
    !
    ! 4.0 Call the radiative transfer routine. (random-maximum overlap)
    ! --------------------------------
    !
    CALL lrtm_rtrnmr(                                                &
         & kproma, klev ,istart       ,iend         ,iout           ,&
         & emis         ,ncbands      ,cldfr        ,taucloud       ,&
         & planklay     ,planklev     ,plankbnd     ,pwvcm          ,&
         & fracs        ,taut         ,totuflux     ,totdflux       ,&
         & fnet         ,totuclfl     ,totdclfl     ,fnetc          ,&
         & idrv         ,dplankbnd_dt ,dtotuflux_dt ,dtotuclfl_dt   )
    !
    ! 5.0 Finalize output
    ! --------------------------------
    !
    DO jk = 0, klev
      DO jl = 1, kproma  ! loop over columns
        uflx(jl,jk+1)  = totuflux(jl,jk)
        dflx(jl,jk+1)  = totdflux(jl,jk)
        uflxc(jl,jk+1) = totuclfl(jl,jk)
        dflxc(jl,jk+1) = totdclfl(jl,jk)
      ENDDO
    ENDDO

    IF (idrv .EQ. 1) THEN
      DO jk = 0, klev
        DO jl = 1, kproma  ! loop over columns
          duflx_dt(jl,jk+1) = dtotuflux_dt(jl,jk)
          duflxc_dt(jl,jk+1) = dtotuclfl_dt(jl,jk)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE lrtm
  !-----------------------------------------------------------------------------
  !>
  !! @brief Copies 2D column height input data into 1D column data
  !!
  SUBROUTINE inatm (                                 &
       & kproma       ,klev         ,psfc           ,&
       & tlev         ,wkl_2d         ,&
       & coldry_2d    ,&
       & tz           ,&
       & wbrodl       ,pwvcm)

    INTEGER, INTENT(in) :: &
         & kproma,       & !< Number of grid points
         & klev            !< Number of model layers

    REAL(wp), INTENT(in) :: &
         & psfc(:),       & !< Surface pressure [hPa, mb] (kbdim)
         & tlev(:,:),     & !< Interface temperatures [K] (kbdim,klev+1)
         & wkl_2d(:,:,:), & !< Gas volume mixing ratios weighted by coldry
         & coldry_2d(:,:)   !< Column dry amount

    REAL(wp), INTENT(out) ::  &
!        & tz(:,0:),      & !< level (interface) temperatures [K] (kproma,0:klev)
         & tz(kproma,0:klev), & !< level (interface) temperatures [K]
         & wbrodl(:,:),   & !< broadening gas column density (mol/cm2)
         & pwvcm(:)         !< precipitable water vapor (cm)

    INTEGER :: jl, jk, imol ! Loop indices
    REAL(wp), DIMENSION(kproma) :: amttl, wvttl, wvsh, summol

!    PRINT *,'LBOUND(tz)',LBOUND(tz) !T.R.
!    PRINT *,'UBOUND(tz)',UBOUND(tz) !T.R.
!    STOP'Stop.' !T.R.

    !
    !  --- Initialize all molecular amounts and cloud properties to zero
    !
    amttl(:)  = 0.0_wp
    wvttl(:)  = 0.0_wp
    !
    !  --- Remap data to 1D column arrays
    !
    tz(1:kproma,0:klev) = tlev(1:kproma,1:klev+1)
    !
    !  --- Calculate auxillary information
    !
    DO jk = 1, klev
      summol(:) = 0.0_wp
      DO jl = 1, kproma
!CDIR EXPAND=nmol
        DO imol = 2, nmol
          summol(jl) = summol(jl) + wkl_2d(jl,imol,jk)
        ENDDO
        wbrodl(jl,jk) = coldry_2d(jl,jk) - summol(jl)
        amttl(jl) = amttl(jl) + coldry_2d(jl,jk)+wkl_2d(jl,1,jk)
        wvttl(jl) = wvttl(jl) + wkl_2d(jl,1,jk)
      ENDDO
    ENDDO

    DO jl = 1, kproma
      wvsh(jl) = (amw * wvttl(jl)) / (amd * amttl(jl))
      pwvcm(jl) = wvsh(jl) * (1.e3_wp * psfc(jl)) / (1.e2_wp * grav)
    ENDDO

  END SUBROUTINE inatm

END MODULE mo_lrtm

