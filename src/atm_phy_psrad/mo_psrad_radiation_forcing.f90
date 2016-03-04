!>
!! Data types and variables used by the psrad radiation package.
!!
!! This module contains
!! <ol>
!! <li> definition of data types for organising psrad diagnostics,
!! <li> the actual variables that are declared of these types, and
!! <li> subroutines for (de-)allocating memory for the variables.
!! </ol>
!!
!! @author Sebastian Rast (MPI-M, 2016-03-04)
!!
!! @par Revision History
!! First version by Sebastian Rast (MPI-M, 2016-03-04)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_psrad_radiation_forcing
  USE mo_kind,                ONLY: wp
  USE mo_psrad_memory,        ONLY: psf=>prm_psrad_forcing
  USE mo_datetime,            ONLY: rdaylen
  USE mo_physical_constants,  ONLY: stbo
  USE mo_radiation_config,    ONLY: lradforcing
  USE mo_psrad_radiation_parameters,ONLY: cemiss

  IMPLICIT NONE
  PRIVATE

  PUBLIC          :: prepare_psrad_radiation_forcing, calculate_psrad_radiation_forcing

  CONTAINS

  !-----------------------------------------------------------------------------
  !>
  !! prepare_forcing: calculate intermediate quantities for forcing calculation
  !!
  !! @par Revision History
  !! original source by J.S.Rast (2016-03-4) based on code by M.A.Thomas and
  !!  S.J.Lorenz
  !-----------------------------------------------------------------------------
  SUBROUTINE prepare_psrad_radiation_forcing(jg                                ,&
       & kproma           ,kbdim             ,klevp1           ,krow           ,&
       & pflx_dnlw        ,pflx_dnsw         ,pflx_dnlw_clr    ,pflx_dnsw_clr   )

    INTEGER, INTENT(in)     :: jg,     & !< domain index 
                             & kproma, & !< block length
                             & kbdim,  & !< maximum block length
                             & klevp1, & !< number of layers + 1
                             & krow      !< block index
    REAL(wp), INTENT(in)    :: &
         pflx_dnlw(kbdim,klevp1),    & !< Net dwnwrd LW flux [Wm2]
         pflx_dnsw(kbdim,klevp1),    & !< Net dwnwrd SW flux [Wm2] for forcing
         pflx_dnlw_clr(kbdim,klevp1),& !< Net dn LW flux (clear sky) [Wm2]
         pflx_dnsw_clr(kbdim,klevp1)   !< Net dn SW flux (clear sky) [Wm2]

    psf(jg)%emter_for(1:kproma,1:klevp1,krow)=pflx_dnlw(1:kproma,1:klevp1)
    psf(jg)%emtef_for(1:kproma,1:klevp1,krow)=pflx_dnlw_clr(1:kproma,1:klevp1)
    psf(jg)%trsol_for(1:kproma,1:klevp1,krow)=pflx_dnsw(1:kproma,1:klevp1)
    psf(jg)%trsof_for(1:kproma,1:klevp1,krow)=pflx_dnsw_clr(1:kproma,1:klevp1)

  END SUBROUTINE prepare_psrad_radiation_forcing
  !-----------------------------------------------------------------------------
  !>
  !! calculate_forcing: calculate radiative forcing 
  !!
  !! @par Revision History
  !! original source by J.S.Rast (2010-04-16) based on code by M.A.Thomas and
  !!  S.J.Lorenz
  !-----------------------------------------------------------------------------
  SUBROUTINE calculate_psrad_radiation_forcing( jg,                &
       &  kproma             ,kbdim               ,klevp1          &
       & ,krow               ,psteplen             ,pi0                 ,pconvfact       &
       & ,pflxs              ,pflxs0              ,pflxt           &
       & ,pflxt0             ,pti                 ,pztsnew         )
    INTEGER, INTENT(in)          :: jg,     & !< domain index 
                                  & kproma, & !< block length
                                  & kbdim,  & !< maximum block length
                                  & klevp1, & !< number of layers + 1
                                  & krow      !< block index
    REAL(wp), INTENT(in)         :: &
         & psteplen,                   & !< time step length
         & pi0(kbdim),                 & !< solar irradiance at top of atmosphere
         & pconvfact(kbdim,klevp1-1),  & !< conversion factor radiation flux -> heating rate
         & pflxs(kbdim,klevp1),        & !< short wave net radiation flux per Watt irradiance (all sky)
         & pflxs0(kbdim,klevp1),       & !< short wave net radiation flux per Watt irradiance (clear sky)
         & pflxt(kbdim,klevp1),        & !< long wave net radiation flux per Watt irradiance (all sky)
         & pflxt0(kbdim,klevp1),       & !< long wave net radiation flux per Watt irradiance (clear sky)
         & pti(kbdim,klevp1),          & !< temperature at interfaces (atmosphere)
         & pztsnew(kbdim)                !< surface temperature

    INTEGER                     :: jk, jl, klev
    REAL(wp)                    :: zemter_for(kbdim,klevp1), zemtef_for(kbdim,klevp1), &
                                 & ztrsol_for(kbdim,klevp1), ztrsof_for(kbdim,klevp1)
    REAL(wp)                    :: zflxs_all_for(kbdim,klevp1), &
                                 & zflxs_clear_for(kbdim,klevp1), &
                                 & zflxt_all_for(kbdim,klevp1), &
                                 & zflxt_clear_for(kbdim,klevp1)
    REAL(wp)                    :: zdtdt_sw, zdtdt_all_for_sw, &
         zdtdt_lw, zdtdt_all_for_lw

    klev=klevp1-1
    zemter_for(1:kproma,1:klevp1)=psf(jg)%emter_for(1:kproma,1:klevp1,krow)
    zemtef_for(1:kproma,1:klevp1)=psf(jg)%emtef_for(1:kproma,1:klevp1,krow)
    ztrsol_for(1:kproma,1:klevp1)=psf(jg)%trsol_for(1:kproma,1:klevp1,krow)
    ztrsof_for(1:kproma,1:klevp1)=psf(jg)%trsof_for(1:kproma,1:klevp1,krow)

    IF (lradforcing(1)) THEN
      !--- Radiation flux forcing:
      zflxs_all_for(1:kproma,1:klevp1)=SPREAD(pi0(1:kproma),2,klevp1)* &
           ztrsol_for(1:kproma,1:klevp1)
      zflxs_clear_for(1:kproma,1:klevp1)=SPREAD(pi0(1:kproma),2,klevp1)* &
           ztrsof_for(1:kproma,1:klevp1)
      !  forcing solar wave length bands:
      DO jk = 1, klev
        DO jl = 1, kproma
          zdtdt_sw=pconvfact(jl,jk)*(pflxs(jl,jk+1)-pflxs(jl,jk))
          zdtdt_all_for_sw=pconvfact(jl,jk)*(zflxs_all_for(jl,jk+1)-zflxs_all_for(jl,jk))
          psf(jg)%netht_sw(jl,jk,krow) = psf(jg)%netht_sw(jl,jk,krow)+rdaylen*(zdtdt_sw-zdtdt_all_for_sw)*psteplen
        ENDDO
      END DO
      DO jk = 1, klevp1
        psf(jg)%d_aflx_sw(1:kproma,jk,krow) = psf(jg)%d_aflx_sw(1:kproma,jk,krow) + &
             (pflxs(1:kproma,jk)-zflxs_all_for(1:kproma,jk))*psteplen
        psf(jg)%d_aflx_swc(1:kproma,jk,krow) = psf(jg)%d_aflx_swc(1:kproma,jk,krow) + &
             (pflxs0(1:kproma,jk)-zflxs_clear_for(1:kproma,jk))*psteplen
      END DO
      psf(jg)%fsw_total_top(1:kproma,krow)=psf(jg)%d_aflx_sw(1:kproma,1,krow)
      psf(jg)%fsw_total_sur(1:kproma,krow)=psf(jg)%d_aflx_sw(1:kproma,klevp1,krow)
      psf(jg)%fsw_clear_top(1:kproma,krow)=psf(jg)%d_aflx_swc(1:kproma,1,krow)
      psf(jg)%fsw_clear_sur(1:kproma,krow)=psf(jg)%d_aflx_swc(1:kproma,klevp1,krow)
    END IF

    IF (lradforcing(2)) THEN
      !--- Radiation flux forcing:
      zflxt_all_for(1:kproma,1:klev)=zemter_for(1:kproma,1:klev)
      ! in the following formulae, the order of calculation has to be as is
      ! if not, for an aerosol free atmosphere the forcing is not exactly 0
      ! because of numeric effects (see the corresponding formulae in radheat)
      zflxt_all_for(1:kproma,klevp1)=zemter_for(1:kproma,klevp1)+ &
           cemiss*stbo*pti(1:kproma,klevp1)**4
      zflxt_all_for(1:kproma,klevp1)=zflxt_all_for(1:kproma,klevp1)- &
           cemiss*stbo*pztsnew(1:kproma)**4
      zflxt_clear_for(1:kproma,1:klev)=zemtef_for(1:kproma,1:klev)
      zflxt_clear_for(1:kproma,klevp1)=zemtef_for(1:kproma,klevp1) + &
           cemiss*stbo*pti(1:kproma,klevp1)**4
      zflxt_clear_for(1:kproma,klevp1)=zflxt_clear_for(1:kproma,klevp1)- &
           cemiss*stbo*pztsnew(1:kproma)**4
      DO jk = 1,klev
        DO jl = 1, kproma
          zdtdt_lw=pconvfact(jl,jk)*(pflxt(jl,jk+1)-pflxt(jl,jk))
          zdtdt_all_for_lw=pconvfact(jl,jk)*(zflxt_all_for(jl,jk+1)-zflxt_all_for(jl,jk))
          psf(jg)%netht_lw(jl,jk,krow) = psf(jg)%netht_lw(jl,jk,krow)+rdaylen*(zdtdt_lw-zdtdt_all_for_lw)*psteplen
        END DO
      END DO
      DO jk = 1, klevp1
        psf(jg)%d_aflx_lw(1:kproma,jk,krow) = psf(jg)%d_aflx_lw(1:kproma,jk,krow) + &
             (pflxt(1:kproma,jk)-zflxt_all_for(1:kproma,jk))*psteplen
        psf(jg)%d_aflx_lwc(1:kproma,jk,krow) = psf(jg)%d_aflx_lwc(1:kproma,jk,krow) + &
             (pflxt0(1:kproma,jk)-zflxt_clear_for(1:kproma,jk))*psteplen
      END DO
      psf(jg)%flw_total_top(1:kproma,krow) = psf(jg)%d_aflx_lw(1:kproma,1,krow)
      psf(jg)%flw_total_sur(1:kproma,krow) = psf(jg)%d_aflx_lw(1:kproma,klevp1,krow)
      psf(jg)%flw_clear_top(1:kproma,krow) = psf(jg)%d_aflx_lwc(1:kproma,1,krow)
      psf(jg)%flw_clear_sur(1:kproma,krow) = psf(jg)%d_aflx_lwc(1:kproma,klevp1,krow)
    END IF
  END SUBROUTINE calculate_psrad_radiation_forcing

END MODULE mo_psrad_radiation_forcing
