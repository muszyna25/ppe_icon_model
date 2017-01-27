!>
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (2011-08)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_surface_diag

  USE mo_kind,              ONLY: wp
  USE mo_physical_constants,ONLY: grav, als, alv, vtmpc2, cpd, rdv, tmelt,    &
                                  vtmpc1, vtmpc2, rd
  USE mo_echam_convect_tables,     ONLY: lookup_ua_list_spline
  USE mo_echam_vdiff_params,ONLY: tpfac2
  USE mo_echam_phy_memory,  ONLY: cdimissval

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wind_stress, surface_fluxes, nsurf_diag

CONTAINS
  !>
  !!
  !!
  SUBROUTINE surface_fluxes( lsfc_heat_flux, psteplen,             &! in
                           & kproma, kbdim, ksfc_type,             &! in
                           & idx_wtr, idx_ice, idx_lnd, ih, iqv,   &! in
                           & pfrc, pcfh_tile, pfac_sfc,            &! in
                           & pcptv_tile, pqsat_tile,               &! in
                           & pca, pcs, bb_btm,                     &! in
                           & plhflx_gbm, pshflx_gbm,               &! out
                           & pevap_gbm,                            &! out
                           & plhflx_tile, pshflx_tile,             &! inout for JSBACH land
                           & pevap_tile,                           &! out
                           & evapotranspiration)

    LOGICAL, INTENT(IN) :: lsfc_heat_flux
    REAL(wp),INTENT(IN) :: psteplen
    INTEGER, INTENT(IN) :: kproma, kbdim, ksfc_type
    INTEGER, INTENT(IN) :: idx_wtr, idx_ice, idx_lnd
    INTEGER, INTENT(IN) :: ih, iqv

    REAL(wp),INTENT(IN) :: pfrc      (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pcfh_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pfac_sfc  (kbdim)
    REAL(wp),INTENT(IN) :: pcptv_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pqsat_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pca       (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pcs       (kbdim,ksfc_type)

    REAL(wp),INTENT(IN) :: bb_btm(kbdim,ksfc_type,ih:iqv)

    REAL(wp),INTENT(OUT) :: plhflx_gbm(kbdim)
    REAL(wp),INTENT(OUT) :: pshflx_gbm(kbdim)
    REAL(wp),INTENT(OUT) :: pevap_gbm(kbdim)

    REAL(wp),INTENT(INOUT) :: plhflx_tile(kbdim,ksfc_type) ! INOUT due to JSBACH land
    REAL(wp),INTENT(INOUT) :: pshflx_tile(kbdim,ksfc_type) ! INOUT due to JSBACH land
    REAL(wp),INTENT(OUT)   :: pevap_tile(kbdim,ksfc_type)

    REAL(wp),OPTIONAL,INTENT(IN)    :: evapotranspiration(kbdim) ! present for JSBACH land

    INTEGER  :: jsfc
    REAL(wp) :: zconst, zdqv(kbdim), zdcptv(kbdim)

    !===================================================================
    ! If surface heat fluxes (incl. latent) are switched off, set
    ! corresponding values to zero and return to the calling subroutine.
    !===================================================================

    DO jsfc = 1,ksfc_type
      IF (jsfc /= idx_lnd .OR. .NOT. PRESENT(evapotranspiration)) THEN  ! for JSBACH land, the fluxes are already in these arrays
        plhflx_tile(1:kproma,jsfc) = 0._wp
        pshflx_tile(1:kproma,jsfc) = 0._wp
      END IF
    END DO

    IF (.NOT.lsfc_heat_flux) THEN
      RETURN
    END IF

    !===================================================================
    ! Otherwise compute diagnostics
    !===================================================================
    zconst = 1._wp/(grav*psteplen)

    !-------------------------------------------------------------------
    ! Moisture fluxes (aka evaporation rates)
    !-------------------------------------------------------------------
    ! Instantaneous moisture flux on each tile

    pevap_tile(1:kbdim,:) = 0._wp

    DO jsfc = 1,ksfc_type

      ! Vertical gradient of specific humidity scaled by factor (1/tpfac1).
      ! Formula: ( qv_{tavg,klev} - qs_tile )/tpfac1
      ! Here qv_{tavg,klev} = tpfac1*qv_klev(t+dt) + (1-tpfac1)*qv_klev(t)
      !                     = tpfac1*bb_qv
      ! where bb_qv is the solution of the linear system at the lowest
      ! model level (i.e., the full level right above surface).
      ! (Formula translated from ECHAM. Question: why using the blended
      ! humidity, not the value on individual surface?)
      ! bb was replaced by bb_btm (according to E. Roeckner), now not blended
      ! quantity used.

      zdqv(1:kproma) =   bb_btm(1:kproma,jsfc,iqv)*pca(1:kproma,jsfc)          &
                     & - tpfac2*pqsat_tile(1:kproma,jsfc)*pcs(1:kproma,jsfc)

      ! Moisture flux ( = evaporation). Formula:
      ! (g*psteplen)**(-1)*[  tpfac1*g*psteplen*(air density)*(exchange coef)
      !                     *(tpfac1)**(-1)*( qv_{tavg,klev} - qs_tile ) ]

      IF (jsfc == idx_lnd) THEN
        pevap_tile(1:kproma,jsfc) = evapotranspiration(1:kproma)
      ELSE
        pevap_tile(1:kproma,jsfc) =  zconst*pfac_sfc(1:kproma) &
                                  & *pcfh_tile(1:kproma,jsfc)  &
                                  & *zdqv(1:kproma)
      END IF
    ENDDO

    ! Compute grid box mean and time integral
    ! The instantaneous grid box mean moisture flux will be passed on
    ! to the cumulus convection scheme.

    pevap_gbm(1:kbdim) = 0._wp   ! "pqhfla" in echam

    DO jsfc = 1,ksfc_type
      pevap_gbm(1:kproma) = pevap_gbm(1:kproma)                           &
        &                 + pfrc(1:kproma,jsfc)*pevap_tile(1:kproma,jsfc)
    ENDDO

    !-------------------------------------------------------------------
    ! Latent heat flux
    !-------------------------------------------------------------------
    ! Instantaneous values

    ! Note: latent heat flux from land is already in plhflx_tile(:,idx_lnd)

    IF (idx_ice<=ksfc_type) &
    plhflx_tile(1:kproma,idx_ice) = als*pevap_tile(1:kproma,idx_ice)

    IF (idx_wtr<=ksfc_type) &
    plhflx_tile(1:kproma,idx_wtr) = alv*pevap_tile(1:kproma,idx_wtr)

    ! Accumulated grid box mean

    plhflx_gbm(1:kbdim) = 0.0_wp

    DO jsfc = 1,ksfc_type
      plhflx_gbm(1:kproma) = plhflx_gbm(1:kproma)                           &
        &                  + pfrc(1:kproma,jsfc)*plhflx_tile(1:kproma,jsfc)
    ENDDO

    !-------------------------------------------------------------------
    ! Sensible heat flux
    !-------------------------------------------------------------------
    ! Instantaneous flux on each tile

    DO jsfc = 1,ksfc_type

      ! Note: sensible heat flux from land is already in pshflx_tile(:,idx_lnd)
      IF (jsfc == idx_lnd) THEN
        ! COMMENT: already done at begin of routine
      ELSE

      ! Vertical gradient of dry static energy.
      ! (Formula translated from ECHAM. Question: why using the blended
      ! dry static energy, not the value on individual surface?)
      ! bb was replaced by bb_btm (according to E. Roeckner), now not blended
      ! quantity used.

        zdcptv(1:kproma) = bb_btm(1:kproma,jsfc,ih) - tpfac2*pcptv_tile(1:kproma,jsfc)

      ! Flux of dry static energy

        pshflx_tile(1:kproma,jsfc) =  zconst*pfac_sfc(1:kproma) &
                                   & *pcfh_tile(1:kproma,jsfc)  &
                                   & *zdcptv(1:kproma)

      END IF

    ENDDO


    ! grid box mean

    pshflx_gbm(1:kbdim) = 0.0_wp

    DO jsfc = 1,ksfc_type
      pshflx_gbm(1:kproma) = pshflx_gbm(1:kproma)                           &
        &                  + pfrc(1:kproma,jsfc)*pshflx_tile(1:kproma,jsfc)
    ENDDO

  END SUBROUTINE surface_fluxes
  !-------------
  !!
  !! Compute wind stress over each surface type
  !!
  SUBROUTINE wind_stress( lsfc_mom_flux, psteplen,              &! in
                        & kproma, kbdim, ksfc_type,             &! in
                        & pfrc, pcfm_tile, pfac_sfc,            &! in
                        & pu_rtpfac1, pv_rtpfac1,               &! in
                        & pu_stress_gbm,  pv_stress_gbm,        &! out
                        & pu_stress_tile, pv_stress_tile        )! out

    LOGICAL, INTENT(IN)    :: lsfc_mom_flux
    REAL(wp),INTENT(IN)    :: psteplen
    INTEGER, INTENT(IN)    :: kproma, kbdim, ksfc_type

    REAL(wp),INTENT(IN)    :: pfrc            (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pcfm_tile       (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pfac_sfc        (kbdim)
    REAL(wp),INTENT(IN)    :: pu_rtpfac1      (kbdim)
    REAL(wp),INTENT(IN)    :: pv_rtpfac1      (kbdim)
    REAL(wp),INTENT(OUT)   :: pu_stress_gbm   (kbdim)
    REAL(wp),INTENT(OUT)   :: pv_stress_gbm   (kbdim)
    REAL(wp),INTENT(OUT)   :: pu_stress_tile  (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pv_stress_tile  (kbdim,ksfc_type)

    INTEGER  :: jsfc
    REAL(wp) :: zconst
    ! Local variables

    INTEGER  :: loidx  (kbdim,ksfc_type) !< counter for masks
    INTEGER  :: is     (ksfc_type)       !< counter for masks
    INTEGER  :: jls, jl, js

    !===================================================================
    ! If surface momentum fluxes is switched off (i.e., using free slip
    ! boundary condition), set wind stress to zero and return to the
    ! calling subroutine.
    !===================================================================
    IF (.NOT.lsfc_mom_flux) THEN

      pu_stress_tile(1:kbdim,:) = 0._wp
      pv_stress_tile(1:kbdim,:) = 0._wp

    !===================================================================
    ! Otherwise do computation
    !===================================================================
    ELSE
      zconst = 1._wp/(grav*psteplen)

      ! Compute wind stress over each surface type, then accumulate
      ! grid box mean. Formula for wind stress:
      !   (grav*psteplen)**(-1)
      !  *[grav*psteplen*tpfac1*(air density)]
      !  *(surface turbulent exchange coeff)
      !  *[(u-/v-wind at lowest model level)/tpfac1]

      pu_stress_tile(1:kbdim,:) = cdimissval
      pv_stress_tile(1:kbdim,:) = cdimissval

      pu_stress_gbm(1:kbdim) = 0.0_wp
      pv_stress_gbm(1:kbdim) = 0.0_wp

      ! check for masks
      !
      DO jsfc = 1,ksfc_type
        is(jsfc) = 0
        DO jl = 1,kproma
          IF(pfrc(jl,jsfc).GT.0.0_wp) THEN
            is(jsfc) = is(jsfc) + 1
            loidx(is(jsfc),jsfc) = jl
          ENDIF
        ENDDO
      ENDDO

      DO jsfc = 1,ksfc_type
        DO jls = 1,is(jsfc)
          ! set index
          js=loidx(jls,jsfc)

          pu_stress_tile(js,jsfc) = zconst*pfac_sfc(js) *pcfm_tile(js,jsfc) &
                                        & *pu_rtpfac1(js)

          pv_stress_tile(js,jsfc) = zconst*pfac_sfc(js) *pcfm_tile(js,jsfc) &
                                        & *pv_rtpfac1(js)
        END DO
      END DO

      DO jsfc = 1,ksfc_type
        DO jls = 1,is(jsfc)
          ! set index
          js=loidx(jls,jsfc)

          pu_stress_gbm(js)       = pu_stress_gbm(js) + pu_stress_tile(js,jsfc) &
                                        &  *pfrc(js,jsfc)

          pv_stress_gbm(js)       = pv_stress_gbm(js) + pv_stress_tile(js,jsfc) &
                                        &  *pfrc(js,jsfc)
        END DO
      END DO

    END IF ! lsfc_mom_flux

  END SUBROUTINE wind_stress
  !-------------
  !!
  !! Compute diagnostics: 10m wind, u and v in 10m,
  !!                      temperature in 2m, dew point temperature in 2m
  !!
  SUBROUTINE nsurf_diag( kproma, kbdim, ksfc_type,        &! in
                       & idx_lnd,                         &! in
                       & pfrc,                            &! in
                       & pqm1,                            &! in humidity 
                       & ptm1,                            &
                       & papm1,     paphm1,               &
                       & pxm1,                            &
                       & pum1,      pvm1,                 &
                       & pocu,      pocv,                 & 
                       & pgeom1,                          &! in geopotential above surface
                       & pcptgz,                          &! in dry static energy
                       & pcpt_tile,                       &! in dry static energy 
                       & pbn_tile,                        &! in for diagnostic
                       & pbhn_tile,                       &! in for diagnostic
                       & pbh_tile,                        &! in for diagnostic
                       & pbm_tile,                        &! in for diagnostic
                       & pri_tile,                        &! in moist Richardson number
                       & psfcWind_gbm,                    &! out 10m windspeed 
                       & ptas_gbm,                        &! out temperature in 2m
                       & pdew2_gbm,                       &! out dew point temperature in 2m
                       & puas_gbm,                        &! out zonal wind in 10m
                       & pvas_gbm,                        &! out meridional wind in 10m
                       & ptasmax,                         &! inout max 2m temperature
                       & ptasmin,                         &! inout min 2m temperature
                       & psfcWind_tile,                   &! out 10m windspeed 
                       & ptas_tile,                       &! out temperature in 2m
                       & pdew2_tile,                      &! out dew point temperature in 2m
                       & puas_tile,                       &! out zonal wind in 10m
                       & pvas_tile                        )! out meridional wind in 10m

    INTEGER, INTENT(IN) :: kproma, kbdim, ksfc_type
    INTEGER, INTENT(IN) :: idx_lnd

    REAL(wp),INTENT(IN) :: pfrc     (kbdim,ksfc_type) !< fraction of the grid box occupied by
                                                      !< each surface type
    REAL(wp), INTENT(in)     :: pqm1(kbdim), pgeom1(kbdim) 
    REAL(wp), INTENT(in)     :: pcptgz(kbdim)         !< dry static energy at surface level
    REAL(wp), INTENT(in)     :: pcpt_tile(kbdim,ksfc_type) !< dry static energy on tiles
    REAL(wp), TARGET, INTENT(in) :: pbn_tile(kbdim,ksfc_type)   !< for diagnostics
    REAL(wp), TARGET, INTENT(in) :: pbhn_tile(kbdim,ksfc_type)  !< for diagnostics
    REAL(wp), INTENT(in)     :: pbh_tile(kbdim,ksfc_type)   !< for diagnostics
    REAL(wp), INTENT(in)     :: pbm_tile(kbdim,ksfc_type)   !< for diagnostics
    REAL(wp), INTENT(in)     :: pri_tile   (kbdim,ksfc_type) !< moist Richardson number
    REAL(wp), INTENT(in)     :: ptm1(kbdim), papm1(kbdim), pxm1(kbdim)
    REAL(wp), INTENT(in)     :: pum1(kbdim), pvm1(kbdim), paphm1(kbdim) ! =paphm1(kbdim, klevp1)
    REAL(wp), INTENT(in)     :: pocu(kbdim), pocv(kbdim)
    REAL(wp), INTENT(out)    :: psfcWind_gbm(kbdim), psfcWind_tile(kbdim,ksfc_type)
    REAL(wp), INTENT(out)    :: ptas_gbm(kbdim),     ptas_tile(kbdim,ksfc_type)
    REAL(wp), INTENT(out)    :: pdew2_gbm(kbdim),    pdew2_tile(kbdim,ksfc_type)
    REAL(wp), INTENT(out)    :: puas_gbm(kbdim),     puas_tile(kbdim,ksfc_type)
    REAL(wp), INTENT(out)    :: pvas_gbm(kbdim),     pvas_tile(kbdim,ksfc_type)
    REAL(wp), INTENT(inout)  :: ptasmax(kbdim),      ptasmin(kbdim)

    ! Local variables
    
    INTEGER  :: loidx  (kbdim,ksfc_type) !< counter for masks
    INTEGER  :: is     (ksfc_type)       !< counter for masks
    INTEGER  :: jls, jl, jsfc, js
    REAL(wp)     :: zhuv, zhtq, zephum, zc2es, zc3les, zc3ies, zc4les, zc4ies
    REAL(wp)     :: zrat, zcbn, zcbs, zcbu, zmerge, zred
    REAL(wp)     :: zh2m(kbdim), zqs1(kbdim), zrh2m(kbdim), zcvm3(kbdim), zcvm4(kbdim)
    REAL(wp)     :: zaph2m(kbdim), zqs2(kbdim), zq2m(kbdim), zfrac(kbdim)
    REAL(wp)     :: ua(kbdim)
    REAL(wp), POINTER :: pbtile(:,:)
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    CONTIGUOUS :: pbtile
#endif
    !CONSTANTS
    zhuv          = 10._wp * grav
    zhtq          = 2._wp * grav
    zephum        = 5.e-2_wp

    zc2es         = 610.78_wp * rdv
    zc3les        = 17.269_wp
    zc3ies        = 21.875_wp
    zc4les        = 35.86_wp
    zc4ies        = 7.66_wp 
   
    ! set total- and tile-fields to zero in order to avoid uninitialised values

    psfcWind_tile(1:kbdim,:) = cdimissval
    puas_tile    (1:kbdim,:) = cdimissval
    pvas_tile    (1:kbdim,:) = cdimissval
    ptas_tile    (1:kbdim,:) = cdimissval
    pdew2_tile   (1:kbdim,:) = cdimissval

    DO jsfc = 1,ksfc_type

    ! check for masks
    !
    is(jsfc) = 0
     DO jl = 1,kproma
       IF(pfrc(jl,jsfc).GT.0.0_wp) THEN
         is(jsfc) = is(jsfc) + 1
         loidx(is(jsfc),jsfc) = jl
       ENDIF
     ENDDO
    ENDDO

    !     Compute new t2m
    !
    DO jsfc = 1,ksfc_type
      IF ( jsfc == idx_lnd ) THEN
        ! land only
        pbtile => pbhn_tile
      ELSE
        ! water and ice
        pbtile => pbn_tile
      END IF
      DO jls=1,is(jsfc)
       jl = loidx(jls,jsfc)
       zrat   = zhtq / pgeom1(jl)
       zcbn   = LOG(1._wp + (EXP(pbtile(jl,jsfc)) - 1._wp) * zrat )
       zcbs   = -(pbtile(jl,jsfc) - pbh_tile(jl,jsfc)) * zrat
       zcbu   = -LOG(1._wp + (EXP(pbtile(jl,jsfc) - pbh_tile(jl,jsfc)) - 1._wp) * zrat)
       zmerge = MERGE(zcbs,zcbu,pri_tile(jl,jsfc) .GT. 0._wp)
       zred   = (zcbn + zmerge) / pbh_tile(jl,jsfc)
       zh2m(jl)   = pcpt_tile(jl,jsfc) + zred * (pcptgz(jl) - pcpt_tile(jl,jsfc))
       ptas_tile(jl,jsfc) = (zh2m(jl) - zhtq ) / (cpd * (1._wp + vtmpc2 * pqm1(jl)))
     ENDDO
    ENDDO
    !
    !           5.96   2M DEW POINT
    !

    DO jsfc = 1,ksfc_type

    CALL lookup_ua_list_spline('nsurf_diag(1)', kproma, is(jsfc), loidx(1,jsfc), ptm1(1), ua(1))

     DO jls=1,is(jsfc)
      jl = loidx(jls,jsfc)
      zqs1(jl)      = ua(jls) / papm1(jl)
      zqs1(jl)      = zqs1(jl) / (1._wp- vtmpc1 * zqs1(jl))
      zrh2m(jl)     = MAX(zephum, pqm1(jl) / zqs1(jl))
 
      zaph2m(jl)     = paphm1(jl) * &  ! = paphm1(1:kproma, klevp1)
          (1._wp - zhtq / ( rd * ptas_tile(jl,jsfc) * (1._wp + vtmpc1 * pqm1(jl) - pxm1(jl))))
     ENDDO

     WHERE(ptas_tile(1:kproma,jsfc) .GT. tmelt)
        zcvm3(1:kproma)   = zc3les
        zcvm4(1:kproma)   = zc4les
     ELSEWHERE
        zcvm3(1:kproma)   = zc3ies
        zcvm4(1:kproma)   = zc4ies
     ENDWHERE

    CALL lookup_ua_list_spline('nsurf_diag(2)', kbdim, is(jsfc), loidx(1,jsfc), ptas_tile(1,jsfc), ua(1))

    DO jls=1,is(jsfc)
      jl = loidx(jls,jsfc)
      zqs2(jl)      = ua(jls) / zaph2m(jl)
      zqs2(jl)      = zqs2(jl) / (1._wp- vtmpc1 * zqs2(jl))
      zq2m(jl)      = zrh2m(jl) * zqs2(jl)
      zfrac(jl)     = LOG(zaph2m(jl) * zq2m(jl) / (zc2es * (1._wp + vtmpc1 * zq2m(jl)))) / zcvm3(jl)
      pdew2_tile(jl,jsfc) = MIN(ptas_tile(jl,jsfc), (tmelt - zfrac(jl) * zcvm4(jl)) / (1._wp - zfrac(jl)))
    ENDDO
   ENDDO

    !
    !*          5.97   10M WIND COMPONENTS
    !
    DO jsfc = 1,ksfc_type
       DO jls=1,is(jsfc)
         jl = loidx(jls,jsfc)
           zrat   = zhuv / pgeom1(jl)
           zcbn   = LOG(1._wp + (EXP (pbn_tile(jl,jsfc)) - 1._wp) * zrat )
           zcbs   = -(pbn_tile(jl,jsfc) - pbm_tile(jl,jsfc)) * zrat
           zcbu   = -LOG(1._wp + (EXP (pbn_tile(jl,jsfc) - pbm_tile(jl,jsfc)) - 1._wp) * zrat)
           zmerge = MERGE(zcbs,zcbu,pri_tile(jl,jsfc) .GT. 0._wp)
           zred   = (zcbn + zmerge) / pbm_tile(jl,jsfc)
           puas_tile(jl,jsfc)    = zred * pum1(jl)
           pvas_tile(jl,jsfc)    = zred * pvm1(jl)
           psfcWind_tile(jl,jsfc)   = zred*SQRT((pum1(jl)-pocu(jl))**2+(pvm1(jl)-pocv(jl))**2)
    ! for ice and land this is identical to
    !      psfcWind_tile(jl,jsfc)   = SQRT(puas_tile(jl,jsfc)**2+pvas_tile(jl,jsfc)**2)
        ENDDO
    ENDDO

    ! Aggregate all diagnostics 
    !
    psfcWind_gbm (1:kbdim)   = 0._wp
    puas_gbm     (1:kbdim)   = 0._wp
    pvas_gbm     (1:kbdim)   = 0._wp
    ptas_gbm     (1:kbdim)   = 0._wp
    pdew2_gbm    (1:kbdim)   = 0._wp
    !
    DO jsfc = 1,ksfc_type
      DO jls = 1,is(jsfc)
! set index
      js=loidx(jls,jsfc)
        psfcWind_gbm(js) = psfcWind_gbm(js) + pfrc(js,jsfc)*psfcWind_tile(js,jsfc)
        puas_gbm    (js) = puas_gbm (js)    + pfrc(js,jsfc)*puas_tile (js,jsfc)
        pvas_gbm    (js) = pvas_gbm (js)    + pfrc(js,jsfc)*pvas_tile (js,jsfc)
        ptas_gbm    (js) = ptas_gbm  (js)   + pfrc(js,jsfc)*ptas_tile  (js,jsfc)
        pdew2_gbm   (js) = pdew2_gbm (js)   + pfrc(js,jsfc)*pdew2_tile (js,jsfc)
      ENDDO
    ENDDO
! 
! find max and min values for 2m temperature
!
    DO jl=1,kproma
        ptasmax  (jl) = MAX(ptasmax(jl),ptas_gbm(jl))
        ptasmin  (jl) = MIN(ptasmin(jl),ptas_gbm(jl))
    ENDDO
    
  END SUBROUTINE nsurf_diag
  !-------------

END MODULE mo_surface_diag
