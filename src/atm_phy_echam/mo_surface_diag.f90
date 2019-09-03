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
  SUBROUTINE surface_fluxes( jcs, kproma, kbdim, ksfc_type,        &! in
                           & idx_wtr, idx_ice, idx_lnd, ih, iqv,   &! in
                           & psteplen,                             &! in
                           & pfrc, lsmask, alake,                  &! in
                           & pcfh_tile, pfac_sfc,                  &! in
                           & pcptv_tile, pqsat_tile,               &! in
                           & pca, pcs, bb_btm,                     &! in
                           & plhflx_lnd, plhflx_lwtr, plhflx_lice, &! in for JSBACH land and lakes
                           & pshflx_lnd, pshflx_lwtr, pshflx_lice, &! in for JSBACH land and lakes
                           & pevap_lnd, pevap_lwtr, pevap_lice,    &! in for JSBACH land and lakes
                           & plhflx_gbm, pshflx_gbm,               &! out
                           & pevap_gbm,                            &! out
                           & plhflx_tile, pshflx_tile,             &! out
                           & pevap_tile )                           ! out

    REAL(wp),INTENT(IN) :: psteplen
    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, ksfc_type
    INTEGER, INTENT(IN) :: idx_wtr, idx_ice, idx_lnd
    INTEGER, INTENT(IN) :: ih, iqv

    REAL(wp),INTENT(IN) :: pfrc(:,:)       ! (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: lsmask(:)       ! (kbdim)
    REAL(wp),INTENT(IN) :: alake(:)        ! (kbdim)
    REAL(wp),INTENT(IN) :: pcfh_tile(:,:)  ! (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pfac_sfc(:)     ! (kbdim)
    REAL(wp),INTENT(IN) :: pcptv_tile(:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pqsat_tile(:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pca(:,:)        ! (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pcs(:,:)        ! (kbdim,ksfc_type)

    REAL(wp),INTENT(IN) :: bb_btm(:,:,ih:) ! (kbdim,ksfc_type,ih:iqv)

    REAL(wp),INTENT(OUT) :: plhflx_gbm(:) ! (kbdim)
    REAL(wp),INTENT(OUT) :: pshflx_gbm(:) ! (kbdim)
    REAL(wp),INTENT(OUT) :: pevap_gbm(:)  ! (kbdim)

    REAL(wp),INTENT(OUT)   :: plhflx_tile(:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pshflx_tile(:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pevap_tile(:,:)  ! (kbdim,ksfc_type)

    ! Input for JSBACH land and lakes
    REAL(wp),INTENT(IN), DIMENSION(:) :: & ! DIMENSION(kbdim) 
      &                         plhflx_lnd, plhflx_lwtr, plhflx_lice, &
      &                         pshflx_lnd, pshflx_lwtr, pshflx_lice, &
      &                         pevap_lnd, pevap_lwtr, pevap_lice

    INTEGER  :: jsfc, jk, jl
    REAL(wp) :: zconst, zdqv, zdcptv

    !$ACC DATA PRESENT( pfrc, lsmask, alake, pcfh_tile, pfac_sfc, pcptv_tile,  &
    !$ACC               pqsat_tile, pca, pcs, bb_btm, plhflx_gbm, plhflx_lice, &
    !$ACC               plhflx_tile, pshflx_gbm, pshflx_lnd, pshflx_lwtr,      &
    !$ACC               pshflx_lice, pevap_gbm, pevap_lnd, pevap_lwtr,         &
    !$ACC               pevap_lice, pevap_tile, plhflx_lnd, plhflx_lwtr,       &
    !$ACC               pshflx_tile )

    !===================================================================
    ! If surface heat fluxes (incl. latent) are switched off, set
    ! corresponding values to zero and return to the calling subroutine.
    !===================================================================

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
      !$ACC LOOP GANG VECTOR
      DO jk = 1, kbdim
        plhflx_tile(jk,jsfc) = 0._wp
        pshflx_tile(jk,jsfc) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL
    ! DO jsfc = 1,ksfc_type
    !   IF (jsfc /= idx_lnd) THEN  ! for JSBACH land, the fluxes are already in these arrays
    !     plhflx_tile(jcs:kproma,jsfc) = 0._wp
    !     pshflx_tile(jcs:kproma,jsfc) = 0._wp
    !   END IF
    ! END DO

    !===================================================================
    ! Otherwise compute diagnostics
    !===================================================================
    zconst = 1._wp/psteplen

    !-------------------------------------------------------------------
    ! Moisture fluxes (aka evaporation rates)
    !-------------------------------------------------------------------
    ! Instantaneous moisture flux on each tile

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
      !$ACC LOOP GANG VECTOR
      DO jk = 1, kbdim
        pevap_tile(jk,jsfc) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
      !$ACC LOOP GANG VECTOR PRIVATE( zdqv )
      DO jl = jcs, kproma
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

        zdqv =   bb_btm(jl,jsfc,iqv)*pca(jl,jsfc)          &
                       & - tpfac2*pqsat_tile(jl,jsfc)*pcs(jl,jsfc)

        ! Moisture flux ( = evaporation). Formula:
        ! (g*psteplen)**(-1)*[  tpfac1*g*psteplen*(air density)*(exchange coef)
        !                     *(tpfac1)**(-1)*( qv_{tavg,klev} - qs_tile ) ]

        ! On solid land (i.e. without lakes)
        !
        IF (jsfc == idx_lnd) THEN
          pevap_tile(jl,jsfc) = pevap_lnd(jl)
        END IF

        ! On open water, ocean and lakes
        !
        IF (jsfc == idx_wtr) THEN
  !!$        WHERE (alake(jl) > 0._wp)
  !!$          pevap_tile(jl,idx_wtr) = pevap_lwtr(jl)
  !!$        ELSE WHERE
  !!$          pevap_tile(jl,idx_wtr) =  zconst*pfac_sfc(jl)   &
  !!$                                       & *pcfh_tile(jl,idx_wtr) &
  !!$                                       & *zdqv
  !!$        END WHERE
          IF (lsmask(jl) < 1._wp) THEN
              pevap_tile(jl,jsfc) = alake(jl)*pevap_lwtr(jl)               & ! lakes
                   &                     +(1._wp-lsmask(jl)-alake(jl))           & ! ocean
                   &                     *zconst*pfac_sfc(jl)*pcfh_tile(jl,jsfc) &
                   &                     *zdqv
              pevap_tile(jl,jsfc) = pevap_tile(jl,jsfc)/(1._wp-lsmask(jl))
          ELSE
              pevap_tile(jl,jsfc) = 0.0_wp
          END IF
        END IF

        ! On ice covered water, ocean and lakes
        !
        IF (jsfc == idx_ice) THEN
  !!$        WHERE (alake(jl) > 0._wp)
  !!$          pevap_tile(jl,idx_ice) = pevap_lice(jl)
  !!$        ELSE WHERE
  !!$          pevap_tile(jl,idx_ice) =  zconst*pfac_sfc(jl)   &
  !!$                                       & *pcfh_tile(jl,idx_ice) &
  !!$                                       & *zdqv
  !!$        END WHERE
           IF (lsmask(jl) < 1._wp) THEN
              pevap_tile(jl,jsfc) = alake(jl)*pevap_lice(jl)               & ! lakes
                   &                     +(1._wp-lsmask(jl)-alake(jl))           & ! ocean
                   &                     *zconst*pfac_sfc(jl)*pcfh_tile(jl,jsfc) &
                   &                     *zdqv
              pevap_tile(jl,jsfc) = pevap_tile(jl,jsfc)/(1._wp-lsmask(jl))
           ELSE
              pevap_tile(jl,jsfc) = 0.0_wp
           END IF
        END IF

      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! Compute grid box mean and time integral
    ! The instantaneous grid box mean moisture flux will be passed on
    ! to the cumulus convection scheme.

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jk = 1, kbdim
      pevap_gbm(jk) = 0._wp   ! "pqhfla" in echam
    END DO
    !$ACC END PARALLEL

    DO jsfc = 1,ksfc_type
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs, kproma
        pevap_gbm(jl) = pevap_gbm(jl) + pfrc(jl,jsfc)*pevap_tile(jl,jsfc)
      END DO
      !$ACC END PARALLEL
    ENDDO

    !-------------------------------------------------------------------
    ! Latent heat flux
    !-------------------------------------------------------------------
    ! Instantaneous values

    IF (idx_lnd <= ksfc_type) THEN
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs, kproma
        plhflx_tile(jl,idx_lnd) = plhflx_lnd(jl)
      END DO
      !$ACC END PARALLEL
    END IF
    IF (idx_wtr <= ksfc_type) THEN
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs, kproma
        IF (alake(jl) > 0._wp) THEN
          plhflx_tile(jl,idx_wtr) = plhflx_lwtr(jl)
        ELSE
          plhflx_tile(jl,idx_wtr) = alv*pevap_tile(jl,idx_wtr)
        END IF
      END DO
      !$ACC END PARALLEL
    END IF
    IF (idx_ice <= ksfc_type) THEN
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs, kproma
        IF (alake(jl) > 0._wp) THEN
          plhflx_tile(jl,idx_ice) = plhflx_lice(jl)
        ELSE
          plhflx_tile(jl,idx_ice) = als*pevap_tile(jl,idx_ice)
        END IF
      END DO
      !$ACC END PARALLEL
    END IF

    ! Accumulated grid box mean

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jk = 1, kbdim
       plhflx_gbm(jk) = 0.0_wp
    END DO
    !$ACC END PARALLEL

    DO jsfc = 1,ksfc_type
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs, kproma
        plhflx_gbm(jl) = plhflx_gbm(jl) + pfrc(jl,jsfc)*plhflx_tile(jl,jsfc)
      END DO
      !$ACC END PARALLEL
    ENDDO

    !-------------------------------------------------------------------
    ! Sensible heat flux
    !-------------------------------------------------------------------
    ! Instantaneous flux on each tile

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
      !$ACC LOOP GANG VECTOR PRIVATE( zdcptv )
      DO jl = jcs, kproma

        ! Vertical gradient of dry static energy.
        ! (Formula translated from ECHAM. Question: why using the blended
        ! dry static energy, not the value on individual surface?)
        ! bb was replaced by bb_btm (according to E. Roeckner), now not blended
        ! quantity used.

        zdcptv = bb_btm(jl,jsfc,ih) - tpfac2*pcptv_tile(jl,jsfc)

        ! Flux of dry static energy

        IF (jsfc == idx_lnd) THEN
          pshflx_tile(jl,jsfc) = pshflx_lnd(jl)
        END IF
        IF (jsfc == idx_wtr) THEN
          IF (alake(jl) > 0._wp) THEN
            pshflx_tile(jl,jsfc) = pshflx_lwtr(jl)
          ELSE
            pshflx_tile(jl,jsfc) =  zconst*pfac_sfc(jl)*pcfh_tile(jl,jsfc)*zdcptv
          END IF
        END IF
        IF (jsfc == idx_ice) THEN
          IF (alake(jl) > 0._wp) THEN
            pshflx_tile(jl,jsfc) = pshflx_lice(jl)
          ELSE
           pshflx_tile(jl,jsfc) =  zconst*pfac_sfc(jl)*pcfh_tile(jl,jsfc)*zdcptv
          END IF
        END IF

      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! grid box mean

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jk = 1, kbdim
      pshflx_gbm(jk) = 0.0_wp
    END DO
    !$ACC END PARALLEL

    DO jsfc = 1,ksfc_type
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs, kproma
        pshflx_gbm(jl) = pshflx_gbm(jl) + pfrc(jl,jsfc)*pshflx_tile(jl,jsfc)
      END DO
      !$ACC END PARALLEL
    ENDDO

    !$ACC END DATA

  END SUBROUTINE surface_fluxes
  !-------------
  !!
  !! Compute wind stress over each surface type
  !!
  SUBROUTINE wind_stress( jcs, kproma, kbdim, ksfc_type,        &! in
                        & psteplen,                             &! in
                        & pfrc, pcfm_tile, pfac_sfc,            &! in
                        & pu_rtpfac1, pv_rtpfac1,               &! in
                        & pu_stress_gbm,  pv_stress_gbm,        &! out
                        & pu_stress_tile, pv_stress_tile        )! out

    REAL(wp),INTENT(IN)    :: psteplen
    INTEGER, INTENT(IN)    :: jcs, kproma, kbdim, ksfc_type

    REAL(wp),INTENT(IN)    :: pfrc            (:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pcfm_tile       (:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pfac_sfc        (:)   ! (kbdim)
    REAL(wp),INTENT(IN)    :: pu_rtpfac1      (:)   ! (kbdim)
    REAL(wp),INTENT(IN)    :: pv_rtpfac1      (:)   ! (kbdim)
    REAL(wp),INTENT(OUT)   :: pu_stress_gbm   (:)   ! (kbdim)
    REAL(wp),INTENT(OUT)   :: pv_stress_gbm   (:)   ! (kbdim)
    REAL(wp),INTENT(OUT)   :: pu_stress_tile  (:,:) ! (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pv_stress_tile  (:,:) ! (kbdim,ksfc_type)

    INTEGER  :: jsfc
    REAL(wp) :: zconst
    ! Local variables

    INTEGER  :: loidx  (kbdim,ksfc_type) !< counter for masks
    INTEGER  :: is     (ksfc_type)       !< counter for masks
    INTEGER  :: jls, jl, js

    zconst = 1._wp/psteplen

    ! Compute wind stress over each surface type, then accumulate
    ! grid box mean. Formula for wind stress:
    !   (grav*psteplen)**(-1)
    !  *[grav*psteplen*tpfac1*(air density)]
    !  *(surface turbulent exchange coeff)
    !  *[(u-/v-wind at lowest model level)/tpfac1]

    !$ACC DATA PRESENT( pu_stress_tile, pv_stress_tile, pfac_sfc, pfrc,        &
    !$ACC               pu_stress_gbm, pv_stress_gbm, pcfm_tile, pu_rtpfac1,   &
    !$ACC               pv_rtpfac1 )  &
    !$ACC      CREATE( is, loidx )

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
      !$ACC LOOP GANG VECTOR
      DO jl = 1, kbdim
        pu_stress_tile(jl,jsfc) = 0.0_wp
        pv_stress_tile(jl,jsfc) = 0.0_wp
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = 1, kbdim
      pu_stress_gbm (jl)   = 0.0_wp
      pv_stress_gbm (jl)   = 0.0_wp
    END DO
    !$ACC END PARALLEL

    ! check for masks
    !
    ! GPU: Compute index list on CPU due to issues with ACC ATOMIC
    !$ACC UPDATE HOST( pfrc )
    DO jsfc = 1,ksfc_type
       is(jsfc) = 0
       DO jl = jcs,kproma
          IF(pfrc(jl,jsfc).GT.0.0_wp) THEN
             is(jsfc) = is(jsfc) + 1
             loidx(is(jsfc),jsfc) = jl
          ENDIF
       ENDDO
    ENDDO
    !$ACC UPDATE DEVICE( is, loidx )

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
       !$ACC LOOP GANG VECTOR PRIVATE( js )
       DO jls = 1,is(jsfc)
          ! set index
          js=loidx(jls,jsfc)

          pu_stress_tile(js,jsfc) = zconst*pfac_sfc(js) *pcfm_tile(js,jsfc)*pu_rtpfac1(js)
          pv_stress_tile(js,jsfc) = zconst*pfac_sfc(js) *pcfm_tile(js,jsfc)*pv_rtpfac1(js)
       END DO
    END DO
    !$ACC END PARALLEL

    DO jsfc = 1,ksfc_type
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR PRIVATE( js )
       DO jls = 1,is(jsfc)
          ! set index
          js=loidx(jls,jsfc)

          pu_stress_gbm(js)       = pu_stress_gbm(js) + pu_stress_tile(js,jsfc)*pfrc(js,jsfc)
          pv_stress_gbm(js)       = pv_stress_gbm(js) + pv_stress_tile(js,jsfc)*pfrc(js,jsfc)
       END DO
       !$ACC END PARALLEL
    END DO

    !$ACC END DATA

  END SUBROUTINE wind_stress
  !-------------
  !!
  !! Compute diagnostics: 10m wind, u and v in 10m,
  !!                      temperature in 2m, dew point temperature in 2m
  !!
  SUBROUTINE nsurf_diag( jcs, kproma, kbdim, ksfc_type,   &! in
                       & idx_lnd,                         &! in
                       & pfrc,                            &! in
                       & pqm1,                            &! in humidity
                       & ptm1,                            &
                       & papm1,     paphm1,               &
                       & pxm1,                            &
                       & pum1,      pvm1,                 &
                       & pocu,      pocv,                 &
                       & pzf,                             &! in height of lowermost full level (m)
                       & pzs,                             &! in height of surface (m)
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

    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, ksfc_type
    INTEGER, INTENT(IN) :: idx_lnd

    REAL(wp),INTENT(IN), DIMENSION(:,:) :: &                !< DIMENSION(kbdim,ksfc_type)
                                pfrc                        !< fraction of the grid box occupied by
                                                            !< each surface type
    REAL(wp), INTENT(in)     :: pqm1(:)                     !< (kbdim)
    REAL(wp), INTENT(in)     :: pzf(:), pzs(:)              !< (kbdim)
    REAL(wp), INTENT(in)     :: pcptgz(:)                   !< (kbdim) dry static energy at surface level
    REAL(wp), INTENT(in)     :: pcpt_tile(:,:)              !< (kbdim,ksfc_type) dry static energy on tiles
    REAL(wp), TARGET, INTENT(in) :: pbn_tile(:,:)           !< (kbdim,ksfc_type) for diagnostics
    REAL(wp), TARGET, INTENT(in) :: pbhn_tile(:,:)          !< (kbdim,ksfc_type) for diagnostics
    REAL(wp), INTENT(in)     :: pbh_tile(:,:)               !< (kbdim,ksfc_type) for diagnostics
    REAL(wp), INTENT(in)     :: pbm_tile(:,:)               !< (kbdim,ksfc_type) for diagnostics
    REAL(wp), INTENT(in)     :: pri_tile(:,:)               !< (kbdim,ksfc_type) moist Richardson number
    REAL(wp), INTENT(in)     :: ptm1(:), papm1(:), pxm1(:)  !< (kbdim)
    REAL(wp), INTENT(in)     :: pum1(:), pvm1(:), paphm1(:) !< (kbdim) =paphm1(kbdim, klevp1)
    REAL(wp), INTENT(in)     :: pocu(:), pocv(:)            !< (kbdim)
    REAL(wp), INTENT(out)    :: psfcWind_gbm(:)             !< (kbdim)
    REAL(wp), INTENT(out)    :: ptas_gbm(:)                 !< (kbdim)
    REAL(wp), INTENT(out)    :: pdew2_gbm(:)                !< (kbdim)
    REAL(wp), INTENT(out)    :: puas_gbm(:)                 !< (kbdim)
    REAL(wp), INTENT(out)    :: pvas_gbm(:)                 !< (kbdim)
    REAL(wp), INTENT(inout)  :: ptasmax(:), ptasmin(:)      !< (kbdim)
    REAL(wp), INTENT(out), DIMENSION(:,:) ::      &         !< DIMENSION(kbdim,ksfc_type)
                                psfcWind_tile,    &
                                ptas_tile,        &
                                pdew2_tile,       &
                                puas_tile,        &
                                pvas_tile

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

    !$ACC DATA PRESENT( pfrc, pqm1, pzf, pzs, pcptgz, pcpt_tile, pbn_tile,     &
    !$ACC               pbhn_tile, pbh_tile, pbm_tile, pri_tile, ptm1, papm1,  &
    !$ACC               pxm1, pum1, pvm1, paphm1, pocu, pocv, ptasmax,         &
    !$ACC               ptasmin )                                              &
    !$ACC      PRESENT( psfcWind_gbm, psfcWind_tile, ptas_gbm, ptas_tile,      &
    !$ACC               pdew2_gbm, pdew2_tile, puas_gbm, puas_tile, pvas_gbm,  &
    !$ACC               pvas_tile )                                            &
    !$ACC      PCREATE( zh2m, zqs1, zrh2m, zcvm3, zcvm4, zaph2m, zqs2, zq2m,   &
    !$ACC               zfrac, ua, is, loidx )

    !CONSTANTS
    zhuv          =  10._wp ! 10m
    zhtq          =   2._wp !  2m
    zephum        = 0.05_wp ! epsilon for rel. humidity

    zc2es         = 610.78_wp * rdv
    zc3les        = 17.269_wp
    zc3ies        = 21.875_wp
    zc4les        = 35.86_wp
    zc4ies        = 7.66_wp

    ! set total- and tile-fields to zero in order to avoid uninitialised values

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
      !$ACC LOOP GANG VECTOR
      DO jl = 1,kbdim
        psfcWind_tile(jl,jsfc) = cdimissval
        puas_tile    (jl,jsfc) = cdimissval
        pvas_tile    (jl,jsfc) = cdimissval
        ptas_tile    (jl,jsfc) = cdimissval
        pdew2_tile   (jl,jsfc) = cdimissval
      END DO
    END DO
    !$ACC END PARALLEL

    ! GPU: Compute index list on CPU due to issues with ACC ATOMIC
    !$ACC UPDATE HOST( pfrc )
    DO jsfc = 1,ksfc_type
      ! check for masks
      !
      is(jsfc) = jcs-1
      DO jl = jcs,kproma
        IF(pfrc(jl,jsfc).GT.0.0_wp) THEN
          is(jsfc) = is(jsfc) + 1
          loidx(is(jsfc),jsfc) = jl
        ENDIF
      ENDDO
    ENDDO
    !$ACC UPDATE DEVICE( is, loidx )

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
      !$ACC DATA PRESENT( pbtile )
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( jl, zrat, zcbn, zcbs, zcbu, zmerge, zred )
      DO jls=jcs,is(jsfc)
       jl = loidx(jls,jsfc)
       zrat   = zhtq / (pzf(jl)-pzs(jl))
       zcbn   = LOG(1._wp + (EXP(pbtile(jl,jsfc)) - 1._wp) * zrat )
       zcbs   = -(pbtile(jl,jsfc) - pbh_tile(jl,jsfc)) * zrat
       zcbu   = -LOG(1._wp + (EXP(pbtile(jl,jsfc) - pbh_tile(jl,jsfc)) - 1._wp) * zrat)
       zmerge = MERGE(zcbs,zcbu,pri_tile(jl,jsfc) .GT. 0._wp)
       zred   = (zcbn + zmerge) / pbh_tile(jl,jsfc)
       zh2m(jl)   = pcpt_tile(jl,jsfc) + zred * (pcptgz(jl) - pcpt_tile(jl,jsfc))
       ptas_tile(jl,jsfc) = (zh2m(jl) - zhtq*grav ) / (cpd * (1._wp + vtmpc2 * pqm1(jl)))
     ENDDO
     !$ACC END PARALLEL
     !$ACC END DATA
    ENDDO
    !
    !           5.96   2M DEW POINT
    !

    DO jsfc = 1,ksfc_type

      CALL lookup_ua_list_spline('nsurf_diag(1)', jcs, kproma, is(jsfc), loidx(:,jsfc), ptm1, ua)

      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( jl )
      DO jls=jcs,is(jsfc)
        jl = loidx(jls,jsfc)
        zqs1(jl)      = ua(jls) / papm1(jl)
        zqs1(jl)      = zqs1(jl) / (1._wp- vtmpc1 * zqs1(jl))
        zrh2m(jl)     = MAX(zephum, pqm1(jl) / zqs1(jl))
 
        zaph2m(jl)     = paphm1(jl) * &  ! = paphm1(jcs:kproma, klevp1)
            (1._wp - zhtq*grav / ( rd * ptas_tile(jl,jsfc) * (1._wp + vtmpc1 * pqm1(jl) - pxm1(jl))))
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs,kproma
        IF(ptas_tile(jl,jsfc) .GT. tmelt) THEN
          zcvm3(jl)   = zc3les
          zcvm4(jl)   = zc4les
        ELSE
          zcvm3(jl)   = zc3ies
          zcvm4(jl)   = zc4ies
        ENDIF
      ENDDO
      !$ACC END PARALLEL

      CALL lookup_ua_list_spline('nsurf_diag(2)', jcs, kbdim, is(jsfc), loidx(:,jsfc), ptas_tile(:,jsfc), ua)

      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( jl )
      DO jls=jcs,is(jsfc)
        jl = loidx(jls,jsfc)
        zqs2(jl)      = ua(jls) / zaph2m(jl)
        zqs2(jl)      = zqs2(jl) / (1._wp- vtmpc1 * zqs2(jl))
        zq2m(jl)      = zrh2m(jl) * zqs2(jl)
        zfrac(jl)     = LOG(zaph2m(jl) * zq2m(jl) / (zc2es * (1._wp + vtmpc1 * zq2m(jl)))) / zcvm3(jl)
        pdew2_tile(jl,jsfc) = MIN(ptas_tile(jl,jsfc), (tmelt - zfrac(jl) * zcvm4(jl)) / (1._wp - zfrac(jl)))
       ENDDO
       !$ACC END PARALLEL
    ENDDO

    !
    !*          5.97   10M WIND COMPONENTS
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
       !$ACC LOOP GANG VECTOR PRIVATE( jl, zrat, zcbn, zcbs, zcbu, zmerge, zred )
       DO jls=jcs,is(jsfc)
         jl = loidx(jls,jsfc)
           zrat   = zhuv / (pzf(jl)-pzs(jl))
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
    !$ACC END PARALLEL

    ! Aggregate all diagnostics
    !
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl = 1,kbdim
      psfcWind_gbm (jl)   = 0._wp
      puas_gbm     (jl)   = 0._wp
      pvas_gbm     (jl)   = 0._wp
      ptas_gbm     (jl)   = 0._wp
      pdew2_gbm    (jl)   = 0._wp
    END DO
    !$ACC END PARALLEL
    !
    DO jsfc = 1,ksfc_type
      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE( js )
      DO jls = jcs,is(jsfc)
! set index
      js=loidx(jls,jsfc)
        psfcWind_gbm(js) = psfcWind_gbm(js) + pfrc(js,jsfc)*psfcWind_tile(js,jsfc)
        puas_gbm    (js) = puas_gbm (js)    + pfrc(js,jsfc)*puas_tile (js,jsfc)
        pvas_gbm    (js) = pvas_gbm (js)    + pfrc(js,jsfc)*pvas_tile (js,jsfc)
        ptas_gbm    (js) = ptas_gbm  (js)   + pfrc(js,jsfc)*ptas_tile  (js,jsfc)
        pdew2_gbm   (js) = pdew2_gbm (js)   + pfrc(js,jsfc)*pdew2_tile (js,jsfc)
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!
! find max and min values for 2m temperature
!
    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jl=jcs,kproma
        ptasmax  (jl) = MAX(ptasmax(jl),ptas_gbm(jl))
        ptasmin  (jl) = MIN(ptasmin(jl),ptas_gbm(jl))
    ENDDO
    !$ACC END PARALLEL

  !$ACC END DATA

  END SUBROUTINE nsurf_diag
  !-------------

END MODULE mo_surface_diag
