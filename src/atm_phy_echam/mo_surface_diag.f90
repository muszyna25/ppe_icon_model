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
  USE mo_exception,         ONLY: finish
  USE mo_physical_constants,ONLY: grav, als, alv, vtmpc2, cpd
  USE mo_echam_vdiff_params,ONLY: tpfac2

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wind_stress, surface_fluxes

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
  !!
  SUBROUTINE surface_fluxes( lsfc_heat_flux, psteplen,             &! in
                           & kproma, kbdim, klev, ksfc_type,       &! in
                           & idx_wtr, idx_ice, idx_lnd, ih, iqv,   &! in
                           & pfrc, pcfh_tile, pfac_sfc,            &! in
                           & pcptv_tile, ptsfc_tile, pqsat_tile,   &! in
                           & pca, pcs, bb,                         &! in
                           & plhflx_gbm, pshflx_gbm,               &! out
                           & pevap_gbm,                            &! out
                           & plhflx_tile, pshflx_tile,             &! out
                           & dshflx_dT_tile,                       &! out
                           & pevap_tile,                           &! out
                           & evapotranspiration)

    LOGICAL, INTENT(IN) :: lsfc_heat_flux
    REAL(wp),INTENT(IN) :: psteplen
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, ksfc_type
    INTEGER, INTENT(IN) :: idx_wtr, idx_ice, idx_lnd
    INTEGER, INTENT(IN) :: ih, iqv

    REAL(wp),INTENT(IN) :: pfrc      (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pcfh_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pfac_sfc  (kbdim)
    REAL(wp),INTENT(IN) :: pcptv_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: ptsfc_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pqsat_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pca       (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pcs       (kbdim,ksfc_type)

    REAL(wp),INTENT(IN) :: bb(kbdim,klev,ih:iqv)

    REAL(wp),INTENT(INOUT) :: plhflx_gbm(kbdim)  ! OUT
    REAL(wp),INTENT(INOUT) :: pshflx_gbm(kbdim)  ! OUT
    REAL(wp),INTENT(INOUT) :: pevap_gbm(kbdim)   ! OUT

    REAL(wp),INTENT(INOUT) :: plhflx_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(INOUT) :: pshflx_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(INOUT) :: pevap_tile(kbdim,ksfc_type)    ! OUT

    REAL(wp),INTENT(INOUT) :: dshflx_dT_tile(kbdim,ksfc_type)  ! OUT

    REAL(wp),OPTIONAL,INTENT(IN)    :: evapotranspiration(kbdim) ! present for JSBACH land

    INTEGER  :: jsfc
    REAL(wp) :: zconst, zdqv(kbdim), zdcptv(kbdim)

    !===================================================================
    ! If surface heat fluxes (incl. latent) are switched off, set
    ! corresponding values to zero and return to the calling subroutine.
    !===================================================================

    ! KF set tile-fields to zero in order to avoid uninitialised values
    ! at diagnostics and output

    pevap_tile(1:kproma,:) = 0._wp
    pevap_gbm (1:kproma)   = 0._wp
    DO jsfc = 1,ksfc_type
      IF (jsfc /= idx_lnd .OR. .NOT. PRESENT(evapotranspiration)) THEN  ! for JSBACH land, the fluxes are already in these arrays
        plhflx_tile(1:kproma,jsfc) = 0._wp
        pshflx_tile(1:kproma,jsfc) = 0._wp
      END IF
    END DO
    ! COMMENT: should span whole array, because there might be dead ends
    dshflx_dT_tile(:,:)= 0._wp

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

    DO jsfc = 1,ksfc_type

      ! Vertical gradient of specific humidity scaled by factor (1/tpfac1).
      ! Formula: ( qv_{tavg,klev} - qs_tile )/tpfac1
      ! Here qv_{tavg,klev} = tpfac1*qv_klev(t+dt) + (1-tpfac1)*qv_klev(t)
      !                     = tpfac1*bb_qv
      ! where bb_qv is the solution of the linear system at the lowest
      ! model level (i.e., the full level right above surface).
      ! (Formula translated from ECHAM. Question: why using the blended
      ! humidity, not the value on individual surface?)

      zdqv(1:kproma) =   bb(1:kproma,klev,iqv)*pca(1:kproma,jsfc)          &
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

    pevap_gbm(1:kproma) = 0._wp   ! "pqhfla" in echam

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

    plhflx_gbm(1:kproma) = 0.0_wp

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
        ! dshflx_dT_tile(1:kproma,jsfc) = 0._wp
      ELSE

      ! Vertical gradient of dry static energy.
      ! (Formula translated from ECHAM. Question: why using the blended
      ! dry static energy, not the value on individual surface?)

        zdcptv(1:kproma) = bb(1:kproma,klev,ih) - tpfac2*pcptv_tile(1:kproma,jsfc)

      ! Flux of dry static energy

        pshflx_tile(1:kproma,jsfc) =  zconst*pfac_sfc(1:kproma) &
                                   & *pcfh_tile(1:kproma,jsfc)  &
                                   & *zdcptv(1:kproma)

      ! Subtract contribution from latent heat
      !  CpTv = CpT(1+vtmpc2*qv)
      !  => CpT = CpTv - CpT*vtmpc2*qv
      ! Question: second term is nonlinear?!

        pshflx_tile(1:kproma,jsfc) =  pshflx_tile(1:kproma,jsfc)       &
                                   & - ptsfc_tile(1:kproma,jsfc)*cpd   &
                                   &  *pevap_tile(1:kproma,jsfc)*vtmpc2

      ! KF: For the Sea-ice model!
      ! attempt to made a first guess of temperature tendency for SHF
      ! over ICE! by assuming only the cp*delta(T) matters. So: d(SHF)/deltaT

        dshflx_dT_tile(1:kproma,jsfc) =  pshflx_tile(1:kproma,jsfc) &
          &                           / zdcptv(1:kproma) !*(1._wp+pevap_tile(1:kproma,jsfc)*vtmpc2)
      END IF

    ENDDO


    ! grid box mean

    pshflx_gbm(1:kproma) = 0.0_wp

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
    REAL(wp),INTENT(INOUT) :: pu_stress_gbm   (kbdim)           ! OUT
    REAL(wp),INTENT(INOUT) :: pv_stress_gbm   (kbdim)           ! OUT
    REAL(wp),INTENT(INOUT) :: pu_stress_tile  (kbdim,ksfc_type) ! OUT
    REAL(wp),INTENT(INOUT) :: pv_stress_tile  (kbdim,ksfc_type) ! OUT

    INTEGER  :: jsfc
    REAL(wp) :: zconst

    !===================================================================
    ! If surface momentum fluxes is switched off (i.e., using free slip
    ! boundary condition), set wind stress to zero and return to the
    ! calling subroutine.
    !===================================================================
    IF (.NOT.lsfc_mom_flux) THEN

      pu_stress_tile(1:kproma,:) = 0._wp
      pv_stress_tile(1:kproma,:) = 0._wp

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

      pu_stress_gbm(1:kproma) = 0.0_wp
      pv_stress_gbm(1:kproma) = 0.0_wp

      DO jsfc = 1,ksfc_type

        pu_stress_tile(1:kproma,jsfc) = zconst*pfac_sfc(1:kproma) &
                                      & *pcfm_tile(1:kproma,jsfc) &
                                      & *pu_rtpfac1(1:kproma)

        pv_stress_tile(1:kproma,jsfc) = zconst*pfac_sfc(1:kproma) &
                                      & *pcfm_tile(1:kproma,jsfc) &
                                      & *pv_rtpfac1(1:kproma)

        pu_stress_gbm(1:kproma)       = pu_stress_gbm(1:kproma)         &
                                      & + pu_stress_tile(1:kproma,jsfc) &
                                      &  *pfrc(1:kproma,jsfc)

        pv_stress_gbm(1:kproma)       = pv_stress_gbm(1:kproma)         &
                                      & + pv_stress_tile(1:kproma,jsfc) &
                                      &  *pfrc(1:kproma,jsfc)
      END DO

    END IF ! lsfc_mom_flux

  END SUBROUTINE wind_stress
  !-------------

END MODULE mo_surface_diag
