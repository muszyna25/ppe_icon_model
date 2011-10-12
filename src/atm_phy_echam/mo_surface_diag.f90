!>
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (2011-08)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_surface_diag

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
#ifdef __ICON__
  USE mo_physical_constants,ONLY: grav, als, alv, vtmpc2, cpd
  USE mo_echam_vdiff_params,ONLY: tpfac2
#else
  USE mo_constants,         ONLY: grav, als, alv, vtmpc2, cpd
  USE mo_physc2,            ONLY: tpfac2
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wind_stress, surface_fluxes

CONTAINS
  !>
  !!
  !!
  SUBROUTINE surface_fluxes( lsfc_heat_flux, pdtime, psteplen,     &! in
                           & kproma, kbdim, klev, ksfc_type,       &! in
                           & idx_wtr, idx_ice, idx_lnd, ih, iqv,   &! in
                           & pfrc, pcfh_tile, pfac_sfc,            &! in
                           & pcptv_tile, ptsfc_tile, pqsat_tile,   &! in
                           & pca, pcs, bb,                         &! in
                           & plhflx_gbm_ac, pshflx_gbm_ac,         &! inout
                           & pevap_gbm_ac,                         &! inout
                           & plhflx_tile, pshflx_tile,             &! out
                           & pevap_tile, pevap_gbm                 )! out

    LOGICAL, INTENT(IN) :: lsfc_heat_flux
    REAL(wp),INTENT(IN) :: pdtime, psteplen
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

    REAL(wp),INTENT(INOUT) :: bb(kbdim,klev,ih:iqv)
   !REAL(wp),INTENT(INOUT) :: bb_btm (kbdim,ksfc_type,ih:iqv)

    REAL(wp),INTENT(INOUT) :: plhflx_gbm_ac(kbdim)
    REAL(wp),INTENT(INOUT) :: pshflx_gbm_ac(kbdim)
    REAL(wp),INTENT(INOUT) ::  pevap_gbm_ac(kbdim)

    REAL(wp),INTENT(OUT)   :: plhflx_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pshflx_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   ::  pevap_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   ::  pevap_gbm (kbdim)

    INTEGER  :: jsfc
    REAL(wp) :: zconst, zdqv(kbdim), zdcptv(kbdim)

    !===================================================================
    ! If surface heat fluxes (incl. latent) are switched off, set
    ! corresponding values to zero and return to the calling subroutine.
    !===================================================================
    pevap_tile(:,:) = 0._wp
    IF (.NOT.lsfc_heat_flux) THEN

       pevap_gbm (1:kproma)   = 0._wp
      plhflx_tile(1:kproma,:) = 0._wp
      pshflx_tile(1:kproma,:) = 0._wp

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

      pevap_tile(1:kproma,jsfc) =  zconst*pfac_sfc(1:kproma) &
                                & *pcfh_tile(1:kproma,jsfc)  &
                                & *zdqv(1:kproma)
    ENDDO

    ! Compute grid box mean and time integral
    ! The instantaneous grid box mean moisture flux will be passed on
    ! to the cumulus convection scheme.

    pevap_gbm(1:kproma) = 0._wp   ! "pqhfla" in echam

    DO jsfc = 1,ksfc_type
      pevap_gbm   (1:kproma) =   pevap_gbm(1:kproma)                         &
                             & + pfrc(1:kproma,jsfc)*pevap_tile(1:kproma,jsfc)

      pevap_gbm_ac(1:kproma) =   pevap_gbm_ac(1:kproma)                      &
                             & + pevap_gbm(1:kproma)*pdtime 
    ENDDO

    !-------------------------------------------------------------------
    ! Latent heat flux
    !-------------------------------------------------------------------
    ! Instantaneous values

    IF (idx_lnd<=ksfc_type) THEN
      CALL finish('','Computation of latent heat flux over land not implemented')
    END IF

    IF (idx_ice<=ksfc_type) &
    plhflx_tile(1:kproma,idx_ice) = als*pevap_tile(1:kproma,idx_ice)

    IF (idx_wtr<=ksfc_type) &
    plhflx_tile(1:kproma,idx_wtr) = alv*pevap_tile(1:kproma,idx_wtr)

    ! Accumulated grid box mean

    DO jsfc = 1,ksfc_type
      plhflx_gbm_ac(1:kproma) =   plhflx_gbm_ac(1:kproma)                        &
                              & + pfrc(1:kproma,jsfc)*plhflx_tile(1:kproma,jsfc) &
                              &  *pdtime
    ENDDO

    !-------------------------------------------------------------------
    ! Sensible heat flux
    !-------------------------------------------------------------------
    ! Instantaneous flux on each tile

    DO jsfc = 1,ksfc_type
      IF (idx_lnd<=ksfc_type) THEN
        CALL finish('','Computation of sensible heat flux over land not implemented')
        ! CYCLE
      END IF

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
    ENDDO

    ! Accumulated grid box mean

    DO jsfc = 1,ksfc_type
      pshflx_gbm_ac(1:kproma) =   pshflx_gbm_ac(1:kproma)                        &
                              & + pfrc(1:kproma,jsfc)*pshflx_tile(1:kproma,jsfc) &
                              &  *pdtime
    ENDDO

  END SUBROUTINE surface_fluxes
  !-------------
  !!
  !! Compute wind stress over each surface type
  !!
  SUBROUTINE wind_stress( lsfc_mom_flux, pdtime, psteplen,      &! in
                        & kproma, kbdim, ksfc_type,             &! in
                        & pfrc, pcfm_tile, pfac_sfc,            &! in
                        & pu_rtpfac1, pv_rtpfac1,               &! in
                        & pu_stress_gbm_ac, pv_stress_gbm_ac,   &! inout
                        & pu_stress_tile,   pv_stress_tile      )! out

    LOGICAL, INTENT(IN)    :: lsfc_mom_flux
    REAL(wp),INTENT(IN)    :: pdtime, psteplen
    INTEGER, INTENT(IN)    :: kproma, kbdim, ksfc_type

    REAL(wp),INTENT(IN)    :: pfrc            (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pcfm_tile       (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pfac_sfc        (kbdim)
    REAL(wp),INTENT(IN)    :: pu_rtpfac1      (kbdim)
    REAL(wp),INTENT(IN)    :: pv_rtpfac1      (kbdim)
    REAL(wp),INTENT(INOUT) :: pu_stress_gbm_ac(kbdim)
    REAL(wp),INTENT(INOUT) :: pv_stress_gbm_ac(kbdim)
    REAL(wp),INTENT(OUT)   :: pu_stress_tile  (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pv_stress_tile  (kbdim,ksfc_type)

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

      DO jsfc = 1,ksfc_type

        pu_stress_tile(1:kproma,jsfc) = zconst*pfac_sfc(1:kproma) &
                                      & *pcfm_tile(1:kproma,jsfc) &
                                      & *pu_rtpfac1(1:kproma)

        pv_stress_tile(1:kproma,jsfc) = zconst*pfac_sfc(1:kproma) &
                                      & *pcfm_tile(1:kproma,jsfc) &
                                      & *pv_rtpfac1(1:kproma)

        pu_stress_gbm_ac(1:kproma) = pu_stress_gbm_ac(1:kproma)      &
                                   & + pu_stress_tile(1:kproma,jsfc) &
                                   &  *pfrc(1:kproma,jsfc)*pdtime

        pv_stress_gbm_ac(1:kproma) = pv_stress_gbm_ac(1:kproma)      &
                                   & + pv_stress_tile(1:kproma,jsfc) &
                                   &  *pfrc(1:kproma,jsfc)*pdtime
      END DO

    END IF ! lsfc_mom_flux

  END SUBROUTINE wind_stress
  !-------------

END MODULE mo_surface_diag
