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
MODULE mo_surface

  USE mo_kind,              ONLY: wp
 !USE mo_exception,         ONLY: finish
  USE mo_surface_diag,      ONLY: wind_stress, surface_fluxes
  USE mo_echam_vdiff_params,ONLY: tpfac2
  USE mo_vdiff_solver,      ONLY: ih, iqv, iu, iv, imh, imqv, imuv, &
                                & nmatrix, nvar_vdiff,              &
                                & matrix_to_richtmyer_coeff
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_surface

CONTAINS
  !>
  !!
  !!
  SUBROUTINE update_surface( lsfc_heat_flux, lsfc_mom_flux,     &! in
                           & pdtime, psteplen,                  &! in
                           & kproma, kbdim, klev, ksfc_type,    &! in
                           & idx_wtr, idx_ice, idx_lnd,         &! in
                           & pfrc, pcfh_tile, pcfm_tile,        &! in
                           & pfac_sfc, pocu, pocv,              &! in
                           & aa, aa_btm, bb, bb_btm,            &! inout
                           & pcpt_tile, pqsat_tile,             &! inout
                           & ptsfc_tile,                        &! inout
                           & pu_stress_gbm_ac, pv_stress_gbm_ac,&! inout
                           & plhflx_gbm_ac, pshflx_gbm_ac,      &! inout
                           & pevap_gbm_ac,  dshflx_dT_ac_tile,  &! inout
                           & pu_stress_tile,   pv_stress_tile,  &! out
                           & plhflx_tile, pshflx_tile,          &! out
                           & dshflx_dT_tile,                    &! out
                           & pevap_tile, pevap_gbm              )! out

    LOGICAL, INTENT(IN)    :: lsfc_heat_flux, lsfc_mom_flux
    REAL(wp),INTENT(IN)    :: pdtime, psteplen
    INTEGER, INTENT(IN)    :: kproma, kbdim, klev, ksfc_type
    INTEGER, INTENT(IN)    :: idx_wtr, idx_ice, idx_lnd

    REAL(wp),INTENT(IN)    :: pfrc      (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pcfh_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pcfm_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pfac_sfc  (kbdim)
    REAL(wp),INTENT(IN)    :: pocu      (kbdim)
    REAL(wp),INTENT(IN)    :: pocv      (kbdim)

    REAL(wp),INTENT(INOUT) :: aa     (kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT) :: aa_btm (kbdim,3,ksfc_type,imh:imqv)
    REAL(wp),INTENT(INOUT) :: bb     (kbdim,klev,nvar_vdiff)
    REAL(wp),INTENT(INOUT) :: bb_btm (kbdim,ksfc_type,ih:iqv)

    REAL(wp),INTENT(INOUT) :: pcpt_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(INOUT) :: pqsat_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(INOUT) :: ptsfc_tile (kbdim,ksfc_type)

    REAL(wp),INTENT(INOUT) :: pu_stress_gbm_ac (kbdim)
    REAL(wp),INTENT(INOUT) :: pv_stress_gbm_ac (kbdim)
    REAL(wp),INTENT(INOUT) ::    plhflx_gbm_ac (kbdim)
    REAL(wp),INTENT(INOUT) ::    pshflx_gbm_ac (kbdim)
    REAL(wp),INTENT(INOUT) ::     pevap_gbm_ac (kbdim)
    REAL(wp),INTENT(INOUT) :: dshflx_dT_ac_tile(kbdim,ksfc_type)

    REAL(wp),INTENT(OUT)   :: pu_stress_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: pv_stress_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   ::    plhflx_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   ::    pshflx_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   :: dshflx_dT_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   ::     pevap_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(OUT)   ::     pevap_gbm  (kbdim)

    INTEGER  :: jsfc, jk, jkm1, im
    REAL(wp) :: se_sum(kbdim), qv_sum(kbdim), wgt_sum(kbdim), wgt(kbdim)
    REAL(wp) :: zca(kbdim,ksfc_type), zcs(kbdim,ksfc_type)
    REAL(wp) :: zfrc_oce(kbdim)

    REAL(wp) :: zen_h (kbdim,ksfc_type)
    REAL(wp) :: zfn_h (kbdim,ksfc_type)
    REAL(wp) :: zen_qv(kbdim,ksfc_type)
    REAL(wp) :: zfn_qv(kbdim,ksfc_type)

    REAL(wp) :: cair_lnd (kbdim)
    REAL(wp) :: csat_lnd (kbdim)

    !===================================================================
    ! BEFORE CALLING land/ocean/ice model
    !===================================================================
    ! Compute wind stress at the old time step.
    ! At this point bb(:,klev,iu) = u_klev(t)/tpfac1 (= udif in echam)
    !               bb(:,klev,iv) = v_klev(t)/tpfac1 (= vdif in echam)

    CALL wind_stress( lsfc_mom_flux, pdtime, psteplen,     &! in
                    & kproma, kbdim, ksfc_type,            &! in
                    & pfrc, pcfm_tile, pfac_sfc,           &! in
                    & bb(:,klev,iu), bb(:,klev,iv),        &! in
                    & pu_stress_gbm_ac, pv_stress_gbm_ac,  &! inout
                    & pu_stress_tile,   pv_stress_tile     )! out

    ! Turbulent transport of moisture:
    ! - finish matrix set up;
    ! - perform bottom level elimination;
    ! - convert matrix entries to Richtmyer-Morton coefficients

    CALL matrix_to_richtmyer_coeff( kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
                                  & cair_lnd, csat_lnd,                      &! in
                                  & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),      &! in
                                  & aa_btm, bb_btm,                          &! inout
                                  & zen_h, zfn_h, zen_qv, zfn_qv             )! out

    ! Set the evapotranspiration coefficients, to be used later in
    ! blending and in diagnoising surface fluxes.

    zca(1:kproma,:) = 1._wp
    zcs(1:kproma,:) = 1._wp

    IF (idx_lnd<=ksfc_type) THEN
      zca(1:kproma,idx_lnd) = cair_lnd(1:kproma)
      zcs(1:kproma,idx_lnd) = csat_lnd(1:kproma)
    END IF

    !===================================================================
    ! CALL surface model(s)
    !===================================================================
    ! Calculate surface temperature over land (and sea ice?)

    ! CALL jsbach_interface
    ! CALL sea_ice_thermodynamics ? 

    !===================================================================
    ! AFTER CALLING land/ocean/ice model
    !===================================================================
    ! Turbulent transport of moisture and dry static energy:
    ! Get solution of the two variables on the lowest model level.
    !-------------------------------------------------------------------
    ! - Over individual tiles
    !   For echam developers: relationship to "update_surface" of echam6:
    !   bb_btm(:,jsfc,ih) : tpfac2*land%ztklevl, tpfac2*ice%ztklevi, tpfac2*ocean%ztklevw
    !   bb_btm(:,jsfc,iqv): tpfac2*land%zqklevl, tpfac2*ice%zqklevi, tpfac2*ocean%zqklevw

    DO jsfc = 1,ksfc_type
       bb_btm(1:kproma,jsfc,ih)  = tpfac2*(    zen_h (1:kproma,jsfc) &
                                 &         *pcpt_tile(1:kproma,jsfc) &
                                 &         +   zfn_h (1:kproma,jsfc) )

       bb_btm(1:kproma,jsfc,iqv) = tpfac2*(    zen_qv(1:kproma,jsfc) &
                                 &        *pqsat_tile(1:kproma,jsfc) &
                                 &        +    zfn_qv(1:kproma,jsfc) )
    END DO

    ! - Grid box mean
    !   For echam developers: relationship to "update_surface" of echam6:
    !   bb(:,klev,ih) : ztdif_new
    !   bb(:,klev,iqv): zqdif_new

     se_sum(1:kproma) = 0._wp    ! sum of weighted solution
     qv_sum(1:kproma) = 0._wp    ! sum of weighted solution
    wgt_sum(1:kproma) = 0._wp    ! sum of weights

    DO jsfc = 1,ksfc_type
           wgt(1:kproma) =  pfrc(1:kproma,jsfc)*pcfh_tile(1:kproma,jsfc)*pfac_sfc(1:kproma)
       wgt_sum(1:kproma) = wgt_sum(1:kproma) + wgt(1:kproma)
        se_sum(1:kproma) =  se_sum(1:kproma) + bb_btm(1:kproma,jsfc,ih)*wgt(1:kproma)
        qv_sum(1:kproma) =  qv_sum(1:kproma) + bb_btm(1:kproma,jsfc,iqv) &
                         &                       *wgt(1:kproma)*zca(1:kproma,jsfc)
    ENDDO

    IF (lsfc_heat_flux) THEN
      bb(1:kproma,klev,ih ) = se_sum(1:kproma)/wgt_sum(1:kproma)
      bb(1:kproma,klev,iqv) = qv_sum(1:kproma)/wgt_sum(1:kproma)
    ELSE
      jsfc = 1
      bb(1:kproma,klev,ih ) = bb_btm(1:kproma,jsfc,ih )
      bb(1:kproma,klev,iqv) = bb_btm(1:kproma,jsfc,iqv)
    END IF

    !-------------------------------------------------------------------
    ! Turbulent transport of u and v: adjust the right-hand side vector,
    ! then perform the bottom level elimination to get the solution
    !-------------------------------------------------------------------
    ! Add additional terms to the r.h.s. of the velocity equations
    ! to take into account ocean currents.
    ! Note that in subroutine rhs_setup the constant tpfac2 has been
    ! multiplied to the r.h.s. array bb. Thus the additional terms here
    ! need to be scaled by the same factor.

    IF (idx_wtr.LE.ksfc_type) THEN   ! Open water is considered
      IF (idx_ice.LE.ksfc_type) THEN ! Sea ice is also considered
        zfrc_oce(1:kproma) = pfrc(1:kproma,idx_wtr)+pfrc(1:kproma,idx_ice)
      ELSE ! only open water
        zfrc_oce(1:kproma) = pfrc(1:kproma,idx_wtr)
      ENDIF
      bb(1:kproma,klev,iu) =   bb(1:kproma,klev,iu)                   &
                           & - pocu(1:kproma)*zfrc_oce(1:kproma)*tpfac2
      bb(1:kproma,klev,iv) =   bb(1:kproma,klev,iv)                   &
                           & - pocv(1:kproma)*zfrc_oce(1:kproma)*tpfac2
    ENDIF

    ! Bottom level elimination

    im   = imuv
    jk   = klev    ! Bottom level index
    jkm1 = jk - 1

    aa(1:kproma,jk,2,im) =  aa(1:kproma,jk,2,im)                      &
                         & -aa(1:kproma,jk,1,im)*aa(1:kproma,jkm1,3,im)
    aa(1:kproma,jk,3,im) =  aa(1:kproma,jk,3,im)/aa(1:kproma,jk,2,im)

    bb(1:kproma,jk,iu) = (bb(1:kproma,jk,iu)                         &
                       & -aa(1:kproma,jk,1,im)*bb(1:kproma,jkm1,iu)) &
                       & /aa(1:kproma,jk,2,im)

    bb(1:kproma,jk,iv) = (bb(1:kproma,jk,iv)                         &
                       & -aa(1:kproma,jk,1,im)*bb(1:kproma,jkm1,iv)) &
                       & /aa(1:kproma,jk,2,im)

   !-------------------------------------------------------------------
   ! Various diagnostics
   !-------------------------------------------------------------------

   CALL surface_fluxes( lsfc_heat_flux, pdtime, psteplen,     &! in
                      & kproma, kbdim, klev, ksfc_type,       &! in
                      & idx_wtr, idx_ice, idx_lnd, ih, iqv,   &! in
                      & pfrc, pcfh_tile, pfac_sfc,            &! in
                      & pcpt_tile, ptsfc_tile, pqsat_tile,    &! in
                      & zca, zcs, bb(:,:,ih:iqv),             &! in
                      & plhflx_gbm_ac, pshflx_gbm_ac,         &! inout
                      & pevap_gbm_ac,  dshflx_dT_ac_tile,     &! inout
                      & plhflx_tile, pshflx_tile,             &! out
                      & dshflx_dT_tile,                       &! out
                      & pevap_tile, pevap_gbm                 )! out

  END SUBROUTINE update_surface
  !-------------

END MODULE mo_surface
