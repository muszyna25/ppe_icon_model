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
  USE mo_exception,         ONLY: finish
#ifdef __ICON__
  USE mo_echam_vdiff_params,ONLY: tpfac2
  USE mo_vdiff_solver,      ONLY: ih, iqv, iu, iv, imh, imqv, imuv, &
                                & nmatrix, nvar_vdiff,              &
                                & matrix_to_richtmyer_coeff
#else
  USE mo_physc2,            ONLY: tpfac2
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_surface

CONTAINS
  !>
  !!
  !!
  SUBROUTINE update_surface( lsfc_heat_flux,                    &! in
                           & kproma, kbdim, klev, ksfc_type,    &! in
                           & idx_wtr, idx_ice, idx_lnd,         &! in
                           & pfrc, pcfh_tile, pfac_sfc,         &! in
                           & pcpt_tile, pqsat_tile, pocu, pocv, &! in
                           & aa, aa_btm, bb, bb_btm             )! inout

    LOGICAL, INTENT(IN)    :: lsfc_heat_flux
    INTEGER, INTENT(IN)    :: kproma, kbdim, klev, ksfc_type
    INTEGER, INTENT(IN)    :: idx_wtr, idx_ice, idx_lnd

    REAL(wp),INTENT(IN)    :: pfrc      (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pcfh_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pfac_sfc  (kbdim)
    REAL(wp),INTENT(IN)    :: pcpt_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pqsat_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pocu      (kbdim)
    REAL(wp),INTENT(IN)    :: pocv      (kbdim)

    REAL(wp),INTENT(INOUT) :: aa     (kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT) :: aa_btm (kbdim,3,ksfc_type,imh:imqv)
    REAL(wp),INTENT(INOUT) :: bb     (kbdim,klev,nvar_vdiff)
    REAL(wp),INTENT(INOUT) :: bb_btm (kbdim,ksfc_type,ih:iqv)

    INTEGER  :: jsfc, jk, jkm1, im, jvar
    REAL(wp) :: zfrc_oce(kbdim)
    REAL(wp) :: zcair(kbdim), zcsat(kbdim), wgt(kbdim)
    REAL(wp) :: se_sum(kbdim), qv_sum(kbdim), wgt_sum(kbdim)

    REAL(wp) :: zen_h (kbdim,ksfc_type)
    REAL(wp) :: zfn_h (kbdim,ksfc_type)
    REAL(wp) :: zen_qv(kbdim,ksfc_type)
    REAL(wp) :: zfn_qv(kbdim,ksfc_type)


    !-------------------------------------------------------------------
    ! For moisture: 
    ! - finish matrix set up;
    ! - perform bottom level elimination; 
    ! - convert matrix entries to Richtmyer-Morton coefficients
    !-------------------------------------------------------------------
    IF (idx_lnd<=ksfc_type) CALL finish('','land surface active in update_surface')
    zcair(:) = 1._wp
    zcsat(:) = 1._wp

    CALL matrix_to_richtmyer_coeff( kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
                                  & zcair, zcsat,                          &! in, dummy
                                  & aa(:,:,:,imh:imqv), bb(:,:,ih:iqv),    &! in
                                  & aa_btm, bb_btm,                        &! inout
                                  & zen_h, zfn_h, zen_qv, zfn_qv           )! out

    !-------------------------------------------------------------------
    ! Calculate surface temperature over land and sea ice;
    ! Provide surface temperature over ocean.
    !-------------------------------------------------------------------

    ! call jsbach and sea ice model here.

    !-------------------------------------------------------------------
    ! For moisture and dry static energy: 
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
        qv_sum(1:kproma) =  qv_sum(1:kproma) + bb_btm(1:kproma,jsfc,iqv)*wgt(1:kproma)
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
    ! For u and v: adjust the right-hand side vector, then perform
    ! the bottom level elimination to get the solution
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

   !!-------------------------------------------------------------------
   !! Various diagnostics
   !!-------------------------------------------------------------------
   !CALL wind_stress( ... )
   !CALL lh_flux( ... )
   !CALL sh_flux( ... )

  END SUBROUTINE update_surface
  !-------------

END MODULE mo_surface
