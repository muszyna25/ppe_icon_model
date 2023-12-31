!>
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (2010-09)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_vdiff_solver

  USE mo_kind,              ONLY: wp
  USE mo_impl_constants,    ONLY: SUCCESS
  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_physical_constants,ONLY: rgrav, cpd, cpv
  USE mo_echam_vdiff_params,ONLY: totte_min, &
    &                             tpfac1, tpfac2, tpfac3, cchar
  USE mo_echam_vdf_config,  ONLY: echam_vdf_config
  USE mo_echam_phy_config,  ONLY: echam_phy_config
  USE mo_echam_vdf_config,  ONLY: echam_vdf_config
  USE mo_nh_testcases_nml,  ONLY: isrfc_type, shflx, lhflx

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_vdiff_solver         !< subroutine
  PUBLIC :: cleanup_vdiff_solver      !< subroutine
  PUBLIC :: matrix_setup_elim         !< subroutine
  PUBLIC :: matrix_to_richtmyer_coeff !< subroutine
  PUBLIC :: rhs_setup, rhs_elim       !< subroutines
  PUBLIC :: rhs_bksub                 !< subroutines
  PUBLIC :: vdiff_tendencies          !< subroutine
  PUBLIC :: nvar_vdiff, nmatrix       !< parameters
  PUBLIC :: ih,iqv,iu,iv              !< parameters
  PUBLIC :: imh,imqv, imuv            !< parameters
  PUBLIC :: matrix_idx,ibtm_var

  ! Module variables

  INTEGER :: nvar_vdiff    !< total number of variables affected by turbulent mixing
  INTEGER :: iu, iv, ih, iqv
  INTEGER :: ixl, ixi, ixv
  INTEGER :: itotte, ithv
  INTEGER :: itrc_start
  INTEGER :: nmatrix
  INTEGER :: imh, imqv, imuv

  INTEGER, ALLOCATABLE :: matrix_idx(:)    !< shape: (nvar_vdiff)
  INTEGER, ALLOCATABLE :: ibtm_var  (:)    !< shape: (nvar_vdiff)
  INTEGER, ALLOCATABLE :: ibtm_mtrx (:)    !< shape: (nmatrix)

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_vdiff_solver'

CONTAINS
  !>
  !!
  !! In this prototype it is assumed that the following variables are subject to
  !! turbulent mixing:
  !!
  !!    variables                              |  # of variables
  !! -----------------------------------------------------------
  !!   u, v, T, qv                             |  4
  !!   all hydrometeors                        |  khydromet
  !!   variance of cloud droplet concentration |  1
  !!   TTE                                     |  1
  !!   variance of theta_v                     |  1
  !!   additional tracers                      |  ktrac
  !! -----------------------------------------------------------
  !!
  SUBROUTINE init_vdiff_solver( khydromet, ktrac, klev )

    INTEGER,INTENT(IN) :: khydromet, ktrac
    INTEGER,INTENT(IN) :: klev
    INTEGER :: ist

    !------------------------------------------
    ! Set up index for prognostic variables
    !------------------------------------------

    nvar_vdiff = 7 + khydromet + ktrac

    iu   = 1;   iv   = 2
    ixl  = 3;   ixi  = 4;  ixv = 5
    itotte= 6;  ithv = 7
    ih   = 8;   iqv  = 9

    !>KF suggestion
    IF((7 + khydromet) > iqv ) &
    CALL finish( TRIM(thismodule),'matrix for vdiff is not properly defined')

    IF(ktrac > 0)  THEN
      itrc_start = 7 + khydromet +1
    ELSE
      itrc_start = 7 + khydromet
    ENDIF
    !<KF

    !-------------------------------------------------------------------
    ! # of vertical levels on which the prognostic equations are solved
    !-------------------------------------------------------------------

    ALLOCATE( ibtm_var(nvar_vdiff),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
      & 'Allocation of ibtm_var failed')

    ! momentum, heat, water substances and tracers are solved on
    ! klev full levels

    ibtm_var(:)    = klev

    ! TTE and the variance of $\theta_v$ are solved at klev-1 half levels.
    ! The upper and lower boundaries of the atmosphere are excluded.

    ibtm_var(itotte) = klev -1
    ibtm_var(ithv) = klev -1

    !------------------------------------------
    ! Set up matrix indices
    !------------------------------------------

    ALLOCATE( matrix_idx(nvar_vdiff),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
      & 'Allocation of matrix_idx failed')

    matrix_idx(iu)   = 1  ; imuv = 1
    matrix_idx(iv)   = 1  ! u and v share the same exchange coeff.
    matrix_idx(ixl)  = 2
    matrix_idx(ixi)  = 2  ! cloud water and ice share the same exchange coeff.
    matrix_idx(ixv)  = 3
    matrix_idx(itotte) = 4
    matrix_idx(ithv) = 5
    matrix_idx(ih)   = 6 ; imh  = 6
    matrix_idx(iqv)  = 7 ; imqv = 7

    IF (ktrac>0) matrix_idx(nvar_vdiff-ktrac+1:nvar_vdiff) = 2

    nmatrix = 7    ! total number of matrices

    !---------------------------------------------------------------------------
    ! # of vertical levels on which elimination of the coefficient will be done
    !---------------------------------------------------------------------------

    ALLOCATE( ibtm_mtrx(nmatrix),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
    & 'Allocation of ibtm_mtrx failed')

    ibtm_mtrx(:)                = klev
    ibtm_mtrx(matrix_idx(itotte)) = klev -1
    ibtm_mtrx(matrix_idx(ithv)) = klev -1

    !$ACC ENTER DATA COPYIN(matrix_idx, ibtm_mtrx, ibtm_var)

  END SUBROUTINE init_vdiff_solver
  !-------------
  !>
  !!
  SUBROUTINE cleanup_vdiff_solver

    INTEGER :: ist

    !$ACC EXIT DATA DELETE(matrix_idx, ibtm_mtrx, ibtm_var)
    DEALLOCATE( matrix_idx,ibtm_mtrx,ibtm_var, STAT=ist)
    IF (ist/=SUCCESS) CALL finish('cleanup_vdiff_solver','Deallocation failed')

  END SUBROUTINE cleanup_vdiff_solver
  !-------------
  !>
  !!
  !! Set up coeffient matrix of the linear algebraic system and 
  !! perform Gauss elimination. For moisture, the last row of the 
  !! matrix (aa_btm) can not be finished yet because the evapotranspiration 
  !! coefficients "cair" and "csat" are not yet available. Thus for this variable 
  !! elimination is performed only till level klev-1. 
  !! 
  SUBROUTINE matrix_setup_elim( jcs, kproma, kbdim, klev, klevm1, &! in
                              & ksfc_type, itop,              &! in
                              & pcfm, pcfh, pcfh_tile, pcfv,  &! in
                              & pcftotte, pcfthv,             &! in
                              & pprfac,                       &! in
                              & prmairm, prmairh, prmrefm,    &! in
                              & aa, aa_btm                    )! out
    ! Arguments

    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, klevm1, ksfc_type
    INTEGER, INTENT(IN) :: itop

    REAL(wp),INTENT(IN) :: pcfm     (:,:)   !< (kbdim,klev) exchange coeff. for u, v
    REAL(wp),INTENT(IN) :: pcfh     (:,:)   !< (kbdim,klevm1) exchange coeff. for heat and tracers
    REAL(wp),INTENT(IN) :: pcfh_tile(:,:)   !< (kbdim,ksfc_type) exchange coeff. for heat and qv, at surface
    REAL(wp),INTENT(IN) :: pcfv     (:,:)   !< (kbdim,klev) exchange coeff. for total water variance
    REAL(wp),INTENT(IN) :: pcftotte (:,:)   !< (kbdim,klev) exchange coeff. for TTE
    REAL(wp),INTENT(IN) :: pcfthv   (:,:)   !< (kbdim,klev) exchange coeff. for variance of theta_v
    REAL(wp),INTENT(IN) :: pprfac   (:,:)   !< (kbdim,klev) prefactor for the exchange coefficients
    REAL(wp),INTENT(IN) :: prmairm  (:,:)   !< (kbdim,klev) reciprocal of layer air mass, full levels
    REAL(wp),INTENT(IN) :: prmairh  (:,:)   !< (kbdim,klevm1) reciprocal of layer air mass, half levels
    REAL(wp),INTENT(IN) :: prmrefm  (:,:)   !< (kbdim,klev) reciprocal of layer ref air mass, full levels

    REAL(wp),INTENT(OUT) :: aa    (:,:,:,:)     !< (kbdim,klev,3,nmatrix) exchange coeff. matrices    out
    REAL(wp),INTENT(OUT) :: aa_btm(:,:,:,imh:)  !< (kbdim,3,ksfc_type,imh:imqv) out
                                                !< last (the klev-th) row of the coeff. matrices 
                                                !< of dry static energy and moisture

    ! Local variables

    REAL(wp) :: zkstar (kbdim,itop-1:klev)     !< scaled exchange coeff on half-levels
    REAL(wp) :: zkh    (kbdim,itop-1:klevm1)   !< scaled exchange doeff on full-levels, 
                                               !< for TTE and variance of theta_v
    INTEGER  :: im             !< index of coefficient matrix
    INTEGER  :: jc, jk, jsfc   !< loop indices
    INTEGER  :: jkm1, jmax

    !$ACC DATA &
    !---- Argument arrays - intent(in)
    !$ACC PRESENT(pcfm,pcfh,pcfh_tile,pcfv,pcftotte,pcfthv,pprfac,prmairm) &
    !$ACC PRESENT(prmairh,prmrefm) &
    !---- Argument arrays - intent(inout)
    !$ACC PRESENT(aa,aa_btm) &
    !---- Local Variables
    !$ACC CREATE(zkstar,zkh) &
    !---- module variable
    !$ACC PRESENT(ibtm_mtrx)

    !-----------------------------------------------------------------------
    ! For all prognostic variables: no turbulent flux at the upper boundary
    !-----------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      zkstar(jc,itop-1) = 0._wp
    ENDDO
    !$ACC END PARALLEL

    !-----------------------------------------------------------------------
    ! For momentum: surface flux is considered
    !-----------------------------------------------------------------------
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klev
      DO jc = jcs,kproma
         zkstar(jc,jk) = pprfac(jc,jk)*pcfm(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    im = matrix_idx(iu)    ! also = matrix_idx(iv)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klev
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmairm(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL


    !---------------------------------------------------------------------
    ! Dry static energy: surface fluxes on different surface types 
    ! are handled separately. 
    !---------------------------------------------------------------------
    im = imh
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
        zkstar(jc,jk) =  pprfac(jc,jk) &
                                 &   *pcfh(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmairm(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! Set the bottom row of the coeff matrix. The same formula applies 
    ! for all surface types (land, water, ice).

    jk = klev
    IF ( isrfc_type == 1 ) THEN
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jc = jcs,kproma
          aa_btm(jc,1,jsfc,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)    ! -K*_{k-1/2}/dm_k
          aa_btm(jc,3,jsfc,im) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prmairm(jc,jk)
          aa_btm(jc,2,jsfc,im) = 1._wp - aa_btm(jc,1,jsfc,im)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jc = jcs,kproma
          aa_btm(jc,1,jsfc,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)    ! -K*_{k-1/2}/dm_k
          aa_btm(jc,3,jsfc,im) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prmairm(jc,jk)
          aa_btm(jc,2,jsfc,im) = 1._wp - aa_btm(jc,1,jsfc,im) - aa_btm(jc,3,jsfc,im)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    END IF

    !---------------------------------------------------------------------
    ! Moisture: different surface types are handled separately.
    !---------------------------------------------------------------------
    im = imqv
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
      zkstar(jc,jk) =  pprfac(jc,jk) &
                                 &   *pcfh(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmrefm(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmrefm(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! Bottom row of the matrix: finish the setup over water and ice;
    ! do part of the computation for land. Later in subroutine
    ! matrix_to_richtmyer_coeff, aa_btm(:,3,idx_land,imqv) will be 
    ! modified, and aa_btm(:,2,idx_land,imqv) re-computed.

    jk = klev
    IF ( isrfc_type == 1 ) THEN
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jc = jcs,kproma
          aa_btm(jc,1,jsfc,im) = -zkstar(jc,jk-1)*prmrefm(jc,jk)    ! -K*_{k-1/2}/dm_k
          aa_btm(jc,3,jsfc,im) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prmrefm(jc,jk)
          aa_btm(jc,2,jsfc,im) = 1._wp - aa_btm(jc,1,jsfc,im)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jc = jcs,kproma
          aa_btm(jc,1,jsfc,im) = -zkstar(jc,jk-1)*prmrefm(jc,jk)    ! -K*_{k-1/2}/dm_k
          aa_btm(jc,3,jsfc,im) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prmrefm(jc,jk)
          aa_btm(jc,2,jsfc,im) = 1._wp - aa_btm(jc,1,jsfc,im) - aa_btm(jc,3,jsfc,im)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    END IF

    !----------------------------------------------------------------------
    ! For all advected tracers except water vapour: no turbulent flux at 
    ! the surface.
    !----------------------------------------------------------------------
    !im = matrix_idx(ixl)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      zkstar(jc,klev) = 0._wp  ! lower boundary, no turbulent flux
    END DO
    !$ACC END PARALLEL

    im = matrix_idx(ixl)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klev
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmrefm(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmrefm(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !----------------------------------------------------------------------
    ! For total water variance: no surface flux. The exchange coefficient
    ! pcfv has been set to to zero in subroutine sfc_exchange_coeff, which
    ! automatically leads to zkstar(:,klev) = 0._wp, thus no additional 
    ! attention is needed here.
    !----------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klev
      DO jc = jcs,kproma
      zkstar(jc,jk) =  pprfac(jc,jk) &
                                 &   *pcfv(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !zkstar(1:kproma,itop:klev) = pprfac(1:kproma,itop:klev) &
    !                           &  *pcfv(1:kproma,itop:klev)

    im = matrix_idx(ixv)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klev
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmrefm(jc,jk)
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmrefm(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !------------------------------------------------------------------------
    ! For TTE: Note that
    ! - Vertical averaging is needed to convert exchange coefficient from
    !   half to full levels, because TTE equation is solved on half levels.
    ! - TTE equation is solved only till array subscript klevm1, which
    !   corresponds to half level (klev - 1/2), i.e., the lowest
    !   interface above surface. Surface value of TTE is (already)
    !   computed in subroutine "sfc_exchange_coeff".
    !------------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klev
      DO jc = jcs,kproma
      zkstar(jc,jk) =  pprfac(jc,jk) &
                                 &   *pcftotte(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
        zkh(jc,jk) = 0.5_wp*(zkstar(jc,jk)+zkstar(jc,jk+1))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      zkh(jc,itop-1) = 0._wp  ! upper boundary, no flux
    ENDDO
    !$ACC END PARALLEL

    im = matrix_idx(itotte)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkh(jc,jk-1)*prmairh(jc,jk)
        aa(jc,jk,3,im) = -zkh(jc,jk  )*prmairh(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !------------------------------------------------
    ! For the variance of theta_v (similar to TTE)
    !------------------------------------------------
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klev
      DO jc = jcs,kproma
      zkstar(jc,jk) =  pprfac(jc,jk) &
                                 &   *pcfthv(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
        zkh(jc,jk) = 0.5_wp*(zkstar(jc,jk)+zkstar(jc,jk+1))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    im = matrix_idx(ithv)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkh(jc,jk-1)*prmairh(jc,jk)
        aa(jc,jk,3,im) = -zkh(jc,jk  )*prmairh(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !-----------------------------------------------------------------------------
    ! Gauss elimination for the coefficient matrices at
    ! - vertical levels [itop,klev-2], for TTE and variance of theta_v;
    ! - vertical levels [itop,klev-1], for all the other variables.
    !-----------------------------------------------------------------------------

#ifndef _OPENACC
    DO im = 1,nmatrix
      DO jc = jcs,kproma
        aa(jc,itop,3,im) = aa(jc,itop,3,im)/aa(jc,itop,2,im)
      ENDDO

      jmax = ibtm_mtrx(im) - 1
      DO jk = itop+1,jmax
        jkm1 = jk - 1
        DO jc = jcs,kproma
          aa(jc,jk,2,im) =  aa(jc,jk,2,im)                       &
                            & -aa(jc,jk,1,im)*aa(jc,jkm1,3,im)
          aa(jc,jk,3,im) =  aa(jc,jk,3,im)/aa(jc,jk,2,im)
        ENDDO
      ENDDO
    END DO
#else
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO im = 1,nmatrix
      DO jc = jcs,kproma

        jmax = ibtm_mtrx(im) - 1
        aa(jc,itop,3,im) = aa(jc,itop,3,im)/aa(jc,itop,2,im)
        !$ACC LOOP SEQ
        DO jk = itop+1,jmax
          jkm1 = jk - 1
          aa(jc,jk,2,im) =  aa(jc,jk,2,im)                       &
                            & -aa(jc,jk,1,im)*aa(jc,jkm1,3,im)
          aa(jc,jk,3,im) =  aa(jc,jk,3,im)/aa(jc,jk,2,im)
        ENDDO
      ENDDO
    END DO
    !$ACC END PARALLEL
#endif


    ! Translation for developers who prefer to think in terms of 
    ! the Richtmyer-Morthon formula and are familiar with the paper by
    ! Polcher et al (1998): after this elimination, 
    !  aa(:,1:ibtm_mtrx(im)-1,2,:) becomes C  (Eqn. 17),
    !  aa(:,1:ibtm_mtrx(im)-1,3,:) becomes -A (Eqn. 19).
    ! See subroutine matrix_to_richtmyer_coeff.

  !$ACC WAIT
  !$ACC END DATA

  END SUBROUTINE matrix_setup_elim

  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  SUBROUTINE rhs_setup( jcs, kproma, kbdim, itop, klev, klevm1,&! in
                      & ksfc_type, ktrac, pdtime,            &! in
                      & pum1, pvm1, pcptgz, pqm1,            &! in
                      & pxlm1, pxim1, pxvar, pxtm1, pxt_emis,&! in
                      & prmrefm, ptottevn, pzthvvar, aa,     &! in
                      & bb, bb_btm                           )! out

    ! Arguments

    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, itop, klev, klevm1
    INTEGER, INTENT(IN) :: ksfc_type, ktrac
    REAL(wp),INTENT(IN) :: pdtime

    REAL(wp),INTENT(IN) :: pum1     (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pvm1     (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pcptgz   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pqm1     (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxlm1    (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxim1    (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxvar    (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxtm1    (:,:,:) !< (kbdim,klev,ktrac)
    REAL(wp),INTENT(IN) :: pxt_emis (:,:)   !< (kbdim,ktrac)
   !REAL(wp),INTENT(IN) :: pxt_emis (:,:,:) ! (kbdim,klev,ktrac) backup for later use
    REAL(wp),INTENT(IN) :: ptottevn (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pzthvvar (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: prmrefm  (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: aa       (:,:,:,:) !< (kbdim,klev,3,nmatrix)

    REAL(wp),INTENT(OUT) :: bb    (:,:,:)   !< (kbdim,klev,nvar_vdiff) OUT
    REAL(wp),INTENT(OUT) :: bb_btm(:,:,ih:) !< (kbdim,ksfc_type,ih:iqv) OUT

    ! Local variables

    REAL(wp) :: ztmp(kbdim,klev)
    INTEGER  :: jsfc, jt, irhs, im, jk, jc

    !$ACC DATA PRESENT( pxtm1, pxt_emis ) IF( ktrac > 0 )
    !$ACC DATA &
    !---- Argument arrays - intent(in)
    !$ACC PRESENT(pum1,pvm1,pcptgz,pqm1,pxlm1,pxim1,pxvar) &
    !$ACC PRESENT(ptottevn,pzthvvar,prmrefm,aa) &
    !---- Argument arrays - intent(inout)
    !$ACC PRESENT(bb,bb_btm) &
    !---- Local Variables
    !$ACC CREATE(ztmp)

    !-------------------------------------------------------------------
    ! First handle variables that are defined on full levels
    !-------------------------------------------------------------------
    ! u and v
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klev
      DO jc = jcs,kproma
        bb(jc,jk,iu) = pum1(jc,jk)
        bb(jc,jk,iv) = pvm1(jc,jk)

    ! Hydrometeors and the variance of cloud droplets

        bb(jc,jk,ixl) = pxlm1(jc,jk)
        bb(jc,jk,ixi) = pxim1(jc,jk)
        bb(jc,jk,ixv) = pxvar(jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL

    ! Other tracers

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( ktrac > 0 )
    !$ACC LOOP SEQ
    DO jt = 1,ktrac
      irhs = jt - 1 + itrc_start
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = itop,klev
          DO jc = jcs,kproma
            bb(jc,jk,irhs) =  pxtm1(jc,jk,jt)
          END DO
        END DO
       !bb(1:kproma,itop:klev,irhs) =  pxtm1(1:kproma,itop:klev,jt)
    ENDDO
    !$ACC END PARALLEL

    ! Heat and moisture

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
        bb(jc,jk,ih ) = pcptgz(jc,jk)
        bb(jc,jk,iqv) = pqm1  (jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jsfc = 1,ksfc_type
      DO jc = jcs,kproma 
        bb_btm(jc,jsfc,ih)  = pcptgz(jc,klev)
        bb_btm(jc,jsfc,iqv) =   pqm1(jc,klev)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !-------------------------------------------------------------------
    ! TTE and the variance of theta_v:
    ! These variables are defined at half levels. Array index jk
    ! correspond to half level k+1/2. Thus klev correspond to the
    ! lower boundary. The linear solver only solves till index klevm1.
    !-------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
        bb(jc,jk,itotte) =  ptottevn(jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL

    im = matrix_idx(itotte)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      bb(jc,     klevm1,itotte) =  bb(jc,klevm1,itotte)   &
                              & -aa(jc,klevm1,3,im)   &
                              & *ptottevn(jc,klev)
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jc = jcs,kproma
        bb(jc,jk,ithv) =  pzthvvar(jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL

    im = matrix_idx(ithv)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      bb(jc,     klevm1,ithv) =  bb(jc,klevm1,ithv)   &
                              & -aa(jc,klevm1,3,im)   &
                              & *pzthvvar(jc,klev)
    ENDDO
    !$ACC END PARALLEL

    !--------------------------------------------------------------------
    ! Apply the implicitness factor
    !--------------------------------------------------------------------
    !bb     = tpfac2*bb
    !bb_btm = tpfac2*bb_btm

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jt = 1, itotte-1
      DO jk = 1,klev
        DO jc = jcs,kproma
          bb(jc,jk,jt)  = tpfac2*bb(jc,jk,jt)
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jt = itotte, iqv
      DO jk = 1,klevm1
        DO jc = jcs,kproma
          bb(jc,jk,jt)  = tpfac2*bb(jc,jk,jt)
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    IF (ktrac>0) THEN
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO jt = itrc_start, nvar_vdiff
        DO jk = 1,klev
          DO jc = jcs,kproma
            bb(jc,jk,jt)  = tpfac2*bb(jc,jk,jt)
          ENDDO
        ENDDO
      ENDDO
    !$ACC END PARALLEL

    ENDIF

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jt = ih,iqv
      DO jk = 1,ksfc_type
        DO jc = jcs,kproma
          bb_btm(jc,jk,jt)  = tpfac2*bb_btm(jc,jk,jt)
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !--------------------------------------------------------------------
    ! Add tracer emissions
    !--------------------------------------------------------------------
    ! Currently we follow ECHAM in which only the surface emission
    ! is treated in "vdiff".
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      ztmp(jc,klev) = prmrefm(jc,klev)*pdtime
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( ktrac > 0 )
    !$ACC LOOP SEQ
    DO jt = 1,ktrac
      irhs = jt - 1 + itrc_start
        !$ACC LOOP GANG VECTOR
        DO jc = jcs,kproma
          bb(jc,klev,irhs) =         bb(jc,klev,irhs) &
                              & + pxt_emis(jc,jt)        &
                              &      *ztmp(jc,klev)
        END DO
    ENDDO
    !$ACC END PARALLEL


    !DO jt = 1,ktrac
    !   irhs = jt - 1 + itrc_start
    !   bb(1:kproma,klev,irhs) =         bb(1:kproma,klev,irhs) &
    !                          & + pxt_emis(1:kproma,jt)        &
    !                          &      *ztmp(1:kproma,klev)
    !ENDDO

    ! Later we may consider treating emission on all vertical levels
    ! in the same way.
    !
    !ztmp(jcs:kproma,itop:klev) = prmrefm(jcs:kproma,itop:klev)*pdtime
    !
    !DO jt = 1,ktrac
    !   irhs = jt - 1 + itrc_start
    !   bb(jcs:kproma,itop:klev,irhs) =         bb(jcs:kproma,itop:klev,irhs) &
    !                               & + pxt_emis(jcs:kproma,itop:klev,jt)   &
    !                               &      *ztmp(jcs:kproma,itop:klev)
    !ENDDO
  !$ACC WAIT
  !$ACC END DATA
  !$ACC END DATA


  END SUBROUTINE rhs_setup

  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  !!
  !! Gauss elimination of the right-hand-side vector 
  !! using coefficients obtained in subroutine "matrix_setup_elim".
  !!
  !! Translation for developers who prefer to think in terms of 
  !! the Richtmyer-Morthon formula and are familiar with the paper by
  !! Polcher et al (1998): after the elimination at the end of 
  !! subroutine matrix_setup_elim, aa(:,1:ibtm_mtrx(im)-1,2,:) 
  !! became the coeff C defined by Eqn. 17 of Polcher et al (1998). 
  !! It is used in this subroutine to convert the variable bb 
  !! into the Richtmyer coeff B (cf Eqn. 19 of Polcher et al 1998).
  !!
  SUBROUTINE rhs_elim( jcs, kproma, kbdim, itop, klev, klevm1, &! in
                     & aa, bb                             )! in, inout

    INTEGER, INTENT(IN)    :: jcs, kproma, kbdim, itop, klev, klevm1
    REAL(wp),INTENT(IN)    :: aa(:,:,:,:) !< (kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT) :: bb(:,:,:)   !< (kbdim,klev,nvar_vdiff)

    REAL(wp) :: znum, zden
    INTEGER  :: jvar, im, jk, jkm1, jmax, jc

    !$ACC DATA &
    !---- Argument arrays - intent(in)
    !$ACC PRESENT(aa) &
    !---- Argument arrays - intent(inout)
    !$ACC PRESENT(bb) &
    !---- Global Variables
    !$ACC COPYIN(matrix_idx,ibtm_var)

    ! 1. Vertical levels [itop+1,klev-2] for TTE and variance of theta_v;
    !    [itop+1,klev-1] for all the other variables.

#ifndef _OPENACC
    DO jvar = 1,nvar_vdiff
      im = matrix_idx(jvar)  ! Index of coefficient matrix
      DO jc = jcs, kproma
        bb(jc,itop,jvar) =  bb(jc,itop,jvar)/aa(jc,itop,2,im)
      ENDDO

      jmax = ibtm_var(jvar) - 1

      DO jk = itop+1,jmax
        jkm1 = jk - 1
        DO jc = jcs,kproma
          znum =  bb(jc,jk  ,jvar)                     &
                   & -bb(jc,jkm1,jvar)*aa(jc,jk,1,im)
          bb(jc,jk,jvar) = znum/aa(jc,jk,2,im)
        ENDDO
      ENDDO
    ENDDO !jvar: variable loop
#else

  !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO jvar = 1,nvar_vdiff
    DO jc = jcs,kproma

      im   = matrix_idx(jvar)
      jmax =   ibtm_var(jvar) - 1

      bb(jc,itop,jvar) =  bb(jc,itop,jvar)/aa(jc,itop,2,im)
      !$ACC LOOP SEQ
      DO jk = itop+1,jmax
        jkm1 = jk - 1

        znum           = bb(jc,jk,jvar) - bb(jc,jkm1,jvar)*aa(jc,jk,1,im)
        bb(jc,jk,jvar) = znum/aa(jc,jk,2,im)
      ENDDO
    ENDDO
  END DO
  !$ACC END PARALLEL

#endif

    ! 2. Bottom level for all variables except u, v, dry static energy
    !    and moisture. After this step the array bb contains the
    !    solution of the linear system.

    DO jvar = 1,nvar_vdiff

      IF (jvar==iu.OR.jvar==iv.OR.jvar==ih.OR.jvar==iqv ) THEN
         CYCLE
      ELSE

      im   = matrix_idx(jvar)  ! Index of coefficient matrix
      jk   = ibtm_var(jvar)    ! Bottom level index
      jkm1 = jk - 1

      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs,kproma
        zden =  aa(jc,jk,2,im)                      &
                      & -aa(jc,jk,1,im)*aa(jc,jkm1,3,im)
        znum =  bb(jc,jk,jvar)                      &
                      & -aa(jc,jk,1,im)*bb(jc,jkm1,jvar)
        bb(jc,jk,jvar) = znum/zden
      ENDDO
      !$ACC END PARALLEL

      END IF
    ENDDO !jvar: variable loop

    ! Note that for TTE and the variance of theta_v, klev-1 is the lowest
    ! level above surface. Now set boundary condition for the variance 
    ! of theta_v.

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      bb(jc,klev,ithv) = bb(jc,klevm1,ithv)
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE rhs_elim

  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  !!
  !! Prepare the Richtmyer-Morton coeffcients for dry static energy and 
  !! moisture, to be used by the surface models (ocean, sea-ice, land).
  !!
  SUBROUTINE matrix_to_richtmyer_coeff( jg, jcs, kproma, kbdim, klev, ksfc_type, idx_lnd, &! in
                                      & aa, bb,                                      &! in
                                      & pdtime, delz,                                &! in
                                      & aa_btm, bb_btm,                              &! inout
                                      & pen_h, pfn_h, pen_qv, pfn_qv,                &! out
                                      & pcair,                                       &! in
                                      & pcsat)                                        ! in

    INTEGER,INTENT(IN)     :: jg, jcs, kproma, kbdim, klev, ksfc_type, idx_lnd
    REAL(wp),INTENT(IN)    :: aa    (:,:,:,imh:) !< (kbdim,klev,3,imh:imqv)
    REAL(wp),INTENT(IN)    :: bb    (:,:,ih:)    !< (kbdim,klev,ih:iqv)
    REAL(wp),INTENT(IN)    :: pdtime 
    REAL(wp),INTENT(IN)    :: delz(:)            !< (kbdim)
    REAL(wp),INTENT(INOUT) :: aa_btm(:,:,:,imh:) !< (kbdim,3,ksfc_type,imh:imqv)
    REAL(wp),INTENT(INOUT) :: bb_btm(:,:,ih:)    !< (kbdim,ksfc_type,ih:iqv)

    REAL(wp),INTENT(OUT) :: pen_h (:,:)  !< (kbdim,ksfc_type) OUT
    REAL(wp),INTENT(OUT) :: pfn_h (:,:)  !< (kbdim,ksfc_type) OUT
    REAL(wp),INTENT(OUT) :: pen_qv(:,:)  !< (kbdim,ksfc_type) OUT
    REAL(wp),INTENT(OUT) :: pfn_qv(:,:)  !< (kbdim,ksfc_type) OUT

    REAL(wp),OPTIONAL,INTENT(IN)    :: pcair(:) !< (kbdim)
    REAL(wp),OPTIONAL,INTENT(IN)    :: pcsat(:) !< (kbdim)

    INTEGER  :: jk, jsfc, klevm1

    !$ACC DATA PRESENT( aa, bb, delz, aa_btm, bb_btm, pen_h, pfn_h, pen_qv, pfn_qv )
    !$ACC DATA PRESENT( pcair, pcsat ) IF( PRESENT(pcair) )

    klevm1 = klev - 1

    !---------------------------------------------------------
    ! Matrix setup and bottom level elimination for moisture
    !---------------------------------------------------------
    ! Evapotranspiration has to be considered over land 

    IF (echam_phy_config(jg)%ljsb .AND. idx_lnd<=ksfc_type) THEN

      jsfc = idx_lnd

      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jk = jcs, kproma
        aa_btm(jk,2,jsfc,imqv) =           1._wp - aa_btm(jk,1,jsfc,imqv) &
                                   & - pcair(jk)*aa_btm(jk,3,jsfc,imqv)
        aa_btm(jk,3,jsfc,imqv) =   pcsat(jk)*aa_btm(jk,3,jsfc,imqv)
      END DO
      !$ACC END PARALLEL

    END IF ! ljsbach

    ! Bottom level elimination for all surface types

    IF ( isrfc_type == 1) THEN
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jk = jcs, kproma
          aa_btm(jk,2,jsfc,imqv) =  aa_btm(jk,2,jsfc,imqv)  &
                                 & -aa_btm(jk,1,jsfc,imqv)  &
                                 & *aa    (jk,klevm1,3,imqv)

          aa_btm(jk,3,jsfc,imqv) =  -lhflx*pdtime/delz(jk) &
                                 & /aa_btm(jk,2,jsfc,imqv)

          bb_btm(jk,jsfc,iqv)    = (bb_btm(jk,jsfc,iqv)    &
                                 & -aa_btm(jk,1,jsfc,imqv) &
                                 & *bb    (jk,klevm1,iqv) )&
                                 & /aa_btm(jk,2,jsfc,imqv)

        END DO
      END DO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jk = jcs, kproma
          aa_btm(jk,2,jsfc,imqv) =  aa_btm(jk,2,jsfc,imqv)  &
                                 & -aa_btm(jk,1,jsfc,imqv)  &
                                 & *aa    (jk,klevm1,3,imqv)

          aa_btm(jk,3,jsfc,imqv) =  aa_btm(jk,3,jsfc,imqv)  &
                                 & /aa_btm(jk,2,jsfc,imqv)

          bb_btm(jk,jsfc,iqv)    = (bb_btm(jk,jsfc,iqv)    &          
                                 & -aa_btm(jk,1,jsfc,imqv) &
                                 & *bb    (jk,klevm1,iqv) )&
                                 & /aa_btm(jk,2,jsfc,imqv)
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

    !---------------------------------------------------------
    ! Bottom level elimination for dry static energy
    !---------------------------------------------------------
    IF ( isrfc_type == 1 ) THEN
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jk = jcs, kproma
          aa_btm(jk,2,jsfc,imh) =  aa_btm(jk,2,jsfc,imh) &
                                      & -aa_btm(jk,1,jsfc,imh) &
                                      & *aa    (jk,klevm1,3,imh)

          aa_btm(jk,3,jsfc,imh) =  -shflx*cpd*pdtime/delz(jk) &
                                      & /aa_btm(jk,2,jsfc,imh)

          bb_btm(jk,jsfc,ih)    = (bb_btm(jk,jsfc,ih)    &
                                      & -aa_btm(jk,1,jsfc,imh) &
                                      & *bb    (jk,klevm1,ih) )&
                                      & /aa_btm(jk,2,jsfc,imh)
        END DO
      END DO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jk = jcs, kproma

          aa_btm(jk,2,jsfc,imh) =  aa_btm(jk,2,jsfc,imh) &
                                      & -aa_btm(jk,1,jsfc,imh) &
                                      & *aa    (jk,klevm1,3,imh)

          aa_btm(jk,3,jsfc,imh) =  aa_btm(jk,3,jsfc,imh) &
                                      & /aa_btm(jk,2,jsfc,imh)

          bb_btm(jk,jsfc,ih)    = (bb_btm(jk,jsfc,ih)    &          
                                      & -aa_btm(jk,1,jsfc,imh) &
                                      & *bb    (jk,klevm1,ih) )&
                                      & /aa_btm(jk,2,jsfc,imh)
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

    !---------------------------------------------------------
    ! Convert matrix entries to Richtmyer-Morton coefficients
    !---------------------------------------------------------

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jsfc = 1,ksfc_type
      DO jk = jcs, kproma
        pen_h (jk,jsfc) = -aa_btm(jk,3,jsfc,imh)
        pen_qv(jk,jsfc) = -aa_btm(jk,3,jsfc,imqv)

        pfn_h (jk,jsfc) =  bb_btm(jk,jsfc,ih )*tpfac1
        pfn_qv(jk,jsfc) =  bb_btm(jk,jsfc,iqv)*tpfac1
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA
    !$ACC END DATA

  END SUBROUTINE matrix_to_richtmyer_coeff
  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  !!
  !! Do Back-substitution to get the solution of the linear system.
  !!
  !! Translation for developers who prefer to think in terms of 
  !! the Richtmyer-Morthon formula and are familiar with the paper by
  !! Polcher et al (1998): on entry bb contains the solution 
  !! at the bottom level and the coeff B on upper levels;
  !! On exit it becomes the solution of the linear system.
  !! aa(:,:,3,:) used here corresponds to -A in the Appendix of 
  !! Polcher et al (1998).
  !! Note that VDIFF uses the implicit time stepping as in IFS
  !! in contrast to Polcher et al (1998). Thus the solution is 
  !! not yet the new value at time step t+dt. 
  !!
  SUBROUTINE rhs_bksub( jcs, kproma, kbdim, itop, klev, aa, bb )

    INTEGER, INTENT(IN)   :: jcs, kproma, kbdim, itop, klev
    REAL(wp),INTENT(IN)   :: aa(:,:,:,:) !< (kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT):: bb(:,:,:)   !< (kbdim,klev,nvar_vdiff)

    INTEGER  :: jvar, im, jk, jkp1, jl, jmax
    !$ACC DATA &
    !---- Argument arrays - intent(in)
    !$ACC PRESENT(aa) &
    !---- Argument arrays - intent(out)
    !$ACC PRESENT(bb) &
    !---- Argument arrays - Module Variables
    !$ACC COPYIN(matrix_idx,ibtm_var)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jvar = 1,nvar_vdiff
      DO jl = jcs,kproma

        im   = matrix_idx(jvar)
        jmax =   ibtm_var(jvar) - 1

        !$ACC LOOP SEQ
        DO jk = jmax,itop,-1
          jkp1 = jk + 1

          bb(jl,jk,jvar) =  bb(jl,jk ,jvar) &
                               & -bb(jl,jkp1,jvar) &
                               & *aa(jl,jk  ,3,im)
        ENDDO
      ENDDO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE rhs_bksub
  !-------------
  !>
  !!
  SUBROUTINE vdiff_tendencies( jg, jcs, kproma, kbdim, itop, klev, klevm1, &! in
                             & ktrac, ksfc_type, idx_wtr,                  &! in
                             & pdtime,                                     &! in
                             & pum1, pvm1, ptm1,                           &! in
                             & pmair,                                      &! in
!!$                             & pmref,                                      &! in
                             & pqm1, pxlm1, pxim1, pxtm1,                  &! in
                             & pgeom1, pcptgz,                             &! in
                             & pztottevn, pzthvvar,                        &! in
                             & pcfm_tile, pfrc, bb,                        &! in
                             & pkedisp,                                    &! out
                             & pxvar, pz0m_tile,                           &! inout
                             & pute_vdf, pvte_vdf, pq_vdf,                 &! out
                             & pqte_vdf, pxlte_vdf, pxite_vdf, pxtte_vdf,  &! out
                             & pz0m, ptotte, pthvvar                       )! out
!!$                             & pz0m, ptotte, pthvvar,                      &! out
!!$                             & psh_vdiff,pqv_vdiff                         )! out

    INTEGER, INTENT(IN) :: jg, jcs, kproma, kbdim, itop, klev, klevm1, ktrac !!$, klevp1
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr
    REAL(wp),INTENT(IN) :: pdtime

    REAL(wp),INTENT(IN)  :: pum1   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pvm1   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: ptm1   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pmair  (:,:)   !< (kbdim,klev) moist air mass [kg/m2]
!!$    REAL(wp),INTENT(IN)  :: pmref  (:,:)   !< (kbdim,klev) dry   air mass [kg/m2]
    REAL(wp),INTENT(IN)  :: pqm1   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxlm1  (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxim1  (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxtm1  (:,:,:) !< (kbdim,klev,ktrac)
    REAL(wp),INTENT(IN)  :: pgeom1 (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pcptgz (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pztottevn(:,:) !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pzthvvar(:,:) !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pcfm_tile     (:,:) !< (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)  :: pfrc          (:,:) !< (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)  :: bb            (:,:,:) !<(kbdim,klev,nvar_vdiff)

    REAL(wp),INTENT(OUT) :: pkedisp(:) !< (kbdim) vertically integrated dissipation
                                       !  of kinetic energy [W/m2]

    REAL(wp),INTENT(INOUT) :: pxvar    (:,:) !< (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pz0m_tile(:,:) !< (kbdim,ksfc_type)

    REAL(wp),INTENT(OUT) :: pute_vdf (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pvte_vdf (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pq_vdf   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pqte_vdf (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxlte_vdf(:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxite_vdf(:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxtte_vdf(:,:,:) !< (kbdim,klev,ktrac)

    REAL(wp),INTENT(OUT) :: pz0m     (:)   !< (kbdim)
    REAL(wp),INTENT(OUT) :: ptotte   (:,:) !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pthvvar  (:,:) !< (kbdim,klev)
!!$    REAL(wp),INTENT(OUT) :: psh_vdiff(:)   !< (kbdim)
!!$    REAL(wp),INTENT(OUT) :: pqv_vdiff(:)   !< (kbdim)

    REAL(wp) :: ztest, zrdt
    REAL(wp) :: zunew, zvnew, zqnew, zsnew, zhnew
    REAL(wp) :: zcp
    REAL(wp) :: zdis  (kbdim,klev)
    REAL(wp) :: z0m_min

    INTEGER  :: jk, jl, jt, irhs, jsfc

    !-------------------------------------------------------------------
    ! Start GPU data region
    !-------------------------------------------------------------------
    !$ACC DATA PRESENT( pxtm1, pxtte_vdf ) IF( ktrac > 0 )
    !$ACC DATA &
    !---- Argument arrays - intent(in)
    !$ACC PRESENT(pum1,pvm1,ptm1,pmair,pqm1,pxlm1,pxim1) &
!!$    !$ACC PRESENT(pum1,pvm1,ptm1,pmair,pmref,pqm1,pxlm1,pxim1) &
    !$ACC PRESENT(pgeom1,pcptgz,pztottevn,pzthvvar,pcfm_tile,pfrc,bb) &
    !---- Argument arrays - intent(inout)
    !$ACC PRESENT(pxvar,pz0m_tile) &
    !---- Argument arrays - intent(out)
    !$ACC PRESENT(pute_vdf,pvte_vdf,pq_vdf,pqte_vdf,pxlte_vdf,pxite_vdf) &
    !$ACC PRESENT(pz0m,ptotte,pthvvar) &
!!$    !$ACC PRESENT(pz0m,ptotte,pthvvar,psh_vdiff,pqv_vdiff) &
    !$ACC PRESENT(pkedisp) &
    !$ACC CREATE(zdis)

    zrdt   = 1._wp/pdtime

    z0m_min = echam_vdf_config(jg)%z0m_min

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, klev
      DO jl = 1, kbdim
        pute_vdf (jl,jk)   = 0._wp
        pvte_vdf (jl,jk)   = 0._wp
        pq_vdf   (jl,jk)   = 0._wp
        pqte_vdf (jl,jk)   = 0._wp
        pxlte_vdf(jl,jk)   = 0._wp
        pxite_vdf(jl,jk)   = 0._wp
        ptotte     (jl,jk)   = 0._wp
        pthvvar  (jl,jk)   = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( ktrac > 0 )
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jt = 1, ktrac
      DO jk = 1, klev
        DO jl = 1, kbdim
          pxtte_vdf(jl,jk,jt) = 0._wp
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jl = 1, kbdim
      pz0m     (jl)     = 0._wp
    END DO
    !$ACC END PARALLEL
    !-------------------------------------------------------------------
    ! Compute TTE at the new time step.
    !-------------------------------------------------------------------
    ztest = 0._wp
    !$ACC PARALLEL DEFAULT(NONE) REDUCTION(+:ztest) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klevm1
      DO jl = jcs,kproma
        ptotte(jl,jk) = bb(jl,jk,itotte) + tpfac3*pztottevn(jl,jk)
        ztest = ztest+MERGE(1._wp,0._wp,ptotte(jl,jk)<0._wp)
      END DO
    END DO
    !$ACC END PARALLEL

    IF( echam_vdf_config(1)%turb == 2 ) THEN
      ztest = 1._wp
    ELSE
      IF(ztest.NE.0._wp) THEN
        CALL finish('vdiff_tendencies','TTE IS NEGATIVE')
      ENDIF
    ENDIF

    !ptotte(jcs:kproma,klev) = pztottevn(jcs:kproma,klev)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs,kproma
      ptotte(jl,klev) = pztottevn(jl,klev)
    END DO
    !$ACC END PARALLEL

    !-------------------------------------------------------------
    ! Variance of virtual potential temperature
    !-------------------------------------------------------------
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = itop,klev
      DO jl = jcs,kproma
        pthvvar(jl,jk) = bb(jl,jk,ithv) + tpfac3*pzthvvar(jl,jk)
        pthvvar(jl,jk) = MAX(totte_min,pthvvar(jl,jk))
      END DO
    END DO
    !$ACC END PARALLEL
    !-------------------------------------------------------------
    ! Tendency of velocity; kinetic energy dissipation
    !-------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jk = 1,kbdim
      pkedisp(jk) = 0._wp   ! initilize the vertical integral
    END DO
    !$ACC END PARALLEL


    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP SEQ
    DO jk = itop,klev
      !$ACC LOOP GANG VECTOR PRIVATE( zunew, zvnew ) 
      DO jl = jcs,kproma
        pute_vdf(jl,jk) = (bb(jl,jk,iu)-tpfac2*pum1(jl,jk))*zrdt
        pvte_vdf(jl,jk) = (bb(jl,jk,iv)-tpfac2*pvm1(jl,jk))*zrdt

        zunew = bb(jl,jk,iu) + tpfac3*pum1(jl,jk)
        zvnew = bb(jl,jk,iv) + tpfac3*pvm1(jl,jk)

        zdis(jl,jk) = 0.5_wp*( pum1(jl,jk)**2 - zunew**2 &
                    &         +pvm1(jl,jk)**2 - zvnew**2 )
        pkedisp(jl)  = pkedisp(jl) + zdis(jl,jk)*pmair(jl,jk)*zrdt
      END DO
    END DO
    !$ACC END PARALLEL
    !-------------------------------------------------------------
    ! Tendency of T and qv, ql, qi; xvar at the new time step
    !-------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE( zqnew, zsnew, zhnew, zcp )
    DO jk=itop,klev
      DO jl=jcs,kproma

        zqnew = bb(jl,jk,iqv) + tpfac3*pqm1(jl,jk)
        pqte_vdf(jl,jk) = (zqnew-pqm1(jl,jk))*zrdt


        ! The computation of the new temperature must be consistent with the computation of
        ! the static energy pcptgz in the subroutine mo_turbulence_diag:atm_exchange_coeff.
        ! The same specific heat must be used.
        !
        zsnew = bb(jl,jk,ih) + tpfac3*pcptgz(jl,jk)
        zhnew = (zsnew + zdis(jl,jk) - pgeom1(jl,jk))
        zcp   = cpd!+(cpv-cpd)*pqm1(jl,jk) ! cp of moist air
        !
        ! Now derive the heating for constant pressure conditions
        ! as needed for provisional updating in the physics.
        ! 
        pq_vdf(jl,jk)   = (zhnew - ptm1(jl,jk)*zcp)*zrdt*pmair(jl,jk)

        pxlte_vdf(jl,jk) = (bb(jl,jk,ixl) - tpfac2*pxlm1(jl,jk))*zrdt
        pxite_vdf(jl,jk) = (bb(jl,jk,ixi) - tpfac2*pxim1(jl,jk))*zrdt

        pxvar(jl,jk) = bb(jl,jk,ixv) + tpfac3*pxvar(jl,jk)
      END DO
    END DO
    !$ACC END PARALLEL

!!$    IF ( get_lebudget() ) THEN
!!$      psh_vdiff(:) = 0._wp
!!$      pqv_vdiff(:) = 0._wp
!!$      DO jk=itop,klev
!!$        ! compute heat budget diagnostic
!!$        psh_vdiff(jcs:kproma) = psh_vdiff(jcs:kproma) + pmref(jcs:kproma,jk) * &
!!$        & (bb(jcs:kproma,jk,ih)  + (tpfac3 - 1._wp)*pcptgz(jcs:kproma,jk)) * zrdt
!!$        ! compute moisture budget diagnostic
!!$        ! ? zdis appears to be dissipation, probably we don't need this for qv??
!!$        pqv_vdiff(jcs:kproma) = pqv_vdiff(jcs:kproma) + pmref(jcs:kproma,jk)* &
!!$        & (bb(jcs:kproma,jk,iqv) + (tpfac3 - 1._wp)*pqm1(jcs:kproma,jk)) * zrdt
!!$      END DO
!!$    END IF
    !-------------------------------------------------------------
    ! Tendency of tracers
    !-------------------------------------------------------------
!   IF (trlist% anyvdiff /= 0) THEN   ! ECHAM
!     DO 577 jt=1,trlist% ntrac       ! ECHAM
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( ktrac > 0 )
        !$ACC LOOP GANG
        DO jt = 1,ktrac
          irhs = itrc_start + jt - 1
!         IF (trlist% ti(jt)% nvdiff /= 1) CYCLE  ! ECHAM
          !$ACC LOOP
          DO jk = itop,klev
            !$ACC LOOP VECTOR
            DO jl = jcs,kproma
              pxtte_vdf(jl,jk,jt) = (bb(jl,jk,irhs)-tpfac2*pxtm1(jl,jk,jt))*zrdt
            ENDDO
          ENDDO
        ENDDO
        !$ACC END PARALLEL

!577  ENDDO
!     END IF

    !----------------------------------------------------------------------------
    ! Update roughness height over open water, then update the grid-box mean
    !----------------------------------------------------------------------------
    IF (idx_wtr<=ksfc_type) THEN  ! water surface exists in the simulation
      !$ACC PARALLEL DEFAULT(NONE) PRESENT(echam_vdf_config) ASYNC(1)
      !acc loop gang vector
      DO jl = 1,kbdim
        pz0m_tile(jl,idx_wtr) = echam_vdf_config(jg)%z0m_oce
      ENDDO
      !$ACC END PARALLEL
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs,kproma
        IF(pfrc(jl,idx_wtr).GT.0._wp) THEN
          pz0m_tile(jl,idx_wtr) = tpfac1*SQRT( bb(jl,klev,iu)**2+bb(jl,klev,iv)**2 ) &
                                & *pcfm_tile(jl,idx_wtr)*cchar*rgrav
          pz0m_tile(jl,idx_wtr) = MAX(z0m_min,pz0m_tile(jl,idx_wtr))
        ENDIF
      ENDDO
      !$ACC END PARALLEL
    ENDIF

    ! Compute grid-box mean 

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jl = 1,kbdim
      pz0m(jl) = 0._wp
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
      !$ACC LOOP GANG VECTOR
      DO jl = jcs,kproma
        pz0m(jl) = pz0m(jl) + pfrc(jl,jsfc)*pz0m_tile(jl,jsfc)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !-------------------------------------------------------------------
    ! End GPU data region
    !-------------------------------------------------------------------
    !$ACC WAIT
    !$ACC END DATA
    !$ACC END DATA


  END SUBROUTINE vdiff_tendencies
  !-------------



END MODULE mo_vdiff_solver
