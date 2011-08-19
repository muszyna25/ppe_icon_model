!>
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (2010-09)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_vdiff_solver

  USE mo_kind,              ONLY: wp
  USE mo_impl_constants,    ONLY: SUCCESS
  USE mo_exception,         ONLY: finish
#ifdef __ICON__
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_physical_constants,ONLY: grav, vtmpc2, cpd, als, alv
  USE mo_echam_vdiff_params,ONLY: clam, da1, tkemin=>tke_min, cons2, cons25, &
                                & tpfac1, tpfac2, tpfac3, cchar, z0m_min
#else
  USE mo_constants,ONLY: grav=>g, vtmpc2, cpd, als, alv
  USE mo_physc2,   ONLY: clam, da1, tkemin, cons2, cons25, &
                       & tpfac1, tpfac2, tpfac3, cchar, z0m_min
  USE mo_time_control,ONLY: lstart
  USE mo_semi_impl,   ONLY: eps
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_vdiff_solver       !< subroutine
  PUBLIC :: cleanup_vdiff_solver    !< subroutine
  PUBLIC :: matrix_setup_elim       !< subroutine
  PUBLIC :: rhs_setup, rhs_elim     !< subroutines
  PUBLIC :: sfc_solve, rhs_bksub    !< subroutines
  PUBLIC :: vdiff_tendencies        !< subroutine
  PUBLIC :: nvar_vdiff, nmatrix     !< parameters
  PUBLIC :: ih,iqv                  !< parameters

  ! Module variables

  INTEGER :: nvar_vdiff    !< total number of variables affected by turbulent mixing
  INTEGER :: iu, iv, ih, iqv
  INTEGER :: ixl, ixi, ixv
  INTEGER :: itke, ithv
  INTEGER :: itrc_start, itrc_end
  INTEGER :: nmatrix

  INTEGER, ALLOCATABLE :: matrix_idx(:)    !< shape: (nvar_vdiff)
  INTEGER, ALLOCATABLE :: ibtm_var  (:)    !< shape: (nvar_vdiff)
  INTEGER, ALLOCATABLE :: ibtm_mtrx (:)    !< shape: (nmatrix)

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_vdiff_solver'
  CHARACTER(len=*), PARAMETER :: version = '$Id$'

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
  !!   TKE                                     |  1
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
    itke = 6;   ithv = 7
    ih   = 8;   iqv  = 9

    !>KF suggestion
    IF((7 + khydromet) > iqv ) &
    CALL finish( TRIM(thismodule),'matrix for vdiff is not properly defined')

    IF(ktrac > 0)  THEN
      itrc_start = 7 + khydromet +1
      itrc_end   = itrc_start + ktrac - 1
    ELSE
      itrc_start = 7 + khydromet
      itrc_end   = itrc_start
    ENDIF
    !<KF

!    itrc_start = 10
!    itrc_end   = itrc_start + ktrac - 1

    !-------------------------------------------------------------------
    ! # of vertical levels on which the prognostic equations are solved
    !-------------------------------------------------------------------

    ALLOCATE( ibtm_var(nvar_vdiff),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
      & 'Allocation of ibtm_var failed')

    ! momentum, heat, water substances and tracers are solved on        
    ! klev full levels 

    ibtm_var(:)    = klev

    ! TKE and the variance of $\theta_v$ are solved at klev-1 half levels.
    ! The upper and lower boundaries of the atmosphere are excluded.  

    ibtm_var(itke) = klev -1
    ibtm_var(ithv) = klev -1

    !------------------------------------------
    ! Set up matrix indices
    !------------------------------------------

    ALLOCATE( matrix_idx(nvar_vdiff),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
      & 'Allocation of matrix_idx failed')

    matrix_idx(iu)   = 1
    matrix_idx(iv)   = 1  ! u and v share the same exchange coeff.
    matrix_idx(ixl)  = 2
    matrix_idx(ixi)  = 2  ! cloud water and ice share the same exchange coeff.
    matrix_idx(ixv)  = 3
    matrix_idx(itke) = 4
    matrix_idx(ithv) = 5
    matrix_idx(ih)   = 6
    matrix_idx(iqv)  = 6
    IF (ktrac>0) matrix_idx(nvar_vdiff-ktrac+1:nvar_vdiff) = 2

    nmatrix = 6    ! total number of matrices

    !---------------------------------------------------------------------------
    ! # of vertical levels on which elimination of the coefficient will be done
    !---------------------------------------------------------------------------

    ALLOCATE( ibtm_mtrx(nmatrix),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
    & 'Allocation of ibtm_mtrx failed')

    ibtm_mtrx(:)                = klev
    ibtm_mtrx(matrix_idx(itke)) = klev -1
    ibtm_mtrx(matrix_idx(ithv)) = klev -1

  END SUBROUTINE init_vdiff_solver
  !-------------
  !>
  !!
  SUBROUTINE cleanup_vdiff_solver

    INTEGER :: ist

    DEALLOCATE( matrix_idx,ibtm_mtrx,ibtm_var, STAT=ist)
    IF (ist/=SUCCESS) CALL finish('cleanup_vdiff_solver','Deallocation failed')

  END SUBROUTINE cleanup_vdiff_solver
  !-------------
  !>
  !!
  SUBROUTINE matrix_setup_elim( kproma, kbdim, klev, klevm1,  &! in
                              & ksfc_type, itop,              &! in
                              & pcfm, pcfh, pcfh_tile, pcfv,  &! in
                              & pcftke, pcfthv,               &! in
                              & pprfac, prdpm, prdph,         &! in
                              & aa, aa_btm                    )! out
    ! Arguments

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1, ksfc_type
    INTEGER, INTENT(IN) :: itop

    REAL(wp),INTENT(IN) :: pcfm     (kbdim,klev)      !< exchange coeff. for u, v
    REAL(wp),INTENT(IN) :: pcfh     (kbdim,klevm1)    !< exchange coeff. for heat and tracers
    REAL(wp),INTENT(IN) :: pcfh_tile(kbdim,ksfc_type) !< exchange coeff. for heat and qv, at surface
    REAL(wp),INTENT(IN) :: pcfv    (kbdim,klev)      !< exchange coeff. for total water variance
    REAL(wp),INTENT(IN) :: pcftke  (kbdim,klev)      !< exchange coeff. for TKE
    REAL(wp),INTENT(IN) :: pcfthv  (kbdim,klev)      !< exchange coeff. for variance of theta_v
    REAL(wp),INTENT(IN) :: pprfac  (kbdim,klev)      !< prefactor for the exchange coefficients
    REAL(wp),INTENT(IN) :: prdpm   (kbdim,klev)      !< reciprocal of layer thickness, full levels
    REAL(wp),INTENT(IN) :: prdph   (kbdim,klevm1)    !< reciprocal of layer thickness, half levels

    REAL(wp),INTENT(OUT) :: aa    (kbdim,klev,3,nmatrix) !< exchange coeff. matrices
    REAL(wp),INTENT(OUT) :: aa_btm(kbdim,3,ksfc_type)    !< last (the klev-th) row of the coeff. 
                                                         !< matrix for sensible heat and moisture

    ! Local variables

    REAL(wp) :: zkstar (kbdim,itop-1:klev)     !< scaled exchange coeff on half-levels
    REAL(wp) :: zkh    (kbdim,itop-1:klevm1)   !< scaled exchange doeff on full-levels, 
                                               !< for TKE and variance of theta_v
    INTEGER  :: im             !< index of coefficient matrix
    INTEGER  :: jc, jk, jsfc   !< loop indices
    INTEGER  :: jkm1, jmax

    !-----------------------------------------------------------------------
    ! For all prognostic variables: no turbulent flux at the upper boundary
    !-----------------------------------------------------------------------

    zkstar(1:kproma,itop-1) = 0._wp

    !-----------------------------------------------------------------------
    ! For momentum: surface flux is considered
    !-----------------------------------------------------------------------
    im = matrix_idx(iu)    ! also = matrix_idx(iv)
    zkstar(1:kproma,itop:klev) = pprfac(1:kproma,itop:klev)  &
                               &  *pcfm(1:kproma,itop:klev)

    DO jk = itop,klev
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prdpm(jc,jk)  ! -K*_{k-1/2}/dp_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prdpm(jc,jk)  ! -K*_{k+1/2}/dp_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    !---------------------------------------------------------------------
    ! For heat and qv: surface flux is considered.
    ! Note that different surface types are treated separately. The shared
    ! part is levels [itop,klevm1]
    !---------------------------------------------------------------------
    im = matrix_idx(ih)              ! also = matrix_idx(iqv)
    zkstar(1:kproma,itop:klevm1) =  pprfac(1:kproma,itop:klevm1) &
                                 &   *pcfh(1:kproma,itop:klevm1)
    DO jk = itop,klevm1
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prdpm(jc,jk)  ! -K*_{k-1/2}/dp_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prdpm(jc,jk)  ! -K*_{k+1/2}/dp_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    ! Lowest model level: consider surface flux for each surface type

    jk = klev
    DO jsfc = 1,ksfc_type
      DO jc = 1,kproma
        aa_btm(jc,1,jsfc) = -zkstar(jc,jk-1)*prdpm(jc,jk)    ! -K*_{k-1/2}/dp_k
        aa_btm(jc,3,jsfc) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prdpm(jc,jk)
        aa_btm(jc,2,jsfc) = 1._wp - aa_btm(jc,1,jsfc) - aa_btm(jc,3,jsfc)
      ENDDO
    ENDDO

    !----------------------------------------------------------------------
    ! For xl, xi, and other tracers: same coefficients as heat and mositure,
    ! but no turbulent flux at the surface
    !----------------------------------------------------------------------
    im = matrix_idx(ixl)
    zkstar(1:kproma,klev) = 0._wp  ! lower boundary, no turbulent flux

    DO jk = itop,klev
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prdpm(jc,jk)  ! -K*_{k-1/2}/dp_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prdpm(jc,jk)  ! -K*_{k+1/2}/dp_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    !----------------------------------------------------------------------
    ! For total water variance: Surface flux already set to zero 
    ! in subroutine sfc_exchange_coeff.
    !----------------------------------------------------------------------
    im = matrix_idx(ixv)
    zkstar(1:kproma,itop:klev) = pprfac(1:kproma,itop:klev) &
                               &  *pcfv(1:kproma,itop:klev)
    DO jk = itop,klev
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prdpm(jc,jk)
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prdpm(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    !------------------------------------------------------------------------
    ! For TKE: Note that
    ! - Vertical averaging is needed to convert exchange coefficient from
    !   half to full levels, because TKE equation is solved on half levels.
    ! - TKE equation is solved only till array subscript klevm1, which
    !   corresponds to half level (klev - 1/2), i.e., the lowest
    !   interface above surface. Surface value of TKE is (already)
    !   computed in subroutine "sfc_exchange_coeff".
    !------------------------------------------------------------------------
    im = matrix_idx(itke)

    zkstar(1:kproma,itop:klev) =  pprfac(1:kproma,itop:klev) &
                               & *pcftke(1:kproma,itop:klev)
    DO jk = itop,klevm1
      zkh(1:kproma,jk) = 0.5_wp*(zkstar(1:kproma,jk)+zkstar(1:kproma,jk+1))
    ENDDO
    zkh(1:kproma,itop-1) = 0._wp  ! upper boundary, no flux

    DO jk = itop,klevm1
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkh(jc,jk-1)*prdph(jc,jk)
        aa(jc,jk,3,im) = -zkh(jc,jk  )*prdph(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    !------------------------------------------------
    ! For the variance of theta_v (similar to TKE)
    !------------------------------------------------
    im = matrix_idx(ithv)
    zkstar(1:kproma,itop:klev) =  pprfac(1:kproma,itop:klev) &
                               & *pcfthv(1:kproma,itop:klev)
    DO jk = itop,klevm1
      zkh(1:kproma,jk) = 0.5_wp*(zkstar(1:kproma,jk)+zkstar(1:kproma,jk+1))
    ENDDO

    DO jk = itop,klevm1
      DO jc = 1,kproma
        aa(jc,jk,1,im) = -zkh(jc,jk-1)*prdph(jc,jk)
        aa(jc,jk,3,im) = -zkh(jc,jk  )*prdph(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO

    !-------------------------------------------------------------------
    ! Gauss elimination for the coefficient matrices at
    ! - vertical levels [itop,klev-2], for TKE and variance of theta_v;
    ! - vertical levels [itop,klev-1], for all the other variables.
    !-------------------------------------------------------------------

    DO im = 1,nmatrix
      aa(1:kproma,itop,3,im) = aa(1:kproma,itop,3,im)/aa(1:kproma,itop,2,im)

      jmax = ibtm_mtrx(im) - 1
      DO jk = itop+1,jmax
        jkm1 = jk - 1
        aa(1:kproma,jk,2,im) =  aa(1:kproma,jk,2,im)                       &
                             & -aa(1:kproma,jk,1,im)*aa(1:kproma,jkm1,3,im)
        aa(1:kproma,jk,3,im) =  aa(1:kproma,jk,3,im)/aa(1:kproma,jk,2,im)
      ENDDO
    END DO

    ! Translation for developers who prefer to think in terms of 
    ! the Richtmyer-Morthon formula and are familiar with the paper by
    ! Polcher et al (1998): after this elimination, 
    !  aa(:,1:ibtm_mtrx(im)-1,2,:) becomes C  (Eqn. 17),
    !  aa(:,1:ibtm_mtrx(im)-1,3,:) becomes -A (Eqn. 19).

  END SUBROUTINE matrix_setup_elim
  !-------------
  !>
  !!
  !!
  SUBROUTINE rhs_setup( kproma, kbdim, itop, klev, klevm1,   &! in
                      & ksfc_type, ktrac, ptpfac2, pstep_len,&! in
                      & pum1, pvm1, pcptgz, pqm1,            &! in
                      & pxlm1, pxim1, pxvar, pxtm1, pxt_emis,&! in
                      & prdpm, ptkevn, pzthvvar, aa,         &! in
                      & bb, bb_btm                           )! out

    ! Arguments

    INTEGER, INTENT(IN) :: kproma, kbdim, itop, klev, klevm1
    INTEGER, INTENT(IN) :: ksfc_type, ktrac
    REAL(wp),INTENT(IN) :: ptpfac2, pstep_len

    REAL(wp),INTENT(IN) :: pum1     (kbdim,klev)
    REAL(wp),INTENT(IN) :: pvm1     (kbdim,klev)
    REAL(wp),INTENT(IN) :: pcptgz   (kbdim,klev)
    REAL(wp),INTENT(IN) :: pqm1     (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxlm1    (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxim1    (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxvar    (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxtm1    (kbdim,klev,ktrac)
    REAL(wp),INTENT(IN) :: pxt_emis (kbdim,ktrac)
   !REAL(wp),INTENT(IN) :: pxt_emis (kbdim,klev,ktrac) ! backup for later use
    REAL(wp),INTENT(IN) :: ptkevn   (kbdim,klev)
    REAL(wp),INTENT(IN) :: pzthvvar (kbdim,klev)
    REAL(wp),INTENT(IN) :: prdpm    (kbdim,klev)
    REAL(wp),INTENT(IN) :: aa       (kbdim,klev,3,nmatrix)

    REAL(wp),INTENT(OUT) :: bb    (kbdim,klev,nvar_vdiff)
    REAL(wp),INTENT(OUT) :: bb_btm(kbdim,ksfc_type,ih:iqv)

    ! Local variables

    REAL(wp) :: ztmp(kbdim,klev)
    INTEGER  :: jsfc, jt, irhs, im

    !-------------------------------------------------------------------
    ! First handle variables that are defined on full levels
    !-------------------------------------------------------------------
    ! u and v

    bb(1:kproma,itop:klev,iu) = pum1(1:kproma,itop:klev)
    bb(1:kproma,itop:klev,iv) = pvm1(1:kproma,itop:klev)

    ! Hydrometeors and the variance of cloud droplets

    bb(1:kproma,itop:klev,ixl) = pxlm1(1:kproma,itop:klev)
    bb(1:kproma,itop:klev,ixi) = pxim1(1:kproma,itop:klev)
    bb(1:kproma,itop:klev,ixv) = pxvar(1:kproma,itop:klev)

    ! Other tracers

    DO jt = 1,ktrac
       irhs = jt - 1 + itrc_start
       bb(1:kproma,itop:klev,irhs) =  pxtm1(1:kproma,itop:klev,jt)
    ENDDO

    ! Heat and moisture

    bb(1:kproma,itop:klevm1,ih ) = pcptgz(1:kproma,itop:klevm1)
    bb(1:kproma,itop:klevm1,iqv) = pqm1  (1:kproma,itop:klevm1)

    DO jsfc = 1,ksfc_type
       bb_btm(1:kproma,jsfc,ih)  = pcptgz(1:kproma,klev)
       bb_btm(1:kproma,jsfc,iqv) =   pqm1(1:kproma,klev)
    ENDDO

    !-------------------------------------------------------------------
    ! TKE and the variance of theta_v:
    ! These variables are defined at half levels. Array index jk
    ! correspond to half level k+1/2. Thus klev correspond to the
    ! lower boundary. The linear solver only solves till index klevm1.
    !-------------------------------------------------------------------
    im = matrix_idx(itke)
    bb(1:kproma,itop:klevm1,itke) =  ptkevn(1:kproma,itop:klevm1)
    bb(1:kproma,     klevm1,itke) =  bb(1:kproma,klevm1,itke)   &
                                  & -aa(1:kproma,klevm1,3,im)   &
                                  & *ptkevn(1:kproma,klev)

    im = matrix_idx(ithv)
    bb(1:kproma,itop:klevm1,ithv) =  pzthvvar(1:kproma,itop:klevm1)
    bb(1:kproma,     klevm1,ithv) =  bb(1:kproma,klevm1,ithv)     &
                                  & -aa(1:kproma,klevm1,3,im)     &
                                  & *pzthvvar(1:kproma,klev)

    !--------------------------------------------------------------------
    ! Apply the implicitness factor
    !--------------------------------------------------------------------
    !bb     = ptpfac2*bb
    !bb_btm = ptpfac2*bb_btm

     bb(1:kproma,1:klev,  1:itke-1) = ptpfac2*bb(1:kproma,1:klev,  1:itke-1)
     bb(1:kproma,1:klevm1,itke:iqv) = ptpfac2*bb(1:kproma,1:klevm1,itke:iqv)

     IF (ktrac>0) THEN
       bb(1:kproma,1:klev,itrc_start:) = ptpfac2*bb(1:kproma,1:klev,itrc_start:)
     ENDIF

     bb_btm(1:kproma,:,:) = ptpfac2*bb_btm(1:kproma,:,:)

    !--------------------------------------------------------------------
    ! Add tracer emissions
    !--------------------------------------------------------------------
    ! Currently we follow ECHAM in which only the surface emission
    ! is treated "vdiff".

    ztmp(1:kproma,klev) = grav*prdpm(1:kproma,klev)*pstep_len

    DO jt = 1,ktrac
       irhs = jt - 1 + itrc_start
       bb(1:kproma,klev,irhs) =         bb(1:kproma,klev,irhs) &
                              & + pxt_emis(1:kproma,jt)        &
                              &      *ztmp(1:kproma,klev)
    ENDDO

    ! Later we may consider treating emission on all vertical levels
    ! in the same way.
    !
    !ztmp(1:kproma,itop:klev) = grav*prdpm(1:kproma,itop:klev)*pstep_len
    !
    !DO jt = 1,ktrac
    !   irhs = jt - 1 + itrc_start
    !   bb(1:kproma,itop:klev,irhs) =         bb(1:kproma,itop:klev,irhs) &
    !                               & + pxt_emis(1:kproma,itop:klev,jt)   &
    !                               &      *ztmp(1:kproma,itop:klev)
    !ENDDO

  END SUBROUTINE rhs_setup
  !-------------
  !>
  !!
  !! Gauss elimination of the right-hand-side vector 
  !! using coefficients obtained in subroutine "matrix_setup_elim".
  !!
  !! Translation for developers who prefer to think in terms of 
  !! the Richtmyer-Morthon formula and are familiar with the paper by
  !! Polcher et al (1998): after the elimination at the end of 
  !! subroutine matrix_setup_elim, aa(:,1:ibtm_mtrx(im)-1,2,:) 
  !! became the coeff C defined by Eqn. 17 of Polcher et al 1998). 
  !! It is used in this subroutine to convert the variable bb 
  !! into the Richtmyer coeff B (cf Eqn. 19 of Polcher et al 1998).
  !!
  SUBROUTINE rhs_elim( kproma, kbdim, itop, klev, klevm1, &! in 
                     & aa, bb                             )! in, inout

    INTEGER, INTENT(IN)    :: kproma, kbdim, itop, klev, klevm1
    REAL(wp),INTENT(IN)    :: aa(kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT) :: bb(kbdim,klev,nvar_vdiff)

    REAL(wp) :: znum(kbdim), zden(kbdim)
    INTEGER  :: jvar, im, jk, jkm1, jmax

    ! 1. Vertical levels [itop+1,klev-2] for TKE and variance of theta_v;
    !    [itop+1,klev-1] for all the other variables.

    DO jvar = 1,nvar_vdiff
      im = matrix_idx(jvar)  ! Index of coefficient matrix
      bb(1:kproma,itop,jvar) =  bb(1:kproma,itop,jvar)/aa(1:kproma,itop,2,im)

      jmax = ibtm_var(jvar) - 1
      DO jk = itop+1,jmax
         jkm1 = jk - 1
         znum(1:kproma) =  bb(1:kproma,jk  ,jvar)                     &
                        & -bb(1:kproma,jkm1,jvar)*aa(1:kproma,jk,1,im)
         bb(1:kproma,jk,jvar) = znum(1:kproma)/aa(1:kproma,jk,2,im)
      ENDDO
    ENDDO !jvar: variable loop

    ! 2. Vertical level klevm1 for TKE and variance of theta_v.

    DO jvar = 1,nvar_vdiff
      IF (jvar==itke.OR.jvar==ithv) THEN

        im   = matrix_idx(jvar)  ! Index of coefficient matrix
        jk   = ibtm_var(jvar)    ! Bottom level index
        jkm1 = jk - 1
        zden(1:kproma) =  aa(1:kproma,jk,2,im)                      &
                       & -aa(1:kproma,jk,1,im)*aa(1:kproma,jkm1,3,im)
        znum(1:kproma) =  bb(1:kproma,jk,jvar)                      &
                       & -aa(1:kproma,jk,1,im)*bb(1:kproma,jkm1,jvar)
        bb(1:kproma,jk,jvar) = znum(1:kproma)/zden(1:kproma)

      END IF
    ENDDO !jvar: variable loop

    ! Note that for TKE and the variance of theta_v, klev-1 is the lowest
    ! level above surface. This means elimination has already finished
    ! for all levels. Now set boundary condition for the variance of theta_v.

    bb(1:kproma,klev,ithv) = bb(1:kproma,klevm1,ithv)

  END SUBROUTINE rhs_elim
  !-------------
  !>
  !! Solve for heat and water vapour in the lowest model layer. 
  !! First treat each surface type separately, then aggregate 
  !! the solutions. Note that heat and water vapour share the 
  !! same coefficient matrix.
  !!
  SUBROUTINE sfc_solve( kproma, kbdim, klev, klevm1, &! in
                      & ksfc_type, idx_wtr, idx_ice, &! in
                      & lsfc_heat_flux, ptpfac2,     &! in
                      & pfrc, pocu, pocv,            &! in
                      & pcpt_tile, pqsat_tile,       &! in
                      & pcfh_tile, pprfac_sfc,       &! in
                      & aa, aa_btm,                  &! in
                      & bb, bb_btm                   )! inout

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, klevm1
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr, idx_ice 
    LOGICAL, INTENT(IN) :: lsfc_heat_flux
    REAL(wp),INTENT(IN) :: ptpfac2
    REAL(wp),INTENT(IN) :: pfrc      (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pocu      (kbdim)
    REAL(wp),INTENT(IN) :: pocv      (kbdim)
    REAL(wp),INTENT(IN) :: pcpt_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pqsat_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pcfh_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN) :: pprfac_sfc(kbdim)

    REAL(wp),INTENT(IN) :: aa    (kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(IN) :: aa_btm(kbdim,3,ksfc_type)

    REAL(wp),INTENT(INOUT) :: bb    (kbdim,klev,nvar_vdiff)
    REAL(wp),INTENT(INOUT) :: bb_btm(kbdim,ksfc_type,ih:iqv)

    REAL(wp) :: zfrc_oce(kbdim)
    REAL(wp) :: znum(kbdim), zden(kbdim)
    REAL(wp) :: se_sum(kbdim), qv_sum(kbdim), wgt_sum(kbdim), wgt(kbdim)
    INTEGER  :: im, jkm1, jsfc, jk, jvar

    !-----------------------------------------------------------------
    ! Add additional terms to the r.h.s. of momentum, static energy
    ! and moisiture equation by taking into account
    ! - ocean currents (for u, v);
    ! - surface sensible heat flux (for static energy);
    ! - surface moisture flux (for water vapor).
    ! Note that in subroutine rhs_setup the constant ptpfac2 has been 
    ! multiplied to the r.h.s. arrays bb and bb_btm. Thus the 
    ! additional terms here need to be scaled by the same factor.
    !-----------------------------------------------------------------

    IF (idx_wtr.LE.ksfc_type) THEN   ! Open water is considered
      IF (idx_ice.LE.ksfc_type) THEN ! Sea ice is also considered
        zfrc_oce(1:kproma) = pfrc(1:kproma,idx_wtr)+pfrc(1:kproma,idx_ice)
      ELSE ! only open water
        zfrc_oce(1:kproma) = pfrc(1:kproma,idx_wtr)
      ENDIF
      bb(1:kproma,klev,iu) =   bb(1:kproma,klev,iu)                     &
                           & - pocu(1:kproma)*zfrc_oce(1:kproma)*ptpfac2
      bb(1:kproma,klev,iv) =   bb(1:kproma,klev,iv)                     &
                           & - pocv(1:kproma)*zfrc_oce(1:kproma)*ptpfac2
    ENDIF

    ! For moisture and heat, each surface type is treated separately

    DO jsfc = 1,ksfc_type
      bb_btm(1:kproma,jsfc,ih)  =     bb_btm(1:kproma,jsfc,ih) & 
                                & -pcpt_tile(1:kproma,jsfc)    &
                                &    *aa_btm(1:kproma,3,jsfc)  &
                                &    *ptpfac2
    ENDDO

    DO jsfc = 1,ksfc_type
      bb_btm(1:kproma,jsfc,iqv) =      bb_btm(1:kproma,jsfc,iqv) & 
                                & -pqsat_tile(1:kproma,jsfc)     &
                                &     *aa_btm(1:kproma,3,jsfc)   &
                                &     *ptpfac2
    ENDDO

    !--------------------------------------------------------------------------
    ! Bottom level elimination for u, v, as well as tracers excluding moisture
    !--------------------------------------------------------------------------
    DO jvar = 1,nvar_vdiff

      IF ( jvar==itke.OR.jvar==ithv.OR. & ! Computation for TKE and thvvar already done
         & jvar==ih  .OR.jvar==iqv ) THEN ! Heat and moisture are handled in the next block
        CONTINUE
      ELSE

        im   = matrix_idx(jvar)  ! Index of coefficient matrix
        jk   = ibtm_var(jvar)    ! Bottom level index
        jkm1 = jk - 1
        zden(1:kproma) =  aa(1:kproma,jk,2,im)                      &
                       & -aa(1:kproma,jk,1,im)*aa(1:kproma,jkm1,3,im)
        znum(1:kproma) =  bb(1:kproma,jk,jvar)                      &
                       & -aa(1:kproma,jk,1,im)*bb(1:kproma,jkm1,jvar)
        bb(1:kproma,jk,jvar) = znum(1:kproma)/zden(1:kproma)
      END IF 
    ENDDO !jvar: variable loop

    !-------------------------------------------------------------------
    ! Bottom level elimination for heat and moisture
    !-------------------------------------------------------------------

    im = matrix_idx(ih)          ! index of coefficient matrix
    jkm1 = klevm1                ! vertical level index klev-1

     se_sum(1:kproma) = 0._wp    ! sum of weighted solution
     qv_sum(1:kproma) = 0._wp    ! sum of weighted solution
    wgt_sum(1:kproma) = 0._wp    ! sum of weights

    DO jsfc = 1,ksfc_type
       ! Weights for aggregation

       wgt(1:kproma) =  pfrc(1:kproma,jsfc)*pcfh_tile(1:kproma,jsfc)*pprfac_sfc(1:kproma)
       wgt_sum(1:kproma) = wgt_sum(1:kproma) + wgt(1:kproma)

       zden(1:kproma) =  aa_btm(1:kproma,2,jsfc)                     &
                      & -aa_btm(1:kproma,1,jsfc)*aa(1:kproma,jkm1,3,im)

       ! Bottom level elimination for dry static energy

       znum(1:kproma) =  bb_btm(1:kproma,jsfc,ih)                   &
                      & -aa_btm(1:kproma,1,jsfc)*bb(1:kproma,jkm1,ih)

       bb_btm(1:kproma,jsfc,ih) = znum(1:kproma)/zden(1:kproma)
       se_sum(1:kproma) = se_sum(1:kproma) + bb_btm(1:kproma,jsfc,ih)*wgt(1:kproma)

       ! Bottom level elimination for water vapour

       znum(1:kproma) =  bb_btm(1:kproma,jsfc,iqv)                    &
                      & -aa_btm(1:kproma,1,jsfc)*bb(1:kproma,jkm1,iqv)

       bb_btm(1:kproma,jsfc,iqv) = znum(1:kproma)/zden(1:kproma)
       qv_sum(1:kproma) = qv_sum(1:kproma) + bb_btm(1:kproma,jsfc,iqv)*wgt(1:kproma)
    ENDDO

    ! Aggregated solutions at the bottom level

    IF (lsfc_heat_flux) THEN
      bb(1:kproma,klev,ih ) = se_sum(1:kproma)/wgt_sum(1:kproma)
      bb(1:kproma,klev,iqv) = qv_sum(1:kproma)/wgt_sum(1:kproma)
    ELSE
      ! If the surface sensible and heat fluxes are switched off, 
      ! we get the same solution on all surface types.
      ! There is thus no need for aggregation. Simply copy.
      jsfc = 1
      bb(1:kproma,klev,ih ) = bb_btm(1:kproma,jsfc,ih )
      bb(1:kproma,klev,iqv) = bb_btm(1:kproma,jsfc,iqv)
    END IF

  END SUBROUTINE sfc_solve
  !-------------
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
  SUBROUTINE rhs_bksub( kproma, kbdim, itop, klev, aa, bb )

    INTEGER, INTENT(IN)   :: kproma, kbdim, itop, klev
    REAL(wp),INTENT(IN)   :: aa(kbdim,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT):: bb(kbdim,klev,nvar_vdiff)

    INTEGER  :: jvar, im, jk, jkp1

    DO jvar = 1,nvar_vdiff
      im = matrix_idx(jvar)
      DO jk = ibtm_var(jvar)-1,itop,-1
         jkp1 = jk + 1
         bb(1:kproma,jk,jvar) =  bb(1:kproma,jk  ,jvar) &
                              & -bb(1:kproma,jkp1,jvar) &
                              & *aa(1:kproma,jk  ,3,im)
      ENDDO
    ENDDO

  END SUBROUTINE rhs_bksub
  !-------------
  !>
  !!
  SUBROUTINE vdiff_tendencies( kproma, kbdim, itop, klev, klevm1, klevp1,  &! in
                             & ktrac, ksfc_type, idx_lnd, idx_wtr, idx_ice,&! in
                             & pdtime, pstep_len,                          &! in
                             & pum1, pvm1, ptm1, pqm1, pxlm1, pxim1,       &! in
                             & pxtm1, pgeom1, pdelpm1, pcptgz,             &! in
#ifdef __ICON__
                             & ptkem1, pztkevn, pzthvvar, prhoh,           &! in
#else
                             & ptkem1, ptkem0, pztkevn, pzthvvar, prhoh,   &! inout, inout, in
#endif
                             & pqshear, ihpbl, pcfh_tile, pqsat_tile, &! in
                             & pcfm_tile, pfrc, bb,                   &! in
                             & pkedisp, pxvar, pz0m_tile,            &! inout
                             & pute, pvte, ptte, pqte,               &! inout
                             & pxlte, pxite, pxtte,                  &! inout
                             & pevap_ac, plhflx_ac, pshflx_ac,       &! inout
                             & pute_vdf, pvte_vdf, ptte_vdf,         &! out
                             & pqte_vdf, pxlte_vdf, pxite_vdf,       &! out
                             & pxtte_vdf, pxvarprod, pz0m,           &! out
                             & ptke, pthvvar, pthvsig, pvmixtau,     &! out
                             & pqv_mflux_sfc                         )! out

    INTEGER, INTENT(IN) :: kproma, kbdim, itop, klev, klevm1, klevp1, ktrac
    INTEGER, INTENT(IN) :: ksfc_type, idx_lnd, idx_wtr, idx_ice
    REAL(wp),INTENT(IN) :: pstep_len, pdtime

    REAL(wp),INTENT(IN)  :: pum1   (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pvm1   (kbdim,klev)
    REAL(wp),INTENT(IN)  :: ptm1   (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pqm1   (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxlm1  (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxim1  (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxtm1  (kbdim,klev,ktrac)
    REAL(wp),INTENT(IN)  :: pgeom1 (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pdelpm1(kbdim,klev)
    REAL(wp),INTENT(IN)  :: pcptgz (kbdim,klev)
#ifdef __ICON__
    REAL(wp),INTENT(IN)  :: ptkem1 (kbdim,klev)
#else
    REAL(wp),INTENT(INOUT)  :: ptkem1 (kbdim,klev)
    REAL(wp),INTENT(INOUT)  :: ptkem0 (kbdim,klev)
#endif
    REAL(wp),INTENT(IN)  :: pztkevn (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pzthvvar(kbdim,klev)
    REAL(wp),INTENT(IN)  :: prhoh   (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pqshear (kbdim,klev)
    INTEGER, INTENT(IN)  :: ihpbl   (kbdim)
    REAL(wp),INTENT(IN)  :: pcfh_tile (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)  :: pqsat_tile(kbdim,ksfc_type)
    REAL(wp),INTENT(IN)  :: pcfm_tile     (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)  :: pfrc          (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)  :: bb            (kbdim,klev,nvar_vdiff)

    REAL(wp),INTENT(INOUT) :: pkedisp(kbdim) !< temporally and vertically
                                             !< integrated dissipation of
                                             !< kinetic energy
    REAL(wp),INTENT(INOUT) :: pxvar    (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pz0m_tile(kbdim,ksfc_type)

    REAL(wp),INTENT(INOUT) :: pute (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pvte (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: ptte (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pqte (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pxlte(kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pxite(kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pxtte(kbdim,klev,ktrac)

    REAL(wp),INTENT(OUT) :: pute_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pvte_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: ptte_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pqte_vdf (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxlte_vdf(kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxite_vdf(kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxtte_vdf(kbdim,klev,ktrac)
    REAL(wp),INTENT(OUT) :: pxvarprod(kbdim,klev) !< "pvdiffp" in echam

    REAL(wp),INTENT(OUT) :: pz0m    (kbdim)
    REAL(wp),INTENT(INOUT) :: ptke    (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pthvvar (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pthvsig (kbdim)
    REAL(wp),INTENT(OUT) :: pvmixtau(kbdim,klev)
    REAL(wp),INTENT(OUT) :: pqv_mflux_sfc(kbdim)  !< surface mass flux of water vapour
                                                  !< "pqhfla" in echam
    REAL(wp)  :: pevap_tile  (kbdim,ksfc_type)
    REAL(wp)  :: plhflx_tile (kbdim,ksfc_type)
    REAL(wp)  :: pshflx_tile (kbdim,ksfc_type)

    REAL(wp),INTENT(INOUT):: pevap_ac    (kbdim)
    REAL(wp),INTENT(INOUT):: plhflx_ac   (kbdim)
    REAL(wp),INTENT(INOUT):: pshflx_ac   (kbdim)

    REAL(wp) :: ztest, zrdt, zconst
    REAL(wp) :: zunew, zvnew, zqnew, zsnew, ztnew, dqv
    REAL(wp) :: zrhodz, zhexp, zlam, zcons23, z2geomf, zz2geo, zmix, ztkesq
    REAL(wp) :: zvidis(kbdim)
    REAL(wp) :: zdis  (kbdim,klev)
    REAL(wp) :: zqflux(kbdim,klevp1)
    REAL(wp) :: zvarpr(kbdim,klevp1)
    REAL(wp) :: zdqtdt(kbdim,klev)
    INTEGER  :: jk, jl, jt, irhs, jsfc

#ifdef __ICON__
#else
   REAL(wp) ::  zeps
#endif

    zrdt   = 1._wp/pstep_len
    zconst = pdtime/(grav*pstep_len)

    IF (itop>1) THEN
      pute_vdf (1:kproma,1:itop-1)   = 0._wp
      pvte_vdf (1:kproma,1:itop-1)   = 0._wp
      ptte_vdf (1:kproma,1:itop-1)   = 0._wp
      pxlte_vdf(1:kproma,1:itop-1)   = 0._wp
      pxite_vdf(1:kproma,1:itop-1)   = 0._wp
      pxtte_vdf(1:kproma,1:itop-1,:) = 0._wp
    END IF

    !-------------------------------------------------------------------
    ! Compute TKE at the new time step.
    !-------------------------------------------------------------------

    DO jk = itop,klevm1
      ztest = 0._wp
      DO jl = 1,kproma
        ptke(jl,jk) = bb(jl,jk,itke) + tpfac3*pztkevn(jl,jk)
        ztest = ztest+MERGE(1._wp,0._wp,ptke(jl,jk)<0._wp)
      END DO
      IF(ztest.NE.0._wp) THEN
        WRITE(message_text,'(a,I4,2E15.5)') 'level, MIN TKE components = ',&
             & jk, MINVAL(bb(:,jk,itke)),MINVAL(pztkevn(:,jk))
        CALL message('', TRIM(message_text))
        CALL finish('vdiff_tendencies','TKE IS NEGATIVE')
      ENDIF
    END DO
    ptke(1:kproma,klev) = pztkevn(1:kproma,klev)

    

#ifdef __ICON__
#else
    !
    ! TIME FILTER FOR TURBULENT KINETIC ENERGY
    !
     IF(.NOT.lstart) THEN
       zeps=eps
     ELSE
       zeps=0._wp
     END IF
     DO 397 jk = itop,klev
       DO 396 jl = 1,kproma
         ptkem1(jl,jk)=ptkem0(jl,jk)+zeps*(ptkem1(jl,jk)-2._wp*ptkem0(jl,jk)+ptke(jl,jk))
         ptkem0(jl,jk)=ptke(jl,jk)
396     END DO
397  END DO
#endif

    !-------------------------------------------------------------
    ! Variance of virtual potential temperature
    !-------------------------------------------------------------
    DO jk = itop,klev
      DO jl = 1,kproma
        pthvvar(jl,jk) = bb(jl,jk,ithv) + tpfac3*pzthvvar(jl,jk)
        pthvvar(jl,jk) = MAX(tkemin,pthvvar(jl,jk))
      END DO
    END DO

    ! STD DEV OF VIRTUAL POT TEMPERATURE AT STANDARD HALF LEVEL KLEV-1/2
    ! (CORRESPONDING TO HALF LEVEL KLEV-1 FOR TKE AND PTHVVAR)

    pthvsig(1:kproma) = SQRT(pthvvar(1:kproma,klev-1))

    !-------------------------------------------------------------
    ! Tendency of velocity; kinetic energy dissipation
    !-------------------------------------------------------------
    zvidis(1:kproma) = 0._wp   ! initilize the vertical integral

    DO jk = itop,klev
      DO jl = 1,kproma
        pute_vdf(jl,jk) = (bb(jl,jk,iu)-tpfac2*pum1(jl,jk))*zrdt
        pvte_vdf(jl,jk) = (bb(jl,jk,iv)-tpfac2*pvm1(jl,jk))*zrdt

        pute(jl,jk) = pute(jl,jk) + pute_vdf(jl,jk)
        pvte(jl,jk) = pvte(jl,jk) + pvte_vdf(jl,jk)

        zunew = bb(jl,jk,iu) + tpfac3*pum1(jl,jk)
        zvnew = bb(jl,jk,iv) + tpfac3*pvm1(jl,jk)

        zdis(jl,jk) = 0.5_wp*( pum1(jl,jk)**2 - zunew**2 &
                    &         +pvm1(jl,jk)**2 - zvnew**2 )
        zvidis(jl)  = zvidis(jl) + zdis(jl,jk)*pdelpm1(jl,jk)
      END DO
    END DO

    ! Save results

    DO jl=1,kproma
     !udif(jl,jrow) = bb(jl,klev,iu)  !for JSBACH, mo_surface
     !vdif(jl,jrow) = bb(jl,klev,iv)  !for JSBACH, mo_surface
      pkedisp(jl) = pkedisp(jl) + zconst*zvidis(jl) ! BLM dissipation
    END DO

    !-------------------------------------------------------------
    ! Tendency of T and qv, ql, qi; xvar at the new time step
    !-------------------------------------------------------------
    DO jk=itop,klev
      DO jl=1,kproma
        zqnew = bb(jl,jk,iqv) + tpfac3*pqm1(jl,jk)
        pqte_vdf(jl,jk) = (zqnew-pqm1(jl,jk))*zrdt
        pqte(jl,jk) = pqte(jl,jk) + pqte_vdf(jl,jk)

        zsnew = bb(jl,jk,ih) + tpfac3*pcptgz(jl,jk)
        ztnew = (zsnew + zdis(jl,jk) - pgeom1(jl,jk)) &
              & /(cpd*(1._wp+vtmpc2*zqnew))
        ptte_vdf(jl,jk) = (ztnew - ptm1(jl,jk))*zrdt
        ptte(jl,jk) = ptte(jl,jk) + ptte_vdf(jl,jk)

        pxlte_vdf(jl,jk) = (bb(jl,jk,ixl) - tpfac2*pxlm1(jl,jk))*zrdt
        pxite_vdf(jl,jk) = (bb(jl,jk,ixi) - tpfac2*pxim1(jl,jk))*zrdt
        zdqtdt     (jl,jk) =   pqte_vdf (jl,jk) &
                           & + pxlte_vdf(jl,jk) &
                           & + pxite_vdf(jl,jk)

        pxlte(jl,jk) = pxlte(jl,jk) + pxlte_vdf(jl,jk)
        pxite(jl,jk) = pxlte(jl,jk) + pxite_vdf(jl,jk)

        pxvar(jl,jk) = bb(jl,jk,ixv) + tpfac3*pxvar(jl,jk)
      END DO
    END DO

    ! When coupled with JSBACH: Correction of tte for snow melt
    ! ptte_vdf(1:kproma,klev) = ptte_vdf(1:kproma,klev)-pztte_corr(1:kproma)

    !-------------------------------------------------------------
    ! Tendency of tracers
    !-------------------------------------------------------------
!   IF (trlist% anyvdiff /= 0) THEN   ! ECHAM
!     DO 577 jt=1,trlist% ntrac       ! ECHAM
        DO jt = 1,ktrac
          irhs = itrc_start + jt - 1
!         IF (trlist% ti(jt)% nvdiff /= 1) CYCLE  ! ECHAM
          DO jk = itop,klev
            DO jl = 1,kproma
              pxtte_vdf(jl,jk,jt) = (bb(jl,jk,irhs)-tpfac2*pxtm1(jl,jk,jt))*zrdt
              pxtte(jl,jk,jt) = pxtte(jl,jk,jt) + pxtte_vdf(jl,jk,jt)
            ENDDO
          ENDDO
        ENDDO
!577  ENDDO
!     END IF

    !------------------------------------------------------------------------------
    ! Derive the production rate of total water variance.
    ! The production rate reads
    !         -\overline{w'q_t'}\frac{\partial q_t}{\partial}
    ! The first multiplicant is obtained by vertically integrate the equation
    !         d(q_t)/dt = - d(\rho w'q_t')/(rho*dz)
    ! The tendency of total water has already been computed above
    ! and stored in variable "zdqtdt";
    ! Air density and vertical shear of q_v have already been computed in
    ! subroutines "atm_exchange_coeff" and "sfc_exchange_coeff";
    !------------------------------------------------------------------------------
    ! Upper boundary condition: no water flux, no variance production.

    zqflux(1:kproma,itop) = 0._wp
    zvarpr(1:kproma,itop) = 0._wp

    ! Note: for zqflux and zvarpr, index jk corresponds to interface k-1/2;
    ! for zdqtdt, jk corresponds to full level jk; for pqshear and prhoh,
    ! index jk corresponds to interface k+1/2.
    ! Note: The vertical loop should not be parallelized!

    DO jk = itop+1,klev+1
      DO jl=1,kproma
        zrhodz = -pdelpm1(jl,jk-1)/grav
        zqflux(jl,jk) = zrhodz*zdqtdt(jl,jk-1) + zqflux(jl,jk-1)
        zvarpr(jl,jk) = pqshear(jl,jk-1)*zqflux(jl,jk)/prhoh(jl,jk-1)
      ENDDO
    ENDDO

    ! Vertical average of variance production rate from half to full levels

    DO jk=itop,klev
     DO jl=1,kproma
       pxvarprod(jl,jk) = 0.5_wp*(zvarpr(jl,jk)+zvarpr(jl,jk+1))
     ENDDO
    ENDDO

    !---------------------------------------------------------------
    ! Compute the vertical mixing time scale, to be used by "cloud".
    ! Note: there is a similar computation in subroutine "atm_exchange_coeff",
    ! but there the computation is done for half levels, while here
    ! it is done at full levels.
    !---------------------------------------------------------------
    DO jk=itop,klev
       DO jl=1,kproma
          zhexp=EXP(1._wp-pgeom1(jl,jk)/pgeom1(jl,ihpbl(jl)))
          zlam=1._wp+(clam-1._wp)*zhexp
          IF(jk.GE.ihpbl(jl)) THEN
             zcons23=cons25
          ELSE
             zcons23=cons2/zlam
          END IF
          z2geomf=2._wp*pgeom1(jl,jk)
          zz2geo=cons2*z2geomf
          zmix=zz2geo/(1._wp+zcons23*z2geomf)
          IF(jk.EQ.1) THEN
             ztkesq=SQRT(MAX(tkemin,ptkem1(jl,1)))
          ELSE
             ztkesq=SQRT(MAX(tkemin,0.5_wp*(ptkem1(jl,jk-1)  &
                                           +ptkem1(jl,jk))))
          END IF
          pvmixtau(jl,jk) = ztkesq/(zmix*da1)
       ENDDO
    ENDDO

    !---------------------------------------------------------------
    ! Derive surface mass flux of moisture (qv, not q_total).
    ! Formula: area-weighted average of
    ! (air density)*(exchange coef)*(qv_{tavg,klev} - qsat_tile).
    ! Here 
    !   qv_{tavg,klev} = tpfac1*qv_klev(t+dt) + (1-tpfac1)*qv_klev(t)
    !                  = tpfac1*bb_qv
    ! where bb_qv is the solution of the linear system at the lowest
    ! model level (i.e., the full level right above surface).
    !---------------------------------------------------------------
    ! Diagnoise instantaneous moisture flux (= evaporation) on each tile

    DO jsfc = 1,ksfc_type
      IF (jsfc==idx_lnd) CYCLE
      DO jl = 1,kproma
        dqv = tpfac1*bb(jl,klev,iqv) - pqsat_tile(jl,jsfc)
        pevap_tile(jl,jsfc) = prhoh(jl,klev)*pcfh_tile(jl,jsfc)*dqv
      ENDDO
    ENDDO

    IF (idx_lnd<=ksfc_type) CALL finish('vdiff_tendencies','land surface not implemented')

    ! Diagnose latent heat flux (need to distinguish ice and water)

    IF (idx_ice<ksfc_type) THEN 
       plhflx_tile(1:kproma,idx_ice) = als*pevap_tile(1:kproma,idx_ice)
    END IF

    IF (idx_wtr<ksfc_type) THEN 
       plhflx_tile(1:kproma,idx_wtr) = alv*pevap_tile(1:kproma,idx_wtr)
    END IF

    ! Compute grid box mean and time average. 
    ! The instantaneous grid box mean moisture flux will be passed on 
    ! to the cumulus convection scheme.

    pqv_mflux_sfc(1:kproma) = 0._wp   ! initialize mass flux ("pqhfla" in echam)

    DO jsfc = 1,ksfc_type

      pqv_mflux_sfc(1:kproma) =  pqv_mflux_sfc(1:kproma) + pfrc(1:kproma,jsfc) &
                              & *pevap_tile(1:kproma,jsfc)

      pevap_ac(1:kproma)      =  pevap_ac(1:kproma) +  pfrc(1:kproma,jsfc)     &
                              & *pevap_tile(1:kproma,jsfc)*pdtime

      plhflx_ac(1:kproma)     =  plhflx_ac(1:kproma) + pfrc(1:kproma,jsfc)     &
                              & *plhflx_tile(1:kproma,jsfc)*pdtime
    ENDDO

    !----------------------------------------------------------------------------
    ! Update roughness height over open water, then update the grid-box mean
    !----------------------------------------------------------------------------
    IF (idx_wtr<=ksfc_type) THEN  ! water surface exists in the simulation
      DO jl = 1,kproma
        pz0m_tile(jl,idx_wtr) = tpfac1*SQRT( bb(jl,klev,iu)**2+bb(jl,klev,iv)**2 ) &
                              & *pcfm_tile(jl,idx_wtr)*cchar/grav
        pz0m_tile(jl,idx_wtr) = MAX(z0m_min,pz0m_tile(jl,idx_wtr))
      ENDDO
    ENDIF

    ! Compute grid-box mean 

    pz0m(1:kproma) = 0._wp
    DO jsfc = 1,ksfc_type
       pz0m(1:kproma) = pz0m(1:kproma) + pfrc(1:kproma,jsfc)*pz0m_tile(1:kproma,jsfc)
    ENDDO

  END SUBROUTINE vdiff_tendencies
  !-------------



END MODULE mo_vdiff_solver
