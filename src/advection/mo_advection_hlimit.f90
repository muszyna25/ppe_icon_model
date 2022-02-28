!>
!! Flux limiter for horizontal tracer transport
!!
!! This module contains flux limiters for horizontal 
!! tracer transport.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2018-11-29)
!! - module mo_advection_limiter has been splitted into two modules, 
!!   in order to separate horizontal and vertical limiter
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
#define LAXFR_UPFLUX_MACRO(PPp_vn,PPp_psi_a,PPp_psi_b) (0.5_wp*((PPp_vn)*((PPp_psi_a)+(PPp_psi_b))-ABS(PPp_vn)*((PPp_psi_b)-(PPp_psi_a))))
#define LAXFR_UPFLUX_V_MACRO(PPp_w,PPp_psi_a,PPp_psi_b) (0.5_wp*((PPp_w)*((PPp_psi_a)+(PPp_psi_b))+ABS(PPp_w)*((PPp_psi_b)-(PPp_psi_a))))

#ifdef __INTEL_COMPILER
#define USE_LAXFR_MACROS
#define laxfr_upflux LAXFR_UPFLUX_MACRO
#endif
!----------------------------
MODULE mo_advection_hlimit

  USE mo_kind,                ONLY: wp, vp
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_fortran_tools,       ONLY: init
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_sync,                ONLY: SYNC_C1, sync_patch_array, &
    &                               sync_patch_array_mult
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_impl_constants,      ONLY: min_rledge_int, min_rlcell_int
#ifndef USE_LAXFR_MACROS
  USE mo_advection_utils,     ONLY: laxfr_upflux
#endif
#ifdef _OPENACC
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif

  IMPLICIT NONE

  PRIVATE


  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: hflx_limiter_mo
  PUBLIC :: hflx_limiter_pd


#if defined( _OPENACC )
#if defined(__ADVECTION_LIMITER_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
  LOGICAL, PARAMETER ::  acc_validate = .FALSE.   ! ONLY SET TO .TRUE. FOR VALIDATION PHASE
#endif

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Flux limiter for horizontal advection
  !!
  !! Zalesak Flux-Limiter (Flux corrected transport)
  !! The corrected flux is a weighted average of the low order flux and the
  !! given high order flux. The high order flux is used to the greatest extent
  !! possible without introducing overshoots and undershoots.
  !! Note: This limiter is positive definite and almost monotone (but not strictly).
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Schaer, C. and P.K. Smolarkiewicz (1996): A synchronous and iterative 
  !!   flux-correction formalism for coupled transport equations. J. comput. Phys., 
  !!   128, 101-120
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2010-03-10)
  !! Modification by Daniel Reinert, DWD (2010-03-25)
  !! - adapted for MPI parallelization
  !! Modification by Daniel Reinert, DWD (2012-09-20)
  !! - possibility for iterative flux correction
  !! Modification by Daniel Reinert, DWD (2016-09-21)
  !! - remove iterative flux correction, since it does not pay off
  !!
  SUBROUTINE hflx_limiter_mo( ptr_patch, ptr_int, p_dtime, p_cc,            &
    &                         p_rhodz_now, p_rhodz_new, p_mass_flx_e,       &
    &                         p_mflx_tracer_h, slev, elev, opt_beta_fct,    &
    &                         opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(inout) ::  &   !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_int_state), TARGET, INTENT(IN) :: & !< pointer to data structure for
      &  ptr_int                               !< interpolation

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable
      &  p_cc(:,:,:)                 !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< density times cell thickness at timestep n
      &  p_rhodz_now(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< density times cell thickness at timestep n+1
      &  p_rhodz_new(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(in) ::     &    !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)         !< (provided by dynamical core)
                                     !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::     &    !< time step
      &  p_dtime

    REAL(wp), INTENT(INOUT) ::  &    !< calculated horizontal tracer mass flux
      &  p_mflx_tracer_h(:,:,:)      !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN) ::      &    !< vertical start level
      &  slev

    INTEGER, INTENT(IN) ::      &    !< vertical end level
      &  elev

    REAL(wp), INTENT(IN), OPTIONAL ::  & !< factor for multiplicative spreading of range 
      &  opt_beta_fct                    !< of permissible values

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)


    REAL(wp) ::                 &    !< first order tracer mass flux
      &  z_mflx_low(nproma,slev:elev,ptr_patch%nblks_e)

    REAL(wp) ::                 &    !< antidiffusive tracer mass flux (F_H - F_L)
      &  z_anti(nproma,slev:elev,ptr_patch%nblks_e)

    REAL(vp) ::                 &    !< antidiffusive tracer mass flux (F_H - F_L)
      &  z_mflx_anti_1,         &    !< (units kg/kg)
      &  z_mflx_anti_2,         &
      &  z_mflx_anti_3

    REAL(vp) ::                 &    !< sum of incoming antidiffusive tracer mass fluxes
      &  z_mflx_anti_in (nproma,slev:elev,ptr_patch%nblks_c) !< (units kg/kg)

    REAL(vp) ::                 &    !< sum of outgoing antidiffusive tracer mass fluxes
      &  z_mflx_anti_out(nproma,slev:elev,ptr_patch%nblks_c) !< (units kg/kg)

    REAL(vp) ::                 &    !< flux divergence at cell center
      &  z_fluxdiv_c(nproma,slev:elev)

    REAL(wp) ::                 &    !< new tracer field after hor. transport,
      &  z_tracer_new_low(nproma,slev:elev,ptr_patch%nblks_c) 
                                     !< if the low order fluxes are used

    REAL(vp) ::                 &    !< local maximum of current tracer value and low
      &  z_tracer_max(nproma,slev:elev,ptr_patch%nblks_c) !< order update

    REAL(vp) ::                 &    !< local minimum of current tracer value and low
      &  z_tracer_min(nproma,slev:elev,ptr_patch%nblks_c) !< order update

    ! remark: single precision would be sufficient for r_m and r_p, but SP-sync is not yet available
    REAL(wp) ::                 &    !< fraction which must multiply all in/out fluxes 
      &  r_p(nproma,slev:elev,ptr_patch%nblks_c),&   !< of cell jc to guarantee
      &  r_m(nproma,slev:elev,ptr_patch%nblks_c)     !< no overshoot/undershoot

    REAL(wp) :: r_frac !< computed minimum fraction which must multiply
                       !< the flux at the edge

    REAL(vp) :: z_min(nproma,slev:elev), & !< minimum/maximum value in cell and neighboring cells
      &         z_max(nproma,slev:elev) 
    REAL(wp) :: z_signum             !< sign of antidiffusive velocity
    REAL(wp) :: beta_fct             !< factor of allowed over-/undershooting in monotonous limiter
    REAL(wp) :: r_beta_fct           !< ... and its reverse value   

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices of two
      &  iilc, iibc                          !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices of three
      &  iilnc, iibnc                        !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iidx, iblk                          !< of edges

    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: i_rlstart_e, i_rlend_e, i_rlstart_c, i_rlend_c
    INTEGER  :: je, jk, jb, jc         !< index of edge, vert level, block, cell

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: z_mflx_low,z_anti,z_mflx_anti_in,z_mflx_anti_out
!DIR$ ATTRIBUTES ALIGN :64 :: z_fluxdiv_c,z_tracer_new_low,z_tracer_max,z_tracer_min
!DIR$ ATTRIBUTES ALIGN :64 :: r_p,r_m,z_min,z_max
#endif
  !-------------------------------------------------------------------------

    ! Set default values
    i_rlstart = grf_bdywidth_e
    i_rlend   = min_rledge_int - 1
    beta_fct  = 1._wp  ! the namelist default is 1.005, but it is passed to the limiter for the Miura3 scheme only

    ! Check for optional arguments
    IF (PRESENT(opt_rlstart)) i_rlstart = opt_rlstart
    IF (PRESENT(opt_rlend)) i_rlend = opt_rlend
    IF (PRESENT(opt_beta_fct)) beta_fct = opt_beta_fct

    r_beta_fct = 1._wp/beta_fct

    ! number of child domains
    i_nchdom = MAX(1,ptr_patch%n_childdom)

    ! Set pointers to index-arrays

    ! line and block indices of two neighboring cells
    iilc => ptr_patch%edges%cell_idx
    iibc => ptr_patch%edges%cell_blk

    ! line and block indices of edges as seen from cells
    iidx => ptr_patch%cells%edge_idx
    iblk => ptr_patch%cells%edge_blk

    ! pointers to line and block indices of three neighbor cells
    iilnc => ptr_patch%cells%neighbor_idx
    iibnc => ptr_patch%cells%neighbor_blk

!$ACC DATA CREATE( z_mflx_low, z_anti, z_mflx_anti_in, z_mflx_anti_out, r_m, r_p ),          &
!$ACC      CREATE( z_tracer_new_low, z_tracer_max, z_tracer_min, z_min, z_max, z_fluxdiv_c ),&
!$ACC      PCOPYIN( p_cc, p_mass_flx_e, p_rhodz_now, p_rhodz_new ), PCOPY( p_mflx_tracer_h ),&
!$ACC      PRESENT( ptr_patch, ptr_int, iilc, iibc, iilnc, iibnc, iidx, iblk ),              &
!$ACC      IF( i_am_accel_node .AND. acc_on )

!$ACC UPDATE DEVICE( p_cc, p_mass_flx_e, p_mflx_tracer_h ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

    IF (p_test_run) THEN
!$ACC KERNELS PRESENT( r_p, r_m) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      r_p = 0._wp
      r_m = 0._wp
!$ACC END KERNELS
    ENDIF

    !
    ! 1. Calculate low (first) order fluxes using the standard upwind scheme and the
    !    antidiffusive fluxes
    !    (not allowed to call upwind_hflux_up directly, due to circular dependency)

    ! loop through all patch edges (and blocks)

    i_rlstart_e  = 5
    i_rlend_e    = min_rledge_int - 2
    i_startblk   = ptr_patch%edges%start_blk(i_rlstart_e,1)
    i_endblk     = ptr_patch%edges%end_blk(i_rlend_e,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart_e, i_rlend_e)

!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif
          !
          ! compute the first order upwind flux; notice
          ! that only the p_cc*p_vn value at cell edge is computed
          ! multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          !
          z_mflx_low(je,jk,jb) =  &
            &  laxfr_upflux(p_mass_flx_e(je,jk,jb),p_cc(iilc(je,jb,1),jk,iibc(je,jb,1)),p_cc(iilc(je,jb,2),jk,iibc(je,jb,2)))


          ! calculate antidiffusive flux for each edge
          ! only correct for i_rlend_e = min_rledge_int - 1, if p_mflx_tracer_h 
          ! is not synchronized. This is sufficient without iterative flux 
          ! correction which turned out to be overly expensive.
          z_anti(je,jk,jb)     = p_mflx_tracer_h(je,jk,jb) - z_mflx_low(je,jk,jb)


        END DO  ! end loop over edges
      END DO  ! end loop over levels
!$ACC END PARALLEL

    END DO  ! end loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(i_rlstart_c,i_rlend_c,i_startblk,i_endblk)

    i_rlstart_c  = grf_bdywidth_c - 1
    i_rlend_c    = min_rlcell_int - 1
    i_startblk   = ptr_patch%cells%start_blk(i_rlstart_c,1)
    i_endblk     = ptr_patch%cells%end_blk(i_rlend_c,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_fluxdiv_c,z_mflx_anti_1,z_mflx_anti_2,z_mflx_anti_3 ) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,        &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
!DIR$ IVDEP,PREFERVECTOR
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          !
          ! 2. Define "antidiffusive" fluxes A(jc,jk,jb,je) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.

          z_mflx_anti_1 =                                                        &
            &     p_dtime * ptr_int%geofac_div(jc,1,jb) / p_rhodz_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,1),jk,iblk(jc,jb,1))

          z_mflx_anti_2 =                                                        &
            &     p_dtime * ptr_int%geofac_div(jc,2,jb) / p_rhodz_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,2),jk,iblk(jc,jb,2))

          z_mflx_anti_3 =                                                        &
            &     p_dtime * ptr_int%geofac_div(jc,3,jb) / p_rhodz_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,3),jk,iblk(jc,jb,3))

          ! Sum of all incoming antidiffusive fluxes into cell jc
          z_mflx_anti_in(jc,jk,jb) = -1._vp * (MIN(0._vp,z_mflx_anti_1) &
            &                                + MIN(0._vp,z_mflx_anti_2) &
            &                                + MIN(0._vp,z_mflx_anti_3) )

          ! Sum of all outgoing antidiffusive fluxes out of cell jc
          z_mflx_anti_out(jc,jk,jb) = MAX(0._vp,z_mflx_anti_1) &
            &                       + MAX(0._vp,z_mflx_anti_2) &
            &                       + MAX(0._vp,z_mflx_anti_3)

          !  compute also divergence of low order fluxes
          z_fluxdiv_c(jc,jk) =  &       ! This optimization works
            & z_mflx_low(iidx(jc,jb,1),jk,iblk(jc,jb,1)) * ptr_int%geofac_div(jc,1,jb) + &
            & z_mflx_low(iidx(jc,jb,2),jk,iblk(jc,jb,2)) * ptr_int%geofac_div(jc,2,jb) + &
            & z_mflx_low(iidx(jc,jb,3),jk,iblk(jc,jb,3)) * ptr_int%geofac_div(jc,3,jb)
!
! TODO:  the datum  z_mflx_low(iidx(jc,jb,3),jk,iblk(jc,jb,3)) yields differences later in z_tracer_new_low
!        The other entries do not cause a problem. 
!        Status 2015_09_07: problem still there in spite of corrections to mo_nonhydro_gpu_types,
!             both iidx(:,:,3) and iblk(:,:,3) possess problem.
!        Status 2015_09_22: this is related to the COLLAPSE directive mentioned above
!
          z_tracer_new_low(jc,jk,jb) =                        &
            &      ( p_cc(jc,jk,jb) * p_rhodz_now(jc,jk,jb)   &
            &      - p_dtime * z_fluxdiv_c(jc,jk) )           &
            &      / p_rhodz_new(jc,jk,jb)

          ! precalculate local maximum of current tracer value and low order
          ! updated value
          z_tracer_max(jc,jk,jb) = MAX(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))

          ! precalculate local minimum of current tracer value and low order
          ! updated value
          z_tracer_min(jc,jk,jb) = MIN(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))

        ENDDO
      ENDDO
!$ACC END PARALLEL

      IF (ptr_patch%id > 1 .OR. l_limited_area) THEN

        ! Due to the lack of dynamic consistency between mass fluxes and cell mass changes
        ! in the boundary interpolation zone, the low-order advected tracer fields may be
        ! nonsense and therefore need artificial limitation

!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG
        DO jc = i_startidx, i_endidx
          IF (ptr_patch%cells%refin_ctrl(jc,jb) == grf_bdywidth_c-1 .OR. &
              ptr_patch%cells%refin_ctrl(jc,jb) == grf_bdywidth_c) THEN
            !$ACC LOOP VECTOR
            DO jk = slev, elev
              z_tracer_new_low(jc,jk,jb) = MAX(0.9_wp*p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
              z_tracer_new_low(jc,jk,jb) = MIN(1.1_wp*p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
              z_tracer_max(jc,jk,jb) = MAX(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
              z_tracer_min(jc,jk,jb) = MIN(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
            ENDDO
          ENDIF
        ENDDO
!$ACC END PARALLEL

      ENDIF

    ENDDO

!$OMP END DO

    ! Additional initialization of lateral boundary points is needed 
    ! for limited-area mode
    IF ( l_limited_area .AND. ptr_patch%id == 1 ) THEN

      i_startblk   = ptr_patch%cells%start_blk(1,1)
      i_endblk     = ptr_patch%cells%end_blk(grf_bdywidth_c-1,1)

      CALL init(r_m(:,:,i_startblk:i_endblk))
      CALL init(r_p(:,:,i_startblk:i_endblk))

!$OMP BARRIER

    ENDIF

    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.

    i_rlstart_c  = grf_bdywidth_c
    i_rlend_c    = min_rlcell_int
    i_startblk   = ptr_patch%cells%start_blk(i_rlstart_c,1)
    i_endblk     = ptr_patch%cells%end_blk(i_rlend_c,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_max,z_min) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          ! max value of cell and its neighbors
          ! also look back to previous time step
          z_max(jc,jk) = beta_fct * MAX( z_tracer_max(jc,jk,jb),               &
            &                 z_tracer_max(iilnc(jc,jb,1),jk,iibnc(jc,jb,1)),  &
            &                 z_tracer_max(iilnc(jc,jb,2),jk,iibnc(jc,jb,2)),  &
            &                 z_tracer_max(iilnc(jc,jb,3),jk,iibnc(jc,jb,3)) )

          ! min value of cell and its neighbors
          ! also look back to previous time step
          z_min(jc,jk) = r_beta_fct * MIN( z_tracer_min(jc,jk,jb),             &
            &                 z_tracer_min(iilnc(jc,jb,1),jk,iibnc(jc,jb,1)),  &
            &                 z_tracer_min(iilnc(jc,jb,2),jk,iibnc(jc,jb,2)),  &
            &                 z_tracer_min(iilnc(jc,jb,3),jk,iibnc(jc,jb,3)) )
        ENDDO
      ENDDO
!$ACC END PARALLEL

!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of q
          r_m(jc,jk,jb) = (z_tracer_new_low(jc,jk,jb) - z_min(jc,jk))/ &
            &             (z_mflx_anti_out(jc,jk,jb) + dbl_eps)

          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of q
          r_p(jc,jk,jb) = (z_max(jc,jk) - z_tracer_new_low(jc,jk,jb))/ &
            &             (z_mflx_anti_in(jc,jk,jb) + dbl_eps)

        ENDDO
      ENDDO
!$ACC END PARALLEL
    ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Synchronize r_m and r_p and determine i_rlstart/i_rlend
    !
    !$ACC WAIT
    CALL sync_patch_array_mult(SYNC_C1, ptr_patch, 2, r_m, r_p, opt_varname='r_m and r_p')

    !
    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !
    !    - at the end, compute new, limited fluxes which are then passed to 
    !      the main program. Note that p_mflx_tracer_h now denotes the 
    !      LIMITED flux.
    !

    i_startblk = ptr_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_patch%edges%end_blk(i_rlend,i_nchdom)



!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,r_frac,z_signum) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk,                &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      !
      ! compute final limited fluxes
      !
!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif

          z_signum = SIGN(1._wp,z_anti(je,jk,jb))

          ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
          ! but is computationally more efficient
          r_frac = 0.5_wp*( (1._wp+z_signum)*                &
             &     MIN(r_m(iilc(je,jb,1),jk,iibc(je,jb,1)),  &
             &         r_p(iilc(je,jb,2),jk,iibc(je,jb,2)))  &
             &     +  (1._wp-z_signum)*                      &
             &     MIN(r_m(iilc(je,jb,2),jk,iibc(je,jb,2)),  &
             &         r_p(iilc(je,jb,1),jk,iibc(je,jb,1)))  )

          ! Limited flux
          p_mflx_tracer_h(je,jk,jb) = z_mflx_low(je,jk,jb)               &
            &                       + MIN(1._wp,r_frac) * z_anti(je,jk,jb)

        ENDDO
      ENDDO
!$ACC END PARALLEL

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC WAIT
!$ACC UPDATE HOST( p_mflx_tracer_h ) IF (acc_validate .AND. i_am_accel_node .AND. acc_on)
!$ACC END DATA

  END SUBROUTINE hflx_limiter_mo




  !-------------------------------------------------------------------------
  !>
  !! Positive definite flux limiter for horizontal advection
  !!
  !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
  !! Only outward fluxes are re-scaled, in order to maintain positive
  !! definiteness.
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
  !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
  !!   the cubed sphere grid. JCP, in press
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2010-10-06)
  !! - Adaption for hexagonal model by Almut Gassmann, MPI-M (2010-11-18)
  !!
  SUBROUTINE hflx_limiter_pd( ptr_patch, ptr_int, p_dtime, p_cc,        &
    &                         p_rhodz_now, p_mflx_tracer_h, slev, elev, &
    &                         opt_rlstart, opt_rlend )

    TYPE(t_patch), TARGET, INTENT(INOUT) ::  &   !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for
      &  ptr_int                               !< interpolation

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable
      &  p_cc(:,:,:)                 !< dim: (nproma,nlev,nblks_c)
                                     !< [kg kg^-1]

    REAL(wp), INTENT(IN) ::    &    !< density times cell thickness at timestep n
      &  p_rhodz_now(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &    !< time step [s]
      &  p_dtime

    REAL(wp), INTENT(INOUT) ::  &    !< calculated horizontal tracer mass flux
      &  p_mflx_tracer_h(:,:,:)      !< dim: (nproma,nlev,nblks_e)
                                     !< [kg m^-2 s^-1]

    INTEGER, INTENT(IN) ::      &    !< vertical start level
      &  slev

    INTEGER, INTENT(IN) ::      &    !< vertical end level
      &  elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

#if defined(__INTEL_COMPILER) || defined(__SX__) || defined(_OPENACC)
    REAL(wp) :: z_mflx1,  z_mflx2, z_mflx3
#else
    REAL(wp) ::                 &    !< tracer mass flux ( total mass crossing the edge )
      &  z_mflx(nproma,slev:elev,3) !< [kg m^-3]
#endif
    ! remark: single precision would be sufficient for r_m, but SP-sync is not yet available
    REAL(wp) ::                 &    !< fraction which must multiply all outgoing fluxes
      &  r_m(nproma,slev:elev,ptr_patch%nblks_c) !< of cell jc to guarantee
                                                      !< positive definiteness

    REAL(wp) :: z_signum                     !< sign of mass flux
                                             !< >0: out; <0: in
#ifndef __INTEL_COMPILER
    REAL(wp) :: p_m                          !< sum of fluxes out of cell jc
                                             !< [kg m^-3]
#else
    REAL(wp) :: p_m(nproma,slev:elev)
#endif
    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices of two
      &  iilc, iibc                          !< neighbor cells (array)

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iidx, iblk                          !< of edges

    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_rlstart_c, i_rlend_c, i_nchdom
    INTEGER  :: je, jk, jb, jc         !< index of edge, vert level, block, cell

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: p_m,r_m
#endif
  !-------------------------------------------------------------------------

    ! set default values
    i_rlstart = grf_bdywidth_e - 1 ! needed for call from miura_cycl scheme, 
                                   ! otherwise grf_bdywidth_e would be sufficient
    i_rlend   = min_rledge_int - 1


    ! Check for optional arguments
    IF (PRESENT(opt_rlstart)) i_rlstart = opt_rlstart
    IF (PRESENT(opt_rlend)) i_rlend = opt_rlend

    ! number of child domains
    i_nchdom = MAX(1,ptr_patch%n_childdom)

    !
    ! Set pointers to index-arrays
    !
    ! line and block indices of two neighboring cells
    iilc => ptr_patch%edges%cell_idx
    iibc => ptr_patch%edges%cell_blk

    ! line and block indices of edges as seen from cells
    iidx => ptr_patch%cells%edge_idx
    iblk => ptr_patch%cells%edge_blk

!$ACC DATA CREATE( r_m ), PCOPYIN( p_cc, p_rhodz_now ), PCOPY( p_mflx_tracer_h ),  &
!$ACC      PRESENT( ptr_patch, ptr_int, iilc, iibc, iidx, iblk ),                          &
!$ACC      IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc, p_mflx_tracer_h ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

    IF (p_test_run) THEN
!$ACC KERNELS PRESENT( r_m ) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      r_m = 0._wp
!$ACC END KERNELS
    ENDIF

!$OMP PARALLEL PRIVATE(i_rlstart_c,i_rlend_c,i_startblk,i_endblk)
    ! Additional initialization of lateral boundary points is needed for limited-area mode
    IF ( l_limited_area .AND. ptr_patch%id == 1) THEN

      i_startblk   = ptr_patch%cells%start_blk(1,1)
      i_endblk     = ptr_patch%cells%end_blk(grf_bdywidth_c-1,1)

      CALL init(r_m(:,:,i_startblk:i_endblk))
!$OMP BARRIER
    ENDIF

    i_rlstart_c = grf_bdywidth_c
    i_rlend_c   = min_rlcell_int
    i_startblk  = ptr_patch%cells%start_blk(i_rlstart_c,1)
    i_endblk    = ptr_patch%cells%end_blk(i_rlend_c,i_nchdom)

    !
    ! 1. Reformulate all fluxes in terms of the total mass [kg m^-3]
    !    that crosses each of the CV-edges and store them in a cell-based structure.
    !
    !    z_mflx > 0: outward
    !    z_mflx < 0: inward
    !

#if defined (__INTEL_COMPILER) || defined (__SX__)
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,p_m, &
!$OMP            z_mflx1,z_mflx2,z_mflx3) ICON_OMP_DEFAULT_SCHEDULE
#else
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,p_m,z_mflx) ICON_OMP_DEFAULT_SCHEDULE
#endif
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,        &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

#if defined (__SX__) || defined ( _OPENACC )
          z_mflx1 = ptr_int%geofac_div(jc,1,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,1),jk,iblk(jc,jb,1))
          z_mflx2 = ptr_int%geofac_div(jc,2,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,2),jk,iblk(jc,jb,2))
          z_mflx3 = ptr_int%geofac_div(jc,3,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,3),jk,iblk(jc,jb,3))

          ! Sum of all outgoing fluxes out of cell jc
          p_m =  MAX(0._wp,z_mflx1) + MAX(0._wp,z_mflx2) + MAX(0._wp,z_mflx3)

          ! fraction which must multiply all fluxes out of cell jc to guarantee no undershoot
          ! Nominator: maximum allowable decrease of \rho q
          r_m(jc,jk,jb) = MIN(1._wp, (p_cc(jc,jk,jb)*p_rhodz_now(jc,jk,jb)) / (p_m + dbl_eps) )

#elif !defined (__INTEL_COMPILER)

          z_mflx(jc,jk,1) = ptr_int%geofac_div(jc,1,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,1),jk,iblk(jc,jb,1))

          z_mflx(jc,jk,2) = ptr_int%geofac_div(jc,2,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,2),jk,iblk(jc,jb,2))
  
          z_mflx(jc,jk,3) = ptr_int%geofac_div(jc,3,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,3),jk,iblk(jc,jb,3))
  
        ENDDO
      ENDDO
!$ACC END PARALLEL

      !
      ! 2. Compute total outward mass
      !
!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx

          ! Sum of all outgoing fluxes out of cell jc
          p_m =  MAX(0._wp,z_mflx(jc,jk,1))  &
            &  + MAX(0._wp,z_mflx(jc,jk,2))  &
            &  + MAX(0._wp,z_mflx(jc,jk,3))

          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of \rho q
          r_m(jc,jk,jb) = MIN(1._wp, (p_cc(jc,jk,jb)*p_rhodz_now(jc,jk,jb)) &
            &                        /(p_m + dbl_eps) )

#else
          z_mflx1 = ptr_int%geofac_div(jc,1,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,1),jk,iblk(jc,jb,1))
          z_mflx2 = ptr_int%geofac_div(jc,2,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,2),jk,iblk(jc,jb,2))
          z_mflx3 = ptr_int%geofac_div(jc,3,jb) * p_dtime &
            &                * p_mflx_tracer_h(iidx(jc,jb,3),jk,iblk(jc,jb,3))
          ! Sum of all outgoing fluxes out of cell jc
          p_m(jc,jk) = MAX(0._wp,z_mflx1) + &
            &          MAX(0._wp,z_mflx2) + &
            &          MAX(0._wp,z_mflx3)
        ENDDO
      ENDDO
      DO jk = slev, elev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          ! fraction which must multiply all fluxes out of cell jc to guarantee
          ! no
          ! undershoot
          ! Nominator: maximum allowable decrease of \rho q
          r_m(jc,jk,jb) = MIN(1._wp, (p_cc(jc,jk,jb)*p_rhodz_now(jc,jk,jb)) &
            &                        /(p_m(jc,jk) + dbl_eps) )

          
#endif
        ENDDO
      ENDDO
!$ACC END PARALLEL
    ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! synchronize r_m
    !
    !$ACC WAIT
    IF(SIZE(r_m)/=0) CALL sync_patch_array(SYNC_C1,ptr_patch,r_m,opt_varname='r_m')

    !
    ! 3. Limit outward fluxes
    !    The inward ones remain untouched.
    !
    i_startblk = ptr_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,z_signum) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk,    &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
#ifdef __LOOP_EXCHANGE
      !$ACC LOOP GANG
      DO je = i_startidx, i_endidx
        ! this is potentially needed for calls from miura_cycl
        IF (ptr_patch%edges%refin_ctrl(je,jb) == grf_bdywidth_e-2) CYCLE
        !$ACC LOOP VECTOR PRIVATE( z_signum )
        DO jk = slev, elev
#else
!$NEC outerloop_unroll(4)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
          IF (ptr_patch%edges%refin_ctrl(je,jb) == grf_bdywidth_e-2) CYCLE
#endif

          ! p_mflx_tracer_h > 0: flux directed from cell 1 -> 2
          ! p_mflx_tracer_h < 0: flux directed from cell 2 -> 1
          z_signum = SIGN(1._wp,p_mflx_tracer_h(je,jk,jb))

          p_mflx_tracer_h(je,jk,jb) = p_mflx_tracer_h(je,jk,jb) * 0.5_wp  &
            & *( (1._wp + z_signum) * r_m(iilc(je,jb,1),jk,iibc(je,jb,1)) &
            &   +(1._wp - z_signum) * r_m(iilc(je,jb,2),jk,iibc(je,jb,2)) )
  
        ENDDO
      ENDDO
!$ACC END PARALLEL
    ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC WAIT
!$ACC UPDATE HOST( p_mflx_tracer_h ) IF (acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE hflx_limiter_pd

END MODULE mo_advection_hlimit

