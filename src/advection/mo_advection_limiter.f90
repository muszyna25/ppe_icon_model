!>
!! Slope and flux limiter for transport scheme
!!
!! This module contains all available slope and flux limiters for both
!! the horizontal and vertical transport.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2010-03-04)
!! Modification by Daniel Reinert, DWD (2010-03-04)
!! - transferred all limiter to this new module
!! - implementation of flux limiter for FFSL-scheme (Miura)
!! Modification by Daniel Reinert, DWD (2010-10-06)
!! - implemented positive definite FCT-limiter for FFSL-scheme
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
MODULE mo_advection_limiter

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_physical_constants,  ONLY: re
  USE mo_model_domain,        ONLY: t_patch
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_sync,                ONLY: SYNC_C, SYNC_C1, sync_patch_array,     &
    &                               sync_patch_array_mult
  USE mo_parallel_config, ONLY: nproma, p_test_run, n_ghost_rows
  USE mo_math_utilities,      ONLY: cc2gc,t_cartesian_coordinates,          &
    &                               t_geographical_coordinates, gvec2cvec,  &
    &                               gnomonic_proj
  USE mo_interpolation,       ONLY: t_int_state
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_impl_constants,      ONLY: min_rledge_int, min_rlcell_int, min_rlcell, &
    &                               min_rledge
  USE mo_math_operators,      ONLY: div
  USE mo_advection_utils,     ONLY: laxfr_upflux, ptr_delp_mc_now,        &
    &                               ptr_delp_mc_new

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: hflx_limiter_mo, hflx_limiter_sm, h_miura_slimiter_mo,        &
    &       h_miura_slimiter_sm, shift_gauss_points, v_muscl_slimiter_mo, &
    &       v_muscl_slimiter_sm, v_ppm_slimiter_mo, v_ppm_slimiter_sm,    &
    &       vflx_limiter_pd, vflx_limiter_pd_ha

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
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2010-03-10)
  !! Modification by Daniel Reinert, DWD (2010-03-25)
  !! - adapted for MPI parallelization
  !!
  SUBROUTINE hflx_limiter_mo( ptr_patch, ptr_int, p_dtime, p_cc, p_mass_flx_e, &
    &                         p_mflx_tracer_h, opt_rlstart, opt_rlend,         &
    &                         opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_int_state), TARGET, INTENT(IN) :: & !< pointer to data structure for
      &  ptr_int                               !< interpolation

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable
      &  p_cc(:,:,:)                 !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(in) ::     &    !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)         !< (provided by dynamical core)
                                     !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::     &    !< time step
      &  p_dtime

    REAL(wp), INTENT(INOUT) ::  &    !< calculated horizontal tracer mass flux
      &  p_mflx_tracer_h(:,:,:)      !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    REAL(wp) ::                 &    !< first order tracer mass flux
      &  z_mflx_low(nproma,ptr_patch%nlev,ptr_patch%nblks_e)

    REAL(wp) ::                 &    !< antidiffusive tracer mass flux (F_H - F_L)
      &  z_anti(nproma,ptr_patch%nlev,ptr_patch%nblks_e)

    REAL(wp) ::                 &    !< antidiffusive tracer mass flux (F_H - F_L)
      &  z_mflx_anti(nproma,ptr_patch%nlev,ptr_patch%nblks_c,3) !< (units kg/kg)

    REAL(wp) ::                 &    !< flux divergence at cell center
      &  z_fluxdiv_c(nproma,ptr_patch%nlev,ptr_patch%nblks_c)

    REAL(wp) ::                 &    !< new tracer field after hor. transport,
      &  z_tracer_new_low(nproma,ptr_patch%nlev,ptr_patch%nblks_c) 
                                     !< if the low order fluxes are used

    REAL(wp) ::                 &    !< local maximum of current tracer value and low
      &  z_tracer_max(nproma,ptr_patch%nlev,ptr_patch%nblks_c) !< order update

    REAL(wp) ::                 &    !< local minimum of current tracer value and low
      &  z_tracer_min(nproma,ptr_patch%nlev,ptr_patch%nblks_c) !< order update

    REAL(wp) ::                 &    !< fraction which must multiply all in/out fluxes 
      &  r_p(nproma,ptr_patch%nlev,ptr_patch%nblks_c),&   !< of cell jc to guarantee
      &  r_m(nproma,ptr_patch%nlev,ptr_patch%nblks_c)     !< no overshoot/undershoot

    REAL(wp) :: r_frac(nproma) !< computed minimum fraction which must multiply
                               !< the flux at the edge

    REAL(wp) :: z_min(nproma), z_max(nproma) !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_signum(nproma)             !< sign of antidiffusive velocity
    REAL(wp) :: p_p(nproma), p_m(nproma)     !< sum of antidiffusive fluxes into and out of cell jc

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices of two
      &  iilc, iibc                          !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices of three
      &  iilnc, iibnc                        !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iidx, iblk                          !< of edges

    INTEGER  :: slev, elev             !< vertical start and end level
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: i_rlstart_e, i_rlend_e, i_rlstart_c, i_rlend_c
    INTEGER  :: je, jk, jb, jc         !< index of edge, vert level, block, cell

  !-------------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = ptr_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_e
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,ptr_patch%n_childdom)


    IF (p_test_run) THEN
      r_p = 0._wp
      r_m = 0._wp
    ENDIF


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
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart_e, i_rlend_e)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=5
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
            &  laxfr_upflux( p_mass_flx_e(je,jk,jb), p_cc(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            &                             p_cc(iilc(je,jb,2),jk,iibc(je,jb,2)) )

          ! calculate antidiffusive flux for each edge
          z_anti(je,jk,jb)     = p_mflx_tracer_h(je,jk,jb) - z_mflx_low(je,jk,jb)

        END DO  ! end loop over edges

      END DO  ! end loop over levels

    END DO  ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL


    !
    ! 2. Compute the updated low order solution z_tracer_new_low
    !
    !  compute divergence of low order fluxes
    CALL div( z_mflx_low, ptr_patch, ptr_int,              &! in
      &       z_fluxdiv_c(:,:,:),                          &! inout
      &       opt_rlstart=3, opt_rlend=min_rlcell_int-1,   &! in
      &       opt_slev=slev, opt_elev=elev                 )! in


!$OMP PARALLEL PRIVATE(i_rlstart_c,i_rlend_c,i_startblk,i_endblk)
    i_rlstart_c  = grf_bdywidth_c - 1
    i_rlend_c    = min_rlcell_int - 1
    i_startblk   = ptr_patch%cells%start_blk(i_rlstart_c,1)
    i_endblk     = ptr_patch%cells%end_blk(i_rlend_c,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,        &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=5
      DO jk = slev, elev

        DO jc = i_startidx, i_endidx
#endif

          z_tracer_new_low(jc,jk,jb) =                           &
            &      ( p_cc(jc,jk,jb) * ptr_delp_mc_now(jc,jk,jb)  &
            &      - p_dtime * z_fluxdiv_c(jc,jk,jb) )           &
            &      / ptr_delp_mc_new(jc,jk,jb)

    !
    ! 3. Define "antidiffusive" fluxes A(jc,jk,jb,je) for each cell. It is the difference
    !    between the high order fluxes (given by the FFSL-scheme) and the low order
    !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
    !    - positive for outgoing fluxes
    !    - negative for incoming fluxes
    !    this sign convention is related to the definition of the divergence operator.

          z_mflx_anti(jc,jk,jb,1) =                                                  &
            &     p_dtime * ptr_int%geofac_div(jc,1,jb) / ptr_delp_mc_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,1),jk,iblk(jc,jb,1))

          z_mflx_anti(jc,jk,jb,2) =                                                  &
            &     p_dtime * ptr_int%geofac_div(jc,2,jb) / ptr_delp_mc_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,2),jk,iblk(jc,jb,2))

          z_mflx_anti(jc,jk,jb,3) =                                                  &
            &     p_dtime * ptr_int%geofac_div(jc,3,jb) / ptr_delp_mc_new(jc,jk,jb)  &
            &   * z_anti(iidx(jc,jb,3),jk,iblk(jc,jb,3))

          ! precalculate local maximum of current tracer value and low order
          ! updated value
          z_tracer_max(jc,jk,jb) = MAX(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))

          ! precalculate local minimum of current tracer value and low order
          ! updated value
          z_tracer_min(jc,jk,jb) = MIN(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO

    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.

    i_rlstart_c  = grf_bdywidth_c
    i_rlend_c    = min_rlcell_int
    i_startblk   = ptr_patch%cells%start_blk(i_rlstart_c,1)
    i_endblk     = ptr_patch%cells%end_blk(i_rlend_c,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_max,p_p,z_min,p_m)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=5
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          ! Sum of all incoming antidiffusive fluxes into cell jc
          p_p(jc) =  -1._wp * (MIN(0._wp,z_mflx_anti(jc,jk,jb,1))   &
                             + MIN(0._wp,z_mflx_anti(jc,jk,jb,2))   &
                             + MIN(0._wp,z_mflx_anti(jc,jk,jb,3)) )

          ! Sum of all outgoing antidiffusive fluxes out of cell jc
          p_m(jc) =  MAX(0._wp,z_mflx_anti(jc,jk,jb,1))  &
                &  + MAX(0._wp,z_mflx_anti(jc,jk,jb,2))  &
                &  + MAX(0._wp,z_mflx_anti(jc,jk,jb,3))

          ! max value of cell and its neighbors
          ! also look back to previous time step
          z_max(jc) = MAX( z_tracer_max(jc,jk,jb),                          &
            &              z_tracer_max(iilnc(jc,jb,1),jk,iibnc(jc,jb,1)),  &
            &              z_tracer_max(iilnc(jc,jb,2),jk,iibnc(jc,jb,2)),  &
            &              z_tracer_max(iilnc(jc,jb,3),jk,iibnc(jc,jb,3)) )

          ! min value of cell and its neighbors
          ! also look back to previous time step
          z_min(jc) = MIN( z_tracer_min(jc,jk,jb),                          &
            &              z_tracer_min(iilnc(jc,jb,1),jk,iibnc(jc,jb,1)),  &
            &              z_tracer_min(iilnc(jc,jb,2),jk,iibnc(jc,jb,2)),  &
            &              z_tracer_min(iilnc(jc,jb,3),jk,iibnc(jc,jb,3)) )

          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of q
          r_m(jc,jk,jb) = (z_tracer_new_low(jc,jk,jb) - z_min(jc))/(p_m(jc) + dbl_eps)

          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of q
          r_p(jc,jk,jb) = (z_max(jc) - z_tracer_new_low(jc,jk,jb))/(p_p(jc) + dbl_eps)

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(SYNC_C1, ptr_patch, 2, r_m, r_p)

    !
    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
    !    program. Note that p_mflx_tracer_h now denotes the LIMITED flux.
    !

    i_startblk = ptr_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = ptr_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,r_frac,z_signum)
    DO jb = i_startblk, i_endblk


      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk,    &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif

          z_signum(je) = SIGN(1._wp,z_anti(je,jk,jb))

          ! This does the same as an IF (z_signum(je) > 0) THEN ... ELSE ... ENDIF,
          ! but is computationally more efficient
          r_frac(je) = 0.5_wp*( (1._wp+z_signum(je))*              &
             &         MIN(r_m(iilc(je,jb,1),jk,iibc(je,jb,1)),    &
             &             r_p(iilc(je,jb,2),jk,iibc(je,jb,2)))    &
             &         -  (z_signum(je)-1._wp)*                    &
             &         MIN(r_m(iilc(je,jb,2),jk,iibc(je,jb,2)),    &
             &             r_p(iilc(je,jb,1),jk,iibc(je,jb,1)))    )

          ! Limited flux
          p_mflx_tracer_h(je,jk,jb) = z_mflx_low(je,jk,jb)                  &
            &                       + MIN(1._wp,r_frac(je)) * z_anti(je,jk,jb)

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

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
  SUBROUTINE hflx_limiter_sm( ptr_patch, ptr_int, p_dtime, p_cc,        &
    &                         p_mflx_tracer_h, opt_rlstart, opt_rlend,  &
    &                         opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for
      &  ptr_int                               !< interpolation

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable
      &  p_cc(:,:,:)                 !< dim: (nproma,nlev,nblks_c)
                                     !< [kg kg^-1]

    REAL(wp), INTENT(IN) ::     &    !< time step [s]
      &  p_dtime

    REAL(wp), INTENT(INOUT) ::  &    !< calculated horizontal tracer mass flux
      &  p_mflx_tracer_h(:,:,:)      !< dim: (nproma,nlev,nblks_e)
                                     !< [kg m^-2 s^-1]

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    REAL(wp) ::                 &    !< tracer mass flux ( total mass crossing the edge )
      &  z_mflx(nproma,ptr_patch%nlev,ptr_patch%nblks_c,ptr_patch%cell_type) !< [kg m^-3]

    REAL(wp) ::                 &    !< fraction which must multiply all outgoing fluxes
      &  r_m(nproma,ptr_patch%nlev,ptr_patch%nblks_c) !< of cell jc to guarantee
                                                      !< positive definiteness

    REAL(wp) :: z_signum(nproma)             !< sign of mass flux
                                             !< >0: out; <0: in
    REAL(wp) :: p_m(nproma)                  !< sum of fluxes out of cell jc
                                             !< [kg m^-3]

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices of two
      &  iilc, iibc                          !< neighbor cells (array)

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iidx, iblk                          !< of edges

    INTEGER  :: slev, elev                   !< vertical start and end level
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_rlstart_c, i_rlend_c, i_nchdom
    INTEGER  :: je, jk, jb, jc         !< index of edge, vert level, block, cell

  !-------------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = ptr_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_e
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,ptr_patch%n_childdom)


    IF (p_test_run) THEN
      r_m = 0._wp
    ENDIF

    !
    ! Set pointers to index-arrays
    !
    ! line and block indices of two neighboring cells
    iilc => ptr_patch%edges%cell_idx
    iibc => ptr_patch%edges%cell_blk

    ! line and block indices of edges as seen from cells
    iidx => ptr_patch%cells%edge_idx
    iblk => ptr_patch%cells%edge_blk


!$OMP PARALLEL PRIVATE(i_rlstart_c,i_rlend_c,i_startblk,i_endblk)

    SELECT CASE (ptr_patch%cell_type)

    CASE(3)

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

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,p_m)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,        &
                           i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = slev, elev
#else
!CDIR UNROLL=5
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
#endif

            z_mflx(jc,jk,jb,1) = ptr_int%geofac_div(jc,1,jb) * p_dtime &
              &                * p_mflx_tracer_h(iidx(jc,jb,1),jk,iblk(jc,jb,1))

            z_mflx(jc,jk,jb,2) = ptr_int%geofac_div(jc,2,jb) * p_dtime &
              &                * p_mflx_tracer_h(iidx(jc,jb,2),jk,iblk(jc,jb,2))
  
            z_mflx(jc,jk,jb,3) = ptr_int%geofac_div(jc,3,jb) * p_dtime &
              &                * p_mflx_tracer_h(iidx(jc,jb,3),jk,iblk(jc,jb,3))
  
          ENDDO
        ENDDO

        !
        ! 2. Compute total outward mass
        !
!CDIR UNROLL=5
        DO jk = slev, elev

          DO jc = i_startidx, i_endidx

            ! Sum of all outgoing fluxes out of cell jc
            p_m(jc) =  MAX(0._wp,z_mflx(jc,jk,jb,1))  &
                  &  + MAX(0._wp,z_mflx(jc,jk,jb,2))  &
                  &  + MAX(0._wp,z_mflx(jc,jk,jb,3))

            ! fraction which must multiply all fluxes out of cell jc to guarantee no
            ! undershoot
            ! Nominator: maximum allowable decrease of \rho q
            r_m(jc,jk,jb) = MIN(1._wp, (p_cc(jc,jk,jb)*ptr_delp_mc_now(jc,jk,jb)) &
              &                        /(p_m(jc) + dbl_eps) )

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

    CASE(6)

      i_rlstart_c = grf_bdywidth_c 
      i_rlend_c   = min_rlcell
      i_startblk  = ptr_patch%cells%start_blk(i_rlstart_c,1)
      i_endblk    = ptr_patch%cells%end_blk(i_rlend_c,i_nchdom)

      !
      ! 1. Reformulate all fluxes in terms of the total mass [kg m^-3]
      !    that crosses each of the CV-edges and store them in a cell-based structure.
      !
      !    z_mflx > 0: outward
      !    z_mflx < 0: inward
      !

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,p_m)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,        &
                           i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = slev, elev
#else
!CDIR UNROLL=5
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx
#endif

            z_mflx(jc,jk,jb,1) = ptr_int%geofac_div(jc,1,jb) * p_dtime &
              &                * p_mflx_tracer_h(iidx(jc,jb,1),jk,iblk(jc,jb,1))

            z_mflx(jc,jk,jb,2) = ptr_int%geofac_div(jc,2,jb) * p_dtime &
              &                * p_mflx_tracer_h(iidx(jc,jb,2),jk,iblk(jc,jb,2))
  
            z_mflx(jc,jk,jb,3) = ptr_int%geofac_div(jc,3,jb) * p_dtime &
              &                * p_mflx_tracer_h(iidx(jc,jb,3),jk,iblk(jc,jb,3))

            z_mflx(jc,jk,jb,4) = ptr_int%geofac_div(jc,4,jb) * p_dtime &
              &                * p_mflx_tracer_h(iidx(jc,jb,4),jk,iblk(jc,jb,4))

            z_mflx(jc,jk,jb,5) = ptr_int%geofac_div(jc,5,jb) * p_dtime &
              &                * p_mflx_tracer_h(iidx(jc,jb,5),jk,iblk(jc,jb,5))

            ! For pentagons, geofac_div=0 and iidx is the 5th edge => no problem
            z_mflx(jc,jk,jb,6) = ptr_int%geofac_div(jc,6,jb) * p_dtime &
              &                * p_mflx_tracer_h(iidx(jc,jb,6),jk,iblk(jc,jb,6))
  
          ENDDO
        ENDDO

        !
        ! 2. Compute total outward mass
        !
!CDIR UNROLL=5
        DO jk = slev, elev

          DO jc = i_startidx, i_endidx

            ! Sum of all outgoing fluxes out of cell jc
            p_m(jc) =  MAX(0._wp,z_mflx(jc,jk,jb,1))  &
                  &  + MAX(0._wp,z_mflx(jc,jk,jb,2))  &
                  &  + MAX(0._wp,z_mflx(jc,jk,jb,3))  &
                  &  + MAX(0._wp,z_mflx(jc,jk,jb,4))  &
                  &  + MAX(0._wp,z_mflx(jc,jk,jb,5))  &
                  &  + MAX(0._wp,z_mflx(jc,jk,jb,6))

            ! fraction which must multiply all fluxes out of cell jc to guarantee no
            ! undershoot
            ! Nominator: maximum allowable decrease of \rho q
            r_m(jc,jk,jb) = MIN(1._wp, (p_cc(jc,jk,jb)*ptr_delp_mc_now(jc,jk,jb)) &
              &                        /(p_m(jc) + dbl_eps) )

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

    END SELECT ! triangles or hexagons
!$OMP END PARALLEL

    ! synchronize r_m
    CALL sync_patch_array(SYNC_C1,ptr_patch,r_m)


    !
    ! 3. Limit outward fluxes
    !    The inward ones remain untouched.
    !
    SELECT CASE (ptr_patch%cell_type)
    CASE(3)

      i_startblk = ptr_patch%edges%start_blk(i_rlstart,1)
      i_endblk   = ptr_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,z_signum)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk,    &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = slev, elev
#else
!CDIR UNROLL=5
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
#endif

            ! p_mflx_tracer_h > 0: flux directed from cell 1 -> 2
            ! p_mflx_tracer_h < 0: flux directed from cell 2 -> 1
            z_signum(je) = SIGN(1._wp,p_mflx_tracer_h(je,jk,jb))

            p_mflx_tracer_h(je,jk,jb) = p_mflx_tracer_h(je,jk,jb) * 0.5_wp      &
              & *( (1._wp + z_signum(je)) * r_m(iilc(je,jb,1),jk,iibc(je,jb,1)) &
              &   +(1._wp - z_signum(je)) * r_m(iilc(je,jb,2),jk,iibc(je,jb,2)) )
  
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    CASE(6)
      i_startblk = ptr_patch%edges%start_blk(i_rlstart,1)
      i_endblk   = ptr_patch%edges%end_blk(min_rledge,i_nchdom)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,z_signum)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = slev, elev
#else
!CDIR UNROLL=5
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
#endif

            ! p_mflx_tracer_h > 0: flux directed from cell 1 -> 2
            ! p_mflx_tracer_h < 0: flux directed from cell 2 -> 1
            z_signum(je) = SIGN(1._wp,p_mflx_tracer_h(je,jk,jb)) &
              & *ptr_patch%edges%system_orientation(je,jb)

            p_mflx_tracer_h(je,jk,jb) = p_mflx_tracer_h(je,jk,jb) * 0.5_wp      &
              & *( (1._wp + z_signum(je)) * r_m(iilc(je,jb,1),jk,iibc(je,jb,1)) &
              &   +(1._wp - z_signum(je)) * r_m(iilc(je,jb,2),jk,iibc(je,jb,2)) )
  
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    END SELECT

  END SUBROUTINE hflx_limiter_sm


  !-------------------------------------------------------------------------
  !>
  !! Positive definite flux limiter for vertical advection
  !!
  !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
  !! for the nonhydrostatic core. Only outward fluxes are re-scaled, in 
  !! order to maintain positive definiteness.
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
  !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
  !!   the cubed sphere grid.  J. Comput. Phys., 230, 1215-1237
  !! - Smolarkiewicz, P. K., 1989: Comment on "A positive definite advection 
  !!   scheme obtained by nonlinear renormalization of the advective fluxes.", 
  !!   Mon. Wea. Rev., 117, 2626-2632
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2011-01-07)
  !!
  SUBROUTINE vflx_limiter_pd( ptr_patch, p_dtime, p_cc, p_mflx_tracer_v, &
    &                         opt_rlstart, opt_rlend, opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  ptr_patch

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable at time (n)
      &  p_cc(:,:,:)                 !< dim: (nproma,nlev,nblks_c)
                                     !< [kg kg^-1]

    REAL(wp), INTENT(IN) ::     &    !< time step [s]
      &  p_dtime

    REAL(wp), INTENT(INOUT) ::  &    !< calculated vertical tracer mass flux
      &  p_mflx_tracer_v(:,:,:)      !< dim: (nproma,nlevp1,nblks_c)
                                     !< [kg m^-2 s^-1]

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp) ::                      & !< fraction which must multiply all
      &  r_m(nproma,ptr_patch%nlev)    !< outgoing fluxes of cell jc
                                       !< to guarantee positive definiteness

    REAL(wp) :: p_m(nproma)            !< sum of fluxes out of cell
                                       !< [kg m^-2]

    REAL(wp) :: z_signum(nproma)       !< sign of mass flux
                                       !< >0: upward; <0: downward

    INTEGER  :: slev, elev             !< vertical start and end level
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: jk, jb, jc     !< index of edge, vert level, block, cell
    INTEGER  :: jkp1, jkm1

  !-------------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = ptr_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,ptr_patch%n_childdom)


    IF (p_test_run) THEN
      r_m = 0._wp
    ENDIF


    i_startblk   = ptr_patch%cells%start_blk(i_rlstart,1)
    i_endblk     = ptr_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,jkp1,p_m,r_m,jkm1,z_signum)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      !
      ! 1. Compute total outward mass
      !
      DO jk = slev, elev
        jkp1 = jk+1

        DO jc = i_startidx, i_endidx

          ! Sum of all outgoing fluxes out of cell jk
          p_m(jc) = p_dtime                                  &
            &     * (MAX(0._wp,p_mflx_tracer_v(jc,jk,jb))    &  ! upper half level
            &      - MIN(0._wp,p_mflx_tracer_v(jc,jkp1,jb)) )   ! lower half level

          ! fraction which must multiply the fluxes out of cell jk to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease \rho^n q^n
          r_m(jc,jk) = MIN(1._wp, (p_cc(jc,jk,jb)*ptr_delp_mc_now(jc,jk,jb)) &
            &         /(p_m(jc) + dbl_eps) )

        ENDDO
      ENDDO

      !
      ! 2. Limit outward fluxes (loop over half levels)
      !    Choose r_m depending on the sign of p_mflx_tracer_v
      !
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev+1, elev
          jkm1 = jk-1
#else
      DO jk = slev+1, elev
        jkm1 = jk-1
        DO jc = i_startidx, i_endidx
#endif
          ! NH:
          ! p_mflx_tracer_v(k-1/2) > 0: flux directed from cell k   -> k-1
          ! p_mflx_tracer_v(k-1/2) < 0: flux directed from cell k-1 -> k
          !
          z_signum(jc) = SIGN(1._wp,p_mflx_tracer_v(jc,jk,jb))

          p_mflx_tracer_v(jc,jk,jb) =  p_mflx_tracer_v(jc,jk,jb)  * 0.5_wp    &
            &                       * ( (1._wp + z_signum(jc)) * r_m(jc,jk)   &
            &                       +   (1._wp - z_signum(jc)) * r_m(jc,jkm1) )
  
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE vflx_limiter_pd



  !-------------------------------------------------------------------------
  !>
  !! Positive definite flux limiter for vertical advection
  !!
  !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
  !! for the hydrostatic core. Only outward fluxes are re-scaled, in 
  !! order to maintain positive definiteness.
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
  !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
  !!   the cubed sphere grid.  J. Comput. Phys., 230, 1215-1237
  !! - Smolarkiewicz, P. K., 1989: Comment on "A positive definite advection 
  !!   scheme obtained by nonlinear renormalization of the advective fluxes.", 
  !!   Mon. Wea. Rev., 117, 2626-2632
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2011-01-07)
  !!
  SUBROUTINE vflx_limiter_pd_ha( ptr_patch, p_dtime, p_cc, p_mflx_tracer_v, &
    &                            opt_rlstart, opt_rlend, opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  ptr_patch

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable at time (n)
      &  p_cc(:,:,:)                 !< dim: (nproma,nlev,nblks_c)
                                     !< [kg kg^-1]

    REAL(wp), INTENT(IN) ::     &    !< time step [s]
      &  p_dtime

    REAL(wp), INTENT(INOUT) ::  &    !< calculated vertical tracer mass flux
      &  p_mflx_tracer_v(:,:,:)      !< dim: (nproma,nlevp1,nblks_c)
                                     !< [kg m^-2 s^-1]

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp) ::                      & !< fraction which must multiply all
      &  r_m(nproma,ptr_patch%nlev)    !< outgoing fluxes of cell jc
                                       !< to guarantee positive definiteness

    REAL(wp) :: p_m(nproma)            !< sum of fluxes out of cell
                                       !< [kg m^-2]

    REAL(wp) :: z_signum(nproma)       !< sign of mass flux
                                       !< >0: upward; <0: downward

    INTEGER  :: slev, elev             !< vertical start and end level
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: jk, jb, jc     !< index of edge, vert level, block, cell
    INTEGER  :: jkp1, jkm1

  !-------------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = ptr_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,ptr_patch%n_childdom)


    IF (p_test_run) THEN
      r_m = 0._wp
    ENDIF


    i_startblk  = ptr_patch%cells%start_blk(i_rlstart,1)
    i_endblk    = ptr_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,jkp1,p_m,r_m,jkm1,z_signum)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,    &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      !
      ! 1. Compute total outward mass
      !
      DO jk = slev, elev
        jkp1 = jk+1

        DO jc = i_startidx, i_endidx

          ! Sum of all outgoing fluxes out of cell jk
          p_m(jc) = p_dtime                                  &
            &     * (MAX(0._wp,p_mflx_tracer_v(jc,jkp1,jb))  &  ! upper half level
            &      - MIN(0._wp,p_mflx_tracer_v(jc,jk,jb)) )     ! lower half level

          ! fraction which must multiply the fluxes out of cell jk to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease \rho^n q^n
          r_m(jc,jk) = MIN(1._wp, (p_cc(jc,jk,jb)*ptr_delp_mc_now(jc,jk,jb)) &
            &         /(p_m(jc) + dbl_eps) )

        ENDDO
      ENDDO

      !
      ! 2. Limit outward fluxes (loop over half levels)
      !    Choose r_m depending on the sign of p_mflx_tracer_v
      !
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev+1, elev
          jkm1 = jk-1
#else
      DO jk = slev+1, elev
        jkm1 = jk-1
        DO jc = i_startidx, i_endidx
#endif

          ! HA:
          ! p_mflx_tracer_v(k-1/2) > 0: flux directed from cell k-1 -> k
          ! p_mflx_tracer_v(k-1/2) < 0: flux directed from cell k   -> k-1
          z_signum(jc) = SIGN(1._wp,p_mflx_tracer_v(jc,jk,jb))

          p_mflx_tracer_v(jc,jk,jb) =  p_mflx_tracer_v(jc,jk,jb)  * 0.5_wp    &
            &                       * ( (1._wp + z_signum(jc)) * r_m(jc,jkm1) &
            &                       +   (1._wp - z_signum(jc)) * r_m(jc,jk) )
  
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE vflx_limiter_pd_ha


  !-------------------------------------------------------------------------
  !>
  !! Monotonous slope limiter for MIURA advection scheme
  !!
  !! Gradient limitation using monotonous Barth-Jesperson limiter
  !! to avoid over/undershoots in the advected field. The limiter is modified
  !! in the sense, that the points where the reconstruction is evaluated
  !! (the Gauss points) have been made time dependent. This wass necessary in
  !! order to achieve positive definiteness.
  !!
  !! @par Literature
  !! - T. J. Barth and D. C. Jespersen. The design and application of upwind
  !!   schemes on unstructured meshes. In 27th Aerospace Sciences Meeting, Jan.
  !!   1989. AIAA paper 89-0366.
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2010-02-26)
  !!
  SUBROUTINE h_miura_slimiter_mo( p_patch, p_cc, ptr_cc, p_distv_gausspoint, p_grad, &
    &                              opt_rlstart, opt_rlend, opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &      !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::    &       !< advected cell centered variable
      &  p_cc(:,:,:)                   !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &       !< different from p_cc if conservative lsq
      &  ptr_cc(:,:,:)                 !< reconstruction is chosen
                                       !< In this case ptr_cc contains coefficient
                                       !< c0 from the lsq reconstruction
                                       !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &       !< distance vectors cell center -->
      &  p_distv_gausspoint(:,:,:,:,:) !< shifted Gauss-points
                                       !< dim: (nproma,nlev,p_patch%nblks_c,3,2)

    REAL(wp), INTENT(INOUT) :: &    !< reconstructed gradients, to be limited
      &  p_grad(:,:,:,:)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    !--- Local variables ---!
    REAL(wp) ::            &
      &  z_min(nproma),    &        !< minimum of cell centered value
      &  z_max(nproma)              !< maximum of cell centered value

    REAL(wp) :: z_lext_val_1,  &    !< linear extrapolation value of cell
      &  z_lext_val_2, z_lext_val_3 !< for edges 1, 2, 3

    REAL(wp) :: z_limit             !< value of gradient limiter

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iilc, iibc

    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: jk, jb, jc         !< index of vert level, block, cell
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

  !-------------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 2
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    ! pointers to line and block indices of three neighbor cells
    iilc => p_patch%cells%neighbor_idx
    iibc => p_patch%cells%neighbor_blk


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_min,z_max,z_lext_val_1,   &
!$OMP            z_lext_val_2,z_lext_val_3,z_limit)

    ! Calculate reconstructed tracer values at Gauss-points provided and limit
    ! the gradient appropriately. Note that one needs to distinguish beween
    ! the prognostic field (p_cc) and the constant c0 (ptr_cc) from the lsq
    ! reconstruction. These fields differ in the case of a conservative least
    ! squares reconstruction.
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev

        DO jc = i_startidx, i_endidx
#endif

          ! min value of cell and its neighbors
          z_min(jc) = MIN( p_cc(jc,jk,jb), p_cc(iilc(jc,jb,1),jk,iibc(jc,jb,1)), &
            &              p_cc(iilc(jc,jb,2),jk,iibc(jc,jb,2)),                 &
            &              p_cc(iilc(jc,jb,3),jk,iibc(jc,jb,3)) )
          ! max value of cell and its neighbors
          z_max(jc) = MAX( p_cc(jc,jk,jb), p_cc(iilc(jc,jb,1),jk,iibc(jc,jb,1)), &
            &              p_cc(iilc(jc,jb,2),jk,iibc(jc,jb,2)),                 &
            &              p_cc(iilc(jc,jb,3),jk,iibc(jc,jb,3)) )


          ! value for linear extrapolation of cell centered value to Gauss-point 1
          z_lext_val_1 = ( ptr_cc(jc,jk,jb)                                      &
            &          + p_grad(jc,jk,jb,1) * p_distv_gausspoint(jc,jk,jb,1,1)   &
            &          + p_grad(jc,jk,jb,2) * p_distv_gausspoint(jc,jk,jb,1,2) ) &
            &          - p_cc(jc,jk,jb)

          ! value for linear extrapolation of cell centered value to Gauss-point 2
          z_lext_val_2 = ( ptr_cc(jc,jk,jb)                                      &
            &          + p_grad(jc,jk,jb,1) * p_distv_gausspoint(jc,jk,jb,2,1)   &
            &          + p_grad(jc,jk,jb,2) * p_distv_gausspoint(jc,jk,jb,2,2) ) &
            &          - p_cc(jc,jk,jb)

          ! value for linear extrapolation of cell centered value to Gauss-point 3
          z_lext_val_3 = ( ptr_cc(jc,jk,jb)                                      &
            &          + p_grad(jc,jk,jb,1) * p_distv_gausspoint(jc,jk,jb,3,1)   &
            &          + p_grad(jc,jk,jb,2) * p_distv_gausspoint(jc,jk,jb,3,2) ) &
            &          - p_cc(jc,jk,jb)



          z_limit = HUGE(1.0_wp)


          ! calculation of limiters
          IF ( z_lext_val_1 > 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_max(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_1) + dbl_eps),z_lext_val_1 ) ) )
          ELSE IF ( z_lext_val_1 < 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_min(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_1) + dbl_eps),z_lext_val_1 ) ) )
          ELSE
            z_limit = MIN( z_limit, 1.0_wp )
          END IF

          IF ( z_lext_val_2 > 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_max(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_2) + dbl_eps),z_lext_val_2 ) ) )
          ELSE IF ( z_lext_val_2 < 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_min(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_2) + dbl_eps),z_lext_val_2 ) ) )
          ELSE
            z_limit = MIN( z_limit, 1.0_wp )
          END IF

          IF ( z_lext_val_3 > 0._wp ) THEN
           z_limit = MIN( z_limit,                          &
              &  MIN( 1.0_wp, (z_max(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_3) + dbl_eps),z_lext_val_3 ) ) )
          ELSE IF ( z_lext_val_3 < 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_min(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_3) + dbl_eps),z_lext_val_3 ) ) )
          ELSE
            z_limit = MIN( z_limit, 1.0_wp )
          END IF


          !
          ! limited value of gradient at cell center
          !
          p_grad(jc,jk,jb,1:2) = z_limit * p_grad(jc,jk,jb,1:2)

          END DO ! end loop over cells

        END DO

      END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE h_miura_slimiter_mo




  !-------------------------------------------------------------------------
  !>
  !! Semi-monotonous slope limiter for MIURA advection scheme
  !!
  !! Gradient limitation using semi-monotonous Barth-Jesperson limiter
  !! to avoid undershoots in the advected field. The limiter is modified
  !! in the sense, that the points where the reconstruction is evaluated
  !! (the Gauss points) have been made time dependent. This was necessary in
  !! order to achieve positive definiteness.
  !!
  !! @par Literature
  !! - T. J. Barth and D. C. Jespersen. The design and application of upwind
  !!   schemes on unstructure meshes. In 27th Aerospace Sciences Meeting, Jan.
  !!   1989. AIAA paper 89-0366.
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2010-02-26)
  !!
  SUBROUTINE h_miura_slimiter_sm( p_patch, p_cc, ptr_cc, p_distv_gausspoint, p_grad, &
    &                             opt_rlstart, opt_rlend, opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::    &   !< advected cell centered variable
      &  p_cc(:,:,:)               !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &   !< different from p_cc if conservative lsq
      &  ptr_cc(:,:,:)             !< reconstruction is chosen
                                   !< In this case ptr_cc contains coefficient
                                   !< c0 from the lsq reconstruction
                                   !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &       !< distance vectors cell center -->
      &  p_distv_gausspoint(:,:,:,:,:) !< shifted Gauss-points
                                       !< dim: (nproma,nlev,p_patch%nblks_c,3,2)

    REAL(wp), INTENT(INOUT) :: &   !< reconstructed gradients, to be limited
      &  p_grad(:,:,:,:)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    !--- Local variables ---!
    REAL(wp) ::            &
      &  z_min(nproma),    &        !< minimum of cell centered value
      &  z_max(nproma)              !< maximum of cell centered value

    REAL(wp) :: z_lext_val_1,  &    !< linear extrapolation value of cell
      &  z_lext_val_2, z_lext_val_3 !< for edges 1, 2, 3

    REAL(wp) :: z_limit             !< value of gradient limiter

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iilc, iibc

    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: jk, jb, jc     !< index of vert level, block, cell
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

  !-------------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 2
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    ! pointers to line and block indices of three neighbor cells
    iilc => p_patch%cells%neighbor_idx
    iibc => p_patch%cells%neighbor_blk


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_min,z_max,z_lext_val_1,   &
!$OMP            z_lext_val_2,z_lext_val_3,z_limit)

    ! Calculate reconstructed tracer values at Gauss-points provided and limit
    ! the gradient appropriately. Note that one needs to distinguish beween
    ! the prognostic field (p_cc) and the constant c0 (ptr_cc) from the lsq
    ! reconstruction. These fields differ in the case of a conservative least
    ! squares reconstruction.
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev

        DO jc = i_startidx, i_endidx
#endif

          ! min value of cell and its neighbors
          z_min(jc) = MIN( p_cc(jc,jk,jb), p_cc(iilc(jc,jb,1),jk,iibc(jc,jb,1)), &
            &              p_cc(iilc(jc,jb,2),jk,iibc(jc,jb,2)),                 &
            &              p_cc(iilc(jc,jb,3),jk,iibc(jc,jb,3)) )
          ! max value of cell and its neighbors
          z_max(jc) = MAX( p_cc(jc,jk,jb), p_cc(iilc(jc,jb,1),jk,iibc(jc,jb,1)), &
            &              p_cc(iilc(jc,jb,2),jk,iibc(jc,jb,2)),                 &
            &              p_cc(iilc(jc,jb,3),jk,iibc(jc,jb,3)) )


          ! value for linear extrapolation of cell centered value to Gauss-point 1
          z_lext_val_1 = ( ptr_cc(jc,jk,jb)                                      &
            &          + p_grad(jc,jk,jb,1) * p_distv_gausspoint(jc,jk,jb,1,1)   &
            &          + p_grad(jc,jk,jb,2) * p_distv_gausspoint(jc,jk,jb,1,2) ) &
            &          - p_cc(jc,jk,jb)

          ! value for linear extrapolation of cell centered value to Gauss-point 2
          z_lext_val_2 = ( ptr_cc(jc,jk,jb)                                      &
            &          + p_grad(jc,jk,jb,1) * p_distv_gausspoint(jc,jk,jb,2,1)   &
            &          + p_grad(jc,jk,jb,2) * p_distv_gausspoint(jc,jk,jb,2,2) ) &
            &          - p_cc(jc,jk,jb)

          ! value for linear extrapolation of cell centered value to Gauss-point 3
          z_lext_val_3 = ( ptr_cc(jc,jk,jb)                                      &
            &          + p_grad(jc,jk,jb,1) * p_distv_gausspoint(jc,jk,jb,3,1)   &
            &          + p_grad(jc,jk,jb,2) * p_distv_gausspoint(jc,jk,jb,3,2) ) &
            &          - p_cc(jc,jk,jb)


          z_limit = HUGE(1.0_wp)

          ! calculation of limiters
          IF ( z_lext_val_1 < 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_min(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_1) + dbl_eps),z_lext_val_1 ) ) )
          ELSE
            z_limit = MIN( z_limit, 1.0_wp )
          END IF

          IF ( z_lext_val_2 < 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_min(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_2) + dbl_eps),z_lext_val_2 ) ) )
          ELSE
            z_limit = MIN( z_limit, 1.0_wp )
          END IF

          IF ( z_lext_val_3 < 0._wp ) THEN
            z_limit = MIN( z_limit,                         &
              &  MIN( 1.0_wp, (z_min(jc) - p_cc(jc,jk,jb))  &
              &  / SIGN( (ABS(z_lext_val_3) + dbl_eps),z_lext_val_3 ) ) )
          ELSE
            z_limit = MIN( z_limit, 1.0_wp )
          END IF

          !
          ! limited value of gradient at cell center
          !
          p_grad(jc,jk,jb,1:2) = z_limit * p_grad(jc,jk,jb,1:2)

          END DO ! end loop over cells

        END DO

      END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE h_miura_slimiter_sm




  !-------------------------------------------------------------------------
  !
  !>
  !! Shift Gauss-points at which the reconstructed tracer fields are
  !! evaluated (for slope-limitation purposes).
  !!
  !! Gauss points (edge midpoints) are shifted by $0.5\Delta t$ times the
  !! reconstructed tangential velocity. Once the shifted Gauss points are known
  !! in geographical coordinates, they are projected (via gnomonic projection)
  !! onto a tangent plane (tangent to cell circumcenter). Then, in this new local
  !! cartesian system the distance vectors between the cell center and the
  !! shifted gauss points are calculated.
  !!
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD, (<2010-03-01>)
  !!
  SUBROUTINE shift_gauss_points( p_patch, p_int, p_dthalf, p_vt_traj,          &
    &                          p_distv_gausspoint, opt_rlstart_c, opt_rlend_c, &
    &                          opt_rlstart_e, opt_rlend_e, opt_slev, opt_elev )


    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), INTENT(IN)  :: p_dthalf  !< \Delta t/2

    REAL(wp), INTENT(IN)  ::  &        !< reconstructed tangential velocity
      &  p_vt_traj(:,:,:)              !< component at edges

    REAL(wp), INTENT(OUT) ::  &        !< distance vectors cell center -->
      &  p_distv_gausspoint(:,:,:,:,:) !< shifted Gauss-points
                                       !< dim: (nproma,nlev,nblks_c,3,2)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart_c                  !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend_c                    !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart_e                  !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend_e                    !< (to avoid calculation of halo points)

    TYPE(t_cartesian_coordinates) :: cc_gausspoint_shifted  !< shifted Gauss-point
                                                          !< in geographical coordinates

    TYPE(t_geographical_coordinates) :: gg_gausspoint_shifted  !< shifted Gauss-point
                                                             !< in geographical coordinates

    REAL(wp) ::  &      !< position of shifted Gauss-points
      &  z_pos_gausspoint_shifted(nproma,p_patch%nlev,p_patch%nblks_e,2)
                        !< (in geographical coordinates)
                        !< dim: (nproma,nlev,p_patch%nblks_e,2)

    REAL(wp) :: xloc,yloc,                   & !< coordinates needed for gnomonic projection
      & xtemp1,xtemp2,xtemp3,ytemp1,ytemp2,  & !< and distance calculation
      & ytemp3

    REAL(wp) :: z_dist_g1(2),  & !< distance vector cell center -> shifted Gauss-point
      &         z_dist_g2(2),  &
      &         z_dist_g3(2)

    REAL(wp) :: z_v_tang_cc(3)  !< reconstructed velocity vector in cartesian coordinates
                                !< (only tangential component)

    REAL(wp) :: v_tang_lon      !< tangential velocity component in lon-direction
    REAL(wp) :: v_tang_lat      !< tangential velocity component in lon-direction

    INTEGER, DIMENSION(:,:,:), POINTER :: & !< Pointer to line and block indices (array)
      &  iile, iibe

    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: je, jk, jb, jc  !< index of edge, vert level, block, cell
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart_c, i_rlend_c, i_rlstart_e, i_rlend_e, i_nchdom

    !-------------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart_c) ) THEN
      i_rlstart_c = opt_rlstart_c
    ELSE
      i_rlstart_c = 2
    ENDIF

    IF ( PRESENT(opt_rlend_c) ) THEN
      i_rlend_c = opt_rlend_c
    ELSE
      i_rlend_c = min_rlcell_int
    ENDIF

    IF ( PRESENT(opt_rlstart_e) ) THEN
      i_rlstart_e= opt_rlstart_e
    ELSE
      i_rlstart_e= 2
    ENDIF

    IF ( PRESENT(opt_rlend_e) ) THEN
      i_rlend_e  = opt_rlend_e
    ELSE
      i_rlend_e  = min_rledge
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    ! pointers to line and block indices of cell edges
    iile => p_patch%cells%edge_idx
    iibe => p_patch%cells%edge_blk


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    i_startblk = p_patch%edges%start_blk(i_rlstart_e,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend_e,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,v_tang_lon,            &
!$OMP            v_tang_lat,z_v_tang_cc,cc_gausspoint_shifted,       &
!$OMP            gg_gausspoint_shifted)

    ! Calculate shifted Gauss-points, by calculating backward trajectories
    ! based on the reconstructed tangential velocity component and p_dthalf.
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,          &
                         i_startidx, i_endidx, i_rlstart_e, i_rlend_e)
!CDIR NOLOOPCHG
      DO jk = slev, elev

        DO je = i_startidx, i_endidx

          ! zonal component of tangential velocity at edge midpoint
          v_tang_lon =                                                   &
            &  p_vt_traj(je,jk,jb) * p_patch%edges%dual_normal(je,jb)%v1

          ! meridional component of tangential velocity at edge midpoint
          v_tang_lat    =                                                 &
            &   p_vt_traj(je,jk,jb) * p_patch%edges%dual_normal(je,jb)%v2

          ! velocity vector: convert from geographical coordinates to
          ! velocity vector in cartesian coordinates
          CALL gvec2cvec( v_tang_lon,v_tang_lat,                       &! in
            &             p_patch%edges%center(je,jb)%lon,             &! in
            &             p_patch%edges%center(je,jb)%lat,             &! in
            &             z_v_tang_cc(1),z_v_tang_cc(2),z_v_tang_cc(3) &! inout
            &            )

          ! calculate shifted Gauss-point in cartesian coordinates
          ! p_int%cart_edge_coord carries cartesian coordinates of edge midpoints
          cc_gausspoint_shifted%x(1) = re*p_int%cart_edge_coord(je,jb,1) &
            &                        - p_dthalf * z_v_tang_cc(1)
          cc_gausspoint_shifted%x(2) = re*p_int%cart_edge_coord(je,jb,2) &
            &                        - p_dthalf * z_v_tang_cc(2)
          cc_gausspoint_shifted%x(3) = re*p_int%cart_edge_coord(je,jb,3) &
            &                        - p_dthalf * z_v_tang_cc(3)

          ! shifted Gauss-point: convert from cartesian to geographical coordinates
          gg_gausspoint_shifted=cc2gc(cc_gausspoint_shifted)

          ! Save shifted Gauss-points in edge based data structure
          z_pos_gausspoint_shifted(je,jk,jb,1) = gg_gausspoint_shifted%lon
          z_pos_gausspoint_shifted(je,jk,jb,2) = gg_gausspoint_shifted%lat

        ENDDO ! loop over cells

      ENDDO   ! loop over vertical levels

    END DO    ! loop over blocks
!$OMP END DO


    i_startblk = p_patch%cells%start_blk(i_rlstart_c,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend_c,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,xloc,yloc,xtemp1,xtemp2,xtemp3,  &
!$OMP            ytemp1,ytemp2,ytemp3,z_dist_g1,z_dist_g2,z_dist_g3)

    ! Calculate distance vectors between cell center and shifted gauss points
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,          &
                         i_startidx, i_endidx, i_rlstart_c, i_rlend_c)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR NOLOOPCHG
      DO jk = slev, elev

        DO jc = i_startidx, i_endidx
#endif
          !
          ! Calculate distance vectors between cell center and shifted
          ! gauss-points at edge 1,2,3
          ! In order to do so, project shifted gauss points onto plane tangent
          ! to cell center
          !

          ! cell center
          xloc = p_patch%cells%center(jc,jb)%lon
          yloc = p_patch%cells%center(jc,jb)%lat

          !
          ! Distance vector between cell center and shifted gauss-point at edge 1
          !
          xtemp1 = z_pos_gausspoint_shifted(iile(jc,jb,1),jk,iibe(jc,jb,1),1)
          ytemp1 = z_pos_gausspoint_shifted(iile(jc,jb,1),jk,iibe(jc,jb,1),2)
          CALL gnomonic_proj( xloc, yloc, xtemp1, ytemp1,  & !in
            &                z_dist_g1(1), z_dist_g1(2) )    !out

          ! multiply with earth radius and save in cell based data structure
          p_distv_gausspoint(jc,jk,jb,1,1:2) = re * z_dist_g1(1:2)

          !
          ! Distance vector between cell center and shifted gauss-point at edge 2
          !
          xtemp2 = z_pos_gausspoint_shifted(iile(jc,jb,2),jk,iibe(jc,jb,2),1)
          ytemp2 = z_pos_gausspoint_shifted(iile(jc,jb,2),jk,iibe(jc,jb,2),2)
          CALL gnomonic_proj( xloc, yloc, xtemp2, ytemp2,  & !in
            &                z_dist_g2(1), z_dist_g2(2) )    !out

          ! multiply with earth radius and save in cell based data structure
          p_distv_gausspoint(jc,jk,jb,2,1:2) = re * z_dist_g2(1:2)

          !
          ! Distance vector between cell center and shifted gauss-point at edge 3
          !
          xtemp3 = z_pos_gausspoint_shifted(iile(jc,jb,3),jk,iibe(jc,jb,3),1)
          ytemp3 = z_pos_gausspoint_shifted(iile(jc,jb,3),jk,iibe(jc,jb,3),2)
          CALL gnomonic_proj( xloc, yloc, xtemp3, ytemp3,  & !in
            &                z_dist_g3(1), z_dist_g3(2) )    !out

          ! multiply with earth radius and save in cell based data structure
          p_distv_gausspoint(jc,jk,jb,3,1:2) = re * z_dist_g3(1:2)

        END DO ! end loop over cells

      END DO

    END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE shift_gauss_points




  !-------------------------------------------------------------------------
  !>
  !! Barth-Jespersen slope limiter for MUSCL (2nd order) vertical advection
  !! (monotone version)
  !!
  !! Limits slope of reconstruced vertical gradients using the
  !! Barth-Jespersen slope limiter (monotone version). Avoids non-physical
  !! over/undershoots in advected fields.
  !!
  !! Literature
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-01-20)
  !!
  SUBROUTINE v_muscl_slimiter_mo( p_patch, p_cc, p_gmoment, p_delp_c1,       &
    &                             p_delp_c2, p_grad, opt_rlstart, opt_rlend, &
    &                             opt_slev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  & !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::    &   !< advected cell centered variable
      &  p_cc(:,:,:)               !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(INOUT) :: &   !< gradient at cell center
      &  p_grad(:,:,:)             !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &   !< pressure difference between half level and
      &  p_delp_c1(:,:,:),     &   !< corresponding upper (p_delp_c1) and lower
      &  p_delp_c2(:,:,:)          !< full level (p_delp_c2).
                                   !< (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::    &   !< geometrical moment which accounts for
      &  p_gmoment(:,:,:)          !< offcentering of mass point at timestep n
                                   !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp) :: z_top, z_mid, z_bot        !< cell values of transported field
    REAL(wp) :: z_min, z_max               !< max and min values of 3 levels
    REAL(wp) :: z_lext_val_1, z_lext_val_2 !< linear extrapolation increment
    REAL(wp) :: z_limit                    !< limiter for gradient at cell center
    INTEGER  :: nlev                       !< number of full levels
    INTEGER  :: slev                       !< vertical start level
    INTEGER  :: jc, jk, jb                 !< index of cell, vertical level and block
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: jkm1, jkm1_ic,  &          !< vertical level minus and plus one
      &         jkp1, jkp1_ic

  !-------------------------------------------------------------------------

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev  = opt_slev
    ELSE
      slev  = 1
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev = p_patch%nlev

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,jkm1_ic,jkm1,jkp1_ic,jkp1,z_top, &
!$OMP            z_mid,z_bot,z_min,z_max,z_lext_val_1,z_lext_val_2,z_limit)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )


      DO jk = slev, nlev

        ! index of top half level
        jkm1_ic = jk - 1
        jkm1 = MAX( jkm1_ic, 1 )
        ! index of bottom half level
        jkp1_ic = jk + 1
        jkp1 = MIN( jkp1_ic, nlev )


        DO jc = i_startidx, i_endidx

          z_top = p_cc(jc,jkm1,jb)
          z_mid = p_cc(jc,jk  ,jb)
          z_bot = p_cc(jc,jkp1,jb)

          ! minimum value of this and surrounding cells
          z_min = MIN( z_top, z_mid, z_bot )
          ! maximum value of this and surrounding cells
          z_max = MAX( z_top, z_mid, z_bot )


          ! linear extrapolated value
          ! first (to face above)
          z_lext_val_2 = p_grad(jc,jk,jb)                              &
            &          * ( p_delp_c2(jc,jk,jb) - p_gmoment(jc,jk,jb) )


          ! calculation of limiter
          IF ( z_lext_val_2 > 0._wp ) THEN
            z_limit = MIN( 1._wp, (z_max-z_mid)/  &
              &       SIGN( (ABS(z_lext_val_2) + dbl_eps),z_lext_val_2 ))

          ELSE IF ( z_lext_val_2 < 0._wp ) THEN
            z_limit = MIN( 1._wp, (z_min-z_mid)/  &
              &       SIGN( (ABS(z_lext_val_2) + dbl_eps),z_lext_val_2 ))
          ELSE
            z_limit = 1._wp
          END IF


          ! linear extrapolated value
          ! second (to face below)
          z_lext_val_1 = p_grad(jc,jk,jb)                                  &
            &          * ( p_delp_c1(jc,jkp1_ic,jb) - p_gmoment(jc,jk,jb) )


          ! calculation of limiter
          IF ( z_lext_val_1 > 0._wp ) THEN
            z_limit = MIN( z_limit,  &
              MIN( 1._wp, (z_max-z_mid)/  &
                &  SIGN( (ABS(z_lext_val_1) + dbl_eps),z_lext_val_1 )))

          ELSE IF ( z_lext_val_1 < 0._wp ) THEN
            z_limit = MIN( z_limit,  &
              MIN( 1._wp, (z_min-z_mid)/  &
                &  SIGN( (ABS(z_lext_val_1) + dbl_eps),z_lext_val_1 )))
          ELSE
            z_limit = MIN( z_limit, 1._wp )
          END IF


          ! limited vertical gradient
          p_grad(jc,jk,jb) = z_limit * p_grad(jc,jk,jb)

        END DO ! end loop over cells

      ENDDO ! end loop over vertical levels

    ENDDO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE v_muscl_slimiter_mo




  !-------------------------------------------------------------------------
  !>
  !! Barth-Jespersen slope limiter for MUSCL (2nd order) vertical advection
  !! (semi-monotonic version)
  !!
  !! Limits slope of reconstruced vertical gradients using the
  !! Barth-Jespersen slope limiter (positive definite version). Avoids
  !! non-physical undershoots in advected fields.
  !!
  !! Literature
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-01-20)
  !!
  SUBROUTINE v_muscl_slimiter_sm( p_patch, p_cc, p_gmoment, p_delp_c1,       &
    &                             p_delp_c2, p_grad, opt_rlstart, opt_rlend, &
    &                             opt_slev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &    !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(INOUT) :: &  !< gradient at cell center
      &  p_grad(:,:,:)            !< dim: (nproma,nlev,p_patch%nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< pressure difference between half level and
      &  p_delp_c1(:,:,:),   &    !< corresponding upper (p_delp_c1) and lower
      &  p_delp_c2(:,:,:)         !< full level (p_delp_c2).
                                  !< (nproma,nlevp1,p_patch%nblks_c)

    REAL(wp), INTENT(IN) :: &     !< geometrical moment which accounts for
      &  p_gmoment(:,:,:)         !< offcentering of mass point at timestep n
                                  !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp) :: z_top, z_mid, z_bot        !< cell values of transported field
    REAL(wp) :: z_min                      !< max and min values of 3 levels
    REAL(wp) :: z_lext_val_1, z_lext_val_2 !< linear extrapolation increment
    REAL(wp) :: z_limit                    !< limiter for gradient at cell center
    INTEGER  :: nlev                       !< number of full levels
    INTEGER  :: slev                       !< vertical start level
    INTEGER  :: jc, jk, jb                 !< index of cell, vertical level and block
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: jkm1, jkm1_ic,  &          !< vertical level minus and plus one
      &         jkp1, jkp1_ic

    !-----------------------------------------------------------------------

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev  = opt_slev
    ELSE
      slev  = 1
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    ! number of vertical levels
    nlev = p_patch%nlev

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,jkm1_ic,jkm1,jkp1_ic,jkp1,z_top,  &
!$OMP            z_mid,z_bot,z_min,z_lext_val_1,z_lext_val_2,z_limit)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )


      DO jk = slev, nlev

        ! index of top half level
        jkm1_ic = jk - 1
        jkm1 = MAX( jkm1_ic, 1 )
        ! index of bottom half level
        jkp1_ic = jk + 1
        jkp1 = MIN( jkp1_ic, nlev )


        DO jc = i_startidx, i_endidx

          z_top = p_cc(jc,jkm1,jb)
          z_mid = p_cc(jc,jk  ,jb)
          z_bot = p_cc(jc,jkp1,jb)

          ! minimum value of this and surrounding cells
          z_min = MIN( z_top, z_mid, z_bot )


          ! linear extrapolated value
          ! first (to face above)
          z_lext_val_2 = p_grad(jc,jk,jb)                              &
            &          * ( p_delp_c2(jc,jk,jb) - p_gmoment(jc,jk,jb) )

          ! calculation of limiter
          IF ( z_lext_val_2 < 0._wp ) THEN
            z_limit = MIN( 1._wp, (z_min-z_mid)/                &
              &       SIGN( (ABS(z_lext_val_2) + dbl_eps),z_lext_val_2 ))
          ELSE
            z_limit = 1._wp
          END IF


          ! linear extrapolated value
          ! second (to face below)
          z_lext_val_1 = p_grad(jc,jk,jb)                                  &
            &          * ( p_delp_c1(jc,jkp1_ic,jb) - p_gmoment(jc,jk,jb) )

          ! calculation of limiter
          IF ( z_lext_val_1 < 0._wp ) THEN
            z_limit = MIN( z_limit, MIN( 1._wp, (z_min-z_mid)/  &
                &  SIGN( (ABS(z_lext_val_1) + dbl_eps),z_lext_val_1 )))
          ELSE
            z_limit = MIN( z_limit, 1._wp )
          END IF


          ! limited vertical gradient
          p_grad(jc,jk,jb) = z_limit * p_grad(jc,jk,jb)

        END DO ! end loop over cells

      ENDDO ! end loop over vertical levels

    ENDDO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE v_muscl_slimiter_sm




  !-------------------------------------------------------------------------
  !>
  !! Limiter for PPM (3rd order) vertical advection (monotone version)
  !!
  !! Removes over- and undershoots in first guess parabola by resetting the
  !! upper or lower interface values.
  !! Avoids non-physical over/undershoots in advected fields.
  !!
  !! Note that this limiter was coded assuming a pressure based vertical
  !! coordinate system. Nevertheless this limiter works for a height based
  !! vertical system, too. This is due to a 'wrong' computation of z_delta
  !! in the case of a height based coordinate system (i.e. z_delta is
  !! implicity multiplied by -1)
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-02-04)
  !!
  SUBROUTINE v_ppm_slimiter_mo( p_patch, p_cc, p_face, p_slope, p_face_up,   &
    &                           p_face_low, opt_rlstart, opt_rlend, opt_slev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &     !< advected cell centered variable
      &  p_cc(:,:,:)               !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &     !< reconstructed face values of the advected field
      &  p_face(:,:,:)             !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &     !< monotonized slope
      &  p_slope(:,:,:)            !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(INOUT) :: &   !< final face value (upper face, height based)
      &  p_face_up(:,:,:)          !< dim: (nproma,nlevp,nblks_c)

    REAL(wp), INTENT(INOUT) :: &   !< final face value (lower face, height based)
      &  p_face_low(:,:,:)         !< dim: (nproma,nlevp,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    INTEGER  :: nlev                      !< number of full levels
    INTEGER  :: slev                      !< vertical start level
    INTEGER  :: jc, jk, jb                !< index of cell, vertical level and block
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: ikp1                      !< vertical level plus one

    REAL(wp) :: z_delta                   !< lower minus upper face value
    REAL(wp) :: z_a6i                     !< curvature of parabola

    !-----------------------------------------------------------------------

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev  = opt_slev
    ELSE
      slev  = 1
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev = p_patch%nlev

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ikp1,z_delta,z_a6i)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )


      DO jk = slev, nlev

        ! index of bottom half level
        ikp1 = jk + 1

        DO jc = i_startidx, i_endidx

          z_delta   = p_face(jc,ikp1,jb) - p_face(jc,jk,jb)
          z_a6i     = 6._wp * (p_cc(jc,jk,jb)                           &
            &       - 0.5_wp * (p_face(jc,jk,jb) + p_face(jc,ikp1,jb)))


          IF ( p_slope(jc,jk,jb) == 0._wp) THEN
            p_face_up(jc,jk,jb)  = p_cc(jc,jk,jb)
            p_face_low(jc,jk,jb) = p_cc(jc,jk,jb)

          ELSE IF (z_delta * z_a6i > z_delta * z_delta) THEN
            p_face_up(jc,jk,jb)  = 3._wp*p_cc(jc,jk,jb) - 2._wp*p_face(jc,ikp1,jb)
            p_face_low(jc,jk,jb) = p_face(jc,ikp1,jb)

          ELSE IF (z_delta * z_a6i < -1._wp * (z_delta * z_delta)) THEN
            p_face_up(jc,jk,jb)  = p_face(jc,jk,jb)
            p_face_low(jc,jk,jb) = 3._wp*p_cc(jc,jk,jb) - 2._wp*p_face(jc,jk,jb)

          ELSE
            p_face_up(jc,jk,jb)  = p_face(jc,jk,jb)
            p_face_low(jc,jk,jb) = p_face(jc,ikp1,jb)
          ENDIF

        END DO

      END DO

    END DO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE v_ppm_slimiter_mo




  !-------------------------------------------------------------------------
  !>
  !! Limiter for PPM (3rd order) vertical advection (semi-monotone version)
  !!
  !! Removes undershoots in first guess parabola by resetting either the
  !! upper or lower interface value.
  !! Avoids non-physical undershoots in advected fields.
  !!
  !! Note that this limiter was coded assuming a pressure based vertical
  !! coordinate system. Nevertheless this limiter works for a height based
  !! vertical system, too - without any modifications.
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-02-04)
  !!
  SUBROUTINE v_ppm_slimiter_sm( p_patch, p_cc, p_face, p_face_up,            &
    &                           p_face_low, opt_rlstart, opt_rlend, opt_slev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &    !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &     !< advected cell centered variable
      &  p_cc(:,:,:)               !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &     !< reconstructed face values of the advected field
      &  p_face(:,:,:)             !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(INOUT) :: &   !< final face value (upper face, height based)
      &  p_face_up(:,:,:)          !< dim: (nproma,nlevp,nblks_c)

    REAL(wp), INTENT(INOUT) :: &   !< final face value (lower face, height based)
      &  p_face_low(:,:,:)         !< dim: (nproma,nlevp,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      & opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp) :: z_delta                   !< lower minus upper face value
    REAL(wp) :: z_a6i                     !< curvature of parabola

    INTEGER  :: nlev                      !< number of full levels
    INTEGER  :: slev                      !< vertical start level
    INTEGER  :: jc, jk, jb                !< index of cell, vertical level and block
    INTEGER  :: ikp1                      !< vertical level plus one
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

  !-------------------------------------------------------------------------

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev  = opt_slev
    ELSE
      slev  = 1
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev = p_patch%nlev

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ikp1,z_delta,z_a6i)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      DO jk = slev, nlev

        ! index of bottom half level
        ikp1 = jk + 1

        DO jc = i_startidx, i_endidx

          z_delta   = p_face(jc,ikp1,jb) - p_face(jc,jk,jb)
          z_a6i     = 6._wp * (p_cc(jc,jk,jb)                           &
            &       - 0.5_wp * (p_face(jc,jk,jb) + p_face(jc,ikp1,jb)))


          IF (ABS(z_delta) < -1._wp*z_a6i) THEN

            IF (p_cc(jc,jk,jb) < MIN(p_face(jc,jk,jb),p_face(jc,ikp1,jb)) ) THEN
              p_face_up(jc,jk,jb)  = p_cc(jc,jk,jb)
              p_face_low(jc,jk,jb) = p_cc(jc,jk,jb)

            ELSE

              IF (p_face(jc,jk,jb) > p_face(jc,ikp1,jb)) THEN
                p_face_up(jc,jk,jb)  = 3._wp*p_cc(jc,jk,jb) - 2._wp*p_face(jc,ikp1,jb)
                p_face_low(jc,jk,jb) = p_face(jc,ikp1,jb)

              ELSE
                p_face_up(jc,jk,jb)  = p_face(jc,jk,jb)
                p_face_low(jc,jk,jb) = 3._wp*p_cc(jc,jk,jb) - 2._wp*p_face(jc,jk,jb)

              ENDIF

            ENDIF

          ELSE
            p_face_up(jc,jk,jb)  = p_face(jc,jk,jb)
            p_face_low(jc,jk,jb) = p_face(jc,ikp1,jb)
          ENDIF

        END DO

      END DO

    END DO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE v_ppm_slimiter_sm


END MODULE mo_advection_limiter

