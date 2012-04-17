!>
!! mo_divergent_modes
!!
!! This module performs the time stepping of the divergent modes -- the fast
!! modes of the nh-system. It comprises the acoustic and gravity wave propa-
!! gation, which are described by the two Poisson brackets: the mass bracket
!! (includes the gradient of the specific mechanical energy and the
!! mass flux divergence) and the theta bracket (includes the pressure gradient
!! term and the flux divergence of the potential temperature).
!! The solving procedure is explicit in the horizontal with the usual forward-
!! backward approach: for approximate energetic consistency a kinetic energy
!! estimate is needed for the future time step. This must be gained in a
!! numerically stable manner. The vertical approach is taken as suggested in
!! Gassmann and Herzog 2008. Clearly, all this is a compromise for a reasonable
!! effiency. Conservation of energy can no longer be expected fully.
!!
!! @author Almut Gassmann, MPI-M
!!
!! @par Revision History
!! Initial release by Almut Gassmann (2009-05-12)
!! Modification by Almut Gassmann, MPI-M (2009-11-17)
!! - solve for w instead for mass flux
!! - call predictor step from outside the routine
!! - add Rayleigh damping
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_divergent_modes

  USE mo_kind,                  ONLY: wp
  USE mo_nonhydrostatic_config, ONLY: gmres_rtol_nh, upstr_beta, ltheta_up_hori
  USE mo_parallel_config,       ONLY: nproma
  USE mo_run_config,            ONLY: dtime
  USE mo_model_domain,          ONLY: t_patch
  USE mo_grid_config,           ONLY: lplane
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp,                  ONLY: cells2edges_scalar, edges2cells_scalar
  USE mo_nonhydro_types,        ONLY: t_nh_state, t_nh_metrics
  USE mo_physical_constants,    ONLY: cpd, rd, cvd, cvd_o_rd, rd_o_cpd
  USE mo_math_gradients,        ONLY: grad_fd_norm, grad_dir_edge
  USE mo_math_divrot,           ONLY: div
  USE mo_math_laplace,          ONLY: directional_laplace
  USE m_gmres,                  ONLY: gmres
  USE mo_exception,             ONLY: finish
  USE mo_nh_init_utils,         ONLY: nflat
  USE mo_sync,                  ONLY: SYNC_E, SYNC_C, sync_patch_array, &
    &                                 sync_patch_array_mult

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: divergent_modes_5band, impl_vert_adv_theta, impl_vert_adv_theta_5diag

  CONTAINS

  !>
  !! divergent_modes_5band
  !!
  !! This is the solver for the divergent modes that solves a 5band matrix
  !! equation in the vertical. Further details are to be bound in the scientific
  !! documentation (file: vertical_solver.pdf(tex))
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2011-01-03)
  !!
  SUBROUTINE divergent_modes_5band (p_nh, p_patch, p_int, know, knew, l_predictor)

    TYPE(t_nh_state),INTENT(INOUT):: p_nh
    TYPE(t_int_state),TARGET,INTENT(IN) :: p_int
    TYPE(t_patch),TARGET,INTENT(IN):: p_patch
    LOGICAL,         INTENT(IN) :: l_predictor
    INTEGER,         INTENT(IN) :: know, knew

    INTEGER  :: nblks_c, nblks_e, npromz_c, npromz_e, nlen, jb, jk
    INTEGER  :: nlev, nlevp1              !< number of full and half levels

    REAL(wp) :: z_theta_v_e     (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_ddxn_emech_e  (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_theta_v_fl_e  (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_grad_corr_e_k (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_jn_vn_e       (nproma,p_patch%nlevp1,p_patch%nblks_e)

    REAL(wp) :: z_mass_fl_div   (nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_exner_mean    (nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_emech_kin2d   (nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_emech_kin3d   (nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_theta_v_fl_div(nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_grad_corr_c_k (nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_grad_corr_c_l (nproma,p_patch%nlevp1,p_patch%nblks_c), &
                z_theta_v_l     (nproma,p_patch%nlevp1,p_patch%nblks_c), &
                z_th_ddz_exner_c(nproma,p_patch%nlevp1,p_patch%nblks_c), &
                z_contr_corr_c  (nproma,p_patch%nlevp1,p_patch%nblks_c)

    REAL(wp) :: z_w_expl        (nproma,p_patch%nlevp1                ), &
                z_kinw_expl_ic  (nproma,p_patch%nlevp1                ), &
                z_rho_expl_ic   (nproma,p_patch%nlevp1                ), &
                z_kinw_expl_mc  (nproma,p_patch%nlev                  ), &
                z_rho_expl      (nproma,p_patch%nlev                  ), &
                z_exner_expl    (nproma,p_patch%nlev                  )

    REAL(wp) :: z_a             (nproma,p_patch%nlevp1                ), &
                z_b             (nproma,p_patch%nlevp1                ), &
                z_c             (nproma,p_patch%nlevp1                ), &
                z_d             (nproma,p_patch%nlevp1                ), &
                z_e             (nproma,p_patch%nlevp1                ), &
                z_q             (nproma,p_patch%nlevp1                ), &
                z_r             (nproma,p_patch%nlevp1                ), &
                z_m             (nproma,p_patch%nlevp1                ), &
                z_u             (nproma,p_patch%nlevp1                ), &
                z_v             (nproma,p_patch%nlevp1                ), &
                z_alpha         (nproma,p_patch%nlevp1                ), &
                z_beta          (nproma,p_patch%nlev                  ), &
                z_gamma         (nproma,p_patch%nlevp1                ), &
                z_delta         (nproma,p_patch%nlevp1                ), &
                z_epsil         (nproma,p_patch%nlevp1                ), &
                z_gaep          (nproma,p_patch%nlevp1                ), &
                z_gaepm         (nproma,p_patch%nlevp1                ), &
                z_gaepd         (nproma,p_patch%nlevp1                )

    REAL(wp) :: z_k             (nproma                       ), &
                z_l             (nproma                       )

    ! for the GMRES solver (needed for the lower boundary condition)
    REAL(wp), PARAMETER :: z_1o2 = 0.5_wp
    REAL(wp) :: z_gmres_rhs     (nproma,p_patch%nblks_e)
    REAL(wp) :: z_vn_surf       (nproma,p_patch%nblks_e)
    INTEGER, PARAMETER :: nmax_iter = 100  ! max. number of allowed iterations
    REAL(wp):: z_residual(nmax_iter)       ! residuum of every iteration
    LOGICAL :: lmaxiter !.TRUE. on output if nmax_iter was encountered
    INTEGER :: niter    ! number of actually needed iterations

    !-------------------------------------------------------------------

    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c
    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e
    nlev      = p_patch%nlev
    nlevp1    = p_patch%nlevp1

    ! 1) What can be done on cell main and interface levels
    !======================================================

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = 1, nlev
        ! specific mechanical energy
        ! - horizontal kinetic energy with implicit vn
        z_emech_kin2d(1:nlen,jk,jb) = p_nh%diag%e_kinh(1:nlen,jk,jb) &
        &                        + p_nh%metrics%geopot(1:nlen,jk,jb)
        ! - horizontal kinetic energy with vector from RK
        z_emech_kin3d(1:nlen,jk,jb) = p_nh%diag%e_kin (1:nlen,jk,jb) &
        &                        + p_nh%metrics%geopot(1:nlen,jk,jb)
      ENDDO
      IF (l_predictor) THEN ! first call
        DO jk = 1, nlev
          ! off-centered exner pressure only for horizontal direction
          ! (see scientific documentation)
          z_exner_mean(1:nlen,jk,jb) = &
          &  (rd_o_cpd-0.5_wp)*p_nh%diag%exner_old  (1:nlen,jk,jb) &
          & +(1.5_wp-rd_o_cpd)*p_nh%prog(know)%exner(1:nlen,jk,jb)
        ENDDO
        p_nh%diag%rho_ic(1:nlen,1,jb) = 0.0_wp
        DO jk = 2,nlev
          ! density at interface levels
          p_nh%diag%rho_ic(1:nlen,jk,jb) =0.5_wp*(p_nh%prog(know)%rho(1:nlen,jk-1,jb) &
          &                                     + p_nh%prog(know)%rho(1:nlen,jk  ,jb))
        ENDDO
        p_nh%diag%rho_ic(1:nlen,nlevp1,jb) = 0.0_wp
      ELSE ! second call
        DO jk = 1, nlev
          ! store old exner if the current call is the main step
          p_nh%diag%exner_old(1:nlen,jk,jb)=p_nh%prog(know)%exner(1:nlen,jk,jb)
        ENDDO
      ENDIF
      ! theta at interfaces
      z_theta_v_l     (1:nlen,1,jb) = 0.0_wp
        DO jk = 2, nlev
          z_theta_v_l(1:nlen,jk,jb) =0.5_wp &
          & *(p_nh%diag%theta_v_impl(1:nlen,jk-1,jb) &
          & + p_nh%diag%theta_v_impl(1:nlen,jk  ,jb) )
        ENDDO
      z_theta_v_l     (1:nlen,nlevp1,jb) = 0.0_wp
      ! vertical measures
      DO jk = 2, nlev
        ! vertical pressure gradient * theta_v (for w-equation)
        z_th_ddz_exner_c(1:nlen,jk,jb) = z_theta_v_l(1:nlen,jk  ,jb) &
        &                   *( p_nh%prog(know)%exner(1:nlen,jk-1,jb) &
        &                     -p_nh%prog(know)%exner(1:nlen,jk  ,jb))&
        &                 / p_nh%metrics%ddqz_z_half(1:nlen,jk  ,jb)
        ! metric correction for horizontal velocity equation
        z_grad_corr_c_l(1:nlen,jk,jb) = (z_emech_kin3d(1:nlen,jk-1,jb) &
        &                               -z_emech_kin3d(1:nlen,jk  ,jb))&
        &                   / p_nh%metrics%ddqz_z_half(1:nlen,jk  ,jb) &
        &                        +cpd*z_th_ddz_exner_c(1:nlen,jk  ,jb)
      ENDDO
      z_grad_corr_c_l (1:nlen,nlevp1,jb) = p_nh%prog(know)%w(1:nlen,nlevp1,jb)/dtime
      ! metric correction for gradient must be fist averaged vertically
      DO jk = nflat+1, nlev
        z_grad_corr_c_k(1:nlen,jk,jb) = 0.5_wp/p_nh%metrics%ddqz_z_full(1:nlen,jk,jb) &
        &  *(z_grad_corr_c_l(1:nlen,jk  ,jb)*p_nh%metrics%ddqz_z_half(1:nlen,jk  ,jb) &
        &   +z_grad_corr_c_l(1:nlen,jk+1,jb)*p_nh%metrics%ddqz_z_half(1:nlen,jk+1,jb))
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! 2) What can be done on edges
    !=============================
    IF (l_predictor) THEN
      ! edge value of theta
      CALL advection_edges (p_nh%prog(know)%theta_v,p_nh%prog(know)%vn,&
      &                     p_patch,p_nh%metrics,p_int,z_theta_v_e)
      ! horizontal gradient of Exner pressure
      CALL grad_fd_norm(z_exner_mean, p_patch, p_nh%diag%horpgrad)
    ENDIF

    ! Compute horizontal gradient of mechanical energy
    CALL grad_fd_norm(z_emech_kin3d, p_patch, z_ddxn_emech_e)

    ! Compute the horizontal average of the metric correction to the gradient
    CALL cells2edges_scalar(z_grad_corr_c_k, p_patch, p_int%c_lin_e, &
         &                  z_grad_corr_e_k, nflat+1, nlev)

    ! 3) Compute horizontal velocity equation
    !========================================

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      IF ( l_predictor ) THEN  ! the pressure gradient term is the same in both calls
        DO jk = 1, nlev        ! note: the metric correction is not, as it contains the
          p_nh%diag%horpgrad(1:nlen,jk,jb) = &                      ! kinetic energy
          &  cpd*z_theta_v_e(1:nlen,jk,jb)*p_nh%diag%horpgrad(1:nlen,jk,jb)
        ENDDO
      ENDIF
      DO jk = 1, nflat
        p_nh%prog(knew)%vn(1:nlen,jk,jb)=p_nh%prog(know)%vn(1:nlen,jk,jb)&
        &                     +dtime*(p_nh%diag%ddt_vn_vort(1:nlen,jk,jb)&
        &                            +p_nh%diag%ddt_vn_phy (1:nlen,jk,jb)&
        &                            -z_ddxn_emech_e       (1:nlen,jk,jb)&
        &                            -p_nh%diag%horpgrad   (1:nlen,jk,jb))
      ENDDO
      DO jk = nflat+1,nlev
        p_nh%prog(knew)%vn(1:nlen,jk,jb)=p_nh%prog(know)%vn(1:nlen,jk,jb)&
        &                     +dtime*(p_nh%diag%ddt_vn_vort(1:nlen,jk,jb)&
        &                            +p_nh%diag%ddt_vn_phy (1:nlen,jk,jb)&
        &                            -z_ddxn_emech_e       (1:nlen,jk,jb)&
        &                            -p_nh%diag%horpgrad   (1:nlen,jk,jb)&
        & +p_nh%metrics%ddxn_z_full(1:nlen,jk,jb)*z_grad_corr_e_k(1:nlen,jk,jb))
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_E,p_patch,p_nh%prog(knew)%vn)

    ! Lower bounadry condition for w: Solve implicit equation
!$OMP PARALLEL WORKSHARE
    z_vn_surf(:,:)   = p_nh%prog(knew)%vn(:,nlev,:)
    z_gmres_rhs(:,:) = p_nh%prog(knew)%vn(:,nlev,:)
!$OMP END PARALLEL WORKSHARE
    CALL gmres(z_vn_surf,                     &! x.Input is the first guess
               lhs_lower_boundary,            &! sbr. calculating l.h.s.
               p_patch,                       &! used for calculating l.h.s.
               p_int,                         &! interpolation state
               p_nh%metrics,                  &! metrics fields
               nblks_e,                       &! number of blocks
               npromz_e,                      &! length of last block
               z_1o2,                         &! used for calculating l.h.s.
               z_gmres_rhs,                   &! right hand side as input
               gmres_rtol_nh,                 &! relative tolerance
               .FALSE.,                       &! NOT absolute tolerance
               nmax_iter,                     &! max. # of iterations to do
               lmaxiter,                      &! out: .true. = not converged
               niter,                         &! out: # of iterations done
               z_residual                     &! out: the residual (array)
               )
    IF (lmaxiter) THEN
      CALL finish('GMRES solver: ','NOT YET CONVERGED !!')
    ENDIF
!$OMP PARALLEL
!$OMP WORKSHARE
    p_nh%prog(knew)%vn(:,nlev,:)=z_vn_surf(:,:)
!$OMP END WORKSHARE
!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      ! This becomes averaged the lower boundary condition for w
      z_jn_vn_e(1:nlen,nlevp1,jb) = p_nh%prog(knew)%vn(1:nlen,nlev,jb) &
      &                     * p_nh%metrics%ddxn_z_full(1:nlen,nlev,jb) &
      &                   * p_nh%metrics%ddqz_z_full_e(1:nlen,nlev,jb)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ! prog(knew)%w(:,nlevp1,:)
    CALL edges2cells_scalar(z_jn_vn_e,p_patch,p_int%e_inn_c,&
         &                  p_nh%prog(knew)%w,nlevp1,nlevp1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1,nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      ! The ddqz_z_full is needed here as vn is defined at full levels (see above)
      p_nh%prog(knew)%w(1:nlen,nlevp1,jb)=p_nh%prog(knew)%w(1:nlen,nlevp1,jb) &
      &                           /p_nh%metrics%ddqz_z_full(1:nlen,nlev,jb)
    ENDDO
!$OMP END DO

    ! 4) finalize horizontal part
    !============================
!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      ! Fluxes at edges
      IF (l_predictor) THEN
        DO jk = 1,nlev
          p_nh%diag%mass_fl_e(1:nlen,jk,jb) = p_nh%diag%rho_e(1:nlen,jk,jb) &
          &                               *p_nh%prog(knew)%vn(1:nlen,jk,jb) &
          &                       *p_nh%metrics%ddqz_z_full_e(1:nlen,jk,jb)
        ENDDO
      ELSE
        DO jk = 1,nlev
          p_nh%diag%mass_fl_e(1:nlen,jk,jb)   =   &
          & 0.5_wp*(p_nh%diag%rho_e(1:nlen,jk,jb)+p_nh%diag%rho_star_e(1:nlen,jk,jb))&
          &        *p_nh%prog(knew)%vn(1:nlen,jk,jb) &
          &*p_nh%metrics%ddqz_z_full_e(1:nlen,jk,jb)
        ENDDO
      ENDIF
      ! metric correction for the divergence
      DO jk = nflat+1,nlev
        z_jn_vn_e(1:nlen,jk,jb)=  p_nh%diag%mass_fl_e(1:nlen,jk,jb) &
        &                   *p_nh%metrics%ddxn_z_full(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL edges2cells_scalar(z_jn_vn_e,p_patch,p_int%e_inn_c,&
                            z_contr_corr_c,nflat+1,nlev)

    IF (.NOT. l_predictor) THEN
      ! edge value of theta
      CALL advection_edges (p_nh%diag%theta_v_ave,p_nh%prog(knew)%vn,&
      &                     p_patch,p_nh%metrics,p_int,z_theta_v_e)
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = 1,nlev
        z_theta_v_fl_e(1:nlen,jk,jb)= p_nh%diag%mass_fl_e(1:nlen,jk,jb) &
        &                                    *z_theta_v_e(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! horizontal divergences
    CALL div(p_nh%diag%mass_fl_e   , p_patch, p_int, z_mass_fl_div   )
    CALL div(z_theta_v_fl_e, p_patch, p_int, z_theta_v_fl_div)

    ! 5) Vertical solution (5band matrix)
    !====================================

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,z_w_expl,z_rho_expl,z_exner_expl,&
!$OMP    z_rho_expl_ic,z_kinw_expl_ic,z_kinw_expl_mc,&
!$OMP    z_alpha,z_beta,z_gamma,z_delta,z_epsil,z_gaep,z_gaepd,z_gaepm,&
!$OMP    z_a,z_b,z_c,z_d,z_e,z_r,z_q,z_m,z_k,z_l,z_u,z_v) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=1,nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF

      ! contravariant vertical velocity times density for explicit part
      z_contr_corr_c(1:nlen,nflat,jb) = 0.0_wp
      DO jk = nflat+1,nlev
        z_contr_corr_c(1:nlen,jk,jb) = z_contr_corr_c(1:nlen,jk,jb) &
        &                   /p_nh%metrics%ddqz_z_full(1:nlen,jk,jb)
      ENDDO
      p_nh%diag%w_concorr_c(1:nlen,1:nflat,jb)= 0.0_wp
      DO jk = nflat+1,nlev
        p_nh%diag%w_concorr_c(1:nlen,jk,jb)= &
        & 0.5_wp*(z_contr_corr_c(1:nlen,jk,jb)+z_contr_corr_c(1:nlen,jk-1,jb))
      ENDDO
      p_nh%diag%w_concorr_c(1:nlen,nlevp1,jb)= 0.0_wp

      ! beta, explicit parts of rho and exner
      DO jk = 1, nlev
        z_beta(1:nlen,jk)=   dtime*rd*p_nh%prog(know)%exner(1:nlen,jk  ,jb) &
        &                  /(cvd*p_nh%prog(know)%rhotheta_v(1:nlen,jk  ,jb) &
        &                         *p_nh%metrics%ddqz_z_full(1:nlen,jk  ,jb))
        z_rho_expl(1:nlen,jk)=          p_nh%prog(know)%rho(1:nlen,jk  ,jb) &
        &                   -dtime/p_nh%metrics%ddqz_z_full(1:nlen,jk  ,jb) &
        &                                   *(z_mass_fl_div(1:nlen,jk  ,jb) &
        &                            -p_nh%diag%w_concorr_c(1:nlen,jk  ,jb) &
        &                            +p_nh%diag%w_concorr_c(1:nlen,jk+1,jb))
        z_exner_expl(1:nlen,jk)=      p_nh%prog(know)%exner(1:nlen,jk  ,jb) &
        &                                           -z_beta(1:nlen,jk     ) &
        &                                *(z_theta_v_fl_div(1:nlen,jk  ,jb) &
        &-z_theta_v_l(1:nlen,jk  ,jb)*p_nh%diag%w_concorr_c(1:nlen,jk  ,jb) &
        &+z_theta_v_l(1:nlen,jk+1,jb)*p_nh%diag%w_concorr_c(1:nlen,jk+1,jb))&
        &                   + dtime*p_nh%diag%ddt_exner_phy(1:nlen,jk  ,jb)
      ENDDO
      ! to be compatible with other tracer advection (until found a solution for that)
      DO jk = nflat+1,nlev
        p_nh%diag%w_concorr_c(1:nlen,jk,jb)= &
        & p_nh%diag%w_concorr_c(1:nlen,jk,jb)/p_nh%diag%rho_ic(1:nlen,jk,jb)
      ENDDO

      ! upper boundary condition for w
      p_nh%prog(knew)%w(1:nlen,1,jb) = 0.0_wp
      z_kinw_expl_ic(1:nlen,1)=0.0_wp
      DO jk = 2, nlev
        ! explicit part for w
        z_w_expl(1:nlen,jk) = p_nh%prog(know)%w(1:nlen,jk  ,jb) &
        &             + dtime*( -(z_emech_kin2d(1:nlen,jk-1,jb) &
        &                        -z_emech_kin2d(1:nlen,jk  ,jb))&
        &            / p_nh%metrics%ddqz_z_half(1:nlen,jk  ,jb) &
        &                 +p_nh%diag%ddt_w_vort(1:nlen,jk  ,jb) &
        &                  +p_nh%diag%ddt_w_phy(1:nlen,jk  ,jb) &
        &                  -rd*z_th_ddz_exner_c(1:nlen,jk  ,jb))
        ! density at interface levels
        z_rho_expl_ic(1:nlen,jk)  = 0.5_wp*(z_rho_expl(1:nlen,jk-1) &
        &                                 + z_rho_expl(1:nlen,jk  ))
        z_kinw_expl_ic(1:nlen,jk) = p_nh%prog(know)%w(1:nlen,jk,jb)**2 &
        & * z_rho_expl_ic(1:nlen,jk) * 0.5_wp / p_nh%diag%rho_ic(1:nlen,jk,jb)
      ENDDO
      ! lower boundary for w is already given 
      ! (note: the minus sign is correct here, do not doubt!)
      z_kinw_expl_ic(1:nlen,nlevp1) = - 0.5_wp* p_nh%prog(know)%w(1:nlen,nlevp1,jb) &
      &                                       * p_nh%prog(knew)%w(1:nlen,nlevp1,jb)

      ! kinetic explicit part at full level
      DO jk = 1,nlev
        z_kinw_expl_mc(1:nlen,jk) = 0.5_wp/p_nh%metrics%ddqz_z_full(1:nlen,jk,jb) &
        & *(z_kinw_expl_ic(1:nlen,jk  )*p_nh%metrics%ddqz_z_half(1:nlen,jk  ,jb) &
        &  +z_kinw_expl_ic(1:nlen,jk+1)*p_nh%metrics%ddqz_z_half(1:nlen,jk+1,jb))
      ENDDO

      ! right hand side, and greek variables
      z_gaep(1:nlen,1     ) = 0.0_wp
      DO jk = 2, nlev
        z_delta(1:nlen,jk) = 0.5_wp*dtime*p_nh%diag%rho_ic(1:nlen,jk,jb) &
        &                        /p_nh%metrics%ddqz_z_half(1:nlen,jk,jb)
        z_alpha(1:nlen,jk) = z_delta(1:nlen,jk)*cvd*z_theta_v_l(1:nlen,jk,jb)
        z_gamma(1:nlen,jk) = 0.25_wp*dtime*p_nh%prog(know)%w(1:nlen,jk,jb)
        z_epsil(1:nlen,jk) = 0.5_wp*p_nh%prog(know)%w(1:nlen,jk,jb) &
        & * p_nh%metrics%ddqz_z_half(1:nlen,jk,jb)/p_nh%diag%rho_ic(1:nlen,jk,jb)

        z_gaep (1:nlen,jk) = z_epsil(1:nlen,jk)*z_gamma(1:nlen,jk)
        z_gaepm(1:nlen,jk) = z_gaep(1:nlen,jk)*p_nh%metrics%mult_1_o_dz(1:nlen,jk,jb)
        z_gaepd(1:nlen,jk) = z_gaep(1:nlen,jk)*p_nh%metrics%diff_1_o_dz(1:nlen,jk,jb)

        z_r(1:nlen,jk) =   0.5_wp*p_nh%diag%rho_ic(1:nlen,jk,jb)*z_w_expl(1:nlen,jk) &
        &           +0.5_wp*z_rho_expl_ic(1:nlen,jk)*p_nh%prog(know)%w(1:nlen,jk,jb) &
        &    -z_alpha(1:nlen,jk)*(z_exner_expl(1:nlen,jk-1)-z_exner_expl(1:nlen,jk)) &
        &+z_delta(1:nlen,jk)*(z_kinw_expl_mc(1:nlen,jk-1)-z_kinw_expl_mc(1:nlen,jk))
      ENDDO
      z_gaep(1:nlen,nlevp1) = 0.0_wp

      ! Setup matrix coefficients
      DO jk = 4, nlev
        z_a(1:nlen,jk)= z_delta(1:nlen,jk)*z_gaepm(1:nlen,jk-1)
      ENDDO
      DO jk = 2, nlev-2
        z_e(1:nlen,jk)= z_delta(1:nlen,jk)*z_gaepm(1:nlen,jk+1)
      ENDDO
      DO jk = 3, nlev
        z_b(1:nlen,jk)= -z_alpha(1:nlen,jk)*z_beta(1:nlen,jk-1)*z_theta_v_l(1:nlen,jk-1,jb) &
        & +(z_gamma(1:nlen,jk)+z_delta(1:nlen,jk)*(z_epsil(1:nlen,jk-1) &
        & -z_gaepd(1:nlen,jk-1)+z_gaepd(1:nlen,jk)))/p_nh%metrics%ddqz_z_full(1:nlen,jk-1,jb)
      ENDDO
      DO jk = 2, nlev-1
        z_d(1:nlen,jk)= -z_alpha(1:nlen,jk)*z_beta(1:nlen,jk)*z_theta_v_l(1:nlen,jk+1,jb) &
        & -(z_gamma(1:nlen,jk)+z_delta(1:nlen,jk)*(z_epsil(1:nlen,jk+1) &
        & +z_gaepd(1:nlen,jk)-z_gaepd(1:nlen,jk+1)))/p_nh%metrics%ddqz_z_full(1:nlen,jk,jb)
      ENDDO
      DO jk = 2, nlev
        z_c(1:nlen,jk)=1.0_wp &
        & +z_alpha(1:nlen,jk)*z_theta_v_l(1:nlen,jk,jb)*(z_beta(1:nlen,jk-1)+z_beta(1:nlen,jk)) &
        & -z_gamma(1:nlen,jk)*p_nh%metrics%diff_1_o_dz(1:nlen,jk,jb) &
        & +z_delta(1:nlen,jk)*(2.0_wp*z_gaepm(1:nlen,jk) &
        &   +z_epsil(1:nlen,jk)*p_nh%metrics%diff_1_o_dz(1:nlen,jk,jb) &
        & -(z_gaep(1:nlen,jk-1)+z_gaep(1:nlen,jk))/(p_nh%metrics%ddqz_z_full(1:nlen,jk-1,jb)**2)&
        & -(z_gaep(1:nlen,jk)+z_gaep(1:nlen,jk+1))/(p_nh%metrics%ddqz_z_full(1:nlen,jk  ,jb)**2))
      ENDDO

      ! Solve pentadiagonal matrix for the mass flux 
      !---------------------------------------------

      z_u(1:nlen,2) = z_c(1:nlen,2)
      z_v(1:nlen,2) = z_d(1:nlen,2)
      z_q(1:nlen,2) = z_r(1:nlen,2)

      z_l(1:nlen)   = z_b(1:nlen,3)/z_u(1:nlen,2)
      z_u(1:nlen,3) = z_c(1:nlen,3)-z_l(1:nlen)*z_v(1:nlen,2)
      z_v(1:nlen,3) = z_d(1:nlen,3)-z_l(1:nlen)*z_e(1:nlen,2)
      z_q(1:nlen,3) = z_r(1:nlen,3)-z_l(1:nlen)*z_q(1:nlen,2)

      DO jk = 4, nlev-1
        z_k(1:nlen)   = z_a(1:nlen,jk)/z_u(1:nlen,jk-2)
        z_l(1:nlen)   =(z_b(1:nlen,jk)-z_k(1:nlen)*z_v(1:nlen,jk-2))/z_u(1:nlen,jk-1)
        z_u(1:nlen,jk)= z_c(1:nlen,jk)-z_k(1:nlen)*z_e(1:nlen,jk-2)-z_l(1:nlen)*z_v(1:nlen,jk-1)
        z_v(1:nlen,jk)= z_d(1:nlen,jk)                             -z_l(1:nlen)*z_e(1:nlen,jk-1)
        z_q(1:nlen,jk)= z_r(1:nlen,jk)-z_k(1:nlen)*z_q(1:nlen,jk-2)-z_l(1:nlen)*z_q(1:nlen,jk-1)
      ENDDO
      z_k(1:nlen)     = z_a(1:nlen,nlev)/z_u(1:nlen,nlev-2)
      z_l(1:nlen)     =(z_b(1:nlen,nlev)-z_k(1:nlen)*z_v(1:nlen,nlev-2))/z_u(1:nlen,nlev-1)
      z_u(1:nlen,nlev)= z_c(1:nlen,nlev)-z_k(1:nlen)*z_e(1:nlen,nlev-2)&
      &                 -z_l(1:nlen)*z_v(1:nlen,nlev-1)
      z_q(1:nlen,nlev)= z_r(1:nlen,nlev)-z_k(1:nlen)*z_q(1:nlen,nlev-2)&
      &                 -z_l(1:nlen)*z_q(1:nlen,nlev-1)

      z_m(1:nlen,nlevp1) = 0.0_wp
      z_m(1:nlen,nlev) = z_q(1:nlen,nlev)/z_u(1:nlen,nlev)
      z_m(1:nlen,nlev-1) = (z_q(1:nlen,nlev-1)-z_v(1:nlen,nlev-1)*z_m(1:nlen,nlev)) &
      &                    /z_u(1:nlen,nlev-1)
      DO jk = nlev-2,2,-1
        z_m(1:nlen,jk) = (z_q(1:nlen,jk)-z_v(1:nlen,jk)*z_m(1:nlen,jk+1)&
        & -z_e(1:nlen,jk)*z_m(1:nlen,jk+2))/z_u(1:nlen,jk)
      ENDDO 
      z_m(1:nlen,1) = 0.0_wp

      ! Results
      DO jk = 1, nlev

        ! density
        p_nh%prog(knew)%rho(1:nlen,jk,jb)=  z_rho_expl(1:nlen,jk)&
        &     -dtime/p_nh%metrics%ddqz_z_full(1:nlen,jk,jb)&
        &     *(z_m(1:nlen,jk)-z_m(1:nlen,jk+1))

        ! exner
        p_nh%prog(knew)%exner(1:nlen,jk,jb)=z_exner_expl(1:nlen,jk)        &
        & -z_beta(1:nlen,jk)*(z_theta_v_l(1:nlen,jk  ,jb)*z_m(1:nlen,jk  ) &
        &                    -z_theta_v_l(1:nlen,jk+1,jb)*z_m(1:nlen,jk+1))

        ! rho*theta
        p_nh%prog(knew)%rhotheta_v(1:nlen,jk,jb)= &
        &    p_nh%prog(know)%rhotheta_v(1:nlen,jk,jb)&
        &      *((p_nh%prog(knew)%exner(1:nlen,jk,jb)&
        &        /p_nh%prog(know)%exner(1:nlen,jk,jb)-1.0_wp)*cvd_o_rd+1.0_wp)

        ! theta
        p_nh%prog(knew)%theta_v(1:nlen,jk,jb) = &
        &       p_nh%prog(knew)%rhotheta_v(1:nlen,jk,jb)&
        &      /p_nh%prog(knew)%rho       (1:nlen,jk,jb)

      ENDDO

      DO jk = 2, nlev
        ! new density at interface levels
        z_l(1:nlen)  = 0.5_wp*(p_nh%prog(knew)%rho(1:nlen,jk-1,jb) &
                              +p_nh%prog(knew)%rho(1:nlen,jk  ,jb)) 
        p_nh%prog(knew)%w(1:nlen,jk,jb) = (2.0_wp*z_m(1:nlen,jk) &
        & -p_nh%prog(know)%w(1:nlen,jk,jb)*z_l(1:nlen))/p_nh%diag%rho_ic(1:nlen,jk,jb)

      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array_mult(SYNC_C,p_patch,5,p_nh%prog(knew)%theta_v,       &
                               p_nh%prog(knew)%rhotheta_v,p_nh%prog(knew)%rho, &
                               p_nh%prog(knew)%exner,p_nh%prog(knew)%w)

  END SUBROUTINE divergent_modes_5band
  !=============================================================================
  !=============================================================================
  !>
  !! lhs_lower_boundary
  !!
  !! This function provides the lhs for the computation of the lower boundary
  !! condition for the horizontal normal wind component. It is called by
  !! the gmres solver.
  !!
  !! @par Revision History
  !! Inital release by Almut Gassmann (2009-05-16)
  !!
  FUNCTION lhs_lower_boundary (p_vn,p_patch,p_int,p_metrics,p_coeff)  RESULT (p_lhs)

    ! input variables
    REAL(wp),INTENT(IN)         :: p_vn(:,:)
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    TYPE(t_int_state), INTENT(IN) :: p_int
    TYPE(t_nh_metrics), INTENT(IN):: p_metrics
    REAL(wp),INTENT(IN)         :: p_coeff

    ! output variables
    REAL(wp) :: p_lhs    (SIZE(p_vn,1),SIZE(p_vn,2))

    REAL(wp) :: z_jn_vn_e(nproma,1,p_patch%nblks_e), &
                w_knew_e (nproma,1,p_patch%nblks_e), &
                w_knew_c (nproma,1,p_patch%nblks_c)

    INTEGER :: nblks_e, npromz_e, nblks_c, npromz_c, nlen, jb
    INTEGER :: nlev, nlevp1            !< number of vertical levels

    !--------------------------------------------------------

    p_lhs=0.0_wp

    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e
    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c

    ! number of vertical levels
    nlev = p_patch%nlev
    nlevp1 = p_patch%nlevp1

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      z_jn_vn_e(1:nlen,1,jb) = p_vn(1:nlen,jb) &
      & *p_metrics%ddxn_z_full(1:nlen,nlev,jb) &
      & *p_metrics%ddqz_z_full_e(1:nlen,nlev,jb)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL edges2cells_scalar(z_jn_vn_e,p_patch,p_int%e_inn_c,w_knew_c,1,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=1,nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      w_knew_c(1:nlen,1,jb)=w_knew_c(1:nlen,1     ,jb)&
      &       /p_metrics%ddqz_z_full(1:nlen,nlev  ,jb)&
      &       /p_metrics%ddqz_z_full(1:nlen,nlev  ,jb)&
      &       *p_metrics%ddqz_z_half(1:nlen,nlevp1,jb)&
      &       *p_coeff
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL cells2edges_scalar(w_knew_c,p_patch,p_int%c_lin_e,w_knew_e,1,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      p_lhs(1:nlen,jb) = p_vn(1:nlen,jb) + w_knew_e(1:nlen,1   ,jb) &
      &                     * p_metrics%ddxn_z_full(1:nlen,nlev,jb)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_E,p_patch,p_lhs)

  END FUNCTION lhs_lower_boundary

  !----------------------------------------------------------------------------
  !>
  !! impl_vert_adv_theta 
  !!
  !! Implicit vertical advection for theta_v
  !! 3 diagonal matrix is solved. (Durran,1999:pp.440-441) 
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2011-01-11)
  !!
  SUBROUTINE impl_vert_adv_theta (p_th,p_mflx_w_con,p_rho,p_patch,p_metrics,p_ti)

    ! Passed variables
    REAL(wp), INTENT(in)              :: p_th(:,:,:)  !< potential temperature
    REAL(wp), INTENT(in)              :: p_mflx_w_con(:,:,:) &
    & ;  !< vertical contravariant mass flux 
    REAL(wp), INTENT(in)              :: p_rho(:,:,:) !< density
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch      !< patch
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics    !< metrics
    REAL(wp), INTENT(inout)           :: p_ti(:,:,:)  &
    &  ;  !< implicit potential temperature: (theta_v(nnow)+theta_v(nnew))/2

    ! Local variables
    REAL(wp):: z_a(nproma,p_patch%nlev), &
      &        z_b(nproma,p_patch%nlev), &
      &        z_c(nproma,p_patch%nlev)  !< matrix coefficients
    REAL(wp):: z_fac(nproma), z_q(nproma,p_patch%nlev)  !< auxiliary vectors
    INTEGER :: nblks_c, npromz_c, jb, jk, nlen
    INTEGER :: nlev            !< number of full levels
    !--------------------------------------------------------------------------

    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c

    ! number of vertical levels
    nlev = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, z_a, z_b, z_c, z_q, z_fac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1,nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = 1, nlev
        p_ti(1:nlen,jk,jb) = p_th(1:nlen,jk,jb)
        z_fac(1:nlen)  = 0.25_wp*dtime/p_metrics%ddqz_z_full(1:nlen,jk,jb) &
        &                /p_rho(1:nlen,jk,jb)
        z_b(1:nlen,jk) = 1.0_wp-z_fac(1:nlen)*(p_mflx_w_con(1:nlen,jk  ,jb)&
        &                                     -p_mflx_w_con(1:nlen,jk+1,jb))
        z_a(1:nlen,jk) =  z_fac(1:nlen)*p_mflx_w_con(1:nlen,jk  ,jb)
        z_c(1:nlen,jk) = -z_fac(1:nlen)*p_mflx_w_con(1:nlen,jk+1,jb)
      ENDDO
      z_q (1:nlen,1   ) = -z_c(1:nlen,1)/z_b(1:nlen,1)
      p_ti(1:nlen,1,jb) = p_ti(1:nlen,1,jb)/z_b(1:nlen,1)
      DO jk = 2, nlev
        z_fac(1:nlen)      = 1.0_wp/(z_b(1:nlen,jk) + z_a(1:nlen,jk)*z_q(1:nlen,jk-1))
        z_q(1:nlen,jk)     = -z_c(1:nlen,jk)*z_fac(1:nlen)
        p_ti(1:nlen,jk,jb) = (p_ti(1:nlen,jk,jb)  &
        &                    -z_a(1:nlen,jk)*p_ti(1:nlen,jk-1,jb))*z_fac(1:nlen)
      ENDDO
      DO jk=nlev-1,1,-1
        p_ti(1:nlen,jk,jb) = p_ti(1:nlen,jk,jb)+z_q(1:nlen,jk)*p_ti(1:nlen,jk+1,jb)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE impl_vert_adv_theta
  !----------------------------------------------------------------------------
  !>
  !! impl_vert_adv_theta_5diag
  !!
  !! Implicit vertical advection for theta_v
  !! 5 diagonal matrix is solved.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2011-05-06)
  !!
  SUBROUTINE impl_vert_adv_theta_5diag (p_th,p_mflx_w_con,p_rho,p_patch,p_metrics,p_ti)

    ! Passed variables
    REAL(wp), INTENT(in)              :: p_th(:,:,:)  !< potential temperature
    REAL(wp), INTENT(in)              :: p_mflx_w_con(:,:,:) &
    & ;  !< vertical contravariant mass flux 
    REAL(wp), INTENT(in)              :: p_rho(:,:,:) !< density
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch      !< patch
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics    !< metrics
    REAL(wp), INTENT(inout)           :: p_ti(:,:,:)  &
    &  ;  !< implicit potential temperature: (theta_v(nnow)+theta_v(nnew))/2

    ! Local variables
    REAL(wp):: z_a(nproma,p_patch%nlev), &
      &        z_b(nproma,p_patch%nlev), &
      &        z_c(nproma,p_patch%nlev), &
      &        z_d(nproma,p_patch%nlev), &
      &        z_e(nproma,p_patch%nlev), &
      &        z_m(nproma,p_patch%nlev), &
      &        z_q(nproma,p_patch%nlev), &
      &        z_r(nproma,p_patch%nlev), &
      &        z_u(nproma,p_patch%nlev), &
      &        z_v(nproma,p_patch%nlev), &
      &        z_k(nproma), &
      &        z_l(nproma)
    REAL(wp):: z_sig_p(nproma,p_patch%nlevp1), &
      &        z_sig_m(nproma,p_patch%nlevp1)
    REAL(wp):: z_fac(nproma,p_patch%nlev)  !< auxiliary vectors
    INTEGER :: nblks_c, npromz_c, jb, jk, nlen
    INTEGER :: nlev,nlevp1
    !--------------------------------------------------------------------------

    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c
    nlev      = p_patch%nlev
    nlevp1    = p_patch%nlevp1

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, z_a, z_b, z_c, z_d, z_e, z_q, z_r, z_fac, z_sig_p, z_sig_m, &
!$OMP            z_k, z_l, z_u, z_v, z_m) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1,nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF

      ! setup matrix coeffs and rhs
      z_sig_p(1:nlen,1)=0.0_wp
      z_sig_m(1:nlen,1)=0.0_wp
      DO jk = 2, nlev
        z_sig_p(1:nlen,jk) = (0.5_wp+upstr_beta*SIGN(0.5_wp,p_mflx_w_con(1:nlen,jk,jb)))/6.0_wp
        z_sig_m(1:nlen,jk) = (0.5_wp-upstr_beta*SIGN(0.5_wp,p_mflx_w_con(1:nlen,jk,jb)))/6.0_wp
      ENDDO
      z_sig_p(1:nlen,nlevp1)= 0.0_wp ! no flux
      z_sig_m(1:nlen,nlevp1)= 0.0_wp
      z_sig_m(1:nlen,2)     = 0.0_wp ! higher order only for upward flux allowed
      z_sig_p(1:nlen,nlev)  = 0.0_wp ! higher order only for downward flux allowed
      DO jk = 1, nlev
        z_r(1:nlen,jk)  = p_th(1:nlen,jk,jb) ! rhs
        z_fac(1:nlen,jk)= 0.5_wp*dtime/p_metrics%ddqz_z_full(1:nlen,jk,jb)/p_rho(1:nlen,jk,jb)
      ENDDO
      DO jk = 3, nlev
        z_a(1:nlen,jk) = -z_fac(1:nlen,jk)*p_mflx_w_con(1:nlen,jk  ,jb)*z_sig_m(1:nlen,jk  ) &
        & * p_metrics%ddqz_z_half(1:nlen,jk  ,jb) / p_metrics%ddqz_z_full(1:nlen,jk  ,jb)
      ENDDO
      DO jk = 1, nlev-2
        z_e(1:nlen,jk) =  z_fac(1:nlen,jk)*p_mflx_w_con(1:nlen,jk+1,jb)*z_sig_p(1:nlen,jk+1) &
        & * p_metrics%ddqz_z_half(1:nlen,jk+1,jb) / p_metrics%ddqz_z_full(1:nlen,jk+1,jb)
      ENDDO
      DO jk = 2, nlev
        z_b(1:nlen,jk) = z_fac(1:nlen,jk)*( &
        &  p_mflx_w_con(1:nlen,jk  ,jb)*(0.5_wp+2.0_wp*z_sig_m(1:nlen,jk  )-z_sig_p(1:nlen,jk  ) &
        & *p_metrics%ddqz_z_half(1:nlen,jk+1,jb)/p_metrics%ddqz_z_full(1:nlen,jk,jb)) &
        & +p_mflx_w_con(1:nlen,jk+1,jb)*z_sig_m(1:nlen,jk+1) &
        & *p_metrics%ddqz_z_half(1:nlen,jk+1,jb)/p_metrics%ddqz_z_full(1:nlen,jk,jb))
      ENDDO
      DO jk = 1, nlev-1
        z_d(1:nlen,jk) = z_fac(1:nlen,jk)*( &
        & -p_mflx_w_con(1:nlen,jk  ,jb)*z_sig_p(1:nlen,jk  ) &
        & *p_metrics%ddqz_z_half(1:nlen,jk  ,jb)/p_metrics%ddqz_z_full(1:nlen,jk,jb) & 
        & -p_mflx_w_con(1:nlen,jk+1,jb)*(0.5_wp+2.0_wp*z_sig_p(1:nlen,jk+1)-z_sig_m(1:nlen,jk+1) &
        & *p_metrics%ddqz_z_half(1:nlen,jk  ,jb)/p_metrics%ddqz_z_full(1:nlen,jk,jb)))
      ENDDO
      z_c(1:nlen,1) = 1.0_wp+z_fac(1:nlen,1)*( &
      & -p_mflx_w_con(1:nlen,2,jb)*(-0.5_wp-z_sig_p(1:nlen,2) &
      & *p_metrics%ddqz_z_half(1:nlen,3,jb)/p_metrics%ddqz_z_full(1:nlen,2,jb)))
      DO jk = 2, nlev-1
        z_c(1:nlen,jk) = 1.0_wp+z_fac(1:nlen,jk)*( &
        &  p_mflx_w_con(1:nlen,jk  ,jb)*(-0.5_wp+2.0_wp*z_sig_p(1:nlen,jk  )-z_sig_m(1:nlen,jk  )&
        & *p_metrics%ddqz_z_half(1:nlen,jk-1,jb)/p_metrics%ddqz_z_full(1:nlen,jk-1,jb)) & 
        & -p_mflx_w_con(1:nlen,jk+1,jb)*(-0.5_wp+2.0_wp*z_sig_m(1:nlen,jk+1)-z_sig_p(1:nlen,jk+1)&
        & *p_metrics%ddqz_z_half(1:nlen,jk+2,jb)/p_metrics%ddqz_z_full(1:nlen,jk+1,jb))) 
      ENDDO
      z_c(1:nlen,nlev) = 1.0_wp+z_fac(1:nlen,nlev)*( &
      &  p_mflx_w_con(1:nlen,nlev,jb)*(-0.5_wp-z_sig_m(1:nlen,nlev) &
      & *p_metrics%ddqz_z_half(1:nlen,nlev-1,jb)/p_metrics%ddqz_z_full(1:nlen,nlev-1,jb)))
      
      ! Solve pentadiagonal matrix 
      !---------------------------------------------

      z_u(1:nlen,1) = z_c(1:nlen,1)
      z_v(1:nlen,1) = z_d(1:nlen,1)
      z_q(1:nlen,1) = z_r(1:nlen,1)

      z_l(1:nlen)   = z_b(1:nlen,2)/z_u(1:nlen,1)
      z_u(1:nlen,2) = z_c(1:nlen,2)-z_l(1:nlen)*z_v(1:nlen,1)
      z_v(1:nlen,2) = z_d(1:nlen,2)-z_l(1:nlen)*z_e(1:nlen,1)
      z_q(1:nlen,2) = z_r(1:nlen,2)-z_l(1:nlen)*z_q(1:nlen,1)

      DO jk = 3, nlev-1
        z_k(1:nlen)   = z_a(1:nlen,jk)/z_u(1:nlen,jk-2)
        z_l(1:nlen)   =(z_b(1:nlen,jk)-z_k(1:nlen)*z_v(1:nlen,jk-2))/z_u(1:nlen,jk-1)
        z_u(1:nlen,jk)= z_c(1:nlen,jk)-z_k(1:nlen)*z_e(1:nlen,jk-2)-z_l(1:nlen)*z_v(1:nlen,jk-1)
        z_v(1:nlen,jk)= z_d(1:nlen,jk)                             -z_l(1:nlen)*z_e(1:nlen,jk-1)
        z_q(1:nlen,jk)= z_r(1:nlen,jk)-z_k(1:nlen)*z_q(1:nlen,jk-2)-z_l(1:nlen)*z_q(1:nlen,jk-1)
      ENDDO
      z_k(1:nlen)     = z_a(1:nlen,nlev)/z_u(1:nlen,nlev-2)
      z_l(1:nlen)     =(z_b(1:nlen,nlev)-z_k(1:nlen)*z_v(1:nlen,nlev-2))/z_u(1:nlen,nlev-1)
      z_u(1:nlen,nlev)= z_c(1:nlen,nlev)-z_k(1:nlen)*z_e(1:nlen,nlev-2)&
      &                 -z_l(1:nlen)*z_v(1:nlen,nlev-1)
      z_q(1:nlen,nlev)= z_r(1:nlen,nlev)-z_k(1:nlen)*z_q(1:nlen,nlev-2)&
      &                 -z_l(1:nlen)*z_q(1:nlen,nlev-1)

      z_m(1:nlen,nlev) = z_q(1:nlen,nlev)/z_u(1:nlen,nlev)
      z_m(1:nlen,nlev-1) = (z_q(1:nlen,nlev-1)-z_v(1:nlen,nlev-1)*z_m(1:nlen,nlev)) &
      &                    /z_u(1:nlen,nlev-1)
      DO jk = nlev-2,1,-1
        z_m(1:nlen,jk) = (z_q(1:nlen,jk)-z_v(1:nlen,jk)*z_m(1:nlen,jk+1)&
        & -z_e(1:nlen,jk)*z_m(1:nlen,jk+2))/z_u(1:nlen,jk)
      ENDDO 

      DO jk = 1, nlev
        p_ti(1:nlen,jk,jb) = z_m(1:nlen,jk)
      ENDDO 
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE impl_vert_adv_theta_5diag
  !----------------------------------------------------------------------------
  !>
  !! advection_edges
  !!
  !! supply the edge value for horizontal advection
  !! (Skamarock and Gassmann, MWR, 2011)
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2011-04-11)
  !!
  SUBROUTINE advection_edges (p_s_c,p_vn,p_patch,p_metrics,p_int,p_s_e)

    ! Passed variables
    REAL(wp), INTENT(in)              :: p_s_c(:,:,:) !< scalar value at centers
    REAL(wp), INTENT(in)              :: p_vn(:,:,:)  !< horizontal velocity
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch      !< patch
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics    !< metrics
    TYPE(t_int_state), INTENT(IN)     :: p_int        !< interpolation state
    REAL(wp), INTENT(inout)           :: p_s_e(:,:,:) !< scalar value at edges

    ! Local variables
    INTEGER :: nblks_c, nblks_e, npromz_c, npromz_e, nlen, jk, jb
    INTEGER :: nlev, nlevp1
    REAL(wp):: z_ddz_s_cl      (nproma,p_patch%nlevp1,p_patch%nblks_c), &
               z_ddz_s_ck      (nproma,p_patch%nlev  ,p_patch%nblks_c), &
               z_ddz_s_ek      (nproma,p_patch%nlev  ,p_patch%nblks_e), &
               z_ddxn_s_e      (nproma,p_patch%nlev  ,p_patch%nblks_e), &
               z_lapl_s_e      (nproma,p_patch%nlev  ,p_patch%nblks_e), &
               z_lapl_terrain_e(nproma,p_patch%nlev  ,p_patch%nblks_e)
    !--------------------------------------------------------------------------

    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c
    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e
    nlev      = p_patch%nlev
    nlevp1    = p_patch%nlevp1

    ! scalar value averaged to edges
    CALL cells2edges_scalar(p_s_c,p_patch,p_int%c_lin_e,p_s_e)

    IF (ltheta_up_hori) THEN
      !include terrain following correction here
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1,nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        !1) vertical gradient of scalar
        DO jk = 2,nlev
          z_ddz_s_cl(1:nlen,jk,jb) = p_s_c(1:nlen,jk-1,jb)-p_s_c(1:nlen,jk,jb) 
        ENDDO
        z_ddz_s_cl(1:nlen,1,jb)      = z_ddz_s_cl(1:nlen,2,jb)
        z_ddz_s_cl(1:nlen,nlevp1,jb) = z_ddz_s_cl(1:nlen,nlev,jb)
        !2) average gradient to main level
        DO jk = 1,nlev
          z_ddz_s_ck(1:nlen,jk,jb)=0.5_wp/p_metrics%ddqz_z_full(1:nlen,jk,jb)&
          &*(z_ddz_s_cl(1:nlen,jk,jb)+z_ddz_s_cl(1:nlen,jk+1,jb))
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      !3) average to edges
      CALL cells2edges_scalar(z_ddz_s_ck,p_patch,p_int%c_lin_e,z_ddz_s_ek)

      IF (lplane) THEN
        ! compute horizontal gradient of potential temperature
        CALL grad_fd_norm(p_s_c,p_patch,z_ddxn_s_e)
        CALL sync_patch_array(SYNC_E,p_patch,z_ddxn_s_e)
        ! compute directional laplace
        CALL grad_dir_edge(p_vn,z_ddxn_s_e,p_patch,p_int,z_lapl_s_e)
        ! terrain following correction
        CALL grad_dir_edge(p_vn,p_metrics%ddxn_z_full,p_patch,p_int,z_lapl_terrain_e)
      ELSE
        CALL directional_laplace(p_vn,p_s_c,p_patch,p_int,upstr_beta,z_lapl_s_e)
        CALL directional_laplace(p_vn,p_metrics%z_mc,p_patch,p_int,upstr_beta,z_lapl_terrain_e)
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO jk = 1, nlev
            ! terrain following correction
            z_lapl_s_e(1:nlen,jk,jb)=z_lapl_s_e(1:nlen,jk,jb) &
            & - z_ddz_s_ek(1:nlen,jk,jb)*z_lapl_terrain_e(1:nlen,jk,jb)    
            ! edge value
            p_s_e(1:nlen,jk,jb) = p_s_e(1:nlen,jk,jb) &
            & - p_patch%edges%dual_edge_length(1:nlen,jb)  &
            & * p_patch%edges%dual_edge_length(1:nlen,jb)  &
            & /6.0_wp * z_lapl_s_e(1:nlen,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF

  END SUBROUTINE advection_edges
  !----------------------------------------------------------------------------

END MODULE mo_divergent_modes

