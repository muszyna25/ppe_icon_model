!>
!! mo_velocity_advection
!!
!! This module contains the subroutine calculating the velocity advection tendencies
!! for the nonhydrostatic dynamical core. Separated from mo_solve_nonhydro in order
!! to speed up compile time
!!
!! @author Guenther Zaengl, DWD
!!
!! @par Revision History
!! Created by Guenther Zaengl, DWD on 2013-09-13)
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

MODULE mo_velocity_advection

  USE mo_kind,                 ONLY: wp, vp
  USE mo_nonhydrostatic_config,ONLY: lextra_diffu, &
                                     lbackward_integr, veladv_offctr
  USE mo_parallel_config,   ONLY: nproma
  USE mo_run_config,        ONLY: lvert_nest, ltestcase
  USE mo_model_domain,      ONLY: t_patch
  USE mo_intp_data_strc,    ONLY: t_int_state
  USE mo_intp,              ONLY: cells2verts_scalar
  USE mo_nonhydro_types,    ONLY: t_nh_state, t_nh_metrics, t_nh_diag, t_nh_prog
  USE mo_math_gradients,    ONLY: grad_green_gauss_cell
  USE mo_math_constants,    ONLY: dbl_eps
  USE mo_math_divrot,       ONLY: rot_vertex
  USE mo_vertical_grid,     ONLY: nrdmax
  USE mo_nh_init_utils,     ONLY: nflatlev
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,    ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int, &
    &                             min_rlcell
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_nh_testcases_nml,  ONLY: nh_test_name
  USE mo_nh_dcmip_gw,       ONLY: fcfugal


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  PUBLIC :: velocity_tendencies

  CONTAINS


  !----------------------------------------------------------------------------
  !>
  !! velocity_tendencies
  !!
  !! Discretization of nonhydrostatic momentum equation similar to hydrostatic core
  !! In particular, the Lamb transformation is applied only to the horizontal
  !! equation of motion, whereas the vertical wind equation is discretized
  !! in advective form
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl (2010-02-03)
  !!
  SUBROUTINE velocity_tendencies (p_prog, p_patch, p_int, p_metrics, p_diag, z_w_concorr_me, z_kin_hor_e, &
                                  z_vt_ie, ntnd, istep, lvn_only, dtime)

    ! Passed variables
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_int_state), TARGET, INTENT(IN):: p_int
    TYPE(t_nh_prog), INTENT(INOUT)       :: p_prog
    TYPE(t_nh_metrics), INTENT(IN)       :: p_metrics
    TYPE(t_nh_diag), INTENT(INOUT)       :: p_diag

    ! Local variables from solve_nh that are passed for efficiency optimization
    REAL(vp), DIMENSION(:,:,:), INTENT(INOUT) :: z_w_concorr_me, z_kin_hor_e, z_vt_ie

    INTEGER, INTENT(IN)  :: ntnd     ! time level of ddt_adv fields used to store tendencies
    INTEGER, INTENT(IN)  :: istep    ! 1: predictor step, 2: corrector step
    LOGICAL, INTENT(IN)  :: lvn_only ! true: compute only vn tendency
    REAL(wp),INTENT(IN)  :: dtime    ! time step

    ! Local variables
    INTEGER :: jb, jk, jc, je
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: i_startblk_2, i_endblk_2, i_startidx_2, i_endidx_2
    INTEGER :: rl_start, rl_end, rl_start_2, rl_end_2
    ! The data type vp (variable precision) is by default the same as wp but reduces
    ! to single precision when the __MIXED_PRECISION cpp flag is set at compile time
    REAL(vp):: z_w_concorr_mc(nproma,p_patch%nlev)
    REAL(vp):: z_w_con_c(nproma,p_patch%nlevp1,p_patch%nblks_c)
    REAL(vp):: z_w_con_c_full(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(vp):: z_ddxn_ekin_e(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(vp):: z_v_grad_w(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp):: z_w_v(nproma,p_patch%nlevp1,p_patch%nblks_v) ! should also become vp

    ! Pointers
    INTEGER, DIMENSION(:,:,:), POINTER   &
#ifdef _CRAYFTN
      , CONTIGUOUS                       &
#endif
      ::                                 &
      icidx, icblk, ieidx, ieblk, iqidx, iqblk, ividx, ivblk, incidx, incblk

    INTEGER  :: nlev, nlevp1          !< number of full and half levels
    ! Local control variable for vertical nesting
    LOGICAL :: l_vert_nested

    INTEGER :: jg

    ! Variables for optional extra diffusion close to the vertical advection stability limit
    REAL(vp) :: cfl_w_limit
    REAL(wp) :: scalfac_exdiff, difcoef

    INTEGER  :: icount(p_patch%nblks_c), iclist(p_patch%nlev*nproma,p_patch%nblks_c), &
                iklist(4*nproma,p_patch%nblks_c), ic, ie, jbe

    !--------------------------------------------------------------------------

    IF ((lvert_nest) .AND. (p_patch%nshift > 0)) THEN  
      l_vert_nested = .TRUE.
    ELSE
      l_vert_nested = .FALSE.
    ENDIF

    !Get patch id
    jg = p_patch%id

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Set pointers to neighbor cells/edges/vertices
    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    incidx => p_patch%cells%neighbor_idx
    incblk => p_patch%cells%neighbor_blk

    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk

    i_nchdom   = MAX(1,p_patch%n_childdom)

    ! Limit on vertical CFL number for applying extra diffusion
    cfl_w_limit = 0.75_wp/dtime   ! this means 75% of the nominal CFL stability limit

    ! Scaling factor for extra diffusion
    scalfac_exdiff = 0.075_wp / ( dtime*(1._wp - cfl_w_limit*dtime) )

    ! Compute w at vertices
    IF (.NOT. lvn_only) CALL cells2verts_scalar(p_prog%w, p_patch, &
      p_int%cells_aw_verts, z_w_v, opt_rlend=min_rlvert_int-1)

    ! Compute vertical vorticity component at vertices
    CALL rot_vertex (p_prog%vn, p_patch, p_int, p_diag%omega_z, opt_rlend=min_rlvert_int-1)


!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)

    IF (istep == 1) THEN ! Computations of velocity-derived quantities that come from solve_nh in istep=2

      rl_start = 5
      rl_end = min_rledge_int - 2

      i_startblk = p_patch%edges%start_blk(rl_start,1)
      i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! RBF reconstruction of tangential wind component
            p_diag%vt(je,jk,jb) = &
              p_int%rbf_vec_coeff_e(1,je,jb) * p_prog%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) + &
              p_int%rbf_vec_coeff_e(2,je,jb) * p_prog%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) + &
              p_int%rbf_vec_coeff_e(3,je,jb) * p_prog%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) + &
              p_int%rbf_vec_coeff_e(4,je,jb) * p_prog%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))
          ENDDO
        ENDDO

        ! Interpolate vn to interface levels and compute horizontal part of kinetic energy on edges
        DO jk = 2, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_diag%vn_ie(je,jk,jb) =                                    &
              p_metrics%wgtfac_e(je,jk,jb)*p_prog%vn(je,jk,jb) +        &
             (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_prog%vn(je,jk-1,jb)
            z_kin_hor_e(je,jk,jb) = 0.5_wp*(p_prog%vn(je,jk,jb)**2 + p_diag%vt(je,jk,jb)**2)
          ENDDO
        ENDDO

        IF (.NOT. lvn_only) THEN ! Interpolate also vt to interface levels
          DO jk = 2, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              z_vt_ie(je,jk,jb) =                                         &
                p_metrics%wgtfac_e(je,jk,jb)*p_diag%vt(je,jk,jb) +        &
               (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_diag%vt(je,jk-1,jb)
            ENDDO
          ENDDO
        ENDIF

        ! Compute contravariant correction for vertical velocity at interface levels
        ! (will be interpolated to cell centers below)
        DO jk = nflatlev(p_patch%id), nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            z_w_concorr_me(je,jk,jb) =                              &
              p_prog%vn(je,jk,jb)*p_metrics%ddxn_z_full(je,jk,jb) + &
              p_diag%vt(je,jk,jb)*p_metrics%ddxt_z_full(je,jk,jb)
          ENDDO
        ENDDO

        IF (.NOT. l_vert_nested) THEN
          ! Top and bottom levels
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            ! Quadratic extrapolation at the top turned out to cause numerical instability in pathological cases,
            ! thus we use a no-gradient condition in the upper half layer
            p_diag%vn_ie(je,1,jb) = p_prog%vn(je,1,jb)
            ! vt_ie(jk=1) is actually unused, but we need it for convenience of implementation
            z_vt_ie(je,1,jb) = p_diag%vt(je,1,jb)
            !
            z_kin_hor_e(je,1,jb) = 0.5_wp*(p_prog%vn(je,1,jb)**2 + p_diag%vt(je,1,jb)**2)
            p_diag%vn_ie(je,nlevp1,jb) =                           &
              p_metrics%wgtfacq_e(je,1,jb)*p_prog%vn(je,nlev,jb) +   &
              p_metrics%wgtfacq_e(je,2,jb)*p_prog%vn(je,nlev-1,jb) + &
              p_metrics%wgtfacq_e(je,3,jb)*p_prog%vn(je,nlev-2,jb)
          ENDDO
        ELSE
          ! vn_ie(jk=1) is extrapolated using parent domain information in this case
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_diag%vn_ie(je,1,jb) = p_diag%vn_ie(je,2,jb) + p_diag%dvn_ie_ubc(je,jb)
            ! vt_ie(jk=1) is actually unused, but we need it for convenience of implementation
            z_vt_ie(je,1,jb) = p_diag%vt(je,1,jb)
            !
            z_kin_hor_e(je,1,jb) = 0.5_wp*(p_prog%vn(je,1,jb)**2 + p_diag%vt(je,1,jb)**2)
            p_diag%vn_ie(je,nlevp1,jb) =                           &
              p_metrics%wgtfacq_e(je,1,jb)*p_prog%vn(je,nlev,jb) +   &
              p_metrics%wgtfacq_e(je,2,jb)*p_prog%vn(je,nlev-1,jb) + &
              p_metrics%wgtfacq_e(je,3,jb)*p_prog%vn(je,nlev-2,jb)
          ENDDO
        ENDIF

      ENDDO
!$OMP END DO
    ENDIF ! istep = 1

    rl_start = 3
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, z_w_concorr_mc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Interpolate horizontal kinetic energy to cell centers
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
!DIR$ IVDEP
        DO jk = 1, nlev
#else
!CDIR UNROLL=6
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif

        p_diag%e_kinh(jc,jk,jb) =  &
          p_int%e_bln_c_s(jc,1,jb)*z_kin_hor_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
          p_int%e_bln_c_s(jc,2,jb)*z_kin_hor_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
          p_int%e_bln_c_s(jc,3,jb)*z_kin_hor_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))

        ENDDO
      ENDDO

      IF (istep == 1) THEN

        ! Interpolate contravariant correction to cell centers ...
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = nflatlev(p_patch%id), nlev
#else
!CDIR UNROLL=6
        DO jk = nflatlev(p_patch%id), nlev
          DO jc = i_startidx, i_endidx
#endif

            z_w_concorr_mc(jc,jk) =  &
              p_int%e_bln_c_s(jc,1,jb)*z_w_concorr_me(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
              p_int%e_bln_c_s(jc,2,jb)*z_w_concorr_me(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
              p_int%e_bln_c_s(jc,3,jb)*z_w_concorr_me(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))

          ENDDO
        ENDDO

        ! ... and to interface levels
        ! Remark: computation of w_concorr_c at nlevp1 is needed in solve_nh only
        ! because this serves solely for setting the lower boundary condition for w
        DO jk = nflatlev(p_patch%id)+1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_diag%w_concorr_c(jc,jk,jb) =                                &
              p_metrics%wgtfac_c(jc,jk,jb)*z_w_concorr_mc(jc,jk) +        &
             (1._wp - p_metrics%wgtfac_c(jc,jk,jb))*z_w_concorr_mc(jc,jk-1) 
          ENDDO
        ENDDO

      ENDIF
    ENDDO
!$OMP END DO

    rl_start = 7
    rl_end = min_rledge_int - 1

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (.NOT. lvn_only) THEN
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Compute v*grad w on edges (level nlevp1 is not needed because w(nlevp1) is diagnostic)
            ! Note: this implicitly includes a minus sign for the gradients, which is needed later on
            z_v_grad_w(je,jk,jb) = p_diag%vn_ie(je,jk,jb) * p_patch%edges%inv_dual_edge_length(je,jb)* &
             (p_prog%w(icidx(je,jb,1),jk,icblk(je,jb,1)) - p_prog%w(icidx(je,jb,2),jk,icblk(je,jb,2))) &
             + z_vt_ie(je,jk,jb) * p_patch%edges%inv_primal_edge_length(je,jb) *                       &
             p_patch%edges%system_orientation(je,jb) *                                                 &
             (z_w_v(ividx(je,jb,1),jk,ivblk(je,jb,1)) - z_w_v(ividx(je,jb,2),jk,ivblk(je,jb,2))) 

            ! Compute horizontal gradient of horizontal kinetic energy
            z_ddxn_ekin_e(je,jk,jb) = z_kin_hor_e(je,jk,jb) *                                     &
             (p_metrics%coeff_gradekin(je,1,jb) - p_metrics%coeff_gradekin(je,2,jb)) +            &
              p_metrics%coeff_gradekin(je,2,jb)*p_diag%e_kinh(icidx(je,jb,2),jk,icblk(je,jb,2)) - &
              p_metrics%coeff_gradekin(je,1,jb)*p_diag%e_kinh(icidx(je,jb,1),jk,icblk(je,jb,1)) 
          ENDDO
        ENDDO
      ELSE ! do not compute w tendency
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
!CDIR UNROLL=6
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Compute horizontal gradient of horizontal kinetic energy
            z_ddxn_ekin_e(je,jk,jb) = z_kin_hor_e(je,jk,jb) *                                     &
             (p_metrics%coeff_gradekin(je,1,jb) - p_metrics%coeff_gradekin(je,2,jb)) +            &
              p_metrics%coeff_gradekin(je,2,jb)*p_diag%e_kinh(icidx(je,jb,2),jk,icblk(je,jb,2)) - &
              p_metrics%coeff_gradekin(je,1,jb)*p_diag%e_kinh(icidx(je,jb,1),jk,icblk(je,jb,1)) 
          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO

    rl_start = 4
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    rl_start_2 = grf_bdywidth_c+1
    rl_end_2   = min_rlcell_int

    i_startblk_2 = p_patch%cells%start_blk(rl_start_2,1)
    i_endblk_2   = p_patch%cells%end_blk(rl_end_2,i_nchdom)

!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, i_startidx_2, i_endidx_2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      z_w_con_c(:,1:nlev,jb) = p_prog%w(:,1:nlev,jb)
      z_w_con_c(:,nlevp1,jb) = 0._wp

!CDIR UNROLL=5
      ! Contravariant vertical velocity on w points and interpolation to full levels
      DO jk = nlev, nflatlev(p_patch%id)+1, -1
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_w_con_c(jc,jk,jb) = z_w_con_c(jc,jk,jb) - p_diag%w_concorr_c(jc,jk,jb)
          z_w_con_c_full(jc,jk,jb) = 0.5_wp*(z_w_con_c(jc,jk,jb)+z_w_con_c(jc,jk+1,jb))
        ENDDO
      ENDDO

      DO jk = 1, nflatlev(p_patch%id)
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          z_w_con_c_full(jc,jk,jb) = 0.5_wp*(z_w_con_c(jc,jk,jb)+z_w_con_c(jc,jk+1,jb))
        ENDDO
      ENDDO

      ! The remaining computations are not needed in vn_only mode and only on prognostic grid points
      IF (lvn_only) CYCLE
      IF (jb < i_startblk_2 .OR. jb > i_endblk_2) CYCLE

      CALL get_indices_c(p_patch, jb, i_startblk_2, i_endblk_2, &
                         i_startidx_2, i_endidx_2, rl_start_2, rl_end_2)

      ! Interpolate horizontal advection of w from edges to cells
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx_2, i_endidx_2
!DIR$ IVDEP
        DO jk = 2, nlev
#else
!CDIR UNROLL=6
      DO jk = 1, nlev ! starting at level 2 would be sufficient, but this improves usage of unrolling
        DO jc = i_startidx_2, i_endidx_2
#endif
          p_diag%ddt_w_adv(jc,jk,jb,ntnd) =                                         &
            p_int%e_bln_c_s(jc,1,jb)*z_v_grad_w(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
            p_int%e_bln_c_s(jc,2,jb)*z_v_grad_w(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
            p_int%e_bln_c_s(jc,3,jb)*z_v_grad_w(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))
        ENDDO
      ENDDO

      ! Sum up remaining terms of vertical wind advection
      DO jk = 2, nlev
!DIR$ IVDEP
        DO jc = i_startidx_2, i_endidx_2
          p_diag%ddt_w_adv(jc,jk,jb,ntnd) = p_diag%ddt_w_adv(jc,jk,jb,ntnd) - z_w_con_c(jc,jk,jb) * &
            (p_prog%w(jc,jk-1,jb)*p_metrics%coeff1_dwdz(jc,jk,jb) -                                 &
             p_prog%w(jc,jk+1,jb)*p_metrics%coeff2_dwdz(jc,jk,jb) +                                 &
             p_prog%w(jc,jk,jb)*(p_metrics%coeff2_dwdz(jc,jk,jb) - p_metrics%coeff1_dwdz(jc,jk,jb)) )
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO

    rl_start = grf_bdywidth_e+1
    rl_end = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Sum up terms of horizontal wind advection: grad(Ekin_h) + vt*(f+relvort_e) + wcon_e*dv/dz
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
!DIR$ IVDEP
        DO jk = 1, nlev
#else
!CDIR UNROLL=2
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          p_diag%ddt_vn_adv(je,jk,jb,ntnd) = - ( z_ddxn_ekin_e(je,jk,jb) +           &
            p_diag%vt(je,jk,jb) * ( p_patch%edges%f_e(je,jb) + 0.5_wp*               &
           (p_diag%omega_z(ividx(je,jb,1),jk,ivblk(je,jb,1))   +                      &
            p_diag%omega_z(ividx(je,jb,2),jk,ivblk(je,jb,2))) ) +                     &
           (p_int%c_lin_e(je,1,jb)*z_w_con_c_full(icidx(je,jb,1),jk,icblk(je,jb,1)) + &
            p_int%c_lin_e(je,2,jb)*z_w_con_c_full(icidx(je,jb,2),jk,icblk(je,jb,2)))* &
           (p_diag%vn_ie(je,jk,jb) - p_diag%vn_ie(je,jk+1,jb))/   &
            p_metrics%ddqz_z_full_e(je,jk,jb) ) 
        ENDDO
      ENDDO

      ! Add centrifugal force for idealized gravity wave test
      IF (ltestcase.AND.nh_test_name=='dcmip_gw_32') THEN
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            p_diag%ddt_vn_adv(je,jk,jb,ntnd) = p_diag%ddt_vn_adv(je,jk,jb,ntnd) &
              &                              + fcfugal(je,jb)
          ENDDO
        ENDDO
      END IF

    ENDDO
!$OMP END DO 
     

    ! Apply some additional diffusion at grid points close to the stability limit for vertical advection.
    ! This happens extremely rarely in practice, the cases encountered so far were breaking gravity waves
    ! over the Andes at resolutions around 10 km.
    ! The diffusion term is added to ddt_vn_adv and ddt_w_adv, respectively, so that it acts also in the
    ! subsequent predictor step if advection is called in the corrector step only
    !
    IF (lextra_diffu .AND. .NOT. lvn_only) THEN
      rl_start = grf_bdywidth_c
      rl_end = min_rlcell_int - 1

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, ic, ie, je, jbe) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        ic = 0

        ! Collect grid points exceeding the specified stability threshold (75% of the theoretical limit) for w_con
        DO jk = MAX(3,nrdmax(jg)-2), nlev-4
          DO jc = i_startidx, i_endidx
            IF (ABS(z_w_con_c(jc,jk,jb)) > cfl_w_limit*p_metrics%ddqz_z_half(jc,jk,jb)) THEN
              ic = ic+1
              iclist(ic,jb) = jc
              iklist(ic,jb) = jk
            ENDIF
          ENDDO
        ENDDO

        icount(jb) = ic

      ENDDO
!$OMP END DO 

      ! If grid points are found, apply second-order diffusion on w and vn
      ! Note that this loop is neither vectorizable nor OpenMP-parallelizable (race conditions!)
!$OMP MASTER
      DO jb = i_startblk, i_endblk

        IF (icount(jb) > 0) THEN
          DO ic = 1, icount(jb)
            jc = iclist(ic,jb)
            jk = iklist(ic,jb)
            difcoef = scalfac_exdiff * MIN(1._wp - cfl_w_limit*dtime,                            &
              ABS(z_w_con_c(jc,jk,jb))*dtime/p_metrics%ddqz_z_half(jc,jk,jb) - cfl_w_limit*dtime )

            ! nabla2 diffusion on w
            p_diag%ddt_w_adv(jc,jk,jb,ntnd) = p_diag%ddt_w_adv(jc,jk,jb,ntnd)        + &
              difcoef * p_patch%cells%area(jc,jb) * (                                  &
              p_prog%w(jc,jk,jb)                          *p_int%geofac_n2s(jc,1,jb) + &
              p_prog%w(incidx(jc,jb,1),jk,incblk(jc,jb,1))*p_int%geofac_n2s(jc,2,jb) + &
              p_prog%w(incidx(jc,jb,2),jk,incblk(jc,jb,2))*p_int%geofac_n2s(jc,3,jb) + &
              p_prog%w(incidx(jc,jb,3),jk,incblk(jc,jb,3))*p_int%geofac_n2s(jc,4,jb)   )

            ! nabla2 diffusion on vn
            ! Note that edge points may be counted up to 4 times. This is intended as this turned out to
            ! be the simplest way to obtain correct parallelization. The diffusion coefficient is multiplied
            ! by 0.25 for this reason.
            DO ie = 1, 3
              je  = ieidx(jc,jb,ie)
              jbe = ieblk(jc,jb,ie)

              jk = iklist(ic,jb)
              p_diag%ddt_vn_adv(je,jk,jbe,ntnd) = p_diag%ddt_vn_adv(je,jk,jbe,ntnd)   +                               &
              0.25_wp * difcoef * p_patch%edges%area_edge(je,jbe) * (                                                 &
              p_int%geofac_grdiv(je,1,jbe)*p_prog%vn(je,jk,jbe)                          +                            &
              p_int%geofac_grdiv(je,2,jbe)*p_prog%vn(iqidx(je,jbe,1),jk,iqblk(je,jbe,1)) +                            &
              p_int%geofac_grdiv(je,3,jbe)*p_prog%vn(iqidx(je,jbe,2),jk,iqblk(je,jbe,2)) +                            &
              p_int%geofac_grdiv(je,4,jbe)*p_prog%vn(iqidx(je,jbe,3),jk,iqblk(je,jbe,3)) +                            &
              p_int%geofac_grdiv(je,5,jbe)*p_prog%vn(iqidx(je,jbe,4),jk,iqblk(je,jbe,4)) +                            &
              p_patch%edges%system_orientation(je,jbe)*p_patch%edges%inv_primal_edge_length(je,jbe) * (               &
              p_diag%omega_z(ividx(je,jbe,2),jk,ivblk(je,jbe,2))-p_diag%omega_z(ividx(je,jbe,1),jk,ivblk(je,jbe,1)) ) ) 

              jk = iklist(ic,jb) - 1
              p_diag%ddt_vn_adv(je,jk,jbe,ntnd) = p_diag%ddt_vn_adv(je,jk,jbe,ntnd)   +                               &
              0.25_wp * difcoef * p_patch%edges%area_edge(je,jbe) * (                                                 &
              p_int%geofac_grdiv(je,1,jbe)*p_prog%vn(je,jk,jbe)                          +                            &
              p_int%geofac_grdiv(je,2,jbe)*p_prog%vn(iqidx(je,jbe,1),jk,iqblk(je,jbe,1)) +                            &
              p_int%geofac_grdiv(je,3,jbe)*p_prog%vn(iqidx(je,jbe,2),jk,iqblk(je,jbe,2)) +                            &
              p_int%geofac_grdiv(je,4,jbe)*p_prog%vn(iqidx(je,jbe,3),jk,iqblk(je,jbe,3)) +                            &
              p_int%geofac_grdiv(je,5,jbe)*p_prog%vn(iqidx(je,jbe,4),jk,iqblk(je,jbe,4)) +                            &
              p_patch%edges%system_orientation(je,jbe)*p_patch%edges%inv_primal_edge_length(je,jbe) * (               &
              p_diag%omega_z(ividx(je,jbe,2),jk,ivblk(je,jbe,2))-p_diag%omega_z(ividx(je,jbe,1),jk,ivblk(je,jbe,1)) ) ) 
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!$OMP END MASTER
    ENDIF

!$OMP END PARALLEL

  END SUBROUTINE velocity_tendencies

END MODULE mo_velocity_advection
