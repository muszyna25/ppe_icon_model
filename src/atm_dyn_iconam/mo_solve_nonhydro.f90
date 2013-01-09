!>
!! mo_solve_nonhydro
!!
!! This module contains the nonhydrostatic dynamical core for the triangular version
!! Its routines were previously contained in mo_divergent_modes and mo_vector_operations
!! but have been extracted for better memory efficiency
!!
!! @author Guenther Zaengl, DWD
!!
!! @par Revision History
!! Initial release by Guenther Zaengl (2010-10-13) based on earlier work
!! by Almut Gassmann, MPI-M
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

MODULE mo_solve_nonhydro

  USE mo_kind,                 ONLY: wp
  USE mo_nonhydrostatic_config,ONLY: itime_scheme,iadv_rhotheta, igradp_method, l_open_ubc, &
                                     kstart_moist, lhdiff_rcf, divdamp_fac, divdamp_order,  &
                                     rayleigh_type, iadv_rcf
  USE mo_dynamics_config,   ONLY: idiv_method
  USE mo_parallel_config,   ONLY: nproma, p_test_run, itype_comm, use_dycore_barrier, &
    & use_icon_comm
  USE mo_run_config,        ONLY: ltimer, timers_level, lvert_nest
  USE mo_model_domain,      ONLY: t_patch
  USE mo_grid_config,       ONLY: l_limited_area
  USE mo_gridref_config,    ONLY: grf_intmethod_e
  USE mo_interpol_config,   ONLY: nudge_max_coeff
  USE mo_intp_data_strc,    ONLY: t_int_state
  USE mo_intp,              ONLY: cells2edges_scalar
  USE mo_intp_rbf,          ONLY: rbf_vec_interpol_edge
  USE mo_nonhydro_types,    ONLY: t_nh_state, t_nh_metrics, t_nh_diag, t_nh_prog, &
                                  t_buffer_memory
  USE mo_physical_constants,ONLY: cpd, rd, cvd, cvd_o_rd, grav, rd_o_cpd, p0ref
  USE mo_math_gradients,    ONLY: grad_green_gauss_cell
  USE mo_math_constants,    ONLY: pi, dbl_eps
  USE mo_math_divrot,       ONLY: div, rot_vertex, div_avg
  USE mo_vertical_grid,     ONLY: nrdmax, nflat_gradp
  USE mo_nh_init_utils,     ONLY: nflatlev
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,    ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int, &
    &                             min_rlcell, RAYLEIGH_CLASSIC, RAYLEIGH_KLEMP
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_advection_hflux,   ONLY: upwind_hflux_miura3
  USE mo_advection_traj,    ONLY: btraj
  USE mo_sync,              ONLY: SYNC_E, SYNC_C, sync_patch_array, sync_patch_array_mult, &
                                  sync_patch_array_gm
  USE mo_mpi,               ONLY: my_process_is_mpi_all_seq, work_mpi_barrier
  USE mo_timer,             ONLY: timer_solve_nh, timer_barrier, timer_start, timer_stop, &
                                  timer_solve_nh_p1, timer_solve_nh_p2, timer_solve_nh_exch
  USE mo_icon_comm_lib,     ONLY: icon_comm_sync

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  REAL(wp), PARAMETER :: rd_o_cvd = 1._wp / cvd_o_rd
  REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd
  REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
  REAL(wp), PARAMETER :: grav_o_cpd = grav / cpd

  PUBLIC :: solve_nh

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
  SUBROUTINE velocity_tendencies (p_prog, p_patch, p_int, p_metrics, p_diag,&
                                  ntnd, istep, lvn_only)

    ! Passed variables
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_int_state), TARGET, INTENT(IN):: p_int
    TYPE(t_nh_prog), INTENT(INOUT)       :: p_prog
    TYPE(t_nh_metrics), INTENT(IN)       :: p_metrics
    TYPE(t_nh_diag), INTENT(INOUT)       :: p_diag

    INTEGER, INTENT(IN)  :: ntnd  ! time level of ddt_adv fields used to store tendencies
    INTEGER, INTENT(IN)  :: istep ! 1: predictor step, 2: corrector step
    LOGICAL, INTENT(IN)  :: lvn_only ! true: compute only vn tendency

    ! Local variables
    INTEGER :: jb, jk, jc, je
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    REAL(wp):: z_w_concorr_me(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp):: z_w_concorr_mc(nproma,p_patch%nlev)
    REAL(wp):: z_w_con_c(nproma,p_patch%nlevp1,p_patch%nblks_c)
    REAL(wp):: z_w_con_c_full(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp):: z_kin_hor_e(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp):: z_ddxn_ekin_e(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp):: z_vnw(nproma,p_patch%nlevp1,p_patch%nblks_e)
    REAL(wp):: z_hadv_w(nproma,p_patch%nlevp1,p_patch%nblks_c)

    INTEGER,  DIMENSION(:,:,:), POINTER :: icidx, icblk, ieidx, ieblk, &
                                           ividx, ivblk, incidx, incblk
    INTEGER  :: nlev, nlevp1          !< number of full and half levels
    ! Local control variable for vertical nesting
    LOGICAL :: l_vert_nested

    !--------------------------------------------------------------------------

    IF ((lvert_nest) .AND. (p_patch%nshift > 0)) THEN  
      l_vert_nested = .TRUE.
    ELSE
      l_vert_nested = .FALSE.
    ENDIF

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

    i_nchdom   = MAX(1,p_patch%n_childdom)

    ! Tangential wind component using RBF reconstruction
    ! Note: vt is also diagnosed in divergent_modes. Thus, computation is needed in predictor step only
    IF (istep == 1) THEN
      CALL rbf_vec_interpol_edge(p_prog%vn, p_patch, p_int, p_diag%vt)
    ENDIF

    ! Compute vertical vorticity component at vertices
    CALL rot_vertex (p_prog%vn, p_patch, p_int, p_diag%omega_z, opt_rlend=min_rlvert_int-1)


!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)

    rl_start = 3
    rl_end = min_rledge_int - 2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Interpolate vn to interface levels and compute horizontal part of kinetic energy on edges
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx
          p_diag%vn_ie(je,jk,jb) =                                    &
            p_metrics%wgtfac_e(je,jk,jb)*p_prog%vn(je,jk,jb) +        &
           (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_prog%vn(je,jk-1,jb)
          z_kin_hor_e(je,jk,jb) = 0.5_wp*(p_prog%vn(je,jk,jb)*p_prog%vn(je,jk,jb) + &
            p_diag%vt(je,jk,jb)*p_diag%vt(je,jk,jb) )
        ENDDO
      ENDDO

      IF (istep == 1) THEN
        ! Compute contravariant correction for vertical velocity at interface levels
        ! (will be interpolated to cell centers below)
        DO jk = nflatlev(p_patch%id), nlev
          DO je = i_startidx, i_endidx
            z_w_concorr_me(je,jk,jb) =                              &
              p_prog%vn(je,jk,jb)*p_metrics%ddxn_z_full(je,jk,jb) + &
              p_diag%vt(je,jk,jb)*p_metrics%ddxt_z_full(je,jk,jb)
          ENDDO
        ENDDO
      ENDIF

      IF (.NOT. l_vert_nested) THEN
        ! Top and bottom levels
        DO je = i_startidx, i_endidx
          p_diag%vn_ie(je,1,jb) =                                &
            p_metrics%wgtfacq1_e(je,1,jb)*p_prog%vn(je,1,jb) +   &
            p_metrics%wgtfacq1_e(je,2,jb)*p_prog%vn(je,2,jb) + &
            p_metrics%wgtfacq1_e(je,3,jb)*p_prog%vn(je,3,jb)
          z_kin_hor_e(je,1,jb) = 0.5_wp*(p_prog%vn(je,1,jb)*p_prog%vn(je,1,jb) + &
            p_diag%vt(je,1,jb)*p_diag%vt(je,1,jb) )
          p_diag%vn_ie(je,nlevp1,jb) =                           &
            p_metrics%wgtfacq_e(je,1,jb)*p_prog%vn(je,nlev,jb) +   &
            p_metrics%wgtfacq_e(je,2,jb)*p_prog%vn(je,nlev-1,jb) + &
            p_metrics%wgtfacq_e(je,3,jb)*p_prog%vn(je,nlev-2,jb)
        ENDDO
      ELSE
        ! vn_ie(jk=1) is extrapolated using parent domain information in this case
        DO je = i_startidx, i_endidx
          p_diag%vn_ie(je,1,jb) = p_diag%vn_ie(je,2,jb) + p_diag%dvn_ie_ubc(je,jb)
          z_kin_hor_e(je,1,jb) = 0.5_wp*(p_prog%vn(je,1,jb)*p_prog%vn(je,1,jb) + &
            p_diag%vt(je,1,jb)*p_diag%vt(je,1,jb) )
          p_diag%vn_ie(je,nlevp1,jb) =                           &
            p_metrics%wgtfacq_e(je,1,jb)*p_prog%vn(je,nlev,jb) +   &
            p_metrics%wgtfacq_e(je,2,jb)*p_prog%vn(je,nlev-1,jb) + &
            p_metrics%wgtfacq_e(je,3,jb)*p_prog%vn(je,nlev-2,jb)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO

    rl_start = 2
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
          DO jc = i_startidx, i_endidx
            p_diag%w_concorr_c(jc,jk,jb) =                                &
              p_metrics%wgtfac_c(jc,jk,jb)*z_w_concorr_mc(jc,jk) +        &
             (1._wp - p_metrics%wgtfac_c(jc,jk,jb))*z_w_concorr_mc(jc,jk-1) 
          ENDDO
        ENDDO

      ENDIF
    ENDDO
!$OMP END DO

    rl_start = 4
    rl_end = min_rledge_int - 2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (.NOT. lvn_only) THEN
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Multiply vn_ie with w interpolated to edges for divergence computation
            z_vnw(je,jk,jb) = p_diag%vn_ie(je,jk,jb)*                            &
             ( p_int%c_lin_e(je,1,jb) * p_prog%w(icidx(je,jb,1),jk,icblk(je,jb,1)) &
             + p_int%c_lin_e(je,2,jb) * p_prog%w(icidx(je,jb,2),jk,icblk(je,jb,2)) )

            ! Compute horizontal gradient of horizontal kinetic energy
            z_ddxn_ekin_e(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb) *  &
             (p_diag%e_kinh(icidx(je,jb,2),jk,icblk(je,jb,2)) -                    &
              p_diag%e_kinh(icidx(je,jb,1),jk,icblk(je,jb,1)) )
          ENDDO
        ENDDO
      ELSE ! do not compute w tendency
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=6
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Compute horizontal gradient of horizontal kinetic energy
            z_ddxn_ekin_e(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb) *  &
             (p_diag%e_kinh(icidx(je,jb,2),jk,icblk(je,jb,2)) -                    &
              p_diag%e_kinh(icidx(je,jb,1),jk,icblk(je,jb,1)) )
          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO

    rl_start = 3
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (.NOT. lvn_only) THEN
        ! Compute horizontal advection of w: -(div(vn*w)-w*div(vn))
        ! (combined into one step for efficiency improvement)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
#endif
            z_hadv_w(jc,jk,jb) = p_prog%w(jc,jk,jb)* ( &
              p_diag%vn_ie(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
              p_diag%vn_ie(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
              p_diag%vn_ie(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%geofac_div(jc,3,jb))- &
             (z_vnw(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))       *p_int%geofac_div(jc,1,jb) + &
              z_vnw(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))       *p_int%geofac_div(jc,2,jb) + &
              z_vnw(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))       *p_int%geofac_div(jc,3,jb))

            z_w_con_c(jc,jk,jb) = p_prog%w(jc,jk,jb)
          ENDDO
        ENDDO

      ELSE ! do not compute w tendency
        z_w_con_c(:,1:nlev,jb) = p_prog%w(:,1:nlev,jb)
      ENDIF

      z_w_con_c(:,nlevp1,jb) = 0._wp

!CDIR UNROLL=5
      ! Contravariant vertical velocity on w points and interpolation to full levels
      DO jk = nlev, nflatlev(p_patch%id)+1, -1
        DO jc = i_startidx, i_endidx
          z_w_con_c(jc,jk,jb) = z_w_con_c(jc,jk,jb) - p_diag%w_concorr_c(jc,jk,jb)
          z_w_con_c_full(jc,jk,jb) = 0.5_wp*(z_w_con_c(jc,jk,jb)+z_w_con_c(jc,jk+1,jb))
        ENDDO
      ENDDO

      DO jk = 1, nflatlev(p_patch%id)
        DO jc = i_startidx, i_endidx
          z_w_con_c_full(jc,jk,jb) = 0.5_wp*(z_w_con_c(jc,jk,jb)+z_w_con_c(jc,jk+1,jb))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO

    IF (.NOT. lvn_only) THEN
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Apply cell averaging to the components of horizontal w advection
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2, nlev
#else
!CDIR UNROLL=4
        DO jk = 1, nlev ! starting at level 2 would be sufficient, but this improves usage of unrolling
          DO jc = i_startidx, i_endidx
#endif
            p_diag%ddt_w_adv(jc,jk,jb,ntnd) =                                         &
                z_hadv_w(jc,jk,jb)                          *p_int%c_bln_avg(jc,1,jb) &
              + z_hadv_w(incidx(jc,jb,1),jk,incblk(jc,jb,1))*p_int%c_bln_avg(jc,2,jb) &
              + z_hadv_w(incidx(jc,jb,2),jk,incblk(jc,jb,2))*p_int%c_bln_avg(jc,3,jb) &
              + z_hadv_w(incidx(jc,jb,3),jk,incblk(jc,jb,3))*p_int%c_bln_avg(jc,4,jb)
          ENDDO
        ENDDO

        ! Sum up remaining terms of vertical wind advection
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx
            p_diag%ddt_w_adv(jc,jk,jb,ntnd) = p_diag%ddt_w_adv(jc,jk,jb,ntnd)   &
              - z_w_con_c(jc,jk,jb)*(p_prog%w(jc,jk-1,jb)-p_prog%w(jc,jk+1,jb)) &
              * p_metrics%inv_ddqz_z_half2(jc,jk,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF

    rl_start = grf_bdywidth_e+1
    rl_end = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Sum up terms of horizontal wind advection:
      ! grad(Ekin_h) + vt*(f+relvort_e) + wcon_e*dv/dz
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
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
    ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  END SUBROUTINE velocity_tendencies


  !>
  !! solve_nh
  !!
  !! Version of divergent_modes used for the triangular grid
  !!
  !! @par Revision History
  !! Based on the initial release of divergent_modes by Almut Gassmann (2009-05-12)
  !! Modified by Guenther Zaengl starting on 2010-02-03
  SUBROUTINE solve_nh (p_nh, p_patch, p_int, bufr, mflx_avg, nnow, nnew, l_init, l_recompute, &
                       lsave_mflx, idyn_timestep, jstep, l_bdy_nudge, dtime)

    TYPE(t_nh_state),  TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state), TARGET, INTENT(IN)    :: p_int
    TYPE(t_patch),     TARGET, INTENT(IN)    :: p_patch
    TYPE(t_buffer_memory),     INTENT(INOUT) :: bufr

    REAL(wp), INTENT(IN)                     :: mflx_avg(:,:,:)

    ! Initialization switch that has to be .TRUE. at the initial time step only (not for restart)
    LOGICAL,                   INTENT(INOUT) :: l_init
    ! Switch to recompute velocity tendencies after a physics call irrespective of the time scheme option
    LOGICAL,                   INTENT(IN)    :: l_recompute
    ! Switch if mass flux needs to be saved for nest boundary interpolation tendency computation
    LOGICAL,                   INTENT(IN)    :: lsave_mflx
    ! Counter of dynamics time step within a large time step (ranges from 1 to iadv_rcf)
    INTEGER,                   INTENT(IN)    :: idyn_timestep
    ! Time step count since last boundary interpolation (ranges from 0 to 2*iadv_rcf-1)
    INTEGER,                   INTENT(IN)    :: jstep
    ! Switch to determine if boundary nudging is executed
    LOGICAL,                   INTENT(IN)    :: l_bdy_nudge
    ! Time levels
    INTEGER,                   INTENT(IN)    :: nnow, nnew
    ! Time step
    REAL(wp),                  INTENT(IN)    :: dtime

    ! Local variables
    INTEGER  :: jb, jk, jc, je
    INTEGER  :: nlev, nlevp1              !< number of full levels
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom, ishift
    INTEGER  :: rl_start, rl_end, istep, ntl1, ntl2, nvar, nshift
    INTEGER  :: ic, ie, ilc0, ibc0, ikp1, ikp2

    REAL(wp) :: z_theta_v_fl_e  (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_theta_v_e     (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_rho_e         (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_gradh_exner   (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_graddiv_vn    (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_graddiv2_vn   (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_w_concorr_me  (nproma,p_patch%nlev,p_patch%nblks_e), &
                z_distv_bary    (nproma,p_patch%nlev  ,p_patch%nblks_e,2)

    INTEGER ::  z_cell_indices  (nproma,p_patch%nlev  ,p_patch%nblks_e,2)

    REAL(wp), TARGET :: z_vn_avg (nproma,p_patch%nlev,p_patch%nblks_e)

    REAL(wp) :: z_mass_fl_div   (nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_theta_v_fl_div(nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_th_ddz_exner_c(nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_dexner_dz_c (2,nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_exner_ex_pr   (nproma,p_patch%nlevp1,p_patch%nblks_c), & ! nlevp1 is intended here
                z_exner_pr      (nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_grad_rth      (nproma,4,p_patch%nlev,p_patch%nblks_c), &
                z_rth_pr        (nproma,2,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) :: z_w_expl        (nproma,p_patch%nlevp1),          &
                z_contr_w_fl_l  (nproma,p_patch%nlevp1),          &
                z_alpha         (nproma,p_patch%nlevp1),          &
                z_beta          (nproma,p_patch%nlev  ),          &
                z_gamma         (nproma,p_patch%nlev  ),          &
                z_a             (nproma,p_patch%nlevp1),          &
                z_b             (nproma,p_patch%nlevp1),          &
                z_c             (nproma,p_patch%nlev  ),          &
                z_g             (nproma,p_patch%nlev  ),          &
                z_q             (nproma,p_patch%nlev  ),          &
                z_rho_expl      (nproma,p_patch%nlev  ),          &
                z_exner_expl    (nproma,p_patch%nlev  ),          &
                z_exner_ic      (nproma,p_patch%nlevp1),          &
                z_theta_v_pr_mc (nproma,p_patch%nlev  ),          &
                z_theta_v_pr_ic (nproma,p_patch%nlevp1),          &
                z_w_concorr_mc  (nproma,p_patch%nlev  ),          &
                z_thermal_exp   (nproma,p_patch%nblks_c),         &
                z_mflx_top      (nproma,p_patch%nblks_c),         &
                z_hydro_corr    (nproma,p_patch%nblks_e)


    REAL(wp):: z_theta1, z_theta2, z_raylfac, wgt_nnow, wgt_nnew, scal_divdamp, z_w_lim, &
               dt_shift, bdy_divdamp
    INTEGER :: nproma_gradp, nblks_gradp, npromz_gradp, nlen_gradp, jk_start
    LOGICAL :: lcompute, lcleanup, lvn_only

    ! Local variables to control vertical nesting
    LOGICAL :: l_vert_nested, l_child_vertnest

    ! Pointers to cell indices
    INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk
    ! Pointers to edge indices
    INTEGER,  DIMENSION(:,:,:),   POINTER :: ieidx, ieblk
    ! Pointers to vertical neighbor indices for pressure gradient computation
    INTEGER,  DIMENSION(:,:,:,:),   POINTER :: ikidx
    ! Pointers to quad edge indices
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iqidx, iqblk
    ! Pointer to velocity field used for mass flux computation
    REAL(wp), DIMENSION(:,:,:),   POINTER :: ptr_vn
    ! Pointers needed for igradp_method = 3
    INTEGER,  DIMENSION(:),   POINTER :: iplev, ipeidx, ipeblk
!     REAL(wp) :: sphere_radius_squared

    !-----------------------------------------------------------------------
!     sphere_radius_squared = grid_sphere_radius * grid_sphere_radius

    !-------------------------------------------------------------------
    IF (use_dycore_barrier) THEN
      CALL timer_start(timer_barrier)
      CALL work_mpi_barrier()
      CALL timer_stop(timer_barrier)
    ENDIF
    !-------------------------------------------------------------------
    
    IF (ltimer) CALL timer_start(timer_solve_nh)
    IF (timers_level > 5) CALL timer_start(timer_solve_nh_p1)

    IF (lvert_nest .AND. (p_patch%nshift_total > 0)) THEN  
      l_vert_nested = .TRUE.
    ELSE
      l_vert_nested = .FALSE.
    ENDIF
    IF (lvert_nest .AND. p_patch%n_childdom > 0 .AND.              &
      (p_patch%nshift_child > 0 .OR. p_patch%nshift_total > 0)) THEN
      l_child_vertnest = .TRUE.
      nshift = p_patch%nshift_child + 1
    ELSE
      l_child_vertnest = .FALSE.
      nshift = 0
    ENDIF

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Set pointers to neighbor cells
    icidx => p_patch%edges%cell_idx
    icblk => p_patch%edges%cell_blk

    ! Set pointers to neighbor edges
    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    ! Set pointer to vertical neighbor indices for pressure gradient
    ikidx => p_nh%metrics%vertidx_gradp

    ! Set pointers to quad edges
    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk

    ! scaling factor for divergence damping: divdamp_fac*delta_x**2
    IF (divdamp_order == 2) THEN
!        scal_divdamp = divdamp_fac*4._wp*pi*sphere_radius_squared &
!                         & / REAL(20*nroot**2*4**(p_patch%level),wp)
      ! get the mean cell area directly from the p_patch%geometry_info
      ! this should work for any (near-uniform) grid geometry
      scal_divdamp = divdamp_fac * p_patch%geometry_info%mean_cell_area 
    ELSE IF (divdamp_order == 4) THEN
!        scal_divdamp = -divdamp_fac*(4._wp*pi*sphere_radius_squared &
!                        & /REAL(20*nroot**2*4**(p_patch%level),wp))**2
      ! get the mean cell area directly from the p_patch%geometry_info
      ! this should work for any (near-uniform) grid geometry
       scal_divdamp = -divdamp_fac * (p_patch%geometry_info%mean_cell_area**2)
    ENDIF

    ! Time increment for backward-shifting of lateral boundary mass flux 
    dt_shift = dtime*(0.5_wp*REAL(iadv_rcf,wp)-0.25_wp)

    ! Coefficient for reduced fourth-order divergence damping along nest boundaries
    bdy_divdamp = 0.75_wp/(nudge_max_coeff + dbl_eps)*ABS(scal_divdamp)

    ! Set pointer to velocity field that is used for mass flux computation
    IF (idiv_method == 1) THEN
      ptr_vn => z_vn_avg
    ELSE
      ptr_vn => p_nh%prog(nnew)%vn
    ENDIF

    IF (p_test_run) THEN
      z_rho_e     = 0._wp
      z_theta_v_e = 0._wp
    ENDIF

    ! Set time levels of ddt_adv fields for call to velocity_tendencies
    IF (itime_scheme >= 4) THEN ! Velocity advection averaging nnow and nnew tendencies
      ntl1 = nnow
      ntl2 = nnew
    ELSE                        ! Velocity advection is taken at nnew only
      ntl1 = 1
      ntl2 = 1
    ENDIF

    ! Weighting coefficients for velocity advection if tendency averaging is used
    ! The off-centering specified here turned out to be beneficial to numerical
    ! stability in extreme situations
    wgt_nnow = 0.25_wp
    wgt_nnew = 1._wp - wgt_nnow

    i_nchdom   = MAX(1,p_patch%n_childdom)

    DO istep = 1, 2

      IF (istep == 1) THEN ! predictor step
        IF (itime_scheme >= 6 .OR. l_init .OR. l_recompute) THEN
          IF (itime_scheme < 6 .AND. .NOT. l_init) THEN
            lvn_only = .TRUE. ! Recompute only vn tendency
          ELSE
            lvn_only = .FALSE.
          ENDIF
          CALL velocity_tendencies(p_nh%prog(nnow),p_patch,p_int,p_nh%metrics,&
                                   p_nh%diag,ntl1,istep,lvn_only)
        ENDIF
        nvar = nnow
      ELSE                 ! corrector step
        lvn_only = .FALSE.
        CALL velocity_tendencies(p_nh%prog(nnew),p_patch,p_int,p_nh%metrics,&
                                 p_nh%diag,ntl2,istep,lvn_only)
        nvar = nnew
      ENDIF

      l_init = .FALSE. ! should be .TRUE. only at initial predictor step

    ! Compute rho and theta at edges
    IF (istep == 1) THEN
      IF (iadv_rhotheta == 2) THEN

        ! Operations from upwind_hflux_miura are inlined in order to process both
        ! fields in one step
        CALL btraj(p_patch, p_int, p_nh%prog(nnow)%vn, p_nh%diag%vt, &
                   0.5_wp*dtime, z_cell_indices, z_distv_bary,       &
                   opt_rlstart=7, opt_rlend=min_rledge_int-1 )

        CALL grad_green_gauss_cell(p_nh%prog(nnow)%rho, p_patch, p_int, z_grad_rth, &
                                   opt_rlend=min_rlcell_int-1, opt_dynmode=.TRUE.,  &
                                   opt_ccin2=p_nh%prog(nnow)%theta_v,               &
                                   opt_ref1=p_nh%metrics%rho_ref_mc,                &
                                   opt_ref2=p_nh%metrics%theta_ref_mc,              &
                                   opt_ccpr=z_rth_pr                                )

        rl_start = 7
        rl_end   = min_rledge_int-1

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

        i_startblk = p_patch%edges%start_blk(min_rledge_int-2,i_nchdom)
        i_endblk   = p_patch%edges%end_blk  (min_rledge_int-3,i_nchdom)

        ! Initialize halo edges with zero in order to avoid access of uninitialized array elements
!$OMP WORKSHARE
        z_rho_e    (:,:,i_startblk:i_endblk) = 0._wp
        z_theta_v_e(:,:,i_startblk:i_endblk) = 0._wp
!$OMP END WORKSHARE

        i_startblk = p_patch%edges%start_blk(rl_start,1)
        i_endblk   = p_patch%edges%end_blk  (rl_end,i_nchdom)

        ! initialize also nest boundary points with zero
        IF (p_patch%id > 1 .OR. l_limited_area) THEN
!$OMP WORKSHARE
          z_rho_e    (:,:,1:i_startblk) = 0._wp
          z_theta_v_e(:,:,1:i_startblk) = 0._wp
!$OMP END WORKSHARE
        ENDIF

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,ilc0,ibc0), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

!CDIR UNROLL=4
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx

              ilc0 = z_cell_indices(je,jk,jb,1)
              ibc0 = z_cell_indices(je,jk,jb,2)  

              ! Calculate "edge values" of rho and theta_v
              ! Note: z_rth_pr contains the perturbation values of rho and theta_v,
              ! and the corresponding gradients are stored in z_grad_rth.
              z_rho_e(je,jk,jb) = p_nh%metrics%rho_ref_me(je,jk,jb)     &
                +                            z_rth_pr(ilc0,1,jk,ibc0)   &
                + z_distv_bary(je,jk,jb,1) * z_grad_rth(ilc0,1,jk,ibc0) &
                + z_distv_bary(je,jk,jb,2) * z_grad_rth(ilc0,2,jk,ibc0)

              z_theta_v_e(je,jk,jb) = p_nh%metrics%theta_ref_me(je,jk,jb) &
                +                            z_rth_pr(ilc0,2,jk,ibc0)     &
                + z_distv_bary(je,jk,jb,1) * z_grad_rth(ilc0,3,jk,ibc0)   &
                + z_distv_bary(je,jk,jb,2) * z_grad_rth(ilc0,4,jk,ibc0)

            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels

        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ELSE IF (iadv_rhotheta == 3) THEN

        lcompute =.TRUE.
        lcleanup =.FALSE.
        ! First call: compute backward trajectory with wind at time level nnow
        CALL upwind_hflux_miura3(p_patch, p_nh%prog(nnow)%rho, p_nh%prog(nnow)%vn,    &
                                p_nh%prog(nnow)%vn, dtime, p_int, lcompute, lcleanup, &
                                0, z_rho_e, opt_rlstart=7, opt_lout_edge=.TRUE.,      &
                                opt_real_vt=p_nh%diag%vt )

        ! Second call: compute only reconstructed value for flux divergence
        lcompute =.FALSE.
        lcleanup =.TRUE.
        CALL upwind_hflux_miura3(p_patch, p_nh%prog(nnow)%theta_v, p_nh%prog(nnow)%vn, &
                                p_nh%prog(nnow)%vn, dtime, p_int, lcompute, lcleanup,  &
                                0, z_theta_v_e, opt_rlstart=7, opt_lout_edge=.TRUE.,   &
                                opt_real_vt=p_nh%diag%vt )

      ELSE

        ! density at edges
        CALL cells2edges_scalar(p_nh%prog(nnow)%rho,p_patch,p_int%c_lin_e,z_rho_e)
        ! virtual potential temperature at edges
        CALL cells2edges_scalar(p_nh%prog(nnow)%theta_v,p_patch,p_int%c_lin_e,z_theta_v_e)

      ENDIF

    ENDIF ! istep = 1

    ! Preparations for igradp_method = 3
    IF (istep == 1 .AND. (igradp_method == 3 .OR. igradp_method == 5)) THEN
      
      iplev  => p_nh%metrics%pg_vertidx
      ipeidx => p_nh%metrics%pg_edgeidx
      ipeblk => p_nh%metrics%pg_edgeblk

      nproma_gradp = MIN(nproma,256)
      nblks_gradp  = INT(p_nh%metrics%pg_listdim/nproma_gradp)
      npromz_gradp = MOD(p_nh%metrics%pg_listdim,nproma_gradp)
      IF (npromz_gradp > 0) THEN
        nblks_gradp = nblks_gradp + 1
      ELSE
        npromz_gradp = nproma_gradp
      ENDIF

    ENDIF

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    rl_start = 3
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,z_exner_ic,z_theta_v_pr_mc,z_theta_v_pr_ic) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (istep == 1) THEN ! to be executed in predictor step only

!CDIR UNROLL=2
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            ! extrapolated perturbation Exner pressure (used for horizontal gradients only)
            z_exner_ex_pr(jc,jk,jb) = - p_nh%metrics%exner_ref_mc(jc,jk,jb) +              &
              (1._wp + p_nh%metrics%exner_exfac(jc,jk,jb))*p_nh%prog(nnow)%exner(jc,jk,jb) &
                     - p_nh%metrics%exner_exfac(jc,jk,jb) *p_nh%diag%exner_old(jc,jk,jb)

            ! non-extrapolated perturbation Exner pressure
            z_exner_pr(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb) - &
              p_nh%metrics%exner_ref_mc(jc,jk,jb)

            ! Now save current time level in exner_old
            p_nh%diag%exner_old(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb)

          ENDDO
        ENDDO

        ! The purpose of the extra level of exner_pr is to simplify coding for
        ! igradp_method=4/5. It is multiplied with zero and thus actually not used
        z_exner_ex_pr(:,nlevp1,jb) = 0._wp

        IF (l_open_ubc .AND. .NOT. l_vert_nested) THEN
          ! Compute contribution of thermal expansion to vertical wind at model top
          ! Isothermal expansion is assumed
          z_thermal_exp(:,jb) = 0._wp
!CDIR UNROLL=4
          DO jk = 2, nlev
            DO jc = i_startidx, i_endidx
              z_thermal_exp(jc,jb) = z_thermal_exp(jc,jb) + cvd_o_rd                &
                * (p_nh%diag%ddt_exner(jc,jk,jb)+p_nh%diag%ddt_exner_phy(jc,jk,jb)) &
                /  p_nh%prog(nnow)%exner(jc,jk,jb)*p_nh%metrics%ddqz_z_full(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF

        IF (igradp_method <= 3) THEN
          ! Perturbation Exner pressure on bottom half level
          DO jc = i_startidx, i_endidx
            z_exner_ic(jc,nlevp1) =                                         &
              p_nh%metrics%wgtfacq_c(jc,1,jb)*z_exner_ex_pr(jc,nlev  ,jb) + &
              p_nh%metrics%wgtfacq_c(jc,2,jb)*z_exner_ex_pr(jc,nlev-1,jb) + &
              p_nh%metrics%wgtfacq_c(jc,3,jb)*z_exner_ex_pr(jc,nlev-2,jb)
          ENDDO

!CDIR UNROLL=3
          DO jk = nlev, MAX(2,nflatlev(p_patch%id)), -1
            DO jc = i_startidx, i_endidx
              ! Exner pressure on remaining half levels for metric correction term
              z_exner_ic(jc,jk) =                                              &
                p_nh%metrics%wgtfac_c(jc,jk,jb)*z_exner_ex_pr(jc,jk,jb) +      &
                (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_exner_ex_pr(jc,jk-1,jb)

              ! First vertical derivative of perturbation Exner pressure
              z_dexner_dz_c(1,jc,jk,jb) =                    &
               (z_exner_ic(jc,jk) - z_exner_ic(jc,jk+1)) *   &
                p_nh%metrics%inv_ddqz_z_full(jc,jk,jb)
            ENDDO
          ENDDO

          IF (nflatlev(p_patch%id) == 1) THEN
            ! Perturbation Exner pressure on top half level
            DO jc = i_startidx, i_endidx
            z_exner_ic(jc,1) =                                          &
              p_nh%metrics%wgtfacq1_c(jc,1,jb)*z_exner_ex_pr(jc,1,jb) + &
              p_nh%metrics%wgtfacq1_c(jc,2,jb)*z_exner_ex_pr(jc,2,jb) + &
              p_nh%metrics%wgtfacq1_c(jc,3,jb)*z_exner_ex_pr(jc,3,jb)

              ! First vertical derivative of perturbation Exner pressure
              z_dexner_dz_c(1,jc,1,jb) =                    &
               (z_exner_ic(jc,1) - z_exner_ic(jc,2)) *   &
                p_nh%metrics%inv_ddqz_z_full(jc,1,jb)
            ENDDO
          ENDIF

        ENDIF

      ENDIF ! istep = 1

      z_theta_v_pr_mc(i_startidx:i_endidx,1) =  0.5_wp *     &
        (p_nh%prog(nnow)%theta_v(i_startidx:i_endidx,1,jb) + &
        p_nh%prog(nvar)%theta_v(i_startidx:i_endidx,1,jb)) - &
        p_nh%metrics%theta_ref_mc(i_startidx:i_endidx,1,jb)

!CDIR UNROLL=8
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          ! density at interface levels for vertical flux divergence computation
          p_nh%diag%rho_ic(jc,jk,jb) = 0.5_wp*(                                     &
            p_nh%metrics%wgtfac_c(jc,jk,jb)*(p_nh%prog(nnow)%rho(jc,jk,jb) +        &
            p_nh%prog(nvar)%rho(jc,jk,jb))+(1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))* &
            (p_nh%prog(nnow)%rho(jc,jk-1,jb)+p_nh%prog(nvar)%rho(jc,jk-1,jb)) )

          ! perturbation virtual potential temperature at main levels
          z_theta_v_pr_mc(jc,jk) = 0.5_wp*(p_nh%prog(nnow)%theta_v(jc,jk,jb) +     &
            p_nh%prog(nvar)%theta_v(jc,jk,jb)) - p_nh%metrics%theta_ref_mc(jc,jk,jb)

          ! perturbation virtual potential temperature at interface levels
          z_theta_v_pr_ic(jc,jk) = &
            p_nh%metrics%wgtfac_c(jc,jk,jb)*z_theta_v_pr_mc(jc,jk) +       &
            (1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*z_theta_v_pr_mc(jc,jk-1)

          ! virtual potential temperature at interface levels
          p_nh%diag%theta_v_ic(jc,jk,jb) = p_nh%metrics%theta_ref_ic(jc,jk,jb) + &
            z_theta_v_pr_ic(jc,jk)

          ! vertical pressure gradient * theta_v
          z_th_ddz_exner_c(jc,jk,jb) = p_nh%metrics%vwind_expl_wgt(jc,jb)* &
            p_nh%diag%theta_v_ic(jc,jk,jb) * (z_exner_pr(jc,jk-1,jb)-      &
            z_exner_pr(jc,jk,jb)) / p_nh%metrics%ddqz_z_half(jc,jk,jb) +   &
            z_theta_v_pr_ic(jc,jk)*p_nh%metrics%d_exner_dz_ref_ic(jc,jk,jb)
        ENDDO
      ENDDO

      ! rho and theta at top level (in case of vertical nesting, upper boundary conditions 
      !                             are set in the vertical solver loop)
      IF (l_open_ubc .AND. .NOT. l_vert_nested) THEN
        DO jc = i_startidx, i_endidx
          p_nh%diag%theta_v_ic(jc,1,jb) = p_nh%metrics%theta_ref_ic(jc,1,jb) + &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*z_theta_v_pr_mc(jc,1) +           &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*z_theta_v_pr_mc(jc,2) +           &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*z_theta_v_pr_mc(jc,3)
          p_nh%diag%rho_ic(jc,1,jb) =  0.5_wp*(                              &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*p_nh%prog(nnow)%rho(jc,1,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*p_nh%prog(nnow)%rho(jc,2,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*p_nh%prog(nnow)%rho(jc,3,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*p_nh%prog(nvar)%rho(jc,1,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*p_nh%prog(nvar)%rho(jc,2,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*p_nh%prog(nvar)%rho(jc,3,jb) )
        ENDDO
      ENDIF

      IF (istep == 1) THEN

        ! Perturbation theta at top and surface levels
        DO jc = i_startidx, i_endidx
          z_theta_v_pr_ic(jc,1)      = 0._wp
          z_theta_v_pr_ic(jc,nlevp1) =                                   &
            p_nh%metrics%wgtfacq_c(jc,1,jb)*z_theta_v_pr_mc(jc,nlev) +   &
            p_nh%metrics%wgtfacq_c(jc,2,jb)*z_theta_v_pr_mc(jc,nlev-1) + &
            p_nh%metrics%wgtfacq_c(jc,3,jb)*z_theta_v_pr_mc(jc,nlev-2)
            p_nh%diag%theta_v_ic(jc,nlevp1,jb) =                                  &
              p_nh%metrics%theta_ref_ic(jc,nlevp1,jb) + z_theta_v_pr_ic(jc,nlevp1)
        ENDDO

        IF (igradp_method <= 3) THEN
!CDIR UNROLL=3
          DO jk = nflat_gradp(p_patch%id), nlev
            DO jc = i_startidx, i_endidx
              ! Second vertical derivative of perturbation Exner pressure (hydrostatic approximation)
              z_dexner_dz_c(2,jc,jk,jb) = -0.5_wp *                                &
               ((z_theta_v_pr_ic(jc,jk) - z_theta_v_pr_ic(jc,jk+1)) *              &
                p_nh%metrics%d2dexdz2_fac1_mc(jc,jk,jb) + z_theta_v_pr_mc(jc,jk)*  &
                p_nh%metrics%d2dexdz2_fac2_mc(jc,jk,jb))
            ENDDO
          ENDDO
        ENDIF

      ENDIF ! istep == 1

    ENDDO
!$OMP END DO

    ! Computations at edge points

    rl_start = grf_bdywidth_e + 1   ! boundary update follows below
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

    IF (istep == 1) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_theta1,z_theta2,ikp1,ikp2) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Store values at nest interface levels
        IF (idyn_timestep == 1 .AND. l_child_vertnest) THEN
          DO je = i_startidx, i_endidx
            p_nh%diag%dvn_ie_int(je,jb) = p_nh%diag%vn_ie(je,nshift,jb) - &
                                          p_nh%diag%vn_ie(je,nshift+1,jb)
          ENDDO
        ENDIF

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nflatlev(p_patch%id)-1
#else
!CDIR UNROLL=6
        DO jk = 1, nflatlev(p_patch%id)-1
          DO je = i_startidx, i_endidx
#endif
            ! horizontal gradient of Exner pressure where coordinate surfaces are flat
            z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)* &
             (z_exner_ex_pr(icidx(je,jb,2),jk,icblk(je,jb,2)) -                  &
              z_exner_ex_pr(icidx(je,jb,1),jk,icblk(je,jb,1)) )
          ENDDO
        ENDDO

        IF (igradp_method <= 3) THEN
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
            DO jk = nflatlev(p_patch%id), nflat_gradp(p_patch%id)
#else
!CDIR UNROLL=6
          DO jk = nflatlev(p_patch%id), nflat_gradp(p_patch%id)
            DO je = i_startidx, i_endidx
#endif
              ! horizontal gradient of Exner pressure, including metric correction
              z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)*         &
               (z_exner_ex_pr(icidx(je,jb,2),jk,icblk(je,jb,2)) -                          &
                z_exner_ex_pr(icidx(je,jb,1),jk,icblk(je,jb,1)) ) -                        &
                p_nh%metrics%ddxn_z_full(je,jk,jb) *                                       &
               (p_int%c_lin_e(je,1,jb)*z_dexner_dz_c(1,icidx(je,jb,1),jk,icblk(je,jb,1)) + &
                p_int%c_lin_e(je,2,jb)*z_dexner_dz_c(1,icidx(je,jb,2),jk,icblk(je,jb,2)))

            ENDDO
          ENDDO

        ! remark: loop exchange is not beneficial here because of 3D indirect addressing
!CDIR UNROLL=5
          DO jk = nflat_gradp(p_patch%id)+1, nlev
            DO je = i_startidx, i_endidx

              ! horizontal gradient of Exner pressure, Taylor-expansion-based reconstruction
              z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)*         &
               (z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2)) +           &
                p_nh%metrics%zdiff_gradp(2,je,jk,jb)*                                      &
               (z_dexner_dz_c(1,icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2)) +         &
                p_nh%metrics%zdiff_gradp(2,je,jk,jb)*                                      &
                z_dexner_dz_c(2,icidx(je,jb,2),ikidx(2,je,jk,jb),icblk(je,jb,2))) -        &
               (z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)) +           &
                p_nh%metrics%zdiff_gradp(1,je,jk,jb)*                                      &
               (z_dexner_dz_c(1,icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)) +         &
                p_nh%metrics%zdiff_gradp(1,je,jk,jb)*                                      &
                z_dexner_dz_c(2,icidx(je,jb,1),ikidx(1,je,jk,jb),icblk(je,jb,1)))))

            ENDDO
          ENDDO
        ELSE IF (igradp_method == 4 .OR. igradp_method == 5) THEN
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
            DO jk = nflatlev(p_patch%id), nlev
#else
          DO jk = nflatlev(p_patch%id), nlev
            DO je = i_startidx, i_endidx
#endif
              ! horizontal gradient of Exner pressure, cubic/quadratic interpolation
              z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)*   &
               (z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)-1,icblk(je,jb,2)) *   &
                p_nh%metrics%coeff_gradp(5,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)  ,icblk(je,jb,2)) *   &
                p_nh%metrics%coeff_gradp(6,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)+1,icblk(je,jb,2)) *   &
                p_nh%metrics%coeff_gradp(7,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,2),ikidx(2,je,jk,jb)+2,icblk(je,jb,2)) *   &
                p_nh%metrics%coeff_gradp(8,je,jk,jb) -                               &
               (z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)-1,icblk(je,jb,1)) *   &
                p_nh%metrics%coeff_gradp(1,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)  ,icblk(je,jb,1)) *   &
                p_nh%metrics%coeff_gradp(2,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)+1,icblk(je,jb,1)) *   &
                p_nh%metrics%coeff_gradp(3,je,jk,jb) +                               &
                z_exner_ex_pr(icidx(je,jb,1),ikidx(1,je,jk,jb)+2,icblk(je,jb,1)) *   &
                p_nh%metrics%coeff_gradp(4,je,jk,jb)) )

            ENDDO
          ENDDO
        ENDIF

        ! compute hydrostatically approximated correction term that replaces downward extrapolation
        IF (igradp_method == 3) THEN
          
          DO je = i_startidx, i_endidx

            z_theta1 = &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb),icblk(je,jb,1)) +  &
              p_nh%metrics%zdiff_gradp(1,je,nlev,jb)*                                       &
             (p_nh%diag%theta_v_ic(icidx(je,jb,1),ikidx(1,je,nlev,jb),  icblk(je,jb,1)) -   &
              p_nh%diag%theta_v_ic(icidx(je,jb,1),ikidx(1,je,nlev,jb)+1,icblk(je,jb,1))) *  &
              p_nh%metrics%inv_ddqz_z_full(icidx(je,jb,1),ikidx(1,je,nlev,jb),icblk(je,jb,1))

            z_theta2 = &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb),icblk(je,jb,2)) +  &
              p_nh%metrics%zdiff_gradp(2,je,nlev,jb)*                                       &
             (p_nh%diag%theta_v_ic(icidx(je,jb,2),ikidx(2,je,nlev,jb),  icblk(je,jb,2)) -   &
              p_nh%diag%theta_v_ic(icidx(je,jb,2),ikidx(2,je,nlev,jb)+1,icblk(je,jb,2))) *  &
              p_nh%metrics%inv_ddqz_z_full(icidx(je,jb,2),ikidx(2,je,nlev,jb),icblk(je,jb,2))

            z_hydro_corr(je,jb) = grav_o_cpd*p_patch%edges%inv_dual_edge_length(je,jb)*    &
              (z_theta2-z_theta1)*4._wp/(z_theta1+z_theta2)**2

          ENDDO
        ELSE IF (igradp_method == 5) THEN
          
          DO je = i_startidx, i_endidx

            ikp1 = MIN(nlev,ikidx(1,je,nlev,jb)+2)
            ikp2 = MIN(nlev,ikidx(2,je,nlev,jb)+2)

            z_theta1 =                                                                       &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb)-1,icblk(je,jb,1)) * &
              p_nh%metrics%coeff_gradp(1,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb)  ,icblk(je,jb,1)) * &
              p_nh%metrics%coeff_gradp(2,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb)+1,icblk(je,jb,1)) * &
              p_nh%metrics%coeff_gradp(3,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikp1                 ,icblk(je,jb,1)) * &
              p_nh%metrics%coeff_gradp(4,je,nlev,jb) 

            z_theta2 =                                                                       &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb)-1,icblk(je,jb,2)) * &
              p_nh%metrics%coeff_gradp(5,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb)  ,icblk(je,jb,2)) * &
              p_nh%metrics%coeff_gradp(6,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb)+1,icblk(je,jb,2)) * &
              p_nh%metrics%coeff_gradp(7,je,nlev,jb) +                                         &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikp2                 ,icblk(je,jb,2)) * &
              p_nh%metrics%coeff_gradp(8,je,nlev,jb) 

            z_hydro_corr(je,jb) = grav_o_cpd*p_patch%edges%inv_dual_edge_length(je,jb)*    &
              (z_theta2-z_theta1)*4._wp/(z_theta1+z_theta2)**2

          ENDDO
        ENDIF

      ENDDO
!$OMP END DO
    ENDIF ! istep = 1

    IF (istep == 1 .AND. (igradp_method == 3 .OR. igradp_method == 5)) THEN

!$OMP DO PRIVATE(jb,je,ie,nlen_gradp,ishift) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_gradp
        IF (jb == nblks_gradp) THEN
          nlen_gradp = npromz_gradp
        ELSE
          nlen_gradp = nproma_gradp
        ENDIF
        ishift = (jb-1)*nproma_gradp
!CDIR NODEP,VOVERTAKE,VOB
        DO je = 1, nlen_gradp
          ie = ishift+je

          z_gradh_exner(ipeidx(ie),iplev(ie),ipeblk(ie))  =              &
            z_gradh_exner(ipeidx(ie),iplev(ie),ipeblk(ie)) +             &
            p_nh%metrics%pg_exdist(ie)*z_hydro_corr(ipeidx(ie),ipeblk(ie))

        ENDDO
      ENDDO
!$OMP END DO
    ENDIF

    ! Update horizontal velocity field
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF ((itime_scheme >= 4) .AND. istep == 2) THEN
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)+ dtime              &
            & *(wgt_nnow*p_nh%diag%ddt_vn_adv(je,jk,jb,ntl1)                                &
            & + wgt_nnew*p_nh%diag%ddt_vn_adv(je,jk,jb,ntl2)+p_nh%diag%ddt_vn_phy(je,jk,jb) &
            & -cpd*z_theta_v_e(je,jk,jb)*z_gradh_exner(je,jk,jb))
          ENDDO
        ENDDO
      ELSE
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)+ dtime     &
            & *(p_nh%diag%ddt_vn_adv(je,jk,jb,ntl1)+p_nh%diag%ddt_vn_phy(je,jk,jb) &
            & -cpd*z_theta_v_e(je,jk,jb)*z_gradh_exner(je,jk,jb))
          ENDDO
        ENDDO
      ENDIF

      IF (lhdiff_rcf .AND. istep == 2) THEN
        ! apply divergence damping if diffusion is not called every sound-wave time step
        IF (divdamp_order == 2) THEN ! standard second-order divergence damping
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)  &
                + scal_divdamp*z_graddiv_vn(je,jk,jb)
            ENDDO
          ENDDO
        ELSE IF (divdamp_order == 4 .AND. (l_limited_area .OR. p_patch%id > 1)) THEN 
          ! fourth-order divergence damping with reduced damping coefficient along nest boundary
          ! (scal_divdamp is negative whereas bdy_divdamp is positive; decreasing the divergence
          ! damping along nest boundaries is beneficial because this reduces the interference
          ! with the increased diffusion applied in nh_diffusion)
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)                       &
                + (scal_divdamp+bdy_divdamp*p_int%nudgecoeff_e(je,jb))*z_graddiv2_vn(je,jk,jb)
            ENDDO
          ENDDO
        ELSE IF (divdamp_order == 4) THEN ! fourth-order divergence damping
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
              p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)  &
                + scal_divdamp*z_graddiv2_vn(je,jk,jb)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      ! Classic Rayleigh damping mechanism for vn (requires reference state !!)
      !
      IF ( rayleigh_type == RAYLEIGH_CLASSIC ) THEN
        DO jk = 1, nrdmax(p_patch%id)
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)       &
              &                          - dtime*p_nh%metrics%rayleigh_vn(jk) &
              &                          * (p_nh%prog(nnew)%vn(je,jk,jb)      &
              &                          - p_nh%ref%vn_ref(je,jk,jb))
          ENDDO
        ENDDO
      ENDIF
    ENDDO
!$OMP END DO

    IF (istep == 2 .AND. l_bdy_nudge) THEN ! apply boundary nudging if requested
!$OMP DO PRIVATE(jb,jk,je,ic) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
      DO ic = 1, p_nh%metrics%nudge_e_dim
        je = p_nh%metrics%nudge_e_idx(ic)
        jb = p_nh%metrics%nudge_e_blk(ic)
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, p_nh%metrics%nudge_e_dim
          je = p_nh%metrics%nudge_e_idx(ic)
          jb = p_nh%metrics%nudge_e_blk(ic)
#endif
          p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnew)%vn(je,jk,jb)  &
            + p_int%nudgecoeff_e(je,jb)*p_nh%diag%grf_tend_vn(je,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF


    ! Boundary update of horizontal velocity
    IF (istep == 1 .AND. (l_limited_area .OR. p_patch%id > 1)) THEN
      rl_start = 1
      rl_end   = grf_bdywidth_e

      i_startblk = p_patch%edges%start_blk(rl_start,1)
      i_endblk   = p_patch%edges%end_blk(rl_end,1)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb) + &
              dtime*p_nh%diag%grf_tend_vn(je,jk,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF

    IF (itype_comm == 2) THEN
      IF (timers_level > 5) THEN
!$OMP MASTER
        CALL timer_stop(timer_solve_nh_p1)
        CALL timer_start(timer_solve_nh_exch)
!$OMP END MASTER
      ENDIF
      ! use OpenMP-parallelized communication using global memory for buffers
      IF (istep == 1) THEN
        IF (idiv_method == 1) THEN
          CALL sync_patch_array_gm(SYNC_E,p_patch,2,bufr%send_e2,bufr%recv_e2,&
               p_nh%prog(nnew)%vn,z_rho_e)
        ELSE
          CALL sync_patch_array_gm(SYNC_E,p_patch,3,bufr%send_e3,bufr%recv_e3, &
               p_nh%prog(nnew)%vn,z_rho_e,z_theta_v_e)
        ENDIF
      ELSE
        CALL sync_patch_array_gm(SYNC_E,p_patch,1,bufr%send_e1,bufr%recv_e1,p_nh%prog(nnew)%vn)
      ENDIF
      IF (timers_level > 5) THEN
!$OMP MASTER
        CALL timer_stop(timer_solve_nh_exch)
        CALL timer_start(timer_solve_nh_p2)
!$OMP END MASTER
      ENDIF
    ENDIF

!$OMP END PARALLEL

    !-------------------------
    ! communication phase
    IF (use_icon_comm) THEN 
      IF (timers_level > 5) THEN
        CALL timer_stop(timer_solve_nh_p1)
        CALL timer_start(timer_solve_nh_exch)
      ENDIF
      IF (istep == 1) THEN
        IF (idiv_method == 1) THEN
          CALL icon_comm_sync(p_nh%prog(nnew)%vn, z_rho_e, p_patch%sync_edges_not_owned, &
            & name="solve_step1_vn")
        ELSE
          CALL icon_comm_sync(p_nh%prog(nnew)%vn, z_rho_e, z_theta_v_e, &
            & p_patch%sync_edges_not_owned, &
            & name="solve_step1_vn")
        ENDIF
      ELSE
        CALL icon_comm_sync(p_nh%prog(nnew)%vn, p_patch%sync_edges_not_owned, &
            & name="solve_step2_vn")
      ENDIF
      IF (timers_level > 5) THEN
        CALL timer_stop(timer_solve_nh_exch)
        CALL timer_start(timer_solve_nh_p2)
      ENDIF
    ELSE
      IF (itype_comm == 1) THEN
        IF (timers_level > 5) THEN
          CALL timer_stop(timer_solve_nh_p1)
          CALL timer_start(timer_solve_nh_exch)
        ENDIF
        IF (istep == 1) THEN
          IF (idiv_method == 1) THEN
            CALL sync_patch_array_mult(SYNC_E,p_patch,2,p_nh%prog(nnew)%vn,z_rho_e)
          ELSE
            CALL sync_patch_array_mult(SYNC_E,p_patch,3,p_nh%prog(nnew)%vn,z_rho_e,z_theta_v_e)
        ENDIF
        ELSE
          CALL sync_patch_array(SYNC_E,p_patch,p_nh%prog(nnew)%vn)
        ENDIF
        IF (timers_level > 5) THEN
          CALL timer_stop(timer_solve_nh_exch)
          CALL timer_start(timer_solve_nh_p2)
        ENDIF
      ENDIF
    ENDIF
    ! end communication phase
    !-------------------------

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    rl_start = 2
    rl_end   = min_rledge_int - 2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (idiv_method == 1 .AND. istep == 1) THEN

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Average normal wind components in order to get nearly second-order accurate divergence
            z_vn_avg(je,jk,jb) = p_int%e_flx_avg(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)        &
              + p_int%e_flx_avg(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%e_flx_avg(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%e_flx_avg(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%e_flx_avg(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

            ! Compute gradient of divergence of vn for divergence damping
            z_graddiv_vn(je,jk,jb) = p_int%geofac_grdiv(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)    &
              + p_int%geofac_grdiv(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%geofac_grdiv(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%geofac_grdiv(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%geofac_grdiv(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

            ! RBF reconstruction of tangential wind component
            p_nh%diag%vt(je,jk,jb) = p_int%rbf_vec_coeff_e(1,je,jb)  &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%rbf_vec_coeff_e(2,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%rbf_vec_coeff_e(3,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%rbf_vec_coeff_e(4,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

           ENDDO
        ENDDO

      ELSE IF (idiv_method == 1 .AND. itime_scheme >= 5) THEN

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Average normal wind components in order to get nearly second-order accurate divergence
            z_vn_avg(je,jk,jb) = p_int%e_flx_avg(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)        &
              + p_int%e_flx_avg(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%e_flx_avg(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%e_flx_avg(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%e_flx_avg(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

            ! RBF reconstruction of tangential wind component
            p_nh%diag%vt(je,jk,jb) = p_int%rbf_vec_coeff_e(1,je,jb)  &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%rbf_vec_coeff_e(2,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%rbf_vec_coeff_e(3,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%rbf_vec_coeff_e(4,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

           ENDDO
        ENDDO

      ELSE IF (idiv_method == 1) THEN

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Average normal wind components in order to get nearly second-order accurate divergence
            z_vn_avg(je,jk,jb) = p_int%e_flx_avg(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)        &
              + p_int%e_flx_avg(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%e_flx_avg(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%e_flx_avg(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%e_flx_avg(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))
           
           ENDDO
        ENDDO

      ELSE IF (istep == 1) THEN ! idiv_method = 2

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Compute gradient of divergence of vn for divergence damping
            z_graddiv_vn(je,jk,jb) = p_int%geofac_grdiv(je,1,jb)*p_nh%prog(nnew)%vn(je,jk,jb)    &
              + p_int%geofac_grdiv(je,2,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%geofac_grdiv(je,3,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%geofac_grdiv(je,4,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%geofac_grdiv(je,5,jb)*p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

            ! RBF reconstruction of tangential wind component
            p_nh%diag%vt(je,jk,jb) = p_int%rbf_vec_coeff_e(1,je,jb)  &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%rbf_vec_coeff_e(2,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%rbf_vec_coeff_e(3,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%rbf_vec_coeff_e(4,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

           ENDDO
        ENDDO

      ELSE IF (itime_scheme >= 5) THEN

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! RBF reconstruction of tangential wind component
            p_nh%diag%vt(je,jk,jb) = p_int%rbf_vec_coeff_e(1,je,jb)  &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%rbf_vec_coeff_e(2,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%rbf_vec_coeff_e(3,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%rbf_vec_coeff_e(4,je,jb)                       &
              * p_nh%prog(nnew)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

           ENDDO
        ENDDO

      ENDIF

      IF (istep == 1 .OR. itime_scheme >= 5) THEN
        ! Compute contravariant correction for vertical velocity at full levels
        DO jk = nflatlev(p_patch%id), nlev
          DO je = i_startidx, i_endidx
            z_w_concorr_me(je,jk,jb) =                                          &
              p_nh%prog(nnew)%vn(je,jk,jb)*p_nh%metrics%ddxn_z_full(je,jk,jb) + &
              p_nh%diag%vt(je,jk,jb)      *p_nh%metrics%ddxt_z_full(je,jk,jb)
          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO

    IF (lhdiff_rcf .AND. istep == 1 .AND. divdamp_order == 4) THEN ! fourth-order divergence damping
      rl_start = grf_bdywidth_e + 1
      rl_end   = min_rledge_int

      i_startblk = p_patch%edges%start_blk(rl_start,1)
      i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif

            ! Compute gradient of divergence of gradient of divergence for fourth-order divergence damping
            z_graddiv2_vn(je,jk,jb) = p_int%geofac_grdiv(je,1,jb)*z_graddiv_vn(je,jk,jb)   &
              + p_int%geofac_grdiv(je,2,jb)*z_graddiv_vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) &
              + p_int%geofac_grdiv(je,3,jb)*z_graddiv_vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) &
              + p_int%geofac_grdiv(je,4,jb)*z_graddiv_vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) &
              + p_int%geofac_grdiv(je,5,jb)*z_graddiv_vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

           ENDDO
        ENDDO
      ENDDO
!$OMP END DO

    ENDIF

    ! It turned out that it is sufficient to compute the contravariant correction in the
    ! predictor step at time level n+1; repeating the calculation in the corrector step
    ! has negligible impact on the results except in very-high resolution runs with extremely steep mountains
    IF (istep == 1 .OR. itime_scheme >= 5) THEN

      rl_start = 3
      rl_end = min_rlcell_int - 1

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,z_w_concorr_mc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Interpolate contravariant correction to cell centers...
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
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
        DO jk = nflatlev(p_patch%id)+1, nlev
          DO jc = i_startidx, i_endidx
            p_nh%diag%w_concorr_c(jc,jk,jb) =                                &
              p_nh%metrics%wgtfac_c(jc,jk,jb)*z_w_concorr_mc(jc,jk) +        &
             (1._wp - p_nh%metrics%wgtfac_c(jc,jk,jb))*z_w_concorr_mc(jc,jk-1) 
          ENDDO
        ENDDO

        DO jc = i_startidx, i_endidx
          p_nh%diag%w_concorr_c(jc,nlevp1,jb) =                         &
            p_nh%metrics%wgtfacq_c(jc,1,jb)*z_w_concorr_mc(jc,nlev) +   &
            p_nh%metrics%wgtfacq_c(jc,2,jb)*z_w_concorr_mc(jc,nlev-1) + &
            p_nh%metrics%wgtfacq_c(jc,3,jb)*z_w_concorr_mc(jc,nlev-2)
        ENDDO

      ENDDO
!$OMP END DO
    ENDIF

    rl_start = 7
    IF (idiv_method == 1) THEN
      rl_end = min_rledge_int - 2
    ELSE
      rl_end = min_rledge_int - 3
    ENDIF

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

    IF (idiv_method == 2 .AND. (p_patch%id > 1 .OR. l_limited_area)) THEN
!$OMP WORKSHARE
      z_theta_v_fl_e(:,:,p_patch%edges%start_blk(5,1):i_startblk) = 0._wp
!$OMP END WORKSHARE
    ENDIF

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Fluxes at edges
      DO jk = 1,nlev
        DO je = i_startidx, i_endidx

          p_nh%diag%mass_fl_e(je,jk,jb) = z_rho_e(je,jk,jb)         &
            * ptr_vn(je,jk,jb) * p_nh%metrics%ddqz_z_full_e(je,jk,jb)

          z_theta_v_fl_e(je,jk,jb)= p_nh%diag%mass_fl_e(je,jk,jb)   &
            * z_theta_v_e(je,jk,jb)

        ENDDO
      ENDDO

      IF (lsave_mflx) THEN
        p_nh%diag%mass_fl_e_sv(i_startidx:i_endidx,:,jb) = p_nh%diag%mass_fl_e(i_startidx:i_endidx,:,jb)
      ENDIF

    ENDDO
!$OMP END DO

    IF (p_patch%id > 1 .AND. grf_intmethod_e >= 5 .AND. idiv_method == 1) THEN
!$OMP DO PRIVATE(ic,je,jb,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO ic = 1, p_nh%metrics%bdy_mflx_e_dim
        je = p_nh%metrics%bdy_mflx_e_idx(ic)
        jb = p_nh%metrics%bdy_mflx_e_blk(ic)

        IF (jstep == 0) THEN
          DO jk = 1, nlev
            p_nh%diag%grf_bdy_mflx(jk,ic,2) = p_nh%diag%grf_tend_mflx(je,jk,jb)
            p_nh%diag%grf_bdy_mflx(jk,ic,1) = mflx_avg(je,jk,jb) - dt_shift*p_nh%diag%grf_bdy_mflx(jk,ic,2)
          ENDDO
        ENDIF

        DO jk = 1, nlev
          p_nh%diag%mass_fl_e(je,jk,jb) = p_nh%diag%grf_bdy_mflx(jk,ic,1) + &
            REAL(jstep,wp)*dtime*p_nh%diag%grf_bdy_mflx(jk,ic,2)
          z_theta_v_fl_e(je,jk,jb) = p_nh%diag%mass_fl_e(je,jk,jb) * z_theta_v_e(je,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
    ENDIF

!$OMP END PARALLEL


    IF (idiv_method == 1) THEN ! use simple divergence based on averaged velocity

      ! horizontal divergences of rho and rhotheta are processed in one step for efficiency
      CALL div(p_nh%diag%mass_fl_e, p_patch, p_int, z_mass_fl_div, opt_in2=z_theta_v_fl_e, &
               opt_out2=z_theta_v_fl_div, opt_rlstart=4, opt_rlend=min_rlcell_int)

    ELSE ! use averaged divergence

      ! horizontal divergences of rho and rhotheta are processed in one step for efficiency
      CALL div_avg(p_nh%diag%mass_fl_e, p_patch, p_int, p_int%c_bln_avg, z_mass_fl_div, &
                   opt_in2=z_theta_v_fl_e, opt_out2=z_theta_v_fl_div, opt_rlstart=4,    &
                   opt_rlend=min_rlcell_int)
    ENDIF

    ! Vertical solution:
!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    IF (l_vert_nested) THEN
      jk_start = 2
    ELSE
      jk_start = 1
    ENDIF

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,z_w_expl,z_contr_w_fl_l,z_rho_expl,z_exner_expl,   &
!$OMP            z_a,z_b,z_c,z_g,z_q,z_alpha,z_beta,z_gamma,ic,z_raylfac,z_w_lim) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! upper boundary conditions for rho_ic and theta_v_ic in the case of vertical nesting
      IF (l_vert_nested .AND. istep == 1) THEN
        DO jc = i_startidx, i_endidx
          p_nh%diag%theta_v_ic(jc,1,jb) = p_nh%diag%theta_v_ic(jc,2,jb) + &
            p_nh%diag%dtheta_v_ic_ubc(jc,jb)
          z_mflx_top(jc,jb)             = p_nh%diag%mflx_ic_ubc(jc,jb,1) + &
            REAL(jstep,wp)*dtime*p_nh%diag%mflx_ic_ubc(jc,jb,2)
          p_nh%diag%rho_ic(jc,1,jb) =  0._wp ! not used in dynamical core in this case, will be set for tracer interface later
        ENDDO
      ELSE IF (l_vert_nested .AND. istep == 2) THEN
        DO jc = i_startidx, i_endidx
          p_nh%diag%theta_v_ic(jc,1,jb) = p_nh%diag%theta_v_ic(jc,2,jb) + &
            p_nh%diag%dtheta_v_ic_ubc(jc,jb)
          z_mflx_top(jc,jb)             = p_nh%diag%mflx_ic_ubc(jc,jb,1) + &
            (REAL(jstep,wp)+0.5_wp)*dtime*p_nh%diag%mflx_ic_ubc(jc,jb,2)
        ENDDO
      ENDIF

      IF (istep == 2 .AND. (itime_scheme >= 4)) THEN
!CDIR UNROLL=5
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx

            ! explicit part for w
            z_w_expl(jc,jk) = p_nh%prog(nnow)%w(jc,jk,jb)  + dtime * &
              (wgt_nnow*p_nh%diag%ddt_w_adv(jc,jk,jb,ntl1) +         &
               wgt_nnew*p_nh%diag%ddt_w_adv(jc,jk,jb,ntl2)           &
              -cpd*z_th_ddz_exner_c(jc,jk,jb) )

            ! contravariant vertical velocity times density for explicit part
            z_contr_w_fl_l(jc,jk) = p_nh%diag%rho_ic(jc,jk,jb)*(-p_nh%diag%w_concorr_c(jc,jk,jb) &
              + p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,jk,jb) )

          ENDDO
        ENDDO
      ELSE
!CDIR UNROLL=5
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx

            ! explicit part for w
            z_w_expl(jc,jk) = p_nh%prog(nnow)%w(jc,jk,jb) + dtime *            &
              (p_nh%diag%ddt_w_adv(jc,jk,jb,ntl1)-cpd*z_th_ddz_exner_c(jc,jk,jb))

            ! contravariant vertical velocity times density for explicit part
            z_contr_w_fl_l(jc,jk) = p_nh%diag%rho_ic(jc,jk,jb)*(-p_nh%diag%w_concorr_c(jc,jk,jb) &
              + p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,jk,jb) )

          ENDDO
        ENDDO
      ENDIF

      ! Solver coefficients
!CDIR UNROLL=3
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          z_beta(jc,jk)=dtime*rd*p_nh%prog(nnow)%exner(jc,jk,jb)         &
          &                   /(cvd*p_nh%prog(nnow)%rhotheta_v(jc,jk,jb) &
          &                   *p_nh%metrics%ddqz_z_full(jc,jk,jb))

          z_alpha(jc,jk)= p_nh%metrics%vwind_impl_wgt(jc,jb)*         &
            &  p_nh%diag%theta_v_ic(jc,jk,jb)*p_nh%diag%rho_ic(jc,jk,jb)
        ENDDO
      ENDDO

      z_alpha(:,nlevp1) = 0.0_wp
!CDIR UNROLL=2
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          z_gamma(jc,jk) =  dtime*cpd*p_nh%metrics%vwind_impl_wgt(jc,jb)* &
            p_nh%diag%theta_v_ic(jc,jk,jb)/p_nh%metrics%ddqz_z_half(jc,jk,jb)

          z_a(jc,jk) = -z_gamma(jc,jk)*z_beta(jc,jk-1)*z_alpha(jc,jk-1)
          z_c(jc,jk) = -z_gamma(jc,jk)*z_beta(jc,jk  )*z_alpha(jc,jk+1)
          z_b(jc,jk) = 1.0_wp+z_gamma(jc,jk)*z_alpha(jc,jk) &
            *(z_beta(jc,jk-1)+z_beta(jc,jk))
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        z_q(jc,2) = -z_c(jc,2)/z_b(jc,2)
      ENDDO

!CDIR UNROLL=4
      DO jk = 3, nlev
        DO jc = i_startidx, i_endidx
          z_g(jc,jk) = 1.0_wp/(z_b(jc,jk)+z_a(jc,jk)*z_q(jc,jk-1))
          z_q(jc,jk) = - z_c(jc,jk)*z_g(jc,jk)
        ENDDO
      ENDDO

      ! upper boundary condition for w (interpolated from parent domain in case of vertical nesting)
      ! Note: the upper b.c. reduces to w(1) = 0 in the absence of diabatic heating
      IF (l_open_ubc .AND. .NOT. l_vert_nested) THEN
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,1,jb) = z_thermal_exp(jc,jb)
          z_contr_w_fl_l(jc,1) = p_nh%diag%rho_ic(jc,1,jb)*p_nh%prog(nnow)%w(jc,1,jb)   &
            * p_nh%metrics%vwind_expl_wgt(jc,jb)
        ENDDO
      ELSE IF (.NOT. l_open_ubc .AND. .NOT. l_vert_nested) THEN
        p_nh%prog(nnew)%w(:,1,jb) = 0._wp
        z_contr_w_fl_l(:,1)       = 0._wp
      ELSE  ! l_vert_nested
        DO jc = i_startidx, i_endidx
          z_contr_w_fl_l(jc,1) = z_mflx_top(jc,jb) * p_nh%metrics%vwind_expl_wgt(jc,jb)
        ENDDO
      ENDIF

      ! lower boundary condition for w, consistent with contravariant correction
      DO jc = i_startidx, i_endidx
        p_nh%prog(nnew)%w(jc,nlevp1,jb) = p_nh%diag%w_concorr_c(jc,nlevp1,jb)
        z_contr_w_fl_l(jc,nlevp1)       = 0.0_wp
      ENDDO


      ! other full level stuff
      ! Top level first
      DO jc = i_startidx, i_endidx
        z_rho_expl(jc,1)=        p_nh%prog(nnow)%rho(jc,1,jb) &
        &        -dtime*p_nh%metrics%inv_ddqz_z_full(jc,1,jb) &
        &                            *(z_mass_fl_div(jc,1,jb) &
        &                            +z_contr_w_fl_l(jc,1   ) &
        &                            -z_contr_w_fl_l(jc,2   ))

        z_exner_expl(jc,1)=             z_exner_pr(jc,1,jb)                    &
        &      -z_beta (jc,1)*(z_theta_v_fl_div(jc,1,jb)                       &
        & +p_nh%diag%theta_v_ic(jc,1,jb)*z_contr_w_fl_l(jc,1)                   &
        & -p_nh%diag%theta_v_ic(jc,2,jb)*z_contr_w_fl_l(jc,2))                  &
        & +dtime*(p_nh%diag%ddt_exner(jc,1,jb)+p_nh%diag%ddt_exner_phy(jc,1,jb))
      ENDDO

      ! Other levels
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          z_rho_expl(jc,jk)=       p_nh%prog(nnow)%rho(jc,jk  ,jb) &
          &        -dtime*p_nh%metrics%inv_ddqz_z_full(jc,jk  ,jb) &
          &                            *(z_mass_fl_div(jc,jk  ,jb) &
          &                            +z_contr_w_fl_l(jc,jk     ) &
          &                             -z_contr_w_fl_l(jc,jk+1   ))

          z_exner_expl(jc,jk)=          z_exner_pr(jc,jk,jb) - z_beta(jc,jk)         &
          &                             *(z_theta_v_fl_div(jc,jk,jb)                 &
          &   +p_nh%diag%theta_v_ic(jc,jk  ,jb)*z_contr_w_fl_l(jc,jk  )               &
          &   -p_nh%diag%theta_v_ic(jc,jk+1,jb)*z_contr_w_fl_l(jc,jk+1))              &
          &   +dtime*(p_nh%diag%ddt_exner(jc,jk,jb)+p_nh%diag%ddt_exner_phy(jc,jk,jb))

          p_nh%prog(nnew)%w(jc,jk,jb) = z_w_expl(jc,jk) - z_gamma(jc,jk)  &
          &      *(z_exner_expl(jc,jk-1)-z_exner_expl(jc,jk))
        ENDDO
      ENDDO

      ! Solve tridiagonal matrix for w
      DO jc = i_startidx, i_endidx
        p_nh%prog(nnew)%w(jc,2,jb)= p_nh%prog(nnew)%w(jc,2,jb)/z_b(jc,2)
      ENDDO

      DO jk = 3, nlev
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,jk,jb) = (p_nh%prog(nnew)%w(jc,jk,jb)  &
            -z_a(jc,jk)*p_nh%prog(nnew)%w(jc,jk-1,jb))*z_g(jc,jk)
        ENDDO
      ENDDO

      DO jk = nlev-1, 2, -1
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnew)%w(jc,jk,jb)&
          &             +p_nh%prog(nnew)%w(jc,jk+1,jb)*z_q(jc,jk)
        ENDDO
      ENDDO

      IF (l_vert_nested) THEN
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,1,jb) = p_nh%prog(nnew)%w(jc,2,jb) + p_nh%diag%dw_ubc(jc,jb)
        ENDDO
      ENDIF

      ! Rayleigh damping mechanism (Klemp,Dudhia,Hassiotis:MWR136,pp.3987-4004)
      !
      IF ( rayleigh_type == RAYLEIGH_KLEMP ) THEN
        DO jk = 2, nrdmax(p_patch%id)
          z_raylfac = 1.0_wp/(1.0_wp+dtime*p_nh%metrics%rayleigh_w(jk))
          DO jc = i_startidx, i_endidx
            p_nh%prog(nnew)%w(jc,jk,jb) = z_raylfac*p_nh%prog(nnew)%w(jc,jk,jb) +    &
                                          (1._wp-z_raylfac)*p_nh%prog(nnew)%w(jc,1,jb)
          ENDDO
        ENDDO
      ! Classic Rayleigh damping mechanism for w (requires reference state !!)
      !
      ELSE IF ( rayleigh_type == RAYLEIGH_CLASSIC ) THEN
        DO jk = 2, nrdmax(p_patch%id)
          DO jc = i_startidx, i_endidx
            p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnew)%w(jc,jk,jb)       &
              &                         - dtime*p_nh%metrics%rayleigh_w(jk) &
              &                         * ( p_nh%prog(nnew)%w(jc,jk,jb)     &
              &                         - p_nh%ref%w_ref(jc,jk,jb) )
          ENDDO
        ENDDO
      ENDIF

      ! Results
      DO jk = jk_start, nlev
        DO jc = i_startidx, i_endidx

          ! density
          p_nh%prog(nnew)%rho(jc,jk,jb) = z_rho_expl(jc,jk)              &
            - p_nh%metrics%vwind_impl_wgt(jc,jb)*dtime                   &
            * p_nh%metrics%inv_ddqz_z_full(jc,jk,jb)                     &
            *(p_nh%diag%rho_ic(jc,jk  ,jb)*p_nh%prog(nnew)%w(jc,jk  ,jb) &
            - p_nh%diag%rho_ic(jc,jk+1,jb)*p_nh%prog(nnew)%w(jc,jk+1,jb))

          ! exner
          p_nh%prog(nnew)%exner(jc,jk,jb) = z_exner_expl(jc,jk) &
            + p_nh%metrics%exner_ref_mc(jc,jk,jb)-z_beta(jc,jk) &
            *(z_alpha(jc,jk  )*p_nh%prog(nnew)%w(jc,jk  ,jb)    &
            - z_alpha(jc,jk+1)*p_nh%prog(nnew)%w(jc,jk+1,jb))

          ! rho*theta
          p_nh%prog(nnew)%rhotheta_v(jc,jk,jb) = p_nh%prog(nnow)%rhotheta_v(jc,jk,jb)   &
            *( (p_nh%prog(nnew)%exner(jc,jk,jb)/p_nh%prog(nnow)%exner(jc,jk,jb)-1.0_wp) &
            *   cvd_o_rd+1.0_wp)

          ! theta
          p_nh%prog(nnew)%theta_v(jc,jk,jb) = &
            p_nh%prog(nnew)%rhotheta_v(jc,jk,jb)/p_nh%prog(nnew)%rho(jc,jk,jb)

        ENDDO
      ENDDO

      ! Special treatment of uppermost layer in the case of vertical nesting
      IF (l_vert_nested) THEN
        DO jc = i_startidx, i_endidx

          ! density
          p_nh%prog(nnew)%rho(jc,1,jb) = z_rho_expl(jc,1)                             &
            - p_nh%metrics%vwind_impl_wgt(jc,jb)*dtime                                &
            * p_nh%metrics%inv_ddqz_z_full(jc,1,jb)                                   &
            *(z_mflx_top(jc,jb) - p_nh%diag%rho_ic(jc,2,jb)*p_nh%prog(nnew)%w(jc,2,jb))

          ! exner
          p_nh%prog(nnew)%exner(jc,1,jb) = z_exner_expl(jc,1)                  &
            + p_nh%metrics%exner_ref_mc(jc,1,jb)-z_beta(jc,1)                  &
            *(p_nh%metrics%vwind_impl_wgt(jc,jb)*p_nh%diag%theta_v_ic(jc,1,jb) &
            * z_mflx_top(jc,jb) - z_alpha(jc,2)*p_nh%prog(nnew)%w(jc,2,jb))

          ! rho*theta
          p_nh%prog(nnew)%rhotheta_v(jc,1,jb) = p_nh%prog(nnow)%rhotheta_v(jc,1,jb)   &
            *( (p_nh%prog(nnew)%exner(jc,1,jb)/p_nh%prog(nnow)%exner(jc,1,jb)-1.0_wp) &
            *   cvd_o_rd+1.0_wp)

          ! theta
          p_nh%prog(nnew)%theta_v(jc,1,jb) = &
            p_nh%prog(nnew)%rhotheta_v(jc,1,jb)/p_nh%prog(nnew)%rho(jc,1,jb)

        ENDDO
      ENDIF

      IF (istep == 2 .AND. l_vert_nested) THEN
        ! Rediagnose appropriate rho_ic(jk=1) for tracer transport
        DO jc = i_startidx, i_endidx
          z_w_lim = p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,1,jb) + &
                    p_nh%metrics%vwind_impl_wgt(jc,jb)*p_nh%prog(nnew)%w(jc,1,jb)
          z_w_lim = SIGN(MAX(1.e-6_wp,ABS(z_w_lim)),z_w_lim)
          p_nh%diag%rho_ic(jc,1,jb) = z_mflx_top(jc,jb)/z_w_lim
        ENDDO
      ENDIF

      IF (istep == 2) THEN
        ! store dynamical part of exner time increment in exner_dyn_incr
        ! the conversion into a temperature tendency is done in the NWP interface
        DO jk = kstart_moist(p_patch%id), nlev
          DO jc = i_startidx, i_endidx
            p_nh%diag%exner_dyn_incr(jc,jk,jb) = p_nh%diag%exner_dyn_incr(jc,jk,jb) + &
              p_nh%prog(nnew)%exner(jc,jk,jb) - p_nh%diag%exner_old(jc,jk,jb) -   &
              dtime*p_nh%diag%ddt_exner_phy(jc,jk,jb)
          ENDDO
        ENDDO
      ENDIF

      IF (istep == 2 .AND. l_child_vertnest) THEN
        ! Store values at nest interface levels
        DO jc = i_startidx, i_endidx

          p_nh%diag%dw_int(jc,jb,idyn_timestep) =                                         &
            0.5_wp*(p_nh%prog(nnow)%w(jc,nshift,jb)   + p_nh%prog(nnew)%w(jc,nshift,jb) - &
                   (p_nh%prog(nnow)%w(jc,nshift+1,jb) + p_nh%prog(nnew)%w(jc,nshift+1,jb)))

          p_nh%diag%mflx_ic_int(jc,jb,idyn_timestep) = p_nh%diag%rho_ic(jc,nshift,jb) * &
            (p_nh%metrics%vwind_expl_wgt(jc,jb)*p_nh%prog(nnow)%w(jc,nshift,jb) + &
             p_nh%metrics%vwind_impl_wgt(jc,jb)*p_nh%prog(nnew)%w(jc,nshift,jb))

          p_nh%diag%dtheta_v_ic_int(jc,jb,idyn_timestep) = p_nh%diag%theta_v_ic(jc,nshift,jb) - &
            p_nh%diag%theta_v_ic(jc,nshift+1,jb)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO

    ! Boundary update in case of nesting
    IF (istep == 1 .AND. (l_limited_area .OR. p_patch%id > 1) &
      & .AND. my_process_is_mpi_all_seq() ) THEN

      rl_start = 1
      rl_end   = grf_bdywidth_c

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

            p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb) + &
              dtime*p_nh%diag%grf_tend_rho(jc,jk,jb)

            p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnow)%theta_v(jc,jk,jb) + &
              dtime*p_nh%diag%grf_tend_thv(jc,jk,jb)

            ! Diagnose rhotheta from rho and theta
            p_nh%prog(nnew)%rhotheta_v(jc,jk,jb) = p_nh%prog(nnew)%rho(jc,jk,jb) * &
              p_nh%prog(nnew)%theta_v(jc,jk,jb)

            ! Diagnose exner from rhotheta
            p_nh%prog(nnew)%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
              p_nh%prog(nnew)%rhotheta_v(jc,jk,jb)))

            p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnow)%w(jc,jk,jb) + &
              dtime*p_nh%diag%grf_tend_w(jc,jk,jb)

          ENDDO
        ENDDO

        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,nlevp1,jb) = p_nh%prog(nnow)%w(jc,nlevp1,jb) + &
            dtime*p_nh%diag%grf_tend_w(jc,nlevp1,jb)
          p_nh%diag%rho_ic(jc,1,jb) =  0.5_wp*(                              &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*p_nh%prog(nnow)%rho(jc,1,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*p_nh%prog(nnow)%rho(jc,2,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*p_nh%prog(nnow)%rho(jc,3,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*p_nh%prog(nnew)%rho(jc,1,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*p_nh%prog(nnew)%rho(jc,2,jb) +  &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*p_nh%prog(nnew)%rho(jc,3,jb) )
        ENDDO
      ENDDO
!OMP END DO

    ELSE IF (istep == 1 .AND. (l_limited_area .OR. p_patch%id > 1)) THEN
      ! In the MPI-parallelized case, only rho and w are updated here,
      ! and theta_v is preliminarily stored on exner in order to save
      ! halo communications

      rl_start = 1
      rl_end   = grf_bdywidth_c

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

            p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb) + &
              dtime*p_nh%diag%grf_tend_rho(jc,jk,jb)

            ! *** Storing theta_v on exner is done to save MPI communications ***
            ! DO NOT TOUCH THIS!
            p_nh%prog(nnew)%exner(jc,jk,jb) = p_nh%prog(nnow)%theta_v(jc,jk,jb) + &
              dtime*p_nh%diag%grf_tend_thv(jc,jk,jb)

            p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnow)%w(jc,jk,jb) + &
              dtime*p_nh%diag%grf_tend_w(jc,jk,jb)

          ENDDO
        ENDDO

        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,nlevp1,jb) = p_nh%prog(nnow)%w(jc,nlevp1,jb) + &
            dtime*p_nh%diag%grf_tend_w(jc,nlevp1,jb)
        ENDDO
      ENDDO
!OMP END DO

    ENDIF

    IF (itype_comm == 2) THEN
      IF (timers_level > 5) THEN
!$OMP MASTER
        CALL timer_stop(timer_solve_nh_p2)
        CALL timer_start(timer_solve_nh_exch)
!$OMP END MASTER
      ENDIF
      ! use OpenMP-parallelized communication using global memory for buffers
      IF (istep == 1) THEN ! Only w is updated in the predictor step
        CALL sync_patch_array_gm(SYNC_C,p_patch,1,bufr%send_c1,bufr%recv_c1,p_nh%prog(nnew)%w)
      ELSE IF (istep == 2) THEN
        ! Synchronize all prognostic variables
        CALL sync_patch_array_gm(SYNC_C,p_patch,3,bufr%send_c3,bufr%recv_c3,                &
                                 p_nh%prog(nnew)%rho,p_nh%prog(nnew)%exner,p_nh%prog(nnew)%w)
      ENDIF
      IF (timers_level > 5) THEN
!$OMP MASTER
        CALL timer_stop(timer_solve_nh_exch)
        IF (istep == 1) CALL timer_start(timer_solve_nh_p1)
!$OMP END MASTER
      ENDIF
    ENDIF

!$OMP END PARALLEL

    !-------------------------
    ! communication phase
    IF (use_icon_comm) THEN 
      IF (timers_level > 5) THEN
        CALL timer_stop(timer_solve_nh_p2)
        CALL timer_start(timer_solve_nh_exch)
      ENDIF
      IF (istep == 1) THEN ! Only w is updated in the predictor step
        CALL icon_comm_sync(p_nh%prog(nnew)%w, p_patch%sync_cells_not_owned, &
            & name="solve_step1_w")
      ELSE IF (istep == 2) THEN
        ! Synchronize all prognostic variables
        CALL icon_comm_sync(p_nh%prog(nnew)%rho, p_nh%prog(nnew)%exner, p_nh%prog(nnew)%w, &
          & p_patch%sync_cells_not_owned, &
          & name="solve_step2_w")
      ENDIF
      IF (timers_level > 5) THEN
        CALL timer_stop(timer_solve_nh_exch)
        IF (istep == 1) CALL timer_start(timer_solve_nh_p1)
      ENDIF
    
    ELSE
      IF (itype_comm == 1) THEN
        IF (timers_level > 5) THEN
          CALL timer_stop(timer_solve_nh_p2)
          CALL timer_start(timer_solve_nh_exch)
        ENDIF
        IF (istep == 1) THEN ! Only w is updated in the predictor step
          CALL sync_patch_array(SYNC_C,p_patch,p_nh%prog(nnew)%w)
        ELSE IF (istep == 2) THEN
          ! Synchronize all prognostic variables
          CALL sync_patch_array_mult(SYNC_C,p_patch,3,p_nh%prog(nnew)%rho,  &
                                    p_nh%prog(nnew)%exner,p_nh%prog(nnew)%w)
        ENDIF
        IF (timers_level > 5) THEN
          CALL timer_stop(timer_solve_nh_exch)
          IF (istep == 1) CALL timer_start(timer_solve_nh_p1)
        ENDIF
      ENDIF
    ENDIF
    ! end communication phase
    !-------------------------

    ENDDO ! istep-loop

    ! The remaining computations are needed for MPI-parallelized applications only
    IF (my_process_is_mpi_all_seq() ) THEN
      IF (ltimer) CALL timer_stop(timer_solve_nh)
      RETURN
    ENDIF 

! OpenMP directives are commented for the NEC because the overhead is too large
#ifndef __SX__
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
#endif
    IF (l_limited_area .OR. p_patch%id > 1) THEN

      ! Index list over halo points lying in the boundary interpolation zone
      ! Note: this list typically contains at most 10 grid points 
#ifndef __SX__
!$OMP DO PRIVATE(jb,ic,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
#endif
      DO ic = 1, p_nh%metrics%bdy_halo_c_dim

        jb = p_nh%metrics%bdy_halo_c_blk(ic)
        jc = p_nh%metrics%bdy_halo_c_idx(ic)

        DO jk = 1, nlev
          p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnew)%exner(jc,jk,jb)

          ! Diagnose rhotheta from rho and theta
          p_nh%prog(nnew)%rhotheta_v(jc,jk,jb) = p_nh%prog(nnew)%rho(jc,jk,jb) * &
            p_nh%prog(nnew)%theta_v(jc,jk,jb)

          ! Diagnose exner from rhotheta
          p_nh%prog(nnew)%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
            p_nh%prog(nnew)%rhotheta_v(jc,jk,jb)))

        ENDDO
      ENDDO
#ifndef __SX__
!$OMP END DO
#endif

      rl_start = 1
      rl_end   = grf_bdywidth_c

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,1)

#ifndef __SX__
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
#endif
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

            p_nh%prog(nnew)%theta_v(jc,jk,jb) = p_nh%prog(nnew)%exner(jc,jk,jb)

            ! Diagnose rhotheta from rho and theta
            p_nh%prog(nnew)%rhotheta_v(jc,jk,jb) = p_nh%prog(nnew)%rho(jc,jk,jb) * &
              p_nh%prog(nnew)%theta_v(jc,jk,jb)

            ! Diagnose exner from rhotheta
            p_nh%prog(nnew)%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
              p_nh%prog(nnew)%rhotheta_v(jc,jk,jb)))

          ENDDO
        ENDDO
      ENDDO
#ifndef __SX__
!$OMP END DO
#endif
    ENDIF

    rl_start = min_rlcell_int - 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

#ifndef __SX__
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
#endif
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          IF (p_nh%metrics%mask_prog_halo_c(jc,jb)) THEN

            p_nh%prog(nnew)%rhotheta_v(jc,jk,jb) = p_nh%prog(nnow)%rhotheta_v(jc,jk,jb)   &
              *( (p_nh%prog(nnew)%exner(jc,jk,jb)/p_nh%prog(nnow)%exner(jc,jk,jb)-1.0_wp) &
              *   cvd_o_rd+1.0_wp)

            p_nh%prog(nnew)%theta_v(jc,jk,jb) = &
              p_nh%prog(nnew)%rhotheta_v(jc,jk,jb)/p_nh%prog(nnew)%rho(jc,jk,jb)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
#ifndef __SX__
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

   IF (ltimer) CALL timer_stop(timer_solve_nh)

  END SUBROUTINE solve_nh

END MODULE mo_solve_nonhydro

