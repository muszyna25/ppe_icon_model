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
MODULE mo_solve_nonhydro

  USE mo_kind,                 ONLY: wp
  USE mo_nonhydrostatic_config,ONLY: itime_scheme,iadv_rhotheta, igradp_method, l_open_ubc
  USE mo_dynamics_config,      ONLY: idiv_method
  USE mo_parallel_config,    ONLY: nproma, p_test_run, itype_comm
  USE mo_run_config,         ONLY: ltimer, lvert_nest
  USE mo_model_domain,       ONLY: t_patch
  USE mo_model_domain_import,ONLY: l_limited_area
  USE mo_interpolation,      ONLY: t_int_state, cells2edges_scalar, rbf_vec_interpol_edge
  USE mo_nonhydro_state,     ONLY: t_nh_state, t_nh_metrics, t_nh_diag, t_nh_prog, &
                                  t_buffer_memory
  USE mo_physical_constants,ONLY: cpd, rd, cvd, cvd_o_rd, grav, rd_o_cpd, p0ref
  USE mo_math_operators,    ONLY: div, rot_vertex, div_avg
  USE mo_vertical_grid,     ONLY: nrdmax, nflat_gradp
  USE mo_nh_init_utils,     ONLY: nflatlev
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,    ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int, min_rlcell
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_advection_hflux,   ONLY: upwind_hflux_miura, upwind_hflux_miura3
  USE mo_sync,              ONLY: SYNC_E, SYNC_C, sync_patch_array, sync_patch_array_mult, &
                                  sync_patch_array_gm
  USE mo_mpi,               ONLY: my_process_is_mpi_all_seq
  USE mo_timer,             ONLY: timer_solve_nh, timer_start, timer_stop

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
                                  ntnd, istep, l_init)

    ! Passed variables
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_int_state), TARGET, INTENT(IN):: p_int
    TYPE(t_nh_prog), INTENT(INOUT)       :: p_prog
    TYPE(t_nh_metrics), INTENT(IN)       :: p_metrics
    TYPE(t_nh_diag), INTENT(inout)       :: p_diag

    INTEGER, INTENT(IN)  :: ntnd  ! time level of ddt_adv fields used to store tendencies
    INTEGER, INTENT(IN)  :: istep ! 1: predictor step, 2: corrector step
    LOGICAL, INTENT(IN)  :: l_init

    ! Local variables
    INTEGER :: jb, jk, jc, je
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    REAL(wp):: z_concorr_e(nproma,p_patch%nlevp1,p_patch%nblks_e)
    REAL(wp):: z_w_con_c(nproma,p_patch%nlevp1,p_patch%nblks_c)
    REAL(wp):: z_w_con_c_full(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp):: z_vt_ie(nproma,p_patch%nlevp1)
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

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx, z_vt_ie)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (istep == 1) THEN

        ! Interpolate vn to interface levels and compute horizontal part of kinetic energy on edges
        DO jk = 2, nlev
          DO je = i_startidx, i_endidx
            p_diag%vn_ie(je,jk,jb) =                                  &
              p_metrics%wgtfac_e(je,jk,jb)*p_prog%vn(je,jk,jb) +        &
             (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_prog%vn(je,jk-1,jb)
            z_kin_hor_e(je,jk,jb) = 0.5_wp*(p_prog%vn(je,jk,jb)*p_prog%vn(je,jk,jb) + &
              p_diag%vt(je,jk,jb)*p_diag%vt(je,jk,jb) )
          ENDDO
        ENDDO

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

      ELSE ! corrector step (istep = 2)

        ! Compute only horizontal kinetic energy
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            z_kin_hor_e(je,jk,jb) = 0.5_wp*                &
             (p_prog%vn(je,jk,jb)*p_prog%vn(je,jk,jb) +    &
              p_diag%vt(je,jk,jb)*p_diag%vt(je,jk,jb) )
          ENDDO
        ENDDO

      ENDIF ! istep == 1

      IF (istep == 1 .AND. l_init) THEN
        ! Interpolate vt to interface levels
        DO jk = nflatlev(p_patch%id)+1, nlev
          DO je = i_startidx, i_endidx
            z_vt_ie(je,jk) =                                     &
              p_metrics%wgtfac_e(je,jk,jb)*p_diag%vt(je,jk,jb) + &
             (1._wp - p_metrics%wgtfac_e(je,jk,jb))*p_diag%vt(je,jk-1,jb)
          ENDDO
        ENDDO

        ! Bottom level
        DO je = i_startidx, i_endidx
          z_vt_ie(je,nlevp1) =                                     &
            p_metrics%wgtfacq_e(je,1,jb)*p_diag%vt(je,nlev,jb) +   &
            p_metrics%wgtfacq_e(je,2,jb)*p_diag%vt(je,nlev-1,jb) + &
            p_metrics%wgtfacq_e(je,3,jb)*p_diag%vt(je,nlev-2,jb)
        ENDDO

        ! Compute contravariant correction for vertical velocity at interface levels
        ! (will be interpolated to cell centers below)
        DO jk = nflatlev(p_patch%id)+1, nlevp1
          DO je = i_startidx, i_endidx
            z_concorr_e(je,jk,jb) = &
              p_diag%vn_ie(je,jk,jb)*p_metrics%ddxn_z_half(je,jk,jb) + &
              z_vt_ie(je,jk)        *p_metrics%ddxt_z_half(je,jk,jb)
          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO

    rl_start = 2
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx)
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

      ! Interpolate contravariant correction to cell centers
      ! To save computing time, the contravariant correction calculated in the preceding
      ! time step is retained here for the predictor step; the impact is that it does
      ! not include the effect of the latest horizontal diffusion call
      IF (istep == 1 .AND. l_init) THEN

#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = nflatlev(p_patch%id)+1, nlevp1
#else
!CDIR UNROLL=6
        DO jk = nflatlev(p_patch%id)+1, nlevp1
          DO jc = i_startidx, i_endidx
#endif

            p_diag%w_concorr_c(jc,jk,jb) =  &
              p_int%e_bln_c_s(jc,1,jb)*z_concorr_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
              p_int%e_bln_c_s(jc,2,jb)*z_concorr_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
              p_int%e_bln_c_s(jc,3,jb)*z_concorr_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))

          ENDDO
        ENDDO

      ENDIF
    ENDDO
!$OMP END DO

    rl_start = 4
    rl_end = min_rledge_int - 2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx)
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

    ENDDO
!$OMP END DO

    rl_start = 3
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

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
           (z_vnw(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))         *p_int%geofac_div(jc,1,jb) + &
            z_vnw(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))         *p_int%geofac_div(jc,2,jb) + &
            z_vnw(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))         *p_int%geofac_div(jc,3,jb))

          z_w_con_c(jc,jk,jb) = p_prog%w(jc,jk,jb)
        ENDDO
      ENDDO

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

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx)
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

    rl_start = grf_bdywidth_e+1
    rl_end = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Sum up terms of horizontal wind advection:
      ! grad(Ekin_h) + vt*(f+relvort_e) + wcon_e*dv/dz
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
! unrolling with 6 would be even better, but only if nlev is a multiple of 6
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
!$OMP END DO

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
  !!
  SUBROUTINE solve_nh (p_nh, p_patch, p_int, bufr, nnow, nnew, l_init, linit_vertnest, &
                       l_bdy_nudge, dtime)

    TYPE(t_nh_state),  TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state), TARGET, INTENT(IN)    :: p_int
    TYPE(t_patch),     TARGET, INTENT(IN)    :: p_patch
    TYPE(t_buffer_memory),     INTENT(INOUT) :: bufr

    ! Initialization switch that has to be .TRUE. at the initial time step only (not for restart)
    LOGICAL,                   INTENT(INOUT) :: l_init
    ! Initialization switch set in dynamics_integration
    LOGICAL,                   INTENT(IN)    :: linit_vertnest(2)
    ! Switch to determine if boundary nudging is executed
    LOGICAL,                   INTENT(IN)    :: l_bdy_nudge
    ! Time levels
    INTEGER,                   INTENT(IN)    :: nnow, nnew
    ! Time step
    REAL(wp),                  INTENT(IN)    :: dtime

    ! Local variables
    INTEGER  :: jb, jk, jc, je
    INTEGER  :: nlev, nlevp1              !< number of full levels
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER  :: rl_start, rl_end, istep, ntl1, ntl2, nvar, nshift
    INTEGER  :: ic, ie

    REAL(wp) :: z_theta_v_fl_e  (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_theta_v_e     (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_rho_e         (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_gradh_exner   (nproma,p_patch%nlev  ,p_patch%nblks_e), &
                z_concorr_e     (nproma,p_patch%nlevp1,p_patch%nblks_e)

    REAL(wp), TARGET :: z_vn_avg (nproma,p_patch%nlev,p_patch%nblks_e)

    REAL(wp) :: z_mass_fl_div   (nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_theta_v_fl_div(nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_th_ddz_exner_c(nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_dexner_dz_c (2,nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_exner_ex_pr   (nproma,p_patch%nlev  ,p_patch%nblks_c), &
                z_exner_pr      (nproma,p_patch%nlev  ,p_patch%nblks_c)

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
                z_vt_ie         (nproma,p_patch%nlevp1),          &
                z_thermal_exp   (nproma,p_patch%nblks_c),         &
                z_hydro_corr    (nproma,p_patch%nblks_e)


    REAL(wp):: z_theta1, z_theta2, z_raylfac
    INTEGER :: nproma_gradp, nblks_gradp, npromz_gradp, nlen_gradp
    LOGICAL :: lcompute, lcleanup

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

    !-------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_solve_nh)

    IF (lvert_nest .AND. (p_patch%nshift > 0)) THEN  
      l_vert_nested = .TRUE.
    ELSE
      l_vert_nested = .FALSE.
    ENDIF
    IF (lvert_nest .AND. (p_patch%nshift_child > 0)) THEN  
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
    IF (itime_scheme == 4 .OR. itime_scheme == 6) THEN ! Velocity advection 2nd order in time
      ntl1 = nnow
      ntl2 = nnew
    ELSE                        ! Velocity advection 1st order in time
      ntl1 = 1
      ntl2 = 1
    ENDIF

    i_nchdom   = MAX(1,p_patch%n_childdom)

    DO istep = 1, 2

      IF (istep == 1) THEN ! predictor step
        IF (.NOT.((itime_scheme == 3 .OR. itime_scheme == 4) .AND. .NOT. l_init)) THEN
          CALL velocity_tendencies(p_nh%prog(nnow),p_patch,p_int,p_nh%metrics,&
                                   p_nh%diag,ntl1,istep,l_init)
        ENDIF
        nvar = nnow
      ELSE                 ! corrector step
        CALL velocity_tendencies(p_nh%prog(nnew),p_patch,p_int,p_nh%metrics,&
                                 p_nh%diag,ntl2,istep,l_init)
        nvar = nnew
      ENDIF

      l_init = .FALSE. ! should be .TRUE. only at initial predictor step

    ! Compute rho and theta at edges
    IF (istep == 1) THEN
      IF (iadv_rhotheta == 2) THEN

        lcompute =.TRUE.
        lcleanup =.FALSE.
        ! First call: compute backward trajectory with wind at time level nnow
        CALL upwind_hflux_miura(p_patch, p_nh%prog(nnow)%rho, p_nh%prog(nnow)%vn,         &
                                p_nh%prog(nnow)%vn, dtime, p_int, lcompute, lcleanup,     &
                                2, 0, 1, z_rho_e, opt_rlstart=7, opt_rlend=min_rledge_int,&
                                opt_lout_edge=.TRUE., opt_real_vt=p_nh%diag%vt )

        ! Second call: compute only reconstructed value for flux divergence
        lcompute =.FALSE.
        lcleanup =.TRUE.
        CALL upwind_hflux_miura(p_patch, p_nh%prog(nnow)%theta_v, p_nh%prog(nnow)%vn, &
                                p_nh%prog(nnow)%vn, dtime, p_int, lcompute, lcleanup, &
                                2, 0, 1, z_theta_v_e, opt_rlstart=7,                  &
                                opt_rlend=min_rledge_int-1, opt_lout_edge=.TRUE.,     &
                                opt_real_vt=p_nh%diag%vt )

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
    IF (istep == 1 .AND. igradp_method == 3) THEN
      
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

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk,jb,i_startidx,i_endidx)

    rl_start = 3
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    ! Computations at cell points; to be executed in predictor step only
    IF (istep == 1) THEN

!$OMP DO PRIVATE(jk,jc,z_exner_ic)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

!CDIR UNROLL=2
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            ! extrapolated perturbation Exner pressure (used for horizontal gradients only)
            z_exner_ex_pr(jc,jk,jb) = p_nh%diag%exner_fphy_incr(jc,jk,jb) -                &
              p_nh%metrics%exner_ref_mc(jc,jk,jb) +                                        &
              (1._wp + p_nh%metrics%exner_exfac(jc,jk,jb))*p_nh%prog(nnow)%exner(jc,jk,jb) &
                     - p_nh%metrics%exner_exfac(jc,jk,jb) *p_nh%diag%exner_old(jc,jk,jb)

            ! non-extrapolated perturbation Exner pressure
            z_exner_pr(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb) + &
              p_nh%diag%exner_fphy_incr(jc,jk,jb) - p_nh%metrics%exner_ref_mc(jc,jk,jb)

            ! Now save current time level in exner_old
            p_nh%diag%exner_old(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb)

            ! and reset the fast-physics increment to zero
            p_nh%diag%exner_fphy_incr(jc,jk,jb) = 0._wp
          ENDDO
        ENDDO

        IF (l_open_ubc .AND. .NOT. l_vert_nested) THEN
          ! Compute contribution of thermal expansion to vertical wind at model top
          ! Isobaric expansion is assumed
          z_thermal_exp(:,jb) = 0._wp
!CDIR UNROLL=4
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              z_thermal_exp(jc,jb) = z_thermal_exp(jc,jb) + cpd_o_rd                &
                * (p_nh%diag%ddt_exner(jc,jk,jb)+p_nh%diag%ddt_exner_phy(jc,jk,jb)) &
                /  p_nh%prog(nnow)%exner(jc,jk,jb)*p_nh%metrics%ddqz_z_full(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF

        ! Perturbation Exner pressure on bottom and top half level
        IF (nflatlev(p_patch%id) == 1) THEN
          DO jc = i_startidx, i_endidx
          z_exner_ic(jc,1) =                                          &
            p_nh%metrics%wgtfacq1_c(jc,1,jb)*z_exner_ex_pr(jc,1,jb) + &
            p_nh%metrics%wgtfacq1_c(jc,2,jb)*z_exner_ex_pr(jc,2,jb) + &
            p_nh%metrics%wgtfacq1_c(jc,3,jb)*z_exner_ex_pr(jc,3,jb)
          ENDDO
        ENDIF
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
      ENDDO
!$OMP END DO
    ENDIF ! istep = 1

    IF (igradp_method == 1) THEN
      rl_end   = min_rlcell_int
      i_endblk = p_patch%cells%end_blk(rl_end,i_nchdom)
    ENDIF


!$OMP DO PRIVATE(jk,jc,z_theta_v_pr_mc,z_theta_v_pr_ic)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      z_theta_v_pr_mc(i_startidx:i_endidx,1) =  0.5_wp *     &
        (p_nh%prog(nnow)%theta_v(i_startidx:i_endidx,1,jb) + &
        p_nh%prog(nvar)%theta_v(i_startidx:i_endidx,1,jb)) - &
        p_nh%metrics%theta_ref_mc(i_startidx:i_endidx,1,jb)

!CDIR UNROLL=8
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          ! density at interface levels for vertical flux divergence computation
          p_nh%diag%rho_ic(jc,jk,jb) = p_nh%metrics%rho_refcorr_ic(jc,jk,jb)+0.5_wp*( &
            p_nh%metrics%wgtfac_c(jc,jk,jb)*(p_nh%prog(nnow)%rho(jc,jk,jb) +          &
            p_nh%prog(nvar)%rho(jc,jk,jb))+(1._wp-p_nh%metrics%wgtfac_c(jc,jk,jb))*   &
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

      ! rho and theta at top level (fields are interpolated from parent domain in case of vertical nesting)
      IF (l_vert_nested) THEN
        DO jc = i_startidx, i_endidx
          p_nh%diag%theta_v_ic(jc,1,jb) = p_nh%diag%theta_v_ic(jc,2,jb) + &
            p_nh%diag%dtheta_v_ic_ubc(jc,jb)
          p_nh%diag%rho_ic(jc,1,jb)     = p_nh%diag%rho_ic(jc,2,jb) +    &
            p_nh%diag%drho_ic_ubc(jc,jb)
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

!CDIR UNROLL=3
        DO jk = nflat_gradp(p_patch%id), nlev
          DO jc = i_startidx, i_endidx
            ! Second vertical derivative of perturbation Exner pressure (hydrostatic approximation)
            z_dexner_dz_c(2,jc,jk,jb) = -0.5_wp / p_nh%prog(nnow)%theta_v(jc,jk,jb) * &
             ((z_theta_v_pr_ic(jc,jk) - z_theta_v_pr_ic(jc,jk+1)) *                   &
              p_nh%metrics%d_exner_dz_ref_mc(jc,jk,jb) + z_theta_v_pr_mc(jc,jk)*      &
              p_nh%metrics%d2_exner_dz2_ref_mc(jc,jk,jb))
          ENDDO
        ENDDO

        ! Store values at nest interface levels
        IF (linit_vertnest(1) .AND. l_child_vertnest) THEN
          DO jc = i_startidx, i_endidx
            p_nh%diag%drho_ic_int    (jc,jb) = p_nh%diag%rho_ic(jc,nshift,jb)     - &
                                               p_nh%diag%rho_ic(jc,nshift+1,jb)
            p_nh%diag%dtheta_v_ic_int(jc,jb) = p_nh%diag%theta_v_ic(jc,nshift,jb) - &
                                               p_nh%diag%theta_v_ic(jc,nshift+1,jb)
          ENDDO
        ENDIF
        IF (linit_vertnest(2) .AND. l_child_vertnest) THEN
          DO jc = i_startidx, i_endidx
            p_nh%diag%dw_int(jc,jb) = p_nh%prog(nnow)%w(jc,nshift,jb)    - &
                                      p_nh%prog(nnow)%w(jc,nshift+1,jb)
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
!$OMP DO PRIVATE(jk,je,z_theta1,z_theta2)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Store values at nest interface levels
        IF (linit_vertnest(1) .AND. l_child_vertnest) THEN
          DO je = i_startidx, i_endidx
            p_nh%diag%dvn_ie_int(je,jb) = p_nh%diag%vn_ie(je,nshift,jb) - &
                                          p_nh%diag%vn_ie(je,nshift+1,jb)
          ENDDO
        ENDIF

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nflatlev(p_patch%id)-1
#else
!CDIR UNROLL=_URD-_URD2
        DO jk = 1, nflatlev(p_patch%id)-1
          DO je = i_startidx, i_endidx
#endif
            ! horizontal gradient of Exner pressure where coordinate surfaces are flat
            z_gradh_exner(je,jk,jb) = p_patch%edges%inv_dual_edge_length(je,jb)* &
             (z_exner_ex_pr(icidx(je,jb,2),jk,icblk(je,jb,2)) -                  &
              z_exner_ex_pr(icidx(je,jb,1),jk,icblk(je,jb,1)) )
          ENDDO
        ENDDO

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

        IF (igradp_method >= 2) THEN
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
        ENDIF

        IF (igradp_method == 3) THEN
        ! compute hydrostatically approximated correction term that replaces downward extrapolation
          
          DO je = i_startidx, i_endidx

            z_theta1 = &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,1),ikidx(1,je,nlev,jb),icblk(je,jb,1)) + &
              p_nh%metrics%zdiff_gradp(1,je,nlev,jb)*                                      &
             (p_nh%diag%theta_v_ic(icidx(je,jb,1),ikidx(1,je,nlev,jb),  icblk(je,jb,1)) -   &
              p_nh%diag%theta_v_ic(icidx(je,jb,1),ikidx(1,je,nlev,jb)+1,icblk(je,jb,1))) *  &
              p_nh%metrics%inv_ddqz_z_full(icidx(je,jb,1),ikidx(1,je,nlev,jb),icblk(je,jb,1))

            z_theta2 = &
              p_nh%prog(nnow)%theta_v(icidx(je,jb,2),ikidx(2,je,nlev,jb),icblk(je,jb,2)) + &
              p_nh%metrics%zdiff_gradp(2,je,nlev,jb)*                                      &
             (p_nh%diag%theta_v_ic(icidx(je,jb,2),ikidx(2,je,nlev,jb),  icblk(je,jb,2)) -   &
              p_nh%diag%theta_v_ic(icidx(je,jb,2),ikidx(2,je,nlev,jb)+1,icblk(je,jb,2))) *  &
              p_nh%metrics%inv_ddqz_z_full(icidx(je,jb,2),ikidx(2,je,nlev,jb),icblk(je,jb,2))

            z_hydro_corr(je,jb) = grav_o_cpd*p_patch%edges%inv_dual_edge_length(je,jb)*    &
              (z_theta2-z_theta1)*4._wp/(z_theta1+z_theta2)**2

            ENDDO
        ENDIF

      ENDDO
!$OMP END DO
    ENDIF ! istep = 1

    IF (istep == 1 .AND. igradp_method == 3) THEN

!$OMP DO PRIVATE(jb,je,ie,nlen_gradp)
      DO jb = 1, nblks_gradp
        IF (jb == nblks_gradp) THEN
          nlen_gradp = npromz_gradp
        ELSE
          nlen_gradp = nproma_gradp
        ENDIF

!CDIR NODEP,VOVERTAKE,VOB
        DO je = 1, nlen_gradp
          ie = (jb-1)*nproma_gradp+je

          z_gradh_exner(ipeidx(ie),iplev(ie),ipeblk(ie))  =              &
            z_gradh_exner(ipeidx(ie),iplev(ie),ipeblk(ie)) +             &
            p_nh%metrics%pg_exdist(ie)*z_hydro_corr(ipeidx(ie),ipeblk(ie))

        ENDDO
      ENDDO
!$OMP END DO
    ENDIF

    ! Update horizontal velocity field
!$OMP DO PRIVATE(jk,je,ic)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF ((itime_scheme == 4 .OR. itime_scheme == 6) .AND. istep == 2) THEN
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)+ dtime     &
            & *(0.5_wp*(p_nh%diag%ddt_vn_adv(je,jk,jb,ntl1)                        &
            & +p_nh%diag%ddt_vn_adv(je,jk,jb,ntl2))+p_nh%diag%ddt_vn_phy(je,jk,jb) &
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
    ENDDO
!$OMP END DO

    IF (istep == 2 .AND. l_bdy_nudge) THEN ! apply boundary nudging if requested
!$OMP DO PRIVATE(jk,je,ic)
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


      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

! OpenMP parallelization is done over jk here because the number of blocks is
! too small for reasonable load balance (may be ifdef'd to go over blocks for other platforms)
!$OMP DO PRIVATE(jk,je)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb) + &
              dtime*p_nh%diag%grf_tend_vn(je,jk,jb)
          ENDDO
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF

    IF (itype_comm == 2) THEN
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
    ENDIF

!$OMP END PARALLEL

    IF (itype_comm == 1) THEN
      IF (istep == 1) THEN
        IF (idiv_method == 1) THEN
          CALL sync_patch_array_mult(SYNC_E,p_patch,2,p_nh%prog(nnew)%vn,z_rho_e)
        ELSE
          CALL sync_patch_array_mult(SYNC_E,p_patch,3,p_nh%prog(nnew)%vn,z_rho_e,z_theta_v_e)
      ENDIF
      ELSE
        CALL sync_patch_array(SYNC_E,p_patch,p_nh%prog(nnew)%vn)
      ENDIF
    ENDIF

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    rl_start = 2
    rl_end   = min_rledge_int - 2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (idiv_method == 1) THEN

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

      ELSE ! idiv_method = 2

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif
            ! Perform only RBF reconstruction of tangential wind component
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
    ENDDO
!$OMP END DO

    rl_start = 3
    rl_end = min_rledge_int - 2

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_vt_ie)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Interpolate vn and vt to interface levels
!CDIR UNROLL=6
      DO jk = nflatlev(p_patch%id), nlev
        DO je = i_startidx, i_endidx
          p_nh%diag%vn_ie(je,jk,jb) = &
            p_nh%metrics%wgtfac_e(je,jk,jb)*p_nh%prog(nnew)%vn(je,jk,jb) +        &
           (1._wp - p_nh%metrics%wgtfac_e(je,jk,jb))*p_nh%prog(nnew)%vn(je,jk-1,jb)
          z_vt_ie(je,jk) =                                                  &
            p_nh%metrics%wgtfac_e(je,jk,jb)*p_nh%diag%vt(je,jk,jb) +        &
           (1._wp - p_nh%metrics%wgtfac_e(je,jk,jb))*p_nh%diag%vt(je,jk-1,jb)
        ENDDO
      ENDDO

      ! Reduced computations required for flat levels
      IF (istep == 1) THEN
        DO jk = 2, nflatlev(p_patch%id)-1
          DO je = i_startidx, i_endidx
            p_nh%diag%vn_ie(je,jk,jb) = &
              p_nh%metrics%wgtfac_e(je,jk,jb)*p_nh%prog(nnew)%vn(je,jk,jb) + &
             (1._wp - p_nh%metrics%wgtfac_e(je,jk,jb))*p_nh%prog(nnew)%vn(je,jk-1,jb)
          ENDDO
        ENDDO
      ENDIF

      IF (istep == 1 .AND. .NOT. l_vert_nested) THEN
        DO je = i_startidx, i_endidx
          p_nh%diag%vn_ie(je,1,jb) =                                       &
            p_nh%metrics%wgtfacq1_e(je,1,jb)*p_nh%prog(nnew)%vn(je,1,jb) + &
            p_nh%metrics%wgtfacq1_e(je,2,jb)*p_nh%prog(nnew)%vn(je,2,jb) + &
            p_nh%metrics%wgtfacq1_e(je,3,jb)*p_nh%prog(nnew)%vn(je,3,jb)
        ENDDO
      ELSE IF (istep == 1 .AND. l_vert_nested) THEN
        DO je = i_startidx, i_endidx
          p_nh%diag%vn_ie(je,1,jb) = p_nh%diag%vn_ie(je,2,jb) + &
            p_nh%diag%dvn_ie_ubc(je,jb)
        ENDDO
      ENDIF

      ! Bottom level
      DO je = i_startidx, i_endidx
        p_nh%diag%vn_ie(je,nlevp1,jb) =                                    &
          p_nh%metrics%wgtfacq_e(je,1,jb)*p_nh%prog(nnew)%vn(je,nlev,jb)   + &
          p_nh%metrics%wgtfacq_e(je,2,jb)*p_nh%prog(nnew)%vn(je,nlev-1,jb) + &
          p_nh%metrics%wgtfacq_e(je,3,jb)*p_nh%prog(nnew)%vn(je,nlev-2,jb)
        z_vt_ie(je,nlevp1) =                                           &
          p_nh%metrics%wgtfacq_e(je,1,jb)*p_nh%diag%vt(je,nlev,jb) +   &
          p_nh%metrics%wgtfacq_e(je,2,jb)*p_nh%diag%vt(je,nlev-1,jb) + &
          p_nh%metrics%wgtfacq_e(je,3,jb)*p_nh%diag%vt(je,nlev-2,jb)
      ENDDO

      ! Compute contravariant correction for vertical velocity at interface levels
      DO jk = nflatlev(p_patch%id)+1, nlevp1
        DO je = i_startidx, i_endidx
          z_concorr_e(je,jk,jb) =                                            &
            p_nh%diag%vn_ie(je,jk,jb)*p_nh%metrics%ddxn_z_half(je,jk,jb) + &
            z_vt_ie(je,jk)*p_nh%metrics%ddxt_z_half(je,jk,jb)
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO

    rl_start = 3
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Interpolate contravariant correction to cell centers
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = nflatlev(p_patch%id)+1, nlevp1
#else
!CDIR UNROLL=6
      DO jk = nflatlev(p_patch%id)+1, nlevp1
        DO jc = i_startidx, i_endidx
#endif

          p_nh%diag%w_concorr_c(jc,jk,jb) =                                          &
            p_int%e_bln_c_s(jc,1,jb)*z_concorr_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) + &
            p_int%e_bln_c_s(jc,2,jb)*z_concorr_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) + &
            p_int%e_bln_c_s(jc,3,jb)*z_concorr_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO

    rl_start = 7
    IF (idiv_method == 1) THEN
      rl_end = min_rledge_int - 2
    ELSE
      rl_end = min_rledge_int - 3
    ENDIF

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je)
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

    ENDDO
!$OMP END DO
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
!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk,jb,i_startidx,i_endidx)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jk,jc,z_w_expl,z_contr_w_fl_l,z_rho_expl,z_exner_expl,z_a,z_b,z_c,&
!$OMP            z_g,z_q,z_alpha,z_beta,z_gamma,ic,z_raylfac)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      IF (istep == 2 .AND. (itime_scheme == 4 .OR. itime_scheme == 6)) THEN
!CDIR UNROLL=5
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx

            ! explicit part for w
            z_w_expl(jc,jk) = p_nh%prog(nnow)%w(jc,jk,jb) + dtime *                          &
              (0.5_wp*(p_nh%diag%ddt_w_adv(jc,jk,jb,ntl1)+p_nh%diag%ddt_w_adv(jc,jk,jb,ntl2))&
              -cpd*z_th_ddz_exner_c(jc,jk,jb))

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
          z_contr_w_fl_l(jc,1) = p_nh%diag%rho_ic(jc,1,jb)*p_nh%prog(nnow)%w(jc,1,jb)   &
            * p_nh%metrics%vwind_expl_wgt(jc,jb)
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
!CDIR ON_ADB(p_nh%prog(nnew)%w)
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
!CDIR ON_ADB(p_nh%prog(nnew)%w)
      DO jc = i_startidx, i_endidx
        p_nh%prog(nnew)%w(jc,2,jb)= p_nh%prog(nnew)%w(jc,2,jb)/z_b(jc,2)
      ENDDO

      DO jk = 3, nlev
!CDIR ON_ADB(p_nh%prog(nnew)%w)
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,jk,jb) = (p_nh%prog(nnew)%w(jc,jk,jb)  &
            -z_a(jc,jk)*p_nh%prog(nnew)%w(jc,jk-1,jb))*z_g(jc,jk)
        ENDDO
      ENDDO

      DO jk = nlev-1, 2, -1
!CDIR ON_ADB(p_nh%prog(nnew)%w)
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
      DO jk = 2, nrdmax(p_patch%id)
        z_raylfac = 1.0_wp/(1.0_wp+dtime*p_nh%metrics%rayleigh_w(jk))
!CDIR ON_ADB(p_nh%prog(nnew)%w)
        DO jc = i_startidx, i_endidx
          p_nh%prog(nnew)%w(jc,jk,jb) = z_raylfac*p_nh%prog(nnew)%w(jc,jk,jb) +    &
                                        (1._wp-z_raylfac)*p_nh%prog(nnew)%w(jc,1,jb)
        ENDDO
      ENDDO

      ! Results
      DO jk = 1, nlev
!CDIR ON_ADB(p_nh%prog(nnew)%w)
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

      IF (istep == 2) THEN
        ! store dynamical part of exner time increment in exner_dyn_incr
        ! the conversion into a temperature tendency is done in the NWP interface
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            p_nh%diag%exner_dyn_incr(jc,jk,jb) = p_nh%diag%exner_dyn_incr(jc,jk,jb) + &
              p_nh%prog(nnew)%exner(jc,jk,jb) - p_nh%diag%exner_old(jc,jk,jb) -   &
              dtime*p_nh%diag%ddt_exner_phy(jc,jk,jb)
          ENDDO
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

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

! OpenMP parallelization is done over jk here because the number of blocks is
! too small for reasonable load balance (may be ifdef'd to go over blocks for other platforms)
!$OMP DO PRIVATE(jk,jc)
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
        ENDDO
!OMP END DO
      ENDDO

    ELSE IF (istep == 1 .AND. (l_limited_area .OR. p_patch%id > 1)) THEN
      ! In the MPI-parallelized case, only rho and w are updated here,
      ! and theta_v is preliminarily stored on exner in order to save
      ! halo communications

      rl_start = 1
      rl_end   = grf_bdywidth_c

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,1)

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

! OpenMP parallelization is done over jk here because the number of blocks is
! too small for reasonable load balance (may be ifdef'd to go over blocks for other platforms)
!$OMP DO PRIVATE(jk,jc)
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
!OMP END DO
      ENDDO

    ENDIF

    IF (itype_comm == 2) THEN
      ! use OpenMP-parallelized communication using global memory for buffers
      IF (istep == 1) THEN ! Only w is updated in the predictor step
        CALL sync_patch_array_gm(SYNC_C,p_patch,1,bufr%send_c1,bufr%recv_c1,p_nh%prog(nnew)%w)
      ELSE IF (istep == 2) THEN
        ! Synchronize all prognostic variables
        CALL sync_patch_array_gm(SYNC_C,p_patch,3,bufr%send_c3,bufr%recv_c3,                &
                                 p_nh%prog(nnew)%rho,p_nh%prog(nnew)%exner,p_nh%prog(nnew)%w)
      ENDIF
    ENDIF

!$OMP END PARALLEL

    IF (itype_comm == 1) THEN
      IF (istep == 1) THEN ! Only w is updated in the predictor step
        CALL sync_patch_array(SYNC_C,p_patch,p_nh%prog(nnew)%w)
      ELSE IF (istep == 2) THEN
        ! Synchronize all prognostic variables
        CALL sync_patch_array_mult(SYNC_C,p_patch,3,p_nh%prog(nnew)%rho,  &
                                   p_nh%prog(nnew)%exner,p_nh%prog(nnew)%w)
      ENDIF
    ENDIF

    ENDDO ! istep-loop

    IF (ltimer) CALL timer_stop(timer_solve_nh)

    ! The remaining computations are needed for MPI-parallelized applications only
    IF (my_process_is_mpi_all_seq() ) RETURN

!$OMP PARALLEL PRIVATE(rl_start,rl_end,jb,i_startblk,i_endblk,i_startidx,i_endidx)
    IF (l_limited_area .OR. p_patch%id > 1) THEN

      ! Index list over halo points lying in the boundary interpolation zone
      ! Note: this list typically contains at most 10 grid points 
!$OMP DO PRIVATE(ic,jk,jc)
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
!$OMP END DO

      rl_start = 1
      rl_end   = grf_bdywidth_c

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,1)

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

!$OMP DO PRIVATE(jk,jc)
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
!$OMP END DO
      ENDDO
    ENDIF

    rl_start = min_rlcell_int - 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

!$OMP DO PRIVATE(jk,jc)
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
!$OMP END DO
    ENDDO
!$OMP END PARALLEL


  END SUBROUTINE solve_nh

END MODULE mo_solve_nonhydro
