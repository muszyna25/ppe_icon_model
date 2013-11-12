!>
!! mo_nh_diffusion
!!
!! Diffusion in the nonhydrostatic model
!!
!! @author Almut Gassmann, MPI-M
!!
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2009-08.25)
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

MODULE mo_nh_diffusion

  USE mo_kind,                ONLY: wp, vp
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: l_limited_area, lfeedback, grid_sphere_radius
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_cell
  USE mo_interpol_config,     ONLY: nudge_max_coeff
  USE mo_intp,                ONLY: verts2edges_scalar, edges2verts_scalar, &
                                    cells2verts_scalar, cells2edges_scalar, &
                                    edges2cells_scalar, verts2cells_scalar, &
                                    edges2cells_vector
  USE mo_nonhydrostatic_config, ONLY: l_zdiffu_t, damp_height, iadv_rcf, lhdiff_rcf
  USE mo_diffusion_config,    ONLY: diffusion_config
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ltimer, iforcing, lvert_nest
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_impl_constants    ,  ONLY: min_rledge, min_rlcell, min_rlvert, &
                                    min_rledge_int, min_rlcell_int, min_rlvert_int, inwp
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_e, grf_bdywidth_c
  USE mo_math_laplace,        ONLY: nabla4_vec
  USE mo_math_constants,      ONLY: dbl_eps, pi
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_gridref_config,      ONLY: denom_diffu_v
  USE mo_parallel_config,     ONLY: p_test_run, itype_comm
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, &
                                    sync_patch_array_mult
  USE mo_physical_constants,  ONLY: cvd_o_rd, cpd, rd, p0ref
  USE mo_timer,               ONLY: timer_nh_hdiffusion, timer_start, timer_stop
  USE mo_vertical_grid,       ONLY: nrdmax
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_exception,           ONLY: finish

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: diffusion

  CONTAINS

  !>
  !! diffusion
  !!
  !! Computes the horizontal diffusion of velocity and temperature
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD (2010-10-13), based on an earlier
  !! version initially developed by Almut Gassmann, MPI-M
  !!
  SUBROUTINE  diffusion(p_nh_prog,p_nh_diag,p_nh_metrics,p_patch,p_int,dtime,linit)

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch    !< single patch
    TYPE(t_int_state),INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_nh_prog), INTENT(inout)    :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_diag), INTENT(inout)    :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    REAL(wp), INTENT(in)            :: dtime      !< time step
    LOGICAL,  INTENT(in)            :: linit      !< initial call or runtime call

    ! local variables - vp means variable precision depending on the __MIXED_PRECISION cpp flag
    REAL(vp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: z_temp
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_nabla2_e
    REAL(vp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: z_nabla2_c
    REAL(vp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_nabla4_e

    REAL(wp):: diff_multfac_vn(p_patch%nlev), diff_multfac_w, diff_multfac_n2w(p_patch%nlev)
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc, je, ic, ishift
    INTEGER :: nlev              !< number of full levels

    ! start index levels and diffusion coefficient for boundary diffusion
    INTEGER :: start_bdydiff_e
    REAL(wp):: fac_bdydiff_v

    ! For Smagorinsky diffusion - vp means variable precision depending on the __MIXED_PRECISION cpp flag
    REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
    REAL(vp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: kh_smag_e
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: u_vert
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: v_vert
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: u_cell
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: v_cell
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt_norm, dvt_tang, smag_offset,   &
                nabv_tang, nabv_norm, rd_o_cvd, nudgezone_diff, bdy_diff, vn_cell1, vn_cell2
    REAL(wp), DIMENSION(p_patch%nlev) :: smag_limit, diff_multfac_smag, enh_smag_fac
    INTEGER  :: nblks_zdiffu, nproma_zdiffu, npromz_zdiffu, nlen_zdiffu

    ! Variables for provisional fix against runaway cooling in local topography depressions
    INTEGER  :: icount(p_patch%nblks_c), iclist(2*nproma,p_patch%nblks_c), iklist(2*nproma,p_patch%nblks_c)
    REAL(wp) :: tdlist(2*nproma,p_patch%nblks_c), tdiff, trefdiff, enh_diffu, thresh_tdiff

    INTEGER,  DIMENSION(:,:,:), POINTER :: icidx, icblk, ieidx, ieblk, ividx, ivblk, &
                                           iecidx, iecblk
    INTEGER,  DIMENSION(:,:),   POINTER :: icell, ilev, iblk, iedge, iedblk
    REAL(wp), DIMENSION(:,:),   POINTER :: vcoef, blcoef, geofac_n2s
    LOGICAL :: ltemp_diffu
    INTEGER :: diffu_type, discr_vn, discr_t
    INTEGER :: jg                 !< patch ID

    !--------------------------------------------------------------------------

    ! The diffusion is an intrinsic part of the NH solver, thus it is added to the timer
    IF (ltimer) CALL timer_start(timer_nh_hdiffusion)

    ! get patch ID
    jg = p_patch%id

    start_bdydiff_e = 5 ! refin_ctrl level at which boundary diffusion starts

    ! number of vertical levels
    nlev = p_patch%nlev

    ! Normalized diffusion coefficient for boundary diffusion
    IF (lhdiff_rcf) THEN
      fac_bdydiff_v = SQRT(REAL(iadv_rcf,wp))/denom_diffu_v
    ELSE
      fac_bdydiff_v = 1._wp/denom_diffu_v
    ENDIF

    ! scaling factor for enhanced diffusion in nudging zone (if present, i.e. for
    ! limited-area runs and one-way nesting)
    nudgezone_diff = 0.04_wp/(nudge_max_coeff + dbl_eps)

    ! scaling factor for enhanced near-boundary diffusion for 
    ! two-way nesting (used with Smagorinsky diffusion only; not needed otherwise)
    bdy_diff = 0.015_wp/(nudge_max_coeff + dbl_eps)

    ! threshold temperature deviation from neighboring grid points 
    ! that activates extra diffusion against runaway cooling
    thresh_tdiff = - 5._wp

    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    icidx => p_patch%cells%neighbor_idx
    icblk => p_patch%cells%neighbor_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    rd_o_cvd = 1._wp/cvd_o_rd

    i_nchdom   = MAX(1,p_patch%n_childdom)

    diffu_type  = diffusion_config(jg)%hdiff_order
    discr_vn    = diffusion_config(jg)%itype_vn_diffu
    discr_t     = diffusion_config(jg)%itype_t_diffu

    IF (linit) THEN ! enhanced diffusion at all levels for initial velocity filtering call
      diff_multfac_vn(:) = diffusion_config(jg)%k4/3._wp*diffusion_config(jg)%hdiff_efdt_ratio
      smag_offset        =  0.0_wp
      diffu_type = 5 ! always combine nabla4 background diffusion with Smagorinsky diffusion for initial filtering call
      smag_limit(:) = 0.125_wp-4._wp*diff_multfac_vn(:)
    ELSE IF (lhdiff_rcf) THEN ! combination with divergence damping inside the dynamical core
      IF (diffu_type == 4) THEN
        diff_multfac_vn(:) = MIN(1._wp/128._wp,diffusion_config(jg)%k4*REAL(iadv_rcf,wp)/ &
                                 3._wp*p_nh_metrics%enhfac_diffu(:))
      ELSE ! For Smagorinsky diffusion, the Smagorinsky coefficient rather than the background
           ! diffusion coefficient is enhanced near the model top (see below)
        diff_multfac_vn(:) = MIN(1._wp/128._wp,diffusion_config(jg)%k4*REAL(iadv_rcf,wp)/3._wp)
      ENDIF
      IF (diffu_type == 3) THEN
        smag_offset   = 0._wp
        smag_limit(:) = 0.125_wp
      ELSE
        smag_offset   = 0.25_wp*diffusion_config(jg)%k4*REAL(iadv_rcf,wp)
        smag_limit(:) = 0.125_wp-4._wp*diff_multfac_vn(:)
      ENDIF
    ELSE           ! enhanced diffusion near model top only
      IF (diffu_type == 4) THEN
        diff_multfac_vn(:) = diffusion_config(jg)%k4/3._wp*p_nh_metrics%enhfac_diffu(:)
      ELSE ! For Smagorinsky diffusion, the Smagorinsky coefficient rather than the background
           ! diffusion coefficient is enhanced near the model top (see below)
        diff_multfac_vn(:) = diffusion_config(jg)%k4/3._wp
      ENDIF
      smag_offset        = 0.25_wp*diffusion_config(jg)%k4
      smag_limit(:)      = 0.125_wp-4._wp*diff_multfac_vn(:)
      ! pure Smagorinsky diffusion does not work without divergence damping
      IF (diffusion_config(jg)%hdiff_order == 3) diffu_type = 5 
    ENDIF

    ! Multiplication factor for nabla4 diffusion on vertical wind speed
    diff_multfac_w = MIN(1._wp/48._wp,diffusion_config(jg)%k4w*REAL(iadv_rcf,wp))

    ! Factor for additional nabla2 diffusion in upper damping zone
    diff_multfac_n2w(:) = 0._wp
    IF (nrdmax(jg) > 1) THEN ! seems to be redundant, but the NEC issues invalid operations otherwise
      DO jk = 2, nrdmax(jg)
        diff_multfac_n2w(jk) = 1._wp/12._wp*((vct_a(jk)-vct_a(nrdmax(jg)+1))/(vct_a(2)-vct_a(nrdmax(jg)+1)))**4
      ENDDO
    ENDIF

    IF (diffu_type == 3 .OR. diffu_type == 5) THEN

      ! temperature diffusion is used only in combination with Smagorinsky diffusion
      ltemp_diffu = diffusion_config(jg)%lhdiff_temp

      ! enhanced factor for Smagorinsky diffusion above the stratopause in order to
      ! properly damp breaking gravity waves
      !
      ! linear increase starting at 25 km, reaching a value of 0.1 at 75 km
      enh_smag_fac(1:nlev) = MIN(0.1_wp,MAX(0._wp,(0.5_wp*(vct_a(1:nlev)+vct_a(2:nlev+1))-25000._wp)/500000._wp))
      ! ... combined with quadratic increase starting at 50 km, reaching a value of 1 at 90 km
      enh_smag_fac(1:nlev) = MIN(1._wp,MAX(enh_smag_fac(1:nlev),                                &
                             MAX(0._wp,(0.5_wp*(vct_a(1:nlev)+vct_a(2:nlev+1))-50000._wp)/40000._wp)**2) )

      ! Smagorinsky coefficient is also enhanced in the six model levels beneath a vertical nest interface
      IF ((lvert_nest) .AND. (p_patch%nshift > 0)) THEN
        enh_smag_fac(1) = MAX(0.333_wp, enh_smag_fac(1))
        enh_smag_fac(2) = MAX(0.25_wp, enh_smag_fac(2))
        enh_smag_fac(3) = MAX(0.20_wp, enh_smag_fac(3))
        enh_smag_fac(4) = MAX(0.16_wp, enh_smag_fac(4))
        enh_smag_fac(5) = MAX(0.12_wp, enh_smag_fac(5))
        enh_smag_fac(6) = MAX(0.08_wp, enh_smag_fac(6))
      ENDIF

      ! empirically determined scaling factor
      diff_multfac_smag(:) = MAX(diffusion_config(jg)%hdiff_smag_fac,enh_smag_fac(:))*dtime

      IF (lhdiff_rcf) diff_multfac_smag(:) = diff_multfac_smag(:)*REAL(iadv_rcf,wp)

    ELSE
      ltemp_diffu = .FALSE.
    ENDIF

    IF (diffu_type == 4) THEN

#ifdef __MIXED_PRECISION
      CALL finish ('', 'hdiff_order=4 not implemented for mixed precision mode; compile without -D__MIXED_PRECISION')
#else
      CALL nabla4_vec( p_nh_prog%vn, p_patch, p_int, z_nabla4_e,  &
                       opt_rlstart=7,opt_nabla2=z_nabla2_e )
#endif

    ELSE IF ((diffu_type == 3 .OR. diffu_type == 5) .AND. discr_vn == 1) THEN

      IF (p_test_run) THEN
        u_vert = 0._wp
        v_vert = 0._wp
      ENDIF

      !  RBF reconstruction of velocity at vertices
      CALL rbf_vec_interpol_vertex( p_nh_prog%vn, p_patch, p_int, &
                                    u_vert, v_vert, opt_rlend=min_rlvert_int )

      rl_start = start_bdydiff_e
      rl_end   = min_rledge_int - 2

      IF (itype_comm == 1 .OR. itype_comm == 3) THEN
        CALL sync_patch_array_mult(SYNC_V,p_patch,2,u_vert,v_vert)
      ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      i_startblk = p_patch%edges%start_blk(rl_start,1)
      i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,vn_vert1,vn_vert2,vn_vert3,vn_vert4,&
!$OMP             dvt_norm,dvt_tang),  ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Computation of wind field deformation

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif

            vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v2

            vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v2

            dvt_tang = p_patch%edges%system_orientation(je,jb)* (   &
                       u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v2 - &
                      (u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v2) )

            vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                       v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%primal_normal_vert(je,jb,3)%v2

            vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                       v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%primal_normal_vert(je,jb,4)%v2

            dvt_norm = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%dual_normal_vert(je,jb,4)%v1 + &
                       v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                       p_patch%edges%dual_normal_vert(je,jb,4)%v2 - &
                      (u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%dual_normal_vert(je,jb,3)%v1 + &
                       v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                       p_patch%edges%dual_normal_vert(je,jb,3)%v2)

            ! Smagorinsky diffusion coefficient
            kh_smag_e(je,jk,jb) = diff_multfac_smag(jk)*SQRT(             (  &
              (vn_vert4-vn_vert3)*p_patch%edges%inv_vert_vert_length(je,jb) -&
              dvt_tang*p_patch%edges%inv_primal_edge_length(je,jb) )**2 + (  &
              (vn_vert2-vn_vert1)*p_patch%edges%system_orientation(je,jb)*   &
              p_patch%edges%inv_primal_edge_length(je,jb) +         &
              dvt_norm*p_patch%edges%inv_vert_vert_length(je,jb))**2 )

            ! The factor of 4 comes from dividing by twice the "correct" length
            z_nabla2_e(je,jk,jb) = 4._wp * (                        &
              (vn_vert4 + vn_vert3 - 2._wp*p_nh_prog%vn(je,jk,jb))  &
              *p_patch%edges%inv_vert_vert_length(je,jb)**2 +       &
              (vn_vert2 + vn_vert1 - 2._wp*p_nh_prog%vn(je,jk,jb))  &
              *p_patch%edges%inv_primal_edge_length(je,jb)**2 )

          ENDDO
        ENDDO

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            ! Subtract part of the fourth-order background diffusion coefficient
            kh_smag_e(je,jk,jb) = MAX(0._wp,kh_smag_e(je,jk,jb) - smag_offset)
            ! Limit diffusion coefficient to the theoretical CFL stability threshold
            kh_smag_e(je,jk,jb) = MIN(kh_smag_e(je,jk,jb),smag_limit(jk))
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ELSE IF ((diffu_type == 3 .OR. diffu_type == 5) .AND. discr_vn >= 2) THEN

      !  RBF reconstruction of velocity at vertices and cells
      CALL rbf_vec_interpol_vertex( p_nh_prog%vn, p_patch, p_int, &
                                    u_vert, v_vert, opt_rlend=min_rlvert_int-1 )

      IF (discr_vn == 2) THEN
        CALL rbf_vec_interpol_cell( p_nh_prog%vn, p_patch, p_int, &
                                    u_cell, v_cell, opt_rlend=min_rlcell_int-1 )
      ELSE
        CALL edges2cells_vector( p_nh_prog%vn, p_nh_diag%vt, p_patch, p_int, &
                                 u_cell, v_cell, opt_rlend=min_rlcell_int-1 )
      ENDIF

      IF (p_test_run) THEN
        z_nabla2_e = 0._wp
      ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,rl_start,rl_end)

      rl_start = start_bdydiff_e
      rl_end   = min_rledge_int - 1

      i_startblk = p_patch%edges%start_blk(rl_start,1)
      i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,vn_vert1,vn_vert2,vn_cell1,vn_cell2,&
!$OMP             dvt_norm,dvt_tang), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ! Computation of wind field deformation

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif

            vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_vert(je,jb,1)%v2

            vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_vert(je,jb,2)%v2

            dvt_tang = p_patch%edges%system_orientation(je,jb)* (   &
                       u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_vert(je,jb,2)%v2 - &
                      (u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_vert(je,jb,1)%v2) )

            vn_cell1 = u_cell(iecidx(je,jb,1),jk,iecblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                       v_cell(iecidx(je,jb,1),jk,iecblk(je,jb,1)) * &
                       p_patch%edges%primal_normal_cell(je,jb,1)%v2

            vn_cell2 = u_cell(iecidx(je,jb,2),jk,iecblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                       v_cell(iecidx(je,jb,2),jk,iecblk(je,jb,2)) * &
                       p_patch%edges%primal_normal_cell(je,jb,2)%v2

            dvt_norm = u_cell(iecidx(je,jb,2),jk,iecblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_cell(je,jb,2)%v1 + &
                       v_cell(iecidx(je,jb,2),jk,iecblk(je,jb,2)) * &
                       p_patch%edges%dual_normal_cell(je,jb,2)%v2 - &
                      (u_cell(iecidx(je,jb,1),jk,iecblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_cell(je,jb,1)%v1 + &
                       v_cell(iecidx(je,jb,1),jk,iecblk(je,jb,1)) * &
                       p_patch%edges%dual_normal_cell(je,jb,1)%v2)

            ! Smagorinsky diffusion coefficient
            kh_smag_e(je,jk,jb) = diff_multfac_smag(jk)*SQRT(             (  &
              (vn_cell2-vn_cell1)*p_patch%edges%inv_dual_edge_length(je,jb) -&
              dvt_tang*p_patch%edges%inv_primal_edge_length(je,jb) )**2 + (  &
              (vn_vert2-vn_vert1)*p_patch%edges%system_orientation(je,jb)*   &
              p_patch%edges%inv_primal_edge_length(je,jb) +         &
              dvt_norm*p_patch%edges%inv_dual_edge_length(je,jb))**2 )

            ! The factor of 4 comes from dividing by twice the "correct" length
            z_nabla2_e(je,jk,jb) = 4._wp * (                        &
              (vn_cell2 + vn_cell1 - 2._wp*p_nh_prog%vn(je,jk,jb))  &
              *p_patch%edges%inv_dual_edge_length(je,jb)**2 +       &
              (vn_vert2 + vn_vert1 - 2._wp*p_nh_prog%vn(je,jk,jb))  &
              *p_patch%edges%inv_primal_edge_length(je,jb)**2 )

          ENDDO
        ENDDO

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            ! Subtract part of the fourth-order background diffusion coefficient
            kh_smag_e(je,jk,jb) = MAX(0._wp,kh_smag_e(je,jk,jb) - smag_offset)
            ! Limit diffusion coefficient to the theoretical CFL stability threshold
            kh_smag_e(je,jk,jb) = MIN(kh_smag_e(je,jk,jb),smag_limit(jk))
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF

    IF (diffu_type == 5) THEN ! Add fourth-order background diffusion

      IF (discr_vn > 1) CALL sync_patch_array(SYNC_E,p_patch,z_nabla2_e)

      ! Interpolate nabla2(v) to vertices in order to compute nabla2(nabla2(v))

      IF (p_test_run) THEN
        u_vert = 0._wp
        v_vert = 0._wp
      ENDIF

      CALL rbf_vec_interpol_vertex( z_nabla2_e, p_patch, p_int, u_vert, v_vert, &
                                    opt_rlstart=4, opt_rlend=min_rlvert_int )

      rl_start = grf_bdywidth_e+1
      rl_end   = min_rledge_int

      IF (itype_comm == 1 .OR. itype_comm == 3) THEN
        CALL sync_patch_array_mult(SYNC_V,p_patch,2,u_vert,v_vert)
      ENDIF

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      i_startblk = p_patch%edges%start_blk(rl_start,1)
      i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,nabv_tang,nabv_norm), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

         ! Compute nabla4(v)
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
            DO jk = 1, nlev
#else
!CDIR UNROLL=5
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
#endif

              nabv_tang = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                          p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                          v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
                          p_patch%edges%primal_normal_vert(je,jb,1)%v2 + &
                          u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                          p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                          v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
                          p_patch%edges%primal_normal_vert(je,jb,2)%v2

              nabv_norm = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                          p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
                          v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
                          p_patch%edges%primal_normal_vert(je,jb,3)%v2 + &
                          u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                          p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
                          v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
                          p_patch%edges%primal_normal_vert(je,jb,4)%v2

              ! The factor of 4 comes from dividing by twice the "correct" length
              z_nabla4_e(je,jk,jb) = 4._wp * (                        &
                (nabv_norm - 2._wp*z_nabla2_e(je,jk,jb))              &
                *p_patch%edges%inv_vert_vert_length(je,jb)**2 +       &
                (nabv_tang - 2._wp*z_nabla2_e(je,jk,jb))              &
                *p_patch%edges%inv_primal_edge_length(je,jb)**2 )

            ENDDO
          ENDDO

        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF

    ! For nested domains, the diffusion contribution is already
    ! included in the time tendency interpolated from the parent
    ! domain. Thus, the boundary zone is skipped here.

    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

    IF (diffu_type == 5) THEN ! Smagorinsky diffusion combined with fourth-order background diffusion
      IF ( jg == 1 .AND. l_limited_area .OR. jg > 1 .AND. .NOT. lfeedback(jg)) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb)  +                    &
                p_patch%edges%area_edge(je,jb) *                                    &
                (MAX(nudgezone_diff*p_int%nudgecoeff_e(je,jb),kh_smag_e(je,jk,jb))* &
                z_nabla2_e(je,jk,jb) - diff_multfac_vn(jk) * z_nabla4_e(je,jk,jb) * &
                p_patch%edges%area_edge(je,jb))
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
      ELSE IF (jg > 1) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb)  +               &
                p_patch%edges%area_edge(je,jb) * (kh_smag_e(je,jk,jb)*         &
                z_nabla2_e(je,jk,jb) - z_nabla4_e(je,jk,jb) *                  &
                MAX(diff_multfac_vn(jk),bdy_diff*p_int%nudgecoeff_e(je,jb)) *  &
                p_patch%edges%area_edge(je,jb))
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
      ELSE
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb)  +                   &
                p_patch%edges%area_edge(je,jb) * (kh_smag_e(je,jk,jb)*             &
                z_nabla2_e(je,jk,jb) - diff_multfac_vn(jk) * z_nabla4_e(je,jk,jb) *&
                p_patch%edges%area_edge(je,jb))
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
      ENDIF
    ELSE IF (diffu_type == 3) THEN ! Only Smagorinsky diffusion
      IF ( jg == 1 .AND. l_limited_area .OR. jg > 1 .AND. .NOT. lfeedback(jg)) THEN
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb)  +                &
                p_patch%edges%area_edge(je,jb) * z_nabla2_e(je,jk,jb) *         &
                MAX(nudgezone_diff*p_int%nudgecoeff_e(je,jb),kh_smag_e(je,jk,jb))
               
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
      ELSE
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          DO jk = 1, nlev
!DIR$ IVDEP
            DO je = i_startidx, i_endidx
              p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb)  +                        &
                p_patch%edges%area_edge(je,jb) * kh_smag_e(je,jk,jb)* z_nabla2_e(je,jk,jb)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
      ENDIF
    ELSE
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlev
!DIR$ IVDEP
          DO je = i_startidx, i_endidx
            p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb)  -    &
              diff_multfac_vn(jk) * z_nabla4_e(je,jk,jb) *        &
              p_patch%edges%area_edge(je,jb)*p_patch%edges%area_edge(je,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF

    IF (l_limited_area .OR. jg > 1) THEN

      ! Lateral boundary diffusion for vn
      i_startblk = p_patch%edges%start_blk(start_bdydiff_e,1)
      i_endblk   = p_patch%edges%end_blk(grf_bdywidth_e,1)

!$OMP DO PRIVATE(jk,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, start_bdydiff_e, grf_bdywidth_e)

        DO jk = 1, nlev
          p_nh_prog%vn(i_startidx:i_endidx,jk,jb) =   &
            p_nh_prog%vn(i_startidx:i_endidx,jk,jb) + &
            z_nabla2_e(i_startidx:i_endidx,jk,jb) * &
            p_patch%edges%area_edge(i_startidx:i_endidx,jb)*fac_bdydiff_v
        ENDDO
      ENDDO
!$OMP END DO

    ENDIF ! vn boundary diffusion

    IF (lhdiff_rcf .AND. diffusion_config(jg)%lhdiff_w) THEN ! add diffusion on vertical wind speed
                     ! remark: the surface level (nlevp1) is excluded because w is diagnostic there

      rl_start = grf_bdywidth_c
      rl_end   = min_rlcell_int-1

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
#endif
            z_nabla2_c(jc,jk,jb) =  &
              p_nh_prog%w(jc,jk,jb)                        *p_int%geofac_n2s(jc,1,jb) + &
              p_nh_prog%w(icidx(jc,jb,1),jk,icblk(jc,jb,1))*p_int%geofac_n2s(jc,2,jb) + &
              p_nh_prog%w(icidx(jc,jb,2),jk,icblk(jc,jb,2))*p_int%geofac_n2s(jc,3,jb) + &
              p_nh_prog%w(icidx(jc,jb,3),jk,icblk(jc,jb,3))*p_int%geofac_n2s(jc,4,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
!DIR$ IVDEP
          DO jk = 1, nlev
#else
!CDIR UNROLL=5
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
#endif
            p_nh_prog%w(jc,jk,jb) = p_nh_prog%w(jc,jk,jb) - diff_multfac_w * p_patch%cells%area(jc,jb)**2 * &
             (z_nabla2_c(jc,jk,jb)                        *p_int%geofac_n2s(jc,1,jb) +                      &
              z_nabla2_c(icidx(jc,jb,1),jk,icblk(jc,jb,1))*p_int%geofac_n2s(jc,2,jb) +                      &
              z_nabla2_c(icidx(jc,jb,2),jk,icblk(jc,jb,2))*p_int%geofac_n2s(jc,3,jb) +                      &
              z_nabla2_c(icidx(jc,jb,3),jk,icblk(jc,jb,3))*p_int%geofac_n2s(jc,4,jb))
          ENDDO
        ENDDO

        ! Add nabla2 diffusion in upper damping layer (if present)
        DO jk = 2, nrdmax(jg)
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_nh_prog%w(jc,jk,jb) = p_nh_prog%w(jc,jk,jb) +                         &
              diff_multfac_n2w(jk) * p_patch%cells%area(jc,jb) * z_nabla2_c(jc,jk,jb)
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO
    ENDIF ! w diffusion

!$OMP END PARALLEL

    IF (itype_comm == 1 .OR. itype_comm == 3) THEN
      CALL sync_patch_array(SYNC_E, p_patch, p_nh_prog%vn)
    ENDIF

    IF (ltemp_diffu) THEN ! Smagorinsky temperature diffusion
      IF (l_zdiffu_t) THEN
        icell      => p_nh_metrics%zd_indlist
        iblk       => p_nh_metrics%zd_blklist
        ilev       => p_nh_metrics%zd_vertidx
        iedge      => p_nh_metrics%zd_edgeidx
        iedblk     => p_nh_metrics%zd_edgeblk
        vcoef      => p_nh_metrics%zd_intcoef
        blcoef     => p_nh_metrics%zd_e2cell
        geofac_n2s => p_nh_metrics%zd_geofac

        nproma_zdiffu = MIN(nproma,256)
        nblks_zdiffu = INT(p_nh_metrics%zd_listdim/nproma_zdiffu)
        npromz_zdiffu = MOD(p_nh_metrics%zd_listdim,nproma_zdiffu)
        IF (npromz_zdiffu > 0) THEN
          nblks_zdiffu = nblks_zdiffu + 1
        ELSE
          npromz_zdiffu = nproma_zdiffu
        ENDIF
      ENDIF

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)

      ! Enhance Smagorinsky diffusion coefficient in the presence of excessive grid-point cold pools
      ! This is restricted to the two lowest model levels
      !
      rl_start = grf_bdywidth_c
      rl_end   = min_rlcell_int-1

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx,ic,tdiff,trefdiff), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        ic = 0

        DO jk = nlev-1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            ! Perturbation potential temperature difference between local point and average of the three neighbors
            tdiff = p_nh_prog%theta_v(jc,jk,jb) -                          &
              (p_nh_prog%theta_v(icidx(jc,jb,1),jk,icblk(jc,jb,1)) +       &
               p_nh_prog%theta_v(icidx(jc,jb,2),jk,icblk(jc,jb,2)) +       &
               p_nh_prog%theta_v(icidx(jc,jb,3),jk,icblk(jc,jb,3)) ) / 3._wp 
            trefdiff = p_nh_metrics%theta_ref_mc(jc,jk,jb) -                       &
              (p_nh_metrics%theta_ref_mc(icidx(jc,jb,1),jk,icblk(jc,jb,1)) +       &
               p_nh_metrics%theta_ref_mc(icidx(jc,jb,2),jk,icblk(jc,jb,2)) +       &
               p_nh_metrics%theta_ref_mc(icidx(jc,jb,3),jk,icblk(jc,jb,3)) ) / 3._wp

            IF (tdiff-trefdiff < thresh_tdiff .AND. trefdiff < 0._wp) THEN
              ic = ic+1
              iclist(ic,jb) = jc
              iklist(ic,jb) = jk
              tdlist(ic,jb) = thresh_tdiff - tdiff
            ENDIF

          ENDDO
        ENDDO

        icount(jb) = ic

      ENDDO
!$OMP END DO

      ! Enhance Smagorinsky coefficients at the three edges of the cells included in the list
      ! Attention: this operation is neither vectorizable nor OpenMP-parallelizable (race conditions!)
!$OMP MASTER
      DO jb = i_startblk,i_endblk

        IF (icount(jb) > 0) THEN
          DO ic = 1, icount(jb)
            jc = iclist(ic,jb)
            jk = iklist(ic,jb)
            enh_diffu = tdlist(ic,jb)*5.e-4_wp
            kh_smag_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) = MAX(enh_diffu,kh_smag_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)))
            kh_smag_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) = MAX(enh_diffu,kh_smag_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)))
            kh_smag_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) = MAX(enh_diffu,kh_smag_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)))
          ENDDO
        ENDIF

      ENDDO
!$OMP END MASTER
!$OMP BARRIER

      IF (discr_t == 1) THEN  ! use discretization K*nabla(theta)

        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int

        i_startblk = p_patch%cells%start_blk(rl_start,1)
        i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! interpolated diffusion coefficient times nabla2(theta)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1, nlev
#else
!CDIR UNROLL=5
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
#endif
              z_temp(jc,jk,jb) =  &
               (kh_smag_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%e_bln_c_s(jc,1,jb)          + &
                kh_smag_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%e_bln_c_s(jc,2,jb)          + &
                kh_smag_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%e_bln_c_s(jc,3,jb))         * &
               (p_nh_prog%theta_v(jc,jk,jb)                        *p_int%geofac_n2s(jc,1,jb) + &
                p_nh_prog%theta_v(icidx(jc,jb,1),jk,icblk(jc,jb,1))*p_int%geofac_n2s(jc,2,jb) + &
                p_nh_prog%theta_v(icidx(jc,jb,2),jk,icblk(jc,jb,2))*p_int%geofac_n2s(jc,3,jb) + &
                p_nh_prog%theta_v(icidx(jc,jb,3),jk,icblk(jc,jb,3))*p_int%geofac_n2s(jc,4,jb))
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO

      ELSE IF (discr_t == 2) THEN ! use conservative discretization div(k*grad(theta))

        rl_start = grf_bdywidth_e
        rl_end   = min_rledge_int - 1

        i_startblk = p_patch%edges%start_blk(rl_start,1)
        i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! compute kh_smag_e * grad(theta) (stored in z_nabla2_e for memory efficiency)
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
!DIR$ IVDEP
            DO jk = 1, nlev
#else
!CDIR UNROLL=6
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
#endif
              z_nabla2_e(je,jk,jb) = kh_smag_e(je,jk,jb) *              &
                p_patch%edges%inv_dual_edge_length(je,jb)*              &
               (p_nh_prog%theta_v(iecidx(je,jb,2),jk,iecblk(je,jb,2)) - &
                p_nh_prog%theta_v(iecidx(je,jb,1),jk,iecblk(je,jb,1)))
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO


        rl_start = grf_bdywidth_c+1
        rl_end   = min_rlcell_int

        i_startblk = p_patch%cells%start_blk(rl_start,1)
        i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          ! now compute the divergence of the quantity above
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1, nlev
#else
!CDIR UNROLL=6
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
#endif
              z_temp(jc,jk,jb) =                                                         &
                z_nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
                z_nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
                z_nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO

      ENDIF


      IF (l_zdiffu_t) THEN ! Compute temperature diffusion truly horizontally over steep slopes
                           ! A conservative discretization is not possible here
!$OMP DO PRIVATE(jb,jc,ic,nlen_zdiffu,ishift) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, nblks_zdiffu
          IF (jb == nblks_zdiffu) THEN
            nlen_zdiffu = npromz_zdiffu
          ELSE
            nlen_zdiffu = nproma_zdiffu
          ENDIF
          ishift = (jb-1)*nproma_zdiffu
!CDIR NODEP,VOVERTAKE,VOB
          DO jc = 1, nlen_zdiffu
            ic = ishift+jc
            z_temp(icell(1,ic),ilev(1,ic),iblk(1,ic)) =                                          &
              z_temp(icell(1,ic),ilev(1,ic),iblk(1,ic)) + p_nh_metrics%zd_diffcoef(ic)*          &
!              MAX(p_nh_metrics%zd_diffcoef(ic),        &
!              kh_smag_e(iedge(1,ic),ilev(1,ic),iedblk(1,ic))* blcoef(1,ic)  +                    &
!              kh_smag_e(iedge(2,ic),ilev(1,ic),iedblk(2,ic))* blcoef(2,ic)  +                    &
!              kh_smag_e(iedge(3,ic),ilev(1,ic),iedblk(3,ic))* blcoef(3,ic) ) *                   &
             (geofac_n2s(1,ic)*p_nh_prog%theta_v(icell(1,ic),ilev(1,ic),iblk(1,ic)) +            &
              geofac_n2s(2,ic)*(vcoef(1,ic)*p_nh_prog%theta_v(icell(2,ic),ilev(2,ic),iblk(2,ic))+&
              (1._wp-vcoef(1,ic))* p_nh_prog%theta_v(icell(2,ic),ilev(2,ic)+1,iblk(2,ic)))  +    &
              geofac_n2s(3,ic)*(vcoef(2,ic)*p_nh_prog%theta_v(icell(3,ic),ilev(3,ic),iblk(3,ic))+&
              (1._wp-vcoef(2,ic))*p_nh_prog%theta_v(icell(3,ic),ilev(3,ic)+1,iblk(3,ic)))  +     &
              geofac_n2s(4,ic)*(vcoef(3,ic)*p_nh_prog%theta_v(icell(4,ic),ilev(4,ic),iblk(4,ic))+&
              (1._wp-vcoef(3,ic))* p_nh_prog%theta_v(icell(4,ic),ilev(4,ic)+1,iblk(4,ic)))  )
          ENDDO
        ENDDO
!$OMP END DO
      ENDIF

!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            p_nh_prog%theta_v(jc,jk,jb) = p_nh_prog%theta_v(jc,jk,jb) + &
              p_patch%cells%area(jc,jb)*z_temp(jc,jk,jb)
             
            p_nh_prog%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
              p_nh_prog%rho(jc,jk,jb)*p_nh_prog%theta_v(jc,jk,jb)))
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! This could be further optimized, but applications without physics are quite rare; 
      ! For LES physics, the sync is always necessary here
      IF ( .NOT. lhdiff_rcf .OR. linit .OR. iforcing /= inwp .OR. atm_phy_nwp_config(jg)%is_les_phy) THEN
        CALL sync_patch_array_mult(SYNC_C,p_patch,2,p_nh_prog%theta_v,p_nh_prog%exner)
      ENDIF

    ENDIF ! temperature diffusion

    IF ( .NOT. lhdiff_rcf .OR. linit .OR. iforcing /= inwp .OR. atm_phy_nwp_config(jg)%is_les_phy) THEN
      IF (diffusion_config(jg)%lhdiff_w) CALL sync_patch_array(SYNC_C,p_patch,p_nh_prog%w)
    ENDIF

    IF (ltimer) CALL timer_stop(timer_nh_hdiffusion)

  END SUBROUTINE diffusion


END MODULE mo_nh_diffusion
