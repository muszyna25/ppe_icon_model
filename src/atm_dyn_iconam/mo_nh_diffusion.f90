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

  USE mo_kind,                ONLY: wp
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_buffer_memory
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: nroot, l_limited_area, lfeedback, grid_sphere_radius
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_cell
  USE mo_interpol_config,     ONLY: nudge_max_coeff
  USE mo_intp,                ONLY: verts2edges_scalar, edges2verts_scalar, &
                                    cells2verts_scalar, cells2edges_scalar, &
                                    edges2cells_scalar, verts2cells_scalar, &
                                    edges2cells_vector
  USE mo_nonhydrostatic_config, ONLY: l_zdiffu_t, damp_height, k2_updamp_coeff, &
                                      iadv_rcf, lhdiff_rcf
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
                                    sync_patch_array_mult, sync_patch_array_gm
  USE mo_physical_constants,  ONLY: cvd_o_rd, cpd, rd, p0ref
  USE mo_timer,               ONLY: timer_nh_hdiffusion, timer_start, timer_stop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: diffusion_tria, diffusion_hex

  CONTAINS

  !>
  !! diffusion_tria
  !!
  !! Computes the horizontal diffusion of velocity and temperature for the triangular grid
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD (2010-10-13), based on an earlier
  !! version initially developed by Almut Gassmann, MPI-M
  !!
  SUBROUTINE  diffusion_tria(p_nh_prog,p_nh_diag,p_nh_metrics,p_patch,p_int,bufr,dtime,linit)

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch    !< single patch
    TYPE(t_int_state),INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_nh_prog), INTENT(inout)    :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_diag), INTENT(inout)    :: p_nh_diag  !< single nh diagnostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    TYPE(t_buffer_memory), INTENT(INOUT) :: bufr
    REAL(wp), INTENT(in)            :: dtime      !< time step
    LOGICAL,  INTENT(in)            :: linit      !< initial call or runtime call

    ! local variables
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c)   :: z_temp
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e)   :: z_nabla2_e
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e)   :: z_nabla4_e

    REAL(wp):: diff_multfac_vn(p_patch%nlev)
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, jc, je, ic

    ! start index levels and diffusion coefficient for boundary diffusion
    INTEGER :: start_bdydiff_e
    REAL(wp):: fac_bdydiff_v

    ! For Smagorinsky diffusion
    REAL(wp), PARAMETER :: rd_o_p0ref = rd / p0ref
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: kh_smag_e
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: u_vert
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: v_vert
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: u_cell
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: v_cell
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt_norm, dvt_tang, smag_offset,   &
                nabv_tang, nabv_norm, rd_o_cvd, nudgezone_diff, bdy_diff, vn_cell1, vn_cell2
    REAL(wp), DIMENSION(p_patch%nlev) :: smag_limit, diff_multfac_smag, enh_smag_fac
    INTEGER  :: nblks_zdiffu, nproma_zdiffu, npromz_zdiffu, nlen_zdiffu
    INTEGER  :: nlev              !< number of full levels

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
    fac_bdydiff_v = 1._wp/denom_diffu_v

    ! scaling factor for enhanced diffusion in nudging zone (if present, i.e. for
    ! limited-area runs and one-way nesting)
    nudgezone_diff = 0.04_wp/(nudge_max_coeff + dbl_eps)

    ! scaling factor for enhanced near-boundary diffusion for 
    ! two-way nesting (used with Smagorinsky diffusion only; not needed otherwise)
    bdy_diff = 0.01_wp/(nudge_max_coeff + dbl_eps)

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

    IF (diffu_type == 3 .OR. diffu_type == 5) THEN

      ! temperature diffusion is used only in combination with Smagorinsky diffusion
      ltemp_diffu = diffusion_config(jg)%lhdiff_temp

      ! enhanced factor for Smagorinsky diffusion above the stratopause in order to
      ! properly damp breaking gravity waves
      enh_smag_fac(1:nlev) = MIN(1._wp,MAX(0._wp,(0.5_wp*(vct_a(1:nlev)+vct_a(2:nlev+1))-50000._wp)/40000._wp)**2)

      ! Smagorinsky coefficient is also enhanced in the three model levels beneath a vertical nest interface
      IF ((lvert_nest) .AND. (p_patch%nshift > 0)) THEN
        enh_smag_fac(1) = MAX(0.25_wp, enh_smag_fac(1))
        enh_smag_fac(2) = MAX(0.175_wp,enh_smag_fac(2))
        enh_smag_fac(3) = MAX(0.10_wp, enh_smag_fac(3))
      ENDIF

      ! empirically determined scaling factor (default of 0.15 for hdiff_smag_fac is somewhat
      ! larger than suggested in the literature)
      diff_multfac_smag(:) = MAX(diffusion_config(jg)%hdiff_smag_fac,enh_smag_fac(:))*dtime

      IF (lhdiff_rcf) diff_multfac_smag(:) = diff_multfac_smag(:)*REAL(iadv_rcf,wp)

    ELSE
      ltemp_diffu = .FALSE.
    ENDIF

    IF (diffu_type == 4) THEN

      CALL nabla4_vec( p_nh_prog%vn, p_patch, p_int, z_nabla4_e,  &
                       opt_rlstart=7,opt_nabla2=z_nabla2_e )

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

      IF (itype_comm == 2) THEN
        ! use OpenMP-parallelized communication using global memory for buffers
        CALL sync_patch_array_gm(SYNC_V,p_patch,2,bufr%send_v2,bufr%recv_v2,u_vert,v_vert)
      ENDIF

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

      rl_start = start_bdydiff_e
      rl_end   = min_rledge_int - 1

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

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

      IF (itype_comm == 2) THEN
        ! use OpenMP-parallelized communication using global memory for buffers
        CALL sync_patch_array_gm(SYNC_V,p_patch,2,bufr%send_v2,bufr%recv_v2,u_vert,v_vert)
      ENDIF

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

    IF (itype_comm == 2) THEN
      ! use OpenMP-parallelized communication using global memory for buffers
      CALL sync_patch_array_gm(SYNC_E,p_patch,1,bufr%send_e1,bufr%recv_e1,p_nh_prog%vn)
    ENDIF

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
!$OMP DO PRIVATE(jb,jc,ic,nlen_zdiffu) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, nblks_zdiffu
          IF (jb == nblks_zdiffu) THEN
            nlen_zdiffu = npromz_zdiffu
          ELSE
            nlen_zdiffu = nproma_zdiffu
          ENDIF
!CDIR NODEP,VOVERTAKE,VOB
          DO jc = 1, nlen_zdiffu
            ic = (jb-1)*nproma_zdiffu+jc
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

      IF (lhdiff_rcf) THEN
!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              p_nh_prog%theta_v(jc,jk,jb) = p_nh_prog%theta_v(jc,jk,jb) + &
                p_patch%cells%area(jc,jb)*z_temp(jc,jk,jb)

              p_nh_prog%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd_o_p0ref* &
                p_nh_prog%rho(jc,jk,jb)*p_nh_prog%theta_v(jc,jk,jb)))
            ENDDO
          ENDDO

          IF (iforcing /= inwp .OR. linit) THEN
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
              p_nh_prog%rhotheta_v(jc,jk,jb) = p_nh_prog%theta_v(jc,jk,jb) * &
                p_nh_prog%rho(jc,jk,jb)
              ENDDO
            ENDDO
          ENDIF

        ENDDO
!$OMP END DO NOWAIT
      ELSE ! diffusion is called at every sound-wave time step
!$OMP DO PRIVATE(jk,jc,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, rl_start, rl_end)

          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              p_nh_diag%ddt_exner(jc,jk,jb) = p_patch%cells%area(jc,jb)* rd_o_cvd / dtime * &
                p_nh_prog%exner(jc,jk,jb)/p_nh_prog%theta_v(jc,jk,jb)*z_temp(jc,jk,jb)
            ENDDO
          ENDDO

        ENDDO
!$OMP END DO NOWAIT
      ENDIF
!$OMP END PARALLEL

      ! This could be further optimized, but applications without physics are quite rare
      IF (lhdiff_rcf .AND. (linit .OR. iforcing /= inwp) )             &
        CALL sync_patch_array_mult(SYNC_C,p_patch,3,p_nh_prog%theta_v, &
          p_nh_prog%rhotheta_v,p_nh_prog%exner)

    ENDIF ! temperature diffusion

    IF (ltimer) CALL timer_stop(timer_nh_hdiffusion)

  END SUBROUTINE diffusion_tria


  !>
  !! diffusion_hex
  !!
  !! Computes the horizontal diffusion of velocity and temperature for the hexagonal grid
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-08.25)
  !! Modification by Guenther Zaengl, DWD (2010-10-13):
  !! Separation of diffusion routines for triangular and hexagonal grid
  !!
  SUBROUTINE  diffusion_hex(p_nh_prog,p_nh_metrics,p_patch,p_int,dtime)

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch    !< single patch
    TYPE(t_int_state),INTENT(in),TARGET :: p_int      !< single interpolation state
    TYPE(t_nh_prog), INTENT(inout)    :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_metrics),INTENT(in),TARGET :: p_nh_metrics !< single nh metric state
    REAL(wp), INTENT(in)            :: dtime      !< time step

    ! local variables
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e)   :: z_nabla4_e

    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: kh_smag_e
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: kh_smag_v
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: kh_smag_c

    REAL(wp):: diff_multfac_vn
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: jk, jb, jc, je, jv, id
    INTEGER :: nlev              !< number of full levels

    LOGICAL :: lsmag_diffu

    ! For Smagorinski diffusion on hexagonal model
    REAL(wp) :: z_rho_v(nproma,p_patch%nlev,p_patch%nblks_v), &
    &           z_rho_e(nproma,p_patch%nlev,p_patch%nblks_e), &
    &           z_shear_def_1 (nproma,p_patch%nlev,p_patch%nblks_e), &
    &           z_shear_def_2 (nproma,p_patch%nlev,p_patch%nblks_e), &
    &           z_strain_def_1(nproma,p_patch%nlev,p_patch%nblks_e), &
    &           z_strain_def_2(nproma,p_patch%nlev,p_patch%nblks_e), &
    &           z_fric_heat_c1(nproma,p_patch%nlev,p_patch%nblks_c), &
    &           z_fric_heat_v (nproma,p_patch%nlev,p_patch%nblks_v), &
    &           z_fric_heat_c (nproma,p_patch%nlev,p_patch%nblks_c), &
    &           z_turb_flx_c1 (nproma,p_patch%nlev,p_patch%nblks_e), &
    &           z_turb_flx_c2 (nproma,p_patch%nlev,p_patch%nblks_e), &
    &           z_turb_flx_v1 (nproma,p_patch%nlev,p_patch%nblks_e), &
    &           z_turb_flx_v2 (nproma,p_patch%nlev,p_patch%nblks_e), &
    &           z_old_rth(nproma)
    REAL(wp) :: zhelp, z_mean_area_edge
    INTEGER, POINTER :: ih1i(:,:,:), ih2i(:,:,:), ih1b(:,:,:), ih2b(:,:,:)
    INTEGER, POINTER :: it1i(:,:,:), it2i(:,:,:), it1b(:,:,:), it2b(:,:,:)
    INTEGER, POINTER :: ici (:,:,:), ivi (:,:,:), icb (:,:,:), ivb (:,:,:)
    INTEGER, POINTER :: icei(:,:,:), ivei(:,:,:), iceb(:,:,:), iveb(:,:,:)
    INTEGER :: jg      !< patch ID
    REAL(wp) :: sphere_radius_squared

    !-----------------------------------------------------------------------
    sphere_radius_squared = grid_sphere_radius * grid_sphere_radius
    !--------------------------------------------------------------------------

    ! get patch ID
    jg = p_patch%id

    IF( diffusion_config(jg)%hdiff_order == 3) THEN
      lsmag_diffu = .TRUE.
    ELSE
      lsmag_diffu = .FALSE.
    ENDIF

    i_nchdom   = MAX(1,p_patch%n_childdom)
    id         = p_patch%id

    ! number of vertical levels
    nlev = p_patch%nlev

    diff_multfac_vn = diffusion_config(id)%k4*3._wp

    IF (.NOT. lsmag_diffu) THEN
      CALL nabla4_vec( p_nh_prog%vn, p_patch, p_int, z_nabla4_e, opt_rlstart=7 )

      i_startblk = p_patch%edges%start_blk(grf_bdywidth_e+1,1)
      i_endblk   = p_patch%edges%end_blk(min_rledge,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, grf_bdywidth_e+1)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb)  -    &
              diff_multfac_vn * z_nabla4_e(je,jk,jb) * &
              p_patch%edges%area_edge(je,jb)*p_patch%edges%area_edge(je,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      CALL sync_patch_array(SYNC_E, p_patch, p_nh_prog%vn)

    ELSE
    ! Smagorinski diffusion for hexagonal model
    !------------------------------------------

      ! a) mean area edge
      z_mean_area_edge=8.0_wp*pi*sphere_radius_squared/REAL(p_patch%n_patch_edges_g,wp)

      ! b) compute density at vertices and edges
      CALL cells2verts_scalar(p_nh_prog%rho, p_patch, p_int%cells_aw_verts, z_rho_v)
      CALL cells2edges_scalar(p_nh_prog%rho, p_patch, p_int%c_lin_e, z_rho_e)

      ! c) abbreviate indices
      ici => p_patch%edges%cell_idx
      icb => p_patch%edges%cell_blk
      ivi => p_patch%edges%vertex_idx
      ivb => p_patch%edges%vertex_blk
      it1i => p_int%dir_gradt_i1
      it2i => p_int%dir_gradt_i2
      it1b => p_int%dir_gradt_b1
      it2b => p_int%dir_gradt_b2
      ih1i => p_int%dir_gradh_i1
      ih2i => p_int%dir_gradh_i2
      ih1b => p_int%dir_gradh_b1
      ih2b => p_int%dir_gradh_b2
      icei => p_patch%cells%edge_idx
      iceb => p_patch%cells%edge_blk
      ivei => p_patch%verts%edge_idx
      iveb => p_patch%verts%edge_blk

!$OMP PARALLEL

      IF (p_test_run) THEN
!$OMP WORKSHARE
        z_turb_flx_v1(:,:,:) = 0.0_wp
        z_turb_flx_v2(:,:,:) = 0.0_wp
        kh_smag_e(:,:,:)     = 0.0_wp
!$OMP END WORKSHARE
      ENDIF

      i_startblk = p_patch%edges%start_blk(3,1)
      i_endblk   = p_patch%edges%end_blk(min_rledge,1)
!$OMP DO PRIVATE(jb,jk,je,i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                  i_startidx, i_endidx, 3, min_rledge)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1,nlev
#else
        DO jk = 1,nlev 
          DO je = i_startidx, i_endidx
#endif
            ! d) Shear deformation at vertices
            z_shear_def_1(je,jk,jb) = &
            &(p_int%shear_def_v1(1,je,jb)*p_nh_prog%vn(it1i(1,je,jb),jk,it1b(1,je,jb)) &
            &+p_int%shear_def_v1(2,je,jb)*p_nh_prog%vn(it1i(2,je,jb),jk,it1b(2,je,jb)) &
            &+p_int%shear_def_v1(3,je,jb)*p_nh_prog%vn(it1i(3,je,jb),jk,it1b(3,je,jb)) &
            &+p_int%shear_def_v1(4,je,jb)*p_nh_prog%vn(it1i(4,je,jb),jk,it1b(4,je,jb)) &
            &+p_int%shear_def_v1(5,je,jb)*p_nh_prog%vn(it1i(5,je,jb),jk,it1b(5,je,jb)) &
            &+p_int%shear_def_v1(6,je,jb)*p_nh_prog%vn(it1i(6,je,jb),jk,it1b(6,je,jb)) &
            &+p_int%shear_def_v1(7,je,jb)*p_nh_prog%vn(it1i(7,je,jb),jk,it1b(7,je,jb)) &
            &+p_int%shear_def_v1(8,je,jb)*p_nh_prog%vn(it1i(8,je,jb),jk,it1b(8,je,jb)) &
            &+p_int%shear_def_v1(9,je,jb)*p_nh_prog%vn(it1i(9,je,jb),jk,it1b(9,je,jb)) )
            z_shear_def_2(je,jk,jb) = &
            &(p_int%shear_def_v2(1,je,jb)*p_nh_prog%vn(it2i(1,je,jb),jk,it2b(1,je,jb)) &
            &+p_int%shear_def_v2(2,je,jb)*p_nh_prog%vn(it2i(2,je,jb),jk,it2b(2,je,jb)) &
            &+p_int%shear_def_v2(3,je,jb)*p_nh_prog%vn(it2i(3,je,jb),jk,it2b(3,je,jb)) &
            &+p_int%shear_def_v2(4,je,jb)*p_nh_prog%vn(it2i(4,je,jb),jk,it2b(4,je,jb)) &
            &+p_int%shear_def_v2(5,je,jb)*p_nh_prog%vn(it2i(5,je,jb),jk,it2b(5,je,jb)) &
            &+p_int%shear_def_v2(6,je,jb)*p_nh_prog%vn(it2i(6,je,jb),jk,it2b(6,je,jb)) &
            &+p_int%shear_def_v2(7,je,jb)*p_nh_prog%vn(it2i(7,je,jb),jk,it2b(7,je,jb)) &
            &+p_int%shear_def_v2(8,je,jb)*p_nh_prog%vn(it2i(8,je,jb),jk,it2b(8,je,jb)) &
            &+p_int%shear_def_v2(9,je,jb)*p_nh_prog%vn(it2i(9,je,jb),jk,it2b(9,je,jb)) )
            ! e) Strain deformation at centers
            z_strain_def_1(je,jk,jb) = &
            &(p_int%strain_def_c1(1,je,jb)*p_nh_prog%vn(ih1i(1,je,jb),jk,ih1b(1,je,jb)) &
            &+p_int%strain_def_c1(2,je,jb)*p_nh_prog%vn(ih1i(2,je,jb),jk,ih1b(2,je,jb)) &
            &+p_int%strain_def_c1(3,je,jb)*p_nh_prog%vn(ih1i(3,je,jb),jk,ih1b(3,je,jb)) &
            &+p_int%strain_def_c1(4,je,jb)*p_nh_prog%vn(ih1i(4,je,jb),jk,ih1b(4,je,jb)) &
            &+p_int%strain_def_c1(5,je,jb)*p_nh_prog%vn(ih1i(5,je,jb),jk,ih1b(5,je,jb)) &
            &+p_int%strain_def_c1(6,je,jb)*p_nh_prog%vn(ih1i(6,je,jb),jk,ih1b(6,je,jb)) )
            z_strain_def_2(je,jk,jb) = &
            &(p_int%strain_def_c2(1,je,jb)*p_nh_prog%vn(ih2i(1,je,jb),jk,ih2b(1,je,jb)) &
            &+p_int%strain_def_c2(2,je,jb)*p_nh_prog%vn(ih2i(2,je,jb),jk,ih2b(2,je,jb)) &
            &+p_int%strain_def_c2(3,je,jb)*p_nh_prog%vn(ih2i(3,je,jb),jk,ih2b(3,je,jb)) &
            &+p_int%strain_def_c2(4,je,jb)*p_nh_prog%vn(ih2i(4,je,jb),jk,ih2b(4,je,jb)) &
            &+p_int%strain_def_c2(5,je,jb)*p_nh_prog%vn(ih2i(5,je,jb),jk,ih2b(5,je,jb)) &
            &+p_int%strain_def_c2(6,je,jb)*p_nh_prog%vn(ih2i(6,je,jb),jk,ih2b(6,je,jb)) )
            ! f) Diffusion coefficient at edge
            !kh_smag_e(je,jk,jb) = diffusion_config(jg)%hdiff_smag_fac &
            !  &                 * z_mean_area_edge
            !                      (Gibt grosse Koeffizienten am Pentagon)
            ! Erst quadrieren und dann mitteln ist besser.
            kh_smag_e(je,jk,jb) = diffusion_config(jg)%hdiff_smag_fac               &
            & *p_patch%edges%area_edge(je,jb)                                       &
            & *SQRT(0.5_wp*(z_strain_def_1(je,jk,jb)**2+z_strain_def_2(je,jk,jb)**2)&
            & +(p_patch%edges%edge_vert_length(je,jb,1)*(z_shear_def_1(je,jk,jb)**2) &
            &  +p_patch%edges%edge_vert_length(je,jb,2)*(z_shear_def_2(je,jk,jb)**2))&
            &  /p_patch%edges%primal_edge_length(je,jb) ) + k2_updamp_coeff*0.5_wp*(1.0_wp-&
            & COS(pi*MAX(0.0_wp,(p_nh_metrics%z_mc_e(je,jk,jb)-damp_height(1)))&
            &                  /(p_nh_metrics%z_mc_e(je, 1,jb)-damp_height(1))))
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      CALL sync_patch_array(SYNC_E, p_patch, kh_smag_e)

      ! g) average the diffusion coefficient to the centers and vertices
      ! to centers
      CALL edges2cells_scalar(kh_smag_e, p_patch, p_int%e_aw_c, kh_smag_c)
      ! to vertices
      CALL edges2verts_scalar(kh_smag_e, p_patch, p_int%e_aw_v, kh_smag_v)
      CALL verts2edges_scalar(kh_smag_v, p_patch, p_int%tria_aw_rhom, kh_smag_e)
      CALL sync_patch_array(SYNC_E, p_patch, kh_smag_e)
      CALL edges2verts_scalar(kh_smag_e, p_patch, p_int%e_1o3_v, kh_smag_v)

!$OMP PARALLEL
      i_startblk = p_patch%edges%start_blk(3,1)
      i_endblk   = p_patch%edges%end_blk(min_rledge,1)
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                  i_startidx, i_endidx, 3, min_rledge)
#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1,nlev
#else
        DO jk = 1,nlev 
          DO je = i_startidx, i_endidx
#endif
            ! h) turbulent fluxes
            z_turb_flx_c1(je,jk,jb) = &
            &    -p_nh_prog%rho(ici(je,jb,1),jk,icb(je,jb,1))  &
            &    *kh_smag_c(ici(je,jb,1),jk,icb(je,jb,1))      &
            &    *z_strain_def_1(je,jk,jb)
            z_turb_flx_c2(je,jk,jb) = &
            &    -p_nh_prog%rho(ici(je,jb,2),jk,icb(je,jb,2))  &
            &    *kh_smag_c(ici(je,jb,2),jk,icb(je,jb,2))      &
            &    *z_strain_def_2(je,jk,jb)
            z_turb_flx_v1(je,jk,jb) = &
            &    -z_rho_v(ivi(je,jb,1),jk,ivb(je,jb,1))        &
            &    *kh_smag_v(ivi(je,jb,1),jk,ivb(je,jb,1))      &
            &    *z_shear_def_1(je,jk,jb)
            z_turb_flx_v2(je,jk,jb) = &
            &    -z_rho_v(ivi(je,jb,2),jk,ivb(je,jb,2))        &
            &    *kh_smag_v(ivi(je,jb,2),jk,ivb(je,jb,2))      &
            &    *z_shear_def_2(je,jk,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF (diffusion_config(jg)%lhdiff_temp) THEN
        CALL sync_patch_array(SYNC_E, p_patch, z_turb_flx_v1)
        CALL sync_patch_array(SYNC_E, p_patch, z_turb_flx_v2)

        i_startblk = p_patch%cells%start_blk(2,1)
        i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,je,zhelp,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                  i_startidx, i_endidx, 2, min_rlcell)
          z_fric_heat_c(:,:,jb)=0.0_wp
          DO je = 1, 6
            DO jk = 1, nlev
              DO jc = i_startidx,i_endidx

                IF (je > p_patch%cells%num_edges(jc,jb)) CYCLE

                zhelp = p_patch%edges%system_orientation(icei(jc,jb,je),iceb(jc,jb,je)) &
                &      *p_patch%cells%edge_orientation(jc,jb,je)

                ! i) frictional heating at centers
                z_fric_heat_c(jc,jk,jb) = z_fric_heat_c(jc,jk,jb)                        &
                &-0.5_wp*((zhelp+1.0_wp)*z_turb_flx_c1(icei(jc,jb,je),jk,iceb(jc,jb,je)) &
                &        -(zhelp-1.0_wp)*z_turb_flx_c2(icei(jc,jb,je),jk,iceb(jc,jb,je)))&
                &*p_nh_prog%vn(icei(jc,jb,je),jk,iceb(jc,jb,je))*p_int%geofac_div(jc,je,jb)

              ENDDO
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
        i_startblk = p_patch%verts%start_blk(2,1)
        i_endblk   = p_patch%verts%end_blk(min_rlvert,1)
!$OMP DO PRIVATE(jb,jk,jv,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk
          CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
          &                  i_startidx, i_endidx, 2, min_rlvert)
#ifdef __LOOP_EXCHANGE
          DO jv = i_startidx, i_endidx
            DO jk = 1,nlev
#else
          DO jk = 1,nlev 
            DO jv = i_startidx, i_endidx
#endif
              ! j) frictional heating at vertices
              z_fric_heat_v(jv,jk,jb) = 0.5_wp* (&
              & ((p_patch%verts%edge_orientation(jv,jb,1)+1.0_wp) &
              &  *z_turb_flx_v1(ivei(jv,jb,1),jk,iveb(jv,jb,1))    &
              & -(p_patch%verts%edge_orientation(jv,jb,1)-1.0_wp) &
              &  *z_turb_flx_v2(ivei(jv,jb,1),jk,iveb(jv,jb,1)))   &
              &*p_nh_prog%vn(ivei(jv,jb,1),jk,iveb(jv,jb,1))*p_int%geofac_rot(jv,1,jb) &
              &+((p_patch%verts%edge_orientation(jv,jb,2)+1.0_wp) &
              &  *z_turb_flx_v1(ivei(jv,jb,2),jk,iveb(jv,jb,2))    &
              & -(p_patch%verts%edge_orientation(jv,jb,2)-1.0_wp) &
              &  *z_turb_flx_v2(ivei(jv,jb,2),jk,iveb(jv,jb,2)))   &
              &*p_nh_prog%vn(ivei(jv,jb,2),jk,iveb(jv,jb,2))*p_int%geofac_rot(jv,2,jb) &
              &+((p_patch%verts%edge_orientation(jv,jb,3)+1.0_wp) &
              &  *z_turb_flx_v1(ivei(jv,jb,3),jk,iveb(jv,jb,3))    &
              & -(p_patch%verts%edge_orientation(jv,jb,3)-1.0_wp) &
              &  *z_turb_flx_v2(ivei(jv,jb,3),jk,iveb(jv,jb,3)))   &
              &*p_nh_prog%vn(ivei(jv,jb,3),jk,iveb(jv,jb,3))*p_int%geofac_rot(jv,3,jb))

            ENDDO
          ENDDO
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ENDIF

      i_startblk = p_patch%edges%start_blk(3,1)
      i_endblk   = p_patch%edges%end_blk(min_rledge,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                  i_startidx, i_endidx, 3, min_rledge)
        DO jk = 1, nlev
          DO je = i_startidx,i_endidx
            ! k) Tendencies to the velocity components
            p_nh_prog%vn(je,jk,jb) = p_nh_prog%vn(je,jk,jb) - dtime  &
            & *((z_turb_flx_c2(je,jk,jb)- z_turb_flx_c1(je,jk,jb)) &
            &   *p_patch%edges%inv_dual_edge_length(je,jb)    &
            &   *p_patch%edges%system_orientation(je,jb)      &
            &  +(z_turb_flx_v2(je,jk,jb)- z_turb_flx_v1(je,jk,jb)) &
            &   *p_patch%edges%inv_primal_edge_length(je,jb)  &
            &  )/z_rho_e(je,jk,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      CALL sync_patch_array(SYNC_E, p_patch, p_nh_prog%vn )

      IF(diffusion_config(jg)%lhdiff_temp) THEN
        ! l) frictional heating at cells averaged from vertices
        CALL verts2cells_scalar(z_fric_heat_v, p_patch, p_int%verts_aw_cells, z_fric_heat_c1)
        i_startblk = p_patch%cells%start_blk(2,1)
        i_endblk   = p_patch%cells%end_blk(min_rlcell,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,z_old_rth,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                  i_startidx, i_endidx, 2, min_rlcell)
          DO jk = 1, nlev
            ! k) update new rho theta and Exner pressure
            z_old_rth(i_startidx:i_endidx) = p_nh_prog%rhotheta_v(i_startidx:i_endidx,jk,jb)
            p_nh_prog%rhotheta_v(i_startidx:i_endidx,jk,jb) = &
            & p_nh_prog%rhotheta_v(i_startidx:i_endidx,jk,jb) &
            & + dtime/cpd/p_nh_prog%exner(i_startidx:i_endidx,jk,jb) &
            & *(z_fric_heat_c(i_startidx:i_endidx,jk,jb)  &
            &  +z_fric_heat_c1(i_startidx:i_endidx,jk,jb))
            p_nh_prog%exner(i_startidx:i_endidx,jk,jb) = &
            & p_nh_prog%exner(i_startidx:i_endidx,jk,jb) &
            & *(1.0_wp+(p_nh_prog%rhotheta_v(i_startidx:i_endidx,jk,jb)&
            & /z_old_rth(i_startidx:i_endidx)-1.0_wp)/cvd_o_rd)
            p_nh_prog%theta_v(i_startidx:i_endidx,jk,jb)= &
            & p_nh_prog%rhotheta_v(i_startidx:i_endidx,jk,jb) &
            & /p_nh_prog%rho(i_startidx:i_endidx,jk,jb)
          ENDDO
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
        CALL sync_patch_array_mult(SYNC_C,p_patch,3,&
                                   p_nh_prog%rhotheta_v,&
                                   p_nh_prog%exner,&
                                   p_nh_prog%theta_v)
      ENDIF
    ENDIF

  END SUBROUTINE diffusion_hex

END MODULE mo_nh_diffusion

