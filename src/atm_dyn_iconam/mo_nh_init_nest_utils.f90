!>
!! This module contains the driver routines needed for initializing nested
!! domains that are started sometime during the integration
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2012-06-06)
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_init_nest_utils

  USE mo_kind,                  ONLY: wp
  USE mo_model_domain,          ONLY: t_patch, p_patch, p_patch_local_parent
  USE mo_nonhydro_types,        ONLY: t_nh_metrics, t_nh_prog
  USE mo_nonhydro_state,        ONLY: p_nh_state
  USE mo_parallel_config,       ONLY: nproma, p_test_run
  USE mo_run_config,            ONLY: ltransport, msg_level, ntracer
  USE mo_dynamics_config,       ONLY: nnow, nnow_rcf
  USE mo_physical_constants,    ONLY: rd, cvd_o_rd, p0ref
  USE mo_impl_constants,        ONLY: min_rlcell, min_rlcell_int, min_rledge_int, &
                                      min_rlvert, min_rlvert_int, MAX_CHAR_LENGTH
  USE mo_grf_nudgintp,          ONLY: interpol_scal_nudging, interpol_vec_nudging
  USE mo_grf_bdyintp,           ONLY: interpol_scal_grf, interpol2_vec_grf
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_mpi,                   ONLY: my_process_is_mpi_parallel
  USE mo_communication,         ONLY: exchange_data, exchange_data_mult
  USE mo_sync,                  ONLY: sync_patch_array, sync_patch_array_mult, &
                                      SYNC_C, SYNC_E, SYNC_V
  USE mo_intp_data_strc,        ONLY: t_int_state, p_int_state, p_int_state_local_parent
  USE mo_grf_intp_data_strc,    ONLY: t_gridref_single_state, t_gridref_state, &
                                      p_grf_state, p_grf_state_local_parent
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e, grf_fbk_start_c, &
                                      grf_nudgintp_start_c, grf_nudgintp_start_e
  USE mo_nwp_lnd_types,         ONLY: t_lnd_prog, t_lnd_diag
  USE mo_lnd_nwp_config,        ONLY: nsfc_subs, lmulti_snow, nlev_snow, nlev_soil
  USE mo_nwp_lnd_state,         ONLY: p_lnd_state
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config
  USE mo_interpol_config,       ONLY: nudge_zone_width
  USE mo_intp,                  ONLY: cells2verts_scalar
  USE mo_ext_data_state,        ONLY: ext_data

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  REAL(wp), PARAMETER :: rd_o_cvd = 1._wp / cvd_o_rd

  PUBLIC :: initialize_nest, topography_blending, topography_feedback

  CONTAINS
  !-------------
  !>
  !! SUBROUTINE initialize_nest
  !!
  !! Driver routine for initializing a nested domain during runtime from the parent domain
  !!
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2012-06-06)
  !!
  !! Note: the interpolation of the land variables does not yet include the
  !! land-water mask and adaptations to make optimal use of tiles
  !! The aggregate_landvars routine must be called before this routine
  !!
  !!
  SUBROUTINE initialize_nest(jg, jgc)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'initialize_nest'


    INTEGER, INTENT(IN) :: jg   ! parent (source) domain ID
    INTEGER, INTENT(IN) :: jgc  ! child  (target) domain ID


    ! local pointers
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog
    TYPE(t_nh_prog),    POINTER     :: p_child_prog
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog_rcf
    TYPE(t_nh_prog),    POINTER     :: p_child_prog_rcf
    TYPE(t_nh_metrics), POINTER     :: p_parent_metrics
    TYPE(t_nh_metrics), POINTER     :: p_child_metrics
    TYPE(t_lnd_prog),   POINTER     :: p_parent_lprog
    TYPE(t_lnd_prog),   POINTER     :: p_child_lprog
    TYPE(t_lnd_diag),   POINTER     :: p_parent_ldiag
    TYPE(t_lnd_diag),   POINTER     :: p_child_ldiag
    TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grfc => NULL()
    TYPE(t_int_state), POINTER      :: p_int => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()

    ! local variables

    ! Indices
    INTEGER :: jb, jc, jk, jk1, jt, je, i_nchdom, i_chidx, i_startblk, i_endblk, &
               i_startidx, i_endidx
    INTEGER :: nlev_c, nlev_p
    INTEGER :: nshift      ! difference between upper boundary of parent or feedback-parent 
                           ! domain and upper boundary of child domain (in terms 
                           ! of vertical levels) 
    INTEGER :: num_lndvars, num_snowvars

    ! Local arrays for variables living on the local parent grid in the MPI case. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: vn_lp, w_lp, thv_pr_lp, rho_pr_lp, &
                                                 lndvars_lp, snowvars_lp
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: tracer_lp

    ! Local arrays on the parent or child grid. These would not have to be allocatable,
    ! but the computational overhead does not matter for an initialization routine
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: thv_pr_par, rho_pr_par, lndvars_par, &
                                                 snowvars_par, lndvars_chi, snowvars_chi
    REAL(wp), ALLOCATABLE :: tsfc_ref_p(:,:), tsfc_ref_c(:,:) ! Reference temperature at lowest level

    ! Soil layer thicknesses (needed for limiting the soil water content)
    REAL(wp) :: dzsoil(9)=(/0.01_wp,0.01_wp,0.02_wp,0.06_wp,0.18_wp,0.54_wp,&
                            1.62_wp,4.86_wp,14.58_wp/)

    LOGICAL :: l_parallel, l_limit(ntracer)
    !-----------------------------------------------------------------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') 'Nest initialization, domain ',jg,' =>',jgc
      CALL message(TRIM(routine),message_text)
    ENDIF

    IF (.NOT. my_process_is_mpi_parallel()) THEN
      l_parallel = .FALSE.
    ELSE
      l_parallel = .TRUE.
    ENDIF

    p_parent_prog     => p_nh_state(jg)%prog(nnow(jg))
    p_child_prog      => p_nh_state(jgc)%prog(nnow(jgc))
    p_parent_prog_rcf => p_nh_state(jg)%prog(nnow_rcf(jg))
    p_child_prog_rcf  => p_nh_state(jgc)%prog(nnow_rcf(jgc))
    p_parent_metrics  => p_nh_state(jg)%metrics
    p_child_metrics   => p_nh_state(jgc)%metrics
    p_parent_lprog    => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
    p_child_lprog     => p_lnd_state(jgc)%prog_lnd(nnow_rcf(jgc))
    p_parent_ldiag    => p_lnd_state(jg)%diag_lnd
    p_child_ldiag     => p_lnd_state(jgc)%diag_lnd
    p_grfc            => p_grf_state(jgc)
    p_pc              => p_patch(jgc)

    IF (l_parallel) THEN
      p_grf => p_grf_state_local_parent(jgc)
      p_int => p_int_state_local_parent(jgc)
      p_pp  => p_patch_local_parent(jgc)
    ELSE
      p_grf => p_grf_state(jg)
      p_int => p_int_state(jg)
      p_pp  => p_patch(jg)
    ENDIF


    i_nchdom = MAX(1,p_patch(jg)%n_childdom)
    i_chidx  = p_pc%parent_child_index

    ! number of full levels of child domain
    nlev_c   = p_pc%nlev

    ! number of full levels of parent domain
    nlev_p   = p_pp%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift

    ! number of land and snow variables to be interpolated
    num_lndvars = 3*(nlev_soil+1)+1+ & ! multi-layer soil variables
                  5+2                  ! single-layer prognostic variables + t_g and freshsnow
    num_snowvars = 5*nlev_snow+1       ! snow fields
    
    ALLOCATE(thv_pr_par  (nproma, nlev_p,      p_patch(jg)%nblks_c), &
             rho_pr_par  (nproma, nlev_p,      p_patch(jg)%nblks_c), &
             lndvars_par (nproma, num_lndvars, p_patch(jg)%nblks_c), &
             snowvars_par(nproma, num_snowvars,p_patch(jg)%nblks_c), &
             lndvars_chi (nproma, num_lndvars, p_patch(jgc)%nblks_c),&
             snowvars_chi(nproma, num_snowvars,p_patch(jgc)%nblks_c),&
             tsfc_ref_p  (nproma,              p_patch(jg)%nblks_c), &
             tsfc_ref_c  (nproma,              p_patch(jgc)%nblks_c) )


    ALLOCATE(vn_lp      (nproma, nlev_p,      p_pp%nblks_e),          &
             w_lp       (nproma, nlev_p+1,    p_pp%nblks_c),          &
             thv_pr_lp  (nproma, nlev_p,      p_pp%nblks_c),          &
             rho_pr_lp  (nproma, nlev_p,      p_pp%nblks_c),          &
             tracer_lp  (nproma, nlev_p,      p_pp%nblks_c, ntracer), &
             lndvars_lp (nproma, num_lndvars, p_pp%nblks_c),          &
             snowvars_lp(nproma, num_snowvars,p_pp%nblks_c)           )

    IF (p_test_run) THEN
      lndvars_chi = 0._wp
      snowvars_chi = 0._wp
    ENDIF

    ! Step 1: boundary interpolation
    !
    ! Note: the sepration between boundary interpolation and the interpolation
    ! of the fields in the prognostic region of the nest is needed because the
    ! existing interpolation routines are runtime-optimized for operations that
    ! are needed regularly. Boundary interpolation is executed from the parent grid
    ! to the child grid, whereas the so-called nudging interpolation works from
    ! the local parent grid to the child grid.


    ! Step 1a: prepare boundary interpolation

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_patch(jg)%cells%start_blk(grf_bdywidth_c+1,1)
    i_endblk   = p_patch(jg)%cells%end_blk(min_rlcell_int-2,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

      ! Compute perturbation quantities for rho and theta at parent level
      DO jk = 1, nlev_p
        DO jc = i_startidx, i_endidx
          rho_pr_par(jc,jk,jb) = &
            p_parent_prog%rho(jc,jk,jb) - p_parent_metrics%rho_ref_mc(jc,jk,jb)

          thv_pr_par(jc,jk,jb) = &
            p_parent_prog%theta_v(jc,jk,jb) - p_parent_metrics%theta_ref_mc(jc,jk,jb)
        ENDDO
      ENDDO

      ! Reference temperature at lowest level (taken as proxy for the surface)
      DO jc = i_startidx, i_endidx
        tsfc_ref_p(jc,jb) = p_parent_metrics%theta_ref_mc(jc,nlev_p,jb) * &
                            p_parent_metrics%exner_ref_mc(jc,nlev_p,jb)
      ENDDO

      IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
        ! Collect soil variables
        DO jk = 1, nlev_soil+1
          jk1 = 3*(jk-1) + 1
          DO jc = i_startidx, i_endidx
            lndvars_par(jc,jk1,jb)   = p_parent_ldiag%t_so(jc,jk,jb) - tsfc_ref_p(jc,jb)
            lndvars_par(jc,jk1+1,jb) = p_parent_ldiag%w_so(jc,jk,jb)
            lndvars_par(jc,jk1+2,jb) = p_parent_ldiag%w_so_ice(jc,jk,jb)
          ENDDO
        ENDDO

        jk1 = 3*(nlev_soil+1)+1
        DO jc = i_startidx, i_endidx
          lndvars_par(jc,jk1,jb)   = p_parent_ldiag%t_so(jc,nlev_soil+2,jb)-tsfc_ref_p(jc,jb)
          lndvars_par(jc,jk1+1,jb) = p_parent_lprog%t_g(jc,jb)    - tsfc_ref_p(jc,jb)
          lndvars_par(jc,jk1+2,jb) = p_parent_ldiag%t_s(jc,jb)    - tsfc_ref_p(jc,jb)
          lndvars_par(jc,jk1+3,jb) = p_parent_ldiag%t_snow(jc,jb) - tsfc_ref_p(jc,jb)
          lndvars_par(jc,jk1+4,jb) = p_parent_ldiag%w_snow(jc,jb)
          lndvars_par(jc,jk1+5,jb) = p_parent_ldiag%rho_snow(jc,jb)
          lndvars_par(jc,jk1+6,jb) = p_parent_ldiag%w_i(jc,jb)
          lndvars_par(jc,jk1+7,jb) = p_parent_ldiag%freshsnow(jc,jb)
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lmulti_snow) THEN
        DO jk = 1, nlev_snow
          jk1 = 5*(jk-1) + 1
          DO jc = i_startidx, i_endidx
            snowvars_par(jc,jk1,jb)   = p_parent_ldiag%t_snow_mult(jc,jk,jb)-tsfc_ref_p(jc,jb)
            snowvars_par(jc,jk1+1,jb) = p_parent_ldiag%rho_snow_mult(jc,jk,jb)
            snowvars_par(jc,jk1+2,jb) = p_parent_ldiag%wliq_snow(jc,jk,jb)
            snowvars_par(jc,jk1+3,jb) = p_parent_ldiag%wtot_snow(jc,jk,jb)
            snowvars_par(jc,jk1+4,jb) = p_parent_ldiag%dzh_snow(jc,jk,jb)
          ENDDO
        ENDDO
        DO jc = i_startidx, i_endidx
          snowvars_par(jc,jk1+5,jb) = p_parent_ldiag%t_snow_mult(jc,nlev_snow+1,jb) - &
                                      tsfc_ref_p(jc,jb)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Step 1b: execute boundary interpolation

    CALL interpol2_vec_grf (p_patch(jg), p_pc, p_int_state(jg), &
      p_grf_state(jg)%p_dom(i_chidx), i_chidx,                  &
      p_parent_prog%vn, p_child_prog%vn)

    CALL interpol_scal_grf (p_patch(jg), p_pc, p_int_state(jg), &
      p_grf_state(jg)%p_dom(i_chidx), i_chidx, 3,               &
      rho_pr_par,      p_child_prog%rho,                        &
      thv_pr_par,      p_child_prog%theta_v,                    &
      p_parent_prog%w, p_child_prog%w                           )

    IF (ltransport) THEN
      l_limit(:) = .TRUE. ! apply positive definite limiter on tracers

      CALL interpol_scal_grf ( p_patch(jg), p_pc, p_int_state(jg),        &
        p_grf_state(jg)%p_dom(i_chidx), i_chidx, ntracer,                 &
        f4din1=p_parent_prog_rcf%tracer, f4dout1=p_child_prog_rcf%tracer, &
        llimit_nneg=l_limit)
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lmulti_snow) THEN
      CALL interpol_scal_grf (p_patch(jg), p_pc, p_int_state(jg),            &
        p_grf_state(jg)%p_dom(i_chidx), i_chidx, 2,                          &
        lndvars_par, lndvars_chi, snowvars_par, snowvars_chi, lnoshift=.TRUE.)
    ELSE IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
      CALL interpol_scal_grf (p_patch(jg), p_pc, p_int_state(jg), &
        p_grf_state(jg)%p_dom(i_chidx), i_chidx, 1,               &
        lndvars_par, lndvars_chi, lnoshift=.TRUE.                 )
    ENDIF

    ! Step 2: Interpolation of fields in the model interior

    ! Step 2a: Copy prognostic variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    IF (l_parallel) THEN

      CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 3, 3*nlev_p+1, &
                              RECV1=rho_pr_lp, SEND1=rho_pr_par,         &
                              RECV2=thv_pr_lp, SEND2=thv_pr_par,         &
                              RECV3=w_lp,      SEND3=p_parent_prog%w     )

      CALL exchange_data(p_pp%comm_pat_glb_to_loc_e, RECV=vn_lp, SEND=p_parent_prog%vn)

      CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, ntracer, ntracer*nlev_p, &
                              RECV4D=tracer_lp, SEND4D=p_parent_prog_rcf%tracer    )

      IF (atm_phy_nwp_config(jg)%inwp_surface == 1) &
        CALL exchange_data(p_pp%comm_pat_glb_to_loc_c, RECV=lndvars_lp, SEND=lndvars_par)

      IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lmulti_snow) &
        CALL exchange_data(p_pp%comm_pat_glb_to_loc_c, RECV=snowvars_lp, SEND=snowvars_par)

    ELSE
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      ! a) cell-based variables
      i_startblk = p_pp%cells%start_blk(grf_nudgintp_start_c+1,i_chidx)
      i_endblk   = p_pp%cells%end_blk(min_rlcell_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
                           grf_nudgintp_start_c+1, min_rlcell_int, i_chidx)

        ! Note: w(nlevp1) is diagnostic and therefore not nudged
        DO jk = 1, nlev_p
          DO jc = i_startidx, i_endidx
            rho_pr_lp(jc,jk,jb) = rho_pr_par(jc,jk,jb)
            thv_pr_lp(jc,jk,jb) = thv_pr_par(jc,jk,jb)
            w_lp(jc,jk,jb)      = p_parent_prog%w(jc,jk,jb)
          ENDDO
        ENDDO

        DO jc = i_startidx, i_endidx
          w_lp(jc,nlev_p+1,jb)      = p_parent_prog%w(jc,nlev_p+1,jb)
        ENDDO

        IF (ltransport) THEN
          DO jt = 1, ntracer
            DO jk = 1, nlev_p
              DO jc = i_startidx, i_endidx
                tracer_lp(jc,jk,jb,jt) = p_parent_prog_rcf%tracer(jc,jk,jb,jt)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
          DO jk = 1, num_lndvars
            DO jc = i_startidx, i_endidx
              lndvars_lp(jc,jk,jb) = lndvars_par(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF

        IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lmulti_snow) THEN
          DO jk = 1, num_snowvars
            DO jc = i_startidx, i_endidx
              snowvars_lp(jc,jk,jb) = snowvars_par(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF

      ENDDO
!$OMP END DO

      ! b) velocity
      i_startblk = p_pp%edges%start_blk(grf_nudgintp_start_e+2,i_chidx)
      i_endblk   = p_pp%edges%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
                           grf_nudgintp_start_e+2, min_rledge_int, i_chidx)

        DO jk = 1, nlev_p
          DO je = i_startidx, i_endidx
            vn_lp(je,jk,jb) = p_parent_prog%vn(je,jk,jb)
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF

    ! Step 2b: Perform interpolation from local parent to child grid

    ! Note: the sync routines cannot be used for the synchronization on the local
    ! parent grid

    IF(l_parallel) CALL exchange_data(p_pp%comm_pat_e, vn_lp)
    CALL interpol_vec_nudging (p_pp, p_pc, p_int, p_grf%p_dom(i_chidx), p_grfc, &
                               i_chidx, nshift, 1, vn_lp, p_child_prog%vn       )
    CALL sync_patch_array(SYNC_E,p_pc,p_child_prog%vn)

    IF(l_parallel) CALL exchange_data_mult(p_pp%comm_pat_c, 3, 3*nlev_p+1, &
                               recv1=thv_pr_lp, recv2=rho_pr_lp, recv3=w_lp)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, nshift, 3, 1, &
                                f3din1=thv_pr_lp, f3dout1=p_child_prog%theta_v,           &
                                f3din2=rho_pr_lp, f3dout2=p_child_prog%rho,               &
                                f3din3=w_lp,      f3dout3=p_child_prog%w                  )
    CALL sync_patch_array_mult(SYNC_C,p_pc,3,p_child_prog%theta_v,p_child_prog%rho,  &
                               p_child_prog%w)

    IF (ltransport) THEN
      IF(l_parallel) CALL exchange_data_mult(p_pp%comm_pat_c, ntracer, ntracer*nlev_p, &
                                             recv4d=tracer_lp)
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, &
                                  nshift, ntracer, 1, f4din=tracer_lp,        &
                                  f4dout=p_child_prog_rcf%tracer              )
      CALL sync_patch_array_mult(SYNC_C,p_pc,ntracer,f4din=p_child_prog_rcf%tracer)
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
      IF(l_parallel) CALL exchange_data(p_pp%comm_pat_c, lndvars_lp)
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, &
                                  1, 1, f3din1=lndvars_lp, f3dout1=lndvars_chi )
      CALL sync_patch_array(SYNC_C,p_pc,lndvars_chi)
    ENDIF

    IF (atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. lmulti_snow) THEN
      IF(l_parallel) CALL exchange_data(p_pp%comm_pat_c, snowvars_lp)
      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, &
                                  1, 1, f3din1=snowvars_lp, f3dout1=snowvars_chi )
      CALL sync_patch_array(SYNC_C,p_pc,snowvars_chi)
    ENDIF

    ! Step 3: Add reference state to thermodynamic variables and copy land fields
    ! from the container arrays to the prognostic variables (for the time being,
    ! all tiles are initialized with the interpolated aggregated variables)
    ! Note that nest boundary points have to be included here to obtain a proper
    ! initialization of all grid points

    i_nchdom = MAX(1,p_pc%n_childdom)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_pc%cells%start_blk(1,1)
    i_endblk   = p_pc%cells%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell)

      ! rho and theta currently contain perturbation fields; thus, the reference
      ! state still needs to be added. Then, rhotheta and exner are diagnosed
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
          p_child_prog%rho(jc,jk,jb) = &
            p_child_prog%rho(jc,jk,jb) + p_child_metrics%rho_ref_mc(jc,jk,jb)

          p_child_prog%theta_v(jc,jk,jb) = &
            p_child_prog%theta_v(jc,jk,jb) + p_child_metrics%theta_ref_mc(jc,jk,jb)

          p_child_prog%rhotheta_v(jc,jk,jb) = &
            p_child_prog%rho(jc,jk,jb) * p_child_prog%theta_v(jc,jk,jb)

          p_child_prog%exner(jc,jk,jb) = EXP(rd_o_cvd*LOG(rd/p0ref* &
            p_child_prog%rhotheta_v(jc,jk,jb)))
        ENDDO
      ENDDO

      ! Reference temperature at lowest level (taken as proxy for the surface)
      DO jc = i_startidx, i_endidx
        tsfc_ref_c(jc,jb) = p_child_metrics%theta_ref_mc(jc,nlev_c,jb) * &
                            p_child_metrics%exner_ref_mc(jc,nlev_c,jb)
      ENDDO

      IF (atm_phy_nwp_config(jgc)%inwp_surface == 1) THEN
        ! Distribute soil variables
        DO jt = 1, nsfc_subs
          DO jk = 1, nlev_soil+1
            jk1 = 3*(jk-1) + 1
            DO jc = i_startidx, i_endidx
              p_child_lprog%t_so_t(jc,jk,jb,jt) = lndvars_chi(jc,jk1,jb) + tsfc_ref_c(jc,jb)
              p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(0._wp,lndvars_chi(jc,jk1+1,jb))
              ! limit w_so_t to pore volume and dryness point - TERRA crashes otherwise
              IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 3) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.012_wp, &
                  MIN(dzsoil(jk)*0.364_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ELSE IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 4) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.030_wp, &
                  MIN(dzsoil(jk)*0.445_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ELSE IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 5) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.035_wp, &
                  MIN(dzsoil(jk)*0.455_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ELSE IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 6) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.060_wp, &
                  MIN(dzsoil(jk)*0.475_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ELSE IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 7) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.065_wp, &
                  MIN(dzsoil(jk)*0.507_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ELSE IF (ext_data(jgc)%atm%soiltyp(jc,jb) == 8) THEN
                p_child_lprog%w_so_t(jc,jk,jb,jt) = MAX(dzsoil(jk)*0.098_wp, &
                  MIN(dzsoil(jk)*0.863_wp,p_child_lprog%w_so_t(jc,jk,jb,jt)))
              ENDIF
              p_child_lprog%w_so_ice_t(jc,jk,jb,jt) = MAX(0._wp,lndvars_chi(jc,jk1+2,jb))
            ENDDO
          ENDDO
        ENDDO

        jk1 = 3*(nlev_soil+1)+1
        DO jt = 1, nsfc_subs
          DO jc = i_startidx, i_endidx
            p_child_lprog%t_so_t(jc,nlev_soil+2,jb,jt) = lndvars_chi(jc,jk1,jb) + &
              tsfc_ref_c(jc,jb)
            ! here we need to initialize t_g rather than t_g_t because t_g is copied to t_g_t
            ! in nwp_surface_init
            p_child_lprog%t_g(jc,jb) = lndvars_chi(jc,jk1+1,jb) + tsfc_ref_c(jc,jb)
            p_child_lprog%t_s_t(jc,jb,jt) = lndvars_chi(jc,jk1+2,jb) + tsfc_ref_c(jc,jb)
            p_child_lprog%t_snow_t(jc,jb,jt) = lndvars_chi(jc,jk1+3,jb) + tsfc_ref_c(jc,jb)
            p_child_lprog%w_snow_t(jc,jb,jt) = MAX(0._wp,lndvars_chi(jc,jk1+4,jb))
            p_child_lprog%rho_snow_t(jc,jb,jt) = lndvars_chi(jc,jk1+5,jb)
            p_child_lprog%w_i_t(jc,jb,jt) = MAX(0._wp,lndvars_chi(jc,jk1+6,jb))
            p_child_ldiag%freshsnow_t(jc,jb,jt) = MAX(0._wp,MIN(1._wp,lndvars_chi(jc,jk1+7,jb)))
          ENDDO
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jgc)%inwp_surface == 1 .AND. lmulti_snow) THEN
        DO jt = 1, nsfc_subs
          DO jk = 1, nlev_snow
            jk1 = 5*(jk-1) + 1
            DO jc = i_startidx, i_endidx
              p_child_lprog%t_snow_mult_t(jc,jk,jb,jt) = snowvars_chi(jc,jk1,jb) + &
                tsfc_ref_c(jc,jb)
              p_child_lprog%rho_snow_mult_t(jc,jk,jb,jt) = snowvars_chi(jc,jk1+1,jb)
              p_child_lprog%wliq_snow_t(jc,jk,jb,jt)     = MAX(0._wp,snowvars_chi(jc,jk1+2,jb))
              p_child_lprog%wtot_snow_t(jc,jk,jb,jt)     = MAX(0._wp,snowvars_chi(jc,jk1+3,jb))
              p_child_lprog%dzh_snow_t(jc,jk,jb,jt)      = MAX(0._wp,snowvars_chi(jc,jk1+4,jb))
            ENDDO
          ENDDO
          DO jc = i_startidx, i_endidx
            p_child_lprog%t_snow_mult_t(jc,nlev_snow+1,jb,jt) = snowvars_chi(jc,jk1+5,jb)+ &
              tsfc_ref_c(jc,jb)
          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DEALLOCATE(thv_pr_par, rho_pr_par, lndvars_par, snowvars_par, lndvars_chi , snowvars_chi,    &
               tsfc_ref_p, tsfc_ref_c, vn_lp, w_lp, thv_pr_lp, rho_pr_lp, tracer_lp, lndvars_lp, &
               snowvars_lp)

  END SUBROUTINE initialize_nest


  !---------------------------------------------------------------------------
  !>
  !! Computes topography blending for nested domains initialized with real 
  !! topography data
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2011-07-05)
  !!
  SUBROUTINE topography_blending(p_pp, p_pc, p_intp, p_intc, p_grf, i_chidx, &
               topo_cp, topo_cc, topo_vc)

    ! patch at parent and child level
    TYPE(t_patch),                TARGET, INTENT(IN) :: p_pp, p_pc
    ! interpolation state at parent and child level
    TYPE(t_int_state),            TARGET, INTENT(IN) :: p_intp, p_intc
    ! grf state is needed at parent level only
    TYPE(t_gridref_single_state), TARGET, INTENT(IN) :: p_grf

    ! child domain index
    INTEGER, INTENT(IN) :: i_chidx

    ! topography on cells and vertices: p = parent level, c = child level
    REAL(wp), INTENT(IN),     DIMENSION(:,:) :: topo_cp
    REAL(wp), INTENT(INOUT),  DIMENSION(:,:) :: topo_cc, topo_vc

    ! auxiliary fields needed for SR calls (SR's expect 3D fields)
    REAL(wp), TARGET :: z_topo_cp(nproma,1,p_pp%nblks_c)
    REAL(wp)         :: z_topo_cc(nproma,1,p_pc%nblks_c)
    REAL(wp)         :: z_topo_vc(nproma,1,p_pc%nblks_v)

    ! another auxiliary needed to handle MPI parallelization, and a related pointer
    REAL(wp), ALLOCATABLE, TARGET :: z_topo_clp (:,:,:)
    REAL(wp), POINTER             :: ptr_topo_cp(:,:,:)

    TYPE(t_gridref_single_state), POINTER :: ptr_grf
    TYPE(t_patch),                POINTER :: ptr_pp
    TYPE(t_int_state),            POINTER :: ptr_int

    INTEGER  :: jgc, jb, jc, jv, nblks_c
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom, i_rlstart, i_rlend
    REAL(wp) :: wfac
    LOGICAL  :: l_parallel

    !-------------------------------------------------------------------------

    jgc = p_pc%id

!    IF (p_nprocs == 1 .OR. p_pe == p_test_pe) THEN
    IF(.NOT. my_process_is_mpi_parallel()) THEN
      l_parallel = .FALSE.
      ptr_grf    => p_grf
      ptr_pp     => p_pp
      ptr_int    => p_intp
    ELSE
      l_parallel = .TRUE.
      ptr_grf    => p_grf_state_local_parent(jgc)%p_dom(i_chidx)
      ptr_pp     => p_patch_local_parent(jgc)
      ptr_int    => p_int_state_local_parent(jgc)
      ALLOCATE(z_topo_clp(nproma,1,ptr_pp%nblks_c))
    ENDIF

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(2(a,i2))') 'topography blending, domain ',&
        p_pp%id,'  => domain ',p_pc%id
      CALL message('', TRIM(message_text))
    ENDIF

    z_topo_cc(:,:,:) = 0._wp
    z_topo_vc(:,:,:) = 0._wp

    nblks_c = p_pp%nblks_c

    ! 1.(a) Copy input topography to auxiliary 3D field
    DO jb = 1, nblks_c
      CALL get_indices_c(p_pp, jb, 1, nblks_c, i_startidx, i_endidx, 1)

      DO jc = i_startidx, i_endidx
        z_topo_cp(jc,1,jb) = topo_cp(jc,jb)
      ENDDO
    ENDDO

    ! 1.(b) Copy this auxiliary field to the local parent in case of MPI parallelization
    IF (l_parallel) THEN
      CALL exchange_data(ptr_pp%comm_pat_glb_to_loc_c, RECV=z_topo_clp, SEND=z_topo_cp)
      ptr_topo_cp => z_topo_clp
      ! synchronization (CALL sync does not work on local parent)
      CALL exchange_data(ptr_pp%comm_pat_c, ptr_topo_cp)
    ELSE
      ptr_topo_cp => z_topo_cp
    ENDIF

    ! 1.(c) Interpolate coarse topography on fine mesh

    ! Comment: one could avoid the awkward procedure of interpolating boundary zone and
    ! prognostic zone separately by introducing appropriate interpolation routines and communication
    ! patterns (both of which are optimized for runtime tasks). However, those routines
    ! would not be needed anywhere else, and terrain blending is not runtime-critical

    ! Lateral boundary zone
    CALL interpol_scal_grf (p_pp, p_pc, p_intp, p_grf, i_chidx, 1,              &
                            f3din1=z_topo_cp, f3dout1=z_topo_cc, lnoshift=.TRUE.)
    ! Prognostic part of the model domain
    ! Note: in contrast to boundary interpolation, nudging expects the input on the local parent grid
    CALL interpol_scal_nudging (ptr_pp, ptr_int, ptr_grf, i_chidx, 0, 1, 1, &
                                f3din1=ptr_topo_cp, f3dout1=z_topo_cc       )


    ! 2. Apply terrain blending to cell points
    ! In the boundary interpolation zone (cell rows 1-4), the interpolated values
    ! are taken; in the subsequent 8 rows, linear blending is applied


    i_nchdom = MAX(1,p_pc%n_childdom)

    ! 2.(a) Pure interpolation in boundary interpolation zone
    i_rlstart  = 1
    i_rlend    = grf_bdywidth_c
    i_startblk = p_pc%cells%start_blk(i_rlstart,1)
    i_endblk   = p_pc%cells%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        topo_cc(jc,jb) = z_topo_cc(jc,1,jb)
      ENDDO
    ENDDO

    ! 2.(b) Blending in prognostic part of the nested domain
    i_rlstart  = grf_bdywidth_c + 1
    i_rlend    = min_rlcell
    i_startblk = p_pc%cells%start_blk(i_rlstart,1)
    i_endblk   = p_pc%cells%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        IF (p_pc%cells%refin_ctrl(jc,jb) > 0 .AND. &
            p_pc%cells%refin_ctrl(jc,jb) <= grf_bdywidth_c + nudge_zone_width) THEN

          wfac = REAL(p_pc%cells%refin_ctrl(jc,jb) - grf_bdywidth_c,wp)/ &
            REAL(nudge_zone_width+1,wp)
          topo_cc(jc,jb) = wfac*topo_cc(jc,jb) + (1._wp-wfac)*z_topo_cc(jc,1,jb)

        ENDIF ! on the interior grid points, topo_cc remains unchanged

      ENDDO
    ENDDO

    CALL sync_patch_array(SYNC_C,p_pc,topo_cc)

    ! 3.(a) Copy blended topography to z_topo_cc to prepare interpolation to vertices
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_pc%cells%start_blk(i_rlstart,1)
    i_endblk   = p_pc%cells%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_pc, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        z_topo_cc(jc,1,jb) = topo_cc(jc,jb) 

      ENDDO
    ENDDO

    ! 3.(b)Interpolate blended topography from cells to vertices
    CALL cells2verts_scalar(z_topo_cc,p_pc,p_intc%cells_aw_verts,z_topo_vc,1,1)


    ! 3.(c) Use interpolated vertex topography in boundary interpolation zone
    i_rlstart  = 1
    i_rlend    = grf_bdywidth_c ! grf_bdywidth_v does not exist, but indexing for vertices is the same as for cells
    i_startblk = p_pc%verts%start_blk(i_rlstart,1)
    i_endblk   = p_pc%verts%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_v(p_pc, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jv = i_startidx, i_endidx
        topo_vc(jv,jb) = z_topo_vc(jv,1,jb)
      ENDDO
    ENDDO

    ! 3.(d) Apply blending in the grid rows for which blending has been applied to cell points
    i_rlstart  = grf_bdywidth_c + 1
    i_rlend    = min_rlvert
    i_startblk = p_pc%verts%start_blk(i_rlstart,1)
    i_endblk   = p_pc%verts%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_v(p_pc, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jv = i_startidx, i_endidx

        IF (p_pc%verts%refin_ctrl(jv,jb) > 0 .AND. &
            p_pc%verts%refin_ctrl(jv,jb) <= grf_bdywidth_c + nudge_zone_width) THEN

          wfac = REAL(p_pc%verts%refin_ctrl(jv,jb) - grf_bdywidth_c,wp)/ &
            REAL(nudge_zone_width+1,wp)
          topo_vc(jv,jb) = wfac*topo_vc(jv,jb) + (1._wp-wfac)*z_topo_vc(jv,1,jb)

        ENDIF ! on the interior grid points, topo_cc remains unchanged

      ENDDO
    ENDDO

    CALL sync_patch_array(SYNC_V,p_pc,topo_vc)

    IF (l_parallel) DEALLOCATE(z_topo_clp)

  END SUBROUTINE topography_blending


  !---------------------------------------------------------------------------
  !>
  !! Computes topography feedback for nested domains initialized with real 
  !! topography data
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2011-07-05)
  !!
  SUBROUTINE topography_feedback(p_pp, p_int, p_grf, i_chidx, &
               topo_cp, topo_cc, topo_vp)

    ! patch at parent level
    TYPE(t_patch),                TARGET, INTENT(IN) :: p_pp
    ! interpolation state at parent level
    TYPE(t_int_state),            TARGET, INTENT(IN) :: p_int
    ! grf state is needed at parent level only
    TYPE(t_gridref_state), TARGET, INTENT(IN) :: p_grf

    ! child domain index
    INTEGER, INTENT(IN) :: i_chidx

    ! topography on cells and vertices: p = parent level, c = child level
    REAL(wp), INTENT(IN),    DIMENSION(:,:)         :: topo_cc
    REAL(wp), INTENT(INOUT), DIMENSION(:,:), TARGET :: topo_cp, topo_vp

    ! auxiliary fields needed for SR calls (SR's expect 3D fields)
    REAL(wp) :: z_topo_cp(nproma,1,p_pp%nblks_c)
    REAL(wp) :: z_topo_vp(nproma,1,p_pp%nblks_v)

    ! another auxiliary to handle MPI parallelization, and related pointer
    REAL(wp), ALLOCATABLE, TARGET :: z_topo_clp(:,:)
    REAL(wp), POINTER             :: ptr_topo_cp(:,:)

    TYPE(t_gridref_state), POINTER :: ptr_grf
    TYPE(t_patch),         POINTER :: ptr_pp

    INTEGER  :: jgc, jb, jc, jv
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom, i_rlstart, i_rlend

    INTEGER, POINTER, DIMENSION(:,:,:) :: iidx, iblk

    LOGICAL :: l_parallel

    !-------------------------------------------------------------------------

    jgc = p_pp%child_id(i_chidx) ! child domain ID

!    IF (p_nprocs == 1 .OR. p_pe == p_test_pe) THEN
      IF(.NOT. my_process_is_mpi_parallel()) THEN
      l_parallel  =  .FALSE.
      ptr_pp      => p_pp
      ptr_grf     => p_grf
      ptr_topo_cp => topo_cp
    ELSE
      l_parallel = .TRUE.
      ALLOCATE(z_topo_clp(nproma,p_patch_local_parent(jgc)%nblks_c))
      ptr_pp      => p_patch_local_parent(jgc)
      ptr_grf     => p_grf_state_local_parent(jgc)
      ptr_topo_cp => z_topo_clp
    ENDIF

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(2(a,i2))') 'topography feedback, domain ',&
        p_pp%child_id(i_chidx),'  => domain ',p_pp%id
        CALL message('', TRIM(message_text))
    ENDIF

    ! Start/End block in the parent domain
    i_rlstart  = grf_fbk_start_c
    i_rlend    = min_rlcell_int

    i_startblk = ptr_pp%cells%start_blk(i_rlstart,i_chidx)
    i_endblk   = ptr_pp%cells%end_blk(i_rlend,i_chidx)

    iidx => ptr_pp%cells%child_idx
    iblk => ptr_pp%cells%child_blk

    ! Perform feedback for cell points in the prognostic part of the nest overlap area
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk,                &
                         i_startidx, i_endidx, i_rlstart, i_rlend, i_chidx)

      DO jc = i_startidx, i_endidx
        ptr_topo_cp(jc,jb) =                                                &
          topo_cc(iidx(jc,jb,1),iblk(jc,jb,1))*ptr_grf%fbk_wgt_c(jc,jb,1) + &
          topo_cc(iidx(jc,jb,2),iblk(jc,jb,2))*ptr_grf%fbk_wgt_c(jc,jb,2) + &
          topo_cc(iidx(jc,jb,3),iblk(jc,jb,3))*ptr_grf%fbk_wgt_c(jc,jb,3) + &
          topo_cc(iidx(jc,jb,4),iblk(jc,jb,4))*ptr_grf%fbk_wgt_c(jc,jb,4)
      ENDDO
    ENDDO

    IF (l_parallel) THEN
      ! Attention: the feedback communication pattern is defined on the local parent (ptr_pp), ...
      CALL exchange_data(ptr_pp%comm_pat_loc_to_glb_c_fbk, RECV=topo_cp, SEND=ptr_topo_cp, &
                         l_recv_exists=.TRUE.)
    ENDIF
    ! ... whereas the subsequent halo communication needs to be done on the parent grid (p_pp)
    CALL sync_patch_array(SYNC_C,p_pp,topo_cp)


    ! Copy fed-back topography to z_topo_cp to prepare interpolation to vertices
    i_nchdom = MAX(1,p_pp%n_childdom)
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_pp%cells%start_blk(i_rlstart,1)
    i_endblk   = p_pp%cells%end_blk(i_rlend,i_nchdom)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,         &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        z_topo_cp(jc,1,jb) = topo_cp(jc,jb) 

      ENDDO
    ENDDO

    ! Interpolate to vertices
    CALL cells2verts_scalar(z_topo_cp,p_pp,p_int%cells_aw_verts,z_topo_vp,1,1)

    ! Copy back to topo_vp in nest overlap area
    i_rlstart  = grf_fbk_start_c
    i_rlend    = min_rlvert_int

    i_startblk = p_pp%verts%start_blk(i_rlstart,i_chidx)
    i_endblk   = p_pp%verts%end_blk(i_rlend,i_chidx)

    DO jb = i_startblk, i_endblk

      CALL get_indices_v(p_pp, jb, i_startblk, i_endblk,                  &
                         i_startidx, i_endidx, i_rlstart, i_rlend, i_chidx)

      DO jv = i_startidx, i_endidx
        topo_vp(jv,jb) = z_topo_vp(jv,1,jb)
      ENDDO
    ENDDO

    CALL sync_patch_array(SYNC_V,p_pp,topo_vp)

    IF (l_parallel) DEALLOCATE(z_topo_clp)

  END SUBROUTINE topography_feedback


END MODULE mo_nh_init_nest_utils

