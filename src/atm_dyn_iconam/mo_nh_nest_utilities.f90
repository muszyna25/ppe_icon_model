!>
!!  This module contains the routines needed for nesting in the nonhydrostatic.
!!  version.
!!
!! @par Revision History
!!  Developed and tested by Guenther Zaengl, DWD (2010-02-10)
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
!!
  
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_nest_utilities
  !
  !
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, message
  USE mo_model_domain,        ONLY: t_patch, t_grid_cells, t_grid_edges, p_patch_local_parent, &
    p_patch
  USE mo_grid_config,         ONLY: n_dom, n_dom_start
  USE mo_intp_data_strc,      ONLY: t_int_state, p_int_state, p_int_state_local_parent
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state, p_grf_state, p_grf_state_local_parent
  USE mo_gridref_config,      ONLY: grf_intmethod_c, grf_intmethod_e, grf_intmethod_ct
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf, interpol_vec_grf, interpol2_vec_grf
  USE mo_grf_nudgintp,        ONLY: interpol_scal_nudging, interpol_vec_nudging
  USE mo_grf_ubcintp,         ONLY: interpol_scal_ubc,interpol_vec_ubc
  USE mo_dynamics_config,     ONLY: nnow, nsav1, nnow_rcf
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: ltransport, msg_level, ntracer, lvert_nest
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydro_state,      ONLY: p_nh_state
  USE mo_nonhydrostatic_config,ONLY: iadv_rcf
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlcell_int, min_rledge_int, &
    &                           MAX_CHAR_LENGTH
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c,                       &
    grf_bdyintp_end_c,                         &
    grf_fbk_start_c,                           &
    grf_bdywidth_c, grf_bdywidth_e,            &
    grf_nudgintp_start_c, grf_nudgintp_start_e,&
    grf_nudge_start_c, grf_nudge_start_e
  USE mo_mpi,                 ONLY: my_process_is_mpi_seq
  USE mo_communication,       ONLY: exchange_data, exchange_data_mult
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, sync_patch_array, &
    global_sum_array3, sync_patch_array_mult
  USE mo_physical_constants,  ONLY: rd, cvd_o_rd, p0ref
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: compute_tendencies, boundary_interpolation, complete_nesting_setup, &
    prep_bdy_nudging, outer_boundary_nudging, nest_boundary_nudging,    &
    prep_rho_bdy_nudging, density_boundary_nudging, prep_outer_bdy_nudging

CONTAINS

  !>
  !! Computes geometric information needed for the conservation correction
  !! in the feedback routine
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD, 2010-05-05
  !!
  SUBROUTINE complete_nesting_setup


    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()
    TYPE(t_patch),      POINTER     :: p_lp => NULL()
    TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()

    REAL(wp), ALLOCATABLE :: cell_volume(:,:,:), z_rho_ref(:,:,:)

    INTEGER :: jg, ji, jgc, jb, jc, jk, jks
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER :: nlev, nlev_c, nshift

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt

    ALLOCATE(p_nh_state(1)%metrics%fbk_dom_volume(p_patch(1)%nlev))
    p_nh_state(1)%metrics%fbk_dom_volume(:) = 0._wp

    DO jg = 1, n_dom-1

      p_pp   => p_patch(jg)
      p_gcp  => p_pp%cells

      nlev = p_pp%nlev

      i_nchdom = p_pp%n_childdom
      IF (i_nchdom == 0) CYCLE

      DO ji = 1, i_nchdom

        jgc    =  p_pp%child_id(ji)
        p_pc   => p_patch(jgc)

        ! Note: the number of levels of fbk_dom_volume is that of the parent grid
        ALLOCATE(p_nh_state(jgc)%metrics%fbk_dom_volume(nlev))

        ! Copy layer thicknesses to local parent of current child

        p_lp => p_patch_local_parent(jgc)
        ALLOCATE(p_lp%cells%ddqz_z_full(nproma, p_lp%nlev, p_lp%n_patch_cells))
        p_lp%cells%ddqz_z_full(:,:,:) = 0._wp ! Safety only
        CALL exchange_data(p_lp%comm_pat_glb_to_loc_c, p_lp%cells%ddqz_z_full, &
          &                p_nh_state(jg)%metrics%ddqz_z_full)

        i_startblk = p_gcp%start_blk(1,1)
        i_endblk   = p_gcp%end_blk(min_rlcell_int,i_nchdom)

        ALLOCATE(cell_volume(nproma, nlev, i_endblk))

        cell_volume = 0._wp

        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            1, min_rlcell_int)

          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              ! Sum must be taken over owned points of nest overlap region only
              IF (p_gcp%refin_ctrl(jc,jb) <= grf_fbk_start_c .AND. p_gcp%child_id(jc,jb) == jgc .AND. &
                p_gcp%decomp_info%owner_mask(jc,jb)) &
                cell_volume(jc,jk,jb) = p_gcp%area(jc,jb)*p_nh_state(jg)%metrics%ddqz_z_full(jc,jk,jb)
            ENDDO
          ENDDO

        ENDDO

        p_nh_state(jgc)%metrics%fbk_dom_volume(:) = global_sum_array3(1,.FALSE.,cell_volume)

        DEALLOCATE(cell_volume)

      ENDDO

      ! Second part: compute correction term needed to use perturbation density for boundary nudging 

      DO ji = 1, i_nchdom

        jgc    =  p_pp%child_id(ji)
        p_pc   => p_patch(jgc)

        nlev_c = p_pc%nlev
        nshift = p_pc%nshift

        p_fbkwgt => p_grf_state_local_parent(jgc)%fbk_wgt_c
        p_gcp => p_patch_local_parent(jgc)%cells
        p_pp  => p_patch_local_parent(jgc)

        iidx  => p_gcp%child_idx
        iblk  => p_gcp%child_blk

        i_startblk = p_gcp%start_blk(grf_nudgintp_start_c+1,ji)
        i_endblk   = p_gcp%end_blk(min_rlcell_int,ji)

        ALLOCATE(p_nh_state(jgc)%metrics%rho_ref_corr(nproma, nlev_c, p_pp%nblks_c), &
          &      z_rho_ref(nproma, nlev, p_pp%nblks_c))
        z_rho_ref(:,:,:) = 0._wp
        p_nh_state(jgc)%metrics%rho_ref_corr(:,:,:) = 0._wp

        CALL exchange_data(p_pp%comm_pat_glb_to_loc_c, RECV=z_rho_ref, &
          &                SEND=p_nh_state(jg)%metrics%rho_ref_mc)

        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            grf_nudgintp_start_c+1, min_rlcell_int, ji)

          DO jk = 1, nlev_c
            jks = jk + nshift
            DO jc = i_startidx, i_endidx

              p_nh_state(jgc)%metrics%rho_ref_corr(jc,jk,jb) = - z_rho_ref(jc,jks,jb) + (           &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1)+ &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2)+ &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3)+ &
                p_nh_state(jgc)%metrics%rho_ref_mc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)  )

            ENDDO
          ENDDO
        ENDDO

        DEALLOCATE(z_rho_ref)

      ENDDO
    ENDDO

  END SUBROUTINE complete_nesting_setup


  !-------------------------------------------------------------------------
  !
  !
  !
  !>
  !! Computes the time tendencies of the prognostic variables needed for
  !! interpolation to the lateral boundaries of the nested domains
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD, 2010-02-10
  !!
  SUBROUTINE compute_tendencies (jg,n_new,n_now,n_new_rcf,n_now_rcf,&
    &                            rdt,rdt_rcf,rdt_mflx,lstep_adv)


    INTEGER,  INTENT(IN) :: jg  ! domain ID
    ! Time levels from which tendencies are computed
    INTEGER,  INTENT(IN) ::  n_new,n_now
    ! Time levels from which tracer-tendencies are computed
    INTEGER,  INTENT(IN) ::  n_new_rcf,n_now_rcf

    ! Inverse value of time step needed for computing the tendencies
    REAL(wp), INTENT(IN) ::  rdt
    ! Inverse value of time step for integration with reduced calling frequency,
    ! needed for computing the tracer-tendencies
    REAL(wp), INTENT(IN) ::  rdt_rcf
    ! Inverse value of time step needed for computing the mass flux tendencies
    REAL(wp), INTENT(IN) ::  rdt_mflx

    LOGICAL :: lstep_adv  ! determines wheter tracer-tendencies should be computed
    ! (.true.) or not (.false.)

    ! local variables

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx,       &
      jb, jc, je, jk, jt, js, i_nchdom, nshift
    INTEGER :: nlev, nlevp1           !< number of full and half levels
    REAL(wp) :: rdt_ubc, dthalf
    ! Switch to control if the child domain is vertically nested and therefore 
    ! needs interpolation of upper boundary conditions
    LOGICAL :: l_child_vertnest

    TYPE(t_nh_state), POINTER :: p_nh
    TYPE(t_nh_prog),  POINTER :: p_prog_now
    TYPE(t_nh_prog),  POINTER :: p_prog_new
    TYPE(t_nh_prog),  POINTER :: p_prog_now_rcf
    TYPE(t_nh_prog),  POINTER :: p_prog_new_rcf

    !-----------------------------------------------------------------------

    i_nchdom = MAX(1,p_patch(jg)%n_childdom)

    p_nh           => p_nh_state(jg)
    p_prog_now     => p_nh%prog(n_now)
    p_prog_new     => p_nh%prog(n_new)
    p_prog_now_rcf => p_nh%prog(n_now_rcf)
    p_prog_new_rcf => p_nh%prog(n_new_rcf)

    ! number of vertical levels
    nlev   = p_patch(jg)%nlev
    nlevp1 = p_patch(jg)%nlevp1

    ! determine if upper boundary interpolation is needed
    IF (lvert_nest .AND. (p_patch(jg)%nshift_child > 0)) THEN  
      l_child_vertnest = .TRUE.
      nshift = p_patch(jg)%nshift_child + 1
    ELSE
      l_child_vertnest = .FALSE.
      nshift = 0
    ENDIF
    rdt_ubc = rdt*REAL(iadv_rcf,wp)/REAL(2*(MAX(1,iadv_rcf-2)),wp)
    dthalf = 1._wp/(2._wp*rdt)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! cell-based variables
    i_startblk = p_patch(jg)%cells%start_blk(grf_bdywidth_c+1,1)
    i_endblk   = p_patch(jg)%cells%end_blk(min_rlcell_int-2,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int-2)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          p_nh%diag%grf_tend_rho(jc,jk,jb) = &
            ( p_prog_new%rho(jc,jk,jb) - p_prog_now%rho(jc,jk,jb) )*rdt

          p_nh%diag%grf_tend_thv(jc,jk,jb) = &
            ( p_prog_new%theta_v(jc,jk,jb) - p_prog_now%theta_v(jc,jk,jb) )*rdt

          ! the div field carries perturbation density for use in SR boundary_interpolation
          p_nh%diag%div(jc,jk,jb) = &
            p_prog_now%rho(jc,jk,jb) - p_nh%metrics%rho_ref_mc(jc,jk,jb)

          ! the dpres_mc field carries perturbation potential temperature for use in SR boundary_interpolation
          p_nh%diag%dpres_mc(jc,jk,jb) = &
            p_prog_now%theta_v(jc,jk,jb) - p_nh%metrics%theta_ref_mc(jc,jk,jb)

          p_nh%diag%grf_tend_w(jc,jk,jb) = &
            ( p_prog_new%w(jc,jk,jb) - p_prog_now%w(jc,jk,jb) )*rdt
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        p_nh%diag%grf_tend_w(jc,nlevp1,jb) = &
          ( p_prog_new%w(jc,nlevp1,jb) - p_prog_now%w(jc,nlevp1,jb) )*rdt
      ENDDO

      IF (l_child_vertnest) THEN ! Compute upper boundary condition and its time derivative for nested domain
        p_nh%diag%dw_int         (:,jb,iadv_rcf+1)            = 0._wp
        p_nh%diag%mflx_ic_int    (:,jb,iadv_rcf+1:iadv_rcf+2) = 0._wp
        p_nh%diag%dtheta_v_ic_int(:,jb,iadv_rcf+1)            = 0._wp
        DO js = 1, iadv_rcf
          DO jc = i_startidx, i_endidx
            p_nh%diag%dw_int(jc,jb,iadv_rcf+1)          = p_nh%diag%dw_int(jc,jb,iadv_rcf+1) +          &
              p_nh%diag%dw_int(jc,jb,js)
            p_nh%diag%mflx_ic_int(jc,jb,iadv_rcf+1)     = p_nh%diag%mflx_ic_int(jc,jb,iadv_rcf+1) +     &
              p_nh%diag%mflx_ic_int(jc,jb,js)
            p_nh%diag%dtheta_v_ic_int(jc,jb,iadv_rcf+1) = p_nh%diag%dtheta_v_ic_int(jc,jb,iadv_rcf+1) + &
              p_nh%diag%dtheta_v_ic_int(jc,jb,js)
          ENDDO
        ENDDO
        DO jc = i_startidx, i_endidx
          p_nh%diag%dw_int(jc,jb,iadv_rcf+1)          = p_nh%diag%dw_int(jc,jb,iadv_rcf+1)          / REAL(iadv_rcf,wp)
          p_nh%diag%mflx_ic_int(jc,jb,iadv_rcf+1)     = p_nh%diag%mflx_ic_int(jc,jb,iadv_rcf+1)     / REAL(iadv_rcf,wp)
          p_nh%diag%dtheta_v_ic_int(jc,jb,iadv_rcf+1) = p_nh%diag%dtheta_v_ic_int(jc,jb,iadv_rcf+1) / REAL(iadv_rcf,wp)
        ENDDO
        IF (iadv_rcf >= 3) THEN
          DO jc = i_startidx, i_endidx
            ! Compute time tendency of mass flux upper boundary condition to obtain second-order accuracy in time
            p_nh%diag%mflx_ic_int(jc,jb,iadv_rcf+2)    = (SUM(p_nh%diag%mflx_ic_int(jc,jb,iadv_rcf-1:iadv_rcf)) -    &
              SUM(p_nh%diag%mflx_ic_int(jc,jb,1:2)))*rdt_ubc
            ! Shift time level of averaged field back to the beginning of the first dynamic substep
            p_nh%diag%mflx_ic_int(jc,jb,iadv_rcf+1)    = p_nh%diag%mflx_ic_int(jc,jb,iadv_rcf+1) -       &
              dthalf*p_nh%diag%mflx_ic_int(jc,jb,iadv_rcf+2)
          ENDDO
        ENDIF
      ENDIF

    ENDDO
!$OMP END DO


    IF (ltransport .AND. lstep_adv) THEN

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int-2)

        DO jt = 1,ntracer
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              p_nh%diag%grf_tend_tracer(jc,jk,jb,jt) =                    &
                &            ( p_prog_new_rcf%tracer(jc,jk,jb,jt)               &
                &            -  p_prog_now_rcf%tracer(jc,jk,jb,jt) )*rdt_rcf
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

    ENDIF


    ! edge-based variables
    i_startblk = p_patch(jg)%edges%start_blk(grf_bdywidth_e+1,1)
    i_endblk   = p_patch(jg)%edges%end_blk(min_rledge_int-3,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch(jg), jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, grf_bdywidth_e+1, min_rledge_int-3)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          p_nh%diag%grf_tend_vn(je,jk,jb) = &
            ( p_prog_new%vn(je,jk,jb) - p_prog_now%vn(je,jk,jb) )*rdt
          p_nh%diag%grf_tend_mflx(je,jk,jb) = &
            ( p_nh%diag%mass_fl_e(je,jk,jb) - p_nh%diag%mass_fl_e_sv(je,jk,jb) )*rdt_mflx
        ENDDO
      ENDDO

      IF (l_child_vertnest) THEN ! Compute tendencies for upper boundary condition
        DO je = i_startidx, i_endidx
          p_nh%diag%dvn_ie_int(je,jb) = 0.5_wp*(p_nh%diag%dvn_ie_int(je,jb) + &
            p_nh%diag%vn_ie(je,nshift,jb) - p_nh%diag%vn_ie(je,nshift+1,jb))
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_tendencies

  !-------------------------------------------------------------------------
  !
  !>
  !! Interpolates time tendencies of prognostic variables to the lateral boundary
  !! of a refined mesh
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2008-07-10
  !!
  SUBROUTINE boundary_interpolation (jg,jgc,ntp_dyn,ntc_dyn,ntp_tr,ntc_tr,lstep_adv, &
    mass_flx_p,mass_flx_c)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_nest_utilities:boundary_interpolation'


    INTEGER, INTENT(IN)     :: jg, jgc      ! domain ID of parent and child grid

    ! Parent and child time levels for dynamical variables and tracers
    INTEGER, INTENT(IN)     :: ntp_dyn, ntc_dyn, ntp_tr, ntc_tr

    LOGICAL, INTENT(IN) :: lstep_adv  ! time tendencies are only interpolated if transport is called

    ! Mass fluxes at parent and child level
    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) :: mass_flx_p, mass_flx_c

    ! local variables

    TYPE(t_nh_diag), POINTER           :: p_diagp => NULL()
    TYPE(t_nh_diag), POINTER           :: p_diagc  => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhp_dyn => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhc_dyn => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhp_tr  => NULL()
    TYPE(t_nh_prog), POINTER           :: p_nhc_tr  => NULL()
    TYPE(t_patch), POINTER             :: p_pp => NULL()
    TYPE(t_patch), POINTER             :: p_pc => NULL()
    TYPE(t_int_state), POINTER         :: p_int => NULL()
    TYPE(t_gridref_state), POINTER     :: p_grf => NULL()
    TYPE(t_gridref_state), POINTER     :: p_grfc => NULL()
    TYPE(t_grid_cells), POINTER        :: p_gcp => NULL()
    TYPE(t_grid_edges), POINTER        :: p_gep => NULL()


    INTEGER :: i_startblk              ! start block
    INTEGER :: i_endblk                ! end index
    INTEGER :: i_startidx              ! start index
    INTEGER :: i_endidx                ! end index

    INTEGER :: jb, jc, jk, jt, ic      ! loop indices

    INTEGER :: nlev_c, nlevp1_c        ! number of full and half levels (child domain)

    INTEGER :: i_chidx, i_nchdom, i_sbc, i_ebc

    REAL(wp) :: aux3dp(nproma,ntracer+4,p_patch(jg)%nblks_c), &
      aux3dc(nproma,ntracer+4,p_patch(jgc)%nblks_c), &
      theta_prc(nproma,p_patch(jgc)%nlev,p_patch(jgc)%nblks_c), &
      rho_prc(nproma,p_patch(jgc)%nlev,p_patch(jgc)%nblks_c)

    ! Pointers to index fields
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

    ! Switch to determine manner of OpenMP parallelization in interpol_scal_grf
    LOGICAL :: lpar_fields=.FALSE.

    ! Switch to control if the child domain is vertically nested and therefore 
    ! needs interpolation of upper boundary conditions
    LOGICAL :: l_child_vertnest

    LOGICAL :: l_limit(2*ntracer)

    !$ INTEGER :: num_threads_omp, omp_get_max_threads
    !-----------------------------------------------------------------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') '========= Interpolate:',jg,' =>',jgc
      CALL message(TRIM(routine),message_text)
    ENDIF

    !$  num_threads_omp = omp_get_max_threads()

    p_diagp       => p_nh_state(jg)%diag
    p_diagc       => p_nh_state(jgc)%diag
    p_nhp_dyn     => p_nh_state(jg)%prog(ntp_dyn)
    p_nhc_dyn     => p_nh_state(jgc)%prog(ntc_dyn)
    p_nhp_tr      => p_nh_state(jg)%prog(ntp_tr)
    p_nhc_tr      => p_nh_state(jgc)%prog(ntc_tr)
    p_int         => p_int_state(jg)
    p_grf         => p_grf_state(jg)
    p_grfc        => p_grf_state(jgc)
    p_pp          => p_patch(jg)
    p_pc          => p_patch(jgc)
    p_gcp         => p_patch(jg)%cells
    p_gep         => p_patch(jg)%edges

    iidx          => p_gcp%child_idx
    iblk          => p_gcp%child_blk

    i_chidx = p_patch(jgc)%parent_child_index
    i_nchdom   = MAX(1,p_pc%n_childdom)

    ! number of vertical levels (child domain)
    nlev_c   = p_pc%nlev
    nlevp1_c = p_pc%nlevp1

    ! determine if upper boundary interpolation is needed
    IF (lvert_nest .AND. (p_pp%nshift_child > 0)) THEN  
      l_child_vertnest = .TRUE.
    ELSE
      l_child_vertnest = .FALSE.
    ENDIF

    ! Interpolation of cell-based variables


    IF (grf_intmethod_c == 1) THEN ! tendency copying for all cell-based variables

      CALL exchange_data(p_pc%comm_pat_interpolation_c, &
        RECV=p_diagc%grf_tend_rho,     &
        SEND=p_diagp%grf_tend_rho)

      CALL exchange_data(p_pc%comm_pat_interpolation_c, &
        RECV=p_diagc%grf_tend_thv,     &
        SEND=p_diagp%grf_tend_thv)

      ! exchange_data should also work for w because it determines the
      ! vertical dimension with UBOUND
      CALL exchange_data(p_pc%comm_pat_interpolation_c, &
        RECV=p_diagc%grf_tend_w,       &
        SEND=p_diagp%grf_tend_w)

      ! grf_intmethod_c = 2, use gradient at cell center for interpolation
    ELSE IF (grf_intmethod_c == 2) THEN

      !$  IF (num_threads_omp <= 14) lpar_fields=.TRUE.

      ! Interpolation of temporal tendencies, full w, perturbation density (stored in div) 
      !  and perturbationvirtual potential temperature (stored in dpres_mc)
      CALL interpol_scal_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 6,         &
        p_diagp%grf_tend_rho, p_diagc%grf_tend_rho,  &
        p_diagp%grf_tend_thv, p_diagc%grf_tend_thv,  &
        p_diagp%grf_tend_w,   p_diagc%grf_tend_w,    &
        p_nhp_dyn%w,          p_nhc_dyn%w,           &
        p_nh_state(jg)%diag%div, rho_prc,            &
        p_nh_state(jg)%diag%dpres_mc, theta_prc,     &
        lpar_fields=lpar_fields )

      ! Start and end blocks for which interpolation is needed
      i_startblk = p_pc%cells%start_blk(1,1)
      i_endblk   = p_pc%cells%end_blk(grf_bdywidth_c,i_nchdom)

      ! This loop is not OpenMP parallelized because the overhead for opening a
      ! parallel section is too large
      DO jb =  i_startblk, i_endblk

        CALL get_indices_c(p_pc, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
          1, grf_bdywidth_c)

        DO jk = 1, nlev_c
          DO jc = i_startidx, i_endidx
            p_nhc_dyn%rho(jc,jk,jb) = rho_prc(jc,jk,jb) + &
              p_nh_state(jgc)%metrics%rho_ref_mc(jc,jk,jb)
            p_nhc_dyn%theta_v(jc,jk,jb) = theta_prc(jc,jk,jb) + &
              p_nh_state(jgc)%metrics%theta_ref_mc(jc,jk,jb)
          ENDDO
        ENDDO
      ENDDO

      ! The following index list contains the halo points of the lateral boundary
      ! cells. These have to be copied as well in order for rho and theta to be
      ! synchronized.
      DO ic = 1, p_nh_state(jgc)%metrics%bdy_halo_c_dim

        jb = p_nh_state(jgc)%metrics%bdy_halo_c_blk(ic)
        jc = p_nh_state(jgc)%metrics%bdy_halo_c_idx(ic)

        DO jk = 1, nlev_c
          p_nhc_dyn%rho(jc,jk,jb) = rho_prc(jc,jk,jb) + &
            p_nh_state(jgc)%metrics%rho_ref_mc(jc,jk,jb)
          p_nhc_dyn%theta_v(jc,jk,jb) = theta_prc(jc,jk,jb) + &
            p_nh_state(jgc)%metrics%theta_ref_mc(jc,jk,jb)

        ENDDO
      ENDDO
    ENDIF

    ! Interpolation of cell based tracer variables
    IF (ltransport .AND. lstep_adv .AND. grf_intmethod_ct == 1) THEN

      ! Start and end blocks for which interpolation is needed
      i_startblk = p_gcp%start_blk(grf_bdyintp_start_c,i_chidx)
      i_endblk   = p_gcp%end_blk(grf_bdyintp_end_c,i_chidx)

      CALL exchange_data_mult(p_pc%comm_pat_interpolation_c, ntracer, ntracer*nlev_c, &
        RECV4D=p_diagc%grf_tend_tracer,SEND4D=p_diagp%grf_tend_tracer)


    ELSE IF (ltransport .AND. lstep_adv .AND. grf_intmethod_ct == 2) THEN

      !$  IF (num_threads_omp <= 4*ntracer+1) lpar_fields=.TRUE.

      ! Apply positive definite limiter on full tracer fields but not on tendencies
      l_limit(1:ntracer) = .FALSE.
      l_limit(ntracer+1:2*ntracer) = .TRUE.

      CALL interpol_scal_grf ( p_pp, p_pc, p_grf%p_dom(i_chidx), 2*ntracer,        &
        f4din1=p_diagp%grf_tend_tracer, f4dout1=p_diagc%grf_tend_tracer,          &
        f4din2=p_nhp_tr%tracer, f4dout2=p_nhc_tr%tracer, lpar_fields=lpar_fields, &
        llimit_nneg=l_limit)

    ENDIF

    ! Interpolation of edge-based variables  (velocity components)
    IF (grf_intmethod_e == 1 .OR. grf_intmethod_e == 2) THEN

      CALL interpol_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), p_diagp%grf_tend_vn, p_diagc%grf_tend_vn)

    ELSE IF (grf_intmethod_e == 3 .OR. grf_intmethod_e == 4) THEN

      CALL interpol2_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 1, p_diagp%grf_tend_vn, p_diagc%grf_tend_vn)

    ELSE IF (grf_intmethod_e == 5 .OR. grf_intmethod_e == 6) THEN

      CALL interpol2_vec_grf (p_pp, p_pc, p_grf%p_dom(i_chidx), 3,        &
        p_diagp%grf_tend_vn, p_diagc%grf_tend_vn, mass_flx_p, mass_flx_c, &
        p_diagp%grf_tend_mflx, p_diagc%grf_tend_mflx )

    ENDIF


    ! Perform interpolations needed for vertical nesting
    IF (l_child_vertnest) THEN

      IF (p_test_run) THEN
        aux3dp = 0._wp
        aux3dc = 0._wp
      ENDIF

      CALL sync_patch_array(SYNC_E,p_pp,p_diagp%dvn_ie_int)

      CALL interpol_vec_ubc (p_pp, p_pc, p_grf%p_dom(i_chidx), &
        p_diagp%dvn_ie_int, p_diagc%dvn_ie_ubc)

      ! Start and end blocks for which interpolation is needed
      i_startblk = p_pp%cells%start_blk(grf_nudgintp_start_c+2,i_chidx)
      i_endblk   = p_pp%cells%end_blk(min_rlcell_int,i_chidx)

      ! For back-copying at child level
      i_nchdom   = MAX(1,p_pc%n_childdom)
      i_sbc      = p_pc%cells%start_blk(grf_nudge_start_c,1)
      i_ebc      = p_pc%cells%end_blk(min_rlcell_int,i_nchdom)

      IF (ltransport .AND. lstep_adv) THEN

        DO jb =  i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            grf_nudgintp_start_c+2, min_rlcell_int, i_chidx)

          DO jc = i_startidx, i_endidx
            aux3dp(jc,1,jb)   = p_diagp%dw_int(jc,jb,iadv_rcf+1)
            aux3dp(jc,2:3,jb) = p_diagp%mflx_ic_int(jc,jb,iadv_rcf+1:iadv_rcf+2)
            aux3dp(jc,4,jb)   = p_diagp%dtheta_v_ic_int(jc,jb,iadv_rcf+1)
          ENDDO
          DO jt = 1, ntracer
            DO jc = i_startidx, i_endidx
              aux3dp(jc,4+jt,jb) = p_diagp%q_int(jc,jb,jt)
            ENDDO
          ENDDO
        ENDDO

        CALL sync_patch_array(SYNC_C,p_pp,aux3dp)

        CALL interpol_scal_ubc (p_pp, p_pc, p_grf%p_dom(i_chidx),  &
          ntracer+4, aux3dp, aux3dc)

        DO jb = i_sbc, i_ebc

          CALL get_indices_c(p_pc, jb, i_sbc, i_ebc, i_startidx, i_endidx, &
            grf_nudge_start_c, min_rlcell_int)

          DO jc = i_startidx, i_endidx
            p_diagc%dw_ubc(jc,jb)          = aux3dc(jc,1,jb)
            p_diagc%mflx_ic_ubc(jc,jb,1:2) = aux3dc(jc,2:3,jb)
            p_diagc%dtheta_v_ic_ubc(jc,jb) = aux3dc(jc,4,jb)
          ENDDO
          DO jt = 1, ntracer
            DO jc = i_startidx, i_endidx
              p_diagc%q_ubc(jc,jb,jt) = aux3dc(jc,jt+4,jb)
            ENDDO
          ENDDO
        ENDDO

      ELSE

        DO jb =  i_startblk, i_endblk

          CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            grf_nudgintp_start_c+2, min_rlcell_int, i_chidx)

          DO jc = i_startidx, i_endidx
            aux3dp(jc,1,jb)   = p_diagp%dw_int(jc,jb,iadv_rcf+1)
            aux3dp(jc,2:3,jb) = p_diagp%mflx_ic_int(jc,jb,iadv_rcf+1:iadv_rcf+2)
            aux3dp(jc,4,jb)   = p_diagp%dtheta_v_ic_int(jc,jb,iadv_rcf+1)
          ENDDO
        ENDDO

        CALL sync_patch_array(SYNC_C,p_pp,aux3dp)

        CALL interpol_scal_ubc(p_pp, p_pc, p_grf%p_dom(i_chidx), 4, aux3dp, aux3dc)

        DO jb = i_sbc, i_ebc

          CALL get_indices_c(p_pc, jb, i_sbc, i_ebc, i_startidx, i_endidx, &
            grf_nudge_start_c, min_rlcell_int)

          DO jc = i_startidx, i_endidx
            p_diagc%dw_ubc(jc,jb)          = aux3dc(jc,1,jb)
            p_diagc%mflx_ic_ubc(jc,jb,1:2) = aux3dc(jc,2:3,jb)
            p_diagc%dtheta_v_ic_ubc(jc,jb) = aux3dc(jc,4,jb)
          ENDDO
        ENDDO

      ENDIF

    ENDIF

  END SUBROUTINE boundary_interpolation


  !>
  !! This routine prepares boundary nudging for use with 1-way nesting.
  !!
  !! The following steps are executed:
  !! 1. Mapping of coarse-grid prognostic variables to intermediate grid sharing
  !!    the domain decomposition and vertical dimension with the child grid
  !! 2. Computation of differences between parent-grid values and averaged child-grid
  !!    variables
  !! 3. Interpolation of difference fields to the child grid
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2010-06-18
  SUBROUTINE prep_bdy_nudging(jgp, jg, lstep_adv)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_nest_utilities:prep_bdy_nudging'


    INTEGER, INTENT(IN) :: jg   ! child grid level
    INTEGER, INTENT(IN) :: jgp  ! parent grid level


    LOGICAL, INTENT(IN) :: lstep_adv  ! Switch if nudging is done for tracers
    ! (only if tracer transport is called)

    ! local variables

    TYPE(t_nh_prog),    POINTER     :: p_parent_prog => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_child_prog  => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_parent_prog_rcf => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_child_prog_rcf  => NULL()
    TYPE(t_nh_diag),    POINTER     :: p_diag        => NULL()
    TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
    TYPE(t_grid_edges), POINTER     :: p_gep => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grfc => NULL()
    TYPE(t_int_state), POINTER      :: p_int => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()

    ! Indices
    INTEGER :: jb, jc, jk, jt, je, js, i_nchdom, i_chidx, i_startblk, i_endblk, &
      i_startidx, i_endidx, istartblk_c, istartblk_e
    INTEGER :: nlev_c, nlev_p
    INTEGER :: nshift      !< difference between upper boundary of parent or feedback-parent 
    !< domain and upper boundary of child domain (in terms 
    !< of vertical levels) 

    ! Local arrays for interpolated parent-grid values, and difference fields. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_thv, diff_thv
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_rho, diff_rho
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_vn , diff_vn
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_w  , diff_w
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: parent_tr , diff_tr

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk, ieidx, ieblk
    LOGICAL :: l_parallel
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt, p_fbkwgt_tr, p_fbkwgt_v
    !-----------------------------------------------------------------------

    IF (msg_level >= 10) THEN
      WRITE(message_text,'(a,i2,a,i2)') '1-way nesting: == Boundary nudging:',jg
      CALL message(TRIM(routine),message_text)
    ENDIF

    IF (my_process_is_mpi_seq()) THEN
      l_parallel = .FALSE.
    ELSE
      l_parallel = .TRUE.
    ENDIF

    p_parent_prog     => p_nh_state(jgp)%prog(nsav1(jgp))
    p_child_prog      => p_nh_state(jg)%prog(nnow(jg))
    p_parent_prog_rcf => p_nh_state(jgp)%prog(nnow_rcf(jgp))
    p_child_prog_rcf  => p_nh_state(jg)%prog(nnow_rcf(jg))
    p_diag            => p_nh_state(jg)%diag
    p_grfc            => p_grf_state(jg)
    p_pc              => p_patch(jg)

    p_grf => p_grf_state_local_parent(jg)
    p_int => p_int_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_gep => p_patch_local_parent(jg)%edges
    p_pp  => p_patch_local_parent(jg)

    i_nchdom = MAX(1,p_pc%n_childdom)
    i_chidx  = p_pc%parent_child_index

    ! number of full levels of child domain
    nlev_c   = p_pc%nlev

    ! number of full/half levels of parent domain
    nlev_p   = p_pp%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift
    js     = nshift

    ! Please note: In the parallel case
    ! - lower bound must be 1 due to synchronization calls
    ! - upper bound must be nblks_c/e to include halo cells/edges
    ! - this doesn't cost extra memory since p_patch_local_parent
    !   only includes the cells/edges really needed

    ! Value of i_startblk needed for subroutine call for parent-to-child interpolation
    istartblk_c = 1

    ALLOCATE(parent_thv  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c),  &
      diff_thv    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c),  &
      parent_rho  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c),  &
      diff_rho    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c),  &
      parent_w    (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c),  &
      diff_w      (nproma, nlev_c+1, p_patch_local_parent(jg)%nblks_c) )

    IF(ltransport .AND. lstep_adv) &
      ALLOCATE(parent_tr(nproma, nlev_p, p_patch_local_parent(jg)%nblks_c, ntracer),&
      &        diff_tr  (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c, ntracer) )

    ! Value of i_startblk needed for subroutine call for parent-to-child interpolation
    istartblk_e = 1

    ALLOCATE(parent_vn  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_e), &
      &      diff_vn    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_e)  )

    ! Set pointers to index and coefficient fields for cell-based variables
    iidx  => p_gcp%child_idx
    iblk  => p_gcp%child_blk
    ieidx => p_gep%child_idx
    ieblk => p_gep%child_blk

    p_fbkwgt    => p_grf%fbk_wgt_c
    p_fbkwgt_tr => p_grf%fbk_wgt_ct
    p_fbkwgt_v  => p_grf%fbk_wgt_e

    ! 1st step: Copy prognostic variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 3, 3*nlev_p, &
      RECV1=parent_rho, SEND1=p_parent_prog%rho,         &
      RECV2=parent_thv, SEND2=p_parent_prog%theta_v,     &
      RECV3=parent_w,   SEND3=p_parent_prog%w            )

    CALL exchange_data(p_pp%comm_pat_glb_to_loc_e,     &
      RECV=parent_vn, SEND=p_parent_prog%vn)

    IF (ltransport .AND. lstep_adv) &
      CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, ntracer, ntracer*nlev_p, &
      RECV4D=parent_tr, SEND4D=p_parent_prog_rcf%tracer)

    ! 2nd step: perform feedback from refined grid to intermediate grid and compute differences

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! a) cell-based variables
    ! Start/End block in the parent domain
    i_startblk = p_gcp%start_blk(grf_nudgintp_start_c+1,i_chidx)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
        grf_nudgintp_start_c+1, min_rlcell_int, i_chidx)

      ! initialize diff_w at surface with zero
      diff_w(:,nlev_c+1,jb) = 0._wp

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
!CDIR UNROLL=5
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
#endif
          
          diff_thv(jc,jk,jb) = parent_thv(jc,jk+js,jb) - (                           &
            p_child_prog%theta_v(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%theta_v(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%theta_v(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%theta_v(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)   )

          diff_rho(jc,jk,jb) = parent_rho(jc,jk+js,jb) - (                       &
            p_child_prog%rho(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%rho(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%rho(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%rho(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4))+ &
            p_nh_state(jg)%metrics%rho_ref_corr(jc,jk,jb)

          diff_w(jc,jk,jb) = parent_w(jc,jk+js,jb) - (                         &
            p_child_prog%w(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%w(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%w(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%w(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)   )

        ENDDO
      ENDDO

      ! Tracers
      IF (ltransport .AND. lstep_adv) THEN

        DO jt = 1, ntracer

#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = 1, nlev_c
#else
!CDIR UNROLL=8
          DO jk = 1, nlev_c
            DO jc = i_startidx, i_endidx
#endif

              diff_tr(jc,jk,jb,jt) = parent_tr(jc,jk+js,jb,jt) - (                                &
                p_child_prog_rcf%tracer(iidx(jc,jb,1),jk,iblk(jc,jb,1),jt)*p_fbkwgt_tr(jc,jb,1) + &
                p_child_prog_rcf%tracer(iidx(jc,jb,2),jk,iblk(jc,jb,2),jt)*p_fbkwgt_tr(jc,jb,2) + &
                p_child_prog_rcf%tracer(iidx(jc,jb,3),jk,iblk(jc,jb,3),jt)*p_fbkwgt_tr(jc,jb,3) + &
                p_child_prog_rcf%tracer(iidx(jc,jb,4),jk,iblk(jc,jb,4),jt)*p_fbkwgt_tr(jc,jb,4)   )
            ENDDO
          ENDDO

        ENDDO

      ENDIF

    ENDDO
!$OMP END DO

    ! b) velocity
    ! Start/End block in the parent domain
    i_startblk = p_gep%start_blk(grf_nudgintp_start_e+2,i_chidx)
    i_endblk   = p_gep%end_blk(min_rledge_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
        grf_nudgintp_start_e+2, min_rledge_int, i_chidx)


#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
!CDIR UNROLL=5
      DO jk = 1, nlev_c
        DO je = i_startidx, i_endidx
#endif
          
          diff_vn(je,jk,jb) = parent_vn(je,jk+js,jb) - (                            &
            p_child_prog%vn(ieidx(je,jb,1),jk,ieblk(je,jb,1))*p_fbkwgt_v(je,jb,1) + &
            p_child_prog%vn(ieidx(je,jb,2),jk,ieblk(je,jb,2))*p_fbkwgt_v(je,jb,2)   )

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! Interpolate differences to child grid; the differences are stored in the grf_tend fields

    ! interpol_vec/scal_nudging needs all fields with full boundaries, so the arrays
    ! have to be sync'd before calling these routines.
    ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
    ! since the arrays don't start with lower bound 1 in the non parallel case!

    ! Synchronization is needed after the interpolation step for cell-based variables because for
    ! those, the nudging tendencies are applied outside the dynamical core for reasons of mass consistency

    IF(l_parallel) CALL exchange_data(p_pp%comm_pat_e, diff_vn)
    CALL interpol_vec_nudging (p_pp, p_pc, p_int, p_grf%p_dom(i_chidx), p_grfc,    &
      &                        i_chidx, 0, istartblk_e, diff_vn,p_diag%grf_tend_vn )

    IF(l_parallel) CALL exchange_data_mult(p_pp%comm_pat_c, 3, 3*nlev_c+1, &
      recv1=diff_thv, recv2=diff_rho, recv3=diff_w)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, 3, istartblk_c, &
      &                         f3din1=diff_thv, f3dout1=p_diag%grf_tend_thv,                  &
      &                         f3din2=diff_rho, f3dout2=p_diag%grf_tend_rho,                  &
      &                         f3din3=diff_w,   f3dout3=p_diag%grf_tend_w                     )
    CALL sync_patch_array_mult(SYNC_C,p_pc,3,p_diag%grf_tend_thv,p_diag%grf_tend_rho,  &
      p_diag%grf_tend_w)

    IF (ltransport .AND. lstep_adv) THEN
      IF(l_parallel) CALL exchange_data_mult(p_pp%comm_pat_c, ntracer, ntracer*nlev_c, recv4d=diff_tr)

      CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx,  &
        &                         0, ntracer, istartblk_c, f4din=diff_tr,      &
        &                         f4dout=p_diag%grf_tend_tracer )
      CALL sync_patch_array_mult(SYNC_C,p_pc,ntracer,f4din=p_diag%grf_tend_tracer)
    ENDIF

    DEALLOCATE(parent_thv, diff_thv, parent_rho, diff_rho, parent_w, diff_w, parent_vn, diff_vn)

    IF(ltransport .AND. lstep_adv) DEALLOCATE(parent_tr, diff_tr)


  END SUBROUTINE prep_bdy_nudging

  !>
  !! This routine prepares boundary nudging for density only (for use with 2-way nesting)
  !!
  !! The following steps are executed:
  !! 1. Mapping of coarse-grid density to intermediate grid sharing
  !!    the domain decomposition and vertical dimension with the child grid
  !! 2. Computation of differences between parent-grid values and averaged child-grid
  !!    values
  !! 3. Interpolation of difference fields to the child grid
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2011-12-08
  SUBROUTINE prep_rho_bdy_nudging(jgp, jg)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_nest_utilities:prep_rho_bdy_nudging'


    INTEGER, INTENT(IN) :: jg   ! child grid level
    INTEGER, INTENT(IN) :: jgp  ! parent grid level


    ! local variables

    TYPE(t_nh_prog),    POINTER     :: p_parent_prog => NULL()
    TYPE(t_nh_prog),    POINTER     :: p_child_prog  => NULL()
    TYPE(t_nh_diag),    POINTER     :: p_diag        => NULL()
    TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
    TYPE(t_gridref_state), POINTER  :: p_grfc => NULL()
    TYPE(t_int_state), POINTER      :: p_int => NULL()
    TYPE(t_patch),      POINTER     :: p_pp => NULL()
    TYPE(t_patch),      POINTER     :: p_pc => NULL()

    ! Indices
    INTEGER :: jb, jc, jk, js, i_nchdom, i_chidx, i_startblk, i_endblk, &
      i_startidx, i_endidx, istartblk_c
    INTEGER :: nlev_c, nlev_p
    INTEGER :: nshift      !< difference between upper boundary of parent or feedback-parent 
    !< domain and upper boundary of child domain (in terms 
    !< of vertical levels) 

    ! Local arrays for interpolated parent-grid values, and difference fields. These have
    ! to be allocatable because their dimensions differ between MPI and non-MPI runs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: parent_rho, diff_rho

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
    LOGICAL :: l_parallel
    REAL(wp), DIMENSION(:,:,:), POINTER :: p_fbkwgt
    !-----------------------------------------------------------------------


    IF (my_process_is_mpi_seq()) THEN
      l_parallel = .FALSE.
    ELSE
      l_parallel = .TRUE.
    ENDIF

    p_parent_prog     => p_nh_state(jgp)%prog(nsav1(jgp))
    p_child_prog      => p_nh_state(jg)%prog(nnow(jg))
    p_diag            => p_nh_state(jg)%diag
    p_grfc            => p_grf_state(jg)
    p_pc              => p_patch(jg)

    p_grf => p_grf_state_local_parent(jg)
    p_int => p_int_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)

    i_nchdom = MAX(1,p_pc%n_childdom)
    i_chidx  = p_pc%parent_child_index

    ! number of full levels of child domain
    nlev_c   = p_pc%nlev

    ! number of full/half levels of parent domain
    nlev_p   = p_pp%nlev

    ! shift between upper model boundaries
    nshift = p_pc%nshift
    js     = nshift

    ! Please note: In the parallel case
    ! - lower bound must be 1 due to synchronization calls
    ! - upper bound must be nblks_c/e to include halo cells/edges
    ! - this doesn't cost extra memory since p_patch_local_parent
    !   only includes the cells/edges really needed

    ! Value of i_startblk needed for subroutine call for parent-to-child interpolation
    istartblk_c = 1

    ALLOCATE(parent_rho  (nproma, nlev_p, p_patch_local_parent(jg)%nblks_c), &
      &      diff_rho    (nproma, nlev_c, p_patch_local_parent(jg)%nblks_c)  )


    ! Set pointers to index and coefficient fields for cell-based variables
    iidx  => p_gcp%child_idx
    iblk  => p_gcp%child_blk

    p_fbkwgt    => p_grf%fbk_wgt_c

    ! 1st step: Copy prognostic variables from parent grid to fields on feedback-parent grid
    ! (trivial without MPI parallelization, but communication call needed for MPI)

    CALL exchange_data(p_pp%comm_pat_glb_to_loc_c,        &
      &                RECV=parent_rho, SEND=p_parent_prog%rho )  

    ! 2nd step: perform feedback from refined grid to intermediate grid and compute differences

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! Start/End block in the parent domain
    i_startblk = p_gcp%start_blk(grf_nudgintp_start_c+1,i_chidx)
    i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
        grf_nudgintp_start_c+1, min_rlcell_int, i_chidx)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
!CDIR UNROLL=8
      DO jk = 1, nlev_c
        DO jc = i_startidx, i_endidx
#endif
          
          diff_rho(jc,jk,jb) = parent_rho(jc,jk+js,jb) - (                       &
            p_child_prog%rho(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
            p_child_prog%rho(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
            p_child_prog%rho(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
            p_child_prog%rho(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4))+ &
            p_nh_state(jg)%metrics%rho_ref_corr(jc,jk,jb)

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL


    ! Interpolate differences to child grid; the differences are stored in the grf_tend fields

    ! interpol_scal_nudging needs all fields with full boundaries, so the arrays
    ! have to be sync'd before calling these routines.
    ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
    ! since the arrays don't start with lower bound 1 in the non parallel case!

    ! Synchronization is needed after the interpolation step for cell-based variables because for
    ! those, the nudging tendencies are applied outside the dynamical core for reasons of mass consistency

    IF(l_parallel) CALL exchange_data(p_pp%comm_pat_c, diff_rho)
    CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, 1, istartblk_c, &
      &                         f3din1=diff_rho, f3dout1=p_diag%grf_tend_rho                   )
    CALL sync_patch_array(SYNC_C,p_pc,p_diag%grf_tend_rho)


    DEALLOCATE(parent_rho, diff_rho)


  END SUBROUTINE prep_rho_bdy_nudging


  !>
  !! This routine prepares outer boundary nudging for the limited-area mode.
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2013-21-10
  SUBROUTINE prep_outer_bdy_nudging (p_patch, p_prog_latbc, p_prog_now, p_metrics, p_diag)

    TYPE(t_patch),   INTENT(IN)    :: p_patch
    TYPE(t_nh_prog), INTENT(IN)    :: p_prog_latbc, p_prog_now
    TYPE(t_nh_metrics), INTENT(IN) :: p_metrics
    TYPE(t_nh_diag), INTENT(INOUT) :: p_diag

    ! local variables
    INTEGER :: jb, jc, jk, je, ic, nlev

    ! number of full levels of child domain
    nlev   = p_patch%nlev

    ! compute differences between lateral boundary data and prognostic variables

!$OMP PARALLEL
!$OMP DO PRIVATE(jk,jc,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_metrics%nudge_c_dim
      jc = p_metrics%nudge_c_idx(ic)
      jb = p_metrics%nudge_c_blk(ic)
      DO jk = 1, nlev
#else
    DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, p_metrics%nudge_c_dim
        jc = p_metrics%nudge_c_idx(ic)
        jb = p_metrics%nudge_c_blk(ic)
#endif
          p_diag%grf_tend_thv(jc,jk,jb) = p_prog_latbc%theta_v(jc,jk,jb) - p_prog_now%theta_v(jc,jk,jb)
          p_diag%grf_tend_rho(jc,jk,jb) = p_prog_latbc%rho    (jc,jk,jb) - p_prog_now%rho    (jc,jk,jb)
          p_diag%grf_tend_w  (jc,jk,jb) = p_prog_latbc%w      (jc,jk,jb) - p_prog_now%w      (jc,jk,jb)

      ENDDO
    ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jk,je,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
      DO ic = 1, p_metrics%nudge_e_dim
        je = p_metrics%nudge_e_idx(ic)
        jb = p_metrics%nudge_e_blk(ic)
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, p_metrics%nudge_e_dim
          je = p_metrics%nudge_e_idx(ic)
          jb = p_metrics%nudge_e_blk(ic)
#endif
          p_diag%grf_tend_vn(je,jk,jb) = p_prog_latbc%vn(je,jk,jb) - p_prog_now%vn(je,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE prep_outer_bdy_nudging



  !>
  !! This routine executes boundary nudging for limited-area simulations.
  !!
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2010-06-18
  SUBROUTINE outer_boundary_nudging(jg, nnow, rcffac)

    INTEGER, INTENT(IN)  :: jg, nnow
    REAL(wp), INTENT(IN) :: rcffac ! Ratio between advective and dynamical time step

    ! Pointers
    TYPE(t_patch),     POINTER ::  p_p
    TYPE(t_nh_state),  POINTER ::  p_nh
    TYPE(t_int_state), POINTER ::  p_int

    ! Indices
    INTEGER :: jb, jc, jk, ic

    INTEGER :: nlev                   ! number of vertical full levels

    REAL(wp) :: rd_o_cvd, rd_o_p0ref

    ! Set pointers
    p_p   => p_patch(jg)
    p_nh  => p_nh_state(jg)
    p_int => p_int_state(jg)

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    ! number of vertical levels
    nlev = p_p%nlev

!$OMP PARALLEL

!$OMP DO PRIVATE(jk,jc,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_nh%metrics%nudge_c_dim
      jc = p_nh%metrics%nudge_c_idx(ic)
      jb = p_nh%metrics%nudge_c_blk(ic)
      DO jk = 1, nlev
#else
    DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, p_nh%metrics%nudge_c_dim
        jc = p_nh%metrics%nudge_c_idx(ic)
        jb = p_nh%metrics%nudge_c_blk(ic)
#endif
        p_nh%prog(nnow)%rho(jc,jk,jb) =                                      &
          p_nh%prog(nnow)%rho(jc,jk,jb) + rcffac*p_int%nudgecoeff_c(jc,jb)*  &
          p_nh%diag%grf_tend_rho(jc,jk,jb)

        p_nh%prog(nnow)%theta_v(jc,jk,jb) =                                     &
          p_nh%prog(nnow)%theta_v(jc,jk,jb) + rcffac*p_int%nudgecoeff_c(jc,jb)* &
          p_nh%diag%grf_tend_thv(jc,jk,jb)

        p_nh%prog(nnow)%rhotheta_v(jc,jk,jb) =                          &
          p_nh%prog(nnow)%rho(jc,jk,jb)*p_nh%prog(nnow)%theta_v(jc,jk,jb)

        p_nh%prog(nnow)%exner(jc,jk,jb) =                                  &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_nh%prog(nnow)%rhotheta_v(jc,jk,jb)))

        p_nh%prog(nnow)%w(jc,jk,jb) =                                      &
          p_nh%prog(nnow)%w(jc,jk,jb) + rcffac*p_int%nudgecoeff_c(jc,jb)*  &
          p_nh%diag%grf_tend_w(jc,jk,jb)

      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE outer_boundary_nudging

  !>
  !! This routine executes boundary nudging for one-way nested domains
  !!
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2010-06-18
  SUBROUTINE nest_boundary_nudging(jg, nnew, nnew_rcf, rcffac)


    INTEGER, INTENT(IN)  :: jg, nnew, nnew_rcf

    REAL(wp), INTENT(IN) :: rcffac ! Ratio between advective and dynamical time step

    ! Pointers
    TYPE(t_nh_state),  POINTER ::  p_nh
    TYPE(t_int_state), POINTER ::  p_int

    ! Indices
    INTEGER :: jb, jc, jk, jt, ic

    INTEGER :: nlev  ! number of vertical full levels

    REAL(wp) :: rd_o_cvd, rd_o_p0ref

    ! Set pointers
    p_nh  => p_nh_state(jg)
    p_int => p_int_state(jg)

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    ! number of vertical levels
    nlev = p_patch(jg)%nlev

!$OMP PARALLEL

!$OMP DO PRIVATE(jk,jc,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_nh%metrics%nudge_c_dim
      jc = p_nh%metrics%nudge_c_idx(ic)
      jb = p_nh%metrics%nudge_c_blk(ic)
      DO jk = 1, nlev
#else
    DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, p_nh%metrics%nudge_c_dim
        jc = p_nh%metrics%nudge_c_idx(ic)
        jb = p_nh%metrics%nudge_c_blk(ic)
#endif
        p_nh%prog(nnew)%rho(jc,jk,jb) =                                      &
          p_nh%prog(nnew)%rho(jc,jk,jb) + rcffac*p_int%nudgecoeff_c(jc,jb)*  &
          p_nh%diag%grf_tend_rho(jc,jk,jb)

        p_nh%prog(nnew)%theta_v(jc,jk,jb) =                                     &
          p_nh%prog(nnew)%theta_v(jc,jk,jb) + rcffac*p_int%nudgecoeff_c(jc,jb)* &
          p_nh%diag%grf_tend_thv(jc,jk,jb)

        p_nh%prog(nnew)%rhotheta_v(jc,jk,jb) =                          &
          p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)

        p_nh%prog(nnew)%exner(jc,jk,jb) =                                  &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_nh%prog(nnew)%rhotheta_v(jc,jk,jb)))

        p_nh%prog(nnew)%w(jc,jk,jb) =                                      &
          p_nh%prog(nnew)%w(jc,jk,jb) + rcffac*p_int%nudgecoeff_c(jc,jb)*  &
          p_nh%diag%grf_tend_w(jc,jk,jb)

      ENDDO
    ENDDO
!$OMP END DO

    IF (ltransport) THEN
!$OMP DO PRIVATE(jk,jc,jb,jt,ic) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
      DO ic = 1, p_nh%metrics%nudge_c_dim
        jc = p_nh%metrics%nudge_c_idx(ic)
        jb = p_nh%metrics%nudge_c_blk(ic)
        DO jt = 1, ntracer
          DO jk = 1, nlev
#else
      DO jt = 1, ntracer
        DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, p_nh%metrics%nudge_c_dim
            jc = p_nh%metrics%nudge_c_idx(ic)
            jb = p_nh%metrics%nudge_c_blk(ic)
#endif

            p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt) = p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt) + &
              rcffac*p_int%nudgecoeff_c(jc,jb)*p_nh%diag%grf_tend_tracer(jc,jk,jb,jt)

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
    ENDIF

!$OMP END PARALLEL


  END SUBROUTINE nest_boundary_nudging

  !>
  !! This routine executes boundary nudging for density (for use with 2-way nesting)
  !!
  !!
  !! @par Revision History
  !! Developed  by Guenther Zaengl, DWD, 2011-12-08
  SUBROUTINE density_boundary_nudging(jg, nnew, rcffac)


    INTEGER, INTENT(IN)  :: jg, nnew

    REAL(wp), INTENT(IN) :: rcffac ! Ratio between advective and dynamical time step

    ! Pointers
    TYPE(t_nh_state),  POINTER ::  p_nh
    TYPE(t_int_state), POINTER ::  p_int

    ! Indices
    INTEGER :: jb, jc, jk, ic

    INTEGER :: nlev  ! number of vertical full levels

    REAL(wp) :: rd_o_cvd, rd_o_p0ref

    ! Set pointers
    p_nh  => p_nh_state(jg)
    p_int => p_int_state(jg)

    ! R/c_v (not present in physical constants)
    rd_o_cvd = 1._wp / cvd_o_rd

    ! R / p0ref
    rd_o_p0ref = rd / p0ref

    ! number of vertical levels
    nlev = p_patch(jg)%nlev

!$OMP PARALLEL

!$OMP DO PRIVATE(jk,jc,jb,ic) ICON_OMP_DEFAULT_SCHEDULE
#ifdef __LOOP_EXCHANGE
    DO ic = 1, p_nh%metrics%nudge_c_dim
      jc = p_nh%metrics%nudge_c_idx(ic)
      jb = p_nh%metrics%nudge_c_blk(ic)
      DO jk = 1, nlev
#else
    DO jk = 1, nlev
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, p_nh%metrics%nudge_c_dim
        jc = p_nh%metrics%nudge_c_idx(ic)
        jb = p_nh%metrics%nudge_c_blk(ic)
#endif
        p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnew)%rho(jc,jk,jb) +  &
          MIN(0.333_wp,3._wp*rcffac*p_int%nudgecoeff_c(jc,jb))*          &
          p_nh%diag%grf_tend_rho(jc,jk,jb)

        p_nh%prog(nnew)%rhotheta_v(jc,jk,jb) =                          &
          p_nh%prog(nnew)%rho(jc,jk,jb)*p_nh%prog(nnew)%theta_v(jc,jk,jb)

        p_nh%prog(nnew)%exner(jc,jk,jb) =                                  &
          EXP(rd_o_cvd*LOG(rd_o_p0ref*p_nh%prog(nnew)%rhotheta_v(jc,jk,jb)))

      ENDDO
    ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  END SUBROUTINE density_boundary_nudging


END MODULE mo_nh_nest_utilities

