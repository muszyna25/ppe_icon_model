!>
!! This module contains utility routines needed for the initialization of the
!! NH model
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-06-29)
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
!----------------------------

MODULE mo_nh_init_utils

  USE mo_kind,                  ONLY: wp
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_metrics, t_nh_state
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_lnd_types,         ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_ext_data_types,        ONLY: t_external_data
  USE mo_parallel_config,       ONLY: nproma
  USE mo_run_config,            ONLY: msg_level, ntracer
  USE mo_grid_config,           ONLY: l_limited_area, n_dom
  USE mo_dynamics_config,       ONLY: iequations, nnow, nnow_rcf
  USE mo_physical_constants,    ONLY: grav, cpd, rd, cvd_o_rd, p0ref
  USE mo_vertical_coord_table,  ONLY: vct_a, vct_b, read_vct
  USE mo_nonhydrostatic_config, ONLY: ivctype
  USE mo_sleve_config,          ONLY: min_lay_thckn, max_lay_thckn, htop_thcknlimit, top_height, &
                                      decay_scale_1, decay_scale_2, decay_exp, flat_height, stretch_fac
  USE mo_impl_constants,        ONLY: max_dom, MAX_CHAR_LENGTH, nclass_aero
  USE mo_math_constants,        ONLY: pi
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_sync,                  ONLY: sync_patch_array, SYNC_C
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp,                  ONLY: edges2cells_scalar
  USE mo_math_laplace,          ONLY: nabla2_scalar
  USE mo_math_gradients,        ONLY: grad_fd_norm
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_initicon_types,        ONLY: t_saveinit_state
  USE mo_initicon_config,       ONLY: type_iau_wgt, is_iau_active, &
    &                                 iau_wgt_dyn, iau_wgt_adv, ltile_coldstart
  USE mo_atm_phy_nwp_config,    ONLY: iprog_aero
  USE mo_lnd_nwp_config,        ONLY: ntiles_total, l2lay_rho_snow, ntiles_water, lmulti_snow, &
                                      nlev_soil, nlev_snow, lsnowtile
!!$  USE mo_util_uuid,             ONLY: t_uuid,  uuid_generate, uuid_parse, &
!!$       &                              uuid_unparse, uuid_string_length
  USE mo_fortran_tools,         ONLY: init, copy

  IMPLICIT NONE

  PRIVATE


  INTEGER:: nflat, nflatlev(max_dom)

  REAL(wp) :: layer_thickness        ! (m)
  INTEGER  :: n_flat_level

  TYPE(t_saveinit_state), ALLOCATABLE  :: saveinit(:)

  PUBLIC :: nflat, nflatlev, n_flat_level, layer_thickness

  PUBLIC :: hydro_adjust, init_hybrid_coord, init_sleve_coord, compute_smooth_topo, &
    &       init_vert_coord, interp_uv_2_vn, init_w, adjust_w, convert_thdvars,     &
    &       convert_omega2w, hydro_adjust_downward, prepare_hybrid_coord,           &
    &       prepare_sleve_coord, save_initial_state, restore_initial_state
  PUBLIC :: compute_iau_wgt

CONTAINS
  !-------------
  !>
  !! SUBROUTINE hydro_adjust
  !! Computes hydrostatically balanced initial condition by bottom-up integration
  !! Virtual temperature is kept constant during the adjustment process
  !!
  !! Input/Output: density, Exner pressure, virtual potential temperature
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-06-29)
  !!
  !!
  !!
  SUBROUTINE hydro_adjust(p_patch, p_nh_metrics, rho, exner, theta_v )


    TYPE(t_patch),      INTENT(IN)       :: p_patch
    TYPE(t_nh_metrics), INTENT(IN)       :: p_nh_metrics

    ! Thermodynamic fields - all defined at full model levels
    REAL(wp), INTENT(INOUT) :: rho(:,:,:)        ! density (kg/m**3)
    REAL(wp), INTENT(INOUT) :: exner(:,:,:)      ! Exner pressure
    REAL(wp), INTENT(INOUT) :: theta_v(:,:,:)    ! virtual potential temperature (K)


    ! LOCAL VARIABLES
    REAL(wp) :: temp_v(nproma,p_patch%nlev) ! virtual temperature
    REAL(wp), DIMENSION(nproma) :: z_fac1, z_fac2, z_fac3, za, zb, zc

    INTEGER :: jb, jk, jc
    INTEGER :: nlen, nlev

    nlev = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,temp_v,z_fac1,z_fac2,z_fac3,za,zb,zc) ICON_OMP_DEFAULT_SCHEDULE

    ! The full model grid including the lateral boundary interpolation zone of
    ! nested domains and MPI-halo points is processed; depending on the setup
    ! of the parallel-read routine, the input fields may need to be synchronized
    ! before entering this routine.

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      ! Compute virtual temperature
      DO jk = 1, nlev
        DO jc = 1, nlen
          temp_v(jc,jk) = theta_v(jc,jk,jb)*exner(jc,jk,jb)
        ENDDO
      ENDDO

      ! Now compute hydrostatically balanced prognostic fields:
      ! The following expressions are derived from the discretized (!) third
      ! equation of motion, assuming dw/dt = 0, and solved for the exner pressure.
      ! Because the vertical discretization differs between the triangular and
      ! hexagonal NH cores, a case discrimination is needed here
      DO jk = nlev-1, 1, -1
        DO jc = 1, nlen
          z_fac1(jc) = p_nh_metrics%wgtfac_c(jc,jk+1,jb)*(temp_v(jc,jk+1) &
            - p_nh_metrics%theta_ref_mc(jc,jk+1,jb)*exner(jc,jk+1,jb))    &
            - (1._wp-p_nh_metrics%wgtfac_c(jc,jk+1,jb))                   &
            * p_nh_metrics%theta_ref_mc(jc,jk,jb)*exner(jc,jk+1,jb)

          z_fac2(jc) = (1._wp-p_nh_metrics%wgtfac_c(jc,jk+1,jb))*temp_v(jc,jk) &
            *exner(jc,jk+1,jb)

          z_fac3(jc) = p_nh_metrics%exner_ref_mc(jc,jk+1,jb)     &
            -p_nh_metrics%exner_ref_mc(jc,jk,jb)-exner(jc,jk+1,jb)

          za(jc) = (p_nh_metrics%theta_ref_ic(jc,jk+1,jb)                     &
            *exner(jc,jk+1,jb)+z_fac1(jc))/p_nh_metrics%ddqz_z_half(jc,jk+1,jb)

          zb(jc) = -(za(jc)*z_fac3(jc)+z_fac2(jc)/p_nh_metrics%ddqz_z_half(jc,jk+1,jb) &
            + z_fac1(jc)*p_nh_metrics%d_exner_dz_ref_ic(jc,jk+1,jb))

          zc(jc) = -(z_fac2(jc)*z_fac3(jc)/p_nh_metrics%ddqz_z_half(jc,jk+1,jb) &
            + z_fac2(jc)*p_nh_metrics%d_exner_dz_ref_ic(jc,jk+1,jb))
        ENDDO !jc

        DO jc = 1, nlen
          exner(jc,jk,jb)      = (zb(jc)+SQRT(zb(jc)**2+4._wp*za(jc)*zc(jc)))/(2._wp*za(jc))
          theta_v(jc,jk,jb)    = temp_v(jc,jk)/exner(jc,jk,jb)
          rho(jc,jk,jb)        = exner(jc,jk,jb)**cvd_o_rd*p0ref/(rd*theta_v(jc,jk,jb))
        ENDDO

      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE hydro_adjust

  !-------------
  !>
  !! SUBROUTINE hydro_adjust_downward
  !! Computes hydrostatically balanced initial condition by top-down integration
  !! In contrast to the above routine, virtual potential temperature is kept constant
  !! during the adjustment, leading to a simpler formula
  !!
  !! Input/Output: density, Exner pressure, virtual potential temperature
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2012-12-28)
  !!
  !!
  !!
  SUBROUTINE hydro_adjust_downward(p_patch, p_nh_metrics, rho, exner, theta_v)


    TYPE(t_patch),      INTENT(IN)       :: p_patch
    TYPE(t_nh_metrics), INTENT(IN)       :: p_nh_metrics

    ! Thermodynamic fields - all defined at full model levels
    REAL(wp), INTENT(INOUT) :: rho(:,:,:)        ! density (kg/m**3)
    REAL(wp), INTENT(INOUT) :: exner(:,:,:)      ! Exner pressure
    REAL(wp), INTENT(INOUT) :: theta_v(:,:,:)    ! virtual potential temperature (K)


    ! LOCAL VARIABLES
    REAL(wp), DIMENSION(nproma) :: theta_v_pr_ic

    INTEGER :: jb, jk, jc
    INTEGER :: nlen, nlev

    nlev = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,theta_v_pr_ic) ICON_OMP_DEFAULT_SCHEDULE

    ! The full model grid including the lateral boundary interpolation zone of
    ! nested domains and MPI-halo points is processed; depending on the setup
    ! of the parallel-read routine, the input fields may need to be synchronized
    ! before entering this routine.

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      ! Now compute hydrostatically balanced prognostic fields:
      ! The following expressions are derived from the discretized (!) third
      ! equation of motion, assuming dw/dt = 0, and solved for the exner pressure.
      DO jk = 2, nlev
        DO jc = 1, nlen
          theta_v_pr_ic(jc) = p_nh_metrics%wgtfac_c(jc,jk,jb) *        &
           (theta_v(jc,jk,jb) - p_nh_metrics%theta_ref_mc(jc,jk,jb)) + &
           (1._wp-p_nh_metrics%wgtfac_c(jc,jk,jb)) *                   &
           (theta_v(jc,jk-1,jb)-p_nh_metrics%theta_ref_mc(jc,jk-1,jb)  )

          exner(jc,jk,jb) = exner(jc,jk-1,jb) + p_nh_metrics%exner_ref_mc(jc,jk,jb) -      &
            p_nh_metrics%exner_ref_mc(jc,jk-1,jb) + p_nh_metrics%ddqz_z_half(jc,jk,jb)*    &
            theta_v_pr_ic(jc)*p_nh_metrics%d_exner_dz_ref_ic(jc,jk,jb)/(theta_v_pr_ic(jc)+ &
            p_nh_metrics%theta_ref_ic(jc,jk,jb))

          rho(jc,jk,jb) = exner(jc,jk,jb)**cvd_o_rd*p0ref/(rd*theta_v(jc,jk,jb))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE hydro_adjust_downward

  !-------------
  !>
  !! SUBROUTINE convert_thdvars
  !! Converts the hydrostatic set of thermodynamic variables into the nonhydrostatic one
  !!
  !! Required input fields: pressure, virtual temperature
  !! Output: density, Exner pressure, virtual potential temperature
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-06-29)
  !!
  !!
  !!
  SUBROUTINE convert_thdvars(p_patch, pres, temp_v, &
                             rho, exner, theta_v    )


    TYPE(t_patch), INTENT(IN) :: p_patch

    ! Input fields - all defined at full model levels
    REAL(wp), INTENT(IN) :: pres  (:,:,:) ! pressure (Pa)
    REAL(wp), INTENT(IN) :: temp_v(:,:,:) ! virtual temperature (K)

    ! Output fields (prognostic model variables) - all defined at full model levels
    REAL(wp), INTENT(OUT) :: rho(:,:,:)        ! density (kg/m**3)
    REAL(wp), INTENT(OUT) :: exner(:,:,:)      ! Exner pressure
    REAL(wp), INTENT(OUT) :: theta_v(:,:,:)    ! virtual potential temperature (K)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, jc
    INTEGER :: nlen, nlev

    nlev = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      DO jk = 1, nlev
        DO jc = 1, nlen
          exner(jc,jk,jb)   = (pres(jc,jk,jb)/p0ref)**(rd/cpd)
          theta_v(jc,jk,jb) = temp_v(jc,jk,jb)/exner(jc,jk,jb)
          rho(jc,jk,jb)     = exner(jc,jk,jb)**cvd_o_rd*p0ref/rd/theta_v(jc,jk,jb)
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE convert_thdvars

  !-------------
  !>
  !! SUBROUTINE convert_omega2w
  !! Converts the hydrostatic vertical velocity (omega, Pa/s)
  !! into physical vertical velocity (m/s)
  !! Note: this routine has to be called on the input grid,
  !! where omega, pressure and temperature are not vertically staggered
  !!
  !! Required input fields: omega, pressure, temperature
  !! Output: vertical wind speed
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-20)
  !!
  !!
  !!
  SUBROUTINE convert_omega2w(omega, w, pres, temp, nblks, npromz, nlev)


    ! Input fields
    REAL(wp), INTENT(IN) :: omega (:,:,:) ! omega (Pa/s)
    REAL(wp), INTENT(IN) :: pres  (:,:,:) ! pressure (Pa)
    REAL(wp), INTENT(IN) :: temp  (:,:,:) ! virtual temperature (K)

    ! Output
    REAL(wp), INTENT(OUT) :: w(:,:,:)  ! vertical velocity (m/s)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlev       ! Number of model levels


    ! LOCAL VARIABLES
    INTEGER :: jb, jk, jc
    INTEGER :: nlen

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
      ENDIF

      DO jk = 1, nlev
        DO jc = 1, nlen
          w(jc,jk,jb) = -rd*omega(jc,jk,jb)*temp(jc,jk,jb)/(grav*pres(jc,jk,jb))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE convert_omega2w


  !-------------
  !>
  !! SUBROUTINE interp_uv_2_vn
  !! Interpolates u and v on cell points to vn on edge points
  !!
  !! Required input fields: u and v on cell points
  !! Output: vn on edge points
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-11)
  !!
  !!
  !!
  SUBROUTINE interp_uv_2_vn(p_patch, p_int, u, v, vn )


    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
    TYPE(t_int_state),     INTENT(IN)   :: p_int

    ! Input fields - all defined at full model levels
    REAL(wp), INTENT(IN) :: u(:,:,:) ! zonal wind component on cell points (m/s)
    REAL(wp), INTENT(IN) :: v(:,:,:) ! meridional wind component on cell points (m/s)

    ! Output field (prognostic model variable) - defined at full model levels
    ! Intent (INOUT) because lateral nest boundaries cannot be filled here
    REAL(wp), INTENT(INOUT) :: vn(:,:,:)  ! edge-normal wind component (m/s)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, je
    INTEGER :: nlev, nblks_e, i_startblk,i_endblk, i_startidx,i_endidx
    REAL(wp) :: z_u, z_v

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

    nlev = p_patch%nlev
    nblks_e = p_patch%nblks_e

    iidx => p_patch%edges%cell_idx
    iblk => p_patch%edges%cell_blk

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    IF (l_limited_area .OR. p_patch%id > 1) THEN ! Fill outermost nest boundary

      i_startblk = p_patch%edges%start_blk(1,1)
      i_endblk   = p_patch%edges%end_blk(1,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk,z_u,z_v) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, 1)

        DO je = i_startidx, i_endidx
          IF (iidx(je,jb,1) >= 1 .AND. iblk(je,jb,1) >= 1) THEN
            DO jk = 1, nlev
              z_u = u(iidx(je,jb,1),jk,iblk(je,jb,1))
              z_v = v(iidx(je,jb,1),jk,iblk(je,jb,1))
              vn(je,jk,jb) =  z_u*p_patch%edges%primal_normal(je,jb)%v1 + &
                              z_v*p_patch%edges%primal_normal(je,jb)%v2
            END DO
          ELSE IF (iidx(je,jb,2) >= 1 .AND. iblk(je,jb,2) >= 1) THEN
            DO jk = 1, nlev
              z_u = u(iidx(je,jb,2),jk,iblk(je,jb,2))
              z_v = v(iidx(je,jb,2),jk,iblk(je,jb,2))
              vn(je,jk,jb) =  z_u*p_patch%edges%primal_normal(je,jb)%v1 + &
                              z_v*p_patch%edges%primal_normal(je,jb)%v2
            END DO
          ENDIF
        END DO

      END DO
!$OMP END DO
    ENDIF

    i_startblk = p_patch%edges%start_blk(2,1)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE

      DO jb = i_startblk, nblks_e

        CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif

            vn(je,jk,jb) = p_int%c_lin_e(je,1,jb)                                              &
              *(u(iidx(je,jb,1),jk,iblk(je,jb,1))*p_patch%edges%primal_normal_cell(je,jb,1)%v1 &
              + v(iidx(je,jb,1),jk,iblk(je,jb,1))*p_patch%edges%primal_normal_cell(je,jb,1)%v2)&
              +            p_int%c_lin_e(je,2,jb)                                              &
              *(u(iidx(je,jb,2),jk,iblk(je,jb,2))*p_patch%edges%primal_normal_cell(je,jb,2)%v1 &
              + v(iidx(je,jb,2),jk,iblk(je,jb,2))*p_patch%edges%primal_normal_cell(je,jb,2)%v2 )

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE interp_uv_2_vn

  !-------------
  !>
  !! SUBROUTINE init_w
  !! Initializes the vertical wind field based on the lower boundary condition
  !! w = v grad h and an empirical vertical decay function
  !! The discretization used here is simpler than that used in the dynamical core
  !! but is sufficient to avoid excessive generation of sound waves during the start phase
  !!
  !! Required input fields: vn, z_ifc
  !! Output: w
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-11)
  !!
  !!
  !!
  SUBROUTINE init_w(p_patch, p_int, vn, z_ifc, w)


    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
    TYPE(t_int_state),     INTENT(IN)   :: p_int

    ! Input fields
    REAL(wp), INTENT(IN) :: vn(:,:,:)    ! edge-normal wind component (m/s)
    REAL(wp), INTENT(IN) :: z_ifc(:,:,:) ! height of half levels (m)

    ! Output field - defined at half model levels
    ! Intent (INOUT) because lateral nest boundaries cannot be filled here
    REAL(wp), INTENT(INOUT) :: w(:,:,:)  ! vertical wind component (m/s)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, je, jc, ktop
    INTEGER :: nlev, nlevp1, nblks_e, nblks_c, nshift, i_startblk, i_startidx, i_endidx

    REAL(wp) :: z_wsfc_e(nproma,1,p_patch%nblks_e) ! w at surface (edge points)
    REAL(wp) :: z_wsfc_c(nproma,1,p_patch%nblks_c) ! w at surface (cell points)
    REAL(wp) :: z_slope_e(nproma,p_patch%nlevp1,p_patch%nblks_e) ! slope at edges

    nlev    = p_patch%nlev
    nlevp1  = p_patch%nlevp1
    nblks_e = p_patch%nblks_e
    nblks_c = p_patch%nblks_c
    nshift  = p_patch%nshift_total

    ! In order to initialize w(1) = 0 except for vertical nesting
    IF (nshift == 0) THEN
      ktop = 2
    ELSE
      ktop = 1
    ENDIF

    ! Compute slope at edges
    CALL grad_fd_norm (z_ifc, p_patch, z_slope_e, 1, nlevp1)

    ! slope cannot be computed at outer boundary edges
    i_startblk = p_patch%edges%start_blk(2,1)

    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

      ! Extrapolation of vn to half levels is neglected here
      DO je = i_startidx, i_endidx
        z_wsfc_e(je,1,jb) = vn(je,nlev,jb)*z_slope_e(je,nlevp1,jb)
      ENDDO
    ENDDO

    CALL edges2cells_scalar(z_wsfc_e,p_patch,p_int%e_inn_c,z_wsfc_c,&
                            1,1,opt_rlstart=2)

    i_startblk = p_patch%cells%start_blk(2,1)

!$OMP PARALLEL
    ! First, initialize w with zero in order to avoid undefined nest boundary points
    CALL init(w(:,:,:))
!$OMP BARRIER

    ! specify a reasonable initial vertical wind speed
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 2)

      DO jc = i_startidx, i_endidx
        w(jc,nlevp1,jb) = z_wsfc_c(jc,1,jb)
      ENDDO
      DO jk = nlev, ktop, -1
        DO jc = i_startidx, i_endidx
          w(jc,jk,jb) = z_wsfc_c(jc,1,jb)*vct_b(jk+nshift)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE init_w


  !-------------
  !>
  !! SUBROUTINE adjust_w
  !! Computes the lower boundary condition for w in a similar way as init_w,
  !! but the result is then merged with an already available vertical wind
  !! field provided by an external data source.
  !!
  !! Required input fields: vn, w, z_ifc
  !! Output: w
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-18)
  !!
  !!
  !!
  SUBROUTINE adjust_w(p_patch, p_int, vn, z_ifc, w)


    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
    TYPE(t_int_state),     INTENT(IN)   :: p_int

    ! Input fields
    REAL(wp), INTENT(IN) :: vn(:,:,:)    ! edge-normal wind component (m/s)
    REAL(wp), INTENT(IN) :: z_ifc(:,:,:) ! height of half levels (m)

    ! INOUT field - defined at half model levels
    REAL(wp), INTENT(INOUT) :: w(:,:,:)  ! vertical wind component (m/s)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, je, jc, ktop
    INTEGER :: nlev, nlevp1, nblks_e, nblks_c, nshift, i_startblk, i_startidx, i_endidx
    REAL(wp):: wfac

    REAL(wp) :: z_wsfc_e(nproma,1,p_patch%nblks_e) ! w at surface (edge points)
    REAL(wp) :: z_wsfc_c(nproma,1,p_patch%nblks_c) ! w at surface (cell points)
    REAL(wp) :: z_slope_e(nproma,p_patch%nlevp1,p_patch%nblks_e) ! slope at edges

    nlev    = p_patch%nlev
    nlevp1  = p_patch%nlevp1
    nblks_e = p_patch%nblks_e
    nblks_c = p_patch%nblks_c
    nshift  = p_patch%nshift_total

    ! In order to initialize w(1) = 0 except for vertical nesting
    IF (nshift == 0) THEN
      ktop = 2
    ELSE
      ktop = 1
    ENDIF

    ! Compute slope at edges
    CALL grad_fd_norm (z_ifc, p_patch, z_slope_e, 1, nlevp1)

    ! slope cannot be computed at outer boundary edges
    i_startblk = p_patch%edges%start_blk(2,1)

    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

      ! Extrapolation of vn to half levels is neglected here
      DO je = i_startidx, i_endidx
        z_wsfc_e(je,1,jb) = vn(je,nlev,jb)*z_slope_e(je,nlevp1,jb)
      ENDDO
    ENDDO

    CALL edges2cells_scalar(z_wsfc_e,p_patch,p_int%e_inn_c,z_wsfc_c,&
                            1,1,opt_rlstart=2)

    i_startblk = p_patch%cells%start_blk(2,1)

    ! specify lower boundary condition and merge with w field provided on input
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,wfac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 2)

      ! Lower boundary condition
      DO jc = i_startidx, i_endidx
        w(jc,nlevp1,jb) = z_wsfc_c(jc,1,jb)
      ENDDO

      ! Merging of lower boundary condition with interpolated data
      DO jk = nlev, ktop, -1
        wfac = vct_b(jk+nshift)**2
        DO jc = i_startidx, i_endidx
          w(jc,jk,jb) = (1._wp-wfac)*w(jc,jk,jb) + wfac*z_wsfc_c(jc,1,jb)
        ENDDO
      ENDDO

      IF (nshift == 0) THEN
        ! Upper boundary condition and smooth transition below
        ! if domain is not vertically nested
        DO jc = i_startidx, i_endidx
          w(jc,1,jb) = 0._wp
          w(jc,2,jb) = 0.33_wp*w(jc,2,jb)
          w(jc,3,jb) = 0.66_wp*w(jc,3,jb)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE adjust_w


  !---------------------------------------------------------------------------
  !>
  !! Initialize hybrid coords by reading the 'a' and 'b'. They are assumed to
  !! be HEIGHT BASED in contrast to the hydrostatic model version. The file
  !! name which contains those data has the same name and structure as its
  !! hydrostatic counterpart, namely 'HYB_PARAMS_XX'.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M, (2009-04-14)
  !!
  SUBROUTINE init_hybrid_coord(nlev, vct_a, vct_b)

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_nh_init_utils:init_hybrid_coord'

    INTEGER,  INTENT(IN)    :: nlev  !< number of full levels
    REAL(wp), INTENT(INOUT) :: vct_a(:), vct_b(:)

    REAL(wp) :: z_height, z_flat
    INTEGER  :: jk
    INTEGER  :: nlevp1            !< number of half levels

    ! number of vertical half levels
    nlevp1 = nlev+1

    ! read hybrid parameters as in the hydrostatic model
    IF ( layer_thickness < 0.0_wp) THEN

      CALL read_vct (iequations,nlev)

    ELSE

      z_flat = REAL(nlevp1-n_flat_level,wp) * layer_thickness
      DO jk = 1, nlevp1
        z_height  = layer_thickness*REAL(nlevp1-jk,wp)
        vct_a(jk) = z_height
        IF ( z_height >= z_flat) THEN
          vct_b(jk) = 0.0_wp
        ELSE
          vct_b(jk) = (z_flat - z_height)/z_flat
        ENDIF
      ENDDO

    ENDIF

  END SUBROUTINE init_hybrid_coord


  !---------------------------------------------------------------------------
  !> Utility routine: computes some values based on vct_a, vct_b
  !
  SUBROUTINE prepare_hybrid_coord(nlev, vct_a, vct_b, vct, nflatlev)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_init_utils:init_hybrid_coord'

    INTEGER,  INTENT(IN)  :: nlev  !< number of full levels
    REAL(wp), INTENT(IN)  :: vct_a(:), vct_b(:)
    REAL(wp), INTENT(OUT) :: vct(:)
    INTEGER,  INTENT(OUT) :: nflatlev(:)

    REAL(wp) :: z_height, z_flat
    INTEGER  :: jk, nflat
    INTEGER  :: nlevp1            !< number of half levels

    ! number of vertical half levels
    nlevp1 = nlev+1

    IF ( layer_thickness < 0.0_wp) THEN

      CALL read_vct (iequations,nlev)

      DO jk = 1, nlevp1
        IF (vct_b(jk) /= 0.0_wp) THEN
          nflat = jk-1
          nflatlev(1) = nflat
          EXIT
        ENDIF
      ENDDO

    ELSE

      nflat  = -1
      z_flat = REAL(nlevp1-n_flat_level,wp) * layer_thickness
      DO jk = 1, nlevp1
        z_height  = layer_thickness*REAL(nlevp1-jk,wp)
        IF ( z_height < z_flat) THEN
          IF (nflat == -1) THEN
            nflat = jk-1
          ENDIF
        ENDIF
      ENDDO
      nflatlev(1) = nflat

      vct(       1:       nlevp1) = vct_a(:)
      vct(nlevp1+1:nlevp1+nlevp1) = vct_b(:)
    ENDIF

    CALL message(TRIM(routine), ' coordinate setup finished')

  END SUBROUTINE prepare_hybrid_coord

  !---------------------------------------------------------------------------
  !>
  !! Initialize SLEVE coordinate for nonhydrostatic model.
  !! In this initial version, the layer distribution is generated by an analytic
  !! formula based on a couple of namelist variables, but an option to read
  !! in a table may be added.
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2010-07-21)
  !!
  SUBROUTINE init_sleve_coord(nlev, vct_a, vct_b)

    INTEGER,  INTENT(IN)    :: nlev  !< number of full levels
    REAL(wp), INTENT(INOUT) :: vct_a(:), vct_b(:)

    REAL(wp) :: z_exp, dvct(nlev), zvcta(nlev+1), stretchfac, zdvct
    INTEGER  :: jk, jk1, jks, jk2
    INTEGER  :: nlevp1        !< number of full and half levels

    ! number of vertical levels
    nlevp1 = nlev+1

    IF (min_lay_thckn > 0.01_wp) THEN
      z_exp = LOG(min_lay_thckn/top_height)/LOG(2._wp/pi*ACOS(REAL(nlev-1,wp)**stretch_fac/&
        &     REAL(nlev,wp)**stretch_fac))

      ! Set up distribution of coordinate surfaces according to the analytical formula
      ! vct = h_top*(2/pi*arccos(jk-1/nlev))**z_exp (taken from the COSMO model, src_artifdata)
      ! z_exp has been calculated above in order to return min_lay_thckn as thickness
      ! of the lowest model layer
      DO jk = 1, nlevp1
        vct_a(jk)      = top_height*(2._wp/pi*ACOS(REAL(jk-1,wp)**stretch_fac/ &
          &              REAL(nlev,wp)**stretch_fac))**z_exp
        vct_b(jk)      = EXP(-vct_a(jk)/5000._wp)
      ENDDO
      ! Apply additional limitation on layer thickness in the middle and upper troposphere if the paramter
      ! max_lay_thckn is specified appropriately
      IF (max_lay_thckn > 2._wp*min_lay_thckn .AND. max_lay_thckn < 0.5_wp*htop_thcknlimit) THEN
        jk1 = 0
        DO jk = 1, nlev
          dvct(jk) = vct_a(jk) - vct_a(jk+1)
          IF (dvct(jk) > max_lay_thckn) jk1 = jk ! lowest layer in which the original layer thickness exceeds
                                                 ! the specified value
        ENDDO
        jks = 0
        jk2 = 0
        zvcta(nlevp1) = 0._wp
        DO jk = nlev, 1, -1
          IF (zvcta(jk+1) < htop_thcknlimit) THEN
            zvcta(jk) = zvcta(jk+1)+MIN(max_lay_thckn,dvct(jk))
          ELSE IF (jk2 == 0) THEN
            jk2 = jk+1
            jks = MAX(0,jk1-jk2)  ! shift layers from which thicknesses are taken downward in order to prevent sudden jumps
            zvcta(jk) = zvcta(jk+1)+dvct(jk+jks)
          ELSE
            zvcta(jk) = zvcta(jk+1)+dvct(jk+jks)
          ENDIF
        ENDDO
        IF (jks == 0) THEN ! either jk1 < htop_thcknlimit, which means that the thickness limiter has nothing to do,
                           ! or htop_thcknlimit is larger than the provisional model top;
          stretchfac = 1   ! in the latter case, the model top height is reset, overriding the top_height parameter
                           !
        ELSE               ! stretch remaining model levels such as to retain the original model top height
          stretchfac = (vct_a(1)-(zvcta(jk2)+REAL(jk2-1,wp)*max_lay_thckn))/&
                       (zvcta(1)-(zvcta(jk2)+REAL(jk2-1,wp)*max_lay_thckn))
        ENDIF
        DO jk = nlev, 1, -1
          IF (vct_a(jk+1) < htop_thcknlimit) THEN
            vct_a(jk) = vct_a(jk+1)+MIN(max_lay_thckn,dvct(jk))
          ELSE
            vct_a(jk) = vct_a(jk+1)+max_lay_thckn+(dvct(jk+jks)-max_lay_thckn)*stretchfac
          ENDIF
        ENDDO

        ! Try to apply additional smoothing on the stretching factor above the constant-thickness layer
        IF (stretchfac /= 1._wp .AND. jk1 < nlev-3) THEN
          DO jk = nlev, 1, -1
            IF (zvcta(jk+1) < htop_thcknlimit) THEN
              zvcta(jk) = vct_a(jk)
            ELSE
              zdvct = MIN(1.025_wp*(vct_a(jk)-vct_a(jk+1)), 1.025_wp*(zvcta(jk1+1)-zvcta(jk1+2))/ &
                         (zvcta(jk1+2)-zvcta(jk1+3))*(zvcta(jk+1)-zvcta(jk+2)) )
              zvcta(jk) = MIN(vct_a(jk),zvcta(jk+1)+zdvct)
            ENDIF
          ENDDO
          IF (zvcta(1) == vct_a(1)) THEN
            vct_a(1:2) = zvcta(1:2)
            vct_a(jk2+1:nlev) = zvcta(jk2+1:nlev)
            DO jk = 3, jk2
              vct_a(jk) = 0.5_wp*(zvcta(jk-1)+zvcta(jk+1))
            ENDDO
          ENDIF
        ENDIF

      ENDIF
    ELSE
     ! Use constant layer thicknesses determined by nlev and top_height
      DO jk = 1, nlevp1
        vct_a(jk) = top_height*(REAL(nlevp1,wp)-REAL(jk,wp))/REAL(nlev,wp)
        vct_b(jk) = EXP(-vct_a(jk)/5000._wp)
      ENDDO
    ENDIF

  END SUBROUTINE init_sleve_coord

  !---------------------------------------------------------------------------
  !> Utility routine: computes some values based on vct_a, vct_b
  !
  SUBROUTINE prepare_sleve_coord(nlev, vct_a, vct_b, vct, nflatlev)

    INTEGER,  INTENT(IN)  :: nlev  !< number of full levels
    REAL(wp), INTENT(IN)  :: vct_a(:), vct_b(:)
    REAL(wp), INTENT(OUT) :: vct(:)
    INTEGER,  INTENT(OUT) :: nflatlev(:)

    INTEGER  :: jk, nflat
    INTEGER  :: nlevp1        !< number of full and half levels

    ! number of vertical levels
    nlevp1 = nlev+1

    !------------------------------
    ! derived quantities: input: vct_a, vct_b -> output: vct, nflat, nflatlev

    IF (min_lay_thckn > 0.01_wp) THEN
      DO jk = 1, nlevp1
        vct(jk)        = vct_a(jk)
        vct(jk+nlevp1) = vct_b(jk)
      ENDDO
    ELSE
      ! Use constant layer thicknesses determined by nlev and top_height
      DO jk = 1, nlevp1
        vct(jk)        = vct_a(jk)
        vct(jk+nlevp1) = vct_b(jk)
      ENDDO
    ENDIF

    ! Determine nflat (first model level for which coordinate surfaces have a
    ! terrain-following component)
    DO jk = 1, nlevp1
      IF (vct_a(jk) < flat_height) THEN
        nflat = jk-1
        EXIT
      ENDIF
    ENDDO

    ! nflat must not be zero for global domain
    nflat = MAX(1,nflat)
    nflatlev(1) = nflat

    CALL message('mo_nh_init_utils: init_sleve_coord', '')

    IF (msg_level >= 7) THEN
     WRITE(message_text,'(a)') 'Nominal heights of coordinate half levels and layer thicknesses (m):'
        CALL message('', TRIM(message_text))

      DO jk = 1, nlevp1
       WRITE(message_text,'(a,i4,2F12.3)') 'jk, vct_a, dvct: ',jk, vct_a(jk), vct_a(jk)-vct_a(MIN(jk+1,nlevp1))
        CALL message('', TRIM(message_text))
      ENDDO
    ENDIF

  END SUBROUTINE prepare_sleve_coord

  !---------------------------------------------------------------------------
  !>
  !! Computes the smoothed topography needed for the SLEVE coordinate.
  !! May be bypassed once an option for reading the smooth topography from data
  !! is available
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2010-07-21)
  !!
  SUBROUTINE compute_smooth_topo(p_patch, p_int, topo_c, topo_smt_c)

    TYPE(t_patch),TARGET,INTENT(IN) :: p_patch
    TYPE(t_int_state), INTENT(IN) :: p_int

    ! Input fields: topography on cells
    REAL(wp), INTENT(IN) :: topo_c(:,:)

    ! Output fields: smooth topography on cells
    REAL(wp), INTENT(OUT) :: topo_smt_c(:,:)

    INTEGER  :: jb, jc, iter, niter
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),nabla2_topo(nproma,1,p_patch%nblks_c)

    !-------------------------------------------------------------------------

    niter = 25 ! 25 smoothing iterations (do we need this to be a namelist variable?)

    ! Initialize auxiliary fields for topography with data and nullify nabla2 field
    z_topo(:,1,:)      = topo_c(:,:)
    nabla2_topo(:,1,:) = 0._wp

    i_startblk = p_patch%cells%start_blk(2,1)
    nblks_c    = p_patch%nblks_c

    CALL sync_patch_array(SYNC_C,p_patch,z_topo)

    ! Apply nabla2-diffusion niter times to create smooth topography
    DO iter = 1, niter

      CALL nabla2_scalar(z_topo, p_patch, p_int, nabla2_topo, 1, 1)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 2)

        DO jc = i_startidx, i_endidx
          z_topo(jc,1,jb) = z_topo(jc,1,jb) + 0.125_wp*nabla2_topo(jc,1,jb) &
            &                               * p_patch%cells%area(jc,jb)
        ENDDO
      ENDDO

      CALL sync_patch_array(SYNC_C,p_patch,z_topo)

    ENDDO

    ! Store smooth topography on output fields
    topo_smt_c(:,:) = z_topo(:,1,:)

  END SUBROUTINE compute_smooth_topo

  !----------------------------------------------------------------------------
  !>
  !! Computes the 3D vertical coordinate fields for the nonhydrostatic model.
  !! (was originally included in subroutine set_nh_metrics but has been
  !! encapsulated because IFS2ICON needs the coordinate fields as input.
  !! Note: this routine is supposed to be used for cells and vertices.
  !! Therefore, field dimensions are passed.
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2011-07-01)
  !!
  SUBROUTINE init_vert_coord(topo, topo_smt, z3d_i,     &
                             nlev, nblks, npromz, nshift, nflat)

    ! Input parameters:
    INTEGER, INTENT(IN) :: nlev, nblks, & ! field dimensions
                           npromz,      & ! length of last block
                           nshift,      & ! shift parameter for vertical nesting
                           nflat          ! index below which model levels are flat

    ! Input fields: "Normal" and smooth topography
    REAL(wp),  INTENT(IN) :: topo    (nproma,nblks), &
                             topo_smt(nproma,nblks)


    ! Output fields: 3D coordinate fields at interface and main levels
    REAL(wp),  INTENT(OUT) :: z3d_i(nproma,nlev+1,nblks)

    INTEGER :: jc, jk, jk1, jb, nlen, nlevp1, ierr(nblks), nerr, ktop_thicklimit(nproma)
    REAL(wp) :: z_fac1, z_fac2, z_topo_dev(nproma), min_lay_spacing, &
                dvct, dvct1, dvct2, minrat1, minrat2, wfac, dz1, dz2, dz3, dzr
    !-------------------------------------------------------------------------

    nlevp1  = nlev+1
    dvct1   = 100._wp
    minrat1 = 1._wp/3._wp ! minimum relative layer thickness for nominal thicknesses <= dvct1 (in m)
    dvct2   = 500._wp
    minrat2 = 0.5_wp  ! minimum relative layer thickness for a nominal thickness of dvct2
    ierr(:) = 0

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jk1, z_fac1, z_fac2, z_topo_dev, min_lay_spacing,&
!$OMP dvct, wfac, ktop_thicklimit, dz1, dz2, dz3, dzr) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1,nblks

      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        z3d_i(nlen+1:nproma,:,jb) = 0._wp
     ENDIF

     z3d_i(1:nlen,nlevp1,jb) = topo(1:nlen,jb)
     ktop_thicklimit(:) = nlevp1

     ! vertical interface height
     IF (ivctype == 1 .OR. decay_scale_1 >= 0.5_wp*top_height) THEN ! hybrid Gal-Chen coordinate
       DO jk = 1, nlev
         jk1 = jk + nshift
         z3d_i(1:nlen,jk,jb) = vct_a(jk1) + vct_b(jk1)*z3d_i(1:nlen,nlevp1,jb)
       ENDDO
     ELSE IF (ivctype == 2) THEN ! SLEVE coordinate (Leuenberger et al. MWR 2010)
       DO jk = 1, nflat
         jk1 = jk + nshift
         z3d_i(1:nlen,jk,jb) = vct_a(jk1)
       ENDDO
       DO jk = nflat + 1, nlev
         jk1 = jk + nshift
         ! Scaling factors for large-scale and small-scale topography
         z_fac1 = SINH((top_height/decay_scale_1)**decay_exp - &
                       (vct_a(jk1)/decay_scale_1)**decay_exp)/ &
                  SINH((top_height/decay_scale_1)**decay_exp)
         z_fac2 = SINH((top_height/decay_scale_2)**decay_exp - &
                       (vct_a(jk1)/decay_scale_2)**decay_exp)/ &
                  SINH((top_height/decay_scale_2)**decay_exp)

         ! Small-scale topography (i.e. full topo - smooth topo)
         z_topo_dev(1:nlen) = topo(1:nlen,jb) - topo_smt(1:nlen,jb)

         z3d_i(1:nlen,jk,jb) = vct_a(jk1) + topo_smt(1:nlen,jb)*z_fac1 + &
            z_topo_dev(1:nlen)*z_fac2
       ENDDO
       ! Ensure that layer thicknesses are not too small; this would potentially cause
       ! instabilities in vertical advection
       DO jk = nlev, 1, -1
         jk1 = jk + nshift
         dvct = vct_a(jk1) - vct_a(jk1+1)
         IF (dvct < dvct1) THEN ! limit layer thickness to minrat1 times its nominal value
           min_lay_spacing = minrat1*dvct
         ELSE IF (dvct < dvct2) THEN ! limitation factor changes from minrat1 to minrat2
           wfac = ((dvct2-dvct)/(dvct2-dvct1))**2
           min_lay_spacing = (minrat1*wfac + minrat2*(1._wp-wfac))*dvct
         ELSE ! limitation factor decreases again
           min_lay_spacing = minrat2*dvct2*(dvct/dvct2)**(1._wp/3._wp)
         ENDIF
         min_lay_spacing = MAX(min_lay_spacing,MIN(50._wp,min_lay_thckn))
         DO jc = 1, nlen
           IF (z3d_i(jc,jk+1,jb)+min_lay_spacing > z3d_i(jc,jk,jb)) THEN
             z3d_i(jc,jk,jb) = z3d_i(jc,jk+1,jb) + min_lay_spacing
             ktop_thicklimit(jc) = jk
           ENDIF
         ENDDO
       ENDDO
       ! Smooth layer thickness ratios in the transition layer of columns where the thickness limiter has been active
       DO jc = 1, nlen
         jk = ktop_thicklimit(jc)
         IF (jk <= nlev-2 .AND. jk >= 4) THEN
           ! TODO : array access with subscript (jk-3)=0
           dz1 = z3d_i(jc,jk+1,jb)-z3d_i(jc,jk+2,jb)
           dz2 = z3d_i(jc,jk-3,jb)-z3d_i(jc,jk-2,jb)
           dzr = (dz2/dz1)**0.25_wp ! stretching factor
           dz3 = (z3d_i(jc,jk-2,jb)-z3d_i(jc,jk+1,jb))/(dzr*(1._wp+dzr*(1._wp+dzr)))
           z3d_i(jc,jk  ,jb) = MAX(z3d_i(jc,jk  ,jb), z3d_i(jc,jk+1,jb) + dz3*dzr)
           z3d_i(jc,jk-1,jb) = MAX(z3d_i(jc,jk-1,jb), z3d_i(jc,jk  ,jb) + dz3*dzr*dzr)
         ENDIF
       ENDDO
       ! Check if level nflat is still flat
       IF (ANY(z3d_i(1:nlen,nflat,jb) /= vct_a(nflat+nshift))) ierr(jb) = 1
       ! Check also if ktop_thicklimit is sufficiently far away from the model top
       IF (nlev > 6 .AND. ANY(ktop_thicklimit(1:nlen) <= 3)) ierr(jb) = ierr(jb) + 1
     ENDIF

   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   nerr = SUM(ierr(1:nblks))
   IF (nerr > 0) CALL finish ('init_vert_coord: ', &
      'flat_height in sleve_nml or model top is too low')

  END SUBROUTINE init_vert_coord

  !----------------------------------------------------------------------------
  !>
  !! Saves the initial state of NWP applications for the IAU iteration mode.
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2016-06-17)
  !!
  SUBROUTINE save_initial_state(p_patch, p_nh, prm_diag, p_lnd, ext_data)

    TYPE(t_patch),             INTENT(IN) :: p_patch(:)
    TYPE(t_nh_state),          INTENT(IN) :: p_nh(:)
    TYPE(t_nwp_phy_diag),      INTENT(IN) :: prm_diag(:)
    TYPE(t_lnd_state), TARGET, INTENT(IN) :: p_lnd(:)
    TYPE(t_external_data),     INTENT(IN) :: ext_data(:)

    INTEGER :: jg, ntl, ntw, nlev, nlevp1, nblks_c, nblks_e

    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog

    ntl = ntiles_total
    ntw = ntiles_total+ntiles_water

    ALLOCATE(saveinit(n_dom))

    DO jg = 1, n_dom

      IF(.NOT. p_patch(jg)%ldom_active) CYCLE

      nlev    = p_patch(jg)%nlev
      nlevp1  = p_patch(jg)%nlevp1
      nblks_c = p_patch(jg)%nblks_c
      nblks_e = p_patch(jg)%nblks_e

      lnd_prog => p_lnd(jg)%prog_lnd(nnow_rcf(jg))
      lnd_diag => p_lnd(jg)%diag_lnd
      wtr_prog => p_lnd(jg)%prog_wtr(nnow_rcf(jg))

      ALLOCATE (saveinit(jg)%fr_seaice(nproma,nblks_c), saveinit(jg)%t_ice(nproma,nblks_c),    &
                saveinit(jg)%h_ice(nproma,nblks_c),     saveinit(jg)%gz0(nproma,nblks_c),      &
                saveinit(jg)%t_mnw_lk(nproma,nblks_c),  saveinit(jg)%t_wml_lk(nproma,nblks_c), &
                saveinit(jg)%h_ml_lk(nproma,nblks_c),   saveinit(jg)%t_bot_lk(nproma,nblks_c), &
                saveinit(jg)%c_t_lk(nproma,nblks_c),    saveinit(jg)%t_b1_lk(nproma,nblks_c),  &
                saveinit(jg)%h_b1_lk(nproma,nblks_c) )

      ALLOCATE (saveinit(jg)%theta_v(nproma,nlev,nblks_c), &
                saveinit(jg)%rho(nproma,nlev,nblks_c),     &
                saveinit(jg)%exner(nproma,nlev,nblks_c),   &
                saveinit(jg)%w(nproma,nlevp1,nblks_c),     &
                saveinit(jg)%tke(nproma,nlevp1,nblks_c),   &
                saveinit(jg)%vn(nproma,nlev,nblks_e),      &
                saveinit(jg)%t_g_t(nproma,nblks_c,ntw),    &
                saveinit(jg)%qv_s_t(nproma,nblks_c,ntw),   &
                saveinit(jg)%freshsnow_t(nproma,nblks_c,ntl), &
                saveinit(jg)%snowfrac_t(nproma,nblks_c,ntl), &
                saveinit(jg)%snowfrac_lc_t(nproma,nblks_c,ntl), &
                saveinit(jg)%w_snow_t(nproma,nblks_c,ntl), &
                saveinit(jg)%w_i_t(nproma,nblks_c,ntl),    &
                saveinit(jg)%h_snow_t(nproma,nblks_c,ntl), &
                saveinit(jg)%t_snow_t(nproma,nblks_c,ntl), &
                saveinit(jg)%rho_snow_t(nproma,nblks_c,ntl), &
                saveinit(jg)%snowtile_flag_t(nproma,nblks_c,ntl), &
                saveinit(jg)%idx_lst_t(nproma,nblks_c,ntl), &
                saveinit(jg)%frac_t(nproma,nblks_c,ntw),    &
                saveinit(jg)%gp_count_t(nblks_c,ntl)        )

      ALLOCATE (saveinit(jg)%tracer(nproma,nlev,nblks_c,ntracer),      &
                saveinit(jg)%w_so_t(nproma,nlev_soil,nblks_c,ntl),     &
                saveinit(jg)%w_so_ice_t(nproma,nlev_soil,nblks_c,ntl), &
                saveinit(jg)%t_so_t(nproma,nlev_soil+1,nblks_c,ntl)    )

      IF (lmulti_snow) THEN
        ALLOCATE (saveinit(jg)%t_snow_mult_t(nproma,nlev_snow+1,nblks_c,ntl), &
                  saveinit(jg)%rho_snow_mult_t(nproma,nlev_snow,nblks_c,ntl), &
                  saveinit(jg)%wtot_snow_t(nproma,nlev_snow,nblks_c,ntl),     &
                  saveinit(jg)%wliq_snow_t(nproma,nlev_snow,nblks_c,ntl),     &
                  saveinit(jg)%dzh_snow_t(nproma,nlev_snow,nblks_c,ntl)       )
      ELSE IF (l2lay_rho_snow) THEN
        ALLOCATE (saveinit(jg)%rho_snow_mult_t(nproma,nlev_snow,nblks_c,ntl))
      ENDIF

      IF (iprog_aero == 1) ALLOCATE (saveinit(jg)%aerosol(nproma,nclass_aero,nblks_c))

!$OMP PARALLEL
      CALL copy(lnd_diag%fr_seaice, saveinit(jg)%fr_seaice)
      CALL copy(wtr_prog%t_ice, saveinit(jg)%t_ice)
      CALL copy(wtr_prog%h_ice, saveinit(jg)%h_ice)
      CALL copy(prm_diag(jg)%gz0, saveinit(jg)%gz0)
      CALL copy(wtr_prog%t_mnw_lk, saveinit(jg)%t_mnw_lk)
      CALL copy(wtr_prog%t_wml_lk, saveinit(jg)%t_wml_lk)
      CALL copy(wtr_prog%h_ml_lk, saveinit(jg)%h_ml_lk)
      CALL copy(wtr_prog%t_bot_lk, saveinit(jg)%t_bot_lk)
      CALL copy(wtr_prog%c_t_lk, saveinit(jg)%c_t_lk)
      CALL copy(wtr_prog%t_b1_lk, saveinit(jg)%t_b1_lk)
      CALL copy(wtr_prog%h_b1_lk, saveinit(jg)%h_b1_lk)

      CALL copy(p_nh(jg)%prog(nnow(jg))%theta_v, saveinit(jg)%theta_v)
      CALL copy(p_nh(jg)%prog(nnow(jg))%rho, saveinit(jg)%rho)
      CALL copy(p_nh(jg)%prog(nnow(jg))%exner, saveinit(jg)%exner)
      CALL copy(p_nh(jg)%prog(nnow(jg))%w, saveinit(jg)%w)
      CALL copy(p_nh(jg)%prog(nnow_rcf(jg))%tke, saveinit(jg)%tke)
      CALL copy(p_nh(jg)%prog(nnow(jg))%vn, saveinit(jg)%vn)
      CALL copy(p_nh(jg)%prog(nnow_rcf(jg))%tracer, saveinit(jg)%tracer)

      CALL copy(lnd_prog%t_g_t, saveinit(jg)%t_g_t)
      CALL copy(lnd_diag%qv_s_t, saveinit(jg)%qv_s_t)
      CALL copy(lnd_diag%freshsnow_t, saveinit(jg)%freshsnow_t)
      CALL copy(lnd_diag%snowfrac_t, saveinit(jg)%snowfrac_t)
      CALL copy(lnd_diag%snowfrac_lc_t, saveinit(jg)%snowfrac_lc_t)
      CALL copy(lnd_prog%w_snow_t, saveinit(jg)%w_snow_t)
      CALL copy(lnd_prog%w_i_t, saveinit(jg)%w_i_t)
      CALL copy(lnd_diag%h_snow_t, saveinit(jg)%h_snow_t)
      CALL copy(lnd_prog%t_snow_t, saveinit(jg)%t_snow_t)
      CALL copy(lnd_prog%rho_snow_t, saveinit(jg)%rho_snow_t)
      CALL copy(lnd_prog%w_so_t, saveinit(jg)%w_so_t)
      CALL copy(lnd_prog%w_so_ice_t, saveinit(jg)%w_so_ice_t)
      CALL copy(lnd_prog%t_so_t, saveinit(jg)%t_so_t)

      IF (ntiles_total > 1 .AND. lsnowtile .AND. .NOT. ltile_coldstart) THEN
        CALL copy(ext_data(jg)%atm%snowtile_flag_t, saveinit(jg)%snowtile_flag_t)
        CALL copy(ext_data(jg)%atm%idx_lst_t, saveinit(jg)%idx_lst_t)
        CALL copy(ext_data(jg)%atm%frac_t, saveinit(jg)%frac_t)
        CALL copy(ext_data(jg)%atm%gp_count_t, saveinit(jg)%gp_count_t)
      ENDIF

      IF (lmulti_snow) THEN
        CALL copy(lnd_prog%t_snow_mult_t, saveinit(jg)%t_snow_mult_t)
        CALL copy(lnd_prog%rho_snow_mult_t, saveinit(jg)%rho_snow_mult_t)
        CALL copy(lnd_prog%wtot_snow_t, saveinit(jg)%wtot_snow_t)
        CALL copy(lnd_prog%wliq_snow_t, saveinit(jg)%wliq_snow_t)
        CALL copy(lnd_prog%dzh_snow_t, saveinit(jg)%dzh_snow_t)
      ELSE IF (l2lay_rho_snow) THEN
        CALL copy(lnd_prog%rho_snow_mult_t, saveinit(jg)%rho_snow_mult_t)
      ENDIF

      IF (iprog_aero == 1)  CALL copy(prm_diag(jg)%aerosol, saveinit(jg)%aerosol)
!$OMP END PARALLEL

    ENDDO

  END SUBROUTINE save_initial_state

  !----------------------------------------------------------------------------
  !>
  !! Restores the initial state of NWP applications for the IAU iteration mode.
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2016-06-17)
  !!
  SUBROUTINE restore_initial_state(p_patch, p_nh, prm_diag, prm_tend, p_lnd, ext_data)

    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state),          INTENT(INOUT) :: p_nh(:)
    TYPE(t_nwp_phy_diag),      INTENT(INOUT) :: prm_diag(:)
    TYPE(t_nwp_phy_tend),      INTENT(INOUT) :: prm_tend(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd(:)
    TYPE(t_external_data),     INTENT(INOUT) :: ext_data(:)

    INTEGER :: jg, ntl, ntw, nlev, nlevp1, nblks_c, nblks_e

    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog

    ntl = ntiles_total
    ntw = ntiles_total+ntiles_water


    DO jg = 1, n_dom

      IF(.NOT. p_patch(jg)%ldom_active) CYCLE

      nlev    = p_patch(jg)%nlev
      nlevp1  = p_patch(jg)%nlevp1
      nblks_c = p_patch(jg)%nblks_c
      nblks_e = p_patch(jg)%nblks_e

      lnd_prog => p_lnd(jg)%prog_lnd(nnow_rcf(jg))
      lnd_diag => p_lnd(jg)%diag_lnd
      wtr_prog => p_lnd(jg)%prog_wtr(nnow_rcf(jg))


!$OMP PARALLEL
      CALL copy(saveinit(jg)%fr_seaice, lnd_diag%fr_seaice)
      CALL copy(saveinit(jg)%t_ice, wtr_prog%t_ice)
      CALL copy(saveinit(jg)%h_ice, wtr_prog%h_ice)
      CALL copy(saveinit(jg)%gz0, prm_diag(jg)%gz0)
      CALL copy(saveinit(jg)%t_mnw_lk, wtr_prog%t_mnw_lk)
      CALL copy(saveinit(jg)%t_wml_lk, wtr_prog%t_wml_lk)
      CALL copy(saveinit(jg)%h_ml_lk, wtr_prog%h_ml_lk)
      CALL copy(saveinit(jg)%t_bot_lk, wtr_prog%t_bot_lk)
      CALL copy(saveinit(jg)%c_t_lk, wtr_prog%c_t_lk)
      CALL copy(saveinit(jg)%t_b1_lk, wtr_prog%t_b1_lk)
      CALL copy(saveinit(jg)%h_b1_lk, wtr_prog%h_b1_lk)

      CALL copy(saveinit(jg)%theta_v, p_nh(jg)%prog(nnow(jg))%theta_v)
      CALL copy(saveinit(jg)%rho, p_nh(jg)%prog(nnow(jg))%rho)
      CALL copy(saveinit(jg)%exner, p_nh(jg)%prog(nnow(jg))%exner)
      CALL copy(saveinit(jg)%w, p_nh(jg)%prog(nnow(jg))%w)
      CALL copy(saveinit(jg)%tke, p_nh(jg)%prog(nnow_rcf(jg))%tke)
      CALL copy(saveinit(jg)%vn, p_nh(jg)%prog(nnow(jg))%vn)
      CALL copy(saveinit(jg)%tracer, p_nh(jg)%prog(nnow_rcf(jg))%tracer)

      CALL copy(saveinit(jg)%t_g_t, lnd_prog%t_g_t)
      CALL copy(saveinit(jg)%qv_s_t, lnd_diag%qv_s_t)
      CALL copy(saveinit(jg)%freshsnow_t, lnd_diag%freshsnow_t)
      CALL copy(saveinit(jg)%snowfrac_t, lnd_diag%snowfrac_t)
      CALL copy(saveinit(jg)%snowfrac_lc_t, lnd_diag%snowfrac_lc_t)
      CALL copy(saveinit(jg)%w_snow_t, lnd_prog%w_snow_t)
      CALL copy(saveinit(jg)%w_i_t, lnd_prog%w_i_t)
      CALL copy(saveinit(jg)%h_snow_t, lnd_diag%h_snow_t)
      CALL copy(saveinit(jg)%t_snow_t, lnd_prog%t_snow_t)
      CALL copy(saveinit(jg)%rho_snow_t, lnd_prog%rho_snow_t)
      CALL copy(saveinit(jg)%w_so_t, lnd_prog%w_so_t)
      CALL copy(saveinit(jg)%w_so_ice_t, lnd_prog%w_so_ice_t)
      CALL copy(saveinit(jg)%t_so_t, lnd_prog%t_so_t)

      IF (ntiles_total > 1 .AND. lsnowtile .AND. .NOT. ltile_coldstart) THEN
        CALL copy(saveinit(jg)%snowtile_flag_t, ext_data(jg)%atm%snowtile_flag_t)
        CALL copy(saveinit(jg)%idx_lst_t, ext_data(jg)%atm%idx_lst_t)
        CALL copy(saveinit(jg)%frac_t, ext_data(jg)%atm%frac_t)
        CALL copy(saveinit(jg)%gp_count_t, ext_data(jg)%atm%gp_count_t)
      ENDIF

      IF (lmulti_snow) THEN
        CALL copy(saveinit(jg)%t_snow_mult_t, lnd_prog%t_snow_mult_t)
        CALL copy(saveinit(jg)%rho_snow_mult_t, lnd_prog%rho_snow_mult_t)
        CALL copy(saveinit(jg)%wtot_snow_t, lnd_prog%wtot_snow_t)
        CALL copy(saveinit(jg)%wliq_snow_t, lnd_prog%wliq_snow_t)
        CALL copy(saveinit(jg)%dzh_snow_t, lnd_prog%dzh_snow_t)
      ELSE IF (l2lay_rho_snow) THEN
        CALL copy(saveinit(jg)%rho_snow_mult_t, lnd_prog%rho_snow_mult_t)
      ENDIF

      IF (iprog_aero == 1)  CALL copy(saveinit(jg)%aerosol, prm_diag(jg)%aerosol)

      ! Fields that need to be reset to zero in order to obtain identical results
      CALL init (p_nh(jg)%diag%ddt_vn_phy)
      CALL init (p_nh(jg)%diag%ddt_tracer_adv)
      CALL init (prm_tend(jg)%ddt_tracer_turb)
      CALL init (prm_tend(jg)%ddt_temp_radsw)
      CALL init (prm_tend(jg)%ddt_temp_radlw)
      CALL init (prm_tend(jg)%ddt_temp_turb)
      CALL init (p_nh(jg)%diag%ddt_temp_dyn)
      CALL init (p_nh(jg)%diag%exner_dyn_incr)
      CALL init (prm_diag(jg)%rain_gsp_rate)
      CALL init (prm_diag(jg)%snow_gsp_rate)
      CALL init (prm_diag(jg)%shfl_s_t)
      CALL init (prm_diag(jg)%qhfl_s_t)
      CALL init (lnd_diag%runoff_s_t)
      CALL init (lnd_diag%runoff_g_t)

!$OMP END PARALLEL



      DEALLOCATE (saveinit(jg)%fr_seaice, saveinit(jg)%t_ice, saveinit(jg)%h_ice, saveinit(jg)%gz0,          &
                  saveinit(jg)%t_mnw_lk, saveinit(jg)%t_wml_lk, saveinit(jg)%h_ml_lk, saveinit(jg)%t_bot_lk, &
                  saveinit(jg)%c_t_lk, saveinit(jg)%t_b1_lk, saveinit(jg)%h_b1_lk )

      DEALLOCATE (saveinit(jg)%theta_v, saveinit(jg)%rho,saveinit(jg)%exner, saveinit(jg)%w, saveinit(jg)%tke,      &
                  saveinit(jg)%vn, saveinit(jg)%t_g_t, saveinit(jg)%qv_s_t, saveinit(jg)%freshsnow_t,               &
                  saveinit(jg)%snowfrac_t, saveinit(jg)%snowfrac_lc_t, saveinit(jg)%w_snow_t,                       &
                  saveinit(jg)%w_i_t, saveinit(jg)%h_snow_t, saveinit(jg)%t_snow_t, saveinit(jg)%rho_snow_t,        &
                  saveinit(jg)%snowtile_flag_t, saveinit(jg)%idx_lst_t, saveinit(jg)%frac_t, saveinit(jg)%gp_count_t)

      DEALLOCATE (saveinit(jg)%tracer, saveinit(jg)%w_so_t, saveinit(jg)%w_so_ice_t, saveinit(jg)%t_so_t)

      IF (lmulti_snow) THEN
        DEALLOCATE (saveinit(jg)%t_snow_mult_t, saveinit(jg)%rho_snow_mult_t, saveinit(jg)%wtot_snow_t, &
                   saveinit(jg)%wliq_snow_t, saveinit(jg)%dzh_snow_t)
      ELSE IF (l2lay_rho_snow) THEN
        DEALLOCATE (saveinit(jg)%rho_snow_mult_t)
      ENDIF

      IF (iprog_aero == 1) DEALLOCATE (saveinit(jg)%aerosol)
    ENDDO

    DEALLOCATE(saveinit)

  END SUBROUTINE restore_initial_state

  !>
  !! Compute weights for incremental analysis update
  !!
  !! Compute weights for incremental analysis update.
  !! 2 weights are provided:
  !! - iau_wgt_dyn can be used for all fields that need to be updated
  !!   every (fast) dynamics time step
  !! - iau_wgt_adv can be used for all fields that need to be updated
  !!   every (slow) advection time step.
  !!
  !! @par Revision History
  !! Initial revision by daniel Reinert, DWD (2014-01-29)
  !!
  SUBROUTINE compute_iau_wgt(sim_time, dt, dt_iau, lreset_wgt_adv)

    REAL(wp)        , INTENT(IN)  :: sim_time          !< Simulation time since model
                                                       !< start
    REAL(wp)        , INTENT(IN)  :: dt                !< time step
    REAL(wp)        , INTENT(IN)  :: dt_iau            !< width of IAU window
    LOGICAL         , INTENT(IN)  :: lreset_wgt_adv    !< If true, reset the accumulated weight for the advective time step

    ! local variables
    REAL(wp)  :: time_iau_elapsed                      !< elapsed time since IAU start [s]
    REAL(wp)  :: fct_eval                              !< result of top-hat or sin2 function

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nh_init_utils:compute_iau_wgt'
    !-------------------------------------------------------------------------

    ! initialize
    IF (lreset_wgt_adv) iau_wgt_adv   = 0._wp

    ! compute elapsed time (in s) since IAU start
    !
    ! trivial so far, however will be changed to mtime when the functionality of
    ! computing the timedelta between two dates becomes available.
    time_iau_elapsed = sim_time


    IF (time_iau_elapsed <= dt_iau) THEN
      is_iau_active = .TRUE.

      SELECT CASE (type_iau_wgt)
        CASE(1)  ! top-hat function
          fct_eval = iau_top_hat(dt_iau,time_iau_elapsed)

        CASE(2)  ! sin2 function
          fct_eval = iau_sin2   (dt_iau,time_iau_elapsed)

        CASE(3)  ! sin function
          fct_eval = iau_sin    (dt_iau,time_iau_elapsed)

        CASE default
          CALL finish(routine,&
                      'Invalid IAU weighting function. Must be 1, 2 or 3.')
      END SELECT

      ! compute weights by multiplying with the time step
      iau_wgt_dyn = fct_eval * dt
      iau_wgt_adv = iau_wgt_adv + iau_wgt_dyn

    ELSE
      is_iau_active = .FALSE.
      iau_wgt_dyn   = 0._wp
      iau_wgt_adv   = 0._wp
    ENDIF

!!$write(0,*) "sim_time, is_iau_active, iau_wgt_dyn, iau_wgt_adv: ", &
!!$  & sim_time, is_iau_active, iau_wgt_dyn, iau_wgt_adv

  END SUBROUTINE compute_iau_wgt


  !>
  !! Evaluates top-hat function at a particular point in time
  !!
  !! Evaluates top-hat function at a particular point in time
  !! Top-hat function is non-zero for 0<=t<=dt and is normalized such that
  !! \int_{t=0}^{t=dt} f(t)\,dt=1
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-29)
  !!
  FUNCTION iau_top_hat (dt,cur_time)  RESULT (fct_eval)

    REAL(wp), INTENT(IN) :: dt                 ! time interval [s]
    REAL(wp), INTENT(in) :: cur_time           ! current time  [s]

    REAL(wp) :: fct_eval
    !-------------------------------------------------------------------------

    IF (cur_time <= dt) THEN
      fct_eval = 1._wp/dt
    ELSE
      fct_eval = 0._wp
    ENDIF

  END FUNCTION iau_top_hat


  !>
  !! Evaluates SIN2 function at a particular point in time
  !!
  !! Evaluates SIN2 function at a particular point in time
  !! SIN2 function is non-zero for 0<=t<=dt and is normalized such that
  !! \int_{t=0}^{t=dt} f(t)\,dt=1
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-01-29)
  !!
  FUNCTION iau_sin2 (dt,cur_time)  RESULT (fct_eval)

    REAL(wp), INTENT(IN) :: dt                 ! time interval [s]
    REAL(wp), INTENT(in) :: cur_time           ! current time  [s]

    REAL(wp) :: fct_eval
    !-------------------------------------------------------------------------

    IF (cur_time <= dt) THEN
      fct_eval = (2._wp/dt) * SIN(pi*cur_time/dt)**2
    ELSE
      fct_eval = 0._wp
    ENDIF

  END FUNCTION iau_sin2

  !>
  !! Evaluates SIN function at a particular point in time
  !!
  !! Evaluates SIN function at a particular point in time
  !! SIN function is non-zero for 0<=t<=dt and is normalized such that
  !! \int_{t=0}^{t=dt} f(t)\,dt=1
  !!
  !! @par Revision History
  !! Initial revision by Harald Anlauf, DWD (2014-04-03)
  !!
  FUNCTION iau_sin (dt, cur_time)  RESULT (fct_eval)

    REAL(wp), INTENT(IN) :: dt                 ! time interval [s]
    REAL(wp), INTENT(in) :: cur_time           ! current time  [s]

    REAL(wp) :: fct_eval
    !-------------------------------------------------------------------------

    IF (cur_time <= dt) THEN
      fct_eval = ((PI/2._wp)/dt) * SIN(PI*cur_time/dt)
    ELSE
      fct_eval = 0._wp
    ENDIF

  END FUNCTION iau_sin

  !----------------------------------------------------------------------------
  ! Notes:
  ! - i refers to the interface between layes and m to the mid of a layer
  ! - two sets are available: for cell centres and vertices.
  ! - what to do on edge points?
  ! - need to be able two write/read grib and netcdf

!!$  SUBROUTINE write_vert_coord (nlev,                                 &
!!$       &                       nproma_c, nblks_c, z3d_ifc,  z3d_mfc, &
!!$       &                       nproma_v, nblks_v, z3d_ifv, z3d_mfv,  &
!!$       &                       hg_uuid, grid_no)
!!$    INTEGER,      INTENT(in) :: nlev, nproma_c, nblks_c, nproma_v, nblks_v
!!$    REAL(wp),     INTENT(in) :: z3d_ifc(nproma_c, nlev+1, nblks_c)
!!$    REAL(wp),     INTENT(in) :: z3d_mfc(nproma_c, nlev,   nblks_c)
!!$    REAL(wp),     INTENT(in) :: z3d_ifv(nproma_v, nlev+1, nblks_v)
!!$    REAL(wp),     INTENT(in) :: z3d_mfv(nproma_v, nlev,   nblks_v)
!!$    TYPE(t_uuid), INTENT(in) :: hg_uuid
!!$    INTEGER,      INTENT(in) :: grid_no
!!$
!!$    TYPE(t_uuid) :: uuid(4)
!!$    CHARACTER(len=uuid_string_length) :: uuid_string(4)
!!$
!!$    INTEGER :: n
!!$
!!$    ! get unique grid file identifier for GRIB2 and updated CF-Convention
!!$
!!$    DO n = 1, 4
!!$      CALL uuid_generate(uuid(n))
!!$      CALL uuid_unparse(uuid(n), uuid_string(n))
!!$    ENDDO
!!$
!!$  END SUBROUTINE write_vert_coord
!!$
!!$  SUBROUTINE read_vert_coord (nlev,                                 &
!!$       &                      nproma_c, nblks_c, z3d_ifc,  z3d_mfc, &
!!$       &                      nproma_v, nblks_v, z3d_ifv, z3d_mfv,  &
!!$       &                      hg_uuid, grid_no)
!!$    INTEGER,      INTENT(in)  :: nlev, nproma_c, nblks_c, nproma_v, nblks_v
!!$    REAL(wp),     INTENT(out) :: z3d_ifc(nproma_c, nlev+1, nblks_c)
!!$    REAL(wp),     INTENT(out) :: z3d_mfc(nproma_c, nlev,   nblks_c)
!!$    REAL(wp),     INTENT(out) :: z3d_ifv(nproma_v, nlev+1, nblks_v)
!!$    REAL(wp),     INTENT(out) :: z3d_mfv(nproma_v, nlev,   nblks_v)
!!$    TYPE(t_uuid), INTENT(out) :: hg_uuid
!!$    INTEGER,      INTENT(out) :: grid_no
!!$
!!$    TYPE(t_uuid) :: uuid(4)
!!$    CHARACTER(len=uuid_string_length) :: uuid_string(4)
!!$
!!$    INTEGER :: n
!!$
!!$    ! get unique grid file identifier for GRIB2 and updated CF-Convention
!!$
!!$    DO n = 1, 4
!!$      CALL uuid_parse(uuid_string(n), uuid(n))
!!$    ENDDO
!!$
!!$  END SUBROUTINE read_vert_coord
  !----------------------------------------------------------------------------

END MODULE mo_nh_init_utils
