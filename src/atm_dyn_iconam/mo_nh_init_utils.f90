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
MODULE mo_nh_init_utils

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_state,      ONLY: t_nh_metrics
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: msg_level
  USE mo_dynamics_config,     ONLY: iequations
  USE mo_physical_constants,  ONLY: grav, cpd, rd, cvd_o_rd, p0ref, vtmpc1
  USE mo_vertical_coord_table,ONLY: vct_a, vct_b, vct, read_vct
  USE mo_nonhydrostatic_config,ONLY: ivctype
  USE mo_sleve_config,        ONLY: min_lay_thckn, top_height, decay_scale_1, &
                                    decay_scale_2, decay_exp, flat_height, stretch_fac
  USE mo_impl_constants,      ONLY: max_dom, SUCCESS, min_rlcell, min_rlcell_int, &
                                    min_rlvert, min_rlvert_int
  USE mo_grf_nudgintp,        ONLY: interpol_scal_nudging
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf
  USE mo_math_constants,      ONLY: pi
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: p_pe,my_process_is_mpi_parallel
  USE mo_communication,       ONLY: exchange_data, exchange_data_mult
  USE mo_sync,                ONLY: sync_patch_array, SYNC_C, SYNC_V
  USE mo_interpolation,       ONLY: t_int_state, cells2verts_scalar, edges2cells_scalar
  USE mo_grf_interpolation,   ONLY: t_gridref_state, t_gridref_single_state 
  USE mo_math_operators,      ONLY: nabla2_scalar, grad_fd_norm
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_v, get_indices_e
  USE mo_interpol_config,     ONLY: nudge_zone_width
  USE mo_impl_constants_grf,  ONLY: grf_fbk_start_c, grf_bdywidth_c
  USE mo_subdivision,         ONLY: p_grf_state_local_parent, p_patch_local_parent, &
                                    p_int_state_local_parent


  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  INTEGER:: nflat, nflatlev(max_dom)

  PUBLIC :: nflat, nflatlev ! diese beiden Variablen muessen aus mo_vertical_grid raus!!!
  PUBLIC :: hydro_adjust, init_hybrid_coord, init_sleve_coord, compute_smooth_topo,    &
            init_vert_coord, topography_blending, topography_feedback, interp_uv_2_vn, &
            init_w, adjust_w, convert_thdvars, virtual_temp, convert_omega2w

CONTAINS
  !-------------
  !>
  !! SUBROUTINE hydro_adjust
  !! Computes hydrostatically balanced initial condition
  !!
  !! Input/Output: density, Exner pressure, virtual potential temperature
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-06-29)
  !!
  !!
  !!
  SUBROUTINE hydro_adjust(p_patch, p_nh_metrics,          &
                          rho, exner, theta_v, rhotheta_v )


    TYPE(t_patch),      INTENT(IN)       :: p_patch
    TYPE(t_nh_metrics), INTENT(IN)       :: p_nh_metrics

    ! Thermodynamic fields - all defined at full model levels
    REAL(wp), INTENT(INOUT) :: rho(:,:,:)        ! density (kg/m**3)
    REAL(wp), INTENT(INOUT) :: exner(:,:,:)      ! Exner pressure
    REAL(wp), INTENT(INOUT) :: theta_v(:,:,:)    ! virtual potential temperature (K)
    REAL(wp), INTENT(INOUT) :: rhotheta_v(:,:,:) ! rho*theta_v


    ! LOCAL VARIABLES
    REAL(wp) :: temp_v(nproma,p_patch%nlev) ! virtual temperature
    REAL(wp), DIMENSION(nproma) :: z_fac1, z_fac2, z_fac3, za, zb, zc  

    INTEGER :: jb, jk, jc
    INTEGER :: nlen, nlev

    nlev = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,temp_v,z_fac1,z_fac2,z_fac3,za,zb,zc)

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
        IF (p_patch%cell_type == 3) THEN
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
        ELSE IF (p_patch%cell_type == 6) THEN
          DO jc = 1, nlen
            z_fac3(jc) = grav/cpd*p_nh_metrics%ddqz_z_half(jc,jk+1,jb)
            z_fac1(jc) = exner(jc,jk+1,jb)
            z_fac2(jc) = exner(jc,jk+1,jb)/temp_v(jc,jk+1)
            za(jc) = 1.0_wp
            zb(jc) = -(z_fac2(jc)*(temp_v(jc,jk)+2.0_wp*z_fac3(jc))-z_fac1(jc))
            zc(jc) = temp_v(jc,jk)*z_fac2(jc)*z_fac1(jc)
          ENDDO
        ENDIF

        DO jc = 1, nlen
          exner(jc,jk,jb)      = (zb(jc)+SQRT(zb(jc)**2+4._wp*za(jc)*zc(jc)))/(2._wp*za(jc))
          theta_v(jc,jk,jb)    = temp_v(jc,jk)/exner(jc,jk,jb)
          rhotheta_v(jc,jk,jb) = exner(jc,jk,jb)**cvd_o_rd*p0ref/rd
          rho(jc,jk,jb)        = rhotheta_v(jc,jk,jb)/theta_v(jc,jk,jb)
        ENDDO

      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE hydro_adjust

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
!$OMP DO PRIVATE(jb,nlen,jk,jc)

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
!$OMP END DO
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
!$OMP DO PRIVATE(jb,nlen,jk,jc)

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
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE convert_omega2w
  !-------------
  !>
  !! SUBROUTINE virtual_temp
  !! Computes virtual temperature
  !!
  !! Required input fields: temperature, specific humidity, cloud and precipitation variables
  !! Output: virtual temperature
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE virtual_temp(p_patch, temp, qv, qc, qi, qr, qs, temp_v)


    TYPE(t_patch), INTENT(IN) :: p_patch

    ! Input fields - all defined at full model levels
    REAL(wp), INTENT(IN)           :: temp(:,:,:) ! temperature (K)
    REAL(wp), INTENT(IN)           :: qv  (:,:,:) ! specific humidity
    REAL(wp), INTENT(IN), OPTIONAL :: qc  (:,:,:) ! specific cloud water
    REAL(wp), INTENT(IN), OPTIONAL :: qi  (:,:,:) ! specific cloud ice
    REAL(wp), INTENT(IN), OPTIONAL :: qr  (:,:,:) ! specific rain water
    REAL(wp), INTENT(IN), OPTIONAL :: qs  (:,:,:) ! specific snow

    REAL(wp), INTENT(OUT) :: temp_v(:,:,:) ! virtual temperature (K)

    INTEGER :: jb, jk, jc
    INTEGER :: nlen, nlev
    LOGICAL :: l_cloud_precip

    nlev = SIZE(temp,2) ! in order to be usable for input and output data

    IF (PRESENT(qc) .AND. PRESENT(qi) .AND. PRESENT(qr) .AND. PRESENT(qs)) THEN
      l_cloud_precip = .TRUE.
    ELSE
      l_cloud_precip = .FALSE.
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc)

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      IF (l_cloud_precip) THEN
        DO jk = 1, nlev
          DO jc = 1, nlen
            temp_v(jc,jk,jb) = temp(jc,jk,jb) * (1._wp + vtmpc1*qv(jc,jk,jb) -      &
                              (qc(jc,jk,jb)+qi(jc,jk,jb)+qr(jc,jk,jb)+qs(jc,jk,jb)) )
          ENDDO
        ENDDO
      ELSE
        DO jk = 1, nlev
          DO jc = 1, nlen
            temp_v(jc,jk,jb) = temp(jc,jk,jb) * (1._wp + vtmpc1*qv(jc,jk,jb))
          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE virtual_temp

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
    INTEGER :: nlev, nblks_e, i_startblk,i_startidx, i_endidx

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

    nlev = p_patch%nlev
    nblks_e = p_patch%nblks_e

    i_startblk = p_patch%edges%start_blk(2,1)

    iidx => p_patch%edges%cell_idx
    iblk => p_patch%edges%cell_blk

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)

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
!$OMP END DO
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
    INTEGER :: jb, jk, je, jc
    INTEGER :: nlev, nlevp1, nblks_e, nblks_c, i_startblk, i_startidx, i_endidx

    REAL(wp) :: z_wsfc_e(nproma,1,p_patch%nblks_e) ! w at surface (edge points)
    REAL(wp) :: z_wsfc_c(nproma,1,p_patch%nblks_c) ! w at surface (cell points)
    REAL(wp) :: z_slope_e(nproma,p_patch%nlevp1,p_patch%nblks_e) ! slope at edges

    nlev    = p_patch%nlev
    nlevp1  = p_patch%nlevp1
    nblks_e = p_patch%nblks_e
    nblks_c = p_patch%nblks_c

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
!$OMP WORKSHARE
    w(:,:,:) = 0._wp
!$OMP END WORKSHARE

    ! specify a reasonable initial vertical wind speed
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 2)

      DO jc = i_startidx, i_endidx
        w(jc,nlevp1,jb) = z_wsfc_c(jc,1,jb)
      ENDDO
      DO jk = nlev, 2, -1
        DO jc = i_startidx, i_endidx
          w(jc,jk,jb) = z_wsfc_c(jc,1,jb)*vct_b(jk)**2
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
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
    INTEGER :: jb, jk, je, jc
    INTEGER :: nlev, nlevp1, nblks_e, nblks_c, i_startblk, i_startidx, i_endidx
    REAL(wp):: wfac

    REAL(wp) :: z_wsfc_e(nproma,1,p_patch%nblks_e) ! w at surface (edge points)
    REAL(wp) :: z_wsfc_c(nproma,1,p_patch%nblks_c) ! w at surface (cell points)
    REAL(wp) :: z_slope_e(nproma,p_patch%nlevp1,p_patch%nblks_e) ! slope at edges

    nlev    = p_patch%nlev
    nlevp1  = p_patch%nlevp1
    nblks_e = p_patch%nblks_e
    nblks_c = p_patch%nblks_c

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
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,wfac)
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 2)

      ! Lower boundary condition
      DO jc = i_startidx, i_endidx
        w(jc,nlevp1,jb) = z_wsfc_c(jc,1,jb)
      ENDDO

      ! Merging of lower boundary condition with interpolated data
      DO jk = nlev, 2, -1
        wfac = vct_b(jk)**6
        DO jc = i_startidx, i_endidx
          w(jc,jk,jb) = (1._wp-wfac)*w(jc,jk,jb) + wfac*z_wsfc_c(jc,1,jb)
        ENDDO
      ENDDO

      ! Upper boundary condition and smooth transition below
      DO jc = i_startidx, i_endidx
        w(jc,1,jb) = 0._wp
        w(jc,2,jb) = 0.33_wp*w(jc,2,jb)
        w(jc,3,jb) = 0.66_wp*w(jc,3,jb)
      ENDDO

    ENDDO
!$OMP END DO
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
  SUBROUTINE init_hybrid_coord(nlev)

    INTEGER, INTENT(IN) :: nlev  !< number of full levels
    INTEGER  :: jk
    INTEGER  :: nlevp1            !< number of half levels
    !-------------------------------------------------------------------------

    ! number of vertical half levels
    nlevp1 = nlev+1

    CALL read_vct (iequations, nlev)

    DO jk = 1, nlevp1
      IF (vct_b(jk) /= 0.0_wp) THEN
        nflat = jk-1
        nflatlev(1) = nflat
        EXIT
      ENDIF
    ENDDO

    CALL message('mo_nh_init_utils: init_hybrid_coord', '')

  END SUBROUTINE init_hybrid_coord

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
  SUBROUTINE init_sleve_coord(nlev)

    INTEGER, INTENT(IN) :: nlev  !< number of full levels
    REAL(wp) :: z_exp
    INTEGER  :: jk, ist
    INTEGER  :: nlevp1        !< number of full and half levels

    !-------------------------------------------------------------------------

    ! number of vertical levels
    nlevp1 = nlev+1


    ! Unlike the hybrid Gal-Chen coordinate, the SLEVE coordinate uses only one
    ! field with coordinate parameters
    ALLOCATE(vct_a(nlevp1), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish ('mo_nh_init_utils:init_sleve_coord', &
                   'allocation of vct_a failed')
    ENDIF

    ! Read namelist for SLEVE coordinate
!    is already done in the all-namelist read-routine
!    CALL sleve_nml_setup

    ! However, vct_b also needs to be defined because it is used elsewhere
    ALLOCATE(vct_b(nlevp1), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish ('mo_nh_init_utils:init_sleve_coord', &
                   'allocation of vct_b failed')
    ENDIF

    ! And vct is allocated because mo_io_vlist accesses it
    ALLOCATE(vct(nlevp1*2), STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish ('mo_nh_init_utils:init_hybrid_coord', &
                   'allocation of vct failed')
    ENDIF

    DO jk = 1, nlevp1
      vct_b(jk) = (REAL(jk-1,wp)/REAL(nlev,wp))**1.25_wp
    ENDDO

    z_exp = LOG(min_lay_thckn/top_height)/LOG(2._wp/pi*ACOS(REAL(nlev-1,wp)**stretch_fac/&
      &     REAL(nlev,wp)**stretch_fac))

    ! Set up distribution of coordinate surfaces according to the analytical formula
    ! vct = h_top*(2/pi*arccos(jk-1/nlev))**z_exp (taken from the COSMO model, src_artifdata)
    ! z_exp has been calculated above in order to return min_lay_thckn as thickness
    ! of the lowest model layer
    DO jk = 1, nlevp1
      vct_a(jk)      = top_height*(2._wp/pi*ACOS(REAL(jk-1,wp)**stretch_fac/ &
        &              REAL(nlev,wp)**stretch_fac))**z_exp
      vct(jk)        = vct_a(jk)
      vct(jk+nlevp1) = vct_b(jk)
    ENDDO

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

    IF (msg_level >= 15) THEN
     WRITE(message_text,'(a)') 'Heights of coordinate half levels (m):'
        CALL message('', TRIM(message_text))

      DO jk = 1, nlevp1
       WRITE(message_text,'(a,i4,F12.3)') 'jk, vct_a: ',jk, vct_a(jk)
        CALL message('', TRIM(message_text))
      ENDDO
    ENDIF

  END SUBROUTINE init_sleve_coord

  !---------------------------------------------------------------------------
  !>
  !! Computes the smoothed topography needed for the SLEVE coordinate.
  !! May be bypassed once an option for reading the smooth topography from data
  !! is available
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD, (2010-07-21)
  !!
  SUBROUTINE compute_smooth_topo(p_patch, p_int, topo_c, topo_smt_c, &
                                 topo_v, topo_smt_v)

    TYPE(t_patch),TARGET,INTENT(IN) :: p_patch
    TYPE(t_int_state), INTENT(IN) :: p_int

    ! Input fields: topography on cells and vertices
    REAL(wp), INTENT(IN) :: topo_c(:,:)
    REAL(wp), INTENT(IN) :: topo_v(:,:)

    ! Output fields: smooth topography on cells and vertices
    REAL(wp), INTENT(OUT) :: topo_smt_c(:,:)
    REAL(wp), INTENT(OUT) :: topo_smt_v(:,:)

    INTEGER  :: jb, jc, iter, niter
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),nabla2_topo(nproma,1,p_patch%nblks_c)
    REAL(wp) :: z_topo_v(nproma,1,p_patch%nblks_v)

    !-------------------------------------------------------------------------

    niter = 25 ! 25 smoothing iterations (do we need this to be a namelist variable?)

    ! Initialize auxiliary fields for topography with data and nullify nabla2 field
    z_topo(:,1,:)      = topo_c(:,:)
    z_topo_v(:,1,:)    = topo_v(:,:)
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

    ! Interpolate smooth topography from cells to vertices
    CALL cells2verts_scalar(z_topo,p_patch,p_int%cells_aw_verts,z_topo_v,1,1)

    CALL sync_patch_array(SYNC_V,p_patch,z_topo_v)

    ! Store smooth topography on output fields
    topo_smt_c(:,:) = z_topo(:,1,:)
    topo_smt_v(:,:) = z_topo_v(:,1,:)

  END SUBROUTINE compute_smooth_topo


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
    IF(my_process_is_mpi_parallel()) THEN
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

    INTEGER  :: jgc, jb, jc, jv, nblks_c
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom, i_rlstart, i_rlend

    INTEGER, POINTER, DIMENSION(:,:,:) :: iidx, iblk

    LOGICAL :: l_parallel

    !-------------------------------------------------------------------------

    jgc = p_pp%child_id(i_chidx) ! child domain ID

!    IF (p_nprocs == 1 .OR. p_pe == p_test_pe) THEN
      IF(my_process_is_mpi_parallel()) THEN
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
  SUBROUTINE init_vert_coord(p_patch, topo, topo_smt, z3d_i, z3d_m,              &
                             nlev, nblks, npromz, nshift, nflat, l_half_lev_centr)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch

    ! Input parameters:
    INTEGER, INTENT(IN) :: nlev, nblks, & ! field dimensions
                           npromz,      & ! length of last block
                           nshift,      & ! shift parameter for vertical nesting
                           nflat          ! index below which model levels are flat

    ! Input fields: "Normal" and smooth topography
    REAL(wp),  INTENT(IN) :: topo    (nproma,nblks), &
                             topo_smt(nproma,nblks)

    LOGICAL, INTENT(IN) :: l_half_lev_centr ! center half levels between full levels if true
 
    ! Output fields: 3D coordinate fields at interface and main levels
    REAL(wp),  INTENT(OUT) :: z3d_i(nproma,nlev+1,nblks), &
                              z3d_m(nproma,nlev,nblks)

    INTEGER :: jc, jk, jk1, jb, nlen, nlevp1
    REAL(wp) :: z_fac1, z_fac2, z_topo_dev(nproma)
    !-------------------------------------------------------------------------

    nlevp1 = nlev+1

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jk1, z_fac1, z_fac2, z_topo_dev)
    DO jb = 1,nblks

      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
     ENDIF

     z3d_i(1:nlen,nlevp1,jb) = topo(1:nlen,jb)

     ! vertical interface height
     IF (ivctype == 1) THEN ! hybrid Gal-Chen coordinate
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
       ! Ensure that layers do not intersect; except for the surface layer, the interface level
       ! distance is limited to min_lay_thckn
       DO jk = nlev-1, 1, -1
         DO jc = 1, nlen
           z3d_i(jc,jk,jb) = MAX(z3d_i(jc,jk,jb),z3d_i(jc,jk+1,jb)+min_lay_thckn)
         ENDDO
       ENDDO
     ENDIF

     DO jk = 1, nlev
       ! geometric height of full levels
       z3d_m(1:nlen,jk,jb) = 0.5_wp*(z3d_i(1:nlen,jk,jb)+z3d_i(1:nlen,jk+1,jb))
     ENDDO

     IF (l_half_lev_centr) THEN
       ! redefine half levels to be in the center between full levels
       DO jk = 2, nlev
         z3d_i(1:nlen,jk,jb) = 0.5_wp*(z3d_m(1:nlen,jk-1,jb)+z3d_m(1:nlen,jk,jb))
       ENDDO
     ENDIF

   ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_vert_coord

END MODULE mo_nh_init_utils
