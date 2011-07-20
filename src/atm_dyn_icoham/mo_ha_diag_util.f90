!>
!! This module contains subroutines that update various prognostic
!! variables.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (MPI-M, 2009-11)
!! Separated from mo_ha_dynamics by Hui Wan (MPI-M, 2010-02-16)
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
MODULE mo_ha_diag_util

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: grav, vtmpc1
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ext_data,           ONLY: t_external_data
  USE mo_math_operators,     ONLY: grad_fd_norm, div, div_avg, rot_vertex
  USE mo_dynamics_config,    ONLY: idiv_method, lshallow_water
  USE mo_ha_dyn_config,      ONLY: ha_dyn_config
  USE mo_io_config,          ONLY: lwrite_omega, l_outputtime
  USE mo_parallel_config,  ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: nlev, nlevp1, iqv, iqc, iqi, iqr, iqs, iforcing
  USE mo_impl_constants,     ONLY: inwp, iecham, ildf_echam
  USE mo_icoham_dyn_types,   ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_interpolation,      ONLY: t_int_state, cells2verts_scalar,        &
                                   cells2edges_scalar, edges2cells_scalar, &
                                   verts2edges_scalar, edges2verts_scalar, &
                                   rbf_vec_interpol_cell, cell_avg,        &
                                   i_cori_method, sick_a, sick_o
  USE mo_eta_coord_diag,     ONLY: half_level_pressure, full_level_pressure, &
                                   auxhyb, geopot
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_sync,               ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_vertical_coord_table, ONLY: delpr, nplev, nplvp1

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  PUBLIC :: update_diag_state
  PUBLIC :: update_dyn_output
  PUBLIC :: update_pres, update_pres_delp_c, update_delp_e
  PUBLIC :: update_tempv_geopot, update_omega

CONTAINS

  !>
  !! Update the diagnostic varibles and calculate some intermediate quantities
  !! from the given diagnostic variables.
  !!
  !! @par Revision History
  !! Separated from subroutine "dyn" and re-written by Hui Wan
  !! (MPI-M, 2009-11-19)
  !!
  SUBROUTINE update_diag_state( pt_prog, pt_patch, pt_int_state, pt_ext_data, &
    &        pt_diag )

  IMPLICIT NONE

  TYPE(t_hydro_atm_prog),INTENT(IN)    :: pt_prog       !< the prognostic variables
  TYPE(t_patch),TARGET,    INTENT(IN)    :: pt_patch      !< grid/patch info.
  TYPE(t_int_state),TARGET,INTENT(IN)    :: pt_int_state  !< horizontal interpolation coeff.
  TYPE(t_external_data),   INTENT(IN)    :: pt_ext_data   !< external data
  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag       !< diagnostic variables

    ! Pressure-related quantities
    CALL update_pres( pt_prog, pt_patch, pt_int_state, pt_diag )

    ! Virtual temperature (hydrostatic model only) and geopotential
    CALL update_tempv_geopot( pt_prog, pt_patch, pt_ext_data, pt_diag )

    ! Relative vorticity at vertices

    SELECT CASE (pt_patch%cell_type)
    CASE (3)
      CALL rot_vertex (pt_prog%vn, pt_patch, pt_int_state, pt_diag%rel_vort)
      CALL sync_patch_array(SYNC_V, pt_patch, pt_diag%rel_vort)
    CASE (6)
      CALL rot_vertex (pt_prog%vn, pt_patch, pt_int_state, pt_diag%rel_vort)
      CALL verts2edges_scalar(pt_diag%rel_vort, pt_patch, pt_int_state%tria_aw_rhom, &
      &                       pt_diag%rel_vort_e)
      CALL sync_patch_array(SYNC_E, pt_patch, pt_diag%rel_vort_e)
      CALL edges2verts_scalar(pt_diag%rel_vort_e, pt_patch, pt_int_state%e_1o3_v, &                      
      &                       pt_diag%rel_vort)
      CALL sync_patch_array(SYNC_V, pt_patch, pt_diag%rel_vort)
    END SELECT

    ! Variables needed for output but not for time stepping

    IF (l_outputtime) THEN
      CALL update_dyn_output( pt_patch, pt_int_state, pt_prog, pt_diag )
    ENDIF

  END SUBROUTINE update_diag_state
  !----------------------

  !>
  !! Compute vertical velocity of the pressure coordinate (omega = dp/dt)
  !!
  SUBROUTINE update_omega( p_vn, p_delp_e, p_pres_mc,  &! in
  &                        p_patch, p_int_state,       &! in
  &                        p_omega                   )  ! inout

    REAL(wp),INTENT(in) :: p_vn     (:,:,:)
    REAL(wp),INTENT(in) :: p_delp_e (:,:,:)
    REAL(wp),INTENT(in) :: p_pres_mc(:,:,:)

    TYPE(t_patch),TARGET,INTENT(in) :: p_patch
    TYPE(t_int_state),INTENT(in) :: p_int_state

    REAL(wp),INTENT(inout) :: p_omega(:,:,:) ! (nproma,nlev,nblks_c)

    REAL(wp) :: z_aux_me  (SIZE(p_vn,1),SIZE(p_vn,2),SIZE(p_vn,3))
    REAL(wp) :: z_aux_mc  (SIZE(p_omega,1),SIZE(p_omega,2),SIZE(p_omega,3))
    REAL(wp) :: z_aux_mc2 (SIZE(p_omega,1),SIZE(p_omega,2),SIZE(p_omega,3))

    INTEGER :: jk, jkp
    INTEGER :: nblks_c, nblks_e, jb, jbs, is, ie

    !------------
    ! Dimension parameters

    nblks_e = p_patch%nblks_int_e
    nblks_c = p_patch%nblks_int_c

    !-------------------------------------------------
    ! Part 1: sum_{j=1}^k div(vn*delp_e)_j
    !-------------------------------------------------
    ! Horizontal mass flux on full levels

    jbs = p_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_e
      CALL get_indices_e(p_patch, jb,jbs,nblks_e, is,ie, 2)
      z_aux_me(is:ie,:,jb) = p_delp_e(is:ie,:,jb)*p_vn(is:ie,:,jb)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Divergence of mass flux

    SELECT CASE(idiv_method)

    CASE(1)

      CALL div( z_aux_me, p_patch, p_int_state, z_aux_mc, opt_rlstart=2 )

    CASE(2)

      CALL div_avg( z_aux_me, p_patch, p_int_state, p_int_state%c_bln_avg, &
                    z_aux_mc, opt_rlstart=2 )

    END SELECT

    ! Vertically integrated mass divergence

    jbs = p_patch%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk,jkp)
    DO jb = jbs,nblks_c
      CALL get_indices_c(p_patch, jb,jbs,nblks_c, is,ie, 2)

      p_omega(is:ie,1,jb) = z_aux_mc(is:ie,1,jb)
      DO jk = 1,nlev-1
        jkp = jk + 1
        p_omega(is:ie,jkp,jb) = p_omega(is:ie,jk,jb) + z_aux_mc(is:ie,jkp,jb)
     ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    !-------------------------------------------------
    ! Part 2: vn*grad(p)
    !-------------------------------------------------

    CALL grad_fd_norm( p_pres_mc, p_patch, z_aux_me, nplvp1 )

    jbs = p_patch%edges%start_blk(4,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk)
    DO jb = jbs,nblks_e
      CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, 4)

      DO jk = nplvp1, nlev
         z_aux_me(is:ie,jk,jb) = z_aux_me(is:ie,jk,jb) * p_vn (is:ie,jk,jb)
      ENDDO

      DO jk = 1, nplev
         z_aux_me(is:ie,jk,jb) = 0._wp  ! no pressure gradient on p-levels
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Interpolation from edges to cell centers

    CALL edges2cells_scalar( z_aux_me, p_patch, p_int_state%e_inn_c, &
                             z_aux_mc, opt_rlstart=3 )

    ! Add smoothing if desired (but only when using the original gradient
    ! operator)

    SELECT CASE (idiv_method)
    CASE(2)
!$OMP PARALLEL
!$OMP WORKSHARE
      z_aux_mc2 = z_aux_mc
!$OMP END WORKSHARE
!$OMP END PARALLEL
      CALL cell_avg( z_aux_mc2, p_patch, p_int_state%c_bln_avg, &
                     z_aux_mc, opt_rlstart=4)
    END SELECT

    !-------------------------------------------------
    ! Combine two parts to get omega
    !-------------------------------------------------

    jbs = p_patch%cells%start_blk(3,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_c
      CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, 3)
      p_omega(is:ie,:,jb) = - p_omega(is:ie,:,jb) + z_aux_mc(is:ie,:,jb)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array( SYNC_C, p_patch, p_omega )

  END SUBROUTINE update_omega


  !>
  !! Compute full- and half-level pressure values, and layerthickess
  !! in pressure coordinate.
  !!
  SUBROUTINE update_pres_delp_c( p_patch, p_pres_sfc,              &! in
                                 p_pres_mc, p_pres_ic, p_delp_mc   )! inout

    TYPE(t_patch),TARGET,INTENT(in) :: p_patch            !< grid info

    REAL(wp),INTENT(in)    :: p_pres_sfc(:,  :)  !< surface pressure
    REAL(wp),INTENT(inout) :: p_pres_mc (:,:,:)  !< full-level pressure
    REAL(wp),INTENT(inout) :: p_pres_ic (:,:,:)  !< half-level pressure
    REAL(wp),INTENT(inout) :: p_delp_mc (:,:,:)  !< $delta p$ - layer thickness

    INTEGER :: jb, nlen, nblks_c, npromz_c

    ! Dimension parameters

    nblks_c  = p_patch%nblks_int_c
    npromz_c = p_patch%npromz_int_c

    ! Diagnose cell based quantities: pressure on full and half levels,
    ! layer thickness, and some auxiliary variables

!$OMP PARALLEL
    IF (.NOT. lshallow_water) THEN
!$OMP DO PRIVATE(jb, nlen)
      DO jb = 1,nblks_c

        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF

        CALL half_level_pressure( p_pres_sfc(:,jb), nproma, nlen, &! in
                                  p_pres_ic(:,:,jb)               )! out

        CALL full_level_pressure( p_pres_ic(:,:,jb),nproma, nlen, &! in
                                  p_pres_mc(:,:,jb)               )! out

        p_delp_mc(1:nlen,:,jb) =   p_pres_ic(1:nlen,2:nlevp1,jb) &
                               & - p_pres_ic(1:nlen,1:nlev  ,jb)

      ENDDO
!$OMP END DO

    ELSE ! shallow water
!$OMP DO PRIVATE(jb, nlen)
      DO jb = 1,nblks_c

        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF

        ! in shallow water mode, the surface pressure takes the role of the
        ! height depth, which is the same as the thickness in hydrostatic mode
        p_delp_mc(1:nlen,1,jb) = p_pres_sfc(1:nlen,jb)

      ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL
  END SUBROUTINE update_pres_delp_c
  !------------------

  !>
  !!
  !!
  SUBROUTINE update_delp_e( p_patch, p_int_state, p_delp_c, &! in
                            p_delp_e                        )! inout

    TYPE(t_patch),TARGET,INTENT(in) :: p_patch      !< grid info
    TYPE(t_int_state),INTENT(in) :: p_int_state  !< horizontal interpolation coeff.

    REAL(wp),INTENT(inout) :: p_delp_c (:,:,:) !< layer thickness at cell centers
    REAL(wp),INTENT(inout) :: p_delp_e (:,:,:) !< layer thickness at edges

    INTEGER :: jk, jb, jbs, nblks_e, is, ie

    REAL(wp) :: z_tmp_v(nproma,nlev,p_patch%nblks_v)
    REAL(wp) :: z_tmp_e(nproma,nlev,p_patch%nblks_e)

    !---

    nblks_e = p_patch%nblks_int_e

    ! Edge-based layer thickness: pressure levels

    jbs = p_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk)
    DO jb = jbs,nblks_e
      CALL get_indices_e(p_patch, jb,jbs,nblks_e, is,ie, 2)

      DO jk = 1,nplev
         p_delp_e(is:ie,jk,jb) = delpr(jk)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Edge-based layer thickness: sigma and transition levels

    CALL cells2edges_scalar (p_delp_c, p_patch,   &! in
                             p_int_state%c_lin_e, &! in
                             p_delp_e,            &! out
                             nplvp1, nlev )        ! optional input

    IF (i_cori_method >= 2) THEN

      ! avoid SICK instability and conserve energy
      IF (p_test_run) z_tmp_v=0.0_wp
      CALL cells2verts_scalar(p_delp_c, p_patch, p_int_state%cells_aw_verts, &
                              z_tmp_v, nplvp1, nlev)
      CALL sync_patch_array(SYNC_V,p_patch,z_tmp_v)
      CALL verts2edges_scalar(z_tmp_v, p_patch, p_int_state%v_1o2_e, &
                              z_tmp_e, nplvp1, nlev)
      jbs = p_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk)
      DO jb = jbs,nblks_e
        CALL get_indices_e(p_patch,jb,jbs,nblks_e,is,ie,2)
        DO jk = nplvp1,nlev
          p_delp_e(is:ie,jk,jb) = sick_a*z_tmp_e(is:ie,jk,jb)+sick_o*p_delp_e(is:ie,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ENDIF

  END SUBROUTINE update_delp_e
  !------------------
  !>
  !! Compute full- and half-level pressure values, layerthickess
  !! in pressure coordinate, and some pressure-related auxiliary quantities.
  !!
  SUBROUTINE update_pres( p_prog, p_patch, p_int_state, &! in
                          p_diag )                       ! inout

    TYPE(t_hydro_atm_prog),INTENT(in)    :: p_prog      !< carrying pres_sfc
    TYPE(t_patch),TARGET,    INTENT(in)    :: p_patch     !< grid info
    TYPE(t_int_state),       INTENT(in)    :: p_int_state !< horizontal interpol. coeff
    TYPE(t_hydro_atm_diag),INTENT(inout) :: p_diag      !< carrying the output

    INTEGER :: jk, jb, jbs, is, ie
    INTEGER :: nlen, nblks_c, npromz_c, nblks_e

    REAL(wp) :: z_tmp_v(nproma,nlev,p_patch%nblks_v)
    REAL(wp) :: z_tmp_e(nproma,nlev,p_patch%nblks_e)

    ! Dimension parameters

    nblks_e  = p_patch%nblks_int_e
    nblks_c  = p_patch%nblks_int_c
    npromz_c = p_patch%npromz_int_c

    ! Diagnose cell based quantities: pressure on full and half levels,
    ! layer thickness, and some auxiliary variables

!$OMP PARALLEL
    IF (.NOT. lshallow_water) THEN
!$OMP DO PRIVATE(jb, nlen)
      DO jb = 1,nblks_c

        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF

        CALL half_level_pressure( p_prog%pres_sfc(:,jb), nproma, nlen, &! in
                                  p_diag%pres_ic(:,:,jb)              ) ! out

        CALL full_level_pressure( p_diag%pres_ic(:,:,jb),nproma, nlen, &! in
                                  p_diag%pres_mc(:,:,jb)              ) ! out

        CALL auxhyb( p_diag%pres_ic(:,:,jb), nproma, nlen,            &! in
                     p_diag%delp_c(:,:,jb), p_diag%rdelp_c(:,:,jb),   &! out
                     p_diag%lnp_ic(:,:,jb), p_diag%rdlnpr_c(:,:,jb),  &! out
                     p_diag%rdalpha_c(:,:,jb)                       )  ! out
      ENDDO
!$OMP END DO

    ELSE ! shallow water
!$OMP DO PRIVATE(jb, nlen)
      DO jb = 1,nblks_c

        IF (jb /= nblks_c) THEN
           nlen = nproma
        ELSE
           nlen = npromz_c
        ENDIF

        ! in shallow water mode, the surface pressure takes the role of the
        ! height depth, which is the same as the thickness in hydrostatic mode
        p_diag%delp_c(1:nlen,1,jb) = p_prog%pres_sfc(1:nlen,jb)

        ! reciprocal layer thickness
        p_diag%rdelp_c(1:nlen,1,jb) = 1._wp/p_diag%delp_c(1:nlen,1,jb)
      ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL

    ! Edge-based layer thickness: pressure levels
    jbs = p_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk)
    DO jb = jbs,nblks_e
      CALL get_indices_e(p_patch, jb,jbs,nblks_e, is,ie, 2)

      DO jk = 1,nplev
         p_diag%delp_e(is:ie,jk,jb) = delpr(jk)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Edge-based layer thickness: sigma and transition levels

    CALL cells2edges_scalar (p_diag%delp_c, p_patch,  &! in
                             p_int_state%c_lin_e,     &! in
                             p_diag%delp_e,           &! out
                             nplvp1, nlev )            ! optional input

    IF (i_cori_method >= 2) THEN

      ! avoid SICK instability and conserve energy
      IF (p_test_run) z_tmp_v=0.0_wp
      CALL cells2verts_scalar(p_diag%delp_c, p_patch, p_int_state%cells_aw_verts, &
                              z_tmp_v, nplvp1, nlev)
      CALL sync_patch_array(SYNC_V,p_patch,z_tmp_v)
      CALL verts2edges_scalar(z_tmp_v, p_patch, p_int_state%v_1o2_e, &
                              z_tmp_e, nplvp1, nlev)
      jbs = p_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk)
      DO jb = jbs,nblks_e
        CALL get_indices_e(p_patch,jb,jbs,nblks_e,is,ie,2)
        DO jk = nplvp1,nlev
          p_diag%delp_e(is:ie,jk,jb) = sick_a*z_tmp_e(is:ie,jk,jb)&
          &                           +sick_o*p_diag%delp_e(is:ie,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ENDIF

  END SUBROUTINE update_pres
  !------------------
  !>
  !! Compute virtual temperature and geopotential.
  !!
  SUBROUTINE update_tempv_geopot( p_prog, p_patch, p_ext_data, &! in
                                  p_diag,                      &! inout
                                  opt_lgeop_wrt_sfc          )  ! optional input

    TYPE(t_hydro_atm_prog),INTENT(in)    :: p_prog
    TYPE(t_patch),           INTENT(in)    :: p_patch
    TYPE(t_external_data),   INTENT(IN)    :: p_ext_data !< external data
    TYPE(t_hydro_atm_diag),INTENT(inout) :: p_diag

    LOGICAL,INTENT(IN),OPTIONAL :: opt_lgeop_wrt_sfc

    ! Local variables

    REAL(wp) :: z_gzs(nproma)      !< surface geopotential
    INTEGER  :: jb, jbs, is, ie
    INTEGER  :: nblks_c
    LOGICAL  :: lgeop_wrt_sfc

    ! check optional input: by default this subroutine diagnoses geopotential
    ! at full and half levels; If opt_lgeop_wrt_sfc is set to .TRUE.,
    ! the values are given as the deviation with respect to
    ! surface geopotential.

    IF ( PRESENT(opt_lgeop_wrt_sfc) ) THEN
       lgeop_wrt_sfc = opt_lgeop_wrt_sfc
    ELSE
       lgeop_wrt_sfc = .FALSE.
    END IF

    ! inquire dimension parameter

    nblks_c  = p_patch%nblks_int_c

    !-------------------------------------------------------
    ! Diagnose virtual temperature (hydrostatic model only)
    !-------------------------------------------------------
!$OMP PARALLEL  PRIVATE(jbs)
    IF (.NOT. lshallow_water) THEN

      IF (ha_dyn_config%ldry_dycore) THEN
!$OMP WORKSHARE
        p_diag% qx       (:,:,:) = 0._wp
        p_diag% virt_incr(:,:,:) = 0._wp
        p_diag% tempv    (:,:,:) = p_prog%temp(:,:,:)
!$OMP END WORKSHARE

      ELSE !moist atmosphere
        jbs = p_patch%cells%start_blk(2,1)
!$OMP DO PRIVATE(jb,is,ie)
        DO jb = jbs,nblks_c
          CALL get_indices_c(p_patch, jb,jbs,nblks_c, is,ie, 2)
          SELECT CASE (iforcing)
          CASE (inwp) !DWD's NWP physics has 4 hydrometeors
            p_diag% qx(is:ie,:,jb) =   p_prog%tracer(is:ie,:,jb,iqc) &! cloud water
                                   & + p_prog%tracer(is:ie,:,jb,iqi) &! cloud ice
                                   & + p_prog%tracer(is:ie,:,jb,iqr) &! rain
                                   & + p_prog%tracer(is:ie,:,jb,iqs)  ! snow

          CASE (iecham,ildf_echam) !ECHAM only considers cloud water and cloud ice
            p_diag% qx(is:ie,:,jb) =   p_prog%tracer(is:ie,:,jb,iqc) &! cloud water
                                   & + p_prog%tracer(is:ie,:,jb,iqi)  ! cloud ice
          END SELECT

          p_diag% virt_incr(is:ie,:,jb) =   vtmpc1*p_prog%tracer(is:ie,:,jb,iqv) &
                                        & - p_diag%qx(is:ie,:,jb)

          p_diag%     tempv(is:ie,:,jb) =  p_prog%temp(is:ie,:,jb)              &
                                        & *( 1._wp+p_diag%virt_incr(is:ie,:,jb) )
        ENDDO
!$OMP END DO
       ENDIF !dry dycore vs moist atmos.

    ENDIF !shallow water vs hydrostatic

    !-----------------------
    ! Diagnose geopotential
    !-----------------------

    IF (lshallow_water) THEN
    ! geopotential = height * gravity
    ! Note: In this version of the shallow water model the thickness
    ! is prognostic and the height is diagnosed. Formerly it was the
    ! other way round.

      jbs = p_patch%cells%start_blk(2,1)
!$OMP DO PRIVATE(jb,is,ie,z_gzs)
      DO jb = jbs,nblks_c
        CALL get_indices_c(p_patch, jb,jbs,nblks_c, is,ie, 2)

        z_gzs(is:ie) = grav*p_ext_data%atm%topography_c(is:ie,jb)

        p_diag%geo_ic(is:ie,nlevp1,jb) = z_gzs(is:ie)
        p_diag%geo_mc(is:ie,1     ,jb) = z_gzs(is:ie)                  &
        &                              + grav*p_prog%pres_sfc(is:ie,jb)
      ENDDO
!$OMP END DO

    ELSE !Integrate the hydrostatic equation

      jbs = p_patch%cells%start_blk(2,1)
!$OMP DO PRIVATE(jb,is,ie,z_gzs)
      DO jb = jbs,nblks_c
        CALL get_indices_c(p_patch, jb,jbs,nblks_c, is,ie, 2)

        IF (lgeop_wrt_sfc) THEN
          z_gzs(is:ie) = 0._wp
        ELSE
          z_gzs(is:ie) = grav*p_ext_data%atm%topography_c(is:ie,jb)
        ENDIF

        CALL geopot( p_diag%tempv(:,:,jb), p_diag%rdlnpr_c(:,:,jb), &! in
        &            p_diag%rdalpha_c(:,:,jb), z_gzs(:),            &! in
        &            nproma, is,ie,                                 &! in
        &            p_diag%geo_mc(:,:,jb), p_diag%geo_ic(:,:,jb) )  ! out
      ENDDO
!$OMP END DO

    ENDIF !shallow water vs hydrostatic
!$OMP END PARALLEL

  END SUBROUTINE update_tempv_geopot
  !-------------------------

  !>
  !! Update variables that are needed for output but not for time integration.
  !!
  SUBROUTINE update_dyn_output( p_patch, p_int_state,  &! in
  &                             p_prog,                &! in
  &                             p_diag               )  ! inout

    TYPE(t_patch),TARGET,INTENT(in) :: p_patch
    TYPE(t_int_state),TARGET,INTENT(in) :: p_int_state

    TYPE(t_hydro_atm_prog),INTENT(in)    :: p_prog
    TYPE(t_hydro_atm_diag),INTENT(inout) :: p_diag

    ! Diagnose divergence

    SELECT CASE(idiv_method)
    CASE(1)
      CALL div( p_prog%vn, p_patch, p_int_state, p_diag%div)
    CASE(2)
      CALL div_avg( p_prog%vn, p_patch, p_int_state, p_int_state%c_bln_avg, &
      &             p_diag%div )
    END SELECT

    ! Reconstruct zonal and meridional wind at cell centers

    SELECT CASE(p_patch%cell_type)
    CASE(3)
      CALL rbf_vec_interpol_cell( p_prog%vn,    &! IN: normal wind
      &                           p_patch,      &! patch
      &                           p_int_state,  &! interpolation state
      &                           p_diag%u,     &! reconstructed u wind
      &                           p_diag%v   )   ! reconstructed y wind
    CASE(6)
      CALL edges2cells_scalar( p_prog%vn,             &! IN: normal wind
      &                        p_patch,               &! patch
      &                        p_int_state%hex_east,  &! interpolation state
      &                        p_diag%u              ) ! reconstr. u wind
      CALL edges2cells_scalar( p_prog%vn,             &! IN: normal wind
      &                        p_patch,               &! patch
      &                        p_int_state%hex_north, &! interpolation state
      &                        p_diag%v              ) ! reconstr. v wind
    END SELECT !triangle vs hexagon

    ! Vertical velocity omega=dp/dt

    IF (lwrite_omega.AND.(.NOT.lshallow_water)) THEN

      CALL update_omega( p_prog%vn, p_diag%delp_e, p_diag%pres_mc,  &! in
      &                  p_patch, p_int_state,                      &! in
      &                  p_diag%wpres_mc                          )  ! inout

    ENDIF

  END SUBROUTINE update_dyn_output
  !----------

END MODULE mo_ha_diag_util
