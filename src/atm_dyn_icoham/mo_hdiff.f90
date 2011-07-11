!>
!!  This module contains parameters, initialization subroutines and.
!!
!!  This module contains parameters, initialization subroutines and
!!  subroutine(s) for applying the horizontal diffusion.
!!
!! @par Revision History
!!  Original version by Hui Wan, MPI-M (2007-09-15)
!!  Modification by Hui Wan, MPI-M (2007-11-15)
!!  - added the option of using different diffusion coefficients
!!    for momentum and temperature.
!!  Modification by Hui Wan, MPI-M (2008-02-04)
!!  - added sub-steps of horizontal diffusion in a single dynamical step.
!!  Modification by Hui Wan, MPI-M (2008-02-09)
!!  - added implicit time stepping for horizontal diffusion.
!!  Modification by Hui Wan, MPI-M (2008-04-05)
!!  - relative tolerance rtol renamed hdiff_rtol.
!!  Modification by Hui Wan, MPI-M (2008-07-01)
!!  - added diagnostics related to the tendencies caused by horizontal diffusion.
!!  Modification by A. Gassmann, MPI-M (2008-09-23)
!!  - lhdiff_verbose removed, the solver should give feedback anyway if something
!!    did not work properly
!!  Restructuing by A. Gassmann, MPI-M (2008-10)
!!  - according to Hui Wan's suggestion:
!!    deleted the implicit treatment of diffusion
!!  Modification by A. Gassmann, MPI-M (2009-03)
!!  - efolding time also for hexagonal model
!!  Modification by A. Gassmann, MPI-M (2010-06)
!!  - Smagorinski diffusion for hexagonal model (hdiff_order=3)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_hdiff

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_model_domain_import, ONLY: nroot
  USE mo_diffusion_nml,       ONLY: k2, k4
  USE mo_diffusion_config,    ONLY: diffusion_config  
  USE mo_ha_dyn_config,       ONLY: ha_dyn_config
  USE mo_parallel_configuration,  ONLY: nproma
  USE mo_run_nml,             ONLY: nlev
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_math_operators,      ONLY: nabla2_vec, nabla2_scalar, &
                                    nabla4_vec, nabla4_scalar
  USE mo_math_constants,      ONLY: pi
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_interpolation,       ONLY: t_int_state, rbf_vec_interpol_cell, &
                                    rbf_vec_interpol_vertex, edges2cells_scalar, &
                                    cells2verts_scalar, edges2verts_scalar, &
                                    verts2cells_scalar
  USE mo_grf_interpolation,   ONLY: denom_diffu_v, denom_diffu_t, grf_intmethod_c
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, sync_patch_array
  USE mo_physical_constants,  ONLY: cpd, re

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: hdiff_expl

 CONTAINS

!--------------------------------------------------------------------
!-------------------------------------------------------------------------
!

 !>
 !!               Apply horizontal diffusion to the prognostic variables.
 !!
 !!
 !! @par Revision History
 !!  Original version by Hui Wan, MPI-M (2007-09)
 !!  Renamed from <i>hdiff</i> to <i>hdiff_expl</i> by Hui Wan (MPI-Met, 2008-2-9)
 !!  Restructuring by Almut Gassmann, MPI-M (2008-10)
 !!
 SUBROUTINE hdiff_expl( pt_patch, pt_new, pt_diag, pt_int, dtime, k_jg )

   IMPLICIT NONE

     TYPE(t_patch), TARGET    :: pt_patch     ! input
     TYPE(t_hydro_atm_prog) :: pt_new       ! in and out
     TYPE(t_hydro_atm_diag) :: pt_diag      ! input
     TYPE(t_int_state),TARGET :: pt_int       ! input

     REAL(wp) :: dtime
     INTEGER  :: k_jg         ! grid ID

     REAL(wp), DIMENSION(nproma,nlev,pt_patch%nblks_e) :: z_edge_val, z_nabla2_e
     REAL(wp), DIMENSION(nproma,nlev,pt_patch%nblks_c) :: z_cell_val, z_nabla2_c
     REAL(wp), DIMENSION(nproma,nlev,pt_patch%nblks_c) :: aux_cell

     INTEGER :: jb, jk, je, jlev, jc, jv

     REAL(wp):: diff_multfac = 0.0_wp ! to avoid compiler warning
     REAL(wp):: min_diffu

     ! index ranges needed for grid refinement
     INTEGER :: i_startblk, i_startidx, i_endidx, i_endblk, i_nchdom

     ! start index levels for boundary diffusion
     INTEGER :: start_bdydiff_c, start_bdydiff_e

     ! diffusion coefficients for boundary diffusion
     REAL(wp) :: fac_bdydiff_c, fac_bdydiff_e

     ! For Smagorinsky diffusion
     REAL(wp), DIMENSION(nproma,nlev,pt_patch%nblks_c) :: kh_smag_c
     REAL(wp), DIMENSION(nproma,nlev,pt_patch%nblks_e) :: kh_smag_e
     REAL(wp), DIMENSION(nproma,nlev,pt_patch%nblks_v) :: kh_smag_v
     REAL(wp), DIMENSION(nproma,nlev,pt_patch%nblks_v) :: u_vert
     REAL(wp), DIMENSION(nproma,nlev,pt_patch%nblks_v) :: v_vert
     REAL(wp) :: dvn_cell, dvt_cell, dvn_vert, dvt_vert

     INTEGER, POINTER :: ih1i(:,:,:), ih2i(:,:,:), ih1b(:,:,:), ih2b(:,:,:)
     INTEGER, POINTER :: it1i(:,:,:), it2i(:,:,:), it1b(:,:,:), it2b(:,:,:)
     INTEGER, POINTER :: ici (:,:,:), ivi (:,:,:), icb (:,:,:), ivb (:,:,:)
     INTEGER, POINTER :: icei(:,:,:), ivei(:,:,:), iceb(:,:,:), iveb(:,:,:)
     REAL(wp)         :: z_delp_v      (nproma,nlev,pt_patch%nblks_v), &
     &                   z_shear_def_1 (nproma,nlev,pt_patch%nblks_e), &
     &                   z_shear_def_2 (nproma,nlev,pt_patch%nblks_e), &
     &                   z_strain_def_1(nproma,nlev,pt_patch%nblks_e), &
     &                   z_strain_def_2(nproma,nlev,pt_patch%nblks_e), &
     &                   z_fric_heat_c1(nproma,nlev,pt_patch%nblks_c), &
     &                   z_fric_heat_v (nproma,nlev,pt_patch%nblks_v), &
     &                   z_fric_heat_c (nproma,nlev,pt_patch%nblks_c), &
     &                   z_turb_flx_c1 (nproma,nlev,pt_patch%nblks_e), &
     &                   z_turb_flx_c2 (nproma,nlev,pt_patch%nblks_e), &
     &                   z_turb_flx_v1 (nproma,nlev,pt_patch%nblks_e), &
     &                   z_turb_flx_v2 (nproma,nlev,pt_patch%nblks_e)
     REAL(wp), PARAMETER :: z_smag_min = 1.0e-20_wp ! security value
     REAL(wp) :: zhelp, z_mean_area_edge
     REAL(wp) :: hdiff_tv_ratio

     LOGICAL :: ltheta_dyn

!--------------------------------------------------------------------
!
    i_nchdom   = MAX(1,pt_patch%n_childdom)
    ltheta_dyn = ha_dyn_config%ltheta_dyn

    ! Parameters for boundary diffusion (relevant for nested domains only)
    start_bdydiff_c = 3 ! refin_ctrl level at which diffusion starts
    start_bdydiff_e = 5

    ! Normalized diffusion coefficients for boundary diffusion
    fac_bdydiff_c = 1._wp/denom_diffu_t  ! for temperature diffusion
    fac_bdydiff_e = 1._wp/denom_diffu_v  ! for momentum diffusion

    ! get hdiff_tv_ratio
    hdiff_tv_ratio = diffusion_config(k_jg)%hdiff_tv_ratio

    ! Pointers to cell and vertex neighbors
    ici => pt_patch%edges%cell_idx
    icb => pt_patch%edges%cell_blk
    ivi => pt_patch%edges%vertex_idx
    ivb => pt_patch%edges%vertex_blk


    IF ((pt_patch%cell_type == 3) .AND.               &
      &  (diffusion_config(k_jg)%hdiff_order == 5)) THEN
       !
       ! Compute diffusion coefficient for Smagorinsky diffusion
       ! ** NOT implemented for cell_type = 6 because      **
       ! ** RBF interpolation is used to compute deformation **

       ! a) Calculation of diffusion coefficient

       ! a.1) RBF reconstruction of velocity at cells and vertices

       CALL rbf_vec_interpol_cell( pt_new%vn, pt_patch, pt_int, &
                                   pt_diag%u, pt_diag%v)

       CALL rbf_vec_interpol_vertex( pt_new%vn, pt_patch, pt_int, &
                                     u_vert, v_vert)

       i_startblk = pt_patch%edges%start_blk(grf_bdywidth_e,1)
       i_endblk   = pt_patch%edges%end_blk(min_rledge,i_nchdom)

       jlev = pt_patch%level
       diff_multfac = MIN(4._wp,0.01375_wp*REAL(nroot*2**jlev,wp))*dtime ! empirically determined scaling factor
       min_diffu  = k2(k_jg)/SQRT(3._wp)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,dvn_cell,dvt_cell,dvn_vert,dvt_vert)
         DO jb = i_startblk,i_endblk

           CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                              i_startidx, i_endidx, grf_bdywidth_e, min_rledge)

           ! a.2 )Computation of wind field deformation

#ifdef __SX__
!CDIR UNROLL=6
#endif
           DO jk = 1, nlev
             DO je = i_startidx, i_endidx

               dvn_cell = pt_diag%u(ici(je,jb,1),jk,icb(je,jb,1)) * &
                          pt_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                          pt_diag%v(ici(je,jb,1),jk,icb(je,jb,1)) * &
                          pt_patch%edges%primal_normal_cell(je,jb,1)%v2 - &
                         (pt_diag%u(ici(je,jb,2),jk,icb(je,jb,2)) * &
                          pt_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                          pt_diag%v(ici(je,jb,2),jk,icb(je,jb,2)) * &
                          pt_patch%edges%primal_normal_cell(je,jb,2)%v2)
               dvt_cell = pt_diag%u(ici(je,jb,1),jk,icb(je,jb,1)) * &
                          pt_patch%edges%dual_normal_cell(je,jb,1)%v1 + &
                          pt_diag%v(ici(je,jb,1),jk,icb(je,jb,1)) * &
                          pt_patch%edges%dual_normal_cell(je,jb,1)%v2 - &
                         (pt_diag%u(ici(je,jb,2),jk,icb(je,jb,2)) * &
                          pt_patch%edges%dual_normal_cell(je,jb,2)%v1 + &
                          pt_diag%v(ici(je,jb,2),jk,icb(je,jb,2)) * &
                          pt_patch%edges%dual_normal_cell(je,jb,2)%v2)

               dvn_vert = u_vert(ivi(je,jb,1),jk,ivb(je,jb,1)) * &
                          pt_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
                          v_vert(ivi(je,jb,1),jk,ivb(je,jb,1)) * &
                          pt_patch%edges%primal_normal_vert(je,jb,1)%v2 - &
                         (u_vert(ivi(je,jb,2),jk,ivb(je,jb,2)) * &
                          pt_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
                          v_vert(ivi(je,jb,2),jk,ivb(je,jb,2)) * &
                          pt_patch%edges%primal_normal_vert(je,jb,2)%v2)
               dvt_vert = u_vert(ivi(je,jb,1),jk,ivb(je,jb,1)) * &
                          pt_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
                          v_vert(ivi(je,jb,1),jk,ivb(je,jb,1)) * &
                          pt_patch%edges%dual_normal_vert(je,jb,1)%v2 - &
                         (u_vert(ivi(je,jb,2),jk,ivb(je,jb,2)) * &
                          pt_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
                          v_vert(ivi(je,jb,2),jk,ivb(je,jb,2)) * &
                          pt_patch%edges%dual_normal_vert(je,jb,2)%v2)

               kh_smag_e(je,jk,jb) = diff_multfac* SQRT(( dvn_cell *    &
                 pt_patch%edges%inv_dual_edge_length(je,jb) -           &
                 dvt_vert*pt_patch%edges%system_orientation(je,jb)*     &
                 pt_patch%edges%inv_primal_edge_length(je,jb) )**2 + (  &
                 dvn_vert*pt_patch%edges%system_orientation(je,jb)*     &
                 pt_patch%edges%inv_primal_edge_length(je,jb) +         &
                 dvt_cell*pt_patch%edges%inv_dual_edge_length(je,jb))**2 )
             ENDDO
           ENDDO

           ! Subtract fourth-order background diffusion coefficient (calculated separately)
           kh_smag_e(i_startidx:i_endidx,:,jb) = &
             MAX(0._wp,kh_smag_e(i_startidx:i_endidx,:,jb) - k4(k_jg))

         ENDDO
!$OMP END DO
!$OMP END PARALLEL

         ! a.3) Interpolate diffusion coefficient to cell midpoints
         CALL edges2cells_scalar( kh_smag_e, pt_patch, pt_int%e_bln_c_s,  &
           kh_smag_c, opt_rlstart=grf_bdywidth_c+1, opt_rlend=min_rlcell)

     ENDIF

    !--------------------------------------------------------------------
    ! compute derivatives (nabla2/nabla4) needed for diffusion
    !--------------------------------------------------------------------

     SELECT CASE (diffusion_config(k_jg)%hdiff_order)
     CASE(-1)
         CONTINUE
     CASE(2) ! 2nd order diffusion ------

       ! Note: later, the diffusion coefficient must be determined by physics

       !---------------------------
       ! velocity
       !
        IF (diffusion_config(k_jg)%lhdiff_vn) THEN

           CALL nabla2_vec( pt_new%vn, pt_patch, pt_int, z_edge_val, &
                           opt_rlstart=5, opt_rlend=min_rledge)

        ENDIF ! if lhdiff_vn

       !---------------------------
       ! temperature
       !---------------------------

        IF (diffusion_config(k_jg)%lhdiff_temp .AND. (.NOT.ltheta_dyn)) THEN

           CALL nabla2_scalar( pt_new%temp, pt_patch, pt_int, z_cell_val, &
                               opt_rlstart=3, opt_rlend=min_rlcell )

        ELSE IF (diffusion_config(k_jg)%lhdiff_temp .AND. ltheta_dyn) THEN

           ! Divide by delta p so as to diffuse the uncoupled theta
           i_startblk = pt_patch%cells%start_blk(1,1)
           i_endblk   = pt_patch%cells%end_blk(min_rlcell,i_nchdom)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk)
           DO jb = i_startblk,i_endblk

             CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, 1, min_rlcell)

             DO jk = 1, nlev
               ! Note: we use delp_c here because convert_t2theta/convert_theta2t
               ! update delp_c but not rdelp_c
               aux_cell(i_startidx:i_endidx,jk,jb) =       &
                 pt_new%theta(i_startidx:i_endidx,jk,jb) / &
                 pt_diag%delp_c(i_startidx:i_endidx,jk,jb)
             ENDDO
           ENDDO
!$OMP END DO
!$OMP END PARALLEL

           CALL nabla2_scalar( aux_cell, pt_patch, pt_int, z_cell_val, &
                               opt_rlstart=3, opt_rlend=min_rlcell )

        ENDIF

     CASE(4,5) ! 4th order diffusion ------

        !---------------------------
        ! velocity
        !---------------------------

        IF (diffusion_config(k_jg)%lhdiff_vn) THEN

           CALL nabla4_vec( pt_new%vn, pt_patch, pt_int, z_edge_val,  &
                            opt_nabla2=z_nabla2_e, opt_rlstart=7, &
                            opt_rlend=min_rledge )

        ENDIF ! if lhdiff_vn

        !---------------------------
        ! temperature
        !---------------------------

        IF (diffusion_config(k_jg)%lhdiff_temp .AND. (.NOT.ltheta_dyn)) THEN

           CALL nabla4_scalar( pt_new%temp, pt_patch, pt_int, z_cell_val, &
                               opt_nabla2=z_nabla2_c, opt_rlstart=4,     &
                               opt_rlend=min_rlcell )

        ELSE IF (diffusion_config(k_jg)%lhdiff_temp .AND. ltheta_dyn) THEN

           ! Divide by delta p so as to diffuse the uncoupled theta
           i_startblk = pt_patch%cells%start_blk(1,1)
           i_endblk   = pt_patch%cells%end_blk(min_rlcell,i_nchdom)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk)
           DO jb = i_startblk,i_endblk

             CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, 1, min_rlcell)

              DO jk = 1, nlev
                ! Note: we use delp_c here because convert_t2theta/convert_theta2t
                ! update delp_c but not rdelp_c
                aux_cell(i_startidx:i_endidx,jk,jb) =       &
                  pt_new%theta(i_startidx:i_endidx,jk,jb) / &
                  pt_diag%delp_c(i_startidx:i_endidx,jk,jb)
              ENDDO
           ENDDO
!$OMP END DO
!$OMP END PARALLEL

           CALL nabla4_scalar( aux_cell, pt_patch, pt_int, z_cell_val, &
                               opt_nabla2=z_nabla2_c, opt_rlstart=4,     &
                               opt_rlend=min_rlcell )
        ENDIF
     END SELECT

    !--------------------------------------------------------------------
    ! apply diffusion
    !--------------------------------------------------------------------

     SELECT CASE (diffusion_config(k_jg)%hdiff_order)
     CASE(-1)
         CONTINUE
     CASE(2) ! 2nd order diffusion ------

       ! Note: later, the diffusion coefficient must be determined by physics

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,i_startidx,i_endidx,diff_multfac)

       !---------------------------
       ! velocity
       !
        IF (diffusion_config(k_jg)%lhdiff_vn) THEN

           ! Diffusion in the domain interior (boundary diffusion follows below)
           i_startblk = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
           i_endblk   = pt_patch%edges%end_blk(min_rledge,i_nchdom)

           diff_multfac = 1._wp/SQRT(3._wp)*k2(k_jg)

!$OMP DO PRIVATE(jb,jk)
           DO jb = i_startblk,i_endblk

             CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, grf_bdywidth_e+1, min_rledge)

             DO jk = 1, nlev
               pt_new%vn(i_startidx:i_endidx,jk,jb) =   &
                 pt_new%vn(i_startidx:i_endidx,jk,jb) + &
                 z_edge_val(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%edges%area_edge(i_startidx:i_endidx,jb)*diff_multfac
             ENDDO
           ENDDO
!$OMP END DO

        ENDIF

        IF (diffusion_config(k_jg)%lhdiff_vn .AND. pt_patch%id > 1) THEN

           ! Lateral boundary diffusion
           i_startblk = pt_patch%edges%start_blk(start_bdydiff_e,1)
           i_endblk   = pt_patch%edges%end_blk(grf_bdywidth_e,1)

           DO jb = i_startblk,i_endblk

             CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, start_bdydiff_e, grf_bdywidth_e)

! OpenMP parallelization is done over jk for boundary diffusion because it affects
! only few blocks
!$OMP DO PRIVATE(jk)
             DO jk = 1, nlev
               pt_new%vn(i_startidx:i_endidx,jk,jb) =   &
                 pt_new%vn(i_startidx:i_endidx,jk,jb) + &
                 z_edge_val(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%edges%area_edge(i_startidx:i_endidx,jb)*fac_bdydiff_e
             ENDDO
!$OMP END DO
           ENDDO

        ENDIF ! vn boundary diffusion

        !---------------------------
        ! temperature
        !---------------------------
        IF (diffusion_config(k_jg)%lhdiff_temp) THEN

          SELECT CASE (pt_patch%cell_type)
          CASE (6)
            ! old
            !diff_multfac = 2._wp/3._wp*SQRT(3._wp)*hdiff_tv_ratio*k2(k_jg)
            diff_multfac = 2._wp/9._wp*hdiff_tv_ratio*k2(k_jg)/SQRT(3._wp)
            ! rejected new
            !!diff_multfac = 2._wp/9._wp*hdiff_tv_ratio*k2(k_jg)
          CASE (3)
            ! old
            diff_multfac = 4._wp/3._wp/SQRT(3._wp)*hdiff_tv_ratio*k2(k_jg)
            ! rejected new
            !!diff_multfac = 4._wp/3._wp*hdiff_tv_ratio*k2(k_jg)
          END SELECT

          IF (.NOT.ltheta_dyn) THEN

            ! Diffusion in the domain interior (boundary diffusion follows below)
            i_startblk = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
            i_endblk   = pt_patch%cells%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,jk)
            DO jb = i_startblk,i_endblk

              CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                 i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

              DO jk = 1, nlev
                pt_new%temp(i_startidx:i_endidx,jk,jb) =   &
                & pt_new%temp(i_startidx:i_endidx,jk,jb) + &
                & z_cell_val(i_startidx:i_endidx,jk,jb) * &
                & pt_patch%cells%area(i_startidx:i_endidx,jb)*diff_multfac
              ENDDO
            ENDDO
!$OMP END DO


          ELSEIF (ltheta_dyn) THEN

            ! Diffusion in the domain interior (boundary diffusion follows below)
            i_startblk = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
            i_endblk   = pt_patch%cells%end_blk(min_rlcell,i_nchdom)

!$OMP DO PRIVATE(jb,jk)
            DO jb = i_startblk,i_endblk

              CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)
              DO jk = 1, nlev
                pt_new%theta(i_startidx:i_endidx,jk,jb) =   &
                  pt_new%theta(i_startidx:i_endidx,jk,jb) + &
                 z_cell_val(i_startidx:i_endidx,jk,jb)*   &
                  pt_patch%cells%area(i_startidx:i_endidx,jb)*diff_multfac* &
                  pt_diag%delp_c(i_startidx:i_endidx,jk,jb)
              ENDDO
            ENDDO
!$OMP END DO

          ENDIF ! if ltheta_dyn
        ENDIF ! lhdiff_temp

        IF (grf_intmethod_c == 1 .AND. diffusion_config(k_jg)%lhdiff_temp &
          & .AND. .NOT.ltheta_dyn .AND. pt_patch%id > 1) THEN

           ! Lateral boundary diffusion
           i_startblk = pt_patch%cells%start_blk(start_bdydiff_c,1)
           i_endblk   = pt_patch%cells%end_blk(grf_bdywidth_c,1)

           DO jb = i_startblk,i_endblk

             CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, start_bdydiff_c, grf_bdywidth_c)

! OpenMP parallelization is done over jk for boundary diffusion because it affects
! only few blocks
!$OMP DO PRIVATE(jk)
             DO jk = 1, nlev
               pt_new%temp(i_startidx:i_endidx,jk,jb) =   &
                 pt_new%temp(i_startidx:i_endidx,jk,jb) + &
                 z_cell_val(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%cells%area(i_startidx:i_endidx,jb)*fac_bdydiff_c
             ENDDO
!$OMP END DO
           ENDDO

        ELSE IF (grf_intmethod_c == 1 .AND. diffusion_config(k_jg)%lhdiff_temp &
          &  .AND. ltheta_dyn .AND. pt_patch%id > 1) THEN

           ! Lateral boundary diffusion
           i_startblk = pt_patch%cells%start_blk(start_bdydiff_c,1)
           i_endblk   = pt_patch%cells%end_blk(grf_bdywidth_c,1)

           DO jb = i_startblk,i_endblk

             CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, start_bdydiff_c, grf_bdywidth_c)

! OpenMP parallelization is done over jk for boundary diffusion because it affects
! only few blocks
!$OMP DO PRIVATE(jk)
             DO jk = 1, nlev
               pt_new%theta(i_startidx:i_endidx,jk,jb) =   &
                 pt_new%theta(i_startidx:i_endidx,jk,jb) + &
                 z_cell_val(i_startidx:i_endidx,jk,jb)*   &
                 pt_patch%cells%area(i_startidx:i_endidx,jb)* &
                 pt_diag%delp_c(i_startidx:i_endidx,jk,jb)*fac_bdydiff_c
             ENDDO
!$OMP END DO
           ENDDO

        ENDIF ! temperature boundary diffusion

!$OMP END PARALLEL

     CASE(3) ! Smagorinski diffusion for hexagonal model

       IF (diffusion_config(k_jg)%lhdiff_vn) THEN

         ! mean area edge
         z_mean_area_edge=8.0_wp*pi*re*re/pt_patch%n_patch_edges_g

         ! c) compute thicknesses at vertices
         CALL cells2verts_scalar(pt_diag%delp_c, pt_patch, pt_int%cells_aw_verts, z_delp_v)

         it1i => pt_int%dir_gradt_i1
         it2i => pt_int%dir_gradt_i2
         it1b => pt_int%dir_gradt_b1
         it2b => pt_int%dir_gradt_b2
         ih1i => pt_int%dir_gradh_i1
         ih2i => pt_int%dir_gradh_i2
         ih1b => pt_int%dir_gradh_b1
         ih2b => pt_int%dir_gradh_b2
         icei => pt_patch%cells%edge_idx
         iceb => pt_patch%cells%edge_blk
         ivei => pt_patch%verts%edge_idx
         iveb => pt_patch%verts%edge_blk

!$OMP PARALLEL

!$OMP WORKSHARE
         z_turb_flx_v1(:,:,:) = 0.0_wp
         z_turb_flx_v2(:,:,:) = 0.0_wp
         kh_smag_e(:,:,:)     = 0.0_wp
!$OMP END WORKSHARE

         i_startblk = pt_patch%edges%start_blk(3,1)
         i_endblk   = pt_patch%edges%end_blk(min_rledge,1)
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
         DO jb = i_startblk, i_endblk
           CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
           &                  i_startidx, i_endidx, 3, min_rledge)
           DO jk = 1, nlev
             DO je = i_startidx, i_endidx
               ! d) Shear deformation at vertices
               z_shear_def_1(je,jk,jb) = &
               &(pt_int%shear_def_v1(1,je,jb)*pt_new%vn(it1i(1,je,jb),jk,it1b(1,je,jb)) &
               &+pt_int%shear_def_v1(2,je,jb)*pt_new%vn(it1i(2,je,jb),jk,it1b(2,je,jb)) &
               &+pt_int%shear_def_v1(3,je,jb)*pt_new%vn(it1i(3,je,jb),jk,it1b(3,je,jb)) &
               &+pt_int%shear_def_v1(4,je,jb)*pt_new%vn(it1i(4,je,jb),jk,it1b(4,je,jb)) &
               &+pt_int%shear_def_v1(5,je,jb)*pt_new%vn(it1i(5,je,jb),jk,it1b(5,je,jb)) &
               &+pt_int%shear_def_v1(6,je,jb)*pt_new%vn(it1i(6,je,jb),jk,it1b(6,je,jb)) &
               &+pt_int%shear_def_v1(7,je,jb)*pt_new%vn(it1i(7,je,jb),jk,it1b(7,je,jb)) &
               &+pt_int%shear_def_v1(8,je,jb)*pt_new%vn(it1i(8,je,jb),jk,it1b(8,je,jb)) &
               &+pt_int%shear_def_v1(9,je,jb)*pt_new%vn(it1i(9,je,jb),jk,it1b(9,je,jb)) )
               z_shear_def_2(je,jk,jb) =  &
               &(pt_int%shear_def_v2(1,je,jb)*pt_new%vn(it2i(1,je,jb),jk,it2b(1,je,jb)) &
               &+pt_int%shear_def_v2(2,je,jb)*pt_new%vn(it2i(2,je,jb),jk,it2b(2,je,jb)) &
               &+pt_int%shear_def_v2(3,je,jb)*pt_new%vn(it2i(3,je,jb),jk,it2b(3,je,jb)) &
               &+pt_int%shear_def_v2(4,je,jb)*pt_new%vn(it2i(4,je,jb),jk,it2b(4,je,jb)) &
               &+pt_int%shear_def_v2(5,je,jb)*pt_new%vn(it2i(5,je,jb),jk,it2b(5,je,jb)) &
               &+pt_int%shear_def_v2(6,je,jb)*pt_new%vn(it2i(6,je,jb),jk,it2b(6,je,jb)) &
               &+pt_int%shear_def_v2(7,je,jb)*pt_new%vn(it2i(7,je,jb),jk,it2b(7,je,jb)) &
               &+pt_int%shear_def_v2(8,je,jb)*pt_new%vn(it2i(8,je,jb),jk,it2b(8,je,jb)) &
               &+pt_int%shear_def_v2(9,je,jb)*pt_new%vn(it2i(9,je,jb),jk,it2b(9,je,jb)) )
               ! e) Strain deformation at centers
               z_strain_def_1(je,jk,jb) = &
               &(pt_int%strain_def_c1(1,je,jb)*pt_new%vn(ih1i(1,je,jb),jk,ih1b(1,je,jb)) &
               &+pt_int%strain_def_c1(2,je,jb)*pt_new%vn(ih1i(2,je,jb),jk,ih1b(2,je,jb)) &
               &+pt_int%strain_def_c1(3,je,jb)*pt_new%vn(ih1i(3,je,jb),jk,ih1b(3,je,jb)) &
               &+pt_int%strain_def_c1(4,je,jb)*pt_new%vn(ih1i(4,je,jb),jk,ih1b(4,je,jb)) &
               &+pt_int%strain_def_c1(5,je,jb)*pt_new%vn(ih1i(5,je,jb),jk,ih1b(5,je,jb)) &
               &+pt_int%strain_def_c1(6,je,jb)*pt_new%vn(ih1i(6,je,jb),jk,ih1b(6,je,jb)) )
               z_strain_def_2(je,jk,jb) = &
               &(pt_int%strain_def_c2(1,je,jb)*pt_new%vn(ih2i(1,je,jb),jk,ih2b(1,je,jb)) &
               &+pt_int%strain_def_c2(2,je,jb)*pt_new%vn(ih2i(2,je,jb),jk,ih2b(2,je,jb)) &
               &+pt_int%strain_def_c2(3,je,jb)*pt_new%vn(ih2i(3,je,jb),jk,ih2b(3,je,jb)) &
               &+pt_int%strain_def_c2(4,je,jb)*pt_new%vn(ih2i(4,je,jb),jk,ih2b(4,je,jb)) &
               &+pt_int%strain_def_c2(5,je,jb)*pt_new%vn(ih2i(5,je,jb),jk,ih2b(5,je,jb)) &
               &+pt_int%strain_def_c2(6,je,jb)*pt_new%vn(ih2i(6,je,jb),jk,ih2b(6,je,jb)) )
               ! f) Diffusion coefficient at edge
               kh_smag_e(je,jk,jb) = diffusion_config(k_jg)%hdiff_smag_fac               &
               & * pt_patch%edges%area_edge(je,jb) * SQRT(z_smag_min                     &
               & + 0.5_wp*(z_strain_def_1(je,jk,jb)**2+z_strain_def_2(je,jk,jb)**2)      &
               & +(pt_patch%edges%edge_vert_length(je,jb,1)*(z_shear_def_1(je,jk,jb)**2) &
               &  +pt_patch%edges%edge_vert_length(je,jb,2)*(z_shear_def_2(je,jk,jb)**2))&
               &  /pt_patch%edges%primal_edge_length(je,jb) )
             ENDDO
           ENDDO
         ENDDO
!$OMP END DO
!$OMP END PARALLEL

         CALL sync_patch_array(SYNC_E, pt_patch, kh_smag_e)
         ! g) average the diffusion coefficient to the centers and vertices
         CALL edges2cells_scalar(kh_smag_e, pt_patch, pt_int%e_aw_c, kh_smag_c)
         CALL edges2verts_scalar(kh_smag_e, pt_patch, pt_int%e_aw_v, kh_smag_v)

!$OMP PARALLEL
         i_startblk = pt_patch%edges%start_blk(3,1)
         i_endblk   = pt_patch%edges%end_blk(min_rledge,1)
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
         DO jb = i_startblk,i_endblk
           CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
           &                  i_startidx, i_endidx, 3, min_rledge)
           DO jk = 1, nlev
             DO je = i_startidx, i_endidx
               ! h) turbulent fluxes
               z_turb_flx_c1(je,jk,jb) = &
               &    -pt_diag%delp_c(ici(je,jb,1),jk,icb(je,jb,1)) &
               &    *kh_smag_c(ici(je,jb,1),jk,icb(je,jb,1))      &
               &    *z_strain_def_1(je,jk,jb)
               z_turb_flx_c2(je,jk,jb) = &
               &    -pt_diag%delp_c(ici(je,jb,2),jk,icb(je,jb,2)) &
               &    *kh_smag_c(ici(je,jb,2),jk,icb(je,jb,2))      &
               &    *z_strain_def_2(je,jk,jb)
               z_turb_flx_v1(je,jk,jb) = &
               &    -z_delp_v(ivi(je,jb,1),jk,ivb(je,jb,1))       &
               &    *kh_smag_v(ivi(je,jb,1),jk,ivb(je,jb,1))      &
               &    *z_shear_def_1(je,jk,jb)
               z_turb_flx_v2(je,jk,jb) = &
               &    -z_delp_v(ivi(je,jb,2),jk,ivb(je,jb,2))       &
               &    *kh_smag_v(ivi(je,jb,2),jk,ivb(je,jb,2))      &
               &    *z_shear_def_2(je,jk,jb)
             ENDDO
           ENDDO
         ENDDO
!$OMP END DO
!$OMP END PARALLEL

         IF (diffusion_config(k_jg)%lhdiff_temp) THEN
           CALL sync_patch_array(SYNC_E, pt_patch, z_turb_flx_v1)
           CALL sync_patch_array(SYNC_E, pt_patch, z_turb_flx_v2)

           i_startblk = pt_patch%cells%start_blk(2,1)
           i_endblk   = pt_patch%cells%end_blk(min_rlcell,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,je,zhelp,i_startidx,i_endidx)
           DO jb = i_startblk,i_endblk
           CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
           &                  i_startidx, i_endidx, 2, min_rlcell)
             z_fric_heat_c(:,:,jb)=0.0_wp
             DO je = 1, 6
               DO jk = 1, nlev
                 DO jc = i_startidx,i_endidx

                   IF (je > pt_patch%cells%num_edges(jc,jb)) CYCLE

                   zhelp = pt_patch%edges%system_orientation(icei(jc,jb,je),iceb(jc,jb,je)) &
                   &      *pt_patch%cells%edge_orientation(jc,jb,je)

                   ! i) frictional heating at centers
                   z_fric_heat_c(jc,jk,jb) = z_fric_heat_c(jc,jk,jb)                       &
                   &-0.5_wp*((zhelp+1.0_wp)*z_turb_flx_c1(icei(jc,jb,je),jk,iceb(jc,jb,je)) &
                   &        -(zhelp-1.0_wp)*z_turb_flx_c2(icei(jc,jb,je),jk,iceb(jc,jb,je)))&
                   &*pt_new%vn(icei(jc,jb,je),jk,iceb(jc,jb,je))*pt_int%geofac_div(jc,je,jb)

                 ENDDO
               ENDDO
             ENDDO
           ENDDO
!$OMP END DO
           i_startblk = pt_patch%verts%start_blk(2,1)
           i_endblk   = pt_patch%verts%end_blk(min_rlvert,1)
!$OMP DO PRIVATE(jb,jk,jv,i_startidx,i_endidx)
           DO jb = i_startblk,i_endblk
           CALL get_indices_v(pt_patch, jb, i_startblk, i_endblk, &
           &                  i_startidx, i_endidx, 2, min_rlvert)
             DO jk = 1, nlev
               DO jv = i_startidx, i_endidx

                 ! j) frictional heating at vertices
                 z_fric_heat_v(jv,jk,jb) = 0.5_wp* (&
                 & ((pt_patch%verts%edge_orientation(jv,jb,1)+1.0_wp) &
                 &  *z_turb_flx_v1(ivei(jv,jb,1),jk,iveb(jv,jb,1))    &
                 & -(pt_patch%verts%edge_orientation(jv,jb,1)-1.0_wp) &
                 &  *z_turb_flx_v2(ivei(jv,jb,1),jk,iveb(jv,jb,1)))   &
                 &*pt_new%vn(ivei(jv,jb,1),jk,iveb(jv,jb,1))*pt_int%geofac_rot(jv,1,jb) &
                 &+((pt_patch%verts%edge_orientation(jv,jb,2)+1.0_wp) &
                 &  *z_turb_flx_v1(ivei(jv,jb,2),jk,iveb(jv,jb,2))    &
                 & -(pt_patch%verts%edge_orientation(jv,jb,2)-1.0_wp) &
                 &  *z_turb_flx_v2(ivei(jv,jb,2),jk,iveb(jv,jb,2)))   &
                 &*pt_new%vn(ivei(jv,jb,2),jk,iveb(jv,jb,2))*pt_int%geofac_rot(jv,2,jb) &
                 &+((pt_patch%verts%edge_orientation(jv,jb,3)+1.0_wp) &
                 &  *z_turb_flx_v1(ivei(jv,jb,3),jk,iveb(jv,jb,3))    &
                 & -(pt_patch%verts%edge_orientation(jv,jb,3)-1.0_wp) &
                 &  *z_turb_flx_v2(ivei(jv,jb,3),jk,iveb(jv,jb,3)))   &
                 &*pt_new%vn(ivei(jv,jb,3),jk,iveb(jv,jb,3))*pt_int%geofac_rot(jv,3,jb))

               ENDDO
             ENDDO
           ENDDO
!$OMP END DO
!$OMP END PARALLEL

         ENDIF

         i_startblk = pt_patch%edges%start_blk(3,1)
         i_endblk   = pt_patch%edges%end_blk(min_rledge,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
         DO jb = i_startblk,i_endblk
           CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
           &                  i_startidx, i_endidx, 3, min_rledge)
           DO jk = 1, nlev
             DO je = i_startidx,i_endidx
               ! k) Tendencies to the velocity components
               pt_new%vn(je,jk,jb) = pt_new%vn(je,jk,jb) - dtime  &
               & *((z_turb_flx_c2(je,jk,jb)- z_turb_flx_c1(je,jk,jb)) &
               &   *pt_patch%edges%inv_dual_edge_length(je,jb)    &
               &   *pt_patch%edges%system_orientation(je,jb)      &
               &  +(z_turb_flx_v2(je,jk,jb)- z_turb_flx_v1(je,jk,jb)) &
               &   *pt_patch%edges%inv_primal_edge_length(je,jb)  &
               &  )/pt_diag%delp_e(je,jk,jb)
             ENDDO
           ENDDO
         ENDDO
!$OMP END DO
!$OMP END PARALLEL

         IF(diffusion_config(k_jg)%lhdiff_temp) THEN
           ! l) frictional heating at cells averaged from vertices
           CALL verts2cells_scalar(z_fric_heat_v, pt_patch, pt_int%verts_aw_cells, z_fric_heat_c1)
           IF (.NOT.ltheta_dyn) THEN
             i_startblk = pt_patch%cells%start_blk(2,1)
             i_endblk   = pt_patch%cells%end_blk(min_rlcell,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx)
             DO jb = i_startblk,i_endblk
               CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
               &                  i_startidx, i_endidx, 2, min_rlcell)
               DO jk = 1, nlev
                 ! k) update new temperature
                 pt_new%temp(i_startidx:i_endidx,jk,jb) = &
                 & pt_new%temp(i_startidx:i_endidx,jk,jb) &
                 & + dtime/cpd*(z_fric_heat_c(i_startidx:i_endidx,jk,jb)  &
                 &             +z_fric_heat_c1(i_startidx:i_endidx,jk,jb))&
                 &             /pt_diag%delp_c(i_startidx:i_endidx,jk,jb)
               ENDDO
             ENDDO
!$OMP END DO
!$OMP END PARALLEL
           ENDIF
         ENDIF
       ENDIF

     CASE(4) ! 4th order diffusion ------

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,i_startidx,i_endidx,diff_multfac)

        !---------------------------
        ! velocity
        !---------------------------

        IF (diffusion_config(k_jg)%lhdiff_vn) THEN

           ! Diffusion in the domain interior (boundary diffusion follows below)
           i_startblk = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
           i_endblk   = pt_patch%edges%end_blk(min_rledge,i_nchdom)

           diff_multfac = k4(k_jg)/3._wp

!$OMP DO PRIVATE(jb,jk)
           DO jb = i_startblk,i_endblk

             CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, grf_bdywidth_e+1, min_rledge)

             DO jk = 1, nlev
               pt_new%vn(i_startidx:i_endidx,jk,jb) =       &
                 pt_new%vn(i_startidx:i_endidx,jk,jb)  -    &
                 z_edge_val(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%edges%area_edge(i_startidx:i_endidx,jb)* &
                 pt_patch%edges%area_edge(i_startidx:i_endidx,jb)*diff_multfac
             ENDDO
           ENDDO
!$OMP END DO

        ENDIF ! if lhdiff_vn

        IF (diffusion_config(k_jg)%lhdiff_vn .AND. pt_patch%id > 1) THEN

           ! Lateral boundary diffusion
           i_startblk = pt_patch%edges%start_blk(start_bdydiff_e,1)
           i_endblk   = pt_patch%edges%end_blk(grf_bdywidth_e,1)

           DO jb = i_startblk,i_endblk

             CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, start_bdydiff_e, grf_bdywidth_e)

! OpenMP parallelization is done over jk for boundary diffusion because it affects
! only few blocks
!$OMP DO PRIVATE(jk)
             DO jk = 1, nlev
               pt_new%vn(i_startidx:i_endidx,jk,jb) =   &
                 pt_new%vn(i_startidx:i_endidx,jk,jb) + &
                 z_nabla2_e(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%edges%area_edge(i_startidx:i_endidx,jb)*fac_bdydiff_e
             ENDDO
!$OMP END DO
           ENDDO

        ENDIF ! vn boundary diffusion

        !---------------------------
        ! temperature
        !---------------------------

        IF (diffusion_config(k_jg)%lhdiff_temp ) THEN


          ! Diffusion in the domain interior (boundary diffusion follows below)
          i_startblk = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
          i_endblk   = pt_patch%cells%end_blk(min_rlcell,i_nchdom)

          SELECT CASE (pt_patch%cell_type)
          CASE (6)
            ! old
            !diff_multfac = hdiff_tv_ratio*k4(k_jg)*4._wp/3._wp
            diff_multfac = hdiff_tv_ratio*k4(k_jg)*4._wp/81._wp/3.0_wp
            ! rejected new
            !diff_multfac = hdiff_tv_ratio*k4(k_jg)*4._wp/81._wp
          CASE (3)
            !old
            diff_multfac = hdiff_tv_ratio*k4(k_jg)*16._wp/27._wp
            !rejected new
            !diff_multfac = hdiff_tv_ratio*k4(k_jg)*16._wp/9._wp
          END SELECT

          IF (.NOT.ltheta_dyn) THEN

!$OMP DO PRIVATE(jb,jk)
            DO jb = i_startblk,i_endblk

              CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                 i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

               DO jk = 1, nlev
                  pt_new%temp(i_startidx:i_endidx,jk,jb) =     &
                    pt_new%temp(i_startidx:i_endidx,jk,jb) -   &
                    z_cell_val(i_startidx:i_endidx,jk,jb) * &
                    pt_patch%cells%area(i_startidx:i_endidx,jb)* &
                    pt_patch%cells%area(i_startidx:i_endidx,jb)*diff_multfac
               ENDDO
            ENDDO
!$OMP END DO

          ELSEIF (ltheta_dyn) THEN

!$OMP DO PRIVATE(jb,jk)
            DO jb = i_startblk,i_endblk

              CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                 i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

              DO jk = 1, nlev
                pt_new%theta(i_startidx:i_endidx,jk,jb) =     &
                  pt_new%theta(i_startidx:i_endidx,jk,jb) -   &
                  z_cell_val(i_startidx:i_endidx,jk,jb) * &
                  pt_patch%cells%area(i_startidx:i_endidx,jb)* &
                  pt_patch%cells%area(i_startidx:i_endidx,jb)*diff_multfac * &
                  pt_diag%delp_c(i_startidx:i_endidx,jk,jb)
              ENDDO
            ENDDO
!$OMP END DO
          ENDIF ! if ltheta_dyn
        ENDIF ! if lhdiff_temp

        IF (grf_intmethod_c == 1 .AND. diffusion_config(k_jg)%lhdiff_temp .AND. &
            .NOT.ltheta_dyn .AND. pt_patch%id > 1) THEN

           ! Lateral boundary diffusion
           i_startblk = pt_patch%cells%start_blk(start_bdydiff_c,1)
           i_endblk   = pt_patch%cells%end_blk(grf_bdywidth_c,1)

           DO jb = i_startblk,i_endblk

             CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, start_bdydiff_c, grf_bdywidth_c)

! OpenMP parallelization is done over jk for boundary diffusion because it affects
! only few blocks
!$OMP DO PRIVATE(jk)
             DO jk = 1, nlev
               pt_new%temp(i_startidx:i_endidx,jk,jb) =   &
                 pt_new%temp(i_startidx:i_endidx,jk,jb) + &
                 z_nabla2_c(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%cells%area(i_startidx:i_endidx,jb)*fac_bdydiff_c
             ENDDO
!$OMP END DO
           ENDDO

        ELSE IF (grf_intmethod_c == 1 .AND. diffusion_config(k_jg)%lhdiff_temp &
          &  .AND. ltheta_dyn .AND. pt_patch%id > 1) THEN

           ! Lateral boundary diffusion
           i_startblk = pt_patch%cells%start_blk(start_bdydiff_c,1)
           i_endblk   = pt_patch%cells%end_blk(grf_bdywidth_c,1)

           DO jb = i_startblk,i_endblk

             CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, start_bdydiff_c, grf_bdywidth_c)

! OpenMP parallelization is done over jk for boundary diffusion because it affects
! only few blocks
!$OMP DO PRIVATE(jk)
             DO jk = 1, nlev
               pt_new%theta(i_startidx:i_endidx,jk,jb) =   &
                 pt_new%theta(i_startidx:i_endidx,jk,jb) + &
                 z_nabla2_c(i_startidx:i_endidx,jk,jb)*   &
                 pt_patch%cells%area(i_startidx:i_endidx,jb)* &
                 pt_diag%delp_c(i_startidx:i_endidx,jk,jb)*fac_bdydiff_c
             ENDDO
!$OMP END DO
           ENDDO

        ENDIF ! temperature boundary diffusion

!$OMP END PARALLEL

     CASE(5) ! 4th order plus 2nd order Smagorinsky diffusion------

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,i_startidx,i_endidx,diff_multfac)

        !---------------------------
        ! velocity
        !---------------------------

        IF (diffusion_config(k_jg)%lhdiff_vn) THEN

           ! Diffusion in the domain interior (boundary diffusion follows below)
           i_startblk = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
           i_endblk   = pt_patch%edges%end_blk(min_rledge,i_nchdom)

           diff_multfac = k4(k_jg)/3._wp

!$OMP DO PRIVATE(jb,jk)
           DO jb = i_startblk,i_endblk

             CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, grf_bdywidth_e+1, min_rledge)

             DO jk = 1, nlev
               pt_new%vn(i_startidx:i_endidx,jk,jb) =       &
                 pt_new%vn(i_startidx:i_endidx,jk,jb)  -    &
                 diff_multfac * z_edge_val(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%edges%area_edge(i_startidx:i_endidx,jb)* &
                 pt_patch%edges%area_edge(i_startidx:i_endidx,jb) + &
                 kh_smag_e(i_startidx:i_endidx,jk,jb) * &
                 z_nabla2_e(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%edges%area_edge(i_startidx:i_endidx,jb)
             ENDDO
           ENDDO
!$OMP END DO

        ENDIF ! if lhdiff_vn

        IF (diffusion_config(k_jg)%lhdiff_vn .AND. pt_patch%id > 1) THEN

           ! Lateral boundary diffusion
           i_startblk = pt_patch%edges%start_blk(start_bdydiff_e,1)
           i_endblk   = pt_patch%edges%end_blk(grf_bdywidth_e,1)

           DO jb = i_startblk,i_endblk

             CALL get_indices_e(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, start_bdydiff_e, grf_bdywidth_e)

! OpenMP parallelization is done over jk for boundary diffusion because it affects
! only few blocks
!$OMP DO PRIVATE(jk)
             DO jk = 1, nlev
               pt_new%vn(i_startidx:i_endidx,jk,jb) =   &
                 pt_new%vn(i_startidx:i_endidx,jk,jb) + &
                 z_nabla2_e(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%edges%area_edge(i_startidx:i_endidx,jb)*fac_bdydiff_e
             ENDDO
!$OMP END DO
           ENDDO

        ENDIF ! vn boundary diffusion

        !---------------------------
        ! temperature
        !---------------------------

        IF (diffusion_config(k_jg)%lhdiff_temp .AND. (.NOT.ltheta_dyn)) THEN

           ! Diffusion in the domain interior (boundary diffusion follows below)
           i_startblk = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
           i_endblk   = pt_patch%cells%end_blk(min_rlcell,i_nchdom)

           SELECT CASE (pt_patch%cell_type)
           CASE (6)
             ! old
             !diff_multfac = hdiff_tv_ratio*k4(k_jg)*4._wp/3._wp
             diff_multfac = hdiff_tv_ratio*k4(k_jg)*4._wp/81._wp/3._wp
             ! rejected new
             !diff_multfac = hdiff_tv_ratio*k4(k_jg)*4._wp/81._wp
           CASE (3)
             ! old
             diff_multfac = hdiff_tv_ratio*k4(k_jg)*16._wp/27._wp
             ! rejected new
             !diff_multfac = hdiff_tv_ratio*k4(k_jg)*16._wp/9._wp
           END SELECT

!$OMP DO PRIVATE(jb,jk)
           DO jb = i_startblk,i_endblk

             CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

              DO jk = 1, nlev
                 pt_new%temp(i_startidx:i_endidx,jk,jb) =     &
                   pt_new%temp(i_startidx:i_endidx,jk,jb) -   &
                   diff_multfac * z_cell_val(i_startidx:i_endidx,jk,jb) * &
                   pt_patch%cells%area(i_startidx:i_endidx,jb)* &
                   pt_patch%cells%area(i_startidx:i_endidx,jb)+ &
                   hdiff_tv_ratio * kh_smag_c(i_startidx:i_endidx,jk,jb) * &
                   z_nabla2_c(i_startidx:i_endidx,jk,jb) * &
                   pt_patch%cells%area(i_startidx:i_endidx,jb)
              ENDDO
           ENDDO
!$OMP END DO

        ELSE IF (diffusion_config(k_jg)%lhdiff_temp .AND. ltheta_dyn) THEN

           ! Diffusion in the domain interior (boundary diffusion follows below)
           i_startblk = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
           i_endblk   = pt_patch%cells%end_blk(min_rlcell,i_nchdom)

           SELECT CASE (pt_patch%cell_type)
           CASE (6)
             ! old
             !diff_multfac = hdiff_tv_ratio*k4(k_jg)*4._wp/3._wp
             diff_multfac = hdiff_tv_ratio*k4(k_jg)*4._wp/81._wp/3._wp 
             ! rejected new
             !diff_multfac = hdiff_tv_ratio*k4(k_jg)*4._wp/81._wp 
           CASE (3)
             ! old
             diff_multfac = hdiff_tv_ratio*k4(k_jg)*16._wp/27._wp
             !diff_multfac = hdiff_tv_ratio*k4(k_jg)*16._wp/9._wp
           END SELECT

!$OMP DO PRIVATE(jb,jk)
           DO jb = i_startblk,i_endblk

             CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell)

             DO jk = 1, nlev
               pt_new%theta(i_startidx:i_endidx,jk,jb) =     &
                 pt_new%theta(i_startidx:i_endidx,jk,jb) +   &
                 pt_patch%cells%area(i_startidx:i_endidx,jb)* &
                 pt_diag%delp_c(i_startidx:i_endidx,jk,jb)* ( &
                 hdiff_tv_ratio * kh_smag_c(i_startidx:i_endidx,jk,jb) * &
                 z_nabla2_c(i_startidx:i_endidx,jk,jb) -   &
                 diff_multfac * z_cell_val(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%cells%area(i_startidx:i_endidx,jb) )
             ENDDO
           ENDDO
!$OMP END DO
        ENDIF ! if lhdiff_temp

        IF (grf_intmethod_c == 1 .AND. diffusion_config(k_jg)%lhdiff_temp &
          &  .AND. .NOT.ltheta_dyn .AND. pt_patch%id > 1) THEN

           ! Lateral boundary diffusion
           i_startblk = pt_patch%cells%start_blk(start_bdydiff_c,1)
           i_endblk   = pt_patch%cells%end_blk(grf_bdywidth_c,1)

           DO jb = i_startblk,i_endblk

             CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, start_bdydiff_c, grf_bdywidth_c)

! OpenMP parallelization is done over jk for boundary diffusion because it affects
! only few blocks
!$OMP DO PRIVATE(jk)
             DO jk = 1, nlev
               pt_new%temp(i_startidx:i_endidx,jk,jb) =   &
                 pt_new%temp(i_startidx:i_endidx,jk,jb) + &
                 z_nabla2_c(i_startidx:i_endidx,jk,jb) * &
                 pt_patch%cells%area(i_startidx:i_endidx,jb)*fac_bdydiff_c
             ENDDO
!$OMP END DO
           ENDDO

        ELSE IF (grf_intmethod_c == 1 .AND. diffusion_config(k_jg)%lhdiff_temp &
          &  .AND. ltheta_dyn .AND. pt_patch%id > 1) THEN

           ! Lateral boundary diffusion
           i_startblk = pt_patch%cells%start_blk(start_bdydiff_c,1)
           i_endblk   = pt_patch%cells%end_blk(grf_bdywidth_c,1)

           DO jb = i_startblk,i_endblk

             CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                                i_startidx, i_endidx, start_bdydiff_c, grf_bdywidth_c)

! OpenMP parallelization is done over jk for boundary diffusion because it affects
! only few blocks
!$OMP DO PRIVATE(jk)
             DO jk = 1, nlev
               pt_new%theta(i_startidx:i_endidx,jk,jb) =   &
                 pt_new%theta(i_startidx:i_endidx,jk,jb) + &
                 z_nabla2_c(i_startidx:i_endidx,jk,jb)*   &
                 pt_patch%cells%area(i_startidx:i_endidx,jb)* &
                 pt_diag%delp_c(i_startidx:i_endidx,jk,jb)*fac_bdydiff_c
             ENDDO
!$OMP END DO
           ENDDO

        ENDIF ! temperature boundary diffusion

!$OMP END PARALLEL

     END SELECT

     CALL sync_patch_array(SYNC_E, pt_patch, pt_new%vn)
     IF(ltheta_dyn) THEN
        CALL sync_patch_array(SYNC_C, pt_patch, pt_new%theta)
     ELSE
        CALL sync_patch_array(SYNC_C, pt_patch, pt_new%temp)
     ENDIF

  END SUBROUTINE hdiff_expl

END MODULE mo_hdiff
