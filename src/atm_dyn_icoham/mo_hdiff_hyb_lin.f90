!>
!! Linear horizontal diffusion 
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Hui Wan, MPI 
!!
!! @par Revision History
!! First version by Hui Wan, MPI (2011-02-28) 
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
MODULE mo_hdiff_hyb_lin

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
  USE mo_diffusion_config,  ONLY: diffusion_config
  USE mo_ha_dyn_config,     ONLY: ha_dyn_config
  USE mo_math_laplace,      ONLY: nabla2_vec, nabla2_scalar, &
                                & nabla4_vec, nabla4_scalar   
  USE mo_model_domain,      ONLY: t_patch
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_icoham_dyn_types,  ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_parallel_config,   ONLY: nproma, p_test_run
  USE mo_impl_constants,    ONLY: min_rlcell, min_rledge
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_gridref_config,    ONLY: denom_diffu_v, denom_diffu_t, grf_intmethod_c
  USE mo_sync,              ONLY: SYNC_C, SYNC_E, sync_patch_array

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: hdiff_hyb_lin 

CONTAINS
  !>
  !!
  SUBROUTINE hdiff_hyb_lin( jg, patch, pint, diag, prog )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: jg      !< grid level

    TYPE(t_patch),    TARGET,INTENT(INOUT) :: patch  !< grid info
    TYPE(t_int_state),TARGET,INTENT(IN)    :: pint   !< interpolation coeff.
    TYPE(t_hydro_atm_diag),  INTENT(IN)    :: diag   !< diagnostic variables
    TYPE(t_hydro_atm_prog),  INTENT(INOUT) :: prog   !< prognostic variables

    ! Temporary arrays
    REAL(wp) :: znabla_c (nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: znabla_e (nproma,patch%nlev,patch%nblks_e)
    REAL(wp) :: ztmp_c   (nproma,patch%nlev,patch%nblks_c)

    ! Temporary arrays for refinement
    REAL(wp) :: znabla2_c(nproma,patch%nlev,patch%nblks_c)
    REAL(wp) :: znabla2_e(nproma,patch%nlev,patch%nblks_e)

    ! Dimension sizes, indices
    INTEGER :: nchilddom
    INTEGER :: jb, jbs, jbe, is, ie, jk
    INTEGER :: start_bdydiff_c, start_bdydiff_e  !< for boundary diffusion

    ! Diffusion coefficients
    REAL(wp) :: zc2temp, zc4temp              !< for domain interia, temp 
    REAL(wp) :: zc2vn,   zc4vn                !< for domain interia, vn
    REAL(wp) :: fac_bdydiff_c, fac_bdydiff_e  !< for boundary of refined region

    LOGICAL :: ltheta_dyn, lhdiff_vn, lhdiff_temp
    INTEGER :: nlev, ik2s, ik2e, ik4s, ik4e
    REAL(wp):: k2, k4, hdiff_tv_ratio

    !===================================================================
    ! 0. Some constants
    !===================================================================
    ! CALL sync_patch_array(SYNC_C, patch, prog%temp)
    ! CALL sync_patch_array(SYNC_E, patch, prog%vn)

    ltheta_dyn = ha_dyn_config%ltheta_dyn
    nlev       = patch%nlev

    hdiff_tv_ratio = diffusion_config(jg)% hdiff_tv_ratio 
    lhdiff_temp    = diffusion_config(jg)% lhdiff_temp
    lhdiff_vn      = diffusion_config(jg)% lhdiff_vn
    ik2s           = diffusion_config(jg)% ik2s
    ik2e           = diffusion_config(jg)% ik2e
    ik4s           = diffusion_config(jg)% ik4s
    ik4e           = diffusion_config(jg)% ik4e
    k2             = diffusion_config(jg)% k2
    k4             = diffusion_config(jg)% k4

    ! Number of child domains
    nchilddom = MAX(1,patch%n_childdom)

    ! Diffusion coefficients for velocity and (potential) temperature

    zc2vn = k2*1._wp/SQRT(3._wp)
    zc4vn = k4/3._wp

    zc2temp = hdiff_tv_ratio*k2*4._wp/3._wp/SQRT(3._wp)
    zc4temp = hdiff_tv_ratio*k4*16._wp/27._wp

    ! For nested domain

    start_bdydiff_c = 3
    start_bdydiff_e = 5

    fac_bdydiff_c = 1._wp/denom_diffu_t  ! normalized coefficient
    fac_bdydiff_e = 1._wp/denom_diffu_v  ! normalized coefficient 

    !===================================================================
    ! 1. Compute Laplacian (nabla2 and/or nabla4)
    !===================================================================
    ! Velocity
    !-----------
    IF (diffusion_config(jg)%lhdiff_vn) THEN
  !    CALL sync_patch_array(SYNC_e, patch, prog%vn)
  !    write(0,*) "start nabla vn"
      IF (p_test_run) znabla_e(:,:,:) = 0.0_wp
      CALL nabla2_vec( prog%vn, patch, pint, znabla_e,     &
                       opt_slev=ik2s,  opt_elev=ik2e,      &
                       opt_rlstart=5, opt_rlend=min_rledge )

      CALL nabla4_vec( prog%vn, patch, pint, znabla_e,     &
                       opt_nabla2=znabla2_e,               &
                       opt_slev=ik4s,  opt_elev=ik4e,      &
                       opt_rlstart=7, opt_rlend=min_rledge )
  !     write(0,*) "start sync znabla_e"
  !    CALL sync_patch_array(SYNC_e, patch, znabla_e)
    END IF
    !-------------
    ! Temperature 
    !-------------
    IF (diffusion_config(jg)%lhdiff_temp.AND.(.NOT.ltheta_dyn)) THEN

      ! CALL sync_patch_array(SYNC_C, patch, prog%temp)
      IF (p_test_run) znabla_c(:,:,:) = 0.0_wp
      CALL nabla2_scalar( prog%temp, patch, pint, znabla_c,   &
                          opt_slev=ik2s,  opt_elev=ik2e,      &
                          opt_rlstart=3, opt_rlend=min_rlcell )

      CALL nabla4_scalar( prog%temp, patch, pint, znabla_c,   &
                          opt_nabla2=znabla2_c,               &
                          opt_slev=ik4s,  opt_elev=ik4e,      &
                          opt_rlstart=4, opt_rlend=min_rlcell )

!       write(0,*) "start sync znabla_c"
!       CALL sync_patch_array(SYNC_C, patch, znabla_c)

   !---------------------------------------------------------------
    ! Potential temperature, if not using temperature as the 
    ! thermodynamic variable. Note that in this case the component
    ! %theta of the prognostic state is in fact theta*delta_p,
    ! defined at cell centers on full levels.
    !---------------------------------------------------------------
    ELSE IF (diffusion_config(jg)%lhdiff_temp .AND. ltheta_dyn) THEN

      ! Divide by delta p so as to diffuse the uncoupled theta
      jbs = patch%cells%start_blk(1,1)
      jbe = patch%cells%end_blk(min_rlcell,nchilddom)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,jbe 
        CALL get_indices_c(patch,jb,jbs,jbe,is,ie,1,min_rlcell)
        ! Note: we use delp_c here because convert_t2theta/convert_theta2t
        ! update delp_c but not rdelp_c
        ztmp_c(is:ie,:,jb) =  prog%theta (is:ie,:,jb) &
                           & /diag%delp_c(is:ie,:,jb)
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      CALL nabla2_scalar( ztmp_c, patch, pint, znabla_c,      &
                          opt_slev=ik2s,  opt_elev=ik2e,      &
                          opt_rlstart=3, opt_rlend=min_rlcell )

      CALL nabla4_scalar( ztmp_c, patch, pint, znabla_c,      &
                          opt_nabla2=znabla2_c,               &
                          opt_slev=ik4s,  opt_elev=ik4e,      &
                          opt_rlstart=4, opt_rlend=min_rlcell )
    ENDIF

!$OMP PARALLEL PRIVATE(jbs,jbe,is,ie)
    !===================================================================
    ! 2. Apply diffusion to the prognostic variables in the domain 
    !    interior. Boundary diffusion follows in step 3.
    !===================================================================
    ! Velocity
    !-----------
    IF (diffusion_config(jg)%lhdiff_vn) THEN

      jbs = patch%edges%start_blk(grf_bdywidth_e+1,1)
      jbe = patch%edges%end_blk(min_rledge,nchilddom)

!$OMP DO PRIVATE(jb,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,jbe
        CALL get_indices_e(patch,jb,jbs,jbe,is,ie,grf_bdywidth_e+1,min_rledge)

        DO jk = ik2s,ik2e
          prog%vn(is:ie,jk,jb) =   prog%vn(is:ie,jk,jb)           &
                               & +znabla_e(is:ie,jk,jb)*zc2vn     &
                               & *patch%edges%area_edge(is:ie,jb)
        ENDDO !jk

        DO jk = ik4s,ik4e
          prog%vn(is:ie,jk,jb) =   prog%vn(is:ie,jk,jb)             &
                               & -znabla_e(is:ie,jk,jb)*zc4vn       &
                               & *patch%edges%area_edge(is:ie,jb)**2
        ENDDO !jk
      ENDDO !jb
!$OMP END DO
    ENDIF ! lhdiff_vn

    !-------------------------
    ! Thermodynamic variable 
    !-------------------------
    IF (diffusion_config(jg)%lhdiff_temp) THEN

      jbs = patch%cells%start_blk(grf_bdywidth_c+1,1)
      jbe = patch%cells%end_blk(min_rlcell,nchilddom)

      IF (.NOT.ltheta_dyn) THEN !---- temperature -----

!$OMP DO PRIVATE(jb,jk) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          CALL get_indices_c(patch,jb,jbs,jbe,is,ie,grf_bdywidth_c+1, min_rlcell)

          DO jk = ik2s,ik2e
            prog%temp(is:ie,jk,jb) =  prog%temp(is:ie,jk,jb)        &
                                   & +znabla_c(is:ie,jk,jb)*zc2temp &
                                   & *patch%cells%area(is:ie,jb)
          ENDDO !jk

          DO jk = ik4s,ik4e
            prog%temp(is:ie,jk,jb) =  prog%temp(is:ie,jk,jb)        &
                                   & -znabla_c(is:ie,jk,jb)*zc4temp & 
                                   & *patch%cells%area(is:ie,jb)**2
          ENDDO !jk
        ENDDO   !jb
!$OMP END DO

      ELSE !---- potential temperature -----

!$OMP DO PRIVATE(jb,jk) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = jbs,jbe
          CALL get_indices_c(patch,jb,jbs,jbe,is,ie,grf_bdywidth_c+1, min_rlcell)

          DO jk = ik2s,ik2e
            prog%theta(is:ie,jk,jb) =  prog%theta(is:ie,jk,jb)       &
                                    & +znabla_c(is:ie,jk,jb)*zc2temp &
                                    & *patch%cells%area(is:ie,jb)    &
                                    & *diag%delp_c(is:ie,jk,jb)
          ENDDO !jk

          DO jk = ik4s,ik4e
            prog%theta(is:ie,jk,jb) =  prog%theta(is:ie,jk,jb)       &
                                    & -znabla_c(is:ie,jk,jb)*zc4temp &
                                    & *patch%cells%area(is:ie,jb)**2 &
                                    & *diag%delp_c(is:ie,jk,jb)
          ENDDO !jk
        ENDDO  !jb
!$OMP END DO

      ENDIF ! if ltheta_dyn
    ENDIF   ! if lhdiff_temp

 
    !===================================================================
    ! 3. Apply 2nd order diffusion for the boundary of refined region
    !===================================================================
    ! NOTE that after the computations in step 2, we have the Laplacian 
    ! (nabla2) of the prognostic variables 
    ! - stored in variable znabla_e/c  for vertical levels [ik2s,ik2e], and  
    ! - stored in variable znabla2_e/c for vertical levels [ik4s,ik4e].
    ! Here for the boundary of refined region, we apply 2nd order 
    ! diffusion for all vertical levels.
    ! TECHNICAL NOTE: OpenMP parallelization is done over jk for 
    ! boundary diffusion because it affects only few blocks.
    !-----------
    ! Velocity
    !-----------
    IF ((patch%id > 1).AND.lhdiff_vn) THEN

      jbs = patch%edges%start_blk(start_bdydiff_e,1)
      jbe = patch%edges%end_blk(grf_bdywidth_e,1)

      DO jb = jbs,jbe
        CALL get_indices_e(patch,jb,jbs,jbe,is,ie,start_bdydiff_e, grf_bdywidth_e)

!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE
        DO jk = ik2s,ik2e
          prog%vn(is:ie,jk,jb) =  prog%vn(is:ie,jk,jb)                &
                               & +znabla_e(is:ie,jk,jb)*fac_bdydiff_e &
                               & *patch%edges%area_edge(is:ie,jb)
        ENDDO !jk
!$OMP END DO
!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE
        DO jk = ik4s,ik4e
          prog%vn(is:ie,jk,jb) =  prog%vn(is:ie,jk,jb)                 &
                               & +znabla2_e(is:ie,jk,jb)*fac_bdydiff_e &
                               & *patch%edges%area_edge(is:ie,jb)
        ENDDO !jk
!$OMP END DO

      ENDDO ! jb
    ENDIF ! patch%id >1 and lhdiff_vn

    !------------------------
    ! Thermodynamic variable
    !------------------------
    ! Necessary only when grf_intmethod_c == 1 (i.e., when the interpolation
    ! method used for the the boundary cells is simple copying).

    IF ((patch%id>1).AND.lhdiff_temp  &
      &  .AND.(grf_intmethod_c==1)) THEN

      jbs = patch%cells%start_blk(start_bdydiff_c,1)
      jbe = patch%cells%end_blk(grf_bdywidth_c,1)

      IF (.NOT.ltheta_dyn) THEN !--- temperature ---

        DO jb = jbs,jbe
          CALL get_indices_c(patch,jb,jbs,jbe,is,ie,start_bdydiff_c,grf_bdywidth_c)

!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE
          DO jk = ik2s,ik2e
            prog%temp(is:ie,jk,jb) =  prog%temp(is:ie,jk,jb)              &
                                   & +znabla_c(is:ie,jk,jb)*fac_bdydiff_c &
                                   & *patch%cells%area(is:ie,jb)
          ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE
          DO jk = ik4s,ik4e
            prog%temp(is:ie,jk,jb) =  prog%temp(is:ie,jk,jb)               &
                                   & +znabla2_c(is:ie,jk,jb)*fac_bdydiff_c &
                                   & *patch%cells%area(is:ie,jb)
          ENDDO
!$OMP END DO NOWAIT

        ENDDO !jb

      ELSE !--- potential temperature ---

        DO jb = jbs,jbe
          CALL get_indices_c(patch,jb,jbs,jbe,is,ie,start_bdydiff_c,grf_bdywidth_c)

!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE
          DO jk = ik2s,ik2e
            prog%theta(is:ie,jk,jb) =  prog%theta(is:ie,jk,jb)             &
                                    & +znabla_c(is:ie,jk,jb)*fac_bdydiff_c &
                                    & *patch%cells%area(is:ie,jb)          &
                                    & *diag%delp_c(is:ie,jk,jb)
          ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE
          DO jk = ik4s,ik4e
            prog%theta(is:ie,jk,jb) =  prog%theta(is:ie,jk,jb)              &
                                    & +znabla2_c(is:ie,jk,jb)*fac_bdydiff_c &
                                    & *patch%cells%area(is:ie,jb)           &
                                    & *diag%delp_c(is:ie,jk,jb)
          ENDDO
!$OMP END DO NOWAIT

        ENDDO
      ENDIF ! temperature or potential temperature
    ENDIF ! lhdiff_temp.AND.(grf_intmethod_c==1)
!$OMP END PARALLEL

    !===================================================================
    ! 4. Computations done. Now synchronize.
    !===================================================================
    CALL sync_patch_array(SYNC_E, patch, prog%vn)
    IF (ltheta_dyn) THEN
      CALL sync_patch_array(SYNC_C, patch, prog%theta)
    ELSE
      CALL sync_patch_array(SYNC_C, patch, prog%temp)
    ENDIF

  END SUBROUTINE hdiff_hyb_lin
  !-------------
END MODULE mo_hdiff_hyb_lin

