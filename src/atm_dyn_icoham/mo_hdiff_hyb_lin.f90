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
MODULE mo_hdiff_hyb_lin

  USE mo_kind,              ONLY: wp
  USE mo_diffusion_nml,     ONLY: lhdiff_vn, lhdiff_temp,    &
                                & hdiff_tv_ratio, k2, k4,    &
                                & k2s, k2e, k4s, k4e
  USE mo_math_operators,    ONLY: nabla2_vec, nabla2_scalar, &
                                & nabla4_vec, nabla4_scalar   
  USE mo_model_domain,      ONLY: t_patch
  USE mo_interpolation,     ONLY: t_int_state
  USE mo_icoham_dyn_types,  ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_parallel_configuration,  ONLY: nproma
  USE mo_run_nml,           ONLY: nlev, i_cell_type, ltheta_dyn
  USE mo_impl_constants,    ONLY: min_rlcell, min_rledge
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_grf_interpolation, ONLY: denom_diffu_v, denom_diffu_t, grf_intmethod_c
  USE mo_sync,              ONLY: SYNC_C, SYNC_E, sync_patch_array

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: hdiff_hyb_lin 
  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
  SUBROUTINE hdiff_hyb_lin( jg, patch, pint, diag, prog )

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: jg      !< grid level

    TYPE(t_patch),    TARGET,INTENT(IN)    :: patch  !< grid info
    TYPE(t_int_state),TARGET,INTENT(IN)    :: pint   !< interpolation coeff.
    TYPE(t_hydro_atm_diag),  INTENT(IN)    :: diag   !< diagnostic variables
    TYPE(t_hydro_atm_prog),  INTENT(INOUT) :: prog   !< prognostic variables

    ! Temporary arrays
    REAL(wp) :: znabla_c (nproma,nlev,patch%nblks_c)
    REAL(wp) :: znabla_e (nproma,nlev,patch%nblks_e)
    REAL(wp) :: ztmp_c   (nproma,nlev,patch%nblks_c)

    ! Temporary arrays for refinement
    REAL(wp) :: znabla2_c(nproma,nlev,patch%nblks_c)
    REAL(wp) :: znabla2_e(nproma,nlev,patch%nblks_e)

    ! Dimension sizes, indices
    INTEGER :: nchilddom
    INTEGER :: jb, jbs, jbe, is, ie, jk
    INTEGER :: start_bdydiff_c, start_bdydiff_e  !< for boundary diffusion

    ! Diffusion coefficients
    REAL(wp) :: zc2temp, zc4temp              !< for domain interia, temp 
    REAL(wp) :: zc2vn,   zc4vn                !< for domain interia, vn
    REAL(wp) :: fac_bdydiff_c, fac_bdydiff_e  !< for boundary of refined region

    !===================================================================
    ! 0. Some constants
    !===================================================================
    ! Number of child domains
    nchilddom = MAX(1,patch%n_childdom)

    ! Diffusion coefficients for velocity and (potential) temperature

    zc2vn = k2(jg)*1._wp/SQRT(3._wp)
    zc4vn = k4(jg)/3._wp

    SELECT CASE (i_cell_type)
    CASE (6)
      zc2temp = hdiff_tv_ratio*k2(jg)*2._wp/9._wp/SQRT(3._wp)
      zc4temp = hdiff_tv_ratio*k4(jg)*4._wp/81._wp/3.0_wp
    CASE (3)
      zc2temp = hdiff_tv_ratio*k2(jg)*4._wp/3._wp/SQRT(3._wp)
      zc4temp = hdiff_tv_ratio*k4(jg)*16._wp/27._wp
    END SELECT

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
    IF (lhdiff_vn) THEN
      CALL nabla2_vec( prog%vn, patch, pint, znabla_e,     &
                       opt_slev=k2s,  opt_elev=k2e,        &
                       opt_rlstart=5, opt_rlend=min_rledge )

      CALL nabla4_vec( prog%vn, patch, pint, znabla_e,     &
                       opt_nabla2=znabla2_e,               &
                       opt_slev=k4s,  opt_elev=k4e,        &
                       opt_rlstart=7, opt_rlend=min_rledge )
    END IF
    !-------------
    ! Temperature 
    !-------------
    IF (lhdiff_temp.AND.(.NOT.ltheta_dyn)) THEN

      CALL nabla2_scalar( prog%temp, patch, pint, znabla_c,   &
                          opt_slev=k2s,  opt_elev=k2e,        &
                          opt_rlstart=3, opt_rlend=min_rlcell )

      CALL nabla4_scalar( prog%temp, patch, pint, znabla_c,   &
                          opt_nabla2=znabla2_c,               &
                          opt_slev=k4s,  opt_elev=k4e,        &
                          opt_rlstart=4, opt_rlend=min_rlcell )

    !---------------------------------------------------------------
    ! Potential temperature, if not using temperature as the 
    ! thermodynamic variable. Note that in this case the component
    ! %theta of the prognostic state is in fact theta*delta_p,
    ! defined at cell centers on full levels.
    !---------------------------------------------------------------
    ELSE IF (lhdiff_temp .AND. ltheta_dyn) THEN

      ! Divide by delta p so as to diffuse the uncoupled theta
      jbs = patch%cells%start_blk(1,1)
      jbe = patch%cells%end_blk(min_rlcell,nchilddom)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
      DO jb = jbs,jbe 
        CALL get_indices_c(patch,jb,jbs,jbe,is,ie,1,min_rlcell)
        ! Note: we use delp_c here because convert_t2theta/convert_theta2t
        ! update delp_c but not rdelp_c
        ztmp_c(is:ie,:,jb) =  prog%theta (is:ie,:,jb) &
                           & /diag%delp_c(is:ie,:,jb)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      CALL nabla2_scalar( ztmp_c, patch, pint, znabla_c,      &
                          opt_slev=k2s,  opt_elev=k2e,        &
                          opt_rlstart=3, opt_rlend=min_rlcell )

      CALL nabla4_scalar( ztmp_c, patch, pint, znabla_c,      &
                          opt_nabla2=znabla2_c,               &
                          opt_slev=k4s,  opt_elev=k4e,        &
                          opt_rlstart=4, opt_rlend=min_rlcell )
    ENDIF

!$OMP PARALLEL PRIVATE(jbs,jbe,is,ie)
    !===================================================================
    ! 2. Apply diffusion to the prognostic variables in the domain 
    !    interior. Boundary diffusion follows in step 3.
    !===================================================================
    ! Velocity
    !-----------
    IF (lhdiff_vn) THEN

      jbs = patch%edges%start_blk(grf_bdywidth_e+1,1)
      jbe = patch%edges%end_blk(min_rledge,nchilddom)

!$OMP DO PRIVATE(jb,jk)
      DO jb = jbs,jbe
        CALL get_indices_e(patch,jb,jbs,jbe,is,ie,grf_bdywidth_e+1,min_rledge)

        DO jk = k2s,k2e
          prog%vn(is:ie,jk,jb) =   prog%vn(is:ie,jk,jb)           &
                               & +znabla_e(is:ie,jk,jb)*zc2vn     &
                               & *patch%edges%area_edge(is:ie,jb)
        ENDDO !jk

        DO jk = k4s,k4e
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
    IF (lhdiff_temp) THEN

      jbs = patch%cells%start_blk(grf_bdywidth_c+1,1)
      jbe = patch%cells%end_blk(min_rlcell,nchilddom)

      IF (.NOT.ltheta_dyn) THEN !---- temperature -----

!$OMP DO PRIVATE(jb,jk)
        DO jb = jbs,jbe
          CALL get_indices_c(patch,jb,jbs,jbe,is,ie,grf_bdywidth_c+1, min_rlcell)

          DO jk = k2s,k2e
            prog%temp(is:ie,jk,jb) =  prog%temp(is:ie,jk,jb)        &
                                   & +znabla_c(is:ie,jk,jb)*zc2temp &
                                   & *patch%cells%area(is:ie,jb)
          ENDDO !jk

          DO jk = k4s,k4e
            prog%temp(is:ie,jk,jb) =  prog%temp(is:ie,jk,jb)        &
                                   & -znabla_c(is:ie,jk,jb)*zc4temp & 
                                   & *patch%cells%area(is:ie,jb)**2
          ENDDO !jk
        ENDDO   !jb
!$OMP END DO

      ELSE !---- potential temperature -----

!$OMP DO PRIVATE(jb,jk)
        DO jb = jbs,jbe
          CALL get_indices_c(patch,jb,jbs,jbe,is,ie,grf_bdywidth_c+1, min_rlcell)

          DO jk = k2s,k2e
            prog%theta(is:ie,jk,jb) =  prog%theta(is:ie,jk,jb)       &
                                    & +znabla_c(is:ie,jk,jb)*zc2temp &
                                    & *patch%cells%area(is:ie,jb)    &
                                    & *diag%delp_c(is:ie,jk,jb)
          ENDDO !jk

          DO jk = k4s,k4e
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
    ! - stored in variable znabla_e/c  for vertical levels [k2s,k2e], and  
    ! - stored in variable znabla2_e/c for vertical levels [k4s,k4e].
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

!$OMP DO PRIVATE(jk)
        DO jk = k2s,k2e
          prog%vn(is:ie,jk,jb) =  prog%vn(is:ie,jk,jb)                &
                               & +znabla_e(is:ie,jk,jb)*fac_bdydiff_e &
                               & *patch%edges%area_edge(is:ie,jb)
        ENDDO !jk
!$OMP END DO
!$OMP DO PRIVATE(jk)
        DO jk = k4s,k4e
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

    IF ((patch%id>1).AND.lhdiff_temp.AND.(grf_intmethod_c==1)) THEN

      jbs = patch%cells%start_blk(start_bdydiff_c,1)
      jbe = patch%cells%end_blk(grf_bdywidth_c,1)

      IF (.NOT.ltheta_dyn) THEN !--- temperature ---

        DO jb = jbs,jbe
          CALL get_indices_c(patch,jb,jbs,jbe,is,ie,start_bdydiff_c,grf_bdywidth_c)

!$OMP DO PRIVATE(jk)
          DO jk = k2s,k2e
            prog%temp(is:ie,jk,jb) =  prog%temp(is:ie,jk,jb)              &
                                   & +znabla_c(is:ie,jk,jb)*fac_bdydiff_c &
                                   & *patch%cells%area(is:ie,jb)
          ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jk)
          DO jk = k4s,k4e
            prog%temp(is:ie,jk,jb) =  prog%temp(is:ie,jk,jb)               &
                                   & +znabla2_c(is:ie,jk,jb)*fac_bdydiff_c &
                                   & *patch%cells%area(is:ie,jb)
          ENDDO
!$OMP END DO

        ENDDO !jb

      ELSE !--- potential temperature ---

        DO jb = jbs,jbe
          CALL get_indices_c(patch,jb,jbs,jbe,is,ie,start_bdydiff_c,grf_bdywidth_c)

!$OMP DO PRIVATE(jk)
          DO jk = k2s,k2e
            prog%theta(is:ie,jk,jb) =  prog%theta(is:ie,jk,jb)             &
                                    & +znabla_c(is:ie,jk,jb)*fac_bdydiff_c &
                                    & *patch%cells%area(is:ie,jb)          &
                                    & *diag%delp_c(is:ie,jk,jb)
          ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jk)
          DO jk = k4s,k4e
            prog%theta(is:ie,jk,jb) =  prog%theta(is:ie,jk,jb)              &
                                    & +znabla2_c(is:ie,jk,jb)*fac_bdydiff_c &
                                    & *patch%cells%area(is:ie,jb)           &
                                    & *diag%delp_c(is:ie,jk,jb)
          ENDDO
!$OMP END DO

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

