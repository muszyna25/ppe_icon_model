!>
!! Contains the implementation of the horizontal tracer transport routines for the ICON ocean model.
!! This comprises vertical advection and diffusion
!! 
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01) 
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
MODULE mo_oce_tracer_transport_vert
!-------------------------------------------------------------------------  
!
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2006
!  
!-------------------------------------------------------------------------  
!
!
!
USE mo_kind,                      ONLY: wp
!USE mo_math_utilities,            ONLY: t_cartesian_coordinates
USE mo_impl_constants,            ONLY: sea_boundary, MIN_DOLIC, min_rlcell !, min_rledge
USE mo_ocean_nml,                 ONLY: n_zlev, expl_vertical_tracer_diff, ab_const, irelax_2d_S,&
  &                                     temperature_relaxation!, ab_gam
USE mo_parallel_config,           ONLY: nproma
USE mo_dynamics_config,           ONLY: nold, nnew
USE mo_run_config,                ONLY: dtime, ltimer
USE mo_timer,                     ONLY: timer_start, timer_stop, timer_adv_vert, timer_ppm_slim, &
  &                                     timer_dif_vert
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, v_base, is_initial_timestep
USE mo_model_domain,              ONLY: t_patch
USE mo_exception,                 ONLY: finish !, message_text, message
USE mo_oce_index,                 ONLY: print_mxmn, ipl_src! , jkc, jkdim
USE mo_loopindices,               ONLY: get_indices_c !, get_indices_e, get_indices_v
USE mo_oce_physics
!USE mo_scalar_product,            ONLY:  map_cell2edges,map_edges2cell,map_edges2cell
!USE mo_oce_math_operators,        ONLY: div_oce, grad_fd_norm_oce, grad_fd_norm_oce_2d
USE mo_advection_utils,           ONLY: laxfr_upflux_v
USE mo_oce_diffusion,             ONLY: tracer_diffusion_vert_expl,&
                                      & tracer_diffusion_vert_impl_hom
USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
USE mo_sync,                      ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

!
! PUBLIC INTERFACE
!
PUBLIC :: advect_diffuse_vertical
! Private implemenation
!
PRIVATE :: upwind_vflux_oce
PRIVATE :: central_vflux_oce
PRIVATE :: mimetic_vflux_oce
PRIVATE :: upwind_vflux_ppm
PRIVATE :: v_ppm_slimiter_mo
PRIVATE :: apply_tracer_flux_top_layer_oce


INTEGER, PARAMETER :: UPWIND = 1
INTEGER, PARAMETER :: CENTRAL= 2
INTEGER, PARAMETER :: MIMETIC= 3
INTEGER, PARAMETER :: MIMETIC_MIURA= 4

CONTAINS
  !-------------------------------------------------------------------------
  !! SUBROUTINE advects vertically the tracers present in the ocean model.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2010).
  !!
  !! mpi parallelized, sync required: trac_out
  SUBROUTINE advect_diffuse_vertical(p_patch, trac_in,       &
                           & p_os,                           &
                           & bc_top_tracer, bc_bot_tracer,   &
                           & A_v,                            &
                           & trac_out, timestep, delta_t,    &
                           & cell_thick_intermed_c,          &
                           & FLUX_CALCULATION_VERT, tracer_id)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    REAL(wp), INTENT(IN)              :: trac_in(nproma,n_zlev, p_patch%nblks_c)
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    REAL(wp)                          :: bc_top_tracer(nproma, p_patch%nblks_c)
    REAL(wp)                          :: bc_bot_tracer(nproma, p_patch%nblks_c)
    REAL(wp)                          :: A_v(:,:,:)                               !vertical mixing coeff
    REAL(wp), INTENT(OUT)             :: trac_out(nproma,n_zlev, p_patch%nblks_c) !new tracer
    INTEGER                           :: timestep                                 ! Actual timestep (to distinghuish initial step from others)
    REAL(wp)                          :: delta_t
    REAL(wp)                          :: cell_thick_intermed_c(nproma,n_zlev, p_patch%nblks_c)
    INTEGER                           :: FLUX_CALCULATION_VERT
    INTEGER, INTENT(IN)               :: tracer_id

    !Local variables
    REAL(wp) :: delta_z, delta_z2
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: jc, jk, jb!, je!jkp1        !< index of edge, vert level, block
    INTEGER  :: z_dolic
    REAL(wp) :: z_adv_flux_v (nproma, n_zlev+1, p_patch%nblks_c)  ! vertical advective tracer flux
    REAL(wp) :: z_div_adv_v (nproma, n_zlev,p_patch%nblks_c)        ! vertical tracer divergence
    REAL(wp) :: z_div_diff_v (nproma, n_zlev,p_patch%nblks_c)        ! vertical tracer divergence
    REAL(wp) :: z_h(nproma,n_zlev, p_patch%nblks_c)
    REAL(wp) :: z_temp(nproma,n_zlev, p_patch%nblks_c)
    LOGICAL  :: ldbg = .FALSE.

    REAL(wp) :: z_diff_flux_v(nproma, n_zlev+1,p_patch%nblks_c)   ! vertical diffusive tracer flux
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_tracer_advection:advect_individual_tracer')
    !-------------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------------
    cells_in_domain => p_patch%cells%in_domain

    z_adv_flux_v   = 0.0_wp
    z_div_adv_v    = 0.0_wp
    z_div_diff_v   = 0.0_wp
    z_diff_flux_v  = 0.0_wp



    ! Initialize timer for horizontal advection
    IF (ltimer) CALL timer_start(timer_adv_vert)

    SELECT CASE(FLUX_CALCULATION_VERT)

    CASE(UPWIND)

      CALL upwind_vflux_oce( p_patch,                    &
                           & trac_in,                    &
                           & p_os%p_diag%w_time_weighted,& 
                           & bc_top_tracer,              &
                           & z_adv_flux_v,tracer_id )
    CASE(CENTRAL)
      CALL central_vflux_oce( p_patch,                   &
                           & trac_in,                    &
                           & p_os%p_diag%w_time_weighted,&
                           & z_adv_flux_v, tracer_id)
    CASE(MIMETIC,MIMETIC_MIURA)
      CALL upwind_vflux_ppm( p_patch, trac_in,           &
                           & p_os%p_diag%w_time_weighted,&
                           & dtime, 1 ,                  & !p_itype_vlimit,  &
                           & cell_thick_intermed_c,      &!p_cellhgt_mc_now, &
                           & z_adv_flux_v, tracer_id)
    END SELECT

    IF (ltimer) CALL timer_stop(timer_adv_vert)

    !divergence is calculated for advective fluxes
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        !interior: from one below surface to the ground
        z_dolic = v_base%dolic_c(jc,jb)
        IF(z_dolic>=MIN_DOLIC)THEN
          DO jk = 1, z_dolic
            ! positive vertical divergence in direction of w (upward positive)
            z_div_adv_v(jc,jk,jb) = z_adv_flux_v(jc,jk,jb) &
                                  &- z_adv_flux_v(jc,jk+1,jb)

          END DO
        ENDIF
      END DO
    END DO

    ipl_src = 5  ! output print level (1-5, fix)
    DO jk = 1, n_zlev
      CALL print_mxmn('adv flux_v',jk,z_adv_flux_v(:,:,:),n_zlev+1, &
        &              p_patch%nblks_c,'trc',ipl_src)
      IF (ldbg) THEN
        write(*,*)'vertical adv:',jk,minval(z_adv_flux_v(:,jk,:)),&
        &maxval(z_adv_flux_v(:,jk,:))
      ENDIF
    END DO
    DO jk = 1, n_zlev
      CALL print_mxmn('div adv-flux_v',jk,z_div_adv_v(:,:,:),n_zlev, &
        &              p_patch%nblks_c,'trc',ipl_src)
     IF (ldbg) THEN
      write(*,*)'vertical div:',jk,minval(z_div_adv_v(:,jk,:)),&
      &maxval(z_div_adv_v(:,jk,:))
     ENDIF
    END DO

    IF (ltimer) CALL timer_start(timer_dif_vert)

    !Case: Implicit Vertical diffusion
    IF(expl_vertical_tracer_diff==1)THEN

      !Add advective part to old tracer
      !surface forcing applied as volume forcing at rhs, i.e.part of explicit term in momentum and tracer eqs.
      !in this case, top boundary ondition of vertical Laplacians are homogeneous
      jk=1
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c
            z_dolic = v_base%dolic_c(jc,jb)
            IF(z_dolic>=MIN_DOLIC)THEN
                delta_z = v_base%del_zlev_m(jk)+p_os%p_prog(nnew(1))%h(jc,jb)

                 z_temp(jc,jk,jb)= (trac_in(jc,jk,jb)*cell_thick_intermed_c(jc,jk,jb)     &
                 & -delta_t*(z_div_adv_v(jc,jk,jb)-bc_top_tracer(jc,jb)))&
                 &/delta_z
            ENDIF
        END DO
      END DO

      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
          DO jc = i_startidx_c, i_endidx_c
            z_dolic = v_base%dolic_c(jc,jb)
            IF(z_dolic>=MIN_DOLIC)THEN
              DO jk = 2, z_dolic
                  delta_z = v_base%del_zlev_m(jk)

                  delta_z = v_base%del_zlev_m(jk)
                  z_temp(jc,jk,jb)= (trac_in(jc,jk,jb)*cell_thick_intermed_c(jc,jk,jb)&
                  & -delta_t*z_div_adv_v(jc,jk,jb))/delta_z
 
               END DO
            ENDIF
        END DO
      END DO

      IF (ldbg) THEN
        DO jk = 1, n_zlev
          write(*,*)'before impl-diff: max/min:',jk,&
          &maxval(z_temp(:,jk,:)), minval(z_temp(:,jk,:))
        END DO
      ENDIF
      ipl_src = 5  ! output print level (1-5, fix)
      DO jk = 1, n_zlev
        CALL print_mxmn('bef impl v-trc:',jk,z_temp(:,:,:),n_zlev, &
        &              p_patch%nblks_c,'trc',ipl_src)
      END DO
      DO jk = 1, n_zlev
        CALL print_mxmn('adv-flux-v',jk,z_adv_flux_v(:,:,:),n_zlev+1, &
        &              p_patch%nblks_c,'trc',ipl_src)
      END DO


      !calculate vert diffusion impicit: result is stored in trac_out
      CALL tracer_diffusion_vert_impl_hom( p_patch,         &
                                   & z_temp,                &
                                   & p_os%p_prog(nnew(1))%h,&
                                   & A_v,                   &
                                   & trac_out(:,:,:))
      IF (ldbg) THEN
        DO jk = 1, n_zlev
          write(*,*)'after impl-diff: max/min:',jk,&
          &maxval(trac_out(:,jk,:)), minval(trac_out(:,jk,:)) 
        END DO
      ENDIF


    !vertival diffusion is calculated explicitely
    ELSEIF(expl_vertical_tracer_diff==0)THEN
    CALL finish("advect_diffuse_vertical", &
    &"xplicit vertical tracer mixing currently not supported -  change namlist option")
 !       CALL tracer_diffusion_vert_expl( p_patch,       &
!                                     & trac_in,        &! z_trac_c,      &
!                                     &  z_h,           &
!                                     &  bc_top_tracer, &
!                                     &  bc_bot_tracer, & 
!                                     &  A_v,           &
!                                     &  z_div_diff_v)
!       DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!         CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
!         DO jc = i_startidx_c, i_endidx_c
!           !interior: from one below surface to the ground
!           z_dolic = v_base%dolic_c(jc,jb)
!           IF(z_dolic>=MIN_DOLIC)THEN
!             DO jk = 1, z_dolic
!                delta_z  = v_base%del_zlev_m(jk)
!                ! positive vertical divergence in direction of w (upward positive)
!                G_n_c_v(jc,jk,jb) = z_div_adv_v(jc,jk,jb)/z_h(jc,jk,jb) - z_div_diff_v(jc,jk,jb)
! 
!             END DO
!           ENDIF
!         END DO
!       END DO
! 
!       ipl_src = 5  ! output print level (1-5, fix)
!       IF (ldbg) THEN
!         DO jk = 1, n_zlev
!           CALL print_mxmn('diff flux_v',jk,z_diff_flux_v(:,:,:),n_zlev+1, &
!             &              p_patch%nblks_c,'trc',ipl_src)
!         END DO
!       ENDIF
!       DO jk = 1, n_zlev
!         CALL print_mxmn('div diff-flux_v',jk,z_div_diff_v(:,:,:),n_zlev, &
!         &              p_patch%nblks_c,'trc',ipl_src)
!       END DO
! 
!       IF( is_initial_timestep(timestep))THEN
!         G_nimd_c_v(:,:,:) = G_n_c_v(:,:,:)
!       ELSE
!         G_nimd_c_v(:,:,:) = (1.5_wp+AB_const)* G_n_c_v(:,:,:)   &
!           &               - (0.5_wp+AB_const)*G_nm1_c_v(:,:,:)
!       ENDIF
! 
! 
!       !Add advective and diffusive part to old tracer
!       DO jb = cells_in_domain%start_block, cells_in_domain%end_block
!         CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
!         DO jc = i_startidx_c, i_endidx_c
!           z_dolic = v_base%dolic_c(jc,jb)
!           IF(z_dolic>=MIN_DOLIC)THEN
!             !top  level
!             jk=1
!             delta_z = (p_os%p_prog(nold(1))%h(jc,jb)+v_base%del_zlev_m(jk))&
!                     &/z_h(jc,jk,jb) 
!     !         z_trac_c(jc,jk,jb) = (trac_in(jc,jk,jb)*delta_z      &
!     !                            & -delta_t*z_div_adv_v(jc,jk,jb)) &
!     !                            &/z_h(jc,jk,jb)
!     !          trac_out(jc,jk,jb) = z_trac_c(jc,jk,jb)               &
!     !                             & +delta_t*z_div_diff_v(jc,jk,jb)  
!     !         trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
!     !                            & -delta_t*(z_div_adv_v(jc,jk,jb) &
!     !                            &/z_h(jc,jk,jb)                  &
!     !                            &-z_div_diff_v(jc,jk,jb))
!             trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
!                                & -delta_t*G_nimd_c_v(jc,jk,jb) 
! 
!             !interior down to bottom
!             DO jk = 2, z_dolic
!               delta_z = v_base%del_zlev_m(jk)/z_h(jc,jk,jb)
!     !           z_trac_c(jc,jk,jb) = (trac_in(jc,jk,jb)*delta_z     &
!     !                              & -delta_t*z_div_adv_v(jc,jk,jb))  &
!     !                              &/z_h(jc,jk,jb) 
!     !          trac_out(jc,jk,jb) = z_trac_c(jc,jk,jb)                &
!     !                               & +delta_t*z_diff_flux_v(jc,jk,jb)
!     !           trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
!     !                              & -delta_t*(z_div_adv_v(jc,jk,jb) &
!     !                              &/z_h(jc,jk,jb)                  &
!     !                              &-z_div_diff_v(jc,jk,jb))
!               trac_out(jc,jk,jb) = trac_in(jc,jk,jb)*delta_z      &
!                                  & -delta_t*G_nimd_c_v(jc,jk,jb) 
!             END DO
!           ELSE
!             trac_out(jc,:,jb) = 0.0_wp
!           ENDIF
!         END DO
!       END DO
    ENDIF!(lvertical_diff_implicit)THEN
    IF (ltimer) CALL timer_stop(timer_dif_vert)

    CALL sync_patch_array(SYNC_C, p_patch, trac_out)

  END SUBROUTINE advect_diffuse_vertical
  !-------------------------------------------------------------------------
  !! First order upwind scheme for vertical tracer advection
  !!
  !!
  !! @par Revision History
  !! Seperated from vertical flux calculation
  !!
  !! mpi parallelized, no sync
  SUBROUTINE apply_tracer_flux_top_layer_oce( ppatch, pvar_c, pw_c,pupflux_i, tracer_id )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch           !< patch on which computation is performed
    REAL(wp), INTENT(IN   )           :: pvar_c(:,:,:)    !< advected cell centered variable
    REAL(wp), INTENT(IN   )           :: pw_c(:,:,:)      !< vertical velocity on cells 
    REAL(wp), INTENT(INOUT)           :: pupflux_i(:,:,:) !< variable in which the upwind flux is stored
                                                          !< dim: (nproma,n_zlev+1,nblks_c)
    INTEGER, INTENT(IN)               :: tracer_id
    ! local variables
    ! height based but reversed (downward increasing depth) coordinate system,
    ! grid coefficient is negative (same as pressure based atmospheric coordinate system
    REAL(wp), PARAMETER :: zcoeff_grid = -1.0_wp
    INTEGER             :: z_dolic
    INTEGER             :: i_startidx_c, i_endidx_c
    INTEGER             :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER             :: jkm1                     !< jk - 1
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------
    cells_in_domain => ppatch%cells%in_domain

    !fluxes at first layer
    !temperature has tracer_id=1 
    IF(tracer_id==1)THEN
      IF(temperature_relaxation/=0)THEN
        pupflux_i(:,1,:) = pvar_c(:,1,:)*pw_c(:,1,:)
      ELSEIF(temperature_relaxation==0)THEN
        pupflux_i(:,1,:) = 0.0
      ENDIF
    !salinity has tracer_id=2 
    ELSEIF(tracer_id==2)THEN
      IF(irelax_2d_S/=0)THEN
        pupflux_i(:,1,:) = pvar_c(:,1,:)*pw_c(:,1,:)
      ELSEIF(irelax_2d_S==0)THEN
        pupflux_i(:,1,:) = 0.0
      ENDIF
    ENDIF

  END SUBROUTINE apply_tracer_flux_top_layer_oce
  !-------------------------------------------------------------------------
  !! First order upwind scheme for vertical tracer advection
  !!
  !! Calculation of vertical tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !! Initial revision by Jochen Foerstner, DWD (2008-05-15)
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems
  !! Modification by Stephan Lorenz, MPI (2010-09-07)
  !! - adapted to hydrostatic ocean core
  !!
  !! mpi parallelized, no sync
  SUBROUTINE upwind_vflux_oce( ppatch, pvar_c, pw_c,top_bc_t, pupflux_i, tracer_id )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch           !< patch on which computation is performed
    REAL(wp), INTENT(IN   )           :: pvar_c(:,:,:)    !< advected cell centered variable
    REAL(wp), INTENT(IN   )           :: pw_c(:,:,:)      !< vertical velocity on cells 
    REAL(wp), INTENT(IN   )           :: top_bc_t(:,:)    !< top boundary condition traver!vertical velocity on cells
    REAL(wp), INTENT(INOUT)           :: pupflux_i(:,:,:) !< variable in which the upwind flux is stored
                                                          !< dim: (nproma,n_zlev+1,nblks_c)
    INTEGER, INTENT(IN)               :: tracer_id
    ! local variables
    ! height based but reversed (downward increasing depth) coordinate system,
    ! grid coefficient is negative (same as pressure based atmospheric coordinate system
    REAL(wp), PARAMETER :: zcoeff_grid = -1.0_wp
    INTEGER             :: z_dolic
    INTEGER             :: i_startidx_c, i_endidx_c
    INTEGER             :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER             :: jkm1                     !< jk - 1
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------
    cells_in_domain => ppatch%cells%in_domain

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        z_dolic = v_base%dolic_c(jc,jb)
        IF(z_dolic>=MIN_DOLIC)THEN
          DO jk = 2, z_dolic
            jkm1 = jk - 1
            ! calculate vertical tracer flux using upwind method
             pupflux_i(jc,jk,jb) =                 &
               &  laxfr_upflux_v( pw_c(jc,jk,jb),  &
               &                  pvar_c(jc,jkm1,jb), pvar_c(jc,jk,jb), zcoeff_grid )
          END DO
          !! no fluxes at bottom boundary
          !pupflux_i(jc,z_dolic+1,jb) = 0.0_wp
        ENDIF
      END DO
    END DO 

    CALL apply_tracer_flux_top_layer_oce( ppatch, pvar_c, pw_c, pupflux_i, tracer_id )

  END SUBROUTINE upwind_vflux_oce
  !-------------------------------------------------------------------------
  !>
  !! First order upwind scheme for vertical tracer advection
  !!
  !! Calculation of vertical tracer fluxes 
  !!
  !! @par Revision History
  !!!! mpi parallelized, no sync
  SUBROUTINE mimetic_vflux_oce( ppatch, pvar_c, pw_c, pupflux_i, tracer_id )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch           !< patch on which computation is performed
    REAL(wp), INTENT(IN)              :: pvar_c(:,:,:)    !< advected cell centered variable
    REAL(wp), INTENT(IN)              :: pw_c(:,:,:)      !< vertical velocity on cells
    REAL(wp), INTENT(INOUT)           :: pupflux_i(:,:,:) !< variable in which the upwind flux is stored
                                                          !< dim: (nproma,n_zlev+1,nblks_c)
    INTEGER, INTENT(IN)               :: tracer_id
    ! local variables
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: jkm1, jkp1, z_dolic                    !< jk - 1
    REAL(wp) :: w_ave(n_zlev)
    REAL(wp) :: w_avep1(n_zlev)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------
    cells_in_domain => ppatch%cells%in_domain

    w_ave(:)  = 0.0_wp
    w_avep1(:)= 0.0_wp

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_dolic = v_base%dolic_c(jc,jb)
          IF(z_dolic>=MIN_DOLIC)THEN
            DO jk = 2, z_dolic
              jkm1 = jk - 1
              jkp1 = jk + 1
              IF(pw_c(jc,jk,jb)>pw_c(jc,jkm1,jb))THEN

                w_ave(jk)=pw_c(jc,jk,jb)*pvar_c(jc,jk,jb)
              ELSE
                w_ave(jk)=pw_c(jc,jk,jb)*pvar_c(jc,jkm1,jb)
              ENDIF
             !w_avep1(jk)=0.5_wp*(pw_c(jc,jk,jb)+pw_c(jc,jkp1,jb))&
             !       &*pvar_c(jc,jk,jb)

            !pupflux_i(jc,jk,jb) = 0.5_wp*(w_ave(jk) +w_avep1(jk))
            !pupflux_i(jc,jk,jb) = (w_ave(jk)*v_base%del_zlev_m(jk) &
            !                    & +w_avep1(jk)*v_base%del_zlev_m(jkp1))&
            !                   & /(v_base%del_zlev_m(jk)+v_base%del_zlev_m(jkp1))
            END DO
            DO jk=1,z_dolic-1
              jkm1 = jk - 1
              jkp1 = jk + 1
              pupflux_i(jc,jk,jb) = (w_ave(jk)*v_base%del_zlev_m(jk) &
                                 & +w_avep1(jk)*v_base%del_zlev_m(jkp1))&
                                 & /(v_base%del_zlev_m(jk)+v_base%del_zlev_m(jkp1))
            END DO
            !w_ave(z_dolic)=0.5_wp*(pw_c(jc,z_dolic,jb)+pw_c(jc,z_dolic-1,jb))&
            !        &*pvar_c(jc,z_dolic,jb)
            !pupflux_i(jc,z_dolic,jb) = w_ave(z_dolic)*v_base%del_zlev_m(z_dolic)&
            !                          &/v_base%del_zlev_i(z_dolic)
            ! no fluxes at bottom boundary
            !pupflux_i(jc,z_dolic+1,jb) = 0.0_wp
        ENDIF
      END DO
    END DO
!     DO jb = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!         DO jc = i_startidx, i_endidx
!           DO jk = 2, n_zlev
!            ! index of top & bottom half level
!            jkm1 = jk - 1
!            jkp1 = jk + 1
!              w_ave(jk)=0.5_wp*(pw_c(jc,jk,jb)+pw_c(jc,jkm1,jb))&
!                    &*pvar_c(jc,jk,jb)
!            END DO ! end cell loop
!        END DO ! end level loop
!      END DO ! end block loop
!     DO jb = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!         DO jc = i_startidx, i_endidx
!          DO jk = 1, n_zlev-1
!            ! index of top & bottom half level
!            jkm1 = jk - 1
!            jkp1 = jk + 1
!             w_avep1(jk)=0.5_wp*(pw_c(jc,jkp1,jb)+pw_c(jc,jk,jb))&
!                    &*pvar_c(jc,jkp1,jb)
!            END DO ! end cell loop
!        END DO ! end level loop
!      END DO ! end block loop
!     DO jb = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!         DO jc = i_startidx, i_endidx 
!          DO jk = 1, n_zlev-1
!            ! index of top & bottom half level
!            jkm1 = jk - 1
!            jkp1 = jk + 1
!            pupflux_i(jc,jk,jb) = 0.5_wp*(w_ave(jk)  *v_base%del_zlev_m(jk) &
!                                &        +w_avep1(jk)*v_base%del_zlev_m(jkp1))/v_base%del_zlev_i(jk)
!            END DO ! end cell loop
!        END DO ! end level loop
!      END DO ! end block loop
    CALL apply_tracer_flux_top_layer_oce( ppatch, pvar_c, pw_c, pupflux_i, tracer_id )
  END SUBROUTINE mimetic_vflux_oce
  !-------------------------------------------------------------------------
  !>
  !!
  !! Calculation of central vertical tracer fluxes
  !!
  !! @par Revision History
  !! Petter Korn, MPI-M
  !!
  !! mpi parallelized, no sync
  SUBROUTINE central_vflux_oce_orig( ppatch, pvar_c, pw_c, pupflux_i, tracer_id )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch          !< patch on which computation is performed
    REAL(wp), INTENT(IN)              :: pvar_c(:,:,:)   !< advected cell centered variable
    REAL(wp), INTENT(IN)              :: pw_c(:,:,:)     !< vertical velocity on cells
    REAL(wp), INTENT(INOUT)           :: pupflux_i(:,:,:)!< variable in which the upwind flux is stored
                                                         !< dim: (nproma,n_zlev+1,nblks_c)
    INTEGER, INTENT(IN)                  :: tracer_id
    ! local variables
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: jkm1                     !< jk - 1
    INTEGER  :: z_dolic
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------
    cells_in_domain => ppatch%cells%in_domain


    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_dolic = v_base%dolic_c(jc,jb)
          IF(z_dolic>=MIN_DOLIC)THEN
            DO jk = 2, z_dolic
              ! index of top half level
              jkm1 = jk - 1

              ! calculate vertical tracer flux using upwind method
              pupflux_i(jc,jk,jb) = 0.5_wp*pw_c(jc,jk,jb) &
                &  * (pvar_c(jc,jkm1,jb)+ pvar_c(jc,jk,jb) )
          END DO
          ! no fluxes at bottom boundary
          pupflux_i(jc,z_dolic+1,jb) = 0.0_wp
        ENDIF
      END DO
    END DO
    CALL apply_tracer_flux_top_layer_oce( ppatch, pvar_c, pw_c, pupflux_i, tracer_id )
  END SUBROUTINE central_vflux_oce_orig
  !-------------------------------------------------------------------------
  !>
  !!
  !! Calculation of central vertical tracer fluxes
  !!
  !! @par Revision History
  !! Petter Korn, MPI-M
  !!
  !! mpi parallelized, no sync
  SUBROUTINE central_vflux_oce( ppatch, pvar_c, pw_c, pupflux_i, tracer_id )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch      !< patch on which computation is performed
    REAL(wp), INTENT(IN   )  :: pvar_c(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(INOUT)  :: pw_c(:,:,:)        !< vertical velocity on cells
    REAL(wp), INTENT(INOUT)  :: pupflux_i(:,:,:)   !< variable in which the upwind flux is stored
                                                   !< dim: (nproma,n_zlev+1,nblks_c)
    INTEGER, INTENT(IN)                  :: tracer_id
    ! local variables
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: jkm1                     !< jk - 1
    INTEGER  :: z_dolic
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-------------------------------------------------------------------------
    cells_in_domain => ppatch%cells%in_domain

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          z_dolic = v_base%dolic_c(jc,jb)
          IF(z_dolic>=MIN_DOLIC)THEN
            DO jk = 2, z_dolic
              ! index of top half level
              jkm1 = jk - 1

              ! calculate vertical tracer flux using upwind method
              pupflux_i(jc,jk,jb) = 0.5_wp*pw_c(jc,jk,jb)             &
                &  * (pvar_c(jc,jkm1,jb)+ pvar_c(jc,jk,jb) )          &
                &          +0.5_wp*pw_c(jc,jk,jb)*pw_c(jc,jk,jb)*dtime&
                &          / v_base%del_zlev_i(jk)                    &
                &          *( pvar_c(jc,jkm1,jb) -pvar_c(jc,jk,jb))

          END DO 
          ! no fluxes at bottom boundary
          pupflux_i(jc,z_dolic+1,jb) = 0.0_wp
        ENDIF
      END DO
    END DO

    CALL apply_tracer_flux_top_layer_oce( ppatch, pvar_c, pw_c, pupflux_i, tracer_id )

  END SUBROUTINE central_vflux_oce
  !------------------------------------------------------------------------
  !! The third order PPM scheme
  !!
  !! Calculation of time averaged vertical tracer fluxes using the third
  !! order PPM scheme.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2009-08-12)
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems. Included
  !!   parameter coeff_grid, in order to apply the same code to either a
  !!   pressure based or height based vertical coordinate system.
  !! Modification by Daniel Reinert, DWD (2011-01-17)
  !! - added optional parameter opt_lout_edge which will provide the 
  !!   reconstructed 'edge' value.
  !!
  !! mpi parallelized, sync required: p_upflux
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201
  ! - Carpenter et al. (1989), MWR, 118, 586-612
  ! - Lin and Rood (1996), MWR, 124, 2046-2070
  !
  SUBROUTINE upwind_vflux_ppm( p_patch, p_cc, &! p_mflx_contra_v,  &
    &                      p_w, p_dtime, p_itype_vlimit,             &
    &                      p_cellhgt_mc_now, &
    &                      p_upflux, tracer_id)

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_advection_vflux: upwind_vflux_ppm'

    TYPE(t_patch), TARGET, INTENT(IN) ::  p_patch  !< patch on which computation is performed
    REAL(wp), INTENT(IN)    :: p_cc(:,:,:)    !< advected cell centered variable
    !INTEGER, INTENT(IN)  :: p_iubc_adv   !< selects upper boundary condition
    !REAL(wp), INTENT(IN) :: p_mflx_contra_v(:,:,:) !< contravariant vertical mass flux dim: (nproma,nlevp1,nblks_c)
    REAL(wp), INTENT(IN)    :: p_w(:,:,:)    !< contravariant vertical velocity
    REAL(wp), INTENT(IN)    :: p_dtime  !< time step
    REAL(wp), INTENT(IN)    :: p_cellhgt_mc_now(:,:,:)    !< layer thickness at cell center at time n
    REAL(wp), INTENT(INOUT) :: p_upflux(:,:,:)    !< output field, containing the tracer mass flux or the reconstructed edge value
    INTEGER, INTENT(IN)     :: p_itype_vlimit  !< parameter to select the limiter for vertical transport
    INTEGER, INTENT(IN)     :: tracer_id
!
!local variables
    REAL(wp) :: z_face(nproma,p_patch%nlevp1,p_patch%nblks_c)   !< face values of transported field
    REAL(wp) :: z_face_up(nproma,p_patch%nlev,p_patch%nblks_c)  !< face value (upper face)
    REAL(wp) :: z_face_low(nproma,p_patch%nlev,p_patch%nblks_c) !< face value (lower face)
    REAL(wp) :: z_lext_1(nproma,p_patch%nlevp1)                 !< linear extrapolation value 1 
    REAL(wp) :: z_lext_2(nproma,p_patch%nlevp1)                 !< linear extrapolation value 2
    REAL(wp) :: z_cfl_m(nproma,p_patch%nlevp1,p_patch%nblks_c)  !< CFL number (weta>0, w<0)
    REAL(wp) :: z_cfl_p(nproma,p_patch%nlevp1,p_patch%nblks_c)  !< CFL number (weta<0, w>0)
    REAL(wp) :: z_slope(nproma,p_patch%nlev,p_patch%nblks_c)    !< monotonized slope
    !REAL(wp) :: zparent_topflx(nproma,p_patch%nblks_c)          !< necessary, to make this routine
    REAL(wp) :: z_slope_u, z_slope_l   !< one-sided slopes
    REAL(wp) :: z_delta_m, z_delta_p   !< difference between lower and upper face value
                                       !< for weta >0 and weta <0
    REAL(wp) :: z_a11, z_a12           !< 1/6 * a6,i (see Colella and Woodward (1984))
    REAL(wp) :: z_weta_dt              !< weta times p_dtime

    INTEGER  :: slev, slevp1           !< vertical start level and start level +1
    INTEGER  :: nlev, nlevp1           !< number of full and half levels
    INTEGER  :: ikm1, ikp1, ikp1_ic, & !< vertical level minus and plus one, plus two
      &  ikp2
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: jc, jk, jb              !< index of cell, vertical level and block
    INTEGER  :: jg                      !< patch ID
    REAL(wp) :: coeff_grid              !< parameter which is used to make the vertical 
                                        !< advection scheme applicable to a height      
                                        !< based coordinate system (coeff_grid=-1)
    !INTEGER :: opt_rlstart
    !INTEGER :: opt_rlend
    LOGICAL  :: opt_lout_edge !< optional: output edge value (.TRUE.),
                              !< or the flux across the edge   !< (.FALSE./not specified)
    !REAL(wp) :: opt_topflx_tra(nproma,p_patch%nblks_c)  !< vertical tracer flux at upper boundary 
    INTEGER, PARAMETER :: islopel_vsm = 1
    !-----------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-----------------------------------------------------------------------
    cells_in_domain => p_patch%cells%in_domain

    slev          = 1
    slevp1        = 2
    opt_lout_edge = .FALSE.
    ! number of vertical levels
    nlev   = n_zlev
    nlevp1 = n_zlev+1

    ! get patch ID
    jg         = p_patch%id
    coeff_grid = -1.0_wp! for height based vertical coordinate system advection_config(jg)%coeff_grid
    ! number of child domains
    i_nchdom   = MAX(1,p_patch%n_childdom)

    ! advection is done with an upwind scheme and a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    ! 3 options:  standard without limiter
    !             standard with semi-monotone or monotone limiter
    !             special version with limiter which handles CFL >1
    !
    ! 1. Calculate Courant number for weta>0 (w<0) and weta<0 (w>0)
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ikm1,z_weta_dt,ikp1_ic,ikp1, &
!$OMP            z_slope_u,z_slope_l,ikp2)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)

      ! Courant number at top
      z_cfl_p(i_startidx:i_endidx,slev,jb)   = 0._wp
      z_cfl_m(i_startidx:i_endidx,slev,jb)   = 0._wp
      ! Courant number at bottom
      z_cfl_p(i_startidx:i_endidx,nlevp1,jb) = 0._wp
      z_cfl_m(i_startidx:i_endidx,nlevp1,jb) = 0._wp

      DO jk = slevp1, nlev
        ! index of top half level
        ikm1 = jk - 1
        DO jc = i_startidx, i_endidx
          ! Calculate local Courant number at half levels
          ! z_cfl_m for weta >0 (w <0)
          ! z_cfl_p for weta <0 (w >0)
          z_weta_dt = ABS(p_w(jc,jk,jb)) * p_dtime * v_base%wet_c(jc,jk,jb)

          !Code modification this was multiplication by inverse height
          ! #slo# division by p_cellhgt_mc_now=zero on dry points
          IF ( v_base%lsm_oce_c(jc,ikm1,jb) <= sea_boundary ) &
            & z_cfl_m(jc,jk,jb) = z_weta_dt / p_cellhgt_mc_now(jc,ikm1,jb)
          IF ( v_base%lsm_oce_c(jc,  jk,jb) <= sea_boundary ) &
            & z_cfl_p(jc,jk,jb) = z_weta_dt / p_cellhgt_mc_now(jc,jk,jb)
          !ELSE
          !  z_weta_dt         = 0.0_wp
          !  z_cfl_m(jc,jk,jb) = 0.0_wp
          !  z_cfl_p(jc,jk,jb) = 0.0_wp
          !ENDIF
        END DO ! end loop over cells
      ENDDO ! end loop over vertical levels
      !
      ! 2. Calculate monotonized slope
      !
      z_slope(i_startidx:i_endidx,slev,jb) = 0._wp

      DO jk = slevp1, nlev
        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1_ic = jk + 1
        ikp1    = MIN( ikp1_ic, nlev )

        DO jc = i_startidx, i_endidx
          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          z_slope_u = 2._wp * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb))
          z_slope_l = 2._wp * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))

          IF ((z_slope_u * z_slope_l) > 0._wp) THEN

            z_slope(jc,jk,jb) = ( p_cellhgt_mc_now(jc,jk,jb)                             &
              &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)            &
              &  + p_cellhgt_mc_now(jc,ikp1,jb)) )                                       &
              &  * ( (2._wp * p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)) &
              &  / (p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                   &
              &  + (p_cellhgt_mc_now(jc,jk,jb) + 2._wp * p_cellhgt_mc_now(jc,ikp1,jb))   &
              &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              &  * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb)) )

            z_slope(jc,jk,jb) = SIGN(                                            &
              &  MIN( ABS(z_slope(jc,jk,jb)), ABS(z_slope_u), ABS(z_slope_l) ),  &
              &    z_slope(jc,jk,jb))
          ELSE
            z_slope(jc,jk,jb) = 0._wp
          ENDIF
          ELSE
            z_slope(jc,jk,jb) = 0._wp
          ENDIF
        END DO ! end loop over cells
      END DO ! end loop over vertical levels

      !
      ! 3. reconstruct face values at vertical half-levels
      !
      ! Boundary values for two highest and lowest half-levels
      !
      ! for faces k=slevp1 and k=nlevp1-1 reconstructed face values are calculated by
      ! interpolating a quadratic (instead of quartic) polynomial through 3
      ! values of the indefinite integral A=\int_{\eta_{0}}^{\eta}q\,\mathrm{d}\eta
      !
      ! for faces k=slev and k=nlevp1 a zero gradient condition is assumed and the
      ! face values are set to the tracer values of the corresponding cell centers
      !
      DO jc = i_startidx, i_endidx
        IF ( v_base%lsm_oce_c(jc,slevp1,jb) <= sea_boundary ) THEN
        z_face(jc,slevp1,jb) = p_cc(jc,slev,jb)*(1._wp - (p_cellhgt_mc_now(jc,slev,jb)&
          &       / p_cellhgt_mc_now(jc,slevp1,jb))) + (p_cellhgt_mc_now(jc,slev,jb)  &
          &       /(p_cellhgt_mc_now(jc,slev,jb) + p_cellhgt_mc_now(jc,slevp1,jb)))   &
          &       * ((p_cellhgt_mc_now(jc,slev,jb) / p_cellhgt_mc_now(jc,slevp1,jb))  &
          &       * p_cc(jc,slev,jb) + p_cc(jc,slevp1,jb))
        ELSE
        z_face(jc,slevp1,jb) = 0.0_wp
        ENDIF 
        IF ( v_base%lsm_oce_c(jc,nlev,jb) <= sea_boundary ) THEN
        z_face(jc,nlev,jb) = p_cc(jc,nlev-1,jb)*( 1._wp                               &
          &       - (p_cellhgt_mc_now(jc,nlev-1,jb) / p_cellhgt_mc_now(jc,nlev,jb)))  &
          &       + (p_cellhgt_mc_now(jc,nlev-1,jb)/(p_cellhgt_mc_now(jc,nlev-1,jb)   &
          &       + p_cellhgt_mc_now(jc,nlev,jb))) * ((p_cellhgt_mc_now(jc,nlev-1,jb) &
          &       / p_cellhgt_mc_now(jc,nlev,jb)) * p_cc(jc,nlev-1,jb)                &
          &       + p_cc(jc,nlev,jb))
        ELSE
          z_face(jc,nlev,jb) = 0.0_wp
        ENDIF
        z_face(jc,slev,jb)   = p_cc(jc,slev,jb)
        z_face(jc,nlevp1,jb) = p_cc(jc,nlev,jb)

      ENDDO


      DO jk = slevp1, nlev-2
        ! index of top half level
        ikm1 = jk - 1
        ! index of bottom half level
        ikp1 = jk + 1
        ikp2 = jk + 2
        DO jc = i_startidx, i_endidx
          IF ( v_base%lsm_oce_c(jc,ikp2,jb) <= sea_boundary ) THEN
          z_face(jc,ikp1,jb) = p_cc(jc,jk,jb) + (p_cellhgt_mc_now(jc,jk,jb)           &
            &  / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb)))         &
            &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                  &
            &  + (1._wp/(p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)    &
            &  + p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,ikp2,jb)))        &
            &  * ( (2._wp * p_cellhgt_mc_now(jc,ikp1,jb) * p_cellhgt_mc_now(jc,jk,jb) &
            &  / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb)))         &
            &  * ( (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))        &
            &  / (2._wp*p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb))    &
            &  - (p_cellhgt_mc_now(jc,ikp2,jb) + p_cellhgt_mc_now(jc,ikp1,jb))        &
            &  / (2._wp*p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb)) )  &
            &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb)) - p_cellhgt_mc_now(jc,jk,jb)     &
            &  * z_slope(jc,ikp1,jb) * (p_cellhgt_mc_now(jc,ikm1,jb)                  &
            &  + p_cellhgt_mc_now(jc,jk,jb)) / (2._wp*p_cellhgt_mc_now(jc,jk,jb)      &
            &  + p_cellhgt_mc_now(jc,ikp1,jb)) + p_cellhgt_mc_now(jc,ikp1,jb)         &
            &  * z_slope(jc,jk,jb) * (p_cellhgt_mc_now(jc,ikp1,jb)                    &
            &  + p_cellhgt_mc_now(jc,ikp2,jb)) / (p_cellhgt_mc_now(jc,jk,jb)          &
            &  + 2._wp*p_cellhgt_mc_now(jc,ikp1,jb)) )
           ELSE
           z_face(jc,ikp1,jb) = 0.0_wp
           ENDIF
        END DO ! end loop over cells
      END DO ! end loop over vertical levels

    END DO
!$OMP END DO
!$OMP END PARALLEL

    !
    ! 4. Limitation of first guess parabola (which is based on z_face)
    ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
    ! the limitation procedure.
    ! Therefore 2 additional fields z_face_up and z_face_low are
    ! introduced.
    !
    IF (p_itype_vlimit == islopel_vsm) THEN
!      ! monotonic (mo) limiter
      IF (ltimer) CALL timer_start(timer_ppm_slim)
      CALL v_ppm_slimiter_mo( p_patch, p_cc, z_face, z_slope, &  !in
        &                   z_face_up, z_face_low)               !inout
      IF (ltimer) CALL timer_stop(timer_ppm_slim)
    ELSE
!      ! simply copy face values to 'face_up' and 'face_low' arrays
!$OMP PARALLEL
!$OMP DO PRIVATE(jk,ikp1,jb,i_startidx,i_endidx)
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
        DO jk = slev, nlev
          ! index of bottom half level
          ikp1 = jk + 1
          !IF ( v_base%lsm_oce_c(jc,ikp1,jb) <= sea_boundary ) THEN

          z_face_up(i_startidx:i_endidx,jk,jb)  = z_face(i_startidx:i_endidx,jk,jb)&
                                                &*v_base%wet_c(jc,jk,jb)
          z_face_low(i_startidx:i_endidx,jk,jb) = z_face(i_startidx:i_endidx,ikp1,jb)&
                                                &*v_base%wet_c(jc,ikp1,jb)
          !ELSE
          !z_face_up(i_startidx:i_endidx,jk,jb)  = 0.0_wp
          !z_face_low(i_startidx:i_endidx,jk,jb) = 0.0_wp
          !ENDIF
        ENDDO
      ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
    ENDIF


!$OMP PARALLEL
    ! 5b. extrapolation using piecewise parabolic approx. of the transported
    ! quantity to the edge and finally, calculation of the upwind fluxes
    !

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_lext_1,z_lext_2,ikm1,z_delta_m, &
!$OMP            z_delta_p,z_a11,z_a12)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)

      z_lext_1(i_startidx:i_endidx,slev)   = p_cc(i_startidx:i_endidx,slev,jb)
      z_lext_2(i_startidx:i_endidx,slev)   = p_cc(i_startidx:i_endidx,slev,jb)
      z_lext_1(i_startidx:i_endidx,nlevp1) = p_cc(i_startidx:i_endidx,nlev,jb)
      z_lext_2(i_startidx:i_endidx,nlevp1) = p_cc(i_startidx:i_endidx,nlev,jb)

      DO jk = slevp1, nlev
        ! index of top half level
        ikm1 = jk -1
        DO jc = i_startidx, i_endidx
        IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN           
          ! linear extrapolated values
          ! for the height based coordinate system multiplication by coeff_grid
          ! is not necessary due to compensating (-) signs.
          ! first (of cell above) (case of w < 0; weta > 0)
          z_delta_m = z_face_low(jc,ikm1,jb) - z_face_up(jc,ikm1,jb)
          z_a11     = p_cc(jc,ikm1,jb)                                  &
            &       - 0.5_wp * (z_face_low(jc,ikm1,jb) + z_face_up(jc,ikm1,jb))

          z_lext_1(jc,jk) = p_cc(jc,ikm1,jb)                            &
            &  + (0.5_wp * z_delta_m * (1._wp - z_cfl_m(jc,jk,jb)))     &
            &  - z_a11*(1._wp - 3._wp*z_cfl_m(jc,jk,jb)                 &
            &  + 2._wp*z_cfl_m(jc,jk,jb)*z_cfl_m(jc,jk,jb))

          ! second (of cell below) (case of w > 0; weta < 0)
          z_delta_p = z_face_low(jc,jk,jb) - z_face_up(jc,jk,jb)
          z_a12     = p_cc(jc,jk,jb)                                    &
            &       - 0.5_wp * (z_face_low(jc,jk,jb) + z_face_up(jc,jk,jb))

          z_lext_2(jc,jk) = p_cc(jc,jk,jb)                              &
            &  - (0.5_wp * z_delta_p * (1._wp - z_cfl_p(jc,jk,jb)))     &
            &  - z_a12*(1._wp - 3._wp*z_cfl_p(jc,jk,jb)                 &
            &  + 2._wp*z_cfl_p(jc,jk,jb)*z_cfl_p(jc,jk,jb))
          !
          ! calculate vertical tracer flux
          !
          p_upflux(jc,jk,jb) =                                  &
            &  laxfr_upflux_v( p_w(jc,jk,jb),       &
            &                z_lext_1(jc,jk), z_lext_2(jc,jk),  &
            &                coeff_grid )
        ELSE
          p_upflux(jc,jk,jb) = 0.0_wp
        ENDIF
        END DO ! end loop over cells
      ENDDO ! end loop over vertical levels

      !
      ! set upper and lower boundary condition
      !
      !p_upflux(:,slev,:)   = p_w(:,1,:)*p_cc(:,1,:)
      !p_upflux(:,nlevp1,:) = 0.0_wp
     CALL apply_tracer_flux_top_layer_oce( p_patch, p_cc, p_w, p_upflux,&
                                         & tracer_id )

    ENDDO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

      !
      ! 6. If desired, apply a flux limiter to limit computed fluxes.
      !    These flux limiters are based on work by Zalesak (1979)
      !
!         IF (p_itype_vlimit == ifluxl_vpd) THEN
!           ! positive-definite (pd) limiter
!           CALL vflx_limiter_pd_ha( p_patch, p_dtime, p_cc, p_upflux,    & !in,inout
!             &                 opt_rlstart=i_rlstart, opt_rlend=i_rlend, & !in
!             &                 opt_slev=slev                             ) !in
!         ENDIF

    CALL sync_patch_array(SYNC_C, p_patch, p_upflux)

  END SUBROUTINE upwind_vflux_ppm
 !-------------------------------------------------------------------------
  !>
  !! Limiter for PPM (3rd order) vertical advection (monotone version)
  !!
  !! Removes over- and undershoots in first guess parabola by resetting the
  !! upper or lower interface values.
  !! Avoids non-physical over/undershoots in advected fields.
  !!
  !! Note that this limiter was coded assuming a pressure based vertical
  !! coordinate system. Nevertheless this limiter works for a height based
  !! vertical system, too. This is due to a 'wrong' computation of z_delta
  !! in the case of a height based coordinate system (i.e. z_delta is
  !! implicity multiplied by -1)
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  !! @par Revision History
  !! Developed by Daniel Reinert, DWD (2010-02-04)
  !!
  !! mpi parallelized, only cells_in_domain are computed, no sync
  SUBROUTINE v_ppm_slimiter_mo( p_patch, p_cc, p_face, p_slope, p_face_up, p_face_low )

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch          !< patch on which computation is performed
    REAL(wp), INTENT(IN)              :: p_cc(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(IN)              :: p_face(:,:,:)    !< reconstructed face values of the advected field
    REAL(wp), INTENT(IN)              :: p_slope(:,:,:)   !< monotonized slope
    REAL(wp), INTENT(INOUT)           :: p_face_up(:,:,:) !< final face value (upper face, height based)
    REAL(wp), INTENT(INOUT)           :: p_face_low(:,:,:)!< final face value (lower face, height based)

    ! locals
    INTEGER  :: nlev                      !< number of full levels
    INTEGER  :: slev                      !< vertical start level
    INTEGER  :: jc, jk, jb                !< index of cell, vertical level and block
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: ikp1                      !< vertical level plus one
    !INTEGER  :: opt_slev, opt_rlstart, opt_rlend
    REAL(wp) :: z_delta                   !< lower minus upper face value
    REAL(wp) :: z_a6i                     !< curvature of parabola
    !-----------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: cells_in_domain
    !-----------------------------------------------------------------------
    cells_in_domain => p_patch%cells%in_domain

    ! check optional arguments
    slev = 1
    nlev = n_zlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx_c,i_endidx_c,ikp1,z_delta,z_a6i)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jk = slev, nlev
        ! index of bottom half level
        ikp1 = jk + 1

        DO jc = i_startidx_c, i_endidx_c
          z_delta   = p_face(jc,ikp1,jb) - p_face(jc,jk,jb)
          z_a6i     = 6._wp * (p_cc(jc,jk,jb)                           &
            &       - 0.5_wp * (p_face(jc,jk,jb) + p_face(jc,ikp1,jb)))

          IF ( p_slope(jc,jk,jb) == 0._wp) THEN
            p_face_up(jc,jk,jb)  = p_cc(jc,jk,jb)
            p_face_low(jc,jk,jb) = p_cc(jc,jk,jb)

          ELSE IF (z_delta * z_a6i > z_delta * z_delta) THEN
            p_face_up(jc,jk,jb)  = 3._wp*p_cc(jc,jk,jb) - 2._wp*p_face(jc,ikp1,jb)
            p_face_low(jc,jk,jb) = p_face(jc,ikp1,jb)

          ELSE IF (z_delta * z_a6i < -1._wp * (z_delta * z_delta)) THEN
            p_face_up(jc,jk,jb)  = p_face(jc,jk,jb)
            p_face_low(jc,jk,jb) = 3._wp*p_cc(jc,jk,jb) - 2._wp*p_face(jc,jk,jb)

          ELSE
            p_face_up(jc,jk,jb)  = p_face(jc,jk,jb)
            p_face_low(jc,jk,jb) = p_face(jc,ikp1,jb)
          ENDIF
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE v_ppm_slimiter_mo

END MODULE mo_oce_tracer_transport_vert
