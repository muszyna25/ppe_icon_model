!>
!! Contains the implementation of the horizontal tracer transport routines for the ICON ocean model.
!! This comprises horizontal advection and diffusion
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/01)
!!
!!  mpi parallelized LL
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
MODULE mo_oce_tracer_transport_horz
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
USE mo_math_utilities,            ONLY: t_cartesian_coordinates, gc2cc
USE mo_math_constants,            ONLY: dbl_eps
USE mo_impl_constants,            ONLY: sea_boundary, sea,&
  &                                     min_rlcell, min_rledge, min_rlcell,MIN_DOLIC
USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer, idisc_scheme,    &
  &                                     irelax_3d_T, relax_3d_mon_T, irelax_3d_S, relax_3d_mon_S,&
  &                                     ab_const, ab_gam, expl_vertical_tracer_diff,l_edge_based,&
  &                                     iswm_oce, UPWIND, CENTRAL,MIMETIC,MIMETIC_MIURA,&
  &                                     FLUX_CALCULATION_HORZ  
USE mo_physical_constants,        ONLY: tf
USE mo_parallel_config,           ONLY: nproma
USE mo_dynamics_config,           ONLY: nold, nnew
USE mo_run_config,                ONLY: dtime, ltimer
USE mo_timer,                     ONLY: timer_start, timer_stop, timer_adv_horz, timer_hflx_lim, &
  &                                     timer_dif_horz
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, v_base, is_initial_timestep
USE mo_model_domain,              ONLY: t_patch
USE mo_exception,                 ONLY: finish !, message_text, message
USE mo_loopindices,               ONLY: get_indices_c, get_indices_e !, get_indices_v
USE mo_oce_boundcond,             ONLY: top_bound_cond_tracer
USE mo_oce_physics
USE mo_sea_ice_types,             ONLY: t_sfc_flx
USE mo_scalar_product,            ONLY:  map_cell2edges_3D,map_edges2cell_3D
USE mo_oce_math_operators,        ONLY:&! div_oce,grad_fd_norm_oce,&
                                       &div_oce_3D, grad_fd_norm_oce_3D!, grad_fd_norm_oce_2d
USE mo_advection_utils,           ONLY: laxfr_upflux
USE mo_oce_diffusion,             ONLY: tracer_diffusion_horz
USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
USE mo_sync,                      ONLY: SYNC_C, SYNC_C1, SYNC_E, sync_patch_array, &
  &                                     sync_patch_array_mult
USE mo_mpi,                       ONLY: my_process_is_mpi_parallel

IMPLICIT NONE

PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'
CHARACTER(len=12)           :: str_module    = 'oceTracHorz '  ! Output of module for 1 line debug
INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

!
! PUBLIC INTERFACE
!
PUBLIC :: advect_diffuse_flux_horz
! Private implemenation
!

PRIVATE :: upwind_hflux_oce
PRIVATE :: central_hflux_oce
PRIVATE :: upwind_hflux_oce_mimetic
PRIVATE :: central_hflux_oce_mimetic
PRIVATE :: mimetic_miura_hflux_oce
PRIVATE :: hflx_limiter_oce_mo
PRIVATE :: elad


 INTEGER, PARAMETER  :: top=1
CONTAINS
!-----------------------------------------------------------------------
!>
!! !  SUBROUTINE advects horizontally the tracers present in the ocean model.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!!
  !!  mpi parallelized LL
SUBROUTINE advect_diffuse_flux_horz( p_patch,             &
                                   & trac_old,            &
                                   & p_os,                &
                                   & p_op_coeff,          &
                                   & K_h,                 &
                                   & flux_horz,           &
                                   & h_e,                 &
                                   & cell_thick_intermed_c)
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch
  REAL(wp)                          :: trac_old(1:nproma,1:n_zlev,1:p_patch%nblks_c)
  TYPE(t_hydro_ocean_state), TARGET :: p_os
  TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
  REAL(wp), INTENT(in)              :: K_h(:,:,:)         !horizontal mixing coeff
  REAL(wp), INTENT(out)             :: flux_horz(1:nproma,1:n_zlev,1:p_patch%nblks_c)
  REAL(wp), INTENT(in)              :: h_e(nproma,p_patch%nblks_e)
  REAL(wp), INTENT(inout)           :: cell_thick_intermed_c(nproma,n_zlev,p_patch%nblks_c)
  !
  !Local variables
  REAL(wp) :: delta_z
  INTEGER  :: i_startidx_c, i_endidx_c
  INTEGER  :: i_startidx_e, i_endidx_e
  INTEGER  :: jc, jk, jb, je
  !REAL(wp) :: max_flux(1:n_zlev),min_flux(1:n_zlev)
  REAL(wp) :: z_adv_flux_h(nproma,n_zlev,p_patch%nblks_e)  ! horizontal advective tracer flux
  REAL(wp) :: z_div_adv_h(nproma,n_zlev,p_patch%nblks_c)   ! horizontal tracer divergence
  REAL(wp) :: z_div_diff_h(nproma,n_zlev,p_patch%nblks_c)  ! horizontal tracer divergence
  REAL(wp) :: z_diff_flux_h(nproma,n_zlev,p_patch%nblks_e) ! horizontal diffusive tracer flux
  TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
  TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain  
  !CHARACTER(len=max_char_length), PARAMETER :: &
  !        & routine = ('mo_tracer_transport_horz:advect_diffuse_flux_horz')
  !-------------------------------------------------------------------------------
  edges_in_domain => p_patch%edges%in_domain
  cells_in_domain => p_patch%cells%in_domain

  z_adv_flux_h(1:nproma,1:n_zlev,1:p_patch%nblks_e) =0.0_wp
  z_div_adv_h(1:nproma,1:n_zlev,1:p_patch%nblks_c)  =0.0_wp
  z_div_diff_h(1:nproma,1:n_zlev,1:p_patch%nblks_c) =0.0_wp
  z_diff_flux_h(1:nproma,1:n_zlev,1:p_patch%nblks_e)=0.0_wp
  flux_horz(1:nproma,1:n_zlev,1:p_patch%nblks_c)    =0.0_wp

  z_vn_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(1)  =0.0_wp
  z_vn_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(2)  =0.0_wp
  z_vn_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(3)  =0.0_wp

  ! Initialize timer for horizontal advection
  IF (ltimer) CALL timer_start(timer_adv_horz)


CALL sync_patch_array(SYNC_E, p_patch, p_os%p_diag%vn_time_weighted)


  !Horizontal advection
  IF( l_edge_based)THEN

    SELECT CASE(FLUX_CALCULATION_HORZ)

    CASE(UPWIND)

      !upwind estimate of tracer flux
      CALL upwind_hflux_oce( p_patch,                      &
                           & trac_old,                     &
                           & p_os%p_diag%vn_time_weighted, &
                           & z_adv_flux_h )
    CASE(CENTRAL)

      !central estimate of tracer flux
      CALL central_hflux_oce( p_patch,                      &
                           &  trac_old,                     &
                           &  p_os%p_diag%vn_time_weighted, &
                           &  z_adv_flux_h )

    CASE(MIMETIC_MIURA)

      !MIURA-type upwind-biased estimate of tracer flux
      CALL mimetic_miura_hflux_oce( p_patch,                    &
                                & trac_old,                     &
                                & p_os%p_diag%vn_time_weighted, &
                                & p_os%p_diag%p_vn_mean,        &
                                & p_op_coeff,                   &
                                & z_adv_flux_h )
    CASE DEFAULT
      CALL finish('TRIM(advect_diffuse_flux_horz)',"This flux option is not supported")

    END SELECT

    !Multiply with height at edges
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      DO jk = 1, n_zlev
        delta_z = v_base%del_zlev_m(jk)
        DO je = i_startidx_e, i_endidx_e
          IF ( v_base%lsm_oce_e(je,jk,jb) == sea ) THEN
            IF(jk==1)THEN
              delta_z=v_base%del_zlev_m(jk) + p_os%p_diag%h_e(je,jb)
            ENDIF
            z_adv_flux_h(je,jk,jb) = z_adv_flux_h(je,jk,jb)*delta_z
          ENDIF
        END DO
      END DO
    END DO

  ELSEIF(.NOT. l_edge_based)THEN

    CALL map_edges2cell_3D( p_patch,                      &
                          & p_os%p_diag%vn_time_weighted, &
                          & p_op_coeff, z_vn_c)

    !Multiply with height at cells
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        DO jk = 1, n_zlev 
          IF (v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
            delta_z = v_base%del_zlev_m(jk)
            IF(jk==1)THEN
              delta_z=v_base%del_zlev_m(jk) + p_os%p_prog(nold(1))%h(jc,jb)
            ENDIF
            z_vn_c(jc,jk,jb)%x = trac_old(jc,jk,jb)*z_vn_c(jc,jk,jb)%x*delta_z
          ENDIF
        END DO
      END DO
    END DO
    CALL sync_patch_array(SYNC_C,p_patch,z_vn_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(1))
    CALL sync_patch_array(SYNC_C,p_patch,z_vn_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(2))
    CALL sync_patch_array(SYNC_C,p_patch,z_vn_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(3))

    SELECT CASE(FLUX_CALCULATION_HORZ)

      CASE(UPWIND)

        CALL upwind_hflux_oce_mimetic( p_patch,                &
                                &trac_old,                     &
                                & p_os%p_diag%vn_time_weighted,&
                                & z_vn_c,                      &
                                & p_op_coeff,                  &
                                & z_adv_flux_h)

      CASE(CENTRAL)

        CALL central_hflux_oce_mimetic( p_patch,               &
                                & trac_old,                    &
                                & p_os%p_diag%vn_time_weighted,&
                                & z_vn_c,                      &
                                & p_op_coeff,                  &
                                & z_adv_flux_h)

      CASE(MIMETIC)

        CALL map_cell2edges_3D( p_patch, &
                              & z_vn_c,  &
                              & z_adv_flux_h,&
                              & p_op_coeff)
    CASE DEFAULT
      CALL finish('TRIM(advect_diffuse_flux_horz)',"This flux option is not supported")

    END SELECT
  ENDIF!l_edge_bases

  CALL sync_patch_array(SYNC_E, p_patch, z_adv_flux_h)

  ! Stop timer for horizontal advection
  IF (ltimer) CALL timer_stop(timer_adv_horz)


  ! Start timer for horizontal flux limitation
  IF (ltimer) CALL timer_start(timer_hflx_lim)

  !Flux limiting process, dependent on tracer configuration
  IF( l_edge_based)THEN  

    IF(FLUX_CALCULATION_HORZ/=UPWIND)THEN

      ! DO jk=1,n_zlev
      ! max_flux(jk)=maxval(z_adv_flux_h(:,jk,:))
      ! min_flux(jk)=minval(z_adv_flux_h(:,jk,:))
      ! END DO
      CALL hflx_limiter_oce_mo( p_patch,                 &
                              & trac_old,                &
                              & p_os%p_diag%mass_flx_e,  &
                              & z_adv_flux_h,            &
                              & p_os%p_diag%depth_c,     &
                              & cell_thick_intermed_c,   &
                              & p_op_coeff,              &
                              & p_os%p_diag%vn_time_weighted)
      ! DO jk=1,n_zlev
      ! write(*,*)'BEFORE/AFTER MAX/MIN LIMITER',jk, &
      ! &max_flux(jk),maxval(z_adv_flux_h(:,jk,:)), min_flux(jk),minval(z_adv_flux_h(:,jk,:))
      ! END DO
    ENDIF

  ELSEIF( .NOT.l_edge_based)THEN  
    !DO jk=1,n_zlev
    !  max_flux(jk)=maxval(z_adv_flux_h(:,jk,:))
    !   min_flux(jk)=minval(z_adv_flux_h(:,jk,:))
    !END DO
    IF(FLUX_CALCULATION_HORZ/=UPWIND)THEN
      CALL hflx_limiter_oce_mo_mimetic( p_patch,           &
                                  & trac_old,               &
                                  & p_os%p_diag%mass_flx_e, &
                                  & z_adv_flux_h,           &
                                  & z_vn_c,                 &
                                  & p_op_coeff,             &
                                  & p_os%p_diag%vn_time_weighted) 
     ENDIF

!  DO jk=1,n_zlev
!  !write(*,*)'BEFORE/AFTER MAX/MIN LIMITER',jk, max_flux(jk),&
!  !&maxval(z_adv_flux_h(:,jk,:)), min_flux(jk),&
!  !& minval(z_adv_flux_h(:,jk,:))
!  END DO
!     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
!       CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
!       DO jk = 1, n_zlev
!         DO je = i_startidx_e, i_endidx_e
!           write(999,*)'adv-flux:old:new:',je,jk,jb,z_adv_flux_h2(je,jk,jb),z_adv_flux_h(je,jk,jb)
!         END DO
!       END DO
!     END DO

 ENDIF!l_edge_based

  IF (ltimer) CALL timer_stop(timer_hflx_lim)

  !The diffusion part: calculate horizontal diffusive flux
  IF (ltimer) CALL timer_start(timer_dif_horz)
  CALL tracer_diffusion_horz( p_patch,      &
                            & trac_old,     &
                            & p_os,         &
                            & K_h,          &
                            & z_diff_flux_h,&
                            & subset_range = edges_in_domain)
  IF (ltimer) CALL timer_stop(timer_dif_horz)



  !Calculate divergence of advective and diffusive fluxes
  CALL div_oce_3D( z_adv_flux_h, p_patch,p_op_coeff%div_coeff, z_div_adv_h,&
     & subset_range=cells_in_domain)

  CALL div_oce_3D( z_diff_flux_h, p_patch,p_op_coeff%div_coeff, z_div_diff_h, &
     & subset_range=cells_in_domain)


  !Final step: calculate sum of advective and diffusive horizontal fluxes
  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
    DO jk = 1, n_zlev
      delta_z = v_base%del_zlev_m(jk)

      DO jc = i_startidx_c, i_endidx_c
        IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          flux_horz(jc,jk,jb) = z_div_diff_h(jc,jk,jb)-z_div_adv_h(jc,jk,jb)
        ENDIF
      END DO
    END DO
  END DO

  CALL sync_patch_array(SYNC_C, p_patch, flux_horz)

END SUBROUTINE advect_diffuse_flux_horz
!-----------------------------------------------------------------------
!-------------------------------------------------------------------------------
  !>
  !! First order upwind scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !!
  !!  mpi parallelized LL
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE upwind_hflux_oce_mimetic( p_patch, pvar_c, pvn_e, z_flux_cc,p_op_coeff,&
                                     & pupflux_e, opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch                        !< patch on which computation is performed
    REAL(wp), INTENT(IN)               :: pvar_c(nproma,n_zlev,p_patch%nblks_c)      !< advected cell centered variable
    REAL(wp), INTENT(IN)               :: pvn_e(nproma,n_zlev,p_patch%nblks_e)       !< normal velocity on edges
    TYPE(t_cartesian_coordinates)      :: z_flux_cc(nproma,n_zlev,p_patch%nblks_c)
    TYPE(t_operator_coeff), INTENT(in) :: p_op_coeff
    !EAL(wp), INTENT(INOUT), OPTIONAL :: ph_e (:,:)                                  !< surface elevation on edges
    REAL(wp), INTENT(INOUT)            :: pupflux_e(nproma,n_zlev,p_patch%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL      :: opt_slev    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL      :: opt_elev    ! optional vertical end level
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc
    INTEGER  :: slev, elev
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: je, jk, jb, il_c, ib_c         !< index of edge, vert level, block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain

    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF
    iilc => p_patch%edges%cell_idx
    iibc => p_patch%edges%cell_blk
    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at leteral walls)
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    ! loop through all patch edges (and blocks)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)

#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN

            il_c=p_op_coeff%upwind_cell_idx(je,jk,jb)
            ib_c=p_op_coeff%upwind_cell_blk(je,jk,jb)

            pupflux_e(je,jk,jb)=&
            & DOT_PRODUCT(z_flux_cc(il_c,jk,ib_c)%x,p_patch%edges%primal_cart_normal(je,jb)%x)

            !    write(999,*)'upwind flux -h',je,jk,jb,z_pub_flux_e_up(je,jk,jb),pupflux_e(je,jk,jb)&
            !    &, pvn_e(je,jk,jb),&
            !    &pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks


  END SUBROUTINE upwind_hflux_oce_mimetic
  !-----------------------------------------------------------------------
  !>
  !! First order upwind scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !!
  !!  mpi parallelized LL
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE central_hflux_oce_mimetic( p_patch, pvar_c, pvn_e, z_flux_cc,p_op_coeff,&
                                     & pupflux_e, opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch                                    !< patch on which computation is performed
    REAL(wp), INTENT(IN)              :: pvar_c (nproma,n_zlev,p_patch%nblks_c)     !< advected cell centered variable
    REAL(wp), INTENT(IN)              :: pvn_e(nproma,n_zlev,p_patch%nblks_e)       !< normal velocity on edges
    TYPE(t_cartesian_coordinates)     :: z_flux_cc(nproma,n_zlev,p_patch%nblks_c)
    TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
    REAL(wp), INTENT(INOUT)           :: pupflux_e(nproma,n_zlev,p_patch%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL     :: opt_slev    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL     :: opt_elev    ! optional vertical end level
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER  :: slev, elev
    INTEGER  :: i_startidx, i_endidx
    REAL(wp) :: z_fac
    INTEGER  :: je, jk, jb      !< index of edge, vert level, block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_cartesian_coordinates):: flux_sum
    REAL(wp) :: orientation
    !-----------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF
    flux_sum%x= 0.0_wp
    !
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at leteral walls)
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    iilc => p_patch%edges%cell_idx
    iibc => p_patch%edges%cell_blk

    ! loop through all patch edges (and blocks)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN

            flux_sum%x=0.5_wp*(z_flux_cc(iilc(je,jb,1),jk,iibc(je,jb,1))%x&
                             &+z_flux_cc(iilc(je,jb,2),jk,iibc(je,jb,2))%x)

            pupflux_e(je,jk,jb)=DOT_PRODUCT(flux_sum%x,&
                               &p_patch%edges%primal_cart_normal(je,jb)%x)
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks

  END SUBROUTINE central_hflux_oce_mimetic

!-------------------------------------------------------------------------------
  !>
  !! First order upwind scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !! Developed by L.Bonaventura  (2004).
  !! Adapted to new grid structure by L. Bonaventura, MPI-M, August 2005.
  !! Modification by Daniel Reinert, DWD (2010-02-09)
  !! - transferred to separate subroutine
  !! Modification by Stephan Lorenz, MPI (2010-09-06)
  !! - adapted to hydrostatic ocean core
  !!
  !!  mpi parallelized LL
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE upwind_hflux_oce( p_patch, pvar_c, pvn_e, pupflux_e, opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch                                    !< patch on which computation is performed
    REAL(wp), INTENT(IN)              :: pvar_c(nproma,n_zlev,p_patch%nblks_c)      !< advected cell centered variable
    REAL(wp), INTENT(IN)              :: pvn_e(nproma,n_zlev,p_patch%nblks_e)       !< normal velocity on edges
    REAL(wp), INTENT(OUT)             :: pupflux_e(nproma,n_zlev,p_patch%nblks_e)   !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL     :: opt_slev    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL     :: opt_elev    ! optional vertical end level
    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER  :: slev, elev
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF
    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at leteral walls)
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    iilc => p_patch%edges%cell_idx
    iibc => p_patch%edges%cell_blk

    ! loop through all patch edges (and blocks)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)

#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            pupflux_e(je,jk,jb) =  &
            &  laxfr_upflux( pvn_e(je,jk,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            &                                 pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
!write(*,*)'upwind flux -h',je,jk,jb,pupflux_e(je,jk,jb), pvn_e(je,jk,jb),&
!&pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))
          ELSE
            pupflux_e(je,jk,jb) = 0.0_wp
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks

  END SUBROUTINE upwind_hflux_oce
  !-----------------------------------------------------------------------
  !>
  !! Central scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using central fluxes
  !!
  !! @par Revision History
  !! Peter korn, MPI-M, 2011
  !!
  !!  mpi parallelized LL
  !!  mpi note: the result is not synced. Should be done in the calling method if required
  SUBROUTINE central_hflux_oce( p_patch, pvar_c, pvn_e, pupflux_e )

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch                                  !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)           :: pvar_c(nproma,n_zlev,p_patch%nblks_c)    !< advected cell centered variable
    REAL(wp), INTENT(INOUT)           :: pvn_e(nproma,n_zlev,p_patch%nblks_e)     !< normal velocity on edges
    REAL(wp), INTENT(INOUT)           :: pupflux_e(nproma,n_zlev,p_patch%nblks_e) !< variable in which the upwind flux is stored

    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: je, jk, jb         
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain

    ! line and block indices of two neighboring cells
    iilc => p_patch%edges%cell_idx
    iibc => p_patch%edges%cell_blk

    ! loop through all patch edges (and blocks)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = 1, n_zlev
        DO je = i_startidx, i_endidx
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            pupflux_e(je,jk,jb) =  0.5_wp*pvn_e(je,jk,jb)             &
              &        *( pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1))      &
              &          +pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)))     !&
!               &          +0.5_wp*pvn_e(je,jk,jb)*pvn_e(je,jk,jb)*dtime&
!               &          * p_patch%edges%inv_dual_edge_length(je,jb)   &
!               &        *( pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))      &
!               &          -pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)))
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
  END SUBROUTINE central_hflux_oce
!-------------------------------------------------------------------------------

 !>
  !! First order scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using Miura'sapproach with first
  !! order reconstruction. Miuras approach is described in H. Miura, Monthly
  !! Weather Review 135, 2007, 4038-4044. The original method has been formulated for
  !! hexagons and is transfered to triangles and uses the Mimetic capabilities.
  !!
  !! @par Revision History
  !! Developed by P. Korn, (2011).
  !!  mpi parallelized, the result is NOT synced. Should be done in the calling method if required
  !!
  SUBROUTINE mimetic_miura_hflux_oce( p_patch, pvar_c, pvn_e, pvn_mean_cc,&
                                    & p_op_coeff, pflux_e,&
                                    & opt_slev, opt_elev )
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    REAL(wp), INTENT(INOUT)           :: pvar_c(nproma,n_zlev,p_patch%nblks_c)!< advected cell centered variable
    REAL(wp), INTENT(INOUT)           :: pvn_e(nproma,n_zlev,p_patch%nblks_e) !< normal velocity on edges
    TYPE(t_cartesian_coordinates)     :: pvn_mean_cc(nproma,n_zlev,p_patch%nblks_c)
    TYPE(t_operator_coeff)            :: p_op_coeff
    REAL(wp), INTENT(INOUT)           :: pflux_e(nproma,n_zlev,p_patch%nblks_e)!< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL     :: opt_slev      ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL     :: opt_elev      ! optional vertical end level

    ! local variables
    INTEGER  :: slev, elev
    INTEGER  :: i_startidx_e, i_endidx_e
    INTEGER  :: je, jk, jb
    !INTEGER  :: il_v1, il_v2, ib_v1, ib_v2
    INTEGER  :: il_c, ib_c
    REAL(wp)                      :: z_gradC(nproma,n_zlev,p_patch%nblks_e)
    TYPE(t_cartesian_coordinates) :: z_gradC_cc(nproma,n_zlev,p_patch%nblks_c)
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-----------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain

    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF
!------------------------------------------------------------------
! ! !-------------upwind comparison
! !     iilc => p_patch%edges%cell_idx
! !     iibc => p_patch%edges%cell_blk
! !     ! loop through all patch edges (and blocks)
! !     DO jb = i_startblk_e, i_endblk_e
! !       CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e,   &
! !         &                i_startidx_e, i_endidx_e, rl_start,rl_end)
! !
! !       DO jk = slev, elev
! !         DO je = i_startidx_e, i_endidx_e          !
! !           IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
! !              pupflux_e(je,jk,jb) =  &
! !              &  laxfr_upflux( pvn_e(je,jk,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
! !              &                                 pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
! ! !write(*,*)'upwind miura',je,jk,jb,pupflux_e(je,jk,jb),pflux_e(je,jk,jb), pvn_e(je,jk,jb)!,&
! ! !&pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2))
! !           ELSE
! !             pupflux_e(je,jk,jb) = 0.0_wp
! !           ENDIF
! !         END DO  ! end loop over edges
! !       END DO  ! end loop over levels
! !     END DO  ! end loop over blocks
! ! !---------------------------------------------------------------

    !Step3: Local linear subgridcell distribution of tracer is calculated
    !3a: calculate tracer gradient
    !3b: map tracer gradient from edge to cell. Result is gradient vector at cell centers
    !3c: project gradient vector at cell center in direction of vector that points from
    !    cell center to upwind point C_i from step 2
    !3a:

    z_gradC(1:nproma,1:n_zlev,1:p_patch%nblks_e) = 0.0_wp

    z_gradC_cc(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(1) = 0.0_wp
    z_gradC_cc(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(2) = 0.0_wp
    z_gradC_cc(1:nproma,1:n_zlev,1:p_patch%nblks_c)%x(3) = 0.0_wp

    CALL grad_fd_norm_oce_3D( pvar_c,                 &
           &                  p_patch,                 &
           &                  p_op_coeff%grad_coeff,  &
           &                  z_gradC)

    CALL sync_patch_array(SYNC_E, p_patch, z_gradC)
    !3b:
    CALL map_edges2cell_3D( p_patch, z_gradC,p_op_coeff, z_gradC_cc)

    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
      DO jk = slev, elev
        DO je =  i_startidx_e, i_endidx_e
          IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary) THEN
            il_c=p_op_coeff%upwind_cell_idx(je,jk,jb)
            ib_c=p_op_coeff%upwind_cell_blk(je,jk,jb)

            pflux_e(je,jk,jb)=pvn_e(je,jk,jb)&
            &*(pvar_c(il_c,jk,ib_c)&
            &+dot_product(p_op_coeff%moved_edge_position_cc(je,jk,jb)%x  &
            &            -p_op_coeff%upwind_cell_position_cc(je,jk,jb)%x,&
            &             z_gradC_cc(il_c,jk,ib_c)%x))

          ENDIF 
        END DO
      END DO
    END DO

  END SUBROUTINE mimetic_miura_hflux_oce
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Flux limiter for horizontal advection
  !!
  !! Zalesak Flux-Limiter (Flux corrected transport)
  !! The corrected flux is a weighted average of the low order flux and the
  !! given high order flux. The high order flux is used to the greatest extent
  !! possible without introducing overshoots and undershoots.
  !! Note: This limiter is positive definite and almost monotone (but not strictly).
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2010-03-10)
  !! Modification by Daniel Reinert, DWD (2010-03-25)
  !! - adapted for MPI parallelization
  !!
  !!  mpi parallelized LL
  !!  mpi note: compute on domain edges. Results is not synced.
  !!
  SUBROUTINE hflx_limiter_oce_mo( p_patch,       &
                                & p_cc,            &
                                & p_mass_flx_e,    &
                                & p_mflx_tracer_h, &
                                & p_thick_old,     &
                                & p_thick_new,     &
                                & p_op_coeff,      &
                                & vn_time_weighted,&
                                & opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch                                        !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)           :: p_cc(nproma,n_zlev,p_patch%nblks_c)            !< advected cell centered variable
    REAL(wp), INTENT(inout)           :: p_mass_flx_e(nproma,n_zlev,p_patch%nblks_e)    !< horizontal mass flux
    REAL(wp), INTENT(INOUT)           :: p_mflx_tracer_h(nproma,n_zlev,p_patch%nblks_e) !< calculated horizontal tracer mass flux
    REAL(wp), INTENT(INOUT)           :: p_thick_old(nproma,n_zlev,p_patch%nblks_c)
    REAL(wp), INTENT(INOUT)           ::  p_thick_new(nproma,n_zlev,p_patch%nblks_c)
    TYPE(t_operator_coeff),INTENT(IN) :: p_op_coeff
    REAL(wp), INTENT(INOUT)           :: vn_time_weighted(nproma,n_zlev,p_patch%nblks_e)
    INTEGER, INTENT(IN), OPTIONAL     :: opt_slev !< optional vertical start level
    INTEGER, INTENT(IN), OPTIONAL     :: opt_elev !< optional vertical end level


!Local variables
    REAL(wp) :: z_mflx_low(nproma,n_zlev,p_patch%nblks_e)      !< first order tracer mass flux
    REAL(wp) :: z_anti(nproma,n_zlev,p_patch%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)
    REAL(wp) :: z_mflx_anti(nproma,n_zlev,p_patch%nblks_c,3) 
    REAL(wp) :: z_fluxdiv_c(nproma,n_zlev,p_patch%nblks_c)     !< flux divergence at cell center
    REAL(wp) :: z_tracer_new_low(nproma,n_zlev,p_patch%nblks_c)!< new tracer field after hor. transport,
                                                               !< if the low order fluxes are used
    REAL(wp) :: z_tracer_max(nproma,n_zlev,p_patch%nblks_c)    !< local maximum of current tracer value and low
                                                               !< order update
    REAL(wp) :: z_tracer_min(nproma,n_zlev,p_patch%nblks_c)    !< local minimum of current tracer value and low
                                                               !< order update
    REAL(wp) :: r_p(nproma,n_zlev,p_patch%nblks_c)             !< fraction which must multiply all in/out fluxes
                                                               !< of cell jc to guarantee
    REAL(wp) :: r_m(nproma,n_zlev,p_patch%nblks_c)             !< no overshoot/undershoot
    REAL(wp) :: r_frac                                         !< computed minimum fraction which must multiply< the flux at the edge
    REAL(wp) :: z_min(nproma,n_zlev)                           !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_max(nproma,n_zlev)
    REAL(wp) :: z_signum                                       !< sign of antidiffusive velocity
    REAL(wp) :: p_p, p_m                                       !< sum of antidiffusive fluxes into and out of cell jc

    INTEGER, DIMENSION(:,:,:), POINTER ::  cell_of_edge_idx, cell_of_edge_blk   !< Pointer to line and block indices of two
                                                                                !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk  !< Pointer to line and block indices of three
                                                                                !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk    !< Pointer to line and block indices (array)
                                                                                !< of edges
    INTEGER  :: slev, elev             !< vertical start and end level
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: je, jk, jb, jc         !< index of edge, vert level, block, cell
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
  !-------------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain
    cells_in_domain => p_patch%cells%in_domain


    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF


    ! Set pointers to index-arrays
    ! line and block indices of two neighboring cells
    cell_of_edge_idx => p_patch%edges%cell_idx
    cell_of_edge_blk => p_patch%edges%cell_blk
    ! line and block indices of edges as seen from cells
    edge_of_cell_idx => p_patch%cells%edge_idx
    edge_of_cell_blk => p_patch%cells%edge_blk
    ! pointers to line and block indices of three neighbor cells
    neighbor_cell_idx => p_patch%cells%neighbor_idx
    neighbor_cell_blk => p_patch%cells%neighbor_blk
    !
    ! 1. Calculate low (first) order fluxes using the standard upwind scheme and the
    !    antidiffusive fluxes
    !    (not allowed to call upwind_hflux_up directly, due to circular dependency)
    z_tracer_new_low(1:nproma,1:n_zlev,1:p_patch%nblks_c) = 0.0_wp
    z_tracer_max(1:nproma,1:n_zlev,1:p_patch%nblks_c)     = 0.0_wp
    z_tracer_min(1:nproma,1:n_zlev,1:p_patch%nblks_c)     = 0.0_wp
    r_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)              = 0.0_wp
    r_p(1:nproma,1:n_zlev,1:p_patch%nblks_c)              = 0.0_wp
    z_fluxdiv_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)      = 0.0_wp
    z_anti(1:nproma,1:n_zlev,1:p_patch%nblks_e)           = 0.0_wp
    z_mflx_low(1:nproma,1:n_zlev,1:p_patch%nblks_e)       = 0.0_wp


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=5
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif
          IF( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
!             z_mflx_low(je,jk,jb) = v_base%wet_e(je,jk,jb)* &
!             & laxfr_upflux( vn_time_weighted(je,jk,jb), &
!             & p_cc(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)), &
!             & p_cc(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)) )
             z_mflx_low(je,jk,jb) = &
             & laxfr_upflux( p_mass_flx_e(je,jk,jb), &
             & p_cc(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)), &
             & p_cc(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)) )

            ! calculate antidiffusive flux for each edge
            z_anti(je,jk,jb)     = p_mflx_tracer_h(je,jk,jb) - z_mflx_low(je,jk,jb)
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!$OMP END DO
   CALL sync_patch_array(SYNC_E, p_patch,z_mflx_low)
   CALL sync_patch_array(SYNC_E, p_patch, z_anti)
   


!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=4
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
          !
          ! 2. Define "antidiffusive" fluxes A(jc,jk,jb,je) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.
         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

           z_mflx_anti(jc,jk,jb,1) =                                                  &
             &     dtime * p_op_coeff%div_coeff(jc,jk,jb,1) / v_base%del_zlev_m(jk)   &!p_thick_new(jc,jk,jb)  &
             &   * z_anti(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1))

           z_mflx_anti(jc,jk,jb,2) =                                                  &
             &     dtime *  p_op_coeff%div_coeff(jc,jk,jb,2) /v_base%del_zlev_m(jk)   &!p_thick_new(jc,jk,jb)  &
             &   * z_anti(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2))

           z_mflx_anti(jc,jk,jb,3) =                                                  &
             &     dtime * p_op_coeff%div_coeff(jc,jk,jb,3) /v_base%del_zlev_m(jk)&  !  p_thick_new(jc,jk,jb)  &
             &   * z_anti(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3))

           !  compute also divergence of low order fluxes
           z_fluxdiv_c(jc,jk,jb) =  &
             & z_mflx_low(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1)) * &
             & p_op_coeff%div_coeff(jc,jk,jb,1) + &
             & z_mflx_low(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2)) * &
             & p_op_coeff%div_coeff(jc,jk,jb,2) + &
             & z_mflx_low(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3)) * &
             & p_op_coeff%div_coeff(jc,jk,jb,3)
         ENDIF
        ENDDO
      ENDDO
      !
      ! 3. Compute the updated low order solution z_tracer_new_low
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

            z_tracer_new_low(jc,jk,jb) = p_cc(jc,jk,jb)&
            & - dtime * z_fluxdiv_c(jc,jk,jb)/v_base%del_zlev_m(jk)!p_thick_new(jc,jk,jb)


          ! precalculate local maximum/minimum of current tracer value and low order
          ! updated value
          z_tracer_max(jc,jk,jb) = MAX(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
          z_tracer_min(jc,jk,jb) = MIN(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
          ELSE
            z_tracer_new_low(jc,jk,jb) = 0.0_wp
            z_tracer_max(jc,jk,jb)     = 0.0_wp
            z_tracer_min(jc,jk,jb)     = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.
    CALL sync_patch_array_mult(SYNC_C1, p_patch, 2, z_tracer_max, z_tracer_min)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_max,p_p,z_min,p_m)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=2
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          ! max value of cell and its neighbors
          ! also look back to previous time step
          z_max(jc,jk) = MAX( z_tracer_max(jc,jk,jb),                          &
            & z_tracer_max(neighbor_cell_idx(jc,jb,1),jk,neighbor_cell_blk(jc,jb,1)),  &
            & z_tracer_max(neighbor_cell_idx(jc,jb,2),jk,neighbor_cell_blk(jc,jb,2)),  &
            & z_tracer_max(neighbor_cell_idx(jc,jb,3),jk,neighbor_cell_blk(jc,jb,3)) )
          ! min value of cell and its neighbors
          ! also look back to previous time step
          z_min(jc,jk) = MIN( z_tracer_min(jc,jk,jb),                          &
            & z_tracer_min(neighbor_cell_idx(jc,jb,1),jk,neighbor_cell_blk(jc,jb,1)),  &
            & z_tracer_min(neighbor_cell_idx(jc,jb,2),jk,neighbor_cell_blk(jc,jb,2)),  &
            & z_tracer_min(neighbor_cell_idx(jc,jb,3),jk,neighbor_cell_blk(jc,jb,3)) )
          ENDIF
        ENDDO
      ENDDO

      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          ! Sum of all incoming antidiffusive fluxes into cell jc
          p_p =  -1._wp * (MIN(0._wp,z_mflx_anti(jc,jk,jb,1))   &
                         + MIN(0._wp,z_mflx_anti(jc,jk,jb,2))   &
                         + MIN(0._wp,z_mflx_anti(jc,jk,jb,3)) )
          ! Sum of all outgoing antidiffusive fluxes out of cell jc
          p_m =  MAX(0._wp,z_mflx_anti(jc,jk,jb,1))  &
            &  + MAX(0._wp,z_mflx_anti(jc,jk,jb,2))  &
            &  + MAX(0._wp,z_mflx_anti(jc,jk,jb,3))

          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of q
          r_m(jc,jk,jb) = (z_tracer_new_low(jc,jk,jb) - z_min(jc,jk))/(p_m + dbl_eps)
          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of q
          r_p(jc,jk,jb) = (z_max(jc,jk) - z_tracer_new_low(jc,jk,jb))/(p_p + dbl_eps)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(SYNC_C1, p_patch, 2, r_m, r_p)

    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
    !    program. Note that p_mflx_tracer_h now denotes the LIMITED flux.
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,r_frac,z_signum)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=3
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif
        IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
          z_signum = SIGN(1._wp,z_anti(je,jk,jb))

          ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
          ! but is computationally more efficient
          r_frac = 0.5_wp*( (1._wp+z_signum)*                &
             &     MIN(r_m(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)),  &
             &         r_p(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)))  &
             &     +  (1._wp-z_signum)*                      &
             &     MIN(r_m(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)),  &
             &         r_p(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)))  )

          ! Limited flux
          p_mflx_tracer_h(je,jk,jb) = z_mflx_low(je,jk,jb)               &
            &                       + MIN(1._wp,r_frac) * z_anti(je,jk,jb)
          ELSE
            p_mflx_tracer_h(je,jk,jb)= 0.0_wp
        ENDIF
        END DO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE hflx_limiter_oce_mo
  !-------------------------------------------------------------------------
  !>
  !! Flux limiter for horizontal advection
  !!
  !! Zalesak Flux-Limiter (Flux corrected transport)
  !! The corrected flux is a weighted average of the low order flux and the
  !! given high order flux. The high order flux is used to the greatest extent
  !! possible without introducing overshoots and undershoots.
  !! Note: This limiter is positive definite and almost monotone (but not strictly).
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !!
  !! @par Revision History
  !! - Inital revision by Daniel Reinert, DWD (2010-03-10)
  !! Modification by Daniel Reinert, DWD (2010-03-25)
  !! - adapted for MPI parallelization
  !!
  !!  mpi parallelized LL
  !!  mpi note: compute on domain edges. Results is not synced.
  !!
  SUBROUTINE hflx_limiter_oce_mo_mimetic( p_patch, p_cc, p_mass_flx_e, &
    &                         adv_tracer_flux_h, flux_cc, p_op_coeff,&
    &                         vn_time_weighted,&
    &                         opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    REAL(wp), INTENT(INOUT)           :: p_cc(1:nproma,1:n_zlev,1:p_patch%nblks_c)         !< advected cell centered variable
    REAL(wp), INTENT(inout)           :: p_mass_flx_e(1:nproma,1:n_zlev,1:p_patch%nblks_e) !horizontal mass flux (from dynamical core
    REAL(wp), INTENT(INOUT)           :: adv_tracer_flux_h(1:nproma,1:n_zlev,1:p_patch%nblks_e)  !< calculated horizontal tracer flux
    TYPE(t_cartesian_coordinates)     :: flux_cc(1:nproma,1:n_zlev,1:p_patch%nblks_c)
    TYPE(t_operator_coeff),INTENT(IN) :: p_op_coeff
    REAL(wp)                          :: vn_time_weighted(1:nproma,1:n_zlev,1:p_patch%nblks_e)
    INTEGER, INTENT(IN), OPTIONAL     :: opt_slev !< optional vertical start level
    INTEGER, INTENT(IN), OPTIONAL     :: opt_elev !< optional vertical end level
!     REAL(wp) :: p_thick_old(:,:,:) !< old surface elevationm, INTENT(in)
!     REAL(wp) ::  p_thick_new    !< new surface elevation INTENT(in)


!Local variables
    REAL(wp) :: z_mflx_low(nproma,n_zlev,p_patch%nblks_e)      !< first order tracer mass flux
    REAL(wp) :: z_anti(nproma,n_zlev,p_patch%nblks_e)          !< antidiffusive tracer mass flux (F_H - F_L)
    REAL(wp) :: z_mflx_anti(nproma,n_zlev,p_patch%nblks_c,3) 
    REAL(wp) :: z_fluxdiv_c(nproma,n_zlev,p_patch%nblks_c)     !< flux divergence at cell center
    REAL(wp) :: z_tracer_new_low(nproma,n_zlev,p_patch%nblks_c)!< new tracer field after hor. transport,
                                                               !< if the low order fluxes are used
    REAL(wp) :: z_tracer_max(nproma,n_zlev,p_patch%nblks_c)    !< local maximum of current tracer value and low
                                                               !< order update
    REAL(wp) :: z_tracer_min(nproma,n_zlev,p_patch%nblks_c)    !< local minimum of current tracer value and low
                                                               !< order update
    REAL(wp) :: r_p(nproma,n_zlev,p_patch%nblks_c)             !< fraction which must multiply all in/out fluxes
                                                               !< of cell jc to guarantee
    REAL(wp) :: r_m(nproma,n_zlev,p_patch%nblks_c)             !< no overshoot/undershoot
    REAL(wp) :: r_frac                                         !< computed minimum fraction which must multiply< the flux at the edge
    REAL(wp) :: z_min(nproma,n_zlev)                           !< minimum/maximum value in cell and neighboring cells
    REAL(wp) :: z_max(nproma,n_zlev)
    REAL(wp) :: z_signum                                       !< sign of antidiffusive velocity
    REAL(wp) :: p_p, p_m                                       !< sum of antidiffusive fluxes into and out of cell jc

    INTEGER, DIMENSION(:,:,:), POINTER ::  cell_of_edge_idx, cell_of_edge_blk   !< Pointer to line and block indices of two
                                                                                !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: neighbor_cell_idx, neighbor_cell_blk  !< Pointer to line and block indices of three
                                                                                !< neighbor cells (array)
    INTEGER, DIMENSION(:,:,:), POINTER :: edge_of_cell_idx, edge_of_cell_blk    !< Pointer to line and block indices (array)
                                                                                !< of edges
    INTEGER  :: slev, elev             !< vertical start and end level
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: je, jk, jb, jc,  il_c, ib_c        !< index of edge, vert level, block, cell
    TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
  !-------------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain
    cells_in_domain => p_patch%cells%in_domain

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = n_zlev
    END IF

    ! these should not be synced here
    CALL sync_patch_array(SYNC_C, p_patch, p_cc)
    CALL sync_patch_array(SYNC_E, p_patch, p_mass_flx_e)
    CALL sync_patch_array(SYNC_E, p_patch, adv_tracer_flux_h)

    ! Set pointers to index-arrays
    ! line and block indices of two neighboring cells
    cell_of_edge_idx => p_patch%edges%cell_idx
    cell_of_edge_blk => p_patch%edges%cell_blk
    ! line and block indices of edges as seen from cells
    edge_of_cell_idx => p_patch%cells%edge_idx
    edge_of_cell_blk => p_patch%cells%edge_blk
    ! pointers to line and block indices of three neighbor cells
    neighbor_cell_idx => p_patch%cells%neighbor_idx
    neighbor_cell_blk => p_patch%cells%neighbor_blk
    !
    ! 1. Calculate low (first) order fluxes using the standard upwind scheme and the
    !    antidiffusive fluxes
    !    (not allowed to call upwind_hflux_up directly, due to circular dependency)
      z_tracer_new_low(1:nproma,1:n_zlev,1:p_patch%nblks_c) = 0.0_wp
      z_tracer_max(1:nproma,1:n_zlev,1:p_patch%nblks_c)     = 0.0_wp
      z_tracer_min(1:nproma,1:n_zlev,1:p_patch%nblks_c)     = 0.0_wp
      r_m(1:nproma,1:n_zlev,1:p_patch%nblks_c)              = 0.0_wp
      r_p(1:nproma,1:n_zlev,1:p_patch%nblks_c)              = 0.0_wp
      z_fluxdiv_c(1:nproma,1:n_zlev,1:p_patch%nblks_c)      = 0.0_wp
      z_anti(1:nproma,1:n_zlev,1:p_patch%nblks_e)           = 0.0_wp
      z_mflx_low(1:nproma,1:n_zlev,1:p_patch%nblks_e)       = 0.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=5
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif
          IF( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
            !il_c=p_op_coeff%upwind_cell_idx(je,jk,jb)
            !ib_c=p_op_coeff%upwind_cell_blk(je,jk,jb)

            z_mflx_low(je,jk,jb) = &
            & laxfr_upflux( p_mass_flx_e(je,jk,jb), &
            & p_cc(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)), &
            & p_cc(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)) )


            !calculate antidiffusive flux for each edge
            z_anti(je,jk,jb) = adv_tracer_flux_h(je,jk,jb) - z_mflx_low(je,jk,jb)
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!$OMP END DO

   CALL sync_patch_array(SYNC_E, p_patch,z_mflx_low)
   CALL sync_patch_array(SYNC_E, p_patch, z_anti)



!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=4
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
          !
          ! 2. Define "antidiffusive" fluxes A(jc,jk,jb,je) for each cell. It is the difference
          !    between the high order fluxes (given by the FFSL-scheme) and the low order
          !    ones. Multiply with geometry factor to have units [kg/kg] and the correct sign.
          !    - positive for outgoing fluxes
          !    - negative for incoming fluxes
          !    this sign convention is related to the definition of the divergence operator.
         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

           z_mflx_anti(jc,jk,jb,1) =                                                  &
             &     dtime * p_op_coeff%div_coeff(jc,jk,jb,1) / v_base%del_zlev_m(jk)   &!p_thick_new(jc,jk,jb)  &
             &   * z_anti(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1))

           z_mflx_anti(jc,jk,jb,2) =                                                  &
             &     dtime *  p_op_coeff%div_coeff(jc,jk,jb,2) /v_base%del_zlev_m(jk)   &!p_thick_new(jc,jk,jb)  &
             &   * z_anti(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2))

           z_mflx_anti(jc,jk,jb,3) =                                                  &
             &     dtime * p_op_coeff%div_coeff(jc,jk,jb,3) /v_base%del_zlev_m(jk)&  !  p_thick_new(jc,jk,jb)  &
             &   * z_anti(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3))

           !  compute also divergence of low order fluxes
           z_fluxdiv_c(jc,jk,jb) =  &
             & z_mflx_low(edge_of_cell_idx(jc,jb,1),jk,edge_of_cell_blk(jc,jb,1)) * &
             & p_op_coeff%div_coeff(jc,jk,jb,1) + &
             & z_mflx_low(edge_of_cell_idx(jc,jb,2),jk,edge_of_cell_blk(jc,jb,2)) * &
             & p_op_coeff%div_coeff(jc,jk,jb,2) + &
             & z_mflx_low(edge_of_cell_idx(jc,jb,3),jk,edge_of_cell_blk(jc,jb,3)) * &
             & p_op_coeff%div_coeff(jc,jk,jb,3)
         ENDIF
        ENDDO
      ENDDO
      !
      ! 3. Compute the updated low order solution z_tracer_new_low
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

            z_tracer_new_low(jc,jk,jb) = p_cc(jc,jk,jb)&
            & - dtime * z_fluxdiv_c(jc,jk,jb)/v_base%del_zlev_m(jk)!p_thick_new(jc,jk,jb)

            ! precalculate local maximum/minimum of current tracer value and low order
            ! updated value
            z_tracer_max(jc,jk,jb) = MAX(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
            z_tracer_min(jc,jk,jb) = MIN(p_cc(jc,jk,jb),z_tracer_new_low(jc,jk,jb))
          ELSE
            z_tracer_new_low(jc,jk,jb) = 0.0_wp
            z_tracer_max(jc,jk,jb)     = 0.0_wp
            z_tracer_min(jc,jk,jb)     = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! 4. Limit the antidiffusive fluxes z_mflx_anti, such that the updated tracer
    !    field is free of any new extrema.
    CALL sync_patch_array_mult(SYNC_C1, p_patch, 2, z_tracer_max, z_tracer_min)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_max,p_p,z_min,p_m)
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=2
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          ! max value of cell and its neighbors
          ! also look back to previous time step
          z_max(jc,jk) = MAX( z_tracer_max(jc,jk,jb),                          &
            & z_tracer_max(neighbor_cell_idx(jc,jb,1),jk,neighbor_cell_blk(jc,jb,1)),  &
            & z_tracer_max(neighbor_cell_idx(jc,jb,2),jk,neighbor_cell_blk(jc,jb,2)),  &
            & z_tracer_max(neighbor_cell_idx(jc,jb,3),jk,neighbor_cell_blk(jc,jb,3)) )
          ! min value of cell and its neighbors
          ! also look back to previous time step
          z_min(jc,jk) = MIN( z_tracer_min(jc,jk,jb),                          &
            & z_tracer_min(neighbor_cell_idx(jc,jb,1),jk,neighbor_cell_blk(jc,jb,1)),  &
            & z_tracer_min(neighbor_cell_idx(jc,jb,2),jk,neighbor_cell_blk(jc,jb,2)),  &
            & z_tracer_min(neighbor_cell_idx(jc,jb,3),jk,neighbor_cell_blk(jc,jb,3)) )
          ENDIF
        ENDDO
      ENDDO

      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
          ! Sum of all incoming antidiffusive fluxes into cell jc
          p_p =  -1._wp * (MIN(0._wp,z_mflx_anti(jc,jk,jb,1))   &
                         + MIN(0._wp,z_mflx_anti(jc,jk,jb,2))   &
                         + MIN(0._wp,z_mflx_anti(jc,jk,jb,3)) )
          ! Sum of all outgoing antidiffusive fluxes out of cell jc
          p_m =  MAX(0._wp,z_mflx_anti(jc,jk,jb,1))  &
            &  + MAX(0._wp,z_mflx_anti(jc,jk,jb,2))  &
            &  + MAX(0._wp,z_mflx_anti(jc,jk,jb,3))

          ! fraction which must multiply all fluxes out of cell jc to guarantee no
          ! undershoot
          ! Nominator: maximum allowable decrease of q
          r_m(jc,jk,jb) = (z_tracer_new_low(jc,jk,jb) - z_min(jc,jk))/(p_m + dbl_eps)
          ! fraction which must multiply all fluxes into cell jc to guarantee no
          ! overshoot
          ! Nominator: maximum allowable increase of q
          r_p(jc,jk,jb) = (z_max(jc,jk) - z_tracer_new_low(jc,jk,jb))/(p_p + dbl_eps)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ! Synchronize r_m and r_p
    CALL sync_patch_array_mult(SYNC_C1, p_patch, 2, r_m, r_p)

    ! 5. Now loop over all edges and determine the minimum fraction which must
    !    multiply the antidiffusive flux at the edge.
    !    At the end, compute new, limited fluxes which are then passed to the main
    !    program. Note that p_mflx_tracer_h now denotes the LIMITED flux.
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,r_frac,z_signum)
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=3
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif
        IF ( v_base%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN
          z_signum = SIGN(1._wp,z_anti(je,jk,jb))

          ! This does the same as an IF (z_signum > 0) THEN ... ELSE ... ENDIF,
          ! but is computationally more efficient
          r_frac = 0.5_wp*( (1._wp+z_signum)*                &
             &     MIN(r_m(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)),  &
             &         r_p(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)))  &
             &     +  (1._wp-z_signum)*                      &
             &     MIN(r_m(cell_of_edge_idx(je,jb,2),jk,cell_of_edge_blk(je,jb,2)),  &
             &         r_p(cell_of_edge_idx(je,jb,1),jk,cell_of_edge_blk(je,jb,1)))  )

          ! Limited flux
          adv_tracer_flux_h(je,jk,jb) = z_mflx_low(je,jk,jb)               &
            &                       + MIN(1._wp,r_frac) * z_anti(je,jk,jb)
        ELSE
          adv_tracer_flux_h(je,jk,jb)= 0.0_wp
        ENDIF
        END DO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE hflx_limiter_oce_mo_mimetic
!-----------------------------------------------------------------------------------------
!>
!! !  SUBROUTINE advects horizontally the tracers present in the ocean model.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!!
SUBROUTINE elad(p_patch, trac_old, trac_new,p_op_coeff, p_os)
!
!
TYPE(t_patch), TARGET, INTENT(in) :: p_patch
REAL(wp), INTENT(IN)  :: trac_old(:,:,:)
REAL(wp), INTENT(OUT) :: trac_new(:,:,:)
TYPE(t_operator_coeff), INTENT(IN)        :: p_op_coeff
TYPE(t_hydro_ocean_state), TARGET :: p_os
!3 arrays for explicit part for tracer in Adams-Bashford  stepping,
!stores information across different timelevels
!REAL(wp) :: trac_out(:,:,:)                              !new tracer
!
!Local variables
INTEGER, PARAMETER :: no_cell_edges = 3
!REAL(wp) :: delta_z, delta_z2
!INTEGER  :: ctr, ctr_total
INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
INTEGER  :: jc, jk, jb!, jkp1        !< index of edge, vert level, block
INTEGER  :: z_dolic
REAL(wp) :: z_in(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_out(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_pred(nproma,n_zlev,p_patch%nblks_c)

REAL(wp) :: z_up(nproma,n_zlev,p_patch%nblks_c,no_cell_edges)
!REAL(wp) :: z_down(nproma,n_zlev,p_patch%nblks_c,no_cell_edges)
REAL(wp) :: z_max(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_min(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_excess(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_grad_excess(nproma,n_zlev,p_patch%nblks_e)
REAL(wp) :: z_diff_excess(nproma,n_zlev,p_patch%nblks_c)
REAL(wp) :: z_K(nproma,n_zlev,p_patch%nblks_e)
!REAL(wp) :: max_val, min_val!, dtime2, z_tmp
INTEGER  :: il_e(no_cell_edges), ib_e(no_cell_edges)
INTEGER  :: il_c1, ib_c1,il_c2, ib_c2, ie
INTEGER  :: stop_ctr
INTEGER            :: iter
INTEGER, PARAMETER :: i_max_iter = 20
!-----------------------------------------------------------------------
  IF (my_process_is_mpi_parallel()) &
    & CALL finish("elad", "is not mpi parallelized")

rl_start_c   = 1
rl_end_c     = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

z_in       = trac_old
z_out      = trac_new
z_pred     = trac_new
z_excess   = 0.0_wp
z_K(:,:,:) = 1.0E12_wp

! DO jk=1,n_zlev
!   write(*,*)'max/min tracer in:',jk,&
!   &maxval(trac_new(:,jk,:)),minval(trac_new(:,jk,:))
! END DO

  !Step 1: determine minimal and maximal permissible values
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)

    DO jc = i_startidx_c, i_endidx_c
      z_dolic = v_base%dolic_c(jc,jb)
      IF(z_dolic>=MIN_DOLIC)THEN
        DO jk = 1, z_dolic
          DO ie=1,no_cell_edges

            !actual edges of cell c1
            il_e(ie) = p_patch%cells%edge_idx(jc,jb,ie)
            ib_e(ie) = p_patch%cells%edge_blk(jc,jb,ie)

            !get neighbor cells of edge
            il_c1 = p_patch%edges%cell_idx(il_e(ie),ib_e(ie),1)
            ib_c1 = p_patch%edges%cell_blk(il_e(ie),ib_e(ie),1)
            il_c2 = p_patch%edges%cell_idx(il_e(ie),ib_e(ie),2)
            ib_c2 = p_patch%edges%cell_blk(il_e(ie),ib_e(ie),2)
!             IF(z_in(il_c1,jk,ib_c1)>z_in(il_c2,jk,ib_c2))THEN
!               z_up(jc,jk,jb,ie)   = z_in(il_c1,jk,ib_c1)
!               z_down(jc,jk,jb,ie) = z_in(il_c2,jk,ib_c2)
!             ELSE
!               z_up(jc,jk,jb,ie)   = z_in(il_c2,jk,ib_c2)
!               z_down(jc,jk,jb,ie) = z_in(il_c1,jk,ib_c1)
!             ENDIF
              IF(p_os%p_diag%ptp_vn(il_e(ie),jk,ib_e(ie))>=0.0_wp)THEN
                z_up(jc,jk,jb,ie) = trac_old(il_c1,jk,ib_c1)
              ELSEIF(p_os%p_diag%ptp_vn(il_e(ie),jk,ib_e(ie))<0.0_wp)THEN
                z_up(jc,jk,jb,ie) = trac_old(il_c2,jk,ib_c2)
              ENDIF
         END DO
         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
           z_max(jc,jk,jb) = maxval(z_up(jc,jk,jb,1:no_cell_edges))
           z_min(jc,jk,jb) = minval(z_up(jc,jk,jb,1:no_cell_edges))
         ELSE
           z_max(jc,jk,jb) = 0.0_wp
           z_min(jc,jk,jb) = 0.0_wp
         ENDIF
!         z_max(jc,jk,jb) = maxval(z_up(jc,jk,jb,1:no_cell_edges))
!         z_min(jc,jk,jb) = minval(z_down(jc,jk,jb,1:no_cell_edges))

        END DO!level-loop
      ENDIF!(z_dolic>0)
    END DO!idx-loop
  END DO !blk-loop

! DO jk=1,n_zlev
! write(*,*)'admissible bound',jk,&
! &maxval(z_max(:,jk,:)), minval(z_max(:,jk,:)),&
! &maxval(z_min(:,jk,:)), minval(z_min(:,jk,:))
! END DO


ITERATION_LOOP: Do iter = 1, i_max_iter
!write(*,*)'iteration----------',iter
!write(1230,*)'iteration----------',iter
  !Step 1: determine excess field
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)

    DO jc = i_startidx_c, i_endidx_c
      z_dolic = v_base%dolic_c(jc,jb)
      IF(z_dolic>MIN_DOLIC)THEN
        DO jk = 1, z_dolic
         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
           z_excess(jc,jk,jb) = max(z_pred(jc,jk,jb)-z_max(jc,jk,jb),0.0_wp)&
                              &+min(z_pred(jc,jk,jb)-z_min(jc,jk,jb),0.0_wp)
          ELSE
            z_excess(jc,jk,jb) = 0.0_wp
          ENDIF

!  IF(z_excess(jc,jk,jb)/=0.0)THEN
!  write(1230,*)'excess field',jc,jk,jb,z_excess(jc,jk,jb),&
! &z_pred(jc,jk,jb), z_max(jc,jk,jb),z_up(jc,jk,jb,1:no_cell_edges)
!  ENDIF
        END DO!level-loop
      ENDIF!(z_dolic>0)
    END DO!idx-loop
  END DO !blk-loop

! DO jk=1,n_zlev
! write(*,*)'max-min excess',jk,&
! &maxval(z_excess(:,jk,:)), minval(z_excess(:,jk,:))
! END DO


   !Step 3: Calculate diffusion of excess field
   !CALL grad_fd_norm_oce( z_excess, p_patch, z_grad_excess)
    CALL grad_fd_norm_oce_3D( z_excess,               &
           &                  p_patch,                &
           &                  p_op_coeff%grad_coeff,  &
           &                  z_grad_excess)
   z_grad_excess = z_K*z_grad_excess
   CALL div_oce_3D( z_grad_excess, p_patch,p_op_coeff%div_coeff, z_diff_excess)

! DO jk=1,n_zlev
! write(*,*)'max-min diffusion',jk,&
! &maxval(z_diff_excess(:,jk,:)), minval(z_diff_excess(:,jk,:))
! END DO

  !Step 4
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                   rl_start_c, rl_end_c)

    DO jc = i_startidx_c, i_endidx_c
      z_dolic = v_base%dolic_c(jc,jb)
      IF(z_dolic>=MIN_DOLIC)THEN
        DO jk = 1, z_dolic
          IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
            z_out(jc,jk,jb) = z_pred(jc,jk,jb) - z_excess(jc,jk,jb)+ z_diff_excess(jc,jk,jb)

!  IF(z_excess(jc,jk,jb)/=0.0)THEN
!  write(1230,*)'alg', z_out(jc,jk,jb),z_pred(jc,jk,jb),&
! &z_excess(jc,jk,jb), z_diff_excess(jc,jk,jb), &
! &z_max(jc,jk,jb),z_up(jc,jk,jb,1:no_cell_edges)
!  ENDIF

          ELSE
            z_out(jc,jk,jb) = 0.0_wp
          ENDIF
!           IF(z_diff_excess(jc,1,jb)/=0.0_wp)THEN
!             write(123,*)'correction',jc,jb,z_in(jc,1,jb),&
!             & z_out(jc,1,jb), z_diff_excess(jc,1,jb)
!           ENDIF
        END DO!level-loop
      ENDIF!(z_dolic>0)
    END DO!idx-loop
  END DO !blk-loop

  !Step 4: Stop criterion
  stop_ctr = 0
  DO jk=1,n_zlev

!      write(*,*)'actual state',jk,&
!       & maxval(z_out(:,jk,:)),minval(z_out(:,jk,:))
!      write(*,*)' pred',jk,&
!       & maxval(z_pred(:,jk,:)),minval(z_pred(:,jk,:))

    IF(   maxval(z_excess(:,jk,:))<1.0E-12_wp&
    &.AND.minval(z_excess(:,jk,:))<1.0E-12_wp)THEN

      stop_ctr = stop_ctr +1

    ELSE
    ENDIF
  END DO

  IF(stop_ctr==n_zlev)THEN
    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
      &                   rl_start_c, rl_end_c)

      DO jc = i_startidx_c, i_endidx_c
        z_dolic = v_base%dolic_c(jc,jb)
        IF(z_dolic>=MIN_DOLIC)THEN
          DO jk = 1, z_dolic
           IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
           trac_new(jc,jk,jb) = z_out(jc,jk,jb)
           ENDIF
!            write(1230,*)'before-after',jc,jk,jb,&
!           &trac_new(jc,jk,jb),trac_old(jc,jk,jb),z_excess(jc,jk,jb)
          END DO
        ENDIF
       ENDDO
     END DO
      exit ITERATION_LOOP
   ELSE
       z_pred(:,:,:) = z_out(:,:,:)
  ENDIF
END DO ITERATION_LOOP

!   DO jb = i_startblk_c, i_endblk_c
!     CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
!     &                   rl_start_c, rl_end_c)
!     DO jc = i_startidx_c, i_endidx_c
! write(123,*)'trac old new',trac_old(jc,1,jb),z_old_new(jc,1,jb),&
! &trac_new(jc,1,jb), z_out(jc,1,jb)
!
!    END DO
!  END DO
END SUBROUTINE elad
! !   !-------------------------------------------------------------------------
! !   !>
! !   !! Positive definite flux limiter for horizontal advection
! !   !!
! !   !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
! !   !! Only outward fluxes are re-scaled, in order to maintain positive
! !   !! definiteness.
! !   !!
! !   !! @par Literature:
! !   !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
! !   !!   Algorithms for Fluids. JCP, 31, 335-362
! !   !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
! !   !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
! !   !!   the cubed sphere grid. JCP, in press
! !   !!
! !   !! @par Revision History
! !   !! - Inital revision by Daniel Reinert, DWD (2010-10-06)
! !   !!
! !   SUBROUTINE hflx_limiter_sm( p_patch, p_dtime, p_cc,  p_op_coeff,  &
! !     &                         p_mflx_tracer_h, opt_slev, opt_elev )
! ! 
! !     TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
! !       &  p_patch
! ! 
! !     REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable
! !       &  p_cc(:,:,:)                 !< dim: (nproma,nlev,nblks_c)
! !                                      !< [kg kg^-1]
! !     TYPE(t_operator_coeff),INTENT(IN) :: p_op_coeff
! ! 
! !     REAL(wp), INTENT(IN) ::     &    !< time step [s]
! !       &  p_dtime
! ! 
! !     REAL(wp), INTENT(INOUT) ::  &    !< calculated horizontal tracer mass flux
! !       &  p_mflx_tracer_h(:,:,:)      !< dim: (nproma,nlev,nblks_e)
! !                                      !< [kg m^-2 s^-1]
! ! 
! ! !     REAL(wp), INTENT(IN), TARGET, OPTIONAL :: &!< density (\rho \Delta z)
! ! !       &  opt_rho(:,:,:)                !< dim: (nproma,nlev,nblks_c)
! ! !                                        !< [kg m^-2]
! ! 
! !     INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
! !       &  opt_slev
! ! 
! !     INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
! !       &  opt_elev
! ! 
! ! 
! !     !
! !     !Local variables
! !     !
! !     !< tracer mass flux ( total mass crossing the edge )
! !     REAL(wp) :: z_mflx(nproma,p_patch%nlev,p_patch%nblks_c,p_patch%cell_type) !< [kg m^-3]
! ! 
! !     !< fraction which must multiply all outgoing fluxes of cell jc to guarantee
! !     !< positive definiteness
! !     REAL(wp) :: r_m(nproma,p_patch%nlev,p_patch%nblks_c) 
! !     !< sign of mass flux 0: out; <0: in
! !     REAL(wp) :: z_signum  
! !     !< sum of fluxes out of cell < [kg m^-3]
! !     REAL(wp) :: p_m(nproma,p_patch%nlev,p_patch%nblks_c)                           
! !     REAL(wp) :: z_mflx_old(nproma,n_zlev,p_patch%nblks_e)
! !     INTEGER, DIMENSION(:,:,:), POINTER :: iilc, iibc
! !     INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
! ! 
! !     INTEGER  :: slev, elev                   
! !     INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
! !     INTEGER  :: i_rlstart, i_rlend, i_rlstart_c, i_rlend_c
! !     INTEGER  :: je, jk, jb, jc         
! !     TYPE(t_subset_range), POINTER :: edges_in_domain,  cells_in_domain
! !   !-------------------------------------------------------------------------
! !     ! Check for optional arguments
! !     IF ( PRESENT(opt_slev) ) THEN
! !       slev = opt_slev
! !     ELSE
! !       slev = 1
! !     END IF
! !     IF ( PRESENT(opt_elev) ) THEN
! !       elev = opt_elev
! !     ELSE
! !       elev = n_zlev
! !     END IF
! !     edges_in_domain=> p_patch%edges%in_domain
! !     cells_in_domain => p_patch%cells%in_domain
! ! 
! ! z_mflx_old=p_mflx_tracer_h
! ! 
! !     !
! !     ! Set pointers to index-arrays
! !     !
! !     ! line and block indices of two neighboring cells
! !     iilc => p_patch%edges%cell_idx
! !     iibc => p_patch%edges%cell_blk
! ! 
! !     ! line and block indices of edges as seen from cells
! !     iidx => p_patch%cells%edge_idx
! !     iblk => p_patch%cells%edge_blk
! ! 
! ! 
! !       ! 1. Reformulate all fluxes in terms of the total mass [kg m^-3]
! !       !    that crosses each of the CV-edges and store them in a cell-based structure.
! !       !    z_mflx > 0: outward
! !       !    z_mflx < 0: inward
! !       !
! ! 
! !    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !       CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
! ! #ifdef __LOOP_EXCHANGE
! !       DO jc = i_startidx, i_endidx
! !         DO jk = slev, elev
! ! #else
! ! !CDIR UNROLL=4
! !       DO jk = slev, elev
! !         DO jc = i_startidx, i_endidx
! ! #endif
! !             z_mflx(jc,jk,jb,1) = p_op_coeff%div_coeff(jc,jk,jb,1) * p_dtime &
! !               & * p_mflx_tracer_h(iidx(jc,jb,1),jk,iblk(jc,jb,1))
! ! 
! !             z_mflx(jc,jk,jb,2) = p_op_coeff%div_coeff(jc,jk,jb,2) * p_dtime &
! !               & * p_mflx_tracer_h(iidx(jc,jb,2),jk,iblk(jc,jb,2))
! ! 
! !             z_mflx(jc,jk,jb,3) = p_op_coeff%div_coeff(jc,jk,jb,3) * p_dtime &
! !               & * p_mflx_tracer_h(iidx(jc,jb,3),jk,iblk(jc,jb,3))
! ! 
! !         ENDDO
! !       ENDDO
! ! 
! !       !
! !       ! 2. Compute total outward mass
! !       !
! !       DO jk = slev, elev
! !         DO jc = i_startidx, i_endidx
! ! 
! !           ! Sum of all outgoing fluxes out of cell jc
! !           p_m(jc,jk,jb) =  MAX(0._wp,z_mflx(jc,jk,jb,1))  &
! !                       &  + MAX(0._wp,z_mflx(jc,jk,jb,2))  &
! !                       &  + MAX(0._wp,z_mflx(jc,jk,jb,3))
! ! 
! !           ! fraction which must multiply all fluxes out of cell jc to guarantee no
! !           ! undershoot
! !           ! Nominator: maximum allowable decrease of \rho q
! ! !           r_m(jc,jk,jb) = MIN(1._wp, (p_cc(jc,jk,jb)*ptr_rho(jc,jk,jb)) &
! ! !             &                        /(p_m + dbl_eps) )
! !            !r_m(jc,jk,jb) = MIN(1._wp, (p_cc(jc,jk,jb)*v_base%del_zlev_m(jk))&
! !            !                           &/(p_m(jc,jk,jb) + dbl_eps) )
! !           IF(p_cc(jc,jk,jb)/=0.0_wp)THEN
! !             IF(p_cc(jc,jk,jb)*v_base%del_zlev_m(jk)<=p_m(jc,jk,jb))THEN
! !               r_m(jc,jk,jb) = (p_cc(jc,jk,jb)*v_base%del_zlev_m(jk))/(p_m(jc,jk,jb) + dbl_eps)
! !             ELSE IF(p_cc(jc,jk,jb)*v_base%del_zlev_m(jk)>p_m(jc,jk,jb))THEN
! !               r_m(jc,jk,jb) = 1.0_wp
! !             ENDIF
! ! 
! !          ELSEIF(p_cc(jc,jk,jb)==0.0_wp)THEN
! !            r_m(jc,jk,jb) = 1.0_wp
! !          ENDIF
! ! ! IF(p_m(jc,jk,jb)/=0.0_wp.OR.p_cc(jc,jk,jb)/=0.0_wp)THEN
! ! ! write(123456,*)'in-out p_m',jc,jk,jb,r_m(jc,jk,jb),&
! ! ! &p_m(jc,jk,jb),p_cc(jc,jk,jb)*v_base%del_zlev_m(jk)
! ! ! ENDIF
! !         ENDDO
! !       ENDDO
! !     ENDDO
! ! !$OMP END DO
! ! 
! ! ! write(*,*)'max/min r_m', maxval(r_m), minval(r_m),&
! ! ! & maxval(p_m), minval(p_m),&
! ! ! & maxval(p_cc(:,1,:)*v_base%del_zlev_m(1)), minval(p_cc(:,1,:)*v_base%del_zlev_m(1))
! !     ! synchronize r_m
! !     CALL sync_patch_array(SYNC_C1,p_patch,r_m)
! ! 
! ! 
! !     !
! !     ! 3. Limit outward fluxes
! !     !    The inward ones remain untouched.
! !     !
! !     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
! ! #ifdef __LOOP_EXCHANGE
! !       DO je = i_startidx, i_endidx
! !         DO jk = slev, elev
! ! #else
! ! !CDIR UNROLL=3
! !       DO jk = slev, elev
! !         DO je = i_startidx, i_endidx
! ! #endif
! !             ! p_mflx_tracer_h > 0: flux directed from cell 1 -> 2
! !             ! p_mflx_tracer_h < 0: flux directed from cell 2 -> 1
! !             z_signum = SIGN(1._wp,p_mflx_tracer_h(je,jk,jb))
! ! 
! ! !             p_mflx_tracer_h(je,jk,jb) = p_mflx_tracer_h(je,jk,jb)                 &
! ! !             &* 0.5_wp *( (1._wp + z_signum) * r_m(iilc(je,jb,1),jk,iibc(je,jb,1)) &
! ! !                     &   +(1._wp - z_signum) * r_m(iilc(je,jb,2),jk,iibc(je,jb,2)) )
! !             IF(p_mflx_tracer_h(je,jk,jb)>=0.0_wp)THEN
! ! 
! !               p_mflx_tracer_h(je,jk,jb)&
! !               & = p_mflx_tracer_h(je,jk,jb)* r_m(iilc(je,jb,1),jk,iibc(je,jb,1))
! ! 
! !             ELSEIF(p_mflx_tracer_h(je,jk,jb)<0.0_wp)THEN
! ! 
! !               p_mflx_tracer_h(je,jk,jb)&
! !               & = p_mflx_tracer_h(je,jk,jb)* r_m(iilc(je,jb,2),jk,iibc(je,jb,2))
! ! 
! !             ENDIF
! ! !p_mflx_tracer_h(je,jk,jb)=p_mflx_tracer_h(je,jk,jb)*v_base%del_zlev_m(jk)
! ! 
! ! ! IF(p_mflx_tracer_h(je,jk,jb)/=0.0_wp.or.z_mflx_old(je,jk,jb)/=0.0_wp)THEN
! ! !  write(123456,*)'in-out',p_mflx_tracer_h(je,jk,jb),z_mflx_old(je,jk,jb),&
! ! !  &(1._wp + z_signum) * r_m(iilc(je,jb,1),jk,iibc(je,jb,1)) &
! ! !  &   +(1._wp - z_signum) * r_m(iilc(je,jb,2),jk,iibc(je,jb,2)),&
! ! !  & (1._wp + z_signum) * r_m(iilc(je,jb,1),jk,iibc(je,jb,1)), &
! ! !  &(1._wp - z_signum) * r_m(iilc(je,jb,2),jk,iibc(je,jb,2)),& 
! ! !  &r_m(iilc(je,jb,1),jk,iibc(je,jb,1)),r_m(iilc(je,jb,1),jk,iibc(je,jb,2))
! ! ! ENDIF
! !           ENDDO
! !         ENDDO
! !       ENDDO
! ! !$OMP END DO
! ! !$OMP END PARALLEL
! ! !write(*,*)'max dev', maxval(p_mflx_tracer_h-z_mflx_old), minval(p_mflx_tracer_h-z_mflx_old)
! ! 
! !   END SUBROUTINE hflx_limiter_sm
! !   !-------------------------------------------------------------------------
! ! !-----------------------------------------------------------------------
! ! !>
! ! !! !  SUBROUTINE advects horizontally the tracers present in the ocean model.
! ! !!
! ! !! @par Revision History
! ! !! Developed  by  Peter Korn, MPI-M (2011).
! ! !!
! !   !!  mpi parallelized LL
! ! SUBROUTINE advect_diffuse_horizontal(p_patch, trac_old,          &
! !                                    & p_os,p_op_coeff,            &
! !                                    & K_h,                        &
! !                                    & trac_new, timestep, delta_t,&
! !                                    & h_e, cell_thick_intermed_c, &
! !                                    & FLUX_CALCULATION_HORZ)
! !   TYPE(t_patch), TARGET, INTENT(in) :: p_patch
! !   REAL(wp)                          :: trac_old(:,:,:)
! !   TYPE(t_hydro_ocean_state), TARGET :: p_os
! !   TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
! !   REAL(wp), INTENT(in) :: K_h(:,:,:)         !horizontal mixing coeff
! !   REAL(wp), INTENT(out):: trac_new(:,:,:)
! !   INTEGER , INTENT(in) :: timestep           ! Actual timestep (to distinghuish initial step from others)
! !   REAL(wp), INTENT(in) :: delta_t
! !   REAL(wp), INTENT(in) :: h_e(nproma,p_patch%nblks_e)
! !   REAL(wp), INTENT(inout) :: cell_thick_intermed_c(nproma,n_zlev,p_patch%nblks_c)!intermediate cell thickness
! !   INTEGER , INTENT(in) :: FLUX_CALCULATION_HORZ
! !   !
! !   !Local variables
! !   REAL(wp) :: delta_z, delta_z2, delta_frac
! !   INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
! !   INTEGER  :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, rl_start_e, rl_end_e
! !   INTEGER  :: jc, jk, jb, je
! !   INTEGER  :: ilc1, ibc1,ilc2, ibc2
! !   REAL(wp) ::trac_content_1,trac_content_2
! !   !INTEGER  :: z_dolic
! !   REAL(wp) :: z_adv_flux_h(nproma,n_zlev,p_patch%nblks_e)  ! horizontal advective tracer flux
! !   REAL(wp) :: z_div_adv_h(nproma,n_zlev,p_patch%nblks_c)   ! horizontal tracer divergence
! !   REAL(wp) :: z_div_diff_h(nproma,n_zlev,p_patch%nblks_c)  ! horizontal tracer divergence
! !   REAL(wp) :: z_diff_flux_h(nproma,n_zlev,p_patch%nblks_e) ! horizontal diffusive tracer flux
! !   REAL(wp) :: z_trac_c(nproma,n_zlev, p_patch%nblks_c)
! !   REAL(wp) :: dummy_h_c2(nproma,n_zlev,p_patch%nblks_c)
! !   TYPE(t_cartesian_coordinates):: z_vn_c(nproma,n_zlev,p_patch%nblks_c)
! !   TYPE(t_subset_range), POINTER :: edges_in_domain, cells_in_domain
! !   !TYPE(t_cartesian_coordinates):: z_vn_c2(nproma,n_zlev,p_patch%nblks_c)
! !   ! CHARACTER(len=max_char_length), PARAMETER :: &
! !   !        & routine = ('mo_tracer_advection:advect_individual_tracer')
! !   !-------------------------------------------------------------------------------
! !   !z_tol= 1.0E-13
! !   !trac_old=10.0_wp
! !   edges_in_domain => p_patch%edges%in_domain
! !   cells_in_domain => p_patch%cells%in_domain
! ! 
! !   z_adv_flux_h  = 0.0_wp
! !   z_div_adv_h   = 0.0_wp
! !   z_div_diff_h  = 0.0_wp
! !   z_diff_flux_h = 0.0_wp
! ! !  z_trac_c      = 0.0_wp
! ! 
! ! ! !   !calculate (dummy) height consistent with divergence of mass fluxx (UPPERMOST LEVEL)
! !   jk = 1
! !   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !     CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !     DO jc = i_startidx_c, i_endidx_c
! !       IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! !         delta_z = v_base%del_zlev_m(jk)+p_os%p_prog(nold(1))%h(jc,jb)
! ! 
! !         p_os%p_diag%depth_c(jc,jk,jb) = v_base%del_zlev_m(jk)+p_os%p_prog(nold(1))%h(jc,jb)
! ! 
! !         cell_thick_intermed_c(jc,jk,jb)=v_base%del_zlev_m(jk)+p_os%p_prog(nnew(1))%h(jc,jb)
! !       ENDIF
! !     END DO
! !   END DO
! ! 
! !   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !     CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !     DO jk = 2, n_zlev
! !       delta_z = v_base%del_zlev_m(jk)
! !       DO jc = i_startidx_c, i_endidx_c
! !         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! !           p_os%p_diag%depth_c(jc,jk,jb) = delta_z
! !           cell_thick_intermed_c(jc,jk,jb)= delta_z
! !         ENDIF
! !       END DO
! !     END DO
! !   END DO
! ! 
! !   ! Initialize timer for horizontal advection
! !   IF (ltimer) CALL timer_start(timer_adv_horz)
! ! 
! ! 
! !   !Step 1) Horizontal advection and diffusion
! ! 
! !   SELECT CASE(FLUX_CALCULATION_HORZ)
! ! 
! !   CASE(UPWIND)!, MIMETIC)
! ! 
! !     !upwind estimate of tracer flux
! !     CALL upwind_hflux_oce( p_patch,        &
! !                          & trac_old,       &!      z_trac_c,       &
! !                          & p_os%p_diag%vn_time_weighted, &
! !                          & z_adv_flux_h )
! !   CASE(CENTRAL)!,MIMETIC)
! ! 
! !     !central estimate of tracer flux
! !     CALL central_hflux_oce( p_patch,       &
! !                          &  trac_old,      &
! !                          &  p_os%p_diag%vn_time_weighted, &
! !                          &  z_adv_flux_h )
! ! 
! ! !   CASE(MIMETIC)
! ! ! 
! ! !   CALL map_edges2cell( p_patch, p_os%p_diag%vn_time_weighted, z_vn_c)
! ! !   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! ! !     CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! ! !     DO jc = i_startidx_c, i_endidx_c
! ! !       DO jk = 1, n_zlev
! ! !           z_vn_c(jc,jk,jb)%x = trac_old(jc,jk,jb)*z_vn_c(jc,jk,jb)%x
! ! !       END DO
! ! !     END DO
! ! !   END DO
! ! !   CALL map_cell2edges( p_patch, z_vn_c, z_adv_flux_h)
! ! 
! !   CASE(MIMETIC_MIURA,MIMETIC)
! ! 
! !     !upwind estimate of tracer flux
! !     CALL mimetic_miura_hflux_oce( p_patch,                      &
! !                                 & trac_old,                     &
! !                                 & p_os%p_diag%vn_time_weighted, &
! !                                 & p_os%p_diag%p_vn_mean,        &
! !                                 & p_op_coeff,                   &
! !                                 & z_adv_flux_h )
! !   END SELECT
! ! !    DO jk = 1, n_zlev
! ! !    write(*,*)'max/min tracer flux horz',jk,&
! ! !    &maxval(z_adv_flux_h(:,jk,:)),minval(z_adv_flux_h(:,jk,:))
! ! !    END DO
! !   IF (ltimer) CALL timer_stop(timer_adv_horz)
! ! 
! !   IF(FLUX_CALCULATION_HORZ/=MIMETIC)THEN
! !     DO jb = edges_in_domain%start_block, edges_in_domain%end_block
! !       CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
! !       DO jk = 1, n_zlev
! !         delta_z = v_base%del_zlev_m(jk)
! !         DO je = i_startidx_e, i_endidx_e
! !           IF ( v_base%lsm_oce_e(je,jk,jb) == sea ) THEN
! !             IF(jk==1)THEN
! !               delta_z=v_base%del_zlev_m(jk) + p_os%p_diag%h_e(je,jb)
! !             ENDIF
! !             z_adv_flux_h(je,jk,jb) = z_adv_flux_h(je,jk,jb)*delta_z
! !           ENDIF
! !         END DO
! !       END DO
! !     END DO
! !   ENDIF
! ! 
! !   IF(FLUX_CALCULATION_HORZ==MIMETIC_MIURA)THEN! .OR. FLUX_CALCULATION_HORZ==MIMETIC)THEN
! !     IF (ltimer) CALL timer_start(timer_hflx_lim)
! ! 
! !     CALL sync_patch_array(SYNC_E, p_patch, z_adv_flux_h)
! !     CALL hflx_limiter_oce_mo( p_patch,                 &
! !                             & trac_old,                &
! !                             & p_os%p_diag%mass_flx_e,  &
! !                             & z_adv_flux_h,            &
! !                             & p_os%p_diag%depth_c,     &
! !                             & cell_thick_intermed_c,   &
! !                             & p_op_coeff,              &
! !                             & p_os%p_diag%vn_time_weighted)
! ! 
! !      IF (ltimer) CALL timer_stop(timer_hflx_lim)
! !   ENDIF
! ! 
! ! 
! !   !The diffusion part: calculate horizontal diffusive flux
! !   IF (ltimer) CALL timer_start(timer_dif_horz)
! !   CALL tracer_diffusion_horz( p_patch,      &
! !                             & trac_old,     &
! !                             & p_os,         &
! !                             & K_h,          &
! !                             & z_diff_flux_h,&
! !                             & subset_range = edges_in_domain)
! !   IF (ltimer) CALL timer_stop(timer_dif_horz)
! ! 
! !   ! LL: the z_diff_flux_h should not be synced in tracer_diffusion_horz
! ! 
! ! 
! !   !Calculate divergence of advective and diffusive fluxes
! !    CALL div_oce_3D( z_adv_flux_h, p_patch,p_op_coeff%div_coeff, z_div_adv_h,&
! !      & subset_range=cells_in_domain)
! ! 
! !    CALL div_oce_3D( z_diff_flux_h, p_patch,p_op_coeff%div_coeff, z_div_diff_h, &
! !      & subset_range=cells_in_domain)
! ! 
! ! 
! ! !---------------------------------------------------------------
! ! !  trac_content_1= 0.0_wp
! ! !   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! ! !     CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! ! !     DO jk = 1, n_zlev
! ! !       delta_z = v_base%del_zlev_m(jk)
! ! ! 
! ! !       DO jc = i_startidx_c, i_endidx_c
! ! !         IF(jk==1)delta_z = v_base%del_zlev_m(jk)+p_os%p_prog(nold(1))%h(jc,jb)
! ! !         !IF( v_base%wet_c(jc,jk,jb)>0.0_wp)THEN
! ! !         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! ! !          trac_content_1=trac_content_1+delta_z*p_patch%cells%area(jc,jb)*trac_old(jc,jk,jb)
! ! !        ENDIF
! ! !         !ENDIF
! ! !       END DO
! ! !     END DO
! ! !   END DO
! ! !---------------------------------------------------------------
! ! 
! !   !Final step: calculate new tracer values
! !     jk = 1
! !     DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !       CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !       DO jc = i_startidx_c, i_endidx_c
! ! 
! !         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! !           delta_z =v_base%del_zlev_m(jk)+p_os%p_prog(nold(1))%h(jc,jb)
! !           delta_z2=v_base%del_zlev_m(jk)+p_os%p_prog(nnew(1))%h(jc,jb)
! ! 
! ! !            trac_new(jc,jk,jb) &
! ! !            &= trac_old(jc,jk,jb)*(delta_z/delta_z2) &
! ! !            & -delta_t*((z_div_adv_h(jc,jk,jb)-z_div_diff_h(jc,jk,jb)*delta_z)&
! ! !            &/delta_z2)
! !            trac_new(jc,jk,jb) &
! !            &= trac_old(jc,jk,jb)*(delta_z/delta_z2) &
! !            & -(delta_t/delta_z2)*(z_div_adv_h(jc,jk,jb)-z_div_diff_h(jc,jk,jb))
! !         ENDIF
! !       END DO
! !     END DO
! ! 
! ! 
! !   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! !     CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! !     DO jk = 2, n_zlev
! !       delta_z = v_base%del_zlev_m(jk)
! ! 
! !       DO jc = i_startidx_c, i_endidx_c
! !         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! !           trac_new(jc,jk,jb) &
! !           &= trac_old(jc,jk,jb)&
! !           & -(delta_t/delta_z)*((z_div_adv_h(jc,jk,jb)-z_div_diff_h(jc,jk,jb)))
! !         ENDIF
! !       END DO
! !     END DO
! !   END DO
! ! !CALL elad(p_patch, trac_old, trac_new, p_os)
! !   CALL sync_patch_array(SYNC_C, p_patch, trac_new)
! ! !---------------------------------------------------------------
! ! !  trac_content_2= 0.0_wp
! ! !   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
! ! !     CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
! ! !     DO jk = 1, n_zlev
! ! !       delta_z = v_base%del_zlev_m(jk)
! ! ! 
! ! !       DO jc = i_startidx_c, i_endidx_c
! ! !         IF(jk==1)delta_z = v_base%del_zlev_m(jk)+p_os%p_prog(nnew(1))%h(jc,jb)
! ! !         !IF( v_base%wet_c(jc,jk,jb)>0.0_wp)THEN
! ! !         IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
! ! !          trac_content_2=trac_content_2+delta_z*p_patch%cells%area(jc,jb)*trac_new(jc,jk,jb)
! ! !        ENDIF
! ! !         !ENDIF
! ! !       END DO
! ! !     END DO
! ! !   END DO
! ! ! write(*,*)'trac:old:new:frac:',trac_content_1,trac_content_2, trac_content_1/trac_content_2
! ! !---------------------------------------------------------------
! ! 
! ! 
! ! END SUBROUTINE advect_diffuse_horizontal
!-----------------------------------------------------------------------

END MODULE mo_oce_tracer_transport_horz
