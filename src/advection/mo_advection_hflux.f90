!>
!! Computation of horizontal tracer flux
!!
!! This routine computes the upwind fluxes for all edges of a uniform
!! resolution patch. In case of different patches on a grid at the
!! same level in the grid hierarchy, no correction is needed
!! (the appropriate values must be then placed in the external halos).
!! The upwind flux function can be replaced by any other flux
!! function (e.g. Engquist-Osher, Godunov), in passive advection
!! case they are all equivalent. These routines compute only
!! the correct edge value of 'c*u' without multiplying by
!! the edge length, so that then the divergence operator
!! can be applied without modifications when computing the
!! flux divergence in the conservative transport formula.
!!
!! Possible options for horizontal tracer flux computation include
!! - first order Godunov method (UP1)
!! - second order MUSCL method
!! - MIURA with second order accurate reconstruction
!! - MIURA with third order accurate reconstruction
!!
!! For MUSCL the piecewise linear approximation is used:
!! See e.g. Lin et al. (1994), Mon. Wea. Rev., 122, 1575-1593
!! An improved, fully 2D version without splitting error is
!! implemeted, too.
!! See Miura, H. (2007), Mon. Wea. Rev., 135, 4038-4044
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Developed by L.Bonaventura  (2004).
!! Adapted to new grid structure by L. Bonaventura, MPI-M, August 2005.
!! Modification by Daniel Reinert, DWD (2009-08-09)
!! - implementation of second order MUSCL
!! Modification by Daniel Reinert, DWD (2009-11-09)
!! - implementation of second order MIURA
!! Modification by Daniel Reinert, DWD (2010-02-23)
!! - swapped slope limiter into new subroutine h_miura_slimiter_mo
!! Modification by Daniel Reinert, DWD (2010-04-23)
!! - moved slope and flux limiter to new module mo_advection_limiter
!! Modification by Daniel Reinert, DWD (2010-05-12)
!! - included MIURA scheme with third order accurate reconstruction
!! Modification by Daniel Reinert, DWD (2010-11-09)
!! - removed MUSCL-type computation of horizontal fluxes
!! Modification by Daniel Reinert, DWD (2011-09-20)
!! - new Miura-type advection scheme with internal time-step subcycling
!! Modification by Daniel Reinert, DWD (2012-05-04)
!! - removed optional slope limiter
!! Modification by Daniel Reinert, DWD (2013-09-30)
!! - new option: FFSL + Miura-type advection with subcycling
!! Modification by Will Sawyer, CSCS (2016-02-26)
!! - added OpenACC support
!! - added temp variables of type t_lsq to circumvent OpenACC compiler bug
!!
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
MODULE mo_advection_hflux

  USE mo_kind,                ONLY: wp, vp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, SUCCESS,                   &
    &                               min_rledge_int, min_rlcell_int,             &
    &                               UP, MIURA, MIURA3, FFSL, FFSL_HYB, MCYCL,   &
    &                               MIURA_MCYCL, MIURA3_MCYCL, FFSL_MCYCL,      &
    &                               FFSL_HYB_MCYCL, UP3, ifluxl_m, ifluxl_sm
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_math_gradients,      ONLY: grad_green_gauss_cell, grad_fe_cell
  USE mo_math_divrot,         ONLY: recon_lsq_cell_l, recon_lsq_cell_q,         &
    &                               recon_lsq_cell_cpoor, recon_lsq_cell_c,     &
    &                               recon_lsq_cell_l_svd, recon_lsq_cell_q_svd, &
    &                               recon_lsq_cell_cpoor_svd,                   &
    &                               recon_lsq_cell_c_svd, recon_lsq_cell_l_consv_svd
  USE mo_interpol_config,     ONLY: llsq_lin_consv, llsq_high_consv, lsq_high_ord, &
    &                               lsq_high_set
  USE mo_intp_data_strc,      ONLY: t_int_state, t_lsq
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_edge
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: timers_level
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c
  USE mo_sync,                ONLY: SYNC_C, SYNC_C1, sync_patch_array,          &
    &                               sync_patch_array_4de1
  USE mo_parallel_config,     ONLY: p_test_run
  USE mo_advection_config,    ONLY: advection_config, lcompute, lcleanup, t_trList
  USE mo_advection_utils,     ONLY: t_list2D
  USE mo_advection_quadrature,ONLY: prep_gauss_quadrature_l,                    &
    &                               prep_gauss_quadrature_l_list,               &
    &                               prep_gauss_quadrature_q,                    &
    &                               prep_gauss_quadrature_q_list,               &
    &                               prep_gauss_quadrature_c,                    &
    &                               prep_gauss_quadrature_c_list,               &
    &                               prep_gauss_quadrature_cpoor
  USE mo_advection_traj,      ONLY: btraj_dreg, t_back_traj,                    &
    &                               btraj_compute_o1, btraj_compute_o2
  USE mo_advection_geometry,  ONLY: divide_flux_area, divide_flux_area_list
  USE mo_advection_hlimit,    ONLY: hflx_limiter_mo, hflx_limiter_pd
  USE mo_timer,               ONLY: timer_adv_horz, timer_start, timer_stop
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_fortran_tools,       ONLY: init, copy
#ifdef _OPENACC
  USE mo_mpi,                 ONLY: i_am_accel_node, my_process_is_work
#endif


  IMPLICIT NONE

  PRIVATE



  PUBLIC :: hor_upwind_flux
  PUBLIC :: upwind_hflux_up
  PUBLIC :: upwind_hflux_miura
  PUBLIC :: upwind_hflux_miura3

#if defined( _OPENACC )
#if defined(__ADVECTION_HFLUX_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
  LOGICAL, PARAMETER ::  acc_validate = .FALSE.   ! ONLY SET TO .TRUE. FOR VALIDATION PHASE
#endif



  !-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Calculation of horizontal upwind flux at triangle edges on half levels
  !!
  !! Calculation of horizontal upwind flux at triangle edges on half levels
  !! using either
  !! - the first order Godunov method (UP1)
  !! - the MIURA method with linear reconstruction (essentially 2D)
  !! - the MIURA method with linear reconstruction and internal subcycling
  !! - the MIURA method with quadr/cubic reconstruction (essentially 2D)
  !! - the FFSL method with quadr/cubic reconstruction (essentially 2D)
  !!
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-02-10)
  !! Modification by Daniel Reinert (2010-11-05)
  !! - tracer loop moved from step_advection to hor_upwind_flux
  !! Modification by Daniel Reinert (2010-11-09)
  !! - removed MUSCL-type horizontal advection
  !! Modification by Daniel Reinert, DWD (2013-09-30)
  !! - new option: FFSL + Miura-type advection with subcycling
  !!
  !
  ! !LITERATURE
  ! MUSCL: Ahmad et al. (2006), Int. J. Num. Meth. Fluids, 50, 1247-1268
  !        Lin et al. (1994), MWR, 122, 1575-1593
  ! MIURA: Miura, H. (2007), Mon. Wea. Rev., 135, 4038-4044
  !
  SUBROUTINE hor_upwind_flux( p_cc, p_rhodz_now, p_rhodz_new, p_mass_flx_e, p_vn,   &
    &                     p_dtime, p_patch, p_int, p_ihadv_tracer, p_igrad_c_miura, &
    &                     p_itype_hlimit, p_iadv_slev, p_iord_backtraj, p_upflux,   &
    &                     opt_rlend )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_hflux: hor_upwind_flux'

    TYPE(t_patch), TARGET, INTENT(INOUT) ::  &     !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  & !< pointer to data structure for interpolation
      &  p_int

    REAL(wp),TARGET, INTENT(IN) ::  & !< advected cell centered variable
      &  p_cc(:,:,:,:)              !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(IN) ::     &   !< density times cell thickness at cell center step (n)
      &  p_rhodz_now(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &   !< density times cell thickness at cell center step (n+1)
      &  p_rhodz_new(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &   !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::     &   !< unweighted velocity field
      &  p_vn(:,:,:)                !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    INTEGER, INTENT(IN) ::      &   !< parameter to select numerical
      &  p_ihadv_tracer(:)          !< scheme for horizontal transport
                                    !< dim: (ntracer)

    INTEGER, INTENT(IN) ::      &   !< parameter to select the gradient
      &  p_igrad_c_miura            !< reconstruction method at cell center

    INTEGER, INTENT(IN) ::      &   !< parameter to select the limiter
      &  p_itype_hlimit(:)          !< for horizontal transport
                                    !< dim: (ntracer)

    INTEGER, INTENT(IN) ::      &   !< vertical start level for advection
      &  p_iadv_slev(:)             !< dim: (ntracer)

    INTEGER, INTENT(IN) ::      &   !< parameter to select the spacial order
      &  p_iord_backtraj            !< of accuracy for the backward trajectory

    REAL(wp), INTENT(INOUT) ::  &   !< variable in which the upwind flux is stored
      &  p_upflux(:,:,:,:)          !< dim: (nproma,nlev,nblks_e,ntracer)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    INTEGER :: jt, nt               !< tracer index and loop index
    INTEGER :: jg                   !< patch ID
    INTEGER :: i_rlend, i_rlend_vt, i_rlend_tr
    INTEGER :: i_rlstart
    INTEGER :: qvsubstep_elev       !< end level for qv-substepping
    INTEGER :: iadv_min_slev        !< scheme specific minimum slev
                                    !< i.e. minimum slev of all tracers which
                                    !< are advected with the given scheme
    INTEGER :: iadv_max_elev        !< scheme specific maximum elev
    INTEGER :: nsubsteps            !< number of substeps in miura_cycl (2 or 3)

    REAL(wp)::   &                  !< unweighted tangential velocity
      &  z_real_vt(nproma,p_patch%nlev,p_patch%nblks_e)!< component at edges

    TYPE(t_back_traj) :: btraj      !< backward trajectories for MIURA, MIURA_MCYCL
    TYPE(t_back_traj) :: btraj_cycl !< backward trajectories for subcycling

    TYPE(t_trList), POINTER :: &    !< pointer to tracer sublist
      &  trAdvect

    REAL(wp) :: z_dthalf            !< 0.5 * pdtime
    REAL(wp) :: z_dthalf_cycl       !< z_dthalf/nsubsteps

#ifdef _OPENACC
    LOGICAL  :: save_i_am_accel_node
#endif
    !-----------------------------------------------------------------------

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_real_vt
#endif
    ! get patch ID
    jg = p_patch%id

    ! tracer fields which are advected
    trAdvect => advection_config(jg)%trAdvect

    ! 
    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

    IF (timers_level > 2) CALL timer_start(timer_adv_horz)
    CALL btraj%construct(nproma,p_patch%nlev,p_patch%nblks_e,2)
    CALL btraj_cycl%construct(nproma,p_patch%nlev,p_patch%nblks_e,2)


!$ACC DATA  PCOPYIN( p_cc, p_mass_flx_e, p_rhodz_now, p_rhodz_new, p_vn ), PCOPYOUT( p_upflux ), &
!$ACC       CREATE( z_real_vt ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc, p_mass_flx_e, p_rhodz_now, p_rhodz_new, p_vn ), &
!$ACC IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

    !*******************************************************************
    !
    ! Tracer-independent part
    !
    !*******************************************************************

    ! In case that different transport schemes (MIURA, MIURA3) are used
    ! for different tracers, the double computation of tangential velocity
    ! vt should be avoided. Instead of computing vt inside each of the
    ! flux-routines, vt is computed only once per timestep prior to the flux
    ! routines. The resulting tangential velocity field is then passed to
    ! the flux routines.

    IF (ANY(p_ihadv_tracer(:)/= UP) .AND. ANY(p_ihadv_tracer(:)/= UP3)) THEN

      i_rlend_vt = MIN(i_rlend, min_rledge_int - 1)

      IF ( p_iord_backtraj /= 1 ) THEN
        i_rlend_vt = min_rledge_int - 3
      ENDIF

      ! reconstruct tangential velocity component at edge midpoints
      CALL rbf_vec_interpol_edge( p_vn, p_patch, p_int,            &! in
        &                         z_real_vt, opt_rlend=i_rlend_vt  )! inout
    ENDIF


    !
    ! Backward trajectory computation for MIURA-scheme with linear 
    ! reconstruction. In that case it is sufficient to compute 
    ! only the barycenter of the departure region (instead of all the vertices).
    i_rlstart  = 5
    i_rlend_tr = MIN(i_rlend, min_rledge_int - 1)
    qvsubstep_elev = advection_config(jg)%iadv_qvsubstep_elev


    IF (advection_config(jg)%isAnyTypeMiura) THEN

      iadv_min_slev = advection_config(jg)%miura_h%iadv_min_slev
      z_dthalf = 0.5_wp * p_dtime

      IF (p_iord_backtraj == 1)  THEN

        ! 1st order backward trajectory
        CALL btraj_compute_o1( btraj     = btraj,            & !inout
          &                  ptr_p       = p_patch,          & !in
          &                  ptr_int     = p_int,            & !in
          &                  p_vn        = p_vn,             & !in
          &                  p_vt        = z_real_vt,        & !in
          &                  p_dthalf    = z_dthalf,         & !in
          &                  opt_rlstart = i_rlstart,        & !in
          &                  opt_rlend   = i_rlend_tr,       & !in
          &                  opt_slev    = iadv_min_slev,    & !in
          &                  opt_elev    = p_patch%nlev      ) !in
      ELSE
        ! 2nd order backward trajectory
        CALL btraj_compute_o2 ( btraj       = btraj,            & !inout
          &                     ptr_p       = p_patch,          & !in
          &                     ptr_int     = p_int,            & !in
          &                     p_vn        = p_vn,             & !in
          &                     p_vt        = z_real_vt,        & !in
          &                     p_dthalf    = z_dthalf,         & !in
          &                     opt_rlstart = i_rlstart,        & !in
          &                     opt_rlend   = i_rlend_tr,       & !in
          &                     opt_slev    = iadv_min_slev,    & !in
          &                     opt_elev    = p_patch%nlev      ) !in
      ENDIF
    ENDIF


    IF (advection_config(jg)%isAnyTypeMcycl) THEN

      iadv_min_slev = advection_config(jg)%mcycl_h%iadv_min_slev
      ! should be moved to advection_config
      iadv_max_elev = MERGE(p_patch%nlev,qvsubstep_elev,ANY(p_ihadv_tracer(:)== MCYCL))

      ! Determine number of substeps in miura_cycl
      ! It is assumed that three substeps are needed if the top of the currently active
      ! model domain is higher than 40 km, otherwise, two are sufficient
      !
      IF (vct_a(iadv_min_slev+p_patch%nshift_total) > 40000._wp) THEN
        nsubsteps = 3
      ELSE
        nsubsteps = 2
      ENDIF
      z_dthalf_cycl = 0.5_wp * p_dtime/REAL(nsubsteps,wp)

      !
      IF (p_iord_backtraj == 1)  THEN
        !
        ! 1st order backward trajectory for subcycled version
        ! The only thing that differs is the time step passed in
        CALL btraj_compute_o1( btraj     = btraj_cycl,       & !inout
          &                  ptr_p       = p_patch,          & !in
          &                  ptr_int     = p_int,            & !in
          &                  p_vn        = p_vn,             & !in
          &                  p_vt        = z_real_vt,        & !in
          &                  p_dthalf    = z_dthalf_cycl,    & !in
          &                  opt_rlstart = i_rlstart,        & !in
          &                  opt_rlend   = i_rlend_tr,       & !in
          &                  opt_slev    = iadv_min_slev,    & !in
          &                  opt_elev    = iadv_max_elev     ) !in
      ELSE
        !
        ! 2nd order backward trajectory for subcycled version
        ! The only thing that differs is the time step passed in
        CALL btraj_compute_o2 ( btraj       = btraj_cycl,       & !inout
          &                     ptr_p       = p_patch,          & !in
          &                     ptr_int     = p_int,            & !in
          &                     p_vn        = p_vn,             & !in
          &                     p_vt        = z_real_vt,        & !in
          &                     p_dthalf    = z_dthalf_cycl,    & !in
          &                     opt_rlstart = i_rlstart,        & !in
          &                     opt_rlend   = i_rlend_tr,       & !in
          &                     opt_slev    = iadv_min_slev,    & !in
          &                     opt_elev    = iadv_max_elev     ) !in
      ENDIF
    ENDIF


    !*******************************************************************
    !
    ! Tracer-specific part
    !
    !*******************************************************************

    DO nt = 1, trAdvect%len ! Tracer loop

      jt = trAdvect%list(nt)

      ! Select desired flux calculation method
      SELECT CASE( p_ihadv_tracer(jt) )

      CASE( UP )      ! ihadv_tracer = 1
        ! CALL first order upwind
        CALL upwind_hflux_up( p_patch      = p_patch,            & !in
          &                   p_cc         = p_cc(:,:,:,jt),     & !in
          &                   p_mass_flx_e = p_mass_flx_e,       & !in
          &                   p_upflux     = p_upflux(:,:,:,jt), & !inout
          &                   opt_slev     = p_iadv_slev(jt),    & !in
          &                   opt_rlend    = i_rlend             ) !in


      CASE( MIURA )   ! ihadv_tracer = 2

        ! CALL MIURA with second order accurate reconstruction
        CALL upwind_hflux_miura( p_patch         = p_patch,            & !in
          &                      p_cc            = p_cc(:,:,:,jt),     & !in
          &                      p_mass_flx_e    = p_mass_flx_e,       & !in
          &                      p_dtime         = p_dtime,            & !in
          &                      p_int           = p_int,              & !in
          &                      btraj           = btraj,              & !in
          &                      p_igrad_c_miura = p_igrad_c_miura,    & !in
          &                      p_itype_hlimit  = p_itype_hlimit(jt), & !in
          &                      p_out_e         = p_upflux(:,:,:,jt), & !inout
          &                      opt_rhodz_now   = p_rhodz_now,        & !in
          &                      opt_rhodz_new   = p_rhodz_new,        & !in
          &                      opt_lconsv      = llsq_lin_consv,     & !in
          &                      opt_rlend       = i_rlend,            & !in
          &                      opt_slev        = p_iadv_slev(jt)     ) !in


      CASE( MIURA3 )  ! ihadv_tracer = 3

        iadv_min_slev = advection_config(jg)%miura3_h%iadv_min_slev

        ! CALL MIURA with third order accurate reconstruction
        CALL upwind_hflux_miura3( p_patch, p_cc(:,:,:,jt), p_mass_flx_e, &! in
          &                p_vn, z_real_vt, p_dtime, p_int,              &! in
          &                lcompute%miura3_h(jt), lcleanup%miura3_h(jt), &! in
          &                p_itype_hlimit(jt),                           &! in
          &                p_upflux(:,:,:,jt),                           &! inout
          &                opt_rhodz_now  = p_rhodz_now,                 &! in
          &                opt_rhodz_new  = p_rhodz_new,                 &! in
          &                opt_rlend      = i_rlend,                     &! in
          &                opt_slev       = p_iadv_slev(jt),             &! in
          &                opt_ti_slev    = iadv_min_slev                )! in

      CASE( FFSL )  ! ihadv_tracer = 4

        iadv_min_slev = advection_config(jg)%ffsl_h%iadv_min_slev

#ifdef _OPENACC
! In GPU mode, copy data to HOST and perform upwind_hflux_ffsl there, then update device
! NOTE: this is only for testing; use upwind_hflux_miura/miura3 for performance
        WRITE(message_text,'(a)') 'GPU mode: performing upwind_hflux_ffsl on host; for performance use upwind_hflux_miura'
        CALL message(TRIM(routine),message_text)
!$ACC UPDATE HOST( p_cc(:,:,:,jt), p_mass_flx_e, p_vn, p_rhodz_now, p_rhodz_new, z_real_vt ), &
!$ACC        IF( i_am_accel_node .AND. acc_on )
        save_i_am_accel_node = i_am_accel_node
        i_am_accel_node = .FALSE.     ! deactivate GPUs throughout upwind_hflux_ffsl
#endif

        ! CALL Flux form semi Lagrangian scheme (extension of MIURA3-scheme)
        ! with second or third order accurate reconstruction
        CALL upwind_hflux_ffsl( p_patch, p_cc(:,:,:,jt),                  &! in
          &                 p_rhodz_now, p_rhodz_new, p_mass_flx_e,       &! in
          &                 p_vn, z_real_vt, p_dtime, p_int,              &! in
          &                 lcompute%ffsl_h(jt), lcleanup%ffsl_h(jt),     &! in
          &                 p_itype_hlimit(jt),                           &! in
          &                 p_upflux(:,:,:,jt),                           &! inout
          &                 opt_lconsv  = llsq_high_consv,                &! in
          &                 opt_rlend   = i_rlend,                        &! in
          &                 opt_slev    = p_iadv_slev(jt),                &! in
          &                 opt_ti_slev = iadv_min_slev                   )! in

#ifdef _OPENACC
        i_am_accel_node =  save_i_am_accel_node    ! reactivate GPUs if appropriate
!$ACC UPDATE DEVICE( p_upflux(:,:,:,jt) ), IF( i_am_accel_node .AND. acc_on )
#endif

      CASE( FFSL_HYB )  ! ihadv_tracer = 5

        iadv_min_slev = advection_config(jg)%ffsl_hyb_h%iadv_min_slev

#ifdef _OPENACC
! In GPU mode, copy data to HOST and perform hflux_ffsl_hybrid there, then update device
! NOTE: this is only for testing; use upwind_hflux_miura/miura3 for performance
        WRITE(message_text,'(a)') 'GPU mode: performing hflux_ffsl_hybrid on host; for performance use upwind_hflux_miura'
        CALL message(TRIM(routine),message_text)
!$ACC UPDATE HOST( p_cc(:,:,:,jt), p_mass_flx_e, p_vn, p_rhodz_now, p_rhodz_new, z_real_vt ), &
!$ACC        IF( i_am_accel_node .AND. acc_on )
        save_i_am_accel_node = i_am_accel_node
        i_am_accel_node = .FALSE.     ! deactivate GPUs throughout hflux_ffsl_hybrid
#endif

        ! CALL hybrid FFSL/miura3 scheme (i.e. if a prescribed CFL threshold
        ! is exceeded, the scheme switches from MIURA3 to FFSL
        CALL hflux_ffsl_hybrid( p_patch, p_cc(:,:,:,jt),                  &! in
          &                 p_rhodz_now, p_rhodz_new, p_mass_flx_e,       &! in
          &                 p_vn, z_real_vt, p_dtime, p_int,              &! in
          &                 lcompute%ffsl_hyb_h(jt),                      &! in
          &                 lcleanup%ffsl_hyb_h(jt),                      &! in
          &                 p_itype_hlimit(jt),                           &! in
          &                 p_upflux(:,:,:,jt),                           &! inout
          &                 opt_lconsv  = llsq_high_consv,                &! in
          &                 opt_rlend   = i_rlend,                        &! in
          &                 opt_slev    = p_iadv_slev(jt),                &! in
          &                 opt_ti_slev = iadv_min_slev                   )! in

#ifdef _OPENACC
        i_am_accel_node =  save_i_am_accel_node    ! reactivate GPUs if appropriate
!$ACC UPDATE DEVICE( p_upflux(:,:,:,jt) ), IF( i_am_accel_node .AND. acc_on )
#endif

      CASE ( MCYCL )   ! ihadv_tracer = 20

        ! CALL MIURA with second order accurate reconstruction and subcycling
        CALL upwind_hflux_miura_cycl( p_patch         = p_patch,            & !in
          &                           p_cc            = p_cc(:,:,:,jt),     & !in
          &                           p_rhodz_now     = p_rhodz_now,        & !in
          &                           p_mass_flx_e    = p_mass_flx_e,       & !in
          &                           p_dtime         = p_dtime,            & !in
          &                           p_ncycl         = nsubsteps,          & !in
          &                           p_int           = p_int,              & !in
          &                           btraj           = btraj_cycl,         & !in
          &                           p_igrad_c_miura = p_igrad_c_miura,    & !in
          &                           p_itype_hlimit  = p_itype_hlimit(jt), & !in
          &                           p_out_e         = p_upflux(:,:,:,jt), & !inout
          &                           elev            = p_patch%nlev,       & !in
          &                           opt_lconsv      = llsq_lin_consv,     & !in
          &                           opt_rlend       = i_rlend             ) !in


      CASE( MIURA_MCYCL )   ! ihadv_tracer = 22

        qvsubstep_elev = advection_config(jg)%iadv_qvsubstep_elev

        ! CALL standard MIURA for lower atmosphere and the subcycling version of
        ! MIURA for upper atmosphere
        CALL upwind_hflux_miura( p_patch         = p_patch,            & !in
          &                      p_cc            = p_cc(:,:,:,jt),     & !in
          &                      p_mass_flx_e    = p_mass_flx_e,       & !in
          &                      p_dtime         = p_dtime,            & !in
          &                      p_int           = p_int,              & !in
          &                      btraj           = btraj,              & !in
          &                      p_igrad_c_miura = p_igrad_c_miura,    & !in
          &                      p_itype_hlimit  = p_itype_hlimit(jt), & !in
          &                      p_out_e         = p_upflux(:,:,:,jt), & !inout
          &                      opt_rhodz_now   = p_rhodz_now,        & !in
          &                      opt_rhodz_new   = p_rhodz_new,        & !in
          &                      opt_lconsv      = llsq_lin_consv,     & !in
          &                      opt_rlend       = i_rlend,            & !in
          &                      opt_slev        = qvsubstep_elev+1,   & !in
          &                      opt_elev        = p_patch%nlev        ) !in


        IF (qvsubstep_elev > 0) THEN

        ! Note that lcompute/lcleanup%mcycl_h is generally used for miura
        ! with substepping. This prevents us from computing the backward
        ! trajectories multiple times when combining the substepping scheme with
        ! different other schemes.

        CALL upwind_hflux_miura_cycl( p_patch         = p_patch,            & !in
          &                           p_cc            = p_cc(:,:,:,jt),     & !in
          &                           p_rhodz_now     = p_rhodz_now,        & !in
          &                           p_mass_flx_e    = p_mass_flx_e,       & !in
          &                           p_dtime         = p_dtime,            & !in
          &                           p_ncycl         = nsubsteps,          & !in
          &                           p_int           = p_int,              & !in
          &                           btraj           = btraj_cycl,         & !in
          &                           p_igrad_c_miura = p_igrad_c_miura,    & !in
          &                           p_itype_hlimit  = p_itype_hlimit(jt), & !in
          &                           p_out_e         = p_upflux(:,:,:,jt), & !inout
          &                           elev            = qvsubstep_elev,     & !in
          &                           opt_lconsv      = llsq_lin_consv,     & !in
          &                           opt_rlend       = i_rlend,            & !in
          &                           opt_slev        = p_iadv_slev(jt)     ) !in
        ENDIF


      CASE( MIURA3_MCYCL )   ! ihadv_tracer = 32

        qvsubstep_elev = advection_config(jg)%iadv_qvsubstep_elev

        ! CALL standard MIURA3 for lower atmosphere and the subcycling version of
        ! MIURA for upper atmosphere
        CALL upwind_hflux_miura3( p_patch, p_cc(:,:,:,jt), p_mass_flx_e, &! in
          &              p_vn, z_real_vt, p_dtime, p_int,                &! in
          &              lcompute%miura3_h(jt), lcleanup%miura3_h(jt),   &! in
          &              p_itype_hlimit(jt),                             &! in
          &              p_upflux(:,:,:,jt),                             &! inout
          &              opt_rhodz_now  = p_rhodz_now,                   &! in
          &              opt_rhodz_new  = p_rhodz_new,                   &! in
          &              opt_rlend      = i_rlend,                       &! in
          &              opt_slev       = qvsubstep_elev+1,              &! in
          &              opt_elev       = p_patch%nlev,                  &! in
          &              opt_ti_slev    = qvsubstep_elev+1,              &! in
          &              opt_ti_elev    = p_patch%nlev                   )! in

        IF (qvsubstep_elev > 0) THEN

        ! Note that lcompute/lcleanup%mcycl_h is generally used for miura
        ! with substepping. This prevents us from computing the backward
        ! trajectories multiple times when combining the substepping scheme with
        ! different other schemes.
        CALL upwind_hflux_miura_cycl( p_patch         = p_patch,            & !in
          &                           p_cc            = p_cc(:,:,:,jt),     & !in
          &                           p_rhodz_now     = p_rhodz_now,        & !in
          &                           p_mass_flx_e    = p_mass_flx_e,       & !in
          &                           p_dtime         = p_dtime,            & !in
          &                           p_ncycl         = nsubsteps,          & !in
          &                           p_int           = p_int,              & !in
          &                           btraj           = btraj_cycl,         & !in
          &                           p_igrad_c_miura = p_igrad_c_miura,    & !in
          &                           p_itype_hlimit  = p_itype_hlimit(jt), & !in
          &                           p_out_e         = p_upflux(:,:,:,jt), & !inout
          &                           elev            = qvsubstep_elev,     & !in
          &                           opt_lconsv      = llsq_lin_consv,     & !in
          &                           opt_rlend       = i_rlend,            & !in
          &                           opt_slev        = p_iadv_slev(jt)     ) !in
        ENDIF


      CASE (FFSL_MCYCL)   ! ihadv_tracer = 42

        qvsubstep_elev = advection_config(jg)%iadv_qvsubstep_elev

#ifdef _OPENACC
! In GPU mode, copy data to HOST and perform upwind_hflux_ffsl there, then update device
! NOTE: this is only for testing; use upwind_hflux_miura/miura3 for performance
        WRITE(message_text,'(a)') 'GPU mode: performing upwind_hflux_ffsl on host; for performance use upwind_hflux_miura'
        CALL message(TRIM(routine),message_text)
!$ACC UPDATE HOST( p_cc(:,:,:,jt), p_mass_flx_e, p_vn, p_rhodz_now, p_rhodz_new, z_real_vt ), &
!$ACC        IF( i_am_accel_node .AND. acc_on )
        save_i_am_accel_node = i_am_accel_node
        i_am_accel_node = .FALSE.     ! deactivate GPUs throughout hflux_ffsl_hybrid
#endif

        ! CALL standard FFSL for lower atmosphere and the subcycling version of
        ! MIURA for upper atmosphere
        CALL upwind_hflux_ffsl( p_patch, p_cc(:,:,:,jt),                  &! in
          &                 p_rhodz_now, p_rhodz_new, p_mass_flx_e,       &! in
          &                 p_vn, z_real_vt, p_dtime, p_int,              &! in
          &                 lcompute%ffsl_h(jt), lcleanup%ffsl_h(jt),     &! in
          &                 p_itype_hlimit(jt),                           &! in
          &                 p_upflux(:,:,:,jt),                           &! inout
          &                 opt_lconsv  = llsq_high_consv,                &! in
          &                 opt_rlend   = i_rlend,                        &! in
          &                 opt_slev    = qvsubstep_elev+1,               &! in
          &                 opt_elev    = p_patch%nlev,                   &! in
          &                 opt_ti_slev = qvsubstep_elev+1,               &! in
          &                 opt_ti_elev = p_patch%nlev                    )! in

#ifdef _OPENACC
        i_am_accel_node =  save_i_am_accel_node    ! reactivate GPUs if appropriate
!$ACC UPDATE DEVICE( p_upflux(:,:,:,jt) ), IF( i_am_accel_node .AND. acc_on )
#endif

        IF (qvsubstep_elev > 0) THEN

        ! Note that lcompute/lcleanup%mcycl_h is generally used for miura
        ! with substepping. This prevents us from computing the backward
        ! trajectories multiple times when combining the substepping scheme with
        ! different other schemes.
        CALL upwind_hflux_miura_cycl( p_patch         = p_patch,            & !in
          &                           p_cc            = p_cc(:,:,:,jt),     & !in
          &                           p_rhodz_now     = p_rhodz_now,        & !in
          &                           p_mass_flx_e    = p_mass_flx_e,       & !in
          &                           p_dtime         = p_dtime,            & !in
          &                           p_ncycl         = nsubsteps,          & !in
          &                           p_int           = p_int,              & !in
          &                           btraj           = btraj_cycl,         & !in
          &                           p_igrad_c_miura = p_igrad_c_miura,    & !in
          &                           p_itype_hlimit  = p_itype_hlimit(jt), & !in
          &                           p_out_e         = p_upflux(:,:,:,jt), & !inout
          &                           elev            = qvsubstep_elev,     & !in
          &                           opt_lconsv      = llsq_lin_consv,     & !in
          &                           opt_rlend       = i_rlend,            & !in
          &                           opt_slev        = p_iadv_slev(jt)     ) !in
        ENDIF


      CASE (FFSL_HYB_MCYCL)   ! ihadv_tracer = 52

        qvsubstep_elev = advection_config(jg)%iadv_qvsubstep_elev

#ifdef _OPENACC
! In GPU mode, copy data to HOST and perform hflux_ffsl_hybrid there, then update device
! NOTE: this is only for testing; use upwind_hflux_miura/miura3 for performance
        WRITE(message_text,'(a)') 'GPU mode: performing hflux_ffsl_hybrid on host; for performance use upwind_hflux_miura'
        CALL message(TRIM(routine),message_text)
!$ACC UPDATE HOST( p_cc(:,:,:,jt), p_mass_flx_e, p_vn, p_rhodz_now, p_rhodz_new, z_real_vt ), &
!$ACC        IF( i_am_accel_node .AND. acc_on )
        save_i_am_accel_node = i_am_accel_node
        i_am_accel_node = .FALSE.     ! deactivate GPUs throughout hflux_ffsl_hybrid
#endif

        ! CALL hybrid FFSL/MIURA3 for lower atmosphere and the subcycling
        ! version of MIURA for upper atmosphere
        CALL hflux_ffsl_hybrid( p_patch, p_cc(:,:,:,jt),                  &! in
          &                 p_rhodz_now, p_rhodz_new, p_mass_flx_e,       &! in
          &                 p_vn, z_real_vt, p_dtime, p_int,              &! in
          &                 lcompute%ffsl_hyb_h(jt),                      &! in
          &                 lcleanup%ffsl_hyb_h(jt),                      &! in
          &                 p_itype_hlimit(jt),                           &! in
          &                 p_upflux(:,:,:,jt),                           &! inout
          &                 opt_lconsv  = llsq_high_consv,                &! in
          &                 opt_rlend   = i_rlend,                        &! in
          &                 opt_slev    = qvsubstep_elev+1,               &! in
          &                 opt_elev    = p_patch%nlev,                   &! in
          &                 opt_ti_slev = qvsubstep_elev+1,               &! in
          &                 opt_ti_elev = p_patch%nlev                    )! in

#ifdef _OPENACC
        i_am_accel_node =  save_i_am_accel_node    ! reactivate GPUs if appropriate
!$ACC UPDATE DEVICE( p_upflux(:,:,:,jt) ), IF( i_am_accel_node .AND. acc_on )
#endif

        IF (qvsubstep_elev > 0) THEN

        ! Note that lcompute/lcleanup%mcycl_h is generally used for miura
        ! with substepping. This prevents us from computing the backward
        ! trajectories multiple times when combining the substepping scheme with
        ! different other schemes.
        CALL upwind_hflux_miura_cycl( p_patch         = p_patch,            & !in
          &                           p_cc            = p_cc(:,:,:,jt),     & !in
          &                           p_rhodz_now     = p_rhodz_now,        & !in
          &                           p_mass_flx_e    = p_mass_flx_e,       & !in
          &                           p_dtime         = p_dtime,            & !in
          &                           p_ncycl         = nsubsteps,          & !in
          &                           p_int           = p_int,              & !in
          &                           btraj           = btraj_cycl,         & !in
          &                           p_igrad_c_miura = p_igrad_c_miura,    & !in
          &                           p_itype_hlimit  = p_itype_hlimit(jt), & !in
          &                           p_out_e         = p_upflux(:,:,:,jt), & !inout
          &                           elev            = qvsubstep_elev,     & !in
          &                           opt_lconsv      = llsq_lin_consv,     & !in
          &                           opt_rlend       = i_rlend,            & !in
          &                           opt_slev        = p_iadv_slev(jt)     ) !in
        ENDIF

      END SELECT

    END DO  ! Tracer loop

!$ACC UPDATE HOST( p_upflux ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

    CALL btraj%destruct()
    CALL btraj_cycl%destruct()

    IF (timers_level > 2) CALL timer_stop(timer_adv_horz)

  END SUBROUTINE hor_upwind_flux




  !-----------------------------------------------------------------------
  !>
  !! The first order Godunov scheme
  !!
  !! Calculation of time averaged horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !! Developed by L.Bonaventura  (2004).
  !! Adapted to new grid structure by L. Bonaventura, MPI-M, August 2005.
  !! Modification by Daniel Reinert, DWD (2010-02-09)
  !! - transferred to separate subroutine
  !!
  SUBROUTINE upwind_hflux_up( p_patch, p_cc, p_mass_flx_e, p_upflux, &
    &                       opt_rlstart, opt_rlend, opt_slev, opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &    !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::     &  !< advected cell centered variable
      &  p_cc(:,:,:)               !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &  !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)       !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(INOUT) ::  &  !< variable in which the upwind flux is stored
      &  p_upflux(:,:,:)           !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
      &  opt_elev

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iilc, iibc
    INTEGER  :: slev, elev           !< vertical start and end level
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block

    !-----------------------------------------------------------------------

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 2
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a picewise constant approx. of the cell centered values
    ! is used.
    !

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    ! line and block indices of two neighboring cells
    iilc => p_patch%edges%cell_idx
    iibc => p_patch%edges%cell_blk

    ! loop through all patch edges (and blocks)

!$ACC DATA  PCOPYIN( p_cc, p_mass_flx_e ), PCOPYOUT( p_upflux ), &
!$ACC       PRESENT( iilc, iibc ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc, p_mass_flx_e ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,   &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=6
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif
          !
          ! compute the first order upwind flux; notice
          ! that only the p_cc*p_mass_flx_e value at cell edge is computed
          ! multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          !
          IF ( p_mass_flx_e(je,jk,jb) .GE. 0.0_wp ) THEN
            p_upflux(je,jk,jb) = p_mass_flx_e(je,jk,jb) * p_cc(iilc(je,jb,2),jk,iibc(je,jb,2))
          ELSE
            p_upflux(je,jk,jb) = p_mass_flx_e(je,jk,jb) * p_cc(iilc(je,jb,1),jk,iibc(je,jb,1))
          ENDIF

        END DO  ! end loop over edges

      END DO  ! end loop over levels
!$ACC END PARALLEL

    END DO  ! end loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_upflux ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE upwind_hflux_up



  !-------------------------------------------------------------------------
  !>
  !! The second order MIURA scheme
  !!
  !! Calculation of time averaged horizontal tracer fluxes at triangle edges
  !! using a Flux form semi-lagrangian (FFSL)-method based on the second order
  !! MIURA-scheme.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2009-11-09)
  !! Modification by Daniel Reinert, DWD (2010-02-10)
  !! - transferred to separate subroutine
  !! Modification by Daniel Reinert, DWD (2010-03-16)
  !! - implemented new computation of backward trajectories on a local tangential
  !!  plane for each edge-midpoint (much faster than the old version)
  !! - moved calculation of backward trajectories and barycenter into subroutine
  !!   back_traj_o1 in module mo_advection_utils. Added second order accurate
  !!   computation of backward trajectories (subroutine back_traj_o2)
  !! Modification by Daniel Reinert, DWD (2016-11-24)
  !! - computation of backward trajectories moved one level up which simplifys 
  !!   flow control. It is now simpler to ensure that backward trajectories are only 
  !!   computed once per time step and that the backward trajectory information is 
  !!   available at all points and levels.
  !!
  !! @par LITERATURE
  !! - Miura, H. (2007), Mon. Weather Rev., 135, 4038-4044
  !! - Skamarock, W.C. (2010), Conservative Transport schemes for Spherical Geodesic
  !!   Grids: High-order Reconstructions for Forward-in-Time Schemes, Mon. Wea. Rev,
  !!   139, 4497-4508
  !!
  SUBROUTINE upwind_hflux_miura( p_patch, p_cc, p_mass_flx_e, p_dtime,      &
    &                      p_int, btraj, p_igrad_c_miura, p_itype_hlimit,   &
    &                      p_out_e, opt_rhodz_now, opt_rhodz_new,           &
    &                      opt_lconsv, opt_rlstart, opt_rlend,              &
    &                      opt_lout_edge, opt_slev, opt_elev )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_hflux: upwind_hflux_miura'

    TYPE(t_patch), TARGET, INTENT(INOUT) ::  &   !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  p_int

    TYPE(t_back_traj), INTENT(IN) :: &      !< information on backward trajectories
      &  btraj

    REAL(wp), TARGET, INTENT(IN) ::     &   !< cell centered variable to be advected
      &  p_cc(:,:,:)                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    INTEGER, INTENT(IN) :: p_igrad_c_miura   !< parameter to select the gradient
                                             !< reconstruction method at cell center

    INTEGER, INTENT(IN) :: p_itype_hlimit    !< parameter to select the limiter
                                             !< for horizontal transport

    REAL(wp), INTENT(INOUT) ::  &   !< output field, containing the upwind flux or the
      &  p_out_e(:,:,:)             !< reconstructed edge value; dim: (nproma,nlev,nblks_e)

    REAL(wp), OPTIONAL, INTENT(IN) ::    &  !< density times cell thickness at cell center step (n)
      &  opt_rhodz_now(:,:,:)               !< dim: (nproma,nlev,nblks_c)
                                            !< needed by flux limiter routines, only

    REAL(wp), OPTIONAL, INTENT(IN) ::    &  !< density times cell thickness at cell center step (n+1)
      &  opt_rhodz_new(:,:,:)               !< dim: (nproma,nlev,nblks_c)
                                            !< needed by flux limiter routines, only

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: if true, conservative reconstruction
      &  opt_lconsv                    !< is used

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: output edge value (.TRUE.),
      &  opt_lout_edge                 !< or the flux across the edge (.FALSE./not specified)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    LOGICAL  :: l_out_edgeval          !< corresponding local variable; default .FALSE.
                                       !< i.e. output flux across the edge

    REAL(vp), TARGET ::    &                   !< reconstructed gradient vector at
      &  z_grad(2,nproma,p_patch%nlev,p_patch%nblks_c)
                                               !< cell center (geographical coordinates)

    REAL(wp), TARGET ::    &                        !< coefficient of the lsq reconstruction
      &  z_lsq_coeff(3,nproma,p_patch%nlev,p_patch%nblks_c)
                                                    !< at cell center (geogr. coordinates)
                                                    !< includes coeff0 and gradients in
                                                    !< zonal and meridional direction

    INTEGER  :: pid
    INTEGER  :: nlev               !< number of full levels
    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    INTEGER  :: ilc0, ibc0         !< line and block index for local cell center
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i, i_rlstart, i_rlend, i_nchdom, i_rlend_c
    LOGICAL  :: l_consv            !< true if conservative lsq reconstruction is used
    LOGICAL  :: use_zlsq           !< true if z_lsq_coeff is used to store the gradients
    TYPE(t_lsq), POINTER :: lsq_lin  !< pointer to p_int_state%lsq_lin
    lsq_lin => p_int%lsq_lin

   !-------------------------------------------------------------------------

!$ACC DATA  PCOPYIN( p_cc, p_mass_flx_e, btraj), PCOPY( p_out_e ), CREATE( z_grad, z_lsq_coeff ), &
!$ACC       PRESENT( btraj%cell_idx, btraj%cell_blk, btraj%distv_bary ), IF( i_am_accel_node .AND. acc_on)
!$ACC DATA  PCOPYIN( opt_rhodz_now ), IF( PRESENT(opt_rhodz_now) .AND. i_am_accel_node .AND. acc_on )
!$ACC DATA  PCOPYIN( opt_rhodz_new ), IF( PRESENT(opt_rhodz_new) .AND. i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc, p_mass_flx_e, btraj, p_out_e ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_grad,z_lsq_coeff
#endif


    ! number of vertical levels
    nlev = p_patch%nlev

    ! get patch ID
    pid = p_patch%id

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = nlev
    END IF

    IF ( PRESENT(opt_lconsv) ) THEN
     l_consv = opt_lconsv
    ELSE
     l_consv = .FALSE. ! non-conservative reconstruction
    ENDIF

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF


    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 5
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF


    i_rlend_c = min_rlcell_int - 1


    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    IF (p_test_run) THEN
!$ACC KERNELS DEFAULT(PRESENT) IF (i_am_accel_node .AND. acc_on)
#ifdef __INTEL_COMPILER
!$OMP PARALLEL DO SCHEDULE(STATIC)
      DO i = 1,SIZE(z_grad,4)
        z_grad(:,:,:,i) = 0._wp
      ENDDO
#else
      z_grad(:,:,:,:) = 0._wp
#endif
!$ACC END KERNELS
    ENDIF

    !
    ! advection is done with an upwind scheme and a piecewise linear
    ! approx. of the cell centered values.
    ! This approx. is evaluated at the barycenter of the area which is
    ! advected through the edge under consideration. This area is
    ! approximated as a rhomboid. The area approximation is based on the
    ! (reconstructed) full 2D velocity field at edge midpoints (at
    ! time t+\Delta t/2) and \Delta t.
    !
    ! 2 options:  without limiter
    !             with    flux limiter following Zalesak (1979)
    !


    use_zlsq = .FALSE. ! default: use z_grad to store the gradient
    !
    ! 2. reconstruction of (unlimited) cell based gradient (lat,lon)
    !
    IF (p_igrad_c_miura == 1) THEN
      ! least squares method
      IF (advection_config(pid)%llsq_svd .AND. l_consv) THEN
        CALL recon_lsq_cell_l_consv_svd( p_cc, p_patch, lsq_lin, z_lsq_coeff,         &
        &                              opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                              opt_lconsv=l_consv)
        use_zlsq = .TRUE.
      ELSE IF (advection_config(pid)%llsq_svd) THEN
        CALL recon_lsq_cell_l_svd( p_cc, p_patch, lsq_lin, z_grad,           &
        &                        opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c)
      ELSE
        CALL recon_lsq_cell_l( p_cc, p_patch, lsq_lin, z_lsq_coeff,         &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_lconsv=l_consv)
        use_zlsq = .TRUE.
      ENDIF

    ELSE IF (p_igrad_c_miura == 2) THEN
      ! Green-Gauss method
      CALL grad_green_gauss_cell( p_cc, p_patch, p_int, z_grad, opt_slev=slev, &
        &                         opt_elev=elev, opt_rlend=i_rlend_c )


    ELSE IF (p_igrad_c_miura == 3) THEN
      ! gradient based on three-node triangular element
      CALL grad_fe_cell( p_cc, p_patch, p_int, z_grad, opt_slev=slev, &
        &                opt_elev=elev, opt_rlend=i_rlend_c )

    ENDIF




    !
    ! 3. Calculate reconstructed tracer value at each barycenter
    !    \Phi_{bary}=\Phi_{circum} + DOT_PRODUCT(\Nabla\Psi,r).
    !    Then calculate the flux v_n*\Delta p*\Phi_{bary}
    !    The fact that the rhomboidal area inevitably overlaps with neighboring
    !    triangles is neglected (local linear approximation instead of piecewise
    !    linear approximation). Only the reconstruction for the local cell
    !    is taken into account.

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! Before starting, preset halo edges that are not processed with zero's in order
    ! to avoid access of uninitialized array elements in subsequent routines
    ! Necessary when called within dycore

    IF ( l_out_edgeval ) THEN
      i_startblk = p_patch%edges%start_blk(i_rlend-1,i_nchdom)
      i_endblk   = p_patch%edges%end_blk(min_rledge_int-3,i_nchdom)

      CALL init(p_out_e(:,:,i_startblk:i_endblk))
!$OMP BARRIER
    ENDIF

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    ! initialize also nest boundary points with zero
    IF ( l_out_edgeval .AND. (p_patch%id > 1 .OR. l_limited_area)) THEN
      CALL init(p_out_e(:,:,1:i_startblk))
!$OMP BARRIER
    ENDIF

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,ilc0,ibc0), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      IF ( l_out_edgeval ) THEN   ! Calculate 'edge value' of advected quantity

!CDIR UNROLL=5
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR PRIVATE(ilc0,ibc0) COLLAPSE(2)
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

            ! Calculate reconstructed tracer value at barycenter of rhomboidal
            ! area which is swept across the corresponding edge.
            ilc0 = btraj%cell_idx(je,jk,jb)
            ibc0 = btraj%cell_blk(je,jk,jb)

            ! Calculate 'edge value' of advected quantity (cc_bary)
            p_out_e(je,jk,jb) = p_cc(ilc0,jk,ibc0)                           &
              &    + btraj%distv_bary(je,jk,jb,1) * z_grad(1,ilc0,jk,ibc0)   &
              &    + btraj%distv_bary(je,jk,jb,2) * z_grad(2,ilc0,jk,ibc0)

          ENDDO ! loop over edges
        ENDDO   ! loop over vertical levels
!$ACC END PARALLEL

      ELSE IF (use_zlsq) THEN

!CDIR UNROLL=5
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR PRIVATE(ilc0,ibc0) COLLAPSE(2)
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

            ! Calculate reconstructed tracer value at barycenter of rhomboidal
            ! area which is swept across the corresponding edge.
            ilc0 = btraj%cell_idx(je,jk,jb)
            ibc0 = btraj%cell_blk(je,jk,jb)

            ! Calculate flux at cell edge (cc_bary*v_{n}* \Delta p)
            p_out_e(je,jk,jb) = ( z_lsq_coeff(1,ilc0,jk,ibc0)                      &
              &    + btraj%distv_bary(je,jk,jb,1) * z_lsq_coeff(2,ilc0,jk,ibc0)    &
              &    + btraj%distv_bary(je,jk,jb,2) * z_lsq_coeff(3,ilc0,jk,ibc0) )  &
              &    * p_mass_flx_e(je,jk,jb)

          ENDDO ! loop over edges
        ENDDO   ! loop over vertical levels
!$ACC END PARALLEL

      ELSE

!CDIR UNROLL=5
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR PRIVATE(ilc0,ibc0), COLLAPSE(2)
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

            ! Calculate reconstructed tracer value at barycenter of rhomboidal
            ! area which is swept across the corresponding edge.
            ilc0 = btraj%cell_idx(je,jk,jb)
            ibc0 = btraj%cell_blk(je,jk,jb)

            ! Calculate flux at cell edge (cc_bary*v_{n}* \Delta p)
            p_out_e(je,jk,jb) = ( p_cc(ilc0,jk,ibc0)                          &
              &    + btraj%distv_bary(je,jk,jb,1) * z_grad(1,ilc0,jk,ibc0)    &
              &    + btraj%distv_bary(je,jk,jb,2) * z_grad(2,ilc0,jk,ibc0) )  &
              &    * p_mass_flx_e(je,jk,jb)

          ENDDO ! loop over edges
        ENDDO   ! loop over vertical levels
!$ACC END PARALLEL

      ENDIF

    ENDDO    ! loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !
    ! 4. If desired, apply a monotonic or positive definite flux limiter 
    !    to limit computed fluxes.
    !    The flux limiter is based on work by Zalesak (1979)
    IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_m) THEN
      IF (.NOT. (PRESENT(opt_rhodz_now) .OR. PRESENT(opt_rhodz_now))) THEN
        CALL finish(TRIM(routine),'Required fields opt_rhodz_now and opt_rhodz_now are not present')
      ENDIF
      CALL hflx_limiter_mo( p_patch, p_int, p_dtime, p_cc,           & !in
        &                   opt_rhodz_now, opt_rhodz_new, p_mass_flx_e,  & !in
        &                   p_out_e, slev, elev, opt_rlend=i_rlend   ) !inout,in

    ELSE IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_sm) THEN
      IF (.NOT. (PRESENT(opt_rhodz_now))) THEN
        CALL finish(TRIM(routine),'Required field opt_rhodz_now not present')
      ENDIF
      ! MPI-sync necessary
      CALL hflx_limiter_pd( p_patch, p_int, p_dtime,                 & !in
        &                   p_cc, opt_rhodz_now, p_out_e,              & !in,inout
        &                   slev, elev, opt_rlend=i_rlend            ) !in
    ENDIF

!$ACC UPDATE HOST( p_out_e ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA

  END SUBROUTINE upwind_hflux_miura




  !-------------------------------------------------------------------------
  !>
  !! The second order MIURA scheme with subcycling-option
  !!
  !! Calculation of time averaged horizontal tracer fluxes at triangle edges
  !! using a Flux form semi-lagrangian (FFSL)-method based on the second order
  !! MIURA-scheme. This particular version of the scheme includes a time
  !! subcyling-option.
  !!
  !! @par Revision History
  !! - Initial revision by Daniel Reinert, DWD (2011-09-20)
  !! Modification by Daniel Reinert, DWD (2012-01-25)
  !! - bug fix for positive definite limiter. Limiter is now called after each
  !!   substep.
  !! Modification by Daniel Reinert, DWD (2013-09-27)
  !! - increased number of subcycling steps from 2 to 3 (hard coded)
  !! Modification by Daniel Reinert, DWD (2013-10-09)
  !! - reduce vertical dimension of local arrays from nlev to elev
  !! Modification by Daniel Reinert, DWD (2016-11-24)
  !! - computation of backward trajectories moved one level up to simplify 
  !!   flow control. It is now simpler to ensure that backward trajectories are only 
  !!   computed once per time step and that the backward trajectory information is 
  !!   available at all points and levels.
  !!
  !! @par LITERATURE
  !! - Miura, H. (2007), Mon. Weather Rev., 135, 4038-4044
  !! - Skamarock, W.C. (2010), Conservative Transport schemes for Spherical Geodesic
  !!   Grids: High-order Reconstructions for Forward-in-Time Schemes, Mon. Wea. Rev,
  !!   139, 4497-4508
  !!
  SUBROUTINE upwind_hflux_miura_cycl( p_patch, p_cc, p_rhodz_now,              &
    &                   p_mass_flx_e, p_dtime,  p_ncycl, p_int, btraj,         &
    &                   p_igrad_c_miura, p_itype_hlimit, p_out_e,              &
    &                   elev, opt_lconsv, opt_rlstart, opt_rlend, opt_slev )


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_hflux: upwind_hflux_miura_cycl'

    TYPE(t_patch), TARGET, INTENT(INOUT) ::  &   !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  p_int

    TYPE(t_back_traj), INTENT(IN) :: &      !< information on backward trajectories
      &  btraj

    REAL(wp), TARGET, INTENT(IN) ::     &   !< cell centered variable to be advected
      &  p_cc(:,:,:)                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< density times cell thickness at timestep n
      &  p_rhodz_now(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)        !< Assumption: constant over p_dtime
                                    !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    INTEGER,  INTENT(IN) ::    &    !< number of sub-timesteps into which p_dtime
      &  p_ncycl                    !< is split (p_ncycl=1 : no subcycling)

    INTEGER, INTENT(IN) :: p_igrad_c_miura   !< parameter to select the gradient
                                             !< reconstruction method at cell center

    INTEGER, INTENT(IN) :: p_itype_hlimit    !< parameter to select the limiter
                                             !< for horizontal transport

    REAL(wp), INTENT(INOUT) ::  &   !< output field, containing the upwind flux or the
      &  p_out_e(:,:,:)             !< reconstructed edge value; dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN)     :: & !< vertical end level; not an optional argument since it
      &  elev                    !< is used to dimension local arrays

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: if true, conservative reconstruction
     &  opt_lconsv                     !< is used

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    REAL(vp), TARGET ::    &                   !< reconstructed gradient vector at
      &  z_grad(2,nproma,p_patch%nlev,p_patch%nblks_c)
                                               !< cell center (geographical coordinates)

    REAL(wp), TARGET ::    &                        !< coefficient of the lsq reconstruction
      &  z_lsq_coeff(3,nproma,p_patch%nlev,p_patch%nblks_c)
                                                    !< at cell center (geogr. coordinates)
                                                    !< includes coeff0 and gradients in
                                                    !< zonal and meridional direction

    REAL(wp) :: z_dtsub                 !< sub timestep p_dtime/p_ncycl
    REAL(wp) ::                     &   !< tracer flux at n + nsub/p_ncycl
      &  z_tracer_mflx(nproma,elev,p_patch%nblks_e,p_ncycl)

    REAL(vp) ::                     &   !< tracer mass flux divergence at cell center
      &  z_rhofluxdiv_c(nproma,elev,p_patch%nblks_c)

    REAL(vp) ::                     &   !< mass flux divergence at cell center
      &  z_fluxdiv_c(nproma,elev)

    REAL(wp), TARGET ::             &   !< 'tracer cell value' at interm.
      &  z_tracer(nproma,elev,p_patch%nblks_c,2) !< old and new timestep

    REAL(wp) ::                     &   !< density (i.e. \rho\Delta z) at interm.
      &  z_rho(nproma,elev,p_patch%nblks_c,2) !< old and new timestep

    INTEGER  :: pid
    INTEGER  :: slev               !< vertical start level
    INTEGER  :: jc, je, jk, jb     !< index of cell, edge, vert level, block
    INTEGER  :: ilc0, ibc0         !< line and block index for local cell center
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom, i_rlend_c
    LOGICAL  :: l_consv            !< true if conservative lsq reconstruction is used
    LOGICAL  :: use_zlsq           !< true if z_lsq_coeff is used to store the gradients
    INTEGER  :: nsub               !< counter for sub-timesteps
    INTEGER  :: nnow, nnew, nsav   !< time indices

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iidx, iblk                          !< of edges
    TYPE(t_lsq), POINTER :: lsq_lin          !< Pointer to p_int_state%lsq_lin
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_grad,z_lsq_coeff,z_tracer_mflx,z_rhofluxdiv_c
!DIR$ ATTRIBUTES ALIGN : 64 :: z_fluxdiv_c,z_tracer,z_rho
#endif
    lsq_lin => p_int%lsq_lin

   !-------------------------------------------------------------------------

    IF (p_ncycl /= 2 .AND. p_ncycl /= 3) &
    CALL finish(TRIM(routine),'current implementation of upwind_hflux_miura_cycl '//&
      &                       'requires 2 or 3 subcycling steps (p_ncycl=2/3)')

    ! get patch ID
    pid = p_patch%id

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_lconsv) ) THEN
     l_consv = opt_lconsv
    ELSE
     l_consv = .FALSE. ! non-conservative reconstruction
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 5
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

!!$    IF (p_igrad_c_miura == 3) THEN
!!$      i_rlend_c = min_rlcell_int
!!$    ELSE
      i_rlend_c = min_rlcell_int - 1
!!$    ENDIF


    ! line and block indices of edges as seen from cells
    iidx => p_patch%cells%edge_idx
    iblk => p_patch%cells%edge_blk

    ! get local sub-timestep
    z_dtsub = p_dtime/REAL(p_ncycl,wp)


    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

!$ACC DATA  PCOPYIN( p_cc, p_rhodz_now, p_mass_flx_e ), PCOPYOUT( p_out_e ), &
!$ACC       CREATE( z_grad, z_lsq_coeff, z_tracer_mflx, z_rhofluxdiv_c, z_fluxdiv_c, z_tracer, z_rho ), &
!$ACC       PRESENT( iidx, iblk, btraj%cell_idx, btraj%cell_blk, btraj%distv_bary ), &
!$ACC       IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc, p_mass_flx_e ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

    IF (p_test_run) THEN
      z_grad(:,:,:,:)   = 0._wp
      z_tracer(:,:,:,:) = 0._wp
    ENDIF

    ! initialize 'now slice' of z_tracer and z_rho
    nnow = 1
    nnew = 2
!$OMP PARALLEL
    CALL copy(p_cc (:,slev:elev,:), z_tracer(:,slev:elev,:,nnow))
    CALL copy(p_rhodz_now(:,slev:elev,:), z_rho   (:,slev:elev,:,nnow))
!$OMP END PARALLEL


    !
    ! advection is done with an upwind scheme and a piecewise linear
    ! approx. of the cell centered values.
    ! This approx. is evaluated at the barycenter of the area which is
    ! advected across the edge under consideration. This area is
    ! approximated as a rhomboid. The area approximation is based on the
    ! (reconstructed) full 2D velocity field at edge midpoints (at
    ! time t+\Delta t/2) and \Delta t.
    !
    ! 2 options:  without limiter
    !             with    flux limiter following Zalesak (1979)
    !

    !
    ! Loop over sub-timesteps (subcycling)
    !
    DO nsub=1, p_ncycl

      use_zlsq = .FALSE.
      !
      ! 2. reconstruction of (unlimited) cell based gradient (lat,lon)
      !
      IF (p_igrad_c_miura == 1) THEN
        ! least squares method
        IF (advection_config(pid)%llsq_svd .AND. l_consv) THEN
          CALL recon_lsq_cell_l_consv_svd( z_tracer(:,:,:,nnow), p_patch, lsq_lin, &
          &                              z_lsq_coeff, opt_slev=slev, opt_elev=elev,        &
          &                              opt_rlend=i_rlend_c, opt_lconsv=l_consv)
          use_zlsq = .TRUE.
        ELSE IF (advection_config(pid)%llsq_svd) THEN
          CALL recon_lsq_cell_l_svd( z_tracer(:,:,:,nnow), p_patch, lsq_lin,       &
          &                    z_grad, opt_slev=slev, opt_elev=elev,             &
          &                    opt_rlend=i_rlend_c)
        ELSE
          CALL recon_lsq_cell_l( z_tracer(:,:,:,nnow), p_patch, lsq_lin,           &
          &                    z_lsq_coeff, opt_slev=slev, opt_elev=elev,        &
          &                    opt_rlend=i_rlend_c, opt_lconsv=l_consv)
          use_zlsq = .TRUE.
        ENDIF

      ELSE IF (p_igrad_c_miura == 2) THEN
        ! Green-Gauss method
        CALL grad_green_gauss_cell( z_tracer(:,:,:,nnow), p_patch, p_int, &
          &                         z_grad, opt_slev=slev, opt_elev=elev, &
          &                         opt_rlend=i_rlend_c )


      ELSE IF (p_igrad_c_miura == 3) THEN
        ! gradient based on three-node triangular element
        CALL grad_fe_cell( z_tracer(:,:,:,nnow), p_patch, p_int, &
          &                z_grad, opt_slev=slev, opt_elev=elev, &
          &                opt_rlend=i_rlend_c )


      ENDIF


      !
      ! 3. Calculate reconstructed tracer value at each barycenter
      !    \Phi_{bary}=\Phi_{circum} + DOT_PRODUCT(\Nabla\Psi,r).
      !    and compute intermediate update of q.
      !

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)


      i_startblk = p_patch%edges%start_blk(i_rlstart,1)
      i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)


    IF ( p_patch%id > 1 .OR. l_limited_area) THEN
      CALL init(z_tracer_mflx(:,:,1:i_startblk,nsub))
!$OMP BARRIER
    ENDIF

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,ilc0,ibc0) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)


        !
        ! 3.1 Compute reconstructed tracer value at barycenter of rhomboidal
        !     area which is swept across the corresponding edge.
        ! 3.2 Compute intermediate tracer mass flux
        !
        IF (use_zlsq) THEN
!$ACC PARALLEL DEFAULT(PRESENT) IF ( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR PRIVATE( ilc0, ibc0 ) COLLAPSE(2)
!CDIR UNROLL=5
          DO jk = slev, elev

            DO je = i_startidx, i_endidx

              ilc0 = btraj%cell_idx(je,jk,jb)
              ibc0 = btraj%cell_blk(je,jk,jb)

              ! compute intermediate flux at cell edge (cc_bary*v_{n}* \Delta p)
              z_tracer_mflx(je,jk,jb,nsub) = ( z_lsq_coeff(1,ilc0,jk,ibc0)        &
                &      + btraj%distv_bary(je,jk,jb,1) * z_lsq_coeff(2,ilc0,jk,ibc0)   &
                &      + btraj%distv_bary(je,jk,jb,2) * z_lsq_coeff(3,ilc0,jk,ibc0) ) &
                &      * p_mass_flx_e(je,jk,jb)

            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels
!$ACC END PARALLEL
        ELSE
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR PRIVATE( ilc0, ibc0 ) COLLAPSE(2)
!CDIR UNROLL=5
          DO jk = slev, elev

            DO je = i_startidx, i_endidx

              ilc0 = btraj%cell_idx(je,jk,jb)
              ibc0 = btraj%cell_blk(je,jk,jb)

              ! compute intermediate flux at cell edge (cc_bary*v_{n}* \Delta p)
              z_tracer_mflx(je,jk,jb,nsub) = ( z_tracer(ilc0,jk,ibc0,nnow)   &
                &      + btraj%distv_bary(je,jk,jb,1) * z_grad(1,ilc0,jk,ibc0)   &
                &      + btraj%distv_bary(je,jk,jb,2) * z_grad(2,ilc0,jk,ibc0) ) &
                &      * p_mass_flx_e(je,jk,jb)

            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels
!$ACC END PARALLEL
        ENDIF

      ENDDO    ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL


      ! 4. limit intermediate tracer fluxes to achieve positive definiteness
      !    The flux limiter is based on work by Zalesak (1979)
      !
      IF ( p_itype_hlimit == ifluxl_sm .OR. p_itype_hlimit == ifluxl_m ) THEN
        !
        CALL hflx_limiter_pd( p_patch, p_int, z_dtsub          , & !in
          &                   z_tracer(:,:,:,nnow)             , & !in
          &                   z_rho(:,:,:,nnow)                , & !in
          &                   z_tracer_mflx(:,:,:,nsub)        , & !inout
          &                   slev, elev, opt_rlend=i_rlend      ) !in
      ENDIF


      ! during the last iteration step, the following computations can be skipped
      IF ( nsub == p_ncycl ) EXIT


      !
      ! 4.1/4.2 compute updated density and tracer fields
      !
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      i_startblk = p_patch%cells%start_blk(3,1)
      i_endblk   = p_patch%cells%end_blk(min_rlcell_int,i_nchdom)


    ! initialize nest boundary points at the second time level
    IF ( nsub == 1 .AND. (p_patch%id > 1 .OR. l_limited_area) ) THEN
      CALL copy(z_tracer(:,slev:elev,1:i_startblk,nnow), &
           z_tracer(:,slev:elev,1:i_startblk,nnew))
!$OMP BARRIER
    ENDIF

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_fluxdiv_c) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, 3, min_rlcell_int)

        ! compute mass flux divergence
        !
        ! This computation needs to be done only once, since the mass flux
        ! p_mass_flx_e is assumed to be constant in time.
        !
        IF ( nsub == 1 ) THEN
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
            DO jc = i_startidx, i_endidx
              DO jk = slev, elev
#else
!CDIR UNROLL=6
            DO jk = slev, elev
              DO jc = i_startidx, i_endidx
#endif

              z_rhofluxdiv_c(jc,jk,jb) =  &
                & p_mass_flx_e(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
                & p_mass_flx_e(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
                & p_mass_flx_e(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)

            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDIF

          ! compute tracer mass flux divergence
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = slev, elev
#else
!CDIR UNROLL=6
          DO jk = slev, elev
            DO jc = i_startidx, i_endidx
#endif

            z_fluxdiv_c(jc,jk) =  &
              & z_tracer_mflx(iidx(jc,jb,1),jk,iblk(jc,jb,1),nsub)*p_int%geofac_div(jc,1,jb) + &
              & z_tracer_mflx(iidx(jc,jb,2),jk,iblk(jc,jb,2),nsub)*p_int%geofac_div(jc,2,jb) + &
              & z_tracer_mflx(iidx(jc,jb,3),jk,iblk(jc,jb,3),nsub)*p_int%geofac_div(jc,3,jb)

          ENDDO
        ENDDO
!$ACC END PARALLEL

! Because the next loop requires z_fluxdiv_c it is crucial to synchronize with END PARALLEL

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slev, elev
          DO jc = i_startidx, i_endidx

            !
            ! 4.1 updated density field for intermediate timestep n + nsub/p_ncycl
            !
            z_rho(jc,jk,jb,nnew) = z_rho(jc,jk,jb,nnow)              &
              &                   - z_dtsub * z_rhofluxdiv_c(jc,jk,jb)

            !
            ! 4.2 updated tracer field for intermediate timestep n + nsub/p_ncycl
            !
            z_tracer(jc,jk,jb,nnew) = ( z_tracer(jc,jk,jb,nnow)       &
              &                      * z_rho(jc,jk,jb,nnow)           &
              &                      - z_dtsub * z_fluxdiv_c(jc,jk) ) &
              &                      / z_rho(jc,jk,jb,nnew)
          ENDDO
        ENDDO
!$ACC END PARALLEL

      ENDDO    ! loop over blocks

!$OMP ENDDO
!$OMP END PARALLEL

      nsav = nnow
      nnow = nnew
      nnew = nsav


      CALL sync_patch_array(SYNC_C,p_patch,z_tracer(:,:,:,nnow),opt_varname='z_tracer')


    ENDDO  ! loop over sub-timesteps



    !
    ! 5. compute averaged tracer mass flux
    !

    ! Before starting, preset halo edges that are not processed with zero's in order
    ! to avoid access of uninitialized array elements in subsequent routines
    ! Necessary when called within dycore

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

        ! Calculate flux at cell edge (cc_bary*v_{n}* \Delta p)
        !
        IF (p_ncycl == 2) THEN
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=5
          DO jk = slev, elev
            DO je = i_startidx, i_endidx
              p_out_e(je,jk,jb) = SUM(z_tracer_mflx(je,jk,jb,1:2))/REAL(p_ncycl,wp)
            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels
!$ACC END PARALLEL
        ELSE IF (p_ncycl == 3) THEN
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=5
          DO jk = slev, elev
            DO je = i_startidx, i_endidx
              p_out_e(je,jk,jb) = SUM(z_tracer_mflx(je,jk,jb,1:3))/REAL(p_ncycl,wp)
            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels
!$ACC END PARALLEL
        ENDIF

    ENDDO    ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_out_e ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

  END SUBROUTINE upwind_hflux_miura_cycl




  !-------------------------------------------------------------------------
  !>
  !! MIURA scheme with third order accurate lsq-reconstruction
  !!
  !! Calculation of time averaged horizontal tracer flux at triangle edges
  !! using a Flux form semi-lagrangian (FFSL)-method based on the MIURA-scheme.
  !! Unlike the standard Miura scheme, a third or fourth order accurate
  !! (i.e. quadratic, cubic) least squares reconstruction is applied.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-05-12)
  !!  Modification by Daniel Reinert, DWD (2010-10-14)
  !! - added possibility of cubic reconstruction
  !!
  !! @par !LITERATURE
  !! - Miura, H. (2007), Mon. Weather Rev., 135, 4038-4044
  !! - Ollivier-Gooch, C. (2002), JCP, 181, 729-752 (for lsq reconstruction)
  !! - Skamarock, W.C. (2010), Conservative Transport schemes for Spherical Geodesic
  !!   Grids: High-order Recosntructions for Forward-in-Time Schemes, Mon. Wea. Rev.
  !!
  SUBROUTINE upwind_hflux_miura3( p_patch, p_cc, p_mass_flx_e, p_vn, p_vt,          &
    &                        p_dtime, p_int, ld_compute, ld_cleanup,                &
    &                        p_itype_hlimit, p_out_e, opt_rhodz_now, opt_rhodz_new, &
    &                        opt_rlstart, opt_rlend, opt_lout_edge, opt_slev,       &
    &                        opt_elev, opt_ti_slev, opt_ti_elev )


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_hflux: upwind_hflux_miura3'

    TYPE(t_patch), INTENT(INOUT) ::  &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), INTENT(IN) ::    &    !< cell centered variable to be advected
      &  p_cc(:,:,:)                !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< contravariant horizontal mass flux at cell edge
      &  p_mass_flx_e(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::    &    !< normal component of velocity field
      &  p_vn(:,:,:)                !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::    &    !< tangential component of velocity field
      &  p_vt(:,:,:)                !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    LOGICAL, INTENT(IN) ::     &    !< flag, if .TRUE. compute geometrical terms
      &  ld_compute

    LOGICAL, INTENT(IN) ::     &    !< flag, if .TRUE. clean up geometrical terms
      &  ld_cleanup

    INTEGER, INTENT(IN) ::     &    !< parameter to select the limiter
      &  p_itype_hlimit             !< for horizontal transport

    REAL(wp), INTENT(INOUT) ::  &   !< output field, containing the tracer mass flux
      &  p_out_e(:,:,:)             !< or the reconstructed edge value depending on
                                    !< opt_lout_edge
                                    !< dim: (nproma,nlev,nblks_e)

    REAL(wp), OPTIONAL, INTENT(IN) :: &!< density times cell thickness at timestep n
      &  opt_rhodz_now(:,:,:)          !< dim: (nproma,nlev,nblks_c)

    REAL(wp), OPTIONAL, INTENT(IN) :: &!< density times cell thickness at timestep n+1
      &  opt_rhodz_new(:,:,:)          !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: output edge value (.TRUE.),
     & opt_lout_edge                   !< or the flux across the edge (.FALSE./not specified)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level (tracer independent part)
      &  opt_ti_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level (tracer independent part)
      &  opt_ti_elev

    LOGICAL  :: l_out_edgeval          !< corresponding local variable; default .FALSE.
                                       !< i.e. output flux across the edge

    REAL(wp) ::   &                     !< coefficients of lsq reconstruction
      &  z_lsq_coeff(lsq_high_set%dim_unk+1,nproma,p_patch%nlev,p_patch%nblks_c)
                                       !< at cell center
                                       !< includes c0 and gradients in zonal and
                                       !< meridional direction

    REAL(vp) ::  &                    !< coordinates of departure region vertices. The origin
      &  z_coords_dreg_v(nproma,4,2,p_patch%nlev,p_patch%nblks_e)
                                      !< of the coordinate system is at the circumcenter of
                                      !< the upwind cell. Unit vectors point to local East
                                      !< and North. (geographical coordinates)
                                      !< dim: (nproma,4,2,nlev,ptr_p%nblks_e)

    REAL(vp), ALLOCATABLE, SAVE ::  & !< gauss quadrature vector
      &  z_quad_vector_sum(:,:,:,:)   !< dim: (nproma,lsq_dim_unk+1,nlev,nblks_e)

    REAL(vp), ALLOCATABLE, SAVE ::  & !< area of departure region [m**2]
      &  z_dreg_area(:,:,:)           !< dim: (nproma,nlev,nblks_e)

    INTEGER, ALLOCATABLE, SAVE, TARGET ::  & !< line indices of upwind cell
      &  z_cell_idx(:,:,:)             !< dim: (nproma,nlev,p_patch%nblks_e)
    INTEGER, ALLOCATABLE, SAVE, TARGET ::  & !< block indices of upwind cell
      &  z_cell_blk(:,:,:)             !< dim: (nproma,nlev,p_patch%nblks_e)

    INTEGER, POINTER ::  &             !< Pointer to line and block indices of the cell
      &  ptr_ilc(:,:,:), ptr_ibc(:,:,:)!< center upstream of the edge

    INTEGER  :: nlev               !< number of full levels
    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: slev_ti, elev_ti   !< vertical start and end level (tracer independent part)
    INTEGER  :: ist                !< status variable
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    INTEGER  :: dim_unk
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_rlend_c, i_nchdom
    INTEGER  :: pid                !< patch ID

    TYPE(t_lsq), POINTER :: lsq_high !< Pointer to p_int_state%lsq_high
    lsq_high => p_int%lsq_high

   !-------------------------------------------------------------------------

    ! number of vertical levels
    nlev = p_patch%nlev

    ! get patch ID
    pid = p_patch%id

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = nlev
    END IF

    IF ( PRESENT(opt_ti_slev) ) THEN
      slev_ti = opt_ti_slev
    ELSE
      slev_ti = 1
    END IF
    IF ( PRESENT(opt_ti_elev) ) THEN
      elev_ti = opt_ti_elev
    ELSE
      elev_ti = nlev
    END IF

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 5
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

    i_rlend_c = min_rlcell_int

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

!$ACC DATA  PCOPYIN( p_cc, p_mass_flx_e, p_vn, p_vt ), PCOPYOUT( p_out_e ), &
!$ACC       CREATE( z_lsq_coeff, z_coords_dreg_v ), IF( i_am_accel_node .AND. acc_on )
!$ACC DATA  PCOPYIN( opt_rhodz_now ), IF( PRESENT(opt_rhodz_now) .AND. i_am_accel_node .AND. acc_on )
!$ACC DATA  PCOPYIN( opt_rhodz_new ), IF( PRESENT(opt_rhodz_new) .AND. i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_cc, p_mass_flx_e, p_vn, p_vt ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

    IF (p_test_run) THEN
!$ACC KERNELS IF (i_am_accel_node .AND. acc_on)
      z_lsq_coeff(:,:,:,:) = 0._wp
!$ACC END KERNELS
    ENDIF

    dim_unk = lsq_high_set%dim_unk+1

    !
    ! advection is done with an upwind scheme and a piecewise quadratic
    ! or cubic approximation of the tracer subgrid distribution.
    ! This approx. is integrated over a rhomboidal approximation of the
    ! departure region which is advected across the edge under consideration.
    ! The approximation is based on the (reconstructed) full 2D velocity
    ! field at edge midpoints (at time t+\Delta t/2) and \Delta t.
    !
    ! 3 options:  without limiter
    !             with monotone flux limiter following Zalesak (1979)
    !             with positive definite flux limiter following Zalesak (1979)
    !


    !
    ! 1. Approximation of the 'departure region'. The coordinates of
    !    all vertices are computed and stored in an edge-based data
    !    structure.
    !    In addition the Gauss-Legendre quadrature is prepared by
    !    calculating some tracer-invariant (i.e. purely geometric) fields.
    !
    IF ( ld_compute ) THEN
      ! allocate temporary arrays for quadrature and upwind cells
      ALLOCATE( z_quad_vector_sum(nproma,dim_unk,nlev,p_patch%nblks_e), &
        &       z_dreg_area(nproma,nlev,p_patch%nblks_e),               &
        &       z_cell_idx(nproma,nlev,p_patch%nblks_e),                &
        &       z_cell_blk(nproma,nlev,p_patch%nblks_e),                &
        &       STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                 &
          &  'allocation for z_quad_vector_sum, z_dreg_area, ' //    &
          &  'z_cell_idx, z_cell_blk' )
      ENDIF

!$ACC ENTER DATA CREATE( z_quad_vector_sum, z_dreg_area, z_cell_idx, z_cell_blk ), &
!$ACC       IF( i_am_accel_node .AND. acc_on )

      ! compute vertex coordinates for the departure region using a first
      ! order accurate (O(\Delta t)) backward trajectory-method
      CALL btraj_dreg( p_patch, p_int, p_vn, p_vt, p_dtime,      &! in
        &              .TRUE.,                                   &! in
        &              z_cell_idx, z_cell_blk, z_coords_dreg_v,  &! out
        &              opt_rlstart=i_rlstart, opt_rlend=i_rlend, &! in
        &              opt_slev=slev_ti, opt_elev=elev_ti        )! in



      ! maps quadrilateral onto the standard rectangle of edge length 2.
      ! provides quadrature points and the corresponding determinant of the
      ! Jacobian for each departure region.
      IF (lsq_high_ord == 2) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a quadratic 2D polynomial
        CALL prep_gauss_quadrature_q( p_patch, z_coords_dreg_v,           &! in
          &                      z_quad_vector_sum, z_dreg_area,          &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

      ELSE IF (lsq_high_ord == 30) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a cubic 2D polynomial without cross derivatives
        CALL prep_gauss_quadrature_cpoor( p_patch, z_coords_dreg_v,       &! in
          &                      z_quad_vector_sum, z_dreg_area,          &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

      ELSE IF (lsq_high_ord == 3) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a full cubic 2D polynomial
        CALL prep_gauss_quadrature_c( p_patch, z_coords_dreg_v,           &! in
          &                      z_quad_vector_sum, z_dreg_area,          &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in
      ENDIF

    END IF ! ld_compute

    !
    ! 2. reconstruction of the tracer subgrid distribution
    !    least squares method
    !    Note: for rlstart=2 we run into a sync-error with nests
    !
    IF (lsq_high_ord == 2) THEN
      ! quadratic reconstruction
      ! (computation of 6 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
        CALL recon_lsq_cell_q_svd( p_cc, p_patch, lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
        CALL recon_lsq_cell_q( p_cc, p_patch, lsq_high, z_lsq_coeff,        &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ELSE IF (lsq_high_ord == 30) THEN
      ! cubic reconstruction without cross derivatives
      ! (computation of 8 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_cpoor_svd( p_cc, p_patch, lsq_high, z_lsq_coeff,&
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_cpoor( p_cc, p_patch, lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ELSE IF (lsq_high_ord == 3) THEN
      ! cubic reconstruction with cross derivatives
      ! (computation of 10 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_c_svd( p_cc, p_patch, lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_c( p_cc, p_patch, lsq_high, z_lsq_coeff,        &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ENDIF


    ! Synchronize polynomial coefficients
    ! Note: a special sync routine is needed here because the fourth dimension
    ! of z_lsq_coeff is (for efficiency reasons) on the third index
    CALL sync_patch_array_4de1(SYNC_C1,p_patch,lsq_high_set%dim_unk+1,z_lsq_coeff,opt_varname='z_lsq_coeff 1')



    ! Pointer to line and block indices of the cell center upstream of the edge
    ptr_ilc => z_cell_idx(:,:,:)
    ptr_ibc => z_cell_blk(:,:,:)



    !
    ! 3. Calculate approximation to the area average \Phi_{avg} of the tracer
    !    in each rhomboidal area.
    !    Then calculate the flux v_n*\Delta p*\Phi_{avg}
    !    The fact that the rhomboidal area inevitably overlaps with neighboring
    !    triangles is neglected (local quadratic/cubic approximation instead
    !    of piecewise quadratic/cubic approximation). Only the reconstruction for
    !    the local cell is taken into account.

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    ! First of all, preset halo edges that are not processed with zero's in order
    ! to avoid access of uninitialized array elements in subsequent routines
    ! Necessary when called within dycore

    IF ( l_out_edgeval ) THEN
      i_startblk = p_patch%edges%start_blk(i_rlend-1,i_nchdom)
      i_endblk   = p_patch%edges%end_blk(min_rledge_int-3,i_nchdom)

      CALL init(p_out_e(:,:,i_startblk:i_endblk))
!$OMP BARRIER
    ENDIF

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    ! initialize also nest boundary points with zero
    IF ( l_out_edgeval .AND. (p_patch%id > 1 .OR. l_limited_area) ) THEN
      CALL init(p_out_e(:,:,1:i_startblk))
!$OMP BARRIER
    ENDIF

 !$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      IF ( l_out_edgeval ) THEN   ! Calculate 'edge value' of advected quantity

      ! Integral over departure region, normalized by departure region area
      ! (equals the tracer area average)
      ! - z_quad_vector_sum : tracer independent part
      ! - z_lsq_coeff       : tracer dependent part (lsq coefficients)

        SELECT  CASE( lsq_high_ord )
        CASE( 2 )  ! quadratic reconstruction

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=4
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=6
            p_out_e(je,jk,jb) =                                                       &
              &  DOT_PRODUCT(z_lsq_coeff(1:6,ptr_ilc(je,jk,jb),jk,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:6,jk,jb) ) / z_dreg_area(je,jk,jb)

          ENDDO
        ENDDO
!$ACC END PARALLEL

        CASE( 30 )  ! cubic reconstruction without third order cross derivatives

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=5
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=8
            p_out_e(je,jk,jb) =                                                       &
              &  DOT_PRODUCT(z_lsq_coeff(1:8,ptr_ilc(je,jk,jb),jk,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:8,jk,jb) ) / z_dreg_area(je,jk,jb)

          ENDDO
        ENDDO
!$ACC END PARALLEL

        CASE( 3 )  ! cubic reconstruction with third order cross derivatives

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=4
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=10
            p_out_e(je,jk,jb) =                                                        &
              &  DOT_PRODUCT(z_lsq_coeff(1:10,ptr_ilc(je,jk,jb),jk,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:10,jk,jb) ) / z_dreg_area(je,jk,jb)

          ENDDO
        ENDDO
!$ACC END PARALLEL

        END SELECT

      ELSE   ! Compute flux at cell edge (edge value * mass_flx)

       ! Calculate flux at cell edge
       !
       ! Integral over departure region, normalized by departure region area
       ! (equals the tracer area average) times the mass flux
       ! - z_quad_vector_sum : tracer independent part
       ! - z_lsq_coeff       : tracer dependent part (lsq coefficients)

        SELECT  CASE( lsq_high_ord )
        CASE( 2 )  ! quadratic reconstruction

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=4
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=6
            p_out_e(je,jk,jb) =                                                       &
              &  DOT_PRODUCT(z_lsq_coeff(1:6,ptr_ilc(je,jk,jb),jk,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:6,jk,jb) ) / z_dreg_area(je,jk,jb)            &
              &  * p_mass_flx_e(je,jk,jb)

          ENDDO
        ENDDO
!$ACC END PARALLEL

        CASE( 30 )  ! cubic reconstruction without third order cross derivatives

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=5
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=8
            p_out_e(je,jk,jb) =                                                       &
              &  DOT_PRODUCT(z_lsq_coeff(1:8,ptr_ilc(je,jk,jb),jk,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:8,jk,jb) ) / z_dreg_area(je,jk,jb)            &
              &  * p_mass_flx_e(je,jk,jb)

          ENDDO
        ENDDO
!$ACC END PARALLEL

        CASE( 3 )  ! cubic reconstruction with third order cross derivatives

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
!CDIR UNROLL=4
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=10
            p_out_e(je,jk,jb) =                                                        &
              &  DOT_PRODUCT(z_lsq_coeff(1:10,ptr_ilc(je,jk,jb),jk,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:10,jk,jb) ) / z_dreg_area(je,jk,jb)            &
              &  * p_mass_flx_e(je,jk,jb)

          ENDDO
        ENDDO
!$ACC END PARALLEL

        END SELECT

      ENDIF
    ENDDO  ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !
    ! 4. If desired, apply a monotonic or positive definite flux limiter 
    !    to limit computed fluxes.
    !    The flux limiter is based on work by Zalesak (1979)
    !
    IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_m) THEN
      IF (.NOT. (PRESENT(opt_rhodz_now) .OR. PRESENT(opt_rhodz_now))) THEN
        CALL finish(TRIM(routine),'Required fields opt_rhodz_now and opt_rhodz_now are not present')
      ENDIF
      !
      CALL hflx_limiter_mo( p_patch, p_int, p_dtime, p_cc,              & !in
        &                   opt_rhodz_now, opt_rhodz_new, p_mass_flx_e, & !in
        &                   p_out_e, slev, elev, opt_rlend=i_rlend,     & !inout,in
        &                   opt_beta_fct=advection_config(pid)%beta_fct ) !in
    ELSE IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_sm) THEN
      IF (.NOT. (PRESENT(opt_rhodz_now))) THEN
        CALL finish(TRIM(routine),'Required field opt_rhodz_now not present')
      ENDIF
      !
      CALL hflx_limiter_pd( p_patch, p_int, p_dtime,                 & !in
        &                   p_cc, opt_rhodz_now, p_out_e,            & !in,inout
        &                   slev, elev, opt_rlend=i_rlend            ) !in
    ENDIF


    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays for quadrature, departure region and
      ! upwind cell indices

!$ACC EXIT DATA DELETE( z_quad_vector_sum, z_dreg_area, z_cell_idx, z_cell_blk ), &
!$ACC           IF( i_am_accel_node .AND. acc_on )

      DEALLOCATE( z_quad_vector_sum, z_dreg_area, z_cell_idx, z_cell_blk, &
        &         STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                 &
          &  'deallocation for z_quad_vector_sum, z_dreg_area, ' //  &
          &  ' z_cell_idx, z_cell_blk failed' )
      ENDIF
    END IF

!$ACC UPDATE HOST( p_out_e ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA

  END SUBROUTINE upwind_hflux_miura3


  !-------------------------------------------------------------------------
  !>
  !! Flux-form semi Lagrangian scheme (extended MIURA3 scheme)
  !!
  !! Flux form semi Lagrangian scheme (extended MIURA3 scheme), where the overlap
  !! between the flux area and the underlying grid cells is taken into account.
  !! The scheme provides the time averaged horizontal tracer fluxes at triangle
  !! edges. A third or fourth order accurate (i.e. quadratic, cubic) least squares
  !! reconstruction can be selected.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-04-12)
  !!
  !! @par !LITERATURE
  !! - Miura, H. (2007), Mon. Weather Rev., 135, 4038-4044
  !! - Ollivier-Gooch, C. (2002), JCP, 181, 729-752 (for lsq reconstruction)
  !! - Skamarock, W.C. (2010), Conservative Transport schemes for Spherical Geodesic
  !!   Grids: High-order Reconstructions for Forward-in-Time Schemes, Mon. Wea. Rev.
  !!
  SUBROUTINE upwind_hflux_ffsl( p_patch, p_cc, p_rhodz_now, p_rhodz_new,     &
    &                      p_mass_flx_e, p_vn, p_vt,                         &
    &                      p_dtime, p_int, ld_compute, ld_cleanup,           &
    &                      p_itype_hlimit, p_out_e, opt_lconsv, opt_rlstart, &
    &                      opt_rlend, opt_lout_edge, opt_slev, opt_elev,     &
    &                      opt_ti_slev, opt_ti_elev  )


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_hflux: upwind_hflux_ffsl'

    TYPE(t_patch), INTENT(INOUT) ::  &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), INTENT(IN) ::    &    !< cell centered variable to be advected
      &  p_cc(:,:,:)                !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< density times cell thickness at timestep n
      &  p_rhodz_now(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< density times cell thickness at timestep n+1
      &  p_rhodz_new(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< contravariant horizontal mass flux at cell edge
      &  p_mass_flx_e(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::    &    !< normal component of velocity field
      &  p_vn(:,:,:)                !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::    &    !< tangential component of velocity field
      &  p_vt(:,:,:)                !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    LOGICAL, INTENT(IN) ::     &    !< flag, if .TRUE. compute geometrical terms
      &  ld_compute

    LOGICAL, INTENT(IN) ::     &    !< flag, if .TRUE. clean up geometrical terms
      &  ld_cleanup

    INTEGER, INTENT(IN) ::     &    !< parameter to select the limiter
      &  p_itype_hlimit             !< for horizontal transport

    REAL(wp), INTENT(INOUT) ::  &   !< output field, containing the tracer mass flux
      &  p_out_e(:,:,:)             !< or the reconstructed edge value
                                    !< dim: (nproma,nlev,nblks_e)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: if true, conservative reconstruction
     &  opt_lconsv                     !< is used (linear reconstruction only)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: output edge value (.TRUE.),
     & opt_lout_edge                   !< or the flux across the edge (.FALSE./not specified)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level (tracer independent part)
      &  opt_ti_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level (tracer independent part)
      &  opt_ti_elev

    LOGICAL  :: l_out_edgeval          !< corresponding local variable; default .FALSE.
                                       !< i.e. output flux across the edge

   REAL(wp) ::   &                     !< coefficients of lsq reconstruction
      &  z_lsq_coeff(lsq_high_set%dim_unk+1,nproma,p_patch%nlev,p_patch%nblks_c)
                                       !< at cell center
                                       !< includes c0 and gradients in zonal and
                                       !< meridional direction

    REAL(vp) ::  &                    !< patch 0,1,2 of subdivided departure region
      &  dreg_patch0(nproma,4,2,p_patch%nlev,p_patch%nblks_e), &  !< coordinates
      &  dreg_patch1(nproma,4,2,p_patch%nlev,p_patch%nblks_e), &
      &  dreg_patch2(nproma,4,2,p_patch%nlev,p_patch%nblks_e)


    REAL(vp), ALLOCATABLE, SAVE ::   & !< gauss quadrature vector for each patch
      &  z_quad_vector_sum0(:,:,:,:),& !< dim: (nproma,lsq_dim_unk+1,nlev,nblks_e)
      &  z_quad_vector_sum1(:,:,:,:),&
      &  z_quad_vector_sum2(:,:,:,:)

    REAL(vp), ALLOCATABLE, SAVE ::  & !< area of each departure region patch
      &  z_dreg_area0(:,:,:),       & !< dim: (nproma,nlev,nblks_e)
      &  z_dreg_area1(:,:,:),       &
      &  z_dreg_area2(:,:,:)

    INTEGER, ALLOCATABLE, SAVE, TARGET ::  & !< line and block indices of underlying cell
      &  patch0_cell_idx(:,:,:), patch0_cell_blk(:,:,:), & !< dim: (nproma,nlev,p_patch%nblks_e)
      &  patch1_cell_idx(:,:,:), patch1_cell_blk(:,:,:), &
      &  patch2_cell_idx(:,:,:), patch2_cell_blk(:,:,:)

    INTEGER, POINTER ::                    & !< Pointer to line and block indices of the cells
      &  ptr_ilc0(:,:,:), ptr_ibc0(:,:,:), & !< to which the departure region patches belong.
      &  ptr_ilc1(:,:,:), ptr_ibc1(:,:,:), &
      &  ptr_ilc2(:,:,:), ptr_ibc2(:,:,:)

    LOGICAL  :: l_consv            !< true if conservative lsq reconstruction is used
    INTEGER  :: nlev               !< number of full levels
    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: slev_ti, elev_ti   !< vertical start and end level (tracer independent part)
    INTEGER  :: ist                !< status variable
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    INTEGER  :: dim_unk
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_rlend_c, i_nchdom
    INTEGER  :: pid                !< patch ID
    TYPE(t_lsq), POINTER :: lsq_high !< Pointer to p_int_state%lsq_high
    lsq_high => p_int%lsq_high

   !-------------------------------------------------------------------------

    ! number of vertical levels
    nlev = p_patch%nlev

    ! get patch ID
    pid = p_patch%id

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = nlev
    END IF

    IF ( PRESENT(opt_ti_slev) ) THEN
      slev_ti = opt_ti_slev
    ELSE
      slev_ti = 1
    END IF
    IF ( PRESENT(opt_ti_elev) ) THEN
      elev_ti = opt_ti_elev
    ELSE
      elev_ti = nlev
    END IF

    IF ( PRESENT(opt_lconsv) ) THEN
     l_consv = opt_lconsv
    ELSE
     l_consv = .FALSE. ! non-conservative reconstruction
    ENDIF

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 5
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

    i_rlend_c = min_rlcell_int

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    IF (p_test_run) THEN
      z_lsq_coeff(:,:,:,:) = 0._wp
    ENDIF

    dim_unk = lsq_high_set%dim_unk+1

    !
    ! advection is done with an upwind scheme and a piecewise linear, quadratic
    ! or cubic approximation of the tracer subgrid distribution.
    ! This approx. is integrated over a rhomboidal approximation of the
    ! departure region which is advected across the edge under consideration.
    ! The approximation is based on the (reconstructed) full 2D velocity
    ! field at edge midpoints (at time t+\Delta t/2) and \Delta t.
    !
    ! 3 options:  without limiter
    !             with monotone flux limiter following Zalesak (1979)
    !             with positive definite flux limiter following Zalesak (1979)
    !

    !
    ! 1. Approximation of the 'departure region'. The coordinates of
    !    all vertices are computed and stored in an edge-based data
    !    structure.
    !    In addition the Gauss-Legendre quadrature is prepared by
    !    calculating some tracer-invariant (i.e. purely geometric) fields.
    !
    IF ( ld_compute ) THEN
      ! allocate temporary arrays for quadrature and upwind cells
      ALLOCATE( z_quad_vector_sum0(nproma,dim_unk,nlev,p_patch%nblks_e), &
        &       z_quad_vector_sum1(nproma,dim_unk,nlev,p_patch%nblks_e), &
        &       z_quad_vector_sum2(nproma,dim_unk,nlev,p_patch%nblks_e), &
        &       z_dreg_area0(nproma,nlev,p_patch%nblks_e),               &
        &       z_dreg_area1(nproma,nlev,p_patch%nblks_e),               &
        &       z_dreg_area2(nproma,nlev,p_patch%nblks_e),               &
        &       patch0_cell_idx(nproma,nlev,p_patch%nblks_e),            &
        &       patch1_cell_idx(nproma,nlev,p_patch%nblks_e),            &
        &       patch2_cell_idx(nproma,nlev,p_patch%nblks_e),            &
        &       patch0_cell_blk(nproma,nlev,p_patch%nblks_e),            &
        &       patch1_cell_blk(nproma,nlev,p_patch%nblks_e),            &
        &       patch2_cell_blk(nproma,nlev,p_patch%nblks_e),            &
        &       STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                        &
          &  'allocation for z_quad_vector_sum0/1/2, z_dreg_area0/1/2, ' // &
          &  'patch0/1/2_cell_idx,  patch0/1/2_cell_blk failed' )
      ENDIF


      ! compute vertex coordinates for the departure region using a first
      ! order accurate (O(\Delta t)) backward trajectory-method
      CALL btraj_dreg( p_patch, p_int, p_vn, p_vt, p_dtime,            &! in
        &              .FALSE.,                                        &! in
        &              patch0_cell_idx, patch0_cell_blk, dreg_patch0,  &! out
        &              opt_rlstart=i_rlstart, opt_rlend=i_rlend,       &! in
        &              opt_slev=slev_ti, opt_elev=elev_ti              )! in


      ! Flux area (aka. departure region) is subdivided according to its overlap
      ! with the underlying grid.
      CALL divide_flux_area(p_patch, p_int, p_vn, p_vt,                &! in
        &                   dreg_patch0, dreg_patch1, dreg_patch2,     &! out
        &                   patch1_cell_idx, patch1_cell_blk,          &! out
        &                   patch2_cell_idx, patch2_cell_blk,          &! out
        &                   opt_rlstart=i_rlstart, opt_rlend=i_rlend,  &! in
        &                   opt_slev=slev_ti, opt_elev=elev_ti         )! in



      ! maps quadrilateral onto the standard rectangle of edge length 2.
      ! provides quadrature points and the corresponding determinant of the
      ! Jacobian for each departure region.
      ! This is done for each of the three patch fragments.
      IF (lsq_high_ord == 1) THEN
        ! Gauss-Legendre quadrature with 1 quadrature point for integrating
        ! a linear 2D polynomial
        CALL prep_gauss_quadrature_l( p_patch, dreg_patch0,               &! in
          &                      z_quad_vector_sum0, z_dreg_area0,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_l( p_patch, dreg_patch1,               &! in
          &                      z_quad_vector_sum1, z_dreg_area1,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_l( p_patch, dreg_patch2,               &! in
          &                      z_quad_vector_sum2, z_dreg_area2,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

      ELSE IF (lsq_high_ord == 2) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a quadratic 2D polynomial
        CALL prep_gauss_quadrature_q( p_patch, dreg_patch0,               &! in
          &                      z_quad_vector_sum0, z_dreg_area0,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_q( p_patch, dreg_patch1,               &! in
          &                      z_quad_vector_sum1, z_dreg_area1,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_q( p_patch, dreg_patch2,               &! in
          &                      z_quad_vector_sum2, z_dreg_area2,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

      ELSE IF (lsq_high_ord == 30) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a cubic 2D polynomial without cross derivatives
        CALL prep_gauss_quadrature_cpoor( p_patch, dreg_patch0,           &! in
          &                      z_quad_vector_sum0, z_dreg_area0,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_cpoor( p_patch, dreg_patch1,           &! in
          &                      z_quad_vector_sum1, z_dreg_area1,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_cpoor( p_patch, dreg_patch2,           &! in
          &                      z_quad_vector_sum2, z_dreg_area2,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

      ELSE IF (lsq_high_ord == 3) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a full cubic 2D polynomial
        CALL prep_gauss_quadrature_c( p_patch, dreg_patch0,               &! in
          &                      z_quad_vector_sum0, z_dreg_area0,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_c( p_patch, dreg_patch1,               &! in
          &                      z_quad_vector_sum1, z_dreg_area1,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_c( p_patch, dreg_patch2,               &! in
          &                      z_quad_vector_sum2, z_dreg_area2,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

      ENDIF

    END IF ! ld_compute


    !
    ! 2. reconstruction of the tracer subgrid distribution
    !    least squares method
    !    Note: for rlstart=2 we run into a sync-error with nests
    !
    IF (lsq_high_ord == 1) THEN
      ! linear reconstruction
      ! (computation of 3 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
        CALL recon_lsq_cell_l_consv_svd( p_cc, p_patch, lsq_high, z_lsq_coeff,    &
          &                              opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
          &                              opt_rlstart=2, opt_lconsv=l_consv )
      ELSE
        CALL recon_lsq_cell_l( p_cc, p_patch, lsq_high, z_lsq_coeff,        &
          &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
          &                    opt_rlstart=2, opt_lconsv=l_consv )
      ENDIF
    ELSE IF (lsq_high_ord == 2) THEN
      ! quadratic reconstruction
      ! (computation of 6 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_q_svd( p_cc, p_patch, lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_q( p_cc, p_patch, lsq_high, z_lsq_coeff,        &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ELSE IF (lsq_high_ord == 30) THEN
      ! cubic reconstruction without cross derivatives
      ! (computation of 8 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_cpoor_svd( p_cc, p_patch, lsq_high, z_lsq_coeff,&
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_cpoor( p_cc, p_patch, lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ELSE IF (lsq_high_ord == 3) THEN
      ! cubic reconstruction with cross derivatives
      ! (computation of 10 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_c_svd( p_cc, p_patch, lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_c( p_cc, p_patch, lsq_high, z_lsq_coeff,        &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ENDIF


    ! Synchronize polynomial coefficients
    ! Note: a special sync routine is needed here because the fourth dimension
    ! of z_lsq_coeff is (for efficiency reasons) on the third index
    CALL sync_patch_array_4de1(SYNC_C,p_patch,lsq_high_set%dim_unk+1,z_lsq_coeff,opt_varname='z_lsq_coeff 2')




    ! Pointer to line and block indices of the cells to which the departure
    ! region patches belong.
    ptr_ilc0 => patch0_cell_idx(:,:,:)
    ptr_ibc0 => patch0_cell_blk(:,:,:)
    ptr_ilc1 => patch1_cell_idx(:,:,:)
    ptr_ibc1 => patch1_cell_blk(:,:,:)
    ptr_ilc2 => patch2_cell_idx(:,:,:)
    ptr_ibc2 => patch2_cell_blk(:,:,:)

    !
    ! 3. Calculate approximation to the area average \Phi_{avg} of the tracer
    !    in each rhomboidal area.
    !    Then calculate the flux v_n*\Delta p*\Phi_{avg}
    !    The fact that the rhomboidal area inevitably overlaps with neighboring
    !    triangles is at least partly taken into account.
    !
    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)


      ! Calculate flux at cell edge
      !
      ! Integral over departure region, normalized by departure region area
      ! (equals the tracer area average) times the mass flux
      ! - z_quad_vector_sum : tracer independent part
      ! - z_lsq_coeff       : tracer dependent part (lsq coefficients)

      SELECT  CASE( lsq_high_ord )
      CASE( 1 )  ! linear reconstruction

!CDIR UNROLL=5
      DO jk = slev, elev
        DO je = i_startidx, i_endidx

!CDIR EXPAND=3
          p_out_e(je,jk,jb) =                                                            &
            &   ( DOT_PRODUCT(z_lsq_coeff(1:3,ptr_ilc0(je,jk,jb),jk,ptr_ibc0(je,jk,jb)), &
            &     z_quad_vector_sum0(je,1:3,jk,jb) )                                     &
            &   + DOT_PRODUCT(z_lsq_coeff(1:3,ptr_ilc1(je,jk,jb),jk,ptr_ibc1(je,jk,jb)), &
            &     z_quad_vector_sum1(je,1:3,jk,jb) )                                     &
            &   + DOT_PRODUCT(z_lsq_coeff(1:3,ptr_ilc2(je,jk,jb),jk,ptr_ibc2(je,jk,jb)), &
            &     z_quad_vector_sum2(je,1:3,jk,jb) ) )                                   &
            &   / (z_dreg_area0(je,jk,jb)+z_dreg_area1(je,jk,jb)+z_dreg_area2(je,jk,jb) )&
            &   * p_mass_flx_e(je,jk,jb)

        ENDDO
      ENDDO

      CASE( 2 )  ! quadratic reconstruction

!CDIR UNROLL=5
      DO jk = slev, elev
        DO je = i_startidx, i_endidx

!CDIR EXPAND=6
          p_out_e(je,jk,jb) =                                                            &
            &   ( DOT_PRODUCT(z_lsq_coeff(1:6,ptr_ilc0(je,jk,jb),jk,ptr_ibc0(je,jk,jb)), &
            &     z_quad_vector_sum0(je,1:6,jk,jb) )                                     &
            &   + DOT_PRODUCT(z_lsq_coeff(1:6,ptr_ilc1(je,jk,jb),jk,ptr_ibc1(je,jk,jb)), &
            &     z_quad_vector_sum1(je,1:6,jk,jb) )                                     &
            &   + DOT_PRODUCT(z_lsq_coeff(1:6,ptr_ilc2(je,jk,jb),jk,ptr_ibc2(je,jk,jb)), &
            &     z_quad_vector_sum2(je,1:6,jk,jb) ) )                                   &
            &   / (z_dreg_area0(je,jk,jb)+z_dreg_area1(je,jk,jb)+z_dreg_area2(je,jk,jb) )&
            &   * p_mass_flx_e(je,jk,jb)

        ENDDO
      ENDDO

      CASE( 30 )  ! cubic reconstruction without third order cross derivatives

      DO jk = slev, elev
        DO je = i_startidx, i_endidx

!CDIR EXPAND=8
          p_out_e(je,jk,jb) =                                                            &
            &   ( DOT_PRODUCT(z_lsq_coeff(1:8,ptr_ilc0(je,jk,jb),jk,ptr_ibc0(je,jk,jb)), &
            &     z_quad_vector_sum0(je,1:8,jk,jb) )                                     &
            &   + DOT_PRODUCT(z_lsq_coeff(1:8,ptr_ilc1(je,jk,jb),jk,ptr_ibc1(je,jk,jb)), &
            &     z_quad_vector_sum1(je,1:8,jk,jb) )                                     &
            &   + DOT_PRODUCT(z_lsq_coeff(1:8,ptr_ilc2(je,jk,jb),jk,ptr_ibc2(je,jk,jb)), &
            &     z_quad_vector_sum2(je,1:8,jk,jb) ) )                                   &
            &   / (z_dreg_area0(je,jk,jb)+z_dreg_area1(je,jk,jb)+z_dreg_area2(je,jk,jb) )&
            &   * p_mass_flx_e(je,jk,jb)

        ENDDO
      ENDDO

      CASE( 3 )  ! cubic reconstruction with third order cross derivatives

      DO jk = slev, elev
        DO je = i_startidx, i_endidx

!CDIR EXPAND=10
          p_out_e(je,jk,jb) =                                                            &
            &   ( DOT_PRODUCT(z_lsq_coeff(1:10,ptr_ilc0(je,jk,jb),jk,ptr_ibc0(je,jk,jb)),&
            &     z_quad_vector_sum0(je,1:10,jk,jb) )                                    &
            &   + DOT_PRODUCT(z_lsq_coeff(1:10,ptr_ilc1(je,jk,jb),jk,ptr_ibc1(je,jk,jb)),&
            &     z_quad_vector_sum1(je,1:10,jk,jb) )                                    &
            &   + DOT_PRODUCT(z_lsq_coeff(1:10,ptr_ilc2(je,jk,jb),jk,ptr_ibc2(je,jk,jb)),&
            &     z_quad_vector_sum2(je,1:10,jk,jb) ) )                                  &
            &   / (z_dreg_area0(je,jk,jb)+z_dreg_area1(je,jk,jb)+z_dreg_area2(je,jk,jb) )&
            &   * p_mass_flx_e(je,jk,jb)


        ENDDO
      ENDDO

      END SELECT

    ENDDO  ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


    !
    ! 4. If desired, apply a monotonic or positive definite flux limiter 
    !    to limit computed fluxes.
    !    The flux limiter is based on work by Zalesak (1979)
    !
    IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_m) THEN
      CALL hflx_limiter_mo( p_patch, p_int, p_dtime, p_cc,              & !in
        &                   p_rhodz_now, p_rhodz_new, p_mass_flx_e,     & !in
        &                   p_out_e, slev, elev, opt_rlend=i_rlend,     & !inout,in
        &                   opt_beta_fct=advection_config(pid)%beta_fct ) !in
    ELSE IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_sm) THEN
      !
      CALL hflx_limiter_pd( p_patch, p_int, p_dtime,                 & !in
        &                   p_cc, p_rhodz_now, p_out_e,              & !in,inout
        &                   slev, elev, opt_rlend=i_rlend            ) !in
    ENDIF


    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays for quadrature, departure region and
      ! upwind cell indices
      DEALLOCATE( z_quad_vector_sum0, z_quad_vector_sum1, z_quad_vector_sum2, &
        &         z_dreg_area0, z_dreg_area1, z_dreg_area2, patch0_cell_idx,  &
        &         patch1_cell_idx, patch2_cell_idx, patch0_cell_blk,          &
        &         patch1_cell_blk, patch2_cell_blk, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                           &
          &  'deallocation for z_quad_vector_sum0/1/2, z_dreg_area0/1/2, '  // &
          &  'patch0/1/2_cell_idx, patch0/1/2_cell_blk failed' )
      ENDIF
    END IF


  END SUBROUTINE upwind_hflux_ffsl




  !-------------------------------------------------------------------------
  !>
  !! Hybrid FFSL/MIURA3-scheme
  !!
  !! Flux form semi Lagrangian scheme (extended MIURA3 scheme), where the overlap
  !! between the flux area and the underlying grid cells is taken into account.
  !! The scheme provides the time averaged horizontal tracer fluxes at triangle
  !! edges. A third or fourth order accurate (i.e. quadratic, cubic) least squares
  !! reconstruction can be selected.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2013-11-01)
  !!
  !! @par !LITERATURE
  !! - Miura, H. (2007), Mon. Weather Rev., 135, 4038-4044
  !! - Ollivier-Gooch, C. (2002), JCP, 181, 729-752 (for lsq reconstruction)
  !! - Skamarock, W.C. (2010), Conservative Transport schemes for Spherical Geodesic
  !!   Grids: High-order Recosntructions for Forward-in-Time Schemes, Mon. Wea. Rev.
  !!
  SUBROUTINE hflux_ffsl_hybrid( p_patch, p_cc, p_rhodz_now, p_rhodz_new,     &
    &                      p_mass_flx_e, p_vn, p_vt,                         &
    &                      p_dtime, p_int, ld_compute, ld_cleanup,           &
    &                      p_itype_hlimit, p_out_e, opt_lconsv, opt_rlstart, &
    &                      opt_rlend, opt_lout_edge, opt_slev, opt_elev,     &
    &                      opt_ti_slev, opt_ti_elev  )


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_hflux: upwind_hflux_hybrid'

    TYPE(t_patch), INTENT(INOUT) ::  &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), INTENT(IN) ::    &    !< cell centered variable to be advected
      &  p_cc(:,:,:)                !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< density times cell thickness at timestep n
      &  p_rhodz_now(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< density times cell thickness at timestep n+1
      &  p_rhodz_new(:,:,:)         !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< contravariant horizontal mass flux at cell edge
      &  p_mass_flx_e(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::    &    !< normal component of velocity field
      &  p_vn(:,:,:)                !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::    &    !< tangential component of velocity field
      &  p_vt(:,:,:)                !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    LOGICAL, INTENT(IN) ::     &    !< flag, if .TRUE. compute geometrical terms
      &  ld_compute

    LOGICAL, INTENT(IN) ::     &    !< flag, if .TRUE. clean up geometrical terms
      &  ld_cleanup

    INTEGER, INTENT(IN) ::     &    !< parameter to select the limiter
      &  p_itype_hlimit             !< for horizontal transport

    REAL(wp), INTENT(INOUT) ::  &   !< output field, containing the tracer mass flux
      &  p_out_e(:,:,:)             !< or the reconstructed edge value
                                    !< dim: (nproma,nlev,nblks_e)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: if true, conservative reconstruction
     &  opt_lconsv                     !< is used (linear reconstruction only)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: output edge value (.TRUE.),
     & opt_lout_edge                   !< or the flux across the edge (.FALSE./not specified)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level (tracer independent part)
      &  opt_ti_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level (tracer independent part)
      &  opt_ti_elev

    LOGICAL  :: l_out_edgeval          !< corresponding local variable; default .FALSE.
                                       !< i.e. output flux across the edge

   REAL(wp) ::   &                     !< coefficients of lsq reconstruction
      &  z_lsq_coeff(lsq_high_set%dim_unk+1,nproma,p_patch%nlev,p_patch%nblks_c)
                                       !< at cell center
                                       !< includes c0 and gradients in zonal and
                                       !< meridional direction

    REAL(vp) ::  &                    !< patch 0,1,2 of subdivided departure region
      &  dreg_patch0(nproma,4,2,p_patch%nlev,p_patch%nblks_e)  !< coordinates

    REAL(vp), ALLOCATABLE ::   & !< dim: (npoints,4,2,nblks_e)
      &  dreg_patch1(:,:,:,:), &
      &  dreg_patch2(:,:,:,:)


    REAL(vp), ALLOCATABLE, SAVE ::   & !< gauss quadrature vector for each patch
      &  z_quad_vector_sum0(:,:,:,:),& !< dim: (nproma,lsq_dim_unk+1,nlev,nblks_e)
      &  z_quad_vector_sum1(:,:,:),  & !< dim: (npoints,lsq_dim_unk+1,nblks_e)
      &  z_quad_vector_sum2(:,:,:)

    REAL(vp), ALLOCATABLE, SAVE ::  & !< sum of area of departure region patches
      &  z_dreg_area(:,:,:)

    INTEGER, ALLOCATABLE, SAVE  ::  & !< line and block indices of underlying cell
      &  patch0_cell_idx(:,:,:), patch0_cell_blk(:,:,:), & !< dim: (nproma,nlev,p_patch%nblks_e)
      &  patch1_cell_idx(:,:),   patch1_cell_blk(:,:),   & !< dim: (npoints,p_patch%nblks_e)
      &  patch2_cell_idx(:,:),   patch2_cell_blk(:,:)

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_lsq_coeff,dreg_patch0,dreg_patch1,dreg_patch2
!DIR$ ATTRIBUTES ALIGN : 64 :: z_quad_vector_sum0,z_quad_vector_sum1,z_quad_vector_sum2
!DIR$ ATTRIBUTES ALIGN : 64 :: z_dreg_area
!DIR$ ATTRIBUTES ALIGN : 64 :: patch0_cell_idx,patch1_cell_idx,patch2_cell_idx
!DIR$ ATTRIBUTES ALIGN : 64 :: patch0_cell_blk,patch1_cell_blk,patch2_cell_blk
#endif

    TYPE(t_list2D), SAVE ::   &    !< list with points for which a local
      &  falist                    !< polynomial approximation is insufficient
                                   !< and a piecewise approximation is needed,
                                   !< instead

    INTEGER  :: nlev               !< number of full levels
    INTEGER  :: npoints            !< number of points per block for ndex list allocation
    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: slev_ti, elev_ti   !< vertical start and end level (tracer independent part)
    INTEGER  :: ist                !< status variable
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    INTEGER  :: ie                 !< index list loop counter
    INTEGER  :: dim_unk
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_rlend_c, i_nchdom
    INTEGER  :: pid                !< patch ID

    TYPE(t_lsq), POINTER :: lsq_high !< Pointer to p_int_state%lsq_high
    lsq_high => p_int%lsq_high

   !-------------------------------------------------------------------------

    ! number of vertical levels
    nlev = p_patch%nlev

    ! get patch ID
    pid = p_patch%id

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = nlev
    END IF

    IF ( PRESENT(opt_ti_slev) ) THEN
      slev_ti = opt_ti_slev
    ELSE
      slev_ti = 1
    END IF
    IF ( PRESENT(opt_ti_elev) ) THEN
      elev_ti = opt_ti_elev
    ELSE
      elev_ti = nlev
    END IF

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 5
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

    i_rlend_c = min_rlcell_int

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    IF (p_test_run) THEN
      z_lsq_coeff(:,:,:,:) = 0._wp
    ENDIF

    dim_unk = lsq_high_set%dim_unk+1

    !
    ! advection is done with an upwind scheme and a piecewise linear, quadratic
    ! or cubic approximation of the tracer subgrid distribution.
    ! This approx. is integrated over a rhomboidal approximation of the
    ! departure region which is advected across the edge under consideration.
    ! The approximation is based on the (reconstructed) full 2D velocity
    ! field at edge midpoints (at time t+\Delta t/2) and \Delta t.
    !
    ! 3 options:  without limiter
    !             with monotone flux limiter following Zalesak (1979)
    !             with positive definite flux limiter following Zalesak (1979)
    !

    !
    ! 1. Approximation of the 'departure region'. The coordinates of
    !    all vertices are computed and stored in an edge-based data
    !    structure.
    !    In addition the Gauss-Legendre quadrature is prepared by
    !    calculating some tracer-invariant (i.e. purely geometric) fields.
    !

    IF ( ld_compute ) THEN
      ! allocate temporary arrays for quadrature and upwind cells
      ALLOCATE( patch0_cell_idx(nproma,nlev,p_patch%nblks_e),            &
        &       patch0_cell_blk(nproma,nlev,p_patch%nblks_e),            &
        &       falist%eidx(nproma*nlev,p_patch%nblks_e),                &
        &       falist%elev(nproma*nlev,p_patch%nblks_e),                &
        &       falist%len(p_patch%nblks_e),                             &
        &       STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                        &
          &  'allocation for patch0_cell_idx,  patch0_cell_blk, falist failed' )
      ENDIF

      falist%len(:) = 0 ! To be safe with the subsequent MAXVAL in the case of nesting

      ! compute vertex coordinates for the departure region using a first
      ! order accurate (O(\Delta t)) backward trajectory-method
      CALL btraj_dreg( p_patch, p_int, p_vn, p_vt, p_dtime,            &! in
        &              .FALSE.,                                        &! in
        &              patch0_cell_idx, patch0_cell_blk, dreg_patch0,  &! out
        &              opt_rlstart=i_rlstart, opt_rlend=i_rlend,       &! in
        &              opt_slev=slev_ti, opt_elev=elev_ti,             &! in
        &              opt_falist=falist                               )! inout

      npoints = MAXVAL(falist%len(:))
      falist%npoints = npoints

      ! allocate temporary arrays for quadrature and upwind cells
      ALLOCATE( z_quad_vector_sum0(nproma,dim_unk,nlev,p_patch%nblks_e), &
        &       z_quad_vector_sum1(npoints,dim_unk,p_patch%nblks_e),     &
        &       z_quad_vector_sum2(npoints,dim_unk,p_patch%nblks_e),     &
        &       z_dreg_area(nproma,nlev,p_patch%nblks_e),                &
        &       patch1_cell_idx(npoints,p_patch%nblks_e),                &
        &       patch2_cell_idx(npoints,p_patch%nblks_e),                &
        &       patch1_cell_blk(npoints,p_patch%nblks_e),                &
        &       patch2_cell_blk(npoints,p_patch%nblks_e),                &
        &       dreg_patch1(npoints,4,2,p_patch%nblks_e),                &
        &       dreg_patch2(npoints,4,2,p_patch%nblks_e),                &
        &       STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                      &
          &  'allocation for z_quad_vector_sum0/1/2, z_dreg_area, ' //    &
          &  'patch1/2_cell_idx,  patch1/2_cell_blk, dreg_patch1/2 failed' )
      ENDIF

      ! Flux area (aka. departure region) is subdivided according to its overlap
      ! with the underlying grid.
      CALL divide_flux_area_list(p_patch, p_int, p_vn, p_vt, falist,   &! in
        &                   dreg_patch0, dreg_patch1, dreg_patch2,     &! out
        &                   patch1_cell_idx, patch1_cell_blk,          &! out
        &                   patch2_cell_idx, patch2_cell_blk,          &! out
        &                   opt_rlstart=i_rlstart, opt_rlend=i_rlend   )! in

      ! maps quadrilateral onto the standard rectangle of edge length 2.
      ! provides quadrature points and the corresponding determinant of the
      ! Jacobian for each departure region.
      ! This is done for each of the three patch fragments.
      IF (lsq_high_ord == 1) THEN
        ! Gauss-Legendre quadrature with 1 quadrature point for integrating
        ! a linear 2D polynomial
        CALL prep_gauss_quadrature_l( p_patch, dreg_patch0,               &! in
          &                      z_quad_vector_sum0, z_dreg_area,         &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_l_list( p_patch, dreg_patch1, falist,  &! in
          &                      z_quad_vector_sum1, z_dreg_area,         &! out/inout
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend )! in

        CALL prep_gauss_quadrature_l_list( p_patch, dreg_patch2, falist,  &! in
          &                      z_quad_vector_sum2, z_dreg_area,         &! out/inout
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend )! in

      ELSE IF (lsq_high_ord == 2) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a quadratic 2D polynomial
        CALL prep_gauss_quadrature_q( p_patch, dreg_patch0,               &! in
          &                      z_quad_vector_sum0, z_dreg_area,         &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_q_list( p_patch, dreg_patch1, falist,  &! in
          &                      z_quad_vector_sum1, z_dreg_area,         &! out/inout
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend )! in

        CALL prep_gauss_quadrature_q_list( p_patch, dreg_patch2, falist,  &! in
          &                      z_quad_vector_sum2, z_dreg_area,         &! out/inout
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend )! in

      ELSE IF (lsq_high_ord == 3) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a full cubic 2D polynomial
        CALL prep_gauss_quadrature_c( p_patch, dreg_patch0,               &! in
          &                      z_quad_vector_sum0, z_dreg_area,         &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev_ti, opt_elev=elev_ti       )! in

        CALL prep_gauss_quadrature_c_list( p_patch, dreg_patch1, falist,  &! in
          &                      z_quad_vector_sum1, z_dreg_area,         &! out/inout
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend )! in

        CALL prep_gauss_quadrature_c_list( p_patch, dreg_patch2, falist,  &! in
          &                      z_quad_vector_sum2, z_dreg_area,         &! out/inout
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend )! in

      ENDIF

      DEALLOCATE(dreg_patch1, dreg_patch2)

    END IF ! ld_compute


    !
    ! 2. reconstruction of the tracer subgrid distribution
    !    least squares method
    !    Note: for rlstart=2 we run into a sync-error with nests
    !
    IF (lsq_high_ord == 1) THEN
      ! linear reconstruction
      ! (computation of 3 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
        CALL recon_lsq_cell_l_consv_svd( p_cc, p_patch, lsq_high, z_lsq_coeff,    &
          &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
          &                    opt_rlstart=2 )
      ELSE
        CALL recon_lsq_cell_l( p_cc, p_patch, lsq_high, z_lsq_coeff,        &
          &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
          &                    opt_rlstart=2 )
      ENDIF
    ELSE IF (lsq_high_ord == 2) THEN
      ! quadratic reconstruction
      ! (computation of 6 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_q_svd( p_cc, p_patch, lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_q( p_cc, p_patch, lsq_high, z_lsq_coeff,        &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ELSE IF (lsq_high_ord == 3) THEN
      ! cubic reconstruction with cross derivatives
      ! (computation of 10 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_c_svd( p_cc, p_patch, lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_c( p_cc, p_patch, lsq_high, z_lsq_coeff,        &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ENDIF


    ! Synchronize polynomial coefficients
    ! Note: a special sync routine is needed here because the fourth dimension
    ! of z_lsq_coeff is (for efficiency reasons) on the third index
    CALL sync_patch_array_4de1(SYNC_C,p_patch,lsq_high_set%dim_unk+1,z_lsq_coeff,opt_varname='z_lsq_coeff 3')

    !
    ! 3. Calculate approximation to the area average \Phi_{avg} of the tracer
    !    in each rhomboidal area.
    !    Then calculate the flux v_n*\Delta p*\Phi_{avg}
    !    The fact that the rhomboidal area inevitably overlaps with neighboring
    !    triangles is at least partly taken into account.
    !
    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,ie,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)


      ! Calculate flux at cell edge
      !
      ! Integral over departure region, normalized by departure region area
      ! (equals the tracer area average) times the mass flux
      ! - z_quad_vector_sum : tracer independent part
      ! - z_lsq_coeff       : tracer dependent part (lsq coefficients)

      SELECT  CASE( lsq_high_ord )
      CASE( 1 )  ! linear reconstruction

!CDIR UNROLL=5
      DO jk = slev, elev
        DO je = i_startidx, i_endidx

!CDIR EXPAND=3
          p_out_e(je,jk,jb) =                                                            &
            &     DOT_PRODUCT(z_lsq_coeff(1:3,patch0_cell_idx(je,jk,jb),jk,patch0_cell_blk(je,jk,jb)), &
            &     z_quad_vector_sum0(je,1:3,jk,jb))
        ENDDO
      ENDDO

      ! Correction for points in index list
!CDIR NODEP,VOVERTAKE,VOB
      DO ie = 1, falist%len(jb)

        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)

!CDIR EXPAND=3
        p_out_e(je,jk,jb) = p_out_e(je,jk,jb)                                            &
          &    + DOT_PRODUCT(z_lsq_coeff(1:3,patch1_cell_idx(ie,jb),jk,patch1_cell_blk(ie,jb)),  &
          &      z_quad_vector_sum1(ie,1:3,jb) )                                      &
          &    + DOT_PRODUCT(z_lsq_coeff(1:3,patch2_cell_idx(ie,jb),jk,patch2_cell_blk(ie,jb)),  &
          &      z_quad_vector_sum2(ie,1:3,jb) )
      ENDDO  ! ie


      CASE( 2 )  ! quadratic reconstruction

!CDIR UNROLL=5
      DO jk = slev, elev
        DO je = i_startidx, i_endidx

!CDIR EXPAND=6
          p_out_e(je,jk,jb) =                                                            &
            &     DOT_PRODUCT(z_lsq_coeff(1:6,patch0_cell_idx(je,jk,jb),jk,patch0_cell_blk(je,jk,jb)), &
            &     z_quad_vector_sum0(je,1:6,jk,jb))
        ENDDO
      ENDDO

      ! Correction for points in index list
!CDIR NODEP,VOVERTAKE,VOB
      DO ie = 1, falist%len(jb)

        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)

!CDIR EXPAND=6
        p_out_e(je,jk,jb) = p_out_e(je,jk,jb)                                            &
          &    + DOT_PRODUCT(z_lsq_coeff(1:6,patch1_cell_idx(ie,jb),jk,patch1_cell_blk(ie,jb)),  &
          &      z_quad_vector_sum1(ie,1:6,jb) )                                      &
          &    + DOT_PRODUCT(z_lsq_coeff(1:6,patch2_cell_idx(ie,jb),jk,patch2_cell_blk(ie,jb)),  &
          &      z_quad_vector_sum2(ie,1:6,jb) )
      ENDDO  ! ie


      CASE( 3 )  ! cubic reconstruction with third order cross derivatives

      DO jk = slev, elev
        DO je = i_startidx, i_endidx

!CDIR EXPAND=10
          p_out_e(je,jk,jb) =                                                            &
            &     DOT_PRODUCT(z_lsq_coeff(1:10,patch0_cell_idx(je,jk,jb),jk,patch0_cell_blk(je,jk,jb)),&
            &     z_quad_vector_sum0(je,1:10,jk,jb))
        ENDDO
      ENDDO

      ! Correction for points in index list
!CDIR NODEP,VOVERTAKE,VOB
      DO ie = 1, falist%len(jb)

        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)

!CDIR EXPAND=10
        p_out_e(je,jk,jb) = p_out_e(je,jk,jb)  &
          &    + DOT_PRODUCT(z_lsq_coeff(1:10,patch1_cell_idx(ie,jb),jk,patch1_cell_blk(ie,jb)),&
          &      z_quad_vector_sum1(ie,1:10,jb) )                                    &
          &    + DOT_PRODUCT(z_lsq_coeff(1:10,patch2_cell_idx(ie,jb),jk,patch2_cell_blk(ie,jb)),&
          &      z_quad_vector_sum2(ie,1:10,jb) )
      ENDDO  ! ie

      END SELECT

      ! Finally compute total flux
!CDIR UNROLL=5
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
          p_out_e(je,jk,jb) = p_mass_flx_e(je,jk,jb) * p_out_e(je,jk,jb) / z_dreg_area(je,jk,jb)
        ENDDO  ! je
      ENDDO  ! jk

    ENDDO  ! loop over blocks
!$OMP END DO
!$OMP END PARALLEL


    !
    ! 4. If desired, apply a monotonic or positive definite flux limiter 
    !    to limit computed fluxes.
    !    The flux limiter is based on work by Zalesak (1979)
    !
    IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_m) THEN
      CALL hflx_limiter_mo( p_patch, p_int, p_dtime, p_cc,              & !in
        &                   p_rhodz_now, p_rhodz_new, p_mass_flx_e,     & !in
        &                   p_out_e, slev, elev, opt_rlend=i_rlend,     & !inout,in
        &                   opt_beta_fct=advection_config(pid)%beta_fct ) !in
    ELSE IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_sm) THEN
      !
      CALL hflx_limiter_pd( p_patch, p_int, p_dtime,                 & !in
        &                   p_cc, p_rhodz_now, p_out_e,              & !in,inout
        &                   slev, elev, opt_rlend=i_rlend            ) !in
    ENDIF


    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays for quadrature, departure region and
      ! upwind cell indices
      DEALLOCATE( z_quad_vector_sum0, z_quad_vector_sum1, z_quad_vector_sum2, &
        &         z_dreg_area, patch0_cell_idx,                               &
        &         patch1_cell_idx, patch2_cell_idx, patch0_cell_blk,          &
        &         patch1_cell_blk, patch2_cell_blk, falist%eidx, falist%elev, &
        &         falist%len, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                           &
          &  'deallocation for z_quad_vector_sum0/1/2, z_dreg_area0/1/2, '  // &
          &  'patch0/1/2_cell_idx, patch0/1/2_cell_blk, falist failed' )
      ENDIF
    END IF


  END SUBROUTINE hflux_ffsl_hybrid


END MODULE mo_advection_hflux

