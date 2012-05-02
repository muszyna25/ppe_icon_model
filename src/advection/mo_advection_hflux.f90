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
!! Modification by Daniel reinert, DWD (2011-09-20)
!! - new Miura-type advection scheme with internal time-step subcycling
!!
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_advection_hflux

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, SUCCESS, TRACER_ONLY,      &
    &                               min_rledge_int, min_rledge, min_rlcell_int, &
    &                               UP, MIURA, MIURA3, FFSL, MCYCL,             &
    &                               MIURA_MCYCL, MIURA3_MCYCL, UP3, islopel_sm, &
    &                               islopel_m, ifluxl_m, ifluxl_sm,             &
    &                               INH_ATMOSPHERE, IHS_ATM_TEMP,               &
    &                               ISHALLOW_WATER, IHS_ATM_THETA
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_math_gradients,      ONLY: grad_green_gauss_cell
  USE mo_math_laplace,        ONLY: directional_laplace
  USE mo_math_divrot,         ONLY: recon_lsq_cell_l, recon_lsq_cell_q,         &
    &                               recon_lsq_cell_cpoor, recon_lsq_cell_c,     &
    &                               recon_lsq_cell_l_svd, recon_lsq_cell_q_svd, &
    &                               recon_lsq_cell_cpoor_svd,                   &
    &                               recon_lsq_cell_c_svd, div
  USE mo_interpol_config,     ONLY: llsq_lin_consv, lsq_high_ord, lsq_high_set
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_edge, rbf_interpol_c2grad                         
  USE mo_intp,                ONLY: cells2edges_scalar
  USE mo_dynamics_config,     ONLY: iequations
  USE mo_ha_dyn_config,       ONLY: ha_dyn_config
  USE mo_nonhydrostatic_config, ONLY: itime_scheme_nh_atm => itime_scheme
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ntracer
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c 
  USE mo_sync,                ONLY: SYNC_C, SYNC_C1, sync_patch_array,          &
    &                               sync_patch_array_mult, sync_patch_array_4de3
  USE mo_parallel_config,     ONLY: p_test_run
  USE mo_advection_config,    ONLY: advection_config, lcompute, lcleanup
  USE mo_advection_utils,     ONLY: laxfr_upflux
  USE mo_advection_quadrature,ONLY: prep_gauss_quadrature_q,                    &
    &                               prep_gauss_quadrature_cpoor,                &
    &                               prep_gauss_quadrature_c
  USE mo_advection_traj,      ONLY: btraj, btraj_o2, btraj_dreg,                &
    &                               btraj_dreg_nosort, divide_flux_area 
  USE mo_advection_limiter,   ONLY: hflx_limiter_mo, hflx_limiter_sm,           &
    &                               h_miura_slimiter_mo, h_miura_slimiter_sm,   &
    &                               shift_gauss_points
  USE mo_df_test,             ONLY: df_distv_barycenter !, df_cell_indices
  USE mo_ha_testcases,        ONLY: ctest_name

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  PUBLIC :: hor_upwind_flux
  PUBLIC :: upwind_hflux_up
  PUBLIC :: upwind_hflux_miura
  PUBLIC :: upwind_hflux_miura3


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
  !!
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-02-10)
  !! Modification by Daniel Reinert (2010-11-05)
  !! - tracer loop moved from step_advection to hor_upwind_flux
  !! Modification by Daniel Reinert (2010-11-09)
  !! - removed MUSCL-type horizontal advection
  !!
  !
  ! !LITERATURE
  ! MUSCL: Ahmad et al. (2006), Int. J. Num. Meth. Fluids, 50, 1247-1268
  !        Lin et al. (1994), MWR, 122, 1575-1593
  ! MIURA: Miura, H. (2007), Mon. Wea. Rev., 135, 4038-4044
  !
  SUBROUTINE hor_upwind_flux( p_cc, p_c0, p_rho, p_mass_flx_e, p_vn, p_dtime,     &
    &                     p_patch, p_int, p_ihadv_tracer, p_igrad_c_miura,        &
    &                     p_itype_hlimit, p_iadv_slev, p_iord_backtraj, p_upflux, &
    &                     opt_rlend )


    TYPE(t_patch), TARGET, INTENT(IN) ::  &     !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  & !< pointer to data structure for interpolation
      &  p_int

    REAL(wp),TARGET, INTENT(IN) ::  & !< advected cell centered variable
      &  p_cc(:,:,:,:)              !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp),TARGET, INTENT(IN) ::  & !< advected cell centered variable step (n)
      &  p_c0(:,:,:,:)              !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp),TARGET, INTENT(IN) ::  & !< density at cell center step (n)
      &  p_rho(:,:,:)               !< dim: (nproma,nlev,nblks_c)

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

    INTEGER :: jt                   !< tracer loop index
    INTEGER :: jg                   !< patch ID
    INTEGER :: i_rlend, i_rlend_vt  
    INTEGER :: qvsubstep_elev       !< end level for qv-substepping

    REAL(wp)::   &                  !< unweighted tangential velocity
      &  z_real_vt(nproma,p_patch%nlev,p_patch%nblks_e)!< component at edges

    !-----------------------------------------------------------------------

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 1
    ENDIF

    ! In case that different transport schemes (MIURA, MIURA3) are used 
    ! for different tracers, the double computation of tangential velocity 
    ! vt should be avoided. Instead of computing vt inside each of the 
    ! flux-routines, vt is computed only once per timestep prior to the flux 
    ! routines. The resulting tangential velocity field is then passed to  
    ! MIURA and MIURA3 as optional argument. 

    IF (ANY(p_ihadv_tracer(:)/= UP) .AND. ANY(p_ihadv_tracer(:)/= UP3)) THEN 

      i_rlend_vt = MIN(i_rlend, min_rledge_int - 1)

      IF (ANY(p_itype_hlimit(:)== islopel_sm) .OR. ANY(p_itype_hlimit(:)== islopel_m)) THEN
        i_rlend_vt = MIN(i_rlend, min_rledge_int - 2)
      ENDIF

      IF ( p_iord_backtraj /= 1 ) THEN
        i_rlend_vt = min_rledge_int - 3
      ENDIF

      ! reconstruct tangential velocity component at edge midpoints
      CALL rbf_vec_interpol_edge( p_vn, p_patch, p_int,            &! in
        &                         z_real_vt, opt_rlend=i_rlend_vt  )! inout
    ENDIF




    DO jt = 1, ntracer ! Tracer loop

      ! Select desired flux calculation method
      SELECT CASE( p_ihadv_tracer(jt) )

      CASE( UP )
        ! CALL first order upwind
        CALL upwind_hflux_up( p_patch, p_cc(:,:,:,jt),                &! in
          &                 p_mass_flx_e, p_upflux(:,:,:,jt),         &! in,inout  
          &                 opt_slev=p_iadv_slev(jt),opt_rlend=i_rlend)! in


      CASE( MIURA )
        ! CALL MIURA with second order accurate reconstruction
        CALL upwind_hflux_miura( p_patch, p_cc(:,:,:,jt), p_mass_flx_e,    &! in
          &                 p_vn, p_dtime, p_int, lcompute%miura_h(jt),    &! in
          &                 lcleanup%miura_h(jt), p_igrad_c_miura,         &! in
          &                 p_itype_hlimit(jt), p_iord_backtraj,           &! in
          &                 p_upflux(:,:,:,jt), opt_lconsv= llsq_lin_consv,&! inout,in
          &                 opt_real_vt=z_real_vt,                         &! in
          &                 opt_slev=p_iadv_slev(jt), opt_rlend=i_rlend    )! in


      CASE( MIURA3 )
        ! CALL MIURA with third order accurate reconstruction
        CALL upwind_hflux_miura3( p_patch, p_cc(:,:,:,jt), p_mass_flx_e, &! in
          &                 p_vn, p_dtime, p_int, lcompute%miura3_h(jt), &! in
          &                 lcleanup%miura3_h(jt), p_itype_hlimit(jt),   &! in
          &                 p_upflux(:,:,:,jt), opt_real_vt=z_real_vt,   &! inout,in
          &                 opt_slev=p_iadv_slev(jt), opt_rlend=i_rlend  )! in

      CASE( FFSL )
        ! CALL Flux form semi lagrangian scheme (extension of MIURA3-scheme) 
        ! with second or third order accurate reconstruction
        CALL upwind_hflux_ffsl( p_patch, p_cc(:,:,:,jt), p_mass_flx_e,    &! in
          &                 p_vn, p_dtime, p_int, lcompute%ffsl_h(jt),    &! in
          &                 lcleanup%ffsl_h(jt), p_itype_hlimit(jt),      &! in
          &                 p_upflux(:,:,:,jt), opt_real_vt=z_real_vt,    &! inout,in
          &                 opt_slev=p_iadv_slev(jt), opt_rlend=i_rlend   )! in


      CASE( UP3 )
        ! CALL 3rd order upwind (only for hexagons, currently)
        CALL upwind_hflux_hex( p_patch, p_int, p_cc(:,:,:,jt), p_c0(:,:,:,jt), &! in
          &                 p_mass_flx_e, p_dtime, p_itype_hlimit(jt),         &! in
          &                 p_upflux(:,:,:,jt), opt_slev=p_iadv_slev(jt)       )! inout


      CASE ( MCYCL )
        ! CALL MIURA with second order accurate reconstruction and subcycling
        CALL upwind_hflux_miura_cycl( p_patch, p_cc(:,:,:,jt), p_rho,     &! in
          &             p_mass_flx_e, p_vn, p_dtime, 2, p_int,            &! in
          &             lcompute%mcycl_h(jt), lcleanup%mcycl_h(jt),       &! in
          &             p_igrad_c_miura, p_itype_hlimit(jt),              &! in
          &             p_iord_backtraj, p_upflux(:,:,:,jt),              &! in,inout
          &             opt_lconsv= llsq_lin_consv, opt_real_vt=z_real_vt,&! in 
          &             opt_slev=p_iadv_slev(jt), opt_rlend=i_rlend       )! in


      CASE( MIURA_MCYCL )

        ! get patch ID
        jg = p_patch%id
        qvsubstep_elev = advection_config(jg)%iadv_qvsubstep_elev

        ! CALL standard MIURA for lower atmosphere and the subcycling version of 
        ! MIURA for upper atmosphere
        CALL upwind_hflux_miura( p_patch, p_cc(:,:,:,jt), p_mass_flx_e,  &! in
          &            p_vn, p_dtime, p_int, lcompute%miura_h(jt),       &! in
          &            lcleanup%miura_h(jt), p_igrad_c_miura,            &! in
          &            p_itype_hlimit(jt), p_iord_backtraj,              &! in
          &            p_upflux(:,:,:,jt), opt_lconsv= llsq_lin_consv,   &! inout,in
          &            opt_real_vt=z_real_vt, opt_slev=qvsubstep_elev+1, &! in
          &            opt_elev=p_patch%nlev, opt_rlend=i_rlend          )! in

        ! Note that lcompute/lcleanup%miura_mcycl_h is only used for miura 
        ! with substepping. This prevents us from computing the backward 
        ! trajectories twice for the standard miura3-scheme.
        CALL upwind_hflux_miura_cycl( p_patch, p_cc(:,:,:,jt), p_rho,    &! in
          &              p_mass_flx_e, p_vn, p_dtime, 2, p_int,          &! in
          &              lcompute%miura_mcycl_h(jt),                     &! in
          &              lcleanup%miura_mcycl_h(jt),                     &! in
          &              p_igrad_c_miura, p_itype_hlimit(jt),            &! in
          &              p_iord_backtraj, p_upflux(:,:,:,jt),            &! in,inout
          &              opt_lconsv= llsq_lin_consv,                     &! in
          &              opt_real_vt=z_real_vt, opt_slev=p_iadv_slev(jt),&! in
          &              opt_elev=qvsubstep_elev, opt_rlend=i_rlend      )! in


      CASE( MIURA3_MCYCL )

        ! get patch ID
        jg = p_patch%id
        qvsubstep_elev = advection_config(jg)%iadv_qvsubstep_elev

        ! CALL standard MIURA3 for lower atmosphere and the subcycling version of 
        ! MIURA for upper atmosphere
        CALL upwind_hflux_miura3( p_patch, p_cc(:,:,:,jt), p_mass_flx_e, &! in
          &           p_vn, p_dtime, p_int, lcompute%miura3_h(jt),       &! in
          &           lcleanup%miura3_h(jt), p_itype_hlimit(jt),         &! in
          &           p_upflux(:,:,:,jt), opt_real_vt=z_real_vt,         &! inout,in
          &           opt_slev=qvsubstep_elev+1, opt_elev=p_patch%nlev,  &! in
          &           opt_rlend=i_rlend                                  )! in

        IF (qvsubstep_elev > 0) &
        ! Note that lcompute/lcleanup%miura3_mcycl_h is only used for miura 
        ! with substepping. This prevents us from computing the backward 
        ! trajectories twice for the standard miura3-scheme.
        CALL upwind_hflux_miura_cycl( p_patch, p_cc(:,:,:,jt), p_rho,    &! in
          &              p_mass_flx_e, p_vn, p_dtime, 2, p_int,          &! in
          &              lcompute%miura3_mcycl_h(jt),                    &! in
          &              lcleanup%miura3_mcycl_h(jt),                    &! in
          &              p_igrad_c_miura, p_itype_hlimit(jt),            &! in
          &              p_iord_backtraj, p_upflux(:,:,:,jt),            &! in,inout
          &              opt_lconsv= llsq_lin_consv,                     &! in
          &              opt_real_vt=z_real_vt, opt_slev=p_iadv_slev(jt),&! in
          &              opt_elev=qvsubstep_elev, opt_rlend=i_rlend      )! in
      END SELECT

    END DO  ! Tracer loop

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
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,   &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)


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
          p_upflux(je,jk,jb) =  &
            &  laxfr_upflux( p_mass_flx_e(je,jk,jb), p_cc(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            &                             p_cc(iilc(je,jb,2),jk,iibc(je,jb,2)) )

        END DO  ! end loop over edges

      END DO  ! end loop over levels

    END DO  ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

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
  !!
  !! @par LITERATURE
  !! - Miura, H. (2007), Mon. Weather Rev., 135, 4038-4044
  !! - Skamarock, W.C. (2010), Conservative Transport schemes for Spherical Geodesic
  !!   Grids: High-order Reconstructions for Forward-in-Time Schemes, Mon. Wea. Rev,
  !!   in Press
  !!
  SUBROUTINE upwind_hflux_miura( p_patch, p_cc, p_mass_flx_e, p_vn, p_dtime,  &
    &                   p_int, ld_compute, ld_cleanup, p_igrad_c_miura,       &
    &                   p_itype_hlimit, p_iord_backtraj, p_out_e, opt_lconsv, &
    &                   opt_rlstart, opt_rlend, opt_lout_edge, opt_real_vt,   &
    &                   opt_slev, opt_elev )


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_hflux: upwind_hflux_miura'

    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), TARGET, INTENT(IN) ::     &   !< cell centered variable to be advected
      &  p_cc(:,:,:)                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::    &    !< unweighted velocity field
      &  p_vn(:,:,:)                !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    LOGICAL, INTENT(IN) ::     &    !< flag, if .TRUE. compute geometrical terms
      &  ld_compute

    LOGICAL, INTENT(IN) ::     &    !< flag, if .TRUE. clean up geometrical terms
      &  ld_cleanup

    INTEGER, INTENT(IN) :: p_igrad_c_miura   !< parameter to select the gradient
                                             !< reconstruction method at cell center

    INTEGER, INTENT(IN) :: p_itype_hlimit    !< parameter to select the limiter
                                             !< for horizontal transport

    INTEGER, INTENT(IN) :: p_iord_backtraj   !< parameter to select the spacial order
                                             !< of accuracy for the backward trajectory

    REAL(wp), INTENT(INOUT) ::  &   !< output field, containing the upwind flux or the
      &  p_out_e(:,:,:)             !< reconstructed edge value; dim: (nproma,nlev,nblks_e)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: if true, conservative reconstruction
     &  opt_lconsv                     !< is used

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: output edge value (.TRUE.),
     & opt_lout_edge                   !< or the flux across the edge (.FALSE./not specified)

    REAL(wp), INTENT(IN), OPTIONAL,  & ! optional: tangential velocity
     & TARGET :: opt_real_vt(:,:,:)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    LOGICAL  :: l_out_edgeval          !< corresponding local variable; default .FALSE.
                                       !< i.e. output flux across the edge

    REAL(wp), TARGET ::    &                   !< reconstructed gradient vector at
      &  z_grad(nproma,p_patch%nlev,p_patch%nblks_c,2) 
                                               !< cell center (geographical coordinates)

    REAL(wp), TARGET ::    &                        !< coefficient of the lsq reconstruction
      &  z_lsq_coeff(nproma,p_patch%nlev,p_patch%nblks_c,3) 
                                                    !< at cell center (geogr. coordinates)
                                                    !< includes coeff0 and gradients in
                                                    !< zonal and meridional direction

    REAL(wp), TARGET ::   &                    !< unweighted tangential velocity
      &  z_real_vt(nproma,p_patch%nlev,p_patch%nblks_e)!< component at edges

    REAL(wp), POINTER :: ptr_real_vt(:,:,:)    !< pointer to z_real_vt or opt_real_vt

    REAL(wp), POINTER :: ptr_cc(:,:,:)         !< ptr to tracer field or coefficient c0
                                               !< of lsq reconstruction (for nonconservative
                                               !< lsq ptr_cc again equals the tracer field).
    REAL(wp), POINTER :: ptr_grad(:,:,:,:)     !< Pointer to reconstructed zonal and
                                               !< meridional gradients.

    REAL(wp), ALLOCATABLE, SAVE ::  &   !< distance vectors cell center -->
      &  z_distv_gausspoint(:,:,:,:,:)  !< shifted Gauss-points
                                        !< dim: (nproma,nlev,p_patch%nblks_c,3,2)

    REAL(wp), ALLOCATABLE, SAVE ::  &   !< distance vectors cell center -->
      &  z_distv_bary(:,:,:,:)          !< barycenter of advected area
                                        !< (geographical coordinates)
                                        !< dim: (nproma,nlev,p_patch%nblks_e,2)

    INTEGER, ALLOCATABLE, SAVE  ::  &   !< line and block indices of cell centers
      &  z_cell_indices(:,:,:,:)        !< in which the calculated barycenters are
                                        !< located
                                        !< dim: (nproma,nlev,p_patch%nblks_e,2)

    REAL(wp) :: z_dthalf                !< \Delta t/2

    INTEGER  :: pid
    INTEGER  :: nlev               !< number of full levels
    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: ist                !< status variable
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    INTEGER  :: ilc0, ibc0         !< line and block index for local cell center
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom, i_rlend_c, i_rlend_tr, i_rlend_vt
    LOGICAL  :: l_consv            !< true if conservative lsq reconstruction is used
    INTEGER  :: itime_scheme

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

    IF ( PRESENT(opt_real_vt) ) THEN
      ptr_real_vt => opt_real_vt
    ELSE
      ptr_real_vt => z_real_vt
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

    IF (p_igrad_c_miura == 3) THEN
      i_rlend_c = min_rlcell_int
    ELSE
      i_rlend_c = min_rlcell_int - 1
    ENDIF

    IF (p_itype_hlimit == islopel_sm .OR. p_itype_hlimit == islopel_m) THEN
      i_rlend_tr = MIN(i_rlend, min_rledge_int - 2)
    ELSE
      i_rlend_tr = MIN(i_rlend, min_rledge_int - 1)
    ENDIF

    IF (p_iord_backtraj == 1)  THEN
      i_rlend_vt = i_rlend_tr 
    ELSE
      i_rlend_vt = MAX(i_rlend_tr - 1, min_rledge)
    ENDIF


    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    IF (p_test_run) THEN
      z_grad(:,:,:,:) = 0._wp
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
    ! 3 options:  without limiter
    !             with    slope limiter (modified Barth-Jesperson limiter)
    !             with    flux limiter following Zalesak (1979)
    !

    IF ( ld_compute ) THEN
      ! allocate temporary arrays for distance vectors and upwind cells
      ALLOCATE( z_distv_bary(nproma,nlev,p_patch%nblks_e,2),             &
        &       z_cell_indices(nproma,nlev,p_patch%nblks_e,2),           &
        &       STAT=ist )
      IF (p_itype_hlimit == islopel_m .OR. p_itype_hlimit == islopel_sm) &
        ALLOCATE( z_distv_gausspoint(nproma,nlev,p_patch%nblks_c,3,2),   &
          &       STAT=ist )

      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                     &
          &  'allocation for z_distv_bary, z_cell_indices, ' //          &
          &  'z_distv_gausspoint failed' )
      ENDIF
    END IF

    !
    ! 1. reconstruction of (unlimited) cell based gradient (lat,lon)
    !
    IF (p_igrad_c_miura == 1) THEN
      ! least squares method
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_l_svd( p_cc, p_patch, p_int%lsq_lin, z_lsq_coeff,     &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_lconsv=l_consv)
      ELSE
      CALL recon_lsq_cell_l( p_cc, p_patch, p_int%lsq_lin, z_lsq_coeff,         &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_lconsv=l_consv)
      ENDIF

      IF (l_consv) THEN
        ptr_cc   => z_lsq_coeff(:,:,:,1)   ! first coefficient of lsq rec.
      ELSE
        ptr_cc   => p_cc(:,:,:)
      ENDIF
      ptr_grad => z_lsq_coeff(:,:,:,2:3) ! gradients of lsq rec.

    ELSE IF (p_igrad_c_miura == 2) THEN
      ! Green-Gauss method
      CALL grad_green_gauss_cell( p_cc, p_patch, p_int, z_grad, opt_slev=slev, &
        &                         opt_elev=elev, opt_rlend=i_rlend_c )

      ptr_cc   => p_cc(:,:,:)
      ptr_grad => z_grad(:,:,:,:)

    ELSE IF (p_igrad_c_miura == 3) THEN
      ! RBF-based method (for rlstart=2 we run into a sync-error with nests)
      CALL rbf_interpol_c2grad( p_cc, p_patch, p_int,                             &
        &                       z_grad(:,:,:,1), z_grad(:,:,:,2),                 &
        &                       opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c,&
        &                       opt_rlstart=4 )

      ptr_cc   => p_cc(:,:,:)
      ptr_grad => z_grad(:,:,:,:)

    ENDIF


    !
    ! 2. Approximation of the 'departure region'. In case of a linear
    !    reconstruction it is sufficient to calculate the barycenter
    !    of the departure region (instead of all the vertices).
    !

    ! one half of current time step
    z_dthalf = 0.5_wp * p_dtime


    IF (ld_compute) THEN

      IF (.NOT. PRESENT(opt_real_vt)) THEN
        ! reconstruct tangential velocity component at edge midpoints
        CALL rbf_vec_interpol_edge( p_vn, p_patch, p_int,             &! in
          &                         ptr_real_vt, opt_rlend=i_rlend_vt )! inout
      ENDIF

      ! compute barycenter of departure region and the distance vector between
      ! the cell center and the barycenter using backward trajectories.
      IF (p_iord_backtraj == 1)  THEN

        ! first order backward trajectory
        CALL btraj   ( p_patch, p_int, p_vn, ptr_real_vt, z_dthalf,   &! in
          &            z_cell_indices, z_distv_bary,                  &! out
          &            opt_rlstart=i_rlstart, opt_rlend=i_rlend_tr,   &! in
          &            opt_slev=slev, opt_elev=elev                   )! in

      ELSE

        ! second order backward trajectory
        CALL btraj_o2( p_patch, p_int, p_vn, ptr_real_vt, z_dthalf,   &! in
          &            z_cell_indices, z_distv_bary,                  &! out
          &            opt_rlstart=i_rlstart, opt_rlend=i_rlend_tr,   &! in
          &            opt_slev=slev, opt_elev=elev                   )! in

      ENDIF


      !
      ! This section has been included for testing purposes
      !
      SELECT CASE (iequations)
      CASE (inh_atmosphere)
        itime_scheme = itime_scheme_nh_atm

      CASE (ihs_atm_temp, ihs_atm_theta, ishallow_water)
        itime_scheme = ha_dyn_config%itime_scheme

      CASE DEFAULT
        CALL finish(TRIM(routine),'cannot get the value of itime_scheme')
      END SELECT

      SELECT CASE (itime_scheme)
       !------------------
       ! Pure advection
       !------------------
       CASE (TRACER_ONLY)

         SELECT CASE ( TRIM(ctest_name) )
           CASE ('DF1', 'DF2', 'DF3', 'DF4') ! deformational flow

!DR           z_cell_indices(:,:,:,:) = df_cell_indices(:,:,:,:)
!DR           z_distv_bary(:,:,:,:)   = df_distv_barycenter(:,:,:,:)

            ! compute difference between distance vector based on 3rd order
            ! Taylor series approximation and operational distance vector of either
            ! first or second order accurracy.
            df_distv_barycenter(:,:,:,:) = df_distv_barycenter(:,:,:,:)      &
              &                          - z_distv_bary(:,:,:,:)

            ! now compute length of distance vectors and store them in the first
            ! position of the df_distv_barycenter field.
            df_distv_barycenter(:,:,:,1) = SQRT(df_distv_barycenter(:,:,:,1) &
              &                          * df_distv_barycenter(:,:,:,1)      &
              &                          + df_distv_barycenter(:,:,:,2)      &
              &                          * df_distv_barycenter(:,:,:,2))

         END SELECT
      END SELECT

    END IF ! ld_compute


    ! 3. If desired, slopes of reconstruction are limited using either a monotonous
    !    or a semi-monotonous version of the Barth-Jespersen slope limiter.
    !
    IF (p_itype_hlimit == islopel_sm) THEN

      IF (ld_compute) THEN
        ! calculate shifted gauss points and distance vector between shifted
        ! gauss point and cell-circumcenter
        CALL shift_gauss_points( p_patch, p_int,                                  & !in
          &                      z_dthalf, ptr_real_vt,                           & !in
          &                      z_distv_gausspoint,                              & !out
          &                      opt_rlend_c=i_rlend_c, opt_rlstart_c=4,          & !in
          &                      opt_rlend_e=i_rlend_vt, opt_rlstart_e=i_rlstart, & !in
          &                      opt_slev=slev, opt_elev=elev                     ) !in
      ENDIF

      ! semi-monotonic (sm) slope limiter
      CALL h_miura_slimiter_sm( p_patch, p_cc, ptr_cc, z_distv_gausspoint, & !in
        &                    ptr_grad, opt_rlend=i_rlend_c, opt_rlstart=4, & !inout,in
        &                    opt_slev=slev, opt_elev=elev                  ) !in

    ELSE IF (p_itype_hlimit == islopel_m) THEN

      IF (ld_compute) THEN
        ! calculate shifted gauss points and distance vector between shifted
        ! gauss point and cell-circumcenter
        CALL shift_gauss_points( p_patch, p_int,                                  & !in
          &                      z_dthalf, ptr_real_vt,                           & !in
          &                      z_distv_gausspoint,                              & !out
          &                      opt_rlend_c=i_rlend_c, opt_rlstart_c=4,          & !in
          &                      opt_rlend_e=i_rlend_vt, opt_rlstart_e=i_rlstart, & !in
          &                      opt_slev=slev, opt_elev=elev                     ) !in
      ENDIF

      ! monotonous (mo) slope limiter
      CALL h_miura_slimiter_mo( p_patch, p_cc, ptr_cc, z_distv_gausspoint, & !in
        &                    ptr_grad, opt_rlend=i_rlend_c, opt_rlstart=4, & !inout,in
        &                    opt_slev=slev, opt_elev=elev                  ) !in

    ENDIF

    ! Synchronize gradient if 10-point RBF-based reconstruction is used.
    ! Only level 1 halo points are synchronized. Thus, independent of 
    ! the reconstruction method, the correct polynomial coefficients 
    ! are available up to halo points of level 1 (after this sync).
    IF (p_igrad_c_miura == 3) &
      & CALL sync_patch_array_mult(SYNC_C1,p_patch,2,f4din=ptr_grad)


    !
    ! 4. Calculate reconstructed tracer value at each barycenter
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

!$OMP WORKSHARE
      p_out_e(:,:,i_startblk:i_endblk) = 0._wp
!$OMP END WORKSHARE
    ENDIF

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    ! initialize also nest boundary points with zero
    IF ( l_out_edgeval .AND. (p_patch%id > 1 .OR. l_limited_area)) THEN
!$OMP WORKSHARE
      p_out_e(:,:,1:i_startblk) = 0._wp
!$OMP END WORKSHARE
    ENDIF

    IF (l_consv .AND. p_igrad_c_miura == 1 .AND.                           &
     & (p_itype_hlimit == islopel_sm .OR. p_itype_hlimit == islopel_m)) THEN

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,ilc0,ibc0), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

        IF ( l_out_edgeval ) THEN   ! Calculate 'edge value' of advected quantity

!CDIR UNROLL=5
          DO jk = slev, elev

            DO je = i_startidx, i_endidx

              ! Calculate reconstructed tracer value at barycenter of rhomboidal
              ! area which is swept across the corresponding edge.
              ilc0 = z_cell_indices(je,jk,jb,1)
              ibc0 = z_cell_indices(je,jk,jb,2)

              ! Calculate 'edge value' of advected quantity (cc_bary)
              p_out_e(je,jk,jb) = p_cc(ilc0,jk,ibc0)                                       &
                &  + ( (z_distv_bary(je,jk,jb,1) - p_int%lsq_lin%lsq_moments(ilc0,ibc0,1)) &
                &  * ptr_grad(ilc0,jk,ibc0,1)                                              &
                &  +   (z_distv_bary(je,jk,jb,2) - p_int%lsq_lin%lsq_moments(ilc0,ibc0,2)) &
                &  * ptr_grad(ilc0,jk,ibc0,2) )

            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels

        ELSE

!CDIR UNROLL=5
          DO jk = slev, elev
            DO je = i_startidx, i_endidx

              ! Calculate reconstructed tracer value at barycenter of rhomboidal
              ! area which is swept across the corresponding edge.
              ilc0 = z_cell_indices(je,jk,jb,1)
              ibc0 = z_cell_indices(je,jk,jb,2)

              ! Calculate flux at cell edge (cc_bary*v_{n}* \Delta p)
              p_out_e(je,jk,jb) = ( p_cc(ilc0,jk,ibc0)                                     &
                &  + ( (z_distv_bary(je,jk,jb,1) - p_int%lsq_lin%lsq_moments(ilc0,ibc0,1)) &
                &  * ptr_grad(ilc0,jk,ibc0,1)                                              &
                &  +   (z_distv_bary(je,jk,jb,2) - p_int%lsq_lin%lsq_moments(ilc0,ibc0,2)) &
                &  * ptr_grad(ilc0,jk,ibc0,2) )) * p_mass_flx_e(je,jk,jb)

            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels

        ENDIF

      END DO    ! loop over blocks
!$OMP END DO
 
    ELSE

      ! If no slope limiter is applied, a more efficient formulation without the
      ! moments can be used.

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,ilc0,ibc0), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

        IF ( l_out_edgeval ) THEN   ! Calculate 'edge value' of advected quantity

!CDIR UNROLL=5
          DO jk = slev, elev
            DO je = i_startidx, i_endidx

              ! Calculate reconstructed tracer value at barycenter of rhomboidal
              ! area which is swept across the corresponding edge.
              ilc0 = z_cell_indices(je,jk,jb,1)
              ibc0 = z_cell_indices(je,jk,jb,2)  

              ! Calculate 'edge value' of advected quantity (cc_bary)
              p_out_e(je,jk,jb) = ptr_cc(ilc0,jk,ibc0)                       &
                &    + z_distv_bary(je,jk,jb,1) * ptr_grad(ilc0,jk,ibc0,1)   &
                &    + z_distv_bary(je,jk,jb,2) * ptr_grad(ilc0,jk,ibc0,2)

            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels

        ELSE

!CDIR UNROLL=5
          DO jk = slev, elev
            DO je = i_startidx, i_endidx

              ! Calculate reconstructed tracer value at barycenter of rhomboidal
              ! area which is swept across the corresponding edge.  
              ilc0 = z_cell_indices(je,jk,jb,1)
              ibc0 = z_cell_indices(je,jk,jb,2)

              ! Calculate flux at cell edge (cc_bary*v_{n}* \Delta p)
              p_out_e(je,jk,jb) = ( ptr_cc(ilc0,jk,ibc0)                        &
                  &    + z_distv_bary(je,jk,jb,1) * ptr_grad(ilc0,jk,ibc0,1)    &
                  &    + z_distv_bary(je,jk,jb,2) * ptr_grad(ilc0,jk,ibc0,2) )  &
                  &    * p_mass_flx_e(je,jk,jb)

            ENDDO ! loop over edges
          ENDDO   ! loop over vertical levels

        ENDIF

      ENDDO    ! loop over blocks
!$OMP END DO NOWAIT

    ENDIF
!$OMP END PARALLEL


    !
    ! 5. If desired, apply a (semi-)monotone flux limiter to limit computed fluxes.
    !    The flux limiter is based on work by Zalesak (1979)
    IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_m) THEN
      CALL hflx_limiter_mo( p_patch, p_int, p_dtime, p_cc, p_mass_flx_e, & !in
        &                p_out_e, opt_rlend=i_rlend, opt_slev=slev,      & !inout,in
        &                opt_elev=elev                                   ) !in

    ELSE IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_sm) THEN
      ! MPI-sync necessary (to be precise: only necessary for
      ! igrad_c_miura /= 3. For simplicity, we perform the sync for
      ! igrad_c_miura = 3 as well.
      CALL hflx_limiter_sm( p_patch, p_int, p_dtime, p_cc, p_out_e, & !in,inout
        &                   opt_rlend=i_rlend, opt_slev=slev,       & !in
        &                   opt_elev=elev                           ) !in
    ENDIF


    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays for velocity, Gauss-points and barycenters
      DEALLOCATE( z_distv_bary, z_cell_indices, STAT=ist )
      IF (p_itype_hlimit == islopel_m .OR. p_itype_hlimit == islopel_sm) &
        & DEALLOCATE( z_distv_gausspoint, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                &
          &  'deallocation for z_distv_bary, z_cell_indices, '  //  &
          &  'z_distv_gausspoint failed' )
      ENDIF
    END IF

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
  !!
  !! @par LITERATURE
  !! - Miura, H. (2007), Mon. Weather Rev., 135, 4038-4044
  !! - Skamarock, W.C. (2010), Conservative Transport schemes for Spherical Geodesic
  !!   Grids: High-order Reconstructions for Forward-in-Time Schemes, Mon. Wea. Rev,
  !!   in Press
  !!
  SUBROUTINE upwind_hflux_miura_cycl( p_patch, p_cc, p_rho, p_mass_flx_e, p_vn,    &
    &                   p_dtime,  p_ncycl, p_int, ld_compute, ld_cleanup,          &
    &                   p_igrad_c_miura, p_itype_hlimit, p_iord_backtraj, p_out_e, &
    &                   opt_lconsv, opt_rlstart, opt_rlend, opt_real_vt,           &
    &                   opt_slev, opt_elev  )


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_hflux: upwind_hflux_miura_cycl'

    TYPE(t_patch), TARGET, INTENT(IN) ::  &   !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), TARGET, INTENT(IN) ::     &   !< cell centered variable to be advected
      &  p_cc(:,:,:)                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< density (i.e. \rho\Delta z) at timestep n
      &  p_rho(:,:,:)               !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)        !< Assumption: constant over p_dtime
                                    !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::    &    !< unweighted velocity field
      &  p_vn(:,:,:)                !< Assumption: constant over p_dtime
                                    !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    INTEGER,  INTENT(IN) ::    &    !< number of sub-timesteps into which p_dtime
      &  p_ncycl                    !< is splitted (p_ncycl=1 : no subcycling)

    LOGICAL, INTENT(IN) ::     &    !< flag, if .TRUE. compute geometrical terms
      &  ld_compute

    LOGICAL, INTENT(IN) ::     &    !< flag, if .TRUE. clean up geometrical terms
      &  ld_cleanup

    INTEGER, INTENT(IN) :: p_igrad_c_miura   !< parameter to select the gradient
                                             !< reconstruction method at cell center

    INTEGER, INTENT(IN) :: p_itype_hlimit    !< parameter to select the limiter
                                             !< for horizontal transport

    INTEGER, INTENT(IN) :: p_iord_backtraj   !< parameter to select the order of
                                             !< spacial accuracy for the backward trajectory

    REAL(wp), INTENT(INOUT) ::  &   !< output field, containing the upwind flux or the
      &  p_out_e(:,:,:)             !< reconstructed edge value; dim: (nproma,nlev,nblks_e)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: if true, conservative reconstruction
     &  opt_lconsv                     !< is used

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    REAL(wp), INTENT(IN), OPTIONAL,  & ! optional: tangential velocity
      & TARGET :: opt_real_vt(:,:,:)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    REAL(wp), TARGET ::    &                   !< reconstructed gradient vector at
      &  z_grad(nproma,p_patch%nlev,p_patch%nblks_c,2) 
                                               !< cell center (geographical coordinates)

    REAL(wp), TARGET ::    &                        !< coefficient of the lsq reconstruction
      &  z_lsq_coeff(nproma,p_patch%nlev,p_patch%nblks_c,3) 
                                                    !< at cell center (geogr. coordinates)
                                                    !< includes coeff0 and gradients in
                                                    !< zonal and meridional direction

    REAL(wp), TARGET ::   &                    !< unweighted tangential velocity
      &  z_real_vt(nproma,p_patch%nlev,p_patch%nblks_e)!< component at edges

    REAL(wp), POINTER :: ptr_real_vt(:,:,:)    !< pointer to z_real_vt or opt_real_vt

    REAL(wp), POINTER :: ptr_cc(:,:,:)         !< ptr to tracer field or coefficient c0
                                               !< of lsq reconstruction (for nonconservative
                                               !< lsq ptr_cc again equals the tracer field).
    REAL(wp), POINTER :: ptr_grad(:,:,:,:)     !< Pointer to reconstructed zonal and
                                               !< meridional gradients.

    REAL(wp), ALLOCATABLE, SAVE ::  &   !< distance vectors cell center -->
      &  z_distv_bary(:,:,:,:)          !< barycenter of advected area
                                        !< (geographical coordinates)
                                        !< dim: (nproma,nlev,p_patch%nblks_e,2)

    INTEGER, ALLOCATABLE, SAVE  ::  &   !< line and block indices of cell centers
      &  z_cell_indices(:,:,:,:)        !< in which the calculated barycenters are
                                        !< located
                                        !< dim: (nproma,nlev,p_patch%nblks_e,2)

    REAL(wp) :: z_dtsub                 !< sub timestep p_dtime/p_ncycl
    REAL(wp) :: z_dthalf                !< z_dtsub/2
    REAL(wp) ::                     &   !< tracer flux at n + nsub/p_ncycl
      &  z_tracer_mflx(nproma,p_patch%nlev,p_patch%nblks_e,p_ncycl)

    REAL(wp) ::                     &   !< tracer mass flux divergence at cell center
      &  z_rhofluxdiv_c(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) ::                     &   !< mass flux divergence at cell center
      &  z_fluxdiv_c(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp), TARGET ::             &   !< 'tracer cell value' at interm.  
      &  z_tracer(nproma,p_patch%nlev,p_patch%nblks_c,2) !< old and new timestep

    REAL(wp), TARGET ::             &   !< density (i.e. \rho\Delta z) at interm.  
      &  z_rho(nproma,p_patch%nlev,p_patch%nblks_c,2) !< old and new timestep

    INTEGER  :: pid
    INTEGER  :: nlev               !< number of full levels
    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: ist                !< status variable
    INTEGER  :: jc, je, jk, jb     !< index of cell, edge, vert level, block
    INTEGER  :: ilc0, ibc0         !< line and block index for local cell center
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom, i_rlend_c, i_rlend_tr, i_rlend_vt
    LOGICAL  :: l_consv            !< true if conservative lsq reconstruction is used
    INTEGER  :: nsub               !< counter for sub-timesteps
    INTEGER  :: nnow, nnew, nsav   !< time indices

   !-------------------------------------------------------------------------

    IF (p_ncycl /= 2) &
    CALL finish(TRIM(routine),'current implementation of upwind_hflux_miura_cycl '//&
      &                       'requires 2 subcycling steps (p_ncycl=2)')

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

    IF ( PRESENT(opt_real_vt) ) THEN
      ptr_real_vt => opt_real_vt
    ELSE
      ptr_real_vt => z_real_vt
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

    i_rlend_tr = MIN(i_rlend, min_rledge_int - 1)

    IF (p_igrad_c_miura == 3) THEN
      i_rlend_c = min_rlcell_int
    ELSE
      i_rlend_c = min_rlcell_int - 1
    ENDIF

    IF (p_iord_backtraj == 1)  THEN
      i_rlend_vt = i_rlend_tr 
    ELSE
      i_rlend_vt = MAX(i_rlend_tr - 1, min_rledge)
    ENDIF


    ! get local sub-timestep
    z_dtsub = p_dtime/REAL(p_ncycl,wp)

    ! one half of current time step
    z_dthalf = 0.5_wp * z_dtsub


    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    IF (p_test_run) THEN
      z_grad(:,:,:,:)   = 0._wp
      z_tracer(:,:,:,:) = 0._wp
    ENDIF

    ! initialize 'now slice' of z_tracer and z_rho
    nnow = 1
    nnew = 2
!$OMP PARALLEL
!$OMP WORKSHARE
    z_tracer(:,slev:elev,:,nnow) = p_cc (:,slev:elev,:)
    z_rho   (:,slev:elev,:,nnow) = p_rho(:,slev:elev,:)
!$OMP END WORKSHARE
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
    IF (ld_compute) THEN

      ! allocate temporary arrays for distance vectors and upwind cells
      ALLOCATE( z_distv_bary(nproma,nlev,p_patch%nblks_e,2),             &
        &       z_cell_indices(nproma,nlev,p_patch%nblks_e,2),           &
        &       STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                     &
          &  'allocation for z_distv_bary, z_cell_indices, failed' )
      ENDIF


      !
      ! 1. Approximation of the 'departure region'. In case of a linear
      !    reconstruction it is sufficient to compute the barycenter
      !    of the departure region (instead of all the vertices).
      !

      IF (.NOT. PRESENT(opt_real_vt)) THEN
        ! reconstruct tangential velocity component at edge midpoints
        CALL rbf_vec_interpol_edge( p_vn, p_patch, p_int,             &! in
          &                         ptr_real_vt, opt_rlend=i_rlend_vt )! inout
      ENDIF

      ! compute barycenter of departure region and the distance vector between
      ! the cell center and the barycenter using backward trajectories.
      IF (p_iord_backtraj == 1)  THEN

        ! first order backward trajectory
        CALL btraj   ( p_patch, p_int, p_vn, ptr_real_vt, z_dthalf,   &! in
          &            z_cell_indices, z_distv_bary,                  &! out
          &            opt_rlstart=i_rlstart, opt_rlend=i_rlend_tr,   &! in
          &            opt_slev=slev, opt_elev=elev                   )! in

      ELSE

        ! second order backward trajectory
        CALL btraj_o2( p_patch, p_int, p_vn, ptr_real_vt, z_dthalf,   &! in
          &            z_cell_indices, z_distv_bary,                  &! out
          &            opt_rlstart=i_rlstart, opt_rlend=i_rlend_tr,   &! in
          &            opt_slev=slev, opt_elev=elev                   )! in

      ENDIF

    END IF ! ld_compute




    !
    ! Loop over sub-timesteps (subcycling)
    !
    DO nsub=1, p_ncycl

      !
      ! 2. reconstruction of (unlimited) cell based gradient (lat,lon)
      !
      IF (p_igrad_c_miura == 1) THEN
        ! least squares method
        IF (advection_config(pid)%llsq_svd) THEN
        CALL recon_lsq_cell_l_svd( z_tracer(:,:,:,nnow), p_patch, p_int%lsq_lin, &
          &                    z_lsq_coeff, opt_slev=slev, opt_elev=elev,        &
          &                    opt_rlend=i_rlend_c, opt_lconsv=l_consv)
        ELSE
        CALL recon_lsq_cell_l( z_tracer(:,:,:,nnow), p_patch, p_int%lsq_lin,     &
          &                    z_lsq_coeff, opt_slev=slev, opt_elev=elev,        &
          &                    opt_rlend=i_rlend_c, opt_lconsv=l_consv)
        ENDIF

        IF (l_consv) THEN
          ptr_cc   => z_lsq_coeff(:,:,:,1)   ! first coefficient of lsq rec.
        ELSE
          ptr_cc   => z_tracer(:,:,:,nnow)
        ENDIF
        ptr_grad => z_lsq_coeff(:,:,:,2:3) ! gradients of lsq rec.

      ELSE IF (p_igrad_c_miura == 2) THEN
        ! Green-Gauss method
        CALL grad_green_gauss_cell( z_tracer(:,:,:,nnow), p_patch, p_int, &
          &                         z_grad, opt_slev=slev, opt_elev=elev, &
          &                         opt_rlend=i_rlend_c )

        ptr_cc   => z_tracer(:,:,:,nnow)
        ptr_grad => z_grad(:,:,:,:)

      ELSE IF (p_igrad_c_miura == 3) THEN
        ! RBF-based method (for rlstart=2 we run into a sync-error with nests)
        CALL rbf_interpol_c2grad( z_tracer(:,:,:,nnow), p_patch, p_int,             &
          &                       z_grad(:,:,:,1), z_grad(:,:,:,2),                 &
          &                       opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c,&
          &                       opt_rlstart=4 )

        ptr_cc   => z_tracer(:,:,:,nnow)
        ptr_grad => z_grad(:,:,:,:)

      ENDIF




      ! Synchronize gradient if 10-point RBF-based reconstruction is used.
      ! Only level 1 halo points are synchronized. Thus, independent of 
      ! the reconstruction method, the correct polynomial coefficients 
      ! are available up to halo points of level 1 (after this sync).
      !
      IF ( p_igrad_c_miura == 3 ) &
        &  CALL sync_patch_array_mult(SYNC_C1,p_patch,2,f4din=ptr_grad)




      !
      ! 3. Calculate reconstructed tracer value at each barycenter
      !    \Phi_{bary}=\Phi_{circum} + DOT_PRODUCT(\Nabla\Psi,r).
      !    and compute intermediate update of q.
      !

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)


      i_startblk = p_patch%edges%start_blk(i_rlstart,1)
      i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)


    IF ( p_patch%id > 1 .OR. l_limited_area) THEN
!$OMP WORKSHARE
      z_tracer_mflx(:,:,1:i_startblk,nsub) = 0._wp
!$OMP END WORKSHARE
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
!CDIR UNROLL=5
        DO jk = slev, elev

          DO je = i_startidx, i_endidx

            ilc0 = z_cell_indices(je,jk,jb,1)
            ibc0 = z_cell_indices(je,jk,jb,2)  

            ! compute intermediate flux at cell edge (cc_bary*v_{n}* \Delta p)
            z_tracer_mflx(je,jk,jb,nsub) = ( ptr_cc(ilc0,jk,ibc0)            &
              &      + z_distv_bary(je,jk,jb,1) * ptr_grad(ilc0,jk,ibc0,1)   &
              &      + z_distv_bary(je,jk,jb,2) * ptr_grad(ilc0,jk,ibc0,2) ) &
              &      * p_mass_flx_e(je,jk,jb)

          ENDDO ! loop over edges
        ENDDO   ! loop over vertical levels

      ENDDO    ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL


      ! 4. limit intermediate tracer fluxes to achieve positive definiteness
      !    The flux limiter is based on work by Zalesak (1979)
      !
      IF ( p_itype_hlimit == ifluxl_sm .OR. p_itype_hlimit == ifluxl_m ) THEN
        ! MPI-sync necessary (to be precise: only necessary for
        ! igrad_c_miura /= 3. For simplicity, we perform the sync for
        ! igrad_c_miura = 3 as well.
        !
        CALL hflx_limiter_sm( p_patch, p_int, z_dtsub          , & !in
          &                   z_tracer(:,:,:,nnow)             , & !in
          &                   z_tracer_mflx(:,:,:,nsub)        , & !inout
          &                   opt_rho = z_rho(:,:,:,nnow)      , & !in 
          &                   opt_rlend=i_rlend, opt_slev=slev , & !in
          &                   opt_elev=elev                      ) !in
      ENDIF


      ! during the last iteration step, the following computations can be skipped
      IF ( nsub == p_ncycl ) EXIT


      ! compute mass flux and tracer mass flux divergence
      !
      ! This computation needs to be done only once, since the mass flux
      ! p_mass_flx_e is assumed to be constant in time.
      IF ( nsub == 1 ) THEN
        CALL div( p_mass_flx_e(:,:,:), p_patch, p_int,    &! in
          &       z_rhofluxdiv_c(:,:,:),                  &! inout
          &       opt_slev=slev, opt_elev=elev,           &! in
          &       opt_rlend=min_rlcell_int                )! in
      ENDIF


      CALL div( z_tracer_mflx(:,:,:,nsub), p_patch, p_int, &! in
        &       z_fluxdiv_c(:,:,:),                        &! inout
        &       opt_slev=slev, opt_elev=elev,              &! in
        &       opt_rlstart=3, opt_rlend=min_rlcell_int    )! in


      !
      ! 4.1/4.2 compute updated density and tracer fields
      !
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

      i_startblk = p_patch%cells%start_blk(3,1)
      i_endblk   = p_patch%cells%end_blk(min_rlcell_int,i_nchdom)


    ! initialize also nest boundary points with zero
    IF ( p_patch%id > 1 .OR. l_limited_area) THEN
!$OMP WORKSHARE
      z_tracer(:,:,1:i_startblk,nnew) = 0._wp
!$OMP END WORKSHARE
    ENDIF

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, 3, min_rlcell_int)

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
            z_tracer(jc,jk,jb,nnew) = ( z_tracer(jc,jk,jb,nnow)          &
              &                      * z_rho(jc,jk,jb,nnow)              &
              &                      - z_dtsub * z_fluxdiv_c(jc,jk,jb) ) &
              &                      / z_rho(jc,jk,jb,nnew)
          ENDDO
        ENDDO

      ENDDO    ! loop over blocks
!$OMP ENDDO
!$OMP END PARALLEL

      nsav = nnow
      nnow = nnew
      nnew = nsav


      CALL sync_patch_array(SYNC_C,p_patch,z_tracer(:,:,:,nnow))


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

!CDIR UNROLL=5
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

            ! Calculate flux at cell edge (cc_bary*v_{n}* \Delta p)
            p_out_e(je,jk,jb) = SUM(z_tracer_mflx(je,jk,jb,1:2))/REAL(p_ncycl,wp)

          ENDDO ! loop over edges
        ENDDO   ! loop over vertical levels

    ENDDO    ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL



    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays for velocity, Gauss-points and barycenters
      DEALLOCATE( z_distv_bary, z_cell_indices, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                 &
          &  'deallocation for z_distv_bary, z_cell_indices, failed' )
      ENDIF
    END IF

  END SUBROUTINE upwind_hflux_miura_cycl




  !-------------------------------------------------------------------------
  !>
  !! MIURA scheme with third order accurate lsq-reconstruction
  !!
  !! Calculation of time averaged horizontal tracer fluxes at triangle edges
  !! using a Flux form semi-lagrangian (FFSL)-method based on the MIURA-scheme.
  !! Unlike the standard Miura scheme, a third or fourth order accurate
  !! (i.e. quadratic, cubic) least squares reconstruction can be applied.
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
  SUBROUTINE upwind_hflux_miura3( p_patch, p_cc, p_mass_flx_e, p_vn, p_dtime,   &
    &                      p_int, ld_compute, ld_cleanup, p_itype_hlimit,       &
    &                      p_out_e, opt_rlstart, opt_rlend, opt_lout_edge,      &
    &                      opt_real_vt, opt_slev, opt_elev )


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_hflux: upwind_hflux_miura3'

    TYPE(t_patch), INTENT(IN)     ::  &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), INTENT(IN) ::    &    !< cell centered variable to be advected
      &  p_cc(:,:,:)                !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< contravariant horizontal mass flux at cell edge
      &  p_mass_flx_e(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::    &    !< normal component of velocity field
      &  p_vn(:,:,:)                !< dim: (nproma,nlev,nblks_e)

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

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: output edge value (.TRUE.),
     & opt_lout_edge                   !< or the flux across the edge (.FALSE./not specified)

    REAL(wp), INTENT(IN), OPTIONAL,  & !< optional: tangential velocity
     & TARGET :: opt_real_vt(:,:,:)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    LOGICAL  :: l_out_edgeval          !< corresponding local variable; default .FALSE.
                                       !< i.e. output flux across the edge

   REAL(wp) ::   &                     !< coefficients of lsq reconstruction
      &  z_lsq_coeff(nproma,p_patch%nlev,lsq_high_set%dim_unk+1,p_patch%nblks_c) 
                                       !< at cell center
                                       !< includes c0 and gradients in zonal and
                                       !< meridional direction

    REAL(wp), TARGET ::   &                    !< tangential component of velocity field
      &  z_real_vt(nproma,p_patch%nlev,p_patch%nblks_e)!< at edges

    REAL(wp), POINTER :: ptr_real_vt(:,:,:)    !< pointer to z_real_vt or opt_real_vt

    REAL(wp) ::  &                    !< coordinates of departure region vertices. The origin
      &  z_coords_dreg_v(nproma,4,2,p_patch%nlev,p_patch%nblks_e)
                                      !< of the coordinate system is at the circumcenter of
                                      !< the upwind cell. Unit vectors point to local East
                                      !< and North. (geographical coordinates)
                                      !< dim: (nproma,4,2,nlev,ptr_p%nblks_e)

    REAL(wp), ALLOCATABLE, SAVE ::  & !< gauss quadrature vector
      &  z_quad_vector_sum(:,:,:,:)   !< dim: (nproma,lsq_dim_unk+1,nlev,nblks_e)

    REAL(wp), ALLOCATABLE, SAVE ::  & !< area departure region [m**2]
      &  z_dreg_area(:,:,:)           !< dim: (nproma,nlev,nblks_e)

    INTEGER, ALLOCATABLE, SAVE, TARGET ::  & !< line and block indices of upwind cell
      &  z_cell_indices(:,:,:,:)             !< dim: (nproma,nlev,p_patch%nblks_e,2)

    INTEGER, POINTER, SAVE ::  &       !< Pointer to line and block indices of the cell
      &  ptr_ilc(:,:,:), ptr_ibc(:,:,:)!< center upstream of the edge

    INTEGER  :: nlev               !< number of full levels
    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: ist                !< status variable
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    INTEGER  :: dim_unk
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_rlend_c, i_nchdom
    INTEGER  :: pid                !< patch ID

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

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF ( PRESENT(opt_real_vt) ) THEN
      ptr_real_vt => opt_real_vt
    ELSE
      ptr_real_vt => z_real_vt
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

    IF ( ld_compute ) THEN
      ! allocate temporary arrays for quadrature and upwind cells
      ALLOCATE( z_quad_vector_sum(nproma,dim_unk,nlev,p_patch%nblks_e), &
        &       z_dreg_area(nproma,nlev,p_patch%nblks_e),               &
        &       z_cell_indices(nproma,nlev,p_patch%nblks_e,2),          &
        &       STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                 &
          &  'allocation for z_quad_vector_sum, z_dreg_area, ' //    &
          &  'z_cell_indices failed' )
      ENDIF
    END IF


    !
    ! 1. reconstruction of the tracer subgrid distribution
    !    least squares method
    !    Note: for rlstart=2 we run into a sync-error with nests
    !
    IF (lsq_high_ord == 2) THEN
      ! quadratic reconstruction
      ! (computation of 6 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_q_svd( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_q( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,        &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ELSE IF (lsq_high_ord == 30) THEN
      ! cubic reconstruction without cross derivatives
      ! (computation of 8 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_cpoor_svd( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,&
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_cpoor( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ELSE IF (lsq_high_ord == 3) THEN
      ! cubic reconstruction with cross derivatives
      ! (computation of 10 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_c_svd( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_c( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,        &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ENDIF


    ! Synchronize polynomial coefficients
    ! Note: a special sync routine is needed here because the fourth dimension
    ! of z_lsq_coeff is (for efficiency reasons) on the third index
    CALL sync_patch_array_4de3(SYNC_C1,p_patch,lsq_high_set%dim_unk+1,z_lsq_coeff)


    !
    ! 2. Approximation of the 'departure region'. The coordinates of
    !    all vertices are computed and stored in an edge-based data
    !    structure.
    !    In addition the Gauss-Legendre quadrature is prepared by
    !    calculating some tracer-invariant (i.e. purely geometric) fields.
    !
    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    IF (ld_compute) THEN

      IF (.NOT. PRESENT(opt_real_vt)) THEN
        ! reconstruct tangential velocity component at edge midpoints
        CALL rbf_vec_interpol_edge( p_vn, p_patch, p_int,          &! in
          &                         ptr_real_vt, opt_rlend=i_rlend )! inout
      ENDIF

      ! compute vertex coordinates for the departure region using a first
      ! order accurate (O(\Delta t)) backward trajectory-method
      CALL btraj_dreg( p_patch, p_int, p_vn, ptr_real_vt, p_dtime, &! in
        &              z_cell_indices, z_coords_dreg_v,            &! out
        &              opt_rlstart=i_rlstart, opt_rlend=i_rlend,   &! in
        &              opt_slev=slev, opt_elev=elev                )! in


!      ! In order to check, whether the vertices are stored in clockwise or
!      ! counterclockwise direction, calculate the cross product between the
!      ! vectors from point 2 to point 3 and point 2 to point 1
!      DO jb = i_startblk, i_endblk

!        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
!                           i_startidx, i_endidx, i_rlstart)
!        DO jk = 1, nlev

!          DO je = i_startidx, i_endidx
!            v23_x = z_coords_dreg_v(je,3,1,jk,jb) - z_coords_dreg_v(je,2,1,jk,jb)
!            v23_y = z_coords_dreg_v(je,3,2,jk,jb) - z_coords_dreg_v(je,2,2,jk,jb)

!            v21_x = z_coords_dreg_v(je,1,1,jk,jb) - z_coords_dreg_v(je,2,1,jk,jb)
!            v21_y = z_coords_dreg_v(je,1,2,jk,jb) - z_coords_dreg_v(je,2,2,jk,jb)
!            ccw_check = (v23_x * v21_y) - (v23_y * v21_x)

!            IF (ccw_check < 0._wp) THEN
!              print*, 'wrong numbering, ccw_check, je,jk,jb ',ccw_check, je,jk,jb

!              z_dummy(1:2) = z_coords_dreg_v(je,2,1:2,jk,jb)
!              z_coords_dreg_v(je,2,1:2,jk,jb) = z_coords_dreg_v(je,4,1:2,jk,jb)
!              z_coords_dreg_v(je,4,1:2,jk,jb) = z_dummy(1:2)
!            ENDIF
!          ENDDO
!        ENDDO
!      ENDDO


!      !
!      ! This section has been included for testing purposes
!      !
!      SELECT CASE (itime_scheme)
!       !------------------
!       ! Pure advection
!       !------------------
!       CASE (TRACER_ONLY)

!         SELECT CASE ( TRIM(ctest_name) )
!           CASE ('DF1', 'DF2', 'DF3', 'DF4') ! deformational flow

!!DR           z_cell_indices(:,:,:,:) = df_cell_indices(:,:,:,:)
!!DR           z_distv_bary(:,:,:,:)   = df_distv_barycenter(:,:,:,:)

!            ! compute difference between distance vector based on 3rd order
!            ! Taylor series approximation and operational distance vector of either
!            ! first or second order accurracy.
!            df_distv_barycenter(:,:,:,:) = df_distv_barycenter(:,:,:,:)      &
!              &                          - z_distv_bary(:,:,:,:)

!            ! now compute length of distance vectors and store them in the first
!            ! position of the df_distv_barycenter field.
!            df_distv_barycenter(:,:,:,1) = SQRT(df_distv_barycenter(:,:,:,1) &
!              &                          * df_distv_barycenter(:,:,:,1)      &
!              &                          + df_distv_barycenter(:,:,:,2)      &
!              &                          * df_distv_barycenter(:,:,:,2))

!         END SELECT
!      END SELECT


      ! maps quadrilateral onto the standard rectangle of edge length 2.
      ! provides quadrature points and the corresponding determinant of the
      ! Jacobian for each departure region.
      IF (lsq_high_ord == 2) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a quadratic 2D polynomial
        CALL prep_gauss_quadrature_q( p_patch, z_coords_dreg_v,           &! in
          &                      z_quad_vector_sum, z_dreg_area,          &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

      ELSE IF (lsq_high_ord == 30) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a cubic 2D polynomial without cross derivatives
        CALL prep_gauss_quadrature_cpoor( p_patch, z_coords_dreg_v,       &! in
          &                      z_quad_vector_sum, z_dreg_area,          &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

      ELSE IF (lsq_high_ord == 3) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a full cubic 2D polynomial
        CALL prep_gauss_quadrature_c( p_patch, z_coords_dreg_v,           &! in
          &                      z_quad_vector_sum, z_dreg_area,          &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in
      ENDIF


    END IF ! ld_compute


    ! Pointer to line and block indices of the cell center upstream of the edge
    ptr_ilc => z_cell_indices(:,:,:,1)
    ptr_ibc => z_cell_indices(:,:,:,2)



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

!$OMP WORKSHARE
      p_out_e(:,:,i_startblk:i_endblk) = 0._wp
!$OMP END WORKSHARE
    ENDIF

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    ! initialize also nest boundary points with zero
    IF ( l_out_edgeval .AND. (p_patch%id > 1 .OR. l_limited_area) ) THEN
!$OMP WORKSHARE
      p_out_e(:,:,1:i_startblk) = 0._wp
!$OMP END WORKSHARE
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

!CDIR UNROLL=4
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=6
            p_out_e(je,jk,jb) =                                                       &
              &  DOT_PRODUCT(z_lsq_coeff(ptr_ilc(je,jk,jb),jk,1:6,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:6,jk,jb) ) / z_dreg_area(je,jk,jb)

          ENDDO
        ENDDO

        CASE( 30 )  ! cubic reconstruction without third order cross derivatives

!CDIR UNROLL=5
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=8
            p_out_e(je,jk,jb) =                                                       &
              &  DOT_PRODUCT(z_lsq_coeff(ptr_ilc(je,jk,jb),jk,1:8,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:8,jk,jb) ) / z_dreg_area(je,jk,jb)

          ENDDO
        ENDDO

        CASE( 3 )  ! cubic reconstruction with third order cross derivatives

!CDIR UNROLL=4
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=10
            p_out_e(je,jk,jb) =                                                        &
              &  DOT_PRODUCT(z_lsq_coeff(ptr_ilc(je,jk,jb),jk,1:10,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:10,jk,jb) ) / z_dreg_area(je,jk,jb)

          ENDDO
        ENDDO

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

!CDIR UNROLL=4
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=6
            p_out_e(je,jk,jb) =                                                       &
              &  DOT_PRODUCT(z_lsq_coeff(ptr_ilc(je,jk,jb),jk,1:6,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:6,jk,jb) ) / z_dreg_area(je,jk,jb)            &
              &  * p_mass_flx_e(je,jk,jb)

          ENDDO
        ENDDO

        CASE( 30 )  ! cubic reconstruction without third order cross derivatives

!CDIR UNROLL=5
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=8
            p_out_e(je,jk,jb) =                                                       &
              &  DOT_PRODUCT(z_lsq_coeff(ptr_ilc(je,jk,jb),jk,1:8,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:8,jk,jb) ) / z_dreg_area(je,jk,jb)            &
              &  * p_mass_flx_e(je,jk,jb)

          ENDDO
        ENDDO

        CASE( 3 )  ! cubic reconstruction with third order cross derivatives

!CDIR UNROLL=4
        DO jk = slev, elev
          DO je = i_startidx, i_endidx

!CDIR EXPAND=10
            p_out_e(je,jk,jb) =                                                        &
              &  DOT_PRODUCT(z_lsq_coeff(ptr_ilc(je,jk,jb),jk,1:10,ptr_ibc(je,jk,jb)), &
              &  z_quad_vector_sum(je,1:10,jk,jb) ) / z_dreg_area(je,jk,jb)            &
              &  * p_mass_flx_e(je,jk,jb)

          ENDDO
        ENDDO

        END SELECT

      ENDIF
    ENDDO  ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    !
    ! 4. If desired, apply a (semi-)monotone flux limiter to limit computed fluxes.
    !    The flux limiter is based on work by Zalesak (1979)
    !
    IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_m) THEN
      CALL hflx_limiter_mo( p_patch, p_int, p_dtime, p_cc, p_mass_flx_e, & !in
        &                p_out_e, opt_rlend=i_rlend, opt_slev=slev,      & !inout,in
        &                opt_elev=elev                                   ) !in
    ELSE IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_sm) THEN
      ! no MPI-sync necessary
      CALL hflx_limiter_sm( p_patch, p_int, p_dtime, p_cc, p_out_e,      & !in,inout
        &                   opt_rlend=i_rlend, opt_slev=slev,            & !in
        &                   opt_elev=elev                                ) !in
    ENDIF


    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays for quadrature, departure region and
      ! upwind cell indices
      DEALLOCATE( z_quad_vector_sum, z_dreg_area, z_cell_indices,    &
        &         STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                                 &
          &  'deallocation for z_quad_vector_sum, z_dreg_area, ' //  &
          &  ' z_cell_indices failed' )
      ENDIF
    END IF


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
  !!   Grids: High-order Recosntructions for Forward-in-Time Schemes, Mon. Wea. Rev.
  !!
  SUBROUTINE upwind_hflux_ffsl( p_patch, p_cc, p_mass_flx_e, p_vn, p_dtime,     &
    &                      p_int, ld_compute, ld_cleanup, p_itype_hlimit,       &
    &                      p_out_e, opt_rlstart, opt_rlend, opt_lout_edge,      &
    &                      opt_real_vt, opt_slev, opt_elev )


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_hflux: upwind_hflux_ffsl'

    TYPE(t_patch), INTENT(IN)     ::  &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) ::  &  !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), INTENT(IN) ::    &    !< cell centered variable to be advected
      &  p_cc(:,:,:)                !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::    &    !< contravariant horizontal mass flux at cell edge
      &  p_mass_flx_e(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) ::    &    !< normal component of velocity field
      &  p_vn(:,:,:)                !< dim: (nproma,nlev,nblks_e)

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

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'edge value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: output edge value (.TRUE.),
     & opt_lout_edge                   !< or the flux across the edge (.FALSE./not specified)

    REAL(wp), INTENT(IN), OPTIONAL,  & !< optional: tangential velocity
     & TARGET :: opt_real_vt(:,:,:)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    LOGICAL  :: l_out_edgeval          !< corresponding local variable; default .FALSE.
                                       !< i.e. output flux across the edge

   REAL(wp) ::   &                     !< coefficients of lsq reconstruction
      &  z_lsq_coeff(nproma,p_patch%nlev,lsq_high_set%dim_unk+1,p_patch%nblks_c) 
                                       !< at cell center
                                       !< includes c0 and gradients in zonal and
                                       !< meridional direction

    REAL(wp), TARGET ::   &                    !< tangential component of velocity field
      &  z_real_vt(nproma,p_patch%nlev,p_patch%nblks_e)!< at edges

    REAL(wp), POINTER :: ptr_real_vt(:,:,:)    !< pointer to z_real_vt or opt_real_vt


    REAL(wp) ::                &       !< coordinates of arrival points. The origin
      &  arrival_pts(nproma,2,2,p_patch%nlev,p_patch%nblks_e)
                                       !< of the coordinate system is at the circumcenter of
                                       !< the upwind cell. Unit vectors point to local East
                                       !< and North. (geographical coordinates)
                                       !< dim: (nproma,nlev,ptr_p%nblks_e,2,2)

    REAL(wp) ::  &                    !< coordinates of departure points. The origin
      &  depart_pts(nproma,2,2,p_patch%nlev,p_patch%nblks_e)
                                      !< of the coordinate system is at the circumcenter of
                                      !< the upwind cell. Unit vectors point to local East
                                      !< and North. (geographical coordinates)
                                      !< dim: (nproma,nlev,ptr_p%nblks_e,2,2)

    REAL(wp) ::  &                    !< patch 0,1,2 of subdivided departure region
      & dreg_patch0(nproma,4,2,p_patch%nlev,p_patch%nblks_e), &  !< coordinates
      & dreg_patch1(nproma,4,2,p_patch%nlev,p_patch%nblks_e), &
      & dreg_patch2(nproma,4,2,p_patch%nlev,p_patch%nblks_e)


    REAL(wp), ALLOCATABLE, SAVE ::   & !< gauss quadrature vector for each patch
      &  z_quad_vector_sum0(:,:,:,:),& !< dim: (nproma,lsq_dim_unk+1,nlev,nblks_e)
      &  z_quad_vector_sum1(:,:,:,:),&
      &  z_quad_vector_sum2(:,:,:,:)

    REAL(wp), ALLOCATABLE, SAVE ::  & !< area of each departure region patch
      &  z_dreg_area0(:,:,:),       & !< dim: (nproma,nlev,nblks_e)
      &  z_dreg_area1(:,:,:),       &
      &  z_dreg_area2(:,:,:)

    INTEGER, ALLOCATABLE, SAVE, TARGET ::  & !< line and block indices of underlying cell
      & patch0_cell_idx(:,:,:), patch0_cell_blk(:,:,:), & !< dim: (nproma,nlev,p_patch%nblks_e)
      & patch1_cell_idx(:,:,:), patch1_cell_blk(:,:,:), &
      & patch2_cell_idx(:,:,:), patch2_cell_blk(:,:,:)

    INTEGER, POINTER ::                    & !< Pointer to line and block indices of the cells
      &  ptr_ilc0(:,:,:), ptr_ibc0(:,:,:), & !< to which the departure region patches belong.
      &  ptr_ilc1(:,:,:), ptr_ibc1(:,:,:), &
      &  ptr_ilc2(:,:,:), ptr_ibc2(:,:,:)

    INTEGER  :: nlev               !< number of full levels
    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: ist                !< status variable
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    INTEGER  :: dim_unk
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_rlend_c, i_nchdom
    INTEGER  :: pid                !< patch ID

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

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF ( PRESENT(opt_real_vt) ) THEN
      ptr_real_vt => opt_real_vt
    ELSE
      ptr_real_vt => z_real_vt
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
    END IF


    !
    ! 1. reconstruction of the tracer subgrid distribution
    !    least squares method
    !    Note: for rlstart=2 we run into a sync-error with nests
    !
    IF (lsq_high_ord == 2) THEN
      ! quadratic reconstruction
      ! (computation of 6 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_q_svd( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_q( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,        &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ELSE IF (lsq_high_ord == 30) THEN
      ! cubic reconstruction without cross derivatives
      ! (computation of 8 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_cpoor_svd( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,&
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_cpoor( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ELSE IF (lsq_high_ord == 3) THEN
      ! cubic reconstruction with cross derivatives
      ! (computation of 10 coefficients -> z_lsq_coeff )
      IF (advection_config(pid)%llsq_svd) THEN
      CALL recon_lsq_cell_c_svd( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,    &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ELSE
      CALL recon_lsq_cell_c( p_cc, p_patch, p_int%lsq_high, z_lsq_coeff,        &
        &                    opt_slev=slev, opt_elev=elev, opt_rlend=i_rlend_c, &
        &                    opt_rlstart=2 )
      ENDIF
    ENDIF


    ! Synchronize polynomial coefficients
    ! Note: a special sync routine is needed here because the fourth dimension
    ! of z_lsq_coeff is (for efficiency reasons) on the third index
    CALL sync_patch_array_4de3(SYNC_C,p_patch,lsq_high_set%dim_unk+1,z_lsq_coeff)


    !
    ! 2. Approximation of the 'departure region'. The coordinates of
    !    all vertices are computed and stored in an edge-based data
    !    structure.
    !    In addition the Gauss-Legendre quadrature is prepared by
    !    calculating some tracer-invariant (i.e. purely geometric) fields.
    !
    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

    IF (ld_compute) THEN

      IF (.NOT. PRESENT(opt_real_vt)) THEN
        ! reconstruct tangential velocity component at edge midpoints
        CALL rbf_vec_interpol_edge( p_vn, p_patch, p_int,          &! in
          &                         ptr_real_vt, opt_rlend=i_rlend )! inout
      ENDIF

      ! compute vertex coordinates for the departure region using a first
      ! order accurate (O(\Delta t)) backward trajectory-method
      CALL btraj_dreg_nosort( p_patch, p_int, p_vn, ptr_real_vt, p_dtime, &! in
        &                     patch0_cell_idx, patch0_cell_blk,           &! out
        &                     arrival_pts, depart_pts,                    &! out
        &                     opt_rlstart=i_rlstart, opt_rlend=i_rlend,   &! in
        &                     opt_slev=slev, opt_elev=elev                )! in



      ! Flux area (aka. departure region) is subdivided according to its overlap 
      ! with the underlying grid.
      CALL divide_flux_area(p_patch, p_int, p_vn, ptr_real_vt,             &! in
        &                   depart_pts, arrival_pts,                       &! in
        &                   dreg_patch0, dreg_patch1, dreg_patch2,         &! out
        &                   patch1_cell_idx, patch1_cell_blk,              &! out
        &                   patch2_cell_idx, patch2_cell_blk,              &! out
        &                   opt_rlstart=i_rlstart, opt_rlend=i_rlend,      &! in
        &                   opt_slev=slev, opt_elev=elev                   )! in



      ! maps quadrilateral onto the standard rectangle of edge length 2.
      ! provides quadrature points and the corresponding determinant of the
      ! Jacobian for each departure region.
      ! This is done for each of the three patch fragments.
      IF (lsq_high_ord == 2) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a quadratic 2D polynomial
        CALL prep_gauss_quadrature_q( p_patch, dreg_patch0,               &! in
          &                      z_quad_vector_sum0, z_dreg_area0,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

        CALL prep_gauss_quadrature_q( p_patch, dreg_patch1,               &! in
          &                      z_quad_vector_sum1, z_dreg_area1,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

        CALL prep_gauss_quadrature_q( p_patch, dreg_patch2,               &! in
          &                      z_quad_vector_sum2, z_dreg_area2,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

      ELSE IF (lsq_high_ord == 30) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a cubic 2D polynomial without cross derivatives
        CALL prep_gauss_quadrature_cpoor( p_patch, dreg_patch0,           &! in
          &                      z_quad_vector_sum0, z_dreg_area0,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

        CALL prep_gauss_quadrature_cpoor( p_patch, dreg_patch1,           &! in
          &                      z_quad_vector_sum1, z_dreg_area1,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

        CALL prep_gauss_quadrature_cpoor( p_patch, dreg_patch2,           &! in
          &                      z_quad_vector_sum2, z_dreg_area2,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

      ELSE IF (lsq_high_ord == 3) THEN
        ! Gauss-Legendre quadrature with 4 quadrature points for integrating
        ! a full cubic 2D polynomial
        CALL prep_gauss_quadrature_c( p_patch, dreg_patch0,               &! in
          &                      z_quad_vector_sum0, z_dreg_area0,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

        CALL prep_gauss_quadrature_c( p_patch, dreg_patch1,               &! in
          &                      z_quad_vector_sum1, z_dreg_area1,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

        CALL prep_gauss_quadrature_c( p_patch, dreg_patch2,               &! in
          &                      z_quad_vector_sum2, z_dreg_area2,        &! out
          &                      opt_rlstart=i_rlstart, opt_rlend=i_rlend,&! in
          &                      opt_slev=slev, opt_elev=elev             )! in

      ENDIF

    END IF ! ld_compute


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
    !    triangles is at least partly taken into account. Depending on the   
    !    The fact that the rhomboidal area inevitably overlaps with neighboring
    !    triangles is neglected (local quadratic/cubic approximation instead
    !    of piecewise quadratic/cubic approximation). Only the reconstruction for
    !    the local cell is taken into account.

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
      CASE( 2 )  ! quadratic reconstruction

!CDIR UNROLL=5
      DO jk = slev, elev
        DO je = i_startidx, i_endidx

!CDIR EXPAND=6
          p_out_e(je,jk,jb) =                                                            &
            &   ( DOT_PRODUCT(z_lsq_coeff(ptr_ilc0(je,jk,jb),jk,1:6,ptr_ibc0(je,jk,jb)), &
            &     z_quad_vector_sum0(je,1:6,jk,jb) )                                     &
            &   + DOT_PRODUCT(z_lsq_coeff(ptr_ilc1(je,jk,jb),jk,1:6,ptr_ibc1(je,jk,jb)), &
            &     z_quad_vector_sum1(je,1:6,jk,jb) )                                     &
            &   + DOT_PRODUCT(z_lsq_coeff(ptr_ilc2(je,jk,jb),jk,1:6,ptr_ibc2(je,jk,jb)), &
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
            &   ( DOT_PRODUCT(z_lsq_coeff(ptr_ilc0(je,jk,jb),jk,1:8,ptr_ibc0(je,jk,jb)), &
            &     z_quad_vector_sum0(je,1:8,jk,jb) )                                     &
            &   + DOT_PRODUCT(z_lsq_coeff(ptr_ilc1(je,jk,jb),jk,1:8,ptr_ibc1(je,jk,jb)), &
            &     z_quad_vector_sum1(je,1:8,jk,jb) )                                     &
            &   + DOT_PRODUCT(z_lsq_coeff(ptr_ilc2(je,jk,jb),jk,1:8,ptr_ibc2(je,jk,jb)), &
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
            &   ( DOT_PRODUCT(z_lsq_coeff(ptr_ilc0(je,jk,jb),jk,1:10,ptr_ibc0(je,jk,jb)),&
            &     z_quad_vector_sum0(je,1:10,jk,jb) )                                    &
            &   + DOT_PRODUCT(z_lsq_coeff(ptr_ilc1(je,jk,jb),jk,1:10,ptr_ibc1(je,jk,jb)),&
            &     z_quad_vector_sum1(je,1:10,jk,jb) )                                    &
            &   + DOT_PRODUCT(z_lsq_coeff(ptr_ilc2(je,jk,jb),jk,1:10,ptr_ibc2(je,jk,jb)),&
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
    ! 4. If desired, apply a (semi-)monotone flux limiter to limit computed fluxes.
    !    The flux limiter is based on work by Zalesak (1979)
    !
    IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_m) THEN
      CALL hflx_limiter_mo( p_patch, p_int, p_dtime, p_cc, p_mass_flx_e, & !in
        &                p_out_e, opt_rlend=i_rlend, opt_slev=slev,      & !inout,in
        &                opt_elev=elev                                   ) !in
    ELSE IF (.NOT. l_out_edgeval .AND. p_itype_hlimit == ifluxl_sm) THEN
      ! no MPI-sync necessary
      CALL hflx_limiter_sm( p_patch, p_int, p_dtime, p_cc, p_out_e,      & !in,inout
        &                   opt_rlend=i_rlend, opt_slev=slev,            & !in
        &                   opt_elev=elev                                ) !in
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


  !-----------------------------------------------------------------------
  !>
  !! The upwind biased 3rd oder advection for the hexagonal grid.
  !! It would be generizable for triangles, too.
  !!
  !! @par Revision History
  !! Developed by Almut Gassmann, MPI-M (2010-11-18)
  !!
  !! @par !LITERATURE
  !! - Skamarock and Gassmann, MWR (to be published)
  !!
  SUBROUTINE upwind_hflux_hex( p_patch, p_int, p_cc, p_c0, p_mass_flx_e,   &
    &                          p_dtime, p_itype_hlimit, p_out_e, opt_slev, &
    &                          opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &    !< patch on which computation is performed
      &  p_patch

    TYPE(t_int_state), TARGET, INTENT(IN) :: & !< pointer to data structure for interpolation
      &  p_int

    REAL(wp), INTENT(IN) ::     &   !< cell centered variable to be advected
      &  p_cc(:,:,:)                !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &   !< advected cell centered variable (step n)
      &  p_c0(:,:,:)                !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::     &   !< contravariant horizontal mass flux
      &  p_mass_flx_e(:,:,:)        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(INOUT) ::  &   !< variable in which the upwind flux is stored
      &  p_out_e(:,:,:)             !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN) :: p_dtime !< time step
 
    INTEGER, INTENT(IN) ::      &   !< parameter to select the limiter
      &  p_itype_hlimit             !< for horizontal transport

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    REAL(wp) :: z_ave_e     (nproma,p_patch%nlev,p_patch%nblks_e), &
      &         z_dir_lapl_e(nproma,p_patch%nlev,p_patch%nblks_e)

    INTEGER  :: slev, elev         !< vertical start and end level
    INTEGER  :: nblks_e, npromz_e
    INTEGER  :: jk, jb, nlen       !< index vert level, block; length of block
    INTEGER  :: jg                 !< patch ID

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

    ! get patch ID
    jg = p_patch%id

    ! compute ordinary average at the edge
    CALL cells2edges_scalar(p_cc,p_patch,p_int%c_lin_e,z_ave_e, &
      &                     opt_slev=slev, opt_elev=elev )

    ! compute directional laplace in edge direction
    CALL directional_laplace(p_mass_flx_e,p_cc,p_patch,p_int,&
      &                      advection_config(jg)%upstr_beta_adv,z_dir_lapl_e,&
      &                      opt_slev=slev, opt_elev=elev )

    ! values for the blocking
    nblks_e  = p_patch%nblks_int_e
    npromz_e = p_patch%npromz_int_e

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = slev, elev
          p_out_e(1:nlen,jk,jb) = p_mass_flx_e(1:nlen,jk,jb) &
          & *(z_ave_e(1:nlen,jk,jb) &
          & - p_patch%edges%dual_edge_length(1:nlen,jb)  &
          & * p_patch%edges%dual_edge_length(1:nlen,jb)  &
          & /6.0_wp * z_dir_lapl_e(1:nlen,jk,jb))
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF ( p_itype_hlimit == ifluxl_sm) THEN
      CALL hflx_limiter_sm( p_patch, p_int, p_dtime, p_c0, p_out_e, & !in,inout
        &                   opt_slev=slev, opt_elev=elev            ) !in
    ENDIF

  END SUBROUTINE upwind_hflux_hex

END MODULE mo_advection_hflux


